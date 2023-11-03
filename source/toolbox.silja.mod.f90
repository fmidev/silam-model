MODULE toolbox
  !
  ! Description:
  ! Contains small useful functions that are needed in several higher
  ! -level modules of silja.
  !
  ! 1. lat-lon-scaling and checking
  ! 2. one-dimensional interpolation
  !
  ! Author: Mika Salonoja email Mika.Salonoja@fmi.fi
  !
  ! All units: SI
  ! 
  ! All degrees: real numbers in degrees (degrees and parts per
  ! hundred/parts per thousand , NOT in degrees and minutes
  ! 
  ! Language: ANSI standard Fortran 90
  ! 
  ! Modules used:
  USE md
  use work_arrays

!  use thermodynamic_tools

  IMPLICIT NONE

  ! The public functions and subroutines available in this module:
  public fu_index
  PUBLIC fu_set_latitude
  PUBLIC fu_scale_latitude
  PUBLIC fu_scale_longitude
  PUBLIC weight_coefficients
  public fu_value_index_in_array
  PUBLIC distance_metres_to_degrees
  PUBLIC distance_degrees_to_metres
  public fu_distance_m
  PUBLIC fu_dx_deg_to_m
  PUBLIC fu_dx_m_to_deg
  PUBLIC fu_dy_deg_to_m
  PUBLIC fu_dy_m_to_deg
  PUBLIC set_array_size
  PUBLIC free_array_memory
  PUBLIC fu_minimum_total_memory_usage
  PUBLIC fu_random_number_center
  PUBLIC fu_random_number_center_rng
  PUBLIC rng_init
  PUBLIC fu_random_number_boundaries
  public random_normal
  PUBLIC next_line_from_input_file
  PUBLIC fu_next_free_unit
  PUBLIC fu_merge_int_arrays
  public fu_sort_real_array
  public fu_2d_interpolation
  public copy_text_file
  public fu_u_case
  public fu_l_case
  public str_u_case
  public str_l_case
  public fu_str_u_case
  public fu_str_l_case
  public compress_int_array ! Remove holes
  public fu_nbrOfWords  ! in the string
  public fu_if_digit
  public fu_solar_zenith_angle_cos
  public fu_day_length_hrs
  public fu_expand_environment
  public init_random_seed
  public fu_merge_integer_to_array
  public matrix_exponent
  public fu_matrix_norm 
  public advect_Galperin_bulk_1d_relat
  public advect_Galperin_bulk_1d_abs
  public fu_process_filepath
  public split_string
  public fu_trim_grads_hat   ! replaces a comon path with a hat "^" symbol
                             ! if dir_slash == '/' (i.e. UNIXes)
  public fu_extend_grads_hat ! extends a hat "^" symbol with ctl file path
  public fu_extend_grads_hat_dir ! extends a hat "^" symbol with given dir
  public fu_if_ends_with
  public replace_string
  public fu_dirname
  public fu_basename
  public linear_regression
  public fu_studnt 
  public fu_lognorm_density
  public get_chunk
  public set_environment_variable
  public fu_get_default_mpi_buf_size
  public remapcon_1d
  public testremapcon_1d
  public set_coastal_value

  public trim_precision
  public fu_trim_cell_fcoord

  public test_sort
  public argMergeSort_char
  public argMergeSort_int
  public argMergeSort_real

  ! The private functions and subroutines of this module:
  private fu_index_of_string
  private fu_index_of_int
  PRIVATE set_arraysize_realpointer
  PRIVATE set_arraysize_realpointer_2d
  PRIVATE free_array_memory_realpointer
  PRIVATE fu_compare_reals
  private fu_compare_double_reals
  private fu_compare_double_to_single
  private fu_compare_single_to_double
  private fu_2d_interp_to_coord
  
  private fu_solar_zenith_angle_cos_par
  private split_string_int
  private split_string_real
  private split_string_char
  private fu_next_non_sep
  private fu_next_sep

  public open_grads_binary_i
  public open_grads_binary_o
  public write_grads_field
  public write_grads_fieldset
  
  interface fu_index
    module procedure fu_index_of_string
    module procedure fu_index_of_int
  end interface

  interface split_string
     module procedure split_string_int
     module procedure split_string_real
     module procedure split_string_char
  end interface

  INTERFACE operator(.eps.)
    MODULE PROCEDURE fu_compare_reals
    MODULE PROCEDURE fu_compare_double_reals
    MODULE PROCEDURE fu_compare_double_to_single
    MODULE PROCEDURE fu_compare_single_to_double
  END INTERFACE

  INTERFACE operator(.deps.)
    MODULE PROCEDURE fu_compare_double_reals
  END INTERFACE

  INTERFACE  set_array_size
    MODULE PROCEDURE set_arraysize_realpointer
    MODULE PROCEDURE set_arraysize_realpointer_2d
  END INTERFACE

  INTERFACE free_array_memory
    MODULE PROCEDURE free_array_memory_realpointer
  END INTERFACE

  interface fu_merge_integer_to_array
    module procedure fu_merge_int_2_arr_with_slave
    module procedure fu_merge_int_to_array
  end interface

  interface fu_2d_interpolation
    module procedure fu_2d_interp_to_coord 
  end interface

  interface fu_solar_zenith_angle_cos
    module procedure fu_solar_zenith_angle_cos_par
  end interface

  interface
#ifdef VS2012
     function setenv(name,value,overwrite) bind(C,name='_putenv')
#else
     function setenv(name,value,overwrite) bind(C,name='setenv')
#endif
        use ISO_C_BINDING
        implicit none
        integer(C_INT) setenv
        character(KIND=C_CHAR), intent(in) :: name(*)
        character(KIND=C_CHAR), intent(in) :: value(*)
        integer(C_INT), value :: overwrite
     end function setenv
  end interface


  REAL, SAVE, PRIVATE :: total_memory_usage = 0

  INTEGER, PARAMETER, PUBLIC :: ascending = 1000 ! Order of sorting
  INTEGER, PARAMETER, PUBLIC :: descending = 2000 ! Order of sorting

  integer, public, parameter :: pressure_unit = 10008  ! Temperature

  integer, private, parameter :: ai = iachar('a')
  integer, private, parameter :: zi = iachar('z')
  integer, private, parameter :: diff = iachar('z') - iachar('Z')


  !Thread-safe random-number generator borrowed from
  !http://jblevins.org/log/openmp
  ! Dimension of the state
  integer, private, parameter :: ns = 4

  ! Default seed vector
  integer, private, parameter, dimension(ns) :: default_seed &
       = (/ 521288629, 362436069, 16163801, 1131199299 /)

  ! A data type for storing the state of the RNG
 
  type :: rng_type
     private
     integer, dimension(ns) :: state
  end type rng_type
  public rng_type

  type fnlen_str
    character(len=fnlen) :: s
  end type fnlen_str

CONTAINS


  ! ****************************************************************
  ! ****************************************************************
  ! 
  !
  !          LATITUDE-LONGITUDE STUFF
  !
  !
  ! ***************************************************************
  ! ***************************************************************

  LOGICAL ELEMENTAL FUNCTION fu_compare_reals(real1, real2)
    !
    ! Compares two real and returns a true value if they're close.
    !
    ! if real1/real2 deviates from unity by less than epsilon, then a
    ! true value is returned.
    !
    IMPLICIT NONE
    !
    ! Imported parameters with intent(in):
    REAL(r4k), INTENT(in) :: real1, real2
    
    ! Local declarations:
    !    REAL, PARAMETER :: eps = 1.0e-10
    REAL(r4k), PARAMETER :: eps = 5.0e-5

    IF (abs(real2) < eps) THEN
      fu_compare_reals = (real1 < eps .and. real1 > -eps)

    ELSE

      IF ((ABS((real1/real2) - 1.)) < eps) THEN
        fu_compare_reals = .true.
      ELSE
        fu_compare_reals = .false.
      END IF

    END IF

  END FUNCTION fu_compare_reals


  ! ****************************************************************

  LOGICAL ELEMENTAL FUNCTION fu_compare_double_reals(real1, real2)
    !
    ! Compares two real and returns a true value if they're close.
    !
    ! if real1/real2 deviates from unity by less than epsilon, then a
    ! true value is returned.
    !
    IMPLICIT NONE
    !
    ! Imported parameters with intent(in):
    real(r8k), INTENT(in) :: real1, real2
    
    ! Local declarations:
    REAL(r8k), PARAMETER :: eps = 1.0e-10
!    real(r8k), PARAMETER :: eps = 5.0e-5

    IF (real2 == 0.) THEN
      fu_compare_double_reals = (real1 < eps .and. real1 > -eps)

    ELSE

      IF ((ABS((real1/real2) - 1.)) < eps) THEN
        fu_compare_double_reals = .true.
      ELSE
        fu_compare_double_reals = .false.
      END IF

    END IF

  END FUNCTION fu_compare_double_reals
  
  ! ****************************************************************

  LOGICAL FUNCTION fu_compare_double_to_single(real_d, real_s)
    !
    ! Compares two real and returns a true value if they're close.
    !
    ! if real1/real2 deviates from unity by less than epsilon, then a
    ! true value is returned.
    !
    IMPLICIT NONE
    !
    ! Imported parameters with intent(in):
    real(r8k), intent(in) :: real_d
    real(r4k), intent(in) :: real_s

    fu_compare_double_to_single = fu_compare_double_reals(real_d, real(real_s,8))

  END FUNCTION fu_compare_double_to_single

  ! ****************************************************************

  LOGICAL FUNCTION fu_compare_single_to_double(real_s_first_arg, real_d_second_arg)
    !
    ! Compares two real and returns a true value if they're close.
    !
    ! if real1/real2 deviates from unity by less than epsilon, then a
    ! true value is returned.
    !
    IMPLICIT NONE
    !
    ! Imported parameters with intent(in):
    real(r8k), intent(in) :: real_d_second_arg
    real(r4k), intent(in) :: real_s_first_arg

    fu_compare_single_to_double = fu_compare_double_reals(real(real_s_first_arg,8), real_d_second_arg)

  END FUNCTION fu_compare_single_to_double



  ! ****************************************************************
  
  REAL FUNCTION  fu_set_latitude(lat, sn_flag)
    !
    ! Description:
    ! Sets a value for a latitude, geographical or modified,
    ! according to the silja-standard: northern latitudes are always
    ! positive.
    ! 
    ! Author: Mika Salonoja, FMI
    ! 
    IMPLICIT NONE
    !
    ! Imported parameters with intent IN:
    REAL, INTENT(in) :: lat
    INTEGER, INTENT(in) :: sn_flag

    SELECT CASE (sn_flag)

    CASE (north_flag)
      fu_set_latitude= lat

    CASE (south_flag)
      fu_set_latitude = -1. * lat

    CASE default

      CALL set_error( 'latitude must be south or north!','fu_set_latitude')
      fu_set_latitude = real_missing
      RETURN

    END SELECT

    IF ((fu_set_latitude > 90.0001) .or. (fu_set_latitude < -90.0001)) THEN
      CALL set_error('the value of latitude out of bounds:' + fu_str(fu_set_latitude),&
                   & 'fu_set_latitude')
      fu_set_latitude = real_missing
    else
      fu_set_latitude = min(90.0,max(-90.0,fu_set_latitude))
    END IF
    
  END FUNCTION fu_set_latitude
  
  
  ! ***************************************************************
  
  
  REAL FUNCTION fu_scale_latitude(lat)
    
    ! Scales a latitude value between -90.0...+90.0 degrees
    ! DUMMY VERSION, DOES NOTHING!
    ! Author: Mika Salonoja, FMI
    ! 
    ! All  units: SI
    ! 
    ! Language : ANSI Fortran 90
    
    IMPLICIT NONE
    !
    ! Imported parameters with intent IN:
    REAL, INTENT(in)  :: lat

    fu_scale_latitude = lat
    
  END FUNCTION fu_scale_latitude
  
  
  ! ************************************************************
  
  REAL FUNCTION fu_scale_longitude(lon_in) result(lon)
    
    ! Scales a longitude value between -180.0...+180.0 degrees
    !
    ! Author: Mika Salonoja, FMI
    ! 
    ! All  units: SI
    ! 
    ! Language : ANSI Fortran 90
    
    IMPLICIT NONE
    !
    ! Imported parameters with intent IN:
    REAL, INTENT(in) :: lon_in

    lon = lon_in ! the initial value
    
    if(abs(lon_in) > 365)then
      call msg('Illegal longitude: ', lon_in)
      call set_error('Illegal longitude given',  'fu_scale_longitude')
      return
    endif
    
    IF (lon > 180.) THEN

      DO
        lon = lon - 360.
        IF (lon <= 180.) EXIT
      END DO

    ELSE IF (lon < -180.) THEN
      
      DO
        lon = lon + 360.
        IF (lon >= -180.) EXIT
      END DO
      
    END IF
    
  END FUNCTION fu_scale_longitude





  ! ****************************************************************
  ! ****************************************************************
  ! 
  !
  !          GENERAL INTERPOLATION 
  !
  !
  ! ***************************************************************
  ! ***************************************************************
  


  !*************************************************************************



  ! ***************************************************************
  
  SUBROUTINE weight_coefficients(diff_a, diff_b, method, coeff_a, coeff_b)
    
    ! Description:
    ! Calculates weight coefficients for 1d-interpolation. The
    ! possible methods are given in globals. In interpolation diff_a
    ! and diff_b have opposite signs, and in
    ! extrapolation same signs.
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
    REAL, INTENT(in) :: diff_a, diff_b
    INTEGER, INTENT(in) :: method
    !
    ! Imported parameters with intent OUT:
    REAL, INTENT(out) :: coeff_a, coeff_b
    
    ! Local declarations:
    REAL :: total

    SELECT CASE (method)

    CASE (linear)

      !----------------------------------------
      !
      ! 1. Linear coefficients.

      IF (diff_a == diff_b) THEN
        coeff_a = 0.5
        coeff_b = 0.5
      ELSE
        coeff_a = diff_b / (diff_b - diff_a)
        coeff_b = diff_a / (diff_a - diff_b)

      END IF


    CASE (nearest_point)

      !----------------------------------------
      !
      ! 3. Nearest point.

      IF (ABS(diff_a) < ABS(diff_b)) THEN
        coeff_a = 1.
        coeff_b = 0.
      ELSE
        coeff_a = 0.
        coeff_b = 1.
      END IF


    CASE default

      CALL set_error('unknown interpolation method',&
                   & 'weight_coefficients')

    END SELECT

  END SUBROUTINE weight_coefficients


  !**************************************************************************
  
  real function fu_value_index_in_array(fValue, arValuesSorted, nValues)
    !
    ! Calculates the relative index of a real value in the sorted array of reals
    ! Distance between each two reals is taken linearly. Index below-the-first
    ! is zero, index above-the-last is nValues+1
    !
    implicit none
    
    ! Imported parameters
    real, intent(in) :: fValue
    real, dimension(:), intent(in) :: arValuesSorted
    integer, intent(in) :: nValues
    
    ! Local variables
    integer :: iTmp, dir, iStart, iEnd, iter
    
    !
    ! Stupidity check
    !
    if(nValues < 2 .or. size(arValuesSorted) < nValues)then
      call set_error('Funny nvalues and/or real array size:' + fu_str(nValues) + ',' + &
                   & fu_str(size(arValuesSorted)),'fu_value_index_in_array')
      return
    endif
    !
    ! Index is computed with regard to sorting
    !
    fu_value_index_in_array = 0  ! if outside from below
    if(arValuesSorted(1) < arValuesSorted(nValues))then               ! Ascending sorting
      dir = 1
    else
      dir = -1
    endif

    if(fValue * dir < arValuesSorted(1) * dir)then
      fu_value_index_in_array = 0.
    elseif(fValue * dir > arValuesSorted(nValues) * dir)then
      fu_value_index_in_array = nValues + 1.0
    else
      iStart = 1
      iEnd = nValues
      !
      ! Cutting the whole array by halfs
      !
      do while(iEnd - iStart > 1)
        iter = nint((iStart + iEnd) / 2.)
        if(fValue * dir > arValuesSorted(iter) * dir)then
          iStart = iter
        else
          iEnd = iter
        endif
      enddo
      if(iEnd == iStart)then
        fu_value_index_in_array = iStart
      else
        fu_value_index_in_array = iStart + (fValue - arValuesSorted(iStart)) / &
                                         & (arValuesSorted(iEnd) - arValuesSorted(iStart))
      endif
    endif
  end function fu_value_index_in_array


  subroutine remapcon_1D(bnd1, bnd2, epsi,  bndout, idx1, idx2, nSeg)
      ! Info for conservative remapping of 1D arrays
      ! Creates a new set of bounds "supermesh" and fills the indices 
      ! increasing sequence of both bnd1 and bnd2 required
      ! Number of out bounds is at most sum of input bounds 
      ! Note, that procedure is assymetric: if bounds are within epsi,
      ! bounds from bnd2 are used
      real, dimension(:), intent(in) :: bnd1, bnd2 ! bound points of arr1
      real, intent(in) ::  epsi !! Consider bunds same if they differ by less than eps
          ! Assumed sorted
      real, dimension(:), intent(out) :: bndout ! Common set of bounds
      integer, dimension(:), intent(out) :: idx1, idx2
      integer, intent(out) ::  nseg

      integer :: i1, i2, iOut, n1, n2

      n1 = size(bnd1)
      n2 = size(bnd2)

      i1 = 0   !!! Counters for segments 
      i2 = 0

      do iOut = 1, n1+n2
        if (bnd1(i1+1) + epsi < bnd2(i2+1)) then
          i1 = i1 + 1
          bndout(iOut) = bnd1(i1)
        else
          if (bnd1(i1+1) < bnd2(i2+1) + epsi ) i1 = i1 + 1 ! if bound within epsi -- count it once
          i2 = i2 + 1
          bndout(iOut) = bnd2(i2)
        endif
        idx1(iOut) = i1
        idx2(iOut) = i2
 !       write (*,*) iOut, bndout(iOut), i1, i2
        if (i1 >= n1) exit
        if (i2 >= n2) exit
      enddo
      do i2 = i2+1, n2
         iOut = iOut + 1
!         write (*,*) iOut, bnd2(i2), i1, i2
         bndout(iOut) = bnd2(i2)
         idx1(iOut) = i1
         idx2(iOut) = i2
      enddo
      do i1 = i1+1, n1
         iOut = iOut + 1
  !       write (*,*) iOut, bnd1(i1), i1, i2
         if (iOut > size(bndout)) call ooops("")
         bndout(iOut) = bnd1(i1)
         idx1(iOut) = i1
         idx2(iOut) = i2
      enddo
      nseg = iOut

  end subroutine remapcon_1D

  subroutine testremapcon_1D()
    integer i, n
    integer, parameter :: n1 = 1801
    integer, parameter :: n2 = 1801

    real, dimension (n1) :: bnd1 
    real, dimension (n2) :: bnd2
    real, dimension (n1+n2) :: bndout
    integer, dimension (n1+n2) :: idx1, idx2
    integer :: nSeg

    do i = 1,n1
      bnd1(i) = -180. + 0.2*(i-1)
    enddo
    do i = 1,n2
      bnd2(i) = -180 + 0.2*(i-1)
    enddo

    call msg("bnd1", bnd1)
    call msg("bnd2", bnd2)
    call remapcon_1D(bnd1, bnd2, 0.001, bndout, idx1, idx2, nSeg)
    call msg("Bndout(1)", bndout(1))
    do i = 1, nSeg-1
        call msg ("idx1,idx2", idx1(i), idx2(i))
        call msg("Bndout(i+1), i", bndout(i+1)) 
    enddo

  endsubroutine testremapcon_1D

  ! ***************************************************************
  ! ***************************************************************
  !
  !
  !        Mapping factors
  !
  !
  ! ***************************************************************
  ! ***************************************************************
  
  SUBROUTINE distance_metres_to_degrees(dx, dy, lat, dx_deg, dy_deg)
    
    ! Description:
    ! Converts a given movement (dx, dy) in metres to a movement in
    ! degrees in spherical coordinates on earth's surface. The given
    ! latitude must be in the same spherical coordinate system than
    ! dx and dy (they have to have the same pole)!!
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
    REAL, INTENT(in) :: dx, dy ! in metres
    REAL, INTENT(in) :: lat ! from where dx and dy are (same pole!!)
    !
    ! Imported parameters with intent OUT:
    REAL, INTENT(out) :: dx_deg, dy_deg

    dy_deg = dy/(earth_radius*degrees_to_radians)
    dx_deg = dx/(earth_radius*degrees_to_radians*COS(degrees_to_radians*lat))
    
  END SUBROUTINE distance_metres_to_degrees


  ! ***************************************************************

  SUBROUTINE distance_degrees_to_metres(dx_deg, dy_deg, lat, dx, dy)
    
    ! Description:
    ! Converts a given movement (dx_deg, dy_deg) in degrees (in
    ! spherical coordinates on earth's surface) to a
    ! movement in metres. The given
    ! latitude must be in the same spherical coordinate system than
    ! dx and dy (they have to have the same pole)!!
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
    REAL, INTENT(in) :: dx_deg, dy_deg ! in degrees
    REAL, INTENT(in) :: lat ! from where dx and dy are (same pole!!)
    
    ! Imported parameters with intent OUT:
    REAL, INTENT(out) :: dx, dy ! in metres
    
    dy = dy_deg*(earth_radius*degrees_to_radians)
    dx = dx_deg*(earth_radius*degrees_to_radians*COS(degrees_to_radians*lat))
    
  END SUBROUTINE distance_degrees_to_metres


  !*******************************************************************
  
  real function fu_distance_m(dx_deg, dy_deg, lat)
    REAL, INTENT(in) :: dx_deg, dy_deg ! in degrees
    REAL, INTENT(in) :: lat ! from where dx and dy are (same pole!!)
    
    real :: dx, dy
    
    call distance_degrees_to_metres(dx_deg, dy_deg, lat, dx, dy)

    fu_distance_m = sqrt(dx*dx + dy*dy)

  end function fu_distance_m


  ! ***************************************************************
  
  REAL FUNCTION fu_dx_deg_to_m(dx_deg, lat)
    
    ! Description:
    ! Converts a given movement dx_deg (east-west-direction)
    ! to a movement in metres. The given
    ! latitude must be in the same spherical coordinate system than
    ! dx (they have to have the same pole)!!
    ! 
    ! All units: SI
    !
    ! Language: ANSI Fortran 90
    !
    ! Author: Mika Salonoja, FMI
    ! 
    IMPLICIT NONE
    !
    ! Imported parameters with intent(in):
    REAL, INTENT(in) :: dx_deg, lat
    real :: dTheta

    dTheta = dx_deg*degrees_to_radians

    fu_dx_deg_to_m = earth_radius*COS(degrees_to_radians*lat)*2.*sin(0.5*dTheta)
!    fu_dx_deg_to_m = earth_radius*COS(degrees_to_radians*lat)* dTheta*(1. - dTheta*dTheta/24.)

!    fu_dx_deg_to_m = earth_radius*COS(degrees_to_radians*lat)* dx_deg*degrees_to_radians

  END FUNCTION fu_dx_deg_to_m

  ! ***************************************************************
  
  REAL FUNCTION fu_dx_m_to_deg(dx_m, lat)
    
    ! Description:
    ! Converts a given movement dx in metres (east-west-direction)
    ! to a movement in degrees. The given
    ! latitude must be in the same spherical coordinate system than
    ! dx (they have to have the same pole)!!
    ! 
    ! All units: SI
    !
    ! Language: ANSI Fortran 90
    !
    ! Author: Mika Salonoja, FMI
    ! 
    IMPLICIT NONE
    !
    ! Imported parameters with intent(in):
    REAL, INTENT(in) :: dx_m, lat

    fu_dx_m_to_deg = dx_m/(earth_radius*degrees_to_radians*COS(degrees_to_radians*lat))

  END FUNCTION fu_dx_m_to_deg


  ! ***************************************************************
  
  REAL FUNCTION fu_dy_deg_to_m(dy_deg)
    
    ! Description:
    ! Converts a given movement dy_deg (north-south-direction)
    ! to a movement in metres. 
    ! All units: SI
    !
    ! Language: ANSI Fortran 90
    !
    ! Author: Mika Salonoja, FMI
    ! 
    IMPLICIT NONE
    !
    ! Imported parameters with intent(in):
    REAL, INTENT(in) :: dy_deg

    fu_dy_deg_to_m = dy_deg*earth_radius*degrees_to_radians

  END FUNCTION fu_dy_deg_to_m


  ! ***************************************************************
  
  REAL FUNCTION fu_dy_m_to_deg(dy_m)
    
    ! Description:
    ! Converts a given movement metric dy_m (north-south-direction)
    ! to a movement in degrees 
    ! All units: SI
    !
    ! Language: ANSI Fortran 90
    !
    ! Author: Mika Salonoja, FMI
    ! 
    IMPLICIT NONE

    ! Imported parameters with intent(in):
    REAL, INTENT(in) :: dy_m

    fu_dy_m_to_deg = dy_m/(earth_radius*degrees_to_radians)

  END FUNCTION fu_dy_m_to_deg










  ! ***************************************************************
  ! ***************************************************************
  !
  !
  !        Dynamic memory stuff.
  !
  !
  ! ***************************************************************
  ! ***************************************************************

  SUBROUTINE set_arraysize_realpointer(array, arraysize, fillvalue)
    
    ! Description:
    ! Makes the 1-d real pointer to have the given size. 
    ! The array is then filled with the optional given value. 
    ! If no value is given, the array values are not touched.
    !
    ! Method:
    ! 1. If the pointer is not associated with a real-array (by
    ! pointing to a target or allocating) then a memory is allocated
    ! and the pointer becomes an alias for the real-array.
    ! 2. If there's already memory allocated, and it is of correct
    ! size, then nothing is done. 
    ! 3. If there's already memory allocated, but not correct
    ! size, then memory is first deallocated and then allocated again.
    ! Addition of MAS: values are kept as much as possible !
    ! 
    ! All units: SI
    !
    ! Language: ANSI Fortran 90
    !
    ! Original code: Mika Salonoja
    ! Author: Mikhail Sofiev, FMI
    ! 
    IMPLICIT NONE
    !
    ! Imported parameters with intent IN:
    INTEGER, INTENT(in) :: arraysize
    REAL, INTENT(in), OPTIONAL :: fillvalue
    
    ! Imported parameters with intent INOUT or POINTER:
    REAL, DIMENSION(:), POINTER :: array
    
    ! Local declarations:
    INTEGER :: status, size_old
    REAL, DIMENSION(:), POINTER :: arrayTmp
    logical :: ifNew

    if(ASSOCIATED(array))then
      ifNew = (size(array) <= 0 .or. size(array) > worksize * 10)
    else
      ifNew = .true.
    endif

    new_or_old: IF (ifNew) THEN

      !
      ! 1. Allocate memory for a new, untouched array pointer.
      !
!            PRINT *, ' Allocating a new array.',arraysize

      ALLOCATE(array(arraysize), stat = status)

      IF (status /= 0) THEN
        call msg(' arraysize: ', arraysize)
        call msg(' Status: ', status)
        CALL set_error('cannot allocate memory for new array','set_realpointer_arraysize')
        RETURN
      END IF
      
      total_memory_usage = total_memory_usage + arraysize

      !
      ! 3. Fill the array with the optional value.
      !
      IF (PRESENT(fillvalue)) array = fillvalue

    ELSE 
      
      size_check: IF (SIZE(array) /= arraysize) THEN      
      !
      ! 2. Allocate memory for an array of wrong size.
      !
!		PRINT *, ' Old array, wrong size. Must allocate....',SIZE(array), arraysize

      ! There's memory allocated alreay, but wrong size:

!        IF(.not.PRESENT(fillvalue))then
          size_old = SIZE(array)
          arrayTmp => fu_work_array(max(worksize,size_old))
          arrayTmp(1:size_old) = array
!        endif
      
        DEALLOCATE(array, stat = status)

        IF(fu_fails(status==0, 'error while deallocating array memory','set_realpointer_arraysize'))return
	
        ALLOCATE(array(arraysize), stat = status)

        IF (status /= 0) THEN
          call msg(' arraysize: ', arraysize)
          call msg(' Status: ', status)
          CALL set_error('cannot allocate memory for old array','set_arraysize_realpointer')
          RETURN
        END IF

        total_memory_usage = total_memory_usage + arraysize

!        IF(.not.PRESENT(fillvalue))then
          array(1:min(size_old,size(array))) = arrayTmp(1:min(size_old,size(array)))
          call free_work_array(arrayTmp)
!        endif
        !
        ! 3. Fill the array with the optional value.
        !
        IF (PRESENT(fillvalue))then
          if(size_old < size(array)) array(size_old+1:size(array)) = fillvalue
        endif

      END IF size_check

!           PRINT *, 'Old array, already same size.'

    END IF new_or_old
        
  END SUBROUTINE set_arraysize_realpointer


  !*********************************************************************************

  SUBROUTINE set_arraysize_realpointer_2d(array, arraysize_1, arraysize_2, fillvalue)
    
    ! Description:
    ! Makes the 1-d real pointer to have the given size. 
    ! The array is then filled with the optional given value. 
    ! If no value is given, the array values are not touched.
    !
    ! Method:
    ! 1. If the pointer is not associated with a real-array (by
    ! pointing to a target or allocating) then a memory is allocated
    ! and the pointer becomes an alias for the real-array.
    ! 2. If there's already memory allocated, and it is of correct
    ! size, then nothing is done. 
    ! 3. If there's already memory allocated, but not correct
    ! size, then memory is first deallocated and then allocated again.
    ! Addition of MAS: values are kept as much as possible !
    ! 
    ! All units: SI
    !
    ! Language: ANSI Fortran 90
    !
    ! Original code: Mika Salonoja
    ! Author: Mikhail Sofiev, FMI
    ! 
    IMPLICIT NONE
    !
    ! Imported parameters with intent IN:
    INTEGER, INTENT(in) :: arraysize_1, arraysize_2
    REAL, INTENT(in), OPTIONAL :: fillvalue

    ! Imported parameters with intent INOUT or POINTER:
    REAL, DIMENSION(:,:), POINTER :: array
    
    ! Local declarations:
    INTEGER :: status, size_old_1, size_old_2
    REAL, DIMENSION(:,:), POINTER :: arrayTmp
    logical :: ifNew

    if(ASSOCIATED(array))then
      ifNew = (size(array) <= 0 .or. size(array,1) > worksize*10 .or. size(array,2) > worksize*10)
    else
      ifNew = .true.
    endif

    new_or_old: IF (ifNew) THEN
      !
      ! 1. Allocate memory for a new, untouched array pointer.
      !
!            PRINT *, ' Allocating a new array.',arraysize

      ALLOCATE(array(arraysize_1, arraysize_2), stat = status)

      IF (status /= 0) THEN
        call msg(' arraysizes: ', arraysize_1, arraysize_2)
        call msg(' Status: ', status)
        CALL set_error('cannot allocate memory for new array','set_arraysize_realpointer_2d')
        RETURN
      END IF
      
      total_memory_usage = total_memory_usage + arraysize_1*arraysize_2

      !
      ! 3. Fill the array with the optional value.
      !
      IF (PRESENT(fillvalue)) array = fillvalue

    ELSE 
      
      size_check: IF (SIZE(array,1) /= arraysize_1 .or. SIZE(array,2) /= arraysize_2) THEN      
        !
        ! 2. Allocate memory for an array of wrong size.
        !
!		PRINT *, ' Old array, wrong size. Must allocate....',SIZE(array), arraysize

        ! There's memory allocated alreay, but wrong size:

!        IF(.not.PRESENT(fillvalue))then
          size_old_1 = SIZE(array,1)
          size_old_2 = SIZE(array,2)
          allocate(arrayTmp(size_old_1,size_old_2), stat=status)
          if(status /= 0)then
            call msg(' arraysizes: ', arraysize_1, arraysize_2)
            call msg(' Status: ', status)
            CALL set_error('cannot allocate memory for new array','set_arraysize_realpointer_2d')
            RETURN
          endif
          arrayTmp(1:size_old_1, 1:size_old_2) = array(1:size_old_1,1:size_old_2)
!        endif

        DEALLOCATE(array, stat = status)
        IF (status /= 0) THEN
          CALL set_error('error while deallocating memory for array','set_arraysize_realpointer_2d')
          RETURN
        END IF
	
        ALLOCATE(array(arraysize_1,arraysize_2), stat = status)
        IF (status /= 0) THEN
          call msg(' arraysizes: ', arraysize_1,arraysize_2)
          call msg(' Status: ', status)
          CALL set_error('cannot allocate memory for old array','set_arraysize_realpointer_2d')
          RETURN
        END IF

        total_memory_usage = total_memory_usage + arraysize_1*arraysize_2

!        IF(.not.PRESENT(fillvalue))then
          array(1:min(size_old_1,size(array,1)),1:min(size_old_2,size(array,2))) = arrayTmp(:,:)
          deallocate(arrayTmp,stat=status)
          IF (status /= 0) THEN
            CALL set_error('error while deallocating memory for temporary array', &
                         & 'set_arraysize_realpointer_2d')
            RETURN
          END IF
!        endif

        !
        ! 3. Fill the array with the optional value.
        !
        IF (PRESENT(fillvalue))then
          if(size_old_1 < size(array,1)) &
                  & array(size_old_1+1:size(array,1), 1:size(array,2)) = fillvalue
          if(size_old_2 < size(array,2)) &
                  & array(1:min(size_old_1+1,size(array,1)), size_old_2+1:size(array,2)) = fillvalue
        endif

      END IF size_check

!           PRINT *, 'Old array, already same size.'

    END IF new_or_old
        
  END SUBROUTINE set_arraysize_realpointer_2d
  


  ! ***************************************************************
  
  SUBROUTINE free_array_memory_realpointer(array)
    
    ! Description:
    ! Frees a memory of a real array, which is of type pointer. 
    !
    ! Language: ANSI Fortran 90
    !
    ! Author: Mika Salonoja, FMI
    ! 
    IMPLICIT NONE
    !
    ! Imported parameters with intent INOUT or POINTER:
    REAL, DIMENSION(:), POINTER :: array
    
    ! Local declarations: 
    INTEGER :: status

    IF (ASSOCIATED(array)) THEN
      !      PRINT *, 'deallocating...'
      DEALLOCATE(array, stat = status)
      IF (status /= 0) THEN
        CALL set_error('error while freeing memory', 'free_memory_real')
        RETURN
      END IF
    END IF
    
  END SUBROUTINE free_array_memory_realpointer


 
  ! ***************************************************************

  INTEGER FUNCTION fu_minimum_total_memory_usage()
    
    ! Description:
    ! Prints the total amount of memory consumed by the set_arraysize
    ! - subroutines during the run.
    !
    ! All units: SI
    !
    ! Language: ANSI Fortran 90
    !
    ! Author: Mika Salonoja, FMI
    ! 
    IMPLICIT NONE
    

    fu_minimum_total_memory_usage = total_memory_usage
  
  END FUNCTION fu_minimum_total_memory_usage




  ! ***************************************************************
  ! ***************************************************************
  ! 
  !
  !        Random numbers 
  !
  !
  ! ***************************************************************
  ! ***************************************************************
  REAL FUNCTION fu_random_number_center_rng(rng, center, width)
    
    ! Returns a random number that is distributed 
    ! uniformly in (center-width) .... (center+width).
    ! given generator is used 
    ! 
    IMPLICIT NONE
    REAL, INTENT(in) :: width, center
    type(rng_type), INTENT(inout) :: rng

    ! Local declarations:
    REAL :: x

    !!!Causes locks in multithread
    !!!CALL RANDOM_NUMBER(x) ! uniform between 0...1 
    x = rng_uniform(rng)
    
    fu_random_number_center_rng = (2.*(x - 0.5)*width) + center
    
  END FUNCTION fu_random_number_center_rng

  subroutine rng_init(self, i)
    ! Inspired by "A Parallel Monte Carlo Experiment"  http://jblevins.org/log/openmp
    type(rng_type), intent(inout) :: self
    integer, intent(in) :: i ! Some integer, usually a function of OMP loop index

    call rng_seed(self, 932117 + i)
  end subroutine rng_init

  ! Seeds the RNG using a single integer and a default seed vector.
  subroutine rng_seed(self, seed)
    ! Marsaglia and Zaman (1994)  from http://jblevins.org/log/openmp
    type(rng_type), intent(inout) :: self
    integer, intent(in) :: seed

    self%state(1) = seed
    self%state(2:ns) = default_seed(2:ns)
  end subroutine rng_seed

  ! Draws a uniform real number on [0,1].
  function rng_uniform(self) result(u)
    ! Marsaglia and Zaman (1994) rng from http://jblevins.org/log/openmp
    type(rng_type), intent(inout) :: self
    real :: u
    integer :: imz

    imz = self%state(1) - self%state(3)

    if (imz < 0) imz = imz + 2147483579

    self%state(1) = self%state(2)
    self%state(2) = self%state(3)
    self%state(3) = imz
    self%state(4) = 69069 * self%state(4) + 1013904243
    imz = imz + self%state(4)
    u = 0.5d0 + 0.23283064d-9 * imz
  end function rng_uniform

  

  REAL FUNCTION fu_random_number_center(center, width)
    
    ! Description:
    ! Returns a random number that is distributed around
    ! the given center. 
    !
    ! The given number is always distributed in range 
    ! (center-width) .... (center+width)
    !
    ! Language: ANSI Fortran 90
    !
    ! Author: Mika Salonoja, FMI
    ! 
    IMPLICIT NONE
    !
    ! The return value of this function:

    ! Imported parameters with intent(in):
    REAL, INTENT(in) :: width, center

    ! Local declarations:
    REAL :: x

    CALL RANDOM_NUMBER(x) ! uniform between 0...1
    
    fu_random_number_center = (2.*(x - 0.5)*width) + center
    
    
  END FUNCTION fu_random_number_center


  ! ***************************************************************

  REAL FUNCTION fu_random_number_boundaries(bottom, top)
    
    ! Description:
    ! Returns a random number that is distributed between given
    ! boundaries.
    !
    ! The given number is always distributed in range 
    ! bottom .... top.
    !
    ! Language: ANSI Fortran 90
    !
    ! Author: Mika Salonoja, FMI
    ! 
    IMPLICIT NONE
    !
    ! The return value of this function:

    ! Imported parameters with intent(in):
    REAL, INTENT(in) :: bottom, top

    ! Local declarations:
    REAL :: x

    CALL RANDOM_NUMBER(x) ! uniform between 0...1
    
    fu_random_number_boundaries = (x * (top - bottom)) + bottom
    
    
  END FUNCTION fu_random_number_boundaries

  !************************************************************************************

  subroutine random_normal(values)
    ! Get normally(0,1) distributed random numbers. Takes array because 2 number are
    ! generated at time.
    implicit none
    real, dimension(:), intent(out) :: values
    
    integer :: ii
    real :: dummy
    
    do ii = 1, size(values)-1, 2
      call box_muller(values(ii), values(ii+1))
    end do

    if (mod(size(values),2) > 0) then
      ! the last element remains if the size was odd or if size(values) == 1
      call box_muller(values(size(values)), dummy)
    end if

  contains
    
    subroutine box_muller(y1, y2)
      ! Generate two normal(0,1) random numbers using the default random number generator
      ! and the Box-Muller transformation.
      implicit none
      real, intent(out) :: y1, y2

      real :: w, rn1, rn2, x1, x2

      w = 2.0
      do while (w >= 1.0)
        ! sample the unit circle
        call random_number(rn1)
        x1 = 2*rn1 - 1.0
        call random_number(rn2)
        x2 = 2*rn2 - 1.0
        w = x1*x1 + x2*x2
      end do
      w = sqrt((-2.0 * log(w)) / w)
      y1 = x1*w
      y2 = x2*w

    end subroutine box_muller

  end subroutine random_normal


  ! ***************************************************************

  SUBROUTINE init_random_seed(const)
    !
    ! Initializes the random-number calculation. Only needs to be
    ! called once for each run. 
    !
    ! Method:
    ! A seed is given from the seconds of wallclock time.
    !
    IMPLICIT NONE
    integer, intent(in), optional :: const

    ! Imported parameters with intent INOUT or POINTER:

    ! Local declarations:
    INTEGER :: i, j
    integer, dimension(:), allocatable :: intAr

!    CALL RANDOM_SEED ( )  ! Processor reinitializes the seed                       
                           ! randomly from the date and time 
    CALL RANDOM_SEED (SIZE = i)  ! i is set to the size of the seed array 
    ALLOCATE (intAr(I))
    if (present(const)) then
      intAr(1:i) = 10*const
    else
      intAr(1:i) = 10
    end if
    
    CALL RANDOM_SEED (PUT = intAr(1:i)) ! Sets seed from array
    
    deallocate(intAr)
        
  END SUBROUTINE init_random_seed






  ! ***************************************************************

  SUBROUTINE next_line_from_input_file(input_file_unit, line, eof)

    ! Description:
    ! Returns next non-comment line fron an open input file (source
    ! term file, control parameter file or trajectory input file).
    ! Also if there's comment in the end of line, it is removed.
    ! All lines beginning with # or ! are considered to be comments.
    !
    IMPLICIT NONE

    ! Imported parameters with intent IN:
    INTEGER, INTENT(in) :: input_file_unit

    ! Imported parameters with intent OUT:
    CHARACTER (LEN=*), INTENT(out) ::  line
    LOGICAL, INTENT(out) :: eof


    ! Local declarations:
    LOGICAL :: file_open
    INTEGER :: status, i

    ! JV CODE

!!$    ! ------------------------------------------------
!!$    !
!!$    ! 1. Check file
!!$    !    ----------
!!$
!!$    !INQUIRE (unit = input_file_unit, opened = file_open)
!!$    !IF (.NOT.file_open) THEN
!!$    !  CALL set_error('file not open','next_line_from_input_file')
!!$    !  RETURN
!!$    !ELSE
!!$      eof = .false.
!!$    !END IF
!!$
!!$    ! ------------------------------------------------
!!$    !
!!$    ! 2. Read a string
!!$    !    -------------
!!$
!!$    DO
!!$      READ(unit = input_file_unit,&
!!$          & fmt = '(A)',&
!!$          & iostat = status) line
!!$      
!!$      IF (status /= 0) THEN
!!$        !        PRINT *, 'next_line_from_input_file status:', status
!!$        eof = .true.
!!$        line(:) = ' '
!!$        RETURN
!!$      END IF
!!$
!!$      ! Non-blanco line beginning with something else than # or !:
!!$      !line = adjustl(line)    ! May be, line starts not from the beginning
!!$      IF (line == '') then
!!$        !print *, 'empty:', line
!!$        CYCLE
!!$      end IF
!!$      i = index(line, '#')
!!$      if (i > 0) then
!!$        line(i:) = ' '
!!$      end if
!!$      
!!$      i = index(line, '!')
!!$      if (i > 0) then
!!$        if (i == 1) then
!!$          line = ''
!!$        else if (line(i-1:i-1) == '') then
!!$          line(i:) = ''
!!$        end if
!!$        
!!$      end if
!!$      if (line == '') then
!!$        cycle
!!$      else
!!$        exit
!!$      end if
!!$      !IF (line(1:1) /= '#' .and. line(1:1) /= '!') EXIT
!!$    END DO
!!$    line = adjustl(line)

    ! ORIGINAL CODE

    INQUIRE (unit = input_file_unit, opened = file_open)
    IF (.NOT.file_open) THEN
      CALL set_error('file not open','next_line_from_input_file')
      RETURN
    ELSE
      eof = .false.
    END IF

    !
    ! 2. Read a string
    !
    DO
      READ(unit = input_file_unit, fmt = '(A)', iostat = status) line

      IF (status /= 0) THEN
        !        PRINT *, 'next_line_from_input_file status:', status
        eof = .true.
        line(:) = ' '
        RETURN
      END IF

      ! Non-blanco line beginning with something else than # or !:
      line = adjustl(line)    ! May be, line starts not from the beginning
      IF (line == ' ') CYCLE
      IF (line(1:1) /= '#' .and. line(1:1) /= '!') EXIT
    END DO

    !
    ! Get rid of tabs
    !
    DO i = 1, LEN_trim(line)
      IF (line(i:i) == achar(09))then
!        call msg('line with tabs:' + line)
        line(i:i) = ' '
!        call msg('line without tabs:' + line)
      endif
    END DO

    ! JV modification

!!$    i = index(line, '#')
!!$    if (i > 0) then
!!$      line(i:) = ''
!!$    end if
!!$    i = index(line, ' !')
!!$    if (i > 0) then
!!$      line(i:) = ''
!!$    end if

    ! original

    DO i = 2, LEN(line)
      IF (line(i:i) == '#' .or. line(i-1:i) == ' !') THEN 
        line(i:) = ' '
        EXIT
      END IF
    END DO

  END SUBROUTINE next_line_from_input_file



  ! ***************************************************************

  INTEGER FUNCTION fu_next_free_unit()
    
    ! Description:
    ! Returns a unit number that is not connected to any opened file.
    ! Thread-safe
    !
    ! 
    IMPLICIT NONE
    !
    ! Local declarations:
    LOGICAL :: fopen
    integer :: ithread, nthreads
    !$ if (omp_in_parallel()) then
    !$    ithread=omp_get_thread_num()
    !$    nthreads=omp_get_num_threads()
    !$ else
        ithread=0
        nthreads=1
    !$ endif

    DO fu_next_free_unit = 7+ithread, 1000, nthreads
      INQUIRE(unit = fu_next_free_unit, opened = fopen)
      IF (.NOT.fopen) EXIT
    END DO

!    call msg('Next free_unit = ', fu_next_free_unit)

  END FUNCTION fu_next_free_unit

  
  !******************************************************************


  FUNCTION fu_merge_int_arrays(ar_in, ar_out, ifHolesAllowed) RESULT(iCountAdded)
  !
  ! Simply merges two arrays ar_in and ar_out putting the output to 
  ! ar_out. The only rule - no repetition of elements. Criteria for 
  ! "non-existing element" is int_missing.
  ! Both arrays can have holes filled with int_missing
  ! Obs. : No sorting is made.
  !
  ! Return value - number of added elements to ar_out
  INTEGER iCountAdded

  ! Import with intention IN:
  INTEGER, DIMENSION(:), INTENT(in) :: ar_in
  logical, intent(in) :: ifHolesAllowed

  ! Import with intent INOUT
  INTEGER, DIMENSION(:), INTENT(inout) :: ar_out

  ! Local declarations:
  INTEGER :: iTmp, iCount

  iCount = 1

  ! Packing the ar_out - removing holes.
  DO iTmp=1,SIZE(ar_out)
    IF(ar_out(iTmp) /= int_missing)THEN
      ar_out(iCount)=ar_out(iTmp)
      iCount=iCount+1
    else
      if(.not. ifHolesAllowed)exit  ! int_missing breaks the cycle
    END IF
  END DO

  iCountAdded = iCount

  ! Adding elements of ar_in to the end of ar_out if they are not already in it
  DO iTmp=1,SIZE(ar_in)
    IF(ar_in(iTmp) /= int_missing)THEN
      IF(.not. ANY(ar_out(1:iCount-1) == ar_in(iTmp))) THEN
        ar_out(iCount)=ar_in(iTmp)
        iCount=iCount+1
        IF(iCount > SIZE(ar_out))THEN
          CALL set_error('too small array size','merge_int_arrays')
          RETURN
        END IF
      END IF
    else
      if(.not. ifHolesAllowed)exit  ! int_missing breaks the cycle
    END IF
  END DO

  iCountAdded = iCount-iCountAdded

  END FUNCTION fu_merge_int_arrays

  !***************************************************************

  integer function fu_merge_int_to_array(int, ar)
    ! 
    ! Checks if this integer is new and adds at the end if it is.
    ! Returns the new size of the array
    !
    implicit none

    ! Imported parameters
    integer, intent(in) :: int
    integer, dimension(:), intent(inout) :: ar

    ! Local variables
    integer :: i
    logical :: ifUnique

    if(int == int_missing) return
    !
    ! Find free place in the array
    !
    fu_merge_int_to_array = -1
    ifUnique = .true.
    do i=1, size(ar)
      if(ar(i) == int_missing) then
        fu_merge_int_to_array = i-1
        exit
      else
        if(ar(i) == int) ifUnique = .false.
      endif
    end do
    if(fu_merge_int_to_array == -1)then ! Array is full
      call msg("Array size:",  size(ar))
      call set_error('Too small receptor array','fu_merge_int_to_array')
      return
    endif
    !
    ! Add the integer in case it is unique.
    !
    if(ifUnique) then
      ar(i) = int
      fu_merge_int_to_array = i
    end if

  end function fu_merge_int_to_array


  !***************************************************************

  integer function fu_merge_int_2_arr_with_slave(int, ar, int_slave, ar_slave)
    ! 
    ! Checks if this integer is new and adds at the end if it is.
    ! Returns the new size of the array
    !
    implicit none

    ! Imported parameters
    integer, intent(in) :: int, int_slave
    integer, dimension(:), intent(inout) :: ar, ar_slave

    ! Local variables
    integer :: i
    logical :: ifUnique

    if(int == int_missing) return
    !
    ! Find free place in the array
    !
    fu_merge_int_2_arr_with_slave = -1
    ifUnique = .true.
    do i=1, size(ar)
      if(ar(i) == int_missing) then
        fu_merge_int_2_arr_with_slave = i-1 ! found the place
        exit
      else
        if(ar(i) == int)then   ! the value is already in array
          ar_slave(i) = max(ar_slave(i), int_slave)
          ifUnique = .false.
        endif
      endif
    end do
    if(fu_merge_int_2_arr_with_slave == -1)then ! Array is full
      call set_error('Too small receptor array','fu_merge_int_2_arr_with_slave')
      return
    endif
    !
    ! Add the integer in case it is unique.
    !
    if(ifUnique) then
      ar(i) = int
      fu_merge_int_2_arr_with_slave = i
      ar_slave(i) = int_slave
    end if

  end function fu_merge_int_2_arr_with_slave

  !*************************************************************************

  logical function fu_sort_real_array(ar, order, ifHolesAllowed)
    !
    ! Just sorts an array ar in an ascending - descending order
    ! The array may have holes filled with real_missing, they will be
    ! eliminated
    ! Method: bubble
    !
    implicit none

    ! Imoprted parameters
    !
    real, dimension(:), intent(inout) :: ar
    integer, intent(in) :: order
    logical, intent(in) :: ifHolesAllowed
    !
    ! Local variables
    integer :: i, iPrev, iCount
    real :: aTmp

    ! Stupidity check
    !
    fu_sort_real_array = .false.
    if(all(ar(1:size(ar)) == real_missing))return
    if((ar(1) .eps. real_missing) .and. (.not. ifHolesAllowed))return ! if holes not allowed return: empty array

    iCount = 1 ! Number of exchanges
    do while(iCount > 0)
      iCount = 0
      iPrev = 1
      do i=2,size(ar)
        if(ar(i).eps.real_missing)then
          if(ifHolesAllowed)cycle     ! if holes allowed, skip the element, exit otherwise
          exit
        endif
        aTmp = ar(i)
        ar(i) = real_missing
        if(order == ascending)then  ! Ascending order
          if(ar(iPrev) > aTmp)then
            ar(iPrev+1) = ar(iPrev)
            ar(iPrev) = aTmp
            iCount = iCount + 1
          else
            ar(iPrev+1) = aTmp
          endif
        elseif(order == descending)then ! Descending order
          if(ar(iPrev) < aTmp)then
            ar(iPrev+1) = ar(iPrev)
            ar(iPrev) = aTmp
            iCount = iCount + 1
          else
            ar(iPrev+1) = aTmp
          endif
        else
          call set_error('Unknown sorting order','fu_sort_real_array')
          return
        endif
        iPrev = iPrev+1
      end do
    end do

    fu_sort_real_array = .true.

  end function fu_sort_real_array


  !***************************************************************************

  REAL FUNCTION fu_2d_interp_to_coord (grid_data, x,y, nx,ny, method, iOutside) result(value)
    !
    ! Interpolates the field given in a grid to given co-ordinates. Most fast
    ! because does not use any function calls and any derived types. 
    ! No grids or projections are known here
    !
    IMPLICIT NONE
    !
    ! Imported parametersn with intent(in):
    real, INTENT(in) :: x, y ! Co-ordinates.
    INTEGER, INTENT(in) :: nx, ny ! grid dimensions
    INTEGER, INTENT(in) :: method ! of interpolation
    integer, intent(in) :: iOutside   ! what to do, if data out of grid (notAllowed, nearestPoint...)
    REAL, DIMENSION(:), intent(in) :: grid_data

    ! Local declarations:
    INTEGER :: southwest_corner, loop, xloop, yloop, ii
    INTEGER, DIMENSION(16)::corner_pointers 
    REAL :: zx, zy 
    REAL, DIMENSION(16)::ww
    REAL:: zxy1, zxy2, zxy3, zxy4, zmin
    REAL:: zx1,zx2,zx3,zx4,zy1,zy2,zy3,zy4     
    INTEGER:: ilat, ilon
    INTEGER :: method_local
    REAL, DIMENSION(4) :: distances
    INTEGER, DIMENSION(1) :: point
    integer, save :: iCount=0


    IF (method == cubic) THEN
      ! If cubic method and position too close to the border of the
      ! grid, use linear instead:
      IF((x<2.).or.(y<2.).or.(x>(nx-1.)).or.(y>(ny-1.))) THEN
        method_local = linear
      ELSE
        method_local = cubic
      END IF
    ELSE
      method_local = method
    END IF

    !
    ! If, due to shifts in the grids, e.g., x=nx+something_small, we must not use
    ! 4-corner interpolation but 2-corner only.
    ! In the most extreme case, we may have just one corner point - e.g. for (0.5,0.5)
    !
    IF((x<=1.).or.(y<=1.).or.(x>=nx).or.(y>=ny)) THEN

      if(iOutside == handleGlobalGrid)then
        call set_error('handleGlobalGrid is not yet supproted','fu_2d_interp_to_coord')
        return
      endif

      
      ! Outside the grid: depending on iOutside set value zero, missing, set error or continue to 
      ! nearest point interpolation

      !!!FIXME Dirty hack Should be cured by double-precision grids 
      !!! Might also hide some severe error
      if((x<0.35).or.(y<0.35).or.(x>nx+0.65).or.(y>ny+0.65)) then
        if(iOutside == setZero)then
          value = 0.
          return
        elseif(iOutside == setMissVal)then
          value = real_missing
          return
        elseif(iOutside == notAllowed)then
          call msg_warning('Out-of-array interpolation not allowed','fu_2d_interp_to_coord')
          call msg('nx, X-co-ordinate:', nx,x)
          call msg('ny, Y-co-ordinate:', ny,y)
          call set_error('Out-of-array interpoaltion not allowed','fu_2d_interp_to_coord')
          value = 0.  ! real_missing is too dangerous in case of wind, etc. Error is enough
          return
        endif
      endif
      ! Nearest point interpolation
      IF(x<=1.)THEN
        !
        ! Near the west border: dependence on x has vanished, left corners disappeared
        !
        if(y<=1)then
          value = grid_data(1)
        elseif(y>=ny)then
          value = grid_data((ny-1)*nx)
        else
          corner_pointers (1) = 1 
          corner_pointers (2) = nx + 1
          southwest_corner = ((int(y) - 1)*nx)
          if(southwest_corner + corner_pointers(2) > nx*ny .or. southwest_corner < 0)then
            call msg_warning('Out-of-array interpolation request 2','fu_2d_interp_to_coord')
            call msg('X-co-ordinate, nx:', nx,x)
            call msg('Y-co-ordinate, ny:', ny,y)
            call set_error('Out-of-array interpoaltion request','fu_2d_interp_to_coord')
            return
          endif
          zy = y - REAL(INT(y))
          ww(1) = (1. - zy)
          ww(2) = zy
          DO loop = 1, 2
            if(grid_data(southwest_corner + corner_pointers(loop)) .eps. real_missing) ww(loop) = 0.0
          end do
          if(sum(ww(1:2)).eps.0.0)then
            value = real_missing
            return
          else
            ww(1:2) = ww(1:2)/sum(ww(1:2))
            value = 0.
            DO loop = 1, 2
              value = value + ww(loop) * grid_data(southwest_corner + corner_pointers(loop))
            END DO
          endif
        endif

      elseIF(y<=1.)THEN
        !
        ! Near the southern border: dependence on y has vanished, south corners disappeared
        !
        if(x<=1)then
          value = grid_data(1)
        elseif(x>=nx)then
          value = grid_data(nx)
        else
          corner_pointers(1) = 0 
          corner_pointers(2) = 1
          southwest_corner = INT(x)
          if(southwest_corner + corner_pointers(2) > nx*ny .or. southwest_corner < 0)then
            call msg_warning('Out-of-array interpolation request 3','fu_2d_interp_to_coord')
            call msg('X-co-ordinate, nx:', nx,x)
            call msg('Y-co-ordinate, ny:', ny,y)
            call set_error('Out-of-array interpoaltion request','fu_2d_interp_to_coord')
            return
          endif
          zx = x - REAL(INT(x))
          ww(1) = 1. - zx
          ww(2) = zx
          DO loop = 1, 2
            if(grid_data(southwest_corner + corner_pointers(loop)) .eps. real_missing) ww(loop) = 0.0
          end do
          if(sum(ww(1:2)).eps.0.0)then
            value = real_missing
            return
          else
            ww(1:2) = ww(1:2)/sum(ww(1:2))
            value = 0.
            DO loop = 1, 2
              value = value + ww(loop) * grid_data(southwest_corner + corner_pointers(loop))
            END DO
          endif
        endif

      elseIF(x>=nx)THEN
        !
        ! Near the eastern border: dependence on x has vanished, east corners disappeared
        !
        if(y<=1)then
          value = grid_data(nx)
        elseif(y>=ny)then
          value = grid_data(nx*ny)
        else
          corner_pointers(1) = 0 
          corner_pointers(2) = nx 
          southwest_corner = INT(y)*nx
          if(southwest_corner + corner_pointers(2) > nx*ny .or. southwest_corner == 0)then
            call msg_warning('Out-of-array interpolation request 4','fu_2d_interp_to_coord')
            call msg('X-co-ordinate, nx:', nx,x)
            call msg('Y-co-ordinate, ny:', ny,y)
            call set_error('Out-of-array interpoaltion request','fu_2d_interp_to_coord')
            return
          endif
          zy = y - REAL(INT(y))
          ww(1) = 1. - zy
          ww(2) = zy
          DO loop = 1, 2
            if(grid_data(southwest_corner + corner_pointers(loop)) .eps. real_missing) ww(loop) = 0.0
          end do
          if(sum(ww(1:2)).eps.0.0)then
            value = real_missing
            return
          else
            ww(1:2) = ww(1:2)/sum(ww(1:2))
            value = 0.
            DO loop = 1, 2
              value = value + ww(loop) * grid_data(southwest_corner + corner_pointers(loop))
            END DO
          endif
        endif

      else
        !
        ! Near the northern border: dependence on y has vanished, north corners disappeared
        !
        if(x<=1)then
          value = grid_data((ny-1)*nx+1)
        elseif(x>=nx)then
          value = grid_data(nx*ny)
        else
          corner_pointers(1) = 0 
          corner_pointers(2) = 1 
          southwest_corner = ((ny - 1)*nx) + INT(x)
          if(southwest_corner + corner_pointers (2) > nx*ny .or. southwest_corner == 0)then
            call msg_warning('Out-of-array interpolation request 5','fu_2d_interp_to_coord')
            call msg('X-co-ordinate, nx:', nx,x)
            call msg('Y-co-ordinate, ny:', ny,y)
            call set_error('Out-of-array interpoaltion request','fu_2d_interp_to_coord')
            return
          endif
          zx = x - REAL(INT(x))
          ww(1) = 1. - zx
          ww(2) = zx
          DO loop = 1, 2
            if(grid_data(southwest_corner + corner_pointers(loop)) .eps. real_missing) ww(loop) = 0.0
          end do
          if(sum(ww(1:2)).eps.0.0)then
            value = real_missing
            return
          else
            ww(1:2) = ww(1:2)/sum(ww(1:2))
            value = 0.
            DO loop = 1, 2
              value = value + ww(loop) * grid_data(southwest_corner + corner_pointers(loop))
            END DO
          endif

        endif

      endif ! types of near-border cases

      return ! near-border interpolation is done

    endif ! if near-border interpolation
    
    !-------------------------------------------------------------------
    !
    ! General case
    !
    method_of_interpolation: SELECT CASE (method_local)

      ! -----------------------------------------------
      ! 
      ! 2. Linear interpolation.
      !    --------------------

      CASE (linear)

      corner_pointers (1) = 0 
      corner_pointers (2) = 1 
      corner_pointers (3) = nx 
      corner_pointers (4) = nx + 1

      southwest_corner = ((INT(y) - 1)*nx) + INT(x)

      if(southwest_corner + corner_pointers (4) > nx*ny .or. southwest_corner == 0)then
        call msg_warning('Out-of-array interpolation request 6','fu_2d_interp_to_coord')
        call msg('X-co-ordinate, nx:', nx,x)
        call msg('Y-co-ordinate, ny:', ny,y)
        call set_error('Out-of-array interpoaltion request','fu_2d_interp_to_coord')
        return
      endif

      zx = x - REAL(INT(x))
      zy = y - REAL(INT(y))

      ww(1) = (1. - zx) * (1. - zy)
      ww(2) = zx * (1. - zy)
      ww(3) = (1. - zx) * zy
      ww(4) = zx * zy

      DO loop = 1, 4
        if(grid_data(southwest_corner + corner_pointers(loop)) .eps. real_missing) ww(loop) = 0.0
      end do
      if(sum(ww(1:4)).eps.0.0)then
        value = real_missing
        return
      else
        ww(1:4) = ww(1:4)/sum(ww(1:4))
        value = 0.
        DO loop = 1, 4
          value = value + ww(loop) * grid_data(southwest_corner + corner_pointers(loop))
        END DO
      endif

      ! -----------------------------------------------
      ! 
      ! 3. Bi-cubic interpolation.
      !    --------------------

      CASE (cubic)

      ii = 0
      ilat = -nx !      

      DO yloop = 1, 4 
        ilon = -1 !
        DO xloop = 1, 4 
          ii = ii + 1
          corner_pointers(ii) = ilon + ilat
          ilon = ilon + 1
        END DO
        ilat = ilat + nx
      END DO

      southwest_corner = ((INT(y) - 1)*nx) + INT(x)

      if(southwest_corner + corner_pointers(ii) > nx*ny .or. southwest_corner == 0)then
        call msg_warning('Out-of-array interpolation request 7','fu_2d_interp_to_coord')
        call msg('X-co-ordinate, nx:', nx,x)
        call msg('Y-co-ordinate, ny:', ny,y)
        call set_error('Out-of-array interpoaltion request','fu_2d_interp_to_coord')
        return
      endif

      zx = x - REAL(INT(x))
      zy = y - REAL(INT(y))

      zx1 = ((-0.5*zx+1.0)*zx-0.5)*zx 
      zx2 = (( 1.5*zx-2.5)*zx    )*zx+1 
      zx3 = ((-1.5*zx+2.0)*zx+0.5)*zx 
      zx4 = (( 0.5*zx-0.5)*zx    )*zx 
      zy1 = ((-0.5*zy+1.0)*zy-0.5)*zy 
      zy2 = (( 1.5*zy-2.5)*zy    )*zy+1 
      zy3 = ((-1.5*zy+2.0)*zy+0.5)*zy 
      zy4 = (( 0.5*zy-0.5)*zy    )*zy

      ww( 1) = zx1*zy1 
      ww( 2) = zx2*zy1
      ww( 3) = zx3*zy1  
      ww( 4) = zx4*zy1 
      ww( 5) = zx1*zy2  
      ww( 6) = zx2*zy2  
      ww( 7) = zx3*zy2  
      ww( 8) = zx4*zy2  
      ww( 9) = zx1*zy3  
      ww(10) = zx2*zy3  
      ww(11) = zx3*zy3 
      ww(12) = zx4*zy3 
      ww(13) = zx1*zy4 
      ww(14) = zx2*zy4   
      ww(15) = zx3*zy4   
      ww(16) = zx4*zy4   

      DO loop = 1, 16
        if(grid_data(southwest_corner + corner_pointers(loop)) .eps. real_missing) ww(loop) = 0.0
      end do
      if(sum(ww(1:16)).eps.0.0)then
        value = real_missing
        return
      else
        ww(1:16) = ww(1:16)/sum(ww(1:16))
        value = 0.
        DO loop = 1, 16
          value = value + ww(loop) * grid_data(southwest_corner + corner_pointers(loop))
        END DO
      endif

      ! -----------------------------------------------
      ! 
      ! 4. Take the nearest point
      !    --------------------

      CASE (nearest_point)

      corner_pointers (1) = 0 
      corner_pointers (2) = 1 
      corner_pointers (3) = nx 
      corner_pointers (4) = nx + 1

      southwest_corner = ((INT(y) - 1)*nx) + INT(x)

      if(southwest_corner + corner_pointers (4) > nx*ny .or. southwest_corner == 0)then
        call msg_warning('Out-of-array interpolation request 8','fu_2d_interp_to_coord')
        call msg('X-co-ordinate, nx:', nx,x)
        call msg('Y-co-ordinate, ny:', ny,y)
        call set_error('Out-of-array interpoaltion request','fu_2d_interp_to_coord')
        return
      endif

      zx = x - REAL(INT(x))
      zy = y - REAL(INT(y))

      distances(1) = zx*zx + zy*zy
      distances(2) = (1.-zx)*(1.-zx) + zy*zy
      distances(3) = zx*zx + (1.-zy)*(1.-zy)
      distances(4) = (1.-zx)*(1.-zx) + (1.-zy)*(1.-zy)

      point = MINLOC(distances)

      value = grid_data(southwest_corner + corner_pointers(point(1)))



    CASE default

      CALL set_error('unknown interpolation method','fu_2d_interp_to_coord')


    END SELECT method_of_interpolation

  END FUNCTION fu_2d_interp_to_coord


  !*************************************************************************

  subroutine copy_text_file(f_from, f_to)
    !
    ! Copies an opened file to another opened file
    !
    implicit none

    ! Imported parameters
    integer, intent(in) :: f_from, f_to

    ! Local variables
    logical :: eof, file_open
    character(len=fnlen) :: line
    integer :: status
    character(len=*), parameter :: sub_name = 'copy_text_file'

    !
    ! Stupidity check first
    !
    INQUIRE (unit = f_from, opened = file_open)
    IF (.NOT.file_open) THEN
      CALL set_error('Input file not open', sub_name)
      RETURN
    ELSE
      rewind(f_from)
      eof = .false.
    END IF

    INQUIRE (unit = f_to, opened = file_open)
    IF (.NOT.file_open) THEN
      CALL set_error('Output file not open', sub_name)
      RETURN
    END IF

    ! ------------------------------------------------
    !
    ! 2. Read a string
    !    -------------
    eof = .false.
    DO while(.not.eof)
      line = ''  ! just in case, if the read fails but does not empties the line
      READ(unit = f_from, fmt = '(A)', iostat = status) line

      IF (status /= 0) THEN
        !        PRINT *, 'next_line_from_input_file status:', status
        eof = .true.
!        RETURN
      END IF

      write(unit = f_to, fmt = '(A)', iostat = status) trim(line)
      IF (status /= 0) THEN
        !        PRINT *, 'next_line_from_input_file status:', status
        call msg(line)
        call set_error('Failed to write the above line', sub_name)
        RETURN
      END IF

    END DO

  end subroutine copy_text_file


  !***************************************************************************

  character function fu_u_case(chA)
    !
    ! Converts the input character to the upper case analog
    !
    implicit none

    ! Imported parameter
    character, intent(in) :: chA

    select case(chA)
      case('a','A')
        fu_u_case = 'A'
      case('b','B')
        fu_u_case = 'B'
      case('c','C')
        fu_u_case = 'C'
      case('d','D')
        fu_u_case = 'D'
      case('e','E')
        fu_u_case = 'E'
      case('f','F')
        fu_u_case = 'F'
      case('g','G')
        fu_u_case = 'G'
      case('h','H')
        fu_u_case = 'H'
      case('i','I')
        fu_u_case = 'I'
      case('j','J')
        fu_u_case = 'J'
      case('k','K')
        fu_u_case = 'K'
      case('l','L')
        fu_u_case = 'L'
      case('m','M')
        fu_u_case = 'M'
      case('n','N')
        fu_u_case = 'N'
      case('o','O')
        fu_u_case = 'O'
      case('p','P')
        fu_u_case = 'P'
      case('q','Q')
        fu_u_case = 'Q'
      case('r','R')
        fu_u_case = 'R'
      case('s','S')
        fu_u_case = 'S'
      case('t','T')
        fu_u_case = 'T'
      case('u','U')
        fu_u_case = 'U'
      case('v','V')
        fu_u_case = 'V'
      case('w','W')
        fu_u_case = 'W'
      case('x','X')
        fu_u_case = 'X'
      case('y','Y')
        fu_u_case = 'Y'
      case('z','Z')
        fu_u_case = 'Z'
      case default
        fu_u_case = chA
    end select

  end function fu_u_case


  !***************************************************************************

  character function fu_l_case(chA)
    !
    ! Converts the input character to the upper case analog
    !
    implicit none

    ! Imported parameter
    character, intent(in) :: chA

    select case(chA)
      case('a','A')
        fu_l_case = 'a'
      case('b','B')
        fu_l_case = 'b'
      case('c','C')
        fu_l_case = 'c'
      case('d','D')
        fu_l_case = 'd'
      case('e','E')
        fu_l_case = 'e'
      case('f','F')
        fu_l_case = 'f'
      case('g','G')
        fu_l_case = 'g'
      case('h','H')
        fu_l_case = 'h'
      case('i','I')
        fu_l_case = 'i'
      case('j','J')
        fu_l_case = 'j'
      case('k','K')
        fu_l_case = 'k'
      case('l','L')
        fu_l_case = 'l'
      case('m','M')
        fu_l_case = 'm'
      case('n','N')
        fu_l_case = 'n'
      case('o','O')
        fu_l_case = 'o'
      case('p','P')
        fu_l_case = 'p'
      case('q','Q')
        fu_l_case = 'q'
      case('r','R')
        fu_l_case = 'r'
      case('s','S')
        fu_l_case = 's'
      case('t','T')
        fu_l_case = 't'
      case('u','U')
        fu_l_case = 'u'
      case('v','V')
        fu_l_case = 'v'
      case('w','W')
        fu_l_case = 'w'
      case('x','X')
        fu_l_case = 'x'
      case('y','Y')
        fu_l_case = 'y'
      case('z','Z')
        fu_l_case = 'z'
      case default
        fu_l_case = chA
    end select

  end function fu_l_case


  !***************************************************************************

  subroutine str_u_case(strA)
    !
    ! Converts the input character to the upper case analog
    !
    implicit none

    ! Imported parameter
    character(len=*), intent(inout) :: strA

    ! local variable
    integer :: i

    do i=1,len_trim(strA)
      strA(i:i) = fu_u_case(strA(i:i))
    enddo

  end subroutine str_u_case


  !***************************************************************************

  subroutine str_l_case(strA)
    !
    ! Converts the input character to the upper case analog
    !
    implicit none

    ! Imported parameter
    character(len=*), intent(inout) :: strA

    ! local variable
    integer :: i

    do i=1,len_trim(strA)
      strA(i:i) = fu_l_case(strA(i:i))
    enddo

  end subroutine str_l_case


  !***************************************************************************

  function fu_str_u_case(strA) result(ucase_A)
    !
    ! Converts the input character string to the upper case analog
    !
    implicit none
    character(len=*), intent(in) :: strA
    character(len=len(strA)) :: ucase_A
    integer :: ci, i

    ucase_A = ''
    do i = 1, len_trim(strA)
      ci = iachar(strA(i:i))
      if (ci >= ai .and. ci <= zi) then
        ucase_A(i:i) = achar(ci - diff)
      else
        ucase_A(i:i) = strA(i:i)
      end if
    end do

  end function fu_str_u_case

    
  !***************************************************************************

  function fu_str_l_case(strA) result(lcase_A)
    !
    ! Converts the input character to the upper case analog
    !
    implicit none
    ! Imported parameter
    character(len=*), intent(in) :: strA

    ! Return value
    character(len=len(strA)) :: lcase_A

    ! local variable
    integer :: i, ci

    lcase_A=''
    do i = 1, len_trim(strA)
      ci = iachar(strA(i:i))
      if (ci >= ai-diff .and. ci <= zi-diff) then
        lcase_A(i:i) = achar(ci + diff)
      else
        lcase_A(i:i) = strA(i:i)
      end if
    end do

  end function fu_str_l_case


  !***************************************************************************

  subroutine compress_int_array(arInt, iMiss, iNbrOfElem)
    !
    ! Removes holes from the array. iMiss defines the hole. 
    !
    implicit none

    ! Imported parameters
    integer, dimension(:), intent(inout) :: arInt
    integer, intent(in) :: iMiss
    integer, intent(inout), optional :: iNbrOfElem ! Nbr of elements before / after compression

    ! Local variables
    integer :: i, j, iMax

    !
    ! First, find iMax - the last meaningfull element. Just to speed-up
    !
    if(present(iNbrOfElem))then
      iMax = iNbrOfElem
    else
      iMax=1
      do i=1, size(arInt)
        if(arInt(i) /= iMiss) iMax = i
      enddo
    end if
    !
    ! Go throught the array as many times as needed
    !
    i=1
    do while (i <= iMax)

      if (arInt(i) == iMiss) then
        do j = i, min(iMax,size(arInt)-1)
          arInt(j) = arInt(j+1)
        end do
        arInt(iMax) = iMiss ! j+1 is dangerous if iMax == size(arInt)
        iMax = iMax -1
        i=i-1 ! May be, the i+1 field was also missing
      end if
      i=i+1

    end do ! Through array

    if(present(iNbrOfElem)) iNbrOfElem = iMax

  end subroutine compress_int_array 


  !*******************************************************************

  integer function fu_nbrOfWords(str)
    !
    ! Counts the number of words in the strings using space as 
    ! a separator
    !
    implicit none

    ! Imported parameters
    character(len=*), intent(in) :: str

    ! Local variables
    integer :: iPos

    fu_nbrOfWords = 0
    iPos = 1

    if(len_trim(str) < 1)return ! Empty string
    do while(str(iPos:iPos) == ' ')
      iPos = iPos + 1
    end do

    do while(iPos <= len_trim(str))
      do while(str(iPos:iPos) /= ' ' .and. iPos <= len_trim(str))
        iPos = iPos + 1
      end do
      fu_nbrOfWords = fu_nbrOfWords + 1
      do while(str(iPos:iPos) == ' ' .and. iPos <= len_trim(str))
        iPos = iPos + 1
      end do
    end do

  end function fu_nbrOfWords



  !*******************************************************************
  !*******************************************************************
  !
  !
  !    UNIT CONVERSION ROUTINES
  !
  !
  !*******************************************************************
  !*******************************************************************


  !***********************************************************************

  logical function fu_if_digit(char)
    !
    ! Checsk whether the given character is a digit
    !
    implicit none

    ! Imported parameters
    character(len=1), intent(in) :: char

    select case(char)
      case('0','1','2','3','4','5','6','7','8','9')
        fu_if_digit = .true.
      case default
        fu_if_digit = .false.
    end select

  end function fu_if_digit


  !***********************************************************************

  real function fu_solar_zenith_angle_cos_par(lon_deg, lat_deg, julian_day, &
                                            & year, mon, day, hour, min, sec)
    !
    ! Computes the solar zenith angle from basic parameters and UTC time.
    ! To account for the local time, we use longitude
    !
    implicit none
    
    ! Imported parameters
    real, intent(in) :: lon_deg, lat_deg
    integer, intent(in) :: julian_day, year, mon, day, hour, min, sec

    ! Local variables
    real :: d, Eqt, Tsolar, w

    !
    ! Declination of the sun
    !
    d = 23.45 * pi / 180. * sin(2. * pi * (284. + julian_day) / 365.)

    !
    ! Equation of time in minutes
    ! http://solardat.uoregon.edu/SolarRadiationBasics.html
    !
    if(julian_day < 106)then
      Eqt = -14.2 * sin(pi * (julian_day + 7.) / 111.)
    elseif(julian_day < 166)then
      Eqt = 4.0 * sin(pi * (julian_day - 106.) / 59.)
    elseif(julian_day < 246)then
      Eqt = -6.5 * sin(pi * (julian_day - 166.) / 80.)
    else
      Eqt = 16.4 * sin(pi * (julian_day - 247.) / 113.)
    endif

    !
    ! Solar time in hours. Longsm -s the longitude of the "standard meridian" for the given time zone,
    ! while longLocal is the actual longitude needed. The differrence is then less than 1 hour.
    ! If we count from Greenwich, Longsm=0, of course
    !
!    Tsolar = Tlocal + Eqt / 60. + (Longsm - Longlocal) / 15
    Tsolar = hour + (min+Eqt) / 60. + sec / 3600. + lon_deg / 15.
    !
    ! Hour angle is:
    !
    w = pi * (12. - Tsolar) / 12.  ! in radians
    !
    ! Cosine of zenith angle
    !
    fu_solar_zenith_angle_cos_par = sin(lat_deg * degrees_to_radians) * sin(d) + &
                                  & cos(lat_deg * degrees_to_radians) * cos(d) * cos(w)

  end function fu_solar_zenith_angle_cos_par

  !***************************************************************************************************


  real function fu_day_length_hrs(fLat_deg, julian_day)
    !
    ! COmputes the day length for the given place at the given Julian day
    ! uses http://solardat.uoregon.edu/SolarRadiationBasics.html
    !
    implicit none
    
    ! Imported paramters
    real, intent(in) :: fLat_deg
    integer, intent(in) :: julian_day
    
    ! Local variables
    real :: d, fTmp
    !
    ! Declination of the sun
    !
    d = 23.45 * pi / 180. * sin(2. * pi * (284. + julian_day) / 365.)
    !
    ! the sunrise / sunset hour angle: positive is sunrise, negative is sunset
    !
    ! Equation of time, minutes - see above
    ! hour angle vs solar time: w = pi * (12 - Tsolar) / 12
    ! Tsolar_sunrise_unset = 12 -_+ w * 12.0 / pi
    ! Solar time vs local time: Tsolar = Tlocal + Eqt / 60 + (Longlocal - Longsm) / 15
    ! daylen = T_sunset - T_sunrise
    ! Now cancel the stupid parts:
    !
    fTmp = tan(fLat_deg * degrees_to_radians) * tan(d)
    ! polar day/night?
    if(abs(fTmp) < 1)then
      fu_day_length_hrs = 2. * acos(-fTmp)  * 12.0 / pi
    else
      fu_day_length_hrs = -1  ! signaling polar day/night
    endif

  end function fu_day_length_hrs

  
  !***************************************************************************************************

  integer function fu_get_day_from_daylen(fLat_deg, fDayLength, ifSearchFromLongToShort)
    !
    ! Searches for a day with the day length close to the one requested
    ! Works for both hemispheres and can search for the first long- and short-length day
    !
    implicit none
    
    ! Imported parameters
    real, intent(in) :: fLat_deg, fDayLength
    logical, intent(in) :: ifSearchFromLongToShort
    
    ! Local variables
    integer :: iDay, iDayStart, iDayEnd  !, iDayTmp, iDayStartTmp
    
    if(ifSearchFromLongToShort)then
      ! Start from the longest day towards short one
      if(fLat_deg > 0)then
        iDayStart = 174        ! Northern Hemisphere, start from 22 June
      else
        iDayStart = -9        ! Southern Hemisphere, start from 22 December
      endif ! hemisphere
    else
      ! Start from the shortest day towards the long one
      if(fLat_deg > 0)then
        iDayStart = -9        ! Northern Hemisphere, start from 22 December
      else
        iDayStart = 174        ! Southern Hemisphere, start from 22 June
      endif ! hemisphere
    endif ! search from long to short or from short to long
    
    iDayEnd = iDayStart + 182
    !
    ! If start/end are from polar day/night, have to shorten the range
    !
    do while(fu_day_length_hrs(fLat_deg, iDayStart) < 0)
      iDayStart = iDayStart + 1
      if(fu_fails(iDayStart < iDayEnd,'Cannot find iDayStart','fu_get_day_from_daylen'))return
    end do
    do while(fu_day_length_hrs(fLat_deg, iDayEnd) < 0)
      iDayEnd = iDayEnd - 1
      if(fu_fails(iDayStart < iDayEnd,'Cannot find iDayEnd','fu_get_day_from_daylen'))return
    end do
    !
    ! Are we inside the possible range?
    ! Remember year-end jump
    !
    
print *, 'iDayStart, iDayEnd, daylen_start, daylen_end', iDayStart, iDayEnd, &
      & fu_day_length_hrs(fLat_deg, iDayStart), fu_day_length_hrs(fLat_deg, iDayEnd)
    
    
    if((fu_day_length_hrs(fLat_deg, iDayStart) - fDayLength) * &
     & (fu_day_length_hrs(fLat_deg, iDayEnd) - fDayLength) > 0)then
      !
      ! the required day length is outside the possible range
      ! where?
      if((fDayLength - fu_day_length_hrs(fLat_deg, iDayStart)) * &
       & (fu_day_length_hrs(fLat_deg, iDayStart) - fu_day_length_hrs(fLat_deg, iDayEnd)) > 0)then
        
        if(iDayStart > 0)then
          fu_get_day_from_daylen = -iDayStart
        else
          fu_get_day_from_daylen = -(iDayStart + 365)
        endif

print *, 'Not meaningful request, start side. Needed, edges daylength:', fDayLength, &
             & fu_day_length_hrs(fLat_deg, iDayStart), fu_day_length_hrs(fLat_deg, iDayEnd)

      else
        fu_get_day_from_daylen = -iDayEnd  ! cannot be negative

print *, 'Not meaningful request, end side. Needed, edges daylength:', fDayLength, &
             & fu_day_length_hrs(fLat_deg, iDayStart), fu_day_length_hrs(fLat_deg, iDayEnd)

      endif
      return
    endif  ! meaningful request
    !
    ! midpoint-split search
    do while(abs(iDayEnd - iDayStart) > 1)
      iDay = nint((iDayEnd + iDayStart) / 2.0)   ! interval midpoint
!print *, 'New midpoint: ', iDay
      !
      ! Cut the edge at the same side with iDay
      if((fu_day_length_hrs(fLat_deg, iDay) - fDayLength) * &
       & (fu_day_length_hrs(fLat_deg, iDayStart) - fDayLength) > 0)then
        iDayStart = iDay
      else
        iDayEnd = iDay
      endif
    end do  ! while iDayEnd - iDayStart > 1
    !
    ! return always positive here
    !
    if(iDay > 0)then
      fu_get_day_from_daylen = iDay
    else
      fu_get_day_from_daylen = iDay + 365
    endif
    
print *, 'Final iDay, daylen, requested daylength', iDay, fu_day_length_hrs(fLat_deg, iDay), fDayLength
    
  end function fu_get_day_from_daylen


  !***************************************************************************************************
  
  function fu_expand_environment(fnm) result(cTmp)
    !
    ! Expand environment variables in a given string. The variables must be
    ! specified like ${VARNAME}. A string can have an arbitrary number of variables,
    ! but the result will not be longer than fnlen.
    !
   
    implicit none
    character(len=*), intent(in) :: fnm
    character(len=fnlen) :: cTmp
    
    ! local
    character(len=fnlen) :: value
    character(len=clen) :: varNm
    integer :: idx, idx_end, lenVarNm, iStat
    
    cTmp = fnm
    idx = index(cTmp, '${')

    do while (idx > 0)
      idx_end = index(cTmp, '}')
      lenVarNm = idx_end - idx - 2

      if (lenVarNm < 1) then
        call set_error('Badly formatted string with environments: ' // trim(fnm), 'fu_expand_environment' )
        return
      end if

      varNm = cTmp(idx+2 : idx_end-1)
      call get_environment_variable(trim(varNm), value, STATUS=iStat)

      if (iStat /= 0) then 
         if (iStat == 1) then
           call msg_warning('Environment variable ' // trim(varNm) // ' is not defined')
         else if (iStat == -1) then 
           call msg_warning('Environment variable ' // trim(varNm) // ' is too long')
         else
           call msg_warning('Other trouble with environment variable ' // trim(varNm))
           call msg ("get_environment_variable returned status", iStat)
         endif
         call msg("The line that caused trouble:")
         call msg(fnm)
         call set_error("Trouble with environment", "fu_expand_environment")
         return
      endif

      cTmp = cTmp(1:idx-1) // trim(value) // cTmp(idx_end+1:fnlen)
      idx = index(cTmp,'${')
    end do

  end function fu_expand_environment

  !************************************************************************************

  function fu_process_filepath(path, must_exist, superfile, superdir) result(processed)
   ! Handles a file path:onverts \ or / to 
    ! the native path separator. If specified, checks that the file exists and sets
    ! an error if not.
    implicit none
    character(len=*), intent(in) :: path
    character(len=fnlen) :: processed
    logical, intent(in), optional ::  must_exist
    character(len=*), optional, intent(in) :: superfile !Filename used to expand "^" symbol
    character(len=*), optional, intent(in) :: superdir !Directory used to expand "^" symbol

    integer :: i
    character(len=32) :: readable
    logical :: exists, must_exist_
    character(len=fnlen) :: superdir_

    if (present(must_exist)) then
      must_exist_ = must_exist
    else
      must_exist_ = .false.
    end if
    if (present(superfile)) then  ! if file given as template, get its directory
      if (present(superdir)) then
        call set_error('superfile and superdir cannot be provided simultaneously','fu_process_filepath')
        return
      else
        superdir_ = fu_dirname(superfile)
      end if
    else
      if (present(superdir)) then  ! directory can also be given straight awaqy
        superdir_ = superdir
      else
        superdir_ = ''
      end if
    end if

    processed = path
    
!    call msg('fu_process_filepath1:' + path + ',' + superdir)
    
!    if (convert_slashes_) then
      do i = 1, len_trim(processed)
        if (processed(i:i) == '\'  .or. processed(i:i) == '/' ) then !'
          processed(i:i) = dir_slash
        end if
      end do
!    end if
    
!    call msg('fu_process_filepath2:' + path)

    processed = fu_extend_grads_hat_dir(path,superdir_)
    
 !   call msg('fu_process_filepath3:' + path)

!    if (convert_slashes_) then
      do i = 1, len_trim(processed)
        if (processed(i:i) == '\'  .or. processed(i:i) == '/' ) then !'
          processed(i:i) = dir_slash
        end if
      end do
!    end if

    if (must_exist_) then
      inquire(file=processed, exist=exists, read=readable)
      if (.not. exists .or. readable=='NO') then
        call set_error('File not readable:' + processed, 'fu_process_filepath')
      end if
    end if
    
  end function fu_process_filepath


  !****************************************************************************************

  function fu_index_of_string(string, string_array) result(ind)
    character(len=*), dimension(:), intent(in) :: string_array
    character(len=*), intent(in) :: string
    
    integer :: ind

    integer :: i
    
    ind = int_missing
    
    do i = 1, size(string_array)
      if (string_array(i) == string) then
        ind = i
        return
      end if
    end do
  end function fu_index_of_string


  !****************************************************************************************

  function fu_index_of_int(int, int_array, nVals) result(ind)
    integer, dimension(:), intent(in) :: int_array
    integer, intent(in) :: int
    integer, intent(in), optional :: nVals
    
    integer :: ind

    integer :: i, n
    
    ind = int_missing
    if (present(nVals)) then 
        n=nVals
    else
        n=size(int_array)
    endif

    
      do i = 1, n
        if (int_array(i) == int) then
          ind = i
          return
        end if
      end do
  end function fu_index_of_int

  !*************************************************************************
  
  subroutine matrix_exponent(matrix_2D, nElem)
    !
    ! Computes matrix exponent via Tailor series with scaling and squaring. 
    ! The basis of the method can be found, for example, in:
    ! Moler, C & van Loan, C. 1978. Nineteen dubious ways to compute
    ! exponential of a matrix. SIAM Review, Vol. 20, N:o 4, p 801-836.
    !
    ! Specifically:
    ! b~ = b/2^k
    ! exp(b~) = I + b + b^2/2! + ... = I + sum_over_j (b^j / j!)
    ! exp(b) = (exp(b~))^(2^k)
    !
    ! where b is the matrix, I is the identity matrix, j tends to infinity, k is
    ! determined from the matrix norm
    !
    implicit none
    
    ! Imported parameters
    real(r8k), dimension(:,:) :: matrix_2D
    integer, intent(in) :: nElem
    
    ! Local declarations:
    real(r8k), DIMENSION(:,:), pointer :: term, part_sum   ! the next term and sum of first n terms of Taylor
    INTEGER :: n, i, j, loop
    real :: norm

    !
    ! Accuracy of the approximation of matrix exponent with Taylor series
    ! Teh extra_scales has an experimental value; adding it to the scaling factor
    ! reduces computing time. The best choice seems to be 3, at least it is
    ! faster than other values between 0 and 6.
    !
    REAL, PARAMETER :: matr_exp_epsilon = 1E-12
    INTEGER, PARAMETER :: extra_scales = 3

    if(nElem < 1 .or. nElem > worksize_2dx)then
      call msg('Strange matrix size:',nElem)
      call set_error('Strange matrix size','matrix_exponent')
      return
    endif

    allocate(term(nElem,nElem),part_sum(nElem,nElem),stat=i)
    if(fu_fails(i==0,'Failed allocation of temporary arrays','matrix_exponent'))return

    part_sum = 0.
    term = 0.

    ! Initialize sum and the term to be the identity matrices:
    DO i = 1, nElem
!      part_sum(i,i) = 1
      term(i,i) = 1
    END DO

    norm = fu_matrix_norm(matrix_2d, nElem)

    ! Find out suitable n: it is the smallest power of 2 so that
    ! norm(b) / n <= 1. extra_scales is an experimental value that makes
    ! computation faster in most cases. Equals to 3.
    !
    ! n = int(LOG(norm) / ln_2) + extra_scales
    n = exponent(norm) + extra_scales
    IF (n < 0) n = 0
    !
    ! Scale the matrix:
    !
    matrix_2d(1:nElem,1:nElem) = matrix_2d(1:nElem,1:nElem) / (2.0**n)
    !!matrix_2d(1:nElem,1:nElem) = matrix_2d(1:nElem,1:nElem) * SET_EXPONENT(1.0, -(n+1))
    !
    ! This loop computes the Taylor:
    !
    DO loop = 1, 100 ! not more than 100 Tyalor terms
      IF (all (abs(term(1:nElem,1:nElem)) < matr_exp_epsilon)) exit
      part_sum(1:nElem,1:nElem) = part_sum(1:nElem,1:nElem) + term(1:nElem,1:nElem)
      term(1:nElem,1:nElem) = MATMUL(term(1:nElem,1:nElem), matrix_2D(1:nElem,1:nElem)) / real(loop)
    END DO  ! Taylor terms
    call msg('number of Taylor terms used:',loop-1)
    !
    ! Return back the scaling
    !
    DO loop = 1, n
      part_sum(1:nElem,1:nElem) = MATMUL(part_sum(1:nElem,1:nElem), part_sum(1:nElem,1:nElem))
    END DO
    matrix_2d(1:nElem,1:nElem) = part_sum(1:nElem,1:nElem)

    deallocate(term,part_sum)

  end subroutine matrix_exponent


  ! ***************************************************************

  real(r8k) FUNCTION fu_max_diagonal (matrix_2D, dim)
    !
    ! Finds the largest diagonal element of a square matrix.

    IMPLICIT NONE

    ! Imported parameters
    real(r8k), DIMENSION(:,:) :: matrix_2D
    INTEGER, INTENT(IN) :: dim       ! dimension of matrix

    ! Local declarations
    INTEGER :: i

    fu_max_diagonal = 0
    DO i = 1, dim
      fu_max_diagonal = MAX(fu_max_diagonal, ABS(matrix_2D(i,i)))
    END DO

  END FUNCTION fu_max_diagonal


  ! ***************************************************************

  real(r8k) FUNCTION fu_matrix_norm (matrix_2D, dim)
    !
    ! This function calculates the norm of an argument matrix. The
    ! method used is described in ORIGEN2 - A Revised and Updated
    ! Version of the Oak Ridge Isotope Generation an d Depletion Code
    ! by A. G. Croff, 1980. Actually the idea is to compute both
    ! the 1-norm and infinity-norm of the matrix and then choose
    ! the smaller of those.

    IMPLICIT NONE

    real(r8k), DIMENSION(:,:) :: matrix_2D
    INTEGER, INTENT(IN) :: dim

    REAL :: max_row, max_col
    INTEGER :: i

    max_row = 0
    max_col = 0
    DO i = 1, dim
      max_row = MAX (SUM(ABS(matrix_2D(i,1:dim))), max_row)
      max_col = MAX (SUM(ABS(matrix_2D(1:dim,i))), max_col)
    END DO
    fu_matrix_norm = MIN (max_row, max_col)

  END FUNCTION fu_matrix_norm


  !*************************************************************************

  subroutine advect_Galperin_bulk_1d_relat(nCells, nSpecies, &
                                         & vMass, vMassCentre, &
                                         & vWind_relative, &
                                         & vMinAdvectedMass, garbage)
    !
    ! A semi-Lagrangian scheme of M.Galperin, as in (Galperin, 2000)
    ! presented at NATO/CCMS and published in IGCE series. 
    !
    ! Advection goes on in relative units, so that wind is in grid cells per time step
    !
    implicit none

    ! Imported parameters
    integer, intent(in) :: nCells, nSpecies
    real, dimension(1:nSpecies, 0:nCells+1), intent(inout) :: vMass
    real, dimension(0:nCells+1), intent(inout) :: vMassCentre
    real, dimension(:), pointer :: vWind_relative, garbage, vMinAdvectedMass

    ! Local variables
    integer :: ix, JR, JL, JC
    real :: RL, RR, fXP, fXC, SS, BR, BL, R, Padv, Qadv, fMTotal
    real, dimension(:), pointer :: vXMomentTmp
    real, dimension(:,:), pointer :: vMassTmp

    !
    ! Preparatory steps
    !
    vMassTmp => fu_work_array_2D()
    vXMomentTmp => fu_work_array()
    if(error)return

    vMassTmp(1:nSpecies, 1:nCells+2) = 0.
    vXMomentTmp(1:nCells+2) = 0.


    !-----------------------------------------------------------------
    !
    ! Cycle over the 1D grid
    !
    DO ix=1,nCells
      !
      ! Check the low-cutting threshold
      !
      fMTotal = sum(vMass(1:nSpecies,ix))

      if(all(vMass(1:nSpecies,ix) < vMinAdvectedMass(:)))then
        !
        ! Mass is small. Drop all to source-related PFB and continue
        !
        garbage(1:nSpecies) = garbage(1:nSpecies) + vMass(1:nSpecies,ix)
        vMass(1:nSpecies,ix) = 0.

      else
        !
        ! Mass is not small, make the advection. 
        !
        ! ATTENTION. Here we store centres of masses rather than the momentum
        ! Reason: with centres of masses all chemistry and deposition are straightforward
        ! Further on, however, we will also involve the momentum because it is much
        ! more convenient in advection routine
        !
        fXC = vMassCentre(ix) ! x-coord
        fXP = fXC * fMTotal
        fXC = fXC + ix

        RL=fXC-ix+0.5  ! ix-0.5 is a left border of the cell. RL is distance from left border
        RR=1.-RL       ! RR is a distance from the right cell border
        R=MIN(RL,RR)

!              if(ix <= nx_dispersion*0.15)then
!                wind = 7.1
!              elseif(ix <= nx_dispersion*0.5)then
!                wind = 7.1 + (ix - nx_dispersion*0.15)*1.3
!              elseif(ix <= nx_dispersion*0.85)then
!                wind = 7.1 + (nx_dispersion*0.85 - ix)*1.3
!              else
!                wind = 7.1
!              endif
!
!              wind = 10.0

        !
        ! Move along x-wind
        !
        fXC=fXC + vWind_relative(ix)

!              fXC=fXC + 0.3 

        BR=fXC+R   ! Right border
        BL=fXC-R   ! Left border
        SS=BR-BL   ! Width of the puff === 2 R
!              JR=INT(BR+0.5)
!              JL=INT(BL+0.5)
!              JC=INT(fXC+0.5)
        JR=min(max(0,INT(BR+0.5)),nCells+1)  ! If Currant number > 1, 1 border-cell
        JL=min(max(0,INT(BL+0.5)),nCells+1)  ! can be not enough, so let's take care of it
        JC=min(max(0,INT(fXC+0.5)),nCells+1)

        !
        ! Here the horizontal diffusion is taken from the passed path. Actually, it may be
        ! computed from actual wind speed and direction but so far we take sigma=0.1*path
        ! With coordinate-wise split it has to be applied to coth dimensions
        !

        IF(JR == JL .OR. SS <= 0.0001) THEN
          !
          ! The whole mass goes to another cell
          !
          vMassTmp(1:nSpecies, JC+1) = vMassTmp(:,JC+1) + vMass(:,ix)
          vXMomentTmp(JC+1) = vXMomentTmp(JC+1) + (fXC-JC) * fMTotal

        ELSE
          !
          ! Only part of the mass goes to another cell
          !
          Padv = (BR - JR + 0.5)/ss ! Fraction of the mass appeared in the right-hand cell
          Qadv = Padv * fMTotal   ! Mass appearing in the right-hand cell

          vMassTmp(1:nSpecies,JR+1) = vMassTmp(:,JR+1) + Padv * vMass(:,ix)
          vXMomentTmp(JR+1) = vXMomentTmp(JR+1) + (BR - JR - 0.5) * 0.5 * Qadv

          PADV=1.-PADV  ! Fraction appearing in the left-hand cell
          Qadv=Padv * fMTotal  ! Mass appearing in the left-hand cell
          vMassTmp(1:nSpecies,JL+1) = vMassTmp(:,JL+1) + Padv * vMass(:,ix)
          vXMomentTmp(JL+1) = vXMomentTmp(JL+1) + (BL - JL + 0.5) * 0.5 * Qadv ! JR is correct!

        END IF

      endif  ! Mass checking
    end do  ! Cycle iy. X-advection is over

    !
    ! Put the tmp arrays back to the main mass ones
    !
    do ix = 0, nCells+1
      vMass(1:nSpecies,ix) = vMassTmp(:,ix+1)
      fMTotal = sum(vMass(1:nSpecies,ix))
      if(fMTotal > 0)then
        vMassCentre(ix) = vXMomentTmp(ix+1) / fMTotal
      else
        vMassCentre(ix) = 0.
      endif
    end do

    call free_work_array(vMassTmp)
    call free_work_array(vXMomentTmp)

  end subroutine advect_Galperin_bulk_1d_relat


  !*************************************************************************

  subroutine advect_Galperin_bulk_1d_abs(nCells, nSpecies, &
                                       & vMass, vMassCentre, &
                                       & vWind_abs, vCellBorder, timestep_sec, &
                                       & vMinAdvectedMass, garbage)
    !
    ! A semi-Lagrangian scheme of M.Galperin, as in (Galperin, 2000)
    ! presented at NATO/CCMS and published in IGCE series. 
    !
    ! Advection goes on in absolute units, so that wind and grid cell sises must be in the same units
    ! Time is then seconds
    !
    ! ALL UNITS: SI
    !
    implicit none

    logical, parameter :: DMAT_advection = .false.

    ! Imported parameters
    integer, intent(in) :: nCells, nSpecies
    real, dimension(1:nSpecies, 0:nCells+1), intent(inout) :: vMass
    real, dimension(0:nCells+1), intent(inout) :: vMassCentre, vCellBorder
    real, dimension(:), pointer :: vWind_abs, & ! must be same linear units!
                                 & garbage, vMinAdvectedMass
    real, intent(in) :: timestep_sec

    ! Local variables
    integer :: ix, JR, JL, JC, iPresent, iFuture, iTime, nTimeSteps
    real :: RL_abs, RR_abs, fXP, fXC, BR_abs, BL_abs, R_abs, Padv, Qadv, fMTotal, fTimeStepSmall, fCourantNbr
    real, dimension(:,:), pointer :: vXMomentTmp
    real, dimension(2,nSpecies,0:nCells+1) :: vMassTmp
    logical :: ifValid, ifGarbage

    !
    ! Preparatory steps
    !
    vXMomentTmp => fu_work_array_2D()
    if(error)return

    vMassTmp(1:2,1:nSpecies, 0:nCells+1) = 0.
    vXMomentTmp(1:2,1:nCells+2) = 0.

    !
    ! The absolute advection is usually taken when the grid cells are substantially uneven.
    ! Then the Courant number can vary largely and the risk of jumping over the small cells
    ! shows up. So, we have to make an effort to keep the max Courant number below 1.
    !
    fCourantNbr = abs(timestep_sec * vWind_abs(1)/(vCellBorder(1) - vCellBorder(0)))
    do ix = 2, nCells
      if(fCourantNbr < abs(timestep_sec * vWind_abs(ix)/(vCellBorder(ix)-vCellBorder(ix-1))))then
        fCourantNbr = abs(timestep_sec * vWind_abs(ix)/(vCellBorder(ix)-vCellBorder(ix-1)))
      endif
    end do
    if(fCourantNbr <= 0.5)then
      fTimeStepSmall = timestep_sec
      nTimeSteps = 1
    else
      nTimeSteps = int(2.*fCourantNbr+1)
      fTimeStepSmall = timestep_sec / real(nTimeSteps)
!        call msg('Max Courant at index:',ix+(iy-1)*nx_meteo,fCourantNbr)
    endif

    !
    ! Cycle over small time steps
    !
    iPresent = 1
    do iTime = 1, nTimeSteps

      iFuture = 3 - iPresent
      ifValid = .false.
      !-----------------------------------------------------------------
      !
      ! Cycle over the 1D grid
      !
      DO ix=1,nCells
        !
        ! Check the low-cutting threshold
        !
        if(iTime == 1)then
          ifGarbage = all(vMass(1:nSpecies,ix) < vMinAdvectedMass(:))
          fMTotal = sum(vMass(1:nSpecies,ix))
        else
          ifGarbage = all(vMassTmp(iPresent,1:nSpecies,ix) < vMinAdvectedMass(:))
          fMTotal = sum(vMassTmp(iPresent,1:nSpecies,ix))
        endif

        if(ifGarbage)then  !fMTotal < fMinAdvectedMass)then
          !
          ! Mass is small. Drop all to source-related PFB and continue
          !
          garbage(1:nSpecies) = garbage(1:nSpecies) + vMass(1:nSpecies,ix)
          vMass(1:nSpecies,ix) = 0.
          cycle

        else  
          !
          ! Mass is not small, make the advection using the small time steps
          !
          ifValid = .true.
          !
          ! fCellBorder and vMassCentre are in absolute units, both measured from the left side 
          ! of the line. fCellBorder(0) is the left end of the line, corresponds to 0.5 in relative
          ! coordinates.
          !
          if(iTime == 1)then
            fXC = vMassCentre(ix)
            fXP = (fXC - 0.5*(vCellBorder(ix-1)+vCellBorder(ix))) * fMTotal
          else
            fXP = vXMomentTmp(iPresent,ix+1) ! z-coord
            fXC = fXP / fMTotal + 0.5*(vCellBorder(ix-1)+vCellBorder(ix))
          endif

          if(fXC < vCellBorder(0) .or. fXC > vCellBorder(nCells+1))then
            call msg('Resetting strange fXC at ix:',fXC,ix)
            fXC = 0.5*(vCellBorder(ix-1)+vCellBorder(ix))
          endif

!          if(if_report_adv)then
!            call msg('abs-adv:ix,fXC:',ix,fXC)
!            call msg('abs-adv:massTotal:',fMTotal)
!          endif

          !
          ! All below computations are going in absolute coordinates, NOT grid indices
          ! Note that for z-coord the absolute levels can be pressure, height, altitude,...
          !
          RL_abs = abs(fXC - vCellBorder(ix-1))
          RR_abs = abs(vCellBorder(ix) - fXC)
          R_abs = min(RL_abs, RR_abs)
          !
          ! Move the centre of mass along the Z-axis. Wind has already been selected 
          ! to have right unit
          !

!                wind = -0.01

          !
          ! Attention: wind can push some part of the species below the surface level.
          ! It means that the 1st level wind is downward. Since so far we do not include the
          ! non-divergent 10m wind with zero vertical component, have to just assume it.
          ! This, in particular, means that the substances can NEVER go underground. We force the
          ! model to keep at least 1% of the first layer thickness above the surface.
          !
          fXC = fXC + vWind_abs(ix) * fTimeStepSmall
          BR_abs = fXC + R_abs ! Right border (upper), absolute
          BL_abs = fXC - R_abs ! Left border (lower), absolute

          !
          ! Find the level where these three points appear
          !
          if(BR_abs < 0.0)then  
            JR = 0
!                  call msg('Underground, JR')
          else
            do JR = 1,nCells
              if((BR_abs - vCellBorder(JR)) * (BR_abs - vCellBorder(JR-1)) <= 0.)exit
            end do
          endif
          if(BL_abs < 0.0)then
            JL = 0
!                  call msg('Underground, JL')
          else
            do JL = 1,nCells
              if((BL_abs - vCellBorder(JL)) * (BL_abs - vCellBorder(JL-1)) <= 0.)exit
            end do
          endif
          if(fXC < 0.0)then
            JC = 0
!                  call msg('Underground, JC')
          else
            do JC = 1,nCells
              if((fXC - vCellBorder(JC)) * (fXC - vCellBorder(JC-1)) <= 0.)exit
            end do
          endif

!                do iSubst = 1, pDispFlds%nSubst
!                  call msg('Z advection',ix,sum(colMassTmp(ix)%subst(iSubst)%mass(1:pDispFlds%nModes(iSubst))))
!                enddo

          if(DMAT_advection)then
            !
            ! In DMAT, all mass always goes/stays in the cell defined by JC
            !
            if(iTime == 1)then
              vMassTmp(iFuture,1:nSpecies,JC) = &
                                 & vMassTmp(iFuture,:,JC) + vMass(:,ix)
            else
              vMassTmp(iFuture,1:nSpecies,JC) = &
                                           & vMassTmp(iFuture,:,JC) + vMassTmp(iPresent,:,ix)
            endif
            vXMomentTmp(iFuture,JC+1) = vXMomentTmp(iFuture,JC+1) + &
                                            & (fXC - 0.5*(vCellBorder(JC)+vCellBorder(JC-1))) * fMTotal
          else
            !
            ! In Galeprin(2000) advection, more options are available
            !

            IF(JR == JL .OR. (BR_abs-BL_abs < 0.001 * (vCellBorder(JL) - vCellBorder(JL-1)))) THEN
              !
              ! The whole mass goes to another cell: thin plume that fits inside borders
              ! or simply the delta-function-type plume
              !
              if(iTime == 1)then
                vMassTmp(iFuture,1:nSpecies,JC) = &
                                   & vMassTmp(iFuture,:,JC) + vMass(:,ix)
              else
                vMassTmp(iFuture,1:nSpecies,JC) = &
                                             & vMassTmp(iFuture,:,JC) + vMassTmp(iPresent,:,ix)
              endif
              vXMomentTmp(iFuture,JC+1) = vXMomentTmp(iFuture,JC+1) + &
                                              & (fXC - 0.5*(vCellBorder(JC)+vCellBorder(JC-1))) * fMTotal
            ELSE
              !
              ! Only part of the mass goes to another cell
              !
              Padv = (BR_abs - vCellBorder(JR-1)) / (BR_abs - BL_abs)  ! Fraction to the right-hand cell, absolute coord

              Qadv=Padv * fMTotal  ! Mass appearing in the right-hand cell
              if(iTime == 1)then
                vMassTmp(iFuture,1:nSpecies,JR) = vMassTmp(iFuture,:,JR) + Padv * vMass(:,ix)
              else
                vMassTmp(iFuture,1:nSpecies,JR) = vMassTmp(iFuture,:,JR) + Padv * vMassTmp(iPresent,:,ix)
              endif
              vXMomentTmp(iFuture,JR+1) = vXMomentTmp(iFuture,JR+1) + &
                                                              & (BR_abs - vCellBorder(JR)) * 0.5 * Qadv

              Padv = 1. - Padv  ! Fraction appearing in the left-hand cell
              Qadv = fMTotal - Qadv
              if(iTime == 1)then
                vMassTmp(iFuture,1:nSpecies,JL) = vMassTmp(iFuture,:,JL) + Padv * vMass(:,ix)
              else
                vMassTmp(iFuture,1:nSpecies,JL) = vMassTmp(iFuture,:,JL) + Padv * vMassTmp(iPresent,:,ix)
              endif

              vXMomentTmp(iFuture,JL+1) = vXMomentTmp(iFuture,JL+1) + &
                                                            & (BL_abs - vCellBorder(JL-1)) * 0.5 * Qadv
            END IF  ! whole mass or split
          endif

          vMassTmp(iPresent,:,ix) = 0.
          vXMomentTmp(iPresent,ix+1) = 0.

        endif  ! Mass checking
      end do  ! over cells
            
      if(.not.ifValid)exit
      iPresent = 3 - iPresent

    end do  ! over small time steps

    !
    ! Put the tmp arrays back to the main mass ones
    !
    do ix = 0, nCells+1
      vMass(1:nSpecies,ix) = vMassTmp(iFuture,:,ix)
      fMTotal = sum(vMass(1:nSpecies,ix))
      if(fMTotal > 0)then
        vMassCentre(ix) = vXMomentTmp(iFuture,ix+1) / fMTotal  + 0.5*(vCellBorder(ix-1)+vCellBorder(ix))
      else
        vMassCentre(ix) =  0.5*(vCellBorder(ix-1)+vCellBorder(ix))
      endif
    end do

    call free_work_array(vXMomentTmp)

  end subroutine advect_Galperin_bulk_1d_abs

!    !***************************************************************************************************                                   
!                                       
!      subroutine avect_lineparams_relative(u, lineParams, lineParamsOut, nCells, time_dir, seconds_abs, iXStart, iXEnd)
!  
!!        real, dimension(-1:), intent(in) :: lineParams  ! in actual units
!!                                                        ! the right/upper boundary of the cell
!        real, dimension(0:), intent(in) :: u            ! wind, actual units
!        real, dimension(-1:), intent(out) :: lineParamsOut ! in actual units
!        integer, intent(in) :: nCells !, iXStart, iXEnd   ! size of line array, start&end of advection
!        real, intent(in) :: time_dir !, seconds_abs       ! absolute timestep and direction in time + or -1
!        !
!        ! Gets locations of cell boundaries after advection
!        ! trick: for backward time call with inverse velocity
!        ! all velocities are forced to zero at ix=0, so no advection
!        ! to underground is possible.
!
!        real ::  alpha, alphat, curpos, nextpos, dz,t
!        integer ::  ix, curix, nextix, ixTo, dir, time_dir
!        logical :: useExp
!
!        time_dir = int(sign(1.0,seconds))
!
!        !
!        ! Scan the line from start to end advecting cell by cell
!        !
!        do ix = iXStart, iXEnd
!          t = abs_seconds
!          curix = ix
!          dir = int(sign(1.0, u(ix))) ! +1 rightwards;  -1  leftwards
!          !
!          ! Where will it land? Search for it
!          !
!          do iCount = n, nCells   ! not more than that - but we should get out before
!            nextix = curix + dir
!            curpos = lineParams(curix)
!            !
!            ! If we do not want to send mass out of the domain, the speed=0 at the border
!            ! and linearly decreasing when approaching it
!            !
!            if ((nextix < 1) .or. (nextix > nCells)) then
!              !
!              ! Off the domain. Global? Closed border?
!              !
!              if(ifGLobal)then
!                !
!                ! Global grid => just turn the line over
!                if(nextix < 1) nextix = nextix + nCells
!                if(nextix > nCells) nextix = nextix - nCells
!                dx = lineParams(nextix) - lineParams(curix)     ??????????
!                alpha = (u(nextix) - u(curix)) / dx ! du/dx
!                alphat= alpha * t
!                useExp = (abs(alphat) > 0.001 )
!                if (useExp) then
!                  nextpos = curpos + u(curix) * (exp(alphat) - 1.0)/alpha     >>??????
!                else
!                  nextpos = curpos + u(curix) * t                              ???????????
!                endif
!
!              elseif(ifClosedBorders)then
!                !
!                ! Closed borders => the mass has to be stopped nicely in front of the border
!                !
!                dx = lineParams(nextix) - lineParams(curix)
!                alpha = (0. - u(curix)) / dx                    ! du/dx
!                alphat = alpha * t
!                if (abs(alphat) > 0.001 ) then
!                  nextpos = curpos + u(curix) * (exp(alphat) - 1.0)/alpha
!                else
!                  nextpos = curpos + u(curix) * t
!                endif
!                lineParamsOut(ix) = nextpos
!                exit   ! border finished the trip
!              else
!                !
!                ! Transparent border: what is gone is gone... but what do we do with the border?
!                ! It cannot be sent further than one grid cell off the domain: no wind data.
!                ! Here the border is allowed to go anywhere inside the 0-th cell, then abrupt stop
!                !
!                dx = lineParams(nextix) - lineParams(curix)
!                alpha = (u(nextix) - u(curix)) /  dx ! du/dx
!                alphat= alpha * t
!                if (abs(alphat) > 0.001) then
!                  move = u(curix) * (exp(alphat) - 1.0) / alpha
!                else
!                  move = u(curix) * t
!                endif
!                if(abs(move) < abs(dx))then
!                  lineParamsOut(ix) = curpos + move  ! move is small, use it
!                else
!                  lineParamsOut(ix) = curpos + dx    ! move sends the border beyond the out-of-domain cell. stop
!                endif
!                exit  ! border finished the trip
!              endif   ! global/closed_borders?
!            else
!              !
!              ! Normal move inside the domain
!              !
!              dx = lineParams(nextix) - lineParams(curix)
!              alpha = (u(nextix) - u(curix)) / dx ! du/dx
!              alphat= alpha * t
!              useExp = (abs(alphat) > 0.001 )
!              if (useExp) then
!                nextpos = curpos + u(curix) * (exp(alphat) - 1.0)/alpha
!              else
!                nextpos = curpos + u(curix) * t
!              endif
!            endif  ! in/out of the domain
!            !
!            ! Next position is suggested. Is it further than the current cell or stops here?
!            !
!            if (nextpos * dir < dir * lineParams(nextix)) then 
!              !
!              ! Move over more than one grid cell
!              ! How long it takes to pass the current cell?
!              !
!              if (useExp) then
!                t = t - log(1.0 + alpha * dx / u(curix)) / alpha
!              else
!                t = t - dx / u(curix)
!              endif
!              curix = nextix
!            else
!              lineParamsOut(ix) = nextpos  ! the trip finishes here
!              exit
!            endif  ! move over 1 or more cells
!          enddo ! going along the cells looking for the destination 
!!#ifdef DEBUG
!  if (.not. lineParamsOut(ix)>=0.)then
!    call msg("Negative cell boundaries after advection!")
!    call msg("Levels from/to:",iXStart,iXEnd)
!    call msg("ix=", ix)
!    call msg("Not all colParamsOut are positive!") 
!    call msg("In ", lineParams(-1:nCells+1))
!    call msg("Out", lineParamsOut(-1:nCells+1))
!    if(useExp)then 
!      call msg("useExp: true")
!    else
!      call msg("useExp: false")
!    endif
!    call msg(" u(curix), u(nextix):", u(curix), u(nextix))
!    call msg("u  ", u(0:nCells))
!    call set_error("Not all lineParamsOut are positive!", "avect_lineparams")
!  endif
!!#endif
!
!        enddo  ! cells
!      end subroutine avect_lineparams_relative
!
!    !========================================================================
!
!      subroutine advect_mass_line_rel(massIn,cmIn, &   ! input mass and its centre. (0:nx+1)
!                                    & massOut,cmOut, & ! output mass and its cenrte, (0:nx+1)
!                                    & passengers, passengersOut, nPassengers, &  ! masses etc transported same way as massIn
!                                    & b_rightAdvected, & ! right-hand borders before/after advection
!                                    & nx, &    ! domain size, valid range, time step
!                                    & fMinAdvectedMass, &  ! low-mass threshold
!                                    & ifAdvAbs)            ! if advecting absolute mass
!        ! 
!        ! Actual advection of massIn using its centre of mass cmIn and the advection itself
!        ! stored as pre-advected borderes b_right and b_rightAdvected. Essentially, slabs are
!        ! projected to the borders and their new position is an interpolation between the
!        ! new position of the borders.
!        ! Passengers is a bunch of variables that can be summed-up and interpolated (masses, moments)
!        ! and which are advected in the same way as the corresponding species in massIn. For instance,
!        ! advection along x-axis requires x-centre of mass for doing it, whereas y- and z-moments are
!        ! transported the same way as mass. Similarly, some other masses can be passengers: e.g., masses
!        ! that constitute aerosols: number conc is advected, masses are passengers.
!        ! 
!        real, dimension(0:), intent(in)  :: massIn, cmIn
!        real, dimension(-1:), intent(in)  :: b_rightAdvected
!        real, dimension(0:), intent(out) :: massOut,cmOut
!        type(silja_rp_1d), dimension(:), pointer :: passengers
!        integer, intent(in) :: nx, nPassengers !, ixStart, ixEnd
!        real, intent(in) :: fMinAdvectedMass !, seconds
!        logical, intent(in) :: ifAdvAbs 
!
!        real :: fC, SS, RL_abs, R_abs, RR_abs, BR_abs, BL_abs, &
!              & ftmp, fracleft, fmass, fmassto, BL_tmp, BR_tmp
!        integer ::  ix, ixTo  
!
!        massOut(0:nx+1) = 0.0
!        cmOut(0:nx+1) = 0.0
!        do iPass = 1, nPassengers
!          passengersOut(iPass)%rp(0:nx+1) = 0.
!        enddo
!!        momXout(0:nx+1) = 0.0
!!        momYout(0:nx+1) = 0.0
!
!! !No advection hack
!!      massOut(0:nlev+1)=massIn(0:nlev+1)
!!      cmOut(0:nlev+1)= ZcmIn(0:nlev+1)
!!      do iPass = 1, nPassengers
!!        passengersOut(iPass)%rp(0:nx+1) = passengersIn(iPass)%rp(0:nx+1)
!!      enddo
!!!      momXout(0:nlev+1)=momXIn(0:nlev+1)
!!!      momYout(0:nlev+1)=momYIn(0:nlev+1)
!!      return
!
!        ixTo = 0
!        cmOut(ixTo) = 0.5 * (b_right(ixTo-1) + b_right(ixTo))
!        !
!        ! Scan all masses in the input vector advecting them one by one.
!        !
!        do ix = ixStart, ixEnd  ! range of cells to be advected
!          fmass = massIn(ix)
!          if (abs(fmass) < fMinAdvectedMass) cycle   ! Garbage should be already collected 
!        
!          fC = cmIn(ix)
!
!          if (ifAdvAbs .and. (fC < b_right(ix) .or. fZC > b_right(ix-1))) then
!            ! Return mass back to cell 
!            fZC = 0.5 * (b_right(ix-1) + b_right(ix))
!          end if
!          ! Left position
!          RL_abs = b_right(ix-1) -fC
!          RR_abs = fC - b_right(ix)
!          R_abs = min(RL_abs, RR_abs)
!  
!          !fraction, of cell occupied by slab
!          ftmp = 2*R_abs / (b_right(ix-1)-b_right(ix))
!
!          ! Slab is not allowed to occupy less than 1e-3 of cell,
!          ! otherwise we get into numerics. Indeed, the target slab should 
!          ! be checked, but it is more complicated...
!          ftmp = max(ftmp, 1e-3)
!
!          if (RL_abs < RR_abs) then
!            ! on the lower side of the layer
!            BL_abs = b_rightAdvected(ix-1)
!            BR_abs = ftmp*b_rightAdvected(ix) + (1.0-ftmp)*BL_abs
!          else
!            BR_abs = b_rightAdvected(ix)
!            BL_abs = ftmp * b_rightAdvected(ix-1) + (1.0-ftmp)*BR_abs
!          end if
!
!          do while (ixTo <= nx .and. BL_abs <= b_right(ixTo))
!            ixTo = ixTo+1 
!            cmOut(ixTo)=0.5*(b_right(ixTo-1)+b_right(ixTo))
!          enddo
!          !now BL_abs is just above ixTo
!          BL_tmp = BL_abs
!          BR_tmp = max(BR_abs, b_right(ixTo)) !lowest of two
!          SS = BL_abs - BR_abs !new size
!          if (SS > 1e-10) then  ! thick slab
!            do while (ixTo <= nx .and. BR_abs < b_right(ixTo))
!              !
!              ! Except for upper part of advected slab
!              ftmp = (BL_tmp - BR_tmp) / SS !Fraction that goes here
!              !
!              ! Passengers: masses and moments advected the same way as the current mass
!              do iPass = 1, nPassengers
!                passengersOut(iPass)%rp(ixTo) = passengersOut(iPass)%rp(ixTo) + &
!                                              & passengersIn(iPass)%rp(ix) * ftmp
!              enddo
!!              momXout(ixTo) = momXout(ixTo) + momXin(ix) * ftmp
!!              momYout(ixTo) = momYout(ixTo) + momYin(ix) * ftmp
!              fmassTo = massOut(ixTo) ! mass here before
!              fTmp = fTmp * fmass ! mass that goes here
!              massOut(ixTo) = fmassTo + ftmp
!              cmOut(ixTo) = (cmOut(ixTo)*fmassTo + 0.5*(BR_tmp+BL_tmp)*ftmp) / massOut(ixTo)
!              ixTo = ixTo+1
!              cmOut(ixTo)=0.5*(b_right(ixTo-1)+b_right(ixTo))
!              BL_tmp = BR_tmp
!              BR_tmp = b_right(ixTo)
!            enddo
!            fTmp =  (BL_tmp - BR_abs) / SS
!            do iPass = 1, nPassengers
!              passengersOut(iPass)%rp(ixTo) = passengersOut(iPass)%rp(ixTo) + &
!                                            & passengersIn(iPass)%rp(ix) * ftmp
!            enddo
!!            momXout(ixTo) = momXout(ixTo) + momXin(ix) * ftmp
!!            momYout(ixTo) = momYout(ixTo) + momYin(ix) * ftmp
!            fmassto = massOut(ixTo)
!            fTmp = fTmp * fmass
!            massOut(ixTo) = fmassto + ftmp
!            cmOut(ixTo) = (cmOut(ixTo)*fmassto + 0.5*(BR_abs+BL_tmp)*fTmp) / massOut(ixTo)
!          else ! zero-thickness slab
!            do iPass = 1, nPassengers
!              passengersOut(iPass)%rp(ixTo) = passengersOut(iPass)%rp(ixTo) + passengersIn(iPass)%rp(ix)
!            enddo
!!            momXout(ixTo) = momXout(ixTo) + momXin(ix)
!!            momYout(ixTo) = momYout(ixTo) + momYin(ix) 
!            fmassto = massOut(ixTo)
!            massOut(ixTo) = fmassto + fmass
!            cmOut(ixTo) = (cmOut(ixTo)*fmassto + BR_abs*fmass )/ massOut(ixTo)
!          endif
!
!        enddo ! ix
!        do ix = ixTo+1, nx+1 
!          cmOut(ix) = 0.5 * (b_right(ix-1)+b_right(ix))
!        enddo
!
!#ifdef DEBUG    
!  if (.not. all(massOut(0:nx+1) >= 0.))then
!    call msg("Not all massOut(0:nlev+1) are positive after advection!" )
!    call msg("b_right, in", b_right(-1:nLev+1))
!    call msg("b_right_advected, out", b_rightAdvected(-1:nLev+1))
!    call msg("Mass, in", massIn(0:nLev+1))
!    call msg("Mass, out", massOut(0:nLev+1))
!    call msg("cm, in", cmIn(0:nLev+1))
!    call msg("cm, out", cmOut(0:nLev+1))
!    call msg("cm, centers", 0.5*(b_right(-1:nLev)+b_right(0:nLev+1)))
!    call set_error("Negative mass","advect_mass_line")
!  endif
!#endif
!
!      end subroutine advect_mass_line_rel

  !************************************************************************************************

  subroutine split_string_int(string, separator, int_array, num_tokens)
    ! 
    ! Split a string consisting of fields separated by a given
    ! character and convert the results to integers. Example:
    ! call split_string(' 12 4 6   8 9 111  ', ' ', array, n)
    ! returns
    ! array = (/12, 4, 6, 8, 9, 111/) 
    ! n = 6
    !
    ! The separator must have length 1, successive separators are merged.
    !
    implicit none
    character(len=*), intent(in) :: string
    character(len=1), intent(in) ::  separator
    integer, dimension(:), intent(out) :: int_array
    integer, intent(out) :: num_tokens

    integer :: len_arr, iostat, value, next_sep, next_non_sep, old_next_non_sep, dummy
    character(len=len(string)) :: token

    len_arr = size(int_array)
    !if (fu_fail(len(separator) /= 1, 'len(separator) must be 1', 'split_string_int')) return
    
    num_tokens = 0
    next_non_sep = fu_next_non_sep(string, separator, 1)

    do while (next_non_sep > 0)
      next_sep = fu_next_sep(string, separator, ind_start=next_non_sep)
      token = string(next_non_sep:next_sep-1)
      !print *, next_non_sep, next_sep, 'tok:', token, 'nt:', num_tokens
      read(unit=token, fmt=*, iostat=iostat) value
      if (fu_fails(iostat == 0, 'Failed to parse:' // token, 'split_string_int')) return
      num_tokens = num_tokens + 1
      if (fu_fails(num_tokens <= len_arr, 'Too many tokens', 'split_string_int')) return
      int_array(num_tokens) = value
      next_non_sep = fu_next_non_sep(string, separator, next_sep)
    end do

  end subroutine split_string_int


  !************************************************************************************************

  subroutine split_string_real(string, separator, real_array, num_tokens)
    ! As above, but for reals.
    !
    implicit none
    character(len=*), intent(in) :: string, separator
    real, dimension(:), intent(out) :: real_array
    integer, intent(out) :: num_tokens

    integer :: len_arr, iostat, next_sep, next_non_sep, old_next_non_sep, dummy
    real :: value
    character(len=len(string)) :: token

    len_arr = size(real_array)
    if (fu_fails(len(separator) == 1, 'len(separator) must be 1', 'split_string_real')) return
    
    num_tokens = 0
    next_non_sep = fu_next_non_sep(string, separator, 1)

    do while (next_non_sep > 0)
      next_sep = fu_next_sep(string, separator, next_non_sep)
      token = string(next_non_sep:next_sep-1)
      read(unit=token, fmt=*, iostat=iostat) value
      if (fu_fails(iostat == 0, 'Failed to parse:' // token, 'split_string_real')) return
      num_tokens = num_tokens + 1
      if (fu_fails(num_tokens <= len_arr, 'Too many tokens', 'split_string_real')) return
      real_array(num_tokens) = value
      next_non_sep = fu_next_non_sep(string, separator, next_sep)
    end do

  end subroutine split_string_real

  !************************************************************************************
  
  subroutine split_string_char(string, separator, char_array, num_tokens)
    ! As above, but for strings.
    !
    implicit none
    character(len=*), intent(in) :: string, separator
    character(len=*), dimension(:), intent(out) :: char_array
    integer, intent(out) :: num_tokens

    integer :: len_arr, next_sep, next_non_sep, old_next_non_sep, dummy
    character(len=len(char_array(1))) :: value

    len_arr = size(char_array)
    if (fu_fails(len(separator) == 1, 'len(separator) must be 1', 'split_string_char')) return
    
    num_tokens = 0
    next_non_sep = fu_next_non_sep(string, separator, 1)

    do while (next_non_sep > 0)
      next_sep = fu_next_sep(string, separator, next_non_sep)
      value = string(next_non_sep:next_sep-1)
      num_tokens = num_tokens + 1
      if (fu_fails(num_tokens <= len_arr, 'Too many tokens', 'split_string_char')) return
      char_array(num_tokens) = value
      next_non_sep = fu_next_non_sep(string, separator, next_sep)
    end do

  end subroutine split_string_char

  !************************************************************************************
  
  ! The following two functions are used by the split_string
  ! subroutines. 
  !
  integer function fu_next_sep(str, separator, ind_start) result(ind)
    implicit none
    character(len=*), intent(in) :: str, separator
    integer, intent(in) :: ind_start

    integer :: i

    ind = 0
    do i = ind_start, len(str)
      if (str(i:i) == separator) then
        ind = i
        return
      end if
    end do
    ind = len(str)+1
  end function fu_next_sep

  integer function fu_next_non_sep(str, separator, ind_start) result(ind)
    implicit none
    character(len=*), intent(in) :: str, separator
    integer, intent(in) :: ind_start

    integer :: i

    ind = 0
    do i = ind_start, len(str)
      if (str(i:i) /= separator) then
        ind = i
        return
      end if
    end do
  end function fu_next_non_sep

  !************************************************************************************

  FUNCTION fu_trim_grads_hat(dataname,ctlname) result(outdataname)
    ! replaces a comon path with a hat "^" symbol
    ! if dir_slash == '/' (i.e. UNIXes)
    ! For Windows, restores a full path if it is relative for ctl file
    !
    ! MAS 17.04.2018. way to handle hat is found for Windows: full(!) switch to UNIX-slash
    ! everywhere. Any presence of Windows slash anywhere in ctl or gs files breaks everything
    ! With this, can get rid of full path

    IMPLICIT NONE
    !
    ! The return value of this function:
    CHARACTER (LEN=fnlen) :: outdataname
    ! Input parameters
    CHARACTER (LEN=*), INTENT(in) :: dataname,ctlname
    ! Local variables
    integer :: iSlash, iTmp
    CHARACTER (LEN=fnlen) :: ctlNameTmp, dataNameTmp

    ctlNameTmp = ctlname
    call check_dir_slash(ctlNameTmp)
    dataNameTmp = dataName
    call check_dir_slash(dataNameTmp)
    
    iSlash=index(ctlnameTmp,dir_slash,.True.) !last slash in ctl file name

!    if (dir_slash == '/') then ! Should be removed when/if writing 
                               ! absolute paths is acceped to be evil everywhere
      if (index(datanameTmp, ctlnameTmp(1:iSlash)) == 1) then !common path
        outdataname = "^" + datanameTmp(iSlash+1 : len_trim(datanameTmp)) !fname
      else
        outdataname = datanameTmp
      endif
!    else
!      !
!      ! Windows: need absolute path. Neither hat nor relative one go
!      !
!      if(index(datanameTmp,':') == 2)then
!        outdataname = datanameTmp
!      else
!        !call get_working_directory(ctlNameTmp)  ! temporary use of the variable
!        ctlNameTmp = fu_current_directory()
!        iTmp = 0
!        do while(index(dataNameTmp(iTmp+1 :),'..') == 1)  ! climb up the tree if needed
!          iTmp = index(dataNameTmp,dir_slash)
!          ctlNameTmp(index(ctlNameTmp,dir_slash,.true.) : ) = ' '
!        end do
!        outdataname = ctlNameTmp + dir_slash + datanameTmp(iTmp+1 : len_trim(dataNameTmp))
!      endif
!    endif
  end FUNCTION fu_trim_grads_hat


!*****************************************************

  FUNCTION fu_extend_grads_hat(dataname,ctlname) result(outdataname)
    ! if dataname starts with hat "^" returns dataname with path of
    ! ctlname

    IMPLICIT NONE
    CHARACTER (LEN=fnlen) :: outdataname
    CHARACTER (LEN=*), INTENT(in) :: dataname,ctlname

    outdataname = fu_extend_grads_hat_dir(dataname,fu_dirname(ctlname))

   end FUNCTION fu_extend_grads_hat

!*****************************************************

  FUNCTION fu_extend_grads_hat_dir(dataname,dirname) result(outdataname)
    ! if dataname starts with hat "^" returns dataname with path of
    ! dirname

    IMPLICIT NONE
    CHARACTER (LEN=fnlen) :: outdataname
    CHARACTER (LEN=*), INTENT(in) :: dataname,dirname
    integer :: iTmp, indHat

    indHat = index(adjustl(dataname), "^")
!    call msg('indHat', indHat)
   
    if (indHat > 0) then ! relative path
      if(len_trim(dirname) > 0)then
        outdataname = adjustl(dirname) + dir_slash + dataname(indHat+1:len(dataname))
      else
        outdataname = dataname(indHat+1:len(dataname))
      endif
    else
      outdataname = adjustl(dataname)
    endif
    if(indHat > 1) outdataname = adjustl(dataname(1:indHat-1)) // outdataname
!    call msg('extend_grads_hat: initial=' + dataname + ', dir=' + dirname + ', final=' + outdataname)
  end FUNCTION fu_extend_grads_hat_dir

  
!*****************************************************

   logical function fu_if_ends_with(string, substring) result(ends_with)
     implicit none
     character(len=*), intent(in) :: string, substring
     
     ends_with = (index(string, substring) > 0 &
             & .and. index(string, substring) == len_trim(string)-len_trim(substring) + 1)
   end function fu_if_ends_with

   
!*******************************************************************
   
   subroutine replace_string(string, replace_what, replace_with, case_sensitive)
     ! Replace in-place all occurances of replace_what with replace_with. By default case
     ! sensitive.
     implicit none
     character(len=*), intent(inout) :: string
     character(len=*), intent(in) :: replace_what, replace_with
     logical, intent(in), optional :: case_sensitive
     
     character(len=len(string)) :: work_string
     integer :: ind_substr_end, ind_repl_end, ind_start
     logical :: replace, case_sensitive_
     
     if (present(case_sensitive)) then
       case_sensitive_ = case_sensitive
     else
       case_sensitive_ = .true.
     end if
     
     replace = .true. ! but exit below if not found 
     
     do while (replace)
       if (case_sensitive_) then
         ind_start = index(string, replace_what)
       else
         ind_start = index(fu_str_l_case(string), fu_str_l_case(replace_what))
       end if
       replace = ind_start > 0
       if (.not. replace) exit
       ! string = what_ever
       ! replace_what = what
       ! replace_with = who
       ! ind_start = 1, ind_repl_end = 3, ind_substr_end = 4
       
       if (ind_start > 1) work_string(1:ind_start-1) = string(1:ind_start-1)
       ind_repl_end = ind_start + len(replace_with) - 1
       if (fu_fails(ind_repl_end <= len(string), 'Cannot replace string', 'replace_string')) return
       ! work_string: what_ever -> whot_ever
       work_string(ind_start:ind_repl_end) = replace_with
       ind_substr_end = ind_start + len(replace_what) - 1
       ! work_string: whot_ever -> who_ever
       work_string(ind_repl_end+1:) = string(ind_substr_end+1:)
       
       string = work_string
     end do
     
   end subroutine replace_string

!******************************************************************************************
   
  character(len=fnlen) function fu_dirname(path) result(dirname)
    ! The directory component from a path.
    implicit none
    character(len=*), intent(in) :: path
     
    integer :: ind_last_slash

    ind_last_slash = max(index(path, '/', back=.true.),index(path, '\', back=.true.))
    if (ind_last_slash <= 1) then
      dirname = ''
      return
    end if
     
    dirname = path(1:ind_last_slash - 1)
     
  end function fu_dirname

  !************************************************************************************
  
  character(len=fnlen) function fu_basename(path) result(basename)
    ! The filename component from a path.
    implicit none
    character(len=*), intent(in) :: path
     
    integer :: ind_last_slash

    ind_last_slash = index(path, dir_slash, back=.true.)
    if (ind_last_slash <= 1) then
      basename = path
      return
    end if
     
    basename = path(ind_last_slash+1:)
     
  end function fu_basename


   
   !****************************************************************************************
   
  subroutine linear_regression(x, y, nP, slope, intercept, r_value, p_value, slope_err, stderr, nMissing)
    !
    ! Calculates the basic parameters of linear regression of two datasets. The algorithm is taken
    ! from Python scipy.stats module.
    ! We also exclude real_missing using direct .eps. function code rather than a call to speedup
    ! 
    implicit none
    
    ! Imported parameters
    real, dimension(:), pointer :: x, y
    real, intent(out) :: slope, intercept, r_value, p_value, slope_err, stderr
    integer, intent(in) :: nP
    integer, intent(out) :: nMissing
    
    ! Local parameters
    real :: xmean, ymean, ssxm, ssym, ssxym, r_den, df, t
    integer :: ii, count
    
    if(fu_fails(nP <= size(x) .and. nP <= size(y),'nP > size(x or y):' + &
              & fu_str(nP) + ',' + fu_str(size(x)) + ',' + fu_str(size(y)),'linear_regression'))return
    nMissing = 0
    xmean = 0.
    ymean = 0.
    ssxm = 0.
    ssym = 0.
    ssxym = 0.
    do ii = 1, nP
      if(abs(x(ii)-real_missing) > 1. .and. abs(y(ii)-real_missing) > 1.)then
        xmean = xmean + x(ii)
        ymean = ymean + y(ii)
        ssxm =  ssxm + x(ii) * x(ii)
        ssym = ssym + y(ii) * y(ii)
        ssxym = ssxym + x(ii) * y(ii)
      else
        nMissing = nMissing + 1
      endif
    end do
    if(fu_fails(nMissing < nP-1, 'Too many missing points:' + fu_str(nMissing), 'linear_regression'))return
    xmean = xmean / (nP - nMissing)
    ymean = ymean / (nP - nMissing)
    ssxm =  ssxm / (nP - nMissing)
    ssym = ssym / (nP - nMissing)
    ssxym = ssxym / (nP - nMissing)

    slope = (ssxym - xmean * ymean) / (ssxm - xmean * xmean)
    intercept = ymean - slope * xmean
    r_den = sqrt((ssxm-xmean*xmean)*(ssym-ymean*ymean))
    if(r_den == 0.0)then
      r_value = 0.0
    else
      r_value = max(-1., min(1., (ssxym - xmean * ymean) / r_den))   ! test for numerical error propagation
    endif
    df = nP - nMissing - 2.
    t = r_value * sqrt(df/((1.0000001 - r_value) * (1.0000001 + r_value)))
    p_value = fu_studnt(abs(t), df) * 2.
!    slope = ssxym / ssxm
!    intercept = ymean - slope * xmean
!    slope_err = sqrt((1. - r_value * r_value) * ssym / ssxm / df)
    slope_err = slope / r_value * sqrt((1. - r_value * r_value) / (nP-nMissing-2))
    stderr = ssym + slope * slope * ssxm + intercept * intercept - &
           & 2.* (slope * ssxym + intercept * ymean - slope * intercept * xmean)

  end subroutine linear_regression

  
  !*****************************************************************************
  
  real function fu_lognorm_density(fD, fDmed, fSigma)            
    !
    ! Returns the value of lognormal distribution density
    !
    implicit none
    
    ! Imported parameters
    real, intent(in) :: fD, fDmed, fSigma
    
    fu_lognorm_density = 1D0 / (fD * log(fSigma)*dsqrt(2.0D0*Pi)) * &
                       & exp(-((log(fD)-log(fDmed))/log(fSigma))**2 / 2.0D0)

  end function fu_lognorm_density
                                     

  !****************************************************************************
  
  REAL FUNCTION fu_studnt (t, doff) RESULT(fn_val)
    !
    !  ALGORITHM AS 27  APPL. STATIST. VOL.19, NO.1
    !  Calculate the upper tail area under Student's t-distribution
    !  Translated from Algol by Alan Miller
    !
    IMPLICIT NONE
    REAL, INTENT(IN)      :: t
    REAL, INTENT(IN)      :: doff

    !     Local variables
    REAL     :: v, x, tt
    REAL, PARAMETER  :: a1 = 0.09979441, a2 = -0.581821, a3 = 1.390993,  &
                        a4 = -1.222452, a5 = 2.151185
    REAL, PARAMETER  :: b1 = 5.537409, b2 = 11.42343
    REAL, PARAMETER  :: c1 = 0.04431742, c2 = -0.2206018, c3 = -0.03317253,  &
                        c4 = 5.679969, c5 = -12.96519
    REAL, PARAMETER  :: d1 = 5.166733, d2 = 13.49862
    REAL, PARAMETER  :: e1 = 0.009694901, e2 = -0.1408854, e3 = 1.88993,  &
                        e4 = -12.75532, e5 = 25.77532
    REAL, PARAMETER  :: f1 = 4.233736, f2 = 14.3963
    REAL, PARAMETER  :: g1 = -9.187228E-5, g2 = 0.03789901, g3 = -1.280346,  &
                        g4 = 9.249528, g5 = -19.08115
    REAL, PARAMETER  :: h1 = 2.777816, h2 = 16.46132
    REAL, PARAMETER  :: i1 = 5.79602E-4, i2 = -0.02763334, i3 = 0.4517029,  &
                        i4 = -2.657697, i5 = 5.127212
    REAL, PARAMETER  :: j1 = 0.5657187, j2 = 21.83269

    !     Check that number of degrees of freedom > 4.
    IF (fu_fails(doff > 4.0,'** Error in Student: degrees of freedom <= 4','studnt'))RETURN

    !   Evaluate series.
    v = 1. / doff
    tt = ABS(t)
    x = 0.5 * (1. +   &
        tt*(((a1 + v*(a2 + v*(a3 + v*(a4 + v*a5)))) / (1. - v*(b1 - v*b2))) +  &
        tt*(((c1 + v*(c2 + v*(c3 + v*(c4 + v*c5)))) / (1. - v*(d1 - v*d2))) +  &
        tt*(((e1 + v*(e2 + v*(e3 + v*(e4 + v*e5)))) / (1. - v*(f1 - v*f2))) +  &
        tt*(((g1 + v*(g2 + v*(g3 + v*(g4 + v*g5)))) / (1. - v*(h1 - v*h2))) +  &
        tt*((i1 + v*(i2 + v*(i3 + v*(i4 + v*i5)))) / (1. - v*(j1 - v*j2))) ))))) ** (-8)
    IF (t >= 0.) THEN
      fn_val = x
    ELSE
      fn_val = 1. - x
    END IF

  END FUNCTION fu_studnt

  
  !****************************************************************************************
  
  subroutine trim_precision(x, factor, missval)
     ! Trims the precision of vector x 
     real (kind=4), dimension(:), intent(inout) :: x
     real (kind=4), optional, intent(in) :: missval
     real, intent(in) :: factor
     real (kind=4)  :: fTmp
     integer (kind=4) :: nbits, bitmask,i !! Must be the same size as X and missval
     

#ifdef NEW_TRIM_PRECISION
     ! FIXME Causes GNU compiler crash on cray
     ! The resulting field is accurate within 0 .. 0.5/factor
     !WARNING!!!  Heavilon_server.silam.mod.f90s relies on 32-bit real 
     !For 64-bit a separate version should be made
       
     ! Get number of bits to keep 
     i = ceiling(factor)
     bitmask= TRANSFER(X'FFFFFFFF',nbits)
     
     do nbits = 1,23
       if (rshift(i,nbits) == 0) exit
     enddo
     bitmask=lshift(bitmask, 23-nbits)

     !NOTE this method is simple, but scales field by 1 to 1 + 0.5/factor
     !          with mean scaling 1 + 0.25./factor
     !
     ! FIXME Could save one more bit and get rid of the offset by rounding instead of cutting
     !
     if (present(missval)) then
       where (x /= missval) 
         x = transfer(IAND(bitmask,transfer(x,1,size(x))), fTmp, size(x))
       end where
     else
       x = transfer(IAND(bitmask,transfer(x,1,size(x))), fTmp, size(x))
     endif
#else
     if (radix(factor) .ne. 2.) then 
       call msg("Some starnge radix=",radix(factor))
       call msg_warning("Not cutting precision","trim_precision")
       return
     endif
     if (present(missval)) then
       where (abs(x) < 1e-37) 
         x = 0. !Kill denormalized to avoid Floating-point exceptions below
       else where (x /= missval) 
         x =  nint(factor * fraction(x))*(2.**exponent(x))/factor
       end where
     else
       where (abs(x) < 1e-37) 
         x = 0. !Denormalized 
       elsewhere
         x =  nint(factor * fraction(x))*(2.**exponent(x))/factor
       endwhere
     endif
#endif

  end subroutine trim_precision

 !*******************************************************

  real function fu_trim_cell_fcoord(f)
    ! adjusts floating-point coordinate of cell in a grid
    ! to avoid incomplete slab in advection
      implicit none
      real, intent(in) :: f
      real :: intf, fracf

      real, parameter :: fZcTrimin = 1.0/6  ! Maximum CM for trapezoid slab
      character(len = *), parameter :: sub_name = 'fu_trim_cell_fcoord'
    ! code

      intf = nint(f)
      fu_trim_cell_fcoord = max(intf - fZcTrimin,min(intf + fZcTrimin, f) )
  end function fu_trim_cell_fcoord

  !***************************************************************************************
  
  SUBROUTINE set_environment_variable(varname,newval,iStat)
        implicit none
        character(len=*), intent(in) :: varname, newval
        character(len=fnlen):: varnamelocal, newvallocal
        integer, intent(out) :: iStat

        varnamelocal = trim(varname)//ACHAR(0)
        newvallocal = trim(newval)//ACHAR(0)        
        
        iStat = setenv(varnamelocal,newvallocal, 1)
           


  end SUBROUTINE set_environment_variable

  !***************************************************************************************
  

  ! Write buffering. The data are normally written into file one field at time. Especially
  ! in MPI runs, each write becomes small, and performance can be increased by first
  ! aggregating data from multiple fields. For serial runs buffering doesn't necessary
  ! improve the performance, however.
  ! 
  ! Current default buffer size is 50 fields, but if the number of MPI tasks is high,
  ! much more buffering can be useful.
  !
  !  SINCE 2020-07-23 also used for MPI NetCDF buffered read/write
  ! 
  
  integer function fu_get_default_mpi_buf_size() result(buf_size_flds)
    implicit none
    integer, parameter :: buf_size_flds_def = 500
    integer :: stat
    character(len=*), parameter :: sub_name = 'fu_get_default_mpi_buf_size'
    character(len=clen) :: strval

    call get_environment_variable('GRADSIO_BUF_SIZE', strval)
    if (error) return
    if (strval /= '') then
      read(unit=strval, fmt=*, iostat=stat) buf_size_flds
      if (fu_fails(stat == 0, 'Failed to parse GRADSIO_BUF_SIZE: ' + strval, sub_name)) return
      if (fu_fails(buf_size_flds < 1000000, 'Suspicious GRADSIO_BUF_SIZE', sub_name)) return
      if (buf_size_flds < 2) call msg_warning('GRADSIO_BUF_SIZE: ' + fu_str(buf_size_flds), sub_name)
      ! size < 1 used for disabling buffer (above)
    else
      buf_size_flds = buf_size_flds_def
    end if
    
  end function fu_get_default_mpi_buf_size

  
  !***************************************************************************************
  
  SUBROUTINE string_to_fnames(string, fnames, iNbrOfNames, ifNonZeroRequired)

    ! Description:
    ! Creates a list of filenames from a given single string. Given
    ! string may contain wildcard character '*' or '?'. In this case this
    ! routine returns the real list of files found through that
    ! filter. 
    !
    ! Method:
    ! If a wildcard charatcter is found then the list of files is
    ! created on shell level to file filelist.silam
    !
    ! ISHELL to execute a command on shell level.
    ! shell command ls, cat, wc
    !
    IMPLICIT NONE

    ! Imported parameters with intent IN:
    CHARACTER (LEN=*), INTENT(in) :: string
    integer, intent(out) :: iNbrOfNames
    logical, intent(in), optional :: ifNonZeroRequired

    ! Imported parameters with intent INOUT or POINTER:
    type(fnlen_str), DIMENSION(:), allocatable, intent(inout) :: fnames

    ! Local declarations:
    LOGICAL :: filelist_exists
    CHARACTER (LEN=fnlen) :: cmd, file_list_name
    character (len=5) :: chTmp
    real :: fTmp
    INTEGER :: i, FUnit, j, iFile, FUnitTmp
    character(len=1) :: chTest
    logical :: ifNonZeroRequiredLocal
    
    if(fu_fails(len_trim(string) > 0, 'Zero-length string given','string_to_fnames'))return
    
    if (present(ifNonZeroRequired)) then
      ifNonZeroRequiredLocal = ifNonZeroRequired
    else
      ifNonZeroRequiredLocal = .false.
    end if
    
    if(allocated(fnames)) then
      deallocate(fnames,stat=i)
      if (fu_fails(i == 0, 'Failed to deallocate string pointers','string_to_fnames'))RETURN
    endif
    !----------------------------------------
    !
    ! 1. No wildcard, just file name
    !
    IF (INDEX(string, '*') == 0 .and. INDEX(string, '?') == 0) THEN
!      CALL set_stringlist_size(fnames, 1)
      allocate(fnames(1),stat=i)
      IF (fu_fails(i == 0, 'Failed to allocate string pointers','string_to_fnames'))RETURN
      iNbrOfNames = 1
      fnames(1)%s = TRIM(string)
      RETURN
    END IF

    !
    ! 2. Wildcard found
    !
    ! Generate the filename
    write(unit=file_list_name, fmt='(A,I0)') 'filelist.silam_', fu_pid()

    !
    ! 2.2. Filter the names using ls:

    call list_files(string, file_list_name)
    IF (error) RETURN

    !
    ! 2.3. How many files found?

    iNbrOfNames = fu_number_of_lines_in_file(file_list_name)
    IF (error) RETURN

    !
    ! 2.4. Read the list-file (first, find free unit)

    funit = fu_next_free_unit()
    if (error) return
    OPEN (unit=FUnit, file=file_list_name, status='old', action='read')

    if(iNbrOfNames > 0)then
      !
      ! Something has been found
      !
      allocate(fnames(iNbrOfNames),stat=i)
      IF(fu_fails(i == 0, 'Failed to allocate string pointers.2','string_to_fnames'))RETURN
      !
      ! If non-zero file size is required, we shall open it and read 1 byte
      !
      if(ifNonZeroRequiredLocal)then
        funitTmp = fu_next_free_unit()
      endif

      iFile = 1
      do i=1,iNbrOfNames
        READ(unit=FUnit, fmt='(A)') fnames(iFile)%s
        fnames(iFile)%s = ADJUSTL(fnames(iFile)%s)
        if(ifNonZeroRequiredLocal)then
          call open_binary(unit=funitTmp, file=fnames(iFile)%s, recl=1, action='read')
          read(FUnitTmp,rec=1,iostat = j)chTest
          close(FUnitTmp)
          if(j==0) iFile=iFile+1
        else
          iFile = ifile + 1
        endif
      END DO  ! iNbrOfNames

      iNbrOfNames = iFile-1
      
      do i=iNbrOfNames+1, size(fnames)
        fnames(i)%s = ''
      end do

    endif   ! non-empty list file

    CLOSE(FUnit,status='DELETE')
    
  END SUBROUTINE string_to_fnames

  
  !************************************************************************************

  INTEGER FUNCTION fu_number_of_lines_in_file(fname)
    
    ! Description:
    ! Returns the number of lines in a file. 
    ! 
    ! Method:
    ! Just read the file till the end
    !
    ! All units: SI
    !
    ! Language: ANSI Fortran 90
    !
    ! Author: M.Sofiev, FMI
    ! 
    IMPLICIT NONE
    !
    ! Imported parameters with intent(in):
    CHARACTER (LEN=*), INTENT(in) :: fname

    ! Local declarations:
    LOGICAL :: file_exists
    integer :: counter, FUnit, status

    ! --------------------------------------
    !
    ! 1. Check the file
    !    --------------

    INQUIRE(file=fname, exist=file_exists)

    IF (.NOT.file_exists) THEN
      CALL set_error(fu_connect_strings('did not find file', fname),&
	               & 'fu_number_of_lines_in_file')
      RETURN
    END IF

    ! --------------------------------------
    !
    ! 2. Read the file till the end
    !    ------------------------
    funit = fu_next_free_unit()
    open(unit=FUnit,file=fname,status='old',action='read')
    counter = -1
    status = 0
    do while(status == 0)
      read(FUnit,fmt=*,iostat = status)
      counter = counter +1
    end do
    
    fu_number_of_lines_in_file = counter

    CLOSE(FUnit)

  END FUNCTION fu_number_of_lines_in_file

  !************************************************************************************

  subroutine open_binary(unit, file, recl, status, action, access, convert, iostat)
    ! 
    ! Open a binary file in machine-independent way.  The actual open
    ! statement is in md, this subroutine gives only some default values for the optional 
    ! arguments.
    ! 
    implicit none 
    integer, intent(in) :: unit, recl 
    character(len=*), intent(in) :: file
    character(len=*), intent(in), optional :: status, convert, action, access 
    integer, intent(out), optional :: iostat
    
    character(len=16) :: convert_, status_, action_, access_
    integer :: iostat_

    if (present(convert)) then
      convert_ = convert
    else
      convert_ = 'native'
    end if
    if (present(status)) then
      status_ = status
    else
      status_ = 'unknown'
    end if
    if (present(action)) then
      action_ = action
    else
      action_ = 'readwrite'
    end if
    if (present(access)) then
      ! Usually we use direct access, but this is inconvenient for the buffered grads
      ! writes, where we'll use access='stream'. This is F2003 standard. Some older Intel
      ! compilers might get the same effect by setting form='binary' and access='sequential'.
      access_ = access
    else
      access_ = 'direct'
      if (fu_fails(recl > 0, 'Bad record lenght for direc access file', 'open_binary')) return
    end if

    select case(trim(fu_str_l_case(convert_)))
      ! Since convert is not a Fortran standard, we define the canonical values here -
      ! they happen to match the current compilers.
    case ('native', 'big_endian', 'little_endian')
      continue
    case default
      call set_error('Illegal value for convert: ' // trim(convert), 'open_binary')
      if (present(iostat)) iostat = -1
      return
    end select
    call open_binary_md(unit, file, recl, status_, action_, convert_, access_, iostat_)
    if (present(iostat)) then
      iostat = iostat_
    else
      if (iostat_ /= 0) then
        call msg('iostat:', iostat_)
        call set_error('Failed to open:'//trim(file), 'open_binary')
      end if
    end if
    
  end subroutine open_binary

  !************************************************************************************
  !
  ! The following grads-related subroutines would be more logical in the grads_io module
  ! but are not there due to some circular dependencies.

  SUBROUTINE open_grads_binary_o(fname, iUnit, NbrOfPoints)
    !
    ! Opens a bimary file for the direct access with an appropriate recordlength
    !
    ! Current code owner Mikhail Sofiev
    !
    IMPLICIT NONE
    !
    ! Imported parameters with intent IN:
    CHARACTER (LEN=*), INTENT(in) :: fname ! File name to open
    INTEGER, INTENT(in) :: iUnit ! Unit to be attributed to the file
    INTEGER :: NbrOfPoints ! Number of values in per record
    integer :: status

    ! Local declarations
    INTEGER :: NbrOfBytes

    IF(iUnit <= 0 .or. iUnit > 100)THEN
      CALL set_error('starnge file unit','fu_open_gradsfile_o_md')
      RETURN
    END IF
    IF(NbrOfPoints < 1)THEN
      CALL set_error('starnge number of elements per record','fu_open_gradsfile_o_md')
      RETURN
    END IF
    IF(fname == '')THEN
      CALL set_error('starnge file name','fu_open_gradsfile_o_md')
      RETURN
    END IF

    NbrOfBytes = NbrOfPoints * sizeof_sp / 8 ! from md, sizeof_sp in bits

    !
    ! GrADS file is a direct-access file, so it is not killed 
    ! while opening and writing. So, we have to do it explicitly.
    !
    OPEN(iUnit,file=fname, iostat = status)
    if(status /= 0)then
      call set_error('Failed to open GrADS file','open_gradsfile_md')
      return
    else
      close(iUnit,status='delete') ! Kill the file
      call open_binary(unit=iUnit, file=fname, recl=nbrOfBytes, action='write')
    endif
  
  END SUBROUTINE open_grads_binary_o


  ! ***************************************************************

  SUBROUTINE open_grads_binary_i(fname, iUnit, NbrOfPoints, ifBigEndian)
    !
    ! Opens a bimary file for the direct access with an appropriate recordlength
    !
    ! Current code owner Mikhail Sofiev
    !
    IMPLICIT NONE
    !
    ! Imported parameters with intent IN:
    CHARACTER (LEN=*), INTENT(in) :: fname ! File name to open
    INTEGER, INTENT(in) :: iUnit ! Unit to be attributed to the file
    INTEGER :: NbrOfPoints ! Number of values in per record
    integer :: status
    logical, intent(in) :: ifBigEndian

    ! Local declarations
    INTEGER :: NbrOfBytes

    IF(iUnit <= 0 .or. iUnit > 100)THEN
      CALL set_error('starnge file unit','fu_open_gradsfile_i_md')
      RETURN
    END IF
    IF(NbrOfPoints < 1)THEN
      CALL set_error('starnge number of elements per record','fu_open_gradsfile_i_md')
      RETURN
    END IF
    IF(fname == '')THEN
      CALL set_error('starnge file name','fu_open_gradsfile_i_md')
      RETURN
    END IF

    NbrOfBytes = sizeof_sp / 8 ! from md  To protect from overflow
    NbrOfBytes = NbrOfPoints * NbrOfBytes ! from md 

    if (.not. NbrOfBytes > 0) then
      call msg('Failed to open GrADS file binary: '//trim(fname))
      
      call set_error("Looks like overflow when opening grads binary", "open_grads_binary_i")
      return
    endif
    
    ! GrADS file is a direct-access file, so it is not killed 
    ! while opening and writing. So, we have to do it explicitly.
    !
    if(ifBigEndian)then
      call open_binary(unit=iUnit, file=fname, recl=nbrOfBytes, action='read', status='old', iostat=status, &
                     & convert='big_endian')
    else
      call open_binary(unit=iUnit, file=fname, recl=nbrOfBytes, action='read', status='old', iostat=status)
    endif
    if(status /= 0)then
      call set_error('Failed to open GrADS file','open_gradsfile_i_md')
    endif
  
  END SUBROUTINE open_grads_binary_i


  ! ****************************************************************

  subroutine write_grads_field(grid_data, npoints, irec, iunit)
    !
    ! Write a field to a direct access file in a machine-independent way.
    !
    use globals, only : r4k
    implicit none
    real, dimension(:), intent(in) :: grid_data
    integer, intent(in) :: npoints, irec, iunit

    integer :: stat
    real(kind=r4k), dimension(:), pointer, save :: data_sp =>null()

    if (fu_fails(npoints <= size(grid_data) .and. npoints > 0, 'Bad nPoints', 'write_grads_field')) return

    if (kind(grid_data) == r4k) then
      ! Normal situation, no conversion.
      !
      write(iunit, rec=irec) grid_data(1:npoints)
    else 
      ! Need conversion, grid_data possibly big.
      if (associated(data_sp)) then
         if (npoints > size(data_sp)) then
            call msg("reallocateing write_grads_field array from, to", data_sp, npoints)
            deallocate(data_sp)
            allocate(data_sp(npoints), stat=stat)
            if (fu_fails(stat == 0, 'Allocate failed', 'write_grads_field')) return
         endif
      else
         allocate(data_sp(npoints), stat=stat)
         if (fu_fails(stat == 0, 'Allocate failed', 'write_grads_field')) return
      endif
      data_sp = grid_data(1:npoints)
      write(iunit, rec=irec) data_sp
    end if
      
  end subroutine write_grads_field

  !************************************************************************************************
  
  subroutine write_grads_fieldset(grid_data_set, npoints, num_flds, iunit)
    !
    ! Write several fields to a direct access file. Must be single precision array. File
    ! must be opened for stream (or maybe sequential) io.
    !
    use globals, only : r4k
    implicit none
    real(r4k), dimension(:,:), intent(in) :: grid_data_set ! ind_fld, ind_point
    integer, intent(in) :: npoints, iunit, num_flds
    integer :: stat

    if (fu_fails(npoints <= size(grid_data_set, 1) .and. npoints > 0, 'Bad nPoints', 'write_grads_field')) return
    if (fu_fails(num_flds <= size(grid_data_set, 2) .and. num_flds > 0, 'Bad num_flds', 'write_grads_field')) return
    
    write(iunit) grid_data_set
      
  end subroutine write_grads_fieldset

  !************************************************************************************

  subroutine get_chunk(count_total, num_chunks, displs, sizes, if_zero_based)
    ! 
    ! Distribute count_total elements into num_chunks contiguous subsets. Example:
    ! [1,2,3,4,5] -> [1,2] [3,4], [5] (count_total = 5, num_chunks = 3).  On output,
    ! displs gives the starting indices (1 or 0 based) in the original array, sizes gives
    ! the chunk size.
    ! 
    implicit none
    integer, intent(in) :: count_total, num_chunks
    integer, dimension(:), intent(out):: displs, sizes
    logical :: if_zero_based

    character(len=*), parameter :: sub_name = 'get_chunk'
    integer :: chunk_size, ind_chunk, remainder, begin, prev_size

    if (fu_fails(num_chunks <= count_total, 'Too many chunks for count', sub_name)) return
    if (fu_fails(size(displs) >= num_chunks, 'displs too small', sub_name)) return
    if (fu_fails(size(sizes) >= num_chunks, 'sizes too small', sub_name)) return

    chunk_size = count_total / num_chunks
    remainder = count_total - num_chunks*chunk_size
    prev_size = 0
    begin = 1
    do ind_chunk = 1, num_chunks
      displs(ind_chunk) = begin + prev_size
      begin = displs(ind_chunk)
      if (ind_chunk <= remainder) then
        sizes(ind_chunk) = chunk_size + 1
      else
        sizes(ind_chunk) = chunk_size
      end if
      prev_size = sizes(ind_chunk)
    end do

    if (if_zero_based) displs(1:num_chunks) = displs(1:num_chunks) - 1

  end subroutine get_chunk

  subroutine set_coastal_value(field, pLandFr, iMeteo, nx_meteo, ny_meteo, &
       & value_in_cell, default_value, if_divide)

    real, dimension(:), pointer, intent(in) :: field, pLandFr
    integer, intent(in) :: iMeteo, nx_meteo, ny_meteo
    real, intent(in) :: default_value
    logical, intent(in) :: if_divide
    real, intent(out) :: value_in_cell
    real :: fTmp
    integer :: iCount, jTmp, iTmp, iyMeteo, ixMeteo, iMeteoTmp

    ! Sets the value of a field in a coastal cell to equal the average of the
    ! mostly land-containing neighbors. In case no mostly land-containing neighbor
    ! is found, the default_value is used. The values can also be divided by the
    ! land area (controlled by if_divide).

    if (pLandFr(iMeteo) < 0.85) then !if (field(iMeteo) < 1e-4) then
      iCount = 0
      fTmp = 0.0
      do jTmp = max(iyMeteo-1,1), min(iyMeteo+1, ny_meteo)
        do iTmp =  max(ixMeteo-1,1), min(ixMeteo+1, nx_meteo)
          iMeteoTmp = iTmp + nx_meteo * (jTmp -1)
          if (iMeteo == iMeteoTmp) cycle
          if (pLandFr(iMeteoTmp) >= 0.85) then
            if (if_divide) then
              fTmp = fTmp + field(iMeteoTmp)/max(0.0, min(1.0, pLandFr(iMeteoTmp)))
            else
              fTmp = fTmp + field(iMeteoTmp)
            end if
            iCount = iCount + 1
          end if
        end do
      end do
      if (iCount > 0) then
        fTmp = fTmp/iCount
        if (fTmp > field(iMeteo)) then
          value_in_cell = fTmp
        else
          if (if_divide) then
            value_in_cell = field(iMeteo)/max(0.0, min(1.0, pLandFr(iMeteo)))
          else
            value_in_cell = field(iMeteo)
          end if
        endif
      else
        value_in_cell = default_value  ! virtually anything: this is a kind-of isolated island
      endif
    else
      if (if_divide) then
        value_in_cell = field(iMeteo)/max(0.0, min(1.0, pLandFr(iMeteo)))
      else
        value_in_cell = field(iMeteo)
      end if
    end if ! value ~0
  end subroutine set_coastal_value


  !!!!!!***********************************************************************************
  !!!!!!
  !!!!!! A bunch of subroutines needed for IS4FIRES.
  !!!!!!
  !!!!!!***********************************************************************************
  !!!!!
  !!!!!real function fu_Planck(Lambda,Temp)
  !!!!!  ! Planck function for wavelangth Lambda and black-body temperature Temp
  !!!!!  implicit none
  !!!!!  real, intent(in) :: Lambda, Temp
  !!!!!  real, parameter :: c1= 3.74E-16, &  ! W.m2
  !!!!!                   & c2 = 1.439E-2    ! m.K
  !!!!!  fu_Planck = c1 / (Lambda**5 * (exp(c2 / (Lambda*Temp))-1.0))
  !!!!!end function fu_Planck
  !!!!!
  !!!!!!===================================================================================
  !!!!!
  !!!!!real function fu_Dozier_fit_error(Trad_21, Trad_31, Tback_21, Tback_31, lambda1, lambda2, fract, Tfire)
  !!!!!  !
  !!!!!  ! Dozier algorithm fitting error: two wavelengths, radiative temperatures of fire and background are given
  !!!!!  ! find a fraction of the pixel burning and temperature of fire
  !!!!!  !
  !!!!!  implicit none
  !!!!!  ! Imported parameters
  !!!!!  real, intent(in) :: Trad_21, Trad_31, Tback_21, Tback_31, lambda1, lambda2, fract, Tfire
  !!!!!  ! Local variables
  !!!!!  real :: eq1, eq2
  !!!!!  eq1 = fu_Planck(lambda1, Trad_21) - (fract * fu_Planck(lambda1, Tf) + (1.-fract) * fu_Planck(lambda1, Tback_21))
  !!!!!  eq2 = fu_Planck(lambda2, Trad_31) - (fract * fu_Planck(lambda2, Tf) + (1.-fract) * fu_Planck(lambda2, Tback_31))
  !!!!!  fu_Dozier_fit_error = eq1**2+eq2**2
  !!!!!end function fu_Dozier_fit_error
  !!!!!
  !!!!!!=====================================================================================
  !!!!!
  !!!!!subroutine get_fire_features_Dozier(Trad_21, Trad_31, Tback_21, Tback_31, lambda1, lambda2, fract, Tfire)
  !!!!!  !
  !!!!!  ! Solver for the Dozier problem
  !!!!!  !
  !!!!!  implicit none
  !!!!!
  !!!!!  ! imported parameters
  !!!!!  real, intent(in) :: Trad_21, Trad_31, Tback_21, Tback_31, lambda1, lambda2
  !!!!!  real, intent(inout) fract, Tfire
  !!!!!
  !!!!!  ! local variables
  !!!!!  real, parameter :: TfireMin = 300., Tfire_max = 2000., fract_min = 0, fract_max = 1., &
  !!!!!                   & lambda1 = 4.0E-6, lambda2 = 11.0E-6  ! m
  !!!!!  real :: err, stepT, stepFract
  !!!!!  !
  !!!!!  ! Temperature of fire should differ from that of the background
  !!!!!  !
  !!!!!  if (Trad_21 - tBack_21 < 2.0 .or. Trad_31 - tBack_31 < 0.05)then
  !!!!!    !
  !!!!!    ! Fire temperature must not be used, neither its fraction
  !!!!!    !
  !!!!!    fraact = -1.
  !!!!!    Tfire = -1.
  !!!!!  else
  !!!!!    !
  !!!!!    ! Solve the two-Planck system of equations
  !!!!!    !
  !!!!!    Tfire = (Tfire_min + tFire_max) * 0.5
  !!!!!    fract = (fract_min + fract_max) * 0.5
  !!!!!    stepT = (tFire_max - Tfire_min) * 0.1
  !!!!!    stepFr = (fract_max - fract_min) * 0.1
  !!!!!    err = fu_Dozier_fit_error(Trad_21, Trad_31, Tback_21, Tback_31, lambda1, lambda2, fract, Tfire)
  !!!!!
  !!!!!    ??????
  !!!!!
  !!!!!
  !!!!!  endif  ! if temperature contrast is sufficient
  !!!!!
  !!!!!end subroutine get_fire_features_Dozier


  !================================ Sorting
  ! Borrowed from Leon Foks https://leonfoks.github.io/coretran/index.html
!!!!BSD 3-Clause License
!!!!
!!!!Copyright (c) 2019, Leon Foks
!!!!All rights reserved.
!!!!
!!!!Redistribution and use in source and binary forms, with or without
!!!!modification, are permitted provided that the following conditions are met:
!!!!
!!!!1. Redistributions of source code must retain the above copyright notice, this
!!!!   list of conditions and the following disclaimer.
!!!!
!!!!2. Redistributions in binary form must reproduce the above copyright notice,
!!!!   this list of conditions and the following disclaimer in the documentation
!!!!   and/or other materials provided with the distribution.
!!!!
!!!!3. Neither the name of the copyright holder nor the names of its
!!!!   contributors may be used to endorse or promote products derived from
!!!!   this software without specific prior written permission.
!!!!
!!!!THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
!!!!AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
!!!!IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
!!!!DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
!!!!FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
!!!!DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
!!!!SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
!!!!CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
!!!!OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
!!!!OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

  ! Only stable sorting

   subroutine argMergeSort_int(this,i)
     integer, dimension(:), intent(in) :: this
     integer, dimension(:), intent(inout) :: i !! Should be 1:N on input
     integer,allocatable :: aux(:)
     integer :: iLeft,iRight, istat

     iLeft=1;iRight=size(this)
     allocate(aux(iRight),stat=istat);
     if  (istat /= 0) then
       call msg('MergeSort:Auxiliary Array',istat)
       call set_error("Allocate failed", "argMergeSort_int")
     endif
     aux=0
     call argSplitSort_i1D(this,i,aux,iLeft,iRight)
     deallocate(aux)
   end subroutine argMergeSort_int
  !====================================================================!
  !====================================================================!
  recursive subroutine argSplitSort_i1D(this,indx,aux,iLeft,iRight)
  !====================================================================!
  integer, intent(in) :: this(:)
  integer :: indx(:)
  integer :: aux(:)
  integer :: iLeft,iRight
  integer :: iMid
  if (iRight <= iLeft) return
  if (iRight <= iLeft + 19) then
    call argInsertionSort_i1D(this,indx,iLeft,iRight)
    return
  end if
  iMid=iLeft+(iRight-iLeft)/2
  call argSplitSort_i1D(this,indx,aux,iLeft,iMid)
  call argSplitSort_i1D(this,indx,aux,iMid+1,iRight)
  call argMergeSorted_i1D(this,indx,aux,iLeft,iMid,iRight)
  end subroutine
  !====================================================================!
  subroutine argMergeSorted_i1D(this,indx,aux,iLeft,iMid,iRight)
  !====================================================================!
  integer :: this(:)
  integer :: indx(:)
  integer :: aux(:)
  integer :: iLeft,iMid,iRight
  integer :: i,j,k

  aux(iLeft:iRight)=indx(iLeft:iRight)
  i=iLeft;j=iMid+1
  do k=iLeft,iRight
    if (i > iMid) then
      indx(k)=aux(j);j=j+1
    elseif (j > iRight) then
      indx(k)=aux(i);i=i+1
    elseif (this(aux(j)) < this(aux(i))) then
      indx(k)=aux(j);j=j+1
    else
      indx(k)=aux(i);i=i+1
    end if
  enddo
  end subroutine
  !====================================================================!
  subroutine argInsertionSort_i1D(this,indx,iLeft,iRight)
  integer :: this(:)
  integer :: indx(:)
  integer :: iLeft,iRight
  integer :: i,j, tmp
  do i=iLeft,iRight
    inner: do j=i,iLeft+1,-1
      if (this(indx(j)) < this(indx(j-1)))then
        tmp = indx(j) ; indx(j) = indx(j-1) ; indx(j-1) = tmp;
        !call swap(indx(j),indx(j-1))
      else
        exit inner
      end if
    end do inner
  end do
  end subroutine argInsertionSort_i1D
  !*************************************************************
  !*************************************************************
  !*************************************************************
  ! Strings just as a mechanical reproduction of integer sorting

   subroutine argMergeSort_char(this,i)
     character(len=*), intent(in) :: this(:)
     integer, dimension(:), intent(inout) :: i !! Should be 1:N on input
     integer, allocatable :: aux(:)
     integer :: iLeft,iRight, istat

     iLeft=1;iRight=size(this)
     allocate(aux(iRight),stat=istat);
     if  (istat /= 0) then
       call msg('MergeSort:Auxiliary Array',istat)
       call set_error("Allocate failed", "argMergeSort_char")
     endif
     aux=0
     call argSplitSort_char(this,i,aux,iLeft,iRight)
     deallocate(aux)
   end subroutine argMergeSort_char
  !====================================================================!
  !====================================================================!
  recursive subroutine argSplitSort_char(this,indx,aux,iLeft,iRight)
  !====================================================================!
  character(len=*), intent(in) :: this(:)
  integer :: indx(:)
  integer :: aux(:)
  integer :: iLeft,iRight
  integer :: iMid
  if (iRight <= iLeft) return
  if (iRight <= iLeft + 19) then
    call argInsertionSort_char(this,indx,iLeft,iRight)
    return
  end if
  iMid=iLeft+(iRight-iLeft)/2
  call argSplitSort_char(this,indx,aux,iLeft,iMid)
  call argSplitSort_char(this,indx,aux,iMid+1,iRight)
  call argMergeSorted_char(this,indx,aux,iLeft,iMid,iRight)
  end subroutine
  !====================================================================!
  subroutine argMergeSorted_char(this,indx,aux,iLeft,iMid,iRight)
  !====================================================================!
  character(len=*), intent(in) :: this(:)
  integer :: indx(:)
  integer :: aux(:)
  integer :: iLeft,iMid,iRight
  integer :: i,j,k

  aux(iLeft:iRight)=indx(iLeft:iRight)
  i=iLeft;j=iMid+1
  do k=iLeft,iRight
    if (i > iMid) then
      indx(k)=aux(j);j=j+1
    elseif (j > iRight) then
      indx(k)=aux(i);i=i+1
    elseif (this(aux(j)) < this(aux(i))) then
      indx(k)=aux(j);j=j+1
    else
      indx(k)=aux(i);i=i+1
    end if
  enddo
  end subroutine
  !====================================================================!
  subroutine argInsertionSort_char(this,indx,iLeft,iRight)
  character(len=*), intent(in) :: this(:)
  integer :: indx(:)
  integer :: iLeft,iRight
  integer :: i,j, tmp
  do i=iLeft,iRight
    inner: do j=i,iLeft+1,-1
      if (this(indx(j)) < this(indx(j-1)))then
        tmp = indx(j) ; indx(j) = indx(j-1) ; indx(j-1) = tmp;
        !call swap(indx(j),indx(j-1))
      else
        exit inner
      end if
    end do inner
  end do
  end subroutine argInsertionSort_char


  !*************************************************************
  !*************************************************************
  !*************************************************************
  ! Real: just as a mechanical reproduction of integer sorting

   subroutine argMergeSort_real(this,i)
     real, intent(in) :: this(:)
     integer, dimension(:), intent(inout) :: i !! Should be 1:N on input
     integer, allocatable :: aux(:)
     integer :: iLeft,iRight, istat

     iLeft=1;iRight=size(this)
     allocate(aux(iRight),stat=istat);
     if  (istat /= 0) then
       call msg('MergeSort:Auxiliary Array',istat)
       call set_error("Allocate failed", "argMergeSort_real")
     endif
     aux=0
     call argSplitSort_real(this,i,aux,iLeft,iRight)
     deallocate(aux)
   end subroutine argMergeSort_real
  !====================================================================!
  !====================================================================!
  recursive subroutine argSplitSort_real(this,indx,aux,iLeft,iRight)
  !====================================================================!
  real, intent(in) :: this(:)
  integer :: indx(:)
  integer :: aux(:)
  integer :: iLeft,iRight
  integer :: iMid
  if (iRight <= iLeft) return
  if (iRight <= iLeft + 19) then
    call argInsertionSort_real(this,indx,iLeft,iRight)
    return
  end if
  iMid=iLeft+(iRight-iLeft)/2
  call argSplitSort_real(this,indx,aux,iLeft,iMid)
  call argSplitSort_real(this,indx,aux,iMid+1,iRight)
  call argMergeSorted_real(this,indx,aux,iLeft,iMid,iRight)
  end subroutine
  !====================================================================!
  subroutine argMergeSorted_real(this,indx,aux,iLeft,iMid,iRight)
  !====================================================================!
  real, intent(in) :: this(:)
  integer :: indx(:)
  integer :: aux(:)
  integer :: iLeft,iMid,iRight
  integer :: i,j,k

  aux(iLeft:iRight)=indx(iLeft:iRight)
  i=iLeft;j=iMid+1
  do k=iLeft,iRight
    if (i > iMid) then
      indx(k)=aux(j);j=j+1
    elseif (j > iRight) then
      indx(k)=aux(i);i=i+1
    elseif (this(aux(j)) < this(aux(i))) then
      indx(k)=aux(j);j=j+1
    else
      indx(k)=aux(i);i=i+1
    end if
  enddo
  end subroutine
  !====================================================================!
  subroutine argInsertionSort_real(this,indx,iLeft,iRight)
  real, intent(in) :: this(:)
  integer :: indx(:)
  integer :: iLeft,iRight
  integer :: i,j, tmp
  do i=iLeft,iRight
    inner: do j=i,iLeft+1,-1
      if (this(indx(j)) < this(indx(j-1)))then
        tmp = indx(j) ; indx(j) = indx(j-1) ; indx(j-1) = tmp;
        !call swap(indx(j),indx(j-1))
      else
        exit inner
      end if
    end do inner
  end do
  end subroutine argInsertionSort_real




  !*************************************************************


  subroutine test_sort()

    integer, parameter ::  N = 160
    integer, dimension(N) :: ivals, indices
    integer :: i
    character(len=8), dimension(160), parameter :: strings= (/&
   & "oo5Aa7ve","taiPai8a","Ob0Yei5n","ahBi8uil","ohV6oovo","Oophah9I","oqu0Ieng","Ufee5lo1", &
   & "kee8Aeri","eish8Ahf","Shei7Viu","ohChi9da","Toaca5ea","xee2iCi4","du5aic1X","daeZahz7", &
   & "ahBi7tha","eeS2vahX","moh7maeM","Shiejah4","eCa9bez2","joQu4zie","Aehee5th","Chie5aep", &
   & "Uu0ge6au","aexa1aeR","ohSei2oo","AiPh2ait","iu3yeiFu","too3RahL","Aaqu6li4","aiYier5S", &
   & "chefah9A","Cahcuju4","Eef1Zie8","ahKah9ei","Iene8ahn","ooL3zoo0","aiMei7ga","aeB0ahva", &
   & "aMeeV4ei","Iu2vohzu","equ9BaH3","wuet7Mee","eiVu2jee","Loop7Eil","as8heeNu","chotoj1A", &
   & "aj7eQuah","eiN4weJi","nooWe2pa","Wughe4de","riesahX0","ahZ3she8","Ioviey3w","xiM5eiqu", &
   & "Uaz2Hael","sheiG4gu","nah2Heiy","hei2Wu4e","Aizae2ca","mohs0Wi6","Thir3tah","eiquuQu7", &
   & "Ahjai9ie","ohCe3Uo4","Ohgh5hoh","Izie3cha","ooPh9rae","leiSeiG1","aeDahk0e","Shai8yie", &
   & "Loo1cies","ieXoina6","toopo4Ie","yahc7Sho","geDu9bo6","Ahz0ohdi","waN3ohGu","Pae5kae0", &
   & "Ieno8ain","quai6aeK","eeQu1ieL","iexeeB6v","DohL1eiD","ieGhoh0c","eJ4guqua","Kiek2ieh", &
   & "Othur4ai","Yah3keel","Vun0eshi","ShohB6ze","Zo4aiche","ohl4theT","ahM9quee","eeYo1oqu", &
   & "Cachai5o","Eevoo9ee","oog2shuG","BahJ2shu","gei2iaNg","Niuj4osh","siLu0hae","oht1Tuc6", &
   & "eeFaej8p","ooh4Aok3","duuho8Ah","deiNgeg6","aiZa0Ohs","giLahd0o","dueCh2zo","Yah6aen8", &
   & "eeBaib8e","queXiex3","ooMie6xu","booKoh3A","Yei5ub8w","sahqu1iW","seeK2ier","ug8Ag3te", &
   & "Kiphahh4","Zae8oogh","eedoh8Ba","kah6LieR","bieLie2a","he3si8Po","eer1Ieme","xahy8Rei", &
   & "iech7aWo","uph5Ku2i","ib7ees2D","ohThu3EP","oopeeNg0","EeveeJ2i","ooV0waip","Ohqua7qu", &
   & "eex4yuJ4","ier8Ungu","yae6Uthe","eeXai8ni","Aihook9z","sasoh8Se","Aepi7Ubi","fa0Mabae", &
   & "ieTh7zi8","Iuyies6y","Ahdai0Oo","pohTii2b","Olu9xo9A","shoh1Cae","chas8Aow","ohQuae4x", &
   & "Shais6Ch","SooH8boo","iV5odee5","fo9Laem1","fuDahng6","shei2Umu","daiti9Ae","ACi2ka8c" /)

    CALL SRAND(0)
    do i=1,N
      ivals(i) = fu_irand()   ! irand is not a standard function of FORTRAN. wrapper is in md
      indices(i) = i
    enddo
    call argMergeSort_int(ivals,indices)
    call msg("Array", ivals)
    call msg("Sorted indices", indices)
    call msg("Sorted:")
    do i=1,N
      call msg("",indices(i), ivals(indices(i)))
    enddo

    !!! Strings
    do i=1,N
      indices(i) = i
    enddo
    call argMergeSort_char(strings,indices)
    call msg("Sorted:")
    do i=1,N
      call msg(strings(indices(i)),indices(i))
    enddo




  end subroutine test_sort



    
    
END MODULE toolbox

