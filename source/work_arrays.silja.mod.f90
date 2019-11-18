MODULE work_arrays

  ! Description:
  ! Here are tools for returning big enough work arrays for various uses.
  ! These can be used for reading big datas, when sizes are not knwon yet,
  ! or when calculating new variables or otherwise filling big arrays.
  !
  ! The sizes of returned array can be set (only once, though), and if
  ! not, the default values given in module max_sizes_and_limits are used.
  !
  ! The arrays returned by this module must be used only for work data,
  ! not storage. All work arrays can be freed
  ! by any subroutine at any time.
  !
  ! The main idea is to save memory (using same arrays in various
  ! places) while not having to allocate/deallocate all the time.
  ! 
  ! Language: ANSI standard Fortran 90
  ! 
  ! Original code: Mika Salonoja, FMI
  ! Author: Mikhail Sofiev, FMI
  ! 
  ! Modules used:
  USE globals

  IMPLICIT NONE

!!!!#define TRACE_ARRAYS  
  
  ! The public functions and subroutines available in this module:
  ! INITIALIZATION:
  PUBLIC initialize_work_arrays

  ! RETURNIG WORK ARRAYS:
  PUBLIC fu_work_array
  PUBLIC fu_work_int_array
  PUBLIC fu_work_array_2d
  PUBLIC fu_work_int_array_2d
  public fu_work_string
  public get_work_arrays_set

  public enlarge_array

  ! freeing or killing (deallocating) work arrays
  PUBLIC free_work_array

  ! PRIVATE functions
  private fu_work_array_worksize
  private fu_work_int_array_worksize
  private fu_big_work_array
  private fu_work_array_2d_worksize
  private fu_work_array_2d_givensize
!  private fu_small_work_array   !
  private fu_big_work_int_array
  private get_work_arrays_set_real
  private get_work_arrays_set_real_2d
  private get_work_arrays_set_real_huge
  private get_work_arrays_set_int
  private get_work_arrays_set_int_huge


!  private enlarge_array_of_strings
  private enlarge_array_of_ints
  private enlarge_array_of_reals
  private enlarge_array_of_stringsPtrs
  private resize_work_arrays_set

  private free_work_array_1d
  private free_work_array_2d
  private free_work_int_array_1d
  private free_work_int_array_2d
  private free_work_string
  private free_work_arrays_set_real
  private free_work_arrays_set_real_2d
  private free_work_arrays_set_int


  interface fu_work_array
    module procedure fu_big_work_array
    module procedure fu_work_array_worksize
  end interface

  interface fu_work_int_array
    module procedure fu_work_int_array_worksize
    module procedure fu_big_work_int_array
  end interface

  interface fu_work_array_2d
    module procedure fu_work_array_2d_worksize
    module procedure fu_work_array_2d_givensize
  end interface

  interface get_work_arrays_set
    module procedure get_work_arrays_set_real_huge
    module procedure get_work_arrays_set_real
    module procedure get_work_arrays_set_real_2d
    module procedure get_work_arrays_set_int_huge
    module procedure get_work_arrays_set_int
  end interface

  INTERFACE free_work_array
    MODULE PROCEDURE free_work_array_1d, free_work_array_2d
    module procedure free_work_int_array_1d, free_work_int_array_2d
    module procedure free_work_string
    module procedure free_work_arrays_set_real
    module procedure free_work_arrays_set_int
    module procedure free_work_arrays_set_real_2d
  END INTERFACE
  
  interface enlarge_array
!    module procedure enlarge_array_of_strings
    module procedure enlarge_array_of_ints
    module procedure enlarge_array_of_reals
    module procedure enlarge_array_of_stringsPtrs
    module procedure resize_work_arrays_set
  end interface

  TYPE silja_rp_1d
    REAL, DIMENSION(:), POINTER :: pp
  END TYPE silja_rp_1d
  
  TYPE silja_ip_1d
    INTEGER, DIMENSION(:), POINTER :: ip
  END TYPE silja_ip_1d

  TYPE silja_rp_2d
    REAL, DIMENSION(:,:), POINTER :: pp
  END TYPE silja_rp_2d
  
  TYPE silja_ip_2d
    INTEGER, DIMENSION(:,:), POINTER :: ip
  END TYPE silja_ip_2d
  
  type silam_sp
    character(len=worksize_string), pointer :: sp
  end type silam_sp


  INTEGER, PRIVATE, SAVE :: ws_1d = worksize
  INTEGER, PRIVATE, SAVE :: ws_int_1d = worksize
  INTEGER, PRIVATE, SAVE :: ws_2d_x = worksize_2dx
  INTEGER, PRIVATE, SAVE :: ws_2d_y = worksize_2dy
  LOGICAL, PRIVATE, SAVE :: size_1d_set = .false.
  LOGICAL, PRIVATE, SAVE :: size_int_1d_set = .false.
  LOGICAL, PRIVATE, SAVE :: size_2d_set = .false.

  INTEGER, PARAMETER, PRIVATE :: empty = 0 ! memory not allocated
  INTEGER, PARAMETER, PRIVATE :: in_use = 1 ! allocated and occupied
  INTEGER, PARAMETER, PRIVATE :: free = 2 ! allocated but free
  integer, parameter, private :: huge_size = 3 ! allocated and huge-size

  INTEGER, DIMENSION(max_work_arrays), PRIVATE, SAVE :: status_1d = empty
  INTEGER, DIMENSION(max_work_arrays), PRIVATE, SAVE :: status_int_1d = empty
  INTEGER, DIMENSION(max_work_arrays), PRIVATE, SAVE :: status_string = empty

  TYPE(silja_rp_1d), DIMENSION(max_work_arrays), PRIVATE, SAVE :: works_1d
  TYPE(silja_ip_1d), DIMENSION(max_work_arrays), PRIVATE, SAVE :: works_int_1d
  TYPE(silam_sp), DIMENSION(max_work_arrays), PRIVATE, SAVE :: works_string

  
  
  
CONTAINS

  ! ***************************************************************


  SUBROUTINE initialize_work_arrays()
    !
    ! Nullifies work array pointer. Does not allocate memory
    ! 
    IMPLICIT NONE

    ! Local declarations:
    INTEGER :: i

    status_1d = empty
    status_int_1d = empty

    DO i = 1, max_work_arrays
      NULLIFY(works_1d(i)%pp)
      NULLIFY(works_int_1d(i)%ip)
    END DO
    
  END SUBROUTINE initialize_work_arrays


  ! ***************************************************************


  FUNCTION fu_work_array_worksize()
    ! 
    ! Returns a pointer to a suitable free real work array.
    !
    IMPLICIT NONE
    !
    ! The return value of this function:
    REAL, DIMENSION(:), POINTER :: fu_work_array_worksize

    ! Local declarations:
    INTEGER :: i, ithread, nthreads, status
!!    character (len=clen) :: strTmp
    
    !$ if (omp_in_parallel()) then
    !$    ithread=omp_get_thread_num()
    !$    nthreads=omp_get_num_threads()
    !$ else
        ithread=0
        nthreads=1
    !$ endif

    DO i = 1+ithread, max_work_arrays, nthreads

      SELECT CASE (status_1d(i))

        CASE (empty) 

        ALLOCATE(works_1d(i)%pp(ws_1d), stat = status)
        IF (status /= 0) THEN
          call msg('status:', status)
          CALL set_error('error while allocating','fu_work_array_worksize')
          exit
        END IF

        fu_work_array_worksize => works_1d(i)%pp
        status_1d(i) = in_use
      !        write (strTmp, '(A,I4,A)') 'Thread ', ithread, ' returning REAL work allocated'
        call msg("returning REAL work allocated (thread, WA No, size)", (/ithread, i, ws_1d/))
        exit

        CASE (free)

        fu_work_array_worksize => works_1d(i)%pp
        status_1d(i) = in_use
        exit

      END SELECT

    END DO


    
    if (i > max_work_arrays) then
       call msg('status_1d',status_1d)
       CALL set_error('all work arrays in use','fu_work_array_worksize')
       return
    endif
#ifdef TRACE_ARRAYS           
   call msg('Real warray:',ithread,i)
#endif

  END FUNCTION fu_work_array_worksize


  ! ***************************************************************


  FUNCTION fu_work_int_array_worksize()
    ! 
    ! Returns a pointer to a suitable free integer work array.
    !
    IMPLICIT NONE
    !
    ! The return value of this function:
    integer, DIMENSION(:), POINTER :: fu_work_int_array_worksize

    ! Local declarations:
    INTEGER :: i, ithread, nthreads, status
    
    !$ if (omp_in_parallel()) then
    !$    ithread=omp_get_thread_num()
    !$    nthreads=omp_get_num_threads()
    !$ else
        ithread=0
        nthreads=1
    !$ endif

    DO i = 1+ithread, max_work_arrays, nthreads

      SELECT CASE (status_int_1d(i))

        CASE (empty) 

          ALLOCATE(works_int_1d(i)%ip(ws_1d), stat = status)
          IF (status /= 0) THEN
            call msg('status:', status)
            CALL set_error('error while allocating','fu_work_int_array')
            exit
          END IF

          fu_work_int_array_worksize => works_int_1d(i)%ip
          status_int_1d(i) = in_use
          !!  call msg('Thread '//trim(fu_str(ithread))//' returning INT work allocated', i, ws_1d)
          call msg("returning INT work allocated (thread, WA No, size)", (/ithread, i, ws_1d/))
          exit

        CASE (free)

          fu_work_int_array_worksize => works_int_1d(i)%ip
          status_int_1d(i) = in_use
#ifdef TRACE_ARRAYS
call msg('Int warray:', ithread, i)
#endif
          exit

      END SELECT

    END DO

    if(i > max_work_arrays)then
      call msg('status_int_1d', status_int_1d)
      CALL set_error('all work arrays in use','fu_work_int_array')
    endif

  END FUNCTION fu_work_int_array_worksize


  ! ***************************************************************

  FUNCTION fu_work_array_2d_worksize()
    ! 
    ! Returns a pointer to a suitable free real 2D work array.
    !
    IMPLICIT NONE
    !
    ! The return value of this function:
    REAL, DIMENSION(:,:), POINTER :: fu_work_array_2d_worksize
    REAL, DIMENSION(:), POINTER :: wa1d

    wa1d => fu_work_array(worksize_2dx*worksize_2dy)

    fu_work_array_2d_worksize(1:worksize_2dx,1:worksize_2dy) => wa1d(1:worksize_2dx*worksize_2dy)

  END FUNCTION fu_work_array_2d_worksize


  ! ***************************************************************

  FUNCTION fu_work_array_2d_givensize(n1,n2)
    ! 
    ! Returns a pointer to a suitable free real 2D work array.
    !
    IMPLICIT NONE
    !
    integer, intent(in) :: n1,n2
    !
    ! The return value of this function:
    REAL, DIMENSION(:,:), POINTER :: fu_work_array_2d_givensize
    REAL, DIMENSION(:), POINTER :: wa1d

    wa1d => fu_work_array(n1*n2)

    fu_work_array_2d_givensize(1:n1,1:n2) => wa1d(1:n1*n2)

  END FUNCTION fu_work_array_2d_givensize


  ! ***************************************************************

  FUNCTION fu_work_int_array_2d()
    ! 
    ! Returns a pointer to a suitable free real 2D work array.
    !
    IMPLICIT NONE
    !
    ! The return value of this function:
    INTEGER, DIMENSION(:,:), POINTER :: fu_work_int_array_2d
    INTEGER, DIMENSION(:), POINTER :: wa1d

    
    wa1d => fu_work_int_array(worksize_2dx*worksize_2dy)

    fu_work_int_array_2d(1:worksize_2dx,1:worksize_2dy) => wa1d(1:worksize_2dx*worksize_2dy)

  END FUNCTION fu_work_int_array_2d


  !****************************************************************

  FUNCTION fu_work_string()
    ! 
    ! Returns a pointer to a suitable free string.
    !
    IMPLICIT NONE
    !
    ! The return value of this function:
    character(len=worksize_string), POINTER :: fu_work_string

    ! Local declarations:
    INTEGER :: i, status

!    call msg('Work strings in use: ',count(status_string(:) == in_use))
!    call msg('Work strings empty: ',count(status_str(:) == empty))
!    call msg('Work strings free: ', count(status_str(:) == free))

    !$OMP CRITICAL(work_arrays_get)

    DO i = 1, max_work_arrays

      SELECT CASE (status_string(i))

        CASE (empty) 
          call msg('allocating string work', i)
          ALLOCATE(works_string(i)%sp, stat = status)
          IF(fu_fails(status==0,'error while allocating string','fu_work_string'))exit

          works_string(i)%sp = ''
          fu_work_string => works_string(i)%sp
          status_string(i) = in_use

          exit

        CASE (free)

          works_string(i)%sp = ''
          fu_work_string => works_string(i)%sp
          status_string(i) = in_use
#ifdef TRACE_ARRAYS
call msg('Wstring:',i)
#endif
          exit

      END SELECT

    END DO

    !$OMP END CRITICAL(work_arrays_get)

    if(i > max_work_arrays)then
      call msg('All work strings are in use', status_string)
      CALL set_error('all work strings in use','fu_work_string')
    endif

  END FUNCTION fu_work_string


  !*****************************************************************

  subroutine get_work_arrays_set_real_huge(nArrays, nElements, pArrays)
    !
    ! Returns a bunch of nArrays work arrays
    !
    implicit none
    
    ! Imported parameters
    integer, intent(in) :: nArrays, nElements
    type(silja_rp_1d), dimension(:), pointer :: pArrays

    ! Local variable
    integer :: iTmp

    allocate(pArrays(nArrays), stat = iTmp)
    if(iTmp == 0)then
      if(nElements > ws_1d)then
        do iTmp = 1, nArrays
          pArrays(iTmp)%pp => fu_work_array(nElements)
          if(error)then
            call set_error('Failed one of huge arrays of work array set','get_work_arrays_set_real_huge')
            exit
          endif
        enddo
      else
        do iTmp = 1, nArrays
          pArrays(iTmp)%pp => fu_work_array()
          if(error)then
            call set_error('Failed one of arrays of work array set','get_work_arrays_set_real_huge')
            exit
          endif
        enddo
      endif
    else
      call set_error('Failed allocation of array set','get_work_arrays_set_real_huge')
    endif

  end subroutine get_work_arrays_set_real_huge


  !*****************************************************************

  subroutine get_work_arrays_set_real(nArrays, pArrays)
    !
    ! Returns a bunch of nArrays work arrays
    !
    implicit none
    
    ! Imported parameters
    integer, intent(in) :: nArrays
    type(silja_rp_1d), dimension(:), pointer :: pArrays

    ! Local variable
    integer :: iTmp

    allocate(pArrays(nArrays), stat = iTmp)
    if(iTmp == 0)then
      do iTmp = 1, nArrays
        pArrays(iTmp)%pp => fu_work_array()
        if(error)then
          call set_error('Failed one of arrays of work array set','get_work_arrays_set_real')
          exit
         endif
      enddo
    else
      call set_error('Failed allocation of array set','get_work_arrays_set_real')
    endif

  end subroutine get_work_arrays_set_real
  

  !*****************************************************************

  subroutine get_work_arrays_set_real_2d(nArrays, pArrays)
    !
    ! Returns a bunch of nArrays work arrays
    !
    implicit none
    
    ! Imported parameters
    integer, intent(in) :: nArrays
    type(silja_rp_2d), dimension(:), pointer :: pArrays

    ! Local variable
    integer :: iTmp

    allocate(pArrays(nArrays), stat = iTmp)
    if(iTmp == 0)then
      do iTmp = 1, nArrays
        pArrays(iTmp)%pp => fu_work_array_2d()
        if(error)then
          call set_error('Failed one of arrays of work array set','get_work_arrays_set_real_2d')
          exit
         endif
      enddo
    else
      call set_error('Failed allocation of array set','get_work_arrays_set_real_2d')
    endif

  end subroutine get_work_arrays_set_real_2d
  

  !*****************************************************************

  subroutine get_work_arrays_set_int_huge(nArrays, nElements, pArrays)
    !
    ! Returns a bunch of nArrays work arrays
    !
    implicit none
    
    ! Imported parameters
    integer, intent(in) :: nArrays, nElements
    type(silja_ip_1d), dimension(:), pointer :: pArrays

    ! Local variable
    integer :: iTmp

    allocate(pArrays(nArrays), stat = iTmp)
    if(iTmp == 0)then
      if(nElements > ws_int_1d)then
        do iTmp = 1, nArrays
          pArrays(iTmp)%ip => fu_work_int_array(nElements)
          if(error)then
            call set_error('Failed one of huge arrays of work array set','get_work_arrays_set_int_huge')
            exit
          endif
        enddo
      else
        do iTmp = 1, nArrays
          pArrays(iTmp)%ip => fu_work_int_array()
          if(error)then
            call set_error('Failed one of arrays of work array set','get_work_arrays_set_int_huge')
            exit
          endif
        enddo
      endif
    else
      call set_error('Failed allocation of array set','get_work_arrays_set_int_huge')
    endif

  end subroutine get_work_arrays_set_int_huge

  
  !*****************************************************************

  subroutine get_work_arrays_set_int(nArrays, pArrays)
    !
    ! Returns a bunch of nArrays work arrays
    !
    implicit none
    
    ! Imported parameters
    integer, intent(in) :: nArrays
    type(silja_ip_1d), dimension(:), pointer :: pArrays

    ! Local variable
    integer :: iTmp

    allocate(pArrays(nArrays), stat = iTmp)
    if(iTmp == 0)then
      do iTmp = 1, nArrays
        pArrays(iTmp)%ip => fu_work_int_array()
        if(error)then
          call set_error('Failed one of arrays of work array set','get_work_arrays_set_int')
          exit
         endif
      enddo
    else
      call set_error('Failed allocation of array set','get_work_arrays_set_int')
    endif

  end subroutine get_work_arrays_set_int
  

  !*********************************************************************************

  SUBROUTINE resize_work_arrays_set(array_set, arraysize_1, arraysize_2)
    !
    ! Resizes the existing work array set of 1D real arrays to have the given size. 
    !
    ! Method:
    ! 1. A new work array set is created with proper sizing.
    ! 2. Information is copied to the new structure
    ! 3. Old work array set is released
    ! 
    IMPLICIT NONE
    !
    ! Imported parameters
    INTEGER, INTENT(in) :: arraysize_1, arraysize_2
    type(silja_rp_1d), dimension(:), pointer :: array_set
    
    ! Local declarations:
    INTEGER :: iAr
    type(silja_rp_1d), DIMENSION(:), POINTER :: array_set_tmp
    logical :: ifNew
    !
    ! Exists?
    !
    ifNew = .not. associated(array_set)

    IF (ifNew) THEN      ! Allocate memory for a new, untouched array pointer.

      call get_work_arrays_set(arraysize_1, arraysize_2, array_set)

    ELSE   ! something already there. Have to copy. Note that we cannot reduce the dimensions

      call get_work_arrays_set(max(arraysize_1, size(array_set)), arraysize_2, array_set_tmp)
      
      do iAr = 1, size(array_set)
        array_set_tmp(iAr)%pp(1:min(size(array_set(iAr)%pp),arraysize_2)) = &
                            & array_set(iAr)%pp(1:min(size(array_set(iAr)%pp),arraysize_2))
      enddo
      
      call free_work_arrays_set_real(array_set)

      array_set => array_set_tmp
      
    END IF ! new_or_old
        
  END SUBROUTINE resize_work_arrays_set


  ! ****************************************************************

  SUBROUTINE free_work_array_1d(work_1d)
    ! 
    ! Frees the given array to be used elsewhere in SILAM
    ! Does not deallocate
    ! memory. Call this after finishing a task.
    !
    ! Method:
    ! Find the work array pointer that given pointer points to,
    ! and set the according index = free.
    !
    IMPLICIT NONE

    ! Imported parameters:
    REAL, DIMENSION(:), POINTER :: work_1d

    ! Local declarations:
    INTEGER :: i
    logical :: found

    found = .false.

    DO i = 1, SIZE(works_1d)
      IF (ASSOCIATED(work_1d, works_1d(i)%pp)) THEN
        found = .true.
        if(status_1d(i) == huge_size)then  ! huge arrays are to be deallocated
#ifdef TRACE_ARRAYS
call msg('Deallocate huge real warray:',i)
#endif
          deallocate(works_1d(i)%pp)
          status_1d(i) = empty
          works_1d(i)%pp => null()
          work_1d => null()
        elseif (status_1d(i) == in_use) then
#ifdef TRACE_ARRAYS
call msg('Free real warray:',i)
#endif
          status_1d(i) = free          ! usual ones are just freed
        else
          call set_error("free_work_array_1d for already not busy work array", &
                & "free_work_array_1d")

        endif
        EXIT
      END IF
    END DO
    if (.not. found) then
      call set_error("Failed to find work array to free", "free_work_array_1d")
    endif

  END SUBROUTINE free_work_array_1d


  ! ***************************************************************

  SUBROUTINE free_work_array_2d(work_2d)

    ! Method:
    ! Find the 1d work array pointer and frees it.
    !
    IMPLICIT NONE

    ! Imported parameters:
    REAL, DIMENSION(:,:), POINTER :: work_2d
    REAL,  DIMENSION(:), POINTER ::  work_1d
    REAL, POINTER :: p2d

    ! Local declarations:
    INTEGER :: i
    logical :: found

    found = .false.
    p2d => work_2d(1,1)

    DO i = 1, SIZE(works_1d)
      if (status_1d(i) == empty) cycle
      IF (ASSOCIATED(p2d, works_1d(i)%pp(1))) THEN
        found=.true.
        work_1d => works_1d(i)%pp
        call free_work_array(work_1d)
        work_2d => null()
        EXIT
      END IF
    END DO
    if (.not. found) then
      call set_error("Failed to find 1d work array to free", "free_work_array_2d")
    endif

  END SUBROUTINE free_work_array_2d


  ! ***************************************************************

  SUBROUTINE free_work_int_array_1d(work_int_1d)
    ! 
    ! Frees the given array to be used elsewhere in Silja.
    ! Does not deallocate
    ! memory. Call this after finishing a task.
    !
    ! Method:
    ! Find the work array pointer that given pointer points to,
    ! and set the according index = free.
    !
    IMPLICIT NONE

    ! Imported parameters:
    integer, DIMENSION(:), POINTER :: work_int_1d

    ! Local declarations:
    INTEGER :: i
    logical :: ifFound

    ifFound = .false.

    DO i = 1, SIZE(works_int_1d)
      IF (ASSOCIATED(work_int_1d, works_int_1d(i)%ip)) THEN
        ifFound = .true.
        if(status_int_1d(i) == huge_size)then  ! huge arrays are to be deallocated

#ifdef TRACE_ARRAYS           
    call msg('Deallocate huge integer warray:',i)
#endif
          deallocate(works_int_1d(i)%ip)
          status_int_1d(i) = empty
          works_int_1d(i)%ip => null()
        else
#ifdef TRACE_ARRAYS
    call msg('Free integer warray:',i)
#endif
          status_int_1d(i) = free
        endif
        work_int_1d => null()
        EXIT
      END IF
    END DO
    if (.not. ifFound) then
            call set_error("Failed to find work array to free", "free_work_int_array_1d")
    endif

  END SUBROUTINE free_work_int_array_1d


  ! ***************************************************************


  SUBROUTINE free_work_int_array_2d(work_2d)
    ! Method:
    ! Find the 1d work array pointer and frees it.
    !
    IMPLICIT NONE

    ! Imported parameters:
    INTEGER, DIMENSION(:,:), POINTER :: work_2d
    INTEGER, DIMENSION(:), POINTER ::  work_1d
    INTEGER, POINTER :: p2d

    ! Local declarations:
    INTEGER :: i
    logical :: found

    found = .false.
    p2d => work_2d(1,1)

    DO i = 1, SIZE(works_int_1d)
      IF (ASSOCIATED(p2d, works_int_1d(i)%ip(1))) THEN
        found=.true.
        work_1d => works_int_1d(i)%ip
        call free_work_array(work_1d)
        work_2d => null()
        EXIT
      END IF
    END DO
    if (.not. found) then
      call set_error("Failed to find 1d int work array to free", "free_work_int_array_2d")
    endif

  END SUBROUTINE free_work_int_array_2d


  ! ***************************************************************


  SUBROUTINE free_work_string(work_string)
    ! 
    ! Frees the given array to be used elsewhere in Silja.
    ! Does not deallocate
    ! memory. Call this after finishing a task.
    !
    ! Method:
    ! Find the work array pointer that given pointer points to,
    ! and set the according index = free.
    !
    IMPLICIT NONE

    ! Imported parameters:
    character(len=*), POINTER :: work_string

    ! Local declarations:
    INTEGER :: i

    DO i = 1, SIZE(works_string)
      IF (ASSOCIATED(work_string, works_string(i)%sp)) THEN
        status_string(i) = free
        work_string => null()
#ifdef TRACE_ARRAYS
    call msg('Free wstring:',i)
#endif
        EXIT
      END IF
    END DO

  END SUBROUTINE free_work_string


  !*****************************************************************

  subroutine free_work_arrays_set_int(pArrays)
    !
    ! Frees the work arrays and deallocates the structure
    !
    implicit none
    
    ! Imported parameters
    type(silja_ip_1d), dimension(:), pointer :: pArrays

    ! Local variable
    integer :: iTmp

    if(.not. associated(pArrays))return

    do iTmp = 1, size(pArrays)
      call free_work_array(pArrays(iTmp)%ip)
    enddo

    deallocate(pArrays)

  end subroutine free_work_arrays_set_int


  !*****************************************************************

  subroutine free_work_arrays_set_real(pArrays)
    !
    ! Frees the work arrays and deallocates the structure
    !
    implicit none
    
    ! Imported parameters
    type(silja_rp_1d), dimension(:), pointer :: pArrays

    ! Local variable
    integer :: iTmp

    if(.not. associated(pArrays))return

    do iTmp = 1, size(pArrays)
      call free_work_array(pArrays(iTmp)%pp)
    enddo

    deallocate(pArrays)
    pArrays => null()

  end subroutine free_work_arrays_set_real


  !*****************************************************************

  subroutine free_work_arrays_set_real_2d(pArrays)
    !
    ! Frees the work arrays and deallocates the structure
    !
    implicit none
    
    ! Imported parameters
    type(silja_rp_2d), dimension(:), pointer :: pArrays

    ! Local variable
    integer :: iTmp

    if(.not. associated(pArrays))return

    do iTmp = 1, size(pArrays)
      call free_work_array(pArrays(iTmp)%pp)
    enddo

    deallocate(pArrays)

  end subroutine free_work_arrays_set_real_2d


  ! ***************************************************************



   
  !******************************************************************************

  subroutine enlarge_array_of_stringsPtrs(ptrs, sizeNew)
    !
    ! Increase the size of the derived-type array of strings.
    !
    implicit none 

    ! Imported parameters
    type(silam_sp), dimension(:), pointer :: ptrs
    integer, intent(in), optional :: sizeNew

    ! Local variables
    integer :: iStatus, i, iSize, iSizeNew
    type(silam_sp), dimension(:), pointer :: vTmp

    if(associated(ptrs))then

      iSize = size(ptrs)
      if(present(sizeNew))then
        iSizeNew = sizeNew
      else
        iSizeNew = 10+iSize
      endif
      if(iSize >= iSizeNew) return

      if(associated(vTmp))then  ! Create temporary array of strings
        if(size(vTmp) <= iSize)then  ! enlarge the size of temporary string pointers
          deallocate(vTmp)
          allocate(vTmp(iSize), stat = iStatus)
        endif
      else
        allocate(vTmp(iSize), stat = iStatus)
      endif
      if(iStatus /= 0)then
        call set_error('Failed to allocate memory for values.0','enlarge_array_of_stringsPtrs')
        return
      endif
      do i=1,iSize
        vTmp(i)%sp => ptrs(i)%sp
      end do

      deallocate(ptrs)
      allocate(ptrs(iSizeNew),stat=iStatus)
      if(iStatus /= 0)then
        call set_error('Failed to allocate memory for values.2','enlarge_array_of_stringsPtrs')
        return
      endif
      do i=1,iSize
        ptrs(i)%sp => vTmp(i)%sp
      end do
      do i=iSize+1, iSizeNew
        allocate(ptrs(i)%sp, stat=iStatus)
        if(iStatus /= 0)then
          call set_error('Failed to allocate memory for values.2.5','enlarge_array_of_stringsPtrs')
          return
        endif
        ptrs(i)%sp = ''
      end do

    else
      if(present(sizeNew))then
        iSizeNew = sizeNew
      else
        iSizeNew = 10
      endif
      allocate(ptrs(iSizeNew),stat=iStatus)
      if(iStatus /= 0)then
        call set_error('Failed to allocate memory for values.3','enlarge_array_of_stringsPtrs')
        return
      endif
      do i=1,iSizeNew
        allocate(ptrs(i)%sp, stat=iStatus)
        if(iStatus /= 0)then
          call set_error('Failed to allocate memory for values.3.5','enlarge_array_of_stringsPtrs')
          return
        endif
        ptrs(i)%sp = ''
      end do
    endif
  end subroutine enlarge_array_of_stringsPtrs



  !******************************************************************************

  subroutine enlarge_array_of_ints(ptrs, sizeNew)
    !
    ! A stupid task - increase the size of the derived-type array
    ! of integers. 
    !
    implicit none 

    ! Imported parameters
    integer, dimension(:), pointer :: ptrs
    integer, intent(in), optional :: sizeNew

    ! Local variables
    integer :: iStatus, i, iSize, iSizeNew
    integer, dimension(:), pointer :: vTmp

    if(associated(ptrs))then
  
      iSize = size(ptrs)
      if(present(sizeNew))then
        iSizeNew = sizeNew
      else
        iSizeNew = 10+iSize
      endif
      if(iSize >= iSizeNew) return

      vTmp => fu_work_int_array()
      if(error)then
        call set_error('Failed to allocate memory for values.1','enlarge_array_of_ints')
        return
      endif
      do i=1,iSize
        vTmp(i) = ptrs(i)
      end do
      deallocate(ptrs)
      allocate(ptrs(iSizeNew),stat=iStatus)
      if(iStatus /= 0)then
        call set_error('Failed to allocate memory for values.2','enlarge_array_of_ints')
        return
      endif
      do i=1,iSize
        ptrs(i) = vTmp(i)
      end do
      ptrs(iSize+1:)=int_missing
      call free_work_array(vTmp)

    else
      if(present(sizeNew))then
        iSizeNew = sizeNew
      else
        iSizeNew = 10
      endif
      allocate(ptrs(iSizeNew),stat=iStatus)
      if(iStatus /= 0)then
        call set_error('Failed to allocate memory for values.3','enlarge_array_of_ints')
        return
      endif
    endif
  end subroutine enlarge_array_of_ints


  !******************************************************************************

  subroutine enlarge_array_of_reals(ptrs, sizeNew)
    !
    ! A stupid task - increase the size of the derived-type array
    ! of the emission cells. We store not the map but the list of cells,
    ! which is better if less than 30% of the map is covered by emission
    !
    implicit none 

    ! Imported parameters
    real, dimension(:), pointer :: ptrs
    integer, intent(in), optional :: sizeNew

    ! Local variables
    integer :: iStatus, i, iSize, iSizeNew
    real, dimension(:), pointer :: vTmp

    if(associated(ptrs))then
  
      iSize = size(ptrs)
      if(present(sizeNew))then
        iSizeNew = sizeNew
      else
        iSizeNew = 10+iSize
      endif
      if(iSize >= iSizeNew) return

      vTmp => fu_work_array()
      if(error)then
        call set_error('Failed to allocate memory for values.1','enlarge_array_of_reals')
        return
      endif
      do i=1,iSize
        vTmp(i) = ptrs(i)
      end do
      deallocate(ptrs)
      allocate(ptrs(iSizeNew),stat=iStatus)
      if(iStatus /= 0)then
        call set_error('Failed to allocate memory for values.2','enlarge_array_of_reals')
        return
      endif
      do i=1,iSize
        ptrs(i) = vTmp(i)
      end do
      ptrs(iSize:) = real_missing
      call free_work_array(vTmp)

    else
      if(present(sizeNew))then
        iSizeNew = sizeNew
      else
        iSizeNew = 10
      endif
      allocate(ptrs(iSizeNew),stat=iStatus)
      if(iStatus /= 0)then
        call set_error('Failed to allocate memory for values.3','enlarge_array_of_reals')
        return
      endif
    endif
  end subroutine enlarge_array_of_reals


  !***********************************************************************

  FUNCTION fu_big_work_array(nPoints)
    ! 
    ! Returns a pointer to a suitable free real work array.
    ! Verifies the size of the work_array. Should it be insufficient, it allocates the
    ! needed space. Such data are allowed to be huge.
    !
    IMPLICIT NONE
    !
    ! The return value of this function:
    REAL, DIMENSION(:), POINTER :: fu_big_work_array

    ! Imported parameter
    integer, intent(in) :: nPoints

    ! Local declarations:
    INTEGER :: i, ithread, nthreads, status
    

    if(nPoints > ws_1d)then  ! Need to get large array?
      !
      ! Large array needed. Find the empty pointer place, allocate it to a
      ! proper size and return as a pointer
      
      !$ if (omp_in_parallel()) then
      !$    ithread=omp_get_thread_num()
      !$    nthreads=omp_get_num_threads()
      !$ else
          ithread=0
          nthreads=1
      !$ endif

    DO i = 1+ithread, max_work_arrays, nthreads

        SELECT CASE (status_1d(i))

          CASE (empty) 

          ALLOCATE(works_1d(i)%pp(nPoints), stat = status)
          IF (status /= 0) THEN
            call msg('status:', status)
            CALL set_error('error while allocating big array','fu_big_work_array')
            exit
          END IF

          fu_big_work_array => works_1d(i)%pp
          status_1d(i) = huge_size
          call msg('Returning HUGE REAL work allocated (thread,No,size)', (/ithread,i, nPoints/))
          exit

        END SELECT

      END DO
      if (i > max_work_arrays) then
         call msg('status_1d',status_1d)
         CALL set_error('all work arrays in use','fu_big_work_array')
      endif

    else

      fu_big_work_array => fu_work_array()

    endif   ! if worksize too small => allocation needed

  end function fu_big_work_array


  !***********************************************************************

  FUNCTION fu_big_work_int_array(nPoints)
    ! 
    ! Returns a pointer to a suitable free real work array.
    ! Verifies the size of the work_array. Should it be insufficient, it allocates the
    ! needed space. Such data are allowed to be huge.
    !
    IMPLICIT NONE
    !
    ! The return value of this function:
    integer, DIMENSION(:), POINTER :: fu_big_work_int_array

    ! Imported parameter
    integer, intent(in) :: nPoints

    ! Local declarations:
    INTEGER :: i, ithread, nthreads, status
    


    if(nPoints > ws_int_1d)then  ! Need to get large array?
      !
      ! Large array needed. Find the empty pointer place, allocate it to a
      ! proper size and return as a pointer
      
      !$ if (omp_in_parallel()) then
      !$    ithread=omp_get_thread_num()
      !$    nthreads=omp_get_num_threads()
      !$ else
          ithread=0
          nthreads=1
      !$ endif
      DO i = 1+ithread, max_work_arrays, nthreads

        IF (status_int_1d(i) == empty) then
          ALLOCATE(works_int_1d(i)%ip(nPoints), stat = status)
          IF (status /= 0) THEN
            call msg('status:', status)
            CALL set_error('error while allocating big array','fu_big_work_array')
            exit
          END IF

          fu_big_work_int_array => works_int_1d(i)%ip
          status_int_1d(i) = huge_size
          call msg('Returning HUGE INT work allocated (thread,No,size)', (/ithread,i, nPoints/))
          exit

          END IF

      END DO

      if(i > max_work_arrays)then
      call msg('status_int_1d',status_int_1d)
      CALL set_error('all work arrays in use','fu_big_work_int_array')
      endif

    else

      fu_big_work_int_array => fu_work_int_array()

    endif   ! if worksize too small => allocation needed

  end function fu_big_work_int_array

END MODULE work_arrays

