MODULE supermarket_of_fields

  ! Description:
  ! This module contains tools for storing and returning field datas
  ! to supermarket data-pool. The fields found are kept in memory
  ! inside this module as a vector of stacks. Each stack is single
  ! time type and single-source , containing all fields from one model
  ! for one observation time.
  !
  ! There are two methods for getting weather data in the supermarket:
  ! fill_supermarket and store_input_to_supermarket.
  !
  ! Original code: Mika Salonoja
  ! Author: Mikhail Sofiev, FMI email mikhail.sofiev@fmi.fi
  !
  ! All units: SI
  ! 
  ! All degrees: real numbers in degrees (degrees and parts per
  ! hundred/parts per thousand , NOT in degrees and minutes
  ! 
  ! Language: ANSI standard Fortran 90
  ! 
  ! Modules used:
  USE optimisation
  USE netcdf_io
  use grads_io
  use grib_api_io
  use shopping_lists
  use ascii_io
!  use grads_templates

  IMPLICIT NONE

  ! The public functions and subroutines available in this module:
  PUBLIC initialize_mini_market
  PUBLIC set_supermarket_empty ! frees all memory
  public nullify_minimarket    ! forces all pointers to null, no checking.
  public set_time_direction_sm
  public fu_get_time_direction_sm
  PUBLIC arrange_supermarket
  public arrange_supermarket_multitime !Same, but only for multitime stacks
  public fill_minimarket_from_namelist
  PUBLIC fill_meteo_market
  PUBLIC store_input_to_supermarket
  public open_input_file
  public close_input_file
  public get_id_list_from_input_file
  public get_input_field
  public write_stack_to_files

  public init_singletime_fields
!  public merge_data_to_supermarket

  public fu_get_field_from_mm_general
  PUBLIC fu_sm_obstime_field
  PUBLIC fu_sm_obstime_3d_field
  PUBLIC fu_sm_obstime_3d_windfield
  PUBLIC fu_sm_simple_field    ! permanent or dispersion field
  public fu_sm_simple_3d_field ! permanent or dispersion field
  public find_field_data_storage_2d
  public find_field_storage_2d

  PUBLIC fu_supermarket_cpu_usage
  PUBLIC supermarket_times
  public report
  PUBLIC supermarket_quantities
  PUBLIC supermarket_variables
  PUBLIC supermarket_3d_quantities
  PUBLIC supermarket_2d_quantities
  PUBLIC fu_field_in_sm
  PUBLIC fu_closest_sm_met_src_time
  public fu_stack
  public fu_name

  public minimarket_initialized
 
  PUBLIC dq_store_2d     

  PUBLIC supermarket_test_fill

  public check_bm_for_times_in_list
  public set_met_srcs_in_sm_single_wdr
  public supermarket_met_srcs
  public check_supermarket_fields_range

  ! The private functions and subroutines not to be used elsewhere:
  private fu_stack_by_indices
  private fu_stack_by_values
  private fu_minimarket_name
  PRIVATE check_sm_for_data_in_list
  PRIVATE find_field_from_hit_list
  PRIVATE put_field_to_hit_list
  PRIVATE put_wind_to_hit_list
  PRIVATE find_wind_from_hit_list
  PRIVATE fu_met_src_storage_index
  private fu_time_storage_index
  PRIVATE update_obstimes
  private check_grid
  private print_supermarket_contents

  interface fu_name
    module procedure fu_minimarket_name
  end interface

  interface fu_stack
    module procedure fu_stack_by_indices
    module procedure fu_stack_by_values
  end interface
  
  interface report
    module procedure print_supermarket_contents
  end interface


  !=================================================================
  !
  ! The main structure of the supermarket is a mini_market_of_fields, which contains
  ! one stackvector and one single stack. A mini-market is then a single largest structure
  ! that can contain the data. Idea is to have several of them whenever needed.
  !
  ! Mini-market contains two sets of stacks: a 2D (met_src,time) vector of stacks 
  ! for multi-time data, and 1D (met_src) vector of stacks for single-time data.
  ! Those include, however, all sorts of times and validities: long-valid, monthly-valid, etc.
  ! Therefore, the single-time stacks NEVER have their own valid_time defined.
  !
  type mini_market_of_stacks
    private
    !
    ! Main data and time structures
    !
    type(silam_stack_ptr), dimension(:,:),pointer :: stacks_multiTime
!    type(silam_stack_ptr), dimension(:), pointer ::  stacks_irregularTime
    type(silam_stack_ptr), dimension(:), pointer ::  stacks_singleTime
    type(silja_time), dimension(:,:), pointer :: obstimes
    !
    ! and some supplementary variables relating to this very mini-market
    !
    character(len=20) :: name
    logical :: stack_multiTime_initialized = .false., &    ! initialization is allowed only once
             & stack_multiTime_exists = .false., &         ! whether time-dependent stacks exist
             & stack_singleTime_initialized = .false., &   ! can be initialised onlyl once
             & stack_singleTime_exists = .false.           ! whether time-INdependent stacks exist
    logical :: minimarket_empty = .true.
    logical :: replace_earliest = .true.
    !
    ! Here are the hit lists: the latest hits among data 
    !
    type(silja_fieldpointer), dimension(30) :: hit_list
    type(silja_windfieldpointer), dimension(20) :: wind_list
    integer :: hit_list_pointer = 0, wind_list_pointer = 0

  end type mini_market_of_stacks
  public mini_market_of_stacks

  !
  ! Types of datasets
  !
  integer, public, parameter :: meteo_dynamic_flag = 51001
  integer, public, parameter :: meteo_single_time_flag  = 51002
  integer, public, parameter :: dispersion_dynamic_flag = 51003
  integer, public, parameter :: dispersion_single_time_flag = 51004
  integer, public, parameter :: multi_time_stack_flag = 51005
  integer, public, parameter :: single_time_stack_flag  = 51006

  !
  ! Public variables for counting fields delivered:
  !
  REAL, PRIVATE :: cpu_usage = 0.
  INTEGER, PUBLIC :: hit_count = 0
  INTEGER, PUBLIC :: hit_miss = 0
  INTEGER, PUBLIC :: fields_delivered = 0
  INTEGER, PUBLIC :: windfields_delivered = 0
  INTEGER, PUBLIC :: fields_3d_delivered = 0
  INTEGER, PUBLIC :: windfields_3d_delivered = 0
  LOGICAL, PUBLIC, SAVE :: supermarket_info = .false.

  ! For direct operation with 2d fields
  ! Similar to one in field_buffer
  ! 
  TYPE Tfld_ptr
    REAL, DIMENSION(:), POINTER :: ptr => null()
    type(silja_field_id), pointer :: idPtr =>null()
  END TYPE Tfld_ptr

CONTAINS

  ! ***************************************************************

  ! ***************************************************************
  !
  !
  !   PUBLIC TOOLS FOR ALTERING THE WHOLE SUPERMARKET
  !
  !    Initialize, organize, make it empty etc.
  !     
  !
  ! ***************************************************************
  ! ***************************************************************

  SUBROUTINE initialize_mini_market(miniMarketPtr, &
                                  & chName, &
                                  & nbr_of_srcs,&      ! first dimension of the stack arrays
                                  & nbr_of_timenodes,& ! if>0, the second dimension in multiTime stack
                                  & nbr_of_fields,&        ! for each timenode, if any
                                  & nbr_of_windfields,&    ! for each timenode, if any
                                  & nbr_of_3d_fields,&     ! for each timenode, if any
                                  & nbr_of_3d_windfields,& ! for each timenode, if any
                                  & replace_earliest_when_full,&
                                  & wdrAr, &
                                  & ifSingleSrc, &
                                  & print_info_to_stdout)
    ! 
    ! Sets the sizes of data-fields kept in memory. Calling this routine is obligatory.
    !
    ! Parameter replace_oldest_when_full defines what to do when
    ! all timenodes are in use. Normally this should be set true
    ! for normal forwards-in-time dispersion modelling, and false
    ! for backward-in-time simulations.
    ! 
    ! In supermarket there is either several single-time-stacks or
    ! just one multi-time-stack, where all data is stored. What is initialized is decided 
    ! from nbr_of_timenodes, which positive value turns the sub to multiTime initialization
    !
    ! If fields_3d_too is set false, then the 3d-fieldsizes have no
    ! meaning.
    !
    IMPLICIT NONE

    ! Imported parameters with intent IN:
    type(mini_market_of_stacks), intent(inout) :: miniMarketPtr
    character(len=*), intent(in) :: chName
    INTEGER, INTENT(in) :: nbr_of_srcs
    INTEGER, INTENT(in) :: nbr_of_timenodes
    INTEGER, INTENT(in) :: nbr_of_fields
    INTEGER, INTENT(in) :: nbr_of_windfields
    INTEGER, INTENT(in) :: nbr_of_3d_fields
    INTEGER, INTENT(in) :: nbr_of_3d_windfields
    LOGICAL, INTENT(in) :: replace_earliest_when_full, ifSingleSrc
    type(wdr_ptr), dimension(:), pointer :: wdrAr
    LOGICAL :: print_info_to_stdout ! if true, then some info is printed
    ! while retrieving data to supermarket

    ! Local declarations:
    INTEGER :: i, j, status !, number_of_timenodes
    CHARACTER(LEN=clen) :: name_of_stack
    type(silja_wdr), pointer :: wdrPtr
    !----------------------------------------
    !
    ! 1. Allocate the stack-vector.
    !
    miniMarketPtr%name = trim(chName)

    if(nbr_of_timenodes > 0)then
      !
      ! Initialize the multiTime stack matrix
      !
      IF (miniMarketPtr%stack_multiTime_initialized) THEN
        call msg('stack_multiTime is already initialized.')
        RETURN
      END IF

      if(error .or. nbr_of_timenodes > max_times)then
        call set_error('Strange number of time nodes:' + fu_str(nbr_of_timenodes), &
                     & 'initialize_mini_market')
        return
      endif

      ALLOCATE(miniMarketPtr%stacks_multiTime(nbr_of_srcs,nbr_of_timenodes),stat=status)
      IF (status /= 0) THEN
        call msg('Allocation status: ', status)
        CALL set_error('Cannot ALLOCATE space for multiTime stack matrix','initialize_mini_market')
        RETURN
      END IF

      ALLOCATE(miniMarketPtr%obstimes(nbr_of_srcs, nbr_of_timenodes), stat=status)
      IF (status /= 0) THEN
        call msg('Allocation status: ', status)
        CALL set_error ( ' Cannot ALLOCATE space for obstimes', 'initialize_mini_market')
        RETURN
      END IF

      miniMarketPtr%obstimes = time_missing
      miniMarketPtr%replace_earliest = replace_earliest_when_full

      !----------------------------------------
      !
      ! 2. Initialize stackvector
      !
!      call msg('Starting multiTime stack vector')
      DO i = 1, nbr_of_srcs
        if(size(wdrAr)==1)then
          wdrPtr => wdrAr(1)%ptr
        elseif(size(wdrAr)==nbr_of_srcs)then
          wdrPtr => wdrAr(i)%ptr
        else
          call msg('Strange nr of met srces',size(wdrAr))
          call set_error('Strange nr of met srces', 'initialize_mini_market')
          return
        endif
!        call report(fu_storage_grid(wdrPtr))
        DO j = 1, nbr_of_timenodes
!          WRITE(unit = name_of_stack, fmt = '(A, I3, A, I3)') trim(miniMarketPtr%name), i, '_', j
          
          name_of_stack = miniMarketPtr%name + fu_str(i,3) + '_' + fu_str(j)

!          call msg(name_of_stack)
          allocate(miniMarketPtr%stacks_multiTime(i, j)%ptr, stat=status)
          IF (status /= 0) THEN
            call msg('Allocation status: ', status)
            CALL set_error('Cannot ALLOCATE stack pointer', 'initialize_mini_market')
            RETURN
          END IF

          call set_defined(miniMarketPtr%stacks_multiTime(i, j)%ptr, silja_false)
          CALL init_stack(nbr_of_fields,&
                        & nbr_of_windfields,&
                        & TRIM(ADJUSTL(name_of_stack)),&
                        & .true., &          ! single time stacks
                        & ifSingleSrc, &          ! single met_src stacks
                        & nbr_of_3d_fields,& ! for each timenode
                        & nbr_of_3d_windfields, &
                        & wdrPtr, &
                        & miniMarketPtr%stacks_multiTime(i, j)%ptr)
          IF (error) RETURN
!call report(fu_storage_grid(fu_wdr(miniMarketPtr%stacks_multiTime(i, j)%ptr)))
        END DO
      END DO

      miniMarketPtr%stack_multiTime_initialized = .true.
      miniMarketPtr%stack_multiTime_exists = .true.
      IF (.not. miniMarketPtr%stack_singleTime_initialized) miniMarketPtr%minimarket_empty  = .true.

    else
      !
      ! Number of time nodes is 0 or negative - initialize the singleTime stack vector
      !
      IF (miniMarketPtr%stack_singleTime_initialized) THEN
        call msg('stack_singleTime is already initialized.')
        RETURN
      END IF

      ALLOCATE(miniMarketPtr%stacks_singleTime(nbr_of_srcs),stat=status)
      IF (status /= 0) THEN
        call msg('Allocation status: ', status)
        CALL set_error('Cannot ALLOCATE space for singleTime stackvector','initialize_mini_market')
        RETURN
      END IF

      DO i = 1, nbr_of_srcs
        WRITE(unit = name_of_stack, fmt = '(A, A, I3)') 'permanent', trim(miniMarketPtr%name), i
        if(size(wdrAr)==1)then
          wdrPtr => wdrAr(1)%ptr
        elseif(size(wdrAr)==nbr_of_srcs)then
          wdrPtr => wdrAr(i)%ptr
        else
          call msg('Strange nr of met srces',size(wdrAr))
          call set_error('Strange nr of met srces', 'initialize_mini_market')
          return
        endif
        
        allocate(miniMarketPtr%stacks_singleTime(i)%ptr, stat = status)
        IF (status /= 0) THEN
          call msg('Allocation status: ', status)
          CALL set_error ( ' Cannot ALLOCATE stack pointer', 'initialize_mini_market')
          RETURN
        END IF
        
        call set_defined(miniMarketPtr%stacks_singleTime(i)%ptr, silja_false)
        CALL init_stack(nbr_of_fields, &       ! number of fields
                      & nbr_of_windfields, &   ! number of windfields
                      & name_of_stack,&        ! name of stack
                      & .false., &             ! single_time
                      & ifSingleSrc, &         ! single_met_src
                      & nbr_of_3d_fields, &       ! number of 3d fields
                      & nbr_of_3d_windfields, &   ! number of 3d windfields
                      & wdrPtr,& 
                      & miniMarketPtr%stacks_singleTime(i)%ptr)  ! stack
        if(error) return
      end do

      miniMarketPtr%stack_singleTime_initialized = .true.
      miniMarketPtr%stack_singleTime_exists = .true.
      IF (.not. miniMarketPtr%stack_multiTime_initialized) miniMarketPtr%minimarket_empty  = .true.

    endif  ! whether multi- or single-time initialization 

    ! Nullify the hit lists:

    do i = 1, size(miniMarketPtr%hit_list)
      nullify(miniMarketPtr%hit_list(i)%fp)
    end do
    do i = 1, size(miniMarketPtr%wind_list)
      nullify(miniMarketPtr%wind_list(i)%fp)
    end do

    supermarket_info = print_info_to_stdout

  END SUBROUTINE initialize_mini_market


  ! ***************************************************************


  SUBROUTINE arrange_supermarket_multitime(miniMarketPtr)

    ! Description:
    ! Forces all arrangements to be done in supermarket.
    ! Updates the valid-time list.
    !
    IMPLICIT NONE
    ! Argument:
    type(mini_market_of_stacks), intent(inout) :: miniMarketPtr


    ! Local declarations:
    INTEGER :: i, j
    TYPE(silja_stack), POINTER :: stack

    IF (miniMarketPtr%stack_multiTime_initialized) then
      DO i = 1, SIZE(miniMarketPtr%stacks_multiTime, 1)
        DO j = 1, SIZE(miniMarketPtr%stacks_multiTime, 2)
            stack => miniMarketPtr%stacks_multiTime(i, j)%ptr
            IF (defined(stack)) CALL arrange_fields_in_stack(stack, .true.)
          IF (error) RETURN
        END DO
      END DO
      CALL update_obstimes(miniMarketPtr)
    end if

  END SUBROUTINE arrange_supermarket_multitime

  ! ***************************************************************


  SUBROUTINE arrange_supermarket(miniMarketPtr, ifLookForPressure)

    ! Description:
    ! Forces all arrangements to be done in supermarket.
    ! Updates the valid-time list.
    !
    IMPLICIT NONE
    ! Argument:
    type(mini_market_of_stacks), intent(inout) :: miniMarketPtr
    logical, optional, intent (in) :: ifLookForPressure ! Look for matching pressure for
                                              ! 3D fields in hybrid coordinates


    ! Local declarations:
    INTEGER :: i, j
    TYPE(silja_stack), POINTER :: stack
    logical MyifLookForPressure

    if ( present(ifLookForPressure)) then 
            MyifLookForPressure = ifLookForPressure
    else
            MyifLookForPressure = .true.
    endif

    IF (miniMarketPtr%stack_multiTime_initialized) then
      DO i = 1, SIZE(miniMarketPtr%stacks_multiTime, 1)
        DO j = 1, SIZE(miniMarketPtr%stacks_multiTime, 2)
            stack => miniMarketPtr%stacks_multiTime(i, j)%ptr
            IF (defined(stack)) CALL arrange_fields_in_stack(stack, MyifLookForPressure)
          IF (error) RETURN
        END DO
      END DO
      CALL update_obstimes(miniMarketPtr)
    end if

    if (miniMarketPtr%stack_singleTime_initialized) then
      do i = 1, size(miniMarketPtr%stacks_singleTime)
        stack => miniMarketPtr%stacks_singleTime(i)%ptr
        if (defined(stack)) then

                call arrange_fields_in_stack(stack, MyifLookForPressure)
                if (error) return
!                call msg("-----------AFTER arrange--------"+ fu_name(stack))
!                call report(fu_sm_simple_field(miniMarketPtr, met_src_missing, &
!                    & ground_pressure_flag, level_missing, single_time_stack_flag))
!                if (error) call unset_error("Tried to report..")
!               !call report(fu_stack(miniMarketPtr,met_src_missing))
!                call msg("-----------AFTER Store--------")
        endif
      end do
    end if

  END SUBROUTINE arrange_supermarket


  ! ***************************************************************


  SUBROUTINE set_supermarket_empty(miniMarketPtr, ifSingleTime, ifMultiTime)

    ! Description:
    ! Frees all memory, all data is lost. After this routine the
    ! supermarket has to be initialized again for memory-allocation.
    ! 
    ! Language: ANSI Fortran 90
    !
    ! Author: Mika Salonoja, FMI
    ! 
    IMPLICIT NONE

    ! Local declarations:
    type(mini_market_of_stacks), intent(inout) :: miniMarketPtr
    INTEGER :: i, j, status
    TYPE(silja_stack), POINTER :: stack
    logical, intent(in) :: ifSingleTime, ifMultiTime

!    PRINT *, 'Setting supermarket empty...'

    !
    ! We start from multi-time supermarket
    !
    if(ifMultiTime)then
      if(miniMarketPtr%stack_multiTime_initialized)then
        DO i = 1, SIZE(miniMarketPtr%stacks_multiTime, 1)
          DO j = 1, SIZE(miniMarketPtr%stacks_multiTime, 2)
            stack => miniMarketPtr%stacks_multiTime(i, j)%ptr
!!!$          PRINT *, 'setting empty:'
!!!$          CALL report(stack)
            CALL set_stack_empty_free_memory(stack)
            
            deallocate(stack, stat=status)
            if (status /= 0) then
              call set_error('Failed to deallocate stack pointer', 'set_supermarket_empty')
              return
            end if
            
            IF (error) RETURN
          END DO
        END DO
        DEALLOCATE(miniMarketPtr%stacks_multiTime)
        DEALLOCATE(miniMarketPtr%obstimes)
        miniMarketPtr%stack_multiTime_initialized = .false.
      endif
    endif
    !
    ! Now - the singleTime one
    !
    if(ifSingleTime)then
      if(miniMarketPtr%stack_singleTime_initialized)then
        DO i = 1, SIZE(miniMarketPtr%stacks_singleTime)
          stack => miniMarketPtr%stacks_singleTime(i)%ptr
          CALL set_stack_empty_free_memory(stack)

          deallocate(stack)
          if (status /= 0) then
            call set_error('Failed to deallocate stack pointer', 'set_supermarket_empty')
            return
          end if
          
          IF (error) RETURN
        END DO
        DEALLOCATE(miniMarketPtr%stacks_singleTime)
        miniMarketPtr%stack_singleTime_initialized = .false.
      endif
    endif
    if(ifSingleTime .and. ifMultiTime)then
      miniMarketPtr%minimarket_empty = .true.
    endif
  END SUBROUTINE set_supermarket_empty


  ! ***************************************************************


  SUBROUTINE nullify_minimarket(miniMarketPtr)
    ! 
    ! Forces all pointers to null.
    ! ATTENTION. in case minimarket is initialised, loss of memory will happen
    ! 
    IMPLICIT NONE

    ! Imported parameters
    type(mini_market_of_stacks), intent(inout) :: miniMarketPtr

    nullify(miniMarketPtr%stacks_multiTime)
    nullify(miniMarketPtr%stacks_singleTime)
    nullify(miniMarketPtr%obstimes)

    miniMarketPtr%minimarket_empty = .true.

  END SUBROUTINE nullify_minimarket


  ! ***************************************************************


  subroutine set_met_srcs_in_sm_single_wdr (miniMarketPtr, ifMultiTime, ifSingleTime, wdr)
    !
    ! Sets the met_srcs to the stackvector from the wdr
    ! Should the size of the stackvector differs from the number
    ! of data sources - the error may be set depending on the 
    ! single_source flag.
    !
    implicit none

    ! Imported parametsr with intent IN
    type(mini_market_of_stacks), intent(inout) :: miniMarketPtr
    type(silja_wdr), intent(in) :: wdr
    logical, intent(in) :: ifMultiTime, ifSingleTime   ! what stacks to set

    ! Local declarations
    type(silja_stack), pointer :: stackPtr
    integer :: i, j, met_src_idx

    !
    ! Some checking first
    !
    if(.not.defined(wdr))then
      call set_error('Undefined wdr','set_met_srcs_in_sm_single_wdr')
      return
    end if

    if(ifMultiTime)then
      !
      ! Set the multiTime stack matrix
      !
      if(.not.miniMarketPtr%stack_multiTime_initialized)then
        call set_error('Non-initialized minimarket multiTime','set_met_srcs_in_sm_single_wdr')
        return
      end if

      ! If multiple met src stack
      ! As if_single_met_src_stack is set for the whole stack vector, only check the first one.
      ! If necessary to have single and multiple met src stacks in the same stackvector, the met src
      ! should be set met_src_missing for multi src stack  
      !
      stackPtr => miniMarketPtr%stacks_multiTime(1,1)%ptr
!      if(.not. fu_single_met_src_stack(stackPtr))then
!        ! set all met-srces to met-src missing.
!        do i=1,size(miniMarketPtr%stacks_multiTime,1)
!          do j=1,size(miniMarketPtr%stacks_multiTime,2)
!            stackPtr => miniMarketPtr%stacks_multiTime(i,j)%ptr
!            call set_met_src(stackPtr, met_src_missing)
!          end do  ! Scan over time dimension
!        end do  ! Scan over met_srcs dimension
!
!      else   ! If multiple met src stack

        if(size(miniMarketPtr%stacks_multiTime,1) /= fu_NbrOfMetSrcs(wdr))then
          call set_error('Stackvector size /= the number of data sources', &
                       & 'set_met_srcs_in_sm_single_wdr')
          return
        end if

      ! Stackvector is OK - just fill the data sources for all single-source stacks
        do i=1,size(miniMarketPtr%stacks_multiTime,1)
          do j=1,size(miniMarketPtr%stacks_multiTime,2)
            stackPtr => miniMarketPtr%stacks_multiTime(i,j)%ptr
            call set_met_src(stackPtr, fu_mds(wdr,i))
          end do  ! Scan over time dimension
        end do  ! Scan over met_srcs dimension
!      endif
    endif  ! if to set the multiTime stack

    !
    ! Now repeat the same exercise for the single time vector, if needed
    !
    if(ifSingleTime)then
      !
      ! Set the singleTime stack matrix
      !
      if(.not.miniMarketPtr%stack_singleTime_initialized)then
        call set_error('Non-initialized minimarket singleTime','set_met_srcs_in_sm_single_wdr')
        return
      end if

!      stackPtr => miniMarketPtr%stacks_singleTime(1)%ptr
!      if(.not. fu_single_met_src_stack(stackPtr))then
!        ! set all met-srces to met-src missing.
!        do i=1,size(miniMarketPtr%stacks_singleTime)
!          stackPtr => miniMarketPtr%stacks_singleTime(i)%ptr
!          call set_met_src(stackPtr, met_src_missing)
!        end do  ! Scan over met_srcs dimension
!
!      else   ! If multiple met src stack

        if(size(miniMarketPtr%stacks_singleTime) /= fu_NbrOfMetSrcs(wdr))then
          call set_error('Stackvector size /= the number of data sources', &
                       & 'set_met_srcs_in_sm_single_wdr')
          return
        end if

      ! Stackvector is OK - just fill the data sources for all single-source stacks
        do i=1,size(miniMarketPtr%stacks_singleTime)
          stackPtr => miniMarketPtr%stacks_singleTime(i)%ptr
          call set_met_src(stackPtr, fu_mds(wdr,i))
        end do  ! Scan over met_srcs dimension
!      endif

    endif  ! if to set the multiTime stack endif  ! if to set the multiTime stack

  end subroutine set_met_srcs_in_sm_single_wdr


  !**********************************************************************

  subroutine set_time_direction_sm(miniMarketPtr, ifReplaceEarliestTime)
    !
    ! Sets the rule for replacement of the stack time vector.
    ! If the run is forward, the earliest field is to be replaced, in adjoint run
    ! the latest one is the first to become outdated
    !
    implicit none

    type(mini_market_of_stacks), intent(inout) :: miniMarketPtr
    logical, intent(in) :: ifReplaceEarliestTime

    miniMarketPtr%replace_earliest = ifReplaceEarliestTime

  end subroutine set_time_direction_sm

  
  !*********************************************************************
  
  integer function fu_get_time_direction_sm(miniMarketptr)
    implicit none
    type(mini_market_of_stacks), intent(in) :: miniMarketPtr
    if(miniMarketPtr%replace_earliest)then
      fu_get_time_direction_sm = forwards
    else
      fu_get_time_direction_sm = backwards
    endif
  end function fu_get_time_direction_sm
  

  ! ***************************************************************

  ! ***************************************************************
  !
  !
  !        Put NWP data from external met_srcs in supermarket
  !
  !
  ! ***************************************************************
  ! ***************************************************************


  subroutine fill_minimarket_from_namelist(miniMarketPtr, &
                                         & nlIn, chItemTitle, & ! namelist and item name
                                         & shop_list, valid_time, st_time_feature, &
                                         & iUpdateType, &   ! create/update/...
                                         & wdr, &
                                         & storage_grid, &
                                         & iAccuracy, ifAdjustGrid, &
                                         & ifGotNewData)
    !
    ! Finds the necessary administrative data and gets the data to the prescribed stack.
    !
    ! ATTENTION.
    ! Stores ONLY to single_time_stack_flag
    !
    IMPLICIT NONE

    ! Imported parameters with intent IN:
    type(mini_market_of_stacks), pointer :: miniMarketPtr
    type(Tsilam_namelist), pointer :: nlIn
    character(len=*), intent(in) :: chItemTitle
    TYPE(silja_shopping_list), INTENT(in) :: shop_list
    type(silja_time), intent(in) :: valid_time  ! Used only for the template!
    type(silja_wdr), intent(in) :: wdr
    integer, intent(in) ::  iUpdateType, iAccuracy  ! indicator of the target stack, ...
    integer, intent(in) ::  st_time_feature  !
    type(silja_grid), intent(in) :: storage_grid
    logical, intent(in), optional :: ifAdjustGrid
    logical, intent(out) :: ifGotNewData

    ! Local variables
    INTEGER :: iFNm, nItems, iItem
    type(silam_sp), dimension(:), save, pointer :: fnames
    logical :: ifNewDataHere
    type(grads_template) :: fname_template
    type(Tsilam_nl_item_ptr), dimension(:), pointer :: ptrItems
    integer :: stack_flag
    type(silam_fformat) :: fform
    character(len=fnlen) :: sp
    logical :: ifAdjustLocal

    ifGotNewData = .false.

    ! Check input
    IF (.NOT.defined(shop_list)) THEN
      CALL set_error('shopping_list undefined','fill_minimarket_from_namelist')
      RETURN
    END IF
    if(.not. associated(nlIn))then
      call set_error('Not associated input namelist','fill_minimarket_from_namelist')
      return
    endif
    
    if(present(ifAdjustGrid))then
      ifAdjustLocal = ifAdjustGrid
    else
      ifAdjustLocal = .false.    ! it will be done this way at the next step anyway
    endif

    if (st_time_feature == dynamic_map) then
      stack_flag = multi_time_stack_flag
    else
      stack_flag = single_time_stack_flag
    endif

    ! Get the items with given title and scan the content to determine the 
    ! format of the input file and its name as decoded from template
    nullify(ptrItems)
    call get_items(nlIn, chItemTitle, ptrItems, nItems)
    if(error .or. nItems < 1)then
      call set_error(fu_connect_strings('Could not extract namelist items:',chItemTitle), &
                   & 'fill_minimarket_from_namelist')
      return
    endif

    ! Scan the items determining the list of input files and their format
    do iItem = 1, nItems
      
      sp = fu_content(ptrItems(iItem))
      fname_template = template_missing
      fform = fu_input_file_format(sp)
      if(error)return
      
      sp=adjustl(sp)
      sp=adjustl(sp(index(sp,' ')+1:))

      if(fform%iFormat == test_field_value_flag)then

          call msg('Storing test field:' + sp)
          CALL store_input_to_supermarket(miniMarketPtr, &
                                    & sp, &
                                    & fform, &
                                    & shop_list, &
                                    & wdr, &
                                    & iUpdateType, &
                                    & stack_flag, &  !
                                    & iAccuracy, &
                                    & ifNewDataHere, &
                                    & storage_grid = storage_grid, &
                                    & ifAdjust=ifAdjustLocal, &
                                    & st_time_feature = static_value, &
                                    & ifforcemetsrcacceptance = .TRUE., & !!Accept test field
                                    & ifVerbose_ = .true.)
          ifGotNewData = ifGotNewData .or. ifNewDataHere

      else
        nullify(fnames)
        call decode_template_string(sp, fname_template)
        call FNm_from_single_template(fname_template, &
                                      & valid_time, fnames, &
                                      & ifAdd = .false., &
                                      & ifStrict = .false., &
                                      & ifAllowZeroFcLen = .false., &
                                      & ifWait = .false., &
                                      & fc_step = fu_obstime_interval(wdr))
        IF (error) THEN
          CALL set_error('error1',  'fill_minimarket_from_namelist')
          return
        END IF

        ! Scan the given files and get the needed stuff
        !
        do iFNm = 1, size(fnames)
          if(fnames(iFNm)%sp=='')exit
          call msg('Reading the file:' + fnames(iFNm)%sp)
          CALL store_input_to_supermarket(miniMarketPtr, &
                                      & fnames(iFNm)%sp, &
                                      & fform, &
                                      & shop_list, &
                                      & wdr,& 
                                      & iUpdateType, &
                                      & stack_flag, &  ! ifUpdateAllowed
                                      & iAccuracy, &
                                      & ifNewDataHere, &
                                      & storage_grid = storage_grid, &
                                      & ifVerbose_ = .True., &
                                      & ifAdjust=ifAdjustLocal, &
                                      & st_time_feature = st_time_feature )
           ifGotNewData = ifGotNewData .or. ifNewDataHere
        end do   ! fnames
      endif   ! if a test field with single or real file
   end do   ! list of input items with their format

  end subroutine fill_minimarket_from_namelist


  !******************************************************************************

  SUBROUTINE fill_meteo_market(miniMarketPtr, wdr, shopping_list, timestep, iAccuracy, ifGotNewData)
    !
    ! 1. Checks if the data on the shopping list are in the supermarket
    ! 2. If not, finds the necessary administrative data and gets the
    ! data in here.
    !
    IMPLICIT NONE

    ! Imported parameters with intent IN:
    type(mini_market_of_stacks), pointer :: miniMarketPtr
    type(silja_wdr), intent(in) :: wdr
    TYPE(silja_shopping_list), INTENT(inout) :: shopping_list
    integer, intent(in) :: iAccuracy
    type(silja_interval), intent(in) :: timestep
    logical, intent(out) :: ifGotNewData

    ! Local declarations:
    TYPE(silja_time), DIMENSION(max_times) :: missing_obstimes, received_obstimes
    INTEGER :: i, j, iTempl, iFNm, iReceivedTimes, nFileNames
    type(meteo_data_source) :: met_src
    TYPE(silja_shopping_list) :: list_of_missing_stuff
    LOGICAL :: data_exists_already, ifNewDataHere
    type(grads_template), dimension(:), pointer :: fname_templates
    logical :: ifAllowZeroFcLen, ifWait
    real :: fSign
    INTEGER ::  direction
    ! Two arrays that remember their size between the calls to avoid reallocation
    type(silam_sp), dimension(:), save, pointer :: fnames
    type(silam_fformat), dimension(:), save, pointer :: fform

    ifGotNewData = .false.

    !
    ! Check input.
    !
    IF (fu_fails(miniMarketPtr%stack_multiTime_initialized .or. &
               & miniMarketPtr%stack_singleTime_initialized, 'minimarket not init','fill_meteo_market'))RETURN
    IF (fu_fails(defined(shopping_list), 'shopping_list undefined','fill_meteo_market'))RETURN

    fSign = SIGN(1.0,fu_sec(timestep))
    if (fSign  > 0) then
      direction = forwards
    else
      direction = backwards
    endif
      
    !
    ! Check if some or all of required data is already in the supermarket.

    CALL check_sm_for_data_in_list(miniMarketPtr, &
                                 & shopping_list, &
                                 & data_exists_already,&
                                 & list_of_missing_stuff,&
                                 & missing_obstimes, &
                                 & wdr)
    IF (error) then
      call msg('*--- list_of_missing_stuff ----*')
      CALL report(list_of_missing_stuff)
      call msg("list_of_missing_stuff")
      return
    endif

    IF (data_exists_already) THEN
!!!$      IF (supermarket_info)&
!!!$          & PRINT *, 'desired data already in supermarket OK!'
      RETURN
    ELSE
!!!$      IF (supermarket_info)&
!!!$          & PRINT *, 'Desired data not in supermarket, retireving..'
    END IF

    !
    ! Get data to supermarket.
    !
    iReceivedTimes = 0
    received_obstimes(iReceivedTimes + 1) = time_missing
    met_src = fu_met_src(list_of_missing_stuff)
    fname_templates => fu_fname_templates(wdr)
    ifAllowZeroFcLen = fu_if_zero_fc_len_allowed(wdr)
    ifWait = fu_ifWait(wdr)
    IF (error) RETURN
    !
    ! Template is allowed to have wildcards or imply test field or whatever.
    ! So, there may be several files for each time and each template.
    !
    times: DO i = 1, SIZE(missing_obstimes)

      IF (.NOT.defined(missing_obstimes(i))) EXIT
      !
      ! For this missing time, try to find the data. Can also take the next time slots 
      ! if the current one is missing (if standard setup namelist allowed this trick)
      !
      call FNm_from_template_array(fu_fname_templates(wdr), &
                                 & fu_file_formats(wdr), &
                                 & missing_obstimes(i), &
                                 & fnames, nFileNames, &     ! actually, output
                                 & fform, &                  ! output too, same size as fnames
                                 & time_missing, &      ! anal_time, &
                                 & .false., &           ! ifStrict = 
                                 & ifAllowZeroFcLen, &  ! ifAllowZeroFcLen = 
                                 & ifWait, &            ! ifWait = 
                                 & fu_obstime_interval(wdr), &
                                 & fu_max_hole_length(wdr))
      !
      ! Did we get anything for this time slot?
      !
      if(nFileNames == 0)then
        !
        ! This time slot has no data. Go to the next one until
        ! the data are found or max length of the hole is exceeded
        !
        ! Note that the missing_obstimes can ebd earlier than the data are found
        !
        if(.not. defined(missing_obstimes(i+1))) then
          missing_obstimes(i+1) = fu_closest_obstime(missing_obstimes(i) + one_second * fSign, &
                                                   & direction, fu_obstime_interval(wdr))
          CALL fix_shopping_time_boundaries(shopping_list, &
                        & fu_earliest_time((/fu_start_time(shopping_list), missing_obstimes(i+1)/)), &
                        & fu_latest_time((/fu_end_time(shopping_list), missing_obstimes(i+1)/)))
          CALL fix_shopping_time_boundaries(list_of_missing_stuff, &
                        & fu_earliest_time((/fu_start_time(list_of_missing_stuff), missing_obstimes(i+1)/)), &
                        & fu_latest_time((/fu_end_time(list_of_missing_stuff), missing_obstimes(i+1)/)))
        endif
        if(error)return
        
        if(iReceivedTimes == 0)then
          !
          ! No times were good this-far. But we know that previous call
          ! succeeded, i.e. it found something: call never ends with the last time empty.
          ! So, just go ahead and find the right time. To make it safe, still check that we do 
          ! not go to infinity 
          !
          if(fu_abs(missing_obstimes(i+1) - missing_obstimes(1)) > fu_max_hole_length(wdr))then
            call msg_warning('Too large hole in the meteo data1','fill_meteo_market')
            call msg('Last tried time:' + fu_str(missing_obstimes(i+1)))
            call set_error('Too large hole in the meteo data1','fill_meteo_market')
            return
          endif
        else
          !
          ! Having the last good received time, check that we do not go too far from it
          !
          if(fu_abs(missing_obstimes(i+1) - received_obstimes(iReceivedTimes)) > &
                                                            & fu_max_hole_length(wdr))then
            call msg_warning('Too large hole in the meteo data2','fill_meteo_market')
            call msg('Last tried time:' + fu_str(missing_obstimes(i)))
            call set_error('Too large hole in the meteo data2','fill_meteo_market')
            return
          endif
        endif  ! if no times were received in this round
      else
        !
        ! Valid data were found. Register this obstime to the list of good ones
        !
        iReceivedTimes = iReceivedTimes + 1
        received_obstimes(iReceivedTimes) = missing_obstimes(i)
        received_obstimes(iReceivedTimes + 1) = time_missing
      endif  ! if any files obtained for this time
        
      
      !!!!
      !!!! Have to scan tempaltes and either create file name or forward the template string
      !!!!
      !!!do iTempl = 1, size(fname_templates)
      !!!  
      !!!  if(fname_templates(iTempl) == template_missing) exit
      !!!
      !!!  fform = fu_file_format(wdr,iTempl)
      !!!
      !!!  if(fu_if_input_file(fform))then
      !!!    !
      !!!    ! File(s) is(are) implied - find their names
      !!!    !
      !!!    call FNm_from_single_template(fname_templates(iTempl), &
      !!!                                & missing_obstimes(i), &
      !!!                                & fnames, &
      !!!                                & ifAdd = .false., &
      !!!                                & fc_step = fu_obstime_interval(wdr))
      !!!    IF (error) THEN
      !!!      CALL msg_warning('no data for this time, but continuing...')
      !!!      CALL unset_error('fill_meteo_market')
      !!!      CYCLE
      !!!    END IF
      !!!  else
      !!!    !
      !!!    ! Template stores a request of some non-file input in its collection
      !!!    !
      !!!    call enlarge_array(fnames,1)
      !!!    fnames(1)%sp = fu_collection(fname_templates(iTempl))
      !!!    if(size(fnames) > 1) fnames(2)%sp = ''
      !!!  endif  ! If this is a file

        do iFNm = 1, nFileNames  !size(fnames) or less

!          if(fnames(iFNm)%sp=='')exit
!          call report(fu_stack(miniMarketPtr,met_src_missing))
!          call report(fu_sm_simple_field(miniMarketPtr, met_src_missing, ground_pressure_flag, level_missing, single_time_stack_flag))
!          if (error) call unset_error("Tried to report..")
!          call msg("-----------BEFORE Store--------"+ fu_name(fu_stack(miniMarketPtr,met_src_missing)))

          CALL store_input_to_supermarket(miniMarketPtr, &
                                        & fnames(iFNm)%sp, &
                                        & fform(iFNm), & ! One format for all these files
                                        & list_of_missing_stuff, &
                                        & wdr, &
                                        & create_if_absent, &    ! ifUpdateAllowed
                                        & multi_time_stack_flag , iAccuracy, &
                                        & ifNewDataHere, &
                                        & ifadjust = .TRUE.,  & !Force storage grid
                                        & storage_grid = meteo_grid)
          ifGotNewData = ifGotNewData .or. ifNewDataHere
!          call report(fu_sm_simple_field(miniMarketPtr, met_src_missing, ground_pressure_flag, level_missing, single_time_stack_flag))
!          if (error) call unset_error("Tried to report..")
!          !call report(fu_stack(miniMarketPtr,met_src_missing))
!          call msg("-----------AFTER Store--------"+ fu_name(fu_stack(miniMarketPtr,met_src_missing)))

        end do   ! size(fnames)
      !!!end do ! Template array
    
      CALL arrange_supermarket(miniMarketPtr, .false.) !Do not look for surface pressure yet..
      IF (error) RETURN
      
    END DO times ! Missing times

#ifdef PREFETCH
    !!!Generate next filename for last defined missing_obstimes
    if (defined(missing_obstimes(i-1)) .and. smpi_global_rank==0) then
      call FNm_from_template_array(fu_fname_templates(wdr), &
                                 & fu_file_formats(wdr), &
                                 & missing_obstimes(i-1)+fu_obstime_interval(wdr), &
                                 & fnames, nFileNames, &     ! actually, output
                                 & fform, &                  ! output too, same size as fnames
                                 & time_missing, &      ! anal_time, &
                                 & .false., &           ! ifStrict = 
                                 & ifAllowZeroFcLen, &  ! ifAllowZeroFcLen = 
                                 & .False., &            ! ifWait = 
                                 & fu_obstime_interval(wdr), &
                                 & fu_max_hole_length(wdr))
        call msg("Calling prefetch for:")
        do iFNm = 1, nFileNames  !size(fnames) or less
                call msg(fnames(iFNm)%sp)
                call EXECUTE_COMMAND_LINE("cat "//trim(fnames(iFNm)%sp)//">/dev/null",&
                                       & wait=.false. )
        end do   ! size(fnames)
        call msg("=====")
    endif
#endif
      
    !
    ! Maybe print some supermarket_info and exit.
    !
    IF (supermarket_info) THEN
      CALL supermarket_times(miniMarketPtr, met_src,missing_obstimes,i) ! missing_obstimes, i are used temporarily
      IF (error) RETURN
      call msg('********************************************************')
      call msg('Data in supermarket now for the following times:')
      CALL report(missing_obstimes)
!      call report(miniMarketPtr)
    END IF

  END SUBROUTINE fill_meteo_market


  ! ***************************************************************


  SUBROUTINE store_input_to_supermarket(miniMarketPtr, &     ! Mini market to put in
                                      & chFName, fformat, & ! File where from to take it
                                      & shopping_list, &     ! The stuff to search
                                      & wdr, &               !Tweak input if needed
                                      & iUpdateType, &       ! Overwrite or not existing fields
                                      & stack_type, &        ! stationary or time-dependent
                                      & iAccuracy, &         ! of reporjection
                                      & ifNewDataConsumed, &  ! OUTPUT flag
                                      & stack_indices_, &  ! optional. OVERWRITES met_src, stacks are picked directly
                                      & storage_grid, &       ! storage grid, grid_missing for stack default
                                      & iBinary, &        ! optional: points at GrADS/NetCDF if already open
                                      & ifForceMetSrcAcceptance, ifadjust, &
                                      & st_time_feature, & ! Tweak time properties of the field ID in singletime stack
                                      &ifVerbose_)

    ! 
    ! Finds fields in the given file and puts them to
    ! supermarket. Decodes GRIB, NetCDF, GrADS, reads ASCII.
    ! 
    ! Only field that are found in shopping list are accepted.
    ! Here it is also ok to use a shopping list with accept_any_met_src,
    ! when the file is a "black box". See module shopping_lists.
    !
    ! Data may be stored in its original grid, user-defined grid
    ! or user-defined sub area of the original grid. For this
    ! use set_supermarket_storage_grid and fu_supermarket_storage_area.
    !
    ! The file may contain fields from several valid times or analysis.
    !
    ! If an unknown field is found, it is skipped (no error set)
    ! 
    use omp_lib_kinds, only : OMP_lock_kind

    IMPLICIT NONE

    ! Imported parameters with intent IN:
    type(mini_market_of_stacks), pointer :: miniMarketPtr
    CHARACTER (LEN=*), INTENT(in) :: chFName
    type(silam_fformat), intent(in) :: fformat
    TYPE(silja_shopping_list), INTENT(in), target :: shopping_list
    type(silja_wdr), intent(in) :: wdr
    integer, INTENT(in) :: stack_type ! dynamic/permanent meteo or dispersion
    integer, intent(in) :: iUpdateType  ! create, overwrite, add, etc - see fields module
    integer, intent(in) :: iAccuracy
    logical, intent(out) :: ifNewDataConsumed
    type(silja_grid), intent(in) :: storage_grid  !Storage grid or grid_missing 
        ! to use   fu_storage_grid(fu_wdr(stack_to_fill)) 
    integer, intent(inout), optional :: iBinary
    integer, intent(in), optional :: st_time_feature
    logical, intent(in), optional :: ifForceMetSrcAcceptance, ifAdjust, ifVerbose_
    integer, dimension(:), intent(in), target, optional :: stack_indices_ ! For boundary market needs, 
                                                              ! allow one field to go to differerent 
                                                              ! stacks with different grids

    ! Local declarations:
    character(len=*), parameter :: sub_name = "store_input_to_supermarket"
    REAL, DIMENSION(:), pointer :: work_array
    INTEGER, DIMENSION(:), pointer :: work_int_array, indexVars, ordermsgs
    INTEGER :: input_unit, i, fields_found, fields_accepted, indFldGrib, indStart, indEnd
    LOGICAL :: eof, ifForceMetSrcAcceptanceLocal, found, ifAdjustLocal, ifVerbose
    TYPE(silja_field_id) :: id, idStore
    
    real ::  scale_factor, fMissingValue
    integer (kind=8) :: before, after
    logical :: ifMapNeeded
    TYPE(silja_field_id), dimension(:), pointer :: idList  
    TYPE(Tfld_ptr), dimension(max_2d_fields) :: MarketFldPtrs
    integer :: iStack, nStacks, number_of_times
    integer :: TweakTime
    integer :: indexVar, nIds, iID, it,  ilev, times_found, io_stat, iTmp, jTmp, iVar, fs_grid, nOrder
    TYPE(silja_shopping_list), pointer :: shpLstPtr
    type(silja_time), dimension(:), pointer ::  timeLst
    type(silja_time), dimension(max_times) ::  valid_times  !Should be more than n of times in multitime stack
    type(silja_time) ::  timeTmp, shopStartTime, time_to_force
    type(silja_interval) :: obstime_interval, intervalTmp
    integer, dimension(:), pointer :: stack_indices ! For boundary market needs,
    integer(omp_lock_kind), dimension(max_threads) :: lockarr
    integer :: old_unit
    integer :: quant, nx, ny, n_out_of_range, n_patched, n_failed

    integer(kind=8) :: count0, count1, count2, count3


    ifNewDataConsumed = .false.

    if(present(st_time_feature))then
      TweakTime=st_time_feature
    else
      TweakTime=dynamic_map !Force (should we  assert it?) IDs to be instant
    endif


    if(present(ifVerbose_))then
      ifVerbose = ifVerbose_
    else
      ifVerbose = .false.
    endif
    shpLstPtr => shopping_list

!    call msg("shopping list for store_input_to_supermarket")
!    call report(shopping_list)
!    call msg("end shopping list for store_input_to_supermarket")


    if(present(ifForceMetSrcAcceptance))then
      ifForceMetSrcAcceptanceLocal = ifForceMetSrcAcceptance
    else
      ifForceMetSrcAcceptanceLocal = .false.
    endif
    if(present(ifAdjust))then
      ifAdjustLocal = ifAdjust
    else
      ifAdjustLocal = .false.
    endif

    if (.not. ifAdjustLocal) then
        call set_error("Grid flexibility not forbidden, surprises might appear", sub_name)
        return
    endif

    call SYSTEM_CLOCK(before)
    
    ! First some checking
    !
    if( present(stack_indices_))then
      call msg('Storing to supermarket '//trim(miniMarketPtr%name)//' with stack_indices: ' // trim(chFName))
      stack_indices => stack_indices_
      do i = 1, size(stack_indices)
        if(stack_indices(i) == int_missing)then 
          nStacks = i-1
          exit
        endif
      enddo
    else
      nStacks = 1
      if(stack_type == single_time_stack_flag) then
        call msg('Storing to supermarket '//trim(miniMarketPtr%name)//' singletime: '//trim(chFName))

        if (.not. miniMarketPtr%stack_singleTime_initialized)then
           call set_error('Single-time stack not initialized for:' + miniMarketPtr%name, sub_name)
           return
        endif
      endif
      if(stack_type == multi_time_stack_flag) then
        call msg('Storing to supermarket '//trim(miniMarketPtr%name)//' multitime: ' // trim(chFName))

        if (.not. miniMarketPtr%stack_multiTime_initialized)then
          call set_error('Multi-time stack not initialized for:' + miniMarketPtr%name, sub_name)
          return
        endif
      endif
      stack_indices => null()
    endif

    !Is this check needed at all?????
    if(nStacks < 1 .or. nStacks > 100)then
      call msg('Funny number of stacks:', nStacks)
      call set_error('Funny number of stacks:', sub_name)
      return
    endif

    !
    ! Start from opening file. 
    ! As one file can include several timesteps, for grads and netcdf files try following:
    ! First check, if it's already open. If not, check if exists. If exists, open; 
    ! if not, look for required fields in binary already open.
    !


    old_unit = int_missing
    if(present(iBinary)) old_unit = iBinary

    obstime_interval =  fu_obstime_interval(wdr)

    call open_input_file(chFName, FFormat, input_unit, old_unit, obstime_interval)
    if(error)return

    if(present(iBinary)) iBinary = input_unit

    !
    ! 2. Go through the file.
    !
    fields_found = 0
    fields_accepted = 0
    eof = .false.

    shopStartTime = fu_start_time(shopping_list)

    select case (FFormat%iformat)
    case (netcdf_file_flag)
      !
      !====================== NETCDF =================================
      !
      ! NetCDF has a list of headers available
                                                        !      ifFirstTime, ifFirstLev
      call id_list_from_netcdf_file(input_unit, idList, fields_found)
      call timelst_from_netcdf_file(input_unit, timeLst, times_found)
      call supermarket_times(miniMarketPtr, met_src_missing, valid_times, number_of_times)
  
      do it = 1, times_found
        if(.not. fu_time_in_list(timeLst(it), interval_missing, shopping_list))cycle
        !!! Shuld be integer number of timesteps from
        if (defined(obstime_interval)) then
          if ( abs(modulo((timeLst(it)-shopStartTime)/obstime_interval + 0.5, 1.) - 0.5) > 0.001 ) cycle
        endif

        if(fu_list_time_indicator(shpLstPtr) == accept_same_month) then   
          found = .false.
          do i = 1, number_of_times    !times already in supermarket
            if(fu_mon(timeLst(it)) == fu_mon(valid_times(i)))then
              found = .true.
              exit
            endif
          enddo
          if(found)cycle
        endif

        do i = 1, fields_found

          call set_valid_time(idList(i), timeLst(it))
          if(fu_accumulated_quantity(fu_quantity(idList(i))))then
            call set_accumulation_length(idList(i), timeLst(it)-fu_analysis_time(idList(i)))
          endif

          IF (fu_field_id_in_list(idList(i), shopping_list, indexVar)) THEN
            !
            ! There can be some files, which have more than worksize elements
            ! For such case, have to allocate a special array
            !
            fs_grid=fu_number_of_gridpoints(fu_grid(idList(i)))
            work_array => fu_work_array(fs_grid)
            if(error)return

            fMissingValue = fu_real_missing_replacement(fu_quantity(idList(i)))
            call field_indices_from_netcdf_file(input_unit, idList(i), idStore, iVar, iTmp, ilev)
            if(ifVerbose)then
              call msg('Reading netcdf field, iVar, it, ilev', (/iVar, iTmp, ilev/))
              call msg("id_requested")
              call report(idList(i))
              call msg("id_read")
              call report(idStore)
            endif
            call read_field_from_netcdf_indices(input_unit, ivar, iTmp, ilev, &
                                           & work_array, &
                                           & fill_value = fMissingValue)

#ifdef DEBUG
            !
            !Dirty hack for bounaries
            if (nstacks> 3) then
                    if (.not. all( work_array(1:fs_grid) >=0)) then
                         do jTmp=1,fs_grid
                            if (.not.  work_array(jTmp) >=0) then
                              call msg("Negatibe value at index:", jTmp, fs_grid)
                              exit
                            endif
                         enddo 
                         call msg("Values nearby",  work_array(max(1,jTmp-5):min(fs_grid,jTmp+5)) )
                          call set_error("Some negatives for boundaries..","store_input_to_supermarket")

                         call read_field_from_netcdf_file(input_unit, &
                                         & idList(i),&
                                         & work_array, fill_value = fMissingValue) !, &
                         return
                   endif
            endif
#endif 
          
           call  put_field_to_sm(miniMarketPtr, &     ! Mini market to put in
                                     & nStacks, &
                                      & idStore, &
                                      & work_array, &
                                      & iUpdateType, &       ! Overwrite or not existing fields
                                      & stack_type, &        ! stationary or time-dependent
                                      & iAccuracy, &         ! of reporjection
                                      & stack_indices, &  ! OVERWRITES met_src, stacks are picked directly
                                      & storage_grid, &       ! storage grid
                                      & ifForceMetSrcAcceptanceLocal, TweakTime, ifAdjustLocal, &
                                      & fu_if_randomise(wdr), &
                                      & shpLstPtr, indexVar, time_to_force)

            if(ifVerbose)call msg('Accepted netcdf field total: ', work_array(1:10))
            if (.not. error) fields_accepted = fields_accepted + 1
            
            call free_work_array(work_array)  

          endif   ! if field is accepted
        enddo   ! fields found
      enddo   ! times found
      deallocate(idList)


    !
    !====================== GRADS =================================
    !
    case (grads_file_flag) !fformat%iformat
      !
      ! GrADS also has similar type of logic for picking the fields
      !
      call supermarket_times(miniMarketPtr, met_src_missing, valid_times, number_of_times)
      if(error)return
      
      indStart = fu_grads_time_index(input_unit, &
                                   & fu_start_time(shopping_list), back_and_forwards, &
                                   & fu_list_time_indicator(shpLstPtr) == accept_same_month)  ! ifAcceptSameMonth
      indEnd = fu_grads_time_index(input_unit, &
                                 & fu_end_time(shopping_list), back_and_forwards, &
                                 & fu_list_time_indicator(shpLstPtr) == accept_same_month)  ! ifAcceptSameMonth
      if(indStart == int_missing .and. indEnd == int_missing)then
        call set_error('Undefined time indices for the requested times. Start=' + &
                                             & fu_str(fu_start_time(shopping_list)) + &
                                             & ', end=' + fu_str(fu_end_time(shopping_list)), &
                      &'store_input_to_supermarket')
        return
      elseif(indStart == int_missing)then
        indStart = indEnd
      elseif(indEnd == int_missing)then
        indEnd = indStart
      endif
      !
      ! Scan the required interval from the sopping list, taking care of accept_same_month
      !
      do iT = indStart, indEnd
        timeTmp = fu_time_of_grads(input_unit, iT)

        !!! Shuld be integer number of timesteps from

        if (defined(obstime_interval)) then
          if ( abs(modulo((timeTmp-shopStartTime)/obstime_interval + 0.5 ,1.) - 0.5) > 0.001 ) cycle
        endif
        !
        ! Check that the existing fields do not cover the needed time - valid check for monthly fields
        !
        if(fu_list_time_indicator(shpLstPtr) == accept_same_month) then   
          found = .false.
          do i = 1, number_of_times    !times already in supermarket => field is already consumed
            if(fu_mon(timeTmp) == fu_mon(valid_times(i)))then
              found = .true.
              exit
            endif
          enddo
          if(found)cycle   ! that needed time is actually not needed
        endif
        
        do iVar = 1, fu_n_gvars(input_unit)
          do iLev = 1, fu_n_gVar_levs(input_unit, iVar)
            call get_grads_var_metadata(input_unit, iVar, iLev, it, id)
            fields_found = fields_found + 1
            if(fu_field_id_in_list(id, shopping_list, &
                                 & indexVar, &
                                 & fu_data_time_features(input_unit), &
                                 & time_to_force))then        ! out: for monthly fields we can force valid time
              if(ifVerbose)then
                call msg('Reading grads field, time index', it)
                call report(id)
              endif

              fs_grid=fu_number_of_gridpoints(fu_grid(id))

              work_array => fu_work_array(fs_grid)
              if(error)return
              call read_field_from_grads_indices(input_unit, &
                                       & iVar, iLev, it, &
                                       & work_array) !, fill_value_ = ???0.) ! the grid-data
              call  put_field_to_sm(miniMarketPtr, &     ! Mini market to put in
                                     & nStacks, &
                                      & id, &
                                      & work_array, &
                                      & iUpdateType, &       ! Overwrite or not existing fields
                                      & stack_type, &        ! stationary or time-dependent
                                      & iAccuracy, &         ! of reporjection
                                      & stack_indices, &  ! OVERWRITES met_src, stacks are picked directly
                                      & storage_grid, &       ! storage grid
                                      & ifForceMetSrcAcceptanceLocal, TweakTime, ifAdjustLocal, &
                                      & fu_if_randomise(wdr), &
                                      & shpLstPtr, indexVar, time_to_force)
              if(ifVerbose)call msg('Accepting grads field total: ', work_array(1:10))

              if (.not. error) fields_accepted = fields_accepted + 1
              call free_work_array(work_array)  ! can be a huge array, so must be freed here
            endif  ! if field in the list
          enddo ! iLev
        enddo  ! iVar
        if(stack_type == single_time_stack_flag .and. fields_accepted > 0)exit  !  just one time needd
      enddo ! iTime

    !
    !======================  TEST FIELD  ================================
    !
     case(test_field_value_flag)        
       !
       ! For the test field, the grid must be selected to be something. 
       !
       work_array => fu_work_array()
       timeTmp = shopStartTime
       call make_test_field(chFName, &  ! chFName contain the field name and value
                          & id, & !idStore, &
                          & work_array, &
                          & shopStartTime, &
                          & zero_interval, & !
                          & storage_grid)
       if(error)return
       eof = .true. ! Exit at the next cycle
!       id = idStore
       
!       call set_grid(id,coarse_geo_global_grid)  ! Dirty hack: shoplist for test field has 
                                  !always coarse_geo_global_grid (see put_test_field_to_content)
                                  ! Thus we have separate id (to match the shoplist)
                                  ! and idStore (to store)
                                  
       if (stack_type == single_time_stack_flag) then
        times_found = 1
        intervalTmp = zero_interval
       else
         if (defined(obstime_interval)) then
           times_found = int((fu_end_time(shopping_list) -  shopStartTime ) /obstime_interval) + 1
           intervalTmp = obstime_interval
         else
           times_found = 1
           intervalTmp = zero_interval
         endif
       endif
       fields_found = times_found
       do iT = 1, times_found
         timeTmp = shopStartTime + intervalTmp*(iT-1)
         call set_valid_time(id,      timeTmp )
!         call set_valid_time(idStore, timeTmp )
         call msg('The following test field has been created (store input)')
         call report(id) !Store)
         call msg('')
         if(fu_field_id_in_list(id, shopping_list,  indexVar))then  
            call  put_field_to_sm(miniMarketPtr, &     ! Mini market to put in
                                & nStacks, &
                                & id, &  !idStore, &
                                 & work_array, &
                                 & iUpdateType, &       ! Overwrite or not existing fields
                                 & stack_type, &        ! stationary or time-dependent
                                 & iAccuracy, &         ! of reporjection
                                 & stack_indices, &  ! OVERWRITES met_src, stacks are picked directly
                                 & storage_grid, &       ! storage grid
                                 & ifForceMetSrcAcceptanceLocal, TweakTime, ifAdjustLocal, &
                                 & fu_if_randomise(wdr), &
                                 & shpLstPtr, indexVar, time_to_force)
           if (.not. error) fields_accepted = fields_accepted + 1
         endif
         if(error)return
       enddo
       call free_work_array(work_array) ! can be a huge array, so must be freed here
       
    !
    !====================== GRIB =================================
    !
    case(grib_file_flag)
        
      CALL SYSTEM_CLOCK(count0)
      call  id_list_from_grib_file(input_unit, idList, nIDs)
      IF (error) THEN
           CALL msg_warning('error during GRIB-decoding',&
                          & 'store_input_to_supermarket')
           return
      endif
      CALL SYSTEM_CLOCK(count1)
      fields_found = nIDs
      work_int_array => fu_work_int_array(2*nIDs) 
      indexVars => work_int_array(1:nIDs)
      orderMsgs => work_int_array(nIDs+1:2*nIDs)
      nOrder = 0
      fs_grid = 0 !Max field size
      ifMapNeeded = .False.

      do iID = 1, nIDs  !Select messages and create storage for them
        if(.not. defined(idList(iID))) cycle

        if(ifVerbose)then
          call msg('Found grib field:')
          call report(idList(iID))
        endif
  
        IF (fu_field_id_in_list(idList(iID), shopping_list, indexVar, data_time_features=TweakTime)) THEN ! Accept the field
          nOrder = nOrder + 1
          indexVars(nOrder) = indexVar
          orderMsgs(nOrder) = iID
          if (indexVar >0) then ! there is an explicit var in the shoplist
            if (fu_if_var_mapped(shopping_list, indexVar)) ifMapNeeded = .true.
          endif
          fs_grid = max(fs_grid, fu_number_of_gridpoints(fu_grid(idList(iID))))

          ! Get final ID FIXME: 1:1 mapping assumed, i.e. not for boundaries
          call  inputID2MarketID(idList(iID), id, &
            & stack_type, storage_grid, TweakTime,  fuIfStaggeredWindsWDR(wdr), shopping_list,  indexVar, time_to_force)
         ! Tfld_ptr
         ! MarketFldPtrs
          call find_field_data_storage_2d(miniMarketPtr, id, stack_type, &
                  & MarketFldPtrs(nOrder)%ptr, MarketFldPtrs(nOrder)%IDptr )

          if(ifVerbose) call msg("Accepting:"+fu_quantity_string(fu_quantity(idList(iID))))
        else
          if(ifVerbose) call msg("NOT Accepting:"+fu_quantity_string(fu_quantity(idList(iID))))
        endif
      enddo
      CALL SYSTEM_CLOCK(count2)

      
      if (defined(storage_grid) .and. (.not. ifMapNeeded)) then
        ! Simple 1:1 mapping, known target grida
#ifdef VS2012 
        ! No parallel GRIB handling at Windows
        !$OMP PARALLEL if (.false.) default(none) &  
#elif defined VOIMA_ECCODES_BUG
        !$OMP parallel if (.false.) default(none) &
#else   
        !$OMP parallel if (Norder > 4) default(none) &
#endif
        !$OMP & shared(nOrder,orderMsgs,indexVars,idList,error, fs_grid,&
        !$OMP input_unit,MarketFldPtrs, wdr, iAccuracy, lockarr) &
        !$OMP & private(work_array,itmp,iID, indexVar, &
        !$OMP & quant, nx, ny, n_out_of_range, n_patched, n_failed)

        ! These locks are needed to generate interpolation structures 
        !$ call OMP_INIT_LOCK(lockarr(omp_get_thread_num()+1))

        !$OMP BARRIER
        work_array => fu_work_array(fs_grid)
        !$OMP DO
        do itmp = 1, nOrder  !Get data
            if (error) cycle
            iID = orderMsgs(iTmp)
            indexVar = indexVars(iTmp)
            ! storage_grid

            call read_field_from_grib_file(input_unit, idList(iID), iID, work_array)
            if (error) then
              call set_error("error after read_field_from_grib_file in par",sub_name)
              cycle
            endif
            quant=fu_quantity(idList(iID)) !!Quant. from grib

            call  remap_field( fu_grid(idList(iID)), &
                        & work_array,&
                        & fu_grid(MarketFldPtrs(iTmp)%IDptr),&
                        & marketFldPtrs(iTmp)%ptr, &
                        & fu_if_randomise(wdr), &
                        & iAccuracy, &
                        & fu_regridding_method(quant), &
                        & fu_outGridInterpType(quant), &
                        & lockarr)
           if(error) then
            call msg_warning("error after remap_field in store_input_to_supermarket. The field:",sub_name)
            call report(idList(iID))
            call msg('target grid:')
            call report(fu_grid(MarketFldPtrs(iTmp)%IDptr))
            call set_error("error after remap_field in par",sub_name)
           endif
           !
           !Final check after remap 
           call grid_dimensions(fu_grid(MarketFldPtrs(iTmp)%IDptr), nx, ny)
           call check_quantity_range(quant, marketFldPtrs(iTmp)%ptr, nx*ny, nx, &
              & .false., .true., &  ! ifRequireValidity, ifSilent 
             & n_out_of_range, n_patched, n_failed)

           if(error) then
            call set_error("error after check_quantity_range in par",sub_name)
           endif
#ifndef VIOMA_BUG
#ifdef DEBUG
           !Triggers voima compiler bug
          if (n_out_of_range > 0) then
            call msg('Out-of-range values removed:' + fu_quantity_string(quant), n_out_of_range)
          end if
#endif
#endif
          
        
                  
    
         END DO ! loop over inputs 
         !$OMP END DO
         call free_work_array(work_array) 

         !$OMP BARRIER
         !$ call OMP_DESTROY_LOCK(lockarr(omp_get_thread_num()+1))
         !$OMP END parallel


         fields_accepted = nOrder


      else
        work_array => fu_work_array(fs_grid)
        do itmp = 1, nOrder  !Get data
            iID = orderMsgs(iTmp)
            indexVar = indexVars(iTmp)

            call read_field_from_grib_file(input_unit, idList(iID), iID, work_array)
            
            call   put_field_to_sm(miniMarketPtr, &     ! Mini market to put in
                               & nStacks, &
                               & idList(iID), &
                               & work_array, &
                               & iUpdateType, &       ! Overwrite or not existing fields
                               & stack_type, &        ! stationary or time-dependent
                               & iAccuracy, &         ! of reporjection
                               & stack_indices, &  ! OVERWRITES met_src, stacks are picked directly
                               & storage_grid, &       ! storage grid
                               & ifForceMetSrcAcceptanceLocal, TweakTime, ifAdjustLocal, &
                               & fu_if_randomise(wdr), &
                               & shpLstPtr, indexVar, &
                               & time_to_force)
            
            if (fu_fails(.not. error, "put_field_to_sm broke", sub_name)) return
              
            fields_accepted = fields_accepted + 1
    
         END DO ! loop over inputs 
         call free_work_array(work_array) !
       endif
       call free_work_array(work_int_array) ! 
      CALL SYSTEM_CLOCK(count3)
      call msg("Storing (IDlist, find storage, stores)  took ms", &
        &(/ int(1e-6*(count1-count0)),  int(1e-6*(count2-count1)), int(1e-6*(count3-count2))/))
      

      case(ascii_file_flag)
      
        DO while(.not.eof)
        !
        ! ASCII file - just erad the next header, together with the data
        !
        work_array => fu_work_array()
        CALL read_next_field_from_ascii_file(input_unit, &
                                               & eof,&
                                               & id,&
             !                                               & fMissingValue, & ! the grid-data
                                               & work_array) ! ASCII is not allowed to have too large array
        if(eof)exit

        IF (error) THEN
          CALL msg_warning('error while getting new field, but continuing.',&
                                 & 'store_input_to_supermarket')
          CALL unset_error('store_input_to_supermarket4')
            IF (eof) EXIT
            CYCLE
          END IF
          IF (fu_field_id_in_list(id, shopping_list, indexVar)) THEN ! Accept the field
          call put_field_to_sm(miniMarketPtr, &     ! Mini market to put in
                             & nStacks, &
                             & id, &
                             & work_array, &
                             & iUpdateType, &       ! Overwrite or not existing fields
                             & stack_type, &        ! stationary or time-dependent
                             & iAccuracy, &         ! of reporjection
                             & stack_indices, &  ! OVERWRITES met_src, stacks are picked directly
                             & storage_grid, &       ! storage grid
                             & ifForceMetSrcAcceptanceLocal, TweakTime, ifAdjustLocal, &
                             & fu_if_randomise(wdr), &
                             & shpLstPtr, indexVar, &
                             & time_to_force)

            if (.not. error) fields_accepted = fields_accepted + 1

            call free_work_array(work_array) ! can be a huge array, so must be freed here

          END if  ! if field is in the list

          if (.not. error) fields_found = fields_found + 1

          IF (error) call unset_error('store_input_to_supermarket5')
  
        END DO ! loop over input
      case default
      call msg('Non-supported file type', FFormat%iformat)
      call set_error('Non-supported file type','store_input_to_supermarket')
    end select

    !
    ! All done, close the input
    !
    if(present(iBinary))then
      call close_input_file(FFormat%iformat, input_unit, iBinary)
    else
      call close_input_file(FFormat%iformat, input_unit)
    endif
    
    IF (fields_accepted > 0)then
      miniMarketPtr%minimarket_empty = .false.
      ifNewDataConsumed = .true.
    endif  ! fields_accepted > 0

    IF (supermarket_info) then
      call msg('Fields found : ',fields_found)
      call msg('Fields accepted: ',fields_accepted)
      call SYSTEM_CLOCK(after)
      call msg('CPU time used: ', 1e-6*(after - before))
    endif
    
!    call msg('coarse_global_geo_grid')
!    call report(coarse_geo_global_grid)
!    call msg('global_geo_grid')
!    call report(geo_global_grid)
!    call msg('End grid reporting')
!    stop

  END SUBROUTINE store_input_to_supermarket

  !************************************************************************************

  subroutine inputID2MarketID(idIn, id, stack_type, storage_grid, TweakTime, ifStaggerWinds, shpLst, indexVar, time_to_force)
      ! Generates idOut as it will appear in stack after  put_field_to_sm 
      ! The idea is to separate operations with  ID and operations with actual data
      ! and to collect all intelligent parts from put_field_to_sm and
      ! put_field_to_stack into one place
      
      

      implicit none

      ! Imported parameters:
      TYPE(silja_field_id), intent(in) :: idIn
      TYPE(silja_field_id), intent(out) :: id
      integer, INTENT(in) :: stack_type ! dynamic/permanent meteo or dispersion
      type(silja_grid), intent(in) :: storage_grid !(can be adjusted for arakawa if ifAdjust is set)
      integer, intent(in) :: TweakTime
      logical, intent(in) :: ifStaggerWinds
      TYPE(silja_shopping_list), intent(in) :: shpLst
      integer, intent(in) :: indexVar  !Not used by now
      type(silja_time),  intent(in)  :: time_to_force
      character(len=*), parameter :: sub_name = 'inputID2MetMarketID'

      integer ::  iFieldTo, iStack, s,t, mon,year, iQ
      TYPE(silja_grid) :: gridTmp !Tmp variables are used to cover 

      id = idIn  !Local copy to play with

      mon  = fu_mon (fu_start_time(shpLst))
      year = fu_year(fu_start_time(shpLst))


      if(defined(time_to_force))call set_valid_time(id, time_to_force)

      if (stack_type == single_time_stack_flag) then 
        !Some tweak needed
        select case (TweakTime)
          case (static_value,static_climatology) !Valid forever
            call set_valid_time(id, really_far_in_past)
            call set_validity_length(id, very_long_interval)
          case (instant_map)
            call set_validity_length(id, zero_interval)
          case (monthly_climatology)
            !Ignore original month in the ID
            call set_valid_time(id, fu_set_time_utc(year,mon, 1, 0, 0, 0.0))   ! start of the month
            call set_validity_length(id, one_day*fu_days_in_month (mon, year)) ! Whole month
          case (period_valid_map) 
              !Live it as is
          case default
            call set_error("Wrong TweakTime for singletime_stack: "//trim(fu_str(TweakTime)), sub_name)
            return
        end select !TweakTime
      else
        if(fu_list_time_indicator(shpLst) == accept_same_month) then
          call set_valid_time(id, fu_set_time_utc(year, mon, 16, 0, 0, 0.0))  ! middle of the month
        endif
      endif

      if (defined(storage_grid)) then
        if (ifStaggerWinds) then
          iQ = fu_quantity(id)
          if (iQ == u_flag) then
            call set_grid(id, fu_staggered_grid('U_grid', storage_grid, .TRUE., .FALSE.))
          elseif (iq == v_flag) then
            call set_grid(id, fu_staggered_grid('V_grid', storage_grid, .FALSE.,.TRUE.))
          else
            call set_grid(id, storage_grid)
          endif
        else
          call set_grid(id, storage_grid)
        endif

      else
        !Missing grid -- supermarket will decide on final grid
        call set_grid(id, grid_missing)
      endif

  end subroutine inputID2MarketID



  
  !************************************************************************************

  subroutine put_field_to_sm(miniMarketPtr, &     ! Mini market to put in
                           & nStacks, &
                           & idIn, &
                           & grid_data, &
                           & iUpdateType, &       ! Overwrite or not existing fields
                           & stack_type, &        ! stationary or time-dependent
                           & iAccuracy, &         ! of reporjection
                           & stack_indices, &  ! OVERWRITES met_src, stacks are picked directly
                           & Grid_in, &       ! storage grid
                           & ifForceMetSrcAcceptance, TweakTime, ifAdjust, &
                           & ifRandomize, &
                           & shpLstPtr, indexVar, time_to_force)
  ! Unification of a piece replicated many times in store_input_to_supermarket
  ! Could it be also unified with dq_store_2d?
  
  !  !Treat missing value ???
  ! At the unification stage t was enabled only for NetCDF  and only for mapped species
  ! Ignoring it for a while
    implicit none

    ! Imported parameters with intent IN:
    type(mini_market_of_stacks), pointer :: miniMarketPtr
    integer, INTENT(in) :: stack_type ! dynamic/permanent meteo or dispersion
    integer, INTENT(in) :: nStacks !There were some checks above, here using as_is
    integer, intent(in) :: iUpdateType  ! create, overwrite, add, etc - see fields module
    integer, intent(in) :: iAccuracy
    type(silja_grid), intent(in) :: Grid_in
    TYPE(silja_field_id), intent(in) :: idIn
    real, dimension(:), intent(in) :: grid_data
    integer, intent(in) :: TweakTime
    logical, intent(in) :: ifForceMetSrcAcceptance, ifAdjust, ifRandomize
    integer, dimension(:), pointer :: stack_indices ! For boundary market needs,
    TYPE(silja_shopping_list), pointer :: shpLstPtr
    integer, intent(in) :: indexVar
    type(silja_time),  intent(in)  :: time_to_force

    !
    type(silja_grid) :: storage_grid
    type(silam_species) :: species
    TYPE(silja_stack), POINTER :: stack_to_fill
    TYPE(silja_field_id), dimension(:), pointer :: map_targets  
    real, dimension(:), pointer :: map_factors
    type(inVar2modVarsMap), pointer :: ptrMapping
    TYPE(silja_field_id) :: id
    integer ::  iFieldTo, iStack, s,t, mon,year
    logical :: found
    character(len=*), parameter :: sub_name = 'put_field_to_sm'
    
    id = idIn  !Local copy to play with

    mon  = fu_mon(fu_start_time(shpLstPtr))
    year = fu_year(fu_start_time(shpLstPtr))


    if(defined(time_to_force))call set_valid_time(id, time_to_force)

    if (stack_type == single_time_stack_flag) then 
      !Some tweak needed
      select case (TweakTime)
        case (static_value,static_climatology) !Valid forever
          call set_valid_time(id, really_far_in_past)
          call set_validity_length(id, very_long_interval)
        case (instant_map)
          call set_validity_length(id, zero_interval)
        case (monthly_climatology)
          !Ignore original month in the ID
          call set_valid_time(id, fu_set_time_utc(year,mon, 1, 0, 0, 0.0))   ! start of the month
          call set_validity_length(id, one_day*fu_days_in_month (mon, year)) ! Whole month
        case (period_valid_map) 
            !Live it as is
        case default
          call set_error("Wrong TweakTime for singletime_stack: "//trim(fu_str(TweakTime)), sub_name)
          return
      end select !TweakTime
    else
      if(fu_list_time_indicator(shpLstPtr) == accept_same_month) then
        call set_valid_time(id, fu_set_time_utc(year, mon, 16, 0, 0, 0.0))  ! middle of the month
      endif
    endif


    !
    ! Now let's find the proper stack
    !
    do iStack = 1, nStacks
      if(.not. associated(stack_indices))then
        s = fu_met_src_storage_index(miniMarketPtr, stack_type, fu_met_src(id))
      else
        s = stack_indices(iStack)
      endif
      IF (error) then
        call unset_error('store_input_to_supermarket6')
        cycle
      endif
      if(stack_type == single_time_stack_flag) THEN
        stack_to_fill => miniMarketPtr%stacks_singleTime(s)%ptr
      elseif(stack_type == multi_time_stack_flag) THEN
        t = fu_time_storage_index(miniMarketPtr, s, &
                & fu_valid_time(id), ifForceMetSrcAcceptance)
        IF (error) then
          call unset_error('store_input_to_supermarket7')
          cycle
        endif
        stack_to_fill => miniMarketPtr%stacks_multiTime(s, t)%ptr
      else  ! Error - unknown stack type
        call msg('Strange stack type switcher:',stack_type)
        call set_error('Wrong stack type switcher','store_input_to_supermarket')
      endif ! Stack type

      if(defined(Grid_in))then
       storage_grid = Grid_in
      else 
        storage_grid = fu_storage_grid(fu_wdr(stack_to_fill))
      endif

      !
      ! If an input id is mapped to a bunch of target ones, use the mapping
      !
      found = .false.
      if(indexVar /= int_missing)then
        call get_var_mapping(fu_shopping_var(shpLstPtr, indexVar), ptrMapping)
        if(associated(ptrMapping))then
          if(defined(ptrMapping))then
            call get_var_map_targets_and_factors(ptrMapping, map_targets, map_factors)
            do iFieldTo = 1, fu_nrTargets(ptrMapping)
              call set_species(species, &
                             & fu_get_material_ptr(fu_substance_name(map_targets(iFieldTo))), &
                             & fu_mode(map_targets(iFieldTo)))
              call set_species(id, species)
              CALL put_field_to_stack(id, &
                                    & grid_data, &
                                    & stack_to_fill, &
                                    ! No adjustment to storage grid => Arakawa grids are allowed:
                                    & ifAdjust, &  
                                    & iUpdateType, &   ! create, overwrite, add
                                    & ifRandomize, &
                                    & storage_grid, &
                                    & area_missing, &
                                    & factor = map_factors(iFieldTo), &
                                    & ifForceMetSrcAcceptance = ifForceMetSrcAcceptance, &
                                    & iAccuracy=iAccuracy, &
                                    & iOutside = fu_outGridInterpType(fu_quantity(id)))
            enddo
            found = .true.
          endif
        endif
      endif  ! associated ptrMapping
      !
      ! In case of no mapping, just put the field
      !
      if(.not. found)then  ! if no species not found, put a general field

        CALL put_field_to_stack(id, &
                              & grid_data, &
                              & stack_to_fill, &
                              & ifAdjust, &  ! No adjustment to storage grid => Arakawa grids are allowed
                              & iUpdateType, &   ! create, overwrite, add
                              & ifRandomize, &
                              & storage_grid, &
                              & area_missing, &
                              & ifForceMetSrcAcceptance = ifForceMetSrcAcceptance, &
                              & iAccuracy=iAccuracy, &
                              & iOutside = fu_outGridInterpType(fu_quantity(id)))
      endif
      IF (error) THEN
        CALL unset_error('store_input_to_supermarket8')
        call msg("ID in trouble")
        call report(id)
        call msg("")

      END IF

    enddo  ! stacks to store the field to     

  end subroutine put_field_to_sm



  ! ***************************************************************

  SUBROUTINE get_input_field(chFName, FFormat, & ! File where from to take it
                           & idRequested, &       ! The stuff to search
                           & data_requested, &
                           & Grid_in, &          ! storage grid
                           & ifAcceptSameMonth, iOutside, idOut, &
                           & iAccuracy, &
                           & wdr, &
                           & fFillValue_, &
                           & iBinary)              ! for GrADS/NetCDF, if file is already open
    ! 
    ! Finds the fields in the given file. Decodes GRIB.
    ! 
    ! Here it is also ok to use accept_any_met_src, when the file is a "black box".
    !
    ! Data may be stored in its original grid, user-defined grid
    ! or user-defined sub area of the original grid. 
    !
    ! The file may contain fields from several valid times or analyses.
    !
    ! If an unknown field is found, then it is skipped (no error set)
    ! 
    use omp_lib_kinds, only : OMP_lock_kind
    IMPLICIT NONE

    ! Imported parameters
    CHARACTER (LEN=*), INTENT(in) :: chFName
    type(silam_fformat), intent(in) :: FFormat
    TYPE(silja_field_id), INTENT(in), target :: idRequested
    real, dimension(:), intent(out), target :: data_requested
    type(silja_grid), intent(in), target :: Grid_in
    logical, intent(in) :: ifAcceptSameMonth
    type(silja_wdr), intent(in) :: wdr
    TYPE(silja_field_id), intent(out) :: idOut
    integer, intent(in) :: iOutside, iAccuracy
    real, intent(in), optional :: fFillValue_
    integer, intent(inout), optional :: iBinary
    
    ! Local declarations:
    INTEGER :: input_unit, iT, iVar, iLev, fields_found, indFldGrib, indGradsTime, fs_grid
    LOGICAL :: eof, found, selection_done,ifZeroNegatives
    TYPE(silja_field_id) :: id
    type(silja_grid), pointer :: pGrid
    type(silja_grid), target :: gridTmp
    TYPE(silja_field_id), dimension(:), pointer :: idList
    type(silja_time), dimension(:), pointer ::  timeLst
    type(silja_time) ::  timeTmp
    type(silja_interval) :: intervalTmp
    real, dimension(:), pointer :: work_array
    logical :: ifWA !!If release WA

    real :: scale_factor, fFillValue
    integer :: old_unit, iID, nIDs, nTimes
    integer(kind = OMP_lock_kind), dimension(1) :: dummyLocks
    integer :: quant, nx, ny, n_out_of_range, n_patched, n_failed
    character(len=*), parameter :: sub_name = 'get_input_field'


    !
    ! Start from opening file. 
    ! As one file can include several timesteps, for grads and netcdf files try following:
    ! First check, if it's already open. If not, check if exists. If exists, open; 
    ! if not, look for required fields in binary already open.
    !
    call msg('Getting a field from input file:' + chFName)
    found = .false.

    if(present(fFillValue_))then
      fFillValue = fFillValue_
    else
      fFillValue = real_missing
    endif

    old_unit = int_missing
    if(present(iBinary)) old_unit = iBinary

    call open_input_file(chFName, FFormat, input_unit, old_unit, fu_obstime_interval(wdr))

    if(error)return
    if(present(iBinary)) iBinary = input_unit

    !
    ! 2. Go through the file.
    !
    fields_found = 0
    eof = .false.
    call set_missing(idOut)


    !!! By default no remapping
    work_array => data_requested
    ifWA = .False.

    select case (FFormat%iformat)
    case (netcdf_file_flag)
      !
      ! NetCDF has own logic for picking the fields - there is a list of headers available
      !

      if(ifAcceptSameMonth) then 
        ! Take whatever with same month of validity
        call timelst_from_netcdf_file(input_unit, timeLst, nTimes)
        do iT=1,nTimes
          if (fu_mon(timeLst(it)) == fu_mon(fu_valid_time(idRequested))) exit
        enddo
        if (iT > nTimes) then
            call set_error("Requested month not found", sub_name)
            return
        endif
        !Substitute valid time
        id = idRequested
        call set_valid_time(id, timeLst(it))
        call field_indices_from_netcdf_file(input_unit, id, idOut, iVar, it, ilev)

      else
          call field_indices_from_netcdf_file(input_unit, idRequested, idOut, iVar, it, ilev)
      endif

          
      IF (fu_field_id_covers_request(idOut, idRequested, .true.)) THEN
 !           call msg("Got it!")

            if(defined(Grid_in) .and.( .not. Grid_in == fu_grid(idOut)) )then
              work_array => fu_work_array(fu_number_of_gridpoints(fu_grid(idOut)))
              ifWA = .True.
              if(error)return
            endif

            call read_field_from_netcdf_indices(input_unit,  ivar, it, ilev, &
                                           & work_array, fFillValue) ! the grid-data
            found = .true.

      endif    ! field id covers the request

    case  (grads_file_flag)
      !
      ! get the correct grads time index. Note that direction is unknown here
      !
      timeTmp = fu_valid_time(idRequested)
      if(defined(timeTmp))then
        indGradsTime = fu_grads_time_index(input_unit, timeTmp, back_and_forwards, ifAcceptSameMonth)
      else
        indGradsTime = 1  ! can be id_simple, then any time is OK
      endif

      var_cycle: do iVar = 1, fu_n_gvars(input_unit)
        do iLev = 1, fu_n_gVar_levs(input_unit, iVar)

          call get_grads_var_metadata(input_unit, iVar, iLev, indGradsTime, idOut)
          if(error)return

!call report(idOut)
          
          fields_found = fields_found + 1


          if(fu_field_id_covers_request(idOut, idRequested, .true.))then

            if(defined(Grid_in) .and. (.not.(fu_grid(idOut) == Grid_in)))then
              work_array => fu_work_array(fu_number_of_gridpoints(fu_grid(idOut)))
              if(error)return
              ifWA = .True.
            endif

            call read_field_from_grads_indices(input_unit, &
                                     & iVar, iLev, indGradsTime, & !idOut, &   !Requested,&
                                     & work_array, fFillValue) !data_requested) ! the grid-data
            found = .true.
            exit var_cycle  ! the field has been found

          endif  ! new ID OK
        enddo   ! iLev
      enddo  var_cycle

    case (test_field_value_flag)
      !
      ! For the test field, the grid must be selected to be something. 
      !
      work_array => data_requested
      if(defined(Grid_in))then 
        pGrid => Grid_in
      else
        if(defined(dispersion_grid))then
          pGrid => dispersion_grid
        else
          gridTmp = geo_global_grid
          pGrid => gridTmp
        endif
      endif
      if(error)return
      
      timeTmp = fu_valid_time(idRequested)
      if (.not. defined(timeTmp)) timeTmp = really_far_in_past
      
      intervalTmp = fu_forecast_length(idRequested)
      if (.not. defined(intervalTmp)) intervalTmp = zero_interval

      call make_test_field(chFName, &  ! chFName contain the field name and value
                               & idOut, &    ! idRequested
                         & work_array, &
                         & timeTmp, &     ! valid time
                         & intervalTmp, & ! froecast lengthh
                         & pGrid)
call msg('The following test field has been created (get_input_field)')
call report(idOut)
call msg('')
!        call ooops("test field")

            found = .true.  ! this field is OK by-definition
            eof = .true.

    case ( grib_file_flag)

       call id_list_from_grib_file(input_unit, idList, nIDs)
       if (fu_fails(.not. error, "after id_list_from_grib_file",'get_input_field')) return

       do iID = 1,nIDs
         if(defined(idList(iID)))then
            IF (fu_field_id_covers_request(idList(iID), idRequested, .true.)) exit
         endif  ! read a resonable id
       ENDDO

       if( iID <= nIDs )then !Loop above terminated

         if(defined(Grid_in))then
           fs_grid=fu_number_of_gridpoints(fu_grid(idList(iID)))
           work_array => fu_work_array(fs_grid)
           if(error)return
           ifWA = .True.
         endif
          
         call read_field_from_grib_file(input_unit, idList(iID), iID, work_array)
         !!!! FIXME No forced missing_replacement in read_field_from_grib_file, 
         ! It is decided from quantity
         if (fu_fails(.not. error, "after read_field_from_grib_file",'get_input_field')) return
         found = .true.
       endif

    case (ascii_file_flag)
      !
      ! Neither of the formats with known systematic list of fields (netCDF, GrADS, test field)
      ! Have to read all through looking for the right field
      !
      DO while(.not.eof)

!      print *, '* Next field...', fields_found

            if(defined(Grid_in))then
              work_array => fu_work_array()
              ifWA=.true.
            endif
            CALL read_next_field_from_ascii_file(input_unit, &
                                               & eof,&
                                               & idOut,&
!                                               & fMissingValue, & ! the grid-data
                                               & work_array, fFillValue) !data_requested)
 
        !
        ! An error may have happened
        !
        IF (error) THEN
          CALL msg_warning('error while getting new field, but continuing.','get_input_field')
          CALL unset_error('get_input_field4')
          IF (eof) THEN
            EXIT
          ELSE
            CYCLE
          END IF
        END IF
        !
        ! Check if the consumed field is what we are after
        !
        if(defined(idOut))then
          fields_found = fields_found + 1
          IF (fu_field_id_covers_request(idOut, idRequested, .true.)) THEN
            found = .true.
            exit  ! main do loop
          END if  ! id covers the request
        endif  ! read a resonable id

        IF (error) call unset_error('get_input_field5')
  
      END DO ! loop over inputs
    case default
            call msg('Non-supported file type',FFormat%iformat)
            call set_error('Non-supported file type','get_input_field')

    end select  ! input file type
    !
    ! If the field has been found, it may need to be regridded
    !
    if(found)then
      if(ifWA)then 

          quant = fu_quantity(idOut)
          call  remap_field( fu_grid(idOut), &
                        & work_array, &
                        & Grid_in, &
                        & data_requested, &
                        !!& ifRandomise, iAccuracy, method, iOutside, omp_locks)
                        & fu_if_randomise(wdr), & ! random smoothing of reprojection
                        & iAccuracy, &
                        & fu_regridding_method(quant), &
                        & iOutside, dummyLocks)

         
           call grid_dimensions(Grid_in, nx, ny)
           call check_quantity_range(quant, data_requested, nx*ny, nx, &
              & .false., .False., &  ! ifRequireValidity, ifSilent 
             & n_out_of_range, n_patched, n_failed)

           if(error) then
            call set_error("error after check_quantity_range ", 'get_input_field')
            return
           endif
!          call grid_data_horizontal_select(fu_grid(idOut), &     ! grid_original
!                                         & work_array, &      ! grid_data
!                                         & selection_done, &  ! selection_done
!                                         & gridTmp,&          ! grid_new
!                                         & data_requested, &  ! selected_grid_data
!                                         & fu_if_randomise(wdr), & ! random smoothing of reprojection
!                                         & iAccuracy, &
!                                         & fu_regridding_method(fu_quantity(id)), &
!                                         & fFillValue, &          ! fMissingValue
!                                         & iOutside)

          call set_grid(idOut,Grid_in)
          call free_work_array(work_array) ! can be a huge array, so release it here
        if(error)return

      endif  ! specific grid requested


    else
      !
      ! The field is not needed. Set ID missing and release the buffer
      !
      call msg('No fields found, fields checked:', fields_found)
#ifdef DEBUG_SRC
   call msg('ID found in the file:')
   call report(idOut)
#endif
      call set_missing(idOut)
    endif  ! found
    
!    idOut = id
    !
    ! All done, close the input
    !
    if(present(iBinary))then
      call close_input_file(FFormat%iformat, input_unit, iBinary)
    else
      call close_input_file(FFormat%iformat, input_unit)
    endif

  END SUBROUTINE get_input_field


  !******************************************************************

  subroutine open_input_file(chFName, FFormat, input_unit, old_unit, obstime_interval_)
    !
    ! Opens the input file taking into account its type. For NetCDF, GRIB and grads, old_unit
    ! can be provided, then no opening happens, just checking that the file name matches one
    ! for already open file
    !
    implicit none

    ! Imported parameters
    character(len=*), intent(in) :: chFName
    type(silam_fformat), intent(in) :: FFormat
    integer, intent(out) :: input_unit
    integer, intent(in)  :: old_unit
    type(silja_interval), optional, intent(in) :: obstime_interval_ !Used only to fix times in GRIB
    type(silja_interval) :: obstime_interval !Used only to fix times in GRIB

    obstime_interval = interval_missing
    if (present(obstime_interval_)) obstime_interval = obstime_interval_




    select case(FFormat%iformat)
      case(grib_file_flag)
        if(old_unit /= int_missing)then
          if(fu_grib_filename(old_unit) ==chFName)then
              input_unit = old_unit
              return
            else
              call msg('Closing grib file i:' + fu_grib_filename(old_unit))
              call close_grib_file_i(old_unit)
          endif
        endif
        CALL open_grib_file_i(chFName, input_unit, obstime_interval)

      case(ascii_file_flag)
        call msg('Opening ASCII file:' + chFName)
        call open_ascii_file_i(chFName, input_unit)

      case(netcdf_file_flag)
        if(old_unit /= int_missing)then
          if(fu_netcdf_filename(old_unit) == chFName)then
              input_unit = old_unit
              return
            else
              call msg('Closing NETCDF file:' + fu_netcdf_filename(old_unit))
              call close_netcdf_file(old_unit)
          endif
        endif
        if (error) return
        call msg('Opening NETCDF file:' + chFName)
        input_unit = open_netcdf_file_i(chFName, FFormat)

      case(grads_file_flag)
!        call msg('Opening GrADS file:' + chFName)
          if(old_unit /= int_missing)then
            if(fu_grads_sctl_filename(old_unit)==chFName)then ! Has to be the super-ctl file name
              input_unit = old_unit
              return
            else
              call close_gradsfile_i(old_unit)
            endif
          endif
          input_unit =  fu_open_gradsfile_i(chFName)

      case(test_field_value_flag)        ! Do nothing - no file
        call msg('Making test field (open_input_file):' + chFName)

      case default
        call msg('Non-supported file type',FFormat%iformat)
        call set_error('Non-supported file type','open_input_file')
    end select
  end subroutine open_input_file


  !****************************************************************************

  subroutine close_input_file(iFFormat, input_unit, iBinary)
    !
    ! Just closes the input file taking into account its format
    !
    implicit none

    integer, intent(in) :: iFFormat, input_unit
    integer, intent(in), optional :: iBinary

    select case(iFFormat)
      case(grib_file_flag)
        CALL close_grib_file_i(input_unit)
      case(ascii_file_flag)
        call close_ascii_file_i(input_unit)
      case(netcdf_file_flag)
        if(.not. present(iBinary))then
          call close_netcdf_file(input_unit)
        endif
      case(grads_file_flag)
        if(.not. present(iBinary))then
          call close_gradsfile_i(input_unit)
        endif
      case(test_field_value_flag)        ! Do nothing - no file
      case default
        call msg('Non-supported file type',iFFormat)
        call set_error('Non-supported file type','close_input_file')
    end select

  end subroutine close_input_file


  !***************************************************************************
  
  subroutine get_id_list_from_input_file(input_unit, fileFormat, idList, number_of_ids)
    !
    ! At least Grads and Netcdf allow for generation of a full ID list for the file. This function
    ! does just that by calling the corresponding subroutines
    !
    implicit none
    
    integer, intent(in) :: input_unit, fileFormat
    type(silja_field_id), dimension(:), pointer :: idList
    integer, intent(out) :: number_of_ids

    select case(fileFormat)
      case(grads_file_flag)
        call get_grads_IDs(input_unit, idList, number_of_ids)
        
      case(netcdf_file_flag)
        call id_list_from_netcdf_file(input_unit, idList, number_of_ids) !ifFirstTime, ifFirstLev)
        
      case default
        call set_error('Only grads and netcdf files support field ID listing','get_id_list_from_input_file')
        nullify(idList)
        number_of_ids = 0
    end select
      
  end subroutine get_id_list_from_input_file 
  
  

  !*********************************************************************
  !
  ! Writing down a stack to the bunch of output files
  !
  ! ********************************************************************

  !********************************************************************************

  subroutine write_stack_to_files(stack, iSource, grads_funit, netcdf_funit, grib_funit, iGribType, &
                                & now, if_regular_output_times)
    !
    ! Drops the whole stack to netcdf file. 
    !
    IMPLICIT NONE

    ! Imported parameters with intent IN:
    TYPE(silja_stack), TARGET, INTENT(in):: stack
    integer, intent(in) :: iSource, grads_funit, netcdf_funit, grib_funit, iGribType
    type(silja_time), intent(in) :: now
    logical, intent(in) :: if_regular_output_times

    ! Local declarations:
    TYPE(silja_stack), POINTER :: stackPtr
    TYPE(silja_field), POINTER :: field_2d
    type(silja_field_id) :: idTmp
    integer :: i
    real, dimension(:), pointer :: pTmp
    !
    ! Checking
    !
    stackPtr => stack
    IF(.not.ASSOCIATED(stackPtr))THEN
      CALL set_error('Non-associated output stack','write_stack_to_files')
      RETURN
    END IF
    IF(.not.defined(stackPtr))THEN
      CALL set_error('Undefined output stack','write_stack_to_files')
      RETURN
    END IF
    !
    ! Simple 2d fields (if any) - the full content of the stack
    !
!    call msg("Writing stack:"+fu_name(stackPtr))
!    call report(stackPtr)

    DO i = 1, fu_number_of_2d_fields(stackPtr)
      field_2d => fu_get_2d_field(stackPtr,i)
      IF(ASSOCIATED(field_2d))THEN
        IF(defined(field_2d)) THEN
          idTmp = fu_id(field_2d)
          if(fu_if_internal_silam_field(idTmp))cycle    ! mds == silam_internal_src
          !
          ! Fix valid time for output: whatever the time is in the stack, we want to output it
          ! here and now. Needed for shifted-meteo runs and for various mid-time-step valid 
          ! and as_is fields
          !
          call set_valid_time(idTmp, now)
          !
          ! GRADS
          
if(.not. defined(fu_level(idTmp)))then
  call report(idTmp)
  call msg('Undefined level')
endif
          if(grads_funit /= int_missing)then
            CALL write_next_field_to_gradsfile(grads_funit, idTmp, fu_grid_data(field_2d), now, &
                                             & if_regular_output_times)
          endif
          !
          ! GRIB
          if(grib_funit /= int_missing)then
            CALL write_next_field_to_gribfile(iGribType, grib_funit, idTmp, fu_grid_data(field_2d))
          endif
          !
          ! NetCDF
          if(netcdf_funit /= int_missing) then
            pTmp => fu_grid_data(field_2d)
            CALL write_next_field_to_netcdf_file(netcdf_funit, idTmp, pTmp)
          endif
          !
          ! went through?
          !
          if(error) then 
            call msg_warning("Failed to write field","write_stack_to_files")
            call report(idTmp)
            call report(stackPtr)
            return
          END IF
        END IF

      END IF
    END DO

  END SUBROUTINE write_stack_to_files



!  !****************************************************************
!
!  subroutine merge_data_to_supermarket(miniMarketPtr, idIn, dataIn, targetId, stack_type)
!    !
!    ! Merges the given field to the supermarket stacks. Uses the merge_data_to_stack
!    ! subroutine from the stacks module. In fact, hardly plays any role
!    ! except for the encapsulation of that function and the stackvector
!    !
!    implicit none
!
!    ! Imported parameters 
!    type(mini_market_of_stacks), intent(inout) :: miniMarketPtr
!    type(silja_field_id), intent(in) :: idIn, targetId
!    real, dimension(:), pointer :: dataIn
!    integer, intent(in) :: stack_type
!
!    ! Local variables
!    type(silja_stack), pointer :: stackPtr
!    integer :: s, t
!    !
!    ! So far, if the field is not in permanent stack, we don't know how to find it
!    !
!
!    s = fu_met_src_storage_index(miniMarketPtr, stack_type, fu_met_src(idIn))
!    if(stack_type == single_time_stack_flag)then
!        stackPtr => miniMarketPtr%stacks_singleTime(s)%ptr
!    elseif(stack_type == multi_time_stack_flag)then
!      t = fu_time_storage_index(miniMarketPtr, s, fu_valid_time(idIn), .false.)
!      stackPtr => miniMarketPtr%stacks_multiTime(s,t)%ptr
!    else
!      call msg('Strange stack type', stack_type)
!      call set_error('Strange stack type', 'merge_data_to_supermarket')
!    endif
!
!    call merge_data_to_stack(idIn, dataIn, stackPtr, targetId)
!
!  end subroutine merge_data_to_supermarket



  ! ***************************************************************

  ! ***************************************************************
  !
  !
  !        GET DATA OUT OF THE SUPERMARKET
  !
  !
  ! ***************************************************************
  ! ***************************************************************

  function fu_get_field_from_mm_general(miniMarket, id, ifMandatory) result(field)
    !
    ! A generic function searching the given mini-market for the given ID.
    ! Some elements can be missing, then the first-met match for the defined parameters is returned
    !
    implicit none
    
    ! Return value of this function:
    TYPE(silja_field), POINTER :: field

    ! Imported parameters with intent IN or pointer:
    type(mini_market_of_stacks), intent(in) :: miniMarket
    type(silja_field_id), intent(inout) :: id
    logical, intent(in) :: ifMandatory     ! whether the set error if not found
    
    ! Local variables
    TYPE(silja_stack), POINTER :: stack
    logical :: found
    integer :: iSrc, iTime
    
    !
    ! Somethig simple...
    !
    NULLIFY(field)
    IF (miniMarket%minimarket_empty)RETURN

    !
    ! If time is defined, look in the multi-time stack
    !
    if(defined(fu_valid_time(id)))then
      IF (.NOT. miniMarket%stack_multiTime_initialized) THEN
        if(ifMandatory) CALL set_error('minimarket not init','fu_get_field_from_mm_general')
        RETURN
      END IF

      ! Where the data shoud be:
      stack => fu_closest_sm_met_src_time(miniMarket, &
                                        & fu_met_src(id), fu_valid_time(id), back_and_forwards, &
                                        & ifMandatory)
      if(error)then
        call unset_error('fu_get_field_from_mm_general1')
        return
      endif

      if(associated(stack))then
        call set_valid_time(id,fu_valid_time(stack))
        call find_field_from_stack(stack, id, field, found)  ! got it?
        if(found)return
        if(error)call unset_error('fu_get_field_from_mm_general22')
      endif

    endif  ! if time is defined
    !
    ! The id valid time is undefined or the field is not found, do the complete search through 
    ! the mini-market
    !
    if(miniMarket%stack_singleTime_initialized)then
      !
      ! Check whether the stuff can be found in permanent stack
      !
      do iSrc = 1, size(miniMarket%stacks_singleTime)
        stack => miniMarket%stacks_singleTime(iSrc)%ptr
        call find_field_from_stack(stack, id, field, found)  ! got it?
        if(found)return
        if(error)call unset_error('fu_get_field_from_mm_general3')
      end do
    endif
      
    if(miniMarket%stack_multiTime_initialized)then
      !
      ! Check whether the stuff can be found in time-related stacks
      !
      do iTime = 1, size(miniMarket%stacks_multiTime,2)
        do iSrc = 1, size(miniMarket%stacks_multiTime,1)
          stack => miniMarket%stacks_multiTime(iSrc, iTime)%ptr
          call find_field_from_stack(stack, id, field, found)  ! got it?
          if(found)return
          if(error)call unset_error('fu_get_field_from_mm_general4')
        end do
      end do
    endif

    nullify(field)  ! No luck, id is not found

  end function fu_get_field_from_mm_general


  !******************************************************************************************

  FUNCTION fu_sm_obstime_field(miniMarketPtr, met_src,quantity,level,time,direction) result(field)
    !
    ! Finds nwpm 2D-fields from memory. Data is returned as a
    ! pointer to a silja_field.
    !
    ! No vertical interpolation is done, so the
    ! quantity has to be available on the given level. If the quantity
    ! defines it's own level (like 2m temperature or snowcover etc.)
    ! the level can be undefined.
    !
    ! If data from any met_src are acceptable - keep met_src as met_src_missing.
    ! Then the first met_src provided the fit to time and quantity 
    ! will be used.
    !
    ! No time interpolation is done, so data is returned to the
    ! closest existing sm observation time.
    ! Direction may have values, backwards, forwards,
    ! back_and_forwards or single_time.
    ! If data is searched
    ! only for one, exact time (not searched in any direction) then
    ! set direction = single_time.
    !
    ! In the field is not found, an undefined field (nullified pointer)
    ! is returned.
    ! 
    IMPLICIT NONE

    ! Return value of this function:
    TYPE(silja_field), POINTER :: field

    ! Imported parameters with intent IN or pointer:
    type(mini_market_of_stacks) :: miniMarketPtr
    type(meteo_data_source), INTENT(in) :: met_src
    INTEGER, INTENT(in) :: quantity
    TYPE(silja_level), INTENT(in)  :: level
    TYPE(silja_time), INTENT(in) :: time
    INTEGER, INTENT(in) :: direction ! in time

    ! Local declarations:
    TYPE(silja_field_id) :: id
    TYPE(silja_time) :: obstime
    LOGICAL :: found_in_memory, found_in_hit_list
    REAL :: before
    TYPE(silja_stack), POINTER :: stack

    before = fu_total_cpu_time_sec()

    ! ------------------------------------------------------
    !
    ! 1. Find correct times and complete id.
    !
    IF (.NOT. miniMarketPtr%stack_multiTime_initialized) THEN
      CALL set_error('minimarket not init','fu_sm_obstime_field')
      RETURN
    END IF

    IF (miniMarketPtr%minimarket_empty) THEN
      NULLIFY(field)
      RETURN
    END IF

    ! Where the data shoud be:
    stack => fu_closest_sm_met_src_time(miniMarketPtr, met_src, time, direction)
    IF (error) RETURN

    obstime = fu_valid_time(stack)
    IF (error) RETURN

    id = fu_set_field_id_simple(met_src, quantity, obstime, level)
    IF (error) RETURN


    ! ------------------------------------------------------
    !
    ! 2. Check the hit list.
    !
    CALL find_field_from_hit_list(miniMarketPtr, id, field, found_in_hit_list)

    IF (found_in_hit_list) THEN
      fields_delivered = fields_delivered + 1
      hit_count = hit_count + 1
      cpu_usage = cpu_usage + fu_total_cpu_time_sec() - before
      RETURN
    ELSE
      hit_miss = hit_miss + 1
    END IF

    ! ----------------------------------------------------
    !
    ! 3. Check the correct time's stack for desired field.
    !    Note that searching the field not via just-set id but via details
    !    is probably faster and also allows tricks with the non-zero validity length
    !
    call find_field_from_stack(met_src, &
                             & quantity,&
                             & obstime,&   ! valid time
                             & stack,&
                             & field,&
                             & found_in_memory) !, &
                             !& chSubstNm, &  ! optional
                             !& fModeVal, &   ! optional
                             ! fWaveLen))     ! optional
!    CALL find_field_from_stack(stack, id, field, found_in_memory)
    IF (error) RETURN

    IF (found_in_memory) THEN
      CALL put_field_to_hit_list(miniMarketPtr, field)
      fields_delivered = fields_delivered + 1
    ELSE
      call msg( 'supermarket contents:')
      CALL report(miniMarketPtr)
      CALL report(id)
      CALL set_error('field not found in supermarket', 'fu_sm_obstime_field')
      NULLIFY(field)
!!!$      PRINT *, '  '
!!!$      PRINT *, '  THIS FIELD IS  being nullified !!!  '
!!!$      CALL report(id)
!!!$      PRINT *, '  '
    END IF
!!!$      PRINT *, 'supermarket contents:'
!!!$      CALL print_supermarket_contents()
    cpu_usage = cpu_usage + fu_total_cpu_time_sec() - before

  END FUNCTION fu_sm_obstime_field


  ! ***************************************************************



  FUNCTION fu_sm_obstime_3d_field(miniMarketPtr, met_src,&
                                & quantity,&
                                & time,&
                                & direction) result(field_3d)

    ! Description:
    ! Finds a 3D scalar field for closest observation time.
    !
    ! If 3d-field contains hybrid-level data, then it always also
    ! contains corresponding surface pressure field (which is needed
    ! in vertical administration).
    !
    ! Direction may have values, backwards, forwards, back_and_forwards or single_time.
    ! If data is searched only for one, exact time (not searched in any direction) then
    ! set direction = single_time.
    ! All units: SI
    !
    IMPLICIT NONE
    !
    ! Return value of this function:
    TYPE(silja_3d_field), POINTER :: field_3d

    ! Imported parameters with intent(in):
    type(mini_market_of_stacks), intent(in) :: miniMarketPtr
    type(meteo_data_source), INTENT(in) :: met_src
    INTEGER, INTENT(in) :: quantity
    TYPE(silja_time), INTENT(in) :: time
    INTEGER, INTENT(in) :: direction

    ! Local declarations:
    TYPE(silja_stack), POINTER :: stack
    TYPE(silja_time) :: obstime
    REAL :: before
    LOGICAL :: found_in_stack

    before = fu_total_cpu_time_sec()

    IF (.NOT. miniMarketPtr%stack_multiTime_initialized) THEN
      CALL set_error('minimarket not init','fu_sm_obstime_3d_field')
      RETURN
    END IF


    IF (miniMarketPtr%minimarket_empty) THEN
      NULLIFY(field_3d)
      RETURN
    END IF

    ! Where the data should be:
    stack => fu_closest_sm_met_src_time(miniMarketPtr, met_src, time, direction)
    IF (error) RETURN

    obstime = fu_valid_time(stack)
    IF (error) RETURN

    CALL find_field_3d_from_stack(met_src, quantity, obstime, stack, field_3d, found_in_stack)

    IF (found_in_stack) THEN
      fields_3d_delivered = fields_3d_delivered + 1
    ELSE
      !      CALL print_supermarket_contents()
      CALL set_error('field 3d not found in supermarket', 'fu_sm_obstime_3d_field')
      call msg('Meteo src:' + fu_name(met_src) + ',' + fu_quantity_string(quantity) + ',' + &
             & fu_str(obstime))
      NULLIFY(field_3d)
    END IF

    cpu_usage = cpu_usage + fu_total_cpu_time_sec() - before

  END FUNCTION fu_sm_obstime_3d_field


  ! ***************************************************************


  FUNCTION fu_sm_obstime_3d_windfield(miniMarketPtr, met_src, time, direction) result(field_3d)

    ! Description:
    ! Finds a 3D windfield for closest observation time.
    !
    ! If 3d-field contains hybrid-level data, then it always also
    ! contains corresponding surface pressure field (which is needed
    ! in vertical administration).
    !
    ! Direction may have values, backwards, forwards,
    ! back_and_forwards or single_time.
    ! If data is searched
    ! only for one, exact time (not searched in any direction) then
    ! set direction = single_time.
    ! All units: SI
    !
    IMPLICIT NONE
    !
    ! Return value of this function:
    TYPE(silja_3d_windfield), POINTER :: field_3d

    ! Imported parameters with intent(in):
    type(mini_market_of_stacks), intent(in) :: miniMarketPtr
    type(meteo_data_source), INTENT(in) :: met_src
    TYPE(silja_time), INTENT(in) :: time
    INTEGER, INTENT(in) :: direction

    ! Local declarations:
    TYPE(silja_stack), POINTER :: stack
    TYPE(silja_time) :: obstime
    REAL :: before
    LOGICAL :: found_in_stack

    before = fu_total_cpu_time_sec()

    IF (.NOT. miniMarketPtr%stack_multiTime_initialized) THEN
      CALL set_error('minimarket not init','fu_sm_obstime_3d_windfield')
      RETURN
    END IF


    IF (miniMarketPtr%minimarket_empty) THEN
      NULLIFY(field_3d)
      RETURN
    END IF

    ! Where the data shoud be:
    stack => fu_closest_sm_met_src_time(miniMarketPtr, met_src, time, direction)
    IF (error) RETURN

    obstime = fu_valid_time(stack)
    IF (error) RETURN

    CALL find_wind_3d_from_stack(met_src, obstime, stack, field_3d, found_in_stack)

    IF (found_in_stack) THEN
      windfields_3d_delivered = windfields_3d_delivered + 1
    ELSE
      call msg('Met source:' + fu_name(met_src) + ',' + fu_str(obstime))
      CALL set_error('wind 3d not found in supermarket','fu_sm_obstime_3d_windfield')
      NULLIFY(field_3d)
    END IF

    cpu_usage = cpu_usage + fu_total_cpu_time_sec() - before

  END FUNCTION fu_sm_obstime_3d_windfield



  ! ***************************************************************


  FUNCTION fu_sm_simple_3d_field(miniMarketPtr, met_src, quantity, stack_type) result(field_3d)

    ! Description:
    ! Finds permanent (time-independent) 2D-fields from memory.
    ! Data is returned as a pointer to a silja_field.
    !
    ! Fields are searched from a separate stack, which congtains
    ! only those fields that were put there with permanent-flag.
    !
    ! In the field is not found, an undefined field (nullified pointer)
    ! is returned.
    ! 
    ! This used to need a level. Don't know why.
    !
    IMPLICIT NONE

    ! Return value of this function:
    TYPE(silja_3d_field), POINTER :: field_3d

    ! Imported parameters with intent IN:
    type(mini_market_of_stacks), intent(in) :: miniMarketPtr
    type(meteo_data_source), INTENT(in) :: met_src
    INTEGER, INTENT(in) :: quantity
    integer, intent(in) :: stack_type

    ! Local declarations:
    TYPE(silja_time) :: obstime
    INTEGER :: ind
    LOGICAL :: found_in_memory, found_in_hit_list
    REAL :: before
    TYPE(silja_stack), POINTER :: stack

    nullify(field_3d)


    IF (.NOT. miniMarketPtr%stack_singleTime_initialized) THEN
      CALL set_error('minimarket not init','fu_sm_simple_3d_field')
      RETURN
    END IF

    IF (miniMarketPtr%minimarket_empty) THEN
      RETURN
    END IF
 
    IF (error) RETURN

    stack => miniMarketPtr%stacks_singleTime(fu_met_src_storage_index(miniMarketPtr, stack_type, met_src))%ptr
    CALL find_field_3d_from_stack(met_src, &
                                & quantity, &
                                & time_missing, &
                                & stack, &
                                & field_3d, &
                                & found_in_memory)
    IF (error) RETURN

    IF (.NOT.found_in_memory) THEN
      CALL set_error('field not found in supermarket','fu_sm_simple_3d_field')
      !CALL report(id)
    END IF

  END FUNCTION fu_sm_simple_3d_field



  ! ***************************************************************


  FUNCTION fu_sm_simple_field(miniMarketPtr, met_src, quantity, level, stack_type) result(field)
    ! 
    ! Finds permanent (time-independent) or dispersion 2D-fields from memory.
    ! Data is returned as a pointer to a silja_field.
    !
    ! In the field is not found, an undefined field (nullified pointer)
    ! is returned.
    ! 
    IMPLICIT NONE

    ! Return value of this function:
    TYPE(silja_field), POINTER :: field

    ! Imported parameters with intent IN:
    type(mini_market_of_stacks), intent(in) :: miniMarketPtr
    type(meteo_data_source), INTENT(in) :: met_src
    INTEGER, INTENT(in) :: quantity
    TYPE(silja_level), INTENT(in) :: level
    integer, intent(in) :: stack_type

    ! Local declarations:
    TYPE(silja_field_id) :: id
    TYPE(silja_time) :: obstime
    INTEGER :: ind
    LOGICAL :: found_in_memory, found_in_hit_list
    REAL :: before
    TYPE(silja_stack), POINTER :: stack

    nullify(field)

    IF (.NOT. (miniMarketPtr%stack_singleTime_initialized)) THEN
      CALL set_error('minimarket not init','fu_sm_simple_field')
      RETURN
    END IF

    IF (miniMarketPtr%minimarket_empty) THEN
      RETURN
    END IF

    id = fu_set_field_id_simple(met_src, quantity,  time_missing, level)
    IF (error) RETURN

    stack => miniMarketPtr%stacks_singleTime(fu_met_src_storage_index(miniMarketPtr, stack_type, met_src))%ptr
    CALL find_field_from_stack(stack, id, field, found_in_memory)
    IF (error) RETURN


    IF (.NOT.found_in_memory) THEN
      CALL set_error('field not found in supermarket','fu_sm_simple_field')
      CALL report(id)
    END IF

  END FUNCTION fu_sm_simple_field


  ! ***************************************************************

  ! ***************************************************************
  !
  !
  !   HERE BEGIN DERIVED QUANTITY-TOOLS
  !
  ! Used for creating new, derived 3D- or 2D fields from data
  ! that is in supermarket. The beginnig "dq" refers to these
  ! fields that are not imported to, but created inside supermarket.
  !
  ! ***************************************************************
  ! ***************************************************************

  
  subroutine find_field_storage_2d(minimarketptr, idin, stack_type, pfield)
    !
    ! Finds a place holder for one time-dependent 2d-field in the supermarket data pool
    !
    ! Obs. the horizontal grid/area of the field stored in memory may NOT change. 
    ! If supermarket storage_grid or storage-area is set, they overrule the grid given by field_id.
    ! This causes inconsistency, which will result in failed search
    !
    IMPLICIT NONE

    ! imported parameters
    type(mini_market_of_stacks), intent(inout) :: miniMarketPtr
    TYPE(silja_field_id), INTENT(in) :: IDin
    integer, INTENT(in) :: stack_type
    type(silja_field), pointer :: pfield !output: pointer to the stored field

    ! Local variables
    TYPE(silja_stack), POINTER :: pStack
    integer :: s, t
    TYPE(silja_field_id) :: ID
      character(len=*), parameter :: sub_name = 'find_field_storage_2d'
    !
    ! Check input parameters
    !

    IF(fu_fails(defined(IDin),'undefined id given',sub_name))RETURN

    ! Find the stack we need to address
    !
      s = fu_met_src_storage_index(miniMarketPtr, stack_type, fu_met_src(IDin))
      if(fu_fails(s >= 1 .and. s <= size(miniMarketPtr%stacks_multiTime,1), &
                & 'Bad met_src, can not find proper stack',sub_name))return
      select case(stack_type)
        case(single_time_stack_flag)
          pStack => miniMarketPtr%stacks_singleTime(s)%ptr
        case(multi_time_stack_flag)
          t = fu_time_storage_index(miniMarketPtr, s, fu_valid_time(IDin), .false.)
          if(fu_fails(t >= 1 .and. t <= size(miniMarketPtr%stacks_multiTime,2), &
                    & 'Bad met_src, can not find proper stack',sub_name))return
          pStack => miniMarketPtr%stacks_multiTime(s, t)%ptr
        case default
          call msg('Strange stack type:', stack_type)
          call set_error('Strange stack type:', sub_name)
      end select

    id = IDin
    if (.not. defined(fu_grid(id))) then
        call set_grid(id, fu_storage_grid(fu_wdr(pStack)))
    endif

    !
    ! Having the stack determined, get the field there - or create it.
    !
    CALL get_field_place_in_stack(pStack, id, pField, .false.) ! ifForceMDS=.false.
    if (error) then
      call set_error('Failed to create the place for this field',sub_name)
      call report(IDin)
      return
    endif

    minimarketptr%minimarket_empty = .false.

  end subroutine find_field_storage_2d
  !**********************************************************************
  
  subroutine find_field_data_storage_2d(minimarketptr, idin, stack_type, pdata, pid)

    !    wrapper around find_field_storage_2d to get data and id pointers directly
    implicit none

    ! imported parameters
    type(mini_market_of_stacks), intent(inout) :: minimarketptr
    type(silja_field_id), intent(in) :: idin
    integer, intent(in) :: stack_type
    real, dimension(:), pointer :: pdata        ! output: pointer to the storage place
    type(silja_field_id), pointer, optional :: pid ! output: poiner to actual idin 

    type(silja_field), pointer :: pfield

    if (present(pid)) pid => null()
    call find_field_storage_2d(minimarketptr, idin, stack_type, pfield)
    pdata => fu_grid_data(pfield)
    if (present(pid)) pid => fu_id(pfield)

  end subroutine find_field_data_storage_2d


  ! ***************************************************************


  SUBROUTINE dq_store_2d(miniMarketPtr, id, grid_data, stack_type, ifRandomise, &
                       & stackptr, iUpdateType, storage_grid, fMissingValue_, iAccuracy_)

    ! Description:
    ! Stores one time-dependent 2d-field into supermarket data pool.
    !
    ! Obs. the horizontal grid/area of the field stored in memory
    ! might actually change. If supermarket
    ! storage_grid is set, they overrule the grid
    ! given by field_id.
    !
    IMPLICIT NONE
    type(mini_market_of_stacks), intent(inout) :: miniMarketPtr
    TYPE(silja_field_id), INTENT(in) :: id
    REAL, DIMENSION(:), INTENT(in) :: grid_data
    integer, INTENT(in) :: stack_type
    logical, intent(in) :: ifRandomise
    integer, intent(in), optional :: iUpdateType
    TYPE(silja_stack), POINTER, optional :: stackptr
    TYPE(silja_grid), INTENT(in), optional :: storage_grid
    real, intent(in), optional :: fMissingValue_
    integer, intent(in), optional :: iAccuracy_

    ! Local declarations:
    INTEGER :: s, t
    integer :: iUpdateTypeLocal, iAccuracy
    TYPE(silja_grid) :: storage_grid_local
    TYPE(silja_stack), POINTER  :: stackptrLocal
    real :: fMissingValue

    if (present(iUpdateType)) then
      iUpdateTypeLocal = iUpdateType
    else
      ! This is used to be ifUpdateAllowed, which defaulted to false...
      iUpdateTypeLocal = create_if_absent
    end if

    if(present(fMissingValue_))then
      fMissingValue = fMissingValue_
    else
      fMissingValue = real_missing
    endif
    if(present(iAccuracy_))then
      iAccuracy = iAccuracy_
    else
      iAccuracy = 5
    endif
    if (present(storage_grid)) then
      storage_grid_local = storage_grid
    else
      ! This is used to be ifUpdateAllowed, which defaulted to false...
      storage_grid_local = grid_missing
    end if
    !
    ! Check input parameters
    !
    IF (.NOT.defined(id)) THEN
      CALL set_error('undefined id given','dq_store_2d')
      RETURN
    END IF

    IF (fu_number_of_gridpoints(fu_grid(id)) > SIZE(grid_data)) THEN
      CALL set_error('number_of_gridpoints bigger than grid_data','dq_store_2d')
      RETURN
    END IF
!call msg('')
!call msg('store_2d quantity:' + fu_quantity_short_string(fu_quantity(id)))
!call report(fu_grid(id))

    ! Get right stack
    if(present(stackptr))then 
        stackptrLocal => stackptr
    else

      s = fu_met_src_storage_index(miniMarketPtr, stack_type, fu_met_src(id))
!        call msg('Store 2',s)
      select case(stack_type)

        case(single_time_stack_flag)
          stackptrLocal => miniMarketPtr%stacks_singleTime(s)%ptr

        case(multi_time_stack_flag)

         t = fu_time_storage_index(miniMarketPtr, s, fu_valid_time(id), .false.)

!        call msg('Store 3',t)

         if(s < 1 .or. s > size(miniMarketPtr%stacks_multiTime,1))then
           call set_error('Bad met_src, can not find proper stack','dq_store_2d')
           return
         end if
         stackptrLocal =>  miniMarketPtr%stacks_multiTime(s, t)%ptr


      case default
        call msg('Strange stack type:', stack_type)
        call set_error('Strange stack type:', 'dq_store_2d')
        return
      end select
    endif
!        call msg('Store 4')
    CALL put_field_to_stack(id,&
                              & grid_data,&
                              & stackptrLocal, &
                              & .True., & ! Force the grid
                              & iUpdateTypeLocal, &
                              & ifRandomise, &
                              & storage_grid = storage_grid_local, &
                              & storage_area = area_missing, &
                              & fMissingValue_ = fMissingValue, &
                              & iAccuracy = iAccuracy, &
                              & iOutside = fu_outGridInterpType(fu_quantity(id))) 

!        call msg('Store 5')


    if(.not. present(stackptr))then
        if(.not.error)miniMarketPtr%minimarket_empty = .false.
    endif 

  END SUBROUTINE dq_store_2d



  ! ***************************************************************

  ! ***************************************************************
  !
  !
  !       Return some information about the present 
  !       contents of the supermarket.
  !
  !
  ! ***************************************************************
  ! ***************************************************************


  ! ***************************************************************


  FUNCTION fu_supermarket_cpu_usage()

    ! Description:
    ! Returns the total cpu-time used by the functions of this module
    ! during the run of program.
    ! 
    IMPLICIT NONE
    !
    ! The return value of this function:
    TYPE(silja_interval) :: fu_supermarket_cpu_usage

    fu_supermarket_cpu_usage = fu_set_interval_sec(cpu_usage)

  END FUNCTION fu_supermarket_cpu_usage


  ! ***************************************************************


  SUBROUTINE supermarket_times(miniMarketPtr, met_src, valid_times, number_of_times, analysis_time)
    !
    ! Finds all observation times (valid times) that fields
    ! are found in supermarket.
    !
    ! The analysis time returned is the first found. Use it when
    ! there is only data from one forecast run (post-processing or
    ! other tools) but not in dispersion model runs.
    !
    IMPLICIT NONE

    ! Imported parameters with intent IN:
    type(mini_market_of_stacks), intent(in) :: miniMarketPtr
    type(meteo_data_source), INTENT(in) :: met_src

    ! Imported parameters with intent OUT:
    TYPE(silja_time), DIMENSION(:), INTENT(inout) :: valid_times
    INTEGER, INTENT(out) :: number_of_times

    ! Optional parameters with intent OUT:
    TYPE(silja_time), INTENT(out), OPTIONAL :: analysis_time

    ! Local declarations:
    TYPE(silja_stack), POINTER :: stack
    INTEGER :: i, n_times, stack_from, stack_to

    valid_times = time_missing

    if (.not. miniMarketPtr%stack_multiTime_exists) then
            call msg("No multitime in:"+ fu_name(miniMarketPtr))
            number_of_times = 0
            return
    endif

    IF (SIZE(valid_times) < SIZE(miniMarketPtr%obstimes,2)) THEN
      CALL set_error('valid_times vector too small','supermarket_times')
      RETURN
    END IF

    if(met_src == met_src_missing)then
      stack_from = 1
      stack_to=size(miniMarketPtr%stacks_multiTime,1)
    else
      stack_from = fu_met_src_storage_index(miniMarketPtr, multi_time_stack_flag, met_src)
      stack_to = stack_from
      if(stack_from < 1 .or.stack_from > size(miniMarketPtr%stacks_multiTime,1))return
    end if

    n_times = 1
    do i=stack_from, stack_to
      stack => miniMarketPtr%stacks_multiTime(i,1)%ptr
      if(.not.defined(stack))cycle

      valid_times(n_times:SIZE(miniMarketPtr%obstimes,2)) &
           & = miniMarketPtr%obstimes(i,1:SIZE(miniMarketPtr%obstimes,2)-n_times+1)

      do while(defined(valid_times(n_times)))
        n_times = n_times+1
      end do

    end do
    !
    ! Sort the gathered time stamps
    !
    CALL times_to_ascending_order(valid_times, number_of_times)

    IF (PRESENT(analysis_time)) THEN
      ! In v5d-usage we assume that all data is from one forecats set,
      ! so one analysis time only:
      stack => miniMarketPtr%stacks_multiTime(1,1)%ptr
      analysis_time = fu_analysis_time(fu_first_field_from_stack(stack))
    END IF

  END SUBROUTINE supermarket_times


    ! ***************************************************************


  SUBROUTINE supermarket_met_srcs(miniMarketPtr, stackType, metSrcs, nMetSrcs)

    ! Description:
    ! Finds all observation times (valid times) that fields
    ! are found in supermarket.
    !
    ! The analysis time returned is the first found. Use it when
    ! there is only data from one forecast run (post-processing or
    ! other tools) but not in dispersion model runs.
    !
    IMPLICIT NONE

    ! Imported parameters with intent IN:
    type(mini_market_of_stacks), intent(in) :: miniMarketPtr
    integer, intent(in) :: stackType

    ! Imported parameters with intent OUT:
    type(meteo_data_source),DIMENSION(:), INTENT(out) :: metSrcs
    INTEGER, INTENT(out) :: nMetSrcs

    ! Local declarations:
    TYPE(silja_stack), POINTER :: stack
    INTEGER :: iMetSrc

    metSrcs = met_src_missing

    select case(stackType)
      case(single_time_stack_flag)

      IF (SIZE(metSrcs) < size(miniMarketPtr%stacks_singletime)) THEN
        CALL set_error('metSrcs vector too small','supermarket_met_srcs')
        RETURN
      END IF

      nMetSrcs = size(miniMarketPtr%stacks_singletime)
      do iMetSrc = 1, nMetSrcs
        stack => miniMarketPtr%stacks_singleTime(iMetSrc)%ptr
        if(.not.defined(stack))cycle
        metSrcs(iMetSrc) = fu_met_src(miniMarketPtr%stacks_singleTime(iMetSrc)%ptr)
      end do

      case(multi_time_stack_flag)

      IF (SIZE(metSrcs) < size(miniMarketPtr%stacks_multitime,1)) THEN
        CALL set_error('metSrcs vector too small','supermarket_met_srcs')
        RETURN
      END IF

      nMetSrcs = size(miniMarketPtr%stacks_multitime,1)
      do iMetSrc = 1, nMetSrcs
        stack => miniMarketPtr%stacks_multiTime(iMetSrc, 1)%ptr
        if(.not.defined(stack))cycle
        metSrcs(iMetSrc) = fu_met_src(miniMarketPtr%stacks_multiTime(iMetSrc,1)%ptr)
      end do
    end select

  END SUBROUTINE supermarket_met_srcs


  ! ***************************************************************


  SUBROUTINE supermarket_quantities(miniMarketPtr, met_src, &
                                  & stack_type,&
                                  & quantities,&
                                  & number_of_quantities)
    ! 
    ! Finds all quantities that fields are found for the 
    ! given met_src in supermarket.
    !
    IMPLICIT NONE

    ! Imported parameters with intent IN:
    type(mini_market_of_stacks), intent(in) :: miniMarketPtr
    type(meteo_data_source), INTENT(in) :: met_src
    integer, intent(in) :: stack_type

    ! Imported parameters with intent OUT:
    INTEGER, DIMENSION(:), INTENT(out) :: quantities
    INTEGER, INTENT(out) :: number_of_quantities

    ! Local declarations:
    INTEGER :: i, j, stack_from, stack_to
    TYPE(silja_stack), POINTER  ::  stack
    INTEGER, dimension(max_quantities) :: qTmp

    quantities(1) = int_missing
    number_of_quantities = 0

    select case(stack_type)

      case(multi_time_stack_flag )

        IF (.NOT. miniMarketPtr%stack_multiTime_initialized) THEN
          CALL set_error('minimarket not init','supermarket_quantities')
          RETURN
        END IF

        if(met_src == met_src_missing)then
          stack_from = 1
          stack_to=size(miniMarketPtr%stacks_multiTime,1)
        else
          stack_from = fu_met_src_storage_index(miniMarketPtr, stack_type, met_src)
          stack_to = stack_from
          if(stack_from < 1 .or.stack_from > size(miniMarketPtr%stacks_multiTime,1))return
        end if

        do i=stack_from, stack_to
          stack => miniMarketPtr%stacks_multiTime(i,1)%ptr
          if(.not.defined(stack))cycle
          CALL stack_quantities(stack, qTmp, j)
          j=fu_merge_int_arrays(qTmp, quantities, .false.)
          number_of_quantities = number_of_quantities + j
        end do

      case(single_time_stack_flag )

        IF (.NOT. miniMarketPtr%stack_singleTime_initialized) THEN
          CALL set_error('minimarket not init','supermarket_quantities')
          RETURN
        END IF

        if(met_src == met_src_missing)then
          stack_from = 1
          stack_to=size(miniMarketPtr%stacks_singleTime)
        else
          stack_from = fu_met_src_storage_index(miniMarketPtr, stack_type, met_src)
          stack_to = stack_from
          if(stack_from < 1 .or.stack_from > size(miniMarketPtr%stacks_singleTime))return
        end if

        do i=stack_from, stack_to
          stack => miniMarketPtr%stacks_singleTime(i)%ptr
          if(.not.defined(stack))cycle
          CALL stack_quantities(stack, qTmp, j)
          j=fu_merge_int_arrays(qTmp, quantities, .false.)
          number_of_quantities = number_of_quantities + j
        end do

      case default
        call msg('Strange stack_type', stack_type)
        call set_error('Strange stack_type', 'supermarket_quantities')

    end select

  END SUBROUTINE supermarket_quantities


  ! ***************************************************************


  SUBROUTINE supermarket_variables(miniMarketPtr, met_src, &
                                 & stack_type,&
                                 & quantities,&
                                 & arSpecies,&
                                 & number_of_variables)
    ! 
    ! Finds all quantities+species combinations that fields are found for the 
    ! given met_src in supermarket.
    !
    IMPLICIT NONE

    ! Imported parameters
    type(mini_market_of_stacks), intent(in) :: miniMarketPtr
    type(meteo_data_source), INTENT(in) :: met_src
    integer, intent(in) :: stack_type
    INTEGER, DIMENSION(:), intent(out) :: quantities
    type(silam_species), dimension(:), intent(out) :: arSpecies
    INTEGER, INTENT(out) :: number_of_variables

    ! Local declarations:
    INTEGER :: i, j, stack_from, stack_to
    TYPE(silja_stack), POINTER  ::  stack
    INTEGER, dimension(max_variables) :: qTmp
    type(silam_species), dimension(max_variables) :: arSpeciesTmp

    quantities(:) = int_missing
    arSpecies(:) = species_missing
    number_of_variables = 0
    !
    ! Let's allow stack_type be undefined
    !
    if(stack_type /= single_time_stack_flag)then                 ! multi-time stack included

      IF (miniMarketPtr%stack_multiTime_initialized) THEN

        if(met_src == met_src_missing)then
          stack_from = 1
          stack_to=size(miniMarketPtr%stacks_multiTime,1)
        else
          stack_from = fu_met_src_storage_index(miniMarketPtr, stack_type, met_src)
          stack_to = stack_from
          if(stack_from < 1 .or.stack_from > size(miniMarketPtr%stacks_multiTime,1))return
        end if

        do i=stack_from, stack_to
          stack => miniMarketPtr%stacks_multiTime(i,1)%ptr
          if(.not.defined(stack))cycle
          CALL stack_variables(stack, qTmp, arSpeciesTmp, j)
          if(j>0)then 
             call merge_variable_arrays(quantities, qTmp, arSpecies, arSpeciesTmp, number_of_variables)
          endif
        end do
      endif  ! if multitime initialised
    endif  ! if not only singletime

    if(stack_type /= multi_time_stack_flag)then                       ! single-time stack included

      IF (miniMarketPtr%stack_singleTime_initialized) THEN

        if(met_src == met_src_missing)then
          stack_from = 1
          stack_to=size(miniMarketPtr%stacks_singleTime)
        else
          stack_from = fu_met_src_storage_index(miniMarketPtr, stack_type, met_src)
          stack_to = stack_from
          if(stack_from < 1 .or.stack_from > size(miniMarketPtr%stacks_singleTime))return
        end if

        do i=stack_from, stack_to
          stack => miniMarketPtr%stacks_singleTime(i)%ptr
          if(.not.defined(stack))cycle
          CALL stack_variables(stack, qTmp, arSpeciesTmp, j)
          if (j>0)  call merge_variable_arrays(quantities, qTmp, arSpecies, arSpeciesTmp, number_of_variables)
        end do

      END IF   ! if singletime initialised
    end if  ! if not inly multitime

    CONTAINS
    
    subroutine merge_variable_arrays(quantities, extraQuantities, arSpecies, arExtraSpecies, &
                                   & number_of_variables)
      !
      ! Merges the extra quantity-species combinations to the main arrays. Duplication of the
      ! combinations is not allowed.
      !    
      implicit none

      ! Imported parameters
      integer, dimension(:), intent(inout) :: quantities
      integer, dimension(:), intent(in) ::  extraQuantities
      type(silam_species), dimension(:),intent(inout) :: arSpecies
      type(silam_species), dimension(:),intent(in) :: arExtraSpecies
      integer, intent(inout) :: number_of_variables
      
      ! Local variables
      integer :: iQ, iQExtra
      logical :: ifAdd
     
      do iQExtra = 1, size(extraQuantities)          ! Scan extra variables
        if(iQExtra > size(arExtraSpecies))exit
        if(extraQuantities(iQExtra) == int_missing)exit
        ifAdd = .true.                               ! Check for duplicate one
        do iQ = 1, size(quantities)
          if(iQ > size(arSpecies))exit
          if(quantities(iQ) == int_missing)exit
          if(quantities(iQ) == extraQuantities(iQExtra))then
            if(arSpecies(iQ) == arExtraSpecies(iQExtra))then
              ifAdd = .false.
              exit
            endif
          endif
        end do   ! iQ
        if(ifAdd)then                                ! Adding unique variable
          if(number_of_variables >= size(quantities))then
            call set_error('Too small quantity array','fu_merge_variable_arrays')
            return
          endif
          number_of_variables = number_of_variables + 1
          arSpecies(number_of_variables) = arExtraSpecies(iQExtra)
          quantities(number_of_variables) = extraQuantities(iQExtra)
        endif
      end do  ! iExtra
      
    end subroutine merge_variable_arrays

  END SUBROUTINE supermarket_variables


  ! ***************************************************************


  SUBROUTINE supermarket_3d_quantities(miniMarketPtr, met_src, stack_type, &
                                     & quantities, number_of_quantities)

    ! Description:
    ! Finds all quantities that 3d scalar fields
    ! are found for the given met_src in supermarket.
    ! Used in v5d-conversion.
    !
    IMPLICIT NONE

    ! Imported parameters with intent IN:
    type(mini_market_of_stacks), intent(in) :: miniMarketPtr
    type(meteo_data_source), INTENT(in) :: met_src
    integer, intent(in) :: stack_type
    INTEGER, intent(out) :: number_of_quantities

    ! Imported parameters with intent OUT:
    INTEGER, DIMENSION(:), INTENT(out) :: quantities

    ! Local declarations:
    INTEGER :: i, j, stack_from, stack_to
    TYPE(silja_stack), POINTER  ::  stack
    INTEGER, dimension(max_quantities) :: qTmp

    quantities(:) = int_missing
    number_of_quantities = 0

    select case(stack_type)

      case(multi_time_stack_flag )

        IF (.NOT. miniMarketPtr%stack_multiTime_initialized) THEN
          CALL set_error('minimarket not init','supermarket_quantities')
          RETURN
        END IF

        if(met_src == met_src_missing)then
          stack_from = 1
          stack_to=size(miniMarketPtr%stacks_multiTime,1)
        else
          stack_from = fu_met_src_storage_index(miniMarketPtr, stack_type, met_src)
          stack_to = stack_from
          if(stack_from < 1 .or.stack_from > size(miniMarketPtr%stacks_multiTime,1))return
        end if

        do i=stack_from, stack_to
          stack => miniMarketPtr%stacks_multiTime(i,1)%ptr
          if(.not.defined(stack))cycle
          CALL stack_3d_quantities(stack, qTmp, j)
          j = fu_merge_int_arrays(qTmp, quantities, .false.)
          number_of_quantities = number_of_quantities + j
        end do

      case(single_time_stack_flag )

        IF (.NOT. miniMarketPtr%stack_singleTime_initialized) THEN
          CALL set_error('minimarket not init','supermarket_quantities')
          RETURN
        END IF

        if(met_src == met_src_missing)then
          stack_from = 1
          stack_to=size(miniMarketPtr%stacks_singleTime)
        else
          stack_from = fu_met_src_storage_index(miniMarketPtr, stack_type, met_src)
          stack_to = stack_from
          if(stack_from < 1 .or.stack_from > size(miniMarketPtr%stacks_singleTime))return
        end if

        do i=stack_from, stack_to
          stack => miniMarketPtr%stacks_singleTime(i)%ptr
          if(.not.defined(stack))cycle
          CALL stack_3d_quantities(stack, qTmp, j)
          j = fu_merge_int_arrays(qTmp, quantities, .false.)
          number_of_quantities = number_of_quantities + j
        end do

      case default
        call msg('Strange stack type:',stack_type)
        call set_error('Strange stack type:','supermarket_3d_quantities')

    end select

  END SUBROUTINE supermarket_3d_quantities


  ! ***************************************************************

  SUBROUTINE supermarket_2d_quantities(miniMarketPtr, met_src, stack_type, &
                                     & quantities_2d, number_of_2d_quantities)

    ! Description:
    ! Finds all possible scalar quantities that are found in supermarket
    ! but only on one, single level (like surface stuff etc.)
    !
    IMPLICIT NONE

    ! Imported parameters 
    type(mini_market_of_stacks), intent(in) :: miniMarketPtr
    type(meteo_data_source), INTENT(in) :: met_src
    integer, intent(in) :: stack_type
    integer, intent(out) :: number_of_2d_quantities
    INTEGER, DIMENSION(:), INTENT(out) :: quantities_2d

    ! Local declarations:
    INTEGER :: i, j
    TYPE(silja_stack), POINTER  ::  stack
    INTEGER :: quantitity
    INTEGER, DIMENSION(max_quantities) :: quantities, quantities_3d
    INTEGER :: number_of_quantities, number_of_3d_quantities

    quantities_2d(:) = int_missing
    number_of_2d_quantities = 0

    !    
    ! 1. Which quantities are found in supermarket?
    ! 
    CALL supermarket_quantities(miniMarketPtr, met_src, stack_type, quantities, number_of_quantities)
    IF (error) RETURN
    IF (number_of_quantities == 0) RETURN

    !    
    ! 2. Which 3D-quantities are found in supermarket?
    !
    CALL supermarket_3d_quantities(miniMarketPtr, met_src, stack_type, &
                                 & quantities_3d, number_of_3d_quantities)
    IF (error) RETURN

    !    
    ! 3. Which quantities are found in for a single 2D level only?
    !
    CALL quantities_single_level_only(quantities,&
                                    & number_of_quantities,&
                                    & quantities_3d,&
                                    & number_of_3d_quantities,&
                                    & quantities_2d,&
                                    & number_of_2d_quantities)
    IF (error) RETURN

  CONTAINS !*****************************************

    SUBROUTINE quantities_single_level_only(quantities,&
                                          & number_of_quantities,&
                                          & quantities_3d,&
                                          & number_of_3d_quantities,&
                                          & quantities_single,&
                                          & number_of_single_quantities)

      IMPLICIT NONE

      ! Imported parameters with intent IN:
      INTEGER, DIMENSION(:), INTENT(in) :: quantities
      INTEGER, DIMENSION(:), INTENT(in) :: quantities_3d
      INTEGER, INTENT(in) :: number_of_quantities
      INTEGER, INTENT(in) :: number_of_3d_quantities

      ! Imported parameters with intent OUT:
      INTEGER, DIMENSION(:), INTENT(out) :: quantities_single
      INTEGER, INTENT(out) :: number_of_single_quantities

      ! Local declarations:
      INTEGER :: i, j
      LOGICAL :: found_in_3d

      quantities_single = int_missing
      number_of_single_quantities = 0

      outer: DO i = 1, number_of_quantities

        found_in_3d = .false.

        inner: DO j = 1, number_of_3d_quantities
          IF (quantities(i) == quantities_3d(j)) found_in_3d = .true.
        END DO inner

        IF (.NOT.found_in_3d) THEN
          number_of_single_quantities = number_of_single_quantities + 1
          quantities_single(number_of_single_quantities) = quantities(i)
        END IF

      END DO outer

    END SUBROUTINE quantities_single_level_only

  END SUBROUTINE supermarket_2d_quantities


  ! ***************************************************************


  LOGICAL FUNCTION fu_field_in_sm(miniMarketPtr, met_src,&
                                & quantity,&
                                & time,&
                                & level,&
                                & look_for_3d,&
                                & permanent, &
                                & species)

    ! Description:
    ! Returns true value if the given field is already in supermarket.
    ! Finds all fieldtypes: scalar and wind, both 2D and 3D.
    !
    IMPLICIT NONE
    !
    ! Imported parameters with intent(in):
    type(mini_market_of_stacks), intent(in) :: miniMarketPtr
    type(meteo_data_source), INTENT(in) :: met_src
    INTEGER, INTENT(in) :: quantity
    TYPE(silja_time), INTENT(in) :: time ! must be exact obstime
    TYPE(silja_level), INTENT(in) :: level ! if missing, then search
    ! only by quantity and time
    LOGICAL, INTENT(in) :: look_for_3d ! search for 3d field only
    LOGICAL, INTENT(in) :: permanent
    type(silam_species), intent(in), optional :: species

    ! Local declarations:
    TYPE(silja_stack), POINTER :: stack
    TYPE(silja_field_id) :: id
    TYPE(silja_field) , POINTER :: field
    TYPE(silja_3d_field), POINTER :: field_3d
    TYPE(silja_windfield) , POINTER :: windfield
    TYPE(silja_3d_windfield), POINTER :: windfield_3d
    integer :: stack_from, stack_to, i, j
    
    !
    ! 1. Basic checks.
    !
    fu_field_in_sm = .false.
    nullify(field)

    IF (miniMarketPtr%minimarket_empty) THEN
      RETURN
    END IF

    !
    ! 2. Look for a permanent field.
    !
    IF (permanent) THEN
      !
      ! Permanent stack
      !
      IF (.NOT. miniMarketPtr%stack_singleTime_initialized) THEN
        RETURN
      END IF
      
      if(met_src == met_src_missing)then
        stack_from = 1
        stack_to=size(miniMarketPtr%stacks_singleTime)
      else
          stack_from = fu_met_src_storage_index(miniMarketPtr, single_time_stack_flag, met_src)
          stack_to = stack_from
          if(stack_from < 1 .or. stack_from > size(miniMarketPtr%stacks_singleTime))return
      end if

      if (time == time_missing) then ! Do not care about time
           do i=stack_from, stack_to
              stack => miniMarketPtr%stacks_singleTime(i)%ptr
              if(.not.defined(stack))cycle
              IF (look_for_3d) THEN
                 CALL find_field_3d_from_stack(met_src,&
                                              & quantity,&
                                              & time,&
                                              & stack,&
                                              & field_3d,&
                                              & fu_field_in_sm)
              ELSE
                  CALL find_field_from_stack(met_src, &
                                        & quantity,&
                                        & time_missing,&
                                        & stack,&
                                        & field,&
                                        & fu_field_in_sm, &
                                        & species)

              ENDIF
!                call msg("Trying: "+fu_quantity_string(quantity) )
!                call report(stack)
!                                DO j = 1, SIZE(stack%fields)
!                        call report(stack%fields(j)%fp)
!                enddo
!                if (fu_field_in_sm) then
!                        call msg("Found: "+fu_quantity_string(fu_quantity(field)))
!                else
!                        call msg("Not Found: "+fu_quantity_string(quantity))
!                endif
              if(fu_field_in_sm)exit
           end do
      else ! respect field time... Is it needed in single-time stack at all??? R. 
              id = fu_set_field_id_simple(met_src, quantity, time_missing, level)
              IF (error) RETURN
              if(present(species))call set_species(id, species)

              do i=stack_from, stack_to
                stack => miniMarketPtr%stacks_singleTime(i)%ptr
                if(.not.defined(stack))cycle
                CALL find_field_from_stack(stack, id, field, fu_field_in_sm)
                if(fu_field_in_sm)exit
              end do
      endif
      return
      
    ELSE 
      !
      ! Non-permanent stack
      !
      IF (.NOT. miniMarketPtr%stack_multiTime_initialized) THEN
         RETURN
      END IF
      !
      ! 3. Where the time-dependent data shoud be
      !
      stack => fu_closest_sm_met_src_time(miniMarketPtr, met_src, time, single_time, .false.) 
                ! !not mandatory -- do not raise error at empty stack
      IF (error) RETURN
      IF (.NOT. associated(stack)) RETURN
      !
      ! 4. Look for a windfield.
      !
      IF (quantity == wind_flag) THEN
 
        IF (look_for_3d) THEN

          CALL find_wind_3d_from_stack(met_src,&
                                     & time,&
                                     & stack,&
                                     & windfield_3d,&
                                     & fu_field_in_sm)
        ELSE
          IF (.NOT.defined(level)) THEN
            CALL set_error('please define level for wind 2D searching',&
                         & 'fu_field_in_sm')
            RETURN
          END IF

          id = fu_set_windfield_id_simple(quantity, met_src, time, level)
          IF (error) RETURN

          CALL find_wind_from_stack(stack, id, windfield, fu_field_in_sm)

        END IF

        RETURN
      END IF ! Wind
      !
      ! 5. Look for a scalar field.
      !
      IF (look_for_3d) THEN

        if(present(species))then
          CALL find_field_3d_from_stack(met_src,&
                                      & quantity,&
                                      & time,&
                                      & stack,&
                                      & field_3d,&
                                      & fu_field_in_sm, &
                                      & species=species)
        else
          CALL find_field_3d_from_stack(met_src,&
                                      & quantity,&
                                      & time,&
                                      & stack,&
                                      & field_3d,&
                                      & fu_field_in_sm)
        endif
      ELSE

        id = fu_set_field_id_simple(met_src, quantity, time, level)
        if(present(species)) call set_species(id, species)
        IF (error) RETURN

        CALL find_field_from_stack(stack, id, field, fu_field_in_sm)

      END IF
    endif  ! non-permanent


  END FUNCTION fu_field_in_sm




  ! ***************************************************************

  ! ***************************************************************
  !
  !
  !      Private functions and subroutines.
  !
  !
  ! ***************************************************************
  ! ***************************************************************

  SUBROUTINE find_field_from_hit_list(miniMarketPtr, field_id, field, found_in_hit_list)

    ! Description:
    ! Checks if the desired field is in the top-20 hit list of last
    ! retrieve fields.
    ! 
    IMPLICIT NONE

    ! Imported parameters with intent IN:
    type(mini_market_of_stacks), intent(in) :: miniMarketPtr
    TYPE(silja_field_id), INTENT(in) :: field_id

    ! Imported parameters with intent out:
    LOGICAL, INTENT(out) :: found_in_hit_list

    ! Imported parameters with intent INOUT or POINTER:
    TYPE(silja_field), POINTER :: field

    ! Local declarations:
    INTEGER :: i

    found_in_hit_list = .false.

    DO i = 1, SIZE(miniMarketPtr%hit_list)

      IF (.NOT.ASSOCIATED(miniMarketPtr%hit_list(i)%fp)) EXIT

      if(.not.defined(miniMarketPtr%hit_list(i)%fp))cycle

      IF (field_id == fu_id(miniMarketPtr%hit_list(i)%fp)) THEN
        !    PRINT *, ' Found in top-20.'
        found_in_hit_list = .true.
        field => miniMarketPtr%hit_list(i)%fp
        EXIT
      END IF

    END DO

  END SUBROUTINE find_field_from_hit_list


  ! ***************************************************************


  SUBROUTINE put_field_to_hit_list(miniMarketPtr, field)

    ! Description:
    ! 
    ! Method:
    ! 
    IMPLICIT NONE

    ! Imported parameters with intent INOUT or POINTER:
    TYPE(silja_field), POINTER :: field
    type(mini_market_of_stacks), intent(inout) :: miniMarketPtr

    ! Local declarations:
    INTEGER :: i

    ! The correct location in top_list:
    miniMarketPtr%hit_list_pointer = miniMarketPtr%hit_list_pointer + 1
    IF (miniMarketPtr%hit_list_pointer > SIZE(miniMarketPtr%hit_list)) miniMarketPtr%hit_list_pointer = 1

    miniMarketPtr%hit_list(miniMarketPtr%hit_list_pointer)%fp => field


  END SUBROUTINE put_field_to_hit_list


  ! ***************************************************************


  SUBROUTINE find_wind_from_hit_list(miniMarketPtr, field_id, field, found_in_hit_list)
    !
    ! Description:
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

    ! Imported parameters with intent IN:
    TYPE(silja_field_id), INTENT(in) :: field_id
    type(mini_market_of_stacks), intent(in) :: miniMarketPtr

    ! Imported parameters with intent out:
    LOGICAL, INTENT(out) :: found_in_hit_list

    ! Imported parameters with intent INOUT or POINTER:
    TYPE(silja_windfield), POINTER :: field

    ! Local declarations:
    INTEGER :: i

    found_in_hit_list = .false.

    DO i = 1, SIZE(miniMarketPtr%wind_list)
      IF (.NOT.ASSOCIATED(miniMarketPtr%wind_list(i)%fp)) EXIT
      field => miniMarketPtr%wind_list(i)%fp
      IF (field_id == fu_id(miniMarketPtr%wind_list(i)%fp)) THEN
        found_in_hit_list = .true.
        field => miniMarketPtr%wind_list(i)%fp
        EXIT
      END IF
    END DO

  END SUBROUTINE find_wind_from_hit_list


  ! ***************************************************************


  SUBROUTINE put_wind_to_hit_list(miniMarketPtr, field)

    ! Description:
    ! 
    ! Method:
    ! 
    IMPLICIT NONE

    ! Imported parameters with intent INOUT or POINTER:
    TYPE(silja_windfield), POINTER :: field
    type(mini_market_of_stacks), intent(inout) :: miniMarketPtr

    ! Local declarations:
    INTEGER :: i

    ! The correct location in top_list:
    miniMarketPtr%wind_list_pointer = miniMarketPtr%wind_list_pointer + 1
    IF (miniMarketPtr%wind_list_pointer > SIZE(miniMarketPtr%wind_list)) miniMarketPtr%wind_list_pointer = 1

    miniMarketPtr%wind_list(miniMarketPtr%wind_list_pointer)%fp => field


  END SUBROUTINE put_wind_to_hit_list


  ! ***************************************************************


  INTEGER FUNCTION fu_met_src_storage_index(miniMarketPtr, stack_type, met_src) result(i)
    ! 
    ! Returns the correct first index of the stack-matrix, that contains
    ! the data for the given met_src.
    !
    ! The stack-array is searched unitl the correct is found, or
    ! the first undefined.
    !
    IMPLICIT NONE
    !
    ! Imported parameters with intent(in):
    type(meteo_data_source), INTENT(in) :: met_src
    type(mini_market_of_stacks), intent(in) :: miniMarketPtr
    integer, intent(in) :: stack_type
    ! Local variables
    type(silja_stack), pointer :: stackPtr
    INTEGER :: j

    ! For time dependent stack
    !
    if(stack_type == multi_time_stack_flag)then

      ! First look for exact match
      !
      DO i = 1, SIZE(miniMarketPtr%stacks_multiTime, 1)
        stackPtr => miniMarketPtr%stacks_multiTime(i,1)%ptr
        IF (fu_single_met_src_stack(stackPtr) .and. met_src == fu_met_src(stackPtr)) RETURN
      END DO
     
      ! No match found. Look for multiple met_source stack or met_source missing stack
      !
      DO i = 1, SIZE(miniMarketPtr%stacks_multiTime, 1)
        stackPtr => miniMarketPtr%stacks_multiTime(i,1)%ptr
        IF ((.not. fu_single_met_src_stack(stackPtr)) .or. fu_met_src(stackPtr) == met_src_missing) RETURN
      END DO
 
      ! Still not found - set error
      !
      call msg('')
      call msg('Sources in the supermarket:')
      DO i = 1, SIZE(miniMarketPtr%stacks_multiTime, 1)
        do j =1, SIZE(miniMarketPtr%stacks_multiTime, 2)
          stackPtr => miniMarketPtr%stacks_multiTime(i,j)%ptr
        call msg(fu_name(fu_met_src(stackPtr))) 
        end do
      END DO

      CALL set_error('Below met_src is not in supermarket','fu_met_src_storage_index')
      call msg(fu_name(met_src))

      i=0

    ! And the same for single-time stack
    !
    elseif(stack_type == single_time_stack_flag)then

      ! First look for exact match
      DO i = 1, SIZE(miniMarketPtr%stacks_singleTime)
        stackPtr => miniMarketPtr%stacks_singleTime(i)%ptr
        IF (fu_single_met_src_stack(stackPtr) .and. met_src == fu_met_src(stackPtr)) RETURN
      END DO

      ! No match found. Look for multiple met_source stack or met_source missing stack
      !
      DO i = 1, SIZE(miniMarketPtr%stacks_singleTime)
        stackPtr => miniMarketPtr%stacks_singleTime(i)%ptr
        IF ((.not. fu_single_met_src_stack(stackPtr)) .or. fu_met_src(stackPtr) == met_src_missing) RETURN
      END DO
 
      ! Still not found - set error
      !
      call msg('')
      call msg('Sources in the supermarket:')
      DO i = 1, SIZE(miniMarketPtr%stacks_singleTime)
        stackPtr => miniMarketPtr%stacks_singleTime(i)%ptr
        call msg(fu_name(fu_met_src(stackPtr))) 
      END DO

      CALL set_error('Below met_src is not in supermarket','fu_met_src_storage_index')
      call msg(fu_name(met_src))

      i=0

    else
      call set_error('Strange stack type','fu_met_src_storage_index')
    endif

  END FUNCTION fu_met_src_storage_index


  ! ***************************************************************


  INTEGER FUNCTION fu_time_storage_index(miniMarketPtr, s, time, ifPreserveMetSrcValue) result(t)
    ! 
    ! Returns a time-index (2.) of the stack, into which data for the
    ! given valid time can be stored.
    ! If such a stack is not found then a index of the first empty one is returned.
    ! 
    ! If correct time or empty stack is not found, then the first stack is set
    ! empty and pointer to it is returned.
    !
    IMPLICIT NONE
    !
    ! Imported parameters with intent(in):
    INTEGER, INTENT(in) :: s ! the met_src index, along which the the time is searched
    TYPE(silja_time), INTENT(in) :: time
    logical, intent(in) :: ifPreserveMetSrcValue
    type(mini_market_of_stacks), intent(inout) :: miniMarketPtr

    ! Local declarations:
    TYPE(silja_time) :: stack_valid
    TYPE(silja_stack), POINTER :: stack
    type(meteo_data_source) :: mds

    !
    ! Find the stack containing the requested time
    !
    DO t = 1, SIZE(miniMarketPtr%stacks_multiTime, 2)

      stack => miniMarketPtr%stacks_multiTime(s, t)%ptr
      stack_valid = fu_valid_time(stack)
      IF (error) RETURN

      ! Skip the empty stacks
      IF (.NOT.defined(stack_valid)) cycle

      ! Correct found?
      IF (time == stack_valid)RETURN

    END DO

    !
    ! No time matching the requested one. Find the empty stack
    !
    DO t = 1, SIZE(miniMarketPtr%stacks_multiTime, 2)

      stack => miniMarketPtr%stacks_multiTime(s, t)%ptr
      stack_valid = fu_valid_time(stack)
      IF (error) RETURN

      ! Empty one?
      IF (.NOT.defined(stack_valid)) RETURN
    END DO

    !
    ! Neither correct-time nor empty are found, so start overwrite the "oldest".
    ! Note that the term "oldest" is not trivial since there can be forward and adjoint runs.
    !
    IF (miniMarketPtr%replace_earliest) THEN
      t = fu_earliest_stack(miniMarketPtr, s)  ! forward run, the earliest are to be overwritten
    ELSE
      t = fu_latest_stack(miniMarketPtr, s)  ! adjoint run, the latest are to be overwritten
    END IF
    stack => miniMarketPtr%stacks_multiTime(s, t)%ptr
    IF (error) RETURN
    call msg('Overwriting the stack with valid time:' + fu_str(fu_valid_time(stack)))
    !
    ! Whle emptying the stack, may need to keep the MDS and restore it afterwards.
    !
    if(ifPreserveMetSrcValue) mds = fu_met_src(miniMarketPtr%stacks_multiTime(s,t)%ptr)
    
    CALL set_stack_empty(miniMarketPtr%stacks_multiTime(s,t)%ptr)
    
    if(ifPreserveMetSrcValue) call set_met_src(miniMarketPtr%stacks_multiTime(s,t)%ptr, mds)

  CONTAINS

    ! ************** private functions *****************


    INTEGER FUNCTION fu_earliest_stack(miniMarketPtr, s)
      !
      ! Returns a index of the stack that contains oldest data along
      ! met_src-index s.
      !
      IMPLICIT NONE
      !
      ! Imported:
      INTEGER, INTENT(in) :: s
      type(mini_market_of_stacks), intent(in) :: miniMarketPtr

      ! Local declarations:
      TYPE(silja_time) :: earliest
      INTEGER :: i
      TYPE(silja_stack), POINTER :: stack

      ! Start value:
      stack => miniMarketPtr%stacks_multiTime(s,1)%ptr
      earliest = fu_valid_time(stack)
      fu_earliest_stack = 1

      DO i = 1, SIZE(miniMarketPtr%stacks_multiTime, 2)
        stack => miniMarketPtr%stacks_multiTime(s, i)%ptr
        IF (.NOT.defined(stack)) CYCLE
        IF (fu_valid_time(stack) < earliest) THEN
          earliest = fu_valid_time(stack)
          fu_earliest_stack = i
        END IF
        IF (error) RETURN
      END DO

    END FUNCTION fu_earliest_stack

    !==========================================================================================

    INTEGER FUNCTION fu_latest_stack(miniMarketPtr, s)
      ! 
      ! Returns a index of the stack that contains latest data along
      ! met_src-index s.
      !
      IMPLICIT NONE
      !
      ! Imported:
      INTEGER, INTENT(in) :: s
      type(mini_market_of_stacks), intent(in) :: miniMarketPtr

      ! Local declarations:
      TYPE(silja_time) :: latest
      INTEGER :: i
      TYPE(silja_stack), POINTER :: stack

      ! Start value:
      stack => miniMarketPtr%stacks_multiTime(s,1)%ptr
      latest = fu_valid_time(stack)
      fu_latest_stack = 1

      DO i = 1, SIZE(miniMarketPtr%stacks_multiTime, 2)
        stack => miniMarketPtr%stacks_multiTime(s, i)%ptr
        IF (.NOT.defined(stack)) CYCLE
        IF (fu_valid_time(stack) > latest) THEN
          latest = fu_valid_time(stack)
          fu_latest_stack = i
        END IF
        IF (error) RETURN
      END DO

    END FUNCTION fu_latest_stack

  END FUNCTION fu_time_storage_index


  ! ***************************************************************


  SUBROUTINE check_sm_for_data_in_list(miniMarketPtr, list,&
                                     & data_exists_already,&
                                     & list_of_missing_stuff,&
                                     & missing_obstimes, &
                                     & wdr)
    !
    ! Cheks if fields required by the shopping list
    ! are already stored in memory in the supermarket
    ! If something on the list is missing, returns a list
    ! with missing items and obstimes.
    !
    ! Also modifies the timelimits so, that the list contains enough
    ! obstimes to cover whole timeperiod of given list. For this
    ! we add the closest obstime from start-time backwards and closest
    ! obstime from end-time forwards to the list. If either happens to be
    ! an obstime, then it is not modified.
    ! 
    ! Method:
    ! 1. If supermarket is empty, the whole given list is returned as it is.
    !
    ! 2. If observation times are found, that are required by the
    ! shopping list, but are not in supermarket obs.times, then
    ! a new list is formed with original quantities, but containing
    ! only the missing times
    !
    ! 3. If there is some data for all required obstimes, but a missing
    ! quantity is found, then the original list is returned. This is
    ! is simplified action.
    !
    ! At the moment there's no vertical level check.
    ! 
    IMPLICIT NONE
    !
    ! Imported parameters with intent(in):
    TYPE(silja_shopping_list), INTENT(in) :: list
    type(mini_market_of_stacks), intent(in) :: miniMarketPtr
    ! Imported parameters with intent(out):
    LOGICAL, INTENT(out) :: data_exists_already
    TYPE(silja_shopping_list), intent(out) ::  list_of_missing_stuff
    TYPE(silja_time), DIMENSION(:), INTENT(out) :: missing_obstimes
    type(silja_wdr), intent(in) :: wdr

    ! Local declarations:
    TYPE(silja_time), DIMENSION(max_times) :: req_obstimes, sm_times
    type(meteo_data_source) :: met_src
    INTEGER :: i, j, miss_count, number_of_sm_times
    INTEGER, DIMENSION(max_quantities) :: quantities
    TYPE(silja_field), POINTER :: field
    TYPE(silja_3d_field), POINTER :: field_3d
    TYPE(silja_time) :: earliest_in_sm, latest_in_sm
    TYPE(silja_time) :: req_earliest, req_latest

    !-------------------------------------------------------------
    !
    ! 1. Modify timelimits.
    !    -----------------

    list_of_missing_stuff = list ! set all the stuff on list as missing
    met_src = fu_met_src(list)
    req_earliest = fu_closest_obstime(fu_start_time(list), backwards, fu_obstime_interval(wdr))
    req_latest = fu_closest_obstime(fu_end_time(list), forwards, fu_obstime_interval(wdr))
    IF (error) RETURN

    CALL fix_shopping_time_boundaries(list_of_missing_stuff, req_earliest, req_latest)
    IF (error) RETURN

    req_obstimes = fu_obstimes(req_earliest, req_latest,.false., fu_obstime_interval(wdr))

    missing_obstimes = req_obstimes ! set all missing first

    !-------------------------------------------------------------
    !
    ! 2. Basic checks
    !
    data_exists_already = .false.

    IF (miniMarketPtr%minimarket_empty .or. &
     & (.NOT.(miniMarketPtr%stack_multiTime_initialized  .or. &
            & miniMarketPtr%stack_singleTime_initialized))) RETURN

!    IF (ALL(met_srcs == int_missing)) RETURN

!    IF (.NOT.ANY(met_srcs == met_src)) RETURN

    quantities = fu_quantities(list_of_missing_stuff)
    IF (error) RETURN

    !-------------------------------------------------------------
    !
    ! 3. Check for times in list that are missing in sm
    !    -----------------------------------------------

    !-------------------------------------------------------------
    !
    ! 3.1. Times that there is already data:

    CALL supermarket_times(miniMarketPtr, met_src, sm_times, number_of_sm_times)
    IF (error) RETURN
    IF (number_of_sm_times == 0) RETURN

    latest_in_sm = fu_latest_time(sm_times)
    earliest_in_sm = fu_earliest_time(sm_times)


    !-------------------------------------------------------------
    !
    ! 3.2. Find which req_obstimes are not within sm-timelimts

    missing_obstimes = time_missing
    miss_count = 0

    DO i = 1, SIZE(req_obstimes)
      IF (.NOT.defined(req_obstimes(i))) EXIT 
      IF (.NOT.fu_between_times(req_obstimes(i), earliest_in_sm, &
                              & latest_in_sm, .true.)) THEN
        !Missing obstime found
        miss_count = miss_count + 1
        missing_obstimes(miss_count) = req_obstimes(i)
      END IF
    END DO

    !-------------------------------------------------------------
    !
    ! 3.3. If missing times found, set new shop-limits and exit

    IF (miss_count > 0) THEN
      CALL fix_shopping_time_boundaries(list_of_missing_stuff,&
                                      & fu_earliest_time(missing_obstimes),&
                                      & fu_latest_time(missing_obstimes))
      RETURN
    END IF


    !-------------------------------------------------------------
    !
    ! 4. Check quantity by quantity
    !    ---------------------------

    missing_obstimes =  req_obstimes

    time_loop: DO i = 1, SIZE(req_obstimes)
      IF (.NOT.defined(req_obstimes(i))) EXIT time_loop

      quantity_loop: DO j = 1, SIZE(quantities)

        IF (.NOT.fu_known_quantity(quantities(j))) EXIT quantity_loop

        IF (fu_multi_level_quantity(quantities(j))) THEN
          IF (.NOT.fu_field_in_sm(miniMarketPtr, met_src,&
                                & quantities(j),&
                                & req_obstimes(i),&
                                & level_missing,&
                                & .true.,&
                                & .false.)) THEN
            call msg_warning(fu_connect_strings('missing multi-level quantity:', &
                                              & fu_quantity_string(quantities(j))))
          END IF
        ENDIF
        IF (.NOT.fu_field_in_sm(miniMarketPtr, met_src,&
                              & quantities(j),&
                              & req_obstimes(i),&
                              & level_missing,& 
                              & .false.,&
                              & .false.)) THEN
          call msg_warning(fu_connect_strings('missing single-level quantity:', &
                                            & fu_quantity_string(quantities(j))))
          RETURN
        END IF

      END DO quantity_loop
    END DO time_loop

    !-------------------------------------------------------------
    !
    ! 5. No times missing, no quantities missing
    !    ---------------------------------------

    data_exists_already = .true.

  END SUBROUTINE check_sm_for_data_in_list


  ! ***************************************************************


  SUBROUTINE check_bm_for_times_in_list(miniMarketPtr, list,&
                                     & data_exists_already,&
                                     & missing_obstimes, &
                                     & observation_interval)
    ! Description:
    ! Cheks if times required by the shopping list
    ! are already stored in memory in the boundary market
    ! Returns a listwith missing obstimes.
    !
    ! Also modifies the timelimits so, that the list contains enough
    ! obstimes to cover whole timeperiod of given list. For this
    ! we add the closest obstime from start-time backwards and closest
    ! obstime from end-time forwards to the list. If either happens to be
    ! an obstime, then it is not modified.
    ! 
    ! Method:
    ! 1. If supermarket is empty, empty timelist is returned
    !
    ! 2. If observation times are found, that are required by the
    ! shopping list, but are not in supermarket obs.times, then
    ! a list is formed containing the missing times

    IMPLICIT NONE
    !
    ! Imported parameters with intent(in):
    type(mini_market_of_stacks), intent(inout) :: miniMarketPtr
    type(silja_shopping_list), intent(inout), target :: list
    ! Imported parameters with intent(out):
    LOGICAL, INTENT(out) :: data_exists_already
    TYPE(silja_time), DIMENSION(:), INTENT(out) :: missing_obstimes
    type(silja_interval), intent(in) :: observation_interval

    ! Local declarations:
    TYPE(silja_time), DIMENSION(max_times) :: req_obstimes, sm_times
    INTEGER :: i,j, miss_count, number_of_sm_times
    TYPE(silja_time) :: earliest_in_sm, latest_in_sm
    TYPE(silja_time) :: req_earliest, req_latest
    type(silja_shopping_list), pointer :: listPtr

    !-------------------------------------------------------------
    !
    ! 1. Basic checks
    !
    listPtr => list

    if(.not. defined(list))then
     call set_error('Boundary shopping list not defined','check_bm_for_times_in_list')
     return
    endif

     if(.NOT.(miniMarketPtr%stack_multiTime_initialized  .or. &
            & miniMarketPtr%stack_singleTime_initialized)) then
       call set_error('boundary market not initialized','check_bm_for_times_in_list')
       RETURN
    endif

    data_exists_already = .false.

    !-------------------------------------------------------------
    !
    ! 2. Modify timelimits.
    !    -----------------

    req_obstimes(:)=time_missing
    if(fu_list_time_indicator(list)==accept_same_month)then
      
      if(fu_day(fu_start_time(list))<16)then
        i = fu_mon(fu_start_time(list))-1
        j = fu_year(fu_start_time(list))
        if(i == 0)then
          i = 12
          j = j-1
        endif
      else
        i = fu_mon(fu_start_time(list))
        j = fu_year(fu_start_time(list))
      endif

      req_earliest = fu_set_time_utc(j, i, 16, 0, 0, 0.0)

      if(fu_day(fu_end_time(list))<16)then
        i = fu_mon(fu_start_time(list))
        j = fu_year(fu_start_time(list))
      else
        i = fu_mon(fu_start_time(list))+1
        j = fu_year(fu_start_time(list))
        if(i == 13)then
          i = 1
          j = j+1
        endif
      endif

      req_latest = fu_set_time_utc(j, i, 16, 0, 0, 0.0)

      CALL fix_shopping_time_boundaries(listPtr, req_earliest, req_latest)
      IF (error) RETURN     

      req_obstimes(1) = req_earliest
      req_obstimes(2) = req_latest

    else

      req_earliest = fu_closest_obstime(fu_start_time(list), backwards, observation_interval)
      req_latest = fu_closest_obstime(fu_end_time(list), forwards, observation_interval)
      IF (error) RETURN
  
      CALL fix_shopping_time_boundaries(listPtr, req_earliest, req_latest)
      IF (error) RETURN
    
      req_obstimes = fu_obstimes(req_earliest, req_latest,.false., observation_interval)

    endif

    missing_obstimes = req_obstimes ! set all missing first

    !-------------------------------------------------------------
    !
    ! 3. Check for times in list that are missing in sm
    !    -----------------------------------------------

    !-------------------------------------------------------------
    !
    ! 3.1. Times that there is already data:
    IF (miniMarketPtr%minimarket_empty)return
    CALL supermarket_times(miniMarketPtr, met_src_missing, sm_times, number_of_sm_times)
    IF (error) RETURN
    IF (number_of_sm_times == 0) RETURN

    latest_in_sm = fu_latest_time(sm_times)
    earliest_in_sm = fu_earliest_time(sm_times)


    !-------------------------------------------------------------
    !
    ! 3.2. Find which req_obstimes are not within sm-timelimts

    missing_obstimes = time_missing
    miss_count = 0

    DO i = 1, SIZE(req_obstimes)
      IF (.NOT.defined(req_obstimes(i))) EXIT 
      IF (.NOT.fu_between_times(req_obstimes(i), earliest_in_sm, &
                              & latest_in_sm, .true.)) THEN
        !Missing obstime found
        miss_count = miss_count + 1
        missing_obstimes(miss_count) = req_obstimes(i)
      END IF
    END DO

    !-------------------------------------------------------------
    !
    ! 3.3. If missing times found, set new shop-limits and exit

    IF (miss_count > 0) THEN
      CALL fix_shopping_time_boundaries(listPtr,&
                                      & fu_earliest_time(missing_obstimes),&
                                      & fu_latest_time(missing_obstimes))     
      RETURN
    END IF

    data_exists_already = .true.

  END SUBROUTINE check_bm_for_times_in_list



  ! ***************************************************************


  SUBROUTINE update_obstimes(miniMarketPtr)
    !
    ! Updates the supermarket's list of observation times for all met_srcs.
    !
    IMPLICIT NONE

    type(mini_market_of_stacks), intent(inout) :: miniMarketPtr

    ! Local declarations:
    INTEGER :: s,t
    TYPE(silja_stack), POINTER :: stack

    miniMarketPtr%obstimes(:,:) = time_missing

    DO s = 1, SIZE(miniMarketPtr%stacks_multiTime,1)
      DO t = 1, SIZE(miniMarketPtr%stacks_multiTime,2)
        stack => miniMarketPtr%stacks_multiTime(s,t)%ptr
        IF (.NOT.defined(stack)) EXIT
        miniMarketPtr%obstimes(s,t) = fu_valid_time(stack)
      END DO
    END DO

  END SUBROUTINE update_obstimes



!  ! ***************************************************************

!
!  SUBROUTINE update_met_srcs()
!
!    !
!    ! Updates the supermarket's list of met_srcs - actually puts it into 
!    ! the wdr object. A trick is that not every stack is necessarily a single-
!    ! met_src object. It may happen that one stack has data from several
!    ! meteo sources. This makes things very complicated.
!    !
!    IMPLICIT NONE
!
!    ! Local declarations:
!    INTEGER :: s
!    TYPE(silja_stack), POINTER :: stack
!
!    DO s = 1, SIZE(stackvector,1)
!      stack => stackvector(s, 1)
!      IF (.NOT.defined(stack)) EXIT
!      stack%met_src_idx = s
!      call set_met_src(stack%wdr_ptr, fu_met_src(stack), s)
!    END DO
!
!  END SUBROUTINE update_met_srcs



  ! ***************************************************************


  FUNCTION fu_closest_sm_met_src_time(miniMarketPtr, met_src, time, direction, ifMandatory) &
                        & result(stack)
    !
    ! Returns a pointer to the stack of the correct met_src and
    ! closest observation time found in supermarket. Used in retrieving
    ! data, not storing, so only existing, defined stacks are taken
    ! into consideration.
    ! 
    IMPLICIT NONE
    !
    ! Return value 
    TYPE(silja_stack), POINTER :: stack

    ! Imported parameters with intent(in):
    type(mini_market_of_stacks), intent(in) :: miniMarketPtr
    TYPE(silja_time), INTENT(in) :: time
    INTEGER, INTENT(in) :: direction
    type(meteo_data_source), intent(in) :: met_src
    logical, intent(in), optional :: ifMandatory

    ! Local declarations:
    INTEGER :: s, t, i
    INTEGER, DIMENSION(max_met_srcs) :: hit
    logical :: ifTimesDefined

    !-------------------------------------------------------------
    !
    ! If time_missing => single-time stack 
    !
    stack => null()
    IF(time == time_missing)THEN
      if(met_src == met_src_missing) then ! Accept all sources
        s = 1
      else
        s = fu_met_src_storage_index(miniMarketPtr, single_time_stack_flag, met_src)
      endif

      IF (error) return
      stack => miniMarketPtr%stacks_singleTime(s)%ptr
      RETURN
    END IF


    !-------------------------------------------------------------
    !
    ! Multi-time stack
    !
    if(met_src == met_src_missing) then ! Accept all sources
      s = 1
    else
      s = fu_met_src_storage_index(miniMarketPtr, multi_time_stack_flag, met_src)
    endif
    if(present(ifMandatory))then
      !
      ! ifMandatory says whether we care whether the time is present or not. 
      ! If it is absent and not mandatory, we shall just return empty stack
      !
      if(.not. ifMandatory)then
        ifTimesDefined = .false.
        do t = 1, size(miniMarketPtr%obstimes,2)
          if(defined(miniMarketPtr%obstimes(s,t)))then
            ifTimesDefined = .true.
            exit
          endif
        end do
        if(.not. ifTimesDefined) return
      endif 
      t = fu_closest_time(time, miniMarketPtr%obstimes(s,:), direction, .not. ifMandatory) ! requires silence
    else
      t = fu_closest_time(time, miniMarketPtr%obstimes(s,:), direction)
    endif

    IF (error.or.(t==0)) return

    stack => miniMarketPtr%stacks_multiTime(s,t)%ptr

  END FUNCTION fu_closest_sm_met_src_time


  !****************************************************************

  function fu_stack_by_indices(miniMarket, indSrc, indTime) result(stack)
    !
    ! Just returns the pointer to stack
    !
    implicit none

    type(mini_market_of_stacks), intent(in) :: miniMarket
    integer, intent(in) :: indSrc
    integer, intent(in), optional :: indTime

    type(silja_stack), pointer :: stack

    nullify(stack)
    !
    ! Depending whether the time index is given, return either the single- or multiTime stack
    !
    if(present(indTime))then
      !
      ! Time dimension is present, return the stack that is dynamic
      !
      if(miniMarket%stack_multiTime_initialized)then
        if(indTime > 0 .and. indTime <= size(miniMarket%stacks_multiTime,2))then
          if(indSrc > 0 .and. indSrc <= size(miniMarket%stacks_multiTime,1))then
            stack => miniMarket%stacks_multiTime(indSrc,indTime)%ptr
          else
            call msg('Source dimension is strange:',indSrc)
            call set_error('Source dimension is strange for miniMarket:' + miniMarket%name, &
                         & 'fu_stack_by_indices')
            nullify(stack)
          endif  ! indSrc OK
          
        else  ! if time index is OK
          call msg('Strange requested time index:', indTime)
          call set_error('Strange requested time index', 'fu_stack_by_indices')
          nullify(stack)
        endif

      else  ! if multiTime stack exists
        call set_error('Time dimension is requested but multitime stack does not exist in miniMarket:' + &
                     & miniMarket%name,'fu_stack_by_indices')
        nullify(stack)
      endif   ! if multiTime stack exists

    else
      !
      ! Absent time dimension, return the requested static stack
      !
      if(miniMarket%stack_singleTime_initialized)then
        if(indSrc > 0 .and. indSrc <= size(miniMarket%stacks_singleTime))then
          stack => miniMarket%stacks_singleTime(indSrc)%ptr
         else
            call msg('Source dimension is strange:',indSrc)
            call set_error('Source dimension is strange for miniMarket:' + miniMarket%name , &
                         & 'fu_stack_by_indices')
            nullify(stack)
          endif  ! indSrc OK
      else  ! if multiTime stack exists
        call set_error('Time dimension is not requested but no single-time stack in miniMarket:' + &
                     & miniMarket%name,'fu_stack_by_indices')
        nullify(stack)
      endif   ! if multiTime stack exists

    endif  ! single- or multi-time

  end function fu_stack_by_indices


  !*************************************************************************************

  function fu_stack_by_values(miniMarket, met_src, time) result(stack)
    !
    ! Just returns the pointer to dispersion stack
    !
    implicit none

    ! Imported parameters
    type(mini_market_of_stacks), intent(in) :: miniMarket
    type(meteo_data_source), INTENT(in) :: met_src
    type(silja_time), intent(in), optional :: time

    ! Return value
    type(silja_stack), pointer :: stack
  
    ! Local variables
    integer :: s, t, stack_type

    if(present(time))then
      stack_type= multi_time_stack_flag
    else
      stack_type= single_time_stack_flag
    endif

    s = fu_met_src_storage_index(miniMarket, stack_type, met_src)
    if(error)return
    
    if(stack_type == multi_time_stack_flag)then
      DO t = 1, SIZE(miniMarket%stacks_multiTime, 2)
        IF (.NOT.defined(fu_valid_time(miniMarket%stacks_multiTime(s, t)%ptr))) cycle
        IF (time == fu_valid_time(miniMarket%stacks_multiTime(s, t)%ptr))then
          stack => fu_stack_by_indices(miniMarket, s, t)
          return
        endif
      END DO
    else
      stack => fu_stack_by_indices(miniMarket, s, 1)
      return
    endif

    nullify(stack)

  end function fu_stack_by_values
  
  
  !*******************************************************************************************
  
  subroutine check_supermarket_fields_range(mm)
    ! Goes stack by stack calling the next-level range check for the stacks.
    implicit none
    type(mini_market_of_stacks), intent(in) :: mm
    integer :: iTmp, jTmp
    call msg('Checking the validity of ranges in minimarket:' + mm%name)

    do iTmp = 1, size(mm%stacks_singleTime)
      call check_stack_fields_ranges(mm%stacks_singleTime(iTmp)%ptr)
    end do
    DO jTmp = 1, SIZE(mm%stacks_multiTime, 2)
      DO iTmp = 1, SIZE(mm%stacks_multiTime, 1)
        call check_stack_fields_ranges(mm%stacks_multiTime(iTmp,jTmp)%ptr)
      end do
    end do
    
  end subroutine check_supermarket_fields_range
  
  
  
  function fu_minimarket_name(mm) result(chNm)
    implicit none
    type(mini_market_of_stacks), intent(in) :: mm
    character(len=20) :: chNm
    chNm = mm%name
  end function fu_minimarket_name

  ! ***************************************************************

  ! ***************************************************************
  ! 
  !
  !                 TESTS
  !
  !
  ! ***************************************************************
  ! ***************************************************************

  SUBROUTINE print_supermarket_contents(miniMarketPtr)

    ! Description:
    ! Print info about the contents of the meteorological fields
    ! in memory.
    !
    IMPLICIT NONE
    type(mini_market_of_stacks), intent(in) :: miniMarketPtr
    ! Local declarations:
    INTEGER :: s,t
    TYPE(silja_stack), POINTER :: stack

    IF (error) RETURN

    IF (.NOT. (miniMarketPtr%stack_multiTime_initialized .or. miniMarketPtr%stack_singleTime_initialized)) THEN
      call msg(' Supermarket untouched')
      RETURN
    END IF

    IF (miniMarketPtr%minimarket_empty ) THEN
      call msg('Supermarket empty')
      RETURN
    END IF

    if (associated(miniMarketPtr%stacks_multiTime)) then
            DO s = 1, SIZE(miniMarketPtr%stacks_multiTime,1)
              call msg(' Here is source number ', s, SIZE(miniMarketPtr%stacks_multiTime,1))
              DO t = 1, SIZE(miniMarketPtr%stacks_multiTime,2)
                call msg(' Here is time slot number ', t, SIZE(miniMarketPtr%stacks_multiTime,2))
                stack => miniMarketPtr%stacks_multiTime(s,t)%ptr
                IF (.NOT.defined(stack)) then
                  call msg('Undefined stack')
                  cycle
                endif
                IF (fu_stack_empty(stack)) then
                  call msg('Empty stack')
                  CYCLE
                endif
                CALL report(stack)
              END DO
            END DO
    endif

    if (associated(miniMarketPtr%stacks_singleTime)) then
      DO s = 1, SIZE(miniMarketPtr%stacks_singleTime,1)
        call msg(' Here is source number in single-time stacks: ', s)
        stack => miniMarketPtr%stacks_singleTime(s)%ptr
        IF (fu_stack_empty(stack)) then
          call msg('Permanent stack empty')
        else
          CALL report(stack)
        endif
      enddo
    endif
  
  END SUBROUTINE print_supermarket_contents



  ! ***************************************************************


  SUBROUTINE supermarket_test_fill(miniMarketPtr) !post_processed)

    ! Description:
    ! Fills supermarket with hirlam data
    ! for few hours around wallclock time.
    ! Used for test purposes.
    ! 
    IMPLICIT NONE

    ! Imported parameters with intent IN:
!    LOGICAL, INTENT(in) :: post_processed
    type(mini_market_of_stacks), intent(inout) :: miniMarketPtr
    ! Local declarations:
    TYPE(silja_shopping_list) :: list
    TYPE(silja_time), DIMENSION(max_times) :: times
    INTEGER :: n, i, j, nt
    type(meteo_data_source) :: met_src
    type(grads_template), dimension(1) :: template
    integer, dimension(1) :: iFFormat
    type(wdr_ptr), dimension(:), pointer :: wdrar
    allocate(wdrar(1))
    met_src = met_src_missing
    wdrar(1)%ptr = wdr_missing
    CALL initialize_mini_market(miniMarketPtr, &
                               & 'test', &
                               & 1, & ! one source
                               & 50,& ! we assume: max 1 timenode per file
                               & 800,& ! for each timenode
                               & 100,& ! for each timenode
                               & 100,& ! for each timenode
                               & 50,& ! for each timenode
                               & .true.,&
                               & wdrar, & 
                               &.true.,&
                               & .true.)
    IF (error) RETURN

    list = fu_set_shopping_list (met_src_missing, &
                               & (/temperature_flag,&
                                 & temperature_2m_flag,&
                                 & dew_point_temp_2m_flag,&
                                 & geopotential_flag, &
                                 & u_flag, &
                                 & v_flag, &
                                 & omega_flag, &
                                 & large_scale_accum_rain_flag,&
                                 & convective_accum_rain_flag,&
                                 & large_scale_rain_int_flag,&
                                 & convective_rain_int_flag,&
                                 & layer_thickness_flag, &
                                 & relative_humidity_flag,&
                                 & specific_humidity_flag,&
                                 & ground_pressure_flag,&
                                 & msl_pressure_flag,&
                                 & cloud_water_flag,&
                                 & total_precipitation_rate_flag /),&
                                 & time_missing,&
                                 & time_missing,&
                                 & level_missing,&
                                 & level_missing)
    IF (error) RETURN

    CALL arrange_supermarket(miniMarketPtr)

    call msg('data now for following times:')
    CALL supermarket_times(miniMarketPtr, met_src, times, nt)
    DO i = 1, nt
      call msg(fu_str(times(i)))
    END DO
    deallocate(wdrar)
  END SUBROUTINE supermarket_test_fill


!************************************************************************

  logical function minimarket_initialized(miniMarketPtr, stacktype)
    implicit none
    type(mini_market_of_stacks), intent(in) :: miniMarketPtr
    integer:: stacktype

    select case(stacktype)
      case(multi_time_stack_flag)
        minimarket_initialized = miniMarketPtr%stack_multiTime_initialized
      case(single_time_stack_flag)
        minimarket_initialized = miniMarketPtr%stack_singleTime_initialized
      case default
        call set_error('Strange stack type','minimarket_initialized')
        return
    end select

  end function minimarket_initialized

!************************************************************************

  subroutine check_grid(miniMarketPtr, i, j)
    implicit none
    type(mini_market_of_stacks), pointer :: miniMarketPtr
    integer::i,j
    call msg('')
    call msg('############### REPORTING GRID OF STACK ################')
    call report(fu_storage_grid(fu_wdr(miniMarketPtr%stacks_multiTime(i, j)%ptr)))
    call msg('')
  end subroutine check_grid

!************************************************************************

  subroutine init_singletime_fields(MarketPtr, pGrid, vertical, qList)
      !
      ! Adds fields to single-time stack of supermarket, fills them with real_missing
      ! sets valid_times to past
      !
      implicit none
      type(mini_market_of_stacks), pointer :: MarketPtr
      type(silja_grid),  intent(in)  :: pGrid
      type(silam_vertical),   intent(in) :: vertical
      TYPE(silja_shopping_list), intent(in) :: qList

      type(silja_field_id) :: idTmp 
      real, dimension(:), pointer :: vals
      logical :: iffound
      type(silja_grid) :: fieldGrid
      character(len=*), parameter :: sub_name = 'init_singletime_fields'
      integer :: fs_grid, iLev, shopQ, iVarLst, ind_grid
     
      do iVarLst = 1, fu_nbr_of_quantities(qList)
        shopQ = fu_quantity(qList,iVarLst)
        if (shopQ < 0) exit
        iffound = fu_field_in_sm(MarketPtr, & ! Already done before?
                              & met_src_missing, & ! metSrcsMultitime(iMetSrc), &
                              & shopQ, &
                              & time_missing, &   ! now,&
                              & level_missing, &
                              & .false., & ! 2d field
                              & permanent=.true.)
        if(error)then
          call msg_warning('Something went wrong with detection of already-done fields but continue')
          call unset_error(sub_name)
        endif

        if (iffound) then ! Already done
           if (fu_realtime_quantity(shopQ)) then
              call msg_warning('Realtime field: ' + fu_quantity_string(shopQ)+&
                      & ' already initialized!' )
           endif
        else ! Create the field
          if (.not. fu_realtime_quantity(shopQ)) then
            call msg("Not realtime quantity:"+ fu_quantity_string(shopQ))
            cycle
          endif
          call msg("initialising instant field:"+fu_quantity_string(shopQ))
          fieldGrid = pGrid
          if (fu_if_staggerX(shopQ)) fieldGrid = fu_staggered_grid("U_stag", pGrid, .true., .false.)
          if (fu_if_staggerY(shopQ)) fieldGrid = fu_staggered_grid("V_stag", pGrid, .false., .true.)
          
          fs_grid = fu_number_of_gridpoints(fieldGrid)
!call msg('fs grid', fs_grid)
          vals => fu_work_array(fs_grid)
          vals(1:fs_grid) = real_missing
          if (fu_multi_level_quantity(shopQ)) then !3D
!call msg('Creating empty 3D field:' + fu_quantity_string(shopQ))
            idtmp = fu_set_field_id(met_src_missing, &                                            
                                  & shopQ, &                     
                                  & really_far_in_past, &                     
                                  & zero_interval, &                                    
                                  & fieldGrid, &
                                  & surface_level) !some defined level
            do iLev=1,fu_NbrOfLevels(vertical)
              call set_level(idtmp, fu_level(vertical, iLev))
              call dq_store_2d(MarketPtr, idtmp, vals, &
                             & single_time_stack_flag, ifRandomise = .false., &
                             & iUpdateType = overwrite_field, storage_grid = fieldGrid)
            enddo
          else !2D
!call msg('Creating empty 2D field:' + fu_quantity_string(shopQ))
            idtmp = fu_set_field_id(met_src_missing, &                                            
                                  & shopQ, &                     
                                  & really_far_in_past, &                     
                                  & zero_interval, &                                    
                                  & fieldGrid, surface_level)
            call dq_store_2d(MarketPtr, idtmp, vals, &
                           & single_time_stack_flag, ifRandomise=.false., &
                           & iUpdateType = overwrite_field, storage_grid = fieldGrid)
          endif !2D/3D
!call report(idTmp)
          call free_work_array(vals)
          if(error)return
          
        endif  ! Already done
        
      enddo  !  Var list

      call arrange_supermarket(MarketPtr)
      if(error)return
  end subroutine init_singletime_fields


END MODULE supermarket_of_fields

