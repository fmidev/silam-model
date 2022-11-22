MODULE stacks 
  !
  ! A pile of met.data, lowest hierarchy level of silja's weatherserver
  !
  ! Description: 
  ! Contains the definition of data-structure silja_stack. In stack
  ! weather-data values are stored both in 2D and 3D-fields at the
  ! same time.
  !
  ! Original code: Mika Salonoja
  ! Author: Mikhail Sofiev, FMI. mikhail.sofiev@fmi.fi
  !
  ! All units: SI
  ! 
  ! Language: ANSI standard Fortran 90
  !
  ! 
  ! Modules used:
  USE fields_3d
  USE windfields_3d
!  USE input_data_rules  !nwpm_administrations

  IMPLICIT NONE

  ! The public functions and subroutines available in this module:
  PUBLIC init_stack
  public set_met_src
  PUBLIC put_field_to_stack
  PUBLIC put_field_to_new_stack
  public get_field_place_in_stack
  public find_receiving_fields_4_merge
  public copy_stack_grid_interpolation
  PUBLIC arrange_fields_in_stack
  public prepare_new_averaging_period
  public update_stack_fields
  PUBLIC find_field_from_stack
  PUBLIC find_wind_from_stack
  PUBLIC find_field_3d_from_stack
  PUBLIC find_wind_3d_from_stack
  public fu_id_in_permanent_stack
  PUBLIC defined
  public set_defined
  PUBLIC fu_size
  PUBLIC fu_valid_time
  PUBLIC fu_met_src
  public fu_name
  PUBLIC report
  public check_stack_fields_ranges
  PUBLIC fu_stack_full
  PUBLIC fu_stack_empty
  PUBLIC set_stack_empty ! make stack empty for refill
  PUBLIC set_stack_empty_free_memory ! deallocates all memory
  PUBLIC fu_one_valid_time_stack
  PUBLIC fu_single_met_src_stack
  PUBLIC fu_number_of_3d_fields
  PUBLIC fu_number_of_wind_3d_fields
  PUBLIC fu_number_of_wind_2d_fields
  PUBLIC fu_number_of_2d_fields
  PUBLIC fu_first_field_from_stack
  PUBLIC fu_first_3d_field_from_stack
  PUBLIC stack_quantities
  PUBLIC stack_3d_quantities
  PUBLIC stack_variables
  PUBLIC fu_field_from_stack_by_quantity
  PUBLIC fu_get_wind_3d_field
  PUBLIC fu_get_wind_2d_field
  PUBLIC fu_get_3d_field
  PUBLIC fu_get_2d_field
  public fu_wdr

  ! The private functions and subroutines not to be used elsewhere:
  PRIVATE fu_stack_defined
  PRIVATE fu_stack_valid_time
  PRIVATE fu_stacksize
  private fu_name_of_stack
  PRIVATE arrange_uvw_to_windfields
  PRIVATE arrange_uv_10m_to_windfields
  PRIVATE arrange_scalar_fields_to_3d
  PRIVATE arrange_windfields_to_3d
  PRIVATE find_field_from_stack_by_id
  PRIVATE find_field_from_stack_direct
  PRIVATE print_stack_report
  private set_stack_met_src
!  private merge_inst_and_aver_2_instant
!  private merge_cumul_2_instant
!  private merge_inst_2_cumulative
!  private merge_cumul_2_cumulative
  private fu_stack_wdr
  private set_stack_defined
  private fu_stack_met_src

  ! Generic names and operator-interfaces of some functions:
  INTERFACE defined
    MODULE PROCEDURE fu_stack_defined
  END INTERFACE

  interface set_defined
    module procedure set_stack_defined
  end interface

  interface set_met_src
    module procedure set_stack_met_src
  end interface

  INTERFACE report
    MODULE PROCEDURE print_stack_report
  END INTERFACE

  INTERFACE fu_size
    MODULE PROCEDURE fu_stacksize
  END INTERFACE

  INTERFACE fu_name
    MODULE PROCEDURE fu_name_of_stack
  END INTERFACE

  INTERFACE fu_valid_time
    MODULE PROCEDURE fu_stack_valid_time
  END INTERFACE
  
  INTERFACE fu_met_src
    MODULE PROCEDURE fu_stack_met_src
  END INTERFACE
  
  INTERFACE fu_wdr
    MODULE PROCEDURE fu_stack_wdr
  END INTERFACE
  
  INTERFACE find_field_from_stack
    MODULE PROCEDURE find_field_from_stack_by_id
    MODULE PROCEDURE find_field_from_stack_direct
  END INTERFACE

  ! Public types defined in this module:  
  TYPE silja_stack
    PRIVATE

    CHARACTER (LEN=clen) :: name ! propably the name of the met_src,
    ! which' data is stored in this stack.

    ! The data in stack: NOTE. fields and winds are a real data, while
    !     3d stuff is just a set of poin%ters helping to arrange the buch of 2d
    !     fields. So, putting a 2d field involves two steps: dropping all
    !     levels to separate 2d fields (data are really copied to the stack)
    !     and then arranging these 2d fields (levels) to a 3d field.
    TYPE(silja_fieldpointer), DIMENSION(:), POINTER :: fields
    TYPE(silja_windfieldpointer), DIMENSION(:), POINTER :: winds
    TYPE(silja_field_3d_pointer), DIMENSION(:), POINTER :: fields_3d
    TYPE(silja_windfield_3d_pointer), DIMENSION(:), POINTER :: winds_3d

    INTEGER :: storepointer ! the index of the last (newest) field,
           ! that has been filled with data. If zero, then empty stack. When
           ! storepoiner reaches the upper limit of stack, then stack is
           ! emptied and storing is started from 1.
    INTEGER :: windstorepointer ! like above
    INTEGER :: storepointer_3d
    INTEGER :: windstorepointer_3d

    LOGICAL :: one_valid_time ! if true, then data from only
                           ! one valid-time is allowed in stack
    TYPE(silja_time) :: valid_time ! if %one_valid_time
    LOGICAL :: single_met_src ! if true, then data from only one met_src
                              ! (one nwp) is allowed in stack
    type(silja_wdr), pointer :: wdr_ptr  ! Pointer to weather_data_rules
    type(meteo_data_source) :: mds
    LOGICAL :: full
    TYPE(silja_logical) :: defined != silja_false
  END TYPE silja_stack
    
  type silam_stack_ptr
    type(silja_stack), pointer :: ptr
  end type silam_stack_ptr




CONTAINS 
  
  ! ***************************************************************


  SUBROUTINE init_stack(number_of_fields,&
                      & number_of_windfields,&
                      & name_of_stack,&
                      & one_valid_time, &
                      & single_met_src, &
                      & number_of_3d_fields,& 
                      & number_of_3d_windfields, &
                      & wdr, & 
                      & stack)
 
    ! Initializes a stack, that contains all the data from one
    ! met_src (eg. one model).
    !
    ! If this routine is calles for an already initialized stack
    ! (propably with data in it) this routine makes it empty! So be
    ! careful. If you wish to make an existing stack larger, use the
    ! routine make_stack_larger.
    !
    IMPLICIT NONE
    !
    ! Imported parameters with intent(in):
    INTEGER, INTENT(in) :: number_of_fields
    INTEGER, INTENT(in) :: number_of_windfields
    LOGICAL, INTENT(in) :: one_valid_time ! if true, then data from only
    ! one valid-time is allowed in stack. The first field inserted
    ! defined the time.
    LOGICAL, INTENT(in) :: single_met_src ! if true, then data from only
    ! one met_src (producer) is allowed in stack. The first field inserted
    ! defined the met_src.
    CHARACTER (LEN=*), INTENT(in) :: name_of_stack
    INTEGER, INTENT(in) :: number_of_3d_fields
    INTEGER, INTENT(in) :: number_of_3d_windfields
    type(silja_wdr), pointer :: wdr

    ! Imported parameters with intent(inout):
    TYPE(silja_stack), INTENT(inout), TARGET :: stack

    ! Local parameters:
    INTEGER :: i ! loop
    INTEGER :: status ! of allocation
    TYPE(silja_stack), POINTER ::sp

    !---------------------------------------------------------------
    !
    !  1. Check if not a new stack.
    !
    sp => stack
    IF (stack%defined == silja_true) THEN
      CALL set_error('cannot initialize stack again','init_stack')
      CALL report(sp)
      CALL set_error('cannot initialize stack again','init_stack')
      RETURN
    END IF

    !---------------------------------------------------------------
    !
    !  2. Allocate space for the vectors containing the fields.
    !     ----------------------------------------------------

    ALLOCATE (stack%fields(number_of_fields), stat = status)

    IF (status /= 0) THEN
      call msg('Allocation status: ', status)
      CALL set_error('Cannot ALLOCATE space for fields',' init_stack ')
      RETURN
    END IF
    IF(.not. associated(sp%fields)) THEN
      call msg('Fields are not associated but allocation status =0')
      CALL set_error('Cannot ALLOCATE space for fields',' init_stack ')
      RETURN
    END IF

    DO i = 1, SIZE(stack%fields)
      allocate(stack%fields(i)%fp, stat=status)
      if (status /= 0) then
        call set_error('Failed to allocate fields',' init_stack ')
        return
      end if
        !      CALL set_field_empty(stack%fields(i), .false.)
      CALL set_new_field(stack%fields(i)%fp)
    END DO


    !---------------------------------------------------------------
    !
    !  3. Allocate space for the vectors containing the windfields.
    !     --------------------------------------------------------

    IF (number_of_windfields > 0) THEN
      ALLOCATE (stack%winds(number_of_windfields), stat = status)
      
      IF (status /= 0) THEN
        PRINT *, 'Allocation status: ', status
        CALL set_error('Cannot ALLOCATE space for windfields',' init_stack ')
        RETURN
      END IF
      
      DO i = 1, SIZE(stack%winds)
        allocate(stack%winds(i)%fp, stat=status)
        if (status /= 0) then
          call set_error('Failed to allocate winds',' init_stack ')
          return
        end if
        CALL set_windfield_empty(stack%winds(i)%fp)
      END DO
    ELSE
      NULLIFY(stack%winds)
    END IF


    !---------------------------------------------------------------
    !
    !  4. Allocate space for the vectors containing the 3d-fields.
    !     ----------------------------------------------------

    IF (number_of_3d_fields > 0) THEN

      ALLOCATE (stack%fields_3d(number_of_3d_fields), stat = status)

      IF (status /= 0) THEN
        PRINT *, 'Allocation status: ', status
        CALL set_error('Cannot ALLOCATE space for 3d-fields','init_stack')
        RETURN
      END IF

      DO i = 1, SIZE(stack%fields_3d)
        allocate(stack%fields_3d(i)%fp, stat=status)
        if (status /= 0) then
          call set_error('Failed to allocate 3D fields',' init_stack ')
          return
        end if
        CALL set_3d_field_empty(stack%fields_3d(i)%fp)
      END DO

    ELSE
      NULLIFY(stack%fields_3d)
    END IF


    !---------------------------------------------------------------
    !
    !  5. Allocate space for the vectors containing the 3d-windfields.
    !     ----------------------------------------------------------

    IF (number_of_3d_windfields > 0) THEN

      ALLOCATE (stack%winds_3d(number_of_3d_windfields), stat=status)

      IF (status /= 0) THEN
        PRINT *, 'Allocation status: ', status
        CALL set_error ( 'Cannot ALLOCATE space for 3d-windfields',' init_stack ')
        RETURN
      END IF

      DO i = 1, SIZE(stack%winds_3d)
        allocate(stack%winds_3d(i)%fp, stat=status)
        if (status /= 0) then
          call set_error('Failed to allocate 3D windfields',' init_stack ')
          return
        end if
        CALL set_3d_windfield_empty(stack%winds_3d(i)%fp)
      END DO

    ELSE
      NULLIFY(stack%winds_3d)
    END IF


    !---------------------------------------------------------------
    !
    !  6. set the rest.
    !     -------------

    stack%storepointer = 0
    stack%windstorepointer = 0
    stack%storepointer_3d = 0
    stack%windstorepointer_3d = 0
    stack%one_valid_time = one_valid_time
    stack%single_met_src = single_met_src
    stack%wdr_ptr => wdr
    if(fu_NbrOfMetSrcs(wdr) == 1) then
      stack%mds = fu_mds(wdr,1)
    else
      stack%mds = met_src_missing
    end if
    stack%full = .false.
    stack%name = TRIM(ADJUSTL(name_of_stack))
    stack%valid_time = time_missing ! if single-time, then the first
                           ! field sets valid time for the whole stack


    !---------------------------------------------------------------
    !
    !  7. Initialization ok.
    !
    stack%defined = silja_true

!    call msg_test('Init_stack: finished')
!    sp => stack
!    call report(sp)

  END SUBROUTINE init_stack


  ! ***************************************************************

  

  subroutine set_stack_met_src(stackPtr, mdsPtr)
    !
    ! Sets the index referring to the wdr as a storage of all 
    ! data sources.
    !
    implicit none

    ! Imported pointers
    type(meteo_data_source), intent(in) :: mdsPtr
    type(silja_stack), pointer :: stackPtr

    if(.not.defined(stackPtr))then
      call set_error('Undefined stack given','set_met_src')
    else
      stackPtr%mds = mdsPtr
    end if


  end subroutine set_stack_met_src

  ! ***************************************************************


  SUBROUTINE put_field_to_stack (id,&
                               & grid_data,&
                               & stack, &
                               & ifMeteo_grid, &
                               & iUpdateType, &
                               & ifRandomise, &
                               & storage_grid, &
                               & storage_area, &
                               & factor, &
                               & ifForceMetSrcAcceptance, &
                               & iAccuracy, &
                               & fMissingValue_, &
                               & iOutside)
    
    ! Description:
    ! Puts a scalar 2D-field into a stack. Finds suitable place for
    ! it. Data is either
    ! 1) interpolated to the given storage (storage_grid defined), or
    ! 2) the given area is selected from the original field
    ! (storage_area defined), or
    ! 3) the field is stored as it is (both storage_grid and
    ! storage_area undefined).
    ! 
    ! If stack is full, the next fields are
    ! written on top of the oldest. 
    !
    IMPLICIT NONE
    !
    ! Imported parameters with intent IN:
    TYPE(silja_grid), INTENT(in) :: storage_grid
    TYPE(silam_area), INTENT(in) :: storage_area
    REAL, DIMENSION(:), INTENT(in), target, optional :: grid_data
    TYPE(silja_field_id), INTENT(in) :: id 
    LOGICAL, INTENT (in) :: ifMeteo_grid ! If grid must correspond to meteo_grid
    ! - If storage_area is defined then ifMeteo_grid determines the nx and ny of the
    !   selected grid - they must correspond to meteo values
    ! - If storage_grid is defined then ifMeteo_grid determines whether the Arakawa-grids
    !   must be interpolated to the meteo one.
    integer, intent(in) :: iUpdateType  ! what to do if the field already exists in stack
!    TYPE(silja_stack), INTENT(inout), TARGET :: stack
    TYPE(silja_stack), pointer :: stack
    real, intent(in), optional  :: factor 
    logical, intent(in) :: ifRandomise
    logical, intent(in), optional :: ifForceMetSrcAcceptance
    real, intent(in), optional :: fMissingValue_
    integer, intent(in) :: iOutside, iAccuracy

    ! Local declarations:
    INTEGER :: place, iCell
    TYPE(silja_field), POINTER :: old
    LOGICAL :: already_in_stack, selection_done, ifForceInputMetSrc
    TYPE(silja_grid) :: gridTmp, grid_new !Tmp variables are used to cover 
    TYPE(silam_area) :: areaTmp           ! presence/absence  of storage grid/area
    REAL, DIMENSION(:), pointer :: selected_grid_data
    TYPE(silja_field_id) :: id_selected, idTmp
    TYPE(silja_stack), POINTER :: stackpointer
    type(silja_wdr), pointer :: wdrPtr
    real, dimension(:), pointer :: dataTmpPtr
    real :: fMissingValue

    TYPE(silja_field), POINTER :: tmpFldPtr

    if(present(fMissingValue_))then
      fMissingValue = fMissingValue_
    else
      fMissingValue = real_missing
    endif

    stackpointer => stack

    !
    ! 1. Check input parameters.
    !
    IF (.NOT.defined(id)) RETURN ! nothing to put

    IF (.NOT.fu_full_id(id)) THEN
      CALL set_error('only fully identified fields to stack','put_field_to_stack')
      RETURN
    END IF

    if(associated(stack))then
      IF (.NOT.fu_true(stack%defined)) THEN
        CALL set_error('Cannot put field into an undefined stack','put_field_to_stack')
        RETURN
      END IF
    else
      call set_error('The forced stack if not associated','put_field_to_stack')
      return
    endif

    IF (stack%full) THEN
      CALL report(stackpointer, .false.)  ! (stackpointer, ifShort)
      CALL  set_error('cannot put data in full stack','put_field_to_stack')
      RETURN
    END IF

    !
    ! Setup the optional input parameters
    !
    if(defined(storage_grid).and.defined(storage_area))then
      call set_error('Both storage grid and area are defined','put_field_to_stack')
      return
    end if
    if(.not.  (defined(storage_grid).or.defined(storage_area)))then
      !'Both storage grid and area are undefined'
      ! use wdr stuff
        gridTmp = fu_storage_grid(stack%wdr_ptr)
        areaTmp = fu_storage_area(stack%wdr_ptr)
     else
      gridTmp = storage_grid
      areaTmp = storage_area
    end if

    !
    ! 2. Single time or not?
    !
    IF (stack%one_valid_time) THEN
      IF (defined(stack%valid_time)) THEN
        !
        ! If valid time already set, check that the new field fits:
        !
        IF (.NOT.fu_field_valid(id, stack%valid_time)) THEN
          CALL set_error('field is not valid at the stack valid time','put_field_to_stack')
          RETURN
        END IF
      ELSE
        !
        ! If no valid time yet (first field going in) then let this field dictate the valid time.
        ! Exception is the always-valid fields, which have undefined valid_time field.
        !
        if(defined(fu_validity_length(id)))then
          if(.not. (very_long_interval ==  fu_validity_length(id)))then
            stack%valid_time = fu_valid_time(id)
          endif
        else
          stack%valid_time = fu_valid_time(id)
        endif
      END IF
    END IF

    !
    ! 3. Single met_src? Careful - we still can force acceptance of any MDS into the stack
    !
    if(present(ifForceMetSrcAcceptance))then
      ifForceInputMetSrc = ifForceMetSrcAcceptance
    else
      ifForceInputMetSrc = .false.
    endif

    idTmp = id

    IF (stack%single_met_src) THEN
      if(ifForceInputMetSrc)then
        call set_met_src(idTmp, stack%mds)
      else
        IF (defined(stack%mds)) THEN
          !
          ! If met_src already set, check that the new field fits:
          !
          IF (.NOT.(fu_met_src(id) == stack%mds)) THEN
            call msg(fu_connect_strings('Field met src:',fu_name(fu_met_src(id))))
            call msg(fu_connect_strings('Stack met src:',fu_name(stack%mds)))
            CALL set_error('field met_src <> stack single met_src','put_field_to_stack')
            RETURN
          END IF
        else
          !
          ! If no met_src yet (first field going in), then let this
          ! field dictate the valid data source:
          !
          stack%mds = fu_mds(stack%wdr_ptr,fu_met_src(id)) ! Pointer must point to wdr element
        END IF  ! definde met_src
      END IF  ! if ifForceInputMetSrc
    else
      if(ifForceInputMetSrc)then
        call set_met_src(idTmp, met_src_missing)
      endif
    endif  ! if single met_src

    !
    ! 4. Check that it isn't already in stack. If it is not, find the right place for the new field
    !    If it is, check whether we should update the field.
    !    Mind the possible any-time relaxed field serach
    !
    if(iUpdateType == overwrite_field_anytime)then
      call find_field_from_stack_direct(fu_met_src(idTmp), &
                                        & fu_quantity(idTmp),&
                                        & time_missing,&
                                        & stackpointer, &
                                        & old, &
                                        & already_in_stack, &
                                        & fu_species(idTmp))
    else
      CALL find_field_from_stack_by_id(stackpointer, idTmp, old, already_in_stack)
    endif

    IF (already_in_stack) then
      !
      ! Field is in the stack. If to update it, grid settings must be forced the same
      !
      if(iUpdateType == create_field)then
        call set_error('Field already in stack:' + fu_quantity_short_string(fu_quantity(id)) + ',' + &
                     & fu_str(fu_species(id)), 'put_field_to_stack')
        RETURN
      elseif(iUpdateType == create_if_absent)then
        call msg_warning('Field already in stack:' + fu_quantity_short_string(fu_quantity(id)) + ',' + &
                       & fu_str(fu_species(id)), &
                       & 'put_field_to_stack')
        RETURN
      else
        areaTmp = area_missing
        gridTmp = fu_storage_grid(stack%wdr_ptr)  !fu_wdr(stackpointer))
      endif

    else
      !
      ! Field is new. Find the correct place for the new field.
      !
      place = stack%storepointer + 1

      IF (place > SIZE(stack%fields)) THEN
        call msg_warning('Stack:' + stack%name + '- now full')
        stack%full = .true.
        RETURN
      END IF
    endif

    if(.not. (defined(areaTmp) .or. defined(gridTmp)))then
      call set_error('Both grid and area are undefined in the stack','put_field_to_stack')
      return
    endif
!    call msg('Put-stack 8')
    
    if(present(grid_data))then
      !
      ! If the grid_data pointer is given, set the field with data and put it to stack
      !
      !  Maybe some memory reductions first.
      !
      ! - If areaTmp is defined then ifMeteo_grid determines the nx and ny of the
      !   selected grid - they must correspond to system values
      ! - If gridTmp is defined - ifMeteo_grid determines whether the Arakawa-grids
      !   must be interpolated to the meteo one.
      ! But here no analysis is done, the ifMeteo_grid is just passed further
      !
      ! ATTENTION.
      ! Here we cannot reserve the space for selected_grid_data since grid_new is unknown
      !
      ! Data selection is done in two steps: (i) get the new grid and allocate space for ther data
      ! (ii) do the reprojection itself
      !
      call grid_data_hor_select_new_grid(fu_grid(idTmp), &  ! grid_original
                                       & gridTmp, &        ! storage_grid
                                       & areaTmp, &       ! storage_area
                                       & ifMeteo_grid, &    ! ifAdjust to system_frid
                                       & grid_new)      ! new grid_new
      if(error)return
      selected_grid_data => fu_work_array(fu_number_of_gridpoints(grid_new))
      if(error)return
      
      CALL grid_data_horizontal_select(fu_grid(idTmp), &
                                     & grid_data,&
                                     & selection_done,&
                                     & grid_new,&
                                     & selected_grid_data, &
                                     & ifRandomise, &
                                     & iAccuracy, &
                                     & fu_regridding_method(fu_quantity(idTmp)), &
                                     & fMissingValue, &
                                     & iOutside)
      IF (error) RETURN

      if(fu_fails(fu_stdSilamGrid(grid_new),'Non-standard grid after data_select','put_field_to_stack'))then
        call msg('Initial requested grid:')
        call report(gridTmp)
        call msg('input-field grid:')
        call report(fu_grid(idTmp))
        call msg('Selected grid:')
        call report(grid_new)
        call msg('')
      endif
    
      !
      ! 7. Put it in!
      !
      IF (selection_done) THEN  ! memory-saving attempt succeeded

        id_selected = idTmp
        CALL set_grid(id_selected, grid_new)
        IF (error) RETURN

        if(present(factor))then
          do iCell = 1, fu_number_of_gridpoints(fu_grid(id_selected))
            selected_grid_data(iCell) = factor * selected_grid_data(iCell)
          enddo
        endif

        if(already_in_stack)then
          call update_field_data(id_selected, selected_grid_data, old, iUpdateType) ! Update existing field
          if(error)return
        else

          CALL set_field(id_selected,&  ! The data and id are copied to the new field
                       & selected_grid_data,&
                       & stack%fields(place)%fp, &
                       & .false.)       ! not necessarily new

        ! *** above is dangerous: if the field is new, an invalid deallocation might happen!
        ! 

          IF (error) RETURN
          stack%storepointer = place
        endif
        !
        ! Mind that: release of the array allocated in the _select sub above
        !
        call free_work_array(selected_grid_data)

      ELSE
        !
        ! Selection is not done. Pointer selected_grid_data is undefined, should be just redirected
        !
        selected_grid_data => grid_data

        if(present(factor))then
          do iCell = 1, fu_number_of_gridpoints(fu_grid(idTmp))
            selected_grid_data(iCell) = factor * selected_grid_data(iCell)
          enddo
        endif

        if(already_in_stack)then
          call update_field_data(idTmp, selected_grid_data, old, iUpdateType)   ! Update existing
          if(error)return
        else
          CALL set_field(idTmp, selected_grid_data, stack%fields(place)%fp, .false.)  ! Add new
          IF (error) RETURN
          stack%storepointer = place
        endif

      END IF  ! selection_done
      
    else
      !
      ! No grid_data is given, i.e. we want the empty void field. Just call the 
      ! set_field without any data pointer, it will understand
      !
      call set_field(idTmp, stack%fields(place)%fp, .false.)       ! not necessarily new
      if(error)return
      stack%storepointer = place
      
    endif  ! if present grid_data

  END SUBROUTINE put_field_to_stack


  ! ***************************************************************


  SUBROUTINE put_field_to_new_stack (id,&
                                   & grid_data,&
                                   & storage_grid,&
                                   & storage_area,&
                                   & newstack, &
                                   & ifMeteo_grid, &
                                   & ifRandomise, &
                                   & iAccuracy, method, & 
                                   & fMissingValue_)
    ! 
    ! Puts a scalar 2D-field into a stack. Finds suitable place for
    ! it. Data is either
    ! 1) interpolated to the given storage (storage_grid defined), or
    ! 2) the given area is selected from the original field
    ! (storage_area defined), or
    ! 3) the field is stored as it is (both storage_grid and
    ! storage_area undefined).
    ! - If storage_area is defined then ifMeteo_grid determines nx and ny of the
    !   selected grid - they must correspond to meteo values
    ! - If storage_grid is defined - ifMeteo_grid determines whether Arakawa-grids
    !   must be interpolated to the meteo one.
    ! 
    ! If stack is full, the next fields are
    ! written on top of the oldest. 
    ! 
    IMPLICIT NONE
    !
    ! Imported parameters with intent IN:
    TYPE(silja_grid), INTENT(in) :: storage_grid
    TYPE(silam_area), INTENT(in) :: storage_area
    REAL, DIMENSION(:), pointer :: grid_data
    TYPE(silja_field_id), INTENT(in) :: id 
    LOGICAL, INTENT (in) :: ifMeteo_grid, ifRandomise
    real, intent(in), optional :: fMissingValue_
    integer, intent(in) :: iAccuracy, method

    
    ! Imported parameters with intent INOUT or POINTER:
    TYPE(silja_stack), INTENT(inout), TARGET :: newstack

    ! Local declarations:
    INTEGER :: place
    TYPE(silja_field), POINTER :: old
    LOGICAL :: already_in_stack
    TYPE(silja_grid) :: gridTmp, grid_new !Tmp variables are used to cover 
    TYPE(silam_area) :: areaTmp           ! presence/absence  of storage grid/area
    REAL, DIMENSION(:), pointer :: selected_grid_data
    LOGICAL :: selection_done
    TYPE(silja_field_id) :: id_selected
    TYPE(silja_stack), POINTER :: stackpointer
    real :: fMissingValue


    !----------------------------------------
    !
    ! 1. Check input parameters.
    !
    if(present(fMissingValue_))then
      fMissingValue = fMissingValue_
    else
      fMissingValue = real_missing
    endif
!    PRINT*,'1'    
    IF (.NOT.defined(id)) RETURN ! nothing to put
    CALL report(id)
    IF (.NOT.fu_full_id(id)) THEN
      CALL set_error('only fully identified fields to stack',&
                   & 'put_field_to_new_stack')
      RETURN
    END IF
!!!$    PRINT*,'fu_true(newstack%defined)',fu_true(newstack%defined)
!!!$    IF (.NOT.fu_true(newstack%defined)) THEN
!!!$      CALL set_error('Cannot put field into an undefined stack',&
!!!$    & 'put_field_to_new_stack')
!!!$      RETURN
!!!$    END IF

    IF (newstack%full) THEN
      stackpointer => newstack
      CALL report(stackpointer, .true.)
      CALL  set_error('cannot put data in full stack','put_field_to_new_stack')
      RETURN
    END IF

    if(defined(storage_grid).and.defined(storage_area))then
      call set_error('Both storage grid and area are defined','put_field_to_new_stack')
      return
    end if
    gridTmp = storage_grid
    areaTmp = storage_area


    ! -----------------------------------------------
    !
    ! 2. Single time or not?
    !    -------------------
!    PRINT*,'2'
    IF (newstack%one_valid_time) THEN
      IF (defined(newstack%valid_time)) THEN
        ! If valid time already set, check that the new field is it
        ! too:
        IF (.NOT.(fu_valid_time(id)==newstack%valid_time)) THEN
          CALL set_error('field time not stack valid time','put_field_to_newstack')
          RETURN
        END IF
      ELSE
        ! If no valid time yet (first field going in), then let this
        ! field dictate the valid time:
        newstack%valid_time = fu_valid_time(id)
      END IF
    END IF


    ! -----------------------------------------------
    !
    ! 3. Single met_src?
    !    --------------
!    PRINT*,'3'
    IF (newstack%single_met_src) THEN
      IF (defined(newstack%mds)) THEN
        ! If met_src already set, check that the new field is it too:
        IF (.NOT.(fu_met_src(id) == newstack%mds)) THEN
          call msg(fu_connect_strings('Field met source:',fu_name(fu_met_src(id))))
          call msg(fu_connect_strings('Stack met source:',fu_name(newstack%mds)))
          CALL set_error('field met_src <> stack single met_src'&
                & ,'put_field_to_new_stack')
          RETURN
        END IF
      else
        ! If no met_src yet (first field going in), then let this
        ! field dictate the valid data source:
        newstack%mds = fu_mds(newstack%wdr_ptr,fu_met_src(id))
      end if
    END IF


    !----------------------------------------
    !
    ! 4. Check that it isn't already in stack.
    !    -------------------------------------
!    PRINT*,'4'   
    stackpointer => newstack
    CALL find_field_from_stack(stackpointer, id, old, already_in_stack)

    IF (already_in_stack) RETURN


    !----------------------------------------
    !
    ! 5. Find the correct place for the new field.
    !    ----------------------------------------
!    PRINT*,'5'    
    place = newstack%storepointer + 1
    
    IF (place > SIZE(newstack%fields)) THEN
      PRINT *, 'Newstack ', TRIM(newstack%name), ' now full.'
      newstack%full = .true.
      RETURN
    END IF

    selected_grid_data => fu_work_array(fu_number_of_gridpoints(fu_grid(id)))

    !-----------------------------------------------------
    !
    ! 6. Maybe some memory reductions first.
    !
    ! Data selection is done in two steps: (i) get the new grid and allocate space for ther data
    ! (ii) do the reprojection itself
    !
    call grid_data_hor_select_new_grid(fu_grid(id), &  ! grid_original
                                     & gridTmp, &      ! storage_grid
                                     & areaTmp, &      ! storage_area
                                     & ifMeteo_grid, & ! ifAdjust to system_frid
                                     & grid_new)       ! new grid_new
    if(error)return
    selected_grid_data => fu_work_array(fu_number_of_gridpoints(grid_new))
    if(error)return

    CALL grid_data_horizontal_select(fu_grid(id), &
                                   & grid_data,&
                                   & selection_done,&
                                   & grid_new,&
                                   & selected_grid_data, &
                                   & ifRandomise, &
                                   & iAccuracy, method, & 
                                   & fMissingValue, &
                                   & notAllowed)
    IF (error) RETURN    

    
    !-----------------------------------------------------
    !
    ! 7. Put it in!
    !    ---------
!    PRINT*,'7'
     IF (selection_done) THEN
  
      id_selected = id
      CALL set_grid(id_selected, grid_new)
      IF (.not. error) CALL set_field(id_selected,&
                   & selected_grid_data,&
                   & newstack%fields(place)%fp, .false.)

    ELSE
      CALL set_field(id, grid_data, newstack%fields(place)%fp, .false.)
    END IF

    call free_work_array(selected_grid_data)

    IF (error) RETURN

    newstack%storepointer = place

  END SUBROUTINE put_field_to_new_stack

  
  !*************************************************************************
  
  subroutine get_field_place_in_stack(pStack, id, pField, ifForceMetSrcAcceptance)
    !
    ! Creates a field storage place in the given stack.
    ! The field is assumed not to be present in the stack
    !
    implicit none
    
    ! Imported parameters
    type(silja_field_id), intent(in) :: id
    type(silja_stack), pointer :: pStack
    type(silja_field), pointer :: pField
    logical, intent(in) :: ifForceMetSrcAcceptance

    ! Local variables
    logical :: ifFound
    !
    ! First, a simple option: does it exist already?
    !
    CALL find_field_from_stack(pStack, id, pField, ifFound) ! available right away?
    if(ifFound)return

    !
    ! Have to create the place holder and then return the pointer to it.
    !
    ! Single-time stack or not?
    !
    IF (pStack%one_valid_time) THEN
      IF (defined(pStack%valid_time)) THEN
        !
        ! If valid time already set, check that the new field fits:
        !
        IF (.NOT.fu_field_valid(id, pStack%valid_time)) THEN
          CALL set_error('field is not valid at the stack valid time','get_field_place_in_stack')
          RETURN
        END IF
      ELSE
        !
        ! If no valid time yet (first field going in) then let this field dictate the valid time.
        ! Exception is the always-valid fields, which have undefined valid_time field.
        !
        if(defined(fu_validity_length(id)))then
          if(.not. (very_long_interval ==  fu_validity_length(id)))then
            pStack%valid_time = fu_valid_time(id)
          endif
        else
          pStack%valid_time = fu_valid_time(id)
        endif
      END IF
    END IF
    !
    ! Single met_src? Careful - we can force acceptance of any MDS into the stack
    !
    IF (pStack%single_met_src) THEN
      IF (defined(pStack%mds)) THEN
        IF (.NOT.(fu_met_src(id) == pStack%mds)) THEN
          if(.not. ifForceMetSrcAcceptance)then
            call msg('Field met src:' + fu_name(fu_met_src(id)))
            call msg('Stack met src:' + fu_name(pStack%mds))
            CALL set_error('field met_src <> stack single met_src','get_field_place_in_stack')
            RETURN
          endif  ! not force MDS
          pStack%mds = fu_mds(pStack%wdr_ptr,fu_met_src(id)) ! Pointer must point to wdr element
        endif  ! same MDS
      else
        pStack%mds = fu_mds(pStack%wdr_ptr,fu_met_src(id)) ! Pointer must point to wdr element
      endif  ! defined stack mds
    endif   ! single-mds stack
    !
    ! Have space for the field? 
    !
    pStack%storepointer = pStack%storepointer + 1
    IF (pStack%storepointer > SIZE(pStack%fields)) THEN
      call msg('');call msg('')
      call msg_warning('Stack full:' + pStack%name,'get_field_place_in_stack')
      call msg("Stack size ", SIZE(pStack%fields))
      call report(pStack)
      call set_error('Stack full:' + pStack%name,'get_field_place_in_stack')
      call msg('');call msg('')
      pStack%full = .true.
      RETURN
    endif
    !
    ! OK, let's finally set the field
    !
    CALL set_field(id, pStack%fields(pStack%storepointer)%fp, .false., ifResetVals=.False.)       ! not necessarily new
    pField => pStack%fields(pStack%storepointer)%fp

  end subroutine get_field_place_in_stack


  !*************************************************************************

  subroutine find_receiving_fields_4_merge(idIn, stack, targetId, &
                                         & found, found_internal, timeDirection, &
                                         & idOut, dataOut, idOut_internal, dataOut_internal)
    !
    ! Looks for the fields in the stack, with which the input Id can be merged.
    ! Place must be reserved in the stack and, possibly, half-filled during previous time steps
    ! So, we should find exact quantity, name, and level, but not time.
    ! The met_src is also a problem - in the input it can be e.g. HIRLAM, while 
    ! in the output it is, of course, SILAM
    ! If the targetId is averaged - we have to find the intermediate field too
    !
    implicit none

    ! Imported parameters 
    type(silja_field_id), intent(in) :: idIn, targetId
    real, dimension(:), pointer :: dataOut, dataOut_internal
    type(silja_stack), intent(in) :: stack
    type(silja_field_id), pointer :: idOut, idOut_internal
    logical, intent(out) :: found, found_internal
    integer, intent(out) :: timeDirection

    ! Local variables
    !
    type(silja_field_id), pointer :: idTmp
    type(silja_field), pointer :: field
    integer :: i, timeDir

    found = .false.
    found_internal = .not. (fu_field_kind(targetId) == averaged_flag .or. &
                          & fu_field_kind(targetId) == accumulated_flag)
    do i = 1, stack%storepointer

      field => stack%fields(i)%fp
      idTmp => fu_id(field)
      
      IF(fu_quantity(idTmp) /= fu_quantity(idIn)) CYCLE

      if(fu_substance_name(idTmp) /= fu_substance_name(idIn)) cycle

      if (.not. fu_mode(idTmp) == fu_mode(idIn))cycle


      if(.not. (fu_optical_wave_length(idIn) .eps. real_missing))then
        if(.not. (1.e6*fu_optical_wave_length(idTmp) .eps. 1.0e6*fu_optical_wave_length(idIn))) cycle
      endif

      !
      ! Basically, grids must be equal, but it may happen that they are 
      ! Arakawa-shifted. This may happen because we did not know this shift
      ! when the stack was initialised. So, during the first step we may need
      ! to adjust the grid to its Arakawa-shifted analog
      !
      if(.not. fu_grid(idTmp) == fu_grid(idIn))then
        if(fu_grids_arakawa_correspond(fu_grid(idTmp), fu_grid(idIn)))then
          !
          ! Arakawa correspondence allows up to a half-cell shift and different nx,ny.
          !
          !One-cell shifted grids are treated in dispersion_server::merge_vector_data_to_stack

          if(fu_number_of_gridpoints(fu_grid(idTmp)) == fu_number_of_gridpoints(fu_grid(idIn)))then
            call msg_warning('Adjusting grid to Arakawa shift','find_receiving_fields_4_merge')
            call set_grid(idTmp, fu_grid(idIn))  ! Adjust grid to Arakawa-shift
!          else
!            call msg_warning('NOT adjusting grid to Arakawa shift. different grids left in data. Trouble ahead','find_receiving_fields_4_merge')
          endif
        else
          cycle   ! Grids do not correspond
        endif
      endif

      ! In the stack there must be no missing-levels, BUT in some
      ! cases we do not know what level will appear for a 2D variable
      ! In this situation, the level is set to missing explicitly. Then,
      ! the first-coming field will overwrite the missing level with its own
      !
      if(fu_cmp_levs_eq(fu_level(idTmp), level_missing))then
        call set_level(idTmp, fu_level(idIn))
      else
        if(.not.fu_cmp_levs_eq(fu_level(idTmp), fu_level(idIn)))cycle
      endif
      !
      ! we found either intermediate or final field - let's store it
      !
      if(fu_met_src(idTmp) == silam_internal_src)then
        found_internal = .true.
        idOut_internal => idTmp
        dataOut_internal => fu_grid_data(field)
      else
        found = .true.
        idOut => idTmp
        dataOut => fu_grid_data(field)
      endif

      if(found_internal .and. found) exit ! Leave when both are found

    end do   ! search the corresponding field in the stack

    if(.not.found)then
      call msg('')
      call msg('Input field to merge: ')
      call report(idIn)
      call msg('')
      call msg('Output stack: ')
      call report(stack)
      call set_error('Cannot find field to merge in given stack','find_receiving_fields_4_merge')
      return
    endif

    if ( (fu_validity_length(idIn)) == zero_interval ) then
      ! only instant fields do need these hacks
      ! Though, i have no clue how to figure out time direction....  R.
      if(fu_valid_time(idOut) < fu_valid_time(targetId))then
        !
        ! Forward in time: target ID is later than existing ID
        ! New field must not be later than target one
        !
        timeDirection = forwards
        if(fu_valid_time(idIn) > fu_valid_time(targetId))then
          call msg('')
          call msg('FORWARD merging. Field that comes:')
          call report(idIn)
          call msg('Target id')
          call report(targetId)
          call set_error('New id valid time is later than target one', 'find_receiving_fields_4_merge')
          return
        endif
      else
        !
        ! Backward in time: target ID is earlier than existing ID
        ! New field must not be earlier than target one
        !
        timeDirection = backwards
        if(fu_valid_time(idIn) < fu_valid_time(targetId))then
          call msg('')
          call msg('INVERSE merging. Field that comes:')
          call report(idIn)
          call msg('Target id')
          call report(targetId)
          call set_error('New id valid time is later than target one', 'find_receiving_fields_4_merge')
          return
        endif
      endif  ! forward vs inverse merging
    endif !instant fields

  end subroutine find_receiving_fields_4_merge



  !**********************************************************************************

  subroutine copy_stack_grid_interpolation(stackFrom, stackTo, gridNew, ifCopyAll, iOutGrid, ifRandomise)
    !
    ! Copies all fields from stackFrom to stackTo, interpolating all fields
    ! to the requested grid gridNew. Attention: the destination stack is cleaned !
    !
    implicit none

    ! Imported parameters
    type(silja_stack), pointer :: stackFrom, stackTo
    type(silja_grid), intent(in) :: gridNew
    logical, intent(in) :: ifCopyAll ! If not - skip internal met_src
    logical, intent(in) :: ifRandomise ! Smoothing the reprojection
    integer, intent(in) :: iOutGrid  ! out of grid interpolation type
    
    ! Local variables
    integer :: iField, qTmp
    type(silja_field), pointer :: fieldPtr

    !
    ! Stupidity check, cleaning the space. Note that empty stack is OK - mass map can still handle
    ! the output
    !
    if(.not.defined(stackFrom))then
      call set_error('Undefined source stack','copy_stack_grid_interpolation')
      return
    endif
    if(.not.defined(stackTo))then
      call set_error('Undefined destination stack','copy_stack_grid_interpolation')
      return
    endif
    if(stackFrom%storepointer < 1)then
!      call set_error('Empty source stack','copy_stack_grid_interpolation')
      return
    endif
    
    call set_stack_empty(stackTo)

!    call msg('copy_stack_grid_interpolation reports new grid')
!    call report(gridNew)

    do iField = 1, stackFrom%storepointer
      
      fieldPtr => stackFrom%fields(iField)%fp
!      call msg(fu_connect_strings('Quantity:',fu_quantity_short_string(fu_quantity(fu_id(fieldPtr)))))
      if((.not.ifCopyAll) .and. fu_met_src(fieldPtr) == silam_internal_src)cycle

!      call msg('Copying...')
      qTmp = fu_quantity(fu_id(fieldPtr))
      CALL put_field_to_stack (fu_id(fieldPtr),& 
                             & fu_grid_data(fieldPtr),&       ! grid data
                             & stackTo, & ! stack to be filled
                             & .not. (fu_if_StaggerX(qTmp) .or. fu_if_StaggerY(qTmp)), &     !.false., & ! ifMeteo_grid adjustment
                             & create_if_absent, & ! ifUpdateAllowed 
                             & ifRandomise, &      ! smoothing the reprojection
                             & gridNew,&  ! storage grid 
                             & area_missing, & ! storage area
                             & iAccuracy = 5, & ! iAccuracy
                             & iOutside = iOutGrid)        ! out of grid value
      if(error)return
    enddo

  end subroutine copy_stack_grid_interpolation


  ! ***************************************************************


  SUBROUTINE arrange_fields_in_stack(stack, ifLookForPressure)
    
    ! Description:
    ! This routine does three things:
    ! 1. Arranges all wind-component fields (u,v,w) into windfields
    ! of type silja_windfield. No data is duplicated, only pointers
    ! are set.
    ! 2. Arranges scalar fields into 3d-scalar fields. A 3d-field
    ! contains all available fields for one quantity, met_src and
    ! valid time, but all levels. No data is duplicated, only pointers
    ! are set.
    ! 3. Arranges the earlier arranged windfields into 3d-windfields,
    ! just like scalar fields described above.
    !
    ! Language: ANSI Fortran 90
    ! 
    ! Author: Mika Salonoja, FMI
    ! 
    IMPLICIT NONE

    ! Imported parameters with intent INOUT or POINTER:
    TYPE(silja_stack), POINTER :: stack
    logical, intent (in) :: ifLookForPressure

    !
    ! 1. Check stack.
    !
    IF (.NOT.defined(stack)) THEN
      CALL set_error('cannot arrange undefined stack','arrange_fields_in_stack')
      RETURN
    END IF
    
    !
    ! 2. Make windfields out of u,v and w.
    !
    CALL arrange_uvw_to_windfields(stack)
    IF (error) RETURN
    
    CALL arrange_uv_10m_to_windfields(stack)
    IF (error) RETURN
    !
    ! 3. Arrange scalar fields to stacked 3d fields.
    !
    IF (ASSOCIATED(stack%fields_3d))then
      CALL arrange_scalar_fields_to_3d(stack, ifLookForPressure)
      IF (error) RETURN
    END IF
    !
    ! 4. Arrange windfields to stacked 3d windfields.
    !
    IF (ASSOCIATED(stack%winds_3d)) CALL arrange_windfields_to_3d(stack)
    
  END SUBROUTINE arrange_fields_in_stack


  ! ***************************************************************


  SUBROUTINE set_stack_empty(stack)
    
    ! Description:
    ! Sets stack empty. All data in stack is lost. Identifications
    ! are set empty, but the actual memory-areas for fields are kept.
    ! After this routine stack is still defined, but empty.
    !
    ! Language: ANSI Fortran 90
    !
    ! Author: Mika Salonoja, FMI
    ! 
    IMPLICIT NONE

    ! Imported parameters with intent INOUT or POINTER:
    TYPE(silja_stack), INTENT(inout) :: stack

    ! Local declarations:
    INTEGER :: status, i
    TYPE(silja_field), POINTER :: field

    IF (ASSOCIATED(stack%fields)) THEN
      DO i = 1, SIZE(stack%fields)
        CALL set_field_empty(stack%fields(i)%fp, .false.)
        IF (error) RETURN
      END DO
    END IF

    IF (ASSOCIATED(stack%winds)) THEN
      DO i = 1, SIZE(stack%winds)
        CALL set_windfield_empty(stack%winds(i)%fp)
        IF (error) RETURN
      END DO
    END IF

    IF (ASSOCIATED(stack%fields_3d)) THEN
      DO i = 1, SIZE(stack%fields_3d)
        CALL set_3d_field_empty(stack%fields_3d(i)%fp)
        IF (error) RETURN
      END DO
    END IF

    IF (ASSOCIATED(stack%winds_3d)) THEN
      DO i = 1, SIZE(stack%winds_3d)
        CALL set_3d_windfield_empty(stack%winds_3d(i)%fp)
        IF (error) RETURN
      END DO
    END IF

    stack%storepointer = 0
    stack%windstorepointer = 0
    stack%storepointer_3d = 0
    stack%windstorepointer_3d = 0
    stack%name  = ' '
    stack%valid_time = time_missing
    stack%mds = met_src_missing
    stack%full = .false.

  END SUBROUTINE set_stack_empty
    

  ! ***************************************************************


  SUBROUTINE set_stack_empty_free_memory(stack)

    ! Description:
    ! Sets stack empty. All data in stack is lost. All memory is
    ! deallocated. After this routine stack is undefined.
    ! Correction: after deallocation, the pointers must be nullified !
    !
    ! Language: ANSI Fortran 90
    !
    ! Original code: Mika Salonoja
    ! Author: Mikhail Sofiev, FMI
    ! 
    IMPLICIT NONE

    ! Imported parameters with intent INOUT or POINTER:
    TYPE(silja_stack), INTENT(inout), TARGET :: stack

    ! Local declarations:
    INTEGER :: status, i
    TYPE(silja_field), POINTER :: field
    TYPE(silja_stack), POINTER :: sp

    IF (ASSOCIATED(stack%winds)) THEN
      DEALLOCATE(stack%winds, stat= status)
      IF (status /= 0) THEN
        sp => stack
        CALL report(sp)
        CALL set_error('cannot deallocate windfields','set_stack_empty_free_memory')
        RETURN    
      END IF
      NULLIFY(stack%winds)
    END IF

    IF (ASSOCIATED(stack%fields_3d)) THEN
      DEALLOCATE(stack%fields_3d, stat= status)
      IF (status /= 0) THEN
        CALL set_error('cannot deallocate 3d-fields','set_stack_empty_free_memory')
        RETURN    
      END IF
      NULLIFY(stack%fields_3d)
    END IF

    IF (ASSOCIATED(stack%winds_3d)) THEN
      DEALLOCATE(stack%winds_3d, stat= status)
      IF (status /= 0) THEN
        CALL set_error('cannot deallocate 3d-windfields','set_stack_empty_free_memory')
        RETURN    
      END IF
      NULLIFY(stack%winds_3d)
    END IF

    IF (ASSOCIATED(stack%fields)) THEN
      DO i = 1, SIZE(stack%fields)
        field => stack%fields(i)%fp
        IF (.NOT.defined(field)) EXIT
        CALL set_field_empty(stack%fields(i)%fp, .true.)
        IF (error) RETURN
      END DO
      
      DEALLOCATE(stack%fields, stat= status)
      
      IF (status /= 0) THEN
        CALL set_error('cannot deallocate fields','set_stack_empty_free_memory')
        RETURN    
      END IF
      NULLIFY(stack%fields)
    END IF
    
    stack%name  = ' '
    stack%defined = silja_undefined

  END SUBROUTINE set_stack_empty_free_memory
    
  
  !***************************************************************************  
  
  subroutine prepare_new_averaging_period(stackPtr)
    !
    ! All fields with field_type "averaged" are set to zero, together with 
    ! their accumulation length
    !
    implicit none

    ! Imported parameter
    type(silja_stack), pointer :: stackPtr

    ! Local variables
    integer :: i
    real, dimension(:), pointer :: dataPtr
    type(silja_field_id), pointer :: idPtr
    type(silja_field), pointer :: fieldPtr

    if(fu_fails(defined(stackPtr), 'Undefined stack given','prepare_new_averaging_period'))return

    do i=1, stackPtr%storepointer
      fieldPtr => stackPtr%fields(i)%fp
      if(fu_met_src(fieldPtr) == silam_internal_src) cycle
      if(fu_field_kind(fu_id(fieldPtr)) == averaged_flag)then
        dataPtr => fu_grid_data(fieldPtr)
        idPtr => fu_id(fieldPtr)
        dataPtr(1:fu_number_of_gridpoints(fu_grid(idPtr)))=0.
        call set_accumulation_length(idPtr, zero_interval)
      endif
    end do

  end subroutine prepare_new_averaging_period


  !******************************************************************************
  
  subroutine update_stack_fields(stackPtr, new_reference_time)
    !
    ! All fields should be updated with new parameters given as input
    !
    implicit none

    ! Imported parameter
    type(silja_stack), pointer :: stackPtr
    type(silja_time), optional, intent(in) :: new_reference_time

    ! Local variables
    integer :: i
    type(silja_field_id), pointer :: idPtr
    type(silja_field), pointer :: fieldPtr
    type(silja_interval) :: fcst_len

    if(fu_fails(defined(stackPtr), 'Undefined stack given','update_stack_fields'))return

    do i=1, stackPtr%storepointer
      fieldPtr => stackPtr%fields(i)%fp
      idPtr => fu_id(fieldPtr)
      if(present(new_reference_time))then  
        ! careful: we want to shift the whole ID, not just redefine one of its times
        fcst_len = fu_forecast_length(idPtr)   ! store
        call set_analysis_time(idPtr, new_reference_time)
        call set_valid_time(idPtr, new_reference_time + fcst_len)
      endif
    end do

  end subroutine update_stack_fields

  
  ! ***************************************************************

  ! ***************************************************************
  !
  !
  !      Get data out of a stack.
  !
  !      These routines do not change stack in any way.
  !
  ! ***************************************************************
  ! ***************************************************************
  
  SUBROUTINE find_field_from_stack_by_id(stack, field_id, field, found)
    
    ! Description:
    ! Finds a field from stack by its identification. If the stack is
    ! not defined at all, no error occurs but found is of course set
    ! to false.
    !
    IMPLICIT NONE
    
    ! Imported parameters with intent IN:
    TYPE(silja_field_id), INTENT(in) :: field_id
    
    ! Imported parameters with intent OUT:
    LOGICAL, INTENT(out) :: found

    ! Imported parameters with intent INOUT or POINTER:
    TYPE(silja_stack), POINTER :: stack
    TYPE(silja_field), POINTER :: field
    type(silja_field_id), pointer :: idTmp

    ! Local declarations:
    INTEGER :: i, start_of_search

    found = .false.

    IF (.NOT.ASSOCIATED(stack)) RETURN

    DO i = 1, SIZE(stack%fields)
      field => stack%fields(i)%fp
      IF (.NOT.defined(field)) EXIT
      idTmp => fu_id(field)
      IF (field_id == idTmp) THEN
        found = .true.
        RETURN
      END IF
    END DO

  END SUBROUTINE find_field_from_stack_by_id


  ! ***************************************************************


  SUBROUTINE find_field_from_stack_direct(met_src, &
                                        & quantity,&
                                        & valid_time,&
                                        & stack,&
                                        & field,&
                                        & found, &
                                        & species)
    !
    ! Returns a field from the stack
    !
    ! Language: ANSI Fortran 90
    !
    ! Author: Mikhail Sofiev, FMI
    ! Code is similar to find_field_3d_from_stack of M.Salonoja
    ! 
    IMPLICIT NONE

    ! Returns value of this function:

    ! Imported parameters with intent IN:
    type(meteo_data_source), INTENT(in) :: met_src
    INTEGER, INTENT(in) :: quantity
    TYPE(silja_time), INTENT(in) :: valid_time
    type(silam_species), intent(in), optional :: species
    ! Imported parameters with intent OUT:
    LOGICAL, INTENT(out) :: found

    ! Imported parameters with intent INOUT or POINTER:
    TYPE(silja_stack), POINTER :: stack
    TYPE(silja_field), POINTER :: field

    ! Local declarations:
    INTEGER :: i

    found = .false.

    IF (.NOT.ASSOCIATED(stack)) RETURN

    IF (.NOT.ASSOCIATED(stack%fields)) THEN
      CALL set_error('no fields in this stack','find_field_from_stack_direct')
      RETURN
    END IF

    DO i = 1, SIZE(stack%fields)
      
      field => stack%fields(i)%fp
      
      IF (.NOT.defined(field)) EXIT

      IF (fu_quantity(field) /= quantity) CYCLE

      if(defined(fu_validity_length(field)))then
        if(valid_time /= time_missing)then
          if(.not.fu_between_times(valid_time, &                              ! time
                                 & fu_valid_time(field), &                    ! limit 1
                                 & fu_valid_time(field) + fu_validity_length(field), &   ! limit2
                                 & .true.)) CYCLE                             ! if accept borders
        endif
      else
        if(valid_time /= time_missing .and. valid_time /= fu_valid_time(field)) CYCLE   
      endif

      IF (.not.(met_src == met_src_missing).and..not.(fu_met_src(field) == met_src)) CYCLE

      if(present(species))then
        if(.not. (fu_species(field) == species))cycle
      endif

      found = .true.
      RETURN

    END DO
  END SUBROUTINE find_field_from_stack_direct


  ! ***************************************************************


  SUBROUTINE find_wind_from_stack(stack, field_id, field, found)
    
    ! Description:
    ! Finds a windfield from stack by its identification. If the
    ! stack is not defined at all, no error occurs but found is of
    ! course set to false.
    !
    IMPLICIT NONE
    
    ! Imported parameters with intent IN:
    TYPE(silja_field_id), INTENT(in) :: field_id
    
    ! Imported parameters with intent OUT:
    LOGICAL, INTENT(out) :: found

    ! Imported parameters with intent INOUT or POINTER:
    TYPE(silja_stack), POINTER :: stack
    TYPE(silja_windfield), POINTER :: field

    ! Local declarations:
    INTEGER :: i, start_of_search

    found = .false.

    IF (.NOT.ASSOCIATED(stack)) RETURN

    DO i = 1, SIZE(stack%winds)
      field => stack%winds(i)%fp
      IF (.NOT.defined(field)) EXIT
      IF (field_id == fu_id(field)) THEN
        found = .true.
        RETURN
      END IF
    END DO
  END SUBROUTINE find_wind_from_stack
  

  ! ***************************************************************


  SUBROUTINE find_field_3d_from_stack(met_src, &
                                    & quantity,&
                                    & valid_time,&
                                    & stack,&
                                    & field_3d,&
                                    & found, &
                                    & species)
    !
    ! Searches for a 3d field from the stack
    ! 
    IMPLICIT NONE

    ! Returns value of this function:

    ! Imported parameters with intent IN:
    type(meteo_data_source), INTENT(in) :: met_src
    INTEGER, INTENT(in) :: quantity
    TYPE(silja_time), INTENT(in) :: valid_time
    type(silam_species), intent(in), optional :: species
    ! Imported parameters with intent OUT:
    LOGICAL, INTENT(out) :: found

    ! Imported parameters with intent INOUT or POINTER:
    TYPE(silja_stack), POINTER :: stack
    TYPE(silja_3d_field), POINTER :: field_3d

    ! Local declarations:
    INTEGER :: i

    found = .false.

    IF (.NOT.ASSOCIATED(stack)) RETURN

    IF (.NOT.ASSOCIATED(stack%fields_3d)) THEN
      CALL set_error('no 3d-fields in this stack','find_field_3d_from_stack')
      RETURN
    END IF

    DO i = 1, SIZE(stack%fields_3d)
      
      field_3d => stack%fields_3d(i)%fp
      
      IF (.NOT.defined(field_3d)) EXIT

!      PRINT *, fu_quantity(field_3d)

      IF (fu_quantity(field_3d) /= quantity) CYCLE

!      CALL report(fu_valid_time(field_3d))

      
      ! time_missing request fits all
      IF ((.not. (valid_time == time_missing)) .and. (.not.(fu_valid_time(field_3d) == valid_time))) CYCLE

!      PRINT *, fu_met_src(field_3d)
    !met_src_missing request fits all
      IF (.not.(met_src == met_src_missing).and..not.(fu_met_src(field_3d) == met_src)) CYCLE

      if(present(species))then
        if(defined(species))then
          if(.not. (fu_species(field_3d) == species))cycle
        endif
      endif

      found = .true.
      RETURN

    END DO
    
  END SUBROUTINE find_field_3d_from_stack


  ! ***************************************************************


  SUBROUTINE find_wind_3d_from_stack(met_src, valid_time, stack, wind_3d, found)
    !
    ! Finds 3d wind field from the stack
    ! 
    IMPLICIT NONE

    ! Returns value of this function:

    ! Imported parameters with intent IN:
    type(meteo_data_source), INTENT(in) :: met_src
    TYPE(silja_time), INTENT(in) :: valid_time

    ! Imported parameters with intent OUT:
    LOGICAL, INTENT(out) :: found

    ! Imported parameters with intent INOUT or POINTER:
    TYPE(silja_stack), POINTER :: stack
    TYPE(silja_3d_windfield), POINTER :: wind_3d

    ! Local declarations:
    INTEGER :: i

    found = .false.

    IF (.NOT.ASSOCIATED(stack)) RETURN

    IF (.NOT.ASSOCIATED(stack%winds_3d)) THEN
      CALL set_error('no 3d-winds in this stack'&
     & ,'find_wind_3d_from_stack')
      RETURN
    END IF

    DO i = 1, SIZE(stack%winds_3d)
      
      wind_3d => stack%winds_3d(i)%fp
      
      IF (.NOT.defined(wind_3d)) EXIT

      if (fu_valid_time(wind_3d) /= valid_time) cycle
      !IF (.not. (valid_time == time_missing) &
      !    .and. fu_valid_time(wind_3d) /= valid_time) CYCLE

!      IF (.not.fu_met_src(wind_3d) == met_src) CYCLE
      IF (.not.(met_src == met_src_missing).and..not.(fu_met_src(wind_3d) == met_src)) CYCLE

      found = .true.
      RETURN
    END DO

    found = .false.

  END SUBROUTINE find_wind_3d_from_stack


  !*****************************************************************

  logical function fu_id_in_permanent_stack(id, stack) result(found)
    !
    ! Checks if the given id is already in the permanent stack.
    ! A trick is that the permanent stack disregards time, so it must
    ! not be taken into account
    !
    implicit none

    ! Imported parameters
    type(silja_field_id), intent(in) :: id
    type(silja_stack), intent(in) :: stack

    ! Local variables
    integer :: i
    type(silja_field), pointer :: field

    found = .false.

    DO i = 1, SIZE(stack%fields)
      
      field => stack%fields(i)%fp
      
      IF (.NOT.defined(field)) EXIT

      IF (fu_quantity(field) /= fu_quantity(id)) CYCLE

      IF (.not.(fu_met_src(id) == met_src_missing .or. &
              & fu_met_src(id) == fu_met_src(field))) CYCLE

      if (.not.fu_cmp_levs_eq(fu_level(id), fu_level(field))) cycle
      found = .true.
      RETURN

    END DO

  end function fu_id_in_permanent_stack


  ! ***************************************************************


  FUNCTION fu_first_field_from_stack(stack) result(field)
    !
    ! Returns the first field found in stack. Used to test the
    ! contents of data.
    !
    ! Language: ANSI Fortran 90
    !
    ! Author: Mika Salonoja, FMI
    ! 
    IMPLICIT NONE
    !
    ! The return value of this function:
    TYPE(silja_field), POINTER :: field

    ! Imported parameters with intent(in):
    TYPE(silja_stack), POINTER :: stack
    
    IF (defined(stack)) THEN

      field => stack%fields(1)%fp

      IF (.NOT.defined(field)) THEN
        CALL set_error('first field not found in stack','fu_first_field_from_stack')
      END IF
    END IF
  
  END FUNCTION fu_first_field_from_stack


  ! ***************************************************************


  FUNCTION fu_field_from_stack_by_quantity(stack, quantity, chName) result(field)
    
    ! Description:
    ! Returns the first field found in stack. Used to test the
    ! contents of data.
    !
    IMPLICIT NONE
    !
    ! The return value of this function:
    TYPE(silja_field), POINTER :: field

    ! Imported parameters with intent(in):
    TYPE(silja_stack), POINTER :: stack
    INTEGER, INTENT(in) :: quantity
    character(len=*),intent(in), optional :: chName

    ! Local:
    INTEGER :: i

    IF (defined(stack)) THEN
      DO i = 1, SIZE(stack%fields)
        field => stack%fields(i)%fp
        IF (.NOT.defined(field)) EXIT
        IF (fu_quantity(field) == quantity) then
          if(present(chName))then
            if(trim(fu_substance_name(field)) == trim(chName))return
          else
            RETURN
          endif
        endif
      END DO

      IF (.NOT.defined(field)) THEN
   CALL set_error('field not found in stack'&
       & ,'fu_field_from_stack_by_quantity')
      END IF
    END IF
  
  END FUNCTION fu_field_from_stack_by_quantity


  ! ***************************************************************


  FUNCTION fu_first_3d_field_from_stack(stack) result(field_3d)
    
    ! Description:
    ! Returns the first field 3d-found in stack. 
    !
    ! Language: ANSI Fortran 90
    !
    ! Author: Mika Salonoja, FMI
    ! 
    IMPLICIT NONE
    !
    ! The return value of this function:
    TYPE(silja_3d_field), POINTER :: field_3d

    ! Imported parameters with intent(in):
    TYPE(silja_stack), POINTER :: stack
    
    IF (defined(stack)) THEN
      field_3d => stack%fields_3d(1)%fp
      IF (.NOT.defined(field_3d)) THEN
        CALL set_error('first field not found in stack','fu_first_3d_field_from_stack')
      END IF
    ELSE
      NULLIFY(field_3d)
    END IF
  
  END FUNCTION fu_first_3d_field_from_stack


  ! ***************************************************************


  SUBROUTINE stack_quantities(stack, quantities, number_of_quantities)

    ! Description:
    ! Finds all possible quantities the  2D scalar fields
    ! are found for in stack. Used in v5d-conversion.
    !
    ! Language: ANSI Fortran 90
    !
    ! Author: Mika Salonoja, FMI
    ! 
    IMPLICIT NONE

    ! Imported parameters with intent OUT:
    INTEGER, DIMENSION(:), INTENT(out) :: quantities
    INTEGER, INTENT(out) :: number_of_quantities
  
    ! Imported parameters with intent POINTER:
    TYPE(silja_stack), POINTER  ::  stack
    
    ! Local declarations:
    INTEGER :: i, j
    INTEGER :: quantity
    TYPE(silja_field), POINTER :: field
    LOGICAL :: found_earlier


    quantities = int_missing
    number_of_quantities = 0

    IF (.NOT.defined(stack)) THEN
      CALL set_error('stack empty','stack_quantities')
      RETURN
    END IF

    DO i = 1, SIZE(stack%fields)

      field => stack%fields(i)%fp

      IF (.NOT.defined(field)) EXIT
      quantity = fu_quantity(field)

      ! a new or existing quantity?
      IF (number_of_quantities == 0) THEN
        quantities(1) = quantity
        number_of_quantities = 1
      ELSE

        found_earlier = .false.
        DO j = 1, number_of_quantities
         IF (quantity == quantities(j)) found_earlier = .true.
        END DO

        IF (found_earlier) CYCLE

        ! A defined time not found earlier in the list:

        number_of_quantities = number_of_quantities + 1
        IF (number_of_quantities > SIZE(quantities)) THEN
          CALL set_error('quantities-vector too small','stack_quantities')
          RETURN
        ELSE
          quantities(number_of_quantities) = quantity
        END IF
      END IF
    END DO

  END SUBROUTINE stack_quantities


  ! ***************************************************************


  SUBROUTINE stack_3d_quantities(stack, quantities, number_of_quantities)

    ! Description:
    ! Finds all possible quantities the 3d scalarfields
    ! are found for in stack. Used in v5d-conversion.
    !
    ! Language: ANSI Fortran 90
    !
    ! Author: Mika Salonoja, FMI
    ! 
    IMPLICIT NONE

    ! Imported parameters with intent OUT:
    INTEGER, DIMENSION(:), INTENT(out) :: quantities
    INTEGER, INTENT(out) :: number_of_quantities
  
    ! Imported parameters with intent POINTER:
    TYPE(silja_stack), POINTER  ::  stack
    
    ! Local declarations:
    INTEGER :: i, j
    INTEGER :: quantity
    TYPE(silja_3d_field), POINTER :: field
    LOGICAL :: found_earlier

    quantities = int_missing
    number_of_quantities = 0

    IF (.NOT.defined(stack)) THEN
      CALL set_error('stack empty','stack_3d_quantities')
      RETURN
    END IF

    IF (.NOT.ASSOCIATED(stack%fields_3d)) RETURN

    DO i = 1, SIZE(stack%fields_3d)

      field => stack%fields_3d(i)%fp

      IF (.NOT.defined(field)) EXIT
      quantity = fu_quantity(field)

      ! A new or existing quantity?
      IF (number_of_quantities == 0) THEN
        quantities(1) = quantity
        number_of_quantities = 1
      ELSE

        found_earlier = .false.
        DO j = 1, number_of_quantities
          IF (quantity == quantities(j)) found_earlier = .true.
        END DO

        IF (found_earlier) CYCLE

        ! A defined time not found earlier in the list:

        number_of_quantities = number_of_quantities + 1
        IF (number_of_quantities > SIZE(quantities)) THEN
          CALL set_error('quantities-vector too small','stack_3d_quantities')
          RETURN
        ELSE
          quantities(number_of_quantities) = quantity
        END IF
      END IF
    END DO

  END SUBROUTINE stack_3d_quantities


  ! ***************************************************************


  SUBROUTINE stack_variables(stack, quantities, arSpecies, number_of_quantities)

    ! Description:
    ! Finds all possible quantities the  2D scalar fields
    ! are found for in stack. Used in v5d-conversion.
    !
    ! Language: ANSI Fortran 90
    !
    ! Author: Mika Salonoja, FMI
    ! 
    IMPLICIT NONE

    ! Imported parameters with intent OUT:
    INTEGER, DIMENSION(:), INTENT(out) :: quantities
    INTEGER, INTENT(out) :: number_of_quantities
    type(silam_species), dimension(:), intent(out) :: arSpecies
    TYPE(silja_stack), intent(in)  ::  stack
    
    ! Local declarations:
    INTEGER :: i, j
    INTEGER :: quantity
    TYPE(silja_field), POINTER :: field
    LOGICAL :: found_earlier
    type(silam_species) :: speciesTmp

    quantities = int_missing
    number_of_quantities = 0

    IF (.NOT.defined(stack)) THEN
      CALL set_error('stack empty','stack_variables')
      RETURN
    END IF

    DO i = 1, SIZE(stack%fields)

      field => stack%fields(i)%fp

      IF (.NOT.defined(field)) EXIT
      quantity = fu_quantity(field)
      speciesTmp = fu_species(field)

      ! a new or existing quantity?
      found_earlier = .false.
      DO j = 1, number_of_quantities
        IF (quantity == quantities(j))then
          if(speciesTmp == arSpecies(j)) found_earlier = .true.
        endif
      END DO

      IF (found_earlier) CYCLE

      ! Add the variable

      IF (number_of_quantities == SIZE(quantities)) THEN
        CALL set_error('quantities-vector too small','stack_variables')
        RETURN
      ELSE
        number_of_quantities = number_of_quantities + 1
        arSpecies(number_of_quantities) = speciesTmp
        quantities(number_of_quantities) = quantity  ! number_of_quantities has just been increased
      END IF
    END DO

  END SUBROUTINE stack_variables


  ! ***************************************************************

  subroutine  check_stack_fields_ranges(stack)
    IMPLICIT NONE
    !
    ! Imported parameters with intent(in):
    TYPE(silja_stack), POINTER :: stack
    
    integer :: iTmp, fs, nx, n_out_of_range, npatched, n_failed, quantity
    TYPE(silja_field), POINTER :: field
    real, dimension(:), pointer :: grid_data

    if(stack%storepointer < 1)return  ! nothing to do in thsi stack
    call msg('Checking stack:' + stack%name)
    fs = fu_number_of_gridpoints(fu_grid(stack%fields(1)%fp))
    
    DO iTmp = 1, stack%storepointer
      field => stack%fields(iTmp)%fp
      if(.not. defined(field))exit
      grid_data => fu_grid_data(field)
      quantity = fu_quantity(field)
      call msg('Checking:' + fu_quantity_string(quantity))
      call check_quantity_range(quantity, grid_data, fs, nx, &
                              & .false., .false., &  ! ifRequireValidity, ifSilent
                              & n_out_of_range, npatched, n_failed)
    end do
    call msg('')
    
  end subroutine  check_stack_fields_ranges

  
  !*****************************************************************

  LOGICAL FUNCTION fu_stack_full(stack)
    !
    ! Returns true value if no more space for scalar 2d-fields.
    ! 
    IMPLICIT NONE
    !
    ! Imported parameters with intent(in):
    TYPE(silja_stack), POINTER :: stack
 
    IF (defined(stack)) THEN
      fu_stack_full = stack%full
    ELSE
      fu_stack_full = .false.
    END IF

  END FUNCTION fu_stack_full


  ! ***************************************************************


  LOGICAL FUNCTION fu_stack_empty(stack)
    
    ! Description:
    ! Returns true value if stack is undefined or empty
    ! 
    ! Language: ANSI Fortran 90
    !
    ! Author: Mika Salonoja, FMI
    ! 
    IMPLICIT NONE
    !
    ! Imported parameters with intent(in):
    TYPE(silja_stack), POINTER :: stack
 
    IF (defined(stack)) THEN
      fu_stack_empty = (stack%storepointer ==0)
    else
      fu_stack_empty = .true.
    END IF

  END FUNCTION fu_stack_empty


  ! ***************************************************************


  LOGICAL FUNCTION fu_one_valid_time_stack(stack)
    ! 
    ! Returns true value if stack can only contain data for single
    ! time
    ! 
    IMPLICIT NONE
    !
    ! Imported parameters with intent(in):
    TYPE(silja_stack), POINTER :: stack
 
    fu_one_valid_time_stack = stack%one_valid_time

  END FUNCTION fu_one_valid_time_stack


  ! ***************************************************************


  LOGICAL FUNCTION fu_single_met_src_stack(stack)
    ! 
    ! Returns true value if stack can only contain data for single
    ! time
    ! 
    IMPLICIT NONE
    !
    ! Imported parameters with intent(in):
    TYPE(silja_stack), POINTER :: stack
 
    fu_single_met_src_stack = stack%single_met_src

  END FUNCTION fu_single_met_src_stack


  ! ***************************************************************


  INTEGER FUNCTION fu_number_of_wind_3d_fields(stack)

    ! Returns the number of defined wind 3d fields in stack.
    !
    IMPLICIT NONE
    !
    ! Imported parameters with intent(in):
    TYPE(silja_stack), POINTER :: stack
    !
    fu_number_of_wind_3d_fields = stack%windstorepointer_3d

  END FUNCTION fu_number_of_wind_3d_fields



  ! ***************************************************************


  INTEGER FUNCTION fu_number_of_wind_2d_fields(stack)

    ! Returns the number of defined wind 3d fields in stack.
    !
    IMPLICIT NONE
    !
    ! Imported parameters with intent(in):
    TYPE(silja_stack), POINTER :: stack

    fu_number_of_wind_2d_fields = stack%windstorepointer

  END FUNCTION fu_number_of_wind_2d_fields




  ! ***************************************************************


  INTEGER FUNCTION fu_number_of_3d_fields(stack)

    ! Returns the number of defined scalar 3d fields in stack.
    !
    IMPLICIT NONE
    !
    ! Imported parameters with intent(in):
    TYPE(silja_stack), POINTER :: stack
    !
    ! Local declarations:
    INTEGER :: i

    fu_number_of_3d_fields = stack%storepointer_3d

  END FUNCTION fu_number_of_3d_fields



  ! ***************************************************************


  INTEGER FUNCTION fu_number_of_2d_fields(stack)

    ! Returns the number of defined scalar fields in stack.
    !
    IMPLICIT NONE
    !
    ! Imported parameters with intent(in):
    TYPE(silja_stack), POINTER :: stack

    fu_number_of_2d_fields = stack%storepointer

  END FUNCTION fu_number_of_2d_fields



  ! ****************************************************************

  ! ****************************************************************
  ! 
  !
  !       Private functions and subroutines of this module.
  !
  !
  ! ***************************************************************
  ! ***************************************************************

  LOGICAL FUNCTION fu_stack_defined(stack)
    
    ! 
    IMPLICIT NONE
    TYPE(silja_stack), intent(in) :: stack
    
    fu_stack_defined = fu_true(stack%defined)

  END FUNCTION fu_stack_defined


  !*****************************************************************

  subroutine set_stack_defined(stack, defined)
    ! 
    ! Stupid FORTRAN sometimes sets stack%defined to silja_true already at
    ! the allocation moment, which then prevents the stack from being initialised
    ! normally. Has to force it to be undefined
    !
    IMPLICIT NONE
    
    ! Imported parameters with intent IN:
    TYPE(silja_stack), intent(inout) :: stack
    type(silja_logical), intent(in) :: defined
    
    stack%defined = defined

  END subroutine set_stack_defined


  ! ***************************************************************


  FUNCTION fu_stack_valid_time(stack)
    
    ! Description:
    ! Returns the valid time of all data kept in the stack. If stack
    ! is not one_valid -type, then undefined time is returned.
    !
    ! Language: ANSI Fortran 90
    !
    ! Author: Mika Salonoja, FMI
    ! 
    IMPLICIT NONE
    !
    ! The return value of this function:
    TYPE(silja_time) :: fu_stack_valid_time
    !
    ! Imported parameters with intent(in):
    TYPE(silja_stack), POINTER :: stack

    IF (defined(stack)) THEN
      IF (stack%one_valid_time) THEN
        fu_stack_valid_time = stack%valid_time
      ELSE
        fu_stack_valid_time = time_missing
      END IF
    ELSE
      fu_stack_valid_time = time_missing
    END IF

  END FUNCTION fu_stack_valid_time


  ! ***************************************************************


  FUNCTION fu_stack_met_src(stack) result(met_src)
    
    ! Description:
    ! Returns the met_src of all data kept in the stack. If stack
    ! is not single-met_src -type, error occurs.
    !
    IMPLICIT NONE
    !
    type(meteo_data_source) :: met_src

    ! Imported parameters
    TYPE(silja_stack), pointer :: stack

    IF (defined(stack)) THEN
      IF (stack%single_met_src) THEN
        met_src = stack%mds
      ELSE
!        CALL set_error('not a single met_src stack','fu_stack_met_src')
         met_src = met_src_missing
      END IF
    ELSE
      met_src = met_src_missing
    END IF

  END FUNCTION fu_stack_met_src

  ! ***************************************************************


  FUNCTION fu_stack_wdr(stack) result(wdr)
    
    ! Description:
    ! Returns the met_src of all data kept in the stack. If stack
    ! is not single-met_src -type, error occurs.
    !
    IMPLICIT NONE
    !
    type(silja_wdr) :: wdr
    ! Imported parameters with intent(in):
    TYPE(silja_stack), POINTER :: stack

    IF (defined(stack)) THEN
      wdr = stack%wdr_ptr
    ELSE
      wdr = wdr_missing
    END IF

  END FUNCTION fu_stack_wdr
  ! ***************************************************************


  FUNCTION fu_name_of_stack(stack) result(name)
    
    IMPLICIT NONE
    !
    character(len=clen) :: name
    ! Imported parameters with intent(in):
    TYPE(silja_stack), POINTER :: stack

    IF (defined(stack)) THEN
      name = stack%name
    ELSE
      name = '***UNDEFINED***'
    END IF

  END FUNCTION fu_name_of_stack



  ! ***************************************************************

  
  SUBROUTINE arrange_uvw_to_windfields(stack)

    ! Description:
    ! Finds corresponding windfield-components, and makes 2d
    ! -windfields out of them, and then stores them in the stack. No
    ! data is duplicated, only pointers are set.
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

    ! Imported parameters with intent INOUT or POINTER:
    TYPE(silja_stack), POINTER :: stack

    ! Local declarations:
    INTEGER :: i
    INTEGER :: v_index, w_index
    LOGICAL :: v_found, w_found, wind_found
    INTEGER :: place
    TYPE(silja_field), POINTER :: u, v, w
    TYPE(silja_windfield), POINTER :: wind

    !----------------------------------------
    !
    ! 1. Loop over scalar fields.
    !    -----------------------

    loop_over_scalars: DO i = 1, SIZE(stack%fields)

      u => stack%fields(i)%fp
      IF (.NOT.defined(u)) then 
        EXIT loop_over_scalars
      end IF
      u_check: IF (fu_quantity(u) == u_flag) THEN
   !----------------------------------------
   !
   ! 2. Start check for one U.
   !    ---------------------

   ! 2.1. First check that the windfield, whose component this u
   ! -field is (same met_src, level and time) isn't already among
   ! the 2d-windfields of this stack:

        CALL corresponding_wind_in_stack(u, place, wind_found)

        !-------------------------------------------
        ! 
        ! 2.2. If windfield already contains w, do nothing

        IF (wind_found) THEN
          wind => stack%winds(place)%fp
          IF (fu_w_available(wind)) CYCLE loop_over_scalars
        END IF

        !--------------------------------------------------------
        !
        ! 2.3. Find corresponding v and w:

        CALL search_match()
        IF (error) RETURN
        !----------------------------------------------
        !
        ! 3. No existing windfield in stack
        !    ------------------------------

        new_windfield: IF (.NOT.wind_found) THEN

          match_found: IF (v_found) THEN

            v => stack%fields(v_index)%fp

            ! -----------------------------------------
            ! 
            ! 3.2. Match found, search in winds a place where to store

            place = stack%windstorepointer + 1
            !    PRINT *, 'place :', place

            IF (place > SIZE(stack%winds)) THEN
              CALL report(stack)
              CALL set_error('no more space for windfields','arrange_uvw_to_windfields')
              RETURN
            END IF


            ! -----------------------------------------
            ! 
            ! 3.3. Place found, set field in.

            w_or_no_w: IF (w_found) THEN
              w => stack%fields(w_index)%fp
              CALL set_windfield_from_uvw(u, v, stack%winds(place)%fp, w)
              IF (error) RETURN

!!!$              CALL set_windfield_from_uvw(&
!!!$                  & u,&
!!!$                  & v, &
!!!$                  & stack%winds(place)%fp)
!!!$              IF (error) RETURN
              !call report(u)
              stack%windstorepointer = stack%windstorepointer + 1

            ELSE
              CALL set_windfield_from_uvw(u, v, stack%winds(place)%fp)
              IF (error) RETURN

              stack%windstorepointer = stack%windstorepointer + 1

            END IF w_or_no_w

          ELSE ! match_found
            call msg("")
            call report(u)
            CALL msg_warning('no matching v found for u')
            call msg("")
            CYCLE loop_over_scalars
          END IF match_found

        ELSE

          !-------------------------------------------------
          !
          ! 4. Windfield exists already, so check it
          !    ------------------------------------

          ! --------------------------------------------------
          !
          ! 4.1. No W in winfield, check for scalar W

          IF (.NOT. w_found) CYCLE loop_over_scalars

          ! -----------------------------------------------
          !
          ! 4.3. Scalar W found, add it to windfield

          w => stack%fields(w_index)%fp
          wind => stack%winds(place)%fp
!          call msg('Adding w to wind')
        !  CALL report(w)
        !  CALL report(wind)

          CALL add_w_to_windfield(w,wind)
          IF (error)THEN
            CALL msg_warning('% No match found for w %')
            RETURN
          END IF

          wind => stack%winds(place)%fp
       !   CALL report(wind)
       !   STOP
        END IF new_windfield

      ELSE
   CYCLE loop_over_scalars
      END IF u_check
      
    END DO loop_over_scalars
    
  CONTAINS
    
    ! *** private functions of arrange_uvw_to_windfields **


    SUBROUTINE search_match()

      ! Searches appropriate v, w fields for given u-field.
      ! ATTENTION. There are many possible types of w-fields:
      ! The routine checks them one-by-one
      !
      ! Language: ANSI Fortran 90
      !
      ! Original code : Mika Salonoja
      ! Author: Mikhail Sofiev, FMI
      ! 
      IMPLICIT NONE
      
      ! Local declarations:
      INTEGER :: j, w_q_tmp
      TYPE(silja_time) :: u_valid
      TYPE(silja_level) :: u_level
      type(meteo_data_source) :: u_met_src
      TYPE(silja_field), POINTER :: field

      v_found = .false.
      w_found = .false.
      u_valid = fu_valid_time(u)
      u_level = fu_level(u)
      u_met_src = fu_met_src(u)

      select case(fu_leveltype(u_level))
        case(constant_altitude)
          w_q_tmp = w_alt_msl_flag
        case(constant_height)
          w_q_tmp = w_height_srf_flag
        case(constant_pressure, hybrid, sigma_level, layer_btw_2_hybrid)
          w_q_tmp = omega_flag
        case default
          call msg_warning('Non-supported level','search_match@arrange_uvw_to_windfields')
          call report(u_level)
          call set_error('Non-supported level','search_match@arrange_uvw_to_windfields')
          return
      end select

      DO j = 1, SIZE(stack%fields)
      
        IF (i==j) CYCLE

        field => stack%fields(j)%fp

        IF (.NOT.defined(field)) EXIT

        ! V found?
        IF (fu_quantity(field) == v_flag) THEN

          IF (.NOT.fu_cmp_levs_eq(fu_level(field), u_level)) CYCLE
          IF (.NOT.(fu_valid_time(field) == u_valid)) CYCLE
          IF (.NOT.(fu_met_src(field) == u_met_src)) CYCLE

          v_found = .true.
          v_index = j
        END IF

        ! Verical velocity found?
        IF (fu_quantity(field) == w_q_tmp) THEN

          IF (.NOT.fu_cmp_levs_eq(fu_level(field), u_level)) CYCLE
          IF (.NOT.(fu_valid_time(field) == u_valid)) CYCLE
          IF (.NOT.(fu_met_src(field) == u_met_src)) CYCLE

          w_found = .true.
          w_index = j
        END IF
   
        IF (v_found.and.w_found) EXIT
      END DO

!      if(.not.v_found) call set_error('Failed to find proper v-wind', &
!                                    & 'search_match@arrange_uvw_to_windfields')
!      if(.not.w_found) call set_error('Failed to find proper vertical wind', &
!                                    & 'search_match@arrange_uvw_to_windfields')
    END SUBROUTINE search_match


    ! ***************************************************************


    SUBROUTINE corresponding_wind_in_stack(u, place, wind_found)
    
      ! Description:
      ! Returns true value if one of the given u-field is already a
      ! component of one of the windfieds in stack.
      !
      ! Language: ANSI Fortran 90
      !
      ! Author: Mika Salonoja, FMI
      ! 
      IMPLICIT NONE
      !
      ! Imported parameters with intent(in):
      TYPE(silja_field), POINTER :: u

      ! Imported parameters with intent(out):
      INTEGER, INTENT(out) :: place
      LOGICAL, INTENT(out) :: wind_found

      ! Local declarations:
      INTEGER :: i
      TYPE(silja_level) :: u_level
      TYPE(silja_time) :: u_valid
      type(meteo_data_source) :: u_met_src
      TYPE(silja_windfield), POINTER :: windfield

      wind_found = .false.
      place = 0

      u_level = fu_level(u)
      u_valid = fu_valid_time(u)
      u_met_src = fu_met_src(u)

      DO i = 1, SIZE(stack%winds)

        windfield => stack%winds(i)%fp
        IF (.NOT.defined(windfield)) EXIT

        IF (.NOT.fu_cmp_levs_eq(fu_level(windfield), u_level)) CYCLE

        IF (.NOT.stack%one_valid_time) THEN
          IF (.NOT.(fu_valid_time(windfield) == u_valid)) CYCLE
        END IF

        ! At the moment wind-components are not arranged:
        IF (.NOT.fu_met_src(windfield) == u_met_src) CYCLE

        place = i
        wind_found = .true.
        EXIT
      END DO

    END SUBROUTINE corresponding_wind_in_stack
    
  END SUBROUTINE arrange_uvw_to_windfields



  ! ***************************************************************

  
  SUBROUTINE arrange_uv_10m_to_windfields(stack)

    ! Description:
    ! Finds corresponding 10m windfield-components, and makes 2d
    ! -windfields out of them, and then stores them in the stack. No
    ! data is duplicated, only pointers are set.
    ! 
    ! All units: SI
    !
    ! Language: ANSI Fortran 90
    !
    ! Original code: Mika Salonoja, FMI
    ! Corrected by M.Sofiev
    ! 
    IMPLICIT NONE

    ! Imported parameters with intent INOUT or POINTER:
    TYPE(silja_stack), POINTER :: stack

    ! Local declarations:
    INTEGER :: i
    INTEGER :: v_index
    LOGICAL :: v_found
    INTEGER :: place
    TYPE(silja_field), POINTER :: u, v
    TYPE(silja_windfield), POINTER :: wind

    !----------------------------------------
    !
    ! 1. Loop over scalar fields.
    !    -----------------------

    loop_over_scalars: DO i = 1, SIZE(stack%fields)
    
      u => stack%fields(i)%fp
      IF (.NOT.defined(u)) EXIT
      
      IF (fu_quantity(u) /= u_10m_flag) CYCLE
     
      !----------------------------------------
      !
      ! 2. Start check for one 10m U.
      !    -------------------------
      
      ! 2.1. First check that the windfield, whose component this u
      ! -field is (same met_src, level and time) isn't already among
      ! the 2d-windfields of this stack:
      IF (fu_corresp_wind_10m_in_stack(u)) CYCLE

      ! 2.2. Find corresponding v:
      CALL search_match()
      IF (error) RETURN
      
      match_found: IF (v_found) THEN
        
        ! -----------------------------------------
        ! 
        ! 2.3. Match found, search in winds a place where to store
        
        place = stack%windstorepointer + 1
        !    PRINT *, 'place :', place
        
        IF (place > SIZE(stack%winds)) THEN
          CALL report(stack)
          CALL set_error('no more space for windfields'&
                & ,'arrange_uv_10m_to_windfields')
          RETURN
        END IF
        
        ! -----------------------------------------
        ! 
        ! 2.4. Place found, set field in.
        
        v => stack%fields(v_index)%fp
        
        CALL set_windfield_from_uvw(&
             & u,&
             & v, &
             & stack%winds(place)%fp)
        IF (error) RETURN
        
        stack%windstorepointer = stack%windstorepointer + 1
        
      ELSE ! match_found
        CALL msg_warning('no matching v10 found for u10')
        CYCLE
      END IF match_found
      
    END DO loop_over_scalars

    
  CONTAINS
    
    ! *** private functions of arrange_uv_10m_to_windfields **


    SUBROUTINE search_match()

      ! Description:
      !
      ! Language: ANSI Fortran 90
      !
      ! Author: Mika Salonoja, FMI
      ! 
      IMPLICIT NONE
      
      ! Local declarations:
      INTEGER :: j
      TYPE(silja_time) :: u_valid
      type(meteo_data_source) :: u_met_src
      TYPE(silja_field), POINTER :: field

      v_found = .false.
      u_valid = fu_valid_time(u)
      u_met_src = fu_met_src(u)

      DO j = 1, SIZE(stack%fields)
      
        IF (i==j) CYCLE

        field => stack%fields(j)%fp

        IF (.NOT.defined(field)) EXIT

        ! V found?
        IF (fu_quantity(field) == v_10m_flag) THEN

          IF (.NOT.(fu_valid_time(field) == u_valid)) CYCLE
          IF (.NOT.(fu_met_src(field) == u_met_src)) CYCLE

          v_found = .true.
          v_index = j
          EXIT ! found ok!!
        END IF

      END DO

    END SUBROUTINE search_match


    ! ***************************************************************


    LOGICAL FUNCTION fu_corresp_wind_10m_in_stack(u)
    
      ! Description:
      ! Returns true value if one of the given u-field is already a
      ! component of one of the windfieds in stack.
      !
      ! Language: ANSI Fortran 90
      !
      ! Author: Mika Salonoja, FMI
      ! 
      IMPLICIT NONE
      !
      ! Imported parameters with intent(in):
      TYPE(silja_field), POINTER :: u

      ! Local declarations:
      INTEGER :: i
      TYPE(silja_level) :: u_level
      TYPE(silja_time) :: u_valid
      type(meteo_data_source) :: u_met_src
      TYPE(silja_windfield), POINTER :: windfield

      fu_corresp_wind_10m_in_stack = .false.

      u_level = fu_level(u)
      u_valid = fu_valid_time(u)
      u_met_src = fu_met_src(u)

      DO i = 1, SIZE(stack%winds)

   windfield => stack%winds(i)%fp
   IF (.NOT.defined(windfield)) EXIT

   IF (.NOT.fu_cmp_levs_eq(fu_level(windfield), u_level)) CYCLE

   IF (.NOT.stack%one_valid_time) THEN
     IF (.NOT.(fu_valid_time(windfield) == u_valid)) CYCLE
   END IF

   IF (.NOT.fu_met_src(windfield) == u_met_src) CYCLE

   fu_corresp_wind_10m_in_stack = .true.
   EXIT

      END DO
     
    END FUNCTION fu_corresp_wind_10m_in_stack
    
  END SUBROUTINE arrange_uv_10m_to_windfields




  ! ***************************************************************


  SUBROUTINE arrange_scalar_fields_to_3d(stack, ifLookForPressure)

    ! Description:
    ! Finds corresponding scalar fields (same quantity, met_src and
    ! time, but different level) and makes 3d-fields out of them, and
    ! then stores them in the stack. No
    ! data is duplicated, only pointers are set.
    !  
    IMPLICIT NONE

    ! Imported parameters with intent INOUT or POINTER:
    TYPE(silja_stack), POINTER :: stack
    logical, intent (in) :: ifLookForPressure

    ! Local declarations:
    INTEGER :: i, j, k
    LOGICAL :: found
    INTEGER :: ind
    TYPE(silja_field), POINTER :: field
    TYPE(silja_3d_field), POINTER :: field_3d
    TYPE(silja_stack), POINTER :: stack_
    character(len=*), parameter :: sub_name="arrange_scalar_fields_to_3d"

    stack_ => stack

    ! -------------------------------------------------------------
    ! 
    ! 1. Loop over 2d scalar fields in stack
    !
    loop_over_2d_fields: DO i = 1, SIZE(stack%fields)

      IF (error) RETURN

      field => stack%fields(i)%fp

      IF (.NOT.defined(field)) EXIT

!      CALL report(field)

      ! -------------------------------------------------------------
      ! 
      ! 2. Check if this field should be added to 3d-fields
      !
      IF (.NOT. fu_multi_level_quantity(fu_quantity(field))) CYCLE


      ! take only stuff of current source's primary levels:
      IF (.NOT.fu_3d_leveltype(fu_leveltype(fu_level(field)))) CYCLE


      ! -------------------------------------------------------------
      ! 
      ! 3. Find the 3d-field into which this field should be added to
      !
      CALL search_match()

      ! If a suitable 3d-field not found, then start filling a new one:
      IF (.NOT.found) THEN
        stack%storepointer_3d = stack%storepointer_3d + 1
        IF (stack%storepointer_3d > SIZE(stack%fields_3d)) THEN
          stack%storepointer_3d = stack%storepointer_3d - 1
          call set_error("Too many fields came", sub_name)
          CALL search_match()


          call msg("Problem fitting field " + fu_quantity_string(fu_quantity(field)) )
          CALL report(field)
          call msg('================================')
          call msg('================================')
          CALL report(stack)
          CALL set_error('no space for fields_3d in stack', sub_name)
          RETURN
        END IF
        ind = stack%storepointer_3d
      END IF


      ! -------------------------------------------------------------
      ! 
      ! 4. Add the current field to the suitable 3d-field.
      !
      CALL add_field_to_3d_field(field, stack%fields_3d(ind)%fp)

    END DO loop_over_2d_fields


    ! -------------------------------------------------------------
    ! 
    ! 5. In each 3d-field organize the fields vertically.
    !
    i=1
    do while (i <= stack%storepointer_3d)

      field_3d => stack%fields_3d(i)%fp

      !
      ! We have to check if this very field_3d is really 3d, which
      ! means that it contains more than one 2d field.
      ! Reason - quantities are considered 3d if they can be defined
      ! at several levels. But this is not necessary that in specific
      ! case they are really 3d
      !
      if(fu_number_of_fields(field_3d) < 2)then
        !
        ! Squeeze the list of fields
        !
        if(i < stack%storepointer_3d)then
          do j=i, stack%storepointer_3d-1

            field_3d => stack%fields_3d(j)%fp
            call set_3d_field_empty(field_3d) ! Kill j-th 3d field

            field_3d => stack%fields_3d(j+1)%fp
            do k=1,fu_number_of_fields(field_3d)
              call add_field_to_3d_field(fu_field_from_3d_field(field_3d,k,.true.), & ! skip order
                                       & stack%fields_3d(j)%fp)
            end do

          end do
        endif
        call set_3d_field_empty(stack%fields_3d(stack%storepointer_3d)%fp)
        stack%storepointer_3d = stack%storepointer_3d - 1

      endif

      if(i <= stack%storepointer_3d) CALL organize_fields_vertically(stack%fields_3d(i)%fp)
      IF (error) RETURN

      i=i+1

    END DO


    ! -------------------------------------------------------------
    ! 
    ! 6. Add surface pressure fields to those 3d-fields that contain
    ! hybrid-level data
    !
    if (ifLookForPressure) then
            DO i = 1, stack%storepointer_3d
              field_3d => stack%fields_3d(i)%fp
              j = fu_leveltype(field_3d) 
              IF ( any(j ==  (/hybrid, sigma_level,layer_btw_2_sigma,layer_btw_2_hybrid/))) THEN
                CALL find_matching_surface_pressure(field_3d, stack)
                IF (error) RETURN
              elseif (j /= layer_btw_2_height) then
                call msg_warning("Unhandled level type for 3d field", sub_name)
                call report(field_3d)
                call set_error("Unhandled level type for 3d field", sub_name)
                return
              END IF
            END DO
    endif


  CONTAINS

    SUBROUTINE search_match()

      ! Description:
      ! Finds fo a field the corresponding 3d-field from a set of 3d
      ! fields.

      IMPLICIT NONE

      ! Local declarations:
      INTEGER :: i
      INTEGER :: quantity
      type(meteo_data_source) :: met_src
      TYPE(silja_time) :: time
      TYPE(silja_3d_field), POINTER :: field_3d
      character(len = substNmLen) :: subst
      real :: mode
      logical :: has_species

      found = .false.
      quantity = fu_quantity(field)
      met_src = fu_met_src(field)
      time = fu_valid_time(field)
      has_species = defined(fu_species(field))

      DO i = 1, SIZE(stack%fields_3d)
        field_3d => stack%fields_3d(i)%fp
        IF (.NOT.defined(field_3d)) EXIT
        IF (quantity /= fu_quantity(field_3d)) CYCLE
        IF (time /= fu_valid_time(field_3d)) CYCLE
        IF (.not.met_src == fu_met_src(field_3d)) CYCLE
        if (.not. has_species .eqv. defined(fu_species(field_3d))) cycle
        if (has_species) then
          if (.not. fu_species(field_3d) == fu_species(field)) cycle
        end if
        found = .true.
        ind = i
        RETURN
      END DO

    END SUBROUTINE search_match

    ! *****************************************************


    SUBROUTINE find_matching_surface_pressure(field_3d, stack)

      ! Description:
      ! For the unfinished 3d-field finds its corresponding surface
      ! pressure field. This is used for hybrid level data, where every
      ! 3d-field  contains also a pointer to to matching surface
      ! pressure field.
      !
      !
      ! Language: ANSI Fortran 90
      !
      ! Author: Mika Salonoja, FMI
      ! 
      IMPLICIT NONE

      ! Imported parameters with intent INOUT or POINTER:
      TYPE(silja_3d_field), POINTER :: field_3d
      TYPE(silja_stack), POINTER :: stack

      ! Local declarations:
      INTEGER :: i
      TYPE(silja_field), POINTER :: field

      !tilap
      INTEGER :: quantity
      TYPE(silja_level) :: level
!call msg('enter:'+ fu_quantity_short_string(fu_quantity(field_3d)))
      loop_over_2d_fields: DO i = 1, SIZE(stack%fields)
!call msg('loop')
      field => stack%fields(i)%fp
      IF (.NOT.defined(field)) EXIT
!call msg('defined:'+fu_quantity_short_string(fu_quantity(field)))


      ! Hack! surface_pressure_flag and ground_pressure_flag should be the same!!!
      IF (.NOT.(fu_quantity(field) == surface_pressure_flag .or. &
              &  fu_quantity(field) == ground_pressure_flag)) CYCLE
!call msg('quantity')
!        IF (.NOT.(fu_level(field) == ground_level)) CYCLE 
      IF (.NOT.stack%one_valid_time) THEN
          IF (.NOT.(fu_valid_time(field) == fu_valid_time(field_3d))) then
                  call msg("fu_valid_time(field):")
                  call msg(fu_str(fu_valid_time(field)))
                  call msg("fu_valid_time(field_3d):")
                  call msg(fu_str(fu_valid_time(field_3d)))
                CYCLE
           endif

      END IF
!call msg('valid time')
!        IF (.NOT.(fu_met_src(field) == fu_met_src(field_3d))) CYCLE
!call msg('src')
        ! Correct field found:
      CALL add_surf_pre_to_3d_field(field_3d, field)
      RETURN

      END DO loop_over_2d_fields
      call msg("Field: ")
      call report(field_3d)
      CALL msg_warning('no surface pressure found for hybrid 3D','find_matching_surface_pressure')
!      stop
      call msg("Stack: "+ fu_name(stack))
!      call report(stack)
!      call msg("")
!      call msg("")
!      call msg("")
!      call msg("")
!      call msg("")

    END SUBROUTINE find_matching_surface_pressure

  END SUBROUTINE arrange_scalar_fields_to_3d


  ! ***************************************************************


  SUBROUTINE arrange_windfields_to_3d(stack)

    ! Description:
    ! Finds corresponding windfields (same met_src and
    ! time, but different level) and makes 3d-fields out of them, and
    ! then stores them in the stack. No
    ! data is duplicated, only pointers are set.
    ! 
    ! Method:
    ! 
    ! All units: SI
    !
    ! Language: ANSI Fortran 90
    !
    IMPLICIT NONE

    ! Imported parameters with intent INOUT or POINTER:
    TYPE(silja_stack), POINTER :: stack

    ! Local declarations:
    INTEGER :: i, j, k
    LOGICAL :: found
    INTEGER :: ind
    TYPE(silja_windfield), POINTER :: windfield
    TYPE(silja_3d_windfield), POINTER :: wind_3d

    type(silja_stack), pointer :: st

    st => stack

    ! -------------------------------------------------------------
    ! 
    ! 1. Loop over 2d scalar windfields in stack
    !
    loop_over_2d_windfields: DO i = 1, SIZE(stack%winds)

      windfield => stack%winds(i)%fp

      ! -------------------------------------------------------------
      ! 
      ! 2. Check if this field should be added to 3d-fields
      !
      IF (.NOT.defined(windfield)) EXIT


      ! take only stuff of current cource's primary levels:
      IF (.NOT.fu_3d_leveltype(fu_leveltype(fu_level(windfield)))) CYCLE


      ! -------------------------------------------------------------
      ! 
      ! 3. Find the 3d-field into which this field should be added to
      !
      CALL search_match()

      ! If a suitable 3d-field not found, then start filling a new one:
      IF (.NOT.found) THEN
        stack%windstorepointer_3d = stack%windstorepointer_3d + 1
        IF (stack%windstorepointer_3d > SIZE(stack%winds_3d)) THEN
           CALL set_error('no space for winds_3d in stack','arrange_windfields_to_3d')
           RETURN
        END IF
        ind = stack%windstorepointer_3d
      END IF


      ! -------------------------------------------------------------
      ! 
      ! 4. Add the current field to the suitable 3d-field.
      !
      CALL add_windfield_to_3d_windfield(windfield,stack%winds_3d(ind)%fp)

    END DO loop_over_2d_windfields


    ! -------------------------------------------------------------
    ! 
    ! 5. In each 3d-windfield organize the fields vertically.
    !
    i = 1
    DO while(i <= stack%windstorepointer_3d)

      wind_3d => stack%winds_3d(i)%fp

      !
      ! We have to check if this very field_3d is really 3d, which
      ! means that it contains more than one 2d field.
      ! Reason - quantities are considered 3d if they can be defined
      ! at several levels. But this is not necessary that in specific
      ! case they are really 3d
      !
      if(fu_number_of_windfields(wind_3d) < 2)then
        !
        ! Squeeze the list of fields
        !
        if(i < stack%windstorepointer_3d)then
          do j=i, stack%windstorepointer_3d-1

            wind_3d => stack%winds_3d(j)%fp
            call set_3d_windfield_empty(wind_3d) ! Kill j-th 3d field

            wind_3d => stack%winds_3d(j+1)%fp
            do k=1,fu_number_of_windfields(wind_3d)
              call add_windfield_to_3d_windfield( &
                        & fu_windfield_from_3d_wind(wind_3d,k,.true.), & ! skip order
                        & stack%winds_3d(j)%fp)
            end do

          end do
        end if
        call set_3d_windfield_empty(stack%winds_3d(stack%windstorepointer_3d)%fp)
        stack%windstorepointer_3d = stack%windstorepointer_3d - 1
        if(i > stack%windstorepointer_3d)exit
      endif

      CALL organize_windfields_vertically(stack%winds_3d(i)%fp)
      IF (error) RETURN

      i=i+1

    END DO


    ! -------------------------------------------------------------
    ! 
    ! 6. Add surface pressure fields to those 3d-fields that contain hybrid-level data
    !
    DO i = 1, stack%windstorepointer_3d

      wind_3d => stack%winds_3d(i)%fp
      IF (fu_leveltype(wind_3d) == hybrid) THEN

   CALL find_matching_surface_pre_wind(wind_3d, stack)

   IF (error) RETURN
      END IF
    END DO


  CONTAINS

    SUBROUTINE search_match()

      ! Description:
      ! Finds fo a windfield the corresponding 3d-windfield from a set of 3d
      ! windfields.
      !
      ! Method:
      ! Check met_src, quantity an time

      IMPLICIT NONE

      ! Local declarations:
      INTEGER :: i
      type(meteo_data_source) :: met_src
      TYPE(silja_time) :: time
      TYPE(silja_3d_windfield), POINTER :: windfield_3d

      found = .false.
      met_src = fu_met_src(windfield)
      time = fu_valid_time(windfield)
      DO i = 1, SIZE(stack%winds_3d)
        windfield_3d => stack%winds_3d(i)%fp
        IF (.NOT.defined(windfield_3d)) EXIT
        IF (time /= fu_valid_time(windfield_3d)) CYCLE
        IF (.not.met_src == fu_met_src(windfield_3d)) CYCLE
        found = .true.
        ind = i
        RETURN
      END DO

    END SUBROUTINE search_match

    ! ***************************************************************


    SUBROUTINE find_matching_surface_pre_wind(wind_3d, stack)

      ! Description:
      ! For the unfinished 3d-windfield finds its corresponding surface
      ! pressure field. This is used for hybrid level data, where every
      ! 3d-windfield  contains also a pointer to to matching surface
      ! pressure field.
      !
      ! Language: ANSI Fortran 90
      !
      ! Author: Mika Salonoja, FMI
      ! 
      IMPLICIT NONE

      ! Imported parameters with intent INOUT or POINTER:
      TYPE(silja_3d_windfield), POINTER :: wind_3d
      TYPE(silja_stack), POINTER :: stack

      ! Local declarations:
      INTEGER :: i
      TYPE(silja_field), POINTER :: field
  
      loop_over_2d_scalar_fields: DO i = 1, SIZE(stack%fields)

        field => stack%fields(i)%fp
        IF (.NOT.defined(field)) EXIT
        IF (fu_quantity(field) /= ground_pressure_flag) CYCLE
!        IF (.NOT.(fu_level(field) == ground_level)) CYCLE  !ehto ei tayty ?!
        IF (error) RETURN
        ! --------------------------------------
        !
        ! Ground pressure found, check the rest: 
        ! Time:
        !
        IF (.NOT.stack%one_valid_time) THEN
          IF (.NOT.(fu_valid_time(field) == fu_valid_time(wind_3d))) CYCLE
        END IF
        ! met_src of data:
!        IF (.NOT.(fu_met_src(field) == fu_met_src(wind_3d))) CYCLE

        ! Correct field found:
        CALL add_surf_pre_to_3d_wind(wind_3d, field)

        RETURN


      END DO loop_over_2d_scalar_fields

      CALL msg_warning('no surface pressure found for hybrid 3D','find_matching_surface_pre_wind')    

    END SUBROUTINE find_matching_surface_pre_wind

  END SUBROUTINE arrange_windfields_to_3d





  ! ***************************************************************


  INTEGER FUNCTION fu_stacksize(stack)
    
    ! Description:
    ! Returns the number of real numbers in a stack.
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
    TYPE(silja_stack), intent(in) :: stack
    !
    ! Local declarations:
    INTEGER :: i    
    TYPE(silja_field), POINTER :: field

    fu_stacksize = 0

    IF (.NOT.defined(stack)) RETURN

    DO i = 1, SIZE(stack%fields)

      field => stack%fields(i)%fp

      IF (.NOT.defined(field)) EXIT

      fu_stacksize = fu_stacksize + fu_size(field)

    END DO

  END FUNCTION fu_stacksize



  ! ***************************************************************


  FUNCTION fu_get_wind_3d_field(stack, iFieldNbr) RESULT(wind_3d_field)
  !
  ! Just returns the wind_3d_field pointed by the iFieldNbr
  !
  ! Code owner Mikhail Sofiev, FMI
  
    IMPLICIT NONE

    ! Imported parameters, intent IN
    TYPE(silja_stack), INTENT(in), TARGET :: stack
    INTEGER, INTENT(in) :: iFieldNbr

    ! Result:
    TYPE(silja_3d_windfield), POINTER :: wind_3d_field

    ! Local declarations
    TYPE(silja_stack), POINTER :: stackpointer

    stackpointer => stack

    IF(.not.defined(stackpointer))THEN
       NULLIFY(wind_3d_field)
    ELSE
      wind_3d_field => stackpointer%winds_3d(iFieldNbr)%fp
    END IF

  END FUNCTION fu_get_wind_3d_field


  ! ***************************************************************


  FUNCTION fu_get_wind_2d_field(stack, iFieldNbr) RESULT(wind_2d_field)
  !
  ! Just returns the windfield pointed by the iFieldNbr
  !
  ! Code owner Mikhail Sofiev, FMI
  
    IMPLICIT NONE

    ! Imported with intent IN
    TYPE(silja_stack), INTENT(in), TARGET :: stack
    INTEGER, INTENT(in) :: iFieldNbr

    ! result
    TYPE(silja_windfield), POINTER :: wind_2d_field

    ! Local declarations
    TYPE(silja_stack), POINTER :: stackpointer

    stackpointer => stack

    IF(.not.defined(stackpointer))THEN
       NULLIFY(wind_2d_field)
    ELSE
      wind_2d_field => stackpointer%winds(iFieldNbr)%fp
    END IF

  END FUNCTION fu_get_wind_2d_field



  ! ***************************************************************


  FUNCTION fu_get_3d_field(stack, iFieldNbr) RESULT(field_3d)
  !
  ! Just returns the 3d_field pointed by the iFieldNbr
  !
  ! Code owner Mikhail Sofiev, FMI
  
  IMPLICIT NONE

  ! Imported IN
  TYPE(silja_stack), INTENT(In), TARGET :: stack
  INTEGER, INTENT(in) :: iFieldNbr

  !result
  TYPE(silja_3d_field), POINTER :: field_3d

  ! Local declarations
  TYPE(silja_stack), POINTER :: stackpointer
!  LOGICAL :: stack_defined

    stackpointer => stack

    IF(.not.defined(stackpointer)) THEN
       NULLIFY(field_3d)
    ELSE
      field_3d => stackpointer%fields_3d(iFieldNbr)%fp
    END IF

  END FUNCTION fu_get_3d_field



  ! ***************************************************************


  FUNCTION fu_get_2d_field(stack, iFieldNbr) RESULT(field_2d)
  !
  ! Just returns the field pointed by the iFieldNbr
  !
  ! Code owner Mikhail Sofiev, FMI
  
  IMPLICIT NONE

  ! Imported IN
  TYPE(silja_stack), INTENT(in), TARGET :: stack
  INTEGER, INTENT(in) :: iFieldNbr

  !result
  TYPE(silja_field), POINTER :: field_2d

  ! Local declarations
  TYPE(silja_stack), POINTER :: stackpointer

    stackpointer => stack

    NULLIFY(field_2d)

    IF(defined(stackpointer))THEN
      IF(fu_field_defined(stackpointer%fields(iFieldNbr)%fp))THEN
        field_2d => stackpointer%fields(iFieldNbr)%fp
      END IF
    END IF

!    IF(.not.defined(stackpointer).or. &
!     & fu_size(stackpointer%fields(iFieldNbr)%fp) <= 0)THEN
!       NULLIFY(field_2d)
!    ELSE
!      field_2d => stackpointer%fields(iFieldNbr)%fp
!    END IF

  END FUNCTION fu_get_2d_field





  ! ***************************************************************

  ! ***************************************************************
  ! 
  !
  !                 TESTS
  !
  !
  ! ***************************************************************
  ! ***************************************************************

  SUBROUTINE print_stack_report(stack, ifShort_)

    ! Description:
    ! Prints contents of stack.
    ! 
    ! Language: ANSI Fortran 90
    !
    ! Author: Mika Salonoja, FMI
    ! 
    IMPLICIT NONE

    ! Imported parameters with intent IN:
    TYPE(silja_stack), intent(in) :: stack
    logical, intent(in), optional :: ifShort_

    ! Local declarations:
    INTEGER :: i, co
    TYPE(silja_field), POINTER :: field
    TYPE(silja_windfield), POINTER :: windfield
    TYPE(silja_3d_field), POINTER :: field_3d
    TYPE(silja_3d_windfield), POINTER :: wind_3d
    INTEGER :: fc, wfc, f3c, wf3c
    logical :: ifShort

    IF (.NOT.defined(stack)) THEN
      call msg(' Undefined stack')
      RETURN
    else
      call msg(fu_connect_strings('%%%%%%%%%%%%_report of stack:', &
                                & stack%name, &
                                & '_%%%%%%%%%%%%%%'))
    end if

    if(present(ifShort_))then
      ifShort = ifShort_
    else
      ifShort = .false.
    endif

    call msg(' Size of stack now ', fu_stacksize(stack))

    call msg(' ********** SCALAR FIELDS ***********')
    fc = 0
    DO i = 1, min(stack%storepointer,size(stack%fields))
      field => stack%fields(i)%fp
      IF (defined(field)) THEN
        fc = fc + 1
        if(ifShort)then
          call msg(fu_quantity_string(fu_quantity(field)))
        else
          CALL report(field)
        endif
        IF (error) RETURN
      END IF
    END DO


    wfc = 0
    IF (ASSOCIATED(stack%winds)) THEN
      call msg(' ********** WINDFIELDS ***********')
      DO i = 1, min(stack%windstorepointer,size(stack%winds))
        windfield => stack%winds(i)%fp
        IF (defined(windfield)) THEN
          wfc = wfc + 1
          if(ifShort)then
            call msg('Wind field')
          else
            CALL report(windfield)
          endif
          IF (error) RETURN
        END IF
      END DO
    ELSE
      call msg(' ********** NO WINDFIELDS ***********')
    END IF


    f3c = 0
    IF (ASSOCIATED(stack%fields_3d)) THEN
      call msg(' ********** 3D-FIELDS ***********, total:', stack%storepointer_3d)
      DO i = 1, min(stack%storepointer_3d,size(stack%fields_3d))
        field_3d => stack%fields_3d(i)%fp
        IF (defined(field_3d)) THEN
          f3c = f3c + 1
          if(ifShort)then
            call msg(fu_quantity_string(fu_quantity(field_3d))+ "_ Field No:"+fu_str(i))
          else
            call msg("Field No:"+fu_str(i))
            CALL report(field_3d)
            call msg("")
          endif
          IF (error) RETURN
        END IF
      END DO
    ELSE
      call msg(' ********** NO 3D-FIELDS ***********')
    END IF


    wf3c = 0
    IF (ASSOCIATED(stack%winds_3d)) THEN
      call msg(' ********** 3D-WINDFIELDS ***********')
      DO i = 1, min(stack%windstorepointer_3d,size(stack%winds_3d))
        wind_3d => stack%winds_3d(i)%fp
        IF (defined(wind_3d)) THEN
          wf3c = wf3c + 1
          if(ifShort)then
            call msg('3D wind field')
          else
            CALL report(wind_3d)
          endif
          IF (error) RETURN
        END IF
      END DO
    ELSE
      call msg(' ********** NO 3D-WINDFIELDS ***********')
    END IF


    call msg('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
    call msg(' fields in total. ', fc)
    call msg(' windfields in total. ', wfc)
    call msg(' 3D-fields in total. ', f3c)
    call msg(' 3D-windfields in total. ', wf3c)
    call msg('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')

  END SUBROUTINE print_stack_report

END MODULE stacks
