MODULE field_buffer
!
! This module contains necessary structures and routines for interface
! between the meteorological part of SILAM and dispersion models themselves
! The main idea is: the structures defined here are known for the model,
! while the meteorological structures are partially known for routines
! here.
! Tasks performed:
! 1. Collecting input needs from the model and then initialising
!    the meteo structures (former PASI_NWP_ini from particle_models)
! 2. At every meteo time step (when the next meteodata should be taken)
!    arranging the collection of pointers to the requested fields.
! 3. 2D fields are worth to be time-interpolated immediately here.
!    They are small enough and the number of operations is less than
!    that of invidual interpolation for each particle
!
! Author: Mikhail Sofiev, FMI, mikhail.sofiev@fmi.fi
!
! All units: SI
! 
! All degrees: real numbers in degrees (degrees and parts per
! hundred/parts per thousand , NOT in degrees and minutes
! 
! Language: ANSI standard Fortran 90
! 

use input_analysis
!use physiographies
!$use omp_lib

implicit none

!  public functions

PUBLIC meteo_init
public init_supplementary_market
public init_supplementary_data_buffer
public init_data_buffer
PUBLIC arrange_buffer
!!PUBLIC set_permanent_buffer_pointers ! arrange_buffer with now=time_missing
                                       ! does this
public fu_dimension ! 2d or 4d field
public fu_2d_field   ! Can be called instead of digging in the structures
public fu_4d_field   ! Can be called instead of digging in the structures
public fu_get_value  ! Get a single value out of meteo buffer - both 4d and 2d
public fu_vert_index  ! Get an index in the certain grid for "closest-cell" interpolation
public interp_field_3d_2_field_3d ! Interpolation 3D for field_3d class
public refine_interp_vert_coefs_v2

public wind_from_buffer_4d

PUBLIC fu_4d_interpolation
public make_2d_from_4d_field
PUBLIC fu_index
PUBLIC fu_make_position_pressure
PUBLIC fu_pressure_index
public fu_vertical_index

! Public functions for the interpolation structures (common, vertical and horizontal specific)
!
public fu_nCoefs      ! common functions start
public fu_interpType
public get_coefs
public fu_ifMeteoDependentInterp ! if the vertical interpolation is meteo-dependent
public fu_vertFrom
public fu_vertTo
public fu_ifMeteoGrd
public fu_grid
public fu_vertical_interp_struct
public get_vertical_range_vert_interp
public set_missing
public report
public refine_all_vert_interp_coefs
public column_from_buffer
public surf_from_buffer


! Private functions

private make_present_variable
PRIVATE fu_index_for_quantity
PRIVATE fu_index_for_quantity_in_buf
private fu_index_for_variable_in_buf
PRIVATE fu_idx4q_in_buf_with_ptr2d
private fu_idx4var_in_buf_with_ptr2d
private fu_dimension_of_meteofield
private fu_4d_interp_to_position
private fu_4d_interp_to_coord
private fu_2d_field_from_buffer
private fu_4d_field_from_buffer
private fu_get_value_from_buffer
private fu_get_value_from_bfr_interp4d
private fu_get_value_from_fldp2d_int3d
private fu_get_value_from_fld2d_int
private fu_get_value_from_bfr_interp
private fu_ind_vert_from_bfr_interp2d
private set_buffer_pointers
private fu_pressure_index_position
private fu_pressure_index_coord
private interp_fld3d_field_2_fld3d
private interp_fld3d_data_2_fld3d

! Private functions of the interpolation structures

private fu_vertFrom_from_interp_struct  ! Vertical interpolation starts
private fu_vertTo_from_interp_struct
private fu_nCoefs_interp_str_vert
private fu_interpType_interp_str_vert
private get_coefs_interp_str_vert
private get_coefs_interp_str_cell_vert
private fu_grid_of_vertInterpStruct
private set_missing_vert_interp_struct

private report_buffer


INTERFACE fu_index
  MODULE PROCEDURE fu_index_for_quantity
  module procedure fu_index_for_quantity_in_buf
  module procedure fu_index_for_variable_in_buf
  module procedure fu_idx4q_in_buf_with_ptr2d
  module procedure fu_idx4var_in_buf_with_ptr2d
END INTERFACE

interface fu_dimension
  module procedure fu_dimension_of_meteofield
end interface

interface fu_2d_field
  module procedure fu_2d_field_from_buffer
end interface

interface fu_4d_field
  module procedure fu_4d_field_from_buffer
end interface

interface fu_get_value
  module procedure fu_get_value_from_buffer
  module procedure fu_get_value_from_bfr_interp4d
  module procedure fu_get_value_from_fldp2d_int3d
  module procedure fu_get_value_from_fld2d_int
  module procedure fu_get_value_from_bfr_interp
end interface

interface fu_vert_index
  module procedure fu_ind_vert_from_bfr_interp2d
end interface

interface fu_pressure_index
   module procedure fu_pressure_index_position
   module procedure fu_pressure_index_coord
end interface

interface fu_4d_interpolation
  module procedure fu_4d_interp_to_position
  module procedure fu_4d_interp_to_coord
end interface

interface report
  module procedure report_buffer
end interface

interface interp_field_3d_2_field_3d
  module procedure interp_fld3d_field_2_fld3d
  module procedure interp_fld3d_data_2_fld3d
end interface

!--------------------- Interface for the interpolation structures

interface fu_vertFrom
  module procedure fu_vertFrom_from_interp_struct
end interface 

interface fu_vertTo
  module procedure fu_vertTo_from_interp_struct
end interface 

interface fu_nCoefs
  module procedure fu_nCoefs_interp_str_vert
end interface 

interface fu_interpType
  module procedure fu_interpType_interp_str_vert
end interface 

interface get_coefs
  module procedure get_coefs_interp_str_vert
  module procedure get_coefs_interp_str_cell_vert
end interface 


interface fu_grid
  module procedure fu_grid_of_vertInterpStruct
end interface 

interface set_missing
  module procedure set_missing_vert_interp_struct
end interface

! procedures for the to-point interpolation


  ! Public types defined here:

  ! For direct operation with 2d fields
  TYPE field_data_ptr
    REAL, DIMENSION(:), POINTER :: ptr => null()
    type(silja_field_id), pointer :: idPtr =>null()
    logical :: ifReady
  END TYPE field_data_ptr

  ! For direct operation with past, present and future 2d fields
  TYPE field_2d_data_ptr
    type(field_data_ptr) :: past, present, future
  END TYPE field_2d_data_ptr

  ! For direct operation with 3d fields
  TYPE field_3d_data_ptr
    TYPE(field_data_ptr), DIMENSION(max_levels) :: p2d
    type(silja_3d_field), pointer  :: field3d => null() !Somewhat redundant, but useful
                                              ! to get 3d f/ield via buffer...
  END TYPE field_3d_data_ptr

  ! For direct operation with past and future of 3d fields
  ! Present field is not necessarily filled-in.
  TYPE field_4d_data_ptr
    TYPE(field_3d_data_ptr) :: past, future, present
  END TYPE field_4d_data_ptr


  !-------------------------------------------------
  ! 
  ! Field_buffer is COMPLETELY OPEN for everybody. The idea:
  ! The pointers to the fields exactly follow the order of quantities 
  ! requested by the model, which are stored in buffer_quantities. 
  ! Important: there is no way to combine 2-d and 3-d field pointers in
  ! one array, so they are split to two arrays, BUT: if the n-th quantity
  ! is 2d then this pointer in p3d is NULL (skipped). So, for each index
  ! only one of arrays p2d and p3d is filled.
  !
  ! Components of the meteo_buffer type:
  ! p2d/p3d - array of pointers to 2d/3d fields
  ! grd_shift - grid shift of the corresponding quantity
  ! The order of quantities is synchronised with buffer_quantities
  !
 TYPE Tfield_buffer
    TYPE(field_2d_data_ptr),DIMENSION(:),POINTER :: p2d => null()
    TYPE(field_4d_data_ptr),DIMENSION(:),POINTER :: p4d => null()
    real :: weight_past
    integer :: nbr_of_levels
    integer, dimension(:), pointer :: buffer_quantities => null()
    type(silam_species), dimension(:), pointer :: buffer_species => null()
    logical, dimension(:), pointer :: ifPointerSet  => null()
    type(silja_time) :: time_past, time_present, time_future
    real, dimension(:,:,:,:,:), pointer :: wings_N => null(), wings_S => null()
    real, dimension(:,:,:,:,:), pointer :: wings_E => null(), wings_W => null()
    ! iPast/iFuture/iRealTime, iMass/iFlux, iLev, ix, iy, iBoundary
    logical :: ifMassFluxBottomUp = .True. !  Quick and dirty
                                           !How the wind was diagnosed: needed for pole vertical advection
  END TYPE Tfield_buffer


  !---------------------------------------------------------------------
  !
  ! Vertial structues of different types can be interpolated between each other
  ! In principle, if the meteodata
  ! are not available, interpolation is strictly one-dimensional. But if meteodata are
  ! actual, there is a map of vertical interpolation coefs, which has to be computed
  ! separately
  !
  type TVertInterpOneCell
    !
    ! Dimensions are all the same and depend on type of interpolation
    !
    integer, dimension(:), pointer :: indLev ! In the vertFrom
    real, dimension(:), pointer :: weight, weight_Z ! of the specific index
  end type TVertInterpOneCell

  type TVertInterpCells
    !
    ! Dimensions are all the same and depend on type of interpolation
    !
    integer, dimension(:,:,:,:), pointer :: indLev ! In the vertFrom
    real, dimension(:,:,:,:), pointer :: weight, weight_Z ! of the specific index
  end type TVertInterpCells

  type TVertInterpStruct
!    private
    character(len=40) :: chName
    type(silam_vertical) :: vertFrom, vertTo
    type(silja_grid) :: grid  ! Grid in which the structure is to be defined
    logical :: ifMeteoGrid    ! whether this grid == meteo_grid
    type(silja_time) :: lastCoefUpdate
    type(silja_interval) :: recommended_update_interval
    integer :: interp_type, nCoefs, iLevStartFrom, iLevEndFrom, iLevStartTo, iLevEndTo
    !integer, dimension(:), pointer :: metQuantities ! needed for fine-tuning interpolation
    integer, dimension(:,:,:,:), pointer :: indLev => null()     ! (nCoefs,nx,ny,nLevsTo)
    real, dimension(:,:,:,:), pointer :: weight =>null(), weight_Z =>null() ! (nCoefs,nx,ny,nLevsTo)
  end type TVertInterpStruct

  type TVertInterpStructPtr
    type(TVertInterpStruct), pointer :: ptr =>null()
  end type TVertInterpStructPtr

  type (TVertInterpStruct), public, parameter :: VertInterpStruct_missing_par = &
        & TVertInterpStruct('',vertical_missing, vertical_missing, grid_missing,&
                         & .false., time_missing, interval_missing, &
                         & int_missing, int_missing, int_missing, int_missing,&
                         & int_missing, int_missing, null(), null(), null())

  type (TVertInterpStruct), public, target :: VertInterpStruct_missing = VertInterpStruct_missing_par
  !
  ! The run will have a pool of interpolation structures, where everyone 
  ! will be free to find the favorite one
  !
  type TVertInterpPool
    type(TVertInterpStructPtr), dimension(40) :: pVIS  ! vertical interpolation structure pointer
    integer :: nVIS = 0
  end type TVertInterpPool
  private TVertInterpPool
  
  type(TVertInterpPool), private, save :: VertInterpPool  ! the bunch of interpolation structures



  logical, private, parameter :: ifDebug = .false.

CONTAINS


  !********************************************************

  SUBROUTINE meteo_init(wdr, &
                      & input_shopLst, full_shopLst, &
                      & input_shopLst_st, full_static_shopLst, &
                      & ifForward, &
                      & out_grid, &
                      & input_quantities, &
                      & meteoMarketPtr, &
                      & nTimeNodesNeeded)
    !
    ! Initializes supermarket and does some other nwp-related
    ! initialization. Sets original value for shopping list
    !
    IMPLICIT NONE

    ! Imported parameters with intent IN:
    LOGICAL, INTENT(in) :: ifForward
    type(silja_grid), intent(in) :: out_grid ! "output" grid to be covered by meteo data
    INTEGER, DIMENSION(:), intent(out) :: input_quantities 
                ! quantities to be used in buffer. Regardless of st/dyn
                ! requested by the model
    integer, intent(in) :: nTimeNodesNeeded

    ! Imported parameters with intent INOUT or OUT:
    TYPE(silja_wdr), intent(inout), target :: wdr
    TYPE(silja_shopping_list), INTENT(inout) :: input_shopLst, input_shopLst_st
    TYPE(silja_shopping_list), INTENT(out) :: full_shopLst, full_static_shopLst
    type(mini_market_of_stacks), pointer :: meteoMarketPtr

    ! Local declarations:
    INTEGER, DIMENSION(max_quantities) :: quantities, qTmp, qTmpStatic, requests, requests_st
    INTEGER :: i, j, iLev, nQ, nQs, fs, iCount, iTempl, iTmp, q
    type(silam_sp), dimension(:), pointer :: fnames
    type(Tinput_content) :: InputContent
    type(wdr_ptr), dimension(:), pointer :: wdrar

    InputContent = input_content_missing

    nullify(fnames)

    ! --------------------------------------------------
    !
    ! List of quantities requested by the model has to be saved. 
    ! Meteo-specific stuff like ABL-generating process may require
    ! some more derived quantites, so the final step 
    ! is to check what the derived_quantities module says.
    !
    ! - full_shopping_list contains ALL quantities needed for the run.
    ! - input_shopping_list will be made so to contain ONLY THOSE to shop
    !   from the input file
    !

    call msg("input_shopLst dyn")
    call report(input_shopLst)
    call msg("input_shopLst st")
    call report(input_shopLst_st)
!    met_buf%nbr_of_levels = 0
    input_quantities(:) = int_missing
    requests = 0
    !---------------------------------------------------
    !
    ! So far, input_shopping_list contains not really input quantities but
    ! rather what model and output required, so most of them are derived ones.
    ! Let's copy them to input_quantities - they will later be
    ! referred by the meteobuffer and by the model.
    !

    !----------------------------------------------------
    !
    ! Now let's create the full_shopping_list, which contains all imaginable
    ! quantities, potentially needed for creation of the input_quantities
    ! This wide list will be asked from the input files for further analysis of
    ! what is available
    !
    full_shopLst = input_shopLst
    full_static_shopLst = input_shopLst_st
    iCount =   fu_nbr_of_quantities(full_shopLst)  &
              &+ fu_nbr_of_quantities(full_static_shopLst)
    DO WHILE (iCount > 0) ! Counter of added quantities
      quantities = fu_quantities(full_shopLst) 
      DO i=1,size(quantities)
        IF(quantities(i) == int_missing) EXIT
        ! Check dynamic quantities
        call quantities_for_derived_one(quantities(i),  fu_request(full_shopLst,i), .false., &
                                      & qTmp, qTmpStatic, requests, requests_st, wdr)
        call add_shopping_quantities(full_shopLst, qTmp, requests)
        call add_shopping_quantities(full_static_shopLst, qTmpStatic, requests_st)
      END DO
      quantities = fu_quantities(full_static_shopLst) 
      DO i=1,size(quantities)
        IF(quantities(i) == int_missing) EXIT
        ! Check only realtime quantities
        if (.not. fu_realtime_quantity(quantities(i))) cycle
        call quantities_for_derived_one(quantities(i), fu_request(full_static_shopLst,i), .true., &
                                      & qTmp, qTmpStatic, requests, requests_st, wdr)
        call add_shopping_quantities(full_shopLst, qTmp, requests)
        call add_shopping_quantities(full_static_shopLst, qTmpStatic, requests_st)
      END DO
      !
      ! If the number of quantities has changed since the last run - reset it and repeat 
      ! the cycle. Otherwise, set the exit condition to true
      !
      iTmp = fu_nbr_of_quantities(full_shopLst)  + fu_nbr_of_quantities(full_static_shopLst)
      if(iCount /= iTmp)then
        iCount = iTmp
      else
        iCount = 0
      endif
    END DO
    if(error)return

    quantities = fu_quantities(full_shopLst) ! restore for the sake of further potential use

!    call msg("Quantities to request from input")
!    call report_list_of_quantities(quantities, "full_shopLst")
!    call report_list_of_quantities(fu_quantities(full_static_shopLst), "full_st_shopLst")

    CALL fix_shopping_levels (full_shopLst, level_missing, level_missing)
    CALL fix_shopping_levels (full_static_shopLst, level_missing, level_missing)
    IF (error) RETURN

    !----------------------------------------
    !
    ! Check what is available in the dynamic input files. In particular, 
    ! meteo_grid is selected here. 
    !
    do iTempl = 1, max_met_files
      if(fu_fname_template(wdr, iTempl) == template_missing) exit
      !
      ! Depending on the format of the input data (file or field), we may need
      ! to expand the file names or simply copy the template output collection
      !
      if(fu_if_input_file(fu_file_format(wdr,iTempl)))then
        call FNm_from_single_template(fu_fname_template(wdr,iTempl),&
                                    & fu_closest_obstime(fu_start_time(wdr) + fu_meteo_time_shift(wdr), &
                                                       & back_and_forwards, &
                                                       & fu_obstime_interval(wdr)), &
                                    & fnames, &
                                    & ifAdd = .false., &
                                    & ifStrict = .false., & 
                                    & ifAllowZeroFcLen = fu_if_zero_fc_len_allowed(wdr), &
                                    & ifWait = .false.)
      else
        call enlarge_array(fnames,1)
        if(error)return
        fnames(1)%sp = fu_collection(fu_fname_template(wdr, iTempl))
        if(size(fnames) > 1) fnames(2)%sp = ''
      endif
      if(error)return

      if((.not. associated(fnames)) .or. size(fnames)<1)then
        call set_error('Seems to be empty input template','meteo_init')
        return
      endif
      if(fnames(1)%sp == '')then
        call set_error('Failed to find any input file from template','meteo_init')
        return
      end if
      !
      ! Read the input file asking all variables stored in full_shopping_list
      !
      do i=1,size(fnames)
        if(fnames(i)%sp == '')exit
        call read_input_content(fnames(i)%sp, wdr, fu_file_format(wdr,iTempl), full_shopLst, InputContent)
        if(error)return
      end do
    end do ! dynamic templates
    !-------------------------------------------------------------------------
    !
    ! Dynamic GRIB content structure is now filled-in, so we can: 
    !   - analyse the GRIB content,
    !   - select the meteo_grid and meteo_vertical.
    !   - re-define the input_shopping_list in accordance with selected data fields
    ! The input_quantities must be filled-in
    !
    !-------------------------------------------------------------------------
    !
    call msg("analyse_input_content got list for  shopping:")
    call report(input_shopLst)
   
   ! Somewhat hackish...
   ! Here we add dynamic quantities needed for real-time quantities to the input_shopLst
   ! Pretend that they were requested...
    quantities = fu_quantities(full_static_shopLst) 
    do i=1, size(quantities)
      if((quantities(i) > 0)  .and. fu_realtime_quantity(quantities(i)))then
        call quantities_for_derived_one(quantities(i), fu_request(full_static_shopLst,i), .true., &
                                      & qTmp, qTmpStatic, requests, requests_st, wdr) 
        call add_shopping_quantities(input_shopLst, qTmp, requests)
      endif
    end do
   
    call analyse_input_content(InputContent, &
                             & out_grid, &
                             & input_shopLst, &
                             & wdr)
    if(error)return


    !----------------------------------------------------------------
    !
    ! The next step is - to re-create the full_shopping_list so that it
    ! contains exactly those variables, which enable the way from 
    ! available quantities in the Grib files to input_shopping_list
    !
    ! Trick: check_quantity assumes truly static quantities available
    !

  !call msg("analyse_input_content returned shopping list:")
  !call report(fu_shopping_list(InputContent))


    call set_deriving_way(fu_shopping_list(InputContent), &  ! all what is available
                        & input_shopLst, full_shopLst, &     ! requested by model, final list (out)
                        & input_shopLst_st, full_static_shopLst, & ! requested by model, final list (out)
                        & wdr)                               ! metadata
    if(error)return

!  call msg("set_deriving_way returned full_shopLst:")
!  call report(full_shopLst)
!  call msg("set_deriving_way returned full_shopLst:")
!  call report(full_static_shopLst)

    iCount=1
    do i=1, size(fu_quantities(full_shopLst))
      q = fu_quantity(full_shopLst,i)
      if( fu_request(full_shopLst, i) > 0  .and. all(input_quantities(1:iCount) /= q) )then
        if (ifCanBeSkippedInDerivation(q)) then
          if( .not. fu_quantity_in_list(q, input_shopLst)) cycle
        endif
        input_quantities(iCount) = q
        iCount = iCount + 1
      endif
    end do
    input_quantities(iCount) = int_missing
!call report_list_of_quantities(input_quantities,"After full_shopLst")

    ! No Intelligent handling of derivation
    input_shopLst_st = full_static_shopLst
!

    !----------------------------------------------------------------
    !
    ! In principle, the same deriving method approach should be applied to static 
    ! fields. However, it will reqire another InputContent and once again
    ! full-range analysis. This seems to be too much for a few permanent fields.
    ! So far, it should be enough to get it straight to the input field.
    ! Physiography will require them or, if needed, will make some derivations.
    ! Here we just add proper quantities to the input_quantities, so that the
    ! model gets what it needs from permanent fields as well.
    !

    do i=1,size(fu_quantities(input_shopLst_st))
      if(fu_request(input_shopLst_st, i) > 0 .and. &
       & .not. any(input_quantities(1:iCount) == fu_quantity(input_shopLst_st,i)))then
        input_quantities(iCount) = fu_quantity(input_shopLst_st, i)
        iCount = iCount + 1
      endif
    end do
    input_quantities(iCount) = int_missing
!call report_list_of_quantities(input_quantities,"After full_shopLst_st")


   !
   ! Now let's expand the static shopping list to include all quantities that might be needed
   !
   iCount = 1
    DO WHILE (iCount > 0) ! Counter of added quantities
      quantities = fu_quantities(input_shopLst_st) 
      DO i=1,size(quantities)
        IF(quantities(i) == int_missing) EXIT

        call quantities_for_derived_one(quantities(i), fu_request(input_shopLst_st,i), .true., &
                                      & qTmp, qTmpStatic, requests, requests_st, wdr) 
        call add_shopping_quantities(full_shopLst, qTmp, requests)
        call add_shopping_quantities(full_static_shopLst, qTmpStatic, requests_st)
      END DO
      !
      ! If the number of quantities has changed since the last cycle - reset it and repeat 
      ! the cycle. Otherwise, set the exit condition to true
      !
      if(iCount /= fu_nbr_of_quantities(full_static_shopLst))then
        iCount = fu_nbr_of_quantities(full_static_shopLst)
      else
        iCount = 0
      endif
    END DO
    if(error)return

    !----------------------------------------------------------------
    !
    ! Now the input_shopping_list converts to what is should be - a list of
    ! quantities shopped from the input meteofiles
    !
    input_shopLst= fu_shopping_list(InputContent)
    if(error)return


    call msg_test('')
    call msg_test('-----------The input shopping list is:-----------')
    if(test_messages) call report(input_shopLst)
    call msg_test('')
    call msg_test('')
    call msg_test('')

    call msg_test('')
    call msg_test('-----------The full shopping list is:-----------')
    if(test_messages) call report(full_shopLst)
    call msg_test('')

    !
    ! Create the maximum possible system_grid, which is still covered
    ! with meteorological fields. It is not really the system grid yet
    ! because there may be some memory reduction attempt, if output grid
    ! and sources are all within some smaller area.
    !
    meteo_grid = fu_meteo_grid(InputContent)
    if(error)return

    call msg('------------------------------------------------------')
    call msg('                 Meteorological grid: ')
    CALL report(meteo_grid)
    call msg('')

    !-----------------------------------------------------------
    !
    ! system_pole can be made from the system_grid
    !
    IF(fu_gridtype(meteo_grid) == lonlat)THEN
      meteo_pole = fu_pole(meteo_grid)
    elseif(fu_gridtype(meteo_grid) == anygrid)THEN
      meteo_pole = pole_geographical
    ELSE
      CALL set_error('Strange grid for system one','meteo_init')
      RETURN
    END IF


    ! --------------------------------------------------
    !
    ! Set the meteo_vertical and vertical level boundaries for the structure
    !
    ! Define the bottom pressure type - either mean sea level pressure
    ! for pressure vertical, or ground level pressure for the terrain-following 
    ! systems
    !
    meteo_vertical = fu_meteo_vertical(InputContent)
    meteo_verticalPtr => meteo_vertical
    call vertical_parameters_ab(meteo_vertical, nz_meteo, a_met, b_met)

    select case (fu_leveltype(meteo_vertical))
      case (constant_pressure, constant_altitude)  ! Pressure co-ordinates
        call set_top_level(wdr, pr_level_10hpa) ! Quite high indeed
        call set_bottom_level(wdr, mean_sea_level)

      case (constant_height, sigma_level, hybrid, layer_btw_2_hybrid) ! Terrain-following co-ordinates
        call set_top_level(wdr, level_missing) ! Accept all
        call set_bottom_level(wdr, level_missing)   ! Accept all

      case default
        call set_error('Unknown system_level type','meteo_init')
    call msg('------------------------------------------------------')
    call msg('                  METEO VERTICAL: ')
    CALL report(meteo_vertical, .true.)
    call msg('------------------------------------------------------')
    call msg('')
        call set_error('Unknown system_level type bis','meteo_init')
        return
    end select

    CALL fix_shopping_levels (full_shopLst, fu_bottom_level(wdr), fu_top_level(wdr))
    IF (error) RETURN

    call msg('------------------------------------------------------')
    call msg('                  METEO VERTICAL: ')
    CALL report(meteo_vertical, .true.)
    call msg('------------------------------------------------------')
    call msg('')


    !------------------------------------------------------------
    !
    ! Now we can also check the number of sources in the current GRIB file set.
    ! Correct place to store them - wdr, of course. Note: the data sources are
    ! set in order so that the first source is the "most important" one.
    !
    ! ATTENTION. There is a dangerous but necessary trick overruling this: if in wdr
    ! the field ifDisregard MDS == .true., there will be only one MDS: met_src_missing.
    ! The reason is simple: ECMWF changes models from time to time and each time sets 
    ! new MDS (actually, new model number) in the GRIB files. Unless the MDS is ignored,
    ! SILAM is unable to pass through the moment of new MDS appearing first time.
    ! A DANGER: if there were several MDSs from the very beginning with 
    ! overlapping content, it will be absolutely impossible to select the proper one
    ! unless their data can be distibguished somehow else - grid, vertical, etc.
    !
    call set_met_srcs(InputContent, wdr) !set the selected met_srcs to wdr
    
    call fix_met_src(input_shopLst, wdr, .true.) ! set proper MDS in main list and in vars

    !--------------------------------------------------------------
    !
    ! Init supermarket. Pay attention to the timenodes. They define
    ! the number of the meteo periods simultaneously stored in memory.
    !
    ! Initialization has to be done separately for time dependent and permanent parts.
    !
    ! Note: minimum number of time nodes is 3 in order to be able to handle the jumps from 
    ! one forecast to another one. Without these three nodes we will be unable to handle the
    ! cumulative quantities if the jump happens not to the start of the forecast
    !
    allocate(wdrar(1))
    wdrar(1)%ptr => wdr
    CALL initialize_mini_market(meteoMarketPtr, &
                              & 'meteo_market', &
                              & fu_NbrOfMetSrcs(wdr), &  ! nbr of MDS
                              & max(2,nTimeNodesNeeded),&   ! timenodes in memory
                              & 3000,& ! for each timenode - fields
                              & 200,&  ! for each timenode - windfields
                              & 50, &  ! for each timenode - 3d fields
                              & 10, &  ! for each timenode  - 3d windfields
                              & ifForward,& ! if replace oldest (true) or latest (false) when full
                              & wdrar, &
                              &  .true. ,&
                              & .true.) ! info to stdout <=> supermarket_info = .true.
    CALL initialize_mini_market(meteoMarketPtr, &
                              & 'meteo_market', &
                              & 1, &  ! nbr of MDS
                              & 0,&   ! timenodes in memory (assume: max 1 timenode per file)
                              & 200,& ! for each timenode - fields
                              & 0,&  ! for each timenode - windfields
                              & 10, &  ! for each timenode - 3d fields
                              & 0, &  ! for each timenode  - 3d windfields
                              & ifForward,& ! if replace oldest (true) or latest (false) when full
                              & wdrar, &
                              & .false., &
                              & .true.) ! info to stdout <=> supermarket_info = .true.
    call msg('meteo market is initialized')
    IF (error) RETURN

    deallocate(wdrar)
    ! Now set the met_src in the supermarket
    !
    call set_met_srcs_in_sm_single_wdr (meteoMarketPtr, .true. , .false., wdr)

  END SUBROUTINE meteo_init


  !*********************************************************************

  subroutine init_supplementary_market(miniMarketPtr, q_dyn, q_stat, &
                                     & wdr, nTimeNodesNeeded, chMarketName)
    !
    ! Initialises a supplementary mini-market. Essentially, encapsulates
    ! the initialize_mini_market from supermarket_of_fields, with some
    ! minor intellect. However, useful for internal markets.
    !
    implicit none

    ! Imported parameters
    type(mini_market_of_stacks), pointer :: miniMarketPtr
    type(silja_wdr), intent(in), target :: wdr
    integer, dimension(:), intent(in) :: q_dyn, q_stat
    integer, intent(in) :: nTimeNodesNeeded
    character(len=*), intent(in) :: chMarketName

    ! Local variables
    type(wdr_ptr), dimension(:), pointer :: wdrPtr
    integer :: nbr_of_fields, nbr_of_windfields, nbr_of_3d_fields, nbr_of_3d_windfields

    allocate(wdrPtr(1))
    allocate(wdrPtr(1)%ptr)
    wdrPtr(1)%ptr = wdr
    !
    ! Here we have to nullify al lpointers in the minimarket. Otherwise
    ! the status of some stack pointer arrays may be undefined: below initilization
    ! does not garantee that they all will be allocated.
    !
    call nullify_minimarket(miniMarketPtr)
    !
    ! Dynamic dispersion market
    !
    call count_fields(q_dyn, nbr_of_fields, nbr_of_windfields, nbr_of_3d_fields, nbr_of_3d_windfields)
    call msg('Fields to mini_market:' + chMarketName + ':', nbr_of_fields)
    call msg('windfields: ', nbr_of_windfields)
    call msg('3d_fields', nbr_of_3d_fields)
    call msg('3d winds', nbr_of_3d_windfields)

    ! We seem to get various failures if no space is reserved for some
    ! type of field, even if there will be none. That's why the
    ! max(nn, 1) everywhere.
    !
    if (nbr_of_fields > 0) then
      ! 
      ! Create time-dependent dispersion stack. Since the addition of
      ! quantities besides the winds, this needs to be based on the
      ! actual quantity lists.
      !
      call initialize_mini_market(miniMarketPtr, &
                                & chMarketName, &
                                & 1,&      ! first dimension of the stack arrays
                                & nTimeNodesNeeded,& ! if>0, the second dimension in multiTime stack
                                & max(nbr_of_fields, 1), &
                                & max(nbr_of_windfields, 1), &
                                & max(nbr_of_3d_fields, 1), &
                                & max(nbr_of_3d_windfields, 1),&
                                & .true.,& ! replace_earliest_when_full
                                & wdrPtr, & 
                                & .false., & ! ifSingleSrc
                                & test_messages)
      if(error)return
    endif
    
    ! Static stacks
    !
    call count_fields(q_stat, nbr_of_fields, nbr_of_windfields, nbr_of_3d_fields, nbr_of_3d_windfields)
    call msg('Fields to static market:' + chMarketName + ':', nbr_of_fields)
    call msg('windfields: ', nbr_of_windfields)
    call msg('3d_fields', nbr_of_3d_fields)
    call msg('3d winds', nbr_of_3d_windfields)
    
    
    if (nbr_of_fields > 0) then
      !
      ! Initialize the dispersion stack
      !
      call initialize_mini_market(miniMarketPtr, &
                                & chMarketName, &
                                & 1,&      ! first dimension of the stack arrays
                                & 0,& ! if>0, the second dimension in multiTime stack
                                & max(nbr_of_fields+40, 70), & !q_disp_stat has only quantities, which might exist for more than one species (pollen)
                                & max(nbr_of_windfields, 1), &
                                & max(nbr_of_3d_fields, 1), &
                                & max(nbr_of_3d_windfields, 1),& 
                                & .false.,&
                                & wdrPtr, &
                                & .false., &
                                & test_messages)
    endif
    deallocate(wdrPtr)

  contains
    
    subroutine count_fields(q_list, nbr_of_fields, nbr_of_windfields, nbr_of_3d_fields, nbr_of_3d_windfields)
      ! Decompose the quantity list into nbr of fields, windfields,
      ! etc. Tries to be general, but assumes that the winds have 3 components.
      implicit none
      integer, dimension(:), intent(in) :: q_list
      integer, intent(out) :: nbr_of_fields, nbr_of_windfields, nbr_of_3d_fields, nbr_of_3d_windfields

      integer :: ind_q, quantity, nbr_of_2d_wind_components, nbr_of_3d_wind_components

      nbr_of_windfields = 0
      nbr_of_3d_fields = 0
      nbr_of_3d_windfields = 0
      nbr_of_fields = 0
      nbr_of_2d_wind_components = 0
      nbr_of_3d_wind_components = 0

      do ind_q = 1, size(q_list)
        quantity = q_list(ind_q)
        if (quantity == int_missing) exit
        if (fu_multi_level_quantity(quantity)) then
          nbr_of_3d_fields = nbr_of_3d_fields + 1
          ! we'll be in trouble if there's to be more than one 3d wind.
          if (fu_wind_quantity(quantity)) nbr_of_3d_wind_components = nbr_of_3d_wind_components + 1
        else
          nbr_of_fields = nbr_of_fields + 1
          if (fu_wind_quantity(quantity)) nbr_of_2d_wind_components = nbr_of_2d_wind_components + 1
        end if
      end do

      if (mod(nbr_of_3d_wind_components, 3) > 0 .or. mod(nbr_of_2d_wind_components, 3) > 0) then
        call set_error('Strange number of wind components (3d, 2d):' + &
                     & fu_str(nbr_of_3d_wind_components) + ',' + &
                     & fu_str(nbr_of_2d_wind_components), 'count_fields')
        return
      end if

      nbr_of_3d_windfields = nbr_of_3d_wind_components / 3 + 1
      nbr_of_windfields = nbr_of_3d_windfields*nz_dispersion + nbr_of_2d_wind_components / 3
      !
      ! Some dq_ functions make several quantities not asking whether the full set is needed.
      ! So, we have to have a reserve for this case
      !
      nbr_of_3d_fields = nbr_of_3d_fields + 2
      nbr_of_fields = nbr_of_fields + 30 + nz_dispersion * nbr_of_3d_fields

    end subroutine count_fields

  end subroutine init_supplementary_market


  !*************************************************************************

  subroutine init_supplementary_data_buffer(pSupplMarket, data_buffer, &
                                          & ifAllocatePresent2D, grid_size, buffer_q)
    !
    ! Initialises the generic data buffer. Essentially, expands the species
    ! and then calls the generic init_data_buffer. This trick is due to dispersion minimarket
    ! uses fields with non-void species definitions. Thus, there can be several different fields
    ! of the same quantity but of different species
    !
    implicit none
    
    ! Imported parameters
    type(mini_market_of_stacks), pointer :: pSupplMarket
    TYPE(Tfield_buffer),POINTER :: data_buffer
    logical, intent(in) :: ifAllocatePresent2D
    integer, intent(in) :: grid_size
    integer, dimension(:), pointer, optional :: buffer_q  ! the demand: needed quantities
                                                          ! if absent, all available are taken
    ! Local variables
    type(silam_species), dimension(max_variables) :: arSpeciesTmp
    integer, dimension(max_variables) :: quantitiesTmp
    integer :: iQdemand, iQ, number_of_variables
    logical :: ifExclude

    !
    ! Get all what exists in the minimarket
    !
    call supermarket_variables(pSupplMarket, met_src_missing, &
                             & int_missing, &         ! stack_type - take all
                             & quantitiesTmp,&
                             & arSpeciesTmp, &
                             & number_of_variables)
    if(error)return

    if(present(buffer_q))then
      !
      ! Scan the minimarket content selecting only variables whose quantities are in buffer_q
      !
      do iQ = 1, number_of_variables
        ifExclude = .true.
        do iQdemand = 1, size(buffer_q)
          if(quantitiesTmp(iQ) == buffer_q(iQdemand))then  ! Needed?
            ifExclude = .false.
            exit
          endif
        end do
        if(ifExclude)then
          quantitiesTmp(iQ) = quantitiesTmp(number_of_variables)
          arSpeciesTmp(iQ) = arSpeciesTmp(number_of_variables)
          quantitiesTmp(number_of_variables) = int_missing
          number_of_variables = number_of_variables - 1
        endif
      end do
      !
      ! Having consumed the content of the minimarket into quantitiesTmp and asSpeciesTmp and
      ! filtered it with the demand array buffer_q, we still have to allow some other
      ! variables to be introduced later on. The trouble is that we do not know their species,
      ! only quantities. Well, nothing to do now - let' s hope that they will come universal - 
      ! with missing species.
      ! So, we turn around and add all quantities in buffer_q but yet-missing from minimarket.
      !
      do iQdemand = 1, size(buffer_q)
        if(buffer_q(iQdemand) == int_missing)exit  ! all done?
        !
        ! Check the so-far introduced quantities. No duplicates now!
        !
        ifExclude = .false.
        do iQ = 1, number_of_variables
          if(buffer_q(iQdemand) == quantitiesTmp(iQ))then
            ifExclude = .true.
            exit
          endif
        end do
        if(ifExclude)cycle
        number_of_variables = number_of_variables + 1
        arSpeciesTmp(number_of_variables) = species_missing
        quantitiesTmp(number_of_variables) = buffer_q(iQdemand)
      end do
    end if  ! present buffer_q demand array
    !
    ! Now we can call the generic buffer initialization
    !
    call init_data_buffer(data_buffer, quantitiesTmp, ifAllocatePresent2D, grid_size, arSpeciesTmp)
    
  end subroutine init_supplementary_data_buffer


  !*************************************************************************

  subroutine init_data_buffer(dat_buf, buffer_q, ifAllocatePresent2D, grid_size, buffer_sp)
    !
    ! Initialises, allocates memory for the buffer and sets its quantities. 
    ! After that they must not be altered
    !
    implicit none

    ! Imported pointer to the main meteoBuffer, 
    TYPE(Tfield_buffer),POINTER::dat_buf
    integer, dimension(:), intent(in) :: buffer_q
    logical, intent(in) :: ifAllocatePresent2D
    integer, intent(in) :: grid_size
    type(silam_species), dimension(:), intent(in), optional :: buffer_sp

    ! Local variables
    integer :: i, iLev, NbrQ, iStatus

    if(.not.defined(meteo_grid))then
      call set_error('Undefined meteo_grid','init_data_buffer')
      return
    endif

    do NbrQ = 0,size(buffer_q)
      if(buffer_q(NbrQ+1) == int_missing)exit
    end do
    if(NbrQ < 1)then
      call set_error('Undefined buffer_quantities','init_data_buffer')
      return
    endif

    if(associated(dat_buf))then
      call set_error('Buffer is already allocated','init_data_buffer')
      return
    endif

    allocate(dat_buf, stat=iStatus)
    if(iStatus /= 0)then
      call set_error('Failed to allocate the buffer','init_data_buffer')
      return
    endif

    ALLOCATE(dat_buf%p2d(NbrQ), dat_buf%p4d(NbrQ), & ! dat_buf%grd_shift(NbrQ), &
           & dat_buf%buffer_quantities(NbrQ), dat_buf%ifPointerSet(NbrQ), stat=iStatus)
    if(iStatus /= 0)then
      call set_error('Failed to allocate the buffer pointers','init_data_buffer')
      return
    endif
    
    if(present(buffer_sp))then
      allocate(dat_buf%buffer_species(NbrQ), stat=iStatus)
    else
      nullify(dat_buf%buffer_species)
    endif
    if(iStatus /= 0)then
      call set_error('Failed to allocate the buffer species','init_data_buffer')
      return
    endif

    DO i=1,NbrQ
      DO iLev = 1, max_levels
        NULLIFY(dat_buf%p4d(i)%past%p2d(iLev)%ptr)
        NULLIFY(dat_buf%p4d(i)%past%p2d(iLev)%idPtr)
        NULLIFY(dat_buf%p4d(i)%future%p2d(iLev)%ptr)
        NULLIFY(dat_buf%p4d(i)%future%p2d(iLev)%idPtr)
        NULLIFY(dat_buf%p4d(i)%present%p2d(iLev)%ptr)
        NULLIFY(dat_buf%p4d(i)%present%p2d(iLev)%idPtr)
        dat_buf%p4d(i)%past%p2d(iLev)%ifReady = .false.
        dat_buf%p4d(i)%present%p2d(iLev)%ifReady = .false.
        dat_buf%p4d(i)%future%p2d(iLev)%ifReady = .false.
      END DO
      !
      ! 2D pointers can either be allocated or nullified - depending on what 
      ! is going to happen later. They can be used in the same manner as
      ! 4D pointers - to point to the fields in some stack, or can be 
      ! physically made via e.g. time interpolation
      !
      if(ifAllocatePresent2D)then
        !
        ! Data will be physically stored in the p2d_pr
        !
        ALLOCATE(dat_buf%p2d(i)%present%ptr(grid_size),stat=iStatus) ! 2d fields are stored here !!!
        if(iStatus /= 0)then
          call set_error('Failed to allocate the 2D grid','init_data_buffer')
          return
        endif
        allocate(dat_buf%p2d(i)%present%idPtr,stat=iStatus) ! field_id will be made with data
        if(iStatus /= 0)then
          call set_error('Failed to allocate the idPtr','init_data_buffer')
          return
        endif
        call set_missing(dat_buf%p2d(i)%present%idPtr)
      else
        !
        ! p2d_pr will be used as a normal pointer
        !
        NULLIFY(dat_buf%p2d(i)%present%ptr, dat_buf%p2d(i)%present%idPtr)
      endif
      !
      ! p2d_past and p2d_fut are always used as normal pointers
      !
      NULLIFY(dat_buf%p2d(i)%past%ptr, dat_buf%p2d(i)%past%idPtr)
      dat_buf%p2d(i)%past%ifReady = .false.
      dat_buf%p2d(i)%present%ifReady = .false.
      dat_buf%p2d(i)%future%ifReady = .false.
    END DO
    !
    ! Now copy the buffer_q into internal buffer_qnaitities
    !
    do i=1, NbrQ
      dat_buf%buffer_quantities(i) = buffer_q(i)
      if(present(buffer_sp))dat_buf%buffer_species(i) = buffer_sp(i)
      dat_buf%ifPointerSet(i) = .false.
    end do

!    dat_buf%max_grd_shift%x_shift =0.
!    dat_buf%max_grd_shift%y_shift =0.
    dat_buf%weight_past = 1.

    dat_buf%time_past = time_missing
    dat_buf%time_future = time_missing
    dat_buf%time_present = time_missing
    
  end subroutine init_data_buffer


  !*******************************************************************

  subroutine arrange_buffer(miniMarketPtr, met_src, now, model_time_step, data_buffer, vert)
    !
    ! Searches appropriate past and future stacks and calls set_buffer_pointers
    ! to actually set the pointers of the buffer. This two-step activity
    ! has to be done because the buffer_pointers are more general than
    ! meteo_pointers
    !
    ! If called with now == time_missing sets only single-time pointers
    ! replaces set_permanent_buffer_pointers
    IMPLICIT NONE

    ! Imported variables 
    type(meteo_data_source), intent(in) :: met_src
    TYPE(silja_time), INTENT (in) :: now
    type(silja_interval), intent(in) :: model_time_step
    TYPE(Tfield_buffer),POINTER::data_buffer
    type(mini_market_of_stacks), pointer :: miniMarketPtr
    type(silam_vertical), INTENT (in) :: vert
!    type(silja_grid), pointer :: gridPtr

    ! Local variables
    TYPE(silja_stack), POINTER :: stack_past, stack_future
    logical :: ifAllFound
    integer :: i,iTmp, iLev
    REAL :: weight_past
    logical :: ifSTonly  ! Set only pointers from single-time stack

    type(mini_market_of_stacks), pointer :: mmPtr

    ifSTonly = (.not. defined(now))

    ! Reset all pointers
    data_buffer%ifPointerSet(:) = .false.
    data_buffer%weight_past = 1. ! Could be nan_missing???

    mmPtr => miniMarketPtr
   
    if ( .not. ifSTonly) then
      !call msg("Arranging multitime buffer for market:"+ fu_name(miniMarketPtr))
      if (.not. minimarket_initialized(miniMarketPtr, multi_time_stack_flag)) then
        call set_error("Trying to set pointers from uninitialized  multi_time_stack", &
                         & "arrange_buffer")
        return
      endif
      !
      ! Find out the meteostacks for past and future
      ! Looking for a really past and really future stacks - they must not be the same
      ! Under a reasonable assumption that the middle of model time step is NEVER the moment
      ! of meteodata validity, the past and future stacks can be found easily.
      !
      stack_past => fu_closest_sm_met_src_time(miniMarketPtr, met_src, &
                                                & now + model_time_step*0.5, backwards, .true.)
      stack_future => fu_closest_sm_met_src_time(miniMarketPtr, met_src, &
                                                & now + model_time_step*0.5, forwards, .true.)
      if(associated(stack_past, stack_future))then
        call msg_warning('stack_past and stack_future are the same for nnow + model_time_step*0.5=' + &
                       & fu_str(now + model_time_step*0.5), 'arrange_buffer')
      endif
      IF(fu_fails(.not.error,'Problem with meteostacks','arrange_buffer'))RETURN
      !
      ! If we found past and future stacks, set the dynamic pointers
      !
      IF(defined(stack_past) .and. defined(stack_future))  THEN
        data_buffer%time_past = fu_valid_time(stack_past)
        data_buffer%time_future = fu_valid_time(stack_future)
        data_buffer%time_present = now
        !
        ! weight_past should point at the centre of the current time step
        !
        if (data_buffer%time_future == data_buffer%time_past) then
           weight_past = 0.                                ! if step mid-point hits the meteo time
        else
           weight_past = ((data_buffer%time_future - now) - model_time_step * 0.5) / &
                    & (data_buffer%time_future - data_buffer%time_past)
        endif

!        data_buffer%ifPointerSet(:) = .false.
!call msg("Setting pointers for multi_time_stack")
!call msg("Stack_past >"+fu_name(stack_past)+"< _"+fu_str(data_buffer%time_past))
!call msg("Stack_future >"+fu_name(stack_future)+"< _"+fu_str(data_buffer%time_future))

        call set_buffer_pointers(met_src,weight_past, &
                               & data_buffer, &
                               & data_buffer%time_past, data_buffer%time_future, &
                               & now, model_time_step, &
                               & stack_past, stack_future, &
!                               & gridPtr, &
                               & vert, &
                               & .true., &  ! If compute present 2D vars
                               & ifAllFound)
        if(error)return
        !
        ! Now force weight_past back to its true value
        !
        data_buffer%weight_past = weight_past
     else
        ifAllFound = .false.
      endif  ! dynamic stacks have been found
    endif

!    call msg("weight_past",weight_past)
    !
    ! Set from singletime stack
    ! 
    if (.not. minimarket_initialized(miniMarketPtr, single_time_stack_flag)) then
      call set_error("Trying to arrange_buffer from uninitialized  single_time_stack", &
                   & "arrange_buffer")
      return
    endif
    stack_past => fu_stack(miniMarketPtr, 1) ! permanent stack is the first
!call msg("Setting pointers for singletime_time_stack")
!call msg("Stack >"+fu_name(stack_past)+"< _")

    call set_buffer_pointers(met_src_missing, 1.0, &
                           & data_buffer, &
                           & time_missing, time_missing, & !time_past, time_future
                           & time_missing, interval_missing, & !now, model time step (irrelevant here)
                           & stack_past, stack_past, &
!                           & gridPtr, &
                           & vert, & 
                           & .false., &  ! If compute present 2D vars
                           & ifAllFound)
    if(error)return
    !
    ! Now all pointers have to be set (if not only singletime)
    !
    if (.not. (ifAllFound .or. ifSTonly)) then 
      call msg_warning('Not all pointers set','arrange_buffer')
      do i=1,size(data_buffer%buffer_quantities)
        if(.not. data_buffer%ifPointerSet(i))then
          if(associated(data_buffer%buffer_species))then
            call msg('Variable not found:' + &
                   & fu_quantity_string(data_buffer%buffer_quantities(i)) + ',' + &
                   & fu_str(data_buffer%buffer_species(i)))
          else
            call msg('Quantity not found:' + fu_quantity_string(data_buffer%buffer_quantities(i)))
          endif
        endif
      end do
      call set_error('Not all pointers set','arrange_buffer')
    endif

  end subroutine arrange_buffer

  !********************************************************

  
  SUBROUTINE set_buffer_pointers(src, weight_past, &
                               & data_buffer, &
                               & time_past, time_future, now, model_time_step, &
                               & stack_past, stack_future, &
                      !         & gridUsed, &
                               & verticalUsed, &
                               & ifMakePresent2DFields, & 
                               & ifAllFound)
    !
    ! Arranges a 5-d set of pointers to the fields (x-y-z-t-quantity) in 
    ! the data_buffer
    ! Called e.g. each time as new weather data are obtained to
    ! supermarket (not so often indeed)
    ! Order of quantities is given in buffer. If any quantity exists in 
    ! that list but not found in supermarket - the ifAllDone is set to false
    ! It also checks correspondence of the grids to the system_grid
    ! and sets an error in case of large discrepancies.
    !
    ! Author: Mikhail Sofiev, FMI
    !
    IMPLICIT NONE

    ! Imported variables 
    type(meteo_data_source), INTENT(in)  :: src
    TYPE(silja_time), INTENT (in) :: time_past, time_future, now
    type(silja_interval), intent(in) :: model_time_step
    REAL, INTENT(in) :: weight_past
    TYPE(Tfield_buffer),POINTER :: data_buffer
    TYPE(silja_stack), POINTER :: stack_past, stack_future
    type(silam_vertical), INTENT(in) :: verticalUsed
!    type(silja_grid), pointer :: gridUsed
    logical, intent(in) :: ifMakePresent2DFields
    logical, intent(out) :: ifAllFound

    ! Local variables
    INTEGER :: iQ, iLev, nQ2d, nQ3d, i
    TYPE(silja_field), POINTER :: field_ptr
    TYPE(silja_3d_field), POINTER :: field_3d_ptr
    LOGICAL ::  ifMultiLevel, ifOK
!    REAL, DIMENSION(2) :: grd_shift ! x_shift,y_shift regarding the system_grid
!    real :: shift_lon, shift_lat
    integer, dimension(max_quantities) :: q2d, q3d
    type(silam_species), dimension(max_quantities) :: Species
    TYPE(Tfield_buffer),POINTER::db
    TYPE(silja_stack), POINTER :: st
    type(silam_species) :: spTmp
    TYPE(silja_field_id), pointer :: idPtr


    db => data_buffer
    st => stack_past
    !
    ! Now, get the complete list of stack 2D variables (variable=quantity+species)
    ! and 3D quantities
    !
    CALL stack_variables(stack_past, q2d, Species, nQ2d) ! In fact, both 2d and 3d quantities
    CALL stack_3d_quantities(stack_past, q3d, nQ3d)
    if(error)return
!    !
!    ! Remove the 3d quantity from 2d list
!    !
!    do iQ = 1, nQ2d
!      if(any(q3d(1:nQ3d) == q2d(iQ))) q2d(iQ)=int_missing
!    end do
!    call compress_int_array(q2d, int_missing, nQ2d)  ! remove missing values

    !-------------------------------------------------
    !
    ! 2. Scan all buffer_quantities and fill-in the pointer
    !    arrays. field_grid >< system_grid is to be checked
    !    Separation of 2d / 3d quantities is also made here
    !
!call msg("set_buffer_pointers: Processing stack:"+fu_name(st))
    data_buffer%nbr_of_levels = fu_NbrOfLevels(verticalUsed)
    DO iQ = 1, SIZE(data_buffer%buffer_quantities)
!      if ( data_buffer%ifPointerSet(iQ) ) cycle ! already set

 
!call msg("Scanning:" + fu_quantity_string(data_buffer%buffer_quantities(iQ)), iQ)

      IF(data_buffer%buffer_quantities(iQ) == int_missing) exit ! All done ?
      if(data_buffer%ifPointerSet(iQ)) cycle ! Do not reset the ready-made pointer
      if(associated(data_buffer%buffer_species))then
        spTmp = data_buffer%buffer_species(iQ)
      else
        spTmp = species_missing
      endif
      !
      ! If vertical-dependent z-componenet of wind, we have to look for a correct
      ! variable
      !
      if(data_buffer%buffer_quantities(iQ) == vertical_velocity_flag) &
                            & data_buffer%buffer_quantities(iQ) = vertical_velocity_pointer

!      print *, iQ, data_buffer%buffer_quantities(iQ),' ',fu_quantity_string(data_buffer%buffer_quantities(iQ))

      !
      ! Check 3D and, if not there, 2D. Note the use of fu_index for simultaneous quantity and species
      ! co-location. See chemical_setup module for the function implementation
      !
      if(fu_quantity_in_quantities(data_buffer%buffer_quantities(iQ),q3d))then
        ifMultiLevel = .true.
      else
        ifMultiLevel = .false.
        
        if (spTmp == species_missing) then
           i = fu_index(data_buffer%buffer_quantities(iQ),Q2d,nQ2d)
        else
           i = fu_index_of_variable(data_buffer%buffer_quantities(iQ), q2d, spTmp, Species, nQ2d)
        endif        
        if(i == int_missing)then
          data_buffer%ifPointerSet(iQ) = .false.
          cycle
        endif
      endif

      !-------------------------------------------------------
      !
      ! Now start actual processing
      !
      IF(ifMultiLevel)THEN
        !
        !  Past 3D pointers
        !
        CALL find_field_3d_from_stack(src, &
                                    & data_buffer%buffer_quantities(iQ),&
                                    & time_past,&
                                    & stack_past,&
                                    & field_3d_ptr,&
                                    & ifOK, &
                                    & spTmp)
        if(error.or..not.ifOK)then 
          call msg("Failed top set past 3d pointer for: " &
                  & +fu_quantity_string(data_buffer%buffer_quantities(iQ))+&
                  & "from stack:"+ fu_name(stack_past))
          data_buffer%ifPointerSet(iQ) = .false.
          if(error) call unset_error('set_buffer_pointers')
          cycle
        endif


        if(fu_number_of_fields(field_3d_ptr) < data_buffer%nbr_of_levels)then
          call msg('')
          call report(field_3d_ptr)
          call set_error('Strange nr of levels in the 3D field','set_buffer_pointers')
          return
        endif
        if(.not. fu_level_belongs_to_vertical(fu_level(fu_field_from_3d_field(field_3d_ptr,1)), &
                                            & verticalUsed))then
          call msg('')
          call msg('Vertical given:')
          call report(verticalUsed)
          call msg('3D field to handle:')
          call report(field_3d_ptr)
          call set_error('Can not handle different level types','set_buffer_pointers')
          return
        endif

        data_buffer%p4d(iQ)%past%field3d => field_3d_ptr
        DO iLev = 1, data_buffer%nbr_of_levels
          data_buffer%p4d(iQ)%past%p2d(iLev)%ptr => fu_grid_data_from_3d(field_3d_ptr,iLev)
          data_buffer%p4d(iQ)%past%p2d(iLev)%idPtr => &
                                 & fu_id(fu_field_from_3d_field(field_3d_ptr,iLev))
          data_buffer%p4d(iQ)%past%p2d(iLev)%ifReady = .not.error
        END DO


        !-----------------------------------------------------
        !
        ! Future 3D pointers. No grid_shift set - it is the same
        ! If last model time step - nothing to search, just set the same
        ! And here the field must be found because otherwise it means that
        ! stacks are not similar.
        !
        CALL find_field_3d_from_stack(src, &
                                    & data_buffer%buffer_quantities(iQ),&
                                    & time_future,&
                                    & stack_future,&
                                    & field_3d_ptr,&
                                    & ifOK, &
                                    & spTmp)
        if(.not.ifOK)then  
          call set_error('No 3d future field in supermarket:' + &
                       & fu_quantity_string(data_buffer%buffer_quantities(iQ)) + ',' + &
                       & fu_str(spTmp), &
                       & 'set_buffer_pointers')
          return
        endif

        data_buffer%p4d(iQ)%future%field3d => field_3d_ptr
        DO iLev = 1, data_buffer%nbr_of_levels
          data_buffer%p4d(iQ)%future%p2d(iLev)%ptr => fu_grid_data_from_3d(field_3d_ptr,iLev)
          data_buffer%p4d(iQ)%future%p2d(iLev)%idPtr =>fu_id(fu_field_from_3d_field(field_3d_ptr,iLev))
          data_buffer%p4d(iQ)%future%p2d(iLev)%ifReady = .not.error
          idPtr => data_buffer%p4d(iQ)%future%p2d(iLev)%idPtr
!                    call msg("")
!          call msg("Lev",iLev)
!          call report(data_buffer%p4d(iQ)%future%p2d(iLev)%idPtr)
        END DO
        !
        !
        if (time_past == time_missing .and. time_future == time_missing) then
          !Single-time quantity: can set present pointer
!call msg('Single-time 4D quantity set:' + fu_quantity_string(data_buffer%buffer_quantities(iQ)))
          do iLev = 1, data_buffer%nbr_of_levels
            data_buffer%p4d(iQ)%present%p2d(iLev)%idPtr   => data_buffer%p4d(iQ)%past%p2d(iLev)%idPtr
            data_buffer%p4d(iQ)%present%p2d(iLev)%ptr     => data_buffer%p4d(iQ)%past%p2d(iLev)%ptr
            data_buffer%p4d(iQ)%present%p2d(iLev)%ifReady =  data_buffer%p4d(iQ)%past%p2d(iLev)%ifReady
            data_buffer%p4d(iQ)%present%field3d => data_buffer%p4d(iQ)%past%field3d
          enddo
        else
          ! Time step is made, so all "present" fields are no longer valid.
          data_buffer%p4d(iQ)%present%p2d(:)%ifReady = .false.
        endif

      ELSE  ! multi_level_quantity

!        print *, 'Enter 2D'
        !------------------------------------
        !
        !   2d quantity - find past pointers
        !
        CALL find_field_from_stack(src, &
                                 & data_buffer%buffer_quantities(iQ),&
                                 & time_past,&
                                 & stack_past,&
                                 & field_ptr,&
                                 & ifOK, &
                                 & spTmp)
        if(.not.ifOK)then
          data_buffer%ifPointerSet(iQ) = .false.
!          call msg("Not found past"+fu_quantity_string(data_buffer%buffer_quantities(iQ)))
          cycle
        endif
!         call msg("Found past"+fu_quantity_string(data_buffer%buffer_quantities(iQ)))

!        CALL grid_chk(fu_grid(field_ptr), gridUsed, grd_shift, ifOK)
!        call grid_shift_indices(fu_grid(field_ptr), gridUsed, shift_lon, shift_lat)
!        IF((shift_lon .eps. real_missing) .or. (shift_lat .eps. real_missing))THEN
!!        IF(.not.ifOK)THEN
!          call msg_warning('Failed grid check')
!          call msg('ID for checking:')
!          call report(fu_id(field_ptr))
!          call msg('reference grid :')
!          call report(gridUsed)
!          CALL set_error('Failed grid check for:' + &
!                       & fu_quantity_string(data_buffer%buffer_quantities(iQ)) + ',' + &
!                       & fu_str(spTmp), &
!                       & 'set_buffer_pointers')
!          RETURN
!        END IF

        !
        ! Store the past pointers
        !
        data_buffer%p2d(iQ)%past%ptr => fu_grid_data(field_ptr)
        data_buffer%p2d(iQ)%past%idPtr => fu_id(field_ptr)
        data_buffer%p2d(iQ)%past%ifReady = .true.

        !--------------------------------------------------------------
        !
        ! 2D quantity. Find and store the future pointers
        !
        CALL find_field_from_stack(src, &
                                 & data_buffer%buffer_quantities(iQ),&
                                 & time_future,&
                                 & stack_future,&
                                 & field_ptr,&
                                 & ifOK, &
                                 & spTmp)
        if(.not.ifOK)then  
          data_buffer%ifPointerSet(iQ) = .false.
!          call msg("Not found future"+fu_quantity_string(data_buffer%buffer_quantities(iQ)))
          cycle
        endif
!         call msg("Found future"+fu_quantity_string(data_buffer%buffer_quantities(iQ)))
        data_buffer%p2d(iQ)%future%ptr => fu_grid_data(field_ptr)
        data_buffer%p2d(iQ)%future%idPtr => fu_id(field_ptr)
        data_buffer%p2d(iQ)%future%ifReady = .true.

        !--------------------------------------
        !
        ! If interpolation is needed - do it. Do not overlook:
        ! p2d_pr may be not allocated yet. Check first
        !
        if (time_past == time_missing .and. time_future == time_missing) then
            !Realtime field, can make present at no cost....
      !      call msg('Single-time 2D quanity set:' + fu_quantity_string(data_buffer%buffer_quantities(iQ)))
            data_buffer%p2d(iQ)%present%idPtr  => data_buffer%p2d(iQ)%past%idPtr
            data_buffer%p2d(iQ)%present%ptr    => data_buffer%p2d(iQ)%past%ptr
            data_buffer%p2d(iQ)%present%ifReady = data_buffer%p2d(iQ)%past%ifReady

        elseif(ifMakePresent2DFields)then
          !
          ! Checking the memory allocation
          !
          if(.not. associated(data_buffer%p2d(iQ)%present%ptr))then
            allocate(data_buffer%p2d(iQ)%present%idPtr, &
                   & data_buffer%p2d(iQ)%present%ptr(fu_number_of_gridpoints(fu_grid( &
                                         & data_buffer%p2d(iQ)%past%idPtr))), stat=i)
            if(i /= 0)then
              call set_error('Failed to allocate memory','set_buffer_pointers')
              return
            endif
          endif

          IF(weight_past .eps. 1.0)THEN
            data_buffer%p2d(iQ)%present%ptr = data_buffer%p2d(iQ)%past%ptr
            data_buffer%p2d(iQ)%present%idPtr = data_buffer%p2d(iQ)%past%idPtr
          ELSE

            call interpolate_fields_in_time(data_buffer%p2d(iQ)%past%idPtr, & 
                                          & data_buffer%p2d(iQ)%past%ptr, &
                                          & data_buffer%p2d(iQ)%future%idPtr, &
                                          & data_buffer%p2d(iQ)%future%ptr, &
                                          & weight_past, &
                                          & data_buffer%p2d(iQ)%present%idPtr, &
                                          & data_buffer%p2d(iQ)%present%ptr, &
                                          & now)
            if(error)then
              call msg_test('NOW is: ')
              call report(now)
              call msg_test('Past id is:')
              call report(data_buffer%p2d(iQ)%past%idPtr)
              call msg_test('Future id is:')
              call report(data_buffer%p2d(iQ)%future%idPtr)
              return
            endif

          END IF  ! Weight_past == 0

          ! Set validity during the current timestep: this is when the values are to be used
          ! Cumulative viariables are of special kind, they are handled by those who use them.
          !
          if (.not. fu_accumulated(data_buffer%p2d(iQ)%present%idPtr))then
            if(fu_interval_positive(model_time_step))then
              call set_valid_time(data_buffer%p2d(iQ)%present%idPtr, now)
              call set_validity_length(data_buffer%p2d(iQ)%present%idPtr, model_time_step)
            else
              call set_valid_time(data_buffer%p2d(iQ)%present%idPtr, now + model_time_step)
              call set_validity_length(data_buffer%p2d(iQ)%present%idPtr, fu_abs(model_time_step))
            endif
          endif
          data_buffer%p2d(iQ)%present%ifReady = .true.
        endif ! Perform time interpolation of 2D quantities

      END IF ! multi-level quantity

      data_buffer%ifPointerSet(iQ) = .true.
!call msg('Pointer set')
      
    END DO ! Quantities in data_buffer

!    print *, 'Done'

    ! Do not do it here!!!! 
    !data_buffer%weight_past = weight_past

    ifAllFound = all(data_buffer%ifPointerSet(1:iQ-1))

  END SUBROUTINE set_buffer_pointers


  ! ****************************************************************


  INTEGER FUNCTION fu_index_for_quantity(q_ptr, quantity)
    !
    ! Finds out the index in the given array q_ptr for the given quantity
    ! 
    IMPLICIT NONE

    ! Imported parameters with the intent IN
    INTEGER, INTENT(in) :: quantity

    ! Imported POINTERS
    INTEGER, DIMENSION(:), intent(in):: q_ptr

    !Local variables
    INTEGER :: i

    fu_index_for_quantity = fu_index(quantity, q_ptr)

   !    DO i=1, SIZE(q_ptr)
   !   IF(q_ptr(i) == int_missing)RETURN
   !   IF(q_ptr(i) == quantity)THEN
   !     fu_index_for_quantity = i
   !     RETURN
   !   END IF
   ! END DO
  END FUNCTION fu_index_for_quantity


  ! ****************************************************************


  INTEGER FUNCTION fu_index_for_quantity_in_buf(buffer, quantity, ifStrict)
    !
    ! Finds out the index in the given array q_ptr for the given quantity
    ! 
    IMPLICIT NONE

    ! Imported parameters 
    INTEGER, INTENT(in) :: quantity
    type(Tfield_buffer), intent(in) :: buffer
    logical, intent(in), optional :: ifStrict

    !Local variables
    INTEGER :: i

    fu_index_for_quantity_in_buf = int_missing

    !call msg('Quantity to find:' + fu_quantity_short_string(quantity), quantity)

    DO i=1, SIZE(buffer%buffer_quantities)
      IF(buffer%buffer_quantities(i) == int_missing)RETURN
      IF(buffer%buffer_quantities(i) == quantity)THEN
        fu_index_for_quantity_in_buf = i
        RETURN
      END IF
    END DO
    ! 
    ! Do we really want this variable?
    !
    if(fu_index_for_quantity_in_buf == int_missing)then
      if(present(ifStrict))then
        if(ifStrict)then
          call msg_warning('Variable >>' + fu_quantity_string(quantity) + &
                         & '<< is really wanted but absent in the buffer','fu_index_for_quantity_in_buf')
          call msg('Buffer content:')
          call report_buffer(buffer,1)  ! just report quantities
          call set_error('Variable >>' + fu_quantity_string(quantity) + &
                       & '<< is really wanted but absent in the buffer','fu_index_for_quantity_in_buf')
        endif
      endif
    endif
  END FUNCTION fu_index_for_quantity_in_buf


  !******************************************************************
  
  INTEGER FUNCTION fu_index_for_variable_in_buf(buffer, quantity, species, ifStrict)
    !
    ! Finds out the index in the given array q_ptr for the given quantity and species
    ! 
    IMPLICIT NONE

    ! Imported parameters 
    INTEGER, INTENT(in) :: quantity
    type(Tfield_buffer), POINTER :: buffer
    type(silam_species), intent(in) :: species
    logical, intent(in) :: ifStrict

    !Local variables
    INTEGER :: iTmp
    type(Tfield_buffer), POINTER :: bf

    fu_index_for_variable_in_buf = int_missing

    bf => buffer
!call msg('Quantity to find:' + fu_quantity_short_string(quantity), quantity)

    DO iTmp = 1, SIZE(buffer%buffer_quantities)
      IF(buffer%buffer_quantities(iTmp) == int_missing)RETURN

!call msg('Quantity available:' + fu_quantity_short_string(bf%buffer_quantities(i)), bf%buffer_quantities(i))

      IF(buffer%buffer_quantities(iTmp) == quantity)THEN
        if(defined(species))then
          if(associated(buffer%p4d(iTmp)%past%p2d(1)%ptr))then
            if(.not. (species == fu_species(buffer%p4d(iTmp)%past%p2d(1)%idPtr)))cycle
          else
            if(.not. (species == fu_species(buffer%p2d(iTmp)%past%idPtr)))cycle
          endif
        endif
        fu_index_for_variable_in_buf = iTmp
        RETURN
      END IF
    END DO
    ! 
    ! Do we really want this variable?
    !
    if(fu_index_for_variable_in_buf == int_missing)then
      if(ifStrict)then
        call msg_warning('Variable >>' + fu_quantity_string(quantity) + &
                       & '<< is really wanted but absent in the buffer','fu_index_for_variable_in_buf')
        call msg('Requested species:')
        call report(species)
        call msg('Buffer content:')
        call report_buffer(buffer,1)  ! just report quantities
        call set_error('Variable >>' + fu_quantity_string(quantity) + &
                     & '<< is really wanted but absent in the buffer','fu_index_for_variable_in_buf')
      endif
    endif
    
  END FUNCTION fu_index_for_variable_in_buf


  ! ****************************************************************


  INTEGER FUNCTION fu_idx4q_in_buf_with_ptr2d(buffer, quantity, ptrVal)
    !
    ! Finds out the index in the given array q_ptr for the given quantity
    ! 
    IMPLICIT NONE

    ! Imported parameters 
    INTEGER, INTENT(in) :: quantity
    type(Tfield_buffer), intent(in) :: buffer
    real, dimension(:), pointer :: ptrVal

    !Local variables
    INTEGER :: i

    fu_idx4q_in_buf_with_ptr2d = int_missing

    !call msg('Quantity to find:' + fu_quantity_short_string(quantity), quantity)

    DO i=1, SIZE(buffer%buffer_quantities)
      IF(buffer%buffer_quantities(i) == int_missing)RETURN
      IF(buffer%buffer_quantities(i) == quantity)THEN
        if(associated(buffer%p4d(i)%past%p2d(1)%ptr))then
          call set_error('4d quantity instead of 2d'+ fu_quantity_string(quantity),'fu_idx4q_in_buf_with_ptr2d')
        else
          fu_idx4q_in_buf_with_ptr2d = i
          ptrVal => buffer%p2d(i)%present%ptr
        endif
        RETURN
      END IF
    END DO
    call set_error('Failed:' + fu_quantity_string(quantity),'fu_idx4q_in_buf_with_ptr2d')
  END FUNCTION fu_idx4q_in_buf_with_ptr2d


  !******************************************************************
  
  INTEGER FUNCTION fu_idx4var_in_buf_with_ptr2d(buffer, quantity, species, ptrVal)
    !
    ! Finds out the index in the given array q_ptr for the given quantity and species
    ! 
    IMPLICIT NONE

    ! Imported parameters 
    INTEGER, INTENT(in) :: quantity
    type(Tfield_buffer), POINTER :: buffer
    type(silam_species), intent(in) :: species
    real, dimension(:), pointer :: ptrVal

    !Local variables
    INTEGER :: iTmp
    type(Tfield_buffer), POINTER :: bf

    fu_idx4var_in_buf_with_ptr2d = int_missing

    bf => buffer
!call msg('Quantity to find:' + fu_quantity_short_string(quantity), quantity)

    DO iTmp = 1, SIZE(buffer%buffer_quantities)
      IF(buffer%buffer_quantities(iTmp) == int_missing)RETURN

!call msg('Quantity available:' + fu_quantity_short_string(bf%buffer_quantities(i)), bf%buffer_quantities(i))

      IF(buffer%buffer_quantities(iTmp) == quantity)THEN
        if(defined(species))then
          if(associated(buffer%p4d(iTmp)%past%p2d(1)%ptr))then
            if(.not. (species == fu_species(buffer%p4d(iTmp)%past%p2d(1)%idPtr)))cycle
          else
            if(.not. (species == fu_species(buffer%p2d(iTmp)%past%idPtr)))cycle
          endif
        endif
        if(associated(buffer%p4d(iTmp)%past%p2d(1)%ptr))then
          call set_error('4d quantity instead of 2d'+ fu_quantity_string(quantity),'fu_idx4var_in_buf_with_ptr2d')
        else
          fu_idx4var_in_buf_with_ptr2d = iTmp
          ptrVal => buffer%p2d(iTmp)%present%ptr
        endif
        RETURN
      END IF
    END DO
  END FUNCTION fu_idx4var_in_buf_with_ptr2d


  !******************************************************************

  integer function fu_dimension_of_meteofield(met_buf, indexFld)
    !
    ! There can be 2d field or 4d field. Reason is that time interpolation
    ! is made only for 2d-spatial fields, not for spatial 3d fields
    ! There is another problem - particular meteovariable can appear at
    ! just one level, being in principle multi-level one. In this case
    ! it will still be referred here as a single-level variable.
    !
    implicit none

    ! Imported parameters
    integer, intent(in) :: indexFld
    TYPE(Tfield_buffer), intent(in)::met_buf

    if(met_buf%p4d(indexFld)%past%p2d(1)%ifReady)then
      fu_dimension_of_meteofield = 4
    else
      fu_dimension_of_meteofield = 2
    endif

  end function fu_dimension_of_meteofield


  ! ****************************************************************


  function fu_2d_field_from_buffer(met_buf, indexFld)result(fieldPtr)
    !
    ! Returns the pointer to 2D field.
    !
    implicit none

    ! Return value - directly array
    type(field_2d_data_ptr), pointer :: fieldPtr

    ! Imported parameters
    TYPE(Tfield_buffer),POINTER::met_buf
    integer, intent(in) :: indexFld

    if(indexFld < 1 .or. indexFld > size(met_buf%p2d))then
      call set_error('Strange idnex','fu_2d_field_from_meteobuffer')
      call msg('Index: ', indexFld)
      nullify(fieldPtr)
      return
    else
      fieldPtr => met_buf%p2d(indexFld)
    endif

  end function fu_2d_field_from_buffer



  ! ****************************************************************


  function fu_4d_field_from_buffer(met_buf, indexFld)result(fieldPtr)
    !
    ! Returns the pointer to 2D or 3D field
    !
    implicit none

    ! Return value - directly array
    type(field_4d_data_ptr), pointer :: fieldPtr

    ! Imported parameters
    integer, intent(in) :: indexFld
    TYPE(Tfield_buffer),POINTER::met_buf

    if(indexFld < 1 .or. indexFld > size(met_buf%p4d))then
      call msg('Index: ', indexFld)
      call set_error('Strange idnex','fu_4d_field_from_meteobuffer')
      nullify(fieldPtr)
      return
    endif
    fieldPtr => met_buf%p4d(indexFld)

  end function fu_4d_field_from_buffer


  ! ****************************************************************


  real function fu_get_value_from_buffer(buf, quantity_index, ix,iy,iz, nx, now)
    !
    ! Returns the value from the given field buffer for the given coordinates
    ! Since the call is expected from within the internal cycle, checking
    ! is minimised
    !
    implicit none

    ! Imported parameters
    type(Tfield_buffer), pointer :: buf
    integer, intent(in) :: ix,iy,iz, nx, quantity_index
    type(silja_time), optional :: now

    ! Local variables
    type(Tfield_buffer), pointer :: b

    !
    ! Find the quantity, get its dimension and then the value.
    ! Beware about the 4d quantities - they are not necessarily ready
    ! for the present time.
    !
#ifdef DEBUG    
    !$ if (omp_in_parallel()) then
    !$ call set_error("Not thread-safe!","fu_get_value_from_buffer")
    !$ 
    !$ endif
#endif    

!    print *, quantity_index, iz, ix, iy

    b => buf

    if(buf%p4d(quantity_index)%past%p2d(1)%ifReady)then
      !
      ! 4d - there is a vertical component.
      !
      if(.not. buf%p4d(quantity_index)%present%p2d(iz)%ifReady)then
        !
        ! Present value is not ready, make it from past and future
        !
        if(present(now))then
          call make_present_variable(buf, quantity_index, iz, &
                                 & size(buf%p4d(quantity_index)%past%p2d(1)%ptr), now)
        else
          call set_error('Now-time absent but the present-time field is not ready', &
                       & 'fu_get_value_from_buffer')
          fu_get_value_from_buffer = real_missing
          return
        endif
      endif

      fu_get_value_from_buffer = buf%p4d(quantity_index)%present%p2d(iz)%ptr(ix+(iy-1)*nx)

    else
      !
      ! 2d - no vertical component, present must be already available
      !
      fu_get_value_from_buffer = buf%p2d(quantity_index)%present%ptr(ix+(iy-1)*nx)

    endif

  end function fu_get_value_from_buffer

  ! ****************************************************************

  ! 
  !
  !          4-d INTERPOLATION of fields
  !
  !
  ! ***************************************************************


  REAL FUNCTION fu_4d_interp_to_position(& !met_buf, &
                                       & data_4d, &
                                       & position, &
                                       & pr_index, &
                                       & past_weight, &
                                       & method_hor, method_vert, &
                                       & grid) &
                & result(value)
    ! 
    ! Interpolates the 4-d field given in a grid cube to a given position.
    !
    ! IMPORTANT. The routine is made as fast as reasonable - it does not 
    ! operate with derived types, but rather takes 4-d cube of data arrays
    ! (x-y-z-t) with size (nx-ny-nLayers-2) and finds the value corresponding
    ! to the position. Since time weight coefficient must be already known -
    ! it is just taken and not computed anyhow. The same is true for vertical
    ! co-ordinate - its relative value is an input parameter.
    !
    ! IMPORTANT. There is faster routine interpolating to co-ordinates. Use it!
    !
    IMPLICIT NONE
    !
    ! Imported parametersn with intent(in):
    TYPE(silam_grid_position), INTENT(in) :: position
    INTEGER, INTENT(in) :: method_hor, method_vert ! of interpolation
    REAL, INTENT(in) :: pr_index ! Relative pressure index 
    REAL, INTENT(in) :: past_weight ! interpolation weight coefficient
    TYPE(field_4d_data_ptr), INTENT(in)::data_4d
!    TYPE(Tfield_buffer),POINTER::met_buf
    TYPE(silja_grid), INTENT(in), OPTIONAL, TARGET :: grid

    ! Local declarations:
    REAL val_future
    TYPE(silja_grid), POINTER :: grid_ptr
    INTEGER :: n_levs

    IF(PRESENT(grid))THEN
      grid_ptr => grid
    ELSE
      grid_ptr => meteo_grid
    END IF

    SELECT CASE(method_vert) ! So far only linear 
    CASE(linear)

      n_levs = SIZE(data_4d%past%p2d)

      IF(REAL(INT(pr_index+0.5)).eps.pr_index)THEN
        !
        ! Take an exact level
        !
        value = fu_2d_interpolation(&
                      & data_4d%past%p2d(max(1,MIN(INT(pr_index+0.5),n_levs)))%ptr, &
                      & grid_ptr, &
                      & position, &
                      & method_hor)

        IF(past_weight.eps.1.0)RETURN !Speed-up if now == time_past

        val_future = fu_2d_interpolation( &
                      & data_4d%future%p2d(max(1,MIN(INT(pr_index+0.5),n_levs)))%ptr, &
                      & grid_ptr, &
                      & position, &
                      & method_hor)
        ! Time interpolation:
        value = past_weight * value + (1.-past_weight) * val_future
      ELSE

        value = fu_2d_interpolation(&
                      & data_4d%past%p2d(max(1,MIN(INT(pr_index),n_levs)))%ptr, &
                      & grid_ptr, &
                      & position, &
                      & method_hor) * (1.-pr_index+INT(pr_index)) + &
                 & fu_2d_interpolation(&
                      & data_4d%past%p2d(MIN(INT(pr_index)+1,n_levs))%ptr, &
                      & grid_ptr, &
                      & position, &
                      & method_hor) * (pr_index-INT(pr_index))

        IF(past_weight.eps.1.0)RETURN !Speed-up if now == time_past

        val_future = fu_2d_interpolation( &
                      & data_4d%future%p2d(max(1,MIN(INT(pr_index),n_levs)))%ptr, &
                      & grid_ptr, &
                      & position, &
                      & method_hor) * (1.-pr_index+INT(pr_index)) + &
                 & fu_2d_interpolation( &
                      & data_4d%future%p2d(MIN(INT(pr_index)+1,n_levs))%ptr, &
                      & grid_ptr, &
                      & position, &
                      & method_hor) * (pr_index-INT(pr_index))

        ! Time interpolation:
        value = past_weight * value + (1.-past_weight) * val_future
      END IF ! if pr_index is integer
      RETURN

    CASE DEFAULT
      CALL set_error('unknown vertical interpolation method','4d_interpolation')
      RETURN
    END SELECT

  END FUNCTION fu_4d_interp_to_position


  !**********************************************************************************

  REAL FUNCTION fu_4d_interp_to_coord(data_4d, &
                                    & x,y,z,nx,ny,nz, &  ! Co-ordinates and grid size
                                    & past_weight, & ! Time index
                                    & method_hor, method_vert, iOut) &
                & result(val)
    ! 
    ! Interpolates the 4-d field given in a grid cube to a given grid co-ordinates.
    ! IMPORTANT. The routine is made as fast as reasonable - it does not 
    ! operate with derived types, but rather takes 4-d cube of data arrays
    ! (x-y-z-t) with size (nx-ny-nLayers-2) and finds the value corresponding
    ! to the co-ordinate. Since time weight coefficient must be already known -
    ! it is just taken and not computed anyhow. The same is true for vertical
    ! co-ordinate - its relative value is an input parameter.
    !
    ! IMPORTANT. Use this routine rather than interpolation to position above - 
    !            it is faster
    ! 
    IMPLICIT NONE
    !
    ! Imported parametersn with intent(in):
    INTEGER, INTENT(in) :: nx,ny,nz,method_hor, method_vert ! of interpolation
    REAL, INTENT(in) :: x,y,z   ! Relative pressure index 
    REAL, INTENT(in) :: past_weight ! interpolation weight coefficient
    TYPE(field_4d_data_ptr), INTENT(in)::data_4d
    integer, intent(in) :: iOut
    
    ! Local declarations:
    REAL val_future

    SELECT CASE(method_vert) ! So far only linear 
    CASE(linear)

      IF(REAL(INT(z+0.5)) .eps. z)THEN
        val = fu_2d_interpolation(data_4d%past%p2d(max(1,MIN(INT(z+0.5),nz)))%ptr, &
                                & x,y, nx,ny, method_hor, iOut)

        IF(past_weight.eps.1.0)RETURN !Speed-up if now == time_past

        val_future = fu_2d_interpolation(data_4d%future%p2d(max(1,MIN(INT(z),nz)))%ptr, &
                                       & x,y, nx,ny, method_hor, iOut)
        ! Time interpolation:
        val = past_weight * val + (1.-past_weight) * val_future
      ELSE

        val = fu_2d_interpolation(data_4d%past%p2d(max(1,MIN(INT(z),nz)))%ptr, &
                                & x,y, nx,ny, method_hor, iOut) * (1.-z+INT(z)) + &
            & fu_2d_interpolation(data_4d%past%p2d(MIN(INT(z)+1,nz))%ptr, &
                                & x,y, nx,ny, method_hor, iOut) * (z-INT(z))

        IF(past_weight.eps.1.0)RETURN !Speed-up if now == time_past

        val_future = fu_2d_interpolation(data_4d%future%p2d(max(1,MIN(INT(z),nz)))%ptr, &
                                       & x,y, nx,ny, method_hor, iOut) * (1.-z+INT(z)) + &
                   & fu_2d_interpolation(data_4d%future%p2d(MIN(INT(z)+1,nz))%ptr, &
                                       & x,y, nx,ny, method_hor, iOut) * (z-INT(z))

        ! Time interpolation:
        val = past_weight * val + (1.-past_weight) * val_future
      END IF ! if pr_index is integer

    CASE DEFAULT
      CALL set_error('unknown vertical interpolation method','4d_interpolation')
      RETURN
    END SELECT

  END FUNCTION fu_4d_interp_to_coord


  !****************************************************************

  subroutine make_2d_from_4d_field (fld_buf, &   ! Field buffer
                                  & varIndex, &  ! Index of variable in buffer to interpolate
                                  & verticalUsed, & ! Vertical structure of the input buffer
                                  & levelTarget, &  ! Level to be interpolated to
                                  & iVerticalTreatment, &  ! may be, integrate/average along the vertical
                                  & now, &          ! target time to interpolate to
                                  & dataPtr, &      ! Output data array
                                  & idRes)          ! Output field id
    !
    ! Interpolates in time and vertical directions, thus making a 2D field from the 4d cube.
    ! Important:
    ! The interpolation is made so that time step is done first for the 
    ! needed levels. The result is stored into the "present" pointer.
    ! However, only needed levels are interpolated. Other levels are
    ! not touched, thus saving time
    !
    implicit none

    ! Imported parameters
    TYPE(Tfield_buffer), intent(inout) :: fld_buf
    integer, intent(in) :: varIndex, iVerticalTreatment  ! variable to be interpolated
    type(silam_vertical), intent(in) :: verticalUsed  ! Vertical in which fields are stored in buffer
    type(silja_level), intent(in) :: levelTarget   ! level which interpoate to
    type(silja_time), intent(in) :: now            ! Time which interpolate to 
    real, dimension(:), intent(out) :: dataPtr         ! Resulting data array
    type(silja_field_id), intent(inout), optional :: idRes ! Resulting id

    ! Local variables
    integer :: i, j, indTmp, iLev, fs
    real :: vP, vF, vertIndex, upWeight, levVal, lowerBorder, upperBorder
    real, dimension(:), pointer :: upPastPtr, downPastPtr, upFutPtr, downFutPtr, upPtr, downPtr
    type(silja_level) :: level

    !
    ! Stupidity check
    !
    if(varIndex > size(fld_buf%p4d) .or. varIndex < 1)then
      call msg('Index:', varIndex)
      call set_error('Strange index given','make_2d_from_4d_field')
      return
    endif
    if(.not.associated(fld_buf%p4d(varIndex)%past%p2d(1)%ptr) .or. &
     & .not.associated(fld_buf%p4d(varIndex)%future%p2d(1)%ptr)) then
      call msg('Index : ', varIndex)
      call set_error('Undefined fields for the index','make_2d_from_4d_field')
      return
    endif


    fs = fu_number_of_gridpoints(fu_grid(fld_buf%p4d(varIndex)%past%p2d(1)%IdPtr))
    if(size(dataPtr) < fs)then
      call set_error('Too small output data array','make_2d_from_4d_field')
      return
    endif

    level = fu_central_level_of_layer(levelTarget)
    if(error)return

    if(iVerticalTreatment == integrate_column_flag)then
      !
      ! Prepare the 3D cube to be integrated
      !
      do iLev = 1, fld_buf%nbr_of_levels
        call make_present_variable(fld_buf, varIndex, iLev, fs, now)
      end do
      
      if(error)return
      !
      ! We should take an integral / average over the whole column available in the buffer.
      ! levelTarget does not matter any more. The height field will give us the thickness of the layers
      !
      select case(fu_leveltype(verticalUsed))

        case(constant_pressure, hybrid, sigma_level, layer_btw_2_hybrid)
          !
          ! Find 3D height of pressure levels
          !
          indTmp = fu_index_for_quantity(fld_buf%buffer_quantities, height_flag) 
          if(indTmp < 1)then
            call set_error('Absent height field','make_2d_from_4d_field')
            return
          endif
          !
          ! Make present-time height
          !
          do iLev = 1, fld_buf%nbr_of_levels
            call make_present_variable(fld_buf, indTmp, iLev, fs, now)
          end do
          if(error)return
          !
          ! Scan the grid and integrate the columns
          !
          do i=1,fs
            dataPtr(i) = 0.0
            lowerBorder = 0.
            do iLev = 1, fld_buf%nbr_of_levels
              upperBorder = lowerBorder + &
                          & 2. * (fld_buf%p4d(indTmp)%present%p2d(iLev)%ptr(i) - lowerBorder)
              dataPtr(i) = dataPtr(i) + fld_buf%p4d(varIndex)%present%p2d(iLev)%ptr(i) * &
                                      & (upperBorder - lowerBorder)
              lowerBorder = upperBorder
            end do

          end do  ! fs

        case(constant_height) ! verticalUsed
          !
          ! Scan the grid and integrate the columns
          !
          do i=1,fs
            dataPtr(i) = 0.0
            lowerBorder = 0.
            do iLev = 1, fld_buf%nbr_of_levels
              upperBorder = lowerBorder + &
                          & 2. * (fu_level_height(fu_level(verticalUsed,iLev)) - lowerBorder)
              dataPtr(i) = dataPtr(i) + fld_buf%p4d(varIndex)%present%p2d(iLev)%ptr(i) * &
                                      & (upperBorder - lowerBorder)
              lowerBorder = upperBorder
            end do

          end do  ! fs

        case(layer_btw_2_height) ! verticalUsed
          !
          ! Scan the grid and integrate the columns
          !
          do i=1,fs
            dataPtr(i) = 0.0
            upperBorder = 0.
            do iLev = 1, fld_buf%nbr_of_levels
              dataPtr(i) = dataPtr(i) + fld_buf%p4d(varIndex)%present%p2d(iLev)%ptr(i) * &
                                      & fu_layer_thickness_m(fu_level(verticalUsed,iLev))
              upperBorder = upperBorder + fu_layer_thickness_m(fu_level(verticalUsed,iLev))
            end do

          end do  ! fs

        case default
          call set_error('Unknown level type of verticalUsed','make_2d_from_4d_field')
          call report(fu_level(verticalUsed,1))
          return
      end select  ! level type of vertical used
      !
      ! Finally, set the output id
      !
      idRes = fld_buf%p4d(varIndex)%present%p2d(fld_buf%nbr_of_levels)%idPtr
      if(iVerticalTreatment == integrate_column_flag)then
        call set_level(idRes, entire_atmosphere_integr_level)
      else
        call set_error('Strange vertical treatment (id):' + fu_str(iVerticalTreatment), &
                     & 'make_2d_from_4d_field')
      endif

    elseif(fu_leveltype(level) == fu_leveltype(verticalUsed) .or. &
         & fu_leveltype(level) == fu_leveltype(fu_level(verticalUsed, 1, .true.)))then
      !
      ! Here we also allow the verticalUsed to consist of layers
      ! between levels with the same type as level.
      !
      ! Level and meteo vertical coincide - one can get a universal coefficient
      ! for the vertical interpolation, which simplifies the whole process
      !
      ! levelTarget can be a single level or a thick layer. So far, this routine 
      ! supports only projection onto a single level and does not support
      ! vertical averaging over a thick layer. So, we have to transform a thick 
      ! layer into a single level placed in the middle of the layer.
      ! For a single level - nothing is to be done
      !
      vertIndex = fu_level_index(level, verticalUsed)
      if(error.or.vertIndex < 0)then
        call set_error('Failed to find the level vertical index','make_2d_from_4d_field')
        return
      endif

      if(vertIndex < 1.) then
        vertIndex = 1.

        ! For vertIndex=1 the time interpolation of level 1 gives immediate answer
        !
        call interpolate_fields_in_time(fld_buf%p4d(varIndex)%past%p2d(1)%idPtr, &
                                      & fld_buf%p4d(varIndex)%past%p2d(1)%ptr, &
                                      & fld_buf%p4d(varIndex)%future%p2d(1)%idPtr, &
                                      & fld_buf%p4d(varIndex)%future%p2d(1)%ptr, &
                                      & fld_buf%weight_past, &
                                      & idRes, & 
                                      & dataPtr, &
                                      & now)
        if(error)then
          call msg_test('NOW is: ')
          call report(now)
          call msg_test('Past id is:')
          call report(fld_buf%p4d(varIndex)%past%p2d(1)%idPtr)
          call msg_test('Future id is:')
          call report(fld_buf%p4d(varIndex)%future%p2d(1)%idPtr)
          return
        endif
        !??
        call set_level(idRes, levelTarget) ! level)

        ! ?? changend
      elseif(vertIndex >= fld_buf%nbr_of_levels)then
        vertIndex = fld_buf%nbr_of_levels

        ! Time interpolation of level nLevs gives immediate answer
        !
        call interpolate_fields_in_time(&
                     & fld_buf%p4d(varIndex)%past%p2d(fld_buf%nbr_of_levels)%idPtr, &
                     & fld_buf%p4d(varIndex)%past%p2d(fld_buf%nbr_of_levels)%ptr, &
                     & fld_buf%p4d(varIndex)%future%p2d(fld_buf%nbr_of_levels)%idPtr, &
                     & fld_buf%p4d(varIndex)%future%p2d(fld_buf%nbr_of_levels)%ptr, &
                     & fld_buf%weight_past, &
                     & idRes, & 
                     & dataPtr, &
                     & now)

        !??
        call set_level(idRes, levelTarget) !level)

        if(error)then
          call msg_test('NOW is: ')
          call report(now)
          call msg_test('Past id is:')
          call report(fld_buf%p4d(varIndex)%past%p2d(1)%idPtr)
          call msg_test('Future id is:')
          call report(fld_buf%p4d(varIndex)%future%p2d(1)%idPtr)
          return
        endif

      else

        ! Perform two time interpolations and then one vertical interpolation
        !
        downPtr => fu_work_array(fu_number_of_gridpoints( &
                                    & fu_grid(fld_buf%p4d(varIndex)%past%p2d(int(vertIndex))%idPtr)))
        call interpolate_fields_in_time(fld_buf%p4d(varIndex)%past%p2d(int(vertIndex))%idPtr, &
                                      & fld_buf%p4d(varIndex)%past%p2d(int(vertIndex))%ptr, &
                                      & fld_buf%p4d(varIndex)%future%p2d(int(vertIndex))%idPtr, &
                                      & fld_buf%p4d(varIndex)%future%p2d(int(vertIndex))%ptr, &
                                      & fld_buf%weight_past, &
                                      & idRes, & 
                                      & downPtr, &
                                      & now)


        call interpolate_fields_in_time(fld_buf%p4d(varIndex)%past%p2d(int(vertIndex+1))%idPtr, &
                                      & fld_buf%p4d(varIndex)%past%p2d(int(vertIndex+1))%ptr, &
                                      & fld_buf%p4d(varIndex)%future%p2d(int(vertIndex+1))%idPtr, &
                                      & fld_buf%p4d(varIndex)%future%p2d(int(vertIndex+1))%ptr, &
                                      & fld_buf%weight_past, &
                                      & idRes, & 
                                      & dataPtr, & ! Just temporarily - we'll overwrite it
                                      & now)
        if(error)then
          call msg_test('NOW is: ')
          call report(now)
          call msg_test('Past id is:')
          call report(fld_buf%p4d(varIndex)%past%p2d(1)%idPtr)
          call msg_test('Future id is:')
          call report(fld_buf%p4d(varIndex)%future%p2d(1)%idPtr)
          return
        endif

        upWeight = vertIndex-int(vertIndex) ! weight coef. for the upper level

        dataPtr(1:fs) = (dataPtr(1:fs)*upWeight + downPtr(1:fs)*(1.-upWeight))

        call set_level(idRes, levelTarget) !level)

        call free_work_array(downPtr)

      endif

    else
      !
      ! Brute-force way. For each grid cell the vertical interpolation
      ! is made for past and future with further time interpolation.
      ! So far, there are three possible level types, and they all are connected 
      ! each-to-each: constant pressure (p-system), constant height (z-system) and 
      ! hybrid. In fact, *-to-hybrid case cannot appear in current setup because
      ! output_server requires hybrid levels in meteoinput to copy the coefficients.
      ! But, when ini files themselves will contain the hybrid coefs - the situation 
      ! changes and this option will become possible.
      !
      ! levelTarget can be a single level or a thick layer. So far, this routine 
      ! supports only projection onto a single level and does not support
      ! vertical averaging over a thick layer. So, we have to transform a thick 
      ! layer into a single level placed in the middle of the layer.
      ! For a single level - nothing is to be done
      !
      level = fu_central_level_of_layer(levelTarget)
      if(error)return

      select case(fu_leveltype(verticalUsed))

        case(constant_pressure)

          select case(fu_leveltype(level))
            case(constant_height) ! From pressure to height (p -> z)
              !
              ! Find 3D height of pressure levels
              !
              indTmp = fu_index_for_quantity(fld_buf%buffer_quantities, height_flag) 
              if(fu_fails(indTmp > 0, 'Absent height field','make_2d_from_4d_field'))return
              do j = 1,fld_buf%nbr_of_levels
                call make_present_variable(fld_buf, indTmp, j, fs, now)
                call make_present_variable(fld_buf, varIndex, j, fs, now)
              end do ! cycle through levels of the present height 3D field
              
              levVal = fu_level_height(level) ! Store the height of the output level
              !
              ! Scan the grid
              !
              do i=1,fs
                !
                ! Find the output level index in the system vertical. 
                ! If corresponding level is not yet time-interpolated - do it
                !
                iLev=fld_buf%nbr_of_levels+1
                do j = 1,fld_buf%nbr_of_levels
                  if(levVal <= fld_buf%p4d(indTmp)%present%p2d(j)%ptr(i)) then
                    iLev = j
                    exit
                  endif
                end do ! cycle through levels of the present height 3D field
                !
                ! Compute the weight for the vertical interpolation and do it.
                !
                if(iLev==1)then ! Below the first level
                  iLev = 2
                  upWeight = 0.
                elseif(iLev > fld_buf%nbr_of_levels)then ! Above the last level
                  iLev = fld_buf%nbr_of_levels
                  upWeight = 1.
                else    ! Above the first, below the last
                  
                  upWeight = (levVal - fld_buf%p4d(indTmp)%present%p2d(ilev-1)%ptr(i)) / &
                           & (fld_buf%p4d(indTmp)%present%p2d(ilev)%ptr(i) - &
                            & fld_buf%p4d(indTmp)%present%p2d(ilev-1)%ptr(i))
                endif
                !
                ! OK, height is found and vertical interpolation is ready. Interpolate!
                !
                dataPtr(i) = fld_buf%p4d(varIndex)%present%p2d(iLev-1)%ptr(i)*(1.-upWeight)+ &
                           & fld_buf%p4d(varIndex)%present%p2d(iLev)%ptr(i)* upWeight
              end do  ! 1..fs

              idRes = fld_buf%p4d(varIndex)%present%p2d(iLev-1)%idPtr
              call set_level(idRes, levelTarget) !level)

            case(hybrid)
              call set_error('pressure-to-hybrid interpolation not supported so far',&
                           & 'make_2d_from_4d_field')
            case default
              call set_error('Unknown level type of level 1','make_2d_from_4d_field')
              call report(level)
              return
          end select ! leveltype of level

        case(constant_height) ! verticalUsed

          select case(fu_leveltype(level))
            case(constant_pressure)
              !
              ! Make present of 3D pressure of the height levels
              !
              indTmp = fu_index_for_quantity(fld_buf%buffer_quantities, pressure_flag) 
              if(fu_fails(indTmp > 0, 'Absent pressure field', 'make_2d_from_4d_field'))return
              do j = 1, fld_buf%nbr_of_levels
                call make_present_variable(fld_buf, indTmp, j, fs, now)
                call make_present_variable(fld_buf, varIndex, j, fs, now)
              end do ! 1..nLevs

              levVal = fu_pr_level_pressure(level) ! Store the height of the output level
              !
              ! Scan the grid
              !
              do i=1,fs
                !
                ! Find the output level index in the vertical for past fields
                ! Note: compare the central points of the levels
                ! Note2: pressure is decreasing with height
                !
                iLev = fld_buf%nbr_of_levels+1
                do j = 1,fld_buf%nbr_of_levels
                  if(levVal >= fld_buf%p4d(indTmp)%present%p2d(j)%ptr(i)) then
                    iLev = j
                    exit
                  endif
                end do ! 1..nLevs

                if(iLev==1)then ! Below the first level
                  iLev = 2
                  upWeight = 0.
                elseif(iLev > fld_buf%nbr_of_levels)then ! Above the last level
                  iLev = fld_buf%nbr_of_levels
                  upWeight = 1.
                else    ! Above the first, below the last
                  upWeight = (levVal - fld_buf%p4d(indTmp)%present%p2d(ilev-1)%ptr(i)) / &
                           & (fld_buf%p4d(indTmp)%present%p2d(ilev)%ptr(i) - &
                            & fld_buf%p4d(indTMp)%present%p2d(ilev-1)%ptr(i))
                endif

                !
                ! vertical interpolation
                !
                dataPtr(i) = fld_buf%p4d(varIndex)%present%p2d(iLev-1)%ptr(i)*(1.-upWeight)+ &
                           & fld_buf%p4d(varIndex)%present%p2d(iLev)%ptr(i) * upWeight
              end do

              idRes = fld_buf%p4d(varIndex)%present%p2d(iLev-1)%idPtr
              call set_level(idRes, levelTarget) !level)

            case(hybrid)
              call set_error('height-to-hybrid interpolation not supported so far',&
                           & 'make_2d_from_4d_field')
            case default
              call set_error('Unknown level type of level 2','make_2d_from_4d_field')
              call report(level)
              return
          end select ! leveltype of level

        case(hybrid, sigma_level, layer_btw_2_hybrid)  ! verticalUsed

          select case(fu_leveltype(level))
            case(constant_pressure)
              !
              ! 3D pressure of the hybrid levels:
              !
              indTmp = fu_index_for_quantity(fld_buf%buffer_quantities, pressure_flag) 
              if(fu_fails(indTmp > 0, 'Absent pressure field','make_2d_from_4d_field')) return
              do j = 1,fld_buf%nbr_of_levels
                call make_present_variable(fld_buf, indTmp, j, fs, now)
                call make_present_variable(fld_buf, varIndex, j, fs, now)
              end do ! 1..nLevs

              levVal = fu_pr_level_pressure(level) ! Store the height of the output level
              !
              ! Scan the grid
              !
              do i=1,fs
                !
                ! Find the output level index in the vertical for past fields
                ! Note: compare the central points of the levels
                ! Note2: pressure is decreasing with height
                !
                iLev = fld_buf%nbr_of_levels + 1
                do j = 1,fld_buf%nbr_of_levels
                  if(levVal >= fld_buf%p4d(indTmp)%present%p2d(j)%ptr(i)) then
                    iLev = j
                    exit
                  endif
                end do ! 1..nLevs

                if(iLev==1)then ! Below the first level
                  iLev = 2
                  upWeight = 0.
                elseif(iLev > fld_buf%nbr_of_levels)then ! Above the last level
                  iLev = fld_buf%nbr_of_levels
                  upWeight = 1.
                else    ! Above the first, below the last
                  upWeight = (levVal - fld_buf%p4d(indTmp)%present%p2d(ilev-1)%ptr(i)) / &
                           & (fld_buf%p4d(indTmp)%present%p2d(ilev)%ptr(i) - &
                            & fld_buf%p4d(indTmp)%present%p2d(ilev-1)%ptr(i))
                endif
                !
                ! vertical interpolation
                !
                dataPtr(i) = fld_buf%p4d(varIndex)%present%p2d(iLev-1)%ptr(i)*(1.-upWeight)+ &
                           & fld_buf%p4d(varIndex)%present%p2d(iLev)%ptr(i) * upWeight
              end do

              idRes = fld_buf%p4d(varIndex)%present%p2d(iLev-1)%idPtr
              call set_level(idRes, levelTarget) !level)

            case(constant_height)
              ! 3D height of the hybrid levels
              indTmp = fu_index_for_quantity(fld_buf%buffer_quantities, height_flag) 
              if(fu_fails(indTmp > 0, 'Absent height field','make_2d_from_4d_field')) return
              do j = 1,fld_buf%nbr_of_levels
                call make_present_variable(fld_buf, indTmp, j, fs, now)
                call make_present_variable(fld_buf, varIndex, j, fs, now)
              end do ! 1..nLevs
              levVal = fu_level_height(level) ! Store the height of the output level
              !
              ! Scan the grid
              !
              do i=1,fs
                !
                ! Find the output level index in the vertical for past fields
                ! Note: compare the central points of the levels
                !
               iLev = fld_buf%nbr_of_levels + 1
                do j = 1,fld_buf%nbr_of_levels
                  if(levVal <= fld_buf%p4d(indTmp)%present%p2d(j)%ptr(i)) then
                    iLev = j
                    exit
                  endif
                end do ! 1..nLevs

                if(iLev==1)then ! Below the first level
                  iLev = 2
                  upWeight = 0.
                elseif(iLev > fld_buf%nbr_of_levels)then ! Above the last level
                  iLev = fld_buf%nbr_of_levels
                  upWeight = 1.
                else    ! Above the first, below the last
                  upWeight = (levVal - fld_buf%p4d(indTmp)%present%p2d(ilev-1)%ptr(i)) / &
                           & (fld_buf%p4d(indTmp)%present%p2d(ilev)%ptr(i) - &
                            & fld_buf%p4d(indTmp)%present%p2d(ilev-1)%ptr(i))
                endif
                !
                ! vertical interpolation
                !
                dataPtr(i) = fld_buf%p4d(varIndex)%present%p2d(iLev-1)%ptr(i)*(1.-upWeight)+ &
                           & fld_buf%p4d(varIndex)%present%p2d(iLev)%ptr(i) * upWeight
              end do

              idRes = fld_buf%p4d(varIndex)%present%p2d(iLev-1)%idPtr
              call set_level(idRes, levelTarget) !level)

            case default ! verticalUsed
              call set_error('Unknown level type of level 3','make_2d_from_4d_field')
              call report(level)
              return
          end select ! leveltype of level

        case default
          call report(level)
          call set_error('Unknown level type of verticalUsed','make_2d_from_4d_field')
          call report(fu_level(verticalUsed,1))
          return
      end select

    endif ! If level types are the same

  end subroutine make_2d_from_4d_field


  ! ***************************************************************

  REAL FUNCTION fu_make_position_pressure(met_buf, &
                                        & pos, & 
                                        & weight_past, &
                                        & h_interp, v_interp, &
                                        & p_, t_, h_) RESULT(pos_pr)
  !
  ! Updates position's pressure if its height above the ground is known either
  ! in metric (position%height) or relative grid co-ordinates (position%z)
  ! Method: 
  ! Compute pressure difference from the closest level of 3D pressure field 
  ! using P1 = P2 * EXP( -1 * H * G / R * T). For that also the mean temperature 
  ! in the layer is needed. NOTE: simple interpolation of pressure field is not
  ! appropriate in case of thick layers.
  ! Order of checks:
  !  1. position%height
  !  2. position%z
  ! If neither OK - set_error. It does not matter if height is up-to-date
  ! 
  ! ATTENTION. Positions are always in the dispersion_grid
  !
  ! NO CHECKING.  Appropriate quantities must exist and the meteofield 
  !               pointers must be set
  !
  ! All units: SI
  !
  ! Author: Mikhail Sofiev
  !
    IMPLICIT NONE

    ! Imported parameters with the intent IN
    TYPE(Tfield_buffer),POINTER::met_buf
    TYPE(silam_grid_position),INTENT(in) :: pos
    INTEGER, INTENT(in) :: h_interp, v_interp ! methods of interpolation
    REAL, INTENT(in) :: weight_past
    TYPE(field_4d_data_ptr),INTENT(in),OPTIONAL:: p_,t_,h_

    ! Local variables
    INTEGER, DIMENSION(:), POINTER :: q_ptr
    INTEGER grd_ind, ind
!    INTEGER, DIMENSION(:), POINTER :: q_arr_ptr
    TYPE(field_data_ptr)::h_surf  
    TYPE(field_4d_data_ptr)::p_dat, t_dat, h_dat
    REAL :: t_mean, z, h, p1, h1
    
    IF(.not.defined(pos))THEN
      CALL set_error('undefined position','fu_update_position_pressure')
      RETURN
    END IF

    grd_ind = INT(fu_x(pos)+0.5) + (INT(fu_y(pos)+0.5)-1)*nx_dispersion
    z = fu_z(pos)
    h = fu_height(pos)
    h_surf%ptr => fu_work_array()
    h_surf%ptr =0.0 ! Represents zero height at the surface level

    IF(PRESENT(p_))THEN
      p_dat=p_
    ELSE
      ind = fu_index_for_quantity(met_buf%buffer_quantities, pressure_flag)
      p_dat = met_buf%p4d(ind)
    END IF

    IF(PRESENT(t_))THEN
      t_dat=t_
    ELSE
      ind = fu_index_for_quantity(met_buf%buffer_quantities, temperature_flag)
      t_dat = met_buf%p4d(ind)
    END IF

    IF(PRESENT(p_))THEN
      h_dat=h_
    ELSE
      ind = fu_index_for_quantity(met_buf%buffer_quantities, height_flag)
      h_dat = met_buf%p4d(ind)
    END IF

    IF (h<0.or.h>100000) THEN ! No height - take z
      IF (z<0.or.z>met_buf%nbr_of_levels) THEN ! Nothing usable - fail
        CALL set_error('Neither z nor height usable','fu_update_position_pressure')
        pos_pr = real_missing
        RETURN
      END IF
      !
      ! Compute height from z
      !
      h = fu_4d_interpolation(h_dat,pos,MAX(z,1.),weight_past,h_interp,v_interp)
      IF(z<1.)THEN ! Below the first layer - involve surface (h=0)
        h = z * h
      END IF
    ELSE ! height is usable - make z
      z = min(fu_vertical_index(fu_x(pos),fu_y(pos), fu_height(pos), &
                              & met_buf%nbr_of_levels, weight_past, &
                              & h_interp, &
                              & .true., & ! growing with index
                              & h_dat, h_surf), real(met_buf%nbr_of_levels))
    END IF
    CALL free_work_array(h_surf%ptr)
    !
    ! Now both height and z are available - compute mean temperature (roughly)
    ! Make it as an arithmetical average of temperature at position and at the
    ! lower(or 1-st if near-surface) layer
    !
    t_mean = 0.5* (fu_4d_interpolation(t_dat, & 
                                     & pos, &
                                     & MAX(z,1.), &
                                     & weight_past, &
                                     & h_interp, v_interp) + &
                 & fu_4d_interpolation(t_dat, & 
                                     & pos, &
                                     & REAL(MAX(INT(z),1)), &
                                     & weight_past, &
                                     & h_interp, v_interp))
    ! Reference values at the layer:
    !
    p1 = fu_4d_interpolation(p_dat, & 
                           & pos, &
                           & REAL(MAX(INT(z),1)), &
                           & weight_past, &
                           & h_interp, v_interp)
    h1 = fu_4d_interpolation(h_dat, & 
                           & pos, &
                           & REAL(MAX(INT(z),1)), &
                           & weight_past, &
                           & h_interp, v_interp)
    !
    ! All known - just take an exponent:
    !
    pos_pr = p1 * EXP((-1. * (h-h1) * g) / (gas_constant_dryair * t_mean))

  END FUNCTION fu_make_position_pressure



  ! ***************************************************************

  REAL FUNCTION fu_pressure_index_position (p_dat, surf_pr_dat, &
                                          & position, nLev, past_coef, interp_method ) !, wdr)
    !
    ! Finds the correct altitude index for a particular position 
    ! from the 4-D pressure cube.
    ! Stupid stuff! Should be used just until the vertical index
    ! is relative as well. To some extent it is similar to sounding,
    ! but hopefully a bit more efficient.
    !
    ! Method: find the interpolated between past and future sounding for 
    !         the exact position and then find vertical index as a linear
    !         combination of corresponding layers.
    !         If the return index <1 or >nLev - truncates.
    ! NOTE. In SILAM the vertical order is always from surface (layer 1) up. 
    !       This is assured by the level-comparison routines. See module levels.
    !
    ! ATTENTION. Positions are always in the dispersion_grid
    !
    IMPLICIT NONE

    ! Imported parameters with intent IN
    TYPE(silam_grid_position), INTENT(in) :: position
    REAL, INTENT(in) :: past_coef
    INTEGER, INTENT (in) ::  nLev, interp_method
!    TYPE(silja_wdr), INTENT(in) :: wdr

    ! Imported parameters with POINTER option
    TYPE(field_4d_data_ptr), INTENT(in) :: p_dat
    TYPE(field_data_ptr), INTENT(in) :: surf_pr_dat

    ! Local declarations
    REAL, DIMENSION(0:max_levels) :: p_vert
    INTEGER :: iLev !, pos_index
    REAL :: pos_pressure, surf_pr_pos

    p_vert = 1.0E10 ! Bigger than any pressure in Pa
    fu_pressure_index_position = -1    ! If search fails - return rabish
    pos_pressure = fu_pressure(position)

    DO iLev = 1, nLev
      !
      ! Find pressure at a particular vertical level but at position x,y
      ! interpoaltion routine requires the grid of the field, which is, by
      ! default, meteo_grid
      !
      IF(past_coef.eps.1.0)THEN
        p_vert(iLev) = fu_2d_interpolation(p_dat%past%p2d(iLev)%ptr, meteo_grid, &
                                               & position, interp_method)
      ELSE
        p_vert(iLev) = fu_2d_interpolation(p_dat%past%p2d(iLev)%ptr, meteo_grid, &
                                               & position, interp_method) * past_coef + &
                     & fu_2d_interpolation (p_dat%future%p2d(iLev)%ptr, meteo_grid, &
                                                & position, interp_method) * (1.-past_coef)
      END IF

      IF(pos_pressure < p_vert(iLev-1).and.pos_pressure >= p_vert(iLev)) THEN

        IF(iLev == 1)THEN ! Position is below 1st level, use surface pressure
          surf_pr_pos = fu_2d_interpolation(surf_pr_dat%ptr, meteo_grid, &
                                                & position, interp_method)
          IF(p_vert(1) >= surf_pr_pos - 1.E-3)THEN ! Problem with 3D pressure
            fu_pressure_index_position = 1.
          ELSE
            fu_pressure_index_position = (surf_pr_pos - MIN(pos_pressure,surf_pr_pos)) / &
                                       & (surf_pr_pos - p_vert(1))
          END IF
        ELSE ! Between some layers
          fu_pressure_index_position = real(iLev) - (pos_pressure - p_vert(iLev)) / &
                                       & (p_vert(iLev-1) - p_vert(iLev))
        END IF

!        ! Actually not needed, but kept so far for debugging purposes
!        IF(fu_pressure_index < 0. .or. fu_pressure_index > nLev)THEN
!          call msg('')
!          call msg('Strange output pressure index:', real_value=fu_pressure_index)
!          CALL set_error('Failed pressure indexing','fu_pressure_index')
!        END IF
        RETURN
      END IF

    END DO

    ! All layers passed, but pressure index not found: very high particle
    ! Still, if it is reasonably close to the top, continue, otherwise set it outside
    !
!    CALL msg_warning('Very high particle','fu_pressure_index')
!    CALL report(position)
    if(pos_pressure >= p_vert(nLev) + 0.5*(p_vert(nLev) - p_vert(nLev-1)))then
      fu_pressure_index_position = nLev  ! still OK
    else
      fu_pressure_index_position = nLev+2.  ! too much, break the story
    endif

  END FUNCTION fu_pressure_index_position


  ! ***************************************************************

  REAL FUNCTION fu_pressure_index_coord(p_dat, surf_pr_dat, fx, fy, pos_pressure, &
                                      & nLev, past_coef, interp_method) !, wdr)
    !
    ! Finds the correct altitude index for a particular position 
    ! from the 4-D pressure cube.
    ! Stupid stuff! Should be used just until the vertical index
    ! is relative as well. To some extent it is similar to sounding,
    ! but hopefully a bit more efficient.
    !
    ! Method: find the interpolated between past and future sounding for 
    !         the exact position and then find vertical index as a linear
    !         combination of corresponding layers.
    !         If the return index <1 or >nLev - truncates.
    ! NOTE. In SILAM the vertical order is always from surface (layer 1) up. 
    !       This is assured by the level-comparison routines. See module levels.
    !
    ! ATTENTION. Positions are always in the dispersion_grid
    !
    IMPLICIT NONE

    ! Imported parameters with intent IN
    real, INTENT(in) :: fx, fy, pos_pressure
    REAL, INTENT(in) :: past_coef
    INTEGER, INTENT (in) ::  nLev, interp_method
!    TYPE(silja_wdr), INTENT(in) :: wdr

    ! Imported parameters with POINTER option
    TYPE(field_4d_data_ptr), INTENT(in) :: p_dat
    TYPE(field_data_ptr), INTENT(in) :: surf_pr_dat

    ! Local declarations
    REAL, DIMENSION(0:max_levels) :: p_vert
    INTEGER :: iLev !, pos_index
    REAL :: surf_pr_pos

    p_vert = 1.0E10 ! Bigger than any pressure in Pa
    fu_pressure_index_coord = -1    ! If search fails - return rabish

    DO iLev = 1, nLev
      !
      ! Find pressure at a particular vertical level but at position x,y
      ! interpoaltion routine requires the grid of the field, which is, by
      ! default, meteo_grid
      !(grid_data, x,y, nx,ny, method, iOutside)
      IF(past_coef.eps.1.0)THEN
        p_vert(iLev) = fu_2d_interpolation(p_dat%past%p2d(iLev)%ptr, fx, fy, &
                                               & nx_meteo, ny_meteo, interp_method, nearestPoint)
      ELSE
        p_vert(iLev) = fu_2d_interpolation(p_dat%past%p2d(iLev)%ptr, fx, fy, &
                                    & nx_meteo, ny_meteo, interp_method, nearestPoint) * past_coef + &
                     & fu_2d_interpolation (p_dat%future%p2d(iLev)%ptr, fx, fy, &
                                    & nx_meteo, ny_meteo, interp_method, nearestPoint) * (1.-past_coef)
      END IF

      IF(pos_pressure < p_vert(iLev-1).and.pos_pressure >= p_vert(iLev)) THEN

        IF(iLev == 1)THEN ! Position is below 1st level, use surface pressure
          surf_pr_pos = fu_2d_interpolation(surf_pr_dat%ptr, fx, fy, &
                                                & nx_meteo, ny_meteo, interp_method, nearestPoint)
          IF(p_vert(1) >= surf_pr_pos - 1.E-3)THEN ! Problem with 3D pressure
            fu_pressure_index_coord = 1.
          ELSE
            fu_pressure_index_coord = (surf_pr_pos - MIN(pos_pressure, surf_pr_pos)) / &
                                    & (surf_pr_pos - p_vert(1))
          END IF
        ELSE ! Between some layers
          fu_pressure_index_coord = real(iLev) - (pos_pressure - p_vert(iLev)) / &
                                  & (p_vert(iLev-1) - p_vert(iLev))
        END IF

!        ! Actually not needed, but kept so far for debugging purposes
!        IF(fu_pressure_index < 0. .or. fu_pressure_index > nLev)THEN
!          call msg('')
!          call msg('Strange output pressure index:', real_value=fu_pressure_index)
!          CALL set_error('Failed pressure indexing','fu_pressure_index')
!        END IF
        RETURN
      END IF

    END DO

    ! All layers passed, but pressure index not found: very high particle
    ! Still, if it is reasonably close to the top, continue, otherwise set it outside
    !
!    CALL msg_warning('Very high particle','fu_pressure_index')
!    CALL report(position)
    if(pos_pressure >= p_vert(nLev) + 0.5*(p_vert(nLev) - p_vert(nLev-1)))then
      fu_pressure_index_coord = nLev  ! still OK
    else
      fu_pressure_index_coord = nLev+2.  ! too much, break the story
    endif

  END FUNCTION fu_pressure_index_coord


  ! ***************************************************************

  REAL FUNCTION fu_vertical_index (x,y, val, &
                                 & nLev, past_coef, &
                                 & interp_method, &
                                 & ifIncr, &
                                 & ptr_4d, ptr_2d)
  !
  ! Finds the correct altitude index for a particular position 
  ! from the 4-D cube for ANY variable. If exists, involves 2d field
  ! field as an underlying (surface) value
  ! To some extent it is similar to sounding, but more efficient.
  !
  ! Method: find the interpolated between past and future incomplete sounding 
  !         for the exact position and then find vertical index as a linear
  !         combination of corresponding layers.
  !         If the return index <1 or >nLev - truncates.
  ! NOTE. In SILAM the vertical order is always from surface (layer 1) up. 
  !       This is assured by the level-comparison routines. See module levels.
  !
  IMPLICIT NONE

  ! Imported parameters with intent IN
!  TYPE(Tfield_buffer),POINTER::met_buf
  REAL, INTENT(in) :: x,y,val ! Actually determines 3D position
  REAL, INTENT(in) :: past_coef
  INTEGER, INTENT (in) ::  nLev, interp_method
  LOGICAL, INTENT(in) :: ifIncr ! Way of sorting .true.=> growing upward

  ! Imported parameters with POINTER option
  TYPE(field_4d_data_ptr), INTENT(in) :: ptr_4d
  TYPE(field_data_ptr), INTENT(in),OPTIONAL :: ptr_2d
  
  ! Local declarations
  REAL, DIMENSION(0:max_levels) :: v_vert
  INTEGER :: iLev
  REAL :: pos_pressure, value_surf
!  TYPE(silam_grid_position) :: pos

!  pos = fu_set_position(x,y,real_missing,time_missing)

!  IF(ifIncr)THEN
!    v_vert = -1.0E20  ! Smaller than everything
!  ELSE
!    v_vert = 1.0E20 ! Bigger than anything
!  END IF

  fu_vertical_index = -1    ! If search fails - return rabish

!  PRINT *, pos_pressure, pos_index, surf_pressure_ptr%ptr(pos_index)

  DO iLev = 1, nLev
    !
    ! Find value at a particular vertical level but at position x,y
    !
    IF(past_coef.eps.1.0)THEN
      v_vert(iLev) = fu_2d_interpolation(ptr_4d%past%p2d(iLev)%ptr, &
                                             & x, y, nx_meteo, ny_meteo, &
                                             & interp_method, nearestPoint)
    ELSE
      v_vert(iLev) = fu_2d_interpolation(ptr_4d%past%p2d(iLev)%ptr, &
                                             & x, y, nx_meteo, ny_meteo, &
                                             & interp_method, nearestPoint) * past_coef + &
                   & fu_2d_interpolation(ptr_4d%future%p2d(iLev)%ptr, &
                                             & x, y, nx_meteo, ny_meteo, &
                                             & interp_method, nearestPoint) * (1.-past_coef)
    END IF

    !
    ! Check if value is between iLev-1 and iLev
    !
!    ifBetween = value < v_vert(iLev) .neqv. ifIncr

    IF((val > v_vert(iLev)) .neqv. ifIncr) THEN

      IF(iLev == 1)THEN ! Position is below 1st level, use 2d field if any
        fu_vertical_index = 1
        IF(PRESENT(ptr_2d))THEN
          value_surf = fu_2d_interpolation(ptr_2d%ptr, x, y, nx_meteo, ny_meteo, &
                                               & interp_method, nearestPoint)
          IF((v_vert(1) < value_surf-1.E-3) .neqv. ifIncr)THEN ! Fields OK
            IF((val < value_surf) .neqv. ifIncr)THEN ! Above the surface
              fu_vertical_index = (value_surf - val) / (value_surf - v_vert(1))
            END IF
          END IF
        END IF

      ELSE ! Between some layers
        fu_vertical_index = real(iLev) - (val - v_vert(iLev)) / &
                               & (v_vert(iLev-1) - v_vert(iLev))
      END IF

      ! Actually not needed, but kept so far for debugging purposes
      IF(fu_vertical_index < 0. .or. fu_vertical_index > nLev)THEN
        CALL msg('Failed vertical indexing:',fu_vertical_index)
        CALL set_error('Failed vertical indexing','fu_vertical_index')
      END IF
      RETURN
    END IF

  END DO

  ! All layers passed, but vertical index not found: very high particle
  !
  CALL msg_warning('Very high particle','fu_vertical_index')
  fu_vertical_index = nLev

  if(val >= v_vert(nLev) + 0.5*(v_vert(nLev) - v_vert(nLev-1)))then
    fu_vertical_index = nLev  ! still OK
  else
    fu_vertical_index = nLev+2.  ! too much, break the story
  endif

  END FUNCTION fu_vertical_index


  !*****************************************************************************

  subroutine make_present_variable(met_buf, varIndex, iLev, iSize, now)
    !
    ! Performs allocation and time interpolation, unless it is already done
    ! Also performs memory allocation, if needed
    !
    implicit none

    ! Imported parameters
    TYPE(Tfield_buffer), intent(inout)::met_buf
    integer, intent(in) :: varIndex, iLev, iSize
    type(silja_time), intent(in) :: now
    integer :: iTmp

    ! Local variables
    integer :: al_status

    !
    ! Stupidity check (future only for weigh_past < 1)
    !
    if(.not.met_buf%p4d(varIndex)%past%p2d(iLev)%ifReady)then
      call msg('')
      call msg('Meteo variable index: ',varIndex)
      call msg('Level:',iLev)
      call set_error('Level does not exist in past','make_present_variable')
      return
    endif
    if(.not.(met_buf%weight_past.eps.1.))then
      if(.not.met_buf%p4d(varIndex)%future%p2d(iLev)%ifReady)then
        call msg('')
        call msg('Meteo variable index: ',varIndex)
        call msg('Level: ',iLev)
        call msg('Quantity:' + &
              & fu_quantity_string(fu_quantity(met_buf%p4d(varIndex)%past%p2d(iLev)%idPtr)))
        call set_error('Level does not exist in future','make_present_variable')
        return
      endif
    endif

    if(.not. met_buf%p4d(varIndex)%present%p2d(iLev)%ifReady)then
      ! One more stupidity check
      ! Check if past and future are the same field (true for single-time stack)
      if(associated(met_buf%p4d(varIndex)%past%p2d(iLev)%ptr, &
                  & target=met_buf%p4d(varIndex)%future%p2d(iLev)%ptr)) then
        call msg("Past and future are the same:")
        call msg('--> Meteo variable index: ',varIndex)
        call msg('--> Level: ',iLev)
        call msg('--> Quantity:' + &
                  & fu_quantity_string(fu_quantity(met_buf%p4d(varIndex)%past%p2d(iLev)%idPtr)))
        if (associated(met_buf%p4d(varIndex)%present%p2d(iLev)%ptr, &
                     & target=met_buf%p4d(varIndex)%future%p2d(iLev)%ptr)) then
          call msg("Setting present ready")
          met_buf%p4d(varIndex)%present%p2d(iLev)%ifReady = .true.
          call msg_warning("Present was not set ready for singletime field", &
                               & "make_present_variable")
        else
          call msg("Present is different. Don't know what to do...")
          met_buf%p4d(varIndex)%present%p2d(iLev)%ifReady = .true.
          call msg_warning("Present was not set ready for singletime field", &
                               & "make_present_variable")
        endif
      else !past and future are different. Make the present...

        !
        !  Create the field if not exists...
        if(.not. associated(met_buf%p4d(varIndex)%present%p2d(iLev)%ptr))then
          allocate(met_buf%p4d(varIndex)%present%p2d(iLev)%idPtr, &
                 & met_buf%p4d(varIndex)%present%p2d(iLev)%ptr(iSize), stat=al_status)
          if(al_status /= 0)then
            call set_error('Failed to allocate memory','make_present_variable')
            return
          endif
        endif

        call interpolate_fields_in_time(met_buf%p4d(varIndex)%past%p2d(iLev)%idPtr, &
                                      & met_buf%p4d(varIndex)%past%p2d(iLev)%ptr, &
                                      & met_buf%p4d(varIndex)%future%p2d(iLev)%idPtr, &
                                      & met_buf%p4d(varIndex)%future%p2d(iLev)%ptr, &
                                      & met_buf%weight_past, &
                                      & met_buf%p4d(varIndex)%present%p2d(iLev)%idPtr, &
                                      & met_buf%p4d(varIndex)%present%p2d(iLev)%ptr, &
                                      & now)
        met_buf%p4d(varIndex)%present%p2d(iLev)%ifReady = .not.error
        if(error)then

          call msg_test('NOW is: ')
          call report(now)
          call msg('')
          call msg('Meteo variable index: ',varIndex)
          call msg('Level: ',iLev)
          call msg(fu_connect_strings('Quantity:', &
                 & fu_quantity_string(fu_quantity(met_buf%p4d(varIndex)%past%p2d(iLev)%idPtr))))
          call set_error('Level does not exist in future','make_present_variable')
          call msg("Level",iLev)
          call msg_test('Past id is:')
          call report(met_buf%p4d(varIndex)%past%p2d(iLev)%idPtr)
          call msg_test('Future id is:')
          call report(met_buf%p4d(varIndex)%future%p2d(iLev)%idPtr)
          return
        endif
      endif ! check for realtime field with unset present...
    endif ! present field ready

  end subroutine make_present_variable


  !*******************************************************************************
  !*******************************************************************************
  !*******************************************************************************
  !
  !
  !  Interpolation structures and operations with them
  !
  !
  !*******************************************************************************
  !*******************************************************************************
  !*******************************************************************************


  !*******************************************************************************
  !
  !  A common routines for both structures
  !
  !*******************************************************************************


  real function fu_get_value_from_bfr_interp(buf, &
                                           & quantityIndex, &
                                           & nxFrom, ixTo, iyTo, iLevTo, &
                                           & pHorizInterpStruct, pVertInterpStruct, &
                                           & ifHorizInterp, ifVertInterp, &
                                           & ifForceWeightPast_) result(value)
    !
    ! The very basic function: returns a value asked to a specific To-grid cell
    ! Made without any checking because it is supposed to be inside main cycles
    !
    implicit none

    ! Imported parameters
    type(Tfield_buffer), pointer :: buf
    type(THorizInterpStruct), intent(in) :: pHorizInterpStruct
    type(TVertInterpStruct), intent(in) :: pVertInterpStruct
    integer, intent(in) :: quantityIndex, nxFrom, ixTo, iyTo, iLevTo
    logical, intent(in) :: ifHorizInterp, ifVertInterp
    logical, intent(in), optional :: ifForceWeightPast_

    ! Local variables
    logical :: ifForceWeightPast

    if(present(ifForceWeightPast_))then
      ifForceWeightPast = ifForceWeightPast_
    else
      ifForceWeightPast = .false.
    endif

    !
    ! Actually, just decides whether this quantity is 2d or 4d and then calls the
    ! appropriate downstream function. Slow but repetition of the extraction looks too ugly
    !
    if(buf%p4d(quantityIndex)%past%p2d(1)%ifReady)then
      value = fu_get_value_from_bfr_interp4d(buf%p4d(quantityIndex), &
                                           & nxFrom, ixTo, iyTo, iLevTo, &
                                           & buf%weight_past, &
                                           & pHorizInterpStruct, pVertInterpStruct, &
                                           & ifHorizInterp, ifVertInterp)
    else
      value = fu_get_value_from_fldp2d_int3d(buf%p2d(quantityIndex), nxFrom, ixTo, iyTo, &
                                           & buf%weight_past, &
                                           & pHorizInterpStruct, ifHorizInterp, ifForceWeightPast)
    endif

  end function fu_get_value_from_bfr_interp


  !*******************************************************************************

  real function fu_get_value_from_bfr_interp4d(field_4d_from, nxFrom, ixTo, iyTo, iLevTo, &
                                             & weight_past, &
                                             & pHorizInterpStruct, pVertInterpStruct, &
                                             & ifHorizInterp, ifVertInterp) result(fValue)
    !
    ! The very basic function: returns a value asked to a specific To-grid cell
    ! Made without any checking because it is supposed to be inside all main cycles
    !
    implicit none

    ! Imported parameters
    type(field_4d_data_ptr), intent(in) :: field_4d_from
    type(THorizInterpStruct), intent(in) :: pHorizInterpStruct
    type(TVertInterpStruct), intent(in) :: pVertInterpStruct
    integer, intent(in) :: nxFrom, ixTo, iyTo, iLevTo
    real, intent(in) :: weight_past
    logical, intent(in) :: ifHorizInterp, ifVertInterp

    ! Local variables
    integer :: iCellIndex, iCoef, iLevFrom, iLev
    real :: fTmp
!    type(TVertInterpOneCell), pointer :: pCoefsVert
    type(THorizInterpOneCell) :: pCoefsHoriz
    real, dimension(max_levels) :: fTmpFrom
    
    !pHor => HorizInterpStruct
    !pVert => VertInterpStruct
    !
    ! ATTENTION. ABSOLUTELY NO CHECKING TO KEEP THE MAX SPEED
    !
    if(ifHorizInterp)then
      !
      ! Full-blown horizontal interpolation
      !
      if(ifVertInterp)then
        !
        ! Full-blown vertical & horizontal interpolations
        !
!        pCoefsVert => pVertInterpStruct%coefs(ixTo,iyTo,iLevTo)
        call get_coefs(pHorizInterpStruct,ixTo,iyTo, pCoefsHoriz)

        do iLevFrom = pVertInterpStruct%indLev(1,ixTo,iyTo,iLevTo), &
                    & pVertInterpStruct%indLev(pVertInterpStruct%nCoefs,ixTo,iyTo,iLevTo)
          if(iLevFrom < 1 .or. iLevFrom > size(fTmpFrom))then
            call msg('Strange levFrom for iLevTo:',iLevFrom,iLevTo)
            call msg('ixTo,iyTo:',ixTo,iyTo)
            call msg('Max number of coefs:',pVertInterpStruct%nCoefs)
            call msg('Levels to scan, from-to:',pVertInterpStruct%indLev(1,ixTo,iyTo,iLevTo), &
                               & pVertInterpStruct%indLev(pVertInterpStruct%nCoefs,ixTo,iyTo,iLevTo))
            call set_error('Strange levFrom', 'fu_get_value_from_bfr_interp4d')
            return
          endif
          fTmpFrom(iLevFrom) = 0.
          do iCoef = 1, pHorizInterpStruct%nCoefs
            iCellIndex = pCoefsHoriz%indX(iCoef) + (pCoefsHoriz%indY(iCoef) - 1.) * nxFrom
            fTmpFrom(iLevFrom) = fTmpFrom(iLevFrom) + pCoefsHoriz%weight(iCoef) * &
                & (field_4d_from%past%p2d(iLevFrom)%ptr(iCellIndex) * weight_past + &
                &  field_4d_from%future%p2d(iLevFrom)%ptr(iCellIndex) * (1.-weight_past))
          end do
        end do ! Filling-up the vertFrom column but for gridTo
        !
        ! The horizontal interpolation made the vertical column in vertFrom
        ! Now it can be interpolated to vertTo
        !
        fValue = 0.
        do iCoef = 1, pVertInterpStruct%nCoefs
          fValue = fValue + pVertInterpStruct%weight(iCoef,ixTo,iyTo,iLevTo) * &
                          & fTmpFrom(pVertInterpStruct%indLev(iCoef,ixTo,iyTo,iLevTo))
        end do

      else
        !
        ! Full horizontal but no vertical interpolation
        ! 
        call get_coefs(pHorizInterpStruct,ixTo,iyTo, pCoefsHoriz)

        fValue = 0.
        do iCoef = 1, pHorizInterpStruct%nCoefs
          iCellIndex = pCoefsHoriz%indX(iCoef) + (pCoefsHoriz%indY(iCoef) - 1.) * nxFrom
          fValue = fValue + pCoefsHoriz%weight(iCoef) * &
              & (field_4d_from%past%p2d(iLevTo)%ptr(iCellIndex) * weight_past + &
              &  field_4d_from%future%p2d(iLevTo)%ptr(iCellIndex) * (1.-weight_past))
        end do
      endif  ! if vertical interpolation needed

    else
      !
      ! No horizontal interpolation
      !
      if(ifVertInterp)then
        !
        ! No horizontal but full vertical interpolation. Have to sum-up along the vertical.
        ! LevelsFrom to be considered are taken from %indLev value, which is sorted
        ! in ascending order
        !
!        if ( iLevTo < 1 .or. iLevTo > size(pVertInterpStruct%weight, 4) ) then
!
!           call msg('fu_get_value: vertFrom:')
!           call report(pVertInterpStruct%vertFrom, .true.)
!           call msg('')
!           call msg('fu_get_value: vertTo:')
!           call report(pVertInterpStruct%vertTo, .true.)
!           call msg('levels from the 4d field:')
!           do iCoef=1,fu_nbrOfLevels(pVertInterpStruct%vertFrom)
!             call report(fu_level(field_4d_from%past%p2d(iCoef)%idPtr))
!           end do
!           call set_error("Gotcha!", "sdgfaerghewtbh")
!           call msg("pVertInterpStruct%weight(:,ixTo,iyTo,iLevTo)",pVertInterpStruct%weight(:,ixTo,iyTo,iLevTo))
!           call msg("pVertInterpStruct%indLev(1,ixTo,iyTo,iLevTo)",pVertInterpStruct%indLev(1,ixTo,iyTo,iLevTo))
!           call msg("pVertInterpStruct%indLev(2,ixTo,iyTo,iLevTo)",pVertInterpStruct%indLev(2,ixTo,iyTo,iLevTo))
!
!        endif

        iCellIndex = ixTo+(iyTo-1)*nxFrom
!        pCoefsVert => pVertInterpStruct%coefs(ixTo,iyTo,iLevTo)
        fValue = 0.
        do iCoef = 1, pVertInterpStruct%nCoefs
          fTmp = pVertInterpStruct%weight(iCoef,ixTo,iyTo,iLevTo)
          iLev = pVertInterpStruct%indLev(iCoef,ixTo,iyTo,iLevTo)
          fValue = fValue + fTmp*( field_4d_from%past%p2d(iLev)%ptr(iCellIndex) * weight_past +&
                    &            field_4d_from%future%p2d(iLev)%ptr(iCellIndex)*(1.-weight_past))
        end do
      else
        !
        ! No interpolation at all: To and From grids and verticals are identical
        !
        iCellIndex = ixTo+(iyTo-1)*nxFrom
        fValue = field_4d_from%past%p2d(iLevTo)%ptr(iCellIndex) * weight_past + &
              & field_4d_from%future%p2d(iLevTo)%ptr(iCellIndex) * (1.-weight_past)
      endif  ! if vertical itnerpolation needed

    endif  ! if horizontal interpolation needed

  end function fu_get_value_from_bfr_interp4d


  !********************************************************************************

  real function fu_get_value_from_fldp2d_int3d(field_p2d_from, nxFrom, ixTo, iyTo, &
                                             & weight_past, &
                                             & pHorizInterpStruct, ifHorizInterp, &
                                             & ifForceWeightPast_) result(fValue)
    !
    ! The very basic function: returns a value asked to a specific To-grid cell
    ! Made without any checking because it is supposed to be inside all main cycles
    !
    implicit none

    ! Imported parameters
    type(field_2d_data_ptr), intent(in) :: field_p2d_from
    type(THorizInterpStruct), intent(in) :: pHorizInterpStruct
    integer, intent(in) :: nxFrom, ixTo, iyTo
    real, intent(in) :: weight_past
    logical, intent(in) :: ifHorizInterp
    logical, intent(in), optional :: ifForceWeightPast_

    ! Local variables
    integer :: iCellIndex, iCoef
    type(THorizInterpOneCell) :: pCoefsHoriz
    logical :: ifForceWeightPast

    !
    ! ATTENTION. ABSOLUTELY NO CHECKING TO KEEP THE MAX SPEED
    !
    if(present(ifForceWeightPast_))then
      ifForceWeightPast = ifForceWeightPast_
    else
      ifForceWeightPast = .false.
    endif

    if(ifHorizInterp)then
      !
      ! Full-blown horizontal interpolation
      !
      call get_coefs(pHorizInterpStruct,ixTo,iyTo, pCoefsHoriz)

      fValue = 0.
      if((.not.ifForceWeightPast) .and. field_p2d_from%present%ifReady)then
        do iCoef = 1, pHorizInterpStruct%nCoefs
          iCellIndex = pCoefsHoriz%indX(iCoef) + (pCoefsHoriz%indY(iCoef) - 1.) * nxFrom
          fValue = fValue + pCoefsHoriz%weight(iCoef) * field_p2d_from%present%ptr(iCellIndex)
        end do
      else
        do iCoef = 1, pHorizInterpStruct%nCoefs
          iCellIndex = pCoefsHoriz%indX(iCoef) + (pCoefsHoriz%indY(iCoef) - 1.) * nxFrom
          fValue = fValue + pCoefsHoriz%weight(iCoef) * &
                & (field_p2d_from%past%ptr(iCellIndex) * weight_past + &
                &  field_p2d_from%future%ptr(iCellIndex) * (1.-weight_past))
        end do
      endif

    else
      !
      ! No horizontal interpolation
      !
      if((.not.ifForceWeightPast) .and. field_p2d_from%present%ifReady)then
        fValue = field_p2d_from%present%ptr(ixTo+(iyTo-1)*nxFrom)
      else
        iCellIndex = ixTo+(iyTo-1)*nxFrom
        fValue = field_p2d_from%past%ptr(iCellIndex) * weight_past + &
               & field_p2d_from%future%ptr(iCellIndex) * (1.-weight_past)
      endif

    endif  ! if horizontal interpolation needed

  end function fu_get_value_from_fldp2d_int3d


  !********************************************************************************

  real function fu_get_value_from_fld2d_int(field_2d_from, nxFrom, ixTo, iyTo, &
                                          & pHorizInterpStruct, ifHorizInterp) result(fValue)
    !
    ! The very basic function: returns a value asked to a specific To-grid cell
    ! Made without any checking because it is supposed to be inside all main cycles
    !
    implicit none

    ! Imported parameters
    type(silja_field), pointer :: field_2d_from
    type(THorizInterpStruct), intent(in) :: pHorizInterpStruct
    integer, intent(in) :: nxFrom, ixTo, iyTo
    logical, intent(in) :: ifHorizInterp

    ! Local variables
    integer :: iCellIndex, iCoef
    type(THorizInterpOneCell) :: pCoefsHoriz
    real, dimension(:), pointer :: pData

    !
    ! ATTENTION. ABSOLUTELY NO CHECKING TO KEEP THE MAX SPEED
    !
    if(ifHorizInterp)then
      !
      ! Full-blown horizontal interpolation
      !
      call get_coefs(pHorizInterpStruct,ixTo,iyTo, pCoefsHoriz)
      pData => fu_grid_data(field_2d_from)

      fValue = 0.
      do iCoef = 1, pHorizInterpStruct%nCoefs
        iCellIndex = pCoefsHoriz%indX(iCoef) + (pCoefsHoriz%indY(iCoef) - 1.) * nxFrom
        fValue = fValue + pCoefsHoriz%weight(iCoef) * pData(iCellIndex)
      end do

    else
      !
      ! No horizontal interpolation
      !
      fValue = pData(ixTo+(iyTo-1)*nxFrom)

    endif  ! if horizontal interpolation needed

  end function fu_get_value_from_fld2d_int


  !***************************************************************************************************
  
  subroutine wind_from_buffer_4d(field_data_u, field_data_v, nxFrom, ixTo, iyTo, iLevTo, weight_past, &
                               & pHorizInterpStructU, pVertInterpStructU, &
                               & pHorizInterpStructV, pVertInterpStructV, &
                               & ifHorizInterp, ifVertInterp, u, v)

    !
    ! Get values (u,v) of a windfield at point ixTo, iyTo,
    ! iLevTo. This subroutine is a wrapper around
    ! fu_value_from_buffer_interp_4d, with the rotation applied as needed.
    ! 
    ! The winds fields u and v can have separate interpolation
    ! structures (due to Arakawa shifts). It is implicitly assumed
    ! that they have the same poles (and hence same rotation coeffecients).

    implicit none
    type(field_4d_data_ptr), intent(in) :: field_data_u, field_data_V
    type(THorizInterpStruct), pointer :: pHorizInterpStructU, pHorizInterpStructV
    type(TVertInterpStruct), pointer :: pVertInterpStructU, pVertInterpStructV
    integer, intent(in) :: nxFrom, ixTo, iyTo, iLevTo
    real, intent(in) :: weight_past
    logical, intent(in) :: ifHorizInterp, ifVertInterp

    real, intent(out) :: u, v
    
    real :: uTmp, vTmp
    logical :: ifRotate

    if (.not. ifHorizInterp) then
      ifRotate = .false.
    else
      ifRotate = pHorizInterpStructU%ifRotation
    end if

    if (ifRotate) then
      if (.not. (associated(pHorizInterpStructU%rotation)) .and. associated(pHorizInterpStructV%rotation)) then
        call set_error('Wind rotation required but not associated', 'wind_from_buffer_4d')
        return
      end if
      uTmp = fu_get_value_from_bfr_interp4d(field_data_u, nxFrom, ixTo, iyTo, iLevTo, weight_past, &
                                          & pHorizInterpStructU, pVertInterpStructU, &
                                          & ifHorizInterp, ifVertInterp)
      vTmp = fu_get_value_from_bfr_interp4d(field_data_v, nxFrom, ixTo, iyTo, iLevTo, weight_past, &
                                          & pHorizInterpStructV, pVertInterpStructV, &
                                          & ifHorizInterp, ifVertInterp)
      u = uTmp*pHorizInterpStructU%rotation(1,1,ixTo,iyTo) + vTmp*pHorizInterpStructU%rotation(1,2,ixTo,iyTo)
      v = uTmp*pHorizInterpStructU%rotation(2,1,ixTo,iyTo) + vTmp*pHorizInterpStructU%rotation(2,2,ixTo,iyTo)
    else
      u = fu_get_value_from_bfr_interp4d(field_data_u, nxFrom, ixTo, iyTo, iLevTo, weight_past, &
                                       & pHorizInterpStructU, pVertInterpStructU, &
                                       & ifHorizInterp, ifVertInterp)
      v = fu_get_value_from_bfr_interp4d(field_data_v, nxFrom, ixTo, iyTo, iLevTo, weight_past, &
                                       & pHorizInterpStructV, pVertInterpStructV, &
                                       & ifHorizInterp, ifVertInterp)
    end if
      

  end subroutine wind_from_buffer_4d


  !*******************************************************************************

  subroutine interp_fld3d_field_2_fld3d(field_3d_in, field_3d_out, &
                                      & ifHorizInterp, ifVertInterp, &
                                      & pHorizInterpStruct, pVertInterpStruct)
    !
    ! A basic task: squeeze one 3d field to another one using interpolation structures
    ! Made as fast as possible, so checking is limited.
    ! No time interpolation
    !
    implicit none

    ! Imported parameters
    type(silja_3d_field), pointer :: field_3d_in
    type(silja_3d_field), pointer :: field_3d_out
    type(TVertInterpStruct), pointer :: pVertInterpStruct
    type(THorizInterpStruct), pointer :: pHorizInterpStruct
    logical, intent(in) :: ifHorizInterp, ifVertInterp
    character(len=*), parameter :: sub_name="interp_fld3d_field_2_fld3d"

    !
    ! Largely resembles the activity of fu_value_from_buffer_interp_3d but for different input type
    !
    ! Local variables
    integer :: iCellIndex, iCoef, iLevFrom, ixTo, iyTo, nxTo, nyTo, iLevTo, nxFrom, nyFrom
!    type(TVertInterpOneCell), pointer :: pCoefsVert
    type(THorizInterpCells) :: pCoefsHoriz
    real, dimension(:), pointer :: fTmpFrom, pDataTo
    type(field_3d_data_ptr) :: pFld3dIn
    real :: fValue

    ! Still some checking
    if(ifHorizInterp)then
      if(.not. associated(pHorizInterpStruct))then
        call set_error('Horizontal interpolation is requested but the structure is undefined', &
                     & sub_name)
        return
      endif
    else  
      ! No horizontal interpolation - the horizontal interpolation structure can be undefined
      !
      if(.not. fu_grid(field_3d_in) == fu_grid(field_3d_out))then
        call set_error('No horizontal interpolation requested but grids mismatch',  sub_name)
      endif
    endif

! The chenck below would be the right thing to do, but too many places cuts this corner making field_3d_in an field_3d_out ooff sync

!     if(.not. fu_quantity(field_3d_in) == fu_quantity(field_3d_out))then
!       call set_error("Quantities do not match", sub_name)
!     endif
! 
!     if(.not. fu_species(field_3d_in) == fu_species(field_3d_out))then
!       call set_error("Species do not match", sub_name)
!     endif

    if (error) then
        call msg('Input 3d field:')
        call report(field_3d_in)
        call msg('Output 3d field:')
        call report(field_3d_out)
        call set_error('No horizontal interpolation requested but grids mismatch',  sub_name)
        return
     endif


    !
    ! Arrange direct access to the input field
    !
    do iLevFrom = 1, fu_NbrOfLevels(pVertInterpStruct%vertFrom)
      pFld3dIn%p2d(iLevFrom)%ptr => fu_grid_data_from_3d(field_3d_in, iLevFrom)
    end do

    ! The actual work is done by a subroutine using the direct access to the input fields

    call interp_fld3d_data_2_fld3d(pFld3dIn, field_3d_out, &
                                      & ifHorizInterp, ifVertInterp, &
                                      & pHorizInterpStruct, pVertInterpStruct)

  end subroutine interp_fld3d_field_2_fld3d


  !*******************************************************************************

  subroutine interp_fld3d_data_2_fld3d(field_3d_data_in, field_3d_out, &
                                     & ifHorizInterp, ifVertInterp, &
                                     & pHorizInterpStruct, pVertInterpStruct)
    !
    ! A basic task: squeeze one 3d field structure to another one using interpolation structures
    ! Made as fast as possible, so checking is limited.
    ! No time interpolation
    !
    implicit none

    ! Imported parameters
    type(field_3d_data_ptr), intent(in) :: field_3d_data_in
    type(silja_3d_field), pointer :: field_3d_out
    type(TVertInterpStruct), pointer :: pVertInterpStruct
    type(THorizInterpStruct), pointer :: pHorizInterpStruct
    logical, intent(in) :: ifHorizInterp, ifVertInterp
    !
    ! Largely resembles the activity of fu_value_from_buffer_interp_3d but for different input type
    !
    ! Local variables
    integer :: iCellIndex, iCoef, iLevFrom, ixTo, iyTo, nxTo, nyTo, iLevTo, nxFrom, nyFrom
    type(THorizInterpCells) :: pCoefsHoriz
    real, dimension(:), pointer :: fTmpFrom, pDataTo
    real :: fValue

    if(ifHorizInterp)then
      if(.not. associated(pHorizInterpStruct))then
        call set_error('Horizontal interpolation is requested but the structure is undefined', &
                     & 'interp_fld3d_data_2_fld3d')
        return
      endif
      call grid_dimensions(pHorizInterpStruct%gridTo, nxTo, nyTo)
      call grid_dimensions(pHorizInterpStruct%gridFrom, nxFrom, nyFrom)
    else  
      ! No horizontal interpolation - the horizontal interpolation structure can be undefined.
      ! Output field defines the grid
      !
      call grid_dimensions(fu_grid(field_3d_out), nxTo, nyTo)
      nxFrom = nxTo
      nyFrom = nyTo
    endif
    if(error)return
    fTmpFrom => fu_work_array()
    if(error)return

    !
    ! Now make the interpolation 
    !
    if(ifHorizInterp)then
      !
      ! Full-blown horizontal interpolation
      !
      if(ifVertInterp)then
        !
        ! Full-blown vertical & horizontal interpolations
        !
        do iLevTo = 1, fu_NbrOfLevels(pVertInterpStruct%vertTo)

          pDataTo => fu_grid_data_from_3d(field_3d_out, iLevTo)
          call get_coefs(pHorizInterpStruct, pCoefsHoriz)
          if(error)return

          
          do iyTo = 1, nyTo
            do ixTo = 1, nxTo
!              pCoefsVert => pVertInterpStruct%coefs(ixTo,iyTo,iLevTo)

              do iLevFrom = pVertInterpStruct%indLev(1,ixTo,iyTo,iLevTo), &
                          & pVertInterpStruct%indLev(pVertInterpStruct%nCoefs,ixTo,iyTo,iLevTo)
                fTmpFrom(iLevFrom) = 0.
                do iCoef = 1, pHorizInterpStruct%nCoefs
                  iCellIndex = pCoefsHoriz%indX(iCoef,ixTo,iyTo) + &
                             & (pCoefsHoriz%indY(iCoef,ixTo,iyTo) - 1.) * nxFrom
                  fTmpFrom(iLevFrom) = fTmpFrom(iLevFrom) + pCoefsHoriz%weight(iCoef,ixTo,iyTo) * &
                                                    & field_3d_data_in%p2d(iLevFrom)%ptr(iCellIndex)
                end do
              end do ! Filling-up the vertFrom column but for gridTo
              !
              ! The horizontal interpolation made the vertical column in vertFrom
              ! Now it can be interpolated to vertTo
              !
              fValue = 0.
              do iCoef = 1, pVertInterpStruct%nCoefs
                fValue = fValue + pVertInterpStruct%weight(iCoef,ixTo,iyTo,iLevTo) * &
                                & fTmpFrom(pVertInterpStruct%indLev(iCoef,ixTo,iyTo,iLevTo))
              end do
              pDataTo(ixTo + (iyTo-1)*nyTo) = fValue

            end do ! ixTo
          end do ! iyTo
        end do ! iLevTo

      else
        !
        ! Full horizontal but no vertical interpolation
        ! 
        do iLevTo = 1, fu_NbrOfLevels(pVertInterpStruct%vertTo)
          pDataTo => fu_grid_data_from_3d(field_3d_out, iLevTo)
          call get_coefs(pHorizInterpStruct, pCoefsHoriz)
          if(error)return

          do iyTo = 1, nyTo
            do ixTo = 1, nxTo
              fValue = 0.
              do iCoef = 1, pHorizInterpStruct%nCoefs
                iCellIndex = pCoefsHoriz%indX(iCoef,ixTo,iyTo) + (pCoefsHoriz%indY(iCoef,ixTo,iyTo) - 1.) * nxFrom
                fValue = fValue + pCoefsHoriz%weight(iCoef,ixTo,iyTo) * &
                                & field_3d_data_in%p2d(iLevTo)%ptr(iCellIndex)
              end do
              pDataTo(ixTo + (iyTo-1)*nyTo) = fValue
            end do ! ixTo
          end do ! iyTo
        end do ! iLevTo
      endif  ! if vertical interpolation needed

    else
      !
      ! No horizontal interpolation
      !
      if(ifVertInterp)then
        !
        ! No horizontal but full vertical interpolation. Have to sum-up along the vertical.
        ! LevelsFrom to be considered are taken from %indLev value, which is sorted
        ! in ascending order
        !

!        call msg('fu_get_value: vertFrom:')
!        call report(pVertInterpStruct%vertFrom, .true.)
!        call msg('')
!        call msg('fu_get_value: vertTo:')
!        call report(pVertInterpStruct%vertTo, .true.)
!        call msg('levels from the 4d field:')
!        do iCoef=1,fu_nbrOfLevels(pVertInterpStruct%vertFrom)
!          call report(fu_level(field_4d_from%past%p2d(iCoef)%idPtr))
!        end do

        do iLevTo = 1, fu_NbrOfLevels(pVertInterpStruct%vertTo)
          pDataTo => fu_grid_data_from_3d(field_3d_out, iLevTo)
          if(error)return
          do iyTo = 1, nyTo
            do ixTo = 1, nxTo
              iCellIndex = ixTo+(iyTo-1)*nxTo
!              pCoefsVert => pVertInterpStruct%coefs(ixTo,iyTo,iLevTo)
              fValue = 0.
              do iCoef = 1, pVertInterpStruct%nCoefs
                fValue = fValue + &
                       & pVertInterpStruct%weight(iCoef,ixTo,iyTo,iLevTo) * &
                       & field_3d_data_in%p2d(pVertInterpStruct%indLev(iCoef,ixTo,iyTo,iLevTo))% &
                                                                                 & ptr(iCellIndex)
              end do
              pDataTo(iCellIndex) = fValue
            end do ! ixTo
          end do ! iyTo
        end do ! iLevTo
      else
        !
        ! No interpolation at all: To and From grids and verticals are identical
        !
        do iLevTo = 1, fu_NbrOfLevels(pVertInterpStruct%vertTo)
          pDataTo => fu_grid_data_from_3d(field_3d_out, iLevTo)
          if(error)return
          do iyTo = 1, nyTo
            do ixTo = 1, nxTo
              iCellIndex = ixTo+(iyTo-1)*nxTo
              pDataTo(iCellIndex) = field_3d_data_in%p2d(iLevTo)%ptr(iCellIndex)
            end do ! ixTo
          end do ! iyTo
        end do ! iLevTo
      endif  ! if vertical itnerpolation needed

    endif  ! if horizontal interpolation needed

    call free_work_array(fTmpFrom)

  end subroutine interp_fld3d_data_2_fld3d



  !*******************************************************************************
  !
  !  Vertical interpolation structure
  !
  !*******************************************************************************


  !******************************************************************

  function fu_vertical_interp_struct(vertFrom, vertTo, gridInterp, interpolation_method, &
                                   & recommended_update_interval, chName, &
                                   & data_buffer, now, &                  ! Meteodata buffer. OPTIONAL
                                   & p_interpHorizStruct4Meteo, p_InterpVert4Meteo) &  ! OPTIONAL
                                   & result(interpVertStruct)
    !
    ! Computes the interpolation defined in the interpVertStruct
    ! ATTENTION! The interpolation method, vertFrom and vertTo
    ! must be given. This routine uses this information to determine
    ! the coefs structure array.
    ! The coefficients must be in accending order, so that max and min indices
    ! are evident without min/max calls
    !
    implicit none

    ! Imported parameters
    type(silam_vertical), intent(in) :: vertFrom, vertTo
    type(silja_grid), intent(in) :: gridInterp ! In which the structure is defined
    integer, intent(in) :: interpolation_method
    type(silja_interval), intent(in) :: recommended_update_interval
    character(len=*), intent(in) :: chName
    type(Tfield_buffer), pointer, optional :: data_buffer  ! optional, if immediate refinement
    type(silja_time), intent(in), optional :: now          ! optional, if immediate refinement
    type(TVertInterpStruct), pointer, optional :: p_InterpVert4Meteo ! optional, if refinement now
    type(THorizInterpStruct), pointer, optional :: p_interpHorizStruct4Meteo ! optional, if refinement now

    ! returned pointer
    type(TVertInterpStruct), pointer :: interpVertStruct

    ! Local variables
    integer :: nLevTo, nLevFrom, status, ix,iy, i, nxInterp, nyInterp, num_quant_req, iCount, &
             & maxCount, iLevFrom, iLevTo
    real :: fLev, zx, zx1, zx2, zx3, zx4
    real, dimension(4) :: weights
    integer, dimension(4) :: indices
    logical :: ifVertsComparable
    type(silja_level) :: levFromTmp ! Level from the vertFrom but projected to VertTo 
    integer, dimension(:), pointer :: quant_req
    real, dimension(:,:), pointer :: weightTmp, weight_Z_Tmp

    !
    ! Stupidity check
    !
    if(.not. (defined(vertFrom) .and. defined(vertTo)) .and. defined(gridInterp))then
      call msg('vertFrom:')
      call report(vertFrom)
      call msg('vertTo')
      call report(vertTo)
      call msg('Reprojection grid:')
      call report(gridInterp)
      call set_error('One of the above is undefined','fu_vertical_interp_struct')
      return
    endif

    !
    ! First, try to find an old structure. If found, check that it is updated frequently enough
    !
!    do iCount = 1, VertInterpPool%nVIS
!      if(trim(VertInterpPool%pVIS(iCount)%ptr%chName) == trim(chName))then
!        interpVertStruct => VertInterpPool%pVIS(iCount)%ptr
!call msg('Returning the requested interpolation structure:' + chName)
!        return
!      endif
!    end do
    do iCount = 1, VertInterpPool%nVIS
      if(VertInterpPool%pVIS(iCount)%ptr%interp_type == interpolation_method)then
        if(VertInterpPool%pVIS(iCount)%ptr%grid == gridInterp)then
          if(fu_cmp_verts_eq(VertInterpPool%pVIS(iCount)%ptr%vertFrom, vertFrom))then
            if(fu_cmp_verts_eq(VertInterpPool%pVIS(iCount)%ptr%vertTo, vertTo))then
              interpVertStruct => VertInterpPool%pVIS(iCount)%ptr
              if(interpVertStruct%recommended_update_interval > recommended_update_interval) &
                        & interpVertStruct%recommended_update_interval = recommended_update_interval
if(trim(VertInterpPool%pVIS(iCount)%ptr%chName) == trim(chName))then
  call msg('Returning the interpolation structure:' + interpVertStruct%chName)
else
  call msg('Returning the interpolation structure:' + interpVertStruct%chName + ', instead of:' + chName)
endif
              return
            endif
          endif
        endif
      endif
    end do
    !
    ! Strcutre is not found. Create new one!
    !
call msg('Making new vertical interpolation structure:' + chName, VertInterpPool%nVIS+1)
call msg('vertFrom:')
call report(vertFrom)
call msg('vertTo')
call report(vertTo)
call msg('Reprojection grid:')
call report(gridInterp)
call msg('')
    if(VertInterpPool%nVIS == size(VertInterpPool%pVIS))then
      call set_error('No more space for vertical interpolation structures','fu_vertical_interp_struct')
      return
    endif
    VertInterpPool%nVIS = VertInterpPool%nVIS + 1
    allocate(VertInterpPool%pVIS(VertInterpPool%nVIS)%ptr, stat = iCount)
    if(iCount /= 0)then
      call set_error('Failed allocation of next horizontal interpolation strcutre','fu_vertical_interp_struct')
      return
    endif
    interpVertStruct => VertInterpPool%pVIS(VertInterpPool%nVIS)%ptr

    interpVertStruct%chName = chName
    interpVertStruct%recommended_update_interval = recommended_update_interval
    interpVertStruct%vertFrom = vertFrom
    interpVertStruct%vertTo = vertTo
    interpVertStruct%interp_type = interpolation_method
    interpVertStruct%grid = gridInterp
    interpVertStruct%ifMeteoGrid = (gridInterp == meteo_grid)
    interpVertStruct%lastCoefUpdate = really_far_in_past
    if(error)return

    !
    ! Allocate the map of coefs
    !
    nLevFrom = fu_NbrOfLevels(interpVertStruct%vertFrom)
    nLevTo = fu_NbrOfLevels(interpVertStruct%vertTo)
    call grid_dimensions(gridInterp,nxInterp, nyInterp)
    if(error)return

    interpVertStruct%iLevStartFrom = nLevFrom
    interpVertStruct%iLevEndFrom = 1
    interpVertStruct%iLevStartTo = nLevTo
    interpVertStruct%iLevEndTo = 1

    ifVertsComparable = fu_verts_comparable(interpVertStruct%vertTo, interpVertStruct%vertFrom)

    !
    ! Allocate the necessary pointers and make up the reprojection coefficients
    ! A trick: vertical interpolation is the same for all grid cells because we do not
    ! know actual meteorology, neither orography.
    !
    select case(interpVertStruct%interp_type)

      case(average, summation)
        !
        ! These are the flux-preserving reprojection ways.
        ! Note that the number of cells affecting the single output one is unknown. Have to use 
        ! temporaries as long as needed.
        ! A relief, however, is that there is no bulky 2D geometry: overlap of layers and its
        ! centre can be computed easily without splitting into small sub-layers.
        !
        ! First, check that we are dealing with layers.
        !
        if(.not. (fu_if_layer(fu_level(interpVertStruct%vertFrom,1)) .and. &
                & fu_if_layer(fu_level(interpVertStruct%vertTo,1))))then
          call msg_warning('Only layers allowed for sum/ave interpolation','fu_vertical_interp_struct')
          call msg('Vertical from:')
          call report(vertFrom)
          call msg('Vertical to:')
          call report(vertTo)
          call set_error('Only layers allowed for sum/ave interpolation','fu_vertical_interp_struct')
          return
        endif
        !
        ! Let's get temporary array that will almost surely incorporate up to 
        !
        weightTmp => fu_work_array_2d(nLevFrom,nLevTo)
        weight_Z_Tmp => fu_work_array_2d(nLevFrom,nLevTo)
        if(error)return
        maxCount = 0
        !
        ! Actual work starts. First do a brute-force: count overlaps of all layers with all layers
        !
        do iLevTo = 1, nLevTo
          !
          ! overlap and momentum of its centre:
          !
          iCount = 0
          do iLevFrom = 1, nLevFrom
            weightTmp(iLevFrom, iLevTo) = fu_vert_overlap_fraction( &
                                                    & fu_level(interpVertStruct%vertFrom,iLevFrom), &
                                                    & fu_level(interpVertStruct%vertTo,iLevTo), &
                                                    & overlap_centre_idx_lyr_with = zx)
            weight_Z_Tmp(iLevFrom, iLevTo) = (zx - 1.) * weightTmp(iLevFrom,iLevTo)

            if(.not. (weightTmp(iLevFrom, iLevTo) .eps. 0.0))then
              iCount = iCount + 1
              interpVertStruct%iLevStartFrom = min(interpVertStruct%iLevStartFrom, iLevFrom)
              interpVertStruct%iLevEndFrom = max(interpVertStruct%iLevEndFrom, iLevFrom)
              interpVertStruct%iLevStartTo = min(interpVertStruct%iLevStartTo, iLevTo)
              interpVertStruct%iLevEndTo = max(interpVertStruct%iLevEndTo, iLevTo)
            endif
          end do  ! iLev To
          maxCount = max(maxCount, iCount)
        enddo  ! iLevFrom

        interpVertStruct%nCoefs = maxCount
        allocate(interpVertStruct%indLev(maxCount,nxInterp,nyInterp,nLevTo), &
               & interpVertStruct%weight(maxCount,nxInterp,nyInterp,nLevTo), &
               & interpVertStruct%weight_Z(maxCount,nxInterp,nyInterp,nLevTo), stat=status)
        if(status /= 0)then
          call msg('Coefs allocation status:',status)
          call set_error('Failed to allocate memory for coefs(iLev)','fu_vertical_interp_struct')
          return
        endif
        !
        ! Now find out the max number of levelsFrom affecting the levelsTo: the coef dimension
        !
        do iy=1,nyInterp
          do ix=1,nxInterp
            do iLevTo = 1, nLevTo
              interpVertStruct%indLev(1:maxCount,ix,iy,iLevTo) = 0
              interpVertStruct%weight(1:maxCount,ix,iy,iLevTo) = 0.0
              interpVertStruct%weight_Z(1:maxCount,ix,iy,iLevTo) = 0.0
              iCount = 1
              do iLevFrom = 1, nLevFrom
                if(weightTmp(iLevFrom,iLevTo) .eps. 0.0)cycle
                interpVertStruct%indLev(iCount,ix,iy,iLevTo) = iLevFrom
                interpVertStruct%weight(iCount,ix,iy,iLevTo) = weightTmp(iLevFrom,iLevTo)
                interpVertStruct%weight_Z(iCount,ix,iy,iLevTo) = weight_Z_Tmp(iLevFrom,iLevTo)
                iCOunt = iCount + 1
              enddo
            end do
          end do    ! ix
        end do   ! iy

        call free_work_array(weightTmp)
        call free_work_array(weight_Z_Tmp)        
        
      case(nearest_point)
        !
        ! coefs vectors are 1-element pointing to the nearest point
        !
        interpVertStruct%nCoefs = 1
        allocate(interpVertStruct%indLev(1,nxInterp,nyInterp,nLevTo), &
               & interpVertStruct%weight(1,nxInterp,nyInterp,nLevTo), &
               & interpVertStruct%weight_Z(1,nxInterp,nyInterp,nLevTo), stat=status)
        if(status /= 0)then
          call msg('Coefs allocation status:',status)
          call set_error('Failed to allocate memory for coefs(iLev)','fu_vertical_interp_struct')
          return
        endif
        interpVertStruct%iLevStartTo = 1
        interpVertStruct%iLevEndTo = nLevTo

        do iLevTo = 1,nLevTo
          !
          ! Project each levelTo to the vertFrom and then check the index value
          !
          if(ifVertsComparable)then
            fLev = fu_level_index(fu_level(interpVertStruct%vertTo,iLevTo), &
                                & interpVertStruct%vertFrom)
          else
            levFromTmp = fu_level_to_vertical_crude(fu_level(interpVertStruct%vertTo, &
                                                           & iLevTo, &   ! level index
                                                           & .true.), &  ! single level needed, not layer
                                                  & interpVertStruct%vertFrom)
            fLev = fu_level_index(levFromTmp, interpVertStruct%vertFrom)
          endif
          if(error)return

          indices(1) = min(nLevFrom, max(1,nint(fLev)))
          interpVertStruct%iLevStartFrom = min(interpVertStruct%iLevStartFrom, indices(1))
          interpVertStruct%iLevEndFrom = max(interpVertStruct%iLevEndFrom, indices(1))
          interpVertStruct%weight(:,:,:,:) = 1.

          do iy=1,nyInterp
            do ix=1,nxInterp
              interpVertStruct%indLev(1,ix,iy,iLevTo) = indices(1)
            end do    ! ix
          end do   ! iy
        end do  ! iLev

      case(linear)
        !
        ! Two surrounding levels elements
        !
        interpVertStruct%nCoefs = 2
        allocate(interpVertStruct%indLev(2,nxInterp,nyInterp,nLevTo), &
               & interpVertStruct%weight(2,nxInterp,nyInterp,nLevTo), &
               & interpVertStruct%weight_Z(2,nxInterp,nyInterp,nLevTo), stat=status)
        if(status /= 0)then
          call msg('Coefs allocation status:',status)
          call set_error('Failed to allocate memory for coefs(iLev)','fu_vertical_interp_struct')
          return
        endif
        interpVertStruct%iLevStartTo = 1
        interpVertStruct%iLevEndTo = nLevTo

        do iLevTo = 1,nLevTo
          !
          ! Project each levelTo to the vertFrom and then check the index value
          !
          if(ifVertsComparable)then
            fLev = fu_level_index(fu_level(interpVertStruct%vertTo,iLevTo), &
                                & interpVertStruct%vertFrom)
          else
            levFromTmp = fu_level_to_vertical_crude(fu_level(interpVertStruct%vertTo, &
                                                           & iLevTo, &  ! level index
                                                           & .true.), &  ! singe level needed, not layer
                                                  & interpVertStruct%vertFrom)
            fLev = fu_level_index(levFromTmp, interpVertStruct%vertFrom)
          endif
          if(error)return

          zx = fLev - REAL(INT(fLev)) ! Exceedance of fLev over lower border

          indices(1) = min(nLevFrom,max(1,int(fLev)))
          weights(1) = 1. - zx
          indices(2) = min(nLevFrom,max(1,int(fLev)+1))
          weights(2) = zx
          interpVertStruct%iLevStartFrom = min(interpVertStruct%iLevStartFrom, indices(1))
          interpVertStruct%iLevEndFrom = max(interpVertStruct%iLevEndFrom, indices(2))

          do iy=1,nyInterp
            do ix=1,nxInterp
              interpVertStruct%indLev(1:2,ix,iy,iLevTo) = indices(1:2)
              interpVertStruct%weight(1:2,ix,iy,iLevTo) = weights(1:2)
            end do
          end do
        end do


      case(cubic)
        !
        ! Four surrounding grid elements
        !
        interpVertStruct%nCoefs = 4
        allocate(interpVertStruct%indLev(4,nxInterp,nyInterp,nLevTo), &
               & interpVertStruct%weight(4,nxInterp,nyInterp,nLevTo), &
               & interpVertStruct%weight_Z(4,nxInterp,nyInterp,nLevTo), stat=status)
        if(status /= 0)then
          call msg('Coefs allocation status:',status)
          call set_error('Failed to allocate memory for coefs(iLev)','fu_vertical_interp_struct')
          return
        endif
        interpVertStruct%iLevStartTo = 1
        interpVertStruct%iLevEndTo = nLevTo

        do iLevTo = 1, nLevTo
          !
          ! Project each levelTo to the vertFrom and then check the index value
          !
          if(ifVertsComparable)then
            fLev = fu_level_index(fu_level(interpVertStruct%vertTo,iLevTo), &
                                & interpVertStruct%vertFrom)
          else
            levFromTmp = fu_level_to_vertical_crude(fu_level(interpVertStruct%vertTo, &
                                                           & iLevTo, &  ! level index
                                                           & .true.), & ! single level needed, not layer
                                                  & interpVertStruct%vertFrom)
            fLev = fu_level_index(levFromTmp, interpVertStruct%vertFrom)
          endif
          zx = fLev - REAL(INT(fLev)) ! Exceedance of fLev over lower border

          indices(1) = min(nLevFrom,max(1,int(fLev)-1))
          weights(1) = ((-0.5*zx+1.0)*zx-0.5)*zx 
          indices(2) = min(nLevFrom,max(1,int(fLev)))
          weights(2) = (( 1.5*zx-2.5)*zx    )*zx+1 
          indices(3) = min(nLevFrom,max(1,int(fLev)+1))
          weights(3) = ((-1.5*zx+2.0)*zx+0.5)*zx 
          indices(4) = min(nLevFrom,max(1,int(fLev)+2))
          weights(4) = (( 0.5*zx-0.5)*zx    )*zx 
          interpVertStruct%iLevStartFrom = min(interpVertStruct%iLevStartFrom, indices(1))
          interpVertStruct%iLevEndFrom = max(interpVertStruct%iLevEndFrom, indices(4))

          do iy=1,nyInterp
            do ix=1,nxInterp
              interpVertStruct%indLev(1:4,ix,iy,iLevTo) = indices(1:4)
              interpVertStruct%weight(1:4,ix,iy,iLevTo) = weights(1:4)
            end do
          end do
        end do

      case default
        call msg('Unknown interpolation method:',interpVertStruct%interp_type)
        call set_error('Unknown interpolation method','fu_vertical_interp_struct')
        return
    end select

    !
    ! A user-friendly piece. The above structure definition uses crude interpolation.
    ! But the user may wish to use it immediately. Then we have refine it right here.
    !
    if(present(data_buffer))then
      if(fu_ifMeteoDependentInterp(interpVertStruct))then
        if(defined(interpVertStruct%recommended_update_interval))then
          if(interpVertStruct%recommended_update_interval < one_day)then
            if(fu_fails(present(now),'time is not present but refinement requested', &
                                                             & 'fu_vertical_interp_struct'))return
            call refine_interp_vert_coefs_v2(interpVertStruct, data_buffer, now)
          endif  ! update interval is shorter than 1 day
        endif  ! defined update itnerval
      endif ! meteo-dependent interpolation
    endif  ! data_buffer pointer given

  end function fu_vertical_interp_struct


  !******************************************************************

  logical function fu_ifMeteoDependentInterp(interpStructVert)
    !
    ! If the verticals are significantly different, the interpolation from
    ! one to another will require actual meteodata, thus becoming dynamic.
    ! This procedure checks this question.
    !
    implicit none
    type(TVertInterpStruct), intent(in) :: interpStructVert

    integer :: q3d, q2d

    !
    ! Is the interpolation trivial?
    !
    if(fu_leveltype(interpStructVert%vertFrom) == fu_leveltype(interpStructVert%vertTo))then
      fu_ifMeteoDependentInterp = .false.
      return
    endif

    ! Can just use the input needs subroutine from levels.
    call projection_input_needs(interpStructVert%vertTo, interpStructVert%vertFrom, q3d, q2d)
    if (q3d == int_missing .and. q2d == int_missing) then
      fu_ifMeteoDependentInterp = .false.
    else
      fu_ifMeteoDependentInterp = .true.
    end if

  end function fu_ifMeteoDependentInterp


  !************************************************************************************

  subroutine refine_all_vert_interp_coefs(meteo_buf_ptr, ifNewMeteoData, now)
    !
    ! Refines coefficients ALL VERTICAL INTERPOLATION STRUCTURES. 
    ! Tries to do it smartly: if not meteo-dependent, skips the whole thing,
    ! sometimes avoid refinement if it is important only when new meteo oarrives, etc.
    !
    implicit none
    
    ! Imported parameters
    type(Tfield_buffer), intent(in) :: meteo_buf_ptr
    logical, intent(in) :: ifNewMeteoData
    type(silja_time), intent(in) :: now

    ! Local vairables
    integer :: iIS

    !
    ! Actually, just goes through the pool of vertical structures and refines them all
    ! by calling below refine_interp_vert_coefs_v2
    !
    do iIS = 1, VertInterpPool%nVIS  ! actually defined ones
      !
      ! There may be no dynamic in relation between the verticals, do nothing if so.
      ! Also, coefficients might have been already updated
      !
      if(.not.fu_ifMeteoDependentInterp(VertInterpPool%pVIS(iIS)%ptr))cycle

      if(defined(VertInterpPool%pVIS(iIS)%ptr%lastCoefUpdate))then
        if(now == VertInterpPool%pVIS(iIS)%ptr%lastCoefUpdate)then
          call msg(VertInterpPool%pVIS(iIS)%ptr%chName + &
                                     & ': interpolation structure coefficients are up to date')
          cycle
        endif
        if(fu_abs(now - VertInterpPool%pVIS(iIS)%ptr%lastCoefUpdate) < &
                                  & VertInterpPool%pVIS(iIS)%ptr%recommended_update_interval)then
          call msg(VertInterpPool%pVIS(iIS)%ptr%chName + &
                                     & ': interpolation structure coeffs are almost up to date')
          cycle
        endif
      endif


!      print *, "Before:"
!      print *, "Ilev", VertInterpPool%pVIS(iIS)%ptr%indLev(1:2,10,10,4)
!      print *, "Coef", VertInterpPool%pVIS(iIS)%ptr%weight(1:2,10,10,4)
      call refine_interp_vert_coefs_v2(VertInterpPool%pVIS(iIS)%ptr, & ! Structure to refine
                                     & meteo_buf_ptr, &    ! Meteo buffer, pointers are ready
                                     & now)
!      print *, "After:"
!      print *, "Ilev", VertInterpPool%pVIS(iIS)%ptr%indLev(1:2,10,10,4)
!      print *, "Coef", VertInterpPool%pVIS(iIS)%ptr%weight(1:2,10,10,4)


!      if(VertInterpPool%pVIS(iIS)%ptr%grid == meteo_grid)then
!        if(fu_cmp_verts_eq(VertInterpPool%pVIS(iIS)%ptr%vertFrom, meteo_vertical))then
!call msg('refining in meteo_grid and meteo_vertical:',iIS)
!          call refine_interp_vert_coefs_v2(VertInterpPool%pVIS(iIS)%ptr, & ! Structure to refine
!                                         & meteo_buf_ptr, &      ! Meteo buffer, pointers are ready
!                                         & now)
!        else
!call msg('refining in meteo_grid but non-meteo_vertical:',iIS)
!          call refine_interp_vert_coefs_v2( &
!                           & VertInterpPool%pVIS(iIS)%ptr, & ! Structure to refine
!                           & meteo_buf_ptr, &      ! Meteo buffer, pointers are ready
!                           & now, &
!                           & p_InterpVert4Meteo = &
!                                  & fu_vertical_interp_struct(meteo_vertical, &
!                                                            & VertInterpPool%pVIS(iIS)%ptr%vertFrom, &
!                                                            & VertInterpPool%pVIS(iIS)%ptr%grid, &
!                                                            & linear, &
!                                                            & VertInterpPool%pVIS(iIS)%ptr%recommended_update_interval, &
!                                                            & 'suppl:'+ VertInterpPool%pVIS(iIS)%ptr%chName))
!        endif
!      else
!        if(fu_cmp_verts_eq(VertInterpPool%pVIS(iIS)%ptr%vertFrom, meteo_vertical))then
!call msg('refining in non-meteo_grid but meteo_vertical:',iIS)
!          call refine_interp_vert_coefs_v2(VertInterpPool%pVIS(iIS)%ptr, & ! Structure to refine
!                                         & meteo_buf_ptr, &      ! Meteo buffer, pointers are ready
!                                         & now, &              ! Current time
!                                         & fu_horiz_interp_struct(meteo_grid, &
!                                                              & VertInterpPool%pVIS(iIS)%ptr%grid, &
!                                                              & linear))
!        else
!call msg('refining in non-meteo_grid and non-meteo vertical:',iIS)
!          call refine_interp_vert_coefs_v2(VertInterpPool%pVIS(iIS)%ptr, & ! Structure to refine
!                                         & meteo_buf_ptr, &      ! Meteo buffer, pointers are ready
!                                         & now, &              ! Current time
!                                         & fu_horiz_interp_struct(meteo_grid, &
!                                                              & VertInterpPool%pVIS(iIS)%ptr%grid, &
!                                                              & linear), &
!                                         & fu_vertical_interp_struct(meteo_vertical, &
!                                                             & VertInterpPool%pVIS(iIS)%ptr%vertFrom, &
!                                                             & VertInterpPool%pVIS(iIS)%ptr%grid, &
!                                                             & linear, &
!                                                             & VertInterpPool%pVIS(iIS)%ptr%recommended_update_interval, &
!                                                             & 'suppl:'+ VertInterpPool%pVIS(iIS)%ptr%chName))
!        endif   ! if meteo vertical
!      endif  ! if meteo grid
      if(error)return
    end do   ! inerpolation structures
  end subroutine refine_all_vert_interp_coefs


  !************************************************************************************

  subroutine refine_interp_vert_coefs_v2(interpVertStruct, & ! Structure to fill
                                       & data_buffer, &      ! Meteo or other buffer, pointers are ready
                                       & now)                ! Current time
    !
    ! Computes the interpolation defined in the interpVertStruct
    ! ATTENTION! The interpolation method, vertFrom and vertTo
    ! must be given. This routine uses this information to determine
    ! the coefs structure array.
    ! The coefficients must be in accending order, so that max and min indices
    ! are evident without min/max calls.
    ! In some cases the interpolation between the verticals require actual
    ! meteorological data. In this case the routine involves the meteo_buffer
    ! dataset in order to find the required quantities. Therefore, the pointers
    ! in that buffer must be set before calling this routine.
    ! For a simplifed crude interpolation use set_interp_vert_coefs_crude from 
    ! the module levels.
    !
    implicit none

    ! Imported parameters
    type(TVertInterpStruct), intent(inout) :: interpVertStruct
    TYPE(Tfield_buffer), intent(in) :: data_buffer
    type(silja_time), intent(in) :: now
    ! Local variables
    integer :: nLevTo, nLevFrom, nCoefs, iLev, status, nxInterp, nyInterp, iQ, nQ, &
             & ixInterp, iyInterp, interp_q_3d, interp_q_2d, ind_q_3d, ind_q_2d, &
             & ind_q_3d_forward, ind_q_2d_forward, nx_buffer, ny_buffer, iLevFrom, iLevTo, iCount
    real :: fLev, zx, zx1, zx2, zx3, zx4, weight_past, interp_surf_data, interp_surf_data_forward
    real, dimension(4) :: weights
    integer, dimension(4) :: indices
    type(THorizInterpStruct), pointer :: p_interpHoriz4Meteo
    type(TVertInterpStruct), pointer :: p_InterpVert4Meteo, p_InterpVert4Meteo_forward
    type(THorizInterpStruct), target :: HinterpMiss
    type(TVertInterpStruct), target :: VinterpMiss
    type(silja_logical) :: ifMetVertFrom
    integer, dimension(:), pointer :: q_arr_ptr
    real, dimension(max_levels) :: interp_col_data_forward, interp_col_data
    real, dimension(:,:), pointer :: weightTmp, weight_Z_tmp
    type(silja_level) :: level, levFromTmp ! Level from the vertFrom but projected to VertTo 
    logical :: if_horiz_interp, if_vert_interp, if_vert_interp_forward, ifOK
    type(silja_grid) :: grid
    type(silam_vertical) :: vertFrom_borders, vertTo_borders
    
    !
    ! Stupidity check
    !
    if(.not. (defined(interpVertStruct%vertFrom) .and. defined(interpVertStruct%vertTo)))then
    call msg("Trouble with interpVertStruct "//interpVertStruct%chName)
      call msg('vertFrom:')
      call report(interpVertStruct%vertFrom)
      call msg('vertTo')
      call report(interpVertStruct%vertTo)
      call set_error('One of verticals is undefined','refine_interp_vert_coefs')
      return
    endif
    interpVertStruct%lastCoefUpdate = time_missing   !Will be set back on 
    !
    ! Checl and get the supplementary interpolation structures if needed: in case if
    ! interpolation in the main structure is not in meteo grid and vertical, we would
    ! need these to bring in meteo data.
    !
    ! Horizontal from meteo grid to the grid of interpolation structure 
    !
    if(interpVertStruct%grid == meteo_grid)then
      if_horiz_interp = .false.
      HinterpMiss = HorizInterpStruct_missing
      p_interpHoriz4Meteo => HinterpMiss
    else
      if_horiz_interp = .true.
      p_interpHoriz4Meteo => fu_horiz_interp_struct(meteo_grid, &             ! from
                                                  & interpVertStruct%grid, &  ! to
                                                  & linear, .true., ifMakeRotation=.False.)  ! method, ifRandomise (irrelevant)
    endif
    !
    ! Vertical interp from meteo vertical to vertFrom
    ! Verticals of borders are needed to get meteo data for the interface points: 
    ! level projection understands this.
    !
    if(fu_cmp_verts_eq(interpVertStruct%vertFrom, meteo_vertical))then
      if_vert_interp = .false.
      VinterpMiss = VertInterpStruct_missing
      p_InterpVert4Meteo => VinterpMiss
    else
      if_vert_interp = .true.
      call make_vertical_of_borders(interpVertStruct%vertFrom, vertFrom_borders, .false.)
      p_InterpVert4Meteo => fu_vertical_interp_struct(meteo_vertical, &
                                                    & vertFrom_borders, & !interpVertStruct%vertFrom, &
                                                    & interpVertStruct%grid, &
                                                    & linear, &
                                                    & interpVertStruct%recommended_update_interval, &
                                                    & 'suppl:'+ interpVertStruct%chName)
    endif

    ! Find the needed meteo quantities and store their indices. Note -
    ! the 'to' and 'from' arguments seem to be reversed because we're
    ! projecting the levels of vertTo into vertFrom.
    !
    ind_q_3d = -1
    ind_q_2d = -1
    call projection_input_needs(interpVertStruct%vertTo, interpVertStruct%vertFrom, &
                              & interp_q_3d, interp_q_2d)
    if (interp_q_3d /= int_missing) then
      ind_q_3d = fu_index(data_buffer, interp_q_3d)
      if(fu_fails(ind_q_3d /= int_missing, 'Quantity not found:' + &
                       & fu_quantity_short_string(interp_q_3d), 'refine_interp_vert_coefs'))return
      grid = fu_grid(data_buffer%p4d(ind_q_3d)%past%p2d(1)%idptr)
      call grid_dimensions(grid, nx_buffer, ny_buffer)
      if (if_horiz_interp) then
        if (.not. (grid == p_interpHoriz4Meteo%gridFrom)) then
          call set_error('Buffer grid does not match horizontal interpolation structure', &
                       & 'refine_interp_vert_coefs')
          return
        end if
      else
        if (.not. fu_grids_match_closely(grid, interpVertStruct%grid)) then
          call set_error('No horizontal interpolation is given but grids do not match', &
                       & 'refine_interp_vert_coefs')
          return
        end if
      end if
    end if ! interp_q_3d /= int_missing

    if (interp_q_2d /= int_missing) then
      ind_q_2d = fu_index(data_buffer, interp_q_2d)
      if(fu_fails(ind_q_2d /= int_missing, 'Quantity not found:' + &
                          & fu_quantity_short_string(interp_q_2d),'refine_interp_vert_coefs'))return
      grid = fu_grid(data_buffer%p2d(ind_q_2d)%past%idPtr)
      
      ! call grid_dimensions also here just in case if no 3d data was needed.
      call grid_dimensions(grid, nx_buffer, ny_buffer)
      if (if_horiz_interp) then
        if (.not. (grid == p_interpHoriz4Meteo%gridFrom)) then
          call set_error('Buffer grid does not match horizontal interpolation structure', &
                       & 'refine_interp_vert_coefs')
          return
        end if
      else
        if (.not. fu_grids_match_closely(grid, interpVertStruct%grid)) then
          call set_error('No horizontal interpolation is given but grids do not match', &
                       & 'refine_interp_vert_coefs')
          return
        end if
      end if
    end if ! interp_q_2d /= int_missing
    
    ! Get all main parameters
    !
    nLevFrom = fu_nbrOfLevels(interpVertStruct%vertFrom)
    nLevTo = fu_nbrOfLevels(interpVertStruct%vertTo)
    nCoefs = fu_nCoefs(interpVertStruct)
    call grid_dimensions(interpVertStruct%grid, nxInterp, nyInterp)
    if(error)return
        
    ! The weight past here might or might not be equal to the one in met_buffer. Better compute it.
    ! 
    if(data_buffer%time_future == data_buffer%time_past)then
      weight_past = 0.0
    else
      weight_past = (data_buffer%time_future - now) / (data_buffer%time_future - data_buffer%time_past)
      if (weight_past > 1.0 .or. weight_past < 0.0) then
        call msg_warning('Interpolation coeffis for invalid time. Weight_past=' + fu_str(weight_past) + &
                       & ', past='+fu_str(data_buffer%time_past) + &
                       & ', future='+fu_str(data_buffer%time_future) + ', now='+fu_str(now), &
                     & 'refine_interp_vert_coefs')
        call set_error('Interpolation coeffis for invalid time. Weight_past=' + fu_str(weight_past), &
                     & 'refine_interp_vert_coefs')
        return
      end if
    endif


    ! Refine the existing coefficients. Crude values are used to determine the 
    ! approximate level from which we take the meteodata
    ! If the grid in the vertical structure is not the meteo_grid, have to find out
    ! the ix_meteo and iy_meteo from the horizontal interpolation structure
    !
    ! ATTENTION. 
    ! There is an inconsistency in the below logic.
    ! interpolation structures and column_from_buffer are all oriented to levels, not layers. 
    ! Exception is summation/average-type structure, which is for layers only.
    ! But fu_project_level distinguishes between these types. In particular,
    ! it assumes that for the layers column_meteo_data are for the layer interfaces, whereas
    ! for levels these are for mid-points.
    ! Solution so far:
    ! the receiving vertical in the fu_project_level (vertFrom) all turned to levels or layers,
    ! depending on the interpolation type.
    !
    select case(interpVertStruct%interp_type)

    case(average, summation)
      call msg(interpVertStruct%chName + ': refining coefficients v.2: Refine average/sum')
      !
      ! Note that the number of cells affecting the single output one may change.
      ! Solution: if too many coefficients - reallocate the arrays and restart the 
      ! whole work. There cannot be too many cases of this kind. After a few time steps
      ! things will come to sense. 
      !
      ! This reprojection requires both inverse-direction data (standard) and 
      ! forward-projection to select the layer range (see field_buffer overlap_fractions)
      !

      ind_q_3d_forward = -1
      ind_q_2d_forward = -1
      call projection_input_needs(interpVertStruct%vertFrom, interpVertStruct%vertTo, &
                                & interp_q_3d, interp_q_2d)
      if (interp_q_3d /= int_missing) then
        ind_q_3d_forward = fu_index(data_buffer, interp_q_3d)
        if(fu_fails(ind_q_3d_forward /= int_missing, 'Quantity not found:' + &
                  & fu_quantity_short_string(interp_q_3d), 'refine_interp_vert_coefs'))return
      end if
      if (interp_q_2d /= int_missing) then
        ind_q_2d_forward = fu_index(data_buffer, interp_q_2d)
        if(fu_fails(ind_q_2d_forward /= int_missing, 'Quantity not found:' + &
                  & fu_quantity_short_string(interp_q_2d),'refine_interp_vert_coefs'))return
      end if
      !
      ! Vertical interp from meteo vertical to vertTo: the "forward projection" needed for
      ! this type of interpolation
      ! Verticals of borders are needed to get meteo data for the interface points: 
      ! level projection understands this.
      !
      if(fu_cmp_verts_eq(interpVertStruct%vertTo, meteo_vertical))then
        if_vert_interp_forward = .false.
        nullify(p_InterpVert4Meteo_forward)
      else
        if_vert_interp_forward = .true.
        call make_vertical_of_borders(interpVertStruct%vertTo, vertTo_borders, .true.)
        p_InterpVert4Meteo_forward => &
                          & fu_vertical_interp_struct(meteo_vertical, &
                                                    & vertTo_borders, & ! need meteo at interface levels
                                                    & interpVertStruct%grid, &
                                                    & linear, &
                                                    & interpVertStruct%recommended_update_interval, &
                                                    & 'suppl_forward:'+ interpVertStruct%chName)
      endif

      weightTmp => fu_work_array_2d()
      weight_Z_Tmp => fu_work_array_2d()
      if(error)return
      !
      ! Start the main cycle
      !
      ifOK = .false.

      ifSufficientSize: do while(.not. ifOK)
        ifOK = .true.
        do iyInterp = 1, nyInterp
          do ixInterp = 1, nxInterp
            call column_from_buffer(data_buffer, ind_q_3d, ixInterp, iyInterp, nx_buffer, &
                                  & interp_col_data, &
                                  & p_interpHoriz4Meteo, p_InterpVert4Meteo, &
                                  & if_horiz_interp, if_vert_interp, weight_past)
            call surf_from_buffer(data_buffer, ind_q_2d, ixInterp, iyInterp, nx_buffer, &
                                & interp_surf_data, &
                                & p_interpHoriz4Meteo, if_horiz_interp, weight_past)
            call column_from_buffer(data_buffer, ind_q_3d_forward, ixInterp, iyInterp, nx_buffer, &
                                  & interp_col_data_forward, &
                                  & p_interpHoriz4Meteo, p_InterpVert4Meteo_forward, &
                                  & if_horiz_interp, if_vert_interp_forward, weight_past)
            call surf_from_buffer(data_buffer, ind_q_2d, ixInterp, iyInterp, nx_buffer, &
                                & interp_surf_data_forward, &
                                & p_interpHoriz4Meteo, if_horiz_interp, weight_past)
            if(error)return
            !
            ! Fill-ing the temporary arrays. They are (nLevsFrom*nLevsTo)
            ! ATTENTION.
            ! Inconsistency.
            ! This function deals exclusively with the layers - but the meteodata are
            ! for midpoints of these layers. Have to make and send there also the
            ! corresponding level-based verticals.
            !
            call overlap_fraction_lyr_in_vert(interpVertStruct%vertFrom, &
                                            & interpVertStruct%vertTo, &
                                            & weightTmp, weight_Z_tmp, &
                                            & interp_col_data_forward, interp_surf_data_forward, &
                                            & interp_col_data, interp_surf_data)
            if(error)return
            !
            ! Copy the non-zero weights checking that the number of coefs is not too large
            !
            do iLevTo = 1, nLevTo
              interpVertStruct%indLev(1:interpVertStruct%nCoefs,ixInterp,iyInterp,iLevTo) = 0
              interpVertStruct%weight(1:interpVertStruct%nCoefs,ixInterp,iyInterp,iLevTo) = 0.0
              interpVertStruct%weight_Z(1:interpVertStruct%nCoefs,ixInterp,iyInterp,iLevTo) = 0.0
              iCount = 1
              do iLevFrom = 1, nLevFrom
                if(weightTmp(iLevFrom,iLevTo) .eps. 0.0)cycle
                if(iCount > interpVertStruct%nCoefs)then
                  call msg_warning('Have to increase the number of coefficients in:' + &
                                 & interpVertStruct%chName,'refine_interp_vert_coefs')
                  ifOK = .false.
                  interpVertStruct%nCoefs = interpVertStruct%nCoefs + 1
                  deallocate(interpVertStruct%indLev, interpVertStruct%weight, &
                           & interpVertStruct%weight_Z)
                  allocate( &
                      & interpVertStruct%indLev(interpVertStruct%nCoefs,nxInterp,nyInterp,nLevTo), &
                      & interpVertStruct%weight(interpVertStruct%nCoefs,nxInterp,nyInterp,nLevTo), &
                      & interpVertStruct%weight_Z(interpVertStruct%nCoefs,nxInterp,nyInterp,nLevTo), &
                      & stat = iCount)
                  if(fu_fails(iCount==0,'Reallocation failed','refine_interp_vert_coefs'))return
                  cycle ifSufficientSize
                endif
                interpVertStruct%indLev(iCount,ixInterp,iyInterp,iLevTo) = iLevFrom
                interpVertStruct%weight(iCount,ixInterp,iyInterp,iLevTo) = weightTmp(iLevFrom,iLevTo)
                interpVertStruct%weight_Z(iCount,ixInterp,iyInterp,iLevTo) = &
                                                                      & weight_Z_Tmp(iLevFrom,iLevTo)
                iCount = iCount + 1
              enddo
            end do  ! iLevTo
          end do   ! ix
        end do  ! iy
      end do ifSufficientSize ! whlie ifOK


    case(nearest_point)
      !
      ! Cannot handle layers...
!      call make_vertical_of_levels(interpVertStruct%vertFrom, vertOfLevels)
!      if(error)return
    call msg(interpVertStruct%chName + ': refining coefficients v.2 :nearest')

      do iyInterp = 1, nyInterp
        do ixInterp = 1, nxInterp
          call column_from_buffer(data_buffer, ind_q_3d, ixInterp, iyInterp, nx_buffer, &
                                & interp_col_data, &
                                & p_interpHoriz4Meteo, p_InterpVert4Meteo, &
                                & if_horiz_interp, if_vert_interp, weight_past)
          call surf_from_buffer(data_buffer, ind_q_2d, ixInterp, iyInterp, nx_buffer, &
                              & interp_surf_data, &
                              & p_interpHoriz4Meteo, if_horiz_interp, weight_past)
          do iLev=1,nLevTo
            ! Project each levelTo to the vertFrom and then check the index value
            !
            level = fu_level(interpVertStruct%vertTo, ilev)
            fLev = fu_project_level(level, interpVertStruct%vertFrom, &
                                  & interp_col_data, interp_surf_data)
            if(error)return

!            flev = max(1.0, min(real(nLevFrom), flev))
!            indices(1) = int(fLev)
!            if(fLev - indices(1) > indices(1)+1 - fLev) indices(1) = indices(1) + 1 ! Which one is closer?
            interpVertStruct%indLev(1,ixInterp,iyInterp,iLev) = min(nLevFrom, max(1,nint(fLev)))
          end do
        end do
      end do

    case(linear)

!      call make_vertical_of_levels(interpVertStruct%vertFrom, vertOfLevels)
!      if(error)return
    call msg(interpVertStruct%chName + ': refining coefficients v.2: linear')
    
      !$OMP PARALLEL DEFAULT(NONE) & 
      !$OMP & SHARED(interpVertStruct,  nLevTo, ifMetVertFrom, data_buffer, nx_meteo, now, &
      !$OMP        & p_interpVert4Meteo, nxInterp, nyInterp, ny_meteo, nq, p_interpHoriz4Meteo, &
      !$OMP        & error, nLevFrom, weight_past, ind_q_3d, ind_q_2d, nx_buffer,if_horiz_interp, &
      !$OMP        & if_vert_interp, vertFrom_borders) &
      !$OMP & PRIVATE(iyInterp, ixInterp, ilev, level, fLev, interp_col_data, interp_surf_data, &
      !$OMP         & zx, weights, indices)

      !$OMP DO COLLAPSE(2)
      do iyInterp = 1, nyInterp
        do ixInterp = 1, nxInterp
          !if (.not. defined(interpVertStruct%vertFrom)) then
          !  call msg('undef')
          !  stop
          !end if
          if (ind_q_3d > 0) then 
            call column_from_buffer(data_buffer, ind_q_3d, ixInterp, iyInterp, nx_buffer, &
                                  & interp_col_data, &
                                  & p_interpHoriz4Meteo, p_InterpVert4Meteo, &
                                  & if_horiz_interp, if_vert_interp, weight_past)
          end if
          if (ind_q_2d > 0) then
            call surf_from_buffer(data_buffer, ind_q_2d, ixInterp, iyInterp, nx_buffer, &
                                & interp_surf_data, &
                                & p_interpHoriz4Meteo, if_horiz_interp, weight_past)
          end if

          do iLev=1,nLevTo
            level = fu_level(interpVertStruct%vertTo, ilev, .true.)
            flev = fu_project_level(level, interpVertStruct%vertFrom, &
                                  & interp_col_data, interp_surf_data)
            flev = max(1.0, min(real(nLevFrom), flev))

            !flev = 1.0
            !flev = min(max(1.0, flev), real(nLevFrom))
            zx = fLev - REAL(INT(fLev)) ! Exceedance of fLev over lower border
            
            indices(1) = int(flev) !min(nLevFrom,max(1,int(fLev)))
            weights(1) = 1. - zx
            indices(2) = min(nLevFrom,int(flev)+1)
            weights(2) = zx

            interpVertStruct%indLev(1:2,ixInterp,iyInterp,iLev) = indices(1:2)
            interpVertStruct%weight(1:2,ixInterp,iyInterp,iLev) = weights(1:2)
          end do ! ilev
          if (error) cycle
        end do ! ix
      end do ! iy
      !$OMP END DO

      !$OMP END PARALLEL

      case(cubic)
!        call make_vertical_of_levels(interpVertStruct%vertFrom, vertOfLevels)
!        if(error)return
    call msg(interpVertStruct%chName + ': refining coefficients v.2 :cubic')

        do iyInterp = 1, nyInterp
          do ixInterp = 1, nxInterp
            call column_from_buffer(data_buffer, ind_q_3d, ixInterp, iyInterp, nx_buffer, &
                                  & interp_col_data, &
                                  & p_interpHoriz4Meteo, p_InterpVert4Meteo, &
                                  & if_horiz_interp, if_vert_interp, weight_past)
            call surf_from_buffer(data_buffer, ind_q_2d, ixInterp, iyInterp, nx_buffer, &
                                & interp_surf_data, &
                                & p_interpHoriz4Meteo, if_horiz_interp, weight_past)

            do iLev=1,nLevTo
              level = fu_level(interpVertStruct%vertTo, ilev, .true.)
              fLev = fu_project_level(level, interpVertStruct%vertFrom, &
                                    & interp_col_data, interp_surf_data)
              flev = max(1.0, min(real(nLevFrom), flev))
              if(error)return

              zx = fLev - REAL(INT(fLev)) ! Exceedance of fLev over lower border

              indices(1) = min(nLevFrom,max(1,int(fLev)-1))
              weights(1) = ((-0.5*zx+1.0)*zx-0.5)*zx 
              indices(2) = min(nLevFrom,max(1,int(fLev)))
              weights(2) = (( 1.5*zx-2.5)*zx    )*zx+1 
              indices(3) = min(nLevFrom,max(1,int(fLev)+1))
              weights(3) = ((-1.5*zx+2.0)*zx+0.5)*zx 
              indices(4) = min(nLevFrom,max(1,int(fLev)+2))
              weights(4) = (( 0.5*zx-0.5)*zx    )*zx 

              interpVertStruct%indLev(1:4,ixInterp,iyInterp,iLev) = indices(1:4)
              interpVertStruct%weight(1:4,ixInterp,iyInterp,iLev) = weights(1:4)
            end do
          end do
        end do

      case default
        call msg('Unknown interpolation method:',interpVertStruct%interp_type)
        call set_error('Unknown interpolation method','refine_interp_vert_coefs')
        return
    end select
 
    interpVertStruct%lastCoefUpdate = now
 
  end subroutine refine_interp_vert_coefs_v2


  !*************************************************************************************
    
  subroutine column_from_buffer(data_buffer, ind_q, ixInterp, iyInterp, nx_buffer, col_data, &
                              & p_interpHorizStruct, p_InterpVert4Meteo, &
                              & if_horiz_interp, if_vert_interp, weight_past)
    implicit none
    type(Tfield_buffer), intent(in) :: data_buffer
    integer, intent(in) :: ind_q, ixInterp, iyInterp, nx_buffer
    real, dimension(:), intent(out) :: col_data
    type(THorizInterpStruct), intent(in) :: p_interpHorizStruct
    type(TVertInterpStruct), intent(in) :: p_InterpVert4Meteo
    logical, intent(in) :: if_horiz_interp, if_vert_interp
    real, intent(in) :: weight_past

    integer :: ilev, ix, iy, i1d
    TYPE(field_3d_data_ptr), pointer :: past, future

    if (.not. if_horiz_interp .and. .not. if_vert_interp) then
      ! Met data is in the interpolation grid and vertical_from -
      ! eg. interpolation from meteo_buffer to meteo_grid. 
      !
      past => data_buffer%p4d(ind_q)%past
      future => data_buffer%p4d(ind_q)%future
      i1d = (iyInterp - 1) * nx_buffer + ixInterp ! nx_buffer is also nx of the interpolation grid.
      do ilev = 1, data_buffer%nbr_of_levels             ! OK for no vertical interpolation
        col_data(ilev) = past%p2d(ilev)%ptr(i1d) * weight_past + &
                       & future%p2d(ilev)%ptr(i1d) * (1.0 - weight_past)
      end do

    elseif (.not. if_vert_interp) then
      !!
      !! Need to interpolate horizontally but not vertically
      !! 
      !ix = int(sum(real(p_interpHorizStruct%indX(1:p_interpHorizStruct%nCoefs,ixInterp,iyInterp)) * &
      !             & p_interpHorizStruct%weight(:,ixInterp,iyInterp))+0.5)
      !iy = int(sum(real(p_interpHorizStruct%indY(1:p_interpHorizStruct%nCoefs,ixInterp,iyInterp)) * &
      !           & p_interpHorizStruct%weight(:,ixInterp,iyInterp))+0.5)
      !past => data_buffer%p4d(ind_q)%past
      !future => data_buffer%p4d(ind_q)%future
      !i1d = (iy-1) * nx_buffer + ix
      !do ilev = 1, data_buffer%nbr_of_levels             ! OK for no vertical interpolation
      !  col_data(ilev) = past%p2d(ilev)%ptr(i1d) * weight_past + &
      !                 & future%p2d(ilev)%ptr(i1d) * (1.0 - weight_past)
      !end do
      do ilev = 1, data_buffer%nbr_of_levels             ! OK for no vertical interpolation
        col_data(ilev) = fu_get_value(data_buffer%p4d(ind_q), nx_buffer, &
                                    & ixInterp, iyInterp, ilev, weight_past, &
                                    & p_interpHorizStruct, p_InterpVert4Meteo, &
                                    & if_horiz_interp, .false.)
      end do
    else
      !
      ! Interpolation in vertical and (may be) horizontal, eg. from boundary
      ! grid & vertical to dispserion grid (/= meteo_grid) and vertical.
      !
      do ilev = 1, fu_nbrOfLevels(p_InterpVert4Meteo%vertTo)
        col_data(ilev) = fu_get_value(data_buffer%p4d(ind_q), nx_buffer, &
                                    & ixInterp, iyInterp, ilev, weight_past, &
                                    & p_interpHorizStruct, p_InterpVert4Meteo, &
                                    & if_horiz_interp, .true.)
      end do
    end if
      
  end subroutine column_from_buffer


  !***************************************************************

  subroutine surf_from_buffer(data_buffer, ind_q, ixInterp, iyInterp, nx_buffer, surf_data, &
                            & p_interpHorizStruct, if_horiz_interp, weight_past)
    implicit none
    TYPE(Tfield_buffer), intent(in) :: data_buffer
    integer, intent(in) :: ind_q, ixInterp, iyInterp, nx_buffer
    real, intent(out) :: surf_data
    type(THorizInterpStruct), pointer :: p_interpHorizStruct
    logical, intent(in) :: if_horiz_interp
    real, intent(in) :: weight_past

    real :: val_past, val_future
    integer :: i1d

    if (.not. if_horiz_interp) then
      i1d = (iyInterp - 1) * nx_buffer + ixInterp
      val_past = data_buffer%p2d(ind_q)%past%ptr(i1d)
      val_future = data_buffer%p2d(ind_q)%future%ptr(i1d)
      surf_data = val_past*weight_past + val_future*(1.0-weight_past)
    else
      surf_data = fu_get_value(data_buffer%p2d(ind_q), nx_buffer, ixInterp, iyInterp, &
                             & weight_past, p_interpHorizStruct, &
                             & ifHorizInterp=.true., ifForceWeightPast_=.true.)
    end if
  end subroutine surf_from_buffer


  !********************************************************************************

  integer function fu_ind_vert_from_bfr_interp2d(ixTo, iyTo, iLevTo, &
                                               & pVertInterpStruct, ifVertInterp) result(iLev)
    !
    ! The very basic function: returns an index in the From grid, which corresponds
    ! to the given ixTo, iyTo. It can later be used for direct addressing the values from the 
    ! gridFrom. However, it is advisable to use then fu_get_value because that function
    ! makes horizontal interpolation rather than picks a single value
    !
    implicit none

    ! Imported parameters
    type(TVertInterpStruct), pointer :: pVertInterpStruct
    integer, intent(in) :: ixTo, iyTo, iLevTo
    logical, intent(in) :: ifVertInterp

    ! Local variables
    real :: fLevFrom
    integer :: iCoef
    !
    ! ATTENTION. ABSOLUTELY NO CHECKING TO KEEP THE MAX SPEED
    !
    if(ifVertInterp)then
      !
      ! Full-blown horizontal interpolation
      !
      fLevFrom = 0.
      do iCoef = 1, pVertInterpStruct%nCoefs
        fLevFrom = fLevFrom + pVertInterpStruct%weight(iCoef,ixTo,iyTo,iLevTo) * &
                            & pVertInterpStruct%indLev(iCoef,ixTo,iyTo,iLevTo)
      end do
      iLev = int(fLevFrom+0.5)
    else
      iLev = iLevTo
    endif  ! if vertical interpolation needed

  end function fu_ind_vert_from_bfr_interp2d



  !******************************************************************
  !
  ! Encapsulation of the vertical interpolation structure
  !
  !******************************************************************

  function fu_vertFrom_from_interp_struct(interpStructVert)result(vertFrom)
    implicit none
    type(silam_vertical) :: vertFrom
    type(TVertInterpStruct), intent(in) :: interpStructVert
    vertFrom = interpStructVert%vertFrom
  end function fu_vertFrom_from_interp_struct

  function fu_vertTo_from_interp_struct(interpStructVert)result(vertTo)
    implicit none
    type(silam_vertical) :: vertTo
    type(TVertInterpStruct), intent(in) :: interpStructVert
    vertTo = interpStructVert%vertTo
  end function fu_vertTo_from_interp_struct

  integer function fu_nCoefs_interp_str_vert(interpStructVert)
    implicit none
    type(TVertInterpStruct), intent(in) :: interpStructVert
    fu_nCoefs_interp_str_vert = interpStructVert%nCoefs
  end function fu_nCoefs_interp_str_vert

  integer function fu_interpType_interp_str_vert(interpStructVert)
    implicit none
    type(TVertInterpStruct), intent(in) :: interpStructVert
    fu_interpType_interp_str_vert = interpStructVert%interp_type
  end function fu_interpType_interp_str_vert

  subroutine get_coefs_interp_str_vert(interpStructVert, coefsInterp)
    implicit none
    type(TVertInterpCells), intent(out) :: coefsInterp ! Just the size of gridTo
    type(TVertInterpStruct), intent(in) :: interpStructVert
    coefsInterp%indLev => interpStructVert%indLev
    coefsInterp%weight => interpStructVert%weight
    coefsInterp%weight_Z => interpStructVert%weight_Z
  end subroutine get_coefs_interp_str_vert

  subroutine get_coefs_interp_str_cell_vert(interpStructVert, ix, iy, iLev, coefsInterp)
    implicit none
    type(TVertInterpOneCell), intent(out) :: coefsInterp
    integer, intent(in) :: ix,iy,iLev
    type(TVertInterpStruct), intent(in), target :: interpStructVert
    coefsInterp%indLev => interpStructVert%indLev(:,ix,iy,iLev)
    coefsInterp%weight => interpStructVert%weight(:,ix,iy,iLev)
    coefsInterp%weight_Z => interpStructVert%weight_Z(:,ix,iy,iLev)
  end subroutine get_coefs_interp_str_cell_vert

  logical function fu_ifMeteoGrd(interpStructVert)
    implicit none
    type(TVertInterpStruct), intent(in) :: interpStructVert
    fu_ifMeteoGrd = interpStructVert%ifMeteoGrid
  end function fu_ifMeteoGrd

  function fu_grid_of_vertInterpStruct(interpStructVert)result(grid)
    implicit none
    type(TVertInterpStruct), intent(in) :: interpStructVert
    type(silja_grid) :: grid
    grid = interpStructVert%grid
  end function fu_grid_of_vertInterpStruct

  !*********************************************************************

  subroutine set_missing_vert_interp_struct(InterpStruct)
    implicit none
    type(TVertInterpStruct), intent(out) :: InterpStruct

    call set_missing(InterpStruct%vertFrom, .true.)
    call set_missing(InterpStruct%vertTo, .true.)
    InterpStruct%grid = grid_missing
    InterpStruct%ifMeteoGrid = .false.
    InterpStruct%interp_type = int_missing
    InterpStruct%nCoefs = int_missing
    nullify(InterpStruct%indLev)
    nullify(InterpStruct%weight)
    nullify(InterpStruct%weight_Z)
  end subroutine set_missing_vert_interp_struct

  !********************************************************************

  subroutine get_vertical_range_vert_interp(interpStructVert, &
                                          & iLevStartFrom, iLevEndFrom, iLevStartTo, iLevEndTo)
    implicit none
    type(TVertInterpStruct), intent(in) :: interpStructVert
    integer, intent(out) :: iLevStartFrom, iLevEndFrom, iLevStartTo, iLevEndTo
    iLevStartFrom = interpStructVert%iLevStartFrom
    iLevEndFrom = interpStructVert%iLevEndFrom
    iLevStartTo = interpStructVert%iLevStartTo
    iLevEndTo = interpStructVert%iLevEndTo
  end subroutine get_vertical_range_vert_interp


  !*******************************************************************************
  !*******************************************************************************
  !
  !
  !   Reporting routines
  !
  !
  !*******************************************************************************
  !*******************************************************************************

  !*******************************************************************************

  subroutine report_buffer(buf, iSwitch)
    !
    ! Depending on the switch value, prints the buffer overview, statistics for 
    ! all the fields in the buffer or dumps the whole content to GrADS-type file
    ! Corresponding values for the switch - 1,2 and 3.
    !
    implicit none

    ! Imported parameters
    type(Tfield_buffer), intent(in) :: buf
    integer, intent(in) :: iSwitch

    ! Local vars
    integer :: iQ, iLev, iCount, i, fs
    type(silam_sp) :: strTmp

    strTmp%sp => fu_work_string()
    if(error)return

    call msg('================== FIELD BUFFER REPORT =================')
    call msg('Number of levels:', buf%nbr_of_levels)
    call msg('Buffer quantities:')

    do iQ=1, size(buf%buffer_quantities)
      !
      ! Scan all quantities one-by-one reporting the field statistics
      ! 
      if(fu_dimension_of_meteofield(buf, iQ) == 2)then
        call msg('2D quantity:',iQ)
        if(buf%p2d(iQ)%present%ifReady)then
          call report(buf%p2d(iQ)%present%idPtr)
        elseif(buf%p2d(iQ)%past%ifReady)then
          call report(buf%p2d(iQ)%past%idPtr)
        elseif(buf%p2d(iQ)%future%ifReady)then
          call report(buf%p2d(iQ)%future%idPtr)
        else
          call msg_warning('None of the past/pres/future IDs is ready','report_buffer')
        endif
      else
        call msg('4D quantity:',iQ)
        if(buf%p4d(iQ)%present%p2d(1)%ifReady)then
          call report(buf%p4d(iQ)%present%p2d(1)%idPtr)
        elseif(buf%p4d(iQ)%past%p2d(1)%ifReady)then
          call report(buf%p4d(iQ)%past%p2d(1)%idPtr)
        elseif(buf%p4d(iQ)%future%p2d(1)%ifReady)then
          call report(buf%p4d(iQ)%future%p2d(1)%idPtr)
        else
          call msg_warning('None of the past/pres/future IDs is ready','report_buffer')
        endif
      endif
      if(buf%ifPointerSet(iQ))then
        call msg('Pointers are set')

        if(iSwitch > 1)then
          !
          ! Statistics for all the fields
          !
          if(fu_dimension_of_meteofield(buf,iQ) == 2)then
            !
            ! 2D quantity past. Fields are printed only ifReady is set
            !
            if(buf%p2d(iQ)%past%ifReady)then
              fs = fu_number_of_gridpoints(fu_grid(buf%p2d(iQ)%past%idPtr))
              iCount=0
              do i=1,fs
                if(buf%p2d(iQ)%past%ptr(i) .eps. real_missing) iCount=iCount+1
              end do
              write(unit=strTmp%sp,fmt=*) &
                  & 'PAST Minimum=',MINVAL(buf%p2d(iQ)%past%ptr(1:fs)), &
                  & 'Average=',(SUM(buf%p2d(iQ)%past%ptr(1:fs))/real(fs)),&
                  & 'Maximum=',MAXVAL(buf%p2d(iQ)%past%ptr(1:fs)), &
                  & ', nMissing=',iCount
              call msg(strTmp%sp)
            else
              call msg('PAST field is not ready')
            endif
            !
            ! 2D PRESENT
            !
            if(buf%p2d(iQ)%present%ifReady)then
              fs = fu_number_of_gridpoints(fu_grid(buf%p2d(iQ)%present%idPtr))
              iCount=0
              do i=1,fs
                if(buf%p2d(iQ)%present%ptr(i) .eps. real_missing) iCount=iCount+1
              end do
              write(unit=strTmp%sp,fmt=*) &
                  & 'PRESENT Min=',MINVAL(buf%p2d(iQ)%present%ptr(1:fs)), &
                  & 'Mean=',(SUM(buf%p2d(iQ)%present%ptr(1:fs))/real(fs)),&
                  & ', Max=',MAXVAL(buf%p2d(iQ)%present%ptr(1:fs)), &
                  & ', nMissing=',iCount
              call msg(strTmp%sp)
            else
              call msg('PRESENT field is not ready')
            endif
            !
            ! 2D FUTURE
            !
            if(buf%p2d(iQ)%future%ifReady)then
              fs = fu_number_of_gridpoints(fu_grid(buf%p2d(iQ)%future%idPtr))
              iCount=0
              do i=1,fs
                if(buf%p2d(iQ)%future%ptr(i) .eps. real_missing) iCount=iCount+1
              end do
              write(unit=strTmp%sp,fmt=*) &
                  & 'FUTURE Minimum=',MINVAL(buf%p2d(iQ)%future%ptr(1:fs)), &
                  & 'Average=',(SUM(buf%p2d(iQ)%future%ptr(1:fs))/real(fs)),&
                  & 'Maximum=',MAXVAL(buf%p2d(iQ)%future%ptr(1:fs)), &
                  & ', nMissing=',iCount
              call msg(strTmp%sp)
            else
              call msg('FUTURE field is not ready')
            endif
          else ! if2d
            !
            ! 4D quantity. Fields are printed only ifReady is set
            !
            do iLev=1,buf%nbr_of_levels
              !
              ! 4D PAST
              !
              if(buf%p4d(iQ)%past%p2d(iLev)%ifReady)then
                fs = fu_number_of_gridpoints(fu_grid(buf%p4d(iQ)%past%p2d(1)%idPtr))
                iCount=0
                do i=1,fs
                  if(buf%p4d(iQ)%past%p2d(1)%ptr(i) .eps. real_missing) iCount=iCount+1
                end do
                write(unit=strTmp%sp,fmt=*) &
                 & 'PAST Minimum=',MINVAL(buf%p4d(iQ)%past%p2d(iLev)%ptr(1:fs)), &
                 & 'Average=',(SUM(buf%p4d(iQ)%past%p2d(iLev)%ptr(1:fs))/real(fs)),&
                 & 'Maximum=',MAXVAL(buf%p4d(iQ)%past%p2d(iLev)%ptr(1:fs)), &
                  & ', nMissing=',iCount
                call msg(strTmp%sp)
              else
                call msg('PAST field is not ready')
              endif
              !
              ! 4D PRESENT
              !
              if(buf%p4d(iQ)%present%p2d(iLev)%ifReady)then
                fs = fu_number_of_gridpoints(fu_grid(buf%p4d(iQ)%present%p2d(1)%idPtr))
                iCount=0
                do i=1,fs
                  if(buf%p4d(iQ)%present%p2d(1)%ptr(i) .eps. real_missing) iCount=iCount+1
                end do
                write(unit=strTmp%sp,fmt=*) &
                  & 'PRESENT Minimum=',MINVAL(buf%p4d(iQ)%present%p2d(iLev)%ptr(1:fs)), &
                  & 'Average=',(SUM(buf%p4d(iQ)%present%p2d(iLev)%ptr(1:fs))/real(fs)),&
                  & 'Maximum=',MAXVAL(buf%p4d(iQ)%present%p2d(iLev)%ptr(1:fs)), &
                  & ', nMissing=',iCount
                call msg(strTmp%sp)
              else
                call msg('PRESENT field is not ready')
              endif
              !
              ! 4D FUTURE
              !
              if(buf%p4d(iQ)%future%p2d(iLev)%ifReady)then
                fs = fu_number_of_gridpoints(fu_grid(buf%p4d(iQ)%future%p2d(1)%idPtr))
                iCount=0
                do i=1,fs
                  if(buf%p4d(iQ)%future%p2d(1)%ptr(i) .eps. real_missing) iCount=iCount+1
                end do
                write(unit=strTmp%sp,fmt=*) &
                  & 'FUTURE Minimum=',MINVAL(buf%p4d(iQ)%future%p2d(iLev)%ptr(1:fs)), &
                  & 'Average=',(SUM(buf%p4d(iQ)%future%p2d(iLev)%ptr(1:fs))/real(fs)),&
                  & 'Maximum=',MAXVAL(buf%p4d(iQ)%future%p2d(iLev)%ptr(1:fs)), &
                  & ', nMissing=',iCount
                call msg(strTmp%sp)
              else
                call msg('FUTURE field is not ready')
              endif
            end do ! iLev
          endif

          if(iSwitch > 2)then
            call set_error('Do not support detail level higher than 2','report_buffer')
            call unset_error('report_buffer')
          endif

        endif ! iSwitch > 1

      else ! ifPointerSet
        call msg('Pointers are not set')
      endif

    end do ! iQ - quantities

    call msg('============== END FIELD BUFFER REPORT =================')

    call free_work_array(strTmp%sp)

  end subroutine report_buffer

  !***************************************************************************************


  subroutine test_project_level()
    implicit none
    
    integer, parameter :: nn = 41
    real, dimension(nn), parameter :: &
         & a = (/0.0, 2006.06, 3996.76, &
              & 5923.68, 7748.44, 9441.12, 10978.87, 12344.79, 13526.95, 14517.65, &
              & 15312.73, 15911.13, 16314.44, 16526.66, 16553.91, 16404.29, &
              & 16087.75, 15615.96, 15002.23, 14261.41, 13409.81, 12465.03, &
              & 11445.84, 10372.0, 9263.99, 8142.78, 7029.45, 5944.81, 4909.01, &
              & 3940.95, 3057.76, 2274.19, 1601.89, 1048.75, 618.02, 307.62, &
              & 109.2, 7.27, 0.0, 0.0, 0.0/), &
         
         & b = (/0.0, 0.0, 0.0, 0.0, 0.0, 0.01, 0.01, 0.02, 0.04, &
              0.05, 0.07, 0.09, 0.12, 0.15, 0.18, 0.22, 0.25, 0.29, 0.33, &
              0.38, 0.42, 0.46, 0.51, 0.55, 0.6, 0.64, 0.68, 0.72, 0.76, &
              0.79, 0.83, 0.86, 0.88, 0.91, 0.93, 0.95, 0.96, 0.97, 0.98, &
              0.99, 1.0/), &
              
         & sigma = (/0.0, 0.0197982940075, 0.0394449964396, &
                  & 0.0584622385405, 0.0764712387565, 0.103176709331, &
                  & 0.118353138057, 0.141833734724, 0.173500759262, &
                  & 0.193278218497, 0.221125056378, 0.247030811507, &
                  & 0.281011175981, 0.313105626772, 0.343374563649, &
                  & 0.381897927482, 0.408773917241, 0.444117706993, &
                  & 0.478060656366, 0.52074932362, 0.552344676114, &
                  & 0.583020412526, 0.622961778553, 0.652363790438, &
                  & 0.691428570284, 0.720363076119, 0.74937535159, &
                  & 0.778670775649, 0.808448213546, 0.828894193977, &
                  & 0.860177777078, 0.882444534186, 0.895809441985, &
                  & 0.920350368803, 0.936099389681, 0.953035976592, &
                  & 0.961077721357, 0.970071749398, 0.98, 0.99, 1.0/), &

         & height = (/0.0, 84.6887171594, 170.072974027, 255.545404244, &
                   & 333.588790243, 403.871373849, 553.477455535, 694.573138588, &
                   & 918.38286538, 1042.36468223, 1252.34487713, 1554.91984704, &
                   & 1757.70821799, 2060.60802189, 2367.90049014, 2681.96431296, &
                   & 3005.55876162, 3460.31113724, 3817.35648327, 4324.91930179, &
                   & 4734.21990025, 5175.50764187, 5807.54162922, 6343.58410412, &
                   & 6938.24363648, 7418.97056005, 8158.22902535, 8787.79381827, &
                   & 9511.66176116, 10355.0912502, 11063.805252, 11917.3727792, &
                   & 12601.9425516, 13879.9469436, 15027.6700619, 15897.9287197, &
                   & 17797.4242392, 19500.3626345, 22004.8436445, 26461.4559563, 100000.0/), &

         & p = (/101324.891894, 100311.642975, 99298.3940561, 98292.4151372, &
              & 97381.0962182, 96566.2672993, 94850.1694614, 93254.4016235, &
              & 90767.7948667, 89413.5970288, 87157.420272, 83987.6145963, &
              & 81915.9278394, 78898.7321637, 75930.3764879, 72990.7108122, &
              & 70058.9251364, 66100.6905417, 63121.5348659, 59074.4802712, &
              & 55966.2645955, 52764.8689197, 48439.444325, 45000.1786493, &
              & 41418.9729735, 38695.7662167, 34792.3905409, 31725.3937841, &
              & 28473.4270273, 25030.3702705, 22405.4724326, 19583.8945947, &
              & 17579.9456758, 14371.2878379, 11992.1189189, 10454.3689189, 7748.44, &
              & 5923.68, 3996.76, 2006.06, 0.0/)

    real, dimension(12) :: hgts = (/0.0, 50.0, 150.0, 300.0, 600.0, 1200.0, &
                                 & 2000.0, 3000.0, 4000.0, 5000.0, 6500.0, 8000.0/)

    real, dimension(12) :: p_small = (/101324.891894, 100725.675485, &
                                    & 99535.850101, 97772.4728, 94321.5833835, 87715.4879072, &
                                    & 79495.1305577, 70108.4696855, 61640.1693196, 54019.8542264, &
                                    & 44034.7988733, 35599.7731748/)

    real, dimension(nn) :: dummy_array
    real :: dummy
    real, parameter :: p_surf = 101324.891894
    integer :: iz, q2d, q3d
    real, dimension(:), pointer :: metdata

    type(silam_vertical) :: vert_hybrid, vert_sigma, vert_press, vert_hgt_small, vert_hgt, vert_press_small
    type(silja_level), dimension(nn) :: levels
    type(silja_level) :: lev
    real :: float
    do iz = 1, nn
      levels(iz) = fu_set_level(hybrid, fval1=a(iz), fval2=b(iz))
    end do
    call set_vertical(levels, vert_hybrid)
    call arrange_levels_in_vertical(vert_hybrid)

    do iz = 1, nn
      levels(iz) = fu_set_level(sigma_level, fval1=sigma(iz))
    end do
    call set_vertical(levels, vert_sigma)
    call arrange_levels_in_vertical(vert_sigma)
    
    do iz = 1, nn
      levels(iz) = fu_set_level(constant_pressure, fval1=p(iz))
    end do
    call set_vertical(levels, vert_press)
    call arrange_levels_in_vertical(vert_press)
    
    do iz = 1, size(hgts) 
      levels(iz) = fu_set_level(constant_height, fval1=hgts(iz))
    end do
    call set_vertical(levels(1:12), vert_hgt_small)
    call arrange_levels_in_vertical(vert_hgt_small)

    do iz = 1, nn
      levels(iz) = fu_set_level(constant_height, fval1=height(iz))
    end do
    call set_vertical(levels, vert_hgt)
    call arrange_levels_in_vertical(vert_hgt)

    do iz = 1, 12
      levels(iz) = fu_set_level(constant_pressure, fval1=p_small(iz))
    end do
    call set_vertical(levels(1:12), vert_press_small)
    print *, 'to arr'
    call arrange_levels_in_vertical(vert_press_small)
    print *, 'from arr'
    !print *, 'Project hybrid(1) to hybrid', &
    !     & fu_project_level(fu_level(vert_hybrid, 1), vert_hybrid, dummy_array, p_surf)
    float = fu_project_level(fu_level(vert_hybrid, 1), vert_hybrid, dummy_array, p_surf)
    print *, 'foo', fu_project_level(fu_level(vert_hybrid, 1), vert_hybrid, dummy_array, p_surf)
    print *, 'Project hybrid(3) to hybrid', &
         & fu_project_level(fu_level(vert_hybrid, 3), vert_hybrid, dummy_array, p_surf)

    
    print *, 'Project height(2) to hybrid', &
         & fu_project_level(fu_level(vert_hgt_small, 2), vert_hybrid, height, dummy)
    ! 1.59 seems to be correct
    
    print *, 'Project height(8) to hybrid', &
         & fu_project_level(fu_level(vert_hgt_small, 8), vert_hybrid, height, dummy)
    ! ~17 is ok
    
    print *, 'Project height(8) to sigma', &
         & fu_project_level(fu_level(vert_hgt_small, 8), vert_sigma, height, dummy)
    
    print *, 'Project 850 hPa to hybrid', &
         & fu_project_level(fu_set_level(constant_pressure, fval1=85000.0), &
                          & vert_hybrid, dummy_array, p_surf)
    ! should be about 11.68
    
    print *, 'Project 850 hPa to big height', &
         & fu_project_level(fu_set_level(constant_pressure, fval1=85000.0), &
                          & vert_hgt, p, dummy)
    ! should be the same

    print *, 'Project hybrid(2) to small height', &
         & fu_project_level(fu_level(vert_hybrid, 2), vert_hgt_small, p_small, p_surf)
    ! hybrid(2) == 84.688717 m => ~2.347
    
    print *, 'Project sigma(2) to small height', &
         & fu_project_level(fu_level(vert_sigma, 2), vert_hgt_small, p_small, p_surf)
    ! should be the same
    
    print *, 'Testing input needs'
    call projection_input_needs(vert_hgt, vert_hybrid, q3d, q2d)
    print *, 'height to hybrid:', fu_quantity_short_string(q3d), ' ', fu_quantity_short_string(q2d)
    
    call projection_input_needs(vert_hybrid, vert_hybrid, q3d, q2d)
    print *, 'hybrid to hybrid:', fu_quantity_short_string(q3d), ' ', fu_quantity_short_string(q2d)

    call projection_input_needs(vert_hybrid, vert_hgt, q3d, q2d)
    print *, 'hybrid to height:', fu_quantity_short_string(q3d), ' ', fu_quantity_short_string(q2d)

    call projection_input_needs(vert_hybrid, vert_sigma, q3d, q2d)
    print *, 'hybrid to sigma:', fu_quantity_short_string(q3d), ' ', fu_quantity_short_string(q2d)

    call projection_input_needs(vert_hgt, vert_sigma, q3d, q2d)
    print *, 'height to sigma:', fu_quantity_short_string(q3d), ' ', fu_quantity_short_string(q2d)

    call projection_input_needs(vert_hgt, vert_press, q3d, q2d)
    print *, 'height to pressure:', fu_quantity_short_string(q3d), ' ', fu_quantity_short_string(q2d)

    call projection_input_needs(vert_hybrid, vert_press, q3d, q2d)
    print *, 'hybrid to pressure:', fu_quantity_short_string(q3d), ' ', fu_quantity_short_string(q2d)

    !.....

    print *, 'Testing crude interpolation'
    metdata => fu_work_array()
    call vert_interp_data_crude(vert_hybrid, vert_press, metdata, dummy)
    print *, 'Surface pressure:', dummy
    call vert_interp_data_crude(vert_hgt, vert_hybrid, metdata, dummy)
    print *, 'Heights:', metdata(1:nn)
    call vert_interp_data_crude(vert_hybrid, vert_hgt, metdata, dummy)
    print *, 'Pressures:', metdata(1:nn)
    call free_work_array(metdata)

  end subroutine test_project_level

  subroutine test_project_level_thick()
    implicit none
    
    integer, parameter :: nn = 41
    real, dimension(nn), parameter :: &
         & a = (/0.0, 2006.06, 3996.76, &
              & 5923.68, 7748.44, 9441.12, 10978.87, 12344.79, 13526.95, 14517.65, &
              & 15312.73, 15911.13, 16314.44, 16526.66, 16553.91, 16404.29, &
              & 16087.75, 15615.96, 15002.23, 14261.41, 13409.81, 12465.03, &
              & 11445.84, 10372.0, 9263.99, 8142.78, 7029.45, 5944.81, 4909.01, &
              & 3940.95, 3057.76, 2274.19, 1601.89, 1048.75, 618.02, 307.62, &
              & 109.2, 7.27, 0.0, 0.0, 0.0/), &
         
         & b = (/0.0, 0.0, 0.0, 0.0, 0.0, 0.01, 0.01, 0.02, 0.04, &
              0.05, 0.07, 0.09, 0.12, 0.15, 0.18, 0.22, 0.25, 0.29, 0.33, &
              0.38, 0.42, 0.46, 0.51, 0.55, 0.6, 0.64, 0.68, 0.72, 0.76, &
              0.79, 0.83, 0.86, 0.88, 0.91, 0.93, 0.95, 0.96, 0.97, 0.98, &
              0.99, 1.0/), &
              
         & sigma = (/0.0, 0.0197982940075, 0.0394449964396, &
                  & 0.0584622385405, 0.0764712387565, 0.103176709331, &
                  & 0.118353138057, 0.141833734724, 0.173500759262, &
                  & 0.193278218497, 0.221125056378, 0.247030811507, &
                  & 0.281011175981, 0.313105626772, 0.343374563649, &
                  & 0.381897927482, 0.408773917241, 0.444117706993, &
                  & 0.478060656366, 0.52074932362, 0.552344676114, &
                  & 0.583020412526, 0.622961778553, 0.652363790438, &
                  & 0.691428570284, 0.720363076119, 0.74937535159, &
                  & 0.778670775649, 0.808448213546, 0.828894193977, &
                  & 0.860177777078, 0.882444534186, 0.895809441985, &
                  & 0.920350368803, 0.936099389681, 0.953035976592, &
                  & 0.961077721357, 0.970071749398, 0.98, 0.99, 1.0/), &

         & height = (/0.0, 84.6887171594, 170.072974027, 255.545404244, &
                   & 333.588790243, 403.871373849, 553.477455535, 694.573138588, &
                   & 918.38286538, 1042.36468223, 1252.34487713, 1554.91984704, &
                   & 1757.70821799, 2060.60802189, 2367.90049014, 2681.96431296, &
                   & 3005.55876162, 3460.31113724, 3817.35648327, 4324.91930179, &
                   & 4734.21990025, 5175.50764187, 5807.54162922, 6343.58410412, &
                   & 6938.24363648, 7418.97056005, 8158.22902535, 8787.79381827, &
                   & 9511.66176116, 10355.0912502, 11063.805252, 11917.3727792, &
                   & 12601.9425516, 13879.9469436, 15027.6700619, 15897.9287197, &
                   & 17797.4242392, 19500.3626345, 22004.8436445, 26461.4559563, 100000.0/), &

         & p = (/101324.891894, 100311.642975, 99298.3940561, 98292.4151372, &
              & 97381.0962182, 96566.2672993, 94850.1694614, 93254.4016235, &
              & 90767.7948667, 89413.5970288, 87157.420272, 83987.6145963, &
              & 81915.9278394, 78898.7321637, 75930.3764879, 72990.7108122, &
              & 70058.9251364, 66100.6905417, 63121.5348659, 59074.4802712, &
              & 55966.2645955, 52764.8689197, 48439.444325, 45000.1786493, &
              & 41418.9729735, 38695.7662167, 34792.3905409, 31725.3937841, &
              & 28473.4270273, 25030.3702705, 22405.4724326, 19583.8945947, &
              & 17579.9456758, 14371.2878379, 11992.1189189, 10454.3689189, 7748.44, &
              & 5923.68, 3996.76, 2006.06, 0.0/)

    real, dimension(12) :: hgts = (/0.0, 50.0, 150.0, 300.0, 600.0, 1200.0, &
                                 & 2000.0, 3000.0, 4000.0, 5000.0, 6500.0, 8000.0/)

    real, dimension(12) :: p_small = (/101324.891894, 100725.675485, &
                                    & 99535.850101, 97772.4728, 94321.5833835, 87715.4879072, &
                                    & 79495.1305577, 70108.4696855, 61640.1693196, 54019.8542264, &
                                    & 44034.7988733, 35599.7731748/)

    real, dimension(nn) :: dummy_array
    real :: dummy
    real, parameter :: p_surf = 101324.891894
    integer :: iz, q2d, q3d
    real, dimension(:), pointer :: metdata

    type(silam_vertical) :: vert_hybrid, vert_sigma, vert_press, vert_hgt_small, vert_hgt, vert_press_small
    type(silja_level), dimension(nn) :: levels
    type(silja_level) :: lev, top, bottom
    real :: float

    do iz = 1, nn-1
      top = fu_set_level(hybrid, fval1=a(iz), fval2=b(iz))      
      bottom = fu_set_level(hybrid, fval1=a(iz+1), fval2=b(iz+1))
      levels(iz) = fu_set_layer_between_two(top, bottom)
    end do
    call set_vertical(levels, vert_hybrid)
    call arrange_levels_in_vertical(vert_hybrid)

    do iz = 1, nn-1
      bottom = fu_set_level(sigma_level, fval1=sigma(iz+1))
      top = fu_set_level(sigma_level, fval1=sigma(iz))
      levels(iz) = fu_set_layer_between_two(top, bottom)
    end do
    call set_vertical(levels, vert_sigma)
    call arrange_levels_in_vertical(vert_sigma)
    
    do iz = 1, nn-1
      top = fu_set_level(constant_pressure, fval1=p(iz+1))
      bottom = fu_set_level(constant_pressure, fval1=p(iz))
      levels(iz) = fu_set_layer_between_two(top, bottom)
    end do
    call set_vertical(levels, vert_press)
    call arrange_levels_in_vertical(vert_press)

    ! some arrays go up->down, some go down->up, sorry.
    do iz = 1, size(hgts) - 1
      top = fu_set_level(constant_height, fval1=hgts(iz+1))
      bottom = fu_set_level(constant_height, fval1=hgts(iz))
      levels(iz) = fu_set_layer_between_two(top, bottom)
    end do
    call set_vertical(levels(1:11), vert_hgt_small)
    call arrange_levels_in_vertical(vert_hgt_small)

    do iz = 1, nn-1
      top = fu_set_level(constant_height, fval1=height(iz+1))
      bottom = fu_set_level(constant_height, fval1=height(iz))
      levels(iz) = fu_set_layer_between_two(top, bottom)
    end do
    call set_vertical(levels, vert_hgt)
    call arrange_levels_in_vertical(vert_hgt)

    do iz = 1, 12-1
      top = fu_set_level(constant_pressure, fval1=p_small(iz+1))
      bottom = fu_set_level(constant_pressure, fval1=p_small(iz))
      levels(iz) = fu_set_layer_between_two(top, bottom)
    end do
    call set_vertical(levels(1:11), vert_press_small)
    print *, 'to arr'
    call arrange_levels_in_vertical(vert_press_small)
    print *, 'from arr, thick'
    !print *, 'Project hybrid(1) to hybrid', &
    !     & fu_project_level(fu_level(vert_hybrid, 1), vert_hybrid, dummy_array, p_surf)
    float = fu_project_level(fu_level(vert_hybrid, 1), vert_hybrid, dummy_array, p_surf)
    print *, 'out'
    print *, 'foo', fu_project_level(fu_level(vert_hybrid, 1), vert_hybrid, dummy_array, p_surf)
    print *, 'Project hybrid(3) to hybrid', &
         & fu_project_level(fu_level(vert_hybrid, 3), vert_hybrid, dummy_array, p_surf)

    
    print *, 'Project height(2) to hybrid', &
         & fu_project_level(fu_level(vert_hgt_small, 2), vert_hybrid, height, dummy)
    ! 1.59 seems to be correct
    
    print *, 'Project height(8) to hybrid', &
         & fu_project_level(fu_level(vert_hgt_small, 8), vert_hybrid, height, dummy)
    ! ~17 is ok
    
    print *, 'Project height(8) to sigma', &
         & fu_project_level(fu_level(vert_hgt_small, 8), vert_sigma, height, dummy)
    
    print *, 'Project 850 hPa to hybrid', &
         & fu_project_level(fu_set_level(constant_pressure, fval1=85000.0), &
                          & vert_hybrid, dummy_array, p_surf)
    ! should be about 11.68
    
    print *, 'Project 850 hPa to big height', &
         & fu_project_level(fu_set_level(constant_pressure, fval1=85000.0), &
                          & vert_hgt, p, dummy)
    ! should be the same

    print *, 'Project hybrid(2) to small height', &
         & fu_project_level(fu_level(vert_hybrid, 2), vert_hgt_small, p_small, p_surf)
    ! hybrid(2) == 84.688717 m => ~2.347
    
    print *, 'Project sigma(2) to small height', &
         & fu_project_level(fu_level(vert_sigma, 2), vert_hgt_small, p_small, p_surf)
    ! should be the same
    
    print *, 'Testing input needs'
    call projection_input_needs(vert_hgt, vert_hybrid, q3d, q2d)
    print *, 'height to hybrid:', fu_quantity_short_string(q3d), ' ', fu_quantity_short_string(q2d)
    
    call projection_input_needs(vert_hybrid, vert_hybrid, q3d, q2d)
    print *, 'hybrid to hybrid:', fu_quantity_short_string(q3d), ' ', fu_quantity_short_string(q2d)

    call projection_input_needs(vert_hybrid, vert_hgt, q3d, q2d)
    print *, 'hybrid to height:', fu_quantity_short_string(q3d), ' ', fu_quantity_short_string(q2d)

    call projection_input_needs(vert_hybrid, vert_sigma, q3d, q2d)
    print *, 'hybrid to sigma:', fu_quantity_short_string(q3d), ' ', fu_quantity_short_string(q2d)

    call projection_input_needs(vert_hgt, vert_sigma, q3d, q2d)
    print *, 'height to sigma:', fu_quantity_short_string(q3d), ' ', fu_quantity_short_string(q2d)

    call projection_input_needs(vert_hgt, vert_press, q3d, q2d)
    print *, 'height to pressure:', fu_quantity_short_string(q3d), ' ', fu_quantity_short_string(q2d)

    call projection_input_needs(vert_hybrid, vert_press, q3d, q2d)
    print *, 'hybrid to pressure:', fu_quantity_short_string(q3d), ' ', fu_quantity_short_string(q2d)

    !.....

    print *, 'Testing crude interpolation'
    metdata => fu_work_array()
    call vert_interp_data_crude(vert_hybrid, vert_press, metdata, dummy)
    print *, 'Surface pressure:', dummy
    call vert_interp_data_crude(vert_hgt, vert_hybrid, metdata, dummy)
    print *, 'Heights:', metdata(1:nn)
    call vert_interp_data_crude(vert_hybrid, vert_hgt, metdata, dummy)
    print *, 'Pressures:', metdata(1:nn)
    call free_work_array(metdata)

  end subroutine test_project_level_thick


END MODULE field_buffer
