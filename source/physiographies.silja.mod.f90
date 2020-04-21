MODULE physiographies ! topography, land-sea-data, albedo & stuff.
  !
  ! This module contains tools for getting physiographic land- and
  ! sea-surface data. Only permanent or semi-permanent (monthly) parameters 
  ! are handled here. 
  ! This module handles monthly meteoMarket, where all slowly changing 
  ! parameters are to be stored: water temperature, leaf area index, biomass volume, ...
  ! 
  ! Author: Mikhail Sofiev, FMI email Mikhail.Sofiev@fmi.fi
  !
  ! All units: SI, unless stated otherwise
  ! 
  ! All degrees: real numbers in degrees (degrees and parts per
  ! hundred/parts per thousand , NOT in degrees and minutes
  ! 
  ! Language: ANSI standard Fortran 90
  ! 
  USE supermarket_of_fields

  IMPLICIT NONE

  ! Init physiography mini-market for monthly and static fields
  !
  public init_physiography
  PUBLIC set_physiography
  public update_physiography
  public set_basic_physiography

  ! The public functions and subroutines available in this module:
  public write_physiography_output
  
  ! Generic function for handling the metadata related to land use
  public read_aggregate_land_use_metadata
  
  public make_time_zone_mapping   ! mapping of time zones as source indices

  ! Public field pointers:
  !
  ! ATTENTION. 
  ! In fact, below fields must not be here, or, at least, must be organised
  ! to some reasonble structure like meteobuffer. 
  ! So far I just allow their use elsewhere in the programm in order
  ! to speed-up the tasks.
  !
  TYPE(silja_field), POINTER, PUBLIC, SAVE :: topography_fld => null()
  TYPE(silja_field), POINTER, PUBLIC, SAVE :: geopotential_sfc_fld => null()
  TYPE(silja_field), POINTER, PUBLIC, SAVE :: albedo_fld => null()
  TYPE(silja_field), POINTER, PUBLIC, SAVE :: fraction_of_land_fld => null()
  TYPE(silja_field), POINTER, PUBLIC, SAVE :: water_eq_snow_depth_fld => null()
  TYPE(silja_field), POINTER, PUBLIC, SAVE :: ground_roughness_meteo_fld => null()
  TYPE(silja_field), POINTER, PUBLIC, SAVE :: ground_roughness_disp_fld => null()
!  TYPE(silja_field), POINTER, PRIVATE, SAVE :: water_roughness
  TYPE(silja_field), POINTER, PUBLIC, SAVE :: meteo_cell_x_size_fld => null() ! in m
  TYPE(silja_field), POINTER, PUBLIC, SAVE :: meteo_cell_y_size_fld => null()
  TYPE(silja_field), POINTER, PUBLIC, SAVE :: meteo_longitude_fld => null() ! in degrees
  TYPE(silja_field), POINTER, PUBLIC, SAVE :: meteo_latitude_fld => null()
  TYPE(silja_field), POINTER, PUBLIC, SAVE :: dispersion_cell_x_size_fld => null() ! in m
  TYPE(silja_field), POINTER, PUBLIC, SAVE :: dispersion_cell_y_size_fld => null()
  TYPE(silja_field), POINTER, PUBLIC, SAVE :: leaf_area_index_fld => null() ! m2/m2  only used with static LAI

  integer, private, dimension(17), parameter :: avail_physiography_quantities = &
                                            & (/relief_height_flag, &
                                            & climatological_albedo_flag, &
                                            & fraction_of_land_flag, &
                                            & fraction_of_ice_flag, &
                                            & water_eq_snow_depth_flag, &
                                            & climatological_albedo_flag, &
                                            & land_roughness_disp_flag, land_roughness_meteo_flag, &
                                            & leaf_area_index_flag, &
                                            & cell_size_x_flag, cell_size_y_flag, &
                                            & soil_sand_mass_fraction_flag, &
                                            & soil_clay_mass_fraction_flag, &
                                            & alluvial_sedim_index_flag, &
                                            & fraction_of_water_flag, &
                                            & longitude_flag, latitude_flag/)
  integer, private, save :: nPhysiographyFields
  TYPE(silja_fieldpointer), POINTER, dimension(:), private, SAVE :: fldPhysiography

  LOGICAL, PRIVATE, SAVE :: physiography_set = .false.

CONTAINS


  !******************************************************************
  
  subroutine init_physiography(meteoMarketPtr, wdr, static_shopping_list, &
                             & physiographyMarketPtr, pdr, monthly_shopping_list)
    !
    ! Initialises the physiography market for monthly fields - if any.
    ! Also sets physiographyDataRules as analogy for wdr, etc.
    ! Makes it ready 
    !
    implicit none

    ! Imported parameters    
    type(mini_market_of_stacks), pointer ::  meteoMarketPtr, physiographyMarketPtr
    type(silja_wdr), pointer :: wdr, pdr
    TYPE(silja_shopping_list), intent(in) :: static_shopping_list, monthly_shopping_list

    ! Local variables
    type(wdr_ptr), dimension(:), pointer :: wdrAr


    allocate(wdrar(1))
    wdrar(1)%ptr => pdr
    CALL initialize_mini_market(physiographyMarketPtr, &
                              & 'physiography_market', &
                              & fu_NbrOfMetSrcs(pdr), &  ! nbr of MDS
                              & 12, &   ! 12 timenodes = 12 months
                              & 100,& ! for each timenode - fields
                              & 0,&  ! for each timenode - windfields
                              & 5, &  ! for each timenode - 3d fields
                              & 0, &  ! for each timenode  - 3d windfields
                              & .true., & ! if replace oldest (true) or latest (false) when full
                              & wdrAr, &
                              & .false., &  ! if single src
                              & .true.) ! info to stdout <=> supermarket_info = .true.
    CALL initialize_mini_market(physiographyMarketPtr, &
                              & 'physiography_market', &
                              & 1, &  ! nbr of MDS
                              & 0,&   ! timenodes in memory (assume: max 1 timenode per file)
                              & 10,& ! for each timenode - fields
                              & 0,&  ! for each timenode - windfields
                              & 0, &  ! for each timenode - 3d fields
                              & 0, &  ! for each timenode  - 3d windfields
                              & .true.,& ! if replace oldest (true) or latest (false) when full
                              & wdrar, &
                              & .false., &
                              & .true.) ! info to stdout <=> supermarket_info = .true.
    call msg('physiography market is initialized')
    IF (error) RETURN

    deallocate(wdrar)
    
  end subroutine init_physiography


  ! ***************************************************************


  SUBROUTINE set_physiography(meteoMarketPtr, wdr, static_shopping_list, full_static_shopping_list, iAccuracy)
    ! 
    ! Finds surface information fields from external met_srcs.
    ! Quantities stored here must be either permanent
    ! (like topographic height) or semi-permanent (like fraction of ice).
    ! 
    ! If supermarket storage-area or storage-grid are defined, then they
    ! are used here as well, since the fields are stored in supermarket.
    !
    ! We make three steps: 
    ! 1. Read the static meteofile to get the requested fields
    ! 2. Make a few basic fields, such as grid cell size
    ! 3. If some fields are still missing, search the files once again
    !    taking any available grid and interpolating it into the meteo_grid
    !
    IMPLICIT NONE

    ! Imported parameters with intent IN:
    type(silja_wdr), intent(in) :: wdr
    TYPE(silja_shopping_list), intent(inout) :: static_shopping_list, full_static_shopping_list
    type(mini_market_of_stacks), pointer :: meteoMarketPtr
    integer, intent(in) :: iAccuracy

    ! Local declarations:
!    type(meteo_data_source) :: met_src
    TYPE(silja_shopping_list) :: shopping_list
    TYPE(silja_time) :: obstime
    TYPE(silja_field_id) :: id
    type(silam_sp), dimension(:), pointer :: fnames
    INTEGER :: i, j,n_req, n_avail, iSrc, n_shop, iTempl
    REAL, DIMENSION(:), POINTER :: arTmp, ptrTopo
    LOGICAL :: OK
    integer, dimension(max_quantities) :: Qs, q_shop, req
    type(silam_vertical) :: vertTmp
    type(silam_fformat)::fform
    character(len=fnlen), dimension(:), pointer :: descriptors
    type(silja_fieldpointer), dimension(:), allocatable :: output_fields

    !----------------------------------------
    !
    ! Check and start
    !
    IF (physiography_set) THEN
      CALL set_supermarket_empty(meteoMarketPtr, .true., .false.)
      call msg('')
      call msg('Re-setting physiography with the following vars:')
      physiography_set = .false.
    ELSE
      call msg('')
      call msg('Setting physiography with the following vars:')
    END IF

    call report(static_shopping_list)
    nullify(fnames)


    !
    ! Apart from the requested quantities, we will create a few below.
    ! In order to have consistent shopping_list and permanent meteo market
    ! have to add them now to the list
    !
    Qs(1:max_quantities) = int_missing
    Qs(1) = cell_size_x_flag
    Qs(2) = cell_size_y_flag
    Qs(3) = longitude_flag
    Qs(4) = latitude_flag
    req(1:max_quantities) = 2

    call add_shopping_quantities(static_shopping_list, Qs, req)
    call add_shopping_quantities(full_static_shopping_list, Qs, req)
    if(error)return

    Qs(1:4) = int_missing

    !
    ! Set the shopping_list twice, so that first physiography is searched
    ! in the meteo_grid and only then - all other sources are taken.
    ! This involves two reading of the input file, but on return ensures that
    ! fields with the meteo_grid will be certainly taken and others
    ! accepted only in case of absence of the meteo_grid based ones.

    call msg('Searching fields with the meteo_grid-corresponding grid')

    !
    ! Since the static_shopping_list may be in unknown grid, we have to copy
    ! only its quantities and create own shopping_list as explained above
    ! Note that we shop from the full list but later will check the availability 
    ! of only those, which are requested. The rest may be needed for the derivation
    !
    n_shop = fu_nbr_of_quantities(full_static_shopping_list)

    Q_shop(1:max_quantities) = fu_quantities(full_static_shopping_list)
    req(1:max_quantities) = fu_requests(full_static_shopping_list)
!    if(n_shop < size(Q_shop)) Q_shop(n_shop+1)=int_missing

    call set_missing(shopping_list)
    call set_missing(vertTmp, .true.)

    do i=1,size(Q_shop)
      if(Q_shop(i) == int_missing) exit
      call add_shopping_variable(shopping_list, &
                               & q_shop(i), species_missing, & ! quaintity, species
                               & meteo_grid, &
                               & vertTmp, &
                               & int_missing, & 
                               & met_src_missing, &
                               & req(i)) ! Requests
    end do
    call msg("Q_shop(1:n_req)",Q_shop(1:i-1))

    call consume_fields(shopping_list)
    


    !----------------------------------------
    !
    ! 2. Prepare some basic fields
    !
    ! Make the grid cell size - 2 more permanent fields - for meteo_grid
    !
!    print *, 'Starting cell sizes'

    id = fu_set_field_id(fu_met_src(wdr,1),&
                       & cell_size_x_flag, &
                       & fu_start_time(wdr),&
                       & zero_interval, &
                       & meteo_grid,&
                       & ground_level)
    call set_validity_length(id, very_long_interval)

    CALL find_field_storage_2d(meteoMarketPtr, id, single_time_stack_flag, meteo_cell_x_size_fld)
    IF(error)RETURN
    arTmp => fu_grid_data(meteo_cell_x_size_fld)

    DO j=1,ny_meteo ! Fill-in temporary array with x-size
     DO i=1,nx_meteo
       arTmp(i+(j-1)*nx_meteo)= fu_dx_cell_m(meteo_grid, i, j)
     END DO
    END DO

    ! meteo y-size  
    !
    CALL set_quantity(id,cell_size_y_flag)
    CALL find_field_storage_2d(meteoMarketPtr, id, single_time_stack_flag, meteo_cell_y_size_fld)
    IF(error)RETURN
    arTmp => fu_grid_data(meteo_cell_y_size_fld)
    DO j=1,ny_meteo   ! Fill-in temporary array with y-size
      DO i=1,nx_meteo
        arTmp(i+(j-1)*nx_meteo)= fu_dy_cell_m(meteo_grid, i, j)
      END DO
    END DO

    !
    ! Longitude
    !
    if(fu_quantity_in_quantities(longitude_flag, Q_shop))then

      CALL set_quantity(id, longitude_flag)
      CALL find_field_storage_2d(meteoMarketPtr, id, single_time_stack_flag, meteo_longitude_fld)
      IF(error)RETURN
      arTmp => fu_grid_data(meteo_longitude_fld)
      DO j=1,ny_meteo   ! Fill-in temporary array with y-size
        DO i=1,nx_meteo
          arTmp(i+(j-1)*nx_meteo) = fu_lon_geographical_from_grid(real(i), real(j), meteo_grid)
        END DO
      END DO
    else
      nullify(meteo_longitude_fld)
    endif
    !
    ! Latitude
    !
    if(fu_quantity_in_quantities(latitude_flag, Q_shop))then
      CALL set_quantity(id, latitude_flag)
      CALL find_field_storage_2d(meteoMarketPtr, id, single_time_stack_flag, meteo_latitude_fld)
      IF(error)RETURN
      arTmp => fu_grid_data(meteo_latitude_fld)
      DO j=1,ny_meteo   ! Fill-in temporary array with y-size
        DO i=1,nx_meteo
          arTmp(i+(j-1)*nx_meteo) = fu_lat_geographical_from_grid(real(i), real(j), meteo_grid)
        END DO
      END DO
    else
      nullify(meteo_latitude_fld)
    endif
    IF(error)RETURN


    !-------------------------------------------------------------------------
    !
    ! OK, now all meteo_grid-based fields are in permanent SM,
    ! the next step is to check if more fields needed
    !
    call supermarket_quantities(meteoMarketPtr, met_src_missing, single_time_stack_flag, Qs, n_avail) 
    !
    ! Now we shall search only those variables, which were not found before
    !
    do i=1,n_shop
      if(any(Qs(1:n_avail)==q_shop(i))) then
        q_shop(i)=int_missing
        req(i) = int_missing
      endif
    end do
    n_req = n_shop
    call compress_int_array(q_shop, int_missing, n_req) 
    n_req = n_shop
    call compress_int_array(req, int_missing, n_req)

    if(n_req > 0)then !n_shop)then

      call msg('')
      call msg('Searching other grids and interpolate to meteo_grid')

!call msg('6,',iTempl)
      !
      ! Shopping list accepts all what is valid inside the computation period - no variables, i.e.
      ! any grid and vertical are fine
      !
      call set_missing(shopping_list)
      call msg("Q_shop(1:n_req)",Q_shop(1:n_req))
      shopping_list = fu_set_shopping_list(met_src_missing, &
                                         & Q_shop(1:n_req), &
                                         & time_missing, & !fu_start_time(wdr), & ! First time boundary
                                         & time_missing, & !fu_start_time(wdr) + fu_period_to_compute(wdr), &
                                         & level_missing, &
                                         & level_missing, &
                                         & requests = req(1:n_req))
      
      call consume_fields(shopping_list)
      call supermarket_quantities(meteoMarketPtr, met_src_missing, single_time_stack_flag, Qs, n_avail) 
      
      
!!!!      do iTempl = 1, max_met_files
!!!!!call msg('1,',iTempl)
!!!!        if(fu_fname_oro_template(wdr, iTempl) == template_missing) exit
!!!!!call msg('2,',iTempl)
!!!!        !
!!!!        ! Depending on the format of the input data (file or field), we may need
!!!!        ! to expand the file names or simply copy the template output collection
!!!!        !
!!!!        if(fu_if_input_file(fu_file_oro_format(wdr,iTempl)))then
!!!!
!!!!!call msg('3,',iTempl)
!!!!          call FNm_from_single_template(fu_fname_oro_template(wdr, iTempl),&
!!!!                                     & fu_closest_obstime(fu_start_time(wdr), &
!!!!                                                        & back_and_forwards, &
!!!!                                                        & fu_obstime_interval(wdr)), &
!!!!                                     & fnames, &
!!!!                                     & ifAdd = .false., &
!!!!                                     & ifStrict = .false., &
!!!!                                     & ifAllowZeroFcLen = fu_if_zero_fc_len_allowed(wdr), &
!!!!                                     & ifWait = .false.)
!!!!!call msg('4,',iTempl)
!!!!          IF (error.or. fnames(1)%sp == '') then
!!!!            call set_error('Problem with one of input files','set_phisiography')
!!!!            RETURN
!!!!          end if
!!!!!call msg('5,',iTempl)
!!!!        else
!!!!          call enlarge_array(fnames,1)
!!!!          if(error)return
!!!!          fnames(1)%sp = fu_collection(fu_fname_oro_template(wdr, iTempl))
!!!!          if(size(fnames) > 1) fnames(2)%sp = ''
!!!!        endif
!!!!
!!!!!call msg('7,',iTempl)
!!!!        do i = 1, size(fnames)
!!!!          if(fnames(i)%sp == '')exit
!!!!          CALL store_input_to_supermarket(meteoMarketPtr, fnames(i)%sp, &
!!!!                                        & fu_file_oro_format(wdr,iTempl), &
!!!!                                        & shopping_list, &
!!!!                                        & wdr, &
!!!!                                        & create_if_absent, &        ! ifUpdateAllowed
!!!!                                        & single_time_stack_flag, 5, &
!!!!                                        & OK, &
!!!!                                        & pGrid_in = meteo_gridPtr, &
!!!!                                        & ifPermanent_ = .true.)
!!!!          IF (error) RETURN
!!!!          fnames(i)%sp = ''
!!!!        end do
!!!!        !
!!!!        ! May be, we get all the missing fields and there is no need to treat the other
!!!!        ! tempaltes ?
!!!!        ! 
!!!!!call msg('8,',iTempl)
!!!!        call supermarket_quantities(meteoMarketPtr, met_src_missing, single_time_stack_flag, Qs, n_avail) 
!!!!!call msg('9,',iTempl)
!!!!!        if(j == n_shop) exit 
!!!!!call msg('10,',iTempl)
!!!!
!!!!      end do ! Array of templates
    end if  ! If ther are needed variables not consumed with a meteo_grid limitation

    !----------------------------------------
    !
    ! Find the fields. A trick: there may be a few sources but they are sorted
    ! according to some rules, so that the first source is the most important
    ! So, we shall scan them from the 1-st until all (if success) fields are set.
    !
    do iSrc = 1, fu_NbrOfMetSrcs(wdr)
      OK = .true.

      if(fu_quantity_in_quantities(climatological_albedo_flag, Qs))then
        albedo_fld => fu_sm_simple_field(meteoMarketPtr, fu_met_src(wdr,iSrc),&
                               & climatological_albedo_flag,&  ! So far - the only way
                               & ground_level, &
                               & single_time_stack_flag)
        IF (error) then 
          CALL unset_error('set_physiography')
          OK = .false.
        else
          id = fu_id(albedo_fld)
        end if
      end if

      if(fu_quantity_in_quantities(land_roughness_meteo_flag, Qs))then
        ground_roughness_meteo_fld => fu_sm_simple_field(meteoMarketPtr, fu_met_src(wdr,iSrc),&
                                                 & land_roughness_meteo_flag,&
                                                 & level_missing, &
                                                 & single_time_stack_flag)
        IF (error) then 
          CALL unset_error('set_physiography')
          OK = .false.
        else
          id = fu_id(ground_roughness_meteo_fld)
        end if
      end if

      if(fu_quantity_in_quantities(land_roughness_disp_flag, Qs))then
        ground_roughness_disp_fld => fu_sm_simple_field(meteoMarketPtr, fu_met_src(wdr,iSrc),&
                                                 & land_roughness_disp_flag,&
                                                 & level_missing, &
                                                 & single_time_stack_flag)
        IF (error) then 
          CALL unset_error('set_physiography')
          OK = .false.
        else
          id = fu_id(ground_roughness_disp_fld)
        end if
      end if

      if(fu_quantity_in_quantities(water_eq_snow_depth_flag, Qs))then
        water_eq_snow_depth_fld => fu_sm_simple_field(meteoMarketPtr, fu_met_src(wdr,iSrc),&
                                           & water_eq_snow_depth_flag,&
                                           & level_missing, &
                                           & single_time_stack_flag)
        IF (error) then 
          CALL unset_error('set_physiography')
          OK = .false.
        else
          id = fu_id(water_eq_snow_depth_fld)
        end if
      end if

      if(fu_quantity_in_quantities(fraction_of_land_flag, Qs))then
        fraction_of_land_fld => fu_sm_simple_field(meteoMarketPtr, fu_met_src(wdr,iSrc),&
                                                 & fraction_of_land_flag,&
                                                 & level_missing, &
                                                 & single_time_stack_flag)
        IF (error) then 
          CALL unset_error('set_physiography')
          OK = .false.
        else
          id = fu_id(fraction_of_land_fld)
        end if
      end if


      !
      ! Funny but static geopotential means relief height, for which
      ! we have own flag.
      ! Procedure: if relief height is available directly, take it 
      ! if not, look for geopotential, from which we should derive the 
      ! relief height and store as a separate field.
      !
      if(fu_quantity_in_quantities(relief_height_flag, Qs))then
        !
        ! Relief height is available directly
        !
        topography_fld => fu_sm_simple_field(meteoMarketPtr, fu_met_src(wdr,iSrc),&
                                           & relief_height_flag,&
                                           & level_missing, &
                                           & single_time_stack_flag)
        IF (error) then 
          CALL unset_error('set_physiography')
          OK = .false.
        else
          id = fu_id(topography_fld)
        end if
      elseif(fu_quantity_in_quantities(geopotential_sfc_flag, Qs))then
        !
        ! Relief height is not available but surface geopotential is
        !
        geopotential_sfc_fld => fu_sm_simple_field(meteoMarketPtr, fu_met_src(wdr,iSrc),&
                                           & geopotential_sfc_flag,&
                                           & ground_level, &
                                           & single_time_stack_flag)
        IF (error) then 
          CALL unset_error('geopotential_sfc set_physiography')
          OK = .false.
        else
          !
          ! Get the geopotential field, scale with g and store it as a new one
          !
          id = fu_id(geopotential_sfc_fld)
          call set_quantity(id, relief_height_flag)
          CALL find_field_storage_2d(meteoMarketPtr, id, single_time_stack_flag, topography_fld)

          arTmp => fu_grid_data(geopotential_sfc_fld)  ! Input
          ptrTopo => fu_grid_data(topography_fld)  ! geopotential has to be scaled with g
          ptrTopo(1:fs_meteo) = arTmp(1:fs_meteo) / g
          call msg('Set topography field')
        end if
      endif

      if(OK) exit

    end do ! scan of the sources

    if(.not.OK)then
      call msg_warning('Not all permanent fields are found','set_physiography')
    end if

    ! If two LAIs come --  sum them up
    !
    if (fu_LAIsrc(wdr) == LAI_static_2) call  make_static_LAI(meteoMarketPtr, fu_if_randomise(wdr)) 

    !--------------------------------------------------------------------------------------------
    !
    ! The last step - read and process the land-use related variables.
    !
    descriptors => fu_land_use_descriptors(wdr)
    do i=1, size(descriptors)
      if(len_trim(descriptors(i)) == 0)exit  ! all done
      call read_aggregate_land_use_metadata(descriptors(i), wdr, shopping_list, output_fields)
      if(error)return
    end do  ! land_use_descriptors

    !----------------------------------------------------------------------------
    !
    ! After all - let's finally check that we got all mandatory fields.
    ! Note that this time we look for the quantities in the requested list, not full one.
    !
    Q_shop(1:max_quantities) = fu_quantities(static_shopping_list)
    req(1:max_quantities) = fu_requests(static_shopping_list)

    ! Here I silently assume that there is only one meteo source. If not, there will 
    ! be an error
    !
    call supermarket_quantities(meteoMarketPtr, met_src_missing, single_time_stack_flag, &
                              & Qs, nPhysiographyFields)
    allocate(fldPhysiography(nPhysiographyFields), stat=j)
    if(fu_fails(j==0,'Failed allocation of internal physiography pointers','set_physiography'))return
    do j = 1, nPhysiographyFields
      fldPhysiography(j)%fp => fu_get_2d_field(fu_stack(meteoMarketPtr, 1), j)
    enddo

    do i=1,size(q_shop)
      if(q_shop(i) == int_missing)exit
      if(fu_realtime_quantity(Q_shop(i))) cycle !No realtime quantities checked here
      OK = .false.
      do j=1, size(qs)
        if(Qs(j) == int_missing)exit
        if(qs(j) == q_shop(i))then
          OK = .true.
          exit
        endif
      end do
      if(.not. OK)then
        if(req(i) < 2)then
          call msg_warning('Desirable quantity not found:' + fu_quantity_short_string(q_shop(i)), &
                         & 'set_physiography')
        else
          call set_error('Mandatory quantity not found:' + fu_quantity_short_string(q_shop(i)), &
                       & 'set_physiography')
          return
        endif
      endif
    enddo
    
    physiography_set = .true.


    if(associated(fnames))then
      do i=1,size(fnames)
        deallocate(fnames(i)%sp)
      end do
      deallocate(fnames)
      nullify(fnames)
    endif

!    print *, 'END PHYSIOGRAPHY'


    CONTAINS

    subroutine consume_fields(shopping_list_local)
      !
      ! Scan all the input files for one time period - the information may be split
      !
      ! Imported parameters
      type(silja_shopping_list), intent(in) :: shopping_list_local

      ! Local variables
      integer :: iTempl
      
      do iTempl = 1, max_met_files
        if(fu_fname_oro_template(wdr, iTempl) == template_missing) exit
        !
        ! Depending on the format of the input data (file or field), we may need
        ! to expand the file names or simply copy the template output collection
        !
        if(fu_if_input_file(fu_file_oro_format(wdr,iTempl)))then

          call FNm_from_single_template(fu_fname_oro_template(wdr, iTempl),&
                                      & fu_closest_obstime(fu_start_time(wdr), &
                                                         & back_and_forwards, &
                                                         & fu_obstime_interval(wdr)), &
                                      & fnames, &
                                      & ifAdd = .false., &
                                      & ifStrict = .false., &
                                      & ifAllowZeroFcLen = fu_if_zero_fc_len_allowed(wdr), &
                                      & ifWait = .false.)
          IF (error) then
            call set_error('Problem with one of input files','set_physiography')
            RETURN
          end if
          if(associated(fnames))then
            if(fnames(1)%sp == '') then
              call set_error('Problem with one of input files','set_physiography')
              return
            end if
          else
            call set_error('names array is not associated','set_physiography')
            return
          endif
        else
          call enlarge_array(fnames,1)
          if(error)return
          fnames(1)%sp = fu_collection(fu_fname_oro_template(wdr, iTempl))
          if(size(fnames) > 1) fnames(2)%sp = ''
        endif  ! 
        !
        ! fnames now contain the set of either file names or instruction for field creation
        !
        do i = 1, size(fnames)
          if(fnames(i)%sp == '')exit
          CALL store_input_to_supermarket(meteoMarketPtr, fnames(i)%sp, fu_file_oro_format(wdr,iTempl), &
                                        & shopping_list_local, &
                                        & wdr, &
                                        & create_if_absent, &
                                        & single_time_stack_flag, &
                                        & iAccuracy, &
                                        & OK, &
                                        & storage_grid = meteo_grid, &
                                        & st_time_feature = static_climatology) !, ifVerbose_ = .true.)
          IF (error) RETURN
          fnames(i)%sp = ''
        end do
      end do  ! Cycle through templates

    end subroutine consume_fields

  END SUBROUTINE set_physiography


  !************************************************************************

  subroutine update_physiography(meteoMarketPtr, wdr, timeStart, timeEnd, iAccuracy )
    !
    ! Searches for the fields that are not valid for the particular time and updates them
    ! from the static_meteo_file
    ! Needed for monthly/weekly etc fields with irregular timing and update frequency
    ! lower than that of normal meteo fields. 
    !
    implicit none
    
    ! Imported parameters with intent IN:
    type(silja_wdr), intent(in) :: wdr
    TYPE(silja_time), intent(in) :: timeStart, timeEnd
    type(mini_market_of_stacks), pointer :: meteoMarketPtr
    integer, intent(in) :: iAccuracy

    ! Local variables
    integer :: iFld, iTempl, i
    type(silja_shopping_list) :: shopping_list_local
    logical :: ifProceed
    type(silam_sp), dimension(:), pointer :: fnames

    IF(.not. physiography_set) THEN
      CALL set_error('Updating the physiography only possible after it is set','')
      return
    END IF
    !
    ! Procedure:
    ! Scan all "static" fields checking if they are still valid for the requested interval
    ! If they are not, add them to the shopping list for updating.
    ! Note that the fields will be UPDATED, not read into a new place. This way we avoid 
    ! problems with space, duplicated field and necessity to reset all static pointers. 
    !
    call set_missing(shopping_list_local)
    if(error)return

    ifProceed = .false.
    do iFld = 1, nPhysiographyFields
      if(.not. fu_between_times(fu_valid_time(fldPhysiography(iFld)%fp), &
                              & timeStart, timeEnd, .true.))then
        !
        ! Add the ID to the shopping list for the update
        !
        call add_shopping_variable(shopping_list_local, fu_id(fldPhysiography(iFld)%fp), 2)
        if(error)return
        ifProceed = .true.
      endif
    end do
    !
    ! Anything to do?
    !
    if(.not. ifProceed)return  ! nothing needs updating here

    ! Let's do the work
    ! Scan all input files - the information may be split
    !
    do iTempl = 1, max_met_files
      if(fu_fname_oro_template(wdr, iTempl) == template_missing) exit
      !
      ! Depending on the format of the input data (file or field), we may need
      ! to expand the file names or simply copy the template output collection
      !
      if(fu_if_input_file(fu_file_oro_format(wdr,iTempl)))then

        call FNm_from_single_template(fu_fname_oro_template(wdr, iTempl),&
                                    & fu_closest_obstime(fu_start_time(wdr), &
                                                       & back_and_forwards, &
                                                       & fu_obstime_interval(wdr)), &
                                    & fnames, &
                                    & ifAdd = .false., &
                                    & ifStrict = .false., &
                                    & ifAllowZeroFcLen = fu_if_zero_fc_len_allowed(wdr), &
                                    & ifWait = .false.)
        IF (error) then
          call set_error('Problem with one of input files','update_physiography')
          RETURN
        end if
        if(associated(fnames))then
          if(fnames(1)%sp == '') then
            call set_error('Problem with one of input files','update_physiography')
            return
          end if
        else
          call set_error('names array is not associated','update_physiography')
          return
        endif
      else
        call enlarge_array(fnames,1)
        if(error)return
        fnames(1)%sp = fu_collection(fu_fname_oro_template(wdr, iTempl))
        if(size(fnames) > 1) fnames(2)%sp = ''
      endif  ! 
      !
      ! fnames now contain the set of either file names or instruction for field creation
      !
      do i = 1, size(fnames)
        if(fnames(i)%sp == '')exit

        CALL store_input_to_supermarket(meteoMarketPtr, fnames(i)%sp, fu_file_oro_format(wdr,iTempl), &
                                      & shopping_list_local, &
                                      & wdr, &
                                      & overwrite_forced, &      ! the field must exist
                                      & single_time_stack_flag, iAccuracy, &
                                      & ifProceed, &
                                      & storage_grid = meteo_grid, &
                                      & st_time_feature = static_climatology)! , ifVerbose_ = .true.)
        IF (error) RETURN
        fnames(i)%sp = ''
      end do
    end do  ! Cycle through templates


  end subroutine update_physiography


  !*****************************************************************

  subroutine write_physiography_output(meteoMarketPtr, chFNm_basic, iGrib, ifGrads, iNetCDF, now, &
                                     & filesWritten, output_grid, ifRandomise)
    !
    ! Writes the permanent stack to the output GRIB/GrADS files. Forces
    ! the field parameters and interpolates the grids if needed. Levels
    ! are copied one-to-one. The first-met field defines the time, the
    ! others use it.
    !
    implicit none

    ! Imported parameters
    character(len=*), intent(in) :: chFNm_basic
    logical, intent(in) :: ifGrads
    integer, intent(in) :: iGrib, iNetCDF
    type(silja_time), intent(in) :: now
    character(len=*), dimension(:), pointer :: filesWritten
    type(silja_grid), intent(in) :: output_grid
    type(mini_market_of_stacks), pointer :: meteoMarketPtr
    logical, intent(in) :: ifRandomise

    ! Local variables
    integer :: i, grib_funit, grads_funit, iNC, nFlds
    type(silja_stack), pointer :: permStack
    type(silja_field), pointer :: field_2d
    type(silja_grid) :: gridTmp
    type(silja_field_id) :: idTmp
    real, dimension(:), pointer :: dataPtr
    logical :: selection_done
    type(silam_vertical) :: vertical
    type(TOutputList) :: MeteoOutLstTmp
    !
    ! Some stupidity check
    !
    if(.not. (ifGrads .or. iNetCDF > 0 .or. iGrib > 0)) then
      call set_error('Neither GrADS nor GRIB nor NetCDF formats are active', 'write_physiography_output')
      return
    endif
    !
    ! Get the permanent stack
    !
    permStack => fu_stack(meteoMarketPtr, 1)
    if(.not.defined(permStack))then
      call set_error('Undefined permanent stack','write_physiography_output')
      return
    endif
    if(fu_number_of_2d_fields(permStack) < 1)then
      call set_error('Empty permanent stack','write_physiography_output')
      return
    endif

    ! Open the files
    !
    if(iGrib > 0)then                                           ! GRIB
      CALL open_gribfile_o('', &  ! directory - not used here
                         & chFNm_basic + '_physiography.grib', grib_funit) 
      if(error)return
    endif
    if(ifGrads)then                                             ! GRADS
      !
      ! For the multi-file output we should store the template rather than a file
      ! name into the GrADS structure. 
      !
      grads_funit = open_gradsfile_o('', &  ! directory not used 
                             & chFNm_basic + '_physiography.grads', output_grid)
      if(error)return
      call msg('Physiography: file opened for the following grid:')
      call report(output_grid)
    endif
    if(iNetCDF > 0)then                                            ! NetCDF
      !
      ! NetCDF requires a full list of output at the opening
      !
      call set_vertical(surface_level, vertical)
      if(error)return

      call expand_output_list(MeteoOutLstTmp, fu_number_of_2d_fields(permStack))
      if(error)return
      nFlds = 0
      DO i = 1, fu_number_of_2d_fields(permStack)
        field_2d => fu_get_2d_field(permStack,i) 
        IF(.not. ASSOCIATED(field_2d))cycle
        IF(.not. defined(field_2d))cycle
        if(.not. fu_quantity_in_quantities(fu_quantity(fu_id(field_2d)), &
                                         & avail_physiography_quantities))cycle
        nFlds = nFlds + 1
        MeteoOutLstTmp%ptrItem(i)%quantity = fu_quantity(fu_id(field_2d))
        MeteoOutLstTmp%ptrItem(i)%request = 2
        MeteoOutLstTmp%ptrItem(i)%AvType = iAsIs
        MeteoOutLstTmp%ptrItem(i)%if3D = .False.
        MeteoOutLstTmp%ptrItem(i)%iVerticalTreatment = do_nothing_flag
      end do ! cycle over permanent stack
      
      iNC = open_netcdf_file_o(chFNm_basic + '_physiography.nc', &  ! name
                             & output_grid, vertical, now, &  ! grid, vertical, valid_time
                             & (/MeteoOutLstTmp, &     ! Meteo variables to write
                              & OutputList_missing, &   !OutDef%Rules%DispOutLst, &      ! Dispersion stack variables
                              & OutputList_missing/), &              ! list of variables
                             & '', .true., iNetCDF, .false., real_missing)  ! chTemplate, ifAllInOne, nnc version, missingVal
    endif ! iNetCDF

    !----------------------------------
    !
    ! Scan the 2d fields (if any) - the full content of the stack, perform the
    ! interpolation to a temporary field and write to the requested files
    ! The first field decides time and all other not quantity-specific parameters
    !
    idTmp = fu_id(fu_get_2d_field(permStack,1))
    call set_grid(idTmp, output_grid)

    DO i = 1, fu_number_of_2d_fields(permStack)
      field_2d => fu_get_2d_field(permStack,i) 
      IF(.not. ASSOCIATED(field_2d))cycle
      IF(.not. defined(field_2d))cycle
      if(.not. fu_quantity_in_quantities(fu_quantity(fu_id(field_2d)), &
                                       & avail_physiography_quantities))cycle
      !
      ! Interpolate the field to the output_grid and form a fictive id
      !
      call set_quantity(idTmp, fu_quantity(fu_id(field_2d)))
      call set_level(idTmp, fu_level(fu_id(field_2d)))
      gridTmp=grid_missing
      !
      ! Data selection is done in two steps: (i) get the new grid and allocate space for ther data
      ! (ii) do the reprojection itself
      !
      call grid_data_hor_select_new_grid(fu_grid(field_2d), &  ! grid_original
                                       & output_grid, &        ! storage_grid
                                       & area_missing, &       ! storage_area
                                       & .false., &    ! ifAdjust to system_frid
                                       & gridTmp)      ! new grid_new
      if(error)return
      dataPtr => fu_work_array(fu_number_of_gridpoints(gridTmp))
      if(error)return
      
      call grid_data_horizontal_select(fu_grid(field_2d), &     ! grid_original
                                     & fu_grid_data(field_2d),& ! grid_data
                                     & selection_done, &  ! selection_done
                                     & gridTmp,&     ! grid_new
                                     & dataPtr, &    ! selected_grid_data
                                     & ifRandomise, &
                                     & 5, fu_regridding_method(fu_quantity(field_2d)), &   ! iAccuracy
                                     & real_missing, & !0.0, &          ! fMissingValue
                                     & setMissVal)
      if(selection_done .and. (gridTmp == output_grid))then
        !
        ! Actually write the fields
        !
        if(iGrib > 0) call write_next_field_to_gribfile(iGrib, grib_funit, idTmp, dataPtr)
        if(error)return
        if(ifGrads) call write_next_field_to_gradsfile(grads_funit, idTmp, dataPtr)
        if(error)return
        if(iNetCDF > 0)call write_next_field_to_netcdf_file(iNC, idTmp, dataPtr)
        if(error)return
      else
        if(gridTmp == output_grid)then
          call msg('Failed to reproject the field')
          call report(fu_id(field_2d))
          call set_error('Failed to reproject the field','write_physiography_output')
        else
          call msg('Failed to reproject the field: new grid is strange')
          call msg('reprojected grid:')
          call report(gridTmp)
          call msg('Requested grid:')
          call report(output_grid)
          call set_error('Failed to reproject the field','write_physiography_output')
        endif
        call unset_error('write_physiography_output')
      endif
      
      call free_work_array(dataPtr)

    END DO  ! 2d fields

    !
    ! Close the file and add the written file to the list
    !
    do i=1,size(filesWritten)
      if(filesWritten(i) == '')then
        nFlds = i
        exit
      endif
    end do
    i = 0
    if(ifGrads) i = i + 2
    if(iGrib /= int_missing) i = i + 2
    if(iNetCDF > 0) i = i + 1
    if(nFlds > size(filesWritten) - i)then
      call msg_warning('List of written files is full','write_physiography_output')
      nFlds = size(filesWritten) - i
    endif
    !
    ! Close files and store the names to the list of written files 
    !
    if(iGrib /= int_missing) then
      call close_gribfile_o(grib_funit)
      filesWritten(nFlds) = chFNm_basic + '_physiography.grib'
      nFlds = nFlds + 1
      filesWritten(nFlds+1) = chFNm_basic + '_physiography.grib.ctl'
      nFlds = nFlds + 1
    endif  ! ifGrib

    if(ifGrads)then
      call close_gradsfile_o(grads_funit,"")
      filesWritten(nFlds) = chFNm_basic + '_physiography.grads'
      nFlds = nFlds + 1
      filesWritten(nFlds+1) = chFNm_basic + '_physiography.grads.ctl'
      nFlds = nFlds + 1
    endif  ! ifGrads

    if(iNetCDF > 0 )then
      call close_netcdf_file(iNC)
      filesWritten(nFlds) = chFNm_basic + '_physiography.nc'
      nFlds = nFlds + 1
    endif   ! netcdf 

  end subroutine write_physiography_output


  !**********************************************************************************

  SUBROUTINE set_basic_physiography(miniMarketPtr, wdr, gridTemplate, verticalTemplate)
    ! 
    ! Make a few basic fields, such as grid cell size
    !
    IMPLICIT NONE

    ! Imported parameters with intent IN:
    type(silja_wdr), intent(in) :: wdr
    type(mini_market_of_stacks), pointer :: miniMarketPtr
    type(silja_grid), intent(in) :: gridTemplate
    type(silam_vertical), intent(in) :: verticalTemplate

    ! Local declarations:
    TYPE(silja_field_id) :: id
    INTEGER :: i,j, nx, ny
    REAL, DIMENSION(:), POINTER :: arTmp

    call grid_dimensions(gridTemplate, nx, ny)
    if(error)return
    arTmp => fu_work_array(nx*ny)
    if(error)return

    !----------------------------------------------------------------------------
    !
    ! Now make dx, dy size fields for the dispersion grid
    !
    DO j=1,ny     ! Fill-in temporary array with x-size
     DO i=1,nx
      arTmp(i+(j-1)*nx) = fu_dx_cell_m(gridTemplate, i,j)
     END DO
    END DO
    id = fu_set_field_id(silam_internal_src, &
                       & cell_size_x_flag, &
                       & fu_start_time(wdr), &
                       & zero_interval, &
                       & gridTemplate,&
                       & ground_level)
    call set_validity_length(id, very_long_interval)
    ! id, grid_data, permanent
    CALL dq_store_2d(miniMarketPtr, id, arTmp, single_time_stack_flag, fu_if_randomise(wdr), &
                   & storage_grid=gridTemplate)
    if(error)return

    DO j=1,ny         ! Fill-in temporary array with y-size
      DO i=1,nx
        arTmp(i+(j-1)*nx) = fu_dy_cell_m(gridTemplate, i,j)
      END DO
    END DO
    CALL set_quantity(id,cell_size_y_flag)

    CALL dq_store_2d(miniMarketPtr, id, arTmp, single_time_stack_flag, fu_if_randomise(wdr), &
                   & storage_grid=gridTemplate)
    if(error)return

    !-----------------------------------------------------------------------------
    !
    ! Depending on the type of vertical, we can set dz as a 3D static or dynamic field.
    ! Strictly speaking, there can also be dz as a 1D array of nz size but field-based
    ! concenpt does not support this trick.
    !
    if(.not. fu_if_level_meteo_dependent(fu_leveltype(verticalTemplate)))then
      !
      ! Static vertical: no need for meteo to get height from it
      !
      CALL set_quantity(id,cell_size_z_flag)

      do i = 1, fu_NbrOfLevels(verticalTemplate)
        arTmp(1:nx*ny) = fu_layer_thickness_m(fu_level(verticalTemplate, i))
        call set_level(id,fu_level(verticalTemplate, i))
        CALL dq_store_2d(miniMarketPtr, id, arTmp, single_time_stack_flag, fu_if_randomise(wdr), &
                       & storage_grid=gridTemplate)
      enddo

    endif

    call free_work_array(arTmp)

  end subroutine set_basic_physiography

  
  !****************************************************************************************
  subroutine make_static_LAI(meteoMarketPtr, ifRandomise) 
    !
    ! Updates static LAI from static LAIhv and LAIlv
    ! The resulting LAI is set valid as long as LAIlv
    ! Roux
    implicit none
    
    ! Imported parameters
    type(mini_market_of_stacks), pointer :: meteoMarketPtr
    logical, intent(in) :: ifRandomise
    
    ! Local variables
    REAL, DIMENSION(:), POINTER ::  lai, LAIhv, LAIlv, hvFrac, lVfrac
    TYPE(silja_field_id), POINTER :: laihv_field, lailv_field
    TYPE(silja_field_id) :: idTmp
    type(silja_field), pointer :: fldLAIhv, fldLAIlv, fldTmp
    integer :: iTmp, iPrDyn, fs, nx,ny

    call msg("make_static_LAI: Making static LAI out of LAI_hv and LAI_lv (with fractions)")

    if (.not. fu_field_in_sm(meteoMarketPtr, met_src_missing, leaf_area_indexhv_flag, &
                     & time_missing,  level_missing, .false., .true.)) then
            call msg("no LAI-hv in single-time stack. doing nothing")
            return
    endif

    if (.not. fu_field_in_sm(meteoMarketPtr, met_src_missing, leaf_area_indexlv_flag, &
                     & time_missing,  level_missing, .false., .true.)) then
            call msg("no LAI-lv in single-time stack. doing nothing")
            return
    endif

    if (.not. fu_field_in_sm(meteoMarketPtr, met_src_missing, fraction_hv_flag, &
                     & time_missing,  level_missing, .false., .true.)) then
            call msg("no frac-hv in single-time stack. doing nothing")
            return
    endif

    if (.not. fu_field_in_sm(meteoMarketPtr, met_src_missing, fraction_lv_flag, &
                     & time_missing,  level_missing, .false., .true.)) then
            call msg("no frac-lv in single-time stack. doing nothing")
            return
    endif


    fldLAIhv => fu_sm_simple_field(meteoMarketPtr, met_src_missing, leaf_area_indexhv_flag, &
                                & level_missing,  single_time_stack_flag)
    if(error)return
    LAIhv => fu_grid_data(fldLAIhv)

    fldLAIlv => fu_sm_simple_field(meteoMarketPtr, met_src_missing, leaf_area_indexlv_flag, &
                                & level_missing,  single_time_stack_flag)
    if(error)return
    LAIlv => fu_grid_data(fldLAIlv)

    fldTmp => fu_sm_simple_field(meteoMarketPtr, met_src_missing, fraction_hv_flag, &
                               & level_missing,  single_time_stack_flag)
    if(error)return
    hvFrac => fu_grid_data(fldTmp)

    fldTmp => fu_sm_simple_field(meteoMarketPtr, met_src_missing, fraction_lv_flag, &
                               & level_missing,  single_time_stack_flag)
    if(error) return
    lVfrac => fu_grid_data(fldTmp)


    idTmp = fu_id(fldLAIhv)

    ! Check that fields match
    call set_quantity(idTmp, leaf_area_indexlv_flag)
    if (.not. idTmp == fu_id(fldLAIlv)) then
       call msg("xxxxxxxxxx ID of LAI_hv xxxxxxxxxxxxxxx")
       call report(fu_id(fldLAIhv))
       call msg("xxxxxxxxxx ID of LAI_lv xxxxxxxxxxxxxxx")
       call report(fu_id(fldLAIlv))
       call set_error("Mismatch id's for LAI_hv and LAI_lv", "make_static_LAI")
       return
    endif


    !Target ID
    call set_quantity(idTmp, leaf_area_index_flag)

    call grid_dimensions(fu_grid(fldLAIlv), nx, ny)
    fs = nx*ny
    lai => fu_work_array(fs)
    lai(1:fs) = LAIhv(1:fs)*hvFrac(1:fs) + LAIlv(1:fs)*lvFrac(1:fs)
    CALL dq_store_2d(meteoMarketPtr, idTmp, lai, single_time_stack_flag, ifRandomise, &
                   & storage_grid=meteo_grid)
    call free_work_array(lai)

  end subroutine make_static_LAI

  
  !************************************************************************************************
  
  subroutine read_aggregate_land_use_metadata(chLandUseFileDescr, wdr, shopping_list, output_fields)
    !
    ! Reads the land-use typeSizes variable from the given file in accordance with its type and sets 
    ! it as a SILAM 2d field
    !
    ! Relates the given land use field_LandUse and aggregates it using the rules described in
    ! the handling nlMetadata namelist.
    ! Two aggregation possibilities are available today:
    ! 1. Lumping the land-use types into a new field. Grid etc are preserved since re-categorization 
    !    does not allow spatial operations
    ! 2. Creation of derived duantities that are averageable, with possible regridding
    !
    ! The metadata file will look like this:
    ! LIST = land_use_file
    !   land_use_file = GRADS <file_name>
    ! END_LIST = land_use_file
    ! LIST = landuse_types_lumping                ! lumping rules, same for all new vars
    !   land_use = 1 forest
    !   land_use = 2  ...
    ! END_LIST = landuse_types_lumping
    ! LIST = aggreagted_variables                 ! list of below namelists
    !   variable_list_name = soil_water_capacity_srf_list
    !   variable_list_name = ...
    ! END_LIST = aggreagted_variables
    ! LIST = soil_water_capacity_srf_list             ! one list per aggregated variable
    !   quantity = water_capac_soil_srf_grav_flag
    !   species =                                     ! can be empty if not needed for this quantity
    !   unit =                                        ! can be empty, then no scaling is applied
    !   grid = same / meteo / dispersion / output     ! forced "same" if aggregation_type == LUMP_CLASSES
    !   variable_value = <lumped_landuse>  <value>
    !   variable_value = forest 10
    ! END_LIST = soil_depth_list
    !
    implicit none
    
    ! Imported parameters
    character(len=*), intent(in) :: chLandUseFileDescr  ! File name with namelist inside
    type(silja_wdr), intent(in) :: wdr
    type(silja_shopping_list), intent(in) :: shopping_list
    type(silja_fieldpointer), dimension(:), allocatable, intent(out) :: output_fields
    
    ! Local variables
    type(Tsilam_namelist_group), pointer :: nlGrpMetadata
    type(Tsilam_namelist), pointer :: nlPtr
    type(Tsilam_nl_item_ptr), dimension(:), pointer :: pItems, pNewVarList
    integer :: iItem, nItems, iStat, iLU, ix, iy, ii, nxFrom, nyFrom, nxTo, nyTo, nLumped, &
             & iVar, nNewVars, iValue, fsFrom, fsTo, quantity, ii_new, iID, nIDs, uIn
    real :: fValue, factor, fX, fY
    type(silam_sp) :: sp
    logical :: ifLumping, ifAverage, ifRegrid, ifFound
    integer, parameter :: nLumpedMax = 1000
    integer, dimension(nLumpedMax) :: indLumping
    character(len=clen), dimension(nLumpedMax) :: chLumpedTypes
    real, dimension(nLumpedMax) :: fValueLumpedLU
    character(len=clen) :: chLumpedLU
    real, dimension(:), pointer :: output_data
    integer, dimension(:), pointer :: iCounter
    type(silam_species) :: species
    type(silja_grid) :: gridFrom, gridTo
    type(silja_field_id) :: idOut
    type(silja_field_id), dimension(:), pointer :: idList
    type(silam_fformat) :: fileFormat
    real, dimension(:), pointer :: dataLU
    type(silja_field_id) :: idLU

    !
    ! Read the file
    !
    call msg('read_aggregate_land_use_metadata is reading the land use data. Can take time!')
    sp%sp => fu_work_string()
    !
    ! Open the file and find the needed land-use file
    !
    ii = fu_next_free_unit()
    open(ii,file=fu_process_filepath(chLandUseFileDescr), iostat = iStat)
    if(fu_fails(iStat==0,'Failed to open file:'+chLandUseFileDescr,'read_aggregate_land_use_metadata'))return
    
    nlGrpMetadata => fu_read_namelist_group(ii,.true.,'END_LAND_USE_DESCRIPTOR')
    if(error)return
    
    close(ii)
    
    nlPtr => fu_namelist(nlGrpMetadata,'land_use_file')
    if(fu_fails(associated(nlPtr),'No land_use_file namelist in:'+chLandUseFileDescr,'read_aggregate_land_use_metadata'))return
    
    sp%sp = fu_content(nlPtr,'land_use_file')  ! get the file format and name 
    
    fileFormat = fu_input_file_format(sp%sp)
    call open_input_file(sp%sp(index(sp%sp,' ')+1:), fileFormat, uIn, int_missing)
    if(error)return
    !
    ! List of field IDs in the input file.
    !
    call get_id_list_from_input_file(uIn, fileFormat%iFormat, idList, nIDs) !! Only first time returned and used
    if(error)return
    !
    ! Scan the input IDs in the file and find the one with the quantity land_use_type_flag
    ! Its grid is then used to request the space
    !
    ifFound = .false.
    do iID = 1, nIDs
      if(fu_quantity(idList(iID)) == land_use_type_flag)then
        dataLU => fu_work_array(fu_number_of_gridpoints(fu_grid(idList(iID))))
        if(error)return
        call get_input_field(sp%sp(index(sp%sp,' ')+1:), fileFormat, & ! File where from to take it
                           & idList(iID), &           ! The stuff to search
                           & dataLU, &
                           & grid_missing, &          ! storage grid
                           & .true., notAllowed, &    ! ifAcceptSameMonth, iOutside, 
                           & idLU, &
                           & 5, &              !iAccuracy, &
                           & wdr, &
                           & fFillValue_=real_missing, &
                           & iBinary=uIn)
        if(error)return
        ifFound = .true.
        exit
      endif  ! found landuse
    end do  ! scan the available IDs looking for the landuse
    if(.not. ifFound) call set_error('Failed to find landuse in:' + chLandUseFileDescr,'read_land_use_data')

    !
    ! Lumping rules for land use types. This is to reduce the numner of effective landuse types 
    ! thus allowing only a few metadata classes, e.g. grass, forest, rest. If this namelist is absent
    ! existing metadata are asssumed for each land use type available from the idLU and dataLU
    !
    nlPtr => fu_namelist(nlGrpMetadata,'landuse_types_lumping')
    if(associated(nlPtr))then
      ifLumping = .true.
      call get_items(nlPtr,'land_use', pItems, nItems)
      if(error .or. fu_fails(nItems > 0, 'No land_use lines in the landuse_types_lumping namelist', &
                                       & 'read_aggregate_land_use_metadata'))return
      if(fu_fails(nItems < nLumpedMax, &
                & 'Insufficient temporary array size. Too many land use typeSizes:' + fu_str(nItems), &
                & 'read_aggregate_land_use_metadata'))return
      chLumpedTypes = ''
      indLumping = int_missing
      !
      ! Scan the lumping lines and fill-in arrays: chLumpedTypes with lumped names and
      ! indLumped with index of the specific landuse (from the main idLU and dataLU) in that lumped-names array
      !
      nLumped = 0
      do iItem = 1, nItems
        sp%sp = fu_content(pItems(iItem))
        read(unit=sp%sp,fmt=*) iLU, chLumpedTypes(nItems+1)  ! safe owing to above check
        do iStat = 1, nLumped
          if(chLumpedTypes(iStat) == chLumpedTypes(nItems+1))then
            indLumping(iLU) = iStat   ! linking the in-file and lumped LUs
            iLU = int_missing
            exit
          endif
        end do  ! search through already existing lumped types
        if(iLU /= int_missing)then  ! have not found this lumped type. Make a new
          nLumped = nLumped + 1
          indLumping(iLU) = nLumped
          chLumpedTypes(nLumped) = chLumpedTypes(nItems+1)  ! add 
        endif
      end do  ! search through all LU types in the LU field
    else
      ifLumping = .false.
    endif  ! if lumping rules namelist is present in the namelist group

    ! Get the list of variables based on the lumped landuse. 
    !
    nlPtr => fu_namelist(nlGrpMetadata,'aggregated_variables')
    if(fu_fails(associated(nlPtr),'Absent aggregated_variables namelist','read_aggregate_land_use_metadata'))return

    call get_items(nlPtr,'variable_list_name', pNewVarList, nNewVars)
    if(fu_fails(nNewVars > 0,'No variable_list_name items in aggregated_variables','read_aggregate_land_use_metadata'))return
    allocate(output_fields(nNewVars), stat=iStat)
    if(fu_fails(iStat==0,'Failed allocation of output field array. Size=' + fu_str(nNewVars),'read_aggregate_land_use_metadata'))return
    do iVar = 1, nNewVars
      allocate(output_fields(iVar)%fp, stat=iStat)
      if(fu_fails(iStat==0,'Failed allocation of output field pointer nrb:'+fu_str(iVar),'read_aggregate_land_use_metadata'))return
    end do
    !----------------------------------------------------------------------------------
    !
    ! The grand cycle over the required variables to be based on the land use lumping
    ! For each, create new variable either in the same grid or with reprojection
    !
    do iVar = 1, nNewVars
      nlPtr => fu_namelist(nlGrpMetadata, fu_content(pNewVarList(iVar)))
      if(fu_fails(associated(nlPtr),'Absent variable list:'+fu_content(pNewVarList(iVar)),'read_aggregate_land_use_metadata'))return
      !
      ! Find the lumping relation and store the lumped land use classification
      !
      call get_items(nlPtr,'variable_value',pItems,nItems)  ! for each lumped landuse type
      if(fu_fails(nItems > 0,'No variable_value items in namelist:'+fu_content(pNewVarList(iVar)),'read_aggregate_land_use_metadata'))return
      fValueLumpedLU(:) = real_missing
      do iItem = 1, nItems
        sp%sp = fu_content(pItems(iItem))
        read(unit=sp%sp, fmt=*) chLumpedLU, fValue
        ifFound = .false.
        do iStat = 1, nLumped
          if(trim(fu_str_u_case(chLumpedLU)) == trim(fu_str_u_case(chLumpedTypes(iStat))))then
            fValueLumpedLU(iStat) = fValue
            ifFound = .true.
            exit
          endif
        end do
        if(.not. ifFound)call msg_warning('Cannot find lumped LU:' + chLumpedLU,'read_aggregate_land_use_metadata')
      end do
      !
      ! All active land-use typeSizes must have the coefficient defined
      !
      do iStat = 1, nLumped
        if(abs(fValueLumpedLU(iStat) - real_missing) < 0.1 * abs(real_missing))then
          call set_error('Undefined coefficient for lumped class:' + chLumpedTypes(iStat) + &
                       & ', variable:' + fu_content(pNewVarList(iVar)),'read_aggregate_land_use_metadata')
        endif
      end do
      if(error)return
      
      !
      ! Action depends on the aggregation_type
      !
      sp%sp = fu_str_u_case(fu_content(nlPtr,'aggregation_type'))
      fsFrom = fu_number_of_gridpoints(fu_grid(idLU))

      if(trim(sp%sp) == 'LUMP_CLASSES')then
        !
        ! Proceed with the output field. No regridding, same ID, quantity, unit, etc.
        !
        if(.not. fu_quantity_in_list(land_use_type_flag, shopping_list))cycle
        
        call set_field(idLU, output_fields(iVar)%fp, .true., ifResetVals=.false.)
        if(error)return
        output_data => fu_grid_data(output_fields(iVar)%fp)
        
        if(ifLumping)then
          do ii = 1, fsFrom
            output_data(ii) = fValueLumpedLU(nint(dataLU(ii))) ! Actually, integers but no such field type
          enddo
        else
          output_data(1:fsFrom) = dataLU(1:fsFrom)
        endif

      elseif(trim(sp%sp) == 'REGRID_SUM_NEW_QUANTITY' .or. &
           & trim(sp%sp) == 'REGRID_AVERAGE_NEW_QUANTITY')then
        !
        ! Get the conversion details: unit, quantity, grid. make the new empty field
        !
        gridFrom = fu_grid(idLU)
        call grid_dimensions(gridFrom, nxFrom, nyFrom)

        quantity = fu_get_silam_quantity(fu_content(nlPtr,'quantity'))
        if(.not. fu_known_quantity(quantity)) &
                & call set_error('Unknown quantity:' + fu_content(nlPtr,'quantity'),'read_aggregate_land_use_metadata')
        if(.not. fu_quantity_in_list(quantity, shopping_list))cycle

        if(len_trim(fu_content(nlPtr,'species')) > 0)then
          species = fu_species_from_short_string(fu_content(nlPtr,'species'))
          if(error .or. .not. defined(species))then
            call set_error('Undefined species:' + fu_content(nlPtr,'species'),'read_aggregate_land_use_metadata')
            return
          endif
          factor = fu_conversion_factor(fu_content(nlPtr,'unit'), fu_quantity_unit(quantity), &
                                      & fu_material(species))
        else
          species = species_missing
          factor = fu_conversion_factor(fu_content(nlPtr,'unit'), fu_quantity_unit(quantity))
        endif
        if(error .or. abs(factor-real_missing) < 0.1 * abs(real_missing))then
          call set_error('Strange conversion factor or species:' + fu_content(nlPtr,'species'),'read_aggregate_land_use_metadata')
          return
        endif
        ifRegrid = .true.
        if(trim(fu_str_u_case(fu_content(nlPtr,'grid'))) == 'SAME_GRID')then
          ifRegrid = .false.
          gridTo = gridFrom
        elseif(trim(fu_str_u_case(fu_content(nlPtr,'grid'))) == 'METEO_GRID')then
          gridTo = meteo_grid
        elseif(trim(fu_str_u_case(fu_content(nlPtr,'grid'))) == 'DISPERSION_GRID')then
          gridTo = dispersion_grid
        elseif(trim(fu_str_u_case(fu_content(nlPtr,'grid'))) == 'OUTPUT_GRID')then
          gridTo = output_grid
        else
          call set_error('Unknown grid type:'+fu_content(nlPtr,'grid') + &
                       & ', allowed only SAME_GRID, METEO_GRID, DISPERSION_GRID, OUTPUT_GRID', &
                       & 'read_aggregate_land_use_metadata')
          return
        endif
        ifAverage = (trim(sp%sp) == 'REGRID_AVERAGE_NEW_QUANTITY')
        !
        ! Set the output field id and data
        !
        idOut = idLU
        call set_quantity(idOut, quantity)
        call set_grid(idOut, gridTo)
        call set_species(idOut, species)
        call set_field(idOut, output_fields(iVar)%fp, .true., ifResetVals=.false.)
        if(error)return
        output_data => fu_grid_data(output_fields(iVar)%fp)
        !
        ! A possibility here is to use the interpolation structure but we do not need any split of 
        ! the input grid to small cells, so it is faster to make it manually.
        !
        ! Cycle over the input grid. Project each point and check where it ends.
        !
        if(ifRegrid)then
          call grid_dimensions(gridTo, nxTo, nyTo)
          fsTo = nxTo * nyTo
          output_data(1:fsTo) = 0.
          if(ifAverage)then
            iCounter = fu_work_int_array(fsTo)  ! of the output grid
            if(error)return
            iCounter(1:fsTo) = 0
          endif
          ii = 0
          do iy = 1, nyFrom
            do ix = 1, nxFrom
              ii = ii + 1        ! counter is faster than ii = ix + (iy-1)*nx
              call project_point_to_grid(gridFrom,real(ix),real(iy), gridTo,fX,fY)
              if(fX < 0.5 .or. fX > nxTo+0.5 .or. fY < 0.5 .or. fY > nyTo+0.5)cycle  ! outside the new grid
              ii_new = nint(fX) + (nint(fY) - 1) * nxTo
              if(ifLumping)then
                output_data(ii_new) = output_data(ii_new) + &
                                    & fValueLumpedLU(indLumping(nint(dataLU(ii)))) * factor
              else
                output_data(ii_new) = output_data(ii_new) + &
                                    & fValueLumpedLU(nint(dataLU(ii))) * factor
              endif
              if(ifAverage) iCounter(ii_new) = iCounter(ii_new) + 1
            end do
          end do
          if(ifAverage)then
            where(iCounter(1:fsTo) > 0) output_data(1:fsTo) = output_data(1:fsTo) / iCounter(1:fsTo)
          endif
        else
          ! No regridding/averaging, just new quantity + lumping. But no vector ops: stack overflow
          if(ifLumping)then
            do ii = 1, fsFrom
              output_data(ii) = fValueLumpedLU(indLumping(nint(dataLU(ii)))) * factor
            end do
          else
            do ii = 1, fsFrom
              output_data(ii) = fValueLumpedLU(nint(dataLU(ii))) * factor
            end do
          endif
        endif   ! if regrid
      else
        
        call set_error('Only LUMP_CLASSES of REGRID_(SUM/AVERAGE)_NEW_QUANTITY actions, not:' + sp%sp,'read_aggregate_land_use_metadata')

      endif  ! type of action for the new variable
    
    end do  ! cycle over new variables

    call free_work_array(sp%sp)
    call free_work_array(dataLU)

  end subroutine read_aggregate_land_use_metadata
  
  
  !****************************************************************************
  
  subroutine make_time_zone_mapping(arSourceIdMapping, chSplitNames, meteoMarketPtr, gridMapping, groups, nGroups)
    !
    ! Maps the time zones and their groups as a filled set of indices.
    ! Called once, so efficiency does not matter
    !
    implicit none
    
    ! imported parameters
    integer, dimension(:,:), intent(inout) :: arSourceIdMapping  ! map of indices
    character(len=10), dimension(:), intent(inout) :: chSplitNames    ! list of names assigned to indices
    type(mini_market_of_stacks), pointer :: meteoMarketPtr
    type(silja_grid), intent(in) :: gridMapping
    type(Tsilam_nl_item_ptr), dimension(:), pointer :: groups  ! just items in the namelist
    integer, intent(in) :: nGroups
    
    ! local variables
    integer :: iGrp, ix, iy, indexMeteo, nx, ny, iTZ, iZone
    type(silja_field), pointer :: TZidxFld
    real, dimension(:), pointer :: tz_data
    type(THorizInterpStruct), pointer :: pHIS
    character(len=10), dimension(:,:), allocatable :: groupItems
    integer, dimension(:), pointer :: nZonesInGrp
    character(len=10) :: cellCode, cellCode2
    logical :: ifFound
    !
    ! Get the time zone map and prepare supplementary variables
    !
    TZidxFld => fu_sm_simple_field(meteoMarketPtr, met_src_missing, timezone_index_flag, &
                                 & level_missing, single_time_stack_flag)
    tz_data => fu_grid_data(TZidxFld)
    ! Note two different nearest-points. One is for interpolation method, one for out-of-grid handling
    pHIS => fu_horiz_interp_struct(meteo_grid, gridMapping, nearest_point, .false., 5, nearestPoint) 
    if(error)return
    ! Handle the zone groups: read the namelist items and decode it
    nZonesInGrp => fu_work_int_array()
    if(error)return
    allocate(groupItems(number_of_time_zones, nGroups), stat = iGrp)
    if(iGrp /= 0)then
      call set_error('Failed allocation of temporary group list, nGrp=' + &
                   & fu_str(nGroups) + ', nTZ=' + fu_str(number_of_time_zones),'make_time_zone_mapping')
      return
    endif
    chSplitNames(:) = ''
    do iGrp = 1, nGroups
      call split_string(fu_content(groups(iGrp)), ' ', groupItems(:,iGrp), nZonesInGrp(iGrp))
      if(error)return
    end do
    !
    ! Scan the given grid and fill it in with indices. Make sure that indices are sequential,
    ! have no holes and account for groupping of time zones
    !
    call grid_dimensions(gridMapping, nx, ny)
    do iy = 1, ny
      do ix = 1, nx
        indexMeteo = fu_grid_index(nx_meteo, ix, iy, pHIS)
        iTZ = nint(tz_data(indexMeteo))   ! index of the time zone
        cellCode = fu_str_u_case(TZ_code(iTZ))   ! country name
        cellCode2 = cellCode + '_' + fu_str(iTZ)  ! country name and the time zone index
        !
        ! is it in any of the groups?
        ! Note that one can use codes (3-char abbreviation of the country, e.g. FIN) or
        ! the code with sub-index, e.g. RUS_25 if the country has nore than one time zone
        !
        ifFound = .false.
        do iGrp = 1, nGroups  
          do iZone = 2, nZonesInGrp(iGrp)  ! first item is the group name
            if(trim(cellCode) == trim(groupItems(iZone,iGrp)) .or. &
             & trim(cellCode2) == trim(groupItems(iZone,iGrp)))then
              cellCode = groupItems(1, iGrp)  ! cell gets the group name as the code
              ifFound = .true.
              exit
            endif
          end do  ! zones in the group
          if(ifFound)exit
        end do  ! groups
        !
        ! Having the code of this grid cell, find its index or create a new
        !
        do iZone = 1, number_of_time_zones
          if(chSplitNames(iZone) == cellCode)then  
            exit                               ! known name, pick the index
          elseif(chSplitNames(iZone) == '')then
            chSplitNames(iZone) = cellCode     ! all known names ended, make new
            exit
          endif
        end do  ! cycle over zones with already registered titles
        !
        ! Got that!
        !
        arSourceIdMapping(ix, iy) = iZone
      enddo  ! ix
    enddo  ! iy
    
  end subroutine make_time_zone_mapping
  
END MODULE physiographies
