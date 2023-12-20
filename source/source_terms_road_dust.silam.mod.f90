MODULE source_terms_road_dust
  !
  ! This module contains description of the road dust emission.
  !
  ! All units: SI, unless otherwise stated
  !
  ! Language: ANSI FORTRAN-90 (or close to it)
  !
  use source_terms_time_params !cocktail_basic
  use cocktail_basic

  implicit none
  private

  !
  ! PUBLIC routines of road dust dust
  !
  public fill_rd_dust_src_from_namelist
  public reserve_rd_dust_source
  public init_emission_rd_dust
  public add_source_species_rd_dust_src
  public create_source_containing_grid
  public add_input_needs
  public link_source_to_species
  public source_2_second_grid
  public compute_emission_for_rd_dust
  public fu_rd_dust_emis_owned_quantity
  public fu_name
  public fu_source_nbr
  public fu_source_id_nbr
  public typical_species_conc
  public report

  !
  ! Private routines of the road dust source
  !
  private create_src_cont_grd_rddust_src
  private add_input_needs_rd_dust_src
  private link_rd_dust_src_to_species
  private fu_species_index_src
  private fu_source_id_nbr_of_rd_dust_src
  private fu_source_nbr_of_rd_dust_src
  private fu_rd_dust_source_name
  private typical_species_cnc_rd_dust
  private report_rd_dust_src
  private update_road_wetness_and_snow

  !
  ! Private subs of the road dust source
  !
  interface create_source_containing_grid
    module procedure create_src_cont_grd_rddust_src
  end interface

  interface add_input_needs
    module procedure add_input_needs_rd_dust_src
  end interface

  interface link_source_to_species
    module procedure link_rd_dust_src_to_species
  end interface

  interface source_2_second_grid
    module procedure project_rddust_src_second_grd
  end interface

  interface fu_source_nbr
    module procedure fu_source_nbr_of_rd_dust_src
  end interface

  interface fu_source_id_nbr
    module procedure fu_source_id_nbr_of_rd_dust_src
  end interface

  interface fu_name
    module procedure fu_rd_dust_source_name
  end interface

  interface typical_species_conc
    module procedure typical_species_cnc_rd_dust
  end interface

  interface report
    module procedure report_rd_dust_src
  end interface

  !
  ! There might be several types of the emission algorithm....
  !
  integer, private, parameter :: simple_road_dust_flag = 6501

  integer :: indSnowdays = int_missing
  integer :: indwetness = int_missing
  
  TYPE silam_road_dust_source
    PRIVATE
    CHARACTER(len=clen) :: src_nm, sector_nm       ! Name of the area source and sector
    character(len=fnlen) :: src_data_dir
    integer :: emisMethod, iSpectrumType, src_nbr, id_nbr  ! A source and id numbers in a WHOLE source list
    integer :: nLevsDispVert, nSpecies
    type(silam_vertical) :: vertLevsDispVert
    real, dimension(:), pointer :: levFractDispVert, fzDisp
    type(Tsilam_namelist), pointer :: nlInputFiles  ! namelist for names of supplementary files
    type(silam_species), dimension(:), pointer :: species
    type(silja_field), pointer :: road_wetness_field  => NULL(), n_snowless_days_field  => NULL()
    type(silja_logical) :: defined
  END TYPE silam_road_dust_source

  type rd_dust_src_ptr
    type(silam_road_dust_source) :: rddust_src
  end type rd_dust_src_ptr
  public rd_dust_src_ptr


CONTAINS


  !*********************************************************************

  subroutine fill_rd_dust_src_from_namelist(nlSetup, srcRDDust, expected_species, src_data_dir)
    
    implicit none

    ! Imported parameters
    type(Tsilam_namelist), pointer :: nlSetup
    type(silam_road_dust_source), intent(inout) :: srcRDDust
    type(silam_species), dimension(:), intent(in), allocatable :: expected_species
    character(len=*), intent(in) :: src_data_dir

    ! Local variables
    integer :: iTmp, jTmp, iSpecies, nFiles
    type(Tsilam_nl_item_ptr), dimension(:), pointer ::  pItems
    type(Taerosol) :: aerosolTmp
    type(silam_species), dimension(1) :: species
    integer, dimension(:), pointer :: indices
    character(len=fnlen) :: chTmp
    logical :: ifFound

    !
    ! Names
    !
    srcRDDust%src_nm = fu_content(nlSetup,'source_name')
    srcRDDust%sector_nm = fu_content(nlSetup,'source_sector_name')
    srcRDDust%src_data_dir = src_data_dir
    srcRDDust%defined = silja_false

    !
    ! Emission index type
    !
    select case(fu_str_u_case(fu_content(nlSetup,'road_dust_emission_method')))

      case ('SIMPLE_ROAD_DUST')
        srcRDDust%emisMethod = simple_road_dust_flag

      case default
        call set_error('Unknown emission method:' + &
                     & fu_content(nlSetup,'road_dust_emission_method'), &
                     & 'fill_rd_dust_src_from_namelist')
        return
    end select

    !
    ! Now the list of aerosol modes that will be emitted. To set them up, we will
    ! use the standard aerosol procedure - for the sake of unification.
    ! Note that there are two potentially concurring definitions: one coming from aerosol dynamics
    ! the other - written in the ini file. The first one prevails, if it exists
    !
    call get_source_aer_species(srcRDDust%species, srcRDDust%nSpecies, &
                              & expected_species, fu_content(nlSetup,'road_dust_substance_name'), &
                              & nlSetup)
    if(error)return

    !
    ! Store the input fields for the source features.
    !
    srcRDDust%nlInputFiles => fu_create_namelist('road_dust_src_supplementary_files')
    if(error)return

    pItems => null()
    call get_items(nlSetup, 'supplementary_file', pItems, nFiles)
    if(error)return
    !
    ! No grads hat expansion or alike: the string contains the file format
    !
    do iTmp = 1, nFiles
      chTmp = fu_process_filepath(fu_content(pItems(iTmp)),superdir=srcRDDust%src_data_dir)
      call add_namelist_item(srcRDDust%nlInputFiles, 'supplementary_file', chTmp)
    end do
    !
    ! Note that the source scaling field cannot be read from the file yet: dispersion grid is undefined
    !
    chTmp = fu_process_filepath(fu_content(nlSetup,'source_scaling_field'),superdir=srcRDDust%src_data_dir)
    call add_namelist_item(srcRDDust%nlInputFiles, 'source_scaling_field', chTmp)

    if(error)return

    srcRDDust%defined = silja_undefined

  end subroutine fill_rd_dust_src_from_namelist


  !*****************************************************************

  subroutine add_input_needs_rd_dust_src(rd_dust_src, q_met_dynamic, q_met_st, &
                                                  & q_disp_dynamic, q_disp_st, wdr)
    !
    ! Returns input needs for the emission/transformation routines. 
    !
    implicit none

    ! Imported parameters
    type(silam_road_dust_source), intent(in) :: rd_dust_src
    integer, dimension(:), intent(inout) :: q_met_dynamic, q_met_st, &
                                          & q_disp_dynamic, q_disp_st
    type(silja_wdr), intent(in), optional :: wdr

    ! Local variables
    integer :: iTmp

    !
    ! Add needed dynamic quantities. Always precipitation, the rest depends on the method
    !
    !iTmp = fu_merge_integer_to_array(total_precipitation_int_flag,  q_met_st)
    !iTmp = fu_merge_integer_to_array(water_eq_snow_depth_flag, q_met_dynamic)

    select case(rd_dust_src%emisMethod)

      case(simple_road_dust_flag)
         iTmp = fu_merge_integer_to_array(temperature_2m_flag, q_met_dynamic) 
         iTmp = fu_merge_integer_to_array(relative_humidity_2m_flag, q_met_dynamic)
         iTmp = fu_merge_integer_to_array(water_eq_snow_depth_flag, q_met_dynamic)
         iTmp = fu_merge_integer_to_array(n_snowless_days_flag, q_disp_dynamic)
         iTmp = fu_merge_integer_to_array(road_wetness_flag, q_disp_dynamic)
         iTmp = fu_merge_integer_to_array(total_precipitation_int_flag, q_met_st)
         iTmp = fu_merge_integer_to_array(road_dust_emis_fact_flag, q_disp_st)
      case default
        call set_error('Unknown emission computation algorithm', 'input_needs_rd_dust_source')
        return
    end select

  end subroutine add_input_needs_rd_dust_src


  !**************************************************************************

  subroutine reserve_rd_dust_source(rd_dust_src, &     ! Src to initialise
                                   & iSrcNbr, &      ! Src number in the grand list
                                   & iSrcIdNbr)      ! SrcID number
    !
    ! Initialises the source:
    ! - stores the reference information: source number and source ID number
    ! - stores the total number of chemical descriptors that will be stored in the source
    ! - nullifies the source dynamic arrays 
    ! - set a few basic internal variables of the source
    !
    implicit none

    ! Imported parameters
    type(silam_road_dust_source), intent(inout) :: rd_dust_src
    integer, intent(in) :: iSrcNbr, iSrcIdNbr

    !
    ! Nullify the basic variables
    !
    rd_dust_src%src_nm = ''
    rd_dust_src%sector_nm = ''
    !
    !
    nullify(rd_dust_src%fZDisp)
    nullify(rd_dust_src%levFractDispVert)
    !
    ! Main source parameters - enough to identify it in the global information list
    !
    rd_dust_src%src_nbr = iSrcNbr
    rd_dust_src%id_nbr = iSrcIdNbr
    !
    ! Finally, mark the source as incomplete
    !
    rd_dust_src%defined = silja_false

  end subroutine reserve_rd_dust_source


  !*********************************************************************

  subroutine init_emission_rd_dust(srcRDDust, dispersionMarketPtr, start_time)
    !
    ! Initializes the road dust source term.
    ! The parameters are read from the given ini file and stored to the returned 
    ! source term.
    !
    implicit none

    ! Imported parameters
    type(Tsilam_namelist), pointer :: nlSetup
    type(silam_road_dust_source), intent(inout) :: srcRDDust
    type(mini_market_of_stacks), pointer :: dispersionMarketPtr
    type(silja_time), intent(in) :: start_time

    ! Local variables
    integer :: iTmp, iSpecies, iFlds, iZ, nQuantities
    real, dimension(:), pointer :: pValues
    type(silam_sp) :: sp
    real :: fTmp
    type(silja_field_id) :: id, idTmp
    real, dimension(:), pointer :: arPtr
    type(silja_shopping_list) :: shop_list
    integer, dimension(:), pointer ::  q_disp_dyn, q_disp_stat, stack_quantities
    type(silam_vertical) :: vertTmp
    logical :: ifOK
    type(silja_field), pointer :: fieldPtr
    type(silja_field), pointer :: field

    q_disp_dyn => fu_work_int_array()
    q_disp_stat => fu_work_int_array()
    stack_quantities => fu_work_int_array()
    q_disp_dyn(1:max_quantities) = int_missing
    q_disp_stat(1:max_quantities) = int_missing
    stack_quantities(1:max_quantities) = int_missing
          
    !
    ! Storing the quantities needed in the dispersion stack
    ! The meteo stack has already been requested
    !
    call add_input_needs_rd_dust_src(srcRDDust, &
                                   & stack_quantities, stack_quantities, & ! meteo quantities, skip
                                   & q_disp_dyn, q_disp_stat)   ! dispersion-buffer quantities, use
    if(error)return
    !
    ! Make the shopping list with all needed quantities. Note that we can have several 
    ! pollen sources with the same or different species emitted - e.g. grass and birch pollen.
    ! The only way to allow them in one run is to use the species names for the dispersion stack 
    ! fields. If the species names are same for some sources, no problem, it is just the same source
    ! written in several parts.
    !
    call set_missing(shop_list)
    call set_missing(vertTmp, .true.)
    if(error)return

    do iFlds = 1, size(q_disp_stat)
      if(q_disp_stat(iFlds) == int_missing)exit
      iSpecies = fu_species_index_src(srcRDDust, q_disp_stat(iFlds))
      if(iSpecies == int_missing)then                             ! Universal quantity
        call add_shopping_variable(shop_list, &
                                 & q_disp_stat(iFlds), &     ! quantity
                                 & species_missing, &
                                 & grid_missing, &
                                 & vertTmp, int_missing, &
                                 & met_src_missing)
      else                                                   ! Species-specific quantity
        call add_shopping_variable(shop_list, &
                                 & q_disp_stat(iFlds), &             ! quantity
                                 & srcRDDust%species(iSpecies), &
                                 & grid_missing, &
                                 & vertTmp, int_missing, &
                                 & met_src_missing)
      endif
      if(error)return
    end do

    !
    ! The supplementary fields have to be taken from supplementary_info files. Note
    ! that due to species names used explicitly in the shopping list, they must be in the files.
    !
    call msg('Filling-in the road dust emission info from supplementary fields')
    call fill_minimarket_from_namelist(dispersionMarketPtr, &
                                     & srcRDDust%nlInputFiles, 'supplementary_file', & ! namelist and item
                                     & shop_list, start_time, &
                                     & static_climatology, &  ! target stack
                                     & create_field, &            ! error if a clash
                                     & wdr_missing, &
                                     & dispersion_gridPtr, &
                                     & 5, .true., & ! iAccuracy, ifAdjustGrid
                                     & ifOK)

    id = fu_set_field_id_simple(met_src_missing,&
         & road_dust_emis_fact_flag, &
         & time_missing, &        ! valid time
         & level_missing)
    call get_field_from_mm_general(dispersionMarketPtr, id, fieldPtr, .false.)
    if(.not. associated(fieldPtr))then
       call set_error('Road dust emission factor map is missing but needed','init_emission_rd_dust')
       return
    endif

    !call msg('finding or creating fields for road dust')
    
    call find_or_create_field(road_wetness_flag, dispersionMarketPtr, species_missing, start_time, surface_level, ifOK, field, pValues, 1.0)
    if(error)return
    
    call find_or_create_field(n_snowless_days_flag, dispersionMarketPtr, species_missing, start_time, surface_level, ifOK, field, pValues, 1000.0)
    if(error)return

    !CALL supermarket_2d_quantities(dispersionMarketPtr, &
    !                             & met_src_missing, multi_time_stack_flag, &
    !                             & stack_quantities, nQuantities)
    
    srcRDDust%defined = silja_true

    call free_work_array(stack_quantities)
    call free_work_array(q_disp_dyn)
    call free_work_array(q_disp_stat)

    call report(srcRDDust)

  end subroutine init_emission_rd_dust


  !****************************************************************************

  subroutine add_source_species_rd_dust_src(rd_dust_src, species_list, nSpecies)
    !
    ! Get the species emitted by this source and make a list of
    ! them. The source must be initialized, of course.
    !
    implicit none
    type(silam_road_dust_source), intent(in) :: rd_dust_src
    type(silam_species), dimension(:), pointer :: species_list
    integer, intent(inout) :: nSpecies

    call addSpecies(species_list, nSpecies, rd_dust_src%species, rd_dust_src%nSpecies)

  end subroutine add_source_species_rd_dust_src


  !**************************************************************************

  subroutine link_rd_dust_src_to_species(species_list, rd_dust_src)
    !
    ! Having the cocktails created, we should establish the shortcut links between the 
    ! source descriptors and cocktail. The link goes via descr%iEmisCocktSpeciesMapping 
    ! and  descr%factor_to_basic_unit./ descr%factor_foreign_2_basic_unit
    ! Note that cocktail is "anything" and the only requirement is that it has all the species
    ! existing in the source inventory.
    ! That has to happen in two steps. Firstly, we establish these links using the 
    ! single descriptor per source. Then, these connections are distributed to each time slot
    !
    implicit none

    ! Imported parameters
    type(silam_road_dust_source), intent(inout) :: rd_dust_src
    type(silam_species), dimension(:), pointer :: species_list
    
  end subroutine link_rd_dust_src_to_species


  !************************************************************************

  integer function fu_species_index_src(rd_dust_src, quantity)
    !
    ! Selects the proper name of the substance for the given quantity - just to be able to 
    ! choose the right input and dispersion-stack fields for each source.
    ! Indices of main pollen species, pollen and free allergen are known by the source
    !
    implicit none

    ! Imported parameters
    type(silam_road_dust_source), intent(in) :: rd_dust_src
    integer, intent(in) :: quantity

    fu_species_index_src = int_missing  ! so far, nothing fancy
        
  end function fu_species_index_src


  !*******************************************************************************
  
  subroutine create_src_cont_grd_rddust_src(rd_dust_src, grid_template, ifVerbose, ifMinimal, ifExtended)
    !
    ! Creates the grid that covers the area with active BVOC emission
    !
    implicit none

    ! Imported parameters
    type(silam_road_dust_source), intent(in) :: rd_dust_src
    type(silja_grid), intent(inout) :: grid_template
    logical, intent(in) :: ifVerbose, ifMinimal
    logical, intent(out) :: ifExtended

    ! So far nothing to do: this source just covers the dispersion grid
    !
    ifExtended = .false.
    return
    
  end subroutine create_src_cont_grd_rddust_src


  !*****************************************************************

  subroutine project_rddust_src_second_grd(rddust_src, grid, vert_disp, vert_proj, iAccuracy)
    !
    ! There is nothing to project in terms of horizontal grid but the vertical has to be made.
    !
    implicit none

    ! Imported parameters
    type(silam_road_dust_source), intent(inout) :: rddust_src
    type(silja_grid), intent(in) :: grid
    type(silam_vertical), intent(in) :: vert_disp, vert_proj
    integer, intent(in) :: iAccuracy

    ! Local variables
    integer :: i, ilev
    type(silam_vertical) :: vertTmp
    real, dimension(:), pointer :: arTmp, dz_disp

    rddust_src%vertLevsDispVert = vert_disp
    allocate(rddust_src%levFractDispVert(fu_NbrOfLevels(vert_disp)), &
           & rddust_src%fzDisp(fu_NbrOfLevels(vert_disp)), stat=i)
    if(fu_fails(i==0, 'Failed to allocate dispersion-vertical level fractions','project_rddust_src_second_grd'))return
    rddust_src%levFractDispVert(:) = 0.0
    rddust_src%fzDisp(:) = 0.0

    arTmp => fu_work_array()
    if(error)return

    call set_vertical(fu_set_layer_between_two(layer_btw_2_height, 0.0, 5.0), vertTmp)
    if(error)return
    arTmp(1) = 1.0

    call reproject_verticals(vertTmp, arTmp, &                   ! vertical from, fractions from
                           & vert_proj, rddust_src%levFractDispVert, &   ! vertical to, fractions to
                           ! mass centres, number of non-zero levels:
                           & rddust_src%fzDisp, rddust_src%nLevsDispVert, &
                           & ifMassCentreInRelUnit=.true.)

    call set_missing(vertTmp, .false.)
    call free_work_array(arTmp)
    
  end subroutine project_rddust_src_second_grd


  !**************************************************************************

  subroutine compute_emission_for_rd_dust(rd_dust_src, &
                                        & met_buf, disp_buf, & 
                                        & now, &      ! current time
                                        & timestep, & ! model time step
                                        & pHorizInterpMet2DispStruct, ifHorizInterp, &
                                        & ifSpeciesMoment, &
                                        & emisMap, mapCoordX, mapCoordY, mapCoordZ, & ! Output
                                        & fMassInjected)                              ! Output
    !
    ! Computes the emission fields for road_dust.
    !
    ! This routine is to be called at each model time step but not inside the main cycle,
    ! therefore its efficiency is of moderate importance.
    !
    implicit none

    ! Imported parameters
    type(silam_road_dust_source), target, intent(in) :: rd_dust_src
    type(Tfield_buffer), pointer ::  met_buf, disp_buf  ! meteo and internal field buffers
    type(silja_time), intent(in) :: now           ! current time
    type(silja_interval), intent(in) :: timestep  ! model time step
    type(THorizInterpStruct), intent(in) :: pHorizInterpMet2DispStruct
    logical, intent(in) :: ifHorizInterp, ifSpeciesMoment
    type(Tmass_map), intent(inout) :: emisMap, mapCoordX, mapCoordY, mapCoordZ
    real(r8k), dimension(:), intent(inout) :: fMassInjected
    character(len=substNmLen) :: substance_name
    
    ! Local variables
    integer :: indLandFr, iDust, iLev,iThresh, iZ, iSp, iSrc, &
         & iMeteo, iDisp, ix, iy, ixMeteo, iyMeteo, iStat, iMode, iTmp, jTmp, iMeteoTmp, &
         & indSnowlessDays, indwetness, indprecrate, indrh2m, indtemp2m, indweqsnowd, indrdemisfact

    real, dimension(:), pointer :: pEmisFact
    real :: fTmp, timestep_sec, fCellTotal, fFluxTotal, snow_depth, scaling_factor, wetness, snowless_days, t2m
    !real :: snow_scaling = 2.5 ! v2
    !real :: snow_scaling = 5.0 !v1
    !real :: snow_scaling = 8.7 !v3,v4
    !real :: snow_scaling = 6.0 !v5
    real :: snow_scaling = 4.4 !v6
    
    real, dimension(:), pointer ::  prh2m, pt2m, pPrecrate, &
                                 & pLandFr, pSnowlessDays, pRoadWetness, &
                                 & pSnowWEQDepth, ptrLonMeteo, ptrLatMeteo
    real*4, dimension(worksize) :: arTmp
    type(silja_field), pointer :: fldMaskPtr
    integer, save :: iCount = 0, iCountTalking = 0
    logical, save :: ifFirst=.true., ifTalking=.true.

    type(silam_sp) :: sp
    
    indweqsnowd = fu_index(met_buf, water_eq_snow_depth_flag, pSnowWEQDepth)
    indrh2m = fu_index(met_buf, relative_humidity_2m_flag, prh2m)
    indtemp2m = fu_index(met_buf, temperature_2m_flag, pt2m)
    indprecrate = fu_index(met_buf, total_precipitation_int_flag, pPrecrate)

    indSnowlessDays = fu_index(disp_buf, n_snowless_days_flag, species_missing, .true.)
    indwetness = fu_index(disp_buf, road_wetness_flag, species_missing, .true.)
    !indrdemisfact = fu_index(disp_buf, road_dust_emis_fact_flag, species_missing, .true.)

    if(fu_fails(fu_index(disp_buf, road_dust_emis_fact_flag, pEmisFact) /= int_missing, &    ! basic map for road dust emission
             & 'Failed road_dust_emis_fact_flag','compute_emission_for_rd_dust'))return

    !call msg('indSnowlessDays', indSnowlessDays)
    !call msg('indwetness', indwetness)
    
    !indwetness = fu_index(disp_buf, road_wetness_flag, pRoadWetness)

    !fu_index(buf, quantity, fu_species_src(src, quantity, ifSpeciesMandatory), .true.)
    !indHS = fu_get_buffer_index(srcPollen, disp_buf, heatsum_flag, .true.)
    
    !ptrXSzDisp => fu_grid_data(dispersion_cell_x_size_fld)  ! get grid cell size in metres
    !ptrYSzDisp => fu_grid_data(dispersion_cell_y_size_fld)

    !ptrLonMeteo => fu_grid_data(meteo_longitude_fld)
    ptrLatMeteo => fu_grid_data(meteo_latitude_fld)

    timestep_sec = abs(fu_sec(timestep))
    !if(error)return

    call update_road_wetness_and_snow(met_buf, disp_buf, indprecrate, indrh2m, indweqsnowd, &
                                 & indSnowlessDays, indwetness, one_day, timestep_sec, now, timestep)

    !indSnowlessDays = fu_index(disp_buf, n_snowless_days_flag, pSnowlessDays)
    !indwetness = fu_index(disp_buf, road_wetness_flag, pRoadWetness)
    
    !
    ! Now scale the emission map
    !
    do iSp = 1, emisMap%nSpecies

       if (.not.  fu_substance_name(emisMap%species(iSp)) == 'road_dust') cycle

       do iSrc = 1, emisMap%nSrc
       
          do iy = 1, ny_dispersion
             do ix = 1, nx_dispersion
                iDisp = ix + (iy-1) * nx_dispersion
        
                iMeteo =  fu_grid_index(nx_meteo, ix, iy, pHorizInterpMet2DispStruct)
                !iyMeteo = int(iMeteo / nx_meteo)
                !ixMeteo = iMeteo - (iyMeteo-1)*nx_meteo

                !if (ptrLatMeteo(iMeteo) > 55) then
                !   latscaling = 0.75
                !else
                !   latscaling = 0.15
                !end if
                !if (pLandFr(iMeteo) < 1.0e-5) cycle
                
                select case(rd_dust_src%emisMethod)
                case(simple_road_dust_flag)

                   !call msg('latitude', ptrLatMeteo(iMeteo))
                   !call msg('road wetness', disp_buf%p2d(indwetness)%present%ptr(iDisp))
                   !call msg('n_snowless_days', disp_buf%p2d(indSnowlessDays)%present%ptr(iDisp))

                   ! v1
                   !if (disp_buf%p2d(indSnowlessDays)%present%ptr(iDisp) > 10 .or. ptrLatMeteo(iMeteo) < 55.0)then
                   !   scaling_factor = 0.1 * (1 - disp_buf%p2d(indwetness)%present%ptr(iDisp))
                   !else
                   !   scaling_factor = snow_scaling * 0.1 * (1 - disp_buf%p2d(indwetness)%present%ptr(iDisp))
                   !end if

                   wetness = disp_buf%p2d(indwetness)%present%ptr(iDisp)
                   t2m = met_buf%p2d(indtemp2m)%present%ptr(iMeteo)
                   snowless_days = disp_buf%p2d(indSnowlessDays)%present%ptr(iDisp)

                   ! v3,v4,v5,v6
                   !if (disp_buf%p2d(indSnowlessDays)%present%ptr(iDisp) > 10 .or. ptrLatMeteo(iMeteo) < 55.0)then
                   if (snowless_days > 10) then
                      if (t2m <= 273) then
                        scaling_factor = 0.65 * (1 - wetness)**10
                      else
                        scaling_factor = 0.65 * (1 - wetness)
                      end if
                   else
                      if (t2m <= 273) then
                         scaling_factor = (1 + pEmisFact(iDisp)*snow_scaling) * 0.65 * (1 - wetness)**10
                      else
                         scaling_factor = (1 + pEmisFact(iDisp)*snow_scaling) * 0.65 * (1 - wetness)
                      end if
                   end if

                   ! v2
                   !if (ptrLatMeteo(iMeteo) < 55.0)then
                   !   scaling_factor = 0.1 * (1 - disp_buf%p2d(indwetness)%present%ptr(iDisp))
                   !else
                   !   if (disp_buf%p2d(indSnowlessDays)%present%ptr(iDisp) > 10) then
                   !      scaling_factor = 0.5 * (1 - disp_buf%p2d(indwetness)%present%ptr(iDisp))
                   !   else
                   !      scaling_factor = snow_scaling * 0.5 * (1 - disp_buf%p2d(indwetness)%present%ptr(iDisp))
                   !   end if
                   !end if

                   !call msg('disp_buf%p2d(indSnowlessDays)%present%ptr(iDisp), disp_buf%p2d(indwetness)%present%ptr(iDisp)', disp_buf%p2d(indSnowlessDays)%present%ptr(iDisp), disp_buf%p2d(indwetness)%present%ptr(iDisp))
                   
                case default
                   call msg('wrong method')
                end select
                !
                ! Fill-in the emission map, not forgetting the number to mass convertion where needed
                !
                do iLev = 1, emisMap%n3D
                  
                   if (emisMap%arM(iSp, iSrc, iLev, ix, iy) < 1e-8) cycle  ! nothing for this dispersion layer
            
                   emisMap%arM(iSp, iSrc, iLev, ix, iy) = scaling_factor * &
                        & emisMap%arM(iSp, iSrc, iLev, ix, iy) !**0.20  !dbg 0.25
                   fMassInjected(iSp) = fMassInjected(iSp) + emisMap%arM(iSp, iSrc, iLev, ix, iy)

                   mapCoordZ%arm(iSp, iSrc, iLev, ix, iy) = emisMap%arM(iSp, iSrc, iLev, ix, iy) * rd_dust_src%fzDisp(iLev)
                   mapCoordX%arm(iSp, iSrc, iLev, ix, iy) = 0.0
                   mapCoordY%arm(iSp, iSrc, iLev, ix, iy) = 0.0
                   
                   !emisMap%ifColumnValid(iSp, iSrc, iLev, ix, iy) = .true.
                   !emisMap%ifGridValid(iSp, iSrc) = .true.

                end do      ! iLev 
             end do      !  nx_dispersion
          end do     ! ny_dispersion
       end do    ! iSrc
    end do    ! iSp

  end subroutine compute_emission_for_rd_dust


  !********************************************************************************************

  logical function fu_rd_dust_emis_owned_quantity(rd_dust_src, quantity)
    !
    ! Check whether the specific quantity is handled by the emission sources exclusively
    ! Such quantities may be related to emission species somewhere else
    !
    implicit none

    ! Imported parameters
    type(silam_road_dust_source), intent(in) :: rd_dust_src
    integer, intent(in) :: quantity
    !
    ! The wind blown dust source does not have own quantities yet
    !
    select case(quantity)
      case default
        fu_rd_dust_emis_owned_quantity = .false.
    end select
  end function fu_rd_dust_emis_owned_quantity


  !*****************************************************************

  integer function fu_source_id_nbr_of_rd_dust_src(rd_dust_src)
    !
    ! Returns the source number. Reason: all sources are enumerated
    ! sequencially, so that the source can, in fact, be refered by its
    ! number without other information. But so far this number is 
    ! copied to the particles in the pollution cloud. 
    !
    ! NOTE. One and only index may be reasonable. The other one MUST be
    ! negative or zero
    !
    implicit none

    ! Imported parameters with intent IN
    type(silam_road_dust_source), intent(in) :: rd_dust_src

    ! Stupidity check
    if(.not.(rd_dust_src%defined == silja_false))then
      call set_error('Undefined source given','fu_source_id_nbr_of_rd_dust_src')
      return
    endif
    fu_source_id_nbr_of_rd_dust_src = rd_dust_src%id_nbr

  end function fu_source_id_nbr_of_rd_dust_src


  !*************************************************************************

  integer function fu_source_nbr_of_rd_dust_src(rd_dust_src)
    !
    ! Returns the source number. Reason: all sources are enumerated
    ! sequencially, so that the source can, in fact, be refered by its
    ! number without other information. But so far this number is 
    ! copied to the particles in the pollution cloud. 
    !
    ! NOTE. One and only index may be reasonable. The other one MUST be
    ! negative or zero
    !
    implicit none

    ! Imported parameters with intent IN
    type(silam_road_dust_source), intent(in) :: rd_dust_src

    ! Stupidity check
    if(.not.(rd_dust_src%defined == silja_false))then
      fu_source_nbr_of_rd_dust_src = rd_dust_src%src_nbr
    else
      fu_source_nbr_of_rd_dust_src = int_missing
      call set_error('Undefined source given','fu_source_nbr_of_rd_dust_src')
      return
    endif

  end function fu_source_nbr_of_rd_dust_src


  !*************************************************************************

  subroutine typical_species_cnc_rd_dust(srcRDDust, species, nSpecies, arConc)
    !
    ! Guesses a typical level of concentration and divides it with the given accuracy factor
    !
    implicit none

    ! Imported parameters
    type(silam_road_dust_source), intent(in) :: srcRDDust
    type(silam_species), dimension(:), pointer :: species
    integer, intent(out) :: nSpecies
    real, dimension(:), pointer :: arConc

    ! Local variables
    integer :: iSpecies
    real :: fTotalFluxNbr, fTotalFluxVol, fMassMeanDiam
    type(Taerosol_mode) :: mode

    real, parameter :: fTypicalMassCnc = 1.0e-9  ! one microgram

    species => srcRDDust%species
    nSpecies = srcRDDust%nSpecies

    do iSpecies = 1, srcRDDust%nSpecies
      arConc(iSpecies) = fTypicalMassCnc
    end do

  end subroutine typical_species_cnc_rd_dust


  !*************************************************************************

  function fu_rd_dust_source_name(rd_dust_src)result(chNm)
    implicit none
    type(silam_road_dust_source), intent(in) :: rd_dust_src
    character(len=clen) :: chNm
    chNm = rd_dust_src%src_nm
  end function fu_rd_dust_source_name


  !*************************************************************************

  subroutine update_road_wetness_and_snow(met_buf, disp_buf, &
       & indprecrate, indrh2m, indweqsnowd, &
       & indSnowlessDays, indwetness, &
       & aver_interval, timestep_sec, mdl_now, mdl_timestep)
      !
      ! Updates the road wetness and days since last snow
    !
    
    implicit none

    ! Imported parameters
    integer, intent(in) :: indprecrate, indrh2m, indweqsnowd
    integer, intent(in) :: indSnowlessDays, indwetness
    type(Tfield_buffer), pointer :: met_buf, disp_buf
    type(silja_interval), intent(in) :: aver_interval, mdl_timestep
    type(real), intent(in) :: timestep_sec
    type(silja_time), intent(in) :: mdl_now
    
    ! Local variables
    integer :: ixDsp, iyDsp, ixMet, iyMet, iDayInYear, iDisp, iMet
    type(Tfield_buffer), pointer :: mb, db
    real :: fRH2m, fprecrate, fsnowdepth, seconds, day_fraction
    !v1,2,3
    !real :: snow_depth_limit = 0.00035, wetrate = 3000./3600., dryrate = 0.17/3600.
    !v4
    real :: snow_depth_limit = 0.00035, wetrate = 6000./3600., dryrate = 0.67/3600.

    !mb => met_buf
    !db => disp_buf
    day_fraction = timestep_sec/(24*3600.)

    !call msg('day_fraction', day_fraction)
      
    do iyDsp = 1, ny_dispersion
       do ixDsp = 1, nx_dispersion
          iDisp = ixDsp + (iyDsp-1) * nx_dispersion
          
          fRH2m = met_buf%p2d(indrh2m)%present%ptr(iDisp)
          fprecrate = met_buf%p2d(indprecrate)%present%ptr(iDisp)
          fsnowdepth = met_buf%p2d(indweqsnowd)%present%ptr(iDisp)

          !!! temporary
          if (disp_buf%p2d(indSnowlessDays)%past%ptr(iDisp) == real_missing) disp_buf%p2d(indSnowlessDays)%past%ptr(iDisp) = 1000
          if (disp_buf%p2d(indwetness)%past%ptr(iDisp) == real_missing) disp_buf%p2d(indwetness)%past%ptr(iDisp) = 1.0
          !!!!!!!!!!!!!!!
          
          if (fsnowdepth < snow_depth_limit) then
             disp_buf%p2d(indSnowlessDays)%present%ptr(iDisp) = &
                  & disp_buf%p2d(indSnowlessDays)%past%ptr(iDisp) + day_fraction
          else
             disp_buf%p2d(indSnowlessDays)%present%ptr(iDisp) = 0.0
          end if

          !v1,2,3
          !disp_buf%p2d(indwetness)%present%ptr(iDisp) = &
          !     & disp_buf%p2d(indwetness)%past%ptr(iDisp) + &
          !     & (wetrate * fprecrate - dryrate * (1-fRH2m))*timestep_sec

          !v4,5,6,7
          disp_buf%p2d(indwetness)%present%ptr(iDisp) = &
               & (disp_buf%p2d(indwetness)%past%ptr(iDisp) + &
               & (wetrate * fprecrate * timestep_sec)) * exp(-dryrate*timestep_sec*(1-fRH2m))

          !call msg('wetness before limit', disp_buf%p2d(indwetness)%present%ptr(iDisp))
          !call msg('exp(-dryrate*timestep_sec*(1-fRH2m))', exp(-dryrate*timestep_sec*(1-fRH2m)))

          !disp_buf%p2d(indwetness)%present%ptr(iDisp) = &
          !     & disp_buf%p2d(indwetness)%past%ptr(iDisp) + &
          !     & (wetrate * fprecrate * timestep_sec) * ( 1.0  - dryrate * timestep_sec * (1.0 - fRH2m))

          if (disp_buf%p2d(indwetness)%present%ptr(iDisp) > 1.0) disp_buf%p2d(indwetness)%present%ptr(iDisp) = 1.0
          if (disp_buf%p2d(indwetness)%present%ptr(iDisp) < 0.0) disp_buf%p2d(indwetness)%present%ptr(iDisp) = 0.0

          !call msg('disp_buf%p2d(indwetness)%past%ptr(iDisp), disp_buf%p2d(indwetness)%present%ptr(iDisp), fprecrate, fRH2m', &
          !     & (/ disp_buf%p2d(indwetness)%past%ptr(iDisp), disp_buf%p2d(indwetness)%present%ptr(iDisp), fprecrate, fRH2m /))
            
       end do   ! ix
    end do   ! iy

    call set_valid_time(disp_buf%p2d(indSnowlessDays)%present%idPtr, mdl_now)
    call set_valid_time(disp_buf%p2d(indwetness)%present%idPtr, mdl_now)
    !call set_validity_length(disp_buf%p2d(indSnowlessDays)%present%idPtr, mdl_timestep)
    !call set_validity_length(disp_buf%p2d(indwetness)%present%idPtr, mdl_timestep)

  end subroutine update_road_wetness_and_snow
  

  subroutine report_rd_dust_src(rd_dust_src)
    implicit none
    type(silam_road_dust_source), intent(in) :: rd_dust_src
    integer :: iSpecies
    call msg('Road dust source'+rd_dust_src%src_nm)
    do iSpecies = 1, rd_dust_src%nSpecies
      call report(rd_dust_src%species(iSpecies))
   end do
  
  end subroutine report_rd_dust_src

  
END MODULE source_terms_road_dust

