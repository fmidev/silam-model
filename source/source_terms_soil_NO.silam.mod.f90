MODULE source_terms_soil_NO
  !
  ! This module contains description of the soil NO emission.
  !
  ! All units: SI, unless otherwise stated
  !
  ! Code owner
  !
  ! Language: ANSI FORTRAN-90 (or close to it)
  !
  use source_terms_time_params !cocktail_basic

  implicit none
  private

  !
  ! PUBLIC routines of soil NO
  !
  public fill_soil_NO_src_from_namelist
  public reserve_soil_NO_source
  public init_emission_soil_NO
  public add_source_species_soil_NO_src
  public add_input_needs
  public link_source_to_species
  public source_2_second_grid
  public compute_emission_for_soil_NO
  public fu_soil_NO_emis_owned_quantity
  public fu_name
  public fu_source_nbr
  public fu_source_id_nbr
  public report
  public create_source_containing_grid

  !
  ! Private routines of the soil NO source
  !
  private create_src_cont_grd_soilno_src
  private add_input_needs_soil_NO_src
  private link_soil_NO_src_to_species
  private fu_source_id_nbr_of_soil_NO_src
  private fu_source_nbr_of_soil_NO_src
  private fu_soil_NO_source_name
  private report_soil_NO_src

  !
  ! Private subs of soil NO source
  !
  interface create_source_containing_grid
    module procedure create_src_cont_grd_soilno_src
  end interface
  
  interface add_input_needs
    module procedure add_input_needs_soil_NO_src
  end interface

  interface link_source_to_species
    module procedure link_soil_NO_src_to_species
  end interface

  interface source_2_second_grid
    module procedure project_soil_NO_src_second_grd
  end interface

  interface fu_source_nbr
    module procedure fu_source_nbr_of_soil_NO_src
  end interface

  interface fu_source_id_nbr
    module procedure fu_source_id_nbr_of_soil_NO_src
  end interface

  interface fu_name
    module procedure fu_soil_NO_source_name
  end interface

  interface report
    module procedure report_soil_NO_src
  end interface

  ! The soil NO source term
  !
  TYPE silam_soil_NO_source
    PRIVATE
    CHARACTER(len=clen) :: src_nm, sector_nm       ! Name of the area source and sector
    character(len=fnlen) :: src_data_dir
    integer :: src_nbr, id_nbr  ! A source and id numbers in a WHOLE source list
    integer :: nLevsDispVert, nSpecies
    type(silam_vertical) :: vertLevsDispVert
    real, dimension(:), pointer :: levFractDispVert, fzDisp
    type(Tsilam_namelist), pointer :: nlInputFiles  ! namelist for names of supplementary files
    type(silam_species), dimension(:), pointer :: species
    type(silja_field), pointer :: source_mask
    type(chemical_adaptor) :: adaptor
    type(silja_logical) :: defined
  END TYPE silam_soil_NO_source

  type soil_NO_src_ptr
    type(silam_soil_NO_source) :: soilNO_src
  end type soil_NO_src_ptr
  public soil_NO_src_ptr


CONTAINS


  !*********************************************************************

  subroutine fill_soil_NO_src_from_namelist(nlSetup, srcsoilNO, expected_species, src_data_dir)
    !
    ! Initializes the soil NO source term.
    !
    implicit none

    ! Imported parameters
    type(Tsilam_namelist), pointer :: nlSetup
    type(silam_soil_NO_source), intent(inout) :: srcsoilNO
    type(silam_species), dimension(:), intent(in), allocatable :: expected_species
    character(len=*), intent(in) :: src_data_dir

    ! Local variables
    integer :: iTmp, jTmp, iSpecies, nFiles, stat
    type(Tsilam_nl_item_ptr), dimension(:), pointer ::  pItems
    type(silam_species) :: species
    type(silam_material), pointer :: material
    integer, dimension(:), pointer :: indices
    character(len=fnlen) :: chTmp
    logical :: ifFound
    character(len=*), parameter :: subname = 'fill_soil_NO_src_from_namelist'

    call msg('dbg filling soil no source from namelist')
    !
    ! Names
    !
    srcsoilNO%src_nm = fu_content(nlSetup,'source_name')
    srcsoilNO%sector_nm = fu_content(nlSetup,'source_sector_name')
    srcsoilNO%src_data_dir = src_data_dir
    srcsoilNO%defined = silja_false

    !
    ! Store the input fields for the source features.
    !
    srcsoilNO%nlInputFiles => fu_create_namelist('soil_NO_src_supplementary_files')
    if(error)return

    pItems => null()
    call get_items(nlSetup, 'supplementary_file', pItems, nFiles)
    if(error)return
    !
    ! No grads hat expansion or alike: the string contains the file format
    !
    do iTmp = 1, nFiles
      chTmp = fu_process_filepath(fu_content(pItems(iTmp)),superdir=srcsoilNO%src_data_dir)
      call add_namelist_item(srcsoilNO%nlInputFiles, 'supplementary_file', chTmp)
    end do
    
    material => fu_get_material_ptr(fu_content(nlSetup, 'emitted_substance'))

    srcsoilNO%nSpecies = 1
    allocate(srcsoilNO%species(1), stat=stat)
    if (fu_fails(stat == 0, 'Allocate failed', subname)) return
    
    call set_species(srcsoilNO%species(1), material, in_gas_phase)
    
    if (fu_fails(fu_content(nlSetup, 'emitted_substance') /= '', 'Missing emitted_substance', subname)) return

    if(error)return

    srcsoilNO%defined = silja_undefined

  end subroutine fill_soil_NO_src_from_namelist


  !*****************************************************************

  subroutine add_input_needs_soil_NO_src(soil_NO_src, q_met_dynamic, q_met_st, &
                                                  & q_disp_dynamic, q_disp_st, wdr)
    !
    ! Returns input needs for the emission/transformation routines. 
    !
    implicit none

    ! Imported parameters
    type(silam_soil_NO_source), intent(in) :: soil_NO_src
    integer, dimension(:), intent(inout) :: q_met_dynamic, q_met_st, &
                                          & q_disp_dynamic, q_disp_st
    type(silja_wdr), intent(in), optional :: wdr

    ! Local variables
    integer :: iTmp
    
    iTmp = fu_merge_integer_to_array(water_eq_snow_depth_flag, q_met_dynamic)
    iTmp = fu_merge_integer_to_array(relative_humidity_2m_flag, q_met_dynamic)
    iTmp = fu_merge_integer_to_array(temperature_2m_flag, q_met_dynamic)

    if (present(wdr)) then ! No wdr -- no meteo requests                                                                                                            
       iTmp = fu_merge_integer_to_array(soil_moisture_vol_frac_nwp_flag, q_met_dynamic) !soil_moisture_content_flag, q_met_dynamic)
       
       if ( wdr == wdr_missing) then  ! Can be undefined  
          iTmp = fu_merge_integer_to_array(leaf_area_index_flag, q_met_dynamic)
       elseif (any (fu_LAIsrc(wdr) == (/LAI_dynamic_1, LAI_dynamic_2/))) then
          iTmp = fu_merge_integer_to_array(leaf_area_index_flag, q_met_dynamic)
       elseif(any (fu_LAIsrc(wdr) == (/LAI_static_1, LAI_static_2/)))then
          iTmp = fu_merge_integer_to_array(leaf_area_index_flag, q_met_st)
       else
          call set_error('Soil NO requires leaf area index', 'input_needs_soil_NO_source')
          return
       endif
    endif
    iTmp = fu_merge_integer_to_array(soil_NO_emis_0_flag, q_disp_st)
    
    ! Finally, fraction of land must always be in the permanent input
    !
    iTmp = fu_merge_integer_to_array(fraction_of_land_flag, q_met_st)

  end subroutine add_input_needs_soil_NO_src


  !**************************************************************************

  subroutine reserve_soil_NO_source(soil_NO_src, &     ! Src to initialise
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
    type(silam_soil_NO_source), intent(inout) :: soil_NO_src
    integer, intent(in) :: iSrcNbr, iSrcIdNbr

    !
    ! Nullify the basic variables
    !
    soil_NO_src%src_nm = ''
    soil_NO_src%sector_nm = ''
    !
    ! A bit of other stuff
    !
    nullify(soil_NO_src%fZDisp)
    nullify(soil_NO_src%levFractDispVert)
    !
    ! Main source parameters - enough to identify it in the global information list
    !
    soil_NO_src%src_nbr = iSrcNbr
    soil_NO_src%id_nbr = iSrcIdNbr
    !
    ! Finally, mark the source as incomplete
    !
    soil_NO_src%defined = silja_false

  end subroutine reserve_soil_NO_source


  !*********************************************************************

  subroutine init_emission_soil_NO(srcsoilNO, dispersionMarketPtr, start_time)
    !
    ! Initializes the soil NO source term.
    !
    implicit none

    ! Imported parameters
    type(Tsilam_namelist), pointer :: nlSetup
    type(silam_soil_NO_source), intent(inout) :: srcsoilNO
    type(mini_market_of_stacks), pointer :: dispersionMarketPtr
    type(silja_time), intent(in) :: start_time

    ! Local variables
    integer :: iTmp, iSpecies, iFlds, iZ, iWind, iThresh
    type(silam_sp) :: sp
    real :: fTmp
    type(silja_field_id) :: id, idTmp
    real, dimension(:), pointer :: arPtr
    type(silja_shopping_list) :: shop_list
    integer, dimension(:), pointer ::  q_disp_dyn, q_disp_stat, stack_quantities
    type(silam_vertical) :: vertTmp
    logical :: ifOK
    type(silja_field), pointer :: fieldPtr

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
    call add_input_needs_soil_NO_src(srcsoilNO, &
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
      iSpecies = int_missing
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
                                 & srcsoilNO%species(iSpecies), &
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
    call msg('Filling-in the soil NO emission info from supplementary fields')
    call fill_minimarket_from_namelist(dispersionMarketPtr, &
                                     & srcsoilNO%nlInputFiles, 'supplementary_file', & ! namelist and item
                                     & shop_list, start_time, &
                                     & static_climatology, &  ! target stack
                                     & create_field, &            ! error if a clash
                                     & wdr_missing, &
                                     & dispersion_gridPtr, &
                                     & 5, .true., & ! iAccuracy, ifAdjustGrid
                                     & ifOK)

    ! Check that all the requested quantities are in the stack.
    ! Some of them might be initialized but if not, can still be added.
    !
    id = fu_set_field_id_simple(met_src_missing,&
         & soil_NO_emis_0_flag, &
         & time_missing, &        ! valid time
         & level_missing)
        fieldPtr => fu_get_field_from_mm_general(dispersionMarketPtr, id, .false.)
        if(.not. associated(fieldPtr))then
          call set_error('Soil NO emission map is missing but needed','init_emission_soil_NO')
          return
        endif

    allocate(srcsoilNO%species(1))
        
    call set_species(srcsoilNO%species(1), fu_get_material_ptr('NO'), in_gas_phase)
        
    srcsoilNO%defined = silja_true

    call free_work_array(sp%sp)
    call free_work_array(stack_quantities)
    call free_work_array(q_disp_dyn)
    call free_work_array(q_disp_stat)

    call report(srcsoilNO)

    !=======================================================================

  end subroutine init_emission_soil_NO


  !****************************************************************************

  subroutine add_source_species_soil_NO_src(soil_NO_src, species_list, nSpecies)
    !
    ! Get the species emitted by this source and make a list of
    ! them. The source must be initialized, of course.
    !
    implicit none
    type(silam_soil_NO_source), intent(in) :: soil_NO_src
    type(silam_species), dimension(:), pointer :: species_list
    integer, intent(inout) :: nSpecies

    call addSpecies(species_list, nSpecies, soil_NO_src%species, soil_NO_src%nSpecies)

  end subroutine add_source_species_soil_NO_src


  !**************************************************************************

  subroutine link_soil_NO_src_to_species(species_list, soil_NO_src)
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
    type(silam_soil_NO_source), intent(inout) :: soil_NO_src
    type(silam_species), dimension(:), pointer :: species_list

    !
    ! Linkage is actually just creation of the chemical adaptor
    !
    call create_adaptor(soil_NO_src%species, species_list, soil_NO_src%adaptor)
    
  end subroutine link_soil_NO_src_to_species

  !*****************************************************************

  subroutine project_soil_NO_src_second_grd(soil_NO_src, grid, vert_disp, vert_proj, iAccuracy)
    !
    ! There is nothing to project in terms of horizontal grid but the vertical has to be made.
    !
    implicit none

    ! Imported parameters
    type(silam_soil_NO_source), intent(inout) :: soil_NO_src
    type(silja_grid), intent(in) :: grid
    type(silam_vertical), intent(in) :: vert_disp, vert_proj
    integer, intent(in) :: iAccuracy

    ! Local variables
    integer :: i, ilev
    type(silam_vertical) :: vertTmp
    real, dimension(:), pointer :: arTmp, dz_disp

    
    if(len_trim(soil_NO_src%sector_nm) > 0)then
      call msg('Re-projecting soil NO source:' + soil_NO_src%src_nm +'_' + soil_NO_src%sector_nm)
    else
      call msg('Re-projecting soil NO source:' + soil_NO_src%src_nm)
    endif

    !
    ! Now re-project the vertical grid of the source to the given vertical
    ! Since vertical distributions are very poorly known and crude, we do not need
    ! any precise meteo-dependent projection, crude will do the job
    !
    soil_NO_src%vertLevsDispVert = vert_disp
    allocate(soil_NO_src%levFractDispVert(fu_NbrOfLevels(vert_disp)), &
           & soil_NO_src%fzDisp(fu_NbrOfLevels(vert_disp)), stat=i)
    if(fu_fails(i==0, 'Failed to allocate dispersion-vertical level fractions','project_soil_NO_src_second_grd'))return
    soil_NO_src%levFractDispVert(:) = 0.0
    soil_NO_src%fzDisp(:) = 0.0

    arTmp => fu_work_array()
    if(error)return
    !
    ! Create the soil NO vertical
    !
    
    ! All emission between 0 m and 5 m
    call set_vertical(fu_set_layer_between_two(layer_btw_2_height, 0.0, 5.0), vertTmp)
    if(error)return
    arTmp(1) = 1.0

    call reproject_verticals(vertTmp, arTmp, &                   ! vertical from, fractions from
                           & vert_proj, soil_NO_src%levFractDispVert, &   ! vertical to, fractions to
                           ! mass centres, number of non-zero levels:
                           & soil_NO_src%fzDisp, soil_NO_src%nLevsDispVert, &
                           & ifMassCentreInRelUnit=.true.)

    call set_missing(vertTmp, .false.)
    call free_work_array(arTmp)

  end subroutine project_soil_NO_src_second_grd


  !**************************************************************************

  subroutine compute_emission_for_soil_NO(soil_NO_src, &
                                        & met_buf, disp_buf, & 
                                        & now, &      ! current time
                                        & timestep, & ! model time step
                                        & pHorizInterpMet2DispStruct, ifHorizInterp, &
                                        & ifSpeciesMoment, &
                                        & emisMap, mapCoordX, mapCoordY, mapCoordZ, & ! Output
                                        & fMassInjected)                              ! Output
    !
    ! This routine is to be called at each model time step but not inside the main cycle,
    ! therefore its efficiency is of moderate importance.
    !
    implicit none

    ! Imported parameters
    type(silam_soil_NO_source), target, intent(in) :: soil_NO_src
    type(Tfield_buffer), pointer ::  met_buf, disp_buf  ! meteo and internal field buffers
    type(silja_time), intent(in) :: now           ! current time
    type(silja_interval), intent(in) :: timestep  ! model time step
    type(THorizInterpStruct), intent(in) :: pHorizInterpMet2DispStruct
    logical, intent(in) :: ifHorizInterp, ifSpeciesMoment
    type(Tmass_map), intent(inout) :: emisMap, mapCoordX, mapCoordY, mapCoordZ
    real(r8k), dimension(:), intent(inout) :: fMassInjected


    ! Local variables
    integer :: iLev, iThresh, iZ, iSoilNO, &
             & iMeteo, iDisp, ix, iy, ixMeteo, iyMeteo, iStat, iTmp, iMeteoTmp
    real :: fTmp, timestep_sec, fCellTotal, fSoilMassWater, &
          & fFluxTotal, fSoilMassWaterThresh, &
          & soil_moisture, snow_depth, lai, temp_2m, relhum_2m
    real, parameter :: a = 1.0, b = 5.0
    real, dimension(:), pointer :: ptrXSzDisp, ptrYSzDisp, pSoilVolumeWater, &
                                 & fMinDArray, fMaxDArray, ptemp2m, &
                                 & prelhum2m, pLAI, pLandFr, pNOEmis0, &
                                 & pSnowWEQDepth, ptrLonMeteo, ptrLatMeteo
    real*4, dimension(worksize) :: arTmp
    type(silja_field), pointer :: fldMaskPtr
    integer, save :: iCount = 0, iCountTalking = 0, uDump, iRec
    logical, save :: ifFirst=.true., ifTalking=.true.

    integer, dimension(10), save :: iArStatistics=0

    type(silam_sp) :: sp

    ! Local parameters 

    real, parameter :: fSnowWEQDepthThreshold = 0.005   !   0.5 cm of water-eq. snow => no soil NO
     
    !
    ! Meteo quantities
    !
    if(fu_fails(fu_index(met_buf, soil_moisture_vol_frac_nwp_flag, pSoilVolumeWater) /= int_missing, &  ! soil water content, m3/m3
    & 'Failed soil moisture input','compute_emission_for_soil_NO'))return
    
    if(fu_fails(fu_index(met_buf, water_eq_snow_depth_flag, pSnowWEQDepth) /= int_missing, &    ! snow water equiv. depth
                                & 'Failed water_eq_snow_depth_flag','compute_emission_for_soil_NO'))return
    if(fu_fails(fu_index(met_buf, leaf_area_index_flag, pLAI) /= int_missing, &                 ! leaf area index
                                & 'Failed leaf area index','compute_emission_for_soil_NO'))return
    if(fu_fails(fu_index(met_buf, fraction_of_land_flag, pLandFr) /= int_missing, &             ! fraction of land                
                                & 'Failed fraction_of_land_flag','compute_emission_for_soil_NO'))return
    if(fu_fails(fu_index(met_buf, temperature_2m_flag, ptemp2m) /= int_missing, &                 
                                & 'Failed 2m temperature index','compute_emission_for_soil_NO'))return
    if(fu_fails(fu_index(met_buf, relative_humidity_2m_flag, prelhum2m) /= int_missing, &                
                                & 'Failed 2m relative humidity index','compute_emission_for_soil_NO'))return
    if(fu_fails(fu_index(disp_buf, soil_NO_emis_0_flag, pNOEmis0) /= int_missing, &    ! basic map for soil NO emission        
             & 'Failed soil_NO_emis_0_flag','compute_emission_for_soil_NO'))return
    
    ptrXSzDisp => fu_grid_data(dispersion_cell_x_size_fld)  ! get grid cell size in metres
    ptrYSzDisp => fu_grid_data(dispersion_cell_y_size_fld)

    ptrLonMeteo => fu_grid_data(meteo_longitude_fld)
    ptrLatMeteo => fu_grid_data(meteo_latitude_fld)

    timestep_sec = abs(fu_sec(timestep))
    if(error)return

    !
    ! Now scan the domain, for each grid cell find the dust emission parameters
    ! and fill-in the emission map
    !
    do iy = 1, ny_dispersion
      do ix = 1, nx_dispersion
        iDisp = ix + (iy-1) * nx_dispersion

        iMeteo =  fu_grid_index(nx_meteo, ix, iy, pHorizInterpMet2DispStruct)
        iyMeteo = int(iMeteo / nx_meteo)
        ixMeteo = iMeteo - (iyMeteo-1)*nx_meteo
        !
        !   ... and separately if meteo-land mask is almost zero
        !
        if(pLandFr(iMeteo) < 1.0e-5)cycle

        soil_moisture = pSoilVolumeWater(iMeteo)
        lai = pLAI(iMeteo)
        snow_depth = pSnowWEQDepth(iMeteo)
        temp_2m = ptemp2m(iMeteo)
        relhum_2m = prelhum2m(iMeteo)

        if(pSnowWEQDepth(iMeteo) > fSnowWEQDepthThreshold)cycle

        temp_2m = min(max(temp_2m-273.0, 0.0), 30.0)

        fFluxTotal = (0.35/0.03) * pNOEmis0(iDisp) * timestep_sec * exp(0.103*(temp_2m)) * a *soil_moisture * exp(-b * soil_moisture**2)
        fFluxTotal = fFluxTotal * (max(0.5, relhum_2m) + 0.5)

        ! Fill-in the emission map
        !
        do iLev = 1, soil_NO_src%nLevsDispVert
          !
          ! First do the check for the overlap: speed-up
          !
          if(abs(soil_NO_src%levFractDispVert(iLev)) < 1e-15)cycle  ! nothing for this dispersion layer

          fCellTotal = 0.0

          do iSoilNO = 1, soil_NO_src%nSpecies
            fTmp = fFluxTotal * &                        ! Total flux
                   & soil_NO_src%levFractDispVert(iLev) * &    ! Vertical overlap
                   & ptrXSzDisp(iDisp) * ptrYSzDisp(iDisp)     ! m2 size

            fCellTotal = fCellTotal + fTmp

            emisMap%arM(soil_NO_src%adaptor%iSp(iSoilNO),soil_NO_src%id_nbr,iLev,ix,iy) = fTmp + &
                       & emisMap%arM(soil_NO_src%adaptor%iSp(iSoilNO),soil_NO_src%id_nbr,iLev,ix,iy)
            fMassInjected(soil_NO_src%adaptor%iSp(iSoilNO)) = &
                                        & fMassInjected(soil_NO_src%adaptor%iSp(iSoilNO)) + fTmp

          end do ! nSpecies
          emisMap%ifColumnValid(soil_NO_src%id_nbr,ix,iy) = .true.
          emisMap%ifGridValid(iLev,soil_NO_src%id_nbr) = .true.

          if (ifSpeciesMoment) then
             mapCoordZ%arm(soil_NO_src%adaptor%iSp(iSoilNO), soil_NO_src%id_nbr, ilev, ix,iy) = &
                   mapCoordZ%arm(soil_NO_src%adaptor%iSp(iSoilNO), soil_NO_src%id_nbr, ilev, ix,iy) + &
                   & fTmp * soil_NO_src%fzDisp(iLev)
          end if
         
          if (.not. ifSpeciesMoment) then
             mapCoordZ%arM(1,soil_NO_src%id_nbr, iLev, ix, iy) = &
                           & soil_NO_src%fzDisp(iLev) * fCellTotal + &
                           & mapCoordZ%arM(1,soil_NO_src%id_nbr, iLev, ix, iy)
          end if
        end do      ! iLev
      end do     ! nx_dispersion
    end do    ! ny_dispersion

  end subroutine compute_emission_for_soil_NO


  !********************************************************************************************

  logical function fu_soil_NO_emis_owned_quantity(soil_NO_src, quantity)
    !
    ! Check whether the specific quantity is handled by the emission sources exclusively
    ! Such quantities may be related to emission species somewhere else
    !
    implicit none

    ! Imported parameters
    type(silam_soil_NO_source), intent(in) :: soil_NO_src
    integer, intent(in) :: quantity
    !
    ! The wind blown dust source does not have own quantities yet
    !
    select case(quantity)
      case default
        fu_soil_NO_emis_owned_quantity = .false.
    end select
  end function fu_soil_NO_emis_owned_quantity


  !*****************************************************************

  integer function fu_source_id_nbr_of_soil_NO_src(soil_NO_src)
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
    type(silam_soil_NO_source), intent(in) :: soil_NO_src

    ! Stupidity check
    if(.not.(soil_NO_src%defined == silja_false))then
      call set_error('Undefined source given','fu_source_id_nbr_of_soil_NO_src')
      return
    endif
    fu_source_id_nbr_of_soil_NO_src = soil_NO_src%id_nbr

  end function fu_source_id_nbr_of_soil_NO_src


  !*************************************************************************

  subroutine create_src_cont_grd_soilno_src(soil_no_src, grid_template, ifVerbose, ifMinimal, ifExtended)
    !
    ! Creates the grid that covers the area with active soil NO emission
    !
    implicit none

    ! Imported parameters
    type(silam_soil_NO_source), intent(in) :: soil_no_src
    type(silja_grid), intent(inout) :: grid_template
    logical, intent(in) :: ifVerbose, ifMinimal
    logical, intent(out) :: ifExtended

    ! So far nothing to do: this source just covers the dispersion grid
    !
    ifExtended = .false.
    return

  end subroutine create_src_cont_grd_soilno_src

  integer function fu_source_nbr_of_soil_NO_src(soil_NO_src)
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
    type(silam_soil_NO_source), intent(in) :: soil_NO_src

    ! Stupidity check
    if(.not.(soil_NO_src%defined == silja_false))then
      fu_source_nbr_of_soil_NO_src = soil_NO_src%src_nbr
    else
      fu_source_nbr_of_soil_NO_src = int_missing
      call set_error('Undefined source given','fu_source_nbr_of_soil_NO_src')
      return
    endif

  end function fu_source_nbr_of_soil_NO_src

  !*************************************************************************

  function fu_soil_NO_source_name(soil_NO_src)result(chNm)
    implicit none
    type(silam_soil_NO_source), intent(in) :: soil_NO_src
    character(len=clen) :: chNm
    chNm = soil_NO_src%src_nm
  end function fu_soil_NO_source_name


  !*************************************************************************

  subroutine report_soil_NO_src(soil_NO_src)
    implicit none
    type(silam_soil_NO_source), intent(in) :: soil_NO_src
    integer :: iSpecies
    call msg('Soil NO source'+soil_NO_src%src_nm)
    do iSpecies = 1, soil_NO_src%nSpecies
      call report(soil_NO_src%species(iSpecies))
    end do
    call msg('')
    call msg('Reporting the soil NO source settings:')
  
  end subroutine report_soil_NO_src

  
END MODULE source_terms_soil_NO

