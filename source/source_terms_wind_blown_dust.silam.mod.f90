MODULE source_terms_wind_blown_dust
  !
  ! This module contains description of the wind-blown dust emission.
  !
  ! Emission is computed for the prescribed size classes onto which the known spectrum is
  ! projected. It needs wind speed and soil parameters described as a series of maps in the
  ! source description file.
  !
  ! Since the dependencies of flux on particle diameter may be tricky, all or most 
  ! of integrals are to be taken numerically. To save this pain during runtime, we 
  ! pre-compute them at the very beginning and will later just pick the values from the array 
  ! defined here.
  !
  ! All units: SI, unless otherwise stated
  !
  ! Code owner Mikhail Sofiev mikhail.sofiev@fmi.fi
  !
  ! Language: ANSI FORTRAN-90 (or close to it)
  !
  use source_terms_time_params !cocktail_basic

  implicit none
  private

  !
  ! PUBLIC routines of wind-blown dust
  !
  public fill_wb_dust_src_from_namelist
  public reserve_wb_dust_source
  public init_emission_wb_dust
  public add_inventory_wb_dust_src
  public create_source_containing_grid
  public add_input_needs
  public link_source_to_species
  public source_2_second_grid
  public compute_emission_for_wb_dust
  public fu_wb_dust_emis_owned_quantity
  public fu_name
  public fu_source_nbr
  public fu_source_id_nbr
  public typical_species_conc
  public report

  !
  ! Private routines of the sea salt source
  !
  private create_src_cont_grd_wbdust_src
  private add_input_needs_wb_dust_src
  private link_wb_dust_src_to_species
  private fu_species_index_src
  private wind_blown_dust_flux4mode
  private fu_source_id_nbr_of_wb_dust_src
  private fu_source_nbr_of_wb_dust_src
  private fu_wb_dust_source_name
  private typical_species_cnc_wb_dust
  private report_wb_dust_src

  !
  ! Private subs of sea salt source
  !
  interface create_source_containing_grid
    module procedure create_src_cont_grd_wbdust_src
  end interface

  interface add_input_needs
    module procedure add_input_needs_wb_dust_src
  end interface

  interface link_source_to_species
    module procedure link_wb_dust_src_to_species
  end interface

  interface source_2_second_grid
    module procedure project_wbdust_src_second_grd
  end interface

  interface fu_source_nbr
    module procedure fu_source_nbr_of_wb_dust_src
  end interface

  interface fu_source_id_nbr
    module procedure fu_source_id_nbr_of_wb_dust_src
  end interface

  interface fu_name
    module procedure fu_wb_dust_source_name
  end interface

  interface typical_species_conc
    module procedure typical_species_cnc_wb_dust
  end interface

  interface report
    module procedure report_wb_dust_src
  end interface

  !
  ! There might be several types of the emission algorithm....
  !
  integer, private, parameter :: emis_sandblasting_v1_flag = 6501
  integer, private, parameter :: emis_Gillette_DMAT = 6502            ! Not reallly implemented
  integer, private, parameter :: simple_dust_flag = 6503
  !
  ! ...and descriptions of the particle spectrum
  !
  integer, private, parameter :: lognorm_4_modes_dust_flag = 6510
  integer, private, parameter :: Kok_Brittle_fragm_theory_flag = 6511
  integer, private, parameter :: lognorm_4_modes_dust_numeric_flag = 6512

  !
  ! Local parameters for tuning the dust production
  !
  real, private, parameter :: fEntrainmentScaling = 0.2e-5, &  ! 2016.10.04.  Zender: 7e-4, EMEP: 2e-5.
!                            & fEntrainmentScaling = 0.3e-5, &  ! 2016.06.21.  Zender: 7e-4, EMEP: 2e-5.
!                            & fEntrainmentScaling = 0.25e-5, &  ! 2016.05.09.  Zender: 7e-4, EMEP: 2e-5.
                     & fSrfRoughSandOnly = 2.5e-4, &   ! 2016.11.02. Wider emission area may be worth trying
!                     & fSrfRoughSandOnly = 2e-4, &   ! 2016.08.24. Wider emission area is still needed
!                     & fSrfRoughSandOnly = 1e-4, &   ! 2016.06.27. Wider emission area is still needed
!                     & fSrfRoughSandOnly = 8e-5, &   ! 2016.06.27. Wider emission area is still needed
!                     & fSrfRoughSandOnly = 6e-5, &   ! 2016.06.27. Wider emission area is still needed
!                     & fSrfRoughSandOnly = 4e-5, &   ! 2016.06.21. Back to this nbr. Simultaneous up-scaling
!                     & fSrfRoughSandOnly = 3.e-5, &   ! 2016.06.03. Just 10% from 1.5e-5, need more
!                     & fSrfRoughSandOnly = 2.5e-5, &   ! 2016.05.24.att2 Just 10% from 1.5e-5, need more
!                     & fSrfRoughSandOnly = 2.0e-5, &   ! 2016.05.24. Need to widen the area
!                     & fSrfRoughSandOnly = 1.5e-5, &   ! 2016.05.13. Need to widen the area
!                     & fSrfRoughSandOnly = 1.2e-5, &   ! 2016.05.09. Simultansous with reduction of scaling
!                     & fSrfRoughSandOnly = 1.1e-5, &   ! 2016.04.17. Too narrow area but high emission inside
!                     & fSrfRoughSandOnly = 3e-5, &   ! still too high emission.
!                     & fSrfRoughSandOnly = 3.5e-5, &   ! too high emission over China, not sure about Africa. 
!                     & fSrfRoughSandOnly = 4e-5, &   ! seems to be too high emission everywhere
                     & fDragPartitioningSandOnly = 1.0/log(0.35*(0.1/fSrfRoughSandOnly)**0.8), &
                     & fLAI_critical = 1.0, &
!                     & u_star_excess_saturation = 0.1   ! m/s
!                     & u_star_excess_saturation = 0.09   ! m/s
!                     & u_star_excess_saturation = 0.06   ! m/s
!                     & u_star_excess_saturation = 0.008   ! m/s  2016.04.17
!                     & u_star_excess_saturation = 0.007   ! m/s  2016.04.19
!                     & u_star_excess_saturation = 0.006   ! m/s  2016.07.21
!                     & u_star_excess_saturation = 0.06   ! m/s  2016.07.21 - after putting saturation to gust
!                     & u_star_excess_saturation = 0.09   ! m/s  2016.07.21
!                     & u_star_excess_saturation = 0.06   ! m/s  2016.08.24
!                     & u_star_excess_saturation = 0.04   ! m/s  2016.09.16
!                     & u_star_excess_saturation = 0.03   ! m/s  2016.09.24
                     & u_star_excess_saturation = 0.02   ! m/s  2016.10.20
  !
  ! The wind-blown dust source term
  !
  TYPE silam_wind_blown_dust_source
    PRIVATE
    CHARACTER(len=clen) :: src_nm, sector_nm       ! Name of the area source and sector
    character(len=fnlen) :: src_data_dir
    integer :: emisMethod, iSpectrumType, src_nbr, id_nbr  ! A source and id numbers in a WHOLE source list
    real, dimension(:), pointer :: fluxPerModeNbr, fluxPerModeVol ! (nModes)
    integer :: nLevsDispVert, nSpecies
    type(silam_vertical) :: vertLevsDispVert
    real, dimension(:), pointer :: levFractDispVert, fzDisp
    type(Tsilam_namelist), pointer :: nlInputFiles  ! namelist for names of supplementary files
    type(silam_species), dimension(:), pointer :: species
    type(silja_field), pointer :: source_mask
    logical, dimension(:), pointer :: ifNumberFlux
    type(chemical_adaptor) :: adaptor
    type(Twind_gust_lookup) :: wgl
    type(silja_logical) :: defined
  END TYPE silam_wind_blown_dust_source

  type wb_dust_src_ptr
    type(silam_wind_blown_dust_source) :: wbdust_src
  end type wb_dust_src_ptr
  public wb_dust_src_ptr

! Do not define it here? please!!!
!!!!!!!!#define DEBUG_DUST 


CONTAINS


  !*********************************************************************

  subroutine fill_wb_dust_src_from_namelist(nlSetup, srcWBDust, expected_species, src_data_dir)
    !
    ! Initializes the sea salt source term.
    ! The parameters are read from the given ini file and stored to the returned 
    ! source term. After the stuff has been read, the main FluxPerMode is to be generated.
    ! This configuration allows for several sea salt sources
    !
    implicit none

    ! Imported parameters
    type(Tsilam_namelist), pointer :: nlSetup
    type(silam_wind_blown_dust_source), intent(inout) :: srcWBDust
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
    srcWBDust%src_nm = fu_content(nlSetup,'source_name')
    srcWBDust%sector_nm = fu_content(nlSetup,'source_sector_name')
    srcWBDust%src_data_dir = src_data_dir
    srcWBDust%defined = silja_false

    !
    ! Emission index type
    !
    select case(fu_str_u_case(fu_content(nlSetup,'wind_blown_dust_emission_method')))

      case ('GILLETTE_DMAT')
        srcWBDust%emisMethod = emis_Gillette_DMAT
        call set_error('GILLETTE_DMAT method does not work yet. use SANDBLASTING_V1','fill_wb_dust_src_from_namelist')
        return

      case ('SANDBLASTING_V1')
        srcWBDust%emisMethod = emis_sandblasting_v1_flag

      case ('SIMPLE_DUST')
        srcWBDust%emisMethod = simple_dust_flag

      case default
        call set_error('Unknown emission method:' + &
                     & fu_content(nlSetup,'wind_blown_dust_emission_method'), &
                     & 'fill_wb_dust_src_from_namelist')
        return
    end select

    !
    ! Spectrum of the particles
    !
    select case(fu_str_u_case(fu_content(nlSetup,'wind_blown_dust_spectrum')))

      case ('LOGNORMAL_FOUR_MODES')
        srcWBDust%iSpectrumType = lognorm_4_modes_dust_numeric_flag

      case ('KOK_BRITTLE_THEORY_SINGLE_MODE')
        srcWBDust%iSpectrumType = Kok_Brittle_fragm_theory_flag

      case default
        call set_error('Unknown spectrum type:' + fu_content(nlSetup,'wind_blown_dust_spectrum'), &
                     & 'fill_wb_dust_src_from_namelist')
        return
    end select

    !
    ! Now the list of aerosol modes that will be emitted. To set them up, we will
    ! use the standard aerosol procedure - for the sake of unification.
    ! Note that there are two potentially concurring definitions: one coming from aerosol dynamics
    ! the other - written in the ini file. The first one prevails, if it exists
    !
    call get_source_aer_species(srcWBDust%species, srcWBDust%nSpecies, &
                              & expected_species, fu_content(nlSetup,'wind_blown_dust_substance_name'), &
                              & nlSetup)
    if(error)return
    !
    ! If number-emission is neeed, it has to be requested via nbr_aer in the species 
    !
    allocate(srcWBDust%ifNumberFlux(srcWBDust%nSpecies), stat=iTmp)
    if(fu_fails(iTmp == 0, 'Failed to allocate number-mass switcher array','fill_wb_dust_src_from_namelist'))return

    do iTmp = 1, srcWBDust%nSpecies
      srcWBDust%ifNumberFlux(iTmp) = (fu_name(fu_material(srcWBDust%species(iTmp))) == 'nbr_aer')
    end do

    !
    ! Store the input fields for the source features.
    !
    srcWBDust%nlInputFiles => fu_create_namelist('wind_blown_dust_src_supplementary_files')
    if(error)return

    pItems => null()
    call get_items(nlSetup, 'supplementary_file', pItems, nFiles)
    if(error)return
    !
    ! No grads hat expansion or alike: the string contains the file format
    !
    do iTmp = 1, nFiles
      chTmp = fu_process_filepath(fu_content(pItems(iTmp)),superdir=srcWBDust%src_data_dir)
      call add_namelist_item(srcWBDust%nlInputFiles, 'supplementary_file', chTmp)
    end do
    !
    ! Note that the source mask cannot be read from the file yet: dispersion grid is undefined
    !
    chTmp = fu_process_filepath(fu_content(nlSetup,'source_area_mask'),superdir=srcWBDust%src_data_dir)
    call add_namelist_item(srcWBDust%nlInputFiles, 'source_area_mask', chTmp)

    if(error)return

    srcWBDust%defined = silja_undefined

  end subroutine fill_wb_dust_src_from_namelist


  !*****************************************************************

  subroutine add_input_needs_wb_dust_src(wb_dust_src, q_met_dynamic, q_met_st, &
                                                  & q_disp_dynamic, q_disp_st, wdr)
    !
    ! Returns input needs for the emission/transformation routines. 
    !
    implicit none

    ! Imported parameters
    type(silam_wind_blown_dust_source), intent(in) :: wb_dust_src
    integer, dimension(:), intent(inout) :: q_met_dynamic, q_met_st, &
                                          & q_disp_dynamic, q_disp_st
    type(silja_wdr), intent(in), optional :: wdr

    ! Local variables
    integer :: iTmp

    !
    ! Add needed dynamic quantities. Always precipitation, the rest depends on the method
    !
    iTmp = fu_merge_integer_to_array(total_precipitation_rate_flag,  q_met_st)
    iTmp = fu_merge_integer_to_array(water_eq_snow_depth_flag, q_met_dynamic)

    select case(wb_dust_src%emisMethod)
      case(emis_Gillette_DMAT)
        iTmp = fu_merge_integer_to_array(friction_velocity_flag, q_met_dynamic)

      case(emis_sandblasting_v1_flag)
        if (present(wdr)) then ! No wdr -- no meteo requests
          iTmp = fu_merge_integer_to_array(friction_velocity_flag, q_met_dynamic)
          iTmp = fu_merge_integer_to_array(convective_velocity_scale_flag, q_met_dynamic)
          iTmp = fu_merge_integer_to_array(soil_moisture_vol_frac_nwp_flag, q_met_dynamic) !soil_moisture_content_flag, q_met_dynamic)
          iTmp = fu_merge_integer_to_array(u_10m_flag, q_met_dynamic)
          iTmp = fu_merge_integer_to_array(v_10m_flag, q_met_dynamic)

          iTmp = fu_merge_integer_to_array(soil_sand_mass_fraction_flag, q_met_st)
          iTmp = fu_merge_integer_to_array(soil_clay_mass_fraction_flag, q_met_st)
          iTmp = fu_merge_integer_to_array(fraction_of_land_flag, q_met_st)
          if ( wdr == wdr_missing) then  ! Can be undefined
             iTmp = fu_merge_integer_to_array(leaf_area_index_flag, q_met_dynamic)
          elseif (any (fu_LAIsrc(wdr) == (/LAI_dynamic_1, LAI_dynamic_2/))) then
            iTmp = fu_merge_integer_to_array(leaf_area_index_flag, q_met_dynamic)
          elseif(any (fu_LAIsrc(wdr) == (/LAI_static_1, LAI_static_2/)))then
            iTmp = fu_merge_integer_to_array(leaf_area_index_flag, q_met_st)
          else
            call set_error('WB dust requires leaf area index', 'input_needs_wb_dust_source')
            return
          endif
          iTmp = fu_merge_integer_to_array(land_roughness_disp_flag, q_met_st)
       endif

        iTmp = fu_merge_integer_to_array(alluvial_sedim_index_flag, q_disp_st)

      case(simple_dust_flag)
          if (present(wdr)) then ! No wdr -- no meteo requests                                                                                                            
            iTmp = fu_merge_integer_to_array(soil_moisture_vol_frac_nwp_flag, q_met_dynamic) !soil_moisture_content_flag, q_met_dynamic)
            iTmp = fu_merge_integer_to_array(u_10m_flag, q_met_dynamic)
            iTmp = fu_merge_integer_to_array(v_10m_flag, q_met_dynamic)
            iTmp = fu_merge_integer_to_array(fraction_of_land_flag, q_met_st)
            if ( wdr == wdr_missing) then  ! Can be undefined  
              iTmp = fu_merge_integer_to_array(leaf_area_index_flag, q_met_dynamic)
            elseif (any (fu_LAIsrc(wdr) == (/LAI_dynamic_1, LAI_dynamic_2/))) then
              iTmp = fu_merge_integer_to_array(leaf_area_index_flag, q_met_dynamic)
            elseif(any (fu_LAIsrc(wdr) == (/LAI_static_1, LAI_static_2/)))then
              iTmp = fu_merge_integer_to_array(leaf_area_index_flag, q_met_st)
            else
              call set_error('WB dust requires leaf area index', 'input_needs_wb_dust_source')
              return
            endif
          endif
          iTmp = fu_merge_integer_to_array(dust_emis_0_flag, q_disp_st)
          iTmp = fu_merge_integer_to_array(u_flag, q_met_dynamic)
          iTmp = fu_merge_integer_to_array(v_flag, q_met_dynamic)
      case default
        call set_error('Unknown emission computation algorithm', 'input_needs_wb_dust_source')
        return
    end select

    ! Finally, fraction of land must always be in the permanent input
    !
    iTmp = fu_merge_integer_to_array(fraction_of_land_flag, q_met_st)

  end subroutine add_input_needs_wb_dust_src


  !**************************************************************************

  subroutine reserve_wb_dust_source(wb_dust_src, &     ! Src to initialise
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
    type(silam_wind_blown_dust_source), intent(inout) :: wb_dust_src
    integer, intent(in) :: iSrcNbr, iSrcIdNbr

    !
    ! Nullify the basic variables
    !
    wb_dust_src%src_nm = ''
    wb_dust_src%sector_nm = ''
    !
    ! A bit of other stuff
    !
    nullify(wb_dust_src%fZDisp)
    nullify(wb_dust_src%levFractDispVert)
    !
    ! Main source parameters - enough to identify it in the global information list
    !
    wb_dust_src%src_nbr = iSrcNbr
    wb_dust_src%id_nbr = iSrcIdNbr
    !
    ! Finally, mark the source as incomplete
    !
    wb_dust_src%defined = silja_false

  end subroutine reserve_wb_dust_source


  !*********************************************************************

  subroutine init_emission_wb_dust(srcWBDust, dispersionMarketPtr, start_time)
    !
    ! Initializes the wind-blown dust source term.
    ! The parameters are read from the given ini file and stored to the returned 
    ! source term. After the stuff has been read, the main FluxPerMode is generated.
    ! This configuration allows for several wind-blown dust sources
    !
    implicit none

    ! Imported parameters
    type(Tsilam_namelist), pointer :: nlSetup
    type(silam_wind_blown_dust_source), intent(inout) :: srcWBDust
    type(mini_market_of_stacks), pointer :: dispersionMarketPtr
    type(silja_time), intent(in) :: start_time

    ! Local variables
    integer :: iTmp, iSpecies, iFlds, iZ, iWind, iThresh
    type(silam_sp) :: sp
    real :: fTmp, fVolFlux, fMassMeanDiam
    type(silja_field_id) :: id
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
    call add_input_needs_wb_dust_src(srcWBDust, &
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
      iSpecies = fu_species_index_src(srcWBDust, q_disp_stat(iFlds))
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
                                 & srcWBDust%species(iSpecies), &
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
    call msg('Filling-in the wind-blown dust emission info from supplementary fields')
    call fill_minimarket_from_namelist(dispersionMarketPtr, &
                                     & srcWBDust%nlInputFiles, 'supplementary_file', & ! namelist and item
                                     & shop_list, start_time, &
                                     & static_climatology, &  ! target stack
                                     & create_field, &            ! error if a clash
                                     & wdr_missing, &
                                     & dispersion_gridPtr, &
                                     & 5, .true., & ! iAccuracy, ifAdjustGrid
                                     & ifOK)
!call msg('wbdust stop')
!stop
    ! Check that all the requested quantities are in the stack.
    ! Some of them might be initialized but if not, can still be added.
    !
    select case(srcWBDust%emisMethod)

      case(emis_sandblasting_v1_flag)
        !                      
        ! First of all, read the source mask
        !                                                                                                                                                                                                                                                                                                                   
        sp%sp => fu_work_string()
        if(error)return
        arPtr => fu_work_array()
        if(error)return

        sp%sp = fu_content(srcWBDust%nlInputFiles, 'source_area_mask')
        if(error .or. len_trim(sp%sp) < 1)then
          call set_error('Source area mask is absent','init_emission_wb_dust')
          return
        endif
        id = fu_set_field_id_simple(met_src_missing, fraction_of_erodible_soil_flag, &
                              & time_missing, level_missing)
        if(error)return
        call set_grid(id, dispersion_grid)
        if(error)return

        call get_input_field(fu_process_filepath(sp%sp(index(adjustl(sp%sp),' ')+1:), &
             & superdir=srcWBDust%src_data_dir), &  ! file name      
             & fu_input_file_format(sp%sp), &          ! file format          
             & id, &                 ! The id to search               
             & arPtr, &              ! data array                                                         
             & dispersion_gridPtr, & ! storage grid                                          
             & .false., &            ! ifAcceptSameMonth                                      
             & nearestPoint, &       ! out of grid interpolation                                                       
             & id, &                 ! idOut                                                   
             & 5, &                  ! iAccuracy                                                                                                     
             & wdr_missing,&
             & 0.0)                  ! Fillvalue                                          
        !                       & iBinary, &       ! for NetCDF                                                                                                                                                                                                                                                                     
        if(error)return

        if(defined(id)) then
          allocate(srcWBDust%source_mask, stat=iTmp)
          if(fu_fails(iTmp == 0,'Failed wind-blown dust source mask allocation','init_emission_wb_dust'))return
          call set_field(id, arPtr, srcWBDust%source_mask, .true.)
        else
          call set_error('Failed to get the source mask','init_emission_wb_dust')
        endif
        if(error)return

        id = fu_set_field_id_simple(met_src_missing,&
                                  & alluvial_sedim_index_flag, &
                                  & time_missing, &        ! valid time
                                  & level_missing)
        fieldPtr => fu_get_field_from_mm_general(dispersionMarketPtr, id, .false.)
        if(.not. associated(fieldPtr))then
          call set_error('Alluvial sediment map is missing but needed','init_emission_wb_dust')
          return
        endif

      case(simple_dust_flag)
        id = fu_set_field_id_simple(met_src_missing,&
             & dust_emis_0_flag, &
             & time_missing, &        ! valid time                                                                                                                                                                                                                                                                        
             & level_missing)
        fieldPtr => fu_get_field_from_mm_general(dispersionMarketPtr, id, .false.)
        if(.not. associated(fieldPtr))then
          call set_error('Dust emission map is missing but needed','init_emission_wb_dust')
          return
        endif

      case default  ! Do nothing. Actually, most of the fields are taken from physiography
    end select

    !
    ! Allocate proper size of other arrays
    !
    allocate(srcWBDust%fluxPerModeNbr(srcWBDust%nSpecies), &
           & srcWBDust%fluxPerModeVol(srcWBDust%nSpecies), stat=iTmp)
    if(iTmp /= 0)then
      call msg('Allocation status:',iTmp)
      call set_error('Failed to allocate dust precomputed flux array', &
                   & 'init_emission_wb_dust')
      return
    endif

    !
    ! Now let's precompute the array of fluxes per size mode
    ! For each mode we have to take care of all temperatures and salinities
    !
    do iSpecies = 1, srcWBDust%nSpecies
      call wind_blown_dust_flux4mode(srcWBDust%species(iSpecies)%mode, &
                                   & srcWBDust%iSpectrumType, &
                                   & srcWBDust%fluxPerModeNbr(iSpecies), &
                                   & srcWBDust%fluxPerModeVol(iSpecies), &
                                   & fMassMeanDiam)
      !
      ! Store the results for this size mode
      !
      call set_massmean_d(srcWBDust%species(iSpecies)%mode, fMassMeanDiam, .false.)
      if(error)return

    end do    ! species
    !
    !
    ! Now, create the lookup table for the gust impact
    !
    ! Make the scanning axes and allocate space
    !
    srcWBDust%wgl%nZ = 50
    srcWBDust%wgl%nWinds = 60
    srcWBDust%wgl%nWindPwr = 50
    allocate(srcWBDust%wgl%fWind(srcWBDust%wgl%nWinds), srcWBDust%wgl%fZr(srcWBDust%wgl%nZ), &
           & srcWBDust%wgl%fWindPwr(srcWBDust%wgl%nWindPwr), stat=iWind)
    if(fu_fails(iWind==0,'Failed gust lolokup axes allocation','init_emission_wb_dust'))return
    do iWind = 1, srcWBDust%wgl%nWinds              ! wind speed 
      if(iWind == 1)then
        srcWBDust%wgl%fWind(iWind) = 0.001   !m/s
      else
        srcWBDust%wgl%fWind(iWind) = real(iWind-1) / 10.0  ! m/s  uStar up to 6m/s...
      endif
    end do
    do iZ = 1, srcWBDust%wgl%nZ                    ! roghness length
      srcWBDust%wgl%fZr(iZ) = 1e-5 * 1.3**iZ    ! z0, m - up to 5m
    end do
    do iThresh = 1, srcWBDust%wgl%nWindPwr      ! wind speed threshold, just borrow windpower 
      if(iThresh == 1)then
        srcWBDust%wgl%fWindPwr(iThresh) = 0.001   !m/s
      else
        srcWBDust%wgl%fWindPwr(iThresh) = real(iThresh-1) / 10.0  ! m/s
      endif
    end do

    allocate(srcWBDust%wgl%pGustVal(srcWBDust%wgl%nWinds, &
                                   & srcWBDust%wgl%nZ, srcWBDust%wgl%nWindPwr), stat=iZ)
    if(fu_fails(iZ==0,'Failed reallocation of windgust lookup table','init_emission_wb_dust'))return
    !
    ! Calculate the lookup table
    !
    do iZ = 1, srcWBDust%wgl%nZ
      do iWind = 1, srcWBDust%wgl%nWinds
        do iThresh = 1, srcWBDust%wgl%nWindPwr
!          srcWBDust%wgl%pGustVal(iWind,iZ,iThresh) = fu_gust_4_dust(srcWBDust%wgl%fWind(iWind), &
!                                                                  & srcWBDust%wgl%fZr(iZ), &
!                                                                  & srcWBDust%wgl%fWindPwr(iThresh))
          srcWBDust%wgl%pGustVal(iWind,iZ,iThresh) = &
                               & fu_gust_4_dust_satur(srcWBDust%wgl%fWind(iWind), &
                                                                  & srcWBDust%wgl%fZr(iZ), &
                                                                  & srcWBDust%wgl%fWindPwr(iThresh))
          if(error)return
        end do  ! threshold
      end do  ! wind-dim
    end do  ! z-dim

!do iZ = 1, wgl%nZ
!write(run_log_funit,fmt='(1000(F5.2,1x))')(wgl%pGustVal(iWind,iZ,wgl%nWindPwr),iWind=1,wgl%nWinds)
!end do
  
    srcWBDust%defined = silja_true

    call free_work_array(sp%sp)
    call free_work_array(stack_quantities)
    call free_work_array(q_disp_dyn)
    call free_work_array(q_disp_stat)

    call report(srcWBDust)

  CONTAINS
  
    !==============================================================================

    real function fu_gust_4_dust(fUstar, z0, fUstarThresh)
      !
      ! Just integrates the normal distribution with the given parameters
      ! gust_integral = S[(u*-uThr)*(u*+uThr)**2 * p(du)] d du,
      ! where u*=u*_input + du
      ! Here du is the deviation of wind from mean value (gust)
      ! u*input is mean u* itself (fUstar)
      ! p is Gaussian function
      ! integral S is taken from minus to plus infinity.
      ! See notebook 11, p.38
      !
      implicit none

      ! Imported parameters
      real, intent(in) :: fUstar, z0, fUstarThresh

      ! Local variables
      integer :: iTmp
      real :: fInt, s, s_2, du, fTmp, fUstarActual
      real, parameter :: f1_sqrt_2Pi = 0.3989422804
      !
      ! Relative standard deviation: sigma/u = 2.185 * k / ln(Z/z0)
      ! k is von Karman = 0.4
      !
      fInt = 0.0
      fTmp = 0.0
      s = 2.185 * fUstar * karmann_c / log(10./z0)   ! sigma_u_star
      s_2 = 0.5 / (s * s)

      do iTmp = -1000, 1000
        du = s * real(iTmp) * 0.01          ! plus-minus 10 sigmas. Should be enough...
        fUstarActual = fUstar + du
        if(fUstarActual > fUstarThresh) fInt = fInt + (fUstarActual-fUstarThresh) * &
                                                    & (fUstarActual + fUstarThresh)**2 * &
                                                    & exp(-du*du*s_2) * 0.01 * s
        fTmp = fTmp + exp(-du*du*s_2) * 0.01 * s
      end do  ! interation cycle
 
      if(abs(fTmp*f1_sqrt_2Pi/s - 1.0) > 1e-5)then
        call msg('PDF integral not unity:',fTmp*f1_sqrt_2Pi/s)
        call set_error('PDF integral not unity','fu_gust_4_dust')
      endif
      fu_gust_4_dust = fInt * f1_sqrt_2Pi / s

! call msg('wind, z/z0, gust forcing factor:' + fu_str(fWind), fZ_z0, fu_gust_integral/fWind**fPwr)

    end function fu_gust_4_dust

    !=======================================================================

    real function fu_gust_4_dust_satur(fUstar, z0, fUstarThresh)
      !
      ! Just integrates the normal distribution with the given parameters
      ! gust_integral = S[(u*-uThr)*(u*+uThr)**2 * p(du)] d du,
      ! where u* = u*_input + du
      ! Here du is the deviation of wind from mean value (gust)
      ! u*input is mean u* itself (fUstar)
      ! Saturation is taken for the excess over the threshold
      ! p is Gaussian function
      ! integral S is taken from minus to plus infinity.
      ! See notebook 11, p.38
      !
      implicit none

      ! Imported parameters
      real, intent(in) :: fUstar, z0, fUstarThresh

      ! Local variables
      integer :: iTmp
      real :: fInt, s, s_2, du, fTmp, fUstarActual
      real, parameter :: f1_sqrt_2Pi = 0.3989422804
      !
      ! Relative standard deviation: sigma/u = 2.185 * k / ln(Z/z0)
      ! k is von Karman = 0.4
      !
      fInt = 0.0
      fTmp = 0.0
      s = 2.185 * fUstar * karmann_c / log(10./z0)   ! sigma_u_star
      s_2 = 0.5 / (s * s)

      do iTmp = -1000, 1000
        du = s * real(iTmp) * 0.01          ! plus-minus 10 sigmas. Should be enough...
        fUstarActual = fUstar + du

        if(fUstarActual > fUstarThresh) then
          !
          ! Saturate the excess
          !
          fUstarActual = fUstarThresh + (fUstarActual - fUstarThresh) * u_star_excess_saturation / &
                                      & (fUstarActual - fUstarThresh + u_star_excess_saturation)
!call msg('Saturation, before and after:',fTmp, uStarSaturated - u_star_min)
          !
          ! ...and add to the integral
          !
          fInt = fInt + (fUstarActual-fUstarThresh) * (fUstarActual + fUstarThresh)**2 * &
                                                    & exp(-du*du*s_2) * 0.01 * s
        endif
        fTmp = fTmp + exp(-du*du*s_2) * 0.01 * s
      end do  ! interation cycle
 
      if(abs(fTmp*f1_sqrt_2Pi/s - 1.0) > 1e-5)then
        call msg('PDF integral not unity:',fTmp*f1_sqrt_2Pi/s)
        call set_error('PDF integral not unity','fu_gust_4_dust')
      endif
      fu_gust_4_dust_satur = fInt * f1_sqrt_2Pi / s

! call msg('wind, z/z0, gust forcing factor:' + fu_str(fWind), fZ_z0, fu_gust_integral/fWind**fPwr)

    end function fu_gust_4_dust_satur

  end subroutine init_emission_wb_dust


  !****************************************************************************

  subroutine add_inventory_wb_dust_src(wb_dust_src, species_list, nSpecies)
    !
    ! Get the species emitted by this source and make a list of
    ! them. The source must be initialized, of course.
    !
    implicit none
    type(silam_wind_blown_dust_source), intent(in) :: wb_dust_src
    type(silam_species), dimension(:), pointer :: species_list
    integer, intent(inout) :: nSpecies

    call addSpecies(species_list, nSpecies, wb_dust_src%species, wb_dust_src%nSpecies)

  end subroutine add_inventory_wb_dust_src


  !**************************************************************************

  subroutine link_wb_dust_src_to_species(species_list, wb_dust_src)
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
    type(silam_wind_blown_dust_source), intent(inout) :: wb_dust_src
    type(silam_species), dimension(:), pointer :: species_list

    !
    ! Linkage is actually just creation of the chemical adaptor
    !
    call create_adaptor(wb_dust_src%species, species_list, wb_dust_src%adaptor)
    
  end subroutine link_wb_dust_src_to_species


  !************************************************************************

  integer function fu_species_index_src(wb_dust_src, quantity)
    !
    ! Selects the proper name of the substance for the given quantity - just to be able to 
    ! choose the right input and dispersion-stack fields for each source.
    ! Indices of main pollen species, pollen and free allergen are known by the source
    !
    implicit none

    ! Imported parameters
    type(silam_wind_blown_dust_source), intent(in) :: wb_dust_src
    integer, intent(in) :: quantity

    fu_species_index_src = int_missing  ! so far, nothing fancy
        
  end function fu_species_index_src


  !*******************************************************************************
  
  subroutine create_src_cont_grd_wbdust_src(wb_dust_src, grid_template, ifVerbose, ifMinimal, ifExtended)
    !
    ! Creates the grid that covers the area with active BVOC emission
    !
    implicit none

    ! Imported parameters
    type(silam_wind_blown_dust_source), intent(in) :: wb_dust_src
    type(silja_grid), intent(inout) :: grid_template
    logical, intent(in) :: ifVerbose, ifMinimal
    logical, intent(out) :: ifExtended

    ! So far nothing to do: this source just covers the dispersion grid
    !
    ifExtended = .false.
    return
    
  end subroutine create_src_cont_grd_wbdust_src


  !*****************************************************************

  subroutine project_wbdust_src_second_grd(wbdust_src, grid, vert_disp, vert_proj, iAccuracy)
    !
    ! There is nothing to project in terms of horizontal grid but the vertical has to be made.
    !
    implicit none

    ! Imported parameters
    type(silam_wind_blown_dust_source), intent(inout) :: wbdust_src
    type(silja_grid), intent(in) :: grid
    type(silam_vertical), intent(in) :: vert_disp, vert_proj
    integer, intent(in) :: iAccuracy

    ! Local variables
    integer :: i, ilev
    type(silam_vertical) :: vertTmp
    real, dimension(:), pointer :: arTmp, dz_disp

    if(iAccuracy < 1 .or. iAccuracy > 10)then
      call msg('Accuracy switch must be from 1 to 10, not:',iAccuracy)
      call msg_warning('Accuracy switch must be from 1 to 10','project_wbdust_src_second_grd')
!      return
    endif

    if(len_trim(wbdust_src%sector_nm) > 0)then
      call msg('Re-projecting wind-blown dust source:' + wbdust_src%src_nm +'_' + wbdust_src%sector_nm)
    else
      call msg('Re-projecting wind-blown dust source:' + wbdust_src%src_nm)
    endif

    !
    ! Now re-project the vertical grid of the source to the given vertical
    ! Since vertical distributions are very poorly known and crude, we do not need
    ! any precise meteo-dependent projection, crude will do the job
    !
    wbdust_src%vertLevsDispVert = vert_disp
    allocate(wbdust_src%levFractDispVert(fu_NbrOfLevels(vert_disp)), &
           & wbdust_src%fzDisp(fu_NbrOfLevels(vert_disp)), stat=i)
    if(fu_fails(i==0, 'Failed to allocate dispersion-vertical level fractions','project_wbdust_src_second_grd'))return
    wbdust_src%levFractDispVert(:) = 0.0
    wbdust_src%fzDisp(:) = 0.0

    arTmp => fu_work_array()
    if(error)return
    !
    ! Create the wind-blown dust vertical, which can depend on the emission method
    !
    select case(wbdust_src%emisMethod)
      
      case(emis_Gillette_DMAT, emis_sandblasting_v1_flag)
        !
        ! So far all the methods play the same profile
        !
        call set_vertical(fu_set_layer_between_two(layer_btw_2_height, 1.0, 1000.0), vertTmp)
        if(error)return
        call add_level(vertTmp, fu_set_layer_between_two(layer_btw_2_height, 1000.0, 2000.0))
        if(error)return
        call add_level(vertTmp, fu_set_layer_between_two(layer_btw_2_height, 2000.0, 3000.0))
        if(error)return
        call add_level(vertTmp, fu_set_layer_between_two(layer_btw_2_height, 3000.0, 5000.0))
        if(error)return
        arTmp(1) = 0.5
        arTmp(2) = 0.25
        arTmp(3) = 0.125
        arTmp(4) = 0.125

      case(simple_dust_flag)
        call set_vertical(fu_set_layer_between_two(layer_btw_2_height, 1.0, 100.0), vertTmp)
        if(error)return
        arTmp(1) = 1.0

      case default
        call set_error('Unknown emission method:' + fu_str(wbdust_src%emisMethod), &
                     & 'project_wbdust_src_second_grd')
        return
    end select

    call reproject_verticals(vertTmp, arTmp, &                   ! vertical from, fractions from
                           & vert_proj, wbdust_src%levFractDispVert, &   ! vertical to, fractions to
                           ! mass centres, number of non-zero levels:
                           & wbdust_src%fzDisp, wbdust_src%nLevsDispVert, &
                           & ifMassCentreInRelUnit=.true.)

    call set_missing(vertTmp, .false.)
    call free_work_array(arTmp)

  end subroutine project_wbdust_src_second_grd


  !**************************************************************************

  subroutine compute_emission_for_wb_dust(wb_dust_src, &
                                        & met_buf, disp_buf, & 
                                        & now, &      ! current time
                                        & timestep, & ! model time step
                                        & pHorizInterpMet2DispStruct, ifHorizInterp, &
                                        & ifSpeciesMoment, &
                                        & emisMap, mapCoordX, mapCoordY, mapCoordZ, & ! Output
                                        & fMassInjected)                              ! Output
    !
    ! Computes the emission fields for wind_blown_dust. Since the dust is a wide-specturm
    ! aerosol, we have to compute the emission flux for each size class separately.
    ! Hence, the flux is a function of the wind speed and type of the surface. 
    ! The latter decides on the spectrum and on the wind speed dependence.
    !
    ! This routine is to be called at each model time step but not inside the main cycle,
    ! therefore its efficiency is of moderate importance.
    !
    implicit none

    ! Imported parameters
    type(silam_wind_blown_dust_source), target, intent(in) :: wb_dust_src
    type(Tfield_buffer), pointer ::  met_buf, disp_buf  ! meteo and internal field buffers
    type(silja_time), intent(in) :: now           ! current time
    type(silja_interval), intent(in) :: timestep  ! model time step
    type(THorizInterpStruct), intent(in) :: pHorizInterpMet2DispStruct
    logical, intent(in) :: ifHorizInterp, ifSpeciesMoment
    type(Tmass_map), intent(inout) :: emisMap, mapCoordX, mapCoordY, mapCoordZ
    real(r8k), dimension(:), intent(inout) :: fMassInjected


    ! Local variables
    integer :: indLandFr, indW10m, indFricVel, iDust, iLev, iWind, iThresh, iZ, &
             & iMeteo, iDisp, ix, iy, ixMeteo, iyMeteo, iStat, iMode, iTmp, jTmp, iMeteoTmp
    real :: fTmp, timestep_sec, fCellTotal, ro_soil_dry, fSoilMassWater, u_star_min, &
          & uStarMerged, fSaltation, fFluxTotal, fRoughness, fSoilMassWaterThresh, fWindSpeed, &
          & zx, zy, zz, fSum, du_star_owen, uStarGustEff !, uStarSaturated
    real, dimension(:), pointer :: ptrXSzDisp, ptrYSzDisp, pSoilVolumeWater, &
                                 & fMinDArray, fMaxDArray, fPartDensity, pSandMassFract, &
                                 & pErodibleFraction, pU10m, pV10m, pU_star, pLAI, &
                                 & pClayMassFraction, pAlluvialSedimentIndex, pLandFr, &
                                 & pSnowWEQDepth, pLandRough, ptrLonMeteo, ptrLatMeteo, &
                                 & pDustEmis0
    real*4, dimension(worksize) :: arTmp
    type(silja_field), pointer :: fldMaskPtr
    integer, save :: iCount = 0, iCountTalking = 0, uDump, iRec
    logical, save :: ifFirst=.true., ifTalking=.true.

#ifdef DEBUG_DUST
    real*4, dimension(:,:,:), pointer, save :: fArDump
    real :: corner_lon_W, corner_lat_S, southpole_lon_E, southpole_lat_N, dx_deg, dy_deg
    logical :: corner_in_geo_lonlat
    integer :: number_of_points_x, number_of_points_y
#endif
    integer, dimension(10), save :: iArStatistics=0

    type(silam_sp) :: sp

    ! Local parameters 
    !
    ! Tuning.
    ! Long way towards tighter smooth roughness: dow to 1e-5 (note: 2.3e-6 turns everything off everywhere)
    ! Also, tightening the empirical saturation: from no limit to 0.1 and then 0.06 m/s of excess u*
    ! 17.4.2016. Much too limited emission area with too strong emission there: turn both parameters, carefully
    ! 19.4.2016 Further tightening the u* excess saturation but it comes to its end. A map of smooth-Z0 is needed
    !
    real, parameter :: fSandDensity = 2600, &
                     & fSandDensity_1_3 = fSandDensity ** 0.3333333, &
                     & fSnowWEQDepthThreshold = 0.001, &     ! 1mm of snow => no dust
                     & sqrt_density_air = sqrt(density_air_288K)
    logical, parameter :: ifWindGust = .true.

    !
    ! Meteo quantities
    !
!    if(fu_fails(fu_index(met_buf, soil_moisture_content_flag, pSoilVolumeWater) /= int_missing, &  ! soil water content, m3/m3
    if(fu_fails(fu_index(met_buf, soil_moisture_vol_frac_nwp_flag, pSoilVolumeWater) /= int_missing, &  ! soil water content, m3/m3
                                & 'Failed soil moisture input','compute_emission_for_wb_dust'))return
    if(fu_fails(fu_index(met_buf, u_10m_flag, pU10m) /= int_missing, &                          ! u10m, m/s
                                & 'Failed 10m u','compute_emission_for_wb_dust'))return
    if(fu_fails(fu_index(met_buf, v_10m_flag, pV10m) /= int_missing, &                          ! v10m, m/s
                                & 'Failed 10m v','compute_emission_for_wb_dust'))return
    if(fu_fails(fu_index(met_buf, water_eq_snow_depth_flag, pSnowWEQDepth) /= int_missing, &    ! snow water equiv. depth
                                & 'Failed water_eq_snow_depth_flag','compute_emission_for_wb_dust'))return
    if(fu_fails(fu_index(met_buf, leaf_area_index_flag, pLAI) /= int_missing, &                 ! leaf area index
                                & 'Failed leaf area index','compute_emission_for_wb_dust'))return
    if(fu_fails(fu_index(met_buf, fraction_of_land_flag, pLandFr) /= int_missing, &             ! fraction of land                
                                & 'Failed fraction_of_land_flag','compute_emission_for_wb_dust'))return


    select case(wb_dust_src%emisMethod)
      case(emis_sandblasting_v1_flag)
        pErodibleFraction => fu_grid_data(wb_dust_src%source_mask)

#ifdef DEBUG_DUST
    call lonlat_grid_parameters(dispersion_grid,&
                       & corner_lon_W, corner_lat_S,corner_in_geo_lonlat, &
                       & number_of_points_x, number_of_points_y, &
                       & southpole_lon_E, southpole_lat_N, &
                       & dx_deg, dy_deg)
    open (50, file='wb_dust_erodible_fr.grads',recl=nx_dispersion*ny_dispersion,form='unformatted',access='direct')
    arTmp(1:nx_dispersion*ny_dispersion) = pErodibleFraction(1:nx_dispersion*ny_dispersion)
    write(50,rec=1)arTmp(1:nx_dispersion*ny_dispersion)
!    write(50,rec=1)(pErodibleFraction(ix),ix=1,nx_dispersion*ny_dispersion)                                                                                                  
    close (50)
    open(50,file='wb_dust_erodible_fr.ctl')
    write(50,*) 'DSET ^wb_dust_erodible_fr.grads'
    write(50,*) 'TITLE SILAM GrADS output'
    write(50,*) 'UNDEF    -0.999999E+15'
    write(50,*) 'XDEF   ',nx_dispersion,' LINEAR   ', corner_lon_W, dx_deg
    write(50,*) 'YDEF  ', ny_dispersion,' LINEAR   ', corner_lat_S, dy_deg
    write(50,*) 'ZDEF 1 LEVELS 24'
    write(50,*) 'TDEF  1  LINEAR  ', fu_time_to_grads_string(now), ' ', fu_interval_to_grads_string(timestep)
    write(50,*) 'VARS   1'
    write(50,*) 'erodible_fr 0 99 99 0 mask'
    write(50,*) 'ENDVARS'
    close(50)
#endif

        if(fu_fails(fu_index(met_buf, friction_velocity_flag, pU_star) /= int_missing, &            ! u*, m/s
             & 'Failed u_star','compute_emission_for_wb_dust'))return
        if(fu_fails(fu_index(met_buf, soil_sand_mass_fraction_flag, pSandMassFract) /= int_missing, &  ! sand mass fraction 
             & 'Failed soil_sand_mass_fraction_flag','compute_emission_for_wb_dust'))return
        if(fu_fails(fu_index(met_buf, soil_clay_mass_fraction_flag, pClayMassFraction) /= int_missing, &  ! clay mass fraction 
             & 'Failed soil_clay_mass_fraction_flag','compute_emission_for_wb_dust'))return
        if(fu_fails(fu_index(met_buf, land_roughness_disp_flag, pLandRough) /= int_missing, &       ! dispersion roughness  
             & 'Failed land_roughness_disp_flag','compute_emission_for_wb_dust'))return
        !                                                      
        ! Alluvial sediments are in the dispersion buffer. Reason: they have to be averaged   
        ! when doing the grid reprojection, not interpolated.                                           
        !                                                                            
        if(fu_fails(fu_index(disp_buf, alluvial_sedim_index_flag, pAlluvialSedimentIndex) /= int_missing, &    ! alluvial sediment
             & 'Failed alluvial_sedim_index_flag','compute_emission_for_wb_dust'))return
      case(simple_dust_flag)
        if(fu_fails(fu_index(disp_buf, dust_emis_0_flag, pDustEmis0) /= int_missing, &    ! basic map for dust emission        
             & 'Failed dust_emis_0_flag','compute_emission_for_wb_dust'))return
      case default
    end select
    ptrXSzDisp => fu_grid_data(dispersion_cell_x_size_fld)  ! get grid cell size in metres
    ptrYSzDisp => fu_grid_data(dispersion_cell_y_size_fld)

    ptrLonMeteo => fu_grid_data(meteo_longitude_fld)
    ptrLatMeteo => fu_grid_data(meteo_latitude_fld)

#ifdef DEBUG_DUST
    if(iCount == 0)then
      call lonlat_grid_parameters(meteo_grid,&
                         & corner_lon_W, corner_lat_S,corner_in_geo_lonlat, &
                         & number_of_points_x, number_of_points_y, &
                         & southpole_lon_E, southpole_lat_N, &
                         & dx_deg, dy_deg)    
      open(50, file='roughness.grads', form='unformatted', access='direct', recl=fs_meteo)
      arTmp(1:fs_meteo) = pLandRough(:)
      write(50,rec=1)arTmp(1:fs_meteo)
      close(50)
!    call msg('Dispersion grid parameters: ', (/corner_lon_E, corner_lat_N, &
!                                             & southpole_lon_E, southpole_lat_N, &
!                                             & dx_deg, dy_deg/))
!    call report(dispersion_grid)
    open(50,file='roughness.ctl')
    write(50,*) 'DSET ^roughness.grads'
    write(50,*) 'TITLE SILAM GrADS output'
    write(50,*) 'UNDEF    -0.999999E+15'
    write(50,*) 'XDEF   ',nx_meteo,' LINEAR   ', corner_lon_W, dx_deg
    write(50,*) 'YDEF  ', ny_meteo,' LINEAR   ', corner_lat_S, dy_deg
    write(50,*) 'ZDEF 1 LEVELS 24'
    write(50,*) 'TDEF  1  LINEAR  ', fu_time_to_grads_string(now), ' ', fu_interval_to_grads_string(timestep)
    write(50,*) 'VARS   1'
    write(50,*) 'roughness 0 99 99 0'
    write(50,*) 'ENDVARS'
    close(50)
    endif
#endif

    timestep_sec = abs(fu_sec(timestep))
    if(error)return

    !-------------------------------------------------------------------
    !
    ! The total number flux released per size class is in the pre-computed array,
    ! from which the values should just be picked.
    !
    ! Particle masses for each size class: !!!! dry particle !!!!
    !
    fPartDensity => fu_work_array()

    do iDust = 1, wb_dust_src%nSpecies
      fPartDensity(iDust) = fu_dry_part_density(wb_dust_src%species(iDust)%material)
    end do

#ifdef DEBUG_DUST
if(ifFirst)then
  allocate(fArDump(nx_dispersion, ny_dispersion, 8))

  call msg('')
  call msg('Reporting the dust source settings:')
  call msg('fEntrainmentScaling',fEntrainmentScaling)
  call msg('fSnowWEQDepthThreshold',fSnowWEQDepthThreshold)
  call msg('fSrfRoughSandOnly',fSrfRoughSandOnly)
  call msg('fDragPartitioningSandOnly = 1.0/log(0.35*(0.1/fSrfRoughSandOnly)**0.8)')
  call msg('fLAI_critical',fLAI_critical)

    call lonlat_grid_parameters(dispersion_grid,&
                       & corner_lon_W, corner_lat_S,corner_in_geo_lonlat, &
                       & number_of_points_x, number_of_points_y, &
                       & southpole_lon_E, southpole_lat_N, &
                       & dx_deg, dy_deg)    

   uDump = fu_next_free_unit()
   open(uDump,file='wind_blown_dump.ctl')
   write(uDump,*) 'DSET ^wind_blown_dump_%m2%d2.grads'
   write(uDump,*) 'TITLE SILAM GrADS output'
   write(uDump,*) 'OPTIONS TEMPLATE'
   write(uDump,*) 'UNDEF    -0.999999E+15'
   write(uDump,*) 'XDEF   ',nx_dispersion,' LINEAR   ', corner_lon_W, dx_deg
   write(uDump,*) 'YDEF  ', ny_dispersion,' LINEAR   ', corner_lat_S, dy_deg
   write(uDump,*) 'ZDEF 1 LEVELS 24'
   write(uDump,*) 'TDEF  1000  LINEAR  ', fu_time_to_grads_string(now), ' ', fu_interval_to_grads_string(timestep)
   write(uDump,*)'VARS 8'
   write(uDump,*)'DragSplit 0 99 99 0 drug split'
   write(uDump,*)'SoilWaterThresh 0 99 99 0 u*min'
   write(uDump,*)'soilwater 0 99 99 0 u*min'
   write(uDump,*)'uStarMin 0 99 99 0 u*min'
   write(uDump,*)'uStarGustEff 0 99 99 0 effective uStar'
   write(uDump,*)'uStarMet 0 99 99 0 u*Owen'
   write(uDump,*)'Saltation 0 99 99 0 saltation'
   write(uDump,*)'fFluxTotal 0 99 99 0 flux total'
   write(uDump,*) 'ENDVARS'
   close(uDump)
   open(uDump,file='wind_blown_dump_' + fu_str(fu_mon(now),2) + fu_str(fu_day(now),2) + '.grads',form='unformatted', &
                       & access='direct', recl=nx_dispersion*ny_dispersion)
   iRec = 1
   ifFirst = .false.
 call msg('fDragPartitioningSandOnly:',fDragPartitioningSandOnly)
 
endif
 
if(fu_hour(now) == 0 .and. fu_min(now) == 0)then
  close(uDump)
  open(uDump,file='wind_blown_dump_' + fu_str(fu_mon(now),2) + fu_str(fu_day(now),2) + '.grads',form='unformatted', &
            & access='direct', recl=nx_dispersion*ny_dispersion)
  iRec = 1
endif
fArDump = 0.0
#endif
 

    !
    ! Now scan the domain, for each grid cell find the dust emission parameters
    ! and fill-in the emission map
    !
    do iy = 1, ny_dispersion
      do ix = 1, nx_dispersion
        iDisp = ix + (iy-1) * nx_dispersion

iArStatistics(1) = iArStatistics(1) + 1
        !
        ! No emission if:
        !   ... no source area fraction in this cell
        !
        select case(wb_dust_src%emisMethod)
          case(emis_sandblasting_v1_flag)
            if(pErodibleFraction(iDisp) < 1.0e-5)cycle
          case(simple_dust_flag)
            if(pDustEmis0(iDisp) < 1.0e-18)cycle
          case default
        end select

        iMeteo =  fu_grid_index(nx_meteo, ix, iy, pHorizInterpMet2DispStruct)
        iyMeteo = int(iMeteo / nx_meteo)
        ixMeteo = iMeteo - (iyMeteo-1)*nx_meteo
        !
        !   ... and separately if meteo-land mask is almost zero
        !
        if(pLandFr(iMeteo) < 1.0e-5)cycle

iArStatistics(2) = iArStatistics(2) + 1

        !
        !   ...and if soil covered with snow
        !
        if(pSnowWEQDepth(iMeteo) > fSnowWEQDepthThreshold)cycle

iArStatistics(3) = iArStatistics(3) + 1

        !
        !   ...or vegetation
        !
        if(pLAI(iMeteo) > fLAI_critical)cycle

iArStatistics(4) = iArStatistics(4) + 1

        select case(wb_dust_src%emisMethod)
          case(emis_sandblasting_v1_flag)
            !   ...or very rough soil or problematic roughness
            !
            fRoughness = max(pLandRough(iMeteo), fSrfRoughSandOnly*2.1)
            if (1. < fDragPartitioningSandOnly * (log(fRoughness/fSrfRoughSandOnly)) )cycle
            iArStatistics(5) = iArStatistics(5) + 1
          case default
            iArStatistics(5) = iArStatistics(5) + 1
        end select
            
        !
        ! Soil moisture can be tricky: for small islands and sea edges it can be set to zero. 
        ! But it can also be zero for deserts. So, for partly sea-cells we will check the surrounding 
        ! mainly land cells if they give any clue on actual moisture in the area
        !
        if(pLandFr(iMeteo) < 0.75)then
          if(abs(pSoilVolumeWater(iMeteo)) < 1.0e-5)then
if(ifTalking)call msg('Checking suspicious soil moisture:',ix,iy)
            iCount = 0
            fTmp = 0.0
            do jTmp = max(iyMeteo-1,1), min(iyMeteo, ny_meteo)
              do iTmp =  max(ixMeteo-1,1), min(ixMeteo, nx_meteo)
                iMeteoTmp = iTmp + nx_meteo * (jTmp -1)
                if(iMeteo == iMeteoTmp)cycle
                if(pLandFr(iMeteoTmp) > 0.75)then
                  fTmp = fTmp + pSoilVolumeWater(iMeteo)
                  iCount = iCount + 1
                endif
              end do
            end do
            if(iCount > 0)then
              fTmp = fTmp / real(iCount)
              if(fTmp > 0.01)then
                 pSoilVolumeWater(iMeteo) = fTmp  ! surrounding is wet: do not believe zero!
if(ifTalking)call msg('Replace zero moisture at (' + fu_str(ix) + ',' + fu_str(iy) + ') with:', fTmp)
              else
if(ifTalking)call msg('Keep zero moisture')
              endif
            else
if(ifTalking)call msg('Replace zero moisture at (' + fu_str(ix) + ',' + fu_str(iy) + ') with default 0.3')
              pSoilVolumeWater(iMeteo) = 0.3  ! virtually anything: this is a kind-of isolated island
            endif
          endif ! soil moisture ~0
        endif  ! significant part of water in cell 
        
        select case(wb_dust_src%emisMethod)
          case(simple_dust_flag)
            
            fWindSpeed = sqrt(pU10m(iMeteo)*pU10m(iMeteo)+pV10m(iMeteo)*pV10m(iMeteo))
            fSaltation = max(1.4, fWindSpeed-4.1)**3 
            fSaltation = fSaltation * max(0.0, ((1 - 5.*pSoilVolumeWater(iMeteo))))**2

            if(fSaltation <= 0) cycle
            
            ! pDustEmis0 refers to a precomputed scaling map that depends on the bare land fraction
            ! (i.e. source area mask), the surface rughness map and the sand fraction map.
            ! Currently it equals scaling_factor * bare_land/roughness**0.5 * factor_based_on_land_features
            ! 1/roughness**0.5 is numerically quite close to 1-log(roughness/some_roughness_scale),
            ! but does not completely cut out the emission in rough areas (e.g. dust storms do occur
            ! in the midwest USA). pDustEmis needed to be preprocessed to cut out mismatching map
            ! pixels mainly along the coasts.

            fFluxTotal = 0.135 * pDustEmis0(iDisp) * fSaltation * timestep_sec * max(0.0, fLAI_critical-pLAI(iMeteo))

            if(fSaltation > 0.0) iArStatistics(6) = iArStatistics(6) + 1
            if(fFluxTotal > 0.0) iArStatistics(7) = iArStatistics(7) + 1
            
          case(emis_sandblasting_v1_flag)
            ro_soil_dry = 2600. * (0.51 + 0.126 * pSandMassFract(iMeteo))
            fSoilMassWater = pSoilVolumeWater(iMeteo) * 1000. / ro_soil_dry  ! gravimentric soil moisture
            !
            ! Saltation and then sandblasting. Modified Malticorena & Bergametti 1985, Zender 2003, etc
            ! Humidity is from Fecan et al, 1999
            !
            u_star_min = 0.03 * &                  ! default fit is 0.0466
                       & fSandDensity_1_3 / &       ! ro_sand ** 1/3
                       & sqrt_density_air / &       ! ro_air
                       & (1. - fDragPartitioningSandOnly * (log(fRoughness/fSrfRoughSandOnly))) 

#ifdef DEBUG_DUST
! Drag partition
fArDump(ix,iy,1) = (1. - fDragPartitioningSandOnly * (log(max(fRoughness,fSrfRoughSandOnly)/fSrfRoughSandOnly)))
#endif
IF(pU_star(iMeteo) > u_star_min)iArStatistics(6) = iArStatistics(6) + 1  ! strong enough wind?

            fSoilMassWaterThresh = 0.14 * pClayMassFraction(iMeteo)**2 + &
                                   & 0.17 * pClayMassFraction(iMeteo)           ! Fecan ea, 1999, as fraction

            if(fSoilMassWater > fSoilMassWaterThresh)then  
              !
              ! Correction only starting from some moisture. Fecan's formula with %->fraction recalc.
              !
              u_star_min = u_star_min * sqrt(1. + 27.7* (fSoilMassWater-fSoilMassWaterThresh)**0.68)
            endif

#ifdef DEBUG_DUST
fArDump(ix,iy,2) = fSoilMassWaterThresh
fArDump(ix,iy,3) = fSoilMassWater
fArDump(ix,iy,4) = u_star_min
#endif

            !
            ! Take Owen's effect:
            ! u* = u*_non_salt + 0.3(U10-U10_thr)**2,
            ! where U10_thr is some velocity. 
            ! So far, I have not got Owen's original work of 1964. From charts in Gillette et al, 1997
            ! JGR, D 22, 2S,989-25,998, it is around 8m/s for U1m, i.e. roughly 20m/s for U10m.
            ! The quadratic dependence and the coefs are taken from Zender, who refers to Gillette, 1998
            ! where, however, nothing like this relation has been obtained. Marticorena, 1995, is another
            ! source but it is French-language thesis, not available from internet. So, until Owen work is
            ! in hands, this relation will stay.
            !
            fWindSpeed = sqrt(pU10m(iMeteo)*pU10m(iMeteo)+pV10m(iMeteo)*pV10m(iMeteo))
!
! Owen effect is commented out following the Kok et al analysis where they claim this is
! wrong idea conflicting with the sandblasting concept.
!
!            if(fWindSpeed > 20.0)then
!              du_star_owen = 0.003 * (fWindSpeed - 20.0)**2
!              if(du_star_owen > 1.5)then
!                call msg('very high Owen for U10: ', du_star_owen, fWindSpeed)
!                du_star_owen = 1.5
!              endif
!            else
              du_star_owen = 0.0
!            endif

IF(pU_star(iMeteo) > u_star_min)iArStatistics(7) = iArStatistics(7) + 1  ! strong wind without gust?



!fTmp = pU_star(iMeteo) / karmann_c * log(10./fRoughness)
!if(fTmp > fWindSpeed * 5.)then
!  ifTalking = .true.
!  iCountTalking = 0
!  call msg('Low wind speed with high uStar:',fWindSpeed, pU_star(iMeteo))
!  call msg('ix, iy:',ix,iy)
!  fWindSpeed = fTmp
!endif

            !
            ! There can be strong excess of the u* over its minimal value, leading to unrealistic 
            ! saltation and flux.
            ! Saturated u* accounts for the fact that very strong wind lifts so much sand that the 
            ! stratification changes, then wind profile adjusts and becomes a limiting factor.
            ! No theory here, just fitting the saturation level. Kok et al showed some more options
            !
!            uStarSaturated = pU_star(iMeteo) + du_star_owen
!
! Saturation has been moved to the gust calculating formula
!
!fTmp = uStarSaturated - u_star_min
!            if(uStarSaturated > u_star_min) then
!              uStarSaturated = u_star_min +(uStarSaturated - u_star_min) * u_star_excess_saturation / &
!                                          & (uStarSaturated - u_star_min + u_star_excess_saturation)
!call msg('Saturation, before and after:',fTmp, uStarSaturated - u_star_min)
!            endif

            if(ifWindGust)then
              !
              ! Effective u* due to gustiness: use the lookup table for (u*-u*min)*(u*+u*min)**2,
              ! which takes into account gust distribution. Note that even for slow wind, there may be 
              ! strong gusts - no hard thresholds here.
              !
              iThresh =  max(1,int(u_star_min * 10.0))+1    ! approx from below
              if(iThresh > wb_dust_src%wgl%nWindPwr-1)cycle  ! above 5 m/s, forget...

              iWind = max(1,int((pU_star(iMeteo)) * 10.0))+1    ! approx from below
!              iWind = max(1,int(uStarSaturated * 10.0))+1    ! approx from below
              if(fu_fails(iWind <= wb_dust_src%wgl%nWinds-1,'Very strong wind', &
                                                     & 'compute_emission_for_wb_dust'))return
              iZ = max(1,int(log(fRoughness/1e-5) / log(1.3)))   ! approx from below
              if(fu_fails(iZ <= wb_dust_src%wgl%nZ-1,'Very high z0', &
                                                             & 'compute_emission_for_wb_dust'))then
                call unset_error('compute_emission_for_wb_dust')
                iZ = wb_dust_src%wgl%nZ-2
              endif
              !
              ! Here we take bilinear interpolation. Not the most-accurate admittedly but still...
              !
              zx = (pU_star(iMeteo) - wb_dust_src%wgl%fWind(iWind)) / &
!              zx = (uStarSaturated - wb_dust_src%wgl%fWind(iWind)) / &
                      & (wb_dust_src%wgl%fWind(iWind+1) - wb_dust_src%wgl%fWind(iWind))
              zy = (fRoughness - wb_dust_src%wgl%fZr(iZ)) / &
                      & (wb_dust_src%wgl%fZr(iZ+1) - wb_dust_src%wgl%fZr(iZ))
              zz = (u_star_min - wb_dust_src%wgl%fWindPwr(iThresh)) / &
                      & (wb_dust_src%wgl%fWindPwr(iThresh+1) - wb_dust_src%wgl%fWindPwr(iThresh))
              fSum = (1.-zx)*(1.-zy)*(1.-zz) + (1.-zx)*(1.-zy)*zz + zx*(1.-zy)*(1.-zz) + zx*(1.-zy)*zz + &
                   & (1.-zx)*zy*(1.-zz) + (1.-zx)*zy*zz + zx*zy*(1.-zz) + zx*zy*zz

              uStarGustEff = ((1.-zx)*(1.-zy)*(1.-zz)* wb_dust_src%wgl%pGustVal(iWind,iZ, iThresh) + &
                          & (1.-zx)*(1.-zy) *   zz * wb_dust_src%wgl%pGustVal(iWind,iZ, iThresh+1) + &
                          & zx *  (1.- zy) *(1.-zz)* wb_dust_src%wgl%pGustVal(iWind+1,iZ,iThresh) + &
                          & zx *  (1. - zy) *  zz  * wb_dust_src%wgl%pGustVal(iWind+1,iZ,iThresh+1) + &
                          & (1. - zx) * zy *(1.-zz)* wb_dust_src%wgl%pGustVal(iWind,iZ+1,iThresh) + &
                          & (1.- zx) *  zy  *  zz  * wb_dust_src%wgl%pGustVal(iWind,iZ+1,iThresh+1) + &
                          &   zx  *  zy * (1. - zz)* wb_dust_src%wgl%pGustVal(iWind+1,iZ+1,iThresh) + &
                          &   zx * zy * zz * wb_dust_src%wgl%pGustVal(iWind+1, iZ+1,iThresh+1)) / fSum
            else
!              if(uStarSaturated > u_star_min)then
!                uStarGustEff = (uStarSaturated-u_star_min) * (uStarSaturated + u_star_min)**2
              if(pU_star(iMeteo) > u_star_min)then
                uStarGustEff = (pU_star(iMeteo)-u_star_min) * (pU_star(iMeteo) + u_star_min)**2
              else
                uStarGustEff = 0.
              endif
            endif   ! ifWindGust

#ifdef DEBUG_DUST
fArDump(ix,iy,5) = uStarGustEff
#endif

            fSaltation = 2.6 * density_air_288K / g * &
                       & (1. - pSnowWEQDepth(iMeteo)/fSnowWEQDepthThreshold) ** 2 * & ! get rid of snowy sand
                       & uStarGustEff

if(fSaltation < 0.)then
 call set_error('Negative saltation','compute_emission_for_wb_dust')
 call msg('fWindSpeed, fWindSPeedCritical',fWindSpeed, fTmp)
 call msg('iThresh,iZ',iThresh,iZ)
 call msg('Owen,iWind',du_star_owen,iWind)
 call msg('zx,zy',zx,zy)
 call msg('zz,fSum',zz,fSum)
 call msg('Gust,Saltation',wb_dust_src%wgl%pGustVal(iWind,iZ, iThresh),fSaltation)
 ifTalking = .true.
 iCountTalking = 0
endif


            if(fSaltation <= 0) cycle  ! even with gust saltation can be negligible
!iArStatistics(5) = iArStatistics(5) + 1

            !
            ! Total flux contains the ratio of horizontal saltation flux to vertical fine-dust flux.
            ! According to MaB95, it is ~10**(linear of clay fraction). However, those folks argue
            ! that it is not exactly so. Their function is based on Gillette, who actually showed
            ! that very high clay fraction means no emission. Max ratio is around clay_fr=0.2.
            ! Hence, we took second-order fit. Optimum has to be checked from the model fit.
            !
            fFluxTotal = fEntrainmentScaling * &
                       & fSaltation * timestep_sec * &
! Original MaB95 with limit at clay_fr=0.2 (Zender argues that it is needed).
!                       & 100. * 10.**(13.4 * min(pClayMassFraction(iMeteo),0.2) - 6.0) * &
! Second-order fit with max at 0.2
                       & 100. * 10.**(-50. * (pClayMassFraction(iMeteo)**2) + &
                                    & 25. * pClayMassFraction(iMeteo) - 6.5) * &
                       & pErodibleFraction(iDisp) * &          ! area fraction of desert: source mask
                       & pAlluvialSedimentIndex(iDisp) * &     ! availability of small particles
                       & (fLAI_critical - max(0.,pLAI(iMeteo)))

#ifdef DEBUG_DUST
fArDump(ix,iy,6) = pU_star(iMeteo)   !+ du_star_owen
fArDump(ix,iy,7) = fSaltation
fArDump(ix,iy,8) = fFluxTotal
#endif


if(fFluxTotal > 1.0e-2)then
  ifTalking = .true.
  call msg('==================Strong flux: saltation, total:',fSaltation, fFluxTotal)
  call msg('ix,iy:', ix,iy)
  call msg('lon,lat:', ptrLonMeteo(iMeteo), ptrLatMeteo(iMeteo))
  call msg('iMeteo, iDisp',iMeteo,iDisp)
  call msg('ustar:', pU_star(iMeteo))
  call msg('ustar_min, sand mass fraction',u_star_min,pSandMassFract(iMeteo))
  call msg('Clay, its function:', pClayMassFraction(iMeteo), 100. * 10.**(13.4 * pClayMassFraction(iMeteo) - 6.0))
  call msg('Erodible fraction, alluvial sediment', pErodibleFraction(iDisp), pAlluvialSedimentIndex(iDisp))
  call msg('')
  iCount = 0
  iCountTalking = 0
elseif(fFluxTotal < 0.)then
  ifTalking = .true.
  call msg('****************************Negative flux: saltation, total:',fSaltation, fFluxTotal)
  call msg('ix,iy:', ix,iy)
  call msg('lon,lat:', ptrLonMeteo(iMeteo), ptrLatMeteo(iMeteo))
  call msg('iMeteo, iDisp',iMeteo,iDisp)
  call msg('ustar:', pU_star(iMeteo))
  call msg('ustar_min, sand mass fraction',u_star_min,pSandMassFract(iMeteo))
  call msg('Clay, its function:', pClayMassFraction(iMeteo), 100. * 10.**(13.4 * pClayMassFraction(iMeteo) - 6.0))
  call msg('Erodible fraction, alluvial sediment', pErodibleFraction(iDisp), pAlluvialSedimentIndex(iDisp))
  call msg('')
  iCount = 0
  iCountTalking = 0
elseif(iCountTalking < 10)then
  call msg('Normal saltation, total:',fSaltation, fFluxTotal)
  call msg('ix,iy:', ix,iy)
  call msg('ustar:', pU_star(iMeteo))
  call msg('ustar_min, sand mass fraction',u_star_min,pSandMassFract(iMeteo))
  call msg('Clay, its function:', pClayMassFraction(iMeteo), 100. * 10.**(13.4 * pClayMassFraction(iMeteo) - 6.0))
  call msg('Erodible fraction, alluvial sediment', pErodibleFraction(iDisp), pAlluvialSedimentIndex(iDisp))
  iCountTalking = iCountTalking + 1
endif
          case default
            call set_error('Unknown emission method','compute_emission_for_wb_dust')
            return
        end select
        !
        ! Fill-in the emission map, not forgetting the number to mass convertion where needed
        !
        do iLev = 1, wb_dust_src%nLevsDispVert
          !
          ! First do the check for the overlap: speed-up
          !
          if(abs(wb_dust_src%levFractDispVert(iLev)) < 1e-6)cycle  ! nothing for this dispersion layer

          fCellTotal = 0.0

          do iDust = 1, wb_dust_src%nSpecies
            if(wb_dust_src%ifNumberFlux(iDust))then
              !
              ! Number flux
              !
              fTmp = fFluxTotal * &                        ! Total flux
                   & wb_dust_src%fluxPerModeNbr(iDust) * &     ! Spectrum fractionation
                   & wb_dust_src%levFractDispVert(iLev) * &    ! Vertical overlap
                   & ptrXSzDisp(iDisp) * ptrYSzDisp(iDisp)     ! m2 size

if(ifTalking)call msg('Number flux final, flux before scaling:',fTmp, fFluxTotal)

            else
              !
              ! Volume flux, converted to mass
              !
              fTmp = fFluxTotal * fPartDensity(iDust) * &
                   & wb_dust_src%fluxPerModeVol(iDust) * &     ! Spectrum fractionation
                   & wb_dust_src%levFractDispVert(iLev) * &    ! Vertical overlap
                   & ptrXSzDisp(iDisp) * ptrYSzDisp(iDisp)     ! m2

if(ifTalking)then
call msg('Volume flux final, flux before scaling:',fTmp, fFluxTotal)
iCountTalking = iCountTalking + 1
if(iCountTalking > 20)ifTalking = .false.
endif
            endif

            fCellTotal = fCellTotal + fTmp

            emisMap%arM(wb_dust_src%adaptor%iSp(iDust),wb_dust_src%id_nbr,iLev,ix,iy) = fTmp + &
                       & emisMap%arM(wb_dust_src%adaptor%iSp(iDust),wb_dust_src%id_nbr,iLev,ix,iy)
            fMassInjected(wb_dust_src%adaptor%iSp(iDust)) = &
                                        & fMassInjected(wb_dust_src%adaptor%iSp(iDust)) + fTmp

!call msg('WB-DUST species index and emission:' + fu_str(iDust), &
!             & wb_dust_src%adaptor%iSp(iDust), fMassInjected(wb_dust_src%adaptor%iSp(iDust)))
if(fMassInjected(wb_dust_src%adaptor%iSp(iDust)) > 1e30)then
  call set_error('Huge emission','compute_emission_for_wb_dust')
  call unset_error('compute_emission_for_wb_dust')
endif


            if (ifSpeciesMoment) then
!              mapCoordX%arm(iSpeciesEmis,a_src%id_nbr,ilev, ix,iy) = &
!                          & mapCoordX%arm(iSpeciesEmis, a_src%id_nbr, ilev, ix,iy) + &
!                          & fMassTmp * 0.0
!              mapCoordY%arm(iSpeciesEmis, a_src%id_nbr, ilev, ix,iy) = &
!                          & mapCoordY%arm(iSpeciesEmis, a_src%id_nbr, ilev, ix,iy) + &
!                          & fMassTmp * 0.0
              mapCoordZ%arm(wb_dust_src%adaptor%iSp(iDust),wb_dust_src%id_nbr, ilev, ix,iy) = &
                   & mapCoordZ%arm(wb_dust_src%adaptor%iSp(iDust),wb_dust_src%id_nbr, ilev, ix,iy) + &
                   & fTmp * wb_dust_src%fzDisp(iLev)
            end if

            emisMap%ifColumnValid(wb_dust_src%id_nbr,ix,iy) = .true.
            emisMap%ifGridValid(iLev,wb_dust_src%id_nbr) = .true.

          end do      ! nSpecies
          if (.not. ifSpeciesMoment) then
!            mapCoordX%arM(1,wbdust_src%id_nbr, iLev, ix, iy) = 0.0 * fCellTotal + &
!                                                & mapCoordX%arM(1,wbdust_src%id_nbr, iLev, ix, iy)
!            mapCoordY%arM(1,wbdust_src%id_nbr, iLev, ix, iy) = 0.0 * fCellTotal + &
!                                                & mapCoordY%arM(1,wbdust_src%id_nbr, iLev, ix, iy)
            mapCoordZ%arM(1,wb_dust_src%id_nbr, iLev, ix, iy) = &
                                                & wb_dust_src%fzDisp(iLev) * fCellTotal + &
                                                & mapCoordZ%arM(1,wb_dust_src%id_nbr, iLev, ix, iy)
          end if
        end do      ! iLev
      end do     ! nx_dispersion
    end do    ! ny_dispersion


#ifdef DEBUG_DUST
!if(fu_min(now) == 0)then
  do iDust = 1, 8
    write(uDump,rec=iRec)fArDump(1:nx_dispersion,1:ny_dispersion,iDust)  !((fArDump(ix,iy,iDust),ix=1,nx_dispersion),iy=1,ny_dispersion)
    iRec = iRec+1
    arTmp(iDust) = sum(fArDump(1:nx_dispersion,1:ny_dispersion,iDust)) / (nx_dispersion * ny_dispersion)
  enddo
  call msg('Dump at:' + fu_str(now),(/arTmp(1:8)/))
!endif
#endif

    call free_work_array(fPartDensity)

do iDust = 1, wb_dust_src%nSpecies
call msg('WB-DUST species index and emission:' + fu_str(iDust), &
             & wb_dust_src%adaptor%iSp(iDust), fMassInjected(wb_dust_src%adaptor%iSp(iDust)))
if(fMassInjected(wb_dust_src%adaptor%iSp(iDust)) < 0.)then
  call set_error('Negative emission','compute_emission_for_wb_dust')
endif
end do
call msg('Out of: ' + fu_str(iArStatistics(1)) + &
       & ', passed through src mask:' + fu_str(iArStatistics(2)) + &
       & ', snow:' + fu_str(iArStatistics(3)) + &
       & ', vegetation:' + fu_str(iArStatistics(4)) + &
       & ', rough soil:' + fu_str(iArStatistics(5)) + &
       & ', uStar dry soil:' + fu_str(iArStatistics(6)) + &
       & ', uStar wet soil:' + fu_str(iArStatistics(7)))

  end subroutine compute_emission_for_wb_dust


  !********************************************************************************************

  logical function fu_wb_dust_emis_owned_quantity(wb_dust_src, quantity)
    !
    ! Check whether the specific quantity is handled by the emission sources exclusively
    ! Such quantities may be related to emission species somewhere else
    !
    implicit none

    ! Imported parameters
    type(silam_wind_blown_dust_source), intent(in) :: wb_dust_src
    integer, intent(in) :: quantity
    !
    ! The wind blown dust source does not have own quantities yet
    !
    select case(quantity)
      case default
        fu_wb_dust_emis_owned_quantity = .false.
    end select
  end function fu_wb_dust_emis_owned_quantity


  !*****************************************************************

  integer function fu_source_id_nbr_of_wb_dust_src(wb_dust_src)
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
    type(silam_wind_blown_dust_source), intent(in) :: wb_dust_src

    ! Stupidity check
    if(.not.(wb_dust_src%defined == silja_false))then
      call set_error('Undefined source given','fu_source_id_nbr_of_wb_dust_src')
      return
    endif
    fu_source_id_nbr_of_wb_dust_src = wb_dust_src%id_nbr

  end function fu_source_id_nbr_of_wb_dust_src


  !*************************************************************************

  integer function fu_source_nbr_of_wb_dust_src(wb_dust_src)
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
    type(silam_wind_blown_dust_source), intent(in) :: wb_dust_src

    ! Stupidity check
    if(.not.(wb_dust_src%defined == silja_false))then
      fu_source_nbr_of_wb_dust_src = wb_dust_src%src_nbr
    else
      fu_source_nbr_of_wb_dust_src = int_missing
      call set_error('Undefined source given','fu_source_nbr_of_wb_dust_src')
      return
    endif

  end function fu_source_nbr_of_wb_dust_src


  !*************************************************************************

  subroutine typical_species_cnc_wb_dust(srcWBDust, species, nSpecies, arConc)
    !
    ! Guesses a typical level of concentration and divides it with the given accuracy factor
    !
    implicit none

    ! Imported parameters
    type(silam_wind_blown_dust_source), intent(in) :: srcWBDust
    type(silam_species), dimension(:), pointer :: species
    integer, intent(out) :: nSpecies
    real, dimension(:), pointer :: arConc

    ! Local variables
    integer :: iSpecies
    real :: fTotalFluxNbr, fTotalFluxVol, fMassMeanDiam
    type(Taerosol_mode) :: mode

    real, parameter :: fTypicalMassCnc = 1.0e-8  ! ten microgram

    species => srcWBDust%species
    nSpecies = srcWBDust%nSpecies

    !
    ! Procedure: compute the flux for each of the dust species, and the total flux
    ! for the whole size range. The total typical dust concentration is a few tens of micrograms
    ! Then the fraction of the mass flux decides on the fraction of the concentration.
    ! Strictly speaking, this is not true due to varying deposition but here we do not have too
    ! many options.
    !
    ! The total flux
    !
    call set_aerosol_mode(mode, '', &                          ! mode, chNm, 
                        & 1.e-8, 1.e-5, &                      ! fp1, fp2, &
                        & 1.e-6, 1.e3, fixed_diameter_flag, 1) ! mass_mean_d, dens, distr_type, solubil
    if(error)return

    call wind_blown_dust_flux4mode(mode, srcWBDust%iSpectrumType, &
                                 & fTotalFluxNbr, fTotalFluxVol, fMassMeanDiam)
    !
    ! With the total flux known, get the typical concentrations
    !
    do iSpecies = 1, srcWBDust%nSpecies
      arConc(iSpecies) = fTypicalMassCnc * &
                       & min(1.0, srcWBDust%fluxPerModeVol(iSpecies) / fTotalFluxVol)
    end do

  end subroutine typical_species_cnc_wb_dust


  !*************************************************************************

  function fu_wb_dust_source_name(wb_dust_src)result(chNm)
    implicit none
    type(silam_wind_blown_dust_source), intent(in) :: wb_dust_src
    character(len=clen) :: chNm
    chNm = wb_dust_src%src_nm
  end function fu_wb_dust_source_name


  !*************************************************************************

  subroutine report_wb_dust_src(wb_dust_src)
    implicit none
    type(silam_wind_blown_dust_source), intent(in) :: wb_dust_src
    integer :: iSpecies
    call msg('Wind-blown dust source'+wb_dust_src%src_nm)
    do iSpecies = 1, wb_dust_src%nSpecies
      call report(wb_dust_src%species(iSpecies))
    end do
    call msg('')
    call msg('Reporting the dust source settings:')
    call msg('fEntrainmentScaling',fEntrainmentScaling)
    call msg('fSrfRoughSandOnly',fSrfRoughSandOnly)
    call msg('fDragPartitioningSandOnly = 1.0/log(0.35*(0.1/fSrfRoughSandOnly)**0.8)',fDragPartitioningSandOnly)
    call msg('fLAI_critical',fLAI_critical)
    call msg('u_star_excess_saturation =',u_star_excess_saturation)
  end subroutine report_wb_dust_src


  !*************************************************************************
  !*************************************************************************
  !
  ! Private functions computing the wind-blown dust emission fluxes
  !
  !*************************************************************************
  !*************************************************************************


  subroutine wind_blown_dust_flux4mode(aerMode, &        ! definition of the spectrum band
                                     & iSpectrumType, &
                                     & fNbrFlux, fVolFlux, fMassMeanDiam)
    !
    ! Computes the number emission flux (shape function) of the wind-blown dust 
    ! for the pre-defined size spectrum.
    !
    implicit none

    ! Imported parameters
    type(Taerosol_mode), intent(inout) :: aerMode
    integer, intent(in) :: iSpectrumType
    real, intent(out) :: fNbrFlux, fVolFlux, fMassMeanDiam

    ! Local variables
    real :: fMinD, fMaxD, fD, fValue, fValue_prev, fIntegrStep, fTmp, factor, &
          & fFluxDens, fFluxDens_prev, fD_um, fIntegrStep_um, fTmpNbr, fTmpVol
    type(Taerosol_mode) :: aerModeTmp
    integer :: iMode, i

    real :: fTotalNumber, fTotalVolume
    !
    ! Parameters for all distribution descriptions
    !
    real, parameter :: fAbsoluteMin = 2.0e-9  ! 2 nm
    real, parameter :: fAbsoluteMax = 1.0e-4  ! 100 um
    !
    ! Parameters for 4-mode spectrum representation (dissertation B.Weinzierl, 2007)
    !
    real, dimension(4), parameter :: nbr_fract_4_modes =     (/0.90455,  0.06552,  0.02842,  0.00151/)
    real, dimension(4), parameter :: nbr_mean_diam_4_modes =(/0.075e-6, 0.36e-6, 0.98e-6, 5.66e-6/)
    real, dimension(4), parameter :: nbr_stdev_4_modes =     (/1.9,       1.6,     2.0,     1.7/)

!    real, dimension(4), parameter :: nbr_fract_4_modes =     (/0.5,  0.5,  0.0,  0.0/)

    type(silam_sp) :: sp, sp2

sp%sp => fu_work_string()
sp2%sp => fu_work_string()
if(error)return

    !
    ! Algorithm: scan each size class and make an integration from min to max diameter of the
    ! number flux. Mean diameter can be then computed too and set to the aerosol back.
    !
    ! Get the parameters of the mode
    !
    fMinD = fu_min_d(aerMode)
    fMaxD = fu_max_d(aerMode)

    if(fMaxD <= fAbsoluteMin .or. fMinD >= fAbsoluteMax) then
      !
      ! The mode range is outside the known range, cannot say anything
      !
      fNbrFlux = 0.
      fVolFlux = 0.
      fMassMeanDiam = 0.5*(fMinD + fMaxD)

    else
      !
      ! The range given has useful part
      !
      if(fMinD < fAbsoluteMin) fMinD = fAbsoluteMin
      if(fMaxD > fAbsoluteMax) fMaxD = fAbsoluteMax

!call msg('Bin requested with the borders, um:', fMinD*1e6, fMaxD*1e6)

sp%sp = 'min:' + fu_str(fMinD*1e6) + ', max:' + fu_str(fMaxD*1e6) + ', number:'
      !
      ! The selector of the dust shape formulations
      !
      select case(iSpectrumType)

        case(lognorm_4_modes_dust_flag)
          !
          ! After Weinzierl dissertation: 4 lognormal modes
          !
call set_error('Analytical approximation does not work. Use numeric','wind_blown_dust_flux4mode')
return
          fNbrFlux = 0.0
          fVolFlux = 0.0
          fMassMeanDiam = 0.0
          fTotalNumber = 0.0
          fTotalVolume = 0.0
          do iMode = 1, 4

            call set_aerosol_mode(aerModeTmp, '', &                         ! mode, chNm, 
                                & nbr_mean_diam_4_modes(iMode), &           ! fp1
                                & nbr_stdev_4_modes(iMode), &  ! fp2
                                & nbr_mean_diam_4_modes(iMode), 2.6e3, &   ! mass_mean_d, dens, 
                                & lognormal_flag, 1)                       ! distr_type, solubil
            if(error)return
            !
            ! Calculate the total-spectrum contribution
            !
            fTotalNumber = fTotalNumber + &
                                & nbr_fract_4_modes(iMode) * &
                                & fu_integrate_number(fAbsoluteMin, fAbsoluteMax, aerModeTmp)
            fTotalVolume = fTotalVolume + &
                                & nbr_fract_4_modes(iMode) * &
                                & fu_integrate_volume(fAbsoluteMin, fAbsoluteMax, aerModeTmp) * &
                                & (fu_mass_mean_diam(fAbsoluteMin, fAbsoluteMax, aerModeTmp))**3
call msg('Cumulative total number and volume, full range:', fTotalNumber, fTotalVolume)
            !
            ! Number-flux fraction is the mode fraction in the distribution * 
            !                             fraction of the mode falling in the output mode
            !
            fTmp = nbr_fract_4_modes(iMode) * fu_integrate_number(fMinD, fMaxD, aerModeTmp)
            fNbrFlux = fNbrFlux + fTmp
sp%sp = sp%sp + ',' + fu_str(real(fTmp))
            !
            ! Volume-flux fraction is the fraction of the mode in the volume
            !
            fTmp = nbr_fract_4_modes(iMode) * fu_integrate_volume(fMinD, fMaxD, aerModeTmp) * &
                 & (fu_mass_mean_diam(fMinD, fMaxD, aerModeTmp))**3
sp2%sp = sp2%sp + ',' + fu_str(real(fTmp))
            fVolFlux = fVolFlux + fTmp
            fMassMeanDiam = fMassMeanDiam + fTmp * fu_mass_mean_diam(fMinD, fMaxD, aerModeTmp)
          end do
          
call msg(sp%sp + ',   volume:' + sp2%sp)

          fNbrFlux = fNbrFlux / fTotalNumber
          fVolFlux = fVolFlux / fTotalVolume

call msg('Total number and volume fluxes:', fTotalNumber, fTotalVolume)
call msg('Normalised number and volume flux for the bin:', fNbrFlux, fVolFlux)
call msg('')

          fMassMeanDiam = fMassMeanDiam / fVolFlux
          
          
        case(Kok_Brittle_fragm_theory_flag)
          !
          ! Kok etal, 2011, (Physics of wind-blown sand and dust) derived fragmentation-theory
          ! based formulas:
          ! dN/dlnD = 1/(cn*D^2)*(1+erf(ln(D/Ds)/(sqrt(2)*ln(Sigma))*ln(-D/L)))
          ! dV/dlnD = D/(cv)*(1+erf(ln(D/Ds)/(sqrt(2)*ln(Sigma))*ln(-D/L)))
          ! cn=0.9539um-2, cv=12.62, Ds=3.4um, SigmaS=3.0, L=12um
          ! That has to be integrated numerically
          !
          ! Get the integration step
          !
!          fIntegrStep = 0.00001*(fMaxD - fMinD)  ! At least 1000 steps over the range
!          fIntegrStep = 0.000001*(fAbsoluteMax - fAbsoluteMin)  ! At least 1000 steps over the range
          factor = (fAbsoluteMax / fAbsoluteMin) ** 0.0001 - 1.
          
          fD = fAbsoluteMin !fMinD
          fNbrFlux = 0.
          fVolFlux = 0.
          fTotalNumber = 0.0
          fTotalVolume = 0.0
          fMassMeanDiam = 0.

          fFluxDens = fu_shape_function_Kok_Brittle_nbr(fD*1e6)
          !
          ! Integration goes over the whole defined range, simultaneously separating the 
          ! given node. Awfull from the point of view of efficiency but fine here since the
          ! sub is called only N-dust-bins times at the initialization phase.
          !
          do while(fD < fAbsoluteMax) !fMaxD)
!            fIntegrStep = min(fIntegrStep, fMaxD - fD)
!            fIntegrStep = min(fIntegrStep, fAbsoluteMax - fD)
!            fD = fD + fIntegrStep

            fIntegrStep = fD * factor

            if(fD < fMinD .and. fD + fIntegrStep > fMinD)then
              fD = fMinD
            elseif(fD < fMaxD .and. fD + fIntegrStep > fMaxD)then
              fD = fMaxD
            else
              fD = fD + fIntegrStep
            endif
            
            fFluxDens_prev = fFluxDens
            fD_um = fD * 1e6
            fIntegrStep_um = fIntegrStep * 1e6
            !
            ! Number-flux density
            !
            fFluxDens = fu_shape_function_Kok_Brittle_nbr(fD_um)
            !
            ! Integrate the number and volume flux. 
            !
            fTmpNbr = fIntegrStep_um * 0.5 * (fFluxDens_prev + fFluxDens)  ! for number flux
            fTmpVol = fIntegrStep_um * 0.5 * 0.5235987756 * (fFluxDens_prev * (fD_um-fIntegrStep_um)**3 + &
                                                           & fFluxDens * fD_um**3)
            fTotalNumber = fTotalNumber + fTmpNbr
            fTotalVolume = fTotalVolume + fTmpVol
            if(fD >= fMinD .and. fD <= fMaxD)then
              fNbrFlux = fNbrFlux + fTmpNbr
              fVolFlux = fVolFlux + fTmpVol
              fMassMeanDiam = fMassMeanDiam + fTmpVol * (fD_um-0.5*fIntegrStep_um)
            endif
          end do  ! cycle over diameter
          
!          fMassMeanDiam = 1.e-6 * fMassMeanDiam / fVolFlux
!          fVolFlux = fVolFlux * 1.0e-18  ! from um^3 to m^3

!call msg('fu_shape_function_Kok_Brittle_nbr: fD range min-max, totalNbr, totalVol, binNbr, binVol', &
!    & (/fMinD, fMaxD, fTotalNumber, fTotalVolume, fNbrFlux, fVolFlux/))

          fMassMeanDiam = fMassMeanDiam * 1e-6 / fVolFlux
          if(fMassMeanDiam < fMinD .or. fMassMeanDiam > fMaxD)then
            call msg('Wrong Kok mean diameter. Min, max, mean:', (/fMinD, fMaxD, fMassMeanDiam/))
            call msg_warning('Resetting mean diameter','wind_blown_dust_flux4mode')
            fMassMeanDiam = sqrt(fMinD * fMaxD)
          endif
          fNbrFlux = fNbrFlux / fTotalNumber
          fVolFlux = fVolFlux / fTotalVolume
        
        case(lognorm_4_modes_dust_numeric_flag)
          !
          ! After Weinzierl dissertation: 4 lognormal modes
          !
!          fIntegrStep = 0.0001*(fMaxD - fMinD)  ! At least 1000 steps over the range
!          fIntegrStep = 0.000001*(fAbsoluteMax - fAbsoluteMin)  ! At least 1000 steps over the range
          factor = (fAbsoluteMax / fAbsoluteMin) ** 0.0001 - 1.

          fD = fAbsoluteMin  !fMinD
          fNbrFlux = 0.
          fVolFlux = 0.
          fTotalNumber = 0.0
          fTotalVolume = 0.0
          fMassMeanDiam = 0.
          
          fFluxDens = 0.0
          do i = 1, 4
            fFluxDens = fFluxDens + nbr_fract_4_modes(i) * &
                                  & fu_lognorm(fD, nbr_mean_diam_4_modes(i), nbr_stdev_4_modes(i))
          end do

          do while(fD < fAbsoluteMax) !fMaxD)
!            fIntegrStep = min(fIntegrStep, fMaxD - fD)
!            fD = fD + fIntegrStep

            fIntegrStep = fD * factor

            if(fD < fMinD .and. fD + fIntegrStep > fMinD)then
              fD = fMinD
            elseif(fD < fMaxD .and. fD + fIntegrStep > fMaxD)then
              fD = fMaxD
            else
              fD = fD + fIntegrStep
            endif
            
            fFluxDens_prev = fFluxDens

            !
            ! Number-flux density
            !
            fFluxDens = 0.0
            do i = 1, 4
              fFluxDens = fFluxDens + nbr_fract_4_modes(i) * &
                                    & fu_lognorm(fD, nbr_mean_diam_4_modes(i), nbr_stdev_4_modes(i))
            end do
            !
            ! Integrate the number and volume flux. 
            !
            fTmpNbr = fIntegrStep * 0.5 * (fFluxDens_prev + fFluxDens)  ! for number flux
            fTmpVol = fIntegrStep * 0.5 * 0.5235987756 * (fFluxDens_prev * (fD-fIntegrStep)**3 + &
                                                        & fFluxDens * fD**3)
            fTotalNumber = fTotalNumber + fTmpNbr
            fTotalVolume = fTotalVolume + fTmpVol
            if(fD >= fMinD .and. fD <= fMaxD)then
              fNbrFlux = fNbrFlux + fTmpNbr
              fVolFlux = fVolFlux + fTmpVol
              fMassMeanDiam = fMassMeanDiam + fTmpVol * (fD-0.5*fIntegrStep)
            endif

          end do  ! cycle over diameter
!call msg('lognorm_4_modes_dust_numeric_flag: fD range min-max, totalNbr, totalVol, binNbr, binVol', &
!    & (/fMinD, fMaxD, fTotalNumber, fTotalVolume, fNbrFlux, fVolFlux/))
          
          fMassMeanDiam = fMassMeanDiam / fVolFlux
          if(fMassMeanDiam < fMinD .or. fMassMeanDiam > fMaxD)then
            call msg('Wrong 4-mode mean diameter. Min, max, mean:', (/fMinD, fMaxD, fMassMeanDiam/))
            call msg_warning('Resetting mean diameter','wind_blown_dust_flux4mode')
            fMassMeanDiam = sqrt(fMinD * fMaxD)
          endif
          fNbrFlux = fNbrFlux / fTotalNumber
          fVolFlux = fVolFlux / fTotalVolume
          
        case default
          call set_error('Unknown type of dust spectrum:'+fu_str(iSPectrumTYpe), &
                       & 'wind_blown_dust_flux4mode')

        end select    ! selector of dust distribution function
    
    endif  ! useful part of the range

    
call free_work_array(sp%sp)
call free_work_array(sp2%sp)
    
    
    CONTAINS

      real function fu_lognorm(fD, fDmed, fSigma)            
        implicit none
        ! Imported parameters
        real, intent(in) :: fD, fDmed, fSigma
        fu_lognorm = 1D0 / (fD * log(fSigma)*dsqrt(2.0D0*Pi)) * exp(-((log(fD)-log(fDmed))/log(fSigma))**2 / 2.0D0)
!call msg('lognorm: D, Dmed,sigma, val',(/fD, fDmed, fSigma, fu_lognorm/))
      end function fu_lognorm
                                     
      real function fu_shape_function_Kok_Brittle_nbr(fD_um)result(fRate)
        ! Kok etal, 2011, (Physics of wind-blown sand and dust) derived fragmentation-theory
        ! based formulas:
        ! dN/dlnD = 1/(cn*D^2)*(1+erf(ln(D/Ds)/(sqrt(2)*ln(Sigma))*ln(-D/L)))
        ! dV/dlnD = D/(cv)*(1+erf(ln(D/Ds)/(sqrt(2)*ln(Sigma))*ln(-D/L)))
        ! cn=0.9539um-2, cv=12.62, Ds=3.4um, SigmaS=3.0, L=12um
        implicit none
        ! Imported parameter
        real, intent(in) :: fD_um
        
        ! Local variables
        real(kind=r8k), parameter :: cn = 0.9539, & ! um-2 
                                   & cv=12.62, &
                                   & Ds=3.4,   &    !um, 
                                   & SigmaS=3.0, ln_sigma_s2 = 1.5536723984242, &  !ln(Sigma) * sqrt(2)
                                   & L=12.0         ! um
        fRate = (1.0D0 + ERF(log(fD_um / Ds) / ln_Sigma_s2)) * exp(-(fD_um / L)**3) / (cn * fD_um**3)
      
      end function fu_shape_function_Kok_Brittle_nbr
    
  end subroutine wind_blown_dust_flux4mode
  
END MODULE source_terms_wind_blown_dust

