MODULE source_terms_dms
  !
  ! A parameterisation for biogenic emission of dimethyl sulfide (DMS) from oceans. 
  ! 
  ! The method is same as described by Kettle et al. (2000) or Tarrason et al. (1995): the
  ! emission is controlled by transfer resistance below the water surface, and bulk
  ! gas-phase concentration is assumed much lower than the Henry's law equilibrium. Then 
  ! 
  ! F = k * cw,
  ! 
  ! where k is the transfer coefficient and cw is concentration in water.
  ! Here, k is evaluated with the formula of Liss and Merilvat (1986) and cw is taken from 
  ! a monthly climatology. 
  ! 
  ! The source is initialized from a setup file, which needs to define
  ! - the species emitted and yield per DMS (allows emitting directly to SO2)
  ! - the path to DMS climatology
  ! - the path to sea mask (similar to the sea salt source).
  !
  ! References: 
  !
  ! Kettle, A.J., Andreae, M.O., 2000. Flux of dimethylsulfide from the oceans
  ! : A comparison of updated data sets and flux models. J. Geophys. Res. 105,
  ! 26,793-26,808. doi:10.1029/2000JD900252
  !
  ! Tarrason, L., Turner, S., Fløisand, I., 1995. Estimation of seasonal dimethyl sulphide
  ! fluxes over the North Atlantic Ocean and their contribution to European pollution
  ! levels. J. Geophys. Res.  Atmos. 100, 11623–11639.
  
  use source_terms_time_params

  implicit none

  private

  public test_dmssrc_no_init
  public test_init_dmssrc

  public fill_next_source_from_namelist
  interface fill_next_source_from_namelist
    module procedure fill_dms_src_from_namelist
  end interface

  public  add_inventory_dmssrc

  public add_input_needs
  interface add_input_needs
    module procedure add_input_needs_dmssrc
  end interface

  public create_source_containing_grid
  interface create_source_containing_grid
    module procedure create_src_cont_grd_dmssrc
  end interface

  public source_2_second_grid
  interface source_2_second_grid
    module procedure project_dmssrc_second_grd
  end interface

  public link_source_to_species
  interface link_source_to_species
    module procedure link_dmssrc_to_species
  end interface

  public fu_source_nbr
  interface fu_source_nbr
    module procedure fu_source_nbr_dmssrc
  end interface

  public fu_source_id_nbr
  interface fu_source_id_nbr
    module procedure fu_source_id_nbr_dmssrc
  end interface

  public fu_name
  interface fu_name
    module procedure fu_name_dmssrc 
  end interface

  public typical_species_conc
  interface typical_species_conc
    module procedure typical_species_cnc_dms
  end interface

  public report
  interface report
    module procedure report_dmssrc
  end interface
  
  public defined
  interface defined
     module procedure defined_dmssrc
  end interface

  public fu_dms_emis_owned_quantity
  public compute_emission_for_dms
  public init_emission_dms
  public fill_dms_src_from_namelist
  public reserve_dms_source

  public silam_dms_source
  TYPE silam_dms_source
    PRIVATE
    ! Source grid. Normally dispersion_grid.
    type(silja_grid) :: grid
    ! Emitted substance
    character(len=substNmLen) :: emis_subst_name = ''
    ! The emitted species determined from emis_subst_name. Only one.
    ! array needed because of the way typical_cnc works:
    type(silam_species), dimension(:), pointer :: species_list => null() ! note: initialization matters
    integer :: ind_emis = int_missing ! in emission species
    
    integer :: src_nbr = int_missing
    integer :: id_nbr = int_missing
    character(len=clen) :: src_name = 'marine_dms'
    character(len=clen) :: src_sector_name = 'biogenic'
    
    real :: yield = real_missing ! mol/mol
    
    ! dms map, in the grid given above: (x, y, month)
    real, dimension(:,:,:), pointer :: map_dms => null()
    ! fraction of water, in the grid given above (x, y)
    real, dimension(:,:), pointer :: water_mask => null()

    character(len=fnlen) :: map_dms_filename, water_mask_filename = ''
    type(silam_fformat) :: map_dms_fileformat, water_mask_fileformat

    ! The year that corresponds to the input file for reading the monthly DMS fields.
    integer :: map_dms_refyear = 2011
    ! A threshold SST below which no emission is evaluated.
    real :: sst_min = 271.15 ! K, -2 C, assumed frozen below that
    ! Not implemented yet.
    logical :: use_ice_fract = .false.
    
    ! Level fractions. Similar to sea salt source.
    real, dimension(:), pointer :: levFractDispVert => null(), fzDisp => null()
    integer :: nLevsDispVert = int_missing
    
    logical :: defined = .false.
 END TYPE silam_dms_source

 public dms_src_ptr
 type dms_src_ptr
    type(silam_dms_source) :: dms_src
 end type dms_src_ptr
 
CONTAINS


  !*********************************************************************

  subroutine fill_dms_src_from_namelist(nlptr, dmssrc)
    implicit none

    ! Imported parameters
    type(Tsilam_namelist), pointer :: nlptr
    type(silam_dms_source), intent(inout) :: dmssrc

    character(len=*), parameter :: subname = 'fill_from_namelist_dmssrc'
    character(len=fnlen) :: filename, fileformat, nl_content
    real :: nl_content_real
    integer :: nl_content_int
    

    if (fu_fails(associated(nlptr), 'nlptr not associated', subname)) return

    ! No defaults set here, the defaults paramaters are in the data type definition.

    nl_content = fu_content(nlptr, 'source_name')
    if (nl_content /= '') dmssrc%src_name = nl_content
    nl_content = fu_content(nlptr, 'source_sector_name')
    if (nl_content /= '') dmssrc%src_sector_name = nl_content
    
    nl_content = fu_content(nlptr, 'dms_map_filename')
    if (fu_fails(nl_content /= '', 'Missing dms_map_filename', subname)) return
    fileformat = nl_content(1:index(nl_content, ' ')-1)
    filename = nl_content(index(nl_content, ' ')+1:)
    dmssrc%map_dms_filename = fu_process_filepath(filename, convert_slashes=.true., must_exist=.true.)
    if (error) return
    dmssrc%map_dms_fileformat = fu_input_file_format(fileformat)
    if (error) return

    nl_content = fu_content(nlptr, 'water_mask_filename')
    if (nl_content == '') nl_content = fu_content(nlptr, 'source_area_mask')
    if (fu_fails(nl_content /= '', 'Missing water_mask_filename', subname)) return
    fileformat = nl_content(1:index(nl_content, ' ')-1)
    filename = nl_content(index(nl_content, ' ')+1:)
    dmssrc%water_mask_filename = fu_process_filepath(filename, convert_slashes=.true., must_exist=.true.)
    if (error) return
    dmssrc%water_mask_fileformat = fu_input_file_format(fileformat)
    if (error) return
    
    nl_content_real = fu_content_real(nlptr, 'sst_min')
    if (.not. (nl_content_real .eps. real_missing)) then
      if (fu_fails(nl_content_real > 250 .and. nl_content_real < 300, 'Strange sst_min', subname)) return
      dmssrc%sst_min = nl_content_real
    end if
    
    nl_content_int = fu_content_int(nlptr, 'dms_map_refyear')
    if (nl_content_int /= int_missing) dmssrc%map_dms_refyear = nl_content_int
    
    nl_content = fu_content(nlptr, 'emitted_substance')
    if (fu_fails(nl_content /= '', 'Missing emitted_substance', subname)) return
    dmssrc%emis_subst_name = nl_content
    nl_content_real = fu_content_real(nlptr, 'yield')
    if (nl_content_real .eps. real_missing) then
      if (dmssrc%emis_subst_name == 'DMS') then 
        dmssrc%yield = 1.0
      else
        call set_error('If emitted_substance != DMS, must specify yield (mol/mol)', subname)
        return
      end if
    else
      if (fu_fails(nl_content_real > 0, 'Strange yield', subname)) return
      dmssrc%yield = nl_content_real
    end if

    call get_species(dmssrc)
    
  end subroutine fill_dms_src_from_namelist


  !*****************************************************************

  subroutine add_input_needs_dmssrc(dmssrc, q_met_dynamic, q_met_static, &
                                  & q_disp_dynamic, q_disp_static)
    !
    ! Returns input needs for the emission/transformation routines. 
    !
    implicit none

    ! Imported parameters
    type(silam_dms_source), intent(in) :: dmssrc
    integer, dimension(:), intent(inout) :: q_met_dynamic, q_met_static, &
                                          & q_disp_dynamic, q_disp_static
    
    integer :: size_new
    
    size_new = fu_merge_integer_to_array(u_10m_flag, q_met_dynamic)
    size_new = fu_merge_integer_to_array(v_10m_flag, q_met_dynamic)
    size_new = fu_merge_integer_to_array(water_surface_temp_flag, q_met_dynamic)
    
  end subroutine add_input_needs_dmssrc


  !**************************************************************************

  subroutine reserve_dms_source(dmssrc, &     ! Src to initialise
                              & iSrcNbr, &      ! Src number in the grand list
                              & iSrcIdNbr)      ! SrcID number
    !
    ! Initialises the source:
    ! - stores the reference information: source number and source ID number
    implicit none

    ! Imported parameters
    type(silam_dms_source), intent(inout) :: dmssrc
    integer, intent(in) :: iSrcNbr, iSrcIdNbr

    dmssrc%src_nbr = iSrcNbr
    dmssrc%id_nbr = iSrcIdNbr

  end subroutine reserve_dms_source


  !*********************************************************************

  subroutine init_emission_dms(dmssrc, grid)
    !
    ! Initializes the sea salt source term.
    ! The parameters are read from the given ini file and stored to the returned 
    ! source term. After the stuff has been read, the main FluxPerMode is to be generated.
    ! This configuration allows for several sea salt sources
    !
    implicit none

    type(silam_dms_source), intent(inout) :: dmssrc
    type(silja_grid), intent(in), optional :: grid ! normally dispersion, but can be overriden for easier testing

    integer :: stat, ix, iy, i1d, nx, ny, month
    type(silja_field_id) :: fid, idout
    type(silja_time) :: time

    character(len=*), parameter :: subname = 'init_emission_dmssrc'
    real, dimension(:), pointer :: values
    type(silam_material), pointer :: material
    type(silam_species) :: species_dms
    real :: max_water

    if (present(grid)) then
      if (fu_fails(defined(grid), 'Grid not defined', subname)) return
      dmssrc%grid = grid
    else
      dmssrc%grid = dispersion_grid
    end if
    call grid_dimensions(dmssrc%grid, nx, ny)
    allocate(dmssrc%map_dms(12, nx, ny), dmssrc%water_mask(nx, ny), stat=stat)
    if (fu_fails(stat == 0, 'Allocate failed', subname)) return
    
    call get_species(dmssrc)
    if (error) return

    ! Read the DMS monthly concentration, interpolate to dispersion grid and copy the values.
    ! 
    values => fu_work_array()
    call set_species(species_dms, fu_get_material_ptr('DMS'), in_gas_phase)
    if (error) return
    do month = 1, 12
      time = fu_set_time_utc(dmssrc%map_dms_refyear, month, 1, 0, 0, 0.0)
      ! need to request DMS - not the emission species which could be ->SO2
      fid = fu_set_field_id_simple(met_src_missing, concentration_flag, time, &
                                & level_missing, species_dms)
      if (error) return
      call get_input_field(dmssrc%map_dms_filename, dmssrc%map_dms_fileformat, &
                         & fid, values, dmssrc%grid, &
                         & iOutside = nearestPoint, &
                         & iAccuracy = 5, &
                         & wdr = wdr_missing, &
                         & ifAcceptSameMonth = .true., &
                         & idout=idout, fFillValue_ =0.) 
      if (error) return
      if (fu_fails(defined(idout), 'Failed to read DMS map', subname)) return
      do iy = 1, ny
        do ix = 1, nx
          i1d = (iy-1)*nx + ix
          dmssrc%map_dms(month, ix, iy) = values(i1d)
        end do
      end do
    end do

    ! Read the water mask
    ! 
    fid = fu_set_field_id_simple(met_src_missing, fraction_of_water_flag, time_missing, level_missing)
    if (error) return
    call get_input_field(dmssrc%water_mask_filename, dmssrc%water_mask_fileformat, &
                       & fid, values, dmssrc%grid, &
                       & iOutside = nearestPoint, &
                       & iAccuracy = 5, &
                       & wdr = wdr_missing, &
                       & ifAcceptSameMonth = .false., &
                       & idout=idout)
    if (error) return
    if (fu_fails(defined(idout), 'Failed to read water fraction', subname)) return
    do iy = 1, ny
      do ix = 1, nx
        i1d = (iy-1)*nx + ix
        dmssrc%water_mask(ix, iy) = values(i1d)
      end do
    end do
    
    call free_work_array(values)
    
    max_water = maxval(dmssrc%water_mask)
    if (fu_fails(max_water < 1.01, 'Too large water fraction', subname)) return
    call msg('Maximum fraction of water:', max_water)
    if (fu_fails(all(dmssrc%water_mask >= 0), 'Negative water fraction', subname)) return
    
    if (fu_fails(all(dmssrc%map_dms >= 0), 'Negative DMS concentration', subname)) return
    
    dmssrc%defined = .true.
    
  end subroutine init_emission_dms

  subroutine get_species(dmssrc)
    ! Set the species list according to the emitted substance.
    implicit none
    type(silam_dms_source), intent(inout) :: dmssrc

    character(len=*), parameter :: subname = 'get_species'
    integer :: stat
    type(silam_material), pointer :: material

    if (associated(dmssrc%species_list)) return
    allocate(dmssrc%species_list(1), stat=stat)
    if (fu_fails(stat == 0, 'Allocate failed', subname)) return
    if (fu_fails(dmssrc%emis_subst_name /= '', 'emis_subst_name is not defined', subname)) return
    material => fu_get_material_ptr(dmssrc%emis_subst_name)
    if (error) return
    call set_species(dmssrc%species_list(1), material, in_gas_phase)
    if (error) return
  end subroutine get_species

  !****************************************************************************

  subroutine add_inventory_dmssrc(dmssrc, species_list, nSpecies)
    !
    ! Get the species emitted by this source and make a list of
    ! them. The source must be initialized, of course.
    !
    implicit none
    type(silam_dms_source), intent(in) :: dmssrc
    type(silam_species), dimension(:), pointer :: species_list
    integer, intent(inout) :: nSpecies
    
    character(len=*), parameter :: subname = 'add_inventory_dmssrc'
    
    if (fu_fails(associated(dmssrc%species_list), 'dmssrc%species_list not associated', subname)) return
    call addSpecies(species_list, nSpecies, dmssrc%species_list, 1)

  end subroutine add_inventory_dmssrc


  !*******************************************************************************
  
  subroutine create_src_cont_grd_dmssrc(dmssrc, grid_template, ifVerbose, ifMinimal, ifExtended)
    implicit none

    ! Imported parameters
    type(silam_dms_source), intent(in) :: dmssrc
    type(silja_grid), intent(inout) :: grid_template
    logical, intent(in) :: ifVerbose, ifMinimal
    logical, intent(out) :: ifExtended

    ifExtended = .false.
    
  end subroutine create_src_cont_grd_dmssrc


  !*****************************************************************

  subroutine project_dmssrc_second_grd(dmssrc, grid, vert_disp, vert_proj, iAccuracy)
    !
    ! There is nothing to project in terms of horizontal grid but the vertical has to be made.
    !
    implicit none

    ! Imported parameters
    type(silam_dms_source), intent(inout) :: dmssrc
    type(silja_grid), intent(in) :: grid
    type(silam_vertical), intent(in) :: vert_disp, vert_proj
    integer, intent(in) :: iAccuracy

    ! Local variables
    integer :: stat
    type(silam_vertical) :: vertTmp
    real, dimension(2) :: fractions
    character(len=*), parameter :: subname = 'project_dmssrc_second_grd'
    !
    ! Now re-project the vertical grid of the source to the given vertical
    ! Since vertical distributions are very poorly known and crude, we do not need
    ! any precise meteo-dependent projection, crude will do the job
    !
    allocate(dmssrc%levFractDispVert(fu_NbrOfLevels(vert_disp)), &
           & dmssrc%fzDisp(fu_NbrOfLevels(vert_disp)), stat=stat)
    if(fu_fails(stat == 0, 'Allocate failed', subname)) return

    ! We use the same vertical as for sea salt.
    ! 
    call set_vertical(fu_set_layer_between_two(layer_btw_2_height, 0.0, 25.0), vertTmp)
    if(error)return
    call add_level(vertTmp, fu_set_layer_between_two(layer_btw_2_height, 25.0, 50.0))
    if(error)return

    if(error)return
    fractions = (/0.6, 0.4/)
    call report(vertTmp, .true.)
    dmssrc%levFractDispVert = 0.0
    dmssrc%fzDisp = 0.0  ! seems that reproject_verticals doesn't initialize zeros.
    call reproject_verticals(vertTmp, fractions, &                    ! vertical from, fractions from
                           & vert_proj, dmssrc%levFractDispVert, &   ! vertical to, fractions to
                           & dmssrc%fzDisp, dmssrc%nLevsDispVert, & ! mass centres, number of non-zero levels
                           & ifMassCentreInRelUnit=.true.)
    call set_missing(vertTmp, ifNew=.false.)
    call msg('Level fractions:', dmssrc%levFractDispVert)
    call msg('number of levels:', dmssrc%nLevsDispVert)
    call report(vert_proj, .true.)
    if (fu_fails(all(dmssrc%levFractDispVert <= 1.0), 'Bad level fractions', subname)) return
    if (fu_fails(all(abs(dmssrc%fzDisp) < 0.5), 'Bad fzDisp', subname)) return
  end subroutine project_dmssrc_second_grd


  !**************************************************************************

  subroutine link_dmssrc_to_species(species_emis, dmssrc)
    implicit none

    ! Imported parameters
    type(silam_dms_source), intent(inout) :: dmssrc
    type(silam_species), dimension(:), intent(in) :: species_emis

    character(len=*), parameter :: subname = 'link_dmssrc_to_species'

    if (fu_fails(associated(dmssrc%species_list), 'species_list not associated', subname)) return
    
    dmssrc%ind_emis = fu_index(dmssrc%species_list(1), species_emis)
    if (fu_fails(dmssrc%ind_emis > 0, 'Failed to link species', subname)) return

  end subroutine link_dmssrc_to_species


  !**************************************************************************

  subroutine compute_emission_for_dms(dmssrc, &
                                    & met_buf, disp_buf, & 
                                    & now, &      ! current time
                                    & timestep, & ! model time step
                                    & pHorizInterpMet2DispStruct, ifHorizInterp, &
                                    & ifSpeciesMoment, &
                                    & emisMap, mapCoordX, mapCoordY, mapCoordZ, & ! Output
                                    & fMassInjected)                              ! output
    implicit none

    ! Imported parameters
    type(silam_dms_source), target, intent(in) :: dmssrc
    type(Tfield_buffer), pointer ::  met_buf, disp_buf  ! meteo and internal field buffers
    type(silja_time), intent(in) :: now           ! current time
    type(silja_interval), intent(in) :: timestep  ! model time step
    type(THorizInterpStruct), intent(in) :: pHorizInterpMet2DispStruct
    logical, intent(in) :: ifHorizInterp, ifSpeciesMoment
    type(Tmass_map), intent(inout) :: emisMap, mapCoordX, mapCoordY, mapCoordZ
    real(r8k), dimension(:),  intent(inout) :: fMassInjected

    character(len=*), parameter :: subname = 'compute_emission_dms'
    integer :: ind_sst, ind_u10m, ind_v10m
    real :: timestep_sec, windspeed, weight_past, emis, fract_water, cnc_dms_water, lev_fract
    integer :: month, ilev, i1d_disp, ix, iy
    real, dimension(:), pointer :: ptrXSzDisp, ptrYSzDisp
    real :: cell_size, sst, u10m, v10m

    if (fu_fails(defined(dmssrc), 'source not defined', subname)) return
    if (fu_fails(ifSpeciesMoment, 'Not implemented for ifSpeciesMoment=false', subname)) return

    ind_sst = fu_index(met_buf, water_surface_temp_flag)
    if (fu_fails(ind_sst > 0, 'Cannot find SST', subname)) return
    ind_u10m = fu_index(met_buf, u_10m_flag)
    if (fu_fails(ind_sst > 0, 'Cannot find u10m', subname)) return
    ind_v10m = fu_index(met_buf, v_10m_flag)
    if (fu_fails(ind_sst > 0, 'Cannot find v10m', subname)) return
    
    timestep_sec = abs(fu_sec(timestep))
    month = fu_mon(now)
    weight_past = met_buf%weight_past
    if (fu_fails(dmssrc%grid == emisMap%gridTemplate, 'Grids don''t match', subname)) return

    ptrXSzDisp => fu_grid_data(dispersion_cell_x_size_fld)  ! get grid cell size in metres
    ptrYSzDisp => fu_grid_data(dispersion_cell_y_size_fld)
    
    do iy = 1, emisMap%ny
      do ix = 1, emisMap%nx
        u10m = fu_get_value(met_buf%p2d(ind_u10m), nx_meteo, ix, iy, weight_past, &
                          & pHorizInterpMet2DispStruct, ifHorizInterp, ifForceWeightPast_=.false.)
        v10m = fu_get_value(met_buf%p2d(ind_v10m), nx_meteo, ix, iy, weight_past, &
                          & pHorizInterpMet2DispStruct, ifHorizInterp, ifForceWeightPast_=.false.)
        windspeed = sqrt(u10m**2 + v10m**2)
        sst = fu_get_value(met_buf%p2d(ind_sst), nx_meteo, ix, iy, weight_past, &
                         & pHorizInterpMet2DispStruct, ifHorizInterp, ifForceWeightPast_=.false.)
        if (sst < dmssrc%sst_min) cycle

        i1d_disp = (iy-1) * emisMap%nx + ix
        cell_size = ptrXSzDisp(i1d_disp) * ptrYSzDisp(i1d_disp)

        cnc_dms_water = dmssrc%map_dms(month, ix, iy)
        if (fu_fails(cnc_dms_water >= 0.0, 'Bad cnc_dms', subname)) return
        fract_water = dmssrc%water_mask(ix,iy)
        if (fu_fails(0 <= fract_water, 'Bad fract_water', subname)) return

        emis = get_flux_dms(windspeed,sst,cnc_dms_water) * cell_size * timestep_sec * fract_water * dmssrc%yield
        
        fMassInjected(dmssrc%ind_emis) = fMassInjected(dmssrc%ind_emis) + emis
        emisMap%ifColumnValid(dmssrc%id_nbr,ix,iy) = .true.

        do ilev = 1, dmssrc%nLevsDispVert
          lev_fract = dmssrc%levFractDispVert(ilev)
          if (lev_fract < 1e-5) cycle
          emisMap%arm(dmssrc%ind_emis, dmssrc%id_nbr, ilev, ix, iy) = &
               & emisMap%arm(dmssrc%ind_emis, dmssrc%id_nbr, ilev, ix, iy) &
               & + emis * lev_fract
          mapCoordZ%arm(dmssrc%ind_emis, dmssrc%id_nbr, ilev, ix, iy) = &
               & mapCoordZ%arm(dmssrc%ind_emis, dmssrc%id_nbr, ilev, ix, iy) &
               & + emis * lev_fract * dmssrc%fzDisp(ilev)
          emisMap%ifGridValid(iLev,dmssrc%id_nbr) = .true.

        end do
        
      end do
    end do

  end subroutine compute_emission_for_dms

  real function get_flux_dms(windspeed, sst, cnc_dms_water) result(flux)
    implicit none
    real, intent(in) :: &
         & windspeed, & ! m/s at 10 m
         & sst ! K
    real, intent(in) :: cnc_dms_water ! mol/l
    
    real, parameter :: schmidt_ref = 660.0, & ! CO2 in 293 K
         & per_sec = 1.0 / 3600.0, to_m = 0.01, per_l_to_per_m3 = 1000.0, p23 = 2.0 / 3.0, p12 = 0.5
    real :: schmidt_dms, kw, sc_ratio, temp_c
    
    temp_c = sst - 273.15
    
    ! Saltzman et al., 1993
    schmidt_dms = 2674.0 - 147.12 * temp_c + 3.726 * temp_c**2 + 0.038*temp_c**3
    sc_ratio = schmidt_ref / schmidt_dms

    ! Liss & Merlivat 1986:
    if (windspeed <= 3.6) then
      kw = 0.17 * windspeed * sc_ratio ** p23
    else if (windspeed <= 13.0) then
      kw = 0.612 * sc_ratio**p23 + (2.85*windspeed - 10.26) * sc_ratio**p12
    else
      kw = 0.612 * sc_ratio**p23 + (5.90*windspeed - 49.90) * sc_ratio**p12
    end if
    kw = per_sec * to_m * kw
    
    flux = kw * cnc_dms_water * per_l_to_per_m3
  end function get_flux_dms

  !********************************************************************************************

  logical function fu_dms_emis_owned_quantity(dmssrc, quantity) result(owned)
    !
    ! Checsk whether the specific quantity is handled by the emission sources exclusively
    ! Such quantities may be related to emission species somewhere else
    !
    implicit none

    ! Imported parameters
    type(silam_dms_source), intent(in) :: dmssrc
    integer, intent(in) :: quantity
    owned = .false.
  end function fu_dms_emis_owned_quantity


  !*****************************************************************

  integer function fu_source_id_nbr_dmssrc(dmssrc) result(id)
    !
    ! Returns the source number. Reason: all sources are enumerated
    ! sequencially, so that the source can, in fact, be refered by its
    ! number without other information. But so far this number is 
    ! copied to the particles in the pollution cloud. 
    !
    implicit none

    ! Imported parameters with intent IN
    type(silam_dms_source), intent(in) :: dmssrc

    if (fu_fails(defined(dmssrc), 'Undefined source', 'fu_source_id_nbr_of_dmssrc')) return
    id = dmssrc%id_nbr

  end function fu_source_id_nbr_dmssrc



  !*************************************************************************

  integer function fu_source_nbr_dmssrc(dmssrc) result(nbr)
    !
    ! Returns the source number. Reason: all sources are enumerated
    ! sequentially, so that the source can, in fact, be refered by its
    ! number without other information. But so far this number is 
    ! copied to the particles in the pollution cloud. 
    !
    implicit none
    type(silam_dms_source), intent(in) :: dmssrc

    if (defined(dmssrc)) then
      nbr = dmssrc%src_nbr
    else
      nbr = int_missing
      call set_error('Undefined source', 'fu_source_nbr_of_dmssrc')
    end if

  end function fu_source_nbr_dmssrc


  !*************************************************************************

  subroutine typical_species_cnc_dms(dmssrc, species, nSpecies, arConc)
    !
    ! Guesses a typical level of concentration.
    !
    implicit none

    ! Imported parameters
    type(silam_dms_source), intent(in) :: dmssrc
    type(silam_species), dimension(:), pointer :: species
    integer, intent(out) :: nSpecies
    real, dimension(:), intent(inout) :: arConc
    
    character(len=*), parameter :: subname = 'typical_species_cnc_dms'
    
    if (fu_fails(defined(dmssrc), 'undefined source', subname)) return
    
    arConc(1) = 1e-8 ! mol/m3, ~0.25 ppb
    nSpecies = 1
    species => dmssrc%species_list

  end subroutine typical_species_cnc_dms


  !*************************************************************************

  function fu_name_dmssrc(dmssrc)result(chNm)
    implicit none
    type(silam_dms_source), intent(in) :: dmssrc
    character(len=clen) :: chNm
    chNm = dmssrc%src_name
  end function fu_name_dmssrc


  !*************************************************************************

  subroutine report_dmssrc(dmssrc)
    implicit none
    type(silam_dms_source), intent(in) :: dmssrc

    call msg('------------------ DMS source report -----------------')
    if (.not. defined(dmssrc)) then
      call msg('Source is undefined.')
      return
    end if
    
    call msg('Source and sector name:' // ' ' // trim(dmssrc%src_name) // ' ' // trim(dmssrc%src_sector_name))
    call msg('Emitted species:')
    call report(dmssrc%species_list(1))
    call msg('Yield (mol/mol):', dmssrc%yield)
    call msg('SST threshold:', dmssrc%sst_min)
    call msg('')
    
  end subroutine report_dmssrc

  !************************************************************************************

  logical function defined_dmssrc(dmssrc)
    implicit none
    type(silam_dms_source), intent(in) :: dmssrc
    defined_dmssrc = dmssrc%defined
  end function defined_dmssrc
  
  subroutine test_dmssrc_no_init()
    implicit none
    
    type(silam_dms_source) :: dmssrc
    type(silja_grid) :: grid
    
    grid = fu_set_lonlat_grid('grid', -10.0, 12.0, .true., 200, 300, pole_geographical, 0.2, 0.2)
    if (error) return
    
    dmssrc%map_dms_filename = 'dms_lana_2011.nc'
    dmssrc%map_dms_fileformat = fu_input_file_format('NETCDF:dms_clim')
    dmssrc%water_mask_filename = 'sslt_emission_global_50km.fld.nc'
    dmssrc%water_mask_fileformat = fu_input_file_format('NETCDF:dms_clim')
    dmssrc%emis_subst_name = 'DMS'
    call msg('Testing DMS source')

    call init_emission_dms(dmssrc, grid)

    call test_dmssrc(dmssrc)

  end subroutine test_dmssrc_no_init
  
  subroutine test_dmssrc(dmssrc)
    implicit none
    type(silam_dms_source), intent(inout) :: dmssrc
    
    integer :: month
    type(silam_vertical) :: vert

    call set_vertical(fu_set_layer_between_two(layer_btw_2_height, 0.0, 12.5), vert)
    call add_level(vert, fu_set_layer_between_two(layer_btw_2_height, 12.5, 25.0))
    call add_level(vert, fu_set_layer_between_two(layer_btw_2_height, 25.0, 50.0))
    
    call project_dmssrc_second_grd(dmssrc, dmssrc%grid, vert, vert, 5)
    if (error) return
    call msg('Level fractions:', dmssrc%levFractDispVert)
    call msg('fzdisp:', dmssrc%fzDisp)

    call report(dmssrc)

    do month = 1, 12
      call msg('month:', month)
      call open_binary(99, file=trim(fu_str(month)) // '.dat', access='stream', recl=0)
      write(99) dmssrc%map_dms(month, :, :)
      call msg('maxval:', maxval(dmssrc%map_dms(month,:,:)))
      call msg('minval:', minval(dmssrc%map_dms(month,:,:)))
    end do
    
    call msg('u sst dms:', (/5.0, 285.0, 5.0e-9/))
    call msg('flux:', get_flux_dms(5.0, 285.0, 5e-9) * 1e6 * 24*3600)
    call msg('flux mgS/year/m2:', get_flux_dms(5.0, 285.0, 5e-9) * 1000 * 32 * 3600*24*365)
    call msg('kw cm/h', get_flux_dms(15.0, 285.0, 1.0) * 3600 * 100 / 1000)
    
  end subroutine test_dmssrc

  subroutine test_init_dmssrc(filename)
    implicit none
    character(len=*), intent(in) :: filename
    
    type(Tsilam_namelist), pointer :: nlptr
    type(silam_dms_source) :: dmssrc
    type(silja_grid) :: grid
    integer :: file_unit, iostat
    character(len=500) :: line
    logical :: eof

    file_unit = fu_next_free_unit()
    open(file_unit, file=filename, action='read', status='old', iostat=iostat)
    if (fu_fails(iostat == 0, 'bad iostat', '')) return
    call next_line_from_input_file(file_unit, line, eof)
    print *, line
    nlptr => fu_read_namelist(file_unit, chStopLine='END_DMS_SOURCE_V5', ifenv=.true.)
    call report(nlptr)
    close(file_unit)
    
    call fill_dms_src_from_namelist(nlptr, dmssrc)
    if (error) return
    grid = fu_set_lonlat_grid('grid', -10.0, 12.0, .true., 200, 300, pole_geographical, 0.2, 0.2)
    call init_emission_dms(dmssrc, grid)
    
    call test_dmssrc(dmssrc)
    
  end subroutine test_init_dmssrc

END MODULE source_terms_dms


