MODULE source_terms_general

  ! This module contains the general description source-term of SILAM
  ! Source-term describes the spatial- and time- distribution of the 
  ! release to atmosphere, and the amount of chemical/radioactive 
  ! materials released.
  !
  ! In SILAM the source term consists of an unlimited number of
  ! release sources, each describing one release
  !
  ! In this module we define:
  !
  ! Words SOURCE and SOURCE TERM REFER TO FULL DESCRIPTION OF THE
  ! RELEASE SITUATION AND CIRCUMSTANCES (location, amounts of
  ! materials released, heat contents, spatial/vertical/time
  ! distribution etc...)
  !
  ! WORD "RELEASE" REFERS ONLY TO THE MATERIALS RELEASED TO THE
  ! ATMOSPHERE (material and their amounts)
  ! 
  ! NOTE. Co-ordinates of the point source (position) are always stored
  !       in the geographical grid as in ini file. System grid is not yet
  !       defined. BUT the cloud from this source will already be in
  !       system grid. This very point is a border between general
  !       geographical grids and the system grid, usually derived by
  !       the meteorological ones.
  !
  ! Currently the module contains descriptions for POINT SOURCE, BOMB SOURCE
  ! and AREA_SOURCE
  ! There can be two point source description formats - simple for the SILAM v.2
  ! and more sophisticated for the SILAM version 3. They called point_source_2 and
  ! point_source_3 correspondingly.
  !
  ! All units: SI
  ! 
  ! All degrees: real numbers in degrees (degrees and parts per
  ! hundred/parts per thousand , NOT in degrees and minutes)
  ! 
  ! Language: ANSI standard Fortran 90
  !
  ! Code owner: Mikhail Sofiev, FMI
  ! 
  ! Modules used:

  !$ use OMP_LIB
  
  use source_terms_area
  use source_terms_bio_voc
  use source_terms_bomb
  use source_term_fires
  use source_terms_point
  use source_terms_pollen
  use source_terms_sea_salt
  use source_terms_volcano
  use source_terms_wind_blown_dust
  use source_terms_dms

  IMPLICIT NONE
  
  private

  ! The public functions and subroutines available in this module:
  
  ! Functions universal for all types of the sources
  !
  public set_source_terms
  public add_source_term_input_needs ! From meteorological modules
  PUBLIC defined
  PUBLIC report
  public store_as_namelist      ! technical dump
  public verify_sources
  public fu_if_dump_emission_flux
  public dump_emission_mass_map ! into time-resolving area source
  public close_emission_dump
  PUBLIC fu_name
  PUBLIC fu_start_time
  PUBLIC fu_end_time
  PUBLIC fu_duration
  public fu_nbr_of_disp_grd_cells
  public link_source_to_species
  public fu_grid
  public fu_source_id_nbr
  public fu_NbrOf_sources_total
  public fu_NbrOf_source_ids
  public fu_get_id_str_from_id_nbr ! Finds the src Id string having the Id number
  public create_src_contain_grd_gen_src
!  public create_src_containing_vert
  public force_source_into_grid
  public source_2_map
  public init_emission_internal_fields
  public link_emission_2_dispersion
  public inject_emission_eulerian
  public inject_emission_lagrangian
  public get_inventory_g_src
  public amounts_from_src_species_unit
  public typical_cnc_from_src_species_unit
  public fu_emission_owned_quantity
  public fu_check_outTemplate ! Checks consistency: outTemplate, ifSplitSources, FileArrangement
  public make_source_id_mapping
  
  ! Parameter assimilation routines
  public assimilation_request_emis
  public observe_params_emis
  public inject_params_emis


  ! The private functions and subroutines for general source

  private count_srcs_in_all_data_files
  private count_srcs_one_data_file
  private allocate_source_space
  private read_sources_from_file
  private set_next_gen_src_from_nl
  PRIVATE fu_source_defined
  PRIVATE report_source
  PRIVATE fu_source_name
  PRIVATE fu_source_start_time
  PRIVATE fu_source_end_time
  PRIVATE fu_earliest_start_time
  PRIVATE fu_latest_end_time
  PRIVATE fu_source_duration
  
  private fu_n_disp_grd_cells_src_inv
  private fu_source_id_nbr_of_source
  private source_2_map_general
  private fu_grid_of_source
  private fu_id_nbr_from_srcNMs
  private fu_id_nbr_from_srcNbr
  private fu_compute_src_id_nbr
  private link_general_src_to_species
  private force_general_src_into_grid
  private fu_src_id
  private store_as_namelist_general_src

  ! Generic names and operator-interfaces of some functions:

  INTERFACE defined
    MODULE PROCEDURE fu_source_defined
  END INTERFACE

  INTERFACE report
    MODULE PROCEDURE report_source
  END INTERFACE

  interface store_as_namelist
    module procedure store_as_namelist_general_src
  end interface

  INTERFACE fu_name
    MODULE PROCEDURE fu_source_name
  END INTERFACE

  interface fu_grid
    module procedure fu_grid_of_source
  end interface

  INTERFACE fu_start_time
    MODULE PROCEDURE fu_source_start_time
    MODULE PROCEDURE fu_earliest_start_time
  END INTERFACE

  INTERFACE fu_end_time
    MODULE PROCEDURE fu_source_end_time
    MODULE PROCEDURE fu_latest_end_time
  END INTERFACE

  INTERFACE fu_duration
    MODULE PROCEDURE fu_source_duration
  END INTERFACE

  interface fu_nbr_of_disp_grd_cells
    module procedure fu_n_disp_grd_cells_src_inv
  end interface

  interface fu_source_id_nbr
    module procedure fu_source_id_nbr_of_source
  end interface

  interface link_source_to_species
    module procedure link_general_src_to_species
  end interface

  interface source_2_map
    module procedure source_2_map_general
  end interface

  interface force_source_into_grid
    module procedure force_general_src_into_grid
  end interface


  ! Public types with private components defined in this module:
  ! The main SILAM source term. Contains pointers to all possible
  ! source types. There may be plenty of the sources of each type
  !
  ! ATTENTION!! Each source (except for bomb one) is uniquely identified
  !             by means of the following parameters: src_nm + sector_nm
  !             or src_nbr. src_nbr is absolutely unique, while src_nm and
  !             sector_nm each may show up several times (but their
  !             combination is unique). This trick allows the following:
  !          1. Complete separation of all source-sector pairs
  !          2. Split of sources, with sectors summed-up
  !          3. Split of sectors, with sources summed-up
  !             This informaion is stored in iSrcIdType field, one for all sources.
  !             
  type silam_source_info
    private
    character(len=clen) :: chSrcNm, chSectorNm, chId
    character(len=clen), dimension(max_descriptors) :: chDescrNames
    integer :: iSrcType, iSrcNbr, iIdNbr, iHorizInterp, iVertInterp, nChemDescr, nCells, iDynamicEnvironment
    logical :: ifNeedsTZindexMap
  end type silam_source_info
  type(silam_source_info) :: silam_source_info_missing= &
          & silam_source_info('','','','', &
          & int_missing,  int_missing, int_missing, int_missing, &
          & int_missing, int_missing, int_missing, int_missing, &
          & .False.)

  integer, parameter :: maxMeteoDepndencies = 5
  TYPE silam_source
    INTEGER :: n_area=0, n_bio_voc=0, n_bomb=0, n_fire=0, n_point=0, &
             & n_pollen=0, n_sea_salt=0, n_volc=0, n_wb_dust=0, n_dms=0
    INTEGER :: iSrcIdType                                ! type of source id
    type(silam_source_info), dimension(:), allocatable :: src_info_lst
    TYPE(silja_logical) :: defined
    TYPE(a_src_ptr), dimension(:), POINTER :: a_ptr => null()
    TYPE(b_src_ptr), dimension(:), POINTER :: b_ptr  => null()
    TYPE(bvoc_src_ptr), dimension(:), POINTER :: bvoc_ptr  => null()
    TYPE(fire_src_ptr), dimension(:), POINTER :: fire_ptr  => null()
    TYPE(p_src_ptr), dimension(:), POINTER :: p_ptr  => null()
    type(pollen_src_ptr), dimension(:), pointer :: pollen_ptr  => null()
    type(sslt_src_ptr), dimension(:), pointer :: sslt_ptr  => null()
    type(wb_dust_src_ptr), dimension(:), pointer :: wbdust_ptr  => null()
    type(dms_src_ptr), dimension(:), pointer :: dms_ptr => null()
    type(volc_src_ptr), dimension(:), pointer :: volc_ptr => null()
    logical :: ifDumpFlux, ifDumpMoment
    type(silja_time) :: FirstDumpTime, LastDumpTime, prevDumpTime, nextDumpTime
    type(silja_interval) :: DumpTimestep
    integer :: DumpFilesArrangement, indexDumpFile
    integer, dimension(:,:), allocatable :: arSourceIdMapping
    character(len=10), dimension(:), allocatable :: chSplitNames
    logical :: ifNeedsTZindexMap, ifUseSourceIdMapping
    integer, dimension(maxMeteoDepndencies) :: arMeteoDepndencies
  END TYPE silam_source
  
  public silam_source

  type Tsource_summary
    private
    character(len=30) :: label
    integer :: source_type
  end type Tsource_summary
    
  type(Tsource_summary), dimension(10), parameter :: known_src_types = (/ &
              & Tsource_summary('AREA_SOURCE_', area_source), &
              & Tsource_summary('BOMB_SOURCE_', bomb_source), &
              & Tsource_summary('DMS_SOURCE_', dms_source), &
              & Tsource_summary('FIRE_SOURCE_', fire_source), &
              & Tsource_summary('POINT_SOURCE_', point_source), &
              & Tsource_summary('POLLEN_SOURCE_', pollen_source), &
              & Tsource_summary('BIOGENIC_VOC_SOURCE_', bio_voc_source), &
              & Tsource_summary('SEA_SALT_SOURCE_', sea_salt_source), &
              & Tsource_summary('VOLCANO_SOURCE_', volcano_source), &
              & Tsource_summary('WIND_BLOWN_DUST_SOURCE_', wind_blown_dust_source)/)
  
  !
  ! A source info - name, sector, number, type, etc. So, some information that is not 
  ! related to the source type and should be filled-in when the source file is consumed
  ! It duplicates the internal source info but is not source-type-specific, which would
  ! much simplify the gathering of across-type information like source id
  ! Alternative to this structue is to kill the type-specific sources, which will make
  ! structure of the module much worse and would significantly complicate additions of new
  ! features to the specific source types or new types.
  ! Second piece of information: 
  ! Indices of horizontal and vertical interpolation structures. Each source may have own
  ! grid but it very likely that for some sources the grids will be the same, as well as
  ! verticals. Therefore, interpolation structures will be made separately and sources will
  ! just know their "own" structures pointed by the indices. Place for computing these
  ! structures is the same as that of set_lagrangian_particles - all major run parameters are 
  ! defined by that time.
  !

  ! Type of the source id - name, sector, or both
  ! 
  integer, public, parameter :: iNoId = 1120  ! Mix all sources
  integer, public, parameter :: iSrcNmId = 1121 ! Separate sources according to names
  integer, public, parameter :: iSrcSectorId = 1122 ! Separate sources according to sectors
  integer, public, parameter :: iSrcNmSectorId = 1123 ! Separate sources according to names and sectors
  integer, public, parameter :: iSrcTimeZoneId = 1124 ! Separate sources according to names and sectors

  real(r8k), private, allocatable, target :: dInjectedMassThread(:,:) !nThreads, nSpeciesEmit


  CONTAINS


  !****************************************************************

  subroutine set_source_terms(em_source, nlEmission, iSrcIdType, ifEmissionGiven, expected_species, &
                            & iSimulationType, ifOldNamelistFormat)
    !
    ! Responsible for the source term initial reading and setting. Essentially, calls the
    ! next-in-line subroutine if there is anything to read at all.
    !
    implicit none

    ! Imported parameters 
    type(silam_source), intent(out) :: em_source
    type(Tsilam_namelist), pointer :: nlEmission
    integer, intent(in) :: iSrcIdType, iSimulationType
    logical, intent(out) :: ifEmissionGiven
    type(silam_species), dimension(:), intent(in), allocatable :: expected_species
    logical, intent(in) :: ifOldNamelistFormat

    ! Local variables
    type(Tsilam_nl_item_ptr), dimension(:), pointer :: pSourceTerms, pDataItems
    integer :: iTmp, jTmp, iSrc, nSourceTerms, nDataFiles, iFile, status
    character(len=fnlen) :: chTmp
    integer, dimension(:), pointer :: nValLinesInventory
    type(Tsilam_namelist), pointer :: nlInventoryDataFilesGlobal
    character(len=*), parameter :: sub_name = 'set_source_terms'

    !
    ! Get all the sources listed in the control file
    !
    em_source%n_area=0
    em_source%n_bio_voc=0
    em_source%n_bomb=0
    em_source%n_fire=0
    em_source%n_point=0
    em_source%n_pollen=0
    em_source%n_sea_salt=0
    em_source%n_volc=0
    em_source%n_wb_dust=0
    em_source%n_dms=0
    nullify(em_source%a_ptr)
    nullify(em_source%bvoc_ptr)
    nullify(em_source%b_ptr)
    nullify(em_source%fire_ptr)
    nullify(em_source%p_ptr)
    nullify(em_source%pollen_ptr)
    nullify(em_source%sslt_ptr)
    nullify(em_source%wbdust_ptr)
    nullify(em_source%volc_ptr)
    nullify(em_source%dms_ptr)
!    nullify(em_source%arSourceIdMapping)  ! not pointer

    em_source%ifNeedsTZindexMap = .False.
    em_source%ifUseSourceIdMapping = .false.

    nullify(pSourceTerms)
    call get_items(nlEmission, 'emission_source', pSourceTerms, nSourceTerms)
    if(error)return

    if(nSourceTerms < 1)then
      call msg_warning('No emission_source lines in the control file',sub_name)
      call report(nlEmission)
      call set_error('No emission_source lines in the control file',sub_name)
      return
    endif

    !-------------------------------------------------------------------------
    !
    ! VOID_SOURCE
    !
    ! Do we have anything to read?
    !
    if(fu_str_u_case(fu_content(pSourceTerms(1))) == 'VOID_SOURCE')then
      if(nSourceTerms == 1)then
        call msg('')
        call msg('Void source term, no dispersion simulations are planned')
        call msg('')
        ifEmissionGiven = .false.
      else
        call msg_warning('VOID_SOURCE is given but there are others',sub_name)
        call report(nlEmission)
        call set_error('VOID_SOURCE is given but there are others',sub_name)
      endif
      return
    else
      ifEmissionGiven = .true.
    endif

    !
    ! A bit of preparations to start reading
    !
    nValLinesInventory => fu_work_int_array()


    !-------------------------------------------------------------------------
    !
    ! Scan the requested source terms and initialise them one by one.
    ! Options are (not all are available so far):
    !
    ! VOID_SOURCE - must be the first line, the others are ignored. Means meteo-processing run
    ! INVENTORY   - data from emission database. Now - point and area sources, bomb so far here as well
    ! SEA_SALT    - meteo-driven emission of sea salt particles
    ! BIOGENIC_VOC- meteo-driven emission of biogenic volatile organic compounds
    ! WILD_LAND_FIRES- source term for wild-land fires
    ! DESERT_DUST - meteo-driven emission of dust from deserts
    ! POLLEN      - pollen and allergen
    !
    call start_count('count_sources')
    call count_srcs_in_all_data_files(em_source, pSourceTerms, nSourceTerms, &
                                    & nlInventoryDataFilesGlobal, & ! list of inventory data files
                                    & nValLinesInventory, &         ! dimensions of inventory val lines
                                    & iSrcIdType, iSimulationType, ifOldNamelistFormat)
    call stop_count('count_sources')
    if(error)return
    !
    ! Allocate the source list and read the data files
    !
    call allocate_source_space(em_source)
    if(error)return

    !--------------------------------------------------------------------------------
    !
    ! Now we can consume the data files and store the structures in the source.
    !
    ! INVENTORY files are all stored in the nlDataFilesGlobal
    ! SEA_SALT, POLLEN, and BIOGENIC_VOC, WIND_BLOWN_DUST files are also there, just different 
    ! procedure of initialization
    !
    nullify(pDataItems)
    call get_items(nlInventoryDataFilesGlobal, 'data', pDataItems, nDataFiles)
    if(error)return
    if(nDataFiles < 1)then
      call report(nlInventoryDataFilesGlobal)
      call set_error('No data files to be read in the above namelist',sub_name)
      return
    endif
    
    call start_count('read_sources')
    do iFile = 1, nDataFiles

      !
      ! Some species can be prescribed by the chemical transformations (actually, by aerosol dynamics)
      ! Then be careful: some sources emit in continuous spectrum and should be projected to this
      ! aerosol. Others will be later re-projected to it.
      !
      call read_sources_from_file(fu_content(pDataItems(iFile)), &
                                & em_source, iSrcIdType, nValLinesInventory, iFile, &
                                & expected_species)
      if(error)return

    end do
    call stop_count('read_sources')

    call msg('Source reading is over')

    call destroy_namelist(nlInventoryDataFilesGlobal)
    call free_work_array(nValLinesInventory)

!    call msg('store as namelist')
!    call store_as_namelist(em_source, 'src.dump2', .true.)


!    !-------------------------------------------------------------------------
!    !
!    ! Other source types:
!    ! BIOGENIC_VOC- meteo-driven emission of biogenic volatile organic compounds
!    ! WILD_LAND_FIRES- source term for wild-land fires
!    ! DESERT_DUST - meteo-driven emission of dust from deserts
!    ! POLLEN      - pollen and allergen
!    !

    !-------------------------------------------------------------------------
    !
    ! Finally, set the source dump if it is needed
    !
    em_source%ifDumpFlux = .false.
    if(len_trim(fu_content(nlEmission,'source_dump_time_step')) > 0)then
      em_source%ifDumpFlux = .true.
      em_source%ifDumpMoment = index(fu_str_u_case(fu_content(nlEmission,'if_dump_emission_moments')), 'YES') > 0
      if(len_trim(fu_content(nlEmission,'source_dump_start_time')) > 0)then
        em_source%FirstDumpTime = &
                          & fu_io_string_to_time(fu_content(nlEmission,'source_dump_start_time'))
        if(len_trim(fu_content(nlEmission,'source_dump_end_time')) > 0)then
          em_source%LastDumpTime = &
                          & fu_io_string_to_time(fu_content(nlEmission,'source_dump_end_time'))
        else
          if(len_trim(fu_content(nlEmission,'source_dump_period'))> 0)then
            em_source%LastDumpTime = em_source%FirstDumpTime + &
                                     & fu_set_named_interval(fu_content(nlEmission,'source_dump_period'))
          else
            em_source%LastDumpTime = time_missing
          endif
        endif
      else
        em_source%FirstDumpTime = time_missing
        em_source%LastDumpTime = time_missing
      endif
      em_source%DumpTimestep = fu_set_named_interval(fu_content(nlEmission,'source_dump_time_step'))
      
      chTmp = fu_str_u_case(fu_content(nlEmission,'source_dump_time_split'))
      if(chTmp== 'ALL_IN_ONE')then
        em_source%DumpFilesArrangement = all_in_one
      elseif(chTmp == 'HOURLY_NEW_FILE')then
        em_source%DumpFilesArrangement = hourly_new_file
      elseif(chTmp == 'DAILY_NEW_FILE')then
        em_source%DumpFilesArrangement = daily_new_file
      elseif(chTmp == 'MONTHLY_NEW_FILE')then
        em_source%DumpFilesArrangement = monthly_new_file
      elseif(chTmp == 'YEARLY_NEW_FILE')then
        em_source%DumpFilesArrangement = yearly_new_file
      else
        call set_error('Unknown source_dump_time_split:' + chTmp, sub_name)
        return
      end if

      em_source%prevDumpTime = em_source%FirstDumpTime
      em_source%nextDumpTime = em_source%FirstDumpTime

    endif  ! if source is to be dumped to the time-resolving area source

    em_source%indexDumpFile = int_missing

    !
    ! When all sources are set, check the meteodependence for those anthropogenic ones that can have it.
    ! Area and point are the ones to check
    !
    iTmp = 1
    em_source%arMeteoDepndencies = int_missing
    do iSrc = 1, em_source%n_area
      jTmp  = fu_meteodep_model(em_source%a_ptr(iSrc)%a_src)
      if (any(em_source%arMeteoDepndencies == jTmp)) cycle !already there
      call msg("adding meteodependence for a_src", jTmp, maxMeteoDepndencies)
      em_source%arMeteoDepndencies(min(iTmp,maxMeteoDepndencies)) = jTmp
      iTmp = iTmp + 1 
    end do
    do iSrc = 1, em_source%n_point
       jTmp  =   fu_meteodep_model(em_source%p_ptr(iSrc)%p_src) 
       if (any(em_source%arMeteoDepndencies == jTmp)) cycle
       call msg("adding meteodependence for p_src", jTmp, maxMeteoDepndencies)
       em_source%arMeteoDepndencies(min(iTmp,maxMeteoDepndencies)) = jTmp
       iTmp = iTmp + 1
    end do
    if (iTmp > maxMeteoDepndencies) then
      call set_error("Too many meteo dependencies!!", sub_name)
      return
    endif
    
    
    !-------------------------------------------------------------------------
    !
    ! The sources have been consumed. Check and report the whole thing
    !
    call msg('')
    em_source%defined = silja_true
    call msg('All sources used for the run (short report):')
    CALL report(em_source, .false.) 
    call msg('')

  end subroutine set_source_terms


  !********************************************************************************************

  subroutine count_srcs_in_all_data_files(em_source, &                  ! the main source
                                        & pSourceTerms, nSourceTerms, & ! from control file
                                        & nlInventoryDataFilesGlobal, & ! list of inventory data files
                                        & nValLinesInventory, &         ! dimensions of inventory val lines
                                        & iSrcIdType, &           ! way the sources are split in the run output
                                        & iSimulationType, &
                                        & ifOldFileFormat)
    !
    ! Scans all the declared source terms and allocates proper space in the em_source
    !
    implicit none

    ! Imported parameters
    TYPE(silam_source), intent(inout) :: em_source  !should be zeroed before
    type(Tsilam_nl_item_ptr), dimension(:), pointer :: pSourceTerms
    integer, intent(in) :: nSourceTerms, iSrcIdType, iSimulationType
    integer, dimension(:), pointer :: nValLinesInventory
    type(Tsilam_namelist), pointer :: nlInventoryDataFilesGlobal
    logical, intent(in) :: ifOldFileFormat

    ! Local variables
    type(silam_source_info), dimension(:), allocatable :: srcInfoTmp
    character(len=fnlen), dimension(:), pointer :: fnames
    character(len=fnlen) :: line, linetmp
    character(len=30) :: chType, chDynamics
    type(Tsilam_namelist), pointer :: nlInventoryDataFiles
    logical :: eof, ifList, ifData
    integer :: file_unit, status, nDataFiles, iFile, iFileGlobal, iDataFileGlobal, &
             & iSrc, iTmp, nSrcRead, nSpecificSources
    type(Tsilam_nl_item_ptr), dimension(:), pointer :: pDataItems
    integer, dimension(:), pointer :: iDynamicsType

    !
    ! We allocate the temporary space for the source IDs - they will be filled in
    ! while counting sources. Note that they will incorporate the possibility
    ! to merge sources, thus saving time and memory. They also will take into account
    ! that merged sources will contain more than one descriptor
    !
    call enlarge_source_info(srcInfoTmp)
    if(error)return

    eof = .false.
    ifList = .false.
    ifData = .false.
    nSrcRead = 0
    iDataFileGlobal = 1
    iDynamicsType => fu_work_int_array()

    nlInventoryDataFilesGlobal => fu_make_missing_namelist()
    if(error)return

    allocate(fnames(nSourceTerms),stat=status)
    if(fu_fails(status == 0,'Failed to allocate temporary file names array','count_srcs_in_all_data_files'))return
    do iSrc = 1, nSourceTerms
      fnames(iSrc) = ''
    enddo

    !-------------------------------------------------------------
    !
    ! Scan the control-file list of sources adding them one-by-one to the proper places.
    ! Note that a single file may be a list of files and probably is a list of sources.
    !
    ! This dictates the procedure:
    ! Iteration 1. Go through the types of the source terms and store "simple" ones
    ! Iteration 2. Explore the files
    !
    nSpecificSources = 0
    do iSrc = 1, nSourceTerms
      line = adjustl(fu_content(pSourceTerms(iSrc)))
      linetmp = fu_str_u_case(line(1:index(line,' ')))
      !
      ! As of 2.11.2015, the words like INVENTORY etc are not needed and can be skipped.
      ! For the sake of backward compatibility, here they are cut out.
      !
      if(index(linetmp,' INVENTORY ') + index(linetmp,' POLLEN ') + index(linetmp,' BIOGENIC_VOC ') +&
       & index(linetmp,' SEA_SALT ') + index(linetmp,' WIND_BLOWN_DUST ') + &
       & index(linetmp,' WILD_LAND_FIRES ') + index(linetmp,' VOLCANO ') + index(linetmp,' DMS ') > 0)then
        call msg_warning('Deprecated INVENTORY/... source type specification. Skipped')
        line = trim(line((index(line,' ')+1):))       ! cut out the type of the source
        line = adjustl(line)
      endif

      nSpecificSources = nSpecificSources + 1

      if(ifOldFileFormat)then
        iDynamicsType(iSrc) = eulerian_flag
        fnames(nSpecificSources) = fu_process_filepath(line)
        if(index(line,'EULERIAN ') > 0 .or. index(line,'LAGRANGIAN ') > 0) &
              & call msg_warning('EULERIAN/LAGRANGIAN attributes are not available in setup_v4_7', &
                               & 'count_srcs_in_all_data_files')
      else
        iTmp = index(line,' ') 
        chDynamics = fu_str_u_case(line(1:iTmp-1))
        line = trim(line((iTmp+1):))             ! cut out the dynamics
        line = adjustl(line)
        if(chDynamics == 'EULERIAN')then
          iDynamicsType(iSrc) = eulerian_flag
        elseif(chDynamics == 'LAGRANGIAN')then
          iDynamicsType(iSrc) = lagrangian_flag
        else
          call msg_warning('Allowed dynamics types: EULERIAN, LAGRANGIAN, not:' + chDynamics, &
                         & 'count_srcs_in_all_data_files')
          call set_error('Unknown dynamics type in the line:' + line,'count_srcs_in_all_data_files')
          return
        endif
!        call msg("chDynamics:"+ chDynamics)
!        call msg("iSimulationType, iSrc", iSimulationType, iSrc)
!        call msg("iDynamicsType(iSrc)", iDynamicsType(iSrc))
        if(fu_fails(iSimulationType == hybrid_flag .or. iSimulationType == iDynamicsType(iSrc), &
                        & 'Dynamics type does not correspond to simulation type:' + line, &
                        & 'count_srcs_in_all_data_files'))return
        fnames(nSpecificSources) = fu_process_filepath(line)
          
      endif  ! old or new file format: EULERIAN / LAGRANGIAN / HYBRID options

    end do  ! nSourceTerms

    !
    ! Cycle over the ini files
    !
    nullify(pDataItems)

    do iFileGlobal = 1, size(fnames)
      if(len_trim(fnames(iFileGlobal)) == 0)exit  ! all files processed
      !
      ! Open file
      !
      file_unit = fu_next_free_unit()
      if(error)return
!call msg(fnames(iFileGlobal))
      OPEN(file=fnames(iFileGlobal), unit=file_unit, action='read', status='old', iostat=status)
      IF (status /= 0) THEN
        call msg('iostat', status)
        call msg('unit', file_unit)
        CALL set_error('cannot open input file:' + fnames(iFileGlobal),'read_all_sources_from_file')
        RETURN
      END IF

      !
      ! Find out the type of the file - a list of other data files or the data-containing file
      !
      do while(.not.eof)
        CALL next_line_from_input_file(file_unit, line, eof)

        IF(index(line,'AREA_SOURCE_') + &
         & index(line,'BOMB_SOURCE_') + &
         & index(line,'FIRE_SOURCE_') + &
         & index(line,'POINT_SOURCE_') + &
         & index(line,'POLLEN_SOURCE_') + &
         & index(line,'BIOGENIC_VOC_SOURCE_') + &
         & index(line,'VOLCANO_SOURCE_') + &
         & index(line,'WIND_BLOWN_DUST_') + &
         & index(line,'SEA_SALT_SOURCE_') + &
         & index(line, 'DMS_SOURCE_') > 0)then
          ifData = .true.
          ifList = .false.
          exit
        endif

        if(index(line,'data_file =') > 0)then
          ifData = .false.
          ifList = .true.
          exit
        endif
      end do
      if(.not.(ifData .or. ifList))then
        call set_error('File contains neither data nor list:' + fnames(iFileGlobal), &
                     & 'read_all_sources_from_file')
        return
      endif

      if(ifData)then
        !
        ! File contains the data or one-step references to the data files. 
        !
        call msg('Counting the simple data file "' + fnames(iFileGlobal)+'" nSrcRead so far:', nSrcRead)

        call count_srcs_one_data_file(file_unit, &
                                    & fu_dirname(fnames(iFileGlobal)), &
                                    & em_source, &
                                    & iSrcIdType, &
                                    & nValLinesInventory, iDataFileGlobal, &
                                    & iDynamicsType(iFileGlobal), &
                                    & srcInfoTmp, nSrcRead)
        if(error)return
        iDataFileGlobal = iDataFileGlobal + 1
        close(file_unit)
        !
        ! This file has been counted. Store its name for the next-step reading.
        ! Clumsy but handy.
        !
        call add_namelist_item(nlInventoryDataFilesGlobal,'data',fnames(iFileGlobal))
        if(error)return

      else
        !
        ! The file is essentially a list of data files. Get its length and content
        !
        call msg('Exploring the multi-data list file:' + fnames(iFileGlobal))
        rewind(file_unit)
        eof = .false.
        nlInventoryDataFiles => fu_read_namelist(file_unit, .true.)
        call get_items(nlInventoryDataFiles, 'data_file', pDataItems, nDataFiles)
        if(error)return
        close(file_unit)
        !
        ! Scan the files one-by-one in order to count their total number
        !
        do iFile = 1, nDataFiles
          linetmp = fu_extend_grads_hat(fu_content(pDataItems(iFile)), fnames(iFileGlobal))

          call msg('Counting the sources for:' + linetmp)
          OPEN(file_unit, &
             & file = linetmp, &
             & action='read', status='old', iostat=status)
          IF (status /= 0) THEN
            CALL set_error('cannot open data file:' + linetmp ,&
                         & 'read_all_sources_from_file')
            RETURN
          END IF
          call count_srcs_one_data_file(file_unit, &
                                      & fu_dirname(linetmp), &
                                      & em_source, &
                                      & iSrcIdType, &
                                      & nValLinesInventory, iDataFileGlobal, &
                                      & iDynamicsType(iFileGlobal), &
                                      & srcInfoTmp, nSrcRead)
          if(error)return
          close(file_unit)
          iDataFileGlobal = iDataFileGlobal + 1
          !
          ! This file has been counted. Store its name for the next-step reading.
          ! Clumsy but handy.
          !
          call add_namelist_item(nlInventoryDataFilesGlobal,'data',linetmp)
          if(error)return

        end do  ! cycle over the low-level data files

        call destroy_namelist(nlInventoryDataFiles)
        nullify(pDataItems)

      end if  ! if the specific inventory file is the data or a list

    end do  ! cycle over the control-file given inventories

    !
    ! Allocate the source info list in the main structure
    !
    if (fu_fails(nSrcRead == fu_num_sources_total(em_source), &
               & 'Counts don''t match', 'count_srcs_in_all_data_files')) return
    allocate(em_source%src_info_lst(nSrcRead), stat = status)
    if(fu_fails(status==0,'Failed to allocate source_info list','count_srcs_in_all_data_files'))return
    !
    ! Copy the source_info structure from the temporary place
    !
    do iSrc = 1, fu_num_sources_total(em_source)
      em_source%src_info_lst(iSrc)   = srcInfoTmp(iSrc)
      em_source%ifNeedsTZindexMap = em_source%ifNeedsTZindexMap .or. srcInfoTmp(iSrc)%ifNeedsTZindexMap
      em_source%src_info_lst(iSrc)%iHorizInterp= int_missing
      em_source%src_info_lst(iSrc)%iVertInterp = int_missing
      if(em_source%src_info_lst(iSrc)%iDynamicEnvironment == lagrangian_flag)then
        if(em_source%src_info_lst(iSrc)%iIdNbr > 99)then
          call set_error('source ID nbr > 99, too much for lagrangian environment', &
                       & 'count_srcs_in_all_data_files')
          return
        endif
      endif
    end do

    if (em_source%ifNeedsTZindexMap) then
            call msg("TZindexMap needed")
    else
            call msg("TZindexMap not needed")
    endif

    if(nSpecificSources < nSourceTerms)then
      call msg('nInitialisedSources < nSourceTerms:', nSpecificSources, nSourceTerms)
      call set_error('nInitialisedSources < nSourceTerms:', 'count_srcs_in_all_data_files')
    endif

    deallocate(srcInfoTmp)
    call free_work_array(iDynamicsType)

  end subroutine count_srcs_in_all_data_files


  !********************************************************************************************

  subroutine count_srcs_one_data_file(file_unit, &
                                    & dirname, &
                                    & em_source, &
                                    & iSrcIdType, & 
                                    & nValLines, iFile, &
                                    & iDynamicsTypeDefault, &
                                    & srcInfoTmp, nSrcRead)
    !
    ! Reads the file through and counts the number of sources
    ! Also counts the number of val lines in all area sources.
    !
    implicit none

    ! Imported parameters
    integer, intent(in) :: file_unit, iFile, iSrcIdType, iDynamicsTypeDefault
    character(len=*), intent(in) :: dirname ! of the file being read
    TYPE(silam_source), target, intent(inout) :: em_source
    integer, dimension(:), pointer :: nValLines
    type(silam_source_info), dimension(:), allocatable :: srcInfoTmp
    integer, intent(inout) :: nSrcRead   ! number of sources read so far

    ! Local variables
    integer :: iTmp, nValLinesTmp, stat, size_from_nc, size_from_nl, nctag, iSrc
    logical :: eof, ifArea, ifFound, ifStartLine
    character(len=fnlen) :: line, nc_entry, ncfilename
    character(len=clen) :: chTmp
    character(len=SubstNmLen) :: chSubstNm
    integer :: iModeNbr, iStarted
    character(len=*), parameter :: sub_name = 'count_srcs_one_data_file'

    !
    ! Scan the input file and check the number of sources of each type to read.
    !
    rewind(file_unit)
    line = ""
    eof = .false.
    nValLines(iFile) = 0
    ifArea = .false.
    iStarted = int_missing
    nValLinesTmp = 0
    chSubstNm = ''

    do while(.not.eof)
      !
      ! Get the type of the source
      !
      CALL next_line_from_input_file(file_unit, line, eof)
      if(error .or. eof)exit
      
      ifStartLine = .false.
      do iSrc = 1, size(known_src_types)
        if(index(line, trim(known_src_types(iSrc)%label)) == 1)then
          ifArea = known_src_types(iSrc)%source_type == area_source
          if(iStarted /= int_missing)then
            call set_error('Start new source without ending previous:'+line, &
                         & sub_name)
            return
          endif
          iStarted = iSrc
          srcInfoTmp(nSrcRead+1)%iSrcType = known_src_types(iSrc)%source_type
          if(ifArea) nValLinesTmp = 0
          ifStartLine = .true.
          exit
        endif
      end do  ! first, check if the new source starts
      if(ifStartLine)cycle
      !
      ! If not a start line, other information may be of use
      !
      if (index(line,'source_name') > 0)then
        if(iStarted /= int_missing)then
          iTmp = index(line,'=')
          srcInfoTmp(nSrcRead+1)%chSrcNm = adjustl(line(iTmp+1:))
        else
          call set_error('source_name is found prior to source start line', &
                       & sub_name)
          return
        endif

      elseif (index(line,'source_sector_name') > 0)then
        if(iStarted /= int_missing)then
          iTmp = index(line,'=')
          srcInfoTmp(nSrcRead+1)%chSectorNm = adjustl(line(iTmp+1:))
        else
          call set_error('source_sector_name is found prior to source start line', &
                       & sub_name)
          return
        endif

      elseif (index(line,'emitted_substance') > 0)then
        if(iStarted /= int_missing)then
          iTmp = index(line,'=')
          chSubstNm = adjustl(line(iTmp+1:))
!          srcInfoTmp(nSrcRead+1)%chSubstNm(1) = adjustl(line(iTmp+1:))
        else
          call set_error('source_sector_name is found prior to source start line', &
                       & sub_name)
          return
        endif

      elseif (index(line,'emitted_size_mode_nbr') > 0)then
        if(iStarted /= int_missing)then
          read(unit=line,fmt=*,iostat=iTmp) chTmp, chTmp, iModeNbr
!          read(unit=line,fmt=*) chTmp, chTmp, srcInfoTmp(nSrcRead+1)%iModeNbr(1)
          if(iTmp /= 0)then
            call set_error('Cannot get mode number from the line:'+line, &
                         & sub_name)
            return
          endif
        else
          call set_error('source_sector_name is found prior to source start line', &
                       & sub_name)
          return
        endif

      elseif (index(line,'par_str') > 0)then
        if(iStarted /= int_missing)then
!call msg('Source:' + srcInfoTmp(nSrcRead+1)%chSrcNm)
          call count_descriptor_names(line, srcInfoTmp(nSrcRead+1)%nChemDescr, &
                                          & srcInfoTmp(nSrcRead+1)%chDescrNames)
          if(error)then
            call set_error('Failed par_str of the source:' + srcInfoTmp(nSrcRead+1)%chSrcNm + '_' + &
                         & srcInfoTmp(nSrcRead+1)%chSectorNm,sub_name)
            return
          endif
        else
          call set_error('par_str is found prior to source start line', &
                       & sub_name)
          return
        endif

      elseif (index(line,'transport_environment_type') > 0)then
        if(iStarted /= int_missing)then
          iTmp = index(line,'=')
          if(index(fu_str_u_case(line(iTmp+1 :)), 'LAGRANGIAN') > 0)then
            srcInfoTmp(nSrcRead+1)%iDynamicEnvironment = lagrangian_flag
          elseif(index(fu_str_u_case(line(iTmp+1 :)), 'EULERIAN') > 0)then
            srcInfoTmp(nSrcRead+1)%iDynamicEnvironment = eulerian_flag
          else
            call set_error('unknown type of transport environment:'+line,sub_name)
            srcInfoTmp(nSrcRead+1)%iDynamicEnvironment = int_missing
            return
          endif
        else
          call set_error('transport_environment_type is found prior to source start line', &
                       & sub_name)
          return
        endif

      elseif (index (line,'val =') > 0)then
        if(iStarted /= int_missing)then
          if(ifArea)then
            nValLinesTmp = nValLinesTmp + 1
            if(nValLines(iFile) < nValLinesTmp) nValLines(iFile) = nValLinesTmp
          endif
        else
          call set_error('val is found prior to source start line', &
                       & sub_name)
          return
        endif

      elseif (index(line, 'netcdf_data') > 0) then
        if (fu_fails(ifArea, 'netcdf_data outside area source', sub_name)) return
        if (srcInfoTmp(nSrcRead+1)%nChemDescr == 0) then
          !!!!We need already parsed par_str by now
          call msg('Source:' + srcInfoTmp(nSrcRead+1)%chSrcNm)
          call set_error("nChemDescr == 0, probably netcdf_data happened before par_str!", &
              sub_name)

          return
        endif
        nc_entry = adjustl(line(index(line,'=')+1:))
        read(unit=nc_entry, iostat=stat, fmt=*) nctag, size_from_nl, ncfilename
        if (fu_fails(stat == 0, 'Failed to parse data item: ' // trim(nc_entry), &
                   & sub_name)) return
        call replace_string(ncfilename, '^', trim(dirname)//dir_slash)
        call unpack_a_src_from_nc('count', ncfilename, nctag, srcInfoTmp(nSrcRead+1)%nChemDescr, size_from_nc)
        if (error) then
          call set_error("after unpack_a_src_from_nc", sub_name)
          return
        endif
        
        ! The size_from_nl is only used to check for consistency:
        if (fu_fails(size_from_nc == size_from_nl, 'Netcdf/namelist sizes don''t match', &
                   & sub_name)) return
        
        nValLinesTmp = nValLinesTmp + size_from_nc
        if(nValLines(iFile) < nValLinesTmp) nValLines(iFile) = nValLinesTmp

      elseif (index(line,'source_timezone') > 0)then
        
        if(iStarted /= int_missing)then
          iTmp = index(line,'=')
          if (trim(adjustl(fu_str_u_case(line(iTmp+1:)))) == 'LOCAL') &
                    & srcInfoTmp(nSrcRead+1)%ifNeedsTZindexMap = .true.
        else
          call set_error('source_timezone is found prior to source start line', &
                       & sub_name)
          return
        endif

        
      elseif (index(line, 'END_') > 0)then
        !
        ! Source is over. Make its ID and check for duplicates. Should there be such,
        ! just increase the number of descriptors of the previously read source - that 
        ! will be it. Otherwise set the new ID and clean the temporary one
        !
        if(iStarted /= int_missing)then
          srcInfoTmp(nSrcRead+1)%chId = fu_src_id(iSrcIdType, & 
                                                & srcInfoTmp(nSrcRead+1)%chSrcNm, &
                                                & srcInfoTmp(nSrcRead+1)%chSectorNm)

          srcInfoTmp(nSrcRead+1)%iIdNbr = fu_compute_src_id_nbr(srcInfoTmp, nSrcRead+1)

          if(srcInfoTmp(nSrcRead+1)%iIdNbr < 0)then ! old source ID, have to add the newly read source info
            !
            ! Old source ID. Update the information and clean the temporary place
            !
            call add_descriptor_names(srcInfoTmp(-1*srcInfoTmp(nSrcRead+1)%iIdNbr)%nChemDescr, &
                                    & srcInfoTmp(-1*srcInfoTmp(nSrcRead+1)%iIdNbr)%chDescrNames, &
                                    & srcInfoTmp(nSrcRead+1)%nChemDescr, &
                                    & srcInfoTmp(nSrcRead+1)%chDescrNames, &
                                    & chSubstNm, &
                                    & iModeNbr)
            if (error) then
              call msg("Failed to add descriptor"+srcInfoTmp(nSrcRead+1)%chId +","&
                                & +srcInfoTmp(nSrcRead+1)%chSrcNm+","+srcInfoTmp(nSrcRead+1)%chSectorNm)
            endif
            srcInfoTmp(-1*srcInfoTmp(nSrcRead+1)%iIdNbr)%nCells = &
                                            & max(nValLinesTmp, &
                                                & srcInfoTmp(-1*srcInfoTmp(nSrcRead+1)%iIdNbr)%nCells)
            srcInfoTmp(nSrcRead+1)%nChemDescr = 0
            srcInfoTmp(nSrcRead+1)%chDescrNames = ''
            srcInfoTmp(nSrcRead+1)%iIdNbr = 0
          else 
            !
            ! new source ID
            !
            nSrcRead = nSrcRead + 1
            call complete_descriptor_names(srcInfoTmp(nSrcRead)%chDescrNames, &
                                         & chSubstNm, &
                                         & iModeNbr, &
                                         & srcInfoTmp(nSrcRead)%nChemDescr)
            srcInfoTmp(nSrcRead)%iSrcNbr = fu_increment_src_counter(known_src_types(iStarted))
            srcInfoTmp(nSrcRead)%nCells = nValLinesTmp
            if(srcInfoTmp(nSrcRead)%iDynamicEnvironment == int_missing) &
                                    & srcInfoTmp(nSrcRead)%iDynamicEnvironment = iDynamicsTypeDefault
            if(size(srcInfoTmp) <= nSrcRead)call enlarge_source_info(srcInfoTmp)
            if(error)return
          endif ! old or new source ID
          iStarted = int_missing
        else
          call set_error('End of the source without start:' + line, sub_name)
          return
        endif
        chSubstNm = ''
      endif  ! content of the line

    end do   ! while .not. eof

    ! Close the possibly open netcdf data file:
    call unpack_a_src_from_nc('close', '', int_missing, int_missing, size_from_nc)

    CONTAINS
                                    
      integer function fu_increment_src_counter(src_summary) result(cnt)
        implicit none
        type(Tsource_summary), intent(in) :: src_summary
        
        select case(src_summary%source_type)
          case(area_source)
            em_source%n_area = em_source%n_area + 1
            cnt = em_source%n_area
          case(bomb_source)
            em_source%n_bomb = em_source%n_bomb + 1
            cnt = em_source%n_bomb
          case(dms_source)
            em_source%n_dms = em_source%n_dms + 1
            cnt = em_source%n_dms
          case(fire_source)
            em_source%n_fire = em_source%n_fire + 1
            cnt = em_source%n_fire
          case(point_source)
            em_source%n_point = em_source%n_point + 1
            cnt = em_source%n_point
          case(pollen_source)
            em_source%n_pollen = em_source%n_pollen + 1
            cnt = em_source%n_pollen
          case(bio_voc_source)
            em_source%n_bio_voc = em_source%n_bio_voc + 1
            cnt = em_source%n_bio_voc
          case(sea_salt_source)
            em_source%n_sea_salt = em_source%n_sea_salt + 1
            cnt = em_source%n_sea_salt
          case(volcano_source)
            em_source%n_volc = em_source%n_volc + 1
            cnt = em_source%n_volc
          case(wind_blown_dust_source)
            em_source%n_wb_dust =em_source%n_wb_dust + 1
            cnt = em_source%n_wb_dust
          case default
            call set_error('Unknown source type:' + src_summary%label, 'fu_increment_src_counter')
            cnt = int_missing
        end select
      end function fu_increment_src_counter
    
  end subroutine count_srcs_one_data_file


  !********************************************************************************************  

  subroutine enlarge_source_info(pSrcInfo)
    !
    ! Enlarges the size of temporary source info structure
    !
    implicit none

    ! Imported parameters
    type(silam_source_info), dimension(:), allocatable, intent(inout) :: pSrcInfo

    ! Local parameters
    integer :: iTmp, nOld, nNew
    type(silam_source_info), dimension(:), allocatable :: srcInfoTmp

    if(allocated(pSrcInfo))then !  ! Have to store the existing information !
      nOld = size(pSrcInfo)
      nNew = nint(nOld*1.5)

      allocate(srcInfoTmp(nNew), stat=iTmp)
      if(fu_fails(iTmp==0,'Failed allocation of temporary source info','enlarge_source_info'))return

      do iTmp = 1, nOld
        srcInfoTmp(iTmp) = pSrcInfo(iTmp)
      enddo

      deallocate(pSrcInfo, stat=iTmp)
      if(fu_fails(iTmp==0,'Failed deallocation of source info','enlarge_source_info')) return

      call move_alloc(srcInfoTmp, pSrcInfo)

    else !  ! Initial allocation !
      nOld = 0
      nNew = 200
      allocate(pSrcInfo(200), stat = iTmp)
      if(fu_fails(iTmp==0,'Failed initial allocation temporary source info','enlarge_source_info'))return

    endif  ! associated pSrcInfo
    !
    ! Nullify the non-existing entries
    do iTmp = nOld+1, nNew
      pSrcInfo(iTmp) = silam_source_info_missing
      pSrcInfo(iTmp)%nChemDescr = 0
      pSrcInfo(iTmp)%nCells = 0
    enddo

  end subroutine enlarge_source_info


  !*********************************************************************

  subroutine allocate_source_space(source)
    !
    ! Allocates the necessary structures for the source. 
    ! The number of specific sources must be filled in
    !
    implicit none

    ! Input and output value of the function
    TYPE(silam_source), intent(inout) :: source

    ! Local variables
    integer :: iTmp, iStatus, i_point, i_area, i_bomb, i_fire, i_sea_salt, i_pollen, &
             & i_bio_voc, i_wb_dust, i_dms, i_volc

    !
    ! Stupidity check
    !
    iTmp = fu_num_sources_total(source)
    if( iTmp <= 0 .or. iTmp> 2000000)then
      call msg('Unrealistic number of sources:', iTmp)
      call set_error('Unrealistic number of sources','allocate_source_space')
      return
    endif

    if(source%n_area > 0) then                               !============= AREA
      call msg('Allocating area sources:',source%n_area)
      allocate(source%a_ptr(source%n_area),stat=iStatus)
      if(fu_fails(iStatus == 0,'Failed to allocate area sources','allocate_source_space'))return
      do iTmp=1,source%n_area
          call set_undef(source%a_ptr(iTmp)%a_src)
      enddo
    else
      nullify(source%a_ptr)
    endif
    if(source%n_bio_voc > 0) then                            !============ BIO-VOC
      call msg('Allocating bio VOC sources:',source%n_bio_voc)
      allocate(source%bvoc_ptr(source%n_bio_voc),stat=iStatus)
      if(fu_fails(iStatus == 0,'Failed to allocate bio VOC sources','allocate_source_space'))return
    else
      nullify(source%bvoc_ptr)
    endif
    if(source%n_bomb > 0) then                               !============= BOMB
      call msg('Allocating bomb sources: ',source%n_bomb)
      allocate(source%b_ptr(source%n_bomb),stat=iStatus)
      if(fu_fails(iStatus == 0,'Failed to allocate bomb sources','allocate_source_space'))return
    else
      nullify(source%b_ptr)
    endif
    if(source%n_fire > 0) then                               !============= FIRE
      call msg('Allocating fire sources: ',source%n_fire)
      allocate(source%fire_ptr(source%n_fire),stat=iStatus)
      if(fu_fails(iStatus == 0,'Failed to allocate fire sources','allocate_source_space'))return
    else
      nullify(source%fire_ptr)
    endif
    if(source%n_point > 0)  then                             !============ POINT
      call msg('Allocating point sources: ',source%n_point)
      allocate(source%p_ptr(source%n_point),stat=iStatus)
      if(fu_fails(iStatus == 0,'Failed to allocate point sources','allocate_source_space'))return
    else
      nullify(source%p_ptr)
    endif
    if(source%n_pollen > 0) then                             !============ POLLEN
      call msg('Allocating pollen sources:',source%n_pollen)
      allocate(source%pollen_ptr(source%n_pollen),stat=iStatus)
      if(fu_fails(iStatus == 0,'Failed to allocate pollen sources','allocate_source_space'))return
    else
      nullify(source%pollen_ptr)
    endif
    if(source%n_sea_salt > 0) then                           !============= SEA SALT
      call msg('Allocating sea salt sources:',source%n_sea_salt)
      allocate(source%sslt_ptr(source%n_sea_salt),stat=iStatus)
      if(fu_fails(iStatus == 0,'Failed to allocate sea salt sources','allocate_source_space'))return
    else
      nullify(source%sslt_ptr)
    endif
    if(source%n_wb_dust > 0) then                           !============= Wind-BLOWN DUST
      call msg('Allocating wind-blown dust sources:',source%n_wb_dust)
      allocate(source%wbdust_ptr(source%n_wb_dust),stat=iStatus)
      if(fu_fails(iStatus == 0,'Failed to allocate wind-blown dust sources','allocate_source_space'))return
    else
      nullify(source%wbdust_ptr)
    endif
    if(source%n_dms > 0) then                           !============= DMS
      call msg('Allocating DMS dust sources:',source%n_dms)
      allocate(source%dms_ptr(source%n_dms),stat=iStatus)
      if(fu_fails(iStatus == 0,'Failed to allocate DMS sources','allocate_source_space'))return
    else
      nullify(source%dms_ptr)
    endif
    if(source%n_volc > 0) then                           !============= VOLCANO
      call msg('Allocating volcano sources:',source%n_volc)
      allocate(source%volc_ptr(source%n_volc),stat=iStatus)
      if(fu_fails(iStatus == 0,'Failed to allocate volcano sources','allocate_source_space'))return
    else
      nullify(source%volc_ptr)
    endif

    !
    ! Go through the sources allocating whatever is needed and setting the 
    ! information fields
    !
    i_area = 0
    i_bio_voc = 0
    i_bomb = 0
    i_fire = 0
    i_point = 0
    i_pollen = 0
    i_sea_salt = 0
    i_wb_dust = 0
    i_dms = 0
    i_volc = 0
    do iTmp = 1, size(source%src_info_lst)  !source%n_point + source%n_area + source%n_bomb
      select case(source%src_info_lst(iTmp)%iSrcType)
        case(area_source)
          i_area = i_area + 1
          call reserve_area_source(source%a_ptr(i_area)%a_src, &     ! Src to initialise
                                 & source%src_info_lst(iTmp)%iSrcNbr, & ! Src number
                                 & source%src_info_lst(iTmp)%iIdNbr, &  ! SrcID number
                                 & source%src_info_lst(iTmp)%nChemDescr, & ! Nbr of chemical descr to reserve
                                 & source%src_info_lst(iTmp)%nCells)    ! max number of cells found so far
        case(bio_voc_source)
          i_bio_voc = i_bio_voc + 1
          call reserve_bio_voc_source(source%bvoc_ptr(i_bio_voc)%bvoc_src, &     ! Src to initialise
                                    & source%src_info_lst(iTmp)%iSrcNbr, & ! Src number
                                    & source%src_info_lst(iTmp)%iIdNbr)  ! SrcID number
        case(bomb_source)
          i_bomb = i_bomb + 1
          call reserve_bomb_source(source%b_ptr(i_bomb)%b_src, &     ! Src to initialise
                                 & source%src_info_lst(iTmp)%iSrcNbr, & ! Src number
                                 & source%src_info_lst(iTmp)%iIdNbr)  ! SrcID number
        case(fire_source)
          i_fire = i_fire + 1
          call reserve_fire_source(source%fire_ptr(i_fire)%fire_src, &     ! Src to initialise
                                 & source%src_info_lst(iTmp)%iSrcNbr, & ! Src number
                                 & source%src_info_lst(iTmp)%iIdNbr)  ! SrcID number
        case(point_source)
          i_point = i_point + 1
          call reserve_point_source(source%p_ptr(i_point)%p_src, &     ! Src to initialise
                                  & source%src_info_lst(iTmp)%iSrcNbr, & ! Src number
                                  & source%src_info_lst(iTmp)%iIdNbr, &  ! SrcID number
                                  & source%src_info_lst(iTmp)%nChemDescr, & ! Nbr of chemical descr to reserve
                                  & source%src_info_lst(iTmp)%iDynamicEnvironment) ! Euleriean or Lagrangian
        case(pollen_source)
          i_pollen = i_pollen + 1
          call reserve_pollen_source(source%pollen_ptr(i_pollen)%pollen_src, &     ! Src to initialise
                                   & source%src_info_lst(iTmp)%iSrcNbr, & ! Src number
                                   & source%src_info_lst(iTmp)%iIdNbr)  ! SrcID number
        case(sea_salt_source)
          i_sea_salt = i_sea_salt + 1
          call reserve_sea_salt_source(source%sslt_ptr(i_sea_salt)%sslt_src, &     ! Src to initialise
                                     & source%src_info_lst(iTmp)%iSrcNbr, & ! Src number
                                     & source%src_info_lst(iTmp)%iIdNbr)  ! SrcID number
        case(wind_blown_dust_source)
          i_wb_dust = i_wb_dust + 1
          call reserve_wb_dust_source(source%wbdust_ptr(i_wb_dust)%wbdust_src, &     ! Src to initialise
                                     & source%src_info_lst(iTmp)%iSrcNbr, & ! Src number
                                     & source%src_info_lst(iTmp)%iIdNbr)  ! SrcID number
        case(dms_source)
          i_dms = i_dms + 1
          call reserve_dms_source(source%dms_ptr(i_dms)%dms_src, &     ! Src to initialise
                                & source%src_info_lst(iTmp)%iSrcNbr, & ! Src number
                                & source%src_info_lst(iTmp)%iIdNbr)  ! SrcID number
        case(volcano_source)
          i_volc = i_volc + 1
          call reserve_volcano_source(source%volc_ptr(i_volc)%volc_src, &     ! Src to initialise
                                    & source%src_info_lst(iTmp)%iSrcNbr, & ! Src number
                                    & source%src_info_lst(iTmp)%iIdNbr, &  ! SrcID number
                                    & source%src_info_lst(iTmp)%nChemDescr, & ! Nbr of chemical descr to reserve
                                    & source%src_info_lst(iTmp)%iDynamicEnvironment) ! Euleriean or Lagrangian
        case(int_missing)
          exit  ! all done

        case default
          call msg('Unknown source type:',source%src_info_lst(iTmp)%iSrcType)
          call set_error('Unknown source type:','allocate_source_space')
          return
      end select
    end do  ! all sources

  end subroutine allocate_source_space


  ! ***************************************************************

  subroutine read_sources_from_file(file_name, source, iSrcIdType, nValLines, iFile, expected_species)
    !
    ! The function reads all SILAM sources from the given file.
    ! First, the file is scanned checking how many sources are to be read,
    ! then they are consumed one-by-one
    ! The source type is determined from the first source line.
    !
    ! When the sources are read, necessary intialisation is performed
    ! for cocktails - chemical or radioactive. Here it is not the best
    ! place for it, but it is virtually the only place where all sources
    ! are available at once and the pollution cloud is not yet initialised.
    !
    IMPLICIT NONE

    ! Imported parameters
    character(len=*), intent(in) :: file_name
    TYPE(silam_source), intent(inout) :: source
    integer, intent(in) :: iSrcIdType, iFile  
    integer, dimension(:), pointer :: nValLines
    type(silam_species), dimension(:), intent(in), allocatable :: expected_species

    ! Local variables
!    INTEGER :: status, i_point, i_area, i_bomb, iGlobalIndex, iStatus
    INTEGER :: iGlobalIndex = 1, status, file_unit, dummy
    LOGICAL :: eof
    CHARACTER(len=clen) :: line
    character(len=80) :: chLabel
    type(Tsilam_namelist), pointer :: nlSrc
    character(len=fnlen) :: chDataDir

    file_unit = fu_next_free_unit()
    if(error)return
    OPEN(file_unit, file=file_name, action='read', status='old', iostat=status)
    IF (status /= 0) THEN
      CALL set_error('cannot open data file:' + file_name,'set_source_terms')
      RETURN
    END IF
    call msg('Reading the data file: ' + file_name)

    !
    ! Read the file through. When each source is consumed, fill-in the src_info_lst element
    ! for the cross-type info freely available afterwards
    !
    !rewind(file_unit)
    line = ""
    eof=.false.
    source%iSrcIdType = iSrcIdType

    FILE_MAIN: do while(.not.eof)

      CALL next_line_from_input_file(file_unit, line, eof)
      if(error .or. eof)exit FILE_MAIN

      do while(index(line,'SOURCE_') ==0)
        CALL next_line_from_input_file(file_unit, line, eof)
        IF (error.or.eof) exit FILE_MAIN
      end do

      chLabel = line
      line = fu_str_u_case(line)

!call msg('Reading the source:'+chLabel)

      IF(index(line,'AREA_SOURCE_') == 1)THEN

!        call msg('Reading area source from file:' + trim(fname))

        nlSrc => fu_read_namelist(file_unit, .false., 'END_'+chLabel, nItems=nValLines(iFIle)+50)
        if(error)return

      ELSEIF (index(line,'BOMB_SOURCE_') == 1 .or. &
            & index(line,'FIRE_SOURCE_') == 1 .or. &
            & index(line,'POINT_SOURCE_') == 1 .or. &
            & index(line,'SEA_SALT_SOURCE_') == 1 .or. &
            & index(line,'WIND_BLOWN_DUST_SOURCE_') == 1 .or. &
            & index(line,'POLLEN_SOURCE_') == 1 .or. &
            & index(line,'BIOGENIC_VOC_SOURCE_') == 1 .or. &
            & index(line, 'DMS_SOURCE_') == 1 .or. &
            & index(line, 'VOLCANO_SOURCE_') == 1)THEN

!        call msg('Reading the:' + chLabel + '- source from file:' + trim(fname))

        nlSrc => fu_read_namelist(file_unit, .false., 'END_'+chLabel)
        if(error)return

      elseIF (index(line,'END_') > 0) THEN

      elseif(eof) then
        exit

      ELSE
        CALL set_error('Unknown sources type found:' + line, 'read_sources_from_file')
        close(file_unit)
        RETURN
      END IF

      if(error)then 
        close(file_unit)
        return
      endif

      !call report(nlSrc)
      !call msg('');call msg('')
      
      ! Area source can point to netcdf file, whose path can be given as relative. By
      ! default, it is relative to the main file but can be overwritten
      !
      if(len_trim(fu_content(nlsrc, 'netcdf_dir')) > 0)then
        call set_error('netcdf_dir is obsolete, use data_dir instead','read_sources_from_file')
        call unset_error('read_sources_from_file')
        chDataDir = fu_process_filepath(fu_content(nlsrc, 'netcdf_dir')) 
      else
        chDataDir = fu_process_filepath(fu_content(nlsrc, 'data_dir'))
      endif
      if(len_trim(fu_content(nlsrc, 'data_dir')) > 0) chDataDir = chDataDir + dir_slash
      if(len_trim(fu_dirname(file_name)) > 0) &
            & chDataDir = chDataDir + fu_process_filepath(fu_dirname(file_name)) + dir_slash
      
      !
      ! Namelist is ready. Set the source
      !
      call set_next_gen_src_from_nl(source, nlSrc, file_name, chDataDir, iGlobalIndex, chLabel, &
                                  & expected_species)

      !
      ! This source is handled, whether successfully or not. Destroy the namelist
      !
      call destroy_namelist(nlSrc)
!      call reset_namelist(nlSrc)
      if(error)return
      
      source%defined = silja_true

    end do FILE_MAIN   ! Scan through the file
    
    close(file_unit)
    
    if(error)return
    
    ! Close the possibly open netcdf data file:
    call unpack_a_src_from_nc('close', '', int_missing, int_missing, dummy)
    
  END subroutine read_sources_from_file


  !*****************************************************************************

  recursive subroutine set_next_gen_src_from_nl(source, nlSrc, chSrcFNm, chDataDir_, iSrcIndex, &
                                              & chLabel, expected_species)
    !
    ! Gets a namelist, sets the source index in the info_list (which is already largely filled-in)
    ! by calling the non-recursive functions for individual sources. The reason for
    ! existence of this function is that the sources can be spread over several files, with
    ! the main file containing just names of these files. Since we have to decide upon the
    ! place for the source BEFORE we actually set it, we have to use this function to do it.
    !
    implicit none

    ! imported parameters
    TYPE(silam_source), intent(inout) :: source
    type(Tsilam_namelist), pointer :: nlSrc
    integer, intent(inout) :: iSrcIndex  ! 1 .. n_point+n_area+n_bomb
    character(len=*), intent(in) :: chLabel, chSrcFNm, chDataDir_
    type(silam_species), dimension(:), intent(in), allocatable :: expected_species

    ! Local variables
    INTEGER :: iLocal, iStatus, iTmp
    character(len=clen) :: chSourceVersion
    character(len=fnlen) :: chDataDir
    type(Tsilam_namelist), pointer :: nlSrcLocal
    type(silam_sp) :: spNm, spSector, spContent

    !
    ! Area source can point to netcdf or grads file, whose path can be given as relative. By
    ! default, it is relative to the main file but can be overwritten via data_dir line
    !
    if(len_trim(fu_content(nlsrc, 'data_dir')) > 0)then
      chDataDir = fu_process_filepath(fu_content(nlsrc,'data_dir')) + dir_slash
    else
      chDataDir = chDataDir_
    endif
    !
    ! A trick: it may happen that the source is, in fact, a separate file
    ! Read it then
    !
    spContent%sp => fu_work_string()
    if(error)return
    spContent%sp = fu_content(nlSrc,'source_file_name')
    if(spContent%sp /= '')then
      iLocal = fu_next_free_unit()
      open(unit=iLocal, file=spContent%sp, status='old', iostat=iStatus)
      if(iStatus == 0)then
        call msg('Reading the source from file:' + spContent%sp + ', data directory:' + chDataDir)
        
        nlSrcLocal => fu_read_namelist(iLocal, .false., 'END_'+chLabel)
        if(error)return
        call set_next_gen_src_from_nl(source, nlSrcLocal, chSrcFNm, chDataDir, &
                                    & iSrcIndex, chLabel, expected_species)
        
        close(iLocal)
        call destroy_namelist(nlSrcLocal)
        call free_work_array(spContent%sp)
        return
      else
        call set_error('Source file does not exist:' + spContent%sp, 'set_next_gen_src_from_nl')
        return
      endif
    endif
    !
    ! OK, the namelist must be reasonable, at least it does not contain the name of a
    ! next-in-chain file. Compute the source ID, its number and then set the source itself
    !
    spNm%sp => fu_work_string()
    spSector%sp => fu_work_string()
    if(error)return

    spNm%sp = fu_content(nlSrc,'source_name')
    spSector%sp = fu_content(nlSrc,'source_sector_name')

    do iTmp = 1, size(source%src_info_lst)
      ! source%n_area + source%n_bomb + source%n_point + source%n_pollen + source%n_sea_salt
      iSrcIndex = int_missing
      if(source%src_info_lst(iTmp)%chSrcNm == spNm%sp .and. &
       & source%src_info_lst(iTmp)%chSectorNm == spSector%sp)then
        iSrcIndex = iTmp
        exit
      endif
    end do
    if(iSrcIndex == int_missing)then
      call set_error('Failed to find the source ID for:' + spNm%sp + ',' + spSector%sp, &
                   & 'set_next_gen_src_from_nl')
      return
    endif

!    call msg('Filling-in to the source:' + &
!           & source%src_info_lst(iSrcIndex)%chSrcNm + ',' + &
!           & source%src_info_lst(iSrcIndex)%chSectorNm)

    iTmp = index(chLabel,"_SOURCE_")+8  !! After it
    chSourceVersion = chLabel(iTmp:)

    !
    ! This index points at the existing point source, to which the information must be added.
    ! Actually, it points to the index in the src_info_lst, which iSrcNbr points to
    ! the source index in the p_src array:
    !
    if(source%src_info_lst(iSrcIndex)%iSrcType == area_source)then
      !
      ! Area source
      !
      call fill_a_src_from_namelist(nlSrc, &
                                  & source%a_ptr(source%src_info_lst(iSrcIndex)%iSrcNbr)%a_src, &
                                  & chSourceVersion, &
                                  & chDataDir)

    elseif(source%src_info_lst(iSrcIndex)%iSrcType == bio_voc_source)then
      !
      ! Bio VOC source. Note somewhat different initialisation procedure
      !
      call fill_bio_voc_src_from_namelist(&
                           & nlSrc, &
                           & source%bvoc_ptr(source%src_info_lst(iSrcIndex)%iSrcNbr)%bvoc_src, &
                           & chDataDir)

    elseif(source%src_info_lst(iSrcIndex)%iSrcType == bomb_source)then
      !
      ! Bomb source
      !
      call fill_b_src_from_namelist( &
                           & nlSrc, &
                           & source%b_ptr(source%src_info_lst(iSrcIndex)%iSrcNbr)%b_src, &
                           & chSourceVersion )

    elseif(source%src_info_lst(iSrcIndex)%iSrcType == fire_source)then
      !
      ! Fire source
      !
      call fill_fire_src_from_namelist( &
                           & nlSrc, &
                           & source%fire_ptr(source%src_info_lst(iSrcIndex)%iSrcNbr)%fire_src, &
                           & expected_species, chDataDir)

    elseif(source%src_info_lst(iSrcIndex)%iSrcType == point_source)then
      !
      ! Point source
      !
      call fill_p_src_from_namelist( &
                           & nlSrc, &
                           & source%p_ptr(source%src_info_lst(iSrcIndex)%iSrcNbr)%p_src, &
                           & chSourceVersion)

    elseif(source%src_info_lst(iSrcIndex)%iSrcType == pollen_source)then
      !
      ! Pollen source. Note somewhat different initialisation procedure
      !
      call fill_pollen_src_from_namelist(&
                           & nlSrc, &
                           & source%pollen_ptr(source%src_info_lst(iSrcIndex)%iSrcNbr)%pollen_src, &
                           & chDataDir)

    elseif(source%src_info_lst(iSrcIndex)%iSrcType == sea_salt_source)then
      !
      ! Sea salt source. Note somewhat different initialisation procedure
      !
      call fill_sslt_src_from_namelist(&
                           & nlSrc, &
                           & source%sslt_ptr(source%src_info_lst(iSrcIndex)%iSrcNbr)%sslt_src, &
                           & expected_species, &
                           & chDataDir)

    elseif(source%src_info_lst(iSrcIndex)%iSrcType == wind_blown_dust_source)then
      !
      ! Wind-blown dust source. Note somewhat different initialisation procedure
      !
      call fill_wb_dust_src_from_namelist(&
                           & nlSrc, &
                           & source%wbdust_ptr(source%src_info_lst(iSrcIndex)%iSrcNbr)%wbdust_src, &
                           & expected_species, &
                           & chDataDir)

    elseif(source%src_info_lst(iSrcIndex)%iSrcType == dms_source)then
      !
      ! DMS source.
      !
      call fill_dms_src_from_namelist(nlSrc, source%dms_ptr(source%src_info_lst(iSrcIndex)%iSrcNbr)%dms_src, chDataDir)

    elseif(source%src_info_lst(iSrcIndex)%iSrcType == volcano_source)then
      !
      ! Volcano source
      !
      call fill_volc_src_from_namelist(nlSrc, &
                                     & source%volc_ptr(source%src_info_lst(iSrcIndex)%iSrcNbr)%volc_src)
      
    else  ! =========== unknown chLabel

      call set_error('Unknown source starting label:' + chLabel, 'set_next_gen_src_from_nl')

    endif  ! content of chLabel

    call free_work_array(spContent%sp)
    call free_work_array(spNm%sp)
    call free_work_array(spSector%sp)

  end subroutine set_next_gen_src_from_nl

  
  !****************************************************************

  subroutine add_source_term_input_needs(src, q_met_dynamic, q_met_static, &
                                            & q_disp_dynamic, q_disp_static, wdr)
    !
    ! Reports the specific needs of the source term concerning the
    ! meteorological data. For example, plume rise will need temperature,
    ! wind speed and ABL parameters. Sea salt needs other stuff, etc.
    !
    implicit none

    ! Imported parameters 
    type(silam_source), intent(in) :: src
    integer, dimension(:), intent(inout) :: q_met_dynamic, q_met_static, &
                                          & q_disp_dynamic, q_disp_static
    type(silja_wdr), intent(in) :: wdr
    ! Local variables
    integer :: iSrc

    !
    ! Area sources might need tz_offset, and, eventually other meteo
    !
    do iSrc = 1, src%n_area
      call add_input_needs_area_source(src%a_ptr(iSrc)%a_src, q_met_dynamic, q_met_static, &
                                     & q_disp_dynamic, q_disp_static)
      if(error)return
    end do
    !
    ! Bio VOC sources depend on actual meteodata
    !
    do iSrc = 1, src%n_bio_voc
      call add_input_needs(src%bvoc_ptr(iSrc)%bvoc_src, q_met_dynamic, q_met_static, &
                                                      & q_disp_dynamic, q_disp_static)
      if(error)return
    end do
    !
    ! Bomb sources needs land mask
    !
    do iSrc = 1, src%n_bomb
      call add_input_needs(src%b_ptr(iSrc)%b_src, q_met_dynamic, q_met_static, &
                                                & q_disp_dynamic, q_disp_static)
      if(error)return
    end do
    !
    ! Fire sources do have this need
    !
    do iSrc = 1, src%n_fire
      call add_input_needs(src%fire_ptr(iSrc)%fire_src, q_met_dynamic, q_met_static, &
                                                      & q_disp_dynamic, q_disp_static)
      if(error)return
    end do
    !
    ! Point sources can have plume rise - at least some of them
    !
    do iSrc = 1, src%n_point
      call add_input_needs_point_source(src%p_ptr(iSrc)%p_src, q_met_dynamic, q_met_static, &
                                                & q_disp_dynamic, q_disp_static)
      if(error)return
    end do
    !
    ! Pollen sources depend on actual meteodata and also requires fields in the dispersion buffer
    !
    do iSrc = 1, src%n_pollen
      call add_input_needs(src%pollen_ptr(iSrc)%pollen_src, q_met_dynamic, q_met_static, &
                                                          & q_disp_dynamic, q_disp_static)
      if(error)return
    end do
    !
    ! Sea salt sources depend on actual meteodata
    !
    do iSrc = 1, src%n_sea_salt
      call add_input_needs(src%sslt_ptr(iSrc)%sslt_src, q_met_dynamic, q_met_static, &
                                                      & q_disp_dynamic, q_disp_static)
      if(error)return
    end do
    !
    ! Wind-blown dust sources depend on actual meteodata
    !
    do iSrc = 1, src%n_wb_dust
      call add_input_needs(src%wbdust_ptr(iSrc)%wbdust_src, q_met_dynamic, q_met_static, &
                                                          & q_disp_dynamic, q_disp_static, wdr)
      if(error)return
    end do
    !
    ! DMS sources depend on actual meteodata
    !
    do iSrc = 1, src%n_dms
      call add_input_needs(src%dms_ptr(iSrc)%dms_src, q_met_dynamic, q_met_static, &
                         & q_disp_dynamic, q_disp_static)
      if(error)return
    end do
    !
    ! Volcano sources depend on actual meteodata 
    !
    do iSrc = 1, src%n_volc
      call add_input_needs_volcano_source(src%volc_ptr(iSrc)%volc_src, q_met_dynamic, q_met_static, &
                                                                     & q_disp_dynamic, q_disp_static)
      if(error)return
    end do

    if(src%iSrcIdType == iSrcTimeZoneId)iSrc = fu_merge_integer_to_array(timezone_index_flag, q_met_Static)

    
  end subroutine add_source_term_input_needs


  !************************************************************************************

  subroutine verify_sources(em_source, timestart, timestep)
    !
    ! Checks that the setup is OK for the source. The only limitation so far is 
    ! area source with time-resolving field in binary. Check that.
    !
    implicit none
    
    ! Imported parameters
    type(silja_time), intent(in) :: timestart
    type(silja_interval), intent(in) :: timestep
    type(silam_source), intent(in) :: em_source
    
    ! Local parameters
    integer :: iSrc

    call msg_warning("check_time_params_a_src missing!", "verify_sources")
    return
    do iSrc = 1, em_source%n_area
    !!      call check_time_params_a_src(em_source%a_ptr(iSrc)%a_src, timestart, timestep)
      if(error)return
    end do
    
  end subroutine verify_sources

  !************************************************************************************

  subroutine get_inventory_g_src(src, dynamics_type, species_list, nspecies)
    !
    ! Collect the lists of emission species from all sources, merge
    ! them and return the full emission species list.
    !
    implicit none

    ! Imported parameters
    type(silam_source), intent(in) :: src
    integer, intent(in) :: dynamics_type     ! Lagrangian / Eulerian / int_missing
    type(silam_species), dimension(:), pointer :: species_list
    integer, intent(out) :: nspecies

    ! Local variables
    integer :: i

    if (.not. defined(src)) then
      call set_error('Source not defined', 'get_inventory_g_src')
      return
    end if

    nSpecies = 0
    nullify(species_list)

    if(dynamics_type == int_missing)then
      do i = 1, src%n_area                                    !======== AREA
        call add_inventory_a_src(src%a_ptr(i)%a_src, species_list, nSpecies)
      end do
      do i = 1, src%n_bio_voc                                  !======== Bio VOC
        call add_inventory_bio_voc_src(src%bvoc_ptr(i)%bvoc_src, species_list, nSpecies)
      end do
      do i = 1, src%n_bomb                                    !======== BOMB
        call add_source_species_b_src(src%b_ptr(i)%b_src, species_list, nSpecies)
      end do
      do i = 1, src%n_fire                                    !======== FIRE
        call add_inventory_fire_src(src%fire_ptr(i)%fire_src, species_list, nSpecies)
      end do
      do i = 1, src%n_point                                   !======== POINT
        call add_inventory_p_src(src%p_ptr(i)%p_src, species_list, nSpecies)
      end do
      do i = 1, src%n_pollen                                  !======== POLLEN
        call add_inventory_pollen_src(src%pollen_ptr(i)%pollen_src, species_list, nSpecies)
      end do
      do i = 1, src%n_sea_salt                                !======== SEA SALT
        call add_inventory_sslt_src(src%sslt_ptr(i)%sslt_src, species_list, nSpecies)
      end do
      do i = 1, src%n_wb_dust                                 !======== WIND-BLOWN DUST
        call add_inventory_wb_dust_src(src%wbdust_ptr(i)%wbdust_src, species_list, nSpecies)
      end do
      do i = 1, src%n_dms
        call add_inventory_dmssrc(src%dms_ptr(i)%dms_src, species_list, nSpecies)
      end do
      do i = 1, src%n_volc
        call add_inventory_volc_src(src%volc_ptr(i)%volc_src, species_list, nSpecies)
      end do
    else
      do i = 1, src%n_area                                    !======== AREA
        if(dynamics_type == &
            & src%src_info_lst(fu_source_nbr(src%a_ptr(i)%a_src))%iDynamicEnvironment) &
                            &  call add_inventory_a_src(src%a_ptr(i)%a_src, species_list, nSpecies)
      end do
      do i = 1, src%n_bio_voc                                  !======== Bio VOC
        if(dynamics_type == &
            & src%src_info_lst(fu_source_nbr(src%bvoc_ptr(i)%bvoc_src))%iDynamicEnvironment) &
                            &  call add_inventory_bio_voc_src(src%bvoc_ptr(i)%bvoc_src, species_list, nSpecies)
      end do
      do i = 1, src%n_bomb                                    !======== BOMB
        if(dynamics_type == &
            & src%src_info_lst(fu_source_nbr(src%b_ptr(i)%b_src))%iDynamicEnvironment) &
                            & call add_source_species_b_src(src%b_ptr(i)%b_src, species_list, nSpecies)
      end do
      do i = 1, src%n_fire                                    !======== FIRE
        if(dynamics_type == &
            & src%src_info_lst(fu_source_nbr(src%fire_ptr(i)%fire_src))%iDynamicEnvironment) &
                            & call add_inventory_fire_src(src%fire_ptr(i)%fire_src, species_list, nSpecies)
      end do
      do i = 1, src%n_point                                   !======== POINT
        if(dynamics_type == &
            & src%src_info_lst(fu_source_nbr(src%p_ptr(i)%p_src))%iDynamicEnvironment) &
                            & call add_inventory_p_src(src%p_ptr(i)%p_src, species_list, nSpecies)
      end do
      do i = 1, src%n_pollen                                  !======== POLLEN
        if(dynamics_type == &
            & src%src_info_lst(fu_source_nbr(src%pollen_ptr(i)%pollen_src))%iDynamicEnvironment) &
                            & call add_inventory_pollen_src(src%pollen_ptr(i)%pollen_src, species_list, nSpecies)
      end do
      do i = 1, src%n_sea_salt                                !======== SEA SALT
        if(dynamics_type == &
            & src%src_info_lst(fu_source_nbr(src%sslt_ptr(i)%sslt_src))%iDynamicEnvironment) &
                            & call add_inventory_sslt_src(src%sslt_ptr(i)%sslt_src, species_list, nSpecies)
      end do
      do i = 1, src%n_wb_dust                                 !======== WIND-BLOWN DUST
        if(dynamics_type == &
            & src%src_info_lst(fu_source_nbr(src%wbdust_ptr(i)%wbdust_src))%iDynamicEnvironment) &
                            & call add_inventory_wb_dust_src(src%wbdust_ptr(i)%wbdust_src, species_list, nSpecies)
      end do
      do i = 1, src%n_dms
        if(dynamics_type == &
            & src%src_info_lst(fu_source_nbr(src%dms_ptr(i)%dms_src))%iDynamicEnvironment) &
                            & call add_inventory_dmssrc(src%dms_ptr(i)%dms_src, species_list, nSpecies)
      end do
      do i = 1, src%n_volc
        if(dynamics_type == &
            & src%src_info_lst(fu_source_nbr(src%volc_ptr(i)%volc_src))%iDynamicEnvironment) &
                            & call add_inventory_volc_src(src%volc_ptr(i)%volc_src, species_list, nSpecies)
      end do
    endif  ! if dynamics type is selective

  end subroutine get_inventory_g_src


  !******************************************************************

  function fu_src_id(iIdType, chSrcNm, chSectorNm) result(id)
    !
    ! Computes the source id from its name and sector
    !
    implicit none

    ! Return variable
    character(len=clen) :: id

    ! Imported parameters
    integer, intent(in) :: iIdType
    character(len=*), intent(in) :: chSrcNm, chSectorNm

    select case(iIdType)
      case(iNoId)
        id = 'ALL_SRCS'
      case(iSrcNmId)
        id = chSrcNm
      case(iSrcSectorId)
        id = chSectorNm
      case(iSrcNmSectorId)
        if(len_trim(chSrcNm) + len_trim(chSectorNm) > clen)then
          call msg('Source name:' + chSrcNm)
          call msg('Sector name:' + chSectorNm)
          call set_error('Too long source and sector names','fu_src_id')
          return
        endif
        id = chSrcNm + '_' + chSectorNm
      case(iSrcTimeZoneId)
        id = 'time_zone'
      case default
        call msg('Strange source id type:',iIdType)
        call set_error('Strange source id type','fu_src_id')
        return
    end select

  end function fu_src_id

  !****************************************************************

  integer function fu_id_nbr_from_srcNbr(src, srcNbr)
    !
    ! Finds out the id number having its number in the total list
    !
    implicit none

    ! Imported parameters
    type(silam_source), intent(in) :: src
    integer, intent(in) :: srcNbr

    if(srcNbr <= 0 .or. srcNbr > size(src%src_info_lst))then
      call msg('Strange source number given:', srcNbr)
      call set_error('Strange source number given','fu_id_nbr_from_srcNbr')
      return
    endif

    fu_id_nbr_from_srcNbr = src%src_info_lst(srcNbr)%iIdNbr

  end function fu_id_nbr_from_srcNbr


  !******************************************************************

  integer function fu_id_nbr_from_srcNms(src, chSrcNm, chSectorNm)
    !
    ! Finds out the id number having its number
    !
    implicit none

    ! Imported parameters
    type(silam_source), intent(in) :: src
    character(len=*), intent(in) :: chSrcNm, chSectorNm

    ! Local variables
    integer :: i
    type(silam_sp) :: sp

    sp%sp => fu_work_string()

    sp%sp = fu_src_id(src%iSrcIdType, chSrcNm, chSectorNm) 

    do i = 1, fu_num_sources_total(src)
      if(src%src_info_lst(i)%chId == sp%sp)then
        fu_id_nbr_from_srcNms = i
        call free_work_array(sp%sp)
        return
      endif
    end do

    call msg_warning('No such source id :' + sp%sp,'fu_id_nbr_from_srcNms')
    fu_id_nbr_from_srcNms = int_missing
    call free_work_array(sp%sp)

  end function fu_id_nbr_from_srcNms


  !**********************************************************************

  integer function fu_compute_src_id_nbr(srcInfoLst, iElementToFill) 
    !
    ! Computes the next id number from the already-filled values
    ! Rules are:
    ! (1) The combination of source name plus source sector name IS UNIQUE
    ! (2) Depending on the requested merging of the emission, the IDs can be
    !     compiled in a different way using these two names. 
    ! (2.1) ID can be one for all sources thus mixing them all in one map
    ! (2.2) ID can be made from source name only
    ! (2.3) ID can be made from source sector name only
    ! (2.4) ID can be made from source name and source sector name combined. In this
    !       case all sources are computed totally separately
    ! (3) Each ID has its own IDnumber, which is then stored into the Lagrangian
    !     particles, maps, etc.
    !
    ! Procedure
    ! (1) Check for uniquness violation - and MERGE the sources with same names
    !     of source and sector
    ! (2) If the name combination is unique, check whether the ID is unique. 
    ! (2.1) If ID is unique - give it a new number
    ! (2.2) If ID is non-unique, use the already-given number for this ID.
    !
    implicit none

    ! Imported parameters
    type(silam_source_info), dimension(:), intent(in) :: srcInfoLst
    integer, intent(in) :: iElementToFill
    
    ! Local variables
    integer :: i, iMaxId

    !
    ! First check the combination of names. It must be unique. Merge the sources
    ! with the same name combinations
    !
    do i=1, iElementToFill-1
      if(srcInfoLst(i)%chSrcNm == srcInfoLst(iElementToFill)%chSrcNm .and. &
       & srcInfoLst(i)%chSectorNm == srcInfoLst(iElementToFill)%chSectorNm)then
        !
        ! Non-unique combination of names. Return the negated index of
        ! the source that has the same name and sector. The new source must be merged
        ! with it.
        !
        fu_compute_src_id_nbr = (-1) * i
        return
      endif
    end do

    !
    ! Name combination is indeed unique. Let's check the ID.
    !
    iMaxId = 0
    do i=1, iElementToFill-1
      if(srcInfoLst(i)%chId == srcInfoLst(iElementToFill)%chId)then
        fu_compute_src_id_nbr = srcInfoLst(i)%iIdNbr
        return
      endif
      iMaxId = max(iMaxId,srcInfoLst(i)%iIdNbr) ! find the last busy id nbr
    end do
    fu_compute_src_id_nbr = iMaxId + 1 ! No such IDs - take the first free number

  end function fu_compute_src_id_nbr


  !****************************************************************

  subroutine create_src_contain_grd_gen_src(em_source, &
                                          & grid_template, &
                                          & ifVerbose, ifMinimal, ifInventoryOnly)
    !
    ! This routine creates the grid similar to the given grid_template
    ! covering all sources of the em_source.
    ! The means the similarity is treated is dictated by the ifMinimal: if true, only
    ! projection is taken, all locations can be completely redefined.
    !
    implicit none

    ! Imported parameters
    type(silam_source), intent(in) :: em_source
    type(silja_grid), intent(inout) :: grid_template
    logical, intent(in) :: ifVerbose, ifMinimal, ifInventoryOnly

    ! Local variables
    integer :: iSrc, nx, ny
    logical :: ifExtended, ifFirstSrc
    type(silja_grid) :: gridTmp

    if(.not. defined(em_source))then
      call set_error('Undefiend source','create_src_cont_grd_gen_src')
      return
    endif
    ifFirstSrc = .true.
    gridTmp = grid_template
    !
    ! For each source we scan its central points and extending the template_grid 
    ! if these points are outside. For point and bomb sources it is only one point,
    ! while area source has plenty of them.
    !
    ! AREA
    !
    do iSrc = 1, em_source%n_area
      call create_src_contain_grd_a_src(em_source%a_ptr(iSrc)%a_src, &
                                      & grid_template, &
                                      & ifVerbose, ifMinimal .and. ifFirstSrc, ifExtended)
      if (error) then
        call msg('create_src_cont_grd_gen_src: Extended the grid for the area source:',iSrc)
        call report(em_source%a_ptr(iSrc)%a_src)
        call unset_error("create_src_contain_grd_gen_src")
        !FIXME Should unset error here or die
      endif                                
      if(ifVerbose .and. ifExtended)then
        call msg('create_src_cont_grd_gen_src: Extended the grid for the area source:',iSrc)
        call report(em_source%a_ptr(iSrc)%a_src)
      endif
      ifFirstSrc = .false.
    end do
    !
    ! NUCLEAR BOMB
    !
    do iSrc = 1, em_source%n_bomb
      call create_source_containing_grid(em_source%b_ptr(iSrc)%b_src, &
                                       & grid_template, &
                                       & ifVerbose, ifMinimal .and. ifFirstSrc, ifExtended)
      ifFirstSrc = .false.
    end do
    !
    ! WILD LAND FIRE
    !
    do iSrc = 1, em_source%n_fire
      call create_src_cont_grd_fire_src(em_source%fire_ptr(iSrc)%fire_src, &
                                      & grid_template, &
                                      & ifVerbose, ifExtended)
      ifFirstSrc = .false.
    end do
    !
    ! POINT
    !
    do iSrc = 1, em_source%n_point
      call create_src_cont_grd_p_src(em_source%p_ptr(iSrc)%p_src, &
                                   & grid_template, &
                                   & ifVerbose, ifMinimal .and. ifFirstSrc, ifExtended)
      if(ifVerbose .and. ifExtended)then
        call msg('create_src_cont_grd_gen_src: POINT source outside the grid:',iSrc)
        call report(em_source%p_ptr(iSrc)%p_src)
      endif   
      ifFirstSrc = .false.
    end do
    !
    ! VOLCANO
    !
    do iSrc = 1, em_source%n_volc
      call create_volc_src_cont_grd(em_source%volc_ptr(iSrc)%volc_src, &
                                   & grid_template, &
                                   & ifVerbose, ifMinimal .and. ifFirstSrc, ifExtended)
      if(ifVerbose .and. ifExtended)then
        call msg('create_src_cont_grd_gen_src: VOLCANO source outside the grid:',iSrc)
        call report(em_source%volc_ptr(iSrc)%volc_src)
      endif   
      ifFirstSrc = .false.
    end do
    !
    ! If we need only grid covering inventory, skip other source terms
    !
    if(.not. ifInventoryOnly)then
      !
      ! BIO-VOC
      !
      do iSrc = 1, em_source%n_bio_voc
        call create_source_containing_grid(em_source%bvoc_ptr(iSrc)%bvoc_src, grid_template, &
                                         & ifVerbose, ifMinimal .and. ifFirstSrc, ifExtended)
        if(ifVerbose .and. ifExtended)then
          call msg('create_src_cont_grd_gen_src: BVOC source outside the grid:',iSrc)
          call report(em_source%bvoc_ptr(iSrc)%bvoc_src)
        endif   
      end do
      !
      ! POLLEN
      !
      do iSrc = 1, em_source%n_pollen
        call create_source_containing_grid(em_source%pollen_ptr(iSrc)%pollen_src, grid_template, &
                                         & ifVerbose, ifMinimal .and. ifFirstSrc, ifExtended)
        if(ifVerbose .and. ifExtended)then
          call msg('create_src_cont_grd_gen_src: POLLEN source outside the grid:',iSrc)
          call report(em_source%pollen_ptr(iSrc)%pollen_src)
        endif   
      end do
      !
      ! SEA SALT
      !
      do iSrc = 1, em_source%n_sea_salt
        call create_source_containing_grid(em_source%sslt_ptr(iSrc)%sslt_src, grid_template, &
                                         & ifVerbose, ifMinimal .and. ifFirstSrc, ifExtended)
        if(ifVerbose .and. ifExtended)then
          call msg('create_src_cont_grd_gen_src: SEA SALT source outside the grid:',iSrc)
          call report(em_source%sslt_ptr(iSrc)%sslt_src)
        endif   
      end do
      !
      ! WIND BLOWN DUST
      !
      do iSrc = 1, em_source%n_wb_dust
        call create_source_containing_grid(em_source%wbdust_ptr(iSrc)%wbdust_src, grid_template, &
                                         & ifVerbose, ifMinimal .and. ifFirstSrc, ifExtended)
        if(ifVerbose .and. ifExtended)then
          call msg('create_src_cont_grd_gen_src: WIND BLOWN DUST source outside the grid:',iSrc)
          call report(em_source%wbdust_ptr(iSrc)%wbdust_src)
        endif   
      end do
      !
      ! DMS
      !
      do iSrc = 1, em_source%n_dms
        call create_source_containing_grid(em_source%dms_ptr(iSrc)%dms_src, grid_template, &
                                         & ifVerbose, ifMinimal .and. ifFirstSrc, ifExtended)
        if(ifVerbose .and. ifExtended)then
          call msg('create_src_cont_grd_gen_src: DMSsource outside the grid:',iSrc)
          call report(em_source%dms_ptr(iSrc)%dms_src)
        endif   
      end do
    endif  ! ifInventoryOnly
    !
    ! Finally, we received the grid, which is either a, possibly, increased grid_template or,
    ! if minimal grid requested, potentially somethig entirely different. Have to check for
    ! stupidity: the suggested grid is to be smaller than the template.
    !
    if(ifMinimal .and. defined(gridTmp))then
      call grid_dimensions(grid_template, nx,ny)
      iSrc = nx*ny
      call grid_dimensions(gridTmp, nx, ny)  ! stored original grid_template
      if(error)return
      if(iSrc > nx*ny)then
        call adjust_grid_dimensions(grid_template, nx, ny)
      endif
    endif  ! ifMinimal

  end subroutine create_src_contain_grd_gen_src


!  !************************************************************************
!  
!  subroutine create_src_containing_vert(source, emisVert, ifVerbose, ifFollowTemplate)
!    !
!    ! Creates a vertical, which covers the emission range. If match to template is forced,
!    ! it can just cut unnecessary layers. If not, the vertical is created as close to that
!    ! of the sources as possible
!    !
!    implicit none
!    
!    ! Imported parameters
!    type(silam_source), intent(in) :: source
!    type(silam_vertical), intent(inout) :: emisVert
!    logical, intent(in) :: ifVerbose, ifFollowTemplate
!    
!    ! Local variables
!    integer :: iLevStart, iLevEnd
!    real, dimension(:), pointer :: fIndices, fIndTmp
!    type(silam_vertical) :: vertTmp
!    
!    fIndices => fu_work_array()
!    fIndTmp => fu_work_array()
!    if(error)return
!    fIndices = real_missing
!    fIndTmp = real_missing
!    
!    if(ifFollowTemplate)then
!      if(defined(emisVert))then
!        vertTmp = emisVert
!      else
!        call set_error('Emission vertical is needed but undefined','create_src_containing_vert')
!        return
!      endif
!    else
!      call set_missing(vertTmp, .true.) ! will be used for lagrangian dynamics
!    endif
!    !
!    ! For Eulerian dynamics: vertical is better to correspond to dispersion_vertical
!    ! But it can have fewer layers
!    !
!    do iSrc = 1, em_source%n_area
!      call project_source_vertical(em_source%a_ptr(iSrc)%a_src, &
!                                 & vertTmp, ifFollowTemplate, fIndTmp, nLevels, &
!                                 & ifVerbose)
!      if(ifVerbose)call msg('Upper level of area source:'+fu_name(em_source%a_ptr(iSrc)%a_src), &
!                          & max(fIndTmp(1:nLevels)))
!      nLevels = fu_merge_real_arrays(fIndTmp, fIndices)
!    end do
!    !
!    ! NUCLEAR BOMB
!    !
!    do iSrc = 1, em_source%n_bomb
!      call project_source_vertical(em_source%b_ptr(iSrc)%b_src, &
!                                 & vertTmp, ifFollowTemplate, fIndTmp, nLevels, &
!                                 & ifVerbose)
!      if(ifVerbose)call msg('Upper level of area source:'+fu_name(em_source%a_ptr(iSrc)%a_src), &
!                          & max(fIndTmp(1:nLevels)))
!      nLevels = fu_merge_real_arrays(fIndTmp, fIndices)
!    end do
!    !
!    ! POINT
!    !
!    do iSrc = 1, em_source%n_point
!      call project_source_vertical(em_source%p_ptr(iSrc)%p_src, &
!                                 & vertTmp, ifFollowTemplate, fIndTmp, nLevels, &
!                                 & ifVerbose)
!      if(ifVerbose)call msg('Upper level of point source:'+fu_name(em_source%p_ptr(iSrc)%p_src), &
!                          & max(fIndTmp(1:nLevels)))
!      nLevels = fu_merge_real_arrays(fIndTmp, fIndices)
!    end do
!    !
!    ! BIO-VOC
!    !
!    do iSrc = 1, em_source%n_bio_voc
!      call project_source_vertical(em_source%bvoc_ptr(iSrc)%bvoc_src, &
!                                 & vertTmp, ifFollowTemplate, fIndTmp, nLevels, &
!                                 & ifVerbose)
!      if(ifVerbose)call msg('Upper level of bio-VOC source:' + &
!                          & fu_name(em_source%bvoc_ptr(iSrc)%bvoc_src), max(fIndTmp(1:nLevels)))
!      nLevels = fu_merge_real_arrays(fIndTmp, fIndices)
!    end do
!    !
!    ! POLLEN
!    !
!    do iSrc = 1, em_source%n_pollen
!      call project_source_vertical(em_source%pollen_ptr(iSrc)%pollen_src, &
!                                 & vertTmp, ifFollowTemplate, fIndTmp, nLevels, &
!                                 & ifVerbose)
!      if(ifVerbose)call msg('Upper level of pollen source:' + &
!                          & fu_name(em_source%pollen_ptr(iSrc)%pollen_src), max(fIndTmp(1:nLevels)))
!      nLevels = fu_merge_real_arrays(fIndTmp, fIndices)
!    end do
!    !
!    ! SEA SALT
!    !
!    do iSrc = 1, em_source%n_sea_salt
!      call project_source_vertical(em_source%sslt_ptr(iSrc)%sslt_src, &
!                                 & vertTmp, ifFollowTemplate, fIndTmp, nLevels, &
!                                 & ifVerbose)
!      if(ifVerbose)call msg('Upper level of sea salt source:' + &
!                          & fu_name(em_source%sslt_ptr(iSrc)%sslt_src), max(fIndTmp(1:nLevels)))
!      nLevels = fu_merge_real_arrays(fIndTmp, fIndices)
!    end do
!      
!    call free_work_array(indices)
!    call free_work_array(indTmp)
!    
!  end subroutine create_src_containing_vert


  !****************************************************************

  subroutine force_general_src_into_grid(em_source, grid_template, ifVerbose)
    !
    ! This routine creates the grid similar to the given grid_template
    ! covering all sources of the em_source
    !
    implicit none

    ! Imported parameters
    type(silam_source), intent(inout) :: em_source
    type(silja_grid), intent(in) :: grid_template
    logical, intent(in) :: ifVerbose

    ! Local variables
    integer :: iSrc
    logical :: ifCut
    type(silja_grid) :: gridTmp

    if(.not. defined(em_source))then
      call set_error('Undefiend source','force_general_src_into_grid')
      return
    endif

    gridTmp = grid_template

    !
    ! For each source we scan its central points and extending the template_grid 
    ! if these points are outside. For point and bomb sources it is only one point,
    ! while area source has plenty of them.
    !
    do iSrc = 1, em_source%n_area
      call force_source_into_grid(em_source%a_ptr(iSrc)%a_src, grid_template, ifVerbose, ifCut)
      if(ifVerbose .and. ifCut)then
        call msg('Cut the area source:'+fu_name(em_source%a_ptr(iSrc)%a_src),iSrc)
!        call report(em_source%a_ptr(iSrc)%a_src)
      endif   
    end do

    do iSrc = 1, em_source%n_bomb
      call create_source_containing_grid(em_source%b_ptr(iSrc)%b_src, gridTmp, &
                                       & ifVerbose, .false., ifCut)
      if(ifCut)then
        call msg('Bomb source is outside the grid:',iSrc)
        call report(grid_template)
        call report(em_source%b_ptr(iSrc)%b_src)
        call set_error('Bomb source is outside the grid:','force_general_src_into_grid')
        call unset_error('force_general_src_into_grid')
      endif
    end do

    do iSrc = 1, em_source%n_fire
      call force_source_into_grid(em_source%fire_ptr(iSrc)%fire_src, grid_template, ifVerbose, ifCut)
      if(ifVerbose .and. ifCut)then
        call msg('Cut the fire source:'+fu_name(em_source%fire_ptr(iSrc)%fire_src),iSrc)
!        call report(em_source%a_ptr(iSrc)%a_src)
      endif   
    end do

    do iSrc = 1, em_source%n_point
      call create_src_cont_grd_p_src(em_source%p_ptr(iSrc)%p_src, gridTmp, &
                                   & ifVerbose, .false., ifCut)
      if(ifCut)then
        call msg('Point source is outside the grid:',iSrc)
        call report(grid_template)
        call msg('It requires the grid:')
        call report(gridTmp)
        call report(em_source%p_ptr(iSrc)%p_src)
        call set_error('Point source is outside the grid:','force_general_src_into_grid')
        call unset_error('force_general_src_into_grid')
      endif
    end do

    do iSrc = 1, em_source%n_volc
      call create_volc_src_cont_grd(em_source%volc_ptr(iSrc)%volc_src, gridTmp, &
                                  & ifVerbose, .false., ifCut)
      if(ifCut)then
        call msg('Vocano source is outside the grid:',iSrc)
        call report(grid_template)
        call msg('It requires the grid:')
        call report(gridTmp)
        call report(em_source%volc_ptr(iSrc)%volc_src)
        call set_error('Volcano source is outside the grid:','force_general_src_into_grid')
        call unset_error('force_general_src_into_grid')
      endif
    end do

!    do iSrc = 1, em_source%n_pollen
!      call create_source_containing_grid(em_source%pollen_ptr(iSrc)%pollen_src, gridTmp, ifVerbose, ifCut)
!      if(ifCut)then
!        call msg('Cut the pollen source:'+fu_name(em_source%pollen_ptr(iSrc)%pollen_src),iSrc)
!      endif
!    end do
!
!    do iSrc = 1, em_source%n_sea_salt
!      call create_source_containing_grid(em_source%sslt_ptr(iSrc)%sslt_src, gridTmp, ifVerbose, ifCut)
!      if(ifCut)then
!        call msg('Cut the sea salt source:'+fu_name(em_source%sslt_ptr(iSrc)%sslt_src),iSrc)
!      endif
!    end do

  end subroutine force_general_src_into_grid


  !********************************************************************************

  subroutine store_as_namelist_general_src(source, chNameTemplate, ifOriginalGrid)
    !
    ! Records the given source as namelist 
    !
    implicit none

    ! Imported parameters
    TYPE(silam_source), INTENT(in) :: source
    character(len=*), intent(in) :: chNameTemplate
    logical, intent(in) :: ifOriginalGrid

    ! Local variables
    integer :: iSrc, uOut

    if(.not. defined(source))then
     call set_error('Undefined source given','store_as_namelist_general_src')
     return
    endif

    uOut = fu_next_free_unit()
    if(error)return
    !
    ! Area sources
    !
    if(source%n_area > 0)then
      open(unit=uOut, file=chNameTemplate+'.sa3', iostat=iSrc)
      if(iSrc /= 0)then
        call set_error('Cannot ope output file:'+chNameTemplate+'_emis.sa3','store_as_namelist_general_src')
        return
      endif

      do iSrc = 1, source%n_area
        call store_as_namelist(source%a_ptr(iSrc)%a_src, uOut, ifOriginalGrid)
        if(error)return
      end do
      close(uOut)
    endif
    !
    ! Point sources
    !
    if(source%n_point > 0)then
      open(unit=uOut, file=chNameTemplate+'.sp5', iostat=iSrc)
      if(iSrc /= 0)then
        call set_error('Cannot ope output file:'+chNameTemplate+'_emis.sp5','store_as_namelist_general_src')
        return
      endif

      do iSrc = 1, source%n_point
        call store_as_namelist(source%p_ptr(iSrc)%p_src, uOut, ifOriginalGrid)
        if(error)return
      end do
      close(uOut)
    endif

    !
    ! Volcano sources
    !
    if(source%n_volc > 0)then
      open(unit=uOut, file=chNameTemplate+'.sv1', iostat=iSrc)
      if(iSrc /= 0)then
        call set_error('Cannot ope output file:'+chNameTemplate+'.sv1','store_as_namelist_general_src')
        return
      endif

      do iSrc = 1, source%n_volc
        call store_volc_src_as_namelist(source%volc_ptr(iSrc)%volc_src, uOut, ifOriginalGrid, .true.)
        if(error)return
      end do
      close(uOut)
    endif

  end subroutine store_as_namelist_general_src


  !************************************************************************************
  
  logical function fu_if_dump_emission_flux(emSrc)
    !
    ! Just checks whether emission dump into time=resolved areas source is needed
    !
    implicit none

    ! Imported parameters
    TYPE(silam_source), INTENT(inout) :: emSrc

    fu_if_dump_emission_flux = emSrc%ifDumpFlux

  end function fu_if_dump_emission_flux


  !************************************************************************************
  
  subroutine dump_emission_mass_map(emSrc, now, chCaseNm, outputFNmTemplate, &
                                  & mapEmis, mapEmis_Px, mapEmis_Py, mapEmis_Pz, &
                                  & timestep)
    !
    ! Sends the full content of emission mass map to the provided GrADS file.
    ! Opens the file, crates grads structure, if needed. Should the indFile index be given
    ! it is used directly, just checked whether the binary should be switched.
    !
    implicit none

    ! Imported parameters
    TYPE(silam_source), INTENT(inout) :: emSrc
    type(silja_time), intent(in) :: now
    character(len=*), intent(in) :: chCaseNm
    type(grads_template), intent(in) :: outputFNmTemplate
    type(TMass_map), target, intent(in) :: mapEmis, mapEmis_Px, mapEmis_Py, mapEmis_Pz
    type(silja_interval), intent(in) :: timestep

    ! Local variables
    integer :: iFileManipulation
    logical, save :: ifFirstTime=.true.
    character(len=fnlen) :: grads_FNm
    type(Tmass_map), pointer :: pMapEmis, pMapPx, pMapPy, pMapPz

    !
    ! Opening files, if requested by the time series allocation.
    ! Number of files to open depends on the number of sources processed
    ! simultaneously
    !
    iFileManipulation = DoNothing
    if(ifFirstTime) then  ! First-time output
      iFileManipulation = OpenFiles
      if(.not. defined(emSrc%FirstDumpTime))then
        emSrc%FirstDumpTime = now
        emSrc%PrevDumpTime = now
        emSrc%NextDumpTime = now
      endif
      ifFirstTime = .false.
    else
      select case(emSrc%DumpFilesArrangement)
        case(all_in_one)
          iFileManipulation = DoNothing
        case(hourly_new_file)
          if(fu_hour(now) /= fu_hour(emSrc%prevDumpTime)) &
                                                & iFileManipulation = SwitchBinary
        case(daily_new_file)
          if(fu_day(now) /= fu_day(emSrc%prevDumpTime)) &
                                                & iFileManipulation = SwitchBinary
        case(monthly_new_file)
          if(fu_mon(now) /= fu_mon(emSrc%prevDumpTime)) &
                                                & iFileManipulation = SwitchBinary
        case(yearly_new_file)
          if(fu_year(now) /= fu_year(emSrc%prevDumpTime)) &
                                                & iFileManipulation = SwitchBinary
        case default
          call set_error('Unknown time series arrangement','dump_emission_mass_map')
          return
      end select
    endif
    if(error)return

    if(.not. (now == emSrc%NextDumpTime))return

    !
    ! Now handle the files
    !
    select case (iFileManipulation)

      case (DoNothing)

!call msg('Do nothing')

      case (OpenFiles) ! Totally new files, including ctls
        !
        ! Close opened files first
        ! 
        if(emSrc%indexDumpFile /= int_missing)then ! Close existing file
          call close_gradsfile_o(emSrc%indexDumpFile,"")
        endif
        !
        ! Get file name from tempaltes
        !
        grads_fNm = fu_FNm(outputFNmTemplate, &
                         & emSrc%FirstDumpTime, &  ! ini_time
                         & emSrc%FirstDumpTime, &  ! anal_time - here it is the same
                         & (now - emSrc%FirstDumpTime), &  ! forecast length
                         & chCaseNm, '') + '_ems_dump.grads'
        if(error)return

!call msg('Open file:' + grads_FNm)
        !
        ! For the multi-file output we should store the template rather than a file
        ! name into the GrADS structure. 
        !
        if(emSrc%DumpFilesArrangement == all_in_one)then
          
          emSrc%indexDumpFile = open_gradsfile_o('', grads_fNm, mapEmis%gridTemplate)
        else

          emSrc%indexDumpFile = open_gradsfile_o('', grads_fNm, mapEmis%gridTemplate, &
                                               & fu_pure_grads_template(outputFNmTemplate, &
                                                                      & chCaseNm, '') + &
                                                                      & '_ems_dump.grads')
        endif

      case(SwitchBinary)
        !
        ! Get new file name from tempaltes
        !
        grads_fNm = fu_FNm(outputFNmTemplate, &
                         & emSrc%FirstDumpTime, &  ! ini_time
                         & emSrc%FirstDumpTime, &  ! anal_time - here it is the same
                         & (now - emSrc%FirstDumpTime), &  ! forecast length
                         & chCaseNm, '') + '_ems_dump.grads'
        if(error)return
        !
        ! Close existing GrADS binary and open a new one - in one routine
        ! ctl file must not be written or closed. May be, some rewinding of the 
        ! structure indices...
        !
!call msg('Switch binary to:'+grads_FNm)
        call switch_grads_binary_o(emSrc%indexDumpFile, '', grads_fNm, 0, now)
        
      case default
        call msg('Unknown file manipulation:', iFileManipulation)
        call set_error('Unknown file manipulation ','dump_emission_mass_map')
    end select
    if(error)return
    !
    ! Now the file is ready, dump the available mass maps to the grads file.
    ! Do not forget: emission mass map has emission masses cumulated during the timestep.
    !
    pMapEmis => mapEmis
call msg('Dumping emission')
    call mass_map_to_grads_file(pMapEmis, int_missing, '', now, &
                              & 1./fu_sec(timestep), emSrc%indexDumpFile)
    if(error)return
    !
    ! Moments may or may not be sent to the file
    !
    if(emSrc%ifDumpMoment)then

      pMapPx => mapEmis_Px
call msg('Dumping Px')
      call mass_map_to_grads_file(pMapPx, int_missing, '', now, &
                                & 1./fu_sec(timestep), emSrc%indexDumpFile)
      if(error)return

      pMapPy => mapEmis_Py
call msg('Dumping Px')
      call mass_map_to_grads_file(pMapPy, int_missing, '', now, &
                                & 1./fu_sec(timestep), emSrc%indexDumpFile)
      if(error)return

      pMapPz => mapEmis_Pz
call msg('Dumping Px')
      call mass_map_to_grads_file(pMapPz, int_missing, '', now, &
                                & 1./fu_sec(timestep), emSrc%indexDumpFile)
      if(error)return
    endif
    !
    ! Final step - close all or prepare target IDs and fields in temporary stack for further
    ! output. Actions with individual fields depend on their types
    !
    if(defined(emSrc%LastDumpTime))then
      if(.not. fu_between_times(now, emSrc%prevDumpTime, emSrc%LastDumpTime, .true.))then
        call close_emission_dump(emSrc, chCaseNm, outputFNmTemplate, pMapEmis, pMapPx, pMapPy, pMapPz)
        emSrc%ifDumpFlux = .false.
      endif
    else
      call msg('Next emission dump time: ' + &
                              & fu_str(emSrc%prevDumpTime + emSrc%DumpTimestep))
      emSrc%prevDumpTime = now
      emSrc%nextDumpTime = now + emSrc%DumpTimestep
    endif

  end subroutine dump_emission_mass_map


  !*****************************************************************************
  
  subroutine close_emission_dump(emSrc, chCaseNm, outputFNmTemplate, &
                               & pMapEmis, pMapPx, pMapPy, pMapPz)
    !
    ! Close the grads file (will create the ctl and super-ctl)
    ! and write down the header of the source term
    !
    implicit none

    TYPE(silam_source), INTENT(in) :: emSrc
    character(len=*), intent(in) :: chCaseNm
    type(grads_template), intent(in) :: outputFNmTemplate
    type(TMass_map), pointer :: pMapEmis, pMapPx, pMapPy, pMapPz

    if(.not. emSrc%ifDumpFlux)return  ! instead of creating fu_if_dump, just return from here
    if (emSrc%indexDumpFile == int_missing) then
      ! if dump was requested but no forecast was run, this could happen if data assimilation 
      ! only runs the analysis iteration.
      call msg_warning('Emission dump file not defined', 'close_emission_dump')
      return
    end if
    call close_gradsfile_o(emSrc%indexDumpFile,"")

    call write_area_src_from_mass_map(pMapEmis, emSrc%FirstDumpTime, emSrc%PrevDumpTime, &
                                    & chCaseNm, '_ems_dump.grads', &
                                    & pMapPx, pMapPy, pMapPz, outputFNmTemplate)
  end subroutine close_emission_dump


  !=========================================================================  
  !=========================================================================  
  !
  !        Private functions for general source
  !
  !=========================================================================  
  !=========================================================================  

  !*****************************************************************

  LOGICAL FUNCTION fu_source_defined(source)

    ! Description:
    ! Returns a true value, if the point_source has been given a value
    ! using correct tools, and no errors occurred.
    !
    IMPLICIT NONE
    
    ! Imported parameters with intent IN:
    TYPE(silam_source), INTENT(in) :: source

    fu_source_defined = fu_true(source%defined)

  END FUNCTION fu_source_defined


  !*****************************************************************

  integer function fu_NbrOf_sources_total(source, iType) result(Nbr)
    !
    ! Returns the number of the sources of the given type or
    ! the total number of sources if iType is omitted
    ! Note: this function does not care on source name-sector grouping,
    ! it returns the total number of separately defined aliases
    !
    ! An update: the function uses the n_point, n_area, n_bomb params that
    ! are supposed to be set during the file reading. Do NOT use this function
    ! before or during that process.
    !
    implicit none

    ! Imported parameters
    type(silam_source), intent(in) :: source
    integer, intent(in), optional :: iType

    if(present(iType))then  ! Compute concrete type
      select case(iType)
        case(area_source)
          Nbr = source%n_area
        case(bio_voc_source)
          Nbr = source%n_bio_voc
        case(bomb_source)
          Nbr = source%n_bomb
        case(fire_source)
          Nbr = source%n_fire
        case(point_source)
          Nbr = source%n_point
        case(pollen_source)
          Nbr = source%n_pollen
        case(sea_salt_source)
          Nbr = source%n_sea_salt
        case(wind_blown_dust_source)
          Nbr = source%n_wb_dust
        case(dms_source)
          nbr = source%n_dms
        case(volcano_source)
          nbr = source%n_volc
        case default 
          call msg('Strange source type:',iType)
          call set_error('Strange source type','fu_NbrOf_sources_total')
          Nbr=0
      end select
    else   ! No iType - sum up all the sources
      Nbr = fu_num_sources_total(source)
    endif  ! if iType present

  end function fu_NbrOf_sources_total


  !******************************************************************

  integer function fu_NbrOf_source_ids(source) result(nbr)
    !
    ! Returns the total number of different source ids 
    ! The value is found via scanning the whole src_info_lst vector
    ! Note that if we use time-zone (or country-) wise split, mapping
    ! is computed specifically and the source indices are ignored
    !
    implicit none

    ! Imported parameters
    type(silam_source), intent(in) :: source

    ! Local variables
    integer :: iTmp

    if(source%ifUseSourceIdMapping)then
      Nbr = maxval(source%arSourceIdMapping)
    else
      Nbr = 0
      do iTmp=1, fu_num_sources_total(source)
        Nbr = max(Nbr, source%src_info_lst(iTmp)%iIdNbr)
      end do
    endif
  end function fu_NbrOf_source_ids


  !*******************************************************************

  function fu_get_id_str_from_id_nbr(src,IdNbr) result(chId)
    !
    ! Find the character Id from its number
    !
    implicit none

    ! Return value of the function
    character(len=clen) :: chId

    ! Imported parameters
    type(silam_source), intent(in) :: src
    integer, intent(in) :: IdNbr

    ! Local variables
    integer :: i

    if(src%ifUseSourceIdMapping)then
      chId = src%chSplitNames(IdNbr)
      return
    else
      do i = 1, fu_num_sources_total(src)
        if(src%src_info_lst(i)%iIdNbr == IdNbr)then
          chId = src%src_info_lst(i)%chId
          return
        endif
      end do
    endif

    call msg('Failed to find the Id number:',IdNbr)
    call set_error('Failed to find the Id number','fu_get_id_str_from_id_nbr')

  end function fu_get_id_str_from_id_nbr


  ! ***************************************************************
  ! Function fu_source_name_all_sources has been removed due to errors and lack of use.


  FUNCTION fu_source_name(src, indexSrc)
    
    ! Returns the concatenated names of all sources included into the structure
 
    IMPLICIT NONE
    !
    ! The return value of this function:
    CHARACTER (LEN = clen) :: fu_source_name
    integer, intent(in) :: indexSrc
    !
    ! Imported parameters with intent(in):
    TYPE(silam_source), INTENT(in) :: src

    ! Local variables
    integer :: i

    IF (.NOT.defined(src)) THEN
      CALL set_error('undefined source given','fu_source_name')
      RETURN
    END IF

    do i=1,src%n_area
      if(fu_source_nbr(src%a_ptr(i)%a_src) == indexSrc) then
        fu_source_name = fu_name(src%a_ptr(i)%a_src)
        return
      endif
    end do
    do i=1,src%n_bio_voc
      if(fu_source_nbr(src%bvoc_ptr(i)%bvoc_src) == indexSrc)then
        fu_source_name = fu_name(src%bvoc_ptr(i)%bvoc_src)
        return
      endif
    end do
    do i=1,src%n_bomb
      if(fu_source_nbr(src%b_ptr(i)%b_src) == indexSrc)then
        fu_source_name = fu_name(src%b_ptr(i)%b_src)
        return
      endif
    end do
    do i=1,src%n_fire
      if(fu_source_nbr(src%fire_ptr(i)%fire_src) == indexSrc)then
        fu_source_name = fu_name(src%fire_ptr(i)%fire_src)
        return
      endif
    end do
    do i=1,src%n_point
      if(fu_source_nbr(src%p_ptr(i)%p_src) == indexSrc) then
        fu_source_name = fu_name(src%p_ptr(i)%p_src)
        return
      endif
    end do
    do i=1,src%n_pollen
      if(fu_source_nbr(src%pollen_ptr(i)%pollen_src) == indexSrc)then
        fu_source_name = fu_name(src%pollen_ptr(i)%pollen_src)
        return
      endif
    end do
    do i=1,src%n_sea_salt
      if(fu_source_nbr(src%sslt_ptr(i)%sslt_src) == indexSrc)then
        fu_source_name = fu_name(src%sslt_ptr(i)%sslt_src)
        return
      endif
    end do
    do i=1,src%n_wb_dust
      if(fu_source_nbr(src%wbdust_ptr(i)%wbdust_src) == indexSrc)then
        fu_source_name = fu_name(src%wbdust_ptr(i)%wbdust_src)
        return
      endif
    end do
    do i=1,src%n_dms
      if(fu_source_nbr(src%dms_ptr(i)%dms_src) == indexSrc)then
        fu_source_name = fu_name(src%dms_ptr(i)%dms_src)
        return
      endif
    end do
    do i=1,src%n_volc
      if(fu_source_nbr(src%volc_ptr(i)%volc_src) == indexSrc)then
        fu_source_name = fu_name(src%volc_ptr(i)%volc_src)
        return
      endif
    end do
    call msg('Source index given:', indexSrc)
    call set_error('index not found in the source list','fu_soruce_name')
    fu_source_name = ''

  END FUNCTION fu_source_name


  ! ***************************************************************


  FUNCTION fu_earliest_start_time(src)
    
    ! Returns the earliest start time of the release of the sources
 
    IMPLICIT NONE
    !
    ! The return value of this function:
    TYPE(silja_time) :: fu_earliest_start_time
    !
    ! Imported parameters with intent(in):
    TYPE(silam_source), INTENT(in) :: src

    ! Local variables
    integer :: i

    IF (.NOT.defined(src)) THEN
      CALL set_error('undefined source given','fu_earliest_start_time')
      RETURN
    END IF

    !
    ! Sources like pollen and sea salt emit always
    !
    if(src%n_bio_voc + src%n_pollen + src%n_sea_salt + src%n_wb_dust + src%n_dms > 0)then
      fu_earliest_start_time = really_far_in_past
      return
    endif

    !
    ! Take something as a starting point
    !
    fu_earliest_start_time = really_far_in_future
    !
    ! Check if there are other times earlier than the first one
    !
    do i=1,src%n_point
      if(.not. defined(src%p_ptr(i)%p_src))cycle
      if(fu_start_time(src%p_ptr(i)%p_src) < fu_earliest_start_time) &
              & fu_earliest_start_time = fu_start_time(src%p_ptr(i)%p_src)
    end do
    do i=1,src%n_area
      if(.not. defined(src%a_ptr(i)%a_src))cycle
      if(fu_start_time(src%a_ptr(i)%a_src) < fu_earliest_start_time) &
              & fu_earliest_start_time = fu_start_time(src%a_ptr(i)%a_src)
    end do
    do i=1,src%n_bomb
      if(.not. defined(src%b_ptr(i)%b_src))cycle
      if(fu_time(fu_start_position(src%b_ptr(i)%b_src)) < fu_earliest_start_time) &
              & fu_earliest_start_time = fu_time(fu_start_position(src%b_ptr(i)%b_src))
    end do
    do i=1,src%n_fire
      if(.not. defined(src%fire_ptr(i)%fire_src))cycle
      if(fu_start_time(src%fire_ptr(i)%fire_src) < fu_earliest_start_time) &
              & fu_earliest_start_time = fu_start_time(src%fire_ptr(i)%fire_src)
    end do
    do i=1,src%n_volc
      if(.not. defined(src%volc_ptr(i)%volc_src))cycle
      if(fu_start_time(src%volc_ptr(i)%volc_src) < fu_earliest_start_time) &
              & fu_earliest_start_time = fu_start_time(src%volc_ptr(i)%volc_src)
    end do

  END FUNCTION fu_earliest_start_time


  ! ***************************************************************


  FUNCTION fu_latest_end_time(src)
    
    ! Returns the earliest start time of the release of the sources
 
    IMPLICIT NONE
    !
    ! The return value of this function:
    TYPE(silja_time) :: fu_latest_end_time
    !
    ! Imported parameters with intent(in):
    TYPE(silam_source), INTENT(in) :: src

    ! Local variables
    integer :: i

    IF (.NOT.defined(src)) THEN
      CALL set_error('undefined source given','fu_latest_end_time')
      RETURN
    END IF
    !
    ! Sources like pollen and sea salt emit always
    !
    if(src%n_bio_voc + src%n_pollen + src%n_sea_salt + src%n_wb_dust + src%n_dms > 0)then
      fu_latest_end_time = really_far_in_future
      return
    else
      fu_latest_end_time = really_far_in_past
    endif
    !
    ! If not, inventory sources have to be taken into account
    !
    do i=1,src%n_area
      if(fu_latest_end_time < fu_end_time(src%a_ptr(i)%a_src)) &
                                 & fu_latest_end_time = fu_end_time(src%a_ptr(i)%a_src)
    end do
    do i=1,src%n_bomb
      if(fu_latest_end_time < fu_time(fu_start_position(src%b_ptr(i)%b_src))) &
                                 & fu_latest_end_time = fu_time(fu_start_position(src%b_ptr(i)%b_src))
    end do
    do i=1,src%n_point
      if(fu_latest_end_time < fu_end_time(src%p_ptr(i)%p_src)) &
                                 & fu_latest_end_time = fu_end_time(src%p_ptr(i)%p_src)
    end do
    do i=1,src%n_fire
      if(fu_latest_end_time < fu_end_time(src%fire_ptr(i)%fire_src)) &
                                 & fu_latest_end_time = fu_end_time(src%fire_ptr(i)%fire_src)
    end do
    do i=1,src%n_volc
      if(fu_latest_end_time < fu_end_time(src%volc_ptr(i)%volc_src)) &
                                 & fu_latest_end_time = fu_end_time(src%volc_ptr(i)%volc_src)
    end do
    if(fu_latest_end_time == really_far_in_past)then
      call set_error('No known sources in the list','fu_latest_end_time')
      return
    endif

  END FUNCTION fu_latest_end_time


  ! ***************************************************************


  FUNCTION fu_source_start_time(src, iType, index)
    
    ! Description:
    ! Returns the end time of the release of point source
 
    IMPLICIT NONE
    !
    ! The return value of this function:
    TYPE(silja_time) :: fu_source_start_time
    !
    ! Imported parameters with intent(in):
    TYPE(silam_source), INTENT(in) :: src
    integer, intent(in) :: iType, index

    IF (.NOT.defined(src)) THEN
      CALL set_error('undefined source given','fu_source_end_time')
      RETURN
    END IF

    SELECT CASE (iType)
      CASE(area_source)
        fu_source_start_time = fu_start_time(src%a_ptr(index)%a_src)
      CASE(bomb_source)
        fu_source_start_time = fu_time(fu_start_position(src%b_ptr(index)%b_src))
      CASE(fire_source)
        fu_source_start_time = fu_start_time(src%fire_ptr(index)%fire_src)
      CASE(point_source)
        fu_source_start_time = fu_start_time(src%p_ptr(index)%p_src)
      CASE(volcano_source)
        fu_source_start_time = fu_start_time(src%volc_ptr(index)%volc_src)
      CASE(bio_voc_source, pollen_source, sea_salt_source, wind_blown_dust_source, dms_source)
        fu_source_start_time = really_far_in_past
      CASE DEFAULT
        CALL set_error('Unknown source type','fu_source_end_time')
    END SELECT

  END FUNCTION fu_source_start_time


  ! ***************************************************************


  FUNCTION fu_source_end_time(src, iType, index)
    
    ! Description:
    ! Returns the end time of the release of point source
 
    IMPLICIT NONE
    !
    ! The return value of this function:
    TYPE(silja_time) :: fu_source_end_time
    !
    ! Imported parameters with intent(in):
    TYPE(silam_source), INTENT(in) :: src
    integer, intent(in) :: iType, index

    IF (.NOT.defined(src)) THEN
      CALL set_error('undefined source given','fu_source_end_time')
      RETURN
    END IF

    SELECT CASE (iType)
      CASE(area_source)
        fu_source_end_time = fu_end_time(src%a_ptr(index)%a_src)
      CASE(bomb_source)
        fu_source_end_time = fu_time(fu_start_position(src%b_ptr(index)%b_src))
      CASE(fire_source)
        fu_source_end_time = fu_end_time(src%fire_ptr(index)%fire_src)
      CASE(point_source)
        fu_source_end_time = fu_end_time(src%p_ptr(index)%p_src)
      CASE(volcano_source)
        fu_source_end_time = fu_end_time(src%volc_ptr(index)%volc_src)
      CASE(bio_voc_source, pollen_source, sea_salt_source, wind_blown_dust_source, dms_source)
        fu_source_end_time = really_far_in_future
      CASE DEFAULT
        CALL set_error('Unknown source type','fu_source_end_time')
    END SELECT

  END FUNCTION fu_source_end_time



  ! ***************************************************************


  FUNCTION fu_source_duration(src, iType, srcIndex)
    
    ! Description:
    ! Returns the duration of release of point source.
 
    IMPLICIT NONE
    !
    ! The return value of this function:
    TYPE(silja_interval) :: fu_source_duration

    ! Imported parameters with intent(in):
    TYPE(silam_source), INTENT(in) :: src
    integer, intent(in) :: iType, srcIndex

    IF (.NOT.defined(src)) THEN
      CALL set_error('undefined source given','fu_source_duration')
      RETURN
    END IF

    SELECT CASE (iType)
    CASE(point_source)
      fu_source_duration = fu_duration(src%p_ptr(srcIndex)%p_src)
    CASE(area_source)
      fu_source_duration = fu_duration(src%a_ptr(srcIndex)%a_src)
    CASE(fire_source)
      fu_source_duration = fu_duration(src%fire_ptr(srcIndex)%fire_src)
    CASE(bomb_source)
      fu_source_duration = interval_missing
    CASE(volcano_source)
      fu_source_duration = fu_duration(src%volc_ptr(srcIndex)%volc_src)
    CASE(bio_voc_source, pollen_source, sea_salt_source, dms_source)
      fu_source_duration = very_long_interval
    CASE DEFAULT
      CALL set_error('Unknown source type','fu_source_duration')
    END SELECT

  END FUNCTION fu_source_duration


!  !*****************************************************************
!
!
!   Lagrangian environment needs substantial rethinking
!   Lagrangian environment needs substantial rethinking
!   Lagrangian environment needs substantial rethinking
!   Lagrangian environment needs substantial rethinking
!   Lagrangian environment needs substantial rethinking
!
!
!  real function fu_source_rate_local_unit(src, iType, index, iTimeSlot) result(rate)
!    !
!    ! Returns the release rate of the specific source of iType
!    ! pointed by the index, at the time moment pointed by iTimeSlot
!    ! Should the time variation indices are to be used - scales the
!    ! rate with them
!    !
!    ! Time variation coefs are NOT HERE. Rate_local_unit is just for 
!    ! the overall slot rate, not for time fluctuations. Slot can be long
!    !
!    ! ATTENTION. Units are NOT SI 
!    !
!    implicit none
!
!    type(silam_source), intent(in) :: src
!    integer, intent(in) :: iType, iIndex, iTimeSlot
!
!    rate = 0.
!    IF (.NOT.defined(src)) THEN
!      CALL set_error('undefined source given','fu_source_rate_local_unit')
!      RETURN
!    END IF
!
!    SELECT CASE (iType)
!    CASE(point_source)
!      if(index > src%n_point)then
!        call set_error('Too large point source index','fu_source_rate_local_unit')
!        return
!      endif
!      rate = fu_rate_local_unit(src%p_ptr(index)%p_src, iTimeSlot)
!
!    CASE(area_source)
!      if(index > src%n_area)then
!        call set_error('Too large area source index','fu_source_rate_local_unit')
!        return
!      endif
!      rate = fu_rate_local_unit(src%a_ptr(index)%a_src, iTimeSlot)
!
!    CASE(bomb_source) ! No time slot is needed
!      rate = fu_yield(src%b_ptr(index)%b_src)
!
!    CASE DEFAULT
!      CALL set_error('Unknown source type','fu_source_rate_local_unit')
!    END SELECT
!
!  end function fu_source_rate_local_unit


  !***************************************************************
!   Lagrangian environment needs substantial rethinking
!   Lagrangian environment needs substantial rethinking
!   Lagrangian environment needs substantial rethinking
!   Lagrangian environment needs substantial rethinking
!   Lagrangian environment needs substantial rethinking
!   Lagrangian environment needs substantial rethinking
!   Lagrangian environment needs substantial rethinking
!   Lagrangian environment needs substantial rethinking

!  function fu_unit_of_source_rate(src, i_point, i_area, i_bomb) result(chUnit)
!    !
!    ! Returns the unit of the rate of the given source - either
!    ! point or area or bomb
!    !
!    implicit none
!
!    ! Return value
!    character(len=unitNmLen) :: chUnit
!
!    ! Imported parameters
!    type(silam_source), intent(in) :: src
!    integer, intent(in), optional :: i_point, i_area, i_bomb
!
!    if(present(i_point))then
!      if(i_point > 0 .and. i_point <= src%n_point)then
!        if(present(i_bomb))then
!          if(i_bomb > 0 .and.i_bomb <= src%n_bomb)then
!            call set_error('Both point and bomb source indices are reasonable', &
!                         & 'fu_unit_of_source_rate')
!            return
!          endif
!        endif
!        if(present(i_area))then
!          if(i_area > 0 .and. i_area <= src%n_area)then
!            call set_error('Both point and area source indices are reasonable', &
!                         & 'fu_unit_of_source_rate')
!            return
!          endif
!        endif
!        chUnit = fu_unit_of_rate(src%p_ptr(i_point)%p_src)
!        return
!      endif  ! Point index reasonable
!    endif ! Point source present
!
!    if(present(i_bomb))then
!      if(i_bomb > 0 .and.i_bomb <= src%n_bomb)then
!        if(present(i_area))then
!          if(i_area > 0 .and. i_area <= src%n_area)then
!            call set_error('Both bomb and area source indices are reasonable', &
!                         & 'fu_unit_of_source_rate')
!            return
!          endif
!        endif
!        chUnit = fu_unit_of_rate(src%b_ptr(i_bomb)%b_src)
!        return
!      endif
!    endif
!
!    if(present(i_area))then
!      if(i_area > 0 .and.i_area <= src%n_area)then
!        chUnit = fu_unit_of_rate(src%a_ptr(i_area)%a_src)
!      else
!        call set_error('Strange area source index','fu_unit_of_source_rate')
!      endif
!    else
!      call set_error('All point, area and bomb source indices are irrelevant', &
!                   & 'fu_unit_of_source_rate')
!    endif
!
!  end function fu_unit_of_source_rate


!   Lagrangian environment needs substantial rethinking
!   Lagrangian environment needs substantial rethinking
!   Lagrangian environment needs substantial rethinking
!   Lagrangian environment needs substantial rethinking
!   Lagrangian environment needs substantial rethinking
!   Lagrangian environment needs substantial rethinking
!   Lagrangian environment needs substantial rethinking
!   Lagrangian environment needs substantial rethinking
!  !*******************************************************************
!
!  subroutine set_factor_to_basic_unit_src(source)
!    !
!    ! Sets factors to the basic units for all sources in the pool
!    !
!    implicit none
!
!    ! Imported parameters
!    type(silam_source), intent(inout) :: source
!
!    ! Local variables
!    integer :: i
!
!    do i=1,source%n_point
!      call set_factor_to_basic_unit(source%p_ptr(i)%p_src)
!    end do
!    do i=1,source%n_area
!      call set_factor_to_basic_unit(source%a_ptr(i)%a_src)
!    end do
!    do i=1,source%n_bomb
!      call set_factor_to_basic_unit(source%b_ptr(i)%b_src)
!    end do
!  end subroutine set_factor_to_basic_unit_src


  !*****************************************************************

  subroutine amounts_from_src_species_unit(src, &
                                         & transport_species, nSpecies,  &
                                         & amounts, &
                                         & start, duration, &
                                         & emission_species, refEmis2Transp_mass)
    !
    ! Returns the total amount of the emission from the inventory-type sources
    ! - of the transported species.
    !
    implicit none

    ! Imported parameters with intent IN
    type(silam_source), intent(in) :: src
    type(silja_time), intent(in), optional :: start
    type(silja_interval), intent(in), optional :: duration
    type(silam_species), dimension(:), pointer :: transport_species, emission_species
    integer, intent(in) :: nSpecies
    type(TspeciesReference), dimension(:), pointer :: refEmis2Transp_mass

    ! This needs to be allocated and sufficiently big.
    real(r8k), dimension(:), intent(out) :: amounts

    ! Local variables
    real, dimension(:), pointer :: fWork
    real, dimension(:), pointer :: amounts_single_src 
    real, dimension(:,:), pointer :: amounts_tmp
    type(silam_species), dimension(:), pointer :: species_single_src
    type(silam_species_arr_ptr), dimension(:), allocatable :: ar_species_single_src
    integer :: nspecies_single_src, i, j, k, isp, iTmp, iThread, nThreads
    logical :: ifTimePresent
    type(chemical_adaptor) :: adaptor
    real(r8k), dimension(:,:), allocatable :: threadamounts
    character(len=*), parameter :: sub_name = 'amounts_from_src_species_unit'
    
    if (present(start) .and. .not. present(duration)) then
      call set_error('If present, must define both start and duration', sub_name)
      return
    end if

    ifTimePresent = present(start) .and. present(duration)

    amounts(1:nSpecies) = 0.0


    
!Glob still fails with parallell, so disable it..
!!    !$OMP PARALLEL IF (.False.)  DEFAULT (NONE) &
!!    !$OMP & PRIVATE(amounts_single_src,  adaptor, i,j,k,isp,iTmp, species_single_src, &
!!    !$OMP & nSpecies_single_src, iThread, nthreads) &
!!    !$OMP & SHARED(src, transport_species, nSpecies, amounts, start, duration, emission_species, threadamounts, &
!!    !$OMP & refEmis2Transp_mass, ifTimePresent, error, amounts_tmp, ar_species_single_src, fWork)

    iThread = 0
!!    !$ iThread = OMP_GET_THREAD_NUM()


!!    !$OMP MASTER
    nthreads = 1
!!    !$ nthreads = omp_get_num_threads()

!!    !$OMP MASTER
!!    !$  call msg(sub_name//"parallel. num threads", nthreads)
      fWork => fu_work_array(nSpecies*nThreads)
      allocate(ar_species_single_src(nthreads), &
               & threadamounts(nSpecies,nThreads), stat = iTmp)
      if(iTmp /= 0)then
         call set_error('Failed allocation of source species for threads', sub_name)
      endif
      threadamounts(:,:) = 0D0
!!    !$OMP END MASTER
!!    !$OMP BARRIER
    amounts_tmp(1:nSpecies,1:nThreads) => fWork(1:nSpecies*nThreads)

    
    nullify(ar_species_single_src(iThread + 1)%pArSp)
    species_single_src => ar_species_single_src(iThread + 1)%pArSp
    amounts_single_src(1:nSpecies) => amounts_tmp(1:nSpecies, iThread + 1)

    if(.not. error)then
      amounts_single_src(1:nSpecies) = 0.0
    endif

    !
    ! Area source
    !
!!    !$OMP MASTER
    if(src%n_area > 0) call msg('Amounts emitted from area sources')
!!    !$OMP END MASTER
!!    !$OMP BARRIER

!!    !$OMP DO
    do i = 1, src%n_area
      !
      ! Get the total source amount for the species emitted from the source
      !
      if(error)cycle
      if (ifTimePresent) then
        call total_amt_species_unit(src%a_ptr(i)%a_src, .false., &  ! src, ifSlotRateOnly
                                  & species_single_src, &
                                  & nspecies_single_src, &
                                  & amounts_single_src, &
                                  & start, duration)
      else
        call total_amt_species_unit(src%a_ptr(i)%a_src, .false., &
                                  & species_single_src, &
                                  & nspecies_single_src, &
                                  & amounts_single_src)
      end if

      if(nspecies_single_src == 0)cycle  ! nothing from this source
      if(sum(amounts_single_src(1:nspecies_single_src)) == 0.0)cycle  ! nothing from this source
      !
      ! Having the source amounts and species, report and store them into the main array
      ! This is quick but affects the shared structures, has to be done in critical section
      !
!!!      !$OMP CRITICAL(asrc)
!     FLUSH (run_log_funit) 
      call report_amounts(fu_name(src%a_ptr(i)%a_src), fu_sector(src%a_ptr(i)%a_src), &
                        & species_single_src, amounts_single_src, nspecies_single_src)
!      if(error)cycle
!     FLUSH (run_log_funit) 
!call msg("iThread, nthreads", iThread, nthreads)
!call msg('Source and species: '// trim(fu_name(src%a_ptr(i)%a_src)) // '_' // trim(fu_sector(src%a_ptr(i)%a_src)))
!do j = 1, nspecies_single_src
!call msg('Amount for species('+fu_str(j)+')' + fu_name(fu_material(species_single_src(j))),amounts_single_src(j))
!end do
!call msg("TRANSP SPECIES:",  nSpecies, size(transport_species))
!do j = 1, nSpecies
! call msg("jTRANSP", j)
! call msg("TRANSP SPECIES: "//trim(fu_str(transport_species(j))))
!enddo
!     FLUSH (run_log_funit) 
      !
      ! Create the link to the transport species
      ! species_single_src(i) = transport_species(adaptor%isp(i))
      !
      if (.not. error) call create_adaptor(species_single_src, transport_species, adaptor)

     ! Collect the emitted amounts to the appropriate places
        do iSp = 1, nspecies_single_src
          threadamounts(adaptor%isp(iSp),iThread+1) = threadamounts(adaptor%isp(iSp),iThread+1) + &
             & real(amounts_single_src(iSp), 8) !Some want to see explicit conversion here
        end do

    end do  ! area sources
!!    !$OMP END DO
    !
    ! Bomb source
    !
!!    !$OMP BARRIER
!!    !$OMP MASTER
    if (src%n_bomb > 0)call msg('Amount from bomb sources')
!!    !$OMP END MASTER
!!    !$OMP BARRIER

!!    !$OMP DO
    do i = 1, src%n_bomb
      if(error)cycle
      call total_bomb_src_species_unit(src%b_ptr(i)%b_src, &
                                     & species_single_src, &
                                     & nspecies_single_src, &
                                     & amounts_single_src)
      if(error)cycle

      if(nspecies_single_src == 0)cycle  ! nothing from this source
      if(sum(amounts_single_src(1:nspecies_single_src)) == 0.0)cycle  ! nothing from this source

      call report_amounts(fu_name(src%b_ptr(i)%b_src), fu_sector(src%b_ptr(i)%b_src), &
                        & species_single_src, amounts_single_src, nspecies_single_src)
      !
      ! Create the link to the transport species
      ! species_single_src(i) = transport_species(adaptor%isp(i))
      !
      if (.not. error)  call create_adaptor(species_single_src, transport_species, adaptor)
!      if(error)cycle

     ! Collect the emitted amounts to the appropriate places
      do iSp = 1, nspecies_single_src
        threadamounts(adaptor%isp(iSp),iThread+1) = threadamounts(adaptor%isp(iSp),iThread+1) + &
           & real(amounts_single_src(iSp), 8) !Some want to see explicit conversion here
      end do


    end do ! bomb sources
!!    !$OMP END DO
    !
    ! Point source
    !
!!    !$OMP BARRIER
!!    !$OMP MASTER
    if(src%n_point > 0)call msg('Amount from point sources')
!!    !$OMP END MASTER
!!    !$OMP BARRIER

!!    !$OMP DO
    do i = 1, src%n_point
      if(error)cycle
      if (ifTimePresent) then
        call total_amt_species_unit(src%p_ptr(i)%p_src, &
                                  & species_single_src, &
                                  & nspecies_single_src, &
                                  & amounts_single_src, &
                                  & start, duration)
      else
        call total_amt_species_unit(src%p_ptr(i)%p_src, &
                                  & species_single_src, &
                                  & nspecies_single_src, &
                                  & amounts_single_src)
      end if
      if(error)cycle

      if(nspecies_single_src == 0)cycle  ! nothing from this source
      if(sum(amounts_single_src(1:nspecies_single_src)) == 0.0)cycle  ! nothing from this source

      call report_amounts(fu_name(src%p_ptr(i)%p_src), fu_sector(src%p_ptr(i)%p_src), &
                        & species_single_src, amounts_single_src, nspecies_single_src)
      
      if (.not. error) call create_adaptor(species_single_src, transport_species, adaptor)
!      if(error)cycle
     ! Collect the emitted amounts to the appropriate places
      do iSp = 1, nspecies_single_src
        threadamounts(adaptor%isp(iSp),iThread+1) = threadamounts(adaptor%isp(iSp),iThread+1) + &
           & real(amounts_single_src(iSp), 8) !Some want to see explicit conversion here
      end do

    end do ! point sources
!!    !$OMP END DO

    !
    ! Fire source
    !
!!    !$OMP BARRIER
!!    !$OMP MASTER
    if(src%n_fire > 0)call msg('Amount from fire sources')
!!    !$OMP END MASTER
!!    !$OMP BARRIER

!!    !$OMP DO
    do i = 1, src%n_fire
      if(error)cycle

      if (ifTimePresent) then
        call total_amt_species_unit(src%fire_ptr(i)%fire_src, &
                                  & species_single_src, &
                                  & nspecies_single_src, &
                                  & amounts_single_src, &
                                  & start, duration)
      else
        call total_amt_species_unit(src%fire_ptr(i)%fire_src, &
                                  & species_single_src, &
                                  & nspecies_single_src, &
                                  & amounts_single_src)
      end if
      if(error)cycle

      if(nspecies_single_src > 0)then
        call report_amounts(fu_name(src%fire_ptr(i)%fire_src), fu_sector(src%fire_ptr(i)%fire_src), &
                          & species_single_src, amounts_single_src, nspecies_single_src)
        !
        if(.not. error) call create_adaptor(species_single_src, transport_species, adaptor)

        ! Collect the emitted amounts to the appropriate places
        do iSp = 1, nspecies_single_src
           threadamounts(adaptor%isp(iSp),iThread+1) = threadamounts(adaptor%isp(iSp),iThread+1) + &
              & real(amounts_single_src(iSp), 8) !Some want to see explicit conversion here
        end do
      else
        call msg('No emission found in the source:' + fu_name(src%fire_ptr(i)%fire_src) + ',' + &
                                                    & fu_sector(src%fire_ptr(i)%fire_src))
      endif

    end do ! fire sources
!!    !$OMP END DO

    !
    ! Volcano source
    !
!!    !$OMP BARRIER
!!    !$OMP MASTER
    if(src%n_volc > 0)call msg('Amount from volcano sources')
!!    !$OMP END MASTER
!!    !$OMP BARRIER

!!    !$OMP DO
    do i = 1, src%n_volc
      if(error)cycle
      if (ifTimePresent) then
        call total_amt_species_unit(src%volc_ptr(i)%volc_src, &
                                  & species_single_src, &
                                  & nspecies_single_src, &
                                  & amounts_single_src, &
                                  & start, duration, level_missing)
      else
        call total_amt_species_unit(src%volc_ptr(i)%volc_src, &
                                  & species_single_src, &
                                  & nspecies_single_src, &
                                  & amounts_single_src, &
                                  & time_missing, interval_missing, level_missing)
      end if
      if(error)cycle

      if(nspecies_single_src == 0)cycle  ! nothing from this source
      if(sum(amounts_single_src(1:nspecies_single_src)) == 0.0)cycle  ! nothing from this source

      call report_amounts(fu_name(src%volc_ptr(i)%volc_src), fu_sector(src%volc_ptr(i)%volc_src), &
                        & species_single_src, amounts_single_src, nspecies_single_src)
      
      if (.not. error) call create_adaptor(species_single_src, transport_species, adaptor)
!      if(error)cycle
     ! Collect the emitted amounts to the appropriate places
      do iSp = 1, nspecies_single_src
        threadamounts(adaptor%isp(iSp),iThread+1) = threadamounts(adaptor%isp(iSp),iThread+1) + &
           & real(amounts_single_src(iSp), 8) !Some want to see explicit conversion here
      end do

    end do ! point sources
!!    !$OMP END DO

!!    !$OMP END PARALLEL



!call msg('Amounts from the inside 1,' + fu_str(nSpecies) + ',' + fu_str(nThreads), amounts(1:nSpecies))


!    do i = 1, nThreads
!      amounts(1:nSpecies) = amounts(1:nSpecies) + threadamounts(1:nSpecies,i)
!    end do

    amounts(1:nSpecies) =  sum( threadamounts(1:nSpecies,1:nThreads), DIM=2)  !Sum over threads


!call msg('Amounts from the inside 2,' + fu_str(nSpecies), amounts(1:nSpecies))

    deallocate(threadamounts)
    deallocate(ar_species_single_src)

    call free_work_array(fWork)


  contains
    
    subroutine report_amounts(name, sector, species_list, amounts, nspecies)
      implicit none
      character(len=*), intent(in) :: name, sector
      type(silam_species), dimension(:), intent(in) :: species_list
      real, dimension(:), intent(in) :: amounts
      integer, intent(in) :: nspecies
      
      character(len=fnlen) :: str1, str2
      INTEGER :: iUnit
      INTEGER, DIMENSION(2) :: fUnits
      
      fUnits(1:2) = (/6, run_log_funit/)
      
      str2=''
      if (nspecies > 30) then
        str1='Source has > 30 emission species, emission not reported'
      elseif(nspecies < 1)then
        str1='No species found in the source'
     else
         write(unit=str1, fmt='(30(A10,1x))') (/(fu_str(species_list(isp)), isp=1, nspecies)/)
         write(unit=str2, fmt='(30(G10.2,1x))') amounts(1:nspecies)
      endif

       
    do iUnit = 1,2
       if (smpi_global_rank /= 0 .and. iUnit==1) cycle !Be quiet at stdut

        write(funits(iUnit),'(A,/,A,/,A)') 'Source: '// trim(name) // '_' // trim(sector), &
             & trim(str1), trim(str2)
    enddo
    end subroutine report_amounts
      
  end subroutine amounts_from_src_species_unit


  !********************************************************************************
  
  subroutine typical_cnc_from_src_species_unit(src, transport_species, nspecies, cnc)
    !
    ! Typical concentrations are known for some source terms, which are related to natural processes
    !
    implicit none

    ! Imported parameters
    type(silam_source), intent(in) :: src
    type(silam_species), dimension(:), pointer :: transport_species
    integer, intent(in) :: nspecies ! transport
    real, dimension(:), intent(out) :: cnc
    
    ! Local variabels
    type(silam_species), dimension(:), pointer :: src_species
    real, dimension(:), pointer :: arCnc
    integer :: nSpeciesSrc, iSrc, j, iSp

    arCnc => fu_work_array()

    cnc(1:nSpecies) = real_missing
    !
    ! Bio VOC knows its typical concentrations
    !
    do iSrc = 1, src%n_bio_voc

      call typical_species_conc(src%bvoc_ptr(iSrc)%bvoc_src, src_species, nSpeciesSrc, arCnc)
      if(error)return

      do j = 1, nSpeciesSrc
        isp = fu_index(src_species(j), transport_species)
        if (isp < 1) then
          call set_error('Strange emitted species', 'typical_cnc_from_src_species_unit')
          return
        end if
        cnc(isp) = min(arCnc(j), cnc(iSp))
      end do
    end do
    !
    ! Pollen knows its typical concentrations
    !
    do iSrc = 1, src%n_pollen

      call typical_species_conc(src%pollen_ptr(iSrc)%pollen_src, src_species, nSpeciesSrc, arCnc)
      if(error)return

      do j = 1, nSpeciesSrc
        isp = fu_index(src_species(j), transport_species)
        if (isp < 1) then
          call set_error('Strange emitted species', 'typical_cnc_from_src_species_unit')
          return
        end if
        cnc(isp) = min(arCnc(j), cnc(iSp))
      end do
    end do
    !
    ! Sea salt knows its typical concentrations
    !
    do iSrc = 1, src%n_sea_salt

      call typical_species_conc(src%sslt_ptr(iSrc)%sslt_src, src_species, nSpeciesSrc, arCnc)
      if(error)return

      do j = 1, nSpeciesSrc
        isp = fu_index(src_species(j), transport_species)
        if (isp < 1) then
          call set_error('Strange emitted species', 'typical_cnc_from_src_species_unit')
          return
        end if
        cnc(isp) = min(arCnc(j), cnc(iSp))
      end do
    end do
    !
    ! Wind-blown dust knows its typical concentrations
    !
    do iSrc = 1, src%n_wb_dust

      call typical_species_conc(src%wbdust_ptr(iSrc)%wbdust_src, src_species, nSpeciesSrc, arCnc)
      if(error)return

      do j = 1, nSpeciesSrc
        isp = fu_index(src_species(j), transport_species)
        if (isp < 1) then
          call set_error('Strange emitted species', 'typical_cnc_from_src_species_unit')
          return
        end if
        cnc(isp) = min(arCnc(j), cnc(iSp))
      end do
    end do
    
    ! DMS knows its typical concentrations (?)
    !
    do iSrc = 1, src%n_dms
      call typical_species_conc(src%dms_ptr(iSrc)%dms_src, src_species, nSpeciesSrc, arCnc)
      if(error)return
      do j = 1, nSpeciesSrc
        isp = fu_index(src_species(j), transport_species)
        if (isp < 1) then
          call set_error('Strange emitted species', 'typical_cnc_from_src_species_unit')
          return
        end if
        cnc(isp) = min(arCnc(j), cnc(iSp))
      end do
    end do

    call free_work_array(arCnc)

  end subroutine typical_cnc_from_src_species_unit


  !********************************************************************************************

  logical function fu_emission_owned_quantity(src, quantity)
    !
    ! Checks whether the specific quantity is handled by the emission sources exclusively
    ! Such quantities may be related to emission species somewhere else
    !
    implicit none

    ! Imported parameters
    type(silam_source), intent(in) :: src
    integer, intent(in) :: quantity
    
    ! Local variables
    integer :: iSrc
    !
    ! A common quantity handled by all sources is emission intensity
    !
    if(quantity == emission_intensity_flag .or. quantity == emission_flux_flag)then
      fu_emission_owned_quantity = .true.
    else
      !
      ! Sources, which have some models inside can handle a few other quantities, have to ask them
      !
      fu_emission_owned_quantity = .false.
      do iSrc = 1, src%n_bio_voc
        fu_emission_owned_quantity = fu_bio_voc_emis_owned_quantity(src%bvoc_ptr(iSrc)%bvoc_src, quantity)
        if(fu_emission_owned_quantity)return
      end do
      do iSrc = 1, src%n_pollen
        fu_emission_owned_quantity = fu_pollen_emis_owned_quantity(src%pollen_ptr(iSrc)%pollen_src, quantity)
        if(fu_emission_owned_quantity)return
      end do
      do iSrc = 1, src%n_sea_salt
        fu_emission_owned_quantity = fu_sslt_emis_owned_quantity(src%sslt_ptr(iSrc)%sslt_src, quantity)
        if(fu_emission_owned_quantity)return
      end do
      do iSrc = 1, src%n_wb_dust
        fu_emission_owned_quantity = fu_wb_dust_emis_owned_quantity(src%wbdust_ptr(iSrc)%wbdust_src, quantity)
        if(fu_emission_owned_quantity)return
      end do
      do iSrc = 1, src%n_dms
        fu_emission_owned_quantity = fu_dms_emis_owned_quantity(src%dms_ptr(iSrc)%dms_src, quantity)
        if(fu_emission_owned_quantity)return
      end do
      do iSrc = 1, src%n_volc
        fu_emission_owned_quantity = fu_volc_emis_owned_quantity(src%volc_ptr(iSrc)%volc_src, quantity)
        if(fu_emission_owned_quantity)return
      end do
    endif

  end function fu_emission_owned_quantity


  !*****************************************************************

  integer function fu_n_disp_grd_cells_src_inv(src) result(nbr)
    !
    ! Returns the total number of cells of dispersion grid that emit something
    ! Since there can be many sources emitting essentially to the same grid cells
    ! we have to take care of the problem by limiting the number with the total grid 
    ! size
    !
    implicit none

    ! Imported parameters with intent IN
    type(silam_source), intent(in) :: src

    ! Local declarations
    integer :: i

    nbr = 0
    !
    ! assume full overlap of fires from different sources
    !
    if(src%n_fire + src%n_area + src%n_point + src%n_bomb + src%n_volc > 0)then
      do i = 1, src%n_fire
        nbr = max(nbr, fu_n_fires(src%fire_ptr(i)%fire_src))
      end do
      !
      !  Add point sources
      !
      nbr = nbr + src%n_point + src%n_bomb + src%n_volc
      !
      !  and all area-sources cells
      !
      do i=1,src%n_area
        nbr = nbr + fu_nbr_of_disp_grd_cells(src%a_ptr(i)%a_src)
        if(error)then
          nbr = int_missing
          return
        endif
        if(nbr > nx_dispersion * ny_dispersion)then
          nbr = nx_dispersion * ny_dispersion
          return
        endif
      end do
      if(nbr < 1 .and. smpi_global_tasks < 2 )then
        call msg('Strange number of dispersion-grid cells after summing all sources:',nbr)
        call msg_warning('Strange number of dispersion-grid cells after summing all sources', &
                       & 'fu_n_disp_grd_cells_src_inv')
      endif
    endif   ! if any sources with defined number of cells are included

  end function fu_n_disp_grd_cells_src_inv


  !*****************************************************************

  integer function fu_source_id_nbr_of_source(src, &
          & i_point, i_area, i_bomb, i_fire, i_sea_salt, i_pollen, i_bio_voc, i_wb_dust, i_dms, i_volc)
    !
    ! Returns the source number. Reason: all sources are enumerated
    ! sequencially, so that the source can, in fact, be refered by its
    ! number without other information. But so far this number is 
    ! copied to the particles in the pollution cloud. 
    !
    ! NOTE. One and only one index may be reasonable. The other one MUST be
    ! negative or zero
    !
    implicit none

    ! Imported parameters with intent IN
    type(silam_source), intent(in) :: src
    integer, intent(in) :: i_point, i_area, i_bomb, i_fire, i_sea_salt, i_pollen, i_bio_voc, i_wb_dust, i_dms, i_volc

    ! Stupidity check
    if(.not.defined(src))then
      call set_error('Undefined source given','fu_source_nbr_of_source')
      return
    endif
    if(i_area > 0)then
      fu_source_id_nbr_of_source = fu_source_id_nbr(src%a_ptr(i_area)%a_src)
      if(fu_fails(i_point+i_bomb+i_fire+i_sea_salt+i_pollen+i_bio_voc+i_wb_dust+i_dms+i_volc == 0, &
                & 'Strange indices given: more than one source type','fu_source_id_nbr_of_source'))return
      return
    endif
    if(i_bio_voc > 0)then
      fu_source_id_nbr_of_source = fu_source_id_nbr(src%bvoc_ptr(i_bio_voc)%bvoc_src)
      if(fu_fails(i_point+i_area+i_bomb+i_fire+i_sea_salt+i_pollen+i_wb_dust+i_dms+i_volc == 0, &
                & 'Strange indices given: more than one source type','fu_source_id_nbr_of_source'))return
      return
    endif
    if(i_bomb > 0)then
      fu_source_id_nbr_of_source = fu_source_id_nbr(src%b_ptr(i_bomb)%b_src)
      if(fu_fails(i_point+i_area+i_fire+i_sea_salt+i_pollen+i_bio_voc+i_wb_dust+i_dms+i_volc == 0, &
                & 'Strange indices given: more than one source type','fu_source_id_nbr_of_source'))return
      return
    endif
    if(i_fire > 0)then
      fu_source_id_nbr_of_source = fu_source_id_nbr(src%fire_ptr(i_fire)%fire_src)
      if(fu_fails(i_point+i_area+i_bomb+i_sea_salt+i_pollen+i_bio_voc+i_wb_dust+i_dms+i_volc == 0, &
                & 'Strange indices given: more than one source type','fu_source_id_nbr_of_source'))return
      return
    endif
    if(i_point > 0)then
      fu_source_id_nbr_of_source = fu_source_id_nbr(src%p_ptr(i_point)%p_src)
      if(fu_fails(i_area+i_bomb+i_fire+i_sea_salt+i_pollen+i_bio_voc+i_wb_dust+i_dms+i_volc == 0, &
                & 'Strange indices given: more than one source type','fu_source_id_nbr_of_source'))return
      return
    endif
    if(i_pollen > 0)then
      fu_source_id_nbr_of_source = fu_source_id_nbr(src%pollen_ptr(i_pollen)%pollen_src)
      if(fu_fails(i_point+i_area+i_bomb+i_fire+i_sea_salt+i_bio_voc+i_wb_dust+i_dms+i_volc == 0, &
                & 'Strange indices given: more than one source type','fu_source_id_nbr_of_source'))return
      return
    endif
    if(i_sea_salt > 0)then
      fu_source_id_nbr_of_source = fu_source_id_nbr(src%sslt_ptr(i_sea_salt)%sslt_src)
      if(fu_fails(i_point+i_area+i_bomb+i_fire+i_pollen+i_bio_voc+i_wb_dust+i_dms+i_volc == 0, &
                & 'Strange indices given: more than one source type','fu_source_id_nbr_of_source'))return
      return
    endif
    if(i_wb_dust > 0)then
      fu_source_id_nbr_of_source = fu_source_id_nbr(src%wbdust_ptr(i_wb_dust)%wbdust_src)
      if(fu_fails(i_point+i_area+i_bomb+i_fire+i_sea_salt+i_pollen+i_bio_voc+i_dms+i_volc == 0, &
                & 'Strange indices given: more than one source type','fu_source_id_nbr_of_source'))return
      return
    endif
    if(i_dms > 0)then
      fu_source_id_nbr_of_source = fu_source_id_nbr(src%dms_ptr(i_dms)%dms_src)
      if(fu_fails(i_point+i_area+i_bomb+i_fire+i_sea_salt+i_pollen+i_bio_voc+i_wb_dust+i_volc == 0, &
                & 'Strange indices given: more than one source type','fu_source_id_nbr_of_source'))return
      return
    endif
    if(i_volc > 0)then
      fu_source_id_nbr_of_source = fu_source_id_nbr(src%volc_ptr(i_volc)%volc_src)
      if(fu_fails(i_point+i_area+i_bomb+i_fire+i_sea_salt+i_pollen+i_bio_voc+i_wb_dust+i_dms == 0, &
                & 'Strange indices given: more than one source type','fu_source_id_nbr_of_source'))return
      return
    endif
    call set_error('None of the source indices is > 0','fu_source_id_nbr_of_source')

  end function fu_source_id_nbr_of_source


  !*******************************************************************************

  logical function fu_check_outTemplate(chOutTemplate, FilesArrangement, iSourceId)
    !
    ! Checks the rules for consistency
    !
    implicit none

    ! Imported parameters with intent IN
    character(len=*), intent(in) :: chOutTemplate
    integer, intent(in) :: FilesArrangement
    integer, intent(in), optional :: iSourceId

    fu_check_outTemplate = .false.

    !
    ! Consistency means that every file will have unique name,so
    ! - if sources are split - their IDs must be in the names
    ! - if every hour/day/month/year a new file is started - corresponding
    !   time must exist in the names
    !
    if(present(iSourceId))then
      select case(iSourceId)
        case(iNoId) ! Mix sources together
        case(iSrcNmId, iSrcSectorId, iSrcNmSectorId, iSrcTimeZoneId)
          if(index(chOutTemplate, trim(fu_silam_source_template_item())) == 0)then
            call msg('Output template:' + chOutTemplate + ', src template:' + &
                                                                & fu_silam_source_template_item())
            call set_error('Source name is not in output template', 'fu_check_outTemplate')
            return
          endif
        case default
          call set_error('Strange source Id switch','fu_check_outTemplate')
          return
      end select
    endif

    !
    ! If one requires a new file every while - he has to ensure that the 
    ! template contains enough varying items. Below checking is a bit too tuff
    ! but it at ensures that there always will be at least one varying parameter.
    ! Example: if we start run from 20:00 for 12 hours only and set 
    ! only valid time hours to vary - we will totally mix the order of files. So, 
    ! it is better not to think about such tricks and ask for full definition.
    ! Exception is - forecast length in hours, which is always enough - it is like 
    ! a counter.
    ! Also, the template may not include parameters, which are varying inside file.
    ! For example, if daily_new_file has several hours inside and template includes
    ! %h2 then GrADS will look for the new file every file, which is not the case.
    ! Therefore, the templates included into the file name must strictly meet the
    ! file arrangements.
    !
    select case(FilesArrangement)

      case(all_in_one)  ! No templates are allowed
        
        if(index(chOutTemplate,trim(fu_forecast_length_templ_str()))>0 .or. &
         & index(chOutTemplate,trim(fu_valid_time_hour_templ_str()))>0.or. &
         & index(chOutTemplate,trim(fu_valid_time_day_templ_str()))>0.or. &
         & index(chOutTemplate,trim(fu_valid_time_month_templ_str()))>0.or. &
         & index(chOutTemplate,trim(fu_valid_time_year_templ_str()))>0)then
          call set_error('ALL_IN_ONE does not allow time templates', &
                       & 'fu_check_outTemplate')
          return
        endif

      case(hourly_new_file) ! Either forecast length or all others together

        if(index(chOutTemplate,trim(fu_forecast_length_templ_str()))==0)then
          if((index(chOutTemplate,trim(fu_valid_time_hour_templ_str()))==0).or. &
           & (index(chOutTemplate,trim(fu_valid_time_day_templ_str()))==0).or. &
           & (index(chOutTemplate,trim(fu_valid_time_month_templ_str()))==0).or. &
           & (index(chOutTemplate,trim(fu_valid_time_year_templ_str()))==0))then
            call set_error('Not enough varying items in template for time series', &
                         & 'fu_check_outTemplate')
            return
          endif
        endif

      case(daily_new_file) ! Must be %y, %m, %d, the rest forbidden

        if((index(chOutTemplate,trim(fu_valid_time_day_templ_str()))==0).or. &
         & (index(chOutTemplate,trim(fu_valid_time_month_templ_str()))==0).or. &
         & (index(chOutTemplate,trim(fu_valid_time_year_templ_str()))==0))then
          call set_error('Not enough varying items in template for time series', &
                       & 'fu_check_outTemplate')
          return
        endif
        if(index(chOutTemplate,trim(fu_forecast_length_templ_str()))>0 .or. &
         &(index(chOutTemplate,trim(fu_valid_time_hour_templ_str()))>0))then
          call set_error('DAILY_NEW_FILE does not allow hourly templates', &
                       & 'fu_check_outTemplate')
          return
        endif

      case(monthly_new_file) ! Must be %y, %m, the rest forbidden

        if((index(chOutTemplate,trim(fu_valid_time_month_templ_str()))==0).or. &
         & (index(chOutTemplate,trim(fu_valid_time_year_templ_str()))==0))then
          call set_error('Not enough varying items in template for time series', &
                       & 'fu_check_outTemplate')
          return
        endif
        if(index(chOutTemplate,trim(fu_forecast_length_templ_str()))>0 .or. &
         &(index(chOutTemplate,trim(fu_valid_time_hour_templ_str()))>0) .or. &
         &(index(chOutTemplate,trim(fu_valid_time_day_templ_str()))>0))then
          call set_error('MONTHLY_NEW_FILE does not allow daily/hourly templates', &
                       & 'fu_check_outTemplate')
          return
        endif

      case(yearly_new_file)

        if((index(chOutTemplate,trim(fu_valid_time_year_templ_str()))==0))then
          call set_error('Not enough varying items in template for time series', &
                       & 'fu_check_outTemplate')
          return
        endif
        if(index(chOutTemplate,trim(fu_forecast_length_templ_str()))>0 .or. &
         &(index(chOutTemplate,trim(fu_valid_time_hour_templ_str()))>0) .or. &
         &(index(chOutTemplate,trim(fu_valid_time_day_templ_str()))>0) .or. &
         &(index(chOutTemplate,trim(fu_valid_time_month_templ_str()))>0))then
          call set_error('YEARLY_NEW_FILE allows only year in templates', &
                       & 'fu_check_outTemplate')
          return
        endif

      case default
        call set_error('Unknown time series file arrangement','fu_check_outTemplate')
        return

    end select

    !
    ! SILAM does not have clear definition of the analysis time because it actually 
    ! does not make the analysis. So, to avoid an ambiguity, the analysis time
    ! is NOT allowed in the output template. Still in place: valid time for the fields,
    ! initial time (start of the computations) and the 
    ! forecast length = valid time - initial time (negative for the inverse run)
    !
    if((index(chOutTemplate,trim(fu_anal_time_hour_templ_str()))>0).or. &
     & (index(chOutTemplate,trim(fu_anal_time_day_templ_str()))>0).or. &
     & (index(chOutTemplate,trim(fu_anal_time_month_templ_str()))>0).or. &
     & (index(chOutTemplate,trim(fu_anal_time_year_templ_str()))>0))then
      call set_error('Analysis time is not allowed in output template', &
                   & 'fu_check_outTemplate')
      return
    endif

    fu_check_outTemplate = .true.

  end function fu_check_outTemplate


  !**********************************************************************
  
  subroutine make_source_id_mapping(em_source, meteoMarketPtr, emis_grid, nlSetup)
    !
    ! If emission is split along whatever rules different from source name and sector
    ! one needs to have a map saying what grid cell belongs to what source index.
    ! This mapping is created by this subroutine
    !
    implicit none
    
    ! Imported parameters
    type(silam_source), intent(inout) :: em_source
    type(mini_market_of_stacks), pointer :: meteoMarketPtr
    type(silja_grid), intent(in) :: emis_grid
    type(Tsilam_namelist), intent(in) :: nlSetup
    
    ! local variables
    type(Tsilam_nl_item_ptr), dimension(:), pointer :: pItems
    integer :: nx, ny, iTmp, nItems
    
    !
    ! Procedure of making the mapping depends on the type of the split
    !
    select case(em_source%iSrcIdType)
      case(iNoId, iSrcNmId, iSrcSectorId, iSrcNmSectorId)
        em_source%ifUseSourceIdMapping = .false.  ! for these types of output split no mapping needed

      case(iSrcTimeZoneId)
        !
        ! Time-zone related split does require mapping: each time zone or timezone group
        ! get their own index. This index is distributed according to the maps of time zones
        !
        call grid_dimensions(emis_grid, nx, ny)
        allocate(em_source%arSourceIdMapping(nx, ny), &
               & em_source%chSplitNames(number_of_time_zones), stat = iTmp)
        if(iTmp /= 0)then
          call report(emis_grid)
          call set_error('Failed allocation of time zone mapping for above emission grid', &
                       & 'make_source_id_mapping')
          return
        endif
        nullify(pItems)
        call get_items(nlSetup, 'time_zone_group', pItems, nItems)
        call make_time_zone_mapping(em_source%arSourceIdMapping, em_source%chSplitNames, &
                                  & meteoMarketPtr, emis_grid, pItems, nItems)
        em_source%ifUseSourceIdMapping = .true.
        
      case default
        call set_error('Unknown type of source split','make_source_id_mapping')
    end select
    
  end subroutine make_source_id_mapping
  
  
  !**********************************************************************
  
  integer function fu_num_sources_total(em_source) result(count_all)
    implicit none
    type(silam_source), intent(in) :: em_source
    
    count_all = (em_source%n_area + em_source%n_bio_voc + em_source%n_bomb + em_source%n_fire + &
               & em_source%n_point + em_source%n_pollen + em_source%n_sea_salt + em_source%n_wb_dust + &
               & em_source%n_dms + em_source%n_volc)
    
  end function fu_num_sources_total

  !**************************************************************************
  !
  ! Data assimilation subroutines for parameter assimilation
  ! Pretty much the intermediates to pass the requests there and back
  !
  !**************************************************************************
  
  subroutine assimilation_request_emis(src, arParams, iLastFilledParam)
    !
    ! Returns the list of parameters to be assimilated for this volcanic source
    ! Also stores their index in the global array of parameters
    !
    implicit none
    
    ! Imported parameters
    type(silam_source), intent(inout) :: src
    type(assimParameter), dimension(:), intent(inout) :: arParams     ! (nParamsToAssim)
    integer, intent(inout) :: iLastFilledParam

    ! local variables
    integer :: iSrc
    
    ! Scan all the sources one-by-one and collect what they want
    !
    !do iSrc=1,src%n_area
    !  call assimilation_request_area(src%a_ptr(iSrc)%a_src, arParams, iLastFilledParam)
    !end do
    !do iSrc=1,src%n_bio_voc
    !  call assimilation_request_bvoc(src%bvoc_ptr(iSrc)%bvoc_src, arParams, iLastFilledParam)
    !end do
    !do iSrc=1,src%n_bomb
    !  call assimilation_request_bomb(src%b_ptr(iSrc)%b_src, arParams, iLastFilledParam)
    !end do
    !do iSrc=1,src%n_dms
    !  call assimilation_request_dms(src%dms_ptr(iSrc)%dms_src, arParams, iLastFilledParam)
    !end do
    !do iSrc=1,src%n_fire
    !  call assimilation_request_fire(src%fire_ptr(iSrc)%fire_src, arParams, iLastFilledParam)
    !end do
    !do iSrc=1,src%n_point
    !  call assimilation_request_point(src%p_ptr(iSrc)%p_src, arParams, iLastFilledParam)
    !end do
    !do iSrc=1,src%n_pollen
    !  call assimilation_request_pollen(src%pollen_ptr(iSrc)%pollen_src, arParams, iLastFilledParam)
    !end do
    !do iSrc=1,src%n_sea_salt
    !  call assimilation_request_sea_salt(src%sslt_ptr(iSrc)%sslt_src, arParams, iLastFilledParam)
    !end do
    do iSrc=1,src%n_volc
      call assimilation_request_volc_src(src%volc_ptr(iSrc)%volc_src, arParams, iLastFilledParam)
    end do
    !do iSrc=1,src%n_wb_dust
    !  call assimilation_request_wb_dust(src%wbdust_ptr(iSrc)%wbdust_src, arParams, iLastFilledParam)
    !end do

  end subroutine assimilation_request_emis

  
  !*******************************************************************************
  
  subroutine observe_params_emis(src, arParams, now)
    !
    ! Provides the values of the current assimilated parameters for the interface array
    ! Obeys their place in the interface array. Mapping already exists and serves as the source of data
    !
    implicit none
    
    ! Imported parameters
    type(silam_source), intent(inout) :: src
    type(assimParameter), dimension(:), intent(inout) :: arParams     ! (nParamsToAssim)
    type(silja_time), intent(in) :: now
  
    ! Local variables
    integer :: iSrc
    !
    ! Scan the sources one0by-one
    !
    !do iSrc=1,src%n_area
    !  call observe_params_ (src%a_ptr(iSrc)%a_src, arParams, now)
    !end do
    !do iSrc=1,src%n_bio_voc
    !  call observe_params_ (src%bvoc_ptr(iSrc)%bvoc_src, arParams, now)
    !end do
    !do iSrc=1,src%n_bomb
    !  call observe_params_ (src%b_ptr(iSrc)%b_src, arParams, now)
    !end do
    !do iSrc=1,src%n_dms
    !  call observe_params_ (src%dms_ptr(iSrc)%dms_src, arParams, now)
    !end do
    !do iSrc=1,src%n_fire
    !  call observe_params_ (src%fire_ptr(iSrc)%fire_src, arParams, now)
    !end do
    !do iSrc=1,src%n_point
    !  call observe_params_ (src%p_ptr(iSrc)%p_src, arParams, now)
    !end do
    !do iSrc=1,src%n_pollen
    !  call observe_params_ (src%pollen_ptr(iSrc)%pollen_src, arParams, now)
    !end do
    !do iSrc=1,src%n_sea_salt
    !  call observe_params_ (src%sslt_ptr(iSrc)%sslt_src, arParams, now)
    !end do
    do iSrc=1,src%n_volc
      call observe_params_volc_src(src%volc_ptr(iSrc)%volc_src, arParams, now)
    end do
    !do iSrc=1,src%n_wb_dust
    !  call observe_params_ (src%wbdust_ptr(iSrc)%wbdust_src, arParams, now)
    !end do
    
  end subroutine observe_params_emis
  
  
  !*******************************************************************************
  
  subroutine inject_params_emis(src, arParams, now)
    !
    ! Stores the new parameter values into the corresponding places of the source.
    ! Note that the actuaql papameters are in the mapping structures, rest is just pointers
    !
    implicit none
    
    ! Imported parameters
    type(silam_source), intent(inout) :: src
    type(assimParameter), dimension(:), intent(inout) :: arParams     ! (nParamsToAssim)
    type(silja_time), intent(in) :: now

    ! Local variables
    integer :: iSrc
    !
    ! Scan the sources one0by-one
    !
    !do iSrc=1,src%n_area
    !  call inject_params_inject_params_ (src%a_ptr(iSrc)%a_src, arParams, now)
    !end do
    !do iSrc=1,src%n_bio_voc
    !  call inject_params_(src%bvoc_ptr(iSrc)%bvoc_src, arParams, now)
    !end do
    !do iSrc=1,src%n_bomb
    !  call inject_params_(src%b_ptr(iSrc)%b_src, arParams, now)
    !end do
    !do iSrc=1,src%n_dms
    !  call inject_params_(src%dms_ptr(iSrc)%dms_src, arParams, now)
    !end do
    !do iSrc=1,src%n_fire
    !  call inject_params_(src%fire_ptr(iSrc)%fire_src, arParams, now)
    !end do
    !do iSrc=1,src%n_point
    !  call inject_params_(src%p_ptr(iSrc)%p_src, arParams, now)
    !end do
    !do iSrc=1,src%n_pollen
    !  call inject_params_(src%pollen_ptr(iSrc)%pollen_src, arParams, now)
    !end do
    !do iSrc=1,src%n_sea_salt
    !  call inject_params_(src%sslt_ptr(iSrc)%sslt_src, arParams, now)
    !end do
    do iSrc=1,src%n_volc
      call inject_params_volc_src(src%volc_ptr(iSrc)%volc_src, arParams, now)
    end do
    !do iSrc=1,src%n_wb_dust
    !  call inject_params_(src%wbdust_ptr(iSrc)%wbdust_src, arParams, now)
    !end do

  end subroutine inject_params_emis
  
  
  
  !*****************************************************************
  
  SUBROUTINE report_source(src, ifDetailed)

    IMPLICIT NONE

    ! Imported parameters with intent IN:
    TYPE(silam_source), INTENT(in) :: src
    logical, intent(in) :: ifDetailed

    integer :: i

    IF (.NOT.defined(src)) THEN
      CALL set_error('undefined source given','report_source')
      RETURN
    END IF

    if(ifDetailed)then
      do i=1,src%n_area
        CALL report(src%a_ptr(i)%a_src)
      end do
      do i=1, src%n_bio_voc
        call report(src%bvoc_ptr(i)%bvoc_src)
      end do
      do i = 1,src%n_bomb
        CALL report(src%b_ptr(i)%b_src)
      end do
      do i = 1,src%n_fire
        CALL report(src%fire_ptr(i)%fire_src)
      end do
      do i=1,src%n_point
        CALL report(src%p_ptr(i)%p_src)
      end do
      do i=1, src%n_pollen
        call report(src%pollen_ptr(i)%pollen_src)
      end do
      do i=1, src%n_sea_salt
        call report(src%sslt_ptr(i)%sslt_src)
      end do
      do i=1, src%n_wb_dust
        call report(src%wbdust_ptr(i)%wbdust_src)
      end do
      do i=1, src%n_dms
        call report(src%dms_ptr(i)%dms_src)
      end do
      do i=1, src%n_volc
        call report(src%volc_ptr(i)%volc_src)
      end do
    else
      do i=1, size(src%src_info_lst)
        if(src%src_info_lst(i)%iSrcNbr > 0)then
          CALL msg('Source:' + src%src_info_lst(i)%chSrcNm + ',' + src%src_info_lst(i)%chSectorNm)
        endif
      end do
    endif

  END SUBROUTINE report_source  


  !*****************************************************************

  subroutine source_2_map_general (em_src, dataPtr, id, ifWholeVertical, iSrcId, iAccuracy, ifRandomise)
    !
    ! General function, which drops the source pointed by the indexSrc to the
    ! map described by id. The field itself is pointed by dataPtr. Paticular 
    ! variable is described by the id as well - by the substance name
    !
    implicit none

    ! Imported parameters
    type(silam_source), intent(inout) :: em_src
    real, dimension(:), pointer :: dataPtr
    type(silja_field_id), intent(in) :: id
    logical, intent(in) :: ifWholeVertical, ifRandomise
    integer, intent(in) :: iSrcId ! Requested source ID number
    integer, intent(in) :: iAccuracy

    ! Local variables
    integer :: iSrc
    logical :: ifFound

    !
    ! Zeroing the output dataPtr has to be done here, not in individual sources.
    ! Otherwise, the sum-up of the sources will fail
    !
    dataPtr(1:fu_number_of_gridpoints(fu_grid(id))) = 0.

    !
    ! Simply find the source pointed by the iSrcId and call its own 
    ! projecting function
    !
    ifFound = .false.
    do iSrc = 1, em_src%n_point + em_src%n_area + em_src%n_bomb
      if(em_src%src_info_lst(iSrc)%iIdNbr == iSrcId)then
        ifFound = .true.
!        call msg('Source:' + em_src%src_info_lst(iSrc)%chId + ',' + &
!               & fu_str(fu_species(id)), em_src%src_info_lst(iSrc)%iIdNbr)
        if(em_src%src_info_lst(iSrc)%iSrcType == area_source)then
          call source_2_map(em_src%a_ptr(em_src%src_info_lst(iSrc)%iSrcNbr)%a_src, &
                          & dataPtr, id, ifWholeVertical, iAccuracy, ifRandomise)
        elseif(em_src%src_info_lst(iSrc)%iSrcType == bomb_source)then
          call source_2_map(em_src%b_ptr(em_src%src_info_lst(iSrc)%iSrcNbr)%b_src, &
                          & dataPtr, id, ifWholeVertical)
        elseif(em_src%src_info_lst(iSrc)%iSrcType == point_source)then
          call source_2_map(em_src%p_ptr(em_src%src_info_lst(iSrc)%iSrcNbr)%p_src, &
                          & dataPtr, id, ifWholeVertical)
        elseif(em_src%src_info_lst(iSrc)%iSrcType == volcano_source)then
          call source_2_map_volc_src(em_src%volc_ptr(em_src%src_info_lst(iSrc)%iSrcNbr)%volc_src, &
                                   & dataPtr, id, ifWholeVertical)
        elseif(em_src%src_info_lst(iSrc)%iSrcType == fire_source .or. &
             & em_src%src_info_lst(iSrc)%iSrcType == sea_salt_source .or. &
             & em_src%src_info_lst(iSrc)%iSrcType == pollen_source .or. &
             & em_src%src_info_lst(iSrc)%iSrcType == bio_voc_source .or. &
             & em_src%src_info_lst(iSrc)%iSrcType == wind_blown_dust_source .or. &
             & em_src%src_info_lst(iSrc)%iSrcType == dms_source)then
          call msg_warning('The type of the source is not supported:' + em_src%src_info_lst(iSrc)%chId)
        else
          call msg('Strange source type:',em_src%src_info_lst(iSrc)%iSrcType)
          call set_error('Strange source type','source_2_map_general')
          return
        endif
      endif
    end do
    
    if(.not. ifFound)then
      call msg('The source ID is not found:', iSrcId)
      call msg('Available IDs:')
      do iSrc = 1, em_src%n_point + em_src%n_area + em_src%n_bomb !+ em_src%n_fire
        call msg('ID:',em_src%src_info_lst(iSrc)%iIdNbr)
      end do
      call set_error('The source ID is not found in the source list','source_2_map_general')
    endif

  end subroutine source_2_map_general


  !*****************************************************************

  function fu_grid_of_source(src) result(grid)
    !
    ! Returns the source grid, which is copied from the area source.
    ! Should the area sources be absent or contain several grids - the 
    ! grid_missing is returned
    !
    implicit none

    ! Return value of the function
    type(silja_grid) :: grid

    ! Imported parameter
    type(silam_source), intent(in) :: src

    ! Local variables
    integer :: i

    if(defined(src)) then
      if(src%n_area >= 1)then
        grid = fu_grid(src%a_ptr(1)%a_src)
        do i=2,src%n_area
          if(.not. (fu_grid(src%a_ptr(i)%a_src) == grid))then
            grid = grid_missing
            return
          endif
        enddo
      endif
    endif
    grid = grid_missing

  end function fu_grid_of_source


  !**************************************************************************

  subroutine link_general_src_to_species(src, species)
    !
    ! Having the cocktails created, we should establish the shortcut links between the 
    ! source descriptors and cocktail. The link goes via descr%iEmisCocktSpeciesMapping 
    ! and  descr%factor_to_basic_unit.
    ! Note that cocktail is "anything" and the only requirement is that it has all the species
    ! existing in the source inventory.
    ! That has to happen in two steps. Firstly, we establish these links using the 
    ! single descriptor per source. Then, these connections are distributed to each time slot
    !
    implicit none

    ! Imported parameters
    type(silam_source), intent(inout) :: src
    type(silam_species), dimension(:), pointer :: species

    ! Local variables
    integer :: iSrc

    !
    ! Area sources
    !
    do iSrc = 1, src%n_area
      call link_source_to_species(species, src%a_ptr(iSrc)%a_src)
      if(error)return
    end do
    !
    ! Bio VOC sources
    !
    do iSrc = 1, src%n_bio_voc
      call link_source_to_species(species, src%bvoc_ptr(iSrc)%bvoc_src)
      if(error)return
    end do
    !
    ! Then bomb sources. Note only one time slot
    !
    do iSrc = 1, src%n_bomb
      call link_source_to_species(species, src%b_ptr(iSrc)%b_src)
      if(error)return
    end do
    !
    ! Fire sources
    !
    do iSrc = 1, src%n_fire
      call link_source_to_species(species, src%fire_ptr(iSrc)%fire_src)
      if(error)return
    end do
    !
    ! Point sources
    !
    do iSrc = 1, src%n_point
      call link_source_to_species(species, src%p_ptr(iSrc)%p_src)
      if(error)return
    end do
    !
    ! Pollen sources
    !
    do iSrc = 1, src%n_pollen
      call link_source_to_species(species, src%pollen_ptr(iSrc)%pollen_src)
      if(error)return
    end do
    !
    ! Sea salt sources
    !
    do iSrc = 1, src%n_sea_salt
      call link_source_to_species(species, src%sslt_ptr(iSrc)%sslt_src)
      if(error)return
    end do
    !
    ! Wind-blown dust sources
    !
    do iSrc = 1, src%n_wb_dust
      call link_source_to_species(species, src%wbdust_ptr(iSrc)%wbdust_src)
      if(error)return
    end do
    !
    ! DMS sources
    !
    do iSrc = 1, src%n_dms
      call link_source_to_species(species, src%dms_ptr(iSrc)%dms_src)
      if(error)return
    end do
    !
    ! VOLCANO sources
    !
    do iSrc = 1, src%n_volc
      call link_volc_src_to_species(species, src%volc_ptr(iSrc)%volc_src)
      if(error)return
    end do

  end subroutine link_general_src_to_species


  !***************************************************************************

  subroutine link_emission_2_dispersion(em_src, emMap_species, meteoMarket, gridMeteo, gridDisp, &
                                      & verticalMeteo, verticalDispersion, iAccuracy, ifRandomise)
    !
    ! Present strategy is that the set_emission calls deal with sources
    ! directly at each time step. Reason: time and composition variations, which
    ! cannot be pre-defined
    ! 
    ! Here we somewhat reduce the cost of the exercise by explicit reprojecting the 
    ! sources to dispersion_grid. Costly action if done carefully! 
    ! The routine is supposed to be called once at the start of the run and thus we
    ! can allow it.
    ! Also, this is the place to initialise the internal source fields.
    !
    implicit none

    ! Imported parameter
    TYPE(silam_source), INTENT(inout) :: em_src
    type(silam_species), dimension(:), pointer :: emMap_species
    type(mini_market_of_stacks), intent(in) :: meteoMarket
    type(silja_grid), intent(in) :: gridMeteo, gridDisp  ! to that all sources should be projected
    type(silam_vertical), intent(in) :: verticalMeteo, verticalDispersion
    integer, intent(in) :: iAccuracy
    logical, intent(in) :: ifRandomise

    ! Local declarations
    integer :: i
    type(silam_vertical) :: vertical_metric
    type(silja_rp_1d), dimension(:), pointer :: tmpArraySet
    type(rng_type) :: rng
    type(THorizInterpStruct), pointer ::  interpCoefMeteo2DispHoriz
    type(silja_field), pointer :: TZidxFld

    !
    ! Stupidity check
    !
    if(.not. defined(gridMeteo))then
      call set_error('undefined meteo grid given','link_emission_2_dispersion')
      return
    endif
    if(.not. defined(gridDisp))then
      call set_error('undefined dispersion grid given','link_emission_2_dispersion')
      return
    endif
    if(.not. defined(em_src))then
      call set_error('undefined source given','link_emission_2_dispersion')
      return
    endif

    !
    ! Link the source terms and the emission mass map - for all source types
    !
    call msg('Linking the source to emission map..')
    call link_source_to_species(em_src, emMap_species)
    if(error)return

    ! Dispersion vertical can be hybrid. For the source we'll define a matching vertical
    ! in meters. This ensures that emission density (per vertical meter) is conserved in
    ! the vertical reprojection.
    !
    call vert_to_metric(verticalDispersion, vertical_metric)
    if (error) return

    !
    ! Scan the sources one-by-one and request them to initialise the internal structures
    ! 
    ! INVENTORY types: store new projections to the given grid and vertical (each source 
    ! has the second position/cell pointer). Strictly speaking, re-projection of verticals 
    ! can be meteo-dependent but here we will use crude parameters of the standard atmosphere. 
    ! Accuracy of source verticals is anyway very crude, so plus-minus a bit does not change 
    ! anything. The only exception will be the plume-rise routine, which will be actually 
    ! called each time step with all dynamic parameters.
    !
    ! AREA sources

    !Interpolation for timezone index map
       interpCoefMeteo2DispHoriz => fu_horiz_interp_struct(gridMeteo, gridDisp, &
                               &  nearest_point, .False.,  0, notAllowed, .False.)
     if (em_src%ifNeedsTZindexMap) then
       TZidxFld =>  fu_sm_simple_field(meteoMarket, met_src_missing, timezone_index_flag, &
                         & level_missing, single_time_stack_flag)
     else
        TZidxFld => null()
     endif


    !
    !$OMP PARALLEL DEFAULT (NONE) &
    !$OMP & PRIVATE(tmpArraySet, i, rng) &
    !$OMP & SHARED(em_src, gridDisp, verticalDispersion, vertical_metric, error, &
    !$OMP        & iAccuracy, ifRandomise, emMap_species, gridMeteo, verticalMeteo,  interpCoefMeteo2DispHoriz, TZidxFld)
    !
    ! Since some of the sources will require pretty much the same temporaries, which need
    ! to be allocated, let's prepare space for them
    !
    call get_work_arrays_set(2, 100, tmpArraySet)    ! any size: the first source will correct this
    if(error) call set_error('Work array set failed','link_emission_2_dispersion')

    !$OMP DO SCHEDULE (DYNAMIC, 1)
    do i=1, em_src%n_area
      if (error) cycle
      call rng_init(rng, i)
      call project_a_src_second_grd(em_src%a_ptr(i)%a_src, gridDisp, verticalDispersion, vertical_metric, &
                              & iAccuracy, tmpArraySet, ifRandomise, rng)
      if(error)call set_error('Source reprojection failed:' + fu_name(em_src%a_ptr(i)%a_src),&
                            & 'link_emission_2_dispersion')
      call prepare_src_vert_params(em_src%a_ptr(i)%a_src, verticalDispersion, vertical_metric)
      if(error)call set_error('Source vertical failed:' + fu_name(em_src%a_ptr(i)%a_src),&
                            & 'link_emission_2_dispersion')

      call init_a_src_TZ_index(em_src%a_ptr(i)%a_src, interpCoefMeteo2DispHoriz, TZidxFld, gridDisp)
      if(error)call set_error('TZ init failed:' + fu_name(em_src%a_ptr(i)%a_src),&
                              & 'link_emission_2_dispersion')
    enddo
    !$OMP END DO

    call free_work_array(tmpArraySet)
    
    !$OMP END PARALLEL
    if(error)return
    
    !
    ! BIOGENIC VOC source. No action needed
    !
    do i=1, em_src%n_bio_voc
!      call msg('bvoc_src -> 2d grd:',i)
      call source_2_second_grid(em_src%bvoc_ptr(i)%bvoc_src, gridDisp, verticalDispersion, &
                              & vertical_metric, iAccuracy)
      if(error)return
!      call msg('done')
    enddo
    !
    ! BOMB sources
    !
    do i=1,em_src%n_bomb
!      call msg('b_src -> 2d grd:',i)
      call source_2_second_grid(em_src%b_ptr(i)%b_src, gridDisp)
!      call msg('done')
      if(error)return
    enddo
    !
    ! FIRE sources
    !
    do i=1, em_src%n_fire
!      call msg('p_src -> 2d grd:',i)
      call source_2_second_grid(em_src%fire_ptr(i)%fire_src, gridDisp, verticalDispersion, &
                              & vertical_metric, iAccuracy)
      if(error)return
!      call msg('Time params')
!      call prepare_src_vert_params(em_src%p_ptr(i)%p_src, vertical_metric)
!      call msg('done')
      if(error)return
    enddo
    !
    ! POINT sources
    !
    do i=1, em_src%n_point
!      call msg('p_src -> 2d grd:',i)
      call project_p_src_2_grids(em_src%p_ptr(i)%p_src, gridMeteo, gridDisp)
      if(error)return
!      call msg('Time params')
      call prepare_src_vert_params_p_src(em_src%p_ptr(i)%p_src, verticalMeteo, verticalDispersion, &
                                       & vertical_metric)
!      call msg('done')
      if(error)return
    enddo
    !
    ! POLLEN source
    !
    do i=1, em_src%n_pollen
!      call msg('a_src -> 2d grd:',i)
      call source_2_second_grid(em_src%pollen_ptr(i)%pollen_src, gridDisp, verticalDispersion, &
                              & vertical_metric, iAccuracy)
      if(error)return
!      call msg('done')
    enddo
    !
    ! SEA SALT source
    !
    do i=1, em_src%n_sea_salt
!      call msg('a_src -> 2d grd:',i)
      call source_2_second_grid(em_src%sslt_ptr(i)%sslt_src, gridDisp, verticalDispersion, &
                              & vertical_metric, iAccuracy)
      if(error)return
!      call msg('done')
    enddo
    !
    ! WIND-BLOWN DUST source
    !
    do i=1, em_src%n_wb_dust
!      call msg('wb_src -> 2d grd:',i)
      call source_2_second_grid(em_src%wbdust_ptr(i)%wbdust_src, gridDisp, verticalDispersion,  &
                              & vertical_metric, iAccuracy)
      if(error)return
!      call msg('done')
    enddo
 !
    ! DMS source
    !
    do i=1, em_src%n_dms
      call source_2_second_grid(em_src%dms_ptr(i)%dms_src, gridDisp, verticalDispersion,  &
                              & vertical_metric, iAccuracy)
      if(error)return
    enddo
    !
    ! VOLCANO sources
    !
    do i=1, em_src%n_volc
!      call msg('p_src -> 2d grd:',i)
      call project_volc_src_2_grids(em_src%volc_ptr(i)%volc_src, gridMeteo, gridDisp)
      if(error)return
!      call msg('Time params')
      call prepare_volc_src_vert_params(em_src%volc_ptr(i)%volc_src, verticalMeteo, verticalDispersion, &
                                      & vertical_metric)
!      call msg('done')
      if(error)return
    enddo

  end subroutine link_emission_2_dispersion


  !**********************************************************************************

  subroutine init_emission_internal_fields(src, pDispersionMarket, start_time)
    !
    ! Initialises the internal fields for the source terms - mainly sending them
    ! to dispersion market. Has to be harmonised with the add_input_needs subroutines
    !
    implicit none

    ! Imported parameters
    TYPE(silam_source), INTENT(inout) :: src
    type(mini_market_of_stacks), pointer :: pDispersionMarket
    type(silja_time), intent(in) :: start_time

    ! Local variables
    integer :: iSrc

    !
    ! INVENTORY sources: point, area, bomb. No action needed
    !

    !
    ! Bio VOC source. No grid but plenty of fields to be initialised
    !
    do iSrc = 1, src%n_bio_voc
      call init_emission_bio_voc(src%bvoc_ptr(iSrc)%bvoc_src, pDispersionMarket, start_time)
      if(error)return
    end do
    !
    ! FIRE source. No grid but plenty of fields to be initialised and stored into dispersion stack
    !
    do iSrc = 1, src%n_fire
      call init_emission_fire(src%fire_ptr(iSrc)%fire_src)
      if(error)return
    end do
    !
    ! POLLEN source. No grid but plenty of fields to be initialised and stored into dispersion stack
    !
    do iSrc = 1, src%n_pollen
      call init_emission_pollen(src%pollen_ptr(iSrc)%pollen_src, pDispersionMarket, start_time)
      if(error)return
    end do
    !
    ! SEA SALT source. Does not have grid to reproject but has internal stuff to be innitialised
    !
    do iSrc = 1, src%n_sea_salt
      call init_emission_sea_salt(src%sslt_ptr(iSrc)%sslt_src)
      if(error)return
    end do
    !
    ! WIND-BLOWN source. Does not have grid to reproject but has internal stuff to be innitialised
    !
    do iSrc = 1, src%n_wb_dust
      call init_emission_wb_dust(src%wbdust_ptr(iSrc)%wbdust_src, pDispersionMarket, start_time)
      if(error)return
    end do

    ! DMS source. Does not have grid to reproject but has internal stuff to be innitialised
    !
    do iSrc = 1, src%n_dms
      call init_emission_dms(src%dms_ptr(iSrc)%dms_src)
      if(error)return
    end do

    ! FIRES source. Complicated: there are many of them (e.g., daily) and they share the land-use,
    !               which is ~1-3km resolution covering up to a globe. That huge map is to be
    !               read once, used when reading the sources, then destroyed. Therefore, all
    !               fire sources have to be initialised at once.
    !
    do iSrc = 1, src%n_fire
      call init_emission_fire(src%fire_ptr(iSrc)%fire_src)
      if(error)return
    enddo

  end subroutine init_emission_internal_fields


  !******************************************************************

  subroutine inject_emission_eulerian(em_src, &
                                    & mapEmis, mapCoordX_emis, mapCoordY_emis, mapCoordZ_emis, &
                                    & met_buf, disp_buf, &
                                    & now, timestep, &
                                    & ifSpeciesMoment, &
                                    & fMassInjected, &
                                    & interpCoefMeteo2DispHoriz, &
                                    & ifMetHorizInterp, &
                                    & interpCoefMeteo2DispVert, &
                                    & ifMetVertInterp)
    !
    ! Part of Eulerian environment.
    ! Computes and injects emission fluxes to the concentration and subgrid-info
    ! fields.
    ! All emission is injected into emisMapPtr (Tmass_map), which afterwards is sent to 
    ! transport mass map. Actually, the only reason for that is to have emission to the output.
    ! There are some other niceties, such as fewer number of vertical layers, which is yet to be 
    ! implemented, etc.
    !
    implicit none

    ! Imported parameters
    TYPE(silam_source), INTENT(inout) :: em_src
    type(Tmass_map), intent(inout) :: mapEmis, mapCoordX_emis, mapCoordY_emis, mapCoordZ_emis
    TYPE(Tfield_buffer), POINTER :: met_buf, disp_buf
    type(silja_time), intent(in) :: now
    type(silja_interval), intent(in) :: timestep
    real(r8k), dimension(:), target, INTENT(inout) :: fMassInjected
    type(THorizInterpStruct), intent(in) ::  interpCoefMeteo2DispHoriz
    type(TVertInterpStruct), intent(in) ::  interpCoefMeteo2DispVert
    logical, intent(in) :: ifMetHorizInterp
    logical, intent(in) :: ifMetVertInterp
    logical, intent(in) :: ifSpeciesMoment

    ! Local variables
    integer :: iSrc, iThread, nThreads, iTmp
    real, dimension(:), pointer :: fWork
    real, dimension(:,:), pointer :: fMassTimeCommon
    integer, dimension(:), pointer :: iSlotStart, iSlotEnd
    real(r8k), dimension(:), pointer :: fMassInjectedThread

    !$ if (.not. allocated(dInjectedMassThread)) then
         nThreads = 1
         !$OMP PARALLEL DEFAULT(shared) 
         !$OMP MASTER
    !$       nThreads = omp_get_num_threads()
         !$OMP END MASTER   
         !$OMP END PARALLEL
         if ( nThreads > 1) then
            call msg("Allocating multithread stuff for OMP emission")
            allocate(dInjectedMassThread(1:mapEmis%nSpecies,1:nThreads-1), stat=iTmp)
            if (iTmp /= 0) then 
               call set_error("Allocation failed!", "inject_emission_eulerian")
            endif
         endif
    !$ endif




!call msg('Grand total of emission map before emission:',sum(mapEmis%arM(:,:,:,:,:)))

    !
    ! Emission is inserted into the emission mass map and then sent to dispersion mass map
    !
    ! Speedup: set the columnValid and gridValid to false before injection.
    !
    mapEmis%ifGridValid = .false.
    mapEmis%ifColumnValid = .false.
    !
    ! Scan the individual sources looking for meteo-independent explicitly given fluxes
    ! Should any is found, add to the concentration map together with the source position
    ! 
    if(em_src%n_area > 0) call prepare_inject_a_src(met_buf)
    !
    ! If some of the sources are meteo-dependent, have to prepare the corresponding pointers
    
    ! Area source can pick stuff only from dispersion buffers
    if(em_src%arMeteoDepndencies(1) /= int_missing) &
                              & call prepare_injection_meteodep(disp_buf, em_src%arMeteoDepndencies)
    
    islotStart => fu_work_int_array()
    islotEnd   => islotStart(em_src%n_area+1:2*em_src%n_area)
    fWork => fu_work_array(em_src%n_area*max_descriptors_in_source)
    fMassTimeCommon(1:max_descriptors_in_source, 1:em_src%n_area) => &
           &  fWork(1:max_descriptors_in_source  * em_src%n_area) 

    !$OMP PARALLEL DEFAULT(none) &
    !$OMP & shared (em_src, mapEmis, mapCoordX_emis, mapCoordY_emis, mapCoordZ_emis, &
    !$OMP &           met_buf, now, timestep, ifSpeciesMoment, fMassInjected, &
    !$OMP &           interpCoefMeteo2DispHoriz, islotEnd,  fMassTimeCommon, &
    !$OMP &           ifMetHorizInterp, interpCoefMeteo2DispVert, ifMetVertInterp, &
    !$OMP &           nThreads, ny_dispersion, error, dinjectedmassthread, islotStart) &
    !$OMP & private (iSrc, iThread, iTmp, fMassInjectedThread) 
    !$OMP 
    !$OMP MASTER
!!    !$ call msg("OMP area source emission: NareaSources Nthreads", em_src%n_area, omp_get_num_threads())
    fMassInjectedThread => fMassInjected ! By default points to global one
    nThreads = 1
    !$ nThreads = omp_get_num_threads()
    !$OMP END MASTER
    !$OMP BARRIER
    iThread = 0
    !$  iThread = omp_get_thread_num()
      
    if (iThread /= 0) then
      fMassInjectedThread => dInjectedMassThread(:,iThread)
      fMassInjectedThread(:) = 0
    endif
    
    ! Count fMassTimeCommon without touchig massmaps.  Sources are independent
    !$OMP DO 
    do iSrc=1,em_src%n_area                                        !========== AREA with no field
      if(error) cycle
      if (.not. em_src%a_ptr(iSrc)%a_src%ifFieldGiven ) then
          call count_emission_area_source_cellist(em_src%a_ptr(iSrc)%a_src, &
                                       & now, timestep, &
                                       & fMassTimeCommon(:,iSrc), &
                                       & islotStart(iSrc), iSlotEnd(iSrc))
      endif
    enddo
    !$OMP END DO
    
    !$OMP BARRIER

    ! Inject stripes
    do iSrc=1,em_src%n_area                                        !========== AREA with no field
      if(error) cycle
      if (.not. em_src%a_ptr(iSrc)%a_src%ifFieldGiven ) then
          if (islotStart(iSrc) >= iSlotEnd(iSrc)) cycle
          call inject_emission_area_source_cellist(em_src%a_ptr(iSrc)%a_src, &
                                     & mapEmis, mapCoordX_emis, mapCoordY_emis, mapCoordZ_emis, &
                                     & met_buf, &
                                     & now, &
                                     & ifSpeciesMoment, &
                                     & fMassInjectedThread, &
                                     & interpCoefMeteo2DispHoriz, &
                                     & ifMetHorizInterp, &
                                     & interpCoefMeteo2DispVert, &
                                     & ifMetVertInterp, fMassTimeCommon(:,iSrc), &
                                     & islotStart(iSrc), iSlotEnd(iSrc), iThread, nThreads)
      endif
    enddo
    !$OMP END PARALLEL
    call free_work_array(islotStart)
    call free_work_array(fWork)

    do iThread = 1, nThreads-1
      fMassInjected(1:mapEmis%nSpecies) = fMassInjected(1:mapEmis%nSpecies) + dInjectedMassThread(:,iThread)
    enddo

!!    !$ call msg("OMP Done with area sources"  )
    if(error)return


    do iSrc=1,em_src%n_area                                        !========== AREA With field
      if ( em_src%a_ptr(iSrc)%a_src%ifFieldGiven ) then
        call inject_emission_area_source_field(em_src%a_ptr(iSrc)%a_src, &
                                     & mapEmis, mapCoordX_emis, mapCoordY_emis, mapCoordZ_emis, &
                                     & met_buf, &
                                     & now, timestep, &
                                     & ifSpeciesMoment, &
                                     & fMassInjected, &
                                     & interpCoefMeteo2DispHoriz, &
                                     & ifMetHorizInterp, &
                                     & interpCoefMeteo2DispVert, &
                                     & ifMetVertInterp)
        if(error) return
      endif
    enddo

    do iSrc=1,em_src%n_bio_voc                                      !========= Bio VOC
      call compute_emission_for_bio_voc(em_src%bvoc_ptr(iSrc)%bvoc_src, &
                                      & met_buf, disp_buf, & 
                                      & now, &      ! current time
                                      & timestep, & ! model time step
                                      & interpCoefMeteo2DispHoriz, ifMetHorizInterp, &
                                      & ifSpeciesMoment, &
                                      & mapEmis, mapCoordX_emis, mapCoordY_emis, mapCoordZ_emis, &  ! Output
                                      & fMassInjected) 
      if(error)return
    end do

    do iSrc=1,em_src%n_bomb                                        !========= BOMB
      call inject_emission_euler_b_src(em_src%b_ptr(iSrc)%b_src, &
                                        & mapEmis, mapCoordX_emis, mapCoordY_emis, mapCoordZ_emis, &
                                        & met_buf, disp_buf, &
                                        & now, timestep, &
                                        & ifSpeciesMoment, &
                                        & fMassInjected, &
                                        & interpCoefMeteo2DispHoriz, &
                                        & ifMetHorizInterp)
      if(error)return
    enddo
    
    if(em_src%n_fire>0)call prepare_inject_fire_src(met_buf)
    do iSrc=1,em_src%n_fire                                       !========== FIRES
      call inject_emission_euler_fire_src(em_src%fire_ptr(iSrc)%fire_src, &
                                        & mapEmis, mapCoordX_emis, mapCoordY_emis, mapCoordZ_emis, &
                                        & met_buf, &
                                        & now, timestep, &
                                        & interpCoefMeteo2DispHoriz, &
                                        & ifMetHorizInterp, &
                                        & interpCoefMeteo2DispVert, &
                                        & ifMetVertInterp, &
                                        & ifSpeciesMoment, &
                                        & fMassInjected)
      if(error)return
    enddo

    if(em_src%n_point > 0)call prepare_inject_p_src(met_buf)
    do iSrc=1,em_src%n_point                                       !========== POINT
      call inject_emission_euler_p_src(em_src%p_ptr(iSrc)%p_src, &
                                     & mapEmis, mapCoordX_emis, mapCoordY_emis, mapCoordZ_emis, &
                                     & met_buf, disp_buf, &
                                     & now, timestep, &
                                     & ifSpeciesMoment, &
                                     & fMassInjected, &
                                     & interpCoefMeteo2DispHoriz, &
                                     & ifMetHorizInterp, &
                                     & interpCoefMeteo2DispVert, &
                                     & ifMetVertInterp)
      if(error)return
    enddo

    do iSrc=1,em_src%n_pollen                                      !========= POLLEN
      call compute_emission_for_pollen(em_src%pollen_ptr(iSrc)%pollen_src, &
                                     & met_buf, disp_buf, & 
                                     & now, &      ! current time
                                     & timestep, & ! model time step
                                     & interpCoefMeteo2DispHoriz, ifMetHorizInterp, &
                                     & ifSpeciesMoment, &
                                     & mapEmis, mapCoordX_emis, mapCoordY_emis, mapCoordZ_emis, &
                                     & fMassInjected, & ! Output
                                     & em_src%arSourceIdMapping, em_src%ifUseSourceIdMapping)
      if(error)return
    end do

    do iSrc=1,em_src%n_sea_salt                                    !========= SEA SALT
      call compute_emission_for_sea_salt(em_src%sslt_ptr(iSrc)%sslt_src, &
                                       & met_buf, disp_buf, & 
                                       & now, &                   ! current time
                                       & timestep, &              ! model time step
                                       & interpCoefMeteo2DispHoriz, ifMetHorizInterp, &
                                       & ifSpeciesMoment, &
                                       & mapEmis, mapCoordX_emis, mapCoordY_emis, mapCoordZ_emis, & ! Output
                                       & fMassInjected)           ! output
      if(error)return
    end do

    do iSrc=1,em_src%n_wb_dust                                    !========= WIND-BLOWN DUST
      call compute_emission_for_wb_dust(em_src%wbdust_ptr(iSrc)%wbdust_src, &
                                      & met_buf, disp_buf, & 
                                      & now, &                   ! current time
                                      & timestep, &              ! model time step
                                      & interpCoefMeteo2DispHoriz, ifMetHorizInterp, &
                                      & ifSpeciesMoment, &
                                      & mapEmis, mapCoordX_emis, mapCoordY_emis, mapCoordZ_emis, & ! Output
                                      & fMassInjected)           ! output
      if(error)return
    end do
    do iSrc=1,em_src%n_dms                                    !========= DMS
      call compute_emission_for_dms(em_src%dms_ptr(iSrc)%dms_src, &
                                  & met_buf, disp_buf, & 
                                  & now, &                   ! current time
                                  & timestep, &              ! model time step
                                  & interpCoefMeteo2DispHoriz, ifMetHorizInterp, &
                                  & ifSpeciesMoment, &
                                  & mapEmis, mapCoordX_emis, mapCoordY_emis, mapCoordZ_emis, & ! Output
                                  & fMassInjected)           ! output
      if(error)return
    end do

    do iSrc=1,em_src%n_volc                                       !========== VOLCANO
      call inject_emission_euler_volc_src(em_src%volc_ptr(iSrc)%volc_src, &
                                        & mapEmis, mapCoordX_emis, mapCoordY_emis, mapCoordZ_emis, &
                                        & met_buf, disp_buf, &
                                        & now, timestep, &
                                        & ifSpeciesMoment, &
                                        & fMassInjected, &
                                        & interpCoefMeteo2DispHoriz, &
                                        & ifMetHorizInterp, &
                                        & interpCoefMeteo2DispVert, &
                                        & ifMetVertInterp)
      if(error)return
    enddo

!call msg('Grand total of emission map after emission:',sum(mapEmis%arM(:,:,:,:,:)))

  end subroutine inject_emission_eulerian


  !****************************************************************************************
  
  subroutine inject_emission_lagrangian(em_src, &
                                      & lpSet, &
                                      & arLowMassThresh, &
                                      & ChemRunSetup, &  ! translate emission species to transport
                                      & met_buf, disp_buf, &
                                      & now, timestep, &
                                      & fMassInjected, &
                                      & interpCoefMeteo2DispHoriz, ifMetHorizInterp, &
                                      & interpCoefMeteo2DispVert, ifMetVertInterp)
    !
    ! Part of Lagrangian environment.
    ! Computes and injects emission fluxes as puffs, which are further transported as lagrngian 
    ! particles.
    ! ATTENTION. Does NOT use emission mass map - for evident reasons. Sources are directly sent to 
    !            lagrangian particles with corresponding transition from source species to transport 
    !            species jumping over the emission species
    !
    implicit none

    ! Imported parameters
    TYPE(silam_source), INTENT(in) :: em_src
    type(Tlagrange_particles_set), INTENT(inout) :: lpSet
    real, dimension(:), INTENT(in) :: arLowMassThresh
    integer, dimension(:), pointer :: arStatus
    type(TchemicalRunSetup), pointer :: ChemRunSetup
    TYPE(Tfield_buffer), POINTER :: met_buf, disp_buf
    type(silja_time), intent(in) :: now
    type(silja_interval), intent(in) :: timestep
    real(r8k), dimension(:),  intent(inout) :: fMassInjected
    type(THorizInterpStruct), pointer ::  interpCoefMeteo2DispHoriz
    type(TVertInterpStruct), pointer ::  interpCoefMeteo2DispVert
    logical, intent(in) :: ifMetHorizInterp
    logical, intent(in) :: ifMetVertInterp

    ! Local variables
    integer :: ix,iy,iz,iSrc, iSubst
    real :: fCellTotal
    character(len=*), parameter :: subname = 'inject_emission_lagrangian'
    !
    ! Scan the individual sources looking for meteo-independent explicitly given fluxes
    ! Should any is found, add to the concentration map together with the source position
    !
    lpSet%nNewPart = 0 ! Reset counter of new particles

    do iSrc=1,em_src%n_area                                        !========== AREA
      call set_error('Lagrangian emission for area sources is not defined',subname)
    enddo

    do iSrc=1,em_src%n_bio_voc                                      !========= Bio VOC
      call set_error('Lagrangian emission for BVOC sources is not defined',subname)
    end do

    do iSrc=1,em_src%n_bomb                                        !========= BOMB
      call inject_emission_lagr_b_src(em_src%b_ptr(iSrc)%b_src, &
                                      & lpSet, arLowMassThresh,  & ! Lagrangian
                                      & ChemRunSetup,  &  ! translate emission species to transport
                                      & met_buf, disp_buf, &
                                      & now, timestep, &
                                      & fMassInjected)
      if(error)return
    enddo
    
    do iSrc=1,em_src%n_fire                                        !========= FIRE
      call set_error('Lagrangian emission for fire sources is not defined',subname)
    enddo
    
    if(em_src%n_point > 0)call prepare_inject_p_src(met_buf)
    do iSrc=1,em_src%n_point                                       !========== POINT
      call inject_emission_lagr_p_src(em_src%p_ptr(iSrc)%p_src, &
                                     & lpSet, arLowMassThresh, & ! Lagrangian
                                     & ChemRunSetup,  &   ! translate emission species to transport
                                     & met_buf, disp_buf, &
                                     & now, timestep, &
                                     & fMassInjected)
      if(error)return
    enddo

    do iSrc=1,em_src%n_volc                                       !========== VOLCANO
      call inject_emission_lagr_volc_src(em_src%volc_ptr(iSrc)%volc_src, &
                                       & lpSet, arLowMassThresh, & ! Lagrangian
                                       & ChemRunSetup,  &   ! translate emission species to transport
                                       & met_buf, disp_buf, &
                                       & now, timestep, &
                                       & fMassInjected)
      if(error)return
    enddo

    do iSrc=1,em_src%n_pollen                                      !========= POLLEN
      call set_error('Lagrangian emission for pollen sources is not defined',subname)
    end do

    do iSrc=1,em_src%n_sea_salt                                    !========= SEA SALT
      call set_error('Lagrangian emission for sea salt sources is not defined',subname)
    end do

    do iSrc=1,em_src%n_wb_dust                                    !========= WIND-BLOWN DUST
      call set_error('Lagrangian emission for wb_dust sources is not defined',subname)
    end do

    if (fu_fails(em_src%n_dms == 0, 'No lagrangian emission for DMS', subname)) return

  end subroutine inject_emission_lagrangian

END MODULE source_terms_general


