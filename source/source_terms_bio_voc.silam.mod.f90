MODULE source_terms_bio_voc
  !
  ! The biogenic dynamic emission sub-module, after Guenther et al, (1993). Obtained from 
  ! the group of Dimitris Melas (University of Thessaloniki): 
  ! Anastasia Poupkou, Thodoris Giannaros, Kostas Markakis, personal communication
  !
  ! All units: SI, unless otherwise stated
  !
  ! Code owner Mikhail Sofiev mikhail.sofiev@fmi.fi
  !
  ! Language: ANSI FORTRAN-90 (or close to it)
  !
  use cocktail_basic

  implicit none
  private

  !
  ! PUBLIC routines of sea salt
  !
  public fill_bio_voc_src_from_namelist
  public reserve_bio_voc_source
  public init_emission_bio_voc
  public add_inventory_bio_voc_src
  public add_input_needs
  public link_source_to_species
  public create_source_containing_grid
  public source_2_second_grid
  public compute_emission_for_bio_voc
  public fu_bio_voc_emis_owned_quantity
  public fu_name
  public fu_source_nbr
  public fu_source_id_nbr
  public typical_species_conc
  public report

  !
  ! Private routines of the sea salt source
  !
  private add_input_needs_bio_voc_src
  private link_bvoc_src_to_species
  private create_src_cont_grd_bvoc_src
  private project_bvoc_src_second_grd
  private fu_source_id_nbr_of_bvoc_src
  private fu_source_nbr_of_bvoc_src
  private fu_bio_voc_source_name
  private typical_species_cnc_bvoc
  private report_bio_voc_src

  !
  ! Private subs of sea salt source
  !

  interface add_input_needs
    module procedure add_input_needs_bio_voc_src
  end interface

  interface link_source_to_species
    module procedure link_bvoc_src_to_species
  end interface

  interface create_source_containing_grid
    module procedure create_src_cont_grd_bvoc_src
  end interface

  interface source_2_second_grid
    module procedure project_bvoc_src_second_grd
  end interface

  interface fu_source_nbr
    module procedure fu_source_nbr_of_bvoc_src
  end interface

  interface fu_source_id_nbr
    module procedure fu_source_id_nbr_of_bvoc_src
  end interface

  interface fu_name
    module procedure fu_bio_voc_source_name
  end interface

  interface typical_species_conc
    module procedure typical_species_cnc_bvoc
  end interface

  interface report
    module procedure report_bio_voc_src
  end interface


  !--------------------------------------------------------------------
  !
  !  The biogenic VOC source term
  !
  TYPE silam_bio_voc_source
    PRIVATE
    CHARACTER(len=clen) :: src_nm, sector_nm        ! Name of the area source and sector
    integer ::  src_nbr, id_nbr          ! A source and id numbers in a WHOLE source list
    integer :: emisMethod, swSource ! Method, source of sw radiation:sw_down or cc
    character(len=fnlen) :: chSrcMaskFile           ! To limit source area
    logical :: ifMasked
    character(len=fnlen), dimension(:), allocatable :: chLandUseMetaDataFNms  ! for land use
    real, dimension(:,:,:,:), allocatable :: arEmisFactor ! (nSpeciesTypes,nx_disp,ny_disp,month)
    integer :: ind_emis_isoprene, ind_emis_monoterpene
    real :: factor_emis_isoprene, factor_emis_monoterpene
    integer :: nLevsDispVert, nSpecies
    type(silam_vertical) :: vertLevsDispVert
    real, dimension(:), allocatable :: levFractDispVert, fzDisp
    type(silam_species), dimension(2) :: species ! Two species for now
    type(chemical_adaptor) :: adaptor
    type(silja_logical) :: defined
  END TYPE silam_bio_voc_source

  type bvoc_src_ptr
    type(silam_bio_voc_source) :: bvoc_src
  end type bvoc_src_ptr
  public bvoc_src_ptr
  
  !
  ! Parameters for biogenic emission
  !
  integer, public, parameter :: emisGuenther_updated_v1 = 601
  
  integer, private, parameter :: totBiomass  = 1             ! indices in the arEmisFactor
  integer, private, parameter :: emsIsoprene = 2
  integer, private, parameter :: emsMonoterpTempr      = 3
  integer, private, parameter :: emsMonoterpTemprLight = 4
  integer, private, parameter :: emsOtherVOC           = 5

CONTAINS

  
  !***********************************************************************

  subroutine fill_bio_voc_src_from_namelist(nlSetup, bvoc_src, src_data_dir)
    !
    ! Creates the set of rules, in accordance with which the computations are to
    ! be made. Input - SILAM namelist
    !
    implicit none

    !Imported parameters
    type(silam_bio_voc_source), intent(inout) :: bvoc_src
    type(Tsilam_namelist), pointer :: nlSetup 
    character(len=*), intent(in) :: src_data_dir

    ! Local variables
    logical :: srcOK
    integer :: nFiles, iFile, iTmp, nSubstEmitted
    type(Tsilam_nl_item_ptr), dimension(:), pointer :: pItems
    type(silam_species) :: speciesTmp
    character(len=substNmLen), dimension(max_species) :: chSubstNm

    bvoc_src%defined = silja_false

    if(.not.defined(nlSetup))then
      call set_error('Undefined namelist given','fill_bio_voc_src_from_namelist')
      return
    endif

    !
    ! Emission index type. It also determines the list of emitted species
    !
    bvoc_src%nSpecies = 0

    select case(fu_str_u_case(fu_content(nlSetup,'bvoc_emission_method')))
      case ('GUENTHER_MODIFIED_V1')
        bvoc_src%emisMethod = emisGuenther_updated_v1
        
      case default
        call set_error('Unknown emission method:' + fu_content(nlSetup,'bvoc_emission_method'), &
                     & 'fill_bio_voc_src_from_namelist')
        return
    end select

    select case(fu_str_l_case(fu_content(nlSetup,'SW_down_method')))
      case ('surf_sw_down_radiation')
        bvoc_src%swSource = surf_sw_down_radiation_flag
        call msg("fill_bio_voc_src_from_namelist uses surf_sw_down_radiation_flag")
      case ('total_cloud')
        bvoc_src%swSource = total_cloud_cover_flag
        call msg("fill_bio_voc_src_from_namelist uses total_cloud_cover_flag")
      case default
        call msg("SW_down_method can be one of surf_sw_down_radiation or total_cloud!")
        call set_error('Unknown SW_down_method:' + fu_content(nlSetup,'SW_down_method'), &
                     & 'fill_bio_voc_src_from_namelist')
        return
    end select

    select case (fu_str_u_case(fu_content(nlsetup, 'if_emit_isoprene')))
      case ('YES')
        call msg('Isoprene emission requested')
        call set_species(speciesTmp, fu_get_material_ptr('C5H8'), in_gas_phase)
        if (error) return
        bvoc_src%nspecies = bvoc_src%nspecies + 1
        bvoc_src%species(bvoc_src%nspecies) = speciesTmp
        bvoc_src%ind_emis_isoprene = bvoc_src%nspecies
        bvoc_src%factor_emis_isoprene = fu_content_real(nlsetup, 'isoprene_emission_factor')
        if(bvoc_src%factor_emis_isoprene == real_missing)then
          bvoc_src%factor_emis_isoprene = 0.5
          call msg('No factor for isoprene emission given - use default:', bvoc_src%factor_emis_isoprene)
        else
          call msg('Isoprene emission factor:', bvoc_src%factor_emis_isoprene)
        endif     
      case ('NO')
        call msg('Isoprene emission not requested')
        bvoc_src%ind_emis_isoprene = int_missing
      case default
        call set_error('Missing if_emit_isoprene in bio VOC source definition', &
                     & 'fill_bio_voc_src_from_namelist')
        return
    end select
    
    select case(fu_str_u_case(fu_content(nlsetup, 'if_emit_monoterpene')))
      case ('YES')
        call set_species(speciesTmp, fu_get_material_ptr('C5H8_2'), in_gas_phase)
        if (error) return
        bvoc_src%nspecies = bvoc_src%nspecies + 1
        bvoc_src%species(bvoc_src%nspecies) = speciesTmp
        bvoc_src%ind_emis_monoterpene = bvoc_src%nspecies ! either 1 or 2
        call msg('Monoterpene emission requested')
        bvoc_src%factor_emis_monoterpene = fu_content_real(nlsetup, 'monoterpene_emission_factor')
        if(bvoc_src%factor_emis_monoterpene .eps. real_missing)then
          bvoc_src%factor_emis_monoterpene = 1.0
          call msg('No factor for monoterpene emission given - use default:', bvoc_src%factor_emis_monoterpene)
        else
          call msg('Monoterpene emission factor:', bvoc_src%factor_emis_monoterpene)
        endif     
      case ('NO')
        bvoc_src%ind_emis_monoterpene = int_missing
        call msg('Monoterpene emission not requested')
      case default
        call set_error('Missing if_emit_monoterpene in bio VOC source definition', &
                     & 'fill_bio_voc_src_from_namelist')  
    end select

    srcOK = .true.

    !
    ! Store the land-use metadata file names
    !
    nullify(pItems)
    call get_items(nlSetup, 'land_use_meta_data_file', pItems, nFiles)
    if(error)return
    if(nFiles < 1)then
      call set_error('No supplementary files are given','fill_bio_voc_src_from_namelist')
      return
    endif
    allocate(bvoc_src%chLandUseMetaDataFNms(nFiles), stat=iFile)
    if(iFile /= 0)then
      call set_error('Failed to allocate land use meta data file names','fill_bio_voc_src_from_namelist')
      return
    endif
    do iFile = 1, nFiles
      bvoc_src%chLandUseMetaDataFNms(iFile) = fu_process_filepath(fu_content(pItems(iFile)), &
                                                  & must_exist = .false., superdir = src_data_dir)
    end do
    !
    ! We can have source mask, which would limit the source area
    !
    bvoc_src%chSrcMaskFile = fu_process_filepath(fu_content(nlSetup,'source_mask_file'), &
                                               & must_exist = .false., superdir = src_data_dir)
    bvoc_src%ifMasked = len_trim(bvoc_src%chSrcMaskFile) > 0
    !
    ! Final checking and cleaning
    !
    if(srcOK)then
      bvoc_src%defined = silja_true
    else
      call set_error('Failed to set the bvoc rules','fill_bio_voc_src_from_namelist')
    endif

  end subroutine fill_bio_voc_src_from_namelist


  !**************************************************************************

  subroutine reserve_bio_voc_source(bvoc_src, &     ! Src to initialise
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
    type(silam_bio_voc_source), intent(inout) :: bvoc_src
    integer, intent(in) :: iSrcNbr, iSrcIdNbr

    !
    ! Nullify the basic variables
    !
    bvoc_src%src_nm = ''
    bvoc_src%sector_nm = ''
    !
    ! Main source parameters - enough to identify it in the global information list
    !
    bvoc_src%src_nbr = iSrcNbr
    bvoc_src%id_nbr = iSrcIdNbr
    !
    ! Finally, mark the source as incomplete
    !
    bvoc_src%defined = silja_false

  end subroutine reserve_bio_voc_source


  !*****************************************************************

  subroutine add_input_needs_bio_voc_src(bvoc_src, q_met_dynamic, q_met_static, &
                                                 & q_disp_dynamic, q_disp_static)
    !
    ! Returns input needs for the transformation routines. 
    !
    implicit none

    ! Imported parameters
    type(silam_bio_voc_source), intent(in) :: bvoc_src
    integer, dimension(:), intent(inout) :: q_met_dynamic, q_met_static, &
                                          & q_disp_dynamic, q_disp_static

    ! Local variables
    integer :: iTmp

    !
    ! Add needed dynamic quantities
    !
    select case(bvoc_src%emisMethod)
      case (emisGuenther_updated_v1)

        iTmp = fu_merge_integer_to_array(temperature_2m_flag, q_met_dynamic)
        iTmp = fu_merge_integer_to_array(bvoc_src%swSource, q_met_dynamic)

      case default
        call msg('Unknown emission method', bvoc_src%emisMethod)
        call set_error('Unknown emission type','add_input_needs_bio_voc_src')
        return
    end select

  end subroutine add_input_needs_bio_voc_src


  !**********************************************************************

  subroutine init_emission_bio_voc(bvoc_src, dispersionMarketPtr, start_time)
    !
    ! Actually creates the fields serving the bvoc emission computations
    ! and drops them to the internal dispersion buffer
    !
    implicit none

    ! Imported parameters
    type(silam_bio_voc_source), intent(inout) :: bvoc_src
    type(mini_market_of_stacks), pointer :: dispersionMarketPtr
    type(silja_time), intent(in) :: start_time

    ! Local variables
    type(Tsilam_namelist), pointer :: nlMetaData
    integer :: uLU, uLUMeta, uCity, io_status, nItems, nx, ny, ix,iy, ixDisp, iyDisp, iItem, &
             & iLUType, iMon, ix_start, ix_end, iy_start, iy_end, ixRec, iyRec
    type(Tsilam_nl_item_ptr), dimension(:), pointer :: ptrItems, ptrItemsLUFeature
    character(len=1), dimension(:), pointer :: chTmp1, chTmp2
    character(len=clen) :: strTmp
    real, dimension(:), pointer :: fWork
    real, dimension(:,:), pointer :: fBiomass, fIsopreneFactor, fMonoterpTemprFactor, &
                                   & fMonoterpTemprLightFactor, fOtherVOCFactor
    real, dimension(2,4) :: fCrdTmp
    real :: xDisp, yDisp, fEarthRadTmp, refLon, refLat, dd, xStartShift, yStartShift, &
          & fCosRefLat, fSinRefLat, fDistX, fDistY, sq1, as1, cosAs1_2, sinAs1_2, fTmp, fTmp2, &
          & lon, lat, sinLatStart, sinLatEnd, cosLatEnd, cosLonStart, sinLonStart, cosLonEnd, &
          & sinLonEnd, cosLatStart, lonStart, lonEnd, latStart, latEnd, &
          & sinLat, cosLat, cosLon, sinLon

    integer, parameter :: city_class_USGS = 13

    fWork => fu_work_array(5*500*12)
    if(error)return
    fWork(1:5*500*12) = 0.0
    fBiomass(1:500,1:12)                  => fWork(1:6000)
    fIsopreneFactor(1:500,1:12)           => fWork(6001:12000)
    fMonoterpTemprFactor(1:500,1:12)      => fWork(12001:18000)
    fMonoterpTemprLightFactor(1:500,1:12) => fWork(18001:24000)
    fOtherVOCFactor(1:500,1:12)           => fWork(24001:30000)
    
    allocate(bvoc_src%arEmisFactor(5,nx_dispersion,ny_dispersion,12), stat=io_status)
    if(io_status /= 0)then
      call set_error('Failed to allocate the biogenic-emission metadata array', &
                   & 'init_emission_bio_voc')
      return
    endif
    bvoc_src%arEmisFactor = 0.
    !
    ! Read each land-use metadata file and then the data file
    !
    do iItem = 1, size(bvoc_src%chLandUseMetaDataFNms)
      !
      ! First basic meta data
      !
      uLUMeta = fu_next_free_unit()
      open(unit=uLUMeta,file=bvoc_src%chLandUseMetaDataFNms(iItem), status='OLD', iostat=io_status)
      if(io_status /= 0)then
        call set_error('Failed to open the metadata file:' + bvoc_src%chLandUseMetaDataFNms(iItem), &
                     & 'init_emission_bio_voc')
        return
      endif
      nlMetaData => fu_read_namelist(uLUMeta, .false., 'vegetation_metadata')
      if(error)return
      close(uLUMeta)
      
      nx = fu_content_int(nlMetaData,'dimension_x')
      ny = fu_content_int(nlMetaData,'dimension_y')
      if(fu_content(nlMetaData,'projection') .ne. 'LAMBERT_AZYMUTHAL_EQUAL_AREA')then
        call set_error(fu_connect_strings('Projection must be LAMBERT_AZYMUTHAL_EQUAL_AREA, not:', &
                                        & fu_content(nlMetaData,'projection')), &
                     & 'init_emission_bio_voc')
        return
      endif
      
      fEarthRadTmp = fu_content_real(nlMetaData,'earth_radius')
      if(fEarthRadTmp .eps. real_missing)then
        call msg_warning('No earth radius, take default','init_emission_bio_voc')
        fEarthRadTmp = earth_radius
      endif
      refLon = fu_content_real(nlMetaData,'reference_point_longitude')
      if(refLon .eps. real_missing)then
        call set_error('No reference longitude','init_emission_bio_voc')
        return
      endif
      refLat = fu_content_real(nlMetaData,'reference_point_latitude')
      if(refLat .eps. real_missing)then
        call set_error('No reference latitude','init_emission_bio_voc')
        return
      endif
      dd = fu_content_real(nlMetaData,'grid_cell_size')
      if(dd .eps. real_missing)then
        call set_error('No grid cell size','init_emission_bio_voc')
        return
      endif
      xStartShift = fu_content_real(nlMetaData,'starting_point_x_shift')
      if(xStartShift .eps. real_missing)then
        call set_error('No start-point x-shift','init_emission_bio_voc')
        return
      endif
      yStartShift = fu_content_real(nlMetaData,'starting_point_y_shift')
      if(yStartShift .eps. real_missing)then
        call set_error('No start-point y-shift','init_emission_bio_voc')
        return
      endif

      fCosRefLat = cos(refLat * degrees_to_radians)
      fSinRefLat = sin(refLat * degrees_to_radians)

      !
      ! Fill-in the temporary arrays with the metadata
      !
      nullify(ptrItemsLUFeature)
      call get_items(nlMetaData, 'lu_feature', ptrItemsLUFeature, ix)
      
      if(ix < 1)then
        call report(nlMetaData)
        call set_error('No lu_feature in the namelist','init_emission_bio_voc')
        return
      endif
      
      do iy = 1, ix
        strTmp = fu_content(ptrItemsLUFeature(iy))
        read(unit=strTmp,fmt=*,iostat=io_status) iLUType, iMon
        read(unit=strTmp,fmt=*,iostat=io_status) iLUType, iMon, &
                                               & fBiomass(iLUType,iMon), &
                                               & fIsopreneFactor(iLUType,iMon), &
                                               & fMonoterpTemprFactor(iLUType,iMon), &
                                               & fMonoterpTemprLightFactor(iLUType,iMon), &
                                               & fOtherVOCFactor(iLUType,iMon)
      end do

      call msg('Reading the main land-use map (will take time!)')

      !
      ! Open the main data files
      !
      uLU = fu_next_free_unit()
      call open_binary(unit=uLU, &
                     & file=fu_process_filepath( fu_content(nlMetaData,'vegetation_land_use_file_name'), &
                                  &  superfile=bvoc_src%chLandUseMetaDataFNms(iItem)), &
                     & recl=nx, status='old', action='read')
      if(error)return

      uCity = fu_next_free_unit()
      call open_binary(unit=uCity, &
                     & file=fu_process_filepath(fu_content(nlMetaData,'city_mask_file_name'), &
                                  &  superfile=bvoc_src%chLandUseMetaDataFNms(iItem)), &
                     & recl=nx, status='old', action='read')
      if(error)return

      !
      ! Allocation of the character array for reading the image depends on the size of image
      !
      allocate(chTmp1(nx), chTmp2(nx), stat=io_status)
      if(io_status /= 0)then
        call set_error('Failed to allocate temporary character field','init_emission_bio_voc')
        return
      endif

      !
      ! Task: define the records to read from the land use file.
      ! Step 1: compute the corners of the dispresion area
      !
      fCrdTmp(1,1) = fu_lon_geographical_from_grid(0.5, 0.5, dispersion_grid)
      fCrdTmp(1,2) = fu_lon_geographical_from_grid(0.5, ny_dispersion+0.5, dispersion_grid)
      fCrdTmp(1,3) = fu_lon_geographical_from_grid(nx_dispersion+0.5, 0.5, dispersion_grid)
      fCrdTmp(1,4) = fu_lon_geographical_from_grid(nx_dispersion+0.5, ny_dispersion+0.5, dispersion_grid)
      fCrdTmp(2,1) = fu_lat_geographical_from_grid(0.5, 0.5, dispersion_grid)
      fCrdTmp(2,2) = fu_lat_geographical_from_grid(0.5, ny_dispersion+0.5, dispersion_grid)
      fCrdTmp(2,3) = fu_lat_geographical_from_grid(nx_dispersion+0.5, 0.5, dispersion_grid)
      fCrdTmp(2,4) = fu_lat_geographical_from_grid(nx_dispersion+0.5, ny_dispersion+0.5, dispersion_grid)

      lonStart = min(fCrdTmp(1,1), fCrdTmp(1,2), fCrdTmp(1,3), fCrdTmp(1,4))
      lonEnd = max(fCrdTmp(1,1), fCrdTmp(1,2), fCrdTmp(1,3), fCrdTmp(1,4))
      latStart = min(fCrdTmp(2,1), fCrdTmp(2,2), fCrdTmp(2,3), fCrdTmp(2,4))
      latEnd = max(fCrdTmp(2,1), fCrdTmp(2,2), fCrdTmp(2,3), fCrdTmp(2,4))

      ix_start = nx
      ix_end = 1
      iy_start = ny
      iy_end = 1

      ! Step 2: adjust the corners to the available land use area
      !
      do iy = 1,ny_dispersion,2
        do ix = 1, nx_dispersion, 2
          lon = fu_lon_geographical_from_grid(real(ix), real(iy), dispersion_grid)
          lat = fu_lat_geographical_from_grid(real(ix), real(iy), dispersion_grid)

          sinLat = sin(lat*degrees_to_radians)
          cosLat = cos(lat*degrees_to_radians)
          cosLon = cos((lon-refLon)*degrees_to_radians)
          sinLon = cos((lon-refLon)*degrees_to_radians)

          ixRec = (nint(-xStartShift + &
                     & (sqrt(2/(1+fSinRefLat*sinLat + fCosRefLat*cosLat*cosLon))) * &
                           & (cosLat*sinLon)*fEarthRadTmp)/dd)
          iyRec = nint((yStartShift - &
                     & ((sqrt(2/(1+fSinRefLat*sinLat + fCosRefLat*cosLat*cosLon))) * &
                              & (fCosRefLat*sinLat - fSinRefLat*cosLat*cosLon)*fEarthRadTmp))/dd)
          if(ix_start > ixRec)ix_Start = ixRec
          if(ix_end < ixRec)ix_End = ixRec
          if(iy_start > iyRec)iy_Start = iyRec
          if(iy_end < iyRec)iy_End = iyRec
        enddo
      enddo

      if(ix_start < 1)ix_Start = 1
      if(ix_end > nx)ix_End = nx
      if(iy_start < 1)iy_Start = 1
      if(iy_end > ny)iy_End = ny


      call msg('Reading records from-to:',iy_start, iy_end)
!      call msg('In each record values from-to:',ix_start,ix_end)
!      call msg('Processing whole record')

      !
      ! Read the image file line-by-line
      !
      do iy = iy_start, iy_End  !1, ny
        read (unit=uLU, rec=iy, iostat=io_status) (chTmp1(ix),ix=1,nx)
        if(io_status /= 0)then
          call msg('failed to read the y-record:',iy)
          call msg('IO-status code:',io_status)
          call set_error('failed to read the land use file','init_emission_bio_voc')
          return
        endif
!        read (unit=uCity, rec=iy) (chTmp2(ix),ix=1,nx)
        read (unit=uCity, rec=iy, iostat=io_status) (chTmp2(ix),ix=1,nx)
        if(io_status /= 0)then
          call msg('failed to read the y-record:',iy)
          call msg('IO-status code:',io_status)
          call set_error('failed to read the city file','init_emission_bio_voc')
          return
        endif

        if(mod(iy,1000) == 0)then
          call msg('Done rec:',iy)
        endif

        !
        ! Process each line by cutting only necessary cells and summing them up into the 
        ! dispersion-grid fields
        !
!        do ix = 1, nx
!        do ix = ix_start, ix_end  !1, nx
        do ix = 1, ix_end  !1, nx

          fDistX = (ix-1) * dd + xStartShift
          fDistY = (-iy+1) * dd + yStartShift

          sq1 = sqrt(fDistX * fDistX + fDistY * fDistY)
          if(sq1 < 1.e-3) sq1 = 1.e-3  ! if less than 1 mm from centre point shift it there

          if(sq1 >= 2.0 * fEarthRadTmp)then
!            call msg_warning('Land-use image point is too far from the centre point. Skip', &
!                           & 'compile_land_use_data')
            cycle
          endif

          as1 = asin(sq1 / (2.0 * fEarthRadTmp))
          cosAs1_2 = cos(2. * as1)
          sinAs1_2 = sin(2. * as1)

          fTmp = fDistX * sinAs1_2 / (sq1 * fCosRefLat * cosAs1_2 - fDistY * fSinRefLat * sinAs1_2)

          lon = atan(fDistX * sinAs1_2 / &
                   & (sq1 * fCosRefLat * cosAs1_2 - fDistY * fSinRefLat * sinAs1_2)) * &
                   & radians_to_degrees + refLon

          fTmp2 =cosAs1_2 * fSinRefLat + (fDistY * sinAs1_2 * fCosRefLat) / sq1
          if(fTmp2 >= 1)then
!            call msg_warning('Numerics at N pole: lat problem, set 90','compile_land_use_data')
            lat = 90.
          elseif(fTmp2 <= -1.)then
!            call msg_warning('Numerics at S pole: lat problem, set -90','compile_land_use_data')
            lat = -90.
          else
            lat=asin(fTmp2) * radians_to_degrees
          endif

          !
          ! Having the geographical coordinates of the small cell, have to project it to 
          ! the dispersion grid
          ! Collection goes on only if the small cell is inside the dispersion grid and
          ! metadata arrays show the actual amounts for the land use type picked from the 
          ! iChar array
          !
          call project_point_to_grid(lon, lat, dispersion_grid, xDisp, yDisp)
          if(error)return

          if (xDisp > 0.5 .and. xDisp < nx_dispersion+0.5 .and. &
            & yDisp > 0.5 .and. yDisp < ny_dispersion+0.5)then
            ixDisp = nint(xDisp)
            iyDisp = nint(yDisp)

            if(ichar(chTmp2(ix)) == city_class_USGS)then
              iLUType = 254                             ! change the land use type to city
            else
              iLUType = ichar(chTmp1(ix))
            endif
            if(iLUType == 0)cycle

            !
            ! Note that average emission factor goes with regard to biomass
            !
            do iMon = 1, 12
              fTmp = fBiomass(iLUType,iMon) * dd * dd ! fbiomass * area
              bvoc_src%arEmisFactor(totBiomass,ixDisp,iyDisp,iMon) = &
                                & bvoc_src%arEmisFactor(totBiomass,ixDisp,iyDisp,iMon) + fTmp
              bvoc_src%arEmisFactor(emsIsoprene,ixDisp,iyDisp,iMon) = &
                                & bvoc_src%arEmisFactor(emsIsoprene,ixDisp,iyDisp,iMon) + &
                                & fIsopreneFactor(iLUType,iMon) * fTmp
              bvoc_src%arEmisFactor(emsMonoterpTempr,ixDisp,iyDisp,iMon) = &
                                & bvoc_src%arEmisFactor(emsMonoterpTempr,ixDisp,iyDisp,iMon) + &
                                & fMonoterpTemprFactor(iLUType,iMon) * fTmp
              bvoc_src%arEmisFactor(emsMonoterpTemprLight,ixDisp,iyDisp,iMon) = &
                                & bvoc_src%arEmisFactor(emsMonoterpTemprLight,ixDisp,iyDisp,iMon) + &
                                & fMonoterpTemprLightFactor(iLUType,iMon) * fTmp
              bvoc_src%arEmisFactor(emsOtherVOC,ixDisp,iyDisp,iMon) = &
                                & bvoc_src%arEmisFactor(emsOtherVOC,ixDisp,iyDisp,iMon) + &
                                & fOtherVOCFactor(iLUType,iMon) * fTmp
            end do  ! month
          endif  ! if the small cell is inside the dispersion grid

        end do  ! single line of the image landuse file

      end do  ! iy: reading the image file

      !
      ! The image file is consumed
      !
      deallocate(chTmp1)
      deallocate(chTmp2)
      close(uLU)
      close(uCity)
      call destroy_namelist(nlMetaData)

    end do  ! cycle over the areas with metadata and land-use

    call free_work_array(fWork)

!    write(unit=strTmp,fmt=*) &
!           & (sum(bvoc_src%arEmisFactor(totBiomass,1:nx_dispersion,1:ny_dispersion,iMon)),iMon=1,12)
!    call msg('Done BVOC mask reading. Total biomass for each month: ' + strTmp)

    !
    ! Some dump in dispersion grid
    !
!    uLU = fu_next_free_unit()
!    call msg('Starting dump of the land use (5param,12mon,nx,ny):',nx_dispersion,ny_dispersion)
!    call report(dispersion_grid)
!    open(unit=uLU, file='landuse.dump',access='direct',form='unformatted',recl=nx_dispersion*ny_dispersion)
!    uCity = 1
!    do iMon = 1, 12 
!      write(uLU,rec=uCity)((arEmisFactor(totBiomass,ixDisp,iyDisp,iMon),ixDisp=1,nx_dispersion), &
!                                                                     & iyDisp=1,ny_dispersion)
!      uCity = uCity + 1
!      write(uLU,rec=uCity)((arEmisFactor(emsIsoprene,ixDisp,iyDisp,iMon),ixDisp=1,nx_dispersion), &
!                                                                      & iyDisp=1,ny_dispersion)
!      uCity = uCity + 1
!      write(uLU,rec=uCity)((arEmisFactor(emsMonoterpTempr,ixDisp,iyDisp,iMon),ixDisp=1,nx_dispersion), &
!                                                                           & iyDisp=1,ny_dispersion)
!      uCity = uCity + 1
!      write(uLU,rec=uCity)((arEmisFactor(emsMonoterpTemprLight,ixDisp,iyDisp,iMon),ixDisp=1,nx_dispersion), &
!                                                                                & iyDisp=1,ny_dispersion)
!      uCity = uCity + 1
!      write(uLU,rec=uCity)((arEmisFactor(emsOtherVOC,ixDisp,iyDisp,iMon),ixDisp=1,nx_dispersion), &
!                                                                      & iyDisp=1,ny_dispersion)
!      uCity = uCity + 1
!    end do
!    close(uLU)

  end subroutine init_emission_bio_voc


  !****************************************************************************

  subroutine add_inventory_bio_voc_src(bvoc_src, species_list, nSpecies)
    !
    ! Get the species emitted by this source and make a list of
    ! them. The source must be initialized, of course.
    !
    implicit none
    type(silam_bio_voc_source), intent(in) :: bvoc_src
    type(silam_species), dimension(:), pointer :: species_list
    integer, intent(inout) :: nSpecies

    call addSpecies(species_list, nSpecies, bvoc_src%species, bvoc_src%nSpecies)

  end subroutine add_inventory_bio_voc_src


  !**************************************************************************

  subroutine link_bvoc_src_to_species(species_list, bvoc_src)
    !
    ! Having the cocktails created, we should establish the shortcut links between the 
    ! source descriptors and species list.
    ! Note that cocktail is "anything" and the only requirement is that it has all the species
    ! existing in the source inventory.
    !
    implicit none

    ! Imported parameters
    type(silam_bio_voc_source), intent(inout) :: bvoc_src
    type(silam_species), dimension(:), pointer :: species_list

    !
    ! Linkage is actually just creation of the chemical adaptor
    !
    call create_adaptor(bvoc_src%species, species_list, bvoc_src%adaptor)
    
  end subroutine link_bvoc_src_to_species


  !**********************************************************************
  
  subroutine create_src_cont_grd_bvoc_src(bvoc_src, grid_template, ifVerbose, ifMinimal, ifExtended)
    !
    ! Creates the grid that covers the area with active BVOC emission
    !
    implicit none

    ! Imported parameters
    type(silam_bio_voc_source), intent(in) :: bvoc_src
    type(silja_grid), intent(inout) :: grid_template
    logical, intent(in) :: ifVerbose, ifMinimal
    logical, intent(out) :: ifExtended

    ! So far nothing to do: this source just covers the dispersion grid
    !
    ifExtended = .false.
    return
    
  end subroutine create_src_cont_grd_bvoc_src


  !*****************************************************************

  subroutine project_bvoc_src_second_grd(bvoc_src, grid, vert_disp, vert_proj, iAccuracy)
    !
    ! There is nothing to project in terms of horizontal grid but the vertical has to be made.
    !
    implicit none

    ! Imported parameters
    type(silam_bio_voc_source), intent(inout) :: bvoc_src
    type(silja_grid), intent(in) :: grid
    type(silam_vertical), intent(in) :: vert_disp, vert_proj
    integer, intent(in) :: iAccuracy

    ! Local variables
    real, dimension(max_levels) :: levFr
   integer :: i, ilev
    type(silam_vertical) :: vertTmp

    if(iAccuracy < 1 .or. iAccuracy > 10)then
      call msg('Accuracy switch must be from 1 to 10, not:',iAccuracy)
      call msg_warning('Accuracy switch must be from 1 to 10','project_bvoc_src_second_grd')
!      return
    endif

    if(len_trim(bvoc_src%sector_nm) > 0)then
      call msg('Re-projecting bio-voc source:' + bvoc_src%src_nm +'_' + bvoc_src%sector_nm)
    else
      call msg('Re-projecting bio-voc source:' + bvoc_src%src_nm)
    endif

    !
    ! Now re-project the vertical grid of the source to the given vertical
    ! Since vertical distributions are very poorly known and crude, we do not need
    ! any precise meteo-dependent projection, crude will do the job
    !
    bvoc_src%vertLevsDispVert = vert_disp
    allocate(bvoc_src%levFractDispVert(fu_NbrOfLevels(vert_disp)), &
           & bvoc_src%fzDisp(fu_NbrOfLevels(vert_disp)), stat=i)
    if(i /= 0)then
      call set_error('Failed to allocate dispersion-vertical level fractions','project_bvoc_src_second_grd')
      return
    endif
    bvoc_src%levFractDispVert(:) = 0.0
    bvoc_src%fzDisp(:) = 0.0

    !
    ! Create the sea salt vertical, which can depend on the emission method
    !
    select case(bvoc_src%emisMethod)
      
      case(emisGuenther_updated_v1)
        !
        ! Guenther presently goes from 10 to 50 m.
        !
        call set_vertical(fu_set_layer_between_two(layer_btw_2_height, 10.0, 25.0), vertTmp)
        if(error)return
        call add_level(vertTmp, fu_set_layer_between_two(layer_btw_2_height, 25.0, 50.0))
        if(error)return
        levFr(1) = 0.6
        levFr(2) = 0.4
      case default
        call set_error('Unknown emission method:'+fu_str(bvoc_src%emisMethod),'project_bvoc_src_second_grd')
        return
    end select

    call reproject_verticals(vertTmp, levFr, &                    ! vertical from, fractions from
                           & vert_proj, bvoc_src%levFractDispVert, &   ! vertical to, fractions to
                           & bvoc_src%fzDisp, bvoc_src%nLevsDispVert, & ! mass centres, number of non-zero levels
                           & ifMassCentreInRelUnit=.true.)
    call set_missing(vertTmp, .false.)

  end subroutine project_bvoc_src_second_grd


  !**************************************************************************

  subroutine compute_emission_for_bio_voc(bvoc_src, &
                                        & met_buf, disp_buf, & 
                                        & now, &      ! current time
                                        & timestep, & ! model time step
                                        & pHorizInterpMet2DispStruct, ifHorizInterp, &
                                        & ifSpeciesMoment, &
                                        & emisMap, mapCoordX, mapCoordY, mapCoordZ, fMassInjected)    ! Output
    !
    ! Compute the bio VOC emission and store in emisMap with indices
    ! ind_isop (isoprene) and ind_mterp (monoterpentine). If index is
    ! int_missing, the species is omitted.
    ! 
    ! The computation follows the main formula:
    !
    ! E_voc = E * D * C * timestep
    ! here:
    ! E is the emission factor = E(landuse,month) [kg / sec kg-dry-weight]
    ! D is the foliar biomass [kg/m2] = D(landuse,month)
    ! C is the depenence on temperature and light [unitless]
    ! timestep is tyhe emission timestep [sec]
    !
    ! This routine is to be called at each model time step but not inside the particle cycle,
    ! therefore its efficiency is of moderate importance.
    !
    implicit none
    ! Arguments
    type(silam_bio_voc_source), intent(in) :: bvoc_src
    type(Tfield_buffer), pointer :: met_buf, disp_buf
    type(silja_time), intent(in) :: now
    type(silja_interval), intent(in) :: timestep
    type(tHorizInterpStruct), intent(in) :: pHorizInterpMet2DispStruct
    logical, intent(in) :: ifHorizInterp, ifSpeciesMoment
    type(Tmass_map), intent(inout), target :: emisMap, mapCoordX, mapCoordY, mapCoordZ
    real(r8k), dimension(:), intent(inout) :: fMassInjected
    
    ! Local variables
    integer :: iMeteo, iDisp, ix, iy, iLev, ixMeteo, iyMeteo, iStat, iMode, iMon 
    integer :: indT2m, indSWRorCC, indP
    type(Tmass_map), pointer :: eM
    real :: fTmp, timestep_sec, T2m, Psfc, SWR, exp_T2m_Ts, exp_T2m_Tm, exp_T2m_Ts_simple, &
          & sqrt_SWR, e_iso, e_mts, e_mtl, e_ter, e_voc, e_iso_total, e_mte_total, &
          & fCellTotal, factor_isop
    real :: lat, lon, cosZen, att
    type(Tfield_buffer), pointer ::  mb, db

    type(silam_sp) :: sp
    character (len=*), parameter :: sub_name = "compute_emission_for_bio_voc"

    ! Local parameters after NATAIR report, also the Guenther notations
    !
    real, parameter :: Ts = 303.
    real, parameter :: Tm = 314.
    real, parameter :: alpha_SWR = 2.7e-3 * 0.45 *4.6 ! PAR=0.45*SWR(W/m2)*4.6(microphotons/s per W)
    real, parameter :: beta_T = 0.09   ! [1/K]
    real, parameter :: C_L_1 = 1.066
    real, parameter :: C_T_1 = 95000   ! [J/mol]
    real, parameter :: C_T_2 = 230000  ! [J/mol]

    !
    ! First, set the output pointer to the locally stored cocktail_map of emission
    ! intensity
    !
    eM => emisMap
    mb => met_buf
    db => disp_buf

    if(.not. emisMap%gridTemplate == dispersion_grid)then
      call msg('Emission map grid:')
      call report(emisMap%gridTemplate)
      call msg('Dispersion grid:')
      call report(dispersion_grid)
      call set_error('Grid of bvoc emission map differs from dispersion grid', &
                   & 'compute_emission_for_bio_voc')
      return
    endif

    !
    ! Some speed-up first: set the whole map set to be void
    ! No - you don't own the emission map!
    !emisMap%ifColumnValid = .false.
    !emisMap%ifGridValid = .false.

    timestep_sec = abs(fu_sec(timestep))
    if(error)return

    e_iso_total = 0.0
    e_mte_total = 0.0

    iMon = fu_mon(now)

    !--------------------------------------------------------------------
    !
    ! Computation of the emission depends on method of emission computation
    !
    select case(bvoc_src%emisMethod)
      
      case(emisGuenther_updated_v1)
        !
        ! Basic case: updated Guenther et al approach
        !
        indT2m = fu_index(met_buf, temperature_2m_flag) 
        if( fu_fails(indT2m > 0,"Failed to find T2m field", sub_name )) return


        indSWRorCC = fu_index(met_buf, bvoc_src%swSource) 
        if ( fu_fails( indSWRorCC > 0, &
              & 'Failed to find '//trim(fu_quantity_string(bvoc_src%swSource))//' field',&
              & sub_name)) return

        if (bvoc_src%swSource == total_cloud_cover_flag) then
          indP = fu_index(met_buf, ground_pressure_flag) ! Must be there anyway
          if( fu_fails(indT2m > 0,"Failed to find Ps field", sub_name )) return
        endif
        if (error) return

        ! 
        ! Computations go grid cell by grid cell, excluding frozen areas and 100% open water
        ! 
        do iy = 1, ny_dispersion
          do ix = 1, nx_dispersion
            iMeteo =  fu_grid_index(nx_meteo, ix, iy, pHorizInterpMet2DispStruct)

            T2m = met_buf%p2d(indT2m)%present%ptr(iMeteo)
            if (T2m > 273.15 .and. bvoc_src%arEmisFactor(totBiomass, ix, iy, imon) > 0.0) then
!!$            if(T2m > 273.15)then
!!$              !
!!$              ! Not frozen: check for biomass
!!$              !
!!$              if(factors(totBiomass,ix,iy,iMon) > 0.0)then
              !
              ! Some green stuff exists. Compute emission! Four species: 
              !   isoprene
              !   temperature- dependent monoterpenes
              !   temperature+light-dependent monoterpenes
              !   other VOCs
              !
              SWR = met_buf%p2d(indSWRorCC)%present%ptr(iMeteo) !Can be cloud_cover
              SWR = max( 0., SWR)
              if (bvoc_src%swSource == total_cloud_cover_flag) then
                Psfc = met_buf%p2d(indP)%present%ptr(iMeteo)
                ! Get SWR from coud cover and pressure
                ! same way as in  g_sto
                lat = fu_lat_geographical_from_grid(real(ix), real(iy), dispersion_grid)
                lon = fu_lon_geographical_from_grid(real(ix), real(iy), dispersion_grid)
                cosZen = max(fu_solar_zenith_angle_cos(lon, lat, now),0.)
                if (cosZen > 1e-5) then ! Avoid arithmetic error
                 att = 1.0 - 0.75*SWR**3.4  !(source: Kasten & Czeplak (1980)) 
                 SWR  = att * Ashrae%a * exp(- Ashrae%b * (Psfc/1e5)/CosZEN) * CosZEN
                else
                 SWR = 0.
               endif
              endif


              exp_T2m_Ts = exp(C_T_1 * (T2m-Ts) / (gas_constant_uni * Ts * T2m))
              exp_T2m_Ts_simple = exp(beta_T * (T2m - Ts))
              exp_T2m_Tm = exp(C_T_2 * (T2m-Tm) / (gas_constant_uni * Ts * T2m)) ! 303 correct!

              sqrt_SWR = sqrt(1 + (alpha_SWR * SWR) * (alpha_SWR * SWR))
              
              ! ISOPRENE
              if (bvoc_src%ind_emis_isoprene /= int_missing) then
                e_iso = bvoc_src%arEmisFactor(emsIsoprene,ix,iy,iMon) * &
                      & C_L_1 * alpha_SWR * SWR / sqrt_SWR * exp_T2m_Ts / (1 + exp_T2m_Tm)
                ! Unit conversion from kg C to moles of C5H8 
                e_iso = e_iso * 16.666666666666667
                e_iso = e_iso * bvoc_src%factor_emis_isoprene
                e_iso = e_iso * timestep_sec
                if(e_iso < 0)then
                  call msg("met_buf%p2d(indSWRorCC)%present%ptr(iMeteo)",met_buf%p2d(indSWRorCC)%present%ptr(iMeteo) )
                  !0.37350E+02  0.44750E+02  0.11783E-01  -.60270E-02 0.10230E+06  0.12007E+04  0.14170E+00

                  call msg("lon,lat,cosZen,att,pSfc,Ashrae%a,Ashrae%b", (/lon,lat,cosZen,att,pSfc,Ashrae%a,Ashrae%b/))
                  !0.10660E+01  0.55890E-02  -.38725E-06  0.10000E+01 0.64622E-01  0.37658E-04
                  call msg("C_L_1, alpha_SWR,SWR,sqrt_SWR, exp_T2m_Ts, exp_T2m_Tm", (/C_L_1, alpha_SWR,SWR,sqrt_SWR, exp_T2m_Ts, exp_T2m_Tm/))
                  call set_error('Negative isoprene emission','compute_emission_for_bio_voc')
                  return
                endif
                e_iso_total = e_iso_total + e_iso
              else
                e_iso = 0.0
              endif
              
              ! MONOTERPENE
              if (bvoc_src%ind_emis_monoterpene /= int_missing) then
                e_mts = bvoc_src%arEmisFactor(emsMonoterpTempr,ix,iy,iMon) * exp_T2m_Ts_simple

                e_mtl = bvoc_src%arEmisFactor(emsMonoterpTemprLight,ix,iy,iMon) * &
                       & C_L_1 * alpha_SWR * SWR / sqrt_SWR * exp_T2m_Ts / (1 + exp_T2m_Tm)
                ! Unit conversion from kg C to moles of C5H8_2
                e_mtl = e_mtl * 8.3333333333333333
                e_mts = e_mts * 8.3333333333333333
                e_mtl = e_mtl * bvoc_src%factor_emis_monoterpene
                e_mts = e_mts * bvoc_src%factor_emis_monoterpene
                e_mtl = e_mtl * timestep_sec
                e_mts = e_mts * timestep_sec
                if(e_mts + e_mtl < 0)then
                  call set_error('Negative monoterpene emission','compute_emission_for_bio_voc')
                  return
                endif
                e_mte_total = e_mte_total + (e_mtl + e_mts) 
              else
                e_mtl = 0.0
                e_mts = 0.0              
              endif
              !
              ! Fill-in the emission map, not forgetting the number to mass convertion where needed
              !
              do iLev = 1, bvoc_src%nLevsDispVert
                !
                ! First do the check for the overlap: speed-up
                !
                if(bvoc_src%levFractDispVert(iLev) == 0.0)cycle  ! nothing for this dispersion layer

                fCellTotal = 0.0
                
                ! ISOPRENE
                if (bvoc_src%ind_emis_isoprene /= int_missing) then
                  if(e_iso > 0.)then
                    emisMap%arM(bvoc_src%adaptor%iSp(bvoc_src%ind_emis_isoprene), &
                              & bvoc_src%id_nbr, iLev,ix,iy) = &
                         & emisMap%arM(bvoc_src%adaptor%iSp(bvoc_src%ind_emis_isoprene), &
                                     & bvoc_src%id_nbr, iLev,ix,iy) + &
                         & e_iso * bvoc_src%levFractDispVert(iLev)
                    emisMap%ifColumnValid(bvoc_src%id_nbr, ix, iy) = .true.
                    emisMap%ifGridValid(iLev, bvoc_src%id_nbr) = .true.
                   
                    fCellTotal = fCellTotal + e_iso * bvoc_src%levFractDispVert(iLev)

                    if (ifSpeciesMoment) then
!                      mapCoordX%arm(iSpeciesEmis,a_src%id_nbr,ilev, ix,iy) = &
!                                  & mapCoordX%arm(iSpeciesEmis, a_src%id_nbr, ilev, ix,iy) + &
!                                  & fMassTmp * 0.0
!                      mapCoordY%arm(iSpeciesEmis, a_src%id_nbr, ilev, ix,iy) = &
!                                  & mapCoordY%arm(iSpeciesEmis, a_src%id_nbr, ilev, ix,iy) + &
!                                  & fMassTmp * 0.0
                      mapCoordZ%arm(bvoc_src%adaptor%iSp(bvoc_src%ind_emis_isoprene), &
                                  & bvoc_src%id_nbr, ilev, ix,iy) = &
                           & mapCoordZ%arm(bvoc_src%adaptor%iSp(bvoc_src%ind_emis_isoprene), &
                                         & bvoc_src%id_nbr, ilev, ix,iy) + &
                           & e_iso * bvoc_src%levFractDispVert(iLev) * bvoc_src%fzDisp(iLev)
                    endif
                  end if
                end if ! isoprene emission requested
                
                ! MONOTERPENE
                if (bvoc_src%ind_emis_monoterpene /= int_missing) then
                  if(e_mtl + e_mts > 0.)then
                    emisMap%arM(bvoc_src%adaptor%iSp(bvoc_src%ind_emis_monoterpene), &
                              & bvoc_src%id_nbr, iLev,ix,iy) = &
                          & emisMap%arM(bvoc_src%adaptor%iSp(bvoc_src%ind_emis_monoterpene), &
                                      & bvoc_src%id_nbr, iLev,ix,iy) + &
                          & (e_mts + e_mtl) * bvoc_src%levFractDispVert(iLev)
                    emisMap%ifColumnValid(bvoc_src%id_nbr, ix, iy) = .true.
                    emisMap%ifGridValid(iLev, bvoc_src%id_nbr) = .true.

                    fCellTotal = fCellTotal + (e_mts + e_mtl) * bvoc_src%levFractDispVert(iLev)

                    if (ifSpeciesMoment) then
!                      mapCoordX%arm(iSpeciesEmis,a_src%id_nbr,ilev, ix,iy) = &
!                                  & mapCoordX%arm(iSpeciesEmis, a_src%id_nbr, ilev, ix,iy) + &
!                                  & fMassTmp * 0.0
!                      mapCoordY%arm(iSpeciesEmis, a_src%id_nbr, ilev, ix,iy) = &
!                                  & mapCoordY%arm(iSpeciesEmis, a_src%id_nbr, ilev, ix,iy) + &
!                                  & fMassTmp * 0.0
                      mapCoordZ%arm(bvoc_src%adaptor%iSp(bvoc_src%ind_emis_monoterpene), &
                                  & bvoc_src%id_nbr, ilev, ix,iy) = &
                           & mapCoordZ%arm(bvoc_src%adaptor%iSp(bvoc_src%ind_emis_monoterpene), &
                                         & bvoc_src%id_nbr, ilev, ix,iy) + &
                           & (e_mts + e_mtl) * bvoc_src%levFractDispVert(iLev) * bvoc_src%fzDisp(iLev)
                    endif
                  end if
                end if

                if (.not. ifSpeciesMoment) then
!                  mapCoordX%arM(1,bvoc_src%id_nbr, iLev, ix, iy) = 0.0 * fCellTotal + &
!                                                      & mapCoordX%arM(1,bvoc_src%id_nbr, iLev, ix, iy)
!                  mapCoordY%arM(1,bvoc_src%id_nbr, iLev, ix, iy) = 0.0 * fCellTotal + &
!                                                      & mapCoordY%arM(1,bvoc_src%id_nbr, iLev, ix, iy)
                  mapCoordZ%arM(iLev, bvoc_src%id_nbr, iLev, ix, iy) = &
                                                      & bvoc_src%fzDisp(iLev) * fCellTotal + &
                                                      & mapCoordZ%arM(1,bvoc_src%id_nbr, iLev, ix, iy)
                end if

              end do  ! iLev

            end if  ! if T2m > 273.15 and biomass in gridpoint

        end do  ! ix
      end do  ! iy

    case default
        call msg('Unknown method for biogenic VOC emission:', bvoc_src%emisMethod)
        call set_error('Unknown method for biogenic VOC emission','compute_emission_for_bio_voc')
        return

    end select  ! emisMethod

    if (bvoc_src%ind_emis_isoprene /= int_missing) then
    
      call msg('Isoprene emission total:', e_iso_total)
      fMassInjected(bvoc_src%adaptor%iSp(bvoc_src%ind_emis_isoprene)) = &
              & fMassInjected(bvoc_src%adaptor%iSp(bvoc_src%ind_emis_isoprene)) + &
              & e_iso_total
    endif
    if (bvoc_src%ind_emis_monoterpene /= int_missing) then
      call msg('Monoterpene emission total:', e_mte_total)
      fMassInjected(bvoc_src%adaptor%iSp(bvoc_src%ind_emis_monoterpene)) = &
            & fMassInjected(bvoc_src%adaptor%iSp(bvoc_src%ind_emis_monoterpene)) + &
            & e_mte_total
    endif
    
    

  end subroutine compute_emission_for_bio_voc


  !********************************************************************************************

  logical function fu_bio_voc_emis_owned_quantity(bvoc_src, quantity)
    !
    ! Checsk whether the specific quantity is handled by the emission sources exclusively
    ! Such quantities may be related to emission species somewhere else
    !
    implicit none

    ! Imported parameters
    type(silam_bio_voc_source), intent(in) :: bvoc_src
    integer, intent(in) :: quantity
    !
    ! The bio-voc source does not have own quantities yet
    !
    select case(quantity)
      case default
        fu_bio_voc_emis_owned_quantity = .false.
    end select

  end function fu_bio_voc_emis_owned_quantity


  !*****************************************************************

  integer function fu_source_id_nbr_of_bvoc_src(bvoc_src)
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
    type(silam_bio_voc_source), intent(in) :: bvoc_src

    ! Stupidity check
    if(.not.(bvoc_src%defined == silja_false))then
      call set_error('Undefined source given','fu_source_id_nbr_of_bvoc_src')
      return
    endif
    fu_source_id_nbr_of_bvoc_src = bvoc_src%id_nbr

  end function fu_source_id_nbr_of_bvoc_src


  !*************************************************************************

  integer function fu_source_nbr_of_bvoc_src(bvoc_src)
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
    type(silam_bio_voc_source), intent(in) :: bvoc_src

    ! Stupidity check
    if(.not.(bvoc_src%defined == silja_false))then
      fu_source_nbr_of_bvoc_src = bvoc_src%src_nbr
    else
      fu_source_nbr_of_bvoc_src = int_missing
      call set_error('Undefined source given','fu_source_nbr_of_bvoc_src')
      return
    endif

  end function fu_source_nbr_of_bvoc_src


  !*************************************************************************

  subroutine typical_species_cnc_bvoc(bvoc_src, &
                                    & species, nSpecies, arConc)  ! output
    !
    ! Guesses a typical level of concentration.
    !
    implicit none

    ! Imported parameters
    type(silam_bio_voc_source), intent(in), target :: bvoc_src
    type(silam_species), dimension(:), pointer :: species
    integer, intent(out) :: nSpecies
    real, dimension(:), pointer :: arConc

    species => bvoc_src%species
    nSpecies = bvoc_src%nSpecies

    arConc(1:nSpecies) = 1.e-9

  end subroutine typical_species_cnc_bvoc


  !*************************************************************************

  function fu_bio_voc_source_name(bvoc_src)result(chNm)
    implicit none
    type(silam_bio_voc_source), intent(in) :: bvoc_src
    character(len=clen) :: chNm
    chNm = bvoc_src%src_nm
  end function fu_bio_voc_source_name


  !*************************************************************************

  subroutine report_bio_voc_src(bvoc_src)
    implicit none
    type(silam_bio_voc_source), intent(in) :: bvoc_src
    integer :: iSpecies
    call msg('bvoc source'+bvoc_src%src_nm)
    do iSpecies = 1, bvoc_src%nSpecies
      call report(bvoc_src%species(iSpecies))
    end do
  end subroutine report_bio_voc_src


END MODULE source_terms_bio_voc



