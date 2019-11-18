MODULE source_term_fires
  !
  ! This module contains description of the fires emission.
  !
  ! Emission is computed for the set of species requested in the source header ini file.
  ! This is an inventory-type emission but with pronounced meteorological influence.
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
  ! PUBLIC routines of fire source
  !
  public fill_fire_src_from_namelist
  public reserve_fire_source
  public init_emission_fire
  public create_src_cont_grd_fire_src
  public force_source_into_grid
  public source_2_second_grid
  public add_inventory_fire_src
  public add_input_needs
  public link_source_to_species
  public prepare_inject_fire_src
  public inject_emission_euler_fire_src
  public fu_name
  public fu_sector
  public fu_source_nbr
  public fu_source_id_nbr
  public total_amt_species_unit
  public fu_start_time
  public fu_end_time
  public fu_duration
  public fu_n_fires
  public report
  public defined
  !
  ! Private routines of the fire source
  !
  private get_FRP_dataset
  private fu_get_fire_metadata
  private clean_fire_metadata
  private add_input_needs_fire_src
  private force_fire_source_into_grid
  private project_fire_src_second_grd
  private remove_fire
  private link_fire_src_to_species
  private fu_source_id_nbr_of_fire_src
  private fu_source_nbr_of_fire_src
  private fu_fires_source_name
  private fu_fires_sector_name
  private tot_amt_species_unit_fire_src
  private fu_start_time_fire_src
  private fu_end_time_fire_src
  private fu_duration_fire_src
  private fu_emission_weighted_flam_smld
  private fires_flux4mode
  private report_fire_src
  private report_metadata
  private defined_fire_src

  !
  ! Private subs of fires source
  !
  interface add_input_needs
    module procedure add_input_needs_fire_src
  end interface

  interface source_2_second_grid
    module procedure project_fire_src_second_grd
  end interface

  interface force_source_into_grid
    module procedure force_fire_source_into_grid
  end interface

  interface link_source_to_species
    module procedure link_fire_src_to_species
  end interface

  interface fu_source_nbr
    module procedure fu_source_nbr_of_fire_src
  end interface

  interface fu_source_id_nbr
    module procedure fu_source_id_nbr_of_fire_src
  end interface

  interface fu_name
    module procedure fu_fires_source_name
  end interface

  interface fu_sector
    module procedure fu_fires_sector_name
  end interface

  interface total_amt_species_unit
    module procedure tot_amt_species_unit_fire_src
  end interface

  interface fu_start_time
    module procedure fu_start_time_fire_src
  end interface

  interface fu_end_time
    module procedure fu_end_time_fire_src
  end interface

  interface fu_duration
    module procedure fu_duration_fire_src
  end interface

  interface report
    module procedure report_fire_src
    module procedure report_metadata
  end interface

  interface defined
    module procedure defined_fire_src
  end interface

  !
  ! The fire metadata. That creature does not have to be one for all sources but
  ! there should be not too many of them. Each of them has a huge land-use map.
  ! A mechanism will be made for several sources to share single metadata file
  !
  ! We allow two types of fire: flaming and smouldering. They emit THE SAME SPECIES,
  ! except for mass-mean diameters in otherwise identical bins - but modes are compared
  ! without this parameter (see chemical setup).
  !
  type Tsilam_fire_metadata
    private
    character(len=fnlen) :: chIniFNm
    integer :: nLU_types, iSpectrumType
    integer :: nSpecies
    type(silja_grid) :: grid_metadata
    character(len=SubstNmLen), dimension(:), pointer :: chLU_names      ! (nLandUseTypes)
    type(silam_species), dimension(:), pointer :: species_flaming=>null(), species_smouldering=>null() ! (nSpecies)
    type(chemical_adaptor) :: adaptor              ! (flaming, smouldering)
    real, dimension(:,:,:), pointer :: pEmsFactors               ! (flam/smld, nLandUseTypes, nSpecies_max)
    real, dimension(:,:), pointer :: pDiurnalVarTot, pDiurnalVarPerFire ! (24, nLandUseTypes)
    real, dimension(:,:), pointer :: flaming_fraction_roughly           ! (iLU, iSp)
    integer*1, dimension(:,:), pointer :: indLU            ! main map of land-use types
    logical, dimension(:), pointer :: ifNumberFlux => null()
    type(silja_logical) :: defined
  end type Tsilam_fire_metadata
  type Tsilam_fire_metadata_ptr
    type(Tsilam_fire_metadata), pointer :: ptr
  end type Tsilam_fire_metadata_ptr
  !
  ! FRP data set. The second component of the fire source term
  ! Apart from basic FRP-based dataset, there can be also analytical description
  ! of the fire. It can be read from ini file instead of FRP dataset or deduced
  ! from the consumed FRP dataset.
  ! Types of fire description are coded via iDescriptorType
  !
  type TFRP_dataset
    private
    integer :: nFires, nMaxObs, indFRPDaily_tot, indFRPDaily_perFire, iDescriptionType
    type(silja_time), dimension(:), allocatable :: day     ! (nFires)
    integer, dimension(:), allocatable :: indLU   ! (nFires)
    real, dimension(:), allocatable :: pLonGeo, fX, fY
    real, dimension(:,:), allocatable :: dx, dy, pFRP, pTA, pT4, pT4b, pT11, pT11b, &
                                   & pMCE, pArea, pHour ! (nFires, nMaxObs)
    real, dimension(:), allocatable :: FPRmax, hStart, hEnd, lon, lat   ! nFires
  end type TFRP_dataset

  !
  ! The fires source term. The idea is to have as many things computed at the
  ! moment of reading as possible. Meteo-dependent thing so far is only injection profile.
  ! Contains general descriptions, pointer to the metadata structure, and
  ! a set of pointers to the FRP datasets.
  !
  TYPE Tsilam_fire_source
    PRIVATE
    CHARACTER(len=clen) :: src_nm, sector_nm  ! Name of the area source and sector
    integer :: src_nbr, id_nbr                ! A source and id numbers in a WHOLE source list
    integer :: nFRPdatasets, nSpecies
    real, dimension(:,:), pointer :: fluxPerModeNbr, fluxPerModeVol ! (nModes, flame/smoulder)
    type(silam_species), dimension(:), pointer :: species_flaming=>null(), species_smouldering=>null()
    logical, dimension(:), pointer :: ifNumberFlux
    type(silja_time), dimension(:), pointer :: first_day, last_day  ! (nFRPsubsets)
    type(silja_grid) :: grid
    type(Tsilam_fire_metadata), pointer :: pFMD
    type(TFRP_dataset), dimension(:), allocatable :: pFRPset  ! (nFRPsubsets)
    logical :: ifGeoCoord             ! can also be grid indices
    type(silja_logical) :: defined
  END TYPE Tsilam_fire_source

  type fire_src_ptr
    type(Tsilam_fire_source) :: fire_src
  end type fire_src_ptr
  public fire_src_ptr

  !
  ! The fire metadata internal collection
  !
  type(Tsilam_fire_metadata_ptr), dimension(10), private, save :: arFMD_glob
  integer, private, save :: nFMD_glob = 0

  !
  ! Some uncertainty for the FRP has to be given. For the time being, let's take
  ! something from a blue sky.
  !
  real, parameter, private :: sigma_FRP_abs = 1.0e7 ! 10 MW...
  real, parameter, private :: sigma_FRP_rel = 0.2 ! plus 20%
  !
  ! Fire plume rise can be obtained either from one- or two-step procedures
  !
  logical, private, parameter :: ifOneStepHeightProcedure = .false.

  !
  !
  !
  logical, private, save :: ifUseHybrid = .false.
  !
  ! The fire plume is split into stem and hat, with some 10% of mass fraction released into stem
  !
  integer, private, parameter :: iPlumeStem = 1, iPlumeHat = 2  ! fire plume is considered to be a mushroom
  real, private :: fStemMassFraction

  !
  ! The fire aerosol is emitted as the following mass spectrum types. Flaming and 
  ! smouldering emit different spectra (same modes, different fractions)
  !
  integer, private, parameter :: lognorm_3_modes_fires_flag = 6120
  integer, private, parameter :: lognorm_3_modes_fires_num_flag = 6121

  integer, parameter :: flaming = 1, smouldering = 2
  character(len=*), dimension(2), parameter :: flamsmold = (/"flaming    ","smouldering"/)
  
  !
  ! Meteodata that will be needed for plume rise
  !
  type(field_4d_data_ptr), private, pointer, save :: fldBVf
  type(field_4d_data_ptr), private, pointer, save :: fldHeight
  type(field_2d_data_ptr), private, save, pointer :: fldAblHeight
  type(field_2d_data_ptr), private, pointer, save :: fldSrfPressure


CONTAINS


  !*********************************************************************

  subroutine fill_fire_src_from_namelist(nlSetup, fs, expected_species, dirname)
    !
    ! Initializes the fires source term.
    ! The parameters are read from the given ini file and stored to the returned 
    ! source term. After the stuff has been read, the main FluxPerMode is to be generated.
    ! This configuration allows for several fires sources
    !
    implicit none

    ! Imported parameters
    type(Tsilam_namelist), intent(in) :: nlSetup
    type(Tsilam_fire_source), intent(inout) :: fs
    type(silam_species), dimension(:), intent(in), allocatable :: expected_species
    character(len=*), intent(in) :: dirname

    ! Local variables
    integer :: iTmp, iUnit
    type(Tsilam_nl_item_ptr), dimension(:), pointer :: pItems
    character (len=fnlen) :: metadata_file
    !
    ! Names
    !
    fs%src_nm = fu_content(nlSetup,'source_name')
    fs%sector_nm = fu_content(nlSetup,'source_sector_name')
    fs%defined = silja_false
    fs%grid = geo_global_grid
    !
    ! Emission factors, cocktails, etc are in the metadata file.
    ! In theory, sources can have different sets - for instance, land-use can be 
    ! made continent-specific. But there is no reason to store same metadata more than once.
    ! The module has a unified metadata repository.
    !
    fs%pFMD => fu_get_fire_metadata(fu_extend_grads_hat_dir( &
                                            &  fu_content(nlSetup,'fire_metadata_file'), dirname), &
                                  & expected_species, nlSetup)

    if(error .or. .not. associated(fs%pFMD))return
    !
    ! Fire data are groupped to FRPset datasets, mostly one-set per day. These are sitting
    ! under frp_dataset namelist items. Get them, allocate the poitners and then
    ! read them one by one.
    !
    nullify(pItems)
    call get_items(nlSetup,'frp_dataset',pItems,fs%nFRPdatasets)
    if(error .or. fs%nFRPdatasets < 1 .or. fs%nFRPdatasets > 100000)then
      call set_error('Failed to get FRP datasets for source:'+fs%src_nm,'fill_fire_src_from_namelist')
      return
    endif
    allocate(fs%first_day(fs%nFRPdatasets), fs%last_day(fs%nFRPdatasets), &
           & fs%pFRPset(fs%nFRPdatasets), stat = iTmp)
    if(iTmp /= 0)then
      call set_error('Failed FRP sets allocation.Source' + fs%src_nm + ', N_sets=' + fu_str(fs%nFRPdatasets), &
                   & 'fill_fire_src_from_namelist')
      return
    endif
    
#ifdef VOIMA_GNU_BUG
    !$OMP PARALLEL if (.False.) default(none), &
#else    
    !$OMP PARALLEL if (.True.) default(none), &
#endif
    !$OMP & private(iTmp, iUnit),  &
    !$OMP & shared (pItems, dirname, fs, error)

    iUnit = fu_next_free_unit()

    !$OMP DO
    do iTmp = 1, fs%nFRPdatasets
      if (error) cycle
      call get_FRP_dataset(fu_process_filepath(fu_extend_grads_hat_dir(fu_content(pItems(iTmp)),dirname)), &
                         & fs%first_day(iTmp), fs%last_day(iTmp), &
                         & fs%pFMD, fs%pFRPset(iTmp), iUnit)
      if(error) call set_error('Failed FRP dataset in source:' + fs%src_nm + '_' + fs%sector_nm, &
                     & 'fill_fire_src_from_namelist')

      if(mod(iTmp,20)==0) call msg('Now reading FRP dataset:' + fs%src_nm + '_' + fs%sector_nm + ', i,n:', &
                                   & iTmp, fs%nFRPdatasets)
    end do  ! FRPdatasets
    !$OMP END DO
    !$OMP END PARALLEL
    
    fs%ifGeoCoord = .true.
    fs%defined = silja_true
    
    deallocate(pitems)

    call report(fs)
    
  end subroutine fill_fire_src_from_namelist


  !***********************************************************************************

  subroutine get_FRP_dataset(chFNm, first_day, last_day, pFMD, pSet, uIn)
    !
    ! Reads one FRP set, usually a daily collection of FRP data
    !
    implicit none
      
    ! Imported parameters
    character(len=*), intent(in) :: chFNm  ! file with FPR set
    type(silja_time), intent(out) :: first_day, last_day  ! range when the set is active
    type(Tsilam_fire_metadata), intent(in) :: pFMD
    type(TFRP_dataset), intent(out) :: pSet     ! main structure to fill-in
    integer, intent(in) :: uIn ! Recycle unit to use

    ! Local variables
    character(len=fnlen) :: lineTmp
    character(len=clen) :: strTmp
    character(len=unitNmLen) :: chSizeUnit, chFRPUnit, chSizeUnit1, chFRPUnit1, strTmp1
    integer :: iItem, nItems, iStat, iStat1, iLU, yr,mon,day,hr,mn, nx, ny, iFire, iTmp, iObs, iCount
    integer :: nFires, nMaxObs
    integer :: iLine
    real :: fSizeConv, fPowerConv !Conversion factors
    real :: fDx, fDy, fpT4, fpT4b, fpT11, fPT11b, fpTA, fFRP, fpMCE, fParea
    real :: sec, fTmp, fX, fY, fLonTmp, fLatTmp, variance_sum, variance_inv, fHour, fWeightPast
    type(silja_time) :: timeTmp
    logical :: ifFound, eof
    character(len=*), parameter :: sub_name = 'get_FRP_dataset'

!call msg('Total variation for LU from get_FRP_dataset',pFMD%pDiurnalVarTot(:,2))
!call msg('Total variation for time from get_FRP_dataset',pFMD%pDiurnalVarTot(12,:))

    !
    ! Read the file
    !
    if(error)return
    open(unit=uIn,file=chFNm,status='old', iostat=iStat)
    if(fu_fails(iStat==0,'Cannot open FRP dataset file:' + chFNm, sub_name))return

    ! Get array sizes from the header
    nFires = int_missing
    nMaxObs = int_missing
    iLine = 0
    do iCount = 1,20  ! Check first 20 lines
      call next_line_from_input_file(uIn, lineTmp, eof)
      if(error .or. eof) exit
      iLine = iLine + 1
      read(unit=lineTmp, fmt =*, iostat=iStat1) strTmp, chSizeUnit, iTmp
      if(iStat1 /= 0) continue ! Wrong format, probably
      if (index(strTmp, "number_of_fires") == 1) nFires = iTmp
      if (index(strTmp, "max_number_of_same_fire_observations") == 1) nMaxObs = iTmp
      if (nFires > 0 .and. nMaxObs > 0) exit
      if (index(strTmp,'fire ') == 1)then
        call set_error('Fire line before the header ends:' + lineTmp + ', line=' + fu_str(iCount), &
                     & 'get_FRP_dataset')
        return
      endif
    enddo
    if(fu_fails(iStat==0,'Failed to parse FRP dataset header: '//trim(chFNm), sub_name))return

    if( nFires < 1 .or. nFires > 1000000) then
       call set_error('Strange number of fires:' + fu_str(nFires), sub_name)
       return
    endif
    if(nMaxObs < 1 .or. nMaxObs > 100) then
       call set_error( 'Strange number of obs:' + fu_str(nMaxObs), sub_name)
       return
    endif

    pSet%nFires = nFires
    pSet%nMaxObs = nMaxObs
    allocate(pSet%pLonGeo(nFires), pSet%fX(nFires), pSet%fY(nFires), &
           & pSet%dx(nFires,nMaxObs),  pSet%dy(nFires,nMaxObs), &
           & pSet%pFRP(nFires,nMaxObs+2), &  ! daily mean
           & pSet%pTA(nFires,nMaxObs),  pSet%pT4(nFires,nMaxObs),  pSet%pT4b(nFires,nMaxObs), &
           & pSet%pT11(nFires,nMaxObs),  pSet%pT11b(nFires,nMaxObs),  pSet%pMCE(nFires,nMaxObs), &
           & pSet%pArea(nFires,nMaxObs),  pSet%day(nFires),  pSet%pHour(nFires,nMaxObs), &
           & pSet%indLU(nFires), stat=iStat)

    if(fu_fails(iStat==0,'Failed fire source allocation','get_FRP_dataset'))return
    pSet%pLonGeo = real_missing
    pSet%fX = real_missing
    pSet%fY = real_missing
    pSet%dx = real_missing
    pSet%dy = real_missing
    pSet%pFRP = real_missing
    pSet%pTA = real_missing
    pSet%pT4 = real_missing
    pSet%pT4b = real_missing
    pSet%pT11 = real_missing
    pSet%pT11b = real_missing
    pSet%pMCE = real_missing
    pSet%pArea = real_missing
    pSet%pHour = real_missing
    pSet%day = time_missing
    pSet%indLU = int_missing
    first_day = time_missing
    last_day = time_missing
    chSizeUnit1=""
    chFRPUnit1=""

    pSet%indFRPDaily_tot = pSet%nMaxObs + 1
    pSet%indFRPDaily_perFire = pSet%nMaxObs + 2

    
!call msg('************* 1',pSet%pFRP(pSet%nFires/2,1:pSet%nMaxObs+2))    
    
    !
    ! Having metadata defined, can associate each fire with the specisic land use.
    ! Strictly speaking, different observations can shift the fire a bit. This may 
    ! mean different land use. In this case we shall take the one that is associated with 
    ! higher FRP. Will do it on the fly.
    !
    call grid_dimensions(pFMD%grid_metadata, nx, ny)

    !
    ! Format of the line:
    !
    ! ind  -----date-----    lon  lat    dx dy unit   FRP    TA  T4  T4b T11 T11b   MCE  area
    !  1  2010 5 1 12 0 0.0  20.00 50.00 1.5 1.5 km   1.5 MW  350 350 300 330 290   0.5  0.05
    !

!call msg('************* 2',pSet%pFRP(pSet%nFires/2,1:pSet%nMaxObs+2))    

    do while (.not. iStat < 0) ! not EOF
      iLine = iLine + 1
      call next_line_from_input_file(uIn, lineTmp, eof)
      if(error .or. eof) exit
      if (index(lineTmp,'fire ') /= 1)then
        if(index(lineTmp,'END_FRP_DATASET_V1') == 1)exit
        call msg(lineTmp, iLine)
        call msg_warning('Non-fire line in the fire file body','get_FRP_dataset')
        cycle
      endif
!      read(uIn,fmt=*,iostat=iStat) strTmp, strTmp1, iFire, yr, mon, day, hr, mn, sec, &
      read(unit=lineTmp,fmt=*,iostat=iStat) strTmp, strTmp1, iFire, yr, mon, day, hr, mn, sec, &
                       & fLonTmp, fLatTmp, &
                       & fDx, fDy, chSizeUnit, fFRP, chFRPUnit, &
                       & fpT4, fpT4b, fpT11, fPT11b, fpTA, fpMCE, fParea
      strTmp=trim(strTmp)
      if (iStat /= 0) then
!        if (strTmp(1:18) == "END_FRP_DATASET_V1") then ! Normal end of the source file
!          iStat = 0
!          exit
!        endif
!        if (iStat > 0) then 
!          if (any(strTmp(1:1) == (/'#','!'/))) cycle
          call msg(trim(chFNm))
          call msg("Trouble inthe above  file, line ", iLine)
!        endif
!        if (strTmp(1:5) /= "fire") then !Last check
!          call msg("Non-fire line in file"//trim(chFNm)//", line ", iLine)
!          iStat = 0
!          cycle
!        endif
!        cycle
        exit
      endif

      if( iFire < 1  .or. iFire > pSet%nFires) then

          call set_error('Strange fire number in file: '//trim(chFNm)//", line "//trim(fu_str(iLine)), &
          & sub_name)
        return
      endif
!      read(unit=strTmp,fmt=*,iostat=iStat) iFire
!call msg('iFire',iFire)
      do iObs = 1, pSet%nMaxObs
        if(abs(pSet%pFRP(iFire,iObs) / real_missing - 1.) < 0.01)exit   ! get the free place
      end do
!call msg('FRP 1:',pSet%pFRP(iFire,1:iObs))

      ! Do some unit magic
      if (chSizeUnit /= chSizeUnit1) then
         fSizeConv = fu_conversion_factor(chSizeUnit, 'm')
         chSizeUnit1 = chSizeUnit
      endif
      if (chFRPUnit /= chFRPUnit1) then
         fPowerConv = fu_conversion_factor(chFRPUnit, 'W')
          chFRPUnit1 = chFRPUnit
      endif

      pSet%dx(iFire,iObs) = fDx * fSizeConv
      pSet%dy(iFire,iObs) = fDy * fSizeConv
      pSet%pFRP(iFire,iObs) = fFRP *  fPowerConv
      pSet%pT4(iFire,iObs) = fpT4
      pSet%pT4b(iFire,iObs) = fpT4b
      pSet%pT11(iFire,iObs) = fpT11
      pSet%pT11b(iFire,iObs) = fPT11b
      pSet%pTA(iFire,iObs) = fpTA
      pSet%pMCE(iFire,iObs) = fpMCE
      pSet%pArea(iFire,iObs) = fParea
      if(pSet%pFRP(iFire,iObs) < 1e-3)then
        pSet%pFRP(iFire,iObs) = real_missing
        cycle
      endif
!call msg('FRP 2:',pSet%pFRP(iFire,1:iObs))
      !
      ! Time operations:
      ! - store the source activity period
      ! - store the actual time of the day of each fire emission
      !   Note that single fire is "valid" only one day. Next day it will be another fire
      !
      timeTmp = fu_set_time_UTC(yr,mon,day,0,0,0.0) ! points at the day start
      if(iObs > 1)then   ! one observation already there
        if(.not. timeTmp == pSet%day(iFire))then
          call set_error('Same-fire data belong to different days:' + fu_str(timeTmp) + &
                             & ',' + fu_str(pSet%day(iFire)),'get_FRP_dataset')
          return
        endif
      else
        pSet%day(iFire) = timeTmp  ! first fire observation decides its day
      endif  ! if one obs already there
      !
      ! Check what obs exist and look if this FRP is the biggest. Store its landuse if yes
      !
      ifFound = .false.
      do iTmp = 1, iObs - 1
        if(pSet%pFRP(iFire,iObs) < pSet%pFRP(iFire,iTmp))then
          ifFound = .true.
          exit
        endif
      end do
      if(.not. ifFound)then ! no one bigger than the new-comer
        call project_point_to_grid(fLonTmp, fLatTmp, pFMD%grid_metadata, fX, fY)
        if(error)return
        if(fX < 0.5 .or. fY < 0.5 .or. fX > nx+0.5 .or. fY > ny+0.5)then
!call msg('Fire out of grid. lon,lat,fXnew,fYnew:', (/fLonTmp, fLatTmp, fX, fY/))
!call report(pFMD%grid_metadata)
          cycle  ! out of grid
        endif  ! out of grid
        pSet%indLU(iFire) = pFMD%indLU(max(min(nint(fX),nx),1), max(min(nint(fY),ny),1))
        pSet%pLonGeo(iFire) = fLonTmp
        pSet%fX(iFire) = fLonTmp
        pSet%fY(iFire) = fLatTmp
!if(pSet%indLU(iFire) < 0 .or. pSet%indLU(iFire) > 10000)then
!call msg('strange land use:',pSet%indLU(iFire))
!endif
      endif  ! no FRP found bigger than new

      pSet%pHour(iFire,iObs) = hr + real(mn)/60. + sec/3600. ! fraction of the day
      if(defined(first_day))then      ! speed-up for the future
        if(timeTmp < first_day)then
          first_day = timeTmp
        elseif(timeTmp > last_day)then
          last_day = timeTmp
        endif
      else
        first_day = timeTmp
        last_day = timeTmp
      endif
    enddo  ! fire file loop
    close(uIn)
    if (iStat /= 0) then

      call set_error("Non-zero status after finished with FRP file: "//trim(chFnm), sub_name)
      return
    endif



    !call msg('************* 3',pSet%pFRP(pSet%nFires/2,1:pSet%nMaxObs+2))    
    !
    ! Clean the mess: may be, some fires are out of grid or corrupted or removed from the list.
    ! This is time to eliminate them.
    !
    iCount = 0
    do iFire = 1, pSet%nFires
      if(pSet%indLU(iFire) == int_missing)then
!        call msg_warning('Removing void fire:'+fu_str(iFire)+', file:'+chFNm,'get_FRP_dataset')
        iCount = iCount + 1
        call remove_fire(pSet, iFire)
      endif
    end do
    if(iCount > 0)call msg('Void fires:' + fu_str(iCount) + ',' + chFNm)
    !
    ! Having collected all fire observations, let's compute daily-mean FRP. 
    ! take unceratinty into account: very low FRPs should not give too much impact if
    ! there are higher values available at other times.
    !
    do iFire = 1, pSet%nFires
!call msg('FRP 3:',pSet%pFRP(iFire,1:pSet%nMaxObs+2))
      if(pSet%indLU(iFire) == int_missing)cycle ! e.g. out of grid
      pSet%pFRP(iFire,pSet%indFRPDaily_tot) = 0
      pSet%pFRP(iFire,pSet%indFRPDaily_perFire) = 0
      variance_sum = 0.
      do iObs = 1, pSet%nMaxObs
        if(pSet%pFRP(iFire,iObs) > 0.)then

          variance_inv = 1. / (sigma_FRP_abs **2 + (pSet%pFRP(iFire,iObs) * sigma_FRP_rel)**2)
          !
          ! Get actual variation
          !
          fTmp = pSet%pHour(iFire,iObs) + pSet%pLonGeo(iFire)/15.0
          if(fTmp < 0.)fTmp = fTmp + 24.
          if(fTmp >= 24.)fTmp = fTmp - 24.
          iTmp = int(fTmp)
          fWeightPast = real(iTmp + 1) - fTmp

          !
          ! Daily value if total-FRP diurnal variation is applied
          !
!call msg('iTmp, iFire, iLU',(/iTmp,iFire,pSet%indLU(iFire)/))
!call msg('Total variation for LU',pFMD%pDiurnalVarTot(:,pSet%indLU(iFire)))
!call msg('Total variation for time',pFMD%pDiurnalVarTot(iTmp+1,:))
          fTmp = pFMD%pDiurnalVarTot(iTmp+1,pSet%indLU(iFire)) * fWeightPast + &
               & pFMD%pDiurnalVarTot(min(24,iTmp+2),pSet%indLU(iFire)) * (1. - fWeightPast)
!call msg('fTmp, variance_inv, fWeightPast:',(/fTmp, variance_inv, fWeightPast/))
          pSet%pFRP(iFire,pSet%indFRPDaily_tot) = pSet%pFRP(iFire,pSet%indFRPDaily_tot) + &
                                                  & pSet%pFRP(iFire,iObs) / fTmp * variance_inv
          !
          ! Daily value if per-fire FRP diurnal variation is applied
          !
          fTmp = pFMD%pDiurnalVarPerFire(iTmp+1,pSet%indLU(iFire)) * fWeightPast + &
               & pFMD%pDiurnalVarPerFire(min(24,iTmp+2),pSet%indLU(iFire)) * (1. - fWeightPast)
          
          pSet%pFRP(iFire,pSet%indFRPDaily_perFire) = pSet%pFRP(iFire,pSet%indFRPDaily_perFire) + &
                                                    & pSet%pFRP(iFire,iObs) / fTmp * variance_inv
          !
          ! Variance is the same for both
          !
          variance_sum = variance_sum + variance_inv
        endif
      end do
!call msg('FRP 4:',pSet%pFRP(iFire,1:iObs))
      pSet%pFRP(iFire,pSet%indFRPDaily_tot) = pSet%pFRP(iFire,pSet%indFRPDaily_tot) / variance_sum
      pSet%pFRP(iFire,pSet%indFRPDaily_perFire) = pSet%pFRP(iFire,pSet%indFRPDaily_perFire) / &
                                                & variance_sum
!call msg('FRP 5:',pSet%pFRP(iFire,1:iObs))
    end do  ! iFire
!call msg('************* 4',pSet%pFRP(pSet%nFires/2,1:pSet%nMaxObs+2))    

    if(pSet%nFires <= 0)then
      call msg_warning('No fires found in:' + chFNm + ', had to disable the set')
      first_day = really_far_in_future
      last_day = really_far_in_past
    endif
!call msg('FRP dataset FRP:',pSet%pFRP(1:pSet%nFires,pSet%indFRPDaily_tot))

  end subroutine get_FRP_dataset


  !*********************************************************************************
    
  function fu_get_fire_metadata(chFMD_FNm, expected_species, nlSetup) result (pFMD)
    !
    ! Searches for the metadata from the same ini file already consumed into fire sources
    !
    implicit none
      
    ! Imported parameter 
    character(len=*), intent(in) :: chFMD_FNm
    type(silam_species), dimension(:), intent(in), allocatable :: expected_species
    type(Tsilam_namelist), intent(in) :: nlSetup

    ! result of the function
    type(Tsilam_fire_metadata), pointer :: pFMD

    integer, parameter :: max_LU = 700 ! max_LU = 100, now it distinguish +6 continents: AU, AS, AF, EU, NA, SA
    integer, parameter :: max_subst = max_species
    ! Local variables
    integer :: uIn, iStat, nItems, iItem, iLU, iFMD, iTmp, nSubst, iSubst, nx, ny, ix, iy, &
             & nSpeciesAerosol, iSp, nSpFlam, nSpSmld, nGrids
    type(Tsilam_namelist), pointer :: nlMD
    type(Tsilam_nl_item_ptr), dimension(:), pointer :: pItems
    character(len=fnlen) :: strTmp! Namelist content
    character(len=substNmLen) :: chSubst, chPhase
    character(len=20) :: chUnit
    character(len=30), dimension(max_LU) :: chLU
    integer, dimension(max_LU) :: arLUFile
    integer, dimension(max_subst) :: iPhase
    real :: fact_flame, fact_smld
    real, dimension(2) :: fMassMeanDiam
    real, dimension(:), pointer :: arMapTmp
    real, dimension(max_LU,max_subst) :: pEmisFactorsTmp_flaming, pEmisFactorsTmp_smouldering !iLUxiSubst
    real, dimension(2,max_species)  :: fluxPerModeNbr, fluxPerModeVol
    type(silam_material_ptr), dimension(max_subst) :: pMaterials
    type(silam_species) :: speciesTmp
    type(silam_species), dimension(:), pointer :: speciesAerosol 
    type(silja_grid), dimension(:), pointer :: grids
    type(silja_time), dimension(:), pointer :: timeLst

    speciesAerosol => null()

    !
    ! Already exists?
    !
    nullify(pFMD)
    do iFMD = 1, nFMD_glob
       if(trim(arFMD_glob(iFMD)%ptr%chIniFNm) == trim(chFMD_FNm))then
         pFMD => arFMD_glob(iFMD)%ptr
        return
      endif
    end do
    !
    ! Have not found anything. Read new!
    !
    allocate(arFMD_glob(nFMD_glob+1)%ptr, stat=iStat)
    if(fu_fails(iStat==0,'Failed Metadata allocation','fu_get_fire_metadata'))return
    nFMD_glob = nFMD_glob + 1
    pFMD => arFMD_glob(nFMD_glob)%ptr

    nullify(pItems)
    uIn = fu_next_free_unit()
    if(error)return
    pFMD%chIniFNm = chFMD_FNm
    
    pEmisFactorsTmp_flaming(:, :) = 0.0
    pEmisFactorsTmp_smouldering(:, :) = 0.0
    
    !
    ! Open and read the content
    !
    open(unit=uIn,file=chFMD_FNm, iostat=iStat)
    if(fu_fails(iStat==0,'Cannot open fire metadata file:' + chFMD_FNm, 'fu_get_fire_metadata'))return
    nlMD => fu_read_namelist(uIn, .false.,'END_FIRE_METADATA_V1')
    close(uIn)
    !
    ! List of lumped land-uses we distinguish and substances, for which we have emission factors
    ! They are mixed together, so have to first sort out the list
    !
    fStemMassFraction = fu_content_real(nlMD,'fraction_of_stem_mass')
    if(error .or. fStemMassFraction < 0. .or. fStemMassFraction > 1.)then
      if(error)call unset_error('fu_get_fire_metadata')
      call set_error('Strange fraction_of_stem_mass in fire metadata file:' + &
                    & fu_str(fStemMassFraction),'fu_get_fire_metadata')
      fStemMassFraction = -1.0
      return
    endif
    !
    ! Emission factors are 3-D: land use, substance and gas/aerosol.
    ! For aerosols, size spectrum is defined further
    !
    call get_items(nlMD,'emission_factor', pItems, nItems)
    if(fu_fails(nItems > 0, 'No emission_factor in metadata file:' + chFMD_FNm, &
                          & 'fu_get_fire_metadata'))return
    pFMD%nLU_types = 0
    nSubst = 0
    do iItem = 1, nItems
      strTmp = fu_content(pItems(iItem))
      read(unit=strTmp,fmt=*)chLU(pFMD%nLU_types+1), chSubst, chPhase, fact_flame, fact_smld     ! read all pieces but unit
      chUnit = trim(strTmp(index(trim(strTmp),' ',.true.)+1:))  ! unit is from last space till the end
      !
      ! Put the stuff to 3 arrays: LU types, materials, 2D emission factor array
      ! Mind the repetitions.
      !
      iLU = pFMD%nLU_types+1
      do iTmp = 1, pFMD%nLU_types
        if(trim(chLU(iTmp)) == trim(chLU(pFMD%nLU_types+1)))then
          iLU = iTmp
          exit
        endif
      end do
      if(iLU == pFMD%nLU_types+1) pFMD%nLU_types = pFMD%nLU_types + 1
      !
      ! Define the species: detect the material and phase
      !
      nSubst = nSubst + 1
      pMaterials(nSubst)%ptr => fu_get_material_ptr(chSubst)   ! store substance and emission phase
      if(fu_str_u_case(chPhase) == 'GAS')then
        iPhase(nSubst) = gas_phase_flag
      elseif(fu_str_u_case(chPhase) == 'FIRE_PM_SPECTRUM')then
        iPhase(nSubst) = fire_mode_flag
      else
        call set_error('Unknown phase:' + chPhase + '. Allowed: GAS or FIRE_PM_SPECTRUM','get_fire_metadata')
        return
      endif
      !
      ! Check for duplicates
      !
      iSubst = nSubst
      do iTmp = 1, nSubst-1                                      ! check for duplicates
        if(associated(pMaterials(nSubst)%ptr, pMaterials(iTmp)%ptr) .and. &
         & iPhase(nSubst) == iPhase(iTmp))then
          iSubst = iTmp
          nSubst = nSubst - 1     ! a duplicate found, remove the last species
          exit
        endif
      end do
      !
      ! Store the emission factors in basic SI unit
      !
      pEmisFactorsTmp_flaming(iLU,iSubst) = fact_flame * &
                  & fu_conversion_factor(chUnit, fu_basic_mass_unit(pMaterials(iSubst)%ptr) + '/J', &
                                       & pMaterials(iSubst)%ptr)
      pEmisFactorsTmp_smouldering(iLU,iSubst) = fact_smld * &
                  & fu_conversion_factor(chUnit, fu_basic_mass_unit(pMaterials(iSubst)%ptr) + '/J', &
                                       & pMaterials(iSubst)%ptr)
    end do  ! emission_factor list
    if(error)return
    !
    ! Spectrum of the particles
    !
    select case(fu_str_u_case(fu_content(nlMD,'fire_aerosol_spectrum')))

      case ('LOGNORMAL_THREE_MODES')
        pFMD%iSpectrumType = lognorm_3_modes_fires_num_flag  !lognorm_3_modes_fires_flag

      case default
        call set_error('Unknown spectrum type:' + fu_content(nlMD,'fire_aerosol_spectrum'), &
                     & 'fu_get_fire_metadata')
        return
      end select
    !
    ! Can now fill-in the species list and store everything.
    ! We can have gaseous and aerosol species. The list is dictated by the emission factors above. 
    ! From the other side, fires emit very specific aerosol size spectrum defined above via its name.
    ! Emission factors therefore cannot refer to individual size bins, ONLY to whole-spectrum.
    ! PM10 is forbidden, only total PM, etc.
    ! The bins onto which the spectrum is projected are defined explicitly - the same for all 
    ! fire aerosol species. To set them, we will use the standard aerosol procedure.
    ! Note that there are two potentially concurring definitions: one coming from aerosol dynamics,
    ! the other - written in the fire ini file. The first one prevails, if exists.
    ! Emission factors for cocktails must be split down to species. The result to produce here is 
    ! the list of species and related emission factors.
    ! The species for flaming and smouldering are the same, except for mass-mean diameter - but that one
    ! is ignored in species comparison (see chemical setup compare_modes_eq).
    !
    pFMD%nSpecies = 0
    nSpFlam = 0
    nSpSmld = 0
    do iSubst = 1, nSubst
      if(iPhase(iSubst) == gas_phase_flag)then
        call set_species(speciesTmp, pMaterials(iSubst)%ptr, in_gas_phase)
        call addSpecies(pFMD%species_flaming, nSpFlam, (/speciesTmp/), 1)
        call addSpecies(pFMD%species_smouldering, nSpSmld, (/speciesTmp/), 1)
      elseif(iPhase(iSubst) == fire_mode_flag)then
        call get_source_aer_species(speciesAerosol, nSpeciesAerosol, &       ! species for given subst, all modes
                                  & expected_species, fu_name(pMaterials(iSubst)%ptr), &
                                  & nlSetup)
        call addSpecies(pFMD%species_flaming, nSpFlam, speciesAerosol, nSpeciesAerosol) ! mean diameter is not set
        call addSpecies(pFMD%species_smouldering, nSpSmld, speciesAerosol, nSpeciesAerosol) ! mean diameter is not set
      else
        call set_error('Unknown phase:'+fu_str(iPhase(iSubst)),'fu_get_fire_metadata')
      endif
    end do
    if(error)return
    !
    ! Having the number of spceies decided, allocate the space and set the emission factors
    ! Later on, the number flaming and smouldering species might become different but sofar it is the same
    !
    if(nSpFlam /= nSpSmld)then
      call msg('Number of flaming and smouldering species:',nSpFlam, nSpSmld)
      call set_error('Different number of flaming and smouldering species','')
      return
    endif
    pFMD%nSpecies = nSpFlam
    allocate(pFMD%ifNumberFlux(pFMD%nSpecies), &
           & pFMD%pEmsFactors(2, pFMD%nLU_types, pFMD%nSpecies), stat=iTmp)
    if(fu_fails(iTmp == 0, 'Failed to allocate species_related arrays','fu_get_fire_metadata'))return
    pFMD%pEmsFactors = 0.0
    !
    ! Finally, set the emission factors and the mean diameters for aerosols
    !
    do iSubst = 1, nSubst
      if(iPhase(iSubst) == gas_phase_flag)then
        !
        ! Gas added as is, emission factor copied. Duplicates have been checked above,
        ! can just use the last-added species index
        !
        call set_species(speciesTmp, pMaterials(iSubst)%ptr, in_gas_phase)
        iSp = fu_index(speciesTmp, pFMD%species_flaming, nSpFlam)
        pFMD%pEmsFactors(flaming, 1:pFMD%nLU_types, iSp) = &
                                              & pEmisFactorsTmp_flaming(1:pFMD%nLU_types, iSubst)
        iSp = fu_index(speciesTmp, pFMD%species_smouldering, nSpSmld)
        pFMD%pEmsFactors(smouldering, 1:pFMD%nLU_types, iSp) = &
                                              & pEmisFactorsTmp_smouldering(1:pFMD%nLU_types, iSubst)
        pFMD%ifNumberFlux(iSp) = .false.
          
      elseif(iPhase(iSubst) == fire_mode_flag)then
        !
        ! Aerosols are blown-up using the prescribed spectrum and given name
        ! Overlay aerosol modes listed in the fire source term with the given substance
        !
        call get_source_aer_species(speciesAerosol, nSpeciesAerosol, &       ! species for given subst, all modes
                                  & expected_species, fu_name(pMaterials(iSubst)%ptr), &
                                  & nlSetup)
        if(error .or. nSpeciesAerosol < 1)then
          call set_error('Failed to get source aerosol modes for substance:' + &
                                               & fu_name(pMaterials(iSubst)%ptr),'get_fire_metadata')
          return
        endif
        do iTmp = 1, nSpeciesAerosol
          call fires_flux4mode(speciesAerosol(iTmp)%mode, &  ! mode to project to (flaming)
                             & pFMD%iSpectrumType, &        ! smoke size distribution 
                             & fluxPerModeNbr(:,iTmp), &     ! number flux for this mode, flames/smld
                             & fluxPerModeVol(:,iTmp), &     ! volume fluxes for this mode, flames/smld
                             & fMassMeanDiam(:))            ! mass-weighted mean diameter, flames/smld
          if(error)return
call msg('mass mean diam flames/smld:',fMassMeanDiam(flaming),fMassMeanDiam(smouldering))
          !
          ! Set the species and emission factors one by one, note size spectrum fractions
          !
          iSp = fu_index(speciesAerosol(iTmp), pFMD%species_flaming, nSpFlam)

          ! This reset caused confusion in naming and identifying of modes
          ! Commenting out until nominal diameter introduced
      !         call msg ("FIXME!!! NOT Resetting flaming mode mass_mean_d from, to ", &
      !            &   fu_mean_d(pFMD%species_flaming(iSp)), fMassMeanDiam(flaming))

          call set_massmean_d(pFMD%species_flaming(iSp)%mode, fMassMeanDiam(flaming), .false.) ! if dynamic diameter
          if(fu_name(pMaterials(iSubst)%ptr) == 'nbr_aer')then
            pFMD%pEmsFactors(flaming, 1:pFMD%nLU_types, iSp) = &
                                              & pEmisFactorsTmp_flaming(1:pFMD%nLU_types, iSubst) * &
                                              & fluxPerModeNbr(flaming,iTmp)           ! NUMBER
            pFMD%ifNumberFlux(iSp) = .true.
          else
            pFMD%pEmsFactors(flaming, 1:pFMD%nLU_types, iSp) = &
                                              & pEmisFactorsTmp_flaming(1:pFMD%nLU_types, iSubst) * &
                                              & fluxPerModeVol(flaming, iTmp)        ! VOLUME
            pFMD%ifNumberFlux(iSp) = .false.
          endif

          iSp = fu_index(speciesAerosol(iTmp), pFMD%species_smouldering, nSpSmld)
          !          call msg ("FIXME!!! NOT Resetting smouldering mode mass_mean_d from, to ", &
          !        &   fu_mean_d(pFMD%species_smouldering(iSp)), fMassMeanDiam(smouldering))
          call set_massmean_d(pFMD%species_smouldering(iSp)%mode, fMassMeanDiam(smouldering), .false.) ! if dynamic diameter
          if(fu_name(pMaterials(iSubst)%ptr) == 'nbr_aer')then
            pFMD%pEmsFactors(smouldering, 1:pFMD%nLU_types, iSp) = &
                                                 & pEmisFactorsTmp_smouldering(1:pFMD%nLU_types, iSubst) * &
                                                 & fluxPerModeNbr(smouldering,iTmp)
          else
            pFMD%pEmsFactors(smouldering, 1:pFMD%nLU_types, iSp) = &
                                                 & pEmisFactorsTmp_smouldering(1:pFMD%nLU_types, iSubst) * &
                                                 & fluxPerModeVol(smouldering, iTmp)
          endif
        end do  ! aerosol species
      else
        call set_error('Unknown phase','fu_get_fire_metadata')
      end if  ! gas / particle
      if(error)return

    end do  ! iSubst  

    !
    ! Deal with land-use types
    !
    allocate(pFMD%chLU_names(pFMD%nLU_types), &
           & pFMD%pDiurnalVarTot(24, pFMD%nLU_types), pFMD%pDiurnalVarPerFire(24,pFMD%nLU_types), &
           & pFMD%flaming_fraction_roughly(pFMD%nLU_types, pFMD%nSpecies), &
           & stat=iStat)
    if(fu_fails(iStat==0,'Failed names & metadata allocation', 'fu_get_fire_metadata'))return
    do iTmp = 1, pFMD%nLU_types
      pFMD%chLU_names(iTmp) = chLU(iTmp)
    end do
    !
    ! Next is the land-use-dependent duirnal variation. Two types considered: 
    ! - total-FRP variation
    !
    call get_items(nlMD,'hour_in_day_index_total',pItems,nItems)
    if(error .or. fu_fails(nItems == pFMD%nLU_types, &
                         & 'Number of total diurnal variations (' + fu_str(nItems) + ') /= number of land use types (' + &
                         & fu_str(pFMD%nLU_types) + ')', 'fu_get_fire_metadata'))return
    do iItem = 1, nItems
      strTmp = fu_content(pItems(iItem))
      iLU = 0
      do iStat = 1, pFMD%nLU_types
        if(index(strTmp,trim(pFMD%chLU_names(iStat))) == 1)then
          iLU = iStat
          exit
        endif
      end do
      if(fu_fails(iLU>0,'Unknown landuse in line:' + strTmp,'fu_get_fire_metadata'))return
      read(unit=strTmp,fmt=*,iostat=iStat)chLU(1),(pFMD%pDiurnalVarTot(ix,iLU),ix=1,24)
      if(fu_fails(iStat==0,'Failed to read the hour_in_day_index_total:'+strTmp,'fu_get_fire_metadata'))return
    end do ! hour_in_day_index_total list
    !
    ! ... and FRP-per-pixel variation
    !
    call get_items(nlMD,'hour_in_day_index_per_fire',pItems,nItems)
    if(error .or. fu_fails(nItems == pFMD%nLU_types, &
                         & 'Number of per-fire diurnal variations (' + fu_str(nItems) + ') /= number of land use types (' + &
                         & fu_str(pFMD%nLU_types) + ')', 'fu_get_fire_metadata'))return
    do iItem = 1, nItems
      strTmp = fu_content(pItems(iItem))
      iLU = 0
      do iStat = 1, pFMD%nLU_types
        if(index(strTmp,trim(pFMD%chLU_names(iStat))) == 1)then
          iLU = iStat
          exit
        endif
      end do
      if(fu_fails(iLU>0,'Unknown landuse in line:' + strTmp,'fu_get_fire_metadata'))return
      read(unit=strTmp,fmt=*,iostat=iStat)chLU(1),(pFMD%pDiurnalVarPerFire(ix,iLU),ix=1,24)
      if(fu_fails(iStat==0,'Failed to read the hour_in_day_index_per_fire:'+strTmp,'fu_get_fire_metadata'))return
    end do
    !
    ! Knowing the metadata, one can roughly calculate the contribution of flaming and smouldering into
    ! the daily emission. This is useful when need quick-and-dirty fire daily total.
    ! Note that it is NOT thre real contribution but a rough estimate of it
    !
    do iLU = 1, pFMD%nLU_types
      do iSp = 1, pFMD%nSpecies
        pFMD%flaming_fraction_roughly(iLU, iSp) = 0.
        do iItem = 1, 24
          pFMD%flaming_fraction_roughly(iLU, iSp) = pFMD%flaming_fraction_roughly(iLU, iSp) + &
                                & fu_emission_weighted_flam_smld(1.0, 1.0, &            ! FRP, seconds
                                                        & pFMD%pEmsFactors(:,iLU,iSp), &
                                                        & pFMD%pDiurnalVarTot(iItem,iLU)) / 24.0
        end do  ! 1..24
        pFMD%flaming_fraction_roughly(iLU, iSp) = &
               & abs(pFMD%flaming_fraction_roughly(iLU, iSp) - pFMD%pEmsFactors(smouldering,iLU,iSp)) / &
               & (abs(pFMD%pEmsFactors(flaming,iLU,iSp) - pFMD%pEmsFactors(smouldering,iLU,iSp))+1e-30)
      end do  ! iSp
    end do  ! iLU
    !
    ! Lumping rules for land use types.
    !
    call get_items(nlMD,'land_use', pItems, nItems)
    if(error .or. fu_fails(nItems > 0, 'No land_use lines in the file:' + chFMD_FNm, &
                                     & 'fu_get_fire_metadata'))return
    do iItem = 1, nItems
      strTmp = fu_content(pItems(iItem))
      read(unit=strTmp,fmt=*) iLU, chLU(iLU)
      do iStat = 1, pFMD%nLU_types
        if(trim(pFMD%chLU_names(iStat)) == trim(chLU(iLU)))then
          arLUFile(iLU) = iStat   ! linking the in-file and lumped LUs
          exit
        endif
      end do
    end do
    !
    ! The last step: read the land use map. Note: it is huge since the resolution 
    ! must be 3 km at least
    !
    call msg('Reading the land use data. Can take time!')
    !
    ! Since the data are large, will do it manually. 
    ! We consider the dataset as coming from GrADS file, which is a large redundancy
    ! since integer*1 is enough for up to 256 land use types
    !
    strTmp = fu_process_filepath(fu_content(nlMD,'land_use_file'), superfile=chFMD_FNm)
    if(index(fu_str_u_case(strTmp),'GRADS') == 1)then
      !
      !  GrADS land use file
      !
      uIn = fu_open_gradsfile_i(fu_extend_grads_hat( strTmp(index(strTmp,' ')+1:), chFMD_FNm))
      if(error .or. fu_fails(uIn>0,'Failed opening of:' + strTmp,'fu_get_fire_metadata'))return

      pFMD%grid_metadata = fu_silamGrid_of_grads(uIn)
      call grid_dimensions(pFMD%grid_metadata, nx, ny)
      if(error)return

      allocate(pFMD%indLU(nx,ny), stat=iStat)
      if(fu_fails(iStat==0,'Failed allocation of landuse map','fu_get_fire_metadata'))return

      arMapTmp => fu_work_array(nx*ny)
      if(error)return

      call read_field_from_grads_id(uIn, &
                                   & fu_set_field_id_simple(met_src_missing, &
                                                          & land_use_type_flag, &
                                                          & fu_time_of_grads(uIn,1), &
                                                          & level_missing), &
                                   & arMapTmp)
      if(error)return
      call close_gradsfile_i(uIn)

    elseif(index(fu_str_u_case(strTmp),'NETCDF') == 1)then
      !
      !  NetCDF land use file
      !
      !strTmp = fu_process_filepath(fu_content(nlMD,'land_use_file'), superfile=chFMD_FNm)
      print *, strTmp
      print *, chFMD_FNm
      print *, fu_content(nlMD,'land_use_file')
      print *, fu_str_u_case(strTmp)
      print *, strTmp(index(strTmp,' ')+1:)
      uIn = open_netcdf_file_i(strTmp(index(strTmp,' ')+1:), fu_input_file_format('NETCDF:LandUseType'))
      if(error .or. fu_fails(uIn>0,'Failed opening of:' + strTmp,'fu_get_fire_metadata'))return

      call timelst_from_netcdf_file(uIn, timeLst, nx) ! nx -- returned length of time list (not to be used here)
      if(error)return

      call get_netcdf_grids(uIn, land_use_type_flag, species_missing, grids, nGrids)
      if(fu_fails(nGrids==1,'Strange number of grids:' + fu_str(nGrids),'fu_get_fire_metadata'))return
      pFMD%grid_metadata = grids(1)
      call grid_dimensions(pFMD%grid_metadata, nx, ny)
      if(error)return


      allocate(pFMD%indLU(nx,ny), stat=iStat)
      if(fu_fails(iStat==0,'Failed allocation of landuse map','fu_get_fire_metadata'))return

      arMapTmp => fu_work_array(nx*ny)
      if(error)return

      call read_field_from_netcdf_file(uIn, &
                                     & fu_set_field_id_simple(met_src_missing, &
                                                            & land_use_type_flag, &
                                                            & timeLst(1), &
                                                            & level_missing), &
                                     & arMapTmp, real_missing)
      if(error)return
      call close_netcdf_file(uIn)
   else
      call set_error('Only GRADS and NetCDF files supported, not:' + strTmp,'fu_get_fire_metadata')
      return
    endif

    !
    ! Find the lumping relation and store lumped land use
    !
    iItem = 1
    do iy = 1, ny
      do ix = 1,nx
        if(mod(iItem,5000000)==0)call msg('Processed:',iItem,nItems)
        if(arMapTmp(iItem) <= 0. .or. arMapTmp(iItem) > size(arLUFile))then
          call set_error('Strange LU type at' + fu_str(iItem) + ':' + fu_str(arMapTmp(iItem)), &
                       & 'fu_get_fire_metadata')
          return
        endif
        iStat = nint(arMapTmp(iItem))
        if(arLUFile(iStat) <= 0 .or. arLUFile(iStat) > pFMD%nLU_types)then
          call set_error('LU type not lumped at:' + fu_str(iItem) + ':' + fu_str(arMapTmp(iItem)), &
                      & 'fu_get_fire_metadata')
          return
        endif
        pFMD%indLU(ix,iy) = arLUFile(iStat)
        iItem = iItem + 1
      enddo
    enddo
call msg('finished processing')
    call free_work_array(arMapTmp)
    
    pFMD%defined = silja_true

  end function fu_get_fire_metadata


  !**********************************************************************

  subroutine clean_fire_metadata(pFMD)
    !
    ! Destroys the huge metadata arrays when these are mo longer needed.
    !
    implicit none
      
    ! Imported parameter 
    type(Tsilam_fire_metadata), pointer :: pFMD

    if(fu_true(pFMD%defined))deallocate(pFMD%indLU)  ! Only the main map
    pFMD%defined = silja_false
    pFMD%chIniFNm = ''

  end subroutine clean_fire_metadata


  !*****************************************************************

  subroutine add_input_needs_fire_src(fire_src, q_met_dynamic, q_met_static, &
                                                  & q_disp_dynamic, q_disp_static)
    !
    ! Returns input needs for the emission/transformation routines. 
    !
    implicit none

    ! Imported parameters
    type(Tsilam_fire_source), intent(in) :: fire_src
    integer, dimension(:), intent(inout) :: q_met_dynamic, q_met_static, &
                                          & q_disp_dynamic, q_disp_static
    ! Local variables
    integer :: iTmp
    !
    ! Add needed dynamic quantities. Injection height needs ABL height and Brunt-Vaisala frequency
    ! Detailed land use will not be requested here, it is a feature of the fire module.
    !
    iTmp = fu_merge_integer_to_array(abl_height_m_flag, q_met_dynamic)
    iTmp = fu_merge_integer_to_array(brunt_vaisala_freq_flag, q_met_dynamic)
    if (fu_leveltype(dispersion_vertical) .ne. layer_btw_2_height) then
      iTmp = fu_merge_integer_to_array(pressure_flag, q_met_dynamic)
      iTmp = fu_merge_integer_to_array(ground_pressure_flag, q_met_dynamic)
    endif

!    iTmp = fu_merge_integer_to_array(emis_factor_fire_flame_flag, q_met_static)
!    iTmp = fu_merge_integer_to_array(emis_factor_fire_smold_flag, q_met_static) 

  end subroutine add_input_needs_fire_src


  !**************************************************************************

  subroutine reserve_fire_source(fire_src, &     ! Src to initialise
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
    type(Tsilam_fire_source), intent(inout) :: fire_src
    integer, intent(in) :: iSrcNbr, iSrcIdNbr

    !
    ! Nullify the basic variables
    !
    fire_src%src_nm = ''
    fire_src%sector_nm = ''

    !
    ! Main source parameters - enough to identify it in the global information list
    !
    fire_src%src_nbr = iSrcNbr
    fire_src%id_nbr = iSrcIdNbr
    !
    ! Finally, mark the source as incomplete
    !
    fire_src%defined = silja_false

  end subroutine reserve_fire_source


  !*********************************************************************

  subroutine init_emission_fire(fire_src)
    !
    implicit none

    ! Imported parameters
    type(Tsilam_fire_source), intent(inout) :: fire_src
    
    !
    ! Free some space from debris
    !
    call clean_fire_metadata(fire_src%pFMD)

  end subroutine init_emission_fire


  !****************************************************************************

  subroutine add_inventory_fire_src(fire_src, species_list, nSpecies)
    !
    ! Get the species emitted by this source and make a list of
    ! them. The source must be initialized, of course.
    !
    implicit none
    type(Tsilam_fire_source), intent(in) :: fire_src
    type(silam_species), dimension(:), pointer :: species_list
    integer, intent(inout) :: nSpecies
    !
    ! Simply add species.
    ! Flaming and smouldering species are identical, no need to duplicate
    !
    call addSpecies(species_list, nSpecies, fire_src%pFMD%species_flaming, fire_src%pFMD%nSpecies)

  end subroutine add_inventory_fire_src


  !*******************************************************************************
  
  subroutine create_src_cont_grd_fire_src(fire_src, grid_template, ifVerbose, ifExtended)
    !
    ! Creates the grid that covers the area with active BVOC emission
    !
    implicit none

    ! Imported parameters
    type(Tsilam_fire_source), intent(in) :: fire_src
    type(silja_grid), intent(inout) :: grid_template
    logical, intent(in) :: ifVerbose
    logical, intent(out) :: ifExtended
    
    ! Local variables
    integer :: iSet, iFire, nx, ny
    real :: x,y, xMin, xMax, yMin, yMax

    
    if(.not. defined(grid_template))then
      !
      ! gridtemplate is undefined. Make it out of this source: just make lon-lat grid covering the fires
      ! with resolution of 0.1 degree
      !
      grid_template = fu_set_grid('', lonlat, pole_geographical, &
                                & fire_src%pFRPset(1)%fX(1)-0.1, fire_src%pFRPset(1)%fY(1)-0.1, &
                                & 3, 3, 0.1, 0.1)   ! nx, ny, dx, dy
    endif
    !
    ! Extend the grid to cover all fires
    !
    call grid_dimensions(grid_template, nx,ny)
    xMin = nx
    xMax = 1
    yMin = ny
    yMax = 1

    do iSet = 1, fire_src%nFRPdatasets
      do iFire = 1, fire_src%pFRPset(iSet)%nFires
        if(fire_src%pFRPset(iSet)%indLU(iFire) == int_missing)cycle
        if(fire_src%ifGeoCoord)then
          call project_point_to_grid(fire_src%pFRPset(iSet)%fX(iFire), &
                                   & fire_src%pFRPset(iSet)%fY(iFire), grid_template, x, y)
        else
          call project_point_to_grid(fire_src%grid, &
                                   & fire_src%pFRPset(iSet)%fX(iFire), &
                                   & fire_src%pFRPset(iSet)%fY(iFire), &
                                   & grid_template, x, y)
        endif
        if(error)return
        if(xMin > x)xMin = x
        if(xMax < x)xMax = x
        if(yMin > y)yMin = y
        if(yMax < y)yMax = y
      end do
    end do
    if(xMin<1 .or. xMax>nx .or. yMin<1 .or. yMax>ny) then
      ifExtended = .true.
      if(ifVerbose)then
        call msg('Grid nx and source xmin/max:' + fu_str(nx),xMin, xMax)
        call msg('Grid ny and source ymin/max:' + fu_str(ny),yMin, yMax)
        call msg('Extending the grid for the fire source:' + fire_src%src_nm)
        call report(fire_src)
      endif
      !
      ! Note that we must first extend the grid towards the high values and only then push 
      ! the starting indices: they are counted from 1, i.e. not changed when nx,ny expand
      !      
      call extend_grid_to_coordinates(grid_template, xMax, yMax)
      call extend_grid_to_coordinates(grid_template, xMin, yMin)
    else
      ifExtended = .false.
    end if
    
  end subroutine create_src_cont_grd_fire_src

  
  !*****************************************************************
  
  subroutine force_fire_source_into_grid(fire_src, grid_template, ifVerbose, ifCut)
    !
    ! Removes the fires not falling inside the grid
    !
    implicit none

    ! Imported parameters
    type(Tsilam_fire_source), intent(inout) :: fire_src
    type(silja_grid), intent(in) :: grid_template
    logical, intent(in) :: ifVerbose
    logical, intent(out) :: ifCut

    ! Local parameters
    integer :: iSet, iFire, nx, ny
    real :: x, y

    call grid_dimensions(grid_template, nx, ny)

    ifCut = .false.
    do iSet = 1, fire_src%nFRPdatasets
      iFire = 1
      do while(iFire <= fire_src%pFRPset(iSet)%nFires)
        if(fire_src%ifGeoCoord)then
          call project_point_to_grid(fire_src%pFRPset(iSet)%fX(iFire), &
                                   & fire_src%pFRPset(iSet)%fY(iFire), grid_template, x, y)
        else
          call project_point_to_grid(fire_src%grid, fire_src%pFRPset(iSet)%fX(iFire), &
                                                  & fire_src%pFRPset(iSet)%fY(iFire), &
                                   & grid_template, x, y)
        endif
        if(x<0.5 .or. x>nx+0.5 .or. y<0.5 .or. y>ny+0.5) then  ! Out of the grid
          call remove_fire(fire_src%pFRPset(iSet), iFire)
          ifCut = .true.
          if(ifVerbose) call msg('Cutting out-of-grid fire,' + fire_src%src_nm + ':',iFire)
          cycle
        endif
        iFire = iFire + 1
      end do
    end do

  end subroutine force_fire_source_into_grid


  !*****************************************************************

  subroutine project_fire_src_second_grd(fire_src, grid, vert_disp, vert_proj, iAccuracy)
    !
    ! There is nothing to project in terms of horizontal grid but the vertical has to be made.
    !
    implicit none

    ! Imported parameters
    type(Tsilam_fire_source), intent(inout) :: fire_src
    type(silja_grid), intent(in) :: grid
    type(silam_vertical), intent(in) :: vert_disp, vert_proj
    integer, intent(in) :: iAccuracy

    ! Local variables
    integer :: nx, ny, iFire, iSet
    real :: x, y

    if(len_trim(fire_src%sector_nm) > 0)then
      call msg('Re-projecting fire source:' + fire_src%src_nm +'_' + fire_src%sector_nm)
    else
      call msg('Re-projecting fire source:' + fire_src%src_nm)
    endif
    !
    ! Reproject the fire locations to the given grid. Use pLon and pLat - they are 
    ! not needed any more
    !
    call grid_dimensions(grid, nx, ny)
    if(error)return

    do iSet = 1, fire_src%nFRPdatasets
      iFire = 1
      do while(iFire <= fire_src%pFRPset(iSet)%nFires)
        if(fire_src%ifGeoCoord)then
          call project_point_to_grid(fire_src%pFRPset(iSet)%fX(iFire), fire_src%pFRPset(iSet)%fY(iFire), &
                                   & grid, x, y)
        else
          call project_point_to_grid(fire_src%grid, &
                                   & fire_src%pFRPset(iSet)%fX(iFire), fire_src%pFRPset(iSet)%fY(iFire), &
                                   & grid, x, y)
        endif
        if (nint(x) < 1 .or. nint(x) > nx .or. nint(y) < 1 .or. nint(y) > ny) then 
          ! Out of the grid - be consistent with rest of the fire source
          call remove_fire(fire_src%pFRPset(iSet), iFire)
          cycle
        else
          fire_src%pFRPset(iSet)%fX(iFire) = x
          fire_src%pFRPset(iSet)%fY(iFire) = y
        endif
        iFire = iFire + 1
      end do ! nFires
      if(fire_src%pFRPset(iSet)%nFires == 0)then
        fire_src%first_day(iSet) = really_far_in_future
        fire_src%last_day(iSet) = really_far_in_past
        cycle
      endif
    end do ! FRP datasets
    fire_src%grid = grid  ! register the new grid
    fire_src%ifGeoCoord = .false.

  end subroutine project_fire_src_second_grd


  !**************************************************************************
  
  subroutine remove_fire(pSet, iFireToRemove)
    !
    ! Removes the fire pointed by the index by copying the last fire over it
    !
    implicit none
    
    ! Imported parameters
    type(TFRP_dataset), intent(inout) :: pSet
    integer, intent(in) :: iFireToRemove

    pSet%pLonGeo(iFireToRemove) = pSet%pLonGeo(pSet%nFires)
    pSet%fX(iFireToRemove) = pSet%fX(pSet%nFires)
    pSet%fY(iFireToRemove) = pSet%fY(pSet%nFires)
    pSet%dx(iFireToRemove,:) = pSet%dx(pSet%nFires,:)
    pSet%dy(iFireToRemove,:) = pSet%dy(pSet%nFires,:)
    pSet%pFRP(iFireToRemove,:) = pSet%pFRP(pSet%nFires,:)
    pSet%pTA(iFireToRemove,:) = pSet%pTA(pSet%nFires,:)
    pSet%pT4(iFireToRemove,:) = pSet%pT4(pSet%nFires,:)
    pSet%pT4b(iFireToRemove,:) = pSet%pT4b(pSet%nFires,:)
    pSet%pT11(iFireToRemove,:) = pSet%pT11(pSet%nFires,:)
    pSet%pT11b(iFireToRemove,:) = pSet%pT11b(pSet%nFires,:)
    pSet%pMCE(iFireToRemove,:) = pSet%pMCE(pSet%nFires,:)
    pSet%pArea(iFireToRemove,:) = pSet%pArea(pSet%nFires,:)
    pSet%pHour(iFireToRemove,:) = pSet%pHour(pSet%nFires,:)
    pSet%day(iFireToRemove) = pSet%day(pSet%nFires)
    pSet%indLU(iFireToRemove) = pSet%indLU(pSet%nFires)
    pSet%nFires = pSet%nFires -1
    
  end subroutine remove_fire


  !**************************************************************************

  subroutine link_fire_src_to_species(species_list, fire_src)
    !
    ! Having the source species created, we should establish the shortcut links. 
    !
    implicit none

    ! Imported parameters
    type(Tsilam_fire_source), intent(inout) :: fire_src
    type(silam_species), dimension(:), pointer :: species_list

    !
    ! Linkage is actually just creation of the chemical adaptors.
    ! Again, species flaming and smouldering are identical in all senses, except for mass mean diameter,
    ! irrelevant here and ignored in all comparisons
    !
    call create_adaptor(fire_src%pFMD%species_flaming, species_list, fire_src%pFMD%adaptor)

  end subroutine link_fire_src_to_species


  ! *************************************************************************
  
  subroutine prepare_inject_fire_src(met_buf)
    !
    ! The subroutine prepares the private module pointers to the fields possibly requested 
    ! for injecting the point sources. As this is called only once before all point sources
    ! their actual needs are unknown. All we can do here is to check the meteo buffer for any 
    ! possibly needed quantity and set the pointer to all existing ones.
    !   
    implicit none

    ! Imported parameters
    type(Tfield_buffer), pointer :: met_buf

    ! Local variables
    integer, dimension(:), pointer :: met_q
    integer :: iQ, iTmp
    

    !Simple switch for Hybrid dispersion
    if (fu_leveltype(dispersion_vertical) .ne. layer_btw_2_height) then
      ifUseHybrid = .true.
    endif

    ! nullify the pointers 
    nullify(fldBVf)
    nullify(fldAblHeight)
    nullify(fldHeight)
    nullify(fldSrfPressure)

    ! Scan the meteo buffer   
    met_q => met_buf%buffer_quantities    
    do iQ = 1, size(met_q)
      if(met_q(iQ) == int_missing)exit
      if(fu_dimension(met_buf, iQ) == 4)then !4D

        select case(met_q(iQ))

          case(brunt_vaisala_freq_flag)
            fldBVf => met_buf%p4d(iQ)

          case(height_flag)
            fldHeight => met_buf%p4d(iQ)

          case default
            cycle
            
        end select
      else !2D

        select case(met_q(iQ))
          
          case(abl_height_m_flag)
            fldAblHeight => met_buf%p2d(iQ)
          
          case(ground_pressure_flag)
            fldSrfPressure => met_buf%p2d(iQ)
          
          case default
            cycle
            
        end select
      endif
    enddo
  end subroutine prepare_inject_fire_src


  !**************************************************************************

  subroutine inject_emission_euler_fire_src(fs, &
                                          & mapEmis, mapCoordX, mapCoordY, mapCoordZ, & ! Output
                                          & met_buf, & 
                                          & now, &      ! current time
                                          & timestep, & ! model time step
                                          & pHorizInterpMet2DispStruct, ifHorizInterp, &
                                          & pVertInterpMet2DispStruct, ifVertInterp, &
                                          & ifSpeciesMoment, &
                                          & fMassInjected)                              ! output
    !
    ! Computes the emission fields for fires. 
    ! The computation is split to two parts: 
    ! 1. Get the emission for each species. So far this is just a multipilication with 
    !    emission factor
    ! 2. Compute the injection profile and project it to dispersion vertical
    !
    ! This routine is to be called at each model time step but not inside the deepest cycle,
    ! therefore its efficiency is of moderate importance.
    !
    implicit none

    ! Imported parameters
    type(Tsilam_fire_source), target, intent(in) :: fs
    type(Tfield_buffer), pointer ::  met_buf  ! meteo and internal field buffers
    type(silja_time), intent(in) :: now           ! current time
    type(silja_interval), intent(in) :: timestep  ! model time step
    type(THorizInterpStruct), intent(in) :: pHorizInterpMet2DispStruct
    type(TVertInterpStruct), intent(in) :: pVertInterpMet2DispStruct
    logical, intent(in) :: ifHorizInterp, ifVertInterp, ifSpeciesMoment
    type(Tmass_map), intent(inout) :: mapEmis, mapCoordX, mapCoordY, mapCoordZ
    real(r8k), dimension(:), intent(inout) :: fMassInjected

    ! Local variable
    integer :: iSet

!call msg('Start injection:', fMassInjected(1))
    
    do iSet = 1, fs%nFRPdatasets
      if(fu_between_times(now, fs%first_day(iSet), fs%last_day(iSet)+one_day-one_minute, .true.))then
        call inject_emission_euler_FRP_set(fs%pFRPset(iSet), &
                                         & fs%pFMD, &
                                         & .not. (fs%first_day(iSet) == fs%last_day(iSet)), &
                                         & fs%id_Nbr)
!call msg('Injection,src:', fMassInjected(1), iSet)
      endif
    enddo
!call msg('End injection')

    CONTAINS

    !**************************************************************************

    subroutine inject_emission_euler_FRP_set(pSet, pFMD, ifManyDays, id_Nbr)
      !
      ! Computes the emission fields for fires. 
      ! The computation is split to two parts: 
      ! 1. Get the emission for each species. So far this is just a multipilication with 
      !    emission factor
      ! 2. Compute the injection profile and project it to dispersion vertical
      !
      ! This routine is to be called at each model time step but not inside the deepest cycle,
      ! therefore its efficiency is of moderate importance.
      !
      implicit none

      ! Imported parameters
      type(TFRP_dataset), intent(in) :: pSet
      type(Tsilam_fire_metadata), pointer :: pFMD
      logical, intent(in) :: ifManyDays
      integer, intent(in) :: id_Nbr

      ! Local variables
      integer :: iLev, ix, iy, iSrc, nLevs, iSpecies, ispeciesEmis, indHrPast, indHrFuture, iFire, &
                & iPlumePart
      real :: fWeightPast, fLevFraction, fCellTotal, timestep_sec, hour_UTC, hour_local, &
            & fPartFraction, fDiurnalVarTot, fDiurnalVarPerFire
      real, dimension(:), pointer :: fMassTimeCommon
      real, dimension(3) :: ptrCoord
      real :: plumeBottomHat, plumeTopHat, plumeBottom, plumeTop, overlapTop, overlapBottom, &
            & p_levBottom, p_levTop, p_srf

      !
      ! First of all, check that we have anything
      !
      fMassTimeCommon => fu_work_array()
      if(error)return
      fMassTimeCommon(1:200) = 0.0

      timestep_sec = abs(fu_sec(timestep))
      hour_UTC = fu_hour(now) + fu_min(now) / 60.

      !
      ! Now, scan the fires one-by-one
      !
      do iFire = 1, pSet%nFires
        if(ifManyDays)then
          if(now < pSet%day(iFire) .or. now > pSet%day(iFire) + one_day)cycle ! speedup
        endif
        !
        ! Since fire intensity changes fast, we take into account the fraction of an hour
        ! that corresponds to "nowLocal". FRP value corresponds to middle of an hour.
        ! Note that hours are 00..23, whereas indices are 1..24.
        !
        hour_local = modulo(hour_UTC + pSet%pLonGeo(iFire) / 15., 24.)
        fWeightPast = 1. - modulo(hour_local - int(hour_local) - 0.5, 1.)  ! centred at mid-hour, (0,1)
        indHrPast = int(modulo(hour_local-0.5, 24.)) + 1     ! (1:24)
        indHrFuture = modulo(indHrPast, 24) + 1           ! +1 and (1:24)
        fDiurnalVarTot = pFMD%pDiurnalVarTot(indHrPast, pSet%indLU(iFire)) * fWeightPast + &
                       & pFMD%pDiurnalVarTot(indHrFuture, pSet%indLU(iFire)) * (1.-fWeightPast)
        fDiurnalVarPerFire = pFMD%pDiurnalVarPerFire(indHrPast, pSet%indLU(iFire))*fWeightPast + &
                           & pFMD%pDiurnalVarPerFire(indHrFuture, pSet%indLU(iFire))*(1.-fWeightPast)
        !
        ! Having all input calculated, get the flaming-smouldering mixture
        !
        do iSpecies = 1, pFMD%nSpecies
          fMassTimeCommon(iSpecies) = &
               & fu_emission_weighted_flam_smld(pSet%pFRP(iFire,pSet%indFRPDaily_tot), &
                                              & timestep_sec, &
                                              & pFMD%pEmsFactors(:,pSet%indLU(iFire),iSpecies), &
                                              & fDiurnalVarTot)
        enddo
        ix = nint(pSet%fX(iFire))
        iy = nint(pSet%fY(iFire))
        !
        ! Plume is described via its top and bottom 
        !
        call fire_plume_vertical_profile(&
                              & pSet%pFRP(iFire,pSet%indFRPDaily_perFire) * fDiurnalVarPerFire, &
                              & ix, iy, met_buf, &
                              & pHorizInterpMet2DispStruct, ifHorizInterp, &
                              & pVertInterpMet2DispStruct, ifVertInterp, &
                              & ifOneStepHeightProcedure, &
                              & plumeBottomHat, plumeTopHat, p_srf, ifUseHybrid)
        if(error)return

        !Prepare pressure
        if (ifUseHybrid) then
          p_levTop = a_half_disp(1) + b_half_disp(1) * p_srf
        endif
        !
        ! We consider the plume stem and hat separately
        !
        do iPlumePart = iPlumeStem, iPlumeHat
          if(iPlumePart == iPlumeStem)then
            if(fStemMassFraction < 1.e-5) cycle  ! if stem is empty, do not waste time
            fPartFraction = fStemMassFraction
            plumeBottom = min(1.0, plumeBottomHat/2.)                 ! 1 m
            plumeTop = plumeBottomHat
          else
            fPartFraction = 1. - fStemMassFraction
            plumeTop = plumeTopHat
            plumeBottom = plumeBottomHat
          endif
          !
          ! We do the injection layer by layer
          !
          do iLev = 1, nz_dispersion
            !
            ! First do the check for the overlap: speed-up
            !
!call msg('Lev:',iLev)
            if (ifUseHybrid) then
               p_levBottom = p_levTop
               p_levTop = a_half_disp(iLev+1)+b_half_disp(iLev+1)*p_srf  
               overlapBottom = min(plumeBottom, p_levBottom) !pressure increases downwards!
               overlapTop = max(plumeTop, p_levTop)
  
               if (overlapBottom <= overlapTop) cycle
               ptrCoord(3) =  (overlapBottom + overLapTop)/2  ! Center of slab
               ptrCoord(3) =  ptrCoord(3) - (p_levTop + p_levBottom)/2 !relative to cell center 
               ptrCoord(3) = - ptrCoord(3) / (p_levTop - p_levBottom) !+-0.5,
                                                                      !positive -up
            else
               overlapBottom = max(plumeBottom, disp_layer_top_m(iLev-1))
               overlapTop = min(plumeTop, disp_layer_top_m(iLev))
  
               if (overlapBottom >= overlapTop) cycle
               ptrCoord(3) = (overlapBottom + overLapTop) / 2
               ptrCoord(3) = (ptrCoord(3) - & 
                          & (disp_layer_top_m(iLev) + disp_layer_top_m(iLev-1))/2) / &  !relative to cell center
                          & (disp_layer_top_m(iLev) - disp_layer_top_m(iLev-1))
            endif
            fLevFraction = (overlapTop - overlapBottom)/(plumeTop - plumeBottom)
            if(abs(ptrCoord(3)) > 0.5)then
              call msg('Relative centre of mass position is strange in layer:',iLev, ptrCoord(3))
              call msg('Plume top, bottom:', plumeTop, plumeBottom)
              if (ifUseHybrid) then
                call msg('Dispersion hybrid layer bottom, top:',p_levBottom, p_levBottom)
              else
                call msg('Dispersion layer top [m]:',disp_layer_top_m(iLev), disp_layer_top_m(iLev-1))
              endif
              call set_error('Wrong calculated vertical plume positon','inject_emission_euler_FRP_set')
              if(abs(ptrCoord(3)) - 0.5 < 0.0001)then
                call unset_error('inject_emission_euler_FRP_set')
                ptrCoord(3)  = sign(0.498, ptrCoord(3))
              else
                return
              endif
            endif

!call msg('ilev, vertical fraction:', ilev, fLevFraction)
!            ptrCoord(1) = max(-0.4999, min(0.4999, pSet%fX(iFire) - ix))
!            ptrCoord(2) = max(-0.4999, min(0.4999, pSet%fY(iFire) - iy))
            ptrCoord(1) = pSet%fX(iFire) - ix
            ptrCoord(2) = pSet%fY(iFire) - iy
            if(abs(ptrCoord(1)) > 0.4999 .or. abs(ptrCoord(2)) > 0.4999) call check_ptrCoord(ptrCoord)
            fCellTotal = 0.0
            if(error)return
            !
            ! Emit the mass species by species
            !
!if(iFire == 10)call msg('fs%fDescr2SpeciesUnit(1,1), fLevFraction, fPartFraction',(/fs%fDescr2SpeciesUnit(1,1), fLevFraction, fPartFraction/))

            do iSpecies = 1, pFMD%nSpecies
              iSpeciesEmis = pFMD%adaptor%iSp(iSpecies)
              mapEmis%arM(iSpeciesEmis, id_Nbr, iLev,ix,iy) = & 
                                        & mapEmis%arM(iSpeciesEmis, id_Nbr, iLev,ix, iy) + &
                                        & fMassTimeCommon(iSpecies) * fLevFraction * fPartFraction
              fCellTotal = fCellTotal + fMassTimeCommon(iSpecies) * fLevFraction * fPartFraction
              fMassInjected(iSpeciesEmis) = fMassInjected(iSpeciesEmis) + &
                                          & fMassTimeCommon(iSpecies) * fLevFraction * fPartFraction
              if (ifSpeciesMoment) then
                mapCoordX%arm(iSpeciesEmis,id_nbr,ilev,ix,iy) = &
                                      & mapCoordX%arm(iSpeciesEmis, id_nbr, ilev, ix, iy) + &
                                      & fMassTimeCommon(iSpecies) * fLevFraction * fPartFraction * ptrCoord(1)
                mapCoordY%arm(iSpeciesEmis, id_nbr, ilev, ix, iy) = &
                                      & mapCoordY%arm(iSpeciesEmis, id_nbr, ilev, ix, iy) + &
                                      & fMassTimeCommon(iSpecies) * fLevFraction * fPartFraction * ptrCoord(2)
                mapCoordZ%arm(iSpeciesEmis, id_nbr, ilev, ix, iy) = &
                                      & mapCoordZ%arm(iSpeciesEmis, id_nbr, ilev, ix, iy) + &
                                      & fMassTimeCommon(iSpecies) * fLevFraction * fPartFraction * ptrCoord(3)
              end if
            end do  ! species

            ! If we use bulk moment, add it of the injected masses and the total injected mass
            !
            if (.not. ifSpeciesMoment) then
              mapCoordx%arM(1,id_nbr, iLev, ix, iy) = &
                           & mapCoordx%arM(1,id_nbr, iLev, ix, iy) + ptrCoord(1) * fCellTotal
              mapCoordy%arM(1,id_nbr, iLev, ix, iy) = &
                           & mapCoordy%arM(1,id_nbr, iLev, ix, iy) + ptrCoord(2) * fCellTotal
              mapCoordz%arM(1,id_nbr, iLev, ix, iy) = &
                           & mapCoordz%arM(1,id_nbr, iLev, ix, iy) + ptrCoord(3) * fCellTotal
            end if

!do iSpeciesEmis = 1, mapEmis%nSpecies
!if(abs(mapCoordX%arM(iSpeciesEmis,id_nbr,iLev,ix,iy) / &
!& mapEmis%arM(iSpeciesEmis,id_nbr,iLev,ix,iy)) >= 0.501)then
!call set_error('Strange X centre of mass','inject_emission_euler_FRP_set')
!call unset_error
!endif
!if(abs(mapCoordY%arM(iSpeciesEmis,id_nbr,iLev,ix,iy) / &
!& mapEmis%arM(iSpeciesEmis,id_nbr,iLev,ix,iy)) >= 0.501)then
!call set_error('Strange Y centre of mass','inject_emission_euler_FRP_set')
!endif
!if(abs(mapCoordZ%arM(iSpeciesEmis,id_nbr,iLev,ix,iy) / &
!& mapEmis%arM(iSpeciesEmis,id_nbr,iLev,ix,iy)) >= 0.501)then
!call set_error('Strange Z centre of mass','inject_emission_euler_FRP_set')
!endif
!if(error)then
!call msg('x,y,lev,src,species:',(/ix,iy,iLev,id_nbr,iSpeciesEmis/))
!call msg('mass, x,y,z moments',(/mapEmis%arM(iSpeciesEmis,id_nbr,iLev,ix,iy), &
!& mapCoordX%arM(iSpeciesEmis,id_nbr,iLev,ix,iy), &
!& mapCoordY%arM(iSpeciesEmis,id_nbr,iLev,ix,iy), &
!& mapCoordZ%arM(iSpeciesEmis,id_nbr,iLev,ix,iy)/))
!endif
!enddo
            mapEmis%ifColumnValid(id_nbr, ix, iy) = .true.
            mapEmis%ifGridValid(ilev, id_nbr) = .true.

          end do  ! iLev dispersion
        end do  ! plume stem, hat
      end do  ! iFire

      call free_work_array(fMassTimeCommon)

    end subroutine inject_emission_euler_FRP_set

    !================================================================================

    subroutine check_ptrCoord(ptrCoord)
      ! check that the apparent error is not a big deal
      implicit none
      real, dimension(:), intent(inout) :: ptrCoord
      integer :: i
      do i = 1, 2
        if(ptrCoord(i) > 0.4999)then
          if(ptrCoord(i) > 0.5001)then
            call set_error('Strange ptrCoord(' + fu_str(i) + '):' + fu_str(ptrCoord(i)),'inject_emission_euler_FRP_set')
            return
          else
            ptrCoord(i) = 0.4999   ! just numerics
          endif
        endif
        if(ptrCoord(i) < -0.4999)then
          if(ptrCoord(i) < -0.5001)then
            call set_error('Strange ptrCoord(' + fu_str(i) + '):' + fu_str(ptrCoord(i)),'inject_emission_euler_FRP_set')
            return
          else
            ptrCoord(i) = -0.4999   ! just numerics
          endif
        endif
      end do
    end subroutine check_ptrCoord
          
    !=======================================================================================

    subroutine fire_plume_vertical_profile(FRP, ix, iy, &
                                         & met_buf, &
                                         & pHorizInterpMet2DispStruct, ifMetHorizInterp, &
                                         & interpCoefMeteo2DispVert, ifMetVertInterp, &
                                         & ifOneStepProcedure, &
                                         & fLevBottom, fLevTop, p_srf, ifPressure)
      implicit none

      ! Imported parameters
      real, intent(in) :: FRP
      integer,  intent(in) :: ix,iy
      type(Tfield_buffer), pointer ::  met_buf
      type(THorizInterpStruct), intent(in) :: pHorizInterpMet2DispStruct
      type(TVertInterpStruct), intent(in) :: interpCoefMeteo2DispVert
      logical, intent(in) :: ifMetHorizInterp, ifMetVertInterp, ifOneStepProcedure
      real, intent(out) :: fLevBottom, fLevTop ! plume bounds
      real, intent(out) :: p_srf !  return only when ifPressure == .true.
      logical, intent(in) :: ifPressure ! If true -- fLevBottom, fLevTop come as pressure
                                     ! otherwise -- as height
      integer :: iLev
      real :: fTmp
      real, dimension(1:max_levels) :: height

      ! Local parameters
      real, parameter :: alpha_full=0.24, alpha_step_1 = 0.15, alpha_step_2 = 0.93, &
                       & beta_full=169,   beta_step_1 = 102.,  beta_step_2 = 298., &
                       & gamma_full=0.35, gamma_step_1 = 0.49, gamma_step_2 = 0.13, &
                       & delta_full=254.,                      delta_step_2 = 259., &
                       & FRP_scale=1.e-6

      ! Local variables
      real :: fABL_height, fBruntVaisFreq
      
      if (.not.(associated(fldBVf) .and. associated(fldAblHeight)))then
        call set_error('Not all meteo found','fire_plume_vertical_profile')
        return
      endif


      fABL_height = fu_get_value(fldAblHeight, nx_meteo, ix, iy, &
                       & met_buf%weight_past, pHorizInterpMet2DispStruct, ifMetHorizInterp)

      if(fABL_height < 10. .or. fABL_height > 6000.)then
        call msg('Funny ABL height at (ix,iy):' + &
               & fu_str(ix) + ',' + fu_str(iy), fABL_height)
      !
      ! Do we need this rubish? Should be checked when ABL is prepared.....
        fABL_height = (fu_get_value(fldAblHeight, nx_meteo, min(ix+1, nx_dispersion), iy, &
                            & met_buf%weight_past, pHorizInterpMet2DispStruct, ifMetHorizInterp) + &
                     & fu_get_value(fldAblHeight, nx_meteo, ix, min(iy+1,ny_dispersion), &
                            & met_buf%weight_past, pHorizInterpMet2DispStruct, ifMetHorizInterp) + &
                     & fu_get_value(fldAblHeight, nx_meteo, max(1,ix-1), iy, &
                            & met_buf%weight_past, pHorizInterpMet2DispStruct, ifMetHorizInterp) + &
                     & fu_get_value(fldAblHeight, nx_meteo, ix, max(1,iy-1), &
                            & met_buf%weight_past, pHorizInterpMet2DispStruct, ifMetHorizInterp)) / 4.0
        fABL_height = max(10.,min(6000.,fABL_height))
        call msg('Funny ABL height, force mean regional:',fABL_height)
      endif

      ! fBV above the ABL 
      fTmp = min(fABL_height * 2.,  fABL_height + 1000.)
      do iLev=1,nz_meteo-1 ! Should be still within meteo vert on finish!!
        height(iLev) = fu_get_value(fldHeight, nx_meteo, ix, iy, iLev,&
                          & met_buf%weight_past, pHorizInterpMet2DispStruct, interpCoefMeteo2DispVert,&
                          & ifMetHorizInterp, .false.)
 !                 call msg("Lev,height"+ fu_int2str(iLev), height(iLev), a_met(iLev)+b_met(iLev)*p_srf)
        if (height(iLev) >= fTmp) exit
      enddo
      if (iLev > nz_meteo - 1) then
        call msg("Fire source:  ABL above met ablh = "//trim(fu_str(fABL_height))//": ix, iy, iLev", (/ ix, iy, iLev/))
      endif
      ! NB!  No vertical interpolation here! Level is found already!
      fBruntVaisFreq = fu_get_value(fldBVf, nx_meteo, ix, iy, iLev, & 
                          & met_buf%weight_past, pHorizInterpMet2DispStruct, interpCoefMeteo2DispVert,&
                          & ifMetHorizInterp, .false.)

      if(ifOneStepProcedure)then
        !
        ! One-step procedure means a unified formula
        !
        fLevTop = alpha_full * fABL_height + &
                & beta_full * ((FRP * FRP_scale) **gamma_full) * exp(-delta_full * fBruntVaisFreq)
        fLevBottom = fLevTop / 3.
      else
        !
        ! For two-steps, first compute the separation value, which is to be compared with ABL
        ! Then - either FT- or unified formula (the later is also OK for ABL)
        !
        fLevTop = alpha_step_1 * fABL_height + beta_step_1 * ((FRP * FRP_scale) **gamma_step_1)
        if(fLevTop > fABL_height)then
          fLevTop = alpha_step_2 * fABL_height + &
                  & beta_step_2 * ((FRP * FRP_scale) **gamma_step_2) * &
                  & exp(-delta_step_2 * fBruntVaisFreq)
          fLevBottom = fLevTop / 2.
        else
          fLevTop = alpha_full * fABL_height + &
                  & beta_full * ((FRP * FRP_scale) **gamma_full) * exp(-delta_full * fBruntVaisFreq)
          fLevBottom = fLevTop / 3.
        endif
      endif

      p_srf = real_missing
      if (ifPressure) then !Recalculate heights to pressures
        ! call msg("Converting to pressure",fLevBottom,fLevTop)

        p_srf = fu_get_value(fldSrfPressure, nx_meteo, ix, iy, &
                       & met_buf%weight_past, pHorizInterpMet2DispStruct, ifMetHorizInterp)
        ! Get missing heights to the column

        do while (height(iLev) < fLevTop) 
          iLev=iLev+1
          if (iLev >= nz_meteo) then
            call msg("Attempt to emit above the meteo grid 1")
            call set_error("Fire source got crazy....", "fire_plume_vertical_profile")
            call unset_error("fire_plume_vertical_profile")
            iLev = nz_meteo
            height(nz_meteo) = fu_get_value(fldHeight, nx_meteo, ix, iy, nz_meteo,&
                         & met_buf%weight_past, pHorizInterpMet2DispStruct, interpCoefMeteo2DispVert,&
                         & ifMetHorizInterp, .false.)
            exit
          endif
          height(iLev) = fu_get_value(fldHeight, nx_meteo, ix, iy, iLev,&
                       & met_buf%weight_past, pHorizInterpMet2DispStruct, interpCoefMeteo2DispVert,&
                       & ifMetHorizInterp, .false.)
  !       call msg("Lev,height"+ fu_int2str(iLev), height(iLev), a_met(iLev)+b_met(iLev)*p_srf)
        enddo   ! while height < fLevTop
        iLev=2
        do while (height(iLev) < fLevBottom) 
          iLev = iLev + 1
          if (iLev >= nz_meteo) then
            call msg("Attempt to emit above the meteo grid 2")
            call set_error("Fire source got crazy....", "fire_plume_vertical_profile")
            call unset_error("fire_plume_vertical_profile")
            iLev = nz_meteo-1   ! bottom must be lower than the top
            exit
          endif
        enddo
        fTmp = (height(iLev) - fLevBottom)/ (height(iLev) - height(iLev-1))
        fLevBottom = (a_met(iLev)+b_met(iLev)*p_srf) * (1. - fTmp) + &
                   & (a_met(iLev-1)+b_met(iLev-1)*p_srf) * fTmp

        do while (height(iLev) < fLevTop) 
          iLev=iLev+1
          if (iLev >= nz_meteo) then
            call msg("Attempt to emit above the meteo grid 3")
            call set_error("Fire source got crazy....", "fire_plume_vertical_profile")
            call unset_error("fire_plume_vertical_profile")
            iLev = nz_meteo
            exit
          endif
        enddo
        fTmp = max((height(iLev) - fLevTop)/ (height(iLev) - height(iLev-1)), 0.0)

        !FIXME! Bullshit dispersion levels should be here!!!
        fLevTop = (a_met(iLev)+b_met(iLev)*p_srf) * (1. - fTmp) + &
                & (a_met(iLev-1)+b_met(iLev-1)*p_srf) * fTmp
!       call msg("Converted to pressure",fLevBottom,fLevTop)
      endif

    end subroutine fire_plume_vertical_profile

  end subroutine inject_emission_euler_fire_src


  !**********************************************************************  
    
  real function fu_emission_weighted_flam_smld(fFRP, timestep_sec, pEmisFactors, fDiurnalVarTot)
    !
    ! Mixes the flaming and smouldering emission fluxes. 
    ! This is the actual place, where the contribution of flaming and smouldering is decided
    ! Emission factors give absolute rate of cocktails. Note: [kg/J]; fire is [W]
    !
    implicit none
    
    ! Imported parameters
    real, intent(in) :: fFRP, timestep_sec, fDiurnalVarTot
    real, dimension(:), intent(in) :: pEmisFactors  ! (flaming/smouldering)

    ! Local variables
    integer :: iType
    real :: flaming_fraction

    !
    ! Here we say a horrible thing: we assume that a fraction of flaming is proportional 
    ! to the diurnal variation. Roughly speaking, during night it will be almost all smouldering,
    ! almost all flaming at the day peak. Specific fractions are completely from the blue sky
    !
    flaming_fraction = 1. - exp(-fDiurnalVarTot)            ! flaming fraction ->1 for peak

    fu_emission_weighted_flam_smld = 0.
    do iType = flaming, smouldering
      fu_emission_weighted_flam_smld = fu_emission_weighted_flam_smld + &
                                     & fFRP * timestep_sec * fDiurnalVarTot * pEmisFactors(iType) * &
                                     & (flaming_fraction * (smouldering - iType) + &     ! flaming = 1
                                      & (1. - flaming_fraction) * (iType - flaming))     ! smouldering = 2
    end do  ! flaming:smouldering

  end function fu_emission_weighted_flam_smld

        
  !*****************************************************************

  integer function fu_source_id_nbr_of_fire_src(fire_src)
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
    type(Tsilam_fire_source), intent(in) :: fire_src

    ! Stupidity check
    if(.not.(fire_src%defined == silja_false))then
      call set_error('Undefined source given','fu_source_id_nbr_of_fire_src')
      return
    endif
    fu_source_id_nbr_of_fire_src = fire_src%id_nbr

  end function fu_source_id_nbr_of_fire_src



  !*************************************************************************

  integer function fu_source_nbr_of_fire_src(fire_src)
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
    type(Tsilam_fire_source), intent(in) :: fire_src

    ! Stupidity check
    if(.not.(fire_src%defined == silja_false))then
      fu_source_nbr_of_fire_src = fire_src%src_nbr
    else
      fu_source_nbr_of_fire_src = int_missing
      call set_error('Undefined source given','fu_source_nbr_of_fire_src')
      return
    endif

  end function fu_source_nbr_of_fire_src


  !*****************************************************************

  subroutine tot_amt_species_unit_fire_src(fire_src, &      ! mandatory, input
                                         & species, nSpecies, amounts, & ! mandatory, output
                                         & start, duration)   ! optional, input
    !
    ! Returns the amount of the released material IN SPECIES UNIT starting from 
    ! start during the duration time interval. Species must be initialised by that moment
    !
    ! Will scan fire by fire. Since each operates one day, no time variation is needed, unless
    ! someone wants exact period. Then will have to see what to do. Note the problem
    ! with local time in that case.
    !
    ! Units: SI basic and derived ones
    !
    implicit none

    ! Imported parameters with intent IN
    type(Tsilam_fire_source), intent(in) :: fire_src
    type(silam_species), dimension(:), pointer :: species
    integer, intent(out) :: nSpecies
    real, dimension(:), pointer :: amounts
    type(silja_time), intent(in), optional :: start
    type(silja_interval), intent(in), optional :: duration

    ! Local variables
    integer :: iSet
    logical :: ifLimitedInTime

    do iSet = 1, fire_src%nFRPdatasets
      if(present(start))then
        if(start > fire_src%last_day(iSet) + one_day) cycle
        ifLimitedInTime = start > fire_src%first_day(iSet)
      else
        ifLimitedInTime = .false.
      endif
      if(present(duration))then
        if(start+duration < fire_src%first_day(iSet))cycle
        ifLimitedInTime = start+duration < fire_src%last_day(iSet)
      else
        ifLimitedInTime = .false.
      endif
      
      nSpecies = 0

      call addSpecies(species, nSpecies, fire_src%pFMD%species_flaming, fire_src%pFMD%nSpecies)
      call addSpecies(species, nSpecies, fire_src%pFMD%species_smouldering, fire_src%pFMD%nSpecies)
    
      amounts(1:fire_src%pFMD%nSpecies) = 0.0

      call tot_amt_species_unit_FRP_set(fire_src%pFRPset(iSet), ifLimitedInTime, fire_src%pFMD)
    
    end do


    CONTAINS
    
    subroutine tot_amt_species_unit_FRP_set(pSet, ifLimitedInTime, pFMD)
      !
      ! Gets the amounts for one subset 
      !
      implicit none
      
      ! Imported parameters
      type(TFRP_dataset), intent(in) :: pSet
      logical, intent(in) :: ifLimitedInTime
      type(Tsilam_fire_metadata), pointer :: pFMD

      ! Local variables
      integer :: iSpecies, iHrLocal, iHr, iFire, iSet
      type(silja_time) :: now, end_time
      real :: seconds
      !
      ! Start from obtaining the amounts in the descriptor units
      !
      if(ifLimitedInTime)then
        !
        ! Limited time: have to go hour by hour accounting for variation
        !
        end_time = start + duration
        do iFire = 1, pSet%nFires
          if(end_time < pSet%day(iFire))cycle
          if(start > pSet%day(iFire)+one_day)cycle
          now = pSet%day(iFire)
          do iHr = 1, 24
            seconds = fu_sec((fu_earliest_time((/end_time, now+one_hour/)) - &
                            & fu_latest_time((/start, now/))))
            if(seconds <= 0.0)cycle
            iHrLocal = iHr + pSet%pLonGeo(iFire) / 15.
            if(iHrLocal < 0)iHrLocal = iHrLocal + 24
            if(iHrLocal >= 24)iHrLocal = iHrLocal - 24
            do iSpecies = 1, pFMD%nSpecies
              amounts(iSpecies) = amounts(iSpecies) + &
                            & fu_emission_weighted_flam_smld(pSet%pFRP(iFire,pSet%indFRPDaily_tot), &
                                                 & seconds, &
                                                 & pFMD%pEmsFactors(:,pSet%indLU(iFire),iSpecies), &
                                                 & pFMD%pDiurnalVarTot(iHrLocal+1, pSet%indLU(iFire)))
            end do   ! iSpecies
          end do  ! iHr
        end do  ! iFire
      else
        !
        ! Whole day for this fire. Take the approximate mixture between the flaming and smouldering
        !
        do iFire = 1, pSet%nFires
          do iSpecies = 1, pFMD%nSpecies
            amounts(iSpecies) = amounts(iSpecies) + &
                              & pSet%pFRP(iFire,pSet%indFRPDaily_tot) * 86400 * &
                                 & (pFMD%pEmsFactors(flaming,pSet%indLU(iFire),iSpecies) * &
                                  & pFMD%flaming_fraction_roughly(pSet%indLU(iFire),iSpecies) + &
                                 & (pFMD%pEmsFactors(smouldering,pSet%indLU(iFire),iSpecies) * &
                                  & (1.-pFMD%flaming_fraction_roughly(pSet%indLU(iFire),iSpecies))))
          enddo   ! nSpecies
        end do  ! nFires
      endif  ! if Limited in time
      if(error)return

    end subroutine tot_amt_species_unit_FRP_set

  end subroutine tot_amt_species_unit_fire_src

  
  !***********************************************************************************************
  
  subroutine fires_flux4mode(aerMode, &   ! definition of the spectrum band - from the run setup
                           & iSpectrumType, &  ! description of the fire-emitted aerosol
                           & fNbrFlux, fVolFlux, fMassMeanD)   ! output
    !
    ! Computes the number emission flux (shape function) of the fire-emitted PM
    ! for the pre-defined bin of the size spectrum of the run.
    !
    implicit none

    ! Imported parameters
    type(Taerosol_mode), intent(inout) :: aerMode
    integer, intent(in) :: iSpectrumType
    real, dimension(2), intent(out) :: fNbrFlux, fVolFlux, fMassMeanD

    ! Local variables
    real :: fMinD, fMaxD, fD, fIntegrStep, fTmp, factor, fFluxDens, fFluxDens_prev, fTmpVol, fTmpNbr
    real, dimension(2) :: fTotalNumber, fTotalVolume
    type(Taerosol_mode) :: aerModeTmp
    integer :: iMode, iType
    !
    ! Parameters for all distribution descriptions
    !
    real, parameter :: fAbsoluteMin = 2.0e-9  ! 3 nm
    real, parameter :: fAbsoluteMax = 1.0e-4  ! 100 um
    !
    ! Parameters for 3-mode spectrum representation (Virkkula et al, ACP 2014) and 
    ! Janhall et al, ACP 2010, and Chubarova et al, Atm.Meas.Technicue 2012.
    ! Ideas: substantial fraction of mass is in mid-size mode, some 150-250 nm. It
    ! gets >70% of the total mass. Coarse fraction, a few um, receives only 25%.
    ! Secondly, smouldering means smaller particles: they show up from the fine side
    ! of the spectrum.
    ! With some modifications, below numbers follow Virkkula - but with adaptations as
    ! described in Notebook 12, p.13-15.
    !
    real, dimension(3,2), parameter :: nbr_fract_3_modes = reshape ((/0.0299, 0.97,   0.0001, & ! flames
                                                                    & 0.1599, 0.84,   0.0001/), &  ! smolder
                                                                  & (/3, 2/) )
    real, dimension(3), parameter :: nbr_mean_diam_3_modes=(/3.2e-8, 8.1e-8, 1.45e-7/)
    real, dimension(3), parameter :: nbr_stdev_3_modes =    (/1.32,   1.58,    2.6/)

!    type(silam_sp) :: sp, sp2

!sp%sp => fu_work_string()
!sp2%sp => fu_work_string()
!if(error)return

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
      fMassMeanD = 0.5*(fMinD + fMaxD)

    else
      !
      ! The range given has useful part
      !
      if(fMinD < fAbsoluteMin) fMinD = fAbsoluteMin
      if(fMaxD > fAbsoluteMax) fMaxD = fAbsoluteMax

call msg('Bin requested with the borders, um:', fMinD*1e6, fMaxD*1e6)

      select case(iSpectrumType)

        case(lognorm_3_modes_fires_flag)
          !
          ! After Weinzierl dissertation: 3 lognormal modes
          !
          fNbrFlux = 0.0
          fVolFlux = 0.0
          fMassMeanD = 0.0
          fTotalNumber = 0.0
          fTotalVolume = 0.0
          
          do iType = flaming, smouldering
            !
            ! Each type of the fire emits into 3 modes, just with different fractionation
            !
if(iType == flaming)call msg('flaming')
if(iType == smouldering)call msg('smouldering')
!sp%sp = 'min_um:' + fu_str(fMinD*1e6) + ', max_um:' + fu_str(fMaxD*1e6) + ', number from 3 modes:'
!sp2%sp = ',   volume from 3 fire modes:'

            do iMode = 1, 3
              call set_aerosol_mode(aerModeTmp, '', &                         ! mode, chNm, 
                                  & nbr_mean_diam_3_modes(iMode), &           ! fp1
                                  & nbr_stdev_3_modes(iMode), &               ! fp2
                                  & nbr_mean_diam_3_modes(iMode), 2.6e3, &   ! mass_mean_d, dens, 
                                  & lognormal_flag, 1)                       ! distr_type, solubil
              if(error)return
              !
              ! Calculate the total-spectrum contribution. fu_integrate_* is a fraction of the number- 
              ! or volume distributions that fall into the given range
              !
              fTotalNumber(iType) = fTotalNumber(iType) + &
                                  & nbr_fract_3_modes(iMode,iType) * &
                                  & fu_integrate_number(fAbsoluteMin, fAbsoluteMax, aerModeTmp)
              fTotalVolume(iType) = fTotalVolume(iType) + &
                                  & nbr_fract_3_modes(iMode,iType) * &
                                  & fu_integrate_volume(fAbsoluteMin, fAbsoluteMax, aerModeTmp) !* &
!                                  & (fu_mass_mean_diam(fAbsoluteMin, fAbsoluteMax, aerModeTmp))**3
call msg('Cumulative total number and volume, full range:', fTotalNumber(iType), fTotalVolume(iType))
              !
              ! Number-flux fraction is the fraction of the mode falling in the output mode
              !
              fTmp = nbr_fract_3_modes(iMode,iType) * fu_integrate_number(fMinD, fMaxD, aerModeTmp)
              fNbrFlux(iType) = fNbrFlux(iType) + fTmp
!sp%sp = sp%sp + ',' + fu_str(real(fTmp))
              !
              ! Volume-flux fraction is the fraction of the mode in the volume
              !
              fTmp = nbr_fract_3_modes(iMode,iType) * fu_integrate_volume(fMinD, fMaxD, aerModeTmp) !* &
!                   & (fu_mass_mean_diam(fMinD, fMaxD, aerModeTmp))**3
!sp2%sp = sp2%sp + ',' + fu_str(real(fTmp))
              fVolFlux(iType) = fVolFlux(iType) + fTmp
!call msg('fMinD, fMaxD, mass_mean_d:',(/fMinD, fMaxD, fu_mass_mean_diam(fMinD, fMaxD, aerModeTmp)/))
              fMassMeanD(iType) = fMassMeanD(iType) + fTmp * fu_mass_mean_diam(fMinD, fMaxD, aerModeTmp)
            end do  ! modes 1..3

!call msg(sp%sp + sp2%sp)

            fMassMeanD(iType) = fMassMeanD(iType) / fVolFlux(iType)  ! get the mean diameter
            fNbrFlux(iType) = fNbrFlux(iType) / fTotalNumber(iType)  ! normalise the nbr fraction
            fVolFlux(iType) = fVolFlux(iType) / fTotalVolume(iType)  ! normalise the volume fraction

call msg('Total number and volume fluxes:', fTotalNumber(iType), fTotalVolume(iType))
call msg('Normalised number and volume flux for the bin:', fNbrFlux(iType), fVolFlux(iType))
call msg('Mass mean diameter:',fMassMeanD(iType))

          end do ! flaming, smouldering
call msg('')
call msg('')

        case(lognorm_3_modes_fires_num_flag)
          !
          ! Same as above but numerical integration. Analytical expressoin is wrong
          !
          factor = (fAbsoluteMax / fAbsoluteMin) ** 0.0001 - 1.

          fNbrFlux = 0.
          fVolFlux = 0.
          fTotalNumber = 0.0
          fTotalVolume = 0.0
          fMassMeanD = 0.
          
          do iType = flaming, smouldering

            fD = fAbsoluteMin  !fMinD
            fFluxDens = 0.0
            do iMode = 1, 3
              fFluxDens = fFluxDens + nbr_fract_3_modes(iMode,iType) * &
                      & fu_lognorm_density(fD, nbr_mean_diam_3_modes(iMode), nbr_stdev_3_modes(iMode))
            end do

            do while(fD < fAbsoluteMax) !fMaxD)

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
              do iMode = 1, 3
                fFluxDens = fFluxDens + nbr_fract_3_modes(iMode, iType) * &
                     & fu_lognorm_density(fD, nbr_mean_diam_3_modes(iMode), nbr_stdev_3_modes(iMode))
              end do
              !
              ! Integrate the number and volume flux. 
              !
              fTmpNbr = fIntegrStep * 0.5 * (fFluxDens_prev + fFluxDens)  ! for number flux
              fTmpVol = fIntegrStep * 0.5 * 0.5235987756 * (fFluxDens_prev * (fD-fIntegrStep)**3 + &
                                                          & fFluxDens * fD**3)
              fTotalNumber(iType) = fTotalNumber(iType) + fTmpNbr
              fTotalVolume(iType) = fTotalVolume(iType) + fTmpVol
              if(fD >= fMinD .and. fD <= fMaxD)then
                fNbrFlux(iType) = fNbrFlux(iType) + fTmpNbr
                fVolFlux(iType) = fVolFlux(iType) + fTmpVol
                fMassMeanD(iType) = fMassMeanD(iType) + fTmpVol * (fD-0.5*fIntegrStep)
              endif

            end do  ! cycle over diameter
!call msg('lognorm_4_modes_dust_numeric_flag: fD range min-max, totalNbr, totalVol, binNbr, binVol', &
!    & (/fMinD, fMaxD, fTotalNumber, fTotalVolume, fNbrFlux, fVolFlux/))
          
            fMassMeanD(iType) = fMassMeanD(iType) / fVolFlux(iType)
            if((.not. fMassMeanD(iType) >= fMinD) .or. (.not. fMassMeanD(iType) <= fMaxD))then
              call msg('Wrong 3-mode mean diameter. Min, max, mean:', (/fMinD, fMaxD, fMassMeanD(iType)/))
              call msg_warning('Resetting mean diameter, fire type:'+fu_str(iType),'fire_flux4mode')
              fMassMeanD(iType) = sqrt(fMinD * fMaxD)
            endif
            fNbrFlux(iType) = fNbrFlux(iType) / fTotalNumber(iType)
            fVolFlux(iType) = fVolFlux(iType) / fTotalVolume(iType)

          end do ! flaming, smouldaering    
          
          
          
          
          
          
          
          
        case default
          call set_error('Unknown type of dust spectrum:' + fu_str(iSPectrumTYpe), 'fires_flux4mode')
      end select
      
!call free_work_array(sp%sp)
!call free_work_array(sp2%sp)
      
    endif  ! useful part of the range

  end subroutine fires_flux4mode

  
  !=========================================================================
  function fu_fires_source_name(fire_src)result(chNm)
    implicit none
    type(Tsilam_fire_source), intent(in) :: fire_src
    character(len=clen) :: chNm
    chNm = fire_src%src_nm
  end function fu_fires_source_name

  !=========================================================================
  function fu_fires_sector_name(fire_src)result(chNm)
    implicit none
    type(Tsilam_fire_source), intent(in) :: fire_src
    character(len=clen) :: chNm
    chNm = fire_src%sector_nm
  end function fu_fires_sector_name

  !=========================================================================
  function fu_start_time_fire_src(fire_src, iSet)
    implicit none
    type(Tsilam_fire_source), intent(in) :: fire_src
    integer, intent(in), optional :: iSet
    type(silja_time) :: fu_start_time_fire_src
    integer :: iTmp
    if(present(iSet))then
      if(iSet < 1 .or. iSet > fire_src%nFRPdatasets)then
        fu_start_time_fire_src = time_missing
        call set_error('Incorrect set index:' + fu_str(iSet),'fu_start_time_fire_src')
      else
        fu_start_time_fire_src = fire_src%first_day(iSet)
      endif
    else
      fu_start_time_fire_src = fire_src%first_day(1)
      do iTmp = 2, fire_src%nFRPdatasets
        if(fu_start_time_fire_src > fire_src%first_day(iTmp)) &
                        & fu_start_time_fire_src = fire_src%first_day(iTmp)
        if (error) exit
      end do
    endif
  end function fu_start_time_fire_src
  
  !=========================================================================
  function fu_end_time_fire_src(fire_src, iSet)
    implicit none
    type(Tsilam_fire_source), intent(in) :: fire_src
    integer, intent(in), optional :: iSet
    type(silja_time) :: fu_end_time_fire_src
    integer :: iTmp
    if(present(iSet))then
      if(iSet < 1 .or. iSet > fire_src%nFRPdatasets)then
        fu_end_time_fire_src = time_missing
        call set_error('Incorrect set index:' + fu_str(iSet),'fu_end_time_fire_src')
      else
        fu_end_time_fire_src = fire_src%last_day(iSet)
      endif
    else
      fu_end_time_fire_src = fire_src%last_day(1)
      do iTmp = 2, fire_src%nFRPdatasets
        if(fu_end_time_fire_src < fire_src%last_day(iTmp)) &
                        & fu_end_time_fire_src = fire_src%last_day(iTmp)
      end do
    endif
  end function fu_end_time_fire_src
  
  !=========================================================================

  function fu_duration_fire_src(fire_src, iSet)
    implicit none
    type(Tsilam_fire_source), intent(in) :: fire_src
    integer, intent(in), optional :: iSet
    type(silja_interval) :: fu_duration_fire_src
    integer :: iTmp
    if(present(iSet))then
      if(iSet < 1 .or. iSet > fire_src%nFRPdatasets)then
        fu_duration_fire_src = interval_missing
        call set_error('Incorrect set index:' + fu_str(iSet),'fu_duration_fire_src')
      else
        fu_duration_fire_src = fire_src%last_day(iSet) + one_day - fire_src%first_day(iSet)
      endif
    else
      fu_duration_fire_src = fu_end_time_fire_src(fire_src) - fu_start_time_fire_src(fire_src)
    endif

  end function fu_duration_fire_src
  
  !========================================================================
  integer function fu_n_fires(fire_src, iSet)
    implicit none
    type(Tsilam_fire_source), intent(in) :: fire_src
    integer, intent(in), optional :: iSet
    type(silja_time) :: fu_end_time_fire_src
    integer :: iTmp
    if(present(iSet))then
      if(iSet < 1 .or. iSet > fire_src%nFRPdatasets)then
        fu_n_fires = int_missing
        call set_error('Incorrect set index:' + fu_str(iSet),'fu_n_fires')
      else
        fu_n_fires = fire_src%pFRPset(iSet)%nFires
      endif
    else
      fu_n_fires = 0
      do iTmp = 1, fire_src%nFRPdatasets
        fu_n_fires = fu_n_fires + fire_src%pFRPset(iTmp)%nFires
      end do
    endif
  end function fu_n_fires


  !*************************************************************************

  subroutine report_fire_src(fire_src)
    implicit none
    type(Tsilam_fire_source), intent(in) :: fire_src
    integer :: iTmp

    call msg('------- Fire source term')
    call msg('Fire source:' + fire_src%src_nm + '_' + fire_src%sector_nm)
    call msg('Number of fires:',fu_n_fires(fire_src))
    call msg('Active time:' + fu_str(fu_start_time_fire_src(fire_src)) + '---' + &
                            & fu_str(fu_end_time_fire_src(fire_src)))
    call report_metadata(fire_src%pFMD)
    call msg('======= End fire source term')

  end subroutine report_fire_src

  !**************************************************************************
  
  subroutine report_metadata(mdata)
    !
    ! Reports the fire source metadata
    !
    implicit none
    
    ! Imported parameter
    type(Tsilam_fire_metadata), pointer :: mdata

    ! Local variables
    character(len=worksize_string) :: sp
    integer :: iSp, iLU, iFlSm


    call msg('Fire source metadata:')
    call msg('Initialization file:' + mdata%chIniFNm)
    call msg('Emitted species and coefficients (kg_or_mole/J) for the land-uses:', &
           & mdata%nSpecies, mdata%nLU_types)

    do iFlSm = 1,2
      call msg("")
      write(unit=sp,fmt='(A3,A11,X,30(A10,1x))') "LUT",flamsmold(iFlSm), &
                                      & (fu_str(mdata%species_flaming(iSp)),iSp=1,mdata%nSpecies)
      call msg(sp)
      do iLU = 1, mdata%nLU_types
          write(unit=sp,fmt='(A14,X,30(E10.3,1x))') mdata%chLU_names(iLU), &
                                      & (mdata%pEmsFactors(flaming,iLU,iSp),iSp=1,mdata%nSpecies)
       
          call msg(sp)
      end do
    enddo


  end subroutine report_metadata


  !**************************************************************************
  
  logical function defined_fire_src(fire_src)
    implicit none
    type(Tsilam_fire_source), intent(in) :: fire_src

    defined_fire_src = fu_true(fire_src%defined)

  end function defined_fire_src

END MODULE source_term_fires


