MODULE source_terms_sea_salt
  !
  ! This module contains description of the sea salt emission.
  !
  ! Emission is computed for the prescribed size classes and needs wind speed and sea water
  ! temperature and salinity
  !
  ! All units: SI, unless otherwise stated
  !
  ! Code owner Mikhail Sofiev mikhail.sofiev@fmi.fi
  !
  ! Language: ANSI FORTRAN-90 (or close to it)
  !
  ! 2018 February: Added halogen emissions by R. Hanninen. Actually only Br2 for now.
  !                Only Br2 has a reasonalbe function for the depletion factor DF(d) as a function
  !                of sea-salt diameter. The species and their amounts should likely be tuned. 
  !   Notes:
  !     - srcSeaSalt%species contains also the halogens:
  !       srcSeaSalt%nSpecies = srcSeaSalt%nSpeciesSeaSalt + srcSeaSalt%nHalogens
  !     - Subroutine init_emission_sea_salt has an extra loop for halogens whose amount
  !       depend on all the sea-salt species. 
  !     - Largest changes are in subroutines sea_salt_flux4mode_hybrid and
  !       sea_salt_flux4mode_spume, where the volume flux is replaced with mass flux
  !       in order to properly calculate the halogen emissions. 
  use source_terms_time_params !cocktail_basic

  implicit none

  !
  ! PUBLIC routines of sea salt
  !
  public fill_sslt_src_from_namelist
  public reserve_sea_salt_source
  public init_emission_sea_salt
  public create_source_containing_grid
  public source_2_second_grid
  public add_inventory_sslt_src
  public add_input_needs
  public link_source_to_species
  public compute_emission_for_sea_salt
  public fu_sslt_emis_owned_quantity
  public fu_name
  public fu_source_nbr
  public fu_source_id_nbr
  public typical_species_conc
  public sea_salt_flux4mode
  public report

  !
  ! Private routines of the sea salt source
  !
  private add_input_needs_sea_salt_src
  private create_src_cont_grd_sslt_src
  private project_sslt_src_second_grd
  private link_sslt_src_to_species
  private sea_salt_flux4mode_hybrid ! general-size flux computation
  private fu_shape_function         ! Merge of Monahan and Martensson
  private fu_salinity_corr          ! flux correction for non-oceanic salinity
  private fu_water_tempr_corr_to_25deg   ! flux correction for non-25C temperature
  private fu_water_tempr_corr_to_20deg   ! ... and to non-20C
  private fu_water_tempr_corr_to_15deg   ! ... and to non-15C
  private sea_salt_flux4mode_spume
  private fu_shape_function_spume
  private fu_source_id_nbr_of_sslt_src
  private fu_source_nbr_of_sslt_src
  private fu_sea_salt_source_name
  private typical_species_cnc_sslt
  private report_sea_salt_src

  !
  ! Private subs of sea salt source
  !
  interface add_input_needs
    module procedure add_input_needs_sea_salt_src
  end interface

  interface create_source_containing_grid
    module procedure create_src_cont_grd_sslt_src
  end interface

  interface source_2_second_grid
    module procedure project_sslt_src_second_grd
  end interface

  interface link_source_to_species
    module procedure link_sslt_src_to_species
  end interface

  interface fu_source_nbr
    module procedure fu_source_nbr_of_sslt_src
  end interface

  interface fu_source_id_nbr
    module procedure fu_source_id_nbr_of_sslt_src
  end interface

  interface fu_name
    module procedure fu_sea_salt_source_name
  end interface

  interface typical_species_conc
    module procedure typical_species_cnc_sslt
  end interface

  interface sea_salt_flux4mode
    module procedure sea_salt_flux4mode_hybrid
  end interface

  interface report
    module procedure report_sea_salt_src
  end interface

  !
  ! There might be several types of the emission algorithm. Therefore, the below list 
  ! of parameters might eventually grow
  !
!  integer, private, parameter :: emisMonMar_w10m = 4101
!  integer, private, parameter :: emisMonahan_w10m = 4102
  integer, private, parameter :: emisHybrid_w10m = 4103
!  integer, private, parameter :: emisMonMar_fricVel = 4102
!  integer, private, parameter :: emisMonahan_fricVel = 4104
  integer, private, parameter :: emisHybrid_spume_w10m = 4105

  !
  ! Since the dependences of flux on diameter, water temperature and salinity are
  ! all tricky, all or most of integrals are to be taken numerically. To save 
  ! this pain during runtime, we pre-compute them at the very beginning and will
  ! later just pick the values from the array defined here.
  ! So far, we pre-define the ranges for temperature from -5 to +25 degrees C
  ! and salinity from 0 to 0.04. Typical values: 5 deg C northern seas and oceans
  ! 20-25 for southern Atlantic, 0.033 salinity for ocean, 0.009 salinity of
  ! the Baltic Sea.
  !
  integer, private, parameter :: nWaterTemps = 30
  integer, private, parameter :: nSalinities = 40
  real, private, parameter :: fMinWaterTemp = 268.15 ! [K]
  real, private, parameter :: fMinSalinity = 0. ! [fraction]
  real, private, parameter :: stepWaterTemp = 1.0 ! [K]
  real, private, parameter :: stepSalinity = 1.e-3 ! fraction

  !
  ! The Sea salt source term
  !
  TYPE silam_sea_salt_source
    PRIVATE
    CHARACTER(len=clen) :: src_nm, sector_nm  ! Name of the area source and sector
    integer :: emisMethod, src_nbr, id_nbr                ! A source and id numbers in a WHOLE source list
    real, dimension(:,:,:), pointer :: fluxPerModeNbr, fluxPerModeMass ! (nModes, nWaterTemps, nSalinities)
    real, dimension(:,:), pointer :: fluxPerModeSpumeNbr, fluxPerModeSpumeMass ! (nModes, nSalinities)
    integer :: nLevsDispVert, nSpecies, nSpeciesSeaSalt, nHalogens !nSpecies=nSpeciesSeaSalt+nHalogens
    logical :: ifWaterSalinityCounts, ifIceFractionCounts
    integer :: iWaterTempHandling
    real :: minOpenWaterFraction, defaultWaterTempr, defaultWaterSalinity
    integer :: uWindFlag, vWindFlag
    type(silam_vertical) :: vertLevsDispVert
    real, dimension(:), pointer :: levFractDispVert, fzDisp
    type(Tsilam_namelist), pointer :: nlInputFiles  ! namelist for names of supplementary files
    type(silam_species), dimension(:), pointer :: species
    type(silja_field), pointer :: source_mask
    logical, dimension(:), pointer :: ifNumberFlux
    type(chemical_adaptor) :: adaptor
    type(silja_logical) :: defined
  END TYPE silam_sea_salt_source

  type sslt_src_ptr
    type(silam_sea_salt_source) :: sslt_src
  end type sslt_src_ptr

                                                    ! default is 20C reference temperature
  logical, private, parameter :: if_25_Degrees_Reference = .false., &  ! if 25C is the reference
                               & if_15_Degrees_Reference = .false., &  ! if 15C is the reference
                               & if_low_spume = .true.


CONTAINS


  !*********************************************************************

  subroutine fill_sslt_src_from_namelist(nlSetup, srcSeaSalt, expected_species, chDataDir)
    !
    ! Initializes the sea salt source term.
    ! The parameters are read from the given ini file and stored to the returned 
    ! source term. After the stuff has been read, the main FluxPerMode is to be generated.
    ! This configuration allows for several sea salt sources
    !
    implicit none

    ! Imported parameters
    type(Tsilam_namelist), intent(in) :: nlSetup
    type(silam_sea_salt_source), intent(inout) :: srcSeaSalt
    type(silam_species), dimension(:), intent(in), allocatable :: expected_species
    character(len=*), intent(in) :: chDataDir

    ! Local variables
    integer :: iTmp, indWaterTemp, indSalinity, iSpecies, nFiles
    type(Tsilam_nl_item_ptr), dimension(:), pointer ::  pItems
    type(Taerosol) :: aerosolTmp
    type(silam_species), dimension(1) :: species
    integer, dimension(:), pointer :: indices
    logical :: ifFound

    !
    ! Names
    !
    srcSeaSalt%src_nm = fu_content(nlSetup,'source_name')
    srcSeaSalt%sector_nm = fu_content(nlSetup,'source_sector_name')
    srcSeaSalt%defined = silja_false

    !
    ! Emission index type
    !
    select case(fu_str_u_case(fu_content(nlSetup,'sea_salt_emission_method')))

!      case ('MONAHAN_MARTENSSON_WIND_10M')
!        srcSeaSalt%emisMethod = emisMonMar_w10m
!
!      case ('MONAHAN_MARTENSSON_FRICTION_VEL')
!        rulesSeaSalt%emisMethod = emisMonMar_fricVel
!
!      case ('MONAHAN_WIND_10M')
!        srcSeaSalt%emisMethod = emisMonahan_w10m

      case ('HYBRID_WIND_10M')
        srcSeaSalt%emisMethod = emisHybrid_w10m

      case ('HYBRID_AND_SPUME_WIND_10M')
        srcSeaSalt%emisMethod = emisHybrid_spume_w10m

!      case ('MONAHAN_FRICTION_VEL')
!        rulesSeaSalt%emisMethod = emisMonahan_fricVel

      case default
        call set_error('Unknown emission method:' + fu_content(nlSetup,'sea_salt_emission_method'), &
                     & 'fill_sslt_src_from_namelist')
        return
    end select

    ! Whether water temperature and salinity are taken into account. Note that
    ! in default case we assume crude emission computation
    !
    if(fu_str_u_case(fu_content(nlSetup,'water_temperature_input_type')) == 'FIXED_VALUE')then
      srcSeaSalt%iWaterTempHandling = static_value

    elseif(fu_str_u_case(fu_content(nlSetup,'water_temperature_input_type')) == 'FIXED_MAP')then
      srcSeaSalt%iWaterTempHandling = static_climatology

    elseif(fu_str_u_case(fu_content(nlSetup,'water_temperature_input_type')) == 'MONTHLY_CLIMATOLOGY')then
      call set_error('Monthly climatology is not yet available for water temparture dependence', &
                   & 'fill_sslt_src_from_namelist')
      return

    elseif(fu_str_u_case(fu_content(nlSetup,'water_temperature_input_type')) == 'DYNAMIC')then
      srcSeaSalt%iWaterTempHandling = dynamic_map

    else
      call set_error('Unknown type of water_temperature_input_type:' + &
                   & fu_content(nlSetup,'water_temperature_input_type'), &
                   & 'fill_sslt_src_from_namelist')
      return
    endif

    srcSeaSalt%ifWaterSalinityCounts = fu_str_u_case( &
                     & fu_content(nlSetup,'sea_salt_emis_depend_on_water_salinity')) == 'YES'
    srcSeaSalt%ifIceFractionCounts = fu_str_u_case( &
                     & fu_content(nlSetup,'sea_salt_emis_depend_on_ice_fraction')) == 'YES'
    !
    ! Some stuff can be taken by-default: salinity, temperature, minimum open-water fraction, 
    !
    srcSeaSalt%defaultWaterSalinity = fu_content_real(nlSetup,'default_water_salinity')
    if(error .or. srcSeaSalt%defaultWaterSalinity < 0. .or. srcSeaSalt%defaultWaterSalinity > 1.)then
      call set_error('No default_water_salinity ??','fill_sslt_src_from_namelist')
      if(error)call unset_error('fill_sslt_src_from_namelist')
      srcSeaSalt%defaultWaterSalinity = 0.033
    endif
    if(srcSeaSalt%iWaterTempHandling == static_value)then
      srcSeaSalt%defaultWaterTempr = fu_content_real(nlSetup,'default_water_temperature')
      if(error .or. srcSeaSalt%defaultWaterTempr < 200. .or. srcSeaSalt%defaultWaterTempr > 400.)then
        call set_error('No default_water_temperature ??','fill_sslt_src_from_namelist')
        call unset_error('fill_sslt_src_from_namelist')
        srcSeaSalt%defaultWaterTempr = 285.  ! Something reasonable for the Atlantic ocean
      endif
    endif
    srcSeaSalt%minOpenWaterFraction = fu_content_real(nlSetup,'min_open_water_area_fraction')
    if(error .or. srcSeaSalt%minOpenWaterFraction < 0. .or. srcSeaSalt%minOpenWaterFraction > 1.)then
      call msg_warning('No min_open_water_area_fraction ??','fill_sslt_src_from_namelist')
      if(error)call unset_error('fill_sslt_src_from_namelist')
      srcSeaSalt%minOpenWaterFraction = 0.  ! No limitations unless stated explicitly
    endif

    !
    ! What do we use as a near-surface wind: wind_10m or the wind at lowest model level?
    !
    if(fu_str_u_case(fu_content(nlSetup,'wind_selection')) == 'WIND_10M')then
      srcSeaSalt%uWindFlag = u_10m_flag
      srcSeaSalt%vWindFlag = v_10m_flag
    elseif(fu_str_u_case(fu_content(nlSetup,'wind_selection')) == 'WIND_LEVEL_1')then
      srcSeaSalt%uWindFlag = u_flag
      srcSeaSalt%vWindFlag = v_flag
    else
      call set_error('wind_selection allows only: WIND_10M, WIND_LEVEL_1, not:' + &
                  & fu_content(nlSetup,'wind_selection'),'fill_sslt_src_from_namelist')
      return
    endif
    !
    ! Now the list of aerosol modes that will be emitted. To set them up, we will
    ! use the standard aerosol procedure - for the sake of unification.
    ! Note that there are two potentially concurring definitions: one coming from aerosol dynamics
    ! the other - written in the ini file. The first one prevails, if it exists
    !
    call get_source_aer_species(srcSeaSalt%species, srcSeaSalt%nSpecies, &
                              & expected_species, fu_content(nlSetup,'sea_salt_substance_name'), &
                              & nlSetup)
    if(error)return
    !
    ! The sea salt source can emit halogens if we ask for it.
    !
    srcSeaSalt%nSpeciesSeaSalt = srcSeaSalt%nSpecies !Save the number of sea-salt species
    if (fu_content(nlSetup,'if_halogen_emission') == '') then
       call msg('No if_halogen_emission option found in sea-salt ini. Assuming no halogen emissions!') 
       srcSeaSalt%nHalogens = 0
    elseif (fu_str_u_case(fu_content(nlSetup,'if_halogen_emission')) == 'YES') then
      call set_species(species(1), fu_get_material_ptr('Br2'), in_gas_phase)
      call addSpecies(srcSeaSalt%species, srcSeaSalt%nSpecies, species, 1)
      srcSeaSalt%nHalogens = 1
      !Omit the chlorine emissions from sea-salt for now. Use GEIA inventory instead for chlorines. 
      !call set_species(species(1), fu_get_material_ptr('HCl'), in_gas_phase) !HCl is better than Cl2
      !call addSpecies(srcSeaSalt%species, srcSeaSalt%nSpecies, species, 1)
      !srcSeaSalt%nHalogens = 2
    else
      srcSeaSalt%nHalogens = 0
    endif
    !
    ! If number-emission is needed, it has to be requested as nbr_aer in the sea salt species list
    !
    allocate(srcSeaSalt%ifNumberFlux(srcSeaSalt%nSpecies), stat=iTmp)
    if(fu_fails(iTmp == 0, 'Failed to allocate number-mass switcher array','fill_sslt_src_from_namelist'))return
    do iTmp = 1, srcSeaSalt%nSpecies
      srcSeaSalt%ifNumberFlux(iTmp) = (fu_name(fu_material(srcSeaSalt%species(iTmp))) == 'nbr_aer')
    end do

    !
    ! Store the input fields for the source features.
    !
    srcSeaSalt%nlInputFiles => fu_create_namelist('sea_salt_src_supplementary_files')
    if(error)return

    nullify(pItems)
    call get_items(nlSetup, 'supplementary_file', pItems, nFiles)
    if(error)return
!    if(nFiles < 1)then
!      call set_error('No supplementary files are given','fill_sslt_src_from_namelist')
!      return
!    endif
    do iTmp = 1, nFiles
      call add_namelist_item(srcSeaSalt%nlInputFiles, 'supplementary_file', &
!                           & fu_extend_grads_hat_dir(fu_content(pItems(iTmp)), chDataDir))
                           & fu_process_filepath(fu_content(pItems(iTmp)), &
                                               & convert_slashes = .true., must_exist = .false., &
                                               & superfile = chDataDir))
    end do

    !
    ! Note that the source mask cannot be read from the file yet: dispersion grid is undefined
    !
    call add_namelist_item(srcSeaSalt%nlInputFiles,  'source_area_mask', &
!                         & fu_extend_grads_hat_dir(fu_content(nlSetup,'source_area_mask'), chDataDir))
                         & fu_process_filepath(fu_content(nlSetup,'source_area_mask'), &
                                             & convert_slashes = .true., must_exist = .false., &
                                             & superfile = chDataDir))
    if(error)return

    srcSeaSalt%defined = silja_undefined

    call report(srcSeaSalt)
    
    
  end subroutine fill_sslt_src_from_namelist


  !*****************************************************************

  subroutine add_input_needs_sea_salt_src(sslt_src, q_met_dynamic, q_met_static, &
                                                  & q_disp_dynamic, q_disp_static)
    !
    ! Returns input needs for the emission/transformation routines. 
    !
    implicit none

    ! Imported parameters
    type(silam_sea_salt_source), intent(in) :: sslt_src
    integer, dimension(:), intent(inout) :: q_met_dynamic, q_met_static, &
                                          & q_disp_dynamic, q_disp_static

    ! Local variables
    integer :: iTmp

    !
    ! Add needed dynamic quantities. Emission needs wind at 10m and, in the future,
    ! friction velocity. We would also need water temperature and ice fraction
    !
    if(sslt_src%ifIceFractionCounts) &
                          & iTmp = fu_merge_integer_to_array(fraction_of_ice_flag, q_met_dynamic)
    !
    ! Depending on the type of the temperature field handling, have to request it to different stacks
    !
    if(sslt_src%iWaterTempHandling == static_climatology)then                   ! fixed map
      iTmp = fu_merge_integer_to_array(water_surface_temp_flag, q_met_static)

    elseif(sslt_src%iWaterTempHandling == dynamic_map)then                          ! dynamic map
      iTmp = fu_merge_integer_to_array(water_surface_temp_flag, q_met_dynamic)

    endif

    ! If water salinity should be considered, one has to include the map
    ! of salinity. 
    !
    if(sslt_src%ifWaterSalinityCounts)  &
                          & iTmp = fu_merge_integer_to_array(water_salinity_flag, q_met_static)

    select case(sslt_src%emisMethod)
      case(emisHybrid_w10m, emisHybrid_spume_w10m)
        iTmp = fu_merge_integer_to_array(sslt_src%uWindFlag, q_met_dynamic)
        iTmp = fu_merge_integer_to_array(sslt_src%vWindFlag, q_met_dynamic)
        if(sslt_src%uWindFlag == u_flag) iTmp = fu_merge_integer_to_array(height_flag, q_met_dynamic)
!!!        iTmp = fu_merge_integer_to_array(windspeed_10m_flag, q_met_dynamic)
        iTmp = fu_merge_integer_to_array(surface_roughness_meteo_flag, q_met_dynamic)
        iTmp = fu_merge_integer_to_array(friction_velocity_flag, q_met_dynamic)

      case default
        call set_error('Unknown emission computation algorythm', &
                     & 'input_needs_sea_salt_source')
        return
    end select

    ! Finally, a land fraction field must always be in the permanent input
    !
    iTmp = fu_merge_integer_to_array(fraction_of_land_flag, q_met_static)


  end subroutine add_input_needs_sea_salt_src


  !**************************************************************************

  subroutine reserve_sea_salt_source(sslt_src, &     ! Src to initialise
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
    type(silam_sea_salt_source), intent(inout) :: sslt_src
    integer, intent(in) :: iSrcNbr, iSrcIdNbr

    !
    ! Nullify the basic variables
    !
    sslt_src%src_nm = ''
    sslt_src%sector_nm = ''

    !
    ! Main source parameters - enough to identify it in the global information list
    !
    sslt_src%src_nbr = iSrcNbr
    sslt_src%id_nbr = iSrcIdNbr
    !
    ! A bit of other stuff
    !
    nullify(sslt_src%fZDisp)
    nullify(sslt_src%levFractDispVert)
    !
    ! Finally, mark the source as incomplete
    !
    sslt_src%defined = silja_false

  end subroutine reserve_sea_salt_source


  !*********************************************************************

  subroutine init_emission_sea_salt(srcSeaSalt)
    !
    ! Initializes the sea salt source term.
    ! The parameters are read from the given ini file and stored to the returned 
    ! source term. After the stuff has been read, the main FluxPerMode is to be generated.
    ! This configuration allows for several sea salt sources
    !
    implicit none

    ! Imported parameters
    type(Tsilam_namelist), pointer :: nlSetup
    type(silam_sea_salt_source), intent(inout) :: srcSeaSalt

    ! Local variables
    integer :: iTmp, indWaterTemp, indSalinity, iSpecies, iHalog
    type(silam_sp) :: sp
    real :: fNbrFlux, fMassFlux, fMassMeanDiam
    type(silja_field_id) :: id, idOut
    real, dimension(:), pointer :: arPtr
    real, dimension(:), pointer :: fPartDensity

    !
    ! First of all, read the source mask
    !
    sp%sp => fu_work_string()
    if(error)return
    arPtr => fu_work_array()
    if(error)return

    sp%sp = fu_process_filepath(fu_content(srcSeaSalt%nlInputFiles, 'source_area_mask'))
    if(error .or. len_trim(sp%sp) < 1)then
      call set_error('Source area mask is absent','init_emission_sea_salt')
      return
    endif
    id = fu_set_field_id_simple(met_src_missing, fraction_of_water_flag, time_missing, level_missing)
    if(error)return
    call set_grid(id, dispersion_grid)
    if(error)return

!    id = fu_set_field_id(met_src_missing, &       ! met_src
!                       & fraction_of_water_flag, & ! quantity
!                       & time_missing, &          ! analysis time
!                       & interval_missing, &      ! forecast length
!                       & dispersion_grid, &       ! grid
!                       & level_missing, &         ! level
!                       & interval_missing, &      ! length of accumulation
!                       & interval_missing, &      ! length of validity
!                       & analysis_flag) !, &         ! type of the field
!!                       & 'sslt', &                ! substance name
!!                       & aerosol_mode_missing)    ! aerosol mode
    if(error)return

    call get_input_field(sp%sp(index(adjustl(sp%sp),' ')+1:), &  ! file name
                       & fu_input_file_format(sp%sp), &          ! file format
                       & id, &                  ! The id to search
                       & arPtr, &               ! data array
                       & dispersion_gridPtr, &  ! storage grid
                       & iOutside = nearestPoint, &         ! out of grid interpolation
                       & iAccuracy = 5, &
                       & wdr = wdr_missing, & 
!                       & iBinary, &           ! for NetCDF
                       & ifAcceptSameMonth = .false., &
                       & idOut = idOut)  ! redefine the id
    if(error)return

    if(defined(idOut)) then
      allocate(srcSeaSalt%source_mask, stat=iTmp)
      if(iTmp /= 0)then
        call set_error('Failed to allocate sea salt source mask','init_emission_sea_salt')
        return
      endif
      call set_field(idOut, arPtr, srcSeaSalt%source_mask, .true.)
    else
      call set_error('Failed to get the source mask','init_emission_sea_salt')
    endif
    if(error)return

    !
    ! Allocate proper size of the other arrays
    !
    allocate(srcSeaSalt%fluxPerModeNbr(srcSeaSalt%nSpecies, nWaterTemps, nSalinities), &
           & srcSeaSalt%fluxPerModeMass(srcSeaSalt%nSpecies, nWaterTemps, nSalinities), &
           & srcSeaSalt%fluxPerModeSpumeNbr(srcSeaSalt%nSpecies, nSalinities), &
           & srcSeaSalt%fluxPerModeSpumeMass(srcSeaSalt%nSpecies, nSalinities), stat=iTmp)
    if(iTmp /= 0)then
      call msg('Allocation status:',iTmp)
      call set_error('Failed to allocate sea salt precomputed flux array', &
                   & 'init_emission_sea_salt')
      return
    endif
   
    !
    ! Particle masses for each size class: !!!! dry particle !!!!
    !
    fPartDensity => fu_work_array()
    if(error)return

    do iSpecies = 1, srcSeaSalt%nSpeciesSeaSalt !NOTE: Omit the halogens
       fPartDensity(iSpecies) = fu_dry_part_density(srcSeaSalt%species(iSpecies)%material) !RISTO NOTE: volume -> mass
    end do
    if(error)return

    !
    ! Now let's precompute the array of fluxes per size mode
    ! For each mode we have to take care of all temperatures and salinities
    !
    select case(srcSeaSalt%emisMethod)
      case(emisHybrid_w10m, emisHybrid_spume_w10m)
        !Initialize these fluxes to zero, since for halogens one needs to sum some of these using all the sea-salt species. 
        srcSeaSalt%fluxPerModeNbr(1:srcSeaSalt%nSpecies, 1:nWaterTemps, 1:nSalinities) = 0.0
        srcSeaSalt%fluxPerModeMass(1:srcSeaSalt%nSpecies, 1:nWaterTemps, 1:nSalinities) = 0.0
        srcSeaSalt%fluxPerModeSpumeNbr(1:srcSeaSalt%nSpecies, 1:nSalinities) = 0.0       
        srcSeaSalt%fluxPerModeSpumeMass(1:srcSeaSalt%nSpecies, 1:nSalinities) = 0.0
        !do iSpecies = 1, srcSeaSalt%nSpecies
        do iSpecies = 1, srcSeaSalt%nSpeciesSeaSalt !When halogens are present this loop must exclude the halogens
          do iHalog = 0, srcSeaSalt%nHalogens !Index 0 is for standard sea-salt, and other indeces correspond to halogens.  
            do indSalinity = 1, nSalinities
              do indWaterTemp = 1, nWaterTemps
                !
                ! Bubble-mediated mechanism
                !
                call sea_salt_flux4mode_hybrid(srcSeaSalt%species(iSpecies)%mode, &
                                           & fMinWaterTemp + (indWaterTemp-1) * stepWaterTemp, &
                                           & fMinSalinity + (indSalinity-1) * stepSalinity, &
                                           & srcSeaSalt%fluxPerModeNbr(iSpecies+iHalog,indWaterTemp,indSalinity), &
                                           & srcSeaSalt%fluxPerModeMass(iSpecies+iHalog,indWaterTemp,indSalinity), &
                                           & fMassMeanDiam, &
                                           & fPartDensity(iSpecies), &
                                           & iHalog)

!                write(unit = sp%sp,fmt=*)'T_water, S_water, flux(all modes):', &
!                                    & fMinWaterTemp + (indWaterTemp-1) * stepWaterTemp, & ! [K]
!                                    & fMinSalinity + (indSalinity-1) * stepSalinity, &
!                                    & (fluxArray(indWaterTemp,indSalinity)%fFluxPerMode(i), i=1,fu_n_modes(cocktail%aerosol))
!                call msg(sp%sp)

               end do  ! nWater temps
               !
               ! If spume mechanism accounted for, add it
               !
               if(srcSeaSalt%emisMethod == emisHybrid_spume_w10m)then
                 call sea_salt_flux4mode_spume(srcSeaSalt%species(iSpecies)%mode, &
                                          & fMinSalinity + (indSalinity-1) * stepSalinity, &
                                          & srcSeaSalt%fluxPerModeSpumeNbr(iSpecies+iHalog,indSalinity), &
                                          & srcSeaSalt%fluxPerModeSpumeMass(iSpecies+iHalog,indSalinity), &
                                          & fMassMeanDiam, &
                                          & fPartDensity(iSpecies), &
                                          & iHalog)
               endif
             end do  ! salinities
          end do   ! halogens
        end do    ! sea-salt species

      case default
        call msg('Unknown emission method:',srcSeaSalt%emisMethod)
        call set_error('Unknown emission method','init_emission_sea_salt')
        return
    end select
    !
    ! Note that now we are setting the mean diameter of aerosol modes. After the above cycle they are set to 
    ! some extreme values of temperature and salinity. Now we call it once again with the "typical" value
    !
    if(if_25_Degrees_Reference)then
      call msg('========> Reference temperature is 25 degrees')
    else
      call msg('========> Reference temperature is 15 degrees')
    endif
    call msg('')
    select case(srcSeaSalt%emisMethod)
      case(emisHybrid_w10m, emisHybrid_spume_w10m)
        call msg('Bubble-mediated mechanism, 15C, 0.03 salinity')
        !do iSpecies = 1, srcSeaSalt%nSpecies
        do iSpecies = 1, srcSeaSalt%nSpeciesSeaSalt !NOTE: Omit the halogens!
          call sea_salt_flux4mode_hybrid(srcSeaSalt%species(iSpecies)%mode, &
                                       & 288.0, &   ! 15 C
                                       & 0.03, &    ! ocean
                                       & fNbrFlux, &   ! number flux
                                       & fMassFlux, &   ! volume flux !RISTO NOTE: volume -> mass
                                       & fMassMeanDiam, &
                                       & fPartDensity(iSpecies))
          write(unit = sp%sp,fmt='(A48,1x,A,F6.1,2(1x,F5.2),A,3(1x,E15.7))') &
                                 & 'T_water,S_water,#-flux, m-flux, m-flux-fake for', &
                                 & trim(fu_str(srcSeaSalt%species(iSpecies))), &
                                 & 288.0, 0.03, &
                                 & fMassMeanDiam*1.e6, ' mkm', &
!                                 & fNbrFlux, fMassFlux, fPartDensity(iSpecies)*fNbrFlux*(Pi/6.)*fMassMeanDiam**3 !RISTO NOTE: volume -> mass
                                 & fNbrFlux, fMassFlux/fPartDensity(iSpecies), fNbrFlux*(Pi/6.)*fMassMeanDiam**3  !NOTE: Output as volume flux (as before)
          call msg(sp%sp)
          !
          ! Store the results for this size mode. Note that we ignore spume mode altogether.
          ! Since wind-speed dependence is different for these mechanisms, the mean diameter 
          ! will be different and thus dynamic. But in most cases at regional scale spume 
          ! particles will not be too important.
          !
          call set_massmean_d(srcSeaSalt%species(iSpecies)%mode, fMassMeanDiam, &
                         & (srcSeaSalt%ifWaterSalinityCounts .or. (srcSeaSalt%iWaterTempHandling /= static_value)))
          if(error)return
        end do
        !
        ! Report spume mechanism
        !
        IF(srcSeaSalt%emisMethod == emisHybrid_spume_w10m)then
          call msg('')
          call msg('Spume-mediated mechanism, 0.03 salinity')
          !do iSpecies = 1, srcSeaSalt%nSpecies
          do iSpecies = 1, srcSeaSalt%nSpeciesSeaSalt !NOTE: Omit the halogens!
            call sea_salt_flux4mode_spume(srcSeaSalt%species(iSpecies)%mode, &
                                        & 0.03, &
                                        & fNbrFlux, &
                                        & fMassFlux, & !RISTO NOTE: volume -> mass
                                        & fMassMeanDiam, &
                                        & fPartDensity(iSpecies))
            write(unit = sp%sp,fmt='(A40,1x,A,2(1x,F5.2),A,3(1x,E15.7))') &
                                   & 'S_water,#-flux, m-flux, m-flux-fake for', &
                                   & trim(fu_str(srcSeaSalt%species(iSpecies))), &
                                   & 0.03, &
                                   & fMassMeanDiam*1.e6, ' mkm', &
!                                   & fNbrFlux, fMassFlux, fPartDensity(iSpecies)*fNbrFlux*(Pi/6.)*fMassMeanDiam**3 !RISTO NOTE: volume -> mass
                                   & fNbrFlux, fMassFlux/fPartDensity(iSpecies), fNbrFlux*(Pi/6.)*fMassMeanDiam**3  !NOTE: Output as volume flux (as before)
            call msg(sp%sp)
          end do    ! species
          call msg('')
        endif

      case default
        call msg('Unknown emission method:',srcSeaSalt%emisMethod)
        call set_error('Unknown emission method','init_emission_sea_salt')
        return
    end select

    call free_work_array(fPartDensity) !RISTO NOTE: volume -> mass
    
    srcSeaSalt%defined = silja_true

    call free_work_array(sp%sp)

  end subroutine init_emission_sea_salt


  !****************************************************************************

  subroutine add_inventory_sslt_src(sslt_src, species_list, nSpecies)
    !
    ! Get the species emitted by this source and make a list of
    ! them. The source must be initialized, of course.
    !
    implicit none
    type(silam_sea_salt_source), intent(in) :: sslt_src
    type(silam_species), dimension(:), pointer :: species_list
    integer, intent(inout) :: nSpecies

    call addSpecies(species_list, nSpecies, sslt_src%species, sslt_src%nSpecies) !RISTO: CHECK Likely nSpecies because lists all emitted species!

  end subroutine add_inventory_sslt_src


  !*******************************************************************************
  
  subroutine create_src_cont_grd_sslt_src(sslt_src, grid_template, ifVerbose, ifMinimal, ifExtended)
    !
    ! Creates the grid that covers the area with active BVOC emission
    !
    implicit none

    ! Imported parameters
    type(silam_sea_salt_source), intent(in) :: sslt_src
    type(silja_grid), intent(inout) :: grid_template
    logical, intent(in) :: ifVerbose, ifMinimal
    logical, intent(out) :: ifExtended

    ! So far nothing to do: this source just covers the dispersion grid
    !
    ifExtended = .false.
    return
    
  end subroutine create_src_cont_grd_sslt_src


  !*****************************************************************

  subroutine project_sslt_src_second_grd(sslt_src, grid, vert_disp, vert_proj, iAccuracy)
    !
    ! There is nothing to project in terms of horizontal grid but the vertical has to be made.
    !
    implicit none

    ! Imported parameters
    type(silam_sea_salt_source), intent(inout) :: sslt_src
    type(silja_grid), intent(in) :: grid
    type(silam_vertical), intent(in) :: vert_disp, vert_proj
    integer, intent(in) :: iAccuracy

    ! Local variables
    real, dimension(:), pointer :: arTmp
    integer :: i
    type(silam_vertical) :: vertTmp

    if(iAccuracy < 1 .or. iAccuracy > 10)then
      call msg('Accuracy switch must be from 1 to 10, not:',iAccuracy)
      call msg_warning('Accuracy switch must be from 1 to 10','project_sslt_src_second_grd')
!      return
    endif

    if(len_trim(sslt_src%sector_nm) > 0)then
      call msg('Re-projecting sea salt source:' + sslt_src%src_nm +'_' + sslt_src%sector_nm)
    else
      call msg('Re-projecting sea salt source:' + sslt_src%src_nm)
    endif

    arTmp => fu_work_array()
    if(error)return

    !
    ! Now re-project the vertical grid of the source to the given vertical
    ! Since vertical distributions are very poorly known and crude, we do not need
    ! any precise meteo-dependent projection, crude will do the job
    !
    sslt_src%vertLevsDispVert = vert_disp
    allocate(sslt_src%levFractDispVert(fu_NbrOfLevels(vert_disp)), &
           & sslt_src%fzDisp(fu_NbrOfLevels(vert_disp)), stat=i)
    if(fu_fails(i == 0, 'Failed allocation dispersion-vertical fractions', &
                                                         & 'project_sslt_src_second_grd'))return
    sslt_src%levFractDispVert(:) = 0.0
    sslt_src%fzDisp(:) = 0.0

    !
    ! Create the sea salt vertical, which can depend on the emission method
    !
    select case(sslt_src%emisMethod)
      
      case(emisHybrid_w10m, emisHybrid_spume_w10m)
        !
        ! Monahan-Martensson-Clarke hybrid presently goes from 5 to 50 m.
        !
        call set_vertical(fu_set_layer_between_two(layer_btw_2_height, 0.0, 25.0), vertTmp)
        if(error)return
        call add_level(vertTmp, fu_set_layer_between_two(layer_btw_2_height, 25.0, 50.0))
        if(error)return
        arTmp(1) = 0.6
        arTmp(2) = 0.4
      case default
        call set_error('Unknown emission method:'+fu_str(sslt_src%emisMethod),'project_sslt_src_second_grd')
        return
    end select

    call reproject_verticals(vertTmp, arTmp, &                    ! vertical from, fractions from
                           & vert_proj, sslt_src%levFractDispVert, &   ! vertical to, fractions to
                           & sslt_src%fzDisp, sslt_src%nLevsDispVert, & ! mass centres, number of non-zero levels
                           & ifMassCentreInRelUnit=.true.)
    call free_work_array(arTmp)
    call set_missing(vertTmp, .false.)

  end subroutine project_sslt_src_second_grd


  !**************************************************************************

  subroutine link_sslt_src_to_species(species_list, sslt_src)
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
    type(silam_sea_salt_source), intent(inout) :: sslt_src
    type(silam_species), dimension(:), pointer :: species_list

    !
    ! Linkage is actually just creation of the chemical adaptor
    !
    call create_adaptor(sslt_src%species, species_list, sslt_src%adaptor)
    
  end subroutine link_sslt_src_to_species


  !**************************************************************************

  subroutine compute_emission_for_sea_salt(sslt_src, &
                                         & met_buf, disp_buf, & 
                                         & now, &      ! current time
                                         & timestep, & ! model time step
                                         & pHorizInterpMet2DispStruct, ifHorizInterp, &
                                         & ifSpeciesMoment, &
                                         & emisMap, mapCoordX, mapCoordY, mapCoordZ, & ! Output
                                         & fMassInjected)                              ! output
    !
    ! Computes the emission fields for sea salt. Since the sea salt is a wide-specturm
    ! aerosol, we have to compute the emission flux for each size class separately.
    ! The computation is split to two parts: 
    ! 1. Computation of the sea area fraction covered with white caps - for each grid cell
    !    but common for all size classes. Here we take ice and land fraction into account.
    ! 2. Computation of the spectrum of the released particles, which depends on
    !    driving parameters selected for the run. It is based on a pre-computed array
    !    that covers all reasonable water temperatures and salinities, and provides
    !    the fluxes for all aerosol ranges of the given run.
    !
    ! This routine is to be called at each model time step but not inside the particle cycle,
    ! therefore its efficiency is of moderate importance.
    !
    ! ALL COMPUTATIONS proceed in METEO GRID
    !
    implicit none

    ! Imported parameters
    type(silam_sea_salt_source), target, intent(in) :: sslt_src
    type(Tfield_buffer), pointer ::  met_buf, disp_buf  ! meteo and internal field buffers
    type(silja_time), intent(in) :: now           ! current time
    type(silja_interval), intent(in) :: timestep  ! model time step
    type(THorizInterpStruct), intent(in) :: pHorizInterpMet2DispStruct
    logical, intent(in) :: ifHorizInterp, ifSpeciesMoment
    type(Tmass_map), intent(inout) :: emisMap, mapCoordX, mapCoordY, mapCoordZ
    real(r8k), dimension(:),  intent(inout) :: fMassInjected

    ! Local variables
    integer :: indIceFr, indSalinityMet, indWaterTemprMet, indFricVel, indZ0, iLev, & !indW10m, &
             & iMeteo, iDisp, ix, iy, ixMeteo, iyMeteo, iStat, iMode, iWaterTemp, iSalinity, &
             & iSeaSalt, indU, indV, indHeight, iHalogen
    real :: fTmp, frOpenWater, timestep_sec, fWhiteCaps, fCellTotal, u_star, u, v, z0, &
          & fSpumeWindFactor, windspeed, height, dtdxdyLevFrac
    real, dimension(:), pointer :: ptrXSzDisp, ptrYSzDisp, fMinDArray, fMaxDArray, pSrcMask
    !real, dimension(:), pointer :: ptrXSzDisp, ptrYSzDisp, fMinDArray, fMaxDArray, fPartDensity, pSrcMask
    logical, save :: ifFirstTime = .true.
    type(silja_field), pointer :: fldMaskPtr
    type(silam_sp) :: sp
    integer, save :: iCount = 0

    ! Local parameters 
    !
    ! For the water area covered with white caps
    real, parameter :: W_scale = 5e-6 !!!1e-4 !!3.84e-6
    real, parameter :: Wind10m_power = 3 !!!2 !! 3.41
!    real, parameter :: fCorrectionFactor = 1.5  ! TEMPORARY, to cope with too strong dry deposition
    !
    ! First, set the output pointer to the locally stored cocktail_map of emission
    ! intensity and get the temporary array for the white caps fraction.
    !
    pSrcMask => fu_grid_data(sslt_src%source_mask)

!open (50, file='sslt_src.dump',recl=nx_dispersion*ny_dispersion,form='unformatted',access='direct')
!write(50,rec=1)(pSrcMask(ix),ix=1,nx_dispersion*ny_dispersion)
!close (50)
!call msg('nx,ny',nx_dispersion, ny_dispersion)
!stop

    if(sslt_src%ifIceFractionCounts)then  ! Ice fraction counts. Find its index
      indIceFr = fu_index(met_buf, fraction_of_ice_flag) 
      if(error .or. indIceFr < 1)then
        call set_error('Failed to find ice fraction field','compute_emission_for_sea_salt')
        return
      endif
    endif

    ptrXSzDisp => fu_grid_data(dispersion_cell_x_size_fld)  ! get grid cell size in metres
    ptrYSzDisp => fu_grid_data(dispersion_cell_y_size_fld)

    timestep_sec = abs(fu_sec(timestep))
    if(error)return

    ! if water salinity counts, find the appropriate indices
    !
    if(sslt_src%ifWaterSalinityCounts)then
      indSalinityMet = fu_index(met_buf, water_salinity_flag) 
      if(error .or. indSalinityMet < 1)then
        call set_error('Failed to find water salinity field','compute_emission_for_sea_salt')
        return
      endif
      !
      ! Salinity map can be damaged with missing values, etc. Clean it if first time
      !
      if(ifFirstTime)then
        do iMeteo = 1, fs_meteo
          if(met_buf%p2d(indSalinityMet)%present%ptr(iMeteo) < 0. .or. &
           & met_buf%p2d(indSalinityMet)%present%ptr(iMeteo) > 1.)then
            met_buf%p2d(indSalinityMet)%present%ptr(iMeteo) = sslt_src%defaultWaterSalinity
          endif
        end do
      endif
    else
      iSalinity = int((sslt_src%defaultWaterSalinity-fMinSalinity) / stepSalinity + 0.5) + 1
      iSalinity = max(1,min(nSalinities,iSalinity))
    endif
    !
    ! If water temperature is not static, find appropriate indices
    !
    if(sslt_src%iWaterTempHandling == static_value)then
      iWaterTemp = int((sslt_src%defaultWaterTempr - fMinWaterTemp) / stepWaterTemp + 0.5) + 1
      iWaterTemp = max(1,min(nWaterTemps,iWaterTemp))
    else
      indWaterTemprMet = fu_index(met_buf, water_surface_temp_flag) 
      if(error .or. indWaterTemprMet < 1)then
        call set_error('Failed to find water temperature field','compute_emission_for_sea_salt')
        return
      endif
    endif

!!$    ! For Halogens the dry part density is now calculated in subr init_emission_sea_salt and where
!!$    ! the volume flux is replaced with massflux
!!$    !
!!$    ! Particle masses for each size class: !!!! dry particle !!!!
!!$    !
!!$    fPartDensity => fu_work_array()
!!$    if(error)return
!!$
!!$    !do iSeaSalt = 1, sslt_src%nSpecies
!!$    do iSeaSalt = 1, sslt_src%nSpeciesSeaSalt !changed due to Halogens
!!$       fPartDensity(iSeaSalt) = fu_dry_part_density(sslt_src%species(iSeaSalt)%material) !RISTO NOTE: volume -> mass
!!$    end do
!!$    if(error)return

    !--------------------------------------------------------------------
    !
    ! Computation of the white caps area depends on method of emission computation
    !
    select case(sslt_src%emisMethod)
      
      case(emisHybrid_w10m, emisHybrid_spume_w10m)
        !
        ! Basic case: a combination of Monahan and Martensson methods with their own
        ! wind-dependent parameterization of white caps on 10m wind speed
        ! Possibly, spume production
        !
!        indW10m = fu_index(met_buf, windspeed_10m_flag) 
!        if(error .or. indW10m < 1)then
!          call set_error('Failed to find 10m windspeed field','compute_emission_for_sea_salt')
!          return
!        endif

        indU = fu_index(met_buf, sslt_src%uWindFlag) 
        indV = fu_index(met_buf, sslt_src%vWindFlag) 
        if(error .or. indU < 1 .or. indV < 1)then
          call set_error('Failed to find near-surface wind fields','compute_emission_for_sea_salt')
          return
        endif

        if(sslt_src%uWindFlag == u_flag)then
          indHeight = fu_index(met_buf, height_flag)
          if(error .or. indHeight < 1)then
            call set_error('Failed to find height field','compute_emission_for_sea_salt')
            return
          endif
        else
          indHeight = int_missing
        endif

        indZ0 = fu_index(met_buf, surface_roughness_meteo_flag) 
        if(error .or. indZ0 < 1)then
          call set_error('Failed to find z0 field','compute_emission_for_sea_salt')
          return
        endif
        indFricVel = fu_index(met_buf, friction_velocity_flag) 
        if(error .or. indFricVel < 1)then
          call set_error('Failed to find 10m windspeed field','compute_emission_for_sea_salt')
          return
        endif
        !
        ! Area of white caps is a fraction of open water (minus the coastal zone, reflected
        ! in min open water fraction). See notebook 7 p.45 or papers of Monahan et al 
        ! (Oceanic Whitecaps, 1986) or Martensson et al (JGR, 2003).
        ! Note: 10m wind <= 38.7 m/s, otherwise white caps area fraction > 1
        !
        do iy = 1, ny_dispersion
          do ix = 1, nx_dispersion

            iDisp = ix+(iy-1)*nx_dispersion

            if(pSrcMask(iDisp) < 0.001)cycle

            iMeteo =  fu_grid_index(nx_meteo, ix, iy,  pHorizInterpMet2DispStruct)

            if(sslt_src%ifIceFractionCounts) then
              frOpenWater = (1. - met_buf%p2d(indIceFr)%present%ptr(iMeteo))
            else 
              frOpenWater = 1.
            end if

            frOpenWater = frOpenWater * pSrcMask(iDisp)

            if (frOpenWater > 0.) then
              u_star = fu_get_value(met_buf%p2d(indFricVel), nx_meteo, ix, iy, 1., &
                                  & pHorizInterpMet2DispStruct, ifHorizInterp)
!              u_10m = fu_get_value(met_buf%p2d(indW10m), nx_meteo, ix, iy, 1., &
!                                 & pHorizInterpMet2DispStruct, ifHorizInterp)
              u = fu_get_value(met_buf, indU, nx_meteo, ix, iy, 1, &
                             & pHorizInterpMet2DispStruct, VertInterpStruct_missing, &
                             & ifHorizInterp, .false.)
              v = fu_get_value(met_buf, indV, nx_meteo, ix, iy, 1, &
                             & pHorizInterpMet2DispStruct, VertInterpStruct_missing,&
                             & ifHorizInterp, .false.)
              if(sslt_src%uWindFlag == u_10m_flag)then
                height = 10.0
              else
                height = fu_get_value(met_buf, indHeight, nx_meteo, ix, iy, 1, &
                                    & pHorizInterpMet2DispStruct, VertInterpStruct_missing, &
                                    & ifHorizInterp, .false.)
              endif
              windspeed = sqrt(u*u + v*v)
              z0 = fu_get_value(met_buf%p2d(indZ0), nx_meteo, ix, iy, 1., &
                              & pHorizInterpMet2DispStruct, ifHorizInterp)
              if(z0 < 1.e-10)then
                z0 = 1.e-6
                if(iCount < 1000)then
                  call msg('gust: strange Z0=',z0)
                  iCount = iCount + 1
                endif
              endif
              fWhiteCaps = min(1.0, W_scale * frOpenWater * &
                                  & fu_effective_wind_pwr_gust(windspeed, height/z0, wind10m_power))
!                                  & fu_effective_wind_pwr_gust(u_10m, 10.0/z0, wind10m_power))
            else
              cycle   ! speedup
            end if

            if(fWhiteCaps < 1.0e-10)cycle  ! speedup...
            !
            ! Find the indices in the arrayFlux for the specific water temperature and 
            ! salinity - either dynamical or default
            !
            if(sslt_src%iWaterTempHandling /= static_value)then
              iWaterTemp = int((met_buf%p2d(indWaterTemprMet)%present%ptr(iMeteo) - fMinWaterTemp) / &
                             & stepWaterTemp + 0.5) + 1
              iWaterTemp = max(1,min(nWaterTemps,iWaterTemp))
            endif
            if(sslt_src%ifWaterSalinityCounts)then
              iSalinity = int((met_buf%p2d(indSalinityMet)%present%ptr(iMeteo)-fMinSalinity) / &
                            & stepSalinity + 0.5) + 1
              !if(iSalinity < 1 .or. iSalinity > nSalinities)then
              !  call msg('Forcing funny iSalinity:',iSalinity)
              if(iSalinity < 0 .or. iSalinity > nSalinities*2)then 
                ! something really wrong 
                call set_error('Strange salinity', 'compute_emission_for_sea_salt' )
                call msg('Strange salinity',iSalinity, met_buf%p2d(indSalinityMet)%present%ptr(iMeteo))
                return
              endif
              iSalinity = max(1,min(nSalinities,iSalinity))

!call msg('Salinity past:',met_buf%p2d(indSalinityMet)%past%ptr(iMeteo),iSalinity)
!call msg('Salinity present:',met_buf%p2d(indSalinityMet)%present%ptr(iMeteo),iSalinity)
!call msg('Salinity future:',met_buf%p2d(indSalinityMet)%future%ptr(iMeteo),iSalinity)
!call msg('Salinity value and index:',met_buf%p2d(indSalinityMet)%present%ptr(iMeteo),iSalinity)

            endif
            if(sslt_src%emisMethod == emisHybrid_spume_w10m)then
!              !
!              ! From Andreas, 1998, with approximation of numerical solution
!              !
!              fSpumeWindFactor = 3.5 * 41 * u_10m / (u_10m - 0.1) * &
!                               & exp((0.32*u_10m-0.5)/(0.008*u_10m+1.))
              !
              ! Andreas, 1998 modified: production for wind <6m/s is zero, <~10m suppressed
              !
              if(windspeed > 6.0)then
                if(if_low_spume)then
                  fSpumeWindFactor = 45 * exp((0.32*windspeed-0.5)/(0.008*windspeed+1.) - &   ! lowspume case
                                                & 20.0/(windspeed-5.)**2)
                else
                  fSpumeWindFactor = 3.5 * 45 * exp((0.32*windspeed-0.5)/(0.008*windspeed+1.) - &
                                                & 20.0/(windspeed-5.)**2)
                endif
              else
                fSpumeWindFactor = 0.0
              endif
            endif
            !
            ! Fill-in the emission map, not forgetting the number to mass convertion where needed
            !
            do iLev = 1, sslt_src%nLevsDispVert
              !
              ! First do the check for the overlap: speed-up
              !
              if(abs(sslt_src%levFractDispVert(iLev)) < 1.0e-5)cycle  ! nothing for this dispersion layer

              fCellTotal = 0.0

              !Common factor for all the fluxes
              dtdxdyLevFrac = timestep_sec * ptrXSzDisp(iDisp) * ptrYSzDisp(iDisp)* sslt_src%levFractDispVert(iLev)
              
              do iSeaSalt = 1, sslt_src%nSpecies
                if(sslt_src%ifNumberFlux(iSeaSalt))then
                  !
                  ! Number flux
                  !
                  if(sslt_src%emisMethod == emisHybrid_spume_w10m)then
                    fTmp = (sslt_src%fluxPerModeNbr(iSeaSalt,iWaterTemp,iSalinity) * fWhiteCaps + &
                          & sslt_src%fluxPerModeSpumeNbr(iSeaSalt,iSalinity) * fSpumeWindFactor) * dtdxdyLevFrac
                  else
                     fTmp = sslt_src%fluxPerModeNbr(iSeaSalt,iWaterTemp,iSalinity) * fWhiteCaps * dtdxdyLevFrac
                  endif
                else
                  !
                  ! Volume flux, converted to mass
                  !
                  if(sslt_src%emisMethod == emisHybrid_spume_w10m)then
                    !
                    ! Sum-up spume and bubble
                    !
                    fTmp = (sslt_src%fluxPerModeMass(iSeaSalt,iWaterTemp,iSalinity) * fWhiteCaps + &
                          & sslt_src%fluxPerModeSpumeMass(iSeaSalt,iSalinity) * fSpumeWindFactor) * &
                          & dtdxdyLevFrac !NOTE: volume -> mass       * fPartDensity(iSeaSalt)
                  else
                    !
                    ! Bubble-only
                    !
                    fTmp = sslt_src%fluxPerModeMass(iSeaSalt,iWaterTemp,iSalinity) * fWhiteCaps * &
                         & dtdxdyLevFrac  !NOTE: volume -> mass       * fPartDensity(iSeaSalt)
                  endif
!call msg('ix,iy,iSalinity,emis:'+fu_str(ix)+','+fu_str(iy),iSalinity,fTmp)
                endif
                
                fCellTotal = fCellTotal + fTmp !* fCorrectionFactor

                emisMap%arM(sslt_src%adaptor%iSp(iSeaSalt),sslt_src%id_nbr,iLev,ix,iy) = &
                      & emisMap%arM(sslt_src%adaptor%iSp(iSeaSalt),sslt_src%id_nbr,iLev,ix,iy) + fTmp
                fMassInjected(sslt_src%adaptor%iSp(iSeaSalt)) = &
                                               & fMassInjected(sslt_src%adaptor%iSp(iSeaSalt)) + fTmp
                
                if (ifSpeciesMoment) then  ! only vertical moment, horizontal ones are irrelevant
                  mapCoordZ%arm(sslt_src%adaptor%iSp(iSeaSalt),sslt_src%id_nbr, ilev, ix,iy) = &
                       & mapCoordZ%arm(sslt_src%adaptor%iSp(iSeaSalt),sslt_src%id_nbr, ilev, ix,iy) + &
                       & fTmp * sslt_src%fzDisp(iLev)
                end if

                !NOTE: The halogens, which are the last nHalogen species in the list, are done in this
                !      same loop because their mass flux (converted to moles) is already calculated in
                !      the tables sslt_src%fluxPerModeMass(:,:,:) and sslt_src%fluxPerModeSpumeMass(:,:)
                !      when these tables were generated in subr init_emission_sea_salt. 
                
                emisMap%ifColumnValid(sslt_src%id_nbr,ix,iy) = .true.
                emisMap%ifGridValid(iLev,sslt_src%id_nbr) = .true.

              end do      ! nSpecies
              if (.not. ifSpeciesMoment) then
                mapCoordZ%arM(1,sslt_src%id_nbr, iLev, ix, iy) = &
                                                    & sslt_src%fzDisp(iLev) * fCellTotal + &
                                                    & mapCoordZ%arM(1,sslt_src%id_nbr, iLev, ix, iy)
              end if
            end do      ! iLev
          end do   ! ix dispersion
        end do   ! iy dispersion

      case default
        call msg('Unknown method for sea salt emission:',sslt_src%emisMethod)
        call set_error('Unknown method for sea salt emission','compute_emission_for_sea_salt')
        return

    end select  ! emisMethod

    !call free_work_array(fPartDensity) !RISTO NOTE: volume -> mass
    ifFirstTime = .false.

!do iSeaSalt = 1, sslt_src%nSpecies
!call msg('SSLT species index and emission:' + fu_str(iSeaSalt), &
!             & sslt_src%adaptor%iSp(iSeaSalt), fMassInjected(sslt_src%adaptor%iSp(iSeaSalt)))
!end do
  end subroutine compute_emission_for_sea_salt


  !********************************************************************************************

  logical function fu_sslt_emis_owned_quantity(sslt_src, quantity)
    !
    ! Checks whether the specific quantity is handled by the emission sources exclusively
    ! Such quantities may be related to emission species somewhere else
    !
    implicit none

    ! Imported parameters
    type(silam_sea_salt_source), intent(in) :: sslt_src
    integer, intent(in) :: quantity
    !
    ! The sea salt source does not have own quantities yet
    !
    select case(quantity)
      case default
        fu_sslt_emis_owned_quantity = .false.
    end select
  end function fu_sslt_emis_owned_quantity


  !*****************************************************************

  integer function fu_source_id_nbr_of_sslt_src(sslt_src)
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
    type(silam_sea_salt_source), intent(in) :: sslt_src

    ! Stupidity check
    if(.not.(sslt_src%defined == silja_false))then
      call set_error('Undefined source given','fu_source_id_nbr_of_sslt_src')
      return
    endif
    fu_source_id_nbr_of_sslt_src = sslt_src%id_nbr

  end function fu_source_id_nbr_of_sslt_src



  !*************************************************************************

  integer function fu_source_nbr_of_sslt_src(sslt_src)
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
    type(silam_sea_salt_source), intent(in) :: sslt_src

    ! Stupidity check: only firmly undefined source returns int_missing
    if(.not. (sslt_src%defined == silja_false))then
      fu_source_nbr_of_sslt_src = sslt_src%src_nbr
    else
      fu_source_nbr_of_sslt_src = int_missing
      call set_error('Undefined source given','fu_source_nbr_of_sslt_src')
      return
    endif

  end function fu_source_nbr_of_sslt_src


  !*************************************************************************

  subroutine typical_species_cnc_sslt(srcSeaSalt, species, nSpecies, arConc)
    !
    ! Guesses a typical level of concentration and divides it with the given accuracy factor
    !
    implicit none

    ! Imported parameters
    type(silam_sea_salt_source), intent(in) :: srcSeaSalt
    type(silam_species), dimension(:), pointer :: species
    integer, intent(out) :: nSpecies
    real, dimension(:), pointer :: arConc
    real, dimension(:), pointer :: fPartDensity

    ! Local variables
    integer :: iSpecies, iWaterTemp, iSalinity
    real :: fTotalFluxNbr, fTotalFluxMass, fMassMeanDiam
    type(Taerosol_mode) :: mode

    real, parameter :: fTypicalMassCnc = 2.0e-9  ! two microgram

    species => srcSeaSalt%species
    nSpecies = srcSeaSalt%nSpecies !RISTO: CHECK

    !
    ! Procedure: compute the flux for each of the sea salt species, and the total flux
    ! for the whole size range. The total typical sea salt concentration is a few micrograms
    ! Then the fraction of the mass flux decides on the fraction of the concentration.
    ! Strictly speaking, this is not true due to varying deposition but here we do not have too
    ! many options.
    !
    ! The total flux
    !
    call set_aerosol_mode(mode, '', 1.e-8, 1.e-5, 1.e-6, 1.e3, fixed_diameter_flag, 1)
    if(error)return

    !
    ! Particle masses for each size class: !!!! dry particle !!!!
    !
    fPartDensity => fu_work_array()
    if(error)return

    do iSpecies = 1, srcSeaSalt%nSpeciesSeaSalt !NOTE: Omit the halogens
       fPartDensity(iSpecies) = fu_dry_part_density(srcSeaSalt%species(iSpecies)%material) !RISTO NOTE: volume -> mass
    end do
    if(error)return

    select case(srcSeaSalt%emisMethod)
      case(emisHybrid_w10m, emisHybrid_spume_w10m)
        call sea_salt_flux4mode_hybrid(mode, &
                                     & 285.0, &
                                     & 0.03, &
                                     & fTotalFluxNbr, &
                                     & fTotalFluxMass, & !RISTO NOTE: volume -> mass
                                     & fMassMeanDiam, &
                                     & fPartDensity(1))
      case default
        call msg('Unknown emission method:',srcSeaSalt%emisMethod)
        call set_error('Unknown emission method','init_emission_sea_salt')
        return
    end select
    !
    ! With the total flux known, get the typical concentrations
    !
    iWaterTemp = int((285 - fMinWaterTemp) / stepWaterTemp + 0.5) + 1
    iWaterTemp = max(1,min(nWaterTemps,iWaterTemp))
    iSalinity = int((0.03-fMinSalinity) / stepSalinity + 0.5) + 1
    iSalinity = max(1,min(nSalinities,iSalinity))

    !do iSpecies = 1, srcSeaSalt%nSpecies
    do iSpecies = 1, srcSeaSalt%nSpeciesSeaSalt  !NOTE: Omit the halogens!
      arConc(iSpecies) = fTypicalMassCnc / 2. * &
                       & min(1.0, srcSeaSalt%fluxPerModeMass(iSpecies,iWaterTemp,iSalinity) / & !RISTO NOTE: volume -> mass
                                & fTotalFluxMass)
    end do

    call free_work_array(fPartDensity) !RISTO NOTE: volume -> mass

  end subroutine typical_species_cnc_sslt


  !*************************************************************************

  function fu_sea_salt_source_name(sslt_src)result(chNm)
    implicit none
    type(silam_sea_salt_source), intent(in) :: sslt_src
    character(len=clen) :: chNm
    chNm = sslt_src%src_nm
  end function fu_sea_salt_source_name


  !*************************************************************************

  subroutine report_sea_salt_src(sslt_src)
    implicit none
    type(silam_sea_salt_source), intent(in) :: sslt_src
    integer :: iSpecies
    call msg('------------------ sea salt source report -----------------')
    if(sslt_src%sector_nm /= '')then
      call msg('Sea salt source' + sslt_src%src_nm + '_' + sslt_src%sector_nm)
    else
      call msg('Sea salt source' + sslt_src%src_nm)
    endif
    call msg('Species:')
    !do iSpecies = 1, sslt_src%nSpecies
    do iSpecies = 1, sslt_src%nSpeciesSeaSalt  !NOTE: Omit the halogens!
      call report(sslt_src%species(iSpecies))
    end do
    select case(sslt_src%emisMethod)
      case(emisHybrid_w10m)
        call msg('Hybrid mechanism, no spume')
      case (emisHybrid_spume_w10m)
        call msg('Hybrid mechanism, with spume')
        if(if_low_spume)then
          call msg('Low-spume sea salt mechanism')
        else
          call msg('Full-spume sea salt mechanism')
        endif
      case default
        call set_error('Unknown emission mechanism:'+fu_str(sslt_src%emisMethod),'report_sea_salt_src')
        return
    end select

    call msg('Wind quantities:' + fu_quantity_short_string(sslt_src%uWindFlag) + ',' + & 
                                & fu_quantity_short_string(sslt_src%vWindFlag)) 
      
    if(sslt_src%ifWaterSalinityCounts)then
      call msg('Water salinity is accounted for')
    else
      call msg('Water salinity is NOT accounted for')
    endif
    call msg('default water salinity:',sslt_src%defaultWaterSalinity)
    if(sslt_src%ifIceFractionCounts)then
      call msg('Ice fraction is accounted for')
    else
      call msg('Ice fraction is NOT accounted for')
    endif
    select case(sslt_src%iWaterTempHandling)
      case(static_value)
        call msg('Water temperature: static value')
        call msg('default awter temperature:', sslt_src%defaultWaterTempr)
      case(static_climatology)
        call msg('Water temperature: static climatology')
      case(dynamic_map)
        call msg('Water temperature: dynamic')
      case default
        call set_error('Unknown water temperature handling:' + fu_str(sslt_src%iWaterTempHandling),'report_sea_salt_src')
        return 
    end select
    
    call msg('Min open water fraction:',sslt_src%minOpenWaterFraction)
    
    call msg('------------------ end sea salt source report -----------------')
    

  end subroutine report_sea_salt_src


  !*************************************************************************
  !*************************************************************************
  !
  ! Private functions computing the sea salt emission fluxes
  !
  !*************************************************************************
  !*************************************************************************


  subroutine sea_salt_flux4mode_hybrid(aerMode, &        ! definition of the spectrum band
                                     & fWaterTempr, &    !  [K]
                                     & fWaterSalinity, & ! [ratio]
                                     & fNbrFlux, fMassFlux, fMassMeanDiam, & !RISTO NOTE: volume -> mass
                                     & dryPartDensity, & ! dry part density kg/m3
                                     & iHalogen)         ! Optional, Halogen index 
    !
    ! Computes the number and mass emission flux of the sea salt for the pre-defined
    ! size spectrum using the combined Monahan-Martensson method. Both water
    ! salinity and temperature are taken into account.
    !
    ! For halogen species calculates the contribution of the corresponding sea-salt
    ! mode to the halogen emissions.
    !
    implicit none

    ! Imported parameters
    type(Taerosol_mode), intent(inout) :: aerMode
    real, intent(in) :: fWaterTempr, fWaterSalinity, dryPartDensity
    !real, intent(out) :: fNbrFlux, fVolFlux, fMassMeanDiam
    real, intent(out) :: fNbrFlux, fMassMeanDiam
    real, intent(inout) :: fMassFlux !RISTO NOTE: volume -> mass
    integer, optional, intent(in) :: iHalogen

    ! Local variables
    real :: fMinD, fMaxD, fD, fD_um, fFluxDens, fFluxDens_prev, fIntegrStep, fIntegrStep_um, fTmp
    real :: fVolFlux, EF, EF_prev
    integer :: iHalog
    type(silam_sp) :: sp

    real, parameter :: fAbsoluteMin = 3.0e-9  ! 3 nm
    real, parameter :: fAbsoluteMax = 1.0e-4  ! 100 um

    !In order to calculate the halogen emission from sea salt (for each mode) we use iHalog to
    !indicate the halogen species.
    if (present(iHalogen)) then
       iHalog = iHalogen
    else !Normal sea-salt case
       iHalog = 0
    end if
    
    !
    ! Algorithm: scan each size class and make an integration from min to max diameter of the
    ! number flux. Mean diameter can be then computed too and set to the aerosol back.
    !
!    sp%sp=> fu_work_string()

    !
    ! Get the parameters of the mode
    !
    fMinD = fu_min_d(aerMode)
    fMaxD = fu_max_d(aerMode)

!if(if_25_Degrees_Reference)then
!call msg('========> Reference temperature is 25 degrees')
!else
!call msg('========> Reference temperature is 15 degrees')
!endif

    !
    ! Now - careful. We have two ranges: Martensson and Monahan. Their border is 1 um. Both are
    ! valid within their own limits. Check all
    !
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

      !
      ! Get the integration step
      !
      fIntegrStep = 0.001*(fMaxD - fMinD)  ! At least 1000 steps over the range

      fD = fMinD
      fNbrFlux = 0.
      fVolFlux = 0.
      fMassMeanDiam = 0.

      if(if_25_Degrees_Reference)then
        fFluxDens = fu_shape_function(fD) * &
                  & fu_water_tempr_corr_to_25deg(fWaterTempr, fD) * &
                  & fu_salinity_corr(fWaterSalinity, fD)
      elseif(if_15_Degrees_Reference)then
        fFluxDens = fu_shape_function(fD) * &
                  & fu_water_tempr_corr_to_15deg(fWaterTempr, fD) * &
                  & fu_salinity_corr(fWaterSalinity, fD)
      else
        fFluxDens = fu_shape_function(fD) * &
                  & fu_water_tempr_corr_to_20deg(fWaterTempr, fD) * &
                  & fu_salinity_corr(fWaterSalinity, fD)
      endif

      !For halogens we use emission factor EF to scale the sea-salt emissions to Br/Cl emissions.
      !Additionally, we convert the emissions to moles.
      !Note that for each sea-salt specie brings it own contribution to every halogen. Therefore,
      !one needs to sum the contribution from every sea-salt specie. 
      fD_um = fD * 1e6 !Diameter in micrometers (needed for fu_depletion_factor_bromine)
      select case(iHalog)
      case(1) !Bromine (currently using Br2, CHECK)
         !fScale = 0.00223 / fu_molar_mass(fu_get_material_ptr('Br2')) ! kg Br / kg SSLT -> mole Br / kg SSLT
         !Note that the bromine mass with respect to mass of NaCl is actually 0.00246 and not 0.00223 as
         !was stated e.g. in Yang et al JGR Vol. 110, D23311 (2008) because the later smaller value is calculated
         !using total mass of Na and Cl. The chlorine mass is more than what the formula NaCl states, i.e. the mole
         !amount of Cl is larger than the mole amount of Na. See Table 1, in Sander et al., ACP 3, 1301 (2003). 
         EF = 0.00246*fu_depletion_factor_bromine(fD_um)/fu_mole_mass(fu_get_material_ptr('Br2')) ! mole(Br)/kg(sea-salt)
      case(2) !Clorine (currently using HCl, CHECK)
         !The total mass fraction of chlorine compared to mass of ideal sea-salt [NaCl] is 0.7069.
         !However, the sea-salt contains more moles of Cl than Na. The mass of this excess Cl is
         !0.1003 times the mass of NaCl. See e.g. Table 1, in Sander et al., ACP 3, 1301 (2003) for
         !the relative amounts in sea-water.
         !The form of chlorine emissions is somewhat unknown but it is more likely that HCl is released
         !by displacement of less-volatile stronger acids, such as nitric and sulfuric acid than
         !the release of Cl2 gas.  
         !The emission of Cl is less efficient than Br. See the above Sander et al. The depletion factor
         !is close to zero but the resolution does not allow a proper value. We take it to be 0.05.
         !EF = 0.1003/fu_mole_mass(fu_get_material_ptr('Cl2'))     !For Cl2 OLD FixMe
         !EF = 0.05*0.7069/fu_mole_mass(fu_get_material_ptr('Cl')) !For HCl FixMe
         EF = 0.0 !OMIT the Chlorines, Use GEIA source for chlorine species
      case default !Normal sea-salt case
         EF = 1.0
      end select

      do while(fD < fMaxD)
        fIntegrStep = min(fIntegrStep, fMaxD - fD)
        fD = fD + fIntegrStep
        fFluxDens_prev = fFluxDens
        EF_prev = EF
        
        fD_um = fD * 1e6
        fIntegrStep_um = fIntegrStep * 1e6
        !
        ! Number-flux density
        !
        if(if_25_Degrees_Reference)then
          fFluxDens = fu_shape_function(fD) * &
                    & fu_water_tempr_corr_to_25deg(fWaterTempr, fD) * &
                    & fu_salinity_corr(fWaterSalinity, fD)
        elseif(if_15_Degrees_Reference)then
          fFluxDens = fu_shape_function(fD) * &
                    & fu_water_tempr_corr_to_15deg(fWaterTempr, fD) * &
                    & fu_salinity_corr(fWaterSalinity, fD)
        else
          fFluxDens = fu_shape_function(fD) * &
                    & fu_water_tempr_corr_to_20deg(fWaterTempr, fD) * &
                    & fu_salinity_corr(fWaterSalinity, fD)
        endif
        !
        ! Integrate the number flux. Note that the shape function is per um.
        !
        fNbrFlux = fNbrFlux + fIntegrStep_um * 0.5 * (fFluxDens_prev + fFluxDens)  ! for number flux
        !
        ! Take into account the emission factor for halogens:
        select case(iHalog)
        case(1) !Bromine (currently using Br2, CHECK)
           EF = 0.00246*fu_depletion_factor_bromine(fD_um)/fu_mole_mass(fu_get_material_ptr('Br2')) ! mole(Br)/kg(sea-salt)
        case(2) !Clorine (currently using HCl, CHECK)
           !EF = 0.1003/fu_mole_mass(fu_get_material_ptr('Cl2'))     !For Cl2 OLD FixMe
           !EF = 0.05*0.7069/fu_mole_mass(fu_get_material_ptr('Cl')) !For HCl FixMe
           EF = 0.0 !OMIT the Chlorines, Use GEIA source for chlorine species
        case default !Normal sea-salt case
           EF = 1.0
        end select
        !
        ! And integrate the mass (volume) flux
        !     
        fTmp = fIntegrStep_um * 0.5 * 0.5235987756 * (EF_prev*fFluxDens_prev * (fD_um-fIntegrStep_um)**3 + &
                                                    & EF*fFluxDens * fD_um**3)
!sp%sp=> fu_work_string()
!if (error) return
!write(unit = sp%sp,fmt=*)'Dp_um, Tw, Sw, Tcorr, Scorr, f1, f2, f3, flux for D, m-2 sec-1 um-1:', &
!                       & fD_um, fMeanD, &
!                       & fWaterTempr, & ! [K]
!                       & fWaterSalinity, &
!                       & fWTcorr, fScorr, f1, f2, f3, fMonahanFunc
!call msg(sp%sp)
!call free_work_array(sp%sp)

        fVolFlux = fVolFlux + fTmp 
        fMassMeanDiam = fMassMeanDiam + fTmp * (fD_um-0.5*fIntegrStep_um)
      end do  ! cycle over diameter

      fMassMeanDiam = 1.e-6 * fMassMeanDiam / fVolFlux  ! from um to m
      fVolFlux = fVolFlux * 1.0e-18                     ! from um^3 to m^3
    endif

    !For halogens one needs to sum the volume flux from each sea-salt species.
    !Also for halogens the number flux and mass mean diameter are not applicable
    if (iHalog > 0) then
       fMassFlux = fMassFlux+dryPartDensity*fVolFlux !RISTO NOTE: volume -> mass
       fNbrFlux = 0.0
       fMassMeanDiam = 0.0
    else !standard sea-salt species
       fMassFlux = dryPartDensity*fVolFlux !RISTO NOTE: volume -> mass
    end if
    
  end subroutine sea_salt_flux4mode_hybrid


  !************************************************************************

  real function fu_shape_function(fD)
    !
    ! The main shape fucntion.
    !
    implicit none
    
    ! Imported parameter
    real, intent(in) :: fD  ![m]
    
    ! Local variables
    real :: fD_um, f1, f2, f3

    ! Locla parameters
    real, parameter :: dFdD_scale = 1.0e6 !689733 
    real, parameter :: power_of_D = 3

    fD_um = fD * 1.e6

    if(if_25_Degrees_Reference)then
!      f1 = 1./ (fD_um+exp(-20.*fD_um)/(20.+exp(-20.*fD_um)))**power_of_D 
!      f2 = dFdD_scale * (1+0.05*fD_um**1.05)
!      f3 = 10**(power_of_10*exp(-(log10(fD_um)/0.8) * (log10(fD_um)/0.8)))

      !
      ! Applied in v.4.5.4 together with Zhang dry deposition. Good results, huge dry dep
      !
!      f1 = exp(-0.12 / fD_um)/(0.05 + exp(- 0.1 / fD_um))
!      f2 = dFdD_scale * (1+0.05*fD_um**1.05) / fD_um ** power_of_D  
!      f3 = 10**(1.6 * exp(-(log10(1.1 * fD_um)/0.8) * (log10(1.1 * fD_um)/0.8)))

!      !
!      ! Fit after correcting the N <-> F error 
!      ! Explanation 30.10.2012: the error was in fitting N-distribution of M03,
!      ! whereas the flux should be fitted. This results in a constant factor for all
!      ! size classes, which was estimated to be 0.2: F=N*f_air_in/area, where
!      ! f_air_in = 6e-5 m3/sec (3.6l/min,M03). With area of emitting surface (white caps) 1e-3 m2
!      ! or less, we get something 0.06 or more. To match the M03 figure 11, 0.2 is good.
!      ! Then refitting the lower range to M03 but keeping the coarser range to M86, we get:
!      !
!      f1 = exp(-0.12 / fD_um)/(0.3 + exp(- 0.15 / fD_um))
!      f2 = dFdD_scale * (1+0.05*fD_um**1.05) / fD_um ** power_of_D  
!      f3 = 10**(1.05 * exp(- (0.1-log10(1.1*fD_um)/0.9)**2 ))
      !
      ! ... and the next iteration
      !
!      f1 = exp(-0.09 / (0.003+fD_um))/(2 + exp(-5. / fD_um))
!      f2 = dFdD_scale * (1+0.05*fD_um**1.05) / fD_um ** power_of_D  
!      f3 = 10**(1.05 * exp(- (0.27-log10(fD_um)/1.1)**2 ))

      !
      ! ... and an attempt to increase AOD by raising submicron part of the spectrum ~40%
      !
      f1 = exp(-0.095 / (0.001+fD_um))/(1.7 + exp(-10. / fD_um))
      f2 = dFdD_scale * (1+0.05*fD_um) / fD_um ** power_of_D  
      f3 = 10**(1.05 * exp(- (0.1-log10(fD_um)/0.95)**2 ))

    elseif(if_15_Degrees_Reference)then
      !
      ! 15 degrees as the reference fit.
      !
      f1 = exp(-0.1 / fD_um)/(0.08 + exp(-0.1 / fD_um))
      f2 = 6.0e5 * (1+0.05*fD_um**1.05) / fD_um ** 3
      f3 = 10**(1.35 * exp(- (0.1-log10(fD_um)/0.95)**2 ))
    else
      !
      ! Default is 20 degrees
      !
! Error: forgot to correct M03 N->flux scaling
!      f1 = exp(-0.1 / fD_um)/(0.15 + exp(-0.1 / fD_um))
!      f2 = 4.7e5 / (0.001 * fD_um**2 + fD_um) ** 3        ! large D, cut-off
!      f3 = 10**(1.4 * exp(- (0.2-log10(fD_um)/1.1)**2 ))
      
      !
      ! This M03 20C fit. However, scaling factor N->F is 3.7 instead of 5.
      ! With such factor it is closer to M86. 20C is then taken as a square root of 15C and 25C
      !
!      f1 = (1. + 0.05 * fD_um) * exp(-0.1 / fD_um)/(0.35 + exp(-0.1 / fD_um))
!      f2 = 6.e5 / (0.0001 * fD_um**2 + fD_um) ** 3        ! large D, cut-off
!      f3 = 10**(1.19 * exp(- (0.25-log10(fD_um)/0.9)**2 ))

      !
      ! And this is M03 20C fit with correct, scaling factor 5 for N->F.
      ! With such factor it is somewhat further from M86 but hopefully closer to observations. 
      ! 20C is then taken as a square root of 15C and 25C
      !
      f1 = (1. + 0.05 * fD_um) * exp(-0.11 / fD_um)/(0.4 + exp(-0.2 / fD_um))
      f2 = 6.e5 / (0.0001 * fD_um**2 + fD_um) ** 3        ! large D, cut-off
      f3 = 10**(1.19 * exp(- (0.35-log10(fD_um)/0.8)**2 ))

    endif  ! if_25_Degrees_Reference
      
    fu_shape_function = f1 * f2 * f3
    
  end function fu_shape_function


  !************************************************************************

  real function fu_salinity_corr(fWaterSalinity, fDPart)
    !
    ! Computes the correction coefficient for aerosol production with regad to 
    ! the water salinity. Correction was taken from the Martensson's paper by
    ! fitting the ratio of fluxes for S=0.033 and S=0.0092 / S=0. The correction 
    ! safe: for large aerosol it tends to zero thus reducing large particle flux
    ! from low-salinity water, while down to 10 nm the correction is reasonable
    ! and does not exceed 2-3. It is better not to go further down due to no
    ! data. Here a limit of the correction of 10 is forced, which correcponds to
    ! well sub-nanometer scale.
    !
    implicit none

    ! Imported parameters
    real, intent(in) :: fWaterSalinity, fDPart ! salinity,[fraction]; aerosol size,[m]

    !
    ! We have a limit correction for fresh water and intermediate coefficient
    ! for low-salinity as in Baltic sea. The rest will be obtained as a linear
    ! combination of these two.
    !
    ! Specific parameters of fitting - see notebook 7, pg. 56
    ! Obs: 42.0168 = 1. / (0.033-0.0092)
    !      108.7 = 1. / 0.0092
    !
    if(fWaterSalinity >= 0.033)then ! very salty water, assume equal to oceanic

      fu_salinity_corr = 1.

    elseif(fWaterSalinity >= 0.0092)then ! in-between ocean (corr=1) and Baltic Sea

      fu_salinity_corr = (fWaterSalinity - 0.0092)*42.0168 + &
                       & (0.033 - fWaterSalinity)*42.0168 * 6.0e-6 * fDPart**(-0.7142)

    elseif(fWaterSalinity >= 0.)then  ! low-saline water, beware of small but negative salinity

      fu_salinity_corr = max(0.,fWaterSalinity)*108.7 * 6.0e-6 * fDPart**(-0.7142) + &
                       & (0.0092 - max(0.,fWaterSalinity))*108.7 * 4.0e-15 * fDPart**(-1.6942)
    else

      fu_salinity_corr = 4.0e-15 * fDPart**(-1.6942)
    endif
    !
    ! Just a precaution:
    !
    if(fu_salinity_corr > 10.) fu_salinity_corr = 10.

  end function fu_salinity_corr


  !************************************************************************

  real function fu_water_tempr_corr_to_25deg(fWaterTempr, fDPart)
    !
    ! Computes the correction coefficient for aerosol production with regard to 
    ! the water temperature. Correction was taken from the Martensson's paper by
    ! fitting the ratio of fluxes for T=15C and 25, 5, -2.  
    !
    implicit none

    ! Imported parameters
    real, intent(in) :: fWaterTempr, fDPart ! salinity,[fraction]; aerosol size,[m]

    ! Local variables
    real :: fD_log

    !
    ! We have a limit correction for -2 C water, intermediate coefficient for 5 and 15C and 
    ! a unity for 25 C. The rest will be obtained as their linear combination.
    !
    ! Specific parameters of fitting - see notebook 7, pg. 56
    ! Obs: 0.1 = 1. / (298 - 288)
    !      0.142857 = 1. / (278 - 271)
    !

    fD_log = log10(fDPart)

    if(fWaterTempr >= 298.)then ! very warm water, assume equal to 25C

      fu_water_tempr_corr_to_25deg = 1.

    elseif(fWaterTempr >= 288.)then ! in-between 25C (corr=1) and 15C

    fu_water_tempr_corr_to_25deg = (fWaterTempr - 288.)*0.1 + &
                                 & (298. - fWaterTempr) * 0.1 * 0.0031 * exp (-0.8353 * fD_log)

    elseif(fWaterTempr >= 278.)then  !in-between 5 and 15

      fu_water_tempr_corr_to_25deg = (fWaterTempr - 278.)*0.1 * &
                                    & 0.0031 * exp (-0.8353 * fD_log) + &
                                   & (288.-fWaterTempr)* 0.1 * 7.99e-7 * exp(-2.02*fD_log)

    elseif(fWaterTempr >= 271.)then

      fu_water_tempr_corr_to_25deg = (fWaterTempr-271.)*0.142857 * 7.99e-7 * exp(-2.02*fD_log) + &
                                   & (278.-fWaterTempr)*0.142857 * 1.6e-7 * exp(-2.21*fD_log)

    else

      fu_water_tempr_corr_to_25deg = 1.6e-7 * exp(-2.21* fD_log)

    endif

    !
    ! Just a precaution:
    !
    if(fu_water_tempr_corr_to_25deg > 10.) fu_water_tempr_corr_to_25deg = 10.

  end function fu_water_tempr_corr_to_25deg


  !************************************************************************

  real function fu_water_tempr_corr_to_15deg(fWaterTempr, fDPart)
    !
    ! Computes the correction coefficient for aerosol production with regard to 
    ! the water temperature. Correction was taken from the Martensson's paper by
    ! fitting the ratio of fluxes for T=15C and 25, 5, -2. The correction 
    ! not very safe: for large aerosol in hot water it grows, thus potentially 
    ! causing trouble. Therefore, it is artificially stopped at 5 um. 
    !
    implicit none

    ! Imported parameters
    real, intent(in) :: fWaterTempr, fDPart ! water temperature,[K]; aerosol size,[m]

    ! Local variables
    real :: fD_um

    !
    ! We have a limit correction for -2 C water, intermediate coefficient for 5 and 25C and 
    ! a unity for 15 C. The rest will be obtained as their linear combination.
    !
    ! Specific parameters of fitting - see refitting_martensson_15_deg_201211.xls and notebook 7, pg. 56
    ! Obs: 0.1 = 1. / (298 - 288)
    !      0.142857 = 1. / (278 - 271)
    !

    fD_um = min(5., fDPart * 1.e6) ! stop at 5 um

    if(fWaterTempr >= 298.)then ! very warm water, assume equal to 25C

      fu_water_tempr_corr_to_15deg = 2. * fD_um ** 0.31

    elseif(fWaterTempr >= 288.)then ! in-between 25C and 15C (where corr=1)

      fu_water_tempr_corr_to_15deg = (fWaterTempr - 288.)*0.1 * 2. * fD_um ** 0.31 + &
                                   & (298. - fWaterTempr)*0.1

    elseif(fWaterTempr >= 278.)then  !in-between 5 and 15

      fu_water_tempr_corr_to_15deg = (fWaterTempr - 278.)*0.1 + &
                                   & (288.-fWaterTempr)*0.1 * 0.35 * (fD_um+0.2) ** (-1.2)

    elseif(fWaterTempr >= 271.)then

      fu_water_tempr_corr_to_15deg = (fWaterTempr-271.)*0.142857 * 0.35 * (fD_um+0.2) ** (-1.2) + &
                                   & (278.-fWaterTempr)*0.142857 * 0.3 * (fD_um+0.3) ** (-1.4)

    else

      fu_water_tempr_corr_to_15deg = 0.3 * (fD_um+0.3) ** (-1.4)

    endif

  end function fu_water_tempr_corr_to_15deg


  !************************************************************************

  real function fu_water_tempr_corr_to_20deg(fWaterTempr, fDPart)
    !
    ! Computes the correction coefficient for aerosol production with regard to 
    ! the water temperature. Correction was taken from the Martensson's paper by
    ! fitting the ratio of fluxes for T=15C and 25, 5, -2. The correction 
    ! not very safe: for large aerosol in hot water it grows, thus potentially 
    ! causing trouble. Therefore, it is artificially stopped at 5 um. 
    !
    implicit none

    ! Imported parameters
    real, intent(in) :: fWaterTempr, fDPart ! water temperature,[K]; aerosol size,[m]

    ! Local variables
    real :: fD_um

    !
    ! We have a limit correction for -2 C water, intermediate coefficient for 5 and 25C and 
    ! a unity for 15 C. The rest will be obtained as their linear combination.
    !
    ! Specific parameters of fitting - see refitting_martensson_15_deg_201211.xls and notebook 7, pg. 56
    ! Obs: 0.1 = 1. / (298 - 288), 
    !      0.2 = 1. / (298 - 293) = 1. / (293 - 288)
    !      0.142857 = 1. / (278 - 271)
    !

    fD_um = min(5., fDPart * 1.e6) ! stop at 5um

    if(fWaterTempr >= 298.)then ! very warm water, assume equal to 25C

      fu_water_tempr_corr_to_20deg = 1.4 * fD_um**2.13 / (fD_um + 0.002)**2

    elseif(fWaterTempr >= 293.)then ! in-between 20C (corr=1) and 25C

      fu_water_tempr_corr_to_20deg = (fWaterTempr - 293.)*0.2 * 1.4 * fD_um**2.13 / (fD_um+0.002)**2 + &
                                   & (298. - fWaterTempr)*0.2

    elseif(fWaterTempr >= 288.)then ! in-between 15C and 20C (corr=1)

      fu_water_tempr_corr_to_20deg = (fWaterTempr - 288.)*0.2 + &
                                   & (293. - fWaterTempr)*0.2 * 0.75 * fD_um**1.85 / (fD_um+0.001)**2

    elseif(fWaterTempr >= 278.)then  !in-between 5 and 15

      fu_water_tempr_corr_to_20deg = (fWaterTempr - 278.)*0.1 * 0.75 * fD_um**1.85 / (fD_um+0.001)**2 + &
                                   & (288.-fWaterTempr)*0.1 * 0.22*fD_um**1.1 / (fD_um+0.012)**2

    elseif(fWaterTempr >= 271.)then  ! in-between -2 and 5

      fu_water_tempr_corr_to_20deg = (fWaterTempr-271.)*0.142857 * 0.22*fD_um**1.1 / (fD_um+0.012)**2 + &
                                   & (278.-fWaterTempr)*0.142857 * 0.13 * fD_um / (fD_um+0.016)**2

    else  ! colder than -2 - take -2.

      fu_water_tempr_corr_to_20deg = 0.13 * fD_um / (fD_um+0.016)**2

    endif

  end function fu_water_tempr_corr_to_20deg


  !**********************************************************************************

  
  subroutine sea_salt_flux4mode_spume(aerMode, &        ! definition of the spectrum band
                                    & fWaterSalinity, &
                                    & fNbrFlux, fMassFlux, fMassMeanDiam, & !RISTO NOTE: volume -> mass
                                    & dryPartDensity, & ! dry part density kg/m3
                                    & iHalogen)         ! Optional, Halogen index 
    !
    ! Computes the number and mass emission flux of the sea salt for the pre-defined
    ! size spectrum using the combined Monahan-Martensson method. Both water
    ! salinity and temperature are taken into account
    !
    implicit none

    ! Imported parameters
    type(Taerosol_mode), intent(inout) :: aerMode
    real, intent(in) :: fWaterSalinity, dryPartDensity
    !real, intent(out) :: fNbrFlux, fVolFlux, fMassMeanDiam
    real, intent(out) :: fNbrFlux, fMassMeanDiam
    real, intent(inout) :: fMassFlux !RISTO NOTE: volume -> mass
    integer, optional, intent(in) :: iHalogen

    ! Local variables
    real :: fMinD, fMaxD, fD, fD_um, fFluxDens, fFluxDens_prev, fIntegrStep, fIntegrStep_um, fTmp
    real :: fVolFlux, EF, EF_prev
    integer :: iHalog
    type(silam_sp) :: sp

    real, parameter :: fAbsoluteMin = 1.0e-6  ! 1 um
    real, parameter :: fAbsoluteMax = 1.0e-4  ! 100 um

    !In order to calculate the halogen emission from sea salt (for each mode) we use iHalog to
    !indicate the halogen species.
    if (present(iHalogen)) then
       iHalog = iHalogen
    else !Normal sea-salt case
       iHalog = 0
    end if
    
    !
    ! Algorithm: scan the size range and make an integration from min to max diameter of the
    ! number flux. Mean diameter can be then computed too and set to the aerosol back.
    ! Get the parameters of the mode
    !
    fMinD = fu_min_d(aerMode)
    fMaxD = fu_max_d(aerMode)

    !
    ! Now - careful. Spume particles strart from 10 um
    !
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

      !
      ! Get the integration step
      !
      fIntegrStep = 0.001*(fMaxD - fMinD)  ! At least 1000 steps over the range

      fD = fMinD
      fNbrFlux = 0.
      fVolFlux = 0.
      fMassMeanDiam = 0.

      fFluxDens = fu_shape_function_spume(fD) * &
                & fu_salinity_corr(fWaterSalinity, fD)

      !For halogens we use emission factor EF to scale the sea-salt emissions to Br/Cl emissions.
      !Additionally, we convert the emissions to moles (the conversion of NaCl from volume to kg is done elsewhere)
      !Note that for each sea-salt specie brings it own contribution to every halogen. Therefore,
      !one needs to sum the contribution from every sea-salt specie. 
      fD_um = fD * 1e6 !Diameter in micrometers (needed for fu_depletion_factor_bromine)
      select case(iHalog)
      case(1) !Bromine (currently using Br2, CHECK)
         !fScale = 0.00223 / fu_molar_mass(fu_get_material_ptr('Br2')) ! kg Br / kg SSLT -> mole Br / kg SSLT
         !Note that the bromine mass with respect to mass of NaCl is actually 0.00246 and not 0.00223 as
         !was stated e.g. in Yang et al JGR Vol. 110, D23311 (2008) because the later smaller value is calculated
         !using total mass of Na and Cl. The chlorine mass is more than what the formula NaCl states, i.e. the mole
         !amount of Cl is larger than the mole amount of Na. See Table 1, in Sander et al., ACP 3, 1301 (2003). 
         EF = 0.00246*fu_depletion_factor_bromine(fD_um)/fu_mole_mass(fu_get_material_ptr('Br2')) ! mole(Br)/kg(sea-salt)
      case(2) !Clorine (currently using HCl, CHECK)
         !The total mass fraction of chlorine compared to mass of ideal sea-salt [NaCl] is 0.7069.
         !However, the sea-salt contains more moles of Cl than Na. The mass of this excess Cl is
         !0.1003 times the mass of NaCl. See e.g. Table 1, in Sander et al., ACP 3, 1301 (2003) for
         !the relative amounts in sea-water.
         !The form of chlorine emissions is somewhat unknown but it is more likely that HCl is released
         !by displacement of less-volatile stronger acids, such as nitric and sulfuric acid than
         !the release of Cl2 gas.  
         !The emission of Cl is less efficient than Br. See the above Sander et al. The depletion factor
         !is close to zero but the resolution does not allow a proper value. We take it to be 0.05.
         !EF = 0.1003/fu_mole_mass(fu_get_material_ptr('Cl2'))     !For Cl2 OLD FixMe
         !EF = 0.05*0.7069/fu_mole_mass(fu_get_material_ptr('Cl')) !For HCl FixMe
         EF = 0.0 !OMIT the Chlorines, Use GEIA source for chlorine species
      case default !Normal sea-salt case
         EF = 1.0
      end select

!sp%sp=> fu_work_string()
!if (error) return
!write(unit = sp%sp,fmt=*)'Dp_um, Tw, Sw, Tcorr, Scorr, f1, f2, f3, flux for D, m-2 sec-1 um-1:', &
!                       & fD_um, fMeanD, &
!                       & fWaterTempr, & ! [K]
!                       & fWaterSalinity, &
!                       & fWTcorr, fScorr, f1, f2, f3, fMonahanFunc
!call msg(sp%sp)
!call free_work_array(sp%sp)

      do while(fD < fMaxD)
        fIntegrStep = min(fIntegrStep, fMaxD - fD)
        fD = fD + fIntegrStep
        fFluxDens_prev = fFluxDens
        EF_prev = EF

        fD_um = fD * 1e6
        fIntegrStep_um = fIntegrStep * 1e6

        fFluxDens = fu_shape_function_spume(fD) * fu_salinity_corr(fWaterSalinity, fD)

        !
        ! Integrate the number flux. Note that the shape function is per um.
        !
        fNbrFlux = fNbrFlux + fIntegrStep_um * 0.5 * (fFluxDens_prev + fFluxDens)  ! for number flux
        !
        ! Take into account the emission factor for halogens:
        select case(iHalog)
        case(1) !Bromine (currently using Br2, CHECK)
           EF = 0.00246*fu_depletion_factor_bromine(fD_um)/fu_mole_mass(fu_get_material_ptr('Br2')) ! mole(Br)/kg(sea-salt)
        case(2) !Clorine (currently using HCl, CHECK)
           !EF = 0.1003/fu_mole_mass(fu_get_material_ptr('Cl2'))     !For Cl2 OLD FixMe
           !EF = 0.05*0.7069/fu_mole_mass(fu_get_material_ptr('Cl')) !For HCl FixMe
           EF = 0.0 !OMIT the Chlorines, Use GEIA source for chlorine species
        case default !Normal sea-salt case
           EF = 1.0
        end select
        !
        ! And integrate the mass (volume) flux
        !        
        fTmp = fIntegrStep_um * 0.5 * 0.5235987756 * (EF_prev*fFluxDens_prev * (fD_um-fIntegrStep_um)**3 + &
                                                    & EF*fFluxDens * fD_um**3)

        !call msg('Diameter and volume flux',fD_um,fVolFlux)
        fVolFlux = fVolFlux + fTmp 
        fMassMeanDiam = fMassMeanDiam + fTmp * (fD_um-0.5*fIntegrStep_um) !Not applicable for halogens
      end do  ! cycle over diameter

      fMassMeanDiam = 1.e-6 * fMassMeanDiam / fVolFlux  ! from um to m
      fVolFlux = fVolFlux * 1.0e-18                     ! from um^3 to m^3

    endif

    !For halogens one needs to sum the volume flux from each sea-salt species.
    !Also for halogens the number flux and mass mean diameter are not applicable
    if (iHalog > 0) then
       fMassFlux = fMassFlux+dryPartDensity*fVolFlux !RISTO NOTE: volume -> mass
       fNbrFlux = 0.0
       fMassMeanDiam = 0.0
    else !standard sea-salt species
       fMassFlux = dryPartDensity*fVolFlux !RISTO NOTE: volume -> mass
    end if
    
  end subroutine sea_salt_flux4mode_spume


  !************************************************************************

  real function fu_shape_function_spume(fD)
    !
    ! The main shape fucntion for spume production
    !
    implicit none
    
    ! Imported parameter
    real, intent(in) :: fD  ![m]
    !
    ! Essentially Andreas, 1998, Journal of Physical Oceanlology, with some speed-up
    ! Note that the wind dependence is in dynamic part. Here only Dp dependence, which is simple
    ! Obs conversion to micrometres
    !
    if(fD < 1.0e-6)then  ! we consider submicron particles outright impossible in the spume
      fu_shape_function_spume = 0.0
    elseif(fD < 37.5e-6)then
      fu_shape_function_spume = 1.e-6 / fD
    elseif(fD < 100.0e-6)then
      fu_shape_function_spume = 681.1712 * (fD*1.e6) ** (-2.8) ! 37.5**1.8
    elseif(fD < 250.0e-6)then
      fu_shape_function_spume = 1.711e13 * (fD*1.e6) ** (-8) ! 37.5**1.8 * 100**5.2
    else
      fu_shape_function_spume = 0.0
    endif    
    
  end function fu_shape_function_spume


  !**********************************************************************************


  real function fu_depletion_factor_bromine(fD_um) result(DF)
    !Calculates the sea-salt aerosol diameter dependent depletion factor DF(d) for bromine.
    !See, X. Yang, J.A. Pyle, R.A. Cox, Geophysical Research Letters, 35 L16815 (2008).
    !The Table 1, contains the values where this fitting fuction is based.
    !This table is based on measuerements by Sander et al., Atmos. Chem. Phys. 3, 1301-1336 (2003).
    !However, Sander typically use enrichment factor EF=1-DF.
    !
    !See also Yang et al., J. Geophys. Research 101, D23311 (2005), where for constant DF
    !the bromine emissions are E = P*Ra*DF, where
    ! P is the sea-salt production rate (kg/m2*s),
    ! Ra is the sea-salt Br/NaCl mass ratio Ra= 0.00223 kg(Br)/kg(NaCl) and
    ! DF is the depletion factor for bromine
    
    implicit none
    real :: fD_um !aerosol diameter, micrometes

    !Fitted shape function using the Table 1 from Yang, et al., GRL, 35 L16815 (2008).
    if (fD_um < 1.2) then
       DF = -2.96088*cos(2.84164*fD_um**0.6)-2.53628
    else
       DF = 0.47512/fD_um**0.63223
    end if
    
  end function fu_depletion_factor_bromine


  













  !**********************************************************************************
  !**********************************************************************************
  !**********************************************************************************
  !   JUST IN CASE: Below are the sea_salt_flux4mode_hybrid/spmume subroutines
  !                 before adding of the halogens
  !                 Note: These old routines use the volume flux and not the mass
  !                       flux as the updated routines above.  
  !**********************************************************************************
  !**********************************************************************************
  !**********************************************************************************

  
  subroutine sea_salt_flux4mode_hybrid_ORIG(aerMode, &        ! definition of the spectrum band
                                     & fWaterTempr, &    !  [K]
                                     & fWaterSalinity, & ! [ratio]
                                     & fNbrFlux, fVolFlux, fMassMeanDiam)
    !
    ! Computes the number emission flux of the sea salt for the pre-defined
    ! size spectrum using the combined Monahan-Martensson method. Both water
    ! salinity and temperature are taken into account
    !
    implicit none

    ! Imported parameters
    type(Taerosol_mode), intent(inout) :: aerMode
    real, intent(in) :: fWaterTempr, fWaterSalinity
    real, intent(out) :: fNbrFlux, fVolFlux, fMassMeanDiam

    ! Local variables
    real :: fMinD, fMaxD, fD, fD_um, fFluxDens, fFluxDens_prev, fIntegrStep, fIntegrStep_um, fTmp
    type(silam_sp) :: sp

    real, parameter :: fAbsoluteMin = 3.0e-9  ! 3 nm
    real, parameter :: fAbsoluteMax = 1.0e-4  ! 100 um

    !
    ! Algorithm: scan each size class and make an integration from min to max diameter of the
    ! number flux. Mean diameter can be then computed too and set to the aerosol back.
    !
!    sp%sp=> fu_work_string()

    !
    ! Get the parameters of the mode
    !
    fMinD = fu_min_d(aerMode)
    fMaxD = fu_max_d(aerMode)

!if(if_25_Degrees_Reference)then
!call msg('========> Reference temperature is 25 degrees')
!else
!call msg('========> Reference temperature is 15 degrees')
!endif

    !
    ! Now - careful. We have two ranges: Martensson and Monahan. Their border is 1 um. Both are
    ! valid within their own limits. Check all
    !
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

      !
      ! Get the integration step
      !
      fIntegrStep = 0.001*(fMaxD - fMinD)  ! At least 1000 steps over the range

      fD = fMinD
      fNbrFlux = 0.
      fVolFlux = 0.
      fMassMeanDiam = 0.

      if(if_25_Degrees_Reference)then
        fFluxDens = fu_shape_function(fD) * &
                  & fu_water_tempr_corr_to_25deg(fWaterTempr, fD) * &
                  & fu_salinity_corr(fWaterSalinity, fD)
      elseif(if_15_Degrees_Reference)then
        fFluxDens = fu_shape_function(fD) * &
                  & fu_water_tempr_corr_to_15deg(fWaterTempr, fD) * &
                  & fu_salinity_corr(fWaterSalinity, fD)
      else
        fFluxDens = fu_shape_function(fD) * &
                  & fu_water_tempr_corr_to_20deg(fWaterTempr, fD) * &
                  & fu_salinity_corr(fWaterSalinity, fD)
      endif

      do while(fD < fMaxD)
        fIntegrStep = min(fIntegrStep, fMaxD - fD)
        fD = fD + fIntegrStep
        fFluxDens_prev = fFluxDens

        fD_um = fD * 1e6
        fIntegrStep_um = fIntegrStep * 1e6
        !
        ! Number-flux density
        !
        if(if_25_Degrees_Reference)then
          fFluxDens = fu_shape_function(fD) * &
                    & fu_water_tempr_corr_to_25deg(fWaterTempr, fD) * &
                    & fu_salinity_corr(fWaterSalinity, fD)
        elseif(if_15_Degrees_Reference)then
          fFluxDens = fu_shape_function(fD) * &
                    & fu_water_tempr_corr_to_15deg(fWaterTempr, fD) * &
                    & fu_salinity_corr(fWaterSalinity, fD)
        else
          fFluxDens = fu_shape_function(fD) * &
                    & fu_water_tempr_corr_to_20deg(fWaterTempr, fD) * &
                    & fu_salinity_corr(fWaterSalinity, fD)
        endif
        !
        ! Integrate the number flux. Note that the shape function is per um.
        !
        fNbrFlux = fNbrFlux + fIntegrStep_um * 0.5 * (fFluxDens_prev + fFluxDens)  ! for number flux
        !
        ! And integrate the mass flux
        !
        fTmp = fIntegrStep_um * 0.5 * 0.5235987756 * (fFluxDens_prev * (fD_um-fIntegrStep_um)**3 + &
                                                    & fFluxDens * fD_um**3)
!sp%sp=> fu_work_string()
!if (error) return
!write(unit = sp%sp,fmt=*)'Dp_um, Tw, Sw, Tcorr, Scorr, f1, f2, f3, flux for D, m-2 sec-1 um-1:', &
!                       & fD_um, fMeanD, &
!                       & fWaterTempr, & ! [K]
!                       & fWaterSalinity, &
!                       & fWTcorr, fScorr, f1, f2, f3, fMonahanFunc
!call msg(sp%sp)
!call free_work_array(sp%sp)

        fVolFlux = fVolFlux + fTmp 
        fMassMeanDiam = fMassMeanDiam + fTmp * (fD_um-0.5*fIntegrStep_um)
      end do  ! cycle over diameter

      fMassMeanDiam = 1.e-6 * fMassMeanDiam / fVolFlux

      fVolFlux = fVolFlux * 1.0e-18  ! from um^3 to m^3

    endif

  end subroutine sea_salt_flux4mode_hybrid_ORIG


  !**********************************************************************************

  subroutine sea_salt_flux4mode_spume_ORIG(aerMode, &        ! definition of the spectrum band
                                    & fWaterSalinity, &
                                    & fNbrFlux, fVolFlux, fMassMeanDiam)
    !
    ! Computes the number emission flux of the sea salt for the pre-defined
    ! size spectrum using the combined Monahan-Martensson method. Both water
    ! salinity and temperature are taken into account
    !
    implicit none

    ! Imported parameters
    type(Taerosol_mode), intent(inout) :: aerMode
    real, intent(in) :: fWaterSalinity
    real, intent(out) :: fNbrFlux, fVolFlux, fMassMeanDiam

    ! Local variables
    real :: fMinD, fMaxD, fD, fD_um, fFluxDens, fFluxDens_prev, fIntegrStep, fIntegrStep_um, fTmp
    type(silam_sp) :: sp

    real, parameter :: fAbsoluteMin = 1.0e-6  ! 1 um
    real, parameter :: fAbsoluteMax = 1.0e-4  ! 100 um

    !
    ! Algorithm: scan the size range and make an integration from min to max diameter of the
    ! number flux. Mean diameter can be then computed too and set to the aerosol back.
    ! Get the parameters of the mode
    !
    fMinD = fu_min_d(aerMode)
    fMaxD = fu_max_d(aerMode)

    !
    ! Now - careful. Spume particles strart from 10 um
    !
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

      !
      ! Get the integration step
      !
      fIntegrStep = 0.001*(fMaxD - fMinD)  ! At least 1000 steps over the range

      fD = fMinD
      fNbrFlux = 0.
      fVolFlux = 0.
      fMassMeanDiam = 0.

      fFluxDens = fu_shape_function_spume(fD) * &
                & fu_salinity_corr(fWaterSalinity, fD)

!sp%sp=> fu_work_string()
!if (error) return
!write(unit = sp%sp,fmt=*)'Dp_um, Tw, Sw, Tcorr, Scorr, f1, f2, f3, flux for D, m-2 sec-1 um-1:', &
!                       & fD_um, fMeanD, &
!                       & fWaterTempr, & ! [K]
!                       & fWaterSalinity, &
!                       & fWTcorr, fScorr, f1, f2, f3, fMonahanFunc
!call msg(sp%sp)
!call free_work_array(sp%sp)

      do while(fD < fMaxD)
        fIntegrStep = min(fIntegrStep, fMaxD - fD)
        fD = fD + fIntegrStep
        fFluxDens_prev = fFluxDens

        fD_um = fD * 1e6
        fIntegrStep_um = fIntegrStep * 1e6

        fFluxDens = fu_shape_function_spume(fD) * fu_salinity_corr(fWaterSalinity, fD)

        !
        ! Integrate the volume flux. Note that the shape function is per um.
        !
        fTmp = fIntegrStep_um * 0.5 * (fFluxDens_prev + fFluxDens)  ! for number flux
        fNbrFlux = fNbrFlux + fTmp
        !
        ! And integrate the number flux
        !
        fTmp = fIntegrStep_um * 0.5 * 0.5235987756 * (fFluxDens_prev * (fD_um-fIntegrStep_um)**3 + &
                                                    & fFluxDens * fD_um**3)
!call msg('Diameter and volume flux',fD_um,fVolFlux)
        fVolFlux = fVolFlux + fTmp 
        fMassMeanDiam = fMassMeanDiam + fTmp * (fD_um-0.5*fIntegrStep_um)
      end do  ! cycle over diameter

      fMassMeanDiam = 1.e-6 * fMassMeanDiam / fVolFlux

      fVolFlux = fVolFlux * 1.0e-18  ! from um^3 to m^3

    endif

  end subroutine sea_salt_flux4mode_spume_ORIG

  
END MODULE source_terms_sea_salt


