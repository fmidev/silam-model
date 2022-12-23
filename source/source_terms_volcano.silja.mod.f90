MODULE source_terms_volcano

  ! This module contains the general description of volcano-type source-term of SILAM
  ! It describes the temporal distribution of the release to atmosphere: the amount of 
  ! SO2 and ash materials released.
  !
  ! Co-ordinates of the source (position) are always stored in the geographical 
  ! grid as in ini file. 
  !
  ! All units: SI
  ! 
  ! All degrees: real numbers in degrees (degrees and parts per
  ! hundred/parts per thousand , NOT in degrees and minutes)
  ! 
  ! Language: ANSI standard Fortran 90
  !
  ! Code author: Mikhail Sofiev, FMI
  ! 
  ! Modules used:

  USE source_terms_time_params   ! derive from this one

  
  IMPLICIT NONE

  private
  
  ! The public functions and subroutines available in this module:
  
  public reserve_volcano_source
  public add_input_needs_volcano_source
  public fill_volc_src_from_namelist
  PUBLIC report
  public defined
  PUBLIC fu_name
  public fu_sector
  PUBLIC fu_start_time
  PUBLIC fu_end_time
  PUBLIC fu_duration
  public total_amt_species_unit
  public fu_source_id_nbr
  public fu_source_nbr
  public link_volc_src_to_species
  public add_source_species_volc_src
  public source_2_map_volc_src
  public project_volc_src_2_grids
  public create_volc_src_cont_grd
  public fu_volc_emis_owned_quantity
  public prepare_volc_src_vert_params
  public inject_emission_euler_volc_src
  public inject_emission_lagr_volc_src
  public store_volc_src_as_namelist

  ! Data assimilation subset
  public assimilation_request_volc_src
  public observe_params_volc_src
  public inject_params_volc_src

  ! Private routines
  
  private total_from_v_src_descr_unit
  private total_from_v_src_species_unit
  private determine_release_params
  private Mastin_height2emis
  private fu_name_v_src
  private fu_sector_v_src
  private fu_source_nbr_of_v_src
  private fu_source_id_nbr_of_v_src
  private fu_NbrOfTimeSlots_of_v_src
  private fu_volcano_source_start_time
  private fu_volcano_source_end_time
  private fu_volcano_source_duration
  private fu_cocktail_descr_of_v_src
  private getTimeSlots_of_volcano_source
  private fu_SlotTime_of_volcano_source
  private report_volcano_source
  private fu_volcano_source_defined

  
  ! Generic names and operator-interfaces of some functions:

  INTERFACE report
    module procedure report_volcano_source
  END INTERFACE

  interface defined
    module procedure fu_volcano_source_defined
  end interface

  INTERFACE fu_name
    module procedure fu_name_v_src
  END INTERFACE

  INTERFACE fu_sector
    module procedure fu_sector_v_src
  END INTERFACE

  INTERFACE fu_start_time
    MODULE PROCEDURE fu_volcano_source_start_time
  END INTERFACE

  INTERFACE fu_end_time
    MODULE PROCEDURE fu_volcano_source_end_time
  END INTERFACE

  INTERFACE fu_duration
    MODULE PROCEDURE fu_volcano_source_duration
  END INTERFACE

  interface total_amt_species_unit
    module procedure total_from_v_src_species_unit
  end interface

  interface fu_source_id_nbr
    module procedure fu_source_id_nbr_of_v_src
  end interface

  interface fu_source_nbr
    module procedure fu_source_nbr_of_v_src
  end interface

  ! Possible values of emissin_vs_vertical: 
  ! For volcano source
  integer, private, parameter :: MastinRelation = 7781
  integer, private, parameter :: fullEmissionMatrix = 7782
  !
  ! For injecting Lagrangian particles
  type(field_4d_data_ptr), private, pointer, save :: fldPressure, fldRdown, fldZ

  !
  ! ATTENTION!! Each source is uniquely identified
  !             src_nm + sector_nm or src_nbr. src_nbr is absolutely unique, while src_nm and
  !             sector_nm each may show up several times (but their
  !             combination is unique). This trick allows the following:
  !          1. Complete separation of all source-sector pairs
  !          2. Split of sources, with sectors summed-up
  !          3. Split of sectors, with sources summed-up
  !             This informaion is stored in iSrcIdType field, one for all sources.

  !===================================================================================
  !
  ! Volcano source is a point source with varying strength, composition and injection height
  ! It can also involve some empirical models for connecting injection height and amount
  ! Some of its parameters can be assimilated
  !
  ! Assimilation map type describes ONE PARAMETER in the assimilation run.
  !
  type TassimMap_volvano
    private
    character(len=clen) :: chAssimNm  ! what parameter to assimilate - as in ini file
    integer :: assimType              ! what do we assimilate
    integer :: nVals, indMain         ! number of values to assimilate / dimention of covariance matrix, index in the main array
    integer, dimension(:), allocatable :: conditions  ! restrictions to the parameter
    real, dimension(:), allocatable :: ini            ! initial values for assimilated scalar parameters
    real, dimension(:,:), allocatable :: ems_zt_ini  ! initial values for assimilated 3d emission
    real, dimension(:,:), allocatable :: cov          ! covariances for all parameters: fixed in time
  end type TassimMap_volvano
  private TassimMap_volvano
  
  integer, private, parameter :: mastin_param_assimilation = 30901
  integer, private, parameter :: emis_z_t_assimilation = 30903
  
  !
  ! Main volcano source
  !
  type silam_volcano_source
    private
    ! General metadata
    CHARACTER(len=clen) :: src_nm, sector_nm ! Name of the point source and sector
    integer :: emission_vs_vertical   ! MASTIN relation or EXPLICIT_PROFILE in par_str
    real :: lon, lat, fXDispGrd, fYDispGrd
    integer :: nDescriptors, src_nbr, id_nbr     ! source and id numbers in a WHOLE source list
    integer :: ixDispGrd, iyDispGrd, izStackDispVert, nzDispVert
    logical :: ifUseTimeVarCoef, ifRateVary, ifEmissionLagrangian, if_inside_domain
    type(Tcocktail_descr), dimension(:), pointer :: cocktail_descr_lst   ! for the whole source
    type(silam_vertical) :: vert4emis, vert4emis_dispVert
    integer, dimension(:), pointer :: nSpeciesInDescr
    real, dimension(:,:), pointer :: fDescr2SpeciesUnit
    integer, dimension(:,:), pointer :: pEmisSpeciesMapping
    real, dimension(:), allocatable :: levFraction, levFractDispVert, fzDisp, dz_m
    ! Mastin model: total emission vs injection height : height as a function of DRE volume flux m3/s H=a*V^b
    real, dimension(:), pointer :: fMastin_cocktail_fraction  ! for each cocktail, (nDescriptors)
    real, pointer :: fMastin_height_pwr, fMastin_height_scaling ! power and scale
    real, pointer :: fMastin_fract_mass_hat, fMastin_fract_height_hat ! fractions of mass and height of mushroom
    ! Data storage
    real, dimension(:,:,:), allocatable :: ems3d, ems3d_dispVert   ! (nDescriptors, nz, nTimes) general
    character(len=clen) :: chEmisUnit                  ! if general matrix is given
    type(silam_release_parameters), dimension(:), allocatable :: params
    ! parameter assimilation
    type(TassimMap_volvano), dimension(:), allocatable :: assimMaps
    type(silja_logical) :: defined
  end type silam_volcano_source
  type volc_src_ptr
    TYPE(silam_volcano_source) :: volc_src
  end type volc_src_ptr
  public volc_src_ptr
  

  CONTAINS


  !*************************************************************************

  subroutine add_input_needs_volcano_source(v_src, q_met_dynamic, q_met_static, &
                                          & q_disp_dynamic, q_disp_static)
    !
    ! Returns input needs for the volcano emission computations 
    !
    implicit none

    ! Imported parameters
    type(silam_volcano_source), intent(in) :: v_src
    integer, dimension(:), intent(inout) :: q_met_dynamic, q_met_static, &
                                            & q_disp_dynamic, q_disp_static
    ! Local variables
    integer :: iTmp
    !
    ! The only need now is for Lagrangian pre-distribution of particles
    !
    if(v_src%ifEmissionLagrangian)then
      iTmp = fu_merge_integer_to_array(abl_height_m_flag, q_met_dynamic)
      iTmp = fu_merge_integer_to_array(height_flag, q_met_dynamic)
      iTmp = fu_merge_integer_to_array(R_down_meteo_flag, q_met_dynamic)
      iTmp = fu_merge_integer_to_array(pressure_flag, q_met_dynamic)
    endif

  end subroutine add_input_needs_volcano_source


  !**************************************************************************

  subroutine reserve_volcano_source(v_src, &        ! Src to initialise
                                  & iSrcNbr, &      ! Src number in the grand list
                                  & iSrcIdNbr, &    ! SrcID number
                                  & nDescriptors, & ! number of chemical descriptors to be reserved
                                  & iDynamicEnvironment)
    !
    ! Initialises the source:
    ! - stores the reference information: source number and source ID number
    ! - stores the total number of chemical descriptors that will be stored in the source
    ! - nullifies the source dynamic arrays 
    ! - set a few basic internal variables of the source
    !
    implicit none

    ! Imported parameters
    type(silam_volcano_source), intent(inout) :: v_src
    integer, intent(in) :: iSrcNbr, iSrcIdNbr, nDescriptors, iDynamicEnvironment
    !
    ! Nullify the basic variables
    !
    v_src%src_nm = ''
    v_src%sector_nm = ''
    allocate(v_src%cocktail_descr_lst(nDescriptors))
    !
    ! Main source parameters - enough to identify it in the global information list
    !
    v_src%src_nbr = iSrcNbr
    v_src%id_nbr = iSrcIdNbr
    v_src%nDescriptors = nDescriptors
    !
    ! Emission goes into Eulerian or Lagrangian clouds. If Lagrangian is available, take it!
    !
    v_src%ifEmissionLagrangian = fu_if_lagrangian_present(iDynamicEnvironment)
    !
    ! Finally, mark the source as incomplete
    !
    v_src%defined = silja_false

  end subroutine reserve_volcano_source


  ! ***************************************************************

  subroutine fill_volc_src_from_namelist(nlSrc, v_src)
    !
    ! Reads and sets one volcano source term from a namelist.
    !
    ! All units: SI, unless otherwise is stated
    !
    IMPLICIT NONE
    !
    ! Imported parameters
    TYPE(silam_volcano_source), intent(inout), target :: v_src
    type(Tsilam_namelist), pointer:: nlSrc

    ! Local variables
    type(silam_sp) :: spContent
    real, dimension(:), pointer :: arTmp
    integer :: nItems, nFields, iTmp, nz, nt, iStat, iTime, iz, iC, nItemsZT
    logical :: ifIjectTop, ifFullEmisMatrix
    logical, dimension(:), allocatable :: ifSpecies
    character(len=clen) :: chEmisUnit, chCocktail
    type(Tsilam_nl_item_ptr), dimension(:), pointer :: ptrItems
    character(len=clen), dimension(100) :: arStr
    type(silja_logical) :: ifMastinUpdated 
    character(len=*), parameter :: sub_name = "fill_volc_src_from_namelist"
    !
    ! Source name, sector, position
    !
    v_src%src_nm = fu_Content(nlSrc, 'source_name')
    v_src%sector_nm = fu_Content(nlSrc, 'source_sector_name')
    v_src%lon = fu_content_real(nlSrc,'source_longitude')
    v_src%lat = fu_content_real(nlSrc,'source_latitude')
    v_src%chEmisUnit = fu_content(nlSrc,'emission_unit')
    !
    ! par_str_volcano: temporal evolution
    !
    ! order of fields in the line: counter y4 m2 d2 h2 m2 sec(real) [<cocktail> <top_height>  ...]
    ! each time stamp shows the height at that moment. In-between, height is interpolated.
    ! The first line starts eruption, the last one ends it.
    ! Note that height is used only if emissin_vs_vertical = MASITN but time stamps and cocktails 
    ! are universal
    !
    call get_items(nlSrc, 'par_str_volcano', ptrItems, nt)
    if(error .or. nt < 1)then
      call report(nlSrc)
      call set_error('No par_str_volcano items in namelist', sub_name )
      return
    endif
    if(.not. allocated(v_src%params))then
      allocate(v_src%params(nt), stat=iTmp)
      if(fu_fails(iTmp == 0,'Failed to allocate parameters for volcano source',  sub_name ))return
      call set_missing(v_src%params, nt)
    endif
    do iTmp = 1, nt
      call decode_time_params_volc_v1(fu_content(ptrItems(iTmp)), fu_content(nlSrc,'vertical_unit'), &
                                    & v_src%params, &
                                    & v_src%cocktail_descr_lst, &    ! list of source cocktail descriptors
                                    & v_src%nDescriptors, &
                                    & ifSpecies)
      if(error)return
    end do
    !
    ! Set the source vertical: whatever the emission definition, the vertical will be the same
    !
    call set_vertical(nlSrc, v_src%vert4emis) ! e.g.L vertical_unit = km   plus   layer_thickness = 1. 1. 2. 2. 2.
                                             ! or a multitude of lines  layer_borders = val_bottom val_top
    if(error)then
      call report(nlSrc)
      call set_error('Failed vertical setup', sub_name )
      return
    endif
    nz = fu_NbrOfLevels(v_src%vert4emis)
    !
    ! The main storage of volcano emission
    !
    allocate(v_src%ems3d(v_src%nDescriptors, nz, nt), stat = iStat)
    if(iStat /= 0)then
      call report_volcano_source(v_src)
      call set_error('Failed allocatoin of z-t emission matrix', sub_name )
      return
    endif
    v_src%ems3d = 0.0
    !
    ! Volcano eruption can be given explicitly or via Mastin equation linking height and mass
    !
    spContent%sp => fu_work_string()
    if(error)return
    select case(fu_str_u_case(fu_Content(nlSrc,'emission_vs_vertical')))
      
      case('MASTIN')                                   ! MASTIN relation
        v_src%emission_vs_vertical = MastinRelation
        v_src%chEmisUnit = 'native'
        !
        ! For MASTIN relation, cloud features are needed (top height is in params%top)
        !
        allocate(v_src%fMastin_cocktail_fraction(v_src%nDescriptors), stat = iStat)
        if(iStat /= 0)then
          call report_volcano_source(v_src)
          call set_error('Failed allocatoin of Mastin cocktail fractions', sub_name )
          return
        endif
        call set_scalar_as_val(v_src%fMastin_height_pwr, fu_content_real(nlSrc,'mastin_height_power'))   ! = 0.241
        call set_scalar_as_val(v_src%fMastin_height_scaling, fu_set_named_value(fu_content(nlSrc,'mastin_height_scaling')))  ! = 2.0 km
        call set_scalar_as_val(v_src%fMastin_fract_mass_hat, fu_content_real(nlSrc,'mastin_mass_fraction_hat'))  ! = 0.75
        call set_scalar_as_val(v_src%fMastin_fract_height_hat, fu_content_real(nlSrc,'mastin_height_fraction_hat'))  ! = 0.25
        spContent%sp = fu_content(nlSrc,'mastin_cocktail_fraction')
        allocate(v_src%fMastin_cocktail_fraction(v_src%nDescriptors), stat=iStat)
        if(fu_fails(iStat==0,'Failed cocktail fractions allocation', sub_name ))return
        read(unit=spContent%sp, fmt=*, iostat=istat) &
                   & (chCocktail, v_src%fMastin_cocktail_fraction(fu_index(chCocktail, v_src%cocktail_descr_lst)), &
                    & iTmp=1, v_src%nDescriptors)
        if(istat /= 0)then
          call report_volcano_source(v_src)
          call set_error('Failed reading mastin_cocktail_fraction', sub_name )
          return
        endif
        
      case('EXPLICIT_PROFILE')                         ! Full profile
        
        v_src%emission_vs_vertical = fullEmissionMatrix 
        !
        ! Set the martrix of the emission rates for each cocktail
        ! The emission line is < time_nbr cocktail emis_z1 emis_z2 ... emis_zn >
        !
        call get_items(nlSrc, 'emission_cocktail_z_t', ptrItems, nItems)
        if(error .or. nItems < 1)then
          call report(nlSrc)
          call set_error('Wrong z-t emission lines', sub_name )
          return
        endif
        !
        ! We shall scan all available lines, search for the specific cocntail and set 
        ! the corresponding items. Note the line format:
        !     emission_cocktail_z_t = 1 COCKTAIL_ASH 0 0 1.0 0 0
        !
        do iTmp = 1, nItems
          spContent%sp = fu_content(ptrItems(iTmp))
          read(unit=spContent%sp,fmt=*, iostat=iStat) iTime, chCocktail, &
                            & (v_src%ems3d(fu_index(chCocktail, v_src%cocktail_descr_lst), iz, iTime), &
                             & iz = 1, fu_NbrOfLevels(v_src%vert4emis)) 
          if(iStat /= 0)then
            call report(nlSrc)
            call set_error('Failed reading emission_cocktail_z_t:' + fu_content(ptrItems(iTmp)), sub_name )
            return
          endif
        end do  ! emission_cocktail_z_t  items

      case default
        call set_error('emission_vs_vertical must be MASTIN or EXPLICIT_PROFILE, not:' + &
                     & fu_Content(nlSrc,'emission_vs_vertical'), sub_name )
        return
    end select   ! if detailed emission matrix is given or Mastin relation is used
    !
    !----------------------------------------------------------------
    !
    ! Assimilation section. Assimilation is considered defined if the standard deviation or covariance
    ! of the parameter is given. Each assimilation request has its name.
    ! Note that Mastin-parameters and z-t emission assimilation are mutually exclusive. 
    !
    ifMastinUpdated = silja_undefined  ! What are we going to assimilate?
    !
    ! First, decide, how many such parameters are given
    !
    call get_items(nlSrc, 'assimilated_param_covariance', ptrItems, nItemsZT)  ! covariances for z-t matrices
    if(nItemsZT > 0) nItemsZT = v_src%nDescriptors                             ! max number of z-t matrices
    call get_items(nlSrc, 'assimilated_param_variance', ptrItems, nItems)   ! Individual parameters
    if(error)return
    
    if(nItems + nItemsZT > 0)then
      allocate(v_src%assimMaps(nItems + nItemsZT), stat=iStat)
      if(fu_fails(iStat==0,'Failed allocation of assimilation mapping', sub_name ))return
      do iTmp = 1, nItems+nItemsZT
        v_src%assimMaps(iTmp)%assimType = int_missing
      end do
      !
      ! Decode the assimilation requests for individual parameters
      ! Note that the assimilated parameter has to be connected with the main one in order to simplify
      ! the "inject" operation.
      ! Solution: The assimMap structure is made the only actual parameter holder. The "named"
      ! parameters are just pointers looking at the specific elements of this structure.
      !
      do iTmp = 1, nitems
        spContent%sp = fu_content(ptrItems(iTmp))
        read(unit=spContent%sp,fmt=*, iostat=iStat) v_src%assimMaps(iTmp)%chAssimNm
        if(fu_fails(iStat==0,'Failed reading assimilation request:'+spContent%sp, sub_name ))return
        spContent%sp = adjustl(spContent%sp(index(spContent%sp,' ')+1 : ))  ! cut out the name
        !
        ! A special case: one line has several cocktails and variances:
        ! assimilated_param_variance = mastin_cocktail_fraction COCKTAIL_SO2 0.1. COCKTAIL_ASH 1.
        !
        if(v_src%assimMaps(iTmp)%chAssimNm == 'mastin_cocktail_fraction')then
          call check_assimilation_setup(ifMastinUpdated, silja_true)
          call split_string(spContent%sp, ' ', arStr, nz)
          if(error)return
          v_src%assimMaps(iTmp)%assimType = mastin_param_assimilation
          allocate(v_src%assimMaps(iTmp)%ini(v_src%nDescriptors), v_src%assimMaps(iTmp)%cov(v_src%nDescriptors, 1), &
                 & v_src%assimMaps(iTmp)%conditions(1), stat=iStat)
          if(fu_fails(iStat==0,'Failed allocation of cocktail fraction for assimilation','set_signle_assim_param'))return
          v_src%assimMaps(iTmp)%conditions(1) = non_negative

          v_src%assimMaps(iTmp)%ini(1:v_src%nDescriptors) = v_src%fMastin_cocktail_fraction(1:v_src%nDescriptors)
          deallocate(v_src%fMastin_cocktail_fraction)
          v_src%fMastin_cocktail_fraction => v_src%assimMaps(iTmp)%ini
          
          do iz = 1, v_src%nDescriptors * 2, 2  ! go through the descriptor-stdev couples one by one
            iC = fu_index(arStr(iz), v_src%cocktail_descr_lst)
            if(fu_fails(iC>0 .and. iC<=v_src%nDescriptors,'strange cocktail name in the line:'+ spContent%sp, sub_name ))return
            read(unit=arStr(iz+1), fmt=*, iostat=iStat) v_src%assimMaps(iTmp)%cov(iC,1)
            if(fu_fails(iStat==0,'Failed reading the cocktail variance line:'+spContent%sp, sub_name ))return
          end do

        elseif(v_src%assimMaps(iTmp)%chAssimNm == 'mastin_height_power')then
          call check_assimilation_setup(ifMastinUpdated, silja_true)
          call set_signle_assim_param(spContent%sp, v_src%assimMaps(iTmp), v_src%fMastin_height_pwr)
          
        elseif(v_src%assimMaps(iTmp)%chAssimNm == 'mastin_height_scaling')then
          call check_assimilation_setup(ifMastinUpdated, silja_true)
          call set_signle_assim_param(spContent%sp, v_src%assimMaps(iTmp), v_src%fMastin_height_scaling)

        elseif(v_src%assimMaps(iTmp)%chAssimNm == 'mastin_mass_fraction_hat')then
          call check_assimilation_setup(ifMastinUpdated, silja_true)
          call set_signle_assim_param(spContent%sp, v_src%assimMaps(iTmp), v_src%fMastin_fract_mass_hat)
        
        elseif(v_src%assimMaps(iTmp)%chAssimNm == 'mastin_height_fraction_hat')then
          call check_assimilation_setup(ifMastinUpdated, silja_true)
          call set_signle_assim_param(spContent%sp, v_src%assimMaps(iTmp), v_src%fMastin_fract_height_hat)
        
        else
          call set_error('Unknown assimilation request:' + spContent%sp, sub_name )
        endif  ! simple or multi-parameter line
        if(error)return
        call msg('Single parameter to assimilate: ' + v_src%assimMaps(iTmp)%chAssimNm)
      end do  ! assimilation standard deviations
      !
      ! The z-t emission 2D matrix can also be assimilated but it has non-zero spatial covariance
      ! And it cn be made only for some cocktails
      !
      call get_items(nlSrc, 'assimilated_param_covariance', ptrItems, nItemsZT)
      if(nItemsZT > 0)then
        !
        ! Scan the lines with covariances keeping in mind that each covariance matrix is described 
        ! by several lines
        !
        do iTmp = 1, nItemsZT
          spContent%sp = fu_content(ptrItems(iTmp))
          !
          ! Actually, only one covariance-involving parameter is defined for volcanoe source
          ! Line format: emission_cocktail_z_t <cocktail_name> <line number> <val1> ... <valn>
          ! Note: this is time-dependent vector with correlated components, same covariance for all times
          !
          if(fu_fails(index(fu_str_l_case(spContent%sp),'emission_cocktail_z_t') == 1, &
                      & 'Unknown assimilation request:'+spContent%sp, sub_name ))return

          call check_assimilation_setup(ifMastinUpdated, silja_false)

          spContent%sp = adjustl(spContent%sp(index(spContent%sp,' '):))  ! cut out assim type
          !
          ! Get the cocktail name, which covariance is in this line and allocate the structures
          !
          iC = fu_index(spContent%sp(1:index(spContent%sp,' ')), v_src%cocktail_descr_lst)  ! get the cocktail number
          nz = fu_nbrOfLevels(v_src%vert4emis)
          nt = size(v_src%params)
          if(.not.allocated(v_src%assimMaps(iC+nItems)%cov))then
            allocate(v_src%assimMaps(iC+nItems)%cov(nz, nz), &
                   & v_src%assimMaps(iC+nItems)%ems_zt_ini(nz, nt), &
                   & v_src%assimMaps(iC+nItems)%conditions(1), stat=iStat)
            v_src%assimMaps(iC+nItems)%cov(:,:) = 0
            v_src%assimMaps(iC+nItems)%ems_zt_ini(:,:) = 0
            v_src%assimMaps(iC+nItems)%assimType = emis_z_t_assimilation
          endif
          read(unit=spContent%sp,fmt=*,iostat=iStat) arStr(1), iz, v_src%assimMaps(iC+nItems)%cov(1:nz, iz)
          if(fu_fails(iStat==0,'Failed reading covariance:'+spContent%sp, sub_name ))return
          v_src%assimMaps(iC+nItems)%ems_zt_ini(:,:) = v_src%ems3d(iC,:,:)  ! store the initial values
          v_src%assimMaps(iC+nItems)%nVals = nz
          v_src%assimMaps(iC+nItems)%conditions(1) = non_negative
        enddo  ! covariane-containing items
      else
        call msg('No z-t emission matrix assimilation')
      endif  ! z-t covariance given
    else
      call msg('No assimilation for volcano source')
    endif  ! if standard deviation is given: assimilation is requested

    call free_work_array(spContent%sp)  
      
    v_src%defined = silja_true

    CONTAINS
  
      !========================================================================

      subroutine set_scalar_as_val(pScalar, fVal)
        !
        ! Allocates the scalar pointer and assigns its value
        !
        implicit none
        ! Imported parameters
        real, pointer :: pScalar
        real, intent(in) :: fVal
        ! local variables
        integer :: iStat
        
        allocate(pScalar, stat=iStat)
        if(fu_fails(iStat==0,'Failed parameter allocation','set_scalar_as_val'))return
        pScalar = fVal
        
      end subroutine set_scalar_as_val
      
      !=========================================================================

      subroutine set_signle_assim_param(string, assimMap, pNamedParam)
        !
        ! sets a single assimilation parameter. It also copies the value of the named parameter
        ! to the structure and reallocates the pointer to look at that structure.
        !
        implicit none
        character(len=*), intent(in) :: string
        type(TassimMap_volvano), intent(inout), target :: assimMap
        real, pointer :: pNamedParam
        ! Locals
        integer :: iStat
        
        ! prepare place and set basics
        assimMap%nVals = 1
        assimMap%assimType = mastin_param_assimilation
        allocate(assimMap%ini(1), assimMap%cov(1,1), assimMap%conditions(1), stat=iStat)
        if(fu_fails(iStat==0,'Failed allocation of assimilation map values','set_signle_assim_param'))return
        !
        ! swap and redirect named parameter
        !
        assimMap%ini(1) = pNamedParam
        deallocate(pNamedParam)
        pNamedParam => assimMap%ini(1)
        !
        ! read covariance: either a single number or a number with a unit
        !
        if(index(trim(string),' ') > 0)then
          assimMap%cov(1,1) = fu_set_named_value(string)
        else
          read(unit=string, fmt=*, iostat=iStat) assimMap%cov(1,1)
        endif
        assimMap%conditions(1) = non_negative
      end subroutine set_signle_assim_param

      !==========================================================
      
      subroutine check_assimilation_setup(ifMastinUpdated, ifNeedMastinUpdate)
        !
        ! If Mastin parameters are assimilated, we would have to recompute the z-t emission matrix
        ! each time the parameters are updated. Evidently, this matrix must not be assimilated in 
        ! this case. This function enforces the exclusivity
        !
        implicit none
        
        ! Imported parameters
        type(silja_logical), intent(inout) :: ifMastinUpdated  ! status of assimilation requests this-far
        type(silja_logical), intent(in) :: ifNeedMastinUpdate    ! new request - what do you want?

        if(fu_true(ifNeedMastinUpdate))then   ! new request assimilates Mastin
          if(fu_false(ifMastinUpdated))then
            call set_error('Both Mastin parameters and z-t emission matrix cannot be assimilated','check_assimilation_setup')
            return
          else
            ifMastinUpdated = silja_true
          endif
        elseif(fu_false(ifNeedMastinUpdate))then    ! new request is z-t matrix, no Mastin assimilation
          if(fu_true(ifMastinUpdated))then
            call set_error('Both z-t emission matrix and Mastin parameters cannot be assimilated','check_assimilation_setup')
            return
          else
            ifMastinUpdated = silja_false
          endif
        else
          ! Do nothing: new request does not bother with Mastin or z-t matrix, i.e. poses no limitations
        endif
      end subroutine check_assimilation_setup
      
  end subroutine fill_volc_src_from_namelist

  
  !*****************************************************************

  subroutine total_from_v_src_descr_unit(v_src, &                ! input
                                       & amounts, &              ! output
                                       & start, duration, layer) ! input
    !
    ! Returns the amount of the released material IN DESCRIPTOR UNIT starting from 
    ! start during the duration time interval.
    !
    ! We have to integrate the release rate, which is given in descriptor basic units - for each
    ! descriptor. Note that the release rates of volcano are CONSTANT within each time slot and given
    ! in volcano_emis matrix, NOT by params 
    !
    ! This subroutine produces a total descriptor-wise emission and returns a vector. Vertical 
    ! dimension is integrated over
    !
    ! Units: SI basic and derived ones
    !
    implicit none

    ! Imported parameters with intent IN
    type(silam_volcano_source), intent(in) :: v_src
    real, dimension(:), intent(out) :: amounts
    type(silja_time), intent(in) :: start
    type(silja_interval), intent(in) :: duration
    type(silja_level), intent(in) :: layer  ! if defined, include only fraction emitted in it

    ! Local variables
    integer :: iSlot, iStartSlot, iEndSlot, iLev, iDescr
    type(silja_time) :: start_integr, end_integr
    real :: fTime_sec

    amounts(1:v_src%nDescriptors) = 0.0
    !
    ! Time period for integration, then start and end slots
    !
    call get_overlapping_slots(v_src%params, start, duration, iStartSlot, iEndSlot)
    if(error)return
    if(iStartSlot == int_missing .or. iEndSlot == int_missing .or. iStartSlot >= iEndSlot)then
      call msg_warning('Source does not overlap with computation period','total_from_v_src_descr_unit')
      return
    endif
    !
    ! Go with time-wise integration
    !
    do iSlot = iStartSlot, iEndSlot-1
      !
      ! Compute the time overlap
      !
      if(v_src%params(iSlot)%time < start)then ! Skip part of slot if needed
        start_integr = start
      else
        start_integr = v_src%params(iSlot)%time
      endif
      if(v_src%params(iSlot+1)%time > start+duration)then ! Skip part of the slot if needed
        end_integr = start+duration
      else
        end_integr = v_src%params(iSlot+1)%time
      endif
      fTime_sec = fu_sec(end_integr - start_integr)
      !
      ! Integrate over time and height for all species
      !
      if(defined(layer))then
        do iLev = 1, fu_NbrOfLevels(v_src%vert4emis_dispVert)
          do iDescr = 1, v_src%nDescriptors
            amounts(iDescr) = amounts(iDescr) + fTime_sec * &               ! time period
                                & fu_vert_overlap_fraction(fu_level(v_src%vert4emis_dispVert,iLev), layer) * &
                                & v_src%ems3d_dispVert(iDescr, iLev, iSlot) ! descr.rate
          end do
        end do
      else
        do iDescr = 1, v_src%nDescriptors
          amounts(iDescr) = amounts(iDescr) + fTime_sec * &               ! time period
                      & sum(v_src%ems3d_dispVert(iDescr, 1:fu_NbrOfLevels(v_src%vert4emis_dispVert), iSlot)) ! rate
        end do
      endif  ! if layer is defined
    end do  ! time slots

  end subroutine total_from_v_src_descr_unit


  !*****************************************************************

  subroutine total_from_v_src_species_unit(v_src, &                      ! input
                                         & species, nSpecies, amounts, & ! output
                                         & start, duration, layer)    ! input
    !
    ! Returns the amount of the released material IN SPECIES UNIT starting from 
    ! start during the duration time interval. Species must be initialised by that moment
    !
    ! Units: SI basic and derived ones
    !
    implicit none

    ! Imported parameters with intent IN
    type(silam_volcano_source), intent(in) :: v_src
    type(silam_species), dimension(:), pointer :: species
    integer, intent(out) :: nSpecies
    real, dimension(:), intent(out) :: amounts
    type(silja_time), intent(in) :: start
    type(silja_interval), intent(in) :: duration
    type(silja_level), intent(in) :: layer  ! if defined, include only fraction emitted in it

    ! Local variables
    integer :: iDescr, iSpecies, nSpeciesDescr
    type(silam_species), dimension(:), pointer :: speciesDescr
    real, dimension(max_descriptors_in_source) :: amountsDescr
    type(chemical_adaptor) :: adaptor

    !
    ! Get the total amounts in descriptor unit
    !
    call total_from_v_src_descr_unit(v_src, amountsDescr, start, duration, layer)
    if(error)return

!do iDescr = 1, v_src%nDescriptors
!call msg('Emission for descriptor: ' + fu_name(v_src%cocktail_descr_lst(iDescr)) + ':',amountsDescr(iDescr))
!end do

    if(sum(amountsDescr(1:v_src%nDescriptors)) == 0.) return
    !
    ! Conversion into the species units is comparatively straightforward.
    !
    nullify(species)
    nSpecies = 0
    call add_source_species_volc_src(v_src, species, nSpecies)
    if(error)return
    amounts(1:nSpecies) = 0.
    !
    ! Now explore the descriptors
    !
    do iDescr = 1, v_src%nDescriptors
      call get_inventory(v_src%cocktail_descr_lst(iDescr), speciesDescr, nSpeciesDescr)
      call create_adaptor(speciesDescr, species, adaptor)
      if (error) return
      do iSpecies = 1, nSpeciesDescr
!call msg('Factor to species unit for:'+fu_name(v_src%cocktail_descr_lst(iDescr)), v_src%fDescr2SpeciesUnit(iSpecies,iDescr))
        amounts(adaptor%iSp(iSpecies)) = amounts(adaptor%iSp(iSpecies)) + &
                                & amountsDescr(iDescr) * v_src%fDescr2SpeciesUnit(iSpecies,iDescr)
!call msg('Amounts for species:' + fu_name(fu_material(species(adaptor%iSp(iSpecies)))),amounts(adaptor%iSp(iSpecies)))
      end do
    end do

  end subroutine total_from_v_src_species_unit


  !*****************************************************************

  subroutine store_volc_src_as_namelist(v_src, uOut, ifOriginalGrd, ifExplicitProfileToo)
    !
    ! Stores the volcano source in the namelist so that it can be read as namelist
    !
    IMPLICIT NONE
    
    ! Imported parameters
    TYPE(silam_volcano_source), intent(in) :: v_src
    integer, intent(in) :: uOut
    logical, intent(in) :: ifOriginalGrd, ifExplicitProfileToo

    ! Local variables
    integer :: iz, it, iDescr
    
    if(fu_fails(fu_true(v_src%defined),'Undefined volcano source','store_volc_src_as_namelist'))return

    write(uOut,fmt='(A)')'VOLCANO_SOURCE_1'
    write(uOut, fmt='(A,A)') 'source_name = ', v_src%src_nm
    write(uOut, fmt='(A,A)') 'source_sector_name = ', v_src%sector_nm
    write(uOut, fmt='(A,F9.4)') 'source_longitude = ', v_src%lon
    write(uOut, fmt='(A,F9.4)') 'source_latitude = ', v_src%lat
    write(uOut, fmt='(A)') 'vertical_unit = m'
    !
    ! Temporal evolution
    !
    do it = 1, size(v_src%params)
      write(uOut,fmt='(A,I2,1x,A,100(A,1x, E9.5))') 'par_str_volcano = ', &
                                      & it, fu_time_to_io_string(v_src%params(it)%time), &
                                      & (fu_name(v_src%cocktail_descr_lst(iDescr)), &
                                       & v_src%params(it)%rate_descr_unit(iDescr), iDescr=1, v_src%nDescriptors)
    end do  ! times in params
    
    !
    ! Emiossion vertical and intensity
    !
    if(v_src%emission_vs_vertical == MastinRelation)then
      !
      ! Generic Mastin-type mushroom
      !
      write(uOut, fmt='(A)')'emission_vs_vertical = MASTIN'
      write(uOut, fmt='(A,F7.3)')'mastin_height_power = ', v_src%fMastin_height_pwr 
      write(uOut, fmt='(A,F7.3)')'mastin_height_scaling = ', v_src%fMastin_height_scaling 
      write(uOut, fmt='(A,F5.3)')'mastin_mass_fraction_hat = ', v_src%fMastin_fract_mass_hat
      write(uOut, fmt='(A,F5.3)')'mastin_height_fraction_hat = ', v_src%fMastin_fract_height_hat
      write(uOut, fmt='(A,100(A,1x,F5.3))')'mastin_cocktail_fraction = ', &
                             & (fu_name(v_src%cocktail_descr_lst(iDescr)), &
                              & v_src%fMastin_cocktail_fraction(iDescr), iDescr=1, v_src%nDescriptors)

    elseif(v_src%emission_vs_vertical == fullEmissionMatrix .or. ifExplicitProfileToo)then
      !
      ! Explicit profile, possibly, computed from Mastin or assimilated
      ! Report the profile and vertical in which it is written: dispersion vertical can be requested
      !
      if(v_src%emission_vs_vertical == fullEmissionMatrix) &
                                    & write(uOut, fmt='(A)')'emission_vs_vertical = EXPLICIT_PROFILE'
      if(ifOriginalGrd)then
        do it = 1, size(v_src%params)
          do iDescr =1, v_src%nDescriptors
            write(uOut,fmt='(A,I2,1x,A,100(E9.3,1x))') 'emission_cocktail_z_t = ', &
                            & it, fu_name(v_src%cocktail_descr_lst(iDescr)), &
                            & (v_src%ems3d(iDescr,iz,it), iz = 1, fu_NbrOfLevels(v_src%vert4emis))
          end do   ! descriptors
        end do  ! times in z_t
        call report_as_namelist(v_src%vert4emis, uOut)
      else
        do it = 1, size(v_src%params)
          do iDescr =1, v_src%nDescriptors
            write(uOut,fmt='(A,1x,I2,1x,A,100(E9.3,1x))') 'emission_cocktail_z_t = ', &
                                    & it, fu_name(v_src%cocktail_descr_lst(iDescr)), &
                                    & (v_src%ems3d_dispVert(iDescr,iz,it), &
                                                      & iz = 1, fu_NbrOfLevels(v_src%vert4emis_dispVert))
          end do   ! descriptors
        end do  ! times in z_t
        call report_as_namelist(v_src%vert4emis_dispVert, uOut)
      endif  ! ifOriginalGrid

    else
      call set_error('Unknown source emission-vertical relation: '+ fu_str(v_src%emission_vs_vertical), 'store_volc_src_as_namelist')
      return
    end if
    
    write(uOut,fmt='(A)')'END_VOLCANO_SOURCE_1'
    write(uOut,fmt=*)
  end subroutine store_volc_src_as_namelist
  
  
  !**************************************************************************

  subroutine inject_emission_euler_volc_src(v_src, &
                                          & mapEmis, mapCoordX, mapCoordY, mapCoordZ, & ! Eulerian
                                          & met_buf, disp_buf, &
                                          & now, timestep, &
                                          & ifSpeciesMoment, &
                                          & fMassInjected, &
                                          & interpCoefMeteo2DispHoriz, ifMetHorizInterp, &
                                          & interpCoefMeteo2DispVert, ifMetVertInterp)
    !
    ! Adds the emission flux to the emission mass map and stores the position
    ! of the injected mass to mapCoord. All maps may contain some data before,
    ! so no overwriting - just adding. Note that for successfull trick one has to 
    ! sum-up momentum, not coordinates. Prior conversion to momentum has been done in
    ! source_general for the whole mass map.
    ! Heavily utilises the precomputed variables: 
    ! - reprojected cells to the dispersion grid
    ! - reprojected vertical to the dispersion vertical
    ! Interpolation goes along these precomputed variables ONLY, no redundancy here
    !
    implicit none

    TYPE(silam_volcano_source), INTENT(in) :: v_src
    type(Tmass_map), intent(inout) :: mapEmis, mapCoordX, mapCoordY, mapCoordZ
    TYPE(Tfield_buffer), pointer :: met_buf, disp_buf
    type(silja_time), intent(in) :: now
    type(silja_interval), intent(in) :: timestep
    logical, intent(in) :: ifSpeciesMoment
    real(r8k), dimension(:), intent(inout) :: fMassInjected
    type(TVertInterpStruct), intent(in) ::  interpCoefMeteo2DispVert
    type(THorizInterpStruct), intent(in) ::  interpCoefMeteo2DispHoriz
    logical, intent(in) :: ifMetVertInterp
    logical, intent(in) :: ifMetHorizInterp

    ! Local variables
    integer :: iSlot, iLev, ix, iy, iSrc, nLevs, iDescr, iSpecies, &
             & iSlotStart, iSlotEnd, ispeciesEmis, iSpeciesDescr
    real :: fWeightPastSrc, fLevFraction, fTimeScale, fCellTotal, factor, timestep_sec
    real, dimension(max_descriptors_in_source) :: fMassTimeCommon
    real, dimension(3) :: ptrCoord
    logical :: ifFound
    type(silam_species), dimension(:), pointer :: species
    real :: plumeBottom, plumeTop, overlapTop, overlapBottom

    !
    ! Have to be very careful: starting and ending of the source slots can be pretty close
    ! especially in case of adjoint run. The same is true for time variation coefficients,
    ! which can vary sharply. Therefore, the total mass injected to the grid will be
    ! integrated along time.
    !
    nullify(species)
    if(.not. v_src%if_inside_domain)return  ! or outside the whole domain if not MPI
    !
    ! General parameters of the release
    !
    call determine_release_params(v_src, now, timestep, mapEmis%nSpecies, &
                                & met_buf, interpCoefMeteo2DispHoriz, ifMetHorizInterp, &
                                & interpCoefMeteo2DispVert, ifMetVertInterp, &
                                & fMassTimeCommon, plumeTop, plumeBottom, fWeightPastSrc, ifFound, &
                                & iSlotStart, iSlotEnd)
    if(error)return
    if(.not. ifFound) return                   ! if any emission found
    !
    ! We do the injection layer by layer (dispersion vertical !)
    !
    do iLev = 1, nz_dispersion
      !
      ! First do the check for the overlap: speed-up
      !
      ptrCoord(3) = v_src%fzDisp(iLev)
      fLevFraction = v_src%levFractDispVert(iLev)

      if(fLevFraction < 1e-5)cycle  ! nothing for this dispersion layer
!   call msg('ilev, vertical fraction:', ilev, fLevFraction)
      ptrCoord(1) = v_src%fXDispGrd - v_src%ixDispGrd
      ptrCoord(2) = v_src%fYDispGrd - v_src%iyDispGrd
      fCellTotal = 0.0
      if(error)return
      !
      ! Emit the mass going species by species
      !
      do iDescr = 1, v_src%nDescriptors
        !
        ! Time variation coefficients can be handled this way:
        !
        do iSpeciesDescr = 1, v_src%nSpeciesInDescr(iDescr)
          iSpeciesEmis = v_src%pEmisSpeciesMapping(iSpeciesDescr,iDescr)
          factor = v_src%fDescr2SpeciesUnit(iSpeciesDescr,iDescr) * fLevFraction 
          mapEmis%arM(iSpeciesEmis, v_src%id_Nbr, iLev,v_src%ixDispGrd,v_src%iyDispGrd) = & 
                 & mapEmis%arM(iSpeciesEmis, v_src%id_Nbr, iLev,v_src%ixDispGrd,v_src%iyDispGrd) + &
                 & fMassTimeCommon(iDescr) * factor
          fCellTotal = fCellTotal + fMassTimeCommon(iDescr) * factor
          fMassInjected(iSpeciesEmis) = fMassInjected(iSpeciesEmis) + fMassTimeCommon(iDescr) * factor

          if (ifSpeciesMoment) then
            mapCoordX%arm(iSpeciesEmis,v_src%id_nbr,ilev,v_src%ixDispGrd,v_src%iyDispGrd) = &
               & mapCoordX%arm(iSpeciesEmis, v_src%id_nbr, ilev, v_src%ixDispGrd, v_src%iyDispGrd) + &
               & fMassTimeCommon(iDescr) * factor * ptrCoord(1)
            mapCoordY%arm(iSpeciesEmis, v_src%id_nbr, ilev, v_src%ixDispGrd, v_src%iyDispGrd) = &
               & mapCoordY%arm(iSpeciesEmis, v_src%id_nbr, ilev, v_src%ixDispGrd, v_src%iyDispGrd) + &
               & fMassTimeCommon(iDescr) * factor * ptrCoord(2)
            mapCoordZ%arm(iSpeciesEmis, v_src%id_nbr, ilev, v_src%ixDispGrd, v_src%iyDispGrd) = &
               & mapCoordZ%arm(iSpeciesEmis, v_src%id_nbr, ilev, v_src%ixDispGrd, v_src%iyDispGrd) + &
               & fMassTimeCommon(iDescr) * factor * ptrCoord(3)
          end if

        end do  ! species inside the descriptor
      end do  ! iDescr

      ! If we use bulk moment, add it of the injected masses and the total injected mass
      !
      if (.not. ifSpeciesMoment) then
        mapCoordx%arM(1,v_src%id_nbr, iLev, v_src%ixDispGrd,v_src%iyDispGrd) = &
                     & mapCoordx%arM(1,v_src%id_nbr, iLev, v_src%ixDispGrd,v_src%iyDispGrd) + &
                     & ptrCoord(1) * fCellTotal
        mapCoordy%arM(1,v_src%id_nbr, iLev, v_src%ixDispGrd,v_src%iyDispGrd) = &
                     & mapCoordy%arM(1,v_src%id_nbr, iLev, v_src%ixDispGrd,v_src%iyDispGrd) + &
                     & ptrCoord(2) * fCellTotal
        mapCoordz%arM(1,v_src%id_nbr, iLev, v_src%ixDispGrd,v_src%iyDispGrd) = &
                     & mapCoordz%arM(1,v_src%id_nbr, iLev, v_src%ixDispGrd,v_src%iyDispGrd) + &
                     & ptrCoord(3) * fCellTotal
      end if

      mapEmis%ifColumnValid(v_src%id_nbr,v_src%ixDispGrd,v_src%iyDispGrd) = .true.
      mapEmis%ifGridValid(ilev, v_src%id_nbr) = .true.

    end do  ! iLev dispersion

  end subroutine inject_emission_euler_volc_src


  !**************************************************************************

  subroutine inject_emission_lagr_volc_src(v_src, &
                                         & lpset, arParticleMass, & ! Lagrangian
                                         & ChemRunSetup, &  ! translate emission species to transport
                                         & met_buf, disp_buf, &
                                         & now, timestep, &
                                         & fMassInjected)
    !
    ! Adds the emission flux to Lagrangian structure by starting new particles.
    ! Note that particles fly in the meteorological grid to utilise max of available 
    ! dynamic information and also to have a simpler connection to pressure and omega-wind.
    !
    ! Therefore, all source variables with "dispersion-grid" meaning here mean "meteo grid"
    ! Refer to the project_volc_src_2_grids, where the selection is done.
    !
    implicit none

    TYPE(silam_volcano_source), INTENT(in) :: v_src
    type(Tlagrange_particles_set), INTENT(inout), target :: lpSet
    real, dimension(:), intent (in) :: arParticleMass
    type(TchemicalRunSetup), pointer :: ChemRunSetup
    TYPE(Tfield_buffer), pointer :: met_buf, disp_buf
    type(silja_time), intent(in) :: now
    type(silja_interval), intent(in) :: timestep
    real(r8k), dimension(:), intent(inout) :: fMassInjected

    ! Local variables
    integer :: iSlot, iLev, ix, iy, iSrc, nLevs, iDescr, nSpEmis, iSpDescr, nP, iParticle, iP, &
             & nLevsToStart, iSpEmis, iSpTo, iSpTransp
    real :: fWeightPastSrc, fLevFraction, fTimeScale, fCellTotal, factor, timestep_sec, fDx, fDy
    integer :: nSpeciesTrn
    real, dimension(:), pointer :: xCellSize, yCellSize
    real, dimension(max_levels) :: xSize, ySize, zSize
    real, dimension(max_species) :: fMassTmp
    real, dimension(max_descriptors_in_source) :: fMassTimeCommon
    integer, dimension(max_levels) :: nPartInLev
    logical :: ifFound
    type(TspeciesReference), dimension(:), pointer :: references
    integer :: iSlotStart, iSlotEnd
    real, dimension(max_levels) :: fLevBottom, fLevTop, meteo_heights
    type(TVertInterpStruct), pointer ::  interpCoefVert_void
    type(THorizInterpStruct), pointer ::  interpCoefHoriz_void
    real, dimension(:,:), pointer :: arDyn, arMass
    integer, dimension(:), pointer :: arStatus

    arDyn    => lpSet%lpDyn
    arMass   => lpSet%lpMassTrn
    arStatus => lpSet%lpStatus
    !
    ! Have to be very careful: starting and ending of the source slots can be pretty close
    ! especially in case of adjoint run. The same is true for time variation coefficients,
    ! which can vary sharply. Therefore, the total mass injected to the grid will be
    ! integrated along time.
    !
    nSpeciesTrn = lpset%nSpeciesTrn
    
    references => chemRunSetup%refEmis2Transp_mass
    nSpEmis = size(references)
    !
    ! General parameters of the release
    !
    call determine_release_params(v_src, now, timestep, nSpEmis, &
                                & met_buf, interpCoefHoriz_void, .false., &
                                & interpCoefVert_void, .false., &
                                & fMassTimeCommon, fLevTop(1), fLevBottom(1), fWeightPastSrc, ifFound, &
                                & iSlotStart, iSlotEnd)
    if(error)return
    if(.not. ifFound) return                   ! if any emission found
    !
    ! Get the amount of each species to be released
    !
    fMassTmp(1:nSpeciesTrn) = 0.0
    do iDescr = 1, v_src%nDescriptors  ! all descriptors
      do iSpDescr = 1, v_src%nSpeciesInDescr(iDescr)    ! all species over descriptor
        iSpEmis = v_src%pEmisSpeciesMapping(iSpDescr,iDescr)
        !
        ! This emission species contributes to nRefSpecies transport species
        !
        do iSpTo = 1, ChemRunSetup%refEmis2Transp_mass(ispEmis)%nRefSpecies  
          factor = ChemRunSetup%refEmis2Transp_mass(ispEmis)%fract(iSpTo)   ! fractionation
          iSpTransp = ChemRunSetup%refEmis2Transp_mass(iSpEmis)%indSpeciesTo(iSpTo)  ! to whom
          fMassTmp(iSpTransp) = fMassTmp(iSpTransp) + &
                 & fMassTimeCommon(iDescr) * v_src%fDescr2SpeciesUnit(iSpDescr,iDescr) * factor
        end do
      end do
    end do  ! descr
    !
    ! Lagrangian injection goes particle by particle. At least one particle is
    ! always released, so the particle masses are NOT identical - but an effort is 
    ! made to have them close. Within one release step they are identical.
    !
    ! Mass of a single particle and number of released particles: take the 
    ! low-mass threshold for each emission species, compare it to the rate of
    ! emission for this species, and take the largest ratio as the number of particles.
    ! Other species are distributed among these particles
    !
    nP = 1
    do iSpTransp = 1, nSpeciesTrn
      nP = max(nP, int(fMassTmp(iSpTransp) / arParticleMass(iSpTransp)))
!call msg('fMassTmp(iSpTransp), arParticleMass(iSpTransp)', fMassTmp(iSpTransp), arParticleMass(iSpTransp))
!call msg('ratio:',fMassTmp(iSpTransp) / arParticleMass(iSpTransp), nP)
    end do

    !
    ! Let's get the dispersion cell size
    !
    xCellSize => fu_grid_data(meteo_cell_x_size_fld)  ! temporary use of xSize variable
    yCellSize => fu_grid_data(meteo_cell_y_size_fld)  ! temporary use of ySize variable
    iP = v_src%ixDispGrd + (v_src%iyDispGrd-1) * nx_meteo
    fDx = 1./xCellSize(iP)
    fDy = 1./yCellSize(iP)     ! end of temporary use of variables
    !
    ! Particles are injected layer-by-layer of the v_src (no connection to dispersion layer)
    ! First determine the parameters and layers
    !
    timestep_sec = fu_sec(timestep)
    if(error)return
    !
    ! many levels. Set the variables to be used further
    !
    nLevsToStart = 0
    do iLev = 1, nz_meteo
      if(v_src%levFractDispVert(iLev) .eps. 0.0)cycle ! Disp is meteo here
      nLevsToStart = nLevsToStart + 1

      xSize(nLevsToStart) = v_src%params(iSlotStart)%xy_size * fWeightPastSrc + &            ! size of the term
                            & v_src%params(iSlotEnd)%xy_size * (1. - fWeightPastSrc) + &
               & sqrt( get_kz(fldZ, fldRdown, iLev, v_src%ixDispGrd, v_src%iyDispGrd, fWeightPastSrc)  * &
                                 & abs(timestep_sec) / 2.)
      ySize(nLevsToStart) = xSize(nLevsToStart)
      zSize(nLevsToStart) = v_src%dz_m(iLev)
      fLevBottom(nLevsToStart) = v_src%fzDisp(iLev)   ! relative, in meteo grid
      fLevTop(nLevsToStart) = v_src%fzDisp(iLev+1)    ! relative, in meteo grid
      nPartInLev(nLevsToStart) = nint(nP * v_src%levFractDispVert(iLev))  ! in meteo grid
    end do
    if(nLevsToStart == 0)then
      call set_error('No layers to emit anything','inject_emission_v_src')
      return
    endif
    !
    ! Ensure that the number of particles to start is exact
    !
    do while(sum(nPartInLev(1:nLevsToStart)) /= nP)
      iP = 1
      do iLev = 1, nLevsToStart   ! find the level with the largest number of particles injected
        if(nPartInLev(iLev) > nPartInLev(iP)) iP = iLev
      end do
      nPartInLev(iP) = max(nPartInLev(iP) + (nP - sum(nPartInLev(1:nLevsToStart))), 0)
    enddo
        
!------------------------------------------------------------------------
!
! UNcomment this in order to get the particles in the absolute pressure coordinates
!    !
!    ! Having computed fLevBottom(1:nLevsToStart) and fLevBottom(1:nLevsToStart) in relative indices
!    ! of meteo vertical, now we have to turn them into absolute presure
!    !
!    do iLev = 1, nLevsToStart
!      fLevBottom(iLev) = fu_4d_interpolation(fldPressure, &
!                                           & v_src%fxDispGrd, v_src%fyDispGrd, fLevBottom(iLev), &
!                                           & nx_meteo, ny_meteo, nz_meteo, &
!                                           & fWeightPast, &
!                                          & linear, linear, notallowed)
!      fLevTop(iLev) =  fu_4d_interpolation(fldPressure, &
!                                         & v_src%fxDispGrd, v_src%fyDispGrd, fLevTop(iLev), &
!                                         & nx_meteo, ny_meteo, nz_meteo, &
!                                         & fWeightPast, &
!                                         & linear, linear, notallowed)
!    end do
    !
    ! Make-up the particles
    !
    iParticle = lpset%iFirstEmptyParticle
    do iLev = 1, nLevsToStart
      do iP = 1, nPartInLev(iLev)
        !
        ! Find free space
        !
        do while(arStatus(iParticle) /= int_missing)
          iParticle = iParticle + 1
          if(iParticle >= lpset%nop)then
call msg('Enlarging the number of particles, 1:',  lpset%nop + max(lpset%nop*5/4,nPartInLev(iLev)-iP))
            call enlarge_lagrange_particles_set(lpset, lpset%nop + max(lpset%nop*5/4,nPartInLev(iLev)-iP))
            exit
          endif
        end do
        !
        ! Put the masses
        ! Emit the mass going species by species. Note multi-step re-indexing: from species index in 
        ! descriptor to species index in emission list, then to species index in transport list.
        ! Simialrly, MassTimeCommon is mass of a descriptor cocktail, which needs to be explored
        !
        arMass(1:nSpeciesTrn, iParticle) = fMassTmp(1:nSpeciesTrn) / real(nP)
        !
        ! Position, turbulent movement, starting size.
        ! Note that vertical has to be converted to relative
        !
        arDyn(lp_x, iParticle) = fu_random_number_center(v_src%fxDispGrd, max(0.1,xSize(iLev)) * fDx)
        arDyn(lp_y, iParticle) = fu_random_number_center(v_src%fyDispGrd, max(0.1,ySize(iLev)) * fDy)
        arDyn(lp_z, iParticle) = fu_random_number_boundaries(fLevBottom(iLev), fLevTop(iLev))
        
        arDyn(lp_uT:lp_wT, iParticle) = 0.0    ! turbulent-wind motion
        arDyn(lp_dx, iParticle) = xSize(iLev)    ! metres
        arDyn(lp_dy, iParticle) = ySize(iLev)    ! metres

        !!!!Zero-sized particels cause troubles
        arDyn(lp_dz, iParticle) = max(zSize(iLev), 100.)   ! 100 meters -- minimum size 
        arStatus(iParticle) = v_src%id_nbr

        lpset%nNewPart = lpset%nNewPart + 1
        lpset%newPartList(lpset%nNewPart) = iParticle

        iParticle = iParticle + 1
        if(iParticle >= lpset%nop)then
call msg('Enlarging lpset, 2:',  lpset%nop + max(lpset%nop*5/4,nPartInLev(iLev)-iP))
          call enlarge_lagrange_particles_set(lpset, lpset%nop + max(lpset%nop*5/4,nPartInLev(iLev)-iP))
          exit
        endif
      end do  ! iP
    end do  ! iLev of the v_src
    !
    ! Fix the starting particle. Not exactly (place can be occupied) but better than starting from 1
    !
    lpset%iFirstEmptyParticle = iParticle

  end subroutine inject_emission_lagr_volc_src


    !*********************************************************************************

  real function get_kz(Z, Rdown, iLev, ix, iy, fweightpast)
    !
    ! Returns Kz at the locatoin of the next Lagrangian particle
    !
    IMPLICIT NONE 
    
    ! Imported parameters
    integer, intent(in) :: iLev, ix, iy
    real, intent(in) :: fweightpast
    type(field_4d_data_ptr), pointer :: Z, Rdown

    ! Local variables
    integer :: iCellIndex
    real :: dz

    iCellIndex = ix+(iy-1)*nx_meteo

    if (iLev==1) then 
      dz = (Z%past%p2d(iLev)%ptr(iCellIndex) * fweightpast + &
               & Z%future%p2d(iLev)%ptr(iCellIndex) * (1.-fweightpast))
    else 
      dz = (  (Z%past%p2d(iLev)%ptr(iCellIndex)- Z%past%p2d(iLev-1)%ptr(iCellIndex)) * fweightpast + &
              & Z%future%p2d(iLev)%ptr(iCellIndex) - Z%future%p2d(iLev-1)%ptr(iCellIndex) * (1.-fweightpast))
    endif

    get_kz =  dz / (Rdown%past%p2d(iLev)%ptr(iCellIndex) * fweightpast + &
                                  & Rdown%future%p2d(iLev)%ptr(iCellIndex) * (1.-fweightpast)) 
  end function get_kz
  
  
  !*******************************************************************************
  !*******************************************************************************
  !
  ! Data assimilation of the volcanic source parameters
  !
  !*******************************************************************************
  !*******************************************************************************

  subroutine assimilation_request_volc_src(v_src, arParams, iLastFilledParam)
    !
    ! Returns the list of parameters to be assimilated for this volcanic source
    ! Also stores their index in the global array of parameters
    !
    implicit none
    
    ! Imported parameters
    type(silam_volcano_source), intent(inout) :: v_src
    type(assimParameter), dimension(:), intent(inout) :: arParams     ! (nParamsToAssim)
    integer, intent(inout) :: iLastFilledParam

    ! local variables
    integer :: iMap, iParam, iStat, it, nVals
    !
    ! Go through the assimilation maps: they have been filled in at the initialization step
    !
    do iMap = 1, size(v_src%assimMaps)
      !
      ! We recognise two types of assimilation for volcanoes: scalars and vectors
      ! Note that they are packed for compactness - may be, unnecessarily...
      !
      select case(v_src%assimMaps(iMap)%assimType)

        case(mastin_param_assimilation)
          !
          ! Simple scalar(s), uncorrelated. But several can be in the single pack - then unpack here
          !
          v_src%assimMaps(iMap)%indMain = iLastFilledParam + 1  ! store the starting place
          do iParam = 1, v_src%assimMaps(iMap)%nVals
            iLastFilledParam = iLastFilledParam + 1
            arParams(iLastFilledParam)%dim = 1
            allocate(arParams(iLastFilledParam)%val(1), arParams(iLastFilledParam)%cov(1,1), &
                   & arParams(iLastFilledParam)%conditions(size(v_src%assimMaps(iMap)%conditions)), stat=iStat)
            if(fu_fails(iStat==0,'Failed assimilation parameter allocation','assimilation_request_volc_src'))return
            arParams(iLastFilledParam)%val(1) = v_src%assimMaps(iMap)%ini(iParam)
            arParams(iLastFilledParam)%cov(1,1) = v_src%assimMaps(iMap)%cov(iParam,1)
            arParams(iLastFilledParam)%conditions(:) = v_src%assimMaps(iMap)%conditions(:)
          enddo  ! bunch of scalars

        case(emis_z_t_assimilation)
          !
          ! A z-t matrix: time-dependent vector with z-correlated elements (uncorrelated in time)
          ! We explore it here to a bunch of z-vectors with correlated components
          !
          nVals = v_src%assimMaps(iMap)%nVals
          do it = 1, size(v_src%params)
            iLastFilledParam = iLastFilledParam + 1
            v_src%assimMaps(iMap)%indMain = iLastFilledParam          ! store the place
            arParams(iLastFilledParam)%dim = nVals
            allocate(arParams(iLastFilledParam)%val(nVals), &
                   & arParams(iLastFilledParam)%cov(nVals, nVals), &
                   & arParams(iLastFilledParam)%conditions(size(v_src%assimMaps(iMap)%conditions)), stat=iStat)
            if(fu_fails(iStat==0,'Failed vector-assimilation allocation','assimilation_request_volc_src'))return
            arParams(iLastFilledParam)%val(1:nVals) = v_src%assimMaps(iMap)%ems_zt_ini(1:nVals,it)
            arParams(iLastFilledParam)%cov(1:nVals, 1:nVals) = v_src%assimMaps(iMap)%cov(1:nVals, 1:nVals)
            arParams(iLastFilledParam)%conditions(:) = v_src%assimMaps(iMap)%conditions(:)
          end do  ! time slots

        case default
          call set_error('Unknown assimilated parameter type:' + fu_str(v_src%assimMaps(iMap)%assimType),'assimilation_request_volc_src')
          return
      end select  ! what to assimilate
      
    end do  ! over assimMap
    
  end subroutine assimilation_request_volc_src

  
  !*******************************************************************************
  
  subroutine observe_params_volc_src(v_src, arParams, now)
    !
    ! Provides the values of the current assimilated parameters for the interface array
    ! Obeys their place in the interface array. Mapping already exists and serves as the source of data
    !
    implicit none
    
    ! Imported parameters
    type(silam_volcano_source), intent(inout) :: v_src
    type(assimParameter), dimension(:), intent(inout) :: arParams     ! (nParamsToAssim)
    type(silja_time), intent(in) :: now

    ! Local variables
    integer :: iMap, iv, iSlot1, iSlot2
    real :: fWeight_past
    
    do iMap = 1, size(v_src%assimMaps)
      !
      ! We recognise two types of assimilation for volcanoes: scalars and vectors
      ! Note that they are packed for compactness - may be, unnecessarily...
      !
      select case(v_src%assimMaps(iMap)%assimType)

        case(mastin_param_assimilation)          ! Simple scalar
          arParams(v_src%assimMaps(iMap)%indMain)%val(1) = v_src%assimMaps(iMap)%ini(1)

        case(emis_z_t_assimilation)
          !
          ! A z-t matrix: time-dependent vector with z-correlated elements (uncorrelated in time)
          !
          call getTimeSlots_from_params(v_src%params, now, iSlot1, iSlot2, fWeight_past)
          if(error) return
          do iv = 1, v_src%assimMaps(iMap)%nVals
            arParams(v_src%assimMaps(iMap)%indMain)%val(iv) = &
                             & v_src%assimMaps(iMap)%ems_zt_ini(iv,iSlot1) * fWeight_past + &
                             & v_src%assimMaps(iMap)%ems_zt_ini(iv,iSlot2) * (1. - fWeight_past)
          end do
        case default
          call set_error('Unknown assimilated parameter type:' + fu_str(v_src%assimMaps(iMap)%assimType),'observe_params_volc_src')
          return
      end select  ! what to assimilate
      
    end do  ! over assimMap
    
  end subroutine observe_params_volc_src
  
  
  !*******************************************************************************

  subroutine inject_params_volc_src(v_src, arParams, now)
    !
    ! Stores the new parameter values into the corresponding places of the source.
    ! Note that the actuaql papameters are in the mapping structures, rest is just pointers
    !
    implicit none
    
    ! Imported parameters
    type(silam_volcano_source), intent(inout) :: v_src
    type(assimParameter), dimension(:), intent(inout) :: arParams     ! (nParamsToAssim)
    type(silja_time), intent(in) :: now

    ! Local variables
    integer :: iMap, iSlot1, iSlot2
    logical :: ifMastinAssim
    
    do iMap = 1, size(v_src%assimMaps)
      !
      ! We recognise two types of assimilation for volcanoes: scalars and vectors
      !
      select case(v_src%assimMaps(iMap)%assimType)

        case(mastin_param_assimilation)          ! Simple scalar
        
          v_src%assimMaps(iMap)%ini(1) = arParams(v_src%assimMaps(iMap)%indMain)%val(1)
          ifMastinAssim = .true.

        case(emis_z_t_assimilation)
          !
          ! A z-t matrix: time-dependent vector with z-correlated elements (uncorrelated in time)
          !
          ifMastinAssim = .false.
          call getTimeSlots_from_params(v_src%params, now, iSlot1, iSlot2)
          if(error) return
          if(iSlot1 /= iSlot2)then
            call msg_warning('Cannot inject z-t emission between the slots','inject_params_volc_src')
            call report_volcano_source(v_src)
            call msg('Time now:' + fu_str(now))
            call set_error('Cannot inject z-t emission between the slots','inject_params_volc_src')
          endif
          
          v_src%assimMaps(iMap)%ems_zt_ini(1:v_src%assimMaps(iMap)%nVals,iSlot1) = &
                        & arParams(v_src%assimMaps(iMap)%indMain)%val(1:v_src%assimMaps(iMap)%nVals)
          
        case default
          call set_error('Unknown assimilated parameter type:' + fu_str(v_src%assimMaps(iMap)%assimType),'inject_params_volc_src')
          return
      end select  ! what to assimilate
      
    end do  ! over assimMap
    !
    ! If Mastin parameters were updated, have to recompute the z-t emission matrix
    !
    if(ifMastinAssim) call Mastin_height2emis(v_src, v_src%vert4emis_dispVert)
    
    
  end subroutine inject_params_volc_src
  
  
  !*****************************************************************
  !
  ! Smaller-scale supplementary routines
  !
  !*****************************************************************

  subroutine link_volc_src_to_species(species_list, v_src)
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
    type(silam_volcano_source), intent(inout) :: v_src
    type(silam_species), dimension(:), pointer :: species_list

    call link_src_species_to_given_list(v_src%nDescriptors, &
                                      & v_src%nSpeciesInDescr, &
                                      & v_src%cocktail_descr_lst, &
                                      & v_src%fDescr2SpeciesUnit, &
                                      & v_src%pEmisSpeciesMapping, &
                                      & species_list)

  end subroutine link_volc_src_to_species


  !*******************************************************************

  subroutine add_source_species_volc_src(v_src, speciesLst, nSpecies)
    !
    ! Fills-in the given list with the own species. Checks for the duplicates
    !
    implicit none

    ! Improted parameters
    type(silam_volcano_source), intent(in) :: v_src
    type(silam_species), dimension(:), pointer :: speciesLst
    integer, intent(inout) :: nSpecies

    ! Local variables
    integer :: iDescr, nSpeciesDescr
    type(silam_species), dimension(:), pointer :: pSpecies

    do iDescr = 1, v_src%nDescriptors
      call get_inventory(v_src%cocktail_descr_lst(iDescr), pSpecies, nSpeciesDescr)
      if(error)return
      call addSpecies(speciesLst, nSpecies, pSpecies, nSpeciesDescr)
      if(error)return
    end do

  end subroutine add_source_species_volc_src

  
  !*************************************************************************
  
  subroutine determine_release_params(v_src, now, timestep, nSpecies, &
                                    & met_buf, interpCoefMeteo2DispHoriz, ifMetHorizInterp, &
                                    & interpCoefMeteo2DispVert, ifMetVertInterp, &
                                    & fMassTimeCommon, fLevTop, fLevBottom, fWeightPastSrc, ifFound, &
                                    & iSlotStart, iSlotEnd)
    !
    ! Determines the main release parameters: spatial profile and absolute release of 
    ! each species
    !
    implicit none
    
    ! Imported parameters
    TYPE(silam_volcano_source), INTENT(in) :: v_src
    type(silja_time), intent(in) :: now
    type(silja_interval), intent(in) :: timestep
    type(silja_time) :: local_now
    integer, intent(in) :: nSpecies
    TYPE(Tfield_buffer), POINTER :: met_buf
    type(TVertInterpStruct), intent(in) ::  interpCoefMeteo2DispVert
    type(THorizInterpStruct), intent(in) ::  interpCoefMeteo2DispHoriz
    logical, intent(in) :: ifMetVertInterp
    logical, intent(in) :: ifMetHorizInterp
    real, dimension(:), intent(out) :: fMassTimeCommon
    real, intent(out) :: fLevTop, fLevBottom, fWeightPastSrc
    logical, intent(out) :: ifFound
    integer, intent(out) :: iSlotStart, iSlotEnd

    ! Local variables
    integer :: iSlot, iLev, nLevs, iDescr, iTmp, iMeteo
    real :: fLevFraction, fTimeScale, fCellTotal, factor, timestep_sec, exhTempr, exhSpeed, stackSize
    type(silja_time) :: timeStart, timeEnd

!call msg('in')

    timeStart = now
    timeEnd = timeStart + fu_abs(timestep)
    ifFound = .false.
    
    ! Check that this source is active during this time period and cut the start-end time if needed
    !
    if(timeEnd <= v_src%params(1)%time .or. timeStart >= v_src%params(size(v_src%params))%time) return

    if(timeStart < v_src%params(1)%time) timeStart = v_src%params(1)%time
    if(timeEnd > v_src%params(size(v_src%params))%time) timeEnd = v_src%params(size(v_src%params))%time
    if(error)return
    !
    ! Find the range of slots to be looked at
    !
    call get_overlapping_slots(v_src%params, timestart, timeEnd-timeStart, iSlotStart, iSlotEnd)
    if(error .or. iSlotStart >= iSlotEnd)return  ! nothing, for whatever reasons

    fWeightPastSrc = (v_src%params(iSlotEnd)%time - now) / &
                   & (v_src%params(iSlotEnd)%time - v_src%params(iSlotStart)%time)
    !
    ! Rate has to be integrated over time. We have such function:
    !
    call total_from_v_src_descr_unit(v_src, &
                                   & fMassTimeCommon, &
                                   & timeStart, timeEnd - timeStart, &
                                   & level_missing)
    !
    ! If zero emission, nothing to emit...
    !
    if(all(fMassTimeCommon(1:nSpecies) <= 0.0))then
      ifFound = .false.
      return  ! nothing, for whatever reasons
    endif
    ifFound = .true.
    timestep_sec = fu_sec(timeEnd - timeStart)
    fLevBottom = real_missing   ! not used for volcano
    fLevTop = real_missing

  end subroutine determine_release_params

  
  !***********************************************************************
  
  subroutine Mastin_height2emis(v_src, verticalTo)
    !
    ! From injection top height to total emission + profile using Mastin relationto for total amount
    !
    implicit none
  
    ! imported parameters
    TYPE(silam_volcano_source), INTENT(inout) :: v_src
    type(silam_vertical), intent(in) :: verticalTo
    
    ! local variables
    type(silam_vertical) :: vertMushroom
    integer :: iDescr, it, nLevsToActive
    real :: fTotalMass
    
    ! A disputable decision: the same formula is applied to all descriptors, with descriptor-specific 
    ! plume top. Strictly speaking, the plume top is for total ash
    !
    do it = 1, size(v_src%params)
      do iDescr = 1, v_src%nDescriptors
        !
        ! Inverse of H = a * V^b and scale with the fraction of specific cocktail
        ! Note that v_src%params(it)%rate_descr_unit is just a convenient place holder for injection top
        !
        fTotalMass = v_src%fMastin_cocktail_fraction(iDescr) * &
                   & (v_src%params(it)%rate_descr_unit(iDescr) / v_src%fMastin_height_scaling) ** &
                                                                    & (1. / v_src%fMastin_height_pwr)
        !
        ! Need to make a vertical for the mushroom. Two layers: hat and stem
        !
        call set_vertical((/fu_set_layer_between_two( &
                               & fu_set_level(constant_height, &
                                      & v_src%params(it)%rate_descr_unit(iDescr) * &
                                                       & (1.- v_src%fMastin_fract_height_hat)), &
                               & fu_set_level(constant_height, 0.0)), &
                          & fu_set_layer_between_two( &
                               & fu_set_level(constant_height, v_src%params(it)%rate_descr_unit(iDescr)), &
                               & fu_set_level(constant_height, &
                                      & v_src%params(it)%rate_descr_unit(iDescr) * &
                                                       & (1.-v_src%fMastin_fract_height_hat))) &
                         & /), &
                        & vertMushroom)
        !
        ! Now, this mushroom needs to be projected to verticalTo. First, reproject fractions, then scale
        !
        call reproject_verticals(vertMushroom, &
                               & (/1.0-v_src%fMastin_fract_mass_hat, v_src%fMastin_fract_mass_hat/), &
                               & verticalTo, v_src%levFractDispVert, &  !v_src%ems3d_dispVert(iDescr,:,it), &
                               & v_src%fzDisp, nLevsToActive, &
                               & ifMassCentreInRelUnit = .true.)
        if(error)return

        v_src%ems3d_dispVert(iDescr,:,it) = v_src%levFractDispVert(:) * fTotalMass

      end do   ! iDescr
    end do   ! it

  end subroutine Mastin_height2emis

  
  !*************************************************************************
  
  subroutine source_2_map_volc_src(v_src, dataPtr, id)
    !
    ! Projects the point source to the map, having the given id as a 
    ! template: the id deterines the grid, level and substance name.
    ! Since the map does not have any aerosol stuff inside - the modes are either
    ! sumed-up or a given one is picked.
    !
    implicit none

    ! Imported parameters
    type(silam_volcano_source), intent(inout) :: v_src
    real, dimension(:), intent(inout) :: dataPtr
    type(silja_field_id), intent(in) :: id

    ! Local variables
    real :: xOut, yOut
    integer :: nSpecies, nx, ny, i, iDescr
    integer, dimension(max_species) :: species_index
    real, dimension(max_species) :: amounts
    type(silam_species), dimension(:), pointer :: species
    logical :: ifFound
    !
    ! Stupidity check
    !
    if(fu_fails(v_src%defined == silja_true, 'Undefined point source given','source_2_map_volc_src'))return
    if(fu_fails(defined(id), 'Undefined id given','source_2_map_volc_src'))return
    if(fu_fails(defined(fu_grid(id)),'Undefined grid given','source_2_map_volc_src'))return
    if(fu_number_of_gridpoints(fu_grid(id)) > size(dataPtr))then
      call msg_warning('Too small data array given for the grid of the id')
      call report(id)
      call msg('Size of the data array: ', size(dataPtr))
      call set_error('Too small data array given','source_2_map_volc_src')
      return
    endif
    !
    ! Find the grid location of the source. A trick: the point source is always
    ! in geographical co-ordinates, while the output grid can be whatever
    !
    call project_point_to_grid(v_src%lon, v_src%lat, fu_grid(id), xOut, yOut)
    call grid_dimensions(fu_grid(id),nx, ny)

    if(xOut <= 0.5 .or. xOut >= nx+0.5)then
      call msg('Volcano source is outside the x-limits.x=', xOut)
      call msg_warning('Volcano source is outside the x-limits.Skipping','source_2_map_volc_src')
      return
    endif
    if(yOut <= 0.5 .or. yOut >= ny+0.5)then
      call msg('Volcano source is outside the y-limits;y=', yOut)
      call msg_warning('Volcano is outside the y-limits.Skipping','source_2_map_volc_src')
      return
    endif
    !
    ! Rate has to be integrated over time. We have three independently varying
    ! parameters - substance mass fraction in the cocktail, aerosol size distribution
    ! (but not the number of modes) and vertical overlap of the layers, all linear in time. 
    !
    ! Rate has to be integrated over time. We have such function:
    ! total_from_p_src_species_unit(a_src, start, duration, species, amounts, ifRateOnly)
    !
    nullify(species)
    nSpecies = 0
    call total_from_v_src_species_unit(v_src, &
                                     & species, nSpecies, amounts, &
                                     & fu_accumulation_start_time(id), &
                                     & fu_accumulation_length(id), &
                                     & fu_level(id))
    if(error)return
    if(sum(amounts(1:nSpecies)) < 1e-10)then  ! quite arbitrary number but 1e-10 is indeed small
      return
    endif
    nullify(species)
    nSpecies = 0
    !
    ! Once the full inventory is returned, have to select the specific species for each descriptor.
    ! In fact, this is almost the worst possible solution but we cannot afford keeping all the
    ! maps for all the sources
    !
    ifFound = .false.
    do iDescr = 1, v_src%nDescriptors
      call get_inventory(v_src%cocktail_descr_lst(iDescr), species, nSpecies)
      species_index(iDescr) = fu_index(fu_species(id), species, nSpecies)
      if(species_index(iDescr) >= 1 .and. species_index(iDescr) <= nSpecies) ifFound = .true.
    end do

    if(.not. ifFound) return
    !
    ! Preparatory work is over
    !
    i = int(int(xOut+0.5) + (int(yOut+0.5)-1)*nx + 0.5)
    do iDescr = 1, v_src%nDescriptors
      if(species_index(iDescr) < 1)cycle    ! Some descriptors may not have the needed species
      dataPtr(i) = dataPtr(i) + amounts(species_index(iDescr)) * &
                              & v_src%fDescr2SpeciesUnit(species_index(iDescr),iDescr)
    end do

  end subroutine source_2_map_volc_src

  
  !*****************************************************************************************

  subroutine project_volc_src_2_grids(vs, gridMeteo, gridDisp)
    !
    ! Just finds the grid-position in the given grid and stores it to the
    ! position_disp_grd
    !
    implicit none

    ! Imported parameters
    type(silam_volcano_source), intent(inout) :: vs
    type(silja_grid), intent(in) :: gridMeteo, gridDisp

    ! Local variables
    real :: fX, fY
    integer :: nx, ny

    !
    ! Just create new position from the original one with modified horizontal coordinates
    ! ATTENTION. If the source is emitting to Lagrangian environment, it takes meteo grid
    ! If it emits in Eulerian environment - take dispersion grid.
    !
    if(vs%ifEmissionLagrangian)then
      call project_point_to_grid(vs%lon, vs%lat, gridMeteo, vs%fXDispGrd, vs%fYDispGrd)
      call grid_dimensions(gridMeteo,nx,ny)
    else
      call project_point_to_grid(vs%lon, vs%lat, gridDisp, vs%fXDispGrd, vs%fYDispGrd)
      call grid_dimensions(gridDisp,nx,ny)
    endif

    if(error)then
      vs%fXDispGrd = real_missing
      vs%fYDispGrd = real_missing
      vs%ixDispGrd = int_missing
      vs%iyDispGrd = int_missing
      call set_error('Failed to project the source','project_volc_src_2_grids')
    else
      vs%ixDispGrd = nint(vs%fXDispGrd)
      vs%iyDispGrd = nint(vs%fYDispGrd)
    endif
    !
    ! Check that it is inside the grid
    !
    vs%if_inside_domain = vs%ixDispGrd >= 0.5 .and. vs%ixDispGrd <= nx + 0.5 .and. &
                        & vs%iyDispGrd >= 0.5 .and. vs%iyDispGrd <= ny + 0.5

  end subroutine project_volc_src_2_grids


  !********************************************************************
  
  subroutine create_volc_src_cont_grd(v_src, grid_template, ifVerbose, ifMinimal, ifExtended)
    !
    ! This routine creates the grid similar to the given grid_template
    ! covering all sources of the em_source
    !
    implicit none

    ! Imported parameters
    type(silam_volcano_source), intent(in) :: v_src
    type(silja_grid), intent(inout) :: grid_template
    logical, intent(in) :: ifVerbose, ifMinimal
    logical, intent(out) :: ifExtended

    ! Local variables
    integer :: iSrc, nx, ny, iCell
    real :: x, y

    if(fu_fails(v_src%defined == silja_true, 'Undefined source','create_volc_src_cont_grd'))return

    if(defined(grid_template))then
      !
      ! For each source we scan its central points and extending the template_grid 
      ! if these points are outside. For point and bomb sources it is only one point,
      ! while area source has plenty of them.
      !
      if(ifMinimal)then
        !
        ! The requested grid is to be minimal-size just covering the current source
        ! We ignore the current size of the grid and its resolution but still keep the 
        ! projection and grid type
        !
        call project_point_to_grid(v_src%lon, v_src%lat, grid_template, x, y)
        if(error)return
        
        call make_minimal_grid(grid_template, x, y)
        
      else
        !
        ! Just check the source to be inside the grid
        !
        call grid_dimensions(grid_template, nx,ny)

        call project_point_to_grid(v_src%lon, v_src%lat, grid_template, x, y)
        if(x<1 .or. x>nx .or. y<1 .or. y>ny) then
          ifExtended = .true.
          if(ifVerbose)then
            call msg('Grid nx and source x:',nx,x)
            call msg('Grid ny and source y:',ny,y)
            call msg('Extending the grid for the point source:',iSrc)
            call report(v_src)
          endif
          call extend_grid_to_coordinates(grid_template, x, y)
        else
          ifExtended = .false.
        end if
      endif  ! if minimal grid is requested
    else
      !
      ! Grid_template is undefined. Create it using this source as the starting point
      !
      grid_template = fu_set_grid('', lonlat, pole_geographical, &
                                & v_src%lon-0.01, v_src%lat-0.01, &
                                & 3, 3, 0.01, 0.01)   ! nx, ny, dx, dy
      if(error)return
    endif  ! if defined grid_template
  end subroutine create_volc_src_cont_grd

  
  !*******************************************************************

  subroutine prepare_volc_src_vert_params(v_src, vertMeteo, vertDisp, vertMetric)
    !
    ! Projects the vertical parameters of the source to the given vertical
    !
    implicit none

    ! Imported parameter
    type(silam_volcano_source), intent(inout) :: v_src
    type(silam_vertical), intent(in) :: vertMeteo, vertDisp, vertMetric

    ! Local variable
    integer :: iz, i
    !
    ! Lagrangian environment works in meteo grid, Eulerian - in dispersion
    !
    if(v_src%ifEmissionLagrangian)then
      v_src%vert4emis_dispVert = vertMeteo
    else
      v_src%vert4emis_dispVert = vertDisp
    endif
    
    v_src%nzDispVert = fu_NbrOfLevels(v_src%vert4emis_dispVert)
    if(v_src%emission_vs_vertical == MastinRelation)then
      !
      ! from injection top height to total emission + profile
      !
      if(.not. allocated(v_src%ems3d_dispVert))then
        allocate(v_src%ems3d_dispVert(v_src%nDescriptors, v_src%nzDispVert, size(v_src%params)), &
               & v_src%fzDisp((v_src%nzDispVert+1)), v_src%levFractDispVert(v_src%nzDispVert), stat=i)
        if(fu_fails(i==0,'Failed allocation of ems3d_dispVert','prepare_volc_src_vert_params'))return
      endif

      call Mastin_height2emis(v_src, vertMetric) !, v_src%vert4emis_dispVert)
      if(error)return
    else
      !
      ! If emis3d given, have to reproject. 
      ! Volcano vertical is fixed in time
      !
      allocate(v_src%levFractDispVert(v_src%nzDispVert), &  ! fractions
             & v_src%fzDisp(v_src%nzDispVert + 1), &        ! bottoms and tops => n+1
             & v_src%dz_m(v_src%nzDispVert), stat=i)        ! thickness, m, of the layers
      if(fu_fails(i==0,'Failed to allocate dispersion-vertical lagrangian level fractions', &
                     & 'prepare_volc_src_vert_params'))return
      v_src%levFractDispVert(1:v_src%nzDispVert) = 0.0
      v_src%fzDisp(1:v_src%nzDispVert + 1) = 0.0
      do iz = 1, v_src%nzDispVert
        v_src%dz_m(iz) = fu_layer_thickness_m(fu_level(v_src%vert4emis_dispVert,iz))
      end do
      !
      ! Dispersion vertical can be hybrid. For the source we defined a matching vertical
      ! in meters. This ensures that emission density (per vertical meter) is conserved in
      ! the vertical reprojection.
      !
      call reproject_verticals(v_src%vert4emis, v_src%levFraction, &   ! vertical from, fractions from
                             & vertMetric, v_src%levFractDispVert, & ! vertical to, fractions to
                             & v_src%fzDisp, v_src%nzDispVert, &  ! mass centres, number of non-zero levels
                             & ifMassCentreInRelUnit=.true.)
    endif  ! Mastin or direct ems3d

  end subroutine prepare_volc_src_vert_params


  !********************************************************************************************

  logical function fu_volc_emis_owned_quantity(v_src, quantity)
    !
    ! Checks whether the specific quantity is handled by the emission sources exclusively
    ! Such quantities may be related to emission species somewhere else
    !
    implicit none

    ! Imported parameters
    type(silam_volcano_source), intent(in) :: v_src
    integer, intent(in) :: quantity
    !
    ! The volcano source does not have own quantities yet
    !
    select case(quantity)
      case default
        fu_volc_emis_owned_quantity = .false.
    end select
  end function fu_volc_emis_owned_quantity

  
  !*******************************************************************************
  !
  ! Encapsulation of various small stuff
  !
  !*******************************************************************************
  !==================================================================
  function  fu_name_v_src(v_src) result(nm)
    implicit none
    character(len=clen) :: nm
    type(silam_volcano_source), intent(in) :: v_src
    nm = v_src%src_nm
  end function fu_name_v_src
  !==================================================================
  function  fu_sector_v_src(v_src) result(nm)
    implicit none
    character(len=clen) :: nm
    type(silam_volcano_source), intent(in) :: v_src
    nm = v_src%src_nm
  end function fu_sector_v_src
  !==================================================================
  integer function fu_source_id_nbr_of_v_src(v_src)
    implicit none
    type(silam_volcano_source), intent(in) :: v_src
    fu_source_id_nbr_of_v_src = v_src%id_nbr
  end function fu_source_id_nbr_of_v_src
  !==================================================================
  integer function fu_source_nbr_of_v_src(v_src)
    implicit none
    type(silam_volcano_source), intent(in) :: v_src
    fu_source_nbr_of_v_src = v_src%src_nbr
  end function fu_source_nbr_of_v_src
  !==================================================================
  SUBROUTINE report_volcano_source(vs)
    IMPLICIT NONE
    ! Imported parameters with intent IN:
    TYPE(silam_volcano_source), INTENT(in) :: vs
    call msg('================= VOLCANO SOURCE =======================')
    call msg('The volcano source description')
    call store_volc_src_as_namelist(vs, run_log_funit, .true., .true.)
    call msg('********************* end of VOLCANO source **********************')
  END SUBROUTINE report_volcano_source
  !==================================================================
  FUNCTION fu_volcano_source_start_time(v_src)
    IMPLICIT NONE
    TYPE(silja_time) :: fu_volcano_source_start_time
    TYPE(silam_volcano_source), INTENT(in) :: v_src
    fu_volcano_source_start_time = v_src%params(1)%time
  END FUNCTION fu_volcano_source_start_time
  !==================================================================
  FUNCTION fu_volcano_source_end_time(v_src)
    IMPLICIT NONE
    TYPE(silja_time) :: fu_volcano_source_end_time
    TYPE(silam_volcano_source), INTENT(in) :: v_src
    fu_volcano_source_end_time = v_src%params(size(v_src%params))%time
  END FUNCTION fu_volcano_source_end_time
  !==================================================================
  FUNCTION fu_volcano_source_duration(v_src)
    IMPLICIT NONE
    TYPE(silja_interval) :: fu_volcano_source_duration
    TYPE(silam_volcano_source), INTENT(in) :: v_src
    fu_volcano_source_duration = v_src%params(size(v_src%params))%time - v_src%params(1)%time
  END FUNCTION fu_volcano_source_duration
  !==================================================================
  subroutine getTimeSlots_of_volcano_source(v_src, time, iSlot1, iSlot2)
    ! Finds the time slots surrounding the given time
    IMPLICIT NONE
    ! Imported parameters
    TYPE(silam_volcano_source), INTENT(in) :: v_src
    type(silja_time), intent(in) :: time
    integer, intent(out) :: iSlot1, iSlot2
    call getTimeSlots_from_params(v_src%params, time, iSlot1, iSlot2)
  END subroutine  getTimeSlots_of_volcano_source
  !=========================================================================  
  function fu_SlotTime_of_volcano_source(v_src, iSlot)
    implicit none
    type(silja_time) :: fu_SlotTime_of_volcano_source
    TYPE(silam_volcano_source), INTENT(in) :: v_src
    integer, intent (in) :: iSlot
    if(fu_fails(iSlot <= size(v_src%params),'Too big slot number','fu_SlotTime_of_volcano_source'))return
    fu_SlotTime_of_volcano_source = v_src%params(iSlot)%time
  end function fu_SlotTime_of_volcano_source
  !=========================================================================  
  function fu_cocktail_descr_of_v_src (v_src) result(descrLst)
    implicit none
    type(Tcocktail_descr), dimension(:), pointer :: descrLst
    type(silam_volcano_source), intent(in), target :: v_src
    descrLst => v_src%cocktail_descr_lst
  end function fu_cocktail_descr_of_v_src
  !==================================================================
  FUNCTION fu_NbrOfTimeSlots_of_v_src(v_src)
    IMPLICIT NONE
    integer :: fu_NbrOfTimeSlots_of_v_src
    TYPE(silam_volcano_source), INTENT(in) :: v_src
    fu_NbrOfTimeSlots_of_v_src = size(v_src%params)
  END FUNCTION fu_NbrOfTimeSlots_of_v_src
  !==================================================================
  logical function fu_volcano_source_defined(v_src)
    implicit none
    type(silam_volcano_source), intent(in) :: v_src
    fu_volcano_source_defined = v_src%defined == silja_true
  end function fu_volcano_source_defined
  
END MODULE source_terms_volcano

