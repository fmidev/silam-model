MODULE source_terms_area

  ! This module contains the general description of area source-term of SILAM
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
  ! All units: SI unless otherwise stated explicitly
  ! 
  ! All degrees: real numbers in degrees (degrees and parts per
  ! hundred/parts per thousand , NOT in degrees and minutes)
  ! 
  ! Language: ANSI standard Fortran 90
  !
  ! Code owner: Mikhail Sofiev, FMI
  ! 
  ! Modules used:

  USE source_terms_time_params

  IMPLICIT NONE
  private

  
  
  ! The public functions and subroutines available in this module:
  
  public reserve_area_source
  public add_input_needs_area_source
  public fill_a_src_from_namelist
  public store_as_namelist
  public write_area_src_from_mass_map

  public defined
  public set_undef
  public report
  public fu_name
  public fu_sector
  public fu_grid
  public fu_start_time
  public fu_end_time
  public fu_duration
  public total_amt_species_unit
  public fu_NbrOfTimeSlots
  public getTimeSlots
  public fu_SlotTime
  public fu_cocktail_descr
  public get_lp_release_height
  public fu_source_id_nbr
  public fu_source_nbr
  public fu_nbr_of_disp_grd_cells
  public source_2_map_area_source
  public project_a_src_second_grd     ! projects to the given grid and stores new cell set
  public prepare_src_vert_params
  public fu_useTimeVarCoef
  public fu_ifMultiLevelFixed
  public fu_n_cells
  public fu_total_cell_value
  public fu_cell_values
  public link_source_to_species
  public add_source_species_a_src
  public get_cell_grid_coords
  public get_cell_geo_coords
  public create_src_contain_grd_a_src
  public inject_emission_area_source_field
  public count_emission_area_source_cellist
  public inject_emission_area_source_cellist
  public check_cells
  public force_source_into_grid
  public check_time_params_a_src
  public unpack_a_src_from_nc
  public init_a_src_TZ_index
  public fu_meteodep_model

  !
  ! PRIVATE routines
  !
  private fill_a_src_from_namelist_v4
  private fill_a_src_from_namelist_v2_3
  private store_a_src_as_namelist
  private report_area_source
  private force_a_src_into_grid
  private fu_source_nbr_of_a_src

  private fu_area_source_defined
  private fu_grid_of_area_source
  private fu_area_source_start_time
  private fu_area_source_end_time
  private fu_area_source_duration
  private fu_name_a_src
  private fu_sector_a_src
  private get_lp_release_height_a_src
  private fu_NbrOfTimeSlots_of_a_src
  private total_from_a_src_descr_unit
  private total_amt_species_unit_a_src
  private getTimeSlots_of_area_source
  private fu_SlotTime_of_area_source
  private fu_cocktail_descr_of_a_src
  private fu_source_id_nbr_of_a_src
  private fu_nbr_of_disp_grd_cells_a_src
  private prepare_src_vert_params_a_src
  private get_cell_grid_coords_a_src
  private get_cell_geo_coords_a_src
  private fu_total_cell_value_a_src
  private fu_cell_values_a_src
  private fu_n_cells_a_Src
  private fu_useTimeVarCoef_a_src
  private fu_ifMultiLevelFixed_a_src
  private link_a_src_to_species

  private fu_meteodep_model_a_src

  ! Generic names and operator-interfaces of some functions:

  INTERFACE defined
    module procedure fu_area_source_defined
  END INTERFACE

  INTERFACE set_undef
    module procedure set_area_source_undefined
  END INTERFACE

  INTERFACE report
    MODULE PROCEDURE report_area_source
  END INTERFACE

  INTERFACE store_as_namelist
    MODULE PROCEDURE store_a_src_as_namelist
  END INTERFACE

  INTERFACE fu_name
    module procedure fu_name_a_src
  END INTERFACE

  INTERFACE fu_sector
    module procedure fu_sector_a_src
  END INTERFACE

  interface fu_grid
    module procedure fu_grid_of_area_source
  end interface

  INTERFACE fu_start_time
    MODULE PROCEDURE fu_area_source_start_time
  END INTERFACE

  INTERFACE fu_end_time
    MODULE PROCEDURE fu_area_source_end_time
  END INTERFACE

  INTERFACE fu_duration
    MODULE PROCEDURE fu_area_source_duration
  END INTERFACE

  interface total_amt_species_unit
    module procedure total_amt_species_unit_a_src
  end interface

  interface fu_NbrOfTimeSlots
    module procedure fu_NbrOfTimeSlots_of_a_src
  end interface 

  interface getTimeSlots
    module procedure getTimeSlots_of_area_source
  end interface

  interface fu_SlotTime
    module procedure fu_SlotTime_of_area_source
  end interface

  interface fu_cocktail_descr
    module procedure fu_cocktail_descr_of_a_src
  end interface

  INTERFACE get_lp_release_height
    MODULE PROCEDURE get_lp_release_height_a_src
  END INTERFACE

  interface fu_source_id_nbr
    module procedure fu_source_id_nbr_of_a_src
  end interface

  interface fu_source_nbr
    module procedure fu_source_nbr_of_a_src
  end interface

  interface get_cell_grid_coords
    module procedure get_cell_grid_coords_a_src
  end interface

  interface get_cell_geo_coords
    module procedure get_cell_geo_coords_a_src
  end interface

  interface  fu_total_cell_value
    module procedure fu_total_cell_value_a_src
  end interface

  interface  fu_cell_values
    module procedure fu_cell_values_a_src
  end interface

  interface  fu_n_cells
    module procedure fu_n_cells_a_Src
  end interface

  interface fu_nbr_of_disp_grd_cells
    module procedure fu_nbr_of_disp_grd_cells_a_src
  end interface

  interface prepare_src_vert_params
    module procedure prepare_src_vert_params_a_src
  end interface

  interface fu_useTimeVarCoef
    module procedure fu_useTimeVarCoef_a_src
  end interface

  interface fu_ifMultiLevelFixed
    module procedure fu_ifMultiLevelFixed_a_src
  end interface

  interface link_source_to_species
    module procedure link_a_src_to_species
  end interface

  interface force_source_into_grid
    module procedure force_a_src_into_grid
  end interface
  
  interface fu_meteodep_model
    module procedure fu_meteodep_model_a_src
  end interface

  integer, private, parameter :: nc_tag_len = 4

  integer, private, parameter :: LocalTimeIndex = -55555 ! Not actual index
  integer, private, parameter :: SolarTimeIndex = -66666
  !
  ! The area source is defined in a very similar way as the point source.
  ! The only difference is that instead of position, the grid is used with 
  ! appropriate vector of x-y-value. Not the matrix - the vector is smaller.
  !

  !==================================================================

  TYPE silam_area_source
    PRIVATE
    !
    ! General environment
    !
    CHARACTER(len=clen) :: src_nm, sector_nm ! Name of the area source and sector
    integer :: src_nbr, id_nbr               ! Just a source and id numbers in a WHOLE source list
    TYPE(silja_grid) :: grid
    type(silam_release_parameters), dimension(:), pointer :: params  => null()
    type(Tcocktail_descr), dimension(:), pointer :: cocktail_descr_lst  => null()   ! for the whole source
    integer :: nCells, nCellsDispGrd, nDescriptors, nLevsDispVert
    logical :: ifUseTimeVarCoef, ifRateVary, ifMultiLevelFixed
    logicaL, dimension(:), allocatable :: ifSpecies
    type(silam_vertical) :: vertLevs, vertLevsDispVert
    integer, dimension(:), pointer :: nSpeciesInDescr => null()
    real, dimension(:,:), pointer :: fDescr2SpeciesUnit  => null()
    integer, dimension(:,:), pointer :: pEmisSpeciesMapping  => null()
    real, dimension(:), pointer :: levFraction  => null() , levFractDispVert  => null(), &
         & fzDisp  => null()
    real, dimension(:,:), pointer :: indHour    => null()          ! Hourly rate indices
    real, dimension(:,:), pointer :: indDay     => null()          ! Daily rate indices
    real, dimension(:,:), pointer :: indMonth    => null()         ! Monthly rate indices
    
    real, dimension(:), pointer :: fDescrTotals  => null()  ! basic descriptor mass unit
    !
    ! Cells for the source given as list of cells
    !
    real, dimension(:), pointer :: cell_fx => null(), cell_fy => null(), &
                                 & cellDispGrd_fx => null(), cellDispGrd_fy => null()
    real, dimension(:,:), pointer :: cell_val => null(), cellDispGrd_val =>null()
    integer, dimension(:), allocatable :: CellTZindex ! Index of timezone used for local and solar time
    !
    ! Environment for field-based emission definition
    !
    logical, public :: ifFieldGiven   ! switch between the cell list and field
    ! GrADS/NetCDF/..., ref to input structures:
    integer, dimension(:), pointer :: indFile4Descr =>null(), uField => null()
    type(silam_fformat), dimension(:), pointer :: FieldFormat => null()
    character(len=fnlen), dimension(:), pointer :: chFieldFNms =>null() ! names/template strings of the binary files 
    type(grads_template), dimension(:), pointer :: chFieldFNmTemplates =>null() ! templates of the binary files 
    type(silja_time), dimension(:), pointer :: times =>null()
!    type(silja_shopping_list), pointer :: pShopLst
    type(silja_field_3d_pointer), dimension(:), pointer :: pFldEmis  =>null()   ! main storage place
    type(silja_logical) :: ifHorizInterp, ifVertInterp  ! if interpolation to dispersion grid is needed
    type(THorizInterpStruct), pointer :: pHorizIS =>null()  ! horizontal interpolation structure
    type(TVertInterpStruct), pointer :: pVertIS  =>null()   ! vertical interpolation strcture
    logical :: ifFluxes = .False. !! True if field icomes as fluxes kg/s/m2 / false if Emission Rates (kg/s)
    integer :: slotOffForFile = int_missing !!! in time-vert-resolving sources offset of time index to be used for a file
                                            !!Grads can handle this on its own
                                            !! in Netcdf files named according to end_of_period slotOffForFile=1 
                                            !! in Netcdf files named according to start or mid of period slotOffForFile=0
                                            !! Checked by name corresponding to the start of first slot in a file and to the end of
                                            !! last slot in the file.
                                            !! Note that it might not correspond to the actual time_label_position in a file
    character(len=tz_name_length) :: tz_name  ! Timezone name 
    integer :: tz_index ! Index of timezone to use with hourly/daily indices
                    ! can be LocalTimeIndex  or SolarTimeIndex
    !
    ! Environment for meteo-dependent variation
    logical, public :: ifTemprDependent  = .FALSE.          ! monthly variation coefficients are NOT ignored
    type(TmeteoDepRules)  :: rulesMeteoDep  = meteoDepRulesMissing
    logical, public :: ems2wholeABL = .false.
    type(silja_3d_field), allocatable :: fldEmis3d  ! if injection into ABL, space for profile
    
    type(silja_logical) :: defined = silja_false
  END TYPE silam_area_source
  public silam_area_source

  !
  ! Pointer to the area source. Needed to encapsulate the area source in the array
  !
  type a_src_ptr
    TYPE(silam_area_source) :: a_src 
  end type a_src_ptr
  public a_src_ptr


CONTAINS

  !**************************************************************************

  subroutine reserve_area_source(a_src, &      ! Src to initialise
                                  & iSrcNbr, &    ! Src number in the grand list
                                  & iSrcIdNbr, &  ! SrcID number
                                  & nDescriptors, & ! number of chemical descriptors to be reserved
                                  & nCells)
    !
    ! Initialises the source:
    ! - stores the reference information: source number and source ID number
    ! - stores the total number of chemical descriptors that will be stored in the source
    ! - nullifies the source dynamic arrays 
    ! - set a few basic internal variables of the source
    !
    implicit none

    ! Imported parameters
    type(silam_area_source), intent(inout) :: a_src
    integer, intent(in) :: iSrcNbr, iSrcIdNbr, nDescriptors, nCells
    
    integer :: idescr

    ! Local variables
    integer :: iTmp
    !
    ! Nullify the basic variables
    !
    a_src%src_nm = ''
    a_src%sector_nm = ''
    nullify(a_src%params)
    nullify(a_src%levFraction)
    nullify(a_src%fzDisp)
    nullify(a_src%levFractDispVert)

    !
    ! Main source parameters - enough to identify it in the global information list
    !
    a_src%src_nbr = iSrcNbr
    a_src%id_nbr = iSrcIdNbr
    a_src%nDescriptors = nDescriptors
    a_src%nCells = 0    ! No cells have been filled-in
    !
    ! Will not need them any time soon
    !
    nullify(a_src%cellDispGrd_fx)
    nullify(a_src%cellDispGrd_fy)
    nullify(a_src%cellDispGrd_val)
    !
    ! Will need these arrays while reading the source. Note that the size may not be final 
    !
    if(nCells > 0)then
      nullify(a_src%cell_fx);    call set_array_size(a_src%cell_fx, nCells, 0.0)
      nullify(a_src%cell_fy);    call set_array_size(a_src%cell_fy, nCells, 0.0)
      nullify(a_src%cell_val);  call set_array_size(a_src%cell_val, a_src%nDescriptors, nCells, 0.0)
    else
      nullify(a_src%cell_fx)
      nullify(a_src%cell_fy)
      nullify(a_src%cell_val)
    endif

    allocate(a_src%cocktail_descr_lst(nDescriptors), a_src%fDescrTotals(nDescriptors), &
           & a_src%indHour(nDescriptors,24), a_src%indDay(nDescriptors,7), &
           & a_src%indMonth(nDescriptors,12), stat=iTmp)
    if(iTmp /= 0)then
      call set_error('Failed to allocate space for basic coefficients','reserve_area_source')
      return
    endif
    a_src%fDescrTotals = real_missing
    do idescr = 1, nDescriptors
      call set_missing(a_src%cocktail_descr_lst(idescr))
    end do
    a_src%indHour = int_missing
    a_src%indDay = int_missing
    a_src%indMonth = int_missing
    a_src%tz_index = 0
    a_src%tz_name = ''
    nullify(a_src%nSpeciesInDescr)
    nullify(a_src%fDescr2SpeciesUnit)
    nullify(a_src%pEmisSpeciesMapping)
!    nullify(a_src%pShopLst)
    
    a_src%ifTemprDependent = .false.
    a_src%ems2wholeABL = .false.
    !
    ! Finally, mark the source as incomplete
    !
    a_src%defined = silja_false

  end subroutine reserve_area_source


  !***************************************************************************

  subroutine add_input_needs_area_source(a_src, q_met_dynamic, q_met_st, q_disp_dynamic, q_disp_st)
    !
    ! Returns input needs for the emission/transformation routines. 
    !
    implicit none

    ! Imported parameters
    type(silam_area_source), intent(in) :: a_src
    integer, dimension(:), intent(inout) :: q_met_dynamic, q_met_st, q_disp_dynamic, q_disp_st
    integer :: iTmp
    !
    ! Meteo data may be needed if the plume rise computations show up
    ! So far, nothing
    !
    if(a_src%tz_index == LocalTimeIndex)then !Use local time
      iTmp = fu_merge_integer_to_array(timezone_index_flag, q_met_st)
    endif
    !
    ! Actual type of dependence and needed input is in source_terms_time_params module
    !
    if(a_src%ifTemprDependent) call add_input_needs_meteo_dependence(a_src%rulesMeteoDep, &
                                                                   & q_met_dynamic, q_met_st, &
                                                                   & q_disp_dynamic, q_disp_st)
    ! and the whole-ABL emission
    if(a_src%ems2wholeABL) iTmp = fu_merge_integer_to_array(abl_height_m_flag, q_met_dynamic)

    return

  end subroutine add_input_needs_area_source


  !****************************************************************
  
  subroutine fill_a_src_from_namelist(nlSrc, a_src, chAreaSrcFileVersion, src_dir_name)
    !
    ! A dispatcher of the reading task to corrspondng subs, 
    ! depending on the source version
    ! "2_1" is fully compatible with "2",  Made to break the run when temperature dependence 
    !                                      is present in source but not supported (e.g. by rarlier Silam)
    !
    IMPLICIT NONE
    !
    ! Imported parameters
    TYPE(silam_area_source), intent(inout) :: a_src
    type(Tsilam_namelist), pointer :: nlSrc
    character(len=*), intent(in) :: chAreaSrcFileVersion, &  !at present, "2", "2_1",  "3", '3_1' or "4"
              & src_dir_name ! for finding the netcdf data file
           

    integer :: iTmp
    character(len=*), parameter :: sub_name = 'fill_a_src_from_namelist'
    
    select case(chAreaSrcFileVersion)
      case('2','2_1')
        call fill_a_src_from_namelist_v2_3(nlSrc, a_src, 2, src_dir_name)
      case('3', '3_1')
        call fill_a_src_from_namelist_v2_3(nlSrc, a_src, 3, src_dir_name)
      case('4')
        call fill_a_src_from_namelist_v4(nlSrc, a_src, src_dir_name)
      case default
        call set_error('Unknown a_src version: "'// trim(chAreaSrcFileVersion)//'"', sub_name)
      end select
      
    !
    ! Whatever the source version is, they all can depend on meteorology
    ! Two possible dependencies: heating degree day (a variety of functions) and emission
    ! into the whole ABL
    !
    a_src%ifTemprDependent = fu_str_u_case(fu_content(nlSrc,'if_temperature_dependent_emission')) == 'YES'
    
    a_src%ems2wholeABL = fu_str_u_case(fu_content(nlSrc,'highest_emission_injection_height')) == 'BLH'
    
    if(a_src%ifTemprDependent)then
      call get_temp_dep_from_namelist(a_src%rulesMeteoDep, nlSrc)
      if(error)then
        call set_error('Failed meteo dependence for:' + a_src%src_nm + '___' + a_src%sector_nm, &
                     & sub_name)
        a_src%defined = silja_false
        return
      endif
      if (any(chAreaSrcFileVersion == (/'2','3'/))) then
          call msg_warning("Meteo-dependent source with no bakward lock", sub_name)
          call msg("Please change a_src")
      endif
      if (chAreaSrcFileVersion == '4') then
        call report_area_source(a_src)
        call set_error("area_source_4 does not support meteo dependence ", sub_name)
        return
      endif

    endif
      
    if(.not. defined(a_src)  .or. a_src%nDescriptors < 1)then
      call report_area_source(a_src)
      call set_error('Undefined or empty area source v'+trim(chAreaSrcFileVersion), &
                                        & sub_name)
      return
    endif

    if(a_src%nCells < 1 )then
      call report_area_source(a_src)
      call msg_warning('Zero-cells area source v'+trim(chAreaSrcFileVersion), &
                                        & sub_name)
    endif
      
  end subroutine fill_a_src_from_namelist


  ! ***************************************************************

  subroutine fill_a_src_from_namelist_v2_3(nlSrc, a_src, iAreaSrcFileVersion, src_dir_name)
    ! Reads and sets one area source term v.2 from an external file.
    ! A trick:
    ! Since the area source is potentially huge, it can be moved to a 
    ! separate file. In this case in the main file we will just have 
    ! a name of the file to read. 
    ! Note also that we rely on the nDescriptors value to be set, as well the 
    ! number of the source in the source info list (available for the general
    ! source only).
    !
    ! All units: SI, unless otherwise is stated
    !
    IMPLICIT NONE
    !
    ! Imported parameters
    TYPE(silam_area_source), intent(inout) :: a_src
    type(Tsilam_namelist), pointer :: nlSrc
    integer, intent(in) :: iAreaSrcFileVersion   ! at present, 2, 3
    character(len=*), intent(in) :: src_dir_name ! for finding the netcdf data file

    ! Local variables
    type(silam_sp) :: spContent
    integer :: iLocal, iTmp, jTmp, nx, ny, nItems, iStat, iNlCell, nCellsInNamelist, iCell, iCurCell, &
             & iDescr, nDescriptorInParLine, iModeNbr, iItem
    integer, dimension(:), pointer :: indDescriptorOfParLine    ! places of descriptors in descriptor list
    type(Tsilam_nl_item_ptr), dimension(:), pointer :: ptrItems
    real :: fLon, fLat, fX, fY
    real, dimension(:), pointer :: fVals
    logical :: ifGeoCoords, ifFound
    type(silja_level) :: level
    character(len=substNmLen) :: chSubstNm
    character(len=*), parameter :: sub_name = 'fill_a_src_from_namelist_v2_3'


    ! A bit of preparations
    !
    !$ if (omp_in_parallel()) call set_error("OMP unsafe",sub_name)
    !$ if (error) return

    spContent%sp => fu_work_string()
    indDescriptorOfParLine => fu_work_int_array()
    fVals => fu_work_array()
    if(error)return

    nullify(ptrItems)
    !
    ! If the source header has not yet been defined, define it. Otherwise check that we 
    ! got the right one
    !
    if(a_src%defined == silja_true)then
      !
      ! The source has already been partly defined. Check namelist to be the right one
      !
      if(.not. fu_same_header())then
        call msg_warning('The header and the namelist are different',sub_name)
        call report(a_Src)
        call report(nlSrc)
        call set_error('The header and the namelist are different',sub_name)
        return
      endif
    else
      !
      ! The source has only been initialised but nothing yet read. Set the header from the namelist
      !
      call set_header()
      if(error)return
    endif

!call report(a_src%vertLevs)

    !
    ! Timeslot parameters
    !
    select case(iAreaSrcFileVersion)
      case(2)
        call get_items(nlSrc, 'par_str', ptrItems, nItems)

      case(3)
        call get_items(nlSrc, 'par_str_area', ptrItems, nItems)
      case default
        call msg('Wrong area source file version (must be 2 or 3):',iAreaSrcFileVersion)
        call set_error('Wrong area source file version',sub_name)
        return
    end select

    if(nItems < 1)then
      call report(nlSrc)
      call set_error('No parameter items in namelist',sub_name)
      return
    endif

    if(.not. a_src%defined == silja_true)then
      !
      ! If nothing has been read yet, param vector is null. Have to allocate it
      !
      allocate(a_src%params(nItems), stat=iTmp)
      if(fu_fails(iTmp == 0, 'Failed to allocate parameters for area source', &
                           & sub_name))return
      call set_missing(a_src%params, nItems)
    endif

    spContent%sp=fu_content(nlSrc,'vertical_unit')
    if(error)return
    indDescriptorOfParLine(1:a_src%nDescriptors) = int_missing
    
    a_src%tz_name = fu_content(nlSrc,'source_timezone')
    if( a_src%tz_name/= '')then 
      if(fu_str_u_case(a_src%tz_name) == 'LOCAL')then
        ! Use local time at the source location
        a_src%tz_index = LocalTimeIndex
      elseif (fu_str_u_case(a_src%tz_name) == 'SOLAR')then
        a_src%tz_index = SolarTimeIndex
      else 
        a_src%tz_index =  timezone_index_by_name(a_src%tz_name)
        if (a_src%tz_index < 1) then
          call set_error("Failed to find Timezone index for '"//trim(a_src%tz_name)//"'", sub_name)
          return
        endif
      endif
    else 
      a_src%tz_name = 'Solar'
      a_src%tz_index = SolarTimeIndex
!      call msg('No timezone for source: '+a_src%src_nm+', '+ a_src%sector_nm +'. Assuming solar time.')
    endif 

    do iTmp = 1, nItems
      select case(iAreaSrcFileVersion)
        case(2)
          chSubstNm = fu_content(nlSrc,'emitted_substance')
          iModeNbr = fu_content_int(nlSrc,'emitted_size_mode_nbr')
          call decode_time_parameter_v4 (fu_content(ptrItems(iTmp)), &  ! line to decode
                                       & iTmp, &                        ! param element index to write to
                                       & spContent%sp, &                ! unit of the vertical
                                       & fu_content(nlSrc,'release_rate_unit'), & ! source release rate unit
                                       & a_src%params, &                ! param array
                                       & a_src%cocktail_descr_lst, &    ! list of source cocktail descriptors
                                       & area_source, &                 ! type of the source
                                       & a_src%nDescriptors, &          ! max nbr of descriptors in the source
                                       & indDescriptorOfParLine, nDescriptorInParLine, & 
                                       & chSubstNm, & 
                                       & iModeNbr, &
                                       & a_src%ifSpecies)
        case(3)
          call decode_time_parameter_v5 (fu_content(ptrItems(iTmp)), &
                                       & iTmp, &
                                       & spContent%sp, &
                                       & fu_content(nlSrc,'release_rate_unit'), & ! source release rate unit
                                       & a_src%params, &
                                       & a_src%cocktail_descr_lst, &    ! list of source cocktail descriptors
                                       & area_source, &   ! type of the source
                                       & a_src%nDescriptors, &
                                       & indDescriptorOfParLine, nDescriptorInParLine, &
                                       & '', fVals, &    ! void here
                                       & a_src%ifSpecies)
        case default
          call msg('Wrong area source file version (must be 2 or 3):',iAreaSrcFileVersion)
          call set_error('Wrong area source file version',sub_name)
          return
      end select
      if(error)then
        call set_error('Failed source:' + a_src%src_nm + '_' + a_src%sector_nm, &
                     & sub_name)
        return
      endif
    end do   ! time slots
    !
    ! Copy the time moments of params to times array.
    !
    allocate(a_src%times(size(a_src%params)), stat=iTmp)
    if(fu_fails(iTmp == 0, 'Failed times array allocation',sub_name))return
    do iTmp = 1, size(a_src%times)
      a_src%times(iTmp) = a_src%params(iTmp)%time
    end do
    !
    ! It is useful to know if the composition varies with time - just a speedup
    !
    SLOTS: do iTmp = 1, size(a_src%params)-1
      do iDescr = 1, a_src%nDescriptors
        a_src%ifRateVary = .not. (a_src%params(iTmp)%rate_descr_unit(iDescr) .eps. &
                                & a_src%params(iTmp+1)%rate_descr_unit(iDescr))
        if(a_src%ifRateVary)exit SLOTS
      end do
    end do SLOTS

    !
    ! If time variation indices are present - read and set them, otherwise - set all to 1
    !
    select case(iAreaSrcFileVersion)
      case(2)
        !
        ! Same species can already be met in other files. Check for the same variation
        !
        if(a_src%indHour(indDescriptorOfParLine(1),1) >= 0.0)then ! something reasonable found
          call set_error('The cocktail:' + fu_name(a_src%cocktail_descr_lst(indDescriptorOfParLine(1))) + &
                 & '- has already been assigned time variation for area source v2:' + a_src%src_nm + &
                 & '_' + a_src%sector_nm,sub_name)
        endif
        !
        ! Put the variation to the corresponding indDescriptorOfParLine
        !
        spContent%sp = fu_content(nlSrc,'hour_in_day_index')
        if(spContent%sp /= '')then
          read(unit=spContent%sp, iostat=iLocal, fmt=*) &
                                 & (a_src%indHour(indDescriptorOfParLine(1),iTmp),iTmp=1,24)
          if(iLocal /= 0)then
            call set_error('Failed to read  hourly variation:' + spContent%sp,'set_header')
            return
          endif
!          call msg('Hourly variation indices are read. sum=',real_value=sum(a_src%indHour(1:24)))
        else
          a_src%indHour(indDescriptorOfParLine(1),:) = 1
        endif
        spContent%sp = fu_content(nlSrc,'day_in_week_index')
        if(spContent%sp /= '')then
          read(unit=spContent%sp, iostat=iLocal, fmt=*) &
                                    & (a_src%indDay(indDescriptorOfParLine(1),iTmp),iTmp=1,7)
          if(iLocal /= 0)then
            call set_error('Failed to read  day-in-week variation:' + spContent%sp,'set_header')
            return
          endif
!          call msg('Daily variation indices are read. sum=',real_value=sum(a_src%indDay(1:24)))
        else
          a_src%indDay(indDescriptorOfParLine(1),:) = 1
        endif
        spContent%sp = fu_content(nlSrc,'month_in_year_index')
        if(spContent%sp /= '')then
          read(unit=spContent%sp, iostat=iLocal, fmt=*) &
                                         & (a_src%indMonth(indDescriptorOfParLine(1),iTmp),iTmp=1,12)
          if(iLocal /= 0)then
            call set_error('Failed to read  monthly variation:' + spContent%sp,'set_header')
            return
          endif
!          call msg('Monthly variation indices are read. sum=',real_value=sum(a_src%indMonth(1:24)))
        else
          a_src%indMonth(indDescriptorOfParLine(1),:) = 1
        endif
        a_src%ifUseTimeVarCoef = a_src%ifUseTimeVarCoef .or. &
                               & any(a_src%indHour(indDescriptorOfParLine(1),1:24) < 0.999) .or. &
                               & any(a_src%indHour(indDescriptorOfParLine(1),1:24) > 1.0001) .or. &
                               & any(a_src%indDay(indDescriptorOfParLine(1),1:7) < 0.999) .or. &
                               & any(a_src%indDay(indDescriptorOfParLine(1),1:7) > 1.0001) .or. &
                               & any(a_src%indMonth(indDescriptorOfParLine(1),1:12) < 0.999) .or. &
                               & any(a_src%indMonth(indDescriptorOfParLine(1),1:12) > 1.0001)
      case(3)
        !
        ! AREA_SOURCE_3: several descriptors in one file.
        ! Same species can already be met in other files. 
        !
        do iDescr = 1, nDescriptorInParLine
          if(a_src%indHour(indDescriptorOfParLine(iDescr),1) >= 0.0)then ! something reasonable found
            call set_error('The species:' + &
                   & fu_name(a_src%cocktail_descr_lst(indDescriptorOfParLine(iDescr))) + &
                   & '- has already been assigned time variation for the area source 3:' + &
                   & a_src%src_nm + '_' + a_src%sector_nm,sub_name)
          endif
        enddo
        !
        ! Put the variation to the corresponding indDescriptorOfParLine
        !
        call get_items(nlSrc,'hour_in_day_index',ptrItems,nItems)
        jTmp = nItems
        do iDescr = 1, nDescriptorInParLine
          a_src%indHour(indDescriptorOfParLine(iDescr),:) = 1.
          do iItem = 1, nItems
            spContent%sp = fu_content(ptrItems(iItem))
            read(unit=spContent%sp,fmt=*)chSubstNm
            if(chSubstNm == trim(fu_name(a_src%cocktail_descr_lst(indDescriptorOfParLine(iDescr))))) then
              read(unit=spContent%sp, iostat=iLocal, fmt=*) chSubstNm, &
                                     & (a_src%indHour(indDescriptorOfParLine(iDescr),iTmp),iTmp=1,24)
              jTmp = jTmp - 1
              if(iLocal /= 0)then
                call set_error('Failed to read  hourly variation:' + spContent%sp,'set_header')
                return
              endif
              exit
            endif
          end do
        end do
        if (fu_fails(jTmp == 0, "Not all hour_in_day_index (asrc3) were used", sub_name)) return 

        call get_items(nlSrc,'day_in_week_index',ptrItems,nItems)
        jTmp = nItems
        do iDescr = 1, nDescriptorInParLine
          a_src%indDay(indDescriptorOfParLine(iDescr),:) = 1.
          do iItem = 1, nItems
            spContent%sp = fu_content(ptrItems(iItem))
            read(unit=spContent%sp,fmt=*)chSubstNm
            if(chSubstNm == trim(fu_name(a_src%cocktail_descr_lst(indDescriptorOfParLine(iDescr))))) then
              read(unit=spContent%sp, iostat=iLocal, fmt=*) chSubstNm, &
                                     & (a_src%indDay(indDescriptorOfParLine(iDescr),iTmp),iTmp=1,7)
              jTmp = jTmp - 1
              if(iLocal /= 0)then
                call set_error('Failed to read  daily-in-week variation:' + spContent%sp,'set_header')
                return
              endif
              exit
            endif
          end do
        end do
        if (fu_fails(jTmp == 0, "Not all day_in_week_index (asrc3) were used", sub_name)) return 

        call get_items(nlSrc,'month_in_year_index',ptrItems,nItems)
        jTmp = nItems
        do iDescr = 1, nDescriptorInParLine
          a_src%indMonth(indDescriptorOfParLine(iDescr),:) = 1.
          do iItem = 1, nItems
            spContent%sp = fu_content(ptrItems(iItem))
            read(unit=spContent%sp,fmt=*)chSubstNm
            if(chSubstNm == trim(fu_name(a_src%cocktail_descr_lst(indDescriptorOfParLine(iDescr))))) then
              read(unit=spContent%sp, iostat=iLocal, fmt=*) chSubstNm, &
                                     & (a_src%indMonth(indDescriptorOfParLine(iDescr),iTmp),iTmp=1,12)
              jTmp = jTmp - 1
              if(iLocal /= 0)then
                call set_error('Failed to read month-in-year variation:' + spContent%sp,'set_header')
                return
              endif
              exit
            endif
          end do
        end do
        if (fu_fails(jTmp == 0, "Not all month_in_year_index (asrc3) were used", sub_name)) return 

        do iDescr = 1, nDescriptorInParLine
          a_src%ifUseTimeVarCoef = a_src%ifUseTimeVarCoef .or. &
                           & any(a_src%indHour(indDescriptorOfParLine(iDescr),1:24) < 0.999) .or. &
                           & any(a_src%indHour(indDescriptorOfParLine(iDescr),1:24) > 1.0001) .or. &
                           & any(a_src%indDay(indDescriptorOfParLine(iDescr),1:7) < 0.999) .or. &
                           & any(a_src%indDay(indDescriptorOfParLine(iDescr),1:7) > 1.0001) .or. &
                           & any(a_src%indMonth(indDescriptorOfParLine(iDescr),1:12) < 0.999) .or. &
                           & any(a_src%indMonth(indDescriptorOfParLine(iDescr),1:12) > 1.0001)
        enddo

      case default
        call set_error('Wrong AREA_SOURCE version',sub_name)
        return
    end select

    !
    ! Now, values
    !
    ! They can be given in the declared grid or normal geographical co-ordinates
    !
    if(fu_content(nlSrc,'coordinate_of_values') == 'GEOGRAPHICAL')then
      ifGeoCoords = .true.
    elseif(fu_content(nlSrc,'coordinate_of_values') == 'GRID_INDICES')then
      ifGeoCoords = .false.
    else
      call set_error('The content of coordinate_of_values must be GEOGRAPHICAL or GRID_INDICES', &
                   & sub_name)
      return
    endif
    
    ! The values can be either included in the namelist or referred into a *packed* netcdf data
    ! file, not to be confused with the files for gridded volume sources. The data is handled below:
    ! 
    call read_values()

    !
    ! Things have been done but do NOT destroy the namelist - it was NOT you who made it, stupid!
    !
    call free_work_array(spContent%sp)
    call free_work_array(indDescriptorOfParLine)
    call free_work_array(fVals)
    call destroy_items(ptrItems)
!    !
!    ! Finally, some reporting and checking
!    !
!    if(a_src%sector_nm == '')then
!      call msg(fu_connect_strings('Area source name:',a_src%src_nm))
!    else
!      call msg(fu_connect_strings('Area source name & sector:',a_src%src_nm,',',a_src%sector_nm))
!    endif
    !
    ! Negative values are not allowed. Check this here.
    !
    do iCell = 1, a_src%nCells
      do iDescr = 1, a_src%nDescriptors
        if(a_src%cell_val(iDescr,iCell) < 0.)then
          call msg('Descriptor, cell index:', iDescr, iCell)
          call msg('X, Y-coord:',a_src%cell_fx(iCell), a_src%cell_fy(iCell))
          call msg('Value:',a_src%cell_fx(iCell), a_src%cell_val(iDescr,iCell))
          call set_error('This source has negative emission rate',sub_name)
          return
        endif
      end do
    end do
    !
    ! Finally, nullify the field-emission structures
    !
    a_src%ifFieldGiven = .false.   ! switch between the cell list and field
!    a_src%ifTimeResolved = .false. ! if field is time-resolving (or use time_param / indXXX)
!    a_src%ifVerticalResolved = .false. ! if field is vertically-resolving (or use vertLevs)
    nullify(a_src%FieldFormat)    ! GrADS / NetCDF / ...
    nullify(a_src%uField)          ! Reference to the input struture
    nullify(a_src%chFieldFNms)     ! names of the binary files 
    nullify(a_src%chFieldFNmTemplates)     ! name templates of the binary files 
    nullify(a_src%indFile4Descr)
    a_src%ifHorizInterp = silja_false
    a_src%ifVertInterp = silja_false   ! if interpolation to dispersion grid is needed
    nullify(a_src%pHorizIS)  ! horizontal interpolation structure
    nullify(a_src%pVertIS)    ! vertical interpolation strcture

    a_src%defined = silja_true

    CONTAINS

    !====================================================================================

    subroutine set_header()
      !
      ! Sets the basic parameters of the source
      !
      implicit none
      !
      ! Source name
      !
      a_src%src_nm = fu_Content(nlSrc, 'source_name')
      a_src%sector_nm = fu_Content(nlSrc, 'source_sector_name')
      a_src%ifUseTimeVarCoef = .false. 

      ! Grid parameters
      !
      a_src%grid = fu_set_grid(nlSrc)
      if(error)return
      if(.not. defined(a_src%grid))then
        call msg_warning('Undefined grid of an area source','set_header')
        call report(nlSrc)
        call set_error('Undefined grid of the area source','set_header')
        return
      endif
      call grid_dimensions(a_src%grid, nx, ny)
      if(error)return
      !
      ! If requested, read the multi-level structure that will be constant in time
      !
      spContent%sp = fu_str_u_case(fu_content(nlSrc,'vertical_distribution'))
      if(spContent%sp == 'SINGLE_LEVEL_DYNAMIC' .or. spContent%sp == '')then
        !
        ! Do nothing - all is in param_strings
        !
        call set_missing(a_src%vertLevs, .true.)
        nullify(a_src%levFraction)
        a_src%ifMultiLevelFixed = .false.

      elseif(spContent%sp == 'MULTI_LEVEL_FIXED')then
        !
        ! Read levels and set the source vertical
        !
        call get_items(nlSrc, 'vert_level', ptrItems, nItems)
        spContent%sp=fu_content(nlSrc,'vertical_unit')
        if(nItems < 1)then
          call set_error('No levels for vertical_distribution = MULTI_LEVEL_FIXED', &
                       & 'set_header')
          return
        else
          allocate(a_src%levFraction(nItems))
          call set_named_level_with_fract(ptrItems(1), level, a_src%levFraction(1), spContent%sp)
          if(error)return
          call set_vertical(level, a_src%vertLevs)
          if(error)return
          do iTmp = 2, nItems
            call set_named_level_with_fract(ptrItems(iTmp), level, a_src%levFraction(iTmp), spContent%sp)
            call add_level(a_src%vertLevs, level)
            if(error)return
          end do
        endif
        if (abs(sum(a_src%levFraction(:)) - 1.) > 1e-3 ) then
          call msg("Fractions",a_src%levFraction(:))
           call report(a_src%vertLevs,.True.)
          call set_error("Strange level fractions", "set_header")
        endif
        a_src%ifMultiLevelFixed = .true.
      else
        call msg('vertical distribution must be MULTI_LEVEL_FIXED or SINGLE_LEVEL_DYNAMIC')
        call set_error(fu_connect_strings('Strange vertical_distribution:',spContent%sp), &
                     & 'set_header')
        return
      endif
      a_src%ifFieldGiven = .false. ! binary used only for source v4
      
    end subroutine set_header


    !==============================================================================

    logical function fu_same_header()
      !
      ! Compares the current source header and the one defined in the namelist
      !
      implicit none

      ! Local variables
      type(silam_release_parameters), dimension(1), target :: paramTmp
      type(silam_release_parameters), dimension(:), pointer :: paramPtr
      
      paramPtr =>null()


      fu_same_header = .false.

      !
      ! Name
      !
      if(.not. (fu_str_u_case(a_src%src_nm) == fu_str_u_case(fu_Content(nlSrc, 'source_name'))))then
        call msg_warning(fu_connect_strings('Source names are different:', &
                                         & a_src%src_nm,',',fu_Content(nlSrc, 'source_name')), &
                       & 'fu_same_header')
        call set_error('Cannot collect species between different sources','fu_same_header')
        return
      endif

      ! Sector
      if(.not. (fu_str_u_case(a_src%sector_nm) == &
                                      & fu_str_u_case(fu_Content(nlSrc, 'source_sector_name'))))then
        call msg_warning(fu_connect_strings('Source sector names are different:', &
                                         & a_src%sector_nm,',',fu_Content(nlSrc, 'source_sector_name')), &
                       & 'fu_same_header')
        call set_error('Cannot collect species between different sources','fu_same_header')
        return
      endif

      ! grid
      if(.not. (fu_set_grid(nlSrc) == a_src%grid))then
        call msg_warning('Grids are different','collect_chemical_2_a_src')
        call msg('New grid from the namelist')
        call report(fu_set_grid(nlSrc))
        call msg('Grid in the area source:')
        call report(a_src%grid)
        call set_error('Cannot collect species between different sources','fu_same_header')
        return
      endif

      ! Timeslot parameters - just check the number of slots. The list of descriptors will be
      ! enlarged later on
      !
      select case(iAreaSrcFileVersion)
        case(2)
          call get_items(nlSrc, 'par_str', ptrItems, nItems)
        case(3)
          call get_items(nlSrc, 'par_str_area', ptrItems, nItems)
        case default
          call msg('Wrong area source file version (must be 2 or 3):',iAreaSrcFileVersion)
          call set_error('Wrong area source file version','set_a_src_from_namelist')
          return
      end select
!      call get_items(nlSrc, 'par_str', ptrItems, nItems)
      if(nItems /= size(a_src%params))then
        call msg('Different number of time slots in a_src and namelist:', nItems, real(size(a_src%params)))
        call set_error('Cannot collect species between different sources','fu_same_header')
        return
      endif

      !
      ! And now - the parameters themselves
      !
      paramPtr => paramTmp
      call set_missing(paramPtr, 1)

      do iTmp = 1, nItems
        select case(iAreaSrcFileVersion)
          case(2)
            call decode_time_parameter_v4(fu_content(ptrItems(iTmp)), &  ! parameter itself
                                        & 1, &                           ! always 1 to get to paramTmp
                                        & 'm', &                         ! arbitrary, vertical is NOT CHECKED
                                        & fu_content(nlSrc,'release_rate_unit'), & ! source release rate unit
                                        & paramPtr, &                    ! receiving temporary
                                        & a_src%cocktail_descr_lst, &    ! list of source cocktail descriptors
                                        & area_source, &                 ! type of the source
                                        & a_src%nDescriptors, &
                                        & indDescriptorOfParLine, nDescriptorInParLine, & 
                                        & fu_content(nlSrc,'emitted_substance'), &
                                        & fu_content_int(nlSrc,'emitted_size_mode_nbr'), &
                                        & a_src%ifSpecies)
          case(3)
            call decode_time_parameter_v5(fu_content(ptrItems(iTmp)), &  ! parameter itself
                                        & 1, &                           ! always 1 to get to paramTmp
                                        & 'm', &                         ! arbitrary, vertical is NOT CHECKED
                                        & fu_content(nlSrc,'release_rate_unit'), & ! source release rate unit
                                        & paramPtr, &                    ! receiving temporary
                                        & a_src%cocktail_descr_lst, &    ! list of source cocktail descriptors
                                        & area_source, &                 ! type of the source
                                        & a_src%nDescriptors, &
                                        & indDescriptorOfParLine, nDescriptorInParLine, &
                                        & '', a_src%fDescrTotals, a_src%ifSpecies)
         case default
            call msg('Wrong area source file version (must be 2 or 3):',iAreaSrcFileVersion)
            call set_error('Wrong area source file version','fu_same_header')
            return
        end select
        if(error)return

        if(.not.(paramTmp(1)%time == a_src%params(iTmp)%time) .or. &
         & .not.(paramTmp(1)%xy_size .eps. a_src%params(iTmp)%xy_size) .or. &
         & .not.(paramTmp(1)%z_velocity .eps. a_src%params(iTmp)%z_velocity) .or. &
         & .not.(paramTmp(1)%tempr .eps. a_src%params(iTmp)%tempr))then
          call msg_warning(fu_connect_strings('Failed comparison for par_str:',fu_content(ptrItems(iTmp))), &
                         & 'fu_same_header')
          call report(a_src)
          call set_error('Different parameters', 'fu_same_header')
          return
        endif

      end do  ! params items

      fu_same_header = .true.

    end function fu_same_header


    !************************************************************************************

    subroutine read_values()
      implicit none
      integer :: num_cells
      
      integer :: nctag, ncsize_in_nl
      character(len=fnlen) :: ncfilename_hat, ncfilename
      logical :: read_from_nc
      real, dimension(:), pointer :: lats, lons, values

       character(len=*), parameter :: sub_name = 'read_values'
       lats =>null()
       lons =>null()
       values =>null()
      
      spcontent%sp = fu_content(nlSrc, 'netcdf_data')
      if (spContent%sp == '') then
        read_from_nc = .false.
        call get_items(nlSrc, 'val', ptrItems, num_cells)
        if(num_cells < 1)then
          call msg('No emission cells in source (fill_a_src_from_namelist_v2_3): ' &
                                & // trim(fu_content(nlSrc,'source_name')))
          return
        endif
      else
        read_from_nc = .true.
        read(unit=spContent%sp, iostat=istat, fmt=*) nctag, ncsize_in_nl, ncfilename_hat
        ncfilename = fu_extend_grads_hat_dir(ncfilename_hat, src_dir_name)
!        call replace_string(ncfilename, '^', trim(src_dir_name)//dir_slash)
        if (fu_fails(istat == 0, 'Failed to parse data item: ' // trim(spcontent%sp), sub_name)) return
        values => fu_work_array()
        lons => fu_work_array()
        lats => fu_work_array()
        call unpack_a_src_from_nc('read', ncfilename, nctag, nDescriptorInParLine, num_cells, values, lons, lats)
        if (error) return
        if (fu_fails(num_cells == ncsize_in_nl, 'Sizes in namelist and netcdf don''t match', sub_name)) return
      end if
          
      call grid_dimensions(a_src%grid,nx,ny)
      if(error)return

      iCurCell = 1

      do iNlCell = 1, num_cells
        if (read_from_nc) then
          flon = lons(iNlcell)
          flat = lats(iNlcell)
          fvals(1) = values(iNlCell)
        else
          !
          ! Careful with reading! Do not destroy the cells, which are not in the list - they can be 
          ! already set from previous files
          !
          spContent%sp = fu_content(ptrItems(iNlCell))
          
          read(unit=spContent%sp, iostat=iStat,fmt=*) fLon, fLat, &
               & (fVals(iDescr), iDescr=1,nDescriptorInParLine)
          if(iStat /= 0)then
            call set_error('Failed to read the line:' + spContent%sp, sub_name)
          endif
          if(any(fVals(1:nDescriptorInParLine) < 0))then
            call msg('Problem:'+spContent%sp,fVals(1), fVals(2))
          endif
        end if
        
        ! Set the read cell into the proper place of the source. Steps:
        ! - find the cell grid indices in the source grid
        ! - if the source is met for the first time, just take the next cell
        ! - if something has been previously read, find the cell to fill in
        !
        if(ifGeoCoords)then
          !
          ! Have to reproject the val line to the source grid
          !
          call project_point_to_grid(fLon, fLat, a_src%grid, fX, fY)
        else
          fX = fLon
          fY = fLat
        endif

        ! Check that the point is inside the source grid. Skip it otherwise
        !
        if(fX < 0.5 .or. fY < 0.5 .or. fX > nx+0.5 .or. fY > ny+0.5)then
          !        call msg('Cell:'+ spContent%sp + ', grid co-ordinates:', fX, fY)
          !        call msg_warning('This cell of area source:' + a_src%src_nm + '-is outside its grid', &
          !                       & sub_name)
          cycle
        endif

        !
        ! Having the cell coordinates in its grid, try to find the corresponding existing cell
        ! or create a new one. Search has to be smart because cells are not sorted anyhow:
        ! (1) if this is the first time we fill the source, take the new cell
        ! (2) if the source has already been partly filled:
        !   2a go from the current position in cell array to the end -may be, just the next one is OK
        !   2b if not found, go from beginning until the current position
        !   2c store the found cell as new current position
        !
        if(a_src%defined == silja_true)then
          !
          ! Try to find something
          !
          ifFound = .false.
          !
          ! Search from curent position until the end
          !
          do iCell = iCurCell, a_src%nCells
            if(a_src%cell_fx(iCell) .eps. fX)then
              if(a_src%cell_fy(iCell) .eps. fY)then
                !
                ! Cell found, update the values. Note:  fVals is zero except for just-read values,
                ! which are in the right positions
                !
                iCurCell = iCell
                ifFound = .true.
                do iDescr = 1, nDescriptorInParLine
                  a_src%cell_val(indDescriptorOfParLine(iDescr),iCell) = &
                       & a_src%cell_val(indDescriptorOfParLine(iDescr),iCell) + fVals(iDescr)
                end do
                exit   ! from search cycle
              endif
            endif
          end do
          !
          ! If not the next cell, may be, it is at the beginning on the list
          !
          if(.not. ifFound)then
            !
            ! Search from the beginning until the current position
            !
            do iCell=1, iCurCell
              if(a_src%cell_fx(iCell) .eps. fX)then
                if(a_src%cell_fy(iCell) .eps. fY)then
                  !
                  ! Cell found, update the values
                  !
                  iCurCell = iCell
                  ifFound = .true.
                  do iDescr = 1, nDescriptorInParLine
                    a_src%cell_val(indDescriptorOfParLine(iDescr),iCell) = &
                         & a_src%cell_val(indDescriptorOfParLine(iDescr),iCell) + fVals(iDescr)
                  end do
                  exit
                endif
              endif
            end do
          endif    ! ifFound

        else
          !
          ! Nothing to search: the cell has not been read before
          !
          ifFound = .false.

        endif  ! if the first-time source reading
        !
        ! If the cell has not yet been stored for whatever reasons, take th new one
        !
        if(.not. ifFound)then
          !
          ! Take entirely new cell
          !
          if(a_src%nCells == size(a_src%cell_fx))then
            call set_array_size(a_src%cell_fx, a_src%nCells+20, 0.0)
            call set_array_size(a_src%cell_fy, a_src%nCells+20, 0.0)
            call set_array_size(a_src%cell_val, a_src%nDescriptors, a_src%nCells+20, 0.0)
          endif
          if(error)return

          a_src%nCells = a_src%nCells+1
          a_src%cell_fx(a_src%nCells) = fX
          a_src%cell_fy(a_src%nCells) = fY
          do iDescr = 1, nDescriptorInParLine
            a_src%cell_val(indDescriptorOfParLine(iDescr),a_src%nCells) = &
                 & a_src%cell_val(indDescriptorOfParLine(iDescr),a_src%nCells) + fVals(iDescr)
          end do
          iCurCell = a_src%nCells

        endif  ! if the cell is new
      end do  ! num_cells

      if (read_from_nc) then
        call free_work_array(values)
        call free_work_array(lons)
        call free_work_array(lats)
      end if

    end subroutine read_values

  end subroutine fill_a_src_from_namelist_v2_3

    !************************************************************************************
    
  subroutine unpack_a_src_from_nc(action, ncfilename, tag, num_descr, src_size, values, lons, lats)
    use netcdf
    ! Read values for a source packed into netcdf.  The source is identified with its
    ! index (tag), which is given in the header namelist. To avoid closing and re-opening
    ! netcdf files repeatedly, the file is cached until a new file is requested, or the
    ! subroutine is called with 'close' as the first argument.
    implicit none
    character(len=*), intent(in) :: action, ncfilename
    ! action is either:
    ! 'read'  : read the values
    ! 'count' : count the size of source
    ! 'close' : close the cached netcdf file, if any.
    integer, intent(out) :: src_size ! how many val lines the source has
    real, dimension(:), intent(out), optional :: values
    real, dimension(:), intent(out), optional :: lons, lats
    integer, intent(in) :: num_descr, tag
    
    character(len=fnlen), save :: ncfilename_prev = ''
    integer, save :: ncid = int_missing, var_id_val, var_id_lon, var_id_lat, var_id_ind_start, &
         & var_id_size, dim_id_src, num_srcs_in_nc
    integer :: ind_start, src_ind_nc, stat

    character(len=*), parameter :: sub_name = 'unpack_a_src_from_nc'

    if (action == 'close') then
      ! close the netcdf if any is open
      if (ncid /= int_missing) then
        stat = nf90_close(ncid)
        if (fu_fails(stat == NF90_NOERR, 'Failed to close nc file: ' // trim(ncfilename_prev), sub_name)) return
        ncfilename_prev = ''
      end if
      ncid = int_missing
      return
    end if

    if (fu_fails(ncfilename /= '', 'Netcdf filename is blank', sub_name)) return

    if (ncfilename /= ncfilename_prev) then
      if (ncid /= int_missing) then
        stat = nf90_close(ncid)
        if (fu_fails(stat == NF90_NOERR, 'Failed to close nc file: ' // trim(ncfilename_prev), sub_name)) return
      end if
      stat = nf90_open(ncfilename, 0, ncid)
      if (fu_fails(stat == NF90_NOERR, 'Failed to open nc file: ' // trim(ncfilename), sub_name)) return
      stat = nf90_inq_dimid(ncid, 'src', dim_id_src)
      ncfilename_prev = ncfilename

      if (fu_fails(stat == NF90_NOERR, 'Failed to inquire src dim', sub_name)) return
      stat = nf90_inquire_dimension(ncid, dim_id_src, len=num_srcs_in_nc)
      if (fu_fails(stat == NF90_NOERR, 'Failed to get source dim size', sub_name)) return
      if (fu_fails(num_srcs_in_nc > 0, 'Bad num_srcs_in_nc', sub_name)) return
      
      stat = nf90_inq_varid(ncid, 'value', var_id_val)
      if (fu_fails(stat == NF90_NOERR, 'Failed to var id for values', sub_name)) return
      stat = nf90_inq_varid(ncid, 'lon', var_id_lon)
      if (fu_fails(stat == NF90_NOERR, 'Failed to var id for lon', sub_name)) return
      stat = nf90_inq_varid(ncid, 'lat ', var_id_lat)
      if (fu_fails(stat == NF90_NOERR, 'Failed to var id for lat', sub_name)) return

      stat = nf90_inq_varid(ncid, 'src_start_index', var_id_ind_start)
      if (fu_fails(stat == NF90_NOERR, 'Failed to var id for start index', sub_name)) return
      stat = nf90_inq_varid(ncid, 'src_size', var_id_size)
      if (fu_fails(stat == NF90_NOERR, 'Failed to var id source size', sub_name)) return
    end if

    if (fu_fails(num_descr == 1, 'Multi-descriptor netcdfs not supported', sub_name)) return

    if (fu_fails(tag >= 0 .and. tag < num_srcs_in_nc, 'Bad source index (tag)', sub_name)) return

    if (action == 'count') then
      stat = nf90_get_var(ncid, var_id_size, src_size, start=(/tag+1/))
      if (fu_fails(stat == NF90_NOERR, 'Failed to get size', sub_name)) return

    else if (action == 'read') then
      if (fu_fails(present(lons) .and. present(lats) .and. present(values),'Missing arguments',sub_name)) return
      ! To read a source, we need 2 indices: the position of the source in the index arrays
      ! (ind_start, src_size), and the index of the first record for this source. The
      ! first-level index is found with the tag, the seconds index is found using the first index.
      src_ind_nc = tag + 1
      
      stat = nf90_get_var(ncid, var_id_ind_start, ind_start, start=(/src_ind_nc/))
      if (fu_fails(stat == NF90_NOERR, 'Failed to get ind_start', sub_name)) return
      
      stat = nf90_get_var(ncid, var_id_size, src_size, start=(/src_ind_nc/))
      if (fu_fails(stat == NF90_NOERR, 'Failed to get size', sub_name)) return
      
      if (fu_fails(size(values) >= src_size, 'values array too small', sub_name)) return
      if (fu_fails(size(lons) >= src_size, 'lons array too small', sub_name)) return
      if (fu_fails(size(lats) >= src_size, 'lats array too small', sub_name)) return
      
      ! ind_start is given zero-based!
      ! 
      stat = nf90_get_var(ncid, var_id_val, values, count=(/src_size/), start=(/ind_start+1/))
      if (fu_fails(stat == NF90_NOERR, 'Failed to read values', sub_name)) return
      stat = nf90_get_var(ncid, var_id_lon, lons, count=(/src_size/), start=(/ind_start+1/))
      if (fu_fails(stat == NF90_NOERR, 'Failed to read lons', sub_name)) return
      stat = nf90_get_var(ncid, var_id_lat, lats, count=(/src_size/), start=(/ind_start+1/))
      if (fu_fails(stat == NF90_NOERR, 'Failed to read lats', sub_name)) return
      
    else 
      call set_error('Strange action: ' // trim(action), sub_name)
      return
    end if

  end subroutine unpack_a_src_from_nc


  ! ***************************************************************

  subroutine fill_a_src_from_namelist_v4(nlSrc, a_src, src_dir_name)
    !
    ! Reads and sets one area source term v.4 from an external file.
    ! This source version consists of a header and external "binary" file,
    ! which contains the 3D emission field. The "binary" fiel can be ASCII, 
    ! GrADS, GRIB, NetCDF, or even TEST_FIELD.
    !
    ! Note: one source must be described with one header. Repetitions are
    ! not allowed.
    !
    ! All units: SI, unless otherwise is stated
    !
    IMPLICIT NONE
    !
    ! Imported parameters
    TYPE(silam_area_source), intent(inout) :: a_src
    type(Tsilam_namelist), pointer :: nlSrc
    character(len=*), intent(in) :: src_dir_name

    ! Local variables
    type(silam_sp) :: spContent
    integer :: iLocal, iTmp, nItems, iStat, iDescr, nDescriptorInParLine, iItem, nIDs, iLev, nTimes, &
         & ind_file, ind_grid, num_grids, ind_vert, num_verts, prev_unit, nx, ny
    integer :: quantity
    integer, dimension(:), pointer :: indDescriptorOfParLine 
    type(Tsilam_nl_item_ptr), dimension(:), pointer :: ptrItems ! 
    type(silja_grid), dimension(:), pointer :: pGrids
    real :: fLon, fLat, fX, fY
    type(silja_level) :: level
    character(len=clen) :: release_rate_unit, vertical_unit, chTmp
    character(len=fnlen) :: fname, fnameTmp
    type(silam_species) :: speciesTmp
    type(silam_species), dimension(:), pointer :: speciesArPtr
    type(silja_time), dimension(:), pointer :: timesTmp 
    type(silja_rp_1d), dimension(:), pointer :: pArTmp 
    logical :: ifFound, ifTimeVertResolved
    type(silja_field_id) :: idTmp, idTmp1
    type(silja_field_id), dimension(:), pointer :: idList
    real, dimension(:), pointer :: pData 
    type(silam_vertical), dimension(:), pointer :: pVerts 
    character(len=*), parameter :: sub_name = "fill_a_src_from_namelist_v4"

    !
    ! Area source v.4 always refers to the separate-field file, below called as "binary"
    !
    a_src%ifFieldGiven = .true.   ! switch between the cell list and field

    !$ if (omp_in_parallel()) call set_error("OMP unsafe",sub_name)
    !$ if (error) return

    ! A bit of preparations
    !
    spContent%sp => fu_work_string()
    indDescriptorOfParLine => fu_work_int_array()
    nullify(ptrItems)
    !
    ! Check that the source header has not yet been defined. Error otherwise.
    ! Note that field-based source cannot be redefined or added-up. All field files
    ! must be defined at once in one header.
    !
    if(a_src%defined == silja_true)then
      !
      ! The source has already been defined. 
      !
      call msg_warning('Area source v.4 cannot be defined twice',sub_name)
      call report(a_Src)
      call report(nlSrc)
      call set_error('Area source v.4 cannot be defined twice',sub_name)
      return
    endif
    !
    ! The source has only been initialised but nothing yet read. Set it from the namelist
    !
    a_src%src_nm = fu_Content(nlSrc, 'source_name')
    a_src%sector_nm = fu_Content(nlSrc, 'source_sector_name')
    !
    ! Do we get time variation from the binary files or it comes from header?
    !
    if(index(fu_str_u_case(fu_Content(nlSrc,'if_time_and_vertical_resolved_field_file')), 'YES') == 1)then

      ifTimeVertResolved = .true.
      
    elseif(index(fu_str_u_case(fu_Content(nlSrc,'if_time_and_vertical_resolved_field_file')), 'NO') == 1)then

      ifTimeVertResolved = .false.

    else
      call set_error('Strange if_time_and_vertical_resolved_field_file line:' + &
                   & fu_Content(nlSrc,'if_time_and_vertical_resolved_field_file'), sub_name)
      return
    endif

    vertical_unit = fu_content(nlSrc,'vertical_unit')
    release_rate_unit = fu_content(nlSrc,'release_rate_unit')

    iTmp = index(release_rate_unit, '/m2', back=.True.)
    if (iTmp > 0) then
        a_Src%ifFluxes = .True.
        if (iTmp /= len_trim(release_rate_unit)-2) then
          call msg("Strange source flux units. only /m2 suffix so far", iTmp, len_trim(release_rate_unit))
          call msg("release_rate_unit = "//trim(release_rate_unit))
          call msg("release_rate_unit = "//trim(release_rate_unit(1:(iTmp-1))))
          call set_error('Strange release_rate_unit', sub_name)
        endif
        release_rate_unit = release_rate_unit(1:(iTmp-1))
        quantity = emission_flux_flag
    else
        a_Src%ifFluxes = .False.
        quantity = emission_intensity_flag
    endif
    !
    ! Read the parameters from the header. Note that
    ! parameter strings are always needed: they describe the speciation and, if needed,
    ! contain scaling of the rate.
    !
    call get_items(nlSrc, 'par_str_area', ptrItems, nItems)
    if(error .or. nItems < 1)then
      call report(nlSrc)
      call set_error('No parameter items in namelist',sub_name)
      return
    endif
    allocate(a_src%params(nItems), stat=iTmp)
    if(fu_fails(iTmp == 0, 'Failed to allocate parameters for area source',sub_name))return
    call set_missing(a_src%params, nItems)

    indDescriptorOfParLine(1:a_src%nDescriptors) = int_missing

    do iTmp = 1, nItems
      call decode_time_parameter_v5(fu_content(ptrItems(iTmp)), &
                                  & iTmp, &
                                  & vertical_unit, &
                                  & release_rate_unit, & ! source release rate unit
                                  & a_src%params, &
                                  & a_src%cocktail_descr_lst, &    ! list of source cocktail descriptors
                                  & area_source, &   ! type of the source
                                  & a_src%nDescriptors, &
                                  & indDescriptorOfParLine, nDescriptorInParLine, &
                                  & fu_content(nlSrc,'total_emission_whole_period'), &
                                  & a_src%fDescrTotals, &      ! where to put decoded emission totals
                                  & a_src%ifSpecies)
      if(error)then
        call set_error('Failed source:' + a_src%src_nm + '_' + a_src%sector_nm, &
                     & sub_name)
        return
      endif
    end do   ! time slots
    !
    ! It is useful to know if the composition varies with time - just a speedup
    !
    SLOTS: do iTmp = 1, size(a_src%params)-1
      do iDescr = 1, a_src%nDescriptors
        a_src%ifRateVary = .not. (a_src%params(iTmp)%rate_descr_unit(iDescr) .eps. &
                                & a_src%params(iTmp+1)%rate_descr_unit(iDescr))
        if(a_src%ifRateVary)exit SLOTS
      end do
    end do SLOTS
    !
    ! Beware: if the binary file is time-resolving, ifRateVary must be .false.
    ! The reason is the times array: it will follow the binary time step, not the
    ! params list.
    !
    if(ifTimeVertResolved .and. a_src%ifRateVary)then
      call set_error('Varying time-slot rate is forbidden with time-resolving binary. Source:' + &
                   & a_src%src_nm + '_' + a_src%sector_nm, sub_name)
      return
    endif
    !
    ! Before opening the files, prepare the indexing: what file contains what descriptor
    !
    allocate(a_src%indFile4Descr(a_src%nDescriptors), stat=iTmp)
    if(fu_fails(iTmp==0,'Failed file4descriptor index allocation',sub_name))return
    a_src%indFile4Descr(1:a_src%nDescriptors) = int_missing

    !
    ! Variation coefficients may be needed if binary contains static field.
    !
    if(ifTimeVertResolved)then
      a_src%tz_index = 1
      a_src%tz_name = 'UTC'
    else
      !
      ! If time variation indices are present - read and set them, otherwise - set all to 1
      !
      call get_items(nlSrc,'hour_in_day_index',ptrItems,nItems)
      do iDescr = 1, nDescriptorInParLine
        a_src%indHour(indDescriptorOfParLine(iDescr),:) = 1.
        do iItem = 1, nItems
          spContent%sp = fu_content(ptrItems(iItem))
          read(unit=spContent%sp,fmt=*)chTmp
          if(chTmp == trim(fu_name(a_src%cocktail_descr_lst(indDescriptorOfParLine(iDescr))))) then
            read(unit=spContent%sp, iostat=iLocal, fmt=*) chTmp, &
                                    & (a_src%indHour(indDescriptorOfParLine(iDescr),iTmp),iTmp=1,24)
            if(iLocal /= 0)then
              call set_error('Failed to read  hourly variation:' + spContent%sp, &
                           & sub_name)
              return
            endif
            exit
          endif
        end do
      end do  ! iDescr

      call get_items(nlSrc,'day_in_week_index',ptrItems,nItems)
      do iDescr = 1, nDescriptorInParLine
        a_src%indDay(indDescriptorOfParLine(iDescr),:) = 1.
        do iItem = 1, nItems
          spContent%sp = fu_content(ptrItems(iItem))
          read(unit=spContent%sp,fmt=*)chTmp
          if(chTmp == trim(fu_name(a_src%cocktail_descr_lst(indDescriptorOfParLine(iDescr))))) then
            read(unit=spContent%sp, iostat=iLocal, fmt=*) chTmp, &
                                   & (a_src%indDay(indDescriptorOfParLine(iDescr),iTmp),iTmp=1,7)
            if(iLocal /= 0)then
              call set_error('Failed to read  daily-in-week variation:' + spContent%sp, &
                           & sub_name)
              return
            endif
            exit
          endif
        end do
      end do

      call get_items(nlSrc,'month_in_year_index',ptrItems,nItems)
      do iDescr = 1, nDescriptorInParLine
        a_src%indMonth(indDescriptorOfParLine(iDescr),:) = 1.
        do iItem = 1, nItems
          spContent%sp = fu_content(ptrItems(iItem))
          read(unit=spContent%sp,fmt=*)chTmp
          if(chTmp == trim(fu_name(a_src%cocktail_descr_lst(indDescriptorOfParLine(iDescr))))) then
            read(unit=spContent%sp, iostat=iLocal, fmt=*) chTmp, &
                                 & (a_src%indMonth(indDescriptorOfParLine(iDescr),iTmp),iTmp=1,12)
            if(iLocal /= 0)then
              call set_error('Failed to read month-in-year variation:' + spContent%sp, &
                           & sub_name)
              return
            endif
            exit
          endif
        end do
      end do

      do iDescr = 1, nDescriptorInParLine
        if(any(a_src%indHour(indDescriptorOfParLine(iDescr),1:24) < 0.999) .or. &
         & any(a_src%indHour(indDescriptorOfParLine(iDescr),1:24) > 1.0001) .or. &
         & any(a_src%indDay(indDescriptorOfParLine(iDescr),1:7) < 0.999) .or. &
         & any(a_src%indDay(indDescriptorOfParLine(iDescr),1:7) > 1.0001) .or. &
         & any(a_src%indMonth(indDescriptorOfParLine(iDescr),1:12) < 0.999) .or. &
         & any(a_src%indMonth(indDescriptorOfParLine(iDescr),1:12) > 1.0001))then
          a_src%ifUseTimeVarCoef = .true.
          exit
        endif
      enddo
      !
      ! ... and the vertical
      !
      ! Data file provides vertically-integrated flux
      ! read the multi-level structure that will be constant in time
      !
      spContent%sp = fu_str_u_case(fu_content(nlSrc,'vertical_distribution'))
      if(spContent%sp == 'SINGLE_LEVEL_DYNAMIC' .or. spContent%sp == '')then
        !
        ! Do nothing - all is in param_strings, which have been read above
        !
        call set_missing(a_src%vertLevs, .true.)
        nullify(a_src%levFraction)
        a_src%ifMultiLevelFixed = .false.

      elseif(spContent%sp == 'MULTI_LEVEL_FIXED')then
        !
        ! Read levels and set the source vertical
        !
        call get_items(nlSrc, 'vert_level', ptrItems, nItems)
        spContent%sp=fu_content(nlSrc,'vertical_unit')
        if(nItems < 1)then
          call set_error('No levels for vertical_distribution = MULTI_LEVEL_FIXED', &
                       & sub_name)
          return
        else
          allocate(a_src%levFraction(nItems))
          call set_named_level_with_fract(ptrItems(1), level, a_src%levFraction(1), spContent%sp)
          if(error)return
          call set_vertical(level, a_src%vertLevs)
          if(error)return
          do iTmp = 2, nItems
            call set_named_level_with_fract(ptrItems(iTmp), level, a_src%levFraction(iTmp), spContent%sp)
            call add_level(a_src%vertLevs, level)
            if(error)return
          end do
        endif
        if (abs(sum(a_src%levFraction(:)) - 1.) > 1e-3 ) then
          call msg("Fractions",a_src%levFraction(:))
           call report(a_src%vertLevs,.True.)
          call set_error("Strange level fractions", sub_name)
        endif

        a_src%ifMultiLevelFixed = .true.
      else
        call msg('vertical distribution must be MULTI_LEVEL_FIXED or SINGLE_LEVEL_DYNAMIC')
        call set_error('Strange vertical_distribution:' + spContent%sp, &
                     & sub_name)
        return
      endif
    
    
      a_src%tz_name = fu_content(nlSrc,'source_timezone')
      if( a_src%tz_name/= '')then

        if(fu_str_u_case(a_src%tz_name) == 'LOCAL')then
          ! Use local time at the source location
          a_src%tz_index = LocalTimeIndex
        elseif (fu_str_u_case(a_src%tz_name) == 'SOLAR')then
          a_src%tz_index = SolarTimeIndex
        else 
          a_src%tz_index =  timezone_index_by_name(a_src%tz_name)
          if (a_src%tz_index < 1) then
            call set_error("Failed to find Timezone index for '"//trim(a_src%tz_name)//"'", sub_name)
            return
          endif
        endif
      else 
        a_src%tz_name = 'Solar'
        a_src%tz_index = SolarTimeIndex
  !      call msg('No timezone for source. Assuming solar time.')
      endif
    endif ! ifTimeVertResolved  ! if vertical is taken from header

    !
    ! Having the basic information stored, get the binary files, get their headers
    ! and get the source grid. Should several files be used, compare the grids.
    ! Different data files are allowed for different species/cocktails. However, all 
    ! files MUST have identical temporal and spatial characteristics. 
    !
    call get_items(nlSrc,'field_emission_file',ptrItems,nItems)
    if(fu_fails(nItems >= 1, 'No field_emission_file lines found',sub_name))return
    allocate(a_src%chFIeldFNms(nItems), a_src%chFIeldFNmTemplates(nItems), &
           & a_src%FieldFormat(nItems), a_src%uField(nItems), stat=iLocal)
    if(fu_fails(iLocal == 0,'Failed binary files names allocation',sub_name))return
    nullify(a_src%times)

    prev_unit=int_missing
    do iItem = 1, nItems
      spContent%sp = adjustl(fu_content(ptrItems(iItem)))
      ! split the input at space:
      a_src%chFIeldFNms(iItem) = fu_process_filepath(spContent%sp(index(spContent%sp,' ')+1:), &  ! file name
                                                   & must_exist = .false., superdir = src_dir_name)
      a_src%FieldFormat(iItem) = fu_input_file_format(spContent%sp)                   ! format
      call decode_template_string(a_src%chFIeldFNms(iItem), a_src%chFIeldFNmTemplates(iItem))

      fname = fu_sample_file_name(a_src%chFIeldFNmTemplates(iItem), &
              & a_src%params(1)%time, a_src%params(size(a_src%params))%time, &
              &  a_src%src_nm, a_src%sector_nm)

      call open_input_file(fname, a_src%FieldFormat(iItem), a_src%uField(iItem), prev_unit) 
      if(error)return
      prev_unit = a_src%uField(iItem) !!Next file might be the same one -- no need to parse it again
      
      select case(a_src%FieldFormat(iItem)%iformat)
        case(grib_file_flag)
          call set_error('GRIB files are not supported',sub_name)

        case(ascii_file_flag)
          call set_error('ASCII files are not supported',sub_name)


        case(grads_file_flag)


          if(iItem == 1)then  
            ! get times and grid
            if(ifTimeVertResolved) &
                          & call get_grads_times(a_src%uField(iItem), a_src%times, nTimes, .true.)
            a_src%grid = fu_silamGrid_of_grads(a_src%uField(iItem))
            a_src%vertLevs = fu_silamVert_of_grads(a_src%uField(iItem))
            a_src%slotOffForFile = 0   ! the slot start is used for the file name generation (well, may be :-) )
            
          else                
            ! check consistency for times
            if(ifTimeVertResolved)then
              call get_grads_times(a_src%uField(iItem), timesTmp, nTimes, .true.)
              call check_times(timesTmp)
              if(associated(timesTmp))deallocate(timesTmp)
            endif
            !
            ! ... and grids
            allocate(pGrids(1),stat=iTmp)
            if(fu_fails(iTmp==0,'Failed temporary grid allocation',sub_name))return

            pGrids(1) = fu_silamGrid_of_grads(a_src%uField(iItem))
            if(.not. (a_src%grid == pGrids(1)))call set_error('new grid is not compatible with source', &
                                                            & sub_name)

            allocate(pVerts(1),stat=iTmp)
            if(fu_fails(iTmp==0,'Failed temporary vertical allocation',sub_name))return
            pVerts(1) = fu_silamVert_of_grads(a_src%uField(iItem))
            if(.not. fu_cmp_verts_eq(a_src%vertLevs, pVerts(1))) &
                                         & call set_error('new grid is not compatible with source', &
                                                        & sub_name)
            call set_missing(pVerts(1), defined(pVerts(1)))
            deallocate(pVerts)
          endif   ! iIitem == 1

!          if(ifTimeVertResolved) 
          call get_grads_IDs(a_src%uField(iItem), idList, nIDs)
         

        case(test_field_value_flag)

          if(iItem == 1)then
            ! Get times and grid
            if(ifTimeVertResolved)then
              allocate(a_src%times(2), stat=iLocal)
              if(iLocal /= 0)then
                call set_error('Failed to allocate times for the area source v4',sub_name)
                return
              endif
              a_src%times(1) = really_far_in_past
              a_src%times(2) = really_far_in_future
            endif
            a_src%grid = fu_set_grid(nlSrc)
            if(fu_fails(defined(a_src%grid),'test_field requires grid given in source namelist', &
                                           & sub_name))return
            call set_vertical(nlSrc, a_src%vertLevs)
            if(fu_fails(defined(a_src%vertLevs),'test_field requires vertical in source namelist', &
                                           & sub_name))return

!          else
!            ! check consistency for time and grid
!            if(ifTimeVertResolved)then
!              allocate(timesTmp(2), stat=iLocal)
!              if(iLocal /= 0)then
!                call set_error('Failed to allocate temporary times for the area source v4',sub_name)
!                return
!              endif
!              timesTmp(1) = really_far_in_past
!              timesTmp(2) = really_far_in_future
!              call check_times(timesTmp)
!              deallocate(timesTmp)
!            endif
!            if(defined(a_src%grid))then
!              call set_error('Test field is not compatible with pre-defined grid',sub_name)
!              return
!            endif
!            gridTmp = fu_set_grid(nlSrc)
!            if(.not. (a_src%grid == gridTmp))call set_error('new grid is not compatible with source', &
!                                                         & sub_name)
          endif  ! iItem == 1

          !
          ! check the list of IDs available from file and fill-in the corresponding descirptors
          !
          if(ifTimeVertResolved)then
            call decode_id_params_from_io_str(a_src%chFIeldFNms(iItem), &   ! test field string to decode
                                    & .false., &      ! ifMultiLevel
                                    & quantity, &   ! decoded quantity
                                    & speciesTmp, &  ! decoded substance name
                                    & .true.)   ! quantity name is really expected in the string
            if     (quantity == emission_intensity_flag) then
              a_src%ifFluxes = .True.
            elseif (quantity == emission_flux_flag) then
              a_src%ifFluxes = .False.
            else
              call set_error('Test emission field must be emission_intensity or emission_flux_flag, not:' + &
                                                                   & fu_quantity_string(iTmp), &
                      & sub_name)
              return
            endif
            allocate(idList(1),stat=iTmp)
            if(fu_fails(iTmp == 0, 'Failed idList allocation',sub_name))return
            nIDs = 1
            idList(1) = fu_set_field_id_simple(met_src_missing, quantity, &
                                             & time_missing, surface_level, chCockt = chTmp)
          endif  ! ifTimeVertResolved

        case(netcdf_file_flag)        !----------------------------------------------  NETCDF

 !         a_src%uField(iItem) = open_netcdf_file_i(a_src%chFIeldFNms(iItem), a_src%FieldFormat(iItem))
 !         if(error)return

          if(iItem == 1)then
            !
            ! Define times
            if(ifTimeVertResolved) then 
              call timelst_from_netcdf_file(a_src%uField(iItem), a_src%times, iTmp, .true., .true.)
              if (fu_fails(.not. error, 'Error after timelst_from_netcdf_file', sub_name)) return
              !! Check if we should use start or end of slot to address the file
              fnameTmp = fu_FNm(a_src%chFIeldFNmTemplates(iItem), &
                           & a_src%times(1), a_src%times(1), zero_interval, &
                           & chSource = a_src%src_nm, chSector = a_src%sector_nm)
              if (fnameTmp == fname) then 
                a_src%slotOffForFile = 0
              else
                !Fname matches end of last slot
                fnameTmp = fu_FNm(a_src%chFIeldFNmTemplates(iItem), &
                           & a_src%times(iTmp), a_src%times(iTmp), zero_interval, &
                           & chSource = a_src%src_nm, chSector = a_src%sector_nm)
                if (fnameTmp == fname) then
                  a_src%slotOffForFile = 1
                else
                  call report(a_src%times)
                  call msg("Template: "//trim(fu_template_string(a_src%chFIeldFNmTemplates(ind_file))) )
                  call set_error("Filename does not match nether beginning nor end of the envelope", sub_name)
                  return
                endif
              endif
            endif
            !
            ! .. .and grids
            call get_netcdf_grids(a_src%uField(iItem), &
                                & quantity, species_missing, pGrids, num_grids)

            if(fu_fails(num_grids >= 1, "No grids for area source",sub_name))return

            do ind_grid = 2, num_grids
              if (fu_fails(pGrids(ind_grid) == pGrids(1), 'Netcdf has several different grids',&
                         & sub_name)) return
            end do
            a_src%grid = pGrids(1)
            deallocate(pGrids); nullify(pGrids)
            !
            ! ... and verticals if resolved
            if (ifTimeVertResolved) then
              call get_netcdf_verticals(a_src%uField(iItem), &
                                      & quantity, species_missing, pVerts, num_verts)
              do ind_vert = 2, num_verts
                if (fu_fails(fu_cmp_verts_eq(pVerts(ind_vert), pVerts(1)), &
                           & 'Netcdf has several different verticals',&
                           & sub_name)) return
              end do
              a_src%vertLevs = pVerts(1)
              deallocate(pVerts); nullify(pVerts)
            end if

          else
            !
            ! Check times
            if(ifTimeVertResolved)then
              call timelst_from_netcdf_file(a_src%uField(iItem), timesTmp, iTmp, .true.,.true.)
              call check_times(timesTmp)
              if(associated(timesTmp))deallocate(timesTmp)
            endif
            !
            ! ... and grids
            call get_netcdf_grids(a_src%uField(iItem), &
                                & quantity, species_missing, pGrids, iTmp)
            if(fu_fails(iTmp == 1,'Strange number of grids in netcdf:' + fu_str(iTmp),&
                                & sub_name))return

            if(fu_fails(a_src%grid == pGrids(1),'new grid is not compatible with source', &
                                              & sub_name))return
            deallocate(pGrids); nullify(pGrids)
            !
            ! ... and verticals
            call get_netcdf_verticals(a_src%uField(iItem), &
                                    & quantity, species_missing, pVerts, iTmp)
            if(fu_fails(iTmp == 1,'Strange number of verticals in netcdf:' + fu_str(iTmp),&
                                & sub_name))return

            if(fu_fails(fu_cmp_verts_eq(a_src%vertLevs, pVerts(1)), &
                      & 'new vertical is not compatible with source',sub_name))then
                call report(a_src%vertLevs)
                do iTmp = 1, fu_NbrOfLevels(a_src%vertLevs)
                    call report(fu_level(a_src%vertLevs,iTmp))
                enddo
                call report(pVerts(1)) 
                do iTmp = 1, fu_NbrOfLevels(a_src%vertLevs)
                    call report(fu_level(pVerts(1),iTmp))
                enddo
                return
            endif
            deallocate(pVerts); nullify(pVerts)
          endif   ! iItem == 1
          !
          ! For the future use - get the list of IDs available
          !
          !if(ifTimeVertResolved) call id_list_from_netcdf_file(a_src%uField(iItem), idList, nIDs, .true.)
          call id_list_from_netcdf_file(a_src%uField(iItem), idList, nIDs)
        case default
          call set_error('Unknown file format:'+fu_str(a_src%FieldFormat(iItem)%iformat), &
                       & sub_name)
      end select  ! file types
      if(error)then
        call set_error('Failed to process the line:' + spContent%sp, sub_name)
        return
      endif
      !
      ! check the list of IDs available from file and fill-in the corresponding descirptors
      !
!      if(ifTimeVertResolved)then
        call fill_file_indices(idList, nIDs, iItem)
        if(error)return
!      endif
    enddo  ! list of data files

    ! Check that all the par-line descriptors are found in the binary
    !
    do iDescr = 1, a_src%nDescriptors
      if(a_src%indFile4Descr(iDescr) == int_missing)then
        call msg('This descriptor is not found in the binary:',iDescr)
        call report(a_src%cocktail_descr_lst(iDescr))
        call msg('Area source consumed this-far:')
        call report(a_src)
        call set_error('Descriptors in source ini and binaries do not meet',sub_name)
        return
      endif
    end do
 
    !
    ! The way to store the fields depends on whether the time- and vertical-related
    ! info is in the binary file or in the source ini.
    ! 
    if(.not.  ifTimeVertResolved)then  !!! resloved fields allocated in project_a_src_second_grd
      !
      ! If the field is neither vertically- nor time-resolved, it can be stored
      ! into the cells array and then used as the structures of source v.3
      !
      call get_work_arrays_set(a_src%nDescriptors, fu_number_of_gridpoints(a_src%grid), pArTmp)
      if(error)return
      
      if (fu_fails(nIds > 0, 'No field ids in list', sub_name)) return
      call msg_warning('Non-time-resolving emission - will take the first valid time from input')
            
      do iDescr = 1, a_src%nDescriptors
        !
        ! Get the field. A trick: different species can easily be in different binary files
        ! So, we scan them all for each descrptor. But it is forbidden to have one
        ! descriptor existing in more than one binary.
        !

        ! Get the field, use to indfile4descr mapping. Do NOT touch the binaries: they are opened above
        !
        ind_file = a_src%indFile4Descr(iDescr)
        fname = fu_FNm(a_src%chFIeldFNmTemplates(ind_file), &
                     & fu_analysis_time(idList(1)), fu_analysis_time(idList(1)), &
                     & fu_forecast_length(idList(1)), &  ! Input file
                     & chSource = a_src%src_nm, chSector = a_src%sector_nm)
        idTmp1 =  fu_set_field_id(met_src_missing,&
                                & quantity, &
                                & fu_analysis_time(idList(1)), & !analysis_time,&
                                & fu_forecast_length(idList(1)), & !forecast_length, &
                                & a_src%grid, &
                                & fu_level(idList(1)), &
                                & chCocktail=fu_name(a_src%cocktail_descr_lst(iDescr)))

        call get_input_field(fname, &
                           & a_src%FieldFormat(ind_file), & ! Input file format
                           & idTmp1, & !!Requested ID
                           & pArTmp(iDescr)%pp, &
                           & a_src%grid, &          ! storage grid
                           & .false., setZero, & !ifAcceptSameMonth, iOutside, 
                           & idTmp, &  ! output
                           & iAccuracy=5, wdr=wdr_missing, &
                           & iBinary = a_src%uField(ind_file))

#ifdef DEBUG                           
       if (any(pArTmp(iDescr)%pp < 0)) then
           call ooops ("Negative field")
           pData => pArTmp(iDescr)%pp
           call grid_dimensions(a_src%grid, nx, ny)
           do iTmp=1,nx*ny
              if (pArTmp(iDescr)%pp(iTmp)<0) call ooops("Here!")
           enddo
       endif
#endif           

      enddo  ! iDescr
      !
      ! All fields are consumed. Transform them into cells, and we're back in the
      ! classical area source.
      !
      call fields_to_cells(pArTmp)  !!!!This guy converts fluxes to rates if needed

      ! The times need to be set: use the par_str lines from the header file.
      ! 
      allocate(a_src%times(size(a_src%params)), stat=itmp)
      if (fu_fails(itmp == 0, 'Allocate failed (times, non-resolving)', sub_name)) return
      do itmp = 1, size(a_src%params)
        a_src%times(itmp) = a_src%params(itmp)%time
      end do
      
      ! Now we should close the file: there can be many hundreds of them, netcdf lib can run
      ! out of buffer space. Since the binary is no longer needed, can close it
      !
      do iItem = 1, nItems
!        call msg('Closing the binary file:' + a_src%chFIeldFNms(iItem))
        call close_input_file(a_src%FieldFormat(iItem)%iformat, a_src%uField(iItem))
      end do
      a_src%ifFieldGiven = .false. ! Now it is the same standard as v3

      call free_work_array(pArTmp)

    endif  ! ifTimeResolved .or. ifVerticalResolved
    !
    ! Things have been done but do NOT destroy the namelist - it was NOT you who made it, stupid!
    !
    call free_work_array(spContent%sp)
    call free_work_array(indDescriptorOfParLine)
    call destroy_items(ptrItems)
    !
    ! We do not yet know the dispersion grid and vertical, so interpolation cannot be set
    !
    a_src%ifHorizInterp = silja_undefined
    a_src%ifVertInterp = silja_undefined
    nullify(a_src%pHorizIS)  ! horizontal interpolation structure
    nullify(a_src%pVertIS)    ! vertical interpolation strcture

    a_src%defined = silja_true


  CONTAINS

    !=============================================================================

    subroutine fields_to_cells(field_list)
      implicit none
      type(silja_rp_1d), dimension(:), intent(in) :: field_list
      
      integer :: ix, iy, nx, ny, field_size, count_nonzero, iDescr, count_values_set, ind_cell, i1d
      real :: fTmp
      integer, dimension(:), pointer :: is_nonzero
      real, dimension(:), pointer :: cellArea
      real, dimension(:), pointer :: field_for_descr 
     
      
      ! Find out how many cells have emission for at least some descriptor.
      !
      cellArea => null()
      call grid_dimensions(a_src%grid, nx, ny)
      field_size = nx*ny
      is_nonzero => fu_work_int_array(field_size)
      is_nonzero(1:field_size) = 0
      do iDescr = 1, a_src%nDescriptors
        field_for_descr => field_list(iDescr)%pp
        do i1d = 1,field_size
          if (field_for_descr(i1d) > 0) is_nonzero(i1d) = 1
        end do
      end do

      if (a_src%ifFluxes) then
        cellArea => fu_work_array(field_size)
        do i1d = 1,field_size
           if (is_nonzero(i1d) > 0) then
             cellArea(i1d) = fu_cell_size(a_src%grid, i1d)
           endif 
        enddo
      endif
      count_nonzero = sum(is_nonzero(1:field_size))
      
      if (count_nonzero == 0) call msg_warning('Area source is empty')

      a_src%nCells = count_nonzero
      call msg ("Allocating a_src%nCells of total grid size",a_src%nCells, field_size)
      nullify(a_src%cell_fx);    call set_array_size(a_src%cell_fx, a_src%nCells, real_missing)
      nullify(a_src%cell_fy);    call set_array_size(a_src%cell_fy, a_src%nCells, real_missing)
      nullify(a_src%cell_val);  call set_array_size(a_src%cell_val, a_src%nDescriptors, a_src%nCells, 0.0)

      if (error) then
          call set_error("Trouble allocating asrc", "fields_to_cells")
          call report(a_src)
          return
      endif

      ! Fill the cell array advancing as indicated by is_nonzero.
      !
      count_values_set = 0
      do iDescr = 1, a_src%nDescriptors
        ind_cell = 0
        field_for_descr => field_list(iDescr)%pp
        do iy = 1, ny
          do ix = 1, nx
            i1d = (iy-1)*nx + ix
            if (is_nonzero(i1d) > 0) then
              ind_cell = ind_cell + 1
              if (idescr == 1) count_values_set = count_values_set + 1
            else
              cycle
            end if
            !if (field_for_descr(i1d) > 0) then
            if (fu_fails(ind_cell <= count_nonzero, 'ind_cell too big', 'fields_to_cells')) return
            
            if (a_src%ifFluxes) then   
              a_src%cell_val(iDescr,ind_cell) = field_for_descr(i1d)*cellArea(i1d)
            else
              a_src%cell_val(iDescr,ind_cell) = field_for_descr(i1d)
            endif
            a_src%cell_fx(ind_cell) = real(ix)
            a_src%cell_fy(ind_cell) = real(iy)
            !end if
          end do
        end do
      end do

      if (fu_fails(count_values_set == count_nonzero, 'Counts mismatch', 'fields_to_cells')) return
      call free_work_array(is_nonzero)
      if (a_src%ifFluxes) then
        call free_work_array(cellArea)
      endif

      if (fu_fails(.not.any(a_src%cell_fx .eps. real_missing), 'real_missing in fx', 'fields_to_cells')) return
    end subroutine fields_to_cells

    logical function fu_if_all_zero(iItemIn)
      implicit none
      integer, intent(in) :: iItemIn
      integer :: iDescrInt
      fu_if_all_zero = .true.
      do iDescrInt = 1, a_src%nDescriptors
        if(pArTmp(iDescrInt)%pp(iItemIn) /= 0.0)then
          fu_if_all_zero = .false.
          return
        endif
      end do
    end function fu_if_all_zero


    !=============================================================================
    
    subroutine check_times(arTimes)
      !
      ! Checks that the given array of times is identical to that written in a_src
      !
      implicit none
      type(silja_time), dimension(:), pointer :: arTimes
      !
      ! Local variables
      logical :: ifOK
      if(associated(arTimes))then
        ifOK = associated(a_src%times)
        if(ifOK) ifOK = (size(a_src%times) == size(arTimes))
        if(ifOK) then
          do iTmp = 1, size(arTimes)
            ifOK = (arTimes(iTmp) == a_src%times(iTmp))
            if(fu_fails(ifOK, 'Source time array is not equal to the new one','check_times'))return
          end do
        endif
        !
        ! Additional trouble: if params have different times from the one here, we will be in trouble.
        ! A way out is to require exact correspondence
        !
        if(a_src%params(1)%time < a_src%times(1) .or. &
         & a_src%params(size(a_src%params))%time > a_src%times(size(a_src%times)))then
          call msg('Start in par_str=' + fu_str(a_src%params(1)%time) + &
                 & ', first source slot time=' + fu_str(a_src%times(1)))
          call msg('Last in par_str=' + fu_str(a_src%params(size(a_src%params))%time) + &
                 & ', first source slot time=' + fu_str(a_src%times(size(a_src%times))))
          call set_error('Inconsistent start and end of the source parameter strings','check_times')
        endif
      else
        if(associated(a_src%times))call set_error('arTimes is associated but a_src%times is not','check_times')
      endif
    end subroutine check_times
    
    !============================================================================
    
    subroutine fill_file_indices(idList, nIDs, iFile)
      !
      ! Determines what descriptor is available from the given file and stores its index
      ! Also overrides the species with one extracted from the file 
      !
      implicit none
      
      ! Imported parameters
      type(silja_field_id), dimension(:), pointer :: idList
      integer, intent(in) :: nIDs, iFile
      
      ! Local variables
      integer :: iDescr, iID, quantity
      character(len = substNmLen) :: chIDDescrNm
      type(silam_species) :: speciesTmp
      type(Tcocktail_descr) :: descrTmp
      
      !
      ! Just scan the IDs checking for quantity and descriptor names. The rest is unimportant.
      !
      quantity = emission_intensity_flag
      if (a_src%ifFluxes) quantity = emission_flux_flag
      do iID = 1, nIDs
        if(fu_quantity(idList(iID)) == quantity)then
          chIDDescrNm = fu_cocktail_name(idList(iID))
          if(len_trim(chIDDescrNm) == 0)then   
            !
            ! there can be still species given, try to make the descriptor out of it
            !
            call msg_warning('cocktail_name is not given in the fieldID, trying species','fill_file_indices')
            speciesTmp = fu_species(idList(iID))
            if(defined(speciesTmp))then
              call set_descriptor_from_species(speciesTmp, descrTmp, fu_str(speciesTmp))
              chIDDescrNm = fu_name(descrTmp)
            else
              call msg_warning('Cannot get emission cocktail name from ID','fill_file_indices')
              call report(idList(iID))
              call set_error('Cannot get emission cocktail name from ID','fill_file_indices')
              return
            endif    ! if descriptor name can be made from species
          endif   ! descriptor name is absent
          !
          ! Find the descriptor in the a_src list of descriptors. Careful: if this is actually substance
          ! or species, we have to compare them, not the name of the fake cocktail
          !
          do iDescr = 1, a_src%nDescriptors
            if(fu_str_u_case(chIDDescrNm) == &
             & fu_str_u_case(fu_name(a_src%cocktail_descr_lst(iDescr))))then
              if(a_src%indFile4Descr(iDescr) == int_missing .or. &
               & a_src%indFile4Descr(iDescr) == iFile)then
                a_src%indFile4Descr(iDescr) = iFile
                if (defined(speciesTmp) ) then 
                  speciesArPtr => fu_species(a_src%cocktail_descr_lst(iDescr))
                  speciesArPtr(1) =  speciesTmp  !!!Force actual species from ID
                endif
              else
                call msg('Descriptor: ' + fu_name(a_src%cocktail_descr_lst(iDescr)) + &
                       & '-is available from two files:', a_src%indFile4Descr(iDescr), iFile)
                call set_error('Descriptor: ' + fu_name(a_src%cocktail_descr_lst(iDescr)) + &
                       & '-is available from two files', 'fill_file_indices@fill_a_src_from_namelist_v4')
              endif
            endif
          end do
        end if
      end do
      
    end subroutine fill_file_indices

  end subroutine fill_a_src_from_namelist_v4


  ! ***************************************************************


  subroutine store_a_src_as_namelist(a_src, uOut, ifOriginalGrd)
    !
    ! Writes one area source term v.3 to an external file. Can do it either in 
    ! original or in dispersion grid and vertical. Will expand all val values
    ! using the emission cocktail descriptors, etc.
    !
    ! All units: SI, unless otherwise is stated
    !
    IMPLICIT NONE
    !
    ! Imported parameters
    TYPE(silam_area_source), intent(in) :: a_src
    integer, intent(in) :: uOut
    logical, intent(in) :: ifOriginalGrd

    ! Local variables
    integer :: iTmp, iDescr

    ! Basic parameters: Source name, sector and grid
    !
    if(a_src%ifFieldGiven)then
      write(uOut,fmt='(A)')'AREA_SOURCE_4'
    else
      write(uOut,fmt='(A)')'AREA_SOURCE_3'
    endif

    write(uOut,fmt='(2A)')'source_name = ', trim(a_src%src_nm)
    write(uOut,fmt='(2A)')'source_sector_name = ', trim(a_src%sector_nm)
    write(uOut,fmt='(2A)')'source_timezone = ', trim(a_src%tz_name)
    call report_as_namelist(a_src%grid, uOut)
    !
    ! Release rate parameters
    ! We will write the whole source in one file, the whole cocktail at once
    !
    write(uOut,fmt='(A)')'release_rate_unit = DESCRIPTOR_DEFAULT'

    !
    ! Time slot and variation parameters
    !
    do iTmp = 1, size(a_src%params)
      call store_time_param_as_namelist(a_src%params(iTmp), a_src%cocktail_descr_lst, uOut, area_source)
    end do
    do iDescr = 1, a_src%nDescriptors
      write(uOut,fmt='(2A,1x,24(F10.6,1x))')'hour_in_day_index = ', &
                               & trim(fu_name(a_src%cocktail_descr_lst(iDescr))),&
                               & (a_src%indHour(iDescr,iTmp),iTmp=1,24)
      write(uOut,fmt='(2A,1x,7(F10.6,1x))')'day_in_week_index = ', &
                               & trim(fu_name(a_src%cocktail_descr_lst(iDescr))),&
                               & (a_src%indDay(iDescr,iTmp),iTmp=1,7)
      write(uOut,fmt='(2A,1x,12(F10.6,1x))')'month_in_year_index = ', &
                               & trim(fu_name(a_src%cocktail_descr_lst(iDescr))),&
                               & (a_src%indMonth(iDescr,iTmp),iTmp=1,12)
    end do
    !
    ! If  multi-level structure is constant in time, write it
    !
    if(a_src%ifMultiLevelFixed)then
      write(uOut,fmt='(A)')'vertical_distribution = MULTI_LEVEL_FIXED'
      if(ifOriginalGrd)then
        call report_as_namelist(uOut, a_src%vertLevs, a_src%levFraction, .true.) ! if skip 0
      else
        call report_as_namelist(uOut, a_src%vertLevsDispVert, a_src%levFractDispVert, .true.)
      end if  ! original grid
    else
      write(uOut,fmt='(A)')'vertical_distribution = SINGLE_LEVEL_DYNAMIC'
    endif  ! if multi level time-fixed emission
    !
    ! If this is the field-based source, there is no reason to store the data
    !
    if(a_src%ifFieldGiven)then

      write(uOut,fmt='(A)')'if_time_and_vertical_resolved_field_file = YES'

      do iTmp = 1, size(a_src%FieldFormat)
        select case(a_src%FieldFormat(iTmp)%iformat)
          case(grib_file_flag)
            write(uOut, fmt='(A)')'field_emission_file = GRIB ', a_src%chFIeldFNms(iTmp)
          case(grads_file_flag)
            write(uOut, fmt='(A)')'field_emission_file = GRADS ', a_src%chFIeldFNms(iTmp)
          case(ascii_file_flag)
            write(uOut, fmt='(A)')'field_emission_file = ASCII ', a_src%chFIeldFNms(iTmp)
          case(netcdf_file_flag)
            write(uOut, fmt='(A)')'field_emission_file = NETCDF ', a_src%chFIeldFNms(iTmp)
          case(test_field_value_flag)
            write(uOut, fmt='(A)')'field_emission_file = TEST_FIELD ', a_src%chFIeldFNms(iTmp)
          case default
            call set_error('Unknown file format:'+fu_str(a_src%FieldFormat(iTmp)%iformat), &
                         & 'store_a_src_as_namelist')
            return
        end select
      end do

    else
      !
      ! Values will be in the above grid
      ! The val array is 2D: (nSubst, nModes)
      !
      write(uOut,fmt='(A)')'coordinate_of_values = GRID_INDICES'

      if(ifOriginalGrd)then
        do iTmp = 1, a_src%nCells
          write(uOut,fmt='(A,100(F15.7,1x))')'val = ', &
                                        & a_src%cell_fx(iTmp), &
                                        & a_src%cell_fy(iTmp), &
                                        & (a_src%cell_val(iDescr,iTmp), iDescr=1, a_src%nDescriptors)
        end do
      else
        do iTmp = 1, a_src%nCellsDispGrd
          write(uOut,fmt='(A,100(F15.7,1x))')'val = ', &
                                       & a_src%cellDispGrd_fx(iTmp), &
                                       & a_src%cellDispGrd_fy(iTmp), &
                                       & (a_src%cellDispGrd_val(iDescr,iTmp), &
                                                              & iDescr=1,a_src%nDescriptors)
        end do
      endif
    endif
    if(a_src%ifFieldGiven)then
      write(uOut,fmt='(A)')'END_AREA_SOURCE_4'
    else
      write(uOut,fmt='(A)')'END_AREA_SOURCE_3'
    endif
    write(uOut,fmt=*)

  end subroutine store_a_src_as_namelist


  !**************************************************************************
  
  subroutine write_area_src_from_mass_map(mapEmis, timeStart, timeEnd, chCase, chSrcNm, &
                                        & pMapPx, pMapPy, pMapPz, gTemplate)
    !
    ! Writes a header of binary-containing area source
    !
    implicit none
    
    ! Imported parameter
    type(Tmass_map), pointer :: mapEmis, pMapPx, pMapPy, pMapPz
    type(silja_time), intent(in) :: timeStart, timeEnd
    character(len=*), intent(in) :: chCase, chSrcNm
    type(grads_template), intent(in) :: gTemplate

    ! Local variables
    integer :: uOut, iTmp
    real, dimension(max_levels) :: fFracs !Dummy level fractions
    character (len=fnlen) :: outName
    character (len=*), parameter :: sub_name= "write_area_src_from_mass_map"

    !$ if (omp_in_parallel()) call set_error("OMP unsafe",sub_name)
    !$ if (error) return

    uOut = fu_next_free_unit()
    if(error)return

    ! File name comes from the pure-GrADS template
    !
    outName = fu_FNm(gTemplate, timeStart, timeStart, timeEnd - timeStart, chCase, chSrcNm) + &
             & '_ems_dump.sa4'
    if(error)return

    open(uOut,file=outName,iostat=iTmp)
    if(fu_fails(iTmp==0,'Failed to open area-source header file',sub_name))return



    if(associated(pMapPx) .or. associated(pMapPy) .or. associated(pMapPz)) &
        & call msg_warning('Cannot make use of moments yet',sub_name)

    ! Basic parameters: Source name, sector and grid
    !
    write(uOut,fmt='(A)')'AREA_SOURCE_4'

    write(uOut,fmt='(A)')'source_name = massMap_source'
    write(uOut,fmt='(A)')'source_sector_name = massMap_source'
    call report_as_namelist(mapEmis%gridTemplate, uOut)
    !
    ! Release rate parameters
    ! We will write the whole source in one file, the whole cocktail at once
    !
    write(uOut,fmt='(A)')'release_rate_unit = DESCRIPTOR_DEFAULT'
    !
    ! Time slot and variation parameters
    call store_time_param_as_namelist(timeStart, 10.0, 0.0, mapEmis%species, mapEmis%nSpecies, &
                                    & uOut, area_source)
    call store_time_param_as_namelist(timeEnd, 10.0, 0.0, mapEmis%species, mapEmis%nSpecies, &
                                    & uOut, area_source)
    !
    ! If  multi-level structure is constant in time, write it
    !
    fFracs(:) = 1.0
    write(uOut,fmt='(A)')'vertical_distribution = MULTI_LEVEL_FIXED'
    call report_as_namelist(uOut, mapEmis%vertTemplate, fFracs, .true.) ! if skip 0

    !
    ! If this is the field-based source, there is no reason to store the data
    !
    write(uOut,fmt='(A)')'if_time_and_vertical_resolved_field_file = YES'

    write(uOut, fmt='(A,1x,A)')'field_emission_file = GRADS ', &
                    & trim(fu_FNm(gTemplate, timeEnd, timeEnd, zero_interval, chCase, chSrcNm) + &
                    & '_ems_dump.grads.super_ctl')
    write(uOut,fmt='(A)')'END_AREA_SOURCE_4'
    write(uOut,fmt=*)
    close(uOut)

  end subroutine write_area_src_from_mass_map


  !**************************************************************************

  subroutine check_time_params_a_src(a_src, timestart, timestep)
    !
    ! If the source is v4 and has time variations active, model time step must be shorter 
    ! than an hour, and exact hour must be the edge of the time steps.
    !
    implicit none

    ! Improted parameters
    type(silam_area_source), intent(in) :: a_src
    type(silja_time), intent(in) :: timestart
    type(silja_interval), intent(in) :: timestep
    
    ! Local variab;es
    integer :: iDescr

    if(a_src%ifUseTimeVarCoef .and. a_src%ifFieldGiven)then
      !
      ! Long time step?
      !
      if(fu_abs(timestep) > one_hour .and. any(abs(a_src%indHour(:,:) - 1.0) > 1e-5))then
        do iDescr = 1, a_src%nDescriptors
          call msg('Hourly coefs for descriptor=' + fu_str(iDescr), a_src%indHour(iDescr,:))
        end do
        call set_error('Model timestep:' + fu_str(timestep) + ', hourly time variation', &
                     & 'check_time_params_a_src')
        return
      endif
      !
      ! Weird run start? We need every hour to be at the edge of the time steps. With
      ! timestep shorter than an hour, it is "almost" enough to check that the start of the 
      ! hour is met exactly from the current starttime using the given timestep. To make it 
      ! absolutely certain, check two sequential hours
      !
      if(.not. fu_next_special_time(timestart, iStartOfHour, forwards, zero_interval) == &
             & fu_next_special_time(timestart, iStartOfHour, forwards, timestep) .or. &
       & abs(real(nint(one_hour / timestep)) - (one_hour / timestep)) > 1e-5 )then
        call set_error('Timestep:' + fu_str(timestep) + &
                & ', is not divisor of one_hour or run steps miss hour edges due to wrong start:' + &
                & fu_str(timestart),'check_time_params_a_src')
        return
      endif
    endif  !  binary and time variatrion

  end subroutine check_time_params_a_src


  !**************************************************************************

  subroutine getTimeSlots_of_area_source(a_src, time, iSlot1, iSlot2)
    !
    ! Finds the time slots surrounding the given time
    ! 
    IMPLICIT NONE

    ! Imported parameters
    TYPE(silam_area_source), INTENT(in) :: a_src
    type(silja_time), intent(in) :: time
    integer, intent(out) :: iSlot1, iSlot2

    call getTimeSlots_from_params(a_src%params, time, iSlot1, iSlot2)

  END subroutine  getTimeSlots_of_area_source


  !**************************************************************************

  subroutine link_a_src_to_species(species_list, a_src)
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
    type(silam_area_source), intent(inout) :: a_src
    type(silam_species), dimension(:), pointer :: species_list

    call link_src_species_to_given_list(a_src%nDescriptors, &
                                      & a_src%nSpeciesInDescr, &
                                      & a_src%cocktail_descr_lst, &
                                      & a_src%fDescr2SpeciesUnit, &
                                      & a_src%pEmisSpeciesMapping, &
                                      & species_list)

  end subroutine link_a_src_to_species


  !*******************************************************************

  subroutine add_source_species_a_src(a_src, speciesLst, nSpecies)
    !
    ! Fills-in the given list with the own species. Checks for the duplicates
    !
    implicit none

    ! Improted parameters
    type(silam_area_source), intent(in) :: a_src
    type(silam_species), dimension(:), pointer :: speciesLst
    integer, intent(inout) :: nSpecies

    ! Local variables
    integer :: iDescr, nSpeciesDescr
    type(silam_species), dimension(:), pointer :: pSpecies
    pSpecies =>null()
    !
    ! If the source species list is not ready - compile it
    !
    do iDescr = 1, a_src%nDescriptors
      call get_inventory(a_src%cocktail_descr_lst(iDescr), pSpecies, nSpeciesDescr)
      if(error)return
      call addSpecies(speciesLst, nSpecies, pSpecies, nSpeciesDescr)
      if(error)return
    end do

  end subroutine add_source_species_a_src


  !*****************************************************************

  integer function fu_source_id_nbr_of_a_src(a_src)
    !
    ! Returns the source number. Reason: all sources are enumerated
    ! sequencially, so that the source can, in fact, be refered by its
    ! number without other information. But so far this number is 
    ! copied to the particles in the pollution cloud. 
    !
    implicit none

    ! Imported parameters with intent IN
    type(silam_area_source), intent(in) :: a_src

    ! Stupidity check
    if(.not.defined(a_src%grid))then
      call set_error('Undefined source given','fu_source_nbr_of_area_source')
      return
    endif
    fu_source_id_nbr_of_a_src = a_src%id_nbr

  end function fu_source_id_nbr_of_a_src


  !*******************************************************************

  subroutine get_lp_release_height_a_src (a_src, iSlot, fWeightPast, z, ifPressure)
    !
    ! Checks the vertical structure of the point source and returns the randomised
    ! vertical position with unit determined by the source definitions - either 
    ! pressure in Pa or height in metres
    !
    implicit none

    ! Imported parameters
    type(silam_area_source), intent(in) :: a_src
    integer, intent(in) :: iSlot
    real, intent(in) :: fWeightPast
    real, intent(out) :: z
    logical, intent(out) :: ifPressure  ! indicator of unit of z

    ! Local variables
    type(silja_level) :: lyr
    real :: fTmp, fTmpSum
    integer :: iTmp

    !
    ! The key question: is the vertical fixed and multi-layer or single-layer and dynamic?
    ! The whole story depends on this
    !
    if(a_src%ifMultiLevelFixed)then
      !
      ! Multi-level fixed layers. Select first the layer by random choice and then
      ! the position inside the layer
      !
      fTmp = fu_random_number_boundaries(0.0,sum(a_src%levFraction(1:size(a_src%levFraction))))
      fTmpSum = 0.
      lyr = level_missing
      do iTmp = 1, size(a_src%levFraction)
        fTmpSum = fTmpSum + a_src%levFraction(iTmp)
        if(fTmp <= fTmpSum)then
          lyr = fu_level(a_src%vertLevs,iTmp)
          exit
        endif
      end do

      if(.not. defined(lyr))then
        call set_error('Could not select the level','get_lp_release_height_a_src')
        call report(a_src)
        return
      endif

    else
      !
      ! Dynamic vertical. Despite strong efficiency requirements, have to check at least basics
      !
      if (iSlot < 1 .or. iSlot > size(a_src%params) .or. &
       & (iSlot == size(a_src%params) .and. fWeightPast < 0.9999))then
        call msg('Strange slot requested',iSlot)
        call set_error('Strange slot requested','get_lp_release_height_a_src')
        return
      endif

      if(fWeightPast < 0.9999)then
        lyr = a_src%params(iSlot)%layerDynamic * fWeightPast + &
            & a_src%params(iSlot+1)%layerDynamic * (1.-fWeightPast)
      else
        lyr = a_src%params(iSlot)%layerDynamic * fWeightPast
      endif
      if(error)then
        call msg('Problem interponlating the slot layers. iSlot, fWeight:',iSlot,fWeightPast)
        call report(a_src)
        return
      endif

    endif ! multi- or singlw layer

    !--------------------------------------------------
    !
    ! Set initial particle conditions with grid transformation.
    !
    SELECT CASE (fu_leveltype(lyr))
      CASE (layer_btw_2_pressure)

        ifPressure = .true.
        SELECT CASE (a_src%params(iSlot)%vertical_distribution)
          CASE (vertically_single_level)
            z = fu_pr_level_pressure(fu_lower_boundary_of_layer(lyr) + fu_upper_boundary_of_layer(lyr)) * 0.5
          CASE (vertically_even_distribution)
            z = fu_random_number_boundaries(fu_pr_level_pressure(fu_upper_boundary_of_layer(lyr)), &
                                          & fu_pr_level_pressure(fu_lower_boundary_of_layer(lyr)))
          CASE default
            CALL set_error('Can not handle this vertical distribution',&
                         & 'get_lp_release_height_a_src')
            RETURN
        END SELECT
        IF (error) RETURN

      CASE (layer_btw_2_height)

        ifPressure = .false.
        SELECT CASE (a_src%params(iSlot)%vertical_distribution)
          CASE (vertically_single_level)
            z = fu_level_height(fu_lower_boundary_of_layer(lyr) + fu_upper_boundary_of_layer(lyr)) * 0.5
          CASE (vertically_even_distribution)
            z = fu_random_number_boundaries(fu_level_height(fu_upper_boundary_of_layer(lyr)), &
                                          & fu_level_height(fu_lower_boundary_of_layer(lyr)))
          CASE default
            CALL set_error('Can not handle this vertical distribution',&
                         & 'get_lp_release_height_a_src')
            RETURN
        END SELECT
        IF (error) RETURN

      CASE default
        CALL set_error('only constant_pressure and constant_height allowed',&
                     & 'get_lp_release_height_a_src')
        RETURN
    END SELECT  ! type of vertical - pressure or height

  end subroutine get_lp_release_height_a_src


  !*****************************************************************

  subroutine total_from_a_src_descr_unit(a_src, ifSlotRateOnly, &     ! mandatory, input
                                       & amounts, &                   ! mandatory, output
                                       & start_, duration_, layer_)   ! optional, input
    !
    ! Returns the amount of the released material IN DESCRIPTOR UNIT starting from 
    ! start during the duration time interval. Descriptors must be initialised by that moment
    !
    ! We have to integrate the release rate, which is given in descriptor basic units - for each
    ! descriptor. 
    !
    ! Four cases are considered (notebook 11 pp.21-23)
    ! 1. slot rates are constant: descriptor rate and vertical fraction
    ! 2. cell-connected time variation coefficients are constant (fHour, fDay, fMonth)
    ! 3. both sets are varying
    ! 4. None of them is time dependent
    !
    ! This subroutine produces a total descriptor-wise emission and returns a vector.
    ! If requested, it can skip the cell-wise spatial integration leaving the space for
    ! further map- or total- aiming integration over cells.
    !
    ! Units: SI basic and derived ones
    !
    implicit none

    ! Imported parameters with intent IN
    type(silam_area_source), intent(in) :: a_src
    logical, intent(in) :: ifSlotRateOnly  ! if cell-wise integration is not neeed, just par_str
    real, dimension(:), intent(out) :: amounts
    type(silja_time), intent(in), optional :: start_
    type(silja_interval), intent(in), optional :: duration_
    type(silja_level), intent(in), optional :: layer_  ! if defined, include only fraction emitted in it

    ! Local variables
    integer :: iSlot, iStartSlot, iEndSlot, iTmp, iDescr, iCell
    real :: fTimeConversion, fTmp, fWeightPast, fVertFraction1, fVertFraction2, fTime_sec, &
          & fTime_sec_start, fTime_sec_end, fSlot_sec
    type(silja_time) :: start, start_integr, end_integr, startTmp, endTmp, timeTmp
    type(silja_interval) :: duration, timeStep
    logical :: ifSlotRateVary, ifFound
    real, dimension(max_descriptors_in_source)  :: fCells
    type(silja_level) :: layer
    type(silja_time) :: StartOfWeek
    integer(4) :: iMyOffsetSec, iSecondsSinceStartOfWeek, iMon, iWeekDay, iHour
!    integer :: iCellX, iCellY
    type(silja_rp_1d), dimension(:), pointer :: pArTmp
    pArTmp =>null()
    !
    ! Time period for integration, then start and end slots. Use times array
    !
    if(present(start_))then
      start = start_
      if(start > a_src%times(size(a_src%times)))then
#ifdef DEBUG                           
        call report(a_src)
        call msg_warning('Source ends before the computation starts','total_from_a_src_descr_unit')
#endif        
        return
      endif
      do iStartSlot = 1, size(a_src%times)-1
        if(a_src%times(iStartSlot+1) > start)exit ! Skip whole slots before the start time
      end do
    else
      start = a_src%times(1)
      iStartSlot = 1
    endif
    if(present(duration_))then
      duration = duration_
      if(start+duration > a_src%times(size(a_src%times))) &
                                      & duration = a_src%times(size(a_src%times)) - start
      if(start+duration < a_src%times(1))then
#ifdef DEBUG                           
        call report(a_src)
        call msg_warning('Source starts after the computation ends','total_from_a_src_descr_unit')
#endif        
        return
      endif
      do iEndSlot = iStartSlot+1, size(a_src%times)
        if(a_src%times(iEndSlot-1) < start+duration)exit ! Skip whole slots after the end time
      end do
      if(iEndSlot > size(a_src%times)) iEndSlot = size(a_src%times)
    else
      iEndSlot = size(a_src%times)
      duration = a_src%times(iEndSlot) - start
    endif
    if(present(layer_))then
      layer = layer_
    else
      layer = level_missing
    endif

    amounts(1:a_src%nDescriptors) = 0.
    !
    ! To have non-zero released amount, there must be at least two slots in the given range
    !
    if(iStartSlot >= iEndSlot)then
      return
    endif


    !--------------------------------------------------------------------------------
    !
    ! If the data are given in the binary file as field, which is by-definition
    ! time- and vertically- resolving, the integration will go over there
    !
    if(a_src%ifFieldGiven)then
      !
      ! The data are in the time- and vertical-resolving datafile. If total 
      ! is given, use it, otherwise call for internal integration procedures
      !
      ifFound = .false.
      do iDescr = 1, a_src%nDescriptors
        if(a_src%fDescrTotals(iDescr) .eps. real_missing)then
          ifFound = .true.
          exit
        endif
      end do
!call msg('11')
      if(ifFound)then 
        !
        ! At least one descriptor has missing total, have to sum-up the binaries
        ! First, prepare the strucutre
        ! 
!call msg('12')
        call get_work_arrays_set(a_src%nDescriptors, pArTmp)
        if(error)return
        do iDescr = 1, a_src%nDescriptors
          pArTmp(iDescr)%pp(1:a_src%nCells) = 0.0
        end do
!call msg('13')
        !
        ! Now use this structure to get maps of time-wise totals over the layer 
        ! for all descriptors from the binaries
        !
        do iSlot = 1, size(a_src%FieldFormat)
          select case(a_src%FieldFormat(iSlot)%iformat)  ! GrADS/NetCDF/..., ref to input structure
            case(grads_file_flag)
              call get_grads_total(a_src%uField(iSlot), a_src%pFldEmis, &
                                 & start, duration, layer, &
                                 & pArTmp)
              if(error)return
            case(netcdf_file_flag)
              call get_netcdf_total(a_src%uField(iSlot), a_src%pFldEmis, &
                                 & start, duration, layer, &
                                 & pArTmp)
              if(error)return
            case(ascii_file_flag) 
              call set_error('ASCII field format is not supported as time-resolving', &
                           & 'total_from_a_src_descr_unit')
              return
            case default
              call set_error('Field format is not supported:' + fu_str(a_src%FieldFormat(iSlot)%iformat), &
                           & 'total_from_a_src_descr_unit')
              return
          end select

          do iDescr = 1, a_src%nDescriptors
            amounts(iDescr) = amounts(iDescr) + sum(pArTmp(iDescr)%pp(1:a_src%nCells)) * &
                                              & a_src%params(1)%rate_descr_unit(iDescr)
          end do

        enddo    ! binary formats
        call free_work_array(pArTmp)

      else
        !
        ! all descriptor totals are reasonable
        !
        amounts(1:a_src%nDescriptors) = amounts(:) + a_src%fDescrTotals(:) * fu_sec(duration) / &
                                      & fu_sec(a_src%times(size(a_src%times)) - a_src%times(1))
      endif  ! if real_missing total found

      return  ! all done

    endif  ! ifFieldGiven


    !--------------------------------------------------------------------------------
    !
    ! The source is given as a list of cells. 
    ! Start the slot integration
    !
    do iSlot = iStartSlot, iEndSlot-1
      !
      ! Check the vertical layer - may be, we do not have to do anything - if there is no emission
      ! into this layer ?
      ! Note that the layer can be dynamic or static.
      !
      if(defined(layer))then
        if(a_src%ifMultiLevelFixed)then 
          !
          ! fixed vertical, decided by vertLevs
          fVertFraction1 = 0.
          do iTmp = 1, fu_NbrOfLevels(a_src%vertLevs)
            fVertFraction1 = fVertFraction1 + a_src%levFraction(iTmp) * &
                                    & fu_vert_overlap_fraction(fu_level(a_src%vertLevs,iTmp), layer)
          end do
          if(fVertFraction1 < 1e-5) return  ! No emission into this output layer
          fVertFraction2 = fVertFraction1
        else                
          !
          ! dynamic vertical, decided by params array
          !
          fVertFraction1 = fu_vert_overlap_fraction(a_src%params(iSlot)%layerDynamic, layer)
          fVertFraction2 = fu_vert_overlap_fraction(a_src%params(iSlot+1)%layerDynamic, layer)
          if((fVertFraction2 < 1e-5 ) .and. (fVertFraction1 < 1e-5)) cycle
        endif
      else
        fVertFraction1 = 1.0
        fVertFraction2 = 1.0
      endif

      !
      ! Compute the time overlap
      !
      if(a_src%times(iSlot) < start)then ! Skip part of slot if needed
        start_integr = start
      else
        start_integr = a_src%times(iSlot)
      endif
      if(a_src%times(iSlot+1) > start+duration)then ! Skip part of the slot if needed
        end_integr = start+duration
      else
        end_integr = a_src%times(iSlot+1)
      endif

      ! Might be not needed, but do it once....
      StartOfWeek = fu_Start_Of_Day_utc(start_integr) - fu_set_interval_d(real(fu_weekday(start_integr)-1))
      !   call msg("Total from area source Starttime:"+fu_str(start_integr))
      !   call msg("Total from area source    StartOfWeek:"+fu_str(StartOfWeek))

      if (a_src%tz_index > 0) then ! One timezone for source
            iMyOffsetSec =  TZ_offset(a_src%tz_index)
      else
            iMyOffsetSec = 0
      endif
      !
      ! Check the rate / composition variation
      !
      ifSlotRateVary = a_src%ifRateVary .or. .not. (fVertFraction1 .eps. fVertFraction2)



      !--------------------------------------------------------------------------------
      !
      ! Here we can finally decide on the integration procedure:
      ! One of the four procedures outlined above. See notebook 11, p.24.
      !
      if(ifSlotRateVary .and. a_src%ifUseTimeVarCoef)then  !-------------------- all dynamic
        !
        ! The most-difficult case. Everything varies, have to go cell-by-cell and hour-by-hour
        ! (or day-by-day, if hourly variation is void). Case 3 of the above function description
        !
        if(ifSlotRateOnly)then
          !
          ! Cells are not involved, then we have a sub for the slot integration
          !
          fCells(1:a_src%nDescriptors) = 1.0
          call integrate_dynamic_slot_descr(a_src%params, &    
                                          & start_integr, end_integr, iSlot, &
                                          & a_src%cocktail_descr_lst, a_src%nDescriptors, &
                                          & fCells, &
                                          & fVertFraction1, fVertFraction2, &
                                          & a_src%ifRateVary, &
                                          & amounts)
                                  
        else
          !
          ! The main and the most-complicated part of the story. Everything varies and all is included.
          ! Some minor simplification might still be possible but probably is not worth the trouble.
          ! Go step-by-step and cell-by-cell
          !
          fSlot_sec = fu_sec(a_src%params(iSlot+1)%time - a_src%params(iSlot)%time)
          fTime_sec_start = fu_sec(start_integr - a_src%params(iSlot)%time)
          fWeightPast = 1. - fTime_sec_start/fSlot_sec
          !
          ! First, cover the start period until the first round hour, then go ahead with hourly step
          !
          startTmp = start_integr
          endTmp = fu_round_closest_hour(startTmp, forwards)
          if(endTmp > end_integr) endTmp = end_integr

          fTime_sec = fu_sec(endTmp - startTmp)

          ! In some cases, start is the round hour, so can skip this step
          !
          do while(startTmp < end_integr)

            fWeightPast = fWeightPast + 0.5*fTime_sec/fSlot_sec


            iSecondsSinceStartOfWeek = nint(fu_sec(startTmp-StartOfWeek),8)
            iMon = fu_mon(startTmp) ! Ignore localtime for month
                                        ! one size fits all
!call msg("Injectarea seconds since StartOfWeek:"+ fu_str(iSecondsSinceStartOfWeek))
!call msg("Injectarea seconds since StartOfWeek:"+ fu_str(iMon))
            do iDescr = 1, a_src%nDescriptors
              fCells(iDescr) = 0.0
              do iCell = 1, a_src%nCellsDispGrd
                ! If cell-wise index -- take it
                if (a_src%tz_index == solarTimeIndex .or. a_src%tz_index == LocalTimeIndex) then
                  iMyOffsetSec = TZ_offset(a_src%CellTZindex(iCell))
                endif
                iHour = modulo(iSecondsSinceStartOfWeek + iMyOffsetSec, i_seconds_in_day) / i_seconds_in_hour + 1
                iWeekDay = modulo(iSecondsSinceStartOfWeek + iMyOffsetSec, i_seconds_in_week)/ i_seconds_in_day + 1
                fCells(iDescr) = fCells(iDescr) + a_src%cellDispGrd_Val(iDescr,iCell) * &
                       & a_src%indHour(iDescr, iHour) * &
                       & a_src%indDay (iDescr, iWeekDay) * &
                       & a_src%indMonth(iDescr,iMon )
              end do
              amounts(iDescr) = amounts(iDescr) + &
                                   & fTime_sec * &                ! time period
                                   & (fWeightPast*fVertFraction1 + (1-fWeightPast)*fVertFraction2) * & !vert
                                   & (fWeightPast * a_src%params(iSlot)%rate_descr_unit(iDescr) + &    ! descr rate
                                       & (1.-fWeightPast) * a_src%params(iSlot+1)%rate_descr_unit(iDescr)) * &
                                   & fCells(iDescr)                       ! sum over cells with time variation incl.
            end do  ! iDescr

            !
            ! The time step is completed. The next one is one hour or shorter if this is the end of integration
            !
            fWeightPast = fWeightPast + 0.5*fTime_sec/fSlot_sec
            startTmp = endTmp

            fTime_sec = 3600.
            endTmp = startTmp + one_hour
            if(endTmp > end_integr)then
              endTmp = end_integr
              fTime_sec = fu_sec(endTmp - startTmp)
            endif

          end do  ! startTmp < end_integr

        endif  ! if Slot Rate only


      elseif(ifSlotRateVary)then  !-------------------------------------------- time variation coef void
        !
        ! Hourly/daily/monthly time variation coefficients =1 but slot rate varies.
        ! Can integrate analytically over the whole slot, then multiply with cell-specific
        ! constant rates (case 2 in the above function description). See notebook 11, p.24.
        !
        if(ifSlotRateOnly)then
          !
          ! Cells are not involved, then we have a sub for the slot integration
          !
          fCells(1:a_src%nDescriptors) = 1.0

        else
          !
          ! Full-rate is requested (still, no time variation coefficients), cells involved. 
          ! A slot-long integration possible
          !
          do iDescr = 1, a_src%nDescriptors
            fCells(iDescr) = 0.0
            do iCell = 1, a_src%nCellsDispGrd
              fCells(iDescr) = fCells(iDescr) + a_src%cellDispGrd_Val(iDescr,iCell)
            end do
          end do

        endif  ! whether the slot only rate is requestde

        call integrate_dynamic_slot_descr(a_src%params, &    
                                        & start_integr, end_integr, iSlot, &
                                        & a_src%cocktail_descr_lst, a_src%nDescriptors, &
                                        & fCells, &
                                        & fVertFraction1, fVertFraction2, &
                                        & a_src%ifRateVary, &
                                        & amounts)
        if(error)return

      elseif(a_src%ifUseTimeVarCoef)then !----------------------------------------- slot rates const
        !
        ! Slot-wise rates are constant but time variation coefficients are non-unity. Case 1 in
        ! the function description. Go hour by hour or day by day
        !
        if(ifSlotRateOnly)then        ! no cells here means no time variation left
          fTime_sec = fu_sec(end_integr - start_integr) 
          do iDescr = 1, a_src%nDescriptors
            amounts(iDescr) = amounts(iDescr) + &
                                & fTime_sec * &            ! time period 
                                & fVertFraction1 * &       ! vertical overlap
                                & a_src%params(iSlot)%rate_descr_unit(iDescr) ! descr.rate
          end do

        else                          ! include cells, i.e. include time variation
          !
          ! Cells are involved. Have to go hour-by-hour. However, speedup is possible because full-day
          ! average of hourly variation is unity. Thus, if the duration is long enough, can count 
          ! days rather than hours - thanks to constant main slot rates. 
          !
          !
          ! First, cover the period until the first round hour, then go with hourly or daily step
          !
          startTmp = start_integr
          endTmp = fu_round_closest_hour(startTmp, forwards)
          if(endTmp > end_integr) endTmp = end_integr

          ! In some cases, the hourly variation is void. Then daily step is possible
          !
          timeStep = one_day
          do iDescr = 1, a_src%nDescriptors
            do iTmp = 1, 24
              if(.not. (a_src%indHour(iDescr,iTmp) .eps. 1.0))then
                timeStep = one_hour
                exit
              endif
            end do
          end do


          fTime_sec = fu_sec(endTmp - startTmp)
          fWeightPast = 0.0
          fSlot_sec = fu_sec(a_src%params(iSlot+1)%time - a_src%params(iSlot)%time)
          ! In some cases, start is the round hour, so can skip this step
          !
          do while(startTmp < end_integr)

            fWeightPast = fWeightPast + 0.5*fTime_sec/fSlot_sec
            iSecondsSinceStartOfWeek = nint(fu_sec(startTmp-StartOfWeek),4)
            iMon = fu_mon(startTmp) ! Ignore localtime for month

            do iDescr = 1, a_src%nDescriptors
              fCells(iDescr) = 0.0
              do iCell = 1, a_src%nCellsDispGrd
                  ! If cell-wise index -- take it
                  if (a_src%tz_index == solarTimeIndex .or. a_src%tz_index == LocalTimeIndex) then
                    iMyOffsetSec = TZ_offset(a_src%CellTZindex(iCell))
                  endif
                  iHour = modulo(iSecondsSinceStartOfWeek + iMyOffsetSec, i_seconds_in_day) / i_seconds_in_hour + 1
                  iWeekDay = modulo(iSecondsSinceStartOfWeek + iMyOffsetSec, i_seconds_in_week)/ i_seconds_in_day + 1

                  fCells(iDescr) = fCells(iDescr) + a_src%cellDispGrd_Val(iDescr,iCell) * &
                                                & a_src%indHour(iDescr, iHour) * &
                                                & a_src%indDay (iDescr, iWeekDay) * &
                                                & a_src%indMonth(iDescr,iMon )
              end do
              amounts(iDescr) = amounts(iDescr) + &
                                   & fTime_sec * &                ! time period
                                   & fVertFraction1 * &           ! vertical overlap
                                   & a_src%params(iSlot)%rate_descr_unit(iDescr) * & ! descr.rate
                                   & fCells(iDescr)                ! sum over cells with time variation incl.
            end do  ! iDescr

            !
            ! The time step is completed. The next one is one hour or
            ! shorter if this is the end of integration
            !
            fWeightPast = fWeightPast + 0.5*fTime_sec/fSlot_sec
            startTmp = endTmp

            fTime_sec = fu_sec(timeStep)
            endTmp = startTmp + timeStep
            if(endTmp > end_integr)then
              endTmp = end_integr
              fTime_sec = fu_sec(endTmp - startTmp)
            endif

          end do  ! startTmp < end_integr

        endif  ! if Slot Rate only

      else  !----------------------------------------------------------------------- all static
        !
        ! All is static, i.e. time-independent. The fastest option (case 4 above)
        !
        fTime_sec = fu_sec(end_integr - start_integr)

        if(ifSlotRateOnly)then        ! no cells
          do iDescr = 1, a_src%nDescriptors
            amounts(iDescr) = amounts(iDescr) + &
                               & fTime_sec * &            ! time period
                               & fVertFraction1 * &       ! vertical overlap
                               & a_src%params(iSlot)%rate_descr_unit(iDescr) ! descr.rate
          end do
        else                          ! include cells
          do iDescr = 1, a_src%nDescriptors
            fCells(iDescr) = 0.0
            do iCell = 1, a_src%nCellsDispGrd
              fCells(iDescr) = fCells(iDescr) + a_src%cellDispGrd_Val(iDescr,iCell)
            end do
            amounts(iDescr) = amounts(iDescr) + &
                                 & fTime_sec * &             ! time period
                                 & a_src%params(iSlot)%rate_descr_unit(iDescr) * & ! descr rate, descr unit
                                 & fVertFraction1 * &        ! vertical overlap
                                 & fCells(iDescr)            ! sum over cells for this descriptor
          end do
        endif
      endif  ! four types of temporal dependence

    end do  ! time slots

!do iDescr = 1, a_src%nDescriptors
!  call msg("Amount from descriotpr "+fu_str(iDescr)+":", amounts(iDescr))
!end do

  end subroutine total_from_a_src_descr_unit


  !*****************************************************************

  subroutine total_amt_species_unit_a_src(a_src, ifSlotRateOnly, &      ! mandatory, input
                                        & species, nSpecies, amounts, & ! mandatory, output
                                        & start_, duration_, layer_)   ! optional, input
    !
    ! Returns the amount of the released material IN SPECIES UNIT starting from 
    ! start during the duration time interval. Species must be initialised by that moment
    !
    ! We have to integrate the release rate, which is given in descriptor basic units - for each
    ! descriptor. Species may have different units, so the return value is an array of species
    !
    ! Four cases are considered (notebook 11 pp.21-23)
    ! 1. slot rates are constant: descriptor rate, species fractions, and vertical fraction
    ! 2. cell-connected time variation coefficients are constant (fHour, fDay, fMonth)
    ! 3. both sets are varying
    ! 4. None of them is time dependent
    !
    ! This subroutine produces a total species-wise emission and returns a vector.
    ! If requested, it can skip the cell-wise spatial integration leaving the space for
    ! further map- or total- aiming integration over cells.
    !
    ! Units: SI basic and derived ones
    !
    implicit none

    ! Imported parameters with intent IN
    type(silam_area_source), intent(in) :: a_src
    logical, intent(in) :: ifSlotRateOnly  ! if cell-wise integration is not neeed, just par_str
    type(silam_species), dimension(:), pointer :: species
    integer, intent(out) :: nSpecies
    real, dimension(:), pointer :: amounts
    type(silja_time), intent(in), optional :: start_
    type(silja_interval), intent(in), optional :: duration_
    type(silja_level), intent(in), optional :: layer_  ! if defined, include only fraction emitted in it

    ! Local variables
    integer :: iDescr, iSpecies, nSpeciesDescr
    type(silam_species), dimension(:), pointer :: speciesDescr
    type(silja_time) :: start
    type(silja_interval) :: duration
    real, dimension(max_descriptors_in_source) :: amountsDescr
    type(silja_level) :: layer
    type(chemical_adaptor) :: adaptor
    !
    ! Time period for integration, then start and end slots
    !
    speciesDescr =>null()

!call msg('1')
    
    nSpecies = 0

    if(present(start_))then
      start = start_
#ifdef DEBUG                                      
      if(start > a_src%params(size(a_src%params))%time)then
        call report(a_src)
        call msg_warning('Source ends before the computation starts','total_amt_species_unit_a_src')
        !return
      endif
#endif        
    else
      start = a_src%params(1)%time
    endif
    if(present(duration_))then
      duration = duration_
      if(start+duration > a_src%params(size(a_src%params))%time) &
                                      & duration = a_src%params(size(a_src%params))%time - start
#ifdef DEBUG                                      
      if(start+duration < a_src%params(1)%time)then
        call report(a_src)
        call msg_warning('Source starts after the computation ends','total_amt_species_unit_a_src')
        !return
      endif
#endif        
    else
      duration = a_src%params(size(a_src%params))%time - start
    endif
    if(present(layer_))then
      layer = layer_
    else
      layer = level_missing
    endif

    !
    ! Start from obtaining the amounts in the descriptor units
    !
    amountsDescr(:) = 0.
    !
    ! ... and the amount themselves
    !
    call total_from_a_src_descr_unit(a_src, ifSlotRateOnly, &
                                   & amountsDescr, &
                                   & start, duration, layer)
    if(error)return
!do iDescr = 1, a_src%nDescriptors
!call msg('Emission for descriptor: ' + fu_name(a_src%cocktail_descr_lst(iDescr)) + ':',amountsDescr(iDescr))
!end do

    !if(sum(amountsDescr(1:a_src%nDescriptors)) == 0.)return

!call msg('3')

    !
    ! Conversion into the species units is comparatively straightforward.
    ! First make the list of species emitted from the source - for the output
    !
    nullify(species)
    nSpecies = 0
    call add_source_species_a_src(a_src, species, nSpecies)
    amounts(1:nSpecies) = 0.
    !
    ! Now explore the descriptors
    !
    do iDescr = 1, a_src%nDescriptors
      call get_inventory(a_src%cocktail_descr_lst(iDescr), speciesDescr, nSpeciesDescr)
      call create_adaptor(speciesDescr, species, adaptor)
      if (error) return
      do iSpecies = 1, nSpeciesDescr
         
!call msg('Factor to species unit for:'+fu_name(a_src%cocktail_descr_lst(iDescr)), a_src%fDescr2SpeciesUnit(iSpecies,iDescr))
        amounts(adaptor%iSp(iSpecies)) = amounts(adaptor%iSp(iSpecies)) + &
                                & amountsDescr(iDescr) * a_src%fDescr2SpeciesUnit(iSpecies,iDescr)
!call msg('Amounts for species:' + fu_name(fu_material(species(adaptor%iSp(iSpecies)))),amounts(adaptor%iSp(iSpecies)))
      end do
    end do

  end subroutine total_amt_species_unit_a_src


  !****************************************************************

  subroutine get_cell_grid_coords_a_src(a_Src, iCell, fCellX, fCellY)
    !
    ! Returns the cell coords, as they are written in the source term.
    ! No grid conversion.
    ! Depending on whether the data are given in field or list of cells,
    ! iCell means different things: 1-D index in the grid or index in the 
    ! list of cells
    !
    implicit none

    ! Imported parameters
    type(silam_area_source), intent(in) :: a_Src
    integer, intent(in) :: iCell
    real, intent(out) :: fCellX, fCellY
    
    integer :: nx, ny, iTmp

    if(a_src%ifFieldGiven)then
      if(iCell <= 0)then
        call set_error('Negative source cell number:'+fu_str(iCell), 'get_cell_coords_a_src')
        return
      endif
      call grid_dimensions(a_src%grid, nx, ny)
      if(iCell > nx*ny)then
        call set_error('Strange source cell number:'+fu_str(iCell), 'get_cell_coords_a_src')
      else
        iTmp = mod(iCell,ny)
        if(iTmp == 0) iTmp = nx
        fCellX = real(iTmp)
        fCellY = real((iCell - iTmp) / (ny-1))
      endif
    else
      if(iCell <=0 .or. iCell > a_src%nCells)then
        call set_error('Strange source cell number:'+fu_str(iCell),'get_cell_coords_a_src')
      else
        fCellX = a_src%cell_fx(iCell)
        fCellY = a_src%cell_fy(iCell)
      endif
    endif  ! ifFieldGiven
  end subroutine get_cell_grid_coords_a_src


  !****************************************************************

  subroutine get_cell_geo_coords_a_src(a_Src, iCell, fCellX, fCellY)
    !
    ! Returns the cell coords, as they are written in the source term.
    ! No grid conversion
    !
    implicit none

    ! Imported parameters
    type(silam_area_source), intent(in) :: a_Src
    integer, intent(in) :: iCell
    real, intent(out) :: fCellX, fCellY
    
    ! local variables
    real :: fX, fY

    if(a_src%ifFieldGiven)then
      call get_cell_grid_coords_a_src(a_src,iCell,fX, fY)
      if(error)return
      fCellX = fu_lon_native_from_grid(fX, fY, a_src%grid)
      fCellY = fu_lat_native_from_grid(fX, fY, a_src%grid)
    else
      if(iCell <=0 .or. iCell > a_src%nCells)then
        call set_error('Strange source cell number','get_cell_coords_a_src')
      else
        fCellX = fu_lon_native_from_grid(a_src%cell_fx(iCell), a_src%cell_fy(iCell), a_src%grid)
        fCellY = fu_lat_native_from_grid(a_src%cell_fx(iCell), a_src%cell_fy(iCell), a_src%grid)
      endif
    endif
  end subroutine get_cell_geo_coords_a_src


  !****************************************************************

  real function fu_total_cell_value_a_Src(a_src, iCell)
    !
    ! Returns the cell value, as they are written in the source term.
    ! No unit conversion, time variation, etc.
    !
    implicit none

    ! Imported parameters
    type(silam_area_source), intent(in) :: a_Src
    integer, intent(in) :: iCell

    if(a_src%ifFieldGiven)then
      call set_error('Does not work yet','fu_total_cell_value_a_Src')
      return
    else
      if(iCell <= 0 .or. iCell > a_src%nCells)then
        call set_error('Strange source cell number','fu_total_cell_value_a_Src')
      else
        fu_total_cell_value_a_Src = sum(a_src%cell_Val(1:a_src%nDescriptors,iCell))
      endif
    endif
  end function fu_total_cell_value_a_Src

  !****************************************************************

  function fu_cell_values_a_Src(a_src, iCell) result(valPtr)
    !
    ! Returns the cell value, as they are written in the source term.
    ! No unit conversion, time variation, etc.
    !
    implicit none

    ! return pointer
    real, dimension(:), pointer :: valPtr

    ! Imported parameters
    type(silam_area_source), intent(in) :: a_Src
    integer, intent(in) :: iCell

    if(a_src%ifFieldGiven)then
      call set_error('Does not work yet','fu_cell_values_a_Src')
      return
    else
      if(iCell <= 0 .or. iCell > a_src%nCells)then
        call set_error('Strange source cell number','fu_cell_values_a_Src')
      else
        valPtr => a_src%cell_Val(:,iCell)
      endif
    endif
  end function fu_cell_values_a_Src


  !****************************************************************

  subroutine create_src_contain_grd_a_src(a_src, grid_template, ifVerbose, ifMinimal, ifExtended)
    !
    ! This routine creates the grid similar to the given grid_template
    ! covering all sources of the em_source
    !
    implicit none

    ! Imported parameters
    type(silam_area_source), intent(in) :: a_src
    type(silja_grid), intent(inout) :: grid_template
    logical, intent(in) :: ifVerbose, ifMinimal
    logical, intent(out) :: ifExtended

    ! Local variables
    integer :: nx, ny, iCell, iType, iMin, iMax, jMin, jMax, ix_start, ix_end, iy_start, iy_end
    real :: x, y, sw_corner_modlon, sw_corner_modlat, dx, dy, pole_x, pole_y
    integer, dimension(:,:), pointer :: arFlag 
    integer, dimension(:), pointer :: iWork 
    character(len=clen) :: chName
    logical :: if_south_pole, if_corner_in_geo_coord



    !$ if (omp_in_parallel()) call set_error("OMP unsafe",'create_src_contain_grd_a_src')
    !$ if (error) return

    if(.not. defined(a_src))then
      call set_error('Undefiend source','create_src_contain_grd_a_src')
      call report(a_src)
      return
    endif

    !
    ! The template has a right to be undefined. Then we just accept the source grid
    ! as the template one but cut it to include this source only
    !
    if(defined(grid_template))then
      !
      ! If defined template, one can try to cut it down to the source area
      !
      if(ifMinimal)then
        ix_start = 1
        iy_start = 1
        call grid_dimensions(grid_template, ix_end, iy_end)
        if(error)return

        iWork => fu_work_int_array(ix_end*iy_end)
        if(error)return
        arFlag(ix_start:ix_end, iy_start:iy_end) => iWork(1:ix_end*iy_end)
        arFlag(:,:) = -1
        !
        ! Cut non-overlapping parts of the grid_template
        !
        call fill_overlap_flag_array(grid_template, &        ! grid to cut
                                   & a_src%grid, &           ! grid area to preserve
                                   & arFlag, &               ! overlap flag array
                                   & cover_the_grid_area, 2) ! filling algorithm, overlap mark
        if(error)then
          call free_work_array(iWork)
          call unset_error('create_src_contain_grd_a_src')  !Why unset????
          return
        endif

        call cut_empty_lines(arFlag, ix_start, iy_start, ix_end, iy_end, &
                           & -1, 1, 2)  !value_to_cut, value_to_stay, value_to_preserve)
        if(error)return

        call cut_grid_size(grid_template, ix_start, iy_start, ix_end, iy_end)

        call free_work_array(iWork)
      endif  ! ifMinimal

    else
      !
      ! Grid_template is undefined: create it using this source grid as a starting point
      ! but cut it down to only existing cells
      ! If gridtype is anygrid, will create some latlon grid 
      !
      if(a_src%nCells < 1)return

      grid_template = a_src%grid

      if(.not. a_src%ifFieldGiven)then 
        !
        ! chacking cells if they exist
        !
        iMin = int(a_src%cell_fx(1)+0.5)
        iMax = iMin
        jMin = int(a_src%cell_fy(1)+0.5)
        jMax = jMin
        do iCell = 2, a_src%nCells
          if(iMin > int(a_src%cell_fx(iCell)+0.5)) iMin = int(a_src%cell_fx(iCell)+0.5)
          if(iMax < int(a_src%cell_fx(iCell)+0.5)) iMin = int(a_src%cell_fx(iCell)+0.5)
          if(jMin > int(a_src%cell_fy(iCell)+0.5)) jMin = int(a_src%cell_fy(iCell)+0.5)
          if(jMax < int(a_src%cell_fy(iCell)+0.5)) jMax = int(a_src%cell_fy(iCell)+0.5)
        end do
        call cut_grid_size(grid_template, iMin, iMax, jMin, jMax)
        if(error)return
      endif  ! if not field in binary
    endif  ! if grid_template defined

    if(a_src%ifFieldGiven)return ! if field given, declared source grid is returned: cannot check all binary

    !
    ! Scan its central points and extend the template_grid if these points are outside. 
    !
    call grid_dimensions(grid_template, nx,ny)

    ifExtended = .false.
    do iCell = 1, a_src%nCells
      call project_point_to_grid(a_src%grid, &
                               & a_src%cell_fx(iCell), &
                               & a_src%cell_fy(iCell), &
                               & grid_template, x, y)
      if(x<1 .or. x>nx .or. y<1 .or. y>ny) then
        call extend_grid_to_coordinates(grid_template, x, y)
        ifExtended = .true.
      endif
    end do

    if(ifVerbose .and. ifExtended)then
      call msg(fu_connect_strings('Extending the grid for the area source:',fu_name(a_src)))
    endif

    if(.not.fu_stdSilamGrid(grid_template))then
      !
      ! Manipulations resulted in non-standard grid. It either jumps over 180/-180 meridian or has
      ! more than 360deg in any direction. However, it can be real if source is global or
      ! itself crosses 180-th meridian. So far, our concern is global source.
      !
      call msg('Area surce:'+ a_src%src_nm + '_' + a_src%sector_nm + &
             & '. Below source-covering grid is not standard. Will try to adjust it')
      call report(grid_template)
      if(fu_ifLonGlobal(a_src%grid))then
        call make_grid_lon_global(grid_template, .true.)  ! adjust resolution if needed
        if(error)return
        !
        ! Note that only lon-global grids can be lat-global
        !
        call msg('The source-cevering grid after lon-global correction:')
        call report(grid_template)
        if(fu_ifLatGlobal(a_src%grid))then
          if(fu_ifLatGlobal(grid_template))then   ! Note that the grid template can be still not lat-global:
                                                  ! declared source grid can be larger than actual data area.
            call make_grid_lat_global(grid_template, .true.)  ! adjust resolution if needed
            if(error)return
          endif
        endif
      endif
      if(.not.fu_stdSilamGrid(grid_template))then
        call msg_warning('Area surce:'+ a_src%src_nm + '_' + a_src%sector_nm + &
                       & '. Resulting grid is not SILAM-standard one','create_src_contain_grd_a_src')
        call report(grid_template)
        call set_error('Resulting grid is not SILAM-standard one','create_src_contain_grd_a_src')
        call msg('')
        return
      endif
    endif  ! if not std SILAM grid

  end subroutine create_src_contain_grd_a_src


  !****************************************************************

  subroutine force_a_src_into_grid(a_src, grid_template, ifVerbose, ifCut)
    !
    ! Checks whether the source is inside the grid and, if some cells appear to be out
    ! deletes them from the cell list
    !
    implicit none

    ! Imported parameters
    type(silam_area_source), intent(inout) :: a_src
    type(silja_grid), intent(in) :: grid_template
    logical, intent(in) :: ifVerbose
    logical, intent(out) :: ifCut

    ! Local variables
    integer :: iSrc, nx, ny, iCell, nCellsOld
    real :: x, y
    logical :: ifReady

    if(.not. defined(a_src))then
      call set_error('Undefiend source','force_a_src_into_grid')
      call unset_error('force_a_src_into_grid')
      return
    endif

    if(a_src%ifFieldGiven)then
      !
      ! If the data are given in field form, just set the "storage grid",
      ! which will be enforced by the reading routine. But do it carefully
      !
      !
      ! Do nothing! smaller grid will be enfirced automatically when reprojecting
      !
      !
!      call cut_grid_size(a_src%grid, &    !grid_to_cut, 
!                       & grid_template, & !grid_area, 
!                       & inside_the_grid_area) !iTypeOfCut, grid_to_preserve)
    else
      !
      ! For each source we scan its central points and extending the template_grid 
      ! if these points are outside. For point and bomb sources it is only one point,
      ! while area source has plenty of them.
      !
      call grid_dimensions(grid_template, nx,ny)

      ifCut = .false.
      iCell = 1
      nCellsOld = a_src%nCells
      do while(iCell <= a_src%nCells)
        call project_point_to_grid(a_src%grid, &
                                 & a_src%cell_fx(iCell), &
                                 & a_src%cell_fy(iCell), &
                                 & grid_template, x, y)
        !
        ! Check the borders and, if not fitting, kill the cell by copying the last cell over
        !
        if(x<1 .or. x>nx .or. y<1 .or. y>ny) then
          if(iCell < a_src%nCells)then
            a_src%cell_fx(iCell) = a_src%cell_fx(a_src%nCells)
            a_src%cell_fy(iCell) = a_src%cell_fy(a_src%nCells)
            a_src%cell_Val(1:a_src%nDescriptors,iCell) &
                 & = a_src%cell_Val(1:a_src%nDescriptors,a_src%nCells)
          endif
          a_src%cell_fx(a_src%nCells) = real_missing
          a_src%cell_fy(a_src%nCells) = real_missing
          a_src%nCells = a_src%nCells - 1
          ifCut = .true.
        else
          iCell = iCell + 1
        endif
      end do
      if(ifVerbose .and. ifCut)then
        call msg('Has cut the cells from ... down to:',nCellsOld,a_src%nCells)
!        call report(a_src)
      endif   

      if(a_src%nCells == 0)then
        call msg('Source is entirely outside the grid',iSrc)
        call report(a_src)
        call report(grid_template)
        call set_error('Source is entirely outside the grid','force_a_src_into_grid')
        call unset_error('force_a_src_into_grid')
      endif

    endif  ! ifFieldGiven

  end subroutine force_a_src_into_grid


  !****************************************************************

  integer function fu_nbr_of_disp_grd_cells_a_src(a_src)
    !
    ! Returns the number of grid cells of the source - for the 
    ! dispersion grid
    !
    implicit none

    ! Imported paramrters
    type(silam_area_source), intent(in) :: a_src

    if(a_src%ifFieldGiven)then
      !
      ! The number of non-zero dispersion-grid cells from interpolation structure
      ! If no interpolation is needed... well, nothing to do, have to accept the whole grid
      !
      if(fu_true(a_src%ifHorizInterp)) then
        fu_nbr_of_disp_grd_cells_a_src = fu_n_non_zero_output_cells(a_src%pHorizIS)
      else
        fu_nbr_of_disp_grd_cells_a_src = fu_number_of_gridpoints(dispersion_grid)
      endif
    else
      !
      ! Just the size of disperion-cell list
      !
      if(a_src%nCellsDispGrd < 0 .or. a_src%nCellsDispGrd >= 1.0e6)then
        call msg('Strange number of dispersion-grid cells in source:' + a_src%src_nm, a_src%nCellsDispGrd)
        call msg_warning('Strange number of dispersion-grid cells in the source:' + a_src%src_nm, &
                       & 'fu_nbr_of_disp_grd_cells_a_src')
        fu_nbr_of_disp_grd_cells_a_src = 0
      else
        fu_nbr_of_disp_grd_cells_a_src = a_src%nCellsDispGrd
      endif
    endif

  end function fu_nbr_of_disp_grd_cells_a_src


  !****************************************************************

  subroutine source_2_map_area_source(a_src, dataPtr, id, iAccuracy, ifRandomise)
    !
    ! Projects the area source (in dispersion grid) to a map, having the given id as a 
    ! template: the id deterines the grid, level, time, and species name.
    ! Reprojection is done from dispersion grid.
    !
    ! ATTENTION. Do not scale or nullify dataPtr: this array can accumulate 
    !            information from many sources.
    !
    ! ATTENTION. This subroutine is tremendously inefficient. Reprojection is done 
    !            as many times as it is called - for each vertical level and for each
    !            descriptor. This is largely forced solution because, in case of split
    !            sources in the emission output file, we will not be able to keep them
    !            all in memory while collecting species. Some simplification comes 
    !            from iAccuracy but still...
    !
    implicit none

    ! Imported parameters
    type(silam_area_source), intent(inout) :: a_src
    real, dimension(:), intent(inout) :: dataPtr
    type(silja_field_id), intent(in) :: id
    integer, intent(in) :: iAccuracy
    logical, intent(in) :: ifRandomise

    ! Local variables
    real :: xOut, yOut, xNew, yNew, fMassSkipped, f_nSmall_2_1, scaling_factor, xL, xR, yL, yU, &
          & fX_ll, fy_ll, fX_ur, fY_ur, fX_save
    integer :: iDescr, i, j, k, iCell, nxSrc, nySrc, nxOut, nyOut, nSmall, nSpecies, iWrap, nWrap
    integer :: nSmallSkipped, ix, iy
    integer, dimension(max_species) :: species_index 
    type(silam_species), dimension(:), pointer :: species
    real, dimension(:), pointer :: amounts 
    logical :: ifFound, ifSimpleReprojection, ifLonGlobalOut
    type(silja_rp_1d), dimension(:), pointer :: pArTmp 
    type(THorizInterpStruct), pointer :: pHorInterp 
    type(silja_grid) :: idGrid
    character(len=*), parameter :: sub_name = 'source_2_map_area_source'

    species =>null()
    amounts  =>null()
    pArTmp =>null()
    pHorInterp =>null()
    idGrid = fu_grid(id)
    !
    ! Stupidity check
    !
    if(.not.defined(id))then
      call set_error('Undefined id given',sub_name)
      return
    endif
    if(.not.defined(idGrid))then
      call set_error('Undefined grid given',sub_name)
      return
    endif
    if(fu_number_of_gridpoints(idGrid) > size(dataPtr))then
      call set_error('Too small data array given',sub_name)
      return
    endif
    if(iAccuracy < 1 .or. iAccuracy > 10)then
      call msg('Accuracy switch must be from 1 to 10, not:',iAccuracy)
      call set_error('Accuracy switch must be from 1 to 10',sub_name)
      return
    endif
if (ifRandomise)then
    call msg('Source reprojection with randomization')
else
    call msg('Source reprojection without randomization')
endif
    !
    ! Compare the source species and the requested template. If the source does not contribute, return
    ! Note that we have to check all descriptors.
    !
    if(error)return
    ifFound = .false.
    do iDescr = 1, a_src%nDescriptors
      call get_inventory(a_src%cocktail_descr_lst(iDescr), species, nSpecies)
      species_index(iDescr) = fu_index(fu_species(id), species, nSpecies)
      if(species_index(iDescr) >= 1 .and. species_index(iDescr) <= nSpecies) ifFound = .true.
    end do
    if(.not. ifFound)then
      if(len_trim(a_src%sector_nm) > 0)then
        call msg('Source:' + a_src%src_nm + '_' + a_src%sector_nm + &
               & '- does not contribute to species:' + fu_str(fu_species(id)))
      else
        call msg('Source:' + a_src%src_nm + &
               & '- does not contribute to species:' + fu_str(fu_species(id)))
      endif
      return
    endif

    call grid_dimensions(idGrid, nxOut, nyOut)
    if(error)return

    if(a_src%ifFieldGiven)then
      !
      ! Have to sum-up the binaries
      ! First, prepare the strucutre
      ! 
      call grid_dimensions(a_src%grid, nxSrc, nySrc)
      call get_work_arrays_set(a_src%nDescriptors, nxSrc*nySrc, pArTmp)
      if(error)return
      do i = 1, a_src%nDescriptors
        pArTmp(i)%pp(1:a_src%nCells) = 0.0
      end do
      !
      ! Now use this structure to get maps of time-wise totals over the layer 
      ! for all descriptors from the binaries
      !
      do j = 1, size(a_src%FieldFormat)
        select case(a_src%FieldFormat(j)%iformat)  ! GrADS/NetCDF/..., ref to input structure
          case(grads_file_flag)
            call get_grads_total(a_src%uField(j), a_src%pFldEmis, &
                               & fu_accumulation_start_time(id), fu_accumulation_length(id), &
                               & fu_level(id), &
                               & pArTmp)
            if(error)return
          case(netcdf_file_flag)
            call get_netcdf_total(a_src%uField(j), a_src%pFldEmis, &
                                & fu_accumulation_start_time(id), fu_accumulation_length(id), &
                                & fu_level(id), &
                                & pArTmp)
            if(error)return
          case(ascii_file_flag) 
            call set_error('ASCII field format is not supported as time-resolving', &
                         & sub_name)
            return
          case default
            call set_error('Field format is not supported:' + fu_str(a_src%FieldFormat(j)%iformat), &
                         & sub_name)
             return
        end select

      enddo    ! binary formats
      !
      ! Now we have summed-up maps for all descriptors. Need to do the actual reprojecting
      ! and splitting of the descriptors geting out the needed species
      !
      pHorInterp => fu_horiz_interp_struct(a_src%grid, idGrid, summation, ifRandomise, iOutside=setZero)
      if(error .or. .not. associated(pHorInterp))then
        call set_error('Failed to get horizontal interpolation structure',sub_name)
        return
      endif

      do iDescr = 1, a_src%nDescriptors
        if(species_index(iDescr) < 1 .or. species_index(iDescr) > nSpecies)cycle ! void descriptor

        scaling_factor = a_src%params(1)%rate_descr_unit(iDescr)* &    ! from the source unit to descriptor unit
                       & a_src%fDescr2SpeciesUnit(species_index(iDescr),iDescr)  ! from descriptor to spcies unit

        do iy = 1, nyOut
          do ix = 1, nxOut
            k = ix + (iy-1)*nxOut
            do i = 1, pHorInterp%nCoefs
              if(pHorInterp%indX(i,ix,iy) == 0)exit  ! all indices are used up
              dataPtr(k) = dataPtr(k) + pArTmp(iDescr)%pp(pHorInterp%indX(i,ix,iy) + &
                                                        & (pHorInterp%indY(i,ix,iy)-1)*nxSrc) * &
                                      & pHorInterp%weight(i,ix,iy) * scaling_factor
            end do
          end do
        end do

      end do ! descriptors
      
      call free_work_array(pArTmp)

    else
      !---------------------------------------------------------------------------------
      !
      ! Source is given as a list of cells. Go hard way.
      !
      if(a_src%nCells < 1)then
        call set_error('No cells in area source:' + a_src%src_nm,sub_name)
        call unset_error(sub_name)
        return
      endif

      nSmallSkipped = 0
      fMassSkipped = 0. 
      !
      ! Get the array of amounts per-descriptor
      !
      amounts => fu_work_array(a_src%nDescriptors)
      if(error)return
      call total_from_a_src_descr_unit(a_src, .true. , &   ! if integrate rates only
                                     & amounts, &
                                     & fu_accumulation_start_time(id), &
                                     & fu_accumulation_length(id), &
                                     & fu_level(id))
      if(error)return

      if(sum(amounts(1:a_src%nDescriptors)) .eps. 0.0)then
        call free_work_array(amounts)
        if(len_trim(a_src%sector_nm) > 0)then
          call msg('Source:' + a_src%src_nm + '_' + a_src%sector_nm + &
                 & '- has zero total for species:' + fu_str(fu_species(id)))
        else
          call msg('Source:' + a_src%src_nm + &
                 & '- has zero total for species:' + fu_str(fu_species(id)))
        endif
        return
      endif

      !
      ! For lonlat grids with the same poles the reprojection can be substantially sped-up
      ifSimpleReprojection = fu_gridtype(dispersion_grid) == lonlat .and. &
                           & fu_gridtype(idGrid) == lonlat
      if(ifSimpleReprojection) ifSimpleReprojection =  fu_pole(dispersion_grid) == fu_pole(idGrid)


      if (idGrid == dispersion_grid) then 
        ! No remapping needed
        if(len_trim(a_src%sector_nm) > 0)then
          call msg('NO mapping for area source:' + a_src%src_nm + '_' + a_src%sector_nm + ',' + &
                                    & fu_str(fu_species(id)))
        else
          call msg('NO mapping for area source:' + a_src%src_nm + ',' + fu_str(fu_species(id)))
        endif
        do iCell = 1, a_src%nCellsDispGrd ! Only these cells are non-zero
            k = nint(a_src%cellDispGrd_fx(iCell)) + (nint(a_src%cellDispGrd_fy(iCell))-1)*nxOut
            do iDescr = 1, a_src%nDescriptors
              if(species_index(iDescr) < 1)cycle     ! Some descriptors may not have the needed species
              dataPtr(k) = dataPtr(k) + a_src%cellDispGrd_val(iDescr,iCell) * amounts(iDescr) * &
                            & a_src%fDescr2SpeciesUnit(species_index(iDescr),iDescr)
            end do
        enddo

      elseif(ifSimpleReprojection)then
        !
        ! Poles are the same, i.e. the grid cell borders are parallel. Then a simple algebra can 
        ! serve the reprojection
        ! This is the direct projection from -> to
        !
        if(len_trim(a_src%sector_nm) > 0)then
          call msg('Simplified mapping area source:' + a_src%src_nm + '_' + a_src%sector_nm + ',' + &
                                    & fu_str(fu_species(id)))
        else
          call msg('Simplified mapping area source:' + a_src%src_nm + ',' + fu_str(fu_species(id)))
        endif
        ifLonGlobalOut = fu_ifLonGlobal(idGrid)

        do iCell = 1, a_src%nCellsDispGrd ! Only these cells are non-zero
          !
          ! Lower left and upper right corners determine the coordinates range in the From grid
          !
          ix = nint(a_src%cellDispGrd_fx(iCell))
          iy = nint(a_src%cellDispGrd_fy(iCell))
          call project_point_to_grid(dispersion_grid, ix-0.5, iy-0.5, idGrid, fX_ll, fY_ll)
          call project_point_to_grid(dispersion_grid, ix+0.5, iy+0.5, idGrid, fX_ur, fY_ur)
          !
          ! Do we yeed to bother?
          !
          if (ifLonGlobalOut) then 
            if(.not.(fY_ur >0.5 .and. fY_ll <nYOut+0.5))cycle
          else
            if(.not.(fX_ur >0.5 .and. fX_ll <nxOut+0.5 .and. fY_ur >0.5 .and. fY_ll <nYOut+0.5))cycle
          endif

          ! we do...
          ! Cycle over the covered To-grid cells sending there a corresponding From mass fraction 
          ! Skip the out the out-of-grid cells
          !
          fX_save = fX_ur
          if (ifLonGlobalOut .and. fX_ll > fX_ur) then ! Handling of global wrapping
            nWrap=2
            fX_ur = nxOut+0.5
          else
            nWrap=1
          endif

          do iWrap = 1,nWrap
            do j = max(1,nint(fY_ll)), min(nyOut,nint(fY_ur))
              do i = max(1,nint(fX_ll)), min(nxOut,nint(fX_ur))
                !
                ! The fraction of mass of the (iSml,jSml) FROM cell that goes into (ix,iy) TO cell
                !
                xL = max(fX_ll, i-0.5)  ! lower left
                xR = min(fX_ur, i+0.5)  ! upper right
                yU = min(fY_ur, j+0.5)
                yL = max(fY_ll, j-0.5)
                f_nSmall_2_1 = (xR - xL) / (fX_ur - fX_ll) * (yU - yL) / (fY_ur - fY_ll)
                if(f_nSmall_2_1 < 1e-10)cycle  ! for the case of nearly-perfect hit of the boundaries

                k = i + (j-1)*nxOut
                do iDescr = 1, a_src%nDescriptors
                  if(species_index(iDescr) < 1)cycle     ! Some descriptors may not have the needed species
                  dataPtr(k) = dataPtr(k) + a_src%cellDispGrd_val(iDescr,iCell) * amounts(iDescr) * &
                                & f_nSmall_2_1 * a_src%fDescr2SpeciesUnit(species_index(iDescr),iDescr)
                end do
              end do  ! ixTo
            end do  ! iyTo
            fX_ur = fX_save 
            fX_ll = 0.5
          enddo !iWrap
        end do  ! iCell
          
      else
        !
        ! Poles are different. No other way, have to go for tough work.
        ! Select the number of small cells, which the original grid cell will be broken into
        !
        ! The recepy for the flux-conserving reprojection is: split each non-zero cell in
        ! fArea to, say, 100 small cells and re-project them as another small grid.
        ! This, of course, is tremendously inefficient once again, but once-per-run can be
        ! allowed.
        !
        ! The number of small emission cells depend on the output grid cell size.
        ! It will be an area ratio of source and output grids * 2000 and then min-max
        ! to limit the variation of the number. Limits are: explicit 1001 and 101, plus
        ! the total number of small cells to proceed: I do not want to deal with tens of 
        ! millions of grid reprojections.
        !
        nSmall = int(sqrt(2000. * fu_cell_size(dispersion_grid) / fu_cell_size(idGrid))/2.)*2+1

        nSmall = min(1001, max(nSmall,101)) ! Hard limits

        nSmall = min(nSmall, int(sqrt(5.e4/a_src%nCells))*2+1) ! nSmall*nSmall*as%nCells < 2*10**6

        nSmall = max(11, nSmall) ! Still not too crude interpolation...

        nSmall = nSmall * iAccuracy / 10
        nSmall = max(nSmall, 3)  ! For reduced accuracy this will be the hard limit

        f_nSmall_2_1 = 1. / real(nSmall * nSmall)

        if(len_trim(a_src%sector_nm) > 0)then
          call msg('Mapping area source:' + a_src%src_nm + '_' + a_src%sector_nm + ',' + &
                                    & fu_str(fu_species(id)) +&
                                    & ', nSmall for emission:', nSmall)
        else
          call msg('Mapping area source:' + a_src%src_nm + ',' + &
                 & fu_str(fu_species(id)) + ', nSmall for emission:', nSmall)
        endif
        !
        ! Finally, can go along the grid cells, not missing the integrated amounts as a general 
        ! multiplier
        !
        do iCell = 1, a_src%nCells ! Only these cells are non-zero

          !
          ! Split the source grid cells to some number of small elements and project
          ! each of them separately. A nice feature is: grid co-ordinates can be real
          ! values
          !
          do i=1,nSmall
            do j=1,nSmall

              ix = nint(a_src%cellDispGrd_fx(iCell))
              iy = nint(a_src%cellDispGrd_fy(iCell))
              if(ifRANDOMISE)then
                 call project_point_to_grid(dispersion_grid, &
                   & ix-0.5 + (real(i)-0.5)/nSmall + fu_random_number_center(0.,1./real(nSmall)), &
                   & iy-0.5 + (real(j)-0.5)/nSmall + fu_random_number_center(0.,1./real(nSmall)), &
                   & idGrid, xNew, yNew)
              else
                 call project_point_to_grid(dispersion_grid, &
                   & ix-0.5 + (real(i)-0.5)/nSmall, &
                   & iy-0.5 + (real(j)-0.5)/nSmall, &
                   & idGrid, xNew, yNew)
              endif
              if(error)return

              if(xNew >= 0.5 .and. xNew <= nxOut+0.5 .and. yNew >= 0.5 .and. yNew <= nyOut+0.5)then
                k = nint(xNew) + (nint(yNew)-1)*nxOut
                do iDescr = 1, a_src%nDescriptors
                  if(species_index(iDescr) < 1)cycle     ! Some descriptors may not have the needed species
                  dataPtr(k) = dataPtr(k) + a_src%cellDispGrd_val(iDescr,iCell) * amounts(iDescr) * &
                                     & f_nSmall_2_1 * a_src%fDescr2SpeciesUnit(species_index(iDescr),iDescr)
                end do
              else
                nSmallSkipped = nSmallSkipped + 1
                do iDescr = 1, a_src%nDescriptors
                  if(species_index(iDescr) < 1)cycle    ! Some descriptors may not have the needed species
                  fMassSkipped = fMassSkipped + a_src%cellDispGrd_val(iDescr,iCell) * amounts(iDescr) * &
                                     & f_nSmall_2_1 * a_src%fDescr2SpeciesUnit(species_index(iDescr),iDescr)
                end do
              endif
            end do   ! nSmall
          end do   ! nSmall

        end do ! Cycle through the source non-zero cells
        
      endif  ! if same pole

      call msg('Total sum for output map:', sum(dataPtr(1:fu_number_of_gridpoints(idGrid))))
      call msg('Emiited to dispersion out of output grid: ', fMassSkipped)

      call free_work_array(amounts)

    endif ! ifFieldGiven

  end subroutine source_2_map_area_source


  !*****************************************************************

  subroutine project_a_src_second_grd(as, grid, vert_disp, vert_proj, iAccuracy, tmpCellDispGrd_val, &
                                    & ifRandomise, rng)
    !
    ! Duplicates the whole set of emission cells re-projecting them to the given new grid
    ! The new cell set is stored as a parallel one to the original list. Reasons for such 
    ! solution are (i) we either store it once and forever or make reprojecting at each 
    ! time step; (ii) having it done once, we can allow the similar approach as above in
    ! source_2_map_area_source and thus ensure the maximum possible coherence.
    !
    implicit none

    ! Imported parameters
    type(silam_area_source), intent(inout) :: as
    type(silja_grid), intent(in) :: grid
    ! vert_proj determines the centre of mass and fractions, vert_disp converts to the
    ! dispersion vertical unit.
    type(silam_vertical), intent(in) :: vert_disp, vert_proj
    integer, intent(in) :: iAccuracy
    type(silja_rp_1d), dimension(:), pointer :: tmpCellDispGrd_val
    logical, intent(in) :: ifRandomise
    type(rng_type), intent(inout) :: rng ! Using random-numbers changes the
                                         ! state of the generator

    ! Local variables
    real :: xNew, yNew, fSmall_2, fMassSml, fMassStored, fX_ll, fY_ll, fX_ur, fY_ur, xL, xR, yL, yU, fTmp
    real :: fxTmp, fyTmp
    integer :: iNewCell, nNewCells, nxOut, nyOut, nSmall, iCell, i,j,k, nSmallSkipped, iDescr
    integer :: iLev, iThread, nSmallTmp, iTry, iTmp
    integer, dimension(:), pointer :: iCellNbrs 
    real, dimension(:), pointer :: tmpCellDispGrd_fx, tmpCellDispGrd_fy
    character(len=clen) :: strTmp
    character(len=fnlen) :: str1Tmp !Longer string for reporting
    type(silam_species) :: speciesTmp
    type(silam_species), dimension(:), pointer :: speciesArPtr
    type(silja_field_id) :: idTmp
    logical :: ifSimpleReprojection

    type(THorizInterpStruct), pointer :: pHIS
    character(len=*), parameter :: sub_name = 'project_a_src_second_grd'

    nSmall=int_missing

   if (as%ifFieldGiven) then  !! Not exactly reprojection, but allocation once dispersion grid known

      iTmp = summation
      if (as%ifFluxes) iTmp = average 
      !$OMP CRITICAL  (HIS_project_a_src_second_grd)
        pHIS => fu_horiz_interp_struct(as%grid, dispersion_grid, summation, ifRandomise, iAccuracy, setzero, ifCreate=.true.)
      !$OMP END CRITICAL (HIS_project_a_src_second_grd)
      if (pHIS%ixStartTo > pHIS%ixEndTo) then  
          call msg("Spource "//trim(as%src_nm)//":"//trim(as%sector_nm) //" outside of our grid. Disablng!")
          as%ifHorizInterp = silja_false
          as%pHorizIS => HorizInterpStruct_missing
          as%ifVertInterp = silja_false
          as%pVertIS => VertInterpStruct_missing
          as%nCellsDispGrd = 0
          return
      endif 
      


      !
      ! Binary contains all info. Use the bunch of field_3d to store the thing
      !
      allocate(as%pFldEmis(as%nDescriptors), stat=i)
      if(fu_fails(i==0,'Failed to allocate emission field_3d',sub_name))return

!      call set_missing(as%pShopLst) 
      do iDescr = 1, as%nDescriptors
        !
        ! Allocate the 3d field pointers
        allocate(as%pFldEmis(iDescr)%fp, stat=i)
        if(fu_fails(i==0,'Failed emission field_3d%fp allocation',sub_name))return
        call set_3d_field_empty(as%pFldEmis(iDescr)%fp)
        !
        ! Species or cocktail?
        !
        if(as%ifSpecies(iDescr))then
          speciesArPtr => fu_species(as%cocktail_descr_lst(iDescr))
          speciesTmp = speciesArPtr(1)
          strTmp = ''
        else
          speciesTmp = species_missing
          strTmp = fu_name(as%cocktail_descr_lst(iDescr))
        endif
        !
        ! And then add levels  
        !
        ! HUOM! Fields stored in dispersion grid!
        do iLev = 1, fu_NbrOfLevels(as%vertLevs)
          idTmp = fu_set_field_id(met_src_missing,&
                                & emission_intensity_flag, &  !!!No matter what, we keep intensity, not fluxes
                                & as%times(1), &          !analysis_time,&
                                & zero_interval, &            !forecast_length, &
                                & dispersion_grid, &
                                & fu_level(as%vertLevs,iLev), &
                                & field_kind = averaged_flag, &  !!!Average fields 
                                & length_of_accumulation = zero_interval, & !!Will be reset on actual read
                                & species = speciesTmp, &
                                & chCocktail=strTmp)
          if(error)return

          call create_empty_field_in_field_3d(as%pFldEmis(iDescr)%fp, idTmp)
          if(error)return

        enddo ! iLev

        call set_3dField_params_from_fieldId(as%pFldEmis(iDescr)%fp, idTmp)
        if(error)return
        call organize_fields_vertically(as%pFldEmis(iDescr)%fp)
        if(error)return

      enddo  ! Descriptors

      as%nCells = int_missing !!! Should not be used
      as%nCellsDispGrd = fu_number_of_gridpoints(dispersion_grid)



    endif
    !
    ! If data are given in the binary file, all what we need is to create 
    ! interpolation structures
    !
    if(as%grid == grid .or. as%ifFieldGiven)then
      as%ifHorizInterp = silja_false
      as%pHorizIS => HorizInterpStruct_missing
    else
      as%ifHorizInterp = silja_true
    endif

    if(fu_cmp_verts_eq(as%vertLevs, vert_proj))then
      as%ifVertInterp = silja_false
      as%pVertIS => VertInterpStruct_missing
    else
      as%ifVertInterp = silja_true
      if(as%ifFieldGiven)then
          if (fu_if_layer(fu_level(as%vertLevs,1))) then
               !$OMP CRITICAL(mk_v_interp_struct)
               as%pVertIS => fu_vertical_interp_struct(&
                                                & as%vertLevs, vert_proj, &
                                                & grid, summation, & ! vertical projection in dispersion grid!
                                                & one_hour*3., &
                                                & 'a_src_to_emis_MassMap_' + fu_str(as%id_nbr))
              !$OMP END CRITICAL(mk_v_interp_struct)
          elseif (as%ems2wholeABL) then
            as%pVertIS => null()
            call msg_warning("Surface-only area source fields..", sub_name)
          else
            call set_error("Non-layer area source fields with 3D interp", sub_name)
          endif
      endif
    endif

    if(as%ifFieldGiven) return ! that's it for the time being

    !------------------------------------------------------------------------
    !
    ! If source is a list of cells, reprojection has to go hard way
    !
    if(iAccuracy < 1 .or. iAccuracy > 10)then
      call msg('Accuracy switch must be from 1 to 10, not:',iAccuracy)
      call set_error('Accuracy switch must be from 1 to 10',sub_name)
      return
    endif

    ! May be, the given grid is the same as in the source?
    ! Then just redirect the pointer and return
    !
    if(as%grid == grid)then
      as%cellDispGrd_fx => as%cell_fx
      as%cellDispGrd_fy => as%cell_fy
      as%cellDispGrd_val => as%cell_val
      as%nCellsDispGrd = as%nCells
      if(.not. as%ifMultiLevelFixed)return
      as%vertLevsDispVert = vert_disp
      allocate(as%levFractDispVert(fu_NbrOfLevels(vert_disp)),as%fzDisp(fu_NbrOfLevels(vert_disp)), stat=i)
      if(fu_fails(i == 0, 'Failed dispersion-vertical level fractions allocation',sub_name))return
      as%levFractDispVert(:) = 0.0
      as%fzDisp(:) = 0.0

      call reproject_verticals(as%vertLevs, as%levFraction, &   ! vertical from, fractions from
                             & vert_proj, as%levFractDispVert, &     ! vertical to, fractions to
                             & as%fzDisp, as%nLevsDispVert, & ! mass centres, number of non-zero levels
                             & ifMassCentreInRelUnit=.true.)

      return
    endif


    call grid_dimensions(grid, nxOut, nyOut)
    ! 
    ! We do not know how many cells will be in the given grid. Have to go carefully
    ! So far, assume the fraction ratio as a scaling for the number of cells
    !
    ! tripled the number - should be sufficient for everything
    !

    !!! RK: This calculation was rewritten to avoid integer overflow
    !!! on very coarse sourcees over very fine dispersion grids
    as%nCellsDispGrd = int( min( nxOut*nyOut*1. ,  &
                    &   3. * as%nCells * fu_cell_size(as%grid) / fu_cell_size(grid) + 1.))


    iCellNbrs => fu_work_int_array(nxOut*nyOut)
    if(error)return
    iCellNbrs(1:nxOut*nyOut) = int_missing

    nullify(as%cellDispGrd_fx);     !call set_array_size(as%cellDispGrd_fx, as%nCellsDispGrd, 0.0)
    nullify(as%cellDispGrd_fy);     !call set_array_size(as%cellDispGrd_fy, as%nCellsDispGrd, 0.0)
    nullify(as%cellDispGrd_val);    !call set_array_size(as%cellDispGrd_val, as%nDescriptors, as%nCellsDispGrd, 0.0)
    if(error)return
    !
    ! Get temporary space for the whole thing. It is better than multiple reallocation. Note also that some
    ! of the stuff may already be allocated
    !
    tmpCellDispGrd_fx => fu_work_array(as%nCellsDispGrd)
    tmpCellDispGrd_fy => fu_work_array(as%nCellsDispGrd)
    if(size(tmpCellDispGrd_val) < as%nDescriptors .or. &
     & size(tmpCellDispGrd_val(1)%pp) < as%nCellsDispGrd)then
      call free_work_array(tmpCellDispGrd_val)
      call get_work_arrays_set(as%nDescriptors, as%nCellsDispGrd, tmpCellDispGrd_val)
      if(error)return
    endif
    nNewCells = 0  ! number of filled-in cells in dispersion_grid, <= nCellsDispGrd
    !
    ! We have to scale each cell with its area, in order to get the
    ! emission density, which can be interpolated with muuta3 (hila.lib).
    ! Apart from this, it is worth of making a bit more accurate re-projecting
    ! than just a simple nearest-point trick: split each cell to N small cells 
    ! and re-project them one-by-one. This, of course, is tremendously inefficient, 
    ! but once-per-run can be allowed.
    !
    ifSimpleReprojection = fu_gridtype(as%grid) == lonlat .and. &
                         & fu_gridtype(grid) == lonlat
    if(ifSimpleReprojection) ifSimpleReprojection = ifSimpleReprojection .and. &
                                                  & fu_pole(as%grid) == fu_pole(grid)
    if(ifSimpleReprojection)then
      !
      ! Poles are the same, i.e. the grid cell borders are parallel. Then a simple algebra can 
      ! serve the reprojection
      ! This is the direct projection from -> to
      !

      do iCell = 1, as%nCells ! Only these cells are non-zero
        !
        ! Lower left and upper right corners determine the coordinates range in the From grid
        !
        call project_point_to_grid(as%grid, as%cell_fx(iCell)-0.5, as%cell_fy(iCell)-0.5, &
                                 & grid,fX_ll,fY_ll)
        call project_point_to_grid(as%grid, as%cell_fx(iCell)+0.5, as%cell_fy(iCell)+0.5, &
                                 & grid,fX_ur,fY_ur)
        !
        ! Do we yeed to bother?
        !
        if(.not.(fX_ur >0.5 .and. fX_ll < nxOut+0.499 .and. fY_ur >0.5 .and. fY_ll < nYOut+0.499))cycle

        ! we do...
        ! Cycle over the covered To-grid cells sending there a corresponding From mass fraction 
        ! Skip the out the out-of-grid cells
        !
        do j = max(1,nint(fY_ll)), min(nyOut,nint(fY_ur))
          do i = max(1,nint(fX_ll)), min(nxOut,nint(fX_ur))
            !
            ! The fraction of mass of the (iSml,jSml) FROM cell that goes into (ix,iy) TO cell
            !
            xL = max(fX_ll, i-0.5)  ! lower left
            xR = min(fX_ur, i+0.5)  ! upper right
            yL = max(fY_ll, j-0.5) 
            yU = min(fY_ur, j+0.5)   
            fSmall_2 = (xR - xL) / (fX_ur - fX_ll) * (yU - yL) / (fY_ur - fY_ll)
            ! Margin must be more than twice machine_epsilon 
            ! Otherwise mass_centers of -0.5  occur
            if(fSmall_2 < 1e-6)cycle  ! for the case of nearly-perfect hit of the boundaries
            fMassSml = sum(as%cell_Val(1:as%nDescriptors,iCell)) * fSmall_2
            if (.not. fMassSml > small_epsilon) cycle
            k = i + (j-1)*nxOut


            ! Move centers towards cell center
            xL = max(fX_ll, i-0.4999)  ! lower left
            xR = min(fX_ur, i+0.4999)  ! upper right Must never round to upper!!
            yL = max(fY_ll, j-0.4999) 
            yU = min(fY_ur, j+0.4999)   
            
            if(iCellNbrs(k) == int_missing)then
              !
              ! The cell is not yet occupied, take the new and set the pointer to the new cell
              !
              nNewCells = nNewCells + 1
              if(nNewCells > size(tmpCellDispGrd_fx))then
                as%nCellsDispGrd = size(tmpCellDispGrd_fx) * 5 / 4 !Enlarge exponentially
                call set_array_size(tmpCellDispGrd_fx, as%nCellsDispGrd)
                call set_array_size(tmpCellDispGrd_fy, as%nCellsDispGrd)
                call enlarge_array(tmpCellDispGrd_val, as%nDescriptors, as%nCellsDispGrd)
                call msg('Enlarged:',as%nCellsDispGrd)
              endif  ! if vector should be extended

              tmpCellDispGrd_fx(nNewCells) = (xL+xR) * 0.5 * fMassSml
              tmpCellDispGrd_fy(nNewCells) = (yL+yU) * 0.5 * fMassSml
              do iDescr = 1, as%nDescriptors
                tmpCellDispGrd_Val(iDescr)%pp(nNewCells) = as%cell_Val(iDescr,iCell) * fSmall_2
              enddo
              iCellNbrs(k) = nNewCells

            else
              !
              ! Adjust the existing cell and add the corresponding mass
              !
              iNewCell = iCellNbrs(k)

              tmpCellDispGrd_fx(iNewCell) = tmpCellDispGrd_fx(iNewCell) + (xL+xR) * 0.5 * fMassSml
              tmpCellDispGrd_fy(iNewCell) = tmpCellDispGrd_fy(iNewCell) + (yL+yU) * 0.5 * fMassSml
              do iDescr = 1, as%nDescriptors
                tmpCellDispGrd_Val(iDescr)%pp(iNewCell) = tmpCellDispGrd_Val(iDescr)%pp(iNewCell) + &
                                                        & as%cell_Val(iDescr,iCell) * fSmall_2
              end do
           endif  ! if cell is already in the disp cells list
          end do  ! ixTo
        end do  ! iyTo
      end do  ! iCell
          
    else
      !
      ! Poles are different. No other way, have to go for tough work.
      ! Select the number of small cells, which the original grid cell will be broken into
      !

      !
      ! The number of small emission cells depend on the output grid cell size.
      ! It will be an area ratio of source and output grids * 1000 and then min-max
      ! to limit the variation of the number. Limits are: explicit 1001 and 101, plus
      ! the total number of small cells to proceed: I do not want to deal with tens of 
      ! millions of grid reprojections.
      !
      nSmall = int(sqrt(2000. * fu_cell_size(as%grid) / fu_cell_size(grid))/2.)*2+1

      nSmall = min(101, max(nSmall,11)) ! Hard limits


!FIXME !!! NONSENSE! nSmall should not depend on as%nCells!!!! R      
!!      if (as%nCells > 0) then
!!        nSmall = min(nSmall, int(sqrt(5.e5/as%nCells))*2+1) ! nSmall*nSmall*as%nCells < 2*10**6
!!      end if
!!
!!      nSmall = max(11, nSmall) ! Still not too crude interpolation...

      nSmall = nSmall * iAccuracy / 10
      if(nSmall < 3)nSmall = 3  ! For reduced accuracy this will be the hard limit
    
      fSmall_2 = 1. / real(nSmall*nSmall)

      do iCell = 1, as%nCells ! Only these cells are non-zero
        fMassSml = sum(as%cell_Val(1:as%nDescriptors,iCell)) * fSmall_2
        if (.not. fMassSml > small_epsilon) cycle

        nSmallTmp=2 !Prepare for trial reprojection

        do iTry=1,2  ! First check if any corner within grid 
                    ! Then do the reprojection or skip the cell
          nSmallSkipped = 0
          do i=1,nSmallTmp
            do j=1,nSmallTmp
              if (iTry==1) then
                ! Corners
                fxTmp = as%cell_fx(iCell) -1.5 + real(i) 
                fyTmp = as%cell_fy(iCell) -1.5 + real(j)
              else
                fxTmp = as%cell_fx(iCell)-0.5 + (real(i)-0.5)/real(nSmall)
                fyTmp = as%cell_fy(iCell)-0.5 + (real(j)-0.5)/real(nSmall)
                if (ifRANDOMISE) then
                  fxTmp = fxTmp + fu_random_number_center_rng(rng, 0.,1./real(nSmall))
                  fyTmp = fyTmp + fu_random_number_center_rng(rng, 0.,1./real(nSmall))
                endif
              endif
               call project_point_to_grid(as%grid, fxTmp, fyTmp, grid, xNew, yNew)
    !call msg('x,y:',xNew,yNew)
              if(error)return
              if(xNew <= 0.501 .or. xNew >= nxOut+0.499)then
    !            call msg('Point is outside the x-limits.iCell;x=',iCell,xNew)
    !            call msg_warning('Point is outside the x-limits.Skipping',sub_name)
                nSmallSkipped = nSmallSkipped + 1
              elseif(yNew <= 0.501 .or. yNew >= nyOut+0.499)then
    !            call msg('Point is outside the y-limits.iCell;y=',iCell,yNew)
    !            call msg_warning('Point is outside the y-limits.Skipping',sub_name)
                nSmallSkipped = nSmallSkipped + 1
              elseif (iTry==2) then !Only now reproject
                !
                ! The emission comes inside the output grid. Check if such cell is already
                ! stored as non-zero. If yes, add new portion of mass, otherwise use new cell
                ! For quick search we use the lookup table, which says which dispersion grid cell
                ! is served by which element in as%cellsDispGrd array
                !
                if(iCellNbrs(nint(xNew)+(nint(yNew)-1)*nxOut) == int_missing)then
                  !
                  ! The cell is not yet occupied, take the new and set the pointer to the new cell
                  !
                  nNewCells = nNewCells + 1
                  if(nNewCells > size(tmpCellDispGrd_fx))then
                    as%nCellsDispGrd = size(tmpCellDispGrd_fx) * 5 / 4 !Enlarge exponentially
                    call set_array_size(tmpCellDispGrd_fx, as%nCellsDispGrd)
                    call set_array_size(tmpCellDispGrd_fy, as%nCellsDispGrd)
                    call enlarge_array(tmpCellDispGrd_val, as%nDescriptors, as%nCellsDispGrd)
                    call msg('Enlarged:',as%nCellsDispGrd)
                  endif  ! if vector should be extended

                  tmpCellDispGrd_fx(nNewCells) = xNew * fMassSml
                  tmpCellDispGrd_fy(nNewCells) = yNew * fMassSml
                  do iDescr = 1, as%nDescriptors
                    tmpCellDispGrd_Val(iDescr)%pp(nNewCells) = as%cell_Val(iDescr,iCell) * fSmall_2
                  enddo
                  iCellNbrs(nint(xNew)+(nint(yNew)-1)*nxOut) = nNewCells

                else
                  !
                  ! Adjust the existing cell and add the corresponding mass
                  !
                  iNewCell = iCellNbrs(nint(xNew)+(nint(yNew)-1)*nxOut)

                  tmpCellDispGrd_fx(iNewCell) = tmpCellDispGrd_fx(iNewCell) + xNew * fMassSml
                  tmpCellDispGrd_fy(iNewCell) = tmpCellDispGrd_fy(iNewCell) + yNew * fMassSml
                  do iDescr = 1, as%nDescriptors
                    tmpCellDispGrd_Val(iDescr)%pp(iNewCell) = tmpCellDispGrd_Val(iDescr)%pp(iNewCell) + &
                                                            & as%cell_Val(iDescr,iCell) * fSmall_2
                  end do
               endif  ! if cell is already in the disp cells list
              endif  ! if new small cell coordinates are outside the grid

            end do  ! nSmall cycle
          end do  ! nSmall cycle

          if (nSmallSkipped == 4) exit ! All corners outside
          nSmallTmp = nSmall !Prepare for reprojection of the cell
        enddo !iTry


      end do ! iCell, Cycle through the source non-zero cells

    endif  ! if simple reprojection    
    if(error)return
    
    as%nCellsDispGrd = nNewCells
    allocate(as%cellDispGrd_Val(as%nDescriptors,nNewCells), as%cellDispGrd_fx(nNewCells), &
           & as%cellDispGrd_fy(nNewCells), stat = iCell)
    if(fu_fails(iCell==0,'Failed allocation of dispersion cells',sub_name))return


    !
    ! Scale/average if needed and compute the longitude of each new grid cell
    ! Averaging is needed if the source is actually a fraction of land or whatever else
    ! but not <mass>/<time>
    !
    do iCell = 1, as%nCellsDispGrd
      fMassStored = 0.
      do iDescr = 1, as%nDescriptors
        fTmp = tmpCellDispGrd_Val(iDescr)%pp(iCell) !! Mass from original emission cell
        as%cellDispGrd_Val(iDescr,iCell) = fTmp
        fMassStored = fMassStored + fTmp
      end do
      as%cellDispGrd_fx(iCell) = tmpCellDispGrd_fx(iCell) / fMassStored
      as%cellDispGrd_fy(iCell) = tmpCellDispGrd_fy(iCell) / fMassStored
#ifdef DEBUG     
      if  (as%cellDispGrd_fx(iCell) >= nxOut+0.5 .or. as%cellDispGrd_fy(iCell) >=  nyOut+0.5 .or.  &
         & as%cellDispGrd_fx(iCell) <= 0.5       .or. as%cellDispGrd_fy(iCell) <= 0.5 ) then 
        call msg("Gotcha! asrc outside the dispersion grid! fx, fy, fMass", &
           & (/as%cellDispGrd_fx(iCell), as%cellDispGrd_fy(iCell), fMassStored/))
        call report(as)
        call set_error("Trouble with src reprojection", "project_a_src_second_grd")
        return
      endif
#endif
    end do  ! iCell

    call free_work_array(iCellNbrs)
    call free_work_array(tmpCellDispGrd_fx)
    call free_work_array(tmpCellDispGrd_fy)
    
    !
    ! Vertical distribution of area source can be dynamic, then reprojection is irrelevant
    !
    if(.not. as%ifMultiLevelFixed)return

    !
    ! Now re-project the vertical grid of the source to the given vertical
    ! Since vertical distributions are very poorly known and crude, we do not need
    ! any precise meteo-dependent projection, crude will do the job
    !
    as%vertLevsDispVert = vert_disp
    allocate(as%levFractDispVert(fu_NbrOfLevels(vert_disp)),as%fzDisp(fu_NbrOfLevels(vert_disp)), stat=i)
    if(fu_fails(i==0,'Failed to allocate dispersion-vertical level fractions', &
                                                                    & sub_name))return
    as%levFractDispVert(:) = 0.0
    as%fzDisp(:) = 0.0

    call reproject_verticals(as%vertLevs, as%levFraction, &   ! vertical from, fractions from
                           & vert_proj, as%levFractDispVert, &     ! vertical to, fractions to
                           & as%fzDisp, as%nLevsDispVert, & ! mass centres, number of non-zero levels
                           & ifMassCentreInRelUnit = .true.)

    if(len_trim(as%sector_nm) > 0)then
      sTrTmp= trim(as%src_nm)//'_'//trim(as%sector_nm)
    else 
      sTrTmp= as%src_nm
    endif
    iThread = 0
    !$  iThread = omp_get_thread_num()

    if (as%nLevsDispVert>0) then  !If things go bad -- can be whatever
       fTmp=sum(as%levFractDispVert(1:as%nLevsDispVert))
    else
       fTmp=real_missing
    endif

    !!
    !! FIXME  Broken in new voima
    !! 
#ifndef VOIMA_GNU_BUG
    if(ifSimpleReprojection)then
      write (str1Tmp, '(A,A30,A,I10,A,F10.6,A,I3,A,I3)') 'Same pole re-projecting area source: ', trim(strTmp(1:30)), &
            & ', nCellsDispGrd=', as%nCellsDispGrd, ', vertFracSum=', fTmp, &
            & ", nDescr=", as%nDescriptors, &
            & ", iThread=", ithread
!      call msg('Same pole re-projecting area source: ' // trim(strTmp) // &
!             & ', nCellsDispGrd='//trim(fu_str(as%nCellsDispGrd)) // &
!             & ', vertFracSum=' //trim(fu_str(fTmp)) // &
!             & ", nDescr=" //trim( fu_str(as%nDescriptors) )// &
!             & ", iThread=" // trim(fu_str(ithread)))
    else
      write (str1Tmp, '(A,A30,A,I8,A,I8,A,F8.6,A,I3,A,I3)') 'Same pole re-projecting area source: ', trim(strTmp(1:30)), &
            &  ': nSmall=', nSmall, &
            & ', nCellsDispGrd=', as%nCellsDispGrd, ', vertFracSum=', fTmp, &
            & ", nDescr=", as%nDescriptors, &
            & ", iThread=", ithread
!      call msg('Re-projecting area source: ' //trim(strTmp) // &
!             & ': nSmall=' //trim( fu_str(nSmall) ) // &
!             & ', nCellsDispGrd='//trim(fu_str(as%nCellsDispGrd)) // &
!             & ', vertFracSum=' //trim(fu_str(fTmp)) // &
!             & ", nDescr=" //trim( fu_str(as%nDescriptors) )// &
!             & ", iThread=" // trim(fu_str(ithread)))
    endif
    call msg(str1Tmp)
#endif    



    !do icell = 1, nz_dispersion
    !  call msg('Source, level: ' // trim(as%src_nm) // '-' // trim(as%sector_nm), icell)
    !  call msg('fzdisp, vertFract', as%fzDisp(icell), as%levFractDispVert(icell))
    !end do
    
    deallocate(as%cell_Val)
    deallocate(as%cell_fX)
    deallocate(as%cell_fY)
    nullify(as%cell_Val)
    nullify(as%cell_fX)
    nullify(as%cell_fY)

  end subroutine project_a_src_second_grd


  !*******************************************************************

  subroutine init_a_src_TZ_index(a_src, interpCoefMeteo2DispHoriz, &
                 &TZidxFld, disp_grid)

    ! Checks if the source uses local or solar time and
    ! Assigns time_zone index for each cell
    implicit none
    type(silam_area_source), intent(inout) :: a_src
    type(THorizInterpStruct), intent(in) ::  interpCoefMeteo2DispHoriz
    type(silja_field), pointer :: TZidxFld
    type(silja_grid), intent(in) :: disp_grid
    
    real :: geolon
    integer :: iCell, indexMeteo
    REAL, DIMENSION(:), POINTER :: tz_data
    
    if (a_src%ifFieldGiven) return
    if (a_src%nCellsDispGrd < 1) return
    if (a_src%tz_index /= solarTimeIndex .and. a_src%tz_index /= LocalTimeIndex) return

    !Need tz_index array in the source
    allocate(a_src%CellTZindex(a_src%nCellsDispGrd), stat = iCell)
    if (iCell /= 0) then
      call msg('Failed allocation status:', iCell)
      call set_error('Failed allocation', 'init_a_src_TZ_index')
      return
    endif

    if (a_src%tz_index == solarTimeIndex) then
          do iCell = 1, a_src%nCellsDispGrd
             geolon = fu_lon_geographical_from_grid( &
                             & a_src%cellDispGrd_fx(iCell), &
                             & a_src%cellDispGrd_fy(iCell), &
                             & disp_grid)
              !Map longitude to 1:24 index
              a_src%CellTZindex(iCell) =  &
              !                 & solar_TZ_index(nint(modulo(geolon+7.5,360.)/15 + 0.5))
                 & solar_TZ_index(modulo(nint(geolon/15),24) + 1 )
           enddo
     else !Local time
        tz_data => fu_grid_data(TZidxFld)
        do iCell = 1, a_src%nCellsDispGrd
           indexMeteo = fu_grid_index(nx_meteo, nint(a_src%cellDispGrd_fx(iCell)), &
                     &  nint(a_src%cellDispGrd_fy(iCell)), interpCoefMeteo2DispHoriz)
           a_src%CellTZindex(iCell) = nint(tz_data(indexMeteo))
        enddo
     endif

  end subroutine init_a_src_TZ_index


  !*******************************************************************

  subroutine prepare_src_vert_params_a_src(a_src, vert_disp, vert_proj)
    !
    ! In case of dynamic vertical, projects it to the given vertical and stores in the 
    ! time params structure
    !
    implicit none

    ! Imported parameter
    type(silam_area_source), intent(inout) :: a_src
    type(silam_vertical), intent(in) :: vert_disp, vert_proj ! projection vertical in meters

    ! Local variable
    integer :: indTime, iz, nz
    real :: fTmp, fTopInd, fBottomInd
    type(silam_vertical) :: vertTmp
    !
    ! If emission is into the whole ABL, prepare the fldEmis3d and make sure that the source 
    ! itself is a single layer
    !
    if(a_src%ems2wholeABL)then
      if(.not. allocated(a_src%fldEmis3d)) allocate(a_src%fldEmis3d)
      return
    endif
    !
    ! vertical structure is enough in case of field
    !
    if(a_src%ifFieldGiven)return 
    !
    ! If a list of cells, have to do it hard way
    ! First, prepare the mass unit for the cocktail conversion
    !
    nz = 0 !fu_NbrOfLevels(vert_disp)
    a_src%nLevsDispVert = 0

    do indTime = 1, size(a_src%params)

      ! Project the dynamic vertical to the given one
      !
      if(a_src%ifMultiLevelFixed)then
        !
        ! Nothing to project - all is fixed and this layer is ignored
        !
        if(allocated(a_src%params(indTime)%levFractDispVert)) deallocate(a_src%params(indTime)%levFractDispVert)
        if(allocated(a_src%params(indTime)%fzDisp)) deallocate(a_src%params(indTime)%fzDisp)
      else
        !
        ! Create the temporary vertical from a single layer of the slot in order to project the 
        ! given vertical onto it
        !
        allocate(a_src%params(indTime)%levFractDispVert(fu_NbrOfLevels(vert_disp)), &
               & a_src%params(indTime)%fzDisp(fu_NbrOfLevels(vert_disp)))
        call set_vertical(a_src%params(indTime)%layerDynamic, vertTmp)
        if(error)return

        a_src%params(indTime)%levFractDispVert = 0.0
        a_src%params(indTime)%fzDisp = 0.0
        
        call reproject_verticals(vertTmp, (/1.0/), &   ! vertical from, fractions from
                                & vert_proj, a_src%params(indTime)%levFractDispVert, &     ! vertical to, fractions to
                                & a_src%params(indTime)%fzDisp, nz, & ! mass centres, number of non-zero levels
                                & ifMassCentreInRelUnit = .true.)
        a_src%nLevsDispVert = max(nz, a_src%nLevsDispVert)

!    call reproject_verticals(as%vertLevs, as%levFraction, &   ! vertical from, fractions from
!                           & vert_proj, as%levFractDispVert, &     ! vertical to, fractions to
!                           & as%fzDisp, as%nLevsDispVert, & ! mass centres, number of non-zero levels
!                           & ifMassCentreInRelUnit = .true.)

        
!!        do iz = 1, nz
!!          !
!!          ! Project the layer of the given vertical to the source vertical 
!!          !
!!          fTopInd = fu_project_level_crude(fu_upper_boundary_of_layer(fu_level(vert_proj,iz)), &
!!                                         & vertTmp)
!!          fBottomInd = fu_project_level_crude(fu_lower_boundary_of_layer(fu_level(vert_proj,iz)), &
!!                                            & vertTmp)
!!          if(error)return
!!
!!          ! Having upper and lower indices of the current given vertical layer with regard to 
!!          ! original source layer, get the fractions that are within this layer
!!          ! 
!!          ! Find out the overlap between the levels (temporary use of xNew variable)
!!          !
!!          fTmp = min(fTopInd,1.5) - max(fBottomInd,0.5)
!!
!!          if(fTmp > 0.0)then
!!            !
!!            ! Mass fraction is the overlap multiplied with mass fraction into the original layer
!!            !
!!            a_src%params(indTime)%levFractDispVert(iz) = fTmp
!!            !
!!            ! Centre of the emitted mass in given vertical requires opposite projection
!!            ! of source layer to given vertical:
!!            !
!!            fTopInd = fu_project_level_crude(fu_upper_boundary_of_layer(a_src%params(indTime)%layerDynamic), &
!!                                           & vert_proj)
!!            fBottomInd = fu_project_level_crude(fu_lower_boundary_of_layer(a_src%params(indTime)%layerDynamic), &
!!                                              & vert_proj)
!!            ! 1st moment of the emitted mass, in relative units [-0.5, 0.5]
!!            !
!!            a_src%params(indTime)%fzDisp(iz) = (0.5 * (max(fBottomInd,real(iz)-0.5) + &
!!                                                     & min(fTopInd,real(iz)+0.5)) - iz) !* &
!!!                                             & fu_layer_thickness_local_unit(fu_level(vert_disp,iz))
!!
!!          else
!!            a_src%params(indTime)%fzDisp(iz) = 0.
!!            a_src%params(indTime)%levFractDispVert(iz) =0.
!!          endif  ! fTmp > 0 - overlap
!!          
!!        end do   ! nz
      end if  ! if single_level_dynamic
    end do  ! time slots

  end subroutine prepare_src_vert_params_a_src

  
  !**************************************************************************

  subroutine inject_emission_area_source_field(a_src, &
                                             & mapEmis, mapCoordX, mapCoordY, mapCoordZ, &
                                             & met_buf, &
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

    TYPE(silam_area_source), target, INTENT(inout) :: a_src ! inout since binary input file can be reopen
    type(Tmass_map), intent(inout) :: mapEmis, mapCoordX, mapCoordY, mapCoordZ 
    type(Tfield_buffer), pointer :: met_buf
    type(silja_time), intent(in) :: now
    type(silja_interval), intent(in) :: timestep
    logical, intent(in) :: ifSpeciesMoment
    real(r8k), dimension(:),  intent(inout) :: fMassInjected
    type(TVertInterpStruct), intent(in) ::  interpCoefMeteo2DispVert
    type(THorizInterpStruct), intent(in) ::  interpCoefMeteo2DispHoriz
    logical, intent(in) :: ifMetVertInterp
    logical, intent(in) :: ifMetHorizInterp

    ! Local variables
    integer :: iSlot, iLev, nLevs, iDescr, iFile, fs, nSlots, uTmp, ix, iy, iEmis, iMeteo, &
             & iSpecies, nSpecies, iSlotStart, iSlotEnd, iSpeciesEmis, indBLH, nx_emis, ny_emis
    real :: fMassTmp, fCellTotal, fWeightPast, fLevFraction, duration_sec, fBLH
    real, dimension(:), pointer :: fMassTimeCommon, pData, cellArea
    type(silja_time), target :: timeStart, timeEnd, timeTmp
    type(silja_time), pointer :: pTimeStart, pTimeEnd
    type(silja_interval) :: timestep_infile
    real, dimension(3) :: ptrCoord
    character (len=fnlen) :: chTmp
    logical :: ifFound, ifVertInterp
    type(silja_field_id), pointer :: idRequest
    type(silja_field_id) :: idFromFile, idTmp
    type(silja_field), pointer :: pField
    type(silja_3d_field), pointer :: pFldEmis3d
    real, dimension(:), pointer :: pBLH, pEmisMain, pDataLev
    character(len=*), parameter :: sub_name = 'inject_emission_area_source_field'

    if (a_src%nCellsDispGrd < 1) return !! No cells to emit

    if(.not. a_src%ifFieldGiven)then
       call set_error("Tried to inject cell_list area source as field",sub_name)
       return
    endif
    pField =>null()
    idRequest =>null()
    pTimeStart =>null()
    pTimeEnd =>null()
    fMassTimeCommon =>null()
    pData =>null()

    if(a_src%ems2wholeABL)then
      if(fu_fails(fu_NbrOfLevels(a_src%vertLevs) == 1, &
       & 'For whole-ABL injection, emission field must have 1 layer',sub_name))return
      indBLH = fu_index(met_buf, abl_height_m_flag, .true.)
      pBLH => met_buf%p2d(indBLH)%present%ptr
    else
      nullify(pBLH)
    endif
    !
    ! Have to be very careful: starting and ending of the source slots can be pretty close
    ! especially in case of adjoint run. The same is true for time variation coefficients,
    ! which can vary sharply. Therefore, the total mass injected to the grid will be
    ! integrated along time. 
    !

#ifdef DEBUG_SRC
    call msg("inject_emission_area_source_field")
#endif
    ! Source time interval is always positive
    if (fu_sec(timestep)>0) then
      timeStart = now
      timeEnd = now + timestep
    else
      timeStart = now + timestep
      timeEnd = now
    endif

    ! Check if source is enabled...
    if (a_src%params(1)%time > timeStart) timeStart = a_src%params(1)%time
    if (a_src%params(2)%time < timeEnd) timeEnd = a_src%params(2)%time
    if (timeStart >= timeEnd) then
      call msg("Do not inject_emission_area_source_field: out of the interval!")
      return
    endif

    ! If there's template in the filename the times in the source can be invalid, as initialization reads 
    ! only some file. Here update the time list for the current file.
    ! Complete time list has probably been gotten from grads ctl, so only netcdf here
    ! HUOM!! One file must serve the whole emission timestep
    !        File a_src%times must be filled with times from some file already
    nSlots = size(a_src%times)

    if(a_src%FieldFormat(a_src%indFile4Descr(1))%iformat==netcdf_file_flag)then
       if (a_src%times(1) > timeEnd .or. a_src%times(nSlots) < timeStart) then 

         !!File outside of the current emission step -- update time slots
         timeTmp=timeStart + (a_src%times(2)-a_src%times(1))*a_src%slotOffForFile 
            !some hint on in-file timestep
         chTmp = fu_FNm(a_src%chFieldFNmTemplates(a_src%indFile4Descr(1)), &    ! Input file name
                                  & timeTmp, timeTmp, &
                                  & zero_interval, &
                                  & chSource = a_src%src_nm, chSector = a_src%sector_nm)
          call msg("Field source switching to new fname: "//trim(chTmp), iSlot)
          call open_input_file(chTmp, & ! Input file name
                          & a_src%FieldFormat(a_src%indFile4Descr(1)), &
                          & a_src%uField(a_src%indFile4Descr(1)), a_src%uField(a_src%indFile4Descr(1)))
          deallocate(a_src%times) ! Must be allocated
          call timelst_from_netcdf_file(a_src%uField(a_src%indFile4Descr(1)), a_src%times, nSlots, .True., .True.)
          if (fu_fails(.not. error, "After switch", sub_name)) return

       endif
    endif
    !
    ! Find the range of slots to be looked at
    !
    iSlotEnd = nSlots+1 !Out of range
    iSlotStart = -1
    do iSlot = 1, nSlots
      if(a_src%times(iSlot) <= timeStart) iSlotStart = iSlot
      if(a_src%times(iSlot) >= timeEnd)then
        iSlotEnd = iSlot
        exit
      endif
    end do
    if(fu_fails(iSlotStart < iSlotEnd, "iSlotStart >= iSlotEnd, should never happen!", sub_name))return
    if(iSlotEnd > nSlots .or. iSlotStart < 1 )then ! time is outside the source validity range
      call msg_warning('Below period is outside the source file validity range',sub_name)
      call msg('Emission step timeStart =' + fu_str(timeStart) + ', timeEnd =' + fu_str(timeEnd))
      call report(a_src%times)
      call msg('Area source has no emission in this period:')
      call report_area_source(a_src)
      return
    endif

    !---------------------------------------------------------------------------------------------
    !
    ! If the data are defined via field, use interpolation structures and careful time integration
    !
    ! Every time slot has to be read separately from the binary file and reprojected to 
    ! the dispersion grid and vertical.
    ! Note that the times are arranged so that the emission valid from one time to another,
    ! i.e. the number of slots is n+1, with emission startign at the first time slot and ending
    ! at the last one.
    !
    ! The most-frequent case is when the current time slot is the same as at previous time step.
    ! No need to read anything then.
    !
    ! If only one timestep is needed from the file, the necessary field can already be read 
    ! All fields in the same source should have same times, so check only the first
    !
    pField => fu_field_from_3d_field(a_src%pFldEmis(1)%fp,1)  
    idRequest => fu_id(pField) 
    if(.not. (iSlotStart == iSlotEnd-1 .and. &
     & fu_valid_time(idRequest) == a_src%times(iSlotStart+1)))then
      !
      ! Reading needed. Get field from file   
      !
      a_src%fDescrTotals(1:a_src%nDescriptors) = 0.0
      nLevs = fu_NbrOfLevels(a_src%vertLevs)
      fs = fu_number_of_gridpoints(dispersion_grid) ! Sic!  Fields stored in dispersion_grid
      if(iSlotStart /= iSlotEnd-1) fMassTimeCommon => fu_work_array(fs) ! intermediate needed for summation
      if (a_src%ifFluxes) then
         cellArea => fu_work_array(fs)
         do ix=1,fs !!Actually i1d
          cellArea(ix) = fu_cell_size(dispersion_grid,ix)
         enddo
      endif

      if(error)return
        
!    call msg("inject_emission_area_source 4")
      do iSlot = iSlotStart, iSlotEnd-1 
        
        if(iSlotStart == iSlotEnd-1)then
          duration_sec = fu_sec(timeEnd - timeStart)
        else
          if(iSlot == iSlotStart)then
            pTimeStart => timeStart
          else
            pTimeStart => a_src%times(iSlot)
          endif
          if(iSlot == iSlotEnd-1)then
            pTimeEnd => timeEnd
          else
            pTimeEnd => a_src%times(iSlot+1)
          endif
          duration_sec = fu_sec(pTimeEnd - pTimeStart)
        endif  ! is more than one slot
        
        timeTmp = a_src%times(iSlot+a_src%slotOffForFile)  !!Time for filename generation
        !
        ! Read the stuff species-by-species and layer by layer. Also file by file.
        !
        do iLev = 1, nLevs
            
          do iDescr = 1, a_src%nDescriptors
            pField => fu_field_from_3d_field(a_src%pFldEmis(iDescr)%fp,iLev)
            idRequest => fu_id(pField)
            
            if(iSlotStart == iSlotEnd-1)then          ! only one time slot to read
              fMassTimeCommon => fu_grid_data(pField) ! take it to main datapointer
            else
            pData => fu_grid_data(pField)    ! these data will be summed-up to intermediate
            if(iSlot == iSlotStart) pData(1:fs) = 0  ! prepare for summation below
            endif
            if(error)return

            call set_valid_time(idRequest, a_src%times(iSlot+1))
            call set_accumulation_length(idRequest, a_src%times(iSlot+1)- a_src%times(iSlot))
            call set_grid(idRequest, a_src%grid) !! Grid from file. Should be reset back 
            if(error)return
            chTmp = fu_FNm(a_src%chFieldFNmTemplates(a_src%indFile4Descr(iDescr)), &    ! Input file name
                          & timeTmp, timeTmp, zero_interval, &
                          & chSource = a_src%src_nm, chSector = a_src%sector_nm) 
            !
            ! Use get_input_field sub providng it with the iBinary pointer to already open file
            !
            if (a_src%ifFluxes) call set_quantity(idRequest, emission_flux_flag)
!              fMassTimeCommon(1:fs) = 0.0   not needed, overwritten anyway

            call get_input_field(chTmp, &
                               & a_src%FieldFormat(a_src%indFile4Descr(iDescr)), &   ! Input file format
                               & idRequest, &         ! The stuff to search
                               & fMassTimeCommon, &   ! data_requested: either field or an intermediate
                               & dispersion_grid, &    ! storage grid. Note, as%grid /= dispersion_grid == storage_grid
                               & .false., setZero, &  ! ifAcceptSameMonth, iOutside, &
                               & idFromFile, &        ! what actually has been read
                               & 5, &                 ! iAccuracy, &
                               & wdr_missing, &
                               & 0., &        ! fFillValue_, &
                               & a_src%uField(a_src%indFile4Descr(iDescr))) ! iBinary of input file

                !!! No totals here! Only fraction of mass might go to the
                !Emission massmap

#ifdef DEBUG_SRC
            if (.not. all(fMassTimeCommon(1:fs)>=0)) then
!              uTmp = fu_next_free_unit()
!              open(unit=uTmp, file='ems_debug.grads',form='binary',recl=fs*4,access='direct')
!              write(uTmp, rec=1) fMassTimeCommon(1:fs)
!              close(uTmp)
              call msg('idRequest')
              call report(idRequest)
              call msg('idFromFile')
              call report(idFromFile)
              call set_error("Negative emission field from file","inject_emission_area_source")
            endif
#endif
            call set_grid(idRequest, dispersion_grid) !! Field is here, reset grid of the stored ID

            !
            ! The array is overwritten when reading, so we cannot have several binaries
            ! contributing to the same descriptor. This is indeed an unlikely situation,
            ! so let's forbid it. Otherwise an extra temporary array is needed
            if(defined(idFromFile))then  ! got something?
              !
              ! if one slot needed, all is read, just scale and go. Otherwise have to sum-up
              !
              if(iSlotStart == iSlotEnd-1)then
                fMassTimeCommon(1:fs) = fMassTimeCommon(1:fs)  * duration_sec
              else
                pData(1:fs) = pData(1:fs) + fMassTimeCommon(1:fs) * duration_sec
              endif
            else
              call msg("")
              call msg("")
              call msg("get_input_field:"+a_src%chFieldFNms(a_src%indFile4Descr(iDescr)))
              call msg("get_input_field:"+a_src%FieldFormat(a_src%indFile4Descr(iDescr))%title )
              call msg("Requested")
              call report(idRequest)
              call msg("Read")
              call report(idFromFile)
              call msg("fMassTimeCommon",fMassTimeCommon(1:100))
              call msg("Source_vertical")
              call report(a_src%vertLevs,.True.)
              call msg("a_src%uField(a_src%indFile4Descr(iDescr)):", a_src%uField(a_src%indFile4Descr(iDescr)))
              call set_error('Did not find any data from source:' + a_src%src_nm, &
                           & 'inject_emission_area_source')
              return
            endif    ! defined idFromFile
          end do  ! iDescriptor
        enddo  ! iLev
      end do  ! time slots
!      call msg("Last ID read from file:")
!      call report(idFromFile)

      if (a_src%ifFluxes) then !!Convert fluxes to rates
        do iLev = 1, nLevs
          do iDescr = 1, a_src%nDescriptors
            pField => fu_field_from_3d_field(a_src%pFldEmis(iDescr)%fp,iLev)
            pData => fu_grid_data(pField)    ! these data will be summed-up to intermediate

            idRequest => fu_id(pField)
            call set_quantity(idRequest, emission_intensity_flag)
            pData(1:fs) = pData(1:fs) * cellArea(1:fs)
          enddo
        enddo
        call free_work_array(cellArea)
      endif
      if(iSlotStart /= iSlotEnd-1) call free_work_array(fMassTimeCommon)

    endif  ! if reading is needed
    !
    ! Source emission array is filled-in. Distribute 
    !
    do iDescr = 1, a_src%nDescriptors
      
      if(a_src%ems2wholeABL)then
        !
        ! Since the actual injection is complicated due to various interpolations and species mapping
        ! just make a temporary 3D field and distribute the emission over it. Then 
        ! the vertical interpolation can be disabled in the big call below
        !
        ! Fill layer by layer adding the next layer if the previous one had something
        !
        ! template for 3d is from the main emission array
        !
        pFldEmis3d => a_src%fldEmis3d
        idTmp = fu_id(fu_field_from_3d_field(a_src%pFldEmis(iDescr)%fp, 1))
        call grid_dimensions(fu_grid(idTmp), nx_emis, ny_emis)
        pEmisMain => fu_grid_data(fu_field_from_3d_field(a_src%pFldEmis(iDescr)%fp, 1))
        !
        ! Scan the layers starting from surface
        !
        do iLev = 1, nz_dispersion
          ! Need to add next layer?
          if(fu_number_of_fields(pFldEmis3d) < iLev)then
            !
            ! next level
            !
            call set_level(idTmp, fu_level(dispersion_vertical, iLev))
            !
            ! create next 2d field, space for data will be allocated inside
            !
            pDataLev => fu_work_array(nx_emis * ny_emis)
            call create_field_in_field_3d(pFldEmis3d, idTmp, pDataLev)
            call free_work_array(pDataLev)
          endif  ! if add next layer
          !
          ! get the pointer to the field data of this layer
          !
          pDataLev => fu_grid_data(fu_field_from_3d_field(pFldEmis3d, iLev))
          pDataLev(:) = 0.0
          !
          ! Scan the horizontal grid computing emission into this layer
          !
          ifFound = .false.
          do iy = 1, ny_dispersion
            do ix = 1, nx_dispersion
              iMeteo =  fu_grid_index(nx_meteo, ix, iy, interpCoefMeteo2DispHoriz)
              iEmis = fu_grid_index(nx_emis, ix, iy, a_src%pHorizIS)
              fBLH = met_buf%p2d(indBLH)%present%ptr(iMeteo)
              if(fBLH > disp_layer_top_m(iLev-1))then
                pDataLev(iEmis) = pEmisMain(iEmis) * &
                                & (min(fBLH, disp_layer_top_m(iLev)) - disp_layer_top_m(iLev-1)) / fBLH
                ifFound = .true.   ! this layer is useful => willl pick next one too
              endif  ! BLH > lev-1
            enddo   ! ix
          enddo   ! iy
          if(.not. ifFound)exit
        enddo  ! iLev
        ifVertInterp = .false.
      else
        ! no ABL dependence
        pFldEmis3d => a_src%pFldEmis(iDescr)%fp
        ifVertInterp = fu_true(a_src%ifVertInterp)
      endif  ! ems2wholeABL
      !
      ! Now interpolate the data with spltitting to species and scaling
      !
#ifdef DEBUG_SRC
    if (fu_true(a_src%ifHorizInterp)) then 
      call set_error("No more field emissions in source grid!", sub_name)
      return
    endif
#endif

      call add_field_to_mass_map_interp(pFldEmis3d, &      ! pointer to valid field_3d
                                      & a_src%pEmisSpeciesMapping(:,iDescr), & ! mapping to emis massmap
                                      & a_src%fDescr2SpeciesUnit(:,iDescr), &  ! fMassFractionTo, 
                                      & a_src%nSpeciesInDescr(iDescr), &       ! nSpeciesTo, &
                                      & fMassInjected, &   ! total amount, inout
                                      & a_src%params(1)%rate_descr_unit(iDescr), &   ! rate scaling
                                      & .False., ifVertInterp, &
                                      & HorizInterpStruct_missing, a_src%pVertIS, &
                                      & a_src%id_nbr, &
                                      & mapEmis, mapCoordX, mapCoordY, mapCoordZ)

       if (any(fMassInjected < 0.)) then 
         call ooops("Neg mass injected")
       endif
!!!!!        !
!!!!!        ! Formality but useful: total mass sent to emission mass map. Note that if we do not
!!!!!        ! read the field the totals still stay
!!!!!        !
!!!!!              !
!!!!!              a_src%fDescrTotals(iDescr) = a_src%fDescrTotals(iDescr) + &
!!!!!                                         & sum(fMassTimeCommon(1:fs)) * duration_sec
!!!!!
#ifdef DEBUG_SRC

      do iSpecies = 1, a_src%nSpeciesInDescr(iDescr)
          iSpeciesEmis = a_src%pEmisSpeciesMapping(iSpecies,iDescr)
!!!!!          fMassInjected(iSpeciesEmis) = fMassInjected(iSpeciesEmis) + &
!!!!!                                                 & a_src%fDescrTotals(iDescr) * &
!!!!!                 & a_src%params(1)%rate_descr_unit(iDescr) * &
!!!!!                 & a_src%fDescr2SpeciesUnit(iSpecies,iDescr)
!!!!!          
          call msg('Time-step emission of:' + fu_substance_name(mapEmis%species(iSpeciesEmis)), & 
                 & fMassInjected(iSpeciesEmis))
      end do  ! iSpecies
#endif      
    end do  ! descriptors

!call msg('Size of area source and field_3d:',sizeof(a_src),sizeof(a_src%pFldEmis))
!call msg('size of horiz & vert structure:',sizeof(a_src%pHorizIS),sizeof(a_src%pVertIS))
!call msg('Size of id_nbr',sizeof(a_src%id_nbr))

  end subroutine inject_emission_area_source_field


  !*******************************************************************************************

  subroutine count_emission_area_source_cellist(a_src, &
                                       & now, timestep, &
                                       & fMassTimeCommon, &
                                       & islotStart, iSlotEnd)
    ! Calculates fMassTimeCommon -- common part of emission 
    ! determined from par_str before applying monthly, daily and weekly caoefficients

    implicit none

    TYPE(silam_area_source), INTENT(inout) :: a_src ! inout since binary input file can be reopen
    type(silja_time), intent(in) :: now
    type(silja_interval), intent(in) :: timestep
    real, dimension(:), intent(out) :: fMassTimeCommon
    integer, intent(out) ::  iSlotStart, iSlotEnd

    ! Local variables
    integer :: iSlot    
    real :: fMassTmp,  duration_sec, fTmp

    type(silja_time), target :: timeStart, timeEnd
    type(silja_time), pointer :: pTimeStart, pTimeEnd
    type(silja_time) :: StartOfWeek
    integer(4) :: iMyOffsetSec, iSecondsSinceStartOfWeek, iMon, iWeekDay, iHour

    real, dimension(3) :: ptrCoord
    logical :: ifFound
    !    integer, save :: iSlotPrev = int_missing  ! Time index in data file
    integer :: iSlotPrev  ! Time index in data file

    
    ! Failsafe values in case we return early
    iSlotStart = 1
    iSlotEnd = 0
    pTimeStart =>null()
    pTimeEnd =>null()

    if(a_src%ifFieldGiven)then
       call set_error("Tried to inject Field area source as cell_list","count_emission_area_source_cellist")
       return
    endif

    iSlotPrev = int_missing  ! Time index in data file

    !FIXME Seems to be wrong!! Wrong time interval is chosen in case of adjoint
    timeStart = now
    timeEnd = timeStart + fu_abs(timestep)

    ! Check that this source is active during this time period and cut the start-end time if needed
    !
    if(timeEnd <= a_src%times(1) .or. timeStart >= a_src%times(size(a_src%times))) return

    if(timeStart < a_src%times(1)) timeStart = a_src%times(1)
    if(timeEnd > a_src%times(size(a_src%times))) timeEnd = a_src%times(size(a_src%times))
    if(error)return
    !
    ! Find the range of slots to be looked at
    !
    do iSlot = 1, size(a_src%times)
      if(a_src%times(iSlot) <= timeStart) iSlotStart = iSlot
      if(a_src%times(iSlot) >= timeEnd)then
        iSlotEnd = iSlot
        exit
      endif
    end do
    if(iSlotStart >= iSlotEnd)return  ! nothing, for whatever reasons

    if(error)return

    !
    ! Rate has to be integrated over time
    !
    call total_from_a_src_descr_unit(a_src, .true. , &   ! if integrate rates only
                                   & fMassTimeCommon, &
                                   & timeStart, timeEnd - timeStart)
    if(error)return
  
  end  subroutine count_emission_area_source_cellist

  !*******************************************************************************************

  subroutine inject_emission_area_source_cellist(a_src, &
                                       & mapEmis, mapCoordX, mapCoordY, mapCoordZ, &
                                       & met_buf, &
                                       & now, &
                                       & ifSpeciesMoment, &
                                       & fMassInjected, &
                                       & interpCoefMeteo2DispHoriz, ifMetHorizInterp, &
                                       & interpCoefMeteo2DispVert, ifMetVertInterp, &
                                       & fMassTimeCommon,  iSlotStart, iSlotEnd, iThread, nThreads)
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

    TYPE(silam_area_source), INTENT(inout) :: a_src ! inout since binary input file can be reopen
    type(Tmass_map),  INTENT(inout) :: mapEmis, mapCoordX, mapCoordY, mapCoordZ
    type(Tfield_buffer), pointer :: met_buf
    type(silja_time), intent(in) :: now
    logical, intent(in) :: ifSpeciesMoment
    real(r8k), dimension(:),  intent(inout) :: fMassInjected
    type(TVertInterpStruct), intent(in) ::  interpCoefMeteo2DispVert
    type(THorizInterpStruct), intent(in) ::  interpCoefMeteo2DispHoriz
    logical, intent(in) :: ifMetVertInterp
    logical, intent(in) :: ifMetHorizInterp
    real, dimension(:), intent(in) :: fMassTimeCommon !(1:a_src%nDescriptors)
    integer, intent(in) ::   iSlotStart, iSlotEnd, iThread, nThreads

    ! Local variables
    integer :: iSlot, iLev, ix, iy, iSrc, iCell, nLevs, iCellX, iCellY, iDescr, iFile, fs, &
             & iSpecies, nSpecies,  iSpeciesEmis
    real :: fMassTmp, fCellTotal, fWeightPast, fLevFraction, duration_sec, fTmp, fMetScaling
    real, pointer :: fPtr 

    type(silja_time), target :: timeStart, timeEnd
    type(silja_time) :: StartOfWeek
    integer(4) :: iMyOffsetSec, iSecondsSinceStartOfWeek, iMon, iWeekDay, iHour, iTmp

    real, dimension(3) :: ptrCoord
    logical :: ifFound, ifCellwiseTz
    !    integer, save :: iSlotPrev = int_missing  ! Time index in data file
    integer :: iSlotPrev  ! Time index in data file

    !---------------------------------------------------------------------------------------------
    !
    ! If the data are defined via field, use interpolation structures and careful time integration
    !
    if (a_src%nCellsDispGrd < 1) return !! No cells to emit

    if(a_src%ifFieldGiven)then
       call set_error("Tried to inject Field area source as cell_list","inject_emission_area_source_cellist")
       return
    endif
    fPtr => null()

    fMetScaling = 1. ! No meteo dependence by default

    ! Do we have any emission this month at all???
    iMon = fu_mon(now) 
    if (all(a_src%indMonth(1:a_src%nDescriptors,iMon) == 0)) return

    if(iSlotStart >= iSlotEnd)return  ! nothing, for whatever reasons
    fWeightPast = (a_src%params(iSlotEnd)%time - now) / &
                & (a_src%params(iSlotEnd)%time - a_src%params(iSlotStart)%time)

    ! Might be not needed, but do it once....
    StartOfWeek = fu_Start_Of_Day_utc(now) - fu_set_interval_d(real(fu_weekday(now)-1))
    iSecondsSinceStartOfWeek = nint(fu_sec(now-StartOfWeek),4)
!call msg("Inject area source Starttime:" + fu_str(now))
!call msg("Injectarea     source StartOfWeek:" + fu_str(StartOfWeek))
!call msg("Injectarea seconds since StartOfWeek:"+ fu_str(iSecondsSinceStartOfWeek))
    if (.not. a_src%ifUseTimeVarCoef) then !All the same
      iHour =1
      iWeekDay = 1
      ifCellwiseTz = .False.
    elseif (a_src%tz_index > 0) then! One timezone for source or constant-rate source
      iMyOffsetSec = TZ_offset(a_src%tz_index)
      iHour = modulo(iSecondsSinceStartOfWeek + iMyOffsetSec, i_seconds_in_day) / i_seconds_in_hour + 1
      iTmp = TZ_special_date_now(a_src%tz_index)
      if (iTmp == 0) then
        iWeekDay = modulo(iSecondsSinceStartOfWeek + iMyOffsetSec, i_seconds_in_week)/ i_seconds_in_day + 1 
      else
        iWeekDay = iTmp
      endif
      ifCellwiseTz = .False.
    else
      ifCellwiseTz = .True.
    endif

    !
    ! We do the injection layer by layer
    !
    do iLev = 1, nz_dispersion
      !
      ! First do the check for the overlap: speed-up
      !
      if(a_src%ifMultiLevelFixed)then
        if( abs(a_src%levFractDispVert(iLev)) < 1e-5)cycle  ! nothing for this dispersion layer
        ptrCoord(3) = a_src%fzDisp(iLev)
        fLevFraction = a_src%levFractDispVert(iLev)
      else
        ifFound = .false.
        do iSlot = iSlotStart, iSlotEnd - 1
          if((a_src%params(iSlot)%levFractDispVert(iLev) +  &
            & a_src%params(iSlot+1)%levFractDispVert(iLev)) > 1e-5) then
            ifFound = .true.
            exit
          endif
        end do
        if(.not. ifFound) cycle
        ptrCoord(3) = a_src%params(iSlotStart)%fzDisp(iLev) * fWeightPast + &
                    & a_src%params(iSlotEnd)%fzDisp(iLev) * (1. - fWeightPast) 
        fLevFraction = a_src%params(iSlotStart)%levFractDispVert(iLev) * fWeightPast + &
                     & a_src%params(iSlotEnd)%levFractDispVert(iLev) * (1. - fWeightPast)
      endif

      ! Now go cell-wise
      do iCell = 1, a_src%nCellsDispGrd

        iCellY = nint(a_src%cellDispGrd_fy(iCell))

        if ( mod(iCellY,nThreads) /= iThread) cycle  ! Is it our stripe at all???

        iCellX = nint(a_src%cellDispGrd_fx(iCell))

        if(ifCellwiseTz)then
          ! If cell-wise index -- take it
          iTmp = a_src%CellTZindex(iCell)
          iMyOffsetSec = TZ_offset(iTmp)
          iHour = modulo(iSecondsSinceStartOfWeek + iMyOffsetSec, i_seconds_in_day) / i_seconds_in_hour + 1
          iTmp = TZ_special_date_now(iTmp)
          if (iTmp == 0) then
            iWeekDay = modulo(iSecondsSinceStartOfWeek + iMyOffsetSec, i_seconds_in_week)/ i_seconds_in_day + 1 
          else
            iWeekDay = iTmp
          endif
        endif
        fCellTotal = 0.

        ptrCoord(1) = a_src%cellDispGrd_fx(iCell) - iCellX
        ptrCoord(2) = a_src%cellDispGrd_fy(iCell) - iCellY

         
        if (a_src%ifTemprDependent)  fMetScaling = &
              & fu_meteodep_emis_scaling(a_src%rulesMeteoDep, iCellX + nx_dispersion * (iCellY-1))

        !
        ! Scan the whole rectangular filtering the non-emitting
        ! subst-mode combinations via coding array
        !
        mapEmis%ifGridValid(ilev,a_src%id_nbr) = .true.
        mapEmis%ifColumnValid(a_src%id_nbr,iCellX,iCellY) = .true.

        do iDescr = 1, a_src%nDescriptors
          !
          ! All is cell-specific. Time-slot cocktail is used for an overall release rate scaling. 
          ! Composition is dictated by cells descriptor-wise. The fractionation of species in the 
          ! descriptors is taken from descriptors themselves, of course. 
          ! Note that the massTmp is the descriptor mass.
          !
          if(a_src%ifUseTimeVarCoef)then
            fMassTmp = fMassTimeCommon(iDescr) * a_src%cellDispGrd_Val(iDescr,iCell) * &
                     & a_src%indHour(iDescr,iHour) * a_src%indDay(iDescr,iWeekDay) * a_src%indMonth(iDescr,iMon)
          else
            fMassTmp = fMassTimeCommon(iDescr) * a_src%cellDispGrd_Val(iDescr,iCell)
          endif

          fMassTmp = fMassTmp * fMetScaling

#ifdef DEBUG          
          if (fMassTmp < 0 ) then
              call msg("a_src%cellDispGrd_Val(iDescr,iCell)", a_src%cellDispGrd_Val(iDescr,iCell))
              call msg("a_src%indHour(iDescr,iHour",a_src%indHour(iDescr,:))
              call msg("a_src%indDay(iDescr,iWeekDay)",a_src%indDay(iDescr,:))
              call msg("a_src%indMonth(iDescr,iMon)", a_src%indMonth(iDescr,:))
              call msg("fMassTimeCommon(iDescr)",fMassTimeCommon(iDescr))
              call msg("iDescr,iCell, iHour, iWeekDay, iMon", (/iDescr,iCell, iHour, iWeekDay, iMon/))
              call set_error("Gotcha -neg emission: "//fu_str(fMassTmp),"inject_emission_area_source_cellist")
          endif
#endif          
          !
          ! Now we can explore the species from the cocktail descriptor
          !
          do iSpecies = 1, a_src%nSpeciesInDescr(iDescr)
            iSpeciesEmis = a_src%pEmisSpeciesMapping(iSpecies,iDescr)

            fTmp = fMassTmp * fLevFraction * a_src%fDescr2SpeciesUnit(iSpecies,iDescr)
            fCellTotal = fCellTotal + fTmp

            fPtr => mapEmis%arM(iSpeciesEmis,a_src%id_Nbr,iLev,iCellX,iCellY) 
            fPtr = fPtr + fTmp

            fMassInjected(iSpeciesEmis) =  fMassInjected(iSpeciesEmis) + fTmp
            
            if (ifSpeciesMoment) then
               fPtr => mapCoordX%arm(iSpeciesEmis,a_src%id_nbr,ilev, iCellX, iCellY)
               fPtr = fPtr + fTmp * ptrCoord(1)

               fPtr => mapCoordY%arm(iSpeciesEmis, a_src%id_nbr, ilev, iCellX,iCellY)  
               fPtr = fPtr + fTmp * ptrCoord(2)

               fPtr => mapCoordZ%arm(iSpeciesEmis, a_src%id_nbr, ilev, iCellX,iCellY)
               fPtr = fPtr + fTmp * ptrCoord(3)

            end if
          
          end do  ! species inside the descriptor
        end do  ! iDescr
        
        ! Add the moment of the injected masses and the total injected mass
        if (.not. ifSpeciesMoment) then
          fPtr => mapCoordX%arM(1,a_src%id_nbr, iLev, iCellX, iCellY)
          fPtr = fPtr + ptrCoord(1) * fCellTotal

          fPtr => mapCoordY%arM(1,a_src%id_nbr, iLev, iCellX, iCellY)
          fPtr = fPtr + ptrCoord(2) * fCellTotal

          fPtr => mapCoordZ%arM(1,a_src%id_nbr, iLev, iCellX, iCellY)
          fPtr = fPtr + ptrCoord(3) * fCellTotal
        end if

      end do  ! iCell
    end do  ! iLev dispersion


  end subroutine inject_emission_area_source_cellist


  !***************************************************************************************

  subroutine check_cells(a_src)
    implicit none
    TYPE(silam_area_source), INTENT(in) :: a_src
    integer :: iCell, iDescr

    do iCell = 1, a_src%nCells
      do iDescr = 1, a_src%nDescriptors
        if(a_src%cell_val(iDescr,iCell) < 0)then
          call msg(fu_connect_strings('Cell<0:',a_src%src_nm,',',a_src%sector_nm), &
                 & iCell, a_src%cell_val(iDescr,iCell))
          call msg('Descriptor:',iDescr)
          call set_error('Cell < 0','check_cells')
        endif
      end do
    end do

  end subroutine check_cells



  !*******************************************************************************
  !
  !  A bit of encapsulation
  !
  !*******************************************************************************

  !=========================================================================  
  logical function fu_area_source_defined(a_src)result(def)
    implicit none 
    type(silam_area_source), intent(in) :: a_src
    
    def = fu_true(a_src%defined)
    
  end function fu_area_source_defined

  !=========================================================================  
  subroutine set_area_source_undefined(a_src)
    implicit none 
    type(silam_area_source), intent(inout) :: a_src
    if (associated(a_src%params)) deallocate(a_src%params)
    a_src%nCells = 0
    a_src%defined = silja_false
  end subroutine set_area_source_undefined


  !=========================================================================  
  function fu_SlotTime_of_area_source(a_src, iSlot)
    implicit none
    type(silja_time) :: fu_SlotTime_of_area_source
    TYPE(silam_area_source), INTENT(in) :: a_src
    integer, intent (in) :: iSlot
    if(iSlot > size(a_src%params))then
      call set_error('Too big slot number','fu_SlotTime_of_area_source')
      return
    endif
    fu_SlotTime_of_area_source = a_src%params(iSlot)%time
  end function fu_SlotTime_of_area_source

  !=========================================================================  
  FUNCTION fu_grid_of_area_source(a_src)
    IMPLICIT NONE
    TYPE(silja_grid) :: fu_grid_of_area_source
    TYPE(silam_area_source), INTENT(in) :: a_src
    fu_grid_of_area_source = a_src%grid
  END FUNCTION fu_grid_of_area_source

  !=========================================================================  
  FUNCTION fu_area_source_start_time(a_src)
    IMPLICIT NONE
    TYPE(silja_time) :: fu_area_source_start_time
    TYPE(silam_area_source), INTENT(in) :: a_src
    fu_area_source_start_time = a_src%params(1)%time
  END FUNCTION fu_area_source_start_time

  !=========================================================================  
  FUNCTION fu_area_source_end_time(a_src)
    IMPLICIT NONE
    TYPE(silja_time) :: fu_area_source_end_time
    TYPE(silam_area_source), INTENT(in) :: a_src
    fu_area_source_end_time = a_src%params(size(a_src%params))%time
  END FUNCTION fu_area_source_end_time

  !=========================================================================  
  FUNCTION fu_area_source_duration(a_src)
    IMPLICIT NONE
    TYPE(silja_interval) :: fu_area_source_duration
    TYPE(silam_area_source), INTENT(in) :: a_src
    fu_area_source_duration = a_src%params(size(a_src%params))%time - a_src%params(1)%time
  END FUNCTION fu_area_source_duration

  !=========================================================================  
  FUNCTION fu_NbrOfTimeSlots_of_a_src(a_src)
    IMPLICIT NONE
    integer :: fu_NbrOfTimeSlots_of_a_src
    TYPE(silam_area_source), INTENT(in) :: a_src
    fu_NbrOfTimeSlots_of_a_src = size(a_src%params)
  END FUNCTION fu_NbrOfTimeSlots_of_a_src

  !=========================================================================  
  logical function fu_ifMultiLevelFixed_a_src(a_src)
    implicit none
    type(silam_area_source), intent(in) :: a_src
    fu_ifMultiLevelFixed_a_src = a_src%ifMultiLevelFixed
  end function fu_ifMultiLevelFixed_a_src

  !=========================================================================  
  integer function fu_n_cells_a_Src(a_src)
    implicit none
    type(silam_area_source), intent(in) :: a_Src
    fu_n_cells_a_Src = a_src%nCells
  end function fu_n_cells_a_Src

  !=========================================================================  
  logical function fu_useTimeVarCoef_a_src(a_src)
    implicit none
    type(silam_area_source), intent(in) :: a_src
    fu_useTimeVarCoef_a_src = a_src%ifUseTimeVarCoef
  end function fu_useTimeVarCoef_a_src

  !=========================================================================  
  function  fu_name_a_src(a_src) result(nm)
    implicit none
    character(len=clen) :: nm
    type(silam_area_source), intent(in) :: a_src
    nm = a_src%src_nm
  end function fu_name_a_src

  !=========================================================================  
  function  fu_sector_a_src(a_src) result(nm)
    implicit none
    character(len=clen) :: nm
    type(silam_area_source), intent(in) :: a_src
    nm = a_src%sector_nm
  end function fu_sector_a_src

  !=========================================================================  
  integer function fu_source_nbr_of_a_src(a_src)
    !
    ! Returns the source number. Reason: all sources are enumerated
    ! sequencially, so that the source can, in fact, be refered by its
    ! number without other information. But so far this number is 
    ! copied to the particles in the pollution cloud. 
    !
    implicit none
    type(silam_area_source), intent(in) :: a_src
    if(.not. (a_src%defined == silja_false))then
      fu_source_nbr_of_a_src = a_src%src_nbr
    else
      fu_source_nbr_of_a_src = int_missing
      call set_error('Undefined source given','fu_source_nbr_of_area_source')
      return
    endif
    
  end function fu_source_nbr_of_a_src

  !=========================================================================  
  function fu_cocktail_descr_of_a_src (a_src) result(descrLst)
    implicit none
    type(Tcocktail_descr), dimension(:), pointer :: descrLst
    type(silam_area_source), intent(in), target :: a_src
    descrLst => a_src%cocktail_descr_lst
  end function fu_cocktail_descr_of_a_src

  !=========================================================================
  integer function fu_meteodep_model_a_src(a_src)
    implicit none
    type(silam_area_source), intent(in) :: a_src
    if(a_src%ifTemprDependent)then
      fu_meteodep_model_a_src = fu_meteodep_model(a_src%rulesMeteoDep)
    else
      fu_meteodep_model_a_src = int_missing
    endif
  end function fu_meteodep_model_a_src
  
  
  ! ***************************************************************

  SUBROUTINE report_area_source(as)

    IMPLICIT NONE

    ! Imported parameters with intent IN:
    TYPE(silam_area_source), INTENT(in) :: as

    integer :: i, j
    real, dimension(:), pointer :: fTmp 

    fTmp => fu_work_array(as%nDescriptors)

    call msg('')
    call msg('********************** AREA source **************************')
    if (defined(as)) then
      call msg('Area source:' + as%src_nm)
    else
      call msg('Area source UNDEFINED:' + as%src_nm)
      return
    endif
    call msg('Grid: ')
    call report(as%grid)
    call msg('start time of release:' + fu_str(as%times(1)))
    call msg('duration of release:' + fu_str(as%times(size(as%times))-as%times(1)))
    call msg('end time of release:' + fu_str(as%times(size(as%times))))
    if(as%ifFieldGiven)then
      call msg('The source is given as time- and vertical-resolving field')
      call msg('Source grid:')
      call report(as%grid)
      call msg('Source vertical:')
      call report(as%vertLevs, .true.)
    else
      call msg(' Parameters at the beginning of the release: ')
      if(as%ifMultiLevelFixed)then
        call msg('Source vertical:')
        call report(as%vertLevs)
      else
        IF (as%params(1)%vertical_distribution == vertically_single_level) THEN
          call msg('single start level: ')
          CALL report(as%params(1)%layerDynamic)
        ELSE
          call msg('bottom & top: ')
          CALL report(as%params(1)%layerDynamic)
        END IF
      endif
      fTmp(1:as%nDescriptors)=0.
      call msg('Number of cells:',as%nCells)
      call msg('Value dimensions (nDescriptors):',as%nDescriptors)
      do i=1,as%nCells
        fTmp(1:as%nDescriptors) = fTmp(1:as%nDescriptors) + as%cell_Val(1:as%nDescriptors,i)
      end do
      do j = 1, as%nDescriptors
        call msg('Starting release rate,' + fu_name(as%cocktail_descr_lst(j)) + ',' + &
                                   & fu_basic_mass_unit(as%cocktail_descr_lst(j)) + '/sec:', &
             & as%params(1)%rate_descr_unit(j) * fTmp(j))
      end do
    endif ! ifFieldGiven
    call msg('********************* end of AREA source **********************')

    call free_work_array(fTmp)

  END SUBROUTINE report_area_source



END MODULE source_terms_area




