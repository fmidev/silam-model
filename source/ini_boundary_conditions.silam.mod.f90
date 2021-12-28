MODULE ini_boundary_conditions
  ! 
  ! Here are tools for managing the initial and boundary conditions from the external
  ! parameters or fields.
  !
  ! Authors: Marje Prank, Mikhail Sofiev, FMI, mikhail.sofiev@fmi.fi
  !
  ! All units: SI
  ! 
  ! All degrees: real numbers in degrees (degrees and parts per
  ! hundred/parts per thousand , NOT in degrees and minutes
  ! 
  ! Language: ANSI standard Fortran 90
  ! 
  ! Modules used:
  use dispersion_server
  use silam_partitioning

  IMPLICIT NONE

  private

  ! Initial conditions
  !
  public set_initial_conditions  ! initialises the MassMap
  public store_restart_file    ! Prepares to the next run

  ! Boundary conditions
  !
  public init_boundary_structures
  public fill_boundary_market
  public boundary_conditions_now
  public fillBoundaryStruct
  public set_missing
 
  public fu_nr_boundary_inFiles
  public fu_shplst 
  public fu_ifHaveInitialConditions
 
  ! Ini & boundary rules stuff
  !
  public set_ini_boundary_rules
  public defined
  public add_ini_boundary_input_needs
  public fu_obstime_interval
  public fu_ifReadBoundaries
  public fu_ifBoundary
  public force_mass_map_cell_mmr
  public update_btypes

  private set_ini_boundary_rules_from_nl
  private fu_obstime_int_bcr
  private if_ini_boundary_rules_defined
  private fu_boundaryStructDefined

  ! INTERFACE
  !
  interface defined
    module procedure fu_boundaryStructDefined
  end interface

  INTERFACE fu_obstime_interval
    MODULE PROCEDURE fu_obstime_int_bcr
  END INTERFACE

  interface set_ini_boundary_rules
    module procedure set_ini_boundary_rules_from_nl
  end interface

  interface defined
    module procedure if_ini_boundary_rules_defined
  end interface

  interface set_missing
     module procedure set_missing_b_rules
  end interface

  !
  ! Structure for storing a single input boundary file, together with BC to apply to
  !
  type Tboundary_file
    character(len=fnlen) :: bHeaderFNm
    type(grads_template), dimension(:), pointer :: bFnameTemplate => null()
    integer, dimension(max_met_files) :: unit_binary
    type(silam_fformat) :: bFileFormat
    logical :: ifStatic
    integer :: fieldType    ! dynamic, static, monthly_climatology
    character(len=6) :: bNames  ! a list of boudnaries to apply to
    logical, dimension(6) :: ifToBoundary ! list of boundaries this header-file serves
    type(silam_vertical) :: vertical
    type(TVertInterpStructPtr), dimension(6) :: interpCoefBndInp2DispVert, interpCoefMet2BndInp
    ! list of field id-s for fields from this header in boundary-stacks
    type(silja_field_id), dimension(max_quantities) :: fieldInStack   
    integer, dimension(max_quantities) :: tSpecies  ! list of transport species this header provides
    integer :: nTSpecies
    type(silja_shopping_list) :: shopLst
    logical :: defined = .false.
  endtype Tboundary_file
  private Tboundary_file

  character, dimension(6), parameter, public :: BoundaryChar = (/'N','S','E','W','T','B'/)
  !
  ! Rules for making the initial and boundary conditions
  !
  type Tini_boundary_rules
    private
    character(len=fnlen), dimension(:), pointer :: chIniFile =>null() ! Initial conditions file
    integer, dimension(:), pointer :: qInitial => null()              ! quantities to initialise
    type(silja_time) :: start_time
    TYPE(silja_interval) :: bInterval
    character, dimension(6) :: bNames = BoundaryChar
    integer, dimension(6) :: bTypes   !Types of boundaries for the whole run,
                            ! reset to the subdomain boundaries by 
    ! a list of files, each with a list of boundaries to apply
    type(Tboundary_file), dimension(:), allocatable :: bFiles 
    logical :: ifBoundary = .false.   !! If any of the boundaries are of dirichlet or ploar type (i.e. requite boundary structures)
    logical :: ifRandomise = .true.        ! if smooth the reprojection aliasing
    type(silja_logical) :: defined = silja_false
  end type Tini_boundary_rules
  public Tini_boundary_rules

  !
  ! The main structure for actual boundary data 
  !
  ! This structure is intended to define a single boundary of the
  ! model cube, although currently only the latitude global closure is implemented. The type
  ! and name of the boundary specify how it should be processed
  ! routine.
  ! 
  type TboundaryStruct
!    private
    integer :: boundaryType = int_missing        ! type of the condition, see below.
    character :: name = char_missing  ! which side of the domain? either N,S,E,W or T, B 
    ! the polar cap, used for spherical computational domains
    integer :: nSp =int_missing , nSrc = int_missing
    real, dimension(:,:,:,:), allocatable :: polarCapMassMom      ! (iSpecies, iSrc, 2, iVertical)
    real, dimension(:), allocatable :: polarCapAirMass            ! Initialize to real_missing, advection will 
                                                                  ! set it to something reasonable....
    real :: polarCapArea
    type(silja_field_3d_pointer), dimension(:,:), pointer :: Fld3dPtr   ! dimension - past/future; 1:nSpTr
    integer :: past_future_switch = 0      !  0 -> [past, future], 1 ->  [future -> past]
    type(THorizInterpStruct), pointer :: meteo2BoundaryInterpHoriz => null()
    type(TVertInterpStruct), pointer :: meteo2BoundaryInterpVert => null()
    type(silam_vertical) :: vertical
    type(silja_wdr), pointer :: bc_wdr
    logical, dimension(:), pointer :: ifBoundSpecies  => null()
    type(silja_logical) :: defined
  end type TboundaryStruct
  public TboundaryStruct

  type TboundaryStructPtr 
    type(TboundaryStruct), pointer :: ptr  => null()
  end type TboundaryStructPtr
  public TboundaryStructPtr

  !
  ! The buffer to the boundary data
  ! The immediate interface between the boundary structure (above) and the advection routine
  ! that takes care of the boundary transport inside the domain.
  ! The above structure cannot be used for this purpose because it is oriented to operations
  ! at every boundary update time, whereas we need the data at each model time step. This
  ! is the function of this buffer.
  ! This is THE ONLY FULLY PUBLIC TYPE here - to speed-up the exchange. pBnd is a small supplementary
  !
  type pBnd
    real, dimension(:,:,:,:), pointer :: ptr
  end type pBnd
  type TboundaryBuffer
    integer, dimension(6) :: iBoundaryType, &        ! replica from boundary structure: zero/polar/...
                           & nBoundarySpecies, &     ! usually less than N transport species
                           & indStep, &              ! 0 or 1 for single-value or full-line boundary
                           & nSrcTrn                 ! number of sources in boundary
    real, dimension(2) :: fPoleCapArea               ! 2: north/south
    real, dimension (:,:), allocatable :: outflowFactor  ! (2, nz) ! Factor limiting outgoing mass from poles
    type(silja_rp_1d), dimension(2) ::  PoleCapAirmass
    real, dimension(:,:,:), allocatable :: wind_to_pole     ! (nx,2, nz)
    real, dimension(:,:,:), allocatable :: poleDryDep    ! (2,iSpecies,iSrc)
    ! index of boundary species in trn mass map (nSpecies,6)
    integer, dimension(:,:), allocatable :: iTransportSpecies, & 
         ! index of transport species in boundary buffer (nTransportSpecies,6)
         & iBoundarySpecies  
    ! Concentrations (per m3)
    real, dimension(:,:,:,:), pointer :: bNorth => null(), bSouth  => null()     ! (nSpecies, nSrc, nx, nz)
    real, dimension(:,:,:,:), pointer :: bEast => null(),  bWest => null()  ! (nSpecies, nSrc, ny, nz)
    real, dimension(:,:,:,:), pointer :: bTop => null(),   bBottom => null() ! (nSpecies, nSrc, nx, ny)
    type(pBnd), dimension(6) :: pBoundaryData
  end type TboundaryBuffer
  public TboundaryBuffer
  public pBnd

  !
  ! Available boundary types and indices of specifis boundaries
  !
  integer, public, parameter :: periodic_boundary_type  = 601, &
                              & polar_boundary_type     = 602, &
                              & dirichlet_boundary_type = 603, &
                              & zero_boundary_type      = 604, &
                              & surface_boundary_type   = 605, &
                              & smpi_comm_boundary_type = 606

  CONTAINS

  !***********************************************************************************
  !***********************************************************************************
  !
  ! INITIAL & BOUNDARY RULES
  !
  !***********************************************************************************
  !***********************************************************************************
  
  subroutine set_ini_boundary_rules_from_nl(nlIni, nlStandardSetup, rulesIniBoundary, start_time)
    !
    ! Sets the boundary rules from the namelist. Does not handle global closure boundaries
    ! as they are not included in the rules.
    !
    implicit none


    ! Imported parameter
    type(Tsilam_namelist), pointer :: nlIni, nlStandardSetup
    type(Tini_boundary_rules) :: rulesIniBoundary
    TYPE(silja_time), intent(in) :: start_time
    ! Local variables
    integer :: iTmp, iStat
    type(Tsilam_nl_item_ptr), dimension(:), pointer :: ptrItems
    character(len = clen) :: chtmp

    character(len=*), parameter :: sub_name="set_ini_boundary_rules_from_nl" 
    !
    ! Stupidity checking
    !
    if(rulesIniBoundary%defined == silja_true)then
      call set_error('Cannot redefine the boundary conditions', sub_name)
      return
    endif
    rulesIniBoundary%defined = silja_false

    if(.not. associated(nlIni))then
      call msg('')
      call msg('')
      call set_error('No initial_and_boundary_conditions list', sub_name)
      return
    endif
    if(.not.defined(nlIni))then
      call set_error('Undefined namelist given', sub_name)
      return
    endif


    !---------------------------------------------------------------------------------
    !
    ! Initial conditions can be of three types:
    ! - initialization of the supplementary quantities
    ! - initialization of concentrations
    ! - a separate creature: restart of the run
    ! In all cases all variables are put into the same namelist in the item "initialize_quantity"
    !
    nullify(ptrItems)
    call get_items(nlIni, 'initialize_quantity', ptrItems, iTmp)
    if(iTmp < 1)then
      call msg('')
      call msg('Empty initialization namelist, start fields are all-zeroes')
      call msg('')
      nullify(rulesIniBoundary%qInitial, rulesIniBoundary%chIniFile)
    else
      !
      ! Note that the list can include concentrations or deposition quantity, which require
      ! substance name for complete definition. In this case the initialization procedure 
      ! will take all substances available in the input files.
      !
      allocate(rulesIniBoundary%qInitial(iTmp), stat=iStat)
      if(fu_fails(iStat == 0, 'Failed to allocate initial quantity list', sub_name))return
      do iStat = 1, iTmp
        rulesIniBoundary%qInitial(iStat) = fu_get_silam_quantity(fu_content(ptrItems(iStat)))
        if(error)then
          call set_error('Strange quantity name:' + fu_content(ptrItems(iStat)), sub_name)
          return
        endif
      end do
      !
      ! Having something to initialize, let's get the input files.
      ! Two types of the initialization are possible:
      ! - from maps stored as GrADS or GRIB files
      ! - concentrations given at a few points - observations-type files
      ! At this level, they are all treated same way: 
      ! initialization_file = GRIB/GRADS/POINT_DATA <fileName>
      !
      call get_items(nlIni, 'initialization_file', ptrItems, iTmp)
      if(iTmp < 1)then
        call set_error('Initialization requested but no files given', sub_name)
        return
      else
        allocate(rulesIniBoundary%chIniFile(iTmp),stat=iStat)
        if(fu_fails(iStat == 0, 'Failed to allocate initial file names', sub_name))return
        do iStat = 1, iTmp
          rulesIniBoundary%chIniFile(iStat) = fu_content(ptrItems(iStat))
        end do
      endif

    endif  ! if any fields to initialize


    !----------------------------------------------------------------------
    !
    ! Boundary conditions for the whole domain parsed here:
    !
    chTmp = fu_str_u_case(fu_content(nlIni,'boundary_type'))
    if (chTmp /= '')  then
      call msg_warning("boundary_type, if_lateral_boundary, if_top_boundary, and if_bottom_boundary depricated", sub_name)
      call msg("Please use lateral_boundary_type and top_boundary_type instead")
      if(chTmp == '' .or. chTmp == 'ZERO')then
        rulesIniBoundary%ifBoundary = .false.         ! Can be set .true. if global, hemispheric or MPI
        rulesIniBoundary%bTypes(1:6) = zero_boundary_type
      elseif(trim(chTmp) == 'DIRICHLET')then
        rulesIniBoundary%ifBoundary = .true.
        rulesIniBoundary%bTypes(1:6) = dirichlet_boundary_type
       else
          call set_error('Unknown boundary type:' + chTmp, sub_name)
          return
        endif

        chTmp = fu_str_u_case(fu_content(nlIni,'if_lateral_boundary'))
        if(trim(chtmp) == 'NO')then
          rulesIniBoundary%bTypes(1:4) = zero_boundary_type
        endif    

        chTmp = fu_str_u_case(fu_content(nlIni,'if_top_boundary'))
        if(trim(chtmp) == 'NO')then
          rulesIniBoundary%bTypes(top_boundary) = zero_boundary_type
        endif

        chTmp = fu_str_u_case(fu_content(nlIni,'if_bottom_boundary'))
        if(trim(chtmp) == 'NO')then
          rulesIniBoundary%bTypes(bottom_boundary) = zero_boundary_type
        endif

    else !!  New BC definition: only lateral_boundary_type and top_boundary_type
       rulesIniBoundary%ifBoundary = .FALSE. 
       chTmp = fu_str_u_case(fu_content(nlIni,'lateral_boundary_type'))
       if (chTmp == 'ZERO') then
         rulesIniBoundary%bTypes(1:4) = zero_boundary_type
       elseif (chTmp == 'PERIODIC') then
         rulesIniBoundary%bTypes(1:4) = periodic_boundary_type
       elseif(trim(chTmp) == 'DIRICHLET')then
         rulesIniBoundary%ifBoundary = .TRUE.
         rulesIniBoundary%bTypes(1:4) = dirichlet_boundary_type
       else
         call set_error('Unknown lateral_boundary_type: "' // trim(chTmp) // &
             & '", can be ZERO, PERIODIC, or DIRICHLET' , sub_name)
         return
       endif

       chTmp = fu_str_u_case(fu_content(nlIni,'top_boundary_type'))
       if (chTmp == 'ZERO') then
         rulesIniBoundary%bTypes(top_boundary) = zero_boundary_type
       elseif(trim(chTmp) == 'DIRICHLET')then
         rulesIniBoundary%ifBoundary = .TRUE.
         rulesIniBoundary%bTypes(top_boundary) = dirichlet_boundary_type
       else
         call set_error('Unknown top_boundary_type: "' // trim(chTmp) // &
             & '", can be ZERO or DIRICHLET' , sub_name)
         return
       endif

       rulesIniBoundary%bTypes(bottom_boundary) = surface_boundary_type
      
    endif  ! Definition of the boundary types
    !
    ! If the boundaries around the whole domain defined, look for headers
    !
    if(rulesIniBoundary%ifBoundary)then

      ! Timestep in boundary files
      rulesIniBoundary%bInterval = fu_set_named_interval(fu_content(nlIni,'boundary_time_step'))
      rulesIniBoundary%start_time =  fu_round_closest(start_time, &
                                                & rulesIniBoundary%bInterval, backwards)  ! Forward run assumed
      !
      ! Having, may be, something to constrain at the boundary, let's get the input files.
      ! Two types of the boundaries are possible:
      ! - from maps stored as GrADS or GRIB or NETCDF files
      ! At this level, they are all treated same way: 
      ! boundary_data_file = GRIB/GRADS/NETCDF/POINT_DATA <fileName>
      !
      call get_items(nlIni, 'boundary_header_filename', ptrItems, iTmp)
      if(iTmp >0 )then
        allocate(rulesIniBoundary%bFiles(iTmp),stat=iStat)
        if(fu_fails(iStat == 0, 'Failed to allocate boundary files structure', sub_name))return !
        ! Set the file format and decode the template
        !
        do iStat = 1, iTmp
          rulesIniBoundary%bFiles(iStat)%bHeaderFNm = fu_process_filepath(fu_content(ptrItems(iStat)), &
                                                          &  must_exist=.true.)
        end do
      endif   ! header file exists
    endif  ! if boundaries
    !
    ! The reprojection aliasing may be smoothed by randomization
    !
    rulesIniBoundary%ifRandomise = .not. fu_str_u_case(fu_content(nlStandardSetup, &
                                                                & 'randomise_reprojection')) == 'NO'
    
    rulesIniBoundary%defined = silja_true

  end subroutine set_ini_boundary_rules_from_nl

  
  !************************************************************************************
  
  subroutine set_missing_b_rules(b_rules)
    implicit none
    type(Tini_boundary_rules), intent(out) :: b_rules
    nullify(b_rules%chIniFile, b_rules%qInitial)
    b_rules%ifBoundary = .false.
    b_rules%defined = silja_false
  end subroutine set_missing_b_rules


  !************************************************************************************

  logical function if_ini_boundary_rules_defined(rulesIniBoundary)
    implicit none
    type(Tini_boundary_rules), intent(in) :: rulesIniBoundary
    if_ini_boundary_rules_defined = rulesIniBoundary%defined == silja_true
  end function if_ini_boundary_rules_defined


  !***********************************************************************************

  subroutine add_ini_boundary_input_needs(rulesIniBoundary, &
                                        & q_met_dyn, q_met_stat, q_disp_dyn, q_disp_stat)
    !
    ! Adds the requested quantities to the provided arrays
    !
    implicit none

    ! Imported parameters
    type(Tini_boundary_rules), intent(in) :: rulesIniBoundary
    integer, dimension(:), intent(inout) :: q_met_dyn, q_met_stat, q_disp_dyn, q_disp_stat

    ! Local variables
    integer :: iTmp

   ! cell mass might be needed to initialize ones or time tracer
   ! Since advection needs it anyhow, just request 
   iTmp = fu_merge_integer_to_array(disp_cell_airmass_flag, q_disp_dyn)

    if(rulesIniBoundary%ifBoundary)then
      if(any(rulesIniBoundary%bTypes(1:size(rulesIniBoundary%bTypes)) == dirichlet_boundary_type))then
        iTmp = fu_merge_integer_to_array(temperature_flag, q_met_dyn)
        iTmp = fu_merge_integer_to_array(pressure_flag, q_met_dyn)
      endif
    endif

  end subroutine add_ini_boundary_input_needs


  !***********************************************************************************
  !***********************************************************************************
  !
  ! INITIAL CONDITIONS
  !
  !***********************************************************************************
  !***********************************************************************************

  subroutine set_initial_conditions(mapMass, mapPx, mapPy, mapPz, mapDD, mapWD, &
                                  & pNorthBoundary, pSouthBoundary, &
                                  & disp_buf, &
                                  & rulesIniBoundary, now, direction, &
                                  & meteoMarketPtr, dispersionMarketPtr)
    !
    ! Actually sets the initial conditions via setting the values of maps of TMassMap object
    ! and, if needed, dispersion minimarket variables. 
    ! Must be called after all these structures are properly initialized
    !
    ! For model run, initializes only concentration and advection moment. For
    ! dispersion minimarket - any field that can be put there.
    !
    implicit none

    ! Imported parameters
    type(Tini_boundary_rules), intent(in) :: rulesIniBoundary
    type(Tmass_map), intent(inout), target  :: mapMass, mapPx, mapPy, mapPz, mapDD, mapWD
    type(TboundaryStruct), intent(inout) :: pNorthBoundary, pSouthBoundary
    type(Tfield_buffer), pointer :: disp_buf
    type(silja_time), intent(in) :: now
    integer, intent(in) :: direction
    type(mini_market_of_stacks), pointer :: meteoMarketPtr, dispersionMarketPtr

    ! Local variables
    type(silja_shopping_list) :: shopLst
    integer :: iFile, iQ, iQDispStackST, iQDispStackMT, iQmassMap, nQmarketST, nQmarketMT, nFields_updated, quantity
    type(Tsilam_namelist), pointer :: nlPtr
    type(silja_time) :: time_tmp
    integer, dimension(:), pointer :: arQDispStackMT, arQDispStackST, arQmassMap, & 
                                     & arQmarketST, arQmarketMT, arIntWork
    type(Tmass_map_ptr), dimension(4) :: ptrMap
    logical :: ifNewDataHere
    character(len=*), parameter :: sub_name="set_initial_conditions" 
    integer, dimension(7), parameter :: MassMapQ = (/concentration_flag, drydep_flag, wetdep_flag, &
      & mass_in_air_flag,  advection_moment_X_flag, advection_moment_Y_flag, advection_moment_Z_flag/)

    !
    ! Stupidity check
    !
    if(associated(rulesIniBoundary%qInitial))then
      if(size(rulesIniBoundary%qInitial) < 1 .or. size(rulesIniBoundary%qInitial) > 1000)then
        call msg('No initialization is requested')
        return
      endif
    else
      call msg('No initialization is requested')
      return
    endif

    if(.not. (defined(mapMass) .and. associated(disp_buf)))then
      call set_error('mapMass and/or dispersion buffer are not associated','set_initial_conditions')
      return
    endif

    arIntWork => fu_work_int_array(5*max_quantities)
    if(error)return
    arQDispStackMT => arIntWork(1:max_quantities)
    arQDispStackST => arIntWork(1*max_quantities+1:2*max_quantities)
    arQmassMap =>     arIntWork(2*max_quantities+1:3*max_quantities)
    arQmarketST =>    arIntWork(3*max_quantities+1:4*max_quantities)
    arQmarketMT =>    arIntWork(4*max_quantities+1:5*max_quantities)

    !
    ! Run over the initialization quantities picking up only those belonging to dispersion_stack
    !
    iQmassMap = 0
    iQDispStackMT = 0
    iQDispStackST = 0


    !!! Figure out which quantities go where
    !!! Everything must be already created in the market
    call supermarket_quantities(dispersionMarketPtr, met_src_missing, &
                                  & single_time_stack_flag,  arQmarketST, nQmarketST )
    call supermarket_quantities(dispersionMarketPtr, met_src_missing, &
                                  & multi_time_stack_flag,  arQmarketMT, nQmarketMT )


    do iQ = 1, size(rulesIniBoundary%qInitial)
      quantity = rulesIniBoundary%qInitial(iQ)
      if (any( quantity == MassMapQ(:) )) then

          iQmassMap = iQmassMap + 1
          arQmassMap(iQmassMap) = quantity
          call msg('Ordering mass-map initialization of:' + fu_quantity_string(quantity))

      elseif (any( quantity == arQmarketST(1:nQmarketST) )) then 

        !! present in dispersion_singletime
        iQDispStackST = iQDispStackST + 1
        arQDispStackST(iQDispStackST) = quantity
        call msg('Ordering dispersion-stack initialization of singetime:' // &
               & fu_quantity_string(quantity), quantity)

      elseif (any( quantity == arQmarketMT(1:nQmarketMT)))then
        !! present in dispersion_singletime
        iQDispStackMT = iQDispStackMT + 1
        arQDispStackMT(iQDispStackMT) = quantity
        call msg('Ordering dispersion-stack initialization of multitime (past):' // &
               & fu_quantity_string(quantity), quantity)

      else
            call msg_warning('Failed initialization of: ' // &
                   & fu_quantity_string(quantity), sub_name)
            call set_error('No initialisation for quantities missing from disp_buffer', sub_name)
            return
      endif
    end do
    arQDispStackST(iQDispStackST+1) = int_missing
    arQDispStackMT(iQDispStackMT+1) = int_missing
    arQmassMap(iQmassMap+1) = int_missing

    if(iQDispStackMT+ iQDispStackST + iQmassMap > 0)then
      !
      ! Scan the input files and consume the values
      ! Have to check two places: dispersion buffer and mass-map structures
      !
      nlPtr => fu_create_namelist()
      if(error)return

      do iFile = 1, size(rulesIniBoundary%chIniFile)
        call add_namelist_item(nlPtr,'initial_file',rulesIniBoundary%chIniFile(iFile))
        if(error)return
      enddo


      if(iQDispStackMT > 0)then
        call msg("Updating Dispersion stack MT:", arQDispStackMT(1:iQDispStackMT))
        time_tmp = disp_buf%time_past 
        shopLst = fu_set_shopping_list(met_src_missing, &
                                     & arQDispStackMT, &
                                     & time_tmp, & ! First time boundary
                                     & time_tmp, &
                                     & level_missing, level_missing)
        if(error)return

        ! overwrite_field_anytime also handles wrong levels coming from the input
        ! 
        call fill_minimarket_from_namelist(dispersionMarketPtr, nlPtr, 'initial_file', &
                                         & shopLst, time_tmp, dynamic_map, overwrite_field_anytime, &
                                         & wdr_missing,&
                                         & dispersion_gridPtr, &
                                         & 5, .true., & ! iAccuracy, ifAdjustGrid
                                         & ifNewDataHere)
!                                         & shopLst, time_missing, single_time_stack_flag)
        if(error)return
        if (fu_fails(ifNewDataHere, "Faild to update Dispersion stack MT", sub_name)) return

      endif ! if some quantities for intialization are in dispersion_stack
      if(iQDispStackST > 0)then
        call msg("Updating Dispersion stack ST:", arQDispStackST(1:iQDispStackST))
        
        shopLst = fu_set_shopping_list(met_src_missing, &
                                     & arQDispStackST, &
                                     & now, & ! First time boundary
                                     & now, &
                                     & level_missing, level_missing)
        if(error)return

        call fill_minimarket_from_namelist(dispersionMarketPtr, nlPtr, 'initial_file', &
                                         & shopLst, now, instant_map, overwrite_field_anytime, &
                                         & wdr_missing,&
                                         & dispersion_gridPtr, &
                                         & 5, .true., & ! iAccuracy, ifAdjustGrid
                                         & ifNewDataHere)
!                                         & shopLst, time_missing, single_time_stack_flag)
        if(error)return
        if (fu_fails(ifNewDataHere, "Faild to update Dispersion stack ST", sub_name)) return

      endif ! if some quantities for intialization are in dispersion_stack


      !----------------------------------------------------------------------
      !
      ! Step 2: initialize the massMap objects - concentration, deposition, and moment fields
      !
      if(iQmassMap > 0)then
        !
        ! Something is found to be initialized in the map of cocktails. Get it!
        !
        call msg('Mass map')
        if(iQmassMap > 4)then
          call msg("iQmassMap", iQmassMap)
          do iQ = 1, iQmassMap
            call msg(fu_quantity_string(arQmassMap(iQ)))
          enddo
          call set_error('Only 4 cocktail maps are known','set_initial_conditions')
          return
        endif
        do iQ = 1, iQmassMap
          select case(arQmassMap(iQ))
            case(concentration_flag, mass_in_air_flag)
              ptrMap(iQ)%ptrMassMap => mapMass

            case(drydep_flag)
              ptrMap(iQ)%ptrMassMap => mapDD

            case(wetdep_flag)
              ptrMap(iQ)%ptrMassMap => mapWD

            case(advection_moment_X_flag)
              ptrMap(iQ)%ptrMassMap => mapPx

            case(advection_moment_Y_flag)
              ptrMap(iQ)%ptrMassMap => mapPy

            case(advection_moment_Z_flag)
              ptrMap(iQ)%ptrMassMap => mapPz

            case default
              call set_error('Do not support massMap intialization for quantity:' + &
                           & fu_quantity_string(arQmassMap(iQ)),'set_initial_conditions')
              return
          end select
        end do        ! quantities to be initialised

        do iQ = iQmassMap+1, 4
          nullify(ptrMap(iQ)%ptrMassMap)
        end do

        call update_mass_map_from_namelist(nlPtr, 'initial_file', ptrMap, iQmassMap, now, &
                                         & rulesIniBoundary%ifRandomise, nFields_updated)
        if(error)return
        !
        ! Check whether this mass map is global and thus poles must be initialised
        !
        if(.not. (defined(pNorthBoundary) .and. defined(pSouthBoundary) ))then
                call set_error("pNorthBoundary and pSouthBoundary nust be set by now!!!",&
                & "set_initial_conditions")
                return
        endif
        do iQ = 1, iQmassMap
          if( any(arQmassMap(iQ) == (/concentration_flag, mass_in_air_flag/))) then
              if (pNorthBoundary%boundaryType == polar_boundary_type) &
                  & call check_poles(ptrMap(iQ)%ptrMassMap, pNorthBoundary, ptrMap(iQ)%ptrMassMap%ny)
              if (pSouthBoundary%boundaryType == polar_boundary_type) &
                  & call check_poles(ptrMap(iQ)%ptrMassMap, pSouthBoundary, 1)
          endif
        end do


        if(fu_fails(nFields_updated > 0, fu_str(iQmassMap) + &
                                      & '- massMap quantities requested but no fields updated', &
                                      & 'set_initial_conditions'))return


      endif ! if some quantities for intialization are in massMap objects

      call destroy_namelist(nlPtr)

      
    endif ! if anything to initialize at all

    call free_work_array(arIntWork)

    CONTAINS

    !=============================================================
    
    subroutine check_poles(pMap, pBoundary, indYEdge)
      !
      ! Poles are not initialised by above mapping, we have to do it manually here
      ! It is simple: concentration over pole must be average over the near surrounding
      !
      implicit none
      
      ! Imported parameters
      type(TMass_map), pointer :: pMap
      type(TboundaryStruct), intent(inout) :: pBoundary
      integer, intent(in) :: indYEdge
      
      ! Local variables
      integer :: ix, iy, iz, iSrc, iSpecies, iPole
      real :: fArea

      !
      ! Get the mean concentration. 
      !
      pBoundary%polarCapMassMom = 0.  !(1:pMap%nSpecies, 1:pMap%nSrc, 1:2, 1:pMap%n3d)
      do iz = 1, pMap%n3d
        fArea = 0.
        do ix = 1, pMap%nx
          fArea = fArea + fu_cell_size(pMap%gridTemplate, ix, indYEdge)
          do iSrc = 1, pMap%nSrc
            do iSpecies = 1, pMap%nSpecies
              pBoundary%polarCapMassMom(iSpecies,iSrc,1,iz) = &
                                          & pBoundary%polarCapMassMom(iSpecies,iSrc,1,iz) + &
                                          & pMap%arM(iSpecies,iSrc,iz,ix,indYEdge)
            end do
          end do  ! isrc
        end do  ! ic
        !
        ! Store volume-mean concentration to the cap
        !
        pBoundary%polarCapMassMom(1:pMap%nSpecies,1:pMap%nSrc,1,iz) = &
                    & pBoundary%polarCapMassMom(1:pMap%nSpecies,1:pMap%nSrc,1,iz) * &
                    & pBoundary%polarCapArea / fArea
        pBoundary%polarCapMassMom(1:pMap%nSpecies,1:pMap%nSrc,2,iz) = 0.    ! moment to zero: cell centre
      end do ! iz
        
    end subroutine check_poles

  end subroutine set_initial_conditions


  !***********************************************************************************

  subroutine store_restart_file()
    !
    ! Prepares to the next run
    !
    implicit none

    call set_error('Does not work yet','store_restart_file')

  end subroutine store_restart_file


  !***********************************************************************************
  !***********************************************************************************
  !
  ! BOUNDARY CONDITIONS
  !
  !***********************************************************************************
  !***********************************************************************************
  !
  ! The idea of the boundary conditions for concentrations is that we have an Eulerian grid
  ! with +1 grid cell around it. This +1 element may serve for the boundary conditions
  ! Basically, we just have to understand how to set this band fast.
  ! The rest is simple:
  ! before the main advection simply inject the appropriate amount (with proper plume size
  ! along the corresponding direction - Courant nbr can be >1)
  !

  !***********************************************************************************

  subroutine update_btypes(grid, r)
    !
    ! Update boundary types accounting for neighbours and grid coverage
    ! Set btypes for the whole domain, and then reset them if a  neighbour exists
    ! 
    TYPE(silja_grid) :: grid !! Should be whole-domain grid
    type(Tini_boundary_rules), intent(inout) :: r

    logical :: ifPeriodicX, ifPeriodicY
    integer :: iStat

    
    ifPeriodicX =  ( r%bTypes(eastern_boundary) ==  periodic_boundary_type)!! From ini file
    if ( fu_ifLonGlobal(grid)) then
       r%bTypes(eastern_boundary:western_boundary) = periodic_boundary_type
       ifPeriodicX = .true.
    endif

    ifPeriodicY =  (r%bTypes(northern_boundary)   ==  periodic_boundary_type)! from ini file

    if (.not. ifPeriodicY) then ! Check if poles needed
      !Note that ALL norterhnmost domains will get polar_boundary_type
      if (fu_ifPolarCapPossible(grid, northern_boundary)) then
          r%bTypes(northern_boundary) = polar_boundary_type
      endif

       if (fu_ifPolarCapPossible(grid, southern_boundary)) then
         r%bTypes(southern_boundary) = polar_boundary_type
       end if
    endif


    call smpi_reset_periodic_topology(ifPeriodicX, ifPeriodicY)

    do iStat = 1, 4
      if(adv_mpi_neighbours(iStat) /= int_missing)then
        r%bTypes(iStat) = smpi_comm_boundary_type
      endif
    end do

    call msg("r%bTypes after update_r%bTypes:", r%bTypes(1:6))

  end subroutine update_btypes

  !***********************************************************************************


  subroutine init_boundary_structures(rulesIniBoundary, &
                                    & wdr, &
                                    & speciesTransp, nSpeciesTransp, &
                                    & mapDisp, &
                                    & boundaryMarketPtr, &
                                    & boundStructArray, &
                                    & pBBuffer)
    ! 
    ! Should analyse the input files and mapping and create boundary structures and boundary stack
    ! 
    implicit none

    ! Imported parameters
    type(Tini_boundary_rules), intent(inout), target :: rulesIniBoundary
    type(silam_species), dimension(:), intent(in) :: speciesTransp
    integer, intent(in) :: nSpeciesTransp
    type(Tmass_Map), intent(in) :: mapDisp
    type(TboundaryStructPtr), dimension(:), pointer :: boundStructArray
    type(mini_market_of_stacks), pointer ::  boundaryMarketPtr
    type(TboundaryBuffer), pointer :: pBBuffer
    type(silja_wdr), pointer :: wdr


    ! Local variables
    character(len=fnlen) :: bFileNm, bFileTp
    type(Tinput_content) :: inputContent
    integer :: iTmp, iStat
    type(Tsilam_namelist_group), pointer :: nlGrpPtr
    type(Tsilam_namelist), pointer :: nlPtr
    type(Tsilam_nl_item_ptr), dimension(:), pointer :: ptrItems
    type(silam_sp) :: spContent
    character(len = fnlen) :: chtmp
    character(len=20) :: ch_bModeVal, ch_tModeVal
    logical :: iffound
    type(silja_shopping_list) :: shopLst
    type(silam_sp), dimension(:), pointer :: fnames
    type(silja_field_id), dimension(:), pointer :: lstid
    type(silja_field_id) ::idTmp
    TYPE(silja_grid), dimension(6) :: grid
    REAL :: corner_x, corner_y, pole_x, pole_y, dx_deg, dy_deg
    INTEGER :: number_of_points_x, number_of_points_y
    LOGICAL :: if_corner_in_geo_coord, if_south_pole
    type(wdr_ptr), dimension(:), pointer :: wdrar
    type(silam_vertical), dimension(6) :: vertical
    type(grads_template), dimension(max_met_files) ::  fnm_multitime_tmpl, fnm_singletime_tmpl
    type(silam_fformat), dimension(max_met_files) :: fFmt_multitime, fFmt_singletime

    integer :: nSingletimeFs, nMultitimeFs, ifile, nlinks, ilink, iTSpecies, iId, &
             & iBSpecies, maxlevs, iHeader, indexvar, nfiles, nSelected
    integer, dimension(:), pointer :: indices, iTransportSpeciesInvolved
    real :: bModeVal, tModeVal, b2tFactor
    character(len=substnmlen) :: bSubstNm, tSubstNm
    type(silam_species) :: species
!    type(Tini_boundary_rules), pointer :: bcrPtr
    type(silja_shopping_list), pointer :: shpLstPtr
    type(Taerosol_mode) :: t_aerosol_mode, b_aerosol_mode
    logical :: included

    character(len=*), parameter :: sub_name="init_boundary_structures" 

!    bcrPtr => rulesIniBoundary

    iTransportSpeciesInvolved => fu_work_int_array()
    indices => fu_work_int_array()
    
    inputContent = input_content_missing
    iTransportSpeciesInvolved(1:nSpeciesTransp) = 0
    !
    ! Stupidity checking
    !
    call msg('Initializing boundary conditions')
    nullify(ptrItems, fnames)



    if(any(rulesIniBoundary%bTypes(1:6) == dirichlet_boundary_type))then

      if(fu_fails(allocated(rulesIniBoundary%bFiles),'Dirichlet boundary requested but no boundary files', &
                                                       & 'init_boundary_structures'))return
      maxlevs = 0
      !-------------------------------------------------------------------------
      !
      ! Go through the header files setting the parameters of the input files in ini_boundary_rules 
      ! and analysing their content
      !
      do iHeader = 1, size(rulesIniBoundary%bFiles)
      
        call set_missing(rulesIniBoundary%bFiles(iHeader)%shopLst) 
        call set_missing(rulesIniBoundary%bFiles(iHeader)%vertical, .true.)
        rulesIniBoundary%bFiles(iHeader)%unit_binary(:) = int_missing
        rulesIniBoundary%bFiles(iHeader)%nTSpecies = 0
        !
        ! Read the header file
        !
        iFile = fu_next_free_unit()
        open(file=rulesIniBoundary%bFiles(iHeader)%bHeaderFNm, &
        & unit=iFile, action='read', status='old', iostat=iStat)
        if(iStat /= 0)then
          call set_error('Failed to open boundary map file:' + &
                       & rulesIniBoundary%bFiles(iHeader)%bHeaderFNm,'init_boundary_structures')
          return
        endif
        nlPtr => fu_read_namelist(iFile, .true.)
        close(iFile)

        if(fu_fails(associated(nlPtr),'Boundary species mapping namelist not associated',&
                                    & 'init_boundary_structures'))return

        rulesIniBoundary%bFiles(iHeader)%bNames = fu_content(nlPtr, 'boundary_names')
        rulesIniBoundary%bFiles(iHeader)%ifToBoundary(:) = .false.    
        do iTmp = 1, 6
          included = index(rulesIniBoundary%bFiles(iHeader)%bNames, trim(rulesIniBoundary%bNames(iTmp)))>0
          if(included .and. rulesIniBoundary%btypes(itmp) == dirichlet_boundary_type) then
            rulesIniBoundary%bFiles(iHeader)%ifToBoundary(iTmp) = .true.  
          endif  
        enddo

        if(fu_content(nlPtr, 'ifClimatology') == 'NO')then
          rulesIniBoundary%bFiles(iHeader)%ifStatic = .false.
          rulesIniBoundary%bFiles(iHeader)%fieldType = dynamic_map
        else
          if(index(fu_content(nlPtr, 'climatologyTimestep'),'MONTHLY') == 1)then
            rulesIniBoundary%bFiles(iHeader)%fieldType = monthly_climatology
            rulesIniBoundary%bFiles(iHeader)%ifStatic = .false.
          elseif(index(fu_content(nlPtr, 'climatologyTimestep'),'STATIC') == 1)then
            rulesIniBoundary%bFiles(iHeader)%fieldType = static_climatology
            rulesIniBoundary%bFiles(iHeader)%ifStatic = .true.
          else
            call set_error('Only monthly and static values accepted for climatological boundaries ', &
                         & 'init_boundary_structures')
            return
          endif
        endif

        rulesIniBoundary%bFiles(iHeader)%bFileFormat = &
                                        & fu_input_file_format(fu_content(nlPtr, 'file_format'))

       ! Get the content of input data-files 

        call get_items(nlPtr, 'boundary_file', ptrItems, nFiles)

        if(nfiles<1)then
          call set_error('No boundary data files given in header file','init_boundary_structures')
          return
        endif

        allocate(rulesIniBoundary%bFiles(iHeader)%bFnameTemplate(nFiles), stat = iStat)
        if(fu_fails(iSTat ==0,'Failed to allocate boundary file names array', &
                            & 'init_boundary_structures'))return

        do iFile = 1, nFiles

          chtmp = fu_content(ptrItems(iFile))
          ! Allow for hat in boundary ini file
          chtmp  =  fu_extend_grads_hat( chtmp, &
                                    & rulesIniBoundary%bFiles(iHeader)%bHeaderFNm )

          call decode_template_string(chtmp, &
                                    & rulesIniBoundary%bFiles(iHeader)%bFnameTemplate(iFile))

        enddo

        !
        ! First get the full list of fields available in the input files
        shopLst = fu_set_shopping_list(met_src_missing, &
                                   & (/accept_all_quantities/), &
                                   & time_missing, & 
                                   & time_missing, &
                                   & level_missing, level_missing)
        if(error)return

        do iFile = 1, nFiles
             
          inputContent = input_content_missing !! Yes, it is a memory leak, but
                        !          content should be different for different
                        !          files otherwise we are in trouble
                        ! Could be fine if content elements would be allocatable

          if(rulesIniBoundary%bFiles(iHeader)%bFnameTemplate(iFile) == template_missing) exit
          if(fu_if_input_file(rulesIniBoundary%bFiles(iHeader)%bFileFormat))then
            call FNm_from_single_template(rulesIniBoundary%bFiles(iHeader)%bFnameTemplate(iFile), &
                                        & rulesIniBoundary%start_time, &
                                        & fnames, &
                                        & ifAdd = .false., & 
                                        & ifStrict = .false., &
                                        & ifAllowZeroFcLen = .true., &
                                        & ifWait = .false., &
                                        & fc_step = fu_obstime_interval(wdr))
            if(error)return
       
          else
            !
            ! Template stores a request of some non-file input in its collection
            !
            call enlarge_array(fnames,1)
            if(error)return
            fnames(1)%sp = fu_collection(rulesIniBoundary%bFiles(iHeader)%bFnameTemplate(iFile))
            if(size(fnames) > 1) fnames(2)%sp = ''
          endif  ! If this is a file
          do iTmp = 1, size(fnames)
            if(fnames(iTmp)%sp=='')exit
            call read_input_content(fnames(iTmp)%sp, wdr_missing, &
                                  & rulesIniBoundary%bFiles(iHeader)%bFileFormat, &
                                  & shopLst, inputContent)
          enddo
        enddo
        !
        ! Make a list of field IDs available from input content
        !
        call id_lst(inputContent, lstId)
        
#ifdef DEBUG
        call msg("Boundary Input content: "+trim(fnames(1)%sp))
        do iId = 1, fu_nbr_of_fields(inputContent)
            call msg(fu_quantity_string(fu_quantity(lstId(iId)))//fu_str(fu_species(lstId(iId))))
        enddo
#endif


        !
        ! Set the mapping between species in input files and boundary
        ! stack. Note that the species in boundary stack are identical
        ! to transport species, quantity and vertical still stay as in
        ! input files.
        ! 
        call get_items(nlPtr, 'par_str', ptrItems, nLinks)
        if(nLinks < 1)then
          call set_error('No mapping given between boundary and transport species', &
                       & 'init_boundary_structures')
          return
        endif

        do iLink = 1, nLinks
          chtmp = fu_content(ptrItems(iLink))
          !
          ! Modes can be named as 'GAS' or via their mean diameters 
          !
!call msg('Resolving link:' + chtmp)
          read(unit=chtmp, fmt=*, iostat=iStat) bSubstNm, tSubstNm, &
                                              & ch_bModeVal, ch_tModeVal, b2tFactor
          if(fu_fails(iStat == 0,'Failed to parse string:'+chtmp,'init_boundary_structures'))return
          !
          ! Boundary mode...
          !
          if(fu_str_u_case(ch_bModeVal) == 'GAS')then      ! Boundary-file mode...
            b_aerosol_mode = in_gas_phase
          else
            read(unit=ch_bModeVal, fmt=*, iostat=iStat) bModeVal
            if (fu_fails(iStat == 0, 'Failed to parse:' // ch_bModeVal, sub_name)) return
            if (fu_fails(bModeVal > 0, 'bmode must be >0 or "GAS", in file: ' // ch_bModeVal, sub_name)) return
            b_aerosol_mode = fu_set_mode(fixed_diameter_flag, bModeVal, bModeVal, bModeVal)
          endif
          !
          ! ...and transported mode
          !
          if(fu_str_u_case(ch_tModeVal) == 'GAS')then     ! Transport-species mode
            t_aerosol_mode = in_gas_phase
            tModeVal = 0
          else
            read(unit=ch_tModeVal, fmt=*, iostat=iStat) tModeVal
            if (fu_fails(iStat == 0, 'Failed to parse:' // ch_tModeVal, sub_name)) return
            if (fu_fails(tModeVal > 0, 'tmode must be >0 or "GAS", in file: ' // ch_tModeVal, sub_name)) return
            t_aerosol_mode = fu_set_mode(fixed_diameter_flag, tModeVal, tModeVal, tModeVal)
          endif
          !
          ! Now find the most-appropriate transport species
          !
          call select_species(speciesTransp, nSpeciesTransp, &
                            & tSubstNm, &
                            & t_aerosol_mode, real_missing, &
                            & indices, nSelected)

          iTSpecies = int_missing

          if (nSelected > 0) iTSpecies = indices(1) !The species selected
          if (nSelected > 1) then
             call msg("More than one transport species matches the boundary")
             call msg("Transp subst: "//trim(tSubstNm)//", mode: " // trim(ch_tModeVal))
             call msg("Matched transport species:")
             do itmp = 1, nSelected
                call report(speciesTransp(indices(itmp)))
             enddo
             call set_error("Non-unique match",sub_name)
             return
          elseif (nSelected == 0 .or. iTSpecies < 1)then
            call msg("nSelected itSpecies", nSelected, iTSpecies)

            call msg('Boundary header links to non-existent transport species: ' + &
                   & tSubstNm + '; modesize:', tModeVal)
            call report(t_aerosol_mode)
            call msg("speciesTransp")
            do iTmp=1,nSpeciesTransp
              call msg("speciesTransp", iTmp)
              call report(speciesTransp(iTmp))
              call msg("")
            enddo
            call set_error('Boundary header links to non-existent transport species', &
                         & 'init_boundary_structures')
            return
          endif
          iTransportSpeciesInvolved(iTSpecies) = 1
          !
          ! Transport species of the link is found. Now get the boundary species linked to it
          !
          iffound = .false. 
!call msg('mapping mode')                         ! from mapping file
!call report(b_aerosol_mode)
          do iId = 1, fu_nbr_of_fields(inputContent)
            if (.not.(fu_quantity(lstId(iId)) == concentration_flag .or. &
                    & fu_ifDiagnosticQuantity(fu_quantity(lstId(iId)), concentration_flag))) cycle

            if (.not.(bSubstNm == fu_substance_name(lstId(iId)))) cycle
!call msg('boundary file mode:')  ! from boundary file
!call report(fu_mode(lstId(iId)))
            if (fu_modes_match(b_aerosol_mode, fu_mode(lstId(iId)))) then
              ifFound = .true.
              exit
            endif
          enddo
          if(.not. iffound)then
            call msg('Boundary header refers to non-existent input-file species: '//trim(bSubstNm) // &
                   & ', modesize: '//trim(ch_bModeVal))
            do iId = 1, fu_nbr_of_fields(inputContent)
              if (.not.(bSubstNm == fu_substance_name(lstId(iId)))) cycle
              call msg('substance:' + fu_substance_name(lstId(iId)))
              call report(fu_mode(lstId(iId)))
            end do
            call set_error('Boundary header refers to non-existent input-file species', &
                         & 'init_boundary_structures')
            return
          endif

          ! If necessary, add the id to boundary shopping list
          !
          idTmp = lstId(iId)
         
          call set_level(idTmp, level_missing)
          !call set_valid_time(idTmp, time_missing)
          if(error)call unset_error('init_boundary_structures')

          call add_shopping_variable(rulesIniBoundary%bFiles(iHeader)%shoplst, idTmp, 2)
          if(.not. fu_field_id_in_list(idTmp, &
                                     & rulesIniBoundary%bFiles(iHeader)%shoplst, &
                                     & indexVar))then
            call set_error('Problems with shopping list','init_boundary_structures')
          endif
          if(error)return

          ! Add the link to mapping
          !
          call set_species(species, fu_get_material_ptr(tSubstNm), t_aerosol_mode) !&
!                         & fu_set_mode(fixed_diameter_flag, tModeVal, tModeVal, tModeVal))
          call set_species(idTmp, species)
          
          shpLstPtr => rulesIniBoundary%bFiles(iHeader)%shoplst
          call add_maplink_to_shopVar(fu_shopping_var(shpLstPtr, &
                                                    & indexVar), idTmp, b2tFactor)

          ifFound = .false.
          do iTmp = 1, rulesIniBoundary%bFiles(iHeader)%nTspecies
            if(rulesIniBoundary%bFiles(iHeader)%tSpecies(iTmp) == iTSpecies)then
              iffound = .true.
              exit
            endif
          enddo
          if(.not. iffound)then
            rulesIniBoundary%bFiles(iHeader)%nTspecies = &
                         & rulesIniBoundary%bFiles(iHeader)%nTSpecies + 1
            rulesIniBoundary%bFiles(iHeader)&
                         & %tSpecies(rulesIniBoundary%bFiles(iHeader)%nTSpecies) = iTSpecies
            rulesIniBoundary%bFiles(iHeader)&
                         & %fieldInStack(rulesIniBoundary%bFiles(iHeader)%nTSpecies) = idTmp  
          endif
        enddo ! cycle over links in mapping


        ! Also have to set the vertical
        ! All the species given by one header file are supposed to have the same vertical, so
        ! cycle over the fields in input_content, adding the level of each field corresponding 
        ! to first shopping variable
        !
        shpLstPtr => rulesIniBoundary%bFiles(iHeader)%shoplst

        do iId = 1, fu_nbr_of_fields(inputContent)
          if(.not. fu_fld_corresponds_to_shop_var(lstId(iId), &
                    & fu_shopping_var(shpLstPtr, 1)))cycle      
          if(.not. defined(rulesIniBoundary%bFiles(iHeader)%vertical))then
            call set_vertical(fu_level(lstId(iId)), rulesIniBoundary%bFiles(iHeader)%vertical)
          elseif(.not. fu_level_belongs_to_vertical(fu_level(lstId(iId)), &
                                                  & rulesIniBoundary%bFiles(iHeader)%vertical))then
            call add_level(rulesIniBoundary%bFiles(iHeader)%vertical, fu_level(lstId(iId)))
          endif
        enddo
        call arrange_levels_in_vertical(rulesIniBoundary%bFiles(iHeader)%vertical)      
        maxlevs = max(fu_NbrOfLevels(rulesIniBoundary%bFiles(iHeader)%vertical),maxlevs)
        call msg('Vertical in boundary file')
        call report(rulesIniBoundary%bFiles(iHeader)%vertical)
        call msg('')
        call msg('Grid in boundary file')
        call report(fu_grid(lstId(1)))
        call msg('')
      enddo !cycle over header files
      if(associated(lstId))deallocate(lstid)

      ! Initialize boundary stacks and structures
      ! First create grids and verticals for each boundary
      !
      do iTmp = 1,6
        grid(iTmp) = fu_boundary_grid(dispersion_grid, rulesIniBoundary%bNames(iTmp))
      enddo

      vertical(1:4) = dispersion_vertical
      call set_vertical(fu_upper_boundary_of_layer(fu_level(dispersion_vertical, nz_dispersion)), vertical(5))
      call set_vertical(fu_lower_boundary_of_layer(fu_level(dispersion_vertical, 1)), vertical(6))
      allocate(wdrar(6))

      ! Cycle over boundaries to create boundary structures
      !
      do iTmp = 1, 6
        ! input files for bc-dr
        nSingletimeFs = 0
        nMultitimeFs = 0
        do iFile = 1, size(rulesIniBoundary%bFiles)
          if(index(rulesIniBoundary%bFiles(iFile)%bNames, rulesIniBoundary%bNames(iTmp))>0)then
            if(rulesIniBoundary%bFiles(iFile)%ifStatic)then
              do iBSpecies = 1, size(rulesIniBoundary%bFiles(iFile)%bFnameTemplate)
                nSingletimeFs = nSingletimeFs + 1
                fnm_singletime_tmpl(nSingletimeFs) = &
                                      & rulesIniBoundary%bFiles(iFile)%bFnameTemplate(iBSpecies)
                fFmt_singletime(nSingletimeFs) = rulesIniBoundary%bFiles(iFile)%bFileFormat
              enddo
            else   
              do iBSpecies = 1, size(rulesIniBoundary%bFiles(iFile)%bFnameTemplate) 
                nMultitimeFs = nMultitimeFs+1  
                fnm_multitime_tmpl(nMultitimeFs) = &
                                      & rulesIniBoundary%bFiles(iFile)%bFnameTemplate(iBSpecies)
                fFmt_multitime(nMultitimeFs) = rulesIniBoundary%bFiles(iFile)%bFileFormat
              enddo
            endif
          endif
        enddo

        boundStructArray(iTmp)%ptr => fu_makeBoundaryStruct(iTmp, &
                                                         & rulesIniBoundary%bNames(iTmp), &
                                                         & rulesIniBoundary%bTypes(iTmp), & 
                                                         & grid(iTmp),  &
                                                         & vertical(iTmp), &
                                                         & fnm_multitime_tmpl, &
                                                         & fnm_singletime_tmpl, &
                                                         & fFmt_multitime, &
                                                         & fFmt_singletime)
        if (error) return

        wdrar(iTmp)%ptr => boundStructArray(iTmp)%ptr%bc_wdr

      end do  ! over boundaries

      ! Initialize boundary minimarket
      call msg('Initiatlizing boundary market multi-time')
      ! Multi time stack initialization
      call initialize_mini_market(boundaryMarketPtr, &
                                & 'boundary_market', &
                                & 6, &      ! first dimension of the stack arrays
                                & 2, &      ! if>0, the second dimension in multiTime stack
                                & nSpeciesTransp * maxlevs+3, & ! for each timenode - fields
                                & 0, &      ! for each timenode - windfields
                                & nSpeciesTransp + 3, &  ! for each timenode - 3d fields
                                & 0, &      ! for each timenode  - 3d windfields
                                & .true., & ! replace_oldest_when_full
                                & wdrAr, &
                                & .true., &  !.false., & - ifSingleMetSrc
                                & .true.)   ! print_info_to_stdout
      ! Single-time stack initialization
      call msg('Initiatlizing boundary market single-time')
      call initialize_mini_market(boundaryMarketPtr, &
                                & 'boundary_market', &
                                & 6, &        ! first dimension of the stack arrays
                                & 0, &        ! if>0, the second dimension in multiTime stack
                                & nSpeciesTransp * maxlevs+3, & ! for each timenode - fields
                                & 0, &        ! for each timenode - windfields
                                & nSpeciesTransp + 3, &  ! for each timenode - 3d fields
                                & 0, &        ! for each timenode  - 3d windfields
                                & .true., &   ! replace_oldest_when_full
                                & wdrAr, &    ! WDR array - as many as stacks
                                & .true., &   ! .false. - ifSingleMetSrc
                                & .true.)     ! print_info_to_stdout
  
      IF (error) RETURN

      IF ((.NOT. minimarket_initialized(boundaryMarketPtr, multi_time_stack_flag)) .or. &
        & (.not. minimarket_initialized(boundaryMarketPtr, single_time_stack_flag))) THEN
        CALL set_error('Boundary market not init','init_boundary_structures')
        RETURN
      END IF
      call msg('Boundary market initialized')

      deallocate(wdrar)  

      ! Now prepare the vertical interpolation structure
      !
      do iHeader = 1, size(rulesIniBoundary%bFiles)
        do iTmp = 1, 6
          ! Interpolate the boundary data to ~dispersion vertical
          !
          rulesIniBoundary%bFiles(iHeader)%interpCoefBndInp2DispVert(iTmp)%ptr => &
                   & fu_vertical_interp_struct(rulesIniBoundary%bFiles(iHeader)%vertical, &
                                             & vertical(iTmp), &
                                             & grid(iTmp), &
                                             & linear, &! a bit crude but fast. Log-linear can be tried
                                             & one_hour*3., &  ! recommended update interval
                                             & 'boundary_to_dispersion') ! just name
          if(error)return
          
          ! Interpolate meteo data to boundary vertical (in order to refine the above):
          !
          rulesIniBoundary%bFiles(iHeader)%interpCoefMet2BndInp(iTmp)%ptr => &
                   & fu_vertical_interp_struct(meteo_vertical, &
                                             & rulesIniBoundary%bFiles(iHeader)%vertical, &
                                             & grid(iTmp), &
                                             & linear, & ! a bit crude but fast. Log-linear can be tried
                                             & one_hour*3., &  ! recommended update interval
                                             & 'meteo_to_boundary_' + fu_str(iTmp)) ! just name
          if(error)return

        enddo
      enddo

    else  ! no dirichlet boundaries but maybe polar
      !
      ! Cycle over boundaries to create non-Dirichlet boundary structures.
      !
      do iTmp = 1, 6
        boundStructArray(iTmp)%ptr => fu_makeBoundaryStruct(iTmp, rulesIniBoundary%bNames(iTmp), &
                                                          & rulesIniBoundary%bTypes(iTmp))
        if (error) return
      end do

    endif ! if any dirichlet boundaries

    !--------------------------------------------------------------------------------
    !
    ! Having the boundary structures created, make the boundary buffer
    !
    call init_boundary_buffer(rulesIniBoundary, boundStructArray, pBBuffer, mapDisp%nSrc, iTransportSpeciesInvolved)
    
    call free_work_array(iTransportSpeciesInvolved)
    call free_work_array(indices)

    
    CONTAINS


    !==================================================================================================

      function fu_makeBoundaryStruct(bIndex, name, boundary_type, &
                                   & bgrid, bvertical, &
                                   & fnm_multitime_tmpl, fnm_singletime_tmpl, &
                                   & fFmt_multitime, fFmt_singletime) result(struct)
      !
      ! Create a new boundary definition. 
      ! Note that the horizontal interpolation is supposed to be done at the moment of reading the 
      ! boundaries, so we need here only vertical interpolation, which is done in the working grid
      ! In most cases, verticalTo == dispersion_vertical, verticalFrom == boundary vertical 
      ! from external file
      !
      implicit none

      ! Imported parameters
      integer, intent(in) :: boundary_type, bIndex
      character, intent(in) :: name
      type(silam_vertical), intent(in), optional :: bvertical 
      type(silja_grid), intent(in), optional :: bgrid 
      type(grads_template), dimension(:), intent(in), optional ::  fnm_multitime_tmpl, fnm_singletime_tmpl
      type(silam_fformat), dimension(:), intent(in), optional :: fFmt_multitime, fFmt_singletime
      ! Return value
      type(TboundaryStruct), pointer :: struct 

      ! Local variables
      integer :: i, iS, iZ
      real :: rCap, arc
      type(silja_grid) :: gridTmp
      type(silja_field_id) :: id
      type(silja_field), pointer :: field
      real, dimension(:), pointer :: work_array

      type(TboundaryStruct), pointer :: structptr

      !structptr => struct

      select case(boundary_type)
        
        case(zero_boundary_type, smpi_comm_boundary_type, periodic_boundary_type, surface_boundary_type)
          !
          ! Zero and mpi boundaries do not require any structures here
          !  
          allocate(struct, stat=i)
          if(fu_fails(i == 0, 'Cannot allocate boundary structure','fu_makeBoundaryStruct'))return
          struct%boundaryType = zero_boundary_type
!          nullify(struct%polarCapConc)
          nullify(struct%Fld3dPtr)
          nullify(struct%meteo2BoundaryInterpHoriz)
          nullify(struct%meteo2BoundaryInterpVert)
          nullify(struct%ifBoundSpecies)
          call set_missing(struct%vertical, .true.)
          allocate(struct%bc_wdr)
          struct%bc_wdr = wdr_missing
      
          struct%defined = silja_true
          return

        case(polar_boundary_type)
          !
          ! Latitude-global closure boundary
          !
          allocate(struct, stat=i)
          if(fu_fails(i == 0, 'Cannot allocate boundary structure','fu_makeBoundaryStruct'))return

          struct%name = name
          struct%boundaryType = polar_boundary_type

          allocate(struct%polarCapMassMom(nSpeciesTransp, mapDisp%nSrc, 2, nz_dispersion), &
                 & struct%polarCapAirMass(nz_dispersion), &
                 & stat=i)
          if(fu_fails(i == 0, 'Cannot allocate polar cap structure','fu_makeBoundaryStruct'))return
          struct%nSrc = mapDisp%nSrc
          struct%nSp  = nSpeciesTransp
          struct%polarCapMassMom(1:nSpeciesTransp, 1:mapDisp%nSrc,1:2, 1:nz_dispersion) = 0.
          struct%polarCapAirMass(1:nz_dispersion) = real_missing

          ! the volume of the cap, needed for computing the fluxes out of it
          !
          if(name == 'S')then
            arc = (90.0 + fu_lat_native_from_grid(1.0, 0.5, dispersion_grid)) / 360.0
          elseif(name == 'N')then
            arc = (90.0 - fu_lat_native_from_grid(1.0, real(ny_dispersion)+0.5, dispersion_grid)) / 360.0
          else
            call set_error('Strange polar boundary name:'+name,'fu_makeBoundaryStruct')
            return
          endif
          rCap = arc * equator_length
          !
          ! Let's limit the arc length to preclude too small polar caps
          !
          if(arc < 2.8e-4)then   ! roughly 0.1 degree
            if(rCap < 0.4 * fu_dy_cell_m(dispersion_grid, 1, 1))then
              call set_error('Border:'+struct%name+'. Pole is too small:' + fu_str(arc*360.) + 'deg','fu_makeBoundaryStruct')
              return
            endif
          endif
          
          struct%polarCapArea = pi * rCap**2
         
          ! Nullify the unnecessary stuff

          nullify(struct%Fld3dPtr)
          nullify(struct%meteo2BoundaryInterpHoriz)
          nullify(struct%meteo2BoundaryInterpVert)
          nullify(struct%ifBoundSpecies)
          call set_missing(struct%vertical, .true.)
          allocate(struct%bc_wdr, stat=i)
          if(fu_fails(i == 0, 'Allocate failed', 'fu_make_boundaryStruct'))return
          struct%bc_wdr = wdr_missing

          struct%defined = silja_true
  
        case(dirichlet_boundary_type)


          if((.not. present(bgrid)).or.(.not. present(bvertical)))then
            call set_error('No boundary grid or vertical given for dirichlet structure', &
                         & 'fu_makeBoundaryStruct')
            return
          endif

          if((.not. present(fnm_multitime_tmpl)) .or. (.not. present(fnm_singletime_tmpl)) .or. &
           & (.not. present(fFmt_multitime)) .or. (.not. present(fFmt_singletime)))then
            call set_error('Problem with input filenames','fu_makeBoundaryStruct')
            return
          endif

          ! Actual boundary conditions given via concentration at upwind boudary
          !
          allocate(struct, stat=i)
          if(fu_fails(i == 0, 'Cannot allocate boundary structure','fu_makeBoundaryStruct'))return

          struct%name = name
          struct%boundaryType = dirichlet_boundary_type
          call msg('Making dirichlet boundary structure:' + struct%name)

         ! First create wdr with corresponding grid
         !
          allocate(struct%bc_wdr)   
          struct%bc_wdr = fu_set_wdr(1, &      ! nr of met srces
                                   & (/fu_set_met_src(bIndex, &       ! Index of the met_src
                                                    & centre_fmi, &   ! int_missing, &  centre
                                                    & bIndex, &   ! sub-centre but we force it to boundary index
                                                    & model_silam, &  ! int_missing, &  model
                                                    & rulesIniBoundary%bNames(bIndex))/), &
                                   & fnm_multitime_tmpl, &
                                   & fnm_singletime_tmpl, &
                                   & fFmt_multitime, &
                                   & fFmt_singletime, &
                                   & fu_top_level(wdr), fu_bottom_level(wdr), &
                                   & fu_start_time(wdr), &
                                   & rulesIniBoundary%bInterval, fu_period_to_compute(wdr), &
                                   & zero_interval, &       ! max hole length
                                   & fu_horizontal_interp_method(wdr), fu_vertical_interp_method(wdr), &
                                   & fu_time_interp_method(wdr), &
                                   & int_missing, int_missing, real_missing, int_missing, &   !  ablh_method, abl_param, minimal_ABLh, Kz_param
                                   & real_missing, &        ! PrecipitationLowLimit
                                   & .true., &              ! ifAllowZeroFcLen
                                   & fu_ifWait(wdr), &      ! ifWait
                                   & .false., &             ! ifDisregardMDS
                                   & storage_grid = bgrid)  ! storage_area, storage_grid (optional)

          !Boundary vertical
          !
          struct%vertical = bvertical

          ! Allocate the 3d-field arrays for the actual boundary values
          allocate(struct%ifBoundSpecies(nSpeciesTransp), struct%Fld3dPtr(2,nSpeciesTransp), stat=i)
          if(fu_fails(i == 0, 'Cannot allocate boundary structures','fu_makeBoundaryStruct'))return
          struct%nSp = nSpeciesTransp
          struct%nSrc = 1 !One source for Dirichlet boundaries
          struct%ifBoundSpecies(:) = .false.

         ! Allocate past and future field structure, fill with 0-s
          work_array => fu_work_array()
          work_array(:) = 0.

          do iS = 1, nSpeciesTransp

            allocate(struct%Fld3dPtr(1,iS)%fp, struct%Fld3dPtr(2,iS)%fp, stat = i)
            if(fu_fails(i == 0, 'Cannot allocate boundary 3d fields','fu_makeBoundaryStruct'))return
            call set_3d_field_empty(struct%Fld3dPtr(1,iS)%fp)
            call set_3d_field_empty(struct%Fld3dPtr(2,iS)%fp)

            do iZ = 1, fu_NbrOfLevels(bvertical)

              id = fu_set_field_id(met_src_missing,&
                                 & concentration_flag, &
                                 & rulesIniBoundary%start_time,&
                                 & zero_interval, &
                                 & bgrid,&
                                 & fu_level(bvertical, iZ, .true.),&
                                 & species = speciesTransp(iS))

              call create_field_in_field_3d(struct%Fld3dPtr(1,iS)%fp, id,work_array)
              if(error)return
              call create_field_in_field_3d(struct%Fld3dPtr(2,iS)%fp, id,work_array)
              if(error)return

!              field => fu_field_from_3d_field(struct%Fld3dPtr(1,iS)%fp, iZ, .true.) 
!              call set_field_empty(field, .true.)
!              call set_field(id,work_array,field)
!              field => fu_field_from_3d_field(struct%Fld3dPtr(2,iS)%fp, iZ, .true.)
!              call set_field_empty(field, .true.)
!              call set_field(id,work_array,field)
            enddo

            call set_3dField_params_from_fieldId(struct%Fld3dPtr(1,iS)%fp, id)
            call set_3dField_params_from_fieldId(struct%Fld3dPtr(2,iS)%fp, id)

!            if(struct%name /= 'T' .and. struct%name /= 'B')then
!              call organize_fields_vertically(struct%Fld3dPtr(1,iS)%fp)
!              call organize_fields_vertically(struct%Fld3dPtr(2,iS)%fp)
!            endif

          enddo
          call free_work_array(work_array)

          ! Now prepare the vertical interpolation structure
          !
          struct%meteo2BoundaryInterpVert => fu_vertical_interp_struct(meteo_vertical, &
                                                                     & bvertical, &
                                                                     & fu_storage_grid(struct%bc_wdr), &
                                                                     & linear, &
                                                                     & one_hour*3., &  ! ~update interval
                                                                     & 'meteo_to_boundary') ! name
          if(error)return
          
          ! Now prepare the horizontal interpolation structure
          !
          struct%meteo2BoundaryInterpHoriz => fu_horiz_interp_struct(meteo_grid, &
                                                                   & fu_storage_grid(struct%bc_wdr), &
                                                                   & linear, fu_if_randomise(wdr), &
                                                                   & iOutside = notAllowed)
          if(error)then
            call set_error('Failed boundary->diepersion horiz.interp coefs for:' + struct%name, &
                         & 'fu_makeBoundaryStruct')
            return
          endif

          struct%defined = silja_true

        case default 

          call set_error('Cannot create boundaries of this type yet.','fu_makeBoundaryStruct')
          return

      end select

    end function fu_makeBoundaryStruct

    !===============================================================================
  
    subroutine init_boundary_buffer(RulesIniBoundary, &
                                  & boundStructArray, &
                                  & BBuffer, &
                                  & nSrcTrn, iTransportSpeciesInvolved)
      !
      ! Allocates the boundary buffer and sets mapping between the boundary and dispoersion species
      !
      implicit none

      ! Imported parameters
      type(Tini_boundary_rules), intent(in) :: RulesIniBoundary
      type(TboundaryStructPtr), dimension(:), pointer :: boundStructArray
      type(TBoundaryBuffer), intent(out) :: BBuffer
      integer, intent(in) :: nSrcTrn
      integer, dimension(:), intent(inout) :: iTransportSpeciesInvolved

      ! Local parameters
      integer :: iB, nPolar, nDirichlet, iSpecies, nPoles

      !
      ! Understand the boundary structure 
      !
      bBuffer%iBoundaryType(1:6) = RulesIniBoundary%bTypes(1:6)
      nPolar = count(bBuffer%iBoundaryType(1:6) ==  polar_boundary_type)
      nDirichlet = count(bBuffer%iBoundaryType(1:6) == dirichlet_boundary_type)
      !
      ! Create the species maping: not all transport species can have boundary conditions
      ! First, count the max number of non-zero species - different borders can have different lists
      ! Second, allocate the index array
      ! Third, store the mapping
      !
      iSpecies = 0
      if(nDirichlet + nPolar > 0)then
        do iB = northern_boundary, bottom_boundary
          if(boundStructArray(iB)%ptr%boundaryType == dirichlet_boundary_type)then
            iSpecies = max(iSpecies,count(iTransportSpeciesInvolved(1:nSpeciesTransp) == 1))
          elseif(boundStructArray(iB)%ptr%boundaryType == polar_boundary_type)then
            iSpecies = nSpeciesTransp
          endif
        end do
      
        allocate(bBuffer%iTransportSpecies(iSpecies,6), bBuffer%iBoundarySpecies(nSpeciesTransp,6), &
               & bBuffer%wind_to_pole(nx_dispersion,2,nz_dispersion), bBuffer%outflowFactor(2,nz_dispersion),  &
               & stat=iB)
        if(fu_fails(iB==0,'Failed species mapping & wind storage allocation','init_boundary_buffer'))return
        bBuffer%iTransportSpecies(:,:) = int_missing
        bBuffer%iBoundarySpecies(:,:) = int_missing
        bBuffer%wind_to_pole(1:nx_dispersion,1:2,:) = 0.
        bBuffer%outflowFactor(:,:) = 1.
        !
        ! ... and now make the mapping itself
        !
        do iB = northern_boundary, bottom_boundary
          bBuffer%nBoundarySpecies(iB) = 0
          if(boundStructArray(iB)%ptr%boundaryType == dirichlet_boundary_type) then
            do iSpecies = 1, nSpeciesTransp
              if(iTransportSpeciesInvolved(iSpecies) == 1)then
                bBuffer%nBoundarySpecies(iB) = bBuffer%nBoundarySpecies(iB) + 1
                bBuffer%iTransportSpecies(bBuffer%nBoundarySpecies(iB),iB) = iSpecies
                bBuffer%iBoundarySpecies(iSpecies,iB) = bBuffer%nBoundarySpecies(iB)
              endif
            end do
          elseif(boundStructArray(iB)%ptr%boundaryType == polar_boundary_type) then
            bBuffer%nBoundarySpecies(iB) = nSpeciesTransp
            do iSpecies = 1, nSpeciesTransp
              bBuffer%iTransportSpecies(iSpecies,iB) = iSpecies
              bBuffer%iBoundarySpecies(iSpecies,iB) = iSpecies
            end do
          endif
        end do
      else
!        nullify(bBuffer%iTransportSpecies)
!        nullify(bBuffer%iBoundarySpecies)
        bBuffer%nBoundarySpecies(1:6) = 0
      endif  ! if Dirichlet boundaries exist
      
      !
      ! Allocate the buffer space. Note that space reserved only for those 
      ! boundaries that actually need it. Since nearly every boundary can be different from the 
      ! others, have to analyse them all one-by-one.
      !
      bBuffer%nSrcTrn = 1
      nPoles = 0
      select case(bBuffer%iBoundaryType(northern_boundary))                                          ! NORTH
        case(polar_boundary_type)
          bBuffer%bNorth => boundStructArray(northern_boundary)%ptr%polarCapMassMom
          bBuffer%PoleCapAirmass(northern_boundary)%pp => boundStructArray(northern_boundary)%ptr%polarCapAirmass
          bBuffer%fPoleCapArea(northern_boundary) = boundStructArray(northern_boundary)%ptr%polarCapArea
          bBuffer%indStep(northern_boundary) = 0
          bBuffer%nSrcTrn(northern_boundary) = nSrcTrn
          nPoles = 1                       ! enough to addresss northern_boundary

        case(dirichlet_boundary_type)
          allocate(bBuffer%bNorth(bBuffer%nBoundarySpecies(northern_boundary), &
                                 & 1, nx_dispersion, nz_dispersion), stat=iB)
          if(fu_fails(iB==0,'Failed northern buffer allocation','init_boundary_buffer'))return
          bBuffer%indStep(northern_boundary) = 1
        case default
          nullify(bBuffer%bNorth)
          bBuffer%indStep(northern_boundary) = int_missing
      end select

      select case(bBuffer%iBoundaryType(southern_boundary))                                         ! SOUTH
        case(polar_boundary_type)
          bBuffer%bSouth => boundStructArray(southern_boundary)%ptr%polarCapMassMom
          bBuffer%PoleCapAirmass(southern_boundary)%pp => boundStructArray(southern_boundary)%ptr%polarCapAirmass
          bBuffer%fPoleCapArea(southern_boundary) = boundStructArray(southern_boundary)%ptr%polarCapArea
          bBuffer%indStep(southern_boundary) = 0
          bBuffer%nSrcTrn(southern_boundary) = nSrcTrn
          nPoles = 2                       ! mandatory to reach southern_boundary

        case(dirichlet_boundary_type)
          allocate(bBuffer%bSouth(bBuffer%nBoundarySpecies(southern_boundary), &
                                 & 1, nx_dispersion, nz_dispersion), stat=iB)
          if(fu_fails(iB==0,'Failed southern buffer allocation','init_boundary_buffer'))return
          bBuffer%indStep(southern_boundary) = 1
        case default
          nullify(bBuffer%bSouth)
          bBuffer%indStep(southern_boundary) = int_missing
      end select
      !
      ! Top and bottom boundaries for poles have to be replicated from the neighbouring cells
      ! This will be done at the boundary preparation step. SO far, allocate the arrays
      !
      if(nPoles > 0)then
        iSPecies = max(bBuffer%nBoundarySpecies(northern_boundary), &
                     & bBuffer%nBoundarySpecies(southern_boundary))
        allocate( bBuffer%poleDryDep(nPoles, iSpecies, nSrcTrn), stat = iB)
        if(fu_fails(iB==0,'Failed allocation of pole dry dep','init_boundary_buffer'))return
        bBuffer%poleDryDep = 0.
      endif

      if(bBuffer%iBoundaryType(eastern_boundary) == dirichlet_boundary_type)then                    ! EAST
        allocate(bBuffer%bEast(bBuffer%nBoundarySpecies(eastern_boundary), &
                              & 1, ny_dispersion, nz_dispersion), stat=iB)
        if(fu_fails(iB==0,'Failed eastern buffer allocation','init_boundary_buffer'))return
        bBuffer%indStep(eastern_boundary) = 1
      else
!        nullify(bBuffer%bEast)
        bBuffer%indStep(eastern_boundary) = int_missing
      end if

      if(bBuffer%iBoundaryType(western_boundary) == dirichlet_boundary_type)then                    ! WEST
        allocate(bBuffer%bWest(bBuffer%nBoundarySpecies(western_boundary), &
                              & 1, ny_dispersion, nz_dispersion), stat=iB)
        if(fu_fails(iB==0,'Failed western buffer allocation','init_boundary_buffer'))return
        bBuffer%indStep(western_boundary) = 1
      else
!        nullify(bBuffer%bWest)
        bBuffer%indStep(western_boundary) = int_missing
      end if

      if(bBuffer%iBoundaryType(top_boundary) == dirichlet_boundary_type)then                        ! TOP
        allocate(bBuffer%bTop(bBuffer%nBoundarySpecies(top_boundary), &
                             & 1, nx_dispersion, ny_dispersion), stat=iB)
        if(fu_fails(iB==0,'Failed top buffer allocation','init_boundary_buffer'))return
        bBuffer%indStep(top_boundary) = 1
      else
!        nullify(bBuffer%bTop)
        bBuffer%indStep(top_boundary) = int_missing
      end if

      if(bBuffer%iBoundaryType(bottom_boundary) == dirichlet_boundary_type)then                     ! BOTTOM
        allocate(bBuffer%bBottom(bBuffer%nBoundarySpecies(bottom_boundary), &
                                & 1, nx_dispersion, ny_dispersion), stat=iB)
        if(fu_fails(iB==0,'Failed bottom buffer allocation','init_boundary_buffer'))return
        bBuffer%indStep(bottom_boundary) = 1
      else
!        nullify(bBuffer%bBottom)
        bBuffer%indStep(bottom_boundary) = int_missing
      end if

      bBuffer%pBoundaryData(northern_boundary)%ptr => bBuffer%bNorth
      bBuffer%pBoundaryData(southern_boundary)%ptr => bBuffer%bSouth
      bBuffer%pBoundaryData(eastern_boundary)%ptr => bBuffer%bEast
      bBuffer%pBoundaryData(western_boundary)%ptr => bBuffer%bWest
      bBuffer%pBoundaryData(top_boundary)%ptr => bBuffer%bTop
      bBuffer%pBoundaryData(bottom_boundary)%ptr => bBuffer%bBottom

    end subroutine init_boundary_buffer

  end subroutine init_boundary_structures


  !******************************************************************************

  SUBROUTINE fill_boundary_market(boundaryMarketPtr, ibr, first_step)
    !
    ! 1. Checks if the data on the shopping list are in the supermarket
    ! 2. If not, finds the necessary administrative data and gets the
    ! data in here.
    !
    IMPLICIT NONE

    ! Imported parameters with intent IN:
    type(mini_market_of_stacks), pointer :: boundaryMarketPtr
    type(Tini_boundary_rules), intent(inout), target :: ibr
    logical, intent(in) :: first_step

    ! Local declarations:
    TYPE(silja_time), DIMENSION(max_times) :: missing_obstimes, miss_obst_tmp
    INTEGER :: i, iTempl, iFNm, iFile, iStat,  iStack, iStack2fill, nStacks, iHeader
    type(silam_fformat) :: fformat
    type(meteo_data_source) :: met_src
    type(silam_sp), dimension(:), save, pointer :: fnames
    LOGICAL :: data_exists_already, ifNewDataHere, ifNewDataAnywhere
    type(grads_template), dimension(:), pointer :: fname_templates
    integer,dimension(7) :: indStacksToFill
    type(silja_shopping_list), pointer :: listPtr


    IF ((.NOT. minimarket_initialized(boundaryMarketPtr, multi_time_stack_flag)) .or. &
      & (.not. minimarket_initialized(boundaryMarketPtr, single_time_stack_flag))) THEN
      CALL set_error('Boundary market not initialized','fill_boundary_market')
      RETURN
    END IF

    indStacksToFill(1:7) = int_missing


    ! Cycle ofer boundary header files in ini-boundary_rules. 
    ! Each of those has separate set of input data files and separate shopping list.

    ! First check the necessary times for every file
    miss_obst_tmp(:) = time_missing
    ifNewDataAnywhere = .false.
    do iHeader = 1, size(ibr%bFiles)

      ! If static boundary, take it only once
      if(ibr%bFiles(iHeader)%ifStatic .and. .not. first_step)cycle

      IF (.NOT.defined(ibr%bFiles(iHeader)%shopLst)) THEN
        CALL set_error('Undefined shopping_list given','fill_boundary_market')
        RETURN
      END IF
      listPtr => ibr%bFiles(iHeader)%shopLst
      missing_obstimes(:) = time_missing
      if(ibr%bFiles(iHeader)%fieldType == monthly_climatology)then
        call set_list_time_indicator(ibr%bFiles(iHeader)%shopLst, accept_same_month)  ! 1 <-> iHeader
      endif


      CALL check_bm_for_times_in_list(boundaryMarketPtr, &
                                    & listPtr, &
                                    & data_exists_already,&
                                    & missing_obstimes, &
                                    & fu_obstime_interval(ibr))

      IF (error) RETURN
      miss_obst_tmp((iHeader-1)*10+1:iHeader*10) = missing_obstimes(1:10)
      
    enddo
     
    do iHeader = 1, size(ibr%bFiles)

      listPtr => ibr%bFiles(iHeader)%shopLst

      missing_obstimes(:) = time_missing
      missing_obstimes(1:10) = miss_obst_tmp((iHeader-1)*10+1:iHeader*10)

      FFormat = ibr%bFiles(iHeader)%bFileFormat


      !
      ! We have to go through all input files, for each go through all stacks to which its 
      ! data belong to.
      ! This means to go through all the input files, compilation of a list of stacks each file
      ! serves and then calling store_input_to_supermarket with this very list.
      !

      !----------------------------------------
      !
      ! 4. Template is allowed to have wildcards or imply test field or whatever.
      !    So, there may be several files for each time and each template. 
        !
      DO i = 1, SIZE(missing_obstimes)

        IF (.NOT.defined(missing_obstimes(i))) EXIT

        CALL fix_shopping_time_boundaries(listPtr,&
                                        & missing_obstimes(i)-one_minute,&
                                        & missing_obstimes(i)+one_minute)

        if(ibr%bFiles(iHeader)%fieldType == monthly_climatology)then
          call set_list_time_indicator(ibr%bFiles(1)%shopLst, accept_same_month)
        endif

        !
        ! Store the indices of the stacks to be filled in with the fields from files. 
        ! The same field can go to more than one boundary
        !

        iStack2fill = 1
        do iStack = 1, 6
          if(ibr%bFiles(iHeader)%ifToBoundary(iStack))then
            indStacksToFill(iStack2fill) = iStack
            iStack2fill = iStack2fill + 1
          endif
        enddo       
 
        !
        ! Have to scan templates and either create file name or forward the template string
        !
        do iFile = 1, size(ibr%bFiles(iHeader)%bFnameTemplate)
         
          if(ibr%bFiles(iHeader)%bFnameTemplate(iFile) == template_missing) exit

          if(fu_if_input_file(FFormat))then
            !
            ! File(s) is(are) implied - find their names
            !
            call FNm_from_single_template(ibr%bFiles(iHeader)%bFnameTemplate(iFile), &
                                        & missing_obstimes(i),&   
                                        & fnames, &
                                        & ifAdd = .false., &
                                        & ifStrict = .false., &
                                        & ifAllowZeroFcLen = .false., &
                                        & ifWait = .false., &
                                        & fc_step = fu_obstime_interval(ibr))
            IF (error) THEN
              CALL msg_warning('no data for this time, but continuing...')
              CALL unset_error('fill_boundary_market')
              exit !!! No file - no love
                   !! This way at least first missing file would pass
            END IF
          

          else
            !
            ! Template stores a request of some non-file input in its collection
            !
            call enlarge_array(fnames,1)
            fnames(1)%sp = fu_collection(ibr%bFiles(iHeader)%bFnameTemplate(iFile))
            if(size(fnames) > 1) fnames(2)%sp = ''
          endif  ! If this is a file

          do iFNm = 1, size(fnames)
            if(fnames(iFNm)%sp=='')exit

              CALL store_input_to_supermarket(boundaryMarketPtr, &
                                            & fnames(iFNm)%sp, &
                                            & FFormat, & 
                                            & listPtr, &
                                            & wdr_missing, &               !Tweak input if needed
                                            & add_field, &    ! ifUpdateAllowed
                                            & multi_time_stack_flag , 5, &
                                            & ifNewDataHere, &
                                            & indStacksToFill, &
                                            & storage_grid = grid_missing, &
                                            & iBinary = ibr%bFiles(iHeader)%unit_binary(iFNm), &
                                            & ifForceMetSrcAcceptance = .true., ifadjust = .true.)
          end do
        end do ! Template array
        ifNewDataAnywhere = .true.

      END DO ! Missing times
      
      IF (error) RETURN
    enddo ! cycle over boundary header files

        !
        ! Finally arrange supermarket, if at least something was found.
        !
    if (ifNewDataAnywhere) then
        call msg_warning('Why to arrange boundary market so many times?','fill_boundary_market')
        CALL arrange_supermarket(boundaryMarketPtr, .false.)
    endif

    !
    ! Maybe print some supermarket_info and exit.
    !
    !IF (supermarket_info) THEN
    IF (.true.) THEN
      ! missing_obstimes, i are used temporarily
      CALL supermarket_times(boundaryMarketPtr, met_src_missing, missing_obstimes, iFile) 
      IF (error) RETURN
      call msg('********************************************************')
      call msg('Data in supermarket now for the following times:')
      CALL report(missing_obstimes)
    END IF

  END SUBROUTINE fill_boundary_market


  !***********************************************************************************

  subroutine fillBoundaryStruct(ibRules, miniMarketBC, now, boundStructArray, &
                              & met_buf, ifFirstTime)
      !
      ! The sub is a bridge between the boundary mini-market and the boundary structure.
      ! In the mini-market the fields are in (up to) 6 stacks for (up to) 6 boudanries with
      ! corresponding grids (bands) but in the original vertical(s). Also, the species are 
      ! identical to transport species but quantities can be 
      ! mixing ratio, concentrations, masses or whatever else (in theory).
      ! The boundary structure keeps the fields exactly around the dispersion grid and 
      ! exactly as concentrations, exactly for this time step
      ! Thus, tasks here are:
      ! - interpolate the vertical from original to dispersion-like ones
      ! - make time interpolation
      ! - convert (diagnose) the input quantity to concentration
      !
      implicit none

      ! Imported parameters
      type(Tini_boundary_rules),  intent(inout), target :: ibRules
      type(TboundaryStructPtr), dimension(:), pointer :: boundStructArray
      type(silja_time ), intent(in) :: now !! Mid-time-step
      type(mini_market_of_stacks), pointer :: miniMarketBC
      type(Tfield_buffer), pointer :: met_buf
      logical, intent(in) :: ifFirstTime

      ! Local variables
      integer :: iBoundary, ifile, ifield 
      real :: weight_past_meteo

      type(silja_stack), pointer :: stackPtr
      type(silja_3d_field), pointer :: field3d
      type(silja_field), pointer :: field2d
      type(silja_field_id), pointer :: idPtr
      logical :: found
      type(TVertInterpStruct), pointer :: interp_met2bnd, interp_bnd2disp 
      integer :: past_fut_switch

      weight_past_meteo = met_buf%weight_past

      ! As different fields in same stack and also same fields in different stacks can originate from
      ! different input files and have different quantities and verticals, cycle over input files, going
      ! through each stack they serve (if2boundary array in bfile in rules) and taking fields originating
      ! from this input file (fieldInStack array in bfile in rules). Bfile structure also keeps the
      ! vertical for these input fields

      ! First time-step have to read the static fields to both past and future structures and also fill 
      ! the past structure with dynamic fields. Later only the future fields will be taken from stack, after 
      ! switching the pointer of past fields to previous future ones. 



      do iBoundary = 1, 6
        if(boundStructArray(iBoundary)%ptr%boundaryType == dirichlet_boundary_type)then

          ! Switch past and future
          !
          boundStructArray(iBoundary)%ptr%past_future_switch = &
                                                 & 1 - boundStructArray(iBoundary)%ptr%past_future_switch
        endif
      enddo


      ! Cycle over the input file sets. 
      !
      do iFile = 1, size(ibRules%bFiles)        ! sets of input files

        if(ibRules%bFiles(iFile)%ifStatic)then
          if(ifFirstTime)then   
            ! read the static field only at the first entrance
            ! their contence has to be put to both past and future fields in boundary structure, 
            !       
            do iBoundary = 1, 6

              if(.not. ibRules%bFiles(iFile)%ifToBoundary(iBoundary))cycle

              ! Stupidity checking
              !
              if(.not.fu_true(boundStructArray(iBoundary)%ptr%defined))then
                call msg('Undefined boundary structure:',iBoundary)
                call set_error('Undefined boundary structure','fillBoundaryStruct')
                return
              endif

              stackPtr => fu_closest_sm_met_src_time(miniMarketBC, & 
                                                   & fu_mds(boundStructArray(iBoundary)%ptr%bc_wdr,1), &
                                                   & time_missing, backwards)
              if(error)then
                call set_error('Failed to set the permanent 3D boundary','fillBoundaryStruct')
                return
              endif

              interp_met2bnd => ibRules%bFiles(iFile)%interpCoefMet2BndInp(iBoundary)%ptr

              interp_bnd2disp => ibRules%bFiles(iFile)%interpCoefBndInp2DispVert(iBoundary)%ptr

              !
              ! Preparatory work over, scan fields from the stack and send them to the boundary
              !
              do iField = 1, ibRules%bFiles(iFile)%nTSpecies

                call find_field_3d_from_stack(met_src_missing, &
                                            & fu_quantity(ibRules%bFiles(iFile)%fieldInStack(iField)),&
                                            & time_missing,&
                                            & stackPtr,&
                                            & Field3d,&
                                            & found, &
                                            & fu_species(ibRules%bFiles(iFile)%fieldInStack(iField)))
                if(found)then

                  ! set the future field in the structure
                  !
                  call field_3d_to_boundary_struct(Field3d, &                 
                                                 & boundStructArray(iBoundary)%ptr, &
                                                 & ibRules%bFiles(iFile)%tSpecies(iField), &
                                                 & interp_bnd2disp, & ! interp rule
                                                 & boundStructArray(iBoundary)%ptr%past_future_switch, &
                                                 &  met_buf, weight_past_meteo,ifFirstTime)
                  if(error)return

                  !point the past field to it
                  !
                  past_fut_switch = boundStructArray(iBoundary)%ptr%past_future_switch
                  boundStructArray(iBoundary)%ptr%Fld3dPtr(2-past_fut_switch, &
                                                         & ibRules%bFiles(iFile)%tSpecies(iField))%fp => &
                  & boundStructArray(iBoundary)%ptr%Fld3dPtr(1+past_fut_switch, &
                                                           & ibRules%bFiles(iFile)%tSpecies(iField))%fp
                else
                  !
                  ! Check 2D, if boundary type allows
                  !
                  if(boundStructArray(iBoundary)%ptr%name == 'T' .or. &
                   & boundStructArray(iBoundary)%ptr%name == 'B')then

                    call find_field_from_stack(stackPtr, ibRules%bFiles(iFile)%fieldInStack(iField), &
                                             & Field2d, found)
                    if(found)then
                      
                      call field_2d_to_boundary_struct(Field2d,&                 
                                                     & boundStructArray(iBoundary)%ptr, &
                                                     & ibRules%bFiles(iFile)%tSpecies(iField), &
                                                     & interp_bnd2disp, &  ! destination
                                                     & boundStructArray(iBoundary)%ptr%past_future_switch, &
                                                      &  met_buf, weight_past_meteo)
                      if(error)return

                      !point the past field to it
                      !
                      past_fut_switch = boundStructArray(iBoundary)%ptr%past_future_switch
                      boundStructArray(iBoundary)%ptr%Fld3dPtr(2-past_fut_switch, &
                                                             & ibRules%bFiles(iFile)%tSpecies(iField))%fp => &
                                               & boundStructArray(iBoundary)%ptr%Fld3dPtr(1+past_fut_switch, &
                                                           & ibRules%bFiles(iFile)%tSpecies(iField))%fp

                    endif ! if found in 2D

                  endif  ! if T or B boundary 

                endif  ! if input 3D field is in permanent stack
                if(.not. found)then
                  call set_error('Permanent boundary input field not found in stack', 'fillBoundaryStruct')
                  return
                endif
                boundStructArray(iBoundary)%ptr%ifBoundSpecies(ibRules%bFiles(iFile)%tSpecies(iField)) = .true.
              enddo  !field
            enddo  ! boundary

          else  ! if first time
            cycle
          endif ! if first time 
      
        else
          !
          ! not a Static file
          !
          do iBoundary = 1, 6

            if(.not. ibRules%bFiles(iFile)%ifToBoundary(iBoundary))cycle

            ! Stupidity checking
            !
            if(.not. boundStructArray(iBoundary)%ptr%defined==silja_true)then
              call msg('Undefined boundary structure:',iBoundary)
              call set_error('Undefined boundary structure','fillBoundaryStruct')
              return
            endif

            !
            ! Update the interpolation structure
            !
            interp_met2bnd => ibRules%bFiles(iFile)%interpCoefMet2BndInp(iBoundary)%ptr
            interp_bnd2disp => ibRules%bFiles(iFile)%interpCoefBndInp2DispVert(iBoundary)%ptr
            
! Handled in diagnostic_quantities for the whole pool     
!            call refine_interp_vert_coefs_v2(interp_bnd2disp, & 
!                                           & met_buf, &
!                                           & now, &
!                                           & boundStructArray(iBoundary)%ptr%meteo2BoundaryInterpHoriz, &
!                                           & interp_met2bnd)
!            if(error)return

            if(ifFirstTime)then
              !
              ! Switch temporarily past and future to fill the past stack
              !
              boundStructArray(iBoundary)%ptr%past_future_switch = &
                                             & 1 -  boundStructArray(iBoundary)%ptr%past_future_switch
              stackPtr => fu_closest_sm_met_src_time(miniMarketBC, &
                                                   & fu_mds(boundStructArray(iBoundary)%ptr%bc_wdr,1), &
                                                   & now, backwards)
              if(error)then
                call set_error('Failed to set the dynamic 3D boundary','fillBoundaryStruct')
                return
              endif

!              call msg('number of transport species', ibRules%bFiles(iFile)%nTSpecies)
              do iField = 1, ibRules%bFiles(iFile)%nTSpecies

                call find_field_3d_from_stack(met_src_missing, &
                                            & fu_quantity(ibRules%bFiles(iFile)%fieldInStack(iField)),&
                                            & time_missing,&
                                            & stackPtr,&
                                            & Field3d,&
                                            & found, &
                                            & fu_species(ibRules%bFiles(iFile)%fieldInStack(iField)))
                if(found)then
                  call field_3d_to_boundary_struct(Field3d,  &                 
                                                 & boundStructArray(iBoundary)%ptr, &
                                                 !species index in structure
                                                 & ibRules%bFiles(iFile)%tSpecies(iField), &
                                                 & interp_bnd2disp, & ! interp rule
                                                 & boundStructArray(iBoundary)%ptr%past_future_switch, &
                                                 &  met_buf, weight_past_meteo, ifFirstTime)
                  if(error)return
                else
                  !
                  ! Check 2D, if boundary type allows
                  !
                  if(boundStructArray(iBoundary)%ptr%name == 'T' .or. &
                   & boundStructArray(iBoundary)%ptr%name == 'B')then

                    call find_field_from_stack(stackPtr, ibRules%bFiles(iFile)%fieldInStack(iField), &
                                             & Field2d, found)
                    if(found)then
             
                      call field_2d_to_boundary_struct(Field2d, &                  
                                                     & boundStructArray(iBoundary)%ptr, &
                                                     & ibRules%bFiles(iFile)%tSpecies(iField), &
                                                     & interp_bnd2disp, &  ! interp rule
                                                     & boundStructArray(iBoundary)%ptr%past_future_switch, &
                                                      &  met_buf, weight_past_meteo)
                      if(error)return
                    endif  ! if boudnary in 2D
                  endif !if T or B
                endif
                if(.not.found)then
                  call msg("Field not found:")
                  idPtr => ibRules%bFiles(iFile)%fieldInStack(iField)
                  call report(idPtr)
                  call set_error('Dynamic boundary input field 1 not found in stack','fillBoundaryStruct')
                  return
                endif
                boundStructArray(iBoundary)%ptr%ifBoundSpecies(ibRules%bFiles(iFile)%tSpecies(iField)) = .true.
              enddo !field
              !
              ! Switch past and future back to fill the future stack
              !
              boundStructArray(iBoundary)%ptr%past_future_switch = &
                                           & 1 - boundStructArray(iBoundary)%ptr%past_future_switch
            endif  ! ifFirstTime



            ! And now to future fields

            !
            ! Find out the right stack
            !
            stackPtr => fu_closest_sm_met_src_time(miniMarketBC, &
                                                 & fu_mds(boundStructArray(iBoundary)%ptr%bc_wdr,1), &
                                                 & now, forwards)
            if(error)then

              call set_error('Failed to set the dynamic 3D boundary','fillBoundaryStruct')
              call unset_error('fillBoundaryStruct')
              call msg_warning('Continuing with last availible boundary fields', 'fillBoundaryStruct')
              do iField = 1, ibRules%bFiles(iFile)%nTSpecies
                past_fut_switch = boundStructArray(iBoundary)%ptr%past_future_switch
                boundStructArray(iBoundary)%ptr%Fld3dPtr(1+past_fut_switch, &
                 & ibRules%bFiles(iFile)%tSpecies(iField))%fp = &
                 & boundStructArray(iBoundary)%ptr%Fld3dPtr(2-past_fut_switch, &
                                                          & ibRules%bFiles(iFile)%tSpecies(iField))%fp
              enddo           
            else


              !
              ! Preparatory work over, scan the stack and send fields to the boundary structure
              !
              do iField = 1, ibRules%bFiles(iFile)%nTSpecies

                ! Now put the fields from stack to future field array
                !
                call find_field_3d_from_stack(met_src_missing, &
                                            & fu_quantity(ibRules%bFiles(iFile)%fieldInStack(iField)),&
                                            & time_missing,&
                                            & stackPtr,&
                                            & Field3d,&
                                            & found, &
                                            & fu_species(ibRules%bFiles(iFile)%fieldInStack(iField)))
                if(found)then

                  !
                  ! Write it down!
                  !

                  call field_3d_to_boundary_struct(Field3d,  &                 
                                                 & boundStructArray(iBoundary)%ptr, &
                                                 !species index in structure:
                                                 & ibRules%bFiles(iFile)%tSpecies(iField), & 
                                                 & interp_bnd2disp, & ! interp rule
                                                 & boundStructArray(iBoundary)%ptr%past_future_switch, &
                                                 &  met_buf, weight_past_meteo,ifFirstTime)
                  if(error)return
                else
                  !
                  ! Check 2D, if boundary type allows
                  !
                  if(boundStructArray(iBoundary)%ptr%name == 'T' .or. &
                   & boundStructArray(iBoundary)%ptr%name == 'B')then
  
                    call find_field_from_stack(stackPtr, ibRules%bFiles(iFile)%fieldInStack(iField), &
                                             & Field2d, found)
                    if(found)then
             
                      call field_2d_to_boundary_struct(Field2d, &                  
                                                     & boundStructArray(iBoundary)%ptr, &
                                                     & ibRules%bFiles(iFile)%tSpecies(iField), &
                                                     & interp_bnd2disp, &! interp rule
                                                     & boundStructArray(iBoundary)%ptr%past_future_switch, &
                                                      &  met_buf, weight_past_meteo)
                      if(error)return

                    endif  ! if boudnary in 2D
                  endif ! if T  or B
                endif ! if found in 3d
                if(.not. found)then
                  call msg("")
                  call msg("Field we were looking for:")
                  call report(ibRules%bFiles(iFile)%fieldInStack(iField))
                  call set_error('Dynamic boundary input field not found in stack',&
                               & 'fillBoundaryStruct')
                  return
                endif
              enddo
            endif  ! if error (stack missing)
          enddo
        endif ! ifStatic

      end do  !input file set
      call msg('Boundary structures ready')

  end subroutine fillBoundaryStruct


        !=======================================================================

        subroutine field_3d_to_boundary_struct(Field_3d, &        ! source field from stack
                                             & boundStruct, &     ! destination structure
                                             & iSpecies, &        ! species index in structure
                                          & interpStructVert, switch, met_buf, weight_past_meteo, ifFirstTime)  ! vertical interpolation rule  
          !
          ! Serves one boundary input field by distributing it from the boundary stack
          ! to the boundary structure
          ! Actions:
          ! - vertical interpolation
          ! - time interpolation
          !
          ! Points the future field ptr at the given field. 
          !
          implicit none

          ! Imported parameters

          type(TboundaryStruct), pointer :: boundStruct 
          type(TVertInterpStruct), pointer :: interpStructVert
          integer :: iSpecies
          TYPE(silja_3d_field), pointer :: Field_3d
       integer, intent(in) :: switch
       type(Tfield_buffer), pointer :: met_buf
       real, intent(in) :: weight_past_meteo
       logical, intent(in) :: ifFirstTime

          ! Local variables
       integer :: iLev, inputQuantity,  nx, ny, iTmp
          real, dimension(:), pointer :: pDat

          !
          ! Steps: 
          ! 1. Scan the 3D boundary structure grid collecting the data from input field via interpolation structure
          ! 2. Decide if conversion is needed and make it using the meteo fields picked via two more interpolation structures
          !    from meteo to boudnary horizontal and vertical. Note that meteo -> boundary vertical is in 4 cases
          !    out of 6 is equal to meteo -> dispersion vertical, so only pointers have to be set and no coefficient
          !    refinement is needed.

          ! 
          ! Interpolate the vertical
          !


          if(fu_valid_time(Field_3d) == fu_valid_time(boundStruct%Fld3dPtr(1+switch, iSpecies)%fp) .and. &
               & (.not. ifFirstTime))then
!            call msg('Field already in boundary structure')
!         call msg('Times in structure: '//fu_str(fu_valid_time(boundStruct%Fld3dPtr(1+switch, iSpecies)%fp)) //&
!                   & fu_str(fu_valid_time(boundStruct%Fld3dPtr(2-switch, iSpecies)%fp)) )
            return
          elseif(fu_valid_time(Field_3d) == fu_valid_time(boundStruct%Fld3dPtr(2-switch, iSpecies)%fp) .and. &
               & (.not. ifFirstTime))then
!            call msg('Slow field')
!         call msg("Field valid time:" //fu_str(fu_valid_time(Field_3d)))
!         call msg('Times in structure: '//fu_str(fu_valid_time(boundStruct%Fld3dPtr(1+switch, iSpecies)%fp)) //&
!                   & fu_str(fu_valid_time(boundStruct%Fld3dPtr(2-switch, iSpecies)%fp)) )


!            call report(Field_3d)
!            call report(boundStruct%Fld3dPtr(2-switch, iSpecies)%fp)
!            call report(boundStruct%Fld3dPtr(1+switch, iSpecies)%fp)


            Field_3d => boundStruct%Fld3dPtr(2-switch, iSpecies)%fp
            boundStruct%Fld3dPtr(2-switch, iSpecies)%fp => boundStruct%Fld3dPtr(1+switch, iSpecies)%fp
            boundStruct%Fld3dPtr(1+switch, iSpecies)%fp => Field_3d

!            call report(boundStruct%Fld3dPtr(2-switch, iSpecies)%fp)
!            call report(boundStruct%Fld3dPtr(1+switch, iSpecies)%fp)

!            call report(fu_valid_time(boundStruct%Fld3dPtr(1+switch, iSpecies)%fp))
!            call report(fu_valid_time(boundStruct%Fld3dPtr(2-switch, iSpecies)%fp))


            return
          
          else
!            call msg('setting new field to b structure')
!            call report(fu_valid_time(Field_3d))
!            call report(fu_valid_time(boundStruct%Fld3dPtr(1+switch, iSpecies)%fp))
          endif

          !
          ! Horizontal interpolation below should not take place, providing that the boundary fields all have
          ! the same grid. In some weird cases it can be violated, so better stay on the safe side.
          !
          if(.not.(fu_grid(Field_3d) == fu_grid(boundStruct%Fld3dPtr(1+switch, iSpecies)%fp)))then
            call msg('Field that comes for the boundary:')
            call report(Field_3d)
            call msg('Boundary structure field:')
            call report(boundStruct%Fld3dPtr(1+switch, iSpecies)%fp)
            call set_error('grids of input field and boundary structure do not correspond', &
                         & 'field_3d_to_boundary_struct')
            return
          endif

          call interp_field_3d_2_field_3d(Field_3d, &                                     ! input
                                        & boundStruct%Fld3dPtr(1+switch, iSpecies)%fp, &  ! output
!                                        & .not.(fu_grid(Field_3d) == &
!                                              & fu_grid(boundStruct%Fld3dPtr(1+switch, iSpecies)%fp)), &  ! if horizontal interp.
                                        & .false., &                                ! if horizontal interp
                                        & .true., &                                 ! if vertical interp
                                        & boundStruct%meteo2BoundaryInterpHoriz, &  ! horizontal interp structure
                                        & interpStructVert)                         ! vertical interp structure
          if(error)return

          call set_valid_time(boundStruct%Fld3dPtr(1+switch, iSpecies)%fp, fu_valid_time(Field_3d))
          if(error)return
!call msg('switch, species', switch, iSpecies)
!call report(fu_valid_time(Field_3d))
!call report(fu_valid_time(boundStruct%Fld3dPtr(1+switch, iSpecies)%fp))

          !
          ! Check whether the quantity conversion is needed
          !
          inputQuantity = fu_quantity(Field_3d)
          if(inputQuantity /= volume_mixing_ratio_flag)then
            !
            ! Have to diagnose from
            !
            if(inputQuantity == concentration_flag)then
              !
              ! Convert
              !
              call make_vmr_from_cnc(met_buf, weight_past_meteo, &
                                   & boundStruct%meteo2BoundaryInterpHoriz, &
                                   & boundStruct%meteo2BoundaryInterpVert, &
                                   & .true., .true., &
                                   & boundStruct%Fld3dPtr(1+switch, iSpecies)%fp)
            else
              call set_error('Concentration cannot be diagnosed from input quantity:' + &
                           & fu_quantity_string(inputQuantity),'field_3d_to_boundary_struct')
              return
            endif
          endif  ! if conversion is needed


        end subroutine field_3d_to_boundary_struct

        !=======================================================================

        subroutine field_2d_to_boundary_struct(Field_2d, &        ! source field from stack
                                             & boundStruct, &     ! destination structure
                                             & iSpecies, &        ! species index in structure
                                          & interpStructVert, switch, met_buf, weight_past_meteo)  ! vertical interpolation rule  
          !
          ! Serves one boundary input field by distributing it from the boundary stack
          ! to the boundary structure
          ! Actions:
          ! - vertical interpolation
          ! - time interpolation
          !
          ! Points the future field ptr at the given field. 
          !
          implicit none

          ! Imported parameters

          type(TboundaryStruct), pointer :: boundStruct 
          type(TVertInterpStruct), pointer :: interpStructVert
       integer, intent(in) :: iSpecies
          TYPE(silja_field), pointer :: Field_2d
       integer, intent(in) :: switch
       type(Tfield_buffer), pointer :: met_buf
       real, intent(in) :: weight_past_meteo

          ! Local variables
       integer :: iLev, inputQuantity


          !
          ! Steps: 
          ! 1. Scan the 3D boundary structure grid collecting the data from input field via 
          !    interpolation structure
          ! 2. Decide if conversion is needed and make it using the meteo fields picked via 
          !    two more interpolation structures from meteo to boudnary horizontal and 
          !    vertical. Note that meteo -> boundary vertical is in 4 cases out of 6 is equal 
          !    to meteo -> dispersion vertical, so only pointers have to be set and no coefficient
          !    refinement is needed.

          ! 
          ! Only for top and bottom boundaries - one level structure
          ! No vertical interpolation, the field goes to the only level of the structure as it is
          ! First set the 3d field empty, then add the new field
          !
          call set_3d_field_empty(boundStruct%Fld3dPtr(1+switch, iSpecies)%fp)
          call add_field_to_3d_field(Field_2d, boundStruct%Fld3dPtr(1+switch, iSpecies)%fp)
   
                   
          ! Check whether the quantity conversion is needed
          !
          inputQuantity = fu_quantity(Field_2d)
          if(inputQuantity /= volume_mixing_ratio_flag)then
            !
            ! Have to diagnose. So far, only volume mixing ratio is known. Here we call directly
            ! but eventually should make it through the generic diagnostic procedure
            !
            if(inputQuantity == concentration_flag)then
              !
              ! Convert
              !
              call msg('Making VMR from cnc...')
              call make_vmr_from_cnc(met_buf, weight_past_meteo, &
                                   & boundStruct%meteo2BoundaryInterpHoriz, &
                                   & boundStruct%meteo2BoundaryInterpVert, &
                                   & .true., .true., &
                                   & boundStruct%Fld3dPtr(1+switch, iSpecies)%fp)
            else
              call set_error(fu_connect_strings('Concentration cannot be diagnosed from input quantity:', &
                                              & fu_quantity_string(inputQuantity)), &
                          & 'field_3d_to_boundary_struct')
              return
            endif
          endif  ! if conversion is needed


        end subroutine field_2d_to_boundary_struct


        !======================================================================

        subroutine make_cnc_from_vmr(met_buf, weight_past, &
                                   & pHorizInterpStruct, pVertInterpStruct, &
                                   & ifHorizInterp, ifVertInterp, &
                                   & vmr_fld_3d)

          ! Creates concentrations from volume mixing ratio
          ! It is assumed, that the concentration field 3d exists and whatever is there, 
          ! is overwritten
          ! 
          implicit none

          ! Imported parameters
          type(Tfield_buffer), intent(in) :: met_buf
          real, intent(in) :: weight_past
          type(THorizInterpStruct), pointer :: pHorizInterpStruct
          type(TVertInterpStruct), pointer :: pVertInterpStruct
          logical, intent(in) :: ifHorizInterp, ifVertInterp
          type(silja_3d_field), pointer :: vmr_fld_3d
          type(silja_field_id), pointer :: id

          ! Local declarations:
          REAL, DIMENSION(:), POINTER :: volume_mixing_ratio, concentration
          real :: temperature, pressure
          INTEGER :: tIndex, pIndex, ix, iy, iz, iCell, nx, ny 

          ! A bit of preparation...
          !
          call grid_dimensions(fu_grid(vmr_fld_3d), nx, ny)

          tIndex = fu_index(met_buf%buffer_quantities, temperature_flag)
          pIndex = fu_index(met_buf%buffer_quantities, pressure_flag)
          do iz = 1, fu_number_of_fields(vmr_fld_3d)
            !call report(fu_level(fu_field_from_3d_field(vmr_fld_3d, iz)))
            volume_mixing_ratio => fu_grid_data_from_3d(vmr_fld_3d, iz)
            concentration => volume_mixing_ratio
            iCell = 0
            do iy = 1, ny
              do ix = 1, nx
                iCell = iCell + 1

                ! Get temperature and pressure
                ! 
                temperature = fu_get_value(met_buf%p4d(tIndex), nx_meteo, ix, iy, iz, &
                                         & weight_past, &
                                         & pHorizInterpStruct, pVertInterpStruct, &
                                         & ifHorizInterp, ifVertInterp)
                IF (error) RETURN
                pressure = fu_get_value(met_buf%p4d(pIndex), nx_meteo, ix, iy, iz, &
                                      & weight_past, &
                                      & pHorizInterpStruct, pVertInterpStruct, &
                                      & ifHorizInterp, ifVertInterp)
                IF (error) RETURN
                !print *, 'ilev, vmr', iz, temperature, pressure, concentration(icell)
                ! Compute the concentration
                !
                concentration(iCell) = volume_mixing_ratio(iCell) * pressure / &
                                     & ( gas_constant_uni * temperature)
                
 
              enddo   ! ix
            enddo  ! iy
            !
            ! Conversion is over, set the new quantity to the field
            !
            id => fu_id(fu_field_from_3d_field(vmr_fld_3d, iz))
            call set_quantity(id, concentration_flag)

          enddo  ! iz
        END SUBROUTINE make_cnc_from_vmr

        !======================================================================

        subroutine make_vmr_from_cnc(met_buf, weight_past, &
                                   & pHorizInterpStruct, pVertInterpStruct, &
                                   & ifHorizInterp, ifVertInterp, &
                                   & vmr_fld_3d)

          ! Creates volume mixing ratio from concentration
          ! It is assumed, that the vmr field 3d exists and whatever is there, 
          ! is overwritten
          ! 
          implicit none

          ! Imported parameters
          type(Tfield_buffer), intent(in) :: met_buf
          real, intent(in) :: weight_past
          type(THorizInterpStruct), intent(in) :: pHorizInterpStruct
          type(TVertInterpStruct), intent(in) :: pVertInterpStruct
          logical, intent(in) :: ifHorizInterp, ifVertInterp
          type(silja_3d_field), pointer :: vmr_fld_3d

          ! Local declarations:
          REAL, DIMENSION(:), POINTER :: vals
          type(silja_field_id), pointer :: id
          real :: temperature, pressure
          INTEGER :: tIndex, pIndex, ix, iy, iz, iCell, nx, ny 

          ! A bit of preparation...
          !
          call grid_dimensions(fu_grid(vmr_fld_3d), nx, ny)

          tIndex = fu_index(met_buf%buffer_quantities, temperature_flag)
          pIndex = fu_index(met_buf%buffer_quantities, pressure_flag)
          do iz = 1, fu_number_of_fields(vmr_fld_3d)
            !call report(fu_level(fu_field_from_3d_field(vmr_fld_3d, iz)))
            vals => fu_grid_data_from_3d(vmr_fld_3d, iz)
            iCell = 0
            do iy = 1, ny
              do ix = 1, nx
                iCell = iCell + 1
                temperature = fu_get_value(met_buf%p4d(tIndex), nx_meteo, ix, iy, iz, &
                                         & weight_past, &
                                         & pHorizInterpStruct, pVertInterpStruct, &
                                         & ifHorizInterp, ifVertInterp)
                IF (error) RETURN
                pressure = fu_get_value(met_buf%p4d(pIndex), nx_meteo, ix, iy, iz, &
                                      & weight_past, &
                                      & pHorizInterpStruct, pVertInterpStruct, &
                                      & ifHorizInterp, ifVertInterp)
                IF (error) RETURN
                vals(iCell) = vals(iCell) * ( gas_constant_uni * temperature) / pressure
              enddo   ! ix
            enddo  ! iy
            !
            ! Conversion is over, set the new quantity to the field
            !
            id => fu_id(fu_field_from_3d_field(vmr_fld_3d, iz))
            call set_quantity(id, volume_mixing_ratio_flag)

          enddo  ! iz
        END SUBROUTINE make_vmr_from_cnc



  ! ***************************************************************

  subroutine boundary_conditions_now(arBoundaries, bBuffer, now)
    !
    ! This subroutine packs the boundary conditions into easily addressable buffer
    !! Also converts SILAM vmr (mass_unit / mole_of_air) into 
    !! mixing ratio needed for advection (mass_unit / kg_of_air)
    ! Copies and interpolates from arBoundaries to bBuffer
    !
    implicit none

    ! Imported parameters
    type(TboundaryStructPtr), dimension(:), pointer :: arBoundaries
    type(TboundaryBuffer), pointer :: bBuffer
    type(silja_time), intent(in) :: now
!    type(Tfield_buffer), intent(in) :: disp_buf

    ! Local parameters
    type(TboundaryStruct), pointer :: bnd
    
    !
    ! Scan the boundaries one by one and get the pointers and values of the data
    ! Note that polar hat will be modified by advection etc, so pointer makes sense
    !
    bnd => arBoundaries(northern_boundary)%ptr                                               ! NORTH
    select case(bnd%boundaryType)
      case(polar_boundary_type)
        bBuffer%bNorth => bnd%polarCapMassMom
        
      case(dirichlet_boundary_type)
        call get_dirichlet_boundary(bnd, &
                                  & bBuffer%bNorth(1:bBuffer%nBoundarySpecies(northern_boundary), &
                                                 & 1:1,1:nx_dispersion,1:nz_dispersion), &
                                  & bBuffer%iTransportSpecies(:,northern_boundary), &
                                  & bBuffer%nBoundarySpecies(northern_boundary), &
                                  & now, &
                                  & 1, nx_dispersion, ny_dispersion, ny_dispersion, 1, nz_dispersion, &
                                  & .false.)
      case default
    end select

    bnd => arBoundaries(southern_boundary)%ptr                                               ! SOUTH
    select case(bnd%boundaryType)
      case(polar_boundary_type)
        bBuffer%bSouth => bnd%polarCapMassMom

      case(dirichlet_boundary_type)
        call get_dirichlet_boundary(bnd, &
                                  & bBuffer%bSouth(1:bBuffer%nBoundarySpecies(southern_boundary), &
                                                 & 1:1,1:nx_dispersion,1:nz_dispersion), &
                                  & bBuffer%iTransportSpecies(:,southern_boundary), &
                                  & bBuffer%nBoundarySpecies(southern_boundary), &
                                  & now, &
                                  & 1,nx_dispersion, 1, 1, 1, nz_dispersion, &
                                  & .false.)
      case default
    end select

    bnd => arBoundaries(eastern_boundary)%ptr                                               ! EAST
    if(bnd%boundaryType == dirichlet_boundary_type) &
        & call get_dirichlet_boundary(bnd, &
                                    & bBuffer%bEast(1:bBuffer%nBoundarySpecies(eastern_boundary), &
                                                  & 1:1, 1:ny_dispersion, 1:nz_dispersion), &
                                    & bBuffer%iTransportSpecies(:,eastern_boundary), &
                                    & bBuffer%nBoundarySpecies(eastern_boundary), &
                                    & now, &
                                    & nx_dispersion, nx_dispersion, 1, ny_dispersion, 1, nz_dispersion, &
                                    & .false.)

    bnd => arBoundaries(western_boundary)%ptr                                               ! WEST
    if(bnd%boundaryType == dirichlet_boundary_type) &
        & call get_dirichlet_boundary(bnd, &
                                    & bBuffer%bWest(1:bBuffer%nBoundarySpecies(western_boundary), &
                                                  & 1:1, 1:ny_dispersion, 1:nz_dispersion), &
                                    & bBuffer%iTransportSpecies(:,western_boundary), &
                                    & bBuffer%nBoundarySpecies(western_boundary), &
                                    & now, &
                                    & 1, 1, 1, ny_dispersion, 1, nz_dispersion, &
                                    & .false.)

    bnd => arBoundaries(top_boundary)%ptr                                                   ! TOP
    if(bnd%boundaryType == dirichlet_boundary_type) &
        & call get_dirichlet_boundary(bnd, &
                                    & bBuffer%bTop(1:bBuffer%nBoundarySpecies(top_boundary), &
                                                 & 1:1, 1:nx_dispersion, 1:ny_dispersion), &
                                    & bBuffer%iTransportSpecies(:,top_boundary), &
                                    & bBuffer%nBoundarySpecies(top_boundary), &
                                    & now, &
                                    & 1, nx_dispersion, 1, ny_dispersion, int_missing, int_missing, &
                                    & .true.)

    bnd => arBoundaries(bottom_boundary)%ptr                                                ! BOTTOM
    if(bnd%boundaryType == dirichlet_boundary_type) &
        & call get_dirichlet_boundary(bnd, &
                                    & bBuffer%bBottom(1:bBuffer%nBoundarySpecies(bottom_boundary), &
                                                  & 1:1, 1:nx_dispersion, 1:ny_dispersion), &
                                    & bBuffer%iTransportSpecies(:,bottom_boundary), &
                                    & bBuffer%nBoundarySpecies(bottom_boundary), &
                                    & now, &
                                    & 1, nx_dispersion, 1, ny_dispersion, int_missing, int_missing, &
                                    & .true.)

    CONTAINS

      !===================================================================
      
      subroutine get_dirichlet_boundary(bnd, buf_data, species_mapping_trn, nSpecies, now, &
                                      & x1, x2, y1, y2, z1, z2, ifPlain)
        !
        ! Makes-up the Dirichlet boundary: actually, just applies time interpolation
        !
        implicit none
        ! Imported parameters
        type(TboundaryStruct), pointer:: bnd
        real, dimension(:,:,:,:), intent(out) :: buf_data  ! (x/y, species, src, z)
        integer, dimension(:), intent(in) :: species_mapping_trn
        integer, intent(in) :: x1, x2, y1, y2, z1, z2, nSpecies
        logical, intent(in) :: ifPlain
        type(silja_time), intent(in) :: now  !!! Should be mid-step here
        
        ! Local variables
        integer :: ix, iy, iz, iCell, izb, iSpecies, iSpeciesTrn !, i1b
        real, dimension(max_species) :: weight_past_boundary
        real, dimension(:), pointer :: val_past, val_future !, &


        !
        ! Get the weight_past_boundary, which can be species-specific
        !
        do iSpecies = 1, nSpecies
          iSpeciesTrn = species_mapping_trn(iSpecies)
          if(fu_valid_time(bnd%Fld3dPtr(2-bnd%past_future_switch, iSpeciesTrn)%fp) == &
           & fu_valid_time(bnd%Fld3dPtr(1+bnd%past_future_switch, iSpeciesTrn)%fp))then
            weight_past_boundary(iSpeciesTrn) = 1.0
          elseif(now >= fu_valid_time(bnd%Fld3dPtr(1+bnd%past_future_switch, iSpeciesTrn)%fp))then
            weight_past_boundary(iSpeciesTrn) = 0.0
          else
            weight_past_boundary(iSpeciesTrn) = &
                   & (fu_valid_time(bnd%Fld3dPtr(1+bnd%past_future_switch, iSpeciesTrn)%fp) - now) / &
                   & (fu_valid_time(bnd%Fld3dPtr(1+bnd%past_future_switch, iSpeciesTrn)%fp) - &
                    & fu_valid_time(bnd%Fld3dPtr(2-bnd%past_future_switch, iSpeciesTrn)%fp))
          endif
          if(weight_past_boundary(iSpeciesTrn) < 0. .or. weight_past_boundary(iSpeciesTrn) > 1.)then
            call msg('Strange weight_past_boundary:', weight_past_boundary(iSpeciesTrn))
            call msg('Now and valid times for 1+switch and 2-switch:')
            call report(now)
            call report(fu_valid_time(bnd%Fld3dPtr(1+bnd%past_future_switch, iSpeciesTrn)%fp))
            call report(fu_valid_time(bnd%Fld3dPtr(2-bnd%past_future_switch, iSpeciesTrn)%fp))
          endif
          if(error)return
        end do  ! iSpecies

        !
        ! The grand injection cycle
        !
        if(ifPlain)then
          !
          ! top or bottom plain boundary: z cycle is void, x and y vary
          !
          do iSpecies = 1, nSpecies
            iSpeciesTrn = species_mapping_trn(iSpecies)
            if(.not. bnd%ifBoundSpecies(iSpeciesTrn))cycle
            val_future => fu_grid_data_from_3d(bnd%Fld3dPtr(1 + bnd%past_future_switch, &
                                                          & iSpeciesTrn)%fp, 1)
            val_past => fu_grid_data_from_3d(bnd%Fld3dPtr(2 - bnd%past_future_switch, &
                                                        & iSpeciesTrn)%fp, 1)
            iCell = 1
            do iy = y1,y2
              do ix = x1,x2
                buf_data(iSpecies,1,ix,iy) = &
                                   & weight_past_boundary(iSpeciesTrn) * val_past(iCell) + &
                                   & (1. - weight_past_boundary(iSpeciesTrn)) * val_future(iCell)
#ifdef DEBUG
                  if (.not.  buf_data(iSpecies,1,ix,iy) >= 0) then
                          call msg ("")
                          call msg ("Negative came from boundaries Plain")
                          call report(now)
                          call msg("Icell "//trim(fu_str(iCell)), val_past(iCell), val_future(iCell))
                          call msg("iSpeciesTrn,iz",iSpeciesTrn,iz)
                          call report(bnd%Fld3dPtr(1 + bnd%past_future_switch, iSpeciesTrn)%fp)
                          call report(bnd%Fld3dPtr(2 - bnd%past_future_switch, iSpeciesTrn)%fp)
                  endif
#endif                  
                iCell = iCell + 1  ! counts move along the plain
              end do !x
            enddo !y
          end do !Species
        else
          !
          ! Vertical walls: either x or y is void, the other and z vary
          !
          do iz = z1, z2
!            z_cell_size_past => pDispBuf%p4d(zSize_ind)%past%p2d(ilev)%ptr
!            z_cell_size_future => pDispBuf%p4d(zSize_ind)%future%p2d(ilev)%ptr
            do iSpecies = 1, nSpecies
              iSpeciesTrn = species_mapping_trn(iSpecies)
              if(.not. bnd%ifBoundSpecies(iSpeciesTrn))cycle
              val_future => fu_grid_data_from_3d(bnd%Fld3dPtr(1 + bnd%past_future_switch, &
                                                            & iSpeciesTrn)%fp, iz)
              val_past => fu_grid_data_from_3d(bnd%Fld3dPtr(2 - bnd%past_future_switch, &
                                                          & iSpeciesTrn)%fp, iz)
              iCell = 1
              do iy = y1,y2
                do ix = x1,x2
                  !
                  ! Even if zeroes: faster than .eps. 0. and buf_data should be nullified anyway
                  !
!                  i1d = ix + (iy - 1) * nx_dispersion
                  buf_data(iSpecies,1,iCell,iz) = &
                                 & (weight_past_boundary(iSpeciesTrn) * val_past(iCell) + &
                                  & (1. - weight_past_boundary(iSpeciesTrn)) * val_future(iCell)) 
#ifdef DEBUG
                  if (.not.  buf_data(iSpecies,1,iCell,iz) >= 0) then
                          call msg ("")
                          call msg ("Negative came from boundaries")
                          call report(now)
                          call msg("Icell "//trim(fu_str(iCell)), val_past(iCell), val_future(iCell))
                          call msg("iSpeciesTrn,iz",iSpeciesTrn,iz)
                          call report(bnd%Fld3dPtr(1 + bnd%past_future_switch, iSpeciesTrn)%fp)
                          call report(bnd%Fld3dPtr(2 - bnd%past_future_switch, iSpeciesTrn)%fp)
                  endif
#endif                  
                  iCell = iCell + 1  ! counts move along x/y
                end do !ix
              end do !iy
            end do ! species
          enddo ! iz
        end if ! if wall or plane

       !!
       !! Finally make the buffer to mass mixing ratio needed for advection
      buf_data(:,:,:,:) = buf_data(:,:,:,:) * (1./ molecular_weight_air)

#ifdef DEBUG
        if (.not. all(buf_data(:,:,:,:)>=0)) then
                call set_error("Non-positive in buf_data", "boundary_conditions_now")
        endif
#endif
        
      end subroutine get_dirichlet_boundary

  end subroutine boundary_conditions_now


  ! ***************************************************************

    function fu_ifBoundary(ibRules)
      implicit none
      type(Tini_boundary_rules), intent(in) :: ibRules
      logical :: fu_ifBoundary

      fu_ifBoundary = ibRules%ifBoundary

    end function fu_ifBoundary
  ! ***************************************************************

    function fu_ifReadBoundaries(bcRules)
      implicit none
      type(Tini_boundary_rules), intent(in) :: bcRules
      logical :: fu_ifReadBoundaries
      integer :: iBoundary
      fu_ifReadBoundaries = .false.

      if(bcRules%ifBoundary)then
        do iBoundary = 1, 6
          if(bcRules%bTypes(iBoundary) == dirichlet_boundary_type)then
            fu_ifReadBoundaries = .true.  
            exit   
          endif
        enddo
      endif

    end function fu_ifReadBoundaries

  ! ***************************************************************


  FUNCTION fu_obstime_int_bcr(bcr) 
    IMPLICIT NONE
    TYPE(silja_interval) :: fu_obstime_int_bcr
    TYPE(Tini_boundary_rules), INTENT(in) :: bcr
    fu_obstime_int_bcr = bcr%bInterval
  END FUNCTION fu_obstime_int_bcr
  
  
  !***********************************************************************************

  function fu_boundaryStructDefined(bStruct) result(ifDef)
    implicit none
    type(TboundaryStruct),  intent(in) :: bStruct
    logical :: ifDef
    ifDef = fu_true(bStruct%defined)
  end function fu_boundaryStructDefined


  !***********************************************************************************
 
   function fu_nr_boundary_inFiles(bcRules)
       IMPLICIT NONE
    TYPE(Tini_boundary_rules), INTENT(in) :: bcRules
    integer :: fu_nr_boundary_inFiles

    fu_nr_boundary_inFiles = size(bcRules%bFiles) 

   end function fu_nr_boundary_inFiles


  !***********************************************************************************
 
   function fu_shplst(bcRules, i)
       IMPLICIT NONE

    TYPE(Tini_boundary_rules), intent(in), target :: bcRules
    integer :: i
    type(silja_shopping_list), pointer ::  fu_shplst

    fu_shplst =>  bcRules%bFiles(i)%shopLst

   end function fu_shplst


  !***********************************************************************************

  function fu_ifHaveInitialConditions(bc_rules) result(ifTrue)
    implicit none
    type(Tini_boundary_rules), intent(in) :: bc_rules
    
    logical :: ifTrue
    
    ifTrue = associated(bc_rules%qInitial)
 
  end function fu_ifHaveInitialConditions

  
  !*******************************************************************************************
  
  subroutine force_mass_map_cell_mmr(mapConc, indSpecies, pBoundaryBuffer, disp_buf)
    !
    ! Sets massmap from cellmasses, initializes airmass in pole
    ! Note that it sets mass mixing ratio and thus needs air mass in the cell.
    ! Do not call before that field is available in dispersion buffer
    !
    implicit none

    type(TMass_map), intent(inout) :: mapConc
    type(TboundaryBuffer), intent(inout) :: pBoundaryBuffer
    TYPE(Tfield_buffer), intent(in) :: disp_buf
    integer, intent(in) :: indSpecies

    !Local
    real :: fTmp, fMass
    integer :: m_ind, ix, iy, iz, i1d, nx, ny, nz
    real, dimension (:), pointer :: pXCellSize, pYCellSize


    nx=mapConc%nx
    ny=mapConc%ny
    nz=mapConc%n3d

    m_ind = fu_index(disp_buf%buffer_quantities, disp_cell_airmass_flag)
    fTmp = 1. !MMR in mol/kg to force

    call msg('Resetting mass map with '//trim(fu_str(fTmp))// ' mol/kg cnc')
    pXCellSize => fu_grid_data(dispersion_cell_x_size_fld)
    pYCellSize => fu_grid_data(dispersion_cell_y_size_fld)
    do iz = 1,nz
      do iy =1,ny
           i1d = (iy-1)*nx+1
           mapConc%arm(indSpecies,1,iz,:,iy) = disp_buf%p4d(m_ind)%past%p2d(iz)%ptr(i1d:i1d+nx-1)*fTmp
      enddo


      if(pBoundaryBuffer%iBoundaryType(northern_boundary) == polar_boundary_type)then
            i1d  = nx*(ny-1)+1 ! Non-staggered index for masses/cell_sizes
            ! Average area density of air in the layer (kg/m2) around * plar_cap_aro $a
            fMass =     sum(disp_buf%p4d(m_ind)%past%p2d(iz)%ptr(i1d:i1d+nx-1))
            fMass = fMass * pBoundaryBuffer%fPoleCapArea(northern_boundary)/ &
                                 & (pYCellSize(i1d)*pXCellSize(i1d)*nx)  
            pBoundaryBuffer%PoleCapAirmass(northern_boundary)%pp(iz) = fMass  ! kg of air
            pBoundaryBuffer%bNorth(indSpecies,1,1,iz) = fMass*fTmp
            pBoundaryBuffer%bNorth(indSpecies,1,2,iz) = 0.
            pBoundaryBuffer%wind_to_pole(1:nx,northern_boundary,:) = 0

       endif

       if(pBoundaryBuffer%iBoundaryType(southern_boundary) == polar_boundary_type)then
            ! Average area density of air in the layer (kg/m2) around * plar_cap_area
            fMass =     sum(disp_buf%p4d(m_ind)%past%p2d(iz)%ptr(1:nx))
            fMass = fMass * pBoundaryBuffer%fPoleCapArea(southern_boundary)/ & 
                                 &( pYCellSize(1)*pXCellSize(1)*nx)  
            pBoundaryBuffer%PoleCapAirmass(southern_boundary)%pp(iz) = fMass  ! kg of air
            pBoundaryBuffer%bSouth(indSpecies,1,1,iz) = fMass*fTmp
            pBoundaryBuffer%bSouth(indSpecies,1,2,iz) = 0.
            pBoundaryBuffer%wind_to_pole(1:nx,southern_boundary,:) = 0
       endif
    enddo
    call msg('Resetting mass map done')

  end subroutine force_mass_map_cell_mmr

  
END MODULE ini_boundary_conditions




