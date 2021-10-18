module io_server
  !
  ! This module handles the input/output features of SILAM. The main tasks are:
  !
  ! 1. to read the list of requested output variables from the output 
  !    configuration file and transform them to some reasonable list of model
  !    variables (both dispersion ones and meteorological ones).
  ! 2. to order necessary input information in case of extra meteorological
  !    variables are requested by the user (input_needs mechanism)
  ! 3. to perform the output itself for all types of the output - GRIB, GrADS
  !    and trajectory. Corresponding io routines will be called from here,
  !    NOT from the particle models. That unit gets just fields and/or 
  !    trajectories as pointers - the rest has to be made here.
  ! 4. If there is no SILAM dispersion quantities in the output list - the 
  !    dispersion model itself must not be run - the signal should also 
  !    come from here
  !
  ! A TRICK. When user asks for e.g. nuclide_cocktail [source_inventory]
  !      it means that the list has to be created somewhere AFTER intialisation
  !      of sources and nuclides. So, this has to be called explicitly 
  !      somewhere in the particle_models. This server can do it by itself
  !      but particle_models is the only place, where correct moment is known
  !
  ! All units: SI
  !
  ! Code language: ANSI FORTRAN-90
  !
  ! Code owner Mikhail Sofiev, FMI  mikhail.sofiev@fmi.fi
  !
  use trajectory_io
  use advection_eulerian
  use advection_lagrangian
  use chemistry_manager  !dispersion_server
  use dispersion_server
  use optical_density
  use source_terms_general
  use pollution_cloud
  use silam_partitioning

  implicit none

  public set_OutDef_basic  ! Basic initialisation of output stack and supplementary vars
  public init_io_files     ! basic initialization of IO files: grads, grib, netcdf
  public global_io_init    ! nain IO initializtion
  public finalise_output   ! Finishes the output, closes files and writes a report
  public collect_output      ! pick the data
  public read_output_configuration ! Reads the file with list of requested variables
  public prepare_meteo
  public report

  ! Encapsulation
  public fu_ifTrajectory_in_output
  public fu_timestep
  private fu_timestep_of_output
  public fu_caseNm
  public fu_ifRunDispersion
  public fu_traj_set
  public get_output_quantities
  public fu_model_output_request
  public fu_species_output_name
  public fu_output_template
  public fu_output_dyn_shopping_list
  public fu_output_st_shopping_list
  public fu_ifDryDep_cumulative
  public fu_ifWetDep_cumulative
  public align_OutDef_initial_time_with_shift

  ! Private functions
  private cut_mesh
  private set_OV_missing
  private output_input_needs  ! Orders the quantities requested by the user for output
  private collect_buffers
  private collect_field_data      ! Takes field-based data to intermediate-IO stacks
  private write_output        ! write the data
  private tune_output_parameters   ! Decodes io settings and initialises intermediate fields
  private init_mass_output  ! Prepares intermediate mass maps
  private explore_item_species     ! Open-up species for the output list item
  private expand_dispersion_output_list  ! Decodes species for some SILAM variables
  private set_emission_output  ! Prepares and, if total requested, writes emission output
  private do_file_manip_grib
  private do_file_manip_grads
  private do_file_manip_netcdf
  private report_output_definition
  private shrink_output_list   ! Delete one quantity from the output list
  private start_new_output_period_lst ! Prepares targetIDs and nullifies averaged fields
  private fu_output_tmp_grid  ! Temporary output grid of the given quantity
  private report_output_list

  interface fu_timestep
    module procedure fu_timestep_of_output
  end interface

  interface report
    module procedure report_output_definition
  end interface
  !
  ! There can be three types of the substance list: SOURCE_INVENTORY, FULL_INVENTORY, and SINGLE_SUBSTANCE
  ! Inventories are given explicitly in the output_config, while the single substance is assumed if
  ! its name is stated
  ! CUSTOM_INVENTORY can be specified as well
  !
  integer, private, parameter :: iSourceInventory = 520
  integer, private, parameter :: iFullInventory = 521
  integer, private, parameter :: iTransportInventory = 522
  integer, private, parameter :: iAerosolInventory = 523
  integer, private, parameter :: iShortLivingInventory = 524
  integer, private, parameter :: iSingleSubstance = 525
  integer, private, parameter :: iNoSubstanceRelation = 526
  ! Not really an inventory  Substituted with  specise on reading  config 
  integer, private, parameter :: iCustomInventory = -1 
  !
  ! The main list of output variables. 
  ! Moved to netcdf_io
  !
!  type TOutputLstItem
!    private
!    integer :: quantity, request, AvType, iMode
!    real :: fWaveLen
!    logical :: if3D
!    type(silja_interval) :: AvPeriod
!    character(len=40) :: chSubstNm ! chemistry
!    character(len=fnlen) :: chSupplem
!    type(silja_field_id) :: targetId
!  end type TOutputLstItem

!  type TOutputList
!    private
!    type(TOutputLstItem), dimension(:), pointer :: ptrItem
!  end type TOutputList

  !
  ! Full description of the output rules and parameters
  !
  TYPE TOutputRules
!    private
    type(TOutputList) :: MeteoOutLst, DispOutLst, MassMapOutLst  ! Lists of output quantities
    type(silja_time) :: ini_time             ! Start of simulations
    TYPE(silja_interval) :: timestep, dump_timestep  ! output timesteps - main and for technical dump
    TYPE(silja_interval) :: rates_dump_timestep !Dump step for reaction rates
    integer :: iOutTimesType                 ! REGULAR or special occasions
!    TYPE(silja_grid) :: grid                 ! horizontal grid of output
!    TYPE(silam_vertical) :: vertical         ! vertical of the output
    logical :: ifRunDispersion               ! May be, just process meteodata...
    logical :: ifAverageMeteo                ! If any of meteo is not AS_IS
    type(silja_logical) :: ifDryDepCumulative, ifWetDepCumulative  ! rates or cumulatives?
    logical :: ifSplitSizeModes              ! If aerosol size modes are summed-up or separate
    LOGICAL :: ifGRADS, ifTrajectory  ! Simple output format switches
    integer :: particles_per_trajectory ! Dump every particles_per_trajectory particle
    integer :: iGrib, iNETCDF         ! 1 or 2 for GRIB, 3 or 4  for NETCDF, int_missing for void
    integer :: OutFilesArrangement           ! The way time is split to the output files
    integer :: cldRepInterv                  ! Reporting interval of cloud mass
    real    :: MMtrimfactor                 ! If positive -- how many discretes per factor 2 to leave in output massmaps
                                            ! Enables easy compression
    type(grads_template) :: outTemplate      ! output file name template without extentions
    character(len=fnlen) :: chFixedNameTemplate
    logical :: ifMeteo2OutHorizInterp, ifMeteo2OutVertInterp, &
             & ifDisp2OutHorizInterp, ifDisp2OutVertInterp
  END TYPE TOutputRules

  TYPE (TOutputRules) , public, parameter :: OutputRules_missing = &
        & TOutputRules(OutputList_missing, OutputList_missing, OutputList_missing,&
              & time_missing, interval_missing,  interval_missing, interval_missing, int_missing,&
              & .False.,.False.,silja_undefined, silja_undefined, &
              & .False.,.False.,.False., int_missing, &
              & int_missing, int_missing, int_missing, int_missing,&
              & real_missing, grads_template_missing, '', &
              & .False.,.False.,.False.,.False.)


  !
  ! Output parameters should also be stored: current units of all files,
  ! their names and the structures describing the GrADS and GRIB files. 
  ! Dimension is set because there can be several outputs processed simulataneously.
  ! The idea is the following. In lagrange framework it is easy 
  ! to attribute each particle to its source. So, it is easy to process several
  ! sources simultaneously writing their results in separate output files. 
  ! Saving is evident - meteodata are processed only once.
  !
  type TOutputParams
!    private
    character(len=clen) :: chCaseNm
    character(len=clen), dimension(:), pointer :: chSrcNm
    integer, dimension(:), pointer :: grib_funit, grads_funit, netcdf_funit
    character(len=fnlen), dimension(:), pointer :: grib_fNm, grads_fNm, &
                                                 & netcdf_fNm, filesWritten
    type(silam_trajectory_set), pointer :: tr_set
    type(grads_file), dimension(:), pointer :: grib_struct, grads_struct, netcdf_struct
  end type TOutputParams
  type(TOutputParams), public, parameter :: OutputParams_missing = &
        &TOutputParams("",null(),null(),null(),null(), &
       & null(),null(),null(),null(),null(),null(),&
       & null(),null())

  !
  ! For output we need several temporary variables, as well as several 
  ! permanent data. Let's define a unique place for them.
  ! Specific variable is SubstMap, which keeps all information where to search
  ! the data for the substance - particle numbers where it is available as
  ! well as the index in the material lists for each particle.
  ! We also will need two meteo and two dispersion output stacks. Tmp ones will be
  ! use for intermediate accumulation of the fields in the system grid (without
  ! interpolation), while the main ones are to be in the output grid. Evidently, 
  ! if system_grid = output_grid the Tmp pointers and the main ones will point 
  ! to the same stacks. Otherwise, there will be interpolation coefficients structures
  ! filled with the reasonable stuff.
  !
  type TOutputVariables
    private
    TYPE(silja_stack), pointer :: meteoStack, dispStack ! Output stacks
    TYPE(silja_stack), pointer :: meteoTmpStack              ! Intermediate output stacks
    TYPE(silja_stack), dimension(:), pointer :: dispTmpStack ! Intermediate output stacks
    type(TMassMapLink) :: MassMapLinks
    type(TLagrangeToMassMapLink) :: Lagr2MMLinks
!    real, dimension(:), pointer :: dx, dy, dz ! Sizes of the output grid cells - permanent data
!    type(field_3d_data_ptr) :: field3d, cnc_air   ! nn 1d array pointers and 1d array of field ids
    type(silja_time) :: lastOutputTime, LastDumpOutputTime, LastRatesDumpOutputTime
    integer :: iNbrCollection ! How many times variables were collected to output stacks
    type(THorizInterpStruct), pointer :: interpCoefMeteo2OutHoriz, & ! meteofields to output ones
                                       & interpCoefDisp2OutHoriz  ! dispersion fields o output ones
    type(TVertInterpStruct), pointer :: interpCoefMeteo2OutVert, interpCoefDisp2OutVert
  end type TOutputVariables
  private TOutputVariables
  !
  ! Output definition uses the above structures:
  !
  type silam_output_definition
!    private
    type(TOutputParams) :: Params 
    type(TOutputRules) :: Rules
!    type(TOutputVariables) :: Vars
  end type silam_output_definition

  type(silam_output_definition) , public, parameter :: output_definition_missing = &
                & silam_output_definition(OutputParams_missing,OutputRules_missing)

  type(TOutputVariables), private, target, save :: OutVars

CONTAINS


  !*****************************************************************************

  subroutine output_input_needs(OutDef, q_st, req_st, q_dyn,  req_dyn, qd_dyn, qd_st, wdr)
    !
    ! Does two things:
    ! 1. Orders necessary quantities for interpolation from system grid/vertical
    !    to the output ones.
    ! 2. Orders the quantities requested by the user for the explicit output
    ! 
    !  Huom! No massmap quantities added here
    !
    ! Typical usage is - if the user wants to have, e.g., ABL height in metres in 
    ! the output, this very function puts it into the shopping_list, which will
    ! later result in creation of this very field, whether it is needed for the 
    ! simulations or not.
    ! Another example - if the simulations are not required at all and just meteo
    ! data are needed. This very function will take care of the stuff.
    !
    implicit none

    ! Imported parameters
    integer, dimension(:), intent(inout) :: q_st, req_st, q_dyn,  req_dyn, qd_dyn, qd_st
    type(silam_output_definition), intent(in) :: OutDef
    type(silja_wdr), pointer :: wdr

    ! Local variables
    integer :: i, iQs, iQd, Qtmp 

    call msg("Making output input needs")

    q_dyn(1) = int_missing
    q_st(1) = int_missing
    req_dyn(1) = int_missing
    req_st(1) = int_missing
    !
    ! Trajectory output may require specific q_dyn 
    !
    if(OutDef%Rules%ifTrajectory)then
      CALL trajectory_input_needs(q_dyn, iQd, q_st, iQs)  !---- Link to trajectory_io module
      req_dyn(1:iQd) = 2
      q_dyn(iQd+1) = int_missing
      req_st(1:iQs) = 2
      q_st(iQs+1) = int_missing
    else
      iQd = 0
      iQs = 0
    endif

    !
    ! q_dyn for recalculation from the meteo vertical to the output one
    ! depend on the type of the output levels. System parameters are not yet 
    ! known. They will be taken into account later by tune_output_parameters
    !
    if(fu_NbrOfLevels(output_vertical) > 1)then
      select case(fu_leveltype(output_vertical))
        case(constant_pressure, layer_btw_2_pressure)
          iQd = fu_merge_integer_to_array(pressure_flag, q_dyn, 2, req_dyn)
        case(constant_height, layer_btw_2_height)
          iQd = fu_merge_integer_to_array(height_flag, q_dyn, 2, req_dyn)
        case(constant_altitude, layer_btw_2_altitude)
          iQd = fu_merge_integer_to_array(height_flag, q_dyn, 2, req_dyn)
        case(hybrid, layer_btw_2_hybrid)
          iQd = fu_merge_integer_to_array(pressure_flag, q_dyn, 2, req_dyn)
          q_dyn(iQd+1) = int_missing
          iQd = fu_merge_integer_to_array(ground_pressure_flag, q_dyn, 2, req_dyn)
        case default ! just 2D output, no vertical interpolation 
          call msg_warning('Possibly strange output vertical','output_input_needs')
          call report(output_vertical)
          return
      end select
    endif
    q_dyn(iQd+1) = int_missing
    !
    ! Add a relief height - this is needed for various conversions, including often
    ! connection between the dispersion and meteorological and dispersion verticals
    !
    iQs = fu_merge_integer_to_array(relief_height_flag, q_st, 2, req_st)
    q_st(iQs+1) = int_missing

    !
    !Surface pressure is needed for many things
    iQd = fu_merge_integer_to_array(ground_pressure_flag, q_dyn, 2, req_dyn)
    q_dyn(iQd+1) = int_missing
    !
    !
    !

    if ( allocated(OutDef%Rules%MeteoOutLst%ptrItem)) then
      do i=1, size(OutDef%Rules%MeteoOutLst%ptrItem)
        Qtmp = OutDef%Rules%MeteoOutLst%ptrItem(i)%quantity

        if(.not.fu_known_quantity(Qtmp))cycle
      
        ! Dirty hack for LAI
        if ( any(Qtmp == (/leaf_area_index_flag, leaf_area_indexhv_flag, leaf_area_indexlv_flag/))) then
          if ( any(fu_LAIsrc(wdr) == (/LAI_static_1, LAI_static_2/) )) then
           iQs = fu_merge_integer_to_array(Qtmp, q_st, &
                        & OutDef%Rules%MeteoOutLst%ptrItem(i)%request, req_st)
           q_st(iQs+1) = int_missing
          else
           iQd = fu_merge_integer_to_array(Qtmp, q_dyn, &
                        & OutDef%Rules%MeteoOutLst%ptrItem(i)%request, req_dyn)
           q_dyn(iQd+1) = int_missing
          endif
          cycle
        endif

        if (fu_realtime_quantity(Qtmp)) then 
          iQs = fu_merge_integer_to_array(Qtmp, q_st, &
                                         & OutDef%Rules%MeteoOutLst%ptrItem(i)%request, req_st)
          q_st(iQs+1) = int_missing
        else
          iQd = fu_merge_integer_to_array(Qtmp, q_dyn, &
                                         & OutDef%Rules%MeteoOutLst%ptrItem(i)%request, req_dyn)
          q_dyn(iQd+1) = int_missing
        endif
      end do
    endif
    !
    ! A stupid workaround: we may need pressure no matter what, if someone will need
    ! z-size, so let's add it here
    !
    iQd = fu_merge_integer_to_array(pressure_flag, q_dyn, 2, req_dyn)



    ! Add requested _Buffer_ quantities
    ! If something does not get to the requests it will not be available from buffer
    !
    if ( allocated(OutDef%Rules%DispOutLst%ptrItem)) then
      do i=1, size(OutDef%Rules%DispOutLst%ptrItem)
        Qtmp = OutDef%Rules%DispOutLst%ptrItem(i)%quantity

        if(.not.fu_known_quantity(Qtmp))cycle
        if ( fu_if_cloud_mass_map_quantity(Qtmp) == silja_true ) cycle 
              !! No need to put it in a buffer: it is a massmap quantity

        ! Not actually needed... I was lazy to push the requests through
        if(OutDef%Rules%DispOutLst%ptrItem(i)%request /= 2) cycle 
     
        ! No requests for dispesion quantities, sorry 
        if (fu_realtime_quantity(Qtmp)) then 
          iQs = fu_merge_integer_to_array(Qtmp, qd_st)
          qd_st(iQs+1) = int_missing
        else
          iQd = fu_merge_integer_to_array(Qtmp, qd_dyn)
          qd_dyn(iQd+1) = int_missing
        endif
      end do
    endif

  end subroutine output_input_needs

  
  !*****************************************************************************
  
  subroutine set_OV_missing(OV)
    type(TOutputVariables), intent(inout) :: OV

    nullify(OV%meteoStack, OV%dispStack, OV%meteoTmpStack, OV%dispTmpStack)
    call set_MML_missing(OV%MassMapLinks)
    call set_LMML_missing(OV%Lagr2MMLinks)
    OV%lastOutputTime = time_missing
    OV%LastDumpOutputTime = time_missing
    OV%LastRatesDumpOutputTime = time_missing
    OV%iNbrCollection = int_missing
    nullify(OV%interpCoefMeteo2OutHoriz, OV%interpCoefDisp2OutHoriz, &
          & OV%interpCoefMeteo2OutVert, OV%interpCoefDisp2OutVert)
  end subroutine set_OV_missing
  

  !*****************************************************************************

  subroutine set_OutDef_basic(nlSetup, nlStandardSetup, OutDef, &
                                          & model_time_step, simulation_type, chCaseNm)
    !
    ! Takes care of initialisation of the output structures - outParams
    ! and Rules. They are consumed from the namelist
    !
    implicit none

    ! return value
    type(silam_output_definition),intent(inout), target :: OutDef

    ! Imported parameters 
    type(Tsilam_namelist), pointer :: nlSetup, nlStandardSetup
    type(silja_interval), intent(in) :: model_time_step
    character(len=*), intent(in) :: chCaseNm
    integer, intent(in) :: simulation_type

    ! Local variables
    integer :: i, j, k, al_status, nx,ny, iSrcId, time_sign, iAvType,iTmp
    logical :: ifFound, ifPrssureNeeded
    real :: dx, dy
    character(len=fnlen) :: strTmp
    type(silja_level), dimension(:), pointer :: levels
    type(silam_species) :: speciesTmp
    type(Tsilam_nl_item_ptr), dimension(:), pointer :: ptrVars

    call msg("Making set_OutDef_basic")
    if(.not.defined(nlSetup))then
      call set_error('Undefined namelist given','set_OutDef_basic')
      return
    endif

    if(model_time_step > zero_interval)then
      time_sign = 1
    else
      time_sign = -1
    endif

    !----------------------------------------------------------------------------
    !
    ! Rules. Store the variables with necessary checking for stupidity.
    ! Note that ifRunDispersion is handled elsewhere
    !
    ! GRID. Undefined grid is possible - it will be then set from system_grid
    ! or somehow else.
    ! Grids available: emission, meteorology, output and internal
    !
    ! Emission grid: 
    ! decided by sources, whcih can be cut down to the dispersion grid. Must be inside
    ! the meteorological and dispersion grids
    !
    ! Meteo grid:
    ! max span is decided by available data, actual size is cut down to smallest
    ! area covering entirely the dispersion grid
    ! 
    ! Output grid:
    ! Defined by the control file. Can be CUSTOM, which is the non-negotiable definition
    ! Can be EMISSION defined by the sources reduced down to meteorological grid, if needed
    ! Can be AREA_BASED constructed around the given area cut down to meteo grid, if needed
    ! Can be METEO copied from meteo grid, which is then built around the sources
    !
    ! Dispersion grid
    ! Must cover the sources and output and be covered by meteorology.
    ! Can be CUSTOM, which is non-negotiable definition.
    ! Can be EMISSION defined by the sources reduced down to meteorological grid, if needed
    ! Can be METEO copied from meteo grid
    ! Can be OUTPUT copied from output grid
    !
    ! At this stage we create only the output grid if possible - i.e. if it is AREA_BASED
    ! or CUSTOM
    !
    if(fu_str_u_case(fu_content(nlSetup,'grid_method')) == 'AREA_BASED')then
      !
      ! AREA_BASED output grid
      !
      strTmp = fu_str_u_case(fu_content(nlSetup,'resolution'))
      read(unit=strTmp,fmt=*,iostat=al_status)dx
      if(al_status /= 0)then
        call set_error('Failed to read grid resolution from:' + strTmp, &
                     & 'set_OutDef_basic')
        return
      endif
      if(index(strTmp,' M') > 0)then
        dy=dx
      elseif(index(strTmp,' KM') > 0)then
        dx=dx*1000.
        dy=dx
      elseif(index(strTmp,' DEG') > 0)then
        dy=dx
      else
        call set_error('Strange resolution:' + strTmp,'set_OutDef_basic')
      endif
      output_grid = fu_lonlat_grid_from_area(fu_set_area(nlSetup), &
                                           & (index(strTmp,'M') > 0), dx, dy)

    elseif(fu_str_u_case(fu_content(nlSetup,'grid_method')) == 'CUSTOM_GRID')then
      !
      ! CUSTOM output grid: non-negotiable
      !
        output_grid = fu_set_grid(nlSetup)
    elseif(fu_str_u_case(fu_content(nlSetup,'grid_method')) == 'METEO_GRID')then
      !
      ! Output grid is a derivative of emission or meteorology
      !
      output_grid = grid_missing  ! Cannot determine it now. Leave for later
    else
      call set_error('Unknown output grid type:' + fu_content(nlSetup,'grid_method'), &
                   & 'set_OutDef_basic')
      return
    endif
    if(error)return

    !
    ! LEVELS. Undefined levels are possible - they will then be copied from
    !         the meteo_vertical
    !
    if(fu_str_u_case(fu_content(nlSetup,'vertical_method')) == 'METEO_LEVELS')then
      call set_missing(output_vertical, .true.)
    else
      call set_vertical(nlSetup, output_vertical)
    endif
    if(error)return

    !
    ! TIME.  Undefined time is possible - then the model output will happen
    !        every meteo time step
    !        It is also possible to ask for output at "special occasions":
    !        REGULAR / START_OF_DAY / START_OF_WEEK / START_OF_MONTH / START_OF_YEAR /
    !                  END_OF_DAY / END_OF_WEEK / END_OF_MONTH / END_OF_YEAR
    !        In case of REGULAR, the output_time_step comes into play, otherwise
    !        the output timestep becomes dynamic.
    !
    ! There can be some dump requested
    !
    call msg('Dump Interval:' + fu_content(nlSetup,'dump_time_step'), &
           & len_trim(fu_content(nlSetup,'dump_time_step')))
    if(len_trim(fu_content(nlSetup,'dump_time_step')) > 0)then
      OutDef%Rules%dump_timestep = fu_set_named_interval(fu_content(nlSetup,'dump_time_step')) * &
                                 & real(time_sign)
    else
      OutDef%Rules%dump_timestep = interval_missing
    endif

    call msg('Reaction rate dump Interval:' + fu_content(nlSetup,'rates_dump_step'), &
           & len_trim(fu_content(nlSetup,'rates_dump_step')))
    if(len_trim(fu_content(nlSetup,'rates_dump_step')) > 0)then
      OutDef%Rules%rates_dump_timestep = fu_set_named_interval(fu_content(nlSetup,'rates_dump_step')) * &
                                 & real(time_sign)
    else
      OutDef%Rules%rates_dump_timestep = interval_missing
    endif
    !
    ! The output for users can be both regular and tight to some specific moments
    !
    strTmp = fu_str_u_case(fu_content(nlSetup,'output_times'))
    if(index(strTmp,'REGULAR') == 1)then
      OutDef%Rules%timestep = fu_set_named_interval(fu_content(nlSetup,'output_time_step')) &
                            & * real(time_sign)
      if(OutDef%Rules%timestep == zero_interval .or. &
       & OutDef%Rules%timestep == interval_missing .or. &
       & error)then
         call msg_warning('Cannot determine the output time step, assume meteo step', &
                        & 'set_OutDef_basic')
         if(error) then 
                 call set_error("Error", 'set_output_def_from_namelist')
                 return
         endif
         OutDef%Rules%timestep = interval_missing
         OutDef%Rules%iOutTimesType = iMeteoTime
      else
         OutDef%Rules%iOutTimesType = iRegular
      end if
      !
      ! A small adjustment: the output time should coinside with the model time step
      !
      OutDef%Rules%timestep = model_time_step * nint(OutDef%Rules%timestep / model_time_step)
      if(defined(OutDef%Rules%dump_timestep))then
        OutDef%Rules%dump_timestep = model_time_step * &
                                   & nint(OutDef%Rules%dump_timestep / model_time_step)
      endif
      if(defined(OutDef%Rules%rates_dump_timestep))then
        OutDef%Rules%rates_dump_timestep = model_time_step * &
                                   & nint(OutDef%Rules%rates_dump_timestep / model_time_step)
      endif
    else
      OutDef%Rules%timestep = interval_missing
      if(index(strTmp,'START_OF_DAY') == 1)then
        OutDef%Rules%iOutTimesType = iStartOfDay
      elseif(index(strTmp,'START_OF_WEEK') == 1)then
        OutDef%Rules%iOutTimesType = iStartOfWeek
      elseif(index(strTmp,'START_OF_MONTH') == 1)then
        OutDef%Rules%iOutTimesType = iStartOfMonth
      elseif(index(strTmp,'START_OF_YEAR') == 1)then
        OutDef%Rules%iOutTimesType = iStartOfYear
      elseif(index(strTmp,'END_OF_DAY') == 1)then
        OutDef%Rules%iOutTimesType = iEndOfDay
      elseif(index(strTmp,'END_OF_WEEK') == 1)then
        OutDef%Rules%iOutTimesType = iEndOfWeek
      elseif(index(strTmp,'END_OF_MONTH') == 1)then
        OutDef%Rules%iOutTimesType = iEndOfMonth
      elseif(index(strTmp,'END_OF_YEAR') == 1)then
        OutDef%Rules%iOutTimesType = iEndOfYear
      else
         call set_error('Cannot determine the output time step, assume meteo step', &
                        & 'set_OutDef_basic')
         call unset_error('set_OutDef_basic')
         call msg_warning('Cannot determine the output time step, assume meteo step', &
                        & 'set_OutDef_basic')
         OutDef%Rules%timestep = interval_missing
         OutDef%Rules%iOutTimesType = iMeteoTime
      endif
    endif

    !
    ! Output formats
    !
    nullify(ptrVars)
    OutDef%Rules%iGrib = int_missing
    OutDef%Rules%iNETCDF = int_missing
    OutDef%Rules%ifGrads = .false.
    OutDef%Rules%ifTrajectory = .false.
    OutDef%Rules%particles_per_trajectory = 500 !Uesd to be hardcoded in trajectory outptut
    
    call get_items(nlSetup, 'output_format', ptrVars, i)
    if(i < 1)then
      strTmp = fu_str_u_case(fu_content(nlSetup,'file_types'))
      if(fu_fails(len_trim(strTmp) > 5, &
                & 'Neither file_types nor output_format lines found in output namelsit', &
                & 'set_OutDef_basic'))return
      if(index(strTmp,'GRIB_YES') > 0)OutDef%Rules%iGrib = 1
      if(index(strTmp,'GRIB1_YES') > 0)OutDef%Rules%iGrib = 1
      if(index(strTmp,'GRIB2_YES') > 0)OutDef%Rules%iGrib = 2
      if(index(strTmp,'NETCDF_YES') > 0)OutDef%Rules%iNetCDF = 3
      if(index(strTmp,'NETCDF3_YES') > 0)OutDef%Rules%iNetCDF = 3
      if(index(strTmp,'NETCDF4_YES') > 0)OutDef%Rules%iNetCDF = 4
      OutDef%Rules%ifGrads = (index(strTmp,'GRADS_YES') > 0)
      OutDef%Rules%ifTrajectory = (index(strTmp,'TRAJECTORY_YES') > 0)
    else
      do j = 1, i
        if(index(fu_str_u_case(fu_content(ptrVars(j))),'GRIB') > 0)then        ! GRIB
          if(index(fu_content(ptrVars(j)),'1') > 0)then
            OutDef%Rules%iGrib = 1
          elseif(index(fu_content(ptrVars(j)),'2') > 0)then
            OutDef%Rules%iGrib = 2
          else
            call set_error('GRIB output requested without version (1 or 2)', &
                         & 'set_OutDef_basic')
          endif
        endif
        if(index(fu_str_u_case(fu_content(ptrVars(j))),'NETCDF') > 0)then        ! NetCDF
          if(index(fu_content(ptrVars(j)),'3') > 0)then
            OutDef%Rules%iNetCDF = 3
          elseif(index(fu_content(ptrVars(j)),'4') > 0)then
            OutDef%Rules%iNetCDF = 4
          else
            call set_error('NETCDF output requested without version (3 or 4)', &
                       & 'set_OutDef_basic')
          endif
        endif
        if(index(fu_str_u_case(fu_content(ptrVars(j))),'GRADS') > 0) &           ! GRADS
                                                 & OutDef%Rules%ifGrads = .true.
        if(index(fu_str_u_case(fu_content(ptrVars(j))),'TRAJECTORY') > 0) &      ! TRAJECTORY
                                                   & OutDef%Rules%ifTrajectory = .true.
      end do
    endif

    iTmp = fu_content_int(nlSetup,"particles_per_trajectory")
    if (iTmp > 0)     OutDef%Rules%particles_per_trajectory = iTmp

    if (OutDef%Rules%ifTrajectory) call msg("Saving trajectory for every "// &
        & trim(fu_str(OutDef%Rules%particles_per_trajectory))//"-th particle") 

    if(OutDef%Rules%iGrib == int_missing .and. OutDef%Rules%iNETCDF == int_missing .and. &
                         & .not. OutDef%Rules%ifGrads .and. .not. OutDef%Rules%ifTrajectory)then
      call set_error('No output format requested','set_OutDef_basic')
      return
    else
      call msg('OutDef%Rules%iGrib, OutDef%Rules%iNetCDF', OutDef%Rules%iGrib, OutDef%Rules%iNetCDF)
    endif
    
    !
    ! Cut precision -- superior for compression
    !
    OutDef%Rules%MMtrimfactor = fu_content_real(nlSetup,'massmap_precision_factor')
    if(OutDef%Rules%MMtrimfactor /= real_missing)then
      call msg("Accuracy cut: massmap_precision_factor=", OutDef%Rules%MMtrimfactor)
      if(OutDef%Rules%MMtrimfactor < 2. .or. OutDef%Rules%MMtrimfactor > 1e6)then
        call set_error('Precision factor must be between 2 and 1e6, or absent', &
                     & 'set_OutDef_basic')
        return
      endif
    endif

    !
    ! Output files arrangement
    !
    strTmp = fu_str_u_case(fu_content(nlSetup,'time_split'))
    if(strTmp== 'ALL_IN_ONE')then
      OutDef%Rules%OutFilesArrangement = all_in_one
    elseif(strTmp == 'HOURLY_NEW_FILE')then
      OutDef%Rules%OutFilesArrangement = hourly_new_file
    elseif(strTmp == 'DAILY_NEW_FILE')then
      OutDef%Rules%OutFilesArrangement = daily_new_file
    elseif(strTmp == 'MONTHLY_NEW_FILE')then
      OutDef%Rules%OutFilesArrangement = monthly_new_file
    elseif(strTmp == 'YEARLY_NEW_FILE')then
      OutDef%Rules%OutFilesArrangement = yearly_new_file
    else
      call set_error('Unknown output file arrangement:' + strTmp, 'set_OutDef_basic')
      return
    end if

    !
    ! GrADS output can require a duplicate ctl file to be created - e.g. "latest.ctl"
    ! Can be empty
    !
    if(OutDef%Rules%ifGrads)then
      outDef%Rules%chFixedNameTemplate = &
                               & fu_content(nlSetup,'grads_metafiles_duplicate_fixed_name_template')
    else
      outDef%Rules%chFixedNameTemplate = ''
    endif
    !
    ! How to handle multi-source case
    !
    if(index(fu_str_u_case(fu_content(nlSetup,'source_id')),'NO_SOURCE_SPLIT') > 0)then
      iSrcId = iNoId  ! Mix sources
    elseif(index(fu_str_u_case(fu_content(nlSetup,'source_id')),'SOURCE_NAME') > 0)then
      iSrcId = iSrcNmId
    elseif(index(fu_str_u_case(fu_content(nlSetup,'source_id')),'SOURCE_SECTOR')>0)then
      iSrcId = iSrcSectorId
    elseif(index(fu_str_u_case(fu_content(nlSetup,'source_id')), 'SOURCE_NAME_AND_SECTOR') > 0)then
      iSrcId = iSrcNmSectorId
    elseif(index(fu_str_u_case(fu_content(nlSetup,'source_id')), 'TIME_ZONE') > 0)then
      iSrcId = iSrcTimeZoneId
    else
      call set_error('Strange source id string:' + fu_content(nlSetup,'source_id'), &
                   & 'set_OutDef_basic')
      return
    endif

    !
    ! Output file template. It has to be checked for consistency with other
    ! rules affecting the number and the content of the output files - 
    ! ifSplitSources and FileArrangement. In particular, we have to confirm
    ! that every file will have unique name
    !
    strtmp = fu_content(nlSetup,'template')
    call check_dir_slash(strtmp)
    if(error)return
    if(fu_check_outTemplate(strTmp, OutDef%Rules%OutFilesArrangement, iSrcId))then
      call decode_template_string(strTmp, OutDef%Rules%outTemplate)
    else
      call set_error('Output template is inconsistent with other rules:' + strTmp, &
                   & 'set_OutDef_basic')
    endif
    if(error)return

    !
    ! List of output variables
    ! There are a few standard output lists - they do not imply any reading
    ! as well as do not imply any meteo output. It is assumed that the output
    ! always takes full information including the aerosol size split
    !
    !Init OutLsts
    call expand_output_list(OutDef%Rules%MeteoOutLst, 1)
    call expand_output_list(OutDef%Rules%DispOutLst, 1)
    call expand_output_list(OutDef%Rules%MassMapOutLst, 1)


    OutDef%Rules%ifRunDispersion = .true.

    strTmp = fu_content(nlSetup,'variable_list')

    if(index(fu_str_u_case(strTmp),'PASSIVE_DISPERSION_LIST') > 0)then

      call set_species(speciesTmp, fu_get_material_ptr('passive'), aerosol_mode_missing)

      if(fu_if_lagrangian_present(simulation_type))then
        call expand_output_list(OutDef%Rules%DispOutLst, 3)
        OutDef%Rules%DispOutLst%ptrItem(1:3)%species = speciesTmp
        OutDef%Rules%DispOutLst%ptrItem(1)%quantity = particle_counter_flag
        OutDef%Rules%DispOutLst%ptrItem(1)%request = 2
        OutDef%Rules%DispOutLst%ptrItem(1)%AvType = iAsIs                ! Instant, in fact
        OutDef%Rules%DispOutLst%ptrItem(2)%quantity = areas_of_risk_flag
        OutDef%Rules%DispOutLst%ptrItem(2)%request = 2
        OutDef%Rules%DispOutLst%ptrItem(2)%AvType = iAsIs                ! Instant, in fact
        OutDef%Rules%DispOutLst%ptrItem(3)%quantity = concentration_flag
        OutDef%Rules%DispOutLst%ptrItem(3)%request = 2
        OutDef%Rules%DispOutLst%ptrItem(3)%AvType = iAsIs                ! Instant, in fact
        OutDef%Rules%ifSplitSizeModes = .false.
      elseif(fu_if_eulerian_present(simulation_type))then
        call expand_output_list(OutDef%Rules%DispOutLst, 1)
        OutDef%Rules%DispOutLst%ptrItem(1)%species = speciesTmp
        OutDef%Rules%DispOutLst%ptrItem(1)%quantity = concentration_flag
        OutDef%Rules%DispOutLst%ptrItem(1)%request = 2
        OutDef%Rules%DispOutLst%ptrItem(1)%AvType = iAsIs                ! Instant, in fact
        OutDef%Rules%ifSplitSizeModes = .false.
      else
        call set_error('Unknown advection method:' + fu_str(simulation_type),'set_OutDef_basic')
        return
      end if

    elseif(index(fu_str_u_case(strTmp),'SOURCE_INVENTORY_DISPERSION_LIST') > 0)then
      call expand_output_list(OutDef%Rules%DispOutLst, 4)
      OutDef%Rules%DispOutLst%ptrItem(1)%quantity = concentration_flag
      OutDef%Rules%DispOutLst%ptrItem(1)%request = 2
      OutDef%Rules%DispOutLst%ptrItem(1)%AvType = iAverage
      OutDef%Rules%DispOutLst%ptrItem(1)%iSpeciesListType = iSourceInventory   !'SOURCE_INVENTORY'
      OutDef%Rules%DispOutLst%ptrItem(2)%quantity = concentration_2m_flag
      OutDef%Rules%DispOutLst%ptrItem(2)%request = 2
      OutDef%Rules%DispOutLst%ptrItem(2)%AvType = iAverage
      OutDef%Rules%DispOutLst%ptrItem(2)%iSpeciesListType = iSourceInventory   !'SOURCE_INVENTORY'
      OutDef%Rules%DispOutLst%ptrItem(3)%quantity = drydep_flag
      OutDef%Rules%DispOutLst%ptrItem(3)%request = 2
      OutDef%Rules%DispOutLst%ptrItem(3)%AvType = iCumulative
      OutDef%Rules%DispOutLst%ptrItem(3)%iSpeciesListType = iSourceInventory   !'SOURCE_INVENTORY'
      OutDef%Rules%DispOutLst%ptrItem(4)%quantity = wetdep_flag
      OutDef%Rules%DispOutLst%ptrItem(4)%request = 2
      OutDef%Rules%DispOutLst%ptrItem(4)%AvType = iCumulative
      OutDef%Rules%DispOutLst%ptrItem(4)%iSpeciesListType = iSourceInventory   !'SOURCE_INVENTORY'

      OutDef%Rules%ifSplitSizeModes = .true.

    elseif(index(fu_str_u_case(strTmp),'FULL_INVENTORY_DISPERSION_LIST') > 0)then
      call expand_output_list(OutDef%Rules%DispOutLst, 4)
      OutDef%Rules%DispOutLst%ptrItem(1)%quantity = concentration_flag
      OutDef%Rules%DispOutLst%ptrItem(1)%request = 2
      OutDef%Rules%DispOutLst%ptrItem(1)%AvType = iAverage
      OutDef%Rules%DispOutLst%ptrItem(1)%iSpeciesListType = iFullInventory   !'FULL_INVENTORY'
      OutDef%Rules%DispOutLst%ptrItem(2)%quantity = concentration_2m_flag
      OutDef%Rules%DispOutLst%ptrItem(2)%request = 2
      OutDef%Rules%DispOutLst%ptrItem(2)%AvType = iAverage
      OutDef%Rules%DispOutLst%ptrItem(2)%iSpeciesListType = iFullInventory   !'FULL_INVENTORY'
      OutDef%Rules%DispOutLst%ptrItem(3)%quantity = drydep_flag
      OutDef%Rules%DispOutLst%ptrItem(3)%request = 2
      OutDef%Rules%DispOutLst%ptrItem(3)%AvType = iCumulative
      OutDef%Rules%DispOutLst%ptrItem(3)%iSpeciesListType = iFullInventory   !'FULL_INVENTORY'
      OutDef%Rules%DispOutLst%ptrItem(4)%quantity = wetdep_flag
      OutDef%Rules%DispOutLst%ptrItem(4)%request = 2
      OutDef%Rules%DispOutLst%ptrItem(4)%AvType = iCumulative
      OutDef%Rules%DispOutLst%ptrItem(4)%iSpeciesListType = iFullInventory   !'FULL_INVENTORY'


      OutDef%Rules%ifSplitSizeModes = .true.

    else
      call read_output_configuration(fu_expand_environment(strTmp), OutDef, model_time_step)
    endif
    if(error)return
    
    !
    ! Request pressure to hybrid-level output
    if (any (fu_leveltype(output_vertical) == &
         & (/layer_btw_2_hybrid, layer_btw_2_sigma, hybrid, sigma_level/))) then
         ! Must have surface pressure in the output if any 3D is there
         ! If any 3d is not instant -- pressure has to be averaged 
         ifPrssureNeeded=.False.
         iavtype = iAsIs
         ! Should check for associated before doing "size"
         if (allocated(OutDef%Rules%MeteoOutLst%ptrItem)) then
           do i = 1, size(OutDef%Rules%MeteoOutLst%ptrItem)
              if (OutDef%Rules%MeteoOutLst%ptrItem(i)%if3D) then 
                  ifPrssureNeeded=.true.
                  if (OutDef%Rules%MeteoOutLst%ptrItem(i)%AvType /= iAsIs) iAvType = iAverage
              endif
              k = OutDef%Rules%MeteoOutLst%ptrItem(i)%quantity
              ifFound = (k == ground_pressure_flag) 
              if (ifFound .or. k == int_missing) exit
           end do
         endif
         if (ifFound) then 
            call msg("Surface pressure already requested")
         elseif (.not. ifPrssureNeeded) then
            call msg("Surface pressure not needed in the output")
         else
            call msg("Forcing surface pressure output")
            if(i > size(OutDef%Rules%MeteoOutLst%ptrItem)) &
                          & call expand_output_list(OutDef%Rules%MeteoOutLst, 2)

            OutDef%Rules%MeteoOutLst%ptrItem(i)%quantity = ground_pressure_flag
            OutDef%Rules%MeteoOutLst%ptrItem(i)%request = 2
            OutDef%Rules%MeteoOutLst%ptrItem(i)%AvType = iAvType
            OutDef%Rules%MeteoOutLst%ptrItem(i)%if3D = .False.
            OutDef%Rules%MeteoOutLst%ptrItem(i)%iVerticalTreatment = do_nothing_flag
         endif
   endif

    OutDef%Params%chCaseNm = chCaseNm

    !
    ! Frequency of cloud reporting - in units of output interval
    !
    i = fu_content_int(nlStandardSetup,'cloud_report_interval')
    if(i > 0 .and. i < 1000)then
      OutDef%Rules%cldRepInterv = i
    else
      OutDef%Rules%cldRepInterv = 1
    endif

  end subroutine set_OutDef_basic
    

  !*****************************************************************************

  subroutine init_io_files(nlStdSetup, nlIOSetup, max_nbr_of_grads_files)
    !
    ! Basic initialization of the input and output file structures
    !
    implicit none
    
    ! Imported parameters
    type(Tsilam_namelist), pointer :: nlStdSetup, nlIOsetup
    integer, intent(in) :: max_nbr_of_grads_files
    
    call init_grib_io(nlStdSetup)
    
    call init_grads_io(max_nbr_of_grads_files) ! Number of files allowed to be opened simultaneously

    call init_netcdf_io(nlStdSetup, nlIOsetup) 

  end subroutine init_io_files
    

  !*****************************************************************************
  !*****************************************************************************

  subroutine global_io_init(input_shopping_list, full_shopping_list, &
                          & static_shopping_list, complete_static_shopping_list, &
                          & disp_dyn_shopping_list, disp_stat_shopping_list, &
                          & meteoMarketPtr, dispersionMarketPtr, outputMarketPtr, &
                          & met_buf, disp_buf, out_buf, &
                          & OutDef, &
                          & traj_set, &  !Container for trajectories
                          & wdr, diagnostic_rules, &
                          & chemRules, dynRules, &
                          & em_source, &
                          & cloud, &
                          & q_disp_dyn, q_disp_stat, &
                          & timestep, &                      ! of advection
                          & ResidenceTime, &                 ! in the grid due to advection
                          & PeriodToCompute, &
                          & iAccuracy, &
                          & DA_time_window, &
                          & nlSetupGrp, nlStdSetup)
    !
    ! Performs the global intialization of all input and output structures.
    ! Tasks:
    ! 1. Meteorology initialization
    ! 2. Conversion of emission source to pollution cloud
    ! 3. Adjustment of output parameters to the run setup and 
    !    intialization of output structures
    !
    implicit none

    ! Imported parameters
    type(silja_shopping_list), intent(inout) :: input_shopping_list, static_shopping_list, &
         & disp_dyn_shopping_list, disp_stat_shopping_list
    type(silja_shopping_list), intent(out) :: full_shopping_list, complete_static_shopping_list
    integer, intent(in) :: iAccuracy
    type(mini_market_of_stacks), pointer :: meteoMarketPtr, dispersionMarketPtr, outputMarketPtr
    type(silam_output_definition), pointer :: OutDef
    type(silja_wdr), pointer :: wdr
    type(Tdiagnostic_rules), intent(inout) :: diagnostic_rules
    type(Tchem_rules), intent(inout) :: chemRules
    type(Tdynamics_rules), intent(in) :: dynRules
    type(silam_source), pointer :: em_source
    type(silam_pollution_cloud), pointer :: cloud
    integer, dimension(:), intent(inout) :: q_disp_dyn, q_disp_stat
    type(silja_interval), intent(in) :: timestep, PeriodToCompute, DA_time_window
    type(Tfield_buffer), pointer :: met_buf, disp_buf, out_buf
    type(silja_interval), intent(out) :: ResidenceTime
    type(Tsilam_namelist_group), pointer :: nlSetupGrp
    type(Tsilam_namelist), pointer :: nlStdSetup
    TYPE(silam_trajectory_set), intent(out), target:: traj_set

    ! Local variables
    integer, dimension(:), pointer ::  q_dyn, q_st, req_dyn, req_st, quantities
    integer :: i,j, iTmp, al_status, nSrcId, nDatTypes,nCtlTypes, nBins, &
             & nTimeNodesNeeded, nTimeNodesNeededDisp, nx, ny 
    type(silam_sp) :: strTmp
    type(silja_grid) :: gridTmp, meteo_grid_reduced
!    type(silja_grid), pointer :: pGrid
    real :: corner_lon_E, corner_lat_N, s_pole_lon_E, s_pole_lat_N, dx_deg, dy_deg, fTmp, max_topo_height
    logical :: corner_in_geographical_latlon, if_south_pole
    real, dimension(:), pointer :: pTmp
    type(Tsilam_namelist), pointer :: nlSetup
    character (len=5) :: chTmp
    character (len=3) :: taskNumber
    character (len=clen) :: chOutGridType, chDispGridType
    type(silam_vertical) :: vertTmp
    type(silam_vertical), pointer :: pVertTmp
!    type(silam_output_definition), pointer :: od
    type(silja_wdr) :: wdrDisp
!    type(wdr_ptr), dimension(:), pointer :: wdrAr
    integer :: iNbrOfDispDynamicQ, iNbrOfDispStaticQ

    !DEBUG only
    integer, dimension(:), pointer :: iWork
    integer, dimension(:,:), pointer :: covermask
    integer ::  nx_gtm, ny_gtm

    strTmp%sp => fu_work_string()
    call set_OV_missing(OutVars)

!    ov => OutVars
!    od => OutDef
!    call msg("**************Global io_init got")
!    call msg ("Input")
!    call report (input_shopping_list)
!    call msg ("Static")
!    call report (static_shopping_list)
!    call msg("************End")

    OutDef%Params%tr_set => traj_set

    !
    ! Model input needs are known. Let's add what is needed for the IO server
    ! and then - be careful with request type. Here it is no longer trivial.
    !
    q_st => fu_work_int_array()
    req_st => fu_work_int_array()
    q_dyn => fu_work_int_array()
    req_dyn => fu_work_int_array()
    q_st = int_missing
    q_dyn = int_missing
    quantities => fu_work_int_array()

    call output_input_needs(OutDef, q_st, req_st, q_dyn, req_dyn, q_disp_dyn, q_disp_stat,  wdr)

    call add_shopping_quantities(input_shopping_list, q_dyn, req_dyn)
    call add_shopping_quantities(static_shopping_list, q_st, req_st)

    !
    ! If any quantity is 3D, we will need surface pressure
    !
    do iTmp = 1, fu_nbr_of_quantities(input_shopping_list)
      if(fu_multi_level_quantity(fu_quantity(input_shopping_list, iTmp)))then
        call add_shopping_quantities(input_shopping_list, (/surface_pressure_flag/), (/2/))
        exit
      endif
    end do

    !
    ! Prepare the log file. It will be opened in the directory of the first output,
    ! which means the first source, if they are split.
    ! Then, the currently collected log information will be copied to the new log file
    ! and the old log will be deleted
    !
    strTmp%sp = fu_FNm(OutDef%Rules%outTemplate, &
                     & fu_start_time(wdr), &    ! ini_time
                     & fu_start_time(wdr), &    ! anal_time
                     & zero_interval, & ! forecast length
                     & OutDef%Params%chCaseNm, &
                     & 'ALL_SRCS')   !log is the same for all sources
                   !OutDef%Params%chSrcNm(1)) cannot use source name - can be timezone-related 
    if(error)return                                                           ! and defined later

    !---------------------------------------------------------------
    !
    ! Having a common log file for MPI parallel runs is not practical,
    ! so for them we open a log file with the task number added just
    ! before the .log suffix.
    !
    iTmp = index(strTmp%sp,dir_slash,.true.)
    call create_directory_tree(strTmp%sp(1:iTmp-1))
    if(error)then
      call msg('Failed creating the directory from:' + strTmp%sp)
      call msg('Note dir_slash:' + dir_slash)
      return
    endif
    run_log_name = strTmp%sp(1:iTmp) + 'run_' + strTmp%sp(iTmp+1:len_trim(strTmp%sp))

    if(smpi_is_mpi_version())then
      call msg("Synchronizing tmp file prefix", i)
      write(taskNumber,'(I3.3)') smpi_global_rank
      i = int(silja_time_to_real8(fu_wallclock()))
      ! Just use what we have to ensure that the string is the same
      call msg("Synchronizing tmp file prefix", i)
      call smpi_allreduce_max_int(i, j, MPI_COMM_WORLD)
      run_log_tmp_name = run_log_name+'_' + &
                & fu_str(real8_to_silja_time(real(j,kind=8)),.false.)+ '_' + taskNumber
      run_log_name = run_log_name + '_' + taskNumber
    else
       run_log_tmp_name=run_log_name+'_'+ fu_str(fu_wallclock(),.false.)
    endif
    run_log_name=run_log_name+'.log'
    run_log_tmp_name=run_log_tmp_name+'.tmp.log'
    call msg("Preparing to swap log file to: "//trim(run_log_tmp_name))

    iTmp = fu_next_free_unit()
    open(iTmp, file = run_log_tmp_name, iostat = al_status)
    if(al_status /= 0)then
      call msg("Can't open for writing: "//trim(run_log_tmp_name))
      call set_error('Failed to open run.log file','global_io_init')
      return
    endif
    call copy_text_file(run_log_funit, iTmp)  ! Copy currently open log file

    close(run_log_funit,status = 'delete',  iostat = al_status)

    if(al_status /= 0)then  !!! Failed to close
      call msg("Can't delete file..., al_status = ", al_status)
      run_log_name=""
      inquire(unit=run_log_funit, NAME=run_log_name) !! going to crash anyway. Reuse variable to get the name..
      call msg("Can't delete file: "//trim(run_log_name))
      call set_error('Failed to remove old run.log file','global_io_init')
      return
    endif
    run_log_funit = iTmp

    !
    ! Time info file can be requested from the control file
    !
    nlSetup => fu_namelist(nlSetupGrp,'general_parameters')
    if(.not. associated(nlSetup))then
      call set_error('Failed to find general_parameters namelist in setup group','global_io_init')
      return
    endif
    strTmp%sp = fu_content(nlSetup,'time_info_file_name')
    if(len_trim(strTmp%sp) > 0)then
      info_funit = fu_next_free_unit()
      open(info_funit, file = strTmp%sp(1:iTmp) + strTmp%sp, iostat = al_status)
      if(al_status /= 0)then
        !
        ! It looks like info_file is open with another application. Create a new one with 
        ! some random value at the end
        !
        call random_number(fTmp)
        if(fTmp < 0.1) fTmp = fTmp + 0.1
        write(unit=chTmp,fmt='(I4)')int(fTmp*10000.)
        open(info_funit, file = strTmp%sp(1:iTmp) + strTmp%sp + chTmp, iostat = al_status)
        if(al_status /= 0)then  ! The problem is more serious, randomised file name does not help
          call set_error('Failed to open info_file_<random_number>','global_io_init')
          return
        endif
      endif
    else
      info_funit = int_missing
    endif  ! time info_file is requested

    !
    ! Grids available: emission, meteorology, output and internal
    !
    ! Emission grid: 
    ! decided by sources, whcih can be cut down to the dispersion grid. Must be inside
    ! the meteorological and dispersion grids
    !
    ! Meteo grid:
    ! max span is decided by available data, actual size is cut down to smallest
    ! area covering entirely the dispersion grid and sources, whatever is defined and 
    ! bigger
    ! 
    ! Output grid:
    ! Defined by the control file. Can be CUSTOM, which is the non-negotiable definition
    ! Can be EMISSION defined by the sources reduced down to meteorological grid, if needed
    ! Can be AREA_BASED constructed around the given area cut down to meteo grid, if needed
    ! Can be METEO copied from meteo grid, which is then built around the sources
    !
    ! Dispersion grid
    ! Must cover the sources and output and be covered by meteorology.
    ! Can be CUSTOM, which is non-negotiable definition.
    ! Can be EMISSION defined by the sources reduced down to meteorological grid, if needed
    ! Can be METEO copied from meteo grid
    ! Can be OUTPUT copied from output grid
    !
    ! Relations of the grids:
    ! Meteo_grid must cover the dispersion_grid
    ! dispersion_grid must cover all the sources
    ! Note that the output grid can be smaller or larger than the others. However,
    ! if it is flexible, we will force it inside the dispersion grid
    !
    ! Determine two grids: dispersion, if it is CUSTOM and non-negotiable and 
    ! the area to be asked from meteo file: envelope of the output grid and sources
    ! If dispersion and output grids are fixed, they have to be checked for 
    ! compatibility
    !
    ! First, output grid type for further references
    !
    nlSetup => fu_namelist(nlSetupGrp,'output_parameters')
    if(error .or. .not. associated(nlSetup))return
    chOutGridType = fu_str_u_case(fu_content(nlSetup,'grid_method'))
    !
    ! Now set the request for the meteo data area and, if possible, dispersion grid
    !
    if(OutDef%rules%ifRunDispersion) then
      if(dynRules%simulation_type == lagrangian_flag .or. dynRules%simulation_type == hybrid_flag)then
        !
        ! For Lagrangian advection, particles fly in meteo grid, hence dispersion_grid = meteo_grid
        ! However, the deposition is stored into the output_grid, which, together with source grid
        ! decides on the request for meteo.
        !
        chDispGridType = 'METEO_GRID'
!        dispersion_grid = output_grid ! whether it is defined or not
        gridTmp = output_grid
        call create_src_contain_grd_gen_src(em_source, gridTmp, &  ! enlarge to cover sources
                                          & .false., .false., .false.) ! ifVerbose, ifMinimal, ifInventoryOnly
        if(error)return

      elseif(dynRules%simulation_type == eulerian_flag)then  ! .or. dynRules%simulation_type == hybrid_flag)then
        !
        ! For Eulerian advection, we have all the niceties of dispersion_grid selection
        !
        nlSetup => fu_namelist(nlSetupGrp,'dispersion_parameters')
        if(error .or. .not.associated(nlSetup))then
          call set_error('dispersion_parameters namelist is absent','global_io_init')
          return
        endif
        chDispGridType = fu_str_u_case(fu_content(nlSetup,'grid_method'))
        !
        ! If dispersion_grid = CUSTOM_GRID, it is fixed
        !
        if(chDispGridType == 'CUSTOM_GRID')then
          !
          ! Dispersion grid is FIXED. 
          !
          dispersion_grid = fu_set_grid(nlSetup)

          if(error)return
          gridTmp = dispersion_grid

          !
          ! Output grid can already be defined by now. Check its relation with dispersion one
          !
          if(defined(output_grid))then
            if(.not. fu_if_grid_covered(output_grid, dispersion_grid))then   ! sml_grid, lrg_grid
              !
              ! Problem with coverage of output with dispersion
              !
              if(chOutGridType /= 'CUSTOM_GRID')then
                !
                ! Output is flexible: cut it.
                ! Call syntax: grid to cut, grid_area_template, cut type, grid_to_preserve
                !
                call cut_grid_size(output_grid, dispersion_grid, inside_the_grid_area)
                if(error)then
                  call set_error('Failed to reduce output grid down to dispersion one','global_io_init')
                  return
                endif
              endif  ! both grids are fixed
            endif  ! output is not covered by dispersion
          endif  ! defined output grid

        elseif(chDispGridType == 'METEO_GRID')then
          !
          ! Flexible dispersion definition, request for meteo is determined by output and source grids.
          !
          gridTmp = output_grid

          call create_src_contain_grd_gen_src(em_source, gridTmp, &  ! enlarge to cover sources
                                            & .false., .false., .false.) ! ifVerbose, ifMinimal, ifInventoryOnly
          if(error)return

        elseif(chDispGridType == 'OUTPUT_GRID')then
          !
          ! Output grid is the only one that decides on dispersion one. But what decides on the output? 
          !
          if(chOutGridType == 'METEO_GRID' )then
            !
            ! The only reasonable grid so far: from the sources
            !
            gridTmp = grid_missing
            call create_src_contain_grd_gen_src(em_source, gridTmp, &  ! enlarge to cover sources
                                              & .false., .false., .false.) ! ifVerbose, ifMinimal, ifInventoryOnly
          elseif(chOutGridType == 'CUSTOM_GRID')then
            !
            ! Output grid should be defined by now
            !
            gridTmp = output_grid

          else
            call set_error('Unknown output grid type:' + chOutGridType,'global_io_init')
            return
          endif

        else
          call set_error('Unknown dispersion grid type:' + chDispGridType,'global_io_init')
          return
        endif  ! if CUSTOM dispersion grid

      endif  ! lagrangian or eulerian advection

      !
      ! Temporal side of the story: data assimilation may require longer period of the data
      ! to be kept in memory
      !
      if(defined(DA_time_window))then
        nlSetup => fu_namelist(nlSetupGrp,'data_assimilation')
        if(error .or. .not.associated(nlSetup))then
          call set_error('dispersion_parameters namelist is absent','global_io_init')
          return
        endif
        
        nTimeNodesNeeded = int(DA_time_window / fu_obstime_interval(wdr) + 0.5) + 6
        nTimeNodesNeeded = max(min(nTimeNodesNeeded, max_times), 2)
      else
        nTimeNodesNeeded = 2   !3  ! no data assimilation
      endif
    else
      !
      ! No dispersion needed - go only with output grid (may still be undefined here)
      !
      gridTmp = output_grid
      nTimeNodesNeeded = 2

    endif  ! ifRunDispersion
    if(error)return

    !----------------------------------------------------------------------------------
    !----------------------------------------------------------------------------------
    !
    ! Initialize weather data usage. This routine from meteobuffer includes 
    ! permanent fields (physiography) and analysis of the GRIB files.
    ! It sets the meteo_grid so that it covers maximum possible part of 
    ! the requested area. However, nobody promised that it will cover 100%.
    !
    call msg('Initializing NWP data usage...')

    strTmp%sp = 'Meteo_initialisation'
    call start_count(chCounterNm = strTmp%sp)
    CALL meteo_init(wdr, &
                  & input_shopping_list, &
                  & full_shopping_list, &
                  & static_shopping_list, &
                  & complete_static_shopping_list, &
                  & fu_interval_positive(timestep), &    ! ifForward                  
                  & gridTmp, &                           ! Grid to cover with meteo (can be undefined)
                  & quantities, &  ! return list of quantities to be stored in meteo_buffer regardless st/dyn
                  & meteoMarketPtr, &   ! The main storage
                  & nTimeNodesNeeded)   ! number of meteo times to be kept in memory
    IF (error) RETURN



    !----------------------------------------------------------------------------------
    !
    ! Now we have the suggested biggest meteo_grid, which can be supported by the 
    ! meteorological data.
    ! The rest depends on the way the other grids are defined
    !
    if(OutDef%rules%ifRunDispersion) then
      !
      ! Run dispersion: set both output and dispersion grids and check/force sources inside
      !
      if(chDispGridType == 'CUSTOM_GRID')then  !----------------------- dispersion = CUSTOM -------------
        !
        ! Fixed dispersion grid. Check coverage with meteorology and set the output one
        !
        if(.not. fu_if_grid_covered(dispersion_grid, meteo_grid))then  ! small_grd, large_grd
          call msg('')
          call msg_warning('Meteo grid does not cover requested custom dispersion grid','grobal_io_init')
          call msg('requested grid: ')
          call report(gridTmp)
          call msg('Available meteo grid:')
          call report(meteo_grid)
          call set_error('Meteo grid does not cover requested custom dispersion grid','global_io_init')
          return
        endif
        !
        ! Fixed dispersion grid is covered by meteo - OK.
        ! Check sources, force inside if needed and allowed
        !
        gridTmp = dispersion_grid
!        call create_source_containing_grid(em_source, gridTmp, .false.)
        call create_src_contain_grd_gen_src(em_source, gridTmp, &  ! enlarge to cover sources
                                          & .false., .false., .false.) ! ifVerbose, ifMinimal, ifInventoryOnly
        if(error)return
        if(.not. gridTmp == dispersion_grid)then
          !
          ! Some sources are outside fixed dispersion grid
          ! 
          nlSetup => fu_namelist(nlSetupGrp,'emission_parameters')
          if(fu_str_u_case(fu_content(nlSetup,'cut_area_source_if_outside_meteo_grid')) == 'YES')then
            call msg('')
            call msg_warning('One of sources is (partially?) outside the fixed disperion grid. Cutting...', &
                           & 'global_io_init')
            !
            ! Careful: it may be not enough to force sources into the dispersion grid.
            ! We must force the output==emission grid into it first and then force sources
            ! Then different projections will not hirt us
            !
            call force_source_into_grid(em_source, dispersion_grid, .false.)
            if(error)return
            call msg('')
          else
            call set_error('One of sources is (partially?) outside the fixed disperion grid', &
                         & 'global_io_init')
            call msg('The disperion grid:')
            call report(dispersion_grid)
            call msg('The source-containing grid:')
            call report(gridTmp)
            return
          endif
        endif

        !
        ! The last grid to set - output (meteo, dispersion and sources are fitting one into another)
        !
        call set_output_grid(chOutGridType, meteo_grid, dispersion_grid, em_source, output_grid)
        if(error)return

      elseif(chDispGridType == 'METEO_GRID')then  !----------------------dispersion=METEO-----------
        !
        ! Dispersion grid is equal to meteorology
        !
        dispersion_grid = meteo_grid

        ! Cut dispersion_grid so any meteo field can be interpolated to stagger of dispersion --
        ! cut one cell on each edge
        call grid_dimensions(dispersion_grid, nx, ny)
        if (fu_ifLonGlobal(dispersion_grid)) then 
           call cut_grid_size(dispersion_grid, 1, 2, nx,   ny-1)
        else
           call cut_grid_size(dispersion_grid, 2, 2, nx-1, ny-1)
        endif

        !Lagrangian has hardcoded assumption that dispersion_grid = meteo_grid
        if(dynRules%simulation_type == lagrangian_flag) then
           meteo_grid = dispersion_grid
        endif

        ! Check if sources fit into dispersion grid
        gridTmp = dispersion_grid
        call create_src_contain_grd_gen_src(em_source, gridTmp, &  ! enlarge to cover sources
                                          & .false., .false., .false.) ! ifVerbose, ifMinimal, ifInventoryOnly
        if(error)return
        if(.not. gridTmp == dispersion_grid)then
          !
          ! Some sources are outside dispersion/meteo grid
          ! 
          nlSetup => fu_namelist(nlSetupGrp,'emission_parameters')
          if(fu_str_u_case(fu_content(nlSetup,'cut_area_source_if_outside_meteo_grid')) == 'YES')then
            call msg('')
            call msg_warning('One of sources is (partially?) outside the meteo grid. Cutting...', &
                           & 'global_io_init')
            call force_source_into_grid(em_source, dispersion_grid, .false.)
            if(error)return
            call msg('')
          else
            call set_error('Source-covering grid is not covered by meteorology', &
                         & 'global_io_init')
            call msg('The meteorological grid:')
            call report(meteo_grid)
            call msg('The source-covering grid:')
            call report(gridTmp)
            return
          endif
        endif
        !
        ! The last grid to set - output (meteo, dispersion and sources are fitting one into another)
        !
        call set_output_grid(chOutGridType, meteo_grid, dispersion_grid, em_source, output_grid)
        if(error)return

      elseif(chDispGridType == 'OUTPUT_GRID')then   !--------------------- dispersion = OUTPUT
        !
        ! Dispersion grid is requested to be the same as output. Hense, first set that grid
        !
        if(chOutGridType == 'CUSTOM_GRID')then
          !
          ! CUSTOM output grid; set above, just check the coverage by meteo
          !
          gridtmp =  fu_staggered_grid('output_stag', output_grid, if_stag_x=.true., if_stag_y=.true.)
          if (.not. fu_if_mesh_intrpolatable(gridtmp, meteo_grid,  ifSpk=.true.)) then
            !!!Only reportig is actually needed
            call cut_mesh(meteo_grid, gridtmp, meteo_grid_reduced, .True.)

            call msg_warning('Dispersion grid is requested = fixed output, which is outside meteo grid')
            call msg('Meteo grid:')
            call report(meteo_grid)
            call msg('Fixed output grid:')
            call report(output_grid)
            call set_error('Dispersion grid is requested = fixed output, which is outside meteo grid', &
                         & 'global_io_init')
            return
          endif

        elseif(chOutGridType == 'METEO_GRID')then
          !
          ! METEO grid used for the output. Most probably, exceeds dispersion_grid but so what?
          !
          output_grid = meteo_grid
          call grid_dimensions(output_grid, nx, ny)
          if (fu_ifLonGlobal(output_grid)) then 
            call cut_grid_size(output_grid, 1, 2, nx,   ny-1)
          else
            call cut_grid_size(output_grid, 2, 2, nx-1, ny-1)
          endif


        elseif(chOutGridType == 'AREA_BASED')then
          !
          ! Flexible output grid. Set above, here just force inside meteo one
          !
          if(.not. fu_if_grid_covered(output_grid, meteo_grid))then  ! small_grd, large_grd
            call cut_grid_size(output_grid, meteo_grid, inside_the_grid_area)
            if(error)return
          endif  ! output is not covered by meteo

        else
          call set_error('Unknown output grid type:' + chOutGridType,'global_io_init')
          return
        endif  ! output grid type

        !
        ! Having output grid set and forced into meteorology, set the dispersion one
        !
        dispersion_grid = output_grid

      else
        call set_error('Unknown dispersion grid type:' + chDispGridType,'global_io_init')
        return
      endif ! type of the dispersion grid definition

    else
      !
      ! No dispersion needed => no dispersion grid. Just output
      !
      dispersion_grid = meteo_grid
      call grid_dimensions(dispersion_grid, nx, ny)
      if (fu_ifLonGlobal(dispersion_grid)) then 
        call cut_grid_size(dispersion_grid, 1, 2, nx,   ny-1)
      else
        call cut_grid_size(dispersion_grid, 2, 2, nx-1, ny-1)
      endif

      if(chOutGridType == 'CUSTOM_GRID')then
        !
        ! CUSTOM output grid; set above, do nothing
        !
      elseif(chOutGridType == 'METEO_GRID')then
        !
        ! METEO grid used for the output. Most probably, exceeds dispersion_grid but so what?
        !
        output_grid = dispersion_grid  ! same as meteo-cut grid

      elseif(chOutGridType == 'AREA_BASED')then
        !
        ! Flexible output grid. Set above
        !
      else
        call set_error('Unknown output grid type:' + chOutGridType,'global_io_init')
        return
      endif  ! output grid type

    endif  ! if dispersion is requested

    !
    ! The output and dispersion grid reports
    !
    call msg('----------------------- Final OUTPUT GRID  -------------------------')
    call report(output_grid)
    call msg(' ')
    output_gridPtr => output_grid
    call grid_dimensions(output_grid, nx_output, ny_output)
    fs_output = nx_output * ny_output

    if(defined(dispersion_grid))then
      call msg('')
      if(smpi_is_mpi_version())then
        call msg('================== Whole-run DISPERSION GRID ===================')
      else
        call msg('======================= DISPERSION GRID ========================')
      endif
      call report(dispersion_grid)
      call msg('----------------------------------------------------------------')
      call msg('')
    else
      call set_error('Dispersion_grid is undefined','global_io_init')
      return
    endif

    ! Having the dispersion grid defined, let's check if it should be decomposed to MPI subdomains
    !
    wholeMPIdispersion_grid = dispersion_grid
    call grid_dimensions(wholeMPIdispersion_grid,  nx_wholeMPIdispersion, ny_wholeMPIdispersion)
    fs_wholeMPIdispersion = nx_wholeMPIdispersion * ny_wholeMPIdispersion

    call smpi_set_grids(wholeMPIdispersion_grid, dispersion_grid, output_grid, timestep)
    ! 
    ! From now on dispersion_grid is for this MPI member only!

    if(smpi_is_mpi_version())then
      call msg('')
      call msg('============== This MPI subdomain DISPERSION GRID ================')
      call report(dispersion_grid)
      call msg('----------------------------------------------------------------')
    end if

      
    ! We can now cut the meteo grid to include just the local dispersion area:
    !
    !In case of Lagrangian run meteo_grid == dispersion_grid the cur)grid will trigger an error
    if (.not. meteo_grid == dispersion_grid) then !! Workaround for lagrangian runs
      call msg('Selecting meteo grid (for this MPI sub domain)') 
      call cut_mesh(meteo_grid, fu_staggered_grid("Tmpstag", dispersion_grid, .True., .True.), &
                               & gridTmp, .false.)
      if (error) return
      meteo_grid = gridTmp
    endif

    if (error) return
    call msg('')
    call msg('============== This MPI subdomain METEOROLOGICAL GRID ================')
    call report(meteo_grid)
    call msg('----------------------------------------------------------------')
    ! repeat this now for possibly changed meteo grid
    call grid_dimensions(meteo_grid, nx_meteo, ny_meteo)
    fs_meteo = nx_meteo * ny_meteo
    meteo_gridPtr => meteo_grid
    !
    ! Time to store the system grid into the weather data rules. So far,
    ! there is some funny area there.
    !
    call set_storage_region(wdr,area_missing,meteo_grid)

    ! Finally, stupidity check: meteo grid must cover the dispersion one, MPI or not.
    !
    
    if(dynRules%simulation_type == eulerian_flag) then
      gridTmp = fu_staggered_grid('disp_stag', dispersion_grid, if_stag_x=.true., if_stag_y=.true.)
    else
      gridTmp = dispersion_grid
    endif

    if(.not.   fu_if_mesh_intrpolatable(gridTmp, meteo_grid)) then
      call set_error('After all, the meteo grid does not cover dispersion one','global_io_init')
      return
    endif

    ! Having the dispersion_grid defined, let's determine the derivatives
    ! of it - pole and dimensions, to speed-up further computations
    !
    dispersion_gridPtr => dispersion_grid
    call grid_dimensions(dispersion_grid,nx_dispersion,ny_dispersion)
    fs_dispersion = nx_dispersion * ny_dispersion
    if(error)return

    !-------------------------------------------------------------------------------
    !
    ! Output timing. Timestep can also coinside with the meteostep
    !
    OutDef%Rules%ini_time = fu_start_time(wdr)
    OutVars%LastDumpOutputTime = time_missing
    OutVars%LastRatesDumpOutputTime = time_missing

    if(defined(OutDef%Rules%dump_timestep)) &
        & OutVars%LastDumpOutputTime = OutDef%Rules%ini_time - OutDef%Rules%dump_timestep
    if(defined(OutDef%Rules%rates_dump_timestep)) &
        & OutVars%LastRatesDumpOutputTime = OutDef%Rules%ini_time - OutDef%Rules%rates_dump_timestep

    select case(OutDef%Rules%iOutTimesType)
      case(iMeteoTime)
        !
        ! Meteo time step
        !
        if(timestep > one_second) then  ! of advection
          OutDef%Rules%timestep = fu_obstime_interval(wdr)
          OutVars%LastOutputTime = fu_round_closest(OutDef%Rules%ini_time, &
                                                  & fu_obstime_interval(wdr), backwards)
        else
          OutDef%Rules%timestep = fu_obstime_interval(wdr) * (-1.)
          OutVars%LastOutputTime = fu_round_closest(OutDef%Rules%ini_time, &
                                                  & fu_obstime_interval(wdr), forwards)
        endif
        if(OutVars%LastOutputTime == OutDef%Rules%ini_time) &    ! Output immediately!
                 & OutVars%LastOutputTime = OutVars%LastOutputTime - OutDef%Rules%timestep
        OutDef%Rules%iOutTimesType = iRegular ! Now it is fixed and regular indeed

      case(iRegular)
        !
        ! Start time and last output time are set so that the output happens
        ! immediately. Timestep is already defined
        !
        OutVars%LastOutputTime = OutDef%Rules%ini_time - OutDef%Rules%timestep

      case default
        !
        ! Some special time => right now may be no output at all. It will happen only if 
        ! the computation start coinsides with this special event Use LastOutputTime as 
        ! temporary variable to determine the right output timestep.
        !
        if(timestep > one_second) then  ! of advection
          OutVars%LastOutputTime = fu_next_special_time(OutDef%Rules%ini_time, & ! - one_minute, & 
                                                      & OutDef%Rules%iOutTimesType, & ! type of occasion
                                                      & forwards, &
                                                      & timestep) !of advection: time quantum
        else
          OutVars%LastOutputTime = fu_next_special_time(OutDef%Rules%ini_time, & ! + one_minute, & 
                                                      & OutDef%Rules%iOutTimesType, & ! type of occasion
                                                      & backwards, &
                                                      & timestep*(-1.)) !of advection: time quantum
        endif
        if(defined(OutVars%LastOutputTime))then
          call msg(fu_connect_strings('First output will happen at:',fu_str(OutVars%LastOutputTime)))
        else
          call set_error('Failed to determine special-occasion output time','global_io_init')
          return
        endif
        if(OutVars%LastOutputTime == OutDef%Rules%ini_time)then
          !
          ! Ini time appeared to be the special occasion time.
          !
          OutDef%Rules%timestep = timestep  ! just for the cirrent moment
          OutVars%LastOutputTime = OutDef%Rules%ini_time - timestep
        else
          !
          ! No output at the beginning of the run!
          !
          OutDef%Rules%timestep = OutVars%LastOutputTime - OutDef%Rules%ini_time
          OutVars%LastOutputTime = OutDef%Rules%ini_time  ! Drop the temporary variable
        endif

    end select

    ! --------------------------------------------------------------
    !
    ! Set physiography (permanent fields) - the data source is not important
    ! here because they are scanned one-by-one starting from the first one,
    ! which is the "most important".
    ! Note that we need the topography field right now - requested it in output_input_needs sub
    !
    CALL set_physiography(meteoMarketPtr, wdr, static_shopping_list, complete_static_shopping_list, iAccuracy)
    if(error)return
    pTmp => fu_grid_data(topography_fld)
    max_topo_height = maxval(pTmp(1:fu_number_of_gridpoints(meteo_grid)))

    !
    ! If the sources are to be split over time zones, this is the right time to define the mapping.
    !
    call make_source_id_mapping(em_source, meteoMarketPtr, dispersion_grid, &
                              & fu_namelist(nlSetupGrp, 'output_parameters'))
    if(error)return

    !-----------------------------------------------------------------
    !
    ! Allocate the space in accordance with the number of sources and 
    ! ifSplitSources
    !
    ! Number of sources
    !
    if(OutDef%rules%ifRunDispersion)then
      nSrcId = fu_NbrOf_source_ids(em_source) ! Number of separate ids
    else
      nSrcId = 1
    endif
    if(error)return

    !
    ! Allocation of file unit numbers, grib/grads ctl structures and file names
    !
!   Netcdf stuff added by Jukkis, 26032008
    allocate(OutDef%Params%grib_funit(nSrcId), OutDef%Params%grads_funit(nSrcId), &
           & OutDef%Params%netcdf_funit(nSrcId), &
           & OutDef%Params%grib_fNm(nSrcId),   OutDef%Params%grads_fNm(nSrcId), &
           & OutDef%Params%netcdf_fNm(nSrcId), &
           & OutDef%Params%grib_struct(nSrcId),OutDef%Params%grads_struct(nSrcId), &
           & OutDef%Params%netcdf_struct(nSrcId), &
           & OutDef%Params%chSrcNm(nSrcId), stat=al_status)
! Netcdf stuff added by Jukkis, 26032008
    if(fu_fails(al_status == 0, 'Cannot allocate output suppl. variables', 'global_io_init'))return

    OutDef%Params%grib_funit = int_missing
    OutDef%Params%grads_funit = int_missing
    OutDef%Params%netcdf_funit = int_missing

    !
    ! Source names are to be connected if their stuff is mixed.
    ! Then, if the case is called after source and there is one source
    ! or the sources are mixed togerther - the case takes source' name 
    ! while source gets nothing
    !
    if(OutDef%rules%ifRunDispersion)then
      do i=1,nSrcId
        OutDef%Params%chSrcNm(i) = fu_get_id_str_from_id_nbr(em_source,i)
      end do
    else
      OutDef%Params%chSrcNm(1) = 'meteo'  ! no sources if no dispersion
    endif
    !
    ! The case name can be after the source ID. In this situation, clear the source ID
    !
    if(OutDef%Params%chCaseNm == '-')then
      OutDef%Params%chCaseNm = OutDef%Params%chSrcNm(1)
      OutDef%Params%chSrcNm(1) = ''
    endif

    !
    ! Number of the output files to be written comes from run_duration, timeStep of 
    ! output, number of separate sources and number of output file formats. 
    !
    nDatTypes=0
    nCtlTypes=0
    if(OutDef%Rules%iGrib > 0) then ! GRIB itself and ctl
      nDatTypes=nDatTypes+1  
      nCtlTypes=nCtlTypes+1
    endif
    if(OutDef%Rules%ifGrads)then ! GrADS itself and ctl
       nDatTypes=nDatTypes+1
       nCtlTypes=nCtlTypes+1
    endif
    if(OutDef%Rules%iNetcdf > 0)then ! netcdf itself and ctl
       nDatTypes=nDatTypes+1
       nCtlTypes=nCtlTypes+1
    endif
    if(OutDef%Rules%ifTrajectory) nDatTypes=nDatTypes+1  ! just trajectory file

    !
    ! Computation of the number of files to be written:
    ! nSrcId source Ids, nDatTypes data types, nCtlTypes ctl types, 
    ! nBins binaries
    ! So: n= nSrc * (nDatTypes*nBinaries + nCtlTypes)
    !
    select case(OutDef%Rules%OutFilesArrangement)
      case(all_in_one)
        nBins = 1  ! One binary only for one output
      case(hourly_new_file)
        nBins = max(1,int(0.5+fu_period_to_compute(wdr)/one_hour))
      case(daily_new_file)
        nBins = max(1,int(0.5+fu_period_to_compute(wdr)/one_day))
      case(monthly_new_file)
        nBins = max(1,int(0.5+fu_period_to_compute(wdr)/one_day/30.))
      case(yearly_new_file)
        nBins = max(1,int(0.5+fu_period_to_compute(wdr)/one_day/365.))
    end select
    allocate(OutDef%Params%filesWritten(nSrcId*(nDatTypes*nBins+nCtlTypes)+nDatTypes*2+1),stat=al_status)
    if(fu_fails(al_status == 0, 'Failed to allocate list of written files', 'global_io_init'))return
    OutDef%Params%filesWritten = ''

    !-----------------------------------------------------------------------------
    !
    ! Vertical definitions
    !
    ! The output levels can be requested from meteo_vertical too
    !
    if(fu_str_u_case(fu_content(fu_namelist(nlSetupGrp, 'output_parameters'),'vertical_method')) &
     & == 'METEO_LEVELS')then
      output_vertical = meteo_vertical
     endif

    !
    ! Dispersion vertical structure can be partly determined by user for Eulerian advection
    ! while for Lagrangian we force selection of the meteo cut to the height of ouput one
    !
    nlSetup => fu_namelist(nlSetupGrp,'dispersion_parameters')
    if(fu_fails(associated(nlSetup), 'dispersion_parameters namelist is absent','global_io_init'))return

    if(dynRules%simulation_type == lagrangian_flag)then

!      dispersion_vertical = output_vertical
      call msg('Lagrange requires dispersion vertical based on meteo levels')
      call levels_to_layers(meteo_vertical, dispersion_vertical)
      ! dispersion vertical top must be below the upper meteo level
      call remove_last_level(dispersion_vertical)

      ! If no output from upper layers -- can cut further
      if(fu_if_cut_vertical_size_from_above(dispersion_vertical, output_vertical, max_topo_height))then
        call msg("Cut done")
      endif

    elseif(dynRules%simulation_type == eulerian_flag .or. dynRules%simulation_type == hybrid_flag)then
      !
      ! Eulerian advection, vertical is selected by the user
      !
      if(fu_str_u_case(fu_content(nlSetup,'vertical_method')) == 'METEO_LEVELS')then
        !dispersion_vertical = meteo_vertical
        call msg('Creating dispersion vertical from meteo levels')
        call levels_to_layers(meteo_vertical, dispersion_vertical)

      elseif(fu_str_u_case(fu_content(nlSetup,'vertical_method')) == 'OUTPUT_LEVELS')then
        dispersion_vertical = output_vertical

      elseif(fu_str_u_case(fu_content(nlSetup,'vertical_method')) == 'CUSTOM_LEVELS' .or. &
           & fu_str_u_case(fu_content(nlSetup,'vertical_method')) == 'CUSTOM_LAYERS')then
        call set_vertical(nlSetup, dispersion_vertical)

      else
        call set_error('Strange dispersion vertical deifinition:' + &
                     & fu_content(nlSetup,'vertical_method'),'global_io_init')
        return
      endif
    endif  ! type of advection
    if(error .or. .not.defined(dispersion_vertical)) then
      call report(dispersion_vertical)
      call set_error('Failed to determine the dispersion_vertical','global_io_init')
      return
    endif

    !
    ! Now, with the domain topography defined, we can reduce the number of vertical levels 
    ! in meteo grid. Possible only if bottom-up wind diagnostic is used, of course.
    !
    if(.not. fu_if_full_meteo_vertical_needed(diagnostic_rules))then
      vertTmp = meteo_vertical
      if(fu_if_cut_vertical_size_from_above(vertTmp, dispersion_vertical, max_topo_height))then
        call replace_vertical_in_list(input_shopping_list, meteo_vertical, vertTmp)
        call replace_vertical_in_list(full_shopping_list, meteo_vertical, vertTmp)
        if(error)then
          call set_error('Failed vertical cut','global_io_init')
          return
        endif
        meteo_vertical = vertTmp
        call vertical_parameters_ab(meteo_vertical, nz_meteo, a_met, b_met)
      endif
    endif

    ! Now some forcing: dispersion_vertical MUST be either hybrid or height layers.
    !
    if (.not. any(fu_leveltype(dispersion_vertical) == (/layer_btw_2_hybrid, layer_btw_2_height/))) then
      call set_error('Only hybrid or height layers are allowed for dispersion vertical', 'global_io_init')
      return
    endif

    dispersion_verticalPtr => dispersion_vertical
    call vertical_parameters(dispersion_vertical, nz_dispersion, &
                       a_half_disp, b_half_disp,  disp_layer_top_m, .true., .true.) 
                       !    & dispersion_layer_z_size, dispersion_layer_boundaries)
    if(error)return

    ! Due to some internal limitations, we require nz_dispersion > 1.
    !
    if(fu_fails(nz_dispersion > 1,'At least 2 levels required in dispersion vertical','global_io_init'))return
    call msg('')
    call msg('========================= METEO VERTICAL ==========================')
    call report(meteo_vertical,.true.)
    call msg('----------------------------------------------------------------')
    call msg('')
    call msg('======================= DISPERSION VERTICAL ========================')
    call report(dispersion_vertical,.true.)
    call msg('----------------------------------------------------------------')
    call msg('')


    call msg('========================= OUTPUT VERTICAL ==========================')
    call report(output_vertical,.true.)
    call msg('----------------------------------------------------------------')
    call msg('')

    output_verticalPtr => output_vertical
    call vertical_parameters(output_vertical, nz_output, &
                       a_half_out, b_half_out, out_layer_top_m, .true., .true.) 

    ! Initalize realtime fields. Permanent ones were already set in physiography
    !
    call init_singletime_fields(meteoMarketPtr, meteo_grid , meteo_vertical, complete_static_shopping_list)

    !
    ! For the dispersion_vertical, we will need coresponding vertical velocity
    !
    select case(fu_leveltype(dispersion_vertical))
      case(constant_pressure, layer_btw_2_pressure, &
         & sigma_level, layer_btw_2_sigma, &
         & hybrid, layer_btw_2_hybrid)

        vertical_velocity_pointer = omega_flag

      case(constant_altitude, layer_btw_2_altitude)

        vertical_velocity_pointer = w_alt_msl_flag

      case(constant_height, layer_btw_2_height, &
         & depth_level, layer_btw_2_depth)

        vertical_velocity_pointer = w_height_srf_flag

      case(entire_atmosphere_single_layer)

        vertical_velocity_pointer = int_missing

      case default
        call set_error('Cannot select proper vertical velocity','global_io_init')
        return
    end select
    if (error) return
    !
    ! Having the vertical type selected, we can refine the selection of the vertical wind
    ! in the full_shopping_list. So far, I just replace teh vertical_velocity_flag
    ! with the vertical_velocity_pointer.
    !
!    call msg('')
!    call msg('Replacing the vertical wind pointer with its value for the selected vertical')
!    call msg('Input shopping list ONCE AGAIN')
!    call report(input_shopping_list)
    call replace_quantity(full_shopping_list, vertical_velocity_flag, vertical_velocity_pointer, .true.)
    call replace_quantity(input_shopping_list, vertical_velocity_flag, vertical_velocity_pointer, .true.)
    if(error)return

!    call msg('... and shopping list after addition/replacement of vertical wind...')
!    call report(full_shopping_list)
!    call msg('')


    call arrange_supermarket(meteoMarketPtr)


    call stop_count(chCounterNm = strTmp%sp)
    

    !--------------------------------------------------------------------
    !
    ! Initialise the dispersion market and buffer
    !
    if(OutDef%rules%ifRunDispersion)then
        ! 
        ! Now include the vertical-dependent quantities to the list of
        ! requested dispersion quantities.
        !
        if (dynRules%simulation_type == eulerian_flag .or. dynRules%simulation_type == hybrid_flag) then
        call euler_adv_input_needs(dynRules%advMethod_Eulerian, &
                                 & .true., &              ! include vertical-dependent fields
                                 & ifUseMassfluxes(diagnostic_rules), &
                                 & q_dyn, q_st, q_disp_dyn, q_disp_stat)
        if(error)return
        endif

        if(fu_if_level_meteo_dependent(fu_leveltype(dispersion_vertical)))then
          iTmp = fu_merge_int_arrays((/cell_size_z_flag, int_missing/), q_disp_dyn, .true.)
          iTmp = fu_merge_int_arrays((/cell_size_x_flag, cell_size_y_flag, int_missing/), &
                                   & q_disp_stat, .true.)
        else
          iTmp = fu_merge_int_arrays((/cell_size_x_flag, cell_size_y_flag, cell_size_z_flag, &
                                     & int_missing/), q_disp_stat, .true.)
        endif


        call add_shopping_quantities(disp_dyn_shopping_list, q_disp_dyn)
        call add_shopping_quantities(disp_stat_shopping_list, q_disp_stat)
!    call msg("") 
!    call msg("*********************************************************")
!    call msg("disp_dyn_shopping_list after euler_adv_input_needs")
!    call report(disp_dyn_shopping_list)
!    call msg("*********************************************************")
!    call msg("") 
!    call msg("") 
!    call msg("*********************************************************")
!    call msg("disp_stat_shopping_list after euler_adv_input_needs")
!    call report(disp_stat_shopping_list)
!    call msg("*********************************************************")
!    call msg("")
!    call ooops("")

      wdrDisp = wdr
      call set_storage_region(wdrDisp, area_missing, dispersion_grid)
      !
      ! How many time nodes are to be kept in the market?
      !
      if(defined(DA_time_window))then
        nTimeNodesNeededDisp = int(DA_time_window / fu_obstime_interval(wdrDisp) + 0.5) + 6
        nTimeNodesNeededDisp = max(min(nTimeNodesNeededDisp,max_times),2)
      else
        nTimeNodesNeededDisp = 2  ! no data assimilation
      endif

      ! Initialize mini market...
      call init_supplementary_market(dispersionMarketPtr, &
                                   & q_disp_dyn, q_disp_stat, &
                                   & wdrDisp, nTimeNodesNeededDisp, 'dispersion minimarket')
      if(error)return

      !
      ! Now, fill-in the dispersion market with the fields that will not change their location
      !
      ! ... fill-in the basic information: dx, dy, dz ...
      !
      call set_basic_physiography(dispersionMarketPtr, wdrDisp, dispersion_grid, dispersion_vertical)
      if(error)return

      ! Allocate structures needed for wind diagnostics
      call init_wind_diag(wholeMPIdispersion_grid, dispersion_grid, dispersion_vertical, &
                  &  meteo_vertical, diagnostic_rules)
      if(error)return

      ! ... create empty fields to store realtime stuff...
      !
      call init_singletime_fields(dispersionMarketPtr, dispersion_grid, &
                                & dispersion_vertical, disp_stat_shopping_list)
      if(error)return

      ! Should this be needed, initialise the water in soil computation procedure
      !
      call init_water_in_soil_model(dispersionMarketPtr, meteoMarketPtr, diagnostic_rules, fu_start_time(wdr))
      if(error)return
      
      ! Initialise emission sources and related fields: they are in dispersionMarket
      !
      call init_emission_internal_fields(em_source, dispersionMarketPtr, fu_start_time(wdr))
      if(error)return
      !
      ! Arrange the market
      !
      call arrange_supermarket(dispersionMarketPtr)
      if(error)return

      ! Finally, initialise the corresponding buffer...
      !
      nullify(disp_buf)
      q_dyn = int_missing
      iTmp = fu_merge_int_arrays(q_disp_dyn, q_dyn, .true.)
      iTmp = fu_merge_int_arrays(q_disp_stat, q_dyn, .true.)

      call init_supplementary_data_buffer(dispersionMarketPtr, disp_buf, &
                                        & .false., fs_dispersion, q_dyn)
      if(error .or. fu_fails(associated(disp_buf),'Non-associated dispersion buffer', &
                                                & 'global_io_init'))return
      ! ... and arrange it
      !
      call arrange_buffer(dispersionMarketPtr, met_src_missing, time_missing, interval_missing, disp_buf,& 
                                       &  dispersion_verticalPtr)
      if(error)return
      !
      ! A couple of other global pointers
      !
      dispersion_cell_x_size_fld => fu_sm_simple_field(dispersionMarketPtr, silam_internal_src,&
                                                     & cell_size_x_flag,&
                                                     & level_missing, &
                                                     & single_time_stack_flag)
      dispersion_cell_y_size_fld => fu_sm_simple_field(dispersionMarketPtr, silam_internal_src,&
                                                     & cell_size_y_flag,&
                                                     & level_missing, &
                                                     & single_time_stack_flag)

      !----------------------------------------------------------------------------
      !
      ! Converting source to initial cloud.
      !
      strTmp%sp = 'Cloud initialisation'
      call start_count(chCounterNm = strTmp%sp)
      call msg(strTmp%sp)

      !
      ! For both Lagrangian and Eulerian advection, we need to know the residence time.
      ! Lagrangian model will guess the single-particle mass, Eulerian will set low-mass threshold
      !
      if (dynRules%simulation_type == eulerian_flag .or. dynRules%simulation_type == hybrid_flag) then
         ResidenceTime = fu_residence_interval(dispersion_grid) ! Positive!
      else
         ResidenceTime = fu_residence_interval(meteo_grid) ! Positive!
      endif

      if(error)return
      !
      ! Whatever the timestep is, the residence interval must be an integer number of hours
      ! to avoid clash with source time variation coefficients
      !
      ResidenceTime = one_hour * max(1., 2. * int(abs(ResidenceTime / one_hour) * 0.5 + 0.5))
      !
      ! Timestep can be < 0 for adjoint runs, so its sign goes to reset interval
      !
      ResidenceTime = timestep * max(1., 2. * int(abs(ResidenceTime / timestep) * 0.5 + 0.5))

      call msg('Estimated residence time is:' + fu_str(ResidenceTime))

      !
      ! Regardless the advection type cloud is responsible for emission, concentration and 
      ! optical density fields. Therefore, all necessary information has to be given.
      ! Emission is in source, chemistry is in source and standard setup namelist,
      ! optical basic data are sitting in the standard namelist, plus the trigger from
      ! output or model input data
      !
      CALL source_to_initial_cloud(em_source, cloud, chemRules, dynRules, &
                                 & OutDef%Rules%ini_time, fu_period_to_compute(wdr), timestep, &
                                 & meteoMarketPtr, &
                                 & output_gridPtr, &
                                 & iAccuracy, fu_if_randomise(wdr))
      IF (error) RETURN

      ! Will we need reaction rates dump?
      if (defined(OutDef%Rules%rates_dump_timestep)) then
         call set_react_rates_map(cloud)
      endif

    endif  ! ifRunDispersion
    IF (error) RETURN
    




    call stop_count(chCounterNm = strTmp%sp)


    !--------------------------------------------------------------------
    !
    ! Initialise the output market
    !
    wdrDisp = wdr
    call set_storage_region(wdrDisp, area_missing, output_grid)
    !
    ! So far, we need only dx, dy, dz fields. 
    !
    q_dyn = int_missing
    q_st = int_missing
    if(fu_if_level_meteo_dependent(fu_leveltype(output_vertical)))then
      q_dyn(1:3) = (/cell_size_z_flag, air_density_flag, int_missing/)
      ! density not really needed but diagnostic quantities will put it there with z size.
      q_st(1:3) = (/cell_size_x_flag, cell_size_y_flag, int_missing/)
    else
      q_dyn(1) = int_missing
      q_st(1:4) = (/cell_size_x_flag, cell_size_y_flag, cell_size_z_flag, int_missing/)
    endif

    ! Initialize mini market...
    call init_supplementary_market(outputMarketPtr, q_dyn, q_st, wdrDisp, 2, 'output minimarket')
    if(error)return

    ! ... and the corresponding buffer...
    nullify(out_buf)
    iTmp = fu_merge_int_arrays(q_st, q_dyn, .true.)
    call init_supplementary_data_buffer(outputMarketPtr, out_buf, &
                                      & .false., fs_output, q_dyn)
    if(error .or. fu_fails(associated(out_buf),'Non-associated output buffer', &
                                             & 'global_io_init'))return

    ! ... fill-in the basic information: dx, dy, dz ...
    call set_basic_physiography(outputMarketPtr, wdrDisp, output_grid, output_vertical)
    if(error)return

    call arrange_supermarket(outputMarketPtr)
    if(error)return

    ! ... and set permanent dispersion pointers
!    CALL set_permanent_buffer_pointers(outputMarketPtr, fu_start_time(wdrDisp), out_buf, &
!                                    & output_gridPtr, output_verticalPtr)
    call arrange_buffer(outputMarketPtr, met_src_missing, time_missing, interval_missing, out_buf,& 
                                       & output_verticalPtr)
    if(error)return

    call free_work_array(q_dyn)
    call free_work_array(q_st)
    call free_work_array(req_dyn)
    call free_work_array(req_st)

    !---------------------------------------------------------------------------
    !
    ! Now all grids seem to be determined, so we can compute all the coefficients
    ! for cross-interpolation.
    !
    call msg('')
    call msg('Setting horizontal and vertical cross-interpolation')
    !
    ! Meteo-to-dispersion grids and verticals. Conversion is set in 
    ! pollution_cloud
    !
    if(OutDef%rules%ifRunDispersion) then
      call msg('Meteo to dispersion grid interpolation')
      call set_meteo2disp_interp(cloud, wdr, timestep) ! from pollution cloud
      if(error)return
    endif

    !
    ! Meteo-to-output grid
    !
    call msg('Meteo to output grid interpolation')
    if(output_grid == meteo_grid)then
      nullify(OutVars%interpCoefMeteo2OutHoriz)
      OutDef%Rules%ifMeteo2OutHorizInterp = .false.
    else
      OutVars%interpCoefMeteo2OutHoriz => fu_horiz_interp_struct(meteo_grid, &    ! grid From
                                                               & output_grid, &  ! grid To
                                                               & fu_horizontal_interp_method(wdr),&
                                                               & fu_if_randomise(wdr), &
                                                               &  iOutside = setMissVal) !method
      if(error)then
        call set_error('Failed meteo to output horizontal interpolation structure','global_io_init')
        return
      endif
      OutDef%Rules%ifMeteo2OutHorizInterp = .true.
    endif
    if(error)return
    !
    ! Meteo-to-output verticals
    !
!call msg('Meteo to output vertical interpolation')
    if(fu_cmp_verts_eq(output_vertical, meteo_vertical))then
      nullify(OutVars%interpCoefMeteo2OutVert)
      OutDef%Rules%ifMeteo2OutVertInterp = .false.
    else
      OutVars%interpCoefMeteo2OutVert => fu_vertical_interp_struct(meteo_vertical, &
                                                                 & output_vertical, &
                                                                 & output_grid, &
                                                                 & fu_vertical_interp_method(wdr), &
                                                                 & one_hour, 'main_meteo_to_output')
      OutDef%Rules%ifMeteo2OutVertInterp = .true.
    endif
    if(error)return

    !
    ! Dispersion-to-output grid
    !
    if(OutDef%rules%ifRunDispersion) then
!call msg('Dispersion to output grid interpolation')
      if(output_grid == dispersion_grid)then
        nullify(OutVars%interpCoefDisp2OutHoriz)
        OutDef%Rules%ifDisp2OutHorizInterp = .false.
      else
        OutVars%interpCoefDisp2OutHoriz => fu_horiz_interp_struct(dispersion_grid, &    ! grid From
                                                               & output_grid, &  ! grid To
                                                               & fu_horizontal_interp_method(wdr), &
                                                               & fu_if_randomise(wdr), &
                                                               &  iOutside = setMissVal) !method
        if(error)then
          call set_error('failed dispersion to output horizontal interpolation structure', &
                       & 'global_io_init')
          return
        endif
        OutDef%Rules%ifDisp2OutHorizInterp = .true.
      endif ! output grid == dispersion grid

      !
      ! Dispersion-to-output verticals
      !
!call msg('Dispersion to output vertical interpolation')
      if(fu_cmp_verts_eq(output_vertical, dispersion_vertical))then
        nullify(OutVars%interpCoefDisp2OutVert)
        OutDef%Rules%ifDisp2OutVertInterp = .false.
      else
        OutVars%interpCoefDisp2OutVert => fu_vertical_interp_struct(dispersion_vertical, &
                                                                  & output_vertical, &
                                                                  & output_grid, &
                                                                  & fu_vertical_interp_method(wdr), &
                                                                  & OutDef%rules%timestep, &
                                                                  & 'disp_to_output')
        OutDef%Rules%ifDisp2OutVertInterp = .true.
      endif ! output vertical == dispersion vertical
      if(error)return
      call report(cloud, .false.)  ! Non-detailed report
      call collect_total_masses(cloud)
      call report_total_masses(cloud, OutVars%iNbrCollection, .true.)
!      call report_inout_mass_cld(cloud,outgoing)
      call msg('')
    endif  ! if run dispersion

    !---------------------------------------------------------------------------
    !
    ! If total emission is expected in the output, let's write it right now 
    ! and remove from the output list before tune_output_parameters allocates
    ! memory in the output stacks.
    ! A trick: e.g. for biological sources we do not know the emission because
    ! it depends on meteorological parameters. However, we do know the forest map
    ! and threshold levels, which may come as surrogates of the emission "totals"
    ! Note: if the emission output is requested in a dynamical form - welcome.
    ! It will stay in the list and set_emission_output will do nothing.
    !
    if(OutDef%rules%ifRunDispersion)then
      do i=1, size(OutDef%Rules%DispOutLst%ptrItem)
        if(OutDef%Rules%DispOutLst%ptrItem(i)%quantity == emission_intensity_flag)then
          call msg('Setting emission output')
          call set_emission_output(OutDef, em_source, cloud, PeriodToCompute, iAccuracy, &
                                 & fu_if_randomise(wdr))
          if(error)return
          exit
        endif
      enddo 
    endif
    !
    ! If permanent fields are requested to the output, let's also write them
    ! here and remove from the list.
    !
    if(allocated(OutDef%Rules%DispOutLst%ptrItem))then
      do i = 1, size(OutDef%Rules%DispOutLst%ptrItem)

        if(OutDef%Rules%DispOutLst%ptrItem(i)%quantity == physiography_field_set_flag) then
          !
          ! Force the permanent stack to the output
          !
          call msg('Setting physiography output')
          call write_physiography_output(meteoMarketPtr, &
                                       & fu_FNm(OutDef%Rules%outTemplate, &
                                              & OutDef%Rules%ini_time, &  ! ini_time
                                              & OutDef%Rules%ini_time, &  ! anal_time
                                              & zero_interval, &      ! forecast length
                                              & OutDef%Params%chCaseNm, &
                                              & OutDef%Params%chSrcNm(1)), &  
                                       & OutDef%Rules%iGrib, &
                                       & OutDef%Rules%ifGrADS, &
                                       & OutDef%Rules%iNetCDF, &
                                       & fu_start_time(wdr), &
                                       & OutDef%Params%filesWritten, &
                                       & output_grid, &
                                       & fu_if_randomise(wdr))
          if(error) call unset_error('global_io_init')
          !
          ! Remove the quantity from the list
          !
          call shrink_output_list(OutDef%Rules%DispOutLst, i)

          if(error)call unset_error('global_io_init')
        endif
      enddo
    endif

    !------------------------------------------------------------------------
    !
    ! Finally, all is set, permanent fields are introduced, 
    ! so we can initiate buffers and try to set permanent pointers. 
    !
    nullify(met_buf)
    call init_data_buffer(met_buf, quantities, .true., fs_meteo)  
    if(error .or. .not. associated(met_buf))then
      call set_error('Non-associated meteo buffer','global_io_init')
      return
    endif
    !CALL set_permanent_buffer_pointers(meteoMarketPtr, fu_start_time(wdr), met_buf, &
     !                                & meteo_gridPtr, meteo_verticalPtr)
    call arrange_buffer(meteoMarketPtr, met_src_missing, time_missing, interval_missing, met_buf,& 
                                       &  meteo_verticalPtr)
    if(error)return

!    call msg('After physiography cloud status=',fu_overall_cloud_status(cloud))

    !-----------------------------------------
    !
    ! Some parameters of the output may depend on the results of the meteo 
    ! data analysis. So, let's adjust these parameters.
    ! Also, the output may contain nuclides, whose availability is known only
    ! when the pollution cloud is created.That's the second task
    !
    call msg('')
    call msg('')
    call msg('Tuning output parameters')

    call tune_output_parameters(wdr, &
                              & full_shopping_list, &
                              & static_shopping_list, &
                              & cloud, &
                              & em_source, &
                              & OutDef, chemRules, dynRules, nlStdSetup, &
                              & DispersionMarketPtr, met_buf, &
                              & timestep) ! Of the model
    if(error)return
    !
    ! All done, stop the counters
    !
    call msg('')
    call msg('============================================================')
    call msg('')
    call msg_test('Done global IO initalization')
    call msg('')
    call msg('============================================================')
    call msg('')



!call report(DispersionMarketPtr)

    
    call free_work_array(quantities)
    call free_work_array(strTmp%sp)

!    call msg('After global IO initialization cloud status=',fu_overall_cloud_status(cloud))
!    call msg('End of global_io_init')

    CONTAINS

    subroutine set_output_grid(chOutGridType, meteoGrd, dispGrd, emSrc, outGrd)
      !
      ! Encapsulates the output grid setup
      !
      implicit none

      ! Imported parameters
      character(len=*), intent(in) :: chOutGridType
      type(silja_grid), intent(in) :: meteoGrd, dispGrd
      type(silam_source), intent(in) :: emSrc
      type(silja_grid), intent(inout) :: outGrd

      if(chOutGridType == 'CUSTOM_GRID')then
        !
        ! CUSTOM output grid; set above, do nothing
        !
      elseif(chOutGridType == 'METEO_GRID')then
        !
        ! METEO grid used for the output. Most probably, exceeds dispersion_grid but so what?
        !
        outGrd = meteoGrd

      elseif(chOutGridType == 'AREA_BASED')then
        !
        ! Flexible output grid. Set above, here just force inside dispersion one
        !
        if(.not. fu_if_grid_covered(outGrd, dispGrd))then  ! small_grd, large_grd
          call cut_grid_size(outGrd, dispGrd, inside_the_grid_area)
          if(error)return
        endif  ! output is not covered by dispersion

      else
        call set_error('Unknown output grid type:' + chOutGridType,'set_output_grid')
        return
      endif  ! output grid type

    call msg("")
    call msg("outdef%Rules%dispoutlst%ptritem at the end global_io_init")
    do i=1, size(outdef%Rules%dispoutlst%ptritem)
      call report(outdef%Rules%dispoutlst%ptritem(i)%targetId)
    enddo

    end subroutine set_output_grid


  end subroutine global_io_init

  !**********************************************************************************
    

  subroutine prepare_meteo( wdr, iAccuracy, DiagRules,&
                          & meteo_input_dyn_shopping_list, meteo_full_dyn_shopping_list, &
                          & disp_dyn_shopping_list, disp_stat_shopping_list, &
                          & output_dyn_shopping_list, output_stat_shopping_list, &
                          & meteo_ptr, disp_buf_ptr, output_buf_ptr, &
                          & now, timestep, &
                          & meteoMarketPtr, dispersionMarketPtr, outputMarketPtr, ifGotNewMeteoData)

    ! Acquire needed meteo and set all buffers

    !
    implicit none
    
    ! Imported parameters
    !
    type(Tdiagnostic_rules) :: DiagRules
    integer, intent(in) :: iAccuracy
    type(silja_time), intent(in) :: now
    type(silja_interval), intent(in) :: timestep
    type(silja_wdr), pointer :: wdr
    type(silja_shopping_list), intent(inout) :: meteo_input_dyn_shopping_list, &
                                              & meteo_full_dyn_shopping_list, &
                                              & disp_dyn_shopping_list, disp_stat_shopping_list
    type(silja_shopping_list), intent(in) :: output_dyn_shopping_list, output_stat_shopping_list
    type(Tfield_buffer), pointer :: meteo_ptr, disp_buf_ptr, output_buf_ptr
    type(mini_market_of_stacks), pointer :: meteoMarketPtr, dispersionMarketPtr, outputMarketPtr
    logical, intent(inout) :: ifGotNewMeteoData



    CHARACTER (LEN=fnlen) :: command_string = ' '
    type(silja_time) :: reftime

    !All stuff needed for mid of next timestep
    reftime = now + timestep*0.5 + fu_meteo_time_shift(wdr)

    if (ifGotNewMeteoData) then
        command_string = 'Meteodata_acquisition'
        call msg(command_string)
        call start_count(chCounterNm = command_string)
        CALL fix_shopping_time_boundaries(meteo_input_dyn_shopping_list, & !pMeteo_input_dyn_shop_list, &
                                        & reftime, &
                                        & reftime) ! + simRules%timestep)
        CALL fill_meteo_market(meteoMarketPtr, wdr, &
                             & meteo_input_dyn_shopping_list, &  !pMeteo_input_dyn_shop_list, &
                             & timestep, &  !! needed only for direction
                             & iAccuracy, &
                             & ifGotNewMeteoData)
        IF (error) return
        call stop_count(chCounterNm = command_string)
     endif
            
      !-------------------------------------------------------------------------
      !
      ! Whether the new data are consumed or not, we should create or, at least, check
      ! the diagnostic quantities and set the pointers in all MiniMarkets.
      ! So far, we handle 
      ! (i) meteomarket using derived_field_quantities modules
      ! (ii) dispersion market using functions in this module
      !
     command_string = 'Meteodata_processing'
     call start_count(chCounterNm = command_string)
     
     call make_all_diagnostic_fields(meteoMarketPtr, meteo_ptr, &         ! meteo market & buffer
                                & dispersionMarketPtr, disp_buf_ptr, & ! disp market & buffer
                                & outputMarketPtr, output_buf_ptr, & ! output market & buffer
                                & meteo_full_dyn_shopping_list, &     ! meteo dynamic shopping list
                                & disp_dyn_shopping_list, &      ! dispersion dynamic shopping list
                                & output_dyn_shopping_list, &  ! full list of output dynamic fields
                                & disp_stat_shopping_list, &     ! dispersion static shopping list
                                & output_stat_shopping_list, &   ! fuopenll list of output static fields
                                & wdr, &                          ! weather data rules
                                & diagRules, &     ! diagnostic rules
                                & ifGotNewMeteoData, &            ! if new meteo data avaialble
                                & now, timestep, & ! now, step & weight of past-meteo time
                                & .True. )       ! if dispersion structures exist

     call stop_count(chCounterNm = command_string)
    end subroutine prepare_meteo

  !**********************************************************************************
    
  subroutine cut_mesh(mesh_to_modify, mesh_to_cover, mesh_new, ifVerbose)
    ! Find minimal sub-mesh of mesh_to_modify needed to provide interpolation
    ! to every point of mesh_to_cover
    !
    implicit none
    type(silja_grid), intent(in) :: mesh_to_modify, mesh_to_cover
    type(silja_grid), intent(out) :: mesh_new
    logical, intent(in) :: ifVerbose
    
    integer :: nx_gtm, ny_gtm, nx_gtc, ny_gtc
    integer :: ix, iy, ix_start_new, ix_end_new, iy_start_new, iy_end_new, nx_new, ny_new
    integer, dimension(:), pointer :: iWork
    integer, dimension(:,:), pointer :: covermask
    logical :: lTmp
    character(len=*), parameter :: sub_name = 'cut_mesh'


    if (fu_fails(defined(mesh_to_modify), 'gtm not defined', sub_name)) return
    if (fu_fails(defined(mesh_to_cover), 'gtc not defined', sub_name)) return
    
    call grid_dimensions(mesh_to_modify, nx_gtm, ny_gtm)

    iWork => fu_work_int_array(nx_gtm*ny_gtm)
    covermask(1:nx_gtm,1:ny_gtm) =>  iWork(1:nx_gtm*ny_gtm)
    if(ifVerbose) call msg("fu_if_mesh_intrpolatable verbose:")
    lTmp = fu_if_mesh_intrpolatable(mesh_to_cover, mesh_to_modify, covermask, ifVerbose)
    if (lTmp) then
      do  iy_start_new = 1, ny_gtm !Bottom
         if (any(covermask(:,iy_start_new) /= 0)  ) exit
      enddo
      do iy_end_new = ny_gtm , iy_start_new, -1 !Top
         if (any(covermask(:,iy_end_new) /= 0)  ) exit
      enddo

      if (fu_ifLonGlobal(mesh_to_cover)) then
        ix_start_new = 1  !Global grid -- no cut in X direction
        ix_end_new = nx_gtm
      else
        do  ix_start_new = 1, nx_gtm !Left
           if (any(covermask(ix_start_new,iy_start_new:iy_end_new) /= 0)  ) exit
        enddo
        do ix_end_new = nx_gtm , ix_start_new, -1 !Right
           if (any(covermask(ix_end_new,iy_start_new:iy_end_new) /= 0)  ) exit
        enddo
      endif

      mesh_new = mesh_to_modify
      if(ifVerbose)then
        call msg("Mesh Before cut")
        call report(mesh_new)
      endif
      call cut_grid_size(mesh_new, ix_start_new, iy_start_new, ix_end_new, iy_end_new)
      if(ifVerbose)then
        call msg("mesh_to_cover")
        call report(mesh_to_cover)
        call msg("Cover mask:")
        do iy = ny_gtm , 1, -1 !Top-bottom, so map looks normal
          call msg("",covermask(:,iy) > 0)
        enddo
        call msg("Mesh After cut")
        call report(mesh_new)
        call msg("stupidity check")
      endif

      ! Final check (Could fail due to numerics)
      lTmp = fu_if_mesh_intrpolatable(mesh_to_cover, mesh_new, covermask, ifVerbose)
      
      if (.not. lTmp .and. iy_end_new < ny_gtm) then !Try to extend one cell up
           call msg("Hmm... Not yet..  Trying to extend north...")
           mesh_new = mesh_to_modify
           iy_end_new = iy_end_new + 1
           call cut_grid_size(mesh_new, ix_start_new, iy_start_new, ix_end_new, iy_end_new)
           lTmp = fu_if_mesh_intrpolatable(mesh_to_cover, mesh_new, covermask, ifVerbose)
     endif

     if (.not. lTmp .and. ix_end_new < nx_gtm) then !Try to extend one cell right
           call msg("Hmm... Not yet..  Trying to extend east...")
           mesh_new = mesh_to_modify
           ix_end_new = ix_end_new + 1
           call cut_grid_size(mesh_new, ix_start_new, iy_start_new, ix_end_new, iy_end_new)
           lTmp = fu_if_mesh_intrpolatable(mesh_to_cover, mesh_new, covermask, ifVerbose)
     endif

     if (.not. lTmp .and. ix_start_new > 1) then !Try to extend one cell left
           call msg("Hmm... Not yet..  Trying to extend west...")
           mesh_new = mesh_to_modify
           ix_start_new = ix_start_new - 1
           call cut_grid_size(mesh_new, ix_start_new, iy_start_new, ix_end_new, iy_end_new)
           lTmp = fu_if_mesh_intrpolatable(mesh_to_cover, mesh_new, covermask, ifVerbose)
     endif

     if (.not. lTmp .and. iy_start_new > 1) then !Try to extend one cell down
           call msg("Hmm... Not yet..  Trying to extend south...")
           mesh_new = mesh_to_modify
           iy_start_new = iy_start_new - 1
           call cut_grid_size(mesh_new, ix_start_new, iy_start_new, ix_end_new, iy_end_new)
           lTmp = fu_if_mesh_intrpolatable(mesh_to_cover, mesh_new, covermask, ifVerbose)
     endif


    endif !Interpolatable??

    if (.not. lTmp .or. ifVerbose) then
       if (.not. mesh_to_cover == mesh_new) then
         !create covermask
         lTmp = fu_if_mesh_intrpolatable(mesh_to_cover, mesh_new, covermask, .True.)
          nx_new =  ix_end_new - ix_start_new -1 
          ny_new =  iy_end_new - iy_start_new -1 
          call msg("Cover mask of final cut attempt:")
          do iy = ny_new , 1, -1 !Top-bottom, so map looks normal
            call msg("",covermask(1:nx_new,iy) > 0)
          enddo
       else
        call msg("mesh_to_cover == mesh_new")
       endif
    endif

    if (.not. lTmp) then



       call msg("Failed to cut mesh: Non-interpolatable")
       call msg("mesh_to_modify")
       call report(mesh_to_modify)
       call msg("mesh_to_cover")
       call report(mesh_to_cover)
       !Call with full verbosity
       lTmp = fu_if_mesh_intrpolatable(mesh_to_cover, mesh_to_modify, covermask, .true.)
       call msg("Cover mask:")
       do iy = ny_gtm , 1, -1 !Top-bottom, so map looks normal
         call msg("",covermask(:,iy) > 0)
       enddo
       call set_error("Failed to cut mesh","cut_mesh")
    endif
    call free_work_array(iWork)
    

  end subroutine cut_mesh


  !*****************************************************************************

  subroutine tune_output_parameters(wdr, meteoVarLst, meteoVarLstST, PCld, em_source, &
                                  & OutDef, chemRules, dynRules, nlStdSetup, &
                                  & pDispersionMarket, met_buf, &
                                  & model_time_step)
    !
    ! After the input meteo files are analysed, we may have to tune the
    ! output definition: adjust the list of output variables, vertical structures, 
    ! etc.
    ! Also, here we explore the generic requests for the output of e.g. 
    ! concentration of source_inventory nuclides end others. Full list of available
    ! nuclides is in the cloud
    !
    implicit none

    ! Imported parameters
    type(silja_wdr), pointer :: wdr
    type(silja_shopping_list), intent(in) :: meteoVarLst, meteoVarLstST
    type(silam_pollution_cloud), pointer :: PCld
    type(silam_source), pointer :: em_source
    type(silam_output_definition), intent(inout), target :: OutDef
    type(silja_interval), intent(in) :: model_time_step
    type(Tchem_rules), intent(in) :: chemRules
    type(Tdynamics_rules), intent(in) :: dynRules
    type(Tsilam_namelist), pointer :: nlStdSetup
    type(mini_market_of_stacks), pointer :: pDispersionMarket
    type(Tfield_buffer), pointer :: met_buf

    ! Local variables
    integer :: i,j, al_status, nSrc, i2D,i3D,iWind2D,iWind3D, tmpFieldKind, iExtra, iTmp
    logical :: found, ifMultiLevel, corner_in_geographical_latlon, if_south_pole, &
             & ifDD_all_cumulative, ifDD_all_rate, ifWD_all_cumulative, ifWD_all_rate
    type(silja_stack), pointer :: stackPtr
    character(len=clen) :: chTmp
    type(ToutputVariables), pointer :: ov
    type(silam_species), dimension(:), pointer :: ptrSpecies

!    call msg("tune_output_parameters got meteoVarLstST:")
!    call report(meteoVarLstST)
!    call msg("tune_output_parameters got meteoVarLst:")
!    call report(meteoVarLst)

    !
    ! Stupidity test first. There was a lot of mess before with definitions. 
    ! Let's check that some basic stuff is OK
    !
    if(.not.defined(meteo_grid))then
      call set_error('Undefined meteo_grid','tune_output_parameters')
      return
    endif
    if(.not.defined(output_grid))then
      call set_error('Undefined output grid','tune_output_parameters')
      return
    endif
    if(.not.defined(meteo_vertical))then
      call set_error('Undefined meteo_vertical','tune_output_parameters')
      return
    endif
    if(.not.defined(output_vertical))then
      call set_error('No output vertical','tune_output_parameters')
      return
    endif

    ov => OutVars
    !
    ! Nullify the MassMapLinks and Lagrangian-mass map links
    !
    OutVars%MassMapLinks%iVerticalTreatment = int_missing
    OutVars%Lagr2MMLinks%iVerticalTreatment = int_missing

    !
    ! We have to check that the meteo and output verticals are compatible
    ! with each other, i.e. we have all necessary data for cross-interpolation
    !
    select case(fu_leveltype(output_vertical))
      case(constant_pressure, constant_height, constant_altitude, &
         & layer_btw_2_pressure, layer_btw_2_height, layer_btw_2_altitude)

        select case(fu_leveltype(meteo_vertical))
          case(constant_pressure, constant_height, sigma_level, hybrid, layer_btw_2_hybrid) ! OK
          case default
            call set_error('Incompatible meteo and output verticals', &
                         & 'tune_output_parameters')
            call report(meteo_vertical)
            call report(output_vertical)
        end select

      case(hybrid)
        !
        ! If output levels are hybrid ones - copy their coefficients from meteodata.
        ! If some coefs do not exist in meteodata - set error
        !
        if(fu_leveltype(meteo_vertical) /= hybrid)then
          call set_error('System vertical is not hybrid, while output one is', &
                       & 'tune_output_parameters')
          return
        end if
        do i =1, fu_nbrOfLevels(output_vertical)
          found = .false.
          do j=1, fu_NbrOfLevels(meteo_vertical)
            if(fu_cmp_levs_eq(fu_level(meteo_vertical,j), fu_level(output_vertical,i)))then
              call set_level(output_vertical,i,fu_set_hybrid_level( &
                              & fu_hybrid_level_number(fu_level(meteo_vertical,j)), &
                              & fu_hybrid_level_coeff_a(fu_level(meteo_vertical,j)), &
                              & fu_hybrid_level_coeff_b(fu_level(meteo_vertical,j))))
              found = .true.
              exit
            end if
          end do
          if(.not.found)then
            call set_error('Failed to find system level for this one:','tune_output_parameters')
            call report(output_vertical)
            return
          endif
        end do

      case default
        if(.not.defined(output_vertical))then
          call set_error('Unsupported output vertical type','tune_output_parameters')
          call report(output_vertical)
          return
        endif
    end select

    !------------------------------------------------------------
    !
    ! Check that all needed variables are in varLst, delete desirable ones
    ! if they are not in the list
    !
    if(allocated(OutDef%Rules%DispOutLst%ptrItem))then
      i=1
      do while(i <= size(OutDef%Rules%DispOutLst%ptrItem))

        ! May be, the list is over...
        !
        if(OutDef%Rules%DispOutLst%ptrItem(i)%quantity == int_missing) exit

!        call msg_test('')
        call msg_test('Checking output dispersion variable:' + &
                    & fu_quantity_short_string(OutDef%Rules%DispOutLst%ptrItem(i)%quantity) + ',' +&
                    & fu_str(OutDef%Rules%DispOutLst%ptrItem(i)%species))
        !
        ! If the dispersion quantity is requested but the model is turned off, check the request and 
        ! proceed apropriately. The rest will be handled by the expand_dispersion_output_list
        !
        if(.not. OutDef%Rules%ifRunDispersion)then
          if(OutDef%Rules%DispOutLst%ptrItem(i)%request == 2)then
            call set_error('SILAM dispersion output requested, but model is off', &
                         & 'tune_output_parameters')
            return
          else
            call msg_warning('Removing variable from output list:' + &
                      & fu_quantity_short_string(OutDef%Rules%DispOutLst%ptrItem(i)%quantity) + ',' + &
                      & fu_str(OutDef%Rules%DispOutLst%ptrItem(i)%species))
            call shrink_output_list(OutDef%Rules%DispOutLst, i)
            cycle
          endif
        endif
        i = i + 1
      end do  ! SILAM dispersion quantities
    endif  ! Allocated dispersion output item

    ! Now - meteorological quantities
    !
    OutDef%Rules%ifAverageMeteo = .false.
    i=1
!call msg("size(OutDef%Rules%MeteoOutLst%ptrItem)", size(OutDef%Rules%MeteoOutLst%ptrItem))
    if (allocated(OutDef%Rules%MeteoOutLst%ptrItem)) then
      do while (i <= size(OutDef%Rules%MeteoOutLst%ptrItem))

        ! May be, the list is over...
        !
        if(OutDef%Rules%MeteoOutLst%ptrItem(i)%quantity == int_missing) exit
        call msg_test('Checking meteo quantity:' + &
                    & fu_quantity_string(OutDef%Rules%MeteoOutLst%ptrItem(i)%quantity))
        !
        ! Non-SILAM quantity => it must be in the meteo shopping list
        ! If the quantity is not in full shopping list - check the request and 
        ! proceed appropriately
        !
        if((.not. fu_quantity_in_list(OutDef%Rules%MeteoOutLst%ptrItem(i)%quantity, meteoVarLst) ).and. &
         & (.not. fu_quantity_in_list(OutDef%Rules%MeteoOutLst%ptrItem(i)%quantity, meteoVarLstST)))then
          if(OutDef%Rules%MeteoOutLst%ptrItem(i)%request == 2)then
    call msg("tune_output_parameters got meteoVarLst:")
    call report(meteoVarLst)
    call msg("tune_output_parameters got meteoVarLst:")
    call report(meteoVarLstST)

            call set_error('Output variable:### "' + &
                         & fu_quantity_string(OutDef%Rules%MeteoOutLst%ptrItem(i)%quantity) + &
                         & '" ### is dropped from run',  &
                         & 'tune_output_parameters')
            return
          else  ! Delete the variable and shrunk the list
            call msg_warning('Removing the quantity from output list:' + &
                           & fu_quantity_string(OutDef%Rules%MeteoOutLst%ptrItem(i)%quantity))
            call shrink_output_list(OutDef%Rules%MeteoOutLst, i)
            cycle
          endif ! request <> 2
        end if  ! quantity is not in shopping list

        ! Check if we should collect meteo at every model timestep....
        if (OutDef%Rules%MeteoOutLst%ptrItem(i)%request > 0 .and. &
             & OutDef%Rules%MeteoOutLst%ptrItem(i)%AvType /= iAsIs ) then
             OutDef%Rules%ifAverageMeteo = .true.
        endif

        i = i + 1

      end do ! list of output meteo quantities
    endif !allocated


    !------------------------------------------------------------------
    !
    ! Expand the output dispersion list. Actually replaces 'SOURCE_INVENTORY'
    ! and 'FULL_INVENTORY' with exact lists of species.
    ! Also takes care of 3D mass-map output
    !
    if(OutDef%rules%ifRunDispersion)then
      !
      ! Finally, check whether ALL deposition fields are cumulative or ALL are rates.
      ! Note that species may or may not be explored, it does not matter here.
      !
      ifDD_all_cumulative = .true.
      ifDD_all_rate = .true.
      ifWD_all_cumulative = .true.
      ifWD_all_rate = .true.
      do i = 1, size(OutDef%Rules%DispOutLst%ptrItem)
        if(OutDef%Rules%DispOutLst%ptrItem(i)%quantity == drydep_flag)then
          ifDD_all_cumulative = ifDD_all_cumulative .and. &
                              & (OutDef%Rules%DispOutLst%ptrItem(i)%avType == iAsIs .or. &
                               & OutDef%Rules%DispOutLst%ptrItem(i)%avType == iCumulative)
          ifDD_all_rate = ifDD_all_rate .and. &
                        & (OutDef%Rules%DispOutLst%ptrItem(i)%avType == iInstant .or. &
                         & OutDef%Rules%DispOutLst%ptrItem(i)%avType == iAverage .or. &
                         & OutDef%Rules%DispOutLst%ptrItem(i)%avType == iMeanLastHrs)
        elseif(OutDef%Rules%DispOutLst%ptrItem(i)%quantity == wetdep_flag)then
          ifWD_all_cumulative = ifWD_all_cumulative .and. &
                              & (OutDef%Rules%DispOutLst%ptrItem(i)%avType == iAsIs .or. &
                               & OutDef%Rules%DispOutLst%ptrItem(i)%avType == iCumulative)
          ifWD_all_rate = ifWD_all_rate .and. &
                        & (OutDef%Rules%DispOutLst%ptrItem(i)%avType == iInstant .or. &
                         & OutDef%Rules%DispOutLst%ptrItem(i)%avType == iAverage .or. &
                         & OutDef%Rules%DispOutLst%ptrItem(i)%avType == iMeanLastHrs)
        endif
        if(ifDD_all_cumulative)then
          OutDef%Rules%ifDryDepCumulative = silja_true
        elseif(ifDD_all_rate)then
          OutDef%Rules%ifDryDepCumulative = silja_false
        else
          OutDef%Rules%ifDryDepCumulative = silja_undefined
        endif
        if(ifWD_all_cumulative)then
          OutDef%Rules%ifWetDepCumulative = silja_true
        elseif(ifWD_all_rate)then
          OutDef%Rules%ifWetDepCumulative = silja_false
        else
          OutDef%Rules%ifWetDepCumulative = silja_undefined
        endif
      end do  ! cycle over dispersion output list
      !
      ! First, initialise the Mass Maps for the 3D dispersion output from
      ! 3D mass maps and lpSets:
      ! concentration_flag, volume_mixing_ratio_flag, 
      ! advection_moment_X_flag, advection_moment_Y_flag, advection_moment_Z_flag
      ! emission_intensity_flag, optical_density_flag, volume_mixing_ratio_flag
      ! dry_deposition_flag, wet_deposition_flag
      ! Note that masses can be in Eulerian and Lagrangian environments
      !
      call init_mass_output(OutDef, PCld, em_source, nlStdSetup, chemRules, dynRules%simulation_type, &
                          & model_time_step)
      if(error)return
      !
      ! Now, expand the remaining pieces of dispersion output: quantities not belonging to mass maps
      !
      call expand_dispersion_output_list(OutDef,PCld, em_source, &
                                       & pDispersionMarket, nlStdSetup, dynRules)
      if(error) return
    endif

    !-----------------------------------------------------------------------------
    !
    ! Let's count the number of fields to be in meteo output 
    !
    nSrc = size(OutDef%Params%chSrcNm)
    if(error)return
    i2D = 0
    i3D = 0
    iWind2D = 0
    iWind3D = 0

    if (allocated(OutDef%Rules%MeteoOutLst%ptrItem)) then
      call msg('')
      call msg('Initializing the meteorological output:')

      do i=1,size(OutDef%Rules%MeteoOutLst%ptrItem)
        if(OutDef%Rules%MeteoOutLst%ptrItem(i)%quantity == int_missing)exit ! Done

        if(OutDef%Rules%MeteoOutLst%ptrItem(i)%avType == iAverage .or. &
         & OutDef%Rules%MeteoOutLst%ptrItem(i)%avType == iCumulative .or. &
         & OutDef%Rules%MeteoOutLst%ptrItem(i)%avType == iMeanLastHrs)then
          iExtra = 1
        else
          iExtra = 0
        endif
        !
        ! Here we have to check how actually the particular quantity is available:
        ! whether is is 2D or 3D. It can be done via meteo shopping list only
        !
        if( OutDef%Rules%MeteoOutLst%ptrItem(i)%iVerticalTreatment == lowest_level_flag .or. &
         & OutDef%Rules%MeteoOutLst%ptrItem(i)%iVerticalTreatment == integrate_column_flag)then
          OutDef%Rules%MeteoOutLst%ptrItem(i)%if3D = .false.
        elseif(fu_mlev_quantity_in_list(meteoVarLst, &
                                 & OutDef%Rules%MeteoOutLst%ptrItem(i)%quantity) == silja_true)then
          OutDef%Rules%MeteoOutLst%ptrItem(i)%if3D = .true.
        elseif(fu_mlev_quantity_in_list(meteoVarLst, &
                                 & OutDef%Rules%MeteoOutLst%ptrItem(i)%quantity) == silja_false)then
          OutDef%Rules%MeteoOutLst%ptrItem(i)%if3D = .false.
        else
          OutDef%Rules%MeteoOutLst%ptrItem(i)%if3D = &
                            & fu_multi_level_quantity(OutDef%Rules%MeteoOutLst%ptrItem(i)%quantity)
        endif

        if(OutDef%Rules%MeteoOutLst%ptrItem(i)%if3D)then
          call msg('3D quantity:' + fu_quantity_string(OutDef%Rules%MeteoOutLst%ptrItem(i)%quantity))
          i2D = i2D + nz_output * (iExtra + 1)
          i3D = i3D + 1 + iExtra
          if(fu_wind_quantity(OutDef%Rules%MeteoOutLst%ptrItem(i)%quantity)) then
            iWind2D = iWind2D + nz_output * (iExtra + 1)
            iWind3D = iWind3D + 1 + iExtra
          endif
        else
          if(OutDef%Rules%MeteoOutLst%ptrItem(i)%iVerticalTreatment == lowest_level_flag)then
            call msg('lowest-level quantity:' + &
                                  & fu_quantity_string(OutDef%Rules%MeteoOutLst%ptrItem(i)%quantity))
          elseif(OutDef%Rules%MeteoOutLst%ptrItem(i)%iVerticalTreatment == integrate_column_flag)then
            call msg('column-integrated quantity:' + &
                                  & fu_quantity_string(OutDef%Rules%MeteoOutLst%ptrItem(i)%quantity))
          else
            call msg('2-D quantity:' + fu_quantity_string(OutDef%Rules%MeteoOutLst%ptrItem(i)%quantity))
          endif
          i2D = i2D + 1 + iExtra
          if(fu_wind_quantity(OutDef%Rules%MeteoOutLst%ptrItem(i)%quantity)) &
                                                                     & iWind2D = iWind2D + 1 + iExtra
        endif

      end do  ! Cycle through the output vars
    endif ! allocated

    !
    ! and intialise the output meteo stack. For all fields - take Nbr+1 to be safe
    !
    if(i2D > 0)then

      allocate(OutVars%meteoStack, stat = al_status)
      if(fu_fails(al_status ==0, 'Failed to allocate output meteo stack','tune_output_parameters'))return
      call set_defined(OutVars%meteoStack,silja_false)  ! To avoid some funny coinsidence

      CALL init_stack(i2D + 1, &         ! Nbr of fields
                    & iWind2D+1, &       ! Nbr of windfields
                    & 'meteo_output', &  ! Stack name
                    & .false., &         ! single time
                    & .false., &         ! sinlge met_src
                    & i3D+1, &           ! Nbr of 3d fields
                    & iWind3D+1, &       ! Nbr of 3d windfields
                    & wdr, &             ! weather data rules
                    & OutVars%meteoStack)  ! Stack itself
      if(error)return

      !
      ! We need an extra tmp stack if the interpolation from the meteo_grid 
      ! to the output grid is needed. Otherwise, the data can be accumulated 
      ! directly to the output stacks without loss of efficiency for 
      ! interpolation
      !
      if(output_grid == meteo_grid)then
        OutVars%MeteoTmpStack => OutVars%meteoStack
      else
        allocate(OutVars%meteoTmpStack, stat = al_status)
        if(fu_fails(al_status == 0,'Failed to allocate output intermediate meteo stack', &
                                 & 'tune_output_parameters'))return
        call set_defined(OutVars%meteoTmpStack,silja_false)  ! To avoid funny coinsidence
        CALL init_stack(i2D + 1, &              ! Nbr of fields
                      & iWind2D, &              ! Nbr of windfields
                      & 'meteo_TMP_output', &  ! Stack name
                      & .false., &               ! single time
                      & .false., &               ! sinlge met_src
                      & i3D, &                  ! Nbr of 3d fields
                      & iWind3D, &              ! Nbr of 3d windfields
                      & wdr, &                  ! weather data rules
                      & OutVars%meteoTmpStack)  ! Stack itself
        IF (error) RETURN
      endif   ! output grid == meteo_grid

      !
      ! The output stack should be filled-in by empty fields but with proper
      ! field_id -s. We shall never add/remove fields from this stack,
      ! they will be just collected in-between the output moments and zeroed,
      ! if needed, at the moments of output.
      ! Use meteo_grid since the fields must be accumulated directly from meteo supermarket
      !
      call fill_stack_with_empty_fields(OutVars%MeteoTmpStack, &  ! stack to fill
                                      & OutDef%Rules%MeteoOutLst, &   ! list of vars
                                      & meteo_grid)
      if(error)return
!      call report(OutVars%MeteoTmpStack)

    else
      nullify(OutVars%MeteoTmpStack, OutVars%MeteoStack)
    endif  ! If any meteo vars for output

    !
    ! Initalization of the model output stack
    !
    if(OutDef%Rules%ifRunDispersion)then
      !
      ! Count the number of fields to be in model output 
      !
      i2D = 0
      i3D = 0
      call msg('')
      call msg('Initializing dispersion output. List size =', size(OutDef%Rules%DispOutLst%ptrItem))

      do i=1,size(OutDef%Rules%DispOutLst%ptrItem)
        
        if(OutDef%Rules%DispOutLst%ptrItem(i)%quantity == int_missing) exit ! List over

        if(OutDef%Rules%DispOutLst%ptrItem(i)%avType == iAverage .or. &
         & OutDef%Rules%DispOutLst%ptrItem(i)%avType == iCumulative .or. &
         & OutDef%Rules%DispOutLst%ptrItem(i)%avType == iMeanLastHrs)then
          iExtra = 1
        else
          iExtra = 0
        endif
        !
        ! We have to take into account that emission has different cocktail from the others
        !
        if(OutDef%Rules%DispOutLst%ptrItem(i)%quantity == emission_intensity_flag)then
          ptrSpecies => fu_species_emission(PCld)
        elseif(OutDef%Rules%DispOutLst%ptrItem(i)%quantity == optical_density_flag .or. &
             & OutDef%Rules%DispOutLst%ptrItem(i)%quantity == optical_density_flag)then
          ptrSpecies => fu_species_optical(PCld)
        else
          ptrSpecies => fu_species_transport(PCld)
        endif

        if(OutDef%Rules%DispOutLst%ptrItem(i)%iVerticalTreatment == lowest_level_flag .or. &
         & OutDef%Rules%DispOutLst%ptrItem(i)%iVerticalTreatment == integrate_column_flag)then
          OutDef%Rules%DispOutLst%ptrItem(i)%if3D = .false.
        else
          OutDef%Rules%DispOutLst%ptrItem(i)%if3D = &
                        & fu_multi_level_quantity(OutDef%Rules%DispOutLst%ptrItem(i)%quantity)
        endif
        if(OutDef%Rules%DispOutLst%ptrItem(i)%if3D) then
          !
          ! Just potentially 3D var, but let's reserve all levels anyway
          !
          i2D = i2D + nz_output * (iExtra + 1)
          i3D = i3D + 1 + iExtra
          call msg('3D quantity:' + &
                 & fu_quantity_string(OutDef%Rules%DispOutLst%ptrItem(i)%quantity) + ',' + &
                 & fu_species_output_name(OutDef%Rules%DispOutLst%ptrItem(i)%species))
        else ! 2D/3D

          i2D = i2D + 1 + iExtra
          if(OutDef%Rules%DispOutLst%ptrItem(i)%iVerticalTreatment == lowest_level_flag)then
            call msg('lowest-level quantity:' + &
                   & fu_quantity_string(OutDef%Rules%DispOutLst%ptrItem(i)%quantity) + ',' + &
                   & fu_species_output_name(OutDef%Rules%DispOutLst%ptrItem(i)%species))
          elseif(OutDef%Rules%DispOutLst%ptrItem(i)%iVerticalTreatment == integrate_column_flag)then
            call msg('column-integrated quantity:' + &
                   & fu_quantity_string(OutDef%Rules%DispOutLst%ptrItem(i)%quantity) + ',' + &
                   & fu_species_output_name(OutDef%Rules%DispOutLst%ptrItem(i)%species))
          else
            call msg('2D quantity:' + &
                   & fu_quantity_string(OutDef%Rules%DispOutLst%ptrItem(i)%quantity) + ',' + &
                   & fu_species_output_name(OutDef%Rules%DispOutLst%ptrItem(i)%species))
          endif
        endif  ! 2D/3D
      end do  ! Cycle through the output vars

      ! Initialization itself. There will be the only output dispersion stack
      ! but possibly nSrc temporary output stacks
      !
      allocate(OutVars%DispTmpStack(nSrc), stat=al_status)
      if(fu_fails(al_status ==0, 'Failed to allocate output intermediate dispersion stack', &
                               & 'tune_output_parameters'))return
      do i=1, nSrc
        call set_defined(OutVars%DispTmpStack(i),silja_false)  ! To avoid funny coinsidence
        stackPtr => OutVars%DispTmpStack(i)
        CALL init_stack(i2D + 1, &
                      & 0, &  ! No windfields
                      & 'dispersion_output_tmp', &
                      & .false., &
                      & .false., &
                      & i3D, &      ! 
                      & 0, &      ! No 3D wind fields
                      & wdr, &
                      & stackPtr)
        IF (error) RETURN
      enddo
      !
      ! Actual output stack is to be separate if grid interpolation is needed,
      ! otherwise it will be just pointing to the tmp stacks one-by-one for
      ! all emission sources.
      !
      if(dispersion_grid == output_grid)then
        OutVars%DispStack => OutVars%DispTmpStack(1)
      else
        allocate(OutVars%DispStack, stat=al_status)
        if(fu_fails(al_status ==0, 'Failed to allocate output dispersion stack','tune_output_parameters'))return
        call set_defined(OutVars%DispStack,silja_false)  ! To avoid funny coinsidence

        CALL init_stack(i2D + 1, &
                      & 0, &      ! No wind fields
                      & 'silam_dispersion_output', &
                      & .false., &
                      & .false., &
                      & i3D, &
                      & 0, &      ! No 3D wind fields
                      & wdr, &
                      & OutVars%DispStack)
        IF (error) RETURN
      endif

      !
      ! Similar to meteo output tmp stacks, all these stacks have to be filled
      ! with empty fields using the dispersion_grid as a template.
      !
      do i=1, nSrc
        stackPtr => OutVars%DispTmpStack(i)
        call fill_stack_with_empty_fields(stackPtr , &        ! stack to fill
                                        & OutDef%Rules%DispOutLst, &  ! list of vars
                                        & dispersion_grid)
        if(error)then
          call msg("************************************************************")
          call msg("iSource:"+fu_str(i))
          call msg("")
          call report(stackPtr)
          call set_error("Failed filling output stacks with place-holders",'tune_output_parameters')
          call msg("************************************************************")
          return
        endif
      end do

    else  ! If Run dispersion

      nullify(OutVars%DispStack,OutVars%DispTmpStack)

    endif  ! If Run dispersion

    OutVars%iNbrCollection = 0
    if(error)return

    !
    ! Finally, initialise trajectories
    !
    IF(OutDef%Rules%ifTrajectory)THEN
      call msg_test('Initalising trajectories...')

      !
      ! !!max steps should use  simRules%residenceInterval !FIXME
      iTmp =  nint(abs(fu_period_to_compute(wdr)/model_time_step))
      iTmp = min(iTmp, max_trajectory_steps)
      call init_trajectory_output(outDef%Params%tr_set, &
            &  iTmp, & !max_steps
            & fu_nbr_of_lagr_particles(PCld,.false.) / OutDef%Rules%particles_per_trajectory, &
            & OutDef%Rules%particles_per_trajectory, & !Stotre every N-th
            & fu_nbr_of_sources(PCld), & 
            & OutDef%Rules%outTemplate, OutDef%Rules%ini_time, &
            & OutDef%Params%chCaseNm, OutDef%Params%chSrcNm(:), &
                                  & met_buf, model_time_step)
      IF(error)return
    END IF

    call msg(''); call msg('========================================='); call msg('')
    call msg('FINAL OUTPUT CONFIGURATION:')
    call msg('')
    call report(OutDef)
    call msg(''); call msg('')


  CONTAINS


    !================================================================================

    subroutine fill_stack_with_empty_fields(stackPtr, VarLst, grid_to_use)
      !
      ! Fills-in the whole stack with empty fields. No re-arrangement
      ! Second task - set the target IDs to the first output time.
      ! So far, there is no difference between the internal and output fields in
      ! terms of forecast length, validity length, accumulation length, etc. - 
      ! they are all zero fields.
      ! A specific type of averaging - total over the whole period. That var will 
      ! not be written at every output time but rather at the very end of the run.
      ! Important: the fields that are stored into the stack do not contain anything,
      ! which has to be reflected in their IDs. 
      ! Target IDs should describe the required situation by the moment of the 
      ! first output, which can be either immediate or after some time.
      !
      implicit none

      ! Imported parameters
      !
      type(silja_stack), pointer :: stackPtr
      type(TOutputList), intent(inout) :: VarLst
      type(silja_grid), intent(in), optional :: grid_to_use

      ! Local variables
      !
      type(silja_field_id) :: idTmp
      type(silja_field_id), pointer :: idPtr
      type(silja_field), pointer :: fieldPtr
      integer :: iQ, iLev, iMetSrc, nMetSrcs, iSubst
      type(silja_grid) :: gridTmp
      type(meteo_data_source), dimension(2), parameter :: met_src = &
                                             & (/fmi_silam_src, silam_internal_src/)
      type(silja_time) :: tmpAnalysisTime 
      type(silja_interval) :: tmpForecastLen, tmpAccumulationLen, tmpValidityLen
      type(silja_level) :: levelTmp
      type(silam_species), dimension(:), pointer :: speciesPtr

      if(.not.defined(stackPtr))then
        call set_error('Undefined stack given','fill_stack_with_empty_fields')
        return
      endif

      do iQ = 1, size(VarLst%ptrItem)

        if(VarLst%ptrItem(iQ)%quantity == int_missing) exit
        !
        ! First output time must be the valid time of the fields. Accumulation /averaging
        ! periods are deterined appropriately. Start of the simulations is in ini_time,
        ! lastOutputTime and Rules%timestep determine the next one.
        ! If this is run with the shifted meteotime, have to take it into account
        !
        tmpAnalysisTime = OutDef%Rules%ini_time + fu_meteo_time_shift(wdr)  !OutVars%LastOutputTime
        tmpForecastLen = zero_interval
        tmpAccumulationLen = zero_interval
        tmpValidityLen = zero_interval

        select case(VarLst%ptrItem(iQ)%AvType)
          case(iAsIs, iInstant)  ! No changes / Instant, valid at next output time
            tmpFieldKind = forecast_flag
            nMetSrcs = 1
          case(iCumulative)  ! Accumulated since start of simulations, valid next output
            tmpFieldKind = accumulated_flag
            nMetSrcs = 2
          case(iAverage) ! Averaged over the timestep, valid by the next 
            tmpFieldKind = averaged_flag
            nMetSrcs = 2
          case(iMeanLastHrs) ! Same as above but averaged over some specific period
            tmpFieldKind = averaged_flag
            nMetSrcs = 2
          case(iTotalWholePeriod)  ! Accumulated over the run, valid at the end
            !
            ! Tricky type of accumulation - it is not time-resolved variable,
            ! so it cannot be put into time-resolving output files. In v3.2.1
            ! this type of averagin is used for the emission set only, which
            ! is handled in a totally separate way. So, be carefull. In fact, here 
            ! it comes only to standardize the io_initialization.
            !
            tmpFieldKind = accumulated_flag
            nMetSrcs = 1
          case(iTechnicalOriginalGrid, iTechnicalDispersionGrid)
            !
            ! That is the technical output of source term. No stack, no output during the run
            !
            cycle

          case default
            call set_error('Unknown averaging type','fill_stack_with_empty_fields')
            call msg('Unknown averaging type', VarLst%ptrItem(iQ)%AvType)
            return
        end select

        ! Depending on the quantity and temporary/final output stack, the grid
        ! can be either dispersion_grid or output grid.
        !
        if(present(grid_to_use))then
          gridTmp = grid_to_use
        else
          gridTmp = fu_output_tmp_grid(VarLst%ptrItem(iQ)%quantity, OutDef, PCld)
        endif

        !
        ! Chemical names etc depend on the type of cocktail: emission or transport.
        ! 
        if(VarLst%ptrItem(iQ)%quantity == emission_intensity_flag)then
          speciesPtr => fu_species_emission(PCld)
        elseif(VarLst%ptrItem(iQ)%quantity == optical_density_flag .or. &
             & VarLst%ptrItem(iQ)%quantity == optical_column_depth_flag)then
          speciesPtr => fu_species_optical(PCld)
        else
          speciesPtr => fu_species_transport(PCld)
        endif
        if(error)return

        !
        ! For averaged variables we have to initialise two fields with
        ! different met_src - the first one is with normal fmi_silam_src
        ! and the second one with silam_internal_src
        !
        do iMetSrc = 1, nMetSrcs
          !
          ! Depending on 2D / 3D type of quantity, we have to initialise 1 or nLevs
          ! number of fields. Again, later the quantity may appear 2D, so we shall
          ! over-reserve the memory. but we do not know this yet
          !
          if(VarLst%ptrItem(iQ)%if3D) then

            do iLev = 1, nz_output

              idTmp = fu_set_field_id(met_src(iMetSrc),&  ! met_src
                                    & VarLst%ptrItem(iQ)%quantity, &
                                    & tmpAnalysisTime, &
                                    & tmpForecastLen, &
                                    & gridTmp, &
                                    & fu_level(output_vertical,iLev),&
!                                    & fu_central_level_of_layer(fu_level(output_vertical,iLev)),&
                                    & tmpAccumulationLen, &
                                    & tmpValidityLen, &
                                    & tmpFieldKind, &
                                    & VarLst%ptrItem(iQ)%species)
              !
              ! We have nothing in the field and have to tell this explicitly
              !
              call set_validity_length(idTmp, interval_missing)
              !
              ! Now store the field to the stack
              !
              call put_field_to_stack (idTmp,&         ! id
                                     & stack = StackPtr, &     ! stack
                                     & ifMeteo_grid = .false., &      ! ifAdjust to meteo_grid
                                     & iUpdateType = create_if_absent, &      ! ifUpdateAllowed
                                     & ifRandomise = fu_if_randomise(wdr), &
                                     & storage_grid = gridTmp, &      ! grid
                                     & storage_area = area_missing, & 
                                     & iAccuracy = 5, &
                                     & iOutside = setZero)        ! out of grid value
!call msg('***:'+fu_quantity_short_string(fu_quantity(idTmp)) + ',' + fu_str(fu_species(idTmp)))
            end do
            !
            ! Reset the level in order to accept all levels in the target id
            !
            call set_level(idTmp, level_missing)

          else  
            !
            ! 2D variables
            !
            if(VarLst%ptrItem(iQ)%iVerticalTreatment == do_nothing_flag)then
              levelTmp = ground_level
            elseif(VarLst%ptrItem(iQ)%iVerticalTreatment == integrate_column_flag)then
              levelTmp = entire_atmosphere_integr_level
            elseif(VarLst%ptrItem(iQ)%iVerticalTreatment == lowest_level_flag)then
              levelTmp = lowest_atmosphere_level
            else
              call set_error('Unknown vertical treatment:' + fu_str(VarLst%ptrItem(iQ)%iVerticalTreatment), &
                           & 'tune_output_parameters')
              return
            endif
            idTmp = fu_set_field_id(met_src(iMetSrc),&
                                  & VarLst%ptrItem(iQ)%quantity, &
                                  & tmpAnalysisTime, &
                                  & tmpForecastLen, &
                                  & gridTmp, &
                                  & levelTmp, &
                                  & tmpAccumulationLen, &
                                  & tmpValidityLen, &
                                  & tmpFieldKind, &
                                  & VarLst%ptrItem(iQ)%species)
            if(error)return
            call put_field_to_stack (idTmp,&      ! id
                                   & stack = StackPtr, &  ! stack
                                   & ifMeteo_grid = .false., &   ! ifAdjust to system_grid
                                   & iUpdateType = create_if_absent, &      ! ifUpdateAllowed
                                   & ifRandomise = fu_if_randomise(wdr), &
                                   & storage_grid = gridTmp, &      ! grid
                                   & storage_area = area_missing, & 
                                   & iAccuracy  = 5, &
                                   & iOutside = setZero)        ! out of grid value

            !
            ! We have nothing in the field and have to tell this explicitly
            !
            call find_field_from_stack(stackPtr, idTmp, fieldPtr, found)
            if(.not.found)then
              call set_error(fu_connect_strings('Problem with stack:',fu_name(stackPtr)), &
                           & 'tune_output_parameters')
              return
            endif
            idPtr => fu_id(fieldPtr)
            call set_level(idPtr, level_missing)
            call set_validity_length(idPtr, interval_missing)

          endif  ! 2D/3D

        end do ! iMetSrc
        !
        ! Set the target id. Level for multi-level quantity is level_missing = acceptt all
        ! There is no target id for the internal field. 
        !
        VarLst%ptrItem(iQ)%targetId = idTmp
        call set_met_src(VarLst%ptrItem(iQ)%targetId, fmi_silam_src)

        select case(VarLst%ptrItem(iQ)%AvType)
          case(iAsIs, iInstant)  ! No changes / Instant, valid at next output time
            call set_valid_time(VarLst%ptrItem(iQ)%targetId, OutVars%LastOutputTime + OutDef%Rules%timeStep)
          
          case(iCumulative, iAverage)  ! Accumulated since start of simulations, average, valid next output
            call set_valid_time(VarLst%ptrItem(iQ)%targetId, OutVars%LastOutputTime + OutDef%Rules%timestep)
            call set_accumulation_length(VarLst%ptrItem(iQ)%targetId, &
                                       & fu_valid_time(VarLst%ptrItem(iQ)%targetId) - OutDef%Rules%ini_time)

          case(iMeanLastHrs) ! Same as above but averaged over some specific period
            call set_valid_time(VarLst%ptrItem(iQ)%targetId, OutVars%LastOutputTime + OutDef%Rules%timeStep)
            call set_accumulation_length(VarLst%ptrItem(iQ)%targetId, VarLst%ptrItem(iQ)%AvPeriod)

          case(iTotalWholePeriod)  ! Accumulated over the run, valid at the end
            !
            ! Tricky type of accumulation - it is not time-resolved variable,
            ! so it cannot be put into time-resolving output files. In v3.2.1
            ! this type of averagin is used for the emission set only, which
            ! is handled in a totally separate way. So, be carefull. In fact, here 
            ! it comes only to standardize the io_initialization.
            !
            call set_valid_time(VarLst%ptrItem(iQ)%targetId, OutVars%LastOutputTime + fu_period_to_compute(wdr))
            call set_accumulation_length(VarLst%ptrItem(iQ)%targetId, fu_period_to_compute(wdr))

          case default
            call set_error('Unknown averaging type','tune_output_parameters')
            call msg('Unknown averaging type', VarLst%ptrItem(iQ)%AvType)
            return
        end select

      end do  ! Cycle through the output vars

    end subroutine fill_stack_with_empty_fields

  end subroutine tune_output_parameters


  !*****************************************************************************

  subroutine set_emission_output(OutDef, em_source, PollutionCloud, PeriodToCompute, &
                               & iAccuracy, ifRandomise)
    !
    ! Prepares emission output. If the total emission is requested, writes
    ! it down with appropriate split. Reason is: total value cannot be put
    ! into the time-resolving files because it is just one field. So, a
    ! separate files would be needed. Otherwise, the emission field is
    ! the same as the particle_counter or other variables, and its dynamic 
    ! values are in the collect_dispersion_data.
    !
    ! Actually, we do not have to prepare anything for the emission output 
    ! so far, so let's check if total is needed and, if yes, perform the 
    ! output and remove it from the list
    !
    !
    implicit none

    ! Imported parameters
    type(silam_output_definition), pointer :: OutDef
    type(silam_source), pointer :: em_source
    type(silam_pollution_cloud), pointer :: PollutionCloud
    type(silja_interval), intent(in) :: PeriodToCompute
    integer, intent(in) :: iAccuracy
    logical, intent(in) :: ifRandomise

    ! Local variables
    integer :: iVar, iSrcId, iLev, iSp, i, iUnitGrib, iUnitGrads, iUnitNetcdf, nSpecies
    type(silam_species), dimension(:), pointer :: pSpecies
    type(silja_field_id) :: idTmp
    real, dimension(:), pointer :: fArPtr
    type(silam_sp) :: sp_grib, sp_grads
    !
    ! Some stupidity check
    !
    if(.not. (OutDef%Rules%ifGrads .or. OutDef%Rules%iGrib == int_missing)) then
      call msg_warning('Neither GrADS nor GRIB formats are active','set_emission_output')
      return
    endif
    !
    ! Cycle through the list of vars looking for the combination of 
    ! emission_intensity_flag and iTotalWholePeriod as averaging type
    !
    do iVar = 1, size(OutDef%Rules%DispOutLst%ptrItem)
     if(OutDef%Rules%DispOutLst%ptrItem(iVar)%quantity == emission_intensity_flag)then
      if(OutDef%Rules%DispOutLst%ptrItem(iVar)%AvType == iTotalWholePeriod) then
        !
        ! In principle, we can have plenty of sources with separate output files
        !
        do iSrcId = 1, fu_NbrOf_source_ids(em_source)
          !
          ! Open GRIB file, if output includes it
          !
          if(OutDef%Rules%iGrib /= int_missing)then
            !
            ! Get binary file name from templates. Temporarily, we can use
            ! grib_fNm from OutDef to store the filename
            !
            sp_grib%sp => fu_work_string()
            sp_grib%sp = fu_FNm(OutDef%Rules%outTemplate, &
                              & OutDef%Rules%ini_time, &  ! ini_time
                              & OutDef%Rules%ini_time, &  ! anal_time, here it's the same
                              & zero_interval, &          ! forecast length
                              & OutDef%Params%chCaseNm, OutDef%Params%chSrcNm(iSrcId)) + &  
                       & '_emis.grib'
            if(error)return

            CALL open_gribfile_o('', sp_grib%sp, iUnitGrib)   ! directory is not used here
                               
            if(error)return
          endif  ! ifGrib
          !
          ! Open GrADS file, if output includes it
          !
          if(OutDef%Rules%ifGrads)then
            !
            ! Get file name from tempaltes
            !
            sp_grads%sp => fu_work_string()
            sp_grads%sp = fu_FNm(OutDef%Rules%outTemplate, &
                               & OutDef%Rules%ini_time, &  ! ini_time
                               & OutDef%Rules%ini_time, &  ! anal_time - here it is the same
                               & zero_interval, &  ! forecast length
                               & OutDef%Params%chCaseNm, OutDef%Params%chSrcNm(iSrcId)) + &
                        & '_emis.grads'
            if(error)return
            !
            ! For the multi-file output we should store the template rather than a file
            ! name into the GrADS structure. 
            !
            iUnitGrads = open_gradsfile_o('', sp_grads%sp, output_grid)  ! directory is not used 
                                                               
            if(error)return
            call msg('Emission GrADS file unit=',iUnitGrads)

          endif  ! ifGrads
          !
          ! NetCDF
          !
          if(OutDef%Rules%iNetcdf /= int_missing)then
            call set_error('EMission output cannot go into NetCDF files','set_emission_output')
            call unset_error('set_emission_output')
!            sp_nc%sp => fu_work_string()
!            sp_nc%sp = fu_FNm(OutDef%Rules%outTemplate, &
!                            & OutDef%Rules%ini_time, &  ! ini_time
!                            & OutDef%Rules%ini_time, &  ! anal_time - here it is the same
!                            & zero_interval, &  ! forecast length
!                            & OutDef%Params%chCaseNm, OutDef%Params%chSrcNm(iSrcId)) + &
!                     & '_emis.nc'
!            if(error)return
!
!            iUnitNetcdf = open_netcdf_file_o( &   
!                                  & sp_nc%sp, &   !fname, 
!                                  & output_grid, &     ! output grid
!                                  & output_vertical, & ! vertical, 
!                                  & OutDef%Rules%ini_time, &  ! here is it a start time
!                                  & OutDef%Rules%MeteoOutLst, OutDef%Rules%DispOutLst, &
!  !                                & chTemplate,               ! optional
!                                & fMissingVal=real_missing)  ! fMissingVal) - optional
          endif
          !
          ! Actual writing the fields. If sources are mixed - they are to be summed up
          !
          call get_cloud_inventory(PollutionCloud, .false., species_emission, pSpecies, nSpecies)
          if(error)return

          fArPtr => fu_work_array()

          do iSp= 1, nSpecies

            do iLev=1, nz_output
              !
              ! Writing will happen in two steps. First, the source projects itself
              ! to the given map. Second, this map is dropped into the output file.
              !
              call msg('Level:', iLev)
              idTmp = fu_set_field_id(fmi_silam_src,&
                                    & emission_intensity_flag, &
                                    & OutDef%Rules%ini_time, &  ! analysis time
                                    & PeriodToCompute, & !zero_interval, &   ! forecast length
                                    & output_grid,&
                                    & fu_level(output_vertical,iLev),&
                                    & PeriodToCompute, & !zero_interval, & ! length_of_accumulation, optional
                                    & zero_interval, & ! length_of_validity, optional
                                    & accumulated_flag, & ! field_kind, optional
                                    & pSpecies(iSp))
!call report(idTmp)
              fArPtr(1:fs_output)=0.

              call source_2_map(em_source, fArPtr, idTmp, .false., iSrcId, iAccuracy, ifRandomise)
              if(error)return

              call msg('Emission sum:',sum(fArPtr(1:fs_output)))

              if(OutDef%Rules%iGrib /= int_missing)then
                call msg_warning('No GRIB code for emission so far','set_emission_output')
                CALL write_next_field_to_gribfile(OutDef%Rules%iGrib, &
                                                & iUnitGrib, &
                                                & idTmp, &
                                                & fArPtr(1:fs_output))
                if(error)return ! call unset_error('set_emission_output')
              endif

              if(OutDef%Rules%ifGrads)then
                CALL write_next_field_to_gradsfile(iUnitGrads, &
                                                 & idTmp, &
                                                 & fArPtr(1:fs_output))
                if(error)return ! call unset_error('set_emission_output')
              endif

!              if(OutDef%Rules%ifNetcdf)then
!                CALL write_next_field_to_netcdf_file(iUnitNetcdf, &
!                                                 & idTmp, &
!                                                 & fArPtr(1:fsOut))
!                if(error) call unset_error('set_emission_output')
!              endif


            end do  ! Cycle over output levels

          end do  ! substances in the source inventory

          call free_work_array(fArPtr)

          !
          ! Close the files 
          !
          if(OutDef%Rules%iGrib /= int_missing) call close_gribfile_o(iUnitGrib)
          if(OutDef%Rules%ifGrads) call close_gradsfile_o(iUnitGrads,"")
!          if(OutDef%Rules%ifNETCDF) call close_netcdf_file(iUnitNetcdf)
          !
          ! Add new file names to the list of written files
          !
          do i=1,size(OutDef%Params%filesWritten)
            if(OutDef%Params%filesWritten(i) == '')exit
          end do
          if(i >= size(OutDef%Params%filesWritten)-2)then
            call msg_warning('List of written files is full','set_emission_output')
          else
            if(OutDef%Rules%iGrib /= int_missing) then 
              OutDef%Params%filesWritten(i) = sp_grib%sp
              i=i+1
              call free_work_array(sp_grib%sp)
            endif
            if(OutDef%Rules%ifGrads) then
              OutDef%Params%filesWritten(i) = sp_grads%sp
              call free_work_array(sp_grads%sp)
            endif
          endif

        end do ! through sources
        !
        ! Output request fulfilled, so the variable can be removed
        !
        call shrink_output_list(OutDef%Rules%DispOutLst, iVar)

      elseif(OutDef%Rules%DispOutLst%ptrItem(iVar)%AvType == iTechnicalOriginalGrid .or. &
           & OutDef%Rules%DispOutLst%ptrItem(iVar)%AvType == iTechnicalDispersionGrid)then
        !
        ! Technical output is requested: store the sources as namelists
        ! The only selection here is the type of grid
        !
        call store_as_namelist(em_source, &
                             & fu_FNm(OutDef%Rules%outTemplate, &
                                   & OutDef%Rules%ini_time, &  ! ini_time
                                   & OutDef%Rules%ini_time, &  ! anal_time, here it's the same
                                   & zero_interval, &          ! forecast length
                                   & OutDef%Params%chCaseNm, ''), & ! case and source names
                             & OutDef%Rules%DispOutLst%ptrItem(iVar)%AvType == iTechnicalOriginalGrid)
        !
        ! Output request fulfilled, so the variable can be removed
        !
        call shrink_output_list(OutDef%Rules%DispOutLst, iVar)

      endif ! type of emission output

     endif ! emission is requested for the output

    end do ! cycle through the output variables

  end subroutine set_emission_output


  !*****************************************************************************

  subroutine collect_output(metBuf, dispBuf, outBuf, &
                          & PCld, now, OutDef, wdr, model_time_step, simulation_type, ifFirstStep, ifLastOutput) 
    !
    ! Collects the instant data and accumulates them into the output stack
    ! Useful for e.g., averaging of some fields.
    ! It does not take care of trajectories - they have to be stored just after emission
    !
    ! Warning! Meteo is for the future time step
    ! Thus the order it to:
    !  1. collect massses
    !  2. collect instant meteo/dispersion quantities
    !  3. output 
    !  4. Collect accumulated/averaged  meteo/dispersion etc....
    ! On first step all fields are AS-IS
    ! On last step do not start new output period
    implicit none

    ! Improted parameters
    type(Tfield_buffer), pointer :: metBuf, dispBuf, outBuf
    logical, intent(in) :: ifFirstStep ! dump everythong AS_IS to the output
                                       ! before collect ing proper averages 
    logical, intent(in) :: ifLastOutput ! put stuff to the outut without collecting
                                      ! averages... Shoulf be always false within dispersion loop
    type(silam_pollution_cloud), pointer :: PCld
    type(silja_time), intent(in) :: now
    type(silam_output_definition), intent(inout), target :: OutDef
    type(silja_wdr), intent(in), optional :: wdr
    type(silja_interval), intent(in) :: model_time_step
    integer, intent(in) :: simulation_type

    ! Local variables
    integer :: iMet, iLev, iSourceId, iQ
    real, dimension(:), pointer :: dataPtr
    type(silja_field_id) :: idTmp
!    type(silam_cocktail), dimension(:), pointer :: cockArPtr
    type(silja_stack), pointer :: stackPtr
    real, dimension(:), pointer :: fTmp
    type(silja_field), pointer :: FldPtr
    logical :: ifOutputTime, ifDumpOutputTime, ifRatesDumpOutputTime
    type(silam_output_definition), pointer :: OD
    type(ToutputVariables), pointer :: OV
    integer, save :: iOutputCounter=0
    character (len=fnlen) :: sp
    character(len=*), parameter :: sub_name = 'collect_output'
    type(silja_interval) :: meteo_time_shift

    !
    ! If meteo is shifted, take this into account.
    ! Meteo should be treated with the shift but dispersion without
    !
    meteo_time_shift = fu_meteo_time_shift(wdr)

    OD => OutDef
    OV => OutVars

    ifOutputTime = now == (OutVars%LastOutputTime + OD%Rules%timestep)

    ifDumpOutputTime = .false.  
    if(defined(OutVars%LastDumpOutputTime))then
      ifDumpOutputTime = now == (OutVars%LastDumpOutputTime + OD%Rules%dump_timestep)
    endif

    ifRatesDumpOutputTime = .false.  
    if(defined(OD%Rules%rates_dump_timestep))then
      ifRatesDumpOutputTime = now == (OutVars%LastRatesDumpOutputTime + OD%Rules%rates_dump_timestep)
    endif
    
    call msg_test('Collecting the output masses')
    !------------------------------------------------------------------
    !
    ! The next step - collect the masses, both in air and deposited, as well as optics, 
    ! for both Eulerian and Lagrangian runs. Note that mass maps are present always: deposition fields
    !
    call start_count('Collect_intermediate_mass_maps')

    call merge_mass_map_2_mass_map(OutVars%MassMapLinks, &
                                 & dispBuf, &          ! data buffer, also has time
                                 & OutVars%LastOutputTime + OutDef%Rules%timestep, & ! next output time
                                 & now, &
                                 & model_time_step, &     ! timestep in mass maps
                                 & ifFirststep, ifLastOutput)
    if(error)return

    if(fu_if_lagrangian_present(simulation_type)) &
        & call merge_lagrPartSet_2_mass_map(OutVars%Lagr2MMLinks, &      ! Structure to merge
                                          & outBuf, &           ! data buffer
                                          & OutVars%LastOutputTime + OutDef%Rules%timestep, & ! next output time
                                          & now, &
                                          & model_time_step)     ! timestep in mass maps
                                        !!! ifLastOutput  should be here

    call stop_count('Collect_intermediate_mass_maps')


    ! Force all AS_IS on the first step to have all outputs defined
    call collect_buffers(metBuf, dispBuf, OD, OV, now + meteo_time_shift, model_time_step, &
        & ifFirstStep, fu_if_randomise(wdr), ifLastOutput, ifOutputTime)

    OutVars%iNbrCollection = OutVars%iNbrCollection + 1
    !
    ! If this is the right time for the output - store the collected data to files
    !
    if(ifOutputTime)then
      !
      ! Write output
      !
      call msg_test('Writing the output')
      call start_count('write_output')
!call check_stack_fields_ranges(OutVars%MeteoTmpStack)
      call write_output(metBuf, PCld, now, OutDef, wdr, model_time_step)
      call stop_count('write_output')
      if (error) return
      call msg_test('Writing finished')
      !
      ! Collect and report total masses
      !
      if(fu_ifRunDispersion(OutDef) .and. mod(iOutputCounter,OutDef%Rules%cldRepInterv) == 0)then
        call collect_total_masses(PCld)
        call report_total_masses(PCld, OutVars%iNbrCollection, .true.) 
        call report_inout_mass_cld(pCld,outgoing)
        call report_inout_mass_cld(pCld,incoming)
      end if
      iOutputCounter = iOutputCounter + 1
      if(error)then
        call msg_warning('Problems with writing',sub_name)
        return
      endif

      !
      ! And the dump. But be careful with the name
      !
      if(ifDumpOutputTime)then

        if(fu_fails(fu_if_eulerian_present(simulation_type),'Dump is not defined for Lagrangian runs',sub_name))return
        if(fu_fails(fu_ifRunDispersion(OutDef),'No dump if no-dispersion run',sub_name))return
        call msg('Dumping the concentration mass map')
        
        do iSourceId = 1, size(OutDef%Params%chSrcNm)

          if(len_trim(OutDef%Params%chSrcNm(iSourceId)) < 1)exit  ! all done

          !
          ! Get binary file name from templates. 
          !
          sp = fu_FNm(OutDef%Rules%outTemplate, &
                       & OutDef%Rules%ini_time, &  ! ini_time
                       & OutDef%Rules%ini_time, &  ! anal_time, here it's the same
                       & now - OutDef%Rules%ini_time, &          ! forecast length
                       & OutDef%Params%chCaseNm, OutDef%Params%chSrcNm(iSourceId)) + '_' + &
                       & trim(fu_str(now,.false.)) + '_dump.grads'
          if(error)return

          !        call mass_map_to_grads_file(fu_concMM_ptr(PCld), iSourceId, sp%sp, now, 1.0)
          call many_mass_maps_to_grads_file((/fu_concMM_ptr(PCld), &
                                            & fu_advection_moment_X_MM_ptr(PCld), &
                                            & fu_advection_moment_Y_MM_ptr(PCld), &
                                            & fu_advection_moment_Z_MM_ptr(PCld)/), &
                                            & iSourceId, sp, now, 1.0)
          if(error)return

        end do  ! iSrc

        OutVars%LastDumpOutputTime = now

      endif  ! if dump output time

      if(ifRatesDumpOutputTime)then
        if (associated(fu_reactRateMM_ptr(PCld)))  then
         call msg('Dumping the reaction rates map')
        else
          call set_error("Rates massmap is not associated...", sub_name) 
          return
        endif
        
        do iSourceId = 1, size(OutDef%Params%chSrcNm)

          if(len_trim(OutDef%Params%chSrcNm(iSourceId)) < 1)exit  ! all done

          !
          ! Get binary file name from templates. 
          sp = fu_FNm(OutDef%Rules%outTemplate, &
                       & OutDef%Rules%ini_time, &  ! ini_time
                       & OutDef%Rules%ini_time, &  ! anal_time, here it's the same
                       & now - OutDef%Rules%ini_time, &          ! forecast length
                       & OutDef%Params%chCaseNm, OutDef%Params%chSrcNm(iSourceId)) + '_' + &
                       & trim(fu_str(now,.false.)) + '_ratesdump.grads'
          if(error)return

          call many_mass_maps_to_grads_file((/fu_reactRateMM_ptr(PCld)/),  iSourceId, sp, now, 1.0)
          if(error)return
        end do  ! iSrc
        OutVars%LastRatesDumpOutputTime = now
      endif  ! if dump output time

      if (ifLastOutput) return
      !
      ! Start new output period
      !
      call msg('Next output time:' + fu_str(OutVars%lastOutputTime + OutDef%Rules%timestep))
      !
      ! Prepare target IDs and fields in temporary stack for further
      ! output. Actions with individual fields depend on their types
      ! Note that meteorology is all taken with the meteo_time_shift - in both duffers
      !
      call start_new_output_period_lst(OutDef%Rules%MeteoOutLst, &
                                     & OutVars%LastOutputTime + OutDef%Rules%timestep + meteo_time_shift, &  ! next output time
                                     & OutDef%Rules%timestep)
      call start_new_output_period_lst(OutDef%Rules%DispOutLst, &
                                     & OutVars%LastOutputTime + OutDef%Rules%timestep + meteo_time_shift, &  ! next output time
                                     & OutDef%Rules%timestep)
      call start_new_output_period_massmap(OutVars%MassMapLinks, &
                                         & dispBuf, &          ! data buffer, also has time
                                         & OutVars%LastOutputTime + OutDef%Rules%timestep, & ! next output time
                                         & now)
      call start_new_output_period_lpSet(OutVars%Lagr2MMLinks, now)
      !
      ! Temporary stacks preparation includes nullifying of averaged fields.
      ! Accumulated fields will continue accumulation, while instant fields will
      ! anyway be overwritten at a proper time.
      !
      if (associated(OutVars%MeteoTmpStack)) then
        if(defined(OutVars%MeteoTmpStack))then
          stackPtr => OutVars%MeteoTmpStack
          call prepare_new_averaging_period(stackPtr)
        endif
      endif

      if(OutDef%Rules%ifRunDispersion)then
        do iQ =1, size(OutDef%Params%chSrcNm)
          stackPtr => OutVars%DispTmpStack(iQ)
          call prepare_new_averaging_period(stackPtr)
        end do
      endif

    endif  ! output time

!    call msg_test('Output Finished')

  end subroutine collect_output


  !***********************************************************************************

  subroutine collect_buffers( metBuf,  dispBuf, OD, OV, now, model_time_step, &
      & ifForceASIS, ifRandomise, ifLastOutput, ifCollectASIS)
    ! 
    ! Should be clled by collect_output before actual output
    !
    ! "AS_IS" values are colected  ifBeforeOut == .True.
    !
    ! All others are collected as "AS_IS" if 
    ! ifBeforeOut == .True. and ifFirststep == .True.
    ! to provide SOME reasonable values for the output
    ! 
    ! Averaged and cumulative values are collected
    ! if  ifBeforeOut == .False.
    !
    type(Tfield_buffer), pointer :: metBuf, dispBuf
    logical, intent(in) :: ifForceASIS, ifRandomise
    logical, intent(in) :: ifLastOutput !! Nothing so far
    logical, intent(in) :: ifCollectASIS !! If we are not going to output
                                         !! ifCollectASIS can be .False.
    type(silam_output_definition), pointer :: OD
    type(ToutputVariables), pointer :: OV
    type(silam_pollution_cloud), pointer :: PCld
    type(silja_time), intent(in) :: now
    type(silja_interval), intent(in) :: model_time_step


    integer :: iQ, iSourceId, avtype, nTmp
    type(silja_logical) :: slTmp
    type(silja_stack), pointer :: stackPtr

      
    !
    ! There are two stacks - meteo and model data. So, let's start from the meteo
    ! It does not depend on splitting of sources, so all meteofields
    ! must be stored to all output files despite the emission sources
    !
    if (allocated(OD%Rules%MeteoOutLst%ptrItem)) then
        ! No averaged meteo quantities  => Skip  the second collection 

        call msg_test('Collecting the output meteo')

        call start_count('collect_met_field_data')
        do iQ=1,size(OD%Rules%MeteoOutLst%ptrItem)

          if(OD%Rules%MeteoOutLst%ptrItem(iQ)%quantity == int_missing) exit

          if (ifForceASIS) then
            avtype = iAsIs
          else
            avtype = OD%Rules%MeteoOutLst%ptrItem(iQ)%AvType
          endif

          if (.not. ifCollectASIS) then
            if (avtype == iAsIs) cycle
          endif
          !
          ! Collect this variable
          !

!call msg_test("Collecting MET with avtype _"+fu_str(avtype)+"_:" &
!                  &+fu_quantity_string(OD%Rules%MeteoOutLst%ptrItem(iQ)%quantity))
!if(OD%Rules%MeteoOutLst%ptrItem(iQ)%quantity == windspeed_10m_flag)then
!  call msg('*')
!endif
          call collect_field_data(metBuf, &             ! Buffer to take the field from
                                & meteo_verticalPtr, &  ! meteorological vertical used in that buffer
                                & now, &                ! Time now
                                & model_time_step, & 
                                & OD%Rules%MeteoOutLst%ptrItem(iQ),&
                                & OV%MeteoTmpStack, &   ! Output stack to merge the field in
                                & avtype, & 
                                & OD, ifRandomise)  ! Overall set of rules and supplementary information
          if(error) exit  !MUST die on error here!
          !call unset_error('collect_buffers')
        end do  ! Scan through the meteo output list
        call stop_count('collect_met_field_data')
        if(error) return
!call msg('done')
!if(error)call unset_error('Let see it further')
    endif ! allocated

    !------------------------------------------------------------------
    !
    ! Now - collect the dispersion quantities, specifically for each source
    !
    call msg_test('Collecting the output dispersion')
    call start_count('Collect_other_dispersion_fields')

    if(associated(OD%Params%chSrcNm))then
      do iSourceId = 1, size(OD%Params%chSrcNm)
      
        if(allocated(OD%Rules%DispOutLst%ptrItem))then
          do iQ=1,size(OD%Rules%DispOutLst%ptrItem)
      
            if(OD%Rules%DispOutLst%ptrItem(iQ)%quantity == int_missing)exit

            if (ifForceASIS) then
              avtype = iAsIs
            else
              avtype = OD%Rules%DispOutLst%ptrItem(iQ)%AvType
            endif
            if (.not. ifCollectASIS) then
              if (avtype == iAsIs) cycle
            endif

    !call msg_test('Collect_buffers:' + fu_quantity_string(OD%Rules%DispOutLst%ptrItem(iQ)%quantity) )

            !
            ! If the quantity is tight to pollution cloud, its output is driven 
            ! by specialised routine collect_dispersion_data.
            ! Otherwise, we are dealing with internal model field stored in 
            ! dispersion_buffer. It can be treated exactly like meteorological one
            !
            slTmp = fu_if_cloud_mass_map_quantity(OD%Rules%DispOutLst%ptrItem(iQ)%quantity) 
            if(slTmp == silja_true)then

call msg('Pollution_cloud quantity cannot be collected here. SOURCE: ', iSourceId)
call msg(fu_quantity_string(OD%Rules%DispOutLst%ptrItem(iQ)%quantity) + ',' + &
                         & fu_str(OD%Rules%DispOutLst%ptrItem(iQ)%species))

!call msg('Collecting data for:' + &
!                     & fu_quantity_short_string(OD%Rules%DispOutLst%ptrItem(iQ)%quantity)))
!          call start_count('collect_disp_data')
!
!          call collect_dispersion_data(PCld, metBuf, dispBuf, OD, iQ, now, model_time_step, iSourceId, &
!                                     & simulation_type)
!          call stop_count('collect_disp_data')

            elseif(slTmp == silja_false)then

              stackPtr => OV%DispTmpStack(iSourceId)

!if(OD%Rules%DispOutLst%ptrItem(iQ)%quantity == 260010 .or. &
! & OD%Rules%DispOutLst%ptrItem(iQ)%quantity == 250132)then
!call msg('Collecting field data for:' + &
!       & fu_quantity_short_string(OD%Rules%DispOutLst%ptrItem(iQ)%quantity))
!endif
              call start_count('collect_field_data')
              call collect_field_data(dispBuf, &  ! Buffer to take the field from
                                    & dispersion_verticalPtr, & ! vertical structure used in the buffer
                                    & now, &         ! Time now
                                    & model_time_step, &
                                    & OD%Rules%DispOutLst%ptrItem(iQ),&
                                    & stackPtr, &  ! Output stack to merge the field in
                                    & avtype, & !
                                    & OD, ifRandomise)   ! OD: all rules and supplementary information
              call stop_count('collect_field_data')

            else
              call msg('Neither mass map nor field quantity',OD%Rules%DispOutLst%ptrItem(iQ)%quantity)
              call set_error('Neither particle- nor field-based quantity','collect_buffers')
              return
            endif
            !
            ! If we fail at this level, have to stop: if this variable was not excluded, model thinks
            ! that it can deliver. If it cannot, something went wrong.
            !
            if(error)then
              call set_error('Cannot collect variable:' + &
                                   & fu_quantity_string(OD%Rules%DispOutLst%ptrItem(iQ)%quantity), &
                             & 'collect_buffers')
              return
            endif

          end do ! Quantities
        endif  ! if associated quantities

      end do ! Emission source ids
    endif  ! associated source IDs

    call stop_count('Collect_other_dispersion_fields')
  end subroutine collect_buffers

  
  !**************************************************************************

  subroutine collect_field_data(datBuffer,  verticalUsed, now, model_time_step, &
                              & listItem, stackOut, avtype, OutDef, ifRandomise)
    !
    ! Collects the data from the given buffer (meteo or dispersion one)
    ! and stores them to the output stack
    !
    implicit none

    ! Imported parameters
    type(Tfield_buffer), intent(inout) :: datBuffer
    type(silam_vertical), intent(in) :: verticalUsed
    type(silja_time), intent(in) :: now ! Neede to assign proper validity to 3d field
    type(TOutputLstItem), target, intent(in) :: listItem
    type(silja_stack), intent(inout) :: stackOut
    integer, intent(in) :: avtype ! might differ from listItem%AvType
                                  ! for the first time step
    type(silam_output_definition), intent(in) :: OutDef
    type(silja_interval), intent(in) :: model_time_step
    logical, intent(in) :: ifRandomise

    ! Local variables
    real, dimension(:), pointer :: dataPtr
    integer :: indQ_arr, iLev, iTmp, i
    type(silja_field_id), target :: idTmp, localIdTarget
    type(silja_field_id), pointer :: idTargetPtr, idInputPtr
    real :: weight_past_out
    logical :: ifWorkArray

    !
    ! Take an index in the buffer
    !
call msg("collect_field_data for:" + fu_quantity_short_string(listItem%quantity) + &
                         & ',' + fu_str(listItem%species))
    if (associated(datBuffer%buffer_species)) then
         indQ_arr = fu_index_of_variable(listItem%quantity, datBuffer%buffer_quantities, &
                      & listItem%species, datBuffer%buffer_species)
    elseif (listItem%species == species_missing) then
                indQ_arr = fu_index(listItem%quantity, datBuffer%buffer_quantities)
    else
            call set_error("No species defined in the buffer", "collect_field_data")
            return
    endif
    !if (listItem%quantity == day_temperature_2m_acc_flag) call ooops("day_temperature_2m_acc_flag")
!    if (listItem%quantity == day_mean_temperature_2m_flag) then
!      call report(datBuffer%p2d(indQ_arr)%present%idPtr)
!      call ooops("day_mean_temperature_2m_flag in collect_field_data")
!    endif
   


    if(indQ_arr < 1 .or. indQ_arr > size(datBuffer%buffer_quantities)) then
      call msg('Cannot find the variable in buffer. Index received:', indQ_arr)
      call msg('variable_requested:' + fu_quantity_short_string(listItem%quantity) + &
             & ',' + fu_str(listItem%species))
      do iLev = 1, size(datBuffer%buffer_quantities)
        call msg('Variable available:' + fu_quantity_short_string(datBuffer%buffer_quantities(iLev)) + &
               & ',' + fu_str(datBuffer%buffer_species(iLev)))
      end do
      call set_error('Cannot find the variable in buffer','collect_field_data')
      return
    endif

    !
    ! Get a field pointer from the meteobuffer - 2D(x-y) or 4D(x-y-z-t)
    ! and drop the output to the stack
    !
    ifWorkArray = .false.
    select case(fu_dimension(datBuffer, indQ_arr))
      case(2)
        !
        ! ATTENTION. 
        ! For output of instant and AS_IS variables we cannot use the present field even if 
        ! it is created. The trouble is that the buffer%weight_past and time tag of the present field 
        ! point at the middle of the model time step, whereas output is at its end. 
        ! Conversely, for averaged and cumulative variables we MUST use the present field because
        ! the mid-step time tag, together with timestep-long validity, ensures correct accumulation
        ! and averaging at the end of the needed period.
        ! Also beware of single-time fields: their past and present are the same pointers, and
        ! future is also pointed at them
        !
        if((avtype == iAsIs .or. avtype == iInstant) .and. &           ! instant / as-is output
         & .not. associated(datBuffer%p2d(indQ_arr)%past%idPtr, &      ! not single-time
                          & datBuffer%p2d(indQ_arr)%present%idPtr))then
          !
          ! Instant / as-is output field, no averaging. Recompute the "true" present field keeping
          ! in mind that we come here only at the output time moment
          !
          weight_past_out = (datBuffer%time_future - now) / (datBuffer%time_future - datBuffer%time_past)
          if(weight_past_out < 1e-5)then                  ! future field fits
            dataPtr => datBuffer%p2d(indQ_arr)%future%ptr
            idInputPtr =>datBuffer%p2d(indQ_arr)%future%idPtr
          elseif(weight_past_out > 0.9999)then            ! past field fits
            dataPtr => datBuffer%p2d(indQ_arr)%past%ptr
            idInputPtr =>datBuffer%p2d(indQ_arr)%past%idPtr
          else
            ! Nothing to do, interpolate
            dataPtr => fu_work_array(fu_number_of_gridpoints(fu_grid(datBuffer%p2d(indQ_arr)%past%idPtr)))
            ifWorkArray = .true.
            call interpolate_fields_in_time(datBuffer%p2d(indQ_arr)%past%idPtr, &   ! fld_id_past, 
                                          & datBuffer%p2d(indQ_arr)%past%ptr, &     ! data_past
                                          & datBuffer%p2d(indQ_arr)%future%idPtr, & ! fld_id_future
                                          & datBuffer%p2d(indQ_arr)%future%ptr, &   ! data_future
                                          & weight_past_out, &
                                          & idTmp, dataPtr, now)             ! fld_id_res, data_res, now
            if(error)return
            idInputPtr => idTmp
          endif  ! weight_past_out is neither 0 nor 1
        else
          !
          ! Averaged or cumulative output variable, or single-time input. Use present field
          !
          !!! Note, that present actually points to the previous time step, since the output goes before any diagnostics!!
          dataPtr => datBuffer%p2d(indQ_arr)%present%ptr
          idInputPtr => datBuffer%p2d(indQ_arr)%present%idPtr
        endif  ! avtype == as-is or instant

        !
        ! 2d field. Should the type of averaging be iAsIs just send two identical ids - use 
        ! the input id instead of target id
        !
        if(avtype == iAsIs)then
          ! grid of idTargetPtr can be staggered... 
          localIdTarget = idInputPtr
          call set_grid(localIdTarget, fu_grid(listItem%targetId))
          idTargetPtr => localIdTarget
        else
          idTargetPtr => listItem%targetId
        endif
        !
        ! Having selected the input data and id, merge them to the output stack
        !
        call merge_vector_data_to_stack(idInputPtr, &
                                      & dataPtr, &
                                      & stackOut, &
                                      & idTargetPtr, &
                                      & ifRandomise)
        if(error)return

      case(4)
        !
        ! 4d field. Scan through the output levels and for each - make time and 
        ! vertical interpolation. Drop obtained 2d fields to the stack level-by-level
        !
        ifWorkArray = .true.
        dataPtr => fu_work_array(fu_number_of_gridpoints( &
                                              & fu_grid(datBuffer%p4d(indQ_arr)%past%p2d(1)%idPtr)))

!      call msg('Found 4D variable in buffer. Index received:', indQ_arr)
!      call msg('variable_requested:' + fu_quantity_short_string(listItem%quantity) + &
!             & ',' + fu_str(listItem%species))

        do iLev = 1, nz_output
          !
          ! Vertical and time interpolations are done in buffer
          !


          if(avtype == iAsIs)then
            call make_2d_from_4d_field (datBuffer, indQ_arr, &
                                      & verticalUsed, &
                                      & fu_level(output_vertical,iLev), &
                                      & listItem%iVerticalTreatment, &
                                      & now, dataPtr, idTmp)
            if(error)return
            ! grid of idTmp can be staggered... 
            localIdTarget = idTmp
            call set_grid(localIdTarget, fu_grid(listItem%targetId))
            idTargetPtr => localIdTarget
            !We do not care about validity here: just "AS_IS"
            !call msg("Taking AS IS")
          else
            !! Dirty hack: here we pretend that the field is for the past time step
            !! The hack introduced due to rearrangement of the dispersion loop
            !! Ideally, the data should be _collected_ for the output
            !! in SILAM they are _pushed_ this causes a lot of trouble, but cost of rewriting
            !! is  big
            call make_2d_from_4d_field (datBuffer, indQ_arr, &
                                      & verticalUsed, &
                                      & fu_level(output_vertical,iLev), &
                                      & listItem%iVerticalTreatment, &
                                      & now - model_time_step*0.5, dataPtr, idTmp)
            if(error)return
            idTargetPtr => listItem%targetId
            !Set valid during the current timestep
            !call msg("Setting valid period")
            if (.not. fu_accumulated(idTmp) )then
              if(fu_interval_positive(model_time_step))then
                call set_valid_time(idTmp, now - model_time_step)
                call set_validity_length(idTmp, model_time_step)
              else
                call set_valid_time(idTmp, now )
                call set_validity_length(idTmp, fu_abs(model_time_step))
              endif
            endif
          endif

          if(listItem%iVerticalTreatment == lowest_level_flag)then 
            call set_level(idTmp, lowest_atmosphere_level)
          endif

          !if (fu_quantity(idTmp) == dispersion_v_flag) then
          !  if (iLev == 1) then
          !    call msg("**************Merging 4D field to output stack****************")
          !    call msg("idTmp (Input)") 
          !    call report(idTmp)
          !    call msg("listItem%targetId") 
          !    call report(idTargetPtr)
          !  endif
          !end if
          call merge_vector_data_to_stack(idTmp, dataPtr, stackOut, idTargetPtr, ifRandomise)
          if(error)return

          if(listItem%iVerticalTreatment == lowest_level_flag) exit !Only one layer

        end do

      case default
        call set_error('Strange dimension of the field:','collect_field_data')
        call msg('Index of the field in buffer: ', indQ_arr)

    end select ! field dimension

    if(ifWorkArray) call free_work_array(dataPtr)

  end subroutine collect_field_data


  !*****************************************************************************

  subroutine write_output(metBuf, PCld, now, OutDef, wdr, model_time_step)
    !
    ! Transfers the collected data to the output stacks, if needed, then
    ! writes the output stacks to files and prepares the stacks for future
    ! use.
    ! Data copying happens if output grid != system_grid. Preparation of the 
    ! fields for future use depends on their type and includes re-setting
    ! the target id, and, if needed, current-id and the data field itself.
    !
    implicit none

    ! Imported parameters
    type(Tfield_buffer), intent(in) :: metBuf
    type(silam_pollution_cloud), intent(in) :: PCld
    type(silja_time), intent(in) :: now
    type(silam_output_definition), intent(inout), target :: OutDef
    type(silja_wdr), intent(in), optional :: wdr
    type(silja_interval), intent(in) :: model_time_step

    ! Local variables
    integer :: iSource, iFileManipulation, iFileManipTmp, iGrads, iNetcdf, iGrib
    type(silja_stack), pointer :: stackPtr
    logical, save :: ifFirstTime = .true.

    !----------------------------------
    !
    ! Opening files, if requested by the time series allocation.
    ! Number of files to open depends on the number of sources processed
    ! simultaneously
    !
    iFileManipulation = DoNothing
    iGrads = int_missing
    iNetcdf = int_missing
    iGrib = int_missing

    if(ifFirstTime) then  ! First-time output
      iFileManipulation = OpenFiles
      ifFirstTime = .false.
    else
      select case(OutDef%Rules%OutFilesArrangement)
        case(all_in_one)
          iFileManipulation = DoNothing
        case(hourly_new_file)
          if(fu_hour(now) /= fu_hour(OutVars%lastOutputTime)) &
                                                & iFileManipulation = SwitchBinary
        case(daily_new_file)
          if(fu_day(now) /= fu_day(OutVars%lastOutputTime)) &
                                                & iFileManipulation = SwitchBinary
        case(monthly_new_file)
          if(fu_mon(now) /= fu_mon(OutVars%lastOutputTime)) &
                                                & iFileManipulation = SwitchBinary
        case(yearly_new_file)
          if(fu_year(now) /= fu_year(OutVars%lastOutputTime)) &
                                                & iFileManipulation = SwitchBinary
        case default
          call set_error('Unknown time series arrangement','write_output')
          return
      end select
    endif

    !
    ! Output consists of two parts: meteostack and dispersion stack, which can be split to sources.
    !
    ! Copy all fields from temporary meteostack to the output one, with interpolation
    !
!call msg_test('Write output 1')
!call msg('Checking stack range before final output collection')
!call check_stack_fields_ranges(OutVars%MeteoStack)
!call check_stack_fields_ranges(OutVars%MeteoTmpStack)
!call msg('Finished')
    if (associated(OutVars%MeteoTmpStack)) then
      if(defined(OutVars%MeteoTmpStack))then
        if(.not. meteo_grid == output_grid)then
          call copy_stack_grid_interpolation(OutVars%MeteoTmpStack, &  !stackFrom
                                           & OutVars%MeteoStack, &  ! stackTo, 
                                           & output_grid, &   !gridNew
                                           & .false., &  ! No copy of internal fields
                                           & setMissVal, &
                                           & fu_if_randomise(wdr))  ! out of grid value
          if(error)return
        endif !meteo_grid == OutDef%Params%grid

  !      call arrange_fields_in_stack(OutVars%MeteoStack)
  !      if(error)return
      endif
    endif

!call msg('Checking stack range after final output collection')
!call check_stack_fields_ranges(OutVars%MeteoStack)
!call check_stack_fields_ranges(OutVars%MeteoTmpStack)
    !
    ! Now - fill-in output dispersion stack(s) with immediate writing
    !
!call msg_test('Write output 2')

    if(OutDef%Rules%ifRunDispersion)then
      !
      ! If sources are split - there will be several names in SrcNm
      ! If they are mixed - there will be just one name in SrcNm. 
      !
      do iSource = 1, size(OutDef%Params%chSrcNm)
        !
        ! Temporary stacks are in dispersion_grid, output stacks are in output_grid
        ! But if the grids are the same, the stack pointers coincide and thus
        ! no data copying needed
        !
        iFileManipTmp = iFileManipulation
        if(dispersion_grid == output_grid)then

          OutVars%DispStack => OutVars%DispTmpStack(iSource) ! No interpolation needed

        else
          !
          ! Copy all fields from temporary stack to the output one with interpolation
          !
          stackPtr => OutVars%DispTmpStack(iSource)
          call copy_stack_grid_interpolation(StackPtr, &  !stackFrom
                                           & OutVars%DispStack, &   ! stackTo, 
                                           & output_grid, &   !gridNew
                                           & .false., &  ! No copy of internal fields
                                           & setMissVal, &
                                           & fu_if_randomise(wdr))  ! out of grid value
          if(error)return

        endif !dispersion_grid == OutDef%Params%grid
        !
        !   Writing MeteoStack to file
        !
        if (associated(OutVars%MeteoStack))then
          if(defined(OutVars%MeteoStack))then

            IF(OutDef%Rules%ifGRADS) then
              CALL do_file_manip_grads(OutDef, now, iSource, iFileManipTmp)
              iGrads = OutDef%params%grads_funit(iSource)
            endif

            IF(OutDef%Rules%iNETCDF /= int_missing) then
              call do_file_manip_netcdf(OutDef, now, iSource, iFileManipTmp, OutDef%Rules%iNETCDF)
              iNetcdf = OutDef%Params%netcdf_funit(iSource)
            endif

            IF(OutDef%Rules%iGRIB /= int_missing) then
              CALL do_file_manip_grib(OutDef, now, iSource, iFileManipTmp)
              iGrib = OutDef%Params%grib_funit(iSource)
            endif
            if(error)return

  !call msg('meteo: write_stack_to_files')
            CALL write_stack_to_files(OutVars%MeteoStack, iSource, &
                                    & iGrads, iNetcdf, iGrib, OutDef%Rules%iGrib, now, &
                                    & OutDef%Rules%iOutTimesType == iRegular)
            if(error)return

            iFileManipTmp = DoNothing

          endif  ! defined output meteo stack
        endif

        !
        !  Writing DispStack to file
        !
        if(defined(OutVars%DispStack))then

          IF(OutDef%Rules%ifGRADS) then
            CALL do_file_manip_grads(OutDef, now, iSource, iFileManipTmp)
            iGrads = OutDef%params%grads_funit(iSource)
          endif

          IF(OutDef%Rules%iNETCDF /= int_missing) then
            call do_file_manip_netcdf(OutDef, now, iSource, iFileManipTmp, OutDef%Rules%iNETCDF)
            iNetcdf = OutDef%Params%netcdf_funit(iSource)
          endif

          IF(OutDef%Rules%iGRIB /= int_missing) then
            CALL do_file_manip_grib(OutDef, now, iSource, iFileManipTmp)
            iGrib = OutDef%Params%grib_funit(iSource)
          endif
          if(error)return

          CALL write_stack_to_files(OutVars%DispStack, iSource, &
                                  & iGrads, iNetcdf, iGrib, OutDef%Rules%iGrib, now, &
                                  & OutDef%Rules%iOutTimesType == iRegular)
          if(error)return

          iFileManipTmp = DoNothing

        endif  ! Defined dispersion stack
        !
        ! Writing the mass map into the file. It will decide inside what and whether to write
        !
!call msg('In write_output, calling:  file manipulation')
        IF(OutDef%Rules%ifGRADS) then
          CALL do_file_manip_grads(OutDef, now, iSource, iFileManipTmp)
          iGrads = OutDef%params%grads_funit(iSource)
        endif

        IF(OutDef%Rules%iNETCDF /= int_missing) then
          call do_file_manip_netcdf(OutDef, now, iSource, iFileManipTmp,OutDef%Rules%iNETCDF)
          iNetcdf = OutDef%Params%netcdf_funit(iSource)
        endif

        IF(OutDef%Rules%iGRIB /= int_missing) then
          CALL do_file_manip_grib(OutDef, now, iSource, iFileManipTmp)
          iGrib = OutDef%Params%grib_funit(iSource)
        endif
        if(error)return

!call msg('writing...')
        CALL write_mass_map_links_to_files(OutVars%MassMapLinks, iSource, &
                                         & output_grid, output_vertical, &
                                         & metBuf, now, OutDef%Rules%timestep, &
                                         & iGrads, iNetcdf, iGrib, OutDef%Rules%iGrib, &
                                         & OutDef%Rules%iOutTimesType == iRegular, OutDef%Rules%MMtrimfactor)
        if(error)return
        CALL write_lpSet_links_to_files(OutVars%Lagr2MMLinks, iSource, &
                                      & output_grid, output_vertical, &
                                      & metBuf, now, OutDef%Rules%timestep, &
                                      & iGrads, iNetcdf, iGrib, OutDef%Rules%iGrib, &
                                      & OutDef%Rules%iOutTimesType == iRegular, OutDef%Rules%MMtrimfactor)
        if(error)return

      end do ! cycle through sources

    else  ! IfRunDispersion
      !
      ! No SILAM dispersion run, so the output stack is made. Drop it to the files. 
      ! Remember - trajectory output and Ensemble output formats are not involved
      ! if there is no dispersion. So, we have to handle two formats - GRIB and GrADS
      !
      if(defined(OutVars%MeteoStack))then

        do iSource = 1, size(OutDef%Params%chSrcNm)
          IF(OutDef%Rules%ifGRADS)THEN !----------------- GrADS
            CALL do_file_manip_grads(OutDef, now, 1, iFileManipulation)
            iGrads = OutDef%params%grads_funit(iSource)
          END IF

          IF(OutDef%Rules%iNETCDF /= int_missing)THEN !----------------- NETCDF
!call msg('In write_output, calling:  stack_to_netcdf_file')
            call do_file_manip_netcdf(OutDef, now, 1, iFileManipulation, OutDef%Rules%iNETCDF)
            iNetcdf = OutDef%Params%netcdf_funit(iSource)
          END IF

          IF(OutDef%Rules%iGRIB /= int_missing)THEN !------------------ GRIB
            CALL do_file_manip_grib(OutDef, now, 1, iFileManipulation)
            iGrib = OutDef%Params%grib_funit(iSource)
          END IF
          if(error)return
          !
          ! Here we call writing for time "now" because the time tags follow the dispersion time
          ! Without meteo_time_shift. Inside the subroutine, time of the id is forced to this "now"
          !
          CALL write_stack_to_files(OutVars%MeteoStack, 1, &
                                  & iGrads, iNetcdf, iGrib, OutDef%Rules%iGrib, now, &
                                  & OutDef%Rules%iOutTimesType == iRegular)
          if(error)return
        enddo
      else
        call set_error('No data to write','write_output')
        return
      endif
    endif ! IfRunDispersion

    !
    ! Fix the current output time and determine the next one
    !
    OutVars%lastOutputTime = now
    call msg('')
    call msg('Starting the new output period')
    call msg('Last output time now: ' + fu_str(OutVars%lastOutputTime))

    if(OutDef%Rules%iOutTimesType /= iRegular)then  
      !
      ! if special occasion, such as start of the next month, the output timestep
      ! has to be recomputed
      !
      if(model_time_step > one_second)then
        OutDef%Rules%timestep = fu_next_special_time(now+one_minute, OutDef%Rules%iOutTimesType, &
                                                   & forwards, model_time_step) - now
      else
        OutDef%Rules%timestep = fu_next_special_time(now-one_minute, OutDef%Rules%iOutTimesType, &
                                                   & backwards, model_time_step*(-1.)) - now
      endif
      !
      ! And adjust the output time to have an integer number of steps...
      !
      OutDef%Rules%timestep = model_time_step * real(int(OutDef%Rules%timestep / model_time_step + 0.5))
      call msg('Nbr of days until new output:',fu_days(OutDef%Rules%timestep))
    endif
!call msg('Checking stack range after starting the new period')
!call check_stack_fields_ranges(OutVars%MeteoStack)
!call check_stack_fields_ranges(OutVars%MeteoTmpStack)

  end subroutine write_output


  !********************************************************************************************

  subroutine start_new_output_period_lst(OutLst, next_output_time, output_time_step)
    !
    ! When output data are written, the new output period starts. This routine
    ! re-sets target IDs for all quantities and calls nullifying actions for
    ! all temporary stacks
    ! It is assumed that last output time is already stored in proper place
    !
    implicit none

    ! Imported parameters
    type(TOutputList), intent(inout) :: OutLst
    type(silja_time), intent(in) :: next_output_time
    type(silja_interval), intent(in) :: output_time_step

    ! Local variables
    integer :: i

    !
    ! In the target IDs, there is always new valid time, plus new accumulation
    ! period for cumulative fields
    !
    ! Meteo stack
    !
    if (allocated(OutLst%ptrItem)) then
      do i=1,size(OutLst%ptrItem)

        if(OutLst%ptrItem(i)%quantity == int_missing)exit
!        call msg("setting Next output time :" + fu_str(next_output_time) + ': for:'+fu_quantity_string(OutLst%ptrItem(i)%quantity))
        call set_valid_time(OutLst%ptrItem(i)%targetId, next_output_time) ! Next output time

        select case(OutLst%ptrItem(i)%AvType)
          case(iAsIs, iInstant, iTotalWholePeriod)
          case(iCumulative) 
            call set_accumulation_length(OutLst%ptrItem(i)%targetId, &
                            & fu_accumulation_length(OutLst%ptrItem(i)%targetId) + output_time_step)
          case(iAverage)
            call set_accumulation_length(OutLst%ptrItem(i)%targetId, output_time_step)
          case(iMeanLastHrs) 
              call set_accumulation_length(OutLst%ptrItem(i)%targetId, OutLst%ptrItem(i)%AvPeriod)
          case default
            call set_error('Unknown averaging type','start_new_output_period')
            return
        end select

!        call msg('new targetID')
!        call report(OutLst%ptrItem(i)%targetId)
      end do
    endif !allocated 

  end subroutine start_new_output_period_lst


  !*****************************************************************************

  subroutine finalise_output(PCld, OutDef, wdr, emSrc)
    !
    ! Writes  ctl files and closes all open output binaries
    ! Dumps remaining trajectories
    ! For Eulerian advection, also reports the out-of-grid transport
    !
    implicit none

    ! Improted parameters
    type(silam_pollution_cloud), intent(in) :: PCld
    type(silam_output_definition), intent(inout) :: OutDef
    type(silja_wdr), intent(in) :: wdr
    TYPE(silam_source), INTENT(in) :: emSrc

    ! Local variables
    integer :: iSource

    !
    ! Write trajectories
    !
    if(OutDef%Rules%ifTrajectory) call  finalize_trajectory_output(OutDef%params%tr_set)

    !
    ! GRIB and GrADS files
    !
    do iSource = 1, size(OutDef%Params%chSrcNm)
      
      IF(OutDef%Rules%ifGRADS)THEN !----------------- GrADS
        if(OutDef%Params%grads_funit(iSource) /= int_missing)then
          if(fu_sec(OutDef%Rules%timestep) < 0)then
            call invert_grads_binary(OutDef%Params%grads_funit(iSource), 2) ! single file so far
          endif
          call msg('Closing GRADS output file')
          if(len_trim(OutDef%Rules%chFixedNameTemplate) > 0)then
            call close_gradsfile_o(OutDef%Params%grads_funit(iSource), &
                                 & OutDef%Rules%chFixedNameTemplate)
          else
            call close_gradsfile_o(OutDef%Params%grads_funit(iSource),"")
          endif
          IF (error) RETURN
        endif
      END IF

      IF(OutDef%Rules%iNETCDF /= int_missing)THEN !------------------ NETCDF
        if(OutDef%Params%netcdf_funit(iSource) /= int_missing)then
          call msg('Closing NETCDF output file')
          call write_ctl_file_4_netcdf(OutDef%Params%netcdf_funit(iSource))
          call close_netcdf_file(OutDef%Params%netcdf_funit(iSource))
          IF (error) RETURN
        endif
      END IF

      IF(OutDef%Rules%iGRIB /= int_missing)THEN !------------------ GRIB
        if(OutDef%Params%grib_funit(iSource) /= int_missing)then
          call msg('Closing GRIB output file')
          call close_gribfile_o(OutDef%Params%grib_funit(iSource))
          IF (error) RETURN
        endif
      END IF

    end do ! nSrc
    
    if(fu_ifRunDispersion(OutDef))then
      !
      ! If there was a source dump, close it
      !
      call close_emission_dump(emSrc, OutDef%params%chCaseNm, OutDef%rules%outTemplate, &
                                                            & fu_emisMM_ptr(pCld), &
                                                            & fu_emission_moment_X_MM_ptr(pCld), &
                                                            & fu_emission_moment_Y_MM_ptr(pCld), &
                                                            & fu_emission_moment_Z_MM_ptr(pCld))
      !
      ! Report detailed out-of-grid transport
      !
      call report_inout_mass_cld(pCld,outgoing)
      call report_inout_mass_cld(pCld,incoming)
    endif

  end subroutine finalise_output


  !*****************************************************************

  function fu_output_tmp_grid(quantity, OutDef, cloud)result(grid)
    !
    ! Different quantities live in different grids during the simulations
    ! In particular, all meteo quantities are in meteo_grid or Arakawa-shifted 
    ! meteo_grid, all Lagrangian dispersion quantities can be projected to any grid, 
    ! while all field dispersion quantities, including depositions are in dispersion_grid.
    !
    ! The task of this function is to select a grid that would be suitable for the specific
    ! quantity. Evidently, it will depend on the quantity and type of advection routine
    !
    implicit none

    ! Return value
    type(silja_grid) :: grid
    
    ! Imported parameter
    integer, intent(in) :: quantity
    TYPE(silam_output_definition), INTENT(in) :: OutDef
    type(silam_pollution_cloud), pointer :: cloud

    if(.not. fu_silam_dispersion_quantity(quantity)) then
      grid = meteo_grid  ! Can be Arakawa-shifted, but here we do not know it
      return
    endif

    select case(quantity)
      case(particle_counter_flag, areas_of_risk_flag)
        grid = output_grid

      case(concentration_flag, drydep_flag, wetdep_flag, concentration_2m_flag,&
         & emission_intensity_flag)
        grid = dispersion_grid

      case default
        
        if(fu_if_cloud_mass_map_quantity(quantity) == silja_true)then 
          ! Indeed error of unknown function
          call msg('')
          call msg('Quantity name:' + fu_quantity_string(quantity) + ', and number:', quantity)
          call set_error('Unknown quantity','fu_output_tmp_grid')
        endif

        grid = meteo_grid  ! Never return missing grid !

    end select

  end function fu_output_tmp_grid


  ! ***************************************************************


  SUBROUTINE do_file_manip_grib(outDef, TimeValid, iSource, iFileManipulation)
    ! 
    ! Drops the whole stack to grib file. Actually it just goes through the
    ! stack and writes the fields to GRIB one-by-one.  
    ! Another trick:
    ! No need to go through the 3d fields - they just contain the pointers to 
    ! 2d ones, which are anyway included into 2d field lists.
    ! GRIB file has a nice feature: an order of the fields does not matter.
    !
    IMPLICIT NONE

    ! Imported parameters with intent IN:
    TYPE(silam_output_definition), INTENT(inout) :: outDef

    ! Imported parameters with intent INOUT or POINTER:
    type(silja_time), intent(in) :: TimeValid
    integer, intent(in) :: iSource, iFileManipulation

    ! Local declarations:
    INTEGER :: i, j, grib_return
    INTEGER :: iLayer, iFieldCount = int_missing
    type(meteo_data_source) :: mdsTmp

    ! Actions with the files: continue to write to currently open, take new binary
    ! or totally start new writing structures
    !
    select case (iFileManipulation)

      case (DoNothing)

      case (OpenFiles) ! Totally new files, including ctls
        !
        ! Close opened files first. Should not be called, in fact, but...
        ! 
        if(OutDef%Params%grib_funit(iSource) /= int_missing)then ! Close existing file
          call close_gribfile_o(OutDef%Params%grib_funit(iSource))
        endif
        !
        ! Get binary file name from templates
        !
        OutDef%Params%grib_fNm(iSource) = &
                        & fu_FNm(OutDef%Rules%outTemplate, &
                               & OutDef%Rules%ini_time, &  ! ini_time
                               & OutDef%Rules%ini_time, &  ! anal_time - here it is the same
                               & (TimeValid-OutDef%Rules%ini_time), & ! forecast length
                               & OutDef%Params%chCaseNm, OutDef%Params%chSrcNm(iSource)) + '.grib'
        if(error)return
        !
        ! For the multi-file output we should store the template rather than a file
        ! name into the GrADS structure. 
        !
        if(OutDef%Rules%OutFilesArrangement == all_in_one)then
          
          CALL open_gribfile_o('', &  ! directory - not used here
                             & OutDef%Params%grib_fNm(iSource), &
                             & OutDef%Params%grib_funit(iSource)) 
        else

          CALL open_gribfile_o('', &  ! directory - not used here
                 & OutDef%Params%grib_fNm(iSource), &
                 & OutDef%Params%grib_funit(iSource), &
                         & fu_pure_grads_template(OutDef%Rules%OutTemplate, &
                                                & OutDef%Params%chCaseNm, &
                                                    & OutDef%Params%chSrcNm(iSource)) + '.grib')
        endif
        if(error)return
        !
        ! Add a new file name to the list of files
        !
        do i=1,size(OutDef%Params%filesWritten)
          if(OutDef%Params%filesWritten(i) == '')exit
        end do
        if(i > size(OutDef%Params%filesWritten))then
          call msg_warning('List of written files is full','do_file_manip_grib')
        else
          OutDef%Params%filesWritten(i) = OutDef%Params%grib_fNm(iSource)
        endif

      case(SwitchBinary) ! Keep ctl but change binary
        !
        ! Get new file name from tempaltes
        !
        OutDef%Params%grib_fNm(iSource) = &
                        & fu_FNm(OutDef%Rules%outTemplate, &
                               & OutDef%Rules%ini_time, &  ! ini_time
                               & OutDef%Rules%ini_time, &  ! anal_time - here it is the same
                               & (timeValid-OutDef%Rules%ini_time), &  ! forecast length
                               & OutDef%Params%chCaseNm, OutDef%Params%chSrcNm(iSource)) + '.grib'
        if(error)return
        !
        ! Close existing GrADS binary and open a new one - in one routine
        ! ctl file must not be written or closed. May be, some rewinding of the 
        ! structure indices...
        !
        call switch_grib_binary_o('', &                   ! directory - not used 
                                & OutDef%Params%grib_fNm(iSource), &
                                & OutDef%Params%grib_funit(iSource))
        if(error)return
        !
        ! Add a new file name to the list of files
        !
        i = count(len_trim(OutDef%Params%filesWritten(:)) == 0) + 1
        if(i > size(OutDef%Params%filesWritten))then
          call msg_warning('List of written files is full','do_file_manip_grib')
        else
          OutDef%Params%filesWritten(i) = OutDef%Params%grib_fNm(iSource)
        endif
        
      case default
        call set_error('Unknown file manipulation ','do_file_manip_grib')
    end select

  END SUBROUTINE do_file_manip_grib


  ! ***************************************************************


  SUBROUTINE do_file_manip_grads(outDef, timeValid, iSource, iFileManipulation)

    ! Description: 
    ! Drops the whole stack to grads file. Actually it just goes through the
    ! stack and writes the fields one-by-one. The file structure is indeed
    ! quite strict, but here no precautions are done. Simply if the field 
    ! is rejected by the routine it is skipped.
    ! Another trick:
    ! No need to go through the 3d fields - they just contain the pointers to 
    ! 2d ones, which are anyway included into 2d field lists.
    !
    IMPLICIT NONE

    ! Imported parameters with intent IN:
    TYPE(silam_output_definition), INTENT(inout) :: outDef
    type(silja_time), intent(in) :: timeValid
    integer, intent(in) :: iSource, iFileManipulation

    ! Local declarations:
    INTEGER :: i, grib_return
    INTEGER :: iLayer, iFieldCount = int_missing, iInvert=0
    !
    ! If needed - open new file for the given source or just switch to another 
    ! binary file keeping the ctl file the same
    !
    select case (iFileManipulation)

      case (DoNothing)

      case (OpenFiles) ! Totally new files, including ctls
        !
        ! Close opened files first
        ! 
        if(OutDef%Params%grads_funit(iSource) /= int_missing)then ! Close existing file
          call close_gradsfile_o(OutDef%Params%grads_funit(iSource),"")
        endif
        !
        ! Get file name from tempaltes
        !
        OutDef%Params%grads_fNm(iSource) = fu_FNm(OutDef%Rules%outTemplate, &
                                                & OutDef%Rules%ini_time, &  ! ini_time
                                                & OutDef%Rules%ini_time, &  ! anal_time
                                                & (timeValid-OutDef%Rules%ini_time), &  ! forecast len
                                                & OutDef%Params%chCaseNm, &
                                                & OutDef%Params%chSrcNm(iSource)) + &
                                         & '.grads'
        if(error)return

        !
        ! For the multi-file output we should store the template rather than a file
        ! name into the GrADS structure. 
        !
        ! The Grads buffering is now switched on whenever mpi-io is on. The defeault
        ! buffering (in grads-io) is 50 fields, but can be set to more with the
        ! GRADSIO_BUF_SIZE environment variable, or disabled by giving it a value < 1.
        if(OutDef%Rules%OutFilesArrangement == all_in_one)then
          
          OutDef%Params%grads_funit(iSource) = open_gradsfile_o('', &  ! directory not used 
                                                              & OutDef%Params%grads_fNm(iSource), &
                                                              & output_grid, ifMPIIO=smpi_use_mpiio_grads, &
                                                              & ifBuffered=smpi_use_mpiio_grads)
        else

          OutDef%Params%grads_funit(iSource) = open_gradsfile_o('', &         ! directory not used 
                                        & OutDef%Params%grads_fNm(iSource), &
                                        & output_grid, &
                                        & fu_pure_grads_template(OutDef%Rules%OutTemplate, &
                                                               & OutDef%Params%chCaseNm, &
                                                               & OutDef%Params%chSrcNm(iSource)) + '.grads', &
                                        & ifMPIIO=smpi_use_mpiio_grads, ifBuffered=smpi_use_mpiio_grads)
        endif
        if(error)return
        !
        ! Add a new file name to the list of files
        !
        do i=1,size(OutDef%Params%filesWritten)
          if(OutDef%Params%filesWritten(i) == '')exit
        end do
        if(i > size(OutDef%Params%filesWritten))then
          call msg_warning('List of written files is full','do_file_manip_grads')
        else
          OutDef%Params%filesWritten(i) = OutDef%Params%grads_fNm(iSource)
        endif

      case(SwitchBinary)
        !
        ! Get new file name from tempaltes
        !
        OutDef%Params%grads_fNm(iSource) = fu_FNm(OutDef%Rules%outTemplate, &
                                                & OutDef%Rules%ini_time, &  ! ini_time
                                                & OutDef%Rules%ini_time, &  ! anal_time
                                                & (timeValid-OutDef%Rules%ini_time), &  ! forecast len
                                                & OutDef%Params%chCaseNm, &
                                                & OutDef%Params%chSrcNm(iSource)) + &
                                         & '.grads'
        if(error)return
        !
        ! Close existing GrADS binary and open a new one - in one routine
        ! ctl file must not be written or closed. May be, some rewinding of the 
        ! structure indices...
        !
        if(fu_sec(OutDef%Rules%timestep) < 0)iInvert = 1
        call switch_grads_binary_o(OutDef%Params%grads_funit(iSource), &
                                 & '', &                   ! directory - not used 
                                 & OutDef%Params%grads_fNm(iSource), iInvert, timeValid)
        if(error)return
        !
        ! Add a new file name to the list of files
        !
        i = count(len_trim(OutDef%Params%filesWritten(:)) == 0) + 1
        if(i > size(OutDef%Params%filesWritten))then
          call msg_warning('List of written files is full','do_file_manip_grads')
        else
          OutDef%Params%filesWritten(i) = OutDef%Params%grads_fNm(iSource)
        endif
        
      case default
        call msg('Unknown file manipulation:', iFileManipulation)
        call set_error('Unknown file manipulation ','do_file_manip_grads')
    end select

  END SUBROUTINE do_file_manip_grads


  ! *****************************************************************************

  subroutine do_file_manip_netcdf(outDef, timeValid, iSource, iFileManipulation, ncversion)
    !
    ! Drops the whole stack to netcdf file. 
    !
    IMPLICIT NONE

    ! Imported parameters with intent IN:
    TYPE(silam_output_definition), INTENT(inout) :: outDef
    type(silja_time), intent(in) :: timeValid
    integer, intent(in) :: iSource, iFileManipulation, ncversion

    ! Local declarations
    integer :: i
    character(len=4),dimension(3:4) :: chTypes = (/'.nc ','.nc4'/)
    character(len=6) :: tasksuff !Topology suffix for filename
    integer :: my_xc, my_yc, size_x, size_y

    ! Actions with the files: continue to write to currently open, take new binary
    ! or totally start new writing structures
    !
    select case (iFileManipulation)

      case (DoNothing)

      case (OpenFiles, SwitchBinary) ! Totally new files, including ctls
        !
        ! Close opened files first. Should not be called, in fact, but...
        ! 
        if(OutDef%Params%netcdf_funit(iSource) /= int_missing)then ! Close existing file

          call close_netcdf_file(OutDef%Params%netcdf_funit(iSource))
          
        endif
        !
        ! Get binary file name from templates
        !
        call smpi_get_process_xy_topology(my_xc, my_yc, size_x, size_y)
        if (size_x*size_y == 1 .or. smpi_use_mpiio_netcdf) then
                
                tasksuff = "      "  !No need for separate filenames on subdomains
        else
                 WRITE(tasksuff, fmt = '(A, I2.2, A, I2.2)') '_',my_xc,'_', my_yc
        endif

        OutDef%Params%netcdf_fNm(iSource) = &
                        & fu_FNm(OutDef%Rules%outTemplate, &
                               & OutDef%Rules%ini_time, &  ! ini_time
                               & OutDef%Rules%ini_time, &  ! anal_time - here it is the same
                               & (TimeValid-OutDef%Rules%ini_time), & ! forecast length
                               & OutDef%Params%chCaseNm, OutDef%Params%chSrcNm(iSource))+trim(tasksuff) + chTypes(ncversion)

        if(error)return
        !
        ! For the multi-file output we should store the template rather than a file
        ! name into the GrADS structure. 
        !
        call msg("Calling open_netcdf_file_o for name:"+ OutDef%Params%netcdf_fNm(iSource))

        OutDef%params%netcdf_funit(iSource) = open_netcdf_file_o( &    
                                & OutDef%Params%netcdf_fNm(iSource), &   !fname, 
                                & output_grid, &     ! output grid
                                & output_vertical, & ! vertical, 
                                & timeValid, &             ! here it is a start time
                                & (/OutDef%Rules%MeteoOutLst, &     ! Meteo variables to write
                                  & OutDef%Rules%DispOutLst, &      ! Dispersion stack variables
                                  & OutDef%Rules%MassMapOutLst/), & ! Mass map variables
                                & fu_pure_grads_template(OutDef%Rules%OutTemplate, &
                                                        & OutDef%Params%chCaseNm, &
                                                        & OutDef%Params%chSrcNm(iSource)) + trim(tasksuff)+chTypes(ncversion), &
                                & OutDef%Rules%OutFilesArrangement == all_in_one, &
                                & ncversion, smpi_use_mpiio_netcdf, &
                                & real_missing)  ! fMissingVal
        if(error)return

        !
        ! Add a new file name to the list of files
        !
        do i=1,size(OutDef%Params%filesWritten)
          if(OutDef%Params%filesWritten(i) == '')exit
        end do
        if(i > size(OutDef%Params%filesWritten))then
          call msg_warning('List of written files is full','do_file_manip_netcdf')
        else
          OutDef%Params%filesWritten(i) = OutDef%Params%netcdf_fNm(iSource)
        endif

      case default
        call set_error('Unknown file manipulation ','do_file_manip_netcdf')
    end select

  END SUBROUTINE do_file_manip_netcdf



  !*****************************************************************************

  subroutine init_mass_output(OutDef, PCld, em_source, nlStdSetup, chemRules, &
                            & simulation_type, model_time_step)
    !
    ! Looks through the dispersion output list, selects the 3D instant Mass-Map based
    ! variables, i.e. concentration_flag, 
    ! advection_moment_X_flag, advection_moment_Y_flag, advection_moment_Z_flag
    ! emission_intensity_flag, optical_density_flag, volume_mixing_ratio_flag
    ! - and stores them into the dedicated mass-map intermediate output structure.
    ! The reason for separate treatment is simply to avoid reordering the indices
    ! at every model time step. With high-resolution calculations and hourly output
    ! the gain can be very large.
    !
    implicit none

    ! Imported parameters
    type(silam_output_definition), intent(inout) :: OutDef
    type(silam_source), pointer :: em_source
    type(silam_pollution_cloud), pointer :: PCld
    type(Tsilam_namelist), pointer :: nlStdSetup
    type(Tchem_rules), intent(in) :: chemRules
    integer, intent(in) :: simulation_type
    type(silja_interval), intent(in) :: model_time_step

    ! Local variables
    integer :: iTmp, iVarIni, nSpeciesOut, iSpecies, qTmp, quantity
    type(Tmass_map), pointer :: pMassMapIn
    type(silam_vertical), pointer :: pVert
    type(Tlagrange_particles_set), pointer :: pLPSetIn
    type(silam_species), dimension(:), pointer :: pSpeciesIn, pSpeciesOut
    logical :: ifQuantityMassMapAdded, ifQuantityMassMapAddedTmp
    type(Toptical_density_rules), pointer :: pOpticalRules
    type(TchemicalRunSetup), pointer :: pChemRunSetup

    ! We shall do it the following way: the dispersion output list will be scanned
    ! searching for the corresponding quantity. Should this be found, we shall take it,
    ! create the intermedate mass map, set the link and remove the list item. 
    ! Note that there can be both Lagrangian and Eulerian links for mass, concentration,
    ! and mixing ratio output.
    !
    ! Split the arrays: copy known mass-related quantities into their own list,
    ! leaving the rest intact in the dispersion output list
    !
!    allocate(OutDef%Rules%MassMapOutLst%ptrItem(size(OutDef%Rules%DispOutLst%ptrItem)), stat=iTmp)
    call expand_output_list(OutDef%Rules%MassMapOutLst, size(OutDef%Rules%DispOutLst%ptrItem ))

    nSpeciesOut = 0
    nullify(pSpeciesOut)
    OutDef%Rules%MassMapOutLst%ptrItem(1)%quantity = int_missing
    iTmp = 1
    iVarIni = 1
    !
    ! If it is from mass map, the below function returns silja_true
    !
    do while(iVarIni <= size(OutDef%Rules%DispOutLst%ptrItem))
      qTmp = OutDef%Rules%DispOutLst%ptrItem(iVarIni)%quantity
      if (qTmp <= 0) then 
        iVarIni = iVarIni + 1  !! fu_if_cloud_mass_map_quantity barks if fed with int_missing
      elseif(fu_true(fu_if_cloud_mass_map_quantity(OutDef%Rules%DispOutLst%ptrItem(iVarIni)%quantity)))then
          OutDef%Rules%MassMapOutLst%ptrItem(iTmp) = OutDef%Rules%DispOutLst%ptrItem(iVarIni)
          iTmp = iTmp + 1
          call shrink_output_list(OutDef%Rules%DispOutLst, iVarIni)
          if(iTmp <= size(OutDef%Rules%MassMapOutLst%ptrItem)) &
                                      & OutDef%Rules%MassMapOutLst%ptrItem(iTmp)%quantity = int_missing
      else
          iVarIni = iVarIni + 1
      endif
    enddo 
    !
    ! The main cycle over the MassMap output list
    !
    do iVarIni = 1, size(OutDef%Rules%MassMapOutLst%ptrItem)

      if(OutDef%Rules%MassMapOutLst%ptrItem(iVarIni)%quantity == int_missing)exit ! cycle ended
      if(OutDef%Rules%MassMapOutLst%ptrItem(iVarIni)%quantity < 0)cycle    ! already treated
      ifQuantityMassMapAdded = .false.
      !
      ! There can be several list items with the same quantity but with different species.
      ! We shall collect them into a single intermediate mass map, of course.
      ! To avoid too much complexity, we shall check and force the same averaging and vertical
      ! treatment for all species of one quantity
      !
      quantity = OutDef%Rules%MassMapOutLst%ptrItem(iVarIni)%quantity
      do iTmp = iVarIni, size(OutDef%Rules%MassMapOutLst%ptrItem)
        if(OutDef%Rules%MassMapOutLst%ptrItem(iTmp)%quantity == int_missing)exit  ! all done
        if(OutDef%Rules%MassMapOutLst%ptrItem(iTmp)%quantity == quantity)then
          if(OutDef%Rules%MassMapOutLst%ptrItem(iTmp)%iVerticalTreatment /= &
                        & OutDef%Rules%MassMapOutLst%ptrItem(iVarIni)%iVerticalTreatment)then
            call set_error('multiple calls with different vertical treatment of:' + &
                         & fu_quantity_string(quantity),'init_mass_output')
            return
          endif
          if(OutDef%Rules%MassMapOutLst%ptrItem(iTmp)%AvType /= &
                        & OutDef%Rules%MassMapOutLst%ptrItem(iVarIni)%AvType)then
            call set_error('multiple calls with different averaging of:' + &
                         & fu_quantity_string(quantity),'init_mass_output')
            return
          endif
          !
          ! Add the species for this item
          !
          call explore_item_species(OutDef%Rules%MassMapOutLst, iTmp, &   ! the list and item to explore
                                  & pCld, pSpeciesOut, nSpeciesOut, nlStdSetup)
          if(error)then 
            call set_error('Failed species for:' + fu_quantity_string(quantity), 'init_mass_output')
            return
          endif
          OutDef%Rules%MassMapOutLst%ptrItem(iTmp)%quantity = -1  ! mark as used
        endif ! repeated quantity call
      end do  ! search for repeating quantities

      if(nSpeciesOut == 0)then
        call msg_warning('No species for intermediate mass map for:' + fu_quantity_string(quantity),'init_mass_output')
        cycle  ! skip the quantity
      else
        call msg('Considering intermediate mass map for:' + fu_quantity_string(quantity))
      endif
      !
      ! Select the parameters for adding the intermediate mass map
      !
      if(quantity == concentration_flag .or. quantity == volume_mixing_ratio_flag)then
        !
        ! This is not trivial: concentrations can be transported, short-living, aerosols,
        ! can be in mass map or lagrangian particle set. Have to handle all carefully.
        !
        ! Add Eulerian links
        !
        if(fu_if_eulerian_present(simulation_type))then
          pMassMapIn => fu_concMM_ptr(PCld)
          pVert  => pMassMapIn%vertTemplate
          call add_intermediate_mass_map(fu_concMM_ptr(PCld), &            ! Transport species
                                       & fu_species_transport(PCld), &
                                       & quantity, &
                                       & pSpeciesOut, &
                                       & OutDef%Rules%MassMapOutLst%ptrItem(iVarIni), &
                                       & OutDef%Rules%ini_time, &
                                       & OutDef%Rules%timestep, &
                                       & model_time_step, &
                                       & ifQuantityMassMapAdded)
          call add_intermediate_mass_map(fu_aerosolMM_ptr(PCld) , &         ! Aerosol species
                                       & fu_species_aerosol(PCld), &
                                       & quantity, &
                                       & pSpeciesOut, &
                                       & OutDef%Rules%MassMapOutLst%ptrItem(iVarIni), &
                                       & OutDef%Rules%ini_time, &
                                       & OutDef%Rules%timestep, &
                                       & model_time_step, &
                                       & ifQuantityMassMapAddedTmp)
          ifQuantityMassMapAdded = ifQuantityMassMapAdded .or. ifQuantityMassMapAddedTmp
          call add_intermediate_mass_map(fu_shortlivedMM_ptr(PCld), &      ! Short-living species
                                       & fu_species_short_lived(PCld), &
                                       & quantity, &
                                       & pSpeciesOut, &
                                       & OutDef%Rules%MassMapOutLst%ptrItem(iVarIni), &
                                       & OutDef%Rules%ini_time, &
                                       & OutDef%Rules%timestep, &
                                       & model_time_step, &
                                       & ifQuantityMassMapAddedTmp)
          ifQuantityMassMapAdded = ifQuantityMassMapAdded .or. ifQuantityMassMapAddedTmp
          if(fu_fails(ifQuantityMassMapAdded,'Failed adding mass map for MassMapIn:' + &
                                       & fu_quantity_string(quantity),'init_mass_output'))return
        endif  ! Eulerian link needed
        !
        ! Add Lagrangian links
        !
        if(fu_if_lagrangian_present(simulation_type))then
          pLPSetIn => fu_lpSet_ptr(PCld)
          pVert  => pLPSetIn%verticalTemplate
          call add_intermed_mass_map_4_lpSet(pLPSetIn, 1, & !pLPSetIn%lpMassTrn, &
                                           & fu_species_transport(PCld), &
                                           & quantity, &
                                           & pSpeciesOut, &
                                           & OutDef%Rules%MassMapOutLst%ptrItem(iVarIni), &
                                           & output_grid, output_vertical, &
                                           & OutDef%Rules%ini_time, &
                                           & OutDef%Rules%timestep, &
                                           & model_time_step, &
                                           & ifQuantityMassMapAdded)
          call add_intermed_mass_map_4_lpSet(pLPSetIn, 2, &  !pLPSetIn%lpMassAer, &
                                           & fu_species_aerosol(PCld), &
                                           & quantity, &
                                           & pSpeciesOut, &
                                           & OutDef%Rules%MassMapOutLst%ptrItem(iVarIni), &
                                           & output_grid, output_vertical, &
                                           & OutDef%Rules%ini_time, &
                                           & OutDef%Rules%timestep, &
                                           & model_time_step, &
                                           & ifQuantityMassMapAddedTmp)
          ifQuantityMassMapAdded = ifQuantityMassMapAdded .or. ifQuantityMassMapAddedTmp
          call add_intermed_mass_map_4_lpSet(pLPSetIn, 3, &    ! pLPSetIn%lpMassSL, &
                                           & fu_species_short_lived(PCld), &
                                           & quantity, &
                                           & pSpeciesOut, &
                                           & OutDef%Rules%MassMapOutLst%ptrItem(iVarIni), &
                                           & output_grid, output_vertical, &
                                           & OutDef%Rules%ini_time, &
                                           & OutDef%Rules%timestep, &
                                           & model_time_step, &
                                           & ifQuantityMassMapAddedTmp)
          ifQuantityMassMapAdded = ifQuantityMassMapAdded .or. ifQuantityMassMapAddedTmp
          if(fu_fails(ifQuantityMassMapAdded,'Failed adding mass map for lpSetIn:' + &
                                           & fu_quantity_string(quantity),'init_mass_output'))return
        endif  ! lagrangian present
      else
        !
        ! Simple mass-map-only quantities
        !
        select case(quantity)

          case(concentration_2m_flag) !! Eulerian so far...
            call set_cnc2m_map(PCld)
            pMassMapIn => fu_concentration2mMM_ptr(PCld)
            pSpeciesIn => fu_species_transport(PCld)

          case(drydep_flag)
            pMassMapIn => fu_drydepMM_ptr(PCld)
            pSpeciesIn => fu_species_transport(PCld)

          case(wetdep_flag)
            pMassMapIn => fu_wetdepMM_ptr(PCld)
            pSpeciesIn => fu_species_transport(PCld)

          case(advection_moment_X_flag)
            pMassMapIn => fu_advection_moment_X_MM_ptr(PCld)
            pSpeciesIn => fu_species_transport(PCld)

          case(advection_moment_Y_flag)
            pMassMapIn => fu_advection_moment_Y_MM_ptr(PCld)
            pSpeciesIn => fu_species_transport(PCld)

          case(advection_moment_Z_flag)
            pMassMapIn => fu_advection_moment_Z_MM_ptr(PCld)
            pSpeciesIn => fu_species_transport(PCld)

          case(emission_intensity_flag, emission_flux_flag)
            pMassMapIn => fu_emisMM_ptr(PCld)
            pSpeciesIn => fu_species_emission(PCld)

          case(optical_density_flag, optical_column_depth_flag)
            !
            ! Make the optical cocktail and link it to the transport one
            !
            pOpticalRules => fu_optical_rules(chemRules)
            pChemRunSetup => fu_ChemRunSetup(chemRules)
            call set_optical_structures(pCld, pOpticalRules, &
                                        & pChemRunSetup, &
                                        & pSpeciesOut, nSpeciesOut, &
                                        & nlStdSetup, quantity)
            if(quantity == optical_density_flag)then
              pMassMapIn => fu_optical_densityMM_ptr(PCld)
            else
              pMassMapIn => fu_optical_column_depthMM_ptr(PCld)
            endif
            if(error)return
            pSpeciesIn => fu_species_optical(PCld)

          case default
            call set_error('Unknown quantity:'+fu_quantity_string(quantity),'init_mass_output')
            return
        end select  ! non-concentration quantities

        pVert  => pMassMapIn%vertTemplate
        !
        ! Having determined the quantity and species, just add the mass map link
        !
        call add_intermediate_mass_map(pMassMapIn, &
                                     & pSpeciesIn, &
                                     & quantity, &
                                     & pSpeciesOut, &
                                     & OutDef%Rules%MassMapOutLst%ptrItem(iVarIni), &
                                     & OutDef%Rules%ini_time, &
                                     & OutDef%Rules%timestep, &
                                     & model_time_step, &
                                     & ifQuantityMassMapAdded)
        if(fu_fails(ifQuantityMassMapAdded,'Failed adding mass map for MassMapIn:' + &
                                         & fu_quantity_string(quantity),'init_mass_output'))return

      endif  ! concentration quantity or all others
      !
      ! Have we added the mass map link? 
      !
      if(ifQuantityMassMapAdded)then
        call msg('Made intermediate mass map:' + fu_quantity_string(quantity))
      else
        call set_error('Failed to make intermediate mass map for:' + &
                     & fu_quantity_string(OutDef%Rules%MassMapOutLst%ptrItem(iVarIni)%quantity), 'init_mass_output')
        return
      endif

      if(nSpeciesOut > 0)then
        deallocate(pSpeciesOut)
        nSpeciesOut = 0
        nullify(pSpeciesOut)
      endif

      !Have to set if3d flag
      ! All species flags are turned to negative in MassMapOutLst, thus using abs
      if (any( OutDef%Rules%MassMapOutLst%ptrItem(iVarIni)%iVerticalTreatment &
            & == (/integrate_column_flag, lowest_level_flag/))) then
         OutDef%Rules%MassMapOutLst%ptrItem(iVarIni)%if3D = .false.
      else
         OutDef%Rules%MassMapOutLst%ptrItem(iVarIni)%if3D = (fu_NbrOfLevels(pVert) > 1)
      endif
      do iTmp = iVarIni, size(OutDef%Rules%MassMapOutLst%ptrItem)
          if(OutDef%Rules%MassMapOutLst%ptrItem(iTmp)%quantity == int_missing)exit  ! all done
          if(abs(OutDef%Rules%MassMapOutLst%ptrItem(iTmp)%quantity) == quantity)then
             OutDef%Rules%MassMapOutLst%ptrItem(iTmp)%if3D =  OutDef%Rules%MassMapOutLst%ptrItem(iVarIni)%if3D 
          endif
      enddo

    end do  ! iVarIni

    !
    ! Finally, we processed the whole list, need to clean the mess and return a clean list
    !
    iVarIni = 1
    do while (iVarIni <= size(OutDef%Rules%MassMapOutLst%ptrItem))
      if(OutDef%Rules%MassMapOutLst%ptrItem(iVarIni)%quantity == int_missing)exit ! cycle ended
      if(OutDef%Rules%MassMapOutLst%ptrItem(iVarIni)%quantity == -1)then  ! remnant of the species exploring
        call shrink_output_list(OutDef%Rules%MassMapOutLst, iVarIni)
        cycle
      endif
      if(OutDef%Rules%MassMapOutLst%ptrItem(iVarIni)%quantity < 0) &     ! species explored, flip the sign
        & OutDef%Rules%MassMapOutLst%ptrItem(iVarIni)%quantity = &
                            & - OutDef%Rules%MassMapOutLst%ptrItem(iVarIni)%quantity
      iVarIni = iVarIni + 1
    end do

    CONTAINS

    !===========================================================================
    
    subroutine add_intermediate_mass_map(pMassMapIn, pSpeciesIn, quantity, &
                                       & pSpeciesOut, pOutItem, &
                                       & ini_time, &
                                       & output_time_step, &
                                       & model_time_step, &
                                       & ifAdded)
      !
      ! Determines the link and sets the link from the pollution-cloud masMap to
      ! intermediate io_server mass map, which is created, of course.
      !
      implicit none

      ! Imported parameters
      type(Tmass_map), pointer :: pMassMapIn
      type(silam_species), dimension(:), pointer :: pSpeciesIn, pSpeciesOut
      integer, intent(in) :: quantity
      type(TOutputLstItem), intent(in) :: pOutItem   
      type(silja_time), intent(in) :: ini_time
      type(silja_interval), intent(in) :: output_time_step, model_time_step

      logical, intent(out) :: ifAdded

      ! Local variables
      integer :: iLink, nSpTmp, iSp
      logical :: ifFound
      type(silam_vertical), pointer :: pVert
      type(silam_species), dimension(:), pointer :: pSpTmp

      ifAdded = .false.

      if(.not. associated(pMassMapIn))return
      if(.not. defined(pMassMapIn))return

      !
      ! Find the free link
      !
      iLink = 1
      ifFound = .false.
      do while(iLink <= size(OutVars%MassMapLinks%iVerticalTreatment))
        if(OutVars%MassMapLinks%iVerticalTreatment(iLink) == int_missing)then
          ifFound = .true.
          exit
        endif
        iLink = iLink + 1
      end do
      if(fu_fails(ifFound,'Failed to find free space for intermediate MassMap', &
                        & 'add_intermediate_mass_map'))return
      !
      ! Questions to be asked:
      ! - inventory => list of species in the output massMap and adaptor
      ! - vertical treatment => reflect in the output massMap size and store
      ! - source split => reflect in the output massMap size and store
      !
      ! Start from making the overlapping set of species: those needed for the output
      ! with those available from the input. Note that the input can be transport, short-lived, etc
      ! In the current construction these are the different mass map links.
      !
      nullify(pSpTmp)
      nSpTmp = 0
      do iSp = 1, size(pSpeciesIn)
        if(fu_index(pSpeciesIn(iSp), pSpeciesOut, ifByName_=.true.) > 0) &
                     & call addSpecies(pSpTmp, nSpTmp, (/pSpeciesIn(iSp)/), 1)
      end do
      if(nSpTmp == 0)then
        call msg_warning('No species found for mass map link output:' + &
                       & fu_quantity_string(pMassMapIn%quantity),'add_intermediate_mass_map')
        return
      endif
      
      call create_adaptor(pSpTmp, pSpeciesIn, OutVars%MassMapLinks%adaptor(iLink), .true.)
      if(error)return
      if(OutVars%MassMapLinks%adaptor(iLink)%nSpecies < 1)return  ! nothing for this mass map
      !
      ! Time averaging: just store
      !
      OutVars%MassMapLinks%iAveragingType(iLink) = pOutItem%AvType
      OutVars%MassMapLinks%AveragingPeriod(iLink) = pOutItem%AvPeriod
      if(fu_accumulated_quantity(quantity))then
        !
        ! Cumulative input quantity.
        !
        select case(pOutItem%AvType)
          case(iAsIs, iAverage, iCumulative)
            OutVars%MassMapLinks%IntegrationStart(iLink) = ini_time
          case(iInstant)
            OutVars%MassMapLinks%IntegrationStart(iLink) = ini_time + output_time_step - model_time_step
          case(iMeanLastHrs)
            OutVars%MassMapLinks%IntegrationStart(iLink) = ini_time + output_time_step - pOutItem%AvPeriod
          case default
            call set_error('Unknown time treatment:' + fu_str(pOutItem%AvType),'add_intermediate_mass_map')
            return
        end select
      else
        !
        ! Instant input quantity.
        !
        select case(pOutItem%AvType)
          case(iAsIs, iInstant)
            OutVars%MassMapLinks%IntegrationStart(iLink) = ini_time
          case(iAverage, iCumulative)
            OutVars%MassMapLinks%IntegrationStart(iLink) = ini_time 
          case(iMeanLastHrs)
            OutVars%MassMapLinks%IntegrationStart(iLink) = ini_time + output_time_step - &
                                                         & pOutItem%AvPeriod - model_time_step
          case default
            call set_error('Unknown time treatment:' + fu_str(pOutItem%AvType),'add_intermediate_mass_map')
            return
        end select
      endif  ! type of input quantity
      !
      ! Vertical treatment options: do_nothing, integrate_column, or pick_lowest_level
      !
      OutVars%MassMapLinks%iVerticalTreatment(iLink) = pOutItem%iVerticalTreatment
      if(pOutItem%iVerticalTreatment == do_nothing_flag)then
        pVert => pMassMapIn%vertTemplate
      elseif(pOutItem%iVerticalTreatment == integrate_column_flag)then
        allocate(pVert)
        call set_vertical(entire_atmosphere_integr_level, pVert)
      elseif(pOutItem%iVerticalTreatment == lowest_level_flag)then
        allocate(pVert)
        call set_vertical(fu_level(pMassMapIn%vertTemplate,1), pVert)
      else
        call set_error('Strange vertical treatment:' + fu_str(pOutItem%iVerticalTreatment), &
                     & 'add_intermediate_mass_map')
        return
      endif
      !
      ! Finally, set the mass map for the intermadiate fields. Note that cumulative-to-whatever
      ! conversion requires intermediate mass map
      !
      OutVars%MassMapLinks%pMMIn(iLink)%ptrMassMap => pMassMapIn
      
      call  set_mass_map(OutVars%MassMapLinks%pMMOut(iLink),&
                     & quantity, &
                     & pMassMapIn%nSrc, &
                     & 0, &       ! border cells
                     & pMassMapIn%gridTemplate, &
                     & pVert, &
                     & pSpTmp, &  !eciesOut, &
                     & val=0.0)               ! val)


      ifAdded = .not. error

    end subroutine add_intermediate_mass_map


    !===========================================================================
    
    subroutine add_intermed_mass_map_4_lpSet(pLPSetIn, iLPSetType, pSpeciesIn, quantity, &
                                           & pSpeciesOut, pOutItem, gridOut, vertOut, &
                                           & ini_time, output_time_step, model_time_step, ifAdded)
      !
      ! Determines the link and sets the link from the pollution-cloud masMap to
      ! intermediate io_server mass map, which is created, of course.
      !
      implicit none

      ! Imported parameters
      type(Tlagrange_particles_set), pointer :: pLPSetIn
      integer, intent(in) :: iLPSetType
      type(silam_species), dimension(:), pointer :: pSpeciesIn, pSpeciesOut
      integer, intent(in) :: quantity
      type(TOutputLstItem), intent(in) :: pOutItem
      type(silja_grid), intent(in) :: gridOut
      type(silam_vertical), target, intent(in) :: vertOut
      type(silja_time), intent(in) :: ini_time
      type(silja_interval), intent(in) :: output_time_step, model_time_step
      logical, intent(out) :: ifAdded

      ! Local variables
      integer :: iLink, nSpTmp, iSp
      logical :: ifFound
      type(silam_vertical), pointer :: pVert
      type(silam_species), dimension(:), pointer :: pSpTmp

      ifAdded = .false.

      if(.not. associated(pLPSetIn) )return
      if(.not. associated(pSpeciesIn))return
      if(.not. (defined(pLPSetIn%gridTemplate) .and. defined(pLPSetIn%verticalTemplate)))return

      !
      ! Find the free link
      !
      iLink = 1
      ifFound = .false.
      do while(iLink <= size(OutVars%Lagr2MMLinks%iVerticalTreatment))
        if(OutVars%Lagr2MMLinks%iVerticalTreatment(iLink) == int_missing)then
          ifFound = .true.
          exit
        endif
        iLink = iLink + 1
      end do
      if(fu_fails(ifFound,'Failed to find free space for intermediate MassMap', &
                        & 'add_intermed_mass_map_4_lpSet'))return
      !
      ! Questions to be asked:
      ! - inventory => list of species in the output massMap and adaptor
      ! - vertical treatment => reflect in the output massMap size and store
      ! - source split => reflect in the output massMap size and store
      !
      ! Start from making the overlapping set of species: those needed for the output
      ! with those available from the input. Note that the input can be transport, short-lived, etc
      ! In the current construction these are the different mass map links.
      !
      nullify(pSpTmp)
      nSpTmp = 0
      do iSp = 1, size(pSpeciesIn)
        if(fu_index(pSpeciesIn(iSp), pSpeciesOut) > 0) &
                     & call addSpecies(pSpTmp, nSpTmp, (/pSpeciesIn(iSp)/), 1)
      end do

!      call create_adaptor(pSpeciesIn, pSpeciesOut, OutVars%Lagr2MMLinks%adaptor(iLink), .true.)
      call create_adaptor(pSpTmp, pSpeciesIn, OutVars%Lagr2MMLinks%adaptor(iLink), .true.)
      if(error)return
      if(OutVars%Lagr2MMLinks%adaptor(iLink)%nSpecies < 1)return  ! nothing for this lpSet
      !
      ! Time averaging: just store
      !
      OutVars%Lagr2MMLinks%iAveragingType(iLink) = pOutItem%AvType
      OutVars%Lagr2MMLinks%AveragingPeriod(iLink) = pOutItem%AvPeriod
      !
      ! Instant input quantity.
      !
      select case(pOutItem%AvType)
        case(iAsIs, iInstant)
          OutVars%Lagr2MMLinks%IntegrationStart(iLink) = ini_time
        case(iAverage, iCumulative)
          OutVars%Lagr2MMLinks%IntegrationStart(iLink) = ini_time 
        case(iMeanLastHrs)
          OutVars%Lagr2MMLinks%IntegrationStart(iLink) = ini_time + output_time_step - &
                                                       & pOutItem%AvPeriod - model_time_step
        case default
          call set_error('Unknown time treatment:' + fu_str(pOutItem%AvType),'add_intermediate_mass_map')
          return
      end select
      !
      ! Vertical treatment options: do_nothing, integrate_column, or pick_lowest_level
      !
      OutVars%Lagr2MMLinks%iVerticalTreatment(iLink) = pOutItem%iVerticalTreatment
      if(pOutItem%iVerticalTreatment == do_nothing_flag)then
        pVert => vertOut             ! lagrangian particles fly in meteo grid but projected to output
      elseif(pOutItem%iVerticalTreatment == integrate_column_flag)then
        allocate(pVert)
        call set_vertical(entire_atmosphere_integr_level, pVert)
      elseif(pOutItem%iVerticalTreatment == lowest_level_flag)then
        allocate(pVert)
        call set_vertical(fu_level(vertOut,1), pVert)
      else
        call set_error('Strange vertical treatment:' + fu_str(pOutItem%iVerticalTreatment), &
                     & 'add_intermed_mass_map_4_lpSet')
        return
      endif
      !
      ! Finally, set the mass map for the intermadiate fields. Note that cumulative-to-whatever
      ! conversion requires intermediate mass map
      !
      OutVars%Lagr2MMLinks%pLpSetIn(iLink)%ptrLpSet => pLPSetIn
      OutVars%Lagr2MMLinks%iLPSetType(iLink) = iLPSetType

      call  set_mass_map(OutVars%Lagr2MMLinks%pMMOut(iLink),&
                       & quantity, &
                       & pLPSetIn%nSrcs, &
                       & 0, &       ! border cells
                       & gridOut, & !!!pLPSetIn%gridTemplate, &
                       & pVert, &
                       & pSpTmp, & !eciesOut, &
                       & val=0.0)            ! val)
              

      ifAdded = .not. error

    end subroutine add_intermed_mass_map_4_lpSet

  end subroutine init_mass_output


  !*********************************************************************************
  
  subroutine explore_item_species(Lst, Item, pCloud, pSpecies, nSpecies, nlStdSetup)
    !
    ! Adds-up the species from the given output list item. One trouble: optics.
    ! No matter what is the species list, for optical quantity part of those may be not
    ! available: simply no optical metadata. Create a separate list of materials+modes 
    ! covering the request and force wave length into it. Then check that the obtained 
    ! species for availability from optical module.
    !
    implicit none
    
    ! Imported parameters
   type(TOutputList), intent(inout) :: Lst
    integer, intent(in) :: item
    type(silam_pollution_cloud), pointer :: pCloud
    type(silam_species), dimension(:), pointer :: pSpecies
    integer, intent (inout) :: nSpecies
    type(Tsilam_namelist), intent(in) :: nlStdSetup
    
    ! Local variables
    integer :: iSp, nSpTmp, nSpInit, iTmp
    type(silam_species), dimension(:), pointer :: pSpTmp
    real :: fWaveLen

    nSpInit = nSpecies   ! store how many were sent in

    if(fu_optics_owned_quantity(pCloud, Lst%ptrItem(item)%quantity))then
      !
      ! Optical quantity. 
      !
      nullify(pSpTmp)
      nSpTmp = 0
      !
      ! A bit ad-hoc solution, which allows giving *_INVENTORY to optical module
      ! and allowing it to through out species with unknown optical features. But
      ! if a specific species is requested with request == 2, the optical data must be available.
      !
      if(Lst%ptrItem(item)%iSpeciesListType /= iSingleSubstance) Lst%ptrItem(item)%request = 1

      if(Lst%ptrItem(item)%iSpeciesListType == iSourceInventory)then                          ! source inventory
        call addSpecies(pSpTmp, nSpTmp, &
                      & fu_species_emission(pCloud), fu_nbr_of_species_emission(pCloud))

      elseif(Lst%ptrItem(item)%iSpeciesListType == iFullInventory)then                        ! full inventory
        call addSpecies(pSpTmp, nSpTmp, &
                      & fu_species_transport(pCloud), fu_nbr_of_species_transport(pCloud))
        call addSpecies(pSpTmp, nSpTmp, &
                      & fu_species_short_lived(pCloud), fu_nbr_of_species_short_lived(pCloud))
        call addSpecies(pSpTmp, nSpTmp, &
                      & fu_species_aerosol(pCloud), fu_nbr_of_species_aerosol(pCloud))

      elseif(Lst%ptrItem(item)%iSpeciesListType == iTransportInventory)then                   ! transported inventory
        call addSpecies(pSpTmp, nSpTmp, &
                      & fu_species_transport(pCloud), fu_nbr_of_species_transport(pCloud))

      elseif(Lst%ptrItem(item)%iSpeciesListType == iAerosolInventory)then                     ! aerosol inventory
        call addSpecies(pSpTmp, nSpTmp, &
                      & fu_species_aerosol(pCloud), fu_nbr_of_species_aerosol(pCloud))

      elseif(Lst%ptrItem(item)%iSpeciesListType == iShortLivingInventory)then                 ! short-living inventory
        call addSpecies(pSpTmp, nSpTmp, &
                      & fu_species_short_lived(pCloud), fu_nbr_of_species_short_lived(pCloud))
      elseif(Lst%ptrItem(item)%iSpeciesListType == iSingleSubstance)then ! specific species
        fWaveLen = fu_optical_wave_length(Lst%ptrItem(item)%species)
        if (.not. fu_fix_matching_species(Lst%ptrItem(item)%chSpecies_string, &
                                        & Lst%ptrItem(item)%species, &
                                        & fu_species_transport(pCloud), &
                                        & fu_nbr_of_species_transport(pCloud))) then
          if (.not. fu_fix_matching_species(Lst%ptrItem(item)%chSpecies_string, &
                                          & Lst%ptrItem(item)%species, &
                                          & fu_species_short_lived(pCloud), &
                                          & fu_nbr_of_species_short_lived(pCloud))) then
            if (.not. fu_fix_matching_species(Lst%ptrItem(item)%chSpecies_string, &
                                            & Lst%ptrItem(item)%species, &
                                            & fu_species_aerosol(pCloud), &
                                            & fu_nbr_of_species_aerosol(pCloud))) then
                call set_error("Could not find  matching species in cloud","explore_item_species")
                call msg(fu_quantity_string(Lst%ptrItem(item)%quantity) +", Species:")
                call report(Lst%ptrItem(item)%species)
                call report(pcloud)
                return
            endif  ! species aerosol
          endif  ! species short lived
        endif  ! species transport
        call set_optical_wave_length(Lst%ptrItem(item)%species, fWaveLen)
        call addSpecies(pSpTmp, nSpTmp, (/Lst%ptrItem(item)%species/), 1)

      elseif(Lst%ptrItem(item)%iSpeciesListType == iNoSubstanceRelation)then
           call set_error(fu_quantity_string(Lst%ptrItem(item)%quantity) + &
                      & ": no substance relation","explore_item_species")
           return
      else
           call set_error(fu_quantity_string(Lst%ptrItem(item)%quantity) + &
                      & ": unknown SpeciesListType:"  + &
                      & fu_str(Lst%ptrItem(item)%iSpeciesListType),"explore_item_species")
           return
      endif  ! type of species list request
      !
      ! Having compiled the list of available material+mode pairs, try them one-by-one for optic output
      !
      do iSp = 1, nSpTmp
        !! For optics output only available species.. 
        call set_optical_wave_length(pSpTmp(iSp), fu_optical_wave_length(Lst%ptrItem(item)%species))
        if( fu_if_species_optics_known(pSpTmp(iSp), real_missing)) then
            call addSpecies(pSpecies, nSpecies, (/pSpTmp(iSp)/), 1)
        endif
        if (error) then 
          call set_error("error after fu_if_species_optics_known", "explore_item_species")
          return
        endif
      end do
      if(nSpTmp > 0)deallocate(pSpTmp)

    else
      !
      ! All other dispersion quantities. Just make the list.
      !
      if(Lst%ptrItem(item)%iSpeciesListType == iSourceInventory)then             ! source inventory
        call addSpecies(pSpecies, nSpecies, &
                      & fu_species_emission(pCloud), fu_nbr_of_species_emission(pCloud))
      elseif(Lst%ptrItem(item)%iSpeciesListType == iFullInventory)then                   ! full inventory
        call addSpecies(pSpecies, nSpecies, &
                      & fu_species_transport(pCloud), fu_nbr_of_species_transport(pCloud))
        call addSpecies(pSpecies, nSpecies, &
                      & fu_species_short_lived(pCloud), fu_nbr_of_species_short_lived(pCloud))
        call addSpecies(pSpecies, nSpecies, &
                      & fu_species_aerosol(pCloud), fu_nbr_of_species_aerosol(pCloud))
      elseif(Lst%ptrItem(item)%iSpeciesListType == iTransportInventory)then              ! transport inventory
        call addSpecies(pSpecies, nSpecies, &
                      & fu_species_transport(pCloud), fu_nbr_of_species_transport(pCloud))
      elseif(Lst%ptrItem(item)%iSpeciesListType == iAerosolInventory)then                ! aerosol inventory
        call addSpecies(pSpecies, nSpecies, &
                      & fu_species_aerosol(pCloud), fu_nbr_of_species_aerosol(pCloud))
      elseif(Lst%ptrItem(item)%iSpeciesListType == iShortLivingInventory)then            ! short living inventory
        call addSpecies(pSpecies, nSpecies, &
                      & fu_species_short_lived(pCloud), fu_nbr_of_species_short_lived(pCloud))
      elseif(Lst%ptrItem(item)%iSpeciesListType == iSingleSubstance)then ! specific species
!        if (.not. fu_fix_matching_species(Lst%ptrItem(item)%species, &
        if (.not. fu_fix_matching_species(Lst%ptrItem(item)%chSpecies_string, &
                                        & Lst%ptrItem(item)%species, &
                                        & fu_species_transport(pCloud), &
                                        & fu_nbr_of_species_transport(pCloud))) then
!          if (.not. fu_fix_matching_species(Lst%ptrItem(item)%species, &
          if (.not. fu_fix_matching_species(Lst%ptrItem(item)%chSpecies_string, &
                                          & Lst%ptrItem(item)%species, &
                                          & fu_species_short_lived(pCloud), &
                                          & fu_nbr_of_species_short_lived(pCloud))) then
!            if (.not. fu_fix_matching_species(Lst%ptrItem(item)%species, &
            if (.not. fu_fix_matching_species(Lst%ptrItem(item)%chSpecies_string, &
                                            & Lst%ptrItem(item)%species, &
                                            & fu_species_aerosol(pCloud), &
                                            & fu_nbr_of_species_aerosol(pCloud))) then
                call set_error("Could not find  matching species in cloud","explore_item_species")
                call msg(fu_quantity_string(Lst%ptrItem(item)%quantity) + &
                       & ", Species string:" + Lst%ptrItem(item)%chSpecies_string + ', species:')
                call report(Lst%ptrItem(item)%species)
                call report(pcloud)
                return
            endif
          endif
        endif
        call addSpecies(pSpecies, nSpecies, (/Lst%ptrItem(item)%species/), 1)      
      elseif(Lst%ptrItem(item)%iSpeciesListType == iNoSubstanceRelation)then
        call set_error(fu_quantity_string(Lst%ptrItem(item)%quantity) + &
                      & ": no substance relation","explore_item_species")
        return
      else
        call set_error(fu_quantity_string(Lst%ptrItem(item)%quantity) + &
                      & ": unknown SpeciesListType:"  + &
                      & fu_str(Lst%ptrItem(item)%iSpeciesListType),"explore_item_species")
        return
      endif  ! type of species request
    endif ! if optical quantity

    !
    ! Having the list of species created, append the corresponding items to the end of the list,
    ! just use -quantity to marke them as "extras"
    !
    if(nSpInit < nSpecies)then         ! got anything?
      !
      ! Find the last filled item
      do iTmp = 1, size(Lst%ptrItem)
        if(Lst%ptrItem(iTmp)%quantity == int_missing)exit
      end do
      if(size(lst%ptrItem) - iTmp + 1 < nSpecies-nSpInit) &
                     & call expand_output_list(Lst, nSpecies - nSpInit)  ! need extra items?
      if(error)return
      do iSp = nSpInit+1, nSpecies
        lst%ptrItem(iTmp)%quantity = - lst%ptrItem(item)%quantity   ! mark as used
        lst%ptrItem(iTmp)%request = lst%ptrItem(item)%request
        lst%ptrItem(iTmp)%AvType = lst%ptrItem(item)%AvType
        lst%ptrItem(iTmp)%iSpeciesListType = int_missing            ! individual species
        lst%ptrItem(iTmp)%iVerticalTreatment = lst%ptrItem(item)%iVerticalTreatment
        lst%ptrItem(iTmp)%if3D = lst%ptrItem(item)%if3D
        lst%ptrItem(iTmp)%AvPeriod = lst%ptrItem(item)%AvPeriod
        lst%ptrItem(iTmp)%species = pSpecies(iSp)                  ! explored species
        lst%ptrItem(iTmp)%targetId = lst%ptrItem(item)%targetId
        iTmp = iTmp + 1
      end do
      if(iTmp <= size(lst%ptrItem))lst%ptrItem(iTmp)%quantity = int_missing ! cut the rest
    endif  ! if new species have been added

  CONTAINS
  
  !!=====================================================================
  !
  !  logical function fu_fix_matching_species(species,speciesLST,nSp) 
  !    !
  !    ! Find speies with matching string and replaces "species"
  !    ! with one from the list. Needed because
  !    ! species generaed from output_config do not have full mode description
  !    !
  !    implicit none
  !    type(silam_species), intent(inout) :: species
  !    type(silam_species), dimension(:), pointer :: speciesLST
  !    integer, intent(in) :: nSp
  !    integer :: i
  !
  !    fu_fix_matching_species = .true.
  !    do i = 1, nSp
  !      if (trim(fu_str(species)) == trim(fu_str(speciesLST(i)))) then
  !        species = speciesLST(i)
  !        return
  !      endif
  !    enddo
  !    fu_fix_matching_species = .false.
  !
  !  end function fu_fix_matching_species
    
  !=====================================================================
  
    logical function fu_fix_matching_species(chSpeciesIn,species,speciesLST,nSp) 
      !
      ! Find speies with matching string and replaces "species"
      ! with one from the list. Needed because
      ! species generaed from output_config do not have full mode description
      ! A peculiarity of this function is that is gets the full output-request string
      ! and tries to find the matching one. Thus, the intermediate decoding of this string
      ! followed by its coding back to string is avoided
      !
      implicit none
      character(len=*), intent(in) :: chSpeciesIn
      type(silam_species), intent(out) :: species
      type(silam_species), dimension(:), pointer :: speciesLST
      integer, intent(in) :: nSp
      integer :: i

      fu_fix_matching_species = .true.
      do i = 1, nSp
        if (trim(chSpeciesIn) == trim(fu_str(speciesLST(i)))) then
          species = speciesLST(i)
          return
        endif
      enddo
      fu_fix_matching_species = .false.

    end function fu_fix_matching_species
    
  end subroutine explore_item_species

  !*****************************************************************************

  subroutine expand_dispersion_output_list(OutDef, PollutionCloud, em_source, pDispersionMarket, &
                                         & nlStdSetup, dynRules)
    !
    ! After the pollution cloud is created we have a complete set of
    ! variables to be used in calculations. Having it, one can create a 
    ! full list of output variables.
    ! A tiny detail, however, is that not every quantity is available for all substances
    ! and vice versa. Have to be very careful
    !
    implicit none

    ! Imported parameters with intent IN
    type(silam_output_definition), intent(inout) :: OutDef
    type(silam_source), pointer :: em_source
    type(silam_pollution_cloud), target, intent(inout) :: PollutionCloud
    type(mini_market_of_stacks), pointer :: pDispersionMarket
    type(Tsilam_namelist), pointer :: nlStdSetup
    type(Tdynamics_rules), intent(in) :: dynRules

    ! Local variables
    integer :: iVarIni, j, iVarOut, k, iTmp, nSpeciesData, nSpeciesSelect, iSpecies, iSpeciesSelect
    type(TOutputList) :: LstTmp
    logical :: ifFound
    type(silam_pollution_cloud), pointer :: ptrCld
    real, dimension(:), pointer :: fTmp
    real :: fModeVal
    type(silam_species), dimension(:), pointer :: pSpeciesData, pSpeciesSelect
    
    ptrCld => PollutionCloud
    !
    ! Expansion goes in four steps:
    ! - copy everything to temporary list
    ! - variable-by-variable copy the items that do not need expansion
    ! - expand species
    ! - expand substances and wavelengths for OD
    !
    ! Allocate temporary array and copy everyting to it
    !
    call expand_output_list(LstTmp, size(OutDef%Rules%DispOutLst%ptrItem))

    do iVarIni = 1, size(OutDef%Rules%DispOutLst%ptrItem)
      LstTmp%ptrItem(iVarIni) = OutDef%Rules%DispOutLst%ptrItem(iVarIni)
      OutDef%Rules%DispOutLst%ptrItem(iVarIni)%quantity = int_missing
    enddo 
    !
    ! Below we re-arrange the output list so that first in the list are specific 
    ! dispersion variables, which do not need expanding for example, area of risk, 
    ! some explicitly requested species, etc. And only at the end of list
    ! there will be expanded lists of source or full inventories.
    !
    ! Collect SILAM dispersion quantities, which do NOT require expanding
    ! No need for list expansion - the number of quantities is not more than the number 
    ! of existing list items
    !
    iVarOut = 1
    do iVarIni = 1, size(LstTmp%ptrItem)
      if(LstTmp%ptrItem(iVarIni)%quantity == int_missing)exit ! cycle ended
      ifFound = .false.

      !
      ! Expansion is needed if the variable for inventory
      !
      if(LstTmp%ptrItem(iVarIni)%iSpeciesListType == iSourceInventory .or. &
       & LstTmp%ptrItem(iVarIni)%iSpeciesListType == iTransportInventory .or. &
       & LstTmp%ptrItem(iVarIni)%iSpeciesListType == iAerosolInventory .or. &
       & LstTmp%ptrItem(iVarIni)%iSpeciesListType == iShortLivingInventory .or. &
       & LstTmp%ptrItem(iVarIni)%iSpeciesListType == iFullInventory) cycle  ! Species need expanding, later

      if(LstTmp%ptrItem(iVarIni)%quantity == concentration_flag .or. &
       & LstTmp%ptrItem(iVarIni)%quantity == drydep_flag .or. & 
       & LstTmp%ptrItem(iVarIni)%quantity == wetdep_flag .or. &
       & LstTmp%ptrItem(iVarIni)%quantity == concentration_2m_flag .or. &
       & LstTmp%ptrItem(iVarIni)%quantity == optical_density_flag .or. &
       & LstTmp%ptrItem(iVarIni)%quantity == optical_column_depth_flag)then 
         if(OutDef%Rules%ifSplitSizeModes) cycle                            ! Species need expanding, later
      endif

      if(LstTmp%ptrItem(iVarIni)%quantity == emission_intensity_flag) then
        LstTmp%ptrItem(iVarIni)%iSpeciesListType = iSourceInventory
        if(fu_if_eulerian_present(dynRules%simulation_type))then
          cycle                                                   ! Species to be expanded
        else
          ifFound = .false.                                       ! Lagrangian does not have it
        endif
      endif

      if(LstTmp%ptrItem(iVarIni)%quantity == advection_moment_X_flag) cycle ! to be expanded
      if(LstTmp%ptrItem(iVarIni)%quantity == advection_moment_Y_flag) cycle ! to be expanded
      if(LstTmp%ptrItem(iVarIni)%quantity == advection_moment_Z_flag) cycle ! to be expanded
      !
      ! If none of the above cases, this is a single variable, can be just copied
      ! after ensuring that it actually exists
      !
      if(fu_if_variable_available(PollutionCloud, pDispersionMarket, &
                                & LstTmp%ptrItem(iVarIni)%quantity, &
                                & LstTmp%ptrItem(iVarIni)%species, nlStdSetup, dynRules, .false.))then
        ifFound = .true.
      else
        !
        ! Problem, possibly. The variable is not available but the dispersion market is already 
        ! initialised by all chemistry- and species-related modules. There is, however, a 
        ! chance that it will be added later on from meteorology - but then no species!
        !
        if(LstTmp%ptrItem(iVarIni)%species == species_missing .and. &
         & LstTmp%ptrItem(iVarIni)%iSpeciesListType == iNoSubstanceRelation)then
!          & LstTmp%ptrItem(iVarIni)%iSpeciesListType == int_missing))then
          ifFound = .true. ! Later: fu_if_species_required(LstTmp%ptrItem(iVarIni)%quantity)
          call msg_warning('Check for species necessity is missing','expand_dispersion_output_list')
        else
          ifFound = .false.
        endif
      endif

      if(ifFound)then
        !
        ! Available: include into the final list
        !
        OutDef%Rules%DispOutLst%ptrItem(iVarOut)%quantity = LstTmp%ptrItem(iVarIni)%quantity
        OutDef%Rules%DispOutLst%ptrItem(iVarOut)%species = LstTmp%ptrItem(iVarIni)%species
        OutDef%Rules%DispOutLst%ptrItem(iVarOut)%request = LstTmp%ptrItem(iVarIni)%request
        OutDef%Rules%DispOutLst%ptrItem(iVarOut)%AvType = LstTmp%ptrItem(iVarIni)%AvType
        OutDef%Rules%DispOutLst%ptrItem(iVarOut)%AvPeriod = LstTmp%ptrItem(iVarIni)%AvPeriod
        OutDef%Rules%DispOutLst%ptrItem(iVarOut)%iVerticalTreatment = &
                                                     & LstTmp%ptrItem(iVarIni)%iVerticalTreatment
        OutDef%Rules%DispOutLst%ptrItem(iVarOut)%iSpeciesListType = &
                                                     & LstTmp%ptrItem(iVarIni)%iSpeciesListType

call msg('Specific variable added:' + fu_quantity_short_string(LstTmp%ptrItem(iVarIni)%quantity) + &
       & fu_species_output_name(LstTmp%ptrItem(iVarIni)%species))

      else
        !
        ! Variable is not available. Check the request and proceed appropriately
        !
        if(LstTmp%ptrItem(iVarIni)%request == 2)then
          call msg('Specific mandatory variable is not available:' + &
                 & fu_quantity_short_string(LstTmp%ptrItem(iVarIni)%quantity) + &
                 & ',' + fu_species_output_name(LstTmp%ptrItem(iVarIni)%species))
          call set_error('Mandatory specific variable is not available','expand_dispersion_output_list')
          return
        else
          call msg('Exclude the specific unavailable variable:' + &
                 & fu_quantity_short_string(LstTmp%ptrItem(iVarIni)%quantity) + &
                 & ',' + fu_species_output_name(LstTmp%ptrItem(iVarIni)%species))
          cycle
        endif
      endif
      iVarOut = iVarOut + 1 
    end do

    !
    ! Finally, collect quantities, which require expanding of species
    ! Note that all the wave lengths have been already handled at the reading stage
    !
    ! Careful: the list might need to be enlarged, which creates more quantities than the
    ! space available.
    !
    do iVarIni = 1, size(LstTmp%ptrItem)

      if(LstTmp%ptrItem(iVarIni)%quantity == int_missing)exit  ! cycle ended

call msg('')
call msg('Expanding quantity:' + fu_quantity_short_string(LstTmp%ptrItem(iVarIni)%quantity))

      !
      ! Set the expansion parameters. First of all, handle the SOURCE_ and FULL_ INVENTORY flags.
      ! They cover most of the problem. Since species have individual aerosol size modes already set
      ! we might need only to merge them, no splitting work is needed.
      ! Two options may show up: all quantities not related to emission, the SOURCE_INVENTORY
      ! actually means an overlap of source and transport inventories. And vice versa.
      ! So, we have to be careful and keep in mind this filtration
      !
      if(LstTmp%ptrItem(iVarIni)%iSpeciesListType == iSourceInventory)then        ! SOURCE_INVENTORY
        pSpeciesData => fu_species_emission(ptrCld)
        nSpeciesData = fu_nbr_of_species_emission(ptrCld)
        if(fu_transport_owned_quantity(ptrCld, LstTmp%ptrItem(iVarIni)%quantity))then
          pSpeciesSelect => fu_species_transport(ptrCld)
          nSpeciesSelect = fu_nbr_of_species_transport(ptrCld)
        elseif(fu_emission_owned_quantity(em_source, LstTmp%ptrItem(iVarIni)%quantity))then
          pSpeciesSelect => fu_species_emission(ptrCld)
          nSpeciesSelect = fu_nbr_of_species_emission(ptrCld)
        elseif(fu_optics_owned_quantity(ptrCld, LstTmp%ptrItem(iVarIni)%quantity))then
          !
          ! Attention. Optical species are not yet created and optical structuers are not initialised.
          ! However, it is still possible to select the output here because they will be copied from
          ! the transport ones. We should just set the correct wave lengths and choose the apprioriate
          ! subsets
          !
          nullify(pSpeciesData)
          nSpeciesData = 0
          call addSpecies(pSpeciesData, nSpeciesData, &
                        & fu_species_emission(ptrCld), fu_nbr_of_species_emission(ptrCld))
          do iSpecies = 1, nSpeciesData
            call set_optical_wave_length(pSpeciesData(iSpecies), &
                                       & fu_optical_wave_length(LstTmp%ptrItem(iVarIni)%species))
          end do
          nullify(pSpeciesSelect)
          nSpeciesSelect = 0
          call addSpecies(pSpeciesSelect, nSpeciesSelect, &
                        & fu_species_transport(ptrCld), fu_nbr_of_species_transport(ptrCld))
          do iSpecies = 1, nSpeciesSelect
            call set_optical_wave_length(pSpeciesSelect(iSpecies), &
                                       & fu_optical_wave_length(LstTmp%ptrItem(iVarIni)%species))
          end do
        else
          call set_error('Quantity:' + fu_quantity_string(LstTmp%ptrItem(iVarIni)%quantity) + &
                       & '- is neither emission nor transport owned but SOURCE_INVENTORY is requested', &
                       & 'expand_dispersion_output_list')
          return
        endif

      elseif(LstTmp%ptrItem(iVarIni)%iSpeciesListType == iTransportInventory .or. &
           & LstTmp%ptrItem(iVarIni)%iSpeciesListType == iFullInventory)then     ! FULL_INVENTORY
        !
        ! Here we should not have anything but transport or full inventory, which
        ! in this content mean the same: all quantities that could be related to aerosol and
        ! short-living mass maps have been treated in mass-map output.
        !
        pSpeciesData => fu_species_transport(ptrCld)
        nSpeciesData = fu_nbr_of_species_transport(ptrCld)
        if(fu_transport_owned_quantity(ptrCld, LstTmp%ptrItem(iVarIni)%quantity))then
          pSpeciesSelect => fu_species_transport(ptrCld)
          nSpeciesSelect = fu_nbr_of_species_transport(ptrCld)
        elseif(fu_emission_owned_quantity(em_source, LstTmp%ptrItem(iVarIni)%quantity))then
          pSpeciesSelect => fu_species_emission(ptrCld)
          nSpeciesSelect = fu_nbr_of_species_emission(ptrCld)
        elseif(fu_optics_owned_quantity(ptrCld, LstTmp%ptrItem(iVarIni)%quantity))then
          !
          ! Attention. Optical species are not yet created and optical structuers are not initialised.
          ! However, it is still possible to select the output here because they will be copied from
          ! the transport ones. We should just set the correct wave lengths.
          !
          nullify(pSpeciesSelect)
          nSpeciesSelect = 0
          call addSpecies(pSpeciesSelect, nSpeciesSelect, &
                        & fu_species_transport(ptrCld), fu_nbr_of_species_transport(ptrCld))
          do iSpecies = 1, nSpeciesSelect
            call set_optical_wave_length(pSpeciesSelect(iSpecies), &
                                       & fu_optical_wave_length(LstTmp%ptrItem(iVarIni)%species))
          end do
          pSpeciesData => pSpeciesSelect
          nSpeciesData = nSpeciesSelect 
        else
          call set_error('Quantity:' + fu_quantity_string(LstTmp%ptrItem(iVarIni)%quantity) + &
                       & '- is neither emission nor transport owned but FULL_INVENTORY is requested', &
                       & 'expand_dispersion_output_list')
          return
        endif
      elseif(LstTmp%ptrItem(iVarIni)%iSpeciesListType == iAerosolInventory)then     ! AEROSOL_INVENTORY
        call set_error('Aerosol inventory is not handled for the quantity:' + &
               & fu_quantity_string(LstTmp%ptrItem(iVarIni)%quantity),'expand_dispersion_output_list')
        return
      elseif(LstTmp%ptrItem(iVarIni)%iSpeciesListType == iShortLivingInventory)then ! SHORT_LIVING_INVENTORY
        call set_error('Short-living inventory is not handled for the quantity:' + &
               & fu_quantity_string(LstTmp%ptrItem(iVarIni)%quantity),'expand_dispersion_output_list')
        return
      else                                                                       ! special quantities
        nSpeciesData = 0
        nSpeciesSelect = 0
      endif    ! type of expansion
      !
      ! Having the list of species defined, just expand the list checking the quantity-species
      ! combinations one-by-one and adding existing ones
      !
      do iSpecies = 1, nSpeciesData
        if(.not. defined(pSpeciesData(iSpecies)))then
          call msg('Warning: strange species at:',iSpecies)
          call report(pSpeciesData(iSpecies))
          call set_error('Strange species in the list','expand_dispersion_output_list')
          return
        endif

call msg('Found species:'  + fu_species_output_name(pSpeciesData(iSpecies)))

        !
        ! Apply the filter. Note that it is void if pSpeciesData and pSpeciesSelect are the same
        !
        ifFound = .false.
        do iSpeciesSelect = 1, nSpeciesSelect
          if(pSpeciesData(iSpecies) == pSpeciesSelect(iSpeciesSelect))then
            ifFound = .true.
            exit
          endif   ! Species name found in selecton list
        enddo  ! Species selection

        if(ifFound )call check_and_add_species(pSpeciesData(iSpecies), iVarIni, iVarOut, &
                                             & OutDef%Rules%ifSplitSizeModes)  ! if use mode
        if(error)then
          call set_error('Failed to add species','expand_dispersion_output_list')
          return
        endif
      end do  ! iSpecies
      !
      ! If the pointers have been allocated for the needds of optical output, free memory here
      !
      if(fu_optics_owned_quantity(ptrCld, LstTmp%ptrItem(iVarIni)%quantity))then
        if(nSpeciesSelect > 0)deallocate(pSpeciesSelect)
        nullify(pSpeciesSelect)
        nSpeciesSelect = 0
        if(LstTmp%ptrItem(iVarIni)%iSpeciesListType == iSourceInventory)then
          if(nSpeciesData > 0)deallocate(pSpeciesData)
          nullify(pSpeciesData)
          nSpeciesData = 0
        endif
      endif
      !
      ! This is not the end of the game: some quantities require further expansion
      !
      ! Trick 1: concentrations of short-living and aerosol species are to be in the output
      !
      if(LstTmp%ptrItem(iVarIni)%quantity == concentration_flag .and. &
       & LstTmp%ptrItem(iVarIni)%iSpeciesListType == iFullInventory)then
        !
        ! Scan the short-living species
        !
        nSpeciesData = fu_nbr_of_species_short_lived(ptrCld)
        pSpeciesData => fu_species_short_lived(ptrCld)
        do iSpecies = 1, nSpeciesData
          if(.not. defined(pSpeciesData(iSpecies)))then
            call msg('Warning: strange short lived species at:',iSpecies)
            call report(pSpeciesData(iSpecies))
            call set_error('Strange species in the list','expand_dispersion_output_list')
            return
          endif

call msg('Found species:' + fu_species_output_name(pSpeciesData(iSpecies)))

          call check_and_add_species(pSpeciesData(iSpecies), iVarIni, iVarOut, &
                                   & OutDef%Rules%ifSplitSizeModes)  ! if use mode 
          if(error)then
            call set_error('Failed to add species','expand_dispersion_output_list')
            return
          endif
        end do  ! iSpecies, short lived
        !
        ! Scan the aerosol species
        !
        nSpeciesData = fu_nbr_of_species_aerosol(ptrCld)
        pSpeciesData => fu_species_aerosol(ptrCld)
        do iSpecies = 1, nSpeciesData
          if(.not. defined(pSpeciesData(iSpecies)))then
            call msg('Warning: strange aerosol species at:',iSpecies)
            call report(pSpeciesData(iSpecies))
            call set_error('Strange species in the list','expand_dispersion_output_list')
            return
          endif

call msg('Found species:' + fu_species_output_name(pSpeciesData(iSpecies)))

          call check_and_add_species(pSpeciesData(iSpecies), iVarIni, iVarOut, &
                                   & OutDef%Rules%ifSplitSizeModes)  ! if use mode
          if(error)then
            call set_error('Failed to add species','expand_dispersion_output_list')
            return
          endif
        end do   ! iSpecies aerosol
      endif    ! if concentration full inventory

    end do ! cycle over expandable variables

    CONTAINS

    !========================================================================================

    subroutine check_and_add_species(species, iVarIni, iVarOut, ifUseSizeMode, fWaveLen)
      !
      ! checks for duplicates and adds the variable to the output list
      !
      implicit none

      ! Imported parameters
      type(silam_species), intent(in) :: species
      integer, intent(in) :: iVarIni
      integer, intent(inout) :: iVarOut
      logical, intent(in) :: ifUseSizeMode
      real, intent(in), optional :: fWaveLen

      ! Local variables
      integer :: iTmp
      logical :: ifDuplicate

!      call msg('Checking for duplicates')

      do iTmp = 1, iVarOut-1
        ifDuplicate = (OutDef%Rules%DispOutLst%ptrItem(iTmp)%quantity == &
                                                         & LstTmp%ptrItem(iVarIni)%quantity) .and. &
                    & (OutDef%Rules%DispOutLst%ptrItem(iTmp)%species == species)
        if(ifDuplicate)then
          call msg_warning('Duplicated var:'+ &
                         & fu_quantity_string(LstTmp%ptrItem(iVarIni)%quantity) + ',' + &
                         & fu_species_output_name(species), &
                         & 'check_and_add_species')
          return
        endif   ! mode is used and also coinsides
      end do  ! check for duplication
      !
      ! OK, complete this variable, ensuring that the place is available
      ! Since iMode is to be taken into account, have to add it to the variable name.
      !

      call msg('Final availability check:' + fu_quantity_string(LstTmp%ptrItem(iVarIni)%quantity) + &
             & ',' + fu_species_output_name(species))

      if(fu_if_variable_available(PollutionCloud, pDispersionMarket, &
                      & LstTmp%ptrItem(iVarIni)%quantity, species, nlStdSetup, dynRules, .true.))then
        !
        ! Finally, all checks passed, the variable can be added
        !
        if(iVarOut > size(OutDef%Rules%DispOutLst%ptrItem)) &
                                           & call expand_output_list(OutDef%Rules%DispOutLst, 2)
        if(error)return
        OutDef%Rules%DispOutLst%ptrItem(iVarOut)%quantity = LstTmp%ptrItem(iVarIni)%quantity
        OutDef%Rules%DispOutLst%ptrItem(iVarOut)%species = species
        OutDef%Rules%DispOutLst%ptrItem(iVarOut)%request = LstTmp%ptrItem(iVarIni)%request
        OutDef%Rules%DispOutLst%ptrItem(iVarOut)%AvType = LstTmp%ptrItem(iVarIni)%AvType
        OutDef%Rules%DispOutLst%ptrItem(iVarOut)%AvPeriod = LstTmp%ptrItem(iVarIni)%AvPeriod
        OutDef%Rules%DispOutLst%ptrItem(iVarOut)%iVerticalTreatment = &
                                                        & LstTmp%ptrItem(iVarIni)%iVerticalTreatment
        if(defined(species))then
          OutDef%Rules%DispOutLst%ptrItem(iVarOut)%iSpeciesListType = iSingleSubstance
        else
          OutDef%Rules%DispOutLst%ptrItem(iVarOut)%iSpeciesListType = iNoSubstanceRelation
        endif

        call msg('The above species is added')

        iVarOut = iVarOut + 1
      else
        !
        ! The quantity-species combination is not available. However, this is possible and error
        ! should not be set even if request is 2. 
        !
        call msg('Exclude the unavailable variable:' + &
               & fu_quantity_short_string(LstTmp%ptrItem(iVarIni)%quantity) + &
               & ',' + fu_species_output_name(species))
      endif

    end subroutine check_and_add_species

  end subroutine expand_dispersion_output_list


  !***********************************************************************

  subroutine read_output_configuration(chFNm, OutDef, model_time_step)
    !
    ! Reads the output configuration file and fills-in the
    ! output lists. It also distributes the requested quantities to meteo list
    ! and dispersion list. 
    ! Do not mix-up the output lists and the shopping list !!! 
    !
    implicit none

    ! Imported variables 
    character(len=*), intent(in) :: chFNm
    type(silam_output_definition), intent(inout) :: OutDef
    type(silja_interval), intent(in) :: model_time_step

    ! Local variables
    !
    integer :: iVar, io_status, iMeteoVar, iDispVar, iAvTmp, iRTmp, nVars, iTmp, k, iWave, nWaves, &
             & iInventory, iInventoryTmp, indSupplem, iVerticalTreatment, iSp
    type(silja_interval) :: intAvTmp
    character(len=clen) :: chVar, chSubstNm, chAveraging ! Not too long
    character(len=fnlen) :: chContent
    real, dimension(:), pointer :: arTmp
    type(Tsilam_namelist), pointer :: nlOutCfg
    type(Tsilam_nl_item_ptr), dimension(:), pointer :: ptrVars
    type(silam_sp) :: spContent, spSupplem, sp
    type(silam_species) :: speciesTmp
    type(silam_species), dimension (max_species) :: CustomList
    integer :: nCustomSpecies
    real :: avUnit

    nullify(ptrVars)
    iVar = fu_next_free_unit()
    iMeteoVar = 1
    iDispVar = 1

    ! Get the namelist
    !
    OPEN(file=chFNm, unit=iVar, action='read', status='old', iostat=io_status)
    IF(fu_fails(io_status == 0,'Cannot open input file:' + chFNm, 'read_output_configuration'))RETURN

    nlOutCfg => fu_read_namelist(iVar, .false., 'END_OUTPUT_CONFIG_3_7')
    if(error)return
    close(iVar)

    !
    ! Check if we have custom inventory specified
    !
    nCustomSpecies = 0
    chContent = fu_content(nlOutCfg,'custom_inventory')
    if (len(trim(chContent)) > 0) then
       k=0
       do while (k < len(chContent)) 
         k = k + 1
         if (chContent(k:k) == ' ') cycle
         iTmp = index(chContent(k:),' ')
         nCustomSpecies = nCustomSpecies + 1
         CustomList(nCustomSpecies) = fu_species_from_short_string(chContent(k:k+iTmp-1))
         k = k + iTmp -1
       enddo
    endif


    !
    ! Set the rule for aerosol size mode reporting (sum-up or separately)
    !
    chVar = fu_content(nlOutCfg,'aerosol_size_modes')
    if(chVar == 'SUM')then
      OutDef%Rules%ifSplitSizeModes = .false.
    elseif(chVar == 'SEPARATE')then
      OutDef%Rules%ifSplitSizeModes = .true.
    else
      call set_error('aerosol_size_modes =' + chVar + ', must be SUM or SEPARATE', &
                   & 'read_output_configuration')
      return
    endif

    !
    ! Get all variables and analyse them
    !
    call get_items(nlOutCfg, 'out_var', ptrVars, nVars)
    if(nVars < 1)then
      call report(nlOutCfg)
      call set_error('No variables in namelist','read_output_configuration')
      return
    endif

    !
    ! The output variable has the following attribuites:
    ! - quantity name, which is to be found now
    ! - substance names or list of substance names: substance will be used to set species, 
    !                                               list saved in iSpeciesListType and expanded later
    ! - averaging type, which is to be saved into avType
    ! - list of wave lengths, which is to be handled now
    !
    ! The only reason why we cannot expand the list of substances is because now we do not know yet
    ! what will be in the run. However, chemicals have been initialised, so species can be set right away.
    !
    OutDef%Rules%ifRunDispersion = .false.
    spContent%sp => fu_work_string()
    spSupplem%sp => fu_work_string()
    arTmp => fu_work_array()
    sp%sp => fu_work_string()
    
    call msg('')
    call msg('The list of requested variables in the output config file')

    do iVar = 1,nVars
      spContent%sp = fu_content(ptrVars(iVar))

      if(index(spContent%sp,'0') == 1) cycle  ! No need in this variable
      call msg(spContent%sp)
      !
      ! Each line is either "Request, quantity, chSubstNm, averaging, supplementary"
      ! or                  "Request, quantity, averaging, supplementary"
      ! chSubstNm is taken into [], so can be easily traced
      ! supplementary information/request starts from %<info_type>, so can be traced too
      !
      ! Let's start from the supplementary info
      ! At present, we understand two: list of wavelengths and column averaging request.
      ! Have to check them one by one
      !
      iVerticalTreatment = do_nothing_flag
      nWaves = 0
      if(index(spContent%sp,'%')> 0)then
        indSupplem = index(fu_str_u_case(spContent%sp),'%INTEGRATE_COLUMN')
        if(indSupplem > 0)then
          iVerticalTreatment = integrate_column_flag
          spContent%sp(indSupplem:indSupplem+17) = ' '
        else
          indSupplem = index(fu_str_u_case(spContent%sp),'%LOWEST_LEVEL')
          if(indSupplem > 0)then
            iVerticalTreatment = lowest_level_flag
            spContent%sp(indSupplem:indSupplem+14) = ' '
          endif
        endif

        indSupplem = index(fu_str_u_case(spContent%sp),'%WAVE_LENGTH')
        if(indSupplem > 0)then
          spSupplem%sp = spContent%sp(indSupplem+1:len_trim(spContent%sp))
          spContent%sp(indSupplem+1:len_trim(spContent%sp)) = ' '
          nWaves = 1  ! at least one must exist
        else
          spSupplem%sp = ''
        endif
      endif  ! supplem info exists
      !
      ! Now let's handle the name fo the substance - or a list of them
      !
      if(index(spContent%sp,'[') > 0)then
        read(unit=spContent%sp, fmt = *) iRTmp, chVar, chSubstNm, chAveraging  
        if(index(chSubstNm,']') /= len_trim(chSubstNm))then
          call set_error('Missing end bracket in the substance (or list of) name:' + chSubstNm, &
                       & 'read_output_configuration')
          return
        else
          !
          ! Remove brackets
          !
          iTmp=len_trim(chSubstNm)
          chSubstNm(1:iTmp-2) = chSubstNm(2:iTmp-1)
          chSubstNm(iTmp-1:iTmp) = ' '
          !
          ! May be a true substance name or one of INVENTORIES, which have to be expanded later
          !
          if(index(chSubstNm,'SOURCE_INVENTORY') == 1)then
            iInventory = iSourceInventory
            speciesTmp = species_missing
          elseif(index(chSubstNm,'FULL_INVENTORY') == 1)then
            iInventory = iFullInventory
            speciesTmp = species_missing
          elseif(index(chSubstNm,'TRANSPORT_INVENTORY') == 1)then
            iInventory = iTransportInventory
            speciesTmp = species_missing
          elseif(index(chSubstNm,'AEROSOL_INVENTORY') == 1)then
            iInventory = iAerosolInventory
            speciesTmp = species_missing
          elseif(index(chSubstNm,'SHORT_LIVING_INVENTORY') == 1)then
            iInventory = iShortLivingInventory
            speciesTmp = species_missing
          elseif(index(chSubstNm,'CUSTOM_INVENTORY') == 1)then
            iInventory = iCustomInventory
            speciesTmp = species_missing
          else
            iInventory = iSingleSubstance
            speciesTmp = fu_species_from_short_string(chSubstNm)
            if(error)then
              call msg_warning('Possible problems with the species string:'+chSubstNm,'read_output_configuration')
              call unset_error('read_output_configuration')
            endif
          endif  ! whether inventory or the substance name
          if(error)return
          if (any(fu_get_SILAM_quantity(chVar) == (/drydep_flag,wetdep_flag/))) then
            if (any(iInventory == (/iFullInventory, iShortLivingInventory/))) then
              call msg_warning("No depositions for short-lived species",  'read_output_configuration')
              call msg("Please correct your output config.")
              call msg("SOURCE_INVENTORY or TRANSPORT_INVENTORY might be good options")
              call set_error("Inventory "//trim(chSubstNm)//" is not allowed for "//trim(chVar), 'read_output_configuration')
              return
            endif
          endif 

        endif
      else
        read(unit=spContent%sp, fmt = *) iRTmp, chVar, chAveraging
        chSubstNm = ''
        iInventory = iNoSubstanceRelation
        speciesTmp = species_missing
      endif    ! substance name or inventory label

      !
      !  Decode the type of averaging
      !
      if(chAveraging == 'AS_IS')then  ! No modifications - just copy of the field
        iAvTmp = iAsIs
      elseif(chAveraging == 'INSTANT')then  ! No averaging - instant fields
        iAvTmp = iInstant
      elseif(chAveraging == 'AVERAGE')then  ! Mean between the output periods
        iAvTmp = iAverage
      elseif(chAveraging == 'CUMULATIVE')then ! Cumulative from the beginning
        iAvTmp = iCumulative
      elseif(chAveraging == 'TOTAL_WHOLE_PERIOD')then ! One value for the whole run period
        iAvTmp = iTotalWholePeriod
      elseif(chAveraging == 'TECHNICAL_ORIGINAL_GRID')then ! One value for the whole run period
        iAvTmp = iTechnicalOriginalGrid
      elseif(chAveraging == 'TECHNICAL_DISPERSION_GRID')then ! One value for the whole run period
        iAvTmp = iTechnicalDispersionGrid
      elseif(index(chAveraging,'MEAN_LAST_') == 1)then ! Mean last X hrs up to output time
        iAvTmp = iMeanLastHrs
        if(index(chAveraging,'_HR') /= 0)then
                avunit=3600.
        elseif (index(chAveraging,'_MN') /= 0) then
                avunit=60.
        else
          call set_error('Must be MEAN_LAST_?_HR or MEAN_LAST_?_MN','read_output_configuration')
          call msg(chAveraging)
          return
        endif
        iTmp=len_trim(chAveraging)
        chAveraging(1:iTmp-10) = chAveraging(11:iTmp)
        chAveraging = chAveraging(1:index(chAveraging,'_')-1)
        if(len_trim(chAveraging) == 0)then
          call set_error('Failed to understand averging','read_output_configuration')
          call msg(fu_connect_strings('Variable:',chVar))
          return
        endif
        if(index(chAveraging, '.') == 0) chAveraging = fu_connect_strings(chAveraging,'.')
        read(unit = chAveraging, fmt=*, iostat=io_status) arTmp(1)
        if(io_status /= 0)then
          call set_error(fu_connect_strings('Can not get averaging period:', &
                                          & chAveraging), &
                       & 'read_output_configuration')
          return
        endif
        intAvTmp = fu_set_interval_sec(arTmp(1) * avunit)
        !
        ! We have to check that this averaging period is an even number of the 
        ! model time steps and adjust it if it is not
        !
        if(.not.(model_time_step * real(abs(nint(intAvTmp / model_time_step))) == &
               & intAvTmp))then
          call report(OutDef%rules%timestep)
          call report(intAvTmp)
          call report(model_time_step)
          call msg("intAvTmp / model_time_step", intAvTmp / model_time_step)
          call report(model_time_step * real(abs(nint(intAvTmp / model_time_step))))
          call msg_warning('Output averaging period is adjusted to model timestep', &
                         & 'read_output_configuration')
          intAvTmp = model_time_step * real(abs(nint(intAvTmp / model_time_step)))
        endif
      else
        call set_error(fu_connect_strings('Unknown averaging type:',chAveraging), &
                     & 'read_output_configuration')
        return
      endif     ! Types of averaging

      !
      ! Put the variable names and request to the appropriate output list. 
      ! If it is SILAM dispersion quantity - set the ifRunDispersion switch
      ! Note that if there are more than one wave length, it can be decoded right here
      !
      ! Procedure: get the list of wavelengths to be considered; then, for each wave length 
      ! expand, if needed, the substance names and modes
      !
      ! Get the supplementary info: wave lengths, their unit and values
      !
      if(nWaves > 0)then
        read(unit=spSupplem%sp, fmt=*, iostat=iTmp) sp%sp
        if(iTmp /= 0)then
          call set_error('Failed to read supplementary string for optical density:' + spSupplem%sp, &
                       & 'read_output_configuration')
          return
        endif
        if(.not. trim(fu_str_u_case(sp%sp)) == 'WAVE_LENGTH')then
          call set_error('Strange supplementary string for optical density:' + spSupplem%sp, &
                       & 'read_output_configuration')
          return
        endif
        !
        ! .. and unit
        !
        read(unit=spSupplem%sp, fmt=*, iostat=io_status) sp%sp, sp%sp
        if(io_status /= 0)then
          call set_error('Failed to read unit for wavelength of optical density:' + spSupplem%sp, &
                       & 'read_output_configuration')
          return
        endif
        !
        ! Scailng factor for wavelength so SI unit [m]
        !
        arTmp(1001) = fu_conversion_factor(sp%sp,'m')
        if(error)return
        !
        ! .. and the values. Re-read twice to get the number of values
        !
        do iTmp = 1, 1000
          read(unit=spSupplem%sp, fmt=*, iostat=io_status) sp%sp, sp%sp, (arTmp(k),k=1,iTmp)
          if(io_status /= 0)then
            nWaves = iTmp-1
            exit
          endif
        end do
        read(unit=spSupplem%sp, fmt=*, iostat=io_status) sp%sp, sp%sp, (arTmp(k),k=1,nWaves)
        !
        ! Scale to SI
        !
        arTmp(1:nWaves) = arTmp(1:nWaves) * arTmp(1001)
      else
        nWaves = 1  ! overwrite the zero number
        arTmp(1) = real_missing
        
      endif  ! nWaves > 0
      !
      ! Having the number of wave lengths, let's create the corresponding species and add 
      ! the resulting output variables to the output list.
      ! Note that the consideration covers also the variables not attributed to any inventory or species
      !

      iInventoryTmp = iInventory
      do iSp = 1, max(1, nCustomSpecies)

         if (iInventory == iCustomInventory) then !Substitute species
           iInventoryTmp = iSingleSubstance
           speciesTmp = CustomList(iSp)
         endif

         do iWave = 1, nWaves

           call set_optical_wave_length(speciesTmp, arTmp(iWave))
           if(error)return

           iTmp =  fu_get_SILAM_quantity(chVar)
           if(fu_SILAM_dispersion_quantity(iTmp))then
             !
             ! Dispersion quantity is to be sent to dispersion list
             !
             if(iDispVar > size(OutDef%Rules%DispOutLst%ptrItem)) &
                         & call expand_output_list(OutDef%Rules%DispOutLst, 2)
             OutDef%Rules%ifRunDispersion = .true.
             OutDef%Rules%DispOutLst%ptrItem(iDispVar)%quantity = iTmp
             OutDef%Rules%DispOutLst%ptrItem(iDispVar)%request = iRTmp
             OutDef%Rules%DispOutLst%ptrItem(iDispVar)%AvType = iAvTmp
             OutDef%Rules%DispOutLst%ptrItem(iDispVar)%AvPeriod = intAvTmp
             OutDef%Rules%DispOutLst%ptrItem(iDispVar)%species = speciesTmp  ! missing if no substance defined
             OutDef%Rules%DispOutLst%ptrItem(iDispVar)%chSpecies_string = chSubstNm  ! just input string 
             OutDef%Rules%DispOutLst%ptrItem(iDispVar)%iSpeciesListType = iInventoryTmp  ! can be iNoSubstanceRelation
             OutDef%Rules%DispOutLst%ptrItem(iDispVar)%iVerticalTreatment = iVerticalTreatment
             OutDef%Rules%DispOutLst%ptrItem(iDispVar)%if3D = &
                & (fu_multi_level_quantity(iTmp) .and. (iVerticalTreatment == do_nothing_flag) )
             iDispVar = iDispVar + 1
             
           else   
             !
             ! non-SILAM dispersion quantity
             !
             if(iMeteoVar > size(OutDef%Rules%MeteoOutLst%ptrItem)) &
                               & call expand_output_list(OutDef%Rules%MeteoOutLst, 2)
             OutDef%Rules%MeteoOutLst%ptrItem(iMeteoVar)%quantity = iTmp
             OutDef%Rules%MeteoOutLst%ptrItem(iMeteoVar)%species = speciesTmp
             OutDef%Rules%MeteoOutLst%ptrItem(iMeteoVar)%chSpecies_string = chSubstNm  ! just input string 
             OutDef%Rules%MeteoOutLst%ptrItem(iMeteoVar)%iSpeciesListType = iInventoryTmp
             OutDef%Rules%MeteoOutLst%ptrItem(iMeteoVar)%request = iRTmp
             OutDef%Rules%MeteoOutLst%ptrItem(iMeteoVar)%AvType = iAvTmp
             OutDef%Rules%MeteoOutLst%ptrItem(iMeteoVar)%AvPeriod = intAvTmp
             OutDef%Rules%MeteoOutLst%ptrItem(iMeteoVar)%iVerticalTreatment = iVerticalTreatment
             OutDef%Rules%MeteoOutLst%ptrItem(iMeteoVar)%if3D = &
                & (fu_multi_level_quantity(iTmp) .and. (iVerticalTreatment == do_nothing_flag) )
             iMeteoVar = iMeteoVar + 1
             
           endif  ! SILAM dispersion quantity

           if(error)return

         end do  ! nWaves

         if (iInventory /= iCustomInventory) exit !No need to loop over custom inventary
       enddo


    end do    ! Cycle through the file

    call msg('done read_output_configuration')
    call msg('')

    call free_work_array(spContent%sp)
    call free_work_array(spSupplem%sp)
    call free_work_array(sp%sp)
    call free_work_array(arTmp)
    call destroy_namelist(nlOutCfg)
    deallocate(ptrVars)

  end subroutine read_output_configuration


  !***************************************************************************

  subroutine shrink_output_list(OutLst, ind)
    !
    ! Deletes one quantity from the output list
    !
    implicit none

    ! Imported parameters
    type(TOutputList), intent(inout) :: OutLst
    integer, intent(in) :: ind

    ! Local variables
    integer :: i

    if(ind < 1)then
      call set_error('Too small index','shrink_output_list')
      return
    endif
    !
    ! Just scan the list from the i-th position shifting all variabls one index up
    !
    do i = ind, size(OutLst%ptrItem)-1
      OutLst%ptrItem(i) = OutLst%ptrItem(i+1)
      if(OutLst%ptrItem(i)%quantity == int_missing)return
    end do
    OutLst%ptrItem(size(OutLst%ptrItem))%quantity = int_missing

  end subroutine shrink_output_list


  !*******************************************************************************

  subroutine report_output_list(list)
    !
    ! Just prints the complete list of the output variables
    !
    implicit none

    type(TOutputList), intent(in) :: list

    ! local variables
    integer :: i, j
    character(len=10),dimension(3), parameter :: chTypes = (/'Not needed','Desirable ','Mandatory '/)
    character(len=2),dimension(3), parameter :: XD = (/'XD','2D','3D'/)


    
    if (.not. allocated(list%ptrItem)) then
       call msg("----Unallocated list-----")
       return
    endif

    call msg('')
    call msg('The list of output variables is:')
    do i=1,size(list%ptrItem)
      if(list%ptrItem(i)%quantity == int_missing .or. list%ptrItem(i)%quantity < 1)cycle
         j=2
         if (list%ptrItem(i)%if3D) j=3
      call msg(chTypes(list%ptrItem(i)%request+1) + ':' + &
             & fu_quantity_string(list%ptrItem(i)%quantity)+"," + XD(j) + ',' + &
             & fu_str(list%ptrItem(i)%species), list%ptrItem(i)%AvType)
    end do
    call msg('End of output list')
    call msg('---------------------------------------------------------------')
    call msg('')
  end subroutine report_output_list


  !**************************************************************************

  function fu_species_output_name(species) result(chNm)
    !
    ! Generates the output name from the substance name, aerosol mode and 
    ! cocktail template
    !
    ! 9/09: The modes now always range from 1 to nModes(iSubst). The
    ! first can be gas, but has to be checked from the species.
    ! 
    implicit none

    ! Output
    character(len=clen) :: chNm

    ! Imported parameters
    type(silam_species), intent(in) :: species
    
    if(.not. defined(species))then
      chNm = 'undefined'
      return
    endif
    !
    ! Component one: substance name
    !
    chNm = fu_substance_name(species)
    !
    ! Component two: aerosol mode
    !
    if(defined(species%mode))then

      if (species%mode == in_gas_phase) then
        chNm = chNm + '_gas'
      elseif (.not. species%mode == no_mode) then
        !
        ! Within the aerosol mode index range, this is an aerosol mode
        !
        chNm = chNm + '_m' + fu_aerosol_mode_size_to_str(fu_nominal_d(species))
      end if

    else
      ! Would be preferable to use an actual mode number. 
      call msg_warning('Why undefined Mode?', 'fu_species_output_name')
    endif

    if(.not. (species%waveLength .eps. real_missing))then
      chNm = chNm + '_w' + fu_optical_wave_length_to_str(species%wavelength)
    endif

  end function fu_species_output_name


  !*****************************************************************************
  
  function fu_output_template(OutDef)
    !
    ! Returns output template
    !
    implicit none
    
    ! Imported parameters
    type(silam_output_definition), pointer :: OutDef

    ! Return value
    type(grads_template), pointer :: fu_output_template
    
    if(associated(OutDef))then
      fu_output_template => OutDef%Rules%outTemplate
    else
      nullify(fu_output_template)
    endif
  
  end function fu_output_template


  !*****************************************************************************
  
  function fu_output_dyn_shopping_list(OutDef, now) result(shList)
    !
    ! Makes-up the shopping list for the output miniMarket. DO NOT MIX IT with
    ! intermediate structures for the various output lines. Sofar, this creature
    ! is just for the dx, dy, dz fields. Later, meteo output can be dumped there.
    !
    implicit none
    
    ! Imported parameter
    type(silam_output_definition), pointer :: OutDef
    type(silja_time), intent(in) :: now
    ! return value
    type(silja_shopping_list) :: shList

    ! So far, only dz, may be.
    !
    if(fu_if_level_meteo_dependent(fu_leveltype(output_vertical)))then
      shList = fu_set_shopping_list(met_src_missing, &
                                  & (/cell_size_z_flag/), &
                                  & now, &
                                  & now, &
                                  & level_missing,&
                                  & level_missing, &
                                  & output_grid, &
                                  & (/2/), &
                                  & (/silja_false/))
    else
      call set_missing(shList)
    endif

  end function fu_output_dyn_shopping_list


  !*****************************************************************************
  
  function fu_output_st_shopping_list(OutDef, now) result(shList)
    !
    ! Makes-up the shopping list for the output miniMarket. DO NOT MIX IT with
    ! intermediate structures for the various output lines. Sofar, this creature
    ! is just for the dx, dy, dz fields. Later, meteo output can be dumped there.
    !
    implicit none
    
    ! Imported parameter
    type(silam_output_definition), pointer :: OutDef
    type(silja_time), intent(in) :: now
    ! return value
    type(silja_shopping_list) :: shList
    integer, parameter, dimension(:) :: &
            rainlist(1:3)=(/large_scale_rain_int_flag, convective_rain_int_flag, & 
                      & total_precipitation_rate_flag/)

    ! So far, only dx, dy
    !
    if(fu_if_level_meteo_dependent(fu_leveltype(output_vertical)))then
      shList = fu_set_shopping_list(met_src_missing, &
                                    & (/cell_size_x_flag, cell_size_y_flag/), &
                                    & now, &
                                    & now, &
                                    & level_missing,&
                                    & level_missing, &
                                    & output_grid, &
                                    & (/2,2/), &
                                    & (/silja_true, silja_true/))
    else
      shList = fu_set_shopping_list(met_src_missing, &
                                    & (/cell_size_x_flag, cell_size_y_flag, cell_size_z_flag/), &
                                    & now, &
                                    & now, &
                                    & level_missing,&
                                    & level_missing, &
                                    & output_grid, &
                                    & (/2,2,2/), &
                                    & (/silja_true, silja_true, silja_false/))
    endif
!    call add_shopping_quantities(shList, rainlist, (/2,2,2/), &
!                        &(/silja_true, silja_true, silja_true/))

  end function fu_output_st_shopping_list


  !*****************************************************************************
  
  subroutine align_OutDef_initial_time_with_shift(OutDef, wdr)
    !
    ! There can be time shift required from the meteo data interface. Output must
    ! comply to it. Here we align these definitions.
    !
    implicit none
    
    ! Imported parameters
    type(silam_output_definition), target, intent(inout) :: OutDef
    type(silja_wdr), intent(in) :: wdr
    
    ! Local variables
    integer :: iStack
    type(silja_time) :: new_initial_time
    type(silja_stack), pointer :: stackPtr
    
    !
    ! First, decide if the action is possible and needed
    !
    if(fu_fails(defined(wdr), 'undefined wdr or OutDef given', 'align_OutDef_initial_time_with_shift'))return
    if(.not. defined(fu_meteo_time_shift(wdr))) return         ! nothing to do
    if(.not. fu_meteo_time_shift(wdr) > zero_interval) return  ! nothing to do
    
    new_initial_time = OutDef%Rules%ini_time + fu_meteo_time_shift(wdr)
    !
    ! Now, do the updates, both dispersion and meteorological stacks
    !    
    do iStack = 1, size(OutVars%DispTmpStack)
      stackPtr => OutVars%DispTmpStack(iStack)
      call update_stack_fields(stackPtr, new_initial_time)
      if(error) return
    enddo
    
    call update_stack_fields(OutVars%MeteoTmpStack, new_initial_time)

  end subroutine align_OutDef_initial_time_with_shift



  !*********************************************************
  !
  ! ENCAPSULATION
  !
  !*********************************************************

  function fu_caseNM(def)
    implicit none
    character(len=clen) :: fu_caseNM
    type(silam_output_definition), intent(in) :: def
    fu_caseNm = def%params%chCaseNm
  end function fu_caseNm

  !-----------------------------------------------------------------------------------------
  logical function fu_ifTrajectory_in_output(def)
    implicit none
    type(silam_output_definition), intent(in) :: def
    fu_ifTrajectory_in_output = def%Rules%ifTrajectory
  end function fu_ifTrajectory_in_output

  !-----------------------------------------------------------------------------------------
  function fu_timestep_of_output(def)result(step)
    implicit none
    type(silja_interval) :: step
    type(silam_output_definition), intent(in) :: def
    step = def%Rules%timestep
  end function fu_timestep_of_output

  !-----------------------------------------------------------------------------------------
  logical function fu_ifRunDispersion(OutDef)
    implicit none
    type(silam_output_definition), intent(in) :: OutDef
    fu_ifRunDispersion= OutDef%RUles%ifRunDispersion
  end function fu_ifRunDispersion

  !-----------------------------------------------------------------------------------------
  function fu_traj_set(OutDef) result(traj_set)
    implicit none
    type (silam_trajectory_set), pointer :: traj_set
    type(silam_output_definition), intent(in), target :: OutDef
    traj_set => OutDef%Params%tr_set
  end function fu_traj_set

  !-----------------------------------------------------------------------------------------
  integer function fu_model_output_request(OutDef)
    implicit none
    type(silam_output_definition), intent(in) :: OutDef
    integer :: i
    fu_model_output_request = 0
    if (allocated(outDef%rules%dispOutLst%ptrItem)) then
      do i=1,size(OutDef%rules%DispOutLst%ptrItem)
        if(OutDef%rules%DispOutLst%ptrItem(i)%quantity == int_missing)exit
        if(fu_silam_dispersion_quantity(OutDef%rules%DispOutLst%ptrItem(i)%quantity)) &
             & fu_model_output_request= max(fu_model_output_request, &
                                          & OutDef%Rules%DispOutLst%ptrItem(i)%request)
      end do
    end if
    if (allocated(outDef%rules%meteoOutLst%ptrItem)) then
      do i=1,size(OutDef%rules%MeteoOutLst%ptrItem)
        if(OutDef%rules%MeteoOutLst%ptrItem(i)%quantity == int_missing)exit
        if(fu_silam_dispersion_quantity(OutDef%rules%MeteoOutLst%ptrItem(i)%quantity)) &
             & fu_model_output_request= max(fu_model_output_request, &
                                          & OutDef%Rules%MeteoOutLst%ptrItem(i)%request)
      end do
    end if

  end function fu_model_output_request

  !-----------------------------------------------------------------------------------------
  subroutine get_output_quantities(OutDef, qArr)
    implicit none
    integer, dimension(:), intent(out) :: qArr  ! return array of quantities
    type(silam_output_definition), target, intent(in) :: OutDef
    integer :: iTmp, jTmp
     jTmp = size(OutDef%rules%DispOutLst%ptrItem)
    do iTmp = 1, jTmp
      qArr(iTmp) = OutDef%rules%DispOutLst%ptrItem(iTmp)%quantity
    end do
    qArr(jTmp+1) = int_missing !Put end marker
  end subroutine get_output_quantities


  !*********************************************************************************

  function fu_ifDryDep_cumulative(OutDef)
    type(silja_logical) :: fu_ifDryDep_cumulative
    type(silam_output_definition), intent(in) :: OutDef
    fu_ifDryDep_cumulative = OutDef%rules%ifDryDepCumulative
  end function fu_ifDryDep_cumulative

  !*********************************************************************************
  
  function fu_ifWetDep_cumulative(OutDef)
    type(silja_logical) :: fu_ifWetDep_cumulative
    type(silam_output_definition), intent(in) :: OutDef
    fu_ifWetDep_cumulative = OutDef%rules%ifWetDepCumulative
  end function fu_ifWetDep_cumulative

  !*********************************************************************************

  subroutine report_output_definition(def)
    !
    ! Prints all output deifnitions
    !
    implicit none

    ! Imported parameters witrh intent IN
    type(silam_output_definition), intent(in) :: def

    ! Local declarations
    integer :: i

    call msg('=============================================================')
    call msg('')
    call msg('Output parameters: ')
    call msg('')

    call report(output_grid)

    call msg('')
    call msg('Output levels')
    call report(output_vertical)

    call msg('Output interval')

    call report(def%Rules%timestep)

    call msg('Output file formats:')
    if(def%Rules%iGrib /= int_missing) call msg('GRIB')
    if(def%Rules%ifGrads) call msg('GRADS')
    if(def%Rules%ifTrajectory ) call msg('TRAJECTORY')
    if(def%Rules%iNetcdf /= int_missing) call msg('NETCDF')

    call msg('')
    call msg('Time series arrangement: ')
    select case(def%Rules%OutFilesArrangement)
      case(all_in_one)
        call msg('ALL_IN_ONE')
      case(hourly_new_file)
        call msg('HOURLY_NEW_FILE')
      case(daily_new_file)
        call msg('DAILY_NEW_FILE')
      case(monthly_new_file)
        call msg('MONTHLY_NEW_FILE')
      case(yearly_new_file)
        call msg('YEARLY_NEW_FILE')
      case default
        call msg('Unknown')
    end select

    call msg('')
    call msg(fu_connect_strings('Output template:',fu_template_string(def%Rules%OutTemplate)))
    call msg("Rules%MeteoOutLst")
     call report_output_list(def%Rules%MeteoOutLst)
    call msg("Rules%DispOutLst")
     call report_output_list(def%Rules%DispOutLst)
    call msg("Rules%MassMapOutLst")
     call report_output_list(def%Rules%MassMapOutLst)

    if(def%Rules%cldRepInterv /= 1)then
      call msg('Cloud mass report every n-th output, n=',def%Rules%cldRepInterv)
    else
      call msg('Cloud mass report every output time')
    endif

    call msg('=============================================================')
    call msg(' end of output parameter report')
    call msg('')
    call msg('')

  end subroutine report_output_definition

end module io_server
