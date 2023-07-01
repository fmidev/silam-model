MODULE pollution_cloud
!
! The module contains type silam_pollution_cloud and necessary
! routines and functions for dealing with it
!
! Author: Mikhail Sofiev, mikhail.sofiev@fmi.fi
!
! All units: SI
! 
! All degrees: real numbers in degrees (degrees and parts per
! hundred/parts per thousand , NOT in degrees and minutes
! 
! Language: ANSI standard Fortran 90
! 
! Modules used:

  use lagrange_particles
  use advection_lagrangian
  use advection_eulerian
  use ini_boundary_conditions
  use chemistry_manager
  use source_apportionment
  use optical_density
  use source_terms_general
  use depositions
  !$ use OMP_LIB
  
  IMPLICIT NONE

  private                    ! Cut out all the mess, which is not needed

  ! Public functions
  PUBLIC init_pollution_cloud
  public add_cloud_driver_input_needs
  public init_cloud_internalFields
  public reset_cloud_data_structures
  public set_meteo2disp_interp
  PUBLIC source_to_initial_cloud
  public set_optical_structures
  public set_aux_disp_massmap
  public set_react_rates_map
  public set_cnc2m_map
  public set_dynamic_emission
  public set_cloud_initial_conditions
  public prepare_transdep
  public check_cloud_species
  PUBLIC defined
  public fu_species_transport
  public fu_nbr_of_species_transport
  public fu_species_emission
  public fu_nbr_of_species_emission
  public fu_species_optical
  public fu_nbr_of_species_optical
  public fu_species_short_lived
  public fu_nbr_of_species_short_lived
  public fu_species_aerosol
  public fu_nbr_of_species_aerosol
  public fu_nbr_of_sources ! There can be many sources IN the pollution cloud
  PUBLIC fu_earliest_start
  public fu_emisMM_ptr
  public fu_reactRateMM_ptr
  public fu_concMM_ptr
  public fu_concentration2mMM_ptr
  public fu_optical_densityMM_ptr
  public fu_optical_column_depthMM_ptr
  public fu_advection_moment_X_MM_ptr
  public fu_advection_moment_Y_MM_ptr
  public fu_advection_moment_Z_MM_ptr
  public fu_emission_moment_X_MM_ptr
  public fu_emission_moment_Y_MM_ptr
  public fu_emission_moment_Z_MM_ptr
  public fu_drydepMM_ptr
  public fu_wetdepMM_ptr
  public fu_aerosolMM_ptr
  public fu_shortlivedMM_ptr
  public fu_concLP_ptr
  public fu_aerosolLP_ptr
  public fu_LP_dynamic_params
  public fu_LP_status
  public fu_lpSet_ptr
  public fu_nbr_of_lagr_particles
  public get_cloud_inventory
  public advect_pollution_cloud_v4 ! Lagrangian advection...
  public transform_pollution_cloud_v5 ! Chemical, radioactive, deposition,...
  public fu_if_cloud_mass_map_quantity ! if this is a particle cloud-related or buffer-field quantity
  public fu_subst_nbr_in_inventory
  public fu_if_variable_available    ! from wherever dispersion related
  public check_masses
  public check_mass_centres
  public collect_total_masses
  PUBLIC report
  PUBLIC report_total_masses
  public report_inout_mass_cld
  public fu_ifMeteo2DispHorizInterp
  public fu_ifMeteo2DispVertInterp
  public fu_interpCoefMeteo2DispHoriz
  public fu_interpCoefMeteo2DispVert
  public set_emission_processor_ptr
  public fu_emission_processor_ptr
  public fu_transport_owned_quantity
  public fu_optics_owned_quantity
  public set_3_sources
  public fu_boundarystructures
  public fu_boundaryBuffer
  
  public check_species

  public init_tla


  ! Private functions

  private prepare_transdep_cloud

  PRIVATE print_pollution_cloud_report
  PRIVATE fu_cloud_defined
  private fu_nbr_of_sources_of_cloud  
  PRIVATE fu_cloud_earliest_start

  private check_cloud_masses
  private check_cloud_mass_centres
  private collect_total_masses_in_cloud
  private collect_mass_budget_eulerian
  private collect_mass_budget_lagrangian

  private report_total_masses_of_cloud 

  ! Generic names and operator-interfaces of some functions:

  INTERFACE report
    MODULE PROCEDURE print_pollution_cloud_report
  END INTERFACE

  INTERFACE report_total_masses
    module procedure report_total_masses_of_cloud 
  END INTERFACE

  interface check_masses
    module procedure check_cloud_masses
  end interface

  interface check_mass_centres
    module procedure check_cloud_mass_centres
  end interface
  
  INTERFACE defined
    MODULE PROCEDURE fu_cloud_defined
  END INTERFACE

  INTERFACE fu_earliest_start
    MODULE PROCEDURE fu_cloud_earliest_start
  END INTERFACE

  interface fu_nbr_of_sources
    module procedure fu_nbr_of_sources_of_cloud
  end interface

  interface collect_total_masses
    module procedure collect_total_masses_in_cloud
  end interface

  interface prepare_transdep
    module procedure prepare_transdep_cloud
  end interface
  
  interface fu_species_emission
     module procedure fu_species_emission_cld
  end interface

  !========================================================================
  !========================================================================
  ! 
  ! The main dispersion bunch of memory: pollution cloud.
  ! It includes information for both running the Lagrangian particles and
  ! Eulerian maps. Therefore, three main structures are defined:
  ! - 1-D array of Lagrangian particles,
  ! - cocktail map (4D object that takes care of Eulerian chemistry)
  ! - maps to keep the Eulerian advection information: coordinates of mass centres
  ! In addition, it keeps the out-of-grid cocktail arrays and a few supplementary
  ! data.
  !
  ! There is an redundancy between the Lagrangian particles and Eulerian cocktail
  ! map. Both contain mass vectors for the species. However, this seems to be a
  ! reasonable compromise at this stage because it:
  ! - allows simple addition of Eulerian mechanism without rewriting the Lagrangian one
  ! - allows simple extension to hybrid model runs due to co-existence of both
  !   data structures.
  !
  !========================================================================
  !========================================================================

  TYPE silam_pollution_cloud
    PRIVATE
    !
    ! General part
    !
    TYPE(silja_time) :: earliest_start
    TYPE(silja_time) :: latest_start
    type(silam_species), dimension(:), pointer :: speciesEmission => null(), &
                                                & speciesTransport => null(), &
                                                & speciesShortLived => null(), &
                                                & speciesAerosol => null(), &
                                                & speciesOpticalDensity => null()
    integer :: nSpEmission, nSpTransport, nSpShortLived, nSpAerosol, nSpOpticalDns, nReactions
    !
    ! Lagrangian chemical and dynamic parts.
    ! Note that Lagrangian components can be sent to eulerian grid with possible lumping/splitting
    !
    type(Tlagrange_particles_set) :: lpSet
    !
    ! Eulerian dynamics: 
    ! moments of centres of masses, interpolation strcutures
    !
    type(Tmass_map) :: mapPx_emis =  mass_map_missing, mapPy_emis =  mass_map_missing, mapPz_emis =  mass_map_missing, &
                              & mapPx_conc =  mass_map_missing, mapPy_conc =  mass_map_missing, mapPz_conc =  mass_map_missing
    type(THorizInterpStruct), pointer :: interpCoefMeteo2DispHoriz => null()
    type(TVertInterpStruct), pointer :: interpCoefMeteo2DispVert => null()
    logical :: ifMeteo2DispHorizInterp, ifMeteo2DispVertInterp ! May be, grids & verticals are same?
    !
    ! Eulerian chemical part: 
    ! Eulerian grid cell masses, dry and wet deposition. 
    ! Note that deposition gridds are common with Lagragian advection.
    type(Tmass_map) :: mapEmis =  mass_map_missing, mapConc =  mass_map_missing, mapConc2m =  mass_map_missing, &
                              & mapDryDep =  mass_map_missing, mapWetDep =  mass_map_missing, &
                              & mapOptDns =  mass_map_missing, mapOptColDepth =  mass_map_missing, &
                              & mapShortlived =  mass_map_missing, mapAerosol =  mass_map_missing, &
                              & mapReactRate =  mass_map_missing
    logical :: ifMakeOpticalDensity, ifMakeOpticalColumnDepth, ifMakeCnc2m, ifMakeReactRates
    !
    ! transport outside(x/y 0/M - ix/iy<1 / ix/iy>xMax/yMax) and 
    ! stopped particles / background (garbage) - for each source
    real, dimension(:,:,:,:), pointer :: x0_mass => null(), y0_mass => null(), &
                                       & xM_mass => null(), yM_mass => null() ! (nSrc,nSpecies,2,nz)
    real, dimension(:,:,:), pointer :: mInAir => null()                       ! (nSrc,nSpecies,nz)
    real, dimension(:,:,:), pointer :: vert_mass_0 => null(), vert_mass_M => null()  ! (nSrc,nSpecies,2)
    real, dimension(:,:), pointer :: garbage_mass => null(), mDryDep => null(), mWetDep => null(), &
                                   & mOut => null(), mIn => null(), mStop => null()  ! (nSrc,nSpecies)
    type(Temission_processor), pointer :: emission_processor => null()
    !
    ! the boundary structures. Array elements are for N S E W boundaries respectively.
    !
    type(TboundaryStructPtr), dimension(:), pointer :: boundaryStructures => null()
    type(TboundaryBuffer), pointer :: pBoundaryBuffer => null()
    
    TYPE(silja_logical) :: defined
  END TYPE silam_pollution_cloud
  
  public silam_pollution_cloud

  !
  ! Structure describing the dynamics of the simulations.
  ! The idea is that the simulations can be based on Eulerian, Lagrangian of hybrid 
  ! ideology, with or without exchange between the grid-based and particle-based parts.
  ! Separation will go along the sources: each source can emit into either of the environments.
  ! Dynamics can be independent or connected into a hybrid simulation, In the later case the 
  ! Lagrangian particles are sent to Eulerian environment when their concentration or size
  ! become comparable with those of Eulerian grid.
  !
  type Tdynamics_rules
    integer :: advMethod_Eulerian ! Eulerian advection is split to vert and horiz
    integer :: advection_variant ! advect_rect, advect_tri, advect_step
    real    :: smoother_factor  ! 1 for no smooth, 0 for "upwind"
    integer :: advMethod_Lagrangian ! methods Lagrangian dynamic
    integer :: diffusionMethod         ! diffusion method
    integer :: advectionType_default   ! If not stated, use this dynamics type
    integer :: simulation_type         ! Eulerian / Lagrangian / hybrid
    logical :: ifMolecDiff ! if use gravitationa separation of gases 
                           ! Works only on hybrid vertical and at pressures upt to tens of Pa
    logical :: ifSubgridDiff ! Diffusion affects also vertical CM
    type(TLagr2Euler_projectionRules) :: projectionRules
  end type Tdynamics_rules
  public Tdynamics_rules

  logical, private :: force_3sources = .false.

  !
  ! Species type selector - for generic calls
  !
  integer, parameter, public :: species_emission = 7101
  integer, parameter, public :: species_transport = 7102
  integer, parameter, public :: species_short_lived = 7103
  integer, parameter, public :: species_aerosol = 7104
  integer, parameter, public :: species_optic_density_3D = 7105
  integer, parameter, public :: species_column_optic_depth_2D = 7106


CONTAINS


  ! ***************************************************************

  SUBROUTINE init_pollution_cloud(cloud, em_source, rulesChem, timestart, timestep, timestep_output)
    !
    ! Initializes the cloud chemistry. NOthing else is defined at this early stage.
    !
    IMPLICIT NONE
    ! Imported parameters
    TYPE(silam_pollution_cloud), INTENT(inout) :: cloud
    TYPE(silam_source), INTENT(inout) :: em_source
    type(Tchem_rules), intent(inout) :: rulesChem
    type(silja_time), intent(in) :: timestart
    type(silja_interval), intent(in) :: timestep, timestep_output

    ! Local declarations:
    INTEGER :: iTmp

    !
    ! Some basic stuff
    !
    cloud%earliest_start = time_missing
    cloud%latest_start = time_missing
    nullify(cloud%mInAir)
    nullify(cloud%mDryDep)
    nullify(cloud%mWetDep)
    nullify(cloud%mOut)
    nullify(cloud%mIn)
    nullify(cloud%mStop)
    nullify(cloud%speciesEmission)
    nullify(cloud%speciesTransport)
    nullify(cloud%speciesShortLived)
    nullify(cloud%speciesAerosol)
    nullify(cloud%speciesOpticalDensity)

    allocate(cloud%emission_processor, stat=iTmp)
    if (fu_fails(iTmp == 0, 'Emission processor allocation failed', 'init_pollution_cloud')) return

    cloud%nSpEmission=0; cloud%nSpTransport=0; cloud%nSpShortLived=0; 
    cloud%nSpAerosol=0;    cloud%nSpOpticalDns=0
    cloud%nReactions=0;
    cloud%ifMakeOpticalDensity = .false.
    cloud%ifMakeOpticalColumnDepth = .false.
    cloud%ifMakeCnc2m = .false.
    cloud%ifMakeReactRates = .false.

    !-----------------------------------------------------------------
    ! 
    ! Chemistry. Get the inventory in order to be able to create the maps of masses etc.
    !
    ! Note that we have to it twice - for Lagrnagian and Eulerian environments
    ! First, let's get all what we have
    !
    call get_inventory_g_src(em_source, int_missing, cloud%speciesEmission, cloud%nSpEmission)
    if (fu_fails(cloud%nSpEmission > 0, 'No emission species', 'init_pollution_cloud'))return
    call msg('')
    call msg('Emission species, all environments:', cloud%nSpEmission)
    do iTmp = 1, cloud%nSpEmission
      call report(cloud%speciesEmission(iTmp))
    end do
    call msg('')
    call global_chemical_init(rulesChem, timestep, timestep_output, &        ! input
                            & cloud%speciesEmission, cloud%nSpEmission, &    ! input
                            & cloud%speciesTransport, cloud%nSpTransport, &    ! output
                            & cloud%speciesShortlived, cloud%nSpShortlived, &  ! output
                            & cloud%speciesAerosol, cloud%nSpAerosol, &
                            & cloud%nReactions)          ! output
    if(error)return

    !
    ! Some of the settings may be incompatible with time steps. Verify it here
    !
    call verify_sources(em_source, timestart, timestep)
    if(error)return
    
    cloud%defined = fu_set_true()

  END SUBROUTINE init_pollution_cloud


  !*****************************************************************

  subroutine add_cloud_driver_input_needs(em_source,chemRules,q_met_dyn, q_met_stat, &
                                                            & q_disp_dyn, q_disp_stat)
    !
    ! Returns the list of quantities needed for the pollution_cloud driving
    ! routines. In fact, most of this is a violation of encapsulation rules
    ! but they sometimes can make the program much faster, so I afford them
    ! and just try not to go too far.
    !
    implicit none

    ! imported parameters
    type(silam_source), intent(in) :: em_source
    type(Tchem_rules), intent(in) :: chemRules
    integer, dimension(:), intent(inout) :: q_met_dyn, q_met_stat, q_disp_dyn, q_disp_stat

    ! Local variabels
    integer :: iTmp

    iTmp = fu_merge_integer_to_array(cell_size_x_flag, q_disp_stat)
    iTmp = fu_merge_integer_to_array(cell_size_y_flag, q_disp_stat)
!    iTmp = fu_merge_integer_to_array(cell_size_z_flag, q_disp_dyn)
!    itmp = fu_merge_integer_to_array(air_density_flag, q_disp_dyn)
  end subroutine add_cloud_driver_input_needs


  !**********************************************************************

  subroutine init_cloud_internalFields(dispersionMarketPtr, Rules, wdr)
    !
    ! This subroutine...just does nothing.
    !
    implicit none

    ! Imported parameter
    type(Tchem_rules), intent(in) :: Rules
    type(silja_wdr), intent(in) :: wdr
    type(mini_market_of_stacks), pointer :: dispersionMarketPtr

  end subroutine init_cloud_internalFields


  !*********************************************************************
  
  subroutine reset_cloud_data_structures(cloud, ifMomentsToo)
    !
    ! Zeroes main data storage, prepares the cloud to re-use. Needed for 
    ! switching between forward and adjoint runs in the data assimilation iteration loop
    !
    implicit none

    ! Imported parameters with intent INOUT or POINTER:
    TYPE(silam_pollution_cloud), intent(inout) :: cloud
    logical, intent(in) :: ifMomentsToo
      
    !
    ! Lagrangian environment
    !
    cloud%lpSet = lagrange_particles_set_missing

    ! Need care: the mass maps seem to be usually allocated but set missing. Still handle
    ! the case of not associated maps.
    !
    call reset_if_possible(cloud%mapEmis)
    call reset_if_possible(cloud%mapConc)
    call reset_if_possible(cloud%mapConc2m)
    call reset_if_possible(cloud%mapDryDep)
    call reset_if_possible(cloud%mapWetDep)
    call reset_if_possible(cloud%mapReactRate)
    call reset_if_possible(cloud%mapOptDns)
    call reset_if_possible(cloud%mapOptColDepth)
    call reset_if_possible(cloud%mapShortlived)
    call reset_if_possible(cloud%mapAerosol)

    if(associated(cloud%x0_mass)) cloud%x0_mass = 0.0
    if(associated(cloud%y0_mass)) cloud%y0_mass = 0.0
    if(associated(cloud%xM_mass)) cloud%xM_mass = 0.0
    if(associated(cloud%yM_mass)) cloud%yM_mass = 0.0
    if(associated(cloud%vert_mass_0)) cloud%vert_mass_0 = 0.0
    if(associated(cloud%vert_mass_M)) cloud%vert_mass_M = 0.0
    if(associated(cloud%garbage_mass)) cloud%garbage_mass = 0.0
    if(associated(cloud%mInAir)) cloud%mInAir = 0.0
    if(associated(cloud%mDryDep)) cloud%mDryDep = 0.0
    if(associated(cloud%mWetDep)) cloud%mWetDep = 0.0
    if(associated(cloud%mOut)) cloud%mOut = 0.0
    if(associated(cloud%mIn)) cloud%mIn = 0.0
    if(associated(cloud%mStop)) cloud%mStop = 0.0

    if(associated(cloud%pBoundaryBuffer)) then
       if (associated(cloud%pBoundaryBuffer%bNorth))  cloud%pBoundaryBuffer%bNorth = 0.0
       if (associated(cloud%pBoundaryBuffer%bSouth))  cloud%pBoundaryBuffer%bSouth = 0.0
       if (associated(cloud%pBoundaryBuffer%bEast))   cloud%pBoundaryBuffer%bEast = 0.0
       if (associated(cloud%pBoundaryBuffer%bWest))   cloud%pBoundaryBuffer%bWest = 0.0
       if (associated(cloud%pBoundaryBuffer%bTop))    cloud%pBoundaryBuffer%bTop  = 0.0
       if (associated(cloud%pBoundaryBuffer%bBottom)) cloud%pBoundaryBuffer%bBottom = 0.0
    endif

    if(ifMomentsToo)then
      call reset_if_possible(cloud%mapPx_emis)
      call reset_if_possible(cloud%mapPy_emis)
      call reset_if_possible(cloud%mapPz_emis)
      call reset_if_possible(cloud%mapPx_conc)
      call reset_if_possible(cloud%mapPy_conc)
      call reset_if_possible(cloud%mapPz_conc)
    endif
    
  contains
    
    subroutine reset_if_possible(p_mass_map)
      implicit none
      type(Tmass_map), intent(inout) :: p_mass_map

      if ( defined(p_mass_map)) p_mass_map%arm = 0.0
    end subroutine reset_if_possible

  end subroutine reset_cloud_data_structures
  
  
  !*******************************************************************

  subroutine set_meteo2disp_interp(cloud, wdr, timestep)
    !
    ! Computes, if needed, the structure for interpolation from meteorological 
    ! to dispersion grid
    !
    implicit none

    ! Imported parameters with intent INOUT or POINTER:
    TYPE(silam_pollution_cloud), pointer :: cloud
    type(silja_wdr), pointer :: wdr
    type(silja_interval), intent(in) :: timestep

    ! local variables
    integer :: al_status
    TYPE(silam_pollution_cloud), pointer :: cld

    cld => cloud
    !
    ! Start from meteo - 2- dispersion grid conversion. It is static and has to be
    ! done once
    !
    cloud%interpCoefMeteo2DispHoriz => fu_horiz_interp_struct(meteo_grid, dispersion_grid, &
                                                              & fu_horizontal_interp_method(wdr), &
                                                              & fu_if_randomise(wdr))
    cloud%ifMeteo2DispHorizInterp = .true. ! Now structure has this field
    !
    ! Set the interpolation coefficients between the meteo and disperison verticals
    ! Since vertical interpolation is position-dependent, we have to define the
    ! grid, in which it happens. Dispersion grid is the only way because horizontal interpolation
    ! must go first as an independent procedure. Vertical coefficients depend on grid cell
    ! and thus have to be done after the grid transformation
    !
    if(fu_cmp_verts_eq(meteo_vertical, dispersion_vertical))then
      nullify(cloud%interpCoefMeteo2DispVert)
      cloud%ifMeteo2DispVertInterp = .false.
    else
      cloud%interpCoefMeteo2DispVert => fu_vertical_interp_struct(meteo_vertical, &
                                                                & dispersion_vertical, &
                                                                & dispersion_grid, &
                                                                & fu_vertical_interp_method(wdr), &
                                                                & one_hour, 'main_meteo_to_disp')
!                                                                & timestep, 'main_meteo_to_disp')
      cloud%ifMeteo2DispVertInterp = .true.
    endif

  end subroutine set_meteo2disp_interp


  ! ***************************************************************

  SUBROUTINE source_to_initial_cloud(source, cloud, &
                                   & chemRules, dynRules, &
                                   & calc_start, calc_dur, timestep, &
                                   & meteoMarket, &
                                   & pOutputGrid, &
                                   & iAccuracy, ifRandomise)
    !
    ! Creates the storage for all masses and moments in pollution cloud to be used 
    ! in dynamics, chemistry, output, etc. Their shapes and content differ for Eulerian 
    ! and Lagrangian dynamics.
    !
    ! If Eulerian advection:
    ! Fields are represented as maps in dispersion grid. 
    !
    ! If Lagrangian advection: 
    ! Particles are stored in the cloud  as  lpset member
    !
    IMPLICIT NONE

    ! Imported parameters with intent IN:
    TYPE(silam_source), INTENT(inout) :: source
    TYPE(silja_time), INTENT(in) :: calc_start
    TYPE(silja_interval), INTENT(in) :: calc_dur, timestep
    type(Tchem_rules), intent(inout) :: chemRules
    type(Tdynamics_rules), intent(in) :: dynRules
    type(mini_market_of_stacks), intent(in) :: meteoMarket
    type(silja_grid), pointer :: pOutputGrid
    integer, intent(in) :: iAccuracy
    logical, intent(in) :: ifRandomise

    ! Imported parameters with intent INOUT or POINTER:
    TYPE(silam_pollution_cloud), intent(inout) :: cloud

    ! Local variables
    integer :: iTmp, al_status, nx, ny, nz, nSrcIds, nMomentSpecies
    type(silam_vertical), target :: verticalTmp
    type(silam_vertical), pointer :: verticalPtr
    type(silja_time) :: srcStartTime, srcEndTime
    integer, dimension(:), pointer :: iPassengers

    IF(fu_fails(defined(source),'Undefined source','source_to_initial_cloud'))RETURN

    iPassengers => fu_work_int_array()
    if(error)return
    !
    ! Set time parameters. We have 2 periods: source_start with source_duration,
    ! calculations_start with forecast_duration. Their common period 
    ! determines the earliest and latest start of particles
    !    
    srcStartTime = fu_start_time(source)
    srcEndTime = fu_end_time(source)

    IF(calc_start < calc_start + calc_dur) THEN 
      !
      !------------------------------------------------- FORWARD calculations
      !
      IF(calc_start < srcStartTime)THEN
        cloud%earliest_start = srcStartTime
      ELSE
        cloud%earliest_start = calc_start
      END IF
      IF(calc_start+calc_dur < srcEndTime)THEN
        cloud%latest_start = calc_start+calc_dur
      ELSE
        cloud%latest_start = srcEndTime
      END IF
    ELSE  
      !------------------------------------------------ INVERSE calculations
      !
      ! Remember: source has always duration >0
      !
      IF(calc_start+calc_dur < srcStartTime)THEN
        cloud%earliest_start = srcStartTime
      ELSE
        cloud%earliest_start = calc_start+calc_dur
      END IF
      IF(calc_start < srcEndTime)THEN
        cloud%latest_start = calc_start
      ELSE
        cloud%latest_start = srcEndTime
      END IF
    END IF
    !
    ! Stupidity check
    !
    if(cloud%latest_start < cloud%earliest_start)then
      call msg('Earliest start:' + fu_str(cloud%earliest_start))
      call msg('Latest start:' + fu_str(cloud%latest_start))
      call set_error('Time inconsistency in the below cloud','source_to_initial_cloud')
      call report(cloud)
      return
    endif

    nSrcIds = fu_NbrOf_source_ids(source)

    call msg('Number of cloud species: ', cloud%nSpTransport)
    call msg('Number of source IDs: ', nSrcIds)

    cloud%ifMakeOpticalDensity = .false.  ! Will be reset later if needed
    cloud%ifMakeOpticalColumnDepth = .false.  ! Will be reset later if needed
    cloud%ifMakeCnc2m = .false.
    cloud%ifMakeReactRates = .false.

    if (force_3sources) then
      nSrcIds = DA_NUM_SOURCES
    else
      nSrcIds = fu_NbrOf_source_ids(source)
    endif

    !-------------------------------------------------------------------------------
    !
    ! Create the pointer array to boundary structs. Whether it will be used decided later
    !
    allocate(cloud%boundaryStructures(6), cloud%pBoundaryBuffer, stat=iTmp)
    if(fu_fails(iTmp==0,'Failed allocation of boundary structures', 'source_to_initial_cloud'))return
    do iTmp = 1,6
      nullify(cloud%boundaryStructures(iTmp)%ptr)
    end do

    !-------------------------------------------------------------------
    !
    ! If any of the species is non-passive, there will be deposition grids
    ! Note the single-level vertical
    !
    call set_vertical(ground_level, verticalTmp) ! for deposition fields
    if(error)return
    verticalPtr => verticalTmp

    !-----------------------------------------------------------------------------
    !
    ! Concentration maps and the corresponding moments are made differently
    ! Eulerian and Lagrangian dynamics.
    !
    ! Advection-type specific structures: secondary grid of emission, deposition maps, etc.
    !
    if(fu_if_Eulerian_present(dynRules%simulation_type))then  !=========================== EULERIAN present

      !
      ! Linking the sources to the selected dispersion grid and vertical 
      ! This should return some memory allocated for temporaries in area sources
      ! which is a good to do before allocating huge massmaps...
      call msg('Linking Eulerian emission to dispersion grid and vertical, memusage kB', fu_system_mem_usage() )
      call link_emission_2_dispersion(source, cloud%speciesEmission, &
                                    & meteoMarket, meteo_grid, dispersion_grid, meteo_vertical, dispersion_vertical, &
                                    & iAccuracy, ifRandomise)

      call msg('Making Eulerian emission map, memusage kB', fu_system_mem_usage() )
      !
      ! Moments are defined using the basic parameters only, so we create them here using the 
      ! above dimensions
      !
      if (fu_if_bulk_eulerian_advection(dynRules%advMethod_Eulerian)) then
        nMomentSpecies = 1
      else
        nMomentSpecies = cloud%nSpEmission
      end if
      iPassengers(1:nMomentSpecies) = 0
      !
      ! Now we can create the emission mass map. Note that this mass map does not have any passengers
      !
      call set_aux_disp_massmap(cloud%mapEmis, nSrcIds,  cloud%speciesEmission,  &
              &cloud%nSpEmission, emission_intensity_flag)

      call set_mass_map_from_basic_par(cloud%mapPx_emis, advection_moment_X_flag, &    ! quantity
                  & nx_dispersion, ny_dispersion, &
                  & nz_dispersion, nSrcIds, & ! fu_NbrOf_source_ids(source), &
                  & nMomentSpecies, &         ! nSpecies
                  & iPassengers, &
                  & 0, 0, 0, 0, &             ! x, y, z, source border cells
!                  & 1, 1, 1, 0, &             ! x, y, z, source border cells
                  & val=0.)                   ! value
      
      call set_mass_map_from_basic_par(cloud%mapPy_emis, advection_moment_Y_flag, &    ! quantity
                  & nx_dispersion, ny_dispersion, &
                  & nz_dispersion, nSrcIds, & ! fu_NbrOf_source_ids(source), &
                  & nMomentSpecies, &         ! nSpecies
                  & iPassengers, &
                   & 0, 0, 0, 0, &             ! x, y, z, source border cells
!                 & 1, 1, 1, 0, &             ! x, y, z, source border cells
                  & val=0.)                   ! value

      call set_mass_map_from_basic_par(cloud%mapPz_emis, advection_moment_Z_flag, &    ! quantity
                                        & nx_dispersion, ny_dispersion, &
                                        & nz_dispersion, nSrcIds, & ! fu_NbrOf_source_ids(source), &
                                        & nMomentSpecies, &         ! nSpecies
                                        & iPassengers, &
                                         & 0, 0, 0, 0, &             ! x, y, z, source border cells
!                                       & 1, 1, 1, 0, &             ! x, y, z, source border cells
                                        & val=0.)                   ! value
      if(error)return
      cloud%mapPx_emis%gridTemplate = dispersion_gridPtr
      cloud%mapPx_emis%vertTemplate = dispersion_verticalPtr
      cloud%mapPy_emis%gridTemplate = dispersion_gridPtr
      cloud%mapPy_emis%vertTemplate = dispersion_verticalPtr
      cloud%mapPz_emis%gridTemplate = dispersion_gridPtr
      cloud%mapPz_emis%vertTemplate = dispersion_verticalPtr
      if(nMomentSpecies == cloud%nSpEmission) then
        call set_map_species(cloud%mapPy_emis, cloud%speciesEmission)
        call set_map_species(cloud%mapPx_emis, cloud%speciesEmission)
        call set_map_species(cloud%mapPz_emis, cloud%speciesEmission)
      endif

      call msg('Making Eulerian concentration map, memusage kB', fu_system_mem_usage() )
      call set_aux_disp_massmap(cloud%mapConc, nSrcIds,  cloud%speciesTransport, &
            & cloud%nSpTransport, mass_in_air_flag)

      if(cloud%nSpShortLived > 0)then
        call msg('Making shortlived map, memusage kB', fu_system_mem_usage() )
        call set_aux_disp_massmap(cloud%mapShortLived, nSrcIds,  cloud%speciesShortLived, &
            & cloud%nSpShortLived, mass_in_air_flag)
      endif
      !
      ! And this strange creature for aerosol dynamics (flipping between diameters and particle numbers)
      !
      if (cloud%nSpAerosol > 0) then
        call msg('Making aerosol map, memusage kB', fu_system_mem_usage() )
        call set_aux_disp_massmap(cloud%mapAerosol, nSrcIds,  cloud%speciesAerosol, &
                                            & cloud%nSpAerosol, aerosol_flag)
      end if

      !
      ! Moments are defined using the basic parameters only, so we create them here using the 
      ! above dimensions
      !
      if(fu_if_bulk_eulerian_advection(dynRules%advMethod_Eulerian)) then
        nMomentSpecies = 1
      else
        nMomentSpecies = cloud%nSpTransport
      end if
      iPassengers(1:nMomentSpecies) = 0
      call msg('Making X moment map, memusage kB', fu_system_mem_usage() )
      call set_mass_map_from_basic_par(cloud%mapPx_conc, advection_moment_X_flag, &  ! quantity
                 & nx_dispersion, ny_dispersion, nz_dispersion, nSrcIds, &
                 & nMomentSpecies, &          ! nSpecies
                 & iPassengers, &
                 & 0, 0, 0, 0, &              ! x, y, z, source border cells
!                 & 1, 1, 1, 0, &              ! x, y, z, source border cells
                 & val=0.)                    ! value
      call msg('Making Y moment map, memusage kB', fu_system_mem_usage() )
      call set_mass_map_from_basic_par(cloud%mapPy_conc, advection_moment_Y_flag, &  ! quantity
                 & nx_dispersion, ny_dispersion, nz_dispersion, nSrcIds, &
                 & nMomentSpecies, &          ! nSpecies
                 & iPassengers, &
                 & 0, 0, 0, 0, &              ! x, y, z, source border cells
!                 & 1, 1, 1, 0, &              ! x, y, z, source border cells
                 & val=0.)                    ! value
      call msg('Making Z moment map, memusage kB', fu_system_mem_usage() )
      call set_mass_map_from_basic_par(cloud%mapPz_conc, advection_moment_Z_flag, &  ! quantity
                                        & nx_dispersion, ny_dispersion, nz_dispersion, nSrcIds, &
                                        & nMomentSpecies, &          ! nSpecies
                                        & iPassengers, &
                                        & 0, 0, 0, 0, &              ! x, y, z, source border cells
!                                        & 1, 1, 1, 0, &              ! x, y, z, source border cells
                                        & val=0.)                    ! value
      cloud%mapPx_conc%gridTemplate = dispersion_gridPtr
      cloud%mapPx_conc%vertTemplate = dispersion_verticalPtr
      cloud%mapPy_conc%gridTemplate = dispersion_gridPtr
      cloud%mapPy_conc%vertTemplate = dispersion_verticalPtr
      cloud%mapPz_conc%gridTemplate = dispersion_gridPtr
      cloud%mapPz_conc%vertTemplate = dispersion_verticalPtr
      if(nMomentSpecies == cloud%nSpTransport) then
        call set_map_species(cloud%mapPx_conc, cloud%speciesTransport)
        call set_map_species(cloud%mapPy_conc, cloud%speciesTransport)
        call set_map_species(cloud%mapPz_conc, cloud%speciesTransport)
      endif
      call msg('Eulerian massmaps done, memusage kB', fu_system_mem_usage() )

    endif  ! if Eulerian present
    
    
    if(fu_if_Lagrangian_present(dynRules%simulation_type))then            !================== LAGRANGIAN present
      !
      call init_lagrange_particles_set(cloud%lpSet, cloud%speciesTransport,  &
                  & cloud%speciesShortLived , cloud%speciesAerosol, meteo_grid,&
                  & meteo_vertical, nSrcIds)

      ! The number of particles. Let the neutral accuracy be 100k, then follow Eulerian algorithm: 
      ! */ by a factor of 200 for the 0-10 range.
      !
      call  enlarge_lagrange_particles_set(cloud%lpset, fu_number_of_lagr_particles(iAccuracy))

      !
      ! Linking the sources to the Lagrangian dispersion grid and vertical. Note that there will be no
      ! emission species
      !
      call msg('Linking Lagrangian emission to dispersion grid and vertical, memusage kB', fu_system_mem_usage() )
      call link_emission_2_dispersion(source, cloud%speciesEmission, meteoMarket, &
                                    & meteo_grid, dispersion_grid, meteo_vertical, dispersion_vertical, &
                                    & iAccuracy, ifRandomise)
    else  ! No Lagrangian dynamics
      call set_lpset_missing(cloud%lpSet)
    endif  ! If Lagrangian present
    if(error)return

    !
    ! Make deposition maps. They are the same for both Lagrangian and Eulerian environments, 
    ! made in dispersion_grid.
    !
    call msg('Making dry deposition, memusage kB', fu_system_mem_usage() )
    call set_aux_disp_massmap(cloud%mapDryDep, nSrcIds,  cloud%speciesTransport, &
                                            & cloud%nSpTransport, drydep_flag)
    call msg('Making wet deposition, memusage kB', fu_system_mem_usage() )
    call set_aux_disp_massmap(cloud%mapWetDep, nSrcIds,  cloud%speciesTransport, &
                                            & cloud%nSpTransport, wetdep_flag)

    if(error)return
    
    !------------------------------------------------------------------------------------
    !
    ! Now, a few general things that are valid for both types of dynamics
    !
    if (cloud%nSpAerosol > 0) then
      !
      ! The aerosol species must be linked to the transported ones to handle the flip between
      ! the particle sizes and their mass-weighted moments (actually, between the single-particle 
      ! volume and number concentrations)
      !
      call msg('Allocating chemRules%ChemRunSetup%mapVolume2NumberCnc, memusage kB', fu_system_mem_usage() )
      allocate(chemRules%ChemRunSetup%mapVolume2NumberCnc)
      call create_mass_2_nbr_mapping(cloud%mapConc, &      ! transport: 
                                   & cloud%mapAerosol, &   ! aerosol
                                   ! mapping of number- and mass- species:
                                   & chemRules%ChemRunSetup%mapVolume2NumberCnc) 
      if(error)return
    else
      nullify(chemRules%ChemRunSetup%mapVolume2NumberCnc)
    endif

    !---------------------------------------------------------------------
    !
    ! Initalise arrays for total masses outside and garbage stuff. 
    !
    if(fu_if_Eulerian_present(dynRules%simulation_type))then  
      nz = nz_dispersion
    else
      nz = 1  ! Lagrangian does not have split over levels
    endif
    call msg('Allocating accounting, memusage kB', fu_system_mem_usage() )
    allocate(cloud%mInAir (nSrcIds,cloud%nSpTransport,nz), &
           & cloud%mDryDep(nSrcIds,cloud%nSpTransport), &
           & cloud%mWetDep(nSrcIds,cloud%nSpTransport), &
           & cloud%mOut   (nSrcIds,cloud%nSpTransport), &
           & cloud%mIn   (nSrcIds,cloud%nSpTransport), &
           & cloud%mStop  (nSrcIds,cloud%nSpTransport), &
           & cloud%x0_mass(nSrcIds,cloud%nSpTransport,2,nz), &
           & cloud%y0_mass(nSrcIds,cloud%nSpTransport,2,nz), &
           & cloud%xM_mass(nSrcIds,cloud%nSpTransport,2,nz), &
           & cloud%yM_mass(nSrcIds,cloud%nSpTransport,2,nz), &
           & cloud%vert_mass_0 (nSrcIds,cloud%nSpTransport,2), &
           & cloud%vert_mass_M(nSrcIds,cloud%nSpTransport,2), &
           & cloud%garbage_mass(nSrcIds,cloud%nSpTransport), stat=al_status)
    if(fu_fails(al_status == 0, 'Cannot ALLOCATE total masses space','source_to_initial_cloud'))return
 
    cloud%mInAir(1:nSrcIds,1:cloud%nSpTransport,1:nz) = 0.
    cloud%mDryDep(1:nSrcIds,1:cloud%nSpTransport) =0.
    cloud%mWetDep(1:nSrcIds,1:cloud%nSpTransport) =0.
    cloud%mOut(1:nSrcIds,1:cloud%nSpTransport) = 0.
    cloud%mIn(1:nSrcIds,1:cloud%nSpTransport) = 0.
    cloud%mStop(1:nSrcIds,1:cloud%nSpTransport) =0.
    cloud%x0_mass(1:nSrcIds,1:cloud%nSpTransport,1:2,1:nz) =0.
    cloud%y0_mass(1:nSrcIds,1:cloud%nSpTransport,1:2,1:nz) = 0.
    cloud%xM_mass(1:nSrcIds,1:cloud%nSpTransport,1:2,1:nz) = 0.
    cloud%yM_mass(1:nSrcIds,1:cloud%nSpTransport,1:2,1:nz) = 0.
    cloud%vert_mass_0(1:nSrcIds,1:cloud%nSpTransport,1:2) = 0.
    cloud%vert_mass_M(1:nSrcIds,1:cloud%nSpTransport,1:2) = 0.
    cloud%garbage_mass(1:nSrcIds,1:cloud%nSpTransport) =0.

    if(.not.error) call msg('Computational environment in pollution cloud is created, memusage kB', fu_system_mem_usage() )

!        call msg('')
!        call msg('Setting trial initial mass')
!        cloud%mapConc%arM(1:cloud%mapConc%nSpecies, &
!                        & 1:cloud%mapConc%nSrc, &
!                        & 1:cloud%mapConc%n3d, &
!                        & 1:cloud%mapConc%nx, &
!                        & 1:cloud%mapConc%ny) = 1.0
!        call msg('Initial mass set')
!        call msg('')

!        call set_advection_moment_mapper(cloud%mapP_emis%mapper, 1, 2, 3)
!        if(error)return

  END SUBROUTINE source_to_initial_cloud


  !*********************************************************************


  subroutine set_aux_disp_massmap(map, nSrc,  pSpecies, nSpecies, iQuantity)
    !
    ! Creates an auxillary massmap
    !
    implicit none

    ! Imported parameters
    type(Tmass_Map), intent(out) :: map
    type(silam_species), dimension(:), intent(in) :: pSpecies

    integer, intent(in) :: nSrc, iQuantity, nSpecies

    ! Local variables
    type(silam_vertical), target :: verticalTmp
    type(silam_vertical), pointer :: verticalPtr


    if (fu_multi_level_quantity(iQuantity)) then
      call msg('Making a multilevel map for '//trim(fu_quantity_string(iQuantity)))
      verticalPtr => dispersion_verticalPtr
    else
      ! and 2D 
       call msg('Making a single-level map for '//trim(fu_quantity_string(iQuantity)))
      call set_vertical(ground_level, verticalTmp) ! for deposition fields
      if(error)return
      verticalPtr => verticalTmp
    endif
    call set_mass_map(map, iQuantity, & 
                   & nSrc, &
                   & 0, &                    ! border cells
                   & dispersion_grid, &
                   & verticalPtr, &
                   & pSpecies(1:nSpecies), & !cocktailPtr%species_list, &
                   & val=0.)   ! value

  end subroutine set_aux_disp_massmap

  !******************************************************************************

  subroutine set_cnc2m_map(cloud)
    !
    ! Creates the _cnc2m_map in cloud
    !
    implicit none

    ! Imported parameters
    TYPE(silam_pollution_cloud), target, INTENT(inout) :: cloud

    call msg("Creating cnc2m  map for "//trim(fu_str(cloud%nSpTransport))//" species")
    call set_aux_disp_massmap(cloud%mapConc2m, cloud%mapConc%nSrc, cloud%speciesTransport, cloud%nSpTransport, concentration_2m_flag)
    cloud%ifMakeCnc2m = .true.

  end subroutine set_cnc2m_map

  !******************************************************************************

  subroutine set_react_rates_map(cloud)
    !
    ! Creates the map for reaction rates
    !
    implicit none

    ! Imported parameters
    TYPE(silam_pollution_cloud), target, INTENT(inout) :: cloud
    type(silam_species), dimension(max_species) :: miss_species

    miss_species(1:cloud%nReactions) = species_missing
    call msg("Creating reaction rates map for "//trim(fu_str(cloud%nReactions))//" reactions")

    call set_aux_disp_massmap(cloud%mapReactRate, cloud%mapConc%nSrc, miss_species, cloud%nReactions, reaction_rate_flag)
    cloud%ifMakeReactRates = .true.

  end subroutine set_react_rates_map

  !******************************************************************************

  subroutine set_optical_structures(cloud, rulesOptical, ChemRunSetup, &
                                  & optSpecies, noptSpecies, nlStdSetup, iQuantity)
    !
    ! Creates the optical cocktail in the cloud
    !
    implicit none

    ! Imported parameters
    TYPE(silam_pollution_cloud), target, INTENT(inout) :: cloud
    type(Toptical_density_rules), intent(inout) :: rulesOptical
    type(silam_species), dimension(:), intent(in) :: optSpecies
    type(TChemicalRunSetup), pointer :: ChemRunSetup
    type(Tsilam_namelist), pointer :: nlStdSetup
    type(silam_vertical), target :: verticalTmp
    type(silam_vertical), pointer :: verticalPtr
    integer, intent(in) :: iQuantity, noptSpecies

    ! Local variables
    integer :: allocstat
    !
    ! Actually, just calls for chemistry_server routine giving it optical rules,
    ! transport species as template, optical metadata for copying and optical species themselves
    !
    if(.not. (cloud%ifMakeOpticalDensity .or. cloud%ifMakeOpticalColumnDepth))then
    
      call addSpecies(cloud%speciesOpticalDensity, cloud%nSpOpticalDns, optSpecies, noptSpecies)

      if(error)return
      !
      ! Establish the link between the transport and optical species
      !
      allocate(chemRunSetup%refsTransp2opt(cloud%nSpTransport), &
             & chemRunSetup%refsOpt2transp(cloud%nSpOpticalDns), &
             & stat=allocstat)
      if(fu_fails(allocstat  == 0, 'Allocate failed', 'set_optical_structures'))return

      call set_speciesReferences(cloud%SpeciesTransport, cloud%SpeciesOpticalDensity, &
                               & chemRunSetup%refsTransp2Opt)
      call set_speciesReferences(cloud%SpeciesOpticalDensity, cloud%SpeciesTransport, &
                               & chemRunSetup%refsOpt2Transp)

      if (error) return

!!$      call link_transp_and_opt_species(cloud%cocktailTransp, cloud%cocktailOpticDns, &
!!$                                     & rulesOptical, ChemRunSetup)
      if(error)return
    endif
    !
    ! With cocktail initialised, can create a optical massMap structure
    !
!    cocktailPtr => cloud%cocktailOpticDns
    if(iQuantity == optical_density_flag)then
      !
      ! 3D optical density...
      !
      call set_aux_disp_massmap(cloud%mapOptDns, cloud%mapConc%nSrc, cloud%speciesOpticalDensity, cloud%nSpOpticalDns, iQuantity)
      if(error)return
      cloud%ifMakeOpticalDensity = .true.

    elseif(iQuantity == optical_column_depth_flag)then

      call set_aux_disp_massmap(cloud%mapOptColDepth, cloud%mapConc%nSrc,  cloud%speciesOpticalDensity, cloud%nSpOpticalDns, iQuantity)
      if(error)return
      cloud%ifMakeOpticalColumnDepth = .true.

    else
      call msg('Unknown quantity to be initialised',iQuantity)
      call set_error('Unknown quantity to be initialised','set_optical_structures')
      return
    endif
    !
    ! Having the cocktail and mass map created, we can order precomputations of the 
    ! mass-to-optical-density scaling coefficients for the specific set of substances
    ! and wave lengths
    !
    if((.not.cloud%ifMakeOpticalDensity) .or. (.not.cloud%ifMakeOpticalColumnDepth))then
      call init_optical_density_data(cloud%speciesOpticalDensity, cloud%nSpOpticalDns, &
                                   & nlStdSetup, ChemRunSetup, rulesOptical)
      if(error)return
    endif

  end subroutine set_optical_structures


  !******************************************************************************

  subroutine set_dynamic_emission(cld, &
                                & em_source, &
                                & met_buf, &
                                & disp_buf, &
                                & chemRules, dynRules, &
                                & now, &
                                & timestep, residence_interval, &
                                & fInjectedMass, dispersionMarketPtr)
    !
    ! Computes the now-time emission rates and translates them to masses of 
    ! lagrangian particles which are about to start during this timestep.
    ! For many cases it is void but for e.g. pollen emission it has all the means.
    ! 
    ! The basic idea is that the particle may be started with unknown mass,
    ! and this very sub sets the actual mass with which the particle starts
    ! Evidently, after start the mass of the particle follows own rules
    !
    implicit none

    ! Imported parameters
    TYPE(silam_pollution_cloud), INTENT(inout) :: cld
    TYPE(Tfield_buffer), POINTER :: met_buf, disp_buf
    TYPE(silam_source), INTENT(inout) :: em_source
    type(Tchem_rules), intent(in) :: chemRules
    type(TDynamics_rules), intent(in) :: dynRules
    type(silja_time), intent(in) :: now
    type(silja_interval), intent(in) :: timestep, residence_interval
    real(r8k), dimension(:), INTENT(inout) :: fInjectedMass
    type(mini_market_of_stacks), pointer :: dispersionMarketPtr

    ! Local variables
    logical :: ifDynamic ! Whether emission is dynamic at all
    type(Tmass_map), pointer :: emisMassPtr ! Emision fields in dispersion_grid
    integer :: iP, ix, iy, iz, iSrc, iVar
    type(silja_time) :: t1, t2
    type(silam_grid_position) :: posTmp
    real :: timestepMdlSec_1, fTmp1, fTmp2, fM1, fM2
    
    !
    ! Actual work with emission fields depends on the type of simulation. There are four such:
    ! (i) data assimilation to emission, (ii) eulerian, (iii) lagrangian, (iv) hybrid
    !
    ! Eulerian evironment:
    ! Upon concluding the emission mass map, it is injected into the transport mass map.
    ! At this point the Eulerian advection means that the mass and moments maps are summed-up.
    !
    ! Lagrangian environment:
    ! Lagrangian particles are started from each source. Hence, no emission map is involved.
    !
!call msg('Set dynamic emission: checking the mass centres 1...')
!call check_mass_centres(cld, 'before_set_dynamic_emission')

    if(defined(cld%emission_processor))then
      call set_dynamic_emission_da()
    else
      select case(dynRules%simulation_type)
        case(eulerian_flag)
          call set_dynamic_emission_std_eul() 
        case(lagrangian_flag)
          call set_dynamic_emission_std_lagr() 
        case(hybrid_flag)
          call set_dynamic_emission_std_eul() 
          call set_dynamic_emission_std_lagr() 
        case default
          call msg('Unknown simulation type:',dynRules%simulation_type)
          call set_error('Unknown simulation type','set_dynamic_emission')
      end select
    endif  ! if data assimilation is defined

    
!call msg('Set dynamic emission: checking the mass centres 2...')
!call check_mass_centres(cld, 'after_set_dynamic_emission')
    
    
  CONTAINS
    
    !=======================================================================
    
    subroutine set_dynamic_emission_std_eul()
      implicit none
      
      call start_count('reset_emission_massmap')
      call reset_map_to_val(cld%mapEmis, 0.0)
      call reset_map_to_val(cld%mapPx_emis, 0.0)
      call reset_map_to_val(cld%mapPy_emis, 0.0)
      call reset_map_to_val(cld%mapPz_emis, 0.0)
      call stop_count('reset_emission_massmap')
      if (error) return

      ! Inject the emission to emission mass map
      !
      call inject_emission_eulerian(em_source, &
                                  & cld%mapEmis, cld%mapPx_emis, cld%mapPy_emis, cld%mapPz_emis, &
                                  & met_buf, disp_buf, now, timestep, &
                                  & .not. fu_if_bulk_eulerian_advection(dynRules%advMethod_Eulerian),&
                                  & fInjectedMass, &
                                  & cld%interpCoefMeteo2DispHoriz, cld%ifMeteo2DispHorizInterp, &
                                  & cld%interpCoefMeteo2DispVert, cld%ifMeteo2DispVertInterp)
      if(error)return


!call report(cld%mapEmis)

      if (fu_if_bulk_eulerian_advection(dynRules%advMethod_Eulerian)) then
        call emis2transpMap_Euler(cld%mapEmis, cld%mapPx_emis, cld%mapPy_emis, cld%mapPz_emis, &
                                & cld%mapConc, cld%mapPx_conc, cld%mapPy_conc, cld%mapPz_conc, &
                                & cld%mapAerosol, &
                                & fu_ChemRunSetup(chemRules))
      else
        call emis2transpMap_Euler_speciesMmt(cld%mapEmis, cld%mapPx_emis, &
                                           & cld%mapPy_emis, cld%mapPz_emis, &
                                           & cld%mapConc, cld%mapPx_conc, &
                                           & cld%mapPy_conc, cld%mapPz_conc, &
                                           & cld%mapAerosol, &
                                           & fu_ChemRunSetup(chemRules))
      end if

!call report(cld%mapConc)

    end subroutine set_dynamic_emission_std_eul
    
    !================================================================================
    
    subroutine set_dynamic_emission_da()
      !
      ! If the emission_processor is defined, it takes over: in forward runs, 
      !
      implicit none
      if (fu_interval_positive(timestep)) then
        call reset_map_to_val(cld%mapEmis, 0.0)
        call reset_map_to_val(cld%mapPx_emis, 0.0)
        call reset_map_to_val(cld%mapPy_emis, 0.0)
        call reset_map_to_val(cld%mapPz_emis, 0.0)
        !
        ! Forward integration: emission is injected normally and processed if applicable.
        !
        call inject_emission_eulerian(em_source, &
                                    & cld%mapEmis, cld%mapPx_emis, cld%mapPy_emis, cld%mapPz_emis, &
                                    & met_buf, disp_buf, now, timestep, &
                                    & .not. fu_if_bulk_eulerian_advection(dynRules%advMethod_Eulerian),&
                                    & fInjectedMass, &
                                    & cld%interpCoefMeteo2DispHoriz, cld%ifMeteo2DispHorizInterp, &
                                    & cld%interpCoefMeteo2DispVert, cld%ifMeteo2DispVertInterp)
        if (defined(cld%emission_processor)) then
          call process_emission_forward(cld%emission_processor,&
                                      & cld%mapEmis, cld%mapPx_emis, cld%mapPy_emis, cld%mapPz_emis, now, timestep)
        end if
        if (fu_if_bulk_eulerian_advection(dynRules%advMethod_Eulerian)) then
          call emis2transpMap_Euler(cld%mapEmis, cld%mapPx_emis, cld%mapPy_emis, cld%mapPz_emis, &
                                  & cld%mapConc, cld%mapPx_conc, cld%mapPy_conc, cld%mapPz_conc, &
                                  & cld%mapAerosol, &
                                  & fu_ChemRunSetup(chemRules))
        else
          call emis2transpMap_Euler_speciesMmt(cld%mapEmis, cld%mapPx_emis, &
                                             & cld%mapPy_emis, cld%mapPz_emis, &
                                             & cld%mapConc, cld%mapPx_conc, &
                                             & cld%mapPy_conc, cld%mapPz_conc, &
                                             & cld%mapAerosol, &
                                             & fu_ChemRunSetup(chemRules))
        end if
        
      else
        ! 
        ! Adjoint integration. Emission field needed only if it is to be adjusted.
        !
        call reset_map_to_val(cld%mapEmis, 0.0)
        if (fu_type(cld%emission_processor) /= processor_void) then
          call inject_emission_eulerian(em_source, &
                                      & cld%mapEmis, cld%mapPx_emis, cld%mapPy_emis, cld%mapPz_emis, &
                                      & met_buf, disp_buf, now, timestep, &
                                      & .not. fu_if_bulk_eulerian_advection(dynRules%advMethod_Eulerian),&
                                      & fInjectedMass, &
                                      & cld%interpCoefMeteo2DispHoriz, cld%ifMeteo2DispHorizInterp, &
                                      & cld%interpCoefMeteo2DispVert, cld%ifMeteo2DispVertInterp)
          call process_emission_adjoint(cld%mapEmis, cld%mapConc, cld%emission_processor, now, timestep)
        end if
        call reset_map_to_val(cld%mapEmis, 0.0)
        call reset_map_to_val(cld%mapPx_emis, 0.0)
        call reset_map_to_val(cld%mapPy_emis, 0.0)
        call reset_map_to_val(cld%mapPz_emis, 0.0)
        !cld%mapEmis%arm = 0.0
      end if
      
    end subroutine set_dynamic_emission_da
    
    !=======================================================================================
    
    subroutine set_dynamic_emission_std_lagr()
      implicit none
      !
      ! Inject the emission to emission mass map
      !
      call inject_emission_lagrangian(em_source, &
                                    & cld%lpSet, &
                                    & chemRules%mass_lagr_particle, &
                                    & fu_ChemRunSetup(chemRules),  &
                                    & met_buf, disp_buf, now, timestep, &
                                    & fInjectedMass, &
                                    & cld%interpCoefMeteo2DispHoriz, cld%ifMeteo2DispHorizInterp, &
                                    & cld%interpCoefMeteo2DispVert, cld%ifMeteo2DispVertInterp)
      if(error)return
      !
      ! If new particles were created, may be, we have now more active ones than before
      !
      if(cld%lpSet%nop < cld%lpSet%iFirstEmptyParticle)then
        cld%lpSet%nop = cld%lpSet%iFirstEmptyParticle - 1
      endif
      !
      ! For hybrid simulations, should project "old" lagrangian particles to eulerian field.
      ! Note that in hybrid type of the run Lagrangian particles fly in Eulerian dispersion grid
      !
      if(dynRules%simulation_type == hybrid_flag)then
        call project_lagr_2_euler_flds(cld%lpSet, &
                                     & cld%mapConc, cld%mapAerosol, &
                                     & cld%mapPx_conc, cld%mapPy_conc, cld%mapPz_conc, &
                                     & dynRules%projectionRules, fu_sec(residence_interval))
      endif

    end subroutine set_dynamic_emission_std_lagr

  end subroutine set_dynamic_emission


  !**********************************************************************************************
  
  subroutine set_cloud_initial_conditions(cloud, disp_buf, &
                                        & rulesIniBoundary, now, timestep, &
                                        & meteoMarketPtr, dispersionMarketPtr)
    !
    ! Sets the initial conditions via setting the values of maps of TMassMap object
    ! and, if needed, dispersion buffer variables. 
    ! Must be called after all these structures are properly initialized
    !
    implicit none

    ! Imported parameters
    TYPE(silam_pollution_cloud), INTENT(inout) :: cloud
    type(Tini_boundary_rules), intent(in) :: rulesIniBoundary
    type(Tfield_buffer), pointer :: disp_buf
    type(silja_time), intent(in) :: now
    type(silja_interval), intent(in) :: timestep
    type(mini_market_of_stacks), pointer :: meteoMarketPtr, dispersionMarketPtr

    ! Local variables
    integer :: direction
    !
    ! Actually, just open-up the cloud structures and passes the request further
    !
    if(fu_interval_positive(timestep))then
      direction = forwards
    else
      direction = backwards
    endif
    
    call set_initial_conditions(cloud%mapConc, &            ! mass map
                              & cloud%mapPx_conc, & !advection moment map
                              & cloud%mapPy_conc, & !advection moment map
                              & cloud%mapPz_conc, & !advection moment map
                              & cloud%mapDryDep, &          ! drydep map
                              & cloud%mapWetDep, &          ! wetdep map
                              & cloud%boundaryStructures(northern_boundary)%ptr, &  ! north
                              & cloud%boundaryStructures(southern_boundary)%ptr, &  ! south
                              & disp_buf, &
                              & rulesIniBoundary, &
                              & now, direction, &
                              & meteoMarketPtr, dispersionMarketPtr)

  end subroutine set_cloud_initial_conditions


  !****************************************************************

  subroutine prepare_transdep_cloud(cloud, meteo_ptr, disp_buf_ptr, chemRules)
    !
    ! Prepares the necessary pointers, meteodata and needed structures and variables
    ! to dynamics computations
    !
    implicit none

    TYPE(silam_pollution_cloud), INTENT(in) :: cloud
    TYPE(Tfield_buffer), POINTER :: meteo_ptr, disp_buf_ptr
    type(Tchem_rules), intent(inout) :: chemRules

    integer :: iTmp

    call msg_warning('Not implemented', 'prepare_transdep_cloud')

!!$    do iTmp = 1, size(cloud%cocktail_types)
!!$      if(cloud%cocktail_types(iTmp) == int_missing)exit
!!$      call prepare_cocktail_transdep(cloud%cocktail_types(iTmp), meteo_ptr, disp_buf_ptr, chemRules, &
!!$                                   & cloud%interpCoefMeteo2DispHoriz, cloud%ifMeteo2DispHorizInterp)
!!$    end do

  end subroutine prepare_transdep_cloud

  
  !****************************************************************
  
  subroutine check_cloud_species(cloud)
    !
    ! Existence of strange constructions, such as no_mode aerosol mode, requires check to
    ! make sure that this garbage is not in the cloud
    !
    implicit none

    ! Imported parameter
    TYPE(silam_pollution_cloud), INTENT(in) :: cloud
    
    ! Local variable
    integer :: iSp
    
    do iSp = 1, cloud%nSpTransport
      select case (fu_mode_type(cloud%speciesTransport(iSp)))
        case(gas_phase_flag, fixed_diameter_flag, gamma_function_flag, moving_diameter_flag, &
           & lognormal_flag, sea_salt_mode_flag, fire_mode_flag)
        case default
          call msg_warning('Strange mode of the below species', 'check_cloud_species')
          call report(cloud%speciesTransport(iSp))
          call set_error('Strange mode of the below species', 'check_cloud_species')
      end select
    end do
           
  end subroutine check_cloud_species
  

  ! ***************************************************************

  subroutine advect_pollution_cloud_v4(cloud,now,timestep, &
                                     & rw_method, &
                                     & have_negatives, &
                                     & IniBoundaryRules, &
                                     & met_buf, &
                                     & disp_buf, &
                                     & wdr, &
                                     & chem_rules, dynamics_rules)
    !
    ! Takes care of Eulerian advection. What has to be done:
    ! 1. Determine the meteo data interpolation coefficients along the
    !    vertical
    ! 2. Preparatory steps like timestep, etc. and then call for advection
    !
    implicit none

    ! Imported parameters with intent IN:
    TYPE(silja_time), INTENT(in) :: now
    TYPE(silja_interval), INTENT(in) :: timestep
    INTEGER, INTENT(in) :: rw_method
    TYPE(silja_wdr), INTENT(in) :: wdr
    logical, intent(in) :: have_negatives
    type(Tchem_rules), intent(in) :: chem_rules
    type(Tdynamics_rules), intent(in) :: dynamics_rules
    type(Tini_boundary_rules), intent(in) :: IniBoundaryRules

    ! Imported parameters with intent INOUT or POINTER
    TYPE(silam_pollution_cloud), INTENT(inout) :: cloud
    TYPE(Tfield_buffer), POINTER :: met_buf, disp_buf
    integer, dimension(:), pointer :: iArBadParticle
    REAL :: weight_past

!call msg_warning('Skipping all transport')
!return
    weight_past = met_buf%weight_past

    if(fu_if_lagrangian_present(dynamics_rules%simulation_type))then
      !
      ! Lagrangian advection-diffusion: only positions are needed. Masses are not affected
      !
      iArBadParticle => fu_work_int_array()
      if(error)return

      call advect_lagrangian_cloud(dynamics_rules%advMethod_Lagrangian, &
                                 & dynamics_rules%diffusionMethod, &
                                 & cloud%lpSet%lpDyn, cloud%lpSet%lpMassTrn, &
                                 & cloud%lpSet%nop, cloud%lpSet%nSpeciesTrn, &
                                 & cloud%lpSet%lpStatus, iArBadParticle, &
                                 & met_buf, & 
                                 & fu_sec(timestep), &
                                 & wdr)
      if(error)return
      !
      ! Handle the out-of-grid transport for limited-area and global grids
      !
      call collect_mass_budget_lagrangian(cloud, iArBadParticle)
      if(error)return

      call free_work_array(iArBadParticle)
    endif

    if(fu_if_eulerian_present(dynamics_rules%simulation_type))then
      !
      ! The main Eulerian advection
      !
      if(fu_ifBoundary(IniBoundaryRules))then
        !! Should be valid for mid-timestep
        call boundary_conditions_now(cloud%boundaryStructures, cloud%pBoundaryBuffer, now + timestep*0.5)
        if(error)return
      endif
#ifdef DEBUG
if(have_negatives)then
  call msg_warning('No mass check if negatives allowed: right before advection','advect_pollutiuon_cloud_v5')
else
  call check_masses(cloud,'right before advection', fu_low_mass_threshold(chem_rules))
  call check_mass_centres(cloud, 'right before advection', 1e-5)
endif
#endif

      call advect_eulerian_cloud(dynamics_rules%advMethod_Eulerian, &
                               & cloud%mapConc, &
                               & cloud%mapPx_conc, cloud%mapPy_conc, cloud%mapPz_conc, &
                               & cloud%mapAerosol, &
                               & cloud%interpCoefMeteo2DispHoriz, cloud%interpCoefMeteo2DispVert, &
                               & cloud%ifMeteo2DispHorizInterp, cloud%ifMeteo2DispVertInterp, &
                               & IniBoundaryRules, cloud%pBoundaryBuffer, &
                               & met_buf, &
                               & disp_buf, &
                               & now, fu_sec(timestep), weight_past, &
                               & cloud%garbage_mass, cloud%mapDryDep, cloud%mapConc2m,&
                               & cloud%x0_mass, cloud%xM_mass, cloud%y0_mass, cloud%yM_mass, &
                               & cloud%vert_mass_0, cloud%vert_mass_M, &
                               & chem_rules, have_negatives, wdr)
      if(error) return
      
#ifdef DEBUG
if(have_negatives)then
  call msg_warning('No mass check if negatives allowed: right after advection','advect_pollutiuon_cloud_v5')
else
  call check_masses(cloud,'right after advection', fu_low_mass_threshold(chem_rules))
  call check_mass_centres(cloud, 'right after  advection', 0.)
endif
#endif
      !
      ! Handle the out-of-grid transport, global closures, etc
      !
!      call collect_mass_budget_eulerian(cloud)
    !  if(error) call unset_error('advect_pollution_cloud_v4')

    endif  ! Simulation types

  end subroutine advect_pollution_cloud_v4


  ! ***************************************************************

  subroutine transform_pollution_cloud_v5(cld, now, timestep, &
                                        & tla_step, &
                                        & met_buf, disp_buf, &
                                        & meteo_input, &
                                        & wdr, &
                                        & chemRules, dynRules, &
                                        & ifDryDep_cumulative_in_output, ifWetDep_cumulative_in_output)
    !
    ! Takes care of Eulerian advection. What has to be done:
    ! 1. Determine the meteo data interpolation coefficients along the
    !    vertical
    ! 2. Preparatory steps like timestep, etc. and then call for advection
    !
    implicit none

    ! Imported parameters with intent IN:
    TYPE(silja_time), INTENT(in) :: now
    TYPE(silja_interval), INTENT(in) :: timestep
    type(t_tla_step), intent(inout) :: tla_step
    type(Tmeteo_input), pointer :: meteo_input
    TYPE(silja_wdr), INTENT(in) :: wdr
    type(Tchem_rules), intent(inout) :: chemRules
    type(Tdynamics_rules), intent(in) :: dynRules
    type(silja_logical), intent(in) :: ifDryDep_cumulative_in_output, ifWetDep_cumulative_in_output

    ! Imported parameters with intent INOUT or POINTER
    TYPE(silam_pollution_cloud), pointer :: cld
    TYPE(Tfield_buffer), POINTER :: met_buf, disp_buf

!call msg_warning('Skip all transformation')
!return

    !
    ! Transformations themselves are the same for both Lagrangian and Eulerian clouds
    !
    if (debug_level > 0) then
      if (defined(cld%mapConc)) &
       & call msg('Grand total of concentration map before transform_maps:',sum(cld%mapConc%arM(:,:,:,:,:)))
    end if

    !----------------------------------------------------------------------------
    !
    ! The main transformation is either Lagrangian or Eulerian, selector is encapsulated because
    ! we either call it here - in forward mode - or after deposition in adjoint mode.
    !
    if(fu_interval_positive(timestep)) call do_transformation()            ! the main transformation

    if (debug_level > 0 .and. fu_interval_positive(timestep)) then
      if (defined(cld%mapConc)) &
       & call msg('Grand total of concentration map before make_wet_deposition_map:', &
             & sum(cld%mapConc%arM(:,:,:,:,:)))
      call check_masses(cld,'after chemistry', fu_low_mass_threshold(chemRules))
    end if

    !
    ! Different cloud types have different deposition procedure and other specifics
    !
    if(dynRules%simulation_type == lagrangian_flag .or. &
     & dynRules%simulation_type == hybrid_flag)then
      !
      ! Both dry and wet deposition have to be handled explicitly
      !
      call start_count("deposition_lagr")
      call make_deposition_lagr(cld%lpSet%lpMassTrn, cld%lpSet%lpDyn, cld%lpSet%lpStatus, &
                              & cld%lpSet%nop, &
                              & cld%mapDryDep, cld%mapWetDep, cld%garbage_mass, &
                              & fu_sec(timestep), &  ! seconds !!
                              & met_buf%weight_past, &   ! meteo time interpolation
                              & chemRules%rulesDeposition)
      call stop_count("deposition_lagr")

    elseif(dynRules%simulation_type == eulerian_flag .or. &
         & dynRules%simulation_type == hybrid_flag)then

      ! Wet deposition has been moved to the main chemistry loop,
      ! so make_wet_deposition_map is no longer called!!!
      !
      ! Call the transformation of the cloud%mapConc.
      ! Since cloud%mapXCoord, cloud%mapYCoord, and cloud%mapZCoord contain the coordinates of
      ! centres of masses, not momentums, chemistry and radioactive decay will not affect them
      !
      !call start_count("wet_deposition_map")
      !call make_wet_deposition_map(cld%mapConc, cld%mapWetDep, cld%garbage_mass, &
      !                           & met_buf, disp_buf, tla_step, &
      !                           & cld%interpCoefMeteo2DispHoriz, cld%interpCoefMeteo2DispVert, &
      !                           & cld%ifMeteo2DispHorizInterp, cld%ifMeteo2DispVertInterp, &
      !                           & met_buf%weight_past, chemRules, fu_sec(timestep), now)
      !call stop_count("wet_deposition_map")
      !if(error)return

!#ifdef DEBUG
!     if (fu_interval_positive(timestep) .and. defined(cld%mapConc)) then
!        call msg('Grand total of concentration map after make_wet_deposition_map:', & 
!               & sum(cld%mapConc%arM(:,:,:,:,:)))
!        call check_masses(cld,'after wet dep', fu_low_mass_threshold(chemRules))
!      end if
!#endif
      !
      ! If the optical density structure is prepared, it has to be filled-in
      !
      if(cld%ifMakeOpticalDensity)then
        call start_count("optical_density")
        call make_optical_dens_map(cld%mapConc, &
                                 & met_buf, &
                                 & cld%interpCoefMeteo2DispHoriz, cld%interpCoefMeteo2DispVert, &
                                 & cld%ifMeteo2DispHorizInterp, cld%ifMeteo2DispVertInterp, &
                                 & cld%mapOptDns, &
                                 & fu_optical_rules(chemRules), &
                                 & fu_chemRunSetup(chemRules))
        call stop_count("optical_density")
        if(error)return
      endif
      if(cld%ifMakeOpticalColumnDepth)then
        !
        ! If optical density is available, column is a sum of density scaled with lyr thickness
        ! Otherwise have to compute everything. Trigger is the defined status of optioncal optical 
        ! density map in the below call
        !
        call start_count("optical_column_depth")
        call make_optical_column_depth_map(cld%mapConc, &
                                           & met_buf, &
                                           & cld%interpCoefMeteo2DispHoriz, &
                                           & cld%interpCoefMeteo2DispVert, &
                                           & cld%ifMeteo2DispHorizInterp, &
                                           & cld%ifMeteo2DispVertInterp, &
                                           & cld%mapOptColDepth, &
                                           & fu_optical_rules(chemRules), &
                                           & fu_chemRunSetup(chemRules), &
                                           & cld%mapOptDns) 
        call stop_count("optical_column_depth")
        if(error)return
      endif

      !call msg('Grand total of concentration map after transformation:',sum(cld%mapConc%arM(:,:,:,:,:)))

!      !
!      ! Having the masses defined, calculate the cloud droplet number concentration
!      !
!      if(fu_quantity_in_quantities(cloud_cond_nucley_nbr_cnc_flag, disp_buf%buffer_quantities)) &
!                      & call make_cld_droplet_nbr_cnc(cld%mapConc, met_buf, disp_buf)
!      if(error)return
      
      
    else
      call set_error('Strange dynamics rules','transform_pollution_cloud_v5')
    end if
    
    !-----------------------------------------------------------
    !
    ! In adjont mode, transformation is called after deposition
    !
    if (.not. fu_interval_positive(timestep)) call do_transformation()   ! adjoint transformation

    CONTAINS

    !===============================================================================
    
    subroutine do_transformation()
      implicit none
      !
      ! Encapsulation for the transformation call. Just a shorter code: in adjoint mode it has
      ! to be called after deposition, whereas in forward mode it is before it
      !
      if(dynRules%simulation_type == lagrangian_flag .or. &
       & dynRules%simulation_type == hybrid_flag)then
        !
        ! Lagrangian transformation: only positions are needed. Masses are not affected
        !
        call transform_lagrangian_part(cld%lpSet, &
                                     & cld%mapDryDep, cld%mapWetDep, &
                                     & cld%garbage_mass, &
                                     & tla_step, &
                                     & met_buf, disp_buf, &
                                     & meteo_input, &
                                     & cld%interpCoefMeteo2DispHoriz, cld%interpCoefMeteo2DispVert, &
                                     & cld%ifMeteo2DispHorizInterp, cld%ifMeteo2DispVertInterp, &
                                     & chemRules, fu_sec(timestep), now)
        if(error)return

      elseif(dynRules%simulation_type == eulerian_flag.or. &
           & dynRules%simulation_type == hybrid_flag)then
        !
        ! The main Eulerian transformation
        !
        call transform_maps(cld%mapConc, cld%mapShortlived, cld%mapAerosol, &
                          & cld%mapDryDep, cld%mapWetDep, &
                          & cld%mapReactRate, &
                          & cld%pBoundaryBuffer, &
                          & cld%garbage_mass, &
                          & tla_step, &
                          & met_buf, disp_buf, &
                          & meteo_input, &
                          & cld%interpCoefMeteo2DispHoriz, cld%interpCoefMeteo2DispVert, &
                          & cld%ifMeteo2DispHorizInterp, cld%ifMeteo2DispVertInterp, &
                          & ifDryDep_cumulative_in_output, ifWetDep_cumulative_in_output, &
                          & chemRules, fu_sec(timestep), now)
        if(error)return
      else
        call set_error('Strange simulation type:' + fu_str(dynRules%simulation_type), &
                     & 'transform_pollution_cloud_v4')
        return
      endif  ! Simulation types

    end subroutine do_transformation


  end subroutine transform_pollution_cloud_v5


  !*****************************************************************

  subroutine collect_mass_budget_eulerian(cloud)
    !
    ! Actually, just takes care of the masses transported outside the grid
    ! by Eulerian advection scheme. 
    ! Lagrangian routine takes care of them by itself
    ! Actual garbage - the lower-than-threshold-mass - is collected automatically
    ! by each routine, so here we do not have to do anything.
    ! The only non-trivial point is to collect the budget for global grid. Here the 
    ! round-the-globe closure for longitude has to be maintained.
    ! Now also latitude closure is available, see ini_boundary_conditions.

    implicit none

    ! Imported parameter
    TYPE(silam_pollution_cloud), intent(inout), target :: cloud

    ! Local variables
    integer :: iTmp, jTmp, iSrc, iz, iSubst, iMode, iSpecies
!    integer, dimension(:), pointer :: nModes
    TYPE(silam_pollution_cloud), pointer :: pc
    logical :: ifLonGlobal, ifNorthPole, ifSouthPole
    logical :: ifVertOutFluxRep
    real, dimension(:), pointer :: fOutN, fOutS, fOutE, fOutW  
   
    ifVertOutFluxRep = .false.
   
    if(ifVertOutFluxRep)then
      fOutN => fu_work_array()
      fOutS => fu_work_array()
      fOutW => fu_work_array()
      fOutE => fu_work_array()
    endif
    
    !
    ! Out-of-grid transport is made via one extra cell around teh computation domain.
    ! Our task is to collect it into 4 variables
    !
    ! Start from x == 0 and x == nx+1
    !
    pc => cloud

!    nSubst = fu_nbr_of_subst(cloud%cocktailTransp)  ! cocktail, ifFull
!    nModes => fu_work_int_array()
!    call n_modes_total(cloud%cocktailTransp, nModes)

    ifNorthPole = (cloud%boundaryStructures(northern_boundary)%ptr%boundaryType == polar_boundary_type)
    ifSouthPole = (cloud%boundaryStructures(southern_boundary)%ptr%boundaryType == polar_boundary_type)
    if(ifNorthPole) call msg('North pole global closure is active')
    if(ifSouthPole) call msg('South pole global closure is active')

    ifLonGlobal = fu_ifLonGlobal(dispersion_grid)

    if(ifLonGlobal) call msg('Longitude global closure is active')

    if(ifVertOutFluxRep)then
      call msg('========= OUT_FLUX_REPORT ========')
    endif
    
    do iSrc = 1, cloud%mapConc%nSrc
      do iz = 1, nz_dispersion
         if(ifVertOutFluxRep)then
           call msg('Flux_report_for_level:', iz)
           call msg('Species mode N S W E')
           fOutN(1:cloud%nSpTransport) = 0.0
           fOutS(1:cloud%nSpTransport) = 0.0
           fOutW(1:cloud%nSpTransport) = 0.0
           fOutE(1:cloud%nSpTransport) = 0.0
         endif
        !
        ! run along the x-borders: garbage or round-the-globe closure.
        ! Strictly speaking, centre of mass must also be taken into account but so far
        ! this seems to be not too bad either.
        !
!        if(ifLonGlobal)then     global clopsure is handled in advection
!          !
!          ! longitude closure
!          !
!          do iTmp = 1, ny_dispersion
!            cloud%mapConc%arM(1:cloud%nSpTransport,iSrc,iz,1,iTmp) = &
!                                          & cloud%mapConc%arM(:,iSrc,iz,1,iTmp) + &
!                                          & cloud%mapConc%arM(:,iSrc,iz,nx_dispersion+1,iTmp)
!            cloud%mapConc%arM(1:cloud%nSpTransport,iSrc,iz,nx_dispersion,iTmp) = &
!                                          & cloud%mapConc%arM(:,iSrc,iz,nx_dispersion,iTmp) + &
!                                          & cloud%mapConc%arM(:,iSrc,iz,0,iTmp)
!            cloud%mapConc%arM(1:cloud%nSpTransport,iSrc,iz,nx_dispersion+1,iTmp) = 0.0
!            cloud%mapConc%arM(1:cloud%nSpTransport,iSrc,iz,0,iTmp) = 0.0
!          end do
!        else
        if(.not. ifLonGlobal)then
          !
          ! garbage collection
          !
          do iTmp = 1, ny_dispersion
            cloud%x0_mass(iSrc,1:cloud%nSpTransport,outgoing,iz) = cloud%x0_mass(iSrc,:,outgoing,iz) + &
                                                    & cloud%mapConc%arM(:,iSrc,iz,0,iTmp)
            if(ifVertOutFluxRep)then
              fOutW(1:cloud%nSpTransport) = fOutW(:) + cloud%mapConc%arM(:,iSrc,iz,0,iTmp) 
            endif
            cloud%xM_mass(iSrc,1:cloud%nSpTransport,outgoing,iz) = cloud%xM_mass(iSrc,:,outgoing,iz) + &
                                                    & cloud%mapConc%arM(:,iSrc,iz,nx_dispersion+1,iTmp)
            if(ifVertOutFluxRep)then
              fOutE(1:cloud%nSpTransport) = fOutE(:) + cloud%mapConc%arM(:,iSrc,iz,nx_dispersion+1,iTmp)
            endif
            cloud%mapConc%arM(1:cloud%nSpTransport,iSrc,iz,nx_dispersion+1,iTmp) = 0.0
            cloud%mapConc%arM(1:cloud%nSpTransport,iSrc,iz,0,iTmp) = 0.0
          end do
        endif  ! fu_ifLonGlobal

!        if (ifNorthPole) then          poles are handled in advection: they now deal with moments etc.
!          !
!          ! latitude: move mass into the polar cell
!          !
!          do iTmp = 1, nx_dispersion
!            cloud%boundaryStructures(northern_boundary)%ptr%polarCapMassMom(1:cloud%nSpTransport,iSrc,1,iz) = &
!                       & cloud%boundaryStructures(northern_boundary)%ptr%polarCapMassMom(:, iSrc, 1, iz) + &
!                       & cloud%mapConc%arM(:,iSrc,iz,iTmp,ny_dispersion+1) !/ &
!!                       & (cloud%boundaryStructures(northern_boundary)%ptr%PolarCapArea * &
!!                        & cloud%boundaryStructures(northern_boundary)%ptr%PolarCapThickness(iz))
!            cloud%mapConc%arM(1:cloud%nSpTransport,iSrc,iz,iTmp,ny_dispersion+1) = 0.0
!          end do
!        else
        if (.not. ifNorthPole) then
          !
          ! or garbage collection: 
          ! run along the y-borders
          !
          do iTmp = 0, nx_dispersion+1
            cloud%yM_mass(iSrc,1:cloud%nSpTransport,outgoing,iz) = cloud%yM_mass(iSrc,:,outgoing,iz) + &
                                                    & cloud%mapConc%arM(:,iSrc,iz,iTmp,ny_dispersion+1)
            if(ifVertOutFluxRep)then                                    
              fOutN(1:cloud%nSpTransport) = fOutN(:) + cloud%mapConc%arM(:,iSrc,iz,iTmp,ny_dispersion+1)                                    
            endif                                    
            cloud%mapConc%arM(1:cloud%nSpTransport,iSrc,iz,iTmp,ny_dispersion+1) = 0.0
          end do
          
        end if
!        if (ifSouthPole) then           poles are handled in advection: they now deal with moments etc.
!          !
!          ! latitude: move mass into the polar cell
!          !
!          do iTmp = 1, nx_dispersion
!            cloud%boundaryStructures(southern_boundary)%ptr%polarCapMassMom(1:cloud%nSpTransport,iSrc,1,iz) =&
!                       & cloud%boundaryStructures(southern_boundary)%ptr%polarCapMassMom(:, iSrc, 1, iz) + &
!                       & cloud%mapConc%arM(:,iSrc,iz,iTmp,0) !/ &
!!                       & (cloud%boundaryStructures(southern_boundary)%ptr%PolarCapArea * &
!!                        & cloud%boundaryStructures(southern_boundary)%ptr%PolarCapThickness(iz))
!            cloud%mapConc%arM(1:cloud%nSpTransport,iSrc,iz,iTmp,0) = 0.0
!          end do
!        else
        if (.not. ifSouthPole) then
          !
          ! or garbage collection: 
          ! run along the y-borders
          !
          do iTmp = 0, nx_dispersion+1
            cloud%y0_mass(iSrc,1:cloud%nSpTransport,outgoing,iz) = cloud%y0_mass(iSrc,:,outgoing,iz) + &
                                          & cloud%mapConc%arM(:,iSrc,iz,iTmp,0)
            if(ifVertOutFluxRep)then                              
              fOutS(1:cloud%nSpTransport) = fOutS(1:cloud%nSpTransport) + &
                                       & cloud%mapConc%arM(1:cloud%nSpTransport,iSrc,iz,iTmp,0)                              
            endif                              
            cloud%mapConc%arM(1:cloud%nSpTransport,iSrc,iz,iTmp,0) = 0.0
          end do
        end if
 
        if(ifVertOutFluxRep)then
          do itmp = 1, cloud%nSpTransport
            if (smpi_global_rank == 0) then
              print('(A15,4(E10.5,1x))'), fu_str(cloud%speciesTransport(itmp)), &
                   & fOutN(itmp), fOutS(itmp), fOutW(itmp), fOutE(itmp)
            end if
            write(run_log_funit, '(A15,4(E10.5,1x))') &
                                     & fu_str(cloud%speciesTransport(itmp)), &
                                     & fOutN(itmp), fOutS(itmp), fOutW(itmp), fOutE(itmp)
          enddo
        endif
      end do  ! iz

      !
      ! Finally, the over-the-top and underground planes
      !
      do jTmp = 0, ny_dispersion+1
        do iTmp = 0, nx_dispersion+1
          cloud%vert_mass_M(iSrc,1:cloud%nSpTransport,outgoing) = cloud%vert_mass_M(iSrc,:,outgoing) + &
                                               & cloud%mapConc%arM(:,iSrc, nz_dispersion+1,iTmp,jTmp)
          cloud%vert_mass_0(iSrc,1:cloud%nSpTransport, outgoing) = cloud%vert_mass_0(iSrc,:,outgoing) + &
                                               & cloud%mapConc%arM(:,iSrc, 0,iTmp,jTmp)
          cloud%mapConc%arM(1:cloud%nSpTransport,iSrc, nz_dispersion+1,iTmp,jTmp) = 0.0
          cloud%mapConc%arM(1:cloud%nSpTransport,iSrc, 0,iTmp,jTmp) = 0.0
        end do
      end do

    end do  ! iSrc
    if(ifVertOutFluxRep)then
      call free_work_array(fOutN)
      call free_work_array(fOutS)
      call free_work_array(fOutW)
      call free_work_array(fOutE)
      call msg('========= END_OUT_FLUX_REPORT ========')
     endif
!    call free_work_array(nModes)

  end subroutine collect_mass_budget_eulerian


  !*************************************************************************

  subroutine collect_mass_budget_lagrangian(cloud, iArBadParticle)
    !
    ! Actually, just takes care of the masses transported outside the grid
    ! by Lagrangian advection scheme. 
    ! Actual garbage - the lower-than-threshold-mass - is assumed for every broken particle
    ! that is still inside the grid.
    ! The non-trivial point is to collect the budget for global grid. Here the 
    ! round-the-globe closure for longitude has to be maintained.
    ! Also the latitude closure is available, see ini_boundary_conditions.

    implicit none

    ! Imported parameter
    TYPE(silam_pollution_cloud), intent(inout), target :: cloud
    integer, dimension(:), pointer :: iArBadParticle

    ! Local variables
    integer :: iPart, iSrc, iy, iz, iSpecies, iTmp
    real, dimension(:,:), pointer :: lpDynLocal, lpMassLocal
    integer, dimension(:), pointer :: lpStatusLocal

    !
    ! Out-of-grid transport is made via one extra cell around teh computation domain.
    ! Our task is to collect it into 4 variables
    !
    ! Start from x == 0 and x == nx+1
    !
    
    lpDynLocal => cloud%lpSet%lpDyn
    lpMassLocal => cloud%lpSet%lpMassTrn
    lpStatusLocal => cloud%lpSet%lpStatus
    !
    ! Loop over bad particles (marked so in advection)
    !
    do iTmp = 1, size(iArBadParticle)

      if(iArBadParticle(iTmp) == int_missing)exit
      iPart = iArBadParticle(iTmp)
    
      iSrc = mod(lpStatusLocal(iPart),100)

      !
      ! In what direction did it go outside?
      !
      if(lpDynLocal(lp_x,iPart) < 0)then
        cloud%x0_mass(iSrc,1:cloud%nSpTransport,outgoing,1) = cloud%x0_mass(iSrc,:,outgoing,1) + &
                                                 & lpMassLocal(1:cloud%nSpTransport,iPart)
      elseif(lpDynLocal(lp_x,iPart) > nx_dispersion + 0.5)then
        cloud%xM_mass(iSrc,1:cloud%nSpTransport,outgoing,1) = cloud%xM_mass(iSrc,:,outgoing,1) + &
                                                 & lpMassLocal(1:cloud%nSpTransport,iPart)
      elseif(lpDynLocal(lp_y,iPart) < 0)then
          cloud%y0_mass(iSrc,1:cloud%nSpTransport,outgoing,1) = cloud%y0_mass(iSrc,1:cloud%nSpTransport,outgoing,1) + &
                                                   & lpMassLocal(1:cloud%nSpTransport,iPart)
      elseif(lpDynLocal(lp_y,iPart) > ny_dispersion + 0.5)then
          cloud%yM_mass(iSrc,1:cloud%nSpTransport,outgoing,1) = cloud%yM_mass(iSrc,1:cloud%nSpTransport,outgoing,1) + &
                                                   & lpMassLocal(1:cloud%nSpTransport,iPart)
      elseif(lpDynLocal(lp_z,iPart) < 0)then
        cloud%vert_mass_0(iSrc,1:cloud%nSpTransport,outgoing) = cloud%vert_mass_0(iSrc,1:cloud%nSpTransport,outgoing) + &
                                                     & lpMassLocal(1:cloud%nSpTransport,iPart)
      elseif(lpDynLocal(lp_z,iPart) > nz_dispersion + 0.5)then
        cloud%vert_mass_M(iSrc,1:cloud%nSpTransport,outgoing) = cloud%vert_mass_M(iSrc,1:cloud%nSpTransport,outgoing) + &
                                                     & lpMassLocal(1:cloud%nSpTransport,iPart)
      else
        cloud%garbage_mass(iSrc,1:cloud%nSpTransport) = cloud%garbage_mass(iSrc,1:cloud%nSpTransport) + &
                                                      & lpMassLocal(1:cloud%nSpTransport,iPart)
      endif

      lpMassLocal(1:cloud%nSpTransport,iPart) = 0.0
      lpStatusLocal(iPart) = int_missing
      if(cloud%lpSet%iFirstEmptyParticle > iPart) cloud%lpSet%iFirstEmptyParticle = iPart
    end do  ! iPart

  end subroutine collect_mass_budget_lagrangian



  !*************************************************************************
  !
  !  Reference functions for the pollution cloud class
  !
  !*************************************************************************


  !******************************************************************

  function fu_if_cloud_mass_map_quantity(iQ) result(ifInCloud)
    !
    ! Returns true if the quantity is attributed to cloud mass map and 
    ! false if it is a field-related parameter. Since an error can cost a lot,
    ! a strict checking is used.
    !
    implicit none

    type(silja_logical) :: ifInCloud

    integer, intent(in) :: iQ ! quantity itself

    select case(iQ)
      case(particle_counter_flag, &
         & emission_intensity_flag, emission_flux_flag, &
         & concentration_flag, &
         & concentration_2m_flag, &
         & areas_of_risk_flag, &
         & mass_in_air_flag, &
         & advection_moment_X_flag, advection_moment_Y_flag, advection_moment_Z_flag, &
         & optical_density_flag, &
         & optical_column_depth_flag, &
         & drydep_flag, &
         & wetdep_flag, &
         & volume_mixing_ratio_flag)

        ifInCloud = silja_true

      case(heatsum_flag, chillsum_flag, start_heatsum_threshold_flag, daily_temp_threshold_flag, &
         & soil_moisture_threshold_flag, &
         & start_calday_threshold_flag, end_calday_threshold_flag, &
         & pollen_rdy_to_fly_flag, allergen_rdy_to_fly_flag, &
         & calday_start_end_diff_flag, heatsum_start_end_diff_flag, temperature_threshold_flag, &
         & growth_season_start_day_flag, heatsum_cutoff_tempr_flag, end_heatsum_threshold_flag, &
         & day_mean_temperature_flag, day_mean_temperature_2m_flag, &
         & day_temperature_acc_flag, day_temperature_2m_acc_flag, &
         & pollen_left_relative_flag, pollen_total_per_m2_flag, pollen_correction_flag, &
         & pollen_potency_flag, plant_growth_flag,  &
         & physiography_field_set_flag, &
         & Vd_correction_DMAT_flag, &
         & interp_met2disp_coef_flag, &
         & interp_met2out_coef_flag, &
         & interp_disp2out_coef_flag,&
         & cell_size_z_flag, &
         & cell_size_x_flag, &
         & cell_size_y_flag, &
         & dispersion_u_flag, &
         & dispersion_v_flag, &
         & dispersion_w_flag, &
         & disp_cell_airmass_flag, &
         & disp_flux_celleast_flag, disp_flux_cellnorth_flag,  disp_flux_celltop_flag, &
         & disp_flux_cellt_rt_flag,disp_flux_celle_rt_flag,disp_flux_celln_rt_flag, &
         & large_scale_rain_int_flag, &
         & convective_rain_int_flag, &
         & total_precipitation_rate_flag, &
         & air_density_flag)

        ifInCloud = silja_false

      case default
        call msg_warning('Unknown quantity:' + fu_quantity_string(iQ), 'fu_if_cloud_mass_map_quantity')
        ifInCloud = silja_undefined

    end select

  end function fu_if_cloud_mass_map_quantity


  !*****************************************************************

  logical function fu_if_variable_available(Cloud, pDispersionMarket, quantity, species, nlStdSetup, &
                                          & dynRules, ifVerbose)
    !
    ! Checks if the quantity can be computed from the cloud data or obtained from dispersion stack
    !
    implicit none

    ! Imported parameters
    type(silam_pollution_cloud), intent(in) :: cloud
    type(mini_market_of_stacks), pointer :: pDispersionMarket
    integer, intent(in) :: quantity
    type(silam_species), intent(in) :: species
    type(Tsilam_namelist), pointer :: nlStdSetup
    type(Tdynamics_rules), intent(in) :: dynRules
    logical, intent(in) :: ifVerbose
    
    ! Local variables
    type(silja_field_id) :: id
    type(silja_field), pointer :: pField
    integer :: iSpecies
    !
    ! First check the cloud.
    ! There are not too many quantities, which can be derived from the pollution cloud data
    !
    select case(quantity)
      case(emission_intensity_flag, emission_flux_flag)          ! Only if species is emitted
        fu_if_variable_available = fu_if_eulerian_present(dynRules%simulation_type) .and. &
                                 & (fu_index(species, Cloud%speciesEmission, Cloud%nSpEmission) > 0)

      case(concentration_flag, concentration_2m_flag, mass_in_air_flag, volume_mixing_ratio_flag)  ! Only if species is transported
        fu_if_variable_available = (fu_index(species, Cloud%speciesTransport, Cloud%nSpTransport) > 0)

      case(optical_density_flag, optical_column_depth_flag)  ! Only if optical params are computed
        !
        ! Careful: optical species might not yet be defined. In this case will check the transport ones for
        ! general availability and then the wave length and optic-specific properties separately
        !
        if(Cloud%nSpOpticalDns > 0)then
          !
          ! Optic species are defined. Easy
          !
          fu_if_variable_available = (fu_index(species, Cloud%speciesOpticalDensity, Cloud%nSpOpticalDns) > 0)
        else
          !
          ! Common case: optics are not yet initialised
          !
          do iSpecies = 1, Cloud%nSpTransport
            if(associated(species%material, Cloud%speciesTransport(iSpecies)%material))then
              if(fu_mode(species) == fu_mode(Cloud%speciesTransport(iSpecies)))then
                ! reference to optic defaults
                fu_if_variable_available = fu_if_species_optics_known(species, real_missing) 
                return
              endif
            endif
          end do
          fu_if_variable_available = .false.
        endif

      case(particle_counter_flag, areas_of_risk_flag)        ! Lagrangian advection only
        fu_if_variable_available = fu_if_lagrangian_present(dynRules%simulation_type)

      case(advection_moment_X_flag, advection_moment_Y_flag, advection_moment_Z_flag)  ! Euler adv only
        fu_if_variable_available = fu_if_eulerian_present(dynRules%simulation_type)

      case(drydep_flag, wetdep_flag)                         ! For aerosol or non-passive gas
        if(fu_index(species, Cloud%speciesTransport, Cloud%nSpTransport) > 0)then
          fu_if_variable_available = (.not. fu_if_gas(fu_material(species))) .or. &  ! not gas
                                   & fu_if_gas_depositing(fu_material(species))      ! depositing gas
        else
          fu_if_variable_available = .false.
        endif

      case default
        !
        ! The quantity is not stored in pollution cloud. Check the dispersion market
        !
        id = fu_set_field_id_simple(met_src_missing, quantity, time_missing, level_missing, species)
        
        pField => fu_get_field_from_mm_general(pDispersionMarket, id, .false.)
        
        if(associated(pField))then
          fu_if_variable_available = defined(pField)
        else
          fu_if_variable_available = .false.
          if(ifVerbose)then
            call msg('Failed to find the field in the market:')
            call report(id)
            call msg('Fields available in the market:')
            call report(pDispersionMarket)
          endif
        endif

    end select

  end function fu_if_variable_available


  !**********************************************************************************************

  logical function fu_transport_owned_quantity(pCld, quantity)
    !
    ! Check whether the specific quantity is used and handled in the emission part of the process.
    !
    implicit none
    
    ! Imported parameters
    type(silam_pollution_cloud), pointer :: pCld
    integer, intent(in) :: quantity
    
    select case(quantity)
      case(particle_counter_flag, areas_of_risk_flag, &
         & concentration_flag, mass_in_air_flag, concentration_2m_flag, volume_mixing_ratio_flag, &
         & drydep_flag, wetdep_flag, &
         & advection_moment_X_flag, advection_moment_Y_flag, advection_moment_Z_flag, &
         & aerosol_flag)
        fu_transport_owned_quantity = .true.
      case default
        fu_transport_owned_quantity = .false.
    end select
    
  end function fu_transport_owned_quantity

  !**********************************************************************************************

  logical function fu_optics_owned_quantity(pCld, quantity)
    !
    ! Check whether the specific quantity is used and handled in the emission part of the process.
    !
    implicit none
    
    ! Imported parameters
    type(silam_pollution_cloud), pointer :: pCld
    integer, intent(in) :: quantity
    
    select case(quantity)
      case(optical_density_flag, optical_column_depth_flag)
        fu_optics_owned_quantity = .true.
      case default
        fu_optics_owned_quantity = .false.
    end select
    
  end function fu_optics_owned_quantity


  !******************************************************************
  !******************************************************************
  !
  !  Some encapsulation
  !
  !******************************************************************
  !******************************************************************

  LOGICAL FUNCTION fu_cloud_defined(cloud)
    !
    ! Returns a true value, if the cloud has been given a value using
    ! one of the setting functions.
    !
    ! Author: Mika Salonoja, FMI
    ! 
    IMPLICIT NONE

    ! Imported parameters with intent IN:
    TYPE(silam_pollution_cloud), INTENT(in) :: cloud

    fu_cloud_defined = fu_true(cloud%defined)

  END FUNCTION fu_cloud_defined


  ! ***************************************************************


  function fu_species_emission_cld(cloud) result(species)
    !
    ! Returns the pointer to the cocktail of the requested particle
    !
    implicit none

    !Return value of the function
    type(silam_species), dimension(:), pointer :: species

    ! Imported parameters with intent IN:
    TYPE(silam_pollution_cloud), intent(in) :: cloud

    species => cloud%speciesEmission

  end function fu_species_emission_cld


  ! ***************************************************************


  integer function fu_nbr_of_species_emission(cloud)
    !
    ! Returns the pointer to the cocktail of the requested particle
    !
    implicit none

    ! Imported parameters with intent IN:
    TYPE(silam_pollution_cloud), intent(in) :: cloud

    fu_nbr_of_species_emission = cloud%nSpEmission

  end function fu_nbr_of_species_emission


  ! ***************************************************************


  function fu_species_transport(cloud) result(species)
    !
    ! Returns the pointer to the cocktail of the requested particle
    !
    implicit none

    !Return value of the function
    type(silam_species), dimension(:), pointer :: species

    ! Imported parameters with intent IN:
    TYPE(silam_pollution_cloud), intent(in) :: cloud

    species => cloud%speciesTransport

  end function fu_species_transport


  ! ***************************************************************


  integer function fu_nbr_of_species_transport(cloud)
    !
    ! Returns the pointer to the cocktail of the requested particle
    !
    implicit none

    ! Imported parameters with intent IN:
    TYPE(silam_pollution_cloud), intent(in) :: cloud

    fu_nbr_of_species_transport = cloud%nSpTransport

  end function fu_nbr_of_species_transport


  ! ***************************************************************


  function fu_species_optical(cloud) result(species)
    !
    ! Returns the pointer to the cocktail of the requested particle
    !
    implicit none

    !Return value of the function
    type(silam_species), dimension(:), pointer :: species

    ! Imported parameters with intent IN:
    TYPE(silam_pollution_cloud), intent(in) :: cloud

    species => cloud%speciesOpticalDensity

  end function fu_species_optical


  ! ***************************************************************


  integer function fu_nbr_of_species_optical(cloud)
    !
    ! Returns the pointer to the cocktail of the requested particle
    !
    implicit none

    ! Imported parameters with intent IN:
    TYPE(silam_pollution_cloud),  intent(in) :: cloud

    fu_nbr_of_species_optical = cloud%nSpOpticalDns

  end function fu_nbr_of_species_optical


  ! ***************************************************************


  function fu_species_short_lived(cloud) result(species)
    !
    ! Returns the pointer to the cocktail of the requested particle
    !
    implicit none

    !Return value of the function
    type(silam_species), dimension(:), pointer :: species

    ! Imported parameters with intent IN:
    TYPE(silam_pollution_cloud), intent(in), target :: cloud

    species => cloud%speciesShortLived

  end function fu_species_short_lived


  ! ***************************************************************


  integer function fu_nbr_of_species_short_lived(cloud)
    !
    ! Returns the pointer to the cocktail of the requested particle
    !
    implicit none

    ! Imported parameters with intent IN:
    TYPE(silam_pollution_cloud), intent(in) :: cloud

    fu_nbr_of_species_short_lived = cloud%nSpShortLived

  end function fu_nbr_of_species_short_lived


  ! ***************************************************************


  function fu_species_aerosol(cloud) result(species)
    !
    ! Returns the pointer to the cocktail of the requested particle
    !
    implicit none

    !Return value of the function
    type(silam_species), dimension(:), pointer :: species

    ! Imported parameters with intent IN:
    TYPE(silam_pollution_cloud), intent(in), target :: cloud

    species => cloud%speciesAerosol

  end function fu_species_aerosol


  ! ***************************************************************


  integer function fu_nbr_of_species_aerosol(cloud)
    !
    ! Returns the pointer to the cocktail of the requested particle
    !
    implicit none

    ! Imported parameters with intent IN:
    TYPE(silam_pollution_cloud), intent(in) :: cloud

    fu_nbr_of_species_aerosol = cloud%nSpAerosol

  end function fu_nbr_of_species_aerosol


  ! ***************************************************************


  INTEGER FUNCTION fu_nbr_of_sources_of_cloud(cloud)
    !
    ! Returns the number of sources in the cloud
    ! 
    IMPLICIT NONE

    ! Imported parameters with intent IN:
    TYPE(silam_pollution_cloud), INTENT(in) :: cloud

    ! local variables
    integer :: i, iSrc

    if(defined(cloud%mapConc))then
      fu_nbr_of_sources_of_cloud = cloud%mapConc%nSrc
    else
      if(cloud%lpSet%nSrcs > 0)then
        fu_nbr_of_sources_of_cloud = cloud%lpSet%nSrcs
      else
        call set_error('Concentration map is not associated and lpSet is empty', &
                     & 'fu_nbr_of_sources_of_cloud')
        fu_nbr_of_sources_of_cloud = int_missing
      endif
    endif

  END FUNCTION fu_nbr_of_sources_of_cloud


  ! ***************************************************************


  FUNCTION fu_cloud_earliest_start(cloud) RESULT(start)
  !
  ! Returns the earliest start time of the cloud particles
  !
    IMPLICIT NONE

    ! Return value
    TYPE(silja_time) :: start
    ! Imported parameters with intent IN:
    TYPE(silam_pollution_cloud), INTENT(in), target :: cloud

    start = cloud%earliest_start

  END FUNCTION fu_cloud_earliest_start


  !******************************************************************

  function fu_emisMM_ptr(cloud)result(ptr)
    !
    ! Returns the pointer to the emission mass map
    !
    IMPLICIT NONE
    ! Return value 
    type(Tmass_map), pointer :: ptr
    ! Imported parameters 
    TYPE(silam_pollution_cloud), INTENT(in), target :: cloud

    ptr => cloud%mapEmis
  end function fu_emisMM_ptr
  

  !******************************************************************
  function fu_reactRateMM_ptr(cloud)result(ptr)
    !
    ! Returns the pointer to the reaction_rate massmap
    !
    IMPLICIT NONE
    ! Return value 
    type(Tmass_map), pointer :: ptr
    ! Imported parameters 
    TYPE(silam_pollution_cloud), INTENT(in), target :: cloud

    ptr => cloud%mapReactRate
  end function fu_reactRateMM_ptr
  

  !******************************************************************

  function fu_concMM_ptr(cloud)result(ptr)
    !
    ! Returns the pointer to the concentration mass map
    !
    IMPLICIT NONE
    ! Return value 
    type(Tmass_map), pointer :: ptr
    ! Imported parameters 
    TYPE(silam_pollution_cloud), INTENT(in), target :: cloud

    ptr => cloud%mapConc
  end function fu_concMM_ptr
  

  !******************************************************************

  function  fu_concentration2mMM_ptr(cloud)result(ptr)
    !
    ! Returns the pointer to the concentration mass map
    !
    IMPLICIT NONE
    ! Return value 
    type(Tmass_map), pointer :: ptr
    ! Imported parameters 
    TYPE(silam_pollution_cloud), INTENT(in), target :: cloud

    ptr => cloud%mapConc2m
  end function  fu_concentration2mMM_ptr
  

  !******************************************************************

  function fu_optical_densityMM_ptr(cloud)result(ptr)
    !
    ! Returns the pointer to the optical density mass map
    !
    IMPLICIT NONE
    ! Return value 
    type(Tmass_map), pointer :: ptr
    ! Imported parameters 
    TYPE(silam_pollution_cloud), INTENT(in), target :: cloud

    ptr => cloud%mapOptDns
  end function fu_optical_densityMM_ptr


  !******************************************************************

  function fu_optical_column_depthMM_ptr(cloud)result(ptr)
    !
    ! Returns the pointer to the optical density mass map
    !
    IMPLICIT NONE
    ! Return value 
    type(Tmass_map), pointer :: ptr
    ! Imported parameters 
    TYPE(silam_pollution_cloud), INTENT(in), target :: cloud

    ptr => cloud%mapOptColDepth
  end function fu_optical_column_depthMM_ptr


  !******************************************************************

  function fu_aerosolMM_ptr(cloud)result(ptr)
    !
    ! Returns the pointer to the concentration mass map
    !
    IMPLICIT NONE
    ! Return value 
    type(Tmass_map), pointer :: ptr
    ! Imported parameters 
    TYPE(silam_pollution_cloud), INTENT(in), target :: cloud

    ptr => cloud%mapAerosol
  end function fu_aerosolMM_ptr
  

  !******************************************************************

  function fu_advection_moment_X_MM_ptr(cloud)result(ptr)
    !
    ! Returns the pointer to the array of the dry deposition cocktails
    !
    IMPLICIT NONE
    ! Return value 
    type(Tmass_map), pointer :: ptr
    ! Imported parameters 
    TYPE(silam_pollution_cloud), INTENT(in), target :: cloud

    ptr => cloud%mapPx_conc
  end function fu_advection_moment_X_MM_ptr
  

  !******************************************************************

  function fu_advection_moment_Y_MM_ptr(cloud)result(ptr)
    !
    ! Returns the pointer to the array of the dry deposition cocktails
    !
    IMPLICIT NONE
    ! Return value 
    type(Tmass_map), pointer :: ptr
    ! Imported parameters 
    TYPE(silam_pollution_cloud), INTENT(in), target :: cloud

    ptr => cloud%mapPy_conc
  end function fu_advection_moment_Y_MM_ptr
  

  !******************************************************************

  function fu_advection_moment_Z_MM_ptr(cloud)result(ptr)
    !
    ! Returns the pointer to the array of the dry deposition cocktails
    !
    IMPLICIT NONE
    ! Return value 
    type(Tmass_map), pointer :: ptr
    ! Imported parameters 
    TYPE(silam_pollution_cloud), INTENT(in), target :: cloud

    ptr => cloud%mapPz_conc
  end function fu_advection_moment_Z_MM_ptr
  
  !******************************************************************

  function fu_emission_moment_X_MM_ptr(cloud)result(ptr)
    !
    ! Returns the pointer to the array of the dry deposition cocktails
    !
    IMPLICIT NONE
    ! Return value 
    type(Tmass_map), pointer :: ptr
    ! Imported parameters 
    TYPE(silam_pollution_cloud), INTENT(in), target :: cloud

    ptr => cloud%mapPx_emis
  end function fu_emission_moment_X_MM_ptr
  

  !******************************************************************

  function fu_emission_moment_Y_MM_ptr(cloud)result(ptr)
    !
    ! Returns the pointer to the array of the dry deposition cocktails
    !
    IMPLICIT NONE
    ! Return value 
    type(Tmass_map), pointer :: ptr
    ! Imported parameters 
    TYPE(silam_pollution_cloud), INTENT(in), target :: cloud

    ptr => cloud%mapPy_emis
  end function fu_emission_moment_Y_MM_ptr
  

  !******************************************************************

  function fu_emission_moment_Z_MM_ptr(cloud)result(ptr)
    !
    ! Returns the pointer to the array of the dry deposition cocktails
    !
    IMPLICIT NONE
    ! Return value 
    type(Tmass_map), pointer :: ptr
    ! Imported parameters 
    TYPE(silam_pollution_cloud), INTENT(in), target :: cloud

    ptr => cloud%mapPz_emis
  end function fu_emission_moment_Z_MM_ptr

  !******************************************************************

  function fu_drydepMM_ptr(cloud)result(ptr)
    !
    ! Returns the pointer to the array of the dry deposition cocktails
    !
    IMPLICIT NONE
    ! Return value 
    type(Tmass_map), pointer :: ptr
    ! Imported parameters 
    TYPE(silam_pollution_cloud), INTENT(in), target :: cloud

    ptr => cloud%mapDryDep
  end function fu_drydepMM_ptr
  
  !******************************************************************

  function fu_wetdepMM_ptr(cloud)result(ptr)
    !
    ! Returns the pointer to the array of the dry deposition cocktails
    !
    IMPLICIT NONE
    ! Return value 
    type(Tmass_map), pointer :: ptr
    ! Imported parameters 
    TYPE(silam_pollution_cloud), INTENT(in), target :: cloud

    ptr => cloud%mapWetDep
  end function fu_wetdepMM_ptr
  
  !************************************************************************************
  
  function fu_shortLivedMM_ptr(cloud) result(ptr)
    implicit none
    type(Tmass_map), pointer :: ptr
    type(silam_pollution_cloud), intent(in), target :: cloud
    ptr => cloud%mapShortlived
  end function fu_shortLivedMM_ptr

  !************************************************************************************
  
  function fu_concLP_ptr(cloud) result(ptr)
    implicit none
    real, dimension(:,:), pointer :: ptr
    type(silam_pollution_cloud), intent(in), target :: cloud
    ptr => cloud%lpSet%lpMassTrn
  end function fu_concLP_ptr

  !************************************************************************************
  
  function fu_aerosolLP_ptr(cloud) result(ptr)
    implicit none
    real, dimension(:,:), pointer :: ptr
    type(silam_pollution_cloud), intent(in), target :: cloud
    ptr => cloud%lpSet%lpMassAer
  end function fu_aerosolLP_ptr

  !************************************************************************************

  function fu_LP_dynamic_params(cloud) result(ptr)
    implicit none
    real, dimension(:,:), pointer :: ptr
    type(silam_pollution_cloud), intent(in), target :: cloud
    ptr => cloud%lpSet%lpDyn
  end function fu_LP_dynamic_params

  !************************************************************************************

  function fu_LP_status(cloud) result(ptr)
    implicit none
    integer, dimension(:), pointer :: ptr
    type(silam_pollution_cloud), intent(in), target :: cloud
    ptr => cloud%lpSet%lpStatus
  end function fu_LP_status

  !************************************************************************************
  
  function fu_lpSet_ptr(cloud)result(ptr)
    implicit none
    type(Tlagrange_particles_set), pointer :: ptr
    type(silam_pollution_cloud), intent(in), target :: cloud
    ptr => cloud%lpSet
  end function fu_lpSet_ptr

  !************************************************************************************
  
  integer function fu_nbr_of_lagr_particles(cloud, ifActive)
    implicit none
    type(silam_pollution_cloud), intent(in) :: cloud
    logical, intent(in) :: ifActive
    if(ifActive)then
      fu_nbr_of_lagr_particles = cloud%lpSet%nop - 1
    else
      fu_nbr_of_lagr_particles = size(cloud%lpSet%lpStatus)
    endif
  end function fu_nbr_of_lagr_particles

  !*****************************************************************

  subroutine get_cloud_inventory(cloud, ifActive, iSpeciesSelector, pSpecies, nSpecies, moles)
    !
    ! Returns the current-time inventory of the cloud.
    ! Asnwer depends on whether transport or emission inventory requested. Note that 
    ! the list of substances can be significantly different!
    !
    ! ATTENTION. Quite slow because scans all particles and merges their 
    !            inventories
    !
    implicit none

    ! Imported parameters with intent IN
    TYPE(silam_pollution_cloud), INTENT(in), target :: cloud
    integer, intent(in) :: iSpeciesSelector  ! Type of the species: emission, transport, ...
    logical, intent(in) :: ifActive     ! If all or only in_air particles considered
    type(silam_species), dimension(:), pointer :: pSpecies
    integer, intent(out) :: nSpecies
    real(r8k), dimension(:), optional, intent(out) :: moles

    ! Local variables
    integer :: i, ix, iy, iz, iSrc, iSubst, iSpecies
    real, dimension(:), pointer :: massVectPtr
    real, dimension(:, :), pointer :: pLPmass
    integer, dimension(:,:,:), pointer :: iSp
    type(Tmass_map), pointer :: pMap
    type(Tlagrange_particles_set), pointer :: pLpSet

    !
    ! Get the names from the cocktail template
    !
    select case(iSpeciesSelector)
      case(species_emission)
        pSpecies => cloud%speciesEmission
        pMap => cloud%mapEmis
        nSpecies = cloud%nSpEmission
        pLpmass => null()
      case(species_transport)
        pSpecies => cloud%speciesTransport
        pMap => cloud%mapConc
        nSpecies = cloud%nSpTransport
        pLpmass => cloud%lpset%lpMassTrn
      case(species_short_lived)
        pSpecies => cloud%speciesShortLived
        pMap => cloud%mapShortlived
        nSpecies = cloud%nSpShortLived
        pLpmass => cloud%lpset%lpMassSL
      case(species_aerosol)
        pSpecies => cloud%speciesAerosol
        pMap => cloud%mapAerosol
        nSpecies = cloud%nSpAerosol
        pLpmass => cloud%lpset%lpMassAER
      case(species_optic_density_3D)
        pSpecies => cloud%speciesOpticalDensity
        pMap => cloud%mapOptDns
        nSpecies = cloud%nSpOpticalDns
        pLpmass => null()
      case(species_column_optic_depth_2D)
        pSpecies => cloud%speciesOpticalDensity
        pMap => cloud%mapOptColDepth
        nSpecies = cloud%nSpOpticalDns
        pLpmass => null()
      case default
        call msg('Unknown type of species:',iSpeciesSelector)
        call set_error('Unknown type of species','get_cloud_inventory')
        return
    end select
    if(error)return
    !
    ! Should the moles be present, get them by summing-up the cloud particles
    ! or integrating the concentration map
    !
    if(present(moles))then 
      moles(1:nSpecies) = 0.
      if (defined(pMap)) then
        if (defined(pMap)) then
      do iy = 1, pMap%ny
       do ix = 1, pMap%nx
        do iz = 1, pMap%n3d
         do iSrc = 1, pMap%nSrc
           do iSpecies = 1, nSpecies
            moles(iSpecies) = moles(iSpecies) + pMap%arM(iSpecies, iSrc, iz, ix, iy)
           end do
         end do   ! iSrc
        end do   ! iz
       end do   ! ix
      end do   ! iy
        endif
      endif

      if (cloud%lpset%defined == silja_true ) then
         if (associated(pLPmass))then
            do ix = 1,cloud%lpset%nop
              do iSpecies = 1, nSpecies
               moles(iSpecies) = moles(iSpecies) + pLpmass(iSpecies, ix)
              end do
            enddo
         endif
      endif

    endif  ! moles are present

  end subroutine get_cloud_inventory


  !***********************************************************************

  integer function fu_subst_nbr_in_inventory(cloud, name)
    !
    ! Checks if the nuclide is available from full inventory
    ! and returns its index there.
    !
    ! ATTENTION. It is believed that all radioactive lagrange particles
    ! have exactly the same list of nuclides in decay_chain (not in the
    ! source list, however).
    !
    implicit none

    ! Imported parameters with intent IN
    TYPE(silam_pollution_cloud), INTENT(in) :: cloud
    character(len=*), intent(in):: name

    ! Local variables
    integer :: iSpecies
!    character(len=clen), dimension(max_nuclides) :: names

!    fu_subst_nbr_in_inventory = fu_iSubst(name, cloud%cocktailTransp%mapper)
    fu_subst_nbr_in_inventory = 0
    do iSpecies = 1, cloud%nSpTransport
      if(fu_substance_name(cloud%speciesTransport(iSpecies)) == name)then
        fu_subst_nbr_in_inventory = iSpecies
      endif
    end do

  end function fu_subst_nbr_in_inventory

  ! Encapsulation for the interpolation coefficients
  !*******************************************************************************************************************

  function fu_ifMeteo2DispHorizInterp(cloud) result(ifTrue)
    implicit none
    type(silam_pollution_cloud) :: cloud
    logical :: ifTrue
    ifTrue = cloud%ifMeteo2DispHorizInterp
  end function fu_ifMeteo2DispHorizInterp

  function fu_ifMeteo2DispVertInterp(cloud) result(ifTrue)
    implicit none
    type(silam_pollution_cloud) :: cloud
    logical :: ifTrue
    ifTrue = cloud%ifMeteo2DispVertInterp
  end function fu_ifMeteo2DispVertInterp
  
  function fu_interpCoefMeteo2DispHoriz(cloud) result(intStruct)
    implicit none
    type(silam_pollution_cloud) :: cloud
    type(THorizInterpStruct), pointer :: intStruct
    intStruct => cloud%interpCoefMeteo2DispHoriz
  end function fu_interpCoefMeteo2DispHoriz

  function fu_interpCoefMeteo2DispVert(cloud) result(intStruct)
    implicit none
    type(silam_pollution_cloud) :: cloud
    type(TVertInterpStruct), pointer :: intStruct
    intStruct => cloud%interpCoefMeteo2DispVert
  end function fu_interpCoefMeteo2DispVert

  ! ***************************************************************
  ! ***************************************************************
  !
  !
  !                  Information on the cloud
  !
  !
  ! ***************************************************************
  ! ***************************************************************


  subroutine check_cloud_masses(pCld, chPlace, fMassThreshold)
    !
    ! Checks all the cloud structures for presence of negative masses
    !
    implicit none

    ! Imported parameters
    TYPE(silam_pollution_cloud), intent(inout) :: pCld
    character(len=*), intent(in) :: chPlace
    real, dimension(:), pointer :: fMassThreshold

    ! Local variables
    integer :: ix, iy, iLev, iSrc, kTmp, iSpecies, iThread, nThreads
    logical :: ifError
    real :: fMass
    real, dimension(:), pointer :: fWork
    real, dimension(:,:), pointer  :: myGarbage

    if (.not. associated(fMassThreshold)) then
      call set_error("low_mass_threshold is not set","check_cloud_masses")
      return
    endif

    !$ fWork => fu_work_array()
    if (error) return

    !$OMP PARALLEL  default(none) &
    !$OMP & PRIVATE (ix, iy, iLev, iSrc, iSpecies, kTmp, ifError, fMass, myGarbage, iThread, nThreads) &
    !$OMP & SHARED (pCld,  chPlace, error,  fMassThreshold, fWork )
    !
    ! Eulerian environment
    !
    iThread = 0
    nThreads=1
    !$ nThreads = omp_get_num_threads()
    !$ iThread = omp_get_thread_num()
    if(defined(pCld%mapConc))then
      if (iThread == 0) then
          call msg("Checking cloud masses at:"+ chPlace &
            !$ &+", OMP nThreads:", nThreads  &
             & )
          myGarbage => pCld%garbage_mass(:,:)
      else
          iSrc=pCld%mapConc%nSrc
          iSpecies=pCld%nSpTransport
          iLev=iThread*iSrc*iSpecies+1
          kTmp=(iThread+1)*iSrc*iSpecies
          myGarbage(1:iSrc,1:iSpecies) => fWork(iLev:kTmp)
          myGarbage(:,:) = 0.0
      endif

      ifError = .false.
      !$OMP DO collapse(4)
      do iy=1,pCld%mapConc%ny
        do ix=1,pCld%mapConc%nx
          do iLev = 1, pCld%mapConc%n3D  ! Vertical levels
           do iSrc = 1, pCld%mapConc%nSrc   ! emission sources
            if (error) cycle
            do iSpecies = 1, pCld%nSpTransport
             fMass = pCld%mapConc%arM(iSpecies,iSrc,iLev,ix,iy)
             if(.not. fMass >= 0.0)then
               if(-fmass < fMassThreshold(iSpecies))then
                 myGarbage(iSrc,iSpecies) = myGarbage(iSrc,iSpecies) + fMass
                 pCld%mapConc%arM(iSpecies,iSrc,iLev,ix,iy) = 0.
                 cycle
               endif
               !$OMP critical (bark_cm)
               call msg('Negative mass found at:' + chPlace + ', ix,iy', ix,iy)
               call msg('Level and source:',iLev,iSrc)
               call msg('Species, index, mass:' + fu_str(pCld%speciesTransport(iSpecies)), &
                                          & iSpecies, pCld%mapConc%arM(iSpecies,iSrc,iLev,ix,iy))
               !$OMP end critical (bark_cm)
               ifError = .true.
             endif
            end do
            !
            ! If negative mass is found in the cell, report it in more details
            !
            if(ifError)then
              !$OMP  critical (bark_cm)
              do iSpecies = 1, pCld%nSpTransport
                do kTmp = 1,pCld%mapConc%n3D
                  call msg('Level:', kTmp)
                  call msg('iSpecies, mass:',iSpecies, &
                           & pCld%mapConc%arM(iSpecies,iSrc,kTmp,ix,iy))
                end do
              end do
              call set_error('Negative mass found at:' + chPlace,'check_cloud_masses')
              call msg('')
              !$OMP end critical (bark_cm)
              cycle
            endif  ! error

          end do  ! iSrc  
         end do  ! iLev
        end do  ! ix
      end do  ! iy
      !$OMP END DO
    endif ! Eulerian structures present
    !
    ! Lagrnagian environment
    !
    if(fu_true(pCld%lpSet%defined))then
      ifError = .false.
      !$OMP DO
      do ix = 1, pCld%lpSet%nop     ! Lagrangian particles
        if (error) cycle
        do iSpecies = 1, pCld%nSpTransport
          
          if(.not. pCld%lpSet%lpMassTrn(iSpecies,ix) >= 0.0)then
            iSrc = mod(pCld%lpSet%lpStatus(ix), 100)  ! WTF? Should be defined in constants...
            if(-pCld%lpSet%lpMassTrn(iSpecies,ix) < fMassThreshold(iSpecies))then
              myGarbage(iSrc,iSpecies) = myGarbage(iSrc,iSpecies) + &
                                               & pCld%lpSet%lpMassTrn(iSpecies,ix)
              pCld%lpSet%lpMassTrn(iSpecies,ix) = 0.
              cycle
            endif
            !$OMP critical (bark_cm)
            call msg('Negative mass found at:' + chPlace + ', Lagr particle', ix)
            call msg('Particle nbr and source:',ix,iSrc)
            call msg('Species, index, mass:' + fu_str(pCld%speciesTransport(iSpecies)), &
                                                  & iSpecies, pCld%lpSet%lpMassTrn(iSpecies,ix))
            !$OMP end critical (bark_cm)
            ifError = .true.
          endif
        end do
        !
        ! If negative mass is found in the cell, report it in more details
        !
        if(ifError)then
          call set_error('Negative mass found at:' + chPlace,'check_cloud_masses')
          cycle
        endif  ! error
      end do  ! particles
      !$OMP END DO  
    endif  ! Lagrangian structures present

    !$OMP BARRIER
    if(defined(pCld%mapConc))then
       !$OMP MASTER
       do iThread = 1,nThreads-1
          iSrc=pCld%mapConc%nSrc
          iSpecies=pCld%nSpTransport
          iLev=iThread*iSrc*iSpecies+1
          kTmp=(iThread+1)*iSrc*iSpecies
          myGarbage(1:iSrc,1:iSpecies) => fWork(iLev:kTmp)
          pCld%garbage_mass(:,:) = pCld%garbage_mass(:,:) +  myGarbage(:,:)
       enddo
       !$OMP end MASTER 
     endif

    !$OMP END PARALLEL
    !$ call free_work_array(fWork)

    if(error)then
      call msg('For references: the list of transport species')
      do ix=1, pCld%nSpTransport
        call msg('The next species number:',ix)
        call report(pCld%speciesTransport(ix))
      end do
    endif

  end subroutine check_cloud_masses

  
  !*************************************************************************
  
  subroutine check_cloud_mass_centres(pCld, chPlace, tolerance_factor)
    !
    ! For Eulerian environment, checks that the centres of masses are in right places
    ! For general case, allow 0.1% of tolerance, i.e. 0.50005 will be considered as a problem
    !
    implicit none
    TYPE(silam_pollution_cloud), intent(inout) :: pCld
    character(len=*), intent(in) :: chPlace
    real, intent(in) :: tolerance_factor

    if(fu_fails(defined(pCld%mapPx_conc),'Not defined Px mass map','check_cloud_mass_centres'))return
    call check_mass_centres_mass_map(pCld%mapConc, pCld%mapPx_conc, chPlace, "Px_concentration", tolerance_factor)
    call check_mass_centres_mass_map(pCld%mapConc, pCld%mapPy_conc, chPlace, "Py_concentration", tolerance_factor)
    call check_mass_centres_mass_map(pCld%mapConc, pCld%mapPz_conc, chPlace, "Pz_concentration", tolerance_factor)

  end subroutine check_cloud_mass_centres
    

  !********************************************************************************************

  subroutine check_species(pCld, chPlace, fMassThreshold, chSp_lowcnc, chSp_highcnc)
    !
    ! Checks all the cloud structures for presence of low masses in higher amounts than the high masses
    !
    implicit none

    ! Imported parameters
    TYPE(silam_pollution_cloud), intent(inout) :: pCld
    character(len=*), intent(in) :: chPlace
    real, dimension(:), intent(in) :: fMassThreshold
    character(len=*), intent(in) :: chSp_lowcnc, chSp_highcnc

    ! Local variables
    integer :: ix, iy, iLev, iSrc, kTmp, iSpecies, iThread, nThreads
    logical :: ifError
    integer ::  iSp_lowcnc, iSp_highcnc
    real :: fMassHi, fMassLo, thHi, thLo
    real, dimension(:), pointer :: fWork
    real, dimension(:,:), pointer  :: myGarbage

    ifError = .false.
    iSp_lowcnc = int_missing
    iSp_highcnc = int_missing
    do iSpecies = 1, pCld%nSpTransport
!      call msg(fu_str(pCld%speciesTransport(iSpecies)) + '--- versus ---' + chSp_lowcnc)
      if(fu_str(pCld%speciesTransport(iSpecies)) == chSp_lowcnc) iSp_lowcnc = iSpecies
      if(fu_str(pCld%speciesTransport(iSpecies)) == chSp_highcnc) iSp_highcnc = iSpecies
    end do
!    call msg('Check-species. low-cnc: '+chSp_lowcnc+', vs high-cnc:' + chSp_highcnc+'. Indices found:',iSp_lowcnc,iSp_highcnc)
    if(iSp_lowcnc == int_missing .or. iSp_highcnc == int_missing)then
      call set_error('Species not found:'+ chSp_lowcnc+'-is:'+fu_str(iSp_lowcnc) + ',' + &
                    & chSp_highcnc + '-is:' + fu_str(iSp_highcnc),'check_species')
      return
    endif
    if(iSp_lowcnc < 1 .or. iSp_lowcnc > size(pCld%mapConc%arM,1) .or. &
     & iSp_highcnc < 1 .or. iSp_highcnc > size(pCld%mapConc%arM,1))then
      call set_error('Strange index:'+ chSp_lowcnc+'-is:'+fu_str(iSp_lowcnc) + ',' + &
                    & chSp_highcnc + '-is:' + fu_str(iSp_highcnc),'check_species')
      return
    endif

    if(size(fMassThreshold) < pCld%nSpTransport)then
      call msg('Strange size of low-mass threshold:',size(fMassThreshold))
      return
    endif

    thHi =  0.1 * fMassThreshold(iSp_highcnc) 
    thLo =  0.1 * fMassThreshold(iSp_lowcnc) 
    
    !$ fWork => fu_work_array()
    if (error) return

    !$OMP PARALLEL  default(none) &
    !$OMP & PRIVATE (ix, iy, iLev, iSrc, iSpecies, kTmp, ifError, fMassHi, fMassLo, myGarbage, iThread, nThreads) &
    !$OMP & SHARED (pCld,  chPlace, error, thHi, thLo, fWork, iSp_highcnc, iSp_lowcnc, &
    !$OMP         & chSp_lowcnc, chSp_highcnc)
    !
    ! Eulerian environment
    !
    iThread = 0
    nThreads=1
    !$ nThreads = omp_get_num_threads()
    !$ iThread = omp_get_thread_num()
    if (iThread == 0) then
!!       call msg("Checking cloud species at:"+ chPlace &
!!         !$ &+", OMP nThreads:", nThreads  &
!!          & )
       myGarbage => pCld%garbage_mass(:,:)
    else
       iSrc=pCld%mapConc%nSrc
       iSpecies=pCld%nSpTransport
       iLev=iThread*iSrc*iSpecies+1
       kTmp=(iThread+1)*iSrc*iSpecies
       myGarbage(1:iSrc,1:iSpecies) => fWork(iLev:kTmp)
       myGarbage(:,:) = 0.0
    endif

    if(defined(pCld%mapConc))then
      ifError = .false.
      !$OMP DO collapse(3)
      do iy=1,pCld%mapConc%ny
        do ix=1,pCld%mapConc%nx
         do iLev = 1, pCld%mapConc%n3D  ! Vertical levels
          if (error) cycle 
          do iSrc = 1, pCld%mapConc%nSrc   ! emission sources
!            call msg('Here (iSrc,iLev,ix,iy)',(/iSrc,iLev,ix,iy/))
            fMassHi  = pCld%mapConc%arM(iSp_highcnc,iSrc,iLev,ix,iy)
            fMassLo = pCld%mapConc%arM(iSp_lowcnc,iSrc,iLev,ix,iy)
!            if(pCld%mapConc%arM(iSp_highcnc,iSrc,iLev,ix,iy) + 0.1 * fMassThreshold(iSp_highcnc) < &
!             & pCld%mapConc%arM(iSp_lowcnc,iSrc,iLev,ix,iy) - 0.1 * fMassThreshold(iSp_lowcnc))then
             if ( fMassHi + thHi < fMassLo - thLo ) then
               !$OMP critical (bark_cm)
               call msg('Low-cnc is higher than high-cnc at:' + chPlace + ', ix,iy,iLev,iSrc', &
                                                                              & (/ix,iy,iLev,iSrc/))
               call msg('mass-low-species(' + chSp_lowcnc + '), threshold_low_cnc, mass-high-species('+ &
                              & chSp_highcnc + '), thresh_high_cnc:', &
                              & (/fMassLo, 10*thLo, fMassHi, 10*thHi/) )
               call msg('Mass excess of low-mass species is sent to garbage:', fMassLo - fMassHi)
               !$OMP end critical (bark_cm)
               myGarbage(iSrc,iSp_lowcnc) = myGarbage(iSrc,iSp_lowcnc) + fMassLo - fMassHi
               pCld%mapConc%arM(iSp_lowcnc,iSrc,iLev,ix,iy) =  fMassHi
!               ifError = .true.
            endif
            !
            ! If wrong relation is found in the cell, report it in more details
            !
            if(ifError)then
              !$OMP critical (bark_cm)
              do iSpecies = 1, pCld%nSpTransport
                do kTmp = 1,pCld%mapConc%n3D
                  call msg('Level:', kTmp)
                  call msg('iSpecies, mass:',iSpecies, &
                           & pCld%mapConc%arM(iSpecies,iSrc,kTmp,ix,iy))
                end do
              end do
              call set_error('Wrong low-high relation found at:' + chPlace,'check_species')
              !$OMP end critical (bark_cm) 
              exit
            endif  ! error
          end do  ! iSrc
         end do  ! iLev
        end do  ! ix
      end do  ! iy
    endif ! Eulerian structures present
    !
    ! Lagrnagian environment
    !
    if(fu_true(pCld%lpSet%defined))then
      !$OMP DO  
      do ix = 1, pCld%lpSet%nop     ! Lagrangian particles
        if (error) cycle
        fMassHi = pCld%lpSet%lpMassTrn(iSp_highcnc,ix)
        fMassLo = pCld%lpSet%lpMassTrn(iSp_lowcnc,ix)
!        if(pCld%lpSet%lpMassTrn(iSp_highcnc,ix) + fMassThreshold(iSp_highcnc) < &
!         & pCld%lpSet%lpMassTrn(iSp_lowcnc,ix) - fMassThreshold(iSp_lowcnc))then
        if ( fMassHi + thHi < fMassLo - thLo ) then
           !$OMP critical (bark_cm)
           call msg('Low-cnc is higher than hgigh-cnc at:' + chPlace + ', ix,iy', ix,iy)
           call msg('Level and source:',iLev,iSrc)
           call msg('mass-low-species (' + chSp_lowcnc + '), mass-high-species (' + &
                  & chSp_highcnc + ')', fMassLo, fMassHi)
           !$OMP end critical (bark_cm)
!           ifError = .true.
        endif
        !
        ! If negative mass is found in the cell, report it in more details
        !
        if(ifError)then
          !$OMP critical (bark_cm)
          do iSpecies = 1, pCld%nSpTransport
            do kTmp = 1,pCld%mapConc%n3D
              call msg('Level:', kTmp)
              call msg('iSpecies, mass:',iSpecies, &
                           & pCld%mapConc%arM(iSpecies,iSrc,kTmp,ix,iy))
            end do
          end do
          !$OMP end critical (bark_cm)
          call set_error('Wrong low-high relation found at:' + chPlace,'check_species')
        endif  ! error
      end do  ! particles
    endif  ! Lagrangian structures present
    !$OMP BARRIER
    !$OMP MASTER
    do iThread = 1,nThreads-1 ! except for the 0-th one
       iSrc=pCld%mapConc%nSrc
       iSpecies=pCld%nSpTransport
       iLev=iThread*iSrc*iSpecies+1
       kTmp=(iThread+1)*iSrc*iSpecies
       myGarbage(1:iSrc,1:iSpecies) => fWork(iLev:kTmp)
       pCld%garbage_mass(:,:) = pCld%garbage_mass(:,:) +  myGarbage(:,:)
    enddo
    !$OMP end MASTER 

    !$OMP END PARALLEL
    !$ call free_work_array(fWork)

    if(error)then
      call msg('For references: the list of transport species')
      do ix=1, pCld%nSpTransport
        call msg('The next species number:',ix)
        call report(pCld%speciesTransport(ix))
      end do
    endif

  end subroutine check_species
  
  
  !*****************************************************************

  SUBROUTINE print_pollution_cloud_report(cloud, ifDetailed)

    IMPLICIT NONE

    ! Imported parameters with intent IN:
    TYPE(silam_pollution_cloud), intent(in) :: cloud
    logical, intent(in), optional :: ifDetailed

    ! Local variable
    integer :: i

    IF (.NOT.defined(cloud)) THEN
      call msg('Undefined cloud')
      RETURN
    END IF

    call msg('')
    call msg('Pollution cloud report:')
    call msg('Emission species:', cloud%nSpEmission)
    do i = 1, cloud%nSpEmission
      call report(cloud%speciesEmission(i))
    end do
    call msg('')
    call msg('Transport species:')
    do i = 1, cloud%nSpTransport
      call report(cloud%speciesTransport(i))
    end do


    if(cloud%lpSet%defined == silja_true)then

      call msg('Lagrangian cloud part:')
      call msg('Earliest start:' + fu_str(cloud%earliest_start))
      call msg('Latest start:' + fu_str(cloud%latest_start))

    elseif(defined(cloud%mapConc))then ! Eulerian advection

      call msg('Eulerian cloud part')
      if(cloud%ifMeteo2DispHorizInterp)then
        call msg('Interpolation meteo-to-dispersion horizontal: YES')
      else
        call msg('Interpolation meteo-to-dispersion horizontal: NO')
      endif
      if(cloud%ifMeteo2DispVertInterp)then
        call msg('Interpolation meteo-to-dispersion vertical: YES')
      else
        call msg('Interpolation meteo-to-dispersion vertical: NO')
      endif
    endif

!    call collect_total_masses_in_cloud(cloud)
!    call msg('Total masses are collected but not reported')
!    call report_total_masses_of_cloud(cloud, 0, .false.)

  END SUBROUTINE print_pollution_cloud_report


  !**********************************************************************

  subroutine collect_total_masses_in_cloud(cloud, pMasses_in_air)
    !
    ! Collects all masses to the cloud totals. Split for sources
    ! and substances is observed. We need all species, including
    ! chemical and radioactive decay chains.
    !
    implicit none

    ! Imported parameter
    TYPE(silam_pollution_cloud), pointer :: cloud
    real, dimension(:), optional, intent(out) :: pMasses_in_air

    ! Local variables
    integer :: iTmp, nGrid, iSrc, ix, iy, iz, nx, ny, nz
    real, dimension(:), pointer :: massVectPtr
    TYPE(silam_pollution_cloud), pointer :: cld


    if(.not.associated(cloud%mInAir))return ! Just stupidity check...

    cld => cloud
!    nSpecies = fu_nbr_of_species(cloud%cocktailTransp)

    cloud%mInAir = 0.
    cloud%mDryDep = 0.
    cloud%mWetDep = 0.
    cloud%mOut= 0.
    cloud%mIn = 0.
    cloud%mStop= 0.

    !
    ! Collect masses from the concentration mass map
    !
    if(defined(cloud%mapConc))then
      do iSrc = 1, cloud%mapConc%nSrc
        do iz=1, cloud%mapConc%n3D
          do iy = 1, cloud%mapConc%ny
            do ix=1, cloud%mapConc%nx
              cloud%mInAir(iSrc,1:cloud%nSpTransport,iz) = cloud%mInAir(iSrc,:,iz) + &
                                                         & cloud%mapConc%arM(:,iSrc,iz,ix,iy)
            end do
          end do
        end do
      end do
    endif  ! if Eulerian environment is active
    !
    ! Masses in Lagrangian lpSet
    !
    cloud%lpSet%nop_active = 0
    if(fu_true(cloud%lpSet%defined))then
      do ix = 1, cloud%lpSet%nop
        if(cloud%lpSet%lpStatus(ix) /= int_missing)then
          iSrc = mod(cloud%lpSet%lpStatus(ix), 100)
          cloud%mInAir(iSrc,1:cloud%nSpTransport,1) = cloud%mInAir(iSrc,:,1) + cloud%lpSet%lpMassTrn(:,ix)
          cloud%lpSet%nop_active = cloud%lpSet%nop_active + 1
        endif
      end do
    endif
    !
    ! Collect deposition grids
    !
    do iSrc = 1, cloud%mapDryDep%nSrc
      do iy = 1, cloud%mapDryDep%ny
        do ix=1, cloud%mapDryDep%nx
          cloud%mDryDep(iSrc,1:cloud%nSpTransport) = cloud%mDryDep(iSrc,:) + &
                                                        & cloud%mapDryDep%arM(:,iSrc,1,ix,iy)
          cloud%mWetDep(iSrc,1:cloud%nSpTransport) = cloud%mWetDep(iSrc,:) + &
                                                        & cloud%mapWetDep%arM(:,iSrc,1,ix,iy)
        end do
      end do
    end do

    ! Collect the outside masses. Have to resolve the global cases
    !

    do iSrc = 1, cloud%mapDryDep%nSrc
      do iz = 1, size(cloud%x0_mass,4)              ! nz_dispersion or 1 - for Eulerian and Lagrangian
        cloud%mOut(iSrc,1:cloud%nSpTransport) = cloud%mOut(iSrc,:) + &
                                                       & cloud%x0_mass(iSrc,:,outgoing,iz) + &
                                                       & cloud%xM_mass(iSrc,:,outgoing,iz) + &
                                                       & cloud%y0_mass(iSrc,:,outgoing,iz) + &
                                                       & cloud%yM_mass(iSrc,:,outgoing,iz)
        cloud%mIn(iSrc,1:cloud%nSpTransport) = cloud%mIn(iSrc,:) + &
                                                       & cloud%x0_mass(iSrc,:,incoming,iz) + &
                                                       & cloud%xM_mass(iSrc,:,incoming,iz) + &
                                                       & cloud%y0_mass(iSrc,:,incoming,iz) + &
                                                       & cloud%yM_mass(iSrc,:,incoming,iz)
      end do   ! iz
      cloud%mOut(iSrc,1:cloud%nSpTransport) = cloud%mOut(iSrc,:) + &
                                                       & cloud%vert_mass_0(iSrc,:,outgoing) + &
                                                       & cloud%vert_mass_M(iSrc,:,outgoing)
      cloud%mIn(iSrc,1:cloud%nSpTransport) = cloud%mIn(iSrc,:) + &
                                                       & cloud%vert_mass_0(iSrc,:,incoming) + &
                                                       & cloud%vert_mass_M(iSrc,:,incoming)
      cloud%mStop(iSrc,1:cloud%nSpTransport) = cloud%mStop(iSrc,:) + cloud%garbage_mass(iSrc,:)

    end do  ! iSrc

    if(present(pMasses_in_air))then
      do iTmp = 1, cloud%nSpTransport
        pMasses_in_air(iTmp) = sum(cloud%mInAir(1:cloud%mapDryDep%nSrc, iTmp, 1:cloud%mapConc%n3D))
      end do
    endif
    
  end subroutine collect_total_masses_in_cloud


  !**********************************************************************

  subroutine report_total_masses_of_cloud(cloud, iStep, ifActive)
    !
    ! Prints total masses of all species for all sources consisting the cloud
    !
    implicit none

    ! Imported parameters
    TYPE(silam_pollution_cloud), INTENT(in) :: cloud
    integer, intent(in) :: iStep
    logical, intent(in) :: ifActive

    ! Local variables
    integer :: iSrc, nSrc, nSp, nop_active, bufsz
    logical :: ifSP, ifNP
    real, dimension(:,:,:), pointer :: vals
    real, dimension(:), pointer :: work
    logical :: ifOk

    if(error)return

    ifSP= .false.
    ifNP= .false.
    if (associated(cloud%boundaryStructures(southern_boundary)%ptr)) &
      & ifSP = (cloud%boundaryStructures(southern_boundary)%ptr%boundaryType == polar_boundary_type)
    
    if (associated(cloud%boundaryStructures(northern_boundary)%ptr)) &
      & ifNP = (cloud%boundaryStructures(northern_boundary)%ptr%boundaryType == polar_boundary_type)
    
    nop_active = int_missing
    if(fu_true(cloud%lpSet%defined)) nop_active = cloud%lpSet%nop_active


    !
    ! What do we report?
    !
    nSrc = size(cloud%mInAir,1)
    nSp = cloud%nSpTransport 
    bufsz = 9*nSrc*nSp

    work => fu_work_array(bufsz*2) !Also reserve for WHOLE_MPI stuff
    work(1:bufsz*2) = 0
    
    vals(1:nSrc,1:nSp,1:9)    => work(1:bufsz)


    !Collect masses
    vals(1:nSrc,1:nSp,1) = sum(cloud%mInAir(1:nSrc,1:nSp,:),3) !In Air
    vals(1:nSrc,1:nSp,2) = cloud%mDryDep(1:nSrc,1:nSp)        !DryDep
    vals(1:nSrc,1:nSp,3) = cloud%mWetDep(1:nSrc,1:nSp)        !WetDep
    if(ifSP)then
      do iSrc=1,nSrc
        vals(iSrc,1:nSp,4) = sum(cloud%boundaryStructures(southern_boundary)%ptr% &
                                & polarcapMassMom(1:nSp,iSrc,1,1:nz_dispersion), 2)
      enddo
    endif
    if(ifNP)then
      do iSrc=1,nSrc
        vals(iSrc,1:nSp,5) = sum(cloud%boundaryStructures(northern_boundary)%ptr% &
                                & polarcapMassMom(1:nSp,iSrc,1,1:nz_dispersion), 2)
      enddo
    endif
    vals(1:nSrc,1:nSp,6) =  cloud%mOut(1:nSrc,1:nSp) !Outgoing
    vals(1:nSrc,1:nSp,7) =  cloud%mIn(1:nSrc,1:nSp) !Incoming
    vals(1:nSrc,1:nSp,8) =  cloud%mStop(1:nSrc,1:nSp) !Garbage

    ! All except for incoming mass: it is already in-air
    vals(1:nSrc,1:nSp,9) = sum(vals(1:nSrc,1:nSp,1:6),3) + vals(1:nSrc,1:nSp,8)

    !
    ! Now write the stuff
    !
    
    if (smpi_adv_tasks > 1) then !Report this subdomain separately
      call report_total_mass("THIS_SUBDOMAIN MASS REPORT", vals, cloud%speciesTransport, &
                  & nSrc, nSp, iStep, nop_active)
      call smpi_reduce_add(work(1:bufsz), work(bufsz+1:2*bufsz), 0, smpi_adv_comm, ifOk)
      vals(1:nSrc,1:nSp,1:9) => work(bufsz+1:2*bufsz) 
    endif
    
    ! Whole mpi report
    if(smpi_adv_rank == 0) then
      call report_total_mass("TOTAL MASS REPORT", vals, cloud%speciesTransport, &
                   & nSrc, nSp, iStep, nop_active)
    endif

    call free_work_array(work)
    
    !call report_total_n_cb4()

  contains
  subroutine report_total_mass(title, vals, species, nSrc, nSp, iStep, nop_active)
    implicit none
    character(len=*), intent(in) :: title
    real, dimension(:,:,:), intent(in) :: vals !(iSrc,iSp,iVal)
    type(silam_species), dimension(:), intent(in) :: species
    integer, intent(in) :: nSrc, nSp
    integer, intent(in) :: iStep, nop_active

    integer :: iUnit, idx1, idx2, iSpecies
    integer, dimension(2) :: fUnits
    character (len=fnlen) :: sp ! Not too short
      
    idx1=3 !If no SP print to WetDep
    idx2=6 !If no NP print from Outgoing

    !Header
    sp = ' Species        InAir        DryDep       WetDep'
    if (sum(sum(vals(:,:,4),1),1) > 0.) then 
      sp = trim(sp) // '       SouthPole'
      idx1 = 4  !Print up to SP
    endif

    if (sum(sum(vals(:,:,5),1),1) > 0.) then
      sp = trim(sp) // '     NorthPole'
      idx2 = 5 !Print from Np
    endif
    sp = trim(sp) // '     Outgoing     -Incoming-    Garbage      Total'


    call msg('')
    !!call msg('========== TOTAL MASS REPORT ============')
    call msg('========== '//title//' ============')
    if (nop_active >= 0) then
      call msg('Model time step, active particles: ',iStep, nop_active)
    else
      call msg('Model time step: ',iStep)
    endif

    fUnits(1:2) = (/run_log_funit, 6/)

    do iUnit = 1,2
      do iSrc=1,nSrc
        if(nSrc > 1) write(fUnits(iUnit), '(A,i3)')'Source: ',iSrc
        write(fUnits(iUnit), '(A)') trim(sp)
        do iSpecies=1, cloud%nSpTransport
          write(fUnits(iUnit), '(A15,9(E12.5,1x))') fu_str(species(iSpecies)),& 
                         & vals(iSrc,iSpecies,1:idx1), vals(iSrc,iSpecies,idx2:9)
        end do
      enddo
      
      if(nSp * nSrc > 1)then
        write(fUnits(iUnit),*)'Grand total:'
        write(fUnits(iUnit), '(15x,9(E12.5,1x))') &
                sum(sum(vals(:,:,1:idx1),1),1), sum(sum(vals(:,:,idx2:9),1),1)
      endif

      if (smpi_global_rank /= 0) exit !No print
    end do  ! iUnits
    call msg('---------- END OF '//title//' ----------')


  end subroutine report_total_mass
!!$    
!!$    subroutine report_total_n_cb4()
!!$      real :: wet, dry, in_air, out, in, garb, factor
!!$      integer, dimension(9), parameter :: &
!!$           & indices = (/ind_no, ind_no2, ind_no3, ind_n2o5, ind_hno3, ind_hono, ind_pna, ind_pan, ind_xo2n/)
!!$      integer :: i, ind
!!$
!!$      call msg('Nitrogen report follows...')
!!$
!!$      wet = 0.0
!!$      dry = 0.0
!!$      in_air = 0.0
!!$      garb = 0.0
!!$      out = 0.0
!!$
!!$      do i = 1, size(indices)
!!$        ind = indices(i)
!!$        if (ind == ind_n2o5) then
!!$          factor = 2.0
!!$        else
!!$          factor = 1.0
!!$        end if
!!$        in_air = in_air + cloud%mInAir(1,ind)*factor
!!$        dry = dry + cloud%mDryDep(1,ind)*factor
!!$        wet = wet + cloud%mWetDep(1,ind)*factor
!!$        out = out + cloud%mOut(1,ind)*factor
!!$        in = in + cloud%mIn(1,ind)*factor
!!$        garb = garb + cloud%mStop(1,ind)*factor
!!$      end do
!!$    
!!$      print '(A)', ' Species      InAir      DryDep      WetDep     Outside    Garbage      Total'
!!$      write(run_log_funit, '(A)') &
!!$           & ' Species      InAir      DryDep      WetDep     Outside    Garbage      Total'
!!$      print '(A8,6(E10.5,1x))', 'total N', in_air, dry, wet, out, garb, in_air+dry+wet+out+garb
!!$      write(run_log_funit, '(A,6(E10.5,1x))') 'Total N', in_air, dry, wet, out, garb, in_air+dry+wet+out+garb
!!$
!!$    end subroutine report_total_n_cb4

  end subroutine report_total_masses_of_cloud

  !*************************************************************************************

  subroutine report_inout_mass_cld(cloud, incout)
    !
    ! Prints the detailed info on incoming/outgoing masses of all species for all sources consisting the cloud
    !
    implicit none

    ! Imported parameters
    TYPE(silam_pollution_cloud), INTENT(in) :: cloud
    integer, intent(in) :: incout

        ! Local variables
    integer :: i, iz, iSrc, nSrc, iUnit, iSpecies, nbr_file_units, nCols
    integer, dimension(2) :: fUnits
    logical :: ifSP, ifNP, ifLG
    real, dimension(6) :: pTmp
!    real :: massSP, massNP, fTmp1, fTmp2
!    character(len=6) :: chTmp1, chTmp2, chTmp3, chTMp4
    character(len=6), dimension(6)  :: colheads

    call msg('')
    select case(incout)
      case(incoming) 
        call msg('========== INCOMING MASS REPORT ============')
      case(outgoing)
        call msg('========== OUTGOING MASS REPORT ============')
      case default
        call set_error("Wrong incoming/outgoing index:"+fu_str(incout),"report_inout_grid_mass_cld" )
        return
    end select

    nSrc = size(cloud%mInAir,1)

    ifSP= .false.
    ifNP= .false.
    if (associated(cloud%boundaryStructures(southern_boundary)%ptr)) &
      & ifSP = (cloud%boundaryStructures(southern_boundary)%ptr%boundaryType == polar_boundary_type)
    
    if (associated(cloud%boundaryStructures(northern_boundary)%ptr)) &
      & ifNP = (cloud%boundaryStructures(northern_boundary)%ptr%boundaryType == polar_boundary_type)

    ifLG = fu_ifLonGlobal(dispersion_grid)

    nCols=0
    if (.not. ifLG) then
       colheads(nCols+1) =  'x<1   '
       colheads(nCols+2) =  'x>xMax'
       nCols=nCols+2
    endif

    if( .not. ifSP)then
      colheads(nCols+1) = 'y<1   '
      nCols = nCols+1
    endif
    if( .not. ifNP)then
      colheads(nCols+1) = 'y>yMax'
      nCols = nCols+1
    endif
    colheads(nCols+1) =  'z<1   '
    colheads(nCols+2) =  'z>zMax'
    nCols=nCols+2

    
    fUnits(1:2) = (/run_log_funit, 6/)
    if (smpi_global_rank == 0) then
      nbr_file_units = 2
    else
      nbr_file_units = 1
    end if
    
    do iUnit = 1, nbr_file_units
      do iSrc=1,nSrc
        if(nSrc > 1) write(fUnits(iUnit), *) 'Source: ',iSrc

!        if(ifSP .or. ifNP)then   ! if any pole exists, longitude closure is mandatory
!          write(fUnits(iUnit), '(5A)')' Species        ', chTmp1, '     ' ,chTmp2, '      z<1       z>zMax'
!        else
!          write(fUnits(iUnit), '(5A)')' Species       x<1       x>xMax      ', chTmp1, '     ' ,chTmp2, '      z<1       z>zMax'
!        endif
         write(fUnits(iUnit), '(A14,'+fu_str(nCols)+'100(1x,A10))') ' Species        ',colheads(1:nCols)

        do iSpecies = 1, cloud%nSpTransport
          
            nCols=0
            if (.not. ifLG) then
               pTmp(nCols+1) =  sum(cloud%x0_mass(iSrc,iSpecies,incout,:))
               pTmp(nCols+2) =  sum(cloud%xM_mass(iSrc,iSpecies,incout,:))
               nCols=nCols+2
            endif

            if( .not. ifSP)then
              pTmp(nCols+1) =  sum(cloud%y0_mass(iSrc,iSpecies,incout,:))
              nCols = nCols+1
            endif
            if( .not. ifNP)then
              pTmp(nCols+1) = sum(cloud%yM_mass(iSrc,iSpecies,incout,:))
              nCols = nCols+1
            endif
            pTmp(nCols+1) =  cloud%vert_mass_0(iSrc,iSpecies,incout)
            pTmp(nCols+2) =  cloud%vert_mass_M(iSrc,iSpecies,incout)
            nCols=nCols+2

            write(fUnits(iUnit), '(A14,'+fu_str(nCols)+'(1x,E10.5))') &
                   & fu_str(cloud%speciesTransport(iSpecies)),  pTmp(1:nCols)
           !! call msg("",pTmp(1:nCols))
            !!write(fUnits(iUnit), '(A14, 100(E10.5,1x),E10.5)')
           !! fu_str(cloud%speciesTransport(iSpecies)),  pTmp(1:nCols)

!            fTmp1 = sum(cloud%y0_mass(iSrc,iSpecies,incout,:))
!            fTmp2 = sum(cloud%yM_mass(iSrc,iSpecies,incout,:))
!  
!  
!            if(ifSP .or. ifNP)then        ! pole closure active only if longitude is global
!              write(fUnits(iUnit), '(A14,3(E10.5,1x),E10.5)')&
!                   & fu_str(cloud%speciesTransport(iSpecies)), &
!                   & fTmp1, &
!                   & fTmp2, &
!                   & cloud%vert_mass_0(iSrc,iSpecies,incout), &
!                   & cloud%vert_mass_M(iSrc,iSpecies,incout)
!            else
!              write(fUnits(iUnit), '(A14,5(E10.5,1x),E10.5)')&
!                   & fu_str(cloud%speciesTransport(iSpecies)), &
!                   & sum(cloud%x0_mass(iSrc,iSpecies,incout,:)), &
!                   & sum(cloud%xM_mass(iSrc,iSpecies,incout,:)), &
!                   & fTmp1, &
!                   & fTmp2, &
!                   & cloud%vert_mass_0(iSrc,iSpecies,incout), &
!                   & cloud%vert_mass_M(iSrc,iSpecies,incout)
!            endif
          end do
!  
!          if(cloud%nSpTransport > 3 .and. nSrc > 1)then
!            pTmp(1:6) = 0.0
!            if(error)return
!            pTmp(1) = pTmp(1) + sum(cloud%x0_mass(iSrc,1:cloud%nSpTransport,incout,:))
!            pTmp(2) = pTmp(2) +  sum(cloud%xM_mass(iSrc,1:cloud%nSpTransport,incout,:))
!            pTmp(3) = pTmp(3) +  sum(cloud%y0_mass(iSrc,1:cloud%nSpTransport,incout,:))
!            pTmp(4) = pTmp(4) +  sum(cloud%yM_mass(iSrc,1:cloud%nSpTransport,incout,:))
!            pTmp(5) = pTmp(5) +  sum(cloud%vert_mass_0(iSrc,1:cloud%nSpTransport,incout))
!            pTmp(6) = pTmp(6) +  sum(cloud%vert_mass_M(iSrc,1:cloud%nSpTransport,incout))
!            if (ifNP) then
!              pTmp(3) = 0.
!              do iz = 1, nz_dispersion
!                pTmp(3) = pTmp(3) + sum(cloud%boundaryStructures(northern_boundary)%ptr%polarCapMassMom(:,iSrc,1,iz)) !* &
!  !                      & cloud%boundaryStructures(northern_boundary)%ptr%PolarCapArea * &
!  !                      & cloud%boundaryStructures(northern_boundary)%ptr%PolarCapThickness(iz)
!              end do
!            endif
!            if (ifSP) then
!              pTmp(4) = 0.
!              do iz = 1, nz_dispersion
!                pTmp(4) = pTmp(4) + sum(cloud%boundaryStructures(southern_boundary)%ptr%polarCapMassMom(:,iSrc,1,iz)) !* &
!  !                       & cloud%boundaryStructures(southern_boundary)%ptr%PolarCapArea * &
!  !                       & cloud%boundaryStructures(southern_boundary)%ptr%PolarCapThickness(iz)
!              end do
!            end if
!            write(fUnits(iUnit),*)'Source total:'
!            write(fUnits(iUnit),'(15x,6(E10.5,1x))') pTmp(1:6)
!          endif
      end do  !iSrc
    end do !iUnit
    if (incout == incoming) then 
      call msg('---------- END OF INCOMING MASS REPORT ----------')
    else
      call msg('---------- END OF OUTGOING MASS REPORT ----------')
    endif
    call msg('')

  end subroutine report_inout_mass_cld



  !==============================================================================================
  !
  ! The emission processors: used with D/A, these operate on the
  ! emission map before it is injected into the tranport mass map.
  ! The actual operations are defined in the source_apportionment module.

  function fu_boundaryStructures(cloud) result(bStructs)
    implicit none
    type(silam_pollution_cloud), intent(in) :: cloud
    type(TboundaryStructPtr), dimension(:), pointer :: bStructs
    bStructs => cloud%boundaryStructures
  end function fu_boundaryStructures

  function fu_boundaryBuffer(cloud) result(bBuffer)
    implicit none
    type(silam_pollution_cloud), intent(in) :: cloud
    type(TboundaryBuffer), pointer :: bBuffer
    bBuffer => cloud%pBoundaryBuffer
  end function fu_boundaryBuffer

  function fu_emission_processor_ptr(cloud) result(proc_ptr)
    implicit none
    type(silam_pollution_cloud), intent(in) :: cloud
    type(Temission_processor), pointer :: proc_ptr
    proc_ptr => cloud%emission_processor
  end function fu_emission_processor_ptr

  subroutine set_emission_processor_ptr(cloud, proc_ptr)
    implicit none
    type(silam_pollution_cloud), intent(inout) :: cloud
    type(Temission_processor), pointer :: proc_ptr
    cloud%emission_processor => proc_ptr
  end subroutine set_emission_processor_ptr

  subroutine set_3_sources(nlptr)
    implicit none
    type(Tsilam_namelist), pointer :: nlptr

    force_3sources = .false.
    
    if (.not. associated(nlptr)) then
      return
    end if

    if (fu_content(nlptr, 'method') == '4D') then
      call set_da_source_indices(1)
      force_3sources = .true.
      call msg('Setting number of sources 4d-var:', DA_NUM_SOURCES)
    else
      force_3sources = .false.
      call set_da_source_indices(1)
    end if

  end subroutine set_3_sources

  subroutine init_tla(cloud, chemrules, traj, num_steps)
    implicit none
    type(silam_pollution_cloud), intent(in) :: cloud
    type(Tchem_rules), intent(in) :: chemrules
    type(t_tla_trajectory), intent(out) :: traj
    integer, intent(in) :: num_steps

    call get_tla_chemistry(chemrules, traj)
    if (error) return

    if (.not. defined(traj)) return ! no need to allocate
    call alloc_tla_trajs(traj, cloud%mapConc%nx, cloud%mapConc%ny, cloud%mapConc%n3d, &
                       & cloud%mapConc%nSpecies, num_steps)
    
  end subroutine init_tla

END MODULE pollution_cloud
