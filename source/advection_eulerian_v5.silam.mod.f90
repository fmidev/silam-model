MODULE advection_eulerian_v5

  ! Description:
  ! Here are tools for calculating Lagrangian and Eulerian advection - transport by
  ! wind.
  !
  ! Authors: Mikhail Sofiev, Julius Vira, Rostislav Kouznetsov, FMI, firstname.lastname@fmi.fi
  !
  ! All units: SI
  ! 
  ! All degrees: real numbers in degrees (degrees and parts per
  ! hundred/parts per thousand , NOT in degrees and minutes
  ! 
  ! Language: ANSI standard Fortran 90
  ! 
  ! Modules used:

  use depositions
  use ini_boundary_conditions
  use chemistry_manager
  use silam_partitioning

  !$ use OMP_LIB

!#define SKIP_VERTICAL
!define SKIP_X_ADVECTION
!define SKIP_Y_ADVECTION

#ifdef MASS_DISTRIBUTOR
#warning 'MASS_DISTRIBUTOR directive is depricated. use "MASS_DISTRIBUTOR" in standard setup instead!'
#endif
!!#define DEBUG
!!!#define DEBUG_V5

  IMPLICIT NONE

  private

  public InitEulerAdvectionFields_v5
  public euler_adv_v5_input_needs
  public report_inout_mass_stuff_v5
  public adv_euler_Galp_xy_v5
  public adv_diffusion_vertical_v5 !both abs and non-abs
      ! z- and eta- verticals, can work with R_down_meteo
      !Diffuses mixing rations
  public test_advect_massMAS
  public test_advect_massR
  public test_advect_mass2
  public test_advect_cell
  public test_advect_mass_many

  private init_molec_diff
  private make_molec_diffusion_passengers
  !Horizontal helpers
  private encode_MPI_moving_masses
  private decode_MPI_moving_masses

  ! Vertical helpers
  private prepareDiffCMvert
  private prepareDiffColumnPressure
  private make_diffusion_passengers

  ! Universal helpers
  private advect_cellboundaries
  private MASS_DISTRIBUTOR
  private advect_mass
  private advect_mass_trislab
  private advect_mass_step
  private get_linestuff_from_boundaries
  private get_linestuff_from_mm
  private count_incoming_cells
  !  private put_linestuff_back

  integer, parameter, public :: advect_rect = 100020
  integer, parameter, public :: advect_tri = 100021
  integer, parameter, public :: advect_step = 100022
 

  !
  ! Actual module variables - metadata, sharable 
  !
 
  logical, private, parameter :: if_report_adv = .false.
  logical, private, save :: do_vert_diff = .true.  ! will be rewritten from fu_kz_param(wdr) if true
  logical, private, save :: do_molec_diff = .false. 
  logical, private, parameter :: report_diff = .false.


  real, private, parameter :: maxmix = 0.005
 ! real, private, parameter :: maxmix = 0.0
  real, private :: CMrelax = 1.0
  integer, private :: distributor_type = int_missing
  logical, private :: diffuse_cm_vert = .false.

  character(10), private :: distributor_name = "Unknown"
 
 
  logical, private, save :: ifXFirst = .true.

  ! Stuff common for v4 vertical and horizontal



  type T_Galp_v5_thread_stuff                            ! will be array 0:nthreads-1
    private
    ! Here nxmax = max(nx_dispersion, ny_dispersion)
    real(r8k), dimension(:), allocatable :: lineParams, lineParamsAdvected !(0:nxmax)
    real, dimension(:), allocatable :: wind_right, scratch, cellmass  !  (0:nxmax)
    logical, dimension(:), allocatable :: ifSkipCell         ! Is there any mass? (0:nxmax)
    real, dimension(:,:,:), allocatable :: cm_line !1:nxy, 1:nspecies, 1:nsources
    real, dimension(:), allocatable ::  moment_line_out !1:nxy
    real, dimension(:,:,:,:), allocatable :: passengers    !0:npassengers, 0:nxMax, 1:nspecies, 1:nsources
    real, dimension(:,:), allocatable :: passengers_out    !0:npassengers, 0:nxMax
    real(r8k), dimension(:,:,:,:), allocatable :: x0_mass, xM_mass, y0_mass, yM_mass   ! 1:nsources, 1:nspecies, 1:2 in/out, 1:nz
    real(r8k), dimension(:,:,:), allocatable :: y0_z_mom, ym_z_mom   ! 1:nsources, 1:nspecies, 1:nz
    real, dimension(:,:), allocatable :: tmp_garbage
    logical, dimension(:,:),  allocatable :: ifSkipSourceSpecies           !1:nspecies, 1:nsources
    real, dimension(:), allocatable :: Cmax1d
    real(r8k), dimension(:), allocatable :: Cmean1d
    integer, dimension(:), allocatable ::Ccount         ! (nz)
    !Vertical-specific 
    real(r8k), dimension(:), allocatable :: PAbove  !Needed for meteo interpolation
    real, dimension(:),   allocatable :: rho_above !(0:nzmax)
    real, dimension(:,:),   allocatable :: vSettling !Used in vertical (1:nSp,0:nzmax)
    real, dimension(:), allocatable :: windRightTmp
    real, dimension(:),   allocatable :: r_down_met, p_met, z_met !(0:nz_meteo+1)
    real, dimension(:),   allocatable :: R_Surf
    real, dimension(:,:),   allocatable :: Vd !(1:nspecies, 1:nsources) 
    real, dimension(:), allocatable   ::   rho_over_R, g_over_deltaP !0:nzMax+1
    real, dimension(:), allocatable :: cm_relax_diff ! (1:nz) CM relaxation factor for subgrid diffusion
    real(r8k), dimension(:), allocatable     ::  A, B, C, P !0:nzMax+1
    real(r8k), dimension(:,:), allocatable   ::  Q !0:npass 0:nzMax+1
    real(r8k), dimension(:,:,:), allocatable ::  bottom_mass, top_mass   ! 1:nsources, 1:nspecies, 1:2 in/out
    type(silja_logical) :: defined = silja_false
  end type T_Galp_v5_thread_stuff
  type T_Galp_v5_vert_shared_stuff
    private
     type(TVertInterpStruct), pointer ::  pVertInterpTop  => null()
     logical,  dimension(:), allocatable :: ifSettles  !nSpecies
     real, dimension(:), pointer :: a_met => null(), b_met => null()
     real, dimension(:), pointer :: a_half_disp  => null(), b_half_disp  => null(), disp_layer_top_m => null()
     logical :: ifMoldiffInitialised = .false.
     real(r8k), dimension(:,:), pointer :: xplus_moldiff  => null(), xminus_moldiff  => null()  !(1:nz, a:nSp), Stuff used for molecular diffusion
  end type T_Galp_v5_vert_shared_stuff

  type T_Eulerian_v5_stuff
    private

   type(T_Galp_v5_thread_stuff), dimension(:), allocatable :: Galp5   !Nthreads
   type(T_Galp_v5_vert_shared_stuff) ::  Galp5zShared !also used in v5
  end type T_Eulerian_v5_stuff 

  type (T_Eulerian_v5_stuff), target, private, save ::  EulerStuff  ! the main internal temporary structures

  real, dimension(:,:,:), allocatable  :: MPI_stuff  ! (nItems, left:right, nThreads)
  integer, dimension(:,:, :), allocatable :: nMPIsend  ! left:right, nThreads, our:their
   integer, dimension(:,:), allocatable :: maxWingUsed ! iWing, iThread

CONTAINS



  ! ****************************************************************

  subroutine euler_adv_v5_input_needs( ifIncludeVerticalDependentFlds, &
                                 & q_met_dynamic, q_met_st, q_disp_dynamic, q_disp_st)
    ! Include the meteorological quantities needed by the eulerian
    ! advection methods.
    implicit none
    logical, intent(in) :: ifIncludeVerticalDependentFlds
    INTEGER, DIMENSION(:), INTENT(inout) :: q_met_dynamic, q_met_st, q_disp_dynamic, q_disp_st
    ! Local variables
    integer :: iTmp


    !call msg("euler_adv_input_needs")    
    ! Check for compatibility: some vertical advection is not compatible with some Kz definitions
    ! Also, some horizontal advections is not compatible with some vertical ones
    !
    ! Wind & pressure:
    ! This peece must stay here until diagnostics_input_needs is done
    iTmp = fu_merge_integer_to_array(u_flag,       q_met_dynamic)
    iTmp = fu_merge_integer_to_array(v_flag,       q_met_dynamic)
!    iTmp = fu_merge_integer_to_array(pressure_flag,q_met_dynamic)
    iTmp = fu_merge_integer_to_array(height_flag,  q_met_dynamic)

    iTmp = fu_merge_integer_to_array(surface_pressure_flag, q_met_dynamic)
    iTmp = fu_merge_integer_to_array(r_down_meteo_flag, q_met_dynamic)
    iTmp = fu_merge_integer_to_array(temperature_flag, q_met_dynamic)
    iTmp = fu_merge_integer_to_array(relative_humidity_flag, q_met_dynamic)

!    if(fu_leveltype(dispersion_vertical) == layer_btw_2_height)then
!         iTmp = fu_merge_integer_to_array(pressure_flag, q_met_dynamic)
!    endif
    iTmp = fu_merge_integer_to_array(disp_flux_celln_rt_flag, q_disp_st)
    iTmp = fu_merge_integer_to_array(disp_flux_cellt_rt_flag, q_disp_st)
    iTmp = fu_merge_integer_to_array(disp_flux_celle_rt_flag, q_disp_st)
    iTmp = fu_merge_integer_to_array(disp_cell_airmass_flag, q_disp_dynamic)
    iTmp = fu_merge_integer_to_array(disp_flux_celltop_flag, q_disp_dynamic)
    iTmp = fu_merge_integer_to_array(disp_flux_celleast_flag, q_disp_dynamic)
    iTmp = fu_merge_integer_to_array(disp_flux_cellnorth_flag, q_disp_dynamic)
!    if(fu_leveltype(dispersion_vertical) .ne. layer_btw_2_hybrid)then
!        !Needed to properly interpolate meteo quantities
        !Always needed 
        iTmp = fu_merge_integer_to_array(height_flag, q_met_dynamic)
!    endif

    !  Not really needed for advection!!!
 !   iTmp = fu_merge_integer_to_array(dispersion_u_flag, q_disp_dynamic)
 !   iTmp = fu_merge_integer_to_array(dispersion_v_flag, q_disp_dynamic)
 !   iTmp = fu_merge_integer_to_array(dispersion_w_flag, q_disp_dynamic)
 !   iTmp = fu_merge_integer_to_array(air_density_flag, q_disp_dynamic)
 !   iTmp = fu_merge_integer_to_array(cell_size_z_flag, q_disp_dynamic)
    if(ifIncludeVerticalDependentFlds)then
      if(.not. defined(dispersion_vertical))then
        call set_error('Dispersion vertical not defined', 'euler_adv_input_needs')
        return
      end if
      if(fu_leveltype(dispersion_vertical) == layer_btw_2_hybrid)then
        iTmp = fu_merge_integer_to_array(cell_size_z_flag, q_disp_dynamic)
        iTmp = fu_merge_integer_to_array(surface_pressure_flag, q_disp_dynamic)
      end if
    end if

   end subroutine euler_adv_v5_input_needs


  !********************************************************************************

  subroutine InitEulerAdvectionFields_v5( nsources, nspecies,  npassengers, &
      & nxy, nz, nthreads, adv_variant, smoother_factor, ifMolecDiffusion, ifSubgridDiffusion )
    !
    ! Initialises the internal advection fields.
    ! It handles only internal fields, which are private for the module.
    !
    implicit none

    ! Imported parameter
    integer, intent(in) :: nsources, nspecies, npassengers, nxy, nz, nthreads, adv_variant
    real, intent(in) :: smoother_factor
    logical, intent(in) :: ifMolecDiffusion, ifSubgridDiffusion

    ! Local variables
    integer :: iLev,iStatus, ithread, nPoles, nxyz
    type (silja_level), dimension(max_levels) :: leveltops
    type (silam_vertical) :: vertdispTop
    
    nxyz = max(nxy, nz)
    !
    ! Stupidity check.
    !
    if(nx_dispersion < 2 .or. ny_dispersion < 2 .or. nx_dispersion > 2000 .or. ny_dispersion > 2000)then
      call msg('nx and ny of dispersion grid:',nx_dispersion, ny_dispersion)
      call msg_warning('Strange dispersion grid parameters','InitEulerAdvectionFields_v5')
    endif
    

    call msg("Allocating EulerStuff%Galp5 for numthreads:", nthreads)
        allocate(EulerStuff%Galp5(0:nthreads-1), stat=iStatus)
        if(fu_fails(iStatus==0,'Failed EulerStuff%Galp5 allocation','InitEulerAdvectionFields_v5'))return

        do ithread = 0, nthreads-1
!          call msg("Allocating for thread:", ithread)
          allocate( &
                & EulerStuff%Galp5(ithread)%lineParams(0:nxyz), &
                & EulerStuff%Galp5(ithread)%lineParamsAdvected(0:nxyz), &
                & EulerStuff%Galp5(ithread)%wind_right(0:nxyz), &
                & EulerStuff%Galp5(ithread)%scratch(0:nxyz), &
                & EulerStuff%Galp5(ithread)%ifSkipCell(0:nxyz+1), &
                & EulerStuff%Galp5(ithread)%cellmass(1:nxyz), &
                & EulerStuff%Galp5(ithread)%cm_line(1:nxyz, 1:nspecies, 1:nsources), &
                & EulerStuff%Galp5(ithread)%moment_line_out(1:nxyz), &
                & EulerStuff%Galp5(ithread)%passengers(0:npassengers+2, 0:nxyz+1, 1:nspecies, 1:nsources), &
                & EulerStuff%Galp5(ithread)%passengers_out(0:npassengers+2, 0:nxyz+1), &
                & EulerStuff%Galp5(ithread)%tmp_garbage(1:nspecies, 1:nsources), &
                & EulerStuff%Galp5(ithread)%x0_mass(1:nsources, 1:nspecies, 1:2, 1:nz_dispersion), &
                & EulerStuff%Galp5(ithread)%xM_mass(1:nsources, 1:nspecies, 1:2, 1:nz_dispersion), &
                & EulerStuff%Galp5(ithread)%y0_mass(1:nsources, 1:nspecies, 1:2, 1:nz_dispersion), &
                & EulerStuff%Galp5(ithread)%yM_mass(1:nsources, 1:nspecies, 1:2, 1:nz_dispersion), &
                & EulerStuff%Galp5(ithread)%y0_z_mom(1:nsources, 1:nspecies, 1:nz_dispersion), &
                & EulerStuff%Galp5(ithread)%yM_z_mom(1:nsources, 1:nspecies, 1:nz_dispersion), &
                & EulerStuff%Galp5(ithread)%ifSkipSourceSpecies(1:nspecies, 1:nsources), &
                & EulerStuff%Galp5(ithread)%Cmax1d(1:nz_dispersion), &
                & EulerStuff%Galp5(ithread)%Cmean1d(1:nz_dispersion), &
                & EulerStuff%Galp5(ithread)%Ccount(1:nz_dispersion), &
                & stat=iStatus)
          if(fu_fails(iStatus==0,'Failed EulerStuff%Galp5('+fu_str(ithread)+') allocation 1',&
                    & 'InitEulerAdvectionFields_v5'))return
          
          EulerStuff%Galp5(ithread)%cm_line(:,:,:) = CONST_NAN
          EulerStuff%Galp5(ithread)%passengers(:,:,:,:) = CONST_NAN
          
             allocate( EulerStuff%Galp5(ithread)%PAbove(0:nz_dispersion), &
                   & EulerStuff%Galp5(ithread)%WindRightTmp(0:nz_dispersion), &
                   & EulerStuff%Galp5(ithread)%rho_above(0:nz_dispersion), &
                   & EulerStuff%Galp5(ithread)%vSettling(1:nspecies,  0:nz_dispersion), &
                   & EulerStuff%Galp5(ithread)%R_surf(1:nspecies), &
                   & EulerStuff%Galp5(ithread)%Vd(1:nspecies, 1:nsources), &
                   & EulerStuff%Galp5(ithread)%r_down_met(0:nz_meteo+1), &
                   & EulerStuff%Galp5(ithread)%p_met(0:nz_meteo+1), &
                   & EulerStuff%Galp5(ithread)%z_met(0:nz_meteo+1), &
                   & EulerStuff%Galp5(ithread)%rho_over_R(0:nz_dispersion+1), &
                   & EulerStuff%Galp5(ithread)%g_over_deltaP(0:nz_dispersion+1), &
                   & EulerStuff%Galp5(ithread)%cm_relax_diff(nz_dispersion), &
                   & EulerStuff%Galp5(ithread)%A(0:nz_dispersion+1), &
                   & EulerStuff%Galp5(ithread)%B(0:nz_dispersion+1), &
                   & EulerStuff%Galp5(ithread)%C(0:nz_dispersion+1), &
                   & EulerStuff%Galp5(ithread)%P(0:nz_dispersion+1), &
                   & EulerStuff%Galp5(ithread)%Q(0:npassengers+2,0:nz_dispersion+1), &
                   & EulerStuff%Galp5(ithread)%bottom_mass(1:nsources, 1:nspecies, 1:2), &
                   & EulerStuff%Galp5(ithread)%top_mass   (1:nsources, 1:nspecies, 1:2), &
                   & stat=iStatus)
             if(fu_fails(iStatus==0,'Failed EulerStuff%Galp5('+fu_str(ithread)+') allocation 2',&
                    & 'InitEulerAdvectionFields_v5'))return
          EulerStuff%Galp5(ithread)%defined = silja_true
        enddo
    !
    ! MPI stuff
    !
    if(any(adv_mpi_neighbours/=int_missing))then
       ! Each thread can fill only 1/nThreads fraction of wings
      allocate(MPI_stuff( & ! 9 data pieces for each species-in-cell
                       & (nsources * nspecies * nz_dispersion * nxy + 1) * &
                       & max(wing_depth_e, wing_depth_n, wing_depth_s, wing_depth_w) * 9 / nThreads, &     
                       & 2, &                         ! left:right
                       & 0:nThreads-1), &              ! OMP threads
             & nMPIsend(2, 0:nThreads-1, 2),& !left:right, OMP threads, our:their
             & maxWingUsed(0:nThreads,4), & !iWing,  OMP threads, (:,nThreads) stores whole-run values
             & stat=iStatus)
             maxWingUsed(:,:) = 0
             
      if(fu_fails(iStatus==0,'MPI stuff allocation failed','InitEulerAdvectionFields_v5'))return
    endif
    !
    ! Init vertical stuff
    !
   if (fu_fails(fu_if_layer(fu_level(dispersion_vertical,1)), &
              & 'Dispersion vertical must be made from layers','InitEulerAdvectionFields_v5'))return
   
   allocate(EulerStuff%Galp5zShared%ifSettles(nSpecies), stat=iStatus)
   if(fu_fails(iStatus==0,'EulerStuff%Galp5zShared%ifSettles allocation failed','InitEulerAdvectionFields_v5'))return


    !Vertical-related arrays
    call vertical_parameters_ab(meteo_vertical, nz_meteo,  &
                       & EulerStuff%Galp5zShared%a_met, EulerStuff%Galp5zShared%b_met)

    call vertical_parameters(dispersion_vertical, nz_dispersion, &
                       & EulerStuff%Galp5zShared%a_half_disp, &
                       & EulerStuff%Galp5zShared%b_half_disp, &
                       & EulerStuff%Galp5zShared%disp_layer_top_m, .true., .true.)

    ! And vertical interpolation structure for level top
    if(.not. associated(EulerStuff%Galp5zShared%pVertInterpTop))then
      call msg('Allocating vertical interpolation for disp_top')
      do iLev = 1, nz_dispersion
        leveltops(iLev) = fu_upper_boundary_of_layer(fu_level(dispersion_vertical, iLev))
      enddo
      call set_vertical(leveltops(1:nz_dispersion), vertdispTop)
      EulerStuff%Galp5zShared%pVertInterpTop => &
               & fu_vertical_interp_struct(meteo_vertical, vertdispTop, dispersion_grid, linear, &
                                             & one_hour, 'meteo_to_disp_leveltop')
      if(error)return
      call set_missing(vertdispTop, .false.) ! Deallocate vertical
    endif
   
    CMrelax = smoother_factor
    distributor_type = adv_variant
    select case (distributor_type)
      case (advect_tri)
            distributor_name="Triangle"
      case (advect_rect)
            distributor_name="Rectangle"
      case (advect_step)
         distributor_name="Step"
      case default
         call set_error("Strange mass distributor type", "InitEulerAdvectionFields_v5")
    end select

    diffuse_cm_vert = ifSubgridDiffusion

    if (ifMolecDiffusion) then
      ! Initialize stuff for molecular diffusion. Currently very primitive:
      ! Use US standard atmosphere and assume diffusion between centers of layers 
      ! As no species here -- allocate only
      if (fu_leveltype(dispersion_vertical) /= layer_btw_2_hybrid) then
        call msg("Can't init molecular diffusion for vertical:")
        call report(dispersion_vertical,.true.)
        call set_error("Molecular diffusion is only iplemented for hybrid layers", &
                  & 'InitEulerAdvectionFields_v5')
         return
      endif
      if(.not. associated(EulerStuff%Galp5zShared%xplus_moldiff))then
        allocate(EulerStuff%Galp5zShared%xplus_moldiff(1:nz, 1:nspecies), &
            & EulerStuff%Galp5zShared%xminus_moldiff(1:nz, 1:nspecies), stat=iStatus)
        if(fu_fails(iStatus==0,'Stuff for molecular diff allocation failed','InitEulerAdvectionFields_v5'))return
        EulerStuff%Galp5zShared%ifMoldiffInitialised = .false. !Yet to be filled with values
      endif
    endif
   
  end subroutine InitEulerAdvectionFields_v5


  !===============================================================================
  !===============================================================================
  !===============================================================================
  !
  !    Eulerian advection.
  !
  !===============================================================================
  !===============================================================================
  !===============================================================================



  
  !************************************************************************************

  subroutine adv_euler_Galp_xy_v5(pDispFlds, &
                                & moment_x, moment_y, moment_z, &
                                & pDispBuf, &
                                & seconds, weight_past, &
                                & garbage, &
                                & x0_mass, xM_mass, y0_mass, yM_mass, &
                                & chem_rules, ifMoments, have_negatives, &
                                & pBBuf, &
                                & ifXfirst)
    !
    ! Each species are taken separately - with own momentum. OMP parallelised.
    !
    use diagnostic_variables, only : iFlux, iMass
    implicit none
    
    ! Imported parameters
    type(Tmass_map), intent(inout) :: pDispFlds
    type(Tmass_map), intent(inout) :: moment_x, moment_y, moment_z 
    type(Tfield_buffer), pointer :: pDispBuf
    real, intent(in) :: seconds, weight_past
    real, dimension(:,:), intent(inout) :: garbage       ! (nSrc,nSpecies)
    real, dimension(:,:,:,:), intent(inout) :: x0_mass, xM_mass, y0_mass, yM_mass  ! (nSrc,nSpecies,nz)
    type(Tchem_rules), intent(in) :: chem_rules
    logical, intent(inout)  :: ifMoments ! In Mass map
    type(TboundaryBuffer), pointer :: pBBuf
    logical, intent(in)  :: ifXfirst, have_negatives


    real(r8k) :: fMTotal   ! Accumulator for mass
    integer :: iDirection, &! Number of cycle over X and Y 
             & iTimesign !Signum of seconds
    integer :: u_ind, v_ind, uadj_ind, vadj_ind, zSize_ind, m_ind ! Field indices in dispersion buffer
    integer ::  isrc, ix, iy, ilev, ispecies, i1d, iPass,  iThread, nThreads, &
             & iCellStartAdv, iCellEndAdv
    integer :: nSrc, nSp, nPass, nLev, iTmp,  &
             & recvCount, recvOffset
    integer, dimension(:), pointer :: iarrptr
    real, dimension(:), pointer :: pXCellSize, pYCellSize,  z_cell_size_past, z_cell_size_future, &
                                 & arMinAdvMass, wind_past, wind_adj, wind_future, cellmass_past, cellmass_future, &
                                 & fptr, exchange_buffer
    real :: fMass, fTmp, fTmp1
    integer, dimension(northern_boundary:western_boundary) :: neighbours, wing_depths ! Communcation neighbours
    logical :: grid_is_lon_glob, ifSkipLine, ifUseMassflux
    integer ::  iExchange, iBorder, iNeighbour, iBoundary, jTmp, wing_off
    
    !Must be variables or parameters to be used in OMP parallel, 
    ! but my_x_coord and my_y_coord can be either depending on mpi or no_mpi
    integer ::  my_x_coord_local, my_y_coord_local 



    type(T_Galp_v5_thread_stuff), pointer :: mystuff, zstuff
    
    integer :: nx_dispersion_mpi, ny_dispersion_mpi
    
    nx_dispersion_mpi = nx_dispersion + wing_depth_e + wing_depth_w
    ny_dispersion_mpi = ny_dispersion + wing_depth_n + wing_depth_s
    !
    !
    wing_depths(northern_boundary:western_boundary) = &
            & (/wing_depth_n, wing_depth_s, wing_depth_e, wing_depth_w/)

    neighbours = adv_mpi_neighbours
    where(neighbours == int_missing) neighbours = smpi_proc_null
    call smpi_get_process_xy_topology(my_x_coord_local, my_y_coord_local, iTmp, jTmp)

    do iBoundary = northern_boundary,western_boundary 
      if ( pBBuf%iBoundaryType(iBoundary) /= smpi_comm_boundary_type &
               & .and. wing_depths(iBoundary) > 0 ) then
         call set_error('Boundary type / halo at '//boundary_name(iBoundary) &
                        &//' mismatch', 'adv_euler_galp_xy_v5')
         return
      end if
    enddo

    iTimesign = sign(1,int(seconds))

    arMinAdvMass => fu_low_mass_threshold(chem_rules)
    nSrc = pDispFlds%nSrc
    nSp = pDispFlds%nSpecies
    nLev = pDispFlds%n3d
    
    iArrPtr => pDispBuf%buffer_quantities
    m_ind = fu_index(iArrPtr, disp_cell_airmass_flag)
    u_ind = fu_index(iArrPtr, disp_flux_celleast_flag)
    v_ind = fu_index(iArrPtr, disp_flux_cellnorth_flag)
    uadj_ind = fu_index(iArrPtr, disp_flux_celle_rt_flag)
    vadj_ind = fu_index(iArrPtr, disp_flux_celln_rt_flag)


    ! Still needed for boundaries with concentrations. 
    ! Should be gone once the boundaries switched to mass mixing ratios
    zSize_ind = fu_index(iArrPtr, cell_size_z_flag)

    if (error .or. any((/u_ind, v_ind, m_ind, uadj_ind, vadj_ind/) < 1)) then
      call msg("u_ind, v_ind, m_ind, uadj_ind, vadj_ind",(/u_ind, v_ind, m_ind, uadj_ind, vadj_ind/))
      call set_error('Failed to get dispersion winds', 'adv_euler_Galp_xy_v5')
      return
    end if

   ! Needed for poles 
    pXCellSize => fu_grid_data(dispersion_cell_x_size_fld)
    pYCellSize => fu_grid_data(dispersion_cell_y_size_fld)

    if(error)then
      call set_error('Cannot get the cell sizes for dispersion grid','adv_euler_Galp_xy_v5')
      return
    endif

    grid_is_lon_glob = fu_ifLonGlobal(pDispFlds%gridTemplate)

    !
    !
    !
    !$OMP PARALLEL DEFAULT(NONE) &
    !$OMP & PRIVATE(ix, iy, ilev, ispecies, i1d,  iPass,  ifSkipLine, mystuff, fTmp, fTmp1, &
    !$OMP         & fMass, wind_adj, wind_past, wind_future, cellmass_past, cellmass_future, &
    !$OMP         & ithread, npass, itmp, iSrc, z_cell_size_past, z_cell_size_future, &
    !$OMP         & iCellStartAdv, iCellEndAdv, fPtr, iExchange, iBorder, iNeighbour, iBoundary, &
    !$OMP         & jTmp, wing_off) &
    !$OMP & SHARED(u_ind, v_ind, uadj_ind, vadj_ind,  zSize_ind, iarrptr, pXCellSize, error, have_negatives, &
    !$OMP         & recvCount, recvoffset, my_x_coord_local, my_y_coord_local, &
    !$OMP        & grid_is_lon_glob, pYCellSize, pDispFlds, moment_x, moment_y, moment_z, CMrelax, &
    !$OMP        & pDispBuf, pBBuf, seconds, weight_past, x0_mass, xM_mass, y0_mass, yM_mass, &
    !$OMP        & garbage, chem_rules, nx_dispersion, ny_dispersion, arMinAdvMass, &
    !$OMP        & ifxfirst, m_ind, ifusemassflux, nsrc, nsp, nLev, ifmoments, EulerStuff, &
    !$OMP        & nThreads, itimesign, MPI_stuff, nx_dispersion_mpi, ny_dispersion_mpi, &
    !$OMP        & wing_depth_w, wing_depth_s, wing_depth_e, wing_depth_n, wing_depth, &
    !$OMP        & nMPIsend, maxWingUsed, neighbours, exchange_buffer, distributor_name, smpi_global_rank)


    iThread = 0
    !$ iThread = OMP_GET_THREAD_NUM()
    mystuff => EulerStuff%Galp5(iThread)

    !$OMP MASTER
    nthreads = 1
    !$ nthreads = omp_get_num_threads()
    if(nthreads == 1)then
      call msg("adv_euler_Galp_xy_v5, " + distributor_name +", CMrelax:",  CMrelax)
    else
      call msg("adv_euler_Galp_xy_v5," + distributor_name +", CMrelax, OMP:",  CMrelax, nThreads)
    endif


    !$OMP END MASTER
    !$OMP BARRIER
    
    mystuff%tmp_garbage = 0.
    mystuff%x0_mass = 0.
    mystuff%xM_mass = 0.
    mystuff%y0_mass = 0.
    mystuff%yM_mass = 0.
    mystuff%y0_z_mom = 0.
    mystuff%yM_z_mom = 0.
    mystuff%Cmax1d = 0.
    mystuff%Cmean1d = 0.
    mystuff%Ccount = 0
    mystuff%lineParamsAdvected(:) = 0. ! If some NaNs there, prediffusion gets mad...

    !Allocated only if neighbours exist
    if (wing_depth >0) maxWingUsed(iThread,:) = 0 !Reset timesteps stats

    do iDirection = 1,2
      !$ if (iDirection == 2) then 
#ifdef DEBUG_V5
      !$ call msg("Barrier") 
#endif
      !$OMP BARRIER !Finish previous round
      !$ endif
      if (error) cycle 

      if(wing_depth > 0) nMPIsend(:, iThread, :) = 0 !Reset send/receive counter

      if ((iDirection == 1) .eqv. ifXfirst ) then    ! if these two logicals are identical
#ifdef  SKIP_X_ADVECTION
        call msg("Skipping X advection")
        cycle
#endif

        if(grid_is_lon_glob)then
            !Make sure that nothing is left in off-grid cells
            mystuff%passengers(:,0,:,:) = 0.
            mystuff%passengers(:,nx_dispersion+1,:,:) = 0.
        endif

        !$OMP DO COLLAPSE (2)
        do iLev = 1, nLev  ! Vertical levels
          do iy = 1, ny_dispersion
            if(error)cycle
            wind_past => pDispBuf%p4d(u_ind)%past%p2d(ilev)%ptr
            wind_future => pDispBuf%p4d(u_ind)%future%p2d(ilev)%ptr
            wind_adj =>  pDispBuf%p4d(uadj_ind)%present%p2d(ilev)%ptr
            cellmass_past   => pDispBuf%p4d(m_ind)%past%p2d(ilev)%ptr
            cellmass_future => pDispBuf%p4d(m_ind)%future%p2d(ilev)%ptr


!call msg('X-adv: iy',iy)
          call get_linestuff_from_mm(pDispFlds, moment_x, moment_y, moment_z, iy, iLev, &
                 &   nx_dispersion, wing_depth_w, nx_dispersion_mpi, &
               &   1, grid_is_lon_glob, western_boundary, eastern_boundary, &
               &   mystuff, nSrc, nSp, have_negatives, arMinAdvMass, pBBuf, ifMoments, ifSkipLine)


          if(ifSkipLine .and. neighbours(western_boundary)==smpi_proc_null .and. & 
             & neighbours(eastern_boundary)==smpi_proc_null) cycle   ! if no masses - skip all

          iCellStartAdv = wing_depth_w + 1
          iCellEndAdv  = wing_depth_w + nx_dispersion

          !
          ! Prepare the line arrays needed for the advect_lineparams, which will move borders of
          ! cells with wind. After that, species masses in slabs will be moved as parts of the cells.
          !
          ix   =  (iy-1) * (nx_dispersion+1) + 1 !staggered base index
          i1d  =  (iy-1) * (nx_dispersion) + 1   ! non-staggered base index

         ! Get line from MPI neighbours
            ! mystuff%cellmass has dimensions of "line"
            if(neighbours(eastern_boundary)/=smpi_proc_null)then
              mystuff%cellmass(wing_depth_w+nx_dispersion+1:nx_dispersion_mpi) = &
                   & pDispBuf%wings_E(iPast,  iMass,iLev,:,iy) *   weight_past + &
                   & pDispBuf%wings_E(iFuture,iMass,iLev,:,iy) * (1.-weight_past)
              mystuff%wind_right(wing_depth_w+nx_dispersion+1:nx_dispersion_mpi) = iTimesign * &
                   & (pDispBuf%wings_E(iPast,    iFlux,iLev,:,iy) *   weight_past + &
                   &  pDispBuf%wings_E(iFuture,  iFlux,iLev,:,iy) * (1.-weight_past) + &
                   &  pDispBuf%wings_E(iRealTime,iFlux,iLev,:,iy))
            end if
            if(neighbours(western_boundary)/=smpi_proc_null)then
              mystuff%cellmass(1:wing_depth_w) = &
                   & pDispBuf%wings_W(iPast,  iMass,iLev,:,iy) *   weight_past + &
                   & pDispBuf%wings_W(iFuture,iMass,iLev,:,iy) * (1.-weight_past)
              mystuff%wind_right(0:wing_depth_w-1) = iTimesign * &
                   & (pDispBuf%wings_W(iPast,    iFlux,iLev,:,iy) *   weight_past + &
                   &  pDispBuf%wings_W(iFuture,  iFlux,iLev,:,iy) * (1.-weight_past) + &
                   &  pDispBuf%wings_W(iRealTime,iFlux,iLev,:,iy))
            end if


            mystuff%cellmass(wing_depth_w+1:nx_dispersion+wing_depth_w) =  & 
                 &   cellmass_past(i1d:i1d+nx_dispersion-1) * weight_past + & !non-staggered
                 & cellmass_future(i1d:i1d+nx_dispersion-1) * (1.-weight_past)

            mystuff%wind_right(wing_depth_w:nx_dispersion+wing_depth_w) = &
                 & iTimesign * (wind_past(ix:ix+nx_dispersion) * weight_past + &
                 & wind_future(ix:ix+nx_dispersion) * (1.-weight_past) + &
                 &    wind_adj(ix:ix+nx_dispersion))
            

         if( .not. grid_is_lon_glob)then
            if (seconds > 0) then
               !Left density
               fTmp = (cellmass_past(i1d) / pDispBuf%p4d(zSize_ind)%past%p2d(ilev)%ptr(i1d)) * weight_past + &
                      & (cellmass_future(i1d) / pDispBuf%p4d(zSize_ind)%future%p2d(ilev)%ptr(i1d)) * (1.-weight_past)
               fTmp =  fTmp / (pXCellSize(i1d)*pYCellSize(i1d)) 
               !Right density
                 iTmp = i1d + nx_dispersion - 1 !Last index of non-staggered line
               fTmp1 = (cellmass_past(iTmp) / pDispBuf%p4d(zSize_ind)%past%p2d(ilev)%ptr(iTmp)) * weight_past + &
                      & (cellmass_future(iTmp) / pDispBuf%p4d(zSize_ind)%future%p2d(ilev)%ptr(iTmp)) * (1.-weight_past)
               fTmp1 =  fTmp1 / (pXCellSize(iTmp)*pYCellSize(iTmp))

               call get_linestuff_from_boundaries(pdispFlds, moment_x, iy, iLev, nx_dispersion_mpi, &
                  & western_boundary, eastern_boundary, iCellStartAdv, iCellEndAdv, &
                  & mystuff%wind_right(0)*seconds,  -mystuff%wind_right(nx_dispersion_mpi)*seconds, &
                  & fTmp, & !airdens_left
                  & fTmp1, &! airdens_right
                  & mystuff, nSrc, nSp, arMinAdvMass, pBBuf, & 
                  & mystuff%x0_mass(:, :, incoming, iLev), mystuff%xM_mass(:, :, incoming, iLev))
              if(error)cycle
           else
               mystuff%passengers(:,0,:,:) = 0.
               mystuff%passengers(:,nx_dispersion_mpi+1,:,:) = 0.
           endif
         else
            ! numerics might break advect_cellboundaries
            if (mystuff%wind_right(0) /= mystuff%wind_right(nx_dispersion_mpi)) then
            !         call msg("wind_right(0) /= wind_right(nx_dispersion_mpi)", &
            !         & mystuff%wind_right(0),  mystuff%wind_right(nx_dispersion_mpi))
            !   call msg_warning("In global domain edge winds must be the same!", "adv_euler_Galp_xy_v5")
            !   call msg("Taking average")
               fTmp = 0.5*(mystuff%wind_right(0) + mystuff%wind_right(nx_dispersion_mpi))
               mystuff%wind_right(0) = fTmp
               mystuff%wind_right(nx_dispersion_mpi) = fTmp
            endif
         endif


           mystuff%ifSkipCell(0:nx_dispersion_mpi+1) = .false. !FIXME hack
          
          
         ! Need these to estmate how much stuff will come from neighbours
         ! Force edges of the domain to be advected
         if (wing_depth_w > 0) mystuff%ifSkipCell(wing_depth_w+1) = .false. 
         if (wing_depth_e > 0) mystuff%ifSkipCell(wing_depth_w + nx_dispersion) = .false. 

          call advect_cellboundaries(mystuff%wind_right(0:nx_dispersion_mpi), &!Staggered
                     & mystuff%cellmass(1:nx_dispersion_mpi),  & ! non-staggered
                     &     mystuff%lineParams(0:nx_dispersion_mpi),  & !Staggered
                     & mystuff%lineParamsAdvected(0:nx_dispersion_mpi),       &!Staggered
                     &     mystuff%ifSkipCell(0:nx_dispersion_mpi+1), & !non-staggered
                     &     mystuff%scratch,   & !staggered
                     &     abs(seconds), grid_is_lon_glob, iCellStartAdv, iCellEndAdv, nx_dispersion_mpi)
           !call msg("")
           !call msg("mystuff%lineParams In (0:)", mystuff%lineParams(0:nx_dispersion_mpi))
           !call msg("mystuff%lineParamsOut ()", mystuff%lineParamsOut(0:nx_dispersion_mpi))
          if(error) then
               call msg("Trouble with X")
               cycle
          endif

          ! Check how many their cells are advected into our domain
          if (wing_depth_e + wing_depth_w > 0)  then
            call  count_incoming_cells(mystuff%lineParamsAdvected(wing_depth_w),  &
                         & mystuff%lineParamsAdvected(wing_depth_w+nx_dispersion), &
                         & mystuff%lineParams(0:nx_dispersion_mpi), &
                         & wing_depth_w,  wing_depth_e, nx_dispersion, iTmp, jTmp)
            if (error) then
             print *, "mystuff%wind_right(0:nx_dispersion_mpi)", mystuff%wind_right(0:nx_dispersion_mpi)
             print *, "mystuff%lineParams In (0:)", mystuff%lineParams(0:nx_dispersion_mpi)
             print *, "mystuff%lineParamsOut (0:)", mystuff%lineParamsAdvected(0:nx_dispersion_mpi)
             cycle
            endif
            nMPIsend(left, iThread, their) = nMPIsend(left, iThread, their) + iTmp
            nMPIsend(right, iThread, their) = nMPIsend(right, iThread, their) + jTmp
            maxWingUsed(iThread, western_boundary) = max(maxWingUsed(iThread,western_boundary),iTmp)
            maxWingUsed(iThread, eastern_boundary) = max(maxWingUsed(iThread,eastern_boundary),jTmp)
          endif

          ! Some Stats (only the current MPI domain)
          !cheap and dirty calculation of courant 
          mystuff%scratch(1:nx_dispersion_mpi) = &
                 & abs(seconds*mystuff%wind_right(1:nx_dispersion_mpi) / mystuff%cellmass(1:nx_dispersion_mpi))
          mystuff%cMax1d(iLev) = max( maxval(mystuff%scratch(1:nx_dispersion_mpi)), &
                                   &  mystuff%cMax1d(iLev))
          mystuff%Cmean1d(iLev) = mystuff%Cmean1d(iLev) &
               & + sum( mystuff%scratch(1:nx_dispersion_mpi))
          mystuff%Ccount(iLev) = mystuff%Ccount(iLev) + nx_dispersion_mpi 
          ! 

          !
          !
          do iSrc = 1, nSrc   ! emission sources
            do ispecies = 1, pdispFlds%nspecies
              if (mystuff%ifSkipSourceSpecies(iSpecies,iSrc)) cycle  ! empty species
              !
              ! Having the mystuff%lineParamsAdvected defined, go species by species "advecting" their slabs following
              ! the cell boundaries template
              !
              nPass = moment_x%passengers(iSpecies)%nPassengers



              call  MASS_DISTRIBUTOR(mystuff%cm_line(1:,iSpecies,iSrc),  mystuff%moment_line_out(1:), &
                          & mystuff%passengers(:,0:,iSpecies,iSrc), mystuff%passengers_out(:,0:),&
                          & nPass+2,  & ! number of passengers + 2 moments
                          & mystuff%lineParams, mystuff%lineParamsAdvected,  mystuff%cellmass,& 
                          & grid_is_lon_glob, iCellStartAdv, iCellEndAdv, nx_dispersion_mpi,&
                          & arMinAdvMass(ispecies), have_negatives, .false., CMrelax) !Don ot convert to CM

              if (error) then 
!                 call  advect_mass_debug(mystuff%cm_line(1:,iSpecies,iSrc),  mystuff%moment_line_out(1:), &
!                          & mystuff%passengers(:,0:,iSpecies,iSrc), mystuff%passengers_out(:,0:),&
!                          & nPass+2,  & ! number of passengers + 2 moments
!                          & mystuff%lineParams, mystuff%lineParamsAdvected,  mystuff%cellmass,& 
!                          & grid_is_lon_glob, iCellStartAdv, iCellEndAdv, nx_dispersion_mpi,&
!                          & arMinAdvMass(ispecies), have_negatives, .false.) !Don ot convert to CM

                call msg("Trouble with advect_mass, x-adv")
                call msg ("iY", iY)
                call msg ("iLev, iSpecies", iLev, iSpecies)
                call msg ("iSrc, iThread", iSrc, iThread)
                cycle
              endif 

              
              !
              ! Return the mass, and moments. Note the MPI shift
              !
              do ix = 1, nx_dispersion
                fMass =  mystuff%passengers_out(0,ix+wing_depth_w)
                if (abs(fMass) <= arMinAdvMass(iSpecies))then
                  !
                  ! Small mass, put to garbage, zero all main arrays
                  !
                  mystuff%tmp_garbage(iSpecies,iSrc) = mystuff%tmp_garbage(iSpecies, iSrc) + fMass
                  pDispFlds%arM(iSpecies,iSrc,iLev,ix,iy) = 0.
                  moment_x%arM(iSpecies,iSrc,iLev,ix,iy) = 0.
                  moment_y%arM(iSpecies,iSrc,iLev,ix,iy) = 0.
                  moment_z%arM(iSpecies,iSrc,iLev,ix,iy) = 0.
                  do iPass = 1, nPass
                    iTmp = moment_x%passengers(iSpecies)%iSpeciesGarbage(iPass)
                    mystuff%tmp_garbage(iTmp,iSrc) =  mystuff%tmp_garbage(iTmp,iSrc) + &
                           &                          mystuff%passengers_out(iPass,ix+wing_depth_w)
                    pDispFlds%arm(iTmp,iSrc,iLev,ix,iy) = 0.
                  end do
                else
                  if (.not. abs(fMass) > 0.) then
                     call msg ("iSpecies,iSrc,iLev,ix,iy,iThread", (/iSpecies,iSrc,iLev,ix,iy,iThread/))
                     call msg ("mass and moments",(/fMass, mystuff%moment_line_out(ix+wing_depth_w), &
                           & mystuff%passengers_out(nPass+1,ix+wing_depth_w), &
                           & mystuff%passengers_out(nPass+2,ix+wing_depth_w)/))
                     call set_error("Gotcha! Wrong mass after X","adv_euler_Galp_xy_v5")
                  endif
                  
                  ! Decent mass. Put it back to main array,  Put all _MOMENTS_ back
                  pDispFlds%arM(iSpecies,iSrc,iLev,ix,iy) = fMass
                  moment_x%arM(iSpecies,iSrc,iLev,ix,iy) = mystuff%moment_line_out(ix+wing_depth_w)
                  moment_y%arM(iSpecies,iSrc,iLev,ix,iy) = mystuff%passengers_out(nPass+1,ix+wing_depth_w) 
                  moment_z%arM(iSpecies,iSrc,iLev,ix,iy) = mystuff%passengers_out(nPass+2,ix+wing_depth_w)

                  do iPass = 1, nPass
                    itmp = moment_x%passengers(iSpecies)%iSpeciesGarbage(iPass)
                    pDispFlds%arm(iTmp,iSrc,iLev,ix,iy) =  mystuff%passengers_out(iPass,ix+wing_depth_w)
                  end do  ! passengers

#ifdef DEBUG_V5
!call msg('After X:('+fu_str(ix)+','+fu_str(iy)+','+fu_str(ilev)+')', &
!       & (/pDispFlds%arm(ispecies, isrc, ilev, ix, iy), moment_x%arm(ispecies, isrc, ilev, ix, iy), &
!         & moment_y%arm(ispecies, isrc, ilev, ix, iy), moment_z%arm(ispecies, isrc, ilev, ix, iy)/))
#endif
                endif  ! low mass
              end do  ! ix
              !
              ! Return the out-of-grid and incoming flows. Do not forget to scale the incoming flow
              ! taking out the mass that was left in the border cell
              !
              if(.not. grid_is_lon_glob)then
                mystuff%x0_mass(iSrc, iSpecies, outgoing, iLev) = &
                       & mystuff%x0_mass(iSrc, iSpecies, outgoing, iLev) + mystuff%passengers_out(0,0)
                mystuff%xM_mass(iSrc, iSpecies, outgoing, iLev) = &
                       & mystuff%xM_mass(iSrc, iSpecies, outgoing, iLev) &
                       & + mystuff%passengers_out(0,nx_dispersion_mpi+1)
              endif  ! global grid
              !
              ! Store the stuff that went out to neighbouring MPI domains
              !
              if(wing_depth_e > 0 .or. wing_depth_w > 0)then
                call encode_MPI_moving_masses(&! masses, x-moment
                                            & mystuff%passengers_out, mystuff%moment_line_out, & 
                                            ! passenger species indexing
                                            & moment_x%passengers(iSpecies)%iSpeciesGarbage, &   
                                            ! what coordinates ?
                                            & iSpecies, iSrc, iLev, int_missing, iy, nPass, &    
                                            & MPI_stuff(:,left,iThread),  nMPIsend(left,iThread, our), &
                                            & MPI_stuff(:,right,iThread), nMPIsend(right,iThread, our))
              endif
            end do ! species
          end do ! iSrc
        end do  ! Cycle iy. X-advection is over
      end do !iLev
      !$OMP END DO
!      !$OMP MASTER
      ifMoments=.true.
!#ifdef DEBUG_V5
!!
!! Nice place to check if something bad is done with centres of masses. Note that they are moments here
!!
!call msg('Checking after X-adv...')
!call check_mass_moments_mass_map(pDispFlds, moment_x, 'after X-adv', 'Px_moment')
!call check_mass_moments_mass_map(pDispFlds, moment_y, 'after X-adv', 'Py_moment')
!call check_mass_moments_mass_map(pDispFlds, moment_z, 'after X-adv', 'Pz_moment')
!call msg('Done checking')
!#endif
!    !$OMP end master
    !
    ! Now or never: X-advection is over, have to communicate MPI stuff to neighbours
    !

!!$    call set_error('After X-advection I want to say to my MPI neighbours','adv_euler_Galp_xy_v4')
!!$    call set_error('Send this:','MPI_stuff(western_boundary, 1:nThreads, 1:maxval(nWest(0:nThreads-1)), 1:9)')
!!$    call set_error('Send this:','MPI_stuff(eastern_boundary, 1:nThreads, 1:maxval(nEast(0:nThreads-1)), 1:9)')
!!$    
!!$    call set_error('After X-advection I want to listen to my MPI neighbours','adv_euler_Galp_xy_v4')
!!$    call set_error('Receive this:','MPI_stuff(western_boundary, 1:nThreads, 1:nWest(1)(intent_in), 1:9)')
!!$    call set_error('Receive this:','MPI_stuff(eastern_boundary, 1:nThreads, 1:nEast(1)(intent_in), 1:9)')
!!$
!!$    call decode_MPI_moving_masses(MPI_stuff(western_boundary,:,:), nWest(1), &
!!$                                & MPI_stuff(eastern_boundary,:,:), nEast(1), &
!!$                                & .true., nThreads, &
!!$                                & pDispFlds, moment_x, moment_y, moment_z)
!!$    stop

       ! error must be synchronous for all threads
       !$OMP BARRIER 
       if (error) cycle  !XY loop
! Do exchange XX
       do iExchange = 1,2
          if (mod(my_x_coord_local+iExchange,2) == 1)  then !Change order -- faster exchange
             if (wing_depth_w < 1 ) cycle
             wing_off = 0
             iBorder = left
             iNeighbour = neighbours(western_boundary)
          else 
             if (wing_depth_e < 1 ) cycle
             wing_off = nx_dispersion-wing_depth_e
             iBorder = right
             iNeighbour = neighbours(eastern_boundary)
          endif

         !$OMP BARRIER 
         !All send/receive numbers must be counted by now
         !$OMP MASTER
           !Worst-scenario receive length

           recvcount  =  nSrc*nSP*9*sum(nMPIsend(iBorder,0:nThreads-1,their)) 
           recvoffset = sum(nMPIsend(iBorder,0:nThreads-1,our)) !Just length to send
           exchange_buffer => fu_work_array(recvoffset + recvcount)
         !$OMP END MASTER
         !$OMP barrier 
          
          iTmp=0 ! Base for this thread in the send buffer
          if (iThread>0) iTmp = sum(nMPIsend(iBorder,0:iThread-1,our))
          jTmp = nMPIsend(iBorder,iThread,our) !Number of items
          exchange_buffer(iTmp+1:iTmp+jTmp) = MPI_stuff(1:jTmp, iBorder, iThread)

          !$OMP barrier
          !Sync before sending
          !$OMP MASTER
!          print *, 'XTASK', smpi_global_rank, 'to', iNeighbour, "recvsize", recvCount
          jTmp = recvOffset !Just length to send
          call smpi_exchange_boundaries(iNeighbour, exchange_buffer, jTmp, recvCount, recvOffset)
          !$OMP END MASTER
          !$OMP barrier
          call  decode_MPI_recv_buffer(exchange_buffer(recvOffset+1:recvCount+recvOffset), &
                               recvCount, wing_off,  0,  iThread, nThreads, &
                                & pDispFlds, moment_x, moment_y, moment_z, ifMoments)
          !$OMP barrier
          !$OMP MASTER
            call free_work_array(exchange_buffer)
          !$OMP END MASTER
        enddo !Exchange


#ifdef DEBUG_V5
   !
   ! Nice place to check if something bad is done with centres of masses. Note that they are moments here
   !
      !$OMP MASTER
   call msg('Checking after X-MPI-decode')
   call check_mass_moments_mass_map(pDispFlds, moment_x, 'after X-adv', 'Px_moment')
   call check_mass_moments_mass_map(pDispFlds, moment_y, 'after X-adv', 'Py_moment')
   call check_mass_moments_mass_map(pDispFlds, moment_z, 'after X-adv', 'Pz_moment')
   call msg('Done checking')
      !$OMP END MASTER
       !$OMP barrier
      !$OMP MASTER
         call msg("Done with X")
      !$OMP END MASTER
#endif

    else !!! Do Y
#ifdef  SKIP_Y_ADVECTION
      call msg("Skipping Y advection")
      cycle
#endif   

!call msg('pBBuf%bNorth before y-adv', pBBuf%bNorth(1,1,1,:))

      !----------------------------------------------------------------------
      !
      ! Y-advection
      ! Never global, possibly with poles or MPI neighbours
      !

      if (pBBuf%iBoundaryType(northern_boundary) == polar_boundary_type .or. &
        &  pBBuf%iBoundaryType(southern_boundary) == polar_boundary_type) then

         !Prepare inflowfactors,  and airmasses and needed for vertical adv
         ! #poles
!call msg('pBBuf%PoleCapThickness(northern_boundary)%pp',pBBuf%PoleCapThickness(northern_boundary)%pp)

        !$OMP DO   
        do iLev = 1, nLev  ! Vertical levels processed independently
          wind_past => pDispBuf%p4d(v_ind)%past%p2d(ilev)%ptr
          wind_future => pDispBuf%p4d(v_ind)%future%p2d(ilev)%ptr
          wind_adj =>  pDispBuf%p4d(vadj_ind)%present%p2d(ilev)%ptr
          cellmass_past   => pDispBuf%p4d(m_ind)%past%p2d(ilev)%ptr
          cellmass_future => pDispBuf%p4d(m_ind)%future%p2d(ilev)%ptr

          !
          ! Treat the poles, if any. Have to ensure that the pole never goes below zero. For that, 
          ! sum-up all outgoing fluxes and attribute them to the amount of stuff in the pole.
          ! From the pole mass > 0 requirement, we get the limitating factor, which varies from 0 (no
          ! mass in the pole) to 1 (unconstrained outflow).
          !
          ! For polar boundary, calculate the limiting factor, otherwise set it to 1.

          do ix = northern_boundary, southern_boundary
             if(pBBuf%iBoundaryType(ix) == polar_boundary_type)then
                i1d = 1   ! Non-staggered index for masses/cell_sizes
                iy = 1    ! Staggered index for v wind
                iTmp = -1 ! to-pole wind sign for v at south
                if (ix == northern_boundary) then 
                   i1d  = nx_dispersion*(ny_dispersion-1)+1 ! Non-staggered index for masses/cell_sizes
                   iy   = nx_dispersion*(ny_dispersion)+1 ! Staggered index for v wind
                   iTmp =  1
                endif

                ! Average area density of air in the layer (kg/m2) around * plar_cap_area
                fMass =     sum(cellmass_past(i1d:i1d+nx_dispersion-1)) *      weight_past + &
                       &  sum(cellmass_future(i1d:i1d+nx_dispersion-1)) *  (1.-weight_past)  !kg

                fMass = fMass * pBBuf%fPoleCapArea(ix)/( pYCellSize(i1d)*pXCellSize(i1d) * nx_dispersion)
                pBBuf%PoleCapAirmass(ix)%pp(iLev) = fMass  ! kg of air

                ! Outflow factor
                fPtr => pBBuf%wind_to_pole(1:nx_dispersion,ix,iLev)
                fPtr = iTmp*(wind_past(iy:iy+nx_dispersion-1)*weight_past + &
                             & wind_future(iy:iy+nx_dispersion-1)*(1.- weight_past) + &
                             & wind_adj(iy:iy+nx_dispersion-1))

                fTmp = - iTimesign * sum(fPtr,  fPtr*iTimesign < 0.) * seconds + 1.e-5  ! to avoid zero
                !Total air flow from pole for the timestep (to pole for backward time) in kg

!                fTmp = max(1.001*fTmp, fMass) 
                 !timesteps / polar_capa
                pBBuf%outflowFactor(ix,iLev) =  fMass / fTmp
             endif

            enddo!Poles
        end do !Levels
        !$OMP END  DO
      !$OMP MASTER
#ifdef DEBUG
   if (pBBuf%iBoundaryType(northern_boundary) == polar_boundary_type) & 
      & call msg("Northern outflow prognosed (timesteps/ cap)",pBBuf%outflowFactor(northern_boundary,:))
   if (pBBuf%iBoundaryType(southern_boundary) == polar_boundary_type) & 
      call msg("Southern outflow prognosed (timesteps/ cap)",pBBuf%outflowFactor(southern_boundary,:))
#endif      
     if (pBBuf%iBoundaryType(northern_boundary) == polar_boundary_type) &
      & pBBuf%outflowFactor(northern_boundary,:) = min(0.9999*pBBuf%outflowFactor(northern_boundary,:),1.)
     if (pBBuf%iBoundaryType(southern_boundary) == polar_boundary_type) &
      & pBBuf%outflowFactor(southern_boundary,:) = min(0.9999*pBBuf%outflowFactor(southern_boundary,:),1.)
      !$OMP END MASTER
      !$OMP BARRIER     !We should finish with polar boundaries before rushing into the main advection

      end if ! #poles
        
      !$OMP DO COLLAPSE (2)
      do iLev = 1, nLev  ! Vertical levels
        do ix = 1, nx_dispersion
          if(error)cycle
          wind_past => pDispBuf%p4d(v_ind)%past%p2d(ilev)%ptr
          wind_future => pDispBuf%p4d(v_ind)%future%p2d(ilev)%ptr
          wind_adj =>  pDispBuf%p4d(vadj_ind)%present%p2d(ilev)%ptr
          z_cell_size_past => pDispBuf%p4d(zSize_ind)%past%p2d(ilev)%ptr
          z_cell_size_future => pDispBuf%p4d(zSize_ind)%future%p2d(ilev)%ptr
          cellmass_past   => pDispBuf%p4d(m_ind)%past%p2d(ilev)%ptr
          cellmass_future => pDispBuf%p4d(m_ind)%future%p2d(ilev)%ptr
        
          call get_linestuff_from_mm(pDispFlds, moment_y, moment_x, moment_z, ix, iLev, &
               &   ny_dispersion, wing_depth_s,  ny_dispersion_mpi, &
               &   2, .false.,  southern_boundary, northern_boundary, &
               &   mystuff, nSrc, nSp, have_negatives, arMinAdvMass, pBBuf, ifMoments, ifSkipLine)

          if(ifSkipLine .and. neighbours(northern_boundary)==smpi_proc_null .and. & 
             & neighbours(southern_boundary)==smpi_proc_null) cycle   ! if no masses - skip all

          !
          ! Prepare the line arrays needed for the advect_lineparams, which will move borders of
          ! cells with wind. After that, species masses in slabs will be moved as parts of the cells.
          !

          iCellStartAdv = wing_depth_s + 1
          iCellEndAdv = wing_depth_s + ny_dispersion

          !
          if(neighbours(northern_boundary)/=smpi_proc_null)then
            mystuff%cellmass(wing_depth_s+ny_dispersion+1:ny_dispersion_mpi) = &
              & pDispBuf%wings_N(iPast,  iMass,iLev,ix,:) * weight_past + &
              & pDispBuf%wings_N(iFuture,iMass,iLev,ix,:) * (1.-weight_past)
            mystuff%wind_right(wing_depth_s+ny_dispersion+1:ny_dispersion_mpi) = iTimesign * &
              & pDispBuf%wings_N(iPast,    iFlux,iLev,ix,:) * weight_past + &
              & pDispBuf%wings_N(iFuture,  iFlux,iLev,ix,:) * (1.-weight_past) + &
              & pDispBuf%wings_N(iRealTime,iFlux,iLev,ix,:)
          end if
          if(neighbours(southern_boundary)/=smpi_proc_null)then
            mystuff%cellmass(1:wing_depth_s) = &
              & pDispBuf%wings_S(iPast,  iMass,iLev,ix,:) * weight_past + &
              & pDispBuf%wings_S(iFuture,iMass,iLev,ix,:) * (1.-weight_past)
            mystuff%wind_right(0:wing_depth_s-1) = iTimesign * &
              & pDispBuf%wings_S(iPast,    iFlux,iLev,ix,:) * weight_past + &
              & pDispBuf%wings_S(iFuture,  iFlux,iLev,ix,:) * (1.-weight_past) + &
              & pDispBuf%wings_S(iRealTime,iFlux,iLev,ix,:)
          end if

          DO iy = 1, ny_dispersion
            i1d = ix + (iy-1) * nx_dispersion
            mystuff%cellmass(wing_depth_s+iy) = &
                 &   cellmass_past(i1d) * weight_past + &
                 & cellmass_future(i1d) * (1.-weight_past)
          end do
          DO iy = 0, ny_dispersion
            i1d = ix + iy * nx_dispersion 
            mystuff%wind_right(wing_depth_s+iy) = iTimesign * &
                 & ( wind_past(i1d) * weight_past + &
                 & wind_future(i1d) * (1.-weight_past) + &
                 &    wind_adj(i1d))
          enddo

!          call msg("wind_past  :   ", wind_past(ix:(nx_dispersion_mpi*(ny_dispersion_mpi+1)):nx_dispersion))
!          call msg("wind_future    :", wind_future(ix:(nx_dispersion_mpi*(ny_dispersion_mpi+1)):nx_dispersion))
!          call msg("wind_adj       :", wind_adj(ix:(nx_dispersion_mpi*(ny_dispersion_mpi+1)):nx_dispersion))
!          call msg("cellmass_past  :", cellmass_past(ix:(nx_dispersion_mpi*ny_dispersion_mpi):nx_dispersion))
!          call msg("cellmass_future:", cellmass_future(ix:(nx_dispersion_mpi*ny_dispersion_mpi):nx_dispersion))
!
!
!          mystuff%wind_right(ny_dispersion_mpi + 1) = -1e27 !Should never be used!!!
!          call msg("Courant number Y:", mystuff%wind_right(1:ny_dispersion_mpi)/mystuff%cellmass(1:ny_dispersion_mpi)*seconds)

           
        if (seconds > 0) then
            !Boundaries might need density
            !Left density
            fTmp = (cellmass_past(ix) / pDispBuf%p4d(zSize_ind)%past%p2d(ilev)%ptr(ix)) * weight_past + &
               (cellmass_future(ix) / pDispBuf%p4d(zSize_ind)%future%p2d(ilev)%ptr(ix)) * (1. -weight_past)
            fTmp =  fTmp / (pXCellSize(ix)*pYCellSize(ix)) 
            !Right density
            iTmp = ix + (ny_dispersion-1)*nx_dispersion  !Last index of non-staggered line
            fTmp1 = (cellmass_past(iTmp) / pDispBuf%p4d(zSize_ind)%past%p2d(ilev)%ptr(iTmp)) * weight_past + &
               (cellmass_future(iTmp) / pDispBuf%p4d(zSize_ind)%future%p2d(ilev)%ptr(iTmp)) * (1. -weight_past)
            fTmp1 =  fTmp1 / (pXCellSize(iTmp)*pYCellSize(iTmp)) !! Air density 

            call get_linestuff_from_boundaries(pdispFlds,moment_y, ix, iLev, ny_dispersion_mpi, &
               & southern_boundary, northern_boundary, iCellStartAdv, iCellEndAdv, &
               &   mystuff%wind_right(0)*seconds,&
               &  -mystuff%wind_right(ny_dispersion_mpi)*seconds, &
               & fTmp, & !airdens_left
               & fTmp1, &! airdens_right
               & mystuff, nSrc, nSp, arMinAdvMass, pBBuf, & 
               & mystuff%y0_mass(:, :, incoming, iLev), mystuff%yM_mass(:, :, incoming, iLev))
         else
            mystuff%passengers(:,0,:,:) = 0.
            mystuff%passengers(:,ny_dispersion_mpi+1,:,:) = 0.
         endif
         if(error)cycle
         mystuff%ifSkipCell(0:ny_dispersion_mpi+1) = .false. !FIXME hack

         ! Need these to estmate how much stuff will come from neighbours
         if (wing_depth_s > 0) mystuff%ifSkipCell(wing_depth_s+1) = .false. 
         if (wing_depth_n > 0) mystuff%ifSkipCell(wing_depth_s + ny_dispersion) = .false. 

          call advect_cellboundaries(mystuff%wind_right(0:ny_dispersion_mpi), &!Staggered
                     & mystuff%cellmass(1:ny_dispersion_mpi),  & ! non-staggered
                     &     mystuff%lineParams(0:ny_dispersion_mpi),  & !Staggered
                     & mystuff%lineParamsAdvected(0:ny_dispersion_mpi),       &!Staggered
                     &     mystuff%ifSkipCell(0:ny_dispersion_mpi+1), & !non-staggered
                     &     mystuff%scratch,   & !staggered
                     &     abs(seconds), .false., iCellStartAdv, iCellEndAdv, ny_dispersion_mpi)

         
          if(error) cycle

          ! Check how many their cells are advected into our domain
         if (wing_depth_s + wing_depth_n > 0)  then
           call  count_incoming_cells(mystuff%lineParamsAdvected(wing_depth_s),  &
                         & mystuff%lineParamsAdvected(wing_depth_s+ny_dispersion), &
                         & mystuff%lineParams(0:ny_dispersion_mpi), &
                         & wing_depth_s,  wing_depth_n, ny_dispersion, iTmp, jTmp)
           nMPIsend(left, iThread, their) = nMPIsend(left, iThread, their) + iTmp
           nMPIsend(right, iThread, their) = nMPIsend(right, iThread, their) + jTmp
           maxWingUsed(iThread, southern_boundary) = max(maxWingUsed(iThread,southern_boundary),iTmp)
           maxWingUsed(iThread, northern_boundary) = max(maxWingUsed(iThread,northern_boundary),jTmp)
         endif

          if(error) then
             call msg("Trouble with Y")
             call msg("wind_past  :   ",        wind_past(ix:(nx_dispersion*(ny_dispersion+1)):nx_dispersion))
             call msg("wind_future    :",     wind_future(ix:(nx_dispersion*(ny_dispersion+1)):nx_dispersion))
             call msg("wind_adj       :",        wind_adj(ix:(nx_dispersion*(ny_dispersion+1)):nx_dispersion))
             call msg("cellmass_past  :",   cellmass_past(ix:(nx_dispersion*ny_dispersion):nx_dispersion))
             call msg("cellmass_future:", cellmass_future(ix:(nx_dispersion*ny_dispersion):nx_dispersion))
             if (wing_depth_s > 0) then
                call msg("Southern wing:")
                call msg("wind_past  :    ",  pDispBuf%wings_S(iPast,    iFlux,iLev,ix,:))
                call msg("wind_future    :", pDispBuf%wings_S(iFuture,  iFlux,iLev,ix,:))
                call msg("wind_adj       :", pDispBuf%wings_S(iRealTime,iFlux,iLev,ix,:))
                call msg("cellmass_past  :", pDispBuf%wings_S(iPast,    iMass,iLev,ix,:))
                call msg("cellmass_future:", pDispBuf%wings_S(iFuture,  iMass,iLev,ix,:))
             endif
             if (wing_depth_n > 0) then
                call msg("Southern wing:")
                call msg("wind_past  :    ",  pDispBuf%wings_N(iPast,    iFlux,iLev,ix,:))
                call msg("wind_future    :",  pDispBuf%wings_N(iFuture,  iFlux,iLev,ix,:))
                call msg("wind_adj       :", pDispBuf%wings_N(iRealTime,iFlux,iLev,ix,:))
                call msg("cellmass_past  :", pDispBuf%wings_N(iPast,    iMass,iLev,ix,:))
                call msg("cellmass_future:", pDispBuf%wings_N(iFuture,  iMass,iLev,ix,:))
             endif
             call msg ("")
             call msg("mystuff%wind_right(0:ny_dispersion_mpi)", mystuff%wind_right(0:ny_dispersion_mpi))
             call msg("mystuff%cellmass(1:ny_dispersion_mpi)", mystuff%cellmass(1:ny_dispersion_mpi))
             call msg("mystuff%lineParams(0:ny_dispersion_mpi)",mystuff%lineParams(0:ny_dispersion_mpi))
             call msg("mystuff%lineParamsAdvected(0:ny_dispersion_mpi)",mystuff%lineParamsAdvected(0:ny_dispersion_mpi))
             cycle
          endif

         
          mystuff%scratch(1:ny_dispersion_mpi) = &
                 & abs(seconds*mystuff%wind_right(1:ny_dispersion_mpi) / mystuff%cellmass(1:ny_dispersion_mpi))

          mystuff%cMax1d(iLev) = max( maxval(mystuff%scratch(1:ny_dispersion_mpi)), &
                                   &  mystuff%cMax1d(iLev))
          mystuff%Cmean1d(iLev) = mystuff%Cmean1d(iLev) &
               & + sum( mystuff%scratch(1:ny_dispersion_mpi))
          mystuff%Ccount(iLev) = mystuff%Ccount(iLev) + ny_dispersion_mpi 
          
          !
          ! Having the mystuff%lineParamsOut defined, go species by species "advecting" their slabs following
          ! the cell boundaries template
          !
!          call msg("")
          do iSrc = 1, nSrc   ! emission sources
            do ispecies = 1, pdispFlds%nspecies
              if (mystuff%ifSkipSourceSpecies(iSpecies,iSrc))cycle  ! empty species
              nPass = moment_y%passengers(iSpecies)%nPassengers

              call  MASS_DISTRIBUTOR(mystuff%cm_line(1:,iSpecies,iSrc),  mystuff%moment_line_out(1:), &
                          & mystuff%passengers(:,0:,iSpecies,iSrc), mystuff%passengers_out(:,0:),&
                          & nPass+2,  & ! number of passengers + 2 moments
                          & mystuff%lineParams, mystuff%lineParamsAdvected,  mystuff%cellmass,& 
                          & .false., iCellStartAdv, iCellEndAdv, ny_dispersion_mpi,&
                          & arMinAdvMass(ispecies), have_negatives, .false., CMrelax) !Don ot convert to CM

              if (error) then 
                call msg("Trouble with advect_mass, y-adv")
                call msg ("iX", iX)
                call msg ("iLev, iSpecies", iLev, iSpecies)
                call msg ("iSrc, iThread", iSrc, iThread)
                call msg("Wind_right", mystuff%wind_right(0:ny_dispersion_mpi))
                cycle
              endif

              !
              ! Return the mass, and moments. Note the MPI shift
              !
              do iy = 1, ny_dispersion
                fMass =  mystuff%passengers_out(0,iy+wing_depth_s)
                if (abs(fMass) <= arMinAdvMass(iSpecies))then
                  !
                  ! Small mass, put to garbage, zero all main arrays
                  !
                  mystuff%tmp_garbage(iSpecies,iSrc) = mystuff%tmp_garbage(iSpecies, iSrc) + fMass
                  pDispFlds%arM(iSpecies,iSrc,iLev,ix,iy) = 0.
                  moment_x%arM(iSpecies,iSrc,iLev,ix,iy) = 0.
                  moment_y%arM(iSpecies,iSrc,iLev,ix,iy) = 0.
                  moment_z%arM(iSpecies,iSrc,iLev,ix,iy) = 0.
                  do iPass = 1, nPass
                    iTmp = moment_y%passengers(iSpecies)%iSpeciesGarbage(iPass)
                    mystuff%tmp_garbage(iTmp,iSrc) = mystuff%tmp_garbage(iTmp,iSrc) + &
                                           & mystuff%passengers_out(iPass,ix)
                    pDispFlds%arm(iTmp,iSrc,iLev,ix,iy) = 0.
                  end do
                else
                  if (.not. abs(fMass) > 0.) then
                    call msg ("iSpecies,iSrc,iLev,ix,iy,iThread", (/iSpecies,iSrc,iLev,ix,iy,iThread/))
                    call msg ("mass and moments",(/fMass, mystuff%moment_line_out(iy+wing_depth_s), &
                           & mystuff%passengers_out(nPass+1,iy+wing_depth_s), &
                           & mystuff%passengers_out(nPass+2,iy+wing_depth_s)/))
                    call set_error("Gotcha! Wrong mass after Y","adv_euler_Galp_xy_v5")
                  endif
                  !
                  ! Decent mass. Put it back to main array, 
                  ! Put all _MOMENTS_ back
                  !
                  pDispFlds%arM(iSpecies,iSrc,iLev,ix,iy) = fMass
                  moment_y%arM(iSpecies,iSrc,iLev,ix,iy) = mystuff%moment_line_out(iy+wing_depth_s)
                  moment_x%arM(iSpecies,iSrc,iLev,ix,iy) = mystuff%passengers_out(nPass+1,iy+wing_depth_s)
                  moment_z%arM(iSpecies,iSrc,iLev,ix,iy) = mystuff%passengers_out(nPass+2,iy+wing_depth_s)

                  do iPass = 1, nPass
                    itmp = moment_y%passengers(iSpecies)%iSpeciesGarbage(iPass)
                    pDispFlds%arm(iTmp,iSrc,iLev,ix,iy) =  mystuff%passengers_out(iPass,iy+wing_depth_s)
                  end do  ! passengers
#ifdef DEBUG_V5
!call msg('After Y:('+fu_str(ix)+','+fu_str(iy)+','+fu_str(ilev)+')', &
!       & (/pDispFlds%arm(ispecies, isrc, ilev, ix, iy), moment_x%arm(ispecies, isrc, ilev, ix, iy), &
!         & moment_y%arm(ispecies, isrc, ilev, ix, iy), moment_z%arm(ispecies, isrc, ilev, ix, iy)/))
#endif
                endif  ! if small mass
              end do  ! iy
              !
              ! Return the out-of-grid flows. Note that for pole we need to return z-moment
              !
              mystuff%y0_mass(iSrc, iSpecies, outgoing, iLev) = &
                                                & mystuff%y0_mass(iSrc, iSpecies, outgoing, iLev) + &
                                                & mystuff%passengers_out(0,0)
              mystuff%yM_mass(iSrc, iSpecies, outgoing, iLev) = &
                                                & mystuff%yM_mass(iSrc, iSpecies, outgoing, iLev) + &
                                                & mystuff%passengers_out(0,ny_dispersion_mpi+1)
              if(pBBuf%iBoundaryType(southern_boundary) == polar_boundary_type) &
                  & mystuff%y0_z_mom(iSrc, iSpecies, iLev) = mystuff%y0_z_mom(iSrc,iSpecies,iLev) + &
                                                & mystuff%passengers_out(nPass+2,0)
              if(pBBuf%iBoundaryType(northern_boundary) == polar_boundary_type) &
                  & mystuff%yM_z_mom(iSrc, iSpecies, iLev) = mystuff%yM_z_mom(iSrc,iSpecies,iLev) + &
                                                & mystuff%passengers_out(nPass+2,ny_dispersion_mpi+1)
              !
              ! Store the stuff that went out to neighbouring MPI domains
              !
              if(ny_dispersion /= ny_dispersion_mpi)then
                call encode_MPI_moving_masses(& ! masses, x-moment
                                            & mystuff%passengers_out, mystuff%moment_line_out, & 
                                            ! passenger species indexing
                                            & moment_y%passengers(iSpecies)%iSpeciesGarbage, &   
                                            ! what coordinates ?
                                            & iSpecies, iSrc, iLev, ix, int_missing, nPass, &    
                                            & MPI_stuff(:,left,iThread),  nMPIsend(left,iThread,our), &
                                            & MPI_stuff(:,right,iThread), nMPIsend(right,iThread,our))
              endif
            end do ! species
          end do ! iSrc
        end do  ! Cycle ix. Y-advection is over
      end do ! Levels : end of Y-cycle of advection
      !$OMP END DO
      
         !$OMP barrier 
            !$OMP MASTER
      ifMoments=.true.
#ifdef DEBUG_V5
!
! Nice place to check if something bad is done with centres of masses. Note that they are moments here
!
call msg('Checking after Y-adv')
call check_mass_moments_mass_map(pDispFlds, moment_x, 'after Y-adv', 'Px_moment')  ! single-thread version should be made...
call check_mass_moments_mass_map(pDispFlds, moment_y, 'after Y-adv', 'Py_moment')
call check_mass_moments_mass_map(pDispFlds, moment_z, 'after Y-adv', 'Pz_moment')
call msg('Done')
#endif
!$OMP END MASTER

       ! error must be synchronous for all threads
       !$OMP BARRIER 
       if (error) cycle

! Do exchange YY
       do iExchange = 1,2
          if (mod(my_y_coord_local+iExchange,2) == 1)  then !Change order -- faster exchange
             if (wing_depth_s < 1 ) cycle
             wing_off = 0
             iBorder = left
             iNeighbour = neighbours(southern_boundary)
          else 
             if (wing_depth_n < 1 ) cycle
             wing_off = ny_dispersion-wing_depth_n
             iBorder = right
             iNeighbour = neighbours(northern_boundary)
          endif

         !$OMP barrier 
         !All send/receive numbers must be counted by now
         !$OMP MASTER
           !Worst-scenario receive length

           recvcount  =  nSrc*nSP*9*sum(nMPIsend(iBorder,0:nThreads-1,their)) 
           recvoffset = sum(nMPIsend(iBorder,0:nThreads-1,our)) !Just length to send
           exchange_buffer => fu_work_array(recvoffset + recvcount)
         !$OMP END MASTER
         !$OMP barrier 
          
          iTmp=0 ! Base for this thread in the send buffer
          if (iThread>0) iTmp = sum(nMPIsend(iBorder,0:iThread-1,our))
          jTmp = nMPIsend(iBorder,iThread,our) ! elements from this thread
          exchange_buffer(iTmp+1:iTmp+jTmp) = MPI_stuff(1:jTmp, iBorder, iThread)

          !$OMP barrier
          !Sync before sending
          !$OMP MASTER
!          print *, 'YTASK', smpi_global_rank, 'to', iNeighbour, "recvOffset", recvOffset, "recvsize", recvCount
          jTmp = recvOffset !Just length to send
          call smpi_exchange_boundaries(iNeighbour, exchange_buffer, jTmp, recvCount, recvOffset)
          !$OMP END MASTER
          !$OMP barrier
          call  decode_MPI_recv_buffer(exchange_buffer(recvOffset+1:recvCount+recvOffset), &
                               recvCount, 0,  wing_off,  iThread, nThreads, &
                                & pDispFlds, moment_x, moment_y, moment_z, ifMoments)
          !$OMP barrier
          !$OMP MASTER
            call free_work_array(exchange_buffer)
          !$OMP END MASTER
        enddo !Exchange

!        write(*,'(A6,I2,A1,6I5)') '+++++ ', smpi_global_rank, ':', lbound(MPI_stuff,1), ubound(MPI_stuff,1), &
!            lbound(MPI_stuff,2), ubound(MPI_stuff,2), &
!            lbound(MPI_stuff,3), ubound(MPI_stuff,3)
!        if(size(MPI_stuff,2) > 1)then
!            call set_error('MPI version does not yet support OpenMP','adv_euler_Galp_xy_v4')
!        end if
!        call smpi_exchange_boundaries(neighbours(southern_boundary), MPI_stuff(southern_boundary,0,:))
#ifdef DEBUG_V5
!
! Nice place to check if something bad is done with centres of masses. Note that they are moments here
!
         !$OMP barrier
          !$OMP MASTER

call msg('Checking after Y-MPI-decode')
call check_mass_moments_mass_map(pDispFlds, moment_x, 'after Y-adv', 'Px_moment')
call check_mass_moments_mass_map(pDispFlds, moment_y, 'after Y-adv', 'Py_moment')
call check_mass_moments_mass_map(pDispFlds, moment_z, 'after Y-adv', 'Pz_moment')
call msg('Done')
      !$OMP END MASTER
#endif

  !   call msg("Done with Y",iCellStartAdv,iCellEndAdv)



    endif !Do Y
    enddo ! Loop over X-Y advection
  
!call msg('pBBuf%bNorth after y-adv', pBBuf%bNorth(1,1,1,:))

  !$OMP CRITICAL (omp_garbage)
  do iSrc = 1, nSrc
    garbage(iSrc, 1:nSp) = garbage(iSrc, 1:nSp) +  mystuff%tmp_garbage(1:nSp,iSrc)
  end do
  !$OMP END CRITICAL (omp_garbage)
  
  
  !$OMP END PARALLEL
    
  ! collect things to the zeroth threads buffer
!call msg('0')
  zstuff => EulerStuff%Galp5(0)
  do  iThread=1,nThreads-1 !Except for the 0-th one
    mystuff => EulerStuff%Galp5(iThread)
!call msg('1',nLev)
    do iLev = 1, nLev
!call msg('2',iLev)
     zstuff%Cmax1d(iLev) = max(mystuff%Cmax1d(iLev), zstuff%Cmax1d(iLev))
    enddo 
!call msg('3')
    zstuff%Cmean1d(:) = zstuff%Cmean1d(:) + mystuff%Cmean1d(:) 
!call msg('4')
    zstuff%Ccount(:) = zstuff%Ccount(:) + mystuff%Ccount(:)
    zstuff%x0_mass(1:nSrc, 1:nSp, 1:2, 1:nLev) = zstuff%x0_mass(1:nSrc, 1:nSp, 1:2, 1:nLev) + & 
                                                      & mystuff%x0_mass(1:nSrc, 1:nSp, 1:2, 1:nLev)
    zstuff%y0_mass(1:nSrc, 1:nSp, 1:2, 1:nLev) = zstuff%y0_mass(1:nSrc, 1:nSp, 1:2, 1:nLev) + & 
                                                      & mystuff%y0_mass(1:nSrc, 1:nSp, 1:2, 1:nLev)
    zstuff%xM_mass(1:nSrc, 1:nSp, 1:2, 1:nLev) = zstuff%xM_mass(1:nSrc, 1:nSp, 1:2, 1:nLev) + & 
                                                      & mystuff%xM_mass(1:nSrc, 1:nSp, 1:2, 1:nLev)
    zstuff%yM_mass(1:nSrc, 1:nSp, 1:2, 1:nLev) = zstuff%yM_mass(1:nSrc, 1:nSp, 1:2, 1:nLev) + & 
                                                      & mystuff%yM_mass(1:nSrc, 1:nSp, 1:2, 1:nLev)
    if(pBBuf%iBoundaryType(southern_boundary) == polar_boundary_type) &
              & zstuff%y0_z_mom(1:nSrc, 1:nSp, 1:nLev) = zstuff%y0_z_mom(1:nSrc, 1:nSp, 1:nLev) + & 
                                                       & mystuff%y0_z_mom(1:nSrc, 1:nSp, 1:nLev)
    if(pBBuf%iBoundaryType(northern_boundary) == polar_boundary_type) &
              & zstuff%yM_z_mom(1:nSrc, 1:nSp, 1:nLev) = zstuff%yM_z_mom(1:nSrc, 1:nSp, 1:nLev) + & 
                                                       & mystuff%yM_z_mom(1:nSrc, 1:nSp, 1:nLev)
    mystuff%x0_mass(1:nSrc, 1:nSp, 1:2, 1:nLev) = 0. ! Keep the budget
    mystuff%xM_mass(1:nSrc, 1:nSp, 1:2, 1:nLev) = 0.
    mystuff%y0_mass(1:nSrc, 1:nSp, 1:2, 1:nLev) = 0.
    mystuff%yM_mass(1:nSrc, 1:nSp, 1:2, 1:nLev) = 0.
    mystuff%y0_z_mom(1:nSrc, 1:nSp, 1:nLev) = 0.
    mystuff%yM_z_mom(1:nSrc, 1:nSp, 1:nLev) = 0.
  enddo

  ! Collect the masses transported in/out during this time step into cumulative ones
  !
  if(.not. grid_is_lon_glob)then
    x0_mass(1:nSrc, 1:nSp, 1:2, 1:nLev) = x0_mass(1:nSrc, 1:nSp, 1:2, 1:nLev) + &
                                         & zstuff%x0_mass(1:nSrc, 1:nSp, 1:2, 1:nLev)
    xM_mass(1:nSrc, 1:nSp, 1:2, 1:nLev) = xM_mass(1:nSrc, 1:nSp, 1:2, 1:nLev) + &
                                         & zstuff%xM_mass(1:nSrc, 1:nSp, 1:2, 1:nLev)
  endif

  if(pBBuf%iBoundaryType(southern_boundary) /= polar_boundary_type)then
      y0_mass(1:nSrc, 1:nSp, 1:2, 1:nLev) = y0_mass(1:nSrc, 1:nSp, 1:2, 1:nLev) + &
                                       & zstuff%y0_mass(1:nSrc, 1:nSp, 1:2, 1:nLev)
  endif
  if(pBBuf%iBoundaryType(northern_boundary) /= polar_boundary_type)then
      yM_mass(1:nSrc, 1:nSp, 1:2, 1:nLev) = yM_mass(1:nSrc, 1:nSp, 1:2, 1:nLev) + &
                                       & zstuff%yM_mass(1:nSrc, 1:nSp, 1:2, 1:nLev)
  endif
  !
  ! Deal with poles: mass and z moment. Note that incoming and outgoing are with regard to grid, not pole
  !
  if(pBBuf%iBoundaryType(southern_boundary) == polar_boundary_type)then

    do iLev = 1, nLev
      do iSrc = 1, nSrc
        do iSpecies = 1, nSp
!          call msg("Southern outflow done (polar timesteps/cap):" + fu_str(iSpecies), iLev  , & 
!               & pBBuf%bSouth(iSpecies,iSrc,1,iLev) * &
!               & pBBuf%fSouthPoleCapArea * pBBuf%SouthPoleCapThickness(iLev)  /  &
!               & zstuff%y0_mass(iSrc, iSpecies, incoming, iLev)  )
          !
          ! cm_new = ((m_ini-d_m_in)*cm_ini + d_mom) / (m_ini - d_m_in + d_m_out)
          ! in and out are with regard to grid, not pole
          !
          fMass = pBBuf%bSouth(iSpecies,iSrc,1,iLev) - zstuff%y0_mass(iSrc, iSpecies, incoming, iLev)
          pBBuf%bSouth(iSpecies,iSrc,1,iLev) = fMass + zstuff%y0_mass(iSrc, iSpecies, outgoing, iLev)

          if(.not. pBBuf%bSouth(iSpecies,iSrc,1,iLev) >= 0.)then
            if(have_negatives)then
              if(.not. pBBuf%bSouth(iSpecies,iSrc,1,iLev) < 0.)then
                call set_error('NaN suspected, southern pole','adv_euler_Galp_xy_v5')
                call msg('pole mass, grid2pole mass, pole2grid mass, limiting factor for species-' &
                              & + fu_str(iSpecies) + &
                              & ', source-' + fu_str(iSrc) + ', level-' + fu_str(iLev), &
                              & (/pBBuf%bSouth(iSpecies,iSrc,1,iLev), &
                                & real(zstuff%y0_mass(iSrc, iSpecies, outgoing, iLev)), &
                                & real(zstuff%y0_mass(iSrc, iSpecies, incoming, iLev)), &
                                & pBBuf%outflowFactor(southern_boundary,iLev) /))
                return
              endif  ! non-negative check: NaN filter
            else
              if(abs(pBBuf%bSouth(iSpecies,iSrc,1,iLev)) > arMinAdvMass(iSpecies))then
                call msg('Negative south pole concentration at iLev, iSpecies, iSrc', (/iLev, iSpecies, iSrc/))
                call msg('Mass, in-pole and out-pole fluxes', (/pBBuf%bSouth(iSpecies,iSrc,1,iLev), &
                                                              & real(zstuff%y0_mass(iSrc, iSpecies, outgoing, iLev)), &
                                                              & real(zstuff%y0_mass(iSrc, iSpecies, incoming, iLev))/))
                call set_error('Negative south pole concentration', 'adv_euler_Galp_xy_v5')
              else
              garbage(iSrc, iSpecies) = garbage(iSrc, iSpecies) + pBBuf%bSouth(iSpecies,iSrc,1,iLev)
              pBBuf%bSouth(iSpecies,iSrc,1,iLev) = 0.
              endif
            endif    ! have_negatives
          endif  ! mass is not positive

          ! Calculate new centre of mass
          !
          if(abs(pBBuf%bSouth(iSpecies,iSrc,1,iLev)) > arMinAdvMass(iSpecies))then
            pBBuf%bSouth(iSpecies,iSrc,2,iLev) = (fMass * pBBuf%bSouth(iSpecies,iSrc,2,iLev) + &
                                                            & zstuff%y0_z_mom(iSrc,iSpecies,iLev)) / &
                                               & pBBuf%bSouth(iSpecies,iSrc,1,iLev)
          else
            garbage(iSrc, iSpecies) = garbage(iSrc, iSpecies) + pBBuf%bSouth(iSpecies,iSrc,1,iLev)
            pBBuf%bSouth(iSpecies,iSrc,1:2,iLev) = 0.0
          endif

        end do  ! iSPecies
      end do ! iSrc
    end do  ! iLev
  endif  ! if south pole boundary

  if(pBBuf%iBoundaryType(northern_boundary) == polar_boundary_type)then
    do iLev = 1, nLev
      do iSrc = 1, nSrc
        do iSpecies = 1, nSp
!          call msg("Northern outflow done (polar timesteps/cap):" + fu_str(iSpecies), iLev  , & 
!               & pBBuf%bnorth(iSpecies,iSrc,1,iLev) * &
!               & pBBuf%fNorthPoleCapArea * pBBuf%NorthPoleCapThickness(iLev)/ &
!               & zstuff%yM_mass(iSrc, iSpecies, incoming, iLev) )
          fMass = pBBuf%bNorth(iSpecies,iSrc,1,iLev) - zstuff%yM_mass(iSrc, iSpecies, incoming, iLev)
          pBBuf%bNorth(iSpecies,iSrc,1,iLev) = fMass + zstuff%yM_mass(iSrc, iSpecies, outgoing, iLev)

          if(.not. pBBuf%bnorth(iSpecies,iSrc,1,iLev) >= 0.)then
            if(have_negatives)then
              if(.not. pBBuf%bnorth(iSpecies,iSrc,1,iLev) < 0.)then
                call set_error('NaN suspected','adv_euler_Galp_xy_v5')
                call msg('pole mass, grid2pole mass, pole2grid mass, limiting factor for species-' &
                              & + fu_str(iSpecies) + &
                              & ', source-' + fu_str(iSrc) + ', level-' + fu_str(iLev), &
                              & (/pBBuf%bnorth(iSpecies,iSrc,1,iLev), &
                                & real(zstuff%yM_mass(iSrc, iSpecies, outgoing, iLev)), &
                                & real(zstuff%yM_mass(iSrc, iSpecies, incoming, iLev)), &
                                & pBBuf%outflowFactor(northern_boundary,iLev) /))
                return
              endif  ! Non-negative check: NaN filter
            else
              if(abs(pBBuf%bNorth(iSpecies,iSrc,1,iLev)) > arMinAdvMass(iSpecies))then
                call msg('Negative north pole concentration at iLev, iSpecies, iSrc', (/iLev, iSpecies, iSrc/))
                call msg('Mass, in-pole and out-pole fluxes', (/pBBuf%bnorth(iSpecies,iSrc,1,iLev), &
                                                              & real(zstuff%yM_mass(iSrc, iSpecies, outgoing, iLev)), &
                                                              & real(zstuff%yM_mass(iSrc, iSpecies, incoming, iLev))/))
                call set_error('Negative north pole concentration', 'adv_euler_Galp_xy_v5')
              else
              garbage(iSrc, iSpecies) = garbage(iSrc, iSpecies) + pBBuf%bnorth(iSpecies,iSrc,1,iLev)
              pBBuf%bNorth(iSpecies,iSrc,1,iLev) = 0.
              endif
            endif  ! have_negatives
          endif   ! mass is not positive

          ! Calculate new centre of mass
          !
          if(abs(pBBuf%bNorth(iSpecies,iSrc,1,iLev)) > arMinAdvMass(iSpecies))then
            pBBuf%bNorth(iSpecies,iSrc,2,iLev) = (fMass * pBBuf%bNorth(iSpecies,iSrc,2,iLev) + &
                                                            & zstuff%yM_z_mom(iSrc,iSpecies,iLev)) / &
                                               & pBBuf%bNorth(iSpecies,iSrc,1,iLev)
          else
            garbage(iSrc, iSpecies) = garbage(iSrc, iSpecies) + pBBuf%bNorth(iSpecies,iSrc,1,iLev)
            pBBuf%bNorth(iSpecies,iSrc,1:2,iLev) = 0.0
          endif
        end do  ! iSPecies
      end do ! iSrc
    end do  ! iLev
  endif  ! if north pole boundary

  if (wing_depth >0) then
!      call msg("Max wing depth used by neighbour This step, NSEW", &
!                &  maxval(maxWingUsed(0:nThreads-1,:), dim=1))
     maxWingUsed(nThreads,:) = maxval(maxWingUsed(0:nThreads,:), dim=1)
      call msg("Max wing depth used by neighbour so far,    NSEW",  maxWingUsed(nThreads,:))
  endif


  call msg('Max Courant for levels (min, mean, max):', (/ minval(zstuff%Cmax1d(:)), sum(zstuff%Cmax1d(:))/nLev, maxval(zstuff%Cmax1d(:)) /))
  call msg('Mean Courant:', sum(zstuff%Cmean1d(:) / (zstuff%Ccount(:)+1e-5))/nLev )  ! to avoid nans
  call msg('')
#ifdef DEBUG_V5
    if(pBBuf%iBoundaryType(northern_boundary) == polar_boundary_type .and. &
         & pBBuf%iBoundaryType(southern_boundary) == polar_boundary_type)then
      do iLev = 1, nLev
        call msg('sum of v-wind to pole at South and North poles at level_ '+fu_str(iLev)+':', &
             & sum(pBBuf%wind_to_pole(1:nx_dispersion, southern_boundary,iLev)), & 
             & sum(pBBuf%wind_to_pole(1:nx_dispersion, northern_boundary,iLev)))
        call msg('sum of abs v-wind at South and North poles at level_ '+fu_str(iLev)+':', &
             & sum(abs(pBBuf%wind_to_pole(1:nx_dispersion, southern_boundary,iLev))), & 
             & sum(abs(pBBuf%wind_to_pole(1:nx_dispersion, northern_boundary,iLev))))
      end do
      call msg('')
    endif
#endif
  
  end subroutine adv_euler_Galp_xy_v5


  !***********************************************************************************

  subroutine get_linestuff_from_mm(pDispFlds, moment_my, moment_other, moment_third, iOther, iThird, &
                                 & nmy_dispersion, shift_myMPI, nmy_dispersion_mpi, &
                                 & advect_axis, ifGlob, left_boundary, right_boundary, &
                                 & mystuff, nSrc, nSp, have_negatives, arMinAdvMass, pBBuf, &
                                 & ifMoments, ifSkipLine)
      !
      ! Replacement for duplicated section in x and y advection
      ! get masses/cm and passengers from massmap and figure out
      ! and fill ifSkipLine, mystuff%ifSkipCell and mystuff%ifSkipSourceSpecies
      ! Should be possible to use in z advection once the stuff is unified
      !

      implicit none

      !Parameters different between axes 
      type(Tmass_map), intent(inout) ::  moment_my, moment_other, moment_third
      integer, intent(in) :: iOther, iThird ! indices in other directions
      integer, intent(in) :: nmy_dispersion, shift_myMPI, nmy_dispersion_mpi !dimensions of advected direction
      integer, intent(in) :: left_boundary, right_boundary ! e.g. east_boundary
      integer, intent(in) :: advect_axis !1 - x, 2 - y, 3 - z
      logical, intent(in) ::  ifGlob, ifMoments
      logical, intent(out) :: ifSkipLine
      ! Parameters same for all axes 
      type (Tmass_map), intent(inout) :: pDispFlds
      type (T_Galp_v5_thread_stuff), intent(inout) :: mystuff
      integer, intent(in) :: nSrc, nSp
      logical, intent(in) :: have_negatives
      real, dimension (:), intent(in) :: arMinAdvMass
      type(TboundaryBuffer), pointer :: pBBuf
            
      !Local
      real :: fMass
      integer, target :: ix, iy, iLev
      integer :: iSrc, iSpecies, iline, iPass, iTmp, nPass
      integer, pointer :: imy
      integer, dimension (:), pointer :: bleftspecies => null(), brightspecies => null()

      integer :: nx_dispersion_mpi, ny_dispersion_mpi

      nx_dispersion_mpi = nx_dispersion + wing_depth_e + wing_depth_w
      ny_dispersion_mpi = ny_dispersion + wing_depth_n + wing_depth_s

      !
      ! Do we need to do this line? Determine iXstart and iXend and flag arrays for species and x
      ifSkipLine = .true.
      mystuff%ifSkipSourceSpecies(1:nSp,1:nSrc) = .true.
      mystuff%ifSkipCell(0:nx_dispersion_mpi) = .true.


      select case  (advect_axis)
          case(1)
             imy => ix
             iy   = iother
             iLev = iThird
           case(2)
             imy => iy
             ix   = iother
             iLev = iThird
           case (3)
             imy => iLev
             ix  = iother
             iy  = iThird
          case default
             call set_error("Wrong axis:"+fu_str(advect_axis),"get_linestuff_from_mm")
!             return
      end select

      DO imy = 1, nmy_dispersion
         iline = imy + shift_myMPI
         do iSrc = 1, nSrc   ! emission sources
           do iSpecies = 1, nSp
             nPass = moment_my%passengers(iSpecies)%nPassengers
             fMass = pDispFlds%arm(ispecies, isrc, ilev, ix, iy)

             if( abs(fMass) <= arMinAdvMass(ispecies))then   !That is bad idea to put negative mass 
                                                           ! to garbage just because it is negative
               !
               ! Small (abs) mass. Send to it garbage together with all  mystuff%passengers, ! nullify moment
               ! 
               mystuff%tmp_garbage(ispecies, iSrc) = mystuff%tmp_garbage(ispecies,isrc) + fMass
               pDispFlds%arM(iSpecies,iSrc,iLev,ix,iy) = 0.
               moment_my%arM(iSpecies,iSrc,iLev,ix,iy) = 0.
               moment_other%arM(iSpecies,iSrc,iLev,ix,iy) = 0.
               moment_third%arM(iSpecies,iSrc,iLev,ix,iy) = 0.
               !
               ! mystuff%passengers are whatever species, same mass map at the end, without moments
               !
               do iPass = 1, nPass
                 iTmp = moment_my%passengers(iSpecies)%iSpeciesGarbage(iPass)
                 mystuff%tmp_garbage(iTmp,iSrc) = &
                           & mystuff%tmp_garbage(iTmp,iSrc) +  pDispFlds%arm(iTmp, iSrc,iLev,ix,iy)
                 pDispFlds%arm(iTmp, iSrc,iLev,ix,iy) = 0.
               end do  ! passengers
               mystuff%cm_line(iline,iSpecies,iSrc) = 0.
               mystuff%passengers(0:nPass+2,iline,iSpecies,iSrc) = 0.
             else
               !
               ! Decent mass. Copy it and its location to line arrays. 
               ! Do not forget: in case of MPI, lines and winds are longer in x-direction than mass map
               ! lines are indexed with iLine
               !
               ifSkipLine = .false. 
               mystuff%ifSkipSourceSpecies(iSpecies,iSrc) = .false.  ! do advection for the species
               mystuff%ifSkipCell(iline) = .false.                      ! do advection for the cell
               mystuff%passengers(0,iline,iSpecies,iSrc) = fMass

               ! Should be centers of mass for x
               ! moments for others are extra  mystuff%passengers
               if (ifMoments) then
                 mystuff%cm_line(iline,iSpecies, iSrc) &
                      & = moment_my%arm(iSpecies, isrc, ilev, ix, iy) / fMass
                 mystuff%passengers(nPass+1,iline,iSpecies,iSrc) &
                      & = moment_other%arM(iSpecies,iSrc,iLev,ix,iy) 
                 mystuff%passengers(nPass+2,iline,iSpecies,iSrc) &
                      & = moment_third%arM(iSpecies,iSrc,iLev,ix,iy)
               else
                 mystuff%cm_line(iline,iSpecies, iSrc) = moment_my%arm(iSpecies, isrc, ilev, ix, iy)
                 mystuff%passengers(nPass+1,iline,iSpecies,iSrc) &
                      & = moment_other%arM(iSpecies,iSrc,iLev,ix,iy) * fMass 
                 mystuff%passengers(nPass+2,iline,iSpecies,iSrc) &
                      & = moment_third%arM(iSpecies,iSrc,iLev,ix,iy) * fMass
               endif
               ! with negative masses centres might get wild...
               if(have_negatives)then
                 if(abs(mystuff%cm_line(iline,iSpecies, iSrc)) >= 0.5) then
                   mystuff%cm_line(iline,iSpecies, iSrc) = 0.0
                 end if
               endif
               !
               ! Get passengers
               do iPass = 1, nPass
                 iTmp = moment_my%passengers(iSpecies)%iSpeciesGarbage(iPass)
                 mystuff%passengers(iPass,iline,iSpecies, iSrc) = pDispFlds%arm(iTmp,iSrc,iLev,ix,iy)
               end do  ! passengers
             endif  ! if mass<threshold

           end do  ! iSpecies
         end do  ! iSrc
      end do   ! imy - taking up the mass

      !
      ! If can get something from boundaries, force calculation of the corresponding edges
      !
      if(.not. ifGlob) then
         !  left (West or South, Bottom)
         iTmp = pBBuf%iBoundaryType(left_boundary)
         if (iTmp == polar_boundary_type) then !All sources
            ifSkipLine = .false.
            mystuff%ifSkipSourceSpecies(:,:) = .false.
            mystuff%ifSkipCell(0) = .false.
         elseif (iTmp == dirichlet_boundary_type) then !Only for iSrc == 1
            where (pBBuf%iBoundarySpecies(:, left_boundary) /= int_missing) &
                                                      & mystuff%ifSkipSourceSpecies(:,1) = .false.
            ifSkipLine = .false.
            mystuff%ifSkipCell(0) = .false.
         endif
          ! right (East or North or Top)
         iTmp = pBBuf%iBoundaryType(right_boundary)
         if (iTmp == polar_boundary_type) then
            ifSkipLine = .false.
            mystuff%ifSkipSourceSpecies(:,:) = .false.
            mystuff%ifSkipCell(nmy_dispersion_mpi + 1) = .false.
         elseif (iTmp == dirichlet_boundary_type) then
            where (pBBuf%iBoundarySpecies(:, right_boundary) /= int_missing) &
                                                      & mystuff%ifSkipSourceSpecies(:,1) = .false.
            ifSkipLine = .false.
            mystuff%ifSkipCell(nmy_dispersion_mpi + 1) = .false.
         endif
      endif  ! global grid

    end subroutine get_linestuff_from_mm
                                 
  !**************************************************************************************
      

  subroutine encode_MPI_moving_masses(passengers, moment_line, iSpeciesIndex, &           ! masses, x/y-moment
                                    & iSpecies, iSrc, iLev, ix, iy, nPass, &      ! what coordinates ?
                                    & MPI_list_left_or_bottom, n_lb, MPI_list_right_or_top, n_rt) ! which domains?
    !
    ! This creature collects the stuff from single x- or y-line to ready-to-MPI-send lists
    !
    implicit none
    
    ! Imported parameters
    real, dimension(0:,0:), intent(in), target :: passengers       ! masses
    real, dimension(1:), intent(in), target :: moment_line         ! x/y-moment
    integer, dimension(:), pointer :: iSpeciesIndex        ! for "true" passengers
    integer, intent(in) :: iSpecies, iSrc, iLev, nPass     ! what coordinates ?
    integer, intent(in), target :: ix, iy                  ! what coordinates ?
    real, dimension(:), intent(out) :: MPI_list_left_or_bottom, MPI_list_right_or_top  ! which domains?      
    integer, intent(inout) :: n_lb, n_rt

    ! Local variables
    integer :: iPass, iLocal, iWing
    
    integer :: nx_dispersion_mpi, ny_dispersion_mpi
    
    nx_dispersion_mpi = nx_dispersion + wing_depth_e + wing_depth_w
    ny_dispersion_mpi = ny_dispersion + wing_depth_s + wing_depth_n
    
    ! Note that either ix or iy in the input is int_missing. Choose axis accordingly
    ! For real passengers, moments do not exist. Use real_missing.
    !
    ! At the left/bottom end we go from zero till iShiftLocal
    ! Note subtraction of the shift, i.e. all indices will be negative. Decoder will add the size 
    ! of the receiving domain.
    ! Note that we have zero-th line cell, which is out-of-domain in principle. But here it is a valid
    ! element - a safeguard against too high Courant.
    ! MPI_list_left_or_bottom(iItem,1:9) = (/ix-shift, iy-shift, iLev, iSrc, iSpecies, m, Px, Py, Pz/)
    !
    if(ix == int_missing)then
      !
      ! Move along x axis: northern or southern border
      !
      !n_lb = 1
      do iWing = 1, wing_depth_w
        if(passengers(0,iWing) == 0.0) cycle
        MPI_list_left_or_bottom(n_lb+1:n_lb+9) = (/real(iWing), &
                                             & real(iy), &
                                             & real(iLev), real(iSrc), real(iSpecies), &
                                             & passengers(0,iWing), &
                                             & moment_line(iWing), passengers(nPass+1,iWing), &
                                             & passengers(nPass+2,iWing)/)
        n_lb = n_lb + 9
        do iPass = 1, nPass
          MPI_list_left_or_bottom(n_lb+1:n_lb+9) = (/real(iWing), &
                                               & real(iy), &
                                               & real(iLev), real(iSrc), real(iSpeciesIndex(iPass)), &
                                               & passengers(iPass,iWing), &
                                               & real_missing, real_missing, real_missing/)
          n_lb = n_lb + 9
        end do  ! passengers
      end do  ! iWing 0:wing_depth_w
      !MPI_list_left_or_bottom(n_lb+1:size(MPI_list_left_or_bottom)) = 0.0
      !n_lb = n_lb - 1
      !
      ! At the right/top end we go from n?_disperion+shift+1 till n?_dispersion_mpi+1 - same logic as above
      ! Note that this time we subtract 
      !
      !n_rt = 1
       do iWing = 1, wing_depth_e
          iLocal = nx_dispersion + wing_depth_w + iWing  !index in input  lines
          if(passengers(0,iLocal) == 0.0) cycle
          MPI_list_right_or_top(n_rt+1:n_rt+9) = (/real(iWing), &
                                             & real(iy), &
                                             & real(iLev), real(iSrc), real(iSpecies), &
                                             & passengers(0,iLocal), &
                                             & moment_line(iLocal), &
                                             & passengers(nPass+1,iLocal),passengers(nPass+2,iLocal)/)
          n_rt = n_rt + 9
          do iPass = 1, nPass
            MPI_list_right_or_top(n_rt+1:n_rt+9) = (/real(iWing), &
                                               & real(iy), &
                                               & real(iLev), real(iSrc), real(iSpeciesIndex(iPass)), &
                                               & passengers(iPass,iLocal), &
                                               & real_missing, real_missing, real_missing/)
            n_rt = n_rt + 9
          end do  ! passengers
        enddo
      !MPI_list_right_or_top(n_rt+1:size(MPI_list_right_or_top)) = 0.0
      !n_rt = n_rt - 1
    else
      !
      ! Move along y dimension: western or eastern border
      !
      !n_lb = 1
      do iWing = 1, wing_depth_s
        if(passengers(0,iWing) == 0.0)cycle
        MPI_list_left_or_bottom(n_lb+1:n_lb+9) = (/real(ix), &
                                             & real(iWing), &
                                             & real(iLev), real(iSrc), real(iSpecies), &
                                             & passengers(0,iWing), &
                                             & passengers(nPass+1,iWing),moment_line(iWing), &
                                             & passengers(nPass+2,iWing)/)
        n_lb = n_lb + 9
        do iPass = 1, nPass
          MPI_list_left_or_bottom(n_lb+1:n_lb+9) = (/real(ix), &
                                              & real(iWing), &
                                              & real(iLev), real(iSrc), real(iSpeciesIndex(iPass)), &
                                              & passengers(iPass,iWing), &
                                              & real_missing, real_missing, real_missing/)
          n_lb = n_lb + 9
        end do  ! passengers
      end do  ! iWing 0:wing_depth_s
      !MPI_list_left_or_bottom(n_lb+1:size(MPI_list_left_or_bottom)) = 0.0
      !n_lb = n_lb - 1
      !
      ! At the right/top end we go from n?_disperion+shift+1 till n?_dispersion_mpi+1 - same logic as above
      ! Note that this time we subtract 
      !
      !n_rt = 1
      do iWing = 1,wing_depth_n
          iLocal = ny_dispersion + wing_depth_s + iWing !index in input  lines
          if(passengers(0,iLocal) == 0.0)cycle
          MPI_list_right_or_top(n_rt+1:n_rt+9) = (/real(ix), &
                                            & real(iWing), &
                                            & real(iLev), real(iSrc), real(iSpecies), &
                                            & passengers(0,iLocal), &
                                            & passengers(nPass+1,iLocal),moment_line(iLocal),&
                                            & passengers(nPass+2,iLocal)/)
          n_rt = n_rt + 9
          do iPass = 1, nPass
            MPI_list_right_or_top(n_rt+1:n_rt+9) = (/real(ix), &
                                              & real(iWing), &
                                              & real(iLev), real(iSrc), real(iSpeciesIndex(iPass)), &
                                              & passengers(iPass,iLocal), &
                                              & real_missing, real_missing, real_missing/)
            n_rt = n_rt + 9
          end do  ! passengers
        enddo
      !MPI_list_right_or_top(n_rt+1:size(MPI_list_right_or_top)) = 0.0
      !n_rt = n_rt - 1
    endif  ! x or y border?

  end subroutine encode_MPI_moving_masses
 
  
  !************************************************************************************
  subroutine decode_MPI_recv_buffer(MPI_list, nRecv, ixoff, iyoff,  iThread, nThreads, &
                                    & pDispFlds, moment_x, moment_y, moment_z, ifMoments)
    !
    ! This creature sends MPI_received list to  mass maps
    !
    implicit none
    
    ! Imported parameters
    real, dimension(:), intent(in) :: MPI_list      
    integer, intent(in) :: nRecv, iThread, nThreads
    integer, intent(in) :: ixoff, iyoff !offset to match list indices to local domain:
                                        !0 for left or second dimension, n_dispersion-wing_depth for right
    logical, intent(in) :: ifMoments
    type(Tmass_map), intent(inout) :: pDispFlds, moment_x, moment_y, moment_z 

    ! Local variables
    integer :: iItem, nxLocal, nyLocal, ix, iy, iLev, iSrc, iSpecies
    !
    ! MPI_list_left_or_bottom(iItem,1:9) = (/ix-shift, iy-shift, iLev, iSrc, iSpecies, m, Px, Py, Pz/)
    !                                            1         2      3      4        5    6   7   8   9
    !
    if (.not. ifMoments) then
         call set_error("decode_MPI_moving_masses deals exclusively with moments", "decode_MPI_moving_masses")
    endif

    do iItem = 1+(iThread*9), nRecv, 9*nThreads !! Simple way to share the job
        
        ix = nint(MPI_list(iItem)) + ixoff
        if(ix > nx_dispersion) call set_error("Gotcha iX","decode_MPI_recv_buffer")
        iy = nint(MPI_list(iItem+1)) + iyoff
        if(iy > ny_dispersion) call set_error("Gotcha iY","decode_MPI_recv_buffer")
        iLev = nint(MPI_list(iItem+2))
        iSrc = nint(MPI_list(iItem+3))
        iSpecies = nint(MPI_list(iItem+4))
        pDispFlds%arM(ispecies, iSrc, iLev, ix, iy) = pDispFlds%arM(ispecies, iSrc, iLev, ix, iy) + &
                                                    & MPI_list(iItem+5)
        if( MPI_list(iItem+6) /= real_missing)then   !Species has moment
          moment_x%arM(ispecies, iSrc, iLev, ix, iy) = moment_x%arM(ispecies, iSrc, iLev, ix, iy) + &
                                                     & MPI_list(iItem+6)
          moment_y%arM(ispecies, iSrc, iLev, ix, iy) = moment_y%arM(ispecies, iSrc, iLev, ix, iy) + &
                                                     & MPI_list(iItem+7)
          moment_z%arM(ispecies, iSrc, iLev, ix, iy) = moment_z%arM(ispecies, iSrc, iLev, ix, iy) + &
                                                     & MPI_list(iItem+8)
        endif
    end do  ! number of pieces
    
  end subroutine decode_MPI_recv_buffer

  subroutine decode_MPI_moving_masses(MPI_list_from_right_or_top, n_rt, MPI_list_from_left_or_bottom, n_lb, &
                                    & ifXAxis, nThreads, &
                                    & pDispFlds, moment_x, moment_y, moment_z, ifMoments)
    !
    ! This creature collects the stuff from single x- or y-line to ready-to-MPI-send lists
    ! Note that we
    !
    implicit none
    
    ! Imported parameters
    real, dimension(0:,:), intent(in) :: MPI_list_from_right_or_top, MPI_list_from_left_or_bottom !which domains?      
    integer, intent(in) :: n_lb, n_rt, nThreads
    logical, intent(in) :: ifXAxis
    logical, intent(in) :: ifMoments

    type(Tmass_map), intent(inout) :: pDispFlds, moment_x, moment_y, moment_z 

    ! Local variables
    integer :: iItem, nxLocal, nyLocal, ix, iy, iLev, iSrc, iSpecies, iThread
    !
    ! At the left/bottom end we go from zero till iShiftLocal
    ! Note subtraction of the shift, i.e. all indices will be negative. Decoder will add the size 
    ! of the receiving domain.
    ! Note that we have zero-th line cell, which is out-of-domain in principle. But here it is a valid
    ! element - a safeguard against too high Courant.
    ! MPI_list_left_or_bottom(iItem,1:9) = (/ix-shift, iy-shift, iLev, iSrc, iSpecies, m, Px, Py, Pz/)
    !                                            1         2      3      4        5    6   7   8   9
    !
    if (.not. ifMoments) then
         call set_error("decode_MPI_moving_masses deals exclusively with moments", "decode_MPI_moving_masses")
    endif

    if(ifXAxis)then
      nxLocal = nx_dispersion - wing_depth_e !!Offset of wings index in local mass map
      nyLocal = 0
    else
      nxLocal = 0
      nyLocal = ny_dispersion - wing_depth_n
    endif
    !
    ! Decode stuff came from left or bottom. Indexing is OK.
    !
    do iThread = 0, nThreads-1
      do iItem = 1, n_lb, 9
        ix = nint(MPI_list_from_left_or_bottom(iThread,iItem))
        iy = nint(MPI_list_from_left_or_bottom(iThread,iItem+1))
        iLev = nint(MPI_list_from_left_or_bottom(iThread,iItem+2))
        iSrc = nint(MPI_list_from_left_or_bottom(iThread,iItem+3))
        iSpecies = nint(MPI_list_from_left_or_bottom(iThread,iItem+4))
        pDispFlds%arM(ispecies, iSrc, iLev, ix, iy) = pDispFlds%arM(ispecies, iSrc, iLev, ix, iy) + &
                                                    & MPI_list_from_left_or_bottom(iThread,iItem+5)
        if( MPI_list_from_left_or_bottom(iThread,iItem+6) /= real_missing)then
          moment_x%arM(ispecies, iSrc, iLev, ix, iy) = moment_x%arM(ispecies, iSrc, iLev, ix, iy) + &
                                                     & MPI_list_from_left_or_bottom(iThread,iItem+6)
          moment_y%arM(ispecies, iSrc, iLev, ix, iy) = moment_y%arM(ispecies, iSrc, iLev, ix, iy) + &
                                                     & MPI_list_from_left_or_bottom(iThread,iItem+7)
          moment_z%arM(ispecies, iSrc, iLev, ix, iy) = moment_z%arM(ispecies, iSrc, iLev, ix, iy) + &
                                                     & MPI_list_from_left_or_bottom(iThread,iItem+8)
        endif
      end do  ! number of pieces
    end do  ! threads from left or bottom
    !
    ! Decode stuff came from right or top. Indices are negative.
    !
    do iThread = 0, nThreads-1
      do iItem = 1, n_rt, 9
        ix = nint(MPI_list_from_right_or_top(iThread,iItem)) + nxLocal
        iy = nint(MPI_list_from_right_or_top(iThread,iItem+1)) + nyLocal
        iLev = nint(MPI_list_from_right_or_top(iThread,iItem+2))
        iSrc = nint(MPI_list_from_right_or_top(iThread,iItem+3))
        iSpecies = nint(MPI_list_from_right_or_top(iThread,iItem+4))
        pDispFlds%arM(ispecies, iSrc, iLev, ix, iy) = pDispFlds%arM(ispecies, iSrc, iLev, ix, iy) + &
                                                    & MPI_list_from_right_or_top(iThread,iItem+5)
        if( MPI_list_from_left_or_bottom(iThread,iItem+6) /= real_missing)then
          moment_x%arM(ispecies, iSrc, iLev, ix, iy) = moment_x%arM(ispecies, iSrc, iLev, ix, iy) + &
                                                     & MPI_list_from_right_or_top(iThread,iItem+6)
          moment_y%arM(ispecies, iSrc, iLev, ix, iy) = moment_y%arM(ispecies, iSrc, iLev, ix, iy) + &
                                                     & MPI_list_from_right_or_top(iThread,iItem+7)
          moment_z%arM(ispecies, iSrc, iLev, ix, iy) = moment_z%arM(ispecies, iSrc, iLev, ix, iy) + &
                                                     & MPI_list_from_right_or_top(iThread,iItem+8)
        endif
      end do  ! number of pieces
    end do  ! threads from left or bottom
    
  end subroutine decode_MPI_moving_masses

  !**********************************************************************************

  subroutine count_incoming_cells(lb_advected, rb_advected, rb, wing_depth_l,  &
                           & wing_depth_r, n_dispersion, NinL, nInR)

         implicit none
         ! Check how many their cells are advected into our domain
         real(r8k), intent(in) :: lb_advected,  rb_advected !domain boundaries
         real(r8k), dimension(0:), intent(in) :: rb !cellboundaries of Whole line
         integer, intent(in) :: wing_depth_l, wing_depth_r, n_dispersion
         integer, intent(out) :: NinL, NinR !Counters to increment
         integer :: iwing
         
         nInL = 0 
         if (wing_depth_l > 0) then
           do iwing = 0, wing_depth_l
                if( lb_advected < rb(wing_depth_l+iwing)) then
                  nInL = iwing
                  exit
               endif
           enddo
           if (iwing > wing_depth_l) call set_error("Incoming from Left beyond the wing", &
                                   &                "count_incoming_cells")
         endif

         nInR = 0 
         if (wing_depth_r > 0) then
           do iwing = 0, wing_depth_r
               if(rb_advected > rb(wing_depth_l + n_dispersion - iwing)) then
                  nInR = iwing
                  exit
               endif
            enddo
            if (iwing > wing_depth_r) call set_error("Incoming from Right beyond the wing", &
                                   &                "count_incoming_cells")
         endif
   end subroutine count_incoming_cells
  
  
  !************************************************************************************
  !************************************************************************************
  !************************************************************************************
  !
  !  Vertical advection routines start
  !
  !************************************************************************************
  !************************************************************************************
  !************************************************************************************



  !************************************************************************************

  subroutine adv_diffusion_vertical_v5(pDispFlds, &
                                           & moment_x, moment_y, moment_z, &
                                           & pAerosolFlds, & ! Not needed here
                                           & pHorizInterpStruct, pVertInterpStruct, &
                                           & ifHorizInterp, ifVertInterp, &
                                           & pMetBuf, pDispBuf, &
                                           & seconds, weight_past, &
                                           & garbage, &
                                           & bottom_mass, top_mass, &
                                           & pDryDep, &
                                           & pCnc2m, &
                                           & pBBuf, &
                                           & chem_rules, &
                                           & have_negatives, wdr, ifAllMoments, ifTalk)
    ! Works with massfluxes
    ! Treats settling as advection except for the first level
    ! Can work with any vertical
    !
    implicit none
    ! Imported parameters
    type(Tmass_map), intent(inout) :: pDispFlds, pDryDep, pAerosolFlds, &
                                    &  moment_x, moment_y, moment_z, pCnc2m  ! Might be unassociated
    type(THorizInterpStruct), pointer :: pHorizInterpStruct
    type(TVertInterpStruct), pointer :: pVertInterpStruct
    logical, intent(in) :: ifHorizInterp, ifVertInterp
    type(Tfield_buffer), intent(in) :: pMetBuf, pDispBuf
    real, intent(in) :: seconds, weight_past
    real, dimension(:,:), intent(inout) :: garbage
    real, dimension(:,:,:), intent(inout) :: bottom_mass, top_mass
    type(Tchem_rules), intent(in) :: chem_rules
    logical, intent(in) :: have_negatives, ifTalk
    logical, intent(inout) :: ifAllMoments
    TYPE(silja_wdr), INTENT(in) :: wdr
    type(TboundaryBuffer), intent(in), pointer :: pBBuf

    ! Local variables
    integer :: isrc, ilev, iLevMet, iPass,  ix, iy, ispecies, ixMeteo, iyMeteo, indexMeteo, indexDisp
    integer ::   nTimeSteps, nSrc, nspecies, nPass, nThreads
    integer ::   m_ind, w_ind, wadj_ind, r_down_met_ind,  temper_ind, &
         & rh_ind, ps_ind, zSize_ind, z_ind
    integer ::  iLevMinDiff, iLevMaxDiff, jTmp, iTmp  
    logical ::  ifDryDep, ifColparamsNotReady, ifBark, ifSkipColumn
    real :: fMass, fMinAdvectedMass, &
         &  fZeroMassThreshold,   seconds_abs, ftmp, ftmp1, weight_up
    real :: dh1, dh2,  cnc2mfrac1, cnc2mfrac2 
    real :: fTemper, fRelHum, XcmTmp, YcmTmp, ZcmTmp
    real(r8k) ::  ps, p2m
    integer :: timeSign
    real :: fVd, iVd2m ! deposition velocity for species according to KS2011
    real, dimension(:), pointer :: tmp_garbage, arMinAdvMass
    real, dimension(:), pointer :: fPtr
    integer, dimension(:), pointer :: iarrPtr
    type(Taerosol_rules), pointer :: rulesAerosol
    type(Tdeposition_rules), pointer :: deprules
    real ::    weightUp 
    integer :: leveltype, n_time_steps, iThread
    integer :: levMinadvection, levMaxadvection
    integer :: bark_count, total_bark
    logical, dimension(2) :: ifPoles
    integer :: spthread = int_missing ! Thread made souh pole advection, saved for reporting 

    logical :: ifSettlingNeeded, ifTunedRs
    !    type(TwetParticle) :: WetParticle

    type(T_Galp_v5_thread_stuff), pointer :: mystuff, zstuff
    type(T_Galp_v5_vert_shared_stuff), pointer :: stuff

    real, dimension(:), pointer :: pXCellSize, pYCellSize 
    
    ifTunedRs = .false.
#ifdef SKIP_VERTICAL         
    call msg("Skipping vertical !!!!")
    return
#endif
    ! Needed for boundaries
    pXCellSize => fu_grid_data(dispersion_cell_x_size_fld)
    pYCellSize => fu_grid_data(dispersion_cell_y_size_fld)

    total_bark = 0

    if (do_vert_diff)   do_vert_diff = (fu_kz_param(wdr) /= zero_kz)
    
    nspecies = pDispFlds%nSpecies
    nSrc = pDispFlds%nSrc
    leveltype = fu_leveltype(dispersion_vertical)

    ! Meteo indices
    iArrPtr => pMetBuf%buffer_quantities
    r_down_met_ind = fu_index(iArrPtr, R_down_meteo_flag)
    temper_ind = fu_index(iArrptr, temperature_flag)
    rh_ind = fu_index(iArrptr, relative_humidity_flag)
    ps_ind = fu_index(iArrptr, surface_pressure_flag)
    z_ind = fu_index(iArrptr, height_flag)

    !Dispersion Indices
    iArrPtr => pDispBuf%buffer_quantities
    m_ind = fu_index(iArrPtr, disp_cell_airmass_flag)
!    v_ind = fu_index(iArrptr, disp_flux_cellnorth_flag)
!    vadj_ind = fu_index(iArrptr, disp_flux_celln_rt_flag)
    w_ind = fu_index(iArrptr, disp_flux_celltop_flag)
    wadj_ind = fu_index(iArrptr, disp_flux_cellt_rt_flag)
    zSize_ind = fu_index(iArrPtr, cell_size_z_flag)


    if (any((/w_ind, wadj_ind, m_ind, ps_ind, r_down_met_ind/) == int_missing)) then
      call msg('Dispersion: w_ind, wadj_ind, m_ind', (/w_ind,  wadj_ind, m_ind/))
      call msg('Meteo: ps_ind, r_down_met_ind', ps_ind, r_down_met_ind)
      call set_error('Not all meteo quantities found', 'adv_diffusion_vertical_v5')
      return
    end if

    depRules => fu_deposition_rules(chem_rules)

    if(.not. fu_ifMeteoGrd(pVertInterpStruct) .and. &
         & .not.(fu_grid(pVertInterpStruct) == dispersion_grid))then
      call set_error('Neither meteo nor dispersion grid in vertical interpolation strcucture', &
                   & 'adv_diffusion_vertical_v5')
      return
    endif

    
    call verify_deposition(pDispFlds%species, nspecies, deprules)
    if (error) return

    fZeroMassThreshold = 0.0 !fu_zero_mass_threshold(chem_rules)
    arMinAdvMass => fu_low_mass_threshold(chem_rules)

    seconds_abs = abs(seconds)
    timeSign = sign(1,int(seconds))

    stuff => EulerStuff%Galp5zShared

    ! Do we need settling velocity?
    do iSpecies= 1,nSpecies
       stuff%ifSettles(iSpecies) = ( fu_massmean_d(pDispFlds%species(iSpecies)%mode) > 1e-9)
    enddo

    


    ifSettlingNeeded = any(stuff%ifSettles)
#ifdef DEBUG_V5
    if (ifSettlingNeeded) then
       call msg("Settling active")
    else
       call msg("No settling calculated")
    endif
#endif
    if ((.not. stuff%ifMoldiffInitialised) .and. associated(stuff%xminus_moldiff)) then
      ! Initialize if needed
      call init_molec_diff(pDispFlds%species, stuff%a_half_disp(1:nz_dispersion+1), stuff%b_half_disp(1:nz_dispersion+1), &
            & stuff%xplus_moldiff, stuff%xminus_moldiff, & 
            & nSpecies, nz_dispersion, seconds_abs)
       stuff%ifMoldiffInitialised = .true.
       do_molec_diff = .true. 
    endif

    ifPoles(northern_boundary) = (pBBuf%iBoundaryType(northern_boundary) == polar_boundary_type)
    ifPoles(southern_boundary) = (pBBuf%iBoundaryType(southern_boundary) == polar_boundary_type)

    !if(error)return
    
    !$OMP PARALLEL DEFAULT(NONE) &
    !$OMP & PRIVATE(ilev, ix, iy,  iLevMinDiff, iLevMaxDiff, iThread, mystuff, &
    !$OMP         &  ispecies,  ixMeteo, iyMeteo, indexMeteo,  iPass, nPass, &
    !$OMP         & jTmp, iTmp, ifSkipColumn, &
    !$OMP         &  ifDryDep, fTemper, fRelHum, &
    !$OMP         & fMass,   dh1, dh2,&
    !$OMP         &  tmp_garbage, ftmp, indexDisp, ftmp1, &
    !$OMP         & fMinAdvectedMass,  fVd, iVd2m, p2m, cnc2mfrac1, cnc2mfrac2,  weight_up, ifColparamsNotReady, &
    !$OMP         & weightUp,  ps, isrc,  N_time_steps, XcmTmp, YcmTmp, ZcmTmp, &
    !$OMP         & levMinadvection, levMaxadvection, ifBark, bark_count) &
    !$OMP & SHARED( stuff, z_ind, r_down_met_ind, m_ind, w_ind, wadj_ind,  &
    !$OMP        & ps_ind,  temper_ind, rh_ind, zSize_ind, ifSettlingNeeded, do_vert_diff, &
    !$OMP        & nspecies,  leveltype, nSrc, top_mass, bottom_mass,  spthread, nthreads,&
    !$OMP        & seconds_abs, timeSign,  pDispFlds, moment_x, moment_y, moment_z, pDryDep, &
    !$OMP        & pHorizInterpStruct, pVertInterpStruct, ifHorizInterp, ifVertInterp, pMetBuf,&
    !$OMP        & pDispBuf, pBBuf, seconds, weight_past, garbage, chem_rules, nx_dispersion, &
    !$OMP        & ny_dispersion, nx_meteo,  fs_meteo, nz_meteo, disp_layer_top_m, a_met, b_met, &
    !$OMP        & arMinAdvMass, depRules, nz_dispersion, have_negatives, wdr,  ifAllMoments, &
    !$OMP        & pXCellSize, pYCellSize,  error, pcnc2m, EulerStuff,  total_bark, &
    !$OMP        & ifPoles,  ifTalk, ifTunedRs, do_molec_diff, diffuse_cm_vert)


    iThread = 0

    !$ iThread = OMP_GET_THREAD_NUM()

    !$OMP MASTER
    nThreads = 1
    !$ nThreads = omp_get_num_threads()
#ifdef DEBUG
    call msg('DEBUG adv_diffusion_vertical_v5')
    if (do_molec_diff) call msg("With molecular diffusion")
    if (diffuse_cm_vert) call msg("With sub-grid diffusion")
#else
      call msg('adv_diffusion_vertical_v5')
#endif
    !$OMP END MASTER
    !$OMP BARRIER

    mystuff => EulerStuff%Galp5(iThread)
    mystuff%tmp_garbage = 0.
    mystuff%bottom_mass = 0.
    mystuff%top_mass = 0.

    bark_count = 0


    ! These will be overwritten by boundaries should there be any
    mystuff%passengers(:,              0,:,:) = 0.
    mystuff%passengers(:,nz_dispersion+1,:,:) = 0.
    mystuff%vSettling(:,:) = 0.  !Will be reset, should it be necessary

    mystuff%cMax1d(1:nz_dispersion) = 0.
    mystuff%Cmean1d(1:nz_dispersion) = 0.
    mystuff%Ccount(1:nz_dispersion) = 0
    
    !$OMP DO  collapse(2)
    DO iy=1, ny_dispersion
      ix_DO: DO ix=1, nx_dispersion
        if (error) cycle

!report_diff = .false. ! (ix == 100 .and. iy == 100)
!call msg('ix', ix)
!call msg('iy', iy)
if(ifTalk)call msg('Starting cycle:',(/ix,iy,iLev,iSpecies/))

        indexDisp = ix + (iy-1)*nx_dispersion !!FIXME should it handle mpi somehow?

        !
        ! Check if column valid and fill local arrays
        !
        call get_linestuff_from_mm(pDispFlds, moment_z, moment_x, moment_y, ix, iy, &
                                 & nz_dispersion, 0, nz_dispersion, &
                                 & 3, .false., bottom_boundary, top_boundary, &
                                 & mystuff, nSrc, nspecies, have_negatives, arMinAdvMass, pBBuf, &
                                 & ifAllMoments, ifSkipColumn)
        !
        ! ifColumnValid is true if ANY of the species has any meaning. Note that sources is the 
        ! second index in ifSkipSourceSpecies
        !
        pDispFlds%ifColumnValid(:,ix,iy) = .not. all(mystuff%ifSkipSourceSpecies(:,:), 1) 

        levMinAdvection = 1
        levMaxAdvection = nz_dispersion
       
        if(ifSkipColumn) cycle

if(ifTalk)then
  call msg('Column valid:',(/ix,iy,iLev,iSpecies/))
  if(ix == 43 .and. iy == 72)then
    call msg('*** here***')
  endif
endif
        
!        !Fixme Hack
!        call msg("Forcing column")
!        mystuff%passengers(0,:,1, 1) = 0.
!        mystuff%cm_line(:,1, 1) = 0.
!        mystuff%passengers(0,levMaxAdvection,1,1) = 1e10
!        levMinAdvection = levMaxAdvection
!        mystuff%passengers(:,:,1, 1) = 0.
!
        !
        ! Get meteo index and settling velocities
        !
          
        indexMeteo = fu_grid_index(nx_meteo, ix, iy, pHorizInterpStruct)
        if(indexMeteo < 1 .or. indexMeteo > fs_meteo)then
            call msg('ncoefs', fu_ncoefs(pHorizInterpStruct))
            call msg('IndexMeteo:',indexMeteo)
            call set_error('Bad meteo index','adv_diffusion_vertical_v5')
            cycle
        endif
      
        !
        ! Set basic meteo for the column 
        ! TODO:  fu_get_value costs a lot. 
        ! Better be replaced with dispersion leveltop quantities 
        ! 

        ! pressure
        !   Make mystuff%PAbove and mystuff%VertMotionAbove

        ps =  fu_get_value(pMetBuf%p2d(ps_ind), nx_meteo, ix, iy, &
                     & weight_past, pHorizInterpStruct, ifHorizInterp)

        ! Meteo is always in hybid-type coordinates
        mystuff%p_met(0:nz_meteo+1) =  stuff%a_met(0:nz_meteo+1) + stuff%b_met(0:nz_meteo+1)*ps



        !   This part is different for eta and z-coordinates
!        mystuff%PAbove(-1) = 2*ps !Does not affect advection, matters for diffusion
        if(leveltype == layer_btw_2_height)then !z-coordinates

          !Getting pressure from meteo fields does not work 
          ! if model levels are outside of meteo level range
          !  !Sea Salt failed with ECHAM input
          ! Then can not rely on Vertical interpolation anymore
          ! Interpolating meteo pressures using meteo heights field
          mystuff%PAbove(0)  = ps

          !Inter/extrapolate pressure
          mystuff%z_met(0)  = 0. 
          call column_from_buffer(pMetBuf, z_ind, ix, iy, nx_meteo, &
                                  & mystuff%z_met(1:nz_meteo), &
                                  & pHorizInterpStruct, stuff%pVertInterpTop, &
                                  & ifHorizInterp, .false., weight_past)
          iTmp=1
          do iLev=1,nz_dispersion
               fTmp=stuff%disp_layer_top_m(iLev)
               do iTmp = iTmp, nz_meteo
                  if (mystuff%z_met(iTmp) > fTmp) exit
               enddo
               fTmp1 = (fTmp-mystuff%z_met(iTmp))/(mystuff%z_met(iTmp-1)-mystuff%z_met(iTmp))
               mystuff%PAbove(iLev) = mystuff%p_met(iTmp-1)*fTmp1 + mystuff%p_met(iTmp)*(1.-fTmp1)
          enddo

          if (any(mystuff%PAbove(0:nz_dispersion-1) <= mystuff%PAbove(1:nz_dispersion)  )) then
                  !$OMP critical(barkv4)
                        
                        call msg("mystuff%PAbove(0:nz_dispersion)",mystuff%PAbove(0:nz_dispersion))
                        call msg('ix, iy',ix, iy)
                        call set_error("non_monotonous pressure at level top...", 'adv_diffusion_vertical_v4')
                 !$OMP END critical(barkv4)
                 continue
          endif
        else  !pressure-coordinates
          ! a and b refer to level bottom
          ! mystuff%PAbove to level top
          mystuff%PAbove(0:nz_dispersion) = stuff%a_half_disp(1:nz_dispersion+1) + &
                                          & stuff%b_half_disp(1:nz_dispersion+1) * ps
        endif  !leveltype
!        mystuff%PAbove(nz_dispersion+1) = mystuff%PAbove(nz_dispersion) - 100  !Pa, virtually anything

        ! Density, vertical motion and settling
        mystuff%wind_right(0) = 0.0


        do ilev = 1, nz_dispersion

          mystuff%wind_right(ilev) = pDispBuf%p4d(w_ind)%past%p2d(ilev)%ptr(indexDisp)*weight_past + &
                     & pDispBuf%p4d(w_ind)%future%p2d(ilev)%ptr(indexDisp) * (1. - weight_past) + &
                     & pDispBuf%p4d(wadj_ind)%present%p2d(ilev)%ptr(indexDisp) 

          mystuff%cellmass(iLev) =  pDispBuf%p4d(m_ind)%past%p2d(ilev)%ptr(indexDisp) *weight_past + &
                           &   pDispBuf%p4d(m_ind)%future%p2d(ilev)%ptr(indexDisp) * (1. - weight_past)

          !mystuff%rho_above(ilev) = mystuff%PAbove(iLev) / (gas_constant_dryair *fTemper)
          ! Not exactly accurate, but should be cheaper than interpolating from meteo
          mystuff%rho_above(ilev) = mystuff%cellmass(iLev) /  &
             & ( pXCellSize(indexDisp) * pYCellSize(indexDisp)  &
             &    * ( pDispBuf%p4d(zSize_ind)%past%p2d(ilev)%ptr(indexDisp) * weight_past + &
             &       pDispBuf%p4d(zSize_ind)%future%p2d(ilev)%ptr(indexDisp) * (1. -weight_past)))
        end do
        mystuff%rho_above(0) = mystuff%rho_above(1) !Used in the diffusion
       
        ! Some Stats (only the current MPI domain)
        !cheap and dirty calculation of courant 
        mystuff%scratch(1:nz_dispersion) = &
                 & abs(seconds*mystuff%wind_right(1:nz_dispersion) / mystuff%cellmass(1:nz_dispersion))
        do ilev = 1, nz_dispersion
             mystuff%cMax1d(iLev) = max(mystuff%scratch(iLev), mystuff%cMax1d(iLev))
             mystuff%Cmean1d(iLev) = mystuff%Cmean1d(iLev) + mystuff%scratch(iLev)
             mystuff%Ccount(iLev) = mystuff%Ccount(iLev) + 1
        enddo
        ! Get settling velocities (if we have particles)
        if (ifSettlingNeeded) then
           do ilev = 1, nz_dispersion
             fTemper = fu_get_value(pMetBuf%p4d(temper_ind), &
                                   & nx_meteo, ix, iy, iLev, weight_past, &
                                   & pHorizInterpStruct, stuff%pVertInterpTop, &
                                   & ifHorizInterp, .true.) !Always interpoate
             fRelHum = fu_get_value(pMetBuf%p4d(rh_ind), &
                                   & nx_meteo, ix, iy, iLev, weight_past, &
                                   & pHorizInterpStruct, stuff%pVertInterpTop, &
                                   & ifHorizInterp, .true.) !Always interpoate

             call get_settling_motion_species(pDispFlds%species, nSpecies, &
                                            & fTemper, real(mystuff%PAbove(iLev)), &
                                            & mystuff%rho_above(ilev), fRelHum, &
                                            & depRules, mystuff%vSettling(:,ilev), &
                                            & .true.) !Get settling in Pa/s
!      call msg("get_settling_motion_species returned in Pa/s for lev:"+fu_str(iLev),  mystuff%vSettling(:,ilev))
           end do

           mystuff%vSettling(:,0) = 0. !  Settlins at surface is treated with Vd

           mystuff%vSettling(:,:) = - mystuff%vSettling(:,:) * &
               & pXCellSize(indexDisp) * pYCellSize(indexDisp) / g !to kg/s positive upwards
              
         endif

        !Metric thicknesses of first two layers needed for cnc2m diagnostics
        if(leveltype == layer_btw_2_height)then !z-coordinates
          dh1 = stuff%disp_layer_top_m(1) - stuff%disp_layer_top_m(0)
          dh2 = stuff%disp_layer_top_m(2) - stuff%disp_layer_top_m(1)
        else  !pressure-coordinates
          dh1 = 2*(mystuff%PAbove(0)-mystuff%PAbove(1)) / ( (mystuff%rho_above(0) + mystuff%rho_above(1)) * g) 
          dh2 = 2*(mystuff%PAbove(1)-mystuff%PAbove(2)) / ( (mystuff%rho_above(1) + mystuff%rho_above(2))* g) 
        endif  !leveltype

        !
        ! Get aerodynamic resistance for diffusion and deposition
!!!      ! Do we need any accuracy here ???
        if(do_vert_diff)then  
             do iLevMet=1,nz_meteo !get column in meteo levels
                 mystuff%r_down_met(iLevMet) = &
                   & pMetBuf%p4d(r_down_met_ind)%past%p2d(iLevMet)%ptr(indexMeteo) * weight_past + &
                   & pMetBuf%p4d(r_down_met_ind)%future%p2d(iLevMet)%ptr(indexMeteo) *  (1.-weight_past)
               if(mystuff%r_down_met(iLevMet) < 0)then
                 call msg('Negative meteo resistance at reading from meteobuffer:', &
                                                              & iLevMet, mystuff%r_down_met(iLevMet))
                 call msg('Meteo and dispersion indices:',(/indexMeteo, ix,iy/))
               endif
             enddo
          mystuff%r_down_met(nz_meteo+1) = mystuff%r_down_met(nz_meteo) !Whatever
          mystuff%r_down_met(0) = 0.


          if ( deprules%RsType == DryD_Rs_standard) then
            call get_rs(pDispFlds%species, nspecies, indexMeteo, weight_past, deprules, mystuff%R_Surf)
          else
            ! Get mmr of species for near-surface cell (sum over sources)
            call get_Rs_2013(pDispFlds%species, sum(mystuff%passengers(0,1,:,:),dim=2)/ mystuff%cellmass(1),  &
                           & nspecies, indexMeteo, weight_past, deprules, &
                           & mystuff%R_Surf, ifTunedRs)
          endif

          ! Can be done once for all species
          mystuff%g_over_deltaP(1:nz_dispersion) =  pXCellSize(indexDisp) * pYCellSize(indexDisp) & 
                                                   & / mystuff%cellmass(1:nz_dispersion)
          mystuff%g_over_deltaP(0) = 0

          ! Prepare sub-grid vertical diffusion
          if (diffuse_cm_vert) then
            call  prepareDiffCMvert(mystuff%cm_relax_diff, mystuff%rho_above, mystuff%PAbove, &
            &  mystuff%r_down_met, mystuff%p_met,  & !Meteo  R1m
            &  seconds_abs) 
          endif

        endif ! do_vert_diff


        ! Boundaries
        if (seconds > 0) then
         !FIXME Mass to inject depends on settling!
         ! Should be fixed at some point !! R
         call get_linestuff_from_boundaries(pdispFlds,moment_z, ix, iy, nz_dispersion, &
               & bottom_boundary, top_boundary, levMinadvection, levMaxadvection, &
               & mystuff%wind_right(0)*seconds,  -mystuff%wind_right(nz_dispersion)*seconds, &
               & mystuff%rho_above(0), & !airdens_left
               & mystuff%rho_above(nz_dispersion), &! airdens_right
               & mystuff, nSrc, nSpecies, arMinAdvMass, pBBuf, & 
               & mystuff%bottom_mass(:, :, incoming), mystuff%top_mass(:, :, incoming))
         else
            mystuff%passengers(:,0,:,:) = 0.
            mystuff%passengers(:,nz_dispersion+1,:,:) = 0.
         endif

        !
        ! Real advection starts  
        !
        ifColparamsNotReady=.True.
        do ispecies = 1, nspecies
          if (error) cycle
          nPass = moment_z%passengers(iSpecies)%nPassengers   
          fMinAdvectedMass = arMinAdvMass(ispecies)
          if (ifColparamsNotReady .or. stuff%ifSettles(ispecies)) then
            !
            ! update advected cell boundaries        

            ! NB: all velosities here are in Pa/s -- settling -positive
            mystuff%windRightTmp(0:nz_dispersion) = &
                & timeSign*(mystuff%wind_right(0:nz_dispersion) &
                          & + mystuff%vSettling(ispecies,0:nz_dispersion))


            !HACK!! 
            !mystuff%VertMotionAboveSpc(:) = 0.      !Hack to zero vertical advection/settling
!           mystuff%VertMotionAboveSpc(:) = 1. !Pa/s     !Constant vertical motion 

!
            mystuff%ifSkipCell(0:nz_dispersion+1) = .false. !FIXME hack
            call advect_cellboundaries(mystuff%windRightTmp(0:nz_dispersion), &!Staggered
                     & mystuff%cellmass(1:nz_dispersion),  & ! non-staggered
                     &     mystuff%lineParams(0:nz_dispersion),  & !Staggered
                     & mystuff%lineParamsAdvected(0:nz_dispersion),       &!Staggered
                     &     mystuff%ifSkipCell(0:nz_dispersion+1), & !non-staggered
                     &     mystuff%scratch,   & !staggered
                     &     abs(seconds), .false., levMinadvection, levMaxadvection, nz_dispersion)


            ! Finite settling -- can't reuse this advecion
            ifColparamsNotReady  = stuff%ifSettles(ispecies)
            if(error) then
               call advect_cellboundaries(mystuff%windRightTmp(0:nz_dispersion), &!Staggered
                     & mystuff%cellmass(1:nz_dispersion),  & ! non-staggered
                     &     mystuff%lineParams(0:nz_dispersion),  & !Staggered
                     & mystuff%lineParamsAdvected(0:nz_dispersion),       &!Staggered
                     &     mystuff%ifSkipCell(0:nz_dispersion+1), & !non-staggered
                     &     mystuff%scratch,   & !staggered
                     &     abs(seconds), .false., levMinadvection, levMaxadvection, nz_dispersion)

               call msg("Trouble with Z")
               cycle
            endif

          endif ! Recycle colparams_advected


          do iSrc=1,nSrc
            
if(ifTalk)call msg('Column valid??',(/ix,iy,iLev,iSpecies, iSrc/))
            
            if(.not. pDispFlds%ifColumnValid(iSrc,ix,iy)) cycle

if(ifTalk)call msg('Prior to advect_mass:',(/ix,iy,iLev,iSpecies/))

            call  MASS_DISTRIBUTOR(mystuff%cm_line(1:,iSpecies,iSrc),  mystuff%moment_line_out(1:), &
                          & mystuff%passengers(:,0:,iSpecies,iSrc), mystuff%passengers_out(:,0:),&
                          & nPass+2,  & ! number of passengers + 2 moments
                          & mystuff%lineParams, mystuff%lineParamsAdvected,  mystuff%cellmass,& 
                          & .false., levMinadvection, levMaxadvection, nz_dispersion,&
                          & arMinAdvMass(ispecies), have_negatives, .true., 1.0) !Convert to CM
!             call msg( " MIn (0:nz+1)", mystuff%passengers(0,0:nz_dispersion+1,iSpecies, iSrc) )
!             call msg("  MOut  (0:nz+1)", mystuff%passengers_out(0,0:nz_dispersion+1))


#ifdef DEBUG
            fTmp = sum(mystuff%passengers(0,0:nz_dispersion+1,iSpecies, iSrc))
            fTmp1 = sum(abs(mystuff%passengers(0,0:nz_dispersion+1,iSpecies, iSrc)))


            if ((.not. all(abs(mystuff%passengers_out(0,0:nz_dispersion+1)) >= 0. )) .or. & !NaNs
               & ( (.not. have_negatives) .and. any(mystuff%passengers_out(0,0:nz_dispersion+1)<0. )) .or. &!Invalid neg
               & (abs(fTmp - sum(mystuff%passengers_out(0,0:nz_dispersion+1)))/ (fTmp1+fMinAdvectedMass)) > 1e-5 ) then
              
              call msg("################ Trouble  After Advection #######################", ix, iy)
              call msg("Column mass 3 before, after", fTmp,   sum(mystuff%passengers_out(0,0:nz_dispersion+1)) )
              call msg("mystuff%PAbove -1:nz+1", mystuff%PAbove(0:nz_dispersion))
              call msg("mystuff%lineParams   :", mystuff%lineParams(0:nz_dispersion))
              call msg("mystuff%lineParamsAdv:", mystuff%lineParamsAdvected(0:nz_dispersion))
              call msg( " MIn (0:nz+1)", mystuff%passengers(0,0:nz_dispersion+1,iSpecies, iSrc) )
              call msg("  MOut  (0:nz+1)", mystuff%passengers_out(0,0:nz_dispersion+1))
              call msg("CMIn (0:nz+1)",mystuff%cm_line(1:nz_dispersion,iSpecies,iSrc) )
              !call msg("CmOut (1:nz)", mystuff%cm_line_out(1:nz_dispersion) )
              call set_error("Trouble after advection1! End", "AdvVertV5") 
              cycle iX_do
            endif
#endif
              
            ! Collect outgoing mass and free zero cell for deposition
            ! Should be done after diffusion ...
            ! ...once it would be capable bringing stuff out of domain...
            if ( abs(mystuff%passengers_out(0,0)) > arMinAdvMass(ispecies)) then
                  call set_error("Gotcha bottom mass", "agergdefbhs")
               endif
            mystuff%bottom_mass(iSrc,iSpecies, outgoing) = mystuff%bottom_mass(nSrc,iSpecies, outgoing) + &
                                                  &  mystuff%passengers_out(0,0)
            mystuff%top_mass(iSrc,iSpecies, outgoing) =  mystuff%top_mass(iSrc,iSpecies, outgoing) + &
                                                  & mystuff%passengers_out(0,nz_dispersion+1) 
            
            do iPass = 1, nPass
               iTmp = moment_z%passengers(iSpecies)%iSpeciesGarbage(iPass)
               mystuff%bottom_mass(iSrc,iTmp, outgoing) =  mystuff%bottom_mass(iSrc,iTmp, outgoing)  + &
                      & mystuff%passengers_out(iPass,0)
               mystuff%passengers_out(iPass,0)  = 0.
               mystuff%top_mass(iSrc,iTmp, outgoing) =  mystuff%top_mass(iSrc,iTmp, outgoing)  + &
                      & mystuff%passengers_out(iPass, nz_dispersion+1)
               mystuff%passengers_out(iPass, nz_dispersion+1)  = 0.
            end do  ! passengers

            ! Free bottom and top for diffusion       
            ! If we do not kill momens here -- we are in trouble
            mystuff%passengers_out(:,0) = 0.
            mystuff%passengers_out(:,nz_dispersion+1) = 0.


#ifdef   DEBUG_V5
!####################################Temp check: FIXME To be removed
            if ( .not. have_negatives) then
               do  iLev = 1, nz_dispersion  ! Vertical levels
                 if (error) cycle
                 fmass = mystuff%passengers_out(0,iLev) 
                 if(abs(fMass) <= arMinAdvMass(iSpecies))then
                   mystuff%passengers_out(nPass+1,iLev) = 0.
                   mystuff%passengers_out(nPass+2,iLev) = 0.
                 else
                   !These guys should be moments already
                   do jTmp = 1,2 
                     XcmTmp = mystuff%passengers_out(nPass+jTmp,iLev) / fMass
                     if(.not. abs(XcmTmp) <= 0.5) then
                         if (jTmp == 1) then
                            call msg('StrangeTmp X position ', XcmTmp)
                         else
                            call msg('StrangeTmp Y position ', XcmTmp)
                         endif
                         call msg('Species:' + fu_str(pDispFlds%species(iSpecies)), iSpecies, arMinAdvMass(iSpecies))
                         call msg('ix,iy',ix,iy)
                         call msg('iLev,Mass,',iLev,fmass)
                         call msg("Levels used in advection:",levMinadvection, levMaxadvection)
                         do iTmp = 1, nz_dispersion
                           fTmp =  max( mystuff%passengers_out(0,iTmp),arMinAdvMass(iSpecies))
                           call msg('iLev,Xcm',iTmp,  mystuff%passengers_out(nPass+1,iTmp) / fTmp)
                           call msg('Ycm,Mass', mystuff%passengers_out(nPass+2,iTmp) / fTmp, &
                                                & mystuff%passengers(0,iTmp,iSpecies, iSrc))
                         end do
                         call set_error("StrangeTmp h-moment position after vert advection", "advV4vert")
                         cycle iX_do
                     endif
                   enddo !jTmp X,Y moments

                   ZcmTmp = mystuff%moment_line_out(iLev) !Center of mass, actually...
                   if(.not.  abs(ZcmTmp) <= 0.5) then 
                       call msg(' Strange Zcm rel.cm',  ZcmTmp)
                       call msg('Species: ', iSpecies)
                       call msg('ix,iy',ix,iy)
                       call msg('iLev,Mass,',iLev,fmass)
                       call msg("Levels used in advection:",levMinadvection, levMaxadvection)
                       call msg("---------------------")
                       do iTmp = 1, nz_dispersion
                         !rel CM position
                         call msg('Before/after adv: iLev, mass.'+fu_str(iTmp), &
                           &  mystuff%passengers(0,iTmp,iSpecies,iSrc) , mystuff%passengers_out(0,iTmp)) 
                         call msg('Before/after adv: iLev, rel. dev.'+fu_str(iTmp), &
                           &  mystuff%cm_line(iTmp,iSpecies,iSrc) , mystuff%moment_line_out(iTmp)) 
                         
                         call msg('mystuff%PAbove',mystuff%PAbove(iTmp) )
                         call msg("-----------")
                       end do
                       call set_error("Trouble after advection2! End", "AdvVertV4") 
                       cycle iX_do
                   endif
                 endif ! Mass Ok
              end do ! levels
            endif
   !##############################################
#endif



            !----------------------------------------------------------------------------------
            !
            !   VERTICAL DIFFUSION 
            !
            

            if(do_vert_diff)then  ! Do vertical diffusion
              
if(ifTalk)call msg('Starting vertical diffusion:',(/ix,iy,iLev,iSpecies/))
          
              
              if (diffuse_cm_vert)  mystuff%moment_line_out(1:nz_dispersion) = &
                   & mystuff%cm_relax_diff(1:nz_dispersion) * mystuff%moment_line_out(1:nz_dispersion)





              call prepareDiffColumnPressure( mystuff%g_over_deltaP, mystuff%rho_over_R, & !g_over_deltaP, rho_over_R_up,
                                     &      mystuff%moment_line_out, &
                                     &      mystuff%rho_above, mystuff%PAbove, & ! Diff column
                                     &      mystuff%r_down_met, mystuff%p_met,  & !Meteo  R1m
                                     &    seconds_abs, N_time_steps) ! Dispersion column

              ! The srf-air interfacing resistance is mystuff%rho_over_R(0) but for matrix it is the A-B-C(1) 
              ! that are at the interface. The 0-th elements are for underground layer
              fVd = 0.

              if(mystuff%R_Surf(iSpecies) < -0.5)then
                iLevMinDiff = 1
                ifDryDep = .false.
              else
                iLevMinDiff = 0
                ifDryDep = .true.
              endif

              if (ifDryDep)  then
                ! Get ref height
!                fTmp = (mystuff%Pabove(0) - mystuff%p_cm_out(1)) / (g * mystuff%rho_above(0))
                fTmp = dh1*(0.5 + mystuff%moment_line_out(1) ) !Positive up!
                ftmp = max(ftmp,0.005)  ! half-centimiter min value

                fVd =   fu_get_vd(fTmp,  & ! Reference height in meters
                                & pDispFlds%species(ispecies), & ! what to deposit
                                & indexMeteo, weight_past, &  ! position in space and time
                                & mystuff%R_Surf(ispecies), &    ! surface resistance 
                                & deprules%DryDepType, timeSign, iVd2m) ! Depo rules
              endif


              mystuff%rho_over_R(0) = fVd * mystuff%rho_above(0)   ! same units as rho_over_R_up
              mystuff%rho_over_R(nz_dispersion) = 0. ! hard cover above

#ifdef DEBUG
            ! save masses before the diffusion for debugging purposes
              mystuff%passengers(:,0:nz_dispersion,iSpecies, iSrc) = mystuff%passengers_out(:,0:nz_dispersion)
#endif              

              call make_diffusion_passengers(mystuff%g_over_deltaP(0:nz_dispersion), &  !g_over_deltaP
                                  & mystuff%rho_over_R(0:nz_dispersion), &              ! rho_over_R_up
                                  & mystuff%passengers_out(:,0:nz_dispersion), nPass+2, &
                                  & iLevMinDiff, nz_dispersion, seconds_abs, &
                                  & N_time_steps, report_diff, &
                                  & mystuff%A, mystuff%B, mystuff%C,  mystuff%P, mystuff%Q) 


#ifdef DEBUG
              if (.not. all(mystuff%passengers_out(0,0:nz_dispersion) >= 0)) then
               if (.not. (have_negatives .and. all(abs(mystuff%passengers_out(0,0:nz_dispersion)) >= 0)))then
                 
                call msg("############# Trouble With diffusion! #########################")
                call msg("Negative masses after diffusion!")
                call msg("mystuff%PAbove  ", mystuff%PAbove(0:nz_dispersion))
                call msg("mystuff%rho      ", mystuff%rho_above(0:nz_dispersion+1))
                call msg( "g_over_deltaP", mystuff%g_over_deltaP(0:nz_dispersion))
                call msg("rho_over_R_up", mystuff%rho_over_R(0:nz_dispersion))
                call msg("  Min   ",  mystuff%passengers(0,:, iSpecies, iSrc))
                call msg("  MOut   ",  mystuff%passengers_out(0,:))

                !restore and run again with full reporting
                mystuff%passengers_out(:,0:nz_dispersion) = mystuff%passengers(:,0:nz_dispersion, iSpecies, iSrc)

                call make_diffusion_passengers(mystuff%g_over_deltaP, mystuff%rho_over_R, & ! g_over_deltaP, rho_over_R_up
                                  & mystuff%passengers_out(:,:), 0, & ! mass and passengers here
                                  & iLevMinDiff, nz_dispersion, seconds_abs, &
                                  & N_time_steps, .true., &
                                  & mystuff%A, mystuff%B, mystuff%C, mystuff%P, mystuff%Q) 
                call msg("  MOut after second-try diffusion  ", mystuff%passengers_out(0,0:nz_dispersion))
                call set_error("Negative mass after diffusion!", "Adv vertical V5") 
               endif
             endif
#endif              
!              call msg("################ After Diffusion #######################", ix, iy)
!              call msg ("iSrc, iSpecies",iSrc, iSpecies)
!              call msg("Column mass 3 before, after", fTmp,  sum(mystuff%passengers_out(0,:)) )
!              call msg("p_above  ", mystuff%PAbove(0:nz_dispersion+1))
!              call msg("rho", mystuff%rho_above(0:nz_dispersion) )
!              call msg( "g_over_deltaP", mystuff%g_over_deltaP(0:nz_dispersion))
!              call msg("rho_over_R_up", mystuff%rho_over_R(0:nz_dispersion))
!              call msg( "   MIn   ", mystuff%passengers(0,:,iSpecies, iSrc) )
!              call msg("  MOut   ", mystuff%passengers_out(0,:))
!              call msg("PCmOut  ", mystuff%p_cm_out(0:nz_dispersion+1))
!              call msg("################## After Diffusion  #####################", ix, iy)
                      
! fTmp = sum(mystuff%mass_line)
! call msg("Column mass 4", ispecies, fTmp)
                if (do_molec_diff) then
                  if (stuff%xplus_moldiff(1,iSpecies) > 0.) then
                    call make_molec_diffusion_passengers(stuff%xplus_moldiff(:,iSpecies), &
                               & stuff%xminus_moldiff(:,iSpecies), mystuff%passengers_out(:,:),  &
                               & nz_dispersion, mystuff%Q)
                  endif
                endif

                if (defined(pCnc2m)) then
                ! Diagnose cnc2m coefficient. output contains mass in 1m-thick layer at 2m
                   ! Hcm1
                  if ( (0.5+mystuff%moment_line_out(1)) * dh1  > 2. ) then ! Extrapolate down..
                    fTmp = 1. 
                    if(ifDryDep .and. fVd > 1e-10) fTmp = iVd2m * fVd
                    cnc2mfrac1 = fTmp/dh1  
                    cnc2mfrac2 = 0.
                    weight_up = real_missing
                  else ! interpolate between the levels  
                    weight_up = (2 -  (0.5+mystuff%moment_line_out(1)) * dh1) / &
                       & (dh1 + (0.5+mystuff%moment_line_out(1)) * dh2  - (0.5+mystuff%moment_line_out(1)) * dh1)
                    cnc2mfrac1 = (1. - weight_up) / dh1
                    cnc2mfrac2 = weight_up / dh2
                  endif
                  pCnc2m%arm(iSpecies,iSrc,1,ix,iy) = &
                      & mystuff%passengers_out(0,1) * cnc2mfrac1  + &
                      & mystuff%passengers_out(0,2) * cnc2mfrac2
                endif
               

            else !   Skipping the diffusion
                if(mod(ix,100) == 1 .and. mod(iy,100) == 1) &
                          & call msg_warning('Skip vertical diffusion','adv_diffusion_vertical_v5')
                if (defined(pCnc2m)) then
                  pCnc2m%arm(iSpecies,iSrc,1,ix,iy) = mystuff%passengers_out(0,1) / dh1 ! cnc(2m) = cnc(first level)
                endif
            endif ! ifVertical diffusion



if(ifTalk)call msg('Return mass:',(/ix,iy,iLev,iSpecies/))

            !####################################
            !  Return column to mass maps 
            !####################################

            pDryDep%arM(iSpecies,iSrc,1,ix,iy) = pDryDep%arM(iSpecies,iSrc,1,ix,iy) + &
                                &  mystuff%passengers_out(0,0)

            if (defined(pCnc2m)) then
!!             ! FIXME Hack: pCnc2m  = 1/Rs
!            if (ifDryDep) then
!                   pCnc2m%arM(iSpecies,iSrc,1,ix,iy) =  1./mystuff%R_Surf(ispecies) * pXCellSize(indexDisp) * pYCellSize(indexDisp)
!                  pCnc2m%arM(iSpecies,iSrc,1,ix,iy) =  fVd * pXCellSize(iTmp) * pYCellSize(iTmp)
!            else
!                    pCnc2m%arM(iSpecies,iSrc,1,ix,iy) = 0
!            endif

              if (.not. (pCnc2m%arM(iSpecies,iSrc,1,ix,iy) >= 0)) then !Achtung!!
                if(.not. have_negatives)then
                  call msg("Strange Cnc2m:  pCnc2m%arM(iSpecies,iSrc,1,ix,iy):", &
                      & pCnc2m%arM(iSpecies,iSrc,1,ix,iy) ) 
                  call msg("Resulting concentration (/m3)", &
                      & pCnc2m%arM(iSpecies,iSrc,1,ix,iy) / (pXCellSize(indexDisp) * pYCellSize(indexDisp) )) 
                  call msg("weightUp", weight_up)
                  call msg("dh1, dh2",  dh1, dh2)
                  call msg("zCM(1),zCM(2)", mystuff%moment_line_out(1), mystuff%moment_line_out(2))
                  call msg("M(1), M(2)",  mystuff%passengers_out(0,1), mystuff%passengers_out(0,2))
                  call msg("vd, vd2m", fVd, 1./iVd2m)
                  call msg("Rs", mystuff%R_Surf(iSpecies))
                  call msg("ix, iy", ix,iy)
                  call set_error("cnc_2m", "advV4vert")
                  cycle
                endif
              endif
            endif !defined (pCnc2m)


            ! Return mases and moments with extensive check 
            do  iLev = 1, nz_dispersion  ! Vertical levels
              fmass = mystuff%passengers_out(0,iLev)

              if(abs(fMass) <= arMinAdvMass(iSpecies))then
                ! Mass is small. Drop it to garbage
                mystuff%tmp_garbage(iSpecies,iSrc) = mystuff%tmp_garbage(iSpecies,iSrc)  + fmass
                pDispFlds%arM(iSpecies,iSrc,iLev,ix,iy) = 0.
                moment_x%arM(iSpecies,iSrc,iLev,ix,iy) = 0.
                moment_y%arM(iSpecies,iSrc,iLev,ix,iy) = 0.
                moment_z%arM(iSpecies,iSrc,iLev,ix,iy) = 0.
              else
                !These guys should be moments already
                XcmTmp = mystuff%passengers_out(nPass+1,iLev) / fMass
                YcmTmp = mystuff%passengers_out(nPass+2,iLev) / fMass
                ZcmTmp = mystuff%moment_line_out(iLev) !This is cm (+-0.5) Positive Up

!                ZcmTmp = 0.5 * (mystuff%Pabove(iLev-1) +  mystuff%Pabove(iLev) - 2.*mystuff%p_cm_out(iLev) ) / &
!                                 & (mystuff%Pabove(iLev-1) - mystuff%Pabove(iLev))
                 
                if (have_negatives) then ! reset crazy CMs
                  if(abs(XcmTmp) > 0.5) XcmTmp = 0.
                  if(abs(YcmTmp) > 0.5) YcmTmp = 0.
                  if(abs(ZcmTmp) > 0.5) ZcmTmp = 0.

                else ! Only positive mass allowed -> moments should not misbehave

                  if(.not. abs(XcmTmp) <= 0.5) then 
                    !$OMP CRITICAL(v4bark)
                    call msg('Strange X position ', XcmTmp)
                    call msg('Species: ', iSpecies)
                    call msg('ix,iy',ix,iy)
                    call msg('iLev,Mass,',iLev,fmass)
                    call msg('Column mass, arMinAdvMass(iSpecies),', sum(mystuff%passengers_out(0,:)), arMinAdvMass(iSpecies))
                    call msg("Levels used in advection:",levMinadvection, levMaxadvection)
                    call msg('iLev       Xcm          Ycm         Mass' &
#ifdef DEBUG                      
                                & // '        XCMin        YCMin       MassIn' &
#endif                      
                    &)
                    do iTmp = 0, nz_dispersion+1
                      fTmp =  max( mystuff%passengers_out(0,iTmp),arMinAdvMass(iSpecies))

#ifdef DEBUG                      
                      fTmp1 = max( mystuff%passengers(0,iTmp,iSpecies, iSrc),arMinAdvMass(iSpecies))
#endif                      
                          call msg( fu_str(iTmp), &
                         (/  mystuff%passengers_out(nPass+1,iTmp) / fTmp, &
                             mystuff%passengers_out(nPass+2,iTmp) / fTmp,  mystuff%passengers_out(0,iTmp) &
#ifdef DEBUG                      
                            , mystuff%passengers(nPass+1,iTmp, iSpecies, iSrc) / fTmp1, &
                             mystuff%passengers(nPass+2,iTmp, iSpecies, iSrc) / fTmp1, &
                             mystuff%passengers(0,iTmp, iSpecies, iSrc) &
#endif                      
                             /))
                    end do
                    !$OMP END CRITICAL(v4bark)
                    XcmTmp = sign(0.499, XcmTmp)
                  endif

                  if(.not.  abs(YcmTmp) <= 0.5) then 
                    !$OMP CRITICAL(v4bark)
                    call msg('Strange Y position ', YcmTmp)
                    call msg('Species: ', iSpecies)
                    call msg('ix,iy',ix,iy)
                    call msg('iLev,Mass,',iLev,fmass)
                    !$OMP END CRITICAL(v4bark)
                    YcmTmp = sign(0.499, YcmTmp)
                  endif

                  if(.not.  abs(ZcmTmp) <= 0.5) then 
                   !$OMP CRITICAL(v4bark)
                    call msg('Strange vertical position , rel.cm ', ZcmTmp)
                    call msg('Species: ', iSpecies)
                    call msg('ix,iy',ix,iy)
                    call msg('iLev,Mass,',iLev,fmass)
                    call msg("Levels used in advection:",levMinadvection, levMaxadvection)
                    call msg("---------------------")
                    call msg('mystuff%PAbove(0)', mystuff%PAbove(0))
                    do iTmp = 1, nz_dispersion
                      call msg('Before adv/diff: iLev, rel. dev.', iTmp, mystuff%cm_line(iLev,iSpecies,iSrc))
                      call msg('After  adv/diff: iLev, rel. dev.', iTmp, mystuff%moment_line_out(iLev)) 
                      call msg('mystuff%PAbove',mystuff%PAbove(iTmp))
                      call msg("-----------")
                    end do
                    !$OMP END CRITICAL(v4bark)
                    ZcmTmp = sign(0.499, ZcmTmp)
                  endif
                endif !if  have_negatives
                pDispFlds%arM(iSpecies,iSrc,iLev,ix,iy) = fMass
                moment_x%arM(iSpecies,iSrc,iLev,ix,iy) = XcmTmp
                moment_y%arM(iSpecies,iSrc,iLev,ix,iy) = YcmTmp
                moment_z%arM(iSpecies,iSrc,iLev,ix,iy) = ZcmTmp !* stuff%dh(iLev)
              endif ! Mass Ok
            enddo ! levels

if(ifTalk)call msg('Mass returned:',(/ix,iy,iLev,iSpecies, nPass/))

            ! Return passengers stuff with no checks
            do iPass = 1, nPass
              iTmp = moment_z%passengers(iSpecies)%iSpeciesGarbage(iPass)

              pDryDep%arM(iTmp,iSrc,1,ix,iy) = pDryDep%arM(iTmp,iSrc,1,ix,iy) + &
                                &  mystuff%passengers(iPass, 0, iSpecies, iSrc)

              if (defined(pCnc2m)) pCnc2m%arm(iTmp,iSrc,1,ix,iy) = &
                     & mystuff%passengers(iPass, 1,iSpecies, iSrc) * cnc2mfrac1  + &
                     & mystuff%passengers(iPass, 2,iSpecies, iSrc) * cnc2mfrac2
              !Masses
              do  iLev = 1, nz_dispersion  
                pDispFlds%arm(iTmp,iSrc,iLev,ix,iy) = mystuff%passengers(iPass,iLev,iSpecies, iSrc)
              enddo
            end do  ! passengers


          enddo !sources loop
        enddo !species loop

      end do ix_DO ! Cycle ix
    end do  ! Cycle iy.                 Z-advection-diffusion is over
    !$OMP END DO

    !-------------------------------------------------------------------------------
    !
    ! Southern and Northern polar caps have to be advected and diffused as well. Above, we stored
    ! all needed parameters, here only the actual function calls are needed
    !
    ! Start from collecting info from all threads into the master
    !
    if(any(ifPoles))then
      !
      ! Having the info collected into one place, can do the advection and diffusion
      !
      !$OMP DO
      do ix = northern_boundary, southern_boundary  ! two poles
        
        if(ifPoles(ix))then    ! pole exists and advection-neighbours are > 0
          !
          ! Get the stuff from boundary structures
          !
          ifSkipColumn = .true.
          mystuff%ifSkipSourceSpecies(:,:) = .true.
          do isrc = 1, nSrc
            do iSpecies = 1, nspecies
              nPass = moment_z%passengers(iSpecies)%nPassengers   
              fMinAdvectedMass = arMinAdvMass(ispecies)
              do iLev = 1, nz_dispersion
                !
                ! Get total mass and check the low-cutting threshold
                !
                fmass = pBBuf%pBoundaryData(ix)%ptr(iSpecies,iSrc,1,iLev)
                if (abs(fmass) < fMinAdvectedMass) then    ! Mass is small. Drop all to source-related PFB and continue
                  garbage(isrc,ispecies) = garbage(isrc,ispecies)  + fmass
                  pBBuf%pBoundaryData(ix)%ptr(iSpecies,iSrc,1:2,iLev) = 0.0   ! both mass and z-moment
                  do iPass = 1, nPass
                    iTmp = moment_z%passengers(iSpecies)%iSpeciesGarbage(iPass)
                    garbage(isrc,iTmp) = garbage(isrc,iTmp) + pBBuf%pBoundaryData(ix)%ptr(iTmp,iSrc,1,iLev)
                    pBBuf%pBoundaryData(ix)%ptr(iTmp,iSrc,1:2,iLev) = 0.
                  end do  !  passengers
                  mystuff%cm_line(iLev,iSpecies,iSrc) = 0.
                  mystuff%passengers(0:nPass,iLev,iSpecies,iSrc) = 0.
                
                else  ! Advectable mass
                  ifSkipColumn = .false.
                  mystuff%ifSkipSourceSpecies(iSpecies,iSrc) = .false.
                  mystuff%passengers(0,iLev,iSpecies,iSrc) = fmass
                  !
                  ! z-coord of pole mass always comes as +-0.5, "+" means upwards
                  !
                  mystuff%cm_line(ilev,iSpecies,iSrc) = pBBuf%pBoundaryData(ix)%ptr(iSpecies,iSrc,2,iLev)
                  !
                  ! If negative mass is allowed, cenres can be wild
                  if(have_negatives)then
                    if(abs(mystuff%cm_line(ilev,iSpecies,iSrc)) >= 0.5) &
                                                          & mystuff%cm_line(ilev,iSpecies,iSrc) = 0.
                  endif
                  ! Get passengers
                  do iPass = 1, nPass
                    iTmp = moment_z%passengers(iSpecies)%iSpeciesGarbage(iPass)
                    mystuff%passengers(iPass,iLev,iSpecies, iSrc) = pBBuf%pBoundaryData(ix)%ptr(iTmp,iSrc,1,iLev)
                  end do  ! passengers
                  levMinAdvection=min(levMinAdvection,iLev)
                  levMaxAdvection=max(levMaxAdvection,iLev)
                  pDispFlds%ifColumnValid(iSrc,1,1) = .true.
                end if
              enddo   ! iLev
            enddo !ISpecies
         enddo  !iSrc

         if(ifSkipColumn) cycle

         ! Minor diagnostics.... Should be ready after y advection
         mystuff%wind_right(0) = 0.
         do iLev = 1, nz_dispersion
             mystuff%wind_right(iLev) = sum(pBBuf%wind_to_pole(:,ix,iLev)) + mystuff%wind_right(iLev-1)
             mystuff%cellmass(iLev) = pBBuf%PoleCapAirmass(ix)%pp(iLev)
         enddo

          ! Some Stats (only the current MPI domain)
          !cheap and dirty calculation of courant 
          mystuff%scratch(1:nz_dispersion) = &
                 & abs(seconds*mystuff%wind_right(1:nz_dispersion) / mystuff%cellmass(1:nz_dispersion))
          do ilev = 1, nz_dispersion
             mystuff%cMax1d(iLev) = max(mystuff%scratch(iLev), mystuff%cMax1d(iLev))
             mystuff%Cmean1d(iLev) = mystuff%Cmean1d(iLev) + mystuff%scratch(iLev)
             mystuff%Ccount(iLev) = mystuff%Ccount(iLev) + 1
         enddo

         if ( .not. pDispBuf%ifMassFluxBottomUp) then !subtract vertical wind at top from the whole column
            mystuff%wind_right(1:nz_dispersion) = mystuff%wind_right(1:nz_dispersion) - mystuff%wind_right(nz_dispersion)
         endif

         ! Trick to make reporting N then S and not vice versa 
         if (ix == northern_boundary) then
#ifdef DEBUG
            call msg("massflux at pole N(kg/s)", mystuff%wind_right(0:nz_dispersion))
#endif
         else
            spthread = iThread !keep it to report South Pole  later
         endif

         mystuff%ifSkipCell(0:nz_dispersion+1) = .false.
         call advect_cellboundaries(mystuff%wind_right(0:nz_dispersion), &!Staggered
                     & mystuff%cellmass(1:nz_dispersion),  & ! non-staggered
                     &     mystuff%lineParams(0:nz_dispersion),  & !Staggered
                     & mystuff%lineParamsAdvected(0:nz_dispersion),       &!Staggered
                     &     mystuff%ifSkipCell(0:nz_dispersion+1), & !non-staggered
                     &     mystuff%scratch,   & !staggered
                     &     abs(seconds), .false., levMinadvection, levMaxadvection, nz_dispersion)
          if(error) then
               call msg("Trouble with Z pole")
               cycle
          endif

          do isrc = 1, nSrc
            do iSpecies = 1, nspecies
              nPass = moment_z%passengers(iSpecies)%nPassengers   
              fMinAdvectedMass = arMinAdvMass(ispecies)

              call  MASS_DISTRIBUTOR(mystuff%cm_line(1:,iSpecies,iSrc),  mystuff%moment_line_out(1:), &
                          & mystuff%passengers(:,0:,iSpecies,iSrc), mystuff%passengers_out(:,0:),&
                          & nPass,  & ! number of passengers. No cross-moments
                          & mystuff%lineParams, mystuff%lineParamsAdvected,  mystuff%cellmass,& 
                          & .false., levMinadvection, levMaxadvection, nz_dispersion,&
                          & arMinAdvMass(ispecies), have_negatives, .true., 1.0) !Convert to CM

            ! Add to and bottom mass
            mystuff%bottom_mass(iSrc,iSpecies, outgoing) = mystuff%bottom_mass(nSrc,iSpecies, outgoing) + &
                                                  &  mystuff%passengers_out(0,0)
            mystuff%top_mass(iSrc,iSpecies, outgoing) =  mystuff%top_mass(iSrc,iSpecies, outgoing) + &
                                                  & mystuff%passengers_out(0,nz_dispersion+1) 
            
            do iPass = 1, nPass
               iTmp = moment_z%passengers(iSpecies)%iSpeciesGarbage(iPass)
               mystuff%bottom_mass(iSrc,iTmp, outgoing) =  mystuff%bottom_mass(iSrc,iTmp, outgoing)  + &
                      & mystuff%passengers_out(iPass,0)
               mystuff%passengers_out(iPass,0)  = 0.
               mystuff%top_mass(iSrc,iTmp, outgoing) =  mystuff%top_mass(iSrc,iTmp, outgoing)  + &
                      & mystuff%passengers_out(iPass, nz_dispersion+1)
               mystuff%passengers_out(iPass, nz_dispersion+1)  = 0.
            end do  ! passengers
              !
              ! Return the stuff back to boundary arrays
              !
              !pBBuf%poleDryDep(ix,iSpecies,iSrc) = pBBuf%poleDryDep(ix,iSpecies,iSrc) + mystuff%passengers_out(0,0)
              do  iLev = 1, nz_dispersion  ! Vertical levels
                fmass = mystuff%passengers_out(0,iLev)
                if(abs(fMass) <= arMinAdvMass(iSpecies))then
                  ! Mass is small. Drop it to garbage
                  garbage(iSrc,iSpecies) = garbage(iSrc,iSpecies)  + fmass
                  pBBuf%pBoundaryData(ix)%ptr(iSpecies,iSrc,1:2,iLev) = 0.
                else
                  ZcmTmp = mystuff%moment_line_out(iLev) !This is cm (+-0.5) Positive Up
                  if (have_negatives) then ! reset crazy CMs
                    if(abs(ZcmTmp) > 0.5) ZcmTmp = 0.
                  else
                    if(.not.  abs(ZcmTmp) <= 0.5) then 
                      call msg('Strange vertical position, rel.cm ', ZcmTmp)
                      call msg('Species: ', iSpecies)
                      call msg('Pole (1=north,2=south):',ix)
                      call msg('iLev,Mass,',iLev,fmass)
                      call msg("Levels used in advection:",levMinadvection, levMaxadvection)
                      call msg("---------------------")
                      do iTmp = 1, nz_dispersion
                        call msg('Before adv/diff: iLev, rel. dev.', iTmp, mystuff%cm_line(iLev,iSpecies,iSrc))
                        call msg('After  adv/diff: iLev, rel. dev.', iTmp, mystuff%moment_line_out(iLev)) 
                        call msg("-----------")
                      end do
                      ZcmTmp = sign(0.499, ZcmTmp)
                    endif
                  endif !if  have_negatives
                  pBBuf%pBoundaryData(ix)%ptr(iSpecies,iSrc,1:2,iLev) = (/fMass, ZcmTmp/)
                endif  ! if decent mass
              end do  ! iLev
            enddo   ! iSpecies
          enddo  ! iSrc
        end if ! if this pole is active
      end do  ! cycle over two poles
    !$OMP END DO
!!    !$OMP MASTER
!!     !$OMP END MASTER
    else
!!      !$OMP MASTER 
!!       call msg("No poles -- no vertical pole advection")
!!      !$OMP END MASTER

    endif  ! if max pole index > 0
    !$OMP END PARALLEL

    if(ifTunedRs) call msg_warning('Tuned Rs values used','adv_diffusion_vertical_v5')
    
#ifdef DEBUG
    if (spthread >= 0) then !Report south pole massflux
        mystuff => EulerStuff%Galp5(spthread)
        call msg("massflux at pole S(kg/s)", mystuff%wind_right(0:nz_dispersion))
    endif
#endif


#ifdef DEBUG_V5
call msg('Checking after Z-adv')   ! no tolerance to error, tolerance_factor = 1.0. Put 0.999 to soften
call check_mass_centres_mass_map(pDispFlds, moment_x, 'after Z-adv', 'Px_centre', 1.0)
call check_mass_centres_mass_map(pDispFlds, moment_y, 'after Z-adv', 'Py_centre', 1.0)
call check_mass_centres_mass_map(pDispFlds, moment_z, 'after Z-adv', 'Pz_centre', 1.0)
call msg('Done')
#endif

    ifAllMoments = .false. ! We have returned centers of mass


    ! collect leftovers to the zeroth thread stuff. Needed for timestep reporting
    zstuff => EulerStuff%Galp5(0)
    do  iThread= 1, nThreads-1 !Except for the 0-th one
      mystuff => EulerStuff%Galp5(iThread)
      do iLev = 1, nz_dispersion
         zstuff%Cmax1d(iLev) = max(mystuff%Cmax1d(iLev), zstuff%Cmax1d(iLev))
      enddo 
      zstuff%Cmean1d(:) = zstuff%Cmean1d(:) + mystuff%Cmean1d(:) 
      zstuff%Ccount(:) = zstuff%Ccount(:) + mystuff%Ccount(:)
      zstuff%top_mass(1:nSrc, 1:nSpecies, 1:2) &
           & = zstuff%top_mass(1:nSrc, 1:nSpecies, 1:2) + mystuff%top_mass(1:nSrc, 1:nSpecies, 1:2)
      mystuff%top_mass(1:nSrc, 1:nSpecies, 1:2) = 0. ! Just to keep the budget
      zstuff%bottom_mass(1:nSrc, 1:nSpecies, 1:2) &
           & = zstuff%bottom_mass(1:nSrc, 1:nSpecies, 1:2) + mystuff%bottom_mass(1:nSrc, 1:nSpecies, 1:2)
      mystuff%bottom_mass(1:nSrc, 1:nSpecies, 1:2) = 0. 
      zstuff%tmp_garbage(1:nspecies,1:nSrc) = zstuff%tmp_garbage(1:nspecies,1:nSrc) + mystuff%tmp_garbage(1:nspecies,1:nSrc)
      mystuff%tmp_garbage(1:nspecies,1:nSrc) = 0
    enddo
    do isrc = 1, nSrc
       garbage(isrc,1:nspecies) = garbage(isrc,1:nspecies) + zstuff%tmp_garbage(1:nspecies,iSrc)
    enddo
    ! Add timestep budget to cumulative one
    top_mass(1:nSrc, 1:nSpecies, 1:2) &
           & = top_mass(1:nSrc, 1:nSpecies, 1:2) + zstuff%top_mass(1:nSrc, 1:nSpecies, 1:2)
    bottom_mass(1:nSrc, 1:nSpecies, 1:2) &
           & = bottom_mass(1:nSrc, 1:nSpecies, 1:2) + zstuff%bottom_mass(1:nSrc, 1:nSpecies, 1:2)
  call msg('Max Courant for levels on vert (min, mean, max):', (/ minval(zstuff%Cmax1d(:)), sum(zstuff%Cmax1d(:))/nz_dispersion, maxval(zstuff%Cmax1d(:)) /))
!  call msg('Max Courant for levels:',  zstuff%Cmax1d(:))
!  call msg('Mean Courant for levels:', zstuff%Cmean1d(:)/(zstuff%Ccount(:)+1e-5))
!  call msg('Mean Courant:', sum(zstuff%Cmean1d(:) / (zstuff%Ccount(:)+1e-5))/nz_dispersion )  ! to avoid nans
  call msg('')


  end subroutine adv_diffusion_vertical_v5


  !**********************************************************************************

  subroutine prepareDiffColumnPressure(g_over_deltaP, rho_over_R_up, &
            &  CMCol, rhoCol, pTopCol, &
            &  Ramet, Pmet,  & !Meteo  R1m
            &  seconds_abs, N_time_steps ) ! Dispersion column 
        ! fills g_over_deltaP, rho_over_R_up
        ! Assumes  all pTopCol, are within [Pmet(0) Pmet(nz_meteo+1)]
        !
        real, dimension(0:), intent(out)  :: rho_over_R_up
        real, dimension(0:), intent(in)   :: g_over_deltaP, rhoCol, Ramet, Pmet 
        real, dimension(1:), intent(inout) :: CMCol
        real(r8k), dimension(0:), intent(in) :: pTopCol
        integer, intent(out)  :: N_time_steps
        real, intent(in) :: seconds_abs
        integer  :: iLevMet,   iLev, iTmp
        real :: omegaMax,omegaLimit,myOmega, fTmp, R, frac_here, RhoPrev, CM
        real(r8k) :: Pprev, pCM1   

        integer, parameter :: MaxSteps = 20  ! Hard maximum allowed diffusion steps

        real, parameter :: RecommendedOmegaPlus = 0.2 ! Guide to select timestep
        ! If no cut occurs timestep is RecommendedOmegaPlus of the fastest link
                                            ! 

        real, parameter :: MaxOmegaPlus = 5. ! maximum allowed timestep/tau_link
        ! to start cutting exchange. MaxSteps is always used then. Hell with
        ! accuracy...

       iLevMet = 1
       R=0
       Pprev=Pmet(0)
       frac_here = 1
       RhoPrev = rhoCol(1)
       do iLev = 1, nz_dispersion

         CM =  CMCol(iLev)


         !Pressure for further evaluation of resistance
         pCM1 =  pTopCol(iLev-1) * (0.5 - CM)  + pTopCol(iLev)*(0.5 + CM)

         do while( Pmet(iLevMet) > pCM1) 
           ! find met level just above dispersion one
           R = R + Ramet(iLevMet)*frac_here  ! add  resistance
           Pprev = Pmet(iLevMet)
           frac_here = 1.
           iLevMet=iLevMet+1
         enddo
         frac_here = (pCM1 - Pprev) / (Pmet(iLevMet) - Pmet(iLevMet -1))
         R = R +  Ramet(iLevMet) * frac_here
         if (.not. (R > 0)) then 
           if (R == 0) then
               R = 0.01
           else
             call msg("Strange R", R)
             call msg("pPrev, pCM",Pprev,pCM1  )
             call msg("iLev, iLevMet", iLev, iLevMet)
             call msg("Pmet(iLevMet), Pmet(iLevMet -1), frac_here:", (/Pmet(iLevMet), Pmet(iLevMet -1), frac_here/))
             call msg('RaMet',Ramet)
             call msg('PMet',Pmet)
             call set_error('resistance < 0','prepareDiffColumnPressure')
             call msg("setting 100 s/m")
             R = 0.01
             call unset_error("prepareDiffColumnPressure/")
         endif
         endif
         rho_over_R_up(iLev-1) = 0.5*(rhoCol(iLev) + RhoPrev) / R !Rho_above
         RhoPrev = rhoCol(iLev)
         R = 0. 
         Pprev = pCM1
         !Remaining fraction 
         frac_here = (Pmet(iLevMet) - Pprev) / (Pmet(iLevMet) - Pmet(iLevMet -1))
      enddo
      rho_over_R_up(nz_dispersion) = 0.
 
      !
      ! Cut too fast exchange if needed and select timestep
      !
      omegaMax = 1./ seconds_abs ! Initial value for 1/Tau -- one diffusion step
      omegaLimit = MaxOmegaPlus * MaxSteps * omegaMax  ! 

      do iLev = 1, nz_dispersion-1 ! loop over interfaces
        fTmp = g_over_deltaP(iLev+1) + g_over_deltaP(iLev) !Inverse capacitance part
        myOmega = rho_over_R_up(iLev) * fTmp ! 1/tay
        if (myOmega > omegaLimit) then !cut  too fast exchange
!                  call msg_warning("Adjusting exchange in the column!")
!                  call msg("Cutting exchange. tau, tau_limit", 1./myOmega, 1./omegaLimit)
                  rho_over_R_up(iLev) = omegaLimit / fTmp
                  myOmega = omegaLimit
        endif
        rho_over_R_up(nz_dispersion) = 0. !hard top for diffusion
        omegaMax = max(omegaMax, myOmega)      
      enddo
      N_time_steps = min(nint(seconds_abs * omegaMax / RecommendedOmegaPlus ), MaxSteps)
!      !FIXME
!      call msg_warning("Hacked diffusion", "prepareDiffColumnSimple")
!      N_time_steps = 11.

  end subroutine prepareDiffColumnPressure
  !**********************************************************************************


  subroutine prepareDiffCMvert(CMfactors, rhoCol, pTopCol, &
            &  Ramet, Pmet,  & !Meteo  R1m
            &  seconds_abs) ! Dispersion column 
        ! Coefficient to diffuse vertical centers of mass
        ! Idea: CM relaxes to zero with tau = dz*dz/ Kz = dz*R,
        ! where dz -- cell thickness, R -- resistance between cell bottom and top
        ! Assumes  all pTopCol, are within [Pmet(0) Pmet(nz_meteo+1)]
        !
        real, dimension(1:), intent(out) :: CMfactors ! <= 1  to multiply CM with
        real, dimension(0:), intent(in)   :: rhoCol, Ramet, Pmet 
        real(r8k), dimension(0:), intent(in) :: pTopCol
        real, intent(in) :: seconds_abs

        integer  :: iLevMet, iLev
        real :: tau, R, frac_here

       iLevMet = 1
       frac_here = 1.
       do iLev = 1, nz_dispersion
         R = 0. 
         do while( Pmet(iLevMet) >  pTopCol(iLev)) 
           ! met. level below disp. layer top
           R = R + Ramet(iLevMet)*frac_here  ! add  resistance
           frac_here = 1.
           iLevMet=iLevMet+1
         enddo
         frac_here = (pTopCol(iLev) - Pmet(iLevMet-1)) / (Pmet(iLevMet) - Pmet(iLevMet -1))
         R = R +  Ramet(iLevMet) * frac_here

         tau = R * (pTopCol(iLev-1) -  pTopCol(iLev)) /(rhoCol(iLev) * g) 
         CMfactors(iLev) = exp(-seconds_abs/tau)
         if (CMfactors(iLev) > 1) then
            call msg("tau, CMfactors(iLev)", tau, CMfactors(iLev))
            call set_error("Gotcha wrong relaxation for CM", "prepareDiffCMvert")
         endif
         frac_here = (Pmet(iLevMet) - pTopCol(iLev)) / (Pmet(iLevMet) - Pmet(iLevMet -1))
      enddo

  end subroutine prepareDiffCMvert
            
            
  !***********************************************************************************

  subroutine make_diffusion_passengers(g_over_deltaP, rho_over_R,  passengers, nPass, &
                          & iLevmin, iLevMax, seconds_abs, N_time_steps,ifReport, &
                          & A, B, C, P, Q)  
      !The routine is supposed to work with any vertical
      real, dimension(0:), intent(in) :: g_over_deltaP, rho_over_R
      real, dimension(0:, 0:), intent(inout) :: passengers
      real, intent(in) :: seconds_abs
      integer, intent(in) :: iLevmin, iLevMax, nPass
      logical, intent(in) :: ifReport
      real(r8k) :: fMAfter, fMbefore, fMBeforeAbs, Tau
      real(r8k) :: xminuscur, xminusnext, xplusprev, xpluscur, fTmp
      integer  :: iLevMet,iLev,  iTmp, jTmp, iPass 
      integer, intent (in) :: N_time_steps
      
      ! No need to output  Just not to torture stack... 
      real(r8k), dimension(0:), intent(out) :: A, B, C, P
      real(r8k), dimension(0:,0:), intent(out) :: Q 

      Tau = seconds_abs / real(max(1,N_time_steps))
      
      ! Fill-in three-diagonal matrix.
      ! See notebook P. 48 (Roux)
      xminuscur = 0.      ! Nothing below underground level
      xplusprev = 0.
      do iLev = iLevMin, iLevMax - 1
        xminusnext =  Tau * rho_over_R(iLev) * g_over_deltaP(iLev+1) 
        xpluscur   =  Tau * rho_over_R(iLev) * g_over_deltaP(iLev)  
        A(iLev) =  - xplusprev
        B(iLev) = 1. + xminuscur + xpluscur
        C(iLev) =  - xminusnext
        if (.not. (B(iLev) > 0)) then
          call msg("Negative B, iLev", B(iLev), iLev) 
          call set_error("Too bad", "make_diffusion_passengers")
          call msg("xminuscur xplusprev", xminuscur, xplusprev)
          call msg("xpluscur, xminusnext", xpluscur, xminusnext)
          call msg("Levels", 0, iLevMax)
          call msg("rho_over_R",  rho_over_R(0:iLevMax))
          call msg("g_over_deltaP",   g_over_deltaP(0:iLevMax))
          return
        endif
        ! prepare for next level
        xminuscur =  xminusnext
        xplusprev = xpluscur
     enddo
     A(iLevMax) = - xplusprev
     B(iLevMax) =   1. + xminuscur
     C(iLevMax) =   0.

     fMBefore =  sum( passengers(0,iLevMin:iLevMax))
     fMBeforeAbs = sum( abs(passengers(0,iLevMin:iLevMax)))

     if (ifReport) then
        call msg('----------------------------------------------------')
        call msg("Timestep, Ntimesteps", real(Tau), 1.0*N_time_steps )
        call msg('Mass Before:', passengers(0,0:iLevMax))
!        call msg('CnC Before []:', passengers(0,0:iLevMax) / g_over_deltaP(0:iLevMax))
        call msg('A:' ,A(0:iLevMax))
        call msg('B:' ,B(0:iLevMax))
        call msg('C:' ,C(0:iLevMax))
     end if
     
     do  jTmp=1, N_time_steps ! Iterate timesteps
          ! Forward part of swift-cycle
          P(iLevMin) = -C(iLevMin) / B(iLevMin)
          Q(:,iLevMin) = passengers(:,iLevMin) / B(iLevMin)
          do iTmp = iLevMin, iLevMax-2
            P(iTmp+1) = -C(iTmp+1) / (B(iTmp+1)+P(iTmp)*A(iTmp+1))
            Q(:,iTmp+1) = (passengers(:,iTmp+1)-Q(:,iTmp)*A(iTmp+1)) / &
                      & (B(iTmp+1) + P(iTmp)*A(iTmp+1))
          end do
          Q(:,iLevMax) = (passengers(:,iLevMax) - &
                & Q(:,iLevMax-1) * A(iLevMax)) / &
                 & (B(iLevMax) + P(iLevMax-1) * A(iLevMax))
          ! Backward (returning) part of swift cycle
          passengers(:,iLevMax) = Q(:,iLevMax)
          do iTmp=iLevMax-1, iLevMin, -1
            passengers(:,iTmp) = Q(:,iTmp) + P(iTmp) * passengers(:,iTmp+1)
          end do
     enddo



     if (ifReport) then
        call msg('MassAfter:', passengers(0,0:iLevMax))
        call msg('CnCAfter:', passengers(0,0:iLevMax) / g_over_deltaP(0:iLevMax))
      end if

   !
   ! Check mass conservation in the column
   !
      fMAfter = sum(passengers(0,iLevMin:iLevMax))

      if(.not. ((abs(fMAfter - fMBefore) <= fMBeforeAbs * 1e-4)))then 
!      if(.not. ((abs(fMAfter - fMBefore) <= abs(fMAfter + fMBefore) * 1e-4)))then 
        !
        ! Something dangerous: set an error. Since the in-surface layer is introduced
        ! at iLevMin (for depositing species only!) the grand total must stay the same.
        ! The dry deposition is then just a change of mass at this in-surface layer.
        !
        call msg('Problem with passenger')
        call msg('Mode, mass before, befreAbs, after diffusion/dd:', (/fMBefore, fMBeforeAbs, fMAfter/))
        call msg('Difference (after-before):',fMAfter - fMBefore)
        do iLev = iLevMin, iLevMax  ! Vertical levels
          call msg('Level and value:',iLev,passengers(0,iLev))
        end do
        call msg('rho/R [kg/m2s]:', rho_over_R(0:iLevMax))
        call msg('g/deltaP [m2/kg]:', g_over_deltaP(0:iLevMax))
        call msg('A:', A(0:iLevMax))
        call msg("B:", B(0:iLevMax))
        call msg('C:', C(0:iLevMax))
        call set_error('Large change of mass during diffusion/drydep', 'make_diffusion_passengers')
        return
      endif    ! if substantial mass gain/loss in the closed column

  


  end subroutine make_diffusion_passengers

      
  !******************************************************************************

  subroutine advect_cellboundaries(u_right, dx_cell, xr, xr_advected,& 
       & ifSkipCell, passtime, abs_seconds, ifLoop, iCellStart, iCellEnd, nCells)

    ! Warning!! If ifLoop==.true. and u_right(0) /= u_right(nCells)  catastrophy happens

    real, dimension(0:), intent(in) :: u_right ! kg/s rightwards flux at right cellside (0:nCells) 
    real, dimension(1:), intent(in) :: dx_cell ! mass of air in cell (1:nCells)
    real, dimension(0:), intent(out) :: passtime ! Actually, scratch area (0:nCells)
    ! passtime(ix) -- how long will it take to right side if ix cell to pass one cell in 
    ! a direction  of u_right(ix)
    logical, dimension(0:), intent(in) :: ifSkipCell  ! zero-mass flag incl. off-domain cell (0:nCells+1)
    real(r8k), dimension(0:), intent(out) :: xr ! pos. of  right border in kg before (0:nCells)
    real(r8k), dimension(0:), intent(out) :: xr_advected  ! pos. of right border in kg after  (0:nCells)
    integer, intent(in) :: iCellStart, iCellEnd, nCells   ! Last index of dx_cell array
    logical, intent(in) :: ifLoop
    real, intent(in) :: abs_seconds

    ! Local variables
    real(r8k) ::  curpos, nextpos, incr, dx, alpha, t, alphat
    integer ::  ix, curix, nextix, dir, iFirstCell, iLastCell
    logical :: useExp, ifSkipPrev

   !
   ! Gets displacement of cell boundaries after advection
   ! trick: for backward time call with opposite-sign  mass-fluxes
   xr_advected(0:nCells) = D_NAN
   passtime(0:nCells) = 2*abs_seconds !Larger than timestep
      
   if (abs_seconds < 0.1) then
      call set_error("abs_seconds should be positive!","advect_cellboundaries")
   endif

   if (ifLoop) then !Advect all cell boundaries
        ifSkipPrev = ifSkipCell(nCells)
        iFirstCell = 0
        iLastCell = nCells-1
      if (u_right(0) /= u_right(nCells)) then
            call msg("u_right(0) /= u_right(nCells) in global domain", u_right(0),u_right(nCells))
        call set_error("Must be the same!", "advect_cellboundaries")
            return
        endif
   else
        ifSkipPrev = .true. ! Current cell decides
        iFirstCell = max(iCellStart-1,0)  
        iLastCell  = min(iCellEnd, nCells)
   endif

   !Get coordinates for right bound whole line is needed
   xr(0)=0.0_r8k
   do ix = 1, nCells
      xr(ix) = xr(ix-1) + dx_cell(ix)
   enddo


    mainloop: do ix = iFirstCell, iLastCell ! Cells, which right side needs advection
      if (ifSkipPrev .and. ifSkipCell(ix)) then
         ifSkipPrev = .true.
        cycle mainloop
      endif
      ifSkipPrev = ifSkipCell(ix)

      t=abs_seconds ! Remaining time
      curix = ix
      dir = int(sign(1.0, u_right(ix))) ! +1 - go right
                                        ! -1 - go left
      curpos = 0.

      innerloop: do while (.true.)
        nextix = curix + dir

        !off-domain left?
        if (nextix < 0) then 
            if (ifLoop) then
               curix = nCells
               nextix = nCells-1
            else
               xr_advected(ix) = u_right(0) * t
            exit innerloop
            endif
         endif

         !off-domain right?
         if (nextix > nCells) then
            if (ifLoop) then
               curix = 0
               nextix = 1
            else
               xr_advected(ix) = xr(nCells) + u_right(nCells) * t
            exit innerloop
            endif
         endif

         dx = xr(nextix) - xr(curix)
         if (passtime(curix) < t) then ! Will not stop here
            t = t - passtime(curix)
            curpos = curpos + dx 
            curix = nextix
          cycle innerloop
         endif

         alpha = (u_right(nextix) - u_right(curix)) / dx ! dw/dz
         alphat= alpha * t
         useExp = (abs(alphat) > 0.00001_r8k )

         if (useExp) then
             incr = u_right(curix) * (exp(alphat) - 1.0_r8k)/alpha
         else
             incr = u_right(curix) * t
         endif

         if (abs(incr - dx) < 1e-8*abs(dx)) then !! Almost never happens, 
                                 !!but causes numerical issues on top-down diagnostics
            xr_advected(ix) = xr(ix) + curpos + dx
            exit
         endif

         if (abs(incr) > abs(dx)) then 
                 ! Move over more than one grid cell
                 ! How long it takes to pass the current cell?
                 if (useExp) then
                      passtime(curix) = log(1.0_r8k + alpha * dx / u_right(curix)) / alpha
                 else
                      passtime(curix) = dx / u_right(curix)
                 endif
                 t = t - passtime(curix)
                 curpos = curpos + dx 
                 curix = nextix
         else
                 xr_advected(ix) = xr(ix) + curpos + incr
          exit innerloop
         endif

      end do innerloop ! expiring timestep
    end do mainloop ! levels

    if(ifLoop) xr_advected(nCells) = xr(nCells) + xr_advected(0)

    if (.not. all(abs(xr_advected(iFirstCell:iLastCell))>=0.)) then
            call msg("")
            call msg("iFirstCell:iLastCell)", iFirstCell, iLastCell)
            call msg("xr(0:nCells)         :",xr(0:nCells))
            call msg("xr_advected(0:nCells):",xr_advected(0:nCells))
            call msg("u_right(0:nCells)*sec:",u_right(0:nCells)*abs_seconds)
            call msg("cellmass(1:nCells)        :",dx_cell(1:nCells))
       call set_error("Gotcha: non-finite xr_advected","here")
    endif

    do ix=iFirstCell+1, iLastCell 
         if (ifSkipCell(ix)) cycle
         if (xr_advected(ix-1) > xr_advected(ix)) then
            call msg("")
            call msg("xr         :",xr(0:nCells))
            call msg("xr_advected:",xr_advected(0:nCells))
            call msg("u_right*sec:",u_right(0:nCells)*abs_seconds)
            call msg("cellmass        :",dx_cell(1:nCells))
            call set_error("Non_monotonous advected cells","")
            exit
         endif
    enddo

!    call msg("")
  end subroutine advect_cellboundaries


!*********************************************************************************************

  subroutine advect_mass(   cmRelIn,  momOut, &
                          & PassengersIn,  PassengersOut, nPass, &
                          & xr, xr_advected, dx,& 
                          & ifLoop, iCellStart, iCellEnd, nCells,&
                          & fMinAdvectedMass, have_negatives, ifCMout, cmrelax)
    ! 
    ! Redistributes masses and moments according to  rightbound_offset.
    ! Advects whatever non-zero mass it finds. Garbage should be collected elsewhere
    ! Warning!! ifLoop loops around line indices as   ..., nCells-1, nCells,  1,  2, ...
    ! ifCmout -- output CM instead of moments
    ! Roux

    real, dimension(1:),    intent(in)  :: cmRelIn ! 1:nCells -- domain only
    real, dimension(1:),    intent(out) :: momOut  ! 1:nCells -- domain only
    real, dimension(1:),    intent(in)  :: dx !(1:nCells) !Cell size (e.g. in kg)
    real(r8k), dimension(0:),    intent(in)  :: xr, xr_advected !(0:nCells)  !same units as above!!!
    real, dimension(0:,0:), intent(in)  :: PassengersIn  !(0:npass, 0:nCells+1) 
    real, dimension(0:,0:), intent(out) :: PassengersOut !(0:npass, 0:nCells+1)
    real, intent(in) :: cmrelax ! Facttor to introduce diffusion bu CM relaxation

! Print full redistribution info: iFrom, iTo , fraction    
!!!#define REPORT(II,frac) call msg("II: ix="+fu_str(ix), ixto, frac)
#define REPORT(II,frac)


#ifdef DEBUG_V5
#define CHECK(what, bark) if (what) call set_error(bark,"advect_mass");  
#define CHECKEXIT(what, bark) if (what) call set_error(bark,"advect_mass"); if (error) exit 
#else
#define CHECK(what, bark)
#define CHECKEXIT(what, bark)
#endif

    ! nlev == nz_dispersion arrays should have one more element
    integer, intent(in) :: iCellStart, iCellEnd, nCells, nPass
    real, intent(in) ::  fMinAdvectedMass ! Used only for conversion to centers of mass
    logical, intent(in) :: have_negatives, ifLoop, ifCMout

    real(r8k) ::  fZC, SS, BR_abs, BL_abs, BL_tmp, BR_tmp, ftmp,  loopx
    real :: fmass, fMass1
    integer ::  ix, ixto
    

    ! !No advection hack
    if (.False.) then
       PassengersOut(:,0:nCells+1)=PassengersIn(:,0:nCells+1)
       if (ifCMout) then
         momOut(1:nCells)= cmRelIn(1:nCells)
       else
         momOut(1:nCells)= cmRelIn(1:nCells)*PassengersIn(0,1:nCells)
       endif
       return
    endif

    momOut(1:nCells)=0.0
    passengersOut(:,0:nCells+1) = 0
    
    ! find first advectable mass
    do ix = iCellStart, iCellEnd
      if (abs(PassengersIn(0,ix)) > 0.) exit
    enddo
    

    !Init distribution loop
    loopx = 0.
    ixto = 0
    if (ix == 0) then !Some left boundary stuff to distribute
       SS  = xr_advected(0) ! Should be positive!!!
       CHECK(.not. (SS > 0), "Can't inject left boundary")
       ixto = 1
       do while (SS > xr(ixto)) 
        PassengersOut(:,ixto) = dx(ixto)/SS * PassengersIn(:,0)
       REPORT (1,dx(ixto)/SS)
        ixto = ixto + 1 ! Moments are zero anyhow
      enddo
      fTmp = xr_advected(0)-xr(ixto-1)  !size of leftover
      PassengersOut(:,ixto) =  fTmp / SS * PassengersIn(:,0) !Remaining part
      REPORT (2, fTmp / SS)
!        call msg("2:ix, frac", 0, fTmp / SS)
      momOut(ixto) =  - 0.5* (1. - fTmp/dx(ixto)) * PassengersOut(0,ixto)  !left-justified slab
      ix = 1
   elseif (xr_advected(ix-1) < 0) then !First mass below the domain
      if (ifLoop) then 
         loopx =  xr(nCells) !Add this part to advected quantities
         do ixto = nCells,1,-1   !Find ixTo for cell ix 
           if (xr_advected(ix-1) + loopx > xr(ixto-1)) exit
         enddo
      endif
   endif


   ! Main distribution loop
   do ix = ix, min(iCellEnd,nCells)
     !if (abs(PassengersIn(0,ix)) < fMinAdvectedMass) cycle
     if (PassengersIn(0,ix) == 0.) cycle

     fZC = cmRelIn(ix) * cmrelax

     if (abs(fZC) >= 0.5 ) then  ! abnormal value check
       if (have_negatives) then 
         call msg('Strange mass centre, return to cell:', fZC)
         fZC = 0.
       else
         call msg('Strange mass centre:', fZC)
         if (abs(fZC) < 0.5001)then ! Still ok
                 fZC  = 0.999*fZC
         else
           call set_error('Strange mass centre, ix=' + fu_str(ix),'vert_adv_column_pass')
           return
         endif 
       endif
     endif
     
     ftmp =  1.0_r8k -  2*abs(fZC) ! rel. size of the slab

     ! Slab is not allowed to occupy less than 1e-3 of cell,
     ! otherwise we get into numerics. Indeed, the target slab should 
     ! be checked, but it is more complicated...
     ftmp = max(ftmp, 1e-6)


     if (fZC < 0.) then
         ! on the lower side of the layer
         bl_abs =  xr_advected(ix-1)
         br_abs = ftmp*xr_advected(ix) + (1.0-ftmp)*bl_abs
     else
         br_abs = xr_advected(ix)
        bl_abs = ftmp * xr_advected(ix-1) + (1.0-ftmp) * br_abs
     end if

     !Still within the domain?
     if (BL_abs + loopx >= xr(nCells))  then 
        if (ifLoop) then
          loopx = loopx - xr(nCells)
          ixto = 1
        else
           !Whole slab and following slabs go above the domain
            ixTo = nCells+1
            passengersOut(0:nPass,ixTo) = passengersOut(0:nPass,ixTo) + &
                              &  passengersIn(0:nPass, ix) 
            REPORT (7, fTmp )
            cycle
        endif
     endif
     
     BL_tmp = BL_abs + loopx

     do  ixto = ixto, nCells
        if (BL_tmp < xr(ixto)) exit
     enddo

     !now BL is just above iLevTo
     SS = BR_abs - BL_abs !new size. can be negative for off-domain target cells

     CHECK(.not. SS>0, "non-positive slab")

     do while (BR_abs + loopx > xr(ixto)) 
           ! Except for rightmost part of advected slab
           BR_tmp = xr(ixto)
           ftmp =  (BR_tmp - BL_tmp) /SS  !Fraction that goes here
           CHECKEXIT( .not. ftmp>=0 , "non-positive slab1")
           passengersOut(0:nPass,ixTo) = passengersOut(0:nPass,ixTo) + &
               &  passengersIn(0:nPass, ix) * ftmp
!             call msg("3:ix, frac", ix, fTmp)
           REPORT (3, fTmp )
           if (ixto > 0) then
               momOut(iXto) = momOut(iXto) +  passengersIn(0, ix) * ftmp * &
                      0.5* (1. - (BR_tmp - BL_tmp)/dx(ixto)) !right-justified slab
               CHECKEXIT(abs(momOut(iXto)) > 0.5*passengersOut(0,ixTo), "Gotcha1")
           endif
           BL_tmp = xr(ixto)
           ixto = ixto + 1
           if (ixto > nCells) then 
              if(.not. ifLoop)  exit
              loopx = loopx - xr(nCells)   !handle loopover
              BL_tmp = 0.0_r8k !  xr(0)  
              ixto = 1
           endif
      enddo
#ifdef DEBUG_V5
      if (error) exit
#endif
      !Leftover
      BR_tmp = BR_abs + loopx
      if (SS > 0.) then
         ftmp =  (BR_tmp - BL_tmp) /SS  !Fraction that goes here
      else
         ftmp = 1.
      endif
      passengersOut(0:nPass,ixTo) = passengersOut(0:nPass,ixTo) + &
                                 & passengersIn(0:nPass, ix) * ftmp
!       call msg("4:ix, frac", ix, fTmp)
       REPORT (4, fTmp )
      if (ixto <= nCells .and.  ixto> 0) then  !non-justified slab
          fTmp =  fTmp * 0.5_r8k* ( BR_tmp + BL_tmp - xr(ixto) - xr(ixto-1))
          fTmp = fTmp / dx(ixto) !!CM * massfraction
          momOut(iXto) = momOut(iXto) + passengersIn(0, ix) * ftmp 
          CHECK(abs(momOut(iXto)) > 0.5*passengersOut(0,ixTo), "Gotcha2")
       endif

    enddo !ix

    !CHECK( (.not. ifLoop) .or. ix == iCellEnd+1, "No leftovers allowed in looped distribution")

    if ( PassengersIn(0,nCells+1) /= 0. &
#ifdef DEBUG_V5
               & .and. .not. error  &
#endif               
                & ) then !get stuff from right boundary
       ! PassengersIn(0,nCells+1) should be non-zero only if something to be injected
       !out_of_domain source cell, inflow, indomain target cell
       SS  = xr(nCells) - xr_advected(nCells) ! Should be positive!!!
       CHECK( .not.(ix == nCells+1 .and. SS > 0.  .and. ixto <= nCells), "Can't inject right boundary")

       do  ixto = ixto, nCells
          if (xr_advected(nCells) < xr(ixto)) exit
       enddo
       
       fTmp = xr(ixto) - xr_advected(nCells)  !first, incomplete slab from boundary
       PassengersOut(:,ixto) =  PassengersOut(:,ixto) + &
              & fTmp / SS * PassengersIn(:,nCells+1)
!        call msg("5:ix, frac", ix, fTmp / SS )
       REPORT (5, fTmp/SS )
       momOut(ixto) = momOut(ixto) + 0.5* (1. - fTmp/dx(ixto)) * & !right-justified slab
            &  fTmp / SS * PassengersIn(0,nCells+1)

       CHECK(abs(momOut(iXto)) > 0.5*passengersOut(0,ixTo), "Gotcha3")
       do ixto = ixto+1,nCells 
         PassengersOut(:,ixto) = dx(ixto)/SS * PassengersIn(:,nCells+1)
!         call msg("6:ix, frac", ix, SS/dx(ixto))
         REPORT (6, 1.0 )
       enddo
    endif


    if (ifCmOut) then 
     ! Turn  Moments to CM within domain
      do ix = 1, nCells 
        fmass = passengersOut(0,ix) 
        if (abs(fmass) >  fMinAdvectedMass) then
          momOut(ix) = momOut(ix) /  fmass  ! CM: +-0.5 , positive -- Up
          if(abs(momOut(ix)) > 0.4999)then 
            if(have_negatives)then      ! Crazy centres are legal...
              momOut(ix) = 0.0
            elseif (abs(momOut(ix)) < 0.5001)then 
               momOut(ix) = sign(0.499,momOut(ix))
            else  !Too bad
               call set_error("Wrong CM at advection output, ix="+fu_str(ix), "advect_mass")
            endif
          endif
        else
          momOut(ix) = 0.0
        endif
      enddo
#ifdef DEBUG_V5 
   else
     ! Turn  Check Moment
      do ix = 1, nCells 
        fmass = passengersOut(0,ix) 
        if (abs(fmass) >  fMinAdvectedMass) then
          fTmp = abs(momOut(ix) /  fmass)  ! CM: +-0.5 , positive -- Up
          if(.not. (fTmp <= 0.5) )then
               call msg("fmass, momOut(ix)",fmass, momOut(ix))
            call set_error("Wrong Moment at advection output, ix="+fu_str(ix), "advect_mass")
          endif
         endif
       enddo
#endif
   endif

#ifdef DEBUG_V5
!    fMass = sum(passengersIn(0,0:nCells+1))
!    fMass1 = sum(passengersOut(0,0:nCells+1))
    fMass = sum(passengersIn(0,iCellStart:iCellEnd))
    fMass1 = sum(passengersOut(0,0:nCells+1))
    if (.not. all(passengersOut(0,0:nCells+1)>=0.) .or. error .or. &
        & (abs(fMass - fMass1) > 1e-7*ncells*(fMass+fMinAdvectedMass)) )then
      !$OMP CRITICAL(v4bark)
       call msg("Trouble after mass advection!" )
       call msg("Line mass in, out", fMass, fMass1)
       call msg("Line mass difference", fmass-fmass1)
       call msg("Mass, in  (0:nCells+1)", passengersIn(0,0:nCells+1))
       call msg("Mass, out (0:nCells+1)", passengersOut(0,0:nCells+1))
       call msg("xr,         (0:nCells):", xr(0:nCells))
       call msg("xr_advected (0:nCells):", xr_advected(0:nCells))
       call msg("cm, in (1:nCells)",   cmRelIn(1:nCells))
       if (ifCmOut) then
          call msg("cm, out (1:nCells)", momOut(1:nCells))
       else
          call msg("mom, out (1:nCells)", momOut(1:nCells))
       endif
       call set_error("Trouble after mass advection!","advect_mass")
      !$OMP END CRITICAL(v4bark)
    endif
#endif
#undef REPORT
#undef CHECK
#undef CHECKEXIT

  end subroutine advect_mass



!*********************************************************************************************

  subroutine advect_mass_step(   cmRelIn,  momOut, &
                          & PassengersIn,  PassengersOut, nPass, &
                          & xr, xr_advected, dx,& 
                          & ifLoop, iCellStart, iCellEnd, nCells,&
                          & fMinAdvectedMass, have_negatives, ifCMout, cmrelax)

    ! Ideal for steps with non-zero background

    ! 
    ! Redistributes masses and moments according to  rightbound_offset.
    ! Advects whatever non-zero mass it finds. Garbage should be collected elsewhere
    ! Warning!! ifLoop loops around line indices as   ..., nCells-1, nCells,  1,  2, ...
    ! ifCmout -- output CM instead of moments
    ! Roux

    real, dimension(1:),    intent(in)  :: cmRelIn ! 1:nCells -- domain only
    real, dimension(1:),    intent(out) :: momOut  ! 1:nCells -- domain only
    real, dimension(1:),    intent(in)  :: dx !(1:nCells) !Cell size (e.g. in kg)
    real*8, dimension(0:),    intent(in)  :: xr, xr_advected !(0:nCells)  !same units as above!!!
    real, dimension(0:,0:), intent(in)  :: PassengersIn  !(0:npass, 0:nCells+1) 
    real, dimension(0:,0:), intent(out) :: PassengersOut !(0:npass, 0:nCells+1)
    real, intent(in) :: cmrelax ! Facttor to introduce diffusion bu CM relaxation

    ! nlev == nz_dispersion arrays should have one more element
    integer, intent(in) :: iCellStart, iCellEnd, nCells, nPass
    real, intent(in) ::  fMinAdvectedMass ! Used only for conversion to centers of mass
    logical, intent(in) :: have_negatives, ifLoop, ifCMout

    real*8 ::  fZC, SS, BR_abs, BL_abs, BL_tmp, BR_tmp, ftmp,  loopx
    real*8 ::  bgmass, bgfrac, exZc, bgdensl, bgdensr
    real*8, dimension(2) :: fracmass
    real*8, dimension(3) :: bounds
    real :: fmass, fMass1
    integer ::  ix, ixto, isubSlab

! Print full redistribution info: iFrom, iTo , fraction    
!#define REPORT(II,frac) print *, "II: ix=", ix, ixto, frac
!#define REPORT(II,frac) call msg("II: ix="+fu_str(ix), ixto, frac)
#define REPORT(II,frac)

#ifdef DEBUG_V5
#define CHECK(what, bark) if (what) call set_error(bark,"advect_mass_step");  
#define CHECKEXIT(what, bark) if (what) call set_error(bark,"advect_mass_step"); if (error) exit 
#else
#define CHECK(what, bark)
#define CHECKEXIT(what, bark)
#endif

    
!print *, ""
!print *, ""
!print *, ""
if(error) call set_error("Secondary error","advect_mass_step")
!call unset_error("Let it crash second_time, advect_mass_step_debug")

    ! !No advection hack
    if (.False.) then
       PassengersOut(:,0:nCells+1)=PassengersIn(:,0:nCells+1)
       if (ifCMout) then
         momOut(1:nCells)= cmRelIn(1:nCells)
       else
         momOut(1:nCells)= cmRelIn(1:nCells)*PassengersIn(0,1:nCells)
       endif
       return
    endif

    momOut(1:nCells)=0.0
    passengersOut(:,0:nCells+1) = 0
    
    ! find first advectable mass
    do ix = iCellStart, iCellEnd
      if (PassengersIn(0,ix) == 0.) cycle
      exit
    enddo
    

    !Init distribution loop
    loopx = 0.
    ixto = 0
    if (ix == 0) then !Some left boundary stuff to distribute
       SS  = xr_advected(0) ! Should be positive!!!
       CHECK(.not. (SS > 0), "Can't inject left boundary")
       ixto = 1
       do while (SS > xr(ixto)) 
        PassengersOut(:,ixto) = dx(ixto)/SS * PassengersIn(:,0)
       REPORT (1,dx(ixto)/SS)
        ixto = ixto + 1 ! Moments are zero anyhow
      enddo
      fTmp = xr_advected(0)-xr(ixto-1)  !size of leftover
      PassengersOut(:,ixto) =  fTmp / SS * PassengersIn(:,0) !Remaining part
      REPORT (2, fTmp / SS)
      momOut(ixto) =  - 0.5* (1. - fTmp/dx(ixto)) * PassengersOut(0,ixto)  !left-justified slab
      ix = 1
   elseif (xr_advected(ix-1) < 0) then !First mass below the domain
      if (ifLoop) then 
         loopx =  xr(nCells) !Add this part to advected quantities
         do ixto = nCells,1,-1   !Find ixTo for cell ix 
           if (xr_advected(ix-1) + loopx > xr(ixto-1)) exit
         enddo
      endif
   endif


   !Some preparation of background info for boundaries
   bgdensl = 0D0
   if (cmRelIn(1) > 0. .and. PassengersIn(0,1) > 0.)  then
      if (ifLoop) then
         bgdensl = PassengersIn(0,nCells) * (1. + 2*cmRelIn(nCells))  / (xr_advected(nCells) - xr_advected(nCells-1))
         !bgdensl = PassengersIn(0,nCells)   / (xr_advected(nCells) - xr_advected(nCells-1))
      else
         if (xr_advected(0)>0D0) bgdensl = PassengersIn(0,0)  / xr_advected(0)
      endif
   endif
   bgdensr = 0D0
   if (cmRelIn(nCells) < 0 .and. PassengersIn(0,nCells) > 0.)  then
      if (ifLoop) then
         bgdensr = PassengersIn(0,1) * (1. + 2*cmRelIn(1))  / (xr_advected(1) - xr_advected(0))
         !bgdensr = PassengersIn(0,1)  / (xr_advected(1) - xr/_advected(0))
      else
         if (xr_advected(nCells) < xr(nCells)) &
            & bgdensr = PassengersIn(0,nCells+1)  / (xr(ncells) - xr_advected(nCells))
      endif
   endif
!   print *, "bgdensl, bgdensr", (/bgdensl, bgdensr/) 


   ! Main distribution loop
IXLOOP:   do ix = ix, min(iCellEnd,nCells)
     if (PassengersIn(0,ix) == 0.) cycle
!     fMass = sum(PassengersIn(0,0:ix-1))
!     fMass1 = sum(PassengersOut(0,0:nCells+1))
!     print *, "ix=", ix, fMass, fMass1, abs(fMass-fMass1)/fMass
     
     fZC = cmRelIn(ix) * cmrelax
     if (abs(fZC) < 1e-5) then 
        fZC = 0.
     elseif (abs(fZC) >= 0.5 ) then  ! abnormal value check
       if (have_negatives) then 
         call msg('Strange mass centre, return to cell:', fZC)
         fZC = 0.
       else
         call msg('Strange mass centre:', fZC)
         if (abs(fZC) < 0.5001)then ! Still ok
                 fZC  = 0.999*fZC
         else
           call set_error('Strange mass centre, ix=' + fu_str(ix),'vert_adv_column_pass')
           return
         endif 
       endif
     endif
     
     bounds(1) = xr_advected(ix-1)
     bounds(3) = xr_advected(ix)

     ! Have to get fraction of each slab and position of their boundary

     if (fZC < 0.) then ! on the left side of the cell
         if (ix==nCells) then
              bgmass = bgdensr * (bounds(3)-bounds(1))
         else
              bgmass = PassengersIn(0,ix+1) * (1. + 2*cmRelIn(ix+1)) * (bounds(3)-bounds(1)) / (xr_advected(ix+1)-bounds(3))
              !bgmass = PassengersIn(0,ix+1)  * (bounds(3)-bounds(1)) / (xr_advected(ix+1)-bounds(3))
         endif

         if (bgmass == 0.) then  !No mass there -- just rectangle with zero bg 
            ftmp =  1.D0 +  2*fZC
            bounds(2) =  ftmp*bounds(3) +  (1.0-ftmp)*bounds(1)  
            fracmass(1) = 1.D0
            ! bl_abs =  xr_advected(ix-1)
            ! br_abs = ftmp*xr_advected(ix) + (1.0-ftmp)*bl_abs
         else
            bgfrac = bgmass/PassengersIn(0,ix)
            exZc = fZC/(1.0D0 - bgfrac) !Mass-centre of over-background part
            if (exZc > -0.5 .and. exZC < 0.) then ! BG consistent with CM
               ftmp =  1.D0 +  2*exZC
               bounds(2) =  ftmp*bounds(3) +  (1.0-ftmp)*bounds(1)
               fracmass(1)  = 1.D0 - bgfrac*(1D0-fTmp)
            else
               bounds(2) = bounds(1) !Maximal fill
               fracmass(1) = -2*fZC  !zero-sized slab on left
            endif
        endif
     else
         if (ix==1) then
            bgmass = bgdensl * (bounds(3)-bounds(1))
         else
            bgmass = PassengersIn(0,ix-1) *(1. + 2*cmRelIn(ix - 1)) *  (bounds(3)-bounds(1)) / (bounds(1)-xr_advected(ix-2))
            !bgmass = PassengersIn(0,ix-1)  *  (bounds(3)-bounds(1)) / (bounds(1)-xr_advected(ix-2))
         endif

         if (bgmass == 0.) then  !No mass there -- just rectangle with zero bg 
            ftmp =  1.D0 -  2*fZC
            bounds(2) =  ftmp*bounds(1) +  (1.0-ftmp)*bounds(3)  
            fracmass(1) = 0.
            fracmass(2) = 1.
            ! br_abs = xr_advected(ix)
            !bl_abs = ftmp * xr_advected(ix-1) + (1.0-ftmp) * br_abs
         else
            bgfrac = bgmass/PassengersIn(0,ix)
            exZc = fZC/(1.0D0 - bgfrac) !Mass-centre of over-background part
            if (exZc < 0.5 .and. exZC > 0.) then ! BG consistent with CM
               ftmp =  1.D0 -  2*exZC
               bounds(2) =  ftmp*bounds(1) +  (1.0-ftmp)*bounds(3)
               fracmass(1) = bgfrac*(1D0-fTmp)
            else
               bounds(2) = bounds(3) !Maximal fill
               fracmass(1) = 1.0D0 - 2*fZC  !zero-sized slab on left
            endif
         endif
     end if
     fracmass(2) = 1.0D0 - fracmass(1)
     CHECK(.not. all(fracmass >= 0), "non-positive fraction")

     do iSubSlab= 1,2
        bl_abs = bounds(iSubSlab)
        br_abs = bounds(iSubSlab+1)
        if (fracmass(iSubSlab) < 1e-20) cycle
        SS = (BR_abs - BL_abs) / fracmass(iSubSlab) 
        CHECK(.not. SS>=0, "non-positive slab")

        !Still within the domain?
        if (BL_abs + loopx >= xr(nCells))  then 
           if (ifLoop) then
             loopx = loopx - xr(nCells)
             ixto = 1
           else
              !Whole slab and following slabs go above the domain
               ixTo = nCells+1
               passengersOut(0:nPass,ixTo) = passengersOut(0:nPass,ixTo) + &
                                 &  passengersIn(0:nPass, ix) * fracmass(iSubSlab)
               REPORT (7,  fracmass(iSubSlab) )
               cycle 
           endif
        endif
        
        BL_tmp = BL_abs + loopx

        do  ixto = ixto, nCells
           if (BL_tmp < xr(ixto)) exit
        enddo

        !now BL_tmp is within ixto cell

        do while (BR_abs + loopx > xr(ixto)) 
              ! Except for rightmost part of advected slab
              BR_tmp = xr(ixto)
              ftmp =  (BR_tmp - BL_tmp) /SS  !Fraction that goes here
              CHECKEXIT( .not. ftmp>=0 , "non-positive slab1")
              passengersOut(0:nPass,ixTo) = passengersOut(0:nPass,ixTo) + &
                  &  passengersIn(0:nPass, ix) * ftmp
   !             call msg("3:ix, frac", ix, fTmp)
              REPORT (3, fTmp )
              if (ixto > 0) then
                  momOut(iXto) = momOut(iXto) +  passengersIn(0, ix) * ftmp * &
                         0.5* (1. - (BR_tmp - BL_tmp)/dx(ixto)) !right-justified slab
                  CHECKEXIT(.not. (abs(momOut(iXto))<= 0.5*passengersOut(0,ixTo)), "Gotcha1")
              endif
              BL_tmp = xr(ixto)
              ixto = ixto + 1
              if (ixto > nCells) then 
                 if(.not. ifLoop)  exit
                 loopx = loopx - xr(nCells)   !handle loopover
                 BL_tmp = 0.D0 !  xr(0)  
                 ixto = 1
              endif
         enddo
#ifdef DEBUG_V5
         if (error) exit
#endif
         !Leftover
         BR_tmp = BR_abs + loopx
         if (SS > 0.) then
            ftmp =  (BR_tmp - BL_tmp) /SS  !Fraction that goes here
         else
            ftmp = fracmass(iSubSlab)
         endif
         passengersOut(0:nPass,ixTo) = passengersOut(0:nPass,ixTo) + &
                                    & passengersIn(0:nPass, ix) * ftmp
   !       call msg("4:ix, frac", ix, fTmp)
         REPORT (4, fTmp )
         if (ixto <= nCells .and.  ixto> 0) then  !non-justified slab
             fTmp =  fTmp * 0.5D0* ( BR_tmp + BL_tmp - xr(ixto) - xr(ixto-1))
             fTmp =  fTmp / dx(ixto) !!CM * massfraction
             momOut(iXto) = momOut(iXto) + passengersIn(0, ix) * ftmp 
             CHECK(.not.(abs(momOut(iXto)) <= 0.5*passengersOut(0,ixTo)), "Gotcha2")
         endif
      enddo !iSubSlab
    enddo IXLOOP !ix

    !CHECK( (.not. ifLoop) .or. ix == iCellEnd+1, "No leftovers allowed in looped distribution")
!     fMass = sum(PassengersIn(0,0:ix-1))
!     fMass1 = sum(PassengersOut(0,0:nCells+1))
!     print *, "ix=", ix, fMass, fMass1, abs(fMass-fMass1)/fMass

    if ( PassengersIn(0,nCells+1) /= 0. &
#ifdef DEBUG_V5
               & .and. .not. error  &
#endif               
                & ) then !get stuff from right boundary
       ! PassengersIn(0,nCells+1) should be non-zero only if something to be injected
       !out_of_domain source cell, inflow, indomain target cell
       SS  = xr(nCells) - xr_advected(nCells) ! Should be positive!!!
       CHECK( .not.(ix == nCells+1 .and. SS > 0.  .and. ixto <= nCells), "Can't inject right boundary")

       do  ixto = ixto, nCells
          if (xr_advected(nCells) < xr(ixto)) exit
       enddo
       
       fTmp = xr(ixto) - xr_advected(nCells)  !first, incomplete slab from boundary
       PassengersOut(:,ixto) =  PassengersOut(:,ixto) + &
              & fTmp / SS * PassengersIn(:,nCells+1)
!        call msg("5:ix, frac", ix, fTmp / SS )
       REPORT (5, fTmp/SS )
       momOut(ixto) = momOut(ixto) + 0.5* (1. - fTmp/dx(ixto)) * & !right-justified slab
            &  fTmp / SS * PassengersIn(0,nCells+1)

       CHECK(abs(momOut(iXto)) > 0.5*passengersOut(0,ixTo), "Gotcha3")
       do ixto = ixto+1,nCells 
         PassengersOut(:,ixto) = dx(ixto)/SS * PassengersIn(:,nCells+1)
!         call msg("6:ix, frac", ix, SS/dx(ixto))
         REPORT (6, 1.0 )
       enddo
    endif


    if (ifCmOut) then 
     ! Turn  Moments to CM within domain
      do ix = 1, nCells 
        fmass = passengersOut(0,ix) 
        if (abs(fmass) >  fMinAdvectedMass) then
          momOut(ix) = momOut(ix) /  fmass  ! CM: +-0.5 , positive -- Up
          if(abs(momOut(ix)) > 0.4999)then 
            if(have_negatives)then      ! Crazy centres are legal...
              momOut(ix) = 0.0
            elseif (abs(momOut(ix)) < 0.5001)then 
               momOut(ix) = sign(0.499,momOut(ix))
            else  !Too bad
               call set_error("Wrong CM at advection output, ix="+fu_str(ix), "advect_mass")
            endif
          endif
        else
          momOut(ix) = 0.0
        endif
      enddo
#ifdef DEBUG_V5 
   else
     ! Turn  Check Moment
      do ix = 1, nCells 
        fmass = passengersOut(0,ix) 
        if (abs(fmass) >  fMinAdvectedMass) then
          fTmp = abs(momOut(ix) /  fmass)  ! CM: +-0.5 , positive -- Up
          if(.not. (fTmp <= 0.5) )then
               call msg("fmass, momOut(ix)",fmass, momOut(ix))
            call set_error("Wrong Moment at advection output, ix="+fu_str(ix), "advect_mass")
          endif
         endif
       enddo
#endif
   endif

#ifdef DEBUG_V5
!    fMass = sum(passengersIn(0,0:nCells+1))
!    fMass1 = sum(passengersOut(0,0:nCells+1))
    fMass = sum(passengersIn(0,iCellStart:iCellEnd))
    fMass1 = sum(passengersOut(0,0:nCells+1))
    if (.not. all(passengersOut(0,0:nCells+1)>=0.) .or. error .or. &
        & (abs(fMass - fMass1) > 1e-7*ncells*(fMass+fMinAdvectedMass)) )then
      !$OMP CRITICAL(v4bark)
       call msg("Trouble after mass advection!" )
       call msg("Line mass in, out", fMass, fMass1)
       call msg("Line mass difference", fmass-fmass1)
       call msg("Mass, in  (0:nCells+1)", passengersIn(0,0:nCells+1))
       call msg("Mass, out (0:nCells+1)", passengersOut(0,0:nCells+1))
       call msg("xr,         (0:nCells):", xr(0:nCells))
       call msg("xr_advected (0:nCells):", xr_advected(0:nCells))
       call msg("cm, in (1:nCells)",   cmRelIn(1:nCells))
       if (ifCmOut) then
          call msg("cm, out (1:nCells)", momOut(1:nCells))
       else
          call msg("mom, out (1:nCells)", momOut(1:nCells))
       endif
       call set_error("Trouble after mass advection!","advect_mass")
      !$OMP END CRITICAL(v4bark)
    endif
#endif
#undef REPORT
#undef CHECK

  end subroutine advect_mass_step





   !*****************************************************************************************

  subroutine advect_mass_trislab(   cmRelIn,  momOut, &
                          & PassengersIn,  PassengersOut, nPass, &
                          & xr, xr_advected, dx,& 
                          & ifLoop, iCellStart, iCellEnd, nCells,&
                          & fMinAdvectedMass, have_negatives, ifCMout, cmrelax)
    ! 
    ! Redistributes masses and moments according to  rightbound_offset.
    ! Advects whatever non-zero mass it finds. Garbage should be collected elsewhere
    ! Warning!! ifLoop loops around line indices as   ..., nCells-1, nCells,  1,  2, ...
    ! ifCmout -- output CM instead of moments
    ! Roux

    real, dimension(1:),    intent(in)  :: cmRelIn ! 1:nCells -- domain only
    real, dimension(1:),    intent(out) :: momOut  ! 1:nCells -- domain only
    real, dimension(1:),    intent(in)  :: dx !(1:nCells) !Cell size (e.g. in kg)
    real(r8k), dimension(0:),    intent(in)  :: xr, xr_advected !(0:nCells)  !same units as above!!!
    real, dimension(0:,0:), intent(in)  :: PassengersIn  !(0:npass, 0:nCells+1) 
    real, dimension(0:,0:), intent(out) :: PassengersOut !(0:npass, 0:nCells+1)
    real, intent(in) :: cmrelax ! Facttor to introduce diffusion bu CM relaxation

#ifdef DEBUG_V5
#define CHECK(what, bark) if (what) call set_error(bark,"advect_mass_trislab"); if (error) return 
! Print full redistribution info: iFrom, iTo , fraction    
!!!#define REPORT(II,frac) if (ifReport) call msg("II: ix="+fu_str(ix), ixto, frac)
#define REPORT(II,frac)
#else
#define CHECK(what, bark)
#define REPORT(II,frac)
#endif

    ! nlev == nz_dispersion arrays should have one more element
    integer, intent(in) :: iCellStart, iCellEnd, nCells, nPass
    real, intent(in) ::  fMinAdvectedMass ! Used only for conversion to centers of mass
    logical, intent(in) :: have_negatives, ifLoop, ifCMout

    real(r8k) ::  fZC, SS, BR_abs, BL_abs, BL_tmp, BR_tmp, ftmp,  loopx,  SS1, alpha, MMRL
    real(r8k), parameter :: fZcTrimax = 1.0_r8k/6  ! Maximum CM for trapezoid slab
    real :: fmass, fMass1
    integer ::  ix, ixto



    ! !No advection hack
    if (.False.) then
       PassengersOut(:,0:nCells+1)=PassengersIn(:,0:nCells+1)
       if (ifCMout) then
         momOut(1:nCells)= cmRelIn(1:nCells)
       else
         momOut(1:nCells)= cmRelIn(1:nCells)*PassengersIn(0,1:nCells)
       endif
       return
    endif

    momOut(1:nCells)=0.0
    passengersOut(:,0:nCells+1) = 0
    
    ! find first advectable mass
    do ix = iCellStart, iCellEnd
      if (PassengersIn(0,ix) == 0.) cycle
      exit
    enddo
    

    !Init distribution loop
    loopx = 0.
    ixto = 0
    if (ix == 0) then !Some left boundary stuff to distribute
       SS  = xr_advected(0) ! Should be positive!!!
       CHECK(.not. (SS > 0), "Can't inject left boundary")
       ixto = 1
       do while (SS > xr(ixto)) 
        PassengersOut(:,ixto) = dx(ixto)/SS * PassengersIn(:,0)
       REPORT (1,dx(ixto)/SS)
        ixto = ixto + 1 ! Moments are zero anyhow
      enddo
      fTmp = xr_advected(0)-xr(ixto-1)  !size of leftover
      PassengersOut(:,ixto) =  fTmp / SS * PassengersIn(:,0) !Remaining part
      REPORT (2, fTmp / SS)
      momOut(ixto) =  - 0.5* (1. - fTmp/dx(ixto)) * PassengersOut(0,ixto)  !left-justified slab
      ix = 1
   elseif (xr_advected(ix-1) < 0) then !First mass below the domain
      if (ifLoop) then 
         loopx =  xr(nCells) !Add this part to advected quantities
         do ixto = nCells,1,-1   !Find ixTo for cell ix 
           if (xr_advected(ix-1) + loopx > xr(ixto-1)) exit
         enddo
      endif
   endif


   ! Main distribution loop
   do ix = ix, min(iCellEnd,nCells)
     if (PassengersIn(0,ix) == 0.) cycle

     fZC = cmRelIn(ix) * cmrelax

     if (abs(fZC) >= 0.5 ) then  ! abnormal value check
       if (have_negatives) then 
         call msg('Strange mass centre, return to cell:', fZC)
         fZC = 0.
       else
         call msg('Strange mass centre:', fZC)
         if (abs(fZC) < 0.5001)then ! Still ok
                 fZC  = 0.999*fZC
         else
           call set_error('Strange mass centre, ix=' + fu_str(ix),'advect_mass_trislab')
           return
         endif 
       endif
     endif
     
     !MMR = MMRL/SS + alpha*(ax-bl_abs)/(SS*SS)
     !MMR(middle_slab) = 1./SS
     if (abs(fZC) < fZcTrimax) then !Trapezoid slab
         ! on the lower side of the layer
         bl_abs =  xr_advected(ix-1)
         br_abs = xr_advected(ix)
         SS = br_abs - bl_abs
         MMRL = (1.0_r8k - 6 * fzc)   ! in fractions of slab mass
         alpha = 12 * fZC 
     else !Triangular slab
        ftmp =  1.0_r8k -  2*abs(fZC) + 2*fZcTrimax ! rel. size of th
        ftmp = max(ftmp, 1e-6)
        if (fZC < 0) then !Triange slab left
           bl_abs =  xr_advected(ix-1)
           br_abs = ftmp*xr_advected(ix) + (1.0-ftmp)*bl_abs
           MMRL =  2.
           alpha = -2.
        else
           br_abs = xr_advected(ix)
           bl_abs = ftmp * xr_advected(ix-1) + (1.0-ftmp) * br_abs
           MMRL =  0.
           alpha = 2.
        end if
        !now BL is just above iLevTo
        SS = BR_abs - BL_abs 

     endif

     !Still within the domain?
     if (BL_abs + loopx >= xr(nCells))  then 
        if (ifLoop) then
          loopx = loopx - xr(nCells)
          ixto = 1
        else
           !Whole slab and following slabs go above the domain
            ixTo = nCells+1
            passengersOut(0:nPass,ixTo) = passengersOut(0:nPass,ixTo) + &
                              &  passengersIn(0:nPass, ix) 
            REPORT (7, fTmp )
            cycle
        endif
     endif
     
     BL_tmp = BL_abs + loopx
     !
     ! Where do we go?
     !
     do  ixto = ixto, nCells
        if (BL_tmp < xr(ixto)) exit
     enddo
     !
     ! If slab is too thin, put it to ixTo and cycle, beware negative one
     !
     if(.not. SS > 1d-5 * dx(ix))then

       if(.not. SS >= 0.0)then
          !$OMP CRITICAL(nonpositive)
          call msg('non-positive slab, advect_mass_trislab, ix=',ix, nCells)
          call msg('xr(ix-2:ix+2)', (/xr(max(0,ix-2):min(ix+2,size(xr)))/))
          call msg('xradvected(ix-2:ix+2)', (/xr_advected(max(0,ix-2):min(ix+2,size(xr_advected)))/))
          call msg('ftmp,bl_abs, br_abs,ss,alpha,MMRL', (/fZC,ftmp,bl_abs, br_abs,ss,alpha,MMRL/))
          call set_error("non-positive slab",'advect_mass_trislab')
          !$OMP END CRITICAL(nonpositive)
       endif
       passengersOut(0:nPass,ixTo) = passengersOut(0:nPass,ixTo) + &
                                   &  passengersIn(0:nPass, ix)
       momOut(iXto) = momOut(iXto) + passengersIn(0, ix) * &
                    & 0.5 * (BR_abs+BL_abs - xr(ixTo-1)-xr(ixTo)) / dx(ixto)
       CHECK(abs(momOut(iXto)) > 0.5*passengersOut(0,ixTo), "Gotcha thin-slab")
       cycle
     endif
     
     do while (BR_abs + loopx > xr(ixto)) 
           ! Except for rightmost part of advected slab
           BR_tmp = xr(ixto)
           SS1 = (BR_tmp - BL_tmp) / SS
           ftmp =  ( MMRL + 0.5* alpha*SS1 ) * SS1   !Fraction that goes here

           CHECK( .not. ftmp>=0 , "non-positive slab1")
           passengersOut(0:nPass,ixTo) = passengersOut(0:nPass,ixTo) + &
               &  passengersIn(0:nPass, ix) * ftmp
!             call msg("3:ix, frac", ix, fTmp)
           REPORT (3, fTmp )
           if (ixto > 0) then
               !momOut(iXto) = momOut(iXto) +  passengersIn(0, ix) * ftmp * &
               !       0.5* (1. - (BR_tmp - BL_tmp)/dx(ixto)) !right-justified slab
!               0.5_r8k* ( BR_tmp + BL_tmp - xr(ixto) - xr(ixto-1))
               fTmp = (BL_tmp - 0.5_r8k*(xr(ixto) + xr(ixto-1)))/SS ! BL Offset from the center
               momOut(iXto) =   momOut(iXto) + passengersIn(0, ix) * &
                  & SS1*(SS1*(SS1*alpha/3.0_r8k + 0.5*MMRL + 0.5*alpha*fTmp) + MMRL*fTmp)*SS/dx(ixto)

               CHECK(abs(momOut(iXto)) > 0.5*passengersOut(0,ixTo), "Gotcha1")
!               if(abs(momOut(iXto)) > 0.5*passengersOut(0,ixTo))then
!                 !$OMP CRITICAL(g1)
!                 call msg_warning("Gotcha1 at:" + fu_str(iXto), 'advect_mass_trislab')
!                 call msg('momOut, passengersOut:', momOut(iXto), passengersOut(0,ixTo))
!                 call set_error("Gotcha1", 'advect_mass_trislab')
!                 !$OMP END CRITICAL(g1)
!                 return
!               endif
           endif

           MMRL = MMRL + alpha*SS1
           BL_tmp = xr(ixto)
           ixto = ixto + 1
           if (ixto > nCells) then 
              if(.not. ifLoop)  exit
              loopx = loopx - xr(nCells)   !handle loopover
              BL_tmp = 0.0_r8k !  xr(0)  
              ixto = 1
           endif
      enddo
      !Leftover
      BR_tmp = BR_abs + loopx
      if (SS == 0) cycle

      SS1 = (BR_tmp - BL_tmp) / SS
      ftmp =  ( MMRL + 0.5* alpha*SS1 ) * SS1 !Fraction that goes here
      
      passengersOut(0:nPass,ixTo) = passengersOut(0:nPass,ixTo) + &
                                  & passengersIn(0:nPass, ix) * ftmp
      REPORT (4, fTmp )
      if (ixto <= nCells .and.  ixto> 0) then  !non-justified slab
!          fTmp =  fTmp * 0.5_r8k* ( BR_tmp + BL_tmp - xr(ixto) - xr(ixto-1))
!          fTmp = fTmp / dx(ixto) !!CM * massfraction
!          momOut(iXto) = momOut(iXto) + passengersIn(0, ix) * ftmp 

           fTmp = (BL_tmp - 0.5_r8k*(xr(ixto) + xr(ixto-1)))/SS ! BL Offset from the center
           momOut(iXto) =   momOut(iXto) + passengersIn(0, ix) * &
               & SS1*(SS1*(SS1*alpha/3.0_r8k + 0.5_r8k*(MMRL + alpha*fTmp)) + MMRL*fTmp)*SS/dx(ixto)
          CHECK(abs(momOut(iXto)) > 0.5*passengersOut(0,ixTo), "Gotcha2")
!          if(abs(momOut(iXto)) > 0.5*passengersOut(0,ixTo))then
!            !$OMP CRITICAL(g1)
!            call msg_warning("Gotcha2 at:" + fu_str(iXto), 'advect_mass_trislab')
!            call msg('momOut, passengersOut:', momOut(iXto), passengersOut(0,ixTo))
!            call set_error("Gotcha2", 'advect_mass_trislab')
!            !$OMP END CRITICAL(g1)
!            return
!          endif
       endif

    enddo !ix

!    CHECK( (.not. ifLoop) .or. ix == iCellEnd+1, "No leftovers allowed in looped distribution")

!    if( abs(PassengersIn(0,nCells+1)) > fMinAdvectedMass) then !get stuff from right boundary
    if ( PassengersIn(0,nCells+1) /= 0.) then !get stuff from right boundary
       ! PassengersIn(0,nCells+1) should be non-zero only if something to be injected
       !out_of_domain source cell, inflow, indomain target cell
       SS  = xr(nCells) - xr_advected(nCells) ! Should be positive!!!
       CHECK( .not.(ix == nCells+1 .and. SS > 0.  .and. ixto <= nCells), "Can't inject right boundary")

       do  ixto = ixto, nCells
          if (xr_advected(nCells) < xr(ixto)) exit
       enddo
       
       fTmp = xr(ixto) - xr_advected(nCells)  !first, incomplete slab from boundary
       PassengersOut(:,ixto) =  PassengersOut(:,ixto) + &
              & fTmp / SS * PassengersIn(:,nCells+1)
!        call msg("5:ix, frac", ix, fTmp / SS )
       REPORT (5, fTmp/SS )
       momOut(ixto) = momOut(ixto) + 0.5* (1. - fTmp/dx(ixto)) * & !right-justified slab
            &  fTmp / SS * PassengersIn(0,nCells+1)

       CHECK(abs(momOut(iXto)) > 0.5*passengersOut(0,ixTo), "Gotcha3")
!       if(abs(momOut(iXto)) > 0.5*passengersOut(0,ixTo))then
!         !$OMP CRITICAL(g1)
!         call msg_warning("Gotcha3 at:" + fu_str(iXto), 'advect_mass_trislab')
!         call msg('momOut, passengersOut:', momOut(iXto), passengersOut(0,ixTo))
!         call set_error("Gotcha3", 'advect_mass_trislab')
!         !$OMP END CRITICAL(g1)
!         return
!       endif
       do ixto = ixto+1,nCells 
         PassengersOut(:,ixto) = dx(ixto)/SS * PassengersIn(:,nCells+1)
!         call msg("6:ix, frac", ix, SS/dx(ixto))
         REPORT (6, 1.0 )
       enddo
    endif


    if (ifCmOut) then 
     ! Turn  Moments to CM within domain
      do ix = 1, nCells 
        fmass = passengersOut(0,ix) 
        if (abs(fmass) >  fMinAdvectedMass) then
          momOut(ix) = momOut(ix) /  fmass  ! CM: +-0.5 , positive -- Up
          if(abs(momOut(ix)) > 0.4999)then 
            if(have_negatives)then      ! Crazy centres are legal...
              momOut(ix) = 0.0
            elseif (abs(momOut(ix)) < 0.5001)then 
               momOut(ix) = sign(0.499,momOut(ix))
            else  !Too bad
               call set_error("Wrong CM at advection output, ix="+fu_str(ix), "advect_mass_trislab")
            endif
          endif
        else
          momOut(ix) = 0.0
        endif
      enddo
#ifdef DEBUG_V5 
   else
     ! Turn  Check Moment
      do ix = 1, nCells 
        fmass = passengersOut(0,ix) 
        if (abs(fmass) >  fMinAdvectedMass) then
          fTmp = abs(momOut(ix) /  fmass)  ! CM: +-0.5 , positive -- Up
          if(.not. (fTmp <= 0.5) )then
               call msg("fmass, momOut(ix)",fmass, momOut(ix))
            call set_error("Wrong Moment at advection output, ix="+fu_str(ix), "advect_mass_trislab")
          endif
         endif
       enddo
#endif
   endif

#ifdef DEBUG_V5
!    fMass = sum(passengersIn(0,0:nCells+1))
!    fMass1 = sum(passengersOut(0,0:nCells+1))
    fMass = sum(passengersIn(0,iCellStart:iCellEnd))
    fMass1 = sum(passengersOut(0,0:nCells+1))
    if (.not. all(passengersOut(0,0:nCells+1)>=0.) .or. error .or. &
        & (abs(fMass - fMass1) > 1e-7*ncells*(fMass+fMinAdvectedMass)) )then
      !$OMP CRITICAL(v4bark)
       call msg("Trouble after mass advection2!" )
       call msg("Line mass in, out", fMass, fMass1)
       call msg("Line mass difference", fmass-fmass1)
       call msg("Mass, in  (0:nCells+1)", passengersIn(0,0:nCells+1))
       call msg("Mass, out (0:nCells+1)", passengersOut(0,0:nCells+1))
       call msg("xr,         (0:nCells):", xr(0:nCells))
       call msg("xr_advected (0:nCells):", xr_advected(0:nCells))
       call msg("cm, in (1:nCells)",   cmRelIn(1:nCells))
       if (ifCmOut) then
          call msg("cm, out (1:nCells)", momOut(1:nCells))
       else
          call msg("mom, out (1:nCells)", momOut(1:nCells))
       endif
       call set_error("Trouble after mass advection!","advect_mass_trislab")
      !$OMP END CRITICAL(v4bark)
    endif
#endif
#undef REPORT
#undef CHECK

  end subroutine advect_mass_trislab
   !*****************************************************************************************

  subroutine MASS_DISTRIBUTOR(   cmRelIn,  momOut, &
                          & PassengersIn,  PassengersOut, nPass, &
                          & xr, xr_advected, dx,& 
                          & ifLoop, iCellStart, iCellEnd, nCells,&
                          & fMinAdvectedMass, have_negatives, ifCMout, cmrelax)

    ! Dispatcher for mass distributors
    ! Roux

    real, dimension(1:),    intent(in)  :: cmRelIn ! 1:nCells -- domain only
    real, dimension(1:),    intent(out) :: momOut  ! 1:nCells -- domain only
    real, dimension(1:),    intent(in)  :: dx !(1:nCells) !Cell size (e.g. in kg)
    real*8, dimension(0:),    intent(in)  :: xr, xr_advected !(0:nCells)  !same units as above!!!
    real, dimension(0:,0:), intent(in)  :: PassengersIn  !(0:npass, 0:nCells+1) 
    real, dimension(0:,0:), intent(out) :: PassengersOut !(0:npass, 0:nCells+1)
    real, intent(in) :: cmrelax ! Facttor to introduce diffusion bu CM relaxation
    integer, intent(in) :: iCellStart, iCellEnd, nCells, nPass
    real, intent(in) ::  fMinAdvectedMass ! Used only for conversion to centers of mass
    logical, intent(in) :: have_negatives, ifLoop, ifCMout

    select case (distributor_type)
      case (advect_tri)
         call advect_mass_trislab(   cmRelIn,  momOut, &
                          & PassengersIn,  PassengersOut, nPass, &
                          & xr, xr_advected, dx,&
                          & ifLoop, iCellStart, iCellEnd, nCells,&
                          & fMinAdvectedMass, have_negatives, ifCMout, cmrelax)

      case (advect_rect)
         call advect_mass(   cmRelIn,  momOut, &
                          & PassengersIn,  PassengersOut, nPass, &
                          & xr, xr_advected, dx,&
                          & ifLoop, iCellStart, iCellEnd, nCells,&
                          & fMinAdvectedMass, have_negatives, ifCMout, cmrelax)
      case (advect_step)
         call advect_mass_step(   cmRelIn,  momOut, &
                          & PassengersIn,  PassengersOut, nPass, &
                          & xr, xr_advected, dx,&
                          & ifLoop, iCellStart, iCellEnd, nCells,&
                          & fMinAdvectedMass, have_negatives, ifCMout, cmrelax)
      case default
         call set_error("Strange mass distributor type", "MASS_DISTRIBUTOR")


    end select


  end subroutine MASS_DISTRIBUTOR



 !**************************************************************************
  subroutine get_linestuff_from_boundaries(pdispFlds, moment_my, ii2, ii3, nCells, &
         & left_boundary, right_boundary, iCellStartAdv, iCellEndAdv, &
         & inflow_left,  inflow_right, airdens_left,  airdens_right, &
         & mystuff, nSrc, nSp, arMinAdvMass, pBBuf,  massIn0,  massInM)

         ! Puts mass_map stuff from boundaries.
         ! Fixme!!! Does not bring _TRUE_ passengers (that are not moments) 

      implicit none

      !Parameters different between axes 
      type(Tmass_map), intent(in) ::  pdispFlds,moment_my
      integer, intent(in) :: ii2, ii3 ! indices in other directions in boundaries
      integer, intent(in) :: nCells  !Boundaries come to cells  0 and nCells+1 
      integer, intent(in) :: left_boundary, right_boundary ! e.g. east_bondary
      integer, intent(out) :: iCellStartAdv, iCellEndAdv ! assigned  0 and nCells+1 
                                                         ! if something comes

      real, intent(in) :: inflow_left,  inflow_right, & ! kg/s incoming positive!!!
                & airdens_left,  airdens_right ! for boundary to (kg|mole|bq)/kg_air conversion
      real(r8k), dimension(:,:), intent(inout) :: massIn0,  massInM ! (1:nSrc 1:nSp) Slice of  x0_mass, xM_mass, 
                                                   !y0_mass, yM_mass  bottom_mass, or top_mass
                                                   ! Accumulates injected mass

      ! Parameters same for all axes 
      type (T_Galp_v5_thread_stuff), intent(inout) :: mystuff
      integer, intent(in) :: nSrc, nSp
      real, dimension (:), intent(in) :: arMinAdvMass
      type(TboundaryBuffer), pointer, intent(in) :: pBBuf
            
      !Local
      real :: fMass, fTmp
      integer :: iSrc, iSpecies, iTmp, nPass
      real, parameter :: lmtfactor = 1.
                                 !!Dirty hack to fight retared coding of
                                !   boundary files that causes negative
                                !   concentrations
                             !Should we put boundaries crap to our nice garbage?



       mystuff%passengers(:,0,:,:) = 0.
       mystuff%passengers(:,nCells+1,:,:) = 0.

       if (inflow_left > 0) then
         select case (pBBuf%iBoundaryType(left_boundary))
           case(zero_boundary_type, smpi_comm_boundary_type)
           case(polar_boundary_type)
             ! left pole is always South, species always mapped 1:1, nothing is skipped
             do iSrc = 1, nSrc   ! emission sources
               do ispecies = 1, pdispFlds%nspecies
                  nPass = moment_my%passengers(iSpecies)%nPassengers
                  fMass = pBBuf%bSouth(iSpecies, iSrc, 1, ii3) /   & ! Mass here!
                               pBBuf%PoleCapAirmass(left_boundary)%pp(ii3) * inflow_left
                  fMass = fMass * pBBuf%outflowFactor(southern_boundary,ii3)
                   iCellStartAdv = 0
                   mystuff%passengers(0,0,iSpecies,iSrc) = fMass         ! mass
                   mystuff%passengers(nPass+2,0,iSpecies,iSrc) = &       ! z moment
                             & fMass * pBBuf%bSouth(iSpecies, iSrc, 2, ii3)
                   massIn0(iSrc,iSpecies) = massIn0(iSrc,iSpecies) + fMass
                enddo
             enddo

           case(dirichlet_boundary_type)
                  iSrc = 1
                  fTmp = 1. / airdens_left * inflow_left ! Boundary to injected_mass factor
                  do ispecies = 1, pdispFlds%nspecies
                    iTmp = pBBuf%iBoundarySpecies(iSpecies, left_boundary)
                    if (iTmp == int_missing) cycle
                    nPass = moment_my%passengers(iSpecies)%nPassengers
                    fMass = pBBuf%pBoundaryData(left_boundary)%ptr(iTmp, iSrc, ii2, ii3) * fTmp
                    if (fMass>0) then 
                      mystuff%passengers(0,0,iSpecies,iSrc) = fMass
                      massIn0(iSrc,iSpecies) = massIn0(iSrc,iSpecies) + fMass
                      iCellStartAdv = 0
                    elseif (fMass >  -lmtfactor*arMinAdvMass(ispecies)) then ! Some numrics, still okay
                      mystuff%tmp_garbage(ispecies, iSrc) =  mystuff%tmp_garbage(ispecies,isrc) + fMass
                    else
                      call msg("iSpecies,iSrc", iSpecies,iSrc)
                      call msg("indices  ii2, ii3", ii2, ii3)
                      call msg("Concentration, iMass", &
                                & pBBuf%pBoundaryData(left_boundary)%ptr(iTmp, iSrc, ii2, ii3), fMass)
                      call set_error("Negative concentration in _"+BoundaryChar(left_boundary)+"_ boundary", &
                                                &"get_linestuff_from_boundaries")
                    endif
                  enddo
             case default
                call set_error("Strange _"+BoundaryChar(left_boundary)+"_ boundary type:"+ &
                   &fu_str(pBBuf%iBoundaryType(left_boundary)), "get_linestuff_from_boundaries")
             endselect
       endif !Left boundary

       if (inflow_right > 0) then
         select case (pBBuf%iBoundaryType(right_boundary))
           case(zero_boundary_type, smpi_comm_boundary_type)
           case(polar_boundary_type)
             ! right pole is always North, species always mapped 1:1, nothing is skipped
             do iSrc = 1, nSrc   ! emission sources
               do ispecies = 1, pdispFlds%nspecies
                  nPass = moment_my%passengers(iSpecies)%nPassengers
                  fMass = pBBuf%bNorth(iSpecies, iSrc, 1, ii3) /   &! Mass here!
                               pBBuf%PoleCapAirmass(right_boundary)%pp(ii3) * inflow_right
                  fMass = fMass * pBBuf%outflowFactor(northern_boundary,ii3)
                  iCellEndAdv = nCells+1
                  mystuff%passengers(0,nCells+1,iSpecies,iSrc) = fMass         ! mass
                  mystuff%passengers(nPass+2,nCells+1,iSpecies,iSrc) = &       ! z moment
                            & fMass * pBBuf%bNorth(iSpecies, iSrc, 2, ii3)
                  massInM(iSrc,iSpecies) = massInM(iSrc,iSpecies) + fMass
                enddo
             enddo
    !         if (inflow_left > 0. ) call msg("MMR * 1e6 from NP, SP", &
    !                & mystuff%passengers(0,nCells+1,1,1)/inflow_right * 1e6, &
    !                & mystuff%passengers(0,0,1,1)/inflow_left * 1e6)

           case(dirichlet_boundary_type)
                  iSrc = 1
                  if (airdens_right < 0.01) call set_error("Gotchasdafaqsdf","fdgsdf")
                  fTmp = 1./airdens_right * inflow_right !Boundary to injected_mass factor 
                  do ispecies = 1, pdispFlds%nspecies
                    iTmp = pBBuf%iBoundarySpecies(iSpecies, right_boundary)
                    if (iTmp == int_missing) cycle
                    nPass = moment_my%passengers(iSpecies)%nPassengers
                    fMass = pBBuf%pBoundaryData(right_boundary)%ptr(iTmp, iSrc, ii2, ii3) *fTmp 
                    if (fMass>0) then 
                      iCellEndAdv = nCells+1
                      mystuff%passengers(0,nCells+1,iSpecies,iSrc) = fMass
                      massInM(iSrc,iSpecies) = massInM(iSrc,iSpecies) + fMass
                    elseif (fMass >  -lmtfactor*arMinAdvMass(ispecies)) then ! Some numrics, still okay
                      mystuff%tmp_garbage(ispecies, iSrc) =  mystuff%tmp_garbage(ispecies,isrc) + fMass
                    else
                      call msg("iSpecies,iSrc", iSpecies,iSrc)
                      call msg("indices  ii2, ii3", ii2, ii3)
                      call msg("Concentration, iMass", &
                                & pBBuf%pBoundaryData(right_boundary)%ptr(iTmp, iSrc, ii2, ii3), fMass)
                      call set_error("Negative concentration in _"+BoundaryChar(right_boundary)+"_ boundary",&
                                                &"get_linestuff_from_boundaries")
                      
                    endif
                  enddo
             case default
                call set_error("Strange _"+BoundaryChar(right_boundary)+"_ boundary type:"+ &
                   &fu_str(pBBuf%iBoundaryType(right_boundary)), "get_linestuff_from_boundaries")
                return
             endselect
          endif !right boundary
  end  subroutine get_linestuff_from_boundaries

  subroutine init_molec_diff(species, a_half, b_half, xplus, xminus, nsp, nz, seconds)
    ! !Roux notebook p.50
    !Fills matrix for molecular diffusion dumb and simple, using standard atmosphere
    ! No fancy things...
    type(silam_species), dimension(:), intent(in) :: species 
    real, dimension (:), intent(in) :: a_half, b_half
    real(r8k), dimension (:, :), intent(out) :: xplus !Exchange factors upwards
    real(r8k), dimension (:, :),   intent(out) :: xminus !Exchange factors downwards
    integer, intent(in) ::  nsp, nz
    real, intent(in) :: seconds !Timestep

    integer :: iLev, iSp
    real :: fStdDiff, mu_over_mu_air_1 ! Diffusivity in air, mu_air/mu_sp - 1
    real :: fDiffusivity, rho_bott, p_bott, TT0bott !Values at layer bottom interface
    real :: p_top !Values at layer top interface
    real :: p_this, dp_this, p_prev ! Values for a layer
    real :: g_over_dp_subst_prev, g_over_dp_subst_this !Adjusted for current species
    real :: rho_over_R_down  !!Adjusted for current species

    ! Hard boundaries
    xminus(1,:) = 0.
    xplus(nz,:) = 0.
    do iSp = 1,  nsp
      
      fStdDiff = fu_gas_molecular_diffusivity_air(species(iSp)%material)
!      call report(species(iSp)%material)
      if (fStdDiff <=0 ) then  ! No diffusion for the species
         xplus(:,iSp) = 0.
         xminus(:,iSp) = 0.
         call msg("No diffusion data for"//trim(fu_str(species(iSp))))
         cycle
      endif
      mu_over_mu_air_1 =  fu_mole_mass(species(iSp)%material)/molecular_weight_air - 1.
      call msg("Filling molecular diffusion matrix for "//trim(fu_str(species(iSp))))

      ! First layer
      p_bott = a_half(1) + b_half(1)*std_pressure_sl
      p_top =  a_half(2) + b_half(2)*std_pressure_sl
      p_this = 0.5 * (p_bott + p_top)
      dp_this = p_bott - p_top
      g_over_dp_subst_this = g * exp(-log(p_this)* mu_over_mu_air_1) / dp_this

      do iLev = 2, nz
        p_prev  = p_this
        p_bott  = p_top
        g_over_dp_subst_prev = g_over_dp_subst_this
        
        p_top = a_half(iLev+1) + b_half(iLev+1)*std_pressure_sl
        p_this = 0.5 * (p_bott + p_top)
        dp_this = p_bott - p_top
        g_over_dp_subst_this = g * exp(-log(p_this)* mu_over_mu_air_1) / dp_this

        !Scale fiffusivity to the bottom interface
        call std_atm_rho_t(p_bott, rho_bott, TT0bott) 
        fDiffusivity = fStdDiff * (std_pressure_sl/p_bott) * (TT0bott)*sqrt(TT0bott)
        rho_over_R_down = fDiffusivity * g * 1.2 *1.2 *rho_bott * rho_bott / (p_prev - p_this) * exp(log(p_bott)* mu_over_mu_air_1) 

        xminus(iLev, iSp)  = g_over_dp_subst_this * rho_over_R_down * seconds !Fraction from this layer down
        xplus(iLev-1, iSp) = g_over_dp_subst_prev * rho_over_R_down * seconds !Fraction from prev layer up
      enddo
     do iLev = 1, nz
        call msg("iLev frac_upward frac_downward (per year)", (/real(iLev), real(xplus(iLev, iSp) / seconds *  31557600.) , real(xminus(iLev, iSp)  / seconds *  31557600.)/))
      enddo
    enddo

  end subroutine init_molec_diff  
!
subroutine make_molec_diffusion_passengers(xplus, xminus,  passengers, nz, Q)
      !Just diffuse with explicit scheme: diffusion assumed very slow
      ! and hardly has any effect except for the upper stratosphere...
      real(r8k), dimension(:), intent(in) :: xplus, xminus
      real, dimension(0:, 0:), intent(inout) :: passengers
      integer, intent(in) ::  nz
      ! No need to output  Just use it as scratch 
      real(r8k), dimension(0:,0:), intent(out) :: Q 

      integer  :: iLev 

    ! Just delta to the masses and moments
    Q(:,1) = passengers(:,2) * xminus(2) - passengers(:,1) * xplus(1)
    do iLev = 2, nz-1
           Q(:,iLev) =  passengers(:,iLev-1) * xplus(iLev-1) &
                      & -  passengers(:,iLev) * (xplus(iLev) + xminus(iLev) ) &
                      & +  passengers(:,iLev+1) * xminus(iLev+1)
     enddo
    Q(:,nz) =   passengers(:,nz-1) * xplus(nz-1)  &
               &  - passengers(:,nz) * xminus(nz)

    passengers(:,1:nz) = passengers(:,1:nz) + Q(:,1:nz)
!    call msg("RelErr", sum(Q(0,1:nz))/sum(abs(Q(0,1:nz))+1D-50))

  end subroutine make_molec_diffusion_passengers


  !*********************************************************************************************

  subroutine report_inout_mass_stuff_v5(incout, MassMap, pBBuf, ifByLevel)
    !
    ! Prints the detailed info on incoming/outgoing masses of all species for all sources consisting the cloud
    !
    implicit none

    ! Imported parameters
    integer, intent(in) :: incout ! incoming/outgoing
    type(Tmass_map), intent(inout) :: MassMap !Whatever mass map that knows species
    type(TboundaryBuffer), pointer :: pBBuf
    logical, intent(in) :: ifByLevel

        ! Local variables
    integer :: i, iz, iSrc, iUnit, iSpecies, iLev
    integer :: nSp, nZ, nSrc, nbr_file_units
    integer, dimension(2) :: fUnits
    logical :: advType
    real, dimension(6) :: pTmp
    real :: massSP, massNP, fTmp1, fTmp2
    character(len=6) :: chTmp1, chTmp2, chTmp3, chTMp4
    type(T_Galp_v5_thread_stuff), pointer :: stuff

  
    ! Everything should be already collected to 0-th thread
    stuff => EulerStuff%Galp5(0) 

    nSrc = size(stuff%x0_mass, 1)
    nSp  = size(stuff%x0_mass, 2)
    nz   = size(stuff%x0_mass, 4)

    call msg('')
    if (incout == incoming) then 
      call msg('========== TIMESTEP GRID-INCOMING MASS REPORT ============')
    else
       call msg('========== TIMESTEP GRID-OUTGOING MASS REPORT ============')
       if (incout /= outgoing) then
          call set_error("Wrong incoming/outgoing index:"+fu_str(incout),"report_inout_mass_stuff" )
          return
       endif
    endif

    if(pbBuf%iBoundaryType(northern_boundary)==polar_boundary_type)then
      chTmp1 = 'Pole S'
    else
      chTmp1 = 'y<1   '
    endif
    if(pbBuf%iBoundaryType(southern_boundary)==polar_boundary_type)then
      chTmp2 = 'Pole N'
    else
      chTmp2 = 'y>yMax'
    endif
    
    fUnits(1:2) = (/run_log_funit, 6/)
    if (smpi_global_rank == 0) then
      nbr_file_units = 2
    else
      nbr_file_units = 1
    end if

    do iUnit = 1,nbr_file_units
      do iSrc=1,nSrc
        if(nSrc > 1) write(fUnits(iUnit), *) 'Source: ',iSrc
        if (ifByLevel) then
           if(fu_ifLonGlobal(dispersion_grid))then
             write(fUnits(iUnit), '(4A)')' Species       iLev  ', chTmp1, '     ' ,chTmp2
           else
             write(fUnits(iUnit), '(4A)')' Species       iLev  x<1       x>xMax      ', chTmp1, '     ' ,chTmp2
           endif
           do iSpecies = 1, nSp
             if(fu_ifLonGlobal(dispersion_grid))then
               do iLev = 1,nz
                 write(fUnits(iUnit), '(A14,1x,I4,2(1x,E10.5))')&
                    & fu_str(MassMap%species(iSpecies)), iLev, &
                    & stuff%y0_mass(iSrc, iSpecies, incout, iLev), &
                    & stuff%yM_mass(iSrc, iSpecies, incout, iLev)
               enddo
             else
               do iLev = 1,nz
                 write(fUnits(iUnit), '(A14,1x,I4,4(1x,E10.5))') &
                    & fu_str(MassMap%species(iSpecies)), iLev, &
                    & stuff%x0_mass(iSrc, iSpecies, incout, iLev), &
                    & stuff%xM_mass(iSrc, iSpecies, incout, iLev), &
                    & stuff%y0_mass(iSrc, iSpecies, incout, iLev), &
                    & stuff%yM_mass(iSrc, iSpecies, incout, iLev)
               enddo 
             endif !Global
             write(fUnits(iUnit), '(1x)')  
           end do
       call msg('---------------Total integrated-----------------')
       endif

      if(fu_ifLonGlobal(dispersion_grid))then
        write(fUnits(iUnit), '(5A)')' Species        ', chTmp1, '     ' ,chTmp2, '      z<1       z>zMax'
      else
        write(fUnits(iUnit), '(5A)')' Species       x<1       x>xMax      ', chTmp1, '     ' ,chTmp2, '      z<1       z>zMax'
      endif

      do iSpecies = 1, nSp
        if(fu_ifLonGlobal(dispersion_grid))then
          write(fUnits(iUnit), '(A14,3(E10.5,1x),E10.5)')&
               & fu_str(MassMap%species(iSpecies)), &
               & sum(stuff%y0_mass(iSrc, iSpecies, incout, 1:nz_dispersion)), &
               & sum(stuff%yM_mass(iSrc, iSpecies, incout, 1:nz_dispersion)), &
               & stuff%bottom_mass(iSrc, iSpecies, incout), &
               & stuff%top_mass(iSrc, iSpecies, incout)
        else
          write(fUnits(iUnit), '(A14,5(E10.5,1x),E10.5)')&
               & fu_str(MassMap%species(iSpecies)), &
               & sum(stuff%x0_mass(iSrc, iSpecies, incout, 1:nz_dispersion)), &
               & sum(stuff%xM_mass(iSrc, iSpecies, incout, 1:nz_dispersion)), &
               & sum(stuff%y0_mass(iSrc, iSpecies, incout, 1:nz_dispersion)), &
               & sum(stuff%yM_mass(iSrc, iSpecies, incout, 1:nz_dispersion)), &
               & stuff%bottom_mass(iSrc, iSpecies, incout), &
               & stuff%top_mass(iSrc, iSpecies, incout)
        endif
      end do

      if(nSp > 3 .and. nSrc > 1)then
        pTmp(1:6) = 0.0
        if(error)return
        pTmp(1) = pTmp(1) + sum(    stuff%x0_mass(iSrc, :, incout, :))
        pTmp(2) = pTmp(2) + sum(    stuff%xM_mass(iSrc, :, incout, :))
        pTmp(3) = pTmp(3) + sum(    stuff%y0_mass(iSrc, :, incout, :))
        pTmp(4) = pTmp(4) + sum(    stuff%yM_mass(iSrc, :, incout, :))
        pTmp(5) = pTmp(5) + sum(stuff%bottom_mass(iSrc, :, incout))
        pTmp(6) = pTmp(6) + sum(   stuff%top_mass(iSrc, :, incout))
        write(fUnits(iUnit),*)'Source total:'
        write(fUnits(iUnit),'(15x,6(E10.5,1x))') pTmp(1:6)
      endif
      end do  !iSrc
    end do
    if (incout == incoming) then 
      call msg('---------- TIMESTEP END OF INCOMING MASS REPORT ----------')
    else
      call msg('---------- TIMESTEP END OF OUT-OF-GRID MASS REPORT ----------')
    endif
   call msg('')


  end subroutine report_inout_mass_stuff_v5

  
  !***********************************************************************************

  subroutine test_advect_mass_many()
    implicit none
    integer :: i
    !
    ! Go through all relaxation strengths
    !
    CMrelax = 1.0
    do i = 1, 11
!      call test_advect_massR() 
      call test_advect_massMAS() 
!      call test_advect_mass2() 
      CMrelax = CMrelax - 0.01
    end do
    CMrelax = 0.5
!    call test_advect_massR() 
    call test_advect_massMAS() 
!    call test_advect_mass2() 
    CMrelax = 0.0
!    call test_advect_massR() 
    call test_advect_massMAS() 
!    call test_advect_mass2() 
    
  end subroutine test_advect_mass_many
  

  !***********************************************************************************


  subroutine test_advect_massR()
    !
    ! Version with output for Roux Gnuplot scripts..
    !
    implicit none
    
     integer, parameter :: nLine = 100 !300 !100
     integer, parameter :: nTests = 9
     integer, parameter :: nSteps = 900
     integer :: W
     real, parameter :: fMax = 100.0
     real*8 :: Cu 
     integer :: delta, output_every = 60 !90 !60
     real, dimension(0:0,0:nLine+1,2) :: passengersInOut = 0
     real, dimension(nLine,0:nSteps) :: Stepresult
     real, dimension(nLine,2) :: cmInOut =0
     real, dimension(nLine) :: inM
     real, dimension(nLine,nTests) :: outM
     real*8, dimension(0:nLine) :: xr, xr_advected =0
     real, dimension(1:nLine) :: dx
     integer :: ix,istep,  iflip, itest, iUnit_rect, iUnit_rectdiff, iUnit_triangle, iUnit_rectdiff_cm, iUnit_step
     character(len=20) :: testname,chMM


     dx(:) = 1.
     xr(0) = 0.
     call msg("#Nline = ", nLine)
!     call msg("#Cu = ", Cu)
     call msg("#W = ", W)

     do ix=1,nLine
       xr(ix) = xr(ix-1) + dx(ix)
     enddo
     
     do  itest=1,nTests
        inm = 0.0
        W = 15
        select case (itest)
           case (1)  
              testname="rectpulse"
              do ix=1,W
!                 inm(ix+nline/2-w/2) = fMax
                 inm(1+ix) = fMax
              enddo
              inm = inm+10
           case (2)
              testname="tripulse"
              W = 20
              do ix=1, W
!                 inm(ix+nline/2-w/2) =  (W - 2*abs(ix- W/2.)) * fMax / real(W)
                 inm(1 + ix) =  (W - 2*abs(ix- W/2.)) * fMax / real(W)
              enddo
              inm = inm+10
           case (3)
              testname="sinpulse"
              do ix=1,W+1
!                 inm(ix+nline/2-w/2) =  1. + 10*(sin(pi*ix/W)**2)
                 inm(1+ix) = fMax*(sin(pi*(ix-1)/W)**2)
              enddo
              inm = inm+10
           case (4)  
              testname="rectgap"
              inm = fMax   !Background cnc
              do ix=1,W
!                 inm(ix+nline/2-w/2) = 0.1 * fMax
                 inm(1+ix) = 0.1 * fMax
              enddo
           case (5)
              testname="trigap"
              inm = fMax   !Background cnc
              do ix=1,W
!                 inm(ix+nline/2-w/2) =  fMax* (0.1 + 1.8  * abs(ix- W/2.) / W)
                 inm(1+ix) =  fMax* (0.1 + 1.8  * abs(ix- W/2.) / W)
              enddo
           case (6)
              testname="singap"
              inm = fMax   !Background cnc
              do ix=1,W
!                 inm(ix+nline/2-w/2) =  (1. - 0.9*(sin(pi*ix/W)**2)) * fMax
                 inm(1+ix) =  (1. - 0.9*(sin(pi*ix/W)**2)) * fMax
              enddo
           case(7)
              testname="delta"
              inm(w/2) = fMax
           case(8)
              inm = 0.1 * fMax
              testname="deltabckgr"
              inm(w/2) = fMax
           !   inm(nLine - 2) = fMax
        endselect
        Stepresult(1:nLine,0)=inm(1:nLine)
   
        cu = nline*1.0D0 / nsteps
        xr_advected(0:nLine) = xr(0:nLine) + Cu

        iUnit_rect = fu_next_free_unit()
        OPEN(iUnit_rect,file="00"+testname+"_rectangle.txt",action='write')
        iUnit_triangle = fu_next_free_unit()
        OPEN(iUnit_triangle,file="00"+testname+"_triangle.txt",action='write')
        iUnit_step = fu_next_free_unit()
        OPEN(iUnit_step,file="00"+testname+"_step.txt",action='write')

        write (iunit_rect,*) "#Test name:", testname
        write (iunit_rect,*) "#Steps per turn:", nsteps
        write (iunit_rect,*) "#Courrant:", cu
!        write (iunit_rect,'(i5,1x,A10,1000(F9.3,1x))') 0, "Initial = ", (inm(ix),ix=1,nLine), sum(inm(:))


        write (iUnit_triangle,*) "#Test name:", testname
        write (iUnit_triangle,*) "#Steps per turn:", nsteps
        write (iUnit_triangle,*) "#Courrant:", cu
!        write (iunit_triangle,'(i5,1x,A10,1000(F9.3,1x))') 0, "Initial = ", (inm(ix),ix=1,nLine), sum(inm(:))

        write (iUnit_step,*) "#Test name:", testname
        write (iUnit_step,*) "#Steps per turn:", nsteps
        write (iUnit_step,*) "#Courrant:", cu
!        write (iunit_step,'(i5,1x,A10,1000(F9.3,1x))') 0, "Initial = ", (inm(ix),ix=1,nLine), sum(inm(:))

        ! Silam default
        !
        passengersInOut(0,1:nLine,1) = inm(:)
        cmInOut(:,:) = 0
        iflip = 1
        do istep=1, nsteps
            call advect_mass(cmInOut(:,iflip), cmInOut(:,3-iflip), &
                           & passengersInOut(:,:,iflip),  passengersInOut(:,:,3-iflip), 0, &
                           & xr, xr_advected, dx,&
                           & .true., 1, nLine, nLine, &! ifLoop, iCellStart, iCellEnd, nCells,
               & 1e-10, .false., .True., 1.0) !& fMinAdvectedMass, have_negatives, ifCMout)
            iflip = 3 -  iflip
            Stepresult(1:nLine,iStep)= passengersInOut(0,1:nLine,iflip)
        enddo
        write (iunit_rect,'(A2,1000(F9.3,1x))')"#",sum(Stepresult(:,0:nSteps:9),dim=1)
        do ix=1,nLine
             write (iunit_rect,'(1000(F9.3,1x))') Stepresult(ix,0:nSteps:9)
        enddo 
            
        ! Pure triangle
        !
        passengersInOut(0,1:nLine,1) = inm(:)
        cmInOut(:,:) = 0
        iflip = 1
        do istep=1, nsteps
           where (passengersInOut(:,:,iflip)< 1e-20) passengersInOut(:,:,iflip) = 0.
           call advect_mass_trislab(cmInOut(:,iflip), cmInOut(:,3-iflip), &
                                  & passengersInOut(:,:,iflip),  passengersInOut(:,:,3-iflip), 0, &
                                  & xr, xr_advected, dx,&
                                  & .true., 1, nLine, nLine, &! ifLoop, iCellStart, iCellEnd, nCells,
                                  & 1e-10, .false., .True., 1.0) !& fMinAdvectedMass, have_negatives, ifCMout)
           iflip = 3 -  iflip
            Stepresult(1:nLine,iStep)= passengersInOut(0,1:nLine,iflip)
        enddo
        write (iunit_triangle,'(A2,1000(F9.3,1x))')"#",sum(Stepresult(:,0:nSteps:9),dim=1)
        do ix=1,nLine
             write (iunit_triangle,'(1000(F9.3,1x))') Stepresult(ix,0:nSteps:9)
        enddo 


        ! Step
        !
        passengersInOut(0,1:nLine,1) = inm(:)
        cmInOut(:,:) = 0
        iflip = 1
        do istep=1, nsteps
!            call msg("")
!           print *,"Step", istep
           where (passengersInOut(:,:,iflip)< 1e-20) passengersInOut(:,:,iflip) = 0.
           call advect_mass_step(cmInOut(:,iflip), cmInOut(:,3-iflip), &
                                  & passengersInOut(:,:,iflip),  passengersInOut(:,:,3-iflip), 0, &
                                  & xr, xr_advected, dx,&
                                  & .true., 1, nLine, nLine, &! ifLoop, iCellStart, iCellEnd, nCells,
                                  & 1e-10, .false., .True., 1.0) !& fMinAdvectedMass, have_negatives, ifCMout)
           iflip = 3 -  iflip
           Stepresult(1:nLine,iStep)= passengersInOut(0,1:nLine,iflip)
        enddo
        write (iunit_step,'(A2,1000(F9.3,1x))')"#",sum(Stepresult(:,0:nSteps:9),dim=1)
        do ix=1,nLine
             write (iunit_step,'(1000(F9.3,1x))') Stepresult(ix,0:nSteps:9)
        enddo 

        close(iUnit_rect)
        close(iUnit_triangle)
        close(iUnit_step)

     enddo ! test patterns

  endsubroutine test_advect_massR


  !***********************************************************************************

  subroutine test_advect_massMAS()
     ! MAS version
    implicit none

    integer, parameter :: nLine = 100
     integer, parameter :: nTests = 8
     integer :: W
     real, parameter :: fMax = 100.0
     real*8 :: Cu 
     integer :: delta, output_every = 60 !90 !60
     real, dimension(0:0,0:nLine+1,2) :: passengersInOut = 0
     real, dimension(nLine,2) :: cmInOut =0
     real, dimension(nLine) :: inM
     real, dimension(nLine,nTests) :: outM
     real*8, dimension(0:nLine) :: xr, xr_advected =0
     real, dimension(1:nLine) :: dx
     integer :: ix,istep, nsteps, iflip, itest, iUnit_rect, iUnit_rectdiff, iUnit_triangle, iUnit_rectdiff_cm
     character(len=20) :: testname,chMM
     logical :: myErr


     dx(:) = 1.
     xr(0) = 0.
     nsteps=9999
     call msg("#Nline = ", nLine)
!     call msg("#Cu = ", Cu)
     call msg("#W = ", W)

     do ix=1,nLine
       xr(ix) = xr(ix-1) + dx(ix)
     enddo
     
     do  itest=0,nTests
        inm = 0.0
        W = 15 ! 25  !15
        select case (itest)
          case(0)
            testname="plain"
            inm(:) = fMax
          case (1)  
              testname="rectpulse"
              do ix=1,W
!                 inm(ix+nline/2-w/2) = fMax
                 inm(1+ix) = fMax
              enddo
          case (2)
              testname="tripulse"
              W = 20
              do ix=1, W
!                 inm(ix+nline/2-w/2) =  (W - 2*abs(ix- W/2.)) * fMax / real(W)
                 inm(1 + ix) =  (W - 2*abs(ix- W/2.)) * fMax / real(W)
              enddo
          case (3)
              testname="sinpulse"
              do ix=1,W+1
!                 inm(ix+nline/2-w/2) =  1. + 10*(sin(pi*ix/W)**2)
                 inm(1+ix) = fMax*(sin(pi*(ix-1)/W)**2)
              enddo
          case (4)  
              testname="rectgap"
              inm = fMax   !Background cnc
              do ix=1,W
!                 inm(ix+nline/2-w/2) = 0.1 * fMax
                 inm(1+ix) = 0.1 * fMax
              enddo
          case (5)
              testname="trigap"
              inm = fMax   !Background cnc
              do ix=1,W
!                 inm(ix+nline/2-w/2) =  fMax* (0.1 + 1.8  * abs(ix- W/2.) / W)
                 inm(1+ix) =  fMax* (0.1 + 1.8  * abs(ix- W/2.) / W)
              enddo
          case (6)
              testname="singap"
              inm = fMax   !Background cnc
              do ix=1,W
!                 inm(ix+nline/2-w/2) =  (1. - 0.9*(sin(pi*ix/W)**2)) * fMax
                 inm(1+ix) =  (1. - 0.9*(sin(pi*ix/W)**2)) * fMax
              enddo
          case(7)
              testname="delta"
              inm(w/2) = fMax
          case(8)
              inm = 0.1 * fMax
              testname="deltabckgr"
              inm(w/2) = fMax
          case(9)
              inm = 0.
              testname="chair"
              do ix=1,W
                 inm(1+ix) =  0.5 * fMax
              enddo
              do ix=W,W*2
                 inm(1+ix) =  fMax
              enddo
        endselect
   
        cu = 0.4   !nline*1.0_r8k / nsteps
        xr_advected(0:nLine) = xr(0:nLine) + Cu

        iUnit_rect = fu_next_free_unit()
        OPEN(iUnit_rect,file="00"+testname+"_rectangle.txt",action='write')
        iUnit_rectdiff = fu_next_free_unit()
        write(unit=chMM,fmt='(F5.3)')maxmix
        print *, '==>',trim(chMM),'<=='
        OPEN(iUnit_rectdiff,file="00"+testname+"_rectangle_diff_"+trim(chMM)+".txt",action='write')
        iUnit_rectdiff_cm = fu_next_free_unit()
        write(unit=chMM,fmt='(F5.3)')CMRelax
        print *, '==>',trim(chMM),'<=='
        OPEN(iUnit_rectdiff_cm,file="00"+testname+"_rectangle_diff_cm_"+trim(chMM)+".txt",action='write')
        iUnit_triangle = fu_next_free_unit()
        OPEN(iUnit_triangle,file="00"+testname+"_triangle.txt",action='write')

        write (iunit_rect,*) "#Test name:", testname
        write (iunit_rect,*) "#Steps per turn:", nsteps
        write (iunit_rect,*) "#Courrant:", cu
        write (iunit_rect,'(i5,1x,A10,1000(F9.3,1x))') 0, "Initial = ", (inm(ix),ix=1,nLine), sum(inm(:))

        write (iunit_rectdiff,*) "#Test name:", testname
        write (iunit_rectdiff,*) "#Steps per turn:", nsteps
        write (iunit_rectdiff,*) "#Courrant, maxmix:", cu, maxmix
        write (iunit_rectdiff,'(i5,1x,A10,1000(F9.3,1x))') 0, "Initial = ", (inm(ix),ix=1,nLine), sum(inm(:))

        write (iunit_rectdiff_cm,*) "#Test name:", testname
        write (iunit_rectdiff_cm,*) "#Steps per turn:", nsteps
        write (iunit_rectdiff_cm,*) "#Courrant, maxmix:", cu, maxmix
        write (iunit_rectdiff_cm,'(i5,1x,A10,1000(F9.3,1x))') 0, "Initial = ", (inm(ix),ix=1,nLine), sum(inm(:))

        write (iUnit_triangle,*) "#Test name:", testname
        write (iUnit_triangle,*) "#Steps per turn:", nsteps
        write (iUnit_triangle,*) "#Courrant:", cu
        write (iunit_triangle,'(i5,1x,A10,1000(F9.3,1x))') 0, "Initial = ", (inm(ix),ix=1,nLine), sum(inm(:))

        ! Silam default
        !
        passengersInOut(0,1:nLine,1) = inm(:)
        cmInOut(:,:) = 0
        iflip = 1
        do istep=1, nsteps
          call advect_mass(cmInOut(:,iflip), cmInOut(:,3-iflip), &
                         & passengersInOut(:,:,iflip),  passengersInOut(:,:,3-iflip), 0, &
                         & xr, xr_advected, dx,&
                         & .true., 1, nLine, nLine, &! ifLoop, iCellStart, iCellEnd, nCells,
                         & 1e-10, .false., .True., 1.0) !& fMinAdvectedMass, have_negatives, ifCMout)
          iflip = 3 -  iflip
          if(mod(istep,output_every) == 0) write (iunit_rect,'(I5,1x,A10,1000(F9.3,1x))') istep, "Values = ", (passengersInOut(0,ix,iflip),ix=1,nLine), sum(passengersInOut(0,1:nline,iflip))
        enddo
        outM(:,1) =  passengersInOut(0,1:nLine,iflip)

        !
        ! Silam default with mass-centre diffusion
        !
        passengersInOut(0,1:nLine,1) = inm(:)
        cmInOut(:,:) = 0
        iflip = 1
        do istep=1, nsteps
          call advect_mass(cmInOut(:,iflip), cmInOut(:,3-iflip), &
                         & passengersInOut(:,:,iflip),  passengersInOut(:,:,3-iflip), 0, &
                         & xr, xr_advected, dx,&
                         & .true., 1, nLine, nLine, &! ifLoop, iCellStart, iCellEnd, nCells,
                         & 1e-10, .false., .True., & !& fMinAdvectedMass, have_negatives, ifCMout)
                         & CMRelax)
          iflip = 3 -  iflip
          if(mod(istep,output_every) == 0) write (iunit_rectdiff_cm,'(I5,1x,A10,1000(F9.3,1x))') istep, "Values = ", (passengersInOut(0,ix,iflip),ix=1,nLine), sum(passengersInOut(0,1:nline,iflip))
        enddo
        outM(:,1) =  passengersInOut(0,1:nLine,iflip)

        ! Pure triangle
        !
        passengersInOut(0,1:nLine,1) = inm(:)
        cmInOut(:,:) = 0
        iflip = 1
        do istep=1, nsteps
 !          print *, ""
 !          print *, "#Step=", istep-1
 !          do ix=1,nLine
 !           print *, ix,  passengersInOut(0,ix,iflip)
 !          enddo
          call advect_mass_trislab(cmInOut(:,iflip), cmInOut(:,3-iflip), &
                                 & passengersInOut(:,:,iflip),  passengersInOut(:,:,3-iflip), 0, &
                                 & xr, xr_advected, dx,&
                                 & .true., 1, nLine, nLine, &! ifLoop, iCellStart, iCellEnd, nCells,
                                 & 1e-10, .false., .True., 1.0) !& fMinAdvectedMass, have_negatives, ifCMout)
          iflip = 3 -  iflip
          if(mod(istep,output_every) == 0) write (iunit_triangle,'(I5,1x,A10,1000(F9.3,1x))') istep, "Values = ", (passengersInOut(0,ix,iflip),ix=1,nLine), sum(passengersInOut(0,1:nline,iflip))
        enddo
!        outM(:,2) =  passengersInOut(0,1:nLine,iflip)

        close(iUnit_rect)
        close(iUnit_rectdiff)
        close(iUnit_rectdiff_cm)
        close(iUnit_triangle)

!        write (iunit,*) "# M0    Standard triangle" ! triangle02 triangle04"
!        do ix=1,nLine
!            write (iunit,*) ix,  inm(ix), outM(ix,1), outM(ix,2) !, outM(ix,3), outM(ix,4) 
!        enddo
!        write (iunit,*), "#" sum(inm), sum(outm(:,1)), sum(outm(:,2)) !, sum(outm(:,3)), sum(outm(:,4))
!        close(iunit)
     enddo ! test patterns

  endsubroutine test_advect_massMAS

  !***********************************************************************************

  subroutine test_advect_mass2()
     implicit none
     integer, parameter :: nLine = 8
     real, dimension(0:0,0:nLine+1,2) :: passengersInOut = 0
     real, dimension(nLine,2) :: cmInOut =0
     real(r8k), dimension(0:nLine) :: xr, xr_advected =0
     real, dimension(1:nLine) :: dx
     integer :: ix,istep, nsteps, iflip, itest, iUnit
     character(len=20) :: testname

  xr = (/0.00000E+00, 0.15247E+12, 0.45545E+12, 0.10485E+13, 0.21797E+13, 0.41732E+13, 0.70728E+13, 0.11237E+14, 0.14675E+14/)
  xr_advected = (/-.00000E+00, 0.15424E+12, 0.46153E+12, 0.10611E+13, 0.22022E+13, 0.41973E+13, 0.70852E+13, 0.11276E+14, 0.14695E+14/)
     dx(1:nLine) = xr(1:nLine)-xr(0:nLine-1)

  passengersInOut(0,0:nLine+1,1) = (/0.00000E+00, 0.32388E+03, 0.75720E+03, 0.17067E+04, 0.29992E+04, 0.39749E+04, 0.37227E+04,  0.00000E+00, 0.00000E+00, 0.11391E+04/)
     cmInOut(:,1) = 0
     iFlip=1
     call advect_mass(cmInOut(:,iflip), cmInOut(:,3-iflip), &
               & passengersInOut(:,:,iflip),  passengersInOut(:,:,3-iflip), 0, &
                & xr, xr_advected, dx,&
               & .false., 1, nLine+1, nLine, &! ifLoop, iCellStart, iCellEnd, nCells,
               & 1e-10, .false., .True., 1.0) !& fMinAdvectedMass, have_negatives, ifCMout)
      print *, passengersInOut(0,1:nLine,2)

  endsubroutine test_advect_mass2
  
  subroutine test_advect_cell()
     implicit none
    integer, parameter :: nLine = 8
     real, dimension(0:0,0:nLine+1,2) :: passengersInOut = 0
     real, dimension(0:0,0:nLine+1,2) :: passengersInOut1 = 0
     real, dimension(nLine,2) :: cmInOut =0
     real, dimension(0:nLine) :: u_right
     real, dimension(0:nLine) :: u_right1
     logical, dimension(0:nLine+1) :: ifSkipCell
     real*8, dimension(0:nLine) :: xr, xr_advected, xr_advected1
     real, dimension(1:nLine) :: dx
     integer :: ix,istep, nsteps, iflip, itest, iUnit
     character(len=20) :: testname
     real, dimension(0:nLine) :: passtime
     real :: seconds
     integer :: i

    seconds=0.7
    do i = 0,nLine
       u_right(0:nLine) = (/2., 1., 0., 1., 2., 3., 2., 1., 0./)
    enddo
    dx(1:nLine) = 1.
    ifSkipCell = .false.

    call advect_cellboundaries(u_right, dx, xr, xr_advected,& 
       & ifSkipCell, passtime, seconds, .false., 1, nLine, nLine)

    print *, xr
    print *, xr_advected  ! Prints advected boundaries

    ! flip wind
    u_right1(0:nLine) = -u_right(nLine:0:-1)
    call advect_cellboundaries(u_right1, dx, xr, xr_advected1,& 
       & ifSkipCell, passtime, seconds, .false., 1, nLine, nLine)
    print *, xr(nLine)-xr_advected1(nLine:0:-1) !must print the same

 print *, ""
  passengersInOut(0,0:nLine+1,1) = (/1., 0., 0., 1., 0., 0., 0.,  0., 0., 0./)
  passengersInOut1(0,0:nLine+1,1) = passengersInOut(0,nLine+1:0:-1,1)
  cmInOut(:,1) = 0.2
     iFlip=1
     call advect_mass_trislab(cmInOut(:,iflip), cmInOut(:,3-iflip), &
               & passengersInOut(:,:,iflip),  passengersInOut(:,:,3-iflip), 0, &
                & xr, xr_advected, dx,&
               & .false., 0, nLine+1, nLine, &! ifLoop, iCellStart, iCellEnd, nCells,
               & 1e-10, .false., .True., 1.0) !& fMinAdvectedMass, have_negatives, ifCMout)
            print *, passengersInOut(0,0:nLine+1,2)
            print *, sum(passengersInOut(0,0:nLine+1,2))
     
  cmInOut(:,1) = -0.2
     call advect_mass_trislab(cmInOut(:,iflip), cmInOut(:,3-iflip), &
               & passengersInOut1(:,:,iflip),  passengersInOut1(:,:,3-iflip), 0, &
                & xr, xr_advected1, dx,&
               & .false., 0, nLine+1, nLine, &! ifLoop, iCellStart, iCellEnd, nCells,
               & 1e-10, .false., .True., 1.0) !& fMinAdvectedMass, have_negatives, ifCMout)

     print *, passengersInOut1(0,nLine+1:0:-1,2)
     print *, sum(passengersInOut1(0,0:nLine+1,2))

  endsubroutine test_advect_cell
  

END MODULE advection_eulerian_v5

