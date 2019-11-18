MODULE advection_lagrangian

  ! Here are the Lagrangian advection and random-walk for a bunch of Lagrangian particles,
  ! which describes the transport with mean wind and turbulent diffusion and other mixing
  ! processes in the atmosphere.
  !
  ! This module is expected to be fast. Tricks used:
  ! 1. There are no Lagrangian partiles themselves - just a bunch of arrays of 
  !    longitude, latitude, pressure, and masses, which together comprise the
  !    cloud of particles
  ! 2. Data are coming in as ready-made form as possible
  ! 3. Lagrangian advection takes place in meteorological grid, no reprojections.
  ! 4. Special preparatory function, which gets the meteofields and prepares
  !    them. It is called once per time step from advection.
  ! 5. Limited number of checking actions
  !
  ! IMPORTANT.
  ! This unit is "physically meaningfull", so the SI units are kept.
  !
  ! Author: Mikhail Sofiev
  !
  ! All units: SI
  ! 
  ! All degrees: real numbers in degrees (degrees and parts per
  ! hundred/parts per thousand), NOT in degrees and minutes
  ! 
  ! Language: ANSI standard Fortran 90
  ! 
  ! Modules used:

  use chemistry_manager
!  use field_buffer
  use lagrange_particles

  IMPLICIT NONE

  private

  ! The public functions and subroutines available in this module:
  public add_lagr_advection_input_needs
  public InitLagrAdvectionFields
  public advect_lagrangian_cloud

  PUBLIC add_random_walk_input_needs
!  PUBLIC random_walk

  ! Private functions and subroutines
  private lagrangian_adv_3d_wind_midpoint
  private lagrangian_adv_2d_wind_midpoint
  private lagrangian_adv_3d_wind_endpoint
  private lagrangian_adv_3d_implicit
  private handle_out_of_grid
  private check_pressure_range
  private prepare_random_walk
  PRIVATE random_walk_gaussian_prof
  PRIVATE bulk_gaussian_rw_turb
  private reflection_rules
  !
  ! Types of lagrangian advection
  !
  integer, parameter, public :: adv_lagrange_3d_wind_midpoint = 100021
  integer, parameter, public :: adv_lagrange_2d_wind_midpoint = 100022
  integer, parameter, public :: adv_lagrange_3d_wind_endpoint = 100023
  integer, parameter, public :: adv_lagrange_3d_implicit = 100024
  !
  ! Types of random walk 
  !
  INTEGER, PARAMETER, PUBLIC :: void_rw = 100030
  INTEGER, PARAMETER, PUBLIC :: fixed_rw = 100031
  INTEGER, PARAMETER, PUBLIC :: fully_mixed_abl_rw = 100032
  INTEGER, PARAMETER, PUBLIC :: bulk_gaussian_rw = 100033
  INTEGER, PARAMETER, PUBLIC :: advanced_gaussian = 100034
  !
  ! Accuracy of iterative advection
  !
  INTEGER, PARAMETER, private :: maxiterations = 5
  INTEGER, PARAMETER, private :: epsilon_distance = 0.001 ! fraction grid cell size
  INTEGER, PARAMETER, private :: epsilon_pressure = 1. ! Pa

  ! Pointers to the fields in meteobuffer. Prepared by prepare_random_walk
  ! and used by above private random-walk functions
  !
  REAL,DIMENSION(:),POINTER,PRIVATE,SAVE :: abl_h_m_ptr, MO_len_inv_ptr, fric_vel_ptr, conv_vel_scale_ptr
  TYPE(field_4d_data_ptr),POINTER,PRIVATE,SAVE :: temp_4d_ptr, pr_4d_ptr, height_4d_ptr
!  type(field_2d_data_ptr), pointer, private, save :: abl_h_pa_ptr, surf_pr_ptr
  integer, private, save :: ind_HABL_Pa, ind_SrfPress

  !
  ! Diffusion is in pressure vertical, so here it is
  !
  real, parameter, private :: pressure_top = 10., pressure_bottom = 110000.
  type(silam_vertical), private, save :: pressure_vertical

  !
  ! Other internal variables
  !
  REAL, PARAMETER, private :: a1=0.25, a2=0.5, b=0.875  ! for filly-mixed ABL
  real, parameter, private :: diffVelocityFT = 0.5 ! Pa/sec, free-troposphere move


CONTAINS


  ! ****************************************************************

  SUBROUTINE add_lagr_advection_input_needs(adv_method, &
                                     & q_met_dynamic, q_met_st, q_disp_dynamic, q_disp_st)
  !
  ! Returns the list of needed quantities for a particular advection
  ! algorithm. 
  ! NOTE. It returns only upper_level parameters (e.g. w_height), but
  !       does not care how they will be computed from the "raw" data.
  !
  IMPLICIT NONE

    ! Input parameters:
    INTEGER, INTENT(in) :: adv_method

    ! Output parameters:
    INTEGER, DIMENSION(:), INTENT(inout) :: q_met_dynamic, q_met_st, q_disp_dynamic, q_disp_st

    ! Local variables
    integer :: iTmp

    ! Wind & pressure:
    iTmp = fu_merge_integer_to_array(u_flag,       q_met_dynamic)
    iTmp = fu_merge_integer_to_array(v_flag,       q_met_dynamic)
    iTmp = fu_merge_integer_to_array(pressure_flag,q_met_dynamic)
    iTmp = fu_merge_integer_to_array(height_flag,  q_met_dynamic)
    iTmp = fu_merge_integer_to_array(temperature_flag,      q_met_dynamic)
    iTmp = fu_merge_integer_to_array(surface_pressure_flag, q_met_dynamic)
!    iTmp = fu_merge_integer_to_array(cell_size_z_flag, q_met_dynamic)

!    iTmp = fu_merge_integer_to_array(cell_size_x_flag, q_met_st)
!    iTmp = fu_merge_integer_to_array(cell_size_y_flag, q_met_st)

    SELECT CASE(adv_method)
      CASE(adv_lagrange_3d_wind_midpoint, adv_lagrange_3d_wind_endpoint, adv_lagrange_3d_implicit)
        call msg_test('3D advection')
        iTmp = fu_merge_integer_to_array(omega_flag, q_met_dynamic)

      CASE(adv_lagrange_2d_wind_midpoint)
        call msg_test('2D advection')

      CASE DEFAULT
        CALL set_error('Unknown advection method','add_lagr_advection_input_needs')
        RETURN
    END SELECT

  END SUBROUTINE add_lagr_advection_input_needs


  !********************************************************************************

  subroutine InitLagrAdvectionFields(adv_method, nspecies, nAerosolSpecies)
    !
    ! Initialises the internal advection fields.
    ! It handles only internal fields, which are private for the module.
    !
    implicit none

    ! Imported parameter
    integer, intent(in) :: adv_method
    integer, intent(in) :: nspecies, nAerosolSpecies

    ! Local variables
    integer :: i,j,iLev,iStatus, iSubst

    !
    ! Anything method-specific?
    !
    select case(adv_method)

      case(adv_lagrange_3d_wind_midpoint, adv_lagrange_2d_wind_midpoint, &
         & adv_lagrange_3d_wind_endpoint, adv_lagrange_3d_implicit)
        !
        ! Lagrangian advection does not have own fields
        !
      case default
        call set_error('Unknown advection','InitLagrAdvectionFields')
    end select
    
  end subroutine InitLagrAdvectionFields



  !===============================================================================
  !===============================================================================
  !===============================================================================
  !
  !    Lagrangian advecion
  !
  !===============================================================================
  !===============================================================================
  !===============================================================================


  ! ***************************************************************

  subroutine advect_lagrangian_cloud(adv_method, rw_method, &
                                   & arDyn, arMassTrn, nop, nSpecies, &
                                   & lpStatus, iArBadParticle, &
                                   & pMetBuf, & !pDispBuf, &
                                   & seconds, weight_past, &
                                   & wdr)
    !
    ! Calls an appropriate advection, depending on the advection_method
    !
    implicit none

    ! Imported parameters
    INTEGER, INTENT(in) :: adv_method, rw_method, nop, nSpecies
    real, dimension(:,:),  INTENT(inout) :: arDyn   ! (9, nParticles)
    real, dimension(:,:), INTENT(inout) :: arMassTrn   ! (nSpecies, nParticles)
    type(Tfield_buffer), INTENT(in) :: pMetBuf !, pDispBuf
    real, intent(in) :: seconds, weight_past
    integer, dimension(:), intent(inout) :: lpStatus
    integer, dimension(:), intent(inout) :: iArBadParticle
    TYPE(silja_wdr), INTENT(in) :: wdr

    select case(adv_method)
      case(adv_lagrange_3d_wind_midpoint)
        Call msg('Doing lagrangian_adv_3d_wind_midpoint')
        call lagrangian_adv_3d_wind_midpoint(rw_method, &
                                           & arDyn, arMassTrn, nop, nSpecies, &
                                           & lpStatus, iArBadParticle, &
                                           & pMetBuf, & !pDispBuf, &
                                           & seconds, weight_past, &
                                           & wdr)
      case(adv_lagrange_2d_wind_midpoint)
        Call msg('Doing lagrangian_adv_2d_wind_midpoint')
        call lagrangian_adv_2d_wind_midpoint(rw_method, &
                                           & arDyn, arMassTrn, nop, nSpecies, &
                                           & lpStatus, iArBadParticle, &
                                           & pMetBuf, & !pDispBuf, &
                                           & seconds, weight_past, &
                                           & wdr)
      case(adv_lagrange_3d_wind_endpoint)
        Call msg('Doing lagrangian_adv_3d_wind_endpoint')
        call lagrangian_adv_3d_wind_endpoint(rw_method, &
                                           & arDyn, arMassTrn, nop, nSpecies, &
                                           & lpStatus, iArBadParticle, &
                                           & pMetBuf, & !pDispBuf, &
                                           & seconds, weight_past, &
                                           & wdr)
      case(adv_lagrange_3d_implicit)
        call lagrangian_adv_3d_implicit(rw_method, &
                                      & arDyn, arMassTrn, nop, nSpecies, &
                                      & lpStatus, iArBadParticle, &
                                      & pMetBuf, & !pDispBuf, &
                                      & seconds, weight_past, &
                                      & wdr)
      case default
        call set_error('Unknown advection method:' + fu_str(adv_method),'advect_lagrangian_cloud')
    end select

  end subroutine advect_lagrangian_cloud


  ! ***************************************************************

  subroutine lagrangian_adv_3d_wind_midpoint(rw_method, &
                                           & arDyn, arMassTrn, nop, nSpecies, &
                                           & lpStatus, iArBadParticle, &
                                           & pMetBuf, & !pDispBuf, &
                                           & seconds, weight_past, &
                                           & wdr)
    !
    ! Advects all Lagrangian particles in the cloud. 
    ! Note that the Lagrangian particles are advected horizontally in the meteo grid realtive 
    ! coordinates, whereas vertical is taken in pressure coordinates.
    ! We also assume that the meteo vertical can be described as hybrid vertical: we use
    ! hybrid coefficients for fast computations.
    !
    implicit none

    ! Imported parameters
    INTEGER, INTENT(in) :: rw_method, nop, nSpecies
    real, dimension(:,:), intent(inout) :: arDyn   ! (9, nParticles)
    real, dimension(:,:), intent(inout) :: arMassTrn   ! (nSpecies, nParticles)
    type(Tfield_buffer), intent(in) :: pMetBuf !, pDispBuf
    real, intent(in) :: seconds, weight_past
    integer, dimension(:), intent(inout) :: lpStatus
    integer, dimension(:), intent(inout) :: iArBadParticle
    TYPE(silja_wdr), INTENT(in) :: wdr

    ! Local declarations:
    INTEGER, DIMENSION(:), POINTER :: q_arr_ptr
    real, dimension(:), pointer :: pXCellSize, pYCellSize, a_full, b_full, a_half, b_half
    INTEGER :: iPart, horiz_interp_method, vert_interp_method, iOutsideAction, ix, iy, iz, &
             & iTmp, p_ind, u_ind, v_ind, w_ind, sp_ind, h_ind, t_ind, grid_index, iter, i1d, &
             & iwp_x, iwp_y, iwp_z, iBadParticle, interp_q_2d, interp_q_3d
    REAL :: move_x, move_y, move_z, move_x_old, move_y_old, move_z_old, fZC, fZC_rel, fTmp, ps, &
          & wind_point_x, wind_point_y, wind_point_z, wind_point_rel_z, wind_u, wind_v, wind_w,pCentre
    logical :: ifOutside
    logical, save :: ifFirst=.true.

    real, parameter :: itertol = 1e-2
    integer, parameter :: niter = 5

    !
    ! If we are here for the first time, store locally the vertical structure of meteo grid
    !
    if(ifFirst)then
      allocate(b_full(nz_meteo+1), a_full(nz_meteo+1), b_half(0:nz_meteo), a_half(0:nz_meteo), &
             & stat=iz)
      if(fu_fails(iz==0,'Failed hybrid coefs allocation','lagrangian_adv_3d_wind_midpoint'))return
      call hybrid_coefs(meteo_vertical, a_full=a_full, b_full=b_full)
      if(error)return
      a_half(0) = 0.
      b_half(0) = 1.
!      a_full(0) = 0.
!      b_full(0) = 2.
      a_full(nz_meteo+1) = a_full(nz_meteo) * (a_full(nz_meteo)/a_full(nz_meteo-1))**1.5
      b_full(nz_meteo+1) = max(0.0 , b_full(nz_meteo) + (b_full(nz_meteo) - b_full(nz_meteo-1)) * 1.5)
      do iz = 1, nz_meteo
        a_half(iz) = (a_full(iz) + a_full(iz+1)) / 2.
        b_half(iz) = (b_full(iz) + b_full(iz+1)) / 2.
      end do
      ifFirst = .false.
    endif  ! ifFirst

    !
    ! Get wind and pressure fields for past and future.
    ! Also prepare indices for needed fields
    !
    q_arr_ptr => pMetBuf%buffer_quantities
    p_ind = fu_index(q_arr_ptr, pressure_flag)
    h_ind = fu_index(q_arr_ptr, height_flag)
    t_ind = fu_index(q_arr_ptr, temperature_flag)
    sp_ind = fu_index(q_arr_ptr,surface_pressure_flag)
    u_ind = fu_index(q_arr_ptr, u_flag)
    v_ind = fu_index(q_arr_ptr, v_flag)
    w_ind = fu_index(q_arr_ptr , omega_flag)
    !
    ! Instead of using 4D cube with pressure in meteo_vertical, we shall use the 
    ! hybrid coefficients and surface pressure. Easier to calculate, surface pressure is always
    ! available as "present".
    ! However, meteo_vertical is thin-layer structure, i.e. only full levels are defined.
    ! We shall have to use their half-sum to get the half-levels.
    !
    ! Below, two versions are presented: for arDyn keeping absolute pressure and relative coord 
    ! in meteo_vertical. Comment/uncomment them and related pieces in other modules to switch.
    !

!    ifPrintDump = .false. !now > tDebugTime !.and. now < tDebugTime + one_day * 2.

    !
    ! Set pointers in the random-walk routines
    !
!call msg('Random walk preparation starts')

    CALL prepare_random_walk(rw_method, pMetBuf)
    if(error)then
      call msg_warning('Random walk preparation failed','lagrangian_adv_3d_wind_midpoint')
      return
    endif

    horiz_interp_method = fu_horizontal_interp_method(wdr)
    vert_interp_method = linear ! So far.
    !
    !  Prepare pointers for the meteo grid cell sizes. Needed in advection
    ! in order to scale the absolute wind into the relative particle motion (and particles
    ! move in advection_grid to minimise reprojections).
    !
    pXCellSize => fu_grid_data(meteo_cell_x_size_fld)
    pYCellSize => fu_grid_data(meteo_cell_y_size_fld)

    if(fu_ifLonGlobal(meteo_grid))then
      iOutsideAction = handleGlobalGrid
    else
      iOutsideAction = nearestPoint !notAllowed
    endif

    if(error)then
      call set_error('Cannot start advection loop','lagrangian_adv_3d_wind_midpoint')
      return
    endif

    iBadParticle = 1
    iArBadParticle(iBadParticle) = int_missing

    !----------------------------------------
    !
    ! Grand loop over particles.
    !
    particle: do iPart = 1, nop
      !
      ! Skip the particle if it is not active
      !
      if(lpStatus(iPart) == int_missing)cycle

!print *, 'Cloud_advection_3'
!call stop_count(chCounterNm = strTmp)

      if(error) return

if(ifPrintDump) then
  call msg('adv x=', arDyn(lp_x,iPart))
  call msg('adv y=', arDyn(lp_y,iPart))
endif
      !
      ! Advection largely follows the Eulerian Galperin's scheme v.3 (implicit)
      ! except for the Eulerian step, which is not needed, of course. Also, it all goes in 
      ! meteo grid
      ! A complexity is: vertical axis is very inhomogeneous, so have to fly in 
      ! absolute pressure. It is also nice since random-walk is also in pressure to keep 
      ! const-mixing-ratio
      !  
      ix = nint(arDyn(lp_x,iPart))
      iy = nint(arDyn(lp_y,iPart))
      i1d = ix + (iy - 1) * nx_meteo
      ps = pMetBuf%p2d(sp_ind)%present%ptr(i1d)

! arDyn pressure--------------------------------
!      ! Turn pressure to relative index in the meteo vertical. Note: omega wind is defined at 
!      ! half levels, i.e at the top of each grid box. See dq_omega in derived_field_quantitiues_2
!      !
!      fZC = arDyn(lp_z,iPart)  ! absolute pressure.
!      do iz = 1, nz_meteo    ! find the level index
!        if((fZC - a_half(iz-1) - b_half(iz-1) * ps) * (fZC - a_half(iz) - b_half(iz) * ps) <= 0.)exit
!      end do
!      fZC_rel = iz + 0.5 + (fZC - a_full(iz) - b_full(iz)*ps) / &
!                                   & (a_half(iz-1) - a_half(iz) + (b_half(iz-1) - b_half(iz)) * ps)
! arDyn index----------------------------------
      ! Get the integer index and absolute pressure
      !
      fZC_rel = arDyn(lp_z,iPart)
      iz = nint(fZC_rel)
if(iPart ==2)then
  call msg('Init, ix,iy',ix,iy)
  call msg('Init, i1d, ps', i1d, ps)
  call msg('Init, fZC_rel, iz', fZC_rel, iz)
endif

if(.not. (iz > 0 .and. iz <= nz_meteo))then
call msg('Funny iz for particle:'+fu_str(iPart), fZC_rel, iz)
call set_error('vertical index outside meteodata range','lagrangian_adv_3d_wind_endpoint')
return
endif

      pCentre = a_full(iz) + b_full(iz) * ps
      if(fZC_rel - iz > 0.)then  ! above the cell centre
        fZC = pCentre + 2. * (fZC_rel - iz) * (a_half(iz) + b_half(iz)*ps - pCentre)
      else
        fZC = pCentre - 2. * (fZC_rel - iz) * (a_half(iz-1) + b_half(iz-1)*ps - pCentre)
      endif

if(iPart ==2)then
  call msg('Init, pCentre,fZC',pCentre,fZC)
endif

if(fZC > ps)then
call msg('Funny fZC')
endif

#ifdef DEBUG_MORE
call msg('X-adv: z:',iz,fZC)
#endif 

      move_x = 0.0
      move_y = 0.0
      move_z = 0.0
      !
      ! Iterations go along 3 axes at once... remembering relative indices for horizontal
      ! and absolute coordinates for vertical position
      ! Move along z is absolute, each time has to find out the new coordinates
      !
      do iter = 1, niter
        !
        ! Inside iterations, push windpoints into the grid boundaries.
        !
if(iPart ==2)then
  call msg('Iter',iter)
endif
        wind_point_x = arDyn(lp_x,iPart) + 0.5 * move_x  ! relative x
        wind_point_y = arDyn(lp_y,iPart) + 0.5 * move_y  ! relative y
        iwp_x = max(1,min(nx_meteo,nint(wind_point_x)))  ! inside the grid
        iwp_y = max(1,min(ny_meteo,nint(wind_point_y)))
        i1d = iwp_x + (iwp_y-1)*nx_meteo

        ps = pMetBuf%p2d(sp_ind)%present%ptr(i1d)
        wind_point_z = min(ps-0.01, fZC + 0.5 * move_z)
        do iwp_z = 1, nz_meteo                                         ! find the level index
          if((wind_point_z - a_half(iwp_z-1) - b_half(iwp_z-1) * ps) * &
           & (wind_point_z - a_half(iwp_z) - b_half(iwp_z) * ps) <= 0.)exit
        end do
        pCentre = a_full(iwp_z) + b_full(iwp_z)*ps
        if(wind_point_z < pCentre)then
          wind_point_rel_z = iwp_z + 0.5 * (wind_point_z - pCentre) / &
                                               & (a_half(iwp_z) + b_half(iwp_z) * ps - pCentre)
        else
          wind_point_rel_z = iwp_z + 0.5 * (wind_point_z - pCentre) / &
                                               & (pCentre - a_half(iwp_z-1) - b_half(iwp_z-1) * ps)
        endif
        wind_point_rel_z = max(0.5, min(nz_meteo+0.5, wind_point_rel_z))
!        iwp_z = nint(wind_point_z)


        !
        ! Note that the wind is so far from the meteo buffer. Vertical is then the omega wind.
        ! Should be properly diagnosed, in principle.
        !
        wind_u = fu_4d_interpolation(pMetBuf%p4d(u_ind), &
                                   & wind_point_x, wind_point_y, wind_point_rel_z, &
                                   & nx_meteo, ny_meteo, nz_meteo, &
                                   & weight_past, &
                                   & horiz_interp_method, vert_interp_method, iOutsideAction)
        wind_v = fu_4d_interpolation(pMetBuf%p4d(v_ind), &
                                   & wind_point_x, wind_point_y, wind_point_rel_z, &
                                   & nx_meteo, ny_meteo, nz_meteo, &
                                   & weight_past, &
                                   & horiz_interp_method, vert_interp_method, iOutsideAction)
        if(w_ind == int_missing)then
          wind_w = 0.0
        else
          !
          ! Note that omega is defined at half-levels, i.e. at the interfaces between the boxes.
          ! I.e., particle with the coordinate 1.5 (at the interface btw 1 and 2 layer) exactly
          ! has wind = omega(z=1)
          !
          wind_w = fu_4d_interpolation(pMetBuf%p4d(w_ind), &
                                     & wind_point_x, wind_point_y, max(0.5,wind_point_rel_z-0.5), &
                                     & nx_meteo, ny_meteo, nz_meteo, &
                                     & weight_past, &
                                     & horiz_interp_method, vert_interp_method, iOutsideAction)
        endif
        !
        ! For the first level, reduce the wind. For a crude and dirty try, make it 
        ! linearly: logarithm would require computations in z-system
        !
if(iPart ==2)then
  call msg('Iter, windpoint x, wind u', wind_point_x, wind_u)
  call msg('Iter, windpoint y, wind v', wind_point_y, wind_v)
  call msg('Iter, windpoint z, wind omega', wind_point_z, wind_w)
endif

        if(wind_point_rel_z < 1.5)then
          wind_u = wind_u * (wind_point_rel_z -0.5)
          wind_v = wind_v * (wind_point_rel_z -0.5)
          wind_w = wind_w * (wind_point_rel_z -0.5)
if(iPart ==2)then
  call msg('Iter, scaling, wind u', wind_point_rel_z -0.5, wind_w)
endif
        endif
        move_x_old = move_x
        move_y_old = move_y
        move_z_old = move_z
        move_x = seconds * wind_u / pXCellSize(i1d)
        move_y = seconds * wind_v / pYCellSize(i1d)
        move_z = seconds * wind_w

if(iPart ==2)then
  call msg('Iter, move_old x, move_x', move_x_old, move_x)
  call msg('Iter, move_old y, move_y', move_y_old, move_y)
  call msg('Iter, move_old z, move_z', move_z_old, move_z)
endif

        if((move_x_old-move_x)**2 + (move_y_old-move_y)**2 + &
         & ((move_z_old-move_z) / (a_half(iz-1)-a_half(iz) + (b_half(iz-1)-b_half(iz))*ps))**2 < &
         & itertol**2) exit

      end do  ! move iterations

      !
      ! Horizontal move and checking the grid
      !
      arDyn(lp_x,iPart) = arDyn(lp_x,iPart) + move_x  ! relative
      arDyn(lp_y,iPart) = arDyn(lp_y,iPart) + move_y  ! relative
if(iPart ==2)then
  call msg('Move, new_x, new_y', arDyn(lp_x,iPart), arDyn(lp_y,iPart))
endif
      if(max((arDyn(lp_x,iPart)-0.5)*(arDyn(lp_x,iPart)-(nx_meteo+0.5)), &
           & (arDyn(lp_y,iPart)-0.5)*(arDyn(lp_y,iPart)-(ny_meteo+0.5))) >= 0.)then !Must not be on the edge
        call handle_out_of_grid(arDyn, iPart, iArBadParticle, iBadParticle, iOutsideAction) !, 'end_iter')
        cycle particle
      endif

      ix = nint(arDyn(lp_x,iPart))
      iy = nint(arDyn(lp_y,iPart))
      i1d = ix + (iy - 1) * nx_meteo
      ps = pMetBuf%p2d(sp_ind)%present%ptr(i1d)

      fZC = fZC + move_z ! absolute
if(iPart ==2)then
  call msg('Move, new_z, ps', fZC, ps)
endif
      if((fZC - ps)*(fZC - a_half(nz_meteo) - b_half(nz_meteo) * ps) > 0.)then
        call check_pressure_range(fZC, ps, a_half(nz_meteo) + b_half(nz_meteo) * ps, &
                                & arDyn, iPart, iArBadParticle, iBadParticle, 'end_iter')
        cycle particle
      endif

      do iz = 1, nz_meteo    ! find the level index
        if((fZC - a_half(iz-1) - b_half(iz-1) * ps) * (fZC - a_half(iz) - b_half(iz) * ps) <= 0.)exit
      end do
      pCentre = a_full(iz) + b_full(iz)*ps
      if(fZC < pCentre)then
        fZC_rel = iz + 0.5 * (fZC - pCentre) / (a_half(iz) + b_half(iz) * ps - pCentre)
      else
        fZC_rel = iz + 0.5 * (fZC - pCentre) / (pCentre - a_half(iz-1) - b_half(iz-1) * ps)
      endif
if(iPart ==2)then
  call msg('Move, new_z_rel, pCentre', fZC_rel, pCentre)
endif

!print *,'7'

      !
      ! Random-walk method selection
      !
      SELECT CASE (rw_method)

        case(void_rw)                        !============================================= VOID


        CASE(fixed_rw)                       !============================================= FIXED 

          arDyn(lp_x,iPart) = arDyn(lp_x,iPart) + &
                            & fu_random_number_center(0., 0.5) * seconds / pXCellSize(i1d) ! m/s
          arDyn(lp_y,iPart) = arDyn(lp_y,iPart) + &
                            & fu_random_number_center(0., 0.5) * seconds / pYCellSize(i1d) ! m/s
          fZC = min(ps-10., fZC + fu_random_number_center(0., (seconds * diffVelocityFT)))  ! Pa/s

        CASE(fully_mixed_abl_rw)             !======================================= FULLY-MIXED ABL
          !
          ! wind-dependent moves along the horizontal axes
          !
          fTmp = sqrt(wind_u*wind_u + wind_v*wind_v)
          arDyn(lp_x,iPart) = arDyn(lp_x,iPart) + &
                     & fu_random_number_center(0.,0.5) * a1 * (fTmp*seconds)**b / pXCellSize(i1d)
          arDyn(lp_y,iPart) = arDyn(lp_y,iPart) + &
                     & fu_random_number_center(0.,0.5) * a1 * (fTmp*seconds)**b / pYCellSize(i1d)
          !
          ! Vertical random walk: in pressure coordinates
          ! Inside ABL - mix-up, above it - fixed RW
          !
          if(fZC > pMetBuf%p2d(ind_HABL_Pa)%present%ptr(i1d))then
            fZC = fu_random_number_boundaries(pMetBuf%p2d(ind_SrfPress)%present%ptr(i1d) - 10., &
                                            & pMetBuf%p2d(ind_HABL_Pa)%present%ptr(i1d))
          else
            fZC = min(ps-10., fZC + fu_random_number_center(0., (seconds * diffVelocityFT)))
          endif

        CASE (bulk_gaussian_rw)                     !=============================== BULK GAUSSIAN

call set_error('Bulk gaussian random walk does not work yet','fu_random_walk')
return

        case default
          call set_error('Unknown random walk method:' + fu_str(rw_method),'lagrangian_adv_3d_wind_midpoint')
          return

      end select  ! Random-walk method
      !
      ! Absolute-pressure vertical position is the outcome of random-walk. Out-of-grid requires 
      ! relative index
      !
if(iPart ==2)then
  call msg('RW, new_x, new_y', arDyn(lp_x,iPart), arDyn(lp_y,iPart))
  call msg('RW, new_z, ps', fZC, ps)
endif
      if(max((arDyn(lp_x,iPart)-0.5)*(arDyn(lp_x,iPart)-(nx_meteo+0.5)), &
           & (arDyn(lp_y,iPart)-0.5)*(arDyn(lp_y,iPart)-(ny_meteo+0.5))) >= 0.)then
        call handle_out_of_grid(arDyn, iPart, iArBadParticle, iBadParticle, iOutsideAction) !, 'end')
        cycle particle
      endif
      if((fZC-ps)*(fZC-a_half(nz_meteo)-b_half(nz_meteo)*ps) > 0.)then
        call check_pressure_range(fZC, ps, a_half(nz_meteo)+b_half(nz_meteo)*ps, &
                                & arDyn, iPart, iArBadParticle, iBadParticle, 'end')
        cycle particle
      endif
      do iz = 1, nz_meteo    ! find the level index
        if((fZC - a_half(iz-1) - b_half(iz-1) * ps) * (fZC - a_half(iz) - b_half(iz) * ps) <= 0.)exit
      end do
      pCentre = a_full(iz) + b_full(iz)*ps
      if(fZC < pCentre)then
        arDyn(lp_z,iPart) = iz + 0.5 * (fZC - pCentre) / (a_half(iz) + b_half(iz) * ps - pCentre)
      else
        arDyn(lp_z,iPart) = iz + 0.5 * (fZC - pCentre) / (pCentre - a_half(iz-1) - b_half(iz-1) * ps)
      endif

if(iPart ==2)then
  call msg('RW, fZC_rel, pCentre', arDyn(lp_z,iPart), pCentre)
endif

!print *, 15
if(arDyn(lp_z,iPart) > nz_meteo + 0.5 .or. arDyn(lp_z,iPart) < 0.5)then
call msg('Funny z-coord:',arDyn(lp_z,iPart),iPart)
endif
      !Random walk could kick the particle  out of the grid
      if(max((arDyn(lp_x,iPart)-0.5)*(arDyn(lp_x,iPart)-(nx_meteo+0.5)), &
           & (arDyn(lp_y,iPart)-0.5)*(arDyn(lp_y,iPart)-(ny_meteo+0.5))) >= 0.)then !Must not be on the edge
        call handle_out_of_grid(arDyn, iPart, iArBadParticle, iBadParticle, iOutsideAction) !, 'end_iter')
        cycle particle
      endif
    end do particle  ! particles

  end subroutine lagrangian_adv_3d_wind_midpoint


  !*****************************************************************

  subroutine lagrangian_adv_2d_wind_midpoint(rw_method, &
                                           & arDyn, arMassTrn, nop, nSpecies, &
                                           & lpStatus, iArBadParticle, &
                                           & pMetBuf, & !pDispBuf, &
                                           & seconds, weight_past, &
                                           & wdr)
    !
    ! Advects all Lagrangian particles in the cloud. 
    ! Note that the Lagrangian particles are advected horizontally in the meteo grid realtive 
    ! coordinates, whereas vertical is taken in pressure coordinates.
    ! We also assume that the meteo vertical can be described as hybrid vertical: we use
    ! hybrid coefficients for fast computations.
    !
    implicit none

    ! Imported parameters
    INTEGER, INTENT(in) :: rw_method, nop, nSpecies
    real, dimension(:,:), intent(inout) :: arDyn   ! (9, nParticles)
    real, dimension(:,:), intent(inout) :: arMassTrn   ! (nSpecies, nParticles)
    type(Tfield_buffer), intent(in) :: pMetBuf !, pDispBuf
    real, intent(in) :: seconds, weight_past
    integer, dimension(:), intent(inout) :: lpStatus
    integer, dimension(:), intent(inout) :: iArBadParticle
    TYPE(silja_wdr), INTENT(in) :: wdr

    ! Local declarations:
    INTEGER, DIMENSION(:), POINTER :: q_arr_ptr
    real, dimension(:), pointer :: pXCellSize, pYCellSize, a_full, b_full, a_half, b_half
    INTEGER :: iPart, horiz_interp_method, vert_interp_method, iOutsideAction, ix, iy, iz, &
             & iTmp, p_ind, u_ind, v_ind, w_ind, sp_ind, h_ind, t_ind, grid_index, iter, i1d, &
             & iwp_x, iwp_y, iBadParticle, interp_q_2d, interp_q_3d
    REAL :: move_x, move_y, move_x_old, move_y_old, fTmp, wind_point_x, wind_point_y, wind_u, wind_v
    logical :: ifOutside
    logical, save :: ifFirst=.true.

    real, parameter :: itertol = 1e-2
    integer, parameter :: niter = 5
    !
    ! If we are here for the first time, store locally the vertical structure of meteo grid
    !
    if(ifFirst)then
      allocate(b_full(nz_meteo+1), a_full(nz_meteo+1), b_half(0:nz_meteo), a_half(0:nz_meteo), &
             & stat=iz)
      if(fu_fails(iz==0,'Failed hybrid coefs allocation','lagrangian_adv_2d_wind_midpoint'))return
      call hybrid_coefs(meteo_vertical, a_full=a_full, b_full=b_full)
      if(error)return
      a_half(0) = 0.
      b_half(0) = 1.
!      a_full(0) = 0.
!      b_full(0) = 2.
      a_full(nz_meteo+1) = a_full(nz_meteo) * (a_full(nz_meteo)/a_full(nz_meteo-1))**1.5
      b_full(nz_meteo+1) = max(0.0 , b_full(nz_meteo) + (b_full(nz_meteo) - b_full(nz_meteo-1)) * 1.5)
      do iz = 1, nz_meteo
        a_half(iz) = (a_full(iz) + a_full(iz+1)) / 2.
        b_half(iz) = (b_full(iz) + b_full(iz+1)) / 2.
      end do
      ifFirst = .false.
    endif  ! ifFirst

    !
    ! Get wind and pressure fields for past and future.
    ! Also prepare indices for needed fields
    !
    q_arr_ptr => pMetBuf%buffer_quantities
    p_ind = fu_index(q_arr_ptr, pressure_flag)
    h_ind = fu_index(q_arr_ptr, height_flag)
    t_ind = fu_index(q_arr_ptr, temperature_flag)
    sp_ind = fu_index(q_arr_ptr,surface_pressure_flag)
    u_ind = fu_index(q_arr_ptr, u_flag)
    v_ind = fu_index(q_arr_ptr, v_flag)
    !
    ! Instead of using 4D cube with pressure in meteo_vertical, we shall use the 
    ! hybrid coefficients and surface pressure. Easier to calculate, surface pressure is always
    ! available as "present".
    ! However, meteo_vertical is thin-layer structure, i.e. only full levels are defined.
    ! We shall have to use their half-sum to get the half-levels.
    !
    ! Below, two versions are presented: for arDyn keeping absolute pressure and relative coord 
    ! in meteo_vertical. Comment/uncomment them and related pieces in other modules to switch.
    !

!    ifPrintDump = .false. !now > tDebugTime !.and. now < tDebugTime + one_day * 2.

    !
    ! Set pointers in the random-walk routines
    !
!call msg('Random walk preparation starts')

    CALL prepare_random_walk(rw_method, pMetBuf)
    if(error)then
      call msg_warning('Random walk preparation failed','lagrangian_adv_2d_wind_midpoint')
      return
    endif

    horiz_interp_method = fu_horizontal_interp_method(wdr)
    vert_interp_method = linear ! So far.
    !
    !  Prepare pointers for the meteo grid cell sizes. Needed in advection
    ! in order to scale the absolute wind into the relative particle motion (and particles
    ! move in advection_grid to minimise reprojections).
    !
    pXCellSize => fu_grid_data(meteo_cell_x_size_fld)
    pYCellSize => fu_grid_data(meteo_cell_y_size_fld)

    if(fu_ifLonGlobal(meteo_grid))then
      iOutsideAction = handleGlobalGrid
    else
      iOutsideAction = nearestPoint !notAllowed
    endif

    if(error)then
      call set_error('Cannot start advection loop','lagrangian_adv_2d_wind_midpoint')
      return
    endif

    iBadParticle = 1
    iArBadParticle(iBadParticle) = int_missing

    !----------------------------------------
    !
    ! Grand loop over particles.
    !
    particle: do iPart = 1, nop
      !
      ! Skip the particle if it is not active
      !
      if(lpStatus(iPart) == int_missing)cycle

!print *, 'Cloud_advection_3'
!call stop_count(chCounterNm = strTmp)

      if(error) return

if(ifPrintDump) then
  call msg('adv x=', arDyn(lp_x,iPart))
  call msg('adv y=', arDyn(lp_y,iPart))
endif
      !
      ! Advection largely follows the Eulerian Galperin's scheme v.3 (implicit)
      ! except for the Eulerian step, which is not needed, of course. Also, it all goes in 
      ! meteo grid
      ! A complexity is: vertical axis is very inhomogeneous, so have to fly in 
      ! absolute pressure. It is also nice since random-walk is also in pressure to keep 
      ! const-mixing-ratio
      !  
      ix = nint(arDyn(lp_x,iPart))
      iy = nint(arDyn(lp_y,iPart))
      iz = nint(arDyn(lp_z,iPart))
      i1d = ix + (iy - 1) * nx_meteo

if(iPart ==2)then
  call msg('Init, ix,iy',ix,iy)
endif

      move_x = 0.0
      move_y = 0.0
      !
      ! Iterations go along 3 axes at once... remembering relative indices for horizontal
      ! and absolute coordinates for vertical position
      ! Move along z is absolute, each time has to find out the new coordinates
      !
      do iter = 1, niter
        !
        ! Inside iterations, push windpoints into the grid boundaries.
        !
if(iPart ==2)then
  call msg('Iter',iter)
endif
        wind_point_x = arDyn(lp_x,iPart) + 0.5 * move_x  ! relative x
        wind_point_y = arDyn(lp_y,iPart) + 0.5 * move_y  ! relative y
        iwp_x = max(1,min(nx_meteo,nint(wind_point_x)))  ! inside the grid
        iwp_y = max(1,min(ny_meteo,nint(wind_point_y)))
        i1d = iwp_x + (iwp_y-1)*nx_meteo

        !
        ! Note that the wind is so far from the meteo buffer. Vertical is then the omega wind.
        ! Should be properly diagnosed, in principle.
        !
        wind_u = fu_4d_interpolation(pMetBuf%p4d(u_ind), &
                                   & wind_point_x, wind_point_y, arDyn(lp_z,iPart), &
                                   & nx_meteo, ny_meteo, nz_meteo, &
                                   & weight_past, &
                                   & horiz_interp_method, vert_interp_method, iOutsideAction)
        wind_v = fu_4d_interpolation(pMetBuf%p4d(v_ind), &
                                   & wind_point_x, wind_point_y, arDyn(lp_z,iPart), &
                                   & nx_meteo, ny_meteo, nz_meteo, &
                                   & weight_past, &
                                   & horiz_interp_method, vert_interp_method, iOutsideAction)
        !
        ! For the first level, reduce the wind. For a crude and dirty try, make it 
        ! linearly: logarithm would require computations in z-system
        !
if(iPart ==2)then
  call msg('Iter, windpoint x, wind u', wind_point_x, wind_u)
  call msg('Iter, windpoint y, wind v', wind_point_y, wind_v)
endif

        if(arDyn(lp_z,iPart) < 1.5)then
          wind_u = wind_u * (arDyn(lp_z,iPart) -0.5)
          wind_v = wind_v * (arDyn(lp_z,iPart) -0.5)
        endif
        move_x_old = move_x
        move_y_old = move_y
        move_x = seconds * wind_u / pXCellSize(i1d)
        move_y = seconds * wind_v / pYCellSize(i1d)

if(iPart ==2)then
  call msg('Iter, move_old x, move_x', move_x_old, move_x)
  call msg('Iter, move_old y, move_y', move_y_old, move_y)
endif

        if((move_x_old-move_x)**2 + (move_y_old-move_y)**2 < itertol**2) exit

      end do  ! move iterations

      !
      ! Horizontal move and checking the grid
      !
      arDyn(lp_x,iPart) = arDyn(lp_x,iPart) + move_x  ! relative
      arDyn(lp_y,iPart) = arDyn(lp_y,iPart) + move_y  ! relative
if(iPart ==2)then
  call msg('Move, new_x, new_y', arDyn(lp_x,iPart), arDyn(lp_y,iPart))
endif
      if(max((arDyn(lp_x,iPart)-0.5)*(arDyn(lp_x,iPart)-(nx_meteo+0.5)), &
           & (arDyn(lp_y,iPart)-0.5)*(arDyn(lp_y,iPart)-(ny_meteo+0.5))) >= 0.)then
        call handle_out_of_grid(arDyn, iPart, iArBadParticle, iBadParticle, iOutsideAction) !, 'end_iter')
        cycle particle
      endif

      ix = nint(arDyn(lp_x,iPart))
      iy = nint(arDyn(lp_y,iPart))
      i1d = ix + (iy - 1) * nx_meteo

!print *,'7'

      !
      ! Random-walk method selection
      !
      SELECT CASE (rw_method)

        case(void_rw)                        !============================================= VOID


        CASE(fixed_rw)                       !============================================= FIXED 

          arDyn(lp_x,iPart) = arDyn(lp_x,iPart) + &
                            & fu_random_number_center(0., 0.5) * seconds / pXCellSize(i1d) ! m/s
          arDyn(lp_y,iPart) = arDyn(lp_y,iPart) + &
                            & fu_random_number_center(0., 0.5) * seconds / pYCellSize(i1d) ! m/s

        CASE(fully_mixed_abl_rw)             !======================================= FULLY-MIXED ABL
          !
          ! wind-dependent moves along the horizontal axes
          !
          fTmp = sqrt(wind_u*wind_u + wind_v*wind_v)
          arDyn(lp_x,iPart) = arDyn(lp_x,iPart) + &
                     & fu_random_number_center(0.,0.5) * a1 * (fTmp*seconds)**b / pXCellSize(i1d)
          arDyn(lp_y,iPart) = arDyn(lp_y,iPart) + &
                     & fu_random_number_center(0.,0.5) * a1 * (fTmp*seconds)**b / pYCellSize(i1d)

        CASE (bulk_gaussian_rw)                     !=============================== BULK GAUSSIAN

call set_error('Bulk gaussian random walk does not work yet','fu_random_walk')
return

        case default
          call set_error('Unknown random walk method:' + fu_str(rw_method),'lagrangian_adv_2d_wind_midpoint')
          return

      end select  ! Random-walk method
      !
      ! Absolute-pressure vertical position is the outcome of random-walk. Out-of-grid requires 
      ! relative index
      !
if(iPart ==2)then
  call msg('RW, new_x, new_y', arDyn(lp_x,iPart), arDyn(lp_y,iPart))
endif
      if(max((arDyn(lp_x,iPart)-0.5)*(arDyn(lp_x,iPart)-(nx_meteo+0.5)), &
           & (arDyn(lp_y,iPart)-0.5)*(arDyn(lp_y,iPart)-(ny_meteo+0.5))) >= 0.)then
        call handle_out_of_grid(arDyn, iPart, iArBadParticle, iBadParticle, iOutsideAction) !, 'end')
        cycle particle
      endif

    end do particle  ! particles

  end subroutine lagrangian_adv_2d_wind_midpoint


  ! ***************************************************************

  subroutine lagrangian_adv_3d_wind_endpoint(rw_method, &
                                           & arDyn, arMassTrn, nop, nSpecies, &
                                           & lpStatus, iArBadParticle, &
                                           & pMetBuf, & !pDispBuf, &
                                           & seconds, weight_past, &
                                           & wdr)
    !
    ! Advects all Lagrangian particles in the cloud. 
    ! Note that the Lagrangian particles are advected horizontally in the meteo grid realtive 
    ! coordinates, whereas vertical is taken in pressure coordinates.
    ! We also assume that the meteo vertical can be described as hybrid vertical: we use
    ! hybrid coefficients for fast computations.
    ! Here the wind iterations are based on start- and end- points of rhe motion vector.
    ! Should be somewhat safer that mid-point-based calculations.
    !
    implicit none

    ! Imported parameters
    INTEGER, INTENT(in) :: rw_method, nop, nSpecies
    real, dimension(:,:), intent(inout) :: arDyn   ! (9, nParticles)
    real, dimension(:,:), intent(inout) :: arMassTrn   ! (nSpecies, nParticles)
    type(Tfield_buffer), intent(in) :: pMetBuf !, pDispBuf
    real, intent(in) :: seconds, weight_past
    integer, dimension(:), intent(inout) :: lpStatus
    integer, dimension(:), intent(inout) :: iArBadParticle
    TYPE(silja_wdr), INTENT(in) :: wdr

    ! Local declarations:
    INTEGER, DIMENSION(:), POINTER :: q_arr_ptr
    real, dimension(:), pointer :: pXCellSize, pYCellSize
    real, dimension(:), pointer, save :: a_full, b_full, a_half, b_half
    INTEGER :: iPart, horiz_interp_method, vert_interp_method, iOutsideAction, ix, iy, iz, &
             & iTmp, p_ind, u_ind, v_ind, w_ind, sp_ind, h_ind, t_ind, grid_index, iter, i1d, &
             & iwp_x, iwp_y, iwp_z, iBadParticle, interp_q_2d, interp_q_3d
    REAL :: move_x, move_y, move_z, move_x_old, move_y_old, move_z_old, fZC, fZC_rel, fTmp, ps, &
          & wind_endpoint_x, wind_endpoint_y, wind_endpoint_z, wind_endpoint_rel_z, pCentre, &
          & wind_u_start, wind_v_start, wind_w_start, wind_u_end, wind_v_end, wind_w_end
real :: ps_ini, fZC_rel_ini
    logical :: ifOutside
    logical, save :: ifFirst=.true.

    real, parameter :: itertol = 1e-2
    integer, parameter :: niter = 5


    call msg("lagrangian_adv_3d_wind_endpoint")
    !
    ! If we are here for the first time, store locally the vertical structure of meteo grid
    !
    if(ifFirst)then
      allocate(b_full(nz_meteo+1), a_full(nz_meteo+1), b_half(0:nz_meteo), a_half(0:nz_meteo), &
             & stat=iz)
      if(fu_fails(iz==0,'Failed hybrid coefs allocation','lagrangian_adv_3d_wind_endpoint'))return
      call hybrid_coefs(meteo_vertical, a_full=a_full, b_full=b_full)
      if(error)return
      a_half(0) = 0.
      b_half(0) = 1.
!      a_full(0) = 0.
!      b_full(0) = 2.
      a_full(nz_meteo+1) = a_full(nz_meteo) * (a_full(nz_meteo)/a_full(nz_meteo-1))**1.5
      b_full(nz_meteo+1) = max(0.0 , b_full(nz_meteo) + (b_full(nz_meteo) - b_full(nz_meteo-1)) * 1.5)
      do iz = 1, nz_meteo
        a_half(iz) = (a_full(iz) + a_full(iz+1)) / 2.
        b_half(iz) = (b_full(iz) + b_full(iz+1)) / 2.
      end do
      ifFirst = .false.
    endif  ! ifFirst

    !
    ! Get wind and pressure fields for past and future.
    ! Also prepare indices for needed fields
    !
    q_arr_ptr => pMetBuf%buffer_quantities
    p_ind = fu_index(q_arr_ptr, pressure_flag)
    h_ind = fu_index(q_arr_ptr, height_flag)
    t_ind = fu_index(q_arr_ptr, temperature_flag)
    sp_ind = fu_index(q_arr_ptr,surface_pressure_flag)
    u_ind = fu_index(q_arr_ptr, u_flag)
    v_ind = fu_index(q_arr_ptr, v_flag)
    w_ind = fu_index(q_arr_ptr , omega_flag)
    !
    ! Instead of using 4D cube with pressure in meteo_vertical, we shall use the 
    ! hybrid coefficients and surface pressure. Easier to calculate, surface pressure is always
    ! available as "present".
    ! However, meteo_vertical is thin-layer structure, i.e. only full levels are defined.
    ! We shall have to use their half-sum to get the half-levels.
    !
    ! Below, two versions are presented: for arDyn keeping absolute pressure and relative coord 
    ! in meteo_vertical. Comment/uncomment them and related pieces in other modules to switch.
    !

!    ifPrintDump = .false. !now > tDebugTime !.and. now < tDebugTime + one_day * 2.

    !
    ! Set pointers in the random-walk routines
    !
!call msg('Random walk preparation starts')

    CALL prepare_random_walk(rw_method, pMetBuf)
    if(error)then
      call msg_warning('Random walk preparation failed','lagrangian_adv_3d_wind_endpoint')
      return
    endif

    horiz_interp_method = fu_horizontal_interp_method(wdr)
    vert_interp_method = linear ! So far.
    !
    !  Prepare pointers for the meteo grid cell sizes. Needed in advection
    ! in order to scale the absolute wind into the relative particle motion (and particles
    ! move in advection_grid to minimise reprojections).
    !
    pXCellSize => fu_grid_data(meteo_cell_x_size_fld)
    pYCellSize => fu_grid_data(meteo_cell_y_size_fld)

    if(fu_ifLonGlobal(meteo_grid))then
      iOutsideAction = handleGlobalGrid
    else
      iOutsideAction = nearestPoint !notAllowed
    endif

    if(error)then
      call set_error('Cannot start advection loop','lagrangian_adv_3d_wind_endpoint')
      return
    endif

    iBadParticle = 1
    iArBadParticle(iBadParticle) = int_missing

    !----------------------------------------
    !
    ! Grand loop over particles.
    !
    particle: do iPart = 1, nop
      !
      ! Skip the particle if it is not active
      !
      if(lpStatus(iPart) == int_missing)cycle

!print *, 'Cloud_advection_3'
!call stop_count(chCounterNm = strTmp)

      if(error) return

!if(ifPrintDump) then
!  call msg('adv x=', arDyn(lp_x,iPart))
!  call msg('adv y=', arDyn(lp_y,iPart))
!endif
      !
      ! Advection largely follows the Eulerian Galperin's scheme v.3 (implicit)
      ! except for the Eulerian step, which is not needed, of course. Also, it all goes in 
      ! meteo grid
      ! A complexity is: vertical axis is very inhomogeneous, so have to fly in 
      ! absolute pressure. It is also nice since random-walk is also in pressure to keep 
      ! const-mixing-ratio
      !  
      ix = nint(arDyn(lp_x,iPart))
      iy = nint(arDyn(lp_y,iPart))
      i1d = ix + (iy - 1) * nx_meteo
      ps = pMetBuf%p2d(sp_ind)%present%ptr(i1d)
ps_ini = ps
      !
      ! Get the integer index and absolute pressure
      !
      fZC_rel = arDyn(lp_z,iPart)
      iz = nint(fZC_rel)
fZC_rel_ini = fZC_rel
!if(iPart ==2)then
!  call msg('Init, ix,iy',ix,iy)
!  call msg('Init, i1d, ps', i1d, ps)
!  call msg('Init, fZC_rel, iz', fZC_rel, iz)
!endif

if(.not. (iz > 0 .and. iz <= nz_meteo))then
call msg('Funny iz for particle:'+fu_str(iPart), fZC_rel, iz)
call set_error('vertical index outside meteodata range','lagrangian_adv_3d_wind_endpoint')
return
endif

      wind_u_start = fu_4d_interpolation(pMetBuf%p4d(u_ind), &
                                       & arDyn(lp_x,iPart), arDyn(lp_y,iPart), arDyn(lp_z,iPart), &
                                       & nx_meteo, ny_meteo, nz_meteo, &
                                       & weight_past, &
                                       & horiz_interp_method, vert_interp_method, iOutsideAction)
      wind_v_start = fu_4d_interpolation(pMetBuf%p4d(v_ind), &
                                       & arDyn(lp_x,iPart), arDyn(lp_y,iPart), arDyn(lp_z,iPart), &
                                       & nx_meteo, ny_meteo, nz_meteo, &
                                       & weight_past, &
                                       & horiz_interp_method, vert_interp_method, iOutsideAction)
      !
      ! Note that omega is defined at half-levels, i.e. at the interfaces between the boxes.
      ! I.e., particle with the coordinate 1.5 (at the interface btw 1 and 2 layer) exactly
      ! has wind = omega(z=1)
      !
      wind_w_start = fu_4d_interpolation(pMetBuf%p4d(w_ind), &
                                       & arDyn(lp_x,iPart), arDyn(lp_y,iPart), max(0.5,arDyn(lp_z,iPart)-0.5), &
                                       & nx_meteo, ny_meteo, nz_meteo, &
                                       & weight_past, &
                                       & horiz_interp_method, vert_interp_method, iOutsideAction)

!if(iPart ==2)then
!  call msg('Iter, windpoint x, wind start u', arDyn(lp_x,iPart), wind_u_start)
!  call msg('Iter, windpoint y, wind start v', arDyn(lp_y,iPart), wind_v_start)
!  call msg('Iter, windpoint z, wind start omega', arDyn(lp_z,iPart), wind_w_start)
!endif

      if(fZC_rel < 1.)then
        wind_u_start = wind_u_start * 2. * (fZC_rel -0.5)
        wind_v_start = wind_v_start * 2. * (fZC_rel -0.5)
        wind_w_start = wind_w_start * 2. * (fZC_rel -0.5)
!if(iPart ==2)then
!  call msg('Iter, scaling, wind u end', wind_endpoint_z -0.5, wind_w_end)
!endif
      endif


      move_x = 0.0
      move_y = 0.0
      move_z = 0.0
      !
      ! Iterations go along 3 axes at once... remembering relative indices for horizontal
      ! and absolute coordinates for vertical position
      ! Move along z is absolute, each time has to find out the new coordinates
      !
      do iter = 1, niter
        !
        ! Inside iterations, push windpoints into the grid boundaries.
        !
!if(iPart ==2)then
!  call msg('Iter',iter)
!endif
        wind_endpoint_x = arDyn(lp_x,iPart) + move_x  ! relative x
        wind_endpoint_y = arDyn(lp_y,iPart) + move_y  ! relative y
        iwp_x = max(1,min(nx_meteo,nint(wind_endpoint_x)))  ! inside the grid
        iwp_y = max(1,min(ny_meteo,nint(wind_endpoint_y)))
        i1d = iwp_x + (iwp_y-1)*nx_meteo
        !
        ! Horizontal move means that surface pressure changed, i.e. absolute fZC at the end of move 
        ! is to be recomputed
        !
        ps = pMetBuf%p2d(sp_ind)%present%ptr(i1d)
        pCentre = a_full(iz) + b_full(iz) * ps
        if(fZC_rel - iz > 0.)then  ! above the cell centre
          fZC = pCentre + 2. * (fZC_rel - iz) * (a_half(iz) + b_half(iz)*ps - pCentre)
        else
          fZC = pCentre - 2. * (fZC_rel - iz) * (a_half(iz-1) + b_half(iz-1)*ps - pCentre)
        endif

        fZC = fZC + move_z      ! Now it is final position for the given iteration

        if(fZC > ps)then
          if(fZC - 0.5 * move_z > ps)then  ! zero wind at the surface => half the initial push
            call msg('Funny fZC end: wind pushes underground ps vs fZC,move_z:' + fu_str(ps), &
                   & fZC-move_z, move_z)
            call check_pressure_range(fZC, ps, a_half(nz_meteo) + b_half(nz_meteo) * ps, &
                                    & arDyn, iPart, iArBadParticle, iBadParticle, 'start_iter')
            cycle particle
          endif
          wind_endpoint_z = min(ps-0.01, fZC)  ! can continue: move seems to be safe
        else
          wind_endpoint_z = fZC
        endif

!if(iPart ==2)then
!  call msg('Init, pCentre,fZC',pCentre,fZC)
!endif

        !
        ! Find the wind_endpoint_z index
        !
        do iwp_z = 1, nz_meteo                                         ! find the level index
          if((wind_endpoint_z - a_half(iwp_z-1) - b_half(iwp_z-1) * ps) * &
           & (wind_endpoint_z - a_half(iwp_z) - b_half(iwp_z) * ps) <= 0.)exit
        end do
        pCentre = a_full(iwp_z) + b_full(iwp_z)*ps
        if(wind_endpoint_z < pCentre)then
          wind_endpoint_rel_z = iwp_z + 0.5 * (wind_endpoint_z - pCentre) / &
                                               & (a_half(iwp_z) + b_half(iwp_z) * ps - pCentre)
        else
          wind_endpoint_rel_z = iwp_z + 0.5 * (wind_endpoint_z - pCentre) / &
                                               & (pCentre - a_half(iwp_z-1) - b_half(iwp_z-1) * ps)
        endif
        wind_endpoint_rel_z = max(0.5, min(nz_meteo+0.5, wind_endpoint_rel_z))

        !
        ! Note that the wind is so far from the meteo buffer. Vertical is then the omega wind.
        ! Should be properly diagnosed, in principle.
        !
        wind_u_end = fu_4d_interpolation(pMetBuf%p4d(u_ind), &
                                       & wind_endpoint_x, wind_endpoint_y, wind_endpoint_rel_z, &
                                       & nx_meteo, ny_meteo, nz_meteo, &
                                       & weight_past, &
                                       & horiz_interp_method, vert_interp_method, iOutsideAction)
        wind_v_end = fu_4d_interpolation(pMetBuf%p4d(v_ind), &
                                       & wind_endpoint_x, wind_endpoint_y, wind_endpoint_rel_z, &
                                       & nx_meteo, ny_meteo, nz_meteo, &
                                       & weight_past, &
                                       & horiz_interp_method, vert_interp_method, iOutsideAction)
        !
        ! Note that omega is defined at half-levels, i.e. at the interfaces between the boxes.
        ! I.e., particle with the coordinate 1.5 (at the interface btw 1 and 2 layer) exactly
        ! has wind = omega(z=1)
        !
        wind_w_end = fu_4d_interpolation(pMetBuf%p4d(w_ind), &
                                       & wind_endpoint_x, wind_endpoint_y, &
                                                                & max(0.5,wind_endpoint_rel_z-0.5), &
                                       & nx_meteo, ny_meteo, nz_meteo, &
                                       & weight_past, &
                                       & horiz_interp_method, vert_interp_method, iOutsideAction)
        !
        ! For the first level, reduce the wind. For a crude and dirty try, make it 
        ! linearly: logarithm would require computations in z-system
        !
!if(iPart ==2)then
!  call msg('Iter, windpoint x, wind end u', wind_endpoint_x, wind_u_end)
!  call msg('Iter, windpoint y, wind end v', wind_endpoint_y, wind_v_end)
!  call msg('Iter, windpoint z, wind end omega', wind_endpoint_rel_z, wind_w_end)
!endif

        if(wind_endpoint_rel_z < 1.)then
          wind_u_end = wind_u_end * 2. * (wind_endpoint_rel_z -0.5)
          wind_v_end = wind_v_end * 2. * (wind_endpoint_rel_z -0.5)
          wind_w_end = wind_w_end * 2. * (wind_endpoint_rel_z -0.5)
!if(iPart ==2)then
!  call msg('Iter, scaling, wind u', wind_endpoint_rel_z -0.5, wind_w_end)
!endif
        endif
        move_x_old = move_x
        move_y_old = move_y
        move_z_old = move_z
        move_x = seconds * 0.5 * (wind_u_start + wind_u_end) / pXCellSize(i1d)
        move_y = seconds * 0.5 * (wind_v_start + wind_v_end) / pYCellSize(i1d)
        move_z = seconds * 0.5 * (wind_w_start + wind_w_end)

!if(iPart ==2)then
!  call msg('Iter, move_old x, move_x', move_x_old, move_x)
!  call msg('Iter, move_old y, move_y', move_y_old, move_y)
!  call msg('Iter, move_old z, move_z', move_z_old, move_z)
!endif

        if((move_x_old-move_x)**2 + (move_y_old-move_y)**2 + &
         & ((move_z_old-move_z) / (a_half(iz-1)-a_half(iz) + (b_half(iz-1)-b_half(iz))*ps))**2 < &
         & itertol**2) exit

      end do  ! move iterations

      !
      ! Horizontal move and checking the grid. All comes from the last iteration, 
      ! including fZC, ps, and pCentre
      !
      arDyn(lp_x,iPart) = wind_endpoint_x
      arDyn(lp_y,iPart) = wind_endpoint_y
      arDyn(lp_z,iPart) = wind_endpoint_rel_z
      ix = iwp_x
      iy = iwp_y
      iz = iwp_z
      fZC = wind_endpoint_z

!if(iPart ==31841)then
!  call msg('Move, new_x, new_y', arDyn(lp_x,iPart), arDyn(lp_y,iPart))
!  call msg('Move(int), new_x, new_y', nint(arDyn(lp_x,iPart)), nint(arDyn(lp_y,iPart)))
!endif
      if(max((arDyn(lp_x,iPart)-0.5)*(arDyn(lp_x,iPart)-(nx_meteo+0.5)), &
           & (arDyn(lp_y,iPart)-0.5)*(arDyn(lp_y,iPart)-(ny_meteo+0.5))) >= 0.)then
        call handle_out_of_grid(arDyn, iPart, iArBadParticle, iBadParticle, iOutsideAction) !, 'end_iter')
        cycle particle
      endif

!if(iPart ==2)then
!  call msg('Move, new_z, ps', fZC, ps)
!endif
      if((fZC - ps)*(fZC - a_half(nz_meteo) - b_half(nz_meteo) * ps) > 0.)then
        call check_pressure_range(fZC, ps, a_half(nz_meteo) + b_half(nz_meteo) * ps, &
                                & arDyn, iPart, iArBadParticle, iBadParticle, 'end_iter')
        cycle particle
      endif

!if(iPart ==2)then
!  call msg('Move, new_z_rel, pCentre', fZC_rel, pCentre)
!endif

!print *,'7'

      !
      ! Random-walk method selection
      !
      SELECT CASE (rw_method)

        case(void_rw)                        !============================================= VOID


        CASE(fixed_rw)                       !============================================= FIXED 

          arDyn(lp_x,iPart) = arDyn(lp_x,iPart) + &
                            & fu_random_number_center(0., 0.5) * abs(seconds) / pXCellSize(i1d) ! m/s
          arDyn(lp_y,iPart) = arDyn(lp_y,iPart) + &
                            & fu_random_number_center(0., 0.5) * abs(seconds) / pYCellSize(i1d) ! m/s
          fZC = min(ps-10., fZC + fu_random_number_center(0., (abs(seconds) * diffVelocityFT)))  ! Pa/s

        CASE(fully_mixed_abl_rw)             !======================================= FULLY-MIXED ABL
          !
          ! wind-dependent moves along the horizontal axes
          !
          fTmp = sqrt(wind_u_start * wind_u_start + wind_v_start * wind_v_start)
          arDyn(lp_x,iPart) = arDyn(lp_x,iPart) + &
                     & fu_random_number_center(0.,0.5) * a1 * (fTmp*abs(seconds))**b / pXCellSize(i1d)
          arDyn(lp_y,iPart) = arDyn(lp_y,iPart) + &
                     & fu_random_number_center(0.,0.5) * a1 * (fTmp*abs(seconds))**b / pYCellSize(i1d)
          !
          ! Vertical random walk: in pressure coordinates
          ! Inside ABL - mix-up, above it - fixed RW.
          ! Note that we do not put particles into the lowest 1m (~10Pa): turbulenc does not 
          ! reach there.
          !
          if(fZC > pMetBuf%p2d(ind_HABL_Pa)%present%ptr(i1d))then
            fZC = fu_random_number_boundaries(pMetBuf%p2d(ind_SrfPress)%present%ptr(i1d) - 10., &
                                            & pMetBuf%p2d(ind_HABL_Pa)%present%ptr(i1d))
          else
            fZC = min(ps-10., fZC + fu_random_number_center(0., (abs(seconds) * diffVelocityFT)))
          endif

        CASE (bulk_gaussian_rw)                     !=============================== BULK GAUSSIAN

call set_error('Bulk gaussian random walk does not work yet','lagrangian_adv_3d_wind_endpoint')
return

        case default
          call set_error('Unknown random walk method:' + fu_str(rw_method),'lagrangian_adv_3d_wind_endpoint')
          return

      end select  ! Random-walk method
      !
      ! Absolute-pressure vertical position is the outcome of random-walk. Out-of-grid requires 
      ! relative index
      !
!if(iPart == 31841)then
!  call msg('RW, new_x, new_y', arDyn(lp_x,iPart), arDyn(lp_y,iPart))
!  call msg('RW, new_x, new_y',nint(arDyn(lp_x,iPart)), nint(arDyn(lp_y,iPart)))
!endif
      if(max((arDyn(lp_x,iPart)-0.5)*(arDyn(lp_x,iPart)-(nx_meteo+0.5)), &
           & (arDyn(lp_y,iPart)-0.5)*(arDyn(lp_y,iPart)-(ny_meteo+0.5))) >= 0.)then
        call handle_out_of_grid(arDyn, iPart, iArBadParticle, iBadParticle, iOutsideAction) !, 'end')
        cycle particle
      endif
      if((fZC-ps)*(fZC-a_half(nz_meteo)-b_half(nz_meteo)*ps) > 0.)then
        call check_pressure_range(fZC, ps, a_half(nz_meteo)+b_half(nz_meteo)*ps, &
                                & arDyn, iPart, iArBadParticle, iBadParticle, 'end')
        cycle particle
      endif
      do iz = 1, nz_meteo    ! find the level index
        if((fZC - a_half(iz-1) - b_half(iz-1) * ps) * (fZC - a_half(iz) - b_half(iz) * ps) <= 0.)exit
      end do
      pCentre = a_full(iz) + b_full(iz)*ps
      if(fZC < pCentre)then
        arDyn(lp_z,iPart) = iz + 0.5 * (fZC - pCentre) / (a_half(iz) + b_half(iz) * ps - pCentre)
      else
        arDyn(lp_z,iPart) = iz + 0.5 * (fZC - pCentre) / (pCentre - a_half(iz-1) - b_half(iz-1) * ps)
      endif

!if(iPart ==2)then
!  call msg('RW, fZC_rel, pCentre', arDyn(lp_z,iPart), pCentre)
!endif

!print *, 15
!if(arDyn(lp_z,iPart) > nz_meteo + 0.5 .or. arDyn(lp_z,iPart) < 0.5)then
!call msg('Funny z-coord:',arDyn(lp_z,iPart),iPart)
!endif
    end do particle  ! particles

  end subroutine lagrangian_adv_3d_wind_endpoint


  ! ***************************************************************

  subroutine lagrangian_adv_3d_implicit(rw_method, &
                                      & arDyn, arMassTrn, nop, nSpecies, &
                                      & lpStatus, iArBadParticle, &
                                      & pMetBuf, & !pDispBuf, &
                                      & seconds, weight_past, &
                                      & wdr)
    !
    ! Advects all Lagrangian particles in the cloud. 
    ! Note that the Lagrangian particles are advected horizontally in the meteo grid realtive 
    ! coordinates, whereas vertical is taken in pressure coordinates.
    ! We also assume that the meteo vertical can be described as hybrid vertical: we use
    ! hybrid coefficients for fast computations.
    ! Here the wind iterations are based on start- and end- points of rhe motion vector.
    ! Should be somewhat safer that mid-point-based calculations.
    !
    implicit none

    ! Imported parameters
    INTEGER, INTENT(in) :: rw_method, nop, nSpecies
    real, dimension(:,:), intent(inout) :: arDyn   ! (9, nParticles)
    real, dimension(:,:), intent(inout) :: arMassTrn   ! (nSpecies, nParticles)
    type(Tfield_buffer), intent(in) :: pMetBuf !, pDispBuf
    real, intent(in) :: seconds, weight_past
    integer, dimension(:), intent(inout) :: lpStatus
    integer, dimension(:), intent(inout) :: iArBadParticle
    TYPE(silja_wdr), INTENT(in) :: wdr

    ! Local declarations:
    INTEGER, DIMENSION(:), POINTER :: q_arr_ptr
    real, dimension(:), pointer :: pXCellSize, pYCellSize
    real, dimension(:), pointer, save :: a_full, b_full, a_half, b_half
    INTEGER :: iPart, horiz_interp_method, vert_interp_method, iOutsideAction, ix, iy, iz, &
             & iTmp, p_ind, u_ind, v_ind, w_ind, sp_ind, h_ind, t_ind, grid_index, iter, i1d, &
             & iwp_x, iwp_y, iwp_z, iBadParticle, interp_q_2d, interp_q_3d, iThread, nThreads
    REAL :: move_x, move_y, move_z, move_x_old, move_y_old, move_z_old, fZC, fZC_rel, fTmp, ps, &
          & wind_endpoint_x, wind_endpoint_y, wind_endpoint_z, wind_endpoint_rel_z, pCentre, &
          & wind_u_end, wind_v_end, wind_w_end
    real :: ps_ini, fZC_rel_ini, fZC_ini
    logical :: ifOutside, ifVerticalSolved
    type(silja_ip_1d), dimension(:), pointer :: pArSetBadParticleTmp
    logical, save :: ifFirst=.true.

    real, parameter :: itertol = 1e-2
    integer, parameter :: niter = 5

    !
    ! If we are here for the first time, store locally the vertical structure of meteo grid
    !
    if(ifFirst)then
      allocate(b_full(nz_meteo+1), a_full(nz_meteo+1), b_half(0:nz_meteo), a_half(0:nz_meteo), &
             & stat=iz)
      if(fu_fails(iz==0,'Failed hybrid coefs allocation','lagrangian_adv_3d_implicit'))return
      call hybrid_coefs(meteo_vertical, a_full=a_full, b_full=b_full)
      if(error)return
      a_half(0) = 0.
      b_half(0) = 1.
!      a_full(0) = 0.
!      b_full(0) = 2.
      a_full(nz_meteo+1) = a_full(nz_meteo) * (a_full(nz_meteo)/a_full(nz_meteo-1))**1.5
      b_full(nz_meteo+1) = max(0.0 , b_full(nz_meteo) + (b_full(nz_meteo) - b_full(nz_meteo-1)) * 1.5)
      do iz = 1, nz_meteo
        a_half(iz) = (a_full(iz) + a_full(iz+1)) / 2.
        b_half(iz) = (b_full(iz) + b_full(iz+1)) / 2.
      end do
      ifFirst = .false.
    endif  ! ifFirst

    !
    ! Get wind and pressure fields for past and future.
    ! Also prepare indices for needed fields
    !
    q_arr_ptr => pMetBuf%buffer_quantities
    p_ind = fu_index(q_arr_ptr, pressure_flag)
    h_ind = fu_index(q_arr_ptr, height_flag)
    t_ind = fu_index(q_arr_ptr, temperature_flag)
    sp_ind = fu_index(q_arr_ptr,surface_pressure_flag)
    u_ind = fu_index(q_arr_ptr, u_flag)
    v_ind = fu_index(q_arr_ptr, v_flag)
    w_ind = fu_index(q_arr_ptr , omega_flag)
    !
    ! Instead of using 4D cube with pressure in meteo_vertical, we shall use the 
    ! hybrid coefficients and surface pressure. Easier to calculate, surface pressure is always
    ! available as "present".
    ! However, meteo_vertical is thin-layer structure, i.e. only full levels are defined.
    ! We shall have to use their half-sum to get the half-levels.
    !
    ! Below, two versions are presented: for arDyn keeping absolute pressure and relative coord 
    ! in meteo_vertical. Comment/uncomment them and related pieces in other modules to switch.
    !

!    ifPrintDump = .false. !now > tDebugTime !.and. now < tDebugTime + one_day * 2.

    !
    ! Set pointers in the random-walk routines
    !
!call msg('Random walk preparation starts')

    CALL prepare_random_walk(rw_method, pMetBuf)
    if(error)then
      call msg_warning('Random walk preparation failed','lagrangian_adv_3d_implicit')
      return
    endif

    horiz_interp_method = fu_horizontal_interp_method(wdr)
    vert_interp_method = linear ! So far.
    !
    !  Prepare pointers for the meteo grid cell sizes. Needed in advection
    ! in order to scale the absolute wind into the relative particle motion (and particles
    ! move in advection_grid to minimise reprojections).
    !
    pXCellSize => fu_grid_data(meteo_cell_x_size_fld)
    pYCellSize => fu_grid_data(meteo_cell_y_size_fld)

    if(fu_ifLonGlobal(meteo_grid))then
      iOutsideAction = handleGlobalGrid
    else
      iOutsideAction = nearestPoint !notAllowed
    endif

    if(error)then
      call set_error('Cannot start advection loop','lagrangian_adv_3d_implicit')
      return
    endif
!call msg('ad1')

    !$OMP PARALLEL DEFAULT(NONE) &
    !$OMP & PRIVATE(iPart, ix, iy, iz, i1d, ps, ps_ini, fZC, fZC_ini, fZC_rel, fZC_rel_ini, iter, &
    !$OMP         & move_x, move_y, move_z, move_x_old, move_y_old, move_z_old, ifVerticalSolved, &
    !$OMP         & wind_u_end, wind_v_end, wind_w_end, fTmp, iBadParticle, iThread, &
    !$OMP         & wind_endpoint_x, wind_endpoint_y, wind_endpoint_z, wind_endpoint_rel_z, &
    !$OMP         & iwp_x, iwp_y, pCentre) &
    !$OMP & SHARED(rw_method, nop, lpStatus, error, arDyn, nx_meteo, ny_meteo, nz_meteo, nThreads, &
    !$OMP         & pMetBuf, sp_ind, a_full, b_full, a_half, b_half, weight_past, u_ind, v_ind, w_ind, &
    !$OMP         & horiz_interp_method, vert_interp_method, iOutsideAction, iArBadParticle, &
    !$OMP         & ind_HAbl_Pa, ind_srfPress, seconds, pXCellSize, pYCellSize, pArSetBadParticleTmp)

    
    iBadParticle = 1

!call msg('ad2')
    !$OMP MASTER
    nThreads = 1
    !$ nThreads = omp_get_num_threads()
     Call msg('Doing lagrangian_adv_3d_implicit, nthreads:', nThreads)
!    call msg('nThreads:',nThreads, nop)
!    call msg('Get array')
   !!FIXME Allocation is wrong if non-static schedule is used
    call get_work_arrays_set(nThreads, nop / nThreads + 2 , pArSetBadParticleTmp)

!     call msg('Got array', size(pArSetBadParticleTmp),size(pArSetBadParticleTmp(1)%ip))
    if(.not. error)then
      do iTmp = 1, nThreads
!        call msg('iTmp',iTmp)
        pArSetBadParticleTmp(iTmp)%ip(1) = int_missing
      end do
    end if
    !$OMP END MASTER
!call msg('ad3')
    !$OMP BARRIER
    
    iThread = 0
    !$ iThread = OMP_GET_THREAD_NUM()

    !----------------------------------------
    !
    ! Grand loop over particles.
    !
    !$OMP DO
    
    particle: do iPart = 1, nop
      !
      ! Skip the particle if it is not active
      !
      if(lpStatus(iPart) == int_missing)cycle

!call msg('ad4')
!print *, 'Cloud_advection_3'
!call stop_count(chCounterNm = strTmp)

      if(error)cycle

!if(ifPrintDump) then
!  call msg('adv x=', arDyn(lp_x,iPart))
!  call msg('adv y=', arDyn(lp_y,iPart))
!endif
      !
      ! Advection largely follows the Eulerian Galperin's scheme v.3 (implicit)
      ! except for the Eulerian step, which is not needed, of course. Also, it all goes in 
      ! meteo grid
      ! A complexity is: vertical axis is very inhomogeneous, so have to fly in 
      ! absolute pressure. It is also nice since random-walk is also in pressure to keep 
      ! const-mixing-ratio
      !  
      ix = nint(arDyn(lp_x,iPart))
      iy = nint(arDyn(lp_y,iPart))
      i1d = ix + (iy - 1) * nx_meteo
      ps = pMetBuf%p2d(sp_ind)%present%ptr(i1d)
ps_ini = ps
      !
      ! Get the integer index and absolute pressure
      !
      fZC_rel = arDyn(lp_z,iPart)
      iz = nint(fZC_rel)
fZC_rel_ini = fZC_rel
!if(iPart ==2)then
!  call msg('Init, ix,iy',ix,iy)
!  call msg('Init, i1d, ps', i1d, ps)
!  call msg('Init, fZC_rel, iz', fZC_rel, iz)
!endif

if(.not. (iz > 0 .and. iz <= nz_meteo))then
call msg('Funny iz for particle:'+fu_str(iPart), fZC_rel, iz)
call set_error('vertical index outside meteodata range','lagrangian_adv_3d_implicit')
cycle
endif

!call msg('ad5')
      move_x = 0.0
      move_y = 0.0
      move_z = 0.0
      ifVerticalSolved = .false.
      !
      ! Iterations go along 3 axes at once... remembering relative indices for horizontal
      ! and absolute coordinates for vertical position
      ! Move along z is absolute, each time has to find out the new coordinates
      !
      do iter = 1, niter
        !
        ! Inside iterations, push windpoints into the grid boundaries.
        !
!if(iPart ==2)then
!  call msg('Iter',iter)
!endif
        wind_endpoint_x = arDyn(lp_x,iPart) + move_x  ! relative x
        wind_endpoint_y = arDyn(lp_y,iPart) + move_y  ! relative y
        iwp_x = max(1,min(nx_meteo,nint(wind_endpoint_x)))  ! inside the grid
        iwp_y = max(1,min(ny_meteo,nint(wind_endpoint_y)))
        i1d = iwp_x + (iwp_y-1)*nx_meteo
        !
        ! Horizontal move means that surface pressure changed, i.e. absolute fZC at the end of move 
        ! is to be recomputed
        !
        ps = pMetBuf%p2d(sp_ind)%present%ptr(i1d)
        pCentre = a_full(iz) + b_full(iz) * ps
        if(.not. ifVerticalSolved)then
          if(fZC_rel - iz > 0.)then  ! above the cell centre
            fZC = pCentre + 2. * (fZC_rel - iz) * (a_half(iz) + b_half(iz)*ps - pCentre)
          else
            fZC = pCentre - 2. * (fZC_rel - iz) * (a_half(iz-1) + b_half(iz-1)*ps - pCentre)
          endif

          if(fZC + move_z > ps)then
            fZC_ini = fZC
            fZC = ps - (ps-fZC_ini)*exp(-move_z/(ps-fZC_ini))  ! this is an exact solution
            move_z = fZC - fZC_ini
            move_z_old = move_z
            ifVerticalSolved = .true.
          else
            fZC = fZC + move_z
          endif
          wind_endpoint_z = fZC
          !
          ! Find the wind_endpoint_z index
          !
          do iwp_z = 1, nz_meteo                                         ! find the level index
            if((wind_endpoint_z - a_half(iwp_z-1) - b_half(iwp_z-1) * ps) * &
             & (wind_endpoint_z - a_half(iwp_z) - b_half(iwp_z) * ps) <= 0.)exit
          end do
          iwp_z = min(nz_meteo, iwp_z)   ! for the case if we gone outside
          pCentre = a_full(iwp_z) + b_full(iwp_z)*ps
          if(wind_endpoint_z < pCentre)then
            wind_endpoint_rel_z = iwp_z + 0.5 * (wind_endpoint_z - pCentre) / &
                                                 & (a_half(iwp_z) + b_half(iwp_z) * ps - pCentre)
          else
            wind_endpoint_rel_z = iwp_z + 0.5 * (wind_endpoint_z - pCentre) / &
                                                 & (pCentre - a_half(iwp_z-1) - b_half(iwp_z-1) * ps)
          endif
          wind_endpoint_rel_z = max(0.5001, min(nz_meteo+0.4999, wind_endpoint_rel_z))
        endif  ! vertical solved

        !
        ! If vertical movement solved to the end, should leave here: iterations along 
        ! horizontal dimension may change surface pressure and make the whole thing inconsistent.
        !
        if(ifVerticalSolved)exit

        wind_u_end = fu_4d_interpolation(pMetBuf%p4d(u_ind), &
                                       & wind_endpoint_x, wind_endpoint_y, wind_endpoint_rel_z, &
                                       & nx_meteo, ny_meteo, nz_meteo, &
                                       & weight_past, &
                                       & horiz_interp_method, vert_interp_method, iOutsideAction)
        wind_v_end = fu_4d_interpolation(pMetBuf%p4d(v_ind), &
                                       & wind_endpoint_x, wind_endpoint_y, wind_endpoint_rel_z, &
                                       & nx_meteo, ny_meteo, nz_meteo, &
                                       & weight_past, &
                                       & horiz_interp_method, vert_interp_method, iOutsideAction)
        !
        ! Note that omega is defined at half-levels, i.e. at the interfaces between the boxes.
        ! I.e., particle with the coordinate 1.5 (at the interface btw 1 and 2 layer) exactly
        ! has wind = omega(z=1).
        ! Also, do not waste time if the problem is solved till the end
        !
        wind_w_end = fu_4d_interpolation(pMetBuf%p4d(w_ind), &
                                       & min(nx_meteo+0.4999,max(0.5001,wind_endpoint_x)), &
                                       & min(ny_meteo+0.4999,max(0.5001,wind_endpoint_y)), &
                                       & min(nz_meteo+0.4999,max(0.5001,wind_endpoint_rel_z-0.5)), &
                                       & nx_meteo, ny_meteo, nz_meteo, &
                                       & weight_past, &
                                       & horiz_interp_method, vert_interp_method, iOutsideAction)
        if(wind_endpoint_rel_z < 1.)then
          wind_u_end = wind_u_end * 2. * (wind_endpoint_rel_z -0.5)
          wind_v_end = wind_v_end * 2. * (wind_endpoint_rel_z -0.5)
          wind_w_end = wind_w_end * 2. * (wind_endpoint_rel_z -0.5)
        endif
        move_x_old = move_x
        move_y_old = move_y
        move_z_old = move_z
        move_x = seconds * wind_u_end / pXCellSize(i1d)
        move_y = seconds * wind_v_end / pYCellSize(i1d)
        move_z = seconds * wind_w_end

        if((move_x_old-move_x)**2 + (move_y_old-move_y)**2 + &
         & ((move_z_old-move_z) / (a_half(iz-1)-a_half(iz) + (b_half(iz-1)-b_half(iz))*ps))**2 < &
         & itertol**2) exit

        !
        ! For the first level, reduce the wind. For a crude and dirty try, make it 
        ! linearly: logarithm would require computations in z-system
        !
!if(iPart ==2)then
!  call msg('Iter, windpoint x, wind end u', wind_endpoint_x, wind_u_end)
!  call msg('Iter, windpoint y, wind end v', wind_endpoint_y, wind_v_end)
!  call msg('Iter, windpoint z, wind end omega', wind_endpoint_rel_z, wind_w_end)
!endif

      end do  ! move iterations

      !
      ! Horizontal move and checking the grid. All comes from the last iteration, 
      ! including fZC, ps, and pCentre
      !
      arDyn(lp_x,iPart) = wind_endpoint_x
      arDyn(lp_y,iPart) = wind_endpoint_y
      arDyn(lp_z,iPart) = wind_endpoint_rel_z
      ix = iwp_x
      iy = iwp_y
      iz = iwp_z
      fZC = wind_endpoint_z

!if(iPart ==2)then
!  call msg('Move, new_x, new_y', arDyn(lp_x,iPart), arDyn(lp_y,iPart))
!endif
      if(max((arDyn(lp_x,iPart)-0.5)*(arDyn(lp_x,iPart)-(nx_meteo+0.5)), &
           & (arDyn(lp_y,iPart)-0.5)*(arDyn(lp_y,iPart)-(ny_meteo+0.5))) >= 0.)then
!        !$OMP CRITICAL
        call handle_out_of_grid(arDyn, iPart, pArSetBadParticleTmp(iThread+1)%ip, iBadParticle, iOutsideAction) !, 'end_iter')
!call msg('Bad particle:'+fu_str(iPart)+','+fu_str(iThread)+':',pArSetBadParticleTmp(iThread+1)%ip(1:10))
!        !$OMP END CRITICAL
        cycle particle
      endif

!if(iPart ==2)then
!  call msg('Move, new_z, ps', fZC, ps)
!endif
      if((fZC - ps)*(fZC - a_half(nz_meteo) - b_half(nz_meteo) * ps) > 0.)then
!        !$OMP CRITICAL
        call check_pressure_range(fZC, ps, a_half(nz_meteo) + b_half(nz_meteo) * ps, &
                                & arDyn, iPart, pArSetBadParticleTmp(iThread+1)%ip, iBadParticle) !, 'end_iter')
!call msg('Bad particle:'+fu_str(iPart)+','+fu_str(iThread)+':',pArSetBadParticleTmp(iThread+1)%ip(1:10))
!        !$OMP END CRITICAL
        cycle particle
      endif

!if(iPart ==2)then
!  call msg('Move, new_z_rel, pCentre', fZC_rel, pCentre)
!endif

!print *,'7'

      !
      ! Random-walk method selection
      !
      SELECT CASE (rw_method)

        case(void_rw)                        !============================================= VOID


        CASE(fixed_rw)                       !============================================= FIXED 

          arDyn(lp_x,iPart) = arDyn(lp_x,iPart) + &
                            & fu_random_number_center(0., 0.5) * abs(seconds) / pXCellSize(i1d) ! m/s
          arDyn(lp_y,iPart) = arDyn(lp_y,iPart) + &
                            & fu_random_number_center(0., 0.5) * abs(seconds) / pYCellSize(i1d) ! m/s
          fZC = fZC + fu_random_number_center(0., (abs(seconds) * diffVelocityFT))  ! Pa/s

        CASE(fully_mixed_abl_rw)             !======================================= FULLY-MIXED ABL
          !
          ! wind-dependent moves along the horizontal axes
          !
          fTmp = sqrt(wind_u_end * wind_u_end + wind_v_end * wind_v_end)
          arDyn(lp_x,iPart) = arDyn(lp_x,iPart) + &
                     & fu_random_number_center(0.,0.5) * a1 * (fTmp*abs(seconds))**b / pXCellSize(i1d)
          arDyn(lp_y,iPart) = arDyn(lp_y,iPart) + &
                     & fu_random_number_center(0.,0.5) * a1 * (fTmp*abs(seconds))**b / pYCellSize(i1d)
          !
          ! Vertical random walk: in pressure coordinates
          ! Inside ABL - mix-up, above it - fixed RW.
          ! Note that we do not put particles into the lowest 1m (~10Pa): turbulenc does not 
          ! reach there.
          !
          fZC_ini = fZC
          if(fZC > pMetBuf%p2d(ind_HABL_Pa)%present%ptr(i1d))then
            fZC = fu_random_number_boundaries(pMetBuf%p2d(ind_SrfPress)%present%ptr(i1d) - 10., &
                                            & pMetBuf%p2d(ind_HABL_Pa)%present%ptr(i1d))
          else
            fZC = min(ps-10.,fZC + fu_random_number_center(0., (abs(seconds) * diffVelocityFT)))
          endif

        CASE (bulk_gaussian_rw)                     !=============================== BULK GAUSSIAN

          call set_error('Bulk gaussian random walk does not work yet','lagrangian_adv_3d_implicit')
          cycle

        case default
          call set_error('Unknown random walk method:' + fu_str(rw_method),'lagrangian_adv_3d_implicit')
          cycle

      end select  ! Random-walk method
      !
      ! Absolute-pressure vertical position is the outcome of random-walk. Out-of-grid requires 
      ! relative index
      !
!if(iPart ==2)then
!  call msg('RW, new_x, new_y', arDyn(lp_x,iPart), arDyn(lp_y,iPart))
!  call msg('RW, new_z, ps', fZC, ps)
!endif
      if(max((arDyn(lp_x,iPart)-0.5)*(arDyn(lp_x,iPart)-(nx_meteo+0.5)), &
           & (arDyn(lp_y,iPart)-0.5)*(arDyn(lp_y,iPart)-(ny_meteo+0.5))) >= 0.)then
        call handle_out_of_grid(arDyn, iPart, pArSetBadParticleTmp(iThread+1)%ip, iBadParticle, iOutsideAction) !, 'end')
        cycle particle
      endif
      if((fZC-ps)*(fZC-a_half(nz_meteo)-b_half(nz_meteo)*ps) > 0.)then
        call check_pressure_range(fZC, ps, a_half(nz_meteo)+b_half(nz_meteo)*ps, &
                                & arDyn, iPart, pArSetBadParticleTmp(iThread+1)%ip, iBadParticle) !, 'end')
        cycle particle
      endif
      do iz = 1, nz_meteo    ! find the level index
        if((fZC - a_half(iz-1) - b_half(iz-1) * ps) * (fZC - a_half(iz) - b_half(iz) * ps) <= 0.)exit
      end do
      pCentre = a_full(iz) + b_full(iz)*ps
      if(fZC < pCentre)then
        arDyn(lp_z,iPart) = iz + 0.5 * (fZC - pCentre) / (a_half(iz) + b_half(iz) * ps - pCentre)
      else
        arDyn(lp_z,iPart) = iz + 0.5 * (fZC - pCentre) / (pCentre - a_half(iz-1) - b_half(iz-1) * ps)
      endif
      
if(arDyn(lp_z,iPart) > nz_meteo + 0.5 .or. arDyn(lp_z,iPart) < 0.5)then
  call msg('Funny z-coord:',arDyn(lp_z,iPart),iPart)
else
  arDyn(lp_z,iPart) = max(0.5001,min(nz_meteo+0.4999,arDyn(lp_z,iPart)))
endif
    end do particle  ! particles

    !$OMP END DO

!call msg('ad4')
    !$OMP END PARALLEL
    
    ix = 1
    do iTmp = 1, nThreads
      iy = 1
!      call msg('Thread and bad particles:' + fu_str(iTmp), pArSetBadParticleTmp(iTmp)%ip(1:10))
      do while(pArSetBadParticleTmp(iTmp)%ip(iy) /= int_missing)
        iArBadParticle(ix) = pArSetBadParticleTmp(iTmp)%ip(iy) 
!        call msg('Bad particle & thread:',iTmp, pArSetBadParticleTmp(iTmp)%ip(iy))
        iy = iy + 1
        ix = ix + 1
      end do
    end do
    iArBadParticle(ix) = int_missing
!    call msg('Bad particles:',iArBadParticle(1:ix))
    call free_work_array(pArSetBadParticleTmp)


  end subroutine lagrangian_adv_3d_implicit


  !********************************************************************************** 

  subroutine handle_out_of_grid(arDyn, iPart, iArBadPart, iBadPart, iOutsideAction,  chPlace)
      !
      ! handles out-of-grid for the given x,y,z. Interestingly, no matter what these are - wind points
      ! or actual particle positions, the result is the same: particle is sent to out-of-grid direction
      !
      implicit none
      integer, intent(in) :: iOutsideAction, iPart
      real, dimension(:,:), intent(inout) :: arDyn
      integer, dimension(:), intent(inout) :: iArBadPart
      integer, intent(inout) :: iBadPart
      character(len=*), intent(in), optional :: chPlace
      logical :: ifHandled, ifBad

      ifHandled=.true.
      ifBad = .True.
      if(arDyn(lp_x,iPart) <= 0.5)then
        if(iOutsideAction == handleGlobalGrid)then   ! global ? Then over-the-x-border jump turns the grid
          arDyn(lp_x,iPart) = arDyn(lp_x,iPart) + nx_meteo
          ifBad = .False.
          if(present(chPlace))call msg(chPlace + ': +nx_meteo', iPart)
        else
          arDyn(lp_x,iPart) = -1
          if(present(chPlace))call msg(chPlace + ': x <= 0.5', iPart)
        endif
      elseif(arDyn(lp_x,iPart) >= nx_meteo + 0.5)then
        if(iOutsideAction == handleGlobalGrid)then   ! global ? Then over-the-x-border jump turns the grid
          arDyn(lp_x,iPart) = arDyn(lp_x,iPart) - nx_meteo
          ifBad = .False.
          if(present(chPlace))call msg(chPlace + ': -nx_meteo', iPart)
        else
          arDyn(lp_x,iPart) = nx_meteo + 1
          if(present(chPlace))call msg(chPlace + ': x >= nx_meteo', iPart)
        endif
      else
            ifHandled=.false.
        endif
      ! handle out-of-grid case
      !
      ! y-out-of-grid
      !
      if(arDyn(lp_y,iPart) <= 0.5)then
        arDyn(lp_y,iPart) = -1
        ifHandled =  .true.
        if(present(chPlace))call msg(chPlace + ': y <= 0.5', iPart)
      elseif(arDyn(lp_y,iPart) >= ny_meteo + 0.5)then
        arDyn(lp_y,iPart) = ny_meteo + 1
        ifHandled =  .true.
        if(present(chPlace))call msg(chPlace + ': y >= ny_meteo', iPart)
      endif  ! handle out-of-grid case

      if (ifHandled) then
          if (ifBad) then
        iArBadPart(iBadPart) = iPart
        iArBadPart(iBadPart+1) = int_missing
        iBadPart = iBadPart + 1
           endif
      else
           
           call msg('Particle No'+fu_str(iPart)+': xy corrds', arDyn(lp_x,iPart), arDyn(lp_y,iPart))
           call msg('nx_meteo, ny_meteo', nx_meteo, ny_meteo)
           call set_error("Called for in-grid particle", "handle_out_of_grid" )
      endif

    end subroutine handle_out_of_grid

    !**************************************************************
    
    subroutine check_pressure_range(fZC, pSrf, pUp, arDyn, iPart, iArBadPart, iBadPart, chPlace)
      implicit none
      real, intent(in) :: fZC, pSrf, pUp
      real, dimension(:,:), intent(inout) :: arDyn
      integer, dimension(:),  intent(inout)  :: iArBadPart
      integer, intent(inout) :: iBadPart
      integer, intent(in) :: iPart
      character(len=*), intent(in), optional :: chPlace

      if(fZC > pSrf)then      ! check the pressure range.
        arDyn(lp_z,iPart) = -1
        iArBadPart(iBadPart) = iPart
        iArBadPart(iBadPart+1) = int_missing
        iBadPart = iBadPart + 1
if(present(chPlace))call msg(chPlace + ': z < 1', iPart)
      elseif(fZC < pUp)then
        arDyn(lp_z,iPart) = nz_meteo + 1
        iArBadPart(iBadPart) = iPart
        iArBadPart(iBadPart+1) = int_missing
        iBadPart = iBadPart + 1
if(present(chPlace))call msg(chPlace + ': z > nz_meteo', iPart)
      endif
    end subroutine check_pressure_range



  ! ****************************************************************
  ! ****************************************************************
  ! ****************************************************************
  !
  !  Random-walk methods start
  !
  ! ****************************************************************
  ! ****************************************************************
  ! ****************************************************************

  SUBROUTINE prepare_random_walk(rw_method, met_buf)
    !
    ! Makes a preparatory steps for the random-walk process. Actually
    ! just sets the pointers to appropriate meteofields. 
    !
    IMPLICIT NONE

    ! Imported parameters
    INTEGER, INTENT(in) :: rw_method
    TYPE(Tfield_buffer), intent(in), target :: met_buf

    ! Local declarations
    INTEGER :: iTmp
    INTEGER, DIMENSION(:), POINTER :: mdl_in_q

    !
    ! Preparation depends on the random-walk method:
    !
    mdl_in_q => met_buf%buffer_quantities

    SELECT CASE(rw_method)
    CASE (void_rw, fixed_rw) ! Do nothing - it does not need any preparations

    CASE (fully_mixed_abl_rw)
      !
      ! ATTENTION! Random-walk is in pressure coord, always
      !
      ind_HABL_Pa = fu_index(mdl_in_q, abl_top_pressure_flag)
      ind_SrfPress = fu_index(mdl_in_q, surface_pressure_flag)

    CASE (bulk_gaussian_rw)

call set_error('Bulk gaussian advection does not work','prepare_random_walk')
return
      iTmp = fu_index(mdl_in_q, temperature_flag)
      temp_4d_ptr => met_buf%p4d(iTmp)

      iTmp = fu_index(mdl_in_q, pressure_flag)
      pr_4d_ptr => met_buf%p4d(iTmp)

      iTmp = fu_index(mdl_in_q, height_flag)
      height_4d_ptr => met_buf%p4d(iTmp)

!      iTmp = fu_index(mdl_in_q, surface_pressure_flag)
!      surf_pr_ptr => met_buf%p2d(iTmp)  !%present%ptr

!      iTmp = fu_index(mdl_in_q, abl_top_pressure_flag)
!      abl_h_pa_ptr => met_buf%p2d(iTmp)  !%present%ptr

      iTmp = fu_index(mdl_in_q, abl_height_m_flag)
      abl_h_m_ptr => met_buf%p2d(iTmp)%present%ptr

      iTmp = fu_index(mdl_in_q, MO_length_inv_flag)
      MO_len_inv_ptr => met_buf%p2d(iTmp)%present%ptr

      iTmp = fu_index(mdl_in_q, friction_velocity_flag)
      fric_vel_ptr => met_buf%p2d(iTmp)%present%ptr

      iTmp = fu_index(mdl_in_q, convective_velocity_scale_flag)
      conv_vel_scale_ptr => met_buf%p2d(iTmp)%present%ptr

    CASE DEFAULT
      CALL set_error('Unknown random-walk method','prepare_random_walk')
      RETURN
    END SELECT

  END SUBROUTINE prepare_random_walk


  ! ***************************************************************^L

  SUBROUTINE add_random_walk_input_needs(rw_method, q_met_dyn, q_met_stat, q_disp_dyn, q_disp_stat)
    !
    ! Returns the list of necessary parameters for calculation of the 
    ! random_walk of particles
    ! NOTE. It returns only upper_level parameters (e.g. abl_height), but
    !       does not care how they will be computed from the "raw" data
    !       That task has to be performed in derived_field_quantities module
    !
    ! Language: FORTRAN-90
    !
    ! Code owner: Mikhail Sofiev, FMI
    !
    IMPLICIT NONE

    ! Input stuff:
    INTEGER, INTENT(in) :: rw_method
!    LOGICAL, INTENT(in) :: post_proc  ! If post-processed

    ! Output - list of quantities
    INTEGER, DIMENSION(:), INTENT(inout) :: q_met_dyn, q_met_stat, q_disp_dyn, q_disp_stat

    ! Local variables
    integer :: iTmp

    SELECT CASE(rw_method)
    CASE (void_rw, fixed_rw) ! Nothing is needed - it is fixed_rw !
    
    CASE(fully_mixed_abl_rw)
      iTmp = fu_merge_integer_to_array(abl_top_pressure_flag, q_met_dyn)
      iTmp = fu_merge_integer_to_array(surface_pressure_flag, q_met_dyn)
      iTmp = fu_merge_integer_to_array(u_flag, q_met_dyn)
      iTmp = fu_merge_integer_to_array(v_flag, q_met_dyn)

    CASE (bulk_gaussian_rw)
      iTmp = fu_merge_integer_to_array(surface_pressure_flag, q_met_dyn)
      iTmp = fu_merge_integer_to_array(MO_length_inv_flag, q_met_dyn)
      iTmp = fu_merge_integer_to_array(friction_velocity_flag, q_met_dyn)
      iTmp = fu_merge_integer_to_array(convective_velocity_scale_flag, q_met_dyn)
      iTmp = fu_merge_integer_to_array(abl_top_pressure_flag, q_met_dyn)
      iTmp = fu_merge_integer_to_array(ABL_height_m_flag, q_met_dyn)
      iTmp = fu_merge_integer_to_array(u_flag, q_met_dyn)
      iTmp = fu_merge_integer_to_array(v_flag, q_met_dyn)
      iTmp = fu_merge_integer_to_array(temperature_flag, q_met_dyn)
      iTmp = fu_merge_integer_to_array(pressure_flag, q_met_dyn)
      iTmp = fu_merge_integer_to_array(height_flag, q_met_dyn)

    CASE DEFAULT
      CALL set_error('Unknown random walk method','add_random_walk_input_needs')
    END SELECT

  END SUBROUTINE add_random_walk_input_needs


  ! ***************************************************************^L
  ! ****************************************************************
  ! 
  !
  !       Private functions and subroutines of this module.
  !
  !
  ! ***************************************************************
  ! ****************************************************************

  subroutine random_walk_gaussian_prof(pos_x, pos_y, pos_z, &
                                     & timestep_sec, xCellSize_1, yCellSize_1, &
                                     & diff_x, diff_y, diff_pa, &
                                     & windspeed,&
                                     & surface_pressure,&
                                     & abl_height, &
                                     & MO_len, u_star, w_star, &
                                     & tempr, dt_dz,&
                                     & pressure, height, dp_dz,&
                                     & rw_method)
    !
    ! A random walk approach is used to parameterize horizontal and
    ! vertical diffusion. Slightly different parametrization is used
    ! for particles located within the ABL and those above the ABL.
    ! Vertical diffusion in the ABL is parametrized as a random
    ! relocation of the particles between surface and the ABL height,
    ! applied after each timestep. A simple method stops the
    ! particles from escaping to space.
    !
    ! Method: From Norwegian Meteorological Institute's SNAP-model
    ! 
    ! All units: SI
    !
    ! Language: ANSI Fortran 90
    !
    ! Author: Mikhail Sofiev, FMI
    ! Original code: Ilkka Valkama and Pilvi Siljamo, FMI
    ! 
  IMPLICIT NONE
    !
    ! Imported parameters
    real, intent(inout) :: pos_x, pos_y, pos_z
    REAL,INTENT(inout) :: diff_x,diff_y,diff_pa !Small-scale turbulent velocity
    real, INTENT(in) :: timestep_sec, xCellSize_1, yCellSize_1
    real, INTENT(in) :: windspeed, surface_pressure, abl_height, &
                      & MO_len, u_star, w_star, tempr, dt_dz, pressure, height, dp_dz
    INTEGER, INTENT(in) :: rw_method
                                                ! of the particle
    ! Local declarations:
    REAL :: t, old_position, new_position, reflection

    REAL, PARAMETER :: a1 = 0.25
    REAL, PARAMETER :: a2 = 0.5 
    REAL, PARAMETER :: b  = 0.875
    REAL, PARAMETER :: delta_pressure = 500 !Pa
    REAL, PARAMETER :: vertical_diffusion_parameter = 2. ! Pa/s

    INTEGER :: levelk
    INTEGER , PARAMETER :: maxlevels=31
    REAL, DIMENSION (maxlevels) :: pressure_level
    REAL, PARAMETER :: critical_richardson = 0.25
    REAL :: sigma_u, sigma_v, sigma_w = 0.1, dsigma_dz, air_dens, density_corr
    REAL ::  lagrangian_tu, lagrangian_tv, lagrangian_tw
        
    
!PRINT *, 'abl_height_pa', abl_height_pa

    !--------------------------------------------
    !
    ! Bottom-top-ABL Thresholds checking:
    !
    !--------------------------------------------
    IF (pos_z <= 10000.) THEN

      ! Above the top threshold. Don't let it climb above 100HPa.
      !
      pos_z = pos_z + 0.2*timestep_sec ! Force it down a bit

    ELSE IF (pos_z >= surface_pressure) THEN
      
      ! Below bottom threshold. Bring it to the surface
      !
      pos_z = surface_pressure - 1. ! [Pa]

    ELSE 
      
      ! Something non-trivial - either free atmosphere, or boundary layer
      !
      IF (pos_z <= abl_height) THEN ! Free atmosphere

        sigma_u = 0.5           ! Velocity standard deviation, unit = m/s
        sigma_v = 0.5           ! unit = m/s
        sigma_w = 0.1           ! unit = m/s

        lagrangian_tu = 300.    ! Lagrangian time scale, unit = s 
        lagrangian_tv = 300.    ! unit = s
        lagrangian_tw = 100.    ! unit = s

        pos_x = pos_x + fu_random_number_center(0.,0.5)*a1*((windspeed*timestep_sec)**b) * xCellSize_1
        pos_y = pos_y + fu_random_number_center(0.,0.5)*a1*((windspeed*timestep_sec)**b) * yCellSize_1
        pos_z = pos_z + fu_random_number_center(0., (timestep_sec*vertical_diffusion_parameter))

      ELSE ! Particle inside ABL :

        SELECT CASE (rw_method)
        CASE (bulk_gaussian_rw)

          CALL bulk_gaussian_rw_turb(height, abl_height, MO_len, u_star, w_star, &
                                   & sigma_u, sigma_v, sigma_w, &
                                   & lagrangian_tu, lagrangian_tv, lagrangian_tw, dsigma_dz)
          IF (error) RETURN

        CASE (advanced_gaussian)
          
!!!$ advanced_gaussian(position,inversion,MO_length,u_star,w_star&
!!!$        & ,sigma_u, sigma_v, sigma_w, lagrangian_tu, lagrangian_tv,&
!!!$        & lagrangian_tw, dsigma_dz)
          
          CALL set_error('Advanced Gaussian RW not ready', 'fu_random_walk_gaussian_prof')
          RETURN
        END SELECT

    
        !-----------------------------------
        !
        ! Compute random displacements due to 3D turbulence field 
        ! NOTE. Time autocorrelation of the turbulent velocity 
        ! fluctuations much shorter than model time step - it is a few
        ! seconds, while usual time step is tens of minutes. So, drift 
        ! component can often be simplified or totally ignored.

        !  horisontal : u-component
        !
        t=timestep_sec/lagrangian_tu !Measure of smallness of correction factors
        
        IF(t > 5.)THEN ! Long time step => skip drift and density correction
          pos_x = pos_x + fu_random_number_center(0.,0.5)*sigma_u * xCellSize_1
        
        ELSEIF (t <= 5. .and. t >=0.5)THEN !Detailed...
          t = EXP(-t)
          pos_x = pos_x + t*diff_x+fu_random_number_center(0.,0.5)*sigma_u*SQRT(1.-t*t) * xCellSize_1

        ELSE  ! Short time step. exponent->Tailor series
          pos_x = pos_x + &
                & (1.-t)*diff_x+fu_random_number_center(0.,0.5)*sigma_u*SQRT(2.*t) * xCellSize_1
        END IF

        ! Horizontal: v-component
        !
        t = timestep_sec/lagrangian_tv
        IF(t > 5.)THEN ! Long time step => skip drift and density correction
          pos_y = pos_y + fu_random_number_center(0.,0.5)*sigma_v * yCellSize_1

        ELSEIF (t <= 5. .and. t >=0.5)THEN !Detailed...
          t = EXP(-t)
          pos_y = pos_y + t*diff_y+fu_random_number_center(0.,0.5)*sigma_v*SQRT(1.-t*t) * yCellSize_1

        ELSE  ! t<0.5. A bit simplified: exponent -> Tailor series
          pos_y = pos_y + &
                & (1.-t)*diff_y+fu_random_number_center(0.,0.5)*sigma_v*SQRT(2.*t) * yCellSize_1
        END IF

        ! Vertical (pressure co-ordinate). 
        !      Density correction = d(density)/dz /density
        !
        t = timestep_sec/lagrangian_tw
        air_dens = pressure/(tempr*gas_constant_dryair)
        density_corr = (pressure*dt_dz-tempr*dp_dz) /  & 
                           & (pressure*pressure*air_dens)
        IF(t > 5.)THEN ! Long time step => skip drift
          pos_z = pos_z + g*air_dens* fu_random_number_center(0.,0.5)*sigma_w + &
                        & lagrangian_tw*(dsigma_dz + density_corr*sigma_w*sigma_w)
        ELSEIF (t <= 5. .and. t >=0.5)THEN !Detailed...

          t = EXP(-t)
          pos_z = pos_z + t*diff_pa + g*air_dens*fu_random_number_center(0.,0.5)* &
                        & sigma_w*SQRT(1. - t*t) + lagrangian_tw*(1. - t)* &
                        & (dsigma_dz + density_corr*sigma_w*sigma_w)
        ELSE  ! t<0.5 A bit simplified: exponent -> Tailor series
          pos_z = pos_z + (1.-t)*diff_pa + g*air_dens* &
                        & fu_random_number_center(0.,0.5)* &
                        & sigma_w*SQRT(2.*t) + lagrangian_tw*t* &
                        & (dsigma_dz + density_corr*sigma_w*sigma_w)
        END IF

      END IF ! Free atmosphere-ABL check
    END IF ! threshold_check

    IF (error) RETURN

  END subroutine random_walk_gaussian_prof


  ! ***************************************************************

  SUBROUTINE reflection_rules(old_position, inversion, new_position, reflection)
   
    ! Description:
    !     RETURNS the new vertical position for a particle.
    !
    !  method :
    !     Total reflection at the bottom and the top of atmosphere
    !     as well as mixing heigth (inversion) 
    
    ! Language: ANSI Fortran 90
    !
    ! Author: Mika Salonoja, FMI
    
    IMPLICIT NONE
    
    ! Local declarations:
    REAL, INTENT(in) :: old_position, inversion
    REAL, INTENT(out) :: new_position, reflection
    ! TYPE(silja_position) :: position
    ! Local declarations
    REAL :: displacement
    
    IF(ABS(old_position) > inversion)THEN
      displacement = old_position
    ELSE
      displacement = MOD(old_position, inversion)
    END IF
    
    IF(displacement > (-old_position))THEN
      ! bottom : reflected upwards 
      reflection = -1.0
      new_position = -old_position -  displacement
    ELSE
      IF(displacement > (inversion - old_position))THEN
        ! top : reflected downwards 
        reflection = -1.0
        new_position = -old_position -  displacement +2.0*inversion
      ELSE
        reflection = 1.0
        new_position = old_position + displacement
      END IF
    END IF
    
  END SUBROUTINE reflection_rules



! ***************************************************************
  
  SUBROUTINE bulk_gaussian_rw_turb(elevation, abl_h_m, MO_length, u_star, w_star, &
                                 & sigma_u, sigma_v, sigma_w, &
                                 & lagrangian_tu, lagrangian_tv, lagrangian_tw, dsigma_dz)
    !
    ! Evaluates turbulence statistics for random walk parametrization :
    !
    !     elevation [m]     altitude of a particel (above ground)   
    !     abl_h [m]         mixing layer height (ABL height)            
    !     MO_lentgh [m]     Monin-Obukhov length scale
    !     u_star [m/s]      friction velocity                                        
    !     w_star [m/s]      convective velocity scale                                               
    !     sigma_u, sigma_v  standard deviations for lateral turbulent velocity fluctuations     
    !     sigma_w           standard deviation for vertical turbulent velocity fluctuations    
    !     dsigma_dz [1/s]   vertical gradient of sigma_w 
    !     lagrangian_tu [s] Lagrangian time scale for the along wind component       
    !     lagrangian_tv [s] Lagrangian time scale for the cross wind component      
    !     lagrangian_tw [s] Lagrangian time scale for the vertical wind component
    !
    ! Method: Bulk schemes from Ryall & Maryon (1997) based on Hanna (1982)
    ! 
    ! All units: SI
    !
    ! Language: ANSI Fortran 90
    !
    ! Author: Ilkka Valkama, FMI
                                                                      
    IMPLICIT NONE

    ! Imported parameters with intent(in/out):
    real, INTENT(in) :: elevation
    REAL, INTENT(in) ::  abl_h_m, MO_length, u_star, w_star
    REAL, INTENT(out) :: sigma_u, sigma_v, sigma_w, dsigma_dz
    REAL, INTENT(out) :: lagrangian_tu,lagrangian_tv, lagrangian_tw

    ! Local declarations:
    REAL :: ref_z, corr_z, corr, coriolis, ustar

    !-----------------------------------
    ! 1. Setting internal parameters
    !-----------------------------------
!PRINT*,'  '
!PRINT*,'Enter bulk_gaussian_rw  '
    !  basic coriolis-force at 60N
    coriolis = 1.0e-4
    ref_z = elevation/abl_h_m
    !  corr_z = coriolis*1.0 - ref_z
    !    corr = -1.3
!    PRINT*,' elevate&refer&MOL :  ', elevation,ref_z,MO_length

    !-----------------------------------
    ! 2. Computing dispersion parameters
    !-----------------------------------

    IF (MO_length <= 0.) THEN

      ! 2.1 Unstably stratified ABL
      !
      sigma_u = u_star*(12.0 - 0.5*abl_h_m/MO_length)**0.33333
      sigma_v = sigma_u
      sigma_w = 1.2*w_star*w_star*(1.0 - 0.9*ref_z)*ref_z**0.66667
      sigma_w =  SQRT(sigma_w + (1.8 - 1.4*ref_z)*u_star**2)
      
      corr_z = -1.4*u_star*u_star + w_star*w_star
      corr_z =  corr_z + (0.8*MAX(ref_z,1.0e-3)**(-0.33333) - 1.8*ref_z**0.66666)
      dsigma_dz = (0.5/sigma_w)/abl_h_m*corr_z
      
      !  Lagrangian time scales   
      
      lagrangian_tu= 0.15*abl_h_m/sigma_u
      lagrangian_tv = lagrangian_tu
      
      IF (elevation < ABS(MO_length)) THEN
        corr_z = sigma_w*(0.55 - 0.38*ABS(elevation/MO_length))
        lagrangian_tw = 0.1*elevation/corr_z
        lagrangian_tw = 0.1*elevation/(sigma_w*(0.55-0.38*ABS(elevation/MO_length)))
      ELSE IF (ref_z < 0.1) THEN
        lagrangian_tw = 0.59*elevation/sigma_w
      ELSE                
        lagrangian_tw = 0.15*abl_h_m/sigma_w*(1.-EXP(-5*ref_z))
      ENDIF
 
    ELSEIF (abl_h_m/MO_length < 1.) THEN

        ! 2.2 Neutrally stratified ABL
        !
        ustar = MAX(1.0e-4,u_star)
        corr_z = elevation/ustar
        corr = coriolis*corr_z
       
        sigma_u = 2.0*ustar*EXP(-3.0*corr)
        sigma_w = 1.3*ustar*EXP(-2.0*corr)
        dsigma_dz = -2.0**coriolis*sigma_w
        sigma_v = sigma_w
        
        lagrangian_tu = 0.5*elevation/sigma_w/(1.0 + 15.0*corr)
        lagrangian_tv = lagrangian_tu
        lagrangian_tw = lagrangian_tu
        
    ELSE

        ! 2.2 Stably stratified ABL
        !
        sigma_u = 2.0*u_star*(1.0 - ref_z)
        sigma_v = 1.3*u_star*(1.0 - ref_z)
        sigma_w = sigma_v
        dsigma_dz= -1.3*u_star/abl_h_m
        
        lagrangian_tu = 0.15*abl_h_m/sigma_u*(SQRT(ref_z))
        lagrangian_tv = 0.467*lagrangian_tu
        lagrangian_tw = 0.1*abl_h_m/sigma_w*(ref_z)**0.8
        
    ENDIF
      
    ! Minimum values for Lagrangian time scales [sec]
    lagrangian_tu = MAX(10.0,lagrangian_tu)
    lagrangian_tv = MAX(10.0,lagrangian_tv)
    lagrangian_tw = MAX(30.0,lagrangian_tw)
      
    IF (dsigma_dz.EQ.0.0) dsigma_dz = 1.0e-10

  END SUBROUTINE bulk_gaussian_rw_turb

END MODULE advection_lagrangian
