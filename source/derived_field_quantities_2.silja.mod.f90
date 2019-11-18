MODULE derived_field_quantities_2

  ! Description:
  ! This module contains tools for calculating some new fields
  ! from the nwp weather data in supermarket. The new fields are
  ! stored to the supermarket along with the original, nwp's own fields.
  ! 
  ! Current code owner: Mikhail Sofiev, FMI
  ! Original code: Pilvi Siljamo, FMI
  !
  ! All units: SI
  ! 
  ! All degrees: real numbers in degrees (degrees and parts per
  ! hundred/parts per thousand , NOT in degrees and minutes
  ! monin_obukhov_flag
  ! Language: ANSI standard Fortran 90
  ! 
  ! Modules used:

  USE physiographies 

  !$use omp_lib
!  use supermarket_of_fields


  IMPLICIT NONE

  ! The public functions and subroutines available in this module:
  PUBLIC dq_eq_potential_temperature
  PUBLIC dq_potential_temperature
  public dq_temperature_from_perturb
  PUBLIC dq_relative_vorticity
  PUBLIC dq_abs_vorticity_advection
  PUBLIC dq_isentropic_pot_vorticity
  PUBLIC dq_thermal_front_parameter
  PUBLIC dq_bulk_thermal_front_parameter
  PUBLIC dq_bulk_richardson_number
  public dq_gradient_richardson_number
  public dq_brunt_vaisala_freq
  PUBLIC dq_ABL_params
  PUBLIC dq_omega
  PUBLIC dq_vertical_velocity_msl
  public dq_vertical_velocity_srf
  public dq_wind_divergence
  public dq_abl_top_pressure
!  PUBLIC dq_albedo
  PUBLIC dq_rad_sw_down_sfc
  public dq_cwcabove_3D
  PUBLIC dq_scavenging_coef
  public dq_scavenging_coef_ls
  PUBLIC dq_pasquill
  public dq_turbulence_profile
  public dq_aerodynamic_resistance
  public dq_hybrid_vertical_wind
  public dq_stomatal_conductance  !!
  public dq_lai_dyn 
  public dq_ref_evapotranspiration
  public fu_g_stomatal
!  public dq_sw_down_sfc

  ! The private functions and subroutines not to be used elsewhere:
  private bulk_richardson_nbr_one_level

  !
  ! The parameters for the similarity theory
  !
  real, public, parameter :: Kz_ref_height = 1.  ! 1 metre

  INTEGER, PRIVATE :: iRec = 1


CONTAINS

  ! ****************************************************************


  SUBROUTINE dq_eq_potential_temperature(meteoMarketPtr, met_src, valid_times)

    ! Description:
    ! Calculates equivalent potential temperature
    !
    ! All units: SI
    !
    ! Language: ANSI Fortran 90
    !
    ! Current code owner: Pilvi Siljamo, FMI
    ! 
    IMPLICIT NONE

    ! Imported parameters with intent IN:
    type(meteo_data_source), INTENT(in) :: met_src
    TYPE(silja_time), DIMENSION(:), INTENT(in) :: valid_times

    type(mini_market_of_stacks), pointer :: meteoMarketPtr
    ! Local declarations:
    TYPE(silja_time) :: time
    INTEGER :: i, tloop
    TYPE(silja_3d_field), POINTER :: t3d, q3d
    TYPE(silja_field_id) :: id
    REAL, DIMENSION(:), POINTER :: t
    REAL, DIMENSION(:), POINTER :: q ! specific humidity, shold be
    !                                  humidity mixing ratio
    TYPE(silja_level), DIMENSION(max_levels) :: levels
    INTEGER :: number_of_levels, fs
    REAL, DIMENSION(:), POINTER :: p, theta

    p => fu_work_array()
!    theta => fu_work_array()
    IF (error) RETURN

    loop_over_times: DO tloop = 1, SIZE(valid_times)
      IF (.NOT.defined(valid_times(tloop))) EXIT loop_over_times
      time = valid_times(tloop)

      IF (fu_field_in_sm(meteoMarketPtr, met_src, & ! Already done before?
                       & eq_pot_temperature_flag,&
                       & time,&
                       & level_missing,&
                       & .true.,&
                       & .false.)) CYCLE loop_over_times

      t3d => fu_sm_obstime_3d_field(meteoMarketPtr, met_src, temperature_flag, time, single_time)
      IF (error) CYCLE loop_over_times

      q3d => fu_sm_obstime_3d_field(meteoMarketPtr, met_src, specific_humidity_flag, time, single_time)
      IF (error) CYCLE loop_over_times

      CALL vertical_levels(t3d, levels, number_of_levels)
      IF (error) CYCLE loop_over_times

      loop_over_levels: DO i = 1, number_of_levels

        t => fu_grid_data_from_3d(t3d, i)

        q => fu_grid_data_from_3d(q3d, i)

        fs = SIZE(t)

        CALL pressure_on_level(t3d,i,p)

        id = fu_set_field_id(fu_met_src(t3d), eq_pot_temperature_flag,&
                           & fu_analysis_time(t3d),&
                           & fu_forecast_length(t3d), fu_grid(t3d), levels(i))
        call find_field_data_storage_2d(meteoMarketPtr, id, multi_time_stack_flag, theta)
        if(fu_fails(.not.error,'Failed theta field data pointer','dq_eq_potential_temperature'))return
        !
        ! From Stull 1988
        !
        theta(1:fs) = (t(1:fs) + &
                     & ((condensation_latentheat/specific_heat_dryair) * q(1:fs)))&
                     & *(1.E5/p(1:fs))**R_per_c_dryair


        !        IF ((SUM(theta(1:fs))/REAL(fs)) > 350.) theta = real_missing

!        CALL dq_store_2d(meteoMarketPtr, id, theta, multi_time_stack_flag )

      END DO loop_over_levels
    END DO loop_over_times


    CALL free_work_array(p)
!    CALL free_work_array(theta)

  END SUBROUTINE dq_eq_potential_temperature



  ! ****************************************************************


  SUBROUTINE dq_potential_temperature(meteoMarketPtr, met_src, valid_times)
    !
    ! Calculates potential temperature for both 3D and 2m temperature fields
    !
    IMPLICIT NONE

    ! Imported parameters with intent IN:
    type(meteo_data_source), INTENT(in) :: met_src
    TYPE(silja_time), DIMENSION(:), INTENT(in) :: valid_times

    type(mini_market_of_stacks), pointer :: meteoMarketPtr

    ! Local declarations:
    TYPE(silja_time) :: time
    INTEGER :: lev, i, tloop, number_of_levels
    TYPE(silja_3d_field), POINTER :: t3d
    TYPE(silja_field), POINTER :: fldTmp
    TYPE(silja_field_id) :: id
    TYPE(silja_level), DIMENSION(max_levels) :: levels
    REAL, DIMENSION(:), POINTER :: p, psrf, theta, t

    p => fu_work_array()
!    theta => fu_work_array()
    IF (error) RETURN

    DO tloop = 1, SIZE(valid_times)
      IF (.NOT.defined(valid_times(tloop))) EXIT ! loop_over_times
      time = valid_times(tloop)

      IF (fu_field_in_sm(meteoMarketPtr, met_src, & ! Already done before?
                       & potential_temperature_flag,&
                       & time,&
                       & level_missing,&
                       & .true.,&
                       & .false.)) CYCLE !loop_over_times
      !
      ! Prepare the fields
      !
      t3d => fu_sm_obstime_3d_field(meteoMarketPtr, met_src,&
                                  & temperature_flag, &
                                  & time,&
                                  & single_time)
      IF (error) CYCLE !loop_over_times

      CALL vertical_levels(t3d, levels, number_of_levels)
      IF (error) CYCLE !loop_over_times

      fldTmp => fu_sm_obstime_field(meteoMarketPtr, met_src, ground_pressure_flag, level_missing, time, forwards)
      if(error)return
      psrf => fu_grid_data(fldTmp)

      fldTmp => fu_sm_obstime_field(meteoMarketPtr, met_src, temperature_2m_flag, level_missing, time, forwards)
      if(error)return
      t => fu_grid_data(fldTmp)

      id = fu_set_field_id(fu_met_src(fldTmp),&
                         & potential_temperature_2m_flag,&
                         & fu_analysis_time(fldTmp),&
                         & fu_forecast_length(fldTmp), &
                         & fu_grid(fldTmp),&
                         & fu_level(fldTmp))
      call find_field_data_storage_2d(meteoMarketPtr, id, multi_time_stack_flag, theta)
      if(fu_fails(.not.error,'Failed theta field data pointer','dq_potential_temperature'))return
      
!!      !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i)
      do i=1,fs_meteo
        theta(i) = t(i)*((std_pressure_sl/psrf(i))**R_per_c_dryair)  
      end do
!!      !$OMP END PARALLEL DO

!      CALL dq_store_2d(meteoMarketPtr, id, theta, multi_time_stack_flag )
!      if(error)return

      !
      ! Now make the whole 3D potential temperature field
      ! 
      DO lev = 1, number_of_levels
        
        t => fu_grid_data_from_3d(t3d, lev)

        CALL pressure_on_level(t3d, lev, p)  ! It really computes the values
        
        id = fu_set_field_id(fu_met_src(t3d),&
                           & potential_temperature_flag,&
                           & fu_analysis_time(t3d),&
                           & fu_forecast_length(t3d), &
                           & fu_grid(t3d),&
                           & levels(lev))
        call find_field_data_storage_2d(meteoMarketPtr, id, multi_time_stack_flag, theta)
        if(fu_fails(.not.error,'Failed theta field data pointer','dq_potential_temperature'))return
        
!!        !$OMP PARALLEL DO DEFAULT(SHARED)
        do i=1,fs_meteo
!          if(i==12379)then
!            call msg('Pot temp: get to i',i)
!          endif
          theta(i) = t(i)*((std_pressure_sl/p(i))**R_per_c_dryair)  
        end do
!!        !$OMP END PARALLEL DO

!        CALL dq_store_2d(meteoMarketPtr, id, theta, multi_time_stack_flag )

      END DO  !, loop over levels
    END DO ! loop_over_times

    CALL free_work_array(p)
!    CALL free_work_array(theta)

  END SUBROUTINE dq_potential_temperature



  ! ****************************************************************


  SUBROUTINE dq_temperature_from_perturb(meteoMarketPtr, met_src, valid_times)
    !
    ! Calculates temperature for potential perturbational temperature
    !
    IMPLICIT NONE

    ! Imported parameters with intent IN:
    type(meteo_data_source), INTENT(in) :: met_src
    TYPE(silja_time), DIMENSION(:), INTENT(in) :: valid_times

    type(mini_market_of_stacks), pointer :: meteoMarketPtr

    ! Local declarations:
    TYPE(silja_time) :: time
    INTEGER :: lev, i, tloop, number_of_levels
    TYPE(silja_3d_field), POINTER :: pert_t3d
    TYPE(silja_field_id) :: id
    TYPE(silja_level), DIMENSION(max_levels) :: levels
    REAL, DIMENSION(:), POINTER :: p, psrf, pert_theta, t

    p => fu_work_array()
!    t => fu_work_array()
    IF (error) RETURN

    DO tloop = 1, SIZE(valid_times)
      IF (.NOT.defined(valid_times(tloop))) EXIT ! loop_over_times
      time = valid_times(tloop)

      IF (fu_field_in_sm(meteoMarketPtr, met_src, & ! Already done before?
                       & temperature_flag,&
                       & time,&
                       & level_missing,&
                       & .true.,&
                       & .false.)) CYCLE !loop_over_times
      !
      ! Prepare the fields
      !
      pert_t3d => fu_sm_obstime_3d_field(meteoMarketPtr, met_src,&
                                       & perturb_pot_temperature_flag, &
                                       & time,&
                                       & single_time)
      IF (error) CYCLE !loop_over_times

      CALL vertical_levels(pert_t3d, levels, number_of_levels)
      IF (error) CYCLE !loop_over_times

      !
      ! Now make the whole 3D potential temperature field
      ! 
      DO lev = 1, number_of_levels
        
        pert_theta => fu_grid_data_from_3d(pert_t3d, lev)

        CALL pressure_on_level(pert_t3d, lev, p)  ! It really computes the values
        if(error)return

        id = fu_set_field_id(fu_met_src(pert_t3d),&
                           & temperature_flag,&
                           & fu_analysis_time(pert_t3d),&
                           & fu_forecast_length(pert_t3d), &
                           & fu_grid(pert_t3d),&
                           & levels(lev))
        call find_field_data_storage_2d(meteoMarketPtr, id, multi_time_stack_flag, t)
        if(fu_fails(.not.error,'Failed t field data pointer','dq_temperature_from_perturb'))return
        
        do i=1,fs_meteo
          t(i) = (pert_theta(i)+300.)/((std_pressure_sl/p(i))**R_per_c_dryair)  
        end do

!        CALL dq_store_2d(meteoMarketPtr, id, t, multi_time_stack_flag )

      END DO  !, loop over levels
    END DO ! loop_over_times

    CALL free_work_array(p)
!    CALL free_work_array(t)

  END SUBROUTINE dq_temperature_from_perturb


  !**************************************************************

  SUBROUTINE dq_relative_vorticity(meteoMarketPtr, met_src, valid_times)

    ! Description:
    ! Calculates relative vorticity
    !
    ! All units: SI
    !
    ! Language: ANSI Fortran 90
    !
    ! Current code owner: Pilvi Siljamo, FMI
    ! 
    IMPLICIT NONE

    ! Imported parameters with intent IN:
    type(meteo_data_source), INTENT(in) :: met_src
    TYPE(silja_time), DIMENSION(:), INTENT(in) :: valid_times
    type(mini_market_of_stacks), pointer :: meteoMarketPtr

    ! Local declarations:
    TYPE(silja_time) :: time
    TYPE(silja_3d_field), POINTER :: u3d, v3d, t3d
    TYPE(silja_field_id) :: id
    TYPE(silja_grid) :: u_grid, v_grid, scalar_grid
    TYPE(silja_level), DIMENSION(max_levels) :: levels
    INTEGER :: number_of_levels
    REAL, DIMENSION(:), POINTER :: u, v, dv_dx, du_dy, vort
!    REAL, DIMENSION(:), POINTER :: u_up, u_down, v_up, v_down
    REAL, DIMENSION(:), POINTER :: pressure
!    REAL, DIMENSION(:), POINTER :: p_up, p_down, dp_dx, dp_dy
    INTEGER :: l, fs, tloop

    dv_dx => fu_work_array()
    du_dy => fu_work_array()
!    dp_dx => fu_work_array()
!    dp_dy => fu_work_array()
    pressure => fu_work_array()
!    p_up => fu_work_array()
!    p_down => fu_work_array()
!    vort => fu_work_array()
    IF (error) RETURN

    loop_over_times: DO tloop = 1, SIZE(valid_times)
      IF (.NOT.defined(valid_times(tloop))) EXIT loop_over_times
      time = valid_times(tloop)

      IF (fu_field_in_sm(meteoMarketPtr, met_src, & ! Already done before?
                       & relative_vorticity_flag,&
                       & time,&
                       & level_missing,&
                       & .true.,&
                       & .false.)) CYCLE loop_over_times

      t3d => fu_sm_obstime_3d_field(meteoMarketPtr, met_src,&
                                  & temperature_flag, &
                                  & time,&
                                  & single_time)
      IF (error) CYCLE loop_over_times

      u3d => fu_sm_obstime_3d_field(meteoMarketPtr, met_src,&
                                  & u_flag, &
                                  & time,&
                                  & single_time)
      IF (error) CYCLE loop_over_times

      v3d => fu_sm_obstime_3d_field(meteoMarketPtr, met_src,&
                                  & v_flag, &
                                  & time,&
                                  & single_time)
      IF (error) CYCLE loop_over_times

      CALL vertical_levels(t3d, levels, number_of_levels)
      IF (error) CYCLE loop_over_times

      scalar_grid = fu_grid(t3d)
      u_grid = fu_grid(u3d)
      v_grid = fu_grid(v3d)
      IF (error) CYCLE loop_over_times

!      vort = 0.
      fs = fu_number_of_gridpoints(scalar_grid)

      loop_over_levels: DO l = 1, number_of_levels

!!!$        IF ((l == 1).or.(l == number_of_levels)) THEN
!!!$          id = fu_set_field_id(&
!!!$              & fu_met_src(t3d),&
!!!$              & relative_vorticity_flag,&
!!!$              & fu_analysis_time(t3d),&
!!!$              & fu_forecast_length(t3d), &
!!!$              & fu_grid(t3d),&
!!!$              & levels(l))
!!!$          vort = real_missing
!!!$          CALL dq_store_2d(meteoMarketPtr, id, vort)
!!!$          IF (error) RETURN
!!!$          CYCLE loop_over_levels
!!!$        END IF

        u => fu_grid_data_from_3d(u3d, l)
        v => fu_grid_data_from_3d(v3d, l)
        IF (error) CYCLE loop_over_levels

!!!$        u_up => fu_grid_data_from_3d(u3d, l+1)
!!!$        v_up => fu_grid_data_from_3d(v3d, l+1)
!!!$        IF (error) RETURN
!!!$
!!!$        u_down => fu_grid_data_from_3d(u3d, l-1)
!!!$        v_down => fu_grid_data_from_3d(v3d, l-1)
!!!$        IF (error) RETURN
!!!$
!!!$        CALL pressure_on_level(t3d, l, pressure) 
!!!$        IF (error) RETURN
!!!$
!!!$        CALL pressure_on_level(t3d, l+1, p_up) 
!!!$        IF (error) RETURN
!!!$        
!!!$        CALL pressure_on_level(t3d, l-1, p_down) 
!!!$        IF (error) RETURN

        ! dv_dx:
        CALL ddx_of_field(v, v_grid, scalar_grid, dv_dx)
        IF (error) CYCLE loop_over_levels

        ! du_dy:
        CALL ddy_of_field(u, u_grid, scalar_grid, du_dy)
        IF (error) CYCLE loop_over_levels

!!!$        ! dp_dx:
!!!$        CALL ddx_of_field(pressure, scalar_grid, scalar_grid, dp_dx)
!!!$        IF (error) RETURN
!!!$        ! dp_dy:
!!!$        CALL ddy_of_field(pressure, scalar_grid, scalar_grid, dp_dy)
!!!$        IF (error) RETURN
!!!$
!!!$        vort(1:fs) = &
!!!$            & (dv_dx(1:fs)-du_dy(1:fs))+ &
!!!$            & (dp_dx(1:fs)*(v_up(1:fs)-v_down(1:fs))- &
!!!$            &  dp_dy(1:fs)*(u_up(1:fs)-u_down(1:fs)))/ &
!!!$            & (p_up(1:fs)-p_down(1:fs))

        id = fu_set_field_id(fu_met_src(t3d),&
                           & relative_vorticity_flag,&
                           & fu_analysis_time(t3d),&
                           & fu_forecast_length(t3d), &
                           & scalar_grid,&
                           & levels(l))
        call find_field_data_storage_2d(meteoMarketPtr, id, multi_time_stack_flag, vort)
        if(fu_fails(.not.error,'Failed vort field data pointer','dq_relative_vorticity'))return

        vort(1:fs) = dv_dx(1:fs)-du_dy(1:fs)

!        PRINT*,'rel_vort max',MAXVAL(vort(1:fs))
!        PRINT*,'rel_vort min',MINVAL( vort(1:fs))
        

!        CALL dq_store_2d(meteoMarketPtr, id, vort, multi_time_stack_flag )

      END DO loop_over_levels
    END DO loop_over_times

!    CALL free_all_work_arrays()
    call free_work_array(dv_dx)
    call free_work_array(du_dy)
    call free_work_array(pressure)
!    call free_work_array(vort)

  END SUBROUTINE dq_relative_vorticity




  !********************************************************************

  SUBROUTINE dq_abs_vorticity_advection(meteoMarketPtr, met_src, valid_times)

    ! Description:
    ! Calculates absolute vorticity advection [m/s2]
    !
    ! All units: SI
    !
    ! Language: ANSI Fortran 90
    !
    ! Current code owner: Pilvi Siljamo, FMI
    ! 
    IMPLICIT NONE

    ! Imported parameters with intent IN:
    type(meteo_data_source), INTENT(in) :: met_src
    TYPE(silja_time), DIMENSION(:), INTENT(in) :: valid_times
    type(mini_market_of_stacks), pointer :: meteoMarketPtr
    ! Local declarations:
    TYPE(silja_time) :: time
    TYPE(silja_3d_field), POINTER :: u3d, v3d, t3d, relvort3d
    TYPE(silja_field_id) :: id
    TYPE(silja_grid) :: u_grid, v_grid, scalar_grid
    TYPE(silja_level), DIMENSION(max_levels) :: levels
    INTEGER :: number_of_levels, nx, ny
    REAL, DIMENSION(:), POINTER :: u, v, coriolis, relvort
    REAL, DIMENSION(:), POINTER :: absolute_vorticity, p 
    REAL, DIMENSION(:), POINTER :: dabsvort_dx,dabsvort_dy,abs_vort_adv
    REAL, DIMENSION(:), POINTER :: dabsvort_dx_f,dabsvort_dy_f
    REAL, DIMENSION(:), POINTER :: dabsvort_dx_b,dabsvort_dy_b
    REAL, DIMENSION(:), POINTER :: ahx, ahxm
    REAL :: ahy
    INTEGER :: i, l, fs, tloop, n

    coriolis => fu_work_array()
    absolute_vorticity => fu_work_array()
    p => fu_work_array()
    dabsvort_dx => fu_work_array()
    dabsvort_dy => fu_work_array() 
!    abs_vort_adv => fu_work_array()

    !    dabsvort_dx_f => fu_work_array()
    !    dabsvort_dy_f => fu_work_array()
    !    dabsvort_dx_b => fu_work_array()
    !    dabsvort_dy_b => fu_work_array()
    !    ahx => fu_work_array()
    !    ahxm => fu_work_array()
    
    IF (error) RETURN

    loop_over_times: DO tloop = 1, SIZE(valid_times)
      IF (.NOT.defined(valid_times(tloop))) EXIT loop_over_times
      time = valid_times(tloop)

      IF (fu_field_in_sm(meteoMarketPtr, met_src, & ! Already done before?
                       & abs_vorticity_advection_flag,&
                       & time,&
                       & level_missing,&
                       & .true.,&
                       & .false.)) CYCLE loop_over_times

      first_loop: IF (tloop == 1) THEN

        IF (.NOT.(fu_field_in_sm(meteoMarketPtr, met_src, & ! Already done before?
                               & relative_vorticity_flag,&
                               & time,&
                               & level_missing,&
                               & .true.,&
                               & .false.))) THEN
          CALL set_error('abs.vort.adv. required but no rel.vort.',&
                       & 'dq_abs_vorticity_advection')
          RETURN
        END IF

        t3d => fu_sm_obstime_3d_field(meteoMarketPtr, met_src,&
                                    & temperature_flag, &
                                    & time,&
                                    & single_time)
        IF (error) RETURN

        scalar_grid = fu_grid(t3d)
        CALL grid_dimensions(scalar_grid, nx, ny)
        IF (error) RETURN

        !        CALL metric_coefficients(scalar_grid, ahx, ahxm, ahy)

        CALL vertical_levels(t3d, levels, number_of_levels)
        IF (error) RETURN

        CALL coriolis_parameters(scalar_grid, coriolis)
        IF (error) RETURN

        CALL grid_dimensions(scalar_grid, nx, ny)
        IF (error) RETURN

        fs = fu_number_of_gridpoints(scalar_grid)

      END IF first_loop


      u3d => fu_sm_obstime_3d_field(meteoMarketPtr, met_src,&
                                  & u_flag, &
                                  & time,&
                                  & single_time)
      IF (error) RETURN

      v3d => fu_sm_obstime_3d_field(meteoMarketPtr, met_src,&
                                  & v_flag, &
                                  & time,&
                                  & single_time)
      IF (error) RETURN

      relvort3d => fu_sm_obstime_3d_field(meteoMarketPtr, met_src,&
                                        & relative_vorticity_flag, &
                                        & time,&
                                        & single_time)
      IF (error) RETURN

      IF (.NOT.(fu_grid(relvort3d) == scalar_grid)) THEN
        CALL report(fu_grid(relvort3d))
        CALL report(scalar_grid)
        CALL set_error('rel.vort. not in scalar grid',&
            & 'dq_abs_vorticity_advection')
        RETURN
      END IF

      u_grid = fu_grid(u3d)
      v_grid = fu_grid(v3d)
      IF (error) RETURN

      abs_vort_adv = 0.

      loop_over_levels: DO l = 1, number_of_levels

        u => fu_grid_data_from_3d(u3d, l)
        v => fu_grid_data_from_3d(v3d, l)
        IF (error) RETURN

        relvort => fu_grid_data_from_3d(relvort3d, l)
        IF (error) RETURN

        absolute_vorticity(1:fs) = relvort(1:fs) + coriolis(1:fs)

          ! dabsvort_dx
          CALL ddx_of_field(absolute_vorticity,&
                          & scalar_grid,&
                          & scalar_grid,&
                          & dabsvort_dx)
          IF (error) RETURN

          ! dabsvort_dy
          CALL ddy_of_field(absolute_vorticity,&
                          & scalar_grid,&
                          & scalar_grid,&
                          & dabsvort_dy)
          IF (error) RETURN

          IF (.NOT.(fu_grid(u3d) == fu_grid(v3d))) THEN ! Arakawa C
            !  CALL report(fu_grid(u3d))
            !  CALL report(fu_grid(v3d))
            DO n=nx, fs 
              u(n) = (u(n-1) + u(n))/2
              v(n) = (v(n)+v(n-nx))/2
            END DO
          END IF

          id = fu_set_field_id(fu_met_src(t3d),&
                             & abs_vorticity_advection_flag,&
                             & fu_analysis_time(t3d),&
                             & fu_forecast_length(u3d), &
                             & fu_grid(t3d),&
                             & levels(l))
          call find_field_data_storage_2d(meteoMarketPtr, id, multi_time_stack_flag, abs_vort_adv)
          if(fu_fails(.not.error,'Failed abs_vort_adv field data pointer','dq_abs_vorticity_advection'))return

          abs_vort_adv(1:fs) = (-1)*(u(1:fs) * dabsvort_dx(1:fs) + &
                                   & v(1:fs) * dabsvort_dy(1:fs))

!!!$        ! Better, but too slow way calculating vorticity advection
!!!$
!!!$        ! dabsvort_dx_forward
!!!$        PRINT*,'ddx_of_field_forward'
!!!$        CALL ddx_of_field_forward(absolute_vorticity,&
!!!$            & scalar_grid,&
!!!$            & u_grid,&
!!!$            & dabsvort_dx_f)
!!!$        IF (error) RETURN
!!!$        PRINT*,'ddy_of_field_forward'
!!!$        ! dabsvort_dy_forward
!!!$        CALL ddy_of_field_forward(absolute_vorticity,&
!!!$            & scalar_grid,&
!!!$            & v_grid,&
!!!$            & dabsvort_dy_f)
!!!$        IF (error) RETURN
!!!$        PRINT*,'ddx_of_field_backward'
!!!$        ! dabsvort_dx_backward
!!!$        CALL ddx_of_field_backward(absolute_vorticity,&
!!!$            & scalar_grid,&
!!!$            & u_grid,&
!!!$            & dabsvort_dx_b)
!!!$        IF (error) RETURN
!!!$        PRINT*,'dabsvort_dy_backward'  
!!!$        ! dabsvort_dy_backward
!!!$        CALL ddy_of_field_backward(absolute_vorticity,&
!!!$            & scalar_grid,&
!!!$            & v_grid,&
!!!$            & dabsvort_dy_b)
!!!$        IF (error) RETURN
!!!$        PRINT*,' CALL pressure_on_level'
!!!$        CALL pressure_on_level(t3d, l, p) 
!!!$        IF (error) RETURN
!!!$        PRINT*,'DO-luuppi'
!!!$        abs_vort_adv(1:fs) = real_missing
!!!$        DO i= 1, fs
!!!$          abs_vort_adv = (0.25/p(i))*&
!!!$              &((u(i)  *(p(i)+p(i+1)) *dabsvort_dx_f(i)+&
!!!$              &  u(i-1)*(p(i)+p(i-1)) *dabsvort_dx_b(i))+&
!!!$              &(v(i)   *(p(i)+p(i+nx))*dabsvort_dy_f *(ahxm(i)/ahx(i))+&
!!!$              & v(i-nx)*(p(i)+p(i-nx))*dabsvort_dy_b)*(ahxm(i-nx)/ahx(i)))
!!!$        END DO
!!!$
!!!$          abs_vort_adv(1:fs) = (-1)*(u(1:fs) * dabsvort_dx(1:fs)) + &
!!!$              & (v(1:fs) * dabsvort_dy(1:fs))
!!!$
!!!$        PRINT*,'end do'
!      END IF
 

!        CALL dq_store_2d(meteoMarketPtr, id, abs_vort_adv, multi_time_stack_flag )
        IF (error) RETURN

      END DO loop_over_levels
    END DO loop_over_times

!    CALL free_all_work_arrays()
    call free_work_array(coriolis)
    call free_work_array(absolute_vorticity)
    call free_work_array(p)
    call free_work_array(dabsvort_dx)
    call free_work_array(dabsvort_dy) 
!    call free_work_array(abs_vort_adv)

  END SUBROUTINE dq_abs_vorticity_advection



  !****************************************************************

  SUBROUTINE dq_isentropic_pot_vorticity(meteoMarketPtr, met_src, valid_times)

    ! Description:
    ! Calculates isentropic potential vorticity (IPV) 
    ! [(K*m2)/(kg*s)]
    !
    ! All units: SI
    !
    ! Language: ANSI Fortran 90
    !
    ! Current code owner: Pilvi Siljamo, FMI
    ! 
    IMPLICIT NONE

    ! Imported parameters with intent IN:
    type(meteo_data_source), INTENT(in) :: met_src
    TYPE(silja_time), DIMENSION(:),INTENT(in) :: valid_times
    type(mini_market_of_stacks), pointer :: meteoMarketPtr
    ! Local declarations:
    TYPE(silja_time) :: time
    TYPE(silja_3d_field), POINTER :: u3d, v3d, theta3d, relvort3d  !, t3d
    TYPE(silja_field_id) :: id
    TYPE(silja_grid) :: u_grid, v_grid, scalar_grid
    TYPE(silja_level), DIMENSION(max_levels) :: levels
    INTEGER :: number_of_levels
    REAL, DIMENSION(:), POINTER :: u, v, coriolis, relvort
    INTEGER :: lev, fs, tloop, n, nx, ny
!    REAL, DIMENSION(:), POINTER :: absolute_vorticity 
!    REAL, DIMENSION(:), POINTER :: U_up, U_down, V_up, V_down
    REAL, DIMENSION(:), POINTER :: pressure_up, pressure_down
    REAL, DIMENSION(:), POINTER :: theta_up, theta_down   !, theta
!    REAL, DIMENSION(:), POINTER :: dtheta_dy, dtheta_dx
    REAL, DIMENSION(:), POINTER :: du_dy, dv_dx, IPV
    logical :: ifMetaData

    coriolis => fu_work_array()
    pressure_up => fu_work_array()
    pressure_down => fu_work_array()
!    dtheta_dx => fu_work_array()
!    dtheta_dy => fu_work_array()
    dv_dx => fu_work_array()
    du_dy => fu_work_array()
!    IPV => fu_work_array()
    ifMetaData = .true.
    IF (error) RETURN
    !
    ! Cycle the requested time making the IPV if it is absent
    !
    DO tloop = 1, SIZE(valid_times)
      IF (.NOT.defined(valid_times(tloop))) EXIT
      time = valid_times(tloop)

      IF (fu_field_in_sm(meteoMarketPtr, met_src, & ! Already done before?
                       & ipv_flag, &
                       & time, &
                       & level_missing, &
                       & .true., &
                       & .false.)) CYCLE

      theta3d => fu_sm_obstime_3d_field(meteoMarketPtr, met_src, potential_temperature_flag, time, single_time)
      IF (error) RETURN

      IF (ifMetadata) THEN
        scalar_grid = fu_grid(theta3d)
        IF (error) RETURN
        CALL vertical_levels(theta3d, levels, number_of_levels)
        CALL coriolis_parameters(scalar_grid, coriolis)
        CALL grid_dimensions(scalar_grid, nx, ny)
        fs = nx*ny
        IF (error) RETURN 
        ifMetadata = .false.
      END IF

      u3d => fu_sm_obstime_3d_field(meteoMarketPtr, met_src, u_flag, time, single_time)
      v3d => fu_sm_obstime_3d_field(meteoMarketPtr, met_src, v_flag, time, single_time)
      IF (error) RETURN

      u_grid = fu_grid(u3d)
      v_grid = fu_grid(v3d)
      IF (error) RETURN
      !
      ! Cycle over vertical levels. Note asymmetric derivatives at the top and bottom
      !
      DO lev = 1, number_of_levels

        u => fu_grid_data_from_3d(u3d, lev)
        v => fu_grid_data_from_3d(v3d, lev)
        IF (error) RETURN

        ! du_dy
        CALL ddy_of_field(u, u_grid, scalar_grid, du_dy)
        IF (error) RETURN

        ! dv_dx 
        CALL ddx_of_field(v, u_grid, scalar_grid, dv_dx)
        IF (error) RETURN

        CALL pressure_on_level(theta3d, min(lev+1,number_of_levels), pressure_up) 
        CALL pressure_on_level(theta3d, max(lev-1,1), pressure_down)
        IF (error) RETURN

        ! potential temperatures
        theta_up => fu_grid_data_from_3d(theta3d, min(lev+1,number_of_levels))
        theta_down => fu_grid_data_from_3d(theta3d, max(lev-1,1))
        IF (error) RETURN

        !
        ! Calculate isentropic potential vorticity
        !
        id = fu_set_field_id(fu_met_src(theta3d),&
                           & ipv_flag,&
                           & fu_analysis_time(theta3d),&
                           & fu_forecast_length(theta3d), &
                           & fu_grid(theta3d),&
                           & levels(lev))
        call find_field_data_storage_2d(meteoMarketPtr, id, multi_time_stack_flag, IPV)
        if(fu_fails(.not.error,'Failed IPV field data pointer','dq_isentropic_pot_vorticity'))return
        
        IPV(1:fs) = -g * (dv_dx(1:fs)-du_dy(1:fs) + coriolis(1:fs)) * &
                       & (theta_up(1:fs)-theta_down(1:fs)) / (pressure_up(1:fs)-pressure_down(1:fs))

!        CALL dq_store_2d(meteoMarketPtr, id, IPV, multi_time_stack_flag )
        IF (error) RETURN
      END DO
    END DO

!    CALL free_all_work_arrays()
    call free_work_array(coriolis)
    call free_work_array(pressure_up)
    call free_work_array(pressure_down)
!    call free_work_array(dtheta_dx)
!    call free_work_array(dtheta_dy)
    call free_work_array(dv_dx)
    call free_work_array(du_dy)
!    call free_work_array(IPV)

  END SUBROUTINE dq_isentropic_pot_vorticity




  !**********************************************************************


  SUBROUTINE dq_thermal_front_parameter(meteoMarketPtr, met_src, valid_times)

    ! Description:
    ! Calculates thermal front parameter (TFP) [K/m2]
    !
    ! TFP = vector product [(change of temperaturegradient) and
    !       (projetion in direction of the temperature gradient)]
    !
    ! All units: SI
    !
    ! Language: ANSI Fortran 90
    !
    ! Current code owner: Pilvi Siljamo, FMI 
    ! 
    IMPLICIT NONE

    ! Imported parameters with intent IN:
    type(meteo_data_source), INTENT(in) :: met_src
    TYPE(silja_time), DIMENSION(:), INTENT(in) :: valid_times
    type(mini_market_of_stacks), pointer :: meteoMarketPtr
    ! Local declarations:
    TYPE(silja_time) :: time
    TYPE(silja_3d_field), POINTER :: t3d
    TYPE(silja_field_id) :: id
    TYPE(silja_grid) :: scalar_grid
    TYPE(silja_level), DIMENSION(max_levels) :: levels
    INTEGER :: number_of_levels
    INTEGER :: l, fs, tloop, i
    REAL,DIMENSION(:), POINTER :: T
    REAL,DIMENSION(:), POINTER :: abs_nabla_T 
    REAL,DIMENSION(:), POINTER :: dT_dy, dT_dx
    REAL,DIMENSION(:), POINTER :: dabsnablaT_dx, dabsnablaT_dy
    REAL,DIMENSION(:), POINTER :: TFP

    abs_nabla_T => fu_work_array() 
!    TFP => fu_work_array()
    dT_dx => fu_work_array()
    dT_dy => fu_work_array()
    dabsnablaT_dx => fu_work_array()
    dabsnablaT_dy => fu_work_array()
    IF (error) RETURN

    loop_over_times: DO tloop = 1, SIZE(valid_times)
      IF (.NOT.defined(valid_times(tloop))) EXIT loop_over_times
      time = valid_times(tloop)

      IF (fu_field_in_sm(meteoMarketPtr, met_src, & ! Already done before?
                       & tfp_flag,&
                       & time,&
                       & level_missing,&
                       & .true.,&
                       & .false.)) CYCLE loop_over_times

      t3d => fu_sm_obstime_3d_field(meteoMarketPtr, met_src, temperature_flag, time, single_time)
      IF (error) RETURN

      CALL vertical_levels(t3d, levels, number_of_levels)
      IF (error) RETURN

      scalar_grid = fu_grid(t3d)
      IF (error) RETURN

      TFP = 0.
      fs = fu_number_of_gridpoints(scalar_grid)

      loop_over_levels: DO l = 1, number_of_levels

        T => fu_grid_data_from_3d(t3d, l)
        IF (error) RETURN

         ! dT_dx:
        CALL ddx_of_field(T, scalar_grid, scalar_grid, dT_dx)
        IF (error) RETURN

        ! dT_dy:
        CALL ddy_of_field(T, scalar_grid, scalar_grid, dT_dy)
        IF (error) RETURN

        abs_nabla_T(1:fs) = SQRT(dT_dx(1:fs)**2. + dT_dy(1:fs)**2.)

        ! dabsnablaT_dx 
        CALL ddx_of_field(abs_nabla_T(1:fs), scalar_grid, scalar_grid,dabsnablaT_dx(1:fs))
        IF (error) RETURN

        ! dabsnablaT_dy 
        CALL ddy_of_field(abs_nabla_T(1:fs), scalar_grid, scalar_grid,dabsnablaT_dy(1:fs))
        IF (error) RETURN


        ! ------------------------------------------
        !
        ! Calculates thermal front parameter
        !
        ! ------------------------------------------

        id = fu_set_field_id(fu_met_src(t3d), &
                           & tfp_flag, &
                           & fu_analysis_time(t3d), &
                           & fu_forecast_length(t3d), &
                           & fu_grid(t3d), &
                           & levels(l))
        call find_field_data_storage_2d(meteoMarketPtr, id, multi_time_stack_flag, TFP)
        if(fu_fails(.not.error,'Failed TFP field data pointer','dq_thermal_front_parameter'))return

        DO i = 1, fs
          IF (abs_nabla_T(i)==0.) THEN
            TFP(i) = 0.
          ELSE
            TFP(i) = -1*((dabsnablaT_dx(i)* dT_dx(i))+&
                   & (dabsnablaT_dy(i)* dT_dy(i)))/abs_nabla_T(i)
          END IF
          IF (error) RETURN
        END DO


!        CALL dq_store_2d(meteoMarketPtr, id, TFP, multi_time_stack_flag )
        IF (error) RETURN

      END DO loop_over_levels
    END DO loop_over_times

!    CALL free_all_work_arrays()
    call free_work_array(abs_nabla_T) 
!    call free_work_array(TFP)
    call free_work_array(dT_dx)
    call free_work_array(dT_dy)
    call free_work_array(dabsnablaT_dx)
    call free_work_array(dabsnablaT_dy)

  END SUBROUTINE dq_thermal_front_parameter


  !********************************************************

  SUBROUTINE dq_bulk_thermal_front_parameter(meteoMarketPtr, met_src, valid_times)

    ! Description:
    ! Calculates thermal front parameter (TFP) [K/m2] from layer thickness
    ! (1000 hPa ->)
    !
    ! TFP = vector product [(change of temperaturegradient) and
    !       (projetion in direction of the temperature gradient)]
    !
    ! All units: SI
    !
    ! Language: ANSI Fortran 90
    !
    ! Current code owner: Pilvi Siljamo, FMI 
    ! 
    IMPLICIT NONE
    
    ! Imported parameters with intent IN:
    type(meteo_data_source), INTENT(in) :: met_src
    TYPE(silja_time), DIMENSION(:), INTENT(in) :: valid_times
    type(mini_market_of_stacks), pointer :: meteoMarketPtr

    ! Local declarations:
    TYPE(silja_time) :: time
    TYPE(silja_3d_field), POINTER :: thickness3d
    TYPE(silja_field_id) :: id
    TYPE(silja_grid) :: scalar_grid
    TYPE(silja_level), DIMENSION(max_levels) :: levels
    INTEGER :: number_of_levels
    INTEGER :: l, fs, tloop, i
    REAL,DIMENSION(:), POINTER :: thickness
    REAL,DIMENSION(:), POINTER :: abs_nabla_Thick 
    REAL,DIMENSION(:), POINTER :: dThick_dy, dThick_dx
    REAL,DIMENSION(:), POINTER :: dabsnablaThick_dx, dabsnablaThick_dy
    REAL,DIMENSION(:), POINTER :: bulkTFP
    REAL, DIMENSION(:), POINTER :: pressure_bottom, pressure_top

    abs_nabla_Thick => fu_work_array() 
!    bulkTFP => fu_work_array()
    dThick_dx => fu_work_array()
    dThick_dy => fu_work_array()
    dabsnablaThick_dx => fu_work_array()
    dabsnablaThick_dy => fu_work_array()
    pressure_bottom  => fu_work_array()
    pressure_top => fu_work_array()
    IF (error) RETURN
    
    loop_over_times: DO tloop = 1, SIZE(valid_times)
      IF (.NOT.defined(valid_times(tloop))) EXIT loop_over_times
      time = valid_times(tloop)

      IF (fu_field_in_sm(meteoMarketPtr, met_src, & ! Already done before?
                       & bulk_tfp_flag, &
                       & time, &
                       & level_missing, &
                       & .true., &
                       & .false.)) CYCLE loop_over_times

      thickness3d => fu_sm_obstime_3d_field(meteoMarketPtr, met_src, &
                                          & temperature_flag, &
                                          & time,&
                                          & single_time)
      IF (error) RETURN

      CALL vertical_levels(thickness3d, levels, number_of_levels)
      IF (error) RETURN

      scalar_grid = fu_grid(thickness3d)
      IF (error) RETURN

      bulkTFP = 0.
      fs = fu_number_of_gridpoints(scalar_grid)

      CALL pressure_on_level(thickness3d, 1, pressure_bottom) !?onko alin?
      IF (error) RETURN

      loop_over_levels: DO l = 1, number_of_levels

        thickness => fu_grid_data_from_3d(thickness3d, l)
        IF (error) RETURN

        ! dThick_dx:
        CALL ddx_of_field(thickness, scalar_grid, scalar_grid, dThick_dx)
        IF (error) RETURN

        ! dThick_dy:
        CALL ddy_of_field(thickness, scalar_grid, scalar_grid, dThick_dy)
        IF (error) RETURN

        abs_nabla_Thick(1:fs) = SQRT(dThick_dx(1:fs)**2. + dThick_dy(1:fs)**2.)
        
        ! dabsnablaThick_dx 
        CALL ddx_of_field(abs_nabla_Thick(1:fs), scalar_grid,&
            & scalar_grid,dabsnablaThick_dx(1:fs))
        IF (error) RETURN

        ! dabsnablaThick_dy 
        CALL ddy_of_field(abs_nabla_Thick(1:fs), scalar_grid,&
            & scalar_grid,dabsnablaThick_dy(1:fs))
        IF (error) RETURN
                
        ! ------------------------------------------
        !
        ! Calculates thermal front parameter
        !
        ! ------------------------------------------


        CALL pressure_on_level(thickness3d, l, pressure_top) 
        IF (error) RETURN

        id = fu_set_field_id(fu_met_src(thickness3d), &
                           & bulk_tfp_flag, &
                           & fu_analysis_time(thickness3d), &
                           & fu_forecast_length(thickness3d), &
                           & fu_grid(thickness3d), &
                           & levels(l))
        call find_field_data_storage_2d(meteoMarketPtr, id, multi_time_stack_flag, bulkTFP)
        if(fu_fails(.not.error,'Failed bulkTFP field data pointer','dq_bulk_thermal_front_parameter'))return

        loop_over_gridpoints: DO i = 1, fs
          IF (abs_nabla_Thick(i)==0.) THEN
            bulkTFP(i) = 0.
          ELSE
            bulkTFP(i) = -1*(((dabsnablaThick_dx(i)* dThick_dx(i))+&
                &(dabsnablaThick_dy(i)* dThick_dy(i)))/abs_nabla_Thick(i)) *&
                & (-g/gas_constant_dryair)/&
                & LOG(pressure_bottom(i)/pressure_top(i))
          END IF
          IF (error) RETURN
        END DO loop_over_gridpoints

!        CALL dq_store_2d(meteoMarketPtr, id, bulkTFP, multi_time_stack_flag )
        IF (error) RETURN

      END DO loop_over_levels
    END DO loop_over_times

!    CALL free_all_work_arrays()
    call free_work_array(abs_nabla_Thick) 
!    call free_work_array(bulkTFP)
    call free_work_array(dThick_dx)
    call free_work_array(dThick_dy)
    call free_work_array(dabsnablaThick_dx)
    call free_work_array(dabsnablaThick_dy)
    call free_work_array(pressure_bottom)
    call free_work_array(pressure_top)

  END SUBROUTINE dq_bulk_thermal_front_parameter

  !************************************************************

  SUBROUTINE dq_bulk_richardson_number(meteoMarketPtr, met_src, valid_times, abl_param)
    !
    ! Returns the bulk Richarson number in Hirlam formulation
    ! with potential temperature (humidity dependent),
    ! directional shear and account for u_*
    !                    g        Theta_v(h1) - Theta_v(h_i)
    !   Rib_hir(h_i) =  ---     -----------------------------  
    !                 theta_v    DealtaU(h1,h_i)^2 + 100*ustar
    ! 
    ! All units: SI
    !
    IMPLICIT NONE
    !
    ! Imported parameters with intent(in):
    type(meteo_data_source), INTENT(in) :: met_src
    TYPE(silja_time), DIMENSION(:), INTENT(in) :: valid_times
    type(mini_market_of_stacks), pointer :: meteoMarketPtr
    integer, intent(in) :: abl_param

    ! Imported parameters with intent(out):

    ! Local declarations:
    TYPE(silja_time) :: time
    INTEGER :: tloop, i, lev
    TYPE(silja_grid) :: scalar_grid
    TYPE(silja_field_id) :: id
    type(silja_field), pointer :: fldTmp
    TYPE(silja_3d_field), POINTER :: theta3d, q3d
    TYPE(silja_3d_field), POINTER :: u3d,v3d, height3d
    REAL, DIMENSION(:), POINTER :: bulk_richardson_number
    REAL, DIMENSION(:), POINTER :: height2, height1
    REAL, DIMENSION(:), POINTER :: theta_1, theta_2, q_1, q_2
    REAL, DIMENSION(:), POINTER :: u_1,v_1, u_2, v_2, ustar
    TYPE(silja_level), DIMENSION(max_levels) :: levels
    INTEGER :: number_of_levels, iTmp, ix, iy, ix1, iy1, iTmp1, iCount


!    bulk_richardson_number => fu_work_array() 
!    IF (error) RETURN

    DO tloop = 1, SIZE(valid_times)

      IF (.NOT.defined(valid_times(tloop))) EXIT   !loop_over_times
      time = valid_times(tloop)

      IF (fu_field_in_sm(meteoMarketPtr, met_src,&   ! Already done before?
                       & bulk_richardson_nbr_flag,&
                       & time,&
                       & level_missing,&
                       & .true.,&
                       & .false.)) CYCLE  ! loop_over_times

      height3d => fu_sm_obstime_3d_field(meteoMarketPtr, met_src, height_flag, time, single_time)
      IF (error) RETURN

      theta3d => fu_sm_obstime_3d_field(meteoMarketPtr, met_src, potential_temperature_flag, &
                                      & time, single_time)
      IF (error) RETURN
      
      q3d => fu_sm_obstime_3d_field(meteoMarketPtr, met_src, specific_humidity_flag, &
                                      & time, single_time)
      IF (error .or. .not. associated(q3d))then
        call unset_error('dq_bulk_richardson_number')
        nullify(q3d)
      endif

      u3d => fu_sm_obstime_3d_field(meteoMarketPtr, met_src, u_flag, time, single_time)
      v3d => fu_sm_obstime_3d_field(meteoMarketPtr, met_src, v_flag, time, single_time)
      IF (error) RETURN

      CALL vertical_levels(theta3d, levels, number_of_levels)
      IF (error) RETURN

      ! lowest model level
      u_1 => fu_grid_data_from_3d(u3d, 1)
      v_1 => fu_grid_data_from_3d(v3d, 1)
      if(abl_param == abl_full_param) q_1 => fu_grid_data_from_3d(q3d, 1)
      theta_1 => fu_grid_data_from_3d(theta3d, 1)
      height1 => fu_grid_data_from_3d(height3d, 1)
      IF (error) RETURN
      !
      ! Finally, u*
      !
      fldTmp => fu_sm_obstime_field(meteoMarketPtr, met_src, friction_velocity_flag, level_missing, time, forwards)
      if(error)return
      ustar => fu_grid_data(fldTmp)
     
      !
      ! Having the initial stuff set, start the layer cycle
      !
      DO lev = 1, number_of_levels

        id = fu_set_field_id(fu_met_src(theta3d),&
                           & bulk_richardson_nbr_flag,&
                           & fu_analysis_time(theta3d),&
                           & fu_forecast_length(theta3d), &
                           & fu_grid(theta3d),&
                           & levels(lev))
        call find_field_data_storage_2d(meteoMarketPtr, id, multi_time_stack_flag, bulk_richardson_number)
        if(fu_fails(.not.error,'Failed bulk_richardson_number field data pointer','dq_bulk_richardson_number'))return
        
        if(lev == 1)then        ! nothing for the first level
          bulk_richardson_number(1:fs_meteo) = 0.
          cycle
        endif
        
        u_2 => fu_grid_data_from_3d(u3d, lev)
        v_2 => fu_grid_data_from_3d(v3d, lev)
        if(abl_param == abl_full_param)q_2 => fu_grid_data_from_3d(q3d, lev)
        theta_2 => fu_grid_data_from_3d(theta3d, lev)
        height2 => fu_grid_data_from_3d(height3d, lev)
        IF (error) RETURN


        call bulk_richardson_nbr_one_level(u_1, u_2, v_1, v_2, q_1, q_2, theta_1, theta_2, &
                                         & height1, height2, uStar, abl_param == abl_full_param, &
                                         & bulk_richardson_number)
!        call msg('level',lev)
!        call msg('Ri min', real_value=MINVAL(bulk_richardson_number(1:fs)))
!        call msg('Ri max', real_value=MAXVAL(bulk_richardson_number(1:fs)))
!        do iTmp=1, fs
!          if(bulk_richardson_number(iTmp) > 1000)then
!            call msg('Problem Ri at i=',iTmp,bulk_richardson_number(iTmp))
!          endif
!        enddo

!        CALL dq_store_2d(meteoMarketPtr, id, bulk_richardson_number, multi_time_stack_flag )
        IF (error) RETURN 


      END DO ! loop_over_levels

    END DO ! loop_over_times

!    CALL free_all_work_arrays()
!    call free_work_array(bulk_richardson_number)

  END SUBROUTINE dq_bulk_richardson_number

  !************************************************************

  SUBROUTINE dq_gradient_richardson_number(meteoMarketPtr, met_src, valid_times, abl_param)

    !
    ! Returns the  Richarson number 
    !  it is 
    !  zero     for the 1st level 
    !  1 -- 2        for the 2nd level
    !  2 -- 4        for the 3-rd level
    !  ...
    !  n-2 -- n-1      for (n-1)-th level
    !  n-1 -- n        for n-th level
    !
    ! with potential temperature (humidity dependent),
    ! directional shear 
    !                    g        Theta(h_{i-1}) - Theta(h_i)
    !        Ri(h_i) = -------    -----------------------------  
    !                  theta_v         DealtaU(h_{i-1},h_i)^2
    ! 
    ! All units: SI
    !
    IMPLICIT NONE
    !
    ! Imported parameters with intent(in):
    type(meteo_data_source), INTENT(in) :: met_src
    TYPE(silja_time), DIMENSION(:), INTENT(in) :: valid_times
    type(mini_market_of_stacks), pointer :: meteoMarketPtr
    integer, intent(in) :: abl_param

    ! Imported parameters with intent(out):

    ! Local declarations:
    TYPE(silja_time) :: time
    INTEGER :: tloop, i, lev
    TYPE(silja_grid) :: scalar_grid
    TYPE(silja_field_id) :: id
    type(silja_field), pointer :: fldTmp
    TYPE(silja_3d_field), POINTER :: theta3d, q3d
    TYPE(silja_3d_field), POINTER :: u3d,v3d, height3d
    REAL, DIMENSION(:), POINTER :: richardson_number
    REAL, DIMENSION(:), POINTER :: height2, height1
    REAL, DIMENSION(:), POINTER :: theta_1, theta_2, q_1, q_2
    REAL, DIMENSION(:), POINTER :: u_1,v_1, u_2, v_2, zero
    TYPE(silja_level), DIMENSION(max_levels) :: levels
    INTEGER :: number_of_levels, iTmp, ix, iy, ix1, iy1, iTmp1, iCount


!    richardson_number => fu_work_array(fs_meteo)
    zero => fu_work_array(fs_meteo) 
    zero(1:fs_meteo) = 0.
    IF (error) RETURN

    DO tloop = 1, SIZE(valid_times)

      IF (.NOT.defined(valid_times(tloop))) EXIT   !loop_over_times
      time = valid_times(tloop)

      IF (fu_field_in_sm(meteoMarketPtr, met_src,&   ! Already done before?
                       & gradient_richardson_nbr_flag,&
                       & time,&
                       & level_missing,&
                       & .true.,&
                       & .false.)) CYCLE  ! loop_over_times

      height3d => fu_sm_obstime_3d_field(meteoMarketPtr, met_src, height_flag, time, single_time)
      IF (error) RETURN

      theta3d => fu_sm_obstime_3d_field(meteoMarketPtr, met_src, potential_temperature_flag, &
                                      & time, single_time)
      IF (error) RETURN
      
      if(abl_param == abl_full_param)then
        q3d => fu_sm_obstime_3d_field(meteoMarketPtr, met_src, specific_humidity_flag, time, single_time)
        IF (error .or. .not. associated(q3d))then
          call unset_error('dq_bulk_richardson_number')
          nullify(q3d)
          nullify(q_2) 
        endif
      endif

      u3d => fu_sm_obstime_3d_field(meteoMarketPtr, met_src, u_flag, time, single_time)
      v3d => fu_sm_obstime_3d_field(meteoMarketPtr, met_src, v_flag, time, single_time)
      IF (error) RETURN

      CALL vertical_levels(theta3d, levels, number_of_levels)
      IF (error) RETURN


      !
      ! Having the initial stuff set, start the layer cycle
      !
      DO lev = 1, number_of_levels
        u_2 => fu_grid_data_from_3d(u3d, lev)
        v_2 => fu_grid_data_from_3d(v3d, lev)
        if (abl_param == abl_full_param) q_2 => fu_grid_data_from_3d(q3d, lev)
        theta_2 => fu_grid_data_from_3d(theta3d, lev)
        height2 => fu_grid_data_from_3d(height3d, lev)
        IF (error) RETURN
        id = fu_set_field_id(fu_met_src(theta3d),&
                           & gradient_richardson_nbr_flag,&
                           & fu_analysis_time(theta3d),&
                           & fu_forecast_length(theta3d), &
                           & fu_grid(theta3d),&
                           & levels(lev))
        call find_field_data_storage_2d(meteoMarketPtr, id, multi_time_stack_flag, richardson_number)
        if(fu_fails(.not.error,'Failed richardson_number field data pointer','dq_gradient_richardson_number'))return

        if (lev > 1) then
           call bulk_richardson_nbr_one_level(u_1, u_2, v_1, v_2, q_1, q_2, theta_1, theta_2, &
                                         & height1, height2, zero, abl_param == abl_full_param, &
                                         & richardson_number)
        else 
            richardson_number(1:fs_meteo) = 0. ! lowest level uses MoninObukhov anyhow
        endif
!        call msg('level',lev)
!        call msg('Ri min', real_value=MINVAL(bulk_richardson_number(1:fs)))
!        call msg('Ri max', real_value=MAXVAL(bulk_richardson_number(1:fs)))
!        do iTmp=1, fs
!          if(bulk_richardson_number(iTmp) > 1000)then
!            call msg('Problem Ri at i=',iTmp,bulk_richardson_number(iTmp))
!          endif
!        enddo

!        CALL dq_store_2d(meteoMarketPtr, id, richardson_number, multi_time_stack_flag )
        IF (error) RETURN 
        u_1 => u_2
        v_1 => v_2
        q_1 => q_2
        theta_1 => theta_2 
        height1 => height2

      END DO ! loop_over_levels

    END DO ! loop_over_times

!    call free_work_array(richardson_number)
    call free_work_array(zero)

  END SUBROUTINE dq_gradient_richardson_number

  !****************************************************************
  
  subroutine bulk_richardson_nbr_one_level(u_1, u_2, v_1, v_2, q_1, q_2, theta_1, theta_2, &
                                         & height_1, height_2, uStar, ifFullABLParam, &
                                         & bulk_richardson_number)
    !
    ! Calculates the bulk Richardson number for the given level
    !
    implicit none

    !
    ! Imported parameters
    real, dimension(:), pointer :: u_1, u_2, v_1, v_2, q_1, q_2, theta_1, theta_2, height_1, height_2
    real, dimension(:), pointer :: uStar, bulk_richardson_number
    logical, intent(in) :: ifFullABLParam

    !
    ! Local variables
    integer :: i, ix, iy, iErrorCount, iCount, ix1, iy1, i1
    real :: du, dv, deltaU2, deltaT, fTmp
    real, parameter :: eps = molecular_weight_air/molecular_weight_water - 1.
    

    iErrorCount = 0
    if( .not. associated(q_1))then
            call msg("Warning! Making Richardson without humidity!!!")
    endif


    do i = 1, fs_meteo

      du = u_2(i) - u_1(i)
      dv = v_2(i) - v_1(i)
      deltaU2 = du*du + dv*dv + 100*ustar(i)*ustar(i) + 0.01 ! as in Hirlam + 0.01
                                        !originally  ustar(i) = amax1(ustar(i),0.01)
                                        ! deltaU2 > 0.01  in any case
                                        ! Check for small deltaU2 removed (RK)
      if(ifFullABLParam)then
        deltaT = theta_2(i)*(1.-eps*q_2(i)) - theta_1(i)*(1.-eps*q_1(i))
      else
        deltaT = theta_2(i) - theta_1(i)
      endif
      bulk_richardson_number(i)= 2. * g * (height_2(i)-height_1(i)) * deltaT / &
                                        & ( (theta_1(i)+theta_2(i))* deltaU2)

      if(bulk_richardson_number(i) > 100.)then
        bulk_richardson_number(i) = 99.99
        iErrorCount = iErrorCount + 1
      elseif(bulk_richardson_number(i) < -100.)then
        bulk_richardson_number(i) = -99.99
        iErrorCount = iErrorCount + 1
      endif
    end do  ! cycle over the grid
    !
    ! Have to check for the cases with extremely large positive and negative Ri.
    ! These are in most cases just the results of some numerics. So, we will set them using the
    ! neighbouring cells.
    !
    if(iErrorCount == 0)return

    call msg("Patching Richardson for Level :", height_2(1))

    do iy = 1, ny_meteo
      do ix = 1, nx_meteo
        i = ix + (iy-1)* nx_meteo
        if(abs(bulk_richardson_number(i)) >= 10. )then
          !
          ! Bad Ri
          !
          fTmp = 0.0
          iCount = 0
          do iy1 = max(iy-1,1),min(ny_meteo,iy+1)
            do ix1 = max(ix-1,1),min(nx_meteo,ix+1)
              i1 = ix1 + (iy1-1)* nx_meteo
              if( abs(bulk_richardson_number(i1)) < 10. )then
                !
                ! Found good one 
                !
                fTmp = fTmp + bulk_richardson_number(i1)
                iCount = iCount + 1
              endif
            end do
          end do
          if(iCount > 0)then
            bulk_richardson_number(i) = fTmp / real(iCount)  ! Set the average
          else
            bulk_richardson_number(i) = sign(1.,bulk_richardson_number(i))  ! Anything...
          endif
          if(iErrorCount == 1)return
          iErrorCount = iErrorCount - 1 
        endif  ! bad Ri
      end do  ! ix
    end do  ! iy

  end subroutine bulk_richardson_nbr_one_level


  !************************************************************

  SUBROUTINE dq_brunt_vaisala_freq(meteoMarketPtr, met_src, valid_times)
    !
    ! Returns Brunt-Vaisala frequency when potential
    ! temperature, air temperature at the two level are
    ! known. Hydrostatic approximation is used for calculating
    ! thickness of atmospheric layer between pressure1 and
    ! pressure2.
    !
    ! Index 1 refers to values at lower level, 2 the higher level and
    ! layer thickness must be positive.
    ! 
    ! All units: SI
    !
    IMPLICIT NONE
    !
    ! Imported parameters with intent(in):
    type(meteo_data_source), INTENT(in) :: met_src
    TYPE(silja_time), DIMENSION(:), INTENT(in) :: valid_times
    type(mini_market_of_stacks), pointer :: meteoMarketPtr

    ! Imported parameters with intent(out):

    ! Local declarations:
    TYPE(silja_time) :: time
    INTEGER :: tloop, i, lev
    TYPE(silja_grid) :: scalar_grid
    TYPE(silja_field_id) :: id
    type(silja_field), pointer :: fldTmp
    TYPE(silja_3d_field), POINTER :: theta3d, height3d
    REAL, DIMENSION(:), POINTER :: BV_fr, &
                                 & height1, height2, height3, &
                                 & theta1, theta2, theta3, &
                                 & dtheta_dz
    TYPE(silja_level), DIMENSION(max_levels) :: levels
    INTEGER :: number_of_levels, iCell
    real :: fTmp

!    BV_fr => fu_work_array() 
    dtheta_dz => fu_work_array() 
    IF (error) RETURN

    DO tloop = 1, SIZE(valid_times)

      IF (.NOT.defined(valid_times(tloop))) EXIT   !loop_over_times
      time = valid_times(tloop)

      IF (fu_field_in_sm(meteoMarketPtr, met_src,&   ! Already done before?
                       & brunt_vaisala_freq_flag,&
                       & time,&
                       & level_missing,&
                       & .true.,&
                       & .false.)) CYCLE  ! loop_over_times

      height3d => fu_sm_obstime_3d_field(meteoMarketPtr, met_src, height_flag, time, single_time)
      IF (error) RETURN

      theta3d => fu_sm_obstime_3d_field(meteoMarketPtr, met_src, potential_temperature_flag, &
                                      & time, single_time)
      IF (error) RETURN

      CALL vertical_levels(theta3d, levels, number_of_levels)
      IF (error) RETURN

      !
      ! Get the screen-level fields for the lowest-level Ri and
      ! make appropriate attributions to the heights
      !
      ! Having the initial stuff set, start the layer cycle
      !
      DO lev = 1, number_of_levels

        if(lev == 1)then
          !
          ! Having height, select the data for determining the profiles: 10m wind 
          ! and 1st model level or just 1st and 2d model levels for everything.
          !
          height2 => fu_grid_data_from_3d(height3d, 1)
          fTmp = sum(height2(1:fs_meteo))/fs_meteo

          if(fTmp <= 3.)then
            !
            ! First model level is dangerously close to 2m, skip the screen-level wind and use 3D only
            ! Accuracy is then 1st order
            !
            theta2 => fu_grid_data_from_3d(theta3d, 1)
            theta3 => fu_grid_data_from_3d(theta3d, 2)
            do iCell=1, fs_meteo
              dtheta_dz(iCell) = (theta3(iCell) - theta2(iCell)) / (height3(iCell) - height2(iCell))
            enddo
          else
            !
            ! Can use 2m and first model level to get the second order of accuracy for theta
            !
            fldTmp => fu_sm_obstime_field(meteoMarketPtr, met_src, potential_temperature_2m_flag, &
                                        & level_missing, time, forwards)
            if(error)return
            theta1 => fu_grid_data(fldTmp)

            theta2 => fu_grid_data_from_3d(theta3d, 1)
            height2 => fu_grid_data_from_3d(height3d, 1)
            theta3 => fu_grid_data_from_3d(theta3d, 2)
            height3 => fu_grid_data_from_3d(height3d, 2)
            do iCell=1, fs_meteo
              dtheta_dz(iCell) = 1. / (2. - height3(iCell)) * &
                               & ((theta2(iCell) - theta3(iCell)) / (height2(iCell) - height3(iCell)) * &
                                & (2. - height2(iCell)) + &
                                & (theta1(iCell) - theta2(iCell)) / (2. - height2(iCell)) * &
                                & (height2(iCell) - height3(iCell)) &
                                & )
            end do

          endif

        elseif(lev == number_of_levels)then
          !
          ! Last model level - again 1st order of accuracy
          !
          height1 => fu_grid_data_from_3d(height3d, lev-1)
          height2 => fu_grid_data_from_3d(height3d, lev)
          theta1 => fu_grid_data_from_3d(theta3d, lev-1)
          theta2 => fu_grid_data_from_3d(theta3d, lev)
          do iCell=1, fs_meteo
            dtheta_dz(iCell) = (theta2(iCell) - theta1(iCell)) / (height2(iCell) - height1(iCell))
          end do

        else
          !
          ! General case: second-order for inhomogenous grid
          !
          height1 => fu_grid_data_from_3d(height3d, lev-1)
          height2 => fu_grid_data_from_3d(height3d, lev)
          height3 => fu_grid_data_from_3d(height3d, lev+1)

          theta1 => fu_grid_data_from_3d(theta3d, lev-1)
          theta2 => fu_grid_data_from_3d(theta3d, lev)
          theta3 => fu_grid_data_from_3d(theta3d, lev+1)

          do iCell=1, fs_meteo
            dtheta_dz(iCell) = 1. / (height1(iCell) - height3(iCell)) * &
                             & ((theta2(iCell) - theta3(iCell)) / (height2(iCell) - height3(iCell)) * &
                              & (height1(iCell) - height2(iCell)) + &
                              & (theta1(iCell) - theta2(iCell)) / (height1(iCell) - height2(iCell)) * &
                              & (height2(iCell) - height3(iCell)) &
                             & )
          end do

        endif  ! Level-dependent wind and temperature gradients computation
        IF (error) RETURN

        id = fu_set_field_id(fu_met_src(theta3d),&
                           & brunt_vaisala_freq_flag,&
                           & fu_analysis_time(theta3d),&
                           & fu_forecast_length(theta3d), &
                           & fu_grid(theta3d),&
                           & levels(lev))
        call find_field_data_storage_2d(meteoMarketPtr, id, multi_time_stack_flag, BV_fr)
        if(fu_fails(.not.error,'Failed BV_fr field data pointer','dq_brunt_vaisala_freq'))return

        do iCell=1, fs_meteo
          BV_fr(iCell) = g * dtheta_dz(iCell) / theta2(iCell)
        end do  ! cycle over the grid

!        call msg('level',lev)
!        call msg('Ri min', real_value=MINVAL(bulk_richardson_number(1:fs)))
!        call msg('Ri max', real_value=MAXVAL(bulk_richardson_number(1:fs)))
!        do iTmp=1, fs
!          if(bulk_richardson_number(iTmp) > 1000)then
!            call msg('Problem Ri at i=',iTmp,bulk_richardson_number(iTmp))
!          endif
!        enddo

!        CALL dq_store_2d(meteoMarketPtr, id, BV_fr, multi_time_stack_flag )
        IF (error) RETURN 

      END DO ! loop_over_levels

    END DO ! loop_over_times

!    CALL free_all_work_arrays()
!    call free_work_array(BV_fr) 
    call free_work_array(dtheta_dz) 

  END SUBROUTINE dq_brunt_vaisala_freq


  !************************************************************

  SUBROUTINE dq_ABL_params(meteoMarketPtr, met_src, valid_times, abl_param, ABL_height_method_switch)
    
    ! This routine for computing of most important ABL and similarity theory parameters :
    ! - Monin Obukhov length (L) [m]
    ! - friction velocity (u*) [m/s]
    ! - turbulent temperature scale (T*) [K]
    ! - surface heat flux (H0)  [W/m2]
    ! - turbulent velocity scale omega* [m/s]
    ! - ABL height
    !
    ! Method follows the Genikhovich's idea for 2-level U, T values.
    ! ABL height calculations are defined in the standard setup file
    !
    ! Index 1 refers to values at 1st model layer. Also 2m temperature and 
    ! 10m wind are used. Since ECMWF has the 1-st level at 10m, in such a case
    ! we have to use 1-st and 2-d model levels. Notations are kept the same.
    ! The same trick is made for temperature (no immediate use, however)
    !
    ! Note that the convective velocity scale alters friction velocity. Therefore,
    ! three steps are made below
    ! 1. Main similarity parameters for surface layer: Kz, u*, T*, L, Q*
    ! 2. ABL height
    ! 3. Convective velocity scale, update of u*
    !
    IMPLICIT NONE
    !
    ! Imported parameters with intent(in):

    type(meteo_data_source), INTENT(in) :: met_src
    TYPE(silja_time), DIMENSION(:), INTENT(in) :: valid_times
    integer, intent(in) :: abl_param, ABL_height_method_switch
    type(mini_market_of_stacks), pointer :: meteoMarketPtr

    ! Imported parameters with intent(out):

    ! Local declarations:
    TYPE(silja_time) :: time
    INTEGER :: tloop, i, j, iu,it, i_tot, i_err, i_err_pos
    TYPE(silja_grid) :: scalar_grid
    TYPE(silja_field_id) :: id
    TYPE(silja_3d_field), POINTER :: height3d,u3d,v3d,t3d, q3d,  theta3d
    TYPE(silja_field), POINTER :: windspeed_10m, field_temp_2m, field_q_2m, ps_2d, &
                                & fraction_of_land_field
    TYPE(silja_level), DIMENSION(max_levels) :: levels
    REAL, DIMENSION(:), POINTER :: monin_obukhov_inv, friction_velocity, & 
     & temperature_scale, sens_heat_flux, height_lev1, height_lev2, &
     & ps, temp2m, q2m, q_1, density, abl_height, wind10m,&
     & temp_1, x1, Kz, u_1, v_1, u_10m, v_10m, land_mask, alpha_Pr, conv_vel_scale, &
     & latent_heat_flux, humid_scale
    REAL :: windspeed_1, fMinAlert, fMinForce, fMaxForce, fMaxAlert
    !
    ! Local parameters:
    ! sigma is unitless constant from Genikhovich' paper
    ! alpha_pr=1./Pr=1/0.74=1.35 
    ! height_K is altitude at which Kz is taken
    ! gamma_dry is dry adiabatic gradient =0.01 K/m
    ! Zilitinkevich: sigma_stable=12, sigma_unstable=5
    ! lev1_height_stable - height of the surface constant-flux layer in stable case
    !
    REAL :: sigma=10.,sigma_stable=12.,sigma_unstable=5., alpha_Pr_neutral=1.35, gamma_dry=0.01
    real :: lev1_Theight_stable = 20. ! thin stable surface layer for T, but do not approach 2m 
    ! Wet adiabatic correction depends on temperature. Below are values from -20 till +40 deg
    real, dimension(7) :: gamma_wet = (/0.00856, 0.00763, 0.00658, 0.00532, 0.00435, 0.00363, 0.00315 /)
    real :: fMeanHeightLev1, height_T2m, height_u10m ! height of the 1st model level, true 2m and 10m heights
    real :: fMinHeightLev1
    logical :: ifWind10m, ifTemp2m

    real :: lnzT, lnzu, Cp_ro, vap_heat, height_1u, height_1T, du, dT, HHH_K, gamma, &
          & fric_vel_cube, arg, Ri, Kz_scaling, y_solver, fTmp, dzu, dzT
          
    integer :: iTmp

!    INTEGER, DIMENSION(100,100)::prob, prob_err, prob_err_pos!
!    REAL, PARAMETER:: du_min=0., du_max=7, dt_min=-3, dt_max=4

!    monin_obukhov_inv => fu_work_array() 
!    friction_velocity => fu_work_array() 
!    temperature_scale => fu_work_array() 
!    sens_heat_flux => fu_work_array() 
    density => fu_work_array()
    x1 =>fu_work_array()
!    Kz =>fu_work_array()
!    alpha_Pr =>fu_work_array()
!    abl_height => fu_work_array()
!    conv_vel_scale => fu_work_array()
    
    if (abl_param == abl_dry_param) then
      q2m => fu_work_array()
      q_1 => fu_work_array()
!    else
!      latent_heat_flux => fu_work_array()
!      humid_scale => fu_work_array()
    end if
    if(error)return

    fraction_of_land_field => fu_sm_simple_field(meteoMarketPtr, met_src_missing,&
                                               & fraction_of_land_flag,&
                                               & level_missing, &
                                               & single_time_stack_flag)
    IF (error) RETURN
    land_mask => fu_grid_data(fraction_of_land_field) ! set in physiography
    IF (error) RETURN

    call quantity_feasible_range(SILAM_sensible_heat_flux_flag, &
                               & fMinAlert, fMinForce, fMaxForce, fMaxAlert)

!    open(90,file ='probability.bin',access='direct',form='unformatted',recl=10000)
!    open(91,file ='probability.txt',access='append',form='formatted')

    !------------------------------------------------------------------------------------
    !
    !  Main loop over times
    !
    loop_over_times: DO tloop = 1, SIZE(valid_times)

!      prob=0
!      prob_err=0
!      prob_err_pos=0
!      i_tot=0
!      i_err=0
!      i_err_pos=0
 
      IF (.NOT.defined(valid_times(tloop))) EXIT loop_over_times
      time = valid_times(tloop)

      IF (fu_field_in_sm(meteoMarketPtr, met_src, & ! Already done before?
                       & MO_length_inv_flag,&
                       & time,&
                       & level_missing,&
                       & .false.,&    ! search for 3D field
                       & .false.)) CYCLE loop_over_times
      
      ps_2d => fu_sm_obstime_field(meteoMarketPtr, met_src,&
                                 & ground_pressure_flag,&
                                 & level_missing,&
                                 & time,&
                                 & single_time)
      IF (error) RETURN
      ps => fu_grid_data(ps_2d)

      !------------------------------------------------------------------------------------
      !
      !  Basic similarity parameters for surface layer
      !

      height3d => fu_sm_obstime_3d_field(meteoMarketPtr, met_src,&
                                       & height_flag, &
                                       & time,&
                                       & single_time)
      IF (error) RETURN
      height_lev1 => fu_grid_data_from_3d(height3d, 1)
      height_lev2 => fu_grid_data_from_3d(height3d, 2)
      IF (error) RETURN
      !
      ! Having height, select the data for determining the profiles: 10m wind and 2m temperature
      ! and 1st model level or just 1st and 2d model levels for everything.
      !
      fMeanHeightLev1 = sum(height_lev1(1:fs_meteo))/fs_meteo
      fMinHeightLev1 = minval(height_lev1(1:fs_meteo))

!      call msg('Mean 1-st level height=',real_value=fMeanHeightLev1)

      ! Need some margin for gradients!!!! 
      if(fMeanHeightLev1 <= 5.)then  
        !
        ! The first model level is too close to 2m, both temperature and wind have to be taken 
        ! from model levels
        !
        ifTemp2m = .false.
        ifWind10m = .false.
      elseif(fMeanHeightLev1 <= 17. .or. fMinHeightLev1 < 15.)then
      !elseif(fMeanHeightLev1 <= 11. .or. fMinHeightLev1 < 10.1)then ! Not suitable for gradients
        ifTemp2m = .true.
        ifWind10m = .false.
        if(any(height_lev1(1:fs_meteo) < 2.))then
          call set_error('First model level is below 2m','dq_monin_obukhov')
          return
        endif
      else
        ifTemp2m = .true.
        ifWind10m = .true.
        if(any(height_lev1(1:fs_meteo) < 10.))then
          call set_error('First model level is below 10m','dq_monin_obukhov')
          return
        endif
      endif
!      print *, 'ifTemp2m, ifWind10m =', ifTemp2m, ifWind10m
      !
      ! Get the temperature from appropriate fields
      !
      t3d => fu_sm_obstime_3d_field(meteoMarketPtr, met_src,&
                                  & temperature_flag, &
                                  & time,&
                                  & single_time)
      IF (error) RETURN

      if(ifTemp2m)then
        field_temp_2m => fu_sm_obstime_field(meteoMarketPtr, met_src,&
                                           & temperature_2m_flag,&
                                           & level_missing,&
                                           & time,&
                                           & single_time)
        if(error)return
        temp2m => fu_grid_data(field_temp_2m)
        temp_1 => fu_grid_data_from_3d(t3d, 1)
      else
        temp2m => fu_grid_data_from_3d(t3d, 1)
        temp_1 => fu_grid_data_from_3d(t3d, 2)
      endif
      IF (error) RETURN
      !
      ! Get the wind from appropriate fields
      !
      u3d => fu_sm_obstime_3d_field(meteoMarketPtr, met_src,&
                                          & u_flag, &
                                          & time,&
                                          & single_time)
      if(error)return
      v3d => fu_sm_obstime_3d_field(meteoMarketPtr, met_src,&
                                          & v_flag, &
                                          & time,&
                                          & single_time)
      if(error)return

      if(ifWind10m)then
        windspeed_10m => fu_sm_obstime_field(meteoMarketPtr, met_src,&
                                           & windspeed_10m_flag,&
                                           & level_missing,&
                                           & time,&
                                           & single_time)
        IF (error) RETURN
        wind10m => fu_grid_data(windspeed_10m)
        u_1 => fu_grid_data_from_3d(u3d, 1)
        v_1 => fu_grid_data_from_3d(v3d, 1)
      else
        u_10m => fu_grid_data_from_3d(u3d, 1)
        v_10m => fu_grid_data_from_3d(v3d, 1)
        u_1 => fu_grid_data_from_3d(u3d, 2)
        v_1 => fu_grid_data_from_3d(v3d, 2)
      endif
      
      theta3d => fu_sm_obstime_3d_field(meteoMarketPtr, met_src,&
                                      & potential_temperature_flag, &
                                      & time,&
                                      & single_time)
      IF (error) RETURN

      !
      ! Find/make places for all fields
      !
      id = fu_id(fu_field_from_3d_field(t3d, 1)) !Some non-staggered 

      !id = fu_id(field_temp_2m) !!!WRONG! field_temp_2m might be missing!
      call set_level(id, ground_level)
      
      call set_quantity(id, MO_length_inv_flag)
      call find_field_data_storage_2d(meteoMarketPtr, id, multi_time_stack_flag, monin_obukhov_inv)
      if(fu_fails(.not.error,'Failed MO_length_inv_flag field data pointer','dq_monin_obukhov'))return

      call set_quantity(id, Prandtl_nbr_flag)
      call find_field_data_storage_2d(meteoMarketPtr, id, multi_time_stack_flag, alpha_Pr)
      if(fu_fails(.not.error,'Failed Prandtl_nbr_flag field data pointer','dq_monin_obukhov'))return

      call set_quantity(id, friction_velocity_flag)
      call find_field_data_storage_2d(meteoMarketPtr, id, multi_time_stack_flag, friction_velocity)
      if(fu_fails(.not.error,'Failed friction_velocity_flag field data pointer','dq_monin_obukhov'))return

      call set_quantity(id, temperature_scale_flag)
      call find_field_data_storage_2d(meteoMarketPtr, id, multi_time_stack_flag, temperature_scale)
      if(fu_fails(.not.error,'Failed temperature_scale_flag field data pointer','dq_monin_obukhov'))return

      call set_quantity(id, SILAM_sensible_heat_flux_flag)
      call find_field_data_storage_2d(meteoMarketPtr, id, multi_time_stack_flag, sens_heat_flux)
      if(fu_fails(.not.error,'Failed SILAM_sensible_heat_flux_flag field data pointer','dq_monin_obukhov'))return

      call set_quantity(id, Kz_scalar_1m_flag)
      call find_field_data_storage_2d(meteoMarketPtr, id, multi_time_stack_flag, Kz)
      if(fu_fails(.not.error,'Failed Kz_scalar_1m_flag field data pointer','dq_monin_obukhov'))return

      call set_quantity(id, convective_velocity_scale_flag)
      call find_field_data_storage_2d(meteoMarketPtr, id, multi_time_stack_flag, conv_vel_scale)
      if(fu_fails(.not.error,'Failed convective_velocity_scale_flag field data pointer','dq_monin_obukhov'))return

      call set_quantity(id, abl_height_m_flag)
      call find_field_data_storage_2d(meteoMarketPtr, id, multi_time_stack_flag, abl_height)
      if(fu_fails(.not.error,'Failed abl_height_m_flag field data pointer','dq_monin_obukhov'))return

      if(abl_param == abl_full_param)then ! Full-blown parameterization, with wetness correction
        call set_quantity(id, silam_latent_heat_flux_flag)
        call find_field_data_storage_2d(meteoMarketPtr, id, multi_time_stack_flag, latent_heat_flux)
        if(fu_fails(.not.error,'Failed silam_latent_heat_flux_flag field data pointer','dq_monin_obukhov'))return

        call set_quantity(id, humidity_scale_flag)
        call find_field_data_storage_2d(meteoMarketPtr, id, multi_time_stack_flag, humid_scale)
        if(fu_fails(.not.error,'Failed humidity_scale_flag field data pointer','dq_monin_obukhov'))return
      end if
      !
      ! ABL parameterization allowed here is abl_full_param and abl_dry_param
      ! Should the dry ABL is used, humidity correction is simply set to zero
      !
      select case(abl_param)

        case(abl_full_param) ! Full-blown parameterization, with wetness correction
          q3d => fu_sm_obstime_3d_field(meteoMarketPtr, met_src,&
                                      & specific_humidity_flag, &
                                      & time,&
                                      & single_time)
          if(ifTemp2m)then
            field_q_2m => fu_sm_obstime_field(meteoMarketPtr, met_src,&
                                            & specific_humidity_2m_flag,&
                                            & level_missing,&
                                            & time,&
                                            & single_time)
            if(error)return
            q2m => fu_grid_data(field_q_2m)
            q_1 => fu_grid_data_from_3d(q3d, 1)
          else
            q2m => fu_grid_data_from_3d(q3d, 1)
            q_1 => fu_grid_data_from_3d(q3d, 2)
          endif

        case(abl_dry_param)  ! Humidity is ignored
          q2m(1:fs_meteo) = 0.
          q_1(1:fs_meteo) = 0.
          nullify(q3d)

        case default
          call set_error('Only full_abl_param and dry_abl_param so far','dq_monin_obukhov')
          return
      end select

      !--------------------------------------------------------
      !
      ! The main cycle over the grid starts

      HHH_K=Kz_ref_height * Kz_ref_height * Kz_ref_height
!      Cp_ro = specific_heat_dryair*density_air

      Kz_scaling = 1.0   !2.0
!      call msg('============= Kz scaling in dq_monin_obukhov:', Kz_scaling)

      DO i=1,fs_meteo

        Cp_ro= specific_heat_dryair*ps(i)/(gas_constant_dryair*temp2m(i))  !Cp*ro

        ! To cope grid-border effects - we have to explicitly cut the wind
        ! difference. In HIRLAM 10m wind can NEVER be stronger than 1-st level 
        ! one. Since it is in denominator, not too small value must be used 
        ! for cut-off. See later the free-convection case
        !
        windspeed_1 = sqrt(u_1(i)*u_1(i)+ v_1(i)*v_1(i))
        if(ifWind10m)then
          du = windspeed_1 - wind10m(i) ! [m/s]
        else
          du  = sqrt((u_1(i)-u_10m(i))**2 + (v_1(i)-v_10m(i))**2)
        endif


        !
        ! Temperature difference. Formally, there should be potential tempr.
        ! difference but we put here just an ordinary one. Extra addition
        ! will hift everything further into stable case a few lines further.
        ! However, thickness of the stable layer should be reduced not immediately
        ! when dQheta > 0 but slightly later - as I do by taking ordinary temperature
        !
        ! Second trick: here I put the wetness correction, from the fields defined
        ! above. If dry_abl_param then this correction will be zero
        !
          ! Workaround for Joana's Ensclim runs: some of  q2m are real_missing
          ! there  !FIXME  Same check below
          !
        if (any((/temp_1(i),temp2m(i), q_1(i), q2m(i)/) == real_missing)) then
           call msg_warning("Crazy input index_meteo="+fu_str(i), "dq_ABL_params")
           call msg("temp_1(i),temp2m(i), q_1(i), q2m(i)", (/temp_1(i),temp2m(i), q_1(i), q2m(i)/))
           dT = 0
        else
        dT = temp_1(i)-temp2m(i) + 0.608*temp2m(i)*(q_1(i)-q2m(i))
           if (.not. abs(dt)>=0.) then 
              call msg_warning("Failed to make dt index_meteo="+fu_str(i), "dq_ABL_params")
              call msg("temp_1(i),temp2m(i), q_1(i), q2m(i)", (/temp_1(i),temp2m(i), q_1(i), q2m(i)/))
              dT = 0
           endif
        endif



        !
        ! Varying value of sigma, plus an update the thickness of the stable 
        ! surface layer. Strictly speaking, there should be gradual thinning
        ! of the layer but so far I put just a straight value. Later a more
        ! sophisticated stuff can come.
        !
        if(ifTemp2m)then
          if(dT > 0.7)then   ! Roughly corresponds to dTheta = 1.0 degree
            height_1T = min(lev1_Theight_stable, height_lev1(i))
          else
            height_1T = height_lev1(i) 
          endif
          height_T2m = 2.0
        else
          if(dT > 0.7)then   ! Roughly corresponds to dTheta = 1.0 degree
            height_1T = min(lev1_Theight_stable, height_lev2(i))
          else
            height_1T = height_lev2(i) 
          endif
          height_T2m = height_lev1(i)
        endif

        if(ifWind10m)then
          height_1u = height_lev1(i) 
          height_u10m = 10.0
        else
          height_1u = height_lev2(i) 
          height_u10m = height_lev1(i)
        endif

        lnzu = LOG(height_1u/height_u10m)
        dzu = height_1u - height_u10m
        lnzT = LOG(height_1T/height_T2m)
        dzT = height_1T - height_T2m

        !
        ! Now switch to potential temperature having height already small
        ! if the case is really stable - in order not to make it too stable.
        ! This all is, of course, a small additions but anyway...
        !
        if(land_mask(i) < 0.5)then ! Less than 50% of land in the grid cell
          gamma = gamma_wet(min(max(int((temp2m(i)-273.15 + 20) / 10 + 1.5),1), size(gamma_wet)))
        else
          gamma = gamma_dry
        endif

        dT = dT + gamma * Kz_ref_height !(height_1T - height_T2m) ! [K] 0.01K/m is dry adibatic correction

        ! 
        ! The real analysis of stability to determine sigma and make actual 
        ! computations of Kz - for the case of free convection and 
        ! a combination of thermal and wind induced turbulence.
        ! Since there may be errors in HIRLAM due to boundary values, 
        ! have to apply below non-elegant way of coding.
        ! However, it ensures that there are no unnecessry computations.
        !
        if(dT > 0.)then
          !
          ! Stable layer: possibly thin layer, stable sigma (large) and at least some
          ! forced wind-induced turbulence (du>=0.001)
          !
          sigma = sigma_stable
          if(du<0.01) du=0.01
          x1(i) = Kz_ref_height * sigma * g * lnzu * lnzu * dT / (temp2m(i)*lnzt*du*du)
          Kz(i) = karmann_c * karmann_c * Kz_ref_height * du * fu_F(x1(i)) / lnzu
        else
          !
          ! Unstable layer: unstable sigma (small) and possibility for free convection
          ! But first apply run-saving patch for dT: reprojection can lead
          ! to arbitrary relations between surface level and 3D temperature. 
          !
          if(dT < -3.0)then
!            call msg_warning('Very strong dT in unstable regime:' + fu_str(dT),'dq_ABL_params')
!            call msg('T2m, T in the grid cell (i1d='+fu_str(i)+')',temp2m(i),temp_1(i))
            dT = -3.0
          endif
          
          sigma = sigma_unstable

          if(windspeed_1<0.01 .or. (du<0.001 .and. windspeed_1<0.1)) then
            !
            ! Free convection: no wind, negative temperature gradient
            ! There will be very large errors, so it is better to use
            ! the free-convection analytical approximation
            !
            Kz(i) = 1.77778 * karmann_c * karmann_c * & 
                  & SQRT(-dT * sigma * g * HHH_K/(temp2m(i)*lnzt))
            du=0.001
          else
            if(du<0.01) du=0.01
            x1(i) = Kz_ref_height * sigma * g * lnzu * lnzu * dT / (temp2m(i)*lnzt*du*du)
            Kz(i) = karmann_c * karmann_c * Kz_ref_height * du * fu_F(x1(i)) / lnzu
          endif
        endif   ! Stable-unstable case
        !
        ! Absolute max: 2 m2/s - taken for Sahara.
        !
        if(Kz(i) > 2.) Kz(i)=2.
        if(Kz(i) < 1.e-3) then
          Kz(i)=1.e-3
        endif

        Kz(i) = Kz(i) * Kz_scaling
        !
        ! If fully log-profile is assumed up to the second level:
        ! friction_velocity(i) = SQRT(Kz(i) * du /(Kz_ref_height*lnzu))  
        !
        ! For log-profile up to the Kz-ref-height (1m)
        !
        friction_velocity(i) = Kz(i) / (karmann_c * Kz_ref_height)
        fric_vel_cube = friction_velocity(i) * friction_velocity(i) * friction_velocity(i)


        ! Heatflux upward is positive (unstable case). Prandtl number depends on stability,
        ! its dependence is just obtained from non-linear regression to Laihtman (unstable) and 
        ! Zilitinkevich (stable) data. The trick is that there have to be Ri number 
        ! (indeed already available) and L (not yet available).
        !
        if(dT > 0.05)then
          !
          ! Stable case, alpha_Pr = F1(Ri) - from Zilitinkevich
          !
          Ri = 2. * g * dzu * dzu * dT / (du * du * dzT * (temp_1(i) + temp2m(i)))
          if(Ri > 10.)Ri = 10.   ! let's limit the Ri value to avoid too low Pr

!          alpha_Pr(i) = alpha_Pr_neutral /(1. + 2.8733 * Ri(i)**1.2653) old stuff 3.6.2
!          alpha_Pr(i) = alpha_Pr_neutral /(1. + 17.76624 * Ri(i)**1.2653) ! fit for limited Ri interval
          alpha_Pr(i) = alpha_Pr_neutral /(1. + 11.1663 * Ri**1.0972) ! fit for wider Ri interval

!if(alpha_Pr(i) < 0.1)then
!  call msg('Large Pr (Ri, Pr):', Ri(i), 1./alpha_Pr(i))
!endif
!if(alpha_Pr(i) > 10.)then
!  call msg('Small Pr (Ri, Pr):', Ri(i), 1./alpha_Pr(i))
!endif

        elseif(dT < -0.05)then
          !
          ! Unstable case, alpha_Pr = F2(z1/L)
          !
          ! A trick: we have alpha_Pr = y/myTmp and non-linear regression arg=y/F(y)
          ! See notebook 7, p.40 for details.
          !
          arg = Kz(i) * karmann_c * g * dT / (fric_vel_cube *  temp2m(i) * lnzt)
          call alpha_Pr_solver(-arg, y_solver, alpha_Pr(i))  ! See the numerical solver below

!if(alpha_Pr(i) < 0.1)then
!  call msg('Large Pr (arg, Pr):', arg, 1./alpha_Pr(i))
!endif
!if(alpha_Pr(i) > 10.)then
!  call msg('Small Pr (arg, Pr):', Ri(i), 1./alpha_Pr(i))
!endif

!          do iTmp = 1,1000
!            call msg('Calling solver with arg=',real(iTmp) * real(iTmp) * 0.00001)
!            call alpha_Pr_solver(real(iTmp) * real(iTmp) * 0.0001, y_solver, alpha_Pr(i))
!            call msg('z/l, alpha_prandtl:',y_solver, alpha_pr(i))
!            call msg('')
!          end do
!          stop

        else
          !
          ! Very small temperature difference, uncertain values. Set neutral 
          !
          alpha_Pr(i) = alpha_Pr_neutral
        endif

        !
        ! Having the alpha_Pr, we can proceed with the other parameters
        !
        sens_heat_flux(i) = - Cp_ro * Kz(i) * alpha_Pr(i) * dT / (Kz_ref_height * lnzt)

        if(.not.((sens_heat_flux(i) < fMaxAlert) .and. sens_heat_flux(i) > fMinAlert)) then
            call msg("Strange sensible heat flux at imeteo="+fu_str(i), sens_heat_flux(i) )
            call msg("Cp_ro,  Kz(i) , alpha_Pr(i) , dT , Kz_ref_height , lnzt):", (/Cp_ro, Kz(i) , alpha_Pr(i) , dT , Kz_ref_height , lnzt/))
            CALL set_error('Strange sensible heatflux','dq_ABL_params')
            return
        endif
        if(.not.((sens_heat_flux(i) < fMaxForce) .and. sens_heat_flux(i) > fMinForce)) then
            call msg("Strange sensible heat flux at imeteo="+fu_str(i), sens_heat_flux(i) )
            call msg("Cp_ro,  Kz(i) , alpha_Pr(i) , dT , Kz_ref_height , lnzt):", (/Cp_ro, Kz(i) , alpha_Pr(i) , dT , Kz_ref_height , lnzt/))
            CALL msg_warning('Strange sensible heatflux, resetting to range: [' + &
                           & fu_str(fMinForce) + ',' + fu_str(fMaxForce) + ']','dq_ABL_params')
            sens_heat_flux(i) = max(min(sens_heat_flux(i),fMinForce),fMaxForce)
        endif

        temperature_scale(i) = -sens_heat_flux(i) / (Cp_ro * friction_velocity(i))

        monin_obukhov_inv(i)= - karmann_c * g * sens_heat_flux(i) / (fric_vel_cube * Cp_ro * temp2m(i))

        if(monin_obukhov_inv(i) > 100)then
          call msg('Large 1/L:',monin_obukhov_inv(i))
        endif

        alpha_Pr(i) = min(1./ alpha_Pr(i), 5.0) ! for the below storage to have real Prandtl nbr, not inverse
                                                ! Also, the upper limit is forced for the very strong stability
                                                ! where the whole methodology is anyway inapplicable
        !
        ! Humidity scale and latent heat flux are computed only in case of full-ABL
        !
        if(abl_param == abl_full_param)then ! Full-blown parameterization, with wetness correction

          vap_heat = ps(i)/(gas_constant_dryair*temp2m(i)) * vaporization_latentheat !Cp*ro
          !
          ! Heatflux upward is positive (unstable case)
          !

          !
          ! Workaround for Joana's Ensclim runs: some of  q2m are real_missing
          ! there  !FIXME
          !
          if (any((/temp_1(i),temp2m(i), q_1(i), q2m(i)/) == real_missing)) then
                  latent_heat_flux(i) = 0.
          else
          latent_heat_flux(i) = - vap_heat  * Kz(i) * alpha_Pr(i) * &
                              & (q_1(i) - q2m(i)) / (Kz_ref_height * lnzt)
          endif
        
          humid_scale(i) = latent_heat_flux(i) / &
                         & (karmann_c * vap_heat * (friction_velocity(i) + 1.0e-5))
        end if ! if dry or wet ABL

      END DO  ! i=1:fs_meteo

      !-------------------------------------------------------------------------------
      !
      ! ABL height is the next step. Note that most methods require full-3D analysis,
      ! which we try to make as fast as possible, each time working with the whole level
      !
      SELECT CASE (ABL_height_method_switch)

        CASE (constant_abl_height)                                           ! constant abl height
          abl_height(1:fs_meteo) = 1000.
          abl_height_m_pointer = abl_height_m_flag

        CASE (parcel_method, richardson_method, combination_method)          ! parcel and/or Richardson

          u3d => fu_sm_obstime_3d_field(meteoMarketPtr, met_src, u_flag, time, single_time)
          v3d => fu_sm_obstime_3d_field(meteoMarketPtr, met_src, v_flag, time, single_time)

          abl_height(1:fs_meteo) = -9.
          if(ABL_height_method_switch == parcel_method .or. &
           & ABL_height_method_switch == combination_method) call apply_parcel(abl_height)
          if(error)return
          if(ABL_height_method_switch == richardson_method) call apply_richardson(abl_height, abl_param == abl_full_param)
          if(error)return
          if(ABL_height_method_switch == combination_method)then
            x1(1:fs_meteo) = -9.
            call apply_richardson(x1, abl_param == abl_full_param)  ! temporary use of the variable
            if(error)return
            DO i = 1, fs_meteo
              IF (x1(i) > abl_height(i)) abl_height(i) = x1(i)
              if(abl_height(i) > 7500.) abl_height(i) = 7500.
            END DO
          endif
          abl_height_m_pointer = abl_height_m_flag

        case (coriolis_method)                                               ! Coriolis-parameter and Kz
          call coriolis_parameters(dispersion_grid, x1)
          IF (error) RETURN 
          do i=1,fs_meteo
            abl_height(i) = Kz(i) / (Kz_ref_height * x1(i))
          end do
          abl_height_m_pointer = abl_height_m_flag

        CASE (nwp_abl)
          abl_height_m_pointer = nwp_abl_height_m_flag

        case default
          call set_error('Unknown ABL height method:' + fu_str(ABL_height_method_switch),'dq_ABL_params')
          return
      end select

      !
      ! Replace ABL height with default if failed
      !
      i=count(abl_height < 0.0)
      if (i>0) then
         call msg_warning("BLH is outside meteo at "//trim(fu_str(i))//"points", "dq_ABL_params")
         call msg("Setting abl height to meteo top....")
         height_lev1 => fu_grid_data_from_3d(height3d, nz_meteo)
         where (abl_height<0) abl_height = height_lev1
      endif

!call msg_warning("Forcing fixed ABL and Kz1m", "dq_ABL_params")
!abl_height(1:fs_meteo) = 1000.
!Kz (1:fs_meteo) = 1e-2
      

      !-------------------------------------------------------------------------------
      !
      !  Convective velocity scale requires abl height and corrects the friction velocity.
      !
      DO i = 1, fs_meteo
        IF(sens_heat_flux(i) > 0.)THEN
          !
          ! Note: w* = (Habl * g * Qsens/(Cp*Rho*T2m))^0.333, where
          ! Cp_ro = specific_heat_dryair*ps(i)/(gas_constant_dryair*temp2m)
          !
          conv_vel_scale(i) = (abl_height(i) * sens_heat_flux(i) * g * gas_constant_dryair / &
                            & (specific_heat_dryair * ps(i))) ** 0.3333333333333333333
          !
          ! Now, we update the friction velocity for convection case, after Beljaars et al, QJRMS 1994:
          ! u10m_effective = sqrt(u10m**2 + (1.2*conv_vel_scale)**2).
          ! Then, u10m = fric_vel/karmann_c * ln(10m/z0) - or close to that. 
          ! For u*, we get:
          ! uStar_eff = sqrt(uStar**2 + (1.2*cnv_vel_scale*karmann_c/ln(10m/z0))**2) - for neutral strat.
          ! For unstable, profile will be different. 
          ! As a crude thing, 0.1 of convective scale was taken. It makes uStar diurnal cycle in desert flat
          !
!          frction_velocity(i) = sqrt(frction_velocity(i)*frction_velocity(i) + &
!                                   & (1.2*conv_vel_scale(i)*karmann_c/log(10./z0(i)))**2)

          friction_velocity(i) = sqrt(friction_velocity(i)**2 + (0.1 * conv_vel_scale(i))**2)

        ELSEIF(sens_heat_flux(i) <= 0.)THEN
          conv_vel_scale(i) = 0.
        ELSE
          if(isNaN(sens_heat_flux(i)))then
            CALL set_error('NaN sensible heatflux','dq_ABL_params')
          else
            CALL set_error('Neither positive nor negative sensible heatflux','dq_ABL_params')
          endif
          call msg('Wrong sensible heatflux, i=', i, sens_heat_flux(i))
          RETURN
        END IF

      END DO   ! 1: fs_meteo


      IF (error) RETURN 

    END DO loop_over_times

    call free_work_array(density)
    call free_work_array(x1)
    if(abl_param == abl_dry_param)then
      call free_work_array(q2m)
      call free_work_array(q_1)
    endif

    CONTAINS

    REAL FUNCTION fu_F(x)
      !
      ! Contains Genikhovich's integral approximated by 3-part algebraic functions
      !
      IMPLICIT NONE
      REAL, INTENT(in) :: x
      integer :: m_stable = 5
      if(x < -1.)then
        fu_F = (1.241+(0.724*((-x-1)**0.9)))*(1.241+(0.724*((-x-1)**0.9)))/(-x)
      elseif(x >= -1. .and. x < 0.)then
        fu_F = 1.+0.54*((-x)**0.8)
      ELSEIF(x >=0. .and. x < 1.)THEN
        fu_F = 1./(1.+0.9*x)
      ELSEIF(x>=1.)THEN
!        fu_F = 0.526316/(x**(m_stable+1))
        fu_F = 0.526316/x
      ELSE
        CALL set_error('Probably NaN-argument','fu_F')
        fu_F=0.
      END IF
    END FUNCTION fu_F

    !===============================================================================

    subroutine alpha_Pr_solver(arg, y, alpha_pr)
      !
      ! Solves the high-order polinomial equation and finds alpha_Pr
      ! Equation is obtained from non-linear fitting of Laihtman's curve:
      ! 0.048x**3 - 0.1756x**2 + 0.4988x
      !
      implicit none
      real, intent(in) :: arg
      real, intent(out) :: y, alpha_pr

      real :: iterator

      ! As seen from the fit, for small my_value, we have a simple linear relation
      !
      if(arg < 0)then

        call msg('Negative argument:',arg)
        call set_error('Negative argument','fu_alpha_Pr_solver')

      elseif(arg < 0.0001)then

        y = 1.3 * arg
        alpha_Pr = fu_alpha_pr(y)  ! for very small value, nothing to solve, nearly neutral

      elseif(arg > 1.0)then

        y = 3.0 * arg
        alpha_Pr = fu_alpha_pr(y)  ! for large value, linear approximation is also OK

      else
        !
        ! Initial point for iterations comes from approximate linear fitting of the same function
        !
        y = 2.0 * arg  ! initial step NOTE 2.0 here is due to alpha_pr_limit = 2
        iterator = y * fu_alpha_pr(y)
        do while(.not. (iterator .eps. arg))
          y = y - (iterator - arg) * 0.5   ! reduce the distance
          iterator = y * fu_alpha_pr(y)  ! still not too bad, although lower than the actual function
        end do
        alpha_Pr = fu_alpha_pr(y)

      endif        
      
    end subroutine alpha_Pr_solver

    !=============================================================================

    real function fu_alpha_pr(y)
      real, intent(in) :: y
      !
      ! The actual fit of Laihtman
      !
!      fu_alpha_pr = 3.0 - 1.65 /(1.0 + 9.28 * y**1.25)  ! the approximation itself
      !
      ! Well, it looks like the max value of 3 for alpha is not the best: both sea 
      ! and land areas demonstrate noticeable (sea area - very strong) over-stating 
      ! the unstable features. Voluntaristically, we put a limit of 2 keeping the neutral 
      ! as 1.35, of course.
      ! ATTENTION. Solver is adapted to this very value
      !
      fu_alpha_pr = 2.0 - 0.65 /(1.0 + 9.28 * y**1.25)
        
    end function fu_alpha_pr

    !==============================================================================
    
    subroutine apply_parcel(abl_height)
      !
      ! Goes along the array and updates all cells, where abl_height is negative
      !
      implicit none
      
      ! Imported parameters
      real, dimension(:), pointer :: abl_height

      ! Local variables
      real, dimension(:), pointer :: ground_pot_temperature, theta_up, theta_down, h2d_up, h2d_down
      integer :: i, iCount, iCount2, iLev
      
      real, parameter :: a = 0.5, b=1.2

      ground_pot_temperature => fu_work_array()

      theta_up  => fu_grid_data_from_3d(theta3d, 1)

      do i = 1, fs_meteo
        if(sens_heat_flux(i) <= 0.0)then
          ground_pot_temperature(i) = temp2m(i)*((std_pressure_sl/ps(i))**R_per_c_dryair) + a
        else
          ground_pot_temperature(i) = temp2m(i)*((std_pressure_sl/ps(i))**R_per_c_dryair) + b
        endif
      end do

      iCount = 0
      iCount2 = 0
      DO iLev = 1, (nz_meteo-1)

        theta_up  => fu_grid_data_from_3d(theta3d, iLev+1)
        theta_down => fu_grid_data_from_3d(theta3d, iLev)
        h2d_up => fu_grid_data_from_3d(height3d, iLev+1)
        h2d_down => fu_grid_data_from_3d(height3d, iLev)
        IF (error) EXIT

        DO i = 1, fs_meteo
          IF (abl_height(i)<= 0.) THEN
            IF (theta_up(i) > ground_pot_temperature(i)) THEN 
              if(theta_down(i) > ground_pot_temperature(i))then
                abl_height(i) = h2d_down(i)
                iCount2 = iCount2 + 1
              else
                abl_height(i) = h2d_down(i) * (theta_up(i) - ground_pot_temperature(i)) &
                            & + h2d_up(i) * (ground_pot_temperature(i) - theta_down(i))
                abl_height(i) = abl_height(i) / (theta_up(i) - theta_down(i))
              endif
              iCount = iCount + 1
            END IF
          END IF
        END DO  !loop 1:fs
        abl_height(1:fs_meteo) = anint(abl_height(1:fs_meteo))+0.2 !! So it can be
                                                       !! distinguished
        IF(iCount >= fs_meteo) EXIT

      END DO ! loop_over_levels

      call free_work_array(ground_pot_temperature)

    end subroutine apply_parcel

    !============================================================================================
    
    subroutine apply_richardson(Habl, ifFullABLParam)
      !
      ! Calculates the ABL height using critical Ri approach
      !
      implicit none

      ! Imported parameters
      real, dimension(:), pointer :: Habl

      ! Local declarations:
      INTEGER :: iLev, i, iCount, iFlip
      REAL, DIMENSION(:), POINTER :: u_1, u_2, v_1, v_2, q_1, q_2, theta_1, theta_2, h2d_1, h2d_2, h2d_2prev
      type(silja_rp_1d), DIMENSION(:), POINTER :: richardson
      logicaL, intent(in) :: ifFullABLParam

      ! Local parameters:
      REAL, PARAMETER :: critical_richardson = 0.25
      !
      ! Get the space and define the flipping pointers
      !
      call get_work_arrays_set(2, fs_meteo, richardson)
      iFlip = 1

      richardson(1)%pp(1:fs_meteo) = 0.0            ! The first-level Ri is always zero
      Habl(1:fs_meteo) = -9.
      !
      !  Loop over levels. Note that we shall stop immediately when ABL is defined in all cells
      !
      u_1 => fu_grid_data_from_3d(u3d, 1)
      v_1 => fu_grid_data_from_3d(v3d, 1)
      theta_1 => fu_grid_data_from_3d(theta3d, 1)
      h2d_1 => fu_grid_data_from_3d(height3d, 1)
      h2d_2prev => h2d_1
      if(ifFullABLParam) q_1 => fu_grid_data_from_3d(q3d, 1)
      iCount = 0
      DO iLev = 2, nz_meteo
        u_2 => fu_grid_data_from_3d(u3d, iLev)
        v_2 => fu_grid_data_from_3d(v3d, iLev)
        theta_2 => fu_grid_data_from_3d(theta3d, iLev)
        h2d_2 => fu_grid_data_from_3d(height3d, iLev)
        if(ifFullABLParam) q_2 => fu_grid_data_from_3d(q3d, iLev)
        IF (error) EXIT
        !
        !  Get bulk_richardson_number for the second flipping array
        !
        call bulk_richardson_nbr_one_level(u_1, u_2, v_1, v_2, q_1, q_2, theta_1, theta_2, &
                                         & h2d_1, h2d_2, friction_velocity, ifFullABLParam, &
                                         & richardson(3-iFlip)%pp)
        IF (error) EXIT
        !
        !  Critical value reached?
        !
        DO i = 1, fs_meteo

          IF (Habl(i) < 0.0) THEN        
            IF (richardson(3-iFlip)%pp(i) > critical_richardson) THEN
              Habl(i) =  h2d_2prev(i) *  (richardson(3-iFlip)%pp(i) - critical_richardson)  &
                       & +   h2d_2(i) *  (critical_richardson - richardson(iFlip)%pp(i))

              Habl(i) = Habl(i) / (richardson(3-iFlip)%pp(i) - richardson(iFlip)%pp(i))
              
              iCount = iCount + 1
            END IF
          END IF
        END DO !loop_over_gridpoints

        IF (iCount >= fs_meteo) EXIT   ! All done?
        !
        ! Flip Ri arrays: the (3-iFlip) array will be used as iFlip in the next cycle.
        !
        h2d_2prev => h2d_2
        iFlip = 3 - iFlip

      END DO ! loop_over_levels
      abl_height(1:fs_meteo) = anint(abl_height(1:fs_meteo))+0.5 !! So it can be
                                                     !! distinguished

      call free_work_array(richardson)

    end subroutine apply_richardson

  END SUBROUTINE dq_ABL_params
  

  !***************************************************************

  SUBROUTINE dq_omega(meteoMarketPtr, met_src, valid_times, ifUpdate)
    !
    ! Calculates vertical velocity omega [Pa/s] using HIRLAM's post prosessing
    ! scheme. 
    !
    ! ONLY for hybrid model level!
    ! Omega is defined at half-levels, i.e. at the interface of the between the
    ! grid boxes.
    !

    IMPLICIT NONE
    
    ! Imported parameters with intent IN:
    type(meteo_data_source), INTENT(in) :: met_src
    TYPE(silja_time), DIMENSION(:),INTENT(in) :: valid_times
    logical, intent(in) :: ifUpdate
    type(mini_market_of_stacks), pointer :: meteoMarketPtr
    
    ! Local declarations:
    TYPE(silja_time) :: time
    TYPE(silja_3d_field), POINTER :: u3d, v3d, t3d !, omega3d
    TYPE(silja_field_id) :: id
    TYPE(silja_grid) :: u_grid, v_grid, w_grid
    TYPE(silja_level), DIMENSION(max_levels) :: levels
    TYPE(silja_field), POINTER :: fieldpointer, fieldpointer_up 
    TYPE(silja_field), POINTER :: ps_2d 
    TYPE(silja_field), POINTER :: ps_past_field, ps_future_field
    REAL :: a_down, b_down, a_up, b_up, ab, db
    INTEGER :: number_of_levels
    REAL, DIMENSION(:), POINTER :: u, v, w  
    INTEGER :: i, lev, tloop
    REAL, DIMENSION(:), POINTER :: pressure, p_up, p_down, dp
    REAL, DIMENSION(:), POINTER :: du_dx, dv_dy
    REAL, DIMENSION(:), POINTER :: ps 
    REAL, DIMENSION(:), POINTER :: dps_dx, dps_dy, dps_dt
    REAL, DIMENSION(:), POINTER :: ps_past, ps_future
    REAL, DIMENSION(:), POINTER :: vertical_omega 
    REAL, DIMENSION(:), POINTER :: div, ps_adv
    REAL, DIMENSION(:), POINTER :: wh, wf, omega_tmp
!    logical :: ifExists
    real :: dt


    p_up => fu_work_array()
    p_down => fu_work_array()
    pressure => fu_work_array()
    dp => fu_work_array()
    div =>fu_work_array()
    ps_adv =>fu_work_array()
    wf => fu_work_array()
    wh => fu_work_array()
    dps_dx => fu_work_array()
    dps_dy => fu_work_array()
    dps_dt => fu_work_array()
    dv_dy => fu_work_array()
    du_dx => fu_work_array()


    IF (error) RETURN

    loop_over_times: DO tloop = 1, SIZE(valid_times)
      IF (.NOT.defined(valid_times(tloop))) EXIT loop_over_times
      time = valid_times(tloop)
      
      IF (fu_field_in_sm(meteoMarketPtr, met_src, & ! Already done before?
                       & omega_flag,&
                       & time,&
                       & level_missing,&
                       & .true.,&
                       & .false.))then
        call msg('... already in SM')
        if(ifUpdate)then
!          omega3d => fu_sm_obstime_3d_field(meteoMarketPtr, met_src, omega_flag, time, single_time)
!          ifExists = .true.
        else
          CYCLE loop_over_times
        endif
!      else
!        ifExists = .false.
!        vertical_omega => fu_work_array()
      endif
      
      t3d => fu_sm_obstime_3d_field(meteoMarketPtr, met_src, temperature_flag, time, single_time)
      IF (error) EXIT loop_over_times
        
      if(fu_leveltype(t3d) /= hybrid)then
        call set_error('Vertical wind is avalilable only for hybrid levels',& 
                     & 'dq_omega')
        exit loop_over_times
      end if
      
      CALL vertical_levels(t3d, levels, number_of_levels)
      IF (error) EXIT loop_over_times

      u3d => fu_sm_obstime_3d_field(meteoMarketPtr, met_src, u_flag, time, single_time)
      v3d => fu_sm_obstime_3d_field(meteoMarketPtr, met_src, v_flag, time, single_time)
      IF (error) EXIT loop_over_times
          
      u_grid = fu_grid(u3d)
      v_grid = fu_grid(v3d)


      ! Surface pressure tendency
      ! Past field
      if(tloop == 1)then
        ps_past_field => fu_sm_obstime_field(meteoMarketPtr, met_src, ground_pressure_flag, &
                                           & level_missing, time, backwards)
      else
        ps_past_field => fu_sm_obstime_field(meteoMarketPtr, met_src, ground_pressure_flag, &
                                           & level_missing, time-ten_seconds, backwards)
        if(error)then
          call unset_error('dq_omega')
          ps_past_field => fu_sm_obstime_field(meteoMarketPtr, met_src, ground_pressure_flag, &
                                             & level_missing, time, single_time)
        endif
      endif
      ps_past  => fu_grid_data(ps_past_field)   

      ! Future field
      if(fu_valid_time(ps_past_field) == time)then 
        ps_future_field => fu_sm_obstime_field(meteoMarketPtr, met_src, ground_pressure_flag, &
                                             & level_missing, time+ten_seconds,  forwards)
      else
        ps_future_field => fu_sm_obstime_field(meteoMarketPtr, met_src, ground_pressure_flag, &
                                             & level_missing, time,  single_time)
      endif

      if(error)then
        call set_error('Only one time in supermarket, cannot compute pressure tendency', 'dq_omega')
        EXIT loop_over_times
      endif

      ps_future  => fu_grid_data(ps_future_field) 
      dt = fu_sec(fu_valid_time(ps_future_field)-fu_valid_time(ps_past_field))
      IF (error) EXIT loop_over_times
     
      do i=1,fs_meteo
       dps_dt(i) = (ps_future(i)-ps_past(i))/dt
      enddo

      ! Surface pressure
      ps_2d => fu_sm_obstime_field(meteoMarketPtr, met_src, ground_pressure_flag, level_missing, &
                                 & time, single_time)
      IF (error) EXIT loop_over_times
        
      ps => fu_grid_data(ps_2d)
      IF (error) EXIT loop_over_times
      
      !
      ! Preparation to surface-pressure advection
      ! dps_dx:
      CALL ddx_of_field(ps, meteo_grid, meteo_grid, dps_dx)
      ! dps_dy:
      CALL ddy_of_field(ps, meteo_grid, meteo_grid, dps_dy)
      IF (error) EXIT loop_over_times

      wh(1:fs_meteo) = 0.

      loop_over_levels: DO lev = 1, number_of_levels          !number_of_levels, 1, -1
        
        !
        ! 1. wind data
        !
        u => fu_grid_data_from_3d(u3d, lev)
        v => fu_grid_data_from_3d(v3d, lev)

        id = fu_set_field_id(fu_met_src(t3d),&
                           & omega_flag,&
                           & fu_analysis_time(t3d),&
                           & fu_forecast_length(t3d), &
                           & meteo_grid, &
                           & levels(lev))
        call find_field_data_storage_2d(meteoMarketPtr, id, multi_time_stack_flag, vertical_omega)
        if(fu_fails(.not.error,'Failed vertical_omega field data pointer','dq_omega'))return

!        if(ifExists)then
!          vertical_omega => fu_grid_data_from_3d(omega3d, lev)
!        endif
        IF (error) EXIT loop_over_levels

        !
        ! 2. hybrid coefficients
        !
        IF (lev==number_of_levels) THEN

          fieldpointer => fu_field_from_3d_field(u3d, lev-1) 
          fieldpointer_up => fu_field_from_3d_field(u3d, lev)      
          IF (error) EXIT loop_over_levels
          
          a_down = 0.5*(  fu_hybrid_level_coeff_a(fu_level(fieldpointer_up)) &
                      & + fu_hybrid_level_coeff_a(fu_level(fieldpointer)))
          b_down = 0.5*(  fu_hybrid_level_coeff_b(fu_level(fieldpointer_up)) &
                      & + fu_hybrid_level_coeff_b(fu_level(fieldpointer)))
    
          a_up = fu_hybrid_level_coeff_a(fu_level(fieldpointer_up))
          b_up = fu_hybrid_level_coeff_b(fu_level(fieldpointer_up))

        elseif (lev==1) THEN

          fieldpointer => fu_field_from_3d_field(u3d, lev) 
          fieldpointer_up => fu_field_from_3d_field(u3d, lev+1)
          IF (error) EXIT loop_over_levels
         
          a_down = 0.
          b_down = 1. 

          a_up = 0.5*(  fu_hybrid_level_coeff_a(fu_level(fieldpointer_up)) &
                    & + fu_hybrid_level_coeff_a(fu_level(fieldpointer)))
          b_up = 0.5*(  fu_hybrid_level_coeff_b(fu_level(fieldpointer_up)) &
                    & + fu_hybrid_level_coeff_b(fu_level(fieldpointer)))

        ELSE    
         
          fieldpointer => fu_field_from_3d_field(u3d, lev) 
          fieldpointer_up => fu_field_from_3d_field(u3d, lev+1) 
          IF (error) EXIT loop_over_levels          
       
          a_up = 0.5*(  fu_hybrid_level_coeff_a(fu_level(fieldpointer_up)) &
                    & + fu_hybrid_level_coeff_a(fu_level(fieldpointer)))
          b_up = 0.5*(  fu_hybrid_level_coeff_b(fu_level(fieldpointer_up)) &
                    & + fu_hybrid_level_coeff_b(fu_level(fieldpointer)))

          fieldpointer => fu_field_from_3d_field(u3d, lev-1) 
          fieldpointer_up => fu_field_from_3d_field(u3d, lev)      
          IF (error) EXIT loop_over_levels

          a_down = 0.5*(  fu_hybrid_level_coeff_a(fu_level(fieldpointer_up)) &
                      & + fu_hybrid_level_coeff_a(fu_level(fieldpointer)))
          b_down = 0.5*(  fu_hybrid_level_coeff_b(fu_level(fieldpointer_up)) &
                      & + fu_hybrid_level_coeff_b(fu_level(fieldpointer)))
   
        END IF
        IF (error) EXIT loop_over_levels

        ab = a_down*b_up-a_up*b_down
        db = b_down - b_up   
      
        !
        ! 3. surface pressure advection
        !
        CALL pressure_on_level(t3d, lev, pressure)
        IF (error) EXIT loop_over_levels

        do i=1,fs_meteo
          p_up(i) = a_up + b_up*ps(i)
          p_down(i) = a_down + b_down*ps(i)
          dp(i) = p_down(i) - p_up(i)
          ps_adv(i) = u(i)*dps_dx(i)+v(i)*dps_dy(i)
        end do

        !
        ! 4. mass flux divergence
        !
        ! du_dx:
        CALL ddx_of_field(u, u_grid, meteo_grid, du_dx)
        ! dv_dy:
        CALL ddy_of_field(v, v_grid, meteo_grid, dv_dy)
        IF (error) EXIT loop_over_levels

        do i=1,fs_meteo
          div(i)= (du_dx(i) + dv_dy(i))*dp(i)
        end do

        !
        ! 5. omega, defined at half-levels, i.e. at the top of each box
        !
        do i=1,fs_meteo
          vertical_omega(i) = wh(i)         ! value on lower boundary
          wh(i) = wh(i) - div(i) - ps_adv(i)*db      ! integratable value on upper boundary
         ! Hirlam full level value actually computed on upper boundary of grid-cell
          wf(i) = pressure(i)*ps_adv(i)*(db + ab*(LOG(p_down(i)/p_up(i)))/dp(i))/dp(i)
!          wf(i) = b*ps_adv(i)    !inclination & change of pressure on level

!          vertical_omega(i) = 0.5 * (vertical_omega(i) + wh(i)) + wf(i) + b*dps_dt(i)
!          vertical_omega(i) = 0.5 * (vertical_omega(i) + wf(i) + wh(i)) + dps_dt(i) ! + b*dps_dt(i)
!          vertical_omega(i) = wf(i) + 0.5*(wh_down(i) + wh_up(i)) + b*dps_dt(i)
!          wh_down(i) = wh_up(i)

          vertical_omega(i) = wf(i) + dps_dt(i) + wh(i) + &
                          & (vertical_omega(i)-wh(i))*(log(pressure(i)/p_up(i))/log(p_down(i)/p_up(i))) 

        end do
        !
        ! 5. store the field
        !

!        CALL dq_store_2d(meteoMarketPtr, id, vertical_omega, multi_time_stack_flag )
        IF (error) EXIT loop_over_levels
     
      END DO loop_over_levels

      if(error)EXIT loop_over_times  ! if error, vertical_omega will be released later
      
!      if(.not. ifExists) call free_work_array(vertical_omega)

    END DO loop_over_times

    call free_work_array(dv_dy)
    call free_work_array(du_dx)
    call free_work_array(p_up)
    call free_work_array(p_down)
    call free_work_array(pressure)
    call free_work_array(dp)
    call free_work_array(div)
    call free_work_array(ps_adv)
    call free_work_array(wf)
    call free_work_array(wh)
    call free_work_array(dps_dx)
    call free_work_array(dps_dy)
    call free_work_array(dps_dt)
!    if(error .and. .not. ifExists) call free_work_array(vertical_omega)

    call msg('Done dq_omega')

  END SUBROUTINE dq_omega


  !****************************************************************

  SUBROUTINE dq_vertical_velocity_msl(meteoMarketPtr, met_src, valid_times)
    !
    ! Calculates vertical velocity w (m/s) from omega (Pa/s) 
    ! using adiabatic approximation w = -omega/(g* density)
    ! and hydrostatic approximation p/density = R* temperature 
    !
    ! All units: SI
    ! 
    IMPLICIT NONE
    
    ! Imported parameters with intent IN:
    type(meteo_data_source), INTENT(in) :: met_src
    TYPE(silja_time), DIMENSION(:), INTENT(in) :: valid_times
    type(mini_market_of_stacks), pointer :: meteoMarketPtr

    ! Local declarations:
    TYPE(silja_time) :: time
    INTEGER :: lev, tloop, i, number_of_levels
    TYPE(silja_3d_field), POINTER :: t3d, tPast3d, tFuture3d, omega3d
    TYPE(silja_field_id) :: id
    TYPE(silja_level), DIMENSION(max_levels) :: levels
    REAL, DIMENSION(:), POINTER :: t, omega, pPast, pFuture
    REAL, DIMENSION(:), POINTER :: p, vertical_velocity
    type(silja_time) :: timeFrw, timeBck
    real :: dt_1

    pPast => fu_work_array()
    pFuture => fu_work_array()
!    vertical_velocity => fu_work_array()
    IF (error) RETURN

    loop_over_times: DO tloop = 1, SIZE(valid_times)
      IF (.NOT.defined(valid_times(tloop))) EXIT loop_over_times
      time = valid_times(tloop)

      IF (fu_field_in_sm(meteoMarketPtr, met_src, & ! Already done before?
                       & w_alt_msl_flag,&
                       & time,&
                       & level_missing,&
                       & .true.,&
                       & .false.)) cycle loop_over_times
      !
      ! Get the basic fields
      !
      t3d => fu_sm_obstime_3d_field(meteoMarketPtr, met_src, temperature_flag, time, single_time)
      if(error)return
      call vertical_levels(t3d, levels, number_of_levels)
      if(error)return

      omega3d => fu_sm_obstime_3d_field(meteoMarketPtr, met_src, omega_flag, time, single_time)
      if(error)return

      !
      ! Get both times - forward and backward and time step between them
      !
      i = fu_closest_time(time+one_second, valid_times, forwards, .true.)
      if(i > 0)then
        timeFrw = valid_times(i)
        timeBck = time
      else
        i = fu_closest_time(time-one_second, valid_times, backwards, .true.)
        if(i > 0)then
          timeFrw = time
          timeBck = valid_times(i)
        else
          call set_error('Failed to find suitable second time in any direction', &
                       & 'dq_vertical_velocity_msl')
          return
        endif
      endif
      if(error)return
      dt_1 = 1. / fu_sec(timeFrw - timeBck)
      !
      ! Temperature is needed only for pressure at its levels
      !
      tPast3d => fu_sm_obstime_3d_field(meteoMarketPtr, met_src, temperature_flag, timeBck, single_time)
      tFuture3d => fu_sm_obstime_3d_field(meteoMarketPtr, met_src, temperature_flag, timeFrw, single_time)
      if(error)return

      DO lev = 1, number_of_levels

        t => fu_grid_data_from_3d(t3d, lev)
        omega => fu_grid_data_from_3d(omega3d, lev)
        if(error)return
        CALL pressure_on_level(tPast3d,lev,pPast)
        CALL pressure_on_level(tFuture3d,lev,pFuture)
        if(error)return

        id = fu_set_field_id(fu_met_src(t3d),&
                           & w_alt_msl_flag,&
                           & fu_analysis_time(t3d),&
                           & fu_forecast_length(t3d), &
                           & fu_grid(t3d),&
                           & levels(lev))
        call find_field_data_storage_2d(meteoMarketPtr, id, multi_time_stack_flag, vertical_velocity)
        if(fu_fails(.not.error,'Failed vertical_velocity field data pointer','dq_vertical_velocity_msl'))return

        do i=1,fs_meteo
          vertical_velocity(i) = ((pFuture(i)-pPast(i))*dt_1 - omega(i)) * &
                                & t(i) * gas_constant_dryair * 2./ (g * (pFuture(i)+pPast(i)))
        end do

        if(error)return

!        CALL dq_store_2d(meteoMarketPtr, id, vertical_velocity, multi_time_stack_flag )
        if(error)return

      END DO
    END DO loop_over_times

!    CALL free_all_work_arrays()
    call free_work_array(pPast)
    call free_work_array(pFuture)
!    call free_work_array(vertical_velocity)
   
  END SUBROUTINE dq_vertical_velocity_msl


  !****************************************************************

  SUBROUTINE dq_vertical_velocity_srf(meteoMarketPtr, met_src, valid_times)
    !
    ! Calculates vertical velocity w (m/s) from omega (Pa/s) 
    ! using adiabatic approximation w = -omega/(g* density)
    ! and hydrostatic approximation p/density = R* temperature.
    ! An additional term to be taken into account is the geopotential 
    ! change along the wind vector. Actually, this is a version of transition
    ! from pressure vertical to z-system
    !
    ! All units: SI
    ! 
    IMPLICIT NONE
    
    ! Imported parameters with intent IN:
    type(meteo_data_source), INTENT(in) :: met_src
    TYPE(silja_time), DIMENSION(:), INTENT(in) :: valid_times
    type(mini_market_of_stacks), pointer :: meteoMarketPtr

    ! Local declarations:
    TYPE(silja_time) :: time
    INTEGER :: lev, tloop, i, number_of_levels
    TYPE(silja_3d_field), POINTER :: t3d, omega3d, u3d, v3d, tPast3d, tFuture3d
    TYPE(silja_field_id) :: id
    TYPE(silja_field), pointer :: ps_2d
    TYPE(silja_level), DIMENSION(max_levels) :: levels
    REAL, DIMENSION(:), POINTER :: t, u, v, omega, pPast, pFuture, pressure
    REAL, DIMENSION(:), POINTER :: p, ps, vertical_velocity, dps_dx, dps_dy
    type(silja_time) :: timeFrw, timeBck
    real :: dt_1, ps_adv, fTmp

!    vertical_velocity => fu_work_array()
    dps_dx => fu_work_array()
    dps_dy => fu_work_array()
    pPast => fu_work_array()
    pFuture => fu_work_array()
    IF (error) RETURN

    !
    ! Grand loop over time
    !
    DO tloop = 1, SIZE(valid_times)
      IF (.NOT.defined(valid_times(tloop))) EXIT  ! loop_over_times
      time = valid_times(tloop)

      IF (fu_field_in_sm(meteoMarketPtr, met_src, & ! Already done before?
                       & w_height_srf_flag,&
                       & time,&
                       & level_missing,&
                       & .true.,&
                       & .false.)) cycle  ! loop_over_times
      t3d => fu_sm_obstime_3d_field(meteoMarketPtr, met_src, temperature_flag, time, single_time)
      omega3d => fu_sm_obstime_3d_field(meteoMarketPtr, met_src, omega_flag, time, single_time)
      if(error)return
      call vertical_levels(t3d, levels, number_of_levels)
      if(error)return
      !
      ! In order to get the pressure trend, have to find out two time moments
      ! and corresponding pressure fields
      !
      !
      ! Get both times - forward and backward and time step between them
      !
      i = fu_closest_time(time+one_second, valid_times, forwards, .true.)
      if(i > 0)then
        timeFrw = valid_times(i)
        timeBck = time
      else
        i = fu_closest_time(time-one_second, valid_times, backwards, .true.)
        if(i > 0)then
          timeFrw = time
          timeBck = valid_times(i)
        else
          call set_error('Failed to find suitable second time in any direction', &
                       & 'dq_vertical_velocity_msl')
          return
        endif
      endif
      if(error)return
      dt_1 = 1. / fu_sec(timeFrw - timeBck)

      ! Temperature is needed only for pressure at its levels
      !
      tPast3d => fu_sm_obstime_3d_field(meteoMarketPtr, met_src, temperature_flag, timeBck, single_time)
      tFuture3d => fu_sm_obstime_3d_field(meteoMarketPtr, met_src, temperature_flag, timeFrw, single_time)
      if(error)return

      !
      ! Preparation to surface-pressure advection
      !
      u3d => fu_sm_obstime_3d_field(meteoMarketPtr, met_src, u_flag, time, single_time)
      v3d => fu_sm_obstime_3d_field(meteoMarketPtr, met_src, v_flag, time, single_time)
      IF (error) RETURN

      ps_2d => fu_sm_obstime_field(meteoMarketPtr, met_src, ground_pressure_flag, level_missing, &
                                 & time, single_time)
      IF (error) RETURN
        
      ps => fu_grid_data(ps_2d)
      IF (error) RETURN

      ! dps_dx:
      CALL ddx_of_field(ps, meteo_grid, meteo_grid, dps_dx)
      IF (error) RETURN

      ! dps_dy:
      CALL ddy_of_field(ps, meteo_grid, meteo_grid, dps_dy)
      IF (error) RETURN

!      call msg_warning('Forcing the global-zero vertical wind')

      !
      ! Loop over levels
      !
      DO lev = 1, number_of_levels

        t => fu_grid_data_from_3d(t3d, lev)
        u => fu_grid_data_from_3d(u3d, lev)
        v => fu_grid_data_from_3d(v3d, lev)
        omega => fu_grid_data_from_3d(omega3d, lev)
        CALL pressure_on_level(tPast3d,lev,pPast)
        CALL pressure_on_level(tFuture3d,lev,pFuture)
        if(error)return

        id = fu_set_field_id(fu_met_src(t3d),&
                           & w_height_srf_flag,&
                           & fu_analysis_time(t3d),&
                           & fu_forecast_length(t3d), &
                           & fu_grid(t3d),&
                           & levels(lev))
        call find_field_data_storage_2d(meteoMarketPtr, id, multi_time_stack_flag, vertical_velocity)
        if(fu_fails(.not.error,'Failed vertical_velocity field data pointer','dq_vertical_velocity_srf'))return

        do i=1,fs_meteo
          ps_adv = u(i)*dps_dx(i)+v(i)*dps_dy(i)
          vertical_velocity(i) = ((pFuture(i)-pPast(i))*dt_1 - omega(i) + &
                                & ps_adv*(pFuture(i)+pPast(i))*0.5/ps(i)) * &
                               & t(i) * gas_constant_dryair * 2. / (g * (pFuture(i)+pPast(i)))

        end do

!        fTmp = sum(vertical_velocity(1:fs_meteo)) / real(fs_meteo)
!        call msg('Level and prior-to-forcing mean wind:',lev,fTmp)
!       do i=1,fs_meteo
!          vertical_velocity(i) = vertical_velocity(i) - fTmp
!        end do
!        call msg('Level and mean wind:',lev,sum(vertical_velocity(1:fs_meteo))/real(fs_meteo))


        if(error)return

!        CALL dq_store_2d(meteoMarketPtr, id, vertical_velocity, multi_time_stack_flag )
        if(error)return

      END DO ! levels
    END DO ! loop_over_times

!    CALL free_work_array(vertical_velocity)
    CALL free_work_array(dps_dx)
    CALL free_work_array(dps_dy)
    CALL free_work_array(pPast)
    CALL free_work_array(pFuture)
   
  END SUBROUTINE dq_vertical_velocity_srf


  !****************************************************************

  subroutine dq_wind_divergence(meteoMarketPtr, met_src, valid_times)
    !
    ! Computes the divergence of the wind field. Idea: the transport equation
    ! is solved under the conditions of non-compressible air, which leads to 
    ! continuity equation, which is essentially div(wind)=0. Here we explicitly 
    ! compute the div(wind) for meteo_grid.
    ! We base the computations on the metrical wind in m/s, which makes the 
    ! computations very simple but not necessarily correct.
    !
    implicit none

    ! Imported parameters with intent(in):
    type(meteo_data_source), INTENT(in) :: met_src
    TYPE(silja_time), DIMENSION(:), INTENT(in) :: valid_times
    type(mini_market_of_stacks), pointer :: meteoMarketPtr

    ! Local variables
    integer :: tloop, iLev, i
    real, dimension(:), pointer :: du_dx, dv_dy, dw_dz, div, u, v, w, h
    TYPE(silja_3d_field), POINTER :: u3d, v3d, w3d, height3d
    type(silja_time) :: time
    type(silja_grid) :: u_grid, v_grid, w_grid
    type(silja_field_id) :: id

    !
    ! div(wind) = du/dx + dv/dy + dw/dz
    !
!    div => fu_work_array()
    du_dx => fu_work_array()
    dv_dy => fu_work_array()
    dw_dz => fu_work_array()
    if(error)return

    !
    ! Grand loop over time
    !
    do tloop = 1, SIZE(valid_times)
      IF (.NOT.defined(valid_times(tloop))) EXIT  ! loop_over_times
      time = valid_times(tloop)

      IF (fu_field_in_sm(meteoMarketPtr, met_src, & ! Already done before?
                       & wind_divergence_flag,&
                       & time,&
                       & level_missing,&
                       & .true.,&
                       & .false.)) cycle  ! loop_over_times
      !
      ! Preparations:
      !
      u3d => fu_sm_obstime_3d_field(meteoMarketPtr, met_src, u_flag, time, single_time)
      v3d => fu_sm_obstime_3d_field(meteoMarketPtr, met_src, v_flag, time, single_time)
      w3d => fu_sm_obstime_3d_field(meteoMarketPtr, met_src, vertical_velocity_pointer, time, single_time)
      height3d => fu_sm_obstime_3d_field(meteoMarketPtr, met_src, height_flag, time, single_time)
      if(error)return

      u_grid = fu_grid(u3d)
      v_grid = fu_grid(v3d)
      w_grid = fu_grid(w3d)

      do iLev = 1, nz_meteo
        !
        ! Get the data
        !
        u => fu_grid_data_from_3d(u3d, iLev)
        v => fu_grid_data_from_3d(v3d, iLev)
        !
        ! du_dx:
        !
        CALL ddx_of_field(u, u_grid, meteo_grid, du_dx)
        IF (error) RETURN
        !
        ! dv_dy:
        !
        CALL ddy_of_field(v, v_grid, meteo_grid, dv_dy)
        IF (error) RETURN
        !
        ! dw_dz:
        !
        CALL ddz_of_field_3d(w3d, height3d, meteo_grid, iLev, dw_dz)
        IF (error) RETURN
        
        !
        ! Divergence itself:
        !
        id = fu_set_field_id(fu_met_src(w3d),&
                           & wind_divergence_flag,&
                           & fu_analysis_time(w3d),&
                           & fu_forecast_length(w3d), &
                           & meteo_grid,&
                           & fu_level(meteo_vertical, iLev))
        call find_field_data_storage_2d(meteoMarketPtr, id, multi_time_stack_flag, div)
        if(fu_fails(.not.error,'Failed div field data pointer','dq_wind_divergence'))return

        do i=1,fs_meteo
          div(i)= (du_dx(i) + dv_dy(i) + dw_dz(i))
        end do

        ! w_grid = scalar_grid
        !

!        call msg('9')

!        CALL dq_store_2d(meteoMarketPtr, id, div, multi_time_stack_flag )
        IF (error) RETURN

      end do  ! loop over levels

    end do  ! time loop

!    call free_work_array(div)
    call free_work_array(du_dx)
    call free_work_array(dv_dy)
    call free_work_array(dw_dz)

  end subroutine dq_wind_divergence


  !****************************************************************

  SUBROUTINE dq_abl_top_pressure(meteoMarketPtr, met_src, valid_times)
    !
    ! This function estimates the abl height in Pa having it in m. A simple gas state equation
    ! is used to compute density
    !
    IMPLICIT NONE
    !
    ! Imported parameters with intent(in):
    type(meteo_data_source), INTENT(in) :: met_src
    TYPE(silja_time), DIMENSION(:), INTENT(in) :: valid_times
    type(mini_market_of_stacks), pointer :: meteoMarketPtr

    ! Local declarations
    REAL, DIMENSION(:), POINTER :: abl_pr, tempr, ps, ABLH
    INTEGER :: tloop, i
    TYPE(silja_3d_field), POINTER :: t3d
    type(silja_field), pointer :: ps_2d, ABL_m_2d
    TYPE(silja_time) :: time
    TYPE(silja_field_id) :: id

!    abl_pr => fu_work_array()

!    call msg('ABL-Pa 1')

    if(error)return

    !
    ! Loop over times
    !
    loop_over_times: DO tloop = 1, SIZE(valid_times)
      IF (.NOT.defined(valid_times(tloop))) EXIT loop_over_times
      time = valid_times(tloop)

      IF (fu_field_in_sm(meteoMarketPtr, met_src, & ! Already done before?
                       & abl_top_pressure_flag,&
                       & time,&
                       & level_missing,&
                       & .false.,&  ! 3D
                       & .false.)) CYCLE loop_over_times

      t3d => fu_sm_obstime_3d_field(meteoMarketPtr, met_src,&         
                                  & temperature_flag, &  
                                  & time,&
                                  & single_time)
      if(error)cycle loop_over_times
      
      tempr => fu_grid_data_from_3d(t3d,1)  ! Lowest model level - surely inside ABL
      if(error)cycle loop_over_times

      ps_2d => fu_sm_obstime_field(meteoMarketPtr, met_src,&
                                 & ground_pressure_flag,&
                                 & level_missing,&
                                 & time,&
                                 & single_time)
      if(error)cycle loop_over_times

      ps => fu_grid_data(ps_2d)
      if(error)cycle loop_over_times

      ABL_m_2d => fu_sm_obstime_field(meteoMarketPtr, met_src,&
                                    & abl_height_m_pointer,&
                                    & level_missing,&
                                    & time,&
                                    & single_time)
      if(error)cycle loop_over_times

      ABLH => fu_grid_data(ABL_m_2d)
      if(error)cycle loop_over_times

      ! State equation: ro=p/(RT), dp/dz=ro*g. Then assume constant temperature equal to that 
      ! at the first level whatever it is. Not the best way but accuracy evidently better than that
      ! of the ABL itself
      !
      id = fu_set_field_id(fu_met_src(t3d),&
                         & abl_top_pressure_flag,&
                         & fu_analysis_time(t3d),&
                         & fu_forecast_length(t3d), &
                         & fu_grid(t3d),&
                         & ground_level)
      call find_field_data_storage_2d(meteoMarketPtr, id, multi_time_stack_flag, abl_pr)
      if(fu_fails(.not.error,'Failed abl_pr field data pointer','dq_abl_top_pressure'))return

      do i = 1, fs_meteo
        abl_pr(i) = ps(i) * (1. - g * ablh(i) / (gas_constant_dryair * tempr(i)))
      end do


!      CALL dq_store_2d(meteoMarketPtr, id, abl_pr, multi_time_stack_flag )

    END DO loop_over_times

!    CALL free_work_array(abl_pr)

!    call msg('ABL-Pa 14')

  END SUBROUTINE dq_abl_top_pressure


!  !******************************************************
!  
!  SUBROUTINE dq_albedo(meteoMarketPtr, met_src, valid_times)
!    ! 
!    ! Calculates albedo from instantaneous downward and net sw radiation fluxes
!    !
!    ! All units: SI
!    !
!    ! Language: ANSI Fortran 90
!    !
!    ! Current code owner: Mikhail Sofiev, FMI
!    ! 
!    IMPLICIT NONE
!
!    ! Imported parameters with intent IN:
!    type(meteo_data_source), INTENT(in) :: met_src
!    TYPE(silja_time), DIMENSION(:), INTENT(in) :: valid_times
!    type(mini_market_of_stacks), pointer :: meteoMarketPtr
!
!    ! Local declarations:
!    TYPE(silja_time) :: time
!    INTEGER :: i, tloop
!    TYPE(silja_field), POINTER :: sw_r_down, sw_r_net
!    TYPE(silja_field_id) :: id
!    REAL, DIMENSION(:), POINTER :: srd, srn, alb
!    INTEGER :: fs
!
!!    alb => fu_work_array()
!    IF (error) RETURN
!
!    loop_over_times: DO tloop = 1, SIZE(valid_times)
!      IF (.NOT.defined(valid_times(tloop))) EXIT loop_over_times
!      time = valid_times(tloop)
!
!      IF (fu_field_in_sm(meteoMarketPtr, met_src, & ! Already done before?
!                       & albedo_flag,&
!                       & time,&
!                       & level_missing,&
!                       & .false.,&   !3D
!                       & .false.)) CYCLE loop_over_times
!
!      sw_r_down => fu_sm_obstime_field(meteoMarketPtr, met_src, &
!                                     & surf_sw_down_radiation_flag, &  
!                                     & level_missing, &
!                                     & time,&
!                                     & single_time)
!      IF (error) CYCLE loop_over_times
!
!      sw_r_net => fu_sm_obstime_field(meteoMarketPtr, met_src,&         
!                                    & surf_sw_net_radiation_flag, &  
!                                    & level_missing, &
!                                    & time,&
!                                    & single_time)
!      IF (error) CYCLE loop_over_times
!
!      srd => fu_grid_data(sw_r_down)
!      srn => fu_grid_data(sw_r_net)
!
!      fs = SIZE(srd)
!
!      id = fu_set_field_id(fu_met_src(sw_r_down),&
!                         & albedo_flag,&
!                         & fu_analysis_time(sw_r_down),&
!                         & fu_forecast_length(sw_r_down), &
!                         & fu_grid(sw_r_down),&
!                         & ground_level)
!      call find_field_data_storage_2d(meteoMarketPtr, id, multi_time_stack_flag, alb)
!      if(fu_fails(.not.error,'Failed abla field data pointer','dq_albedo'))return
!
!      alb(1:fs) = MIN(MAX(srd - srn, 0.) / MAX(srd,0.00001), 1.000)
!
!
!!      CALL dq_store_2d(meteoMarketPtr, id, alb, multi_time_stack_flag )
!
!    END DO loop_over_times
!
!!    CALL free_work_array(alb)
!
!  END SUBROUTINE dq_albedo

  !****************************************************

    
  SUBROUTINE dq_rad_sw_down_sfc(meteoMarketPtr, met_src, valid_times)
    ! 
    ! Used if net radiation given, but we need downwards component (e.g. some HIRLAM data)
    ! Climatological albedo flag here can mean actually the dynamic meteo albedo 
    ! Cumulative fluxes are used because only one derivation path is allowed per quantity 
    ! and in most cases instant flux is computed from cumulative itself. 
    ! !!! This derivation is actually wrong !!!
    ! Albedo being an instant changing quantity in the files can lead to all kinds of 
    ! "interesting" results when the accumulation period for fluxes is more than one meteo step. 
    ! Hopefully still OK for most runs ..
    ! 
    IMPLICIT NONE

    ! Imported parameters with intent IN:
    type(meteo_data_source), INTENT(in) :: met_src
    TYPE(silja_time), DIMENSION(:), INTENT(in) :: valid_times
    type(mini_market_of_stacks), pointer :: meteoMarketPtr

    ! Local declarations:
    TYPE(silja_time) :: time
    INTEGER :: i, tloop
    TYPE(silja_field), POINTER :: sw_r_net, albedo
    TYPE(silja_field_id) :: id
    REAL, DIMENSION(:), POINTER :: srd, srn, alb
    INTEGER :: fs

    loop_over_times: DO tloop = 1, SIZE(valid_times)
      IF (.NOT.defined(valid_times(tloop))) EXIT loop_over_times
      time = valid_times(tloop)

      IF (fu_field_in_sm(meteoMarketPtr, met_src, & ! Already done before?
                       & surf_sw_down_radiation_ac_flag,&
                       & time,&
                       & level_missing,&
                       & .false.,&   !3D
                       & .false.)) CYCLE loop_over_times

      albedo => fu_sm_obstime_field(meteoMarketPtr, met_src, &
                                     & climatological_albedo_flag, &  
                                     & level_missing, &
                                     & time,&
                                     & single_time)
      IF (error) CYCLE loop_over_times

      sw_r_net => fu_sm_obstime_field(meteoMarketPtr, met_src,&         
                                    & surf_sw_net_radiation_ac_flag, &  
                                    & level_missing, &
                                    & time,&
                                    & single_time)
      IF (error) CYCLE loop_over_times

      alb => fu_grid_data(albedo)
      srn => fu_grid_data(sw_r_net)

      fs = SIZE(srn)
      
      id = fu_set_field_id(fu_met_src(sw_r_net),&
                         & surf_sw_down_radiation_ac_flag,&
                         & fu_analysis_time(sw_r_net),&
                         & fu_forecast_length(sw_r_net), &
                         & fu_grid(sw_r_net),&
                         & ground_level,&
                         & length_of_accumulation=fu_accumulation_length(sw_r_net), & 
                         & length_of_validity=fu_validity_length(sw_r_net), &     
                         & field_kind=fu_field_kind(fu_id(sw_r_net)))             

      call find_field_data_storage_2d(meteoMarketPtr, id, multi_time_stack_flag, srd)
      if(fu_fails(.not.error,'Failed field data pointer','dq_rad_sw_down_sfc'))return

      do i = 1, fs
        if(alb(i)<1.0)then
          srd(i) = srn(i) / (1.0-alb(i))
        else ! actually no idea what it was, net should be 0
          srd(i) = srn(i)
          call msg('Mirror earth: ', alb(i), srd(i))
        endif
        if(srd(i)<0.0)then
          call msg('Shiny shiny earth: ', alb(i), srd(i))
        endif
      enddo

    END DO loop_over_times


  END SUBROUTINE dq_rad_sw_down_sfc

  !****************************************************
  
  
  SUBROUTINE dq_cwcabove_3D(meteoMarketPtr, met_src, valid_times, ifMakeCWC, ifMakePWC)
    !
    ! Cloud water column above bottom of meteo layers.
    ! Used for wet deposition calculation
    !
    ! Code owner: Rostislav Kouznetsov, FMI
    !
    ! All units: SI
    
    IMPLICIT NONE

    ! Imported parameters with intent IN:
    type(meteo_data_source), INTENT(in) :: met_src
    TYPE(silja_time), DIMENSION(:), INTENT(in) :: valid_times
    type(mini_market_of_stacks), pointer :: meteoMarketPtr
    logical, intent(in) :: ifMakeCWC, ifMakePWC

    ! Local declarations:
    TYPE(silja_time) :: time
    INTEGER :: i, tloop
    TYPE(silja_3d_field), POINTER :: cwc3d, cic3d
    type(silja_field), pointer :: fldTmp
    TYPE(silja_field_id) :: id
    TYPE(silja_grid) :: grid
    REAL, DIMENSION(:), POINTER :: cwc_above, pCWC_prev, pwc_above, pPWC_prev, psrf, cwc, cic

    REAL :: dag, dbg  !! sfc pressure to level airmass coefffs
    INTEGER :: fs, lev,  number_of_levels
    LOGICAL :: ifOnecwc, ifCWCinSM, ifPWCinSM
    real, dimension(max_levels) :: a_met, b_met, a_half_met, b_half_met

    REAL :: fTmp, fWe1, fWe2   ! tmp, weights for heights for level thickness

    IF (error) RETURN



    loop_over_times: DO tloop = 1, SIZE(valid_times)
      IF (.NOT.defined(valid_times(tloop))) EXIT loop_over_times
      time = valid_times(tloop)

      ifCWCinSM = fu_field_in_sm(meteoMarketPtr, met_src, &                       ! Already done before?
                       & cwcabove_3D_flag,&
                       & time,&
                       & level_missing,&
                       & .true.,&   ! 3D
                       & .false.)
      ifPWCinSM = fu_field_in_sm(meteoMarketPtr, met_src, &                       ! Already done before?
                       & pwcabove_3D_flag,&
                       & time,&
                       & level_missing,&
                       & .true.,&   ! 3D
                       & .false.)


      if ( ( (.not. ifMakeCWC) .or. ifCWCinSM) .and. ((.not. ifMakeCWC) .or. ifCWCinSM)) CYCLE loop_over_times

      !------------------------------------------------------
      !
      ! 1. Get necessary fields and their parameters
      !
      fldTmp => fu_sm_obstime_field(meteoMarketPtr, met_src, surface_pressure_flag, level_missing, &
                                  & time, single_time)
      if(error)return
      psrf => fu_grid_data(fldTmp)


      IF (fu_field_in_sm(meteoMarketPtr, met_src, &    !Only one cloud water?
                       & cloud_cond_water_flag, &
                       & time,&
                       & level_missing,&
                       & .true.,&   ! 3D
                       & .false.))  then
        cwc3d => fu_sm_obstime_3d_field(meteoMarketPtr, met_src,&         
                                    & cloud_cond_water_flag, &  
                                    & time,&
                                    & single_time)
        IF (error) CYCLE loop_over_times
        ifOnecwc = .true.
      else
        cwc3d => fu_sm_obstime_3d_field(meteoMarketPtr, met_src,&         
                                    & cloud_water_flag, &  
                                    & time,&
                                    & single_time)
        IF (error) CYCLE loop_over_times
       
        cic3d => fu_sm_obstime_3d_field(meteoMarketPtr, met_src,&         
                                    & cloud_ice_flag, &  
                                    & time,&
                                    & single_time)
        IF (error) CYCLE loop_over_times
        ifOnecwc = .false.
      endif

      IF (error) CYCLE loop_over_times

        
      grid = fu_grid(cwc3d) ! Any non-shifted field
      
      number_of_levels = fu_NbrOfLevels(meteo_vertical)
      IF (error) RETURN

      call hybrid_coefs(meteo_vertical, a_full=a_met, b_full=b_met) 
       !Indeed, a_half, b_half) should be used here
        
      fs = fu_number_of_gridpoints(grid)
      
    
      
      !-----------------------------------------------------
      !
      ! 2. Calculate cumulative water column above the 
      !    bottom of 'layers' and store it as precipitation rate for now
      !    'layer' boundaries are taken as mid-heights between levels
      !    or ground
      !
!      cwc_above = 0.

!!! Invert  meteo to half-level coefficients
        !prepare da db calvulation: attempt to reconstruct half-level coefficients
        ! FIXME Meteo vertical could have half-levels: grib files do have them
      a_half_met(1) = 0.
      b_half_met(1) = 1.
      DO lev = 1,number_of_levels
          a_half_met(lev+1) = a_half_met(lev) + 2*(a_met(lev)-a_half_met(lev))
          b_half_met(lev+1) = b_half_met(lev) + 2*(b_met(lev)-b_half_met(lev))
      enddo


      DO lev = number_of_levels, 1, -1 ! top->bottom

        if (.not. ifOnecwc) cic => fu_grid_data_from_3d(cic3d,lev)

        cwc => fu_grid_data_from_3d(cwc3d, lev)


        dag = a_half_met(lev) - a_half_met(lev+1) 
        dbg = b_half_met(lev) - b_half_met(lev+1)
        dag = dag / g  !! dag+dbg*psrf gives airmass  of meteo layer in kg/m2
        dbg = dbg / g
#ifdef DEBUG        
        if( dag+dbg*40000 < 0) then 
           call msg("minval(psrf)", minval(psrf(1:fs)))
           call set_error("Negative meteo layer mass", "dq_cwcabove_3D")
        endif
#endif        

        if (ifMakeCWC) then 
          id = fu_set_field_id(fu_met_src(cwc3d),&
                             & cwcabove_3D_flag,&
                             & fu_analysis_time(cwc3d),&
                             & fu_forecast_length(cwc3d), &
                             & grid,&
                             & fu_level(meteo_vertical, lev))
          call find_field_data_storage_2d(meteoMarketPtr, id, multi_time_stack_flag, cwc_above)
          if(fu_fails(.not.error,'Failed cwc_above field data pointer','dq_cwcabove_3D'))return

          if (lev == number_of_levels) then !! Init accumulation
             pCWC_prev => cwc_above
             cwc_above(1:fs) = 0.
          endif
          if (ifOnecwc) then  !Only cwc used
              !$omp workshare
                cwc_above(1:fs) = pCWC_prev(1:fs) + cwc(1:fs) * (dag + psrf(1:fs)*dbg)
              !$omp end workshare
          else ! cwc+cic
              !$omp workshare
                cwc_above(1:fs) = pCWC_prev(1:fs) + (cwc(1:fs)+cic(1:fs)) *  (dag + psrf(1:fs)*dbg)
              !$omp end workshare
          endif !ifOnecwc
          pCWC_prev => cwc_above
        endif
        if (ifMakePWC) then
          !!Set the pointer once
          id = fu_set_field_id(fu_met_src(cwc3d),&
                             & pwcabove_3D_flag,&
                             & fu_analysis_time(cwc3d),&
                             & fu_forecast_length(cwc3d), &
                             & grid,&
                             & fu_level(meteo_vertical, lev))
          call find_field_data_storage_2d(meteoMarketPtr, id, multi_time_stack_flag, pwc_above)
          if(fu_fails(.not.error,'Failed pwc_above field data pointer','dq_cwcabove_3D'))return

          if (lev == number_of_levels) then !! Init accumulation
             pPWC_prev => pwc_above
             pwc_above(1:fs) = 0.
          endif
          if (ifOnecwc) then  !Only cwc used
              !$omp workshare
                pwc_above(1:fs) = pPWC_prev(1:fs) + max(cwc(1:fs)-cloud_saturation, 0.) * (dag + psrf(1:fs)*dbg)
              !$omp end workshare
          else ! cwc+cic
              !$omp workshare
                pwc_above(1:fs) = pPWC_prev(1:fs) + max(cwc(1:fs)+cic(1:fs)-cloud_saturation, 0.) *  (dag + psrf(1:fs)*dbg)
              !$omp end workshare
          endif !ifOnecwc
          pPWC_prev => pwc_above
        endif



      END DO !!levels

      !!!Total column -- just copy the lowest level
      if (ifMakeCWC) then
          id = fu_set_field_id(fu_met_src(cwc3d),&
                             & cwcolumn_flag,&
                             & fu_analysis_time(cwc3d),&
                             & fu_forecast_length(cwc3d), &
                             & grid,&
                             & surface_level)
          call find_field_data_storage_2d(meteoMarketPtr, id, multi_time_stack_flag, cwc_above)
          if(fu_fails(.not.error,'Failed cwc_above_sfc field data pointer','dq_cwcabove_3D'))return
          cwc_above(1:fs) = pCWC_prev(1:fs)
      endif
      if (ifMakePWC) then
          id = fu_set_field_id(fu_met_src(cwc3d),&
                             & pwcolumn_flag,&
                             & fu_analysis_time(cwc3d),&
                             & fu_forecast_length(cwc3d), &
                             & grid,&
                             & surface_level)
          call find_field_data_storage_2d(meteoMarketPtr, id, multi_time_stack_flag, pwc_above)
          if(fu_fails(.not.error,'Failed cwc_above_sfc field data pointer','dq_cwcabove_3D'))return
          pwc_above(1:fs) = pPWC_prev(1:fs)
      endif


    END DO loop_over_times


!    CALL free_work_array(cwc_above)

  END SUBROUTINE dq_cwcabove_3D


  !****************************************************

  SUBROUTINE dq_scavenging_coef(meteoMarketPtr, met_src, valid_times, ifRandomise)
! Deprecated!

          !
    ! First-guess rough estimate of the scavenging coefficient regardless the
    ! reatures of the scavenged substance. 
    !
    ! The in-cloud/sub-cloud division is determined for a 'reference' cloud,
    ! with constant cloud top and bottom heights. The water/snow
    ! phase is chosen according to the air temperature at the particle
    ! elevation.
    !     
    ! "Instatenous" HIRLAM precipitation values are used
    ! whenever they are available. Otherwise the accumulated values are taken
    !
    ! The water/snow phase are split regarding the temperature at the particle height
    ! The scavenging coefficients are computed from the NAME-II type methods
    !
    ! Code owner Mikhail Sofiev
    ! Implemented algorithm is taken from diy_weather_tools of I.Valkama
    !
    ! All units: SI
    
    IMPLICIT NONE

    ! Imported parameters with intent IN:
    type(meteo_data_source), INTENT(in) :: met_src
    TYPE(silja_time), DIMENSION(:), INTENT(in) :: valid_times
    type(mini_market_of_stacks), pointer :: meteoMarketPtr
    logical, intent(in) :: ifRandomise

    ! Local declarations:
    TYPE(silja_time) :: time
    INTEGER :: i, tloop
    TYPE(silja_3d_field), POINTER :: t3d, p3d
    TYPE(silja_field), POINTER :: prec_ls, prec_conv
    TYPE(silja_field_id) :: id
    TYPE(silja_grid) :: grid
    TYPE(silja_level), DIMENSION(max_levels) :: levels
    REAL, DIMENSION(:), POINTER :: scav, pr_ls, pr_conv, t, p
    INTEGER :: fs, l, ix, number_of_levels
    REAL :: ls_scav, conv_scav

!!Deprecated
!

    call set_error("Deprecated code!","dq_scavenging_coef")
    return

    scav => fu_work_array()
    IF (error) RETURN

    loop_over_times: DO tloop = 1, SIZE(valid_times)
      IF (.NOT.defined(valid_times(tloop))) EXIT loop_over_times
      time = valid_times(tloop)

      IF (fu_field_in_sm(meteoMarketPtr, met_src, &                       ! Already done before?
                       & scavenging_coefficient_flag,&
                       & time,&
                       & level_missing,&
                       & .true.,&   ! 3D
                       & .false.)) CYCLE loop_over_times

      !------------------------------------------------------
      !
      ! 1. Get necessary fields and their parameters
      !
      t3d => fu_sm_obstime_3d_field(meteoMarketPtr, met_src,&         
                                  & temperature_flag, &  
                                  & time,&
                                  & single_time)
      IF (error) CYCLE loop_over_times

      p3d => fu_sm_obstime_3d_field(meteoMarketPtr, met_src,&         
                                  & pressure_flag, &  
                                  & time,&
                                  & single_time)
      IF (error) CYCLE loop_over_times

      prec_ls => fu_sm_obstime_field(meteoMarketPtr, met_src, &  ! Instantaneous field is required
                                   & large_scale_rain_int_flag, &  
                                   & level_missing, &
                                   & time,&
                                   & single_time)
      IF (error) CYCLE loop_over_times

      pr_ls => fu_grid_data(prec_ls)

      prec_conv => fu_sm_obstime_field(meteoMarketPtr, met_src, &  ! Instantaneous field is required
                                     & convective_rain_int_flag, &  
                                     & level_missing, &
                                     & time,&
                                     & single_time)
      IF (error) CYCLE loop_over_times

      pr_conv => fu_grid_data(prec_conv)
        
      grid = fu_grid(t3d) ! Any non-shifted field
      
      CALL vertical_levels(t3d, levels, number_of_levels)
      IF (error) RETURN
        
      fs = fu_number_of_gridpoints(grid)
 
      scav = real_missing
    
      !-----------------------------------------------------
      !
      ! 2. Compute scavenging for each layer
      !
         
      loop_over_levels: DO l = 1, number_of_levels

        t => fu_grid_data_from_3d(t3d, l)
        p => fu_grid_data_from_3d(p3d, l)

        loop_grid: DO ix = 1,fs

          !---------- Ensure no scavenging for fog or above the cloud top

          IF(pr_ls(ix)+pr_conv(ix) < 1.e-06 .or. p(ix) < 70000.)THEN 
            scav(ix)=0.0
            CYCLE loop_grid
          END IF

          !---------- If something exists - check the in/below cloud options
          !           Note precipitation unit: mm s-1
          
          IF(p(ix) > 90000.)THEN 

            !---  below the cloud base convective and dynamic efficiencies equal
            !
            IF (t(ix) > 273.16)THEN  !---rain
              conv_scav = pr_conv(ix)**0.79 * 0.05417  !8.40E-5*(3600**0.79)
              ls_scav = pr_ls(ix)**0.79 * 0.05417

            ELSE   !-----snow
              conv_scav = pr_conv(ix)**0.305 * 9.722e-04  ! 8.00E-5* 3600.**0.305
              ls_scav = pr_ls(ix)**0.305 * 9.722e-04
            END IF

          ELSE   
            !
            !------- for in-cloud scavenging (700 hPa < height < 900 hpa) 
            !
            ! 1. Convective rain efficiencies of water and snow are equal 
            !
            conv_scav = pr_conv(ix)**0.79 * 0.2167 ! 3.36E-4*(3600**0.79)

            ! 2. Dynamic (large-scale) rain/snow efficiencies
            !
            IF(t(ix) > 273.16)THEN  !------ rain
              ls_scav = pr_ls(ix)**0.79 * 0.05417  ! 8.40E-5*(3600**0.79)

            ELSE !----------- snow
              ls_scav = pr_ls(ix)**0.305 * 9.722e-04  ! 8.00E-5* 3600.**0.305
            END IF

          END IF ! In/sub cloud scavenging

          scav(ix) = ls_scav + conv_scav

        END DO loop_grid

        id = fu_set_field_id(fu_met_src(prec_ls),&
                           & scavenging_coefficient_flag,&
                           & fu_analysis_time(prec_ls),&
                           & fu_forecast_length(prec_ls), &
                           & fu_grid(prec_ls),&
                           & levels(l))
        
        CALL dq_store_2d(meteoMarketPtr, id, scav, multi_time_stack_flag, ifRandomise)

      END DO loop_over_levels
    END DO loop_over_times

    CALL free_work_array(scav)

  END SUBROUTINE dq_scavenging_coef

    !****************************************************

  SUBROUTINE dq_scavenging_coef_ls(meteoMarketPtr, met_src, valid_times, ifRandomise)
 !Deprecated
    !
    ! First-guess rough estimate of the scavenging coefficient regardless the
    ! reatures of the scavenged substance. 
    !
    ! This subroutine uses only large-scale rain fields.
    !
    ! The in-cloud/sub-cloud division is determined for a 'reference' cloud,
    ! with constant cloud top and bottom heights. The water/snow
    ! phase is chosen according to the air temperature at the particle
    ! elevation.
    !     
    ! "Instatenous" HIRLAM precipitation values are used
    ! whenever they are available. Otherwise the accumulated values are taken
    !
    ! The water/snow phase are split regarding the temperature at the particle height
    ! The scavenging coefficients are computed from the NAME-II type methods
    !
    ! Code owner Mikhail Sofiev
    ! Implemented algorithm is taken from diy_weather_tools of I.Valkama
    !
    ! All units: SI
    
    IMPLICIT NONE

    ! Imported parameters with intent IN:
    type(meteo_data_source), INTENT(in) :: met_src
    TYPE(silja_time), DIMENSION(:), INTENT(in) :: valid_times
    type(mini_market_of_stacks), pointer :: meteoMarketPtr
    logical, intent(in) :: ifRandomise

    ! Local declarations:
    TYPE(silja_time) :: time
    INTEGER :: i, tloop
    TYPE(silja_3d_field), POINTER :: t3d, p3d
    TYPE(silja_field), POINTER :: prec_ls, prec_conv
    TYPE(silja_field_id) :: id
    TYPE(silja_grid) :: grid
    TYPE(silja_level), DIMENSION(max_levels) :: levels
    REAL, DIMENSION(:), POINTER :: scav, pr_ls, t, p
    INTEGER :: fs, l, ix, number_of_levels
    REAL :: ls_scav, conv_scav

!!Deprecated
!

    call set_error("Deprecated code!","dq_scavenging_coef")
    return


    scav => fu_work_array()
    IF (error) RETURN

    loop_over_times: DO tloop = 1, SIZE(valid_times)
      IF (.NOT.defined(valid_times(tloop))) EXIT loop_over_times
      time = valid_times(tloop)

      IF (fu_field_in_sm(meteoMarketPtr, met_src, &                       ! Already done before?
                       & scavenging_coefficient_flag,&
                       & time,&
                       & level_missing,&
                       & .true.,&   ! 3D
                       & .false.)) CYCLE loop_over_times

      !------------------------------------------------------
      !
      ! 1. Get necessary fields and their parameters
      !
      t3d => fu_sm_obstime_3d_field(meteoMarketPtr, met_src,&         
                                  & temperature_flag, &  
                                  & time,&
                                  & single_time)
      IF (error) CYCLE loop_over_times

      p3d => fu_sm_obstime_3d_field(meteoMarketPtr, met_src,&         
                                  & pressure_flag, &  
                                  & time,&
                                  & single_time)
      IF (error) CYCLE loop_over_times

      prec_ls => fu_sm_obstime_field(meteoMarketPtr, met_src, &  ! Instantaneous field is required
                                   & large_scale_rain_int_flag, &  
                                   & level_missing, &
                                   & time,&
                                   & single_time)
      IF (error) CYCLE loop_over_times

      pr_ls => fu_grid_data(prec_ls)
        
      grid = fu_grid(t3d) ! Any non-shifted field
      
      CALL vertical_levels(t3d, levels, number_of_levels)
      IF (error) RETURN
        
      fs = fu_number_of_gridpoints(grid)
 
      scav = real_missing
    
      !-----------------------------------------------------
      !
      ! 2. Compute scavenging for each layer
      !
         
      loop_over_levels: DO l = 1, number_of_levels

        t => fu_grid_data_from_3d(t3d, l)
        p => fu_grid_data_from_3d(p3d, l)

        loop_grid: DO ix = 1,fs

          !---------- Ensure no scavenging for fog or above the cloud top

          IF(pr_ls(ix) < 1.e-06 .or. p(ix) < 70000.)THEN 
            scav(ix)=0.0
            CYCLE loop_grid
          END IF

          !---------- If something exists - check the in/below cloud options
          !           Note precipitation unit: mm s-1
          
          IF(p(ix) > 90000.)THEN 

            !---  below the cloud base convective and dynamic efficiencies equal
            !
            IF (t(ix) > 273.16)THEN  !---rain
              ls_scav = pr_ls(ix)**0.79 * 0.05417

            ELSE   !-----snow
              ls_scav = pr_ls(ix)**0.305 * 9.722e-04

            END IF

          ELSE   
            !
            !------- for in-cloud scavenging (700 hPa < height < 900 hpa) 
            !
            ! Dynamic (large-scale) rain/snow efficiencies
            !
            IF(t(ix) > 273.16)THEN  !------ rain
              ls_scav = pr_ls(ix)**0.79 * 0.05417  ! 8.40E-5*(3600**0.79)

            ELSE !----------- snow
              ls_scav = pr_ls(ix)**0.305 * 9.722e-04  ! 8.00E-5* 3600.**0.305
            END IF

          END IF ! In/sub cloud scavenging

          scav(ix) = ls_scav

        END DO loop_grid

        id = fu_set_field_id(fu_met_src(prec_ls),&
                           & scavenging_coefficient_flag,&
                           & fu_analysis_time(prec_ls),&
                           & fu_forecast_length(prec_ls), &
                           & fu_grid(prec_ls),&
                           & levels(l))
  
        CALL dq_store_2d(meteoMarketPtr, id, scav, multi_time_stack_flag, ifRandomise)

      END DO loop_over_levels
    END DO loop_over_times

    CALL free_work_array(scav)

  END SUBROUTINE dq_scavenging_coef_ls


  ! ***************************************************************


  SUBROUTINE dq_pasquill(meteoMarketPtr, met_src, valid_times)
    ! 
    ! Makes Pasquill stability class field. Requires:
    ! Monin-Obukhov length and surface roughness for wind
    ! 
    ! Method: Golder (1977), modified by Pentti Vaajama, FMI
    ! 
    ! All units: SI
    !
    ! Language: ANSI Fortran 90
    !
    ! Current code owner: Mikhail Sofiev,
    ! Elements are taken from diu_weather_tools of I.Valkama

    IMPLICIT NONE

    ! Imported parameters with intent IN:
    type(meteo_data_source), INTENT(in) :: met_src
    TYPE(silja_time), DIMENSION(:), INTENT(in) :: valid_times
    type(mini_market_of_stacks), pointer :: meteoMarketPtr

    ! Local declarations:
    TYPE(silja_time) :: time
    INTEGER :: i, tloop
    TYPE(silja_field), POINTER :: MO_l_inv_f, sr_f
    TYPE(silja_field_id) :: id
    REAL, DIMENSION(:), POINTER :: MO_l_inv, sr, stab_data
    INTEGER :: ix
    REAL :: local_MO, local_pasquill

!    stab_data => fu_work_array()
    IF (error) RETURN

    loop_over_times: DO tloop = 1, SIZE(valid_times)
      IF (.NOT.defined(valid_times(tloop))) EXIT loop_over_times
      time = valid_times(tloop)

      IF (fu_field_in_sm(meteoMarketPtr, met_src, & ! Already done before?
                       & pasquill_class_flag,&
                       & time,&
                       & level_missing,&
                       & .false.,&   ! 3D
                       & .false.)) CYCLE loop_over_times

      !------------------------------------------------------
      !
      ! 1. Get necessary fields and their parameters
      !
      MO_l_inv_f => fu_sm_obstime_field(meteoMarketPtr, met_src, &  ! Monin Obukhov_length
                                      & MO_length_inv_flag, &  
                                      & level_missing, &
                                      & time,&
                                      & single_time)
      IF (error) THEN
        CALL unset_error('dq_pasquill')
        CYCLE loop_over_times
      END IF

      MO_l_inv => fu_grid_data(MO_l_inv_f)

      !-- Does dynamic / permanent roughness exist ?
      !
      IF(fu_field_in_sm(meteoMarketPtr, met_src, surface_roughness_meteo_flag,& !dynamical roughness ?
                      & time, level_missing,&
                      & .false.,.false.)) THEN    ! look for 3d; if_permanent

        sr_f => fu_sm_obstime_field(meteoMarketPtr, &  ! Surface roughness for wind
                                  & met_src,&         
                                  & surface_roughness_meteo_flag, &  
                                  & level_missing, &
                                  & time,&
                                  & single_time)
        sr => fu_grid_data(sr_f) ! Dynamical roughness

      ELSEIF(fu_field_in_sm(meteoMarketPtr, met_src, surface_roughness_meteo_flag,& ! permanent ?
                          & time, level_missing,&
                          & .false.,.true.)) THEN     ! look for 3d; if_permanent

        sr_f => fu_sm_simple_field(meteoMarketPtr, met_src,&
                                 & surface_roughness_meteo_flag,&
                                 & level_missing, &
                                 & single_time_stack_flag)
        sr => fu_grid_data(sr_f) ! permanent roughness

      else ! No roughness at all - use something default
        IF(error) CALL unset_error('dq_pasquill')
        CALL msg_warning('No roughness at all. Take default','dq_pasquill')

        sr => fu_work_array()
        sr = 0.05
      END IF


      IF (error) RETURN

      !------------------------------------------------------
      !
      ! 2. Now - compute the stability class
      !
      id = fu_set_field_id(fu_met_src(MO_l_inv_f),&
                         & pasquill_class_flag,&
                         & fu_analysis_time(MO_l_inv_f),&
                         & fu_forecast_length(MO_l_inv_f), &
                         & meteo_grid,&
                         & ground_level)
      call find_field_data_storage_2d(meteoMarketPtr, id, multi_time_stack_flag, stab_data)
      if(fu_fails(.not.error,'Failed stab_data field data pointer','dq_pasquill'))return

      loop_grid: DO ix = 1,fs_meteo

        local_MO = ABS(MO_l_inv(ix))

        local_pasquill = local_MO*(14.6 - 0.167&
              & *sr(ix)) + (1.5979 + 0.21625*LOG(sr(ix)))&
              & *(1.0 - EXP(local_MO*(-47.843 - 178.46*sr(ix))))

        local_pasquill = SIGN(1.0,MO_l_inv(ix))*local_pasquill

        IF(local_pasquill <= -2.5)THEN
          stab_data(ix) = 1
        ELSE IF(local_pasquill <= -1.5)THEN
          stab_data(ix) = 2
        ELSE IF(local_pasquill <= -0.5)THEN
          stab_data(ix) = 3
        ELSE IF(local_pasquill <= 0.5)THEN
          stab_data(ix) = 4
        ELSE IF(local_pasquill <= 1.5)THEN
          stab_data(ix) = 5
        ELSE IF(local_pasquill <= 2.5)THEN
          stab_data(ix) = 6
        ELSE
          stab_data(ix) = 7
        END IF

      END DO loop_grid


!      CALL dq_store_2d(meteoMarketPtr, id, stab_data, multi_time_stack_flag )

    END DO loop_over_times

!    CALL free_work_array(stab_data)

  END SUBROUTINE dq_pasquill


  !*****************************************************************************
  
  SUBROUTINE dq_lscale_profile(meteoMarketPtr, met_src, valid_times, Kz_method, ABLh_min)
    ! Lscale for mixing
    ! Follows ECMWF procedure. Intent to be used for output only
    ! Repeats  dq_Rdown_profile
    !
    IMPLICIT NONE
    !
    ! Imported parameters with intent(in):
    type(meteo_data_source), INTENT(in) :: met_src
    TYPE(silja_time), DIMENSION(:), INTENT(in) :: valid_times
    INTEGER, INTENT(in) :: Kz_method
    real, INTENT(in) :: ABLh_min
    type(mini_market_of_stacks), pointer :: meteoMarketPtr

    ! Imported parameters with intent(out):

    ! Local declarations:
    TYPE(silja_time) :: time
    INTEGER :: tloop, lev, i, iFlip
    TYPE(silja_field_id) :: id
    type(silja_field), pointer :: fldTmp
    TYPE(silja_3d_field), POINTER ::  height3d,  u3d, v3d, t3d, q3d !3D input
    REAL, DIMENSION(:), POINTER :: ptrHabl, ptrKz_1m                       !2D input
    REAL, DIMENSION(:), POINTER :: height1, height2, u_1, u_2, v_1, v_2, theta_1, theta_2
    REAL, DIMENSION(:), POINTER :: t2, q2, psrf, pml, R_tmp
    type(silja_rp_1d), DIMENSION(:), POINTER :: fliparr
    TYPE(silja_level), DIMENSION(max_levels) :: levels
    INTEGER :: number_of_levels
    real ::  fTmp, Habl_01, Habl, Ri, du, dv, dz, deltaU2, deltaT, l, PhiMPhiH, z, zprev 
    real:: extrashear,  Rscale
    real, parameter :: eps = molecular_weight_air/molecular_weight_water - 1.

!    call set_error("Not implemented yet","dq_R1m_profile")
!    return
      !
      ! Get the space and define the flipping pointers
      !
    call get_work_arrays_set(2, fs_meteo, fliparr)
    iFlip = 1
!    R_tmp => fu_work_array(fs_meteo)
    pml => fu_work_array(fs_meteo)

    nullify(height2)
    nullify(u_2)
    nullify(v_2)
    nullify(theta_1)

    IF (error) RETURN

    DO tloop = 1, SIZE(valid_times)

      IF (.NOT.defined(valid_times(tloop))) EXIT   !loop_over_times
      time = valid_times(tloop)

      IF (fu_field_in_sm(meteoMarketPtr, met_src,&   ! Already done before?
                       & turb_length_scale_flag,&
                       & time,&
                       & level_missing,&
                       & .true.,&
                       & .false.)) CYCLE  ! loop_over_times

      ! 2D fields

      fldTmp => fu_sm_obstime_field(meteoMarketPtr, met_src, surface_pressure_flag, level_missing, &
                                  & time, single_time)
      if(error)return
      psrf => fu_grid_data(fldTmp)


      ! 3D fields
      height3d => fu_sm_obstime_3d_field(meteoMarketPtr, met_src, height_flag, time, single_time)
      IF (error) RETURN

      CALL vertical_levels(height3d, levels, number_of_levels)
      IF (error) RETURN
      

      u3d => fu_sm_obstime_3d_field(meteoMarketPtr, met_src, u_flag, time, single_time)
      IF (error) RETURN

      v3d => fu_sm_obstime_3d_field(meteoMarketPtr, met_src, v_flag, time, single_time)
      IF (error) RETURN

      t3d => fu_sm_obstime_3d_field(meteoMarketPtr, met_src, temperature_flag, time, single_time)
      IF (error) RETURN

      
      q3d => fu_sm_obstime_3d_field(meteoMarketPtr, met_src, specific_humidity_flag, &
                                      & time, single_time)
      IF (error .or. .not. associated(q3d))then
        call msg_warning("Can't get humidity, using dry atmosphere","dq_lscale_profile")

        call unset_error('dq_lscale_profile')
        nullify(q3d)
      endif
     
       
      call msg("Making dq_lscale_profile. ECMWF method")
      DO lev = 1, number_of_levels
        height1 => height2
        u_1 => u_2
        v_1 => v_2
        theta_1 => theta_2
        
        height2 => fu_grid_data_from_3d(height3d, lev)
        u_2 => fu_grid_data_from_3d(u3d, lev)
        v_2 => fu_grid_data_from_3d(v3d, lev)
        if (associated(q3d)) q2 => fu_grid_data_from_3d(q3d, lev)
        t2 => fu_grid_data_from_3d(t3d, lev)
        CALL pressure_on_level(t3d, lev, pml)  ! It really computes the values

        theta_2 => fliparr(iFlip)%pp
        if (associated(q3d)) then
              theta_2(1:fs_meteo) = t2(1:fs_meteo) * &
                   & ((std_pressure_sl/pml(1:fs_meteo))**R_per_c_dryair) &
                   & * (1.-eps*q2(1:fs_meteo))
        else
              theta_2(1:fs_meteo) = t2(1:fs_meteo) * &
                   & ((std_pressure_sl/pml(1:fs_meteo))**R_per_c_dryair)
        endif
        
        id = fu_set_field_id(fu_met_src(height3d),&
                           & turb_length_scale_flag,&
                           & fu_analysis_time(height3d),&
                           & fu_forecast_length(height3d), &
                           & fu_grid(height3d),&
                           & levels(lev))
        call find_field_data_storage_2d(meteoMarketPtr, id, multi_time_stack_flag, R_tmp)
        if(fu_fails(.not.error,'Failed R_tmp field data pointer','dq_lscale_profile'))return

        if (lev == 1) then ! Just some scale
            R_tmp(1:fs_meteo) = height2(1:fs_meteo)
        else

           do i=1,fs_meteo
                   
               z = height2(i)

               du = u_2(i) - u_1(i)
               dv = v_2(i) - v_1(i)
               dz = height2(i) - height1(i)
               deltaT = theta_2(i) - theta_1(i)

               l = 1./(1./150. + 1./height2(i)) ! Magic ECMWF scale

               ! Unresolved shear as in ECMWF IFS Cy38r1, part 4, eq. 3.57.
               ! there extrashear is added to sqrt(deltaU2),
               ! we add its square to squared shear 
               fTmp = (1. - pml(i)/psrf(i))
               extrashear = 5e-3 * 27.2 * fTmp * exp(-10*fTmp) * dz

               deltaU2 = du*du + dv*dv + extrashear*extrashear +  0.01  
                                 !ensure deltaU2 > 0.01  in any case

               Ri = 2. * g * dz * deltaT / ( (theta_1(i)+theta_2(i))* deltaU2)
               if (Ri > 0.) then
                       ! Approximation of ECMWF gradient stability functions
                       ! PhiMPhiH(Ri=0) = 1
                       ! PhiMPhiH(Ri -> \infty) ~ Ri^3 
                       ! Kz = l^2 |dU/dz| / (PhiMPhiH) 
                       ! p. 47
                  PhiMPhiH = ( (305.*Ri + 78.)*Ri + 13.)*Ri + 1.
               else ! Neutral/unstable...
                   PhiMPhiH = 1./sqrt(sqrt(1 - 16.* Ri)) !PhiMPhiH = (1 - 16. Ri)^-0.75
                   PhiMPhiH = PhiMPhiH * PhiMPhiH * PhiMPhiH
               endif

               R_tmp(i) = l / sqrt(PhiMPhiH)  ! ECMWF length scaled with stratification
           end do ! loop meteo
        endif !(lev == 1)
        IF (error) RETURN 

!        CALL dq_store_2d(meteoMarketPtr, id, R_tmp, multi_time_stack_flag )

        iFlip = 3 - iFlip
        IF (error) RETURN 
      END DO ! loop_over_levels

    END DO ! loop_over_times

!    call free_work_array(R_tmp)
    call free_work_array(pml)
    call free_work_array(fliparr)

  END SUBROUTINE dq_lscale_profile

  !*****************************************************************************
  !*****************************************************************************
  
  SUBROUTINE dq_Rdown_profile(meteoMarketPtr, met_src, valid_times, Kz_method, ABLh_min)
    ! Aerodynamic resistane between this and previous level
    ! flux between levels -- (C1 - C2) / (R2)
    ! Resisratnce to 1 m according to SILAM standard Kz for the lowest level
    ! The idea is to calculate everything from the basic parameters,
    ! So only one 3D field is needed for Kz
    !
    IMPLICIT NONE
    !
    ! Imported parameters with intent(in):
    type(meteo_data_source), INTENT(in) :: met_src
    TYPE(silja_time), DIMENSION(:), INTENT(in) :: valid_times
    INTEGER, INTENT(in) :: Kz_method
    real, INTENT(in) :: ABLh_min
    type(mini_market_of_stacks), pointer :: meteoMarketPtr

    ! Imported parameters with intent(out):

    ! Local declarations:
    TYPE(silja_time) :: time
    INTEGER :: tloop, lev, i, iFlip
    TYPE(silja_field_id) :: id
    type(silja_field), pointer :: fldTmp
    TYPE(silja_3d_field), POINTER ::  height3d,  u3d, v3d, t3d, q3d !3D input
    REAL, DIMENSION(:), POINTER :: ptrHabl, ptrKz_1m, ptrPs, ptrMOs, ptrWstar         !2D input
    REAL, DIMENSION(:), POINTER :: height1, height2, u_1, u_2, v_1, v_2, theta_prev
    REAL, DIMENSION(:), POINTER :: t2, q2, psrf, R_tmp
    REAL, DIMENSION(:), POINTER :: a_full, b_full, fWork
    type(silja_rp_1d), DIMENSION(max_levels) :: Rptr
    TYPE(silja_level), DIMENSION(max_levels) :: levels
    INTEGER :: number_of_levels, ithread, nthreads,  iTmp
    real ::  fTmp, Habl_01, Habl, Ri, du, dv, dz, pml, deltaU2, deltaT, l, PhiMPhiH, z, zprev, theta 
    real:: extrashear,  Rscale, kZ_tab
    real, parameter :: eps = molecular_weight_air/molecular_weight_water - 1.
    integer, parameter  :: KZtabLen = 50
    real, dimension(KZtabLen) :: HuntenKZ 
    real, parameter :: HuntenKZpar(KZtabLen) = (/10.0,10.0,10.0,10.0,10.0,10.0,10.0,10.0,6.5,3.0,3.0,3.0,1.6,.23,.23,.27,.30,.33,.37,.42,.45,.52,.56,.64,.72,.80,.90,1.0,1.1,1.2,1.3,1.5,1.6,1.8,2.0,2.2,2.4,2.8,3.0,3.4,3.8,4.1,4.7,5.2,5.9,6.3,7.2,8.0,9.0,10.0/)
    ! Hunten 1975 Index -- height in km
    ! JOURNAL OF GEOPHYSICAL RESEARCH, VOL. 86, NO. CI0, PAGES 9859-9868, OCTOBER 20, 1981
    ! Stratospheric Eddy Diffusion Coefficients From   Tracer  Data
    !                 S. T. MASSIE AND  D. M. HUNTEN
    HuntenKZ(:) = HuntenKZpar(:) 

!    call set_error("Not implemented yet","dq_R1m_profile")
!    return

    nullify(height2)
    nullify(u_2)
    nullify(v_2)
    nullify(theta_prev)

    IF (error) RETURN

    DO tloop = 1, SIZE(valid_times)

      IF (.NOT.defined(valid_times(tloop))) EXIT   !loop_over_times
      time = valid_times(tloop)

      IF (fu_field_in_sm(meteoMarketPtr, met_src,&   ! Already done before?
                       & R_down_meteo_flag,&
                       & time,&
                       & level_missing,&
                       & .true.,&
                       & .false.)) CYCLE  ! loop_over_times

      if (Kz_method /= zero_kz) then ! some meteo needed

              ! 2D fields
              fldTmp => fu_sm_obstime_field(meteoMarketPtr, met_src, abl_height_m_pointer, level_missing, &
                                          & time, single_time)
              if(error)return
              ptrHabl => fu_grid_data(fldTmp)

              fldTmp => fu_sm_obstime_field(meteoMarketPtr, met_src, Kz_scalar_1m_flag, level_missing, &
                                                 & time, single_time)
              if(error)return
              ptrKz_1m => fu_grid_data(fldTmp)

              fldTmp => fu_sm_obstime_field(meteoMarketPtr, met_src, surface_pressure_flag, level_missing, &
                                          & time, single_time)
              if(error)return
              psrf => fu_grid_data(fldTmp)

              fldTmp => fu_sm_obstime_field(meteoMarketPtr, met_src, convective_velocity_scale_flag, level_missing, &
                                          & time, single_time)
              if(error)return
              ptrWstar => fu_grid_data(fldTmp)

              fldTmp => fu_sm_obstime_field(meteoMarketPtr, met_src, MO_length_inv_flag, level_missing, &
                                          & time, single_time)
              if(error)return
              ptrMOs => fu_grid_data(fldTmp) 
      endif


      ! 3D fields
      height3d => fu_sm_obstime_3d_field(meteoMarketPtr, met_src, height_flag, time, single_time)
      IF (error) RETURN

      CALL vertical_levels(height3d, levels, number_of_levels)
      IF (error) RETURN
      
      if ( any(Kz_method == (/silam_abl_ec_ft_kz, ec_kz, simple_abl_ec_ft_kz/))) then ! Additional fields needed

        u3d => fu_sm_obstime_3d_field(meteoMarketPtr, met_src, u_flag, time, single_time)
        IF (error) RETURN

        v3d => fu_sm_obstime_3d_field(meteoMarketPtr, met_src, v_flag, time, single_time)
        IF (error) RETURN

        t3d => fu_sm_obstime_3d_field(meteoMarketPtr, met_src, temperature_flag, time, single_time)
        IF (error) RETURN

        q3d => fu_sm_obstime_3d_field(meteoMarketPtr, met_src, specific_humidity_flag, &
                                              & time, single_time)
        IF (error .or. .not. associated(q3d))then
          call msg_warning("Can't get humidity, using dry atmosphere","dq_Rdown_profile")
          call unset_error('dq_Rto1m_profile')
          nullify(q3d)
        endif
      endif

      DO lev = 1, number_of_levels
        id = fu_set_field_id(fu_met_src(height3d),&
                           & R_down_meteo_flag,&
                           & fu_analysis_time(height3d),&
                           & fu_forecast_length(height3d), &
                           & fu_grid(height3d),&
                           & levels(lev))
        call find_field_data_storage_2d(meteoMarketPtr, id, multi_time_stack_flag, Rptr(lev)%pp)
        if(fu_fails(.not.error,'Failed R_tmp field data pointer','dq_Rdown_profile'))return
      enddo

      select case (Kz_method)
          case (zero_kz)
            if(lev==1)call msg("Making Rdown for zero exchange. Method:"+fu_str(Kz_method))
            !$OMP PARALLEL shared(Rptr,number_of_levels,fs_meteo) &
            !$OMP & private (lev)
            !$OMP DO
            DO lev = 1, number_of_levels
               Rptr(lev)%pp(1:fs_meteo) = 1e20 ! Some very big number
            ENDDO
            !$OMP ENDDO
            !$OMP END PARALLEL

          case (silam_abl_ec_ft_kz, ec_kz, simple_abl_ec_ft_kz)
            theta_prev => fu_work_array(fs_meteo)
            fWork => fu_work_array(2*number_of_levels)
            a_full(1:number_of_levels) => fWork(1:number_of_levels)
            b_full(1:number_of_levels) => fWork(number_of_levels+1:2*number_of_levels)
            call hybrid_coefs(meteo_vertical, a_full=a_full, b_full=b_full)

            call msg("Making Rdown from Ri and surface scheme. Method:"+fu_str(Kz_method)+"ABL limit:",ABLh_min)

            !$OMP PARALLEL IF (.TRUE.) DEFAULT(NONE) &
            !$OMP & SHARED(Kz_method, height3d, fs_meteo,&
            !$OMP &    number_of_levels,ablh_min, &
            !$OMP &    u3d,v3d,q3d,t3d,error,Rptr, &
            !$OMP &    a_full,b_full,psrf,ptrhabl,ptrkz_1m, ptrMOs, ptrWstar, theta_prev) &
            !$OMP & PRIVATE(height1, height2,u_1,u_2,v_1,v_2,R_tmp ,i,z,Habl,Habl_01, lev,  pml,theta,ftmp,du,dv,dz,deltat,l,extrashear, &
            !$OMP &      t2,q2,deltau2,ri,phimphih,ithread, nthreads, zprev)
            ithread = 0
            nthreads = 1
            q2 => null()
            !$ ithread  = omp_get_thread_num()
            !$ nthreads = omp_get_num_threads()
            zprev= Kz_ref_height ! For the first level
            do lev=1, number_of_levels
               !set pointers
               ! These pointers must not be used at the lowest level
               height1 => height2
               u_1 => u_2
               v_1 => v_2
                
               height2 => fu_grid_data_from_3d(height3d, lev)
               u_2 => fu_grid_data_from_3d(u3d, lev)
               v_2 => fu_grid_data_from_3d(v3d, lev)
               if (associated(q3d)) q2 => fu_grid_data_from_3d(q3d, lev)
               t2 => fu_grid_data_from_3d(t3d, lev)

               R_tmp => Rptr(lev)%pp

               do i=1+ithread,fs_meteo, nthreads
                  if (error) exit
                  pml = a_full(lev) + b_full(lev)*psrf(i)
                  theta = t2(i) *  ((std_pressure_sl/pml)**R_per_c_dryair)

                  if (associated(q3d)) theta = theta  * (1.-eps*q2(i))
                      
                  z = height2(i)
                  if (lev > 1) zprev=height1(i)
                  Habl = max(ptrHabl(i), ABLh_min)
                  Habl_01 = Habl * 0.1
                  if ( lev == 1 ) then ! use standard SILAM PBL for lowest level

                      if(z < Habl_01)then      !  within 10% of ABL
                        fTmp = log(z/Kz_ref_height)
                      elseif (z < Habl) then ! upper within ABL
                        fTmp = (log(Habl_01/Kz_ref_height) + (z - Habl_01) / Habl_01)
                      else ! outside ABL -- only for 1st level....
                        fTmp =  (log(Habl_01/Kz_ref_height) + 9. & !(Habl - Habl_01) / Habl_01
                            & +      (z - Habl)/Habl_01  * Kz_FT_factor )
                      endif
                      fTmp = fTmp * Kz_ref_height / ptrKz_1m(i)

                  elseif (z < Habl) then 
                      if ( Kz_method == silam_abl_ec_ft_kz) then  !old  routine within ABL
                        if(z < Habl_01)then      ! both  within 10% of ABL
                          fTmp = log(z/height1(i))  !dimensionless resistance
                        elseif (height1(i) < Habl_01) then
                          fTmp = (log(Habl_01/height1(i)) + (z - Habl_01) / Habl_01)
                        else  !both  within 10% to 100% of ABL
                          fTmp = (z - height1(i)) / Habl_01
                        endif
                        !Scale it with ptrKz_1m
                        fTmp = fTmp * Kz_ref_height / ptrKz_1m(i)
                      else  !  Kz_method == simple_abl_ec_ft_kz
                       if ( ptrMOs(i) > 0) then
                         ! Kz for Stable BL (m2/s) 
                         fTmp = ptrKz_1m(i)/ Kz_ref_height / (1+4.7*z*ptrMOs(i))  *z* (1-z/Habl)**2 
                       else
                         !Unstable BL  
                         !(0.1*2.8) **1/3
                         fTmp = max(ptrKz_1m(i), karmann_c*0.65*ptrWstar(i)) *z* (1-z/Habl)**2
                       endif
                       fTmp = (z - zprev) / max(0.1,min(fTmp,500.)) !Resistance Dz/Kz

                      endif


                  else ! Upper layers -- ECMWF approach
                    du = u_2(i) - u_1(i)
                    dv = v_2(i) - v_1(i)
                    dz = height2(i) - height1(i)
                    deltaT = theta - theta_prev(i)

                    l = 1./(1./150. + 1./height2(i)) ! Magic ECMWF scale

                   ! Unresolved shear as in ECMWF IFS Cy38r1, part 4, eq. 3.57.
                   ! there extrashear is added to sqrt(deltaU2),
                   ! we add its square to squared shear 
                    fTmp = (1. - pml/psrf(i))
                    extrashear = 5e-3 * 27.2 * fTmp * exp(-10*fTmp) * dz

                    deltaU2 = du*du + dv*dv + extrashear*extrashear +  0.01  
                                      !ensure deltaU2 > 0.01  in any case

                    Ri = 2. * g * dz * deltaT / ( (theta + theta_prev(i))* deltaU2)
                    if (Ri > 0.) then
                            ! Approximation of ECMWF gradient stability functions
                            ! PhiMPhiH(Ri=0) = 1
                            ! PhiMPhiH(Ri -> \infty) ~ Ri^3 
                            ! Kz = l^2 |dU/dz| / (PhiMPhiH) 
                            ! p. 47
                       PhiMPhiH = ( (305.*Ri + 78.)*Ri + 13.)*Ri + 1.
                    else ! Neutral/unstable...
                        PhiMPhiH = 1./sqrt(sqrt(1 - 16.* Ri)) !PhiMPhiH = (1 - 16. Ri)^-0.75
                        PhiMPhiH = PhiMPhiH * PhiMPhiH * PhiMPhiH
                    endif

                     fTmp = PhiMPhiH * dz /(l*l*sqrt(deltaU2)) !1/kZ
                     ftmp = max( min(fTmp, 300.), 0.01)
                     !Minimum value Kz =0.003 m^2/s applied
                     ! 10^6/0.003 = 3*10^8s ~ 1 year turnover for 1-km layer,
                     ! 100y for 10km....

                     !11.9.2013 Reduced to 100 m^2/s so diffusion does not break
                     !! Max value Kx = 1000 m^2/s
                     !! 10^6/1000 = 10^3 s ~ 15 min turnover for 1-km layer
                     fTmp = dz *fTmp !dz/Kz

                     !!FIXME HACK
                     !!fTmp = dz * 0.1 ! Constant 10m^2/s Kz

                  endif ! ECMWF Kz
                  R_tmp(i) = fTmp  !
                  theta_prev(i) = theta !Save for the next level
                  if(r_tmp(i) <0)then
                    call msg('negative R silamAbl-ecFT: level, i',lev,i)
                    call msg('z,dz,Habl,Habl_01',(/z,dz,Habl,Habl_01/))
                    call msg('fTmp, ptrKz_1m(i), Kz_ref_height, R_tmp(i)', &
                           & (/fTmp,ptrKz_1m(i),Kz_ref_height,R_tmp(i)/))
                    call set_error('Negative resistance','dq_Rto1m_profile')
                  endif
              end do ! loop meteo
              if (error) exit
            end do !levels
            !$OMP END PARALLEL 
            call free_work_array(theta_prev)
            call free_work_array(fWork)

          case (simple_kz) !Impementation from eq 11-13 from GMD
            call msg("Making Rdown with simple ABL. Method:"+fu_str(Kz_method)+"ABL limit:",ABLh_min)
            ! Simple and stupid
            ! Take Kz at the level and divide it with dz to the previous level
            !$OMP PARALLEL IF (.TRUE.) DEFAULT(NONE) &
            !$OMP & SHARED(height3d, fs_meteo, number_of_levels,ablh_min, &
            !$OMP &   error,Rptr,ptrhabl,ptrkz_1m, ptrMOs, ptrWstar) &
            !$OMP & PRIVATE(height1, height2, R_tmp ,i,z,Habl, lev,ftmp,ithread, nthreads,  zprev )
            ithread = 0
            nthreads = 1
            !$ ithread  = omp_get_thread_num()
            !$ nthreads = omp_get_num_threads()
            zprev= Kz_ref_height ! For the first level
            do lev=1, number_of_levels
              R_tmp => Rptr(lev)%pp
              height1 => height2 !Warning! uninitialized at first level
              height2 => fu_grid_data_from_3d(height3d, lev)
              do i=1+ithread,fs_meteo, nthreads
                 if (error) cycle
                 z = height2(i)
                 if (lev > 1) zprev=height1(i)
                  
                 Habl = max(ptrHabl(i), ABLh_min)
                 if (z >= Habl) then
                   fTmp = 0.1  !Kz = 0.1 m2/s
                 elseif ( ptrMOs(i) > 0) then
                   ! Kz for Stable BL
                   fTmp = ptrKz_1m(i)/(1+4.7*z*ptrMOs(i))  *z* (1-z/Habl)**2 
                 else
                   !Unstable BL 
                   !(0.1*2.8) **1/3
                   fTmp = max(ptrKz_1m(i), karmann_c*0.65*ptrWstar(i)) *z* (1-z/Habl)**2
                 endif
                 R_tmp(i) = (z - zprev) / max(0.1,min(fTmp,500.))
              end do ! loop meteo
            END DO ! loop_over_levels

            !$OMP END PARALLEL
            IF (error) RETURN 



          case (silam_kz_emulator, hunten_kz) ! Just the same as old KZ from Kz1m and PBL height
            call msg("Making Rdown from Kz and ABL only. Method:"+fu_str(Kz_method)+"ABL limit:",ABLh_min)
            !$OMP PARALLEL IF (.TRUE.) DEFAULT(NONE) &
            !$OMP & SHARED(Kz_method, height3d, fs_meteo, &
            !$OMP &    number_of_levels,ablh_min, &
            !$OMP &    u3d,v3d,q3d,t3d,error,Rptr, &
            !$OMP &    a_full,b_full,psrf,ptrhabl,ptrkz_1m,theta_prev, HuntenKZ) &
            !$OMP & PRIVATE(height1, height2,u_1,u_2,v_1,v_2,R_tmp ,i,z,Habl,Habl_01, lev,  pml,theta,ftmp,du,dv,dz,deltat,l,extrashear, &
            !$OMP &      t2,q2,deltau2,ri,phimphih,ithread, nthreads, iTmp, Kz_tab, zprev )
            ithread = 0
            nthreads = 1
            q2 => null()
            !$ ithread  = omp_get_thread_num()
            !$ nthreads = omp_get_num_threads()
            do lev=1, number_of_levels
              R_tmp => Rptr(lev)%pp
              zprev= Kz_ref_height
              height1 => height2 !Warning! uninitialized at first level
              height2 => fu_grid_data_from_3d(height3d, lev)
              do i=1+ithread,fs_meteo, nthreads
                 if (error) cycle
                 z = height2(i)
                 Habl = max(ptrHabl(i), ABLh_min)
                 Habl_01 = Habl * 0.1
                 if (lev > 1) zprev=height1(i)
                 if(z < Habl_01)then      !  within 10% of ABL
                   fTmp =  log(z/zPrev)
                 elseif(z < Habl) then
                    if (zprev < Habl_01) then
                       fTmp =  (log(Habl_01/zPrev) + (z - Habl_01) / Habl_01)
                    else
                       fTmp =  (z - zprev) / Habl_01
                    endif
                 else
                    if (zprev < Habl_01) then ! upper within ABL
                       fTmp = (log(Habl_01/zPrev) + 9. &
                                & + (z - Habl)/Habl_01  * Kz_FT_factor)
                    elseif(zprev < Habl)then
                       fTmp = ( Habl - zprev) / Habl_01 &
                               & + (z - Habl)/Habl_01  * Kz_FT_factor 
                    else
                       fTmp = (z - zprev)/Habl_01  * Kz_FT_factor 
                    endif
                 endif
                 R_tmp(i) = fTmp *  Kz_ref_height / ptrKz_1m(i)
                 if (Kz_method == hunten_kz) then
                    z  = 0.5*(height2(i) + zprev)
                    dz = height2(i) - zprev
                    iTmp = max(1,nint(z/1000.)) ! in km
                    if (iTmp <= KZtabLen) then
                       kZ_tab = HuntenKZ(iTmp) 
                    else
                       kZ_tab = HuntenKZ(KZtabLen) + &  !Extrapolate it
                          &(HuntenKZ(KZtabLen) - HuntenKZ(KZtabLen-1))*(iTmp - KZtabLen)
                    endif
                    fTmp = max(0., min(1.,  (z - 10000. )/5000.))  !Transition from 10lm to 15 km
                    R_tmp(i) = R_tmp(i) * (1. - fTmp) +  (dz/kZ_tab) * fTmp
                  endif


                 if(.not. r_tmp(i) > 0. .or. R_tmp(i) > 1e20 )then
                   call msg('negative R old Kz: level, i',lev,i)
                   call msg('z,zprev,Habl,Habl_01',(/z,zprev,Habl,Habl_01/))
                   call msg('level, fTmp, ptrKz_1m(i), Kz_ref_height, R_tmp(i)', &
                          & (/fTmp,ptrKz_1m(i),Kz_ref_height,R_tmp(i)/))
                   call set_error('Negative resistance','dq_Rto1m_profile')
                 endif
              end do ! loop meteo
            END DO ! loop_over_levels

            !$OMP END PARALLEL
            IF (error) RETURN 
          case default
            call msg("Unknown Kz_method:",Kz_method)
            call set_error("Can't make R_down_meteo_flag", "dq_Rto1m_profile")
            return
        end select


    END DO ! loop_over_times

!    call free_work_array(R_tmp)

  END SUBROUTINE dq_Rdown_profile

  !*****************************************************************************
  
  SUBROUTINE dq_turbulence_profile(meteoMarketPtr, met_src, valid_times, list)
    !
    IMPLICIT NONE
    !
    ! Imported parameters with intent(in):
    type(meteo_data_source), INTENT(in) :: met_src
    TYPE(silja_time), DIMENSION(:), INTENT(in) :: valid_times
    type(mini_market_of_stacks), pointer :: meteoMarketPtr
    TYPE(silja_shopping_list), INTENT(IN) :: list

    ! Imported parameters with intent(out):

    ! Local declarations:
    TYPE(silja_time) :: time
    INTEGER :: tloop, lev, iCell
    TYPE(silja_grid) :: scalar_grid
    TYPE(silja_field_id) :: id
    type(silja_field), pointer :: fldTmp
    TYPE(silja_3d_field), POINTER :: windspeed3d, height3d, Ri3d
    REAL, DIMENSION(:), POINTER :: Ez, Km,Kh, Kd, S, lz
    REAL, DIMENSION(:), POINTER :: height2, height1, height3, windspeed1, windspeed2, windspeed3, Ri, abl
    TYPE(silja_level), DIMENSION(max_levels) :: levels
    INTEGER :: number_of_levels
    real :: Rif, Psi, Psi_tau, fTmp
    logical :: ifTKE, ifKzMom, ifKzHeat, ifKzScal, ifS, ifTLS

    ifTKE = fu_quantity_in_list(turb_kinetic_energy_SILAM_flag, list)
    ifKzMom = fu_quantity_in_list(Kz_momentum_3d_flag, list)
    ifKzHeat = fu_quantity_in_list(Kz_heat_3d_flag, list)
    ifKzScal = fu_quantity_in_list(Kz_scalar_3d_flag, list)
    ifS = fu_quantity_in_list(wind_vertical_shear_flag, list)
    ifTLS = fu_quantity_in_list(turb_length_scale_flag, list)
    
    if(.not. ifS) S => fu_work_array() 
    if(.not. ifTLS) lz => fu_work_array() 
    if(.not. ifTKE) Ez => fu_work_array() 
    if(.not. ifKzMom) Km => fu_work_array() 
    if(.not. ifKzHeat) Kh => fu_work_array() 
    if(.not. ifKzScal) Kd => fu_work_array() 
    IF (error) RETURN

    DO tloop = 1, SIZE(valid_times)

      IF (.NOT.defined(valid_times(tloop))) EXIT   !loop_over_times
      time = valid_times(tloop)

      IF (fu_field_in_sm(meteoMarketPtr, met_src,&   ! Already done before?
                       & Kz_momentum_3d_flag,&
                       & time,&
                       & level_missing,&
                       & .false.,&
                       & .false.)) CYCLE  ! loop_over_times

      height3d => fu_sm_obstime_3d_field(meteoMarketPtr, met_src, height_flag, time, single_time)
      IF (error) RETURN

      Ri3d => fu_sm_obstime_3d_field(meteoMarketPtr, met_src, gradient_richardson_nbr_flag, &
                                      & time, single_time)
      IF (error) RETURN

      windspeed3d => fu_sm_obstime_3d_field(meteoMarketPtr, met_src, windspeed_flag, time, single_time)
      IF (error) RETURN

      CALL vertical_levels(windspeed3d, levels, number_of_levels)
      IF (error) RETURN

      !
      ! Get the screen-level fields for the lowest-level and
      ! make appropriate attributions to the heights
      !
      ! ABL
      !
      fldTmp => fu_sm_obstime_field(meteoMarketPtr, met_src, abl_height_m_pointer, level_missing, &
                                  & time, forwards)
      if(error)return
      abl => fu_grid_data(fldTmp)

      !
      ! Having the initial stuff set, start the layer cycle
      ! For vertical derivative, point 2 is always the height at which the derivative is taken, 
      ! point 1 is lower, point 3 is upper
      !
      DO lev = 1, number_of_levels

        Ri => fu_grid_data_from_3d(Ri3d, lev)

        if(ifS)then
          id = fu_set_field_id(fu_met_src(windspeed3d),&
                             & wind_vertical_shear_flag,&
                             & fu_analysis_time(windspeed3d),&
                             & fu_forecast_length(windspeed3d), &
                             & fu_grid(windspeed3d),&
                             & levels(lev))
          call find_field_data_storage_2d(meteoMarketPtr, id, multi_time_stack_flag, S)
          if(fu_fails(.not.error,'Failed S field data pointer','dq_turbulence_profile'))return
        endif

        if(ifTLS)then
          id = fu_set_field_id(fu_met_src(windspeed3d),&
                             & turb_length_scale_flag,&
                             & fu_analysis_time(windspeed3d),&
                             & fu_forecast_length(windspeed3d), &
                             & fu_grid(windspeed3d),&
                             & levels(lev))
          call find_field_data_storage_2d(meteoMarketPtr, id, multi_time_stack_flag, lz)
          if(fu_fails(.not.error,'Failed lz field data pointer','dq_turbulence_profile'))return
        endif
        if(ifTKE)then
          id = fu_set_field_id(fu_met_src(windspeed3d),&
                             & turb_kinetic_energy_SILAM_flag,&
                             & fu_analysis_time(windspeed3d),&
                             & fu_forecast_length(windspeed3d), &
                             & fu_grid(windspeed3d),&
                             & levels(lev))
          call find_field_data_storage_2d(meteoMarketPtr, id, multi_time_stack_flag, Ez)
          if(fu_fails(.not.error,'Failed Ez field data pointer','dq_turbulence_profile'))return
        endif
        if(ifKzMom)then
          id = fu_set_field_id(fu_met_src(windspeed3d),&
                             & Kz_momentum_3d_flag,&
                             & fu_analysis_time(windspeed3d),&
                             & fu_forecast_length(windspeed3d), &
                             & fu_grid(windspeed3d),&
                             & levels(lev))
          call find_field_data_storage_2d(meteoMarketPtr, id, multi_time_stack_flag, Km)
          if(fu_fails(.not.error,'Failed Km field data pointer','dq_turbulence_profile'))return
        endif
        if(ifKzHeat)then
          id = fu_set_field_id(fu_met_src(windspeed3d),&
                             & Kz_heat_3d_flag,&
                             & fu_analysis_time(windspeed3d),&
                             & fu_forecast_length(windspeed3d), &
                             & fu_grid(windspeed3d),&
                             & levels(lev))
          call find_field_data_storage_2d(meteoMarketPtr, id, multi_time_stack_flag, Kh)
          if(fu_fails(.not.error,'Failed Kh field data pointer','dq_turbulence_profile'))return
        endif
        if(ifKzScal)then
          id = fu_set_field_id(fu_met_src(windspeed3d),&
                             & Kz_scalar_3d_flag,&
                             & fu_analysis_time(windspeed3d),&
                             & fu_forecast_length(windspeed3d), &
                             & fu_grid(windspeed3d),&
                             & levels(lev))
          call find_field_data_storage_2d(meteoMarketPtr, id, multi_time_stack_flag, Kd)
          if(fu_fails(.not.error,'Failed Kd field data pointer','dq_turbulence_profile'))return
        endif
        
        if(lev == 1)then
          !
          ! Having height, select the data for determining the profiles: 10m wind 
          ! and 1st model level or just 1st and 2d model levels for everything.
          !
          height2 => fu_grid_data_from_3d(height3d, 1)

          if(sum(height2(1:fs_meteo))/fs_meteo <= 11.)then
            !
            ! First model level is dangerously close to 10m, skip the screen-level wind and use 3D only
            ! Accuracy is then 1st order
            !
            windspeed2 => fu_grid_data_from_3d(windspeed3d, 1)
            windspeed3 => fu_grid_data_from_3d(windspeed3d, 2)
            height3 => fu_grid_data_from_3d(height3d, 2)
            do iCell=1, fs_meteo
              S(iCell) = (windspeed3(iCell) - windspeed2(iCell)) / (height3(iCell) - height2(iCell))
            enddo
          else
            !
            ! Can use 10m and first model level to get the second order of accuracy
            !
            fldTmp => fu_sm_obstime_field(meteoMarketPtr, met_src, windspeed_10m_flag, level_missing, time, forwards)
            if(error)return
            windspeed1 => fu_grid_data(fldTmp)
            windspeed2 => fu_grid_data_from_3d(windspeed3d, 1)
            height2 => fu_grid_data_from_3d(height3d, 1)
            windspeed3 => fu_grid_data_from_3d(windspeed3d, 2)
            height3 => fu_grid_data_from_3d(height3d, 2)
            do iCell=1, fs_meteo
              S(iCell) = 1. / (10. - height3(iCell)) * &
                       & ((windspeed2(iCell) - windspeed3(iCell)) / (height2(iCell) - height3(iCell)) * &
                        & (10. - height2(iCell)) + &
                        & (windspeed1(iCell) - windspeed2(iCell)) / (10. - height2(iCell)) * &
                        & (height2(iCell) - height3(iCell)) &
                       & )
            end do
          endif

        elseif(lev == number_of_levels)then
          !
          ! Last model level - again 1st order of accuracy
          !
          windspeed1 => fu_grid_data_from_3d(windspeed3d, lev-1)
          height1 => fu_grid_data_from_3d(height3d, lev-1)
          windspeed2 => fu_grid_data_from_3d(windspeed3d, lev)
          height2 => fu_grid_data_from_3d(height3d, lev)
          do iCell=1, fs_meteo
            S(iCell) = (windspeed2(iCell) - windspeed1(iCell)) / (height2(iCell) - height1(iCell))
          end do

        else
          !
          ! General case: second-order for inhomogenous grid
          !
          windspeed1 => fu_grid_data_from_3d(windspeed3d, lev-1)
          height1 => fu_grid_data_from_3d(height3d, lev-1)
          windspeed2 => fu_grid_data_from_3d(windspeed3d, lev)
          height2 => fu_grid_data_from_3d(height3d, lev)
          windspeed3 => fu_grid_data_from_3d(windspeed3d, lev+1)
          height3 => fu_grid_data_from_3d(height3d, lev+1)
          do iCell=1, fs_meteo
            S(iCell) = 1. / (height1(iCell) - height3(iCell)) * &
                     & ((windspeed2(iCell) - windspeed3(iCell)) / (height2(iCell) - height3(iCell)) * &
                      & (height1(iCell) - height2(iCell)) + &
                      & (windspeed1(iCell) - windspeed2(iCell)) / (height1(iCell) - height2(iCell)) * &
                      & (height2(iCell) - height3(iCell)) &
                     & )
          end do

        endif  ! Level-dependent shear computation
        IF (error) RETURN

        !
        ! The rest of algebra for Ri, TKE and Km, Kh, Kd
        !
        do iCell=1, fs_meteo

          if(Ri(iCell) > 0)then
            Rif = Ri(iCell) * 1.25 * (1 + 36 * Ri(iCell))**1.7 / (1 + 19 * Ri(iCell))**2.7

            Psi_tau = 0.228 - 0.208 * Rif
            fTmp = 1 - 2.25 * Rif  !Psi3
            Psi =  0.55 * fTmp * Psi_tau * (1 - (1 + 1 / fTmp) * Rif)

            !
            ! A discussion of lz: therre has to be a large-scale limits for the size of eddies.
            ! Options are:
            ! 1) 0.1*Habl - in case of strongly convectional ABL is indeed a sort of measure
            !    But what impact does it have at e.g. tropopause?
            ! 2) A sum of squared reciprocals of several limiting factors: distance to the surface z
            !    rotation of Earth via Coriolis and wind shear
            !
            ! Option 1 would look like this:
            !
!            lz(iCell) = min(0.1 * abl(iCell), height2(iCell)) * (1 - Rif / 0.195)**(4. / 3.)

            ! Option 2 looks like this
            lz(iCell) = (1 - Rif / 0.195)**(4. / 3.) / sqrt( &
                                & 1. / (height2(iCell)*height2(iCell)) + &
                                & (basic_coriolis * basic_coriolis + S(iCell) * S(iCell)) / &
                                       & (0.01 * windspeed2(iCell) * windspeed2(iCell))  &
                                & )

            fTmp = S(iCell) * lz(iCell)
            Ez(iCell) = Psi * fTmp * fTmp



            Km(iCell) = 2 * Psi_tau * sqrt(Ez(iCell)) * lz(iCell)

            Kh(iCell) = Km(iCell) * Rif / Ri(iCell)

            fTmp = 1. / Psi
            Kd(iCell) =  Km(iCell)  * (1 - 0.3 * Rif * Psi_tau * fTmp) / ((1 + 0.3 * Ri(iCell) * fTmp) * Psi_tau)

          else
            Km(iCell) = 0. !real_missing
            Kh(iCell) = 0. !real_missing
            Kd(iCell) = 0. !real_missing
            Ez(iCell) = 0. !real_missing
            lz(iCell) = 0. !real_missing
          endif

        end do  ! cycle over the grid


!        call msg('')
!        call msg('level',lev)
!        call msg('Ri grad min', MINVAL(ri(1:fs_meteo)))
!        call msg('Ri grad max', MAXVAL(ri(1:fs_meteo)))
        do iCell=1, fs_meteo
          if(ri(iCell) > 1000)then
            call msg('Problem Ri at i=',iCell,ri(iCell))
          endif
        enddo

      

      END DO ! loop_over_levels

    END DO ! loop_over_times

    if(.not. ifS) call free_work_array(S)
    if(.not. ifTLS)call free_work_array(lz) 
    if(.not. ifTKE)call free_work_array(Ez) 
    if(.not. ifKzMom)call free_work_array(Km) 
    if(.not. ifKzHeat)call free_work_array(Kh) 
    if(.not. ifKzScal)call free_work_array(Kd) 

  END SUBROUTINE dq_turbulence_profile



  !*******************************************************************************
  subroutine dq_lai_dyn(meteoMarketPtr, met_src, valid_times, LAIsrc, ifRandomise)
     ! Just sums up leaf-area indices

    implicit none
    
    ! Imported parameters with intent(in):
    type(meteo_data_source), INTENT(in) :: met_src
    TYPE(silja_time), DIMENSION(:), INTENT(in) :: valid_times
    type(mini_market_of_stacks), pointer :: meteoMarketPtr
    integer, intent(in) :: LAIsrc
    logical, intent(in) :: ifRandomise

    ! Local declarations:
    TYPE(silja_time) :: time
    INTEGER :: tloop, iCell, ix, iy
    TYPE(silja_field_id) :: id
    type(silja_field), pointer :: fldTmp
    REAL, DIMENSION(:), POINTER :: LAI, LAIhv, LAIlv, hvFrac, lVfrac
    real :: lat,lon,cosZen, fTmp


    IF (error) RETURN

    DO tloop = 1, SIZE(valid_times)

      IF (.NOT.defined(valid_times(tloop))) EXIT   !loop_over_times
      time = valid_times(tloop)

      IF (fu_field_in_sm(meteoMarketPtr, met_src,&   ! Already done before?
                       & leaf_area_index_flag,&
                       & time,&
                       & surface_level,&
                       & .false.,&
                       & .false.)) CYCLE  ! loop_over_times
      !
      !
      if (LAIsrc == LAI_dynamic_2) then
             !Static fields
            fldTmp => fu_sm_obstime_field(meteoMarketPtr, met_src, fraction_hv_flag, &
                                        & level_missing,  time_missing, single_time)
            if(error)exit
            hvFrac => fu_grid_data(fldTmp)

            fldTmp => fu_sm_obstime_field(meteoMarketPtr, met_src, fraction_lv_flag, &
                                        & level_missing,  time_missing, single_time)
            if(error) exit
            lVfrac => fu_grid_data(fldTmp)

            !Dynamic fields
            fldTmp => fu_sm_obstime_field(meteoMarketPtr, met_src, leaf_area_indexhv_flag, &
                                        & level_missing,  time, single_time)
            if(error)exit
            LAIhv => fu_grid_data(fldTmp)

            fldTmp => fu_sm_obstime_field(meteoMarketPtr, met_src, leaf_area_indexlv_flag, &
                                        & level_missing,  time, single_time)
            if(error) exit
            LAIlv => fu_grid_data(fldTmp)
      else
            call msg("LAIsrc", LAIsrc)
            call set_error("Tried to make LAI form LAI... Should never happen", "dq_lai_dyn")
      endif



      id = fu_set_field_id(fu_met_src(fldTmp),&
                         & leaf_area_index_flag, &
                         & fu_analysis_time(fldTmp),&
                         & fu_forecast_length(fldTmp), &
                         & fu_grid(fldTmp),&
                         & surface_level)

      call find_field_data_storage_2d(meteoMarketPtr, id, multi_time_stack_flag, LAI)

      LAI(1:fs_meteo) = max(LAIhv(1:fs_meteo) * hvFrac(1:fs_meteo), 0.) &
                     &+ max(LAIlv(1:fs_meteo) * lvFrac(1:fs_meteo), 0.)

      !CALL dq_store_2d(meteoMarketPtr, id, lai, multi_time_stack_flag, ifRandomise)
      IF (error) exit

    end do  ! times
  end subroutine dq_lai_dyn
  
  
  !********************************************************************************
  
  subroutine dq_stomatal_conductance(meteoMarketPtr, met_src, valid_times, LAIsrc)
    !
    ! Computes the aerodynamic resistance Ra for dry deposition - from the given altitude down to 
    ! laminar layer
    !
    implicit none
    
    ! Imported parameters with intent(in):
    type(meteo_data_source), INTENT(in) :: met_src
    TYPE(silja_time), DIMENSION(:), INTENT(in) :: valid_times
    type(mini_market_of_stacks), pointer :: meteoMarketPtr
    integer, intent(in) :: LAIsrc

    ! Local declarations:
    TYPE(silja_time) :: time
    INTEGER :: tloop, iCell, ix, iy
    TYPE(silja_field_id) :: id
    type(silja_field), pointer :: fldTmp
    REAL, DIMENSION(:), POINTER :: h_canopy, LAI,  cloud, pressure, T2m, rh2m,  G_sto
    real :: lat,lon,cosZen

!    G_sto => fu_work_array() 
!    IF (error) RETURN

    DO tloop = 1, SIZE(valid_times)

      IF (.NOT.defined(valid_times(tloop))) EXIT   !loop_over_times
      time = valid_times(tloop)

      IF (fu_field_in_sm(meteoMarketPtr, met_src,&   ! Already done before?
                       & stomatal_conductance_flag,&
                       & time,&
                       & surface_level,&
                       & .false.,&
                       & .false.)) CYCLE  ! loop_over_times
      !
      !
      select case(LAIsrc)
         case (LAI_dynamic_1, LAI_dynamic_2)
            fldTmp => fu_sm_obstime_field(meteoMarketPtr, met_src, leaf_area_index_flag, &
                                        & level_missing, time, single_time)
         case (LAI_static_1, LAI_static_2)
            fldTmp => fu_sm_simple_field(meteoMarketPtr, met_src, leaf_area_index_flag, &
                                        & level_missing,  single_time_stack_flag)
         case default
            call msg("LAIsrc", LAIsrc)
            call set_error("Unknown LAIsrc", "dq_stomatal_conductance")
      end select
      if(error)return
      LAI => fu_grid_data(fldTmp)

      fldTmp => fu_sm_obstime_field(meteoMarketPtr, met_src, surface_pressure_flag, level_missing, & 
                                  & time, single_time)
      if(error)return
      pressure => fu_grid_data(fldTmp)

      fldTmp => fu_sm_obstime_field(meteoMarketPtr, met_src, total_cloud_cover_flag, level_missing, &  
                                  & time, single_time)
      if(error)return
      cloud => fu_grid_data(fldTmp)

      fldTmp => fu_sm_obstime_field(meteoMarketPtr, met_src, temperature_2m_flag, level_missing, &  
                                  & time, single_time)
      if(error)return
      t2m => fu_grid_data(fldTmp)

      fldTmp => fu_sm_obstime_field(meteoMarketPtr, met_src, relative_humidity_2m_flag, level_missing, &  
                                  & time, single_time)
      if(error)return
      rh2m => fu_grid_data(fldTmp)
      if(error)return
      !
      ! The main cycle itself
      !
      id = fu_set_field_id(fu_met_src(fldTmp),&
                         & stomatal_conductance_flag, &
                         & fu_analysis_time(fldTmp),&
                         & fu_forecast_length(fldTmp), &
                         & fu_grid(fldTmp),&
                         & surface_level)
      call find_field_data_storage_2d(meteoMarketPtr, id, multi_time_stack_flag, G_sto)
      if(fu_fails(.not.error,'Failed G_sto field data pointer','dq_stomatal_conductance'))return

      !call msg_warning('Three-fold lowered stomatal conductance','dq_stomatal_conductance')
!call msg('Full stomatal conductance','dq_stomatal_conductance')
      
      call SolarSetup(time)
      do iy = 1, ny_meteo
        do ix = 1, nx_meteo
          iCell = ix + (iy-1)* nx_meteo
          if (LAI(iCell) > 0.05) then
             lat = fu_lat_geographical_from_grid(real(ix), real(iy), meteo_grid)
             lon = fu_lon_geographical_from_grid(real(ix), real(iy), meteo_grid)
             cosZen = max(fu_solar_zenith_angle_cos(lon, lat, time),0.)
             G_sto(iCell) = fu_g_stomatal(LAI(iCell), T2m(iCell), rh2m(iCell), pressure(iCell), cloud(iCell), cosZEN)
!             G_sto(iCell) = 0.3 * fu_g_stomatal(LAI(iCell), T2m(iCell), rh2m(iCell), pressure(iCell), cloud(iCell), cosZEN)
          else
            G_sto(iCell) = 0.
          endif
        enddo
      enddo

    end do  ! times
    


  end subroutine dq_stomatal_conductance
  
  !************************************************************************

  real  function fu_g_stomatal(LAI, T2m, rh2m, p, cloud, cosZEN)
    !
    ! Get stomatal resistance 
    ! As done in EMEP model with some corner-cut...
        !!    Calculates stomatal conductance g_sto based upon methodology from 
        !!    EMEP MSC-W Note 6/00 and Mapping Manual (2004), and revisions (Simpson
        !!    and Emberson, Chapter 5, EMEP Rep 1/2006, Mapping Manual revisions, 2007,
        !!    and l. Emberson Notes from Forest group, Dec. 2007):
        !!
        !!    g_sto = [g_max * f_pot * f_light * f_temp * f_vpd * f_swp ]/41000.0
        !!
    ! Added nighttime g_sto 

    !
    real, intent(in) ::  cosZEN,  T2m, rh2m, p, cloud
    real, intent(in)  :: LAI       ! leaf area index (m^2/m^2), one-sided

    real :: Idirect, Idiffuse, Idrctn
    real :: LAIsun    ! sunlit LAI
    real :: ftmp, bt, dts, dg, att
    real :: vpd, vpsat, PARshade, PARsun, LAIsunfrac  
    real :: f_env,  f_vpd,  f_swp, f_phen, f_sun,   f_shade, f_light, f_temp
    real :: mmol2sm

    real,  parameter :: PARfrac = 0.45 ! approximation to fraction (0.45 to 0.5) of total 
                          ! radiation in PAR waveband (400-700nm)
    real,  parameter ::  Wm2_uE  = 4.57 ! converts from W/m^2 to umol/m^2/s
    real,  parameter ::  Wm2_2uEPAR= PARfrac * Wm2_uE ! converts from W/m^2 to umol/m^2/s PAR
    real,  parameter ::  do3se_f_light = 8e-3 ! LUC-dependent parameter in EMEP 0.005--0.013
                          ! except 0.002 for root_crop...
                          ! see Inputs_DO3SE.csv
    real,  parameter :: do3se_T_opt = 23 + 273   ! ! Some average... about med_needle
    real,  parameter :: do3se_T_min = 5  + 273   ! 
    real,  parameter :: do3se_T_max = 40 + 273   !

    !vapour pressure defficit
    real,  parameter :: do3se_VPD_min = 1. 
    real,  parameter :: do3se_VPD_max = 3. ! 2.5--3.5 
    real,  parameter :: do3se_f_min =  1e-1 
    real,  parameter :: do3se_f_min_lowveg =  1e-2
    real,  parameter :: do3se_g_max = 200  !!! mmol/m2/s


   !!
   !! ROUX: Nighttime g_sto does not get to zero! 
   ! Mairgareth A. Caird, James H. Richards, Lisa A. Donovan
   ! "Nighttime Stomatal Conductance and Transpiration in C3 and C4 Plant"
   !  DOI: https://doi.org/10.1104/pp.106.092940
   ! Fig1 suggests ~0.02-0.06 mol/m2/s for majority of species
   ! thus 20-60 mmol/m2/s
   !
    real,  parameter :: g_night = 40 !! mmol/m2/s
            !! 40 seems to be okay (VRA2016)

    real, parameter :: cosA    = 0.5   ! A = mean leaf inclination (60 deg.), 
     ! where it is assumed that leaf inclination has a spherical distribution


    ! Calculated by ClearSkyRadn in EMEP model
   if ( CosZEN > 1e-5 ) then
      att = 1.0 - 0.75*cloud**3.4  !(source: Kasten & Czeplak (1980)) 
      Idrctn     = att * Ashrae%a * exp(- Ashrae%b * (p/1e5)/CosZEN)
      Idiffuse   = Ashrae%c * Idrctn
      Idirect    = Idrctn *  CosZEN

      ! from subroutine CanopyPAR of EMEP model
      LAIsun = (1.0 - exp(-0.5*LAI/cosZEN) ) * cosZEN/cosA
      LAIsunfrac = LAIsun/LAI
  
  
      !! PAR flux densities evaluated using method of
      !! Norman (1982, p.79): 
      !! "conceptually, 0.07 represents a scattering coefficient"  
  
      PARshade = Idiffuse * exp(-0.5*LAI**0.7) +  &
                 0.07 * Idirect  * (1.1-0.1*LAI)*exp(-cosZEN)
  
      PARsun = Idirect *cosA/cosZEN + PARshade
  
      !.. Convert units, and to PAR fraction
  
      PARshade = PARshade * Wm2_2uEPAR
      PARsun   = PARsun   * Wm2_2uEPAR
!!---------------------------------------
!!    Calculate f_light, using methodology as described in Emberson et 
!!    al. (1998), eqns. 31-35, based upon sun/shade method of  
!!    Norman (1979,1982)

      f_sun   = (1.0 - exp (-do3se_f_light*PARsun  ) )
      f_shade = (1.0 - exp (-do3se_f_light*PARshade) )
      f_light = LAIsunfrac * f_sun + (1.0 - LAIsunfrac) * f_shade
    else
!      Idirect    = 0.0
!      Idiffuse   = 0.0
!      PARshade = 0.0
!      PARsun   = 0.0
!      f_sun  = 0.0
      f_light = 0.0

    end if
    f_light=max(g_night/do3se_g_max,f_light) 
            ! SILAM invention to leave some conductance 
            ! for night  
            ! For some vegetation can be up to 0.1
            

!!..3) Calculate  f_temp
!!---------------------------------------
!! Asymmetric  function from Mapping Manual
!! NB _ much more efficient to tabulate this - do later!

  dg  =    do3se_T_opt - do3se_T_min 
  bt  =    (do3se_T_max - do3se_T_opt ) / dg
  dTs = max( do3se_T_max - T2m, 0.0 )      !CHECK why max?
  f_temp = dTs / ( do3se_T_max - do3se_T_opt )
  f_temp = ( T2m - do3se_T_min ) / dg *  f_temp**bt
  f_temp = max( f_temp, 0.001 )  !Roux: Reduced  from 0.01 to 0.001 to supress the deposition in winter
  ! Revised usage of min value during 2007



!!..4) Calculate f_vpd ! Water pressure defficit
  !! Re-written from what is in EMEP


      ! Saturated pressure from Tetens equation https://en.wikipedia.org/wiki/Tetens_equation
  fTmp   = 17.67 * (T2m-273.15)/(T2m-29.65)
  vpSat = 611.2 * exp(fTmp) !in Pa
  vpd   = (1.0 - min(rh2m,0.95)) * vpSat  ! in Pa 
  
      ! Wang-VPD model, Journal of Hydrometeorology 13, 239 (Wang, 2012)
      !     https://doi.org/10.1175/JHM-D-11-043.1
  f_vpd = 51.5*vpd**(-.65)   !! Coefficient to match unity at +15

  !!!! double-counting for vpd improves O3 scores on summer evenings in Europe
  !! Suppresses deposition of O3
  !! I believe, HONO storage should take care of this in some future

!  f_vpd = f_vpd*f_vpd

!!..5) Calculate f_swp
!!---------------------------------------
!  !/  Use SWP_Mpa to get f_swp. We just need this updated
!  !   once per day, but for simplicity we do it every time-step.

!  !ds     f_swp = do3se(iLC)%f_min + &
!  !ds            (1-do3se(iLC)%f_min)*(do3se(iLC)%PWP-L%SWP)/ &
!  !ds                                (do3se(iLC)%PWP-do3se(iLC)%SWP_max)
!  !ds     f_swp = min(1.0,f_swp)
!  ! Aug 2010: use HIRLAM's SW, and simple "DAM" function
!  ! ************************************
!   if ( USE_SOILWATER ) f_swp =  L%fSW
!  ! ************************************       
   f_swp = 1.  !!! Disabled in original EMEP code
   f_phen = 1. 
!!.. And finally,
!!---------------------------------------
!!  ( with revised usage of min value for f_temp during 2007)

   f_env = f_vpd * f_swp
   f_env = max( f_env, do3se_f_min )
   f_env = f_temp * f_env

   f_env = f_phen * f_env * f_light  ! Canopy average

! From mmol O3/m2/s to s/m given in Jones, App. 3, gives 41000 for 20 deg.C )
!   (should we just use P=100 hPa?)

   mmol2sm = 8.3144e-8 * t2m       ! 0.001 * RT/P

   fu_g_stomatal = do3se_g_max * f_env * mmol2sm


   !g_sun = g_sto * f_sun/f_light       ! sunlit part

  endfunction fu_g_stomatal

  
  !*******************************************************************************

  subroutine dq_aerodynamic_resistance(meteoMarketPtr, met_src, valid_times, zRef)
    !
    ! Computes the aerodynamic resistance Ra for dry deposition - from the given altitude down to 
    ! laminar layer
    !
    implicit none
    
    ! Imported parameters with intent(in):
    type(meteo_data_source), INTENT(in) :: met_src
    TYPE(silja_time), DIMENSION(:), INTENT(in) :: valid_times
    type(mini_market_of_stacks), pointer :: meteoMarketPtr
    real, intent(in) :: zRef

    ! Local declarations:
    TYPE(silja_time) :: time
    INTEGER :: tloop, iCell
    TYPE(silja_field_id) :: id
    type(silja_field), pointer :: fldTmp
    REAL, DIMENSION(:), POINTER :: z0, L_inv, u_star, R_a

    IF (error) RETURN

    DO tloop = 1, SIZE(valid_times)

      IF (.NOT.defined(valid_times(tloop))) EXIT   !loop_over_times
      time = valid_times(tloop)

      IF (fu_field_in_sm(meteoMarketPtr, met_src,&   ! Already done before?
                       & R_a_flag,&
                       & time,&
                       & surface_level,&
                       & .false.,&
                       & .false.)) CYCLE  ! loop_over_times
      !
      ! Get the screen-level fields for the lowest-level and
      ! make appropriate attributions to the heights
      !
      fldTmp => fu_sm_obstime_field(meteoMarketPtr, met_src, friction_velocity_flag, level_missing, & ! u*
                                  & time, forwards)
      if(error)return
      u_star => fu_grid_data(fldTmp)

      fldTmp => fu_sm_obstime_field(meteoMarketPtr, met_src, surface_roughness_meteo_flag, level_missing, & ! z0
                                  & time, forwards)
      if(error)return
      z0 => fu_grid_data(fldTmp)

      fldTmp => fu_sm_obstime_field(meteoMarketPtr, met_src, MO_length_inv_flag, level_missing, &  ! L_inv
                                  & time, forwards)
      if(error)return
      L_inv => fu_grid_data(fldTmp)

      !
      ! The main cycle itself
      !
      id = fu_set_field_id(fu_met_src(fldTmp),&
                         & R_a_flag, &
                         & fu_analysis_time(fldTmp),&
                         & fu_forecast_length(fldTmp), &
                         & fu_grid(fldTmp),&
                         & surface_level)
      call find_field_data_storage_2d(meteoMarketPtr, id, multi_time_stack_flag, R_a)
      if(fu_fails(.not.error,'Failed R_a field data pointer','dq_aerodynamic_resistance'))return

      do iCell=1, fs_meteo
        R_a(iCell) = fu_R_a(zRef, z0(iCell), L_inv(iCell), u_star(iCell))
      enddo

!      CALL dq_store_2d(meteoMarketPtr, id, R_a, multi_time_stack_flag )
!      IF (error) RETURN 

    end do  ! times
    

  end subroutine dq_aerodynamic_resistance


  !**************************************************************************************************

  subroutine dq_hybrid_vertical_wind(meteoMarketPtr, met_src, valid_times)
    use vertical_motion
    implicit none

    ! Compute the vertical motion with respect to a hybrid vertical
    ! coordinate. 

    type(meteo_data_source), intent(in) :: met_src
    type(mini_market_of_stacks), intent(inout) :: meteoMarketPtr
    type(silja_time), dimension(:), intent(in) :: valid_times

    real, dimension(:), pointer :: u, v, massflux_u, massflux_v, dps_dt, dp_dt, &
         & div_x, div_y, div_integral, dp, eta_dot_2d
    
    type(silam_vertical), save :: meteo_layer_vertical
    logical :: first_time = .true.
    integer :: tloop, nlevs, fs, ilev
    type(silja_time) :: time
    real(vm_p), dimension(:,:), allocatable :: u3d, v3d, eta_dot_3d
    real(vm_p), dimension(:), allocatable :: ps_vm_p, spt_vm_p
    if (first_time) then
      call levels_to_layers(meteo_vertical, meteo_layer_vertical)
      !call report(meteo_layer_vertical, .true.)
      first_time = .false.
    end if
    
    fs = fs_meteo
     
    nlevs = fu_nbrOfLevels(meteo_layer_vertical)
      
    call init()

    do tloop = 1, size(valid_times)
      time = valid_times(tloop)
      if (.not. defined(time)) exit
      if (fu_field_in_sm(meteoMarketPtr, met_src, eta_dot_flag, &
                       & time, level_missing, &
                       & look_for_3d=.true., &
                       & permanent=.false.)) then
        call msg('Done already')
        cycle
      end if
      call process_valid_time(time)
    end do

    call cleanup()
    
  contains
    subroutine init()
      implicit none
      allocate(u3d(fs, nlevs), v3d(fs, nlevs), eta_dot_3d(fs, nlevs+1), ps_vm_p(fs), spt_vm_p(fs))
      dps_dt => fu_work_array()
!      eta_dot_2d => fu_work_array()
    end subroutine init
    
    subroutine cleanup
      call free_work_array(dps_dt)
!      call free_work_array(eta_dot_2d)
    end subroutine cleanup

    subroutine process_valid_time(time)
      use vertical_motion
      implicit none
      type(silja_time) :: time

      type(silja_3d_field), pointer :: u3d_fld, v3d_fld
      type(silja_field), pointer :: ps_fld
      real, dimension(:), pointer :: ps
      type(silja_level) :: lev_full
      type(silja_field_id) :: fid

      integer :: i

      u3d_fld => fu_sm_obstime_3d_field(meteoMarketPtr, met_src, u_flag, time, single_time)
      v3d_fld => fu_sm_obstime_3d_field(meteoMarketPtr, met_src, v_flag, time, single_time)
      ps_fld => fu_sm_obstime_field(meteoMarketPtr, met_src, ground_pressure_flag, &
                                  & level_missing, time, single_time)
      if (.not. associated(ps_fld)) then
        call set_error('Surface pressure field is not associated', 'process_valid_time')
        return
      end if
      ps => fu_grid_data(ps_fld)
      call get_tendency(ground_pressure_flag, level_missing, time, dps_dt)
      
      call start_count('process_valid_time_in_out')
      do ilev = 1, nlevs
        u => fu_grid_data_from_3d(u3d_fld, ilev)
        v => fu_grid_data_from_3d(v3d_fld, ilev)
        u3d(:,ilev) = u
        v3d(:,ilev) = v
      end do
      do i = 1, fs
        ps_vm_p(i) = real(ps(i), vm_p)
        spt_vm_p(i) = real(dps_dt(i), vm_p)
      end do
      call stop_count('process_valid_time_in_out')

      call diagnose_eta_dot(u3d, v3d, ps_vm_p, spt_vm_p, &
                          & meteo_layer_vertical, &
                          & fu_grid(u3d_fld), fu_grid(v3d_fld), meteo_grid, &
                          & .false., .false., eta_dot_3d) ! half levels, no top-down

      call start_count('process_valid_time_in_out')

      do ilev = 1, nlevs
        lev_full = fu_level(meteo_vertical, ilev)
        fid = fu_set_field_id(met_src, eta_dot_flag, &
                            & fu_analysis_time(u3d_fld), fu_forecast_length(u3d_fld), &
                            & meteo_grid, lev_full)
        call find_field_data_storage_2d(meteoMarketPtr, fid, multi_time_stack_flag, eta_dot_2d)
        if(fu_fails(.not.error,'Failed eta_dot_2d field data pointer','process_valid_time'))return
        
        eta_dot_2d(1:fs) = 0.5 * (eta_dot_3d(1:fs,ilev) + eta_dot_3d(1:fs,ilev+1))
!        call dq_store_2d(meteoMarketPtr, fid, eta_dot_2d, multi_time_stack_flag)
      end do
      call stop_count('process_valid_time_in_out')

    end subroutine process_valid_time

    subroutine get_tendency(quantity, level, time, tendency)
      implicit none
      integer, intent(in) :: quantity
      type(silja_level) :: level
      type(silja_time) :: time
      real, dimension(:), intent(out) :: tendency
      
      real, dimension(:), pointer :: past, future
      real :: dt
      logical :: is_first_time
      type(silja_field), pointer :: past_field, future_field
      integer :: fs

      fs = fs_meteo

      is_first_time = time == valid_times(1)
      
      past => fu_work_array()
      future => fu_work_array()

      if (is_first_time) then
        past_field => fu_sm_obstime_field(meteoMarketPtr, met_src, quantity, &
                                        & level, time, backwards)
      else
        past_field => fu_sm_obstime_field(meteoMarketPtr, met_src, quantity, &
                                        & level, time-ten_seconds, backwards)
        if (error) then
          call msg('Why?')
          return
        end if
        
      end if
      
      if (fu_valid_time(past_field) == time) then
        future_field => fu_sm_obstime_field(meteoMarketPtr, met_src, quantity, &
                                         & level, time+ten_seconds,  forwards)
      else
        future_field => fu_sm_obstime_field(meteoMarketPtr, met_src, quantity, &
                                         & level, time,  single_time)
      end if

      if (error) return
      
      past => fu_grid_data(past_field)
      future => fu_grid_data(future_field)
      dt = fu_sec(fu_valid_time(future_field)-fu_valid_time(past_field))
      tendency(1:fs) = (future(1:fs) - past(1:fs)) / dt
      
      call free_work_array(past)
      call free_work_array(future)

    end subroutine get_tendency
    
  end subroutine dq_hybrid_vertical_wind

  
  !**********************************************************************************************

  subroutine dq_ref_evapotranspiration(meteoMarketPtr, met_src, valid_times)
    !
    ! Computes the referece evapotranspiration ET0 following FAO-1998 standard model
    ! Allen,R.G., Pereira,L.S., Raes,D., Smith,M. (1998) Crop evapotranspiration –
    ! Guidelines for computing crop water requirements. FAO Irrigation and drainage paper 56.
    ! Food and Agriculture Organization of the United Nations, Rome
    ! Reference ET0 is evapotranspiration from a reference surface. The reference surface is
    ! 0.12cm dense grass, fixed surface resistance of 70 s/m and albedo 0.23. The reference surface
    ! is well watered, green, actively growing and completely shadowing the underlying soil.
    !
    implicit none
    
    ! Imported prameters
    type(meteo_data_source), intent(in) :: met_src
    type(mini_market_of_stacks), intent(inout) :: meteoMarketPtr
    type(silja_time), dimension(:), intent(in) :: valid_times
    
    ! Local parameters
    real, parameter :: r_bulk_stom_one_leaf = 100.0         ! sec/m, a typical value

            
    ! Local declarations:
    TYPE(silja_time) :: time
    INTEGER :: tloop, iCell
    TYPE(silja_field_id) :: id
    type(silja_field), pointer :: fldTmp
    REAL, DIMENSION(:), POINTER :: T, net_sw_rad_srf, net_lw_rad_srf,LAI, Psrf, rh, r_a, ET0  !, air_density
    real :: net_radiation, T_C, r_s, dPsw_dT, psychrometric_const, saturWaterVapourPr

    IF (error) RETURN

    DO tloop = 1, SIZE(valid_times)

      IF (.NOT.defined(valid_times(tloop))) EXIT   !loop_over_times
      time = valid_times(tloop)

      IF (fu_field_in_sm(meteoMarketPtr, met_src,&   ! Already done before?
                       & ref_evapotranspiration_flag,&
                       & time,&
                       & surface_level,&
                       & .false.,&
                       & .false.)) CYCLE  ! loop_over_times
      !
      ! Get the screen-level fields for the lowest-level and
      ! make appropriate attributions to the heights
      !
      fldTmp => fu_sm_obstime_field(meteoMarketPtr, met_src, temperature_2m_flag, level_missing, time, forwards)
      if(error)return
      T => fu_grid_data(fldTmp)

      fldTmp => fu_sm_obstime_field(meteoMarketPtr, met_src, surf_sw_net_radiation_flag, level_missing, time, forwards)
      if(error)return
      net_sw_rad_srf => fu_grid_data(fldTmp)
      
      fldTmp => fu_sm_obstime_field(meteoMarketPtr, met_src, surf_lw_net_radiation_flag, level_missing, time, forwards)
      if(error)return
      net_lw_rad_srf => fu_grid_data(fldTmp)
      
      fldTmp => fu_sm_obstime_field(meteoMarketPtr, met_src, leaf_area_index_flag, level_missing, time, forwards)
      if(error)return
      LAI => fu_grid_data(fldTmp)
      
      fldTmp => fu_sm_obstime_field(meteoMarketPtr, met_src, surface_pressure_flag, level_missing, time, forwards)
      if(error)return
      Psrf => fu_grid_data(fldTmp)
      
!      fldTmp => fu_sm_obstime_field(meteoMarketPtr, met_src, air_density_flag, level_missing, time, forwards)
!      if(error)return
!      air_density => fu_grid_data(fldTmp)
      
      fldTmp => fu_sm_obstime_field(meteoMarketPtr, met_src, relative_humidity_flag, level_missing, time, forwards)
      if(error)return
      rh => fu_grid_data(fldTmp)
      
      fldTmp => fu_sm_obstime_field(meteoMarketPtr, met_src, r_a_flag, level_missing, time, forwards)
      if(error)return
      r_a => fu_grid_data(fldTmp)
      !
      ! The main cycle itself
      !
      id = fu_set_field_id(fu_met_src(fldTmp),&
                         & ref_evapotranspiration_flag, &
                         & fu_analysis_time(fldTmp),&
                         & fu_forecast_length(fldTmp), &
                         & fu_grid(fldTmp),&
                         & surface_level)
      call find_field_data_storage_2d(meteoMarketPtr, id, multi_time_stack_flag, ET0)
      if(fu_fails(.not.error,'Failed ET0 field data pointer','dq_ref_evapotranspiration'))return
      !
      ! The main cycle over grid
      !
      do iCell = 1, fs_meteo
        net_radiation = net_sw_rad_srf(iCell) + net_lw_rad_srf(iCell)
        if(net_radiation > 0.)then   ! daytime => plants breath but soil takes roughly 10% of energy for heating-up
          net_radiation = 0.9 * net_radiation
        else
          ET0(iCell) = 0.0
          cycle
        endif
        T_C = T(iCell) - 273.15
        saturWaterVapourPr = 610.78 * exp(17.2694 * T_C / (T_C + 238.3))      ! Pa
        dPsw_dT = 4098. * saturWaterVapourPr / ((T_C + 237.3)*(T_C + 237.3))  ! Pa/K, just derivative
        r_s = r_bulk_stom_one_leaf * 2. / LAI(iCell)                  ! only upper leaf side => LAI/2.
        psychrometric_const = specific_heat_dryair * Psrf(iCell) / (condensation_latentheat * 0.622) ! 0.622=18g/mole / 29g/mole
          
!        ET0(iCell) = (dPsw_dT * net_radiation + air_density(iCell) * specific_heat_dryair * &
        ET0(iCell) = (dPsw_dT * net_radiation + density_air * specific_heat_dryair * &
                                          & saturWaterVapourPr * (1. - rh(iCell)) / r_a(iCell)) / &
                   & (dPsw_dT + psychrometric_const * (1. + r_s / r_a(iCell)))

      end do! fs_meteo
      
    end do   ! tloop

    end subroutine dq_ref_evapotranspiration
  

!
!subroutine dq_rain_intensity(meteoMarketPtr, & 
!                &field,  validtime,  desiredQ)
!  type(Tfield_buffer), pointer :: met_buf              ! input meteo data: buffer!
!  type(mini_market_of_stacks), pointer :: dispMarketPtr   ! working mini-market
!  type(silja_field), pointer :: field
!  type(silja_time), intent(in) :: validtime
!  integer, intent(in) :: desiredQ                 ! input and output quantities wanted
!  integer :: accQ, indAcc, ix, iy
!  real, dimension(:), pointer :: instant_rain
!  real :: acc_next, acc_prev, aclen
!  logical :: ifHorizInterp, previous_same_start
!  type(THorizInterpStruct), pointer :: pHorizInterpStruct ! meteo 2 dispersion horizontal
!  type(silja_field_id), pointer :: idFuture, idPast, idPtr 
!  type(silja_field_id) :: idTmp
!  
!
!    select case (desiredQ)
!      case (large_scale_rain_int_flag)
!              accQ =  large_scale_accum_rain_flag
!      case (convective_rain_int_flag)
!              accQ = convective_accum_rain_flag
!      case (total_precipitation_rate_flag)
!              accQ = total_precipitation_acc_flag
!      case default
!              call set_error("Can't make rain out of "+ fu_quantity_string(desiredQ), &
!                      "df_make_rain_intensity")
!              return
!    end select
!    
!    ifHorizInterp = .not. (meteo_grid == dispersion_grid)
!    if(ifHorizInterp) &
!       & pHorizInterpStruct => fu_horiz_interp_struct(meteo_grid, dispersion_grid, linear)
!
!!    call msg("Filling field with data...")
!!    call msg(fu_quantity_string(desiredQ))
!
!    instant_rain => fu_grid_data(field)
!    indAcc = fu_index(met_buf%buffer_quantities, accQ)
!
!    idFuture => met_buf%p2d(indAcc)%future%idPtr
!    idPast  => met_buf%p2d(indAcc)%past%idPtr
!
!
!    previous_same_start = fu_accumulation_start_time(idFuture) == &
!                            & fu_accumulation_start_time(idPast)
!
!    if(previous_same_start) then 
!        ! Previous field has same start of accumulation.
!        aclen = fu_sec(fu_accumulation_length(idFuture) - fu_accumulation_length(idPast))
!
!        IF (aclen > 0.) THEN
!           if (ifHorizInterp)  then
!              do iy = 1, ny_dispersion
!                do ix = 1, nx_dispersion
!                    acc_next = fu_get_value(met_buf%p2d(indAcc), &
!                               & nx_meteo, ix, iy, &
!                               & 0.0, &
!                               & pHorizInterpStruct, ifHorizInterp, .true.)
!                    acc_prev = fu_get_value(met_buf%p2d(indAcc), &
!                               & nx_meteo, ix, iy, &
!                               & 1.0, &
!                               & pHorizInterpStruct, ifHorizInterp, .true.)
!                    instant_rain(ix + (iy-1)*nx_dispersion) =  &
!                               & (acc_next - acc_prev) / aclen
!               enddo !ix
!              enddo !iy
!           else ! no interpolation
!              instant_rain(1:fs_dispersion) = &
!                      & ( met_buf%p2d(indAcc)%future%ptr(1:fs_meteo) &
!                      &  - met_buf%p2d(indAcc)%past%ptr(1:fs_meteo) ) / aclen
!           endif
!
!        ELSE
!          instant_rain(1:fs_dispersion) = 0.
!          CALL report(met_buf%p2d(indAcc)%future%idPtr)
!          CALL report(met_buf%p2d(indAcc)%past%idPtr)
!          CALL set_error('zero accumulation between obstimes',&
!              & 'dq_average_between_obstimes')
!          RETURN
!        END IF
!
!      ELSE ! No previous field exists, or following fields
!
!        aclen = fu_sec(fu_accumulation_length(idFuture))
!
!        IF (aclen > 0.) THEN
!           if (ifHorizInterp)  then
!              do iy = 1, ny_dispersion
!                do ix = 1, nx_dispersion
!                    acc_next = fu_get_value(met_buf%p2d(indAcc), &
!                               & nx_meteo, ix, iy, &
!                               & 0.0, &
!                               & pHorizInterpStruct, ifHorizInterp, .true.)
!                    instant_rain(ix + (iy-1)*nx_dispersion) =  acc_next  / aclen
!               enddo !ix
!              enddo !iy
!           else ! no interpolation
!              instant_rain(1:fs_dispersion) =     &
!                      & met_buf%p2d(indAcc)%future%ptr(1:fs_meteo) /aclen
!           endif
!        ELSE
!          instant_rain(1:fs_dispersion) = 0.
!          CALL report(met_buf%p2d(indAcc)%future%idPtr)
!          CALL msg_warning('zero accumulation', 'dq_average_between_obstimes')
!        END IF
!
!      END IF  ! previous exist
!
!
!      !
!      ! Almost done. Now set proper validity in field ID
!      !
!      idPtr => fu_id(field)
!      call set_valid_time(idPtr, fu_valid_time(idPast))
!      call set_analysis_time(idPtr, fu_analysis_time(idPast))
!      call set_validity_length(idPtr, fu_set_interval_sec(aclen))
!
!end subroutine dq_rain_intensity
!
!

    
      !****************************************************************

  SUBROUTINE dq_cloudcvr_from_3d(meteoMarketPtr, met_src, valid_times)
    !
    ! Make total cloud cover from 3d
    ! Take max over the levels
    !
    IMPLICIT NONE
    !
    ! Imported parameters with intent(in):
    type(meteo_data_source), INTENT(in) :: met_src
    TYPE(silja_time), DIMENSION(:), INTENT(in) :: valid_times
    type(mini_market_of_stacks), pointer :: meteoMarketPtr

    ! Local declarations
    REAL, DIMENSION(:), POINTER :: cc_lev, cc_tot 
    INTEGER :: tloop, i, nLevs, iLev
    TYPE(silja_3d_field), POINTER :: cc3d
    TYPE(silja_level), DIMENSION(max_levels) :: levels
    TYPE(silja_time) :: time
    TYPE(silja_field_id) :: id

    loop_over_times: DO tloop = 1, SIZE(valid_times)
      IF (.NOT.defined(valid_times(tloop))) EXIT loop_over_times
      time = valid_times(tloop)

      IF (fu_field_in_sm(meteoMarketPtr, met_src, & ! Already done before?
                       & total_cloud_cover_flag,&
                       & time,&
                       & level_missing,&
                       & .false.,&  ! 3D
                       & .false.)) CYCLE loop_over_times

      cc3d => fu_sm_obstime_3d_field(meteoMarketPtr, met_src,&         
                                  & cloud_cover_flag, &  
                                  & time,&
                                  & single_time)
      if(error)cycle loop_over_times
      
      CALL vertical_levels(cc3d, levels, nLevs)
      IF (error) CYCLE loop_over_times
      
      id = fu_set_field_id(fu_met_src(cc3d),&
                         & total_cloud_cover_flag,&
                         & fu_analysis_time(cc3d),&
                         & fu_forecast_length(cc3d), &
                         & fu_grid(cc3d),&
                         & ground_level)
      call find_field_data_storage_2d(meteoMarketPtr, id, multi_time_stack_flag, cc_tot)
      if(fu_fails(.not.error,'Failed cc_tot field data pointer','dq_cloudcvr_from_3d'))return

      cc_tot(1:fs_meteo) = 0.0
      do iLev = 1, nLevs
        cc_lev => fu_grid_data_from_3d(cc3d,iLev)  
        if(error)cycle loop_over_times
        do i = 1, fs_meteo
          cc_tot(i) = max(cc_tot(i), cc_lev(i))
        end do
      enddo
    END DO loop_over_times
  END SUBROUTINE dq_cloudcvr_from_3d
    
    
    
END MODULE derived_field_quantities_2

