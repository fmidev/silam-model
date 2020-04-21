module enkf
  ! 
  ! Ensemble Kalman Filter
  !
  ! Below is implementation of EnKF analysis algorithm adapted from the TOPAZ assimilation
  ! system of P. Sakov et al. The analysis step is rearranged to run on single pass,
  ! in-memory, and linear algebra has been changed to use BLAS subroutines.
  !
  ! Two versions of EnKF are available. With enkf_flavor == enkf_enkf, the classical EnKF
  ! with perturbed observations is used. With enkf_flavor == enkf_denkf, the
  ! "deterministic EnKF" [1] is used.
  ! 
  ! Localisation is performed unless turned off with loc_type == loc_none. Either step or
  ! Gaspari-Cohn localisation functions are available.
  !
  ! References: 
  ! [1] Sakov, P., Oke, P.R., 2008. A deterministic formulation of the
  ! ensemble Kalman filter: an alternative to ensemble square root filters. Tellus A 60,
  ! 361–371. doi:10.1111/j.1600-0870.2007.00299.x

!  use toolbox
  use geography_tools  !, only : fu_gc_distance
  use optimisation, only : start_count, stop_count
  implicit none

  private
  
  public enkf_update
  public enkf_update_openmp

  public get_row

  integer, parameter, public :: loc_none = 0, loc_step = 1, loc_gaspari_cohn = 2
  integer, parameter, public :: enkf_enkf = 0, enkf_denkf = 1
  real, private, parameter :: degrad = pi / 180.0

contains
  
  subroutine enkf_update(ens, ens_obs, obs_data, obs_var, obs_loc, mdl_loc, &
                       & loc_dist_m, loc_type, enkf_flavor, rfactor, diagn_out_dir)
    implicit none
    real, dimension(:,:) :: ens     ! E, (ind_mdl, ind_ens)
    real, dimension(:,:) :: ens_obs ! HE (ind_obs, ind_ens)
    real, dimension(:) :: obs_data  ! d
    real, dimension(:) :: obs_var   ! diag(R)
    ! localisation of model and obs points
    real, dimension(:,:) :: obs_loc ! (lon/lat,ind_obs)
    real, dimension(:,:) :: mdl_loc ! (lon/lat, ind_mdl)
    real, intent(in) :: loc_dist_m 
    integer, intent(in) :: loc_type ! see above
    integer, intent(in) :: enkf_flavor ! see above
    real, intent(in) :: rfactor ! inflation factor, multiplies R for updating the ensemble anomalies.
    ! directory for some diagnostic output, char_missing -> no output made
    character(len=*), intent(in) :: diagn_out_dir 

    character(len=*), parameter :: sub_name = 'run_enkf'
    real, dimension(:), pointer :: obs_ens_mean, stdev_inv, innov_mean, tmp
    integer :: ens_size, ii, obs_size, mdl_size
    integer :: stat, ind_obs, ind_ens, ind_mdl, pseudoinv_count
    real, dimension(:,:), allocatable :: obs_anom, &  ! S, (ind_end, ind_obs) &
         & subS, &! local subset of S
         & X1, &  ! either SS^T + I or S^TS + I depending on which is smaller
         & Gmatr, Gtransp,&   ! gain matrix 
         & X5tmp, &     ! local subset of X5
         & obs_pert, &     ! observation perturbation in ENKF
         & subD            ! local subset of obs_pert
    real, dimension(:), allocatable :: subdy, & ! local subset of innov_mean
         & lfactors, & ! localization factors
         & obs_dist    ! distance to the local observations
    integer, dimension(:), allocatable :: obs_ind ! indices of locally selected observations
    real :: sqrtminv, lon_prev, lat_prev
    integer :: iter, loc_obs_found, num_loc_obs, skip_count, skip_count_obs
    logical :: use_obs_space
    integer, dimension(:), pointer :: loc_obs_used
    real, dimension(size(ens,2)) :: grpt_ens

    !call set_error('Crash!', sub_name)
    !return

    call start_count('enkf_prepare')
    ens_size = size(ens, 2)
    mdl_size = size(ens, 1)
    obs_size = size(obs_data)
    if (fu_fails(ens_size > 1, 'ens_size < 2', sub_name)) return
    if (fu_fails(size(obs_var) == obs_size, 'obs_var and obs_data don''tt conform', sub_name)) return
    if (.not. fu_index(loc_type, (/loc_none, loc_step, loc_gaspari_cohn/)) > 0) then 
      call set_error('Bad loc_type', sub_name)
    end if
    if (.not. (enkf_flavor == enkf_enkf .or. enkf_flavor == enkf_denkf)) then
      call set_error('Bad enkf_flavor', sub_name)
    end if
    if (.not. loc_dist_m > 0.0 .and. loc_dist_m < 2e7) then
      call msg('loc_dist_m:', loc_dist_m)
      call set_error('Strange loc_dist_m', sub_name)
    end if
    if (fu_fails(rfactor > 0, 'Strange rfactor', sub_name)) return
    if (error) return

    call msg('ENKF: Dimension of observation vector:', obs_size)
    call msg('ENKF: Dimension of state vector:', mdl_size)

    if (fu_fails(obs_size > 0, 'No observations given', sub_name)) return
     
    allocate(obs_anom(obs_size, ens_size), obs_ind(obs_size), obs_dist(obs_size), &
           & X5tmp(ens_size, ens_size), &
           & stat=stat)
    if (fu_fails(stat == 0, 'Allocate failed', sub_name)) return

    ! calculate S = R^(-1/2) . (HE - mean(HE, dim=ens)) / sqrt(m-1), m = ensemble size
    sqrtminv = 1/sqrt(real(ens_size) - 1)
    stdev_inv => fu_work_array(obs_size)
    stdev_inv(1:obs_size) = 1 / sqrt(obs_var(1:obs_size))
    obs_ens_mean => fu_work_array(obs_size)
    if (error) return
    do ind_obs = 1, obs_size
      obs_ens_mean(ind_obs) = sum(ens_obs(ind_obs, :)) / ens_size
      !print *, 'OEM', obs_ens_mean(ind_obs)
      do ind_ens = 1, ens_size
        obs_anom(ind_obs, ind_ens) &
             & = sqrtminv * stdev_inv(ind_obs) * (ens_obs(ind_obs, ind_ens) - obs_ens_mean(ind_obs)) 
      end do
    end do
    
    innov_mean => fu_work_array(obs_size)
    if (error) return
    ! calculate s = R^(-1/2) . (d - Hx) / sqrt(m-1), m = ensemble size, 
    ! d = data = obs, s is d or dy in TOPAZ code
    do ind_obs = 1, obs_size
      innov_mean(ind_obs) = sqrtminv * stdev_inv(ind_obs) * (obs_data(ind_obs) - obs_ens_mean(ind_obs))
    end do
    
    if (enkf_flavor == enkf_enkf) then
      allocate(obs_pert(obs_size,ens_size), stat=stat)
      if (fu_fails(stat == 0, 'Allocate failed', sub_name)) return
      do ind_obs = 1, obs_size
        call random_normal(obs_pert(ind_obs,:))
        ! topaz code has here rfactor without sqrt ??
        obs_pert(ind_obs,:) = obs_pert(ind_obs,:) * sqrtminv / sqrt(rfactor)
      end do
    end if
    
    ! preparation done.
    call stop_count('enkf_prepare')
    ! pre-allocate local arrays for 100 obs (enlarge if needed)
    call allocate_loc_arr(2000, ens_size, subS, subdy, lfactors, Gmatr, Gtransp, subD, X1)
    !tmp => fu_work_array()

    ! go cell by cell, calculate local X5 and update
    loc_obs_found = 0
    loc_obs_used => fu_work_int_array(mdl_size)
    loc_obs_used(1:mdl_size) = 0
    lon_prev = real_missing
    lat_prev = real_missing
    call start_count('enkf_update')
    
    ! ****
    ! Loop over gridpoints
    pseudoinv_count = 0
    skip_count = 0
    skip_count_obs = 0
    do ind_mdl = 1, mdl_size
      if (fu_str(int(ind_mdl/real(mdl_size)*10)) /= fu_str(int((ind_mdl-1)/real(mdl_size)*10))) then
        call msg('ENKF: gridpoints done, total', ind_mdl, mdl_size)
      end if
      grpt_ens = ens(ind_mdl,:)
      
      ! If the spread is very small (possibly because all values are near-zero), skip the analysis:
      if (maxval(grpt_ens) - minval(grpt_ens) < max(1e-12, 1e-6*sum(grpt_ens)/ens_size)) then
        skip_count = skip_count + 1
        cycle
      end if
      
      ! If new horizontal gridpoint, get new local observations and ensemble transformation X5.
      !
      if (.not. ((mdl_loc(1,ind_mdl) .eps. lon_prev) .and. (mdl_loc(2,ind_mdl) .eps. lat_prev))) then
        lon_prev = mdl_loc(1,ind_mdl)
        lat_prev = mdl_loc(2,ind_mdl)
        call get_local_obs(loc_type, mdl_loc(1,ind_mdl), mdl_loc(2,ind_mdl), obs_loc, &
                         & loc_dist_m, obs_ind, obs_dist, num_loc_obs)
        if (num_loc_obs == 0) then
          skip_count_obs = skip_count_obs + 1
          cycle
        end if
                
        loc_obs_used(ind_mdl) = num_loc_obs
        use_obs_space = num_loc_obs < ens_size
        !use_obs_space = .false.
        !use_obs_space = .true.
        if (size(subdy) < num_loc_obs) then 
        !if (.true.) then 
          deallocate(subS, subdy, lfactors, Gmatr, Gtransp, X1, subD)
          loc_obs_found = num_loc_obs
        end if
        if (.not. allocated(subS)) then
          call allocate_loc_arr(num_loc_obs, ens_size, subS, subdy, lfactors, Gmatr, Gtransp, subD, X1)
        end if
        if (error) return
        ! Extract local sub-matrices of S, dy
        subS(1:num_loc_obs,:) = obs_anom(obs_ind(1:num_loc_obs),:)
        subdy(1:num_loc_obs) = innov_mean(obs_ind(1:num_loc_obs))
        if (enkf_flavor == enkf_enkf) subD(1:num_loc_obs,:) = obs_pert(obs_ind(1:num_loc_obs), :)

        call apply_locfun(loc_type, enkf_flavor, loc_dist_m, subS, subdy, subD, obs_dist, num_loc_obs, lfactors)
        if (error) return
        
        call get_x5_blas(subS, X1, Gtransp, Gmatr, subdy, X5tmp, subD, rfactor, ens_size, num_loc_obs, &
                       & pseudoinv_count, enkf_flavor)
        !call get_x5_f90()
      else if(num_loc_obs == 0) then
        cycle
      end if
      ! now have local X5 -> update the ensemble
      ! y = alpha*Ax + beta*y
      ! SUBROUTINE SGEMV ( TRANS, M, N, ALPHA, A, LDA, X, INCX, BETA, Y, INCY )
      call sgemv('t', ens_size, ens_size, 1.0, X5tmp, ens_size, grpt_ens, 1, 0.0, ens(ind_mdl,1), mdl_size)
    
    end do ! over state vector
    call stop_count('enkf_update')
    call msg('ENKF: number of pseudoinverses:', pseudoinv_count)
    
    ! Note: if the filter is used with domain decomposition, most or all gridpoints might end up skipped!
    call msg('ENKF: number of local analyses skipped:', skip_count)
    call msg('ENKF: number of local analyses skipped due to no observations:', skip_count_obs)
    
    if (diagn_out_dir /= char_missing) then
      ! Write the number of observations used for each analysis point.
      call dump_int_arr(loc_obs_used(1:mdl_size), diagn_out_dir + dir_slash + 'loc_obs.dat')
    end if

    call free_work_array(loc_obs_used)
    call free_work_array(innov_mean)
    call free_work_array(stdev_inv)
    call free_work_array(obs_ens_mean)

  end subroutine enkf_update

  !************************************************************************************

  subroutine enkf_update_openmp(ens, ens_obs, obs_data, obs_var, obs_loc, mdl_loc, &
                              & loc_dist_m, loc_type, enkf_flavor, rfactor, diagn_out_dir)
    ! 
    ! OpenMP version of enkf_update. For parallelisation details, see below.
    ! 
    implicit none
    real, dimension(:,:) :: ens     ! E, (ind_mdl, ind_ens)
    real, dimension(:,:) :: ens_obs ! HE (ind_obs, ind_ens)
    real, dimension(:) :: obs_data  ! d
    real, dimension(:) :: obs_var   ! diag(R)
    ! localisation of model and obs points
    real, dimension(:,:) :: obs_loc ! (lon/lat,ind_obs)
    real, dimension(:,:) :: mdl_loc ! (lon/lat, ind_mdl)
    real, intent(in) :: loc_dist_m 
    integer, intent(in) :: loc_type ! see above
    integer, intent(in) :: enkf_flavor ! see above
    real, intent(in) :: rfactor ! inflation factor, multiplies R for updating the ensemble anomalies.
    ! directory for some diagnostic output, char_missing -> no output made
    character(len=*), intent(in) :: diagn_out_dir 

    character(len=*), parameter :: sub_name = 'run_enkf'
    real, dimension(:), pointer :: obs_ens_mean, stdev_inv, innov_mean, tmp
    integer :: ens_size, obs_size, mdl_size
    integer :: stat, ind_obs, ind_ens, ind_mdl, pseudoinv_count
    real, dimension(:,:), allocatable :: obs_anom, &  ! S, (ind_ens, ind_obs) &
         & subS, &! local subset of S
         & X1, &  ! either SS^T + I or S^TS + I depending on which is smaller
         & Gmatr, Gtransp,&   ! gain matrix 
         & X5tmp, &     ! local subset of X5
         & obs_pert, &     ! observation perturbation in ENKF
         & subD            ! local subset of obs_pert
    real, dimension(:), allocatable :: subdy, & ! local subset of innov_mean
         & lfactors, & ! localization factors
         & obs_dist    ! distance to the local observations
    integer, dimension(:), allocatable :: obs_ind ! indices of locally selected observations
    real :: sqrtminv, lon_prev, lat_prev
    integer :: loc_obs_found, num_loc_obs, skip_count, ind_first, ind_last, ind_chunk, &
         & num_omp_chunks, chunk_target_size
    integer, dimension(:), pointer :: loc_obs_used
    real, dimension(size(ens,2)) :: grpt_ens
    integer, dimension(:), pointer :: chunk_starts, chunk_sizes
    
    !call set_error('Crash!', sub_name)
    !return

    call start_count('enkf_prepare')
    ens_size = size(ens, 2)
    mdl_size = size(ens, 1)
    obs_size = size(obs_data)
    if (fu_fails(ens_size > 1, 'ens_size < 2', sub_name)) return
    if (fu_fails(size(obs_var) == obs_size, 'obs_var and obs_data don''tt conform', sub_name)) return
    if (.not. fu_index(loc_type, (/loc_none, loc_step, loc_gaspari_cohn/)) > 0) then 
      call set_error('Bad loc_type', sub_name)
    end if
    if (.not. (enkf_flavor == enkf_enkf .or. enkf_flavor == enkf_denkf)) then
      call set_error('Bad enkf_flavor', sub_name)
    end if
    if (.not. loc_dist_m > 0.0 .and. loc_dist_m < 2e7) then
      call msg('loc_dist_m:', loc_dist_m)
      call set_error('Strange loc_dist_m', sub_name)
    end if
    if (fu_fails(rfactor > 0, 'Strange rfactor', sub_name)) return
    if (error) return

    call msg('ENKF: Dimension of observation vector:', obs_size)
    call msg('ENKF: Dimension of state vector:', mdl_size)

    if (fu_fails(obs_size > 0, 'No observations given', sub_name)) return
     
    allocate(obs_anom(obs_size, ens_size), loc_obs_used(mdl_size), stat=stat)
    if (fu_fails(stat == 0, 'Allocate failed', sub_name)) return
    loc_obs_used(1:mdl_size) = 0

    ! calculate S = R^(-1/2) . (HE - mean(HE, dim=ens)) / sqrt(m-1), m = ensemble size
    sqrtminv = 1/sqrt(real(ens_size) - 1)
    stdev_inv => fu_work_array(obs_size)
    stdev_inv(1:obs_size) = 1 / sqrt(obs_var(1:obs_size))
    obs_ens_mean => fu_work_array(obs_size)
    if (error) return
    do ind_obs = 1, obs_size
      obs_ens_mean(ind_obs) = sum(ens_obs(ind_obs, :)) / ens_size
      !print *, 'OEM', obs_ens_mean(ind_obs)
      do ind_ens = 1, ens_size
        obs_anom(ind_obs, ind_ens) &
             & = sqrtminv * stdev_inv(ind_obs) * (ens_obs(ind_obs, ind_ens) - obs_ens_mean(ind_obs)) 
      end do
    end do
    
    innov_mean => fu_work_array(obs_size)
    if (error) return
    ! calculate s = R^(-1/2) . (d - Hx) / sqrt(m-1), m = ensemble size, 
    ! d = data = obs, s is d or dy in TOPAZ code
    do ind_obs = 1, obs_size
      innov_mean(ind_obs) = sqrtminv * stdev_inv(ind_obs) * (obs_data(ind_obs) - obs_ens_mean(ind_obs))
    end do
    
    if (enkf_flavor == enkf_enkf) then
      allocate(obs_pert(obs_size,ens_size), stat=stat)
      if (fu_fails(stat == 0, 'Allocate failed', sub_name)) return
      do ind_obs = 1, obs_size
        call random_normal(obs_pert(ind_obs,:))
        ! topaz code has here rfactor without sqrt ??
        obs_pert(ind_obs,:) = obs_pert(ind_obs,:) * sqrtminv / sqrt(rfactor)
      end do
    end if
    
    ! preparation done.
    call stop_count('enkf_prepare')
    call start_count('enkf_update')

    pseudoinv_count = 0
    skip_count = 0
    
    ! The OpenMP parallelisation is along the mdl points. This is done by splitting the
    ! model points to a number of chunks and looping in OpenMP over those, with the
    ! following rationale:
    
    ! - It is is useful to keep the order the points are analysed: when the points have
    ! the same localisation, the ensemble transform doesn't need to be recomputed. 
    ! - On the other hand, some grid points might be skipped due to having no local
    ! observations, or due to not having ensemble spread, which results in load imbalance.
    ! 
    ! Thus, we make the number of chunks higher than number of threads, and try to deal
    ! with the load imbalance by using SCHEDULE(dynamic).
    
    chunk_target_size = 100
    chunk_starts => fu_work_int_array()
    chunk_sizes => fu_work_int_array()
    num_omp_chunks = mdl_size / chunk_target_size + 1
    call msg('ENKF: number of OpenMP chunks', num_omp_chunks)
    call get_chunk(mdl_size, num_omp_chunks, chunk_starts, chunk_sizes, if_zero_based=.false.)
    if (error) return
    

    !$OMP PARALLEL DEFAULT(NONE) SHARED(ens, &
    !$OMP & ens_obs, obs_data, obs_var, obs_loc, mdl_loc, loc_dist_m, loc_type, enkf_flavor, rfactor, &
    !$OMP & diagn_out_dir, obs_ens_mean, stdev_inv, innov_mean, ens_size, obs_size, mdl_size, obs_anom, &
    !$OMP & obs_pert, loc_obs_used, error, chunk_starts, chunk_sizes, num_omp_chunks) &
    !$OMP & PRIVATE(stat, ind_obs, ind_ens, ind_mdl, subS, X1, Gmatr, Gtransp, X5tmp, subD, &
    !$OMP & lfactors, obs_dist, obs_ind, lon_prev, lat_prev, loc_obs_found, subdy, num_loc_obs, &
    !$OMP & grpt_ens, ind_first, ind_last, ind_chunk) &
    !$OMP & REDUCTION(+:skip_count,pseudoinv_count)

    allocate(obs_ind(obs_size), obs_dist(obs_size), &
           & X5tmp(ens_size, ens_size), &
           & stat=stat)
    if (fu_fails(stat == 0, 'Allocate failed', sub_name)) continue

    ! pre-allocate local arrays for 100 obs (enlarge if needed)
    call allocate_loc_arr(2000, ens_size, subS, subdy, lfactors, Gmatr, Gtransp, subD, X1)

    !tmp => fu_work_array()

    ! go cell by cell, calculate local X5 and update
    loc_obs_found = 0
    lon_prev = real_missing
    lat_prev = real_missing
    
    ! ****
    ! Loop over gridpoints
    !$OMP DO SCHEDULE(dynamic)
    do ind_chunk = 1, num_omp_chunks
      ind_first = chunk_starts(ind_chunk)
      ind_last = ind_first + chunk_sizes(ind_chunk) - 1
      !print *, 'thread chunk', omp_get_thread_num(), ind_chunk, ind_first, ind_last
      do ind_mdl = ind_first, ind_last

        if (error) cycle

        !if (fu_str(int(ind_mdl/real(mdl_size)*10)) /= fu_str(int((ind_mdl-1)/real(mdl_size)*10))) then
        !  call msg('ENKF: gridpoints done, total', ind_mdl, mdl_size)
        !end if
        grpt_ens = ens(ind_mdl,:)

        ! If the spread is very small (possibly because all values are near-zero), skip the analysis:
        if (maxval(grpt_ens) - minval(grpt_ens) < max(1e-12, 1e-6*sum(grpt_ens)/ens_size)) then
          skip_count = skip_count + 1
          cycle
        end if

        ! If new horizontal gridpoint, get new local observations and ensemble transformation X5.
        !
        if (.not. ((mdl_loc(1,ind_mdl) .eps. lon_prev) .and. (mdl_loc(2,ind_mdl) .eps. lat_prev))) then
          lon_prev = mdl_loc(1,ind_mdl)
          lat_prev = mdl_loc(2,ind_mdl)
          call get_local_obs(loc_type, mdl_loc(1,ind_mdl), mdl_loc(2,ind_mdl), obs_loc, &
               & loc_dist_m, obs_ind, obs_dist, num_loc_obs)
          if (num_loc_obs == 0) then
            cycle
          end if

          loc_obs_used(ind_mdl) = num_loc_obs
          
          !use_obs_space = .false.
          !use_obs_space = .true.
          if (size(subdy) < num_loc_obs) then 
            !if (.true.) then 
            deallocate(subS, subdy, lfactors, Gmatr, Gtransp, X1, subD)
            loc_obs_found = num_loc_obs
          end if
          if (.not. allocated(subS)) then
            call allocate_loc_arr(num_loc_obs, ens_size, subS, subdy, lfactors, Gmatr, Gtransp, subD, X1)
          end if
          if (error) cycle
          ! Extract local sub-matrices of S, dy
          subS(1:num_loc_obs,:) = obs_anom(obs_ind(1:num_loc_obs),:)
          subdy(1:num_loc_obs) = innov_mean(obs_ind(1:num_loc_obs))
          if (enkf_flavor == enkf_enkf) subD(1:num_loc_obs,:) = obs_pert(obs_ind(1:num_loc_obs), :)

          call apply_locfun(loc_type, enkf_flavor, loc_dist_m, subS, subdy, subD, obs_dist, num_loc_obs, lfactors)
          if (error) cycle

          call get_x5_blas(subS, X1, Gtransp, Gmatr, subdy, X5tmp, subD, rfactor, ens_size, num_loc_obs, &
               & pseudoinv_count, enkf_flavor)

          if (error) cycle

        else if(num_loc_obs == 0) then
          cycle
        end if
        ! now have local X5 -> update the ensemble
        ! y = alpha*Ax + beta*y
        ! SUBROUTINE SGEMV ( TRANS, M, N, ALPHA, A, LDA, X, INCX, BETA, Y, INCY )
        call sgemv('t', ens_size, ens_size, 1.0, X5tmp, ens_size, grpt_ens, 1, 0.0, ens(ind_mdl,1), mdl_size)

      end do ! over state vector
    end do ! over OpenMP chunks
    !$OMP END DO

    deallocate(subS, subdy, lfactors, Gmatr, Gtransp, X1, subD)

    !$OMP END PARALLEL
    deallocate(loc_obs_used)

    call stop_count('enkf_update')
    call msg('ENKF: number of pseudoinverses:', pseudoinv_count)
    
    ! Note: if the filter is used with domain decomposition, most or all gridpoints might end up skipped!
    call msg('ENKF: number of local analyses skipped:', skip_count)
    
    if (diagn_out_dir /= char_missing) then
      ! Write the number of observations used for each analysis point.
      call dump_int_arr(loc_obs_used(1:mdl_size), diagn_out_dir + dir_slash + 'loc_obs.dat')
    end if

    !call free_work_array(loc_obs_used)
    call free_work_array(innov_mean)
    call free_work_array(stdev_inv)
    call free_work_array(obs_ens_mean)
    call free_work_array(chunk_starts)
    call free_work_array(chunk_sizes)

  end subroutine enkf_update_openmp
  
  !************************************************************************************

  subroutine get_x5_blas(subS, X1, Gtransp, Gmatr, subdy, X5tmp, subD, rfactor, ens_size, &
                       & num_loc_obs, pseudoinv_count, &
                       & enkf_flavor)
    !
    ! Internal routine for evaluating the local ensemble transformation X5.
    ! 
    implicit none
    real, dimension(:), intent(in) :: subdy
    real, dimension(:,:), intent(inout) :: X5tmp, subD, subS, X1, Gtransp, Gmatr
    real, intent(in) :: rfactor
    integer, intent(inout) :: pseudoinv_count
    integer, intent(in) :: ens_size, num_loc_obs, enkf_flavor

    integer :: iter, ind_ens, ind_obs, stat
    character(len=*), parameter :: sub_name = 'get_x5_blas'
    logical :: use_obs_space

    use_obs_space = num_loc_obs < ens_size

    ! calculate the ensemble transformation matrix (X5 in Sakov & Oke 2008 paper) using
    ! BLAS subroutines.
    do iter = 1, 2
      if (iter == 2 .and. .not. (rfactor .eps. 1.0)) then
        subS = subS / sqrt(rfactor)
      end if

      ! The analysis will use observation or ensemble space depending on which is
      ! smaller. Depending on the choice, it is convenient to calculate either G or G^T. 
      if (use_obs_space) then
        ! G = S^T (I + SS^T)^-1
        ! X1 = SS^T + I
        ! SSYRK - perform one of the symmetric rank k operations   C := alpha*A*A' + beta*C
        ! SUBROUTINE SSYRK ( UPLO, TRANS, N, K, ALPHA, A, LDA, BETA, C, LDC )
        call ssyrk('u', 'n', num_loc_obs, ens_size, 1.0, subS, size(subS,1), 0.0, X1, size(X1,1))
        ! note: only the upper triangular part of X1 is stored.
        do ind_obs = 1, num_loc_obs
          X1(ind_obs, ind_obs) = X1(ind_obs, ind_obs) + 1.0
        end do
        Gtransp(1:num_loc_obs, 1:ens_size) = subS(1:num_loc_obs, 1:ens_size)
        call multiply_inv(X1, Gtransp, num_loc_obs, ens_size, &
             & size(X1,1), size(Gtransp,1), 'cholesky', stat)
        if (stat /= 0) then
          ! X1 is close to singular. Try with pseudoinverse, but need to recompute X1 first.
          call ssyrk('u', 'n', num_loc_obs, ens_size, 1.0, subS, size(subS,1), 0.0, X1, size(X1,1))
          do ind_obs = 1, num_loc_obs
            X1(ind_obs, ind_obs) = X1(ind_obs, ind_obs) + 1.0
          end do
          call multiply_inv(X1, Gtransp, num_loc_obs, ens_size, &
               & size(X1,1), size(Gtransp,1), 'pseudoinv', stat)
          pseudoinv_count = pseudoinv_count + 1
        end if
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !if (error .or. fu_fails(stat == 0, 'ENKF: Problem with X1, stat = '//fu_str(stat), sub_name)) then
        if (error .or. fu_fails(stat == 0, 'ENKF: Problem with X1', sub_name)) then
          write (*,*) stat
          call msg('ENKF: dumping matrices, num_loc_obs = ', num_loc_obs)
          !call msg('X1:', reshape(X1(1:num_loc_obs,1:num_loc_obs),(/num_loc_obs*num_loc_obs/)))
          call matrixdump(subS(1:num_loc_obs, 1:ens_size), 'subS.dat')
          call matrixdump(X1(1:num_loc_obs, 1:num_loc_obs), 'X1.dat')
        end if

      else ! model space
        ! G = (I + S^T S)^-1 S^T,
        ! X1 = I + S^T S
        X1(1:ens_size, 1:ens_size) = matmul(transpose(subS(1:num_loc_obs,:)), subS(1:num_loc_obs,:))
        do ind_ens = 1, ens_size
          X1(ind_ens, ind_ens) = X1(ind_ens,ind_ens) + 1.0
        end do
        Gmatr(1:ens_size, 1:num_loc_obs) = transpose(subS(1:num_loc_obs, 1:ens_size))
        call multiply_inv(X1, Gmatr, ens_size, num_loc_obs, size(X1,1), size(Gmatr, 1), 'cholesky', stat)
        if (stat /= 0) then
          ! X1 is close to singular. Try with pseudoinverse, but need to recompute X1 first.
          X1(1:ens_size, 1:ens_size) = matmul(transpose(subS(1:num_loc_obs,:)), subS(1:num_loc_obs,:))
          do ind_ens = 1, ens_size
            X1(ind_ens, ind_ens) = X1(ind_ens,ind_ens) + 1.0
          end do
          call multiply_inv(X1, Gmatr, ens_size, num_loc_obs, size(X1,1), size(Gmatr, 1), &
               & 'pseudoinv', stat)
          pseudoinv_count = pseudoinv_count + 1
        end if
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !if (fu_fails(stat == 0, 'ENKF: Problem inverting X1, stat = '//fu_str(stat), sub_name)) return
        if (fu_fails(stat == 0, 'ENKF: Problem inverting X1', sub_name)) then
          write(*,*) stat
          return
        end if
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      end if ! model/obs space inversion

      if (iter == 1) then
        ! evaluate Gs1^T in equation below 17 in the Sakov 2008 paper
        ! we'll store it to X5tmp. Second round goes (possibly) with different r facrtor, 
        ! so Gmatr and X1 will change.
        if (use_obs_space) then
          do ind_ens = 1, ens_size
            X5tmp(ind_ens, :) = sum(Gtransp(1:num_loc_obs,ind_ens) * subdy(1:num_loc_obs))
          end do
        else
          do ind_ens = 1, ens_size
            X5tmp(ind_ens, :) = sum(Gmatr(ind_ens,1:num_loc_obs) * subdy(1:num_loc_obs))
          end do
        end if
      end if

      ! calculate DFS at iteration 1, SRF at iteration 2 ??

    end do ! iter

    if (use_obs_space) then
      ! note Gmatr vs Gtransp
      if (enkf_flavor == enkf_enkf) then
        ! classic enkf: T = G (D - S)
        !X5tmp = X5tmp + matmul(transpose(Gtransp(1:num_loc_obs,1:ens_size)), &
        !                     & subD(1:num_loc_obs,:) - subS(1:num_loc_obs,:))
        subD(1:num_loc_obs,:) = subD(1:num_loc_obs,:) - subS(1:num_loc_obs,:)
        call sgemm('t', 'n', ens_size, ens_size, num_loc_obs, 1.0, &
             & Gtransp, size(Gtransp, 1), subD, size(subD,1), 1.0, X5tmp, size(X5tmp, 1))
      else
        ! denkf scheme: T = -(1/2) * GS
        ! SGEMM - perform one of the matrix-matrix operations   C := alpha*op( A )*op( B ) + beta*C,
        ! SUBROUTINE SGEMM ( TRANSA, TRANSB, M, N, K, ALPHA, A, LDA, B, LDB, BETA, C, LDC )
        call sgemm('t', 'n', ens_size, ens_size, num_loc_obs, -0.5, &
             & Gtransp, size(Gtransp,1), subS, size(subS, 1), 1.0, X5tmp, size(X5tmp, 1))
      end if
    else
      if (enkf_flavor == enkf_enkf) then
        subD(1:num_loc_obs,:) = subD(1:num_loc_obs,:) - subS(1:num_loc_obs,:)
        call sgemm('n', 'n', ens_size, ens_size, num_loc_obs, 1.0, &
             & Gmatr, size(Gmatr, 1), subD, size(subD,1), 1.0, X5tmp, size(X5tmp, 1))
        !X5tmp = X5tmp + matmul(Gmatr(1:ens_size,1:num_loc_obs), &
        !                     & subD(1:num_loc_obs,:) - subS(1:num_loc_obs,:))
      else
        ! denkf scheme: add T = -(1/2) * GS
        call sgemm('n', 'n', ens_size, ens_size, num_loc_obs, -0.5, &
             & Gmatr, size(Gmatr,1), subS, size(subS, 1), 1.0, X5tmp, size(X5tmp, 1))
      end if

    end if

    do ind_ens = 1, ens_size
      X5tmp(ind_ens, ind_ens) = X5tmp(ind_ens,ind_ens) + 1.0
    end do

  end subroutine get_x5_blas

  subroutine multiply_inv(A, B, n, m, lda, ldb, method, stat)
    ! evaluate B -> A^-1 B using LAPACK subroutines, with A symmetric positively
    ! (semi-)definite.
    implicit none
    real, dimension(lda,*), intent(in) :: A ! n x n matrix
    real, dimension(ldb,*), intent(inout) :: B ! n x m matrix
    ! Need to give leading dimensions to comply with LAPACK interface and avoid
    ! temporary copies.
    integer, intent(in) :: n, m ! n = order of A, m = number of RHSs = number of columns in B
    integer, intent(in) :: lda, ldb
    character(len=*), intent(in) :: method
    integer, intent(out) :: stat

    real :: ssyev_work(3*n), eigval(n)
    real, parameter :: min_eigval_ratio = 1e-5
    real :: eigval_ratio
    integer :: npos, ind_ev, ii
    character(len=*), parameter :: sub_name = 'multiply_inv'

    select case(method)
    case ('cholesky')
      ! SUBROUTINE SPOTRF( UPLO, N, A, LDA, INFO )
      call spotrf('U', n, A, lda, stat)
      if (stat /= 0) return
      ! SUBROUTINE SPOTRS( UPLO, N, NRHS, A, LDA, B, LDB, INFO )
      call spotrs('U', n, m, A, lda, B, ldb, stat)
      if (fu_fails(stat == 0, 'ENKF: Error with spotrs', sub_name)) return
    case ('pseudoinv')
      ! SUBROUTINE SSYEV( JOBZ, UPLO, N, A, LDA, W, WORK, LWORK, INFO )
      call ssyev('v', 'u', n, A, lda, eigval, ssyev_work, size(ssyev_work), stat)
      if (fu_fails(stat == 0, 'ENKF: Error with ssyev', sub_name)) then
        call msg('ssyev stat', stat)
        return
      end if
      if (eigval(1) < eigval(n)*min_eigval_ratio) then
        do ind_ev = 1, n
          if (eigval(ind_ev) >= eigval(n)*min_eigval_ratio) exit
        end do
      else
        ind_ev = 1
      end if
      if (fu_fails(ind_ev < n+1, 'ENKF: Failed to truncate eigenvalues', sub_name)) return
      if (fu_fails(all(eigval(ind_ev:) >= 0.0), 'ENKF: Negative eigv', sub_name)) return
      npos = n - ind_ev + 1 ! number of positive eigenvalues
      !call msg('n, npos', n, npos)

      ! now A = USU^T, A^-1 = US^-1U^T
      ! set B -> U^T B, npos <= n
      B(1:npos, 1:m) = matmul(transpose(A(1:n, ind_ev:n)), B(1:n,1:m))
      ! B -> pinv(S)*B
      do ii = 1, npos
        B(ii,1:m) = B(ii,1:m) * (1.0/eigval(ind_ev+ii-1))
      end do
      ! B -> UB
      B(1:n,1:m) = matmul(A(1:n,ind_ev:n), B(1:npos,1:m))

    case default 
      call set_error('Unsupported inversion method', sub_name)
    end select

  end subroutine multiply_inv
  
  !************************************************************************************
  
  subroutine apply_locfun(loc_type, enkf_flavor, loc_dist, subS, subdy, subD, obs_dist, num_loc_obs, lfactors)
    implicit none
    integer, intent(in) :: loc_type
    real, intent(in) :: loc_dist
    integer, intent(in) :: enkf_flavor
    real, dimension(:,:), intent(inout) :: subS
    real, dimension(:), intent(inout) :: subdy
    real, dimension(:,:), intent(inout) :: subD
    real, dimension(:), intent(in) :: obs_dist
    integer, intent(in) :: num_loc_obs
    real, dimension(:), intent(out) :: lfactors

    integer :: ind_obs
    real :: dist_ratio, xx, xx2, xx3

    select case(loc_type)
    case (loc_none)
      lfactors(1:num_loc_obs) = 1.0
      return
    case (loc_step)
      where (obs_dist(1:num_loc_obs) < loc_dist) 
        lfactors(1:num_loc_obs) = 1.0
      elsewhere
        lfactors(1:num_loc_obs) = 0.0
      end where
    case (loc_gaspari_cohn)
      do ind_obs = 1, num_loc_obs
        dist_ratio = obs_dist(ind_obs) / loc_dist
        if (dist_ratio > 1) then
          lfactors(ind_obs) = 0.0
        else
          xx = dist_ratio * 2.0
          xx2 = xx**2
          xx3 = xx2 * xx
          if (xx < 1) then
            lfactors(ind_obs) = 1.0 + xx2 * (- xx3 / 4.0 + xx2 / 2.0)&
                 + xx3 * (5.0 / 8.) - xx2 * (5.0 / 3.0)
          else
            lfactors(ind_obs) = xx2 * (xx3 / 12.0 - xx2 / 2.0)&
                 + xx3 * (5.0 / 8.0) + xx2 * (5.0 / 3.0)&
                 - xx * 5.0 + 4.0 - (2.0 / 3.0) / xx
          end if
        end if ! dist_ratio > 1
      end do

    case default
      call set_error('loc_type no supported', 'apply_locfun')
      return

    end select

    do ind_obs = 1, num_loc_obs
      subS(ind_obs,:) = subS(ind_obs,:) * lfactors(ind_obs)
      subdy(ind_obs) = subdy(ind_obs) * lfactors(ind_obs)
    end do
    if (enkf_flavor == enkf_enkf) then
      do ind_obs = 1, num_loc_obs
        subD(ind_obs,:) = subD(ind_obs,:) * lfactors(ind_obs)
      end do
    end if

  end subroutine apply_locfun

  subroutine allocate_loc_arr(alloc_size, ens_size, subS, subdy, lfactors, Gmatr, Gtransp, subD, X1)
    implicit none
    real, dimension(:,:), allocatable :: subS, Gmatr, Gtransp, subD, X1
    real, dimension(:), allocatable :: lfactors, subdy
    integer, intent(in) :: alloc_size, ens_size

    integer :: size_x1_req, stat
    character(len=*), parameter :: sub_name = 'allocate_loc_arr'

    if (.not. allocated(subS)) then
      !call msg('Allocating: num_loc_obs, ens_size', (/alloc_size, ens_size, omp_get_thread_num()/))
      allocate(subS(alloc_size, ens_size), &
           & subdy(alloc_size), &
           & lfactors(alloc_size), &
           & Gmatr(ens_size, alloc_size), Gtransp(alloc_size, ens_size), &
           & subD(alloc_size, ens_size), stat=stat)
      if (fu_fails(stat == 0, 'Allocate subS etc failed', sub_name)) return

      size_x1_req = min(ens_size, alloc_size)
      if (allocated(X1)) then
        ! X1 needs to be either ensemble or observation sized. The smaller will be used.
        if (size(X1, 1) < size_x1_req) then
          deallocate(X1)
          allocate(X1(size_x1_req, size_x1_req), stat=stat)
        end if
      else
        allocate(X1(size_x1_req, size_x1_req), stat=stat)
      end if
      if (fu_fails(stat == 0, 'Allocate X1 failed', sub_name)) return
    end if

  end subroutine allocate_loc_arr

  subroutine get_local_obs(loc_type, lon, lat, obs_loc, &
                         & loc_dist_m, obs_ind, &
                         & obs_dist, num_loc_obs)
    ! Find the observations within the localisation radius and store their indices and
    ! distance from the analysis point.
    !
    implicit none
    real, intent(in) :: lon, lat
    real, dimension(:,:), intent(in) :: obs_loc
    real, intent(in) :: loc_dist_m
    integer, dimension(:), intent(out) :: obs_ind
    real, dimension(:), intent(out) :: obs_dist
    integer, intent(in) :: loc_type
    integer, intent(out) :: num_loc_obs

    integer :: ind_obs, num_obs
    real :: dxmin, loc_dist_km2, distance, dlat, dlon, obslat, obslon
    real, parameter :: dy = earth_radius*degrad*1e-3

    num_obs = size(obs_loc, 2)
    loc_dist_km2 = (1e-3*loc_dist_m)**2
    ! if loc_none is used, this sub will collect all observations!
    if (loc_type == loc_none) loc_dist_km2 = huge(kind(loc_dist_km2))

    num_loc_obs = 0
    do ind_obs = 1, num_obs
      obslon = obs_loc(1,ind_obs)
      obslat = obs_loc(2,ind_obs)
      dxmin = fu_dx_deg_to_m(1.0, max(abs(lat), abs(obslat)))*1e-3
      dlon = abs(lon-obslon)
      dlat = abs(lat-obslat)
      if ((dlon*dxmin)**2 + (dlat*dy)**2 > loc_dist_km2) cycle
      distance = fu_gc_distance(lon, obslon, lat, obslat)
      if (distance > loc_dist_m .and. loc_type /= loc_none) cycle
      num_loc_obs = num_loc_obs + 1
      obs_ind(num_loc_obs) = ind_obs
      obs_dist(num_loc_obs) = distance
    end do

  end subroutine get_local_obs

  !************************************************************************************

  subroutine dump_int_arr(arr, filename)
    implicit none
    integer, dimension(:), intent(in) :: arr
    character(len=*), intent(in) :: filename
    
    integer :: file_unit, iostat

    file_unit = fu_next_free_unit()
    open(file_unit, file=filename, access='stream', form='unformatted', status='replace', iostat=iostat)
    if (fu_fails(iostat == 0, 'Failed to open: ' // trim(filename), 'dump_int_arr')) return
    
    write(file_unit) arr
    close(file_unit)

  end subroutine dump_int_arr

  subroutine get_row(ens, ind_row, row, stdev, mean, corr_or_cov)
    ! Computes a row from the ensemble covariance matrix, along with the mean and stdev.
    ! Either covariance or correlation can be requested.
    implicit none
    real, dimension(:,:), intent(in) :: ens
    integer, intent(in) :: ind_row
    real, dimension(:), intent(out) :: row, stdev, mean
    character(len=*), intent(in) :: corr_or_cov

    integer :: ens_size, mdl_size, ind_ens, ind_mdl
    character(len=*), parameter :: sub_name = 'get_row'
    
    ens_size = size(ens,2)
    mdl_size = size(ens,1)

    if (fu_fails(size(row) >= mdl_size, 'row too small', sub_name)) return
    if (fu_fails(size(stdev) == mdl_size, 'stdev size doesn''t match', sub_name)) return
    if (fu_fails(size(stdev) == size(mean), 'stdev and mean size don''t match', sub_name)) return
    if (fu_fails(ind_row <= mdl_size .and. ind_row > 0, 'bad ind_mdl', sub_name)) return

    mean = sum(ens, dim=2) / ens_size
    stdev = 0.0
    do ind_ens = 1, ens_size
      stdev = stdev + (ens(:,ind_ens) - mean)**2
    end do
    stdev = sqrt(stdev / ens_size)

    select case(corr_or_cov)
    case('cov')
      row = 0.0
      do ind_ens = 1, ens_size
        do ind_mdl = 1, mdl_size
          row(ind_mdl) &
               & = sum( (ens(ind_mdl,:)-mean(ind_mdl)) * (ens(ind_row,:)-mean(ind_row)) ) / ens_size
        end do
      end do

    case ('corr')
      row = 0.0
      do ind_ens = 1, ens_size
        do ind_mdl = 1, mdl_size
          row(ind_mdl) &
               & = sum( (ens(ind_mdl,:)-mean(ind_mdl)) * (ens(ind_row,:)-mean(ind_row)) ) &
               & / (ens_size*stdev(ind_mdl)*stdev(ind_row))
        end do
      end do
      
    case default
      call set_error('Bad corr_or_cov', sub_name)
      return
    end select

  end subroutine get_row
  
  subroutine matrixdump(matrix, filename)
    implicit none
    real, dimension(:,:), intent(in) :: matrix
    character(len=*), intent(in) :: filename

    integer :: file_unit, iostat

    file_unit = fu_next_free_unit()
    open(file_unit, file=filename, access='stream', form='unformatted', iostat=iostat, status='replace')
    if (fu_fails(iostat == 0, 'Failed to open: '//trim(filename), 'matrixdump')) return
    write(file_unit) matrix
    close(file_unit)

  end subroutine matrixdump

end module enkf
