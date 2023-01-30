MODULE silam_partitioning
  !
  ! This module contains the higher level domain decomposition and MPI routines
  !
  USE grids_geo
  USE silam_times
  !$ use OMP_LIB

#ifdef PUHTI_BUG
!!!#pragma message("PUHTI_BUG enabled")
    use,intrinsic :: ieee_arithmetic
#endif
  IMPLICIT NONE

  PRIVATE

  ! This type has the decomposition info for parallel Grads output
  TYPE decomposition
     INTEGER :: size_x = int_missing, size_y = int_missing     ! Size of the decomposed domain
     INTEGER :: offset_x = int_missing, offset_y = int_missing ! Offsets related to the global domain
     INTEGER :: globalNx = int_missing, globalNy = int_missing ! Size of the global domain in cells
     INTEGER :: maxsize_x = int_missing, maxsize_y = int_missing ! Max disp. grid size among domains
     INTEGER :: subarray_mpi_type = int_missing                ! MPI datatype for file IO
  END TYPE decomposition

  TYPE(decomposition), PRIVATE :: output_decomposition
  LOGICAL :: decomposed_output_grid

  ! A few things are public
  INTEGER, PUBLIC :: adv_mpi_neighbours(4), &  ! Ranks of neighbors (N S E W)
       & wing_depth_e, wing_depth_w, &
       & wing_depth_n, wing_depth_s, &
       & wing_depth ! depth of the buffer zone, this is global
                    ! max of all wing depths
       
  INTEGER, PRIVATE :: buffer_size_x, buffer_size_y

  LOGICAL :: is_io_task
  INTEGER, PUBLIC :: smpi_ens_size, smpi_ens_ind
  INTEGER :: nof_iotasks, x_divisions, y_divisions
  REAL, DIMENSION (max_divisions) :: division_hint = real_missing
  INTEGER, PUBLIC, dimension(max_divisions):: smpi_x_sizes, smpi_y_sizes, smpi_x_offsets, smpi_y_offsets !Only (x|y)-master knows
  INTEGER :: my_x_coord, my_y_coord !Domain number in X or Y dimension 0:n_divisions-1 

  LOGICAL, PUBLIC, SAVE :: smpi_use_mpiio_grads
  LOGICAL, PUBLIC, SAVE :: smpi_use_mpiio_netcdf

  ! Functions here and in silam_mpi:
  !
  PUBLIC smpi_abort
  PUBLIC smpi_advection_barrier
  PUBLIC smpi_finalize
  PUBLIC smpi_init
  PUBLIC smpi_is_mpi_version
  PUBLIC smpi_setup_parameters
  PUBLIC smpi_get_process_xy_topology
  PUBLIC smpi_reset_periodic_topology
  PUBLIC fu_decompose_grid_mpi
  PUBLIC smpi_write_grads_field_parallel
  PUBLIC smpi_write_grads_fieldset_parallel
  PUBLIC smpi_read_grads_fieldset_parallel
  PUBLIC smpi_set_grids
  PUBLIC smpi_get_decomposition
  PUBLIC smpi_get_maxsubdomain_size
  public smpi_gather_field
  public smpi_bcast_aray
  public smpi_bcast_char
  public smpi_bcast_int_aray
  public smpi_bcast_int8_aray
  public smpi_allgather_char
  public smpi_bcast_string
  public smpi_send
  public smpi_recv
  public smpi_gatherv_real
  public smpi_scatterv_real
  PUBLIC smpi_close_gradsfile_mpiio
  PUBLIC smpi_exchange_boundaries
  PUBLIC smpi_exchange_wings
  PUBLIC smpi_open_gradsfile_mpiio_r
  PUBLIC smpi_open_gradsfile_mpiio_w
  PUBLIC smpi_global_barrier
  PUBLIC smpi_reduce_add
  PUBLIC smpi_reduce_max
  PUBLIC smpi_allreduce_add
  PUBLIC smpi_allreduce_max_int


CONTAINS

  !************************************************************************************

  SUBROUTINE smpi_get_maxsubdomain_size(nx, ny)
    !
    ! Return info needed for output chunking
    !
    IMPLICIT NONE
    INTEGER, INTENT(out) :: nx, ny
    nx   = output_decomposition%maxsize_x
    ny   = output_decomposition%maxsize_y

  END SUBROUTINE smpi_get_maxsubdomain_size

  !************************************************************************************

  SUBROUTINE smpi_get_decomposition(nx, ny, offx, offy, gnx, gny)
    !
    ! Return the information about the dispersion grid decomposition
    !
    IMPLICIT NONE
    INTEGER, INTENT(out) :: nx, ny, offx, offy, gnx, gny

    nx   = output_decomposition%size_x
    ny   = output_decomposition%size_y
    offx = output_decomposition%offset_x
    offy = output_decomposition%offset_y
    gnx  = output_decomposition%globalNx
    gny  = output_decomposition%globalNy

  END SUBROUTINE smpi_get_decomposition

  !************************************************************************************

  SUBROUTINE smpi_get_process_xy_topology(my_xc, my_yc, size_x, size_y)
    !
    ! Return the information about the domain decomposition process layout
    !
    IMPLICIT NONE
    INTEGER, INTENT(out) :: my_xc, my_yc, size_x, size_y

#ifdef SILAM_MPI    
    IF(is_io_task)THEN
       CALL set_error('IO tasks do not have y-coordinate in advection communicator!', &
            & 'smpi_get_process_xy_topology')
       RETURN
    END IF

    my_xc = my_x_coord
    my_yc = my_y_coord
    size_x = x_divisions
    size_y = y_divisions
#else    
    my_xc = 0
    my_yc = 0
    size_x = 1
    size_y = 1
#endif
  END SUBROUTINE smpi_get_process_xy_topology

  !************************************************************************************

  SUBROUTINE smpi_set_grids(wholeMPI_grid_disp, grid_disp, grid_out, timestep)
    !
    ! With the whole-domain dispersion, output and meteo grids defined, re-define them
    ! according to the domain decomposition.
    !
    ! - if output grid == dispersion grid initially, both are decomposed. If they differ,
    !   the subdomains will keep the full-domain output grid!
    !
    IMPLICIT NONE
    TYPE(silja_grid), INTENT(in) :: wholeMPI_grid_disp
    TYPE(silja_grid), INTENT(inout) :: grid_disp, grid_out
    TYPE(silja_interval), INTENT(in) :: timestep

    INTEGER ::  nx, ny, ix, iy, stat
    integer :: iMPI, jMPI,nxMPI,nyMPI
    real, dimension(2*max_divisions) :: sendbuf, recvbuf
    logical :: ifOk

#ifdef SILAM_MPI

    CALL msg('Setting MPI grids')


    IF (wholeMPI_grid_disp == grid_out) THEN
      ! First save the global dimensions to output_decomposition, then decompose
      ! the grid
      call grid_dimensions(grid_out, output_decomposition%globalNx, output_decomposition%globalNy)

      decomposed_output_grid = .TRUE.
      grid_out = fu_decompose_grid_mpi(grid_out, timestep)
    ELSE
      CALL msg_warning('Output grid not decomposed', 'smpi_set_grids')
    END IF

    IF (fu_gridtype(wholeMPI_grid_disp) /= lonlat) THEN
       ! anygrids can be decomposed without problem, but the corresponding fields would
       ! need to be handled as well.
       call msg("disp. grid:")
       call report(wholeMPI_grid_disp)
       CALL set_error('Need lonlat grid for MPI', 'smpi_set_grids')
       RETURN
    END IF
#endif

    ! Save the original dispersion grid for later use and then decompose
    grid_disp = fu_decompose_grid_mpi(wholeMPI_grid_disp, timestep)


    CALL grid_dimensions(grid_disp, nx, ny)
    ! Get the cell sizes of dispersion grid
    ALLOCATE(disp_grid_size_x(nx, ny), &
         & disp_grid_size_y(nx, ny), stat=stat)
    IF (fu_fails(stat == 0, 'Allocate failed', 'smpi_set_grid')) RETURN
    DO iy = 1, ny
       DO ix = 1, nx
          disp_grid_size_x(ix,iy) = fu_dx_cell_m(grid_disp, ix, iy)
          disp_grid_size_y(ix,iy) = fu_dy_cell_m(grid_disp, ix, iy)
       END DO
    END DO
#ifdef SILAM_MPI
    !Max mpi subdomain sizes needed for collective exchange
    CALL smpi_allreduce_max_int(nx, output_decomposition%maxsize_x, smpi_adv_comm)
    CALL smpi_allreduce_max_int(ny, output_decomposition%maxsize_y, smpi_adv_comm)

    !!!! Get whole MPI grid structure to master
     sendbuf(1:2*max_divisions) = 0.
     recvbuf(1:2*max_divisions) = float(int_missing)
     sendbuf(1 + my_x_coord) = nx
     sendbuf(1 + max_divisions + my_x_coord) = output_decomposition%offset_x
     call smpi_reduce_max(sendbuf, recvbuf, 2*max_divisions, 0, smpi_adv_x_comm, ifOk)
     smpi_x_sizes(1:max_divisions)=nint(recvbuf(1:max_divisions))
     smpi_x_offsets(1:max_divisions)=nint(recvbuf(1+max_divisions:2*max_divisions))

     sendbuf(1:2*max_divisions) = 0.
     recvbuf(1:2*max_divisions) = float(int_missing)
     sendbuf(1 + my_y_coord) = ny
     sendbuf(1 + max_divisions + my_y_coord) = output_decomposition%offset_y
     call smpi_reduce_max(sendbuf, recvbuf, 2*max_divisions, 0, smpi_adv_y_comm, ifOk)
     smpi_y_sizes(1:max_divisions)=nint(recvbuf(1:max_divisions))
     smpi_y_offsets(1:max_divisions)=nint(recvbuf(1+max_divisions:2*max_divisions))

     CALL smpi_get_process_xy_topology(iMPI, jMPI, nxMPI, nyMPI)
     call msg("My rank X,Y", my_x_coord, my_y_coord)
     if (my_x_coord == 0) call msg("I'm X-master. Xsizes  :", smpi_x_sizes(1:nxMPI))
     if (my_x_coord == 0) call msg("I'm X-master. Xoffsets:", smpi_x_offsets(1:nxMPI))
     if (my_y_coord == 0) call msg("I'm Y-master. Ysizes  :", smpi_y_sizes(1:nyMPI))
     if (my_y_coord == 0) call msg("I'm Y-master. Yoffsets:", smpi_y_offsets(1:nyMPI))


    CALL msg('')

    CALL create_subarray_datatype()
#else
     smpi_x_sizes(1)=nx
     smpi_x_offsets(1)=0
     smpi_y_sizes(1)=ny
     smpi_y_offsets(1)=0
     output_decomposition%globalNx = nx
     output_decomposition%globalNy = ny

     output_decomposition%offset_x = 0
     output_decomposition%offset_y = 0
     output_decomposition%size_x = nx
     output_decomposition%size_y = ny
     output_decomposition%maxsize_x = nx
     output_decomposition%maxsize_y = ny


#endif     


  END SUBROUTINE smpi_set_grids

  ! ************************************************************************************

  FUNCTION fu_decompose_grid_mpi(whole_dispersion_grid, timestep)RESULT(grid)
    !
    ! Determines the MPI topology and calls a generic grid decomposition.
    !
    IMPLICIT NONE

    ! The return value of this function:
    TYPE(silja_grid) :: grid

    ! Imported parameter with intent IN
    TYPE(silja_grid), INTENT(in) :: whole_dispersion_grid
    TYPE(silja_interval), INTENT(in) :: timestep

#ifdef SILAM_MPI    
    ! Local parameters
    INTEGER :: iMPI, jMPI, nxMPI, nyMPI, nxMY, nyMY, ierr
    INTEGER :: offset_x, offset_y, globalNx, globalNy, local_max
    real, dimension (2*max_divisions) :: sendbuf, recvbuf !exchange geometry
    logical :: ifOk

    ! Then get the decomposition rules
    CALL smpi_get_process_xy_topology(iMPI, jMPI, nxMPI, nyMPI)
    IF(error)RETURN

    ! Pick the small grid out of the big one
    grid = fu_pick_subgrid(whole_dispersion_grid, &
       & iMPI, jMPI, nxMPI, nyMPI, division_hint, offset_x, offset_y)

    ! Save the offset for output data type definitions, if they are used
    output_decomposition%offset_x = offset_x
    output_decomposition%offset_y = offset_y

    ! We need to have the buffer zone for advection. Assume max wind speed and calculate the size
    ! then store it to global variables.
    CALL grid_dimensions(grid, nxMY, nyMY)
    output_decomposition%size_x = nxMY
    output_decomposition%size_y = nyMY
    if (nxMPI > 1) then
     buffer_size_x = CEILING(max_wind_speed * abs(fu_sec(timestep)) / &
         & MIN(fu_dx_cell_m(grid, 1, 1), fu_dx_cell_m(grid, nxMY, nyMY)))
    else 
      buffer_size_x = 0
    endif

    if (nyMPI > 1) then
       buffer_size_y = CEILING(max_wind_speed * abs(fu_sec(timestep)) / &
           & MIN(fu_dy_cell_m(grid, 1, 1), fu_dy_cell_m(grid, nxMY, nyMY)))
    else
        buffer_size_y = 0
    endif

    local_max = max(buffer_size_x, buffer_size_y)
    
    ! Determine the maximum buffer width of all ranks
    CALL smpi_allreduce_max_int(local_max, wing_depth, smpi_adv_comm)
#else
    ! Simply copy the grid...
    !
    grid = whole_dispersion_grid
    wing_depth = 0
#endif    
    
    !
    ! Where to add this buffer?
    !
    wing_depth_e = 0
    wing_depth_w = 0
    wing_depth_n = 0
    wing_depth_s = 0

#ifdef SILAM_MPI
    IF(adv_mpi_neighbours(western_boundary) /= int_missing)THEN
       wing_depth_w = wing_depth
    ENDIF
    IF(adv_mpi_neighbours(eastern_boundary) /= int_missing) THEN
       wing_depth_e = wing_depth
    END IF

    IF(adv_mpi_neighbours(southern_boundary) /= int_missing)THEN
       wing_depth_s = wing_depth
    ENDIF
    IF(adv_mpi_neighbours(northern_boundary) /= int_missing) THEN
       wing_depth_n = wing_depth
    END IF

    !Single-task MPI 
    if (all( (/wing_depth_w, wing_depth_e, wing_depth_s, wing_depth_n/) ==0)) wing_depth = 0

    call msg("SMPI wing depth, MyReqX, MyReqY", (/wing_depth, buffer_size_x, buffer_size_y/))
    if (wing_depth > nxMY .or. wing_depth > nyMY) then
         call msg("Depth, nx, ny", (/wing_depth, nxMY, nyMY/))
         call set_error("Wing depth exceeds subdomain dimensions.","fu_decompose_grid_mpi")
         return
    endif
    !PRINT *, 'DEBUG_SI:', smpi_global_rank, output_decomposition%globalNx, output_decomposition%globalNy, &
    !     &                output_decomposition%offset_x, output_decomposition%offset_y, &
    !     &                output_decomposition%size_x,   output_decomposition%size_y

    ! Now we should have everything that is needed for the datatype
#endif    

  END FUNCTION fu_decompose_grid_mpi

  ! ************************************************************************************

  SUBROUTINE smpi_setup_parameters(nlPtr)
    !
    ! Set up the basic parameters for MPI runs from namelist
    !
    IMPLICIT NONE

    TYPE(Tsilam_namelist_group), POINTER, INTENT(in) :: nlPtr

    INTEGER :: ierr, iTmp
    CHARACTER(len=fnlen) :: nl_content

    CHARACTER(len=*), PARAMETER :: sub_name='smpi_setup_parameters'

    !! Both MPI and non-MPI
    itmp = fu_content_int(nlPtr, 'omp_num_threads')
    !$ if (itmp > 0) call omp_set_num_threads(iTmp)   


#ifdef SILAM_MPI    


    nof_iotasks = fu_content_int(nlPtr, 'io_tasks')
    IF (nof_iotasks == int_missing) nof_iotasks = 0
    IF (fu_fails(nof_iotasks == 0, 'IO tasks not implemented', sub_name)) RETURN

    if (fu_str_l_case(fu_content(nlptr, 'ensemble_tasks')) == 'automatic') then
      smpi_ens_size = smpi_global_tasks
    else
      smpi_ens_size = fu_content_int(nlptr, 'ensemble_size')
      IF (smpi_ens_size == int_missing) smpi_ens_size = 1
      IF (fu_fails(smpi_ens_size > 0, 'Bad ensemble_size', sub_name)) RETURN
    end if

    nl_content = fu_str_l_case(fu_content(nlptr, 'decomposition_method'))
    SELECT CASE (nl_content)
    CASE ('automatic')
       IF (smpi_ens_size > 1) THEN
          CALL set_error('ensemble_tasks > 1 requires manual decomposition', sub_name)
          RETURN
       END IF
       y_divisions = smpi_global_tasks - nof_iotasks
       x_divisions = 1
    CASE ('manual', '')
       x_divisions = fu_content_int(nlPtr, 'x_divisions')
       y_divisions = fu_content_int(nlPtr, 'y_divisions')
       nl_content = fu_content(nlPtr, 'division_hint')  ! Stripe width proportional to 
       if (len(trim(nl_content)) > 0) then
         call split_string(nl_content, ' ', division_hint, iTmp)
         if ( x_divisions /= 1 .or. iTmp /= y_divisions )  then
           call msg("division_hint = ",division_hint(1:iTmp))
           call msg("length of division_hint must be equal to y_divisions, but equals", iTmp)
           call msg("x_divisions,y_divisions", x_divisions,y_divisions)
           call set_error("Can't apply division hint to the split", sub_name)
           return
         endif
       endif

       IF (fu_fails(y_divisions /= int_missing .AND. y_divisions > 0, 'invalid y_divisions', sub_name)) RETURN
    CASE default
       CALL set_error('Strange decomposition_method', sub_name)
       RETURN
    END SELECT

    ! Currently the ensemble and domain decomposition can not used together
    IF (smpi_ens_size > 1 .AND. (x_divisions > 1 .OR. y_divisions > 1)) THEN
       CALL set_error('Using domain decomposition and ensemble together is not yet supported', &
            & 'smpi_setup_parameters')
       RETURN
    END IF

    IF (smpi_ens_size*(x_divisions*y_divisions + nof_iotasks) /= smpi_global_tasks) THEN
       CALL msg('Processes required for decomposition:', smpi_ens_size*(x_divisions*y_divisions + nof_iotasks))
       CALL set_error('Number of processes does not match request in control', 'smpi_setup_parameters')
       RETURN
    END IF

    IF (fu_content(nlPtr, 'use_mpiio') == 'YES') THEN
       IF (nof_iotasks > 0) THEN
          CALL set_error('gradsMPIIO can not be used with iotasks', 'smpi_setup_parameters')
          RETURN
       ELSE
          smpi_use_mpiio_grads = .TRUE.
       END IF
    ELSE
       smpi_use_mpiio_grads = .FALSE.
    END IF

    IF (fu_content(nlPtr, 'use_mpiio_netcdf') == 'YES') THEN
       IF (nof_iotasks > 0) THEN
          CALL set_error('ncMPIIO can not be used with iotasks', 'smpi_setup_parameters')
          RETURN
       ELSE
          smpi_use_mpiio_netcdf = .TRUE.
       END IF
    ELSE
       smpi_use_mpiio_netcdf = .FALSE.
    END IF

#else
    CALL msg('No-MPI run')
#endif 

    CALL init_communicators()

  END SUBROUTINE smpi_setup_parameters


  ! ************************************************************************************

  SUBROUTINE init_communicators()
    !
    ! Internal procedure to initialize the communicators. The grand idea is as follows:
    ! - MPI_COMM_WORLD is split to smpi_ensmember_comm such that the processes belonging to the
    ! same ensemble member are share the smpi_ensmember_comm.
    ! - smpi_ensmember_comm is split the smpi_adv_comm and smpi_io_comm if IO tasks are defined; 
    ! otherwise smpi_ensmember_comm == smpi_adv_comm
    ! - the smpi_enkf_comm is created for communications across the ensemble members as needed 
    ! for the EnKF analysis.
    ! 
    ! However, the reality is simpler:
    ! - IO tasks are not used
    ! - parallel advection is not allowed with EnKF
    ! so in practical runs, the only non-trivial commnunicators will be smpi_enkf_comm and
    ! smpi_adv_comm, which will be identical to MPI_COMM_WORLD. Although the idea was to
    ! eventually allow using MPI advection with EnKF, the communicators may need to be
    ! rethought, since the EnKF is now parallel in itself (by decomposing the state
    ! vector).
    
    IMPLICIT NONE

#ifdef SILAM_MPI    
    INTEGER :: io_color, adv_color, key, ierr
    INTEGER, DIMENSION(2) :: coords
    CHARACTER(len=*), PARAMETER :: sub_name = 'init_communicators'
    INTEGER :: smpi_adv_x_rank, smpi_adv_y_rank 

    IF (nof_iotasks < 0 .OR. x_divisions < 0 .OR. y_divisions < 0) THEN
       CALL set_error('init_communicators called without setting correct parameters', &
            & sub_name)
       RETURN
    END IF

    ! Check that the number of tasks and dimensions match
    IF (smpi_ens_size * (x_divisions * y_divisions + nof_iotasks) /= smpi_global_tasks) THEN
       CALL set_error('Number of required tasks job launch parameters do not match', &
            & sub_name)
       RETURN
    END IF

    ! First create the domain communicator
    !!!!CALL create_ens_comm()

    ! MPI errors should not kill Silam instantly. A chance for last breath needed.
    call MPI_Errhandler_set(MPI_COMM_WORLD, MPI_ERRORS_RETURN, ierr)
    if (smpi_error(ierr, "after MPI_Errhandler_set", sub_name)) return


    smpi_ens_ind = smpi_global_rank / (x_divisions * y_divisions + nof_iotasks)  + 1! index of the ensemble member
                                                            ! same for all procs within a memeber
    CALL mpi_comm_split(MPI_COMM_WORLD, smpi_ens_ind, smpi_global_rank, smpi_ensmember_comm, ierr)
    if (smpi_error(ierr, "after world to ensmember split", sub_name)) return
    call mpi_comm_rank(smpi_ensmember_comm, smpi_ens_rank, ierr)
    if (smpi_error(ierr, "after ensmember rank", sub_name)) return
    !MPI_COMM_WORLD 
    !      smpi_ensmember_comm (communication within ensemble member), smpi_ens_size, smpi_ens_rank 
    !     
    
    
    if (smpi_ens_size > 1) then
      call mpi_comm_split(MPI_COMM_WORLD, smpi_ens_rank, smpi_ens_ind, smpi_enkf_comm, ierr)
      if (smpi_error(ierr, "after  ens to enkf split", sub_name)) return
      call mpi_comm_rank(smpi_enkf_comm, smpi_ens_rank, ierr)
    if (smpi_error(ierr, "after smpi_enkf_comm rank", sub_name)) return
    else
      smpi_ens_rank = 0
      smpi_enkf_comm = MPI_UNDEFINED
    end if
    !!!!end create_ens_comm
    IF (error) RETURN

    CALL mpi_comm_rank(smpi_ensmember_comm, key, ierr) !Use same order for ranking
    if (smpi_error(ierr, "after smpi_ensmember_comm rank", sub_name)) return
    ! Generate the communicators using comm_split
    IF (key < x_divisions * y_divisions - nof_iotasks) THEN
       is_io_task = .FALSE.
       io_color   = MPI_UNDEFINED
       adv_color  = 1
    ELSE
       is_io_task = .TRUE.
       io_color   = 1
       adv_color  = MPI_UNDEFINED
    END IF

    CALL mpi_comm_split(smpi_ensmember_comm, adv_color, key, smpi_adv_comm, ierr)
    if (smpi_error(ierr, "after smpi_ensmember_comm to smpi_adv_comm split", sub_name)) return
    CALL mpi_comm_size(smpi_adv_comm, smpi_adv_tasks, ierr)
    if (smpi_error(ierr, "after smpi_adv_comm size", sub_name)) return
    CALL mpi_comm_split(smpi_ensmember_comm, io_color, key, smpi_io_comm, ierr)
    if (smpi_error(ierr, "after smpi_ensmember_comm to smpi_io_comm split", sub_name)) return

    ! Basic cartesian topology
    IF (is_io_task) THEN
       CALL mpi_comm_rank(smpi_io_comm, smpi_io_rank, ierr)
       if (smpi_error(ierr, "after  smpi_io_comm rank", sub_name)) return

    ELSE
       CALL mpi_comm_rank(smpi_adv_comm, smpi_adv_rank, ierr)
       if (smpi_error(ierr, "after  smpi_adv_comm rank", sub_name)) return

       ! Periodicity only in x-direction, poles in y-direction
       CALL mpi_cart_create(smpi_adv_comm, 2, (/x_divisions, y_divisions/), &
         & (/x_divisions /= 1,  y_divisions/=1 /), .TRUE., smpi_adv_cart_comm, ierr) 
            ! Set all BC periodic for a while, except if one task per direction
       ! They will be reset after grids settled
       if (smpi_error(ierr, "after mpi_cart_create", sub_name)) return
       CALL mpi_cart_shift(smpi_adv_cart_comm, 1, 1, &
            & adv_mpi_neighbours(southern_boundary), &
            & adv_mpi_neighbours(northern_boundary), ierr)
       if (smpi_error(ierr, "after mpi_cart_shift 1 1", sub_name)) return
       CALL mpi_cart_shift(smpi_adv_cart_comm, 0, 1, &
            & adv_mpi_neighbours(western_boundary), &
            & adv_mpi_neighbours(eastern_boundary), ierr)
       if (smpi_error(ierr, "after mpi_cart_shift 0 1", sub_name)) return
       ! Determine own position in the grid
       CALL mpi_comm_rank(smpi_adv_cart_comm, smpi_adv_cart_rank, ierr)
       if (smpi_error(ierr, "after smpi_adv_cart_comm rank", sub_name)) return
       CALL mpi_cart_coords(smpi_adv_cart_comm, smpi_adv_cart_rank, 2, coords, ierr)
       if (smpi_error(ierr, "after mpi_cart_coords", sub_name)) return

       my_x_coord = coords(1)
       my_y_coord = coords(2)
       WHERE (adv_mpi_neighbours == smpi_proc_null) adv_mpi_neighbours = int_missing
       CALL msg('all-periodic MPI topology created (will be reset later)')
       CALL report_neighbours()

       
       CALL mpi_comm_split(smpi_adv_comm, my_y_coord, key, smpi_adv_x_comm, ierr) !rows have same y
       if (smpi_error(ierr, "after smpi_adv_comm, to x split", sub_name)) return
       CALL mpi_comm_rank(smpi_adv_x_comm, smpi_adv_x_rank, ierr)
       if (smpi_error(ierr, "after smpi_adv_x rank", sub_name)) return
       if (fu_fails(my_x_coord == smpi_adv_x_rank, "my_x_coord /= smpi_adv_x_rank",sub_name)) return

       CALL mpi_comm_split(smpi_adv_comm, my_x_coord, key, smpi_adv_y_comm, ierr) !cols have same x
       if (smpi_error(ierr, "after smpi_adv_comm, to y split", sub_name)) return
       CALL mpi_comm_rank(smpi_adv_y_comm, smpi_adv_y_rank, ierr)
       if (smpi_error(ierr, "after smpi_adv_y rank", sub_name)) return
       if (fu_fails(my_y_coord == smpi_adv_y_rank, "my_y_coord /= smpi_adv_y_rank",sub_name)) return


    END IF

    call msg('Communicators initialized')
    call msg('Global rank: ', smpi_global_rank)
    call msg('Rank in ensemble_member: ', smpi_ens_rank)
    call msg('Advect rank: ', smpi_adv_cart_rank)
#else
  adv_mpi_neighbours = int_missing

#endif    
    
  END SUBROUTINE init_communicators

  ! ************************************************************************************

  SUBROUTINE smpi_reset_periodic_topology(ifXperiodic, ifYperiodic)
    
    logical, intent(in) :: ifXperiodic, ifYperiodic
    
    if (x_divisions == 1) then
      adv_mpi_neighbours(eastern_boundary) = int_missing
      adv_mpi_neighbours(western_boundary) = int_missing
    elseif (.NOT. ifXperiodic) then
      if (my_x_coord == 0) adv_mpi_neighbours(western_boundary) = int_missing
      if (my_x_coord == x_divisions - 1 ) adv_mpi_neighbours(eastern_boundary) = int_missing
    endif

    if (y_divisions == 1) then
      adv_mpi_neighbours(southern_boundary) = int_missing
      adv_mpi_neighbours(northern_boundary) = int_missing
    elseif (.NOT. ifYperiodic) then
      if (my_y_coord == 0) adv_mpi_neighbours(southern_boundary) = int_missing
      if (my_y_coord == y_divisions - 1) adv_mpi_neighbours(northern_boundary) = int_missing
    endif

    !reset wing depths (Could have been done more elegantly)

    IF(adv_mpi_neighbours(western_boundary) == int_missing)  wing_depth_w = 0
    IF(adv_mpi_neighbours(eastern_boundary) == int_missing) wing_depth_e = 0
    IF(adv_mpi_neighbours(southern_boundary) == int_missing) wing_depth_s = 0
    IF(adv_mpi_neighbours(northern_boundary) == int_missing) wing_depth_n = 0
    

    CALL msg('MPI topology adter smpi_reset_periodic_topology:')
    CALL report_neighbours()

  END SUBROUTINE smpi_reset_periodic_topology

  ! ************************************************************************************

  SUBROUTINE report_neighbours()

       CALL msg('my_x_coord:', my_x_coord)
       CALL msg('my_y_coord:', my_y_coord)
       CALL msg('Neighbours:')
       CALL msg('-> east:', adv_mpi_neighbours(eastern_boundary))
       CALL msg('-> west:', adv_mpi_neighbours(western_boundary))
       CALL msg('-> south:', adv_mpi_neighbours(southern_boundary))
       CALL msg('-> north:', adv_mpi_neighbours(northern_boundary))

  END SUBROUTINE report_neighbours

  ! ************************************************************************************

  SUBROUTINE smpi_write_grads_field_parallel(grid_data, npoints, irec, iunit)
    !
    ! Write a single field to a grads file using mpiio
    ! File has to be opened with mpiio routines
    !
    IMPLICIT NONE
    REAL, DIMENSION(:), INTENT(in) :: grid_data
    INTEGER, INTENT(in) :: npoints, irec, iunit
#ifdef SILAM_MPI    
    INTEGER :: ierr, dtypesize, local_datasize, fs_global
    INTEGER(kind=MPI_OFFSET_KIND) :: byteoffset
    CHARACTER(len=*), PARAMETER :: sub_name = 'smpi_write_grads_field_parallel'

    REAL(kind=r4k), ALLOCATABLE, DIMENSION(:) :: write_buffer

    local_datasize = output_decomposition%size_x * output_decomposition%size_y
    fs_global = output_decomposition%globalNx * output_decomposition%globalNy

    IF (fu_fails(local_datasize == npoints, 'Datasizes do not match!', sub_name)) RETURN

    ! Compute the byte offset for this 'record' allowing for offsets bigger than fits
    ! into integer*4.
    CALL mpi_type_size(MPI_REAL, dtypesize, ierr)
    if (smpi_error(ierr, "after mpi_type_size", sub_name)) return
    byteoffset = int(irec-1, MPI_OFFSET_KIND) &
         & * int(fs_global, MPI_OFFSET_KIND) * int(dtypesize, MPI_OFFSET_KIND)

    CALL mpi_file_set_view(iunit, byteoffset, MPI_REAL, &
         & output_decomposition%subarray_mpi_type, &
         & 'native', MPI_INFO_NULL, ierr)
    if (smpi_error(ierr, "after mpi_file_set_view", sub_name)) return
    
#ifdef DOUBLE_PRECISION
    ! Copy the data to single precision buffer
    allocate(write_buffer(size(grid_data)))
    write_buffer = REAL(grid_data, kind=r4k)
    CALL mpi_file_write_all(iunit, write_buffer, local_datasize, MPI_REAL, &
         & MPI_STATUS_IGNORE, ierr)
    IF (smpi_error(ierr, "after mpi_file_write_all double",sub_name)) return
#else
    CALL mpi_file_write_all(iunit, grid_data, local_datasize, MPI_REAL, &
          & MPI_STATUS_IGNORE, ierr)
    IF (smpi_error(ierr, "after mpi_file_write_all",sub_name)) return
#endif

#endif 

  END SUBROUTINE smpi_write_grads_field_parallel

  ! ************************************************************************************

  SUBROUTINE smpi_write_grads_fieldset_parallel(grid_data, npoints, nflds, irec, iunit)
    !$ use omp_lib
    !
    ! Write a set of fields to a grads file using mpiio.  The file has to be opened with
    ! the mpiio routines. The values must be given in single precision (MPI_REAL); the
    ! conversion is probably better done before.
    !
    IMPLICIT NONE
    REAL(r4k), DIMENSION(:,:), INTENT(in) :: grid_data
    INTEGER, INTENT(in) :: npoints, nflds, irec, iunit
    character(len=*), parameter :: sub_name = 'smpi_write_grads_fieldset_parallel'
#ifdef SILAM_MPI    
    INTEGER :: ierr, dtypesize, local_datasize, fs_global, tmp, reslen
    INTEGER(kind=8) :: local_datasize8
    INTEGER(kind=MPI_OFFSET_KIND) :: byteoffset
    integer, save :: dtype_fieldset = MPI_UNDEFINED, nflds_prev = 0
    integer :: fieldsize_decomp
    character(len=200) :: errstr
    integer, dimension(3) :: sizes, starts, subsizes
    integer :: ndims = 3

    ! The static variables make this not threadsafe.
    !$if(fu_fails(.not. omp_in_parallel(), 'This subroutine is not OpenMP compatible', sub_name)) return
    
    fieldsize_decomp = output_decomposition%size_x * output_decomposition%size_y
    fs_global = output_decomposition%globalNx * output_decomposition%globalNy
    local_datasize8 = 1
    local_datasize8 = local_datasize8 * fieldsize_decomp * nflds

    IF (fu_fails(  local_datasize8 < MAX_INT32,'(local_datasize < MAX_INT32) failed', sub_name)) RETURN
    local_datasize = local_datasize8

    IF (fu_fails(fieldsize_decomp == npoints, 'Fieldsizes do not match!', sub_name)) RETURN
    if (fu_fails(size(grid_data,2) == nflds, 'Number of fields does not match', sub_name)) return

    ! The MPI datatype for writing is cached, create if needed.
    if (nflds_prev /= nflds) then
      if (dtype_fieldset /= MPI_UNDEFINED) then
        call mpi_type_free(dtype_fieldset, ierr)
        if (smpi_error(ierr, "after mpi_type_free", sub_name)) return
      end if
      ! Create a datatype corresponding to catenated sub-arrays.
      ! call mpi_type_contiguous(nflds, output_decomposition%subarray_mpi_type, dtype_fieldset, ierr)
      ! the above  should also work but now we're using what is below...

      sizes = (/output_decomposition%globalNx, output_decomposition%globalNy, nflds/)
      subsizes = (/output_decomposition%size_x, output_decomposition%size_y, nflds/)
      starts = (/output_decomposition%offset_x, output_decomposition%offset_y, 0/)

      call mpi_type_create_subarray(ndims, sizes, subsizes, starts, &
                                  & MPI_ORDER_FORTRAN, MPI_REAL, dtype_fieldset, ierr)
      if (smpi_error(ierr, "after mpi_type_create_subarray", sub_name)) return

      call mpi_type_commit(dtype_fieldset, ierr)
      if (smpi_error(ierr, "after mpi_type_commit", sub_name)) return
      nflds_prev = nflds
    end if

    ! Compute the byte offset for this 'record' allowing for offsets bigger than fits
    ! into integer*4.
    CALL mpi_type_size(MPI_REAL, dtypesize, ierr)
    if (smpi_error(ierr, "after mpi_type_size", sub_name)) return
    byteoffset = int(irec-1, MPI_OFFSET_KIND) &
         & * int(fs_global, MPI_OFFSET_KIND) * int(dtypesize, MPI_OFFSET_KIND)
   
    CALL mpi_file_set_view(iunit, byteoffset, MPI_REAL, &
         & dtype_fieldset, &
         & 'native', MPI_INFO_NULL, ierr)
    if (smpi_error(ierr, "after mpi_file_set_view", sub_name)) return
    
    call mpi_file_write_all(iunit, grid_data, local_datasize, MPI_REAL, MPI_STATUS_IGNORE, ierr)
    if (smpi_error(ierr, "mpi_file_write_all", sub_name)) return
#else    
    call set_error('MPI subroutine called in no-mpi version', sub_name)
#endif    

  END SUBROUTINE smpi_write_grads_fieldset_parallel

  !************************************************************************************

  SUBROUTINE smpi_read_grads_fieldset_parallel(grid_data, fieldsize, iUnit, iFirstFld, nflds)
    !$ use omp_lib
    !
    ! Read a set of fields to a grads file using mpiio.  The file has to be opened with
    ! the mpiio routines. The values must be given in single precision (MPI_REAL); the
    ! conversion is probably better done before.
    !
    IMPLICIT NONE
    REAL(r4k), DIMENSION(:,:), INTENT(in) :: grid_data
    INTEGER, INTENT(in) :: fieldsize, nflds, iUnit, iFirstFld
    CHARACTER(len=*), PARAMETER :: sub_name = 'smpi_read_grads_fieldset_parallel'
#ifdef SILAM_MPI    
    INTEGER :: ierr, iTmp, reslen, dtypesize, local_datasize, fs_global
    integer(kind=8) :: local_datasize8
    INTEGER(kind=MPI_OFFSET_KIND) :: byteoffset  !!! This is 8 in both MPICH and OpenMPI
    integer  :: dtype_fieldset 
    integer(kind=8) :: fieldsize_decomp
    character(len=200) :: errstr
    integer, dimension(3) :: sizes, starts, subsizes
    integer :: ndims = 3

    ! The static variables make this not threadsafe.
    !$ if(fu_fails(.not. omp_in_parallel(), 'This subroutine is not OpenMP compatible', sub_name)) return
   
    fieldsize_decomp = output_decomposition%size_x * output_decomposition%size_y
    fs_global = output_decomposition%globalNx * output_decomposition%globalNy

    local_datasize8 = 1
    local_datasize8 = local_datasize8 * fieldsize_decomp * nflds

    IF (fu_fails(local_datasize8  < MAX_INT32,'(local_datasize < MAX_INT32) failed', sub_name)) RETURN
    local_datasize = local_datasize8

    IF (fu_fails(fieldsize_decomp == fieldsize, 'Fieldsizes do not match!', sub_name)) RETURN
    if (fu_fails(size(grid_data,2) >= nflds, 'Number of fields in buffer is not sufficient', sub_name)) return

    ! Create a datatype corresponding to catenated sub-arrays.
    ! call mpi_type_contiguous(nflds, output_decomposition%subarray_mpi_type, dtype_fieldset, ierr)
    ! the above  should also work but now we're using what is below...

    sizes = (/output_decomposition%globalNx, output_decomposition%globalNy, nflds/)
    subsizes = (/output_decomposition%size_x, output_decomposition%size_y, nflds/)
    starts = (/output_decomposition%offset_x, output_decomposition%offset_y, 0/)

    call mpi_type_create_subarray(ndims, sizes, subsizes, starts, &
                                & MPI_ORDER_FORTRAN, MPI_REAL, dtype_fieldset, ierr)

    if (fu_fails(ierr == MPI_SUCCESS, 'Failed to create fieldset datatype', sub_name)) return
    call mpi_type_commit(dtype_fieldset, ierr)
    if (fu_fails(ierr == MPI_SUCCESS, 'Failed to commit fieldset datatype', sub_name)) return

    ! Compute the byte offset for this 'record' allowing for offsets bigger than fits
    ! into integer*4.
    CALL mpi_type_size(MPI_REAL, dtypesize, ierr)
    IF (fu_fails(ierr == MPI_SUCCESS, 'MPI_TYPE_SIZE failed', sub_name)) then
      call mpi_error_string(ierr, errstr, reslen, itmp)
      call msg(errstr)
      return
    endif

    byteoffset = int(iFirstFld-1, MPI_OFFSET_KIND) &
         & * int(fs_global, MPI_OFFSET_KIND) * int(dtypesize, MPI_OFFSET_KIND)
   
    CALL mpi_file_set_view(iunit, byteoffset, MPI_REAL, dtype_fieldset, &
         & 'native', MPI_INFO_NULL, ierr)
    IF (fu_fails(ierr == MPI_SUCCESS, 'MPI_FILE_SET_VIEW failed', sub_name)) then
      call mpi_error_string(ierr, errstr, reslen, itmp)
      call msg(errstr)
      return
    endif
    call mpi_file_read_all(iunit, grid_data, local_datasize, MPI_REAL, MPI_STATUS_IGNORE, ierr)
    IF (ierr /= MPI_SUCCESS) THEN
      call msg("local_datasize,",int(local_datasize/1000), int(mod(local_datasize,1000)))
      call mpi_error_string(ierr, errstr, reslen, itmp)
      call msg(errstr)
      call set_error('mpi_file_read_all failed', sub_name)
      return
    endif

    call mpi_type_free(dtype_fieldset, ierr)
    if (fu_fails(ierr == MPI_SUCCESS, 'Failed to free fieldset datatype', sub_name)) return

#else    
    call set_error('MPI subroutine called in no-mpi version', sub_name)
#endif    

  END SUBROUTINE smpi_read_grads_fieldset_parallel

  !************************************************************************************

#ifdef SILAM_MPI  
  SUBROUTINE create_subarray_datatype
    !
    ! This routine creates the datatypes that are needed when storing the parts
    ! of the output grid to a file using MPIIO.
    !
    ! This routine is internal for this module and should be called after the
    ! domain decomposition is done.
    !
    INTEGER :: ierr
    INTEGER :: ndims = 2
    INTEGER, DIMENSION(2) :: sizes, subsizes, starts

    sizes = (/output_decomposition%globalNx,output_decomposition%globalNy/)
    subsizes = (/output_decomposition%size_x, output_decomposition%size_y/)
    starts = (/output_decomposition%offset_x, output_decomposition%offset_y/)

    ! Note that here we have to use the datatype used in the FILE, not in memory
    CALL mpi_type_create_subarray(ndims, sizes, subsizes, starts, MPI_ORDER_FORTRAN, MPI_REAL, &
         & output_decomposition%subarray_mpi_type, ierr)
    IF (fu_fails(ierr == MPI_SUCCESS, 'MPI_TYPE_CREATE_SUBARRAY failed', 'create_subarray_datatype')) RETURN
    CALL mpi_type_commit(output_decomposition%subarray_mpi_type, ierr)
    IF (fu_fails(ierr == MPI_SUCCESS, 'MPI_TYPE_COMMIT failed', 'create_subarray_datatype')) RETURN

  END SUBROUTINE create_subarray_datatype
#endif  

  ! ************************************************************************************!
  


  ! ************************************************************************************

  SUBROUTINE smpi_gather_field(localFld, globalFLD, ifStagX, ifStagY)
    !
    ! Gathers globalFLD to master (rank 0)
    ! 
    !
    REAL, dimension (:), pointer :: localFld
    REAL, dimension (:,:), pointer ::  globalFLD
    logical, intent(in) :: ifStagX, ifStagY

#ifdef SILAM_MPI    
    REAL, dimension (:),pointer  :: recvbuf, sendbuf
    REAL, dimension (1), target  :: dummy
    REAL, dimension (:,: ),pointer ::ptr2d

    INTEGER  :: nx, ny, offx, offy, gnx, gny, maxmsgsize, itask, ioff, iTmp
    integer :: addx, addy
    INTEGER :: ierr

    addx = 0
    addy = 0
    if (ifStagX) addx = 1
    if (ifStagY) addy = 1

    call smpi_get_maxsubdomain_size(nx, ny)
    maxmsgsize = (nx+addx)*(ny+addy)+5
    call  smpi_get_decomposition(nx, ny, offx, offy, gnx, gny)
    nx  = nx  + addx 
    gnx = gnx + addx
    ny  = ny  + addy
    gny = gny + addy 
    
    recvbuf => dummy
    if (smpi_global_rank == 0) then ! Distribute our local stuff, prepare receive buffer
       ptr2d(1:nx,1:ny) =>  localFld(1:nx*ny)
       globalFLD(offx+1:offx+nx, offy+1:offy+ny) = ptr2d(1:nx,1:ny)
       recvbuf => fu_work_array(smpi_global_tasks*maxmsgsize)
    endif
    sendbuf => fu_work_array(maxmsgsize)
    sendbuf(1:5)=(/real(smpi_global_rank),real(nx), real(ny), real(offx), real(offy)/)
    sendbuf(6:5+nx*ny) = localFld(1:nx*ny)

    call msg("Gathering"+fu_str(smpi_global_rank),(/size(sendbuf), size(recvbuf), maxmsgsize/))
   iTmp=0
    !   do while(iTmp == 0)
    !           call sleep(10)
    !           exit
!
!    enddo
    call MPI_Gather(sendbuf, maxmsgsize, smpi_real_type,&
                  & recvbuf, maxmsgsize, smpi_real_type, 0,  smpi_adv_comm, ierr)
    IF (fu_fails(ierr == MPI_SUCCESS, 'MPI_Gather failed', 'smpi_gather_field')) RETURN

    call free_work_array(sendbuf) ! no need anymore

    !get stuff
    if (smpi_global_rank == 0) then
       do iTask = 1,x_divisions*y_divisions-1
            ioff=iTask*maxmsgsize
            nx   = int(recvbuf(ioff+2))
            ny   = int(recvbuf(ioff+3))
            offx = int(recvbuf(ioff+4))
            offy = int(recvbuf(ioff+5))
            call msg("Got from rank:"+fu_str(int(recvbuf(ioff+1)))+", (/nx, ny, offx, offy/)", (/nx, ny, offx, offy/))
            ptr2d(1:nx,1:ny) =>  recvbuf(ioff+6:ioff+5+nx*ny)
            globalFLD(offx+1:offx+nx, offy+1:offy+ny) = ptr2d(1:nx,1:ny)
       enddo
       call free_work_array(recvbuf)
    endif
#else    
    call set_error("Should not be called ever", "smpi_gather_field")
#endif    
  end SUBROUTINE smpi_gather_field

  !***************************************************************************************

  SUBROUTINE smpi_bcast_aray(array, len, comm,rank)

    REAL, dimension (:), intent(inout) :: array
    INTEGER, intent(in)  :: len, comm,rank
    INTEGER :: ierr
#ifdef SILAM_MPI
    call MPI_Bcast( array,  len, smpi_real_type, rank,  comm, ierr)
    IF (fu_fails(ierr == MPI_SUCCESS, 'MPI_Bcast failed', 'smpi_bcast_aray')) RETURN
#else
    call set_error("Should not be called ever", "smpi_bcast_aray")
#endif

  end SUBROUTINE smpi_bcast_aray

  !***************************************************************************************
  SUBROUTINE smpi_bcast_int8_aray(array, len, comm, rank)

     integer(kind=8), dimension (:), intent(inout) :: array
     INTEGER, intent(in)  :: len, comm, rank
    INTEGER :: ierr
#ifdef SILAM_MPI
     call MPI_Bcast( array,  len, MPI_INTEGER8, rank,  comm, ierr)
    IF (fu_fails(ierr == MPI_SUCCESS, '', '')) RETURN
#else
      call set_error("Should not be called ever", "smpi_bcast_int_aray")
#endif

  end SUBROUTINE smpi_bcast_int8_aray

  !***************************************************************************************

  
  SUBROUTINE smpi_bcast_int_aray(array, len, comm, rank)

     integer, dimension (:), intent(inout) :: array
     INTEGER, intent(in)  :: len, comm, rank
    INTEGER :: ierr
#ifdef SILAM_MPI
     call MPI_Bcast( array,  len, MPI_INTEGER, rank,  comm, ierr)
    IF (fu_fails(ierr == MPI_SUCCESS, '', '')) RETURN
#else
      call set_error("Should not be called ever", "smpi_bcast_int_aray")
#endif

  end SUBROUTINE smpi_bcast_int_aray

  !***************************************************************************************

  SUBROUTINE smpi_allgather_char(myarr, mysize, fullarr)

    character, dimension (:), pointer :: myarr, fullarr
    INTEGER, intent(in)  :: mysize
    INTEGER :: ierr
#ifdef SILAM_MPI    
!       PRINT *, "smpi_allgather_char", smpi_global_rank, mysize, size(myarr), size(fullarr)
    call MPI_Allgather( myarr, mysize, MPI_CHARACTER, &
                       & fullarr, mysize, MPI_CHARACTER, smpi_adv_comm, ierr)
    IF (fu_fails(ierr == MPI_SUCCESS, '', '')) RETURN
!    PRINT *, "Done", smpi_global_rank, ierr
#else    
    call set_error("Should not be called ever", "smpi_allgather_char")
#endif    
  end SUBROUTINE smpi_allgather_char

  !***************************************************************************************

  SUBROUTINE smpi_bcast_char(myarr, mysize, rank)

    character, dimension (:), intent(inout) :: myarr
    INTEGER(kind=8), intent(in)  :: mysize
    integer, intent(in) :: rank
    INTEGER(kind=8) ::  i, j
    INTEGER(kind=8) ::  size_new
    INTEGER :: ierr, cnt
#ifdef SILAM_MPI
         ! Exchaneg by MAX_INT32 chunks
         do i = 0,mysize,MAX_INT32  
            cnt = int(min((i+1)*MAX_INT32, mysize) - i)
            !           call msg("Broadcasting chars", cnt, MAX_INT32)
            call MPI_bcast( myarr(i+1), cnt, MPI_CHARACTER, rank, smpi_adv_comm, ierr)
            IF (fu_fails(ierr == MPI_SUCCESS, '', '')) RETURN
         enddo
#else    
    call set_error("Should not be called ever", "smpi_bcast_char")
#endif    
  end SUBROUTINE smpi_bcast_char

  !***************************************************************************************

  SUBROUTINE smpi_bcast_string(str)

     CHARACTER (*), intent(inout)  :: str
#ifdef SILAM_MPI    
     INTEGER :: clen,ierr
     clen=len(str)

     call MPI_Bcast( str,  clen, MPI_CHARACTER, 0,  smpi_adv_comm, ierr)
     IF (fu_fails(ierr == MPI_SUCCESS, '', '')) RETURN
#else    
    call set_error("Should not be called ever", "smpi_bcast_string")
#endif    

  end SUBROUTINE smpi_bcast_string
 
  !***************************************************************************************

  SUBROUTINE smpi_send(myarr, mysize, dest_rank)

    real, dimension (:), intent(in) :: myarr
    INTEGER, intent(in)  :: mysize, dest_rank
    INTEGER :: ierr, tag
    
#ifdef SILAM_MPI    
    tag = 0
!    print *, "Sending"
    
!    print *, myarr(1:5), "...", myarr(mysize-5:mysize)
!    PRINT *, "smpi_send", smpi_global_rank, "to", dest_rank, "size", mysize, size(myarr)
    call MPI_Send( myarr, mysize, smpi_real_type, dest_rank, tag, smpi_adv_comm, ierr)
    IF (fu_fails(ierr == MPI_SUCCESS, '', '')) RETURN
!    print *, myarr(1:5), "...", myarr(mysize-5:mysize)
!    PRINT *, "Done", smpi_global_rank, mysize, ierr
#else    
    call set_error("Should not be called ever", "smpi_send")
#endif    
  end SUBROUTINE smpi_send

  !***************************************************************************************

  SUBROUTINE smpi_recv(myarr, mysize, src_rank)

    real, dimension (:), intent(out) :: myarr
    INTEGER, intent(in)  :: mysize, src_rank
#ifdef SILAM_MPI    
    integer, dimension(MPI_STATUS_SIZE) :: stat
    INTEGER :: ierr, tag, rank, gotcount

   
    tag=MPI_ANY_TAG
    rank=MPI_ANY_SOURCE
!    print *, "Receiving"
!    call backtrace()

!    print *, myarr(1:5), "...", myarr(mysize-5:mysize)
!    PRINT *, "smpi_recv", smpi_global_rank, "from", src_rank, "maxsize", mysize, size(myarr)
    call MPI_Recv( myarr, mysize, smpi_real_type, src_rank, tag, smpi_adv_comm, stat, ierr)
    IF (fu_fails(ierr == MPI_SUCCESS, '', '')) RETURN
!    print *, stat(:)

    call MPI_Get_count( stat, smpi_real_type, gotcount,  ierr )
    IF (fu_fails(ierr == MPI_SUCCESS, '', '')) RETURN
!   print *, myarr(1:5), "...", myarr(mysize-5:mysize)
!    PRINT *, "Done", smpi_global_rank, "Got size:", ierr
#else    
    call set_error("Should not be called ever", "smpi_recv")
#endif    
  end SUBROUTINE smpi_recv

  !***************************************************************************************
  
  SUBROUTINE smpi_scatterv_real(sendbuf, counts, offsets, recvbuf, recvcount,  comm)

    real,    dimension (:), intent(in) :: sendbuf
    integer, dimension (:), intent(in) :: counts, offsets
    real, dimension (:), intent(out) :: recvbuf
    INTEGER, intent(in)  :: recvcount, comm
    INTEGER :: ierr
#ifdef SILAM_MPI
    call MPI_Scatterv(sendbuf, counts, offsets, &
             &  smpi_real_type, recvbuf, recvcount, smpi_real_type, 0, comm, ierr)
    IF (fu_fails(ierr == MPI_SUCCESS, '', '')) RETURN
#else    
    call set_error("Should not be called ever", "smpi_scatterv_real")
#endif    
           
   
  end SUBROUTINE smpi_scatterv_real

  !***************************************************************************************

  SUBROUTINE smpi_gatherv_real(sendbuf, sendcount, recvbuf, recvcounts, offsets, comm)

    real,    dimension (:), intent(in) :: sendbuf
    integer, intent(in) :: sendcount
    real, dimension (:), intent(out) :: recvbuf
    integer, dimension (:), intent(in) :: recvcounts, offsets
    INTEGER, intent(in)  ::  comm
    INTEGER :: ierr

#ifdef SILAM_MPI
    call MPI_Gatherv(sendbuf, sendcount, smpi_real_type,&
              &       recvbuf, recvcounts, offsets, &
              &        smpi_real_type, 0, comm, ierr)
    IF (fu_fails(ierr == MPI_SUCCESS, '', '')) RETURN
#else    
    call set_error("Should not be called ever", "smpi_gatherv_real")
#endif    
   
  end SUBROUTINE smpi_gatherv_real
 
  ! ***********************************************************************************
#ifdef SILAM_MPI
 logical function smpi_error(ierr, place, sub)
   INTEGER, INTENT(in) :: ierr
   CHARACTER(len=*), INTENT(in) :: place, sub
   INTEGER :: res_len, err2
   CHARACTER(len=MPI_MAX_ERROR_STRING) :: errormsg
    
    smpi_error = .False.
    IF (ierr /= MPI_SUCCESS) THEN
       CALL mpi_error_string(ierr, errormsg, res_len, err2)
       call set_error("MPI ERROR after "//trim(place)//": "//trim(errormsg), sub)
       smpi_error = .True.
    END IF
   end function smpi_error
#endif    

  !******************************************************************

  SUBROUTINE smpi_exchange_boundaries(direction_rank, buffer, send_count, rcv_count, rcv_off)
    ! This routine exchanges the boundary data between two tasks
    ! Argumnents:
    !   direction_rank(in) - MPI rank of the task that we are exchanging data with
    !   buffer(inout)  - work array that has in the beginning the data to be sent; end of that
    !                    buffer is used as receive buffer
    !   in_count(in)   - number of elements to send
    !   out_count(out) - number of received elements
    !   out_index(out) - starting index of the receive buffer
    IMPLICIT NONE
    INTEGER, INTENT(in) :: direction_rank, send_count, rcv_off
    INTEGER, INTENT(inout) :: rcv_count !Size of buffer on in, received count on out
    REAL, DIMENSION(:), POINTER :: buffer
    CHARACTER(len=*), PARAMETER :: sub_name='smpi_exchange_boundaries'
#ifdef SILAM_MPI
    INTEGER :: ierr, sr_status(MPI_STATUS_SIZE), res_len
    CHARACTER(len=MPI_MAX_ERROR_STRING) :: errormsg
    CALL mpi_sendrecv(buffer, send_count, smpi_real_type, direction_rank, 1, &
         & buffer(1+rcv_off), rcv_count, smpi_real_type, direction_rank, 1, &
         & smpi_adv_comm, sr_status, ierr)
    IF (smpi_error(ierr, "after mpi_sendrecv",sub_name)) return

    CALL mpi_get_count(sr_status, smpi_real_type, rcv_count, ierr)
    IF (smpi_error(ierr, "after mpi_get_count",sub_name)) return
    !PRINT *, 'task', smpi_global_rank, "to", direction_rank, 'send_count', send_count, 'received', rcv_count
#else
    RETURN
#endif    

  END SUBROUTINE smpi_exchange_boundaries

  ! ***********************************************************************************

  SUBROUTINE smpi_exchange_wings(direction_rank, send_buf, send_count, recv_buf, recv_count)
    ! This routine exchanges the boundary data between two tasks
    ! Argumnents:
    !   direction_rank(in) - MPI rank of the task that we are exchanging data with
    !   send_buf(in)       - work array that has the data to be sent
    !   send_count(in)     - number of elements to send
    !   recv_buf(out)      - work array that will hold the received elements
    !   recv_count(inout)  - input value is the size of the receive buffer, on return the
    !                        variable will hold the number of received elements
    IMPLICIT NONE
    INTEGER, INTENT(in) :: direction_rank, send_count
    INTEGER, INTENT(out) :: recv_count
    REAL, DIMENSION(:), POINTER :: send_buf, recv_buf
    CHARACTER(len=*), PARAMETER :: sub_name='smpi_exchange_wings'

#ifdef SILAM_MPI 
    INTEGER :: ierr, received_elements, sr_status(MPI_STATUS_SIZE), res_len
    CHARACTER(len=MPI_MAX_ERROR_STRING) :: errormsg

    ! Exchange implies that the send and receive buffers are of same size
    recv_count = send_count
    CALL mpi_sendrecv(send_buf, send_count, smpi_real_type, direction_rank, 1, &
         &            recv_buf, recv_count, smpi_real_type, direction_rank, 1, &
         &            smpi_adv_comm, sr_status, ierr)
    IF (smpi_error(ierr, "after mpi_sendrecv",sub_name)) return

    CALL mpi_get_count(sr_status, smpi_real_type, recv_count, ierr)
    IF (smpi_error(ierr, "after mpi_get_count",sub_name)) return
#else
    RETURN
#endif    
    
  END SUBROUTINE smpi_exchange_wings

  ! ***********************************************************************************

  SUBROUTINE smpi_open_gradsfile_mpiio_w(fname, filehandle)
    ! Opens a file using mpi_file_open, this is needed for mpiio.
    ! Currently this uses the ADVECTION communicator because the mpiio
    ! version does not use separate io tasks, but all tasks that participate
    ! in the advection computation do also write the data.
    IMPLICIT NONE
    CHARACTER(len=*), INTENT(in) :: fname
    INTEGER, INTENT(inout) :: filehandle
    INTEGER :: ierr, res_len
#ifdef SILAM_MPI 
    CHARACTER(len=MPI_MAX_ERROR_STRING) :: errormsg

    CALL mpi_file_open(smpi_adv_comm, fname, MPI_MODE_WRONLY+MPI_MODE_CREATE, MPI_INFO_NULL, &
         & filehandle, ierr)
    IF (smpi_error(ierr, "after mpi_file_open","smpi_open_gradsfile_mpiio_w")) return
#else
    RETURN
#endif    

  END SUBROUTINE smpi_open_gradsfile_mpiio_w

  ! ***********************************************************************************

  SUBROUTINE smpi_open_gradsfile_mpiio_r(fname, filehandle)
    ! Opens a file using mpi_file_open, this is needed for mpiio.
    ! Currently this uses the ADVECTION communicator because the mpiio
    ! version does not use separate io tasks, but all tasks that participate
    ! in the advection computation do also write the data.
    IMPLICIT NONE
    CHARACTER(len=*), INTENT(in) :: fname
    INTEGER, INTENT(inout) :: filehandle
    INTEGER :: ierr, res_len
#ifdef SILAM_MPI 
    CHARACTER(len=MPI_MAX_ERROR_STRING) :: errormsg

    CALL mpi_file_open(smpi_adv_comm, fname, MPI_MODE_RDONLY, MPI_INFO_NULL, &
         & filehandle, ierr)
    IF (smpi_error(ierr, "after mpi_file_open","smpi_open_gradsfile_mpiio_r")) return
#else
    RETURN
#endif    

  END SUBROUTINE smpi_open_gradsfile_mpiio_r

  ! ***********************************************************************************

  SUBROUTINE smpi_close_gradsfile_mpiio(filehandle)
    ! Closes a file that has been opened using mpi_file_open
    INTEGER, INTENT(inout) :: filehandle
    INTEGER :: ierr
#ifdef SILAM_MPI
    CALL mpi_file_close(filehandle, ierr)
    IF (smpi_error(ierr, "after mpi_file_close","smpi_close_gradsfile_mpiio")) return
#else
    RETURN
#endif    

  END SUBROUTINE smpi_close_gradsfile_mpiio

  ! ***********************************************************************************

  SUBROUTINE smpi_reduce_add(sendbuf, recvbuf, root, comm, success)
    IMPLICIT NONE
    REAL, DIMENSION(:), INTENT(in) :: sendbuf
    REAL, DIMENSION(:), INTENT(out) :: recvbuf
    INTEGER, INTENT(in) :: root, comm
    LOGICAL, INTENT(out) :: success

#ifdef SILAM_MPI    
    INTEGER :: ierr
    CALL mpi_reduce(sendbuf, recvbuf, SIZE(sendbuf), smpi_real_type, MPI_SUM, root, comm, ierr)
    IF (smpi_error(ierr, "","smpi_reduce_add")) return
    success = ierr == MPI_SUCCESS
#else
    success = .FALSE.
#endif    

  END SUBROUTINE smpi_reduce_add

  ! **********************************************************************************
  SUBROUTINE smpi_reduce_max(sendbuf, recvbuf, nvals, root, comm, success)
    IMPLICIT NONE
    REAL, DIMENSION(:), INTENT(in) :: sendbuf
    REAL, DIMENSION(:), INTENT(out) :: recvbuf
    INTEGER, INTENT(in) :: root, comm,nvals
    LOGICAL, INTENT(out) :: success

    INTEGER :: ierr, errstr_len
    CHARACTER(len=500) :: errstr
#ifdef SILAM_MPI
    CALL mpi_reduce(sendbuf, recvbuf, nvals, smpi_real_type, MPI_MAX, root, comm, ierr)
    if (smpi_error(ierr, "","smpi_reduce_max")) return
    success = ierr == MPI_SUCCESS
#else
    recvbuf(1:nvals) = sendbuf(1:nvals)
    success = .True.
#endif    

  END SUBROUTINE smpi_reduce_max

  ! **********************************************************************************
  SUBROUTINE smpi_allreduce_add(sendbuf, recvbuf, comm, success)
    IMPLICIT NONE
    REAL, DIMENSION(:), INTENT(in) :: sendbuf
    REAL, DIMENSION(:), INTENT(out) :: recvbuf
    INTEGER, INTENT(in) :: comm
    LOGICAL, INTENT(out) :: success

    INTEGER :: ierr, errstr_len, nvals
#ifdef SILAM_MPI
    CALL mpi_allreduce(sendbuf, recvbuf, SIZE(sendbuf), smpi_real_type, MPI_SUM, comm, ierr)
    if (smpi_error(ierr, "","smpi_allreduce_add")) return
    success = ierr == MPI_SUCCESS
#else
    nvals = SIZE(sendbuf)
    recvbuf(1:nvals) = sendbuf(1:nvals)
    success = .True.
#endif    

  END SUBROUTINE smpi_allreduce_add

  ! **********************************************************************************

  SUBROUTINE smpi_allreduce_max_int(local_max, global_max, comm)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: local_max
    INTEGER, INTENT(OUT) :: global_max
    INTEGER, INTENT(IN) :: comm
    INTEGER :: ierr
    
#ifdef SILAM_MPI
    CALL mpi_allreduce(local_max, global_max, 1, MPI_INTEGER, MPI_MAX, comm, ierr)
    if (smpi_error(ierr, "","smpi_allreduce_max_int")) return
#else
    global_max = local_max
#endif    
    
  END SUBROUTINE smpi_allreduce_max_int

  ! **********************************************************************************

  SUBROUTINE smpi_advection_barrier()
    !
    ! Barrier for the advection computation tasks
    !
    IMPLICIT NONE
    INTEGER :: ierr
#ifdef SILAM_MPI
    CALL mpi_barrier(smpi_adv_comm, ierr)
    if (smpi_error(ierr, "","smpi_advection_barrier")) return
#else
    return
#endif

  END SUBROUTINE smpi_advection_barrier

  ! **********************************************************************************

  SUBROUTINE smpi_global_barrier()
    !
    ! Barrier for the advection computation tasks
    !
    IMPLICIT NONE
    INTEGER :: ierr
#ifdef SILAM_MPI
    CALL mpi_barrier(MPI_COMM_WORLD, ierr)
    if (smpi_error(ierr, "","smpi_global_barrier")) return
#else
    return
#endif

  END SUBROUTINE smpi_global_barrier

  ! **********************************************************************************

  SUBROUTINE smpi_finalize()
    !
    ! Wrapper call fro MPI_Abort
    !
    IMPLICIT NONE
    INTEGER :: ierr

    ! TODO: add somewhere the calls to release the initialized datatypes and
    ! communicators
#ifdef SILAM_MPI    
    CALL mpi_finalize(ierr)
    if (smpi_error(ierr, "","smpi_finalize")) return

#else
    return
#endif

  END SUBROUTINE smpi_finalize

  ! **********************************************************************************

  SUBROUTINE smpi_init()
    !
    ! Initialize the MPI communication using MPI_FUNNELED mode and
    ! check the number of tasks and own global rank.
    !
    IMPLICIT NONE
    CHARACTER(len=*), PARAMETER :: sub_name='smpi_init'
#ifdef SILAM_MPI    
    INTEGER :: required, provided, ierr

#ifdef PUHTI_BUG
    LOGICAL :: FLAGZ, FLAGO,FLAGI

    ! save
    CALL IEEE_GET_HALTING_MODE(IEEE_DIVIDE_BY_ZERO, FLAGZ)
    CALL IEEE_GET_HALTING_MODE(IEEE_OVERFLOW, FLAGO)
    CALL IEEE_GET_HALTING_MODE(IEEE_INVALID, FLAGI)
    !mask
    CALL IEEE_SET_HALTING_MODE(IEEE_DIVIDE_BY_ZERO, .FALSE.)
    CALL IEEE_SET_HALTING_MODE(IEEE_OVERFLOW, .FALSE.)
    CALL IEEE_SET_HALTING_MODE(IEEE_INVALID, .FALSE.)
    print *, "FPE errors masked"
#endif

    !required = MPI_THREAD_FUNNELED
    required = MPI_THREAD_SINGLE
    
    CALL mpi_init_thread(required, provided, ierr)
#ifdef PUHTI_BUG
   !restore
    CALL IEEE_SET_HALTING_MODE(IEEE_DIVIDE_BY_ZERO, FLAGZ)
    CALL IEEE_SET_HALTING_MODE(IEEE_OVERFLOW, FLAGO)
    CALL IEEE_SET_HALTING_MODE(IEEE_INVALID, FLAGI)
    print *, "FPE errors mask restored"
#endif
    if (smpi_error(ierr, "after mpi_init_thread",sub_name)) return
    IF (provided < required) THEN
       call msg('The MPI library does not support MPI_THREAD_SINGLE')
       return
    END IF

    CALL mpi_comm_size(MPI_COMM_WORLD, smpi_global_tasks, ierr)
    if (smpi_error(ierr, "after mpi_comm_size",sub_name)) return
    CALL mpi_comm_rank(MPI_COMM_WORLD, smpi_global_rank, ierr)
    if (smpi_error(ierr, "after mpi_comm_rank",sub_name)) return
#else 
    smpi_global_tasks = 1
    smpi_global_rank = 0
    smpi_ens_rank = 0 !FIXME Should be set elsewhere
  
#endif

  END SUBROUTINE smpi_init

  ! **********************************************************************************

  SUBROUTINE smpi_abort(status)
    IMPLICIT NONE
    INTEGER, INTENT(in) :: status
    INTEGER :: ierr

#ifdef SILAM_MPI
     CALL mpi_abort(MPI_COMM_WORLD, status, ierr)
#else
    call exit(status)
#endif    

  END SUBROUTINE smpi_abort

END MODULE silam_partitioning
