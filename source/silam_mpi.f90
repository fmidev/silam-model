MODULE silam_mpi
  !
  ! Low level MPI routines
  !
  ! This module has all parts that need the low level MPI layer, that is,
  ! this module does not use any higher level modules (like globals.mod).
  !

#ifdef SILAM_MPI
  USE mpi
#else

!!!STUBS and checks for NON-MPI version
#ifdef  WITH_PNETCDF
#error "WITH_PNETCDF defined for non-MPI version. forgot to define SILAM_MPI?"
#endif
#define MPI_DOUBLE_PRECISION 0

#endif



  IMPLICIT NONE

  PRIVATE


  ! These routines are public
  PUBLIC smpi_is_mpi_version
  PUBLIC smpi_abort
  PUBLIC smpi_advection_barrier
  PUBLIC smpi_close_gradsfile_mpiio
  PUBLIC smpi_exchange_boundaries
  PUBLIC smpi_exchange_wings
  PUBLIC smpi_init
  PUBLIC smpi_open_gradsfile_mpiio_r
  PUBLIC smpi_open_gradsfile_mpiio_w
  PUBLIC smpi_global_barrier
  PUBLIC smpi_finalize
  PUBLIC smpi_reduce_add
  PUBLIC smpi_reduce_max
  PUBLIC smpi_allreduce_add
  PUBLIC smpi_allreduce_max_int





#ifdef SILAM_MPI   
  INTEGER, PARAMETER, PUBLIC :: smpi_proc_null = MPI_PROC_NULL
  LOGICAL, PARAMETER :: is_mpi_silam_version = .TRUE.

  ! These mpi constants are exported to next level
  PUBLIC MPI_LOGICAL, MPI_LOR, MPI_STATUS_IGNORE, MPI_SUCCESS, MPI_COMM_WORLD, MPI_PROC_NULL, &
       & MPI_UNDEFINED, MPI_OFFSET_KIND, MPI_INFO_NULL, MPI_ORDER_FORTRAN, MPI_REAL, MPI_INTEGER, &
       & MPI_DOUBLE_PRECISION, MPI_CHARACTER, MPI_STATUS_SIZE, MPI_ANY_TAG, MPI_ANY_SOURCE, &
       & MPI_INTEGER8, MPI_MODE_RDONLY
#else
   public mpi_gather
   public mpi_scatter
   public mpi_allreduce
   public mpi_bcast
   public mpi_allgather
   public mpi_alltoallv
  INTEGER, PUBLIC, PARAMETER :: MPI_UNDEFINED = -1, MPI_REAL = -1, &
       & MPI_SUCCESS = -1, MPI_LOGICAL = -1, MPI_LOR = -1, &
       & MPI_COMM_WORLD = -1,  MPI_INFO_NULL = -1
  INTEGER, PUBLIC, PARAMETER :: smpi_proc_null = -1
  LOGICAL, PRIVATE, PARAMETER :: is_mpi_silam_version = .FALSE.
#endif

#ifdef DOUBLE_PRECISION
  INTEGER, PARAMETER, PUBLIC :: smpi_real_type = MPI_DOUBLE_PRECISION
#else
  INTEGER, PARAMETER, PUBLIC :: smpi_real_type = MPI_REAL
#endif

  INTEGER, PUBLIC :: smpi_global_tasks = 1, smpi_global_rank = 0
  INTEGER, PUBLIC :: smpi_adv_tasks = 1, smpi_adv_rank =0
  INTEGER, PUBLIC :: smpi_adv_cart_rank = 0, smpi_io_rank = 0, smpi_ens_rank = 0
  ! MPI communicators:
  INTEGER, PUBLIC :: smpi_io_comm, smpi_ensmember_comm, smpi_enkf_comm 
  INTEGER, PUBLIC :: smpi_adv_comm, smpi_adv_cart_comm, smpi_adv_x_comm, smpi_adv_y_comm 

  ! Global error status variable is set here
  LOGICAL, SAVE, PUBLIC :: error

  


CONTAINS

  ! **********************************************************************************

  FUNCTION smpi_is_mpi_version() RESULT(val)
    !
    ! This routine can be used to check if the binary was compiled
    ! with MPI support
    !
    LOGICAL :: val
    val = is_mpi_silam_version
  END FUNCTION smpi_is_mpi_version

  ! **********************************************************************************

  SUBROUTINE smpi_init()
    !
    ! Initialize the MPI communication using MPI_FUNNELED mode and
    ! check the number of tasks and own global rank.
    !
    IMPLICIT NONE
#ifdef SILAM_MPI    
    INTEGER :: required, provided, ierr

    !required = MPI_THREAD_FUNNELED
    required = MPI_THREAD_SINGLE
    CALL mpi_init_thread(required, provided, ierr)
    IF (provided < required) THEN
       WRITE(0,*) 'The MPI library does not support MPI_THREAD_SINGLE'
       STOP
    END IF

    CALL mpi_comm_size(MPI_COMM_WORLD, smpi_global_tasks, ierr)
    CALL mpi_comm_rank(MPI_COMM_WORLD, smpi_global_rank, ierr)
#else 
    smpi_global_tasks = 1
    smpi_global_rank = 0
    smpi_ens_rank = 0 !FIXME Should be set elsewhere
  
#endif

  END SUBROUTINE smpi_init

  ! **********************************************************************************

  SUBROUTINE smpi_advection_barrier()
    !
    ! Barrier for the advection computation tasks
    !
    IMPLICIT NONE
    INTEGER :: ierr
#ifdef SILAM_MPI
    CALL mpi_barrier(smpi_adv_comm, ierr)
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
#else
    return
#endif

  END SUBROUTINE smpi_finalize

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

  ! ***********************************************************************************

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
#ifdef SILAM_MPI
    INTEGER :: ierr, sr_status(MPI_STATUS_SIZE)
    CALL mpi_sendrecv(buffer, send_count, smpi_real_type, direction_rank, 1, &
         & buffer(1+rcv_off), rcv_count, smpi_real_type, direction_rank, 1, &
         & smpi_adv_comm, sr_status, ierr)

    CALL mpi_get_count(sr_status, smpi_real_type, rcv_count, ierr)

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

#ifdef SILAM_MPI 
    INTEGER :: ierr, received_elements, sr_status(MPI_STATUS_SIZE)

    ! Exchange implies that the send and receive buffers are of same size
    recv_count = send_count
    CALL mpi_sendrecv(send_buf, send_count, smpi_real_type, direction_rank, 1, &
         &            recv_buf, recv_count, smpi_real_type, direction_rank, 1, &
         &            smpi_adv_comm, sr_status, ierr)

    CALL mpi_get_count(sr_status, smpi_real_type, recv_count, ierr)
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
    IF (ierr /= MPI_SUCCESS) THEN
       WRITE(0,'(A,A)') 'Failed to open file for mpi write: ', fname
       CALL mpi_error_string(ierr, errormsg, res_len, ierr)
       PRINT *, errormsg
       CALL smpi_abort(1)
    END IF
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
    IF (ierr /= MPI_SUCCESS) THEN
       WRITE(0,'(A,A)') 'Failed to open file for mpi read: ', fname
       CALL mpi_error_string(ierr, errormsg, res_len, ierr)
       PRINT *, errormsg
       CALL smpi_abort(1)
    END IF
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
    IF (ierr /= MPI_SUCCESS) THEN
       WRITE(0,*) 'Failed to close file'
       CALL smpi_abort(1)
    END IF
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
    INTEGER :: ierr, errstr_len
    CHARACTER(len=500) :: errstr
    CALL mpi_reduce(sendbuf, recvbuf, SIZE(sendbuf), smpi_real_type, MPI_SUM, root, comm, ierr)
    success = ierr == MPI_SUCCESS
    IF (.NOT. success) THEN
       CALL mpi_error_string(INT(ierr), errstr, errstr_len, ierr)
       PRINT *, 'MPI ERROR:', errstr
    END IF
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
    success = ierr == MPI_SUCCESS
    IF (.NOT. success) THEN
       CALL mpi_error_string(INT(ierr), errstr, errstr_len, ierr)
       PRINT *, 'MPI ERROR:', errstr
    END IF
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
    CHARACTER(len=500) :: errstr
#ifdef SILAM_MPI
    CALL mpi_allreduce(sendbuf, recvbuf, SIZE(sendbuf), smpi_real_type, MPI_SUM, comm, ierr)
    success = ierr == MPI_SUCCESS
    IF (.NOT. success) THEN
       CALL mpi_error_string(INT(ierr), errstr, errstr_len, ierr)
       PRINT *, 'MPI ERROR:', errstr
    END IF
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
#else
    global_max = local_max
#endif    
    
  END SUBROUTINE smpi_allreduce_max_int

#ifndef SILAM_MPI 
  subroutine mpi_gather(sendbuf, sendcount, sendtype, recvbuf, recvcount, recvtype, root, comm, ierr)
    implicit none
    real  :: sendbuf, recvbuf
    integer :: sendcount, sendtype, recvcount, recvtype, root, comm, ierr
    ierr = 9999 
  end subroutine mpi_gather

  subroutine mpi_scatter(sendbuf, sendcount, sendtype, recvbuf, recvcount, recvtype, root, comm, ierr)
    implicit none
    real  :: sendbuf, recvbuf
    integer :: sendcount, sendtype, recvcount, recvtype, root, comm, ierr
    ierr = 9999
  end subroutine mpi_scatter

  subroutine mpi_allreduce(sendbuf, recvbuf, count, datatype, op, comm, ierr)
    implicit none
    logical :: sendbuf, recvbuf
    integer :: count, datatype, op, comm, ierr
    ierr = 9999
  end subroutine mpi_allreduce

  subroutine mpi_bcast(buffer, count, datatype, root, comm, ierr)
    implicit none
    real :: buffer
    integer :: count, datatype, root, comm, ierr
    ierr = 99990
  end subroutine mpi_bcast

  subroutine mpi_allgather(sendbuf, sendcount, sendtype, recvbuf, recvcount, recvtype, comm, ierr)
    implicit none
    real :: sendbuf, recvbuf
    integer :: sendcount, sendtype, recvcount, recvtype, comm, ierr
    ierr = 9999
  end subroutine mpi_allgather

  subroutine mpi_alltoallv(sendbuf, sendcounts, sdispls, sendtype, &
                         & recvbuf, recvcounts, rdispls, recvtype, &
                         & comm, ierror)
    implicit none
    real :: sendbuf, recvbuf
    integer :: sendcounts, sdispls, sendtype, recvcounts, rdispls, recvtype, comm, ierror
    ierror = 9999
  end subroutine mpi_alltoallv

#endif

! ***********************************************************************************


END MODULE silam_mpi
