module tangent_linear

  ! A module for handling tangent linear storage.
  ! 
  ! Usage of adjoint for a non-linear operation requires storing the forward
  ! trajectory. The values need to be stored just before each process, so several steps
  ! may be needed inside each timestep. The values are kept in memory if possible. If not,
  ! the values are flushed to temporary files in directory pointed by the tmpdir element
  ! in the trajectory structure.
  !
  ! Throughout this module, "trajectory" is a series of states of the model along time axis - these
  ! states are stored/retrieved for tangent linear computations. 
  ! Do not mix them with Lagrangian trajectores.
  ! 
  ! The tractories are used as follows:
  ! 
  ! - On initialization, each process requiring linearisation requests storage, which is identified by a tag
  ! - After requests are collected, the storage is allocated
  ! - On each timestep, the pointer to the storage for the current step is requested (get_tla_step)
  ! - Each process requests pointer to its storage (get_tla_point)
  ! - The values are either stored (forward) or consumed (adjoint).
  ! 
  !use globals, only : msg, set_error
  use toolbox, only : fu_fails, msg, set_error, int_missing, error, fnlen, dir_slash, fu_pid, fu_next_free_unit
  use optimisation, only : start_count, stop_count

  !use globals, only : fnlen, error
  implicit none
  
  private
  
  public add_tla_traj
  public alloc_tla_trajs
  public get_tla_step
  public fu_get_tla_point

  public t_tla_trajectory
  public t_tla_step
  
  public test_tl_structs

  public defined
  interface defined
     module procedure fu_defined_tla_step
     module procedure fu_defined_tla_traj
  end interface

  integer, parameter, private :: max_trajs = 8

  type rp_5d
     real, dimension(:,:,:,:,:), pointer :: rp ! (species, z, x, y, step)
  end type rp_5d

  type rp_4d
     real, dimension(:,:,:,:), pointer :: rp ! (species, z, x, y)
  end type rp_4d

  type t_tla_trajectory
     private
     ! the trajectories themselves
     type(rp_5d), dimension(max_trajs) :: trajs
     ! the tag of the requesting process
     integer, dimension(max_trajs) :: tags = int_missing
     ! dimensionality: if < 4, dimensions are made singleton starting from left
     integer, dimension(max_trajs) :: dimensions = int_missing
     integer :: num_trajs = 0
     ! the range of timesteps currently loaded in memory
     integer :: in_mem_first = int_missing, in_mem_last = int_missing
     integer :: num_steps_mem = int_missing, num_steps_total = int_missing
     ! the directory for temporary files if needed
     character(len=fnlen) :: tmpdir='.'
     logical :: defined = .false.
     logical :: allocated = .false.
  end type t_tla_trajectory
  
  type t_tla_step
     private
     type(rp_4d), dimension(max_trajs) :: steps
     integer, dimension(max_trajs) :: tags = int_missing     
     integer :: num_steps = 0
     logical :: defined = .false.
  end type t_tla_step

  ! Default size (in bytes) allowed for the TL trajectory storage in memory
  integer(8), parameter, private :: MAX_TRAJ_SIZE_DEF = 10000000000_8

contains
  
  logical function fu_defined_tla_step(step) result(if_defined)
    implicit none
    type(t_tla_step), intent(in) :: step
    if_defined = step%defined
  end function fu_defined_tla_step

  logical function fu_defined_tla_traj(traj) result(if_defined)
    implicit none
    type(t_tla_trajectory), intent(in) :: traj
    if_defined = traj%defined
  end function fu_defined_tla_traj
  
  !*******************************************************************************
  
  subroutine add_tla_traj(tag, traj, dimensions)
    implicit none
    integer, intent(in) :: tag
    type(t_tla_trajectory), intent(inout) :: traj
    integer, intent(in), optional :: dimensions
    integer :: ind_traj

    do ind_traj = 1, traj%num_trajs
      if (traj%tags(ind_traj) == tag) then
        call msg('Tag requested:', tag)
        call set_error('Tag already requested', 'add_tla_trj')
        return
      end if
    end do
    
    traj%num_trajs = traj%num_trajs + 1
    if (fu_fails(traj%num_trajs <= max_trajs, 'Too many TLA trajectories requested', 'add_tla_traj')) return
    traj%tags(traj%num_trajs) = tag
    if (present(dimensions)) then
      if (fu_fails(dimensions < 5 .and. dimensions >= 0, 'Bad dimensionality', 'add_tla_traj')) return
      traj%dimensions(traj%num_trajs) = dimensions
    else
      traj%dimensions = 4
    end if

    traj%defined = .true.
        
    ! allocation later.

  end subroutine add_tla_traj
  
  !*******************************************************************************
  
  subroutine alloc_tla_trajs(traj, nx, ny, nlevs, nspecies_transp, n_timesteps, max_size_in_mem_opt) 
    implicit none
    type(t_tla_trajectory), intent(inout) :: traj
    integer, intent(in) :: nx, ny, nlevs, nspecies_transp, n_timesteps
    integer(8), intent(in), optional :: max_size_in_mem_opt ! bytes

    integer :: stat, numvals, ind_traj, steps_to_mem, values_per_step
    integer :: dimensions, n1d, n2d, n3d, n4d
    integer(8) :: max_size
    
    if (fu_fails(traj%defined, 'TLA trajectory not defined', 'alloc_tla_trajs')) return
       
    ! first count the size
    values_per_step = 0
    do ind_traj = 1, max_trajs
      nullify(traj%trajs(ind_traj)%rp)
      if (ind_traj > traj%num_trajs) exit
      if (fu_fails(traj%tags(ind_traj) /= int_missing, 'Invalid tag', 'alloc_tla_trajs')) return
      n1d = nspecies_transp
      n2d = nlevs
      n3d = nx
      n4d = ny
      dimensions = traj%dimensions(ind_traj)
      if (dimensions < 4) n1d = 1
      if (dimensions < 3) n2d = 1
      if (dimensions < 2) n3d = 1
      if (dimensions < 1) n4d = 1
      numvals = n1d*n2d*n3d*n4d
      call msg('Init tangent linearization, nbr of values per step, real(tag):', &
             & numvals, real(traj%tags(ind_traj)))
      values_per_step = values_per_step + n1d*n2d*n3d*n4d
    end do
    
    call msg('Total size required per timestep (reals):', values_per_step)
    call msg('Timesteps required:', n_timesteps)
    if (present(max_size_in_mem_opt)) then
      max_size = max_size_in_mem_opt
    else
      max_size = MAX_TRAJ_SIZE_DEF
    end if
    steps_to_mem = max_size / 4  / values_per_step
    if (fu_fails(steps_to_mem > 0, 'Cannot fit a TL step in memory', 'alloc_tla_trajs')) return 
    traj%num_steps_mem = min(steps_to_mem, n_timesteps)
    call msg('Steps stored in memory:', traj%num_steps_mem)
    if (n_timesteps <= steps_to_mem) then
      call msg('-> disk storage not used')
    else
      call msg('-> disk storage required, directory: ' //  trim(traj%tmpdir))
    end if

    do ind_traj = 1, max_trajs
      nullify(traj%trajs(ind_traj)%rp)
      if (ind_traj > traj%num_trajs) exit
      
      n1d = nspecies_transp
      n2d = nlevs
      n3d = nx
      n4d = ny
      dimensions = traj%dimensions(ind_traj)
      if (dimensions < 4) n1d = 1
      if (dimensions < 3) n2d = 1
      if (dimensions < 2) n3d = 1
      if (dimensions < 1) n4d = 1
      allocate(traj%trajs(ind_traj)%rp(n1d, n2d, n3d, n4d, traj%num_steps_mem), stat=stat)
      if (fu_fails(stat == 0, 'Allocate failed', 'alloc_tla_trajs')) return
     
    end do

    traj%allocated = .true.
    traj%in_mem_first = 1
    traj%in_mem_last = traj%num_steps_mem
    traj%num_steps_total = n_timesteps

  end subroutine alloc_tla_trajs
  
  
  !*******************************************************************************

  subroutine handle_storage(traj, ind_step)
    !
    ! Handles the TLA disk io if needed. The data are stored into flat binary files each
    ! containing a chunk small enough to fit in memory. Full files are always written and
    ! read.
    ! It both stores and reads the data - depending on whether we go forward of backward in time
    !
    implicit none
    type(t_tla_trajectory), intent(inout) :: traj
    integer, intent(in) :: ind_step
    
    integer :: in_mem_first, in_mem_last

    call start_count('tla_io')
      
    if (traj%in_mem_first <= ind_step .and. ind_step <= traj%in_mem_last) then
      continue
    end if

    if (ind_step < traj%in_mem_first .and. ind_step == 1) then
      ! restarted forward run -> no disk oper but set limits to the first chunk
      traj%in_mem_first = 1
      traj%in_mem_last = min(traj%num_steps_mem, traj%num_steps_total)

    else if (traj%in_mem_last < ind_step) then
      if (fu_fails(ind_step == traj%in_mem_last+1, 'Bad ind_step', 'handle_storage')) return
      ! forward run -> store prev to disk  
      call traj_to_disk(traj)
      if (error) return
      traj%in_mem_first = ind_step
      traj%in_mem_last = min(traj%in_mem_first + traj%num_steps_mem-1, traj%num_steps_total)

    else if (ind_step < traj%in_mem_first) then
      if (fu_fails(ind_step == traj%in_mem_first-1, 'Bad ind_step', 'handle_storage')) return
      ! adjoint run -> load next from disk
      in_mem_last = ind_step
      in_mem_first = max(in_mem_last - traj%num_steps_mem+1, 1)
      call traj_from_disk(traj, in_mem_first, in_mem_last)
      if (error) return

    end if
    
    call stop_count('tla_io')

  contains
    
    subroutine traj_to_disk(traj)
      implicit none
      type(t_tla_trajectory), intent(in) :: traj

      character(len=fnlen) :: filename
      integer :: file_unit, ind_traj, stat, pid, ind_step
      logical :: exists

      pid = fu_pid()
      filename = get_filename(traj%in_mem_first, traj%in_mem_last)
      file_unit = fu_next_free_unit()
      
      inquire(file=filename, exist=exists)
      if (exists) then
        open(unit=file_unit, file=filename)
        close(unit=file_unit, status='delete')
      end if
      open(file_unit, file=filename, form='unformatted', action='write', iostat=stat)
      if (fu_fails(stat == 0, 'Failed to open TL trajectory file', 'traj_to_disk')) return

      call msg('Writing TL traj to file: ' // trim(filename))
      do ind_traj = 1, traj%num_trajs
        do ind_step = 1, traj%num_steps_mem
          write(file_unit, iostat=stat) traj%trajs(ind_traj)%rp(:,:,:,:,ind_step)
          if (fu_fails(stat == 0, 'Failed write', 'traj_to_disk')) return
        end do
      end do
      close(file_unit)
      
    end subroutine traj_to_disk

    subroutine traj_from_disk(traj, in_mem_first, in_mem_last)
      implicit none
      type(t_tla_trajectory), intent(inout) :: traj
      integer, intent(in) :: in_mem_first, in_mem_last ! to be read

      character(len=fnlen) :: filename
      integer :: file_unit, ind_traj, stat, ind_step

      filename = get_filename(in_mem_first, in_mem_last)

      file_unit = fu_next_free_unit()
      open(file_unit, file=filename, form='unformatted', action='read', iostat=stat)
      if (fu_fails(stat == 0, 'Failed to open TL trajectory file', 'traj_from_disk')) return
      
      call msg('Reading TL traj from file: ' // trim(filename))
      do ind_traj = 1, traj%num_trajs
        do ind_step = 1, traj%num_steps_mem
          read(file_unit, iostat=stat) traj%trajs(ind_traj)%rp(:,:,:,:,ind_step)
          if (fu_fails(stat == 0, 'Failed read', 'traj_from_disk')) return
        end do
      end do
      close(file_unit)
      traj%in_mem_first = in_mem_first
      traj%in_mem_last = in_mem_last
    end subroutine traj_from_disk

    function get_filename(in_mem_first, in_mem_last) result(name)
      implicit none
      integer, intent(in) :: in_mem_first, in_mem_last
      character(len=fnlen) :: name
      
      !name = 'filename'

      write(name, fmt='("TL", I0, "_", I0, "_", I0)') fu_pid(), in_mem_first, in_mem_last
      
    end function get_filename

  end subroutine handle_storage
  
  !*******************************************************************************

  subroutine get_tla_step(traj, step, ind_step)
    !
    ! finds and reads stored data or creates and, if time, stores the new chunk
    !
    implicit none
    type(t_tla_trajectory), intent(inout) :: traj
    type(t_tla_step), intent(out) :: step
    integer, intent(in) :: ind_step ! the model timestep, counted from begin of forward run

    integer :: ind_traj, ind_in_arr

    if (fu_fails(traj%allocated, 'TLA trajectory not allocated', 'get_t_tla_step')) return

    call handle_storage(traj, ind_step)
    if (error) return
    ind_in_arr = ind_step - traj%in_mem_first + 1

    do ind_traj = 1, traj%num_trajs
      step%steps(ind_traj)%rp => traj%trajs(ind_traj)%rp(:,:,:,:,ind_in_arr)
      step%tags(ind_traj) = traj%tags(ind_traj)
    end do
    
    step%num_steps = traj%num_trajs
    step%defined = .true.
    
  end subroutine get_tla_step
  
  
  !*******************************************************************************
  
  function fu_get_tla_point(step, tag) result(array_ptr)
    !
    ! Searches for the TLA point of the given tag. It better exist
    !
    implicit none
    type(t_tla_step), intent(in) :: step
    integer, intent(in) :: tag
    
    real, dimension(:,:,:,:), pointer :: array_ptr

    integer :: ind_traj
    
    if (fu_fails(step%defined, 'Undefined tla_step', 'fu_get_tla_point')) return

    do ind_traj = 1, step%num_steps
      if (step%tags(ind_traj) == tag) then
        array_ptr => step%steps(ind_traj)%rp
        return
      end if
    end do
    
    call msg('Looking for tag:', tag)
    call set_error('Tag not found', 'get_tla_point')
    
  end function fu_get_tla_point
  
  !*******************************************************************************

  subroutine test_tl_structs()
    implicit none

    type(t_tla_trajectory) :: traj
    type(t_tla_step) :: step
    integer, parameter :: nx = 2, ny = 3, nz = 4, nspecies = 3, nsteps = 20
    integer, parameter :: dim1 = 2, dim2=3, dim3=4
    real, dimension(:,:,:,:), pointer :: tlapoint
    integer :: ind_step
    
    call add_tla_traj(dim1, traj, dim1)
    call add_tla_traj(dim2, traj, dim2)
    call add_tla_traj(dim3, traj, dim3)
    
    if (error) return

    call alloc_tla_trajs(traj, nx, ny, nz, nspecies, nsteps, 6000_8)

    if (error) return
    
    do ind_step = 1, nsteps
      print *, 'Forward:', ind_step
      call get_tla_step(traj, step, ind_step)
      if (error) return
      tlapoint => fu_get_tla_point(step, dim1)
      if (fu_fails(associated(tlapoint), 'tlapoint not associated', 'test_tl_structs')) return
      if (error) return

      tlapoint = dim1*ind_step
      !print *, lbound(tlapoint), ubound(tlapoint)

      tlapoint => fu_get_tla_point(step, dim2) 
      tlapoint = dim2*ind_step
      !print *, lbound(tlapoint), ubound(tlapoint)

      tlapoint => fu_get_tla_point(step, dim3) 
      tlapoint = dim3*ind_step
    end do

    do ind_step = nsteps, 1, -1
      print *, 'Backwards:', ind_step
      call get_tla_step(traj, step, ind_step)
      if (error) return

      tlapoint => fu_get_tla_point(step, dim1) 
      if (error) return

      print *, 'dim1', tlapoint == dim1*ind_step

      tlapoint => fu_get_tla_point(step, dim2) 
      print *, 'dim2', tlapoint == dim2*ind_step
      !tlapoint = dim2*ind_step

      tlapoint => fu_get_tla_point(step, dim3) 
      print *, 'dim3', tlapoint == dim3*ind_step
       !tlapoint = dim3*ind_step
    end do

    do ind_step = 1, nsteps
      print *, 'Forward:', ind_step
      call get_tla_step(traj, step, ind_step)
      if (error) return
      tlapoint => fu_get_tla_point(step, dim1)
      if (fu_fails(associated(tlapoint), 'tlapoint not associated', 'test_tl_structs')) return
      if (error) return

      tlapoint = dim1*ind_step+1
      !print *, lbound(tlapoint), ubound(tlapoint)

      tlapoint => fu_get_tla_point(step, dim2) 
      tlapoint = dim2*ind_step+2
      !print *, lbound(tlapoint), ubound(tlapoint)

      tlapoint => fu_get_tla_point(step, dim3) 
      tlapoint = dim3*ind_step+3
    end do

    do ind_step = 1, nsteps
      print *, 'Forward:', ind_step
      call get_tla_step(traj, step, ind_step)
      if (error) return
      tlapoint => fu_get_tla_point(step, dim1)
      if (fu_fails(associated(tlapoint), 'tlapoint not associated', 'test_tl_structs')) return
      if (error) return

      tlapoint = dim1*ind_step+3
      !print *, lbound(tlapoint), ubound(tlapoint)

      tlapoint => fu_get_tla_point(step, dim2) 
      tlapoint = dim2*ind_step+4
      !print *, lbound(tlapoint), ubound(tlapoint)

      tlapoint => fu_get_tla_point(step, dim3) 
      tlapoint = dim3*ind_step+5
    end do

    do ind_step = nsteps, 1, -1
      print *, 'Backwards:', ind_step
      call get_tla_step(traj, step, ind_step)
      if (error) return

      tlapoint => fu_get_tla_point(step, dim1) 
      if (error) return

      print *, 'dim1', tlapoint == dim1*ind_step+3

      tlapoint => fu_get_tla_point(step, dim2) 
      print *, 'dim2', tlapoint == dim2*ind_step+4
      !tlapoint = dim2*ind_step

      tlapoint => fu_get_tla_point(step, dim3) 
      print *, 'dim3', tlapoint == dim3*ind_step+5
       !tlapoint = dim3*ind_step
    end do

  end subroutine test_tl_structs

end module tangent_linear
