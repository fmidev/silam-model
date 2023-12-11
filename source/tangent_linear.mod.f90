module tangent_linear

  ! A module for handling tangent linear storage.
  ! 
  ! Usage of adjoint for a non-linear operation requires storing the forward
  ! trajectory. The values need to be stored just before each process, so several steps
  ! may be needed inside each timestep. The values are kept in memory if possible. If not,
  ! the values are flushed to temporary files in directory pointed by the tmpdir element
  ! in the trajectory structure.
  !
  ! Throughout this module, "trajectory" is a series of states of the model along time axis - 
  ! these states are stored in forward and retrieved for tangent linear computations. 
  ! Do not mix them with Lagrangian trajectores.
  ! 
  ! The tractories are used as follows:
  ! 
  ! - On initialization, each process requiring linearisation requests storage, identified by a tag
  ! - After requests are collected, the storage is allocated
  ! - On each timestep, the pointer to the storage for the current step is requested (get_tla_step)
  ! - Each process requests pointer to its storage (get_tla_point)
  ! - The values are either stored (forward) or consumed (adjoint).
  ! 
  !use globals, only : msg, set_error
  use toolbox, only : fu_fails, msg, set_error, int_missing, error, fnlen, dir_slash, fu_pid, &
           & fu_next_free_unit, proc_ID_string, smpi_adv_comm, fu_str, fu_index, smpi_adv_rank,&
           & create_directory_tree, F_NAN, ooops
  use optimisation, only : start_count, stop_count
  use silam_partitioning, only : smpi_allreduce_max_int, smpi_advection_barrier
  !$ use omp_lib

  !use globals, only : fnlen, error
  implicit none
  
  private
  
  public add_tla_traj
  public alloc_tla_trajs
  public prepare_tla_step

  public fu_get_tla_point
  public fu_get_tla_column
  public fu_get_tla_map

  public set_store_forward  !!Set store-forward flag
  public if_store_forward !!get store-forward flag

  public t_tla_trajectory
  
  public test_tl_structs

  public defined
  interface defined
     module procedure fu_defined_tla_traj
  end interface

  integer, parameter, private :: max_trajs = 8

  type rp_5d
     real, dimension(:,:,:,:,:), allocatable :: rp ! (species, z, x, y, step)
  end type rp_5d

  !!Main storage for TLA with 
  type t_tla_trajectory
     private
     ! the trajectories themselves
     type(rp_5d), dimension(max_trajs) :: trajs
     ! the tag of the requesting process
     integer, dimension(max_trajs) :: tags = int_missing
     ! dimensionality: if < 4, dimensions are made singleton starting from left
     integer, dimension(max_trajs) :: nSpecies = int_missing
     integer :: num_tags = 0
     ! the range of timesteps currently loaded in memory
     integer :: in_mem_first = int_missing, in_mem_last = int_missing
     integer :: num_steps_mem = int_missing, num_steps_total = int_missing
     integer :: in_mem_index_current = int_missing !! Index set by prepare_tla_step
     ! the directory for temporary files if needed
     character(len=fnlen) :: tmpdir='.'
     logical :: store_forward = .FALSE.
     logical :: defined = .false.
     logical :: allocated = .false.
     character(len=fnlen) :: store_directory ! Place where to store the disk stuff (not "here")
  end type t_tla_trajectory
  
  ! Default size (in bytes) allowed for the TL trajectory storage in memory
  integer(8), parameter, private :: MAX_TRAJ_SIZE_DEF = 10000000000_8  !! Per whole MPI family

contains
  
  logical function fu_defined_tla_traj(traj) result(if_defined)
    implicit none
    type(t_tla_trajectory), intent(in) :: traj
    if_defined = traj%defined
  end function fu_defined_tla_traj
  
  !*******************************************************************************
  
  subroutine add_tla_traj(tag, traj, nSpecies) !! Actually, nSpecies is just a dimension 
    implicit none
    integer, intent(in) :: tag
    type(t_tla_trajectory), intent(inout) :: traj
    integer, intent(in) :: nSpecies
    integer :: ind_traj

    if(fu_fails(tag /= int_missing, 'strange TLA tag', 'add_tla_tag'))return
    
    do ind_traj = 1, traj%num_tags
      if (traj%tags(ind_traj) == tag) then
        call msg('Tag requested:', tag)
        call set_error('Tag already requested', 'add_tla_trj')
        return
      end if
    end do
    
    traj%num_tags = traj%num_tags + 1
    if (fu_fails(traj%num_tags <= max_trajs, 'Too many TLA trajectories requested', 'add_tla_traj')) return
    traj%tags(traj%num_tags) = tag
    traj%nSpecies(traj%num_tags)= nSpecies

    traj%defined = .true.
        
    ! allocation later.

  end subroutine add_tla_traj
  
  !*******************************************************************************
  
  subroutine alloc_tla_trajs(traj, nx, ny, nlevs, n_timesteps, max_size_in_mem_opt, store_directory) 
    implicit none
    type(t_tla_trajectory), intent(inout) :: traj
    integer, intent(in) :: nx, ny, nlevs, n_timesteps
    integer(8), intent(in), optional :: max_size_in_mem_opt ! bytes
    character(len=*), intent(in) :: store_directory 

    integer :: stat, numvals, ind_traj, steps_to_mem, values_per_step, values_per_step_global
    integer :: dimensions, n1d, n2d, n3d, n4d
    integer(8) :: max_size
    
    if (fu_fails(traj%defined, 'TLA trajectory not defined', 'alloc_tla_trajs')) return
       
    ! first count the size
    values_per_step = 0
    do ind_traj = 1, traj%num_tags
      if (fu_fails(traj%tags(ind_traj) /= int_missing, 'Invalid tag', 'alloc_tla_trajs')) return
      numvals = nx*ny*nlevs*traj%nSpecies(ind_traj)
      call msg('Init tangent linearization, nbr of values per step, (tag='&
                  & //trim(fu_str(traj%tags(ind_traj)))//'):',  numvals)
      values_per_step = values_per_step + numvals
    end do
    
    !! These needed to synchronized TLA sorage among MPI members, so disk io goes synchronously
    call smpi_allreduce_max_int(values_per_step, values_per_step_global, smpi_adv_comm)  
    call msg('Total size required per timestep (reals), ours, globalmax:', values_per_step, values_per_step_global)
    values_per_step = values_per_step_global
    call msg('Timesteps required:', n_timesteps)
    if (present(max_size_in_mem_opt)) then
      max_size = max_size_in_mem_opt
    else
      max_size = MAX_TRAJ_SIZE_DEF !! Keep total size grows with more workers
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

    do ind_traj = 1,  traj%num_tags
      allocate(traj%trajs(ind_traj)%rp(traj%nSpecies(ind_traj), nlevs, nx, ny, traj%num_steps_mem), stat=stat)
      if (fu_fails(stat == 0, 'Allocate failed', 'alloc_tla_trajs')) return
    end do

    traj%store_directory = store_directory

    traj%allocated = .true.
    traj%num_steps_total = n_timesteps

    call reset_tla_traj(traj, 1, traj%num_steps_mem)

  end subroutine alloc_tla_trajs
  
  
  !*******************************************************************************

  subroutine prepare_tla_step(traj, ind_step)
    !
    ! Handles the TLA disk io if needed  to make sure the needed step is in memory
    ! Shold be called before dealing with a given step
    !
    ! The data are stored into flat binary files each
    ! containing a chunk small enough to fit in memory. Full files are always written and
    ! read.
    ! It both stores and reads the data - depending on whether we go forward of backward in time
    !
    implicit none
    type(t_tla_trajectory), intent(inout) :: traj
    integer, intent(in) :: ind_step
    
    integer :: in_mem_first, in_mem_last

    character(len=*), parameter :: sub_name = 'prepare_tla_step'


    !$ if (omp_in_parallel()) then
    !$    call set_error("Should not be called in parallel!", sub_name)
    !$    return
    !$ endif

    if (traj%in_mem_first <= ind_step .and. ind_step <= traj%in_mem_last) then
      traj%in_mem_index_current = ind_step - traj%in_mem_first + 1
      return
    endif

    if (ind_step < 1 .or. ind_step > traj%num_steps_total) then
        call msg ("ind_step, traj%num_steps_total", ind_step, traj%num_steps_total)
        call set_error("Requested step is beyond the trajectory", sub_name)
        return
    endif

    call start_count('tla_io')
      
    if (fu_fails(traj%allocated, 'TLA trajectory not allocated', sub_name)) return

    if (ind_step < traj%in_mem_first .and. ind_step == 1) then
      ! restarted forward run -> no disk oper but set limits to the first chunk
      in_mem_first = 1
      in_mem_last = min(traj%num_steps_mem, traj%num_steps_total)
      if (traj%store_forward) then
        call set_error("Reset of store_forward trajectories not implemented", sub_name)
        !! This would work only if assimilation window equals
        !! assimilation interval, and would need complete reset of trajectory (incl files)
        !! a check needed above then, and reset here
      else
        call msg("prepare_tla_step start")
      endif
      call reset_tla_traj(traj, in_mem_first, in_mem_last)

    else if (traj%in_mem_last < ind_step) then !! Going forward
      call msg("ind_step, traj%in_mem_first, traj%in_mem_last", (/ind_step, traj%in_mem_first, traj%in_mem_last/))
      call msg("traj%num_steps_mem, traj%num_steps_total", traj%num_steps_mem, traj%num_steps_total)
      if (fu_fails(ind_step == traj%in_mem_last+1, 'Bad ind_step fwd', 'handle_storage')) return
      ! forward run -> store prev to disk  
      call traj_to_disk(traj)
      if (error) return
      in_mem_first = ind_step
      in_mem_last = min(in_mem_first + traj%num_steps_mem-1, traj%num_steps_total)
      !! Load TL from  previous forward if exists
      if (traj%store_forward) then
        call traj_from_disk(traj, in_mem_first, in_mem_last, .TRUE.) !! Can be missing
      else
        call reset_tla_traj(traj, in_mem_first, in_mem_last)
      endif
    else if (ind_step < traj%in_mem_first) then
      call msg("ind_step, traj%in_mem_first, traj%in_mem_last", (/ind_step, traj%in_mem_first, traj%in_mem_last/))
      if (fu_fails(ind_step == traj%in_mem_first-1, 'Bad ind_step adj', 'handle_storage')) return
      if (traj%store_forward) then !!Store prev to disk, normally not needed
        call traj_to_disk(traj)
        if (error) return
      endif
      ! adjoint run -> load next from disk
      in_mem_last = ind_step
      in_mem_first = max(in_mem_last - traj%num_steps_mem+1, 1)
      call traj_from_disk(traj, in_mem_first, in_mem_last, .FALSE.) !!Must be present
      if (error) return
    end if

    traj%in_mem_index_current = ind_step - traj%in_mem_first + 1
    
    call stop_count('tla_io')

  contains
    
    subroutine traj_to_disk(traj)
      implicit none
      type(t_tla_trajectory), intent(in) :: traj

      character(len=fnlen) :: filename
      integer :: file_unit, ind_traj, stat, ind_step
      logical :: exists
      character(len=*), parameter :: sub_name = 'traj_to_disk'

      !Only one creates directory, otherwie the call might get stuck (did at Mahti)
      if (smpi_adv_rank == 0) call create_directory_tree(trim(traj%store_directory))
      call smpi_advection_barrier()

      filename = get_filename(traj%in_mem_first, traj%in_mem_last, traj%store_directory)
      file_unit = fu_next_free_unit()
      
      open(file_unit, file=filename, form='unformatted', action='write',  ACCESS='STREAM',  STATUS='REPLACE',  iostat=stat)
      if (fu_fails(stat == 0, 'Failed to open TL trajectory file', sub_name)) return

      call msg('Writing TL traj to file: ' // trim(filename))
      do ind_traj = 1, traj%num_tags
        do ind_step = 1, traj%num_steps_mem
          write(file_unit, iostat=stat) traj%trajs(ind_traj)%rp(:,:,:,:,ind_step)
          if (fu_fails(stat == 0, 'Failed write', sub_name)) return
        end do
      end do
      close(file_unit)
      
    end subroutine traj_to_disk

    subroutine traj_from_disk(traj, in_mem_first, in_mem_last, ifAllowMissing)
      implicit none
      type(t_tla_trajectory), intent(inout) :: traj
      integer, intent(in) :: in_mem_first, in_mem_last ! to be read
      logical, intent(in) :: ifAllowMissing

      character(len=fnlen) :: filename
      integer :: file_unit, ind_traj, stat, ind_step
      character(len=*), parameter :: sub_name = 'traj_from_disk'

      filename = get_filename(in_mem_first, in_mem_last, traj%store_directory)

      file_unit = fu_next_free_unit()
      open(file_unit, file=filename, form='unformatted', action='read', ACCESS='STREAM',iostat=stat)
      
      if (stat == 0) then
        call msg('Reading TL traj from file: ' // trim(filename))
        do ind_traj = 1, traj%num_tags
          do ind_step = 1, traj%num_steps_mem
            read(file_unit, iostat=stat) traj%trajs(ind_traj)%rp(:,:,:,:,ind_step)
            if (fu_fails(stat == 0, 'Failed read', sub_name)) return
          end do
        end do
        close(file_unit)
        traj%in_mem_first = in_mem_first
        traj%in_mem_last = in_mem_last
      else
        if (fu_fails(ifAllowMissing, 'Failed to open TL trajectory file', sub_name)) return
        call msg('No file present, resetting trajectory (would be from file): ' // trim(filename))
        call reset_tla_traj(traj, in_mem_first, in_mem_last)
      endif
    end subroutine traj_from_disk

    function get_filename(in_mem_first, in_mem_last, dirname) result(name)
      implicit none
      integer, intent(in) :: in_mem_first, in_mem_last
      character(len=fnlen),intent(in) :: dirname
      character(len=fnlen) :: name

      character(len=fnlen) :: strtmp
       
      write(strtmp, fmt='("TL", A, "_", I4.4, "_", I4.4)') trim(proc_ID_string), in_mem_first, in_mem_last
      name = trim(dirname)//dir_slash//trim(strtmp)
      
    end function get_filename

  end subroutine prepare_tla_step


  
  !*******************************************************************************
  
  function fu_get_tla_point(traj, tag, iz, ix, iy) result(array_ptr)
    !
    ! Gets a pointer to TLA vector for a single cell  (1:nSp)
    !
    implicit none
    type(t_tla_trajectory), target, intent(in) :: traj
    integer, intent(in) :: tag, ix, iy, iz
    
    real, dimension(:), pointer :: array_ptr

    integer ::  indStepMem, ind_traj
    character(len = *), parameter :: sub_name = 'fu_get_tla_point'
    
    ind_traj = fu_index(tag, traj%tags)
    if (ind_traj > 0) then
      array_ptr => traj%trajs(ind_traj)%rp(:,iz,ix,iy, traj%in_mem_index_current)
    else
      array_ptr => null()
    endif
    
  end function fu_get_tla_point
  
  !*******************************************************************************
  
  function fu_get_tla_column(traj, tag, ix, iy) result(array_ptr)
    !
    ! Gets a pointer to TLA vector for a single column (1:nsp, 1:nz)
    !
    implicit none
    type(t_tla_trajectory), target, intent(in) :: traj
    integer, intent(in) ::   tag, ix, iy
    
    real, dimension(:,:), pointer :: array_ptr 

    integer ::  indStepMem, ind_traj
    character(len = *), parameter :: sub_name = 'fu_get_tla_column'
    
    ind_traj = fu_index(tag, traj%tags)
    
    if (ind_traj > 0) then
      array_ptr => traj%trajs(ind_traj)%rp(:,:,ix,iy, traj%in_mem_index_current)
    else
      array_ptr => null()
    endif
    
  end function fu_get_tla_column
  !*******************************************************************************
  
  function fu_get_tla_map(traj, tag) result(array_ptr)
    !
    ! Gets a pointer to TLA vector for a 4D thing (nSp,nz,nx,ny)
    !
    implicit none
    type(t_tla_trajectory), target, intent(in) :: traj
    integer, intent(in) :: tag
    
    real, dimension(:,:,:,:), pointer :: array_ptr

    integer ::  indStepMem, ind_traj
    character(len = *), parameter :: sub_name = 'fu_get_tla_map'
    
    ind_traj = fu_index(tag, traj%tags)
    if (ind_traj > 0) then
      array_ptr => traj%trajs(ind_traj)%rp(:,:,:,:, traj%in_mem_index_current)
    else
      array_ptr => null()
    endif
    
  end function fu_get_tla_map
  
  !*******************************************************************************

  subroutine test_tl_structs()
    implicit none

    type(t_tla_trajectory) :: traj
    integer, parameter :: nx = 2, ny = 3, nz = 4,  nsteps = 200
    integer, parameter :: dim1 = 2, dim2=3, dim3=4
    integer :: ind_step, itry
    
    call add_tla_traj(dim1, traj, dim1)
    call add_tla_traj(dim2, traj, dim2)
    call add_tla_traj(dim3, traj, dim3)
    
    if (error) return

    call alloc_tla_trajs(traj, nx, ny, nz, nsteps, 6000_8, '.')

    if (error) return
    
    do itry=1,2
      do ind_step = 1, nsteps
        print *, 'Forward:', ind_step
        call prepare_tla_step(traj, ind_step)
        if (error) return
        call storeTLAbypoint(traj, ind_step, dim1)
        call storeTLAbypoint(traj, ind_step, dim2)
        call storeTLAbypoint(traj, ind_step, dim3)
      end do

      do ind_step = nsteps, 1, -1
        print *, 'Backwards:', ind_step
        call prepare_tla_step(traj, ind_step)
        if (error) return
        call checkTLAbycol(traj, ind_step, dim3)
        call checkTLAbycol(traj, ind_step, dim2)
        call checkTLAbycol(traj, ind_step, dim1)
      end do

    end do
    contains
       !*******************************************************
      
      subroutine storeTLAbypoint(traj, indstep, tag)
          implicit none
          type(t_tla_trajectory), intent(in) :: traj
          integer, intent(in) :: indstep, tag
          character(len = *), parameter :: sub_name = 'storeTLA'
          integer :: ix, iy, iz
          real, dimension(:), pointer :: storevec

          do iy = 1,ny
             do ix = 1,nx
               do iz = 1,nz
                  storevec => fu_get_tla_point(traj, tag, iz, ix, iy)
                  storevec(:) = real(tag + 10*(ix +10*(iy +10*(iz + 1000*indstep))))
               enddo
             enddo
          enddo
      end subroutine storeTLAbypoint

      subroutine checkTLAbycol(traj, indstep, tag)
          implicit none
          type(t_tla_trajectory), intent(in) :: traj
          integer, intent(in) :: indstep, tag
          character(len = *), parameter :: sub_name = 'storeTLA'
          integer :: ix, iy, iz
          real :: val
          real, dimension(:,:), pointer :: storevec

          do iy = 1,ny
             do ix = 1,nx
               storevec => fu_get_tla_column(traj, tag, ix, iy)
               do iz = 1,nz
                 val = real(tag + 10*(ix +10*(iy +10*(iz + 1000*indstep))))
                 if (any( storevec(:, iz) /= val)) then
                    print *, "Failed  traj, istep, tag, ix, iy" , indstep, tag, ix, iy
                    print *, storevec(:, iz)
                    print *, val
                 endif
               enddo
             enddo
          enddo
      end subroutine checkTLAbycol
  end subroutine test_tl_structs

  !*******************************************************
 
  subroutine reset_tla_traj(traj, in_mem_first, in_mem_last)
     implicit none
     type(t_tla_trajectory), target, intent(inout) :: traj
     integer, intent(in) :: in_mem_first, in_mem_last
     integer :: ind_traj

      do ind_traj = 1, traj%num_tags
            traj%trajs(ind_traj)%rp(:,:,:,:,:) = F_NAN
      end do
      traj%in_mem_first = in_mem_first
      traj%in_mem_last = in_mem_last

 end subroutine reset_tla_traj

  !*******************************************************
 
 subroutine set_store_forward(traj, ifStore)
     implicit none
     type(t_tla_trajectory), target, intent(inout) :: traj
     logical, intent(in) :: ifStore
     character(len = *), parameter :: sub_name = 'set_store_forward'

     traj%store_forward = ifStore
 end subroutine set_store_forward

  !*******************************************************

 function if_store_forward(traj) result(ifStore)
     implicit none
     type(t_tla_trajectory), target, intent(in) :: traj
     logical :: ifStore
     character(len = *), parameter :: sub_name = 'if_store_forward'

     ifStore = traj%store_forward
 end function if_store_forward

end module tangent_linear
