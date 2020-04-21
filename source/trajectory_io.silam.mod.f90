MODULE trajectory_io
!
! Current module is made for the needs of the back-compatibility with old
! output format of the results stored as a set of trajectories (TRADOS style)
!
! Here a minimum set of routines for storing these trajectories is provided.
! Trajectories are selected from the main particle cloud. They are put in one
! huge dynamical structure and stored during the simulation followed by 
! writing the output file.
! Searching for the weather data for each particle is an
! expensive procedure, which slows down the program.
!
! BE CAREFUL, DO NOT ORDER TOO MANY TRAJECTORIES.
!
! Author: Mikhail Sofiev mikhail.sofiev@fmi.fi
! Author: Roux aka rostislav.kouznetsov@fmi.fi
!
! All units: SI
!
! Language: FORTRAN-90
!
use lagrange_particles
use md
!use field_buffer


  implicit none

private

  PUBLIC trajectory_input_needs ! Just a list of quantities to be written
  PUBLIC init_trajectory_output ! Allocate memory and fill in the initial values
                                ! open output files
  PUBLIC store_traj_time        ! Drops particle positions to the array
                                ! 
  PUBLIC finalize_trajectory_output     ! Flushes the trajectory array to a given file
  PUBLIC defined

  private dump_trajectory


  INTERFACE defined
    MODULE PROCEDURE fu_trajectory_set_defined
  END INTERFACE
  !-------------------------------------------------------------
  !
  TYPE silam_trajectory_set
    PRIVATE
    INTEGER :: nTraj = 0, max_steps = 0, nSrcs = 0 !sizes
    TYPE(silja_interval) :: time_step = interval_missing
    INTEGER :: traj_part = 500 ! Store every traj_part particle

    !Bookkeeping
    integer :: nFree=0, nReady=0, iTrajMax=0 !Stack indices, Highest index of valid trajectory
    integer, DIMENSION(:), ALLOCATABLE :: ReadyToDumpList, FreeList !Stacks of iTraj indices
     
    ! Dump-related stuff
    integer, DIMENSION(:), ALLOCATABLE :: traj_dumped ! (1:nSrc) Counter for dumped trajectories
    integer,  DIMENSION(:), ALLOCATABLE :: ounit !(1:nSrc) ! File unut
    character(len=fnlen),  DIMENSION(:), ALLOCATABLE :: fname_tmp, fname !(1:nSrc)
    integer, DIMENSION(:), ALLOCATABLE  :: NoTraj_offset ! (1:nSrc) Offset to write number of trajectories
    ! Data
    integer,  DIMENSION(:), ALLOCATABLE ::nbr_of_steps, iParticle, iSrc !nbr_of_traj
    TYPE(silja_time), DIMENSION(:), ALLOCATABLE :: time_start !(1:nbr_of_traj)
    REAL, DIMENSION(:,:,:), ALLOCATABLE :: traj_data !  i_x:i_abl,1:max_nbr_of_steps,1:nbr_of_traj
    type(silja_logical) :: defined = silja_false
  END TYPE silam_trajectory_set

  PUBLIC silam_trajectory_set

  type(silam_trajectory_set), public, parameter :: trajectory_set_missing = &
                    & silam_trajectory_set(0,0,0,interval_missing,0,0,0,0, &
                    & null(),null(),null(),null(),null(),null(), &
                    & null(),null(),null(),null(),null(),null(), silja_false )

!  REAL,DIMENSION(:),POINTER,PRIVATE,SAVE :: abl_h_m_ptr => null(), pasquill_ptr => null()
!  integer, private, save :: ind_abl_h_m = int_missing, ind_pasquill = int_missing,&
!           &  ind_pr_4d = int_missing, ind_height_4d = int_missing, ind_scav_4d = int_missing
   
   integer, dimension(1:7), private, save :: idx_meteo = int_missing ! Indices in
                                                                    ! meteo buffer

   integer, parameter, private :: i_x = 1, i_y = 2, i_pressure = 3, i_height = 4, i_abl_h = 5,&
                        & i_scav = 6, i_stability = 7

!                        
!   integer, parameter, private :: traj_part = 500 ! Number of traj =  NoP/traj_part


  !-------------------------------------------------------------
  !
  ! THE MAIN ARRAY WITH TRAJECTORIES. The number of sources determine
  ! the dimension of the array
  !
  TYPE(silam_trajectory_set), private, target, save:: tr_set = &
         &  trajectory_set_missing 


CONTAINS


  ! ***************************************************************

  SUBROUTINE trajectory_input_needs(q_met_dyn, iQ, q_met_st, iQs)
    !
    ! Returns the list of needed quantities for creating of the TRADOS-like output
    !
    IMPLICIT NONE

    ! Imported parameter
    INTEGER, DIMENSION(:), intent(inout) :: q_met_dyn, q_met_st
    integer, intent(out) :: iQ, iQs

    iQ =  fu_merge_int_arrays((/height_flag, &
                              & pressure_flag, &
                              & abl_height_m_flag, &
                              & pasquill_class_flag, &
                              & int_missing/), q_met_dyn, .false.)
!    iQs = fu_merge_int_arrays((/scavenging_coefficient_flag, &
!                              & int_missing/), q_met_st, .false.)
    iQs = fu_merge_int_arrays((/int_missing/), q_met_st, .false.) 
    !! missing scavenging_coefficient_flag should not crash the run

  END SUBROUTINE trajectory_input_needs


  ! ***************************************************************

  SUBROUTINE init_trajectory_output(tr_set, max_nbr_of_steps, nTraj, every, nSrcs,  &
                            & Template, ini_time, chCaseNm, chSrcNm, &
                                  & pMetBuf, time_step)
  !
  ! Allocates memory for the main structure silam_trajectory_set, opens input
  ! files  and fills a few initial values.
  !
    IMPLICIT NONE

    ! Imported parameters with the intent IN
    TYPE(silam_trajectory_set), intent(out):: tr_set
    type(grads_template), intent(in) :: Template
    character (len=*), intent(in) :: chCaseNm
    character (len=*), dimension(:), intent(in) :: chSrcNm
    type (silja_time), intent(in) :: ini_time
    INTEGER, INTENT(in) :: max_nbr_of_steps,  nTraj, every, nSrcs
    type(silja_interval),  INTENT(in) :: time_step
    type(Tfield_buffer), intent(in) :: pMetBuf

    ! Local variables

    integer, dimension(:), pointer:: q_arr_ptr
    character (len=clen) :: tstr
    integer :: istat, iSrc, unit

    character (len=*), parameter :: sub_name="init_trajectory_output"

    q_arr_ptr => pMetBuf%buffer_quantities
    idx_meteo(i_pressure)  = fu_index(q_arr_ptr, pressure_flag)
    idx_meteo(i_height)    = fu_index(q_arr_ptr, height_flag)
    idx_meteo(i_abl_h)     = fu_index(q_arr_ptr, abl_height_m_flag)
    idx_meteo(i_scav)      = fu_index(q_arr_ptr, scavenging_coefficient_flag)
    idx_meteo(i_stability) = fu_index(q_arr_ptr, pasquill_class_flag)

    if (any(idx_meteo(i_pressure:i_stability)<0)) then
       if (idx_meteo(i_scav) > 0) then  !!! scavenging_coefficient_flag might be missing....
         call msg ("idx_pressure, idx_height, idx_abl_h, idx_scav, idx_stability", &
                  idx_meteo(i_pressure:i_stability))
         call set_error("Missing input for trajectory output",sub_name)
       endif
    endif


    !Init to zero
    tr_set = trajectory_set_missing

    ! Set Basic parameters
    tr_set%time_step = time_step
    tr_set%nSrcs = nSrcs
    tr_set%max_steps = max_nbr_of_steps
    tr_set%traj_part = every
    tr_set%nFree = 0

    allocate(tr_set%traj_dumped(tr_set%nSrcs), tr_set%ounit(tr_set%nSrcs), &
           & tr_set%fname_tmp(tr_set%nSrcs), tr_set%fname(tr_set%nSrcs), &
           & tr_set%NoTraj_offset(tr_set%nSrcs), stat=istat )
    IF(istat /= 0) then
       call msg("Requested: nSrcs, nTraj, max_steps", (/ tr_set%nSrcs, tr_set%nTraj, tr_set%max_steps/))
       call set_error('Cannot allocate trajectory set',sub_name)
       return
    endif

    ! Output fiels and per-source stuff
    tstr=fu_computer_time_string()
    do iSrc = 1, nSrcs
      tr_set%traj_dumped(iSrc) = 0
      tr_set%fname(iSrc) = fu_FNm(template, &
                          & ini_time, & ! ini_time
                          & ini_time, & ! anal_time
                          & zero_interval, &   ! forecast length
                          & chCaseNm, &
                          & chSrcNm(iSrc), '')

      tr_set%fname_tmp(iSrc) = tr_set%fname(iSrc) +'.tmp-'+tstr+'.traj'
      tr_set%fname(iSrc) = tr_set%fname(iSrc) +'.traj'
      unit=fu_next_free_unit()

      open(unit, file = tr_set%fname_tmp(iSrc), action="write", status="replace", iostat = istat)
      if (istat /= 0) then
        call set_error('Failed to open: '//trim(tr_set%fname_tmp(iSrc)),sub_name)
        return
      endif
      tr_set%ounit(iSrc) = unit

      WRITE(unit,'(A,A)')'#Old trajectory output(TRADOS style) written by ', revision_str
      if (nSrcs>1) WRITE(unit,'(A,A)')'#Source name: ', chSrcNm(iSrc)
      WRITE(unit,'(A,A)')'#Case name: ',chCaseNm
      WRITE(unit,'(A)')'# add. quantity:  1: boundary layer height from gound [m]'
      WRITE(unit,'(A)')'# add. quantity:   2: scavenging coefficient [1/s]'
      WRITE(unit,'(A)')'# add. quantity:   3: pasquill stability class'
      istat = ftell(unit) !Save position to write final number of trajs
      tr_set%NoTraj_offset(iSrc) = istat !Save position to write final number of trajs
      !!! WARNING If changing line below, adjust corresponding line in "finalize_trajectory_output"
      WRITE(unit,'(A25)') 'XXXXXXXXXX trajectories' ! Placeholder for trajectory number in file
   enddo

      
   tr_set%defined = silja_true

   ! And adjust main arrays to proper size
   ! Must reserve something, so increase by 20% would be an increas
   call resize_tr_set(tr_set, max(nTraj, 100)) 


  END SUBROUTINE init_trajectory_output

 subroutine resize_tr_set(tr_set, nTraj) 

    ! Allocate only 

    type(silam_trajectory_set), intent(inout) :: tr_set
    integer, intent(in) :: nTraj !New number of trajectories

    type(silam_trajectory_set) :: tr_tmp !Dummy container for temporary arrays
    integer :: iTmp, istat
    character (len=*), parameter :: sub_name="resize_tr_set"
    

    if ( nTraj == tr_set%nTraj)  return

    if (tr_set%iTrajMax > nTraj) then 
        call msg("Active trajectories, new size", tr_set%iTrajMax, nTraj)
        call set_error("Active trajectories will not fit new size",sub_name)
    endif

    if (ntraj == 0) then !Reset
      tr_tmp = trajectory_set_missing
      return
    endif

    if (.not. all((/tr_set%nSrcs, tr_set%max_steps/) > 0)) then
      call set_error("Can resize only initialized silam_trajectory_set",sub_name)
      return
    endif

    ! Allocate stuff

    call msg ("Resizing trajectory set from. to, max_steps, nSrc", &
         & (/tr_set%nTraj, nTraj,  tr_set%max_steps, tr_set%nSrcs/))
    
    allocate(tr_tmp%FreeList(nTraj),&
           & tr_tmp%ReadyToDumpList(nTraj), &
           & tr_tmp%nbr_of_steps(nTraj),&
           & tr_tmp%iParticle(nTraj),&
           & tr_tmp%iSrc(nTraj),&
           & tr_tmp%time_start(nTraj),&
           & tr_tmp%traj_data(size(idx_meteo), tr_set%max_steps, nTraj), & 
                STAT=istat)
    IF(istat /= 0) then
       call msg("Requested: nSrcs, nTraj, max_steps",  (/tr_set%nSrcs, nTraj, tr_set%max_steps/))
       call set_error('Cannot allocate trajectory set',sub_name)
       return
    endif

    if (tr_set%nTraj > 0) then !Not the initial allocation: Copy tr_set to tr_tmp and free tr_set
       !copy valuable content to new
       tr_tmp%ReadyToDumpList(1:tr_set%nTraj) =  tr_set%ReadyToDumpList(1:tr_set%nTraj)  
       tr_tmp%nbr_of_steps(1:tr_set%nTraj)    =  tr_set%nbr_of_steps(1:tr_set%nTraj)
       tr_tmp%iParticle(1:tr_set%nTraj)       =  tr_set%iParticle(1:tr_set%nTraj)
       tr_tmp%iSrc(1:tr_set%nTraj)            =  tr_set%iSrc(1:tr_set%nTraj)
       tr_tmp%time_start(1:tr_set%nTraj)      =  tr_set%time_start(1:tr_set%nTraj)
       tr_tmp%traj_data(1:size(idx_meteo),1:tr_set%max_steps,1:tr_set%nTraj)  & 
                               &=  tr_set%traj_data(1:size(idx_meteo),1:tr_set%max_steps,1:tr_set%nTraj)
       
       !Adjust bookkeeping - put new particle numbers to the bottom of the free stack
       do iTmp=1, nTraj-tr_set%nTraj ! New trajectories
          tr_tmp%FreeList(iTmp) = nTraj + 1 - iTmp !First In Last out
       enddo
       ! And put an old free list on top of the new one
       tr_tmp%FreeList(iTmp:iTmp+tr_set%nFree-1) = tr_set%FreeList(1:tr_set%nFree)
       tr_set%nFree = tr_set%nFree + (nTraj - tr_set%nTraj)
     else
       !Init bookkeeping
       tr_set%nFree = nTraj
       do iTmp = 1,nTraj
         tr_tmp%FreeList(iTmp) = nTraj + 1 - iTmp !First In Last out
       enddo
    endif

    tr_tmp%iParticle(tr_set%nTraj+1:nTraj) =  int_missing ! Void newly allocated trajectories

    ! Move reassign new stuff to tr_set
    CALL MOVE_ALLOC(tr_tmp%ReadyToDumpList, tr_set%ReadyToDumpList)
    CALL MOVE_ALLOC(tr_tmp%FreeList,        tr_set%FreeList)
    CALL MOVE_ALLOC(tr_tmp%nbr_of_steps,    tr_set%nbr_of_steps)
    CALL MOVE_ALLOC(tr_tmp%iParticle,       tr_set%iParticle)
    CALL MOVE_ALLOC(tr_tmp%iSrc,            tr_set%iSrc)
    CALL MOVE_ALLOC(tr_tmp%time_start,      tr_set%time_start)
    CALL MOVE_ALLOC(tr_tmp%traj_data,       tr_set%traj_data)

    tr_set%nTraj = nTraj

    call msg("Resizing rajectory set to "//fu_str(nTraj)// " done!")

  end subroutine resize_tr_set




  ! ***************************************************************

  SUBROUTINE store_traj_time(tr_set, met_buf, lpSet, now)
  !
  ! Drops current particle positions to the trajectory structure.
  ! Dumps trajectories that are ready to files
  !
  !
  IMPLICIT NONE
    
    ! Imported parameters with the intent IN
    type(silam_trajectory_set), intent(inout) :: tr_set
    type(Tfield_buffer), pointer :: met_buf
    type(Tlagrange_particles_set), INTENT(in), target :: lpSet
    TYPE(silja_time), INTENT(in):: now

    ! Local variables
    INTEGER :: iTraj, iStep, iTmp, iPart, nTrajBefore
    REAL :: fTmp
    real, dimension(:,:),  pointer :: lpSetDyn
    integer, dimension(:), pointer :: lpStatus
    character (len=*), parameter :: sub_name="store_traj_time"

    lpSetDyn => lpSet%lpDyn
    lpStatus => lpSet%lpStatus

    if (.not. defined(tr_set)) then
         call set_error("tr_set undefined", sub_name)
         return
    endif
    
    !
    ! Assign (some of) new particles to trajectories
    !
    !    call msg("tr_set%FreeList", tr_set%FreeList(max(1,tr_set%nFree-20):tr_set%nFree))
    
    nTrajBefore =  tr_set%Ntraj - tr_set%nFree
    do iTmp=1, lpset%nNewPart
      if (lpSet%newPartList(iTmp) == int_missing) exit !All new particles accounted
      CALL RANDOM_NUMBER(fTmp) 
      if (fTmp*tr_set%traj_part > 1.) cycle ! No traj for this particle
      
      if  (tr_set%nFree < 1) & ! Enlarge by 25%
            & call  resize_tr_set(tr_set, tr_set%Ntraj * 4 / 3 + 100) 
      
      iTraj = tr_set%FreeList(tr_set%nFree)  ! Get next free trajectory
      tr_set%nFree = tr_set%nFree - 1

      ! Init the trajectory
      tr_set%nbr_of_steps(iTraj) = 0
      iPart = lpSet%newPartList(iTmp)  !New particcle number in lpset
      tr_set%iParticle(iTraj) = iPart
      tr_set%iSrc(iTraj) = mod(lpstatus(iPart), 100)
      tr_set%time_start(iTraj) = now
      tr_set%iTrajMax = max(iTraj,tr_set%iTrajMax)
    enddo
    call msg("store_traj: nNewPart, nTraj, NtrajNew, iFreeTraj", &
     & (/lpset%nNewPart, nTrajBefore, tr_set%Ntraj - tr_set%nFree, tr_set%FreeList(tr_set%nFree)/))

    !
    ! Scan through the trajectories checking for each of them, whether the corresponding particle
    ! is active. Store it if yes. If not, Set the trajectory done
    !
    do iTraj = 1, tr_set%iTrajMax ! Scan only active part of trajectories
      !
      if(tr_set%iParticle(iTraj) == int_missing)cycle !void trajectory

      if(lpStatus(tr_set%iParticle(iTraj)) == int_missing .or.& ! Particle is dead
         & tr_set%nbr_of_steps(iTraj) >= tr_set%max_steps )then !Too long trajectory
        ! trajectory done
        tr_set%nReady = tr_set%nReady + 1
        tr_set%ReadyToDumpList(tr_set%nReady)=iTraj
        if (tr_set%nbr_of_steps(iTraj) >= tr_set%max_steps) then
                call msg("Trajectory ran out of steps. iTraj, iPart", iTraj, tr_set%iParticle(iTraj))
        endif
        cycle
      endif
      
      ! All fine.  Store the time step.
      iStep = tr_set%nbr_of_steps(iTraj) + 1
      iPart = tr_set%iParticle(iTraj)

      tr_set%traj_data(i_x,iStep,iTraj) = lpSetDyn(lp_x, iPart)
      tr_set%traj_data(i_y,iStep,iTraj) = lpSetDyn(lp_y, iPart)
      ! It would be great to get height and pressure properly....

      do iTmp = i_pressure,i_stability 
        if (idx_meteo(iTmp) > 0) then  !! only for non-missing data
         ! Nearest-neighbour
         tr_set%traj_data(iTmp,iStep,iTraj) =  fu_get_value(met_buf, idx_meteo(iTmp),&
                                                & nint(lpSetDyn(lp_x, iPart)), &
                                                & nint(lpSetDyn(lp_y, iPart)), &
                                                & nint(lpSetDyn(lp_z, iPart)), &
                                                & nx_meteo, now)
        else
          tr_set%traj_data(iTmp,iStep,iTraj) = real_missing
        endif
      enddo
      tr_set%nbr_of_steps(iTraj) = iStep
    end do  ! iTraj

    ! Dump to output files  whatever is ready by now
    do iTmp = 1, tr_set%nReady
        iTraj = tr_set%ReadyToDumpList(iTmp)
        call dump_trajectory(tr_set, iTraj)
        ! Release the trajectory
        tr_set%nFree = tr_set%nFree + 1
        tr_set%FreeList(tr_set%nFree) = iTraj
        tr_set%iParticle(iTraj) = int_missing
    enddo
    tr_set%nReady = 0

    ! Adjust the highest valid  trajectory index
    do iTmp = tr_set%iTrajMax,1,-1
        if (tr_set%iParticle(iTmp) > 0) exit
    enddo
    tr_set%iTrajMax = iTmp

  END SUBROUTINE store_traj_time
 
 subroutine dump_trajectory(tr_set, iTraj)

  !Just dumps given trajectory to the corresponding file 
  ! Does not care about book keeping


  ! More or less follows this format:
  !
  !#/tmp/hirdata//trajectories_sobo1_1999_02_18_04.27UTC.silja
  !# produced by silja_trajectory on:cray using weather data:Atlantic HIRLAM FMIml
  !# case: sobo1
  !# comment: Sosnovyj Bor 29.-31.7.1998
  !# comment: 3 partikkelia/min*12h
  !# comment: Pilvi
  !# add. quantity:   1: boundary layer height from gound [m]
  !# add. quantity:   2: scavenging coefficient [1/s]
  !# add. quantity:   3: pasquill stability class
  ! 2160 trajectories
  !************  Trajectory number     1  ********************
  !  97 steps: time(6)UTC       lat     lon     pre;HPa  h;m  additional
  !    1 1998  7 29  6  0  0.0  59.900  29.080   994.0    32  357.327   0.328E-05    4.000 
  !    2 1998  7 29  6 30  0.0  59.957  28.823   964.9   285  319.136   0.102E-05    4.000 
  !    3 1998  7 29  7  0  0.0  60.064  28.489   986.0   124  336.109   0.102E-05    4.000 
  !
    IMPLICIT NONE
   type(silam_trajectory_set), intent(inout) :: tr_set
   integer, intent(in) :: iTraj
  
    INTEGER :: unit, iTr, iStep, yr,mon,day,hr,min, iSrc
    real :: x, y, sec

    iSrc = tr_set%iSrc(iTraj)
    unit = tr_set%ounit(iSrc)
    iTr = tr_set%traj_dumped(iSrc)+1  ! Counter of trajectories in a file
    
   

      WRITE(unit,'(A,I5,A)')'************  Trajectory number', iTr, ' ********************'
      WRITE(unit,'(I5,A)') tr_set%nbr_of_steps(iTraj), &
                         & ' steps: step, time(6)UTC, lat,    lon,  pre(HPa), height(m), extra'

      DO iStep = 1, tr_set%nbr_of_steps(iTraj)
        x = tr_set%traj_data(i_x, iStep, iTraj)
        y = tr_set%traj_data(i_y, iStep, iTraj)
        CALL get_time_components_utc(tr_set%time_start(iTraj) + tr_set%time_step*(iStep-1), &
                                 & yr,mon,day,hr,min,sec)
        WRITE(unit,'(I4,I5,4I3,f5.1,2(f9.3),2f8.0,f9.3,1x,e10.3,I3)') &
            & iStep, &
            & yr, mon, day, hr, min, sec, &
            & fu_lat_geographical_from_grid(x, y, meteo_grid), &
            & fu_lon_geographical_from_grid(x, y, meteo_grid), &
            & 0.01 * tr_set%traj_data(i_pressure, iStep, iTraj), &
            & tr_set%traj_data(i_height, iStep, iTraj), &
            & tr_set%traj_data(i_abl_h, iStep, iTraj), &
            & tr_set%traj_data(i_scav, iStep, iTraj), &
            & nint(tr_set%traj_data(i_stability, iStep, iTraj))
      END DO
      tr_set%traj_dumped(iSrc) = iTr
 END SUBROUTINE dump_trajectory


 subroutine finalize_trajectory_output(tr_set)
    ! Flushes rest of the trajrctories and resets the main tr_set
    IMPLICIT NONE
    type(silam_trajectory_set), intent(inout) :: tr_set
    
    integer :: iTraj, iSrc, iStat, unit
    character (len=*), parameter :: sub_name="finalize_trajectory_output"
    character (len=25) :: tmpstr

    if (.not. defined(tr_set)) then
         call set_error("tr_set undefined", sub_name)
         return
    endif
    call msg("finalize_trajectory_output: tr_set%ounit()", tr_set%ounit(:))
    call msg("finalize_trajectory_output: NoTraj_offset()", tr_set%NoTraj_offset(:))

    do iTraj = 1, tr_set%iTrajMax
       if (tr_set%iParticle(iTraj) == int_missing) cycle
       if (tr_set%nbr_of_steps(iTraj) > 2) call dump_trajectory(tr_set, iTraj)
       tr_set%iParticle(iTraj) = int_missing
       tr_set%nFree = tr_set%nFree + 1
       tr_set%FreeList(tr_set%nFree) = iTraj 
    enddo
    tr_set%iTrajMax = 0

    ! Put trajectories count and close the files
    do iSrc = 1,tr_set%nSrcs
      unit = tr_set%ounit(iSrc)
      tr_set%ounit(iSrc) = int_missing
      close(unit) !!Intel compiler requires unformatted access to write in the middle of the file
      !WARNING! Should be in sync with placeholder line in init_trajectory_output
      WRITE(tmpstr,'(I10,A15)') tr_set%traj_dumped(iSrc), ' trajectories'
      open(unit, file = tr_set%fname_tmp(iSrc), status="old", access="stream", iostat = istat)
      CALL FSEEK(unit, tr_set%NoTraj_offset(iSrc),  0, iStat) !status 
      write (unit) tmpstr
      close(unit)
      if ( RENAME(tr_set%fname_tmp(iSrc), tr_set%fname(iSrc)) /=0) then
         call set_error("Failed to rename "//trim(tr_set%fname_tmp(iSrc)//&
                     & "to" // trim(tr_set%fname(iSrc))),sub_name)
      endif
    enddo
    call msg("Done with trajectories output for "//trim(fu_str(tr_set%nSrcs))//" sources.")
    call msg(" Trajectries dumped for each source", tr_set%traj_dumped(1:tr_set%nSrcs))

    tr_set = trajectory_set_missing !Clean things up

 end SUBROUTINE finalize_trajectory_output

  LOGICAL FUNCTION fu_trajectory_set_defined(tr_set)
    !
    IMPLICIT NONE
    type(silam_trajectory_set), intent(in) :: tr_set

    if (tr_set%defined == silja_true) then
       fu_trajectory_set_defined = .true.
    else
       fu_trajectory_set_defined = .false.
       if (.not. tr_set%defined == silja_false) then
          call set_error("OOps!", "fu_trajectory_set_defined")
       endif
    endif
  END FUNCTION fu_trajectory_set_defined

END MODULE trajectory_io

