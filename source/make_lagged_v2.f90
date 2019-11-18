program make_lagged_v2
  ! 
  ! This program generates sets of grib files with a random time lag. The settings are
  ! given in a configuration file. The lag is fixed for each set and assumed gaussian.
  ! 
  ! For each valid time obtained from the config file, the corresponding met file is
  ! found. A file in the output directory is created with the same name, fields and metadata. 
  ! The data is replaced by interpolating (reading new input files as needed) into valid_time + lag.
  !
  ! The output appears in subdirectories numbered by the sample.
  !
  ! This program is really inefficient. Too bad.
  !
  use grib_api
  use silam_namelist
  use silam_times
  use grads_templates
  use input_data_rules, only : FNm_from_single_template, fu_closest_obstime
  use toolbox, only : random_normal
  implicit none

  type t_setup
     character(len=fnlen) :: output_dir
     type(grads_template) :: input_templ
     type(silja_interval) :: lag_stdev
     integer :: nbr_samples
     type(silja_time) :: first_time, last_time
     type(silja_interval) :: timestep
  end type t_setup

  integer :: stat, ind_sample, ind_nl
  character(len=fnlen) :: filename_conf, arg
  type(t_setup) :: setup
  type(Tsilam_namelist_group), pointer :: nlgrp
  type(Tsilam_namelist), pointer :: nlptr
  integer :: my_task, num_tasks
  logical :: is_shared

  call get_command_argument(1, arg)
  call get_workshare(arg, my_task, num_tasks, is_shared)
  if (error) return
  if (is_shared) then
    call get_command_argument(2, filename_conf)
  else
    call get_command_argument(1, filename_conf)
  end if

  if (filename_conf == '') then
    print *, 'Usage: make_lagged config_file_name'
    stop
  end if

  run_log_funit = fu_next_free_unit()
  open(run_log_funit, file=trim(filename_conf)//'.log')
  
  call get_nlgrp(filename_conf, nlgrp)
  if (error) return

  do ind_nl = 1, fu_nbr_of_namelists(nlgrp)
    nlptr => fu_namelist(nlgrp, ind_nl)
    if (error) return
    call msg('Reporing namelist:')
    call report(nlptr)

    call init(nlptr, setup)
    if (error) return
    call destroy_namelist(nlptr)

    call generate(setup)

  end do

  !************************************************************************************
  
contains

  subroutine get_workshare(arg, my_task, num_tasks, is_shared)
    implicit none
    character(len=*), intent(in) :: arg
    integer, intent(out) :: my_task, num_tasks
    logical, intent(out) :: is_shared

    character(len=8) :: wrk
    integer :: iostat, num_tokens
    integer, dimension(3) :: tokens
    
    wrk = trim(arg)
    ! ' -w1:1'
    if (wrk(1:2) == '-w') then
      call split_string(trim(wrk(3:)), ':', tokens, num_tokens)
      if (error) return
      if (fu_fails(num_tokens == 2, 'something wrong with splitting string', 'get_workshare')) return
      my_task = tokens(1)
      num_tasks = tokens(2)
      is_shared = .true.
    else
      my_task = 1
      num_tasks = 1
      is_shared = .false.
    end if

  end subroutine get_workshare

  subroutine get_nlgrp(filename_conf, nlgrp)
    implicit none
    character(len=*), intent(in) :: filename_conf
    type(Tsilam_namelist_group), pointer :: nlgrp

    integer :: file_unit, stat
    
    file_unit = fu_next_free_unit()
    open(file_unit, file=filename_conf, action='read', status='old', iostat=stat)
    if (fu_fails(stat == 0, 'Failed to open:'//trim(filename_conf), 'get_nlgrp')) return
    nlgrp => fu_read_namelist_group(file_unit, ifenv=.true.)
    if (error) return
    if (fu_fails(associated(nlgrp), 'No namelists read from '//trim(filename_conf), 'get_nlgrp')) return
    close(file_unit)

  end subroutine get_nlgrp
  
  subroutine init(nlptr, setup)
    implicit none
    type(Tsilam_namelist), pointer :: nlptr
    type(t_setup), intent(out) :: setup

    character(*), parameter :: sub_name = 'init'
    character(len=fnlen) :: content
    integer :: content_int, file_unit, stat
    
    if (fu_fails(associated(nlptr), 'No namelist found', sub_name)) return

    call chk_from_nl(nlptr, 'output_dir', setup%output_dir, sub_name)
    if (error) return
    call chk_from_nl(nlptr, 'input_templ', content, sub_name)
    if (error) return
    call decode_template_string(content, setup%input_templ)
    if (error) return
    call chk_from_nl(nlptr, 'lag_stdev', content, sub_name)
    if (error) return
    setup%lag_stdev = fu_set_named_interval(content)
    if (error) return
    content_int = fu_content_int(nlptr, 'nbr_samples')
    if (error .or. fu_fails(content_int /= int_missing, 'Missing nbr_samples', sub_name)) return
    setup%nbr_samples = content_int

    call chk_from_nl(nlptr, 'first_time', content, sub_name)
    if (error) return
    setup%first_time = fu_io_string_to_time(content)
    if (error) return
    call chk_from_nl(nlptr, 'last_time', content, sub_name)
    if (error) return
    setup%last_time = fu_io_string_to_time(content)
    if (error) return
    call chk_from_nl(nlptr, 'timestep', content, sub_name)
    if (error) return
    setup%timestep = fu_set_named_interval(content)
    
  end subroutine init
  
  subroutine chk_from_nl(nlptr, key, store, caller, default)
    implicit none
    type(Tsilam_namelist), pointer :: nlptr
    character(*), intent(in) :: key, caller
    character(*), intent(out) :: store
    character(*), intent(in), optional :: default
    
    character(len=fnlen) :: content

    content = fu_content(nlptr, key)
    if (content == '') then
      if (present(default)) then
        store = default
      else
        call set_error('Missing ' // trim(key), caller)
        return
      end if
    else
      store = content
    end if
        
  end subroutine chk_from_nl

  subroutine get_lags(stdev, lags)
    implicit none
    type(silja_interval), intent(in) :: stdev
    type(silja_interval), dimension(:), intent(out) :: lags

    real, dimension(size(lags)) :: sample_arr
    integer :: ii

    call random_normal(sample_arr)
    sample_arr = sample_arr * fu_sec(stdev)
    
    if (.not. fu_sort_real_array(sample_arr, ascending, ifHolesAllowed=.false.)) then
      call set_error('Sorting returned false', 'get_lags')
      return
    end if
    
    call msg('Reporting lags')
    do ii = 1, size(lags)
      call msg('seconds:', sample_arr(ii))
      lags(ii) = fu_set_interval_sec(sample_arr(ii))
    end do
    
  end subroutine get_lags
  
  subroutine generate(setup)
    implicit none
    type(t_setup), intent(in) :: setup

    type(silja_interval), dimension(setup%nbr_samples) :: lags

    integer, parameter :: max_grib_ids = 800
    character(len=*), parameter :: sub_name = 'generate'

    integer, dimension(:), allocatable :: grib_ids_out, levels_req, leveltypes_req, params_req
    real, dimension(:,:), allocatable :: fields_next, fields_prev, fields_out
    integer :: mesg_size, mesg_count, mesg_size_new
    character(len=fnlen) :: filename_now, filename_out
    integer :: grib_id_now, gribf_now, stat, gribf_out, ind_mesg, grib_ed
    type(silja_time) :: now, now_lag, time_next, time_prev, time_next_new, time_prev_new
    real :: weight_next, weight_prev
    logical :: alloc_done

    if (fu_fails(setup%output_dir /= '', 'no output directory', sub_name)) return
    do ind_sample = 1, setup%nbr_samples
      ! use sample-1 in naming so it can be accessed with %task
      call create_directory_tree(trim(setup%output_dir) // dir_slash // fu_str(ind_sample-1))
    end do
    if (error) return

    ! generate the sorted lags
    call init_random_seed()
    call get_lags(setup%lag_stdev, lags)
    if (error) return
    
    allocate(grib_ids_out(max_grib_ids), levels_req(max_grib_ids), leveltypes_req(max_grib_ids), &
           & params_req(max_grib_ids))
    alloc_done = .false.

    ! loop over valid times
    now = setup%first_time
    do while(now <= setup%last_time)
      call msg('Now:' // trim(fu_str(now)))
      ! ++ open, clone all messages from the valid time
      call get_input_filename(setup%input_templ, now, filename_now)
      if (error) return
      call grib_open_file(gribf_now, filename_now, 'r')
      mesg_count = 0
      mesg_size = -1
      do
        call grib_new_from_file(gribf_now, grib_id_now, stat)
        if (stat == GRIB_END_OF_FILE) then
          exit
        else if (stat /= GRIB_SUCCESS) then
          call set_error('Failed to read valid-time file', sub_name)
        end if
        mesg_count = mesg_count + 1
        if (fu_fails(mesg_count <= max_grib_ids, 'Too many messages', sub_name)) return
        call grib_get(grib_id_now, 'editionNumber', grib_ed)
        if (grib_ed == 1) then
          call grib_get(grib_id_now, 'indicatorOfParameter', params_req(mesg_count))
          call grib_get(grib_id_now, 'indicatorOfTypeOfLevel', leveltypes_req(mesg_count))
        else
          call grib_get(grib_id_now, 'parameterNumber', params_req(mesg_count))
          call grib_get(grib_id_now, 'typeOfFirstFixedSurface', leveltypes_req(mesg_count))
        end if
        call grib_get(grib_id_now, 'level', levels_req(mesg_count))
          
        call grib_get_size(grib_id_now, 'values', mesg_size_new)
        if (mesg_size > 0 .and. mesg_size_new /= mesg_size) then
          call set_error('Message size changed', sub_name)
          return
        else
          mesg_size = mesg_size_new
        end if
        call grib_clone(grib_id_now, grib_ids_out(mesg_count))
        call grib_release(grib_id_now)
      end do

      call msg('Number of messages:', mesg_count)
      call msg('Size of messages:', mesg_size)

      if (.not. alloc_done) then
        allocate(fields_next(mesg_size, mesg_count), &
               & fields_prev(mesg_size, mesg_count), fields_out(mesg_size, mesg_count))
        alloc_done = .true.
      end if

      ! ++ loop over lags
      do ind_sample = 1, setup%nbr_samples
        call msg('Sample:', ind_sample)
        now_lag = now + lags(ind_sample)
        ! ++ ++ get past, future times
        time_prev_new = fu_closest_obstime(now_lag, backwards, setup%timestep)
        ! ++ ++ if past == prev future, move data
        if (time_prev_new == time_prev) then
          continue
        else if (time_prev_new == time_next) then 
          fields_prev = fields_next
          time_prev = time_prev_new
        else
          ! ++ ++ else read data
          time_prev = time_prev_new
          call read_fields(time_prev, params_req, levels_req, leveltypes_req, fields_prev)
          if (error) return
        end if
        time_next_new = fu_closest_obstime(now_lag, forwards, setup%timestep)
        if (time_next_new == time_next) then
          continue
        else if (time_next_new == time_prev) then
          fields_next = fields_prev
          time_next = time_next_new
        else
          time_next = time_next_new
          call read_fields(time_next, params_req, levels_req, leveltypes_req, fields_next)
          if (error) return
        end if
        call msg('Lagged time:' // trim(fu_str(now_lag)))
        call msg('Previous met time:' // trim(fu_str(time_prev)))
        call msg('Next met time:' // trim(fu_str(time_next)))
        if (.not. (time_next == time_prev)) then
          ! for zero-lag times coincide
          if (fu_fails(time_next - time_prev == setup%timestep, 'steps don''t match', sub_name)) return
          weight_next = fu_sec(now_lag-time_prev) / fu_sec(time_next - time_prev)
        else
          weight_next = 1.0
        end if
      
        if (fu_fails(weight_next >= 0.0, 'Bad weight_prev', sub_name)) return
        weight_prev = 1 - weight_next
      
        call msg('Weight previous:', weight_prev)
        call msg('Weight next:', weight_next)

        ! ++ ++ generate new data with weights
        
        fields_out = fields_prev*weight_prev + fields_next*weight_next
        ! ++ ++ set message data, open, write, close
        filename_out = trim(setup%output_dir) // dir_slash // trim(fu_str(ind_sample-1)) &
             & // dir_slash // trim(fu_basename(filename_now)) // '.shifted'
        call grib_open_file(gribf_out, filename_out, 'w')
        do ind_mesg = 1, mesg_count
          call grib_set(grib_ids_out(ind_mesg), 'values', fields_out(1:mesg_size, ind_mesg))
          call grib_write(grib_ids_out(ind_mesg), gribf_out)
        end do
        call grib_close_file(gribf_out)
      end do ! sample
      do ind_mesg = 1, mesg_count
        call grib_release(grib_ids_out(ind_mesg))
      end do
    ! ++ ++ else read data for past, future, check match with valid
    ! ++ close valid-time input
      now = now + setup%timestep
    end do

  end subroutine generate
    
  subroutine read_fields(valid_time, params_req, levels_req, leveltypes_req, data_out)
    implicit none
    type(silja_time) :: valid_time
    integer, dimension(:), intent(in) :: params_req, levels_req, leveltypes_req
    real, dimension(:,:), intent(out) :: data_out

    integer :: nbr_mesg, size_mesg, grib_ed
    character(len=fnlen) :: filename
    integer :: gribf, gribid, curr_size, param, ind_mesg

    nbr_mesg = size(data_out, 2)
    size_mesg = size(data_out, 1)
    call get_input_filename(setup%input_templ, valid_time, filename)
    if (error) return
    call msg('Opening grib file: ' // trim(filename))
    call grib_open_file(gribf, filename, 'r')

    do ind_mesg = 1, nbr_mesg
      call grib_new_from_file(gribf, gribid, stat)
      ! if error, either the file has different content than expected, or a real error
      if (fu_fails(stat == GRIB_SUCCESS, 'Failed to read everything from input', 'read_fields')) return
      call grib_get(gribid, 'editionNumber', grib_ed)
      if (grib_ed == 1) then
        call grib_get(gribid, 'indicatorOfParameter', param)
      else
        call grib_get(gribid, 'parameterNumber', param)
      end if
      if (fu_fails(param == params_req(ind_mesg), 'Surprising parameter ' // trim(fu_str(ind_mesg)), 'read_fields')) return
      call grib_get(gribid, 'level', param)
      if (fu_fails(param == levels_req(ind_mesg), 'Surprising level', 'read_fields')) return
      if (grib_ed == 1) then
        call grib_get(gribid, 'indicatorOfTypeOfLevel', param)
      else
        call grib_get(gribid, 'typeOfFirstFixedSurface', param)
      end if
      if (fu_fails(param == leveltypes_req(ind_mesg), 'Surprising level type', 'read_fields')) return
      call grib_get_size(gribid, 'values', curr_size)
      if (fu_fails(curr_size == size_mesg, 'Surprising size of message', 'read_fields')) return
      ! finally clear with checking
      call grib_get(gribid, 'values', data_out(:,ind_mesg))
      call grib_release(gribid)
    end do
    call grib_close_file(gribf)

  end subroutine read_fields

  subroutine get_input_filename(templ, valid_time, filename)
    implicit none
    type(grads_template), intent(in) :: templ
    type(silja_time), intent(in) :: valid_time
    character(len=*), intent(out) :: filename
    
    character(len=*), parameter :: sub_name = 'get_input_filename'
    type(silam_sp), dimension(:), pointer :: filenames
    
    nullify(filenames)
    call fnm_from_single_template(templ, valid_time, filenames, &
                                & ifAdd = .false., &
                                & ifStrict = .false., &
                                & ifAllowZeroFcLen = .false., &
                                & ifWait = .false.)
    if (error) return
    if (fu_fails(associated(filenames), 'filenames not associated', sub_name)) return
    if (fu_fails(size(filenames) > 0, 'No filenames', sub_name)) return
    filename = filenames(1)%sp
    
    if (associated(filenames)) deallocate(filenames)
    
  end subroutine get_input_filename

  
end program make_lagged_v2
