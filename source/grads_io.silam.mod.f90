MODULE grads_io

  ! This module contains all necessary information and routine for 
  ! reading/writing of GrADS and GRIB files.
  ! Two types of binaries are considered: GrADS and GRIB. GrADS binaries
  ! are written here, while the GRIB ones are in the grib_io module.
  ! However, all information needed for GrADS to display both types
  ! of binaries is here.
  !
  ! Limitations for GrADS output are: 
  ! 1. Exact order of cycling: grid_indices -> levels -> variables -> time
  ! 2. Time step must be fixed in the output but not input, where e.g. 1 month is allowed
  ! 3. The list of variables and levels must exactly coinside between times
  ! 4. Horizontal grids must be exactly the same for all fields
  !
  ! As a result, during the first-time writing ALL information but time step 
  ! is collected and then must not be modified. The second writing completely 
  ! defines the file, except for the number of time steps.
  !
  ! Since there may be several GrADS + GRIB files written by a program, all
  ! vars are arrays
  !
  ! GRIB limitations assumed here: 
  ! - Time dimension is regular and motion along it is regular too at least 
  !   for the first and the second time periods. The rest can be in any order, 
  !   but after the first time period. 
  ! - Grid is unique or does not matter. It is taken for ctl from the first
  !   field and then is not checked at all.
  ! The rest is absolutely arbitrary.
  !
  ! Author: Mikhail Sofiev, FMI mikhail.sofiev@fmi.fi
  !
  ! All units: SI
  ! 
  ! Language: ANSI standard Fortran 90
  ! 
  !
  use stacks !ascii_io
  use silam_partitioning
  use optimisation, only : start_count, stop_count
  IMPLICIT NONE

  ! The public functions and subroutines available in this module:

  public init_grads_io
  public fu_next_free_grads_structure ! Finds a free structure
  PUBLIC fu_open_gradsfile_i      !-- file with GrADS input - via super-ctl file
  PUBLIC open_gradsfile_o         !-- file for GrADS output + ctl file
  public switch_grads_binary_o
  public read_field_from_grads_id    ! reads the requested field from GrADS file
  public read_field_from_grads_indices    ! reads the requested field from GrADS file
  PUBLIC write_next_field_to_gradsfile   ! writing command + ctl list enlarging
  PUBLIC close_gradsfile_i
  PUBLIC close_gradsfile_o          ! with writing the ctl file
  public invert_grads_binary      ! For negative time steps in adjoint runs
  PUBLIC fill_structures_for_grib ! fill structures for ctl file for GRIB binary
  PUBLIC init_ctl_for_grib   ! Called when GRIB file is opened
  PUBLIC write_ctl_for_grib  ! Same as for GrADS below, but more sophisticated

  PUBLIC fu_unit_bin       !-- returns the file unit pointed by unit_bin
  PUBLIC release_index     !-- Frees-up the structure
  public set_grib_binary_unit  ! Sets new value for the binary GRIB unit
  public fu_grads_time_index
  public fu_n_gvars
  public fu_n_gVar_levs
  public fu_n_glevs
  public fu_n_gtimes
  public fu_silamVert_of_grads
  public fu_silamGrid_of_grads
  public fu_time_of_grads
  public get_grads_times
  public get_grads_total
  public fu_validity_length
  public fu_data_time_features
  public get_grads_var_metadata
  public get_grads_IDs
  public report
  public fu_grads_filename
  public fu_grads_sctl_filename
  public fu_grads_missing_value
  public flush_buffer

  ! Private functions of this module
  private fu_glevel_type
  PRIVATE write_ctl_file  ! Uses the structures filled-in during writing
  PRIVATE fill_structures ! Fills the gfile structure during the first time period
  PRIVATE fu_get_glevel   ! Translates silja_level to glevel
  private wr_nxt_field_to_gf_from_fld
  private wr_n_fld_2_gf_from_id_and_data
  private report_grads_file
  private fu_find_grads_structure_i
  private fu_validity_length_from_grads
  private fu_data_time_features_grads
  private field_to_buffer
  private init_grads_buffer
  private free_grads_buffer
  private fu_get_default_buf_size
  private grads_str_2_times

  interface fu_validity_length
    module procedure fu_validity_length_from_grads
  end interface

  interface fu_data_time_features
    module procedure fu_data_time_features_grads
  end interface


  interface write_next_field_to_gradsfile
    module procedure wr_nxt_field_to_gf_from_fld
    module procedure wr_n_fld_2_gf_from_id_and_data
  end interface

  interface report
    module procedure report_grads_file
  end interface

  !-----------------------------------------------------------------------
  !
  !  GrADS-defining structures. 
  !
  TYPE grads_grid
    PRIVATE
    REAL :: x_start=-1., y_start=-1., x_step=-1., y_step=-1.
    INTEGER :: nx=-1,ny=-1
    TYPE(silja_logical) :: defined = silja_false
  END TYPE grads_grid
  type (grads_grid), private, parameter :: grads_grid_missing = &
                  &   grads_grid(-1.,-1.,-1.,-1.,-1,-1,silja_false)

  TYPE grads_levels
    PRIVATE
    REAL, DIMENSION(max_levels) :: levels=-1.
    TYPE(silja_logical) :: defined = silja_false
  END TYPE grads_levels
  type (grads_levels), private, parameter :: grads_levels_missing = &
          & grads_levels(-1.,silja_false)

  TYPE grads_variable
    PRIVATE
    INTEGER :: quantity = -1, grib_quantity = -1, n_levs = -1, grib_lev_typ = -1, &
             & grib_time_ind = -999
    character(len=substNmLen) :: chCocktailNm = ''  ! chSubstNm = '', 
    character(len=15) :: chGradsNm = ''
    type(silam_species) :: species
    type(silja_level) :: SilamLevel  ! may or may not be useful: 3D var has no definite level
    REAL :: grib_lev_val  = -1.0
    real :: fScalingFactor = 1.0
    integer :: iVerticalFeature = int_missing  ! do_nothing_flag, integrate_column_flag, lowest_level_flag, level_3D_type_flag, level_2d_type_flag
    TYPE(silja_logical) :: defined = silja_false
  END TYPE grads_variable
    type (grads_variable), private, parameter :: grads_variable_missing = & 
        & grads_variable(-1,-1,-1, -1,-999,"","", &
        & species_missing, level_missing, -1., 1., int_missing,   silja_false)

  TYPE grads_time
    PRIVATE
    TYPE(silja_time) :: start = time_missing, start_bin = time_missing
    TYPE(silja_interval) :: step = interval_missing
    type(silja_interval) :: validity_length = interval_missing
    logical :: ifVaryingStep = .false.   ! e.g. 1 month
    type(silja_time), dimension(:), pointer :: arTimes => null()  ! all times in the file
    TYPE(silja_logical) :: defined = silja_false
  END TYPE grads_time
  type (grads_time), private, parameter :: grads_time_missing = &
                          & grads_time(time_missing, time_missing, interval_missing, interval_missing, &
                                     & .false., null(), silja_false)

  TYPE grads_buffer
     real(r4k), dimension(:,:), pointer :: buf !fieldsize, max_flds
     logical :: allocated = .false.
     integer :: num_flds = 0, fieldsize = int_missing, max_flds = int_missing
  END TYPE grads_buffer

  !---------------------------------------------------------------------
  ! Structure completely defining a single GrADS file  
  !
  TYPE grads_file
    PRIVATE
    CHARACTER (LEN=fnlen) :: super_ctl_fname='', fname='', fname_initial = '', chTemplate=''
    type(grads_template) :: grTemplate
    INTEGER :: unit_bin=-1 !, unit_ctl=-1
    TYPE(grads_grid) :: ggrid
    type(silja_grid) :: silamGrid
    TYPE(grads_levels) :: glevs
    type(silam_vertical) :: silamVertical
    TYPE(grads_variable), DIMENSION(max_variables) :: gvars
    TYPE(grads_time) :: gtime
    type(grads_buffer) :: gradsbuf  
    integer, DIMENSION(max_variables) :: first_field_offset_in_tstep  !!  Only for input, offsets in fields
    INTEGER :: n_times=1, n_vars=1, n_levs=1, n_times_bin=1 ! Numbers of aliases
    INTEGER :: time_nbr=1, var_nbr=1, lev_nbr=1, rec=1 ! Expected field indices
    INTEGER :: levType = int_missing       ! Type of 3D levels in the file
    integer :: time_label_position = int_missing
    integer :: data_time_features = dynamic_map
    real :: missing_value
    logical :: ifBigEndian, ifYInverse, ifZInverse
    logical :: ifMPIIO=.false., ifBuffered=.false.
    TYPE(silja_logical) :: defined = silja_false
 END TYPE grads_file

  type grads_file_ptr
    private
    type(grads_file), pointer :: ptr
  end type grads_file_ptr

  !--------------------------------------------------------------------
  !  Further structure will be filled and used during writing a particular
  !  GrADS file
  !
  INTEGER, public, PARAMETER :: max_nbr_of_grads_files = 92 ! Start from unit=20, should not exceed 99
  INTEGER, private :: nbr_grads_files = 92 ! Start from unit=20, should not exceed 99^M 

  TYPE(grads_file_ptr), DIMENSION(:), PRIVATE, POINTER, SAVE :: gfile   


CONTAINS


!***************************************************************************
!
!   GrADS related stuff, including private functions
!
!***************************************************************************

!***********************************************************************

  subroutine init_grads_io(nFiles)
    !
    ! Just allocates the GrADS file structure
    !
    implicit none

    ! Imported variables
    integer, intent(in), optional :: nFiles

    ! Local variables
    integer :: iStat, iTmp

    if(present(nFiles))then
      if(nFiles > max_nbr_of_grads_files)then
        call set_error('Too large number of requested GrADS files','init_grads_io')
        return
      else
        allocate(gfile(nFiles), stat=iStat)
        nbr_grads_files = nFiles
      endif
    else
      allocate(gfile(max_nbr_of_grads_files), stat=iStat)
    endif

    if(iStat /= 0)then
      call set_error('Failed to allocate GrADS file structures','init_grads_io')
      return
    endif

    do iTmp = 1, size(gfile)
      allocate(gfile(iTmp)%ptr, stat=iStat)
      if(iStat /= 0)then
        call set_error('Failed to allocate GrADS file pointer','init_grads_io')
        return
      endif
    end do

  end subroutine init_grads_io


  !*************************************************************************

  integer function fu_next_free_grads_structure() result(index)
    !
    ! Finds the first free GrADS structure and returns its index
    !
    implicit none

    ! Local variables
    integer :: i

    do i=1, nbr_grads_files-1
      if(gfile(i)%ptr%defined == silja_false)then
        if(gfile(i)%ptr%unit_bin == -1)then
          index = i
          return
        endif
      endif
    end do
    call set_error('No free structures left','fu_next_free_grads_structure')

  end function fu_next_free_grads_structure 


  !*************************************************************************

  integer function fu_find_grads_structure_i(chCtlFNm)
    !
    ! Scans the filled GrADS structures looking for the ctl file name.
    !
    implicit none

    ! Imported parameters
    character(len=*), intent(in) :: chCtlFNm

    ! Local variables
    integer :: i

    fu_find_grads_structure_i = int_missing
    if(.not. associated(gfile))return
!call msg('Looking for the file:' + chCtlFNm)
    do i=1, nbr_grads_files-1
      if(gfile(i)%ptr%defined == silja_true)then
!call msg('Occupied GrADS structure:' + gfile(i)%ptr%fname, i)
        if(trim(gfile(i)%ptr%fname) == trim(chCtlFNm))then
          fu_find_grads_structure_i = i
          return
        endif
      endif
    end do

  end function fu_find_grads_structure_i


  !*****************************************************************

  integer function fu_open_gradsfile_i(chSuperCtlFNm)
    !
    ! Opens the GrADS file for input. Note: since the GrADS ctl file does not provide
    ! all necessary information about the file, there has to be a super-ctl file
    ! where this information is stored. This function reads this super-ctl file
    ! and then opens and consumes the basic ctl.
    ! The result of this exercise is stored into the GrADS structure. The function returns
    ! the place of the file in the gfile_ptr structure
    !
    !  ! restart from same-grid dump can be done with MPI read
    !  In case of #ifdef SILAM_MPI and smpi_use_mpiio_grads == .True. : check that
    !    1. domain split is used 
    !    and 
    !    2. The file has only one timestep
    !    and
    !    3. Grid is whole_mpi_dispersion_grid
    !
    !    In this case:
    !       1. MPI-read the whole binary (our subgrid) to the buffer
    !       2. Report this subdomain dispersion grid as the file grid
    !       3. set ifBuffered to .True.
    !       4. On read report the buffer 
    !       
    !
    implicit none

    ! Imported parameters
    character(len=*), intent(in) :: chSuperCtlFNm

    ! Local variables
    integer :: iUnit, iFile, iStatus, iVar, iLev, nPoints
    type(Tsilam_namelist_group), pointer :: nlGrpSCtlPtr
    type(Tsilam_namelist), pointer :: nlPtr
    logical :: eof, if_corner_in_geo_coord, if_south_pole, ifTemplate
    type(silam_sp) :: sp, sp_u_case, spTmp
    real :: pole_x, pole_y, fTmp
    real, dimension(:), pointer :: fPtr
    type(grads_file), pointer :: gf
    type(grads_variable), pointer :: gvar
    CHARACTER(LEN = fnlen) :: chTmp
    integer :: iTmp, nx, ny, offx, offy, gnx, gny, nFlds
    logical :: EndianSpecified, ifOK
    type(silam_sp), dimension(:), save, pointer :: fnames

    !
    ! Open and read the super-ctl file, which is in a namelist format
    ! It should contain just three parts: the complete description of grid and 
    ! a complete description of vertical. The rest is well presented by the ctl itself
    !
    sp%sp => fu_work_string()
    iUnit = fu_next_free_unit()
    if(error)return
    open(iUnit, file=chSuperCtlFNm, status='old', action='read', iostat= iStatus)
    if(iStatus /= 0)then
      call set_error('Super-ctl file does not exist:' + chSuperCtlFNm, 'grads_file_open_i')
      return
    endif

    nlGrpSCtlPtr => fu_read_namelist_group(iUnit, .false.)
    if(error)return

    nlPtr => fu_namelist(nlGrpSCtlPtr,'general')
    if(error .or. .not. associated(nlPtr))then
      call set_error('Cannot find general namelist in the super-ctl','fu_open_gradsfile_i')
      return
    endif

    chTmp = trim(fu_content(nlPtr,'ctl_file_name'))
    if(len(chTmp) < 1)then
      call set_error('The super-ctl file does not have ctl_file_name item', 'grads_file_open_i')
      return
    endif
    !
    ! Extend hat in ctl file name with superctl path
    !
    call replace_namelist_item(nlPtr, 'ctl_file_name', 'ctl_file_name', &
                             & fu_extend_grads_hat(chTmp,chSuperCtlFNm) )
    close(iUnit)

    !
    ! Check if this ctl file has already been opened
    !
    iFile = fu_find_grads_structure_i(fu_content(nlPtr,'ctl_file_name'))
    if(iFile /= int_missing)then
      fu_open_gradsfile_i = iFile
      call destroy_namelist_group(nlGrpSCtlPtr)
      return
    endif
    !
    !----------------------------------------------------------------------------
    !
    ! Real work starts
    !
    ! 1. Get the super-ctl content and translate it to GrADS terms
    ! 2. Read the ctl file and compare its content with the stored one
    ! 3. Fill-in the rest of the gfile structure
    !
    iFile = fu_next_free_grads_structure()
    if(error)return

    gf => gfile(iFile)%ptr
    !
    ! Some default stuff
    !
    gf%ifYInverse = .false.
    gf%ifZInverse = .false.
    !
    ! Consume the basic super-ctl information
    !
    gf%super_ctl_fname = trim(chSuperCtlFNm)
    gf%fname = fu_content(nlPtr,'ctl_file_name')
    gf%silamGrid = fu_set_grid(nlPtr)
    call set_vertical(nlPtr, gf%silamVertical)
!    call arrange_levels_in_vertical(gf%silamVertical)
    if(error)return

    sp%sp = fu_content(nlPtr, 'time_label_position')

    if(index(trim(fu_str_l_case(sp%sp)), 'start_of_period') > 0)then
      gf%time_label_position = start_of_period
    elseif(index(trim(fu_str_l_case(sp%sp)), 'mid_of_period') > 0)then
      gf%time_label_position = mid_of_period
    elseif(index(trim(fu_str_l_case(sp%sp)), 'end_of_period') > 0)then
      gf%time_label_position = end_of_period
    elseif(index(trim(fu_str_l_case(sp%sp)), 'instant_fields') > 0)then
      gf%time_label_position = instant_fields
    else
      call set_error('strange time_label_position in superctl:' + sp%sp, 'fu_open_gradsfile_i')
      call msg(trim(fu_str_l_case(sp%sp)), index(trim(fu_str_l_case(sp%sp)), 'end_of_period'))
    endif

    if(len_trim(fu_content(nlPtr,'data_time_features')) > 0)then
      if(index(fu_str_l_case(fu_content(nlPtr,'data_time_features')),'dynamic_data') > 0)then
        gf%data_time_features = dynamic_map
      elseif(index(fu_str_l_case(fu_content(nlPtr,'data_time_features')),'monthly_data') > 0)then
        gf%data_time_features = monthly_climatology
      elseif(index(fu_str_l_case(fu_content(nlPtr,'data_time_features')),'static_data') > 0)then
        gf%data_time_features = static_climatology
      else
        call msg('data_time_features can be dynamic_data/monthly_data/static_data, not:' + &
               & fu_content(nlPtr,'data_time_features'))
        call set_error('Wrong data_time_features line','fu_open_gradsfile_i')
      endif
      !
      ! Backward compatibility causes trouble: contradicting lines
      !
      if(len_trim(fu_content(nlPtr,'if_permanent_data')) > 0)then
        call set_error('if_permanent_data is obsolete, use data_time_features','fu_open_gradsfile_i')
        call msg('data_time_features can be dynamic_data/monthly_data/static_data')
        return
      endif
    else
      !
      ! Backward compatibility - may be, obsolete way of introducing the permanent data?
      !
      if(fu_str_u_case(fu_content(nlPtr,'if_permanent_data')) == 'YES')then
        call msg('data_time_features can be dynamic_data/monthly_data/static_data')
        call set_error('if_permanent_data is obsolete, use data_time_features','fu_open_gradsfile_i')
        call unset_error('Compatibility patch applied by fu_open_gradsfile_i')
        gf%data_time_features = static_climatology
      else
        gf%data_time_features = dynamic_map
      endif
    endif

    if(len_trim(fu_content(nlPtr, 'validity_length')) > 0)then
      gf%gtime%validity_length = &
                                     & fu_set_named_interval(fu_content(nlPtr, 'validity_length'))
    else
      if(gf%data_time_features == static_climatology)then
        gf%gtime%validity_length = very_long_interval
      else
        gf%gtime%validity_length = zero_interval
      endif
    endif

    !
    ! Here we translate the normal grid and vertical to GrADS terms but only for checking with 
    ! the ctl file. Actual reading will use the normal SILAM grid and vertical
    !
    gf%n_levs = fu_NbrOfLevels(gf%silamVertical)
    do iLev = 1, gf%n_levs
      iStatus = fu_glevel_type(fu_leveltype(gf%silamVertical))
      gf%glevs%levels(iLev) = fu_get_glevel(fu_level(gf%silamVertical, iLev), &
                                                    & iStatus)
    end do
    if(error)return

    SELECT CASE(fu_gridtype(gf%silamGrid))
    CASE(lonlat)
      CALL lonlat_grid_parameters(gf%silamGrid, &
                         & gf%ggrid%x_start, gf%ggrid%y_start, &
                         & if_corner_in_geo_coord,&
                         & gf%ggrid%nx, gf%ggrid%ny, &
                         & pole_x, pole_y,  & 
                         & gf%ggrid%x_step, gf%ggrid%y_step)
    case (anygrid)
      call grid_dimensions(gf%silamGrid, gf%ggrid%nx, gf%ggrid%ny)
      gf%ggrid%y_start = fu_lat_geographical_from_grid(1., 1., gf%silamGrid)
      gf%ggrid%x_start = fu_lon_geographical_from_grid(1., 1., gf%silamGrid)
      gf%ggrid%x_step = fu_dx_cell_deg(gf%silamGrid, 1, 1)
      gf%ggrid%y_step = fu_dy_cell_deg(gf%silamGrid, 1, 1)

    CASE DEFAULT
      CALL set_error('Only latlon and any-grids so far','fu_open_gradsfile_i')
      RETURN
    END SELECT
    IF(error)RETURN

    !
    ! Step 2: read the ctl file and compare with the above parameters
    !
    iUnit = fu_next_free_unit()
    IF(error)RETURN
    
    open(iUnit, file=gf%fname, action='read', status='old', iostat=iStatus)
    if(iStatus /= 0)then
      call set_error('Failed to open the ctl file:' + gf%fname, 'fu_open_gradsfile_i')
      return
    endif

    sp_u_case%sp => fu_work_string()
    spTmp%sp => fu_work_string()
    fPtr => fu_work_array()
    eof = .false.
    ifTemplate = .false.

    EndianSpecified = .false. 
    !
    ! Start the grand cycle through the ctl file
    !
    do while(.not.eof)

      call next_line_from_input_file(iUnit, sp%sp, eof)
      if(error.or.eof)exit
      sp_u_case%sp = fu_str_u_case(sp%sp)
      iStatus = 0

      if(index(sp_u_case%sp,'DSET') == 1)then
        chTmp = adjustl(trim(sp%sp(6:)))    
        gf%chTemplate = fu_extend_grads_hat(chTmp,gf%fname)
!        if (index(chTmp, "^") == 1) then ! path relative to ctl location
!                iTmp = index(gf%fname,dir_slash,.True.) !Last slash
!                gf%chTemplate = &
!                   & fu_connect_strings(gf%fname(1:iTmp),chTmp(2:len(chTmp))) !fname
!        else
!                gf%chTemplate = chTmp 
!        endif

!        gf%chTemplate = adjustl(gf%chTemplate)
        call decode_template_string(gf%chTemplate, gf%grTemplate)

      elseif(index(sp_u_case%sp,'TITLE') == 1)then

      elseif(index(sp_u_case%sp,'OPTIONS') == 1)then

        ifOK = .false.
        if(index(sp_u_case%sp,'ENDIAN') > 0)then
          gf%ifBigEndian = (index(sp_u_case%sp,'BIG_') > 0)
          EndianSpecified = .true. 
          ifOK = .true.
        endif
        if(index(sp_u_case%sp,'TEMPLATE') > 1)then
          ifTemplate = .true.
          ifOK = .true.
        endif
        if(index(sp_u_case%sp,'YREV') > 0)then
          gf%ifYInverse = .true.
          ifOK = .true.
        endif
        if(index(sp_u_case%sp,'ZREV') > 0)then
          gf%ifZInverse = .true.
          ifOK = .true.
        endif
        if(.not. ifOK)then
          call set_error('Unknown option in ctl file:' + sp_u_case%sp, 'fu_open_gradsfile_i')
          return
        endif
        

      elseif(index(sp_u_case%sp,'UNDEF') == 1)then

        read(unit=sp%sp,fmt=*,iostat=iStatus) spTmp%sp, gf%missing_value

      elseif(index(sp_u_case%sp,'XDEF') == 1)then

        if(index(sp_u_case%sp,'LINEAR') > 0)then
          read(unit=sp%sp,fmt=*,iostat=iStatus) spTmp%sp, nPoints, spTmp%sp, fPtr(1), fPtr(2) ! to check
          if(nPoints /= gf%ggrid%nx .or. &
           & .not. (fPtr(1) .eps. gf%ggrid%x_start) .or. &
           & .not. (fPtr(2) .eps. gf%ggrid%x_step))then
            call msg('X-dimension of grid is incorrect:' + sp%sp)
            call report(gf%silamGrid)
            call set_error('X-dimension of grid is incorrect','fu_open_gradsfile_i')
            return
          endif
        else
          call set_error('Strange XDEF line: only linear x-axe allowed' + sp%sp, &
                       & 'fu_open_gradsfile_i')
         return
        endif

      elseif(index(sp_u_case%sp,'YDEF') == 1)then

        if(index(sp_u_case%sp,'LINEAR') > 0)then
          read(unit=sp%sp,fmt=*,iostat=iStatus) spTmp%sp, nPoints, spTmp%sp, fPtr(1), fPtr(2)  ! to check
          if(nPoints /= gf%ggrid%ny .or. &
           & .not. (fPtr(1) .eps. gf%ggrid%y_start) .or. &
           & .not. (fPtr(2) .eps. gf%ggrid%y_step))then
            call msg('Y-dimension of grid is incorrect:' + sp%sp)
            call report(gf%silamGrid)
            call set_error('Y-dimension of grid is incorrect','fu_open_gradsfile_i')
            return
          endif
        else
          call set_error('Strange XDEF line: only linear x-axe allowed' + sp%sp, &
                       & 'fu_open_gradsfile_i')
          return
        endif

      elseif(index(sp_u_case%sp,'ZDEF') == 1)then

        if(index(sp_u_case%sp,'LEVELS') > 0)then
          read(unit=sp%sp,fmt=*,iostat=iStatus) spTmp%sp, nPoints, spTmp%sp, &
                                              & (fPtr(iLev), iLev=1,gf%n_levs)
          if(nPoints /= gf%n_levs)then
            call msg('Z-dimension is incorrect:' + sp%sp + ', number of levels mismatch:', &
                   & nPoints, gf%n_levs)
            call report(gf%silamVertical)
            call set_error('Z-dimension size is incorrect (number of levels)','fu_open_gradsfile_i')
            exit
          endif
          if (gf%n_levs == 1)then       ! Single-level can be funny, allow it
            if(((gf%glevs%levels(1) .eps. 0.0) .and. fPtr(1) < 0.) .or. &
             & (gf%glevs%levels(1) .eps. fPtr(1))) then
              continue
            else
              call msg('Single-layer Z-dimension is incorrect:' + sp%sp)
              call report(gf%silamVertical, .true.)
              call set_error('Z-dimension size is incorrect','fu_open_gradsfile_i')
              exit
            endif
          else                                        ! multi-level structures must match
            do iLev = 1, gf%n_levs
              if(.not. (fPtr(iLev) .eps. gf%glevs%levels(iLev)))then
                call msg('Z-dimension is incorrect:' + sp%sp + ', Level number and values:' + &
                       & fu_str(iLev),fPtr(iLev), gf%glevs%levels(iLev))
                call msg('Super-ctl vertical:')
                call report(gf%silamVertical)
                call msg('Super-ctl vertical translated to grads:', gf%glevs%levels)
                call set_error('Z-dimension is incorrect','fu_open_gradsfile_i')
                return
              endif
            end do
          end if
        else
          call set_error('Strange ZDEF line: only levls z-axe allowed' + sp%sp, 'fu_open_gradsfile_i')
          return
        endif

      elseif(index(sp_u_case%sp,'TDEF') == 1)then

        call grads_str_2_times(sp%sp, &  ! the TDEF line
                             & gf%n_times, & ! time dimension
                             & gf%gtime%start, &  ! first GrADS time
                             & gf%gtime%step, & ! defined if not varying step
                             & gf%gtime%ifVaryingStep, &   ! e.g. 1 month?
                             & gf%gtime%arTimes)  ! ifVarying stores all times
        if(error)then
          call msg_warning('Failed to read TDEF time line:' + sp%sp)
          return
        endif

      elseif(index(sp_u_case%sp,'VARS') == 1)then

        read(unit=sp%sp,fmt=*,iostat=iStatus) spTmp%sp, gf%n_vars

        do iVar = 1, gf%n_vars
          call next_line_from_input_file(iUnit, sp%sp, eof)
          if(error .or.eof)then
            call set_error('Failed to read all ctl variables','fu_open_gradsfile_i')
            return
          endif
          gvar => gf%gvars(iVar)
          read(sp%sp,fmt=*,iostat=iStatus) spTmp%sp, &
                                         & gvar%n_levs, &
                                         & gvar%grib_quantity, &
                                         & gvar%grib_lev_typ
          if(iStatus /= 0)then
            call set_error('Failed to read ctl variable:' + sp%sp,'fu_open_gradsfile_i')
            return
          endif

          !
          ! In many cases the GrADS variable is shorter than needed for complete description of SILAM 
          ! output var. This is what the super-ctl is made for - it contains a complete description 
          ! of the whole thing.
          !
          nlPtr => fu_namelist(nlGrpSCtlPtr,spTmp%sp)
          if(error .or. .not. associated(nlPtr))then
            call set_error('Cannot find variables in the super-ctl:' + spTmp%sp,'fu_open_gradsfile_i')
            return
          endif
          !
          ! Quantity (other parameters are not present in the string):
          !
          call decode_id_params_from_io_str(fu_content(nlPtr,'quantity_short_name'), &
                                          & gvar%n_levs > 1, &
                                          & gvar%quantity, &
                                          & gvar%species, &
                                          & .true.)  ! ifPush
          !
          ! Vertical features. Note that for some levels info is in super-ctl
          !
          gvar%SilamLevel = level_missing

          if(fu_str_u_case(fu_content(nlPtr,'vertical_feature')) == 'FIELD_2D')then
            gvar%iVerticalFeature = level_2D_type_flag
            if(gvar%n_levs > 1)then
              call set_error('More than one level for 2D variable:' + &
                           & fu_quantity_string(gvar%quantity), &
                           & 'fu_open_gradsfile_i')
            endif
            if(fu_content(nlPtr,'vert_level') /= '') &
                 & call set_named_level_with_fract(fu_get_item(nlPtr,'vert_level'), &
                                                 & gvar%SilamLevel, fTmp)

          elseif(fu_str_u_case(fu_content(nlPtr,'vertical_feature')) == 'FIELD_3D')then
            gvar%iVerticalFeature = level_3D_type_flag
            if(gvar%n_levs < 1)then
              call set_error('Less than one level for 3D variable:' + &
                           & fu_quantity_string(gvar%quantity), &
                           & 'fu_open_gradsfile_i')
            endif

          elseif(fu_str_u_case(fu_content(nlPtr,'vertical_feature')) == 'COLUMN_INTEGRATED')then
            gvar%iVerticalFeature = integrate_column_flag
            if(gvar%n_levs > 1)then
              call set_error('More than one level for column-integrated variable:' + &
                           & fu_quantity_string(gvar%quantity), &
                           & 'fu_open_gradsfile_i')
            endif
            if(fu_content(nlPtr,'vert_level') /= '') &
                 & call set_named_level_with_fract(fu_get_item(nlPtr,'vert_level'), &
                                          & gvar%SilamLevel, fTmp)

          elseif(fu_str_u_case(fu_content(nlPtr,'vertical_feature')) == 'LOWEST_LEVEL')then
            gvar%iVerticalFeature = lowest_level_flag
            if(gvar%n_levs > 1)then
              call set_error('More than one level for lowest-level variable:' + &
                           & fu_quantity_string(gvar%quantity), &
                           & 'fu_open_gradsfile_i')
            endif
            if(fu_content(nlPtr,'vert_level') /= '') &
                 & call set_named_level_with_fract(fu_get_item(nlPtr,'vert_level'), &
                                          & gvar%SilamLevel, fTmp)

          else
            if(gvar%n_levs == 0)then  ! backward compatibility
              gvar%iVerticalFeature = level_2D_type_flag
            elseif(gvar%n_levs == gf%n_levs)then
              gvar%iVerticalFeature = level_3D_type_flag
            else
              call set_error('Strange number fo levels:' + sp%sp,'fu_open_gradsfile_i')
            endif
            if(fu_content(nlPtr,'vert_level') /= '') &
                 & call set_named_level_with_fract(fu_get_item(nlPtr,'vert_level'), &
                                          & gvar%SilamLevel, fTmp)
          endif  ! vertical features type
          if(error)return
          !
          ! Correction factor, if any. Either deduce it from unit or from explicit declaration
          !
          if(len_trim(fu_content(nlPtr,'unit')) > 0)then
            gvar%fScalingFactor = &
                                    & fu_set_named_value('1.0  ' // fu_content(nlPtr,'unit'))
            if(len_trim(fu_content(nlPtr,'scaling_factor')) > 0)then
              call set_error('Both scaling_factor and unit are given','fu_open_gradsfile_i')
              call report(nlPtr)
              return
            endif
          elseif(len_trim(fu_content(nlPtr,'scaling_factor')) > 0)then
            gvar%fScalingFactor = fu_content_real(nlPtr,'scaling_factor')
          else
            gvar%fScalingFactor = 1.0
          endif
          if(error .or. (gvar%fScalingFactor .eps. real_missing))then
            call set_error('cannot get the scaling factor from the namelist', 'fu_open_gradsfile_i')
            call report(nlPtr)
            return
          endif
          !
          ! And the species, if any
          !
          if(len_trim(fu_content(nlPtr,'substance_name')) > 0)then
            if(len_trim(fu_content(nlPtr,'cocktail_name')) > 0) then
              call set_error('substance_name: ' + fu_content(nlPtr,'substance_name') + &
                           & ', and cocktail_name: ' + fu_content(nlPtr,'cocktail_name') + &
                           & ' cannot be given simultaneously', &
                           & 'fu_open_grads_file_i')
              return
            endif
            call set_species(gvar%species, &
                           & fu_get_material_ptr(fu_content(nlPtr,'substance_name')), &
                           & fu_set_mode(nlPtr), &
                           & fu_set_named_value(fu_content(nlPtr,'optical_wavelength'),.true.)) ! keep silent
          else   ! if no substance_name, concktail_name can be given
            gvar%species = species_missing
            if(len_trim(fu_content(nlPtr,'cocktail_name')) > 0)then
              gvar%chCocktailNm = fu_content(nlPtr,'cocktail_name')
            else
              gvar%chCocktailNm = ''
            endif
          endif   ! is substance_name is given
          if(error)then
            call set_error('Failed the variable:' + sp%sp,'fu_open_gradsfile_i')
            return
          endif
        end do  ! cycle over variables

      elseif(index(sp_u_case%sp,'ENDVARS') == 1)then
        exit            ! all done
      else
        call msg_warning('Unknown ctl line:' + sp%sp,'fu_open_gradsfile_i')
      endif

      if(iStatus /= 0)then
        call set_error('Something strange in ctl line:' + sp%sp,'fu_open_gradsfile_i')
        return
      endif

    end do  ! reading the ctl file

    close(iUnit)

    gf%first_field_offset_in_tstep(1) = 0
    do iVar = 1, gf%n_vars
       gvar => gf%gvars(iVar)
       gf%first_field_offset_in_tstep(iVar+1) = gf%first_field_offset_in_tstep(iVar) + &
                            & max(1,gvar%n_levs) ! grads allows 0 levels for 2D vars but here we do not
    enddo
   !!    gf%first_field_offset_in_tstep(gf%n_vars+1) has total number of fields per time step



    
    ! little_endian by default: makes sense when there is no "OPTIONS" line in .ctl
    if (.not. EndianSpecified) then
            gf%ifBigEndian = .false. ! Force little endian
      call msg("Reading .ctl file:" + gf%fname)
            call msg_warning("No endian specified, forcing LITTLE_ENDIAN!","fu_open_gradsfile_i")
    endif

    call free_work_array(sp%sp)
    call free_work_array(sp_u_case%sp)
    call free_work_array(spTmp%sp)
    call free_work_array(fPtr)

    if(error)return

    fu_open_gradsfile_i = iFile
    gfile(iFile)%ptr%defined = silja_true

    gf%ifBuffered = .false.
#ifdef SILAM_MPI 
    if (gf%silamGrid == wholeMPIdispersion_grid .and. smpi_use_mpiio_grads)then
      if (.not. gf%ifYInverse  .and.  gf%n_times == 1) then
         call smpi_get_decomposition(nx, ny, offx, offy, gnx, gny)
         if (nx /= gnx .or. ny /= gny) then !! Makes sense to read fields to the buffer
           gf%ifBuffered = .true.
           if (gf%gradsbuf%allocated) call free_grads_buffer(gf%gradsbuf)
           !Just for the whole timestep
           nflds = gf%first_field_offset_in_tstep(gf%n_vars+1)
           call init_grads_buffer(gf%gradsbuf, nx*ny, nFlds)

           call FNm_from_single_template(gf%grTemplate, &
                                        & fu_time_of_grads(iFile, 1), &
                                        & fnames, &
                                        & anal_time = fu_time_of_grads(iFile, 1), &  ! analysis time
                                        & ifStrict = .false., &
                                        & ifAdd = .false., &  ! to the fnames array
                                        & ifAllowZeroFcLen = .true., &
                                        & ifWait = .false.)
            if(error)return


            !  call open_grads_binary_i(fnames(1)%sp, gfile(igf)%ptr%unit_bin, nPoints, &
            !                         & gfile(igf)%ptr%ifBigEndian)

            call smpi_read_grads_fieldset_parallel(gf%gradsbuf%buf, nx*ny, nFlds, fnames(1)%sp)
            !!! Pretend that our subdomain grid is Grads GRID
            gf%silamGrid = dispersion_grid
            CALL lonlat_grid_parameters(gf%silamGrid, &
                         & gf%ggrid%x_start, gf%ggrid%y_start, &
                         & if_corner_in_geo_coord,&
                         & gf%ggrid%nx, gf%ggrid%ny, &
                         & pole_x, pole_y,  & 
                         & gf%ggrid%x_step, gf%ggrid%y_step)
         endif
       endif
    endif
#endif
    call msg("Grid_from_grads")
    call report(gf%silamGrid)


!
!   DEBUG MISSION ONLY
!
!gfile(iFile)%ptr%fname = fu_connect_strings(gfile(iFile)%ptr%fname,'_debug')
!call write_ctl_file(gfile(iFile)%ptr)  ! It will write the ctl file
!gfile(iFile)%ptr%fname = gfile(iFile)%ptr%fname(1:len_trim(gfile(iFile)%ptr%fname)-6)
  end function fu_open_gradsfile_i

  !*********************************************************************************

  subroutine grads_str_2_times(chLine, &  ! the TDEF line
                             & n_times, & ! time dimension
                             & timeStart, &  ! first GrADS time
                             & intervStep, & ! defined if not varying step
                             & ifVaryingStep, &   ! e.g. 1 month?
                             & arTimes)  ! ifVarying stores all times
    !
    ! Converts the GrADS-type string to the SILAM type structure
    ! The input template example is: TDEF 2 LINEAR 00:00Z01jan1900  1dy
    ! For start time the template is: hh:mmZddmmmyyyy
    ! For example 25.4.1996 at 06.00 UTC for 06:00Z25apr1996 
    ! 
    implicit none

    ! Imported parameters
    character(len=*), intent(in) :: chLine  ! the TDEF line
    integer, intent(out) :: n_times          ! time dimension
    type(silja_time), intent(out) ::  timeStart     ! first GrADS time
    type(silja_interval), intent(out) :: intervStep ! defined if not varying step
    logical, intent(out) :: ifVaryingStep           ! e.g. 1 month?
    type(silja_time), dimension(:), pointer :: arTimes ! if varying step

    ! Local variables
    integer :: iTmp, year, month, day, hour, minute, iStatus
    real :: second
    character(len=3) :: chTmp
    character(len=fnlen) :: strTmp
    !
    ! Stupidity check: only LINEAR time series are allowed
    !
    if(index(fu_str_u_case(chLine),'LINEAR') == 0)then
      call set_error('Strange TDEF line: only linear t-axe allowed ' + chLine, &
                   & 'grads_str_2_times')
      return
    endif

    second = 0.
    minute = 0
    hour = 0
    n_times = -1

    !
    ! Get the number of time steps and prepare to get the start time 
    !
    read(unit=chLine,fmt=*,iostat=iStatus) strTmp, n_times, strTmp, strTmp
    if(iStatus /= 0)then
      call set_error('Failed to get the number of times from:' + chLine,'grads_str_2_times')
      return
    endif
    strTmp = fu_str_u_case(strTmp)

    !
    ! Determine the start time
    !
    iTmp = index(strTmp,'Z')
    if(iTmp > 0)then  ! time exists
      if(index(strTmp(1:iTmp),':') == 0)then ! no minutes
        read(unit=strTmp,fmt='(I2)', iostat=iStatus) hour
      else
        read(unit=strTmp,fmt='(I2,A1,I2)', iostat=iStatus) hour, chTmp, minute
      endif
    endif
    if(iStatus /= 0)then
      call set_error('Failed to get time from:' + strTmp, 'grads_str_2_times')
      return
    endif
    !
    ! Determine the start date
    !
    read(unit=strTmp(iTmp+1:), fmt='(I2,A3,I4)', iostat=iStatus) day, chTmp, year
    if(iStatus /= 0)then
      call set_error('Failed to get date from:' + strTmp, 'grads_str_2_times')
      return
    endif

    month=0
    do iTmp=1,12
      if(fu_str_u_case(chTmp) == fu_str_u_case(chMonthNames_3chr(iTmp)))then
        month=iTmp
        exit
      end if
    end do
    if(month == 0)then
      call set_error('Strange month name:' + chTmp, 'grads_str_2_times')
      return
    endif

    timeStart = fu_set_time_utc(year, month, day, hour, minute, second)

    !
    ! Get the interval and set either regular time step or fill-in the 
    ! whole bunch of GrADS times
    !
    read(unit=chLine,fmt=*,iostat=iStatus) strTmp, n_times, strTmp, strTmp, strTmp

    if(index(strTmp,'dy') > 0)then
      read(unit=strTmp(1:index(strTmp,'dy')-1), fmt=*, iostat=iStatus) iTmp
      intervStep = fu_set_interval_h(iTmp*24)
      nullify(arTimes)
      ifVaryingStep = .false.

    elseif(index(strTmp,'hr') > 0)then
      read(unit=strTmp(1:index(strTmp,'hr')-1), fmt=*, iostat=iStatus) iTmp
      intervStep = fu_set_interval_h(iTmp)
      nullify(arTimes)
      ifVaryingStep = .false.

    elseif(index(strTmp,'mn') > 0)then
      read(unit=strTmp(1:index(strTmp,'mn')-1), fmt=*, iostat=iStatus) iTmp
      intervStep = fu_set_interval_min(iTmp)
      nullify(arTimes)
      ifVaryingStep = .false.

    elseif(index(strTmp,'mo') > 0 .or. &
         & index(strTmp,'yr') > 0)then
      !
      ! Irregular interval! Set full array of times
      !
      intervStep = interval_missing
      ifVaryingStep = .true.
      if(n_times > 10000)then
        call msg('Very large number of irregular GrADS times:', n_times)
        call set_error('Very large number of irregular GrADS times:','grads_str_2_times')
        return
      endif
      allocate(arTimes(0:n_times+1), stat = iStatus)
      if(iStatus /= 0)then
        call set_error('Failed to allocate times array','grads_str_2_times')
        return
      endif

      if(index(strTmp,'mo') > 0)then
        read(unit=strTmp(1:index(strTmp,'mo')-1), fmt=*, iostat=iStatus) iTmp
      else
        read(unit=strTmp(1:index(strTmp,'yr')-1), fmt=*, iostat=iStatus) iTmp
      endif
      !
      ! Set one step backwards and the starting time
      !
      if(month == 1)then
        arTimes(0) = fu_set_time_utc(year-1, 12, day, hour, minute, second)
      else
        arTimes(0) = fu_set_time_utc(year, month-1, day, hour, minute, second)
      endif
      arTimes(1) = timeStart              ! first time is the start
      !
      ! And now all other times, plus one step forwards
      !
      do iStatus = 2, n_times+1              ! fill-in the rest
        if(index(strTmp,'mo') > 0)then
          month = month + iTmp
          if(month > 12)then
            year = year + 1
            month = 1
          endif
        else
          year = year + iTmp
        endif
        arTimes(iStatus) = fu_set_time_utc(year, month, day, hour, minute, second)
        if(error)return
      end do  ! n_times-1

    else
      call set_error('Cannot handle the string:' + strTmp, 'grads_str_2_times')
      return
    endif

  end subroutine grads_str_2_times


  
  !*********************************************************************************
  
  integer function open_grads_file_o_full(dir, fname, chTemplate, fMissingVal)
    !
    ! Opens the GrADS and ctl files, defines grads_grid (mandatory for binary 
    ! opening) and returns the number of the file in the GrADS file list gfile_ptr
    ! Uses binary opening from md module.
    ! Also, receives a complete list of field IDs to be stored in the file later on.
    ! This allows a full definition of the GrADS file structure, so that later on it only 
    ! finds out the correct record where the field is to be put and writes it down.
    !
    ! Author: M.Sofiev, FMI
    !
    IMPLICIT NONE
    !
    ! Imported paraneters with intent IN
    CHARACTER (LEN=*), INTENT(in) :: dir, fname ! Directory and name of the GrADS file
    CHARACTER (LEN=*), INTENT(in) :: chTemplate
    real, intent(in) :: fMissingVal

    call set_error('Not ready yet','open_grads_file_o_full')
    open_grads_file_o_full = int_missing
  end function open_grads_file_o_full

  !************************************************************************************

  INTEGER FUNCTION open_gradsfile_o(dir, fname, grid, chTemplate, fMissingVal, ifMPIIO, ifBuffered, &
                                  & time_label_position)
    !
    ! Opens the GrADS and ctl files, defines grads_grid (mandatory for binary 
    ! opening) and returns the number of the file in the GrADS file list gfile_ptr
    ! Uses binary opening from md module.
    !
    ! Author: M.Sofiev, FMI
    !
    IMPLICIT NONE
    !
    ! Imported paraneters with intent IN
    CHARACTER (LEN=*), INTENT(in) :: dir, fname ! Directory and name of the GrADS file
    TYPE(silja_grid), INTENT(in) :: grid ! will be stored to ggrid
    CHARACTER (LEN=*), INTENT(in), optional :: chTemplate
    real, intent(in), optional :: fMissingVal
    logical, intent(in), optional :: ifMPIIO, ifBuffered
    integer, intent(in), optional :: time_label_position
    ! Note, ifBuffered allows buffering but ultimately the choice dependes on
    ! GRADSIO_BUF_SIZE environment variable.

    ! Local declarations
    INTEGER :: NbrOfPoints ! in the grid
    INTEGER :: iFile,nx,ny
    LOGICAL :: if_corner_in_geo_coord, if_south_pole, if_parallel
    REAL :: pole_x, pole_y, x,y,xs,ys

    open_gradsfile_o = -1

    !---------------------------------------------------------
    !
    ! Find a spare file counter. The last structure is always kept free!!
    !
    iFile = fu_next_free_grads_structure()
    if(error)return

    if(present(ifMPIIO)) then
      if(ifMPIIO.and.(.not.smpi_is_mpi_version())) then
        call set_error('Trying to open file in parallel mode - this is &
            &a serial version', 'open_grads_file_o')
        return
      end if
      if_parallel = ifMPIIO
    else
      if_parallel = .false.
    end if
    
    if (present(ifBuffered)) then
      gfile(iFile)%ptr%ifBuffered = ifBuffered .and. fu_get_default_buf_size() > 0
    else
      gfile(iFile)%ptr%ifBuffered = .false.
    end if

    if(present(fMissingVal)) then
      gfile(iFile)%ptr%missing_value = fMissingVal
    else
      gfile(iFile)%ptr%missing_value = real_missing
    endif

    if(present(time_label_position))then
      gfile(iFile)%ptr%time_label_position = time_label_position
    else
      gfile(iFile)%ptr%time_label_position = end_of_period
    endif
    
    gfile(iFile)%ptr%defined = silja_undefined

    !---------------------------------------------------------
    !
    ! Define the grid for GrADS from the silja_grid
    !
    gfile(iFile)%ptr%silamGrid = grid

    SELECT CASE(fu_gridtype(grid))
    CASE(lonlat)
      CALL lonlat_grid_parameters(grid, &
            & gfile(iFile)%ptr%ggrid%x_start, gfile(iFile)%ptr%ggrid%y_start, &
            & if_corner_in_geo_coord,&
            & gfile(iFile)%ptr%ggrid%nx, gfile(iFile)%ptr%ggrid%ny, &
            & pole_x, pole_y, & 
            & gfile(iFile)%ptr%ggrid%x_step, gfile(iFile)%ptr%ggrid%y_step)

      call msg('open_gradsfile_o reports lonlat grid:')
      call msg('GrADS grid xStart:', gfile(iFile)%ptr%ggrid%x_start)
      call msg('GrADS grid yStart:', gfile(iFile)%ptr%ggrid%y_start)
      call msg('GrADS grid nx:', gfile(iFile)%ptr%ggrid%nx)
      call msg('GrADS grid ny:', gfile(iFile)%ptr%ggrid%ny)
      call msg('GrADS grid xStep:', gfile(iFile)%ptr%ggrid%x_step)
      call msg('GrADS grid yStep:', gfile(iFile)%ptr%ggrid%y_step)

    case (anygrid)
      call grid_dimensions(grid, gfile(iFile)%ptr%ggrid%nx, gfile(iFile)%ptr%ggrid%ny)
      gfile(iFile)%ptr%ggrid%y_start = fu_lat_geographical_from_grid(1., 1., grid)
      gfile(iFile)%ptr%ggrid%x_start = fu_lon_geographical_from_grid(1., 1., grid)
      gfile(iFile)%ptr%ggrid%x_step = fu_dx_cell_deg(grid, 1, 1)
      gfile(iFile)%ptr%ggrid%y_step = fu_dy_cell_deg(grid, 1, 1)
      call msg('open_gradsfile_o reports anygrid:')
      call msg('GrADS grid nx:', gfile(iFile)%ptr%ggrid%nx)
      call msg('GrADS grid ny:', gfile(iFile)%ptr%ggrid%ny)

    CASE DEFAULT
      CALL set_error('Strange grid type','open_grads_file_o')
      RETURN
    END SELECT
    IF(error)RETURN


    NbrOfPoints = gfile(iFile)%ptr%ggrid%nx * gfile(iFile)%ptr%ggrid%ny
    IF(NbrOfPoints < 1)THEN
      CALL set_error('Problem with grid dimensions','open_grads_file_o')
      RETURN
    END IF

    gfile(iFile)%ptr%ggrid%defined = fu_set_true()

    !---------------------------------------------------------
    !
    ! Open files. Attention: NbrOfPoints has to be translated to the 
    ! number of bytes for the binary file, which is machine-dependent
    !
    if(.not.if_parallel) gfile(iFile)%ptr%unit_bin = fu_next_free_unit() !------------ binary file
    if(len_trim(dir) > 0)then
      gfile(iFile)%ptr%fname = fu_connect_strings(dir,dir_slash,fname) !fname
    else
      gfile(iFile)%ptr%fname = fname
    endif
    !
    ! Store the fname - that will be the ctl with full length of the forecast
    !
    gfile(iFile)%ptr%fname_initial = gfile(iFile)%ptr%fname

    if(present(chTemplate))then
      gfile(iFile)%ptr%chTemplate = chTemplate
    else
      gfile(iFile)%ptr%chTemplate = gfile(iFile)%ptr%fname
    endif

    ! Note that currently the mpi version will exit if the file creation fails (low-level
    ! mpi parts are not able to set correct error flags)
    call check_create_dir(gfile(iFile)%ptr%fname)
    if (error) return

    if(if_parallel)then
      call smpi_open_gradsfile_mpiio(gfile(iFile)%ptr%fname, gfile(iFile)%ptr%unit_bin)
      gfile(iFile)%ptr%ifMPIIO = .true.
    else if (gfile(ifile)%ptr%ifBuffered) then
      ! direct access not used
      call open_binary(gfile(iFile)%ptr%unit_bin, gfile(iFile)%ptr%fname, recl=int_missing, &
                     & action='write', access='stream', status='replace')
    else
      CALL open_grads_binary_o(gfile(iFile)%ptr%fname, gfile(iFile)%ptr%unit_bin, NbrOfPoints)
      gfile(iFile)%ptr%ifMPIIO = .false.
    end if

    open_gradsfile_o = iFile
      
  END FUNCTION open_gradsfile_o  

  !************************************************************************************

  subroutine switch_grads_binary_o(Findex, dir, fname, iInvert, now)
    !
    ! Closes currently open GrADS file pointed by the index and opens the new GrADS 
    ! binary file.
    ! 
    ! Be careful to use MPI-aware routines when needed.
    !
    ! Author: M.Sofiev, FMI
    !
    IMPLICIT NONE
    !
    ! Imported paraneters with intent IN
    integer, intent(in) :: Findex ! File index in the array of grads structures
    CHARACTER (LEN=*), INTENT(in) :: dir, fname ! Directory and name of the GrADS file
    ! 0 - do nothing, 1 - invert, 2 - invert and change start time in gFile structure:
    integer, intent(in) :: iInvert 
    type(silja_time), intent(in) :: now

    ! Local declarations
    INTEGER :: NbrOfPoints ! in the grid
    logical :: dir_exists
    character(len=*), parameter :: sub_name = 'switch_grads_binary_o'

    ! If needed (adjoint run), have the leading MPI rank invert the previous binary.
    !
    if(iInvert > 0 .and. (.not. smpi_is_mpi_version() .or. smpi_adv_rank == 0)) then
      ! If we're running some fancy MPI configuration, then we should be checking within
      ! the ensemble, or within the IO group, or something...
      if (smpi_is_mpi_version()) then
        if (fu_fails(smpi_io_rank == 0, 'Can''t invert binary in this mpi conf', sub_name)) return
        if (fu_fails(smpi_ens_rank == 0, 'Can''t invert binary in this mpi conf', sub_name)) return
      end if
      call invert_grads_binary(Findex, iInvert)
    end if

    if (error) return
    
    if (gfile(Findex)%ptr%ifBuffered) call flush_buffer(findex)

    if (gfile(Findex)%ptr%ifMPIIO) then
      call smpi_close_gradsfile_mpiio(gfile(Findex)%ptr%unit_bin)
    else
      close(gfile(Findex)%ptr%unit_bin) 
    end if
    !
    ! Seize this moment to write the ctl and super_ctl for the current file. Note the number of time steps!
    !
    if (.not. gfile(Findex)%ptr%ifMPIIO .or. smpi_adv_rank==0)then
      call write_ctl_file(gfile(Findex)%ptr, '', .true.)
    endif
    gfile(Findex)%ptr%n_times_bin = 0
    gfile(Findex)%ptr%gtime%start_bin = now

    NbrOfPoints = gfile(Findex)%ptr%ggrid%nx * gfile(Findex)%ptr%ggrid%ny
    
    IF(NbrOfPoints < 1)THEN
      CALL set_error('Problem with grid dimensions', sub_name)
      RETURN
    END IF

    !---------------------------------------------------------
    !
    ! Open files. Attention: NbrOfPoints has to be translated to the 
    ! number of bytes for the binary file, which is machine-dependent
    !
    gfile(Findex)%ptr%unit_bin = fu_next_free_unit() !------------ binary file
    if(len_trim(dir) > 0)then
      gfile(Findex)%ptr%fname = fu_connect_strings(dir,dir_slash,fname) !fname
    else
      gfile(Findex)%ptr%fname = fname
    endif
    
    call check_create_dir(gfile(Findex)%ptr%fname) 
    if (error) return

    if (gfile(Findex)%ptr%ifMPIIO) then
      call smpi_open_gradsfile_mpiio(gfile(Findex)%ptr%fname, gfile(Findex)%ptr%unit_bin)
    else
      CALL open_grads_binary_o(gfile(Findex)%ptr%fname, gfile(Findex)%ptr%unit_bin, NbrOfPoints)
    end if
    if (error) return
    
    gfile(FIndex)%ptr%rec = 1

  END subroutine switch_grads_binary_o

  subroutine check_create_dir(filename)
    ! Check if filename is in a non-existent directory, and create if needed.
    implicit none
    character(len=*), intent(in) :: filename

    character(len=len(filename)) :: dirname
    logical :: dir_exists

    dirname = fu_dirname(filename)
    if (dirname /= '') then
      inquire(file=dirname, exist=dir_exists)
      if (.not. dir_exists) call create_directory_tree(dirname)
    end if
    
  end subroutine check_create_dir

  !******************************************************************************

  subroutine read_field_from_grads_id(igf, field_id, grid_data, fill_value_, direction_)
    !
    ! Reads the requested field from the GrADS file.
    ! ATTENTION. To get the field, one has to order it. This is entirely different
    ! from the GRIB file processing where the field location is unknown. Here we can
    ! do the job properly by reading exactly those fields, which we need
    !
    implicit none

    ! Imported parameters
    integer, intent(in) :: igf
    TYPE(silja_field_id), INTENT(in) :: field_id ! Field id 
    real, dimension(:), intent(inout) :: grid_data  ! must exist
    real, intent(in), optional :: fill_value_
    integer, intent(in), optional :: direction_

    ! Local variables
    integer :: iTmp, iVar, iLev, iVarRequested, iLevRequested, iRec, nTimeSteps, indTime, ix, iy
!    integer :: nPoints
    logical :: ifFound
    type(silam_sp), dimension(:), save, pointer :: fnames
    type(silja_time) :: timeFileBeg, time_in_file
    type(grads_file), pointer :: gf
    integer :: direction
    real :: fill_value

    fill_value = real_missing
    if(present(fill_value_)) fill_value = fill_value_

    direction = forwards
    if(present(direction_)) direction = direction_

    !
    ! Find the field index in the GrADS file
    !
    gf => gfile(igf)%ptr

    iVarRequested = 1
    ifFound = .false.
    do while(iVarRequested <= gfile(igf)%ptr%n_vars)
      if(fu_quantity(field_id) == gfile(igf)%ptr%gvars(iVarRequested)%quantity)then
        if(fu_species(field_id) == gfile(igf)%ptr%gvars(iVarRequested)%species)then
!        if(fu_str_u_case(fu_substance_name(field_id)) == &
!                               & fu_str_u_case(gfile(igf)%ptr%gvars(iVarRequested)%chSubstNm))then
!          if(fu_mode(field_id) == gfile(igf)%ptr%gvars(iVarRequested)%aerosolMode)then
!            if(fu_optical_wave_length(field_id) .eps. &
!                                    & gfile(igf)%ptr%gvars(iVarRequested)%fOpticalWaveLength)then
              iLevRequested = 1

              if(gfile(igf)%ptr%gvars(iVarRequested)%n_levs > 1)then

                do while(iLevRequested <= gfile(igf)%ptr%n_levs)
                  iTmp = fu_glevel_type(fu_leveltype(fu_level(field_id)))
                  if(fu_get_glevel(fu_level(field_id),iTmp) .eps. &
                                               & gfile(igf)%ptr%glevs%levels(iLevRequested))then
                    ifFound = .true.
                    exit
                  endif
                  iLevRequested = iLevRequested + 1
                enddo    ! cycle over levels for the multi-level variable

              else
                ifFound = .true.
              end if    ! if single-level variable

              if(ifFound)exit
!            endif   ! optical wave length is the same
!          endif   ! mode size is the same
!        endif   ! substance is the same
        endif   ! species is the same
      endif  ! quantity is right
      iVarRequested = iVarRequested + 1
    end do  ! cycle over GrADS variables

    if(.not. ifFound)then
      call msg_warning('Field not found','read_field_from_grads_id')
      call report(field_id)
      call set_error('Field not found','read_field_from_grads_id')
      return
    else
!call msg("Found field with valid time:" + fu_str(fu_valid_time(field_id)))
    endif
!    !FIXME
!call msg("=============================")
!call msg("Getting field:")
!call report(field_id)
!call msg("=============================")

    !
    ! Determine the right time stamp among those available in the file
    ! Note that the id can be either generated by metadata generator in this module or created
    ! elsewhere. It may not follow the rules of validity periods - but we must stay in agreement 
    ! with the metadata generator.
    !
    select case (gfile(igf)%ptr%data_time_features)
!!!      case(dynamic_map)
!!!        if(defined(fu_validity_length(field_id)))then
!!!          if(fu_validity_length(field_id) > one_second)then
!!!            !
!!!            ! If the field is period-valid, a chance for ambiguity exists, have to check it
!!!            !
!!!            ix = fu_grads_time_index(igf, fu_valid_time(field_id), forwards, .false.)
!!!            iy = fu_grads_time_index(igf, fu_valid_time(field_id) + fu_validity_length(field_id), forwards, .false.)
!!!            if(ix == int_missing)ix = iy   ! Stupid, I know
!!!            if(iy == int_missing)iy = ix
!!!!call msg('Indices for valid and valid+period:',ix,iy)
!!!            if(ix /= iy)then
!!!              if(abs(ix-iy)>1)then
!!!                call msg('Period-valid field requested allowing many grads times. From/to:', ix, iy)
!!!                time_in_file = time_missing
!!!              else
!!!              !
!!!              ! Ambiguity: the validity length allows two grads moments to fall into the interval. 
!!!              ! Most-probably, these are just two ends of the interval but resolving the problem is
!!!              ! not possible without knowledge of direction of the calculations. 
!!!              ! The rule is: take the index that is the closest to the "end" of the interval, i.e. the
!!!              ! latest time for forward run but the earliest time for inverse run.
!!!              !
!!!              if(present(direction))then
!!!                if(direction == forwards)then
!!!                  time_in_file = fu_time_of_grads(igf, iy)
!!!                else
!!!                  time_in_file = fu_time_of_grads(igf, ix)
!!!                endif
!!!              else
!!!                  !
!!!                  ! The last resort: if the label shows at the ambiguity point, let's take this very moment,
!!!                  ! essentially reducing the length of validity a little bit.
!!!                  !
!!!                call msg_warning('Period-valid field requested allowing more than one grads time', &
!!!                             & 'read_field_from_grads_id')
!!!                  select case(gfile(igf)%ptr%time_label_position)
!!!                    case (start_of_period)     ! take the latest
!!!                      time_in_file = fu_time_of_grads(igf, ix)
!!!                      call msg('Resolve start_of_period. Took:' + fu_str(time_in_file), ix)
!!!                    case (end_of_period)
!!!                      time_in_file = fu_time_of_grads(igf, iy)
!!!                      call msg('Resolve end_of_period. Took:' + fu_str(time_in_file), iy)
!!!                    case default
!!!                      call msg_warning('Cannot resolve the ambiguity','read_field_from_grads_id')
!!!                      time_in_file = time_missing
!!!                  end select
!!!                endif  ! present direction
!!!              endif  ! |ix-iy| > 1
!!!            else
!!!              time_in_file = fu_time_of_grads(igf, ix)
!!!            endif
!!!          else
!!!            time_in_file = fu_time_of_grads(igf, fu_grads_time_index(igf, fu_valid_time(field_id), forwards,.false.))
!!!          endif
!!!        else
!!!          time_in_file = fu_time_of_grads(igf, fu_grads_time_index(igf, fu_valid_time(field_id), forwards,.false.))
!!!        endif  ! defined validity length
      
      case(dynamic_map)
        iTmp = fu_grads_time_index(igf, fu_valid_time(field_id),  direction, .false.)

      case(monthly_climatology)
        iTmp = fu_grads_time_index(igf, fu_valid_time(field_id),  direction, .true.)

      case(static_climatology)
         iTmp = 1

      case default
        call set_error('Unknown grads time featires:' + fu_str(gfile(igf)%ptr%data_time_features), &
                     & 'read_field_from_grads_id')
        return
    end select
    if (iTmp > 0) then
      time_in_file = fu_time_of_grads(igf, iTmp)
    else
      call msg("index of fu_valid_time(field_id) "+fu_str(fu_valid_time(field_id)), iTmp)
      call set_error("Strange index from grads file","read_field_from_grads_id")
    endif
    !
    ! The above clumsy logic can fail for many reasons
    !
    if(error .or. .not. defined(time_in_file))then
      call report(field_id)
      call msg('Possible times and grads-file time indices:' + &
                         & fu_str(fu_time_of_grads(igf,ix))+','+fu_str(fu_time_of_grads(igf,iy)), ix, iy)
      call set_error('Failed search for the time in grads file', 'read_field_from_grads_id')
      return
    endif

! call msg('grads: time to search:' + fu_str(time_in_file))
    !
    ! Get the right binary file name. Note that static-field file may not have time step
    !
    call FNm_from_single_template(gfile(igf)%ptr%grTemplate, &
                                & time_in_file, &
                                & fnames, &
                                & anal_time = fu_analysis_time(field_id), &
                                & ifStrict = .false., &
                                & ifAdd = .false., &  ! to the fnames array
                                & ifAllowZeroFcLen = .true., &
                                & ifWait = .false.)
    if(error)return

!    !FIXME
!    Call msg("Grades file:" + fnames(1)%sp)
    !
    ! There must be just one file
    !
    if(size(fnames) > 1)then
      if(len_trim(fnames(2)%sp) > 0)then
        call msg_warning('More than one GrADS file contains the field','read_field_from_grads_id')
        do iVar = 1, size(fnames)
          call msg('File name:' + fnames(iVar)%sp)
        end do
        call set_error('More than one GrADS file contains the field','read_field_from_grads_id')
        return
      endif
    endif

    !!!nPoints = INT(REAL(gfile(igf)%ptr%ggrid%nx) * REAL(gfile(igf)%ptr%ggrid%ny) + 0.0001)
    !!!
    !!!!
    !!!! Open it if needed
    !!!!
    !!!if(fnames(1)%sp /= gfile(igf)%ptr%fname)then
    !!!
    !!!  if(gfile(igf)%ptr%unit_bin < 0)then
    !!!    gfile(igf)%ptr%unit_bin = fu_next_free_unit()
    !!!  else
    !!!    close(gfile(igf)%ptr%unit_bin, stat = iTmp)  ! attempt tp close but if fails so be it
    !!!  endif
    !!!  call open_grads_binary_i(fnames(1)%sp, gfile(igf)%ptr%unit_bin, nPoints, &
    !!!                         & gfile(igf)%ptr%ifBigEndian)
    !!!  if(error)return
    !!!  gfile(igf)%ptr%fname = fnames(1)%sp
    !!!endif

    nTimeSteps = nint((time_in_file - gfile(igf)%ptr%gtime%start) / gfile(igf)%ptr%gtime%step)
    !
    ! Having the indices defined, get the field
    !
    call read_field_from_grads_indices(igf, iVarRequested, iLevRequested, &
                                     & nTimeSteps+1, &
                                     & grid_data, fill_value)

  end subroutine read_field_from_grads_id


  !*******************************************************************
  
  subroutine read_field_from_grads_indices(igf, indVar_, indLev_, indTime_, grid_data, fill_value_)
    !
    ! Reads the requested field from the GrADS file using exact indices
    ! ATTENTION. To get the field, one has to order it. This is entirely different
    ! from the GRIB file processing where the field location is unknown. Here we can
    ! do the job properly by reading exactly those fields, which we need.
    !
    implicit none

    ! Imported parameters
    integer, intent(in) :: igf, indVar_, indLev_, indTime_
    real, dimension(:), intent(out) :: grid_data  ! must exist
    real, intent(in), optional :: fill_value_

    ! Local variables
    integer :: iRec, iVTmp, iTTmp, iLTmp, nPoints, ix, iy, indTimeLocal
    real :: fill_value
    type(silam_sp), dimension(:), save, pointer :: fnames
    type(grads_file), pointer :: gf
    type(grads_variable), pointer :: gvar

    fill_value = real_missing
    if(present(fill_value_))fill_value = fill_value_



    gf => gfile(igf)%ptr
    nPoints = gf%ggrid%nx * gf%ggrid%ny

    
    iRec = 1 + gf%first_field_offset_in_tstep(indVar_)
    if(gf%ifZInverse)then
      iRec = iRec + gf%gvars(indVar_)%n_levs - indLev_
    else
      iRec = iRec + indLev_ - 1
    endif
    !!!IRec is now number of the field in the timestep

    if (gf%ifBuffered) then
      grid_data(1:nPoints) = gf%gradsbuf%buf(1:nPoints,iRec)
    else

      !
      ! Get the right binary file name. Note that static-field file may not have time step
      !
      if(error)return
      call FNm_from_single_template(gf%grTemplate, &
                                  & fu_time_of_grads(igf, indTime_), &
                                  & fnames, &
                                  & anal_time = fu_time_of_grads(igf, 1), &  ! analysis time
                                  & ifStrict = .false., &
                                  & ifAdd = .false., &  ! to the fnames array
                                  & ifAllowZeroFcLen = .true., &
                                  & ifWait = .false.)
      if(error)return
      !
      ! Open it if needed
      !
      if(fnames(1)%sp /= gf%fname)then
        if(gf%unit_bin < 0)then
          gf%unit_bin = fu_next_free_unit()
        else
          close(gf%unit_bin, iostat = ix)
        endif
        call open_grads_binary_i(fnames(1)%sp, gf%unit_bin, nPoints, &
                               & gf%ifBigEndian)
        if(error)return
        gf%fname = fnames(1)%sp
      endif
      !
      ! Current binary starts from some time step
      !
      indTimeLocal = fu_time_index_in_grads_binary(igf, indTime_)

      iRec = iRec + (indTimeLocal - 1)*gf%first_field_offset_in_tstep(gf%n_vars+1)
      !!gf%first_field_offset_in_tstep(gf%n_vars+1) -- number of fields per timestep
      ! iRec is number of the field in the binary file
      !
      !
      ! Finally, READ, allowing for flipped y axis
      !
      if(fu_fails(size(grid_data) >= nPoints,'Too small array for field:' + fu_str(size(grid_data)) + &
                                              & ',' + fu_str(nPoints),'read_field_from_grads_id'))return
      if(gf%ifYInverse)then
        read(gf%unit_bin,rec=iRec, iostat=iVTmp) &
                   & ((grid_data(ix+(iy-1)*gf%ggrid%nx), &
                                 & ix = 1,gf%ggrid%nx), iy = gf%ggrid%ny,1,-1)
      else
        read(gf%unit_bin,rec=iRec, iostat = iVTmp)(grid_data(ix), ix = 1, nPoints)
      endif
      if(iVTmp /= 0)then
        call msg('failed to read the record,size,',iRec, nPoints)
        call msg('iostat = ', iVTmp)
        call msg('file unit: ', gf%unit_bin)
        call report_open_files()
        call set_error('failed to read the record','read_field_from_grads_id')
      endif
    endif !gf%buffered

    !!!Need for scaling and missing_value handling
    if (gf%missing_value /= fill_value .or. gf%gvars(indVar_)%fScalingFactor /= 1.) then
      !Replace the missing values
      iRec = 0
      do ix = 1, nPoints
        if(grid_data(ix) == gf%missing_value)then
          grid_data(ix) = fill_value
          iRec = iRec + 1
        else
          grid_data(ix) = grid_data(ix) * gf%gvars(indVar_)%fScalingFactor
        endif
      enddo
      if (iRec > 0)then
        call msg(fu_quantity_string(gf%gvars(indVar_)%quantity) + ':' + fu_str(iRec) + &
               & '- missing values found and replaced with', fill_value)
  !      call report(field_id)
      endif
    endif
!call msg('iRec,sum(value)',iRec,sum(grid_data(1:nPoints)))

    
  contains
  
    !===================================================
  
    subroutine report_open_files()
      implicit none
      integer :: ind_file
      character(len=fnlen) :: filename, access, form
      logical :: opened
      ind_file = 8

      do ind_file = 1, 20
        inquire(unit=ind_file, opened=opened)
        if (ind_file > 20) exit
        if (.not. opened) then
          call msg('Not open: ', ind_file)
          cycle
        end if
        !if (.not. opened) then
        !  exit
        inquire(unit=ind_file, access=access, name=filename, form=form)
        !inquire(unit=ind_file, name=filename)
        
        call msg('File unit: ', ind_file)
        call msg('  access: ' // trim(access))
        call msg('  form: ' // trim(form))
        call msg('  file: ' // trim(filename))
        
        !ind_file = ind_file + 1
      end do
  
    end subroutine report_open_files

    !========================================================================
    
    integer function fu_time_index_in_grads_binary(igf, indTimeGlob)
      !
      ! Computes time index in a specific grads binary from its lgobal time index
      !
      implicit none
      
      ! Imported parameters
      integer, intent(in) :: igf, indTimeGlob
    
      ! Local variables
      type(silja_time) :: timeFileBeg, time_in_file
      integer :: nTimeSteps, iTmp
      
      !
      ! Get the start time for the specific binary that contains the given global time index
      !
      time_in_file = fu_time_of_grads(igf, indTimeGlob)
      timeFileBeg = fu_same_template_start_time(gfile(igf)%ptr%grTemplate, &
                                              & time_in_file, &
                                              & gfile(igf)%ptr%gtime%start)
!call msg("timeFileBeg:"+ fu_str(timeFileBeg))
!call msg("time_in_file:"+ fu_str(time_in_file))
!call msg("gfile(igf)%ptr%gtime%start:"+ fu_str(gfile(igf)%ptr%gtime%start))
!call msg("gfile(igf)%ptr%gtime%step:"+ fu_str(gfile(igf)%ptr%gtime%step))
      ! 
      ! if this is the first binary file, and if the run did not start at midnight, the first time
      ! could be later than the first possible by template:
      !
      if (timeFileBeg < gfile(igf)%ptr%gtime%start) timeFileBeg = gfile(igf)%ptr%gtime%start

      if(error .or. .not.defined(timeFileBeg))return

      if(gfile(igf)%ptr%gtime%ifVaryingStep)then
        nTimeSteps = -1
        do iTmp = 1, gfile(igf)%ptr%n_times
          if(time_in_file == gfile(igf)%ptr%gtime%arTimes(iTmp))then
            nTimeSteps = iTmp-1
            exit
          endif
        enddo
        if(nTimeSteps < 0)then
          call msg_warning('Failed to find the following time in GrADS file records:' + &
                         & fu_time_to_io_string(time_in_file), 'fu_time_index_in_grads_binary')
          call msg('Times available from GrADS file:')
          do iTmp = 1, gfile(igf)%ptr%n_times
            call report(gfile(igf)%ptr%gtime%arTimes(iTmp))
          enddo
          call set_error('Failed to find the following index/time in GrADS fiel records:' + &
                       & fu_str(indTimeGlob) + ',' + &
                       & fu_time_to_io_string(time_in_file), 'fu_time_index_in_grads_binary')
          return
        endif
      else
        nTimeSteps = nint((time_in_file - timeFileBeg) / gfile(igf)%ptr%gtime%step)
      endif
      
      fu_time_index_in_grads_binary = nTimeSteps + 1

    end function fu_time_index_in_grads_binary
    
  end subroutine read_field_from_grads_indices
  
  
  ! ****************************************************************


  SUBROUTINE  wr_n_fld_2_gf_from_id_and_data (igf, field_id, grid_data, forced_valid_time, &
                                            & if_regular_output_times)
    !
    ! Checks and writes the requested field to the GrADS file with
    ! appropriate filling of the ctl file structures, if needed.
    ! Field must follow the rules outlined at the top of the module.
    !
    ! Specific tricks have to be applied to the 2D variables, which are 
    ! all treated here as those with "no levels", i.e. having the number 
    ! of levels 0 (see GrADS manual). Thus, for such a variables no level
    ! comparison is performed, they do not affect the number of layers in
    ! the file, but such a variable may be written only once per time period
    ! This makes the routine EXTREMELY CLUMSY, but be careful in corrections:
    ! any possible combinations have to be considered
    !
    ! Author: M.Sofiev
    !
    IMPLICIT NONE

    ! Imported variables with intent IN
    TYPE(silja_field_id), INTENT(in) :: field_id ! Field id 
    real, dimension(:), intent(in) :: grid_data
    INTEGER, INTENT(in) :: igf      ! Index of the GrADS File
    type(silja_time), intent(in), optional :: forced_valid_time
    logical, intent(in), optional :: if_regular_output_times

    ! Local variables
    INTEGER :: nx,ny
    integer(kind=4) :: nPoints
    REAL :: x_corner, y_corner, x_step, y_step, pole_x, pole_y
    LOGICAL :: if_corner_in_geo_coord, if_south_pole, ifTimeOK, ifLevelOK
!    TYPE(silja_field), POINTER :: fieldptr ! Field to write to the file
    TYPE(grads_file), POINTER :: gf ! Pointer to the selected GrADS file

!    fieldptr => field

    if(igf == int_missing)then
      call set_error('Undefined index of GrADS file','wr_n_fld_2_gf_from_id_and_data')
      return
    endif
    if(igf < 1 .or. igf > size(gfile))then
      call msg('Index of GrADS file:',igf)
      call set_error('Strange index of GrADS file','wr_n_fld_2_gf_from_id_and_data')
      return
    endif

    gf => gfile(igf)%ptr


!    call msg('')
!    call msg('Field for the GrADS file:' + &
!           & fu_quantity_string(fu_quantity(field_id)) + fu_str(fu_species(field_id)))
!    call report(field_id)
!    call msg('------------------GrADS structure so far:----------------------------')
!    call report(gf)
!    call msg('')



    !----------------------------------------------------
    !
    ! Check trivial stuff: that the file exists and the grids match
    !
    IF(gf%unit_bin < 1)THEN
      CALL set_error('Wrong GrADS file index',' wr_n_fld_2_gf_from_id_and_data')
      RETURN
    END IF
      
    SELECT CASE(fu_gridtype(fu_grid(field_id))) 
    CASE(lonlat)
      CALL lonlat_grid_parameters(fu_grid(field_id), &
                       & x_corner, y_corner, if_corner_in_geo_coord,&
                       & nx, ny, &
                       & pole_x, pole_y, & 
                       & x_step, y_step)

      IF(error)RETURN
      IF(.not.((x_corner .eps. gf%ggrid%x_start) .and. &
             & (y_corner .eps. gf%ggrid%y_start) .and. &
             & (nx == gf%ggrid%nx) .and. &
             & (ny == gf%ggrid%ny) .and. &
             & (x_step .eps. gf%ggrid%x_step) .and. &
             & (y_step .eps. gf%ggrid%y_step))) THEN
        !
        ! Grids do not match exactly but, may be, they are Arakawa-corresponding?
        ! FORBIDDEN now. 2017.04.03
        !
!        call msg_warning('Grids do not match, trying Arakawa','wr_n_fld_2_gf_from_id_and_data')
        !if(.not. fu_grids_arakawa_correspond(fu_grid(field_id), &
        !                                   & fu_set_lonlat_grid ('', &
        !                                                       & gf%ggrid%x_start,&
        !                                                       & gf%ggrid%y_start, &
        !                                                       & if_corner_in_geo_coord, &
        !                                                       & gf%ggrid%nx,&
        !                                                       & gf%ggrid%ny, &
        !                                                       & fu_pole(fu_grid(field_id)), & 
        !                                                       & gf%ggrid%x_step,&
        !                                                       & gf%ggrid%y_step)))then
          call msg('Input field grid: ')
          call report(fu_grid(field_id))
          call msg('GrADS grid xStart:', gf%ggrid%x_start)
          call msg('GrADS grid yStart:', gf%ggrid%y_start)
          call msg('GrADS grid nx:', gf%ggrid%nx)
          call msg('GrADS grid ny:', gf%ggrid%ny)
          call msg('GrADS grid xStep:', gf%ggrid%x_step)
          call msg('GrADS grid yStep:', gf%ggrid%y_step)
          CALL set_error('Grids do not match',' wr_n_fld_2_gf_from_id_and_data')
          RETURN
!        endif
      END IF
    case(anygrid)
      !
      ! Just check that nx and ny are the same
      !
      call grid_dimensions(fu_grid(field_id), nx, ny)
      IF(error)RETURN
      IF(.not.((nx == gf%ggrid%nx) .and. &
             & (ny == gf%ggrid%ny)))then
        call msg('Input field grid: ')
        call report(fu_grid(field_id))
        call msg('GrADS grid xStart:', gf%ggrid%x_start)
        call msg('GrADS grid yStart:', gf%ggrid%y_start)
        call msg('GrADS grid nx:', gf%ggrid%nx)
        call msg('GrADS grid ny:', gf%ggrid%ny)
        call msg('GrADS grid xStep:', gf%ggrid%x_step)
        call msg('GrADS grid yStep:', gf%ggrid%y_step)
        CALL set_error('Grids do not match',' wr_n_fld_2_gf_from_id_and_data')
        RETURN
      endif

    case default
      call set_error('Only for lonlat  and anygrids this far','wr_n_fld_2_gf_from_id_and_data')
    end select


    !----------------------------------------------------
    !
    ! If the structures for this file are not finished - continue to fill them.
    ! Note. After fill_structures the values of the structures
    ! are set equal to the field, so checking must pass 
    !
    IF (.not.(fu_true(gf%gtime%defined).and. &
            & fu_true(gf%gvars(gf%var_nbr)%defined).and. &
            & fu_true(gf%glevs%defined))) THEN
!      print *, gf%var_nbr, &
!             & fu_true(gf%gtime%defined), &
!             & fu_true(gf%gvars(gf%var_nbr)%defined), &
!             & fu_true(gf%glevs%defined)

      if(present(forced_valid_time))then
        CALL fill_structures(field_id, gf, forced_valid_time)
      else
        CALL fill_structures(field_id, gf)
      endif
      
      IF(error)RETURN
    END IF
    !-----------------------------------------------------
    !
    ! If expected field corresponds to what has come - write, increment counters
    ! and return.
    ! Note that for irregular time steps we shall force the time tag and accept it no metter what
    !
    ifTimeOK = .not. gf%gtime%ifVaryingStep  ! for irregular output no time can be rejected

    if(.not. ifTimeOK .and. present(if_regular_output_times))then
      gf%gtime%ifVaryingStep = gf%gtime%ifVaryingStep .or. if_regular_output_times
      ifTimeOK = .not. gf%gtime%ifVaryingStep  ! for irregular output no time can be rejected
    endif  ! forced (ir)regular output

    if(.not. ifTimeOK)then
      IF(gf%time_nbr > 1)THEN
        ifTimeOK = fu_between_times(gf%gtime%start + gf%gtime%step * REAL(gf%time_nbr-1), &
                                  & fu_valid_time(field_id), &
                                  & fu_valid_time(field_id) + fu_validity_length(field_id), &
                                  & .true.)  ! if accept boundaries
      ELSE
        ifTimeOK = fu_between_times(gf%gtime%start, &           ! Time
                                  & fu_valid_time(field_id), &  ! Limit 1
                                  & fu_valid_time(field_id) + fu_validity_length(field_id), & !Lim2
                                  & .true.)  ! if accept boundaries
      END IF
    endif  ! if forced timeOK

    if(.not.ifTimeOK) call set_error('Time mis-matched','wr_n_fld_2_gf_from_id_and_data')
    ! Check time
    IF(ifTimeOK)THEN

      ! Check quantity
      IF(gf%gvars(gf%var_nbr)%quantity == fu_quantity(field_id) .and. &
      & gf%gvars(gf%var_nbr)%species == fu_species(field_id) .and. &
      & gf%gvars(gf%var_nbr)%chCocktailNm == fu_cocktail_name(field_id))THEN

        ! Check level
        IF(gf%gvars(gf%var_nbr)%iVerticalFeature /= level_3d_type_flag)THEN
          ifLevelOK = .true.
        ELSE
          ifLevelOK = gf%glevs%levels(gf%lev_nbr) == &
                    & fu_get_glevel(fu_level(field_id),gf%levType)
          if(.not.ifLevelOK)then
            !
            ! May be, this variable is 2d with just potential for being 3D in 
            ! some other cases ? Necesary condition for that - current level number
            ! must be 1. There is no sufficient condition to check. If 2D assumption 
            ! is wrong - the variable will fail at the next cycle.
            !
            if(gf%lev_nbr == 1)then
              ifLevelOK = .true.
              gf%gvars(gf%var_nbr)%iVerticalFeature = level_2d_type_flag
            endif
          endif
        END IF
        if(.not.ifLevelOK) call set_error('Level mis-matched','wr_n_fld_2_gf_from_id_and_data')
        IF(ifLevelOK)THEN

          !----------------------------------------------------------------
          !
          !                      WRITE THE FIELD 
          ! and set defined and increment counters
          !
          gf%n_times=gf%time_nbr
          !
          nPoints = gf%ggrid%nx * gf%ggrid%ny
          if (gf%ifBuffered) then
            if (fu_buffer_full(gf%gradsbuf)) call flush_buffer(igf)
            if (error) return
            call field_to_buffer(gf%gradsbuf, grid_data, nPoints)
            ! record counter incremented by flush_buffer
          else
            if(gf%ifMPIIO)then
              call smpi_write_grads_field_parallel(grid_data, nPoints, gf%rec, gf%unit_bin)
            else
              CALL write_grads_field(grid_data, nPoints, gf%rec, gf%unit_bin)  ! Writing
            end if
            gf%rec=gf%rec+1
          end if
          if (error) return
          gf%defined = silja_true

          !
          ! Increment counters 
          !

          IF(gf%gvars(gf%var_nbr)%iVerticalFeature == level_3d_type_flag) THEN
            !
            ! 3D variable
            !
            gf%lev_nbr = gf%lev_nbr + 1  
            IF(gf%time_nbr == 1) RETURN !- For the first time period lists are incomplete
!            IF(gf%glevs%levels(gf%lev_nbr) < 0.)THEN
            IF(gf%lev_nbr > gf%n_levs)THEN
              gf%lev_nbr=1             !--- The next variable
              gf%var_nbr=gf%var_nbr+1
              IF(gf%var_nbr > gf%n_vars)THEN ! May be, next time period
!              IF(gf%gvars(gf%var_nbr)%quantity < 0)THEN ! May be, next time period
                gf%var_nbr=1
                gf%time_nbr=gf%time_nbr+1
                gf%n_times_bin = gf%n_times_bin + 1
!      call msg('New time period2:',gf%n_times_bin)
!      call msg('New time period2 total time:',gf%time_nbr)
              END IF
            END IF
          ELSE
            !
            ! One of 2D-type variables
            !
            gf%var_nbr = gf%var_nbr + 1
            gf%lev_nbr = 1
            IF(gf%time_nbr == 1) RETURN !- For the first time period lists are incomplete
            IF(gf%gvars(gf%var_nbr)%quantity < 0)THEN ! The next time period
              gf%var_nbr = 1
              gf%time_nbr=gf%time_nbr + 1
              gf%n_times_bin = gf%n_times_bin + 1
!      call msg('New time period1:',gf%n_times_bin)
!      call msg('New time period1 total time:',gf%time_nbr)
            END IF
          END IF ! if2D variable

          RETURN  ! Normal return from this sub

        END IF ! If levels match
      else
           call set_error("Quantity or species mismatched", &
                        & "wr_n_fld_2_gf_from_id_and_data")
      END IF ! If quantities match
    END IF ! If time periods match

    call msg('')
    call msg('Field came:')
    call report(field_id)
    call msg('Field expected:')
    IF(gf%time_nbr > 1)THEN
      call report(gf%gtime%start + gf%gtime%step * REAL(gf%time_nbr-1))
    ELSE
      call report(gf%gtime%start)
    END IF

    IF(gf%gvars(gf%var_nbr)%iVerticalFeature == level_3d_type_flag)THEN
      call msg('Level value:', gf%glevs%levels(gf%lev_nbr))
    ELSEIF(gf%gvars(gf%var_nbr)%iVerticalFeature == integrate_column_flag)THEN
      call msg('Column-integrated var')
    ELSEIF(gf%gvars(gf%var_nbr)%iVerticalFeature == lowest_level_flag)THEN
      call msg('Lowest-level var')
    ELSEIF(gf%gvars(gf%var_nbr)%iVerticalFeature == level_2d_type_flag)THEN
      call msg('2D var')
    ELSE
      call set_error('Unknown vertical features of the variable' &
         & + fu_str(gf%gvars(gf%var_nbr)%iVerticalFeature),'wr_n_fld_2_gf_from_id_and_data')
    END IF
    call msg(fu_quantity_string(gf%gvars(gf%var_nbr)%quantity) + '_' + &
           & fu_str(gf%gvars(gf%var_nbr)%species))

    CALL set_error('Field does not match',' wr_n_fld_2_gf_from_id_and_data')

  END SUBROUTINE  wr_n_fld_2_gf_from_id_and_data


  ! ****************************************************************


  SUBROUTINE wr_nxt_field_to_gf_from_fld(igf, field)
    !
    ! Actually, refers to the above routine and set here ONLY for backward compatibility
    !
    implicit none

    TYPE(silja_field), INTENT(in), target :: field ! Field to write
    integer, intent(in) :: igf

    ! Local variables
    TYPE(silja_field), pointer:: fPtr ! Field to write

    fPtr => field
    call wr_n_fld_2_gf_from_id_and_data(igf,fu_id(fPtr),fu_grid_data(fPtr))

  END SUBROUTINE wr_nxt_field_to_gf_from_fld


  ! ****************************************************************


  ! ****************************************************************


  SUBROUTINE fill_structures(field_id, gf, forced_valid_time)
    !
    ! Fills-in the gfile structures during the first time step
    ! Special attention is paid to the 2D variables
    !
    ! Author: M.Sofiev, FMI
    !
    IMPLICIT NONE
    !
    ! Imported parameters with intent IN
    TYPE(grads_file), INTENT(inout) :: gf
    TYPE(silja_field_id), intent(in):: field_id
    type(silja_time), intent(in), optional :: forced_valid_time

    ! Local variables
    INTEGER :: iTmp
    LOGICAL :: ifSurface, ifTimeStart
    !-----------------------------------------------------
    ! If something has not matched during the second time period => error
    !
    IF(gf%time_nbr > 1)THEN
      CALL set_error('Wrong field','fill_structures')
      RETURN
    END IF

!    call report(field_id)

    !-----------------------------------------------------
    !
    ! If totally new GrADS file - start all lists.
    ! Note. If current field describes a 2D variable - there will be no levels
    !
    IF(gf%gtime%defined == silja_false)THEN ! Totally new GrADS file
      gf%gtime%start = fu_valid_time(field_id)
      gf%gtime%defined = silja_undefined
      gf%time_nbr = 1
      gf%n_times = 1
      gf%gtime%start_bin = gf%gtime%start

      gf%gvars(1)%quantity = fu_quantity(field_id)
      gf%gvars(1)%species = fu_species(field_id)
      gf%gvars(1)%chCocktailNm = fu_cocktail_name(field_id)
      gf%gvars(1)%SilamLevel = fu_level(field_id)    ! might or might not be used: can be 3D grads var
      gf%gvars(1)%n_levs = 1
      gf%gvars%defined = silja_undefined
      gf%var_nbr = 1
      gf%n_vars = 1

      call set_missing(gf%silamVertical, .true.)

      IF(fu_multi_level_quantity(gf%gvars(1)%quantity))THEN
        !
        ! Potentially possible to have several levels for this variable
        !
        if(fu_cmp_levs_eq(fu_level(field_id),entire_atmosphere_integr_level))then
          !
          ! Actually, this is column-integrated field, i.e. 2D from GrADS point of view
          !
          gf%gvars(1)%iVerticalFeature = integrate_column_flag

        elseif(fu_cmp_levs_eq(fu_level(field_id), lowest_atmosphere_level))then
          !
          ! Actually, this is lowest-level field, i.e. 2D from GrADS point of view
          !
          gf%gvars(1)%iVerticalFeature = lowest_level_flag

        else
          call set_vertical(fu_level(field_id), gf%silamVertical)
          if(error)return
          gf%glevs%levels(1) = fu_get_glevel(fu_level(field_id),gf%levType)
          gf%glevs%defined = silja_undefined
          gf%lev_nbr = 1
          gf%n_levs = 1
          gf%gvars(1)%iVerticalFeature = level_3d_type_flag ! This is a 3D variable
        endif
      ELSE
        gf%gvars(1)%iVerticalFeature = level_2d_type_flag ! This is a 2D variable, no levels
      END IF
      RETURN
    END IF  ! totally new file

    !-----------------------------------------------------
    !
    ! May be, the new field starts the second time period ?
    ! Careful here: there may be a tricky field with long validity and valid time start in past
    ! To handle this, we'll check the forced_valid_time if it is available
    !
    if(present(forced_valid_time))then
      ifTimeStart = gf%gtime%start == forced_valid_time
    else
      ifTimeStart = fu_between_times(gf%gtime%start, fu_valid_time(field_id), &
                      & fu_valid_time(field_id) + fu_validity_length(field_id), .true.)
    endif
    
    if(ifTimeStart)THEN 

      !---------------------------------------------------
      !
      ! Still the first time period => check if new variable
      ! "New" here means that the previous field was not the same. "Old" variable does NOT mean, 
      ! of course, that this is still the first one - can be the second level of the third
      ! variable
      !
      IF(gf%gvars(gf%var_nbr)%quantity == fu_quantity(field_id) .and. &
       & gf%gvars(gf%var_nbr)%species == fu_species(field_id) .and. &
        & gf%gvars(gf%var_nbr)%chCocktailNm == fu_cocktail_name(field_id))THEN
        
        if(gf%glevs%defined == silja_true)then
          !
          ! The list of levels is closed. Thus, this must be the 2nd/3rd/4th/... level
          ! of some variable. Basically, we do not have to do anything - just to check
          ! that this variable fits into the list of levels and gets correct place in
          ! the overall list of fields
          !
          if(.not. fu_get_glevel(fu_level(field_id),gf%levType) == gf%glevs%levels(gf%lev_nbr))then
            call msg('Expected level in structure:',gf%lev_nbr)
            call report(gf%silamVertical,.true.)
            call report(field_id)
            call set_error('Next variable has arrived but levels do not match','fill_structures')
            return
          endif
          gf%gvars(gf%var_nbr)%n_levs = gf%lev_nbr

        else
          !--------------------------------------
          !
          ! New level of the very first variable
          !
          IF(gf%gvars(gf%var_nbr)%iVerticalFeature /= level_3d_type_flag)THEN
            call msg('')
            call report(field_id)
            CALL set_error('2D var claims to have several levels','fill_structures')
            RETURN
          END IF
          IF(fu_get_glevel(fu_level(field_id),gf%levType) <= gf%glevs%levels(gf%lev_nbr-1).and. &
           & fu_get_glevel(fu_level(field_id),gf%levType) > gf%glevs%levels(1))THEN
            CALL set_error('Intermediate level found','fill_structures')
            RETURN
          END IF
!          call msg('Adding level to vertical:')
!          call report(gf%silamVertical)

          call add_level(gf%silamVertical, fu_level(field_id))  ! fill-in the SILAM vertical

          gf%glevs%levels(gf%lev_nbr) = fu_get_glevel(fu_level(field_id),gf%levType)
          IF(error)RETURN
          gf%gvars(gf%var_nbr)%n_levs = gf%lev_nbr
          gf%n_levs = MAX(gf%n_levs, gf%lev_nbr)
        endif  ! next multi-level variable

      ELSE

        !----------------------------------------
        !
        ! New variable => must be the first level or have a level indicator 0 
        ! (2D variable). See GrADS manual for such variables.
        !
        ! But first check if it is not in the list already.
        !
        DO iTmp = 1, gf%n_vars

          IF(gf%gvars(iTmp)%quantity == fu_quantity(field_id) .and. &
           & gf%gvars(iTmp)%species == fu_species(field_id) .and. &
           & gf%gvars(gf%var_nbr)%chCocktailNm == fu_cocktail_name(field_id))THEN
            call msg('Field that comes')
            call report(field_id)
            call msg('GrADS var N:',iTmp)
            call msg(fu_quantity_string(gf%gvars(iTmp)%quantity) + '_' + &
                   & fu_str(gf%gvars(iTmp)%species))
            CALL set_error('Duplicated variable:' + &
                         & fu_quantity_short_string(gf%gvars(iTmp)%quantity)  + '_' + &
                         & fu_str(gf%gvars(iTmp)%species), &
                         & 'fill_structures')
            RETURN
          END IF
        END DO
        !
        ! Check if the previous variable has several levels, set iVerticalFeature
        ! and close the level definition if yes. Criteria for the previous var to be 2D is:
        ! - the lev_nbr is 1 or 2. It is 1 if the variable by-definition is 2D and it is 2
        ! with levels not yet closed if the variable is theoretically 3d but appeared to be
        ! 2d now.
        !
        if(gf%lev_nbr <= 2)then
          if(gf%gvars(gf%n_vars)%iVerticalFeature == level_3d_type_flag) &
                                     & gf%gvars(gf%n_vars)%iVerticalFeature = level_2d_type_flag
        else
          gf%glevs%defined = silja_true ! Close the level definition
        endif
        !
        ! Start new variable
        !
        gf%n_vars = gf%n_vars + 1
        gf%var_nbr = gf%n_vars
        gf%lev_nbr = 1
        gf%gvars(gf%var_nbr)%quantity = fu_quantity(field_id)
        gf%gvars(gf%var_nbr)%species = fu_species(field_id)
        gf%gvars(gf%var_nbr)%chCocktailNm = fu_cocktail_name(field_id)
        gf%gvars(gf%var_nbr)%SilamLevel = fu_level(field_id)  ! might or might not be used: can be 3D grads var
        !
        ! Can it potentially be multi-level ?
        !
        if(fu_multi_level_quantity(gf%gvars(gf%var_nbr)%quantity))then
          !
          ! Potentially possible to have several levels for this variable
          !
          if(fu_cmp_levs_eq(fu_level(field_id),entire_atmosphere_integr_level))then
            !
            ! Actually, this is column-integrated field, i.e. 2D from GrADS point of view
            !
            gf%gvars(gf%var_nbr)%iVerticalFeature = integrate_column_flag

          elseif(fu_cmp_levs_eq(fu_level(field_id), lowest_atmosphere_level))then
            !
            ! Actually, this is lowest-level field, i.e. 2D from GrADS point of view
            !
            gf%gvars(gf%var_nbr)%iVerticalFeature = lowest_level_flag

          else
            gf%gvars(gf%var_nbr)%iVerticalFeature = level_3d_type_flag ! This is a 3D variable
            !
            ! Have level structure been already defined ?
            ! If yes - do nothing, the level will be checked before writing
            ! If not - start level definition
            !
            if(.not. gf%glevs%defined == silja_true)then
              call set_vertical(fu_level(field_id), gf%silamVertical)  ! fill-in the SILAM vertical
              if(error)return
              gf%glevs%levels(1) = fu_get_glevel(fu_level(field_id),gf%levType)
              gf%glevs%defined = silja_undefined
              gf%lev_nbr = 1
              gf%n_levs = 1
            endif
          endif  ! level in the field id
        else
          gf%gvars(gf%var_nbr)%iVerticalFeature = level_2d_type_flag
        endif

      END IF ! IF new variable

    ELSE

      !-------------------------------------------------------
      !
      ! The second time period starts - reset all other counters. All is defined,
      ! except for the number of steps.
      !
      ! We still need to check whether the last variable has several levels and set iVerticalFeature
      ! and close the level definition if yes. Criteria for the previous var to be 2D is:
      ! - the lev_nbr is 1 or 2. It is 1 if the variable by-definition is 2D and it is 2
      ! with levels not yet closed if the variable is theoretically 3d but appeared to be
      ! 2d now.
      !
      if(gf%lev_nbr <= 2)then      ! it is still set to the value of last variable
        if(gf%gvars(gf%n_vars)%iVerticalFeature == level_3d_type_flag) &
                                     & gf%gvars(gf%n_vars)%iVerticalFeature = level_2d_type_flag
      endif

      ! Note that the increment of time counter time_nbr is a fake offsetting the
      ! no-increment during the first writing. We have to ensure that counter always
      ! point at the field that is to come - future time, next var, next level, etc.
      !
      call msg('All GrADS variables fixed: ',gf%n_vars)
      if (.not. defined(gf%silamVertical)) call set_vertical(fu_level(field_id), gf%silamVertical)
!      call report(gf%silamVertical)
      gf%time_nbr = gf%time_nbr + 1 
!      gf%n_times_bin = gf%n_times_bin + 1  
!      call msg('Keep time period3:',gf%n_times_bin)
!      call msg('New time period3 total time:',gf%time_nbr)
      gf%var_nbr = 1
      gf%lev_nbr = 1
      if(present(forced_valid_time))then
        gf%gtime%step = forced_valid_time - gf%gtime%start
      else
        gf%gtime%step = fu_valid_time(field_id) - gf%gtime%start
      endif
      gf%gtime%ifVaryingStep = .false.   ! actually, do not know but better to initialise
      gf%gvars(1:gf%n_vars)%defined = silja_true
      gf%glevs%defined = silja_true
      gf%gtime%defined = silja_true
!      do iTmp=1,gf%n_vars
!        print *, fu_true(gf%gvars(iTmp)%defined)
!      end do
!      print *, (fu_true(gf%gvars(iTmp)%defined),iTmp=1,gf%n_vars)
      !
      ! Take this moment to write down the ctl file: for further references if the run crashes
      !
      if (.not. gf%ifMPIIO .or. smpi_adv_rank==0)then
        CALL write_ctl_file(gf, '')
      end if

    END IF ! IF new time period

!    call msg('Variable: grads Nbr & quantity:',gf%var_nbr, real(gf%gvars(gf%var_nbr)%quantity))
!    call msg('End fill_structures')
!    call msg('')

  END SUBROUTINE fill_structures


  ! ****************************************************************


  REAL FUNCTION fu_get_glevel(lev, levTypeGlobal)
    !
    ! Extracts the level value depending on the type of level. Handles
    ! exotic levels like surface-level.
    ! ATTENTION. Unit of levels depends on the type of the level. Variable
    ! levTypeGlobal helps to avoid mixture of units.
    !
    ! Code owner: M.Sofiev
    !
    IMPLICIT NONE
    !
    ! Imported parameters with intent IN
    TYPE(silja_level), INTENT(in) :: lev
    INTEGER, INTENT (inout) :: levTypeGlobal

    ! Local variable
    LOGICAL :: ifCheck ! If leveltype of lev is meaningful (any of 3D-type) => true

    SELECT CASE(fu_leveltype(lev))
    CASE(surface)
      fu_get_glevel = 0.0 ! Whatever, e.g. metres
      ifCheck = .false.
    CASE(constant_height)
      fu_get_glevel = fu_level_height(lev)  ! Metres
      ifCheck = .true.
    CASE(constant_altitude)
      fu_get_glevel = fu_level_altitude(lev) ! Metres
      ifCheck = .true.
    CASE(constant_pressure)
      fu_get_glevel = fu_pr_level_pressure(lev) ! Pa
      ifCheck = .true.
    CASE(layer_btw_2_height, layer_btw_2_altitude, layer_btw_2_pressure)
      fu_get_glevel = fu_layer_centre_value(lev) ! Pa
      ifCheck = .true.
    CASE(hybrid, layer_btw_2_hybrid)
      !fu_get_glevel = fu_hybrid_level_number(lev) ! UNITLESS
      fu_get_glevel = fu_height_for_press(fu_hybrid_level_pressure(lev, std_pressure_sl))
      ifCheck = .true.
    CASE DEFAULT
      if(fu_leveltype(lev) == entire_atmosphere_single_layer)then
        fu_get_glevel = real_missing
        ifCheck = .false.
      else
        CALL set_error('Can not handle level','fu_get_glevel')
        RETURN
      endif
    END SELECT

    IF(ifCheck)THEN !----------------------- Meaningful 3D level in lev
      IF(levTypeGlobal == int_missing)THEN ! Nothing is stored so far
        levTypeGlobal = fu_glevel_type(fu_leveltype(lev))
      ELSE  !----------- Something is stored in levTypeGlobal => make checking
        IF(levTypeGlobal /= fu_glevel_type(fu_leveltype(lev)))THEN
          CALL set_error('mixture of 3D level types','fu_get_glevel')
          fu_get_glevel = real_missing
          RETURN
        END IF
      END IF
    END IF
  END FUNCTION fu_get_glevel


  !**********************************************************************

  function fu_silamVert_of_grads(igf) result(vert)
    type(silam_vertical), pointer :: vert
    integer, intent(in) :: igf

    vert => gfile(igf)%ptr%silamVertical

  end function fu_silamVert_of_grads


  !**********************************************************************

  function fu_silamGrid_of_grads(igf) result(grid)
    type(silja_grid) :: grid
    integer, intent(in) :: igf

    grid = gfile(igf)%ptr%silamGrid

  end function fu_silamGrid_of_grads

  !**********************************************************************

  function fu_time_of_grads(igf, iTime) result(time)
    type(silja_time) :: time
    integer, intent(in) :: igf, iTime

    if(gfile(igf)%ptr%gtime%ifVaryingStep)then
      time = gfile(igf)%ptr%gtime%arTimes(iTime)
    else
      time = gfile(igf)%ptr%gtime%start + gfile(igf)%ptr%gtime%step * (iTime-1)
    endif

  end function fu_time_of_grads

  !**********************************************************************

  subroutine get_grads_times(igf, times, nTimes, ifEnvelope)
    !
    ! This subroutine returns the ENVELOPE of grads time period:
    ! - the first time is when the first field validity starts
    ! - the last time is when the validity of the last field ends
    ! - total number of times returned is n_times+1, of course
    !
    implicit none

    ! Imported parameters
    type(silja_time), dimension(:), pointer :: times
    integer, intent(in) :: igf
    logical, intent(in) :: ifEnvelope
    integer, intent(out) :: nTimes

    ! Local variables
    integer :: iTmp, iShift
    real :: fShift

    if(ifEnvelope)then
      allocate(times(gfile(igf)%ptr%n_times+1),stat=iTmp)
    else
      allocate(times(gfile(igf)%ptr%n_times),stat=iTmp)
    endif
    if(fu_fails(iTmp==0, 'Failed times array allocation,size=' + fu_str(gfile(igf)%ptr%n_times), &
                       & 'get_grads_times'))return

    if(gfile(igf)%ptr%gtime%ifVaryingStep)then
      !
      ! Varying time step. Note that the arTimes is from 0 to n_times+1. Use it!
      !
      if(ifEnvelope)then
        select case (gfile(igf)%ptr%time_label_position)
          case(start_of_period)
            iShift = 0
          case(mid_of_period)
            call set_error('Mid-period time label is not supported for varying time step', &
                         & 'get_grads_times')
          case(end_of_period)
            iShift = 1
          case default
            call set_error('Unknown time_label_position:'+fu_str(gfile(igf)%ptr%time_label_position), &
                         & 'get_grads_times')
        end select
        do iTmp = 1, gfile(igf)%ptr%n_times+1
          times(iTmp) = gfile(igf)%ptr%gtime%arTimes(iTmp-iShift)  ! get times with possile one-value shift
        end do
        nTimes = gfile(igf)%ptr%n_times + 1
      else
        do iTmp = 1, gfile(igf)%ptr%n_times
          times(iTmp) = gfile(igf)%ptr%gtime%arTimes(iTmp)  ! get times with one-value shift
        end do
        nTimes = gfile(igf)%ptr%n_times
      endif
    else
      !
      ! Fixed time step.
      !
      if(ifEnvelope)then
        select case (gfile(igf)%ptr%time_label_position)
          case(start_of_period)
            fShift = 1.0
          case(mid_of_period)
            fShift = 1.5
          case(end_of_period)
            fShift = 2.0
          case default
            call set_error('Unknown time_label_position:'+fu_str(gfile(igf)%ptr%time_label_position), &
                         & 'get_grads_times')
            return
        end select
        do iTmp = 1, gfile(igf)%ptr%n_times+1
          times(iTmp) = gfile(igf)%ptr%gtime%start + gfile(igf)%ptr%gtime%step * (real(iTmp)-fShift)
        end do
        nTimes = gfile(igf)%ptr%n_times+1
      else
        do iTmp = 1, gfile(igf)%ptr%n_times
          times(iTmp) = gfile(igf)%ptr%gtime%start + gfile(igf)%ptr%gtime%step * (real(iTmp)-1)
        end do
        nTimes = gfile(igf)%ptr%n_times
      endif  ! ifEnvelope

    endif  ! if varying time step

  end subroutine get_grads_times

!  !**********************************************************************
!
!  function fu_time_step_of_grads(igf) result(step)
!    type(silja_interval) :: step
!    integer, intent(in) :: igf
!
!    step = gfile(igf)%ptr%gtime%step
!
!  end function fu_time_step_of_grads


  !**********************************************************************

  function fu_validity_length_from_grads(igf) result(duration)
    type(silja_interval) :: duration
    integer, intent(in) :: igf

    duration = gfile(igf)%ptr%gtime%validity_length

  end function fu_validity_length_from_grads

  !**********************************************************************

  integer function fu_data_time_features_grads(igf)
    integer, intent(in) :: igf

    fu_data_time_features_grads = gfile(igf)%ptr%data_time_features

  end function fu_data_time_features_grads

  !**********************************************************************

  integer function fu_grads_time_index(igFile, timeMoment, direction, ifAcceptSameMonth_) result(ind)
    !
    ! Computes the actual index of the time moment in the given GrADS structure.
    ! Note that in case of non-integer time, the given interval must cover the required time
    ! rather than just be the closest possible one.
    ! However, exact hit into a grads slot leads to this very index returned.
    !
    implicit none

    ! Imported parameters
    integer, intent(in) :: igFile, direction
    type(silja_time), intent(in) :: timeMoment
    logical, intent(in) :: ifAcceptSameMonth_

    ! Local variables
    logical :: ifAmbiguity, ifAcceptSameMonth
    integer :: iTmp, iShift, iEnvelope, iMonthNeeded
    real :: fIndex
    !
    ! First of all, if we are dealing with static fields, nothing depends on the input time
    !
    if(gfile(igFile)%ptr%data_time_features == static_climatology)then
      ind = 1
      return
    endif
    !
    ! In case of monthly climatology, reduced dependence on time too:
    !
    ifAcceptSameMonth = ifAcceptSameMonth_ .or. &
                      & (gfile(igFile)%ptr%data_time_features == monthly_climatology)

    if(fu_fails(defined(timeMoment),'Undefined time given','fu_grads_time_index'))return
    !
    ! The shift depends on position of the grads time stamp: the first time is actually in grads
    ! only if the stamp is at the start of validity period - or the file contains instsnt fields
    !
    select case (gfile(igfILE)%ptr%time_label_position)
      case(start_of_period)
        iShift = 0
        iEnvelope = 1
      case(mid_of_period)
        call set_error('Mid-period time label is not supported for varying time step', &
                     & 'fu_grads_time_index')
      case(end_of_period)
        iShift = 1
        iEnvelope = 1
      case(instant_fields)
        iShift = 0
        iEnvelope = 0
      case default
        call set_error('Unknown time_label_position:'+fu_str(gfile(igfILE)%ptr%time_label_position), &
                     & 'fu_grads_time_index')
    end select
      
    iMonthNeeded = fu_mon(timeMoment)  ! might need it
    
    !
    ! Now, the procedure depends on whether the file has fixed time step
    !
    if(gfile(igFile)%ptr%gtime%ifVaryingStep)then
      !
      ! Find the first time smaller than the given one. 
      ! Return the previous index, whihc will be either exact hit or within the validity interval 
      ! from that index
      !
      do iTmp = 1, gfile(igfILE)%ptr%n_times + iEnvelope
        
        if(ifAcceptSameMonth)then
          if(fu_mon(gfile(igfile)%ptr%gtime%arTimes(iTmp) - gfile(igfile)%ptr%gtime%step*iShift) == iMonthNeeded)then
            ind = iTmp
            exit
          endif
        else
          if (gfile(igfile)%ptr%gtime%arTimes(iTmp) - gfile(igfile)%ptr%gtime%step*iShift + one_second > timeMoment) then
            ind = iTmp
            exit
          endif
        endif  ! if same month
      end do
      
    else
      if(ifAcceptSameMonth)then
        !
        ! Have to scan until hit the right month
        !
        do iTmp = 1, gfile(igfILE)%ptr%n_times + iEnvelope
          if(fu_mon(gfile(igFile)%ptr%gtime%start + (gfile(igFile)%ptr%gtime%step * real(iTmp-1))) == iMonthNeeded)then
            ind = iTmp
            exit
          endif
        end do
        
      else
        !
        ! For regular files, first, dumb search of the index taking care of not hitting the edge
        !
        fIndex = (timeMoment - gfile(igFile)%ptr%gtime%start) / gfile(igFile)%ptr%gtime%step


        if(abs(fIndex - real(nint(fIndex))) < 0.1)then
          !
          ! If exact hit, just return the index
          !
          !!! Adding shift here breaks restore form dump
          ind = nint(fIndex) + 1 ! + iShift

          ! Dirty hack to still enable reading the first time step on exact hit and end_of_period tag
          !
          if(ind == iShift .and. gfile(igFile)%ptr%time_label_position == end_of_period)then
            ind = 1 
            call set_error('grads index hack triggered','fu_grads_time_index')
            call unset_error('fu_grads_time_index')
          endif
        else
          !
          ! No hit into grads time stamp, need to consider validity period
          !
          select case(gfile(igFile)%ptr%time_label_position)
            case (start_of_period)
              ind = int(fIndex - 0.01) + 1
            case (mid_of_period)
              ind = int(fIndex + 0.49) + 1
            case (end_of_period)
              ind = int(fIndex + 0.99) + 1
            case default
              call set_error('Strange time_label_position','fu_grads_time_index')
          end select
        endif  ! if regular file
      endif  ! if accept same month
      !!! Checking the ambiguity: are we at the edge? We can then use either of the neighbouring
      !!! intervals. Let's do it this way: 
      !!! - if only one of the intervals is present, give it
      !!! - if both intervals present, give the one that does not overlap with the next model time step
      !!!
      !!select case(gfile(igFile)%ptr%time_label_position)
      !!  case (start_of_period)
      !!    ifHit = (timeMoment == (gfile(igFile)%ptr%gtime%start + &
      !!                                & gfile(igFile)%ptr%gtime%step * real(ind)))
      !!  case (mid_of_period)
      !!    ifHit = (timeMoment == (gfile(igFile)%ptr%gtime%start + &
      !!                                & gfile(igFile)%ptr%gtime%step * (real(ind) - 0.5)))
      !!  case (end_of_period)
      !!    ifHit = (timeMoment == (gfile(igFile)%ptr%gtime%start + &
      !!                   a            & gfile(igFile)%ptr%gtime%step * (real(ind) - 1.)))
      !!  case default
      !!    call set_error('Strange time_label_position','fu_grads_time_index')
      !!end select
      !!!
      !!! If we have two intervals covering the required time, i.e. it is at their edges,
      !!! we shall use the direction information to choose the right index, if it exists, of course.
      !!! If only one of the indices exists, take it - formally it is valid.
      !!!
      !!if(ifAmbiguity)then
      !!    if(direction == backwards)then
      !!      if(ind + 1 <= gfile(igFile)%ptr%n_times) ind = ind + 1  ! usable
      !!    elseif(direction == forwards)then
      !!      if(ind < 1) ind = ind + 1 ! no coverage with default index => +1
      !!    else
      !!      call set_error('Unknown direction:'+fu_str(direction),'fu_grads_time_index')
      !!    endif
      !!endif  ! if ambiguity in time coverage
    
  !    call msg('time moment:' + fu_str(timeMoment) + &
  !           & ', grads start:' + fu_str(gfile(igFile)%ptr%gtime%start) + &
  !           & ', time step [hr] real index:', fu_hour(gfile(igFile)%ptr%gtime%step), &
  !           & ((timeMoment - gfile(igFile)%ptr%gtime%start) / &
  !                                 & gfile(igFile)%ptr%gtime%step) + 1.499)
  !    call msg('Resulting time & index:' + &
  !            & fu_str(gfile(igFile)%ptr%gtime%start + &
  !                              & (gfile(igFile)%ptr%gtime%step * real(fu_grads_time_index-1))), &
  !            & fu_grads_time_index)
  !    call msg('Residual resulting_time - target[hr]:', &
  !           & fu_hour(gfile(igFile)%ptr%gtime%start + &
  !                   & (gfile(igFile)%ptr%gtime%step * real(fu_grads_time_index-1)) - timeMoment))

    endif   ! if ifVaryingStep
    
    if(error .or. ind < 1 .or. ind > gfile(igFile)%ptr%n_times) ind = int_missing

  end function fu_grads_time_index


  !******************************************************************

  integer function fu_n_gvars(igFile)
    implicit none
    integer, intent(in) :: igFile
    if(gfile(igFile)%ptr%defined == silja_true)then
      fu_n_gvars = gfile(igFile)%ptr%n_vars
    else
      fu_n_gvars = int_missing
    endif
  end function fu_n_gvars


  !******************************************************************

  integer function fu_n_glevs(igFile)
    implicit none
    integer, intent(in) :: igFile
    if(gfile(igFile)%ptr%defined == silja_true)then
      fu_n_glevs = gfile(igFile)%ptr%n_levs
    else
      fu_n_glevs = int_missing
    endif
  end function fu_n_glevs
  !******************************************************************

  integer function fu_n_gVar_levs(igFile, iVar)
    implicit none
    integer, intent(in) :: igFile, iVar
    if (iVar < 1 .or. iVar > gfile(igFile)%ptr%n_vars)then
      call set_error('Strange var index','fu_n_gVar_levs')
      return
    endif
    if(gfile(igFile)%ptr%defined == silja_true)then
      if (gfile(igFile)%ptr%gvars(iVar)%n_levs == 0) then
        fu_n_gVar_levs = 1
      else
        fu_n_gVar_levs = gfile(igFile)%ptr%gvars(iVar)%n_levs
      endif
    else
      fu_n_gVar_levs = int_missing
    endif
  end function fu_n_gVar_levs

  !******************************************************************

  integer function fu_n_gtimes(igFile)
    implicit none
    integer, intent(in) :: igFile
    if(gfile(igFile)%ptr%defined == silja_true)then
      fu_n_gtimes = gfile(igFile)%ptr%n_times
    else
      fu_n_gtimes = int_missing
    endif
  end function fu_n_gtimes


  !********************************************************************

  subroutine get_grads_var_metadata(igFile, iVar, iLev, indTime, id) ! indices: gfile, gvar, glev, gtime; SILAM-id
    !
    ! Returns the metadata of the GrADS variable in a form of SILAM field id
    !
    implicit none

    ! Imported parameters
    integer, intent(in) :: igFile, iVar, iLev, indTime ! indices of gfile, gvar,, glev, gtime; 
    type(silja_field_id), intent(out) :: id            ! SILAM-id
    type(silja_interval) :: fcast_len, validity_len

    ! Local variables
    type(silja_level) :: level

    call set_missing(id)

    if(.not. (gfile(igFile)%ptr%defined == silja_true)) return

    if(gfile(igFile)%ptr%gtime%ifVaryingStep)then
      select case(gfile(igFile)%ptr%time_label_position)
        case(start_of_period)
          fcast_len = gfile(igFile)%ptr%gtime%arTimes(indTime) - gfile(igFile)%ptr%gtime%start
        case(mid_of_period)
          fcast_len = ((gfile(igFile)%ptr%gtime%arTimes(indTime-1) - &
                                                          & gfile(igFile)%ptr%gtime%start) + &
                    &  (gfile(igFile)%ptr%gtime%arTimes(indTime) - &
                                                          & gfile(igFile)%ptr%gtime%start)) / 2.0
        case(end_of_period)
          fcast_len = (gfile(igFile)%ptr%gtime%arTimes(indTime-1) - gfile(igFile)%ptr%gtime%start)
        case default
          call set_error('strange time_label_position','get_grads_var_metadata')
          return
      end select
    else
      select case(gfile(igFile)%ptr%time_label_position)
        case(start_of_period, instant_fields)
          fcast_len = gfile(igFile)%ptr%gtime%step * (indTime - 1)
        case(mid_of_period)
          fcast_len = gfile(igFile)%ptr%gtime%step * (real(indTime) - 1.5)
        case(end_of_period)
          fcast_len = gfile(igFile)%ptr%gtime%step * (indTime - 2) !+ one_minute
        case default
          call set_error('strange time_label_position','get_grads_var_metadata')
          return
      end select
    endif
    
    select case(gfile(igFile)%ptr%data_time_features)
      case(dynamic_map)
        validity_len = gfile(igFile)%ptr%gtime%step
      case(monthly_climatology)
        validity_len = interval_missing
      case(static_climatology)
        validity_len = very_long_interval
      case default
        call set_error('Unknown data_time_features:' + &
                     & fu_str(gfile(igFile)%ptr%data_time_features),'get_grads_var_metadata')
        return
    end select
    !
    ! SilamLevel, if defined, has priority over iVerticalFeature
    !
    select case(gfile(igFile)%ptr%gvars(iVar)%iVerticalFeature)
      case(level_2d_type_flag)
        if(defined(gfile(igFile)%ptr%gvars(iVar)%SilamLevel))then
          level = gfile(igFile)%ptr%gvars(iVar)%SilamLevel
        else
          level = surface_level
        endif
      case(integrate_column_flag)
        if(defined(gfile(igFile)%ptr%gvars(iVar)%SilamLevel))then
          level = gfile(igFile)%ptr%gvars(iVar)%SilamLevel
        else
          level = entire_atmosphere_integr_level
        endif
      case(lowest_level_flag)
        ! level = fu_level(gfile(igFile)%ptr%silamVertical, 1)
        if(defined(gfile(igFile)%ptr%gvars(iVar)%SilamLevel))then
          level = gfile(igFile)%ptr%gvars(iVar)%SilamLevel
        else
          level = lowest_atmosphere_level
        endif
      case default
        level = fu_level(gfile(igFile)%ptr%silamVertical, iLev)  ! 3d variable
    end select
    id = fu_set_field_id(met_src_missing,&
                       & gfile(igFile)%ptr%gvars(iVar)%quantity , &
                       & gfile(igFile)%ptr%gtime%start, &             ! analysis_time,&
                       & fcast_len, &                                 ! forecast_length, &
                       & gfile(igFile)%ptr%silamGrid,&                ! grid
                       & level, &                             ! level
                       & interval_missing, &                    ! length_of_accumulation, &
                       & validity_len, &                        ! length_of_validity, &
                       & forecast_flag, &                       ! field_kind, &
                       & species = gfile(igFile)%ptr%gvars(iVar)%species, &   ! species
                       & chCocktail = gfile(igFile)%ptr%gvars(iVar)%chCocktailNm) ! of cocktail
    if(error)call set_missing(id)

  end subroutine get_grads_var_metadata


  !*****************************************************************
  
  subroutine get_grads_IDs(igf, idList, nIDs)
    !
    ! Creates a full ID list for gards file. All IDs have first time as the analysis time
    ! and zero forecast length.
    ! Algorithm is based on numerous calls of get_grads_var_metadata
    !
    implicit none
    
    ! Imported parameters
    integer, intent(in) :: igf
    type(silja_field_id), dimension(:), pointer :: idList
    integer, intent(out) :: nIDs
    
    ! Local variables
    integer :: iVar, iLev, iID
    
    !
    ! First, count the number of independent IDs
    !
    nIDs = 0
    do iVar = 1, gfile(igf)%ptr%n_vars
      if(gfile(igf)%ptr%gvars(iVar)%iVerticalFeature == level_3d_type_flag)then
        nIDs = nIDs + gfile(igf)%ptr%n_levs
      else
        nIDs = nIDs + 1
      endif
    end do
    allocate(idList(nIDs), stat = iVar)
    if(fu_fails(iVar==0,'Failed allocation of ID list, size='+fu_str(nIDs),'get_gards_IDs'))return

    iID = 1
    do iVar = 1, gfile(igf)%ptr%n_vars
      if(gfile(igf)%ptr%gvars(iVar)%iVerticalFeature == level_3d_type_flag)then
        do iLev = 1, gfile(igf)%ptr%n_levs
          call get_grads_var_metadata(igf, iVar, iLev, 1, idList(iID))
          iID = iID + 1
        end do
      else
        call get_grads_var_metadata(igf, iVar, 1, 1, idList(iID))
        iID = iID + 1
      endif
    end do

  end subroutine get_grads_IDs


  !******************************************************************
  
  subroutine get_grads_total(ugf, flds_3D_Requested, start, duration, layer, pOutput)
    !
    ! Searches for the given quantity and species, if defined, limits time period and 
    ! gets the grand totals for all species (or just one, if it is defined)
    !
    implicit none
    
    ! Imported parameters
    integer, intent(in) :: ugf
    type(silja_field_3d_pointer), dimension(:), pointer :: flds_3D_Requested
    type(silja_time), intent(in) :: start
    type(silja_interval), intent(in) :: duration
    type(silja_level), intent(in) :: layer
    type(silja_rp_1d), dimension(:), pointer :: pOutput
    
    ! Local parameters
    type(grads_file), pointer :: gf
    integer :: iVar, iTime, iFldRequested, iFldQuantity, iLev, fs, field_kind, nTimes, iTmp
    type(silam_species) :: FldSpecies
    character(len=substNmLen) :: chFldCocktailNm
    real, dimension(:), pointer :: pData
    type(silja_field_id) :: idTmp
    type(silja_interval) :: length_of_accumulation
    real :: duration_sec, fVertFract
    type(silja_time) :: startTmp, endTmp
    type(silja_time), dimension(:), pointer :: envelope_times
    !
    ! Some preparations first
    !
    if(.not. associated(gfile))then
      call set_error('gf not associated','get_grads_total')
      return
    endif
    if(.not. associated(gfile(ugf)%ptr))then
      call set_error('gf ptr not associated','get_grads_total')
      return
    endif
    gf => gfile(ugf)%ptr
    do iVar = 1, size(flds_3d_Requested)
      pOutput(iVar)%pp(1 : fu_number_of_gridpoints(fu_grid( & 
                              & fu_id(fu_field_from_3d_field(flds_3D_Requested(iVar)%fp,1))))) = 0.0
    end do
    call get_grads_times(ugf, envelope_times, nTimes, .true.) ! ifEnvelope_
    if(error)return
    !
    ! Store the required 
    !
    do iTime = 1, nTimes-1
      !
      ! Scan the list of variables, set the required ID and read the stuff
      !
      if(envelope_times(iTime+1) <= start)then
        cycle    ! too early
      elseif(envelope_times(iTime) >= start + duration)then
        return  ! all done
      endif
      if(envelope_times(iTime) < start)then
        startTmp = start
      else
        startTmp = envelope_times(iTime)
      endif
      if(envelope_times(iTime+1) > start + duration)then
        endTmp = start + duration
      else
        endTmp = envelope_times(iTime+1)
      endif
      duration_sec = fu_sec(endTmp - startTmp)

      do iFldRequested = 1, size(flds_3d_Requested)
        idTmp = fu_id(fu_field_from_3d_field(flds_3D_Requested(iFldRequested)%fp,1))
        pData => fu_grid_data(fu_field_from_3d_field(flds_3D_Requested(iFldRequested)%fp,1))
        iFldQuantity = fu_quantity(idTmp)
        FldSpecies = fu_species(idTmp)
        chFldCocktailNm = fu_cocktail_name(idTmp)
        fs = fu_number_of_gridpoints(fu_grid(idTmp))

        do iVar = 1, gf%n_vars
          !
          ! Needed variable? It is defined by quantity and species / cocktail_name
          !
          if(gf%gVars(iVar)%quantity == int_missing)cycle
          if(iFldQuantity /= int_missing)then
            if(iFldQuantity /= gf%gVars(iVar)%quantity)cycle
          endif
          if(defined(FldSpecies))then
            if(.not. (FldSpecies == gf%gVars(iVar)%species))cycle
          endif
          if(len_trim(chFldCocktailNm) > 0)then
            if(.not. (chFldCocktailNm == gf%gVars(iVar)%chCocktailNm))cycle
          endif
          !
          ! OK, variable found that meets the requirements of this field_3d
          ! Go layer by layer summing-up the maps
          !
          if(fu_accumulated_quantity(gf%gVars(iVar)%quantity))then
            field_kind = accumulated_flag
            length_of_accumulation = envelope_times(iTime) - envelope_times(1)
          else
            field_kind = forecast_flag
            length_of_accumulation = zero_interval 
          endif

          do iLev = 1, gf%gVars(iVar)%n_levs
            fVertFract =  fu_vert_overlap_fraction(fu_level(gf%silamVertical, iLev), layer)
            if(error)return
            if(fVertFract < 1e-5)cycle

            idTmp = fu_set_field_id(met_src_missing, &
                                  & gf%gVars(iVar)%quantity, &
                                  & envelope_times(1),&
                                  & envelope_times(iTime) - envelope_times(1), &
                                  & gf%silamGrid,&
                                  & fu_level(gf%silamVertical, iLev),&
                                  & length_of_accumulation, & ! optional
                                  & zero_interval, &     !optional
                                  & field_kind, &             ! optional
                                  & species = gf%gVars(ivar)%species, &
                                  & chCocktail = gf%gVars(iVar)%chCocktailNm)
            if(error)return

            call read_field_from_grads_id(ugf, idTmp, pData, 0.0)
            if(error)return

            do iTmp = 1, fs
              pOutput(iFldRequested)%pp(iTmp) = pOutput(iFldRequested)%pp(iTmp) + &
                                              & pData(iTmp) * duration_sec * fVertFract
            end do
          enddo   ! iLev
        enddo  ! iVar
      enddo  ! iFldRequested
    enddo  ! nf%n_times

  end subroutine get_grads_total

  !*****************************************************************

  integer function fu_glevel_type(leveltype)
    !
    ! Translates SILAM level types to the GrADS level types.
    ! Actually, the only purpose of the routine is to translate
    ! the thick layer types to their central points
    !
    implicit none

    ! Imported parameter
    integer, intent(in) :: leveltype

    select case(leveltype)
      case(surface, &
         & top_of_the_atmosphere, &
         & constant_pressure, & 
         & mean_sea, &
         & constant_altitude, &
         & constant_height, &
         & sigma_level, &
         & hybrid, &
         & depth_level)
        fu_glevel_type = leveltype ! single levels do not change

      case(layer_btw_2_pressure)
        fu_glevel_type = constant_pressure

      case(layer_btw_2_altitude)
        fu_glevel_type = constant_altitude

      case(layer_btw_2_height)
        fu_glevel_type = constant_height

      case(layer_btw_2_depth)
        fu_glevel_type = depth_level

      case(layer_btw_2_sigma)
        fu_glevel_type = sigma_level

      case(layer_btw_2_hybrid)
        fu_glevel_type = hybrid

    case default
      call msg('Strange leveltype',leveltype)
      call set_error('Non-supproted leveltype','fu_glevel_type')
      fu_glevel_type = int_missing
    end select

  end function fu_glevel_type

  ! ****************************************************************


  SUBROUTINE close_gradsfile_i(igf)
    !
    ! Closes the GrADS binary file and then cleans up the structure
    !
    ! Author: M.Sofiev
    !
    IMPLICIT NONE
    !
    ! Imported parameter with intent IN
    INTEGER, INTENT(in) :: igf ! counter of the GrADS file to close
    TYPE(grads_file), POINTER :: gf

    if(igf == int_missing .or. igf < 1)then
      call set_error('Undefined grads file index','close_grads_file_i')
      return
    endif

    gf => gfile(igf)%ptr

    !!! Free reads buffer
    if (gf%gradsbuf%allocated) call free_grads_buffer(gf%gradsbuf)

    if (gf%unit_bin > 0) close(gf%unit_bin) ! Close the binary (if any)

    CALL release_index(igf) ! Free-up the structure.

  END SUBROUTINE close_gradsfile_i


  ! ****************************************************************


  SUBROUTINE close_gradsfile_o(igf, chFixedNameTemplate)
    !
    ! Closes the GrADS binary file, writes its ctl file using the gfile data,
    ! and then cleans up the structure
    !
    ! Author: M.Sofiev
    !
    IMPLICIT NONE
    !
    ! Imported parameter with intent IN
    INTEGER, INTENT(in) :: igf ! counter of the GrADS file to close
    character(len=*), intent(in) :: chFixedNameTemplate
!    TYPE(grads_file), POINTER :: gfileptr

    if(igf == int_missing .or. igf < 1)then
      call set_error('Undefined grads file index','close_grads_file_o')
      return
    endif
!    gfileptr => gfile(igf)%ptr

    if (gfile(igf)%ptr%ifBuffered) call flush_buffer(igf)
    if(gfile(igf)%ptr%ifMPIIO)then
      call smpi_close_gradsfile_mpiio(gfile(igf)%ptr%unit_bin)
    else
      close(gfile(igf)%ptr%unit_bin) ! Close the binary
    end if

    ! Only single process should write the control file for MPIIO
    if(.not.gfile(igf)%ptr%ifMPIIO.or.smpi_adv_rank==0)then
      !
      ! Write two ctl/super_ctl file sets: one for the just-closed binary, one for the whole series
      ! Note that the binary-specific ctl will cover only this binary time range
      ! The fname_initial ctl file (actually, the first file name used to open this igf) will have 
      ! the whole length
      !
      CALL write_ctl_file(gfile(igf)%ptr, chFixedNameTemplate, .true.)
      
      gfile(igf)%ptr%fname = gfile(igf)%ptr%fname_initial
      
      CALL write_ctl_file(gfile(igf)%ptr, chFixedNameTemplate)
    end if

    CALL release_index(igf) ! Free-up the structure.

  END SUBROUTINE close_gradsfile_o


  !*******************************************************************

  subroutine invert_grads_binary(igf, iInvert)
    !
    ! Used for time invesion of the data fields in the GrADS binary
    ! So far handles only single-binary GrADS file.
    ! Suitable for inverse runs when timestep is negative, which is 
    ! forbidden in GrADS. As a result of this subroutine, the last
    ! time in the binary becomes the first one.
    ! 
    implicit none

    ! Imported parameters
    integer, intent(in) :: igf, iInvert

    ! Local variables
    TYPE(grads_file), POINTER :: gf
    real, dimension(:), pointer :: dataPtr
    integer :: iSize, iRecRead, iRecWrite, i, iVar, iLev, indTime, uFTmp, nRecsPerTime
    character(len=clen) :: chTmp
      
    gf => gfile(igf)%ptr
    !
    ! Stupidity check
    !
    if(gf%unit_bin < 0)then
      call set_error('Undefined GrADS binary','invert_grads_binary')
      return
    endif
    !
    ! Internal preparations
    !
    iSize = gf%ggrid%nx * gf%ggrid%ny

    uFTmp = fu_next_free_unit()
    if(error)return

    write(unit=chTmp,fmt='(A10,I0,A4)') 'grads_tmp_', fu_pid(), '.bin'
    call open_grads_binary_o(chTmp, uFTmp, iSize)    
    !call open_gradsfile_o_md('grads.bin.tmp', uFTmp, iSize)

    nRecsPerTime = 0
    do iVar = 1, gf%n_vars
      nRecsPerTime = nRecsPerTime + max(gf%gvars(iVar)%n_levs, 1)
    end do

    call msg (fu_connect_strings('Inverting GrADS binary: ', gf%fname))
    call msg('Number of timesteps in binary file: ', gf%n_times_bin)
    call msg('Number of records per time: ', nRecsPerTime)

    dataPtr => fu_work_array()
    if(error)return

    !
    ! Preliminary cycle - copying the main file to the temporary one
    !     
    do iRecRead = 1, nRecsPerTime * gf%n_times_bin
      read(gf%unit_bin,rec=iRecRead)(dataPtr(i),i=1,iSize)
      write(uFTmp, rec=iRecRead)(dataPtr(i),i=1,iSize)
    end do

    !
    ! The main cycle inverting the time order but keeping variable order inside 
    ! one time
    !
    iRecWrite = 1
    do indTime = gf%n_times_bin, 1, -1
      iRecRead = nRecsPerTime * (indTime-1) + 1
      do iVar = 1, gf%n_vars
        do iLev = 1, max(gf%gvars(iVar)%n_levs,1)
          read(uFTmp,rec=iRecRead)(dataPtr(i),i=1,iSize)
          write(gf%unit_bin, rec=iRecWrite)(dataPtr(i),i=1,iSize)
          iRecWrite = iRecWrite + 1
          iRecRead = iRecRead + 1
        end do
      end do
    end do

    close(uFTmp, status='delete')

    !
    ! Now we have to adjust the time parameter in the file structure
    !
    if(iInvert == 2)gf%gtime%start = gf%gtime%start + (gf%gtime%step * (real(gf%n_times)-1.))
    gf%gtime%start_bin = gf%gtime%start_bin + (gf%gtime%step * (real(gf%n_times_bin)-1.))

    call free_work_array(dataPtr)

  end subroutine invert_grads_binary


  ! ****************************************************************


  SUBROUTINE write_ctl_file(gf, chFixedNameTemplate, if_current_binary_)
    !
    ! Just writes the ctl file having the complete information in the grads_file 
    ! structure.
    ! After that writes the super-ctl file for further reading by SILAM itself
    ! GrADS output can require a duplicate ctl file to be created - e.g. "latest.ctl",
    ! this is handled via chFixedNameTemplate.
    !
    USE silam_partitioning
    IMPLICIT NONE
    !
    ! Parameters imported with intent IN
    TYPE(grads_file), INTENT(inout) :: gf
    character(len=*), intent(in) :: chFixedNameTemplate
    logical, intent(in), optional :: if_current_binary_

    ! Local stuff
    INTEGER :: iTmp, iCtlUnit, globalNX, globalNY
    REAL :: globalStartX, globalStartY, stepX, stepY
    type(silam_sp) :: sp
    CHARACTER(LEN = fnlen) :: chTmp
    logical :: if_current_binary
    REAL :: southpole_lon_E, southpole_lat_N
    LOGICAL :: corner_in_geo_lonlat

    if(present(if_current_binary_))then
      if_current_binary = if_current_binary_
    else
      if_current_binary = .false.
    endif
    
    !------------------------------------------------------------------
    !
    ! Start the GrADS ctl file
    !
    iCtlUnit = fu_next_free_unit()
    if(error)then
      call set_error('Failed to get unit for ctl file','write_ctl_file')
      return
    endif
    open(iCtlUnit, file = trim(gf%fname +'.ctl'))

!    if(.not. gf%defined == silja_true)then
!      call msg_warning('Fully undefined ctl file','write_ctl_file')
!      write(iCtlUnit,*)'FULLY UNDEFINED FILE'
!      return
!    endif

    sp%sp => fu_work_string()

!Bad thing with absolute path
    !WRITE(iCtlUnit,'(A5,A)') 'DSET ',gf%chTemplate(1:len_trim(gf%chTemplate)) 
!This is better...
    chTmp = fu_trim_grads_hat(gf%chTemplate,gf%fname)
    call replace_string(chtmp, '%task', trim(fu_str(smpi_adv_rank)))
    WRITE(iCtlUnit,'(A6,A)') 'DSET  ',trim(chTmp)
    WRITE(iCtlUnit,'(A24)')  'TITLE SILAM GrADS output'
    ! Binaray output no longer converted to big endian
    WRITE(iCtlUnit,'(A21)')  'OPTIONS LITTLE_ENDIAN'
    if(index(gf%chTemplate, '%') > 0) write(iCtlUnit,'(A16)') 'OPTIONS TEMPLATE'
    WRITE(iCtlUnit,'(A6,1X,E15.6)') 'UNDEF ', gf%missing_value
    
    if(gf%ifMPIIO)then ! Write values of global grid to the ctl file
    !  call smpi_get_global_dispgrid_dims(globalStartX, globalNX, stepX, &
    !            & globalStartY, globalNY, stepY)

        CALL lonlat_grid_parameters(wholeMPIdispersion_grid, globalStartX,  globalStartY, corner_in_geo_lonlat, &
         & globalNX, globalNY, southpole_lon_E, southpole_lat_N, stepX, stepY)


      WRITE(iCtlUnit,'(A5,I5,A8,2(F15.7,1x))') 'XDEF ',globalNX,' LINEAR ', &
                                                  & globalStartX, &
                                                  & stepX
      WRITE(iCtlUnit,'(A5,I5,A8,2(F15.7,1x))') 'YDEF ',globalNY,' LINEAR ', &
                                                  & globalStartY, &
                                                  & stepY
    else ! Otherwise write out the parameters from local grid
      WRITE(iCtlUnit,'(A5,I5,A8,2(F15.7,1x))') 'XDEF ',gf%ggrid%nx,' LINEAR ', &
                                                  & gf%ggrid%x_start, &
                                                  & gf%ggrid%x_step
      WRITE(iCtlUnit,'(A5,I5,A8,2(F15.7,1x))') 'YDEF ',gf%ggrid%ny,' LINEAR ', &
                                                  & gf%ggrid%y_start, &
                                                  & gf%ggrid%y_step
      !
      ! Attention with levels: may be too long string. Write compact
      !
    end if
    
    IF(gf%n_levs == 0)THEN
      WRITE(iCtlUnit,'(A15)') 'ZDEF 1 LEVELS 0'
    ELSE
      WRITE(unit=sp%sp,fmt='(I4,A7)') gf%n_levs,' LEVELS'
      do iTmp = 1, gf%n_levs
        if(abs(gf%glevs%levels(iTmp)) < 10)then
          write(unit=chTmp,fmt='(F8.5)') gf%glevs%levels(iTmp)
        elseif(abs(gf%glevs%levels(iTmp)) < 1e3)then
          write(unit=chTmp,fmt='(F8.3)') gf%glevs%levels(iTmp)
        elseif(abs(gf%glevs%levels(iTmp)) < 1e5)then
          write(unit=chTmp,fmt='(F8.1)') gf%glevs%levels(iTmp)
        else
          write(unit=chTmp,fmt='(F8.0)') gf%glevs%levels(iTmp)
        endif
        sp%sp = fu_connect_with_space(sp%sp, chTmp)
      end do
      WRITE(iCtlUnit,'(A5,A)') 'ZDEF ',trim(sp%sp)
    END IF
    if(gf%n_times == 1)then
      WRITE(iCtlUnit,'(A16,A20,A4)') 'TDEF  1  LINEAR ', &
                             & fu_time_to_grads_string(gf%gtime%start), ' 1hr'
    else
      if(if_current_binary)then
        WRITE(iCtlUnit,'(A5,I7,A8,A20,A10)') 'TDEF ',gf%n_times_bin,' LINEAR ', &
                               & fu_time_to_grads_string(gf%gtime%start_bin), &
                               & fu_interval_to_grads_string(gf%gtime%step)
      else
        WRITE(iCtlUnit,'(A5,I7,A8,A20,A10)') 'TDEF ',gf%n_times,' LINEAR ', &
                               & fu_time_to_grads_string(gf%gtime%start), &
                               & fu_interval_to_grads_string(gf%gtime%step)
      endif
    endif
    WRITE(iCtlUnit,'(A5,I5)')'VARS ',gf%n_vars


    DO iTmp =1, gf%n_vars
      if(defined(gf%gvars(iTmp)%species))then
        sp%sp = fu_quantity_short_string(gf%gvars(iTmp)%quantity) + '_' + &
              & fu_str(gf%gvars(iTmp)%species)
      elseif(gf%gvars(iTmp)%chCocktailNm /= '') then
        sp%sp = fu_quantity_short_string(gf%gvars(iTmp)%quantity) + '_' + &
              & gf%gvars(iTmp)%chCocktailNm
      else
        sp%sp = fu_quantity_short_string(gf%gvars(iTmp)%quantity)
      endif
      !
      ! If the variable name appears too long, we have to select the shorter one, still obeying 
      ! uniqueness of the names
      !
      if(len_trim(sp%sp) > 15)then
      !
      ! The full variable name is too long. Cut it down to 11 digits and put _XXX at the end
      !
          gf%gvars(iTmp)%chGrADSNm = sp%sp(1:11) + '_' + fu_str(iTmp,3)

      else

        gf%gvars(iTmp)%chGrADSNm = sp%sp  ! name is short enough

      endif  ! if too long variable name

      if(gf%gvars(iTmp)%iVerticalFeature == level_3d_type_flag)then
        WRITE(iCtlUnit,'(A,I4,A9,A,2x,A)')trim(gf%gvars(iTmp)%chGrADSNm), &
                                        & max(gf%gvars(iTmp)%n_levs,0),' 99 99 0 ', &
                                        & TRIM(fu_quantity_string(gf%gvars(iTmp)%quantity)), &
                                        & trim(fu_substance_name(gf%gvars(iTmp)%species))
      else
        WRITE(iCtlUnit,'(A,A12,A,2x,A)') trim(gf%gvars(iTmp)%chGrADSNm), ' 0 99 99 0 ',&
                                       & trim(fu_quantity_string(gf%gvars(iTmp)%quantity)), &
                                       & trim(fu_substance_name(gf%gvars(iTmp)%species))
      END IF
    END DO   ! over variables

    WRITE(iCtlUnit,'(A7)') 'ENDVARS'
    close(iCtlUnit)

    !-------------------------------------------------------------------------
    !
    ! Now write the super-ctl file for SILAM
    ! It contains just three items: name of the ctl file, and complete grid and 
    ! vertical namelists
    !
    open(iCtlUnit, file = gf%fname+'.super_ctl')

    write(iCtlUnit,fmt='(A)')'# This is the super-ctl file allowing SILAM to read the GrADS format'
    !
    ! Generic list contains the universal definitions for the file: grid, vertical, time labeling
    !
    write(iCtlUnit,fmt='(A)')'LIST = general'
    !    write(iCtlUnit,fmt='(2x,2A)')'ctl_file_name = ',trim(gf%fname + '.ctl')
    chTmp=trim(gf%fname + '.ctl')
    write(iCtlUnit,fmt='(2x,2A)')'ctl_file_name = ',trim(fu_trim_grads_hat(chTmp,chTmp))
    if (gf%ifMPIIO) then
      call report_as_namelist(wholeMPIdispersion_grid, iCtlUnit, gf%fname(1:(index(gf%fname, '.grads')-1)))
    else
      call report_as_namelist(gf%silamGrid, iCtlUnit, gf%fname(1:(index(gf%fname, '.grads')-1)))
    end if
    if (.not. defined(gf%silamVertical)) then
      call msg_warning('Undefined vertical, will used surface level', 'write_ctl_file')
      call set_vertical(surface_level, gf%silamVertical)
    end if
    call report_as_namelist(gf%silamVertical, iCtlUnit)
    select case(gf%time_label_position)
      case(end_of_period)
        write(iCtlUnit,*)'time_label_position  =  end_of_period'
      case(start_of_period)
        write(iCtlUnit,*)'time_label_position  =  start_of_period'
      case(mid_of_period)
        write(iCtlUnit,*)'time_label_position  =  mid_of_period'
      case(instant_fields)
        write(iCtlUnit,*)'time_label_position  =  instant_fields'
      case default
        call set_error('Unknown time label position:' + fu_str(gf%time_label_position),'write_ctl_file')
        return
    end select
    write(iCtlUnit,fmt='(A)')'END_LIST = general'
    !
    ! Individual lists for the variables contain the GrADS variable name as the list name
    ! and then the exact definition of these variables one-by-one
    !
    do iTmp =1, gf%n_vars
      write(iCtlUnit,fmt='(2A)')'LIST = ',gf%gvars(iTmp)%chGrADSNm
      write(iCtlUnit,fmt='(3x,2A)')'quantity_short_name = ', &
                                 & trim(fu_quantity_short_string(gf%gvars(iTmp)%quantity))
      if(defined(gf%gvars(iTmp)%species))then
        call report(gf%gvars(iTmp)%species, iCtlUnit)
        !
        ! Molar mass and half-life time, if needed
        !
        write(iCtlUnit,fmt='(3A)')'molar_mass = ', &
                                    & fu_str(fu_mole_mass(fu_material(gf%gvars(iTmp)%species))), ' kg/mole'
        if(fu_if_radioactive(fu_material(gf%gvars(iTmp)%species)))then
           write(iCtlUnit,fmt='(3A)')'half_life = ', &
                                    & fu_str(fu_half_life(fu_material(gf%gvars(iTmp)%species))), ' sec'
        endif
      else
        if(len_trim(gf%gvars(iTmp)%chCocktailNm) > 0) &
                & write(iCtlUnit,fmt='(3x,2A)')'cocktail_name = ',trim(gf%gvars(iTmp)%chCocktailNm)
      endif  ! subst name OK
      !
      ! What to report as level is a tricky question. Grads variables can be 3D, then th ewhole verticacl is the 
      ! choice. But the level of the 2D field should be stored. Can be 2m etc.
      !
      select case(gf%gvars(iTmp)%iVerticalFeature)
        case(level_2D_type_flag, do_nothing_flag)
          write(iCtlUnit,fmt='(3x,A)')'vertical_feature = FIELD_2D'
          call report_as_namelist(iCtlUnit, gf%gvars(iTmp)%SilamLevel, 1.)    ! Useful info
        case(level_3D_type_flag)
          write(iCtlUnit,fmt='(3x,A)')'vertical_feature = FIELD_3D'
        case(integrate_column_flag)
          write(iCtlUnit,fmt='(3x,A)')'vertical_feature = COLUMN_INTEGRATED'
          call report_as_namelist(iCtlUnit, gf%gvars(iTmp)%SilamLevel, 1.)    ! trivial but still, let's store
        case(lowest_level_flag)
          write(iCtlUnit,fmt='(3x,A)')'vertical_feature = LOWEST_LEVEL'
          call report_as_namelist(iCtlUnit, gf%gvars(iTmp)%SilamLevel, 1.)    ! trivial but still, let's store
        case default
          call set_error('Unknown type of vertical feature','write_ctl_file')
      end select

      write(iCtlUnit,fmt='(2A)')'END_LIST = ',gf%gvars(iTmp)%chGrADSNm
    end do  ! over GrADS variables

    close(iCtlUnit)
    !
    ! If we need a duplicate ctl file, the easiest is to copy this one
    !
    if(len_trim(chFixedNameTemplate) > 1)then
      open(iCtlUnit, file = trim(chFixedNameTemplate), iostat = iTmp)
      if(fu_fails(iTmp==0,'Failed to open duplicate ctl file:' + chFixedNameTemplate, &
                        & 'write_ctl_file'))return
      iTmp = fu_next_free_unit()
      if(error)return
      open(iTmp, file = trim(gf%fname +'.ctl')) ! Just written file, should be available
      call copy_text_file(iTmp, iCtlUnit)
      close(iTmp)
      close(iCtlUnit)
    endif

    call free_work_array(sp%sp)

  END SUBROUTINE write_ctl_file

  ! ***************************************************************




!******************************************************************************
!******************************************************************************
!
!   GRIB related stuff
!
!******************************************************************************
!******************************************************************************


  INTEGER FUNCTION init_ctl_for_grib (grib_unit, dir, grib_fname, chTemplate)
    !
    ! Called when GRIB file is opened. Makes the same supplementary work as 
    ! open_grads_file_o: finds empty gfile structure, fills names and units, 
    ! opens ctl file
    !
    ! Author: M.Sofiev
    !
    IMPLICIT NONE
    !
    ! Imported variables with intent IN
    INTEGER, INTENT(in) :: grib_unit
    CHARACTER(LEN=*), INTENT(in) :: dir, grib_fname
    CHARACTER(LEN=*), INTENT(in), optional :: chTemplate

    ! Local variables
    INTEGER :: iFile

    init_ctl_for_grib = -1

    !---------------------------------------------------------
    ! Find a spare file counter. The last structure is always kept free!!
    !
    DO iFile=1, nbr_grads_files
      IF(gfile(iFile)%ptr%defined == silja_false)EXIT
    END DO
    IF(iFile >= nbr_grads_files)THEN
      CALL set_error('Too many grads files','init_ctl_for_grib')
      RETURN
    END IF
    !
    ! So far it is just occupied but not yet defined structure.
    ! It will become silja_true when at least one variable is stored there
    !
    gfile(iFile)%ptr%defined = silja_undefined

    !---------------------------------------------------------
    !
    ! Set up the file names and open the ctl file
    !
    gfile(iFile)%ptr%unit_bin = grib_unit !------------ binary file
    if(len_trim(dir) > 0)then
      gfile(iFile)%ptr%fname = dir + dir_slash + grib_fname
    else
      gfile(iFile)%ptr%fname = grib_fname !grib_fname
    endif

    if(present(chTemplate))then
      gfile(iFile)%ptr%chTemplate = chTemplate
    else
      gfile(iFile)%ptr%chTemplate = gfile(iFile)%ptr%fname
    endif

!    gfile(iFile)%unit_ctl = fu_next_free_unit() !------------ text ctl file
!    OPEN(gfile(iFile)%unit_ctl,file=fu_connect_strings(gfile(iFile)%fname,'.ctl'))

    init_ctl_for_grib = iFile

  END FUNCTION init_ctl_for_grib 


  !******************************************************************

  subroutine set_grib_binary_unit(gIndex, UnitBin)
    !
    ! Just sets new value for the binary unit in the GRIB-related structure
    !
    implicit none

    integer, intent(in) :: gIndex, UnitBin

    gfile(gIndex)%ptr%unit_bin = UnitBin

  end subroutine set_grib_binary_unit


  ! ****************************************************************


  SUBROUTINE fill_structures_for_grib(id,grib_quantity,lev_type,lev_val,t_range,igf)
    !
    ! Fill structures for ctl file describing the GRIB binary. Specifics
    ! are that GRIB can actually handle absolutely any set of the fields. 
    ! As a result, virtually no checking should be done. Algorithm is simple -
    ! if something (time moment, quantity or level) is not in the list - add !
    !
    ! ASSUMPTIONS/LIMITATIONS coming from necessity to show this file in GrADS
    !                          and from the logic of the program.
    ! 1.Time must be regular and motion along time axis must be regular at 
    !    least for the first and the second time peiords.
    ! 2.All 3D variables are treated as multi-layer ones with all layers
    !    defined in a particular file. 
    ! 3.If some quantity appeared at some specific level (e.g. surface) - it
    !    never comes with another one. So, e.g. ground pressure and 3D pressure
    !    are different
    !
    ! This routine has to be called from-inside grib_io since it requires
    ! GRIB internal parameters to be defined.
    !
    ! Author: Mikhail Sofiev
    !
    IMPLICIT NONE

    ! Imported parameters with intent IN
    TYPE(silja_field_id), INTENT(in) :: id
    INTEGER, INTENT(in) :: grib_quantity, & !See grib_code_table module, ksec1(6)
               & lev_type, & ! ksec1(7) - code table 3 WMO GRIB manual
!               & lev_val, &  ! ksec1(8,9) -code table 3 WMO GRIB manual
               & t_range, &  ! time range indicator, ksec1(18), code table 5
               & igf        ! Number of the structure to be filled in
    real, intent(in) :: lev_val

    ! Local variables
    TYPE(grads_file), POINTER :: gf
    INTEGER :: iTmp, iVar
    LOGICAL :: ifTimeExist,ifVarExist,ifLevExist, &
             & if_corner_in_geo_coord,if_south_pole
    REAL :: fTmp, pole_x, pole_y

    gf => gfile(igf)%ptr

    !
    ! If totally new GRIB file - start all lists and fill-in grid (different
    ! from GrADS !!!)
    !
    IF(gf%gtime%defined == silja_false)THEN ! Totally new GRIB file
      gf%gtime%start = fu_valid_time(id)
      gf%gtime%defined = silja_undefined
      gf%time_nbr = 1
      gf%n_times = 1
      gf%gtime%start_bin = gf%gtime%start

      gf%gvars(1)%quantity = fu_quantity(id)
      gf%gvars(1)%species = fu_species(id)
      gf%gvars(1)%ChCocktailNm = fu_cocktail_name(id)
      gf%gvars(1)%grib_quantity = grib_quantity
      gf%gvars(1)%grib_lev_typ = lev_type
      gf%gvars(1)%grib_time_ind = t_range
      gf%gvars(1)%n_levs = 0   !1
      gf%gvars%defined = silja_undefined
      gf%var_nbr = 1
      gf%n_vars = 1

      IF(fu_multi_level_quantity(gf%gvars(1)%quantity))THEN
        !
        ! Potentially multi-level variable
        !
        if(fu_cmp_levs_eq(fu_level(id),entire_atmosphere_integr_level))then
          !
          ! Actually, this is column-integrated field, i.e. 2D from GrADS point of view
          !
          gf%gvars(1)%iVerticalFeature = integrate_column_flag
          gf%gvars(1)%grib_lev_val = lev_val ! level value

        elseif(fu_cmp_levs_eq(fu_level(id), lowest_atmosphere_level))then
          !
          ! Actually, this is lowest-level field, i.e. 2D from GrADS point of view
          !
          gf%gvars(1)%iVerticalFeature = lowest_level_flag
          gf%gvars(1)%grib_lev_val = lev_val ! level value

        else
          gf%glevs%levels(1) = fu_get_glevel(fu_level(id),gf%levType)
          gf%gvars(1)%grib_lev_val = lev_val  ! May be, it is still 2D...
          gf%glevs%defined = silja_undefined
          gf%lev_nbr = 1
          gf%n_levs = 1
        endif
      ELSE
        gf%gvars(1)%iVerticalFeature = level_2d_type_flag ! This is a 2D variable, no levels
        gf%gvars(1)%grib_lev_val = lev_val ! level value
      END IF

      SELECT CASE(fu_gridtype(fu_grid(id))) ! Define the grid for ctl file
      CASE(lonlat)
        CALL lonlat_grid_parameters(fu_grid(id), &
                           & gf%ggrid%x_start, gf%ggrid%y_start,  &
                           & if_corner_in_geo_coord,&
                           & gf%ggrid%nx, gf%ggrid%ny, &
                           & pole_x, pole_y, & 
                           & gf%ggrid%x_step, gf%ggrid%y_step)
      CASE DEFAULT
        CALL set_error('Only latlon grids so far','open_grads_file_o')
        RETURN
      END SELECT
      gf%defined = silja_true
      RETURN
    END IF


    !-------------------------------------------------
    !
    ! Does time period exist in the list ?
    !
    DO iTmp = 1, gf%n_times
      IF(iTmp > 1)THEN
        ifTimeExist = (gf%gtime%start + gf%gtime%step * REAL(iTmp-1)) ==  &
                    &  fu_valid_time(id)
      ELSE
        ifTimeExist = gf%gtime%start == fu_valid_time(id)
      END IF
      IF(ifTimeExist) EXIT
    END DO
    IF(.not. ifTimeExist) THEN
      gf%n_times =gf%n_times + 1
      IF(gf%n_times == 2) gf%gtime%step = fu_valid_time(id) - gf%gtime%start
    END IF


    !-------------------------------------------------
    !
    ! Does variable exist in the list ?
    !
    DO iVar = 1, gf%n_vars
      ifVarExist = gf%gvars(iVar)%quantity == fu_quantity(id) .and. &
                 & gf%gvars(iVar)%species == fu_species(id) .and. &
                 & gf%gvars(iVar)%ChCocktailNm == fu_cocktail_name(id) 
      IF(ifVarExist)EXIT
    END DO

    IF(ifVarExist)THEN
      if(gf%gvars(iVar)%grib_time_ind /= t_range) gf%gvars(iVar)%grib_time_ind = -999
    ELSE
      gf%n_vars = gf%n_vars + 1
      iVar = gf%n_vars
      gf%gvars(gf%n_vars)%quantity = fu_quantity(id)
      gf%gvars(gf%n_vars)%species = fu_species(id)
      gf%gvars(gf%n_vars)%ChCocktailNm = fu_cocktail_name(id) 
      gf%gvars(gf%n_vars)%grib_quantity = grib_quantity
      gf%gvars(gf%n_vars)%grib_lev_typ = lev_type
      gf%gvars(gf%n_vars)%grib_time_ind = t_range

      gf%gvars(gf%n_vars)%grib_lev_val = lev_val

      IF(fu_multi_level_quantity(gf%gvars(gf%n_vars)%quantity))THEN
        if(fu_cmp_levs_eq(fu_level(id),entire_atmosphere_integr_level))then
          gf%gvars(gf%n_vars)%iVerticalFeature = integrate_column_flag
        elseif(fu_cmp_levs_eq(fu_level(id), lowest_atmosphere_level))then
          gf%gvars(gf%n_vars)%iVerticalFeature = lowest_level_flag
        else
          gf%gvars(gf%n_vars)%iVerticalFeature = level_3d_type_flag
        endif
      ELSE
        gf%gvars(gf%n_vars)%iVerticalFeature = level_2d_type_flag
      END IF
    END IF

    !---------------------------------------------------
    !
    ! Does level exist in the list ?
    !
    IF(fu_multi_level_quantity(fu_quantity(id)))THEN !3D => do levels
      fTmp = fu_get_glevel(fu_level(id),gf%levType)
      IF(error)RETURN
      ifLevExist = .false.
      DO iTmp = 1, gf%n_levs
        ifLevExist = gf%glevs%levels(iTmp) .eps. fTmp
        IF(ifLevExist)EXIT
      END DO

      IF(.not. ifLevExist) THEN 
        !
        ! We have to insert the level to the place where it should be
        ! Problem is that nobody promised the order - ascending or descending
        ! Solution: First two layers follow each other and define the order
        !
        IF(gf%n_levs == 0)THEN
          gf%glevs%levels(1) = fTmp
        ELSEIF(gf%n_levs == 1)THEN
          gf%glevs%levels(2) = fTmp
        ELSE
          IF(gf%glevs%levels(1) < gf%glevs%levels(gf%n_levs))THEN
            ! Ascending order
            DO iTmp = gf%n_levs, 1, -1
              IF(gf%glevs%levels(iTmp) < fTmp)THEN
                gf%glevs%levels(iTmp+1) = fTmp
                EXIT
              ELSE
                gf%glevs%levels(iTmp+1) = gf%glevs%levels(iTmp)
              END IF
            END DO
          ELSE
            ! Descending order
            DO iTmp = gf%n_levs, 1, -1
              IF(gf%glevs%levels(iTmp) > fTmp)THEN
                gf%glevs%levels(iTmp+1) = fTmp
                EXIT
              ELSE
                gf%glevs%levels(iTmp+1) = gf%glevs%levels(iTmp)
              END IF
            END DO
          END IF
        END IF
        gf%n_levs = gf%n_levs + 1
      END IF  !--- insert level
    END IF !-- if 3D

  END SUBROUTINE fill_structures_for_grib



  ! ****************************************************************


  SUBROUTINE write_ctl_for_grib(igf)
    !
    ! Writes the ctl file describing the GRIB binary. There are a few
    ! tricks:
    ! 1. GRIB codes are written for each variable
    ! 2. It is stated that all 3D variables have maximum number of
    !    levels. This is not necessary the case, but corrections are
    !    to be done by gribmap routine
    !
    ! Author: M.Sofiev
    !
    IMPLICIT NONE
    !
    ! Imported parameters with the intent IN
    INTEGER, INTENT(in) :: igf

    ! Local variables
    INTEGER :: iTmp, iCtlUnit
    CHARACTER(LEN = fnlen) :: chTmp
    TYPE(grads_file), POINTER :: gf

    gf => gfile(igf)%ptr

    call set_error('Ctl for grib is outdated','write_ctl_for_grib')
    return

    iCtlUnit = fu_next_free_unit()
    if(error)then
      call set_error('Failed to get unit for ctl file','write_ctl_file_for_grib')
      return
    endif
    open(iCtlUnit,file=fu_connect_strings(gf%fname,'.ctl'))

    if(.not. gf%defined == silja_true)then
      WRITE(iCtlUnit,*) 'ALL STRUCTURES UNDEFINED'
      call msg_warning('All GRIB structures are undefined','write_ctl_for_grib')
      return
    endif
        
    ! WRITE(iCtlUnit,'(A5,A)')'DSET ',gf%chTemplate(1:len_trim(gf%chTemplate))

    chTmp = trim(gf%chTemplate)
    WRITE(iCtlUnit,'(A5,A)')'DSET ',fu_trim_grads_hat(chTmp,gf%fname)
    
    WRITE(iCtlUnit,'(A23)') 'TITLE SILAM GRIB output'
    WRITE(iCtlUnit,'(A10)') 'DTYPE grib'
    if(index(gf%chTemplate, '%') > 0) write(iCtlUnit,'(A16)') 'OPTIONS TEMPLATE'
    chTmp = fu_connect_strings(gf%fname,'.gmp')
    WRITE(iCtlUnit,'(A6,A)') 'INDEX ', chTmp(1:LEN_TRIM(chTmp))
    WRITE(iCtlUnit,'(A6,1X,E15.6)') 'UNDEF ', &
                                 & fu_grib_missing_real(gf%gvars(1)%quantity)
    WRITE(iCtlUnit,'(A5,I5,A8,2(F15.7,1x))') 'XDEF ',gf%ggrid%nx,' LINEAR ', &
                                                & gf%ggrid%x_start, &
                                                & gf%ggrid%x_step
    WRITE(iCtlUnit,'(A5,I5,A8,2(F15.7,1x))') 'YDEF ',gf%ggrid%ny,' LINEAR ', &
                                                & gf%ggrid%y_start, &
                                                & gf%ggrid%y_step
    IF(gf%n_levs == 0)THEN
      WRITE(iCtlUnit,'(A15)') 'ZDEF 1 LEVELS 0'
    ELSE
      WRITE(iCtlUnit,'(A5,I4,A8,100(F15.7,1x))') 'ZDEF ',gf%n_levs,' LEVELS ', &
                                  & (gf%glevs%levels(iTmp),iTmp=1,gf%n_levs)
    END IF
    if(gf%n_times == 1)then
      WRITE(iCtlUnit,'(A16,A20,A4)') 'TDEF  1  LINEAR ', &
                            & fu_time_to_grads_string(gf%gtime%start), ' 1hr'
    else
      WRITE(iCtlUnit,'(A5,I7,A8,A20,A10)') 'TDEF ',gf%n_times,' LINEAR ', &
                            & fu_time_to_grads_string(gf%gtime%start), &
                            & fu_interval_to_grads_string(gf%gtime%step)
    endif

    WRITE(iCtlUnit,'(A5,I5)')'VARS ',gf%n_vars
    DO iTmp =1, gf%n_vars
      if(defined(gf%gvars(iTmp)%species))then
        chTmp = fu_quantity_short_string(gf%gvars(iTmp)%quantity) + &
                                 & '_' + fu_substance_name(gf%gvars(iTmp)%species)
      else
        chTmp = fu_quantity_short_string(gf%gvars(iTmp)%quantity)
      endif
      if(len_trim(chTmp) > 15)then
        call msg_warning('Cutting too long variable name','write_ctl_file')
        chTmp(15:)=' '
      endif
      IF(gf%gvars(iTmp)%iVerticalFeature == level_3d_type_flag)THEN
        WRITE(iCtlUnit,'(A,A4,2(I6,A2),F9.2,A2,I6,2x,A,2x,A)') &
            & chTmp(1:LEN_TRIM(chTmp)), &            ! Var name
            & gf%n_levs, &                         ! Nbr of levels (MAX available)
            & gf%gvars(iTmp)%grib_quantity,', ', &  ! GRIB code of quantity
            & gf%gvars(iTmp)%grib_lev_typ,', ', &   ! GRIB level type
            & gf%gvars(iTmp)%grib_lev_val,', ', &   ! GRIB level value 
            & gf%gvars(iTmp)%grib_time_ind, &       ! GRIB time range indicator
            & TRIM(fu_quantity_string(gf%gvars(iTmp)%quantity)), & ! complete name
            & trim(fu_str(gf%gvars(iTmp)%species))
      ELSE
        WRITE(iCtlUnit,'(A,I4,2(I6,A2),A5,I5,2x,A,2x,A)') &
             & chTmp(1:LEN_TRIM(chTmp)), &          ! Var name
             & '0 ', &                               ! 2D => Nbr of levels =0
             & gf%gvars(iTmp)%grib_quantity,', ', & ! GRIB code of quantity
             & gf%gvars(iTmp)%grib_lev_typ,', ', &  ! GRIB level type
             & '0, ', &                             ! 3D => GRIB level value = 0
             & gf%gvars(iTmp)%grib_time_ind, &     ! GRIB time range indicator
             & TRIM(fu_quantity_string(gf%gvars(iTmp)%quantity)), & ! complete name
             & trim(fu_str(gf%gvars(iTmp)%species))
      END IF
    END DO
    WRITE(iCtlUnit,'(A7)') 'ENDVARS'

    close(iCtlUnit)

  END SUBROUTINE write_ctl_for_grib


  !********************************************************************************

  subroutine report_grads_file(gf)
    !
    ! Simply prints the main content of the GrADS file structure
    !
    implicit none

    ! Imported parameter
    TYPE(grads_file), POINTER :: gf

    ! Local variable
    integer :: iVar

    call msg('============ report for GrADS file:' + gf%fname + ', vars:', gf%n_Vars)
    do iVar = 1, gf%n_Vars
      call msg('Species and quantity:' + fu_str(gf%gvars(iVar)%species), &
             & gf%gvars(iVar)%quantity)
    end do
    call msg('Variable current index:',gf%var_nbr)
    call msg('Grid nx, x_start:',gf%ggrid%nx,gf%ggrid%x_start)
    call msg('Grid ny, y_start:',gf%ggrid%ny,gf%ggrid%y_start)
    call msg('Grid x-step:', gf%ggrid%x_step)
    call msg('Grid y-step:', gf%ggrid%y_step)
    do iVar = 1, gf%n_Levs
      call msg('Level:',iVar,gf%glevs%levels(iVar))
    end do
    call msg('Current level index:',gf%lev_nbr)
    call msg('Start time & current index:' + fu_time_to_io_string(gf%gtime%start),gf%time_nbr)
    call msg('============ End of report')

  end subroutine report_grads_file


!********************************************************************************
!********************************************************************************
!
!   Encapsulation stuff
!
!********************************************************************************
!********************************************************************************

  INTEGER FUNCTION fu_unit_bin(ifg)
    IMPLICIT NONE
    INTEGER, INTENT(in) :: ifg
    fu_unit_bin = gfile(ifg)%ptr%unit_bin
  END FUNCTION fu_unit_bin


  SUBROUTINE release_index(gIndex)
    !
    ! Frees-up the grads_file structure in gfile array
    !
    IMPLICIT NONE

    INTEGER, INTENT(in) :: gIndex

    gfile(gIndex)%ptr = gfile(nbr_grads_files)%ptr ! The last structure is always undefined

  END SUBROUTINE release_index


  function fu_grads_filename(iFile)

    implicit none

    integer, intent(in) :: iFile
    CHARACTER (LEN=fnlen) :: fu_grads_filename
    
    fu_grads_filename = gfile(iFile)%ptr%fname

  end function fu_grads_filename
 

  function fu_grads_sctl_filename(iFile)

    implicit none

    integer, intent(in) :: iFile
    CHARACTER (LEN=fnlen) :: fu_grads_sctl_filename

    fu_grads_sctl_filename = gfile(iFile)%ptr%super_ctl_fname

  end function fu_grads_sctl_filename
 
  real function fu_grads_missing_value(iFile)
    implicit none
    integer, intent(in) :: iFile
    fu_grads_missing_value = gfile(iFile)%ptr%missing_value
  end function fu_grads_missing_value

  ! Write buffering. The data are normally written into file one field at time. Especially
  ! in MPI runs, each write becomes small, and performance can be increased by first
  ! aggregating data from multiple fields. For serial runs buffering doesn't necessary
  ! improve the performance, however.
  ! 
  ! Current default buffer size is 50 fields, but if the number of MPI tasks is high,
  ! much more buffering can be useful.
  ! 
  
  integer function fu_get_default_buf_size() result(buf_size_flds)
    implicit none
    integer, parameter :: buf_size_flds_def = 50
    integer :: stat
    character(len=*), parameter :: sub_name = 'fu_get_default_buf_size'
    character(len=clen) :: strval

    call get_environment_variable('GRADSIO_BUF_SIZE', strval)
    if (error) return
    if (strval /= '') then
      read(unit=strval, fmt=*, iostat=stat) buf_size_flds
      if (fu_fails(stat == 0, 'Failed to parse GRADSIO_BUF_SIZE: ' + strval, sub_name)) return
      if (fu_fails(buf_size_flds < 1000000, 'Suspicious GRADSIO_BUF_SIZE', sub_name)) return
      if (buf_size_flds < 2) call msg_warning('GRADSIO_BUF_SIZE: ' + fu_str(buf_size_flds), sub_name)
      ! size < 1 used for disabling buffer (above)
    else
      buf_size_flds = buf_size_flds_def
    end if
    
  end function fu_get_default_buf_size

  subroutine init_grads_buffer(buf, fieldsize, buf_size_flds)
    implicit none
    type(grads_buffer) :: buf
    integer, intent(in) :: fieldsize, buf_size_flds

    integer :: stat
    character(len=*), parameter :: sub_name = 'init_grads_buffer'
    
    if (fu_fails(buf_size_flds > 0, 'Bad buf_size_flds', sub_name)) return
        
    call msg('Allocating grads buffer, num flds:', buf_size_flds)
    allocate(buf%buf(fieldsize, buf_size_flds), stat=stat)
    if (fu_fails(stat == 0, 'Failed to allocate write buffer', sub_name)) return
    buf%max_flds  = buf_size_flds
    buf%fieldsize = fieldsize
    buf%allocated = .true.
    buf%num_flds = 0

    call msg('Grads buffer allocated, fieldsize:', fieldsize)

  end subroutine init_grads_buffer

  subroutine free_grads_buffer(buf)
    implicit none
    type(grads_buffer), intent(inout) :: buf
    
    if (fu_fails(buf%allocated, 'buffer not allocate', 'free_grads_buffer')) return
    deallocate(buf%buf)
    buf%max_flds = 0
    buf%num_flds = 0
    buf%allocated = .false.
  end subroutine free_grads_buffer
 
  subroutine field_to_buffer(buf, grid_data, fieldsize)
    implicit none
    type(grads_buffer), intent(inout) :: buf
    integer, intent(in) :: fieldsize
    real, dimension(:), intent(in) :: grid_data
    
    integer :: ind_buffer
    character(len=*), parameter :: sub_name = 'field_to_buffer'

    if (.not. buf%allocated) then
      call init_grads_buffer(buf, fieldsize, fu_get_default_buf_size())
      if (error) return
    end if

    if (fu_fails(buf%fieldsize == fieldsize, 'Wrong fieldsize in buffer', sub_name)) return
    if (fu_fails(.not. fu_buffer_full(buf), 'Buffer full', sub_name)) return
    
    buf%buf(:,buf%num_flds+1) = grid_data(1:buf%fieldsize)
    buf%num_flds = buf%num_flds + 1

    !call msg('Field stored in buffer, index, count', ind_buffer, buf%num_flds)
  end subroutine field_to_buffer
  
  logical function fu_buffer_full(buf) result(if_full)
    implicit none
    type(grads_buffer), intent(in) :: buf
    if_full = buf%num_flds == buf%max_flds
  end function fu_buffer_full

  subroutine flush_buffer(igf)
    ! Flush the write buffer. The buffer needs not to be full. The is incremented by the
    ! number of fields written.
    ! 
    implicit none
    integer, intent(in) :: igf
    
    character(len=*), parameter :: sub_name = 'flush_buffer'
    type(grads_file), pointer :: gf
    integer :: ind_start, ind_end

    call start_count('flush_buffer')

    if (fu_fails(igf > 0 .and. igf <= size(gfile), 'Invalid igf', sub_name)) return
    gf => gfile(igf)%ptr
    if (fu_fails(associated(gf), 'gf not associated', sub_name)) return
    
    if (.not. gf%gradsbuf%allocated) return
    if (gf%gradsbuf%num_flds == 0) return
    
    !ind_start = gf%gradsbuf%max_flds - gf%gradsbuf%num_flds + 1
    !ind_end = gf%gradsbuf%max_flds
    ind_start = 1
    ind_end = gf%gradsbuf%num_flds
!       call msg('Flushing write buffer, indices', ind_start, ind_end)
    if (gf%ifMPIIO) then
      call smpi_write_grads_fieldset_parallel(gf%gradsbuf%buf(:,ind_start:ind_end), gf%gradsbuf%fieldsize, &
                                            & gf%gradsbuf%num_flds, gf%rec, gf%unit_bin)
    else
      call write_grads_fieldset(gf%gradsbuf%buf(:,ind_start:ind_end), &
                              & gf%gradsbuf%fieldsize, gf%gradsbuf%num_flds, &
                              & gf%unit_bin)
    end if
    if (error) return

    !NB. When buffering is used, the record counter is used (only) by MPI-IO.
    !
    gf%rec = gf%rec + gf%gradsbuf%num_flds
    gf%gradsbuf%num_flds = 0
!    call msg('Done')
    call stop_count('flush_buffer')
  end subroutine flush_buffer

END MODULE grads_io
