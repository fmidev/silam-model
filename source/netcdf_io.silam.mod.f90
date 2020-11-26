MODULE netcdf_io

  ! This module contains all necessary information and routines for 
  ! reading/writing of NETCDF files.
  !
  ! Writing:
  ! Grid, vertical and variables are all defined while opening the file for writing.
  ! Two vertical axis' are defined: normal one for 3D variables and 1 level for 2D variables.
  ! Time is treated as record dimension with unlimited length. 
  ! Timesteps don't have to be regular, time value is written to file at each new timestep.
  ! For inverse runs the binary is not inverted and times go backwards. 
  ! If output time split is not ALL_IN_ONE, something like ctl is created to enable GrADS to 
  ! handle the multiple files with template filenames like one. In that case regular timestep
  ! is assumed.
  ! Conventions followed in writing netcdf files are as close to CF as we could manage.
  !
  ! Author: Marje Prank, FMI marje.prank@fmi.fi
  !
  ! All units: SI
  ! 
  ! Language: ANSI standard Fortran 90
  !

! Following macros allow to enable different features  
! WITH_PNETCDF



  !USE grads_io
!  use grib_api_io
  use netcdf
  use stacks !fields_3d !_identifications
  use silam_partitioning

#ifdef WITH_PNETCDF
  use pnetcdf, only : nf90mpi_buffer_attach,nf90mpi_buffer_detach,nf90mpi_create,nf90mpi_def_dim,nf90mpi_def_var,nf90mpi_enddef,nf90mpi_put_att
  use pnetcdf, only : nf90mpi_put_var,nf90mpi_bput_var,nf90mpi_strerror,nf90mpi_wait_all,nf90mpi_begin_indep_data,nf90mpi_end_indep_data, nf90mpi_close
#endif


  implicit none

  PRIVATE   ! cut out all the lengthy NetCDF libraries

  public init_netcdf_io
  public open_netcdf_file_o     ! Opens a netcdf file for writing
  public close_netcdf_file      ! Closes a netcdf file
  public write_next_field_to_netcdf_file  ! Writes one 2D variable to a netcdf file
  public fu_next_free_netcdf_structure 
  public write_ctl_file_4_netcdf

  public open_netcdf_file_i
  public id_list_from_netcdf_file
  public get_netcdf_grids   !Used in sources and grads2io
  public get_netcdf_verticals ! Used in sources and grads2io and init_massmap
  public get_netcdf_total   !  Used in sources 
  public read_field_from_netcdf_file
  public field_indices_from_netcdf_file
  public read_field_from_netcdf_indices

  public timelst_from_netcdf_file

  public fu_netcdf_filename
  public expand_output_list

  private make_staggered_fld
  private set_missing_netcdf_dim
  private fu_netcdf_str_to_silam_time
  private put_var_nc
  private put_var_nc_buf
  private def_dim_nc
  private def_var_nc
  private put_chatt_nc
  private put_fatt_nc
  private create_nc
  
  interface set_missing
     module procedure set_missing_netcdf_dim
  end interface

  private set_missing


#ifdef DEBUG_NC   
    logical, private, parameter :: ifMsgs = .True.
#else
    logical, private, parameter :: ifMsgs = .False.
#endif

  !----------------------------------------------------------------------------
  !  netcdf-defining structures. 
  !


  TYPE netcdf_grid
    PRIVATE
    REAL :: x_start=-1., y_start=-1., x_step=-1., y_step=-1.
    INTEGER :: nx=-1,ny=-1
    character (len = nf90_max_name) :: chYdimName, chXdimName
    integer :: latDimId = int_missing, latVarId = int_missing, lonDimId = int_missing, lonVarId = int_missing
    integer :: dxMVarId = int_missing, dyMVarId = int_missing, &
         & sinMapRotVarId = int_missing, cosMapRotVarId = int_missing
    real :: sPole_lon = 0.0, sPole_lat = -90.0 ! For rotated lon-lat grid
    type(silja_grid) :: sGrid = grid_missing
    TYPE(silja_logical) :: defined = silja_false
 
  END TYPE netcdf_grid

  TYPE netcdf_levels
    PRIVATE
    real, dimension(max_levels) :: levels
    real, dimension(max_levels) :: std_z ! Needed for non-z heights
    character (len = nf90_max_name) :: CHZDIMNAME="", zUnit="", positive="", long_name="", standard_name=""
    character (len = nf90_max_name) :: CHZhalfDIMNAME=""
    integer :: levVarId = int_missing, levDimId = int_missing
    integer :: aVarId = int_missing,  bVarId = int_missing   !For hybrid
    integer :: daVarId = int_missing, dbVarId = int_missing   !For hybrid layers
    integer :: dzVarID = int_missing
    integer :: ahalfVarId = int_missing,  bhalfVarId = int_missing, levhalfDimId =int_missing !For hybrid layers
    TYPE(silja_logical) :: defined = silja_false
  END TYPE netcdf_levels

  TYPE netcdf_time
    PRIVATE
    TYPE(silja_time) :: first_valid_time=time_missing, last_valid_time = time_missing, &
                     & tDimStart = time_missing, analysis_time = time_missing
    TYPE(silja_time), dimension(:), pointer :: times => null()
    TYPE(silja_interval) :: step = interval_missing
    character (len = nf90_max_name) :: CHTDIMNAME = ''
    logical :: ifRegular = .false.
    integer :: tDimId, tVarId, time_label_position, analysis_t_from
    TYPE(silja_logical) :: defined = silja_false
  END TYPE netcdf_time

  TYPE netcdf_variable
    PRIVATE
    INTEGER :: quantity=int_missing, n_levs=int_missing, n_dims=int_missing
    integer, dimension(nf90_max_var_dims):: dimIds = int_missing

    type(silam_species) :: species
    character(len=nf90_max_name) :: chVarNm = ''
    character(len=substNmLen) :: chCocktailNm = ''
    LOGICAL :: if3D      = .false., ifTimeDep = .false.
    integer :: varid = int_missing
    integer :: xtype = int_missing
    type(netcdf_grid), pointer :: nGrid => null() !Pointer to files grid
    type(silam_vertical):: sVert
    TYPE(silja_logical) :: defined = silja_false
    character(len=nf90_max_name) ::  unit = "", silam_amount_unit = ""
    real :: nc_molar_mass = -1.
    real :: offset=0., scaleFactor=1.
    real :: offset_nc = 0., scaleFactor_nc = 1.
    character(len = nf90_max_name) :: stagger = '', axis='', positive_direction='', calendar=''

    real :: missing_value=real_missing, valid_min=real_missing, valid_max=real_missing

    character(len = nf90_max_name) :: standard_name='', long_name='', coordinates='', bounds='', compress=''
 
  END TYPE netcdf_variable
  
  type netcdf_dim
    PRIVATE
    character :: axis='u', positive='m'  ! x, y, z, t, u-unknown; positive: u-up, d-down, m-missing  
    character(len=nf90_max_name) :: dimName='', varName='', chType=''
    character(len=nf90_max_name) :: PsVar=''
    real, dimension(1) :: p0=1
    real(8), dimension(:), pointer :: values => null()
    real, dimension(:), pointer :: a=>NULL(), b=>NULL(), a_half=>NULL(), b_half=>NULL(), dz=>null()
    real, dimension(:, :), pointer :: values2d =>null()
    type(silam_vertical) :: svert 
    character(len = nf90_max_name) :: A_var='', B_var='', P0_var='', PS_var='', A_half_var='', &
       & B_half_var='', dz_var = ''
    integer :: dimId=int_missing, varId=int_missing, dimLen=int_missing
    logical :: defined = .false., ifTimeDep=.false., if2d=.false., ifReverse=.false.

  end type netcdf_dim




  !---------------------------------------------------------------------
  ! Structure completely defining a single netcdf file  
  !
  INTEGER, public, PARAMETER :: max_nbr_of_grids_in_file = 5
  TYPE netcdf_file
    PRIVATE
    type(netcdf_grid), dimension(max_nbr_of_grids_in_file) :: nGrids
    CHARACTER (LEN=fnlen) :: fname='', chTemplate='', title='',tmpfname=''
    type(meteo_data_source) :: met_src = met_src_missing
    INTEGER :: unit_bin=-1
    type(netcdf_dim), dimension(:), pointer :: nDims =>null()
    TYPE(netcdf_levels) :: nlevs
    type(silam_vertical) :: silamVertical  
    TYPE(netcdf_variable), DIMENSION(:), pointer :: nvars =>null()
    TYPE(silja_field_id), dimension(max_2d_fields) :: id
    TYPE(netcdf_time) :: ntime 
    INTEGER :: n_times=0, n_vars=0, n_levs=0, n_Dims=0 ! Numbers of aliases
    real :: missing_value = real_missing
#ifdef WITH_PNETCDF
    integer, dimension(:), allocatable ::  nf90mpiReq,  nf90mpiStat
    integer :: nf90mpifieldCnt, nf90mpiBUFfdlds
#endif
    logical :: ifCocktails = .false.
    logical :: ifMPIIO = .false.
    integer :: ncver = int_missing
    TYPE(silja_logical) :: defined = silja_false

  END TYPE netcdf_file



  !----------------------------------------------------------------------
  ! Structure for ctl-like creature for netcdf
  !
  type netcdf_ctl
    private
    TYPE(silja_time) :: start_time , last_valid_time 
    TYPE(silja_interval) :: step   
    integer :: nTimes 
    logical :: defined
  end type netcdf_ctl

  type(netcdf_ctl), private, save :: nCtl = netcdf_ctl(time_missing, time_missing, interval_missing, int_missing,.false.)

  !----------------------------------------------------------------------
  !
  !Some parameters
  !
  INTEGER, public, PARAMETER :: max_nbr_of_netcdf_files = 200 

  integer, public, parameter :: dim_start = 11101
  integer, public, parameter :: first_value = 11102

  INTEGER, private, PARAMETER :: nc4_deflate_level = 5

  TYPE(netcdf_file), DIMENSION(max_nbr_of_netcdf_files), PRIVATE, target, SAVE :: nfile

  !---------------------------------------------------------------------
  !
  ! The main list of output variables. Moved here from io_server
  !
  type TOutputLstItem
!    private
    integer :: quantity, request, AvType, iSpeciesListType, iVerticalTreatment
    logical :: if3D
    type(silja_interval) :: AvPeriod = interval_missing
    type(silam_species) :: species 
    character(len=clen) :: chSpecies_string
    type(silja_field_id) :: targetId 
  end type TOutputLstItem

  type (TOutputLstItem), parameter, public :: OutputLstItem_missing = &
        & TOutputLstItem(int_missing, int_missing, int_missing, int_missing, int_missing,&
        & .False., interval_missing, species_missing, '', field_id_missing)

  type TOutputList
!    private
    type(TOutputLstItem), dimension(:), allocatable ::ptrItem
  end type TOutputList

  type (TOutputList), parameter, public :: OutputList_missing = TOutputList(null())

  public TOutputLstItem
  public TOutputList

!--------------------------------------------------------------------------
  type(Tsilam_namelist_group), private, pointer, save :: nmTblNlGrp => null()

  integer, parameter :: buffer_timesteps = 1 !!! 24 !!!How many timesteps to buffer
                                             !!  for pnetcdf


CONTAINS

  !
  ! Wrappers-dispatchers
  !
  !==================================================================================
  subroutine put_chatt_nc(nf, varID, varname, attNname, att)
     ! A wrapper to nf90_put_att 
    type (netcdf_file), intent(in) :: nf
    integer, intent(in) :: varID
    character(len=*), intent(in) :: varname, attNname, att

    character (len=nf90_max_name) ::  chAttName
    character (len=fnlen) :: chAtt
    integer :: iStat
    chAttName = attNname
    chAtt = att
    if (nf%ncver==3 .and. nf%ifmpiio) then
#ifdef  WITH_PNETCDF
      iStat = nf90mpi_put_att(nf%unit_bin, varId, chAttName, att) 
      if (iStat /= 0) call msg("nf90mpi_put_att: "//nf90mpi_strerror(iStat))
#else
      call set_error("Compiled without pNetCDF",'put_chatt_nc')
      return
#endif
    else
      iStat = nf90_put_att(nf%unit_bin, varId, chAttName, chAtt) 
      if (iStat /= 0) call msg("nf90_put_att: "//nf90_strerror(iStat))
    endif
    if (iStat /= 0)then 
      call set_error('Failed to put "'+attNname+'" attribute to "'+&
                 &  varname+'" variable','put_chatt_nc')
    endif
  end subroutine put_chatt_nc
   
  !========================================================================================
  subroutine put_fatt_nc(nf, varID, varname, attNname, att)
     ! A wrapper to nf90_put_att
    type (netcdf_file), intent(in) :: nf
    integer, intent(in) :: varID
    character(len=*), intent(in) :: varname, attNname
    real, intent(in) :: att

    character (len=nf90_max_name) ::  chAttName
#ifdef  WITH_PNETCDF
    real (kind =  r4k) :: fAtt4
#endif
    real :: fAtt
    integer :: iStat
    chAttName = attNname
    if (nf%ncver==3 .and. nf%ifmpiio) then
#ifdef  WITH_PNETCDF
      fAtt4 = att
      iStat = nf90mpi_put_att(nf%unit_bin, varID, chAttName, fAtt4)

      if (iStat /= 0) call msg("nf90mpi_put_att: "//nf90mpi_strerror(iStat))
#else
      call set_error("Compiled without pNetCDF",'put_chatt_nc')
      return
#endif
    else
      fAtt = att
      iStat = nf90_put_att(nf%unit_bin, varID, chAttName, fAtt) 
      if (iStat /= 0) call msg("nf90_put_att: "//nf90_strerror(iStat))
    endif
    if (iStat /= 0)then 
      call set_error('Failed to put "'+attNname+'" attribute to "'+&
                 &  varname+'" variable','put_fhatt_nc')
    endif
  end subroutine put_fatt_nc

  !===========================================================================================
  subroutine create_nc(fname, ounit, ncver, ifMPIIO)
     character (len=*), intent (in) :: fname
     integer, intent(inout) :: ounit
     integer, intent(in) :: ncver
     logical, intent(in) :: ifmpiio
     

     integer :: iStat
     character (len=*), parameter ::sub_name = "create_nc"

     if (ifmpiio) then
       call msg("Creating mpiio ncver:", ncver)
     else
       call msg("Creating non-mpiio ncver:", ncver)
     endif

    iStat = 0

    if (ncver == 3) then
      if (ifMPIIO) then
#ifdef WITH_PNETCDF
      iStat = nf90mpi_create(smpi_adv_comm, fname, NF90_64BIT_OFFSET, MPI_INFO_NULL, ounit)
      if (iStat/=0) call set_error("nf90mpi_create error: "//trim(nf90mpi_strerror(iStat)), sub_name)
      
#else
      call set_error("Tried to use MPIIO  nc3 without PNETCDF support", sub_name)
#endif
      else
        iStat = nf90_create(fname, NF90_64BIT_OFFSET, ounit)
        if (iStat/=0) call set_error("nf90_create1 error: "//trim(nf90_strerror(iStat)), sub_name)
      endif
   elseif (ncver == 4) then
      if (ifMPIIO) then
        call set_error('Tried to use NETCDF4_MPIIO while compiled without it',sub_name)
      else
        iStat = nf90_create(fname, NF90_NETCDF4, ounit)
      endif
      if (iStat/=0) call set_error("nf90_create2 error: "//trim(nf90_strerror(iStat)), sub_name)
    else
       call set_error('Strange ncver: '//trim(fu_str(ncver)),sub_name)
    endif
  end subroutine create_nc

  !=================================================================================

  subroutine def_dim_nc(nf, dim_name, dimlen, dim_id)
     type (netcdf_file), intent(in) :: nf
     character (len=*), intent (in) :: dim_name
     integer, intent(in) :: dimlen
     integer, intent(out) :: dim_id

     integer :: iStat
     character (len=*), parameter ::sub_name = "def_dim_nc"

    iStat = int_missing 
    if (nf%ncver==3 .and. nf%ifmpiio) then
#ifdef WITH_PNETCDF      
        iStat = nf90mpi_def_dim(nf%unit_bin, dim_name, int(dimlen, kind=8), dim_id) 
        if (iStat /= 0) call msg_warning("nf90mpi_def_dim: "//nf90mpi_strerror(iStat), sub_name)
#endif
    else
        iStat = nf90_def_dim(nf%unit_bin, dim_name, dimlen, dim_id) 
        if (iStat /= 0) call msg_warning("nf90_def_dim: "//nf90_strerror(iStat), sub_name)
    endif
    if (iStat /= 0)then 
      call set_error('Failed to define "'+dim_name+'" dimension',sub_name)
    endif
  end subroutine def_dim_nc

  !=================================================================================

  subroutine def_var_nc(nf, var_name, xType, dim_ids, chunks, varId, ifCollective)
     type (netcdf_file), intent(in) :: nf
     character (len=*), intent (in) :: var_name
     integer,  intent(in) :: xType
     integer, dimension(:), intent(in) :: dim_ids, chunks
     logical, intent(in) :: ifCollective
     integer, intent(out) :: varId

     integer :: iStat
     character (len=*), parameter ::sub_name = "def_var_nc"

    iStat = int_missing 
    if (nf%ncver==3 .and. nf%ifmpiio) then
#ifdef WITH_PNETCDF
     iStat = nf90mpi_def_var(nf%unit_bin, var_name, xType,  dim_ids, varId) 
     if (iStat /= 0) call set_error("nf90mpi_def_var: "//nf90mpi_strerror(iStat), sub_name)
#else             
       call set_error("Compiled without PNETCDF", sub_name)
#endif
    else
      if (nf%ncver== 3 .or. size(chunks) < 3) then
        iStat = nf90_def_var(nf%unit_bin, var_name, xType,  dim_ids, varId) 
      else   
#ifdef VS2012
        iStat = nf90_def_var(nf%unit_bin, var_name, xType, dim_ids, varId)
#else
        if (nf%ifmpiio) then !No compression
          iStat = nf90_def_var(nf%unit_bin, var_name, xType, dim_ids, varId, chunksizes = chunks)
        else !Make chunked, compressed variable
          iStat = nf90_def_var(nf%unit_bin, var_name, xType, dim_ids, varId, &
                             & chunksizes = chunks, shuffle=.false., deflate_level = nc4_deflate_level)
        endif
#endif
      endif
      if (iStat /= 0) call set_error("nf90_def_var: "//nf90_strerror(iStat), sub_name)
    endif
    if (iStat /= 0)then 
      call set_error('Failed to create variable "'+var_name+'"',sub_name)
    endif
  end subroutine def_var_nc

  !=================================================================================

  subroutine put_var_nc(nf, varId, outData, starts, counts)
      !
      ! Put simple variables -- dimensions etc....
      !
     type (netcdf_file),  intent(in) :: nf
     integer, intent(in) :: varId
     real(kind=4), dimension(:), intent(inout) :: outData
     integer, dimension (:), intent(in) :: starts, counts

     integer(kind=8), dimension (nf90_max_var_dims) :: mystarts, mycounts
     integer :: ndims

     integer :: iStat
     character (len=*), parameter ::sub_name = "put_var_nc"

     if (nf%ncver == 3  .and. nf%ifmpiio) then
#ifdef WITH_PNETCDF
        ndims=size(starts)
        mystarts(:ndims) = starts(:)
        mycounts(:ndims) = counts(:)
        iStat = nf90mpi_put_var(nf%unit_bin, VarId, outData, mystarts(1:ndims), mycounts(1:ndims))
        if (iStat /= 0) call set_error("nf90mpi_put_var: "//nf90mpi_strerror(iStat), sub_name)
#else
       call set_error("Compiled without PNETCDF", sub_name)
#endif
     else
        iStat = nf90_put_var(nf%unit_bin, VarId, outData, start=starts, count=counts)
        if (iStat /= 0) call set_error("nf90_put_var: "//nf90_strerror(iStat), sub_name)
     endif
  end subroutine put_var_nc

  !=================================================================================

  subroutine put_var_nc_buf(nf, varId, outData, starts, counts)
    ! Uses buffered feature of nf90mpi if availablae
     type (netcdf_file), intent(inout) :: nf
     integer, intent(in) :: varId
    ! For some reason nf90mpi_bput_var uses intent(inout) for data
     real(kind=4), dimension(:), intent(inout) :: outData
     integer, dimension(:), intent(in) :: starts, counts

     integer(kind=8), dimension (nf90_max_var_dims) :: mystarts, mycounts
     integer :: ndims
     integer :: iStat, iTmp

     character (len=*), parameter ::sub_name = "put_var_nc_buf"

     if (nf%ncver == 3  .and. nf%ifmpiio) then
#ifdef WITH_PNETCDF 
      nf%nf90mpifieldCnt =  nf%nf90mpifieldCnt + 1
      ndims=size(starts)
      mystarts(:ndims) = starts(:)
      mycounts(:ndims) = counts(:)
!      call msg("size(outData)",size(outData))
!      call msg("mystarts(1:ndims)",starts(1:ndims))
!      call msg("mycounts(1:ndims)",counts(1:ndims))

      if (ifMsgs) call msg("nf%nf90mpifieldCnt nf%nf90mpiBUFfdlds",nf%nf90mpifieldCnt, nf%nf90mpiBUFfdlds)
      iStat = nf90mpi_bput_var(nf%unit_bin, VarId, outData, nf%nf90mpiReq(nf%nf90mpifieldCnt), mystarts(1:ndims), mycounts(1:ndims))
      if (iStat /= 0) call set_error("nf90mpi_put_var: "//nf90mpi_strerror(iStat), sub_name)

      ! Was it the last field? Flush it
      if (nf%nf90mpifieldCnt == nf%nf90mpiBUFfdlds ) then ! Extra for time variable
        call flush_nf90mpi_buf(nf, sub_name)
        if (fu_fails(.not. error,"Error set", sub_name)) return
      endif
#else
       call set_error("Compiled without PNETCDF", sub_name)
#endif
     else
        iStat = nf90_put_var(nf%unit_bin, VarId, outData,  starts, counts)
        if (iStat /= 0) call set_error("nf90_put_var: "//nf90_strerror(iStat), sub_name)
     endif
  end subroutine put_var_nc_buf


  !*******************************************


#ifdef WITH_PNETCDF
  subroutine flush_nf90mpi_buf(nf, place)
        type (netcdf_file), intent(inout) :: nf
        character (len=*), intent(in) :: place
        character (len=*), parameter ::sub_name = "flush_nf90mpi_buf"

        integer :: iStat, iTmp

        if (ifMsgs)call msg("Dumping all",nf%nf90mpifieldCnt, nf%nf90mpiBUFfdlds)

        iStat = nf90mpi_wait_all(nf%unit_bin, nf%nf90mpifieldCnt, nf%nf90mpiReq, nf%nf90mpiStat)

        if (iStat /= 0) call set_error("nf90mpi_wait_all in "//place//": "//nf90mpi_strerror(iStat), sub_name)

        if (any(nf%nf90mpiStat(1:nf%nf90mpifieldCnt) /= 0)) then ! Trouble writing something
          do iTmp = 1, nf%nf90mpifieldCnt
            call msg("nf90mpi_status after nf90mpi_wait_all in "//place//" ("//trim(fu_str(iTmp))//"): "// &
                    & trim(nf90mpi_strerror(nf%nf90mpiStat(iTmp))))
          enddo
          call set_error("Trouble with status", sub_name)
          return
        endif
        nf%nf90mpifieldCnt = 0
        if(ifMsgs) call msg("Dump complete")


  end subroutine flush_nf90mpi_buf

#endif
 !
 !
 ! Main routines
 !




  !*******************************************

  integer function open_netcdf_file_o(fname, myOutGrid, vertical, timeValid, &
                                    & lstsOutVars, &
                                    & chTemplate, ifAllInOne, ncver, ifMPIIO, fMissingVal)
    
    ! Opens the NETCDF files, fills nFile structure, defines dimensions and variables
    ! Returns the number of the file in the netcdf file list nfile_ptr
    ! Called for every output file
    !

    IMPLICIT NONE
    !
    ! Imported parameters with intent IN
    CHARACTER (LEN=*), INTENT(in) :: fname ! Directory and name of the GrADS file
    TYPE(silja_grid), INTENT(in) :: myOutGrid ! This MPI domain output grid
                                        ! We will get proper one to store in
                                        ! netcdf structures
    type(silam_vertical), INTENT(in) :: vertical
    type(silja_time), INTENT(in) :: timeValid
    type(TOutputList), dimension(3), INTENT(in) :: lstsOutVars
    CHARACTER (LEN=*), INTENT(in) :: chTemplate
    logical, intent(in) :: ifAllInOne, ifMPIIO
    integer, intent(in) :: ncver
    real, intent(in) :: fMissingVal

    ! Local declarations
    INTEGER :: iFile, iStat,iVal, iTmp, outvar_quantity, iLst, jTmp
    character (len=nf90_max_name) ::  chAttName
    TYPE(silja_grid) :: grid !Full output grid
    character (len=fnlen) :: chAtt
    real(r4k) :: fAtt
    integer, dimension(4), target :: dimids3d, chunks3d
    integer, dimension(3),target :: dimids2d,  chunks2d
    integer ::  nx, ny, offx, offy, gnx, gny
    integer, dimension(:), pointer :: dimidsXd, chunksXd !Pointer to one of above
    integer :: NoFields2d ! total 2d fields in the file
    LOGICAL :: if_corner_in_geo_coord, if_rotated_pole
    logical :: if_ab_needed, if_ab_half_needed, if_dz_needed
    REAL :: pole_x, pole_y
    real, dimension(:), pointer :: varData
    real, dimension(:), pointer :: fPtr,fPtra,fPtrb,fPtrda,fPtrdb,fPtrahalf,fPtrbhalf
    real :: fTmp

    type(Tsilam_namelist), pointer :: nl
    type(Tsilam_nl_item_ptr),pointer :: nlitem
    TYPE(netcdf_file),POINTER :: nf 
    integer, dimension(0) :: zeroints 
    character (len=*), parameter :: sub_name='open_netcdf_file_o'

    nf => null()
    dimidsXd => null()
    chunksXd => null()
    NoFields2d = 0
    !
    ! Open the file, find structure for it
    !
    if (ifMPIIO) then 
      if (.not.smpi_is_mpi_version()) then
        call set_error('Trying to open file in parallel mode - this is &
                     & a serial version', 'open_netcdf_file_o')
        return
      endif
      grid = wholeMPIdispersion_grid
    else
      grid = myOutGrid
    end if
    

    open_netcdf_file_o = int_missing
    
    ! Find a spare file counter. The last structure is always kept free!!
    iFile = fu_next_free_netcdf_structure()
    nf => nfile(iFile)
    if(error)return

    nf%defined = silja_undefined

    nf%fname = fname
    nf%tmpfname = nf%fname+".tmp"



    if(ifAllInOne)then                           
      nf%chTemplate = nf%fname
    else
      nf%chTemplate = chTemplate
    endif

    call create_directory_tree(nf%fname(1:index(nfile(iFile)%fname,dir_slash,.true.)-1))

    nf%ncver=ncver 
    nf%ifMPIIO = ifMPIIO
    call create_nc(nf%tmpfname, nf%unit_bin, nf%ncver, nf%ifMPIIO)
    if (error) then
      call set_error('Failed to create netcdf file','open_netcdf_file_o')
      return
    endif
    open_netcdf_file_o = iFile
  

    if(.not. nCtl%defined) then
       ! if first time write start time to nctl structure -- this is the start time of the run
        nCtl%nTimes = 0
        nCtl%start_time = timeValid
        nCtl%step = interval_missing
        nCtl%defined = .true.
   endif
   
   nl => fu_create_namelist()
    
   ! And now write the headers (attributes, axes, variables)
   ! First global attributes:
    call  put_chatt_nc(nf, nf90_global, "global", "title", "SILAM_OUTPUT")
    call  put_chatt_nc(nf, nf90_global, "global", "Conventions",  "CF-1.3")
    call  put_chatt_nc(nf, nf90_global, "global", "source",  revision_str)
    call  put_chatt_nc(nf, nf90_global, "global", "_CoordinateModelRunDate", &
           &  fu_time_to_thredds_string(nCtl%start_time))
    call  put_chatt_nc(nf, nf90_global, "global", "SIMULATION_START_DATE", &
           &  fu_time_to_thredds_string(nCtl%start_time))

    nf%missing_value = fMissingVal
    
    ! Define dimensions:
    ! Grid:

    nf%nGrids(1)%sGrid = grid

    call grid_dimensions(grid, nf%nGrids(1)%nx, nf%nGrids(1)%ny)

    if(fu_gridtype(grid) == lonlat)then
      CALL lonlat_grid_parameters(grid, &
            & nf%nGrids(1)%x_start, nf%nGrids(1)%y_start, if_corner_in_geo_coord,&
            & nf%nGrids(1)%nx, nf%nGrids(1)%ny, &
            & pole_x, pole_y, & 
            & nf%nGrids(1)%x_step, nf%nGrids(1)%y_step)

      ! Silam grid rotation spec
      call  put_chatt_nc(nf, nf90_global, "global", "grid_projection",  "lonlat")
      call  put_fatt_nc(nf, nf90_global, "global", "pole_lat", pole_y)
      call  put_fatt_nc(nf, nf90_global, "global", "pole_lon", pole_x)
      
      !CF grid rotation spec
      if   (pole_y .eps. -90.)  then 
          nf%nGrids(1)%chYdimName = "lat" 
          nf%nGrids(1)%chXdimName = "lon" 
          if_rotated_pole = .false.
      else
         nf%nGrids(1)%chYdimName = "rlat" 
         nf%nGrids(1)%chXdimName = "rlon" 
          if_rotated_pole = .true.
         call def_var_nc(nf, "rp", nf90_char, zeroints, zeroints , iTmp, .false.) !varId=iTmp, no collective
         if (error) return

         
         call  put_chatt_nc(nf, iTmp, "rp", "grid_mapping_name",  "rotated_latitude_longitude")
         call  put_fatt_nc(nf, iTmp, "rp", "grid_north_pole_latitude", -pole_y)
         call  put_fatt_nc(nf, iTmp, "rp", "grid_north_pole_longitude", mod(pole_x+360., 360.) - 180.)
         call  put_chatt_nc(nf, iTmp, "rp", "_CoordinateTransformType",  "Projection")
         call  put_chatt_nc(nf, iTmp, "rp", "_CoordinateAxisTypes",  "GeoX GeoY")
      endif
    else
      nf%nGrids(1)%chYdimName = "Y" 
      nf%nGrids(1)%chXdimName = "X" 
      call  put_chatt_nc(nf, nf90_global, "global", "grid_projection",  "unknown")
    endif
    
    nf%nGrids(1)%defined = silja_true

    call def_dim_nc(nf, nf%nGrids(1)%chXdimName, nf%nGrids(1)%nx, nf%nGrids(1)%lonDimId)
    if (fu_fails(.not. error, "def_dim_nc X failed", sub_name)) return

    call def_dim_nc(nf, nf%nGrids(1)%chYdimName, nf%nGrids(1)%ny, nf%nGrids(1)%latDimId) 
    if (fu_fails(.not. error, "def_dim_nc Y failed", sub_name)) return


    ! Vertical:
    nf%silamVertical = vertical
    nf%n_levs = fu_NbrOfLevels(vertical)
    if_ab_needed = .false.
    if_ab_half_needed = .false.
    if_dz_needed = .false.

    SELECT CASE(fu_leveltype(vertical))
    CASE(surface)
      nf%nlevs%chZdimName = "surface"
      nf%nlevs%long_name = "surface level"
      nf%nlevs%levels(1) = 0.0 ! Whatever, e.g. metres
      nf%nlevs%zUnit = "level"
      nf%nlevs%positive = "up"
    CASE(constant_height)
      nf%nlevs%chZdimName = "height"
      nf%nlevs%long_name = "constant height from surface"
      nf%nlevs%zUnit = "m"
      nf%nlevs%positive = "up"
      do iVal = 1, nf%n_levs
        nf%nlevs%levels(iVal) = fu_level_height(fu_level(vertical, iVal))  ! Metres
      enddo
    CASE(constant_altitude)
      nf%nlevs%chZdimName = "altitude"
      nf%nlevs%long_name = "constant altitude from mean sea level"
      nf%nlevs%zUnit = "m"
      nf%nlevs%positive = "up"
      do iVal = 1, nf%n_levs
        nf%nlevs%levels(iVal) = fu_level_altitude(fu_level(vertical, iVal)) ! Metres
      enddo
    CASE(constant_pressure)
      nf%nlevs%chZdimName = "pressure"
      nf%nlevs%long_name = "constant pressure"
      nf%nlevs%zUnit = "Pa"
      nf%nlevs%positive = "down"
      do iVal = 1, nf%n_levs
        nf%nlevs%levels(iVal) = fu_pr_level_pressure(fu_level(vertical, iVal)) ! Pa
      enddo
    CASE(layer_btw_2_height)
      nf%nlevs%chZdimName = "height"
      nf%nlevs%long_name = "layer midpoint constant height from surface"
      nf%nlevs%standard_name = 'layer_midpoint_height_above_ground'
      nf%nlevs%zUnit = "m"
      nf%nlevs%positive = "up"
      if_dz_needed = .True.
      do iVal = 1, nf%n_levs
        nf%nlevs%levels(iVal) = fu_layer_centre_value(fu_level(vertical, iVal)) ! Metres
      enddo
    CASE(layer_btw_2_altitude)
      nf%nlevs%chZdimName = "altitude"
      nf%nlevs%long_name = "layer midpoint constant altitude from mean sea level"
      nf%nlevs%zUnit = "m"
      nf%nlevs%positive = "up"
      if_dz_needed = .True.
      do iVal = 1, nf%n_levs
        nf%nlevs%levels(iVal) = fu_layer_centre_value(fu_level(vertical, iVal)) ! Pa
      enddo
    CASE(layer_btw_2_pressure)
      nf%nlevs%chZdimName = "pressure"
      nf%nlevs%long_name = "layer midpoint constant pressure"
      nf%nlevs%zUnit = "Pa"
      nf%nlevs%positive = "down"
      do iVal = 1, nf%n_levs
        nf%nlevs%levels(iVal) = fu_layer_centre_value(fu_level(vertical, iVal)) ! Pa
      enddo
    CASE(hybrid)
      nf%nlevs%chZdimName = "hybrid"
      nf%nlevs%long_name = "Standard height of hybrid level"
      nf%nlevs%zUnit = "m"
      nf%nlevs%positive = "up"
      if_ab_needed = .true.
      do iVal = 1, nf%n_levs
        nf%nlevs%levels(iVal) = fu_hybrid_level_number(fu_level(vertical, iVal)) ! UNITLESS
      enddo
    CASE(layer_btw_2_hybrid)
      nf%nlevs%chZdimName = "hybrid"
      nf%nlevs%chZhalfdimName = "hybrid_half"
      nf%nlevs%long_name = "standard height of hybrid layer midpoint"
      nf%nlevs%standard_name = "hybrid_layer_midpoint_standard_height_above_ground"
      nf%nlevs%zUnit = "m"
      nf%nlevs%positive = "up"
      if_ab_needed = .true.
      if_ab_half_needed = .true.
      do iVal = 1, nf%n_levs
        ! These values will be used to match field that comes to existing netcdf level
        nf%nlevs%levels(iVal) = fu_hybrid_level_number(fu_level(vertical, iVal)) ! UNITLESS
      enddo
    CASE DEFAULT
      call msg("Leveltype", fu_leveltype(vertical)) 
      call report(vertical)
      CALL set_error('Cant handle vertical','open_netcdf_file_o')
      RETURN
    END SELECT

    nf%nlevs%defined = silja_true

    call def_dim_nc(nf, nf%nlevs%chZdimName, nf%n_levs, nf%nlevs%levDimId) 
    if (fu_fails(.not. error, "def_dim_nc failed z", sub_name)) return

    if (if_ab_half_needed) then  ! One more dim for half-levels
       call def_dim_nc(nf, nf%nlevs%chZhalfdimName,  nf%n_levs+1, nf%nlevs%levhalfDimId) 
       if (fu_fails(.not. error, "def_dim_nc failed ab_half", sub_name)) return
    endif

    ! Time:
    nf%ntime%first_valid_time = timeValid
    nf%ntime%last_valid_time = time_missing
    if (nf%ntime%first_valid_time == time_missing ) then
      call set_error('Failed to get start time','open_netcdf_file_o')
      return
    endif

    nf%ntime%defined = silja_true    

    nf%ntime%chTdimName = "time"  

    call def_dim_nc(nf, nf%ntime%chTdimName, nf90_unlimited, nf%ntime%tDimId)
    if (fu_fails(.not. error, "def_dim_nc failed time", sub_name)) return

    dimIds2d =  (/nf%ngrids(1)%lonDimId, nf%ngrids(1)%latDimId, nf%ntime%tDimId/)
    dimIds3d =  (/nf%ngrids(1)%lonDimId, nf%ngrids(1)%latDimId, nf%nlevs%levDimId ,nf%ntime%tDimId/)


    ! Makes difference only for NetCDF4
    if ( ncver == 4 .and.  ifMPIIO) then
      !NC4 and MPIIO is very slow otherwise... Makes nearly unreadable files on large grids
      chunks3d =  (/nf%ngrids(1)%nx, nf%ngrids(1)%ny, nf%n_levs , 1/)
      chunks2d =  (/nf%ngrids(1)%nx, nf%ngrids(1)%ny, 1 /)
    else !Matters only for NetCDF4
      !Smaller chunks -- faster reading of fields afterwards
      !Minimal chunk size -- 90 kbytes, max 360 kBytes: large enough for
      !efficient compression
      iTmp = (nf%ngrids(1)%nx - 1) / 300  + 1 !Parts to divide to 300 points per chunk max
      iTmp = (nf%ngrids(1)%nx - 1) / iTmp + 1 !chunk size
      jTmp = (nf%ngrids(1)%ny - 1) / 300  + 1 
      jTmp = (nf%ngrids(1)%ny - 1) / jTmp + 1
      chunks3d =  (/iTmp, jTmp, 1 , 1/)
      chunks2d =  (/iTmp, jTmp, 1 /)
    endif
    
       

    !
    ! Define variables, assign variable attributes
    ! Dimension variables:
    

    SELECT CASE(fu_gridtype(grid))

    CASE(lonlat)
      call def_var_nc(nf, nf%ngrids(1)%chXdimName, nf90_float, (/nf%ngrids(1)%lonDimId/), (/nf%ngrids(1)%nx/), nf%ngrids(1)%lonVarId, .false.)
      if (fu_fails(.not. error, "def_var_nc failed x1d", sub_name)) return

      call def_var_nc(nf, nf%ngrids(1)%chYdimName, nf90_float, (/nf%ngrids(1)%latDimId/), (/nf%ngrids(1)%ny/), nf%ngrids(1)%latVarId, .false.)
      if (fu_fails(.not. error, "def_var_nc failed y1d", sub_name)) return

    case(anygrid)
      call def_var_nc(nf, 'lon2d', nf90_float, (/nf%ngrids(1)%lonDimId, nf%ngrids(1)%latDimId/),  (/nf%ngrids(1)%nx,nf%ngrids(1)%ny/), nf%ngrids(1)%lonVarId, .false.)
      if (fu_fails(.not. error, "def_var_nc failed x2d", sub_name)) return

      call def_var_nc(nf, 'lat2d', nf90_float, (/nf%ngrids(1)%lonDimId, nf%ngrids(1)%latDimId/), (/nf%ngrids(1)%nx,nf%ngrids(1)%ny/),  nf%ngrids(1)%latVarId, .false.)
      if (fu_fails(.not. error, "def_var_nc failed y2d", sub_name)) return
    case default
      call set_error('unknown gridtype', sub_name)
    end select

    call put_chatt_nc(nf, nf%ngrids(1)%lonVarId, "x dim ",  "units", "degrees_east")
    call put_chatt_nc(nf, nf%ngrids(1)%lonVarId, "x dim ",  "axis", "X")
    call put_chatt_nc(nf, nf%ngrids(1)%latVarId, "y dim ",  "units", "degrees_north")
    call put_chatt_nc(nf, nf%ngrids(1)%latVarId, "y dim ",  "axis", "Y")

    if (if_rotated_pole) then
       call put_chatt_nc(nf, nf%ngrids(1)%lonVarId, "x dim ",  "long_name", "longitude in rotated pole grid")
       call put_chatt_nc(nf, nf%ngrids(1)%lonVarId, "x dim ",  "standard_name", "grid_longitude")
       call put_chatt_nc(nf, nf%ngrids(1)%lonVarId, "x dim ",  "_CoordinateAxisType", "GeoX")
       call put_chatt_nc(nf, nf%ngrids(1)%latVarId, "y dim ",  "long_name", "latitude in rotated pole grid")
       call put_chatt_nc(nf, nf%ngrids(1)%latVarId, "y dim ",  "standard_name", "grid_latitude")
       call put_chatt_nc(nf, nf%ngrids(1)%latVarId, "y dim ",  "_CoordinateAxisType", "GeoY")
    else
       call put_chatt_nc(nf, nf%ngrids(1)%lonVarId, "x dim ",  "long_name", "longitude")
       call put_chatt_nc(nf, nf%ngrids(1)%lonVarId, "x dim ",  "standard_name", "longitude")
       call put_chatt_nc(nf, nf%ngrids(1)%lonVarId, "x dim ",  "_CoordinateAxisType", "Lon")

       call put_chatt_nc(nf, nf%ngrids(1)%latVarId, "y dim ",  "long_name", "latitude")
       call put_chatt_nc(nf, nf%ngrids(1)%latVarId, "y dim ",  "standard_name", "latitude")
       call put_chatt_nc(nf, nf%ngrids(1)%latVarId, "y dim ",  "_CoordinateAxisType", "Lat")
    endif

    !
    ! Grid of anygrid type needs some extra grid variables
    !
    if(fu_gridtype(grid)==anygrid)then

      call def_var_nc(nf, 'dx', nf90_float, (/nf%ngrids(1)%lonDimId, nf%ngrids(1)%latDimId/), (/nf%ngrids(1)%nx,nf%ngrids(1)%ny/), nf%ngrids(1)%dxMVarId, .false.)
      if (fu_fails(.not. error, "def_var_nc failed dx", sub_name)) return

      call def_var_nc(nf, 'dy', nf90_float, (/nf%ngrids(1)%lonDimId, nf%ngrids(1)%latDimId/), (/nf%ngrids(1)%nx,nf%ngrids(1)%ny/), nf%ngrids(1)%dyMVarId, .false.)
      if (fu_fails(.not. error, "def_var_nc failed dy", sub_name)) return

      call def_var_nc(nf, 'sin_map_rot', nf90_float, (/nf%ngrids(1)%lonDimId, nf%ngrids(1)%latDimId/), (/nf%ngrids(1)%nx,nf%ngrids(1)%ny/), nf%ngrids(1)%sinMapRotVarId, .false.)
      if (fu_fails(.not. error, "def_var_nc failed sin_map_rot", sub_name)) return

      call def_var_nc(nf, 'cos_map_rot', nf90_float, (/nf%ngrids(1)%lonDimId, nf%ngrids(1)%latDimId/), (/nf%ngrids(1)%nx,nf%ngrids(1)%ny/), nf%ngrids(1)%cosMapRotVarId, .false.)
      if (fu_fails(.not. error, "def_var_nc failed cos_map_rot", sub_name)) return


      call put_chatt_nc(nf, nf%ngrids(1)%dxMVarId, "dx",  "units", "m")
      call put_chatt_nc(nf, nf%ngrids(1)%dyMVarId, "dy",  "units", "m")

      call put_chatt_nc(nf, nf%ngrids(1)%sinMapRotVarId, &
            & "sin of local map rotation",  "long_name", "sine of local map rotation")

      call put_chatt_nc(nf, nf%ngrids(1)%cosMapRotVarId, &
            & "cos of local map rotation",  "long_name", "cosine of local map rotation")
      endif

    !
    ! Z variable
    !

    call def_var_nc(nf, nf%nlevs%chZdimName, nf90_float, (/nf%nlevs%levDimId/), (/nf%n_levs/), nf%nlevs%levVarId, .false.)
    if (fu_fails(.not. error, "def_var_nc failed z dim", sub_name)) return

    call put_chatt_nc(nf, nf%nlevs%levVarId, "z dim",  "units", nf%nlevs%zUnit)
    call put_chatt_nc(nf, nf%nlevs%levVarId, "z dim",  "positive", nf%nlevs%positive)
    call put_chatt_nc(nf, nf%nlevs%levVarId, "z dim",  "long_name", nf%nlevs%long_name)
    call put_chatt_nc(nf, nf%nlevs%levVarId, "z dim",  "axis", "Z")
    call put_chatt_nc(nf, nf%nlevs%levVarId, "z dim",  "standard_name",nf%nlevs%standard_name )

    if (if_dz_needed) then
       call def_var_nc(nf, "dz", nf90_float, (/nf%nlevs%levDimId/), (/nf%n_levs/), nf%nlevs%dzVarId, .false.)
       if (fu_fails(.not. error, "def_var_nc failed dz", sub_name)) return

       call put_chatt_nc(nf, nf%nlevs%dzVarId, "dz var",  "long_name", "Layer thickness" )
       call put_chatt_nc(nf, nf%nlevs%dzVarId, "dz var",  "standard_name", "layer_thickness" )
       call put_chatt_nc(nf, nf%nlevs%dzVarId, "dz var",  "units", "m" )
    endif

    
    if (if_ab_needed) then
       call def_var_nc(nf, "a", nf90_float, (/nf%nlevs%levDimId/), (/nf%n_levs/), nf%nlevs%aVarId, .false.)
       if (fu_fails(.not. error, "def_var_nc failed a", sub_name)) return

       call def_var_nc(nf, "b", nf90_float, (/nf%nlevs%levDimId/), (/nf%n_levs/), nf%nlevs%bVarId, .false.)
       if (fu_fails(.not. error, "def_var_nc failed b", sub_name)) return

       call put_chatt_nc(nf, nf%nlevs%aVarId, "a var",  "long_name", "layer mid-point A(N) coefficients" )
       call put_chatt_nc(nf, nf%nlevs%aVarId, "a var",  "units", "Pa" )
       call put_chatt_nc(nf, nf%nlevs%bVarId, "b var",  "long_name", "layer mid-point B(N) coefficients" ) 
       call put_chatt_nc(nf, nf%nlevs%bVarId, "b var",  "units", "1" )
       if (if_ab_half_needed) then

          call def_var_nc(nf, "a_half", nf90_float, (/nf%nlevs%levhalfDimId/), (/nf%n_levs+1/), nf%nlevs%ahalfVarId, .false.)
          if (fu_fails(.not. error, "def_var_nc failed a_half", sub_name)) return

          call def_var_nc(nf, "b_half", nf90_float, (/nf%nlevs%levhalfDimId/), (/nf%n_levs+1/), nf%nlevs%bhalfVarId, .false.) 
          if (fu_fails(.not. error, "def_var_nc failed b_half", sub_name)) return

          call def_var_nc(nf, "da", nf90_float, (/nf%nlevs%levDimId/), (/nf%n_levs/), nf%nlevs%daVarId, .false.)
          if (fu_fails(.not. error, "def_var_nc failed da", sub_name)) return

          call def_var_nc(nf, "db", nf90_float, (/nf%nlevs%levDimId/), (/nf%n_levs/), nf%nlevs%dbVarId, .false.) 
           if (fu_fails(.not. error, "def_var_nc failed db", sub_name)) return

          call put_chatt_nc(nf, nf%nlevs%daVarId, "da var",  "long_name", "layer increment of 'A' coefficients" )
          call put_chatt_nc(nf, nf%nlevs%daVarId, "da var",  "units", "Pa" )
          call put_chatt_nc(nf, nf%nlevs%dbVarId, "db var",  "long_name", "layer increment of 'B' coefficients" )
          call put_chatt_nc(nf, nf%nlevs%dbVarId, "db var",  "units", "1" )
          call put_chatt_nc(nf, nf%nlevs%ahalfVarId, "a_half var",  "long_name", "layer top A(N) coefficients" )
          call put_chatt_nc(nf, nf%nlevs%ahalfVarId, "a_half var",  "units", "Pa" )
          call put_chatt_nc(nf, nf%nlevs%bhalfVarId, "b_half var",  "long_name", "layer top B(N) coefficients" )
          call put_chatt_nc(nf, nf%nlevs%bhalfVarId, "b_half var",  "units", "1" )
        endif !if_ab_half_needed
    endif !if_ab_needed


    !
    ! Time
    !

!    netcdf4 mpiio must have collective access to the record dimension variable
    call def_var_nc(nf, nf%ntime%chTdimName, nf90_int, (/nf%ntime%tDimId/), (/1/), nf%ntime%tVarId, .true.) 
    if (fu_fails(.not. error, "def_var_nc failed t dim", sub_name)) return


    call put_chatt_nc(nf, nf%ntime%tVarId, "time var",  "units", &
            & "seconds since " //fu_time_to_netcdf_string(nCtl%start_time) )
    call put_chatt_nc(nf, nf%ntime%tVarId, "time var",  "long_name", "time")
    call put_chatt_nc(nf, nf%ntime%tVarId, "time var",  "axis",  "T")
    call put_chatt_nc(nf, nf%ntime%tVarId, "time var",  "calendar",  "standard")
    call put_chatt_nc(nf, nf%ntime%tVarId, "time var",  "standard_name",  "time")


      
    ! Output variables:
    nf%n_vars = 0
    !
    ! There can be three lists: meteorological variables, dispersion variables from stack
    ! and dispersion mass-map variables. Let's scan all three lists
    !
    ! No separate vars for dimensions

    allocate(nf%nvars(size(lstsOutVars(1)%ptrItem) &
                    & + size(lstsOutVars(2)%ptrItem)&
                    & + size(lstsOutVars(3)%ptrItem)), stat=istat)
    if (istat /= 0)then 
      call msg('Failed to allocate NetCDF file structure No:', iFile)
      call msg ("size(lstsOutVars(1)%ptrItem)+ size(lstsOutVars(2)%ptrItem)+size(lstsOutVars(3)%ptrItem)", &
              & (/size(lstsOutVars(1)%ptrItem), size(lstsOutVars(2)%ptrItem), size(lstsOutVars(3)%ptrItem)/))
      !msg puts stars if there is something strange
      print *, size(lstsOutVars(1)%ptrItem), size(lstsOutVars(2)%ptrItem), size(lstsOutVars(3)%ptrItem)    
      call set_error('failed to allocate nf%ndims(nf%n_Dims), nf%nvars(nf%n_Vars)',sub_name)
      return
    endif  

    do iLst = 1, 3   ! meteo, dispersion, mass map
      if(.not. allocated(lstsOutVars(iLst)%ptrItem))cycle
      if  (size(lstsOutVars(iLst)%ptrItem) < 1) cycle
      do iTmp=1, size(lstsOutVars(iLst)%ptrItem)
        if(lstsOutVars(iLst)%ptrItem(iTmp)%quantity == int_missing .or. &
         & lstsOutVars(iLst)%ptrItem(iTmp)%quantity < 1)exit
        outvar_quantity = lstsOutVars(iLst)%ptrItem(iTmp)%quantity
        if(outvar_quantity == int_missing .or. outvar_quantity < 1) exit
        nf%n_vars = nf%n_vars + 1
        nf%nvars(nf%n_vars)%quantity = outvar_quantity
        nf%nvars(nf%n_vars)%if3D = lstsOutVars(iLst)%ptrItem(iTmp)%if3D
        if (defined(lstsOutVars(iLst)%ptrItem(iTmp)%species)) then
          nf%nvars(nf%n_vars)%species = lstsOutVars(iLst)%ptrItem(iTmp)%species
          nf%nvars(nf%n_vars)%chVarNm = fu_quantity_short_string(outvar_quantity) + &
                                      & '_' + fu_str(nf%nvars(nf%n_vars)%species)
        else
          nf%nvars(nf%n_vars)%species = species_missing
          nf%nvars(nf%n_vars)%chVarNm = fu_quantity_short_string(outvar_quantity)
        end if

        if (nf%nvars(nf%n_vars)%if3D) then
           dimidsXd => dimids3d
           chunksXd => chunks3d
           NoFields2D = NoFields2D + nf%n_levs
           if (ifMsgs) call msg('Adding 3D variable:' + fu_quantity_string(outvar_quantity) + '_' + fu_str(nf%nvars(nf%n_vars)%species))
        else
           dimidsXd => dimids2d
           chunksXd => chunks2d
           NoFields2D = NoFields2D + 1
           if (ifMsgs) call msg('Adding 2D variable:' + fu_quantity_string(outvar_quantity) + '_' + fu_str(nf%nvars(nf%n_vars)%species))
        endif

        call def_var_nc(nf, nf%nvars(nf%n_vars)%chVarNm, nf90_float, &
           & dimidsXd, chunksXd, nf%nvars(nf%n_vars)%varId, nf%ifmpiio)
        if (fu_fails(.not. error, "def_var_nc main vars", sub_name)) return

        call put_fatt_nc(nf, nf%nvars(nf%n_vars)%varId, nf%nvars(nf%n_vars)%chVarNm, &
                       & "_FillValue", nf%missing_value)

!        if (.not. nf%nvars(nf%n_vars)%if3D) then
!            call level_to_short_string(fu_level(lstsOutVars(iLst)%ptrItem(iTmp)%targetId), chAtt)
!            call put_chatt_nc(nf, nf%nvars(nf%n_vars)%varId, nf%nvars(nf%n_vars)%chVarNm, &
!                     & 'vert_level'  , trim(chAtt))
!        endif


        if(defined(nf%nvars(nf%n_vars)%species)) then
          call reset_namelist(nl)
          call put_species_to_namelist(nf%nvars(nf%n_vars)%species, nl)
          do jTmp=1,fu_nbr_of_items(nl)
            nlitem =>  fu_get_item(nl, jTmp)
            call put_chatt_nc(nf, nf%nvars(nf%n_vars)%varId, nf%nvars(nf%n_vars)%chVarNm, &
                     & fu_name(nlitem) , fu_content(nlitem))
          enddo
          ! Some unit handling: declare the unit. Note that scale_factor is always 1.0
          ! SILAM provides the output in the same units it keeps the computations
          chAtt = fu_content(nl,"silam_amount_unit")
          nf%nvars(nf%n_vars)%unit = fu_quantity_unit(outvar_quantity, chAtt)
          !
          ! Molar mass and half-life time, if needed
          !
          fTmp = fu_mole_mass(fu_material(nf%nvars(nf%n_vars)%species))
          if (fTmp > 0 .and. fTmp < 10000) then
             chAtt = fu_str(fTmp) // ' kg/mole'
             call put_chatt_nc(nf, nf%nvars(nf%n_vars)%varId, nf%nvars(nf%n_vars)%chVarNm, &
                             & "molar_mass", chAtt)
          endif

          if(fu_if_radioactive(fu_material(nf%nvars(nf%n_vars)%species)))then
            chAtt = fu_str(fu_half_life(fu_material(nf%nvars(nf%n_vars)%species))) // ' sec'
            call put_chatt_nc(nf, nf%nvars(nf%n_vars)%varId, nf%nvars(nf%n_vars)%chVarNm, &
                            & "half_life", chAtt)
          endif
          !
          ! long_name
          !
          chAtt=trim(fu_quantity_string(outvar_quantity))//' '//trim(fu_str(nf%nvars(nf%n_vars)%species))

        else
          ! Just use default units
          nf%nvars(nf%n_vars)%unit = fu_quantity_unit(outvar_quantity)
          chAtt=trim(fu_quantity_string(outvar_quantity)) !long_name
        endif
        call put_chatt_nc(nf, nf%nvars(nf%n_vars)%varId, nf%nvars(nf%n_vars)%chVarNm, &
                        & "units", nf%nvars(nf%n_vars)%unit)
        call put_chatt_nc(nf, nf%nvars(nf%n_vars)%varId, nf%nvars(nf%n_vars)%chVarNm, &
                        & "long_name", chAtt)

        if (if_rotated_pole) then
          call put_chatt_nc(nf, nf%nvars(nf%n_vars)%varId, nf%nvars(nf%n_vars)%chVarNm, "grid_mapping", "rp")
        endif

        !We have no clue on standard names....
!        call put_chatt_nc(nf, nf%nvars(nf%n_vars)%varId, nf%nvars(nf%n_vars)%chVarNm, "standard_name", fu_quantity_standard_name(outvar_quantity))
      end do  ! iList item
    end do  ! iLst

    if (nf%n_vars == 0) then
      call set_error('No variables in output list','open_netcdf_file_o')
      return
    endif

    call grid_dimensions(myOutGrid,nx,ny)
    ! Close define mode
    if (nf%ncver==3 .and. nf%ifmpiio) then
#ifdef WITH_PNETCDF
      !!! NoFields2d
      nf%nf90mpiBUFfdlds = fu_get_default_mpi_buf_size()
      !nf%nf90mpiBUFfdlds = min(100, NoFields2d * buffer_timesteps)
      !      nf%nf90mpiNfieldsPerTimestep = NoFields2d

      nf%nf90mpifieldCnt = 0
      call msg("Allocating output buffer of  fields2D, fields per time step", nf%nf90mpiBUFfdlds, NoFields2d)
      allocate ( nf%nf90mpiReq(nf%nf90mpiBUFfdlds), &
              & nf%nf90mpiStat(nf%nf90mpiBUFfdlds), stat=iStat)
      if (fu_fails(iStat == 0, "Failed to allocate nf90MpiStuff", sub_name)) return
      call msg("nf90mpi_buffer_attach size, kB", nx*ny*(nf%nf90mpiBUFfdlds)*4/1024)
      iStat =nf90mpi_buffer_attach(nf%unit_bin, nx*ny*int((nf%nf90mpiBUFfdlds)*4, kind=8))
      if (fu_fails(iStat == 0, "nf90mpi_buffer_attach: "//nf90mpi_strerror(iStat), sub_name)) return 
      iStat =nf90mpi_enddef(nf%unit_bin, V_ALIGN=int(iTmp*4,kind=8), R_ALIGN=int(iTmp*4,kind=8))
      if (fu_fails(iStat == 0, "nf90mpi_enddef: "//nf90mpi_strerror(iStat), sub_name)) return 
#else             
      call set_error("Compiled without PNETCDF enddef", sub_name)
#endif
    else
       iTmp = chunks2d(1)*chunks2d(2)*4 ! 
       iStat = nf90_enddef(nf%unit_bin, V_ALIGN=iTmp*4, R_ALIGN=iTmp*4) 
       if (fu_fails(iStat == 0, "nf90_enddef: "//nf90_strerror(iStat), sub_name)) return 
    endif

    
    ! Write dimension variables (only one of the processes in case of MPIIO)
    
#ifdef WITH_PNETCDF
    if (nf%ncver==3 .and. nf%ifmpiio) then
        iStat =  nf90mpi_begin_indep_data(nf%unit_bin)
        if (fu_fails(iStat == 0, "nf90mpi_end_indep_data: "//nf90mpi_strerror(iStat), sub_name)) return
    endif
#endif
    
    if (smpi_global_rank == 0 .or. .not. ifMPIIO ) then 
      iTmp = maxval((/nf%ngrids(1)%nx,  nf%ngrids(1)%ny, 6*nf%n_levs+2/))
      varData => fu_work_array(iTmp)
      SELECT CASE(fu_gridtype(grid))
        CASE(lonlat)
          do iVal = 1, nf%ngrids(1)%nx
           varData(iVal) = nf%ngrids(1)%x_start + (iVal - 1) * nf%ngrids(1)%x_step
          end do
          call put_var_nc(nf, nf%ngrids(1)%lonVarId, varData, (/1/), (/nf%ngrids(1)%nx/))
          if (error) then
            call set_error('Failed to write lon values','open_netcdf_file_o')
            return
          endif
          do iVal = 1, nf%ngrids(1)%ny
           varData(iVal) = nf%ngrids(1)%y_start + (iVal - 1) * nf%ngrids(1)%y_step
          end do
          call put_var_nc(nf, nf%ngrids(1)%latVarId, varData,  (/1/),  (/nf%ngrids(1)%ny/))
          if (error) then 
            call set_error('Failed to write lat values','open_netcdf_file_o')
            return
          endif
        case(anygrid)
          fPtr => fu_geolats_fld(grid)
          call put_var_nc(nf, nf%ngrids(1)%latVarId, fPtr, &
              & (/1,1/), (/nf%ngrids(1)%nx, nf%ngrids(1)%ny/))
          if (error)then 
            call set_error('Failed to write 2d lat values','open_netcdf_file_o')
            return
          endif
          fPtr => fu_geolons_fld(grid)
          call put_var_nc(nf, nf%ngrids(1)%lonVarId, fPtr, &
                         &  (/1,1/), (/nf%ngrids(1)%nx, nf%ngrids(1)%ny/))
          if (error) then
            call set_error('Failed to write 2d lon values','open_netcdf_file_o')
            return
          endif
   
          fPtr => fu_cos_map_rot_fld(grid)
          call put_var_nc(nf, nf%ngrids(1)%cosMapRotVarId, fPtr, &
                         &  (/1,1/),  (/nf%ngrids(1)%nx, nf%ngrids(1)%ny/)) 
          if (error)then 
            call set_error('Failed to write cosine values','open_netcdf_file_o')
            return
          endif
          fPtr => fu_sin_map_rot_fld(grid)
          call put_var_nc(nf, nf%ngrids(1)%sinMapRotVarId, fPtr, &
                         & (/1,1/), (/nf%ngrids(1)%nx, nf%ngrids(1)%ny/)) 
          if (error)then 
            call set_error('Failed to write sine values','open_netcdf_file_o')
            return
          endif
          fPtr => fu_dy_fld_m(grid)
          call put_var_nc(nf, nf%ngrids(1)%dyMVarId, fPtr, &
                         & (/1,1/), (/nf%ngrids(1)%nx, nf%ngrids(1)%ny/)) 
          if (error)then 
            call set_error('Failed to write dy values','open_netcdf_file_o')
            return
          endif
          fPtr => fu_dx_fld_m(grid)
          call put_var_nc(nf, nf%ngrids(1)%dxMVarId, fPtr, &
                         & (/1,1/), (/nf%ngrids(1)%nx, nf%ngrids(1)%ny/)) 
          if (error)then 
            call set_error('Failed to write dx values','open_netcdf_file_o')
            return
          endif
        CASE DEFAULT
          CALL set_error('Unknown gridtype','open_netcdf_file_o')
          RETURN
      END SELECT
   
   
      iTmp = nf%n_levs

      if (if_dz_needed) then
         do jTmp=1,iTmp !Level loop
                varData(jTmp) = fu_layer_thickness_m(fu_level(vertical, jTmp))
         enddo
         call put_var_nc(nf, nf%nlevs%dzVarId,varData(1:iTmp), (/1/),(/iTmp/))
         if (error) then
           call set_error('Failed to write "dz" values','open_netcdf_file_o')
           return
         endif
      endif

      if (if_ab_needed) then
         fPtra=>varData(1:iTmp)
         fPtrb=>varData(iTmp+1:2*iTmp)
         fPtrda=>varData(2*iTmp+1:3*iTmp)
         fPtrdb=>varData(3*iTmp+1:4*iTmp)
         fPtrahalf=>varData(4*iTmp+1:5*iTmp+1)
         fPtrbhalf=>varData(5*iTmp+2:6*iTmp+2)
         call hybrid_coefs(vertical, a_full=fPtra, b_full=fPtrb, a_half=fPtrahalf, b_half = fPtrbhalf)

   
         call put_var_nc(nf, nf%nlevs%aVarId,fPtra(1:iTmp), (/1/),(/iTmp/))
         if (error) then
           call set_error('Failed to write "a" values','open_netcdf_file_o')
           return
         endif
         call put_var_nc(nf, nf%nlevs%bVarId,fPtrb(1:iTmp), (/1/),(/iTmp/))
         if (error) then
           call set_error('Failed to write "b" values','open_netcdf_file_o')
           return
         endif
         if (if_ab_half_needed)then
            fPtrda = fPtrahalf(2:itmp+1)- fPtrahalf(1:itmp)
            fPtrdb = fPtrbhalf(2:itmp+1)- fPtrbhalf(1:itmp)
            call put_var_nc(nf, nf%nlevs%ahalfVarId,fPtrahalf(1:iTmp+1), (/1/),(/iTmp+1/))
            if (error) then
              call set_error('Failed to write "a_half" values','open_netcdf_file_o')
              return
            endif
            call put_var_nc(nf, nf%nlevs%bhalfVarId,fPtrbhalf(1:iTmp+1), (/1/),(/iTmp+1/))
            if (error) then
              call set_error('Failed to write "b_half" values','open_netcdf_file_o')
              return
            endif
            call put_var_nc(nf, nf%nlevs%daVarId,fPtrda(1:iTmp), (/1/),(/iTmp/))
            if (error) then
              call set_error('Failed to write "da" values','open_netcdf_file_o')
              return
            endif
            call put_var_nc(nf, nf%nlevs%dbVarId,fPtrdb(1:iTmp), (/1/),(/iTmp/))
            if (error) then
              call set_error('Failed to write "db" values','open_netcdf_file_o')
              return
            endif
         endif
         do iVal = 1,iTmp !Approximate height of levels -- make GrADS happy
              nf%nlevs%std_z(iVal) = fu_height_for_press(fPtra(iVal) + fPtrb(iVal)*std_pressure_sl)
         enddo
      else  !Height levels
         ! std_z = z 
         nf%nlevs%std_z(1:iTmp) = nf%nlevs%levels(1:iTmp)
      endif
   
      !
      ! Put vertical dimension variable
      call put_var_nc(nf, nf%nlevs%levVarId, nf%nlevs%std_z(1:iTmp), (/1/),(/iTmp/))
      if (error) then
        call set_error('Failed to write lev values','open_netcdf_file_o')
        return
      endif
      call free_work_array(varData)

    endif !smpi_rank == 0
#ifdef WITH_PNETCDF
    if (nf%ncver==3 .and. nf%ifmpiio) then
        iStat =  nf90mpi_end_indep_data(nf%unit_bin)
        if (fu_fails(iStat == 0, "nf90mpi_end_indep_data: "//nf90mpi_strerror(iStat), sub_name)) return
    endif
#endif


    nf%defined = silja_true
    call destroy_namelist(nl)

    call msg(('Netcdf file open for output: "'//trim(nf%tmpfname))//'"')


  END FUNCTION open_netcdf_file_o  


  !********************************************
   
  integer function open_netcdf_file_i(fName, fFormat)
    !
    ! Opens the netcdf file for input; looks, what's inside; tries to fill the netcdf file structure:
    ! dimensions, variable names, etc 
    ! The function returns the place of the file in the nfile_ptr structure.
    !
    ! A bit more complicated than output - unknown number of  dimensions, each variable can come on it's
    ! own grid (staggered grid for winds, for instance). 
    ! Some dimensions can be used for string length. 
    ! 
    ! Variables and dimensions can either be specified in the name table file or the file should
    ! follow CF conventions. 
    !

    implicit none

    ! Imported parameters
    character(len=*), intent(in) :: fName 
    type(silam_fformat) :: fformat

    TYPE(netcdf_file),POINTER :: nf
    integer :: iStat, iFile, formatNum, nAtts, uLimDimId, iTmp, attLen, xType, iVar, iVal, iNl, &
             & nVals, idim, it, days_year, ilev, nDims, nMon, nDay, nHr, nMin, nxStag, nyStag, &
             & nxRef, nyRef, jTmp, nx, ny
    integer, dimension(:), pointer ::  iAtt
!    integer, dimension(4) :: dimIds
    character (len=nf90_max_name) :: aVar, bVar, aHalfVar, bHalfVar,  P0Var, PsVar, ChTmp2
    character (len=fnlen) :: chAtt, chTmp
    character (len=nf90_max_name) :: attName
    character(len=255) :: substance_name
    real :: P0, fSec, sPoleLat_tmp, sPoleLon_tmp
    character (len=nf90_max_name), dimension(:), allocatable :: lstVarName, lstDimName, &
                                       &lstDimVar,  gVarNmNf, gVarNm,  refVarNm
    character (len=clen), dimension(:), allocatable :: lstAxis, lstSilamQ, lstVType,  &
                                                   & lstDimType, lstSubst
    character, dimension(:), allocatable :: arrChTmp
    real, dimension(:), pointer :: fAtt, vals
    real, dimension(:), allocatable :: lstlevValue, lstMode, lstWavelen, lstFactor, lstOffset
    real(kind=8), dimension(:), pointer :: f8arrPtr
    real, dimension(:,:), pointer :: lons2Dptr, lats2Dptr
    integer, dimension(:), allocatable :: gVarId, refVarId
    integer, parameter :: nMaxAttVals = 15
    type(Tsilam_namelist), pointer :: nlPtr
    type(Tsilam_nl_item_ptr), dimension(:), pointer :: ptrItems
    character (len=fnlen) :: strContent
    logical :: ifFound, ifDImFound, ifNew          !, tFromFnm
    character (len=nf90_max_name), dimension(4) :: lstChTmp
    real ::  fTmp, x1, y1, x2, y2, dx, dy, wavelength, mean_diameter, bottom, thickness, wrfdxm, wrfdym
    real(r8k) :: jDate
    type(silja_level), dimension(:), pointer :: lstLevs
    type(silja_time) :: ttmp
    type(Tsilam_namelist_group), pointer::nlgrpptr
    type(Taerosol_mode) :: aerosol_mode
    type(silam_species) :: species
    type(Tsilam_namelist), pointer :: nl !!!For species handling
    logical :: ifSpeciesFromAtts, if2DLatLon

   character (len=*), parameter :: sub_name='open_netcdf_file_i'
    
    nl => fu_create_namelist()

    ! Defaults to be overriden with  the nametable if needed
    sPoleLat_tmp = -90.0
    sPoleLon_tmp = 0.0
    aVar = 'a'
    bVar = 'b'
    aHalfVar = 'a_half'
    bHalfVar = 'b_half'
    !
    ! Open the file, find structure for it
    !
    if(.not. associated(nmTblNlGrp))then
      call set_error('NETCDF name table not associated', sub_name)
    endif
    if(error)return
    nlgrpptr => nmTblNlGrp
    open_netcdf_file_i = -1
    
    ! Find a spare file counter. The last structure is always kept free!!
    iFile = fu_next_free_netcdf_structure()
    if(error)return

    open_netcdf_file_i = iFile
    nf => nfile(iFile)

    if (.not. associated(nf))then
      call set_error('nf not associated', sub_name)
    endif
    nf%defined = silja_undefined
    nf%nTime%first_valid_time = time_missing
    nf%nTime%last_valid_time = time_missing
    nf%nTime%tDimStart = time_missing
    nf%nTime%analysis_time   = time_missing  

    if(len_trim(fName) > 0)then
      nf%fname = fName
    else
      call set_error('Filename missing', sub_name)
      return
    endif

    
    istat = nf90_open(fname, NF90_NOWRITE, nf%unit_bin)
    if (istat /= 0)then 
      call msg_warning(fu_connect_strings('Netcdf error1:', nf90_strerror(iStat)),sub_name)
      call set_error(fu_connect_strings('Failed to open input netcdf file: ', fName),sub_name)
      return
    endif
    call msg('Opened NetCDF file for reading: '// fName,  iFile)

    istat = nf90_inquire(nf%unit_bin, nf%n_Dims, nf%n_Vars, nAtts, ulimDimId)  !, formatNum)
    if (istat /= 0)then 
      call msg_warning(fu_connect_strings('Netcdf error2:', nf90_strerror(iStat)),sub_name)
      call set_error(fu_connect_strings('Failed to inquire input netcdf file: ', fName),sub_name)
      return
    endif    

    allocate (nf%ndims(nf%n_Dims), nf%nvars(nf%n_Vars), stat=istat)

    if(ifMsgs .or. istat /= 0)then
      call msg('Number of dimensions in file:', nf%n_Dims)
      call msg('Number of variables:', nf%n_Vars)
      call msg('Number of global attributes:', nAtts)
      call msg('Unlimited dim id:', uLimDimId)
    endif
    if (istat /= 0)then 
         call msg('Failed to allocate NetCDF file structure No:', iFile)
         call set_error('failed to allocate nf%ndims(nf%n_Dims), nf%nvars(nf%n_Vars)',sub_name)
         return
    endif  

    
    ! Global attributes
    ! Report them, get file title
    ! The file will be recognized and name table chosen by file title (not nice, but don't see any other way)
    ! 
    fAtt => fu_work_array()
    iAtt => fu_work_int_array() 
    ! WRF thingies, multiplied with map factor map give real dx and dy in meters
    dx=1
    dy=1

    do iTmp = 1, nAtts
       iStat = nf90_inq_attname(nf%unit_bin, nf90_global, iTmp, attName)
       if (istat /= 0)then 
         call msg('Attr',iTmp)
         call msg_warning(fu_connect_strings('Netcdf error3:', nf90_strerror(iStat)),sub_name)
         call msg_warning('Failed to inquire global attribute name',sub_name)
       endif  
       iStat = nf90_inquire_attribute(nf%unit_bin, nf90_global, attName, xType, attLen)
       if (istat /= 0)then 
         call msg_warning(fu_connect_strings('Netcdf error4:', nf90_strerror(iStat)),sub_name)
         call msg_warning(fu_connect_strings('Failed to inquire global attribute: ', attName),sub_name)
       endif 

       select case(xType)
       case(nf90_char)
         if(attlen > len(chAtt))then
          chAtt = '===>> too long, sorry <<==='
          call msg('Global attribute ('//trim(fu_str(iTmp))//'): ' // trim(attname) // ' :: too long, sorry')
         else
           iStat = nf90_get_att(nf%unit_bin, nf90_global, attName, chAtt)
           if (istat /= 0)then 
             call msg_warning(fu_connect_strings('Netcdf error5:', nf90_strerror(iStat)),sub_name)
             call msg_warning(fu_connect_strings('Failed to get global attribute (',fu_str(iTmp),'): ', attName),sub_name)
           endif 
           if(ifMsgs)then
            call msg('Global attribute ('//trim(fu_str(iTmp))//') str: ' // trim(attname) // '  :: ' // trim(chAtt))
           endif

           select case (fu_str_u_case(trim(attName)))
             case('TITLE')
                nf%title = trim(chAtt)
              case ('SIMULATION_START_DATE') ! WRF thingies
                nf%nTime%analysis_time = fu_netcdf_str_to_silam_time(chAtt)
              case('_COORDINATEMODELRUNDATE')
                nf%nTime%analysis_time = fu_netcdf_str_to_silam_time(chAtt)
              case default
                if (ifMsgs) then
                  do jTmp = 1, attlen
                    if (ichar(chAtt(jTmp:jTmp)) == 0) exit !Null character 
                  enddo
                  call msg_warning('Unrecognized global string att (' // &
                      & trim(fu_str(iTmp)) //') "'// fu_str_u_case(trim(attName))//'" = '//chAtt(1:jTmp-1), sub_name)
                endif
            end select
         endif

       case(nf90_int) 
         iStat = nf90_get_att(nf%unit_bin, nf90_global, attName, iAtt)
         if (istat /= 0)then 
           call msg_warning(fu_connect_strings('Netcdf error6:', nf90_strerror(iStat)),sub_name)
           call msg_warning(fu_connect_strings('Failed to get global attribute (',fu_str(iTmp),'): ', attName),sub_name)
         endif

         if(ifMsgs) call msg(fu_connect_strings('Global attribute (',fu_str(iTmp),') int:', attname, '::'), iAtt(1:attLen))

       case(nf90_float, nf90_double)
         iStat = nf90_get_att(nf%unit_bin, nf90_global, attName, fAtt)
         if (istat /= 0)then 
           call msg_warning(fu_connect_strings('Netcdf error7:', nf90_strerror(iStat)),sub_name)
           call msg_warning(fu_connect_strings('Failed to get global attribute (',fu_str(iTmp),'): ', attName),sub_name)
         endif 
         if(ifMsgs) &
            &call msg(fu_connect_strings('Global attribute (',fu_str(iTmp),') float:', attname, '::'), fAtt(1:attlen))
         
         select case(fu_str_u_case(trim(attName)))
           case ('DX') ! WRF thingies
             wrfdxm = fAtt(1)
           case ('DY') ! WRF thingies
             wrfdym = fAtt(1)
           case ('POLE_LAT') !SILAM thingies
             sPoleLat_tmp = fAtt(1)
           case ('POLE_LON') !SILAM thingies
             sPoleLon_tmp = fAtt(1)
           case default
              if (ifMsgs) then
               call msg_warning('Unrecognized global float att ' // &
                  & trim(fu_str(iTmp)) // ' = '//trim(fu_str(fAtt(1))), sub_name)
              endif
          end select


         if(fu_str_u_case(trim(adjustl(attName))) == 'DY') dy = fAtt(1)

       case default
         call msg_warning('Attribute of unknown type',sub_name)
       end select
    enddo 


    ! Looking for an exact match for file title in name table namelist
    ! Note that this namelist group might not even exist is fiel is self-defining
    !
    ifFound = .false.
    nullify(ptrItems)
    if(len(trim(fformat%title)) > 0)then
      call msg(fu_connect_strings('File title taken from format string:', fformat%title))
      nf%title = fformat%title
    endif
    if(len(trim(nf%title)) < 1)then
      call set_error('No empty nametable titles allowed anymore',sub_name)
      return
    endif
    if(associated(nlgrpptr))then
      do iNl = 1, fu_nbr_of_namelists(nlgrpptr)
        nlPtr => fu_namelist(nlgrpptr, iNl)
        call get_items(nlPtr, 'title', ptrItems, nVals)
        if(nVals < 1)then
           call set_error('No titles in namelist',sub_name)
           return
        endif
        do iTmp = 1, nVals
            if(len(trim(fu_content(ptrItems(iTmp)))) < 1)cycle
            if(adjustl(trim(fu_content(ptrItems(iTmp)))) == adjustl(trim(nf%title)))then  
            !if((index(trim(adjustl(fu_content(ptrItems(iTmp)))),trim(adjustl(nf%title))) > 0) &
            ! & .or. (index(trim(adjustl(nf%title)),trim(adjustl(fu_content(ptrItems(iTmp))))) > 0))then  
              ifFound = .true.
              call msg('Table in standard_name_table: ', iNl)
              exit
            endif
        enddo
        if(ifFound)exit
      enddo  ! iNl
    endif

    if(ifFound)then !Proceed with the nametable
      chTmp = fu_str_u_case(trim(nf%title))
      if(chTmp == 'SILAM_OUTPUT_UG' .or. chTmp == 'SILAM_OUTPUT' .or. chTmp=='SILAM_EMIS')then
        call msg('Silam now reads its output without nametable')
        call msg('Please remove SILAM_OUTPUT nametable in your netcdf_nametable')
        call set_error('SILAM_OUTPUT nametable is obsolete', sub_name)
        return
      endif

      ! Dimensions. Can be empty e.g. for SILAM output - neither even has to exist.
      ! If the file is self-explantory, nothing is needed here.
      !
      call get_items(nlPtr, 'dim', ptrItems, nDims)
      if(nDims > 0)then
        allocate(lstDimName(nDims), lstAxis(nDims), lstDimVar(nDims), lstDimType(nDims), stat=istat)
        if (istat /= 0)then 
             call msg('Failed to allocate NetCDF file structure No:', iFile)
             call set_error('failed to allocate nf%ndims(nf%n_Dims), nf%nvars(nf%n_Vars)',sub_name)
             return
        endif  

        do iTmp = 1, nDims
           strContent = fu_content(ptrItems(iTmp))
           read(unit=strContent, iostat=iStat, fmt=*) lstDimName(iTmp), lstAxis(iTmp), lstDimVar(iTmp), lstDimType(iTmp)
        enddo
      endif

      do iDim = 1, nf%n_Dims
        iStat = nf90_inquire_dimension(nf%unit_bin, iDim, nf%nDims(iDim)%dimName, nf%nDims(iDim)%dimLen)
        nf%nDims(iDim)%dimId = iDim
        if (ifMsgs) call msg("Dimension in file "//trim(nf%nDims(iDim)%dimName), idim, nf%nDims(iDim)%dimLen)

        do iTmp = 1, nDims
          if(fu_str_u_case(adjustl(lstDimName(iTmp)))==fu_str_u_case(adjustl(nf%nDims(iDim)%dimName)))then
            nf%nDims(iDim)%axis = lstAxis(iTmp)
            nf%nDims(iDim)%varName = lstDimVar(iTmp)
            nf%nDims(iDim)%chType = lstDimType(iTmp)
            nf%nDims(iDim)%defined = .true.
          endif
        enddo
      enddo  !nf%n_Dims

      !
      ! Checking explicitly coefficients for vertical.
      !


      call get_items(nlPtr, 'a', ptrItems, nVals)
      if (nVals > 0)then
        do iTmp = 1, nVals
          strContent = fu_content(ptrItems(iTmp))
          jTmp = index(strContent, " ")
          do iDim = 1, nf%n_dims
            if (fu_str_u_case(adjustl(nf%nDims(iDim)%dimName)) == fu_str_u_case(strContent(1:jTmp-1)))then
             call split_string(strContent(jTmp+1:), ' ', fAtt, iVar)
             if (iVar == nf%nDims(iDim)%dimlen) then
                  allocate(nf%nDims(iDim)%a(iVar), stat=istat)
                  if (istat /= 0)then 
                       call msg('Failed to allocate NetCDF file structure No:', iFile)
                       call set_error('failed to allocate nf%ndims(nf%n_Dims), nf%nvars(nf%n_Vars)',sub_name)
                       return
                  endif  

                  nf%nDims(iDim)%a(1:iVar) = fAtt(1:iVar)
             else
               call msg("Dimension "//trim(nf%nDims(iDim)%dimName)//" has length ", nf%nDims(iDim)%dimlen)
               call msg("Got  coefficients", iVar)
               call set_error('Failed to read hybrid a from name table:', sub_name)
               call msg(strContent)
               return
             endif
             exit
            endif
          enddo
        enddo
      endif

      call get_items(nlPtr, 'b', ptrItems, nVals)
      if (nVals > 0)then
        do iTmp = 1, nVals
          strContent = fu_content(ptrItems(iTmp))
          jTmp = index(strContent, " ")
          do iDim = 1, nf%n_dims
            if (fu_str_u_case(adjustl(nf%nDims(iDim)%dimName)) == fu_str_u_case(strContent(1:jTmp-1)))then
             call split_string(strContent(jTmp+1:), ' ', fAtt, iVar)
             if (iVar == nf%nDims(iDim)%dimlen) then
                  allocate(nf%nDims(iDim)%b(iVar), stat=istat)
                  if (istat /= 0)then 
                       call msg('Failed to allocate NetCDF file structure No:', iFile)
                       call set_error('failed to allocate nf%ndims(nf%n_Dims), nf%nvars(nf%n_Vars)',sub_name)
                       return
                  endif  

                     nf%nDims(iDim)%b(1:iVar) = fAtt(1:iVar)
             else
               call msg("Dimension "//trim(nf%nDims(iDim)%dimName)//" has length ", nf%nDims(iDim)%dimlen)
               call msg("Got  coefficients", iVar)
               call set_error('Failed to read hybrid b from name table:', sub_name)
               call msg(strContent)
               return
             endif
             exit
            endif
          enddo
        enddo
      endif

      call get_items(nlPtr, 'a_half', ptrItems, nVals)
      if (nVals > 0)then
        do iTmp = 1, nVals
          strContent = fu_content(ptrItems(iTmp))
          jTmp = index(strContent, " ")
          do iDim = 1, nf%n_dims
            if (fu_str_u_case(adjustl(nf%nDims(iDim)%dimName)) == fu_str_u_case(strContent(1:jTmp-1)))then
             call split_string(strContent(jTmp+1:), ' ', fAtt, iVar)
             if (iVar == nf%nDims(iDim)%dimlen) then
                  allocate(nf%nDims(iDim)%a_half(iVar), stat=istat)
                  if (istat /= 0)then 
                       call msg('Failed to allocate NetCDF file structure No:', iFile)
                       call set_error('failed to allocate nf%ndims(nf%n_Dims), nf%nvars(nf%n_Vars)',sub_name)
                       return
                  endif  
                     nf%nDims(iDim)%a_half(1:iVar) = fAtt(1:iVar)
             else
               call msg("Dimension "//trim(nf%nDims(iDim)%dimName)//" has length ", nf%nDims(iDim)%dimlen)
               call msg("Got  coefficients", iVar)
               call set_error('Failed to read hybrid a_half from name table:', sub_name)
               call msg(strContent)
               return
             endif
             exit
            endif
          enddo
        enddo
      endif

      call get_items(nlPtr, 'b_half', ptrItems, nVals)
      if (nVals > 0)then
        do iTmp = 1, nVals
          strContent = fu_content(ptrItems(iTmp))
          jTmp = index(strContent, " ")
          do iDim = 1, nf%n_dims
            if (fu_str_u_case(adjustl(nf%nDims(iDim)%dimName)) == fu_str_u_case(strContent(1:jTmp-1)))then
             call split_string(strContent(jTmp+1:), ' ', fAtt, iVar)
             if (iVar == nf%nDims(iDim)%dimlen) then
                  allocate(nf%nDims(iDim)%a_half(iVar), stat=istat)
                  if (istat /= 0)then 
                       call msg('Failed to allocate NetCDF file structure No:', iFile)
                       call set_error('failed to allocate nf%ndims(nf%n_Dims), nf%nvars(nf%n_Vars)',sub_name)
                       return
                  endif  
                     nf%nDims(iDim)%b_half(1:iVar) = fAtt(1:iVar)
             else
               call msg("Dimension "//trim(nf%nDims(iDim)%dimName)//" has length ", nf%nDims(iDim)%dimlen)
               call msg("Got  coefficients", iVar)
               call set_error('Failed to read hybrid b_half from name table:', sub_name)
               call msg(strContent)
               return
             endif
             exit
            endif
          enddo
        enddo
      endif



      call get_items(nlPtr, 'P0', ptrItems, nVals)
      if (nVals > 0)then
        do iTmp = 1, nVals
          strContent = fu_content(ptrItems(iTmp))
          read(unit=strContent, iostat=iStat, fmt=*) chTmp, fTmp
          do iDim = 1, nf%n_dims
            if (fu_str_u_case(adjustl(nf%nDims(iDim)%dimName)) == fu_str_u_case(adjustl(chTmp)))then
              nf%nDims(iDim)%P0 = fTmp
              exit
            endif
          enddo
        enddo
      endif

      call get_items(nlPtr, 'aHalfVar', ptrItems, nVals)
      if (nVals > 0)then
        do iTmp = 1, nVals
          strContent = fu_content(ptrItems(iTmp))
          read(unit=strContent, iostat=iStat, fmt=*) chTmp, chTmp2
          do iDim = 1, nf%n_dims
            if (fu_str_u_case(adjustl(nf%nDims(iDim)%dimName)) == fu_str_u_case(adjustl(chTmp)))then
              nf%nDims(iDim)%a_half_Var = chTmp2
              exit
            endif
          enddo
        enddo
      endif
      call get_items(nlPtr, 'aVar', ptrItems, nVals)
      if (nVals > 0)then
        do iTmp = 1, nVals
          strContent = fu_content(ptrItems(iTmp))
          read(unit=strContent, iostat=iStat, fmt=*) chTmp, chTmp2
          do iDim = 1, nf%n_dims
            if (fu_str_u_case(adjustl(nf%nDims(iDim)%dimName)) == fu_str_u_case(adjustl(chTmp)))then
              nf%nDims(iDim)%a_Var = chTmp2
              exit
            endif
          enddo
        enddo
      endif
      call get_items(nlPtr, 'bHalfVar', ptrItems, nVals)
      if (nVals > 0)then
        do iTmp = 1, nVals
          strContent = fu_content(ptrItems(iTmp))
          read(unit=strContent, iostat=iStat, fmt=*) chTmp, chTmp2
          do iDim = 1, nf%n_dims
            if (fu_str_u_case(adjustl(nf%nDims(iDim)%dimName)) == fu_str_u_case(adjustl(chTmp)))then
              nf%nDims(iDim)%b_half_Var = chTmp2
              exit
            endif
          enddo
        enddo
      endif
      call get_items(nlPtr, 'bVar', ptrItems, nVals)
      if (nVals > 0)then
        do iTmp = 1, nVals
          strContent = fu_content(ptrItems(iTmp))
          read(unit=strContent, iostat=iStat, fmt=*) chTmp, chTmp2
          do iDim = 1, nf%n_dims
            if (fu_str_u_case(adjustl(nf%nDims(iDim)%dimName)) == fu_str_u_case(adjustl(chTmp)))then
              nf%nDims(iDim)%b_Var = chTmp2
              exit
            endif
          enddo
        enddo
      endif
      call get_items(nlPtr, 'P0Var', ptrItems, nVals)
      if (nVals > 0)then
        do iTmp = 1, nVals
          strContent = fu_content(ptrItems(iTmp))
          read(unit=strContent, iostat=iStat, fmt=*) chTmp, chTmp2
          do iDim = 1, nf%n_dims
            if (fu_str_u_case(adjustl(nf%nDims(iDim)%dimName)) == fu_str_u_case(adjustl(chTmp)))then
              nf%nDims(iDim)%P0_Var = chTmp2
              exit
            endif
          enddo
        enddo
      endif
      call get_items(nlPtr, 'PsVar', ptrItems, nVals)
      if (nVals > 0)then
        do iTmp = 1, nVals
          strContent = fu_content(ptrItems(iTmp))
          read(unit=strContent, iostat=iStat, fmt=*) chTmp, chTmp2
          do iDim = 1, nf%n_dims
            if (fu_str_u_case(adjustl(nf%nDims(iDim)%dimName)) == fu_str_u_case(adjustl(chTmp)))then
              nf%nDims(iDim)%Ps_Var = chTmp2
              exit
            endif
          enddo
        enddo
      endif


      call get_items(nlPtr, 'gridvar', ptrItems, nVals)
      if(nVals > 0)then
        allocate(gVarNmNf(nVals)) 
        allocate(gVarNm(nVals))
        allocate(gVarId(nVals))
        allocate(refVarNm(nVals))
        allocate(refVarId(nVals))
        do iTmp = 1, nVals
          strContent = fu_content(ptrItems(iTmp))
          read(unit=strContent, iostat=iStat, fmt=*) gVarNmNf(iTmp), gVarNm(iTmp), refVarNm(iTmp)
          gVarId(iTmp) = int_missing
          refVarId(iTmp) = int_missing
        enddo
      endif

      ! If the map projection is rotated lon-lat, the rotation has to be defined in the nametable
      call get_items(nlPtr, 'south_pole_lat', ptrItems, nVals)
      if(nVals > 0)then
          sPoleLat_tmp  = fu_content_real(ptrItems(1))
          call get_items(nlPtr, 'south_pole_lon', ptrItems, nVals)
          if(nVals > 0)then
             sPoleLon_tmp  = fu_content_real(ptrItems(1))
          else
             call set_error('Pole lat given in nc name table but lon missing',sub_name)
          endif
      endif

      ! If analysis time is the start of the time dimension or the first tiem in file
      call get_items(nlPtr, 'analysis_time', ptrItems, nVals)
      if(nVals > 1)then
        call set_error('Too many analysis_time lines in name table', sub_name)
        return
      elseif(nVals == 1)then
        if(fu_str_u_case(adjustl(trim(fu_content(ptrItems(1))))) == 'TIME_DIMENSION_START')then
          nf%nTime%analysis_t_from = dim_start
        elseif(fu_str_u_case(adjustl(trim(fu_content(ptrItems(1))))) == 'FIRST_TIME_IN_FILE')then
          nf%nTime%analysis_t_from = first_value
        else
          nf%nTime%analysis_time = fu_io_string_to_time(fu_content(ptrItems(1)))
          if(error)return
          nf%nTime%analysis_t_from = int_missing
        endif
      else
        if(ifMsgs) call msg('No analysis_time line in name table - taking the first time in the file')
        nf%nTime%analysis_t_from = first_value
      endif

      ! Variables
      ! Either species or cocktails ..
      call get_items(nlPtr, 'ifCocktails', ptrItems, nVals)
      if(nVals > 1)then
        call set_error('Too many ifCocktails lines in name table', sub_name)
        return
      elseif(nVals == 1)then
        if(fu_str_u_case(adjustl(trim(fu_content(ptrItems(1))))) == 'YES')then
          nf%ifCocktails = .true.
          if(ifMsgs) call msg('File includes cocktails')
        elseif(fu_str_u_case(adjustl(trim(fu_content(ptrItems(1))))) == 'NO')then
          nf%ifCocktails = .false.
          if(ifMsgs) call msg('File includes species')
        else
            call set_error('Strange ifCocktails in nameTable:'+fu_content(ptrItems(1)),sub_name)
            return
        endif
      else
        nf%ifCocktails = .false.
        if(ifMsgs) call msg('No ifCocktails line in name table - assuming species')
      endif

        !
        ! For averaged variables, the time label can be at the 
        ! beginning, middle and end of the period. No info of this kind in NetCDF, have to
        ! use name table
        !
        if(index(fu_str_l_case(fu_content(nlPtr, 'time_label_position')), 'start_of_period') > 0)then
          nf%nTime%time_label_position = start_of_period
        elseif(index(fu_str_l_case(fu_content(nlPtr, 'time_label_position')), 'mid_of_period') > 0)then
          nf%nTime%time_label_position = mid_of_period
        elseif(index(fu_str_l_case(fu_content(nlPtr, 'time_label_position')), 'end_of_period') > 0)then
          nf%nTime%time_label_position = end_of_period
        elseif(index(fu_str_l_case(fu_content(nlPtr, 'time_label_position')), 'instant') > 0)then
          nf%nTime%time_label_position = instant_fields
        elseif(fu_content(nlPtr, 'time_label_position') == '')then
          call msg_warning('No time_label_position in netcdf name table - assuming instant', sub_name)
          nf%nTime%time_label_position = instant_fields
        else
          call set_error('strange time_label_position in superctl:' + &
                       & fu_content(nlPtr, 'time_label_position'), sub_name)
          return
        endif
      

      call get_items(nlPtr, 'var', ptrItems, nVals)
      if(nVals > 0)then
        allocate(lstVarName(nVals)) 
        allocate(lstSilamQ(nVals))
        allocate(lstVType(nVals))
        allocate(lstlevValue(nVals))
        allocate(lstSubst(nVals))
        allocate(lstMode(nVals))
        allocate(lstWavelen(nVals))
        allocate(lstFactor(nVals))
        allocate(lstOffset(nVals))
        do iTmp = 1, nVals
          strContent = fu_content(ptrItems(iTmp))
          read(unit=strContent, iostat=iStat, fmt=*) lstVarName(iTmp), lstSilamQ(iTmp), &
                                       & lstVType(iTmp), lstlevValue(iTmp), lstSubst(iTmp), &
                                       & lstMode(iTmp), lstWavelen(iTmp),  lstFactor(iTmp), lstOffset(iTmp) 
        enddo

      endif ! if nVals > 0

    else

      chTmp = fu_str_u_case(trim(nf%title))
      if(chTmp == 'SILAM_OUTPUT' .or. chTmp=='SILAM_EMIS')then
        nf%ntime%time_label_position = end_of_period  !!!Silam averages this way
      else
        nf%nTime%time_label_position = instant_fields  !! Others must be instant
        call msg_warning(fu_connect_strings('File not in nametable: ', nf%title), sub_name)
        call msg('Proceed without nametable')
      endif

      do iDim = 1, nf%n_Dims
        iStat = nf90_inquire_dimension(nf%unit_bin, iDim, nf%nDims(iDim)%dimName, nf%nDims(iDim)%dimLen)
        if (fu_fails (iStat == 0 , "CF nf90_inquire_dimension :"//nf90_strerror(iStat), sub_name)) return
        nf%nDims(iDim)%dimId = iDim
        nf%nDims(iDim)%defined = .false.
      enddo  !nf%n_Dims

    endif  ! if associated namelist group from NetCDF code table 

    !-------------------------------------------------------------------
    !
    ! Main cycle getting variable names from the file and their attributes
    !
    do iVar = 1, nf%n_Vars
      nf%nVars(iVar)%varid = iVar
      iStat = nf90_inquire_variable(nf%unit_bin, iVar, nf%nVars(iVar)%chVarNm, nf%nVars(iVar)%xtype, &
                                  & nf%nVars(iVar)%n_dims, nf%nVars(iVar)%dimids, nAtts)
      if(iStat /= 0)then
        call msg_warning(fu_connect_strings('Cannot get the variable,',nf90_strerror(iStat)))
        cycle
      endif
      
      if(allocated(gVarNmNf))then
        do iTmp = 1, size(gVarNmNf)
          if(fu_str_u_case(adjustl(trim(nf%nVars(iVar)%chVarNm))) == fu_str_u_case(adjustl(trim(gVarNmNf(iTmp)))))then
            gVarId(iTmp) = iVar
          endif
          if(fu_str_u_case(adjustl(trim(nf%nVars(iVar)%chVarNm))) == fu_str_u_case(adjustl(trim(refVarNm(iTmp)))))then
            refVarId(iTmp) = iVar
          endif
        enddo
      endif

      ! First the simple case of SILAM output. Needs to know, if 2d or 3d variable.
      !
      chTmp = fu_str_u_case(trim(nf%title))

      ifSpeciesFromAtts = .False.

      if( chTmp == 'SILAM_OUTPUT_UG' .or. chTmp == 'SILAM_OUTPUT' .or. chTmp=='SILAM_EMIS')then
        !
        nf%nVars(iVar)%if3D = .false.
        nf%nVars(iVar)%n_levs = int_missing
        do iTmp = 1, nf%nVars(iVar)%n_dims
           !Now only 3D vars have z dimension
           if ( nf%nDims(nf%nVars(iVar)%dimids(iTmp))%axis == 'z')then
             nf%nVars(iVar)%if3D = .true.
             exit
           endif
        enddo

        call decode_id_params_from_io_str(nf%nVars(iVar)%chVarNm, &
                                          & nf%nVars(iVar)%if3D, &
                                          & nf%nVars(iVar)%quantity, &
                                          & species, &
                                          & .false.)  ! can be non-SILAM quantity

        ifSpeciesFromAtts = defined(species)


        if(error)return

        if(ifMsgs) call msg("Variable "//trim(nf%nVars(iVar)%chVarNm)//" "//fu_quantity_short_string(nf%nVars(iVar)%quantity)+':'+fu_str(species))
      else
        !
        ! Not SILAM output
        ! Compare the variable name to the name table and find the SILAM quantity
        !
        if (allocated(lstVarName)) then
          do iTmp = 1, size(lstVarName)
            if(trim(nf%nVars(iVar)%chVarNm) == trim(lstVarName(iTmp)))then
               nf%nVars(iVar)%quantity = fu_get_silam_quantity(lstSilamQ(iTmp))
               if(error)then
                 call set_error("Failed to parse variable listed in the nametable",sub_name)
                 return
               endif

               if(ifMsgs) call msg(fu_connect_strings('Found variable:', &
                                            & nf%nVars(iVar)%chVarNm,  ', silam quantity: '), &
                          & nf%nVars(iVar)%quantity)


               jTmp = fu_str2leveltype(lstVType(iTmp))
               if(jTmp /= any_level .and. jTmp /= no_level) then
                 nf%nVars(iVar)%n_levs = 1
                 if (jTmp == int_missing) then
                    call set_error('Strange type of level:'//trim(lstVType(iTmp)),  sub_name)
                    return
                 elseif (jTmp == layer_btw_2_height) then ! 'HEIGHT_LYR_FROM_SURF'
                   call set_vertical(&
                       & fu_set_layer_between_two(layer_btw_2_height, lstlevValue(iTmp)*2.0, 0.0), &
                       & nf%nVars(iVar)%sVert)
                 else
                    call set_vertical(fu_set_level(jTmp,lstlevValue(iTmp)), nf%nVars(iVar)%sVert)
                 endif
                    
                  if(error)return
               endif

               substance_name = ''
               mean_diameter = real_missing
               wavelength = real_missing
              if(lstSubst(iTmp) /= 'XXX')substance_name = lstSubst(iTmp)
              if(lstMode(iTmp) /= -1.)mean_diameter = lstMode(iTmp)
              if(lstWavelen(iTmp) /= -1.) wavelength = lstWavelen(iTmp)
               if ((mean_diameter > 0.0) .and. (mean_diameter < 1.0)) then
                 aerosol_mode = fu_set_mode(fixed_diameter_flag, &
                                          & mean_diameter*0.99, mean_diameter*1.01, mean_diameter)               
               else
              aerosol_mode = in_gas_phase
               end if
               if (substance_name /= '') then
                if(nf%ifCocktails)then
                 species = species_missing
                  nf%nVars(iVar)%chCocktailNm = substance_name
                else
                  call set_species(species, fu_get_material_ptr(substance_name), aerosol_mode, wavelength)
                endif
              else
                species = species_missing
              endif
               if (error) return
               if (ifMsgs) call msg("Substance name:"+ substance_name)
               exit
             endif
          enddo
        endif !associated lstVarName
      endif !Mot a silam output


      ! Variable attributes
      ! Have to know the attribute names we're looking for
      !
!      call msg("Variabe attributes, quantity:"+ fu_quantity_short_string( nf%nVars(iVar)%quantity))

      nf%nVars(iVar)%offset_nc = 0.
      nf%nVars(iVar)%ScaleFactor_nc= 1.
      nf%nVars(iVar)%offset = 0.
      nf%nVars(iVar)%ScaleFactor= 1.
      do iTmp = 1, nAtts
        iStat = nf90_inq_attname(nf%unit_bin, iVar, iTmp, attName)
        if (istat /= 0)then 
          call msg_warning(fu_connect_strings('Netcdf error8:', nf90_strerror(iStat)),sub_name)
          call msg_warning('Failed to inquire attribute name',sub_name)
        else  

         ! Report all the attributes anyway
          iStat = nf90_inquire_attribute(nf%unit_bin, iVar, attName, xType, attLen)
          if (istat /= 0)then 
            call msg_warning(fu_connect_strings('Netcdf error9:', nf90_strerror(iStat)),sub_name)
            call msg_warning(fu_connect_strings('Failed to inquire attribute: ', attName),sub_name)
          else
 
            select case(xType)
            case(nf90_char)
              iStat = nf90_get_att(nf%unit_bin, iVar, attName, chAtt)
              if (istat /= 0)then 
                call msg_warning(fu_connect_strings('Netcdf error10:', nf90_strerror(iStat)),sub_name)
                call msg_warning(fu_connect_strings('Failed to get attribute: ', attName),sub_name)
              endif 
              if(ifMsgs)then
                call msg(fu_connect_strings('Attribute:', attname, '::' , chAtt))
              endif

              if (ifSpeciesFromAtts) call add_namelist_item(nl,attName, chAtt) !!Save attribute as a namelist
              
            case(nf90_int, NF90_USHORT,  NF90_SHORT, NF90_UINT, NF90_INT64, NF90_UINT64) 
              iStat = nf90_get_att(nf%unit_bin, iVar, attName, iAtt)
              if (istat /= 0)then 
                call msg_warning(fu_connect_strings('Netcdf error11:', nf90_strerror(iStat)),sub_name)
                call msg_warning(fu_connect_strings('Failed to get global attribute: ', attName),sub_name)
              endif
              do iVal = 1, attLen 
                if(ifMsgs) call msg(fu_connect_strings('Attribute:', attname, '::'), int_value = iAtt(iVal))
                fAtt(iVal) =  real(iAtt(iVal)) ! Make float, just in case
              enddo
            case(nf90_float, nf90_double)
              iStat = nf90_get_att(nf%unit_bin, iVar, attName, fAtt)
              if (istat /= 0)then 
                call msg_warning(fu_connect_strings('Netcdf error12:', nf90_strerror(iStat)),sub_name)
                call msg_warning(fu_connect_strings('Failed to get attribute: ', attName),sub_name)
              endif 
              if(ifMsgs)then
                do iVal = 1, attlen
                  call msg(fu_connect_strings('Attribute:', attname, '::'), fAtt(iVal))
                enddo
              endif
            case default
              call msg_warning('Attribute of unknown type',sub_name)
            end select

       ! And now check, if we know, what to do with it. For now just put the values of those to the variable structure 
       !
            select case(attname)          
            case('units')
               nf%nVars(iVar)%unit = chAtt
            case('add_offset')
               nf%nVars(iVar)%offset_nc = fAtt(1)
            case('scale_factor')
               nf%nVars(iVar)%scaleFactor_nc =  fAtt(1)
            case('silam_amount_unit')
               nf%nVars(iVar)%silam_amount_unit =  chAtt
            case('molar_mass')
              if(fu_str_u_case(trim(nf%title)) == 'SILAM_OUTPUT')then
                 fTmp = fu_set_named_value(chAtt) ! kg/mole
                 ! Silam tends to put some rubbish as molar mass sometimes
                 if (fTmp > 0 .and. fTmp < 100) then 
                   nf%nVars(iVar)%nc_molar_mass =  fTmp
                 endif
              endif
            case('_CoordinateAxisType') !Fall-back option
               if (nf%nVars(iVar)%axis == '') nf%nVars(iVar)%axis = chAtt(1:1)         
            case('axis')
               nf%nVars(iVar)%axis = chAtt         
            case('positive')
               nf%nVars(iVar)%positive_direction = chAtt
            case('calendar')
               nf%nVars(iVar)%calendar = chAtt    
            case('missing_value', '_FillValue')
               if(any (xType == (/nf90_int, NF90_USHORT, NF90_SHORT, NF90_UINT, NF90_INT64, &
                                & NF90_UINT64, nf90_float, nf90_double, NF90_BYTE/)))then
                 nf%nVars(iVar)%missing_value = fAtt(1) !fatt must already contain converted value
               else
                 call set_error('Strange data type for missing value attribute:'+fu_str(fatt(1)), sub_name)
               endif                
            case('valid_min')
               if(xType == nf90_float .or. xType == nf90_double)then
                 nf%nVars(iVar)%valid_min = fAtt(1)
               elseif(xType == nf90_int)then
                 nf%nVars(iVar)%valid_min = real(iAtt(1)) 
               else
                 call msg_warning('Strange data type for valid_min attribute')
               endif
            case('valid_max')
               if(xType == nf90_float .or. xType == nf90_double)then
                 nf%nVars(iVar)%valid_max = fAtt(1)
               elseif(xType == nf90_int)then
                 nf%nVars(iVar)%valid_max = real(iAtt(1)) 
               else
                 call msg_warning('Strange data type for valid_max attribute')
               endif
            case('valid_range')
               if(xType == nf90_float .or. xType == nf90_double)then
                 nf%nVars(iVar)%valid_min = fAtt(1)
                 nf%nVars(iVar)%valid_max = fAtt(2)
               elseif(xType == nf90_int)then
                 nf%nVars(iVar)%valid_min = real(iAtt(1))
                 nf%nVars(iVar)%valid_max = real(iAtt(2)) 
               else
                 call msg_warning('Strange data type for valid_range attribute')
               endif
            case('A_var')
               aVar = chAtt
            case('B_var')
               bVar = chAtt
            case('A_half_var')
               aHalfVar = chAtt
            case('B_half_var')
               bHalfVar = chAtt
            case('P0_var')
               P0Var = chAtt
            case('PS_var')
               PsVar = chAtt
            case('standard_name')
               nf%nVars(iVar)%standard_name = chAtt
            !   if(nf%nVars(iVar)%quantity == int_missing)then
            !     nf%nVars(iVar)%quantity = fu_silam_quantity_from_std_name(chAtt)
            !     if(error)then
            !        call unset_error(sub_name)
            !        nf%nVars(iVar)%quantity = int_missing
            !     endif
            !   endif
            case('long_name', 'description')
               nf%nVars(iVar)%long_name = chAtt
            case('coordinates')
               nf%nVars(iVar)%coordinates = chAtt
            case('stagger')
               nf%nVars(iVar)%stagger = chAtt
            case('bounds')
               nf%nVars(iVar)%bounds = chAtt
            case('compress')
              call msg_warning('Variable with compressed dimensions', sub_name)
               nf%nVars(iVar)%compress = chAtt
               nf%nVars(iVar)%quantity = int_missing   ! Not usable
            case default
            end select
          endif
        endif
      enddo !attribute

      if(nf%nVars(iVar)%quantity == int_missing)then
        if(ifMsgs)then
          call msg_warning(fu_connect_strings('Not a SILAM quantity:',nf%nVars(iVar)%chVarNm))
        endif
      endif

      
      !!! Final setup of species and conversions
      if (ifSpeciesFromAtts) then 
        !!Restore needed parameters from attributes
        substance_name =  fu_content(nl, "substance_name")
        if (len_trim(substance_name) > 0) then
           aerosol_mode = fu_set_mode(nl)
           chAtt = fu_content(nl, "optical_wavelength")

           wavelength = real_missing
           if (len_trim(chatt) > 0) wavelength = fu_set_named_value(fu_content(nlPtr,'optical_wavelength'),.true.)

           call set_species(nf%nVars(iVar)%species, fu_get_material_ptr(substance_name), aerosol_mode, wavelength)
           if (error) then
              call set_error("trouble with nc variable at SpeciesFromAtts:"+trim(nf%nVars(iVar)%chVarNm),sub_name)
              return
           endif
             
           !Silam file might have been converted to ug/m3 for thredds
           if (nf%nVars(iVar)%unit(1:2) == 'ug') then
             chAtt=fu_content(nl, "silam_amount_unit")
             if (chatt == 'kg') then
               nf%nVars(iVar)%ScaleFactor_nc =  nf%nVars(iVar)%ScaleFactor_nc * 1e-9
               nf%nVars(iVar)%unit = 'kg'//nf%nVars(iVar)%unit(3:)
             else if (chAtt == 'mole') then
               fTmp = fu_set_named_value(fu_content(nl,'molar_mass'),.true.)
               nf%nVars(iVar)%ScaleFactor_nc =  nf%nVars(iVar)%ScaleFactor_nc * 1e-9/fTmp
               nf%nVars(iVar)%unit = 'mole'//nf%nVars(iVar)%unit(5:)
             else
               call set_error("Got strange silam_amount_unit for silam input in ug..", sub_name)
               return
             endif
           endif
        else
           call msg_warning("substance_name attribute missing, at SpeciesFromAtts:"+trim(nf%nVars(iVar)%chVarNm),sub_name)
           call msg("Resetting species!")
           nf%nVars(iVar)%species = species_missing
        endif
        call reset_namelist(nl)

      else
        nf%nVars(iVar)%species = species  
      endif

      ! Dimension var-s. Now try CF: same name as dim, right units, t - unlimited dim, 
      ! standard name attribute for z to get the vertical type
      do iDim = 1, nf%n_Dims
        if( nf%nDims(iDim)%defined)then
          if (ifMsgs) call msg("Dimension already defined", iDim)
        else
!          call msg("Trying CF for dimension " // trim(nf%nDims(iDim)%dimName) // ':', iDim )
!          call msg('"'//fu_str_u_case(adjustl(nf%nVars(iVar)%chVarNm)) //'" "' //fu_str_u_case(adjustl(nf%nDims(iDim)%dimName))//'"')
          if(fu_str_u_case(adjustl(nf%nVars(iVar)%chVarNm)) == fu_str_u_case(adjustl(nf%nDims(iDim)%dimName)))then
!            call msg("Got dimension var")
            nf%nDims(iDim)%varName = nf%nVars(iVar)%chVarNm
            nf%nDims(iDim)%varId = iVar
            ! Axis
            if(iDim == uLimDimId .or. (fu_str_l_case(nf%nVars(iVar)%axis) =='t') )then
              nf%nDims(iDim)%axis = 't'
              nf%nDims(iDim)%defined = .true.

            elseif(fu_str_l_case(adjustl(nf%nVars(iVar)%unit)) == 'degrees_north' .or. &
                 & fu_str_l_case(adjustl(nf%nVars(iVar)%unit)) == 'degree_north' .or. &
                 & fu_str_l_case(adjustl(nf%nVars(iVar)%unit)) == 'degrees_n' .or. &
                 & fu_str_l_case(adjustl(nf%nVars(iVar)%unit)) == 'degree_n' .or. &
                 & fu_str_l_case(adjustl(nf%nVars(iVar)%unit)) == 'degreesn' .or. &
                 & fu_str_l_case(adjustl(nf%nVars(iVar)%unit)) == 'degreen')then
                nf%nDims(iDim)%axis = 'y'
                nf%nDims(iDim)%defined = .true.
 
            elseif(fu_str_l_case(adjustl(nf%nVars(iVar)%unit)) == 'degrees_east' .or. &
                 & fu_str_l_case(adjustl(nf%nVars(iVar)%unit)) == 'degree_east' .or. &
                 & fu_str_l_case(adjustl(nf%nVars(iVar)%unit)) == 'degrees_e' .or. &
                 & fu_str_l_case(adjustl(nf%nVars(iVar)%unit)) == 'degree_e' .or. &
                 & fu_str_l_case(adjustl(nf%nVars(iVar)%unit)) == 'degreese' .or. &
                 & fu_str_l_case(adjustl(nf%nVars(iVar)%unit)) == 'decleneee')then
                nf%nDims(iDim)%axis = 'x'
                nf%nDims(iDim)%defined = .true.
 
            elseif(nf%nVars(iVar)%standard_name /= '')then
               ! Defaults for the case below  
              nf%nDims(iDim)%axis = 'z'
              nf%nDims(iDim)%defined = .true.                
              select case(fu_str_l_case(adjustl(nf%nVars(iVar)%standard_name)))
                case('layer_midpoint_height_above_ground') !Silam output goes 
                   nf%nDims(iDim)%chType = 'HEIGHT_LYR_FROM_SURF'
                case('height') !Silam output goes
                      nf%nDims(iDim)%chType = 'HEIGHT_FROM_SURF'
                 case('altitude')   
                   nf%nDims(iDim)%chType = 'ALTITUDE_FROM_SEA_LEVEL'
                 case('air_pressure')
                   nf%nDims(iDim)%chType = 'PRESSURE'
                 case('atmosphere_sigma_coordinate')
                   nf%nDims(iDim)%chType = 'sigma_level'
                 case('hybrid_layer_midpoint_standard_height_above_ground')
                   nf%nDims(iDim)%chType = 'HYBRID_LYR'
                   nf%nDims(iDim)%a_half_Var =  aHalfVar
                   nf%nDims(iDim)%b_half_Var =  bHalfVar
                 case('atmosphere_hybrid_sigma_pressure_coordinate')
                   nf%nDims(iDim)%chType = 'HYBRID'
                   nf%nDims(iDim)%a_Var = aVar
                   nf%nDims(iDim)%b_Var = bVar
                 case('surface')
                   nf%nDims(iDim)%chType = 'SURFACE_LEVEL'
                 case default
                   call msg_warning(fu_connect_strings('Strange standard name for vertical axes: ', nf%nVars(iVar)%standard_name),sub_name)
                   nf%nDims(iDim)%axis = ''
                   nf%nDims(iDim)%defined = .false.             
              end select
              if (ifMsgs) call msg("Z-dimension type (from standard_name): "//trim(nf%nDims(iDim)%chType))
            elseif(nf%nVars(iVar)%long_name /= '')then  !CDO kills standard_name
               ! Defaults for the case below  
              nf%nDims(iDim)%axis = 'z'
              nf%nDims(iDim)%defined = .true.                
              select case(fu_str_l_case(adjustl(nf%nVars(iVar)%long_name)))
                 case('layer midpoint constant height from surface') !Silam output goes 
                   nf%nDims(iDim)%chType = 'HEIGHT_LYR_FROM_SURF'
                 case('height above surface') !OLD Silam output...
                   nf%nDims(iDim)%chType = 'HEIGHT_LYR_FROM_SURF'
                 case('standard height of hybrid layer midpoint')
                   nf%nDims(iDim)%chType = 'HYBRID_LYR'
                   nf%nDims(iDim)%a_half_Var =  aHalfVar
                   nf%nDims(iDim)%b_half_Var =  bHalfVar
                 case default
                   call msg_warning('Strange long name for vertical axes: '//trim(nf%nVars(iVar)%long_name),sub_name)
                   call msg("Faild to recognize dimension from long_name"// trim(nf%nDims(iDim)%dimName) // ':', iDim )

                   call set_error("Failed with vertical", sub_name)
                   return
                  ! nf%nDims(iDim)%axis = ''
                  ! nf%nDims(iDim)%defined = .false.             
              end select
              if (ifMsgs) call msg("Z-dimension type (from long_name): "//trim(nf%nDims(iDim)%chType))
            else
              call msg("Faild to recognize dimension"// trim(nf%nDims(iDim)%dimName) // ':', iDim )
              call set_error("Failed", sub_name)
            endif
          endif


          if(nf%nDims(iDim)%defined) exit   ! go for the next variable

        endif  ! if nf%nDims(iDim)%defined

      enddo  ! cycle over dimensions

! FIXME Seems to be not needed at all....      
!      ! z dim variable with coefficients as attributes
!      !
!      iTmp=iDim  ! Save coeffs
!      if (iDim< nf%n_Dims) then ! The variable is dimension that jut had beed defined
!         do iDim = 1, nf%n_Dims
!           if (ifMsgs) call msg(trim(nf%nDims(iDim)%varName)//" dimension type: "//trim(nf%nDims(iDim)%chType))
!           if(nf%nVars(iVar)%chVarNm == nf%nDims(iTmp)%varName)then
!             if(nf%nDims(iTmp)%axis == 'z' .and. nf%nDims(iTmp)%chType == 'HYBRID')then
!                if(nf%nDims(iDim)%a_Var == '')nf%nDims(iDim)%a_Var =  aVar
!                if(nf%nDims(iDim)%b_Var == '')nf%nDims(iDim)%b_Var =  bVar
!                if(nf%nDims(iDim)%a_Half_Var == '')nf%nDims(iDim)%a_half_Var =  aHalfVar
!                if(nf%nDims(iDim)%b_Half_Var == '')nf%nDims(iDim)%b_half_Var =  bHalfVar
!                if(nf%nDims(iDim)%P0_Var == '')nf%nDims(iDim)%P0_Var =  P0Var
!                if(nf%nDims(iDim)%Ps_Var == '')nf%nDims(iDim)%Ps_Var = PsVar
!             endif
!           endif  
!         enddo
!      endif
!      avar = ''
!      bVar = ''
!      ahalfvar = ''
!      bhalfVar = ''
!      P0Var = ''
!      PsVar = ''

      ! Scaling and adding offset from name-table file for strange unit conversion etc
      !
      if(allocated(lstOffset))then
        do iTmp = 1, size(lstOffset)
          if(adjustl(nf%nVars(iVar)%chVarNm) == adjustl(lstVarName(iTmp)))then
!            call msg("Increasing scalefactor nametable from, by", nf%nVars(iVar)%scaleFactor , lstFactor(iTmp))
            nf%nVars(iVar)%scaleFactor =  lstFactor(iTmp)
            nf%nVars(iVar)%offset =  lstOffset(iTmp)
          endif
        enddo
     ! else
      ! FIXME Handle unit conversion on our own
     !        nf%nVars(iVar)%quantity
     !   if (nf%nVars(iVar)%unit == '') then 
     !           call set_error("Unknown unit")
     !
     !
     !
      endif
      
    enddo !var


    ! Varid-s for dimvars
    do iDim = 1, nf%n_Dims
      ifFound = .false.
      do iVar = 1, nf%n_Vars
        if(fu_str_u_case(adjustl(nf%nDims(iDim)%varName)) == fu_str_u_case(adjustl(nf%nVars(iVar)%chVarNm)))then
          nf%nDims(iDim)%varId = nf%nVars(iVar)%varId
          ifFound = .true.
          nf%nDims(iDim)%defined =  .true.
          exit
        endif
      enddo
      if(.not. ifFound)then
        if(ifMsgs)then
          call msg_warning(fu_connect_strings('Cannot find varId for dimvar: ', nf%nDims(iDim)%varName, ' for dimension: ', nf%nDims(iDim)%dimName),sub_name)
        endif
        nf%nDims(iDim)%defined = .false.
      endif
    enddo


    ! Read dimension variables. For now all dimensions should be defined and all dim-variables known

    do iTmp = 1, nf%n_Dims
      if(nf%nDims(iTmp)%defined)then
        if(ifMsgs) call msg("Processong dimension "//trim(nf%nDims(iTmp)%varName)//". Axis:"//trim(nf%nDims(iTmp)%axis))
        !TIME
        !
        if(fu_str_l_case(adjustl(nf%nDims(iTmp)%axis)) == 't')then
          nf%n_times = nf%nDims(iTmp)%dimlen
          if(ifMsgs) call msg('Time values:')
          allocate(nf%nTime%times(nf%n_times))
          if(nf%nVars(nf%nDims(iTmp)%varId)%xType == nf90_char)then 
            ! assume standard netcdf string (yyyy-mm-dd hh:mm:ss)
            jTmp = nf%nDims(nf%nVars(nf%nDims(iTmp)%varId)%dimIds(1))%dimLen
            allocate(arrChTmp(jTmp), stat = istat)
            if (istat /= 0)then 
                call msg_warning('Failed to allocate arrChTmp of size: '//fu_str(jTmp), sub_name)
            endif

            do iT = 1, nf%nDims(iTmp)%dimLen
              iStat = NF90_get_var(nf%unit_bin, &
                                 & nf%nVars(nf%nDims(iTmp)%varId)%varId, &
                                 & arrChTmp, &
                                 & start=(/1, iT/), &
                                 & count=(/nf%nDims(nf%nVars(nf%nDims(iTmp)%varId)%dimIds(1))%dimLen, 1/)) 
              write(chTmp, fmt=*)arrChTmp(1:4), '-', arrChTmp(6:7), '-', arrChTmp(9:10), &
                   & '_',arrChTmp(12:13), ':',arrChTmp(15:16), ':', arrChTmp(18:19)  
              nf%nTime%times(iT) =  fu_netcdf_str_to_silam_time(adjustl(chTmp))    
              if(ifMsgs) call report(nf%nTime%times(iT))
              if(error)return          
            enddo
            deallocate(arrChTmp)
          else ! not char format
            
            nf%nVars(nf%nDims(iTmp)%varId)%unit = fu_str_l_case(adjustl(trim(nf%nVars(nf%nDims(iTmp)%varId)%unit)))
            iT = index(nf%nVars(nf%nDims(iTmp)%varId)%unit, 'since')
            if(iT > 0)then
              chTmp = nf%nVars(nf%nDims(iTmp)%varId)%unit(iT+5:)
              nf%nTime%tDimStart = fu_netcdf_str_to_silam_time(adjustl(trim(chTmp)))
            else
              call set_error('Cannot parse time unit string', sub_name)
            endif
            if(ifMsgs)then
                call msg('Time dimension start: ')
                call report(nf%nTime%tDimStart)
            endif
            
            read(unit=nf%nVars(nf%nDims(iTmp)%varId)%unit, iostat=iStat, fmt=*) chTmp
            select case(chTmp)
              case('day', 'days', 'd')
                nf%nVars(nf%nDims(iTmp)%varid)%Unit = 'DAY'
              case('hour', 'hours', 'h')
                nf%nVars(nf%nDims(iTmp)%varid)%Unit = 'HR'
              case('minute', 'minutes', 'm')
                nf%nVars(nf%nDims(iTmp)%varid)%Unit = 'MIN'
              case('second', 'seconds', 's')
                nf%nVars(nf%nDims(iTmp)%varid)%Unit = 'SEC'
              case default
                call set_error(fu_connect_strings('Unsupported time unit: ', &
                                                  & nf%nVars(nf%nDims(iTmp)%varId)%Unit) ,sub_name)
                return
            end select
            if(ifMsgs)then
              call msg('Time unit:' + nf%nVars(nf%nDims(iTmp)%varid)%Unit)
            endif                                      
            
!              nf%nDims(iTmp)%values => fu_work_array()
              allocate(nf%nDims(iTmp)%values(nf%nDims(iTmp)%dimLen))
              iStat = NF90_get_var(ncid=nf%unit_bin, &
                                 & varid=nf%nDims(iTmp)%varId, &
                                 & values=nf%nDims(iTmp)%values, &
                                 & start=(/1/), &
                                 & count=(/nf%nDims(iTmp)%dimLen/)) 

              if(ifMsgs) call msg('Time values are:', nf%nDims(iTmp)%values(1:nf%nDims(iTmp)%dimLen))
              if (istat /= 0)then 
                call msg_warning('Failed to get time values',sub_name)
                call msg_warning(fu_connect_strings('Netcdf error13:', nf90_strerror(iStat)),sub_name)
                nf%nDims(iTmp)%defined = .false.       ! undefine the dimension
              endif 
              if(ifMsgs)then
                call msg(fu_connect_strings('Calendar:',nf%nVars(nf%nDims(iTmp)%varid)%calendar))
              endif
              if(ifMsgs)then
                do iT = 1, nf%nDims(iTmp)%dimLen
                  call msg('Index & value:',iT,nf%nDims(iTmp)%values(iT))
                end do
              endif
              !
              ! Determine the type of calendar
              !
              nf%nVars(nf%nDims(iTmp)%varid)%calendar = &
                   & fu_str_l_case(adjustl(nf%nVars(nf%nDims(iTmp)%varid)%calendar))
              if(nf%nVars(nf%nDims(iTmp)%varid)%calendar == '' .or. &
               & nf%nVars(nf%nDims(iTmp)%varid)%calendar == 'none')then
                call msg_warning('Calendar not known, assuming standard',sub_name)
                nf%nVars(nf%nDims(iTmp)%varId)%calendar = 'standard'
              endif
              if(index(nf%nVars(nf%nDims(iTmp)%varid)%calendar, 'standard') > 0 &
               & .or.  index('standard', nf%nVars(nf%nDims(iTmp)%varid)%calendar) > 0) then
                nf%nVars(nf%nDims(iTmp)%varId)%calendar = 'standard'
              endif
              
              call parse_times(nf%nTime%tDimStart, nf%nDims(iTmp)%values, nf%nDims(iTmp)%dimLen, &
                             & nf%nVars(nf%nDims(iTmp)%varId)%unit, nf%nVars(nf%nDims(iTmp)%varid)%calendar, &
                             & nf%nTime%times)
              
              if(ifMsgs)then
                call msg('Times found in the NetCDF file:')
                do iT = 1, nf%nDims(iTmp)%dimLen             
                  call report(nf%nTime%times(iT))
                end do
              endif
          endif ! char format

          nf%nTime%first_valid_time = nf%nTime%times(1)
          nf%nTime%last_valid_time = nf%nTime%times(nf%nDims(iTmp)%dimLen)
          if(nf%nDims(iTmp)%dimLen > 1) then
            nf%nTime%step = nf%nTime%times(2) - nf%nTime%times(1)
          else
            nf%nTime%step = zero_interval  !Still valid in arithmetics
          endif



          if(nf%nTime%analysis_time == time_missing)then
            if(nf%nTime%analysis_t_from == dim_start)then
              nf%nTime%analysis_time = nf%nTime%tDimStart
            elseif(nf%nTime%analysis_t_from == first_value)then
              nf%nTime%analysis_time = nf%nTime%first_valid_time
              !!! Dirty hack to avoid 00:30 analysis on mid_priod hourly averages
              if (nf%nTime%time_label_position == mid_of_period) then
                nf%nTime%analysis_time = nf%nTime%analysis_time - nf%nTime%step*0.5
              endif
            else
              call set_error('Failed to set analysis time', sub_name)
            endif
          endif
          call msg('Analysis time:')
          call report(nf%nTime%analysis_time)
        
        !
        !LEVEL
        !
        elseif(fu_str_l_case(adjustl(nf%nDims(iTmp)%axis)) == 'z')then

!          nf%nDims(iTmp)%values => fu_work_array()
          allocate(nf%nDims(iTmp)%values(nf%nDims(iTmp)%dimLen))
          iStat = NF90_get_var(ncid=nf%unit_bin, &
                             & varid=nf%nDims(iTmp)%varId, &
                             & values=nf%nDims(iTmp)%values, &
                             & start=(/1/), &
                             & count=(/nf%nDims(iTmp)%dimLen/)) 
          if (istat /= 0)then 
            call msg_warning('Failed to get level values',sub_name)
            call msg_warning(fu_connect_strings('Netcdf error14:', nf90_strerror(iStat)),sub_name)
            nf%nDims(iTmp)%defined = .false.       ! undefine the dimension
          endif
          if (ifMsgs) call msg("Z-dimension type: "//trim(nf%nDims(iTmp)%chType))

          allocate(lstLevs(nf%nDims(iTmp)%dimlen+1))
          chTmp = fu_str_u_case(adjustl(nf%nDims(iTmp)%chType)) 
          if( chTmp == 'PRESSURE')then
            do iLev = 1, nf%nDims(iTmp)%dimLen
              lstLevs(ilev) = fu_set_pressure_level(real(nf%nDims(iTmp)%values(iLev)))
            enddo
            lstLevs(nf%nDims(iTmp)%dimlen+1) = level_missing
            call set_vertical(lstLevs, nf%nDims(iTmp)%sVert)
          elseif( chTmp == 'SIGMA_LEVEL')then
            do iLev = 1, nf%nDims(iTmp)%dimLen
              lstLevs(ilev) = fu_set_sigma_level(real(nf%nDims(iTmp)%values(iLev)))
            enddo
            lstLevs(nf%nDims(iTmp)%dimlen+1) = level_missing
            call set_vertical(lstLevs, nf%nDims(iTmp)%sVert)
          elseif( chTmp == 'HEIGHT_FROM_SURF')then
            do iLev = 1, nf%nDims(iTmp)%dimLen
              lstLevs(ilev) = fu_set_constant_height_level(real(nf%nDims(iTmp)%values(iLev)))
            enddo
            lstLevs(nf%nDims(iTmp)%dimlen+1) = level_missing
            call set_vertical(lstLevs, nf%nDims(iTmp)%sVert)
          elseif( chTmp == 'HEIGHT_LYR_FROM_SURF')then
            bottom = 0
            do iLev = 1, nf%nDims(iTmp)%dimLen
              thickness = 2.0 * (real(nf%nDims(iTmp)%values(iLev)) - bottom)
              lstLevs(ilev) = fu_set_layer_between_two(layer_btw_2_height, bottom+thickness, bottom)
              bottom = bottom + thickness
            enddo
            lstLevs(nf%nDims(iTmp)%dimlen+1) = level_missing
            call set_vertical(lstLevs, nf%nDims(iTmp)%sVert)
           elseif( chTmp == 'HEIGHT_LYR_BY_TOP')then
            bottom = 0
            do iLev = 1, nf%nDims(iTmp)%dimLen
              thickness = real(nf%nDims(iTmp)%values(iLev)) - bottom
              lstLevs(ilev) = fu_set_layer_between_two(layer_btw_2_height, bottom+thickness, bottom)
              bottom = bottom + thickness
            enddo
            lstLevs(nf%nDims(iTmp)%dimlen+1) = level_missing
            call set_vertical(lstLevs, nf%nDims(iTmp)%sVert)
          elseif( chTmp == 'ALTITUDE_FROM_SEA_LEVEL')then
            do iLev = 1, nf%nDims(iTmp)%dimLen
              lstLevs(ilev) = fu_set_constant_altitude_level(real(nf%nDims(iTmp)%values(iLev)))
            enddo
            lstLevs(nf%nDims(iTmp)%dimlen+1) = level_missing
            call set_vertical(lstLevs, nf%nDims(iTmp)%sVert)
          elseif( chTmp == 'HYBRID')then
            if (.not. associated(nf%nDims(iTmp)%a)) then
              ! Get hybrid coefficients
              allocate(nf%nDims(iTmp)%a(nf%nDims(iTmp)%dimlen))
              allocate(nf%nDims(iTmp)%b(nf%nDims(iTmp)%dimlen))
              nf%nDims(iTmp)%a(:) = int_missing
              nf%nDims(iTmp)%b(:) = int_missing         
              do iVar = 1, nf%n_Vars 
                if(fu_str_l_case(adjustl(nf%nVars(iVar)%chVarNm)) == nf%nDims(iTmp)%a_Var)then
                  iStat = NF90_get_var(nf%unit_bin, nf%nVars(iVar)%varId, nf%nDims(iTmp)%a)
                elseif(fu_str_l_case(adjustl(nf%nVars(iVar)%chVarNm)) == nf%nDims(iTmp)%b_Var)then
                  iStat = NF90_get_var(nf%unit_bin, nf%nVars(iVar)%varId, nf%nDims(iTmp)%b)
                elseif(fu_str_l_case(adjustl(nf%nVars(iVar)%chVarNm)) == nf%nDims(iTmp)%P0_Var)then
                  iStat = NF90_get_var(nf%unit_bin, nf%nVars(iVar)%varId, nf%nDims(iTmp)%P0)
                endif
                if(iStat /= 0)then
                  if(ifMsgs)then
                    call msg_warning(fu_connect_strings('Failed to read the dimension variable:', &
                                                    & nf%nDims(iTmp)%varName), &
                                                    & sub_name)
                  endif
                  nf%nDims(iTmp)%defined = .false.       ! undefine the dimension
                  cycle
                endif
              enddo
            end if
            if(nf%nDims(iTmp)%a(1)==int_missing)then 
              if(nf%nDims(iTmp)%b(1)==int_missing)then
                call set_error('Hybrid coefficients missing', sub_name)
                return
              else
                !
                ! Hybrid levels can be defined either by one or two coefficients.
                ! Here handle the case where only one coefficient is given:
                ! pressure = coef/1000*surface_pressure
                ! a = 0; b = coef/1000
                !
                nf%nDims(iTmp)%a(:) = 0. 
                do iStat = 1, size(nf%nDims(iTmp)%a)
                  nf%nDims(iTmp)%b(iStat) = nf%nDims(iTmp)%b(iStat)/1000
                enddo
                nf%nDims(iTmp)%P0(1) = 0.
              endif
            endif
            do iLev = 1, nf%nDims(iTmp)%dimLen
              lstLevs(ilev) = fu_set_hybrid_level(iLev, nf%nDims(iTmp)%a(iLev)*nf%nDims(iTmp)%P0(1), &
                                                & nf%nDims(iTmp)%b(iLev))
            enddo  
            lstLevs(nf%nDims(iTmp)%dimlen+1) = level_missing
            call set_vertical(lstLevs, nf%nDims(iTmp)%sVert)
          elseif( chTmp == 'HYBRID_LYR')then
            if (.not. associated(nf%nDims(iTmp)%a_half)) then
              ! Get hybrid coefficients
              allocate(nf%nDims(iTmp)%a_half(nf%nDims(iTmp)%dimlen+1))
              allocate(nf%nDims(iTmp)%b_half(nf%nDims(iTmp)%dimlen+1))
              nf%nDims(iTmp)%a_half(:) = int_missing
              nf%nDims(iTmp)%b_half(:) = int_missing         
              call msg("Searching a_half"+nf%nDims(iTmp)%a_half_Var)
              call msg("Searching b_half"+nf%nDims(iTmp)%b_half_Var)
              do iVar = 1, nf%n_Vars
                call msg("Trying:"+fu_str_l_case(adjustl(nf%nVars(iVar)%chVarNm)))
                if(fu_str_l_case(adjustl(nf%nVars(iVar)%chVarNm)) == nf%nDims(iTmp)%a_half_Var)then
                  call msg("Found:"+nf%nDims(iTmp)%a_half_Var)
                  iStat = NF90_get_var(nf%unit_bin, nf%nVars(iVar)%varId, nf%nDims(iTmp)%a_half)
                  
                elseif(fu_str_l_case(adjustl(nf%nVars(iVar)%chVarNm)) == nf%nDims(iTmp)%b_half_Var)then
                  call msg("Found:"+nf%nDims(iTmp)%b_half_Var)
                  iStat = NF90_get_var(nf%unit_bin, nf%nVars(iVar)%varId, nf%nDims(iTmp)%b_half)
                elseif(fu_str_l_case(adjustl(nf%nVars(iVar)%chVarNm)) == nf%nDims(iTmp)%P0_Var)then
                  iStat = NF90_get_var(nf%unit_bin, nf%nVars(iVar)%varId, nf%nDims(iTmp)%P0)
                endif
                if(iStat /= 0)then
                  if(ifMsgs)then
                    call msg_warning(fu_connect_strings('Failed to read the dimension variable:', &
                                                    & nf%nDims(iTmp)%varName), &
                                                    & sub_name)
                  endif
                  nf%nDims(iTmp)%defined = .false.       ! undefine the dimension
                  cycle
                endif
              enddo
            end if
            if ((nf%nDims(iTmp)%a_half(1)==int_missing) .or. &
               & (nf%nDims(iTmp)%b_half(1)==int_missing) ) then
                call set_error('Hybrid half coefficients missing', sub_name)
                return
            endif

            ! Check only first layer pressure drop, so layers are always bottom-top
            fTmp = (nf%nDims(iTmp)%a_half(2) - nf%nDims(iTmp)%a_half(1) )  &
                  &  + (nf%nDims(iTmp)%b_half(2) - nf%nDims(iTmp)%b_half(1))*std_pressure_sl
!            call msg("a_half", nf%nDims(iTmp)%a_half(1:nf%nDims(iTmp)%dimLen+1))
!            call msg("b_half", nf%nDims(iTmp)%b_half(1:nf%nDims(iTmp)%dimLen+1))
!            call ooops("ooops")
            
            do iLev = 1, nf%nDims(iTmp)%dimLen
              if (fTmp < 0) then 
                ! pressure decreases with lev number: bottopm-up
                lstLevs(ilev) = fu_set_layer_between_two(layer_btw_2_hybrid, ilev, ilev-1, &
                        & nf%nDims(iTmp)%a_half(iLev+1), nf%nDims(iTmp)%b_half(iLev+1), &  !top
                        & nf%nDims(iTmp)%a_half(iLev), nf%nDims(iTmp)%b_half(iLev))        !bottom
              else
                ! top-down
                lstLevs(ilev) = fu_set_layer_between_two(layer_btw_2_hybrid, ilev-1, ilev, &
                        & nf%nDims(iTmp)%a_half(iLev), nf%nDims(iTmp)%b_half(iLev), &    !top
                        & nf%nDims(iTmp)%a_half(iLev+1), nf%nDims(iTmp)%b_half(iLev+1))  !bottom

              endif
            enddo  
            lstLevs(nf%nDims(iTmp)%dimlen+1) = level_missing
            call set_vertical(lstLevs, nf%nDims(iTmp)%sVert)
          elseif( chTmp == 'SURFACE_LEVEL')then
            call set_vertical(surface_level, nf%nDims(iTmp)%sVert)
          elseif( chTmp == 'TOP_ATMOSPHERE_LEVEL')then
            call set_vertical(top_atmosphere_level, nf%nDims(iTmp)%sVert)
          elseif( chTmp == 'MEAN_SEA_LEVEL')then
            call set_vertical(mean_sea_level, nf%nDims(iTmp)%sVert)
          elseif( chTmp == 'ENTIRE_ATMOSPHERE_MEAN_LAYER')then
            call set_vertical(entire_atmosphere_mean_level, nf%nDims(iTmp)%sVert)
          elseif( chTmp == 'ENTIRE_ATMOSPHERE_INTEGR_LAYER')then
            call set_vertical(entire_atmosphere_integr_level, nf%nDims(iTmp)%sVert)
          else
            call set_error(fu_connect_strings('Strange type of level:',nf%nDims(iTmp)%chType), &
                                            & sub_name)
            return
          endif
          deallocate(lstLevs)

          !FIXME Levels here are as they are in file. No reshuffling!!!!!
!          call arrange_levels_in_vertical(nf%nDims(iTmp)%sVert, ifChanged)
          nf%nDims(iTmp)%ifReverse = .False.
           
          ! assuming, that it was upside down, not a random mess ...

        !LON
        ! LAT and LON coordinates can have either 1d or 2d values. If also has time or z dimension, the last dimensions will be ignored. 
        ! Here just check the units and read the values.
        !
        elseif(fu_str_l_case(adjustl(nf%nDims(iTmp)%axis)) == 'x')then
          if(nf%nVars(nf%nDims(iTmp)%varid)%n_dims == 1)then ! easy case - 1d
            allocate(nf%nDims(iTmp)%values(nf%nDims(iTmp)%dimLen))            
            iStat = NF90_get_var(ncid=nf%unit_bin, &
                               & varid=nf%nDims(iTmp)%varId, &
                               & values=nf%nDims(iTmp)%values, &
                               & start=(/1/), &
                               & count=(/nf%nDims(iTmp)%dimLen/)) 
            if (istat /= 0)then 
              call msg_warning('Failed to get lon values',sub_name)
              call msg_warning(fu_connect_strings('Netcdf error15:', nf90_strerror(iStat)),sub_name)
              nf%nDims(iTmp)%defined = .false.       ! undefine the dimension
            endif
            if(nf%nDims(iTmp)%values(2) < nf%nDims(iTmp)%values(1))then
              nf%nDims(iTmp)%ifReverse = .true.
              call msg('Reversed x dimension')
            else
              nf%nDims(iTmp)%ifReverse = .false.
            endif
          else ! more than 1d. As here it might yet not be known, which of those are important, leave for later when reading gridvars

!            call msg(fu_connect_strings('Nr of dims for: ', nf%nDims(iTmp)%varName), nf%nVars(nf%nDims(iTmp)%varid)%n_dims)

          endif
        
        !LAT
        !
        elseif(fu_str_l_case(adjustl(nf%nDims(iTmp)%axis)) == 'y')then
          if(nf%nVars(nf%nDims(iTmp)%varid)%n_dims == 1)then ! easy case - 1d

            allocate(nf%nDims(iTmp)%values(nf%nDims(iTmp)%dimLen))            
            iStat = NF90_get_var(ncid=nf%unit_bin, &
                               & varid=nf%nDims(iTmp)%varId, &
                               & values=nf%nDims(iTmp)%values, &
                               & start=(/1/), &
                               & count=(/nf%nDims(iTmp)%dimLen/)) 
            if (istat /= 0)then 
              call msg_warning('Failed to get lat values',sub_name)
              call msg_warning(fu_connect_strings('Netcdf error16:', nf90_strerror(iStat)),sub_name)
              nf%nDims(iTmp)%defined = .false.       ! undefine the dimension
            endif
            if(nf%nDims(iTmp)%values(2) < nf%nDims(iTmp)%values(1))then
              nf%nDims(iTmp)%ifReverse = .true.
              call msg('Reversed y dimension')
            else
              nf%nDims(iTmp)%ifReverse = .false.
            endif
          else ! more than 1d. As here it might yet not be known, which of those are important, leave for later when reading gridvars
            
!            call msg(fu_connect_strings('Nr of dims for: ', nf%nDims(iTmp)%varName), nf%nVars(nf%nDims(iTmp)%varid)%n_dims)
          
          endif
        else
          call msg_warning(fu_connect_strings('Failed to attribute the dimension variable:', &
                                            & nf%nDims(iTmp)%varName), &
                         & sub_name)
        endif !t, z, x, y
      else        
        ! Try to identify undefined z dim vithout dim var
        ! Could also be string length or something else funny
        ! Only thing to check is dim name
        if(fu_str_l_case(adjustl(nf%nDims(iTmp)%dimName)) == 'lev' .or. &
         & fu_str_l_case(adjustl(nf%nDims(iTmp)%dimName)) == 'ilev' .or. &
         & fu_str_l_case(adjustl(nf%nDims(iTmp)%dimName)) == 'level' .or. &
         & fu_str_l_case(adjustl(nf%nDims(iTmp)%dimName)) == 'height' .or. &
         & fu_str_l_case(adjustl(nf%nDims(iTmp)%dimName)) == 'hybrid' .or. &
         & fu_str_l_case(adjustl(nf%nDims(iTmp)%dimName)) == 'bottom_top' .or. &
         & fu_str_l_case(adjustl(nf%nDims(iTmp)%dimName)) == 'bottom_top_stag' .or. &
         & fu_str_l_case(adjustl(nf%nDims(iTmp)%dimName)) == 'ilevel' .or. &
         & fu_str_l_case(adjustl(nf%nDims(iTmp)%dimName)) == 'z')then
          nf%nDims(iTmp)%axis = 'z'
          nf%nDims(iTmp)%varId = int_missing
          nf%nDims(iTmp)%varName = ''
          nf%nDims(iTmp)%defined = .true.
        else
          if(ifMsgs)then
            call msg_warning(fu_connect_strings('Unidentified dimension: ',nf%nDims(iTmp)%dimName))
          endif
        endif 

      endif !if defined
      
    enddo  


   ! Check the availibility of vars (availible if silam quantity and dims defined)
   ! If not availible, set quantity to int_missing
   ! Set silam verticals and grid pointers for variables

    do iVar = 1, nf%n_Vars
      ifFound = .false. 

      !Set 2D not initialised earlier. Will be reset for 3D variables
      ! The only chance to get initialised for variables with no Z dimension
      if (nf%nVars(iVar)%n_levs < 0) then
              nf%nVars(iVar)%if3d = .false.
              nf%nVars(iVar)%n_levs = 1
              call set_vertical(surface_level, nf%nVars(iVar)%sVert) 
      endif
       
      do iDim = 1, nf%nVars(iVar)%n_Dims
        if(.not. nf%nDims(nf%nVars(iVar)%dimIds(iDim))%defined)then
          if(ifMsgs)then
            call msg_warning(fu_connect_strings('Variable with undefined dimension:', nf%nDims(nf%nVars(iVar)%dimIds(iDim))%dimName),sub_name)
          endif
          nf%nVars(iVar)%quantity = int_missing
          exit
        endif
!        print *, iVar, iDim, fu_str_l_case(adjustl(trim(nf%nDims(nf%nVars(iVar)%dimIds(iDim))%axis)))
        if(fu_str_l_case(adjustl(trim(nf%nDims(nf%nVars(iVar)%dimIds(iDim))%axis))) == 't')then
          nf%nVars(iVar)%ifTimeDep = .true.
        elseif(fu_str_l_case(adjustl(trim(nf%nDims(nf%nVars(iVar)%dimIds(iDim))%axis))) == 'z')then
          nf%nVars(iVar)%n_levs = nf%nDims(nf%nVars(iVar)%dimIds(iDim))%dimlen
          nf%nVars(iVar)%sVert = nf%nDims(nf%nVars(iVar)%dimIds(iDim))%sVert
          if(nf%nVars(iVar)%n_levs > 1)then
            nf%nVars(iVar)%if3d = .true.        
          endif
       
        elseif(fu_str_l_case(adjustl(trim(nf%nDims(nf%nVars(iVar)%dimIds(iDim))%axis))) == 'x')then  
       
          do iTmp = 1, nf%nVars(iVar)%n_Dims
            if(fu_str_l_case(adjustl(trim(nf%nDims(nf%nVars(iVar)%dimIds(iTmp))%axis))) == 'y')then

              ! First check, if required grid is already available in nGrids array, if yes, set the pointer  
              do iVal = 1, size(nf%nGrids)
                if(nf%nGrids(iVal)%defined == silja_false)exit
                if(nf%nGrids(iVal)%latDimId == nf%nVars(iVar)%dimIds(iTmp) .and. &
                 & nf%nGrids(iVal)%lonDimId == nf%nVars(iVar)%dimIds(iDim))then
                  ifFound = .true.
                  nf%nVars(iVar)%nGrid => nf%nGrids(iVal)
                  exit
                endif
              enddo   
              if(ifFound)exit
              ! Grid not there yet - have to define a new one
              ! If 1d dim vars, assume lonlat grid
              nf%nGrids(iVal)%latDimId = nf%nVars(iVar)%dimIds(iTmp)
              nf%nGrids(iVal)%lonDimId = nf%nVars(iVar)%dimIds(iDim)
              nf%nGrids(iVal)%latVarId = nf%nDims(nf%nVars(iVar)%dimIds(iTmp))%varId
              nf%nGrids(iVal)%lonVarId = nf%nDims(nf%nVars(iVar)%dimIds(iDim))%varId
              nf%nGrids(iVal)%nx = nf%nDims(nf%nVars(iVar)%dimIds(iDim))%dimLen
              nf%nGrids(iVal)%ny = nf%nDims(nf%nVars(iVar)%dimIds(iTmp))%dimLen
              nf%nGrids(iVal)%sPole_lon = sPoleLon_tmp
              nf%nGrids(iVal)%sPole_lat = sPoleLat_tmp
              nf%nGrids(iVal)%defined = silja_true
              nf%nVars(iVar)%nGrid => nf%nGrids(iVal)
              ifFound = .true.
              exit
            endif 
          enddo
        endif
      enddo
      if(.not. ifFound)then ! grid missing
         nf%nVars(iVar)%quantity = int_missing
         if(ifMsgs)then
           call msg_warning('Variable with undefined grid',sub_name)
         endif
      endif
      if(nf%nVars(iVar)%quantity == int_missing)then
        if(ifMsgs)then
          call msg_warning(fu_connect_strings('Unavailible variable: ', nf%nVars(iVar)%chVarNm))
        endif
        nf%nVars(iVar)%defined = silja_false        
      else
        if(ifMsgs)call msg(fu_connect_strings('Availible variable: ', nf%nVars(iVar)%chVarNm))
      endif
    enddo 

    !
    ! Now all actually used grids should be in nGrids array, all variables pointing at correct grid.
    ! Corresponding SILAM grids have to be defined 
    ! For not lonlat grids, now both x and y dimensions are known, so lat and lon fields can be read.
    ! Some other variables might also be available (map factor, rotation); those have to be defined 
    ! in nc name-table as gridvars
    !
    vals => fu_work_array()
    do iTmp = 1, size(nf%nGrids)
 
      if(.not. nf%nGrids(iTmp)%defined == silja_true)cycle
      if(nf%nGrids(iTmp)%lonVarId == int_missing .or. nf%nGrids(iTmp)%latVarId == int_missing)then
        call set_error('grid trouble 1',sub_name)
        return
      endif
      if(nf%nGrids(iTmp)%lonDimId == int_missing .or. nf%nGrids(iTmp)%latDimId == int_missing)then
        call set_error('grid trouble 2',sub_name)
        return
      endif


      nx = nf%nGrids(iTmp)%nx  
      ny = nf%nGrids(iTmp)%ny
      if(nf%nVars(nf%nGrids(iTmp)%latVarId)%n_dims == 1 .and. nf%nVars(nf%nGrids(iTmp)%lonVarId)%n_dims == 1)then 

          !!! Handling of wrapped X dimensions like 358 359 0 1 2
          f8arrPtr => nf%nDims(nf%nGrids(iTmp)%lonDimId)%values(1:nx)
          if (nf%nDims(nf%nGrids(iTmp)%lonDimId)%ifReverse) then
            x1 = f8arrPtr(nx)
            dx = f8arrPtr(1) - f8arrPtr(nx)
          else
            x1 = f8arrPtr(1)
            dx = f8arrPtr(nx) - f8arrPtr(1)
          endif 
          if (dx < 0) dx = dx + 360.  !!Add wrap
          dx = dx / (nx - 1)

          f8arrPtr => nf%nDims(nf%nGrids(iTmp)%latDimId)%values(1:ny)
          if (nf%nDims(nf%nGrids(iTmp)%latDimId)%ifReverse) then
            y1 = f8arrPtr(ny)
            dy = (f8arrPtr(1) - f8arrPtr(ny)) / (ny - 1)
          else
            y1 = f8arrPtr(1)
            dy = (f8arrPtr(ny) - f8arrPtr(1)) / (ny - 1)
          endif 

          nf%nGrids(iTmp)%sGrid &
               & = fu_set_lonlat_grid('gridfromnc',  x1, y1,  .false., nx, ny, &
                                    & fu_set_pole(south_flag, nf%nGrids(iTmp)%sPole_lat,  nf%nGrids(iTmp)%sPole_lon), &
                                    & dx, dy)
      else ! 2d lon and lat fields, other dimensions will be ignored
        !
        ! need the values of lats and lons, to set the anygrid, for checking if it already exists and not 
        ! eating all the memory in the world allocating new anygrid parameters each time new file is read
        !
        vals(1) = real_missing  !!Temp storage for lon-lat
        fAtt(1) = real_missing
        if(allocated(gVarId))then
          do iVar = 1, size(gVarId)
            if (gVarId(iVar) < 1) cycle ! Not all listed quantities are bound to exist in the file
            if(.not. (any(nf%nvars(gVarId(iVar))%dimIds == nf%nGrids(iTmp)%lonDimId) .and. &
             & any(nf%nvars(gVarId(iVar))%dimIds == nf%nGrids(iTmp)%latDimId)))cycle
            if(gVarNm(iVar)=="lon")then
              call read_2d_gridvar(nf, gVarId(iVar), gVarNm(iVar), refVarId(iVar), refVarNm(iVar), vals)
            elseif(gVarNm(iVar)=="lat")then
              call read_2d_gridvar(nf, gVarId(iVar), gVarNm(iVar), refVarId(iVar), refVarNm(iVar), fAtt)
            endif
          enddo

          if((vals(1) ==  real_missing) .or. (fAtt(1) == real_missing)) then 
            call set_error('lats or lons unfound',sub_name)
            return
          endif
          lats2Dptr(1:nx,1:ny) => fAtt(1:nx*ny)
          lons2Dptr(1:nx,1:ny) => vals(1:nx*ny) 
          if2DLatLon = all(lons2Dptr(1:nx,1) == lons2Dptr(1:nx,ny))
          if (if2DLatLon) then
            !Few additional checks
            if2DLatLon = all(lats2Dptr(1,1:ny) == lats2Dptr(1,1:ny))
            if (if2DLatLon) then
              if2DLatLon = all(lats2Dptr(1,1) == lats2Dptr(1:nx,1))
              if (if2DLatLon) then
                if2DLatLon = all(lons2Dptr(1,1) == lons2Dptr(1,1:ny))
                if (if2DLatLon) then
                  !!Check for even spacing of latitudes beginning, within 1% of mean spacing (Mercator grids should fail here)
                  dy = (lats2Dptr(1,ny) - lats2Dptr(1,1))/(ny-1)
                  if2DLatLon = all( (/abs(lats2Dptr(1,2) - lats2Dptr(1,1) - dy),  &   
                                   & abs(lats2Dptr(1,ny) - lats2Dptr(1,ny-1) - dy) /) < 0.01*abs(dy)) 

                endif
              endif
            endif
          endif

          if (if2DLatLon) then
            dx =  abs(lons2Dptr(1,1)-lons2Dptr(nx,1))/ (nx - 1)
            dy =  abs(lats2Dptr(1,1)-lats2Dptr(1,ny))/ (ny - 1)
            nf%nGrids(iTmp)%sGrid &
               & = fu_set_lonlat_grid('gridfromnc2d', &
                                    & real(minval(lons2Dptr(1:nx,1)),DEFAULT_REAL_KIND), &
                                    & real(MINVAL(lats2Dptr(1,1:ny)),DEFAULT_REAL_KIND), &
                                    & .false., &
                                    & nf%nGrids(iTmp)%nx,  nf%nGrids(iTmp)%ny, &
                                    & pole_geographical, & !! WRF has POLE_LAT attrinute, but for NORHERN pole
                                    & dx, dy)
!!                                  nf%nGrids(iTmp)%sPole_lat and nf%nGrids(iTmp)%sPole_lon are wrong here!




          else
          ! True ANYGRID
            nf%nGrids(iTmp)%sGrid = fu_set_any_grid('anygridformnc', nx, ny)

            if(error)return
            call setAgLonlatFlds(nf%nGrids(iTmp)%sGrid, vals, fAtt, ifnew)

            !
            ! And now the rest of available grid parameters
            !
            if(.not. ifNew)cycle  ! if grid already exists, umnlikely that we have to overwrite the parameters
            do iVar = 1, size(gVarId)
              if(gVarId(iVar) == int_missing)then
                call set_error('gridvar unfound','')
                return
              endif
              if(any(nf%nvars(gVarId(iVar))%dimIds == nf%nGrids(iTmp)%lonDimId) .and. &
               & any(nf%nvars(gVarId(iVar))%dimIds == nf%nGrids(iTmp)%latDimId))then

                call read_2d_gridvar(nf, gVarId(iVar), gVarNm(iVar), refVarId(iVar), refVarNm(iVar), vals)
                if(gVarNm(iVar) == "dx") vals(1:nx*ny) = vals(1:nx*ny)*wrfdxm
                if(gVarNm(iVar) == "dy") vals(1:nx*ny) = vals(1:nx*ny)*wrfdym
                call setAnygridParam(nf%nGrids(iTmp)%sGrid, gVarNm(iVar), vals)
                if(error)return 
   
              endif
            enddo
            call completeAnygrid3DParam(nf%nGrids(iTmp)%sGrid)
          endif !!!Need for true anygrid
        endif !!1d gridvar
        call report(nf%nGrids(iTmp)%sGrid) 
      endif
    enddo !!Grids

    call free_work_array(iAtt) 
    call free_work_array(fAtt)
    call free_work_array(vals)
    if (defined(nl)) call destroy_namelist(nl) !

    if(allocated(lstDimName))deallocate(lstDimName)
    if(allocated(lstAxis))   deallocate(lstAxis)
    if(allocated(lstDimVar)) deallocate(lstDimVar)
    if(allocated(lstDimType))deallocate(lstDimType)

    if(allocated(lstVarName)) deallocate(lstVarName) 
    if(allocated(lstSilamQ)) deallocate(lstSilamQ)
    if(allocated(lstVType))  deallocate(lstVType)
    if(allocated(lstlevValue))deallocate(lstlevValue)
    if(allocated(lstSubst)) deallocate(lstSubst)
    if(allocated(lstMode)) deallocate(lstMode)
    if(allocated(lstWavelen))deallocate(lstWavelen)
    if(allocated(lstFactor))deallocate(lstFactor)
    if(allocated(lstOffset))deallocate(lstOffset)

    if(allocated(gVarNmNf))deallocate(gVarNmNf)
    if(allocated(gVarNm))deallocate(gVarNm)
    if(allocated(gVarId))deallocate(gVarId)
    if(allocated(refVarNm))deallocate(refVarNm)
    if(allocated(refVarId))deallocate(refVarId)

    if(ifMsgs)then
       call msg('nc file open_i: '+trim(nf%fname), iFile)
    endif
    nfile(iFile)%defined = silja_true


  contains
! ===============================================================================
    subroutine read_2d_gridvar(nf, varId, varNm, refVarId, refVarNm, vals)
      implicit none
       
      TYPE(netcdf_file),POINTER :: nf
      integer :: varId, refVarId
      character (len=nf90_max_name) :: varNm, refVarNm
      integer, dimension(:), pointer :: arrStart, arrCount
      real, dimension(:), pointer :: vals, refVals
      integer :: nxStag, nyStag, iDim, iStat, nxRef, nyRef
        
      allocate(arrStart(nf%nVars(varId)%n_dims), stat = iStat)
      allocate(arrCount(nf%nVars(varId)%n_dims), stat = iStat)

      do iDim = 1, nf%nVars(varId)%n_dims
        arrStart(iDim) = 1
        if(nf%nDims(nf%nvars(varId)%dimIds(iDim))%axis == 'x')then
          arrCount(iDim) = nf%nDims(nf%nvars(varId)%dimIds(iDim))%dimlen
          nxStag = nf%nDims(nf%nvars(varId)%dimIds(iDim))%dimlen
        elseif(nf%nDims(nf%nvars(varId)%dimIds(iDim))%axis == 'y') then 
          arrCount(iDim) = nf%nDims(nf%nvars(varId)%dimIds(iDim))%dimlen
          nyStag = nf%nDims(nf%nvars(varId)%dimIds(iDim))%dimlen
        else
          arrCount(iDim) = 1
        endif
      enddo
      iStat = NF90_get_var(nf%unit_bin, varId, vals, &
                         & start=arrStart, count=arrCount )        
      deallocate(arrStart, arrCount)
      if(iStat /= 0)then
        call msg_warning(fu_connect_strings('Failed to read gridvar:', &
                                          & varNm), &
                                          & sub_name)
        call msg_warning(nf90_strerror(iStat))
        return
      endif
      
      !
      ! fixing 0 lats and lons in WRF staggered grids
      !
      if(fu_str_u_case(adjustl(trim(refVarNm))) /= 'XXX')then
        if(refVarId == int_missing)then
          call set_error('refvar unfound',sub_name)
          return
        endif
        allocate(arrStart(nf%nVars(refVarId)%n_dims), stat = iStat)
        allocate(arrCount(nf%nVars(refVarId)%n_dims), stat = iStat)
        do iDim = 1, nf%nVars(refVarId)%n_dims
          arrStart(iDim) = 1
          if(nf%nDims(nf%nvars(refVarId)%dimIds(iDim))%axis == 'x')then
            arrCount(iDim) = nf%nDims(nf%nvars(refVarId)%dimIds(iDim))%dimlen
            nxRef = nf%nDims(nf%nvars(refVarId)%dimIds(iDim))%dimlen
          elseif(nf%nDims(nf%nvars(refVarId)%dimIds(iDim))%axis == 'y') then 
            arrCount(iDim) = nf%nDims(nf%nvars(refVarId)%dimIds(iDim))%dimlen
            nyRef = nf%nDims(nf%nvars(refVarId)%dimIds(iDim))%dimlen
          else
            arrCount(iDim) = 1
          endif
        enddo
        refVals => fu_work_array() 
        iStat = NF90_get_var(nf%unit_bin, refVarId, refVals, &
                           & start=arrStart, count=arrCount )        
        deallocate(arrStart, arrCount)
        if(iStat /= 0)then
          call msg_warning(fu_connect_strings('Failed to read refvar:', &
                                            & refVarNm), &
                                            & sub_name)
            call msg_warning(nf90_strerror(iStat))
            return
          endif
          call make_staggered_fld(refVals, nxRef, nyRef, vals, nxStag, nyStag)
          call free_work_array(refVals)
        endif
        
      end subroutine read_2d_gridvar
   
      subroutine parse_times(time_start, increments, num_times, incr_unit, calendar_type, times_parsed)
        implicit none
        type(silja_time), intent(in) :: time_start
        real(r8k), dimension(:), intent(in) :: increments
        integer, intent(in) :: num_times
        character(len=*), intent(in) :: incr_unit, calendar_type
        type(silja_time), dimension(:), intent(out) :: times_parsed
        
        integer :: days_year, ref_year, incr_years
        real(r8k) :: reftime_sec, incr_unit_sec, jul_date_start, incr_sec, incr_days, incr_unit_days
        real :: fSec
        character(len=nf90_max_name) :: chInterval
        integer :: nMod, nday, nHr, nMin, iT
        
        select case(calendar_type)
        case('standard','gregorian','proleptic_gregorian','proleptic')
          !
          ! Normal calendar
          !
          write(chInterval, fmt=*) '1 ', trim(fu_str_l_case(incr_unit))
          incr_unit_sec = fu_sec8(fu_set_named_interval(chInterval))
          reftime_sec = silja_time_to_real8(time_start)


          do iT = 1, nf%nDims(iTmp)%dimLen
            !force 1-sec precision 
            times_parsed(iT) = real8_to_silja_time(1D0*nint(reftime_sec + incr_unit_sec * increments(iT),kind=8))
          end do
          
        case('noleap', '365_day', 'all_leap', '366_day', '365_days')
          !
          ! All years are considered to have 365 or 366 days. Below julian dates are taken
          ! wrt this. The algorithm is subject to roundoff errors, but using seconds
          ! avoids this in most practical cases.
          !
          if (calendar_type == 'noleap' .or. calendar_type == '365_day' .or. calendar_type == '365_days') then 
            days_year = 365
            ref_year = 2007
          else
            days_year = 366
            ref_year = 2008
          end if
          
          tTmp = fu_set_time_utc(ref_year, &
                           & fu_mon(time_start), fu_day(time_start), fu_hour(time_start), &
                           & fu_min(time_start), fu_sec(time_start))
          jul_date_start = fu_julian_date_real(tTmp)
          incr_unit_sec = fu_conversion_factor(fu_str_l_case(incr_unit), 'sec')
          
          do iT = 1, nf%nDims(iTmp)%dimLen             
            incr_sec = increments(iT) * real(incr_unit_sec, r8k)
            incr_years = int((jul_date_start + incr_sec / (3600*24)) / days_year)
            ! correct day and month, wrong year:
            times_parsed(iT) = fu_julian_date_to_time(jul_date_start+(incr_sec/(24*3600)) &
                                                    & - incr_years*days_year, ref_year)
            times_parsed(iT) = fu_set_time_utc(fu_year(time_start) + incr_years, &
                                         & fu_mon(times_parsed(iT)), &
                                         & fu_day(times_parsed(iT)), &
                                         & fu_hour(times_parsed(iT)), &
                                         & fu_min(times_parsed(iT)), &
                                         & fu_sec(times_parsed(iT)))
          end do
          if(error)return
          
        case ('360_day')
          !
          ! Stupid year with 30 days per month. CAREFUL
          !
          
          incr_unit_days = fu_conversion_factor(fu_str_l_case(incr_unit), 'day')
          days_year = 360
          
          incr_years = int(increments(iT) / days_year * incr_unit_days) ! years passed
          incr_days = increments(iT) * incr_unit_days - incr_years * days_year  ! fraction of the year, [day]
          nMon = int(incr_days / 30)
          nDay = int(incr_days - nMon * 30)
          nHr = int((incr_days - nMon * 30 - nDay) * 24)
          nMin = int((incr_days - nMon * 30 - nDay - nHr/24) * 1440)
          fSec = (incr_days - nMon * 30 - nDay - nHr/24 - nMin/1440) * 86400
          nMon = fu_mon(time_start) + nMon
          nDay = fu_day(time_start) + nDay
          nHr = fu_hour(time_start) + nHr
          nMin = fu_min(time_start) + nMin
          fSec = fu_mon(time_start) + fSec
          if(fSec >= 60)then 
            fSec = fSec - 60
            nMin = nMin + 1
          endif
          if(nMin >= 60)then
            nMin = nMin - 60
            nHr = nHr + 1
          endif
          if(nHr >= 24)then
            nHr = nHr - 24
            nDay = nDay + 1
          endif
          if(nDay > 30)then
            nDay = nDay - 30
            nMon = nMon + 1
          endif
          if(nMon > 12)then
            nMon = nMon - 1
            incr_years = incr_years + 1
          endif
          
          times_parsed(iT) = fu_set_time_utc(fu_year(time_start) + incr_years, &
                                       & nMon, &
                                       & nDay, &
                                       & nHr, &
                                       & nMin, &
                                       & fSec)
          
        case default
          call set_error(fu_connect_strings('Unsupported calendar type: ', calendar_type), sub_name)
          return
        end select  ! type of calendar
        
      end subroutine parse_times

  end function open_netcdf_file_i

  !******************************************** 

  subroutine read_field_from_netcdf_file(iFile, field_id, grid_data, fill_value)
  ! Reads the requested field(quantity, level, time) from netcdf file iFile to grid_data array
  ! Bad in general, because actual ID found from file is ignored

   implicit none
    ! Imported parameters
    integer, intent(in) :: iFile
    TYPE(silja_field_id), INTENT(in) :: field_id ! Field id requestd
    real, dimension(:), intent(out) :: grid_data  
    real, intent(in) :: fill_value
    ! Local parameters
    integer :: iVar, iLev, iT
    TYPE(silja_field_id) :: id_found  !!!Ignored here


    call field_indices_from_netcdf_file(ifile, field_id, id_found, ivar, it, ilev)
    if (error) return
    call read_field_from_netcdf_indices(ifile, ivar, it, ilev, grid_data, fill_value)

  end subroutine read_field_from_netcdf_file

  !******************************************** 

  subroutine read_field_from_netcdf_indices(iFile, ivar, it, ilev, grid_data, fill_value)
  ! Reads the requested field(quantity, level, time) from netcdf file iFile to grid_data array

   implicit none
    ! Imported parameters
    integer, intent(in) :: iFile, ivar, it, ilev
    real, dimension(:), intent(out) :: grid_data  ! must exist
    real, intent(in) :: fill_value
    ! Local parameters
    integer ::  nx, ny, varId, iTmp, jTmp, idimX, idimY, iQty
    type(netcdf_file), pointer :: nf 
    logical :: ifFound
    integer, dimension(nf90_max_var_dims) :: arrStart, arrCount
    real, dimension(:), pointer :: data_tmp
    character(len=fnlen) :: fnm
    real :: sc1, off1, sc2, off2, missval, fTmp

    
    nf => nfile(iFile)
    varId = nf%nvars(iVar)%varid
    iQty =  nf%nvars(iVar)%quantity


   ! And finally read the requested field
   ! first need to define start and count. For that we need to know the order of the dimensions in the file.
   call grid_dimensions(nf%nvars(iVar)%nGrid%sGrid, nx, ny)
   
   do jTmp = 1, nf%nvars(iVar)%n_dims
      if(nf%nDims(nf%nvars(iVar)%dimIds(jTmp))%axis == 'x')then 
        arrStart(jTmp) = 1
        arrCount(jTmp) = nx
        iDimX = jTmp
      elseif(nf%nDims(nf%nvars(iVar)%dimIds(jTmp))%axis == 'y')then 
        arrStart(jTmp) = 1
        arrCount(jTmp) = ny
        iDimY = jTmp
      elseif(nf%nDims(nf%nvars(iVar)%dimIds(jTmp))%axis == 'z')then
        if(nf%nDims(nf%nvars(iVar)%dimIds(jTmp))%ifReverse)then
          arrStart(jTmp) = nf%nvars(iVar)%n_levs - iLev + 1
        else
          arrStart(jTmp) = iLev
        endif


        arrCount(jTmp) = 1
      elseif(nf%nDims(nf%nvars(iVar)%dimIds(jTmp))%axis == 't')then
        arrStart(jTmp) = iT
        arrCount(jTmp) = 1
      else !unknown
        call set_error(fu_connect_strings('Strange dimension', nf%nDims(nf%nvars(iVar)%dimIds(jTmp))%axis),'read_field_from_netcdf_file')
        return
      endif
   enddo
   
#ifdef DEBUG_NC
    call msg('Reading variable:' + fu_quantity_string(iQty))
!    if ( trim(nf%nvars(iVar)%chVarNm) == "cnc_SO2_gas") call ooops("cnc_SO2_gas")
    call msg("Variable: '"+nf%nvars(iVar)%chVarNm+"'" )
    call msg("Arrstart: ", arrStart(1:nf%nvars(iVar)%n_dims))
    call msg("Arrcount: ", arrCount(1:nf%nvars(iVar)%n_dims))
    call msg("")

#endif

   ! Read the variable:
   if(iDimX > idimY .or. nf%nDims(nf%nvars(iVar)%dimIds(iDimX))%ifReverse .or. nf%nDims(nf%nvars(iVar)%dimIds(iDimY))%ifReverse)then
     data_tmp => fu_work_array(product(arrCount(1:nf%nvars(iVar)%n_dims)))  
     jTmp = nf90_get_var(nf%unit_bin, varId, data_tmp, start = arrStart, count=arrCount)
     if (jTmp /= 0) then
       call msg_warning('Failed to read field from netcdf file', 'read_field_from_netcdf_file') 
       call msg_warning(fu_connect_strings('Netcdf error17: ',nf90_strerror(jTmp)), 'read_field_from_netcdf_file')
       call set_error('Failed to read field from netcdf file', 'read_field_from_netcdf_file')
       return
     endif
     
     ! column--, row--, column-consequtive
     if(iDimX > idimY .and. nf%nDims(nf%nvars(iVar)%dimIds(iDimX))%ifReverse .and. nf%nDims(nf%nvars(iVar)%dimIds(iDimY))%ifReverse)then
!call msg('column--, row--, column-consequtive')
          do jTmp=1,ny
            do iTmp=1,nx
              grid_data(iTmp+(jTmp-1)*nx) = data_tmp(ny-jTmp+1 + (nx-iTmp)*ny)
            end do
          end do

     ! column++, row--, column-consequtive
     elseif(iDimX > idimY .and. nf%nDims(nf%nvars(iVar)%dimIds(iDimY))%ifReverse)then
!call msg('column++, row--, column-consequtive')
          do jTmp=1,ny
            do iTmp=1,nx
              grid_data(iTmp+(jTmp-1)*nx) = data_tmp(ny-jTmp+1 + (iTmp-1)*ny)
            end do
          end do

     ! column--, row++, column-consequtive
     elseif(iDimX > idimY .and. nf%nDims(nf%nvars(iVar)%dimIds(iDimX))%ifReverse)then
 !call msg('column--, row++, column-consequtive')
         do jTmp=1,ny
            do iTmp=1,nx
              grid_data(iTmp+(jTmp-1)*nx) = data_tmp(ny-jTmp+1 + (iTmp-1)*ny)
            end do
          end do

     ! column--, row--, row-consequtive
     elseif(nf%nDims(nf%nvars(iVar)%dimIds(iDimX))%ifReverse .and. nf%nDims(nf%nvars(iVar)%dimIds(iDimY))%ifReverse)then
 !call msg('column--, row--, row-consequtive')
         do jTmp=1,ny
            do iTmp=1,nx
              grid_data(iTmp+(jTmp-1)*nx) = data_tmp(nx-iTmp+1 + (ny-jTmp)*nx)
            end do
          end do

     ! column++, row++, column-consequtive
     elseif(iDimX > idimY)then
! call msg('column++, row++, column-consequtive')
         do jTmp=1,ny
            do iTmp=1,nx
              grid_data(iTmp+(jTmp-1)*nx) = data_tmp(jTmp + (iTmp-1)*ny)
            end do
          end do

     ! column++, row--, row-consesqutive
     elseif(nf%nDims(nf%nvars(iVar)%dimIds(iDimY))%ifReverse)then 
! call msg('column++, row--, row-consesqutive')
       do jTmp=1,ny
          do iTmp=1,nx
            grid_data(iTmp+(jTmp-1)*nx) = data_tmp(iTmp + (ny-jTmp)*nx)
          end do
        end do
! call msg('Val 1, 2', grid_data(1), grid_data(2))
     ! column--, row++, row-consequtive
     elseif(nf%nDims(nf%nvars(iVar)%dimIds(iDimX))%ifReverse)then
!  call msg('column--, row++, row-consequtive')
        do jTmp=1,ny
            do iTmp=1,nx
              grid_data(iTmp+(jTmp-1)*nx) = data_tmp(nx-iTmp+1 + (jTmp-1)*nx)
            end do
          end do
     endif

     call free_work_array(data_tmp)
   
   else ! Normal scanning mode:  010 column++, row++, row-consequtive (copied from grib_api_io)
! call msg('Normal scanning mode')
    jTmp = nf90_get_var(nf%unit_bin, varId, grid_data, start = arrStart, count=arrCount)
     if (jTmp /= 0) then
       call msg_warning('Failed to read field from netcdf file', 'read_field_from_netcdf_file') 
       call msg_warning(fu_connect_strings('Netcdf error18: ',nf90_strerror(jTmp)), 'read_field_from_netcdf_file')
       call set_error('Failed to read field from netcdf file', 'read_field_from_netcdf_file')
       return
     endif
   endif
    

 
   ! Offset and scaleFactor might exist in the file or come from unit conversion. 
   ! Reducing two transforms into one results in numerical errros
   ! that lead to catastrophic results, e.g.negative values of vmr  from Mozart boundaries  (Roux)

   sc1=nf%nVars(iVar)%scaleFactor_nc
   off1=nf%nvars(iVar)%offset_nc
   sc2=nf%nVars(iVar)%scaleFactor
   off2=nf%nvars(iVar)%offset
   

#ifdef DEBUG_NC                      
   call msg('The following factor-offset sets are applied:',(/sc1,off1,sc2,off2/))
#endif

   iTmp = 1 ! No of offset-scales to apply
   if (sc1 /= 1. .or.  off1 /= 0.) then
       if (sc2 /= 1. .or. off2 /= 0.) iTmp = 2
   else
       if (sc2 /= 1. .or.  off2 /= 0.) then
           sc1  = sc2
           off1 = off2
           sc2 = 1.
           off2 = 0.
        else
           iTmp = 0
       endif
   endif

  
   if (nf%nVars(iVar)%missing_value /= real_missing) then ! Full-blown masking
     iTmp = 0
     missval = nf%nVars(iVar)%missing_value
     do jTmp = 1, nx*ny
       if(grid_data(jTmp) == missval)then
         grid_data(jTmp) = fill_value
         iTmp = iTmp + 1
       else
         grid_data(jTmp) = (grid_data(jTmp) * sc1 + off1) * sc2 + off2
       endif
     enddo 

     if (iTmp > 0)then
        call msg(fu_quantity_string(iQty) + ':' + fu_str(iTmp) + &
               & '- missing values found and replaced with', fill_value)
     endif
   else ! No masking
     if (iTmp == 1) grid_data(1:nx*ny) = grid_data(1:nx*ny) * sc1 + off1
     if (iTmp == 2) grid_data(1:nx*ny) = (grid_data(1:nx*ny) * sc1 + off1) * sc2 + off2
#ifdef DEBUG_NC                      
      call msg("No of scalings, nx, ny", (/iTmp, nx,ny/))
      if (iTmp>0) call msg("Field read with  sc1  off1  sc2  off2", (/ sc1, off1, sc2, off2/))
#endif
   endif
#ifdef DEBUG_NC
   call msg("min(field), max(field)", minval(grid_data(1:nx*ny)), maxval(grid_data(1:nx*ny)))
#endif

!!! Dirty hack for concentrations packed as integers
  if (any(iQty == (/concentration_flag, volume_mixing_ratio_flag, emission_intensity_flag, emission_flux_flag/) )) then
    if (any(nf%nVars(iVar)%xtype == (/nf90_int, NF90_SHORT/) )) then
       ! Need to take care of possible errors due to numerics of integer scaling
       ! If the result deviates from zero by less than the discrete of the
       ! variable, set it to zero.
      fTmp = 2*abs(sc1*sc2) !
      if (any(grid_data(1:nx*ny) < -fTmp)) then
        call msg(fu_quantity_string(iQty))
        call msg("Field read with  sc1  off1  sc2  off2", (/ sc1, off1, sc2, off2/))
        call msg("min(field), max(field)", minval(grid_data(1:nx*ny)), maxval(grid_data(1:nx*ny)))
        call set_error("Irrecoverable negatives in the data", "read_field_from_netcdf_file")
        return
      endif
      where(grid_data(1:nx*ny) < fTmp) grid_data(1:nx*ny) = 0.
#ifdef DEBUG_NC
   call msg("min(field), max(field) after hack1", minval(grid_data(1:nx*ny)), maxval(grid_data(1:nx*ny)))
#endif

       
    elseif ( nint(nf%nVars(iVar)%missing_value) == -32767) then
       ! Too bad... The source was unpacked by someone who was unable to adjust minimal value...
       ! Just remove negatives and hope for best
       where(grid_data(1:nx*ny) < 0) grid_data(1:nx*ny) = 0
#ifdef DEBUG_NC
   call msg("min(field), max(field) after hack 2", minval(grid_data(1:nx*ny)), maxval(grid_data(1:nx*ny)))
#endif
    endif
  endif

 
  end subroutine read_field_from_netcdf_indices


  !********************************************

  subroutine close_netcdf_file(inf)
    !
    ! Close the netcdf file
    !
    implicit none

    integer, intent(in)::inf
    integer :: i, iStat, iTmp
    character (len=*), parameter :: sub_name='close_netcdf_file'

    TYPE(netcdf_file),POINTER :: nf
    
    nf => nfile(inf) 

    call msg("Closing nc file", inf)
    if (nf%ncver == 3  .and. nf%ifmpiio) then
#ifdef WITH_PNETCDF

!!!       if (fu_fails(modulo(nf%nf90mpifieldCnt, nf%nf90mpiNfieldsPerTimestep + 1) == 0, "Incomplete timestep was written", sub_name)) return

       if (nf%nf90mpifieldCnt /= 0) then !Have to finish the file
        call flush_nf90mpi_buf(nf, sub_name)
        if (fu_fails(.not. error,"Error set", sub_name)) return
       endif

       i = nf90mpi_buffer_detach(nf%unit_bin)
       if (i /= 0) then
         call set_error("nf90mpi_buffer_detach: "// nf90mpi_strerror(i), sub_name)
       endif
       deallocate (nf%nf90mpiReq, nf%nf90mpiStat, stat=i)
       if (fu_fails(i == 0, "Failed to deallocate nf90MpiStuff", sub_name)) return 
       i =  nf90mpi_close(nf%unit_bin)
       if (fu_fails(i == 0, "Failed to close netcdf file", sub_name)) return 
#else
       call set_error("Compiled without PNETCDF", sub_name)
#endif
    else  
      i =  nf90_close(nf%unit_bin)
      if (i /= 0)then 
        call msg("NF90 Error: "//trim(nf90_strerror(i)),i)
        call set_error('Failed to close netcdf file','close_netcdf_file')
        return
      endif
    endif
    if (  len(trim(nf%tmpfname)) > 0 .and. &
      &  (smpi_global_rank == 0 .or. .not. nf%ifMPIIO ) ) then
      !FIXME Should not rename on error
      if ( RENAME(nf%tmpfname, nf%fname) /=0) then
        call msg("TMP    file: "//trim(nf%tmpfname))
        call msg("target file: "//trim(nf%fname))
        call set_error("Rename from tmp failed",'close_netcdf_file')
        return
      endif
    endif
 
    ! Zero the structure
    nf%nlevs%levels = real_missing
    if(associated(nf%ntime%times))deallocate(nf%ntime%times)
    if (associated(nf%nDims)) then
       do i = 1,nf%n_Dims
         call set_missing(nf%ndims(i))
       enddo
       deallocate(nf%ndims)
      nf%n_Dims=0
    endif

      do i = 1, size(nf%nGrids)
        nf%nGrids(i)%defined = silja_false
      enddo
    if (associated(nf%nvars)) deallocate(nf%nvars)

    nf%fname=''
    nf%tmpfname=''
    nf%chTemplate=''
    nf%title = ''
    nf%met_src = met_src_missing
    nf%unit_bin=-1
     call set_missing(nfile(inf)%silamVertical, .false.)
    nf%n_times=0
    nf%n_vars=0
    nf%n_levs=0
    nf%defined = silja_false


  end subroutine close_netcdf_file


   !************************************************************************

  subroutine write_ctl_file_4_netcdf(inf)

  !
  ! Actually only necessary for GrADS to be able to connect files with time templates
  ! in filename if file opened with xdfopen. Only *.nc file name and time definition is
  ! written in the ctl, the rest will be read from the netcdf file itself. 
  !

    IMPLICIT NONE
    !
    ! Parameters imported with intent IN
    integer, intent(in)::inf

    ! Local stuff
    TYPE(netcdf_file),POINTER :: nf
    INTEGER :: iCtlUnit, iVar, iOff, iStat, iTmp, jTmp, nVar
    integer :: nx, ny
    character (len=nf90_max_name) :: chTmp, chTmp1, chTmp2

    if (smpi_global_rank /= 0) return !Only master writes the file
    nf => nfile(inf)

    iCtlUnit = fu_next_free_unit()
    if(error)then
      call set_error('Failed to get unit for ctl file','write_ctl_file_4_netcdf')
      return
    endif

    open(iCtlUnit,file=fu_connect_strings(nf%fname,'.ctl'))

    if(.not. nf%defined == silja_true)then
      call msg_warning('Fully undefined ctl file','write_ctl_file')
      write(iCtlUnit,*)'FULLY UNDEFINED FILE'
      return
    endif

    !WRITE(iCtlUnit,'(A5,A)') 'DSET ',nf%chTemplate(1:len_trim(nf%chTemplate))
    chTmp=nf%chTemplate(1:len_trim(nf%chTemplate))
    WRITE(iCtlUnit,'(A5,A)') 'DSET ',trim(fu_trim_grads_hat(chTmp,nf%fname))
    WRITE(iCtlUnit,'(A)') 'DTYPE NETCDF'
    WRITE(iCtlUnit,'(A,F0.0)') 'UNDEF ', nf%missing_value

    if(index(nf%chTemplate, '%') > 0)then
      if(nCtl%defined)then 
        write(iCtlUnit,'(A16)') 'OPTIONS TEMPLATE'
        WRITE(iCtlUnit,'(A10,I7,A8,A20,A10)') 'TDEF time ',nCtl%nTimes,' LINEAR ', &
                         & fu_time_to_grads_string(nCtl%start_time), &
                         & fu_interval_to_grads_string(nCtl%step)
      else
        call set_error('nCtl not defined','write_ctl_file_4_netcdf')
        return
      endif
    endif
    
    if(fu_gridtype(nf%nGrids(1)%sGrid)==anygrid)then   ! x and y have to be defined
      call grid_dimensions(nf%nGrids(1)%sGrid, nx, ny)
      WRITE(iCtlUnit,'(A10,I7,A20)') 'XDEF lon  ',nx,'  LINEAR   1.0   1.0'
      WRITE(iCtlUnit,'(A10,I7,A20)') 'YDEF lat  ',ny,'  LINEAR   1.0   1.0'
    else
      WRITE(iCtlUnit,'(A4,1X,A,1X,I7,1X,A,1X,F9.4,1X,F9.4)') 'XDEF',trim(nf%nGrids(1)%chXdimName), &
                    & nf%nGrids(1)%nx,'LINEAR', nf%nGrids(1)%x_start, nf%nGrids(1)%x_step
      WRITE(iCtlUnit,'(A4,1X,A,1X,I7,1X,A,1X,F9.4,1X,F9.4)') 'YDEF',trim(nf%nGrids(1)%chYdimName), &
                    & nf%nGrids(1)%ny,'LINEAR', nf%nGrids(1)%y_start, nf%nGrids(1)%y_step
    endif

    WRITE(iCtlUnit,'(A,1X,A,1X,I5,1X,A,1X,100(F8.1,1x))') 'ZDEF',trim(nf%nlevs%chZdimName), nf%n_levs,'LEVELS', nf%nlevs%std_z(1:nf%n_levs)
    iOff = ftell(iCtlUnit) !Save position to write final number of trajs
    WRITE(iCtlUnit,'(A6,A3)') "VARS  ","XXX"
    iVar = 0
    if (nf%nlevs%dzVarId > 0) then
        WRITE(iCtlUnit,"(A,X,I3,X,A)") 'dz=>dz ', nf%n_levs, ' z  Layer thickness [m]'
        iVar = iVar + 1
    endif
    if (nf%nlevs%daVarId > 0) then
        WRITE(iCtlUnit,"(A,X,I3,X,A)") 'da=>da ', nf%n_levs, ' Layer increment of A coefficients [Pa]'
        iVar = iVar + 1
    endif
    if (nf%nlevs%dbVarId > 0) then
        WRITE(iCtlUnit,"(A,X,I3,X,A)") 'db=>db ', nf%n_levs, ' Layer increment of B coefficients []'
        iVar = iVar + 1
    endif
    if (nf%nlevs%aVarId > 0) then
        WRITE(iCtlUnit,"(A,X,I3,X,A)") 'a=>a ', nf%n_levs, ' Layer mid A coefficients [Pa]'
        iVar = iVar + 1
    endif
    if (nf%nlevs%bVarId > 0) then
        WRITE(iCtlUnit,"(A,X,I3,X,A)") 'b=>b ', nf%n_levs, ' Layer mid B coefficients []'
        iVar = iVar + 1
    endif

    ! "Normal variables"
    do nVar = 1,nf%n_vars
      chTmp = nf%nvars(nVar)%chVarNm
      chTmp1 = ""
      if(defined(nf%nvars(nVar)%species)) chTmp1=fu_str(nf%nvars(nVar)%species)

      !Try to construct reasonable GRADS name
         if(len(trim(chTmp1))>0) then !Species exist
               chTmp =  trim(fu_quantity_short_string(nf%nvars(nVar)%quantity))
               jTmp  = min(len(chTmp),10)  !Cut quantity length
               iTmp = index(chTmp1,"POLLEN_") 
               if (iTmp>0) then !Get taxon name
                  chTmp1=chTmp1(iTmp+7:)
                  !Cut if needed short_name + taxon 
                  chTmp = chTmp(1:jTmp)+'_p'+chTmp1(1:3)
            elseif (index(chTmp1,"APHIDS")>0) then !!Aphids. Yes. Ugly solution for ugly problem
                   chTmp = chTmp(1:jTmp)+'_pAPH'  !! Pretend it is also pollen
               else !No pollen
                  iTmp = index(chTmp1,"NH415SO4")
                  if (iTmp>0) then
                      chTmp = chTmp(1:jTmp)+'_AmS'+chTmp1(9:) !_mode
                  else
                      ! By default -- cut quantity and hope for the best
                      chTmp = chTmp(1:jTmp)+'_'+chTmp1 
                  endif
               endif
         endif

      if (len(trim(chTmp))> 15) then 
         chTmp =  chTmp(1:15) !Jus cut it in any case
      endif

      
      if (nf%nvars(nVAr)%if3D) then
        WRITE(iCtlUnit,'(A60,X,I3,X,A,X,A,X,A)') nf%nvars(nVar)%chVarNm+'=>'+chTmp, nf%n_levs, " t,z,y,x ", &
               & fu_quantity_string(nf%nvars(nVar)%quantity), trim(chTmp1)
      else
        WRITE(iCtlUnit,'(A60,A,A,X,A,X,A)') nf%nvars(nVar)%chVarNm + '=>'+chtmp ,'   0   t,y,x ', &
               & fu_quantity_string(nf%nvars(nVar)%quantity), trim(chTmp1)
      endif
      iVar = iVar + 1
    enddo
    WRITE(iCtlUnit,'(A)') 'ENDVARS'

    !Write vars string
    CALL FSEEK(iCtlUnit, iOff,  0, iStat) !status 
    WRITE(iCtlUnit,'(A6,I3)') "VARS  ", iVar

    close(iCtlUnit)

  end subroutine write_ctl_file_4_netcdf

  !************************************************************************

  subroutine write_next_field_to_netcdf_file(iFile, field_id, grid_data)
    !
    ! Puts the field to an open netcdf-file, to which
    ! a unit-number netcdf_unit points. The field must have two parts:
    ! the identification (silja_field_id) and grid-data (real array).
    !
    ! Fills the nCtl structure when necessary
    ! 
   implicit none

   type (silja_field_id), intent(in) :: field_id
   INTEGER, INTENT(in) :: iFile
   REAL, DIMENSION(:), intent(inout) :: grid_data !!nf90mpi wants it so...

   integer, dimension(4)::iCount3d, iStart3d
   integer, dimension(3)::iCount2d, iStart2d
   integer :: varid, iLev, iStat, iTmp2, iVar, nx, ny, offx, offy, gnx, gny

   TYPE(netcdf_file),POINTER :: nf

   ! Buffering of 3D fields
   real, dimension(:), pointer, save :: buf3D => null()
   integer, save :: bufvarid = -1,  buflevel = 0, bufsize = 0
   integer(kind=4), dimension(1), target :: timeval
   character (len=*), parameter :: sub_name='write_next_field_to_netcdf_file '

   if (fu_fails(ifile > 0 .and. ifile <= size(nfile), 'Bad ifile: '//fu_str(ifile), sub_name)) return
   nf => nfile(iFile)

   if (nf%ifMPIIO) then
     call smpi_get_decomposition(nx, ny, offx, offy, gnx, gny)
   else
     offx=0
     offy=0
     nx = nf%ngrids(1)%nx
     ny = nf%ngrids(1)%ny
   endif

   if (nf%ifmpiio .and. nf%ncver== 4) then
     iTmp2 = nx*ny*nf%n_levs
     if ( bufsize < iTmp2 ) then
       if (associated(buf3D)) deallocate(buf3D)
       allocate(buf3D(iTmp2), stat=iStat)
       bufsize = iTmp2
     endif
   endif


   ! Find out timestep, write new time value if necessary
   if (fu_valid_time(field_id) /= nf%ntime%last_valid_time) then
!      call report(field_id)
     nf%ntime%last_valid_time = fu_valid_time(field_id)
     timeval(1) = nint(fu_sec8(nf%ntime%last_valid_time - nCtl%start_time))
     nf%n_times = nf%n_times + 1
     call msg("Making a new time: "//fu_str(fu_valid_time(field_id))//" No in file, Value", nf%n_times, timeval(1))
     !Stupid enough, but all ranks should write the dimension variable
     if (nf%ncver == 3  .and. nf%ifmpiio) then
#ifdef WITH_PNETCDF 
      
      nf%nf90mpifieldCnt =  nf%nf90mpifieldCnt + 1
      istat = nf90mpi_bput_var(nf%unit_bin, nf%ntime%tVarId, &
                        & timeval, &
                        &  nf%nf90mpiReq(nf%nf90mpifieldCnt),  int((/nf%n_times/),kind=8), int((/1/),kind=8)) 
      if (iStat /= 0) then 
        call set_error("nf90mpi_put_var time: "//nf90mpi_strerror(iStat), sub_name)
        return
      endif
      if (nf%nf90mpifieldCnt == nf%nf90mpiBUFfdlds ) then ! Extra for time variable
        call flush_nf90mpi_buf(nf, sub_name)
        if (fu_fails(.not. error,"Error set", sub_name)) return
      endif
#else
       call set_error("Compiled without PNETCDF", sub_name)
#endif
     else
        istat = NF90_PUT_VAR(nf%unit_bin, nf%ntime%tVarId, &
                        & timeval, (/nf%n_times/), (/1/)) 
       if (istat /= nf90_noerr)then
         call msg("NF90_PUT_VAR error:"+nf90_strerror(istat), iStat) 
         call set_error('Failed to write t dimension variable',sub_name)
         return
       endif
     endif
   endif 

   ! Find out variable id
   varid = int_missing
   do iVar = 1, nf%n_Vars
     if((fu_quantity(field_id) == nf%nvars(iVar)%quantity) .and. &
      & (fu_species(field_id) == nfile(ifile)%nvars(ivar)%species)) then
       varid = nf%nvars(iVar)%varid
       exit
     endif
   enddo
   if (varid == int_missing)then
      call set_error('Quantities dont match',sub_name)
      call report(field_id)
      call msg("Variables in the file (quantity, species)")
           do iVar = 1, nf%n_Vars
             call msg( fu_str(ivar)+":" & 
                & +fu_quantity_string(nf%nvars(iVar)%quantity))
             call report(nfile(ifile)%nvars(ivar)%species) 
           enddo
      call msg("")
      return 
   endif

   ! Write the data. If only part of the data is wanted, define start/count/stride
   ! for nf90_put_var call

   if(nf%nvars(iVar)%if3D) then
     ! Find out level
     do iLev = 1, nf%n_levs !Dirty hack: NetCDF should be able to output proper layers
       if (fu_level_value(fu_central_level_of_layer(fu_level(field_id))) .eps. nf%nlevs%levels(iLev)) exit
     enddo

     if (iLev > nf%n_levs)then
       call msg('Field level value:', fu_level_value(fu_level(field_id)))
       call report(field_id)
       do iTmp2 = 1, nf%n_levs 
         call msg('Output level: ', iTmp2, fu_level_value(fu_level(nf%silamVertical, iTmp2)))
       enddo
       call set_error(fu_connect_strings('Levels dont match for: ',nf%nvars(iTmp2)%chVarNm),sub_name) 
       return
     endif

     if (ifMsgs) call msg('3D: variable:' + fu_quantity_string(nf%nvars(iVar)%quantity) + '_' + fu_str(nf%nvars(iVar)%species), iLev)

     if (nf%ifmpiio .and. nf%ncver== 4) then
       if (bufvarid /= varid .and. (buflevel /= 0 .or. iLev /= 1)) then !New var, dirty buffer
         call msg("buflevel, iLev, bufvarid, varid", (/buflevel, iLev, bufvarid, varid/))
         call set_error("buffer must be clean for the next 3D field",sub_name)
       elseif (iLev-1 /= buflevel) then 
         call msg("buflevel, iLev, bufvarid, varid", (/buflevel, iLev, bufvarid, varid/))
         call set_error("buffer must have exactly iLev-1 levels", sub_name)
       else ! continue accumulation 
         bufvarid = varid
         iTmp2 = nx*ny
         buf3D(iTmp2*buflevel+1:iTmp2*(buflevel+1)) = grid_data(1:iTmp2)
         buflevel = iLev
         if (iLev == nf%n_levs) then !Dump stuff to netcdf file
           iStart3d = (/1+offx, 1+offy, 1, nf%n_Times/)
           iCount3d = (/nx, ny, nf%n_levs, 1/)
           iStat = NF90_PUT_VAR(nf%unit_bin, varid, buf3D, istart3d, icount3d)
           if (iStat /= 0)then 
             call set_error('Failed to write 3d variable:'+nf%nvars(iVar)%chVarNm,sub_name)
             return
           endif
           buflevel = 0 !Reset buffer
           bufvarid = -1
  !         call msg("write_next_field_to_netcdf_file 3D: "+nf%nvars(iVar)%chVarNm)
         endif
       endif
     else ! No buffer, plain 2d fields
       iStart3d = (/1+offx, 1+offy, iLev, nf%n_Times/)
       iCount3d = (/nx, ny, 1, 1/)
       call PUT_VAR_nc_buf(nf, varid, grid_data, istart3d, icount3d)
       if (error)then 
         call set_error('Failed to write 3d variable:'+nf%nvars(iVar)%chVarNm,sub_name)
         return
       endif
     endif
   else  ! 2D variable -- do not care about levels
     if (ifMsgs) call msg('2D: variable:' + fu_quantity_string(nf%nvars(iVar)%quantity) + '_' + fu_str(nf%nvars(iVar)%species))
     iStart2d = (/1+offx, 1+offy, nf%n_Times/)
     iCount2d = (/nx, ny, 1/)
     call PUT_VAR_nc_buf(nf, varid, grid_data, istart2d, icount2d)
     if (error)then 
       call set_error('Failed to write 2d variable:  '+nf%nvars(iVar)%chVarNm,sub_name)
       return
     endif
   endif

   if (error) then
       call msg("NF90_PUT_VAR error:"+nf90_strerror(iStat), iStat)
       call msg("iStart2d", iStart2d)
       call msg("iCount2d", iCount2d)
       call msg("nf%ngrids(1)%nx nf%ngrids(1)%ny", nf%ngrids(1)%nx, nf%ngrids(1)%ny)
       call msg('')
   endif

   !------------------------------------------------------------------------------------------
   ! update in the nCtl structure
   if(.not. nCtl%defined) then
     call set_error('nCtl structure must be defined by now',sub_name)
   endif
   if(nCtl%last_valid_time /= nf%ntime%last_valid_time)then
     nCtl%last_valid_time = nf%ntime%last_valid_time
     nCtl%nTimes = nCtl%nTimes + 1
   endif
   ! if second timestep, write the step to nctl structure
   if(nCtl%step == interval_missing) then
     if(nCtl%nTimes == 2)then 
        nCtl%step = nCtl%last_valid_time - nCtl%start_time
     elseif(nCtl%nTimes > 2)then
        call set_error('Failed to fill the nCtl structure',sub_name)
     endif
   endif
  END SUBROUTINE write_next_field_to_netcdf_file


 !***************************************************************

  integer function fu_next_free_netcdf_structure()
    !
    ! Finds the first free NETCDF structure and returns its index
    !
    implicit none

    ! Local variables
    integer :: i

    fu_next_free_netcdf_structure = int_missing
    !$OMP critical (next_free_netcdf)
    do i=1, max_nbr_of_netcdf_files-1
      if(nfile(i)%defined == silja_false)then
        if(nfile(i)%unit_bin == -1)then
          fu_next_free_netcdf_structure = i
          exit
        endif
      endif
    end do
    !$OMP end critical(next_free_netcdf)
    if (fu_next_free_netcdf_structure == int_missing) &
                     & call set_error('No free structures left','fu_next_free_netcdf_structure')

  end function fu_next_free_netcdf_structure 


  !*****************************************************************


  subroutine init_netcdf_io(nlStdSetup, nlOut)
    !
    ! Just allocates the NETCDF file structure
    !
    implicit none

    ! Imported variables
    type(Tsilam_namelist), pointer :: nlStdSetup, nlOut


    ! Local variables
    integer :: iStat, iTmp

    ! Read the name table file. 
    iTmp = fu_next_free_unit()
    open(file=fu_expand_environment(fu_content(nlStdSetup,'netcdf_name_table_fnm')), &
      & unit=iTmp, action='read', status='old', iostat=iStat)
    if(iStat /= 0)then
      call set_error(fu_connect_strings('Failed to open netcdf name table:', &
                   & fu_content(nlStdSetup,'netcdf_name_table_fnm')), 'init_netcdf_io')
      return
    endif

    nmTblNlGrp => fu_read_namelist_group(iTmp, .false., 'END_NETCDF_NAME_TABLE')
    close(iTmp)
    IF (error) RETURN

    if(.not. associated(nmTblNlGrp))then
      call msg_warning('NetCDF coding namelist group is not associated','init_netcdf_io')
    endif

    if(len_trim(fu_content(nlOut,'netcdf_preferred_mass_unit')) > 0)then
      call set_error('netcdf_preferred_mass_unit has been disabled','init_netcdf_io')
      return
    endif

  end subroutine init_netcdf_io


!***************************************************************************

  subroutine id_list_from_netcdf_file(iUnit, idList, iCount)
    !
    ! Get ID list
    ! Always all levels, first time

    implicit none
    integer, intent(in) :: iUnit
    integer, intent(out) :: iCount
  
    type(silja_field_id),dimension(:), pointer :: idList
    type(netcdf_file),pointer :: nf
    type(silja_interval) :: length_of_accumulation, fc_length
    integer :: iVar, iLev,  field_kind, nFields, ntimes, nlevs
    type(silja_level) :: levTmp

    if(.not. nfile(iUnit)%defined == silja_true )then
      call set_error('nf ptr not defined','id_list_from_netcdf_file')
      return
    endif
    nf => nfile(iUnit)

    
    !! Count fields
    nFields = 0
    do iVar = 1, nf%n_vars
      if(nf%nVars(iVar)%quantity == int_missing)cycle
      nFields = nFields + nf%nVars(iVar)%n_levs
    enddo


    iCount = 0
    allocate(idList(nFields), stat = iVar)
    if(iVar /= 0)then
      call set_error('Failed to allocate id list', 'id_list_from_netcdf_file')
    endif

    do iVar = 1, nf%n_vars
      if(nf%nVars(iVar)%quantity == int_missing)cycle
      nlevs = nf%nVars(iVar)%n_levs
      do iLev = 1, nlevs
         iCount = iCount + 1
         call id_for_nc_index(iUnit, idList(iCount), iVar, 1, iLev)
        !! call report(idList(iCount) )
      enddo   ! iLev
    enddo  ! iVar
  end subroutine id_list_from_netcdf_file


  subroutine id_for_nc_index(iUnit, id_out, iVar, iT, iLev)

    !
    ! This supposed to be the only place that is supposed to generate field_id
    ! from netcdf file
    !  Not much checks
    !
    ! Loose copy from id_list_from_netcdf_file
    !

    implicit none
    integer, intent(in) :: iUnit
    type(silja_field_id), intent(out) :: id_out
    integer, intent(in) ::  iVar, iT, iLev
  
    type(netcdf_file),pointer :: nf
    type(silja_interval) :: length_of_accumulation, fc_length
    integer ::  field_kind, nFields, ntimes, nlevs
    type(silja_level) :: levTmp
     character (len=*), parameter :: sub_name = "id_for_nc_index"

    if(.not. nfile(iUnit)%defined == silja_true )then
      call set_error('nf ptr not defined',sub_name)
      return
    endif
    nf => nfile(iUnit)

    if(fu_accumulated_quantity(nf%nVars(iVar)%quantity))then
      field_kind = accumulated_flag
      length_of_accumulation = nf%nTime%times(iT) - nf%nTime%analysis_time
      fc_length = length_of_accumulation
    else
      if (nf%nTime%time_label_position == instant_fields) then
         field_kind = forecast_flag
         length_of_accumulation = zero_interval
         fc_length = nf%nTime%times(iT) - nf%nTime%analysis_time
      else
         length_of_accumulation =  nf%nTime%step
         field_kind = averaged_flag

         select case (nf%nTime%time_label_position)
            case (end_of_period) 
               fc_length = nf%nTime%times(iT) - nf%nTime%analysis_time
            case (mid_of_period)
               fc_length = nf%nTime%times(iT) - nf%nTime%analysis_time + nf%nTime%step*0.5 
            case (start_of_period)
               fc_length = nf%nTime%times(iT) - nf%nTime%analysis_time + nf%nTime%step

            case default
               call set_error("nf%nTime%time_label_position="//fu_str(nf%nTime%time_label_position)&
                  &//"  handling not implemented", sub_name)
               return
         end select
       endif
    endif

    !
    ! The variable may not have any meaningful vertical: deposition, for isntance
    !
    if(defined(nf%nVars(iVar)%sVert))then
      levTmp = fu_level(nf%nVars(iVar)%sVert, iLev)
    else
      levTmp = surface_level  !!!!Level from nametable is ignored FIXME
    endif

    id_out = fu_set_field_id(met_src_missing, &
                                   & nf%nVars(iVar)%quantity, &
                                   & nf%nTime%analysis_time,&
                                   & fc_length, &
                                   & nf%nVars(iVar)%nGrid%sGrid,&
                                   & levTmp,&
                                   & length_of_accumulation, & ! optional
                                   & zero_interval, &     ! ## length_of_validity,
                                   & field_kind, &             ! optional
                                   & species = nf%nvars(ivar)%species, &
                                   & chCocktail = nf%nvars(ivar)%chCocktailNm)



  end subroutine id_for_nc_index



  !**************************************************************************

  subroutine get_netcdf_grids(incf, quantity, species, grids, nGrids)
    !
    ! Finds out the grids of the given quantity and species - if defined.
    ! Otherwise provides all in the file
    !
    implicit none
    
    ! Imported parameters
    integer, intent(in) :: incf, quantity
    type(silam_species), intent(in) :: species
    type(silja_grid), dimension(:), pointer :: grids
    integer, intent(out) :: nGrids

    ! Local parameters
    type(netcdf_file),pointer :: nf
    integer :: iVar
    !
    ! Some preparations first
    !
    if(.not. nfile(incf)%defined == silja_true )then
      call set_error('nf ptr not associated','get_netcdf_grids')
      return
    endif
    nf => nfile(incf)
    !
    ! Scan the list of variables twice (SORRY!) in order to understand how many grids we have
    !
    nGrids = 0
    do iVar = 1, nf%n_vars
      !
      ! filtering unnecesary stuff
      !
      if(nf%nVars(iVar)%quantity == int_missing)cycle
      if(quantity /= int_missing)then
        if(quantity /= nf%nVars(iVar)%quantity)cycle
      endif
      if(defined(species))then
        if(.not. (species == nf%nVars(iVar)%species))cycle
      endif
      nGrids = nGrids + 1
    enddo
    if(nGrids == 0)then
      call msg_warning('No grids found for the given quantity:' + fu_quantity_string(quantity) + &
                     & ', and species:' + fu_str(species),'get_netcdf_grids')
      nullify(grids)
      return
    endif
    !
    ! Store the grids to the given array
    !
    allocate(grids(nGrids), stat=iVar)
    if(fu_fails(iVar==0,'Failed grid array allocation','get_netcdf_grids'))return
    nGrids = 0
    do iVar = 1, nf%n_vars
       if(nf%nVars(iVar)%quantity == int_missing)cycle
       if(quantity /= int_missing)then
         if(quantity /= nf%nVars(iVar)%quantity)cycle
       endif
       if(defined(species))then
        if(.not. (species == nf%nVars(iVar)%species))cycle
       endif
       nGrids = nGrids + 1
       grids(nGrids) = nf%nVars(iVar)%nGrid%sGrid
    enddo

  end subroutine get_netcdf_grids




  !**************************************************************************

  subroutine get_netcdf_verticals(incf, quantity, species, verticals, nVerts)
    !
    ! Finds out the grids of the given quantity and species - if defined.
    ! Otherwise provides all in the file
    !
    implicit none
    
    ! Imported parameters
    integer, intent(in) :: incf, quantity
    type(silam_species), intent(in) :: species
    type(silam_vertical), dimension(:), pointer :: verticals
    integer, intent(out) :: nVerts

    ! Local parameters
    type(netcdf_file),pointer :: nf
    integer :: iVar
    !
    ! Some preparations first
    !
    if(.not. nfile(incf)%defined == silja_true )then
      call set_error('nf ptr not associated','get_netcdf_verticals')
      return
    endif
    nf => nfile(incf)
    !
    ! Scan the list of variables twice (SORRY!) in order to understand how many grids we have
    !
    nVerts = 0
    do iVar = 1, nf%n_vars
      !
      ! filtering unnecesary stuff
      !
      if(nf%nVars(iVar)%quantity == int_missing)cycle
      if(quantity /= int_missing)then
        if(quantity /= nf%nVars(iVar)%quantity)cycle
      endif
      if(defined(species))then
        if(.not. (species == nf%nVars(iVar)%species))cycle
      endif
      nVerts = nVerts + 1
    enddo
    if(nVerts == 0)then
      call msg_warning('No verticals found for the given quantity:' + fu_quantity_string(quantity) + &
                     & ', and species:' + fu_str(species),'get_netcdf_verticals')
      nullify(verticals)
      return
    endif
    !
    ! Store the grids to the given array
    !
    allocate(verticals(nVerts), stat=iVar)
    if(fu_fails(iVar==0,'Failed grid array allocation','get_netcdf_verticals'))return
    nVerts = 0
    do iVar = 1, nf%n_vars
       if(nf%nVars(iVar)%quantity == int_missing)cycle
       if(quantity /= int_missing)then
         if(quantity /= nf%nVars(iVar)%quantity)cycle
       endif
       if(defined(species))then
        if(.not. (species == nf%nVars(iVar)%species))cycle
       endif
       nVerts = nVerts + 1
       verticals(nVerts) = nf%nVars(iVar)%sVert
    enddo

  end subroutine get_netcdf_verticals


  !***************************************************************************
 
 
  function fu_netcdf_filename(iFile)

    implicit none

    integer, intent(in) :: iFile
    CHARACTER (LEN=fnlen) :: fu_netcdf_filename
    
    fu_netcdf_filename = nfile(iFile)%fname

  end function fu_netcdf_filename
  
  
!***************************************************************************
 
  subroutine timelst_from_netcdf_file(input_unit, timeLst, nTimes, ifCopy_, ifEnvelope_, an_time)
    !
    ! Returns the times of NetCDF. Can copy them to time list or just give a link,
    ! can also return an envelope of the field time ranges, which is needed for 
    ! averaged and/or period-valid fields. Envelope starts at the beginning of the first
    ! time range and ends at the end of the last one.
    !
    implicit none

    integer, intent(in) :: input_unit
    integer, intent(out) :: nTimes
    logical, intent(in), optional :: ifCopy_, ifEnvelope_
    type(silja_time), intent(out), optional :: an_time 
    type(silja_time), dimension(:), pointer :: timeLst
    
    logical :: ifCopy, ifEnvelope
    integer :: iTmp

    if(present(ifCopy_))then
      ifCopy = ifCopy_
    else
      ifCopy = .false.
    endif
    if(present(ifEnvelope_))then
      ifEnvelope = ifEnvelope_
    else
      ifEnvelope = .false.
    endif


    if (present(an_time)) an_time =  nfile(input_unit)%ntime%analysis_time

    if(ifCopy)then      ! Copy the time array
      if(ifEnvelope)then    ! times should cover the range of validity of the fields
        if(fu_fails(defined(nfile(input_unit)%nTime%step), &
                  & 'Time step must be defined to ask for envelope time array', &
                  & 'timelst_from_netcdf_file'))return
        nTimes = nfile(input_unit)%n_times + 1
        allocate(timeLst(nTimes), stat = iTmp)
        if(fu_fails(iTmp==0,'Failed time array allocation','timelst_from_netcdf_file'))return
        !
        ! Envelope is arranged differently depending on the time label position
        !
        select case(nfile(input_unit)%ntime%time_label_position)
          case(start_of_period)
            do iTmp = 1, nTimes-1
              timeLst(iTmp) = nfile(input_unit)%ntime%times(iTmp)
            end do
          case(mid_of_period)
            do iTmp = 1, nTimes-1
              timeLst(iTmp) = nfile(input_unit)%ntime%times(iTmp) - &
                            & nfile(input_unit)%nTime%step * 0.5
            end do
          case(end_of_period)
            do iTmp = 1, nTimes-1
              timeLst(iTmp) = nfile(input_unit)%ntime%times(iTmp) - &
                            & nfile(input_unit)%nTime%step
            end do
          case default
            call set_error('Unknown time label position:' + &
                         & fu_str(nfile(input_unit)%ntime%time_label_position), &
                         & 'timelst_from_netcdf_file')
            return
        end select
        timeLst(nTimes) = timeLst(nTimes-1) + nfile(input_unit)%nTime%step
      else   ! not envelope times, just copy the array
        nTimes = nfile(input_unit)%n_times
        allocate(timeLst(nTimes), stat = iTmp)
        if(fu_fails(iTmp==0,'Failed time array allocation','timelst_from_netcdf_file'))return
        do iTmp = 1, nTimes
          timeLst(iTmp) = nfile(input_unit)%ntime%times(iTmp)
        end do
      endif  ! if envelope time array
    else   ! no copy, just set pointer
      if(ifEnvelope)then
        call set_error('Cannot make envelope without copying time array','timelst_from_netcdf_file')
        nullify(timeLst)
        nTimes = int_missing
        return
      else
        nTimes = nfile(input_unit)%n_times
        timeLst => nfile(input_unit)%ntime%times
      endif   ! ifEnvelope time array
    endif  ! if copy times
  end subroutine timelst_from_netcdf_file

!***************************************************************************
  

  subroutine field_indices_from_netcdf_file(iUnit, id_req, id_found, iVar, it, ilev)
    !
    ! Tries to figure out indices that cover given id and set them according to the 
    ! times in the file.
    ! One should request instant, accumulated or 

    !
    implicit none

    integer, intent(in) :: iUnit
    type(silja_field_id), intent(in) :: id_req
    type(silja_field_id), intent(out) :: id_found
    integer, intent(out) :: it, iVar, ilev
    
    

    type(netcdf_file),pointer :: nf
    integer :: label_pos, ntimes, itmp
    type(silja_time) :: valtime_req, timetmp, valtime_nc
    type(silja_interval) :: valint_req, tstep, valint_nc, offset  !!!netcdf_time = valtime + offset
     character (len=*), parameter ::sub_name = "field_indices_from_netcdf_file"

    id_found = field_id_missing
    nf => nfile(iUnit)
    it = int_missing
    nTimes = nfile(iUnit)%n_times


    ! First find out the variable id
    do iVar = 1,nf%n_vars
      if(fu_quantity(id_req) /= nf%nvars(iVar)%quantity)cycle
      if (.not. fu_species(id_req) == nf%nvars(ivar)%species) cycle
      if (.not. fu_cocktail_name(id_req) == nf%nvars(ivar)%chCocktailNm) cycle
      exit      
    end do  ! cycle over variables

    if(iVar > nf%n_vars)then
      call msg_warning('Requested variable not found',sub_name)
      call report(id_req)
      call set_error('Variabl not found',sub_name)
      return
    endif    

    ! Check if the requested level exists
    if(nf%nvars(iVar)%n_levs > 1)then
      do iLev = 1,nf%nvars(iVar)%n_levs
        if(fu_cmp_levs_eq(fu_level(id_req),fu_level(nf%nVars(iVar)%sVert,iLev))) exit
      enddo    ! cycle over levels for the multi-level variable

      if( iLev > nf%nvars(iVar)%n_levs)then
        call msg_warning('Level not found',sub_name)
        call report(id_req)
        call set_error('Level not found',sub_name)
        return
      endif
    else
      iLev = 1
    end if  ! if single-level variable

    if (defined(fu_grid(id_req))) then
      if (.not. (nf%nvars(iVar)%nGrid%sGrid == fu_grid(id_req))) then
        call msg_warning("Grids missmatch", sub_name)
        call msg("requested grid")
        call report(fu_grid(id_req))
        call msg("available grid")
        call report(nf%nvars(iVar)%nGrid%sGrid)
        call set_error("Grids missmatch", sub_name)
      endif
    endif


    !
    ! Times
    ! 
    if( fu_sec(fu_validity_length(id_req)) > 0 ) then
      call msg("Requested ID")
      call report(id_req)
      call set_error("No period-valid fields from NetCDF", sub_name)
    endif



    label_pos = nf%ntime%time_label_position
    select case(label_pos)
      case(end_of_period, instant_fields)
         offset = zero_interval
       case(mid_of_period)
         offset = nf%nTime%step * 0.5
       case(start_of_period)
         offset = nf%nTime%step
       case default
         call set_error('Unknown time label position:' + &
                      & fu_str(nf%ntime%time_label_position), &
                      & sub_name)
         return
    end select

    timetmp = fu_valid_time(id_req)
    if (fu_field_kind(id_req) == averaged_flag ) then
       timetmp = timetmp - offset
       if (.not.  fu_accumulation_length(id_req) == nf%nTime%step ) then
             call msg("NCfile timestep: "//trim(fu_str(nf%nTime%step)) )
             call msg("id_requested:")
             call report(id_req)
             call set_error("Averaging time missmatch",sub_name)
             return
       endif
    endif

    !Just find matching time
    do iT = 1,nTimes
      if (timetmp ==  nf%ntime%times(iT)) exit
    enddo
    if  (iT > nTimes) then
      call set_error("No matching timestamp" ,sub_name)
      return
    endif

    call id_for_nc_index(iUnit, id_found, iVar, iT, iLev)

  end subroutine field_indices_from_netcdf_file




!***************************************************************************

  subroutine make_staggered_fld(refFld, nxRef, nyRef, stagFld, nxStag, nyStag, stagDirX, stagDirY)

    !
    ! So here we fix 0 lons and lats for staggered wrf grids, interpolating/extrapolating the reference one.
    ! number of gridcells for staggered field has to be either one bigger or one smaller than the ref one
    ! or the staggering direction has to be defined
    !
    
    implicit none
    integer, intent(in) :: nxRef, nyRef, nxStag, nyStag
    real, dimension(:), intent(in):: refFld
    real, dimension(:), intent(inout):: stagFld
    integer, intent(in), optional :: stagDirX, stagDirY
    integer :: ix, iy



    if(abs(nxRef-nxStag)>1)then
      call set_error('nx too different','make_staggered_fld')
      return
    endif
    if(abs(nyRef-nyStag)>1)then
      call set_error('ny too different','make_staggered_fld')
      return
    endif
    if(nxRef==nxStag .and. nyRef==nyStag .and. .not. (present(stagDirX) .or. present(stagDirY)))then
      call set_error('confused','make_staggered_fld')
      return
    endif

    if(nxStag-nxRef==1 .and. nyStag==nyRef)then ! staggered grid bigger by 1 cell in x dir

      do iy = 1, nyStag
        do ix = 1, nxStag
          if(ix == 1)then
            stagFld((iy-1)*nxStag+1) =  refFld((iy-1)*nxRef+1) - 0.5*(refFld((iy-1)*nxRef+2) - refFld((iy-1)*nxRef+1))
          elseif(ix == nxStag)then
            stagFld(iy*nxStag) = refFld(iy*nxRef) + 0.5*(refFld(iy*nxRef) - refFld(iy*nxRef-1))
          else
            stagFld((iy-1)*nxStag+ix) = 0.5 * (refFld((iy-1)*nxRef+ix-1) + refFld((iy-1)*nxRef+ix))
          endif
        enddo
      enddo

      return
    endif

    if(nyStag-nyRef==1 .and. nxStag==nxRef)then ! staggered grid bigger by 1 cell in y dir

      do iy = 1, nyStag
        do ix = 1, nxStag
          if(iy == 1)then
            stagFld(ix) =  refFld(ix) - 0.5*(refFld(nxRef+ix) - refFld(ix))
          elseif(iy == nyStag)then
            stagFld((nyStag-1)*nxStag+ix) &
                 & = refFld((nyRef-1)*nxRef+ix) + 0.5*(refFld((nyRef-1)*nxRef+ix) - refFld((nyRef-2)*nxRef+ix))
          else
            stagFld((iy-1)*nxStag+ix) = 0.5 * (refFld((iy-2)*nxRef+ix) + refFld((iy-1)*nxRef+ix))
          endif
        enddo
      enddo

      return
    endif

    call set_error('not implemented','make_staggered_fld')

  end subroutine make_staggered_fld


  !****************************************************************************************

  subroutine set_missing_netcdf_dim(ncdim)
    implicit none
    type(netcdf_dim), intent(inout) :: ncdim
    type(netcdf_dim) :: blank_dim
   
    if (associated(ncdim%values)) deallocate(ncdim%values)
    if (associated(ncdim%a_half)) deallocate(ncdim%a_half)
    if (associated(ncdim%b_half)) deallocate(ncdim%b_half)
    if (associated(ncdim%a)) deallocate(ncdim%a)
    if (associated(ncdim%b)) deallocate(ncdim%b)
    if (associated(ncdim%values2d)) deallocate(ncdim%values2d)
    call set_missing(ncdim%svert, .false.)
    
    ! The variable blank_dim has default initialization, which we copy...
    ncdim = blank_dim
    ! ...but still need to nullify since this is not in the default initialization.
    nullify(ncdim%values, ncdim%a, ncdim%b,ncdim%a_half,ncdim%b_half, ncdim%values2d)
  end subroutine set_missing_netcdf_dim


  !*****************************************************************

  function fu_netcdf_str_to_silam_time(chTime) result(time)
    !
    ! Converts the netcdf-type string to the SILAM type structure
    ! The input template is: yyyy-mm-dd hh:mm:ss
    ! For example 0000-01-01 00:00:00
    !
    implicit none

    ! Return value
    type(silja_time) :: time

    ! Imported parameter
    character(len=*), intent(in) :: chTime
    character :: chtmp

    ! Local variables
    integer :: year, month, day, hour, minute, iStatus
    integer :: iywidth,imwidth,idwidth, i, ind
    real :: second
    character(len=80) :: chF
    character(len=1), dimension(2) :: delimiters
    year=int_missing; month=int_missing; day=int_missing 

    second = 0.
    minute = 0
    hour = 0
    time = time_missing

    delimiters = (/'T', '_' /)

    iywidth=index(chTime,'-')-1
    imwidth=index(chTime(iywidth+2:),'-') -1
    idwidth=index(chTime(iywidth+imwidth+3:),' ') -1
    if (all(idwidth /= (/1,2/)) ) then !Try to look for T -- delimiter in thredds dates
      do i=1, size(delimiters)
        ind = index(chTime(iywidth+imwidth+3:), delimiters(i))
        if (ind > 0) then
          idwidth = ind - 1
          exit
        end if
      end do
        

      if (idwidth /= 2) then
        call set_error('Failed to get time from string1"'//trim(chTime)//'"' , &
             & 'fu_netcdf_str_to_silam_time')
        return
      endif
    endif
    chF='(I'+fu_str(iywidth)+', A, I'+fu_str(imwidth)+', A, I'+fu_str(idwidth)+', A, I2, A, I2, A, F2.0)'



      !
    read(unit=chTime,fmt=chF, iostat=iStatus) &
                       & year, chtmp, month, chtmp, day, chtmp, hour, chtmp, minute, chtmp, second
    if (iStatus /= 0)then
      call set_error('Failed to get time from netcdf  string "'//trim(chTime)//'"' , &
                  & 'fu_netcdf_str_to_silam_time')
      return
    endif

    time = fu_set_time_utc(year, month, day, hour, minute, second)
    
  end function fu_netcdf_str_to_silam_time


  !***********************************************************************************
  
  subroutine get_netcdf_total(incf, flds_3D_Requested, start, duration, layer, pOutput)
    !
    ! Searches for the given quantity and species, if defined, limits time period and 
    ! gets the grand totals for all species (or just one, if it is defined)
    !
    implicit none
    
    ! Imported parameters
    integer, intent(in) :: incf
    type(silja_field_3d_pointer), dimension(:), pointer :: flds_3D_Requested
    type(silja_time), intent(in) :: start
    type(silja_interval), intent(in) :: duration
    type(silja_level), intent(in) :: layer
    type(silja_rp_1d), dimension(:), pointer :: pOutput
    
    ! Local parameters
    type(netcdf_file),pointer :: nf
    integer :: iVar, iTime, iFldRequested, iFldQuantity, iLev, fs, field_kind, nTimes
    type(silam_species) :: FldSpecies
    character(len=substNmLen) :: chFldCocktailNm
    real, dimension(:), pointer :: pData
    type(silja_field_id) :: idTmp
    type(silja_interval) :: length_of_accumulation
    real :: duration_sec
    type(silja_time) :: startTmp, endTmp
    type(silja_time), dimension(:), pointer :: envelope_times
    !
    ! Some preparations first
    !
    if(.not. nfile(incf)%defined == silja_true )then
      call set_error('nf ptr not associated','get_netcdf_total')
      return
    endif
    nf => nfile(incf)
    do iVar = 1, size(flds_3d_Requested)
      pOutput(iVar)%pp(:) = 0.0
    end do

    call timelst_from_netcdf_file(incf, envelope_times, nTimes, .true., .true.) !ifCopy_, ifEnvelope_
    if(error)return
    !
    ! Store the required 
    !
    do iTime = 1, nTimes-1
      !
      ! Scan the list of variables, set the required ID and read the stuff
      !
      if(envelope_times(iTime+1) <= start .or. envelope_times(iTime) >= start + duration)cycle
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

        do iVar = 1, nf%n_vars
          !
          ! Needed variable? It is defined by quantity and species / cocktail_name
          !
          if(nf%nVars(iVar)%quantity == int_missing)cycle
          if(iFldQuantity /= int_missing)then
            if(iFldQuantity /= nf%nVars(iVar)%quantity)cycle
          endif
          if(defined(FldSpecies))then
            if(.not. (FldSpecies == nf%nVars(iVar)%species))cycle
          endif
          if(len_trim(chFldCocktailNm) > 0)then
            if(.not. (chFldCocktailNm == nf%nVars(iVar)%chCocktailNm))cycle
          endif
          !
          ! OK, variable found that meets the requirements of this field_3d
          ! Go layer by layer summing-up the maps
          !
          if(fu_accumulated_quantity(nf%nVars(iVar)%quantity))then
            field_kind = accumulated_flag
            length_of_accumulation = nf%nTime%times(iTime) - nf%nTime%analysis_time
          else
            field_kind = forecast_flag
            length_of_accumulation = zero_interval 
          endif

          do iLev = 1, nf%nVars(iVar)%n_levs
 
            idTmp = fu_set_field_id(met_src_missing, &
                                  & nf%nVars(iVar)%quantity, &
                                  & nf%nTime%analysis_time,&
                                  & nf%nTime%times(iTime) - nf%nTime%analysis_time, &
                                  & nf%nVars(iVar)%nGrid%sGrid,&
                                  & fu_level(nf%nVars(iVar)%sVert, iLev),&
                                  & length_of_accumulation, & ! optional
                                  & zero_interval, &     !optional
                                  & field_kind, &             ! optional
                                  & species = nf%nvars(ivar)%species, &
                                  & chCocktail = nf%nVars(iVar)%chCocktailNm)
            if(error)return

            call read_field_from_netcdf_file(incf, idTmp, pData, 0.0)
            if(error)return

            pOutput(iFldRequested)%pp(1:fs) = pOutput(iFldRequested)%pp(1:fs) + &
                                            & pData(1:fs) * duration_sec
          enddo   ! iLev
        enddo  ! iVar
      enddo  ! iFldRequested
    enddo  ! nf%n_times

  end subroutine get_netcdf_total


  !****************************************************************

  subroutine expand_output_list(Lst, nNewItems)
    !
    ! Enlarges the list of output variables by nNewItems items
    !
    implicit none

    type(TOutputList), intent(inout) :: Lst
    integer, intent(in) :: nNewItems

    ! Local variables
    type(TOutputList) :: LstTmp
    integer :: iItem, nItems

    !
    ! Steps: allocate temporary; copy the existing stuff to safety, reallocate bigger main list, 
    ! copy the stuff back
    !
    if(allocated(Lst%PtrItem))then
      !
      ! Have to save the stuff to temporary and reallocate
      !
      nItems = size(Lst%PtrItem)

      allocate(LstTmp%ptrItem(nItems+nNewItems), stat = iItem)
      if(fu_fails(iItem == 0, 'Failed to allocate temporary','expand_output_list'))return

      LstTmp%ptrItem(1:nItems) = Lst%ptrItem(1:nItems)
      
      !Does it not happen automatically?
      LstTmp%ptrItem(nItems+1:nItems+nNewItems) =  OutputLstItem_missing

      call move_alloc(LstTmp%ptrItem, Lst%ptrItem)

    else
      !
      ! No problem at all - just allocate the space for a couple of variables
      !
      allocate(Lst%ptrItem(nNewItems), stat = iItem)  ! plus new variables
      if(fu_fails(iItem == 0, 'Failed to allocarte temporary','expand_output_list'))return

      !Does it not happen automatically?
      Lst%ptrItem(1:nNewItems) = OutputLstItem_missing         ! kill the new variables
    endif

  end subroutine expand_output_list

  


END MODULE netcdf_io


