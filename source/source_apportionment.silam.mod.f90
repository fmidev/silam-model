module source_apportionment
  
  ! This module contains routines for handling the source term in data
  ! assimilation. This is implemented through "emission processors"
  ! which is intended to offer a generic interface for any
  ! adjustements done to the emission mass map before adding it to the tranport map.
  !
  ! In forward integrations the emission processor corresponds to application of the
  ! control variable into the model, ie. scaling or forcing the emission flux. In adjoint
  ! integration (4D-Var) sources become receptors: the processor accumulates values from
  ! the adjoint variable (with weighting and/or masking), and in the end of integration,
  ! the gradient with respect to the control variable can be collected from the processor.
  !
  ! The emission adjustment can be either
  ! - 2D map of scaling factors for each species, applied to every level
  ! - 2D scaling factors resolving time and vertical, applied to every gridcell
  ! - 2D forcings resolving time and vertical, applied to every gridcell with nonzero emission,
  !   replacing the existing emission map.
  ! - 2D forcings as above but weighted with square root of level thickness.
  ! 
  ! The weighted 2D forcing aims to make the background term in the cost function
  ! independent of the vertical discretization - even if the background term is not
  ! explicitly included, but provided by truncated iteration.
  ! 
  ! Let E_i be the emission (kg/s) to ith layer. Then define the control variable as
  ! 
  ! E_i' = E_i / sqrt(dz_i), 
  ! 
  ! where dz_i is the layer thickness.
  ! 
  ! Now 
  ! 
  ! sum(E_i'^2) = sum(E_i^2 / dz_i) = sum(rho^2*dz_i), 
  ! 
  ! where rho is the emission density in kg / (m*s). This way, the 
  ! squared sum, as used in the cost function, becomes only
  ! dependent on the (continuous) distribution of rho, not dz.

  ! Using the weighted mode generally recommended. Without it, the inversion will favour
  ! putting the emission into thinner levels.

  use dispersion_server
  implicit none
  private

  public defined
  public fu_data
  public set_emission_processor
  public process_emission_forward
  public process_emission_adjoint
  public zero_processor_grad
  public fu_type
  public fu_grid
  public fu_vertical
  public set_da_source_indices
  public fu_species_list
  public report

  private fu_defined_processor
  private fu_type_processor
  private fu_vertical_proc
  private fu_grid_proc
  private report_emission_processor
  
  interface set_emission_processor
     module procedure set_emission_processor_from_params
     module procedure set_emission_processor_from_input
  end interface

  interface defined
     module procedure fu_defined_processor
  end interface

  interface fu_type
     module procedure fu_type_processor
  end interface

  interface fu_grid
     module procedure fu_grid_proc
  end interface

  interface fu_vertical
     module procedure fu_vertical_proc
  end interface

  interface set_missing
     module procedure set_missing_processor
  end interface
  
  interface report
    module procedure report_emission_processor
  end interface

  type Temission_processor
     private
     real, dimension(:), allocatable :: scale_factor_fwd ! (nx*ny*nsp) or (nz*num_times*nsp)
     real, dimension(:), allocatable :: scale_factor_grad ! (nx*ny*nsp) or (nz*num_times*nsp)
     type(silja_grid) :: grid = grid_missing
     type(silam_vertical) :: vertical
     type(silam_species), dimension(:), allocatable :: species_list ! currently species_emission
     type(silja_time), dimension(:,:), allocatable :: time_slots ! for time/height mode 
     integer :: num_times 
     type(chemical_adaptor) :: adapt_transp_to_emis 
     logical :: defined = .false.
     integer :: emProc_type, nspecies
  end type Temission_processor
  public Temission_processor

  integer, parameter, public :: processor_scaling_flag = 91, &
                              & processor_void_flag = 90, &
                              & processor_tm_hgt_scale_flag = 92, &
                              & processor_tm_hgt_force_flag = 93, &
                              & processor_tm_hgt_force_wgt_flag = 94, &
                              & processor_time_height_flag = 95, &
                              & processor_log_scaling_flag = 96
                             
  ! Indices in the 1th dimension of time_slots array:
  integer, parameter, public :: slot_end = 2, slot_start = 1
  real, parameter, public :: tm_hgt_force_scaling = 1e3

contains

  !************************************************************

  subroutine set_da_source_indices(num_sources)
    implicit none
    integer, intent(in) :: num_sources
    
    select case (num_sources)
    case (1)
      DA_POSITIVE = 1
      DA_NEGATIVE = 1
      DA_ZERO = 1
      DA_ZERO_COEF_GRAD = 0.0
      DA_NEGT_COEF_OBS = 1.0
      DA_NEGT_COEF_GRAD = 0.0
    case (2)
      DA_POSITIVE = 1
      DA_NEGATIVE = 2
      DA_ZERO = 1
      DA_ZERO_COEF_GRAD = 0.0
      DA_NEGT_COEF_OBS = -1.0
      DA_NEGT_COEF_GRAD = -1.0
    case default
      call set_error('DA source index mode not available', 'set_da_source_indices')
    end select

    DA_NUM_SOURCES = num_sources
  end subroutine set_da_source_indices

  !************************************************************

  subroutine set_missing_processor(proc)
    implicit none
    type(Temission_processor), intent(out) :: proc
    proc%defined = .false.
    !
    ! Default initialization will occur.
    !
  end subroutine set_missing_processor

  !************************************************************

  logical function fu_defined_processor(proc) result(is_defined)
    implicit none
    type(Temission_processor), intent(in) :: proc
    is_defined = proc%defined
  end function fu_defined_processor

  !************************************************************

  function fu_data(proc, flag) result(data_ptr)
    implicit none
    type(Temission_processor), intent(in), target :: proc
    real, dimension(:), pointer :: data_ptr
    integer, intent(in) :: flag
    character(len=*), parameter :: sub_name = 'fu_data'
    
    if (proc%emProc_type == processor_void_flag) then
      nullify(data_ptr)
    elseif (flag == forwards) then
      data_ptr => proc%scale_factor_fwd
    elseif (flag == backwards) then
      data_ptr => proc%scale_factor_grad
    else
      nullify(data_ptr)
      call set_error("Wrong flag for emission_processor data", sub_name)
    end if

  end function fu_data

  !************************************************************

  integer function fu_type_processor(proc) result(proc_type)
    implicit none
    type(Temission_processor), intent(in) :: proc
    proc_type = proc%emProc_type
  end function fu_type_processor

  !************************************************************

  function fu_grid_proc(proc) result(grd)
    implicit none
    type(Temission_processor), intent(in) :: proc
    type(silja_grid) :: grd
    grd = proc%grid
  end function fu_grid_proc

  !************************************************************

  function fu_vertical_proc(proc) result(vert)
    implicit none
    type(Temission_processor), intent(in) :: proc
    type(silam_vertical) :: vert
    vert = proc%vertical
  end function fu_vertical_proc

  !************************************************************

  function fu_species_list(proc) result(species_list)
    implicit none
    type(Temission_processor), intent(in), target :: proc
    type(silam_species), dimension(:), pointer :: species_list
    species_list => proc%species_list
  end function fu_species_list


  !************************************************************************************

  subroutine set_emission_processor_from_params(map_emis, species_transport, proc_type, proc, &
                                              & time_slots)
    implicit none
    type(tMass_map), intent(in) :: map_emis
    type(silam_species), dimension(:), intent(in) :: species_transport
    type(Temission_processor), intent(out) :: proc
    ! time-resolving types only: dimension (2,:) for start and end times
    type(silja_time), dimension(:,:), intent(in), optional :: time_slots 

    integer, intent(in) :: proc_type

    integer :: nx, ny, stat, nsp, nsp_transp, num_times, nz
    type(silam_species), dimension(:), pointer :: species_emission

    select case (proc_type)
    case(processor_scaling_flag, processor_log_scaling_flag)
      if (fu_fails(.not. present(time_slots), 'Time dimension not supported', 'set_emission_processor')) return
      nx = map_emis%nx
      ny = map_emis%ny
      nsp = map_emis%nSpecies
      nsp_transp = size(species_transport)
      species_emission => map_emis%species
      allocate(proc%scale_factor_fwd(nx*ny*nsp), proc%scale_factor_grad(nx*ny*nsp), proc%species_list(nsp), stat=stat)
      if (fu_fails(stat == 0, 'Allocate failed', 'set_Temission_processor')) return


      if (proc_type == processor_scaling_flag) then
        proc%scale_factor_fwd(:) = 1. !! unity scaling by deafult. If defined is set below, the scale_factor must be meaningful.
      else
        proc%scale_factor_fwd(:) = 0. !! exp(0) = 1
      endif
      proc%scale_factor_grad(:) = 0.

      proc%grid = map_emis%gridTemplate
      call set_missing(proc%vertical, ifNew=.true.)
      proc%species_list = species_emission
      proc%nspecies = nsp
      proc%num_times = int_missing

      call create_adaptor(species_transport, species_emission, proc%adapt_transp_to_emis, &
                        & allow_missing=.true.)
      if (error) return

    case(processor_tm_hgt_scale_flag, processor_tm_hgt_force_flag, processor_tm_hgt_force_wgt_flag)
      if (fu_fails(present(time_slots), 'time_slots required', 'set_emission_processor')) return
      if (fu_fails(size(time_slots, 1) == 2, 'Strange shape for time_slots', 'set_emission_processor')) return
      nsp = map_emis%nspecies
      nsp_transp = size(species_transport)
      species_emission => map_emis%species
      num_times = size(time_slots, 2)
      proc%num_times = num_times
      nz = map_emis%n3d
      allocate(proc%scale_factor_fwd(nz*num_times*nsp),proc%scale_factor_grad(nz*num_times*nsp), &
            & proc%species_list(nsp), proc%time_slots(2, num_times),  stat=stat)
      if (fu_fails(stat == 0, 'Allocate failed', 'set_Temission_processor')) return
      proc%scale_factor_fwd(:) = 1.
      proc%scale_factor_grad(:) = 0.
      proc%grid = grid_missing
      proc%vertical = map_emis%vertTemplate
      proc%species_list = species_emission
      proc%nspecies = nsp
      proc%time_slots = time_slots
      call create_adaptor(species_transport, species_emission, proc%adapt_transp_to_emis, &
                        & allow_missing=.true.)
      if (error) return

    case(processor_void_flag)
      proc%grid = grid_missing

    case default
      call set_error('Strange proc_type', 'set_emission_processor')
      return
    end select

    proc%defined = .true.
    proc%emProc_type = proc_type

  end subroutine set_emission_processor_from_params

  
  !*******************************************************************
  
  subroutine set_emission_processor_from_input(disp_market, valid_time, map_emis, &
                                             & species_transport, proc_type, proc)
    implicit none
    type(mini_market_of_stacks), intent(in) :: disp_market
    type(tMass_map), intent(in) :: map_emis
    type(silja_time), intent(in) :: valid_time
    type(silam_species), dimension(:), intent(in) :: species_transport
    type(Temission_processor), pointer :: proc
    integer, intent(in) :: proc_type
    
    integer :: nx, ny, stat, nsp, nsp_transp, ind_start, ind_end, stride, isp_emis
    type(silam_species), dimension(:), pointer :: species_emission
    real, dimension(:), pointer :: emis_data_for_species
    type(silja_field_id) :: fid
    type(silja_field), pointer :: field_emis
    logical :: have_data

    have_data = .false.

    select case (proc_type)
      
    case(processor_scaling_flag, processor_log_scaling_flag)
      species_emission => map_emis%species
      do isp_emis = 1, size(species_emission)
        fid = fu_set_field_id(met_src_missing, ln_emission_scaling_flag, valid_time, zero_interval, &
                            & dispersion_grid, surface_level, species=species_emission(isp_emis))
        if (error) return
        call get_field_from_mm_general(disp_market, fid, field_emis, ifMandatory=.false.)
        if (error) then
          call set_error('Failed to set emission_processor from input', 'set_emission_processor_from_input')
          return
        end if

        if (.not. associated(field_emis)) cycle

        if (.not. have_data) then
          nx = map_emis%nx
          ny = map_emis%ny
          nsp = map_emis%nSpecies
          nsp_transp = size(species_transport)
          
          allocate(proc%scale_factor_fwd(nx*ny*nsp), proc%scale_factor_grad(nx*ny*nsp), &
                   &  proc%species_list(nsp), stat=stat)
          if (fu_fails(stat == 0, 'Allocate failed', 'set_Temission_processor')) return
          proc%scale_factor_fwd(:) = 1.
          proc%scale_factor_grad(:) = 0.
          proc%grid = map_emis%gridTemplate
          proc%species_list = species_emission
          proc%nspecies = nsp
        end if
        have_data = .true.
        emis_data_for_species => fu_grid_data(field_emis)
        ind_start = (isp_emis-1)*fs_dispersion + 1
        ind_end = ind_start + fs_dispersion - 1
        stride = size(species_emission)
        proc%scale_factor_fwd(ind_start:ind_end:stride) = emis_data_for_species(1:fs_dispersion)
      end do

    case(processor_tm_hgt_scale_flag, processor_tm_hgt_force_flag, processor_tm_hgt_force_wgt_flag)
      call set_error('Cannot set this processor type from disp market:' + &
                   & fu_str(processor_scaling_flag), 'set_emission_processor')
      return

    case default
      call set_error('Cannot set this type of emission_processor', 'set_emission_processor_from_input')
      return
    end select
    
    if (.not. have_data) call set_missing(proc)

  end subroutine set_emission_processor_from_input


  !*******************************************************************

  subroutine process_emission_forward(proc, map_emis, map_emis_px, map_emis_py, map_emis_pz, now, timestep)
    ! 
    ! Apply the emission adjustment in forward mode. Modifies the emission mass map
    ! according to the mode and content of the processor.
    ! Moments are assumed given, not centre of mass.
    ! If using centres of mass remove modifications to map_emis_px etc.
    ! 
    implicit none
    type(Temission_processor), intent(in) :: proc
    type(tmass_map), intent(inout) :: map_emis, map_emis_px, map_emis_py, map_emis_pz
    type(silja_time), intent(in) :: now
    type(silja_interval), intent(in) :: timestep

    integer :: ix, iy, isp, nx, ny, nspecies, ind, ilev, nz, ind_time
    integer, parameter :: isrc=1
    real :: scaling, scaling_max, seconds, mass
    real, dimension(map_emis%n3d) :: lev_wgt
    
    if (.not. proc%defined) return
    if (map_emis%nspecies /= map_emis_px%nspecies) then
      call set_error('Bulk advection not working with emission processors', 'process_emission_forward')
      return
    end if

    seconds = fu_sec(timestep)
    
    select case(proc%emProc_type)
      case(processor_void_flag)
        return

      case(processor_scaling_flag, processor_log_scaling_flag)
        call grid_dimensions(proc%grid, nx, ny)
        nspecies = proc%nspecies
        scaling_max = -1
        map_emis%ifGridValid = .true.

        !$OMP PARALLEL default(none) shared(nx, ny, nspecies, proc, &
        !$OMP   & map_emis, map_emis_px, map_emis_py, map_emis_pz) &
        !$OMP  private(ix,iy,isp,ind,scaling)
        !$OMP DO collapse(2)
        do iy = 1, ny
          do ix = 1, nx
            do isp = 1, nspecies
              ind = (iy-1)*nspecies*nx + (ix-1)*nspecies + isp
              if (proc%emProc_type == processor_scaling_flag) then
                scaling = proc%scale_factor_fwd(ind)
              else
                scaling = exp(proc%scale_factor_fwd(ind))
              endif
              map_emis%arm(isp,isrc,:,ix,iy) = map_emis%arm(isp,isrc,:,ix,iy)*scaling 
              map_emis_px%arm(isp,isrc,:,ix,iy) = map_emis_px%arm(isp,isrc,:,ix,iy)*scaling
              map_emis_py%arm(isp,isrc,:,ix,iy) = map_emis_py%arm(isp,isrc,:,ix,iy)*scaling
              map_emis_pz%arm(isp,isrc,:,ix,iy) = map_emis_pz%arm(isp,isrc,:,ix,iy)*scaling
            end do
          end do
        end do
        !$OMP END DO
        !$OMP END PARALLEL
      case(processor_tm_hgt_scale_flag, processor_tm_hgt_force_flag, processor_tm_hgt_force_wgt_flag)
        do ind_time = 1,  proc%num_times
          if (proc%time_slots(slot_start, ind_time) <= now .and. &
            & now < proc%time_slots(slot_end, ind_time)) exit
        end do
        if (ind_time > proc%num_times) then
          call msg_warning('No time slot found for processor')
          return
        else
          call msg('time_height_emission forward: slot start: ' &
                 & // fu_str(proc%time_slots(slot_start, ind_time)))
        end if

        nz = map_emis%n3d
        nspecies = proc%nspecies
        if (proc%emProc_type == processor_tm_hgt_scale_flag) then
          do ilev = 1, nz
            do isp = 1, nspecies
              ind = (ind_time-1)*nspecies*nz + (ilev-1)*nspecies + isp
              scaling = proc%scale_factor_fwd(ind)
              map_emis%arm(isp,isrc,ilev,:,:) = map_emis%arm(isp,isrc,ilev,:,:) * scaling
              map_emis_px%arm(isp,isrc,ilev,:,:) = map_emis_px%arm(isp,isrc,ilev,:,:)*scaling
              map_emis_py%arm(isp,isrc,ilev,:,:) = map_emis_py%arm(isp,isrc,ilev,:,:)*scaling
              map_emis_pz%arm(isp,isrc,ilev,:,:) = map_emis_pz%arm(isp,isrc,ilev,:,:)*scaling
              if (scaling > scaling_max) scaling_max = scaling
            end do
          end do
          call msg('Max scaling', scaling_max)

        else
          ! The forcing mode. The sqrt-weighting can be applied, see top of the module. 
          ! The forced emission is applied with zero centre of mass.
          if (proc%emProc_type == processor_tm_hgt_force_wgt_flag) then
            lev_wgt(1:nz) = (/(sqrt(fu_layer_thickness_m(fu_level(map_emis%vertTemplate, ilev))), ilev=1, nz)/)
          else
            lev_wgt(1:nz) = 1.0
          end if
          do ilev = 1, nz
            do isp = 1, nspecies
              ind = (ind_time-1)*nspecies*nz + (ilev-1)*nspecies + isp
              mass = proc%scale_factor_fwd(ind)*seconds
              where (map_emis%arm(isp,isrc,ilev,:,:) > 0)
                map_emis%arm(isp,isrc,ilev,:,:) = mass * tm_hgt_force_scaling  * lev_wgt(ilev)
              end where
              map_emis_px%arm(isp,isrc,ilev,:,:) = 0.0!map_emis_px%arm(isp,isrc,ilev,:,:)*scaling
              map_emis_py%arm(isp,isrc,ilev,:,:) = 0.0!map_emis_py%arm(isp,isrc,ilev,:,:)*scaling
              map_emis_pz%arm(isp,isrc,ilev,:,:) = 0.0!map_emis_pz%arm(isp,isrc,ilev,:,:)*scaling
            end do
          end do

          !call msg('sum of emission map after forcing:', sum(map_emis%arm(:,:,:,:,:)))
        end if
      
      case default
        call set_error('Strange processor type', 'process_emission_forward')
        return
    end select
    
  contains
    
    subroutine report_map(map)
      implicit none
      type(tmass_map), intent(in) :: map
      
      do isp = 1, nspecies
        call msg('Species, total emission', isp, sum(map%arm(isp,isrc,:,:,:)))
      end do
      
    end subroutine report_map
    
  end subroutine process_emission_forward

  
  !*******************************************************************

  subroutine process_emission_adjoint(map_emis, map_transp, cell_metric, buf, proc, now, timestep)
    !!   "observation on" adjoint massmap
    ! 
    ! Apply the emission processor in adjoint mode. In the end of assimilation/inversion
    ! window, the action of the adjoint processor recovered from the processor data array.
    ! In practice, this means collecting the emissin masses, with appropriate scalings and
    ! integrating in the dimensions not resolved by the processor.
    ! 
    implicit none
    type(tmass_map), intent(in) :: map_emis
    type(tmass_map), intent(in) :: map_transp
    integer, intent(in) :: cell_metric
    TYPE(Tfield_buffer), intent(in) :: buf  !! dispersion_buffer
    type(Temission_processor), intent(inout) :: proc
    type(silja_time), intent(in) :: now
    type(silja_interval), intent(in) :: timestep

    integer :: ix, iy, isp_transp, isp_emis, ind, ilev, nx, ny, nlevs, nspecies, nspecies_transp, ind_time
    real :: transp_val, dEmis_dScale, volume, seconds, weight_past, cell_area
    type(chemical_adaptor) :: adaptor
    real, dimension(map_transp%n3d) :: thickness, lev_wgt
    type(field_3d_data_ptr), pointer  :: airmass3dp, airmass3df
    integer :: iSpOnes
    character(len=*), parameter :: sub_name = 'process_emission_adjoint'

    if (.not. proc%defined) return
    if (proc%emProc_type == processor_void_flag) return

    if (cell_metric /= cloud_metric_geometry_flag) then
      if ( .not. any( proc%emProc_type == (/processor_scaling_flag,processor_log_scaling_flag/))) then
        call set_error("Unimplemented combination of cell_metric and emProc_type", sub_name)
        return
      ! RK: Should be trivial to implement here, but I have no test case to debug it
      ! RK: processor_scaling_flag can be used as an example
      endif
    endif

    adaptor = proc%adapt_transp_to_emis
    nx = map_emis%nx
    ny = map_emis%ny
    nlevs = map_emis%n3d
    nspecies_transp = map_transp%nspecies
    nspecies = proc%nspecies

    ind = int_missing
    airmass3dp => null()
    airmass3df => null()
    select case (cell_metric)
      case (cloud_metric_geometry_flag)
         thickness(1:nlevs) = (/(fu_layer_thickness_m(fu_level(map_emis%vertTemplate, ilev)), ilev=1, nlevs)/)
       case (cloud_metric_cellmass_flag)
         ind = fu_index(buf, disp_cell_airmass_flag)
         if (ind < 1) then                              
           call set_error('Failed to find the air mass field', sub_name)
           return
         end if
         airmass3dp => buf%p4d(ind)%past
         airmass3df => buf%p4d(ind)%future
         weight_past = (buf%time_future - now) / (buf%time_future - buf%time_past)
       case(cloud_metric_ones_flag)
         iSpOnes = select_single_species(map_transp%species, nspecies_transp, &
                           & 'ones', in_gas_phase, real_missing)
         if (error .or. (iSpOnes < 1)) then
           call set_error("Couldn't find 'ones' species", sub_name)
           return
         endif
         !! Stupidity check
         volume = map_transp%arm(iSpOnes, 1, 1, 1, 1)    ! kg or whatever
         if( fu_fails(volume > 1, "Strange ones mass", sub_name)) return
      case default
        call set_error("Strange cloud_metric_flag: "//trim(fu_str(cell_metric)), sub_name)
         return
    end select 
    
    select case(proc%emProc_type) 
    case (processor_scaling_flag, processor_log_scaling_flag)
      
      call msg("Processing adjoint emissions")
      
      if (fu_fails(map_transp%nsrc == DA_NUM_SOURCES, &
           & 'Wrong number of sources in transport map', sub_name)) return

      if (fu_fails(map_transp%gridTemplate == proc%grid, 'Grids do not match', sub_name)) return
      
      !$OMP parallel if (.TRUE.) default(none) &
      !$OMP & shared(map_emis,map_transp,nx,ny,nlevs, nspecies_transp, nspecies,cell_metric,thickness, &
      !$OMP & airmass3dp,airmass3df,proc,weight_past, ispones, adaptor, DA_POSITIVE, DA_ZERO, DA_NEGATIVE,&
      !$OMP & DA_ZERO_COEF_GRAD, DA_NEGT_COEF_GRAD) &
      !$OMP & private(cell_area,ix,iy,ilev,isp_transp,volume,ind,isp_emis,dEmis_dScale,transp_val)
      !$OMP DO  collapse (2)
      do iy = 1, ny
        do ix = 1, nx
          if(cell_metric == cloud_metric_geometry_flag) cell_area = fu_cell_size(map_emis%gridTemplate, ix, iy)
          do ilev = 1, nlevs
            if (cell_metric == cloud_metric_geometry_flag) then
              volume = thickness(ilev) * cell_area
            elseif (cell_metric == cloud_metric_cellmass_flag) then
              ind = ix + nx * (iy - 1) !index in 2d field
              volume = airmass3dp%p2d(ilev)%ptr(ind) * weight_past + &  !! Actually cell mass
                       airmass3df%p2d(ilev)%ptr(ind) * (1. - weight_past)
            elseif (cell_metric == cloud_metric_ones_flag) then
              volume = map_transp%arm(iSpOnes, 1, ilev, ix, iy)    ! kg or whatever
            else 
              call set_error("Strange cloud_metric_flag: "//trim(fu_str(cell_metric)), sub_name)
            endif
            do isp_transp = 1, nspecies_transp
              isp_emis = adaptor%isp(isp_transp)
              if (isp_emis == int_missing) cycle
              ind = (iy-1)*nspecies*nx + (ix-1)*nspecies + isp_emis

              dEmis_dScale = map_emis%arm(isp_emis,1,ilev,ix,iy) !! just emission for now
              !! dEmis_dScale = dEmis / dScaling
              if (proc%emProc_type == processor_log_scaling_flag) then
                  dEmis_dScale = dEmis_dScale * exp(proc%scale_factor_fwd(ind))
              endif
              transp_val = map_transp%arm(isp_transp,DA_POSITIVE,ilev,ix,iy) &
                   & + DA_ZERO_COEF_GRAD * map_transp%arm(isp_transp,DA_ZERO,ilev,ix,iy) &
                   & + DA_NEGT_COEF_GRAD * map_transp%arm(isp_transp,DA_NEGATIVE,ilev,ix,iy)

              proc%scale_factor_grad(ind) = proc%scale_factor_grad(ind) + dEmis_dScale*transp_val / volume
            end do
          end do
        end do
      end do
      !$OMP ENDDO
      !$OMP END parallel

    case (processor_tm_hgt_scale_flag, processor_tm_hgt_force_flag, processor_tm_hgt_force_wgt_flag)

      do ind_time = 1,  proc%num_times
        if (proc%time_slots(slot_start, ind_time) <= now .and. now < proc%time_slots(slot_end, ind_time)) exit
      end do
      if (ind_time > proc%num_times) then
        call msg_warning('No time slot found for processor')
        return
      else
        call msg('time_height_emission adjoint: slot start: ' &
               & // fu_str(proc%time_slots(slot_start, ind_time)))
      end if
      ! summation over x, y (and timeslot), ilev, isp resolved

      if (proc%emProc_type == processor_tm_hgt_scale_flag) then
        do iy = 1, ny
          do ix = 1, nx
            do ilev = 1, nlevs
              volume = thickness(ilev) * fu_cell_size(map_emis%gridTemplate, ix, iy)
              do isp_transp = 1, nspecies_transp
                isp_emis = adaptor%isp(isp_transp)
                if (isp_emis == int_missing) cycle
                ind = (ind_time-1)*nspecies*nlevs + (ilev-1)*nspecies + isp_emis
                dEmis_dScale = map_emis%arm(isp_emis,1,ilev,ix,iy)
                transp_val = map_transp%arm(isp_transp,DA_POSITIVE,ilev,ix,iy) &
                     & + DA_ZERO_COEF_GRAD * map_transp%arm(isp_transp,DA_ZERO,ilev,ix,iy) &
                     & + DA_NEGT_COEF_GRAD * map_transp%arm(isp_transp,DA_NEGATIVE,ilev,ix,iy)
                proc%scale_factor_grad(ind) = proc%scale_factor_grad(ind) + dEmis_dScale*transp_val / volume
              end do
            end do
          end do
        end do

      else
        if (proc%emProc_type == processor_tm_hgt_force_wgt_flag) then
          lev_wgt(1:nlevs) = sqrt(thickness(1:nlevs))
        else
          lev_wgt(1:nlevs) = 1.0
        end if
        seconds = abs(fu_sec(timestep))
        do iy = 1, ny
          do ix = 1, nx
            do ilev = 1, nlevs
              do isp_transp = 1, nspecies_transp
                isp_emis = adaptor%isp(isp_transp)
                if (isp_emis == int_missing) cycle
                ind = (ind_time-1)*nspecies*nlevs + (ilev-1)*nspecies + isp_emis
                dEmis_dScale = map_emis%arm(isp_emis,1,ilev,ix,iy)
                if ((dEmis_dScale .eps. 0.0)) cycle
                volume = thickness(ilev) * fu_cell_size(map_emis%gridTemplate, ix, iy)
                transp_val = map_transp%arm(isp_transp,DA_POSITIVE,ilev,ix,iy) &
                     & + DA_ZERO_COEF_GRAD * map_transp%arm(isp_transp,DA_ZERO,ilev,ix,iy) &
                     & + DA_NEGT_COEF_GRAD * map_transp%arm(isp_transp,DA_NEGATIVE,ilev,ix,iy)
                !proc%scale_factor(ind) = proc%scale_factor(ind) + transp_val/volume! * seconds
                proc%scale_factor_grad(ind) = proc%scale_factor_grad(ind) &
                     & + transp_val/volume * seconds * tm_hgt_force_scaling * lev_wgt(ilev)
                
              end do
            end do
          end do
        end do

      end if
    case default
      call set_error('Strange processor type', sub_name)
      return
      
    end select
  end subroutine process_emission_adjoint


  !************************************************************

  subroutine zero_processor_grad(proc)
    implicit none
    type(Temission_processor), intent(inout) :: proc
    if (allocated(proc%scale_factor_grad)) proc%scale_factor_grad(:) = 0.0
  end subroutine zero_processor_grad

  
  !************************************************************

  subroutine report_emission_processor(proc, message_)
    !
    ! Emission processor report
    !
    implicit none
    
    ! Imported parameters
    type(Temission_processor), intent(in) :: proc
    character(len=*), intent(in), optional :: message_
    character(len=fnlen) :: message

    
    ! Local variables
    integer :: iTmp
    
    if(present(message_)) then 
        message = message_
    else
       message = ""
    endif

    if(proc%defined)then
      if (proc%emProc_type == processor_void_flag) then
        call msg(trim(message) // ": VOID emission processor")
      else
        call msg('')
        call msg(message)
        call msg('----------------- emission processor report ----------------')
        call msg('Processor type:' + fu_str(proc%emProc_type))
        call report(proc%grid)
        call report(proc%vertical)
        call msg('Number of species and times:', proc%nspecies, proc%num_times)
        do iTmp = 1, proc%nspecies
          call msg('Species:' + fu_str(proc%species_list(iTmp)))
        end do
        if(proc%num_times > 0 .and. allocated(proc%time_slots))then
          call msg('Starting time:' + fu_str(proc%time_slots(1,1)))
          call msg('Ending time:' + fu_str(proc%time_slots(2,proc%num_times)))
        endif
        
        if(allocated(proc%scale_factor_fwd))then
          call msg('scale_factor_fwd min / max:', minval(proc%scale_factor_fwd), maxval(proc%scale_factor_fwd))
          call msg('scale_factor_fwd minloc:', minloc(proc%scale_factor_fwd, KIND=4))
          call msg('scale_factor_fwd maxloc:', maxloc(proc%scale_factor_fwd, KIND=4))
        else
          call msg('Scale_factor_fwd is not allocated')
        endif

        if(allocated(proc%scale_factor_grad))then
          call msg('scale_factor_grad norm:', sqrt(dot_product(proc%scale_factor_grad, proc%scale_factor_grad)))
        else
          call msg('Scale_factor_grad is not allocated')
        endif

  !      call report(proc%adapt_transp_to_emis) no such function
        call msg('----------------- END emission processor report ----------------')
      endif
    else
      call msg(trim(message) // ": UNDEFINED emission processor")
    endif
    
  end subroutine report_emission_processor

end module source_apportionment

