program test_perturb
  use perturbations
  use silam_namelist
  use da_interface
  implicit none
  
  integer, parameter :: nx = 200, ny = 200, nz = 5, fs = nx*ny
  
  real, dimension(:), allocatable :: xx
  real :: xx_mean, xx_var

  allocate(xx(fs))

  call random_normal(xx)

  xx_mean = sum(xx) / size(xx)
  print *, 'Mean:', xx_mean
  xx_var = sum((xx - xx_mean)**2) / size(xx)
  print *, 'Var:', xx_var

  ! test the last element
  deallocate(xx)
  allocate(xx(11))
  call random_normal(xx)
  print *, xx(11)

  !call test_with_rules()
  call test_volc()
  
  print *, 'Bye'

contains

  subroutine test_volc()
    integer :: file_unit 
    type(Tsilam_namelist), pointer :: nl
    type(da_rules) :: darules
    type(silja_time) :: start_time
    type(silja_interval) :: period_to_compute
    
    type(silja_grid) :: grd
    type(silam_vertical) :: vert
    type(silam_species) :: species_o3
    type(t_perturbation), dimension(8) :: perturbations
    real, dimension(:), allocatable :: pert_values
    type(Tmass_map), pointer :: mass_map, pmm
    type(silam_pollution_cloud) :: cloud
    type(Tfield_buffer) :: buf
    type(silja_time), dimension(:,:), pointer :: time_slots
    type(Temission_processor) :: proc
    type(silja_interval) :: timestep 
    integer :: ii, jj
    integer, parameter :: num_samples = 1000, num_times = 72*4

    file_unit = fu_next_free_unit()
    call init_random_seed(44)

    open(file='test_perturb.ini', unit=file_unit)
    nl => fu_read_namelist(file_unit)
    if (.not. associated(nl)) then
      print *, 'Not associated nl'
      stop
    end if
    close(file_unit)

    open(run_log_funit, file='test_perturb.log')

    period_to_compute = fu_set_named_interval(fu_content(nl, 'computed_period'))
    start_time = fu_io_string_to_time(fu_content(nl, 'start_time'))

    call set_da_rules(nl, darules, start_time, period_to_compute)

    grd = fu_set_grid(nl)
    if (error) return
    call set_vertical(nl, vert)
    if (error) return
    
    call init_chemical_materials(nl)
    call set_species(species_o3, fu_get_material_ptr('O3'), in_gas_phase)

    perturbations(:)%target = int_missing
    perturbations(1)%target = perturb_concentration
    call set_perturbations(perturbations, (/species_o3/), (/species_o3/), darules, grd, vert, start_time)
    
    mass_map => fu_set_mass_map(concentration_flag, 1, 0, grd, vert, (/species_o3/), val=0.0)
    
    call get_emis_time_slots(darules, time_slots)
    call set_emission_processor(mass_map, (/species_o3/), processor_tm_hgt_force, proc, time_slots)

    timestep = fu_set_interval_sec(900.0)
    
    call msg('n. timeslots', size(time_slots,2))
    
    call msg('proc data:', fu_data(proc))
    
    write(71, fmt=*) (/(fu_level_height(fu_level(vert, ii, .true.)), ii=1, fu_nbrOfLevels(vert))/)
    
    do ii = 1, num_samples
      print *, 'sample', ii
      perturbations(1)%first_time = .true.
      do jj = 1, num_times
        call perturb_emis_proc_volc(proc, perturbations(1), timestep)
        write(72, fmt=*) fu_data(proc)
      end do
    end do

  end subroutine test_volc

  subroutine test_with_rules()
    implicit none
    integer :: file_unit 
    type(Tsilam_namelist), pointer :: nl
    type(da_rules) :: darules
    type(silja_time) :: start_time
    type(silja_interval) :: period_to_compute
    
    type(silja_grid) :: grd
    type(silam_vertical) :: vert
    type(silam_species) :: species_o3
    type(t_perturbation), dimension(8) :: perturbations
    real, dimension(:), allocatable :: pert_values
    type(Tmass_map), pointer :: mass_map, pmm
    type(silam_pollution_cloud) :: cloud
    type(Tfield_buffer) :: buf

    file_unit = fu_next_free_unit()

    open(file='test_perturb.ini', unit=file_unit)
    nl => fu_read_namelist(file_unit)
    if (.not. associated(nl)) then
      print *, 'Not associated nl'
      stop
    end if
    close(file_unit)

    period_to_compute = fu_set_named_interval(fu_content(nl, 'computed_period'))
    start_time = fu_io_string_to_time(fu_content(nl, 'start_time'))

    call set_da_rules(nl, darules, start_time, period_to_compute)

    grd = fu_set_grid(nl)
    if (error) return
    call set_vertical(nl, vert)
    if (error) return
    
    call init_chemical_materials(nl)
    call set_species(species_o3, fu_get_material_ptr('O3'), in_gas_phase)
    
    perturbations(:)%target = int_missing
    perturbations(1)%target = perturb_concentration
    call set_perturbations(perturbations, (/species_o3/), (/species_o3/), darules, grd, vert, start_time)

    allocate(pert_values(perturbations(1)%pert_size))
    call get_perturbation(perturbations(1), pert_values)

    xx_mean = sum(pert_values) / size(pert_values)
    print *, 'Mean:', xx_mean
    xx_var = sum((pert_values - xx_mean)**2) / size(pert_values)
    print *, 'Var:', xx_var
    
    print *, 'Size of pert:', size(pert_values)
    open(file_unit, access='stream', file='pert.out')
    write(file_unit) pert_values
    close(file_unit)
    
    mass_map => fu_set_mass_map(concentration_flag, 1, 0, grd, vert, (/species_o3/), val=0.0)

    call perturb_mass_map(mass_map, perturbations(1))
    
    open(file_unit, access='stream', file='pert_lev2.out')
    write(file_unit) mass_map%arm(1,1,2,:,:)
    close(file_unit)
    

    call report(grd)
    call report(vert)
  end subroutine test_with_rules

  
end program test_perturb
