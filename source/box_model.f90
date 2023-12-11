program box
  use cbm4_Precision
  use chemistry_manager
  !use io_server, only : fu_species_output_name
  use chemical_setup
  use cocktail_basic
  use photolysis 

  implicit none
  
  type Tbox_setup
     type(silja_time) :: start_time, end_time
     type(silja_interval) :: timestep
     character (len=fnlen) :: case_name
     real :: lat = 45.0
     real :: lon = 0.0
     integer :: output_unit_cnc, output_unit_nbr
  end type Tbox_setup

  type(silam_species), dimension(:), pointer :: species_normal, species_shortlived, &
       species_aerosol, species_emission
  integer :: nspecies, nspecies_shortlived = 0, nspecies_aerosol = 0, nspecies_emission = 0
  integer :: nreact_rates = 0
  integer :: retcode = 0
  
  type(Tmeteo_input) :: meteo_input
  type(Tchem_rules), pointer :: pChemrules
  type(Tbox_setup) :: box_setup

  real(sp), dimension(:), allocatable, save :: cnc_shortlived, cnc_normal, cnc_aerosol
  real, dimension(:,:), allocatable, save :: metdat
  real(sp), dimension(:,:), allocatable, save :: tan_lin_traj
  real(sp), dimension(:), allocatable :: h_start_arr

  nullify(species_normal, species_shortlived, species_aerosol, species_emission)
  allocate(pChemRules)

  call setup(pChemrules)
  if ( .not. error) then 

    call msg('')
    call adjoint_test(box_setup, pchemrules)
 !   stop
 !    call run(box_setup, pChemrules)
    close(box_setup%output_unit_cnc)
    close(box_setup%output_unit_nbr)
  endif


  if (error) retcode = 10 !! Some nonzero
  call msg("Returning", retcode)
  close(run_log_funit)
  call exit_with_status(retcode)



  !************************************************************************************

contains

  subroutine initial_conditions(nlInit, species_list, values, species_list_size)
 
    type(Tsilam_namelist), intent(in) :: nlInit
    type(silam_species), dimension(:), pointer :: species_list
    real(sp), dimension(:) :: values
    integer, intent(in) :: species_list_size

    integer :: iSpecies, iostat
    character(len=clen) :: ident, unit
    character(len=fnlen) ::  line, content, default
    real :: val, conversion
    type(silam_material), pointer :: material
    integer :: ind_pres, ind_tempr
    logical :: first_ppb

    default = fu_content(nlInit, 'default')

    first_ppb = .true.
    do iSpecies = 1, species_list_size 
      ident = fu_str(species_list(ispecies))
      val = fu_content_real(nlInit, ident)
      content = fu_content(nlInit, ident)
      if (content == '') then
        if (default == '') then
          call msg('Looking for:')
          call report(species_list(iSpecies))
          call msg('With name ' // trim(ident))
          call set_error('Initial condition not defined', 'initial_conditions')
          return
        else
          content = default
        end if
      end if
      read(content, fmt=*, iostat=iostat) val, unit
      if (fu_fails(iostat==0, 'Failed to parse:'+ident+'='+content, 'initial_conditions'))then
        return
      endif
      material => fu_material(species_list(ispecies))
      if (fu_fails(associated(material), 'Material not associated', 'initial_conditions')) return
      if (unit == 'ppb') then
        ind_pres = fu_index_in_meteo_input(meteo_input, pressure_flag)
        ind_tempr = fu_index_in_meteo_input(meteo_input, temperature_flag)
        if (fu_fails(ind_pres /= int_missing, 'Where is pressure', 'initial_conditions')) return
        if (fu_fails(ind_tempr /= int_missing, 'Where is temperature', 'initial_conditions')) return
        conversion = metdat(ind_pres,1) / (gas_constant_uni * metdat(ind_tempr,1)) * 1e-9
        if (first_ppb) call msg('PPB to concentration conversion:', conversion)
        first_ppb = .false.
      else
        conversion = fu_conversion_factor(unit, fu_basic_mass_unit(material), material)
      end if
      values(iSpecies) = val*conversion
      write(line, fmt='(I3,2x, A20, G20.3, A10, A3)') iSpecies, trim(ident), val*conversion, &
           & trim(fu_basic_mass_unit(material)), '/m3'
      call msg(line)
    end do

  end subroutine initial_conditions
  
  !************************************************************************************
 
  subroutine init_aerosol_distr(nlInit, species_mass, values_mass, mass_list_size, &
                                & species_nbr, values_nbr, nbr_list_size)
 
    type(Tsilam_namelist), pointer :: nlInit
    type(silam_species), dimension(:), pointer :: species_mass, species_nbr
    real(sp), dimension(:) :: values_mass, values_nbr
    integer, intent(in) :: mass_list_size, nbr_list_size


    integer :: iTmp, nModes, i, iType, iStatus, S
    real :: f1, f2, cnc
    character(len=clen) ::  chMatNm, chType
    type(silam_species), dimension(:), pointer :: speciesTmp  
    type(TspeciesReference), dimension(:), pointer :: mapMass,  mapNbr
    type(Tsilam_nl_item_ptr), dimension(:), pointer :: ptrItems
    type(silam_sp) :: spContent

    allocate(speciesTmp(1))
    nullify(ptrItems)
    call get_items(nlInit, 'aerosol_mode', ptrItems, nModes)
    allocate(spContent%sp)

    do iTmp = 1, nModes
      
    spContent%sp = fu_content(ptrItems(iTmp))
      read(unit=spContent%sp, iostat=iStatus, fmt=*) chMatNm, chType, f1, f2, s, cnc

      select case(trim(chType))
        case('lognormal')
          iType = lognormal_flag
        case('fixed_diameter')
          iType = fixed_diameter_flag
        case('moving_diameter')
          iType = moving_diameter_flag
        case default
          call set_error('strange distr type','init_aerosol_distr')
          return
      end select


!      f1 = f1 / exp(log(f2)**2. / 2.)  ! mean d to median d

!      f1 = exp(log(f1) - 3.0 * log(f2)**2)                             ! from vol med to nbr med

      
      call set_species(speciesTmp(1), fu_get_material_ptr(chMatNm), &
                     & fu_set_mode(iType, f1, f2, solubility = s))

      call create_mode_projection(speciesTmp, 1, &   
                              & species_mass, mass_list_size, &   
                              & species_nbr, nbr_list_size, &     
                              & mapMass, mapNbr, .true.) 

      !  Mass
      do i = 1, mapMass(1)%nRefSpecies
        values_mass(mapMass(1)%indSpeciesTo(i)) = values_mass(mapMass(1)%indSpeciesTo(i)) + &
                                              & cnc * mapMass(1)%fract(i)
      enddo

      ! Number
      do i = 1, mapNbr(1)%nRefSpecies
        values_nbr(mapNbr(1)%indSpeciesTo(i)) = values_nbr(mapNbr(1)%indSpeciesTo(i)) + &
                                              & cnc * mapNbr(1)%fract(i)
      enddo

    enddo

  end subroutine init_aerosol_distr
  
  !************************************************************************************

  subroutine get_meteo_params(metInput, nlMet, values)
    implicit none
    type(Tmeteo_input), intent(in) :: metInput
    real, dimension(:,:), intent(out) :: values
    type(Tsilam_namelist), pointer :: nlMet

    integer :: i
    real :: val
    character(len=clen) :: quantname

    if (.not. associated(nlmet)) then
      call set_error('meteo namelist missing', 'get_meteo_params')
      return
    end if

    do i = 1, metInput%nQuantities
      quantname = fu_quantity_short_string(metInput%quantity(i))
      val = fu_content_real(nlMet, quantname)
      if (val .eps. real_missing) then
        call set_error('Missing quantity ' // trim(quantname) // ' in meteo namelist', &
                     & 'get_meteo_params')
      end if
      values(i,1) = val 
    end do
    

  end subroutine get_meteo_params

  !************************************************************************************

  subroutine setup(chemrules)
    implicit none

    type(Tchem_rules), intent(inout) :: chemrules

    character(len=fnlen) :: inifilename, filename, emis_cockt_name, inidir
    integer :: stat, inifile, stdsetupfile, iTmp
    type(Tsilam_namelist_group), pointer :: nlSetupGrp
    type(Tsilam_namelist), pointer :: nl, nlStdSetup
    integer, dimension(:), pointer :: q_met_dyn, q_met_stat, q_disp_dyn, q_disp_stat
    type(Tcocktail_descr) :: descr_emis_cockt
    logical :: ifSpecies
    type(Tsilam_nl_item_ptr), dimension(:), pointer :: pItems
    integer :: nItems
    
    run_log_funit = fu_next_free_unit()
    open(run_log_funit, file='box.log', iostat=stat, action='write')
    if (stat /= 0) then
      call set_error('Failed to open log file', 'setup')
    end if

    !!! Own file name
    CALL GET_COMMAND_ARGUMENT(0, filename, STATUS=stat)
    if (stat/=0) then
      call set_error('Failed to get own filename', 'setup')
      return
    end if


    CALL GET_COMMAND_ARGUMENT(1, inifilename, STATUS=stat)
    if (stat/=0) then
      call set_error('Usage: '// trim(filename)//' box_ini_file.ini', 'setup')
      return
    end if

    inifile = fu_next_free_unit()
    open(inifile, file=inifilename, iostat=stat, action='read', status='old')
    if (stat /= 0) then
      call set_error('Failed to open ini file: ' // inifilename, 'setup')
      return
    end if

    nlSetupGrp => fu_read_namelist_group(inifile, ifEnv=.true.)
    if (error) return
    close(inifile)

    !nlTransf => fu_read_namelist(inifile)
    
    nl => fu_namelist(nlSetupGrp, 'setup')
    if (.not.  associated(nl)) then
      call set_error('Missing setup namelist', 'setup')
      return
    end if
    
    box_setup%start_time = fu_io_string_to_time(fu_content(nl, 'start_time'))
    box_setup%end_time = fu_io_string_to_time(fu_content(nl, 'end_time'))
    box_setup%timestep = fu_set_named_interval(fu_content(nl, 'timestep'), box_setup%start_time)
    box_setup%case_name = fu_content(nl, 'case_name')
    if (error) return
    
    box_setup%lat = fu_content_real(nl, 'latitude')
    if (box_setup%lat == real_missing) then
      call set_error('Missing latitude in setup namelist', 'setup')
      return
    end if
    box_setup%lon = fu_content_real(nl, 'longitude')
    if (box_setup%lon == real_missing) then
      call set_error('Missing longitude in setup namelist', 'setup')
      return
    end if

    nlStdSetup => fu_namelist(nlSetupGrp, 'STANDARD_SETUP')

    
    nl => fu_namelist(nlSetupGrp, 'transformation_parameters')
    if (.not. associated(nl)) then
      call set_error('Missing transformation_parameters namelist', 'setup')
      return
    end if

    inidir = fu_content(nlStdSetup,'standard_setup_directory')
    call expand_paths(nlStdSetup, inidir)
    call init_chemical_materials(fu_content(nlStdSetup,'chemical_database_fnm'), &
                               & fu_content(nlStdSetup,'nuclide_database_fnm'))
    
!    call set_standard_cocktails_fname(fu_expand_environment(fu_content(nlStdSetup,'standard_cocktail_fnm')))
    ! Set the standard cocktail descriptions
    !
    call get_items(nlStdSetup, 'standard_cocktail_fnm', pItems, nItems)
    if(error)return
    call read_standard_cocktail_descrs(pItems, nItems)
    if(error)return
    
    
    call set_chemistry_rules(nl, nlStdSetup, chemrules, .true., .False.)
    if (error) return
    
    ! Some transformations depend on emission species. We can fake them with the
    ! emission_cocktail input field.
    emis_cockt_name = fu_content(nl, 'emission_cocktail')
    if (emis_cockt_name /= '') then
      call msg('Forced emission cocktail: ' // trim(emis_cockt_name))
      call set_cocktail_description(emis_cockt_name, descr_emis_cockt, ifSpecies)
      if (error) return
      call get_inventory(descr_emis_cockt, species_emission, nspecies_emission)
      if (error) return
      ! Do not deallocate the descriptor, emission species now point into it.
    end if


    q_met_dyn => fu_work_int_array()
    q_disp_dyn => fu_work_int_array()
    q_met_stat => fu_work_int_array()
    q_disp_stat => fu_work_int_array()

    q_met_dyn(:) = int_missing
    q_disp_dyn(:) = int_missing
    q_met_stat(:) = int_missing
    q_disp_stat(:) = int_missing

    
    call add_transformation_input_needs(chemrules, q_met_dyn, q_met_stat, q_disp_dyn, &
                                      & q_disp_stat, meteo_input)
    if (error) return
    call free_work_array(q_disp_dyn)
    allocate(metdat(meteo_input%nquantities, 1))
    call get_meteo_params(meteo_input, fu_namelist(nlSetupGrp, 'meteo'), metdat)
    if (error) return
    
    
    call global_chemical_init(chemrules, box_setup%timestep, box_setup%timestep, &
                            & species_emission, nspecies_emission,&
                            & species_normal, nspecies, &
                            & species_shortlived, nspecies_shortlived, &
                            & species_aerosol, nspecies_aerosol, nreact_rates)


    call free_work_array(q_disp_stat)
    call free_work_array(q_met_dyn)
    call free_work_array(q_met_stat)
    if (error) return

    nl => fu_namelist(nlSetupGrp, 'initial_conditions')
    
        
    if (nspecies > 0) then
      allocate(cnc_normal(nspecies), stat=stat)
      if (stat /= 0) then
        call set_error('Allocate failed', 'setup')
        return
      end if
      call msg('')
      call msg('Setting initial conditions for normal species')
      call msg('')
      call initial_conditions(nl, species_normal, cnc_normal, nspecies)
    else
      call set_error('No species', 'setup')
      return
    end if
    
    if (error) return

    if (nspecies_shortlived > 0) then
      allocate(cnc_shortlived(nspecies_shortlived), stat=stat)
      if (stat /= 0) then
        call set_error('Allocate failed', 'setup')
        return
      end if
      call msg('')
      call msg('Setting initial conditions for shortlived species')
      call msg('')
      call initial_conditions(nl, species_shortlived, cnc_shortlived, nspecies_shortlived)
    end if

    if (error) return

    if (nspecies_aerosol > 0) then
      allocate(cnc_aerosol(nspecies_aerosol), stat=stat)
      if (stat /= 0) then
        call set_error('Allocate failed', 'setup')
        return
      end if
      call msg('')
      call msg('Setting initial conditions for aerosol species')
      call msg('')
      call initial_conditions(nl, species_aerosol, cnc_aerosol, nspecies_aerosol)

      call init_aerosol_distr(nl, species_normal, cnc_normal, nspecies, species_aerosol, &
                            & cnc_aerosol, nspecies_aerosol)


    end if

    if (error) return

    !call free_work_array(q_met_dyn)
    !call free_work_array(q_met_stat)

    allocate(chemRules%low_mass_trsh(nspecies))
    chemRules%low_mass_trsh = 0.

    box_setup%output_unit_cnc = fu_next_free_unit()
    filename = 'box_cnc'//trim(box_setup%case_name)//'.out'
    open(box_setup%output_unit_cnc, file=filename, iostat=stat)
    if (stat /= 0) then
      call set_error('Failed to open output file', 'setup')
      return
    end if

    write(box_setup%output_unit_cnc, fmt='(1000(A,1x))') &
                      & ' . ',(trim(fu_str(species_normal(iTmp))), iTmp=1,nSpecies), &
                      & (trim(fu_str(species_shortlived(iTmp))), iTmp=1,nSpecies_shortlived)

    box_setup%output_unit_nbr = fu_next_free_unit()
    filename = 'box_nbr'//trim(box_setup%case_name)//'.out'
    open(box_setup%output_unit_nbr, file=filename, iostat=stat)
    if (stat /= 0) then
      call set_error('Failed to open output file', 'setup')
      return
    end if
    if (nSpecies_aerosol > 0) then
      write(box_setup%output_unit_nbr,fmt='(15x,1000G15.6)') &
           !& (fu_mean_D(species_aerosol(iTmp))*1e6, iTmp=1,nSpecies_aerosol)    !fu_mean_d is obsolite
           & (fu_massmean_D(species_aerosol(iTmp))*1e6, iTmp=1,nSpecies_aerosol) !RH: CHECK: fu_massmean_d or fu_nominal_d
      write(box_setup%output_unit_nbr,fmt='(15x,1000G15.6)') &
           & ((fu_max_D(species_aerosol(iTmp)) - fu_min_D(species_aerosol(iTmp)))*1e6, &
            & iTmp=1,nSpecies_aerosol)
    end if

    if (.not. fu_interval_positive(box_setup%timestep)) call init_tangent_linear(box_setup)

  end subroutine setup

    
  subroutine expand_paths(nlptr, setupdir)
    ! Process the (possible) filepaths in the setup namelist:
    ! replace environment variables
    ! replace ^ with the directory setupfile is in
    implicit none
    type(Tsilam_namelist), pointer :: nlptr
    character(len=*), intent(in) :: setupdir

    integer :: ind_item, num_items
    type(Tsilam_nl_item_ptr), pointer :: item
    character(len=worksize_string) :: content

    if (fu_fails(associated(nlptr), 'nlptr not associated', 'expand_paths')) return
    call msg('Expanding paths: ' // trim(setupdir))

    do ind_item = 1, fu_nbr_of_items(nlptr)
      item => fu_get_item(nlptr, ind_item)
      content = fu_content(item)
      if (error) return
      content = fu_extend_grads_hat_dir(content, setupdir)
      if (error) return
      call replace_namelist_item(item, "", fu_name(item), content)
      if (error) return
    end do

  end subroutine expand_paths


  !************************************************************************************

  subroutine init_tangent_linear(setup)
    implicit none
    type(Tbox_setup), intent(in) :: setup
    
    integer :: num_steps, stat

    num_steps = (setup%end_time - setup%start_time) / setup%timestep + 1
  
    if (fu_fails(num_steps > 0, 'num_steps not positive', 'init_tangent_linear')) return

    allocate(tan_lin_traj(nspecies, num_steps), h_start_arr(num_steps), stat=stat)
    if (fu_fails(stat == 0, 'Allocate failed', 'init_tangent_linear')) return
    tan_lin_traj = 0.0
  end subroutine init_tangent_linear

  !************************************************************************************

  subroutine run(box_setup, chemRules)
    implicit none
    type(Tbox_setup) :: box_setup
    type(Tchem_rules), intent(in) :: chemrules

    type(silja_time) :: now
    integer :: iTransf, step_count, num_steps, ind_press
    real :: seconds, cell_volume, lat, lon, zenith_cos, h_start, fixed_albedo
    real, dimension(:), pointer :: garb_array, low_mass_thres
    real(sp), dimension(size(cnc_normal)) :: lin_pt
    logical, parameter :: no_reporting = .false.
    logical :: forward, print_it
    real, dimension(1000,1):: photorates

    lat = box_setup%lat
    lon = box_setup%lon
    cell_volume = 1.0
    h_start = -1
    now = box_setup%start_time
    seconds = fu_sec(box_setup%timestep)

    garb_array => fu_work_array()
    low_mass_thres => fu_work_array()
    low_mass_thres(1:nspecies) = 1e-15
    garb_array(1:nspecies) = 0.0
    forward = fu_interval_positive(box_setup%timestep)

    step_count = 0
    num_steps = size(tan_lin_traj, 2)
    h_start = seconds
    

    if (chemRules%useDynamicAlbedo) then 
       fixed_albedo = real_missing  !!Use it from meteo input
    else
       fixed_albedo = chemRules%defaultStaticAlbedo
    end if  
    
    if (allocated(cnc_shortlived)) cnc_shortlived = 0.0
    !do while (now <= box_setup%end_time)
    do while (fu_between_times(now, box_setup%start_time, box_setup%end_time, accept_boundaries_too=.true.))
      call output(box_setup%output_unit_cnc, cnc_normal, nspecies, cnc_shortlived, nSpecies_shortlived, &
                & fu_sec(now-box_setup%start_time))
     
      call SolarSetup(now)
      step_count = step_count + 1
      if (nspecies_aerosol > 0) then
        call msg('write nbr')
        call output(box_setup%output_unit_nbr, cnc_aerosol, nspecies_aerosol, cnc_shortlived, nSpecies_shortlived, &
                  & fu_sec(now-box_setup%start_time))
      end if
      zenith_cos = fu_solar_zenith_angle_cos(lon, lat, now)
      
     
      ! Aerosol dynamics, adjoint
      !
!!$      if (.not.forward) then
!!$        do iTransf = 1, chemRules%nAerosolDynamics
!!$          select case(chemrules%iAerosolDynTypes(iTransf))
!!$          case(aerosol_dynamics_basic)
!!$            call transform_AerDynBasic(chemRules%rulesAerDynBasic, &
!!$                                     & cnc_normal, &
!!$                                     & cnc_shortlived, &
!!$                                     & cnc_aerosol, &
!!$                                     & chemRules%low_mass_trsh, &
!!$                                     & garb_array(:), &
!!$                                     & chemRules%ChemRunSetup%mapVolume2NumberCnc, &
!!$                                     & metdat, &
!!$                                     & seconds)
!!$            
!!$            
!!$          case(aerosol_dynamics_simple)
!!$            call transform_AerDynSimple(cnc_normal, &
!!$                                      & cnc_shortlived, &
!!$                                      & garb_array, &
!!$                                      & chemRules%rulesAerDynSimple, &
!!$                                      & low_mass_thres, &
!!$                                      & metdat, seconds)
!!$            
!!$            
!!$          case default
!!$            continue
!!$          end select
!!$        enddo
!!$      end if
      
      ! Gas phase
      !

      do iTransf = 1, chemrules%nTransformations
        !exit
        select case(chemrules%iTransformTypes(itransf))
        case(transformation_passive)
          call transform_passive(cnc_normal, &
                               & chemRules%rulesPassive, &
                               & metdat(:,1), &
                               & 1, &
                               & seconds, &
                               & print_it)
          
        case (transformation_sulphur_dmat)
          call transform_dmat(cnc_normal, &
                            & cnc_shortlived, &
                            & chemRules%rulesSulphurDMAT, &
                            & metdat(:,1), &
                            & zenith_cos, &
                            & now, lat, lon, &
                            & seconds, low_mass_thres, garb_array(:), &
                            & print_it)

        case (transformation_acid_basic)
          
          call transform_acid_Basic(cnc_normal, &
                                  & cnc_shortlived, &
                                  & chemRules%rulesAcidBasic, &
                                  & metdat(:,1),&
                                  & seconds, &
                                  & garb_array(:), &
                                  & zenith_cos, &
                                  & species_normal, low_mass_thres, &
                                  & no_reporting, 1, 1, 1, &
                                  & lat, lon, now, &
                                  & print_it)
          
        case (transformation_cbm4)
          ind_press = fu_index_in_meteo_input(meteo_input, pressure_flag)
          if (fu_fails(ind_press > 0, 'Failed to get pressure', 'run')) return
          zenith_cos = fu_solar_zenith_angle_cos(lon, lat, now)
          !call get_photorates_column(lon, lat, now, zenith_cos, (/metdat(ind_press)/), photorates)
          call prepare_step_cbm4() !.false.)
          if (forward) then
            if (allocated(tan_lin_traj)) then 
              tan_lin_traj(:,step_count) = cnc_normal
              h_start_arr(step_count) = h_start
            end if
            call msg('hstart::', h_start)
            call transform_cbm4(cnc_normal, &
!                              & photorates(:,1), &
                              & chemRules%rulesCBM4, &
                              & metdat(:,1), &
                              & seconds, &
                              & garb_array(:), &
                              & zenith_cos, h_start, print_it)
          else
            !print *, 'tlt_adj', tan_lin_traj(:,step_count)
            !print *, 'tlt_ajd_sum', sum(tan_lin_traj(:,step_count))
            lin_pt(1:nspecies) = tan_lin_traj(:,num_steps - step_count + 1)
            h_start = h_start_arr(step_count)
            call transform_cbm4_adj(cnc_normal, lin_pt(1:nspecies), &
!                                  & photorates(:,1), &
                                  & chemRules%rulesCBM4, &
                                  & metdat(:,1), &
                                  & seconds, &
                                  & garb_array(:), &
                                  & zenith_cos, h_start, print_it)
          end if

        case (transformation_radioactive)
          if (forward) then
            call transform_radioactive(cnc_normal, &
                                     & chemRules%rulesRadioactive, &
                                     & metdat(:,1), &
                                     & seconds, &
                                     & print_it)
          else
            call set_error('Adjoint radioactive decay is not available', 'run')
          end if

        case (transformation_cbm42_strato)
          ind_press = fu_index_in_meteo_input(meteo_input, pressure_flag)
          if (fu_fails(ind_press > 0, 'Failed to get pressure', 'run')) return
          zenith_cos = fu_solar_zenith_angle_cos(lon, lat, now)
          !call get_photorates_column(lon, lat, now, zenith_cos, (/metdat(ind_press)/), photorates) !Old buggy
          !call get_photorates_column(metdat_col, zenith_cos, fixed_albedo, now, aodext, & !copy from chemistry_manager
          !       & aodscat, photorates, chemRules%ifPhotoAOD)
          !metdat should be two dimensions array:
          call get_photorates_column(metdat, zenith_cos, fixed_albedo, now, (/1.0/), &
                 & (/0.9/), (/300./), photorates, chemRules%ifPhotoAOD, standard_atmosphere, (/0./), 0.999,fake_cloud) !CHECK

              !  subroutine  get_photorates_column(metdat_col, zenith_cos, alb_sfc_fixed, now, &
         !!& aod_ext, aod_scat, o3_col, rates, ifPhotoAOD, PhotoO3colType, tau_above_bott, ssa, cloud_model)


          if (error) return
          call prepare_step_cbm42_strato(.false., photorates(:,1))
          if (forward) then
            if (allocated(tan_lin_traj)) then 
              tan_lin_traj(:,step_count) = cnc_normal
              h_start_arr(step_count) = h_start
            end if
            call transform_cbm42_strato(cnc_normal, &
                                      & photorates(:,1), &
                                      & chemRules%rulesCBM42_strato, &
                                      & metdat(:,1), &
                                      & seconds, &
                                      & garb_array(:), &
                                      & zenith_cos, lat, &
                                      & h_start, now, print_it)
                        
          else
            lin_pt(1:nspecies) = tan_lin_traj(:,num_steps - step_count + 1)
            h_start = h_start_arr(num_steps - step_count + 1)
            call transform_cbm42_strato_adj(cnc_normal, lin_pt(1:nspecies), &
                                          & photorates(:,1), &
                                          & chemRules%rulesCBM42_strato, &
                                          & metdat(:,1), &
                                          & seconds, &
                                          & garb_array(:), &
                                          & zenith_cos, lat, h_start, now, print_it)
          end if
          
          
          
        case (transformation_cbm5_SOA)
            zenith_cos = fu_solar_zenith_angle_cos(lon, lat, now)
            ind_press = fu_index_in_meteo_input(meteo_input, pressure_flag)
            !call get_photorates_column(lon, lat, now, zenith_cos, (/metdat(ind_press)/), photorates) !Old buggy
            !call get_photorates_column(metdat_col, zenith_cos, fixed_albedo, now, aodext, & !copy from chemistry_manager
            !     & aodscat, photorates, chemRules%ifPhotoAOD)
            !metdat should be two dimensions array:
            call get_photorates_column(metdat, zenith_cos, fixed_albedo, now, (/1.0/), &
                   & (/0.9/), (/300./), photorates, chemRules%ifPhotoAOD, standard_atmosphere, (/0./), 0.999,fake_cloud) !CHECK
            call prepare_step_cbm5_soa(.false., photorates(:,1))  !! Sets constant rates, reports others
            if (forward) then
                if (allocated(tan_lin_traj)) then 
                  tan_lin_traj(:,step_count) = cnc_normal
                  h_start_arr(step_count) = h_start
                end if
               call msg('hstart::', h_start)

               call transform_cbm5_SOA(cnc_normal, &
                                     & photorates(:,1), &
                                     & chemRules%rulesCBM5_SOA, &
                                     & metdat(:,1), &
                                     & seconds, &
                                     & garb_array(:), &
                                     & zenith_cos, lat, &
                                     & h_start, &
                                     & now, &
                                     & print_it)
            else
              lin_pt(1:nspecies) = tan_lin_traj(:,num_steps - step_count + 1)
              h_start = h_start_arr(num_steps - step_count + 1)
              call transform_cbm5_SOA_adj(cnc_normal, lin_pt(1:nspecies), &
                                            & photorates(:,1), &
                                            & chemRules%rulesCBM5_SOA, &
                                            & metdat(:,1), &
                                            & seconds, &
                                            & garb_array(:), &
                                            & zenith_cos, lat, h_start, now, print_it)
            end if
        end select
      end do

      if (forward) then
        do iTransf = 1, chemRules%nAerosolDynamics
          select case(chemrules%iAerosolDynTypes(iTransf))
          case(aerosol_dynamics_basic)
            call transform_AerDynBasic(chemRules%rulesAerDynBasic, &
                                     & cnc_normal, &
                                     & cnc_shortlived, &
                                     & cnc_aerosol, &
                                     & chemRules%low_mass_trsh, &
                                     & garb_array(:), &
                                     & chemRules%ChemRunSetup%mapVolume2NumberCnc, &
                                     & metdat(:,1), &
                                     & seconds, print_it)
                      
          case(aerosol_dynamics_simple)
            call transform_AerDynSimple(cnc_normal, &
                                      & cnc_shortlived, &
                                      & garb_array, &
                                      & chemRules%rulesAerDynSimple, &
                                      & low_mass_thres, &
                                      & metdat(:,1), seconds, zenith_cos, print_it)
            
          case(aerosol_dynamics_Mid_Atmosph)
            call transform_AerDynMidAtm(cnc_normal, &
                                      & cnc_shortlived, &
                                      & cnc_aerosol, &
                                      & garb_array, &
                                      & chemRules%rulesAerDynMidAtmosph, &
                                      & low_mass_thres, &
                                      & metdat(:,1), &
                                      & seconds, &
                                      & nSpecies, species_normal, 1,1,1, &
                                      & print_it)
            
          case(aerosol_dynamics_VBS)
              call transform_AerDynVBS(cnc_normal, &
                                    & cnc_shortlived, &
                                    & garb_array, &
                                    & chemRules%rulesAerDynVBS, low_mass_thres, &
                                    & metdat(:,1), seconds, print_it)
            
          case default
            continue
          end select
        enddo
      end if
      
      if(error)return

      now = now + box_setup%timestep
    end do
    call free_work_array(low_mass_thres)
    call free_work_array(garb_array)
  end subroutine run

  !************************************************************************************

  subroutine output(unit, values, nvalues, values_sl, nvalues_sl, seconds)
    implicit none
    integer, intent(in) :: unit
    integer, intent(in) :: nvalues, nValues_sl
    real(sp), dimension(:), intent(in) :: values, values_sl
    real, intent(in) :: seconds

    write(unit, fmt='(f15.1, 1000G15.6)') seconds, values(1:nvalues), values_sl(1:nvalues_sl)
        
  end subroutine output

  !************************************************************************************
  
  subroutine adjoint_test(setup_fwd, chemrules)
    implicit none
    type(Tbox_setup), intent(in) :: setup_fwd
    type(Tchem_rules), intent(in) :: chemrules

    type(Tbox_setup) :: setup_adj
    integer :: ispecies
    real :: delta = 1e-9
    real(sp), dimension(:), allocatable :: cnc1, cnc2, cnc_plus, cnc_minus
    real(sp) :: cost1, cost0, cost_apr
    real(sp) :: deriv, adj_deriv
    integer :: ispecies1, ispecies2

    if(fu_fails(fu_interval_positive(setup_fwd%timestep), 'Need forward setup', 'adjoint_test')) return

    setup_adj = setup_fwd
    setup_adj%timestep = fu_opposite(setup_fwd%timestep)
    setup_adj%end_time = setup_fwd%start_time
    setup_adj%start_time = setup_fwd%end_time
    
    allocate(cnc1(nspecies), cnc_plus(nspecies), cnc_minus(nspecies))
    
    cnc1 = cnc_normal + delta
    call init_tangent_linear(setup_adj)
    
    do while (delta > 1e-18)
      do ispecies1 = 1, nspecies
        !if (ispecies1 /= 6) cycle
        delta = max(1e-10, 1e-3*cnc1(ispecies1))
        cnc_normal = cnc1
        cnc_normal(ispecies1) = cnc1(ispecies1) + delta
        call run(setup_fwd, chemrules)
        cnc_plus = cnc_normal
        
        cnc_normal = cnc1
        cnc_normal(ispecies1) = cnc1(ispecies1) - delta
        call run(setup_fwd, chemrules)
        cnc_minus = cnc_normal
        
        cnc_normal = cnc1
        call run(setup_fwd, chemrules) ! linearization
        do ispecies2 = 1, nspecies
          !if (ispecies2 /= 14) cycle
          deriv = (cnc_plus(ispecies2) - cnc_minus(ispecies2)) / (2*delta)

            cnc_normal = 0.0
            cnc_normal(ispecies2) = 1.0
            call run(setup_adj, chemrules)
            adj_deriv = cnc_normal(ispecies1)

            write (*,"(3(A10, X), 3(E13.5, X))") "adj_test", trim(fu_str(species_normal(ispecies1))), &
                & trim(fu_str(species_normal(ispecies2))), deriv, adj_deriv, cnc_plus(ispecies2)
        end do
      end do
      exit
    end do
    
    
    
  end subroutine adjoint_test

  !************************************************************************************
  
  ! Replicated from io_server to speed up compilation!
  function fu_species_output_name(species) result(chNm)
    !
    ! Generates the output name from the substance name, aerosol mode and 
    ! cocktail template
    !
    ! 9/09: The modes now always range from 1 to nModes(iSubst). The
    ! first can be gas, but has to be checked from the species.
    ! 
    implicit none

    ! Output
    character(len=clen) :: chNm

    ! Imported parameters
    type(silam_species), intent(in) :: species
    
    if(.not. defined(species))then
      chNm = 'undefined'
      return
    endif
    !
    ! Component one: substance name
    !
    chNm = fu_substance_name(species)
    !
    ! Component two: aerosol mode
    !
    if(defined(species%mode))then

      if (species%mode == in_gas_phase) then
        chNm = chNm + '_gas'
      else 
        !
        ! Within the aerosol mode index range, this is an aerosol mode
        !
        !chNm = chNm + '_m' + fu_aerosol_mode_size_to_str(fu_mean_d(species%mode))    !fu_mean_d is obsolite
        chNm = chNm + '_m' + fu_aerosol_mode_size_to_str(fu_massmean_d(species%mode)) !RH: CHECK: fu_massmean_d or fu_nominal_d
      end if

    else
      ! Would be preferable to use an actual mode number. 
      call msg_warning('Why undefined Mode?', 'fu_species_output_name')
    endif

    if(.not. (species%waveLength .eps. real_missing))then
      chNm = chNm + '_w' + fu_optical_wave_length_to_str(species%wavelength)
    endif

  end function fu_species_output_name

end program box
