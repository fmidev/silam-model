MODULE source_terms_NH3_FANv2
  !
  ! This module links the NH3 emission model of J.Vira, Flow of Agriculture Nitrogen,
  ! FANv2. The model is described by Vira et al, GMD 2020 and Vira et al, 2022 ACP.
  ! The current module is an interface to FAN, which is a separate module FanMod.f90.
  ! In-essence, this module is a modified copy of Fan2CTSMMod.f90 interface module.
  !
  ! All units: SI, unless otherwise stated
  !
  ! Original code: J.Vira, Fan2CTSMMod.f90
  ! Author of the SILAM adaptation: Mikhail Sofiev mikhail.sofiev@fmi.fi
  !
  ! Language: ANSI FORTRAN-90 (or close to it)
  !
  use source_terms_time_params !cocktail_basic
  use FanMod 

  implicit none
  private

  !
  ! PUBLIC routines of sea salt
  !
  public fill_NH3_src_from_namelist
  public reserve_NH3_source
  public init_emission_NH3
  public create_source_containing_grid
  public source_2_second_grid
  public add_source_species_NH3_src
  public add_input_needs
  public link_source_to_species
  public prepare_inject_NH3_src
  public compute_emission_NH3_src
  public fu_NH3_emis_owned_quantity
  public fu_name
  public fu_source_nbr
  public fu_source_id_nbr
  public typical_species_conc
  public report

  !
  ! Private routines of the sea salt source
  !
  private add_input_needs_NH3_src
  private create_src_cont_grd_NH3_src
  private project_NH3_src_second_grd
  private link_NH3_src_to_species
  private fu_source_id_nbr_of_NH3_src
  private fu_source_nbr_of_NH3_src
  private fu_source_name_NH3_src
  private typical_species_cnc_NH3_src
  private report_NH3_src

  !
  ! Private subs of sea salt source
  !
  interface add_input_needs
    module procedure add_input_needs_NH3_src
  end interface

  interface create_source_containing_grid
    module procedure create_src_cont_grd_NH3_src
  end interface

  interface source_2_second_grid
    module procedure project_NH3_src_second_grd
  end interface

  interface link_source_to_species
    module procedure link_NH3_src_to_species
  end interface

  interface fu_source_nbr
    module procedure fu_source_nbr_of_NH3_src
  end interface

  interface fu_source_id_nbr
    module procedure fu_source_id_nbr_of_NH3_src
  end interface

  interface fu_name
    module procedure fu_source_name_NH3_src
  end interface

  interface typical_species_conc
    module procedure typical_species_cnc_NH3_src
  end interface

  interface report
    module procedure report_NH3_src
  end interface

  !
  ! There might be several types of the emission algorithm. Therefore, the below list 
  ! of parameters might eventually grow
  !
  integer, private, parameter :: emis_FANv2 = 4120

  !
  ! The NH3 source term
  !
  TYPE silam_NH3_source
    PRIVATE
    CHARACTER(len=clen) :: src_nm, sector_nm  ! Name of the area source and sector
    character(len=fnlen) :: dataDir           ! main directory for NH3 source metadata
    integer :: emisMethod, src_nbr, id_nbr    ! A source and id numbers in a WHOLE source list
    integer :: nLevsDispVert, nSpecies, nSpeciesNH3
    type(silam_vertical) :: vertLevsDispVert
    real, dimension(:), pointer :: levFractDispVert, fzDisp
    type(Tsilam_namelist), pointer :: nlInputFiles  ! namelist for names of supplementary files
    type(silam_species), dimension(:), pointer :: species
    type(silja_field) :: source_mask
    type(chemical_adaptor) :: adaptor
    
    ! Reduction factor for fertilizer due to mechanical incorporation.
    ! N available for volatilization becomes multiplied by (1-fert_incorp_reduct).
    real(r8) :: fert_incorp_reduct = 0.25_r8
  
    type(silja_logical) :: defined
  END TYPE silam_NH3_source

  type NH3_src_ptr
    type(silam_NH3_source) :: NH3_src
  end type NH3_src_ptr
  public NH3_src_ptr


CONTAINS


  !*********************************************************************

  subroutine fill_NH3_src_from_namelist(nlSetup, srcNH3, expected_species, chDataDir)
    !
    ! Initializes the sea salt source term.
    ! The parameters are read from the given ini file and stored to the returned 
    ! source term. After the stuff has been read, the main FluxPerMode is to be generated.
    ! This configuration allows for several sea salt sources
    !
    implicit none

    ! Imported parameters
    type(Tsilam_namelist), intent(in) :: nlSetup
    type(silam_NH3_source), intent(inout) :: srcNH3
    type(silam_species), dimension(:), intent(in), allocatable :: expected_species
    character(len=*), intent(in) :: chDataDir

    ! Local variables
    integer :: iTmp, iSpecies, nFiles
    type(Tsilam_nl_item_ptr), dimension(:), pointer ::  pItems
    type(silam_species), dimension(1) :: species
    integer, dimension(:), pointer :: indices
    logical :: ifFound

    !
    ! Names
    !
    srcNH3%src_nm = fu_content(nlSetup,'source_name')
    srcNH3%sector_nm = fu_content(nlSetup,'source_sector_name')
    srcNH3%defined = silja_false

    !
    ! Emission index type
    !
    select case(fu_str_u_case(fu_content(nlSetup,'NH3_emission_method')))

      case ('FAN_v2')
        srcNH3%emisMethod = emis_FANv2

      case default
        call set_error('Unknown emission method:' + fu_content(nlSetup,'NH3_emission_method'), &
                     & 'fill_NH3_src_from_namelist')
        return
    end select
    !
    ! Now the list of aerosol modes that will be emitted. To set them up, we will
    ! use the standard aerosol procedure - for the sake of unification.
    ! Note that there are two potentially concurring definitions: one coming from aerosol dynamics
    ! the other - written in the ini file. The first one prevails, if it exists
    !
    call get_source_aer_species(srcNH3%species, srcNH3%nSpecies, &
                              & expected_species, fu_content(nlSetup,'NH3_substance_name'), &
                              & nlSetup)
    if(error)return
    !
    ! Store the input fields for the source features.
    !
    srcNH3%nlInputFiles => fu_create_namelist('NH3_src_supplementary_files')
    if(error)return

    nullify(pItems)
    call get_items(nlSetup, 'supplementary_file', pItems, nFiles)
    if(error)return
    do iTmp = 1, nFiles
      call add_namelist_item(srcNH3%nlInputFiles, 'supplementary_file', fu_content(pItems(iTmp)))
    end do
    !
    ! Note that the source mask cannot be read from the file yet: dispersion grid is undefined
    ! We have to store the datadir instead
    !
    srcNH3%dataDir = chDataDir
    call add_namelist_item(srcNH3%nlInputFiles,  'source_area_mask', &
                         & fu_content(nlSetup,'source_area_mask'))
    if(error)return

    srcNH3%defined = silja_undefined

    call report(srcNH3)

  end subroutine fill_NH3_src_from_namelist


  !*****************************************************************

  subroutine add_input_needs_NH3_src(NH3_src, q_met_dynamic, q_met_static, &
                                            & q_disp_dynamic, q_disp_static)
    !
    ! Returns input needs for the emission/transformation routines. 
    !
    implicit none

    ! Imported parameters
    type(silam_NH3_source), intent(in) :: NH3_src
    integer, dimension(:), intent(inout) :: q_met_dynamic, q_met_static, &
                                          & q_disp_dynamic, q_disp_static

    ! Local variables
    integer :: iTmp

    !
    ! Add needed dynamic quantities. Emission needs wind at 10m and, in the future,
    ! friction velocity. We would also need water temperature and ice fraction
    !
    iTmp = fu_merge_integer_to_array(windspeed_10m_flag, q_met_dynamic)

    ! Finally, a land fraction field must always be in the permanent input
    !
    iTmp = fu_merge_integer_to_array(fraction_of_land_flag, q_met_static)

  end subroutine add_input_needs_NH3_src


  !**************************************************************************

  subroutine reserve_NH3_source(NH3_src, &     ! Src to initialise
                              & iSrcNbr, &      ! Src number in the grand list
                              & iSrcIdNbr)      ! SrcID number
    !
    ! Initialises the source:
    ! - stores the reference information: source number and source ID number
    ! - stores the total number of chemical descriptors that will be stored in the source
    ! - nullifies the source dynamic arrays 
    ! - set a few basic internal variables of the source
    !
    implicit none

    ! Imported parameters
    type(silam_NH3_source), intent(inout) :: NH3_src
    integer, intent(in) :: iSrcNbr, iSrcIdNbr

    !
    ! Nullify the basic variables
    !
    NH3_src%src_nm = ''
    NH3_src%sector_nm = ''

    !
    ! Main source parameters - enough to identify it in the global information list
    !
    NH3_src%src_nbr = iSrcNbr
    NH3_src%id_nbr = iSrcIdNbr
    !
    ! A bit of other stuff
    !
    nullify(NH3_src%fZDisp)
    nullify(NH3_src%levFractDispVert)
    !
    ! Finally, mark the source as incomplete
    !
    NH3_src%defined = silja_false

  end subroutine reserve_NH3_source


  !*********************************************************************

  subroutine init_emission_NH3(srcNH3)
    !
    ! Initializes the sea salt source term.
    ! The parameters are read from the given ini file and stored to the returned 
    ! source term. After the stuff has been read, the main FluxPerMode is to be generated.
    ! This configuration allows for several sea salt sources
    !
    implicit none

    ! Imported parameters
    type(silam_NH3_source), intent(inout) :: srcNH3

    ! Local variables
    integer :: iTmp, iSpecies
    character(len=fnlen) :: strTmp, fName
    type(silja_field_id) :: id
    type(silja_field_id), pointer :: idPtr
    real, dimension(:), pointer :: arPtr
    integer, dimension(2), parameter :: mask_quantities = (/fraction_of_water_flag, fraction_of_land_flag/)

    !
    ! First of all, read the source mask
    !
    strTmp = adjustl(fu_content(srcNH3%nlInputFiles, 'source_area_mask'))
    if(error .or. len_trim(strTmp) < 1)then
      call set_error('Source area mask is absent','init_emission_NH3')
      return
    endif
    id = fu_set_field_id_simple(met_src_missing, int_missing, time_missing, level_missing)
    if(error)return
    call set_grid(id, dispersion_grid)
    if(error)return
    call set_field(id, srcNH3%source_mask, .true.)
    if(error)return

    arPtr => fu_grid_data(srcNH3%source_mask)
    idPtr => fu_id(srcNH3%source_mask)

    fname = fu_process_filepath(strTmp(index(strTmp,' ')+1:),superdir=srcNH3%dataDir)
    do iTmp = 1,2
       call set_quantity(id, mask_quantities(iTmp))
       call get_input_field(fname, &  ! file name
                       & fu_input_file_format(strTmp), &          ! file format
                       & id, &                  ! The id to search
                       & arPtr, &               ! data array
                       & dispersion_grid, &  ! storage grid
                       & iOutside = nearestPoint, &         ! out of grid interpolation
                       & iAccuracy = 5, &
                       & wdr = wdr_missing, & 
                       & ifAcceptSameMonth = .false., &
                       & idOut = idPtr)  ! redefine the id
        if(defined(idPtr)) exit
    enddo
    if(error)return

    if (iTmp > 2) then
      call set_error('Failed to get the source mask','init_emission_NH3')
      return
    endif

    !! Turn fraction_of_land_flag to fraction_of_water_flag if needed
    if (fu_quantity(idPtr) == fraction_of_land_flag) then
        arPtr(:) = 1. - arPtr(:)
       call set_quantity(id, fraction_of_water_flag)
    endif

    srcNH3%defined = silja_true

  end subroutine init_emission_NH3


  !****************************************************************************

  subroutine add_source_species_NH3_src(NH3_src, species_list, nSpecies)
    !
    ! Get the species emitted by this source and make a list of
    ! them. The source must be initialized, of course.
    !
    implicit none
    type(silam_NH3_source), intent(in) :: NH3_src
    type(silam_species), dimension(:), pointer :: species_list
    integer, intent(inout) :: nSpecies

    call addSpecies(species_list, nSpecies, NH3_src%species, NH3_src%nSpecies)

  end subroutine add_source_species_NH3_src


  !*******************************************************************************
  
  subroutine create_src_cont_grd_NH3_src(NH3_src, grid_template, ifVerbose, ifExtended)
    !
    ! Creates the grid that covers the area with active NH3 emission
    ! Since the source is global, this is a void routine.
    !
    implicit none

    ! Imported parameters
    type(silam_NH3_source), intent(in) :: NH3_src
    type(silja_grid), intent(inout) :: grid_template
    logical, intent(in) :: ifVerbose
    logical, intent(out) :: ifExtended

    ! So far nothing to do: this source just covers the dispersion grid
    !
    ifExtended = .false.
    return
    
  end subroutine create_src_cont_grd_NH3_src


  !*****************************************************************

  subroutine project_NH3_src_second_grd(NH3_src, grid, vert_disp, vert_proj, iAccuracy)
    !
    ! There is nothing to project in terms of horizontal grid but the vertical has to be made.
    !
    implicit none

    ! Imported parameters
    type(silam_NH3_source), intent(inout) :: NH3_src
    type(silja_grid), intent(in) :: grid
    type(silam_vertical), intent(in) :: vert_disp, vert_proj
    integer, intent(in) :: iAccuracy

    ! Local variables
    real, dimension(2) :: arTmp
    integer :: i
    type(silam_vertical) :: vertTmp

    if(iAccuracy < 1 .or. iAccuracy > 10)then
      call msg('Accuracy switch must be from 1 to 10, not:',iAccuracy)
      call msg_warning('Accuracy switch must be from 1 to 10','project_sslt_src_second_grd')
!      return
    endif

    if(len_trim(NH3_src%sector_nm) > 0)then
      call msg('Re-projecting NH3 source:' + NH3_src%src_nm +'_' + NH3_src%sector_nm)
    else
      call msg('Re-projecting NH3 source:' + NH3_src%src_nm)
    endif

    !
    ! Now re-project the vertical grid of the source to the given vertical
    ! Since vertical distributions are very poorly known and crude, we do not need
    ! any precise meteo-dependent projection, crude will do the job
    !
    NH3_src%vertLevsDispVert = vert_disp
    allocate(NH3_src%levFractDispVert(fu_NbrOfLevels(vert_disp)), &
           & NH3_src%fzDisp(fu_NbrOfLevels(vert_disp)), stat=i)
    if(fu_fails(i == 0, 'Failed allocation dispersion-vertical fractions', &
                                                         & 'project_NH3_src_second_grd'))return
    NH3_src%levFractDispVert(:) = 0.0
    NH3_src%fzDisp(:) = 0.0

    !
    ! Create the NH3 vertical, which can depend on the emission method
    !
    select case(NH3_src%emisMethod)
      
      case(emis_FANv2)
        !
        ! FANv2 goes from surface to 50 m.
        !
        call set_vertical(fu_set_layer_between_two(layer_btw_2_height, 0.0, 25.0), vertTmp)
        if(error)return
        call add_level(vertTmp, fu_set_layer_between_two(layer_btw_2_height, 25.0, 50.0))
        if(error)return
        arTmp(1) = 0.6
        arTmp(2) = 0.4
      case default
        call set_error('Unknown emission method:'+fu_str(NH3_src%emisMethod),'project_sslt_src_second_grd')
        return
    end select

    call reproject_verticals(vertTmp, arTmp, &                    ! vertical from, fractions from
                           & vert_proj, NH3_src%levFractDispVert, &   ! vertical to, fractions to
                           & NH3_src%fzDisp, NH3_src%nLevsDispVert, & ! mass centres, number of non-zero levels
                           & ifMassCentreInRelUnit=.true.)
    call set_missing(vertTmp, .false.)

  end subroutine project_NH3_src_second_grd


  !**************************************************************************

  subroutine link_NH3_src_to_species(species_list, NH3_src)
    !
    ! Having the cocktails created, we should establish the shortcut links between the 
    ! source descriptors and cocktail. The link goes via descr%iEmisCocktSpeciesMapping 
    ! and  descr%factor_to_basic_unit./ descr%factor_foreign_2_basic_unit
    ! Note that cocktail is "anything" and the only requirement is that it has all the species
    ! existing in the source inventory.
    ! That has to happen in two steps. Firstly, we establish these links using the 
    ! single descriptor per source. Then, these connections are distributed to each time slot
    !
    implicit none

    ! Imported parameters
    type(silam_NH3_source), intent(inout) :: NH3_src
    type(silam_species), dimension(:), pointer :: species_list

    !
    ! Linkage is actually just creation of the chemical adaptor
    !
    call create_adaptor(NH3_src%species, species_list, NH3_src%adaptor)
    
  end subroutine link_NH3_src_to_species


  ! *************************************************************************
  
  subroutine prepare_inject_NH3_src(met_buf)
    !
    ! The subroutine prepares the private module pointers to the fields requested 
    ! for injecting the NH3 sources. 
    !
    implicit none

    ! Imported parameters
    type(Tfield_buffer), pointer :: met_buf

!    ! Local variables
!    integer, dimension(:), pointer :: met_q
!    integer :: iQ, iTmp
    
    call set_error('Not implemented','prepare_inject_NH3_src')
    return
    
!    ! nullify the pointers 
!    nullify(fldBVf)
!    nullify(fldAblHeight)
!    nullify(fldHeight)
!    nullify(fldSrfPressure)
!    nullify(fldT)
!    nullify(fldQ)
!
!    ! Scan the meteo buffer   
!    met_q => met_buf%buffer_quantities    
!    do iQ = 1, size(met_q)
!      if(met_q(iQ) == int_missing)exit
!      if(fu_dimension(met_buf, iQ) == 4)then !4D
!
!        select case(met_q(iQ))
!
!          case(brunt_vaisala_freq_flag)
!            fldBVf => met_buf%p4d(iQ)
!
!          case(height_flag)
!            fldHeight => met_buf%p4d(iQ)
!
!          case(temperature_flag)
!            fldT => met_buf%p4d(iQ)
!
!          case(specific_humidity_flag)
!            fldQ => met_buf%p4d(iQ)
!
!          case default
!            cycle
!            
!        end select
!      else !2D
!
!        select case(met_q(iQ))
!          
!          case(abl_height_m_flag)
!            fldAblHeight => met_buf%p2d(iQ)
!          
!          case(ground_pressure_flag)
!            fldSrfPressure => met_buf%p2d(iQ)
!          
!          case default
!            cycle
!            
!        end select
!      endif
!    enddo
  end subroutine prepare_inject_NH3_src

  
  !**************************************************************************

  subroutine compute_emission_NH3_src(NH3_src, &
                                    & met_buf, disp_buf, & 
                                    & now, &      ! current time
                                    & timestep, & ! model time step
                                    & pHorizInterpMet2DispStruct, ifHorizInterp, &
                                    & ifSpeciesMoment, &
                                    & emisMap, mapCoordX, mapCoordY, mapCoordZ, & ! Output
                                    & fMassInjected)                              ! output
    !
    ! Computes the emission fields for NH3 by calling FAN.
    !
    implicit none

    ! Imported parameters
    type(silam_NH3_source), target, intent(in) :: NH3_src
    type(Tfield_buffer), pointer ::  met_buf, disp_buf  ! meteo and internal field buffers
    type(silja_time), intent(in) :: now           ! current time
    type(silja_interval), intent(in) :: timestep  ! model time step
    type(THorizInterpStruct), intent(in) :: pHorizInterpMet2DispStruct
    logical, intent(in) :: ifHorizInterp, ifSpeciesMoment
    type(Tmass_map), intent(inout) :: emisMap, mapCoordX, mapCoordY, mapCoordZ
    real(r8k), dimension(:),  intent(inout) :: fMassInjected

    ! Local variables
    integer :: indFricVel, indZ0, iLev, indW10m, iNH3, &
             & iMeteo, iDisp, ix, iy, ixMeteo, iyMeteo, iStat, iMode, iWaterTemp, iSalinity, &
             & indU, indV, indHeight
    real :: fTmp, timestep_sec, fCellTotal, u_star, u, v, z0, &
          & windspeed, height, dtdxdyLevFrac
    real, dimension(:), pointer :: ptrXSzDisp, ptrYSzDisp, fMinDArray, fMaxDArray, pSrcMask
    logical, save :: ifFirstTime = .true.
    type(silja_field), pointer :: fldMaskPtr
    type(silam_sp) :: sp
    integer, save :: iCount = 0

    !
    ! First, set the output pointer to the locally stored cocktail_map of emission
    ! intensity and get the temporary array for the white caps fraction.
    !
    pSrcMask => fu_grid_data(NH3_src%source_mask)

!open (50, file='NH3_src.dump',recl=nx_dispersion*ny_dispersion,form='unformatted',access='direct')
!write(50,rec=1)(pSrcMask(ix),ix=1,nx_dispersion*ny_dispersion)
!close (50)
!call msg('nx,ny',nx_dispersion, ny_dispersion)
!stop

    ptrXSzDisp => fu_grid_data(dispersion_cell_x_size_fld)  ! get grid cell size in metres
    ptrYSzDisp => fu_grid_data(dispersion_cell_y_size_fld)

    timestep_sec = abs(fu_sec(timestep))
    if(error)return

    !
    ! Computation depends on method of emission computation
    !
    select case(NH3_src%emisMethod)
      
      case(emis_FANv2)
        !
        ! Basic case
        !
        call set_error('Not implemented','compute_emission_for_NH3')
        return

        indZ0 = fu_index(met_buf, surface_roughness_meteo_flag) 
        if(error .or. indZ0 < 1)then
          call set_error('Failed to find z0 field','compute_emission_for_NH3')
          return
        endif
        indFricVel = fu_index(met_buf, friction_velocity_flag) 
        if(error .or. indFricVel < 1)then
          call set_error('Failed to find friction velocity field','compute_emission_for_NH3')
          return
        endif
        !
        ! Basically, a series of calls of FAN for each grid cell
        !
        do iy = 1, ny_dispersion
          do ix = 1, nx_dispersion

            iDisp = ix+(iy-1)*nx_dispersion

            if(pSrcMask(iDisp) < 0.001)cycle

            iMeteo =  fu_grid_index(nx_meteo, ix, iy,  pHorizInterpMet2DispStruct)

            u_star = fu_get_value(met_buf%p2d(indFricVel), nx_meteo, ix, iy, 1., &
                                & pHorizInterpMet2DispStruct, ifHorizInterp)
            z0 = fu_get_value(met_buf%p2d(indZ0), nx_meteo, ix, iy, 1., &
                            & pHorizInterpMet2DispStruct, ifHorizInterp)
            if(z0 < 1.e-10)then
              z0 = 1.e-6
              if(iCount < 1000)then
                call msg('gust: strange Z0=',z0)
                iCount = iCount + 1
              endif
            endif

            !
            ! Fill-in the emission map
            !
            do iLev = 1, NH3_src%nLevsDispVert
              !
              ! First do the check for the overlap: speed-up
              !
              if(abs(NH3_src%levFractDispVert(iLev)) < 1.0e-5)cycle  ! nothing for this dispersion layer

              fCellTotal = 0.0

              !Common factor for all the fluxes
              dtdxdyLevFrac = timestep_sec * ptrXSzDisp(iDisp) * ptrYSzDisp(iDisp)* NH3_src%levFractDispVert(iLev)
              
              do iNH3 = 1, NH3_src%nSpecies

!                fTmp = ...
                call set_error('Not implemented 2','compute_emission_for_NH3')
                return
                
                fCellTotal = fCellTotal + fTmp !* fCorrectionFactor

                emisMap%arM(NH3_src%adaptor%iSp(iNH3),NH3_src%id_nbr,iLev,ix,iy) = &
                      & emisMap%arM(NH3_src%adaptor%iSp(iNH3),NH3_src%id_nbr,iLev,ix,iy) + fTmp
                fMassInjected(NH3_src%adaptor%iSp(iNH3)) = &
                                               & fMassInjected(NH3_src%adaptor%iSp(iNH3)) + fTmp
                
                if (ifSpeciesMoment) then  ! only vertical moment, horizontal ones are irrelevant
                  mapCoordZ%arm(NH3_src%adaptor%iSp(iNH3),NH3_src%id_nbr, ilev, ix,iy) = &
                       & mapCoordZ%arm(NH3_src%adaptor%iSp(iNH3),NH3_src%id_nbr, ilev, ix,iy) + &
                       & fTmp * NH3_src%fzDisp(iLev)
                end if

                emisMap%ifColumnValid(NH3_src%id_nbr,ix,iy) = .true.
                emisMap%ifGridValid(iLev,NH3_src%id_nbr) = .true.

              end do      ! nSpecies
              if (.not. ifSpeciesMoment) then
                mapCoordZ%arM(1,NH3_src%id_nbr, iLev, ix, iy) = &
                                                    & NH3_src%fzDisp(iLev) * fCellTotal + &
                                                    & mapCoordZ%arM(1,NH3_src%id_nbr, iLev, ix, iy)
              end if
            end do      ! iLev
          end do   ! ix dispersion
        end do   ! iy dispersion

      case default
        call msg('Unknown method for NH3 emission:',NH3_src%emisMethod)
        call set_error('Unknown method for NH3 emission','compute_emission_for_NH3')
        return

    end select  ! emisMethod

    ifFirstTime = .false.

!do iNH3 = 1, NH3_src%nSpecies
!call msg('NH3 species index and emission:' + fu_str(iNH3), &
!             & NH3_src%adaptor%iSp(iNH3), fMassInjected(NH3_src%adaptor%iSp(iNH3)))
!end do
  end subroutine compute_emission_NH3_src


  !********************************************************************************************

  logical function fu_NH3_emis_owned_quantity(NH3_src, quantity)
    !
    ! Checks whether the specific quantity is handled by the emission sources exclusively
    ! Such quantities may be related to emission species somewhere else
    !
    implicit none

    ! Imported parameters
    type(silam_NH3_source), intent(in) :: NH3_src
    integer, intent(in) :: quantity
    !
    ! The NH3 source does not have own quantities yet
    !
    select case(quantity)
      case default
        fu_NH3_emis_owned_quantity = .false.
    end select
  end function fu_NH3_emis_owned_quantity


  !*****************************************************************

  integer function fu_source_id_nbr_of_NH3_src(NH3_src)
    !
    ! Returns the source number. Reason: all sources are enumerated
    ! sequencially, so that the source can, in fact, be refered by its
    ! number without other information. But so far this number is 
    ! copied to the particles in the pollution cloud. 
    !
    ! NOTE. One and only index may be reasonable. The other one MUST be
    ! negative or zero
    !
    implicit none

    ! Imported parameters with intent IN
    type(silam_NH3_source), intent(in) :: NH3_src

    ! Stupidity check
    if(.not.(NH3_src%defined == silja_false))then
      call set_error('Undefined source given','fu_source_id_nbr_of_NH3_src')
      return
    endif
    fu_source_id_nbr_of_NH3_src = NH3_src%id_nbr

  end function fu_source_id_nbr_of_NH3_src



  !*************************************************************************

  integer function fu_source_nbr_of_NH3_src(NH3_src)
    !
    ! Returns the source number. Reason: all sources are enumerated
    ! sequencially, so that the source can, in fact, be refered by its
    ! number without other information. But so far this number is 
    ! copied to the particles in the pollution cloud. 
    !
    ! NOTE. One and only index may be reasonable. The other one MUST be
    ! negative or zero
    !
    implicit none

    ! Imported parameters with intent IN
    type(silam_NH3_source), intent(in) :: NH3_src

    ! Stupidity check: only firmly undefined source returns int_missing
    if(.not. (NH3_src%defined == silja_false))then
      fu_source_nbr_of_NH3_src = NH3_src%src_nbr
    else
      fu_source_nbr_of_NH3_src = int_missing
      call set_error('Undefined source given','fu_source_nbr_of_NH3_src')
      return
    endif

  end function fu_source_nbr_of_NH3_src


  !*************************************************************************

  subroutine typical_species_cnc_NH3_src(NH3_src, species, nSpecies, arConc)
    !
    ! Guesses a typical level of concentration and divides it with the given accuracy factor
    !
    implicit none

    ! Imported parameters
    type(silam_NH3_source), intent(in) :: NH3_src
    type(silam_species), dimension(:), pointer :: species
    integer, intent(out) :: nSpecies
    real, dimension(:), pointer :: arConc

    ! Local variables
    integer :: iSpecies

    real, parameter :: fTypicalMassCnc = 1.0e-9  ! one microgram

    species => NH3_src%species
    nSpecies = NH3_src%nSpecies

  end subroutine typical_species_cnc_NH3_src


  !*************************************************************************

  function fu_source_name_NH3_src(NH3_src)result(chNm)
    implicit none
    type(silam_NH3_source), intent(in) :: NH3_src
    character(len=clen) :: chNm
    chNm = NH3_src%src_nm
  end function fu_source_name_NH3_src


  !*************************************************************************

  subroutine report_NH3_src(NH3_src)
    implicit none
    type(silam_NH3_source), intent(in) :: NH3_src
    integer :: iSpecies
    call msg('------------------ NH3 source report -----------------')
    if(NH3_src%sector_nm /= '')then
      call msg('NH3 source' + NH3_src%src_nm + '_' + NH3_src%sector_nm)
    else
      call msg('NH3 source' + NH3_src%src_nm)
    endif
    call msg('Species:')
    do iSpecies = 1, NH3_src%nSpeciesNH3
      call report(NH3_src%species(iSpecies))
    end do
    select case(NH3_src%emisMethod)
      case(emis_FANv2)
        call msg('FANv2 mechanism')
      case default
        call set_error('Unknown emission mechanism:'+fu_str(NH3_src%emisMethod),'report_NH3_src')
        return
    end select

    
    call msg('------------------ end NH3 source report -----------------')
    

  end subroutine report_NH3_src


  !*************************************************************************
  !*************************************************************************
  !
  ! Private functions computing the NH3 emission fluxes
  !
  !*************************************************************************
  !*************************************************************************


  
END MODULE source_terms_NH3_FANv2


