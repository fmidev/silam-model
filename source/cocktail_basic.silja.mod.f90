MODULE cocktail_basic
  ! 
  ! 9/2009: this module now contains mainly the types Tmass_map and
  ! silam_cocktail, which is essentially an encapsulation of the
  ! chemical structure of the mass_map. The type depends strongly on
  ! the chemical_setup module.
  !
  ! The aerosol type is moved into chemical_setup, while deposition
  ! and scavenging of chemically inactive aerosols should be defined
  ! in the deposition module. Some routines are still here,
  ! commented. Also some thermodiffusion routines are still located here.
  ! 
  ! A number of unused Tmasses and mass_vector procedures has been
  ! removed.
  !
  ! Code owner Mikhail Sofiev, FMI, mikhail.sofiev@fmi.fi
  ! Modifications Sep 2009: Julius Vira, FMI
  !
  use field_buffer
!!!!!  use lagrange_particles
  use tangent_linear

!  use chemical_setup

  implicit none
  

  ! public routines for MASS_MAP
  !
  public set_mass_map
  public set_mass_map_from_basic_par
  public dealloc_mass_map
  public set_missing
  public fu_nbr_of_species
  public emis2transpMap_Euler_speciesMmt
  public emis2transpMap_Euler
  public reset_map_to_val
  public set_map_species
  public fu_material
  public fu_species
  public fu_quantity
  public check_mass_vector
  public flip_moment_and_basic_var_Euler
  public flip_moment_and_basic_var_Lagr
  public many_mass_maps_to_grads_file
  public mass_map_to_grads_file
  public mass_to_mixing_ratio
  public mixing_ratio_to_mass
  public report
  public check_mass_moments_mass_map
  public check_mass_centres_mass_map

  ! public routines for meteo_input
  !
  public define_meteo_input   ! allocates memory and fills-in the quantities and refs to markets
  public pick_meteo_pointers  ! picks the refs to input fields
  public fill_in_meteo_input     ! gets the values for the given point from the input fields
  public fill_in_meteo_input_column  !!!Same for column 
!  public fu_index

  public aerosol_processes_input_needs

  PUBLIC fu_OH_cnc_forced
  
  !--------------------------------------------------------------------
  !
  ! private routines for MASS_MAP
  !
  private fu_nbr_of_species_of_massMap
  private fu_quantity_of_mass_map
  private report_mass_map
  private fu_mass_map_defined  

  private fu_index_in_meteo_input

  ! Interface part
  !
  interface fu_nbr_of_species
    module procedure fu_nbr_of_species_of_massMap
  end interface
 
  interface report
    module procedure report_mass_map
  end interface

  interface fu_quantity
     module procedure fu_quantity_of_mass_map
  end interface

!  interface fu_index
!    module procedure fu_index_in_meteo_input
!  end interface

  interface defined
     module procedure fu_mass_map_defined
  end interface

  type TPassenger
    integer :: nPassengers
    integer, dimension(:), pointer :: iSpeciesGarbage, iSpeciesPass
  end type TPassenger

  !------------------------------------------------------------------------------
  !
  ! The mass map is a 1-4D map, which can be used in ANY application type.
  ! - pollution cloud: (nLagrageanParticles, 1, 1, 1); 
  ! - a single-source 3D grid: (nx,ny,nLevels,1), 
  ! - multi-source 2D grid it comes with (nx,ny,nSources,1) or (nx,ny,1,nSources)
  ! - dry-wet deposition it is (nx,ny,2,nSources) or (nGridCells,2,nSources,1),...
  !
  ! In some cases, this can also carry the out-of-storage information via n*BorderCells,
  ! which tells how many cells are reserved OUTSIDE the main storage cube. E.g., if
  ! nXBorderCells == 1, true dimensions are 0:nx+1, etc.
  !
  ! So, it is like a memory bunch, which can be linked to any 1-, 2- 3- or 4-D
  ! structure keeping the main chemical info via the cocktail-structured arMass
  !
  ! ATTENTION. 
  ! 12.08.2007 in attempt to reduce the memory consumption of the model, the Tmasses
  ! type is eliminated from the arM storage. Instead, the 5-th dimension is introduced: species.
  ! A single species is a substance in specific phase - gas or aerosol of specific size.
  ! Conversion of species to substance-size 2D structure is done via iSp array.
  !
  ! Changes 9/2009: the map currently has no cocktail_template,
  ! although a pointer could be added if needed. A chemical_mapper is
  ! provided for indexing, where the 1d indices correspond to list of
  ! species. That element fundamentally defines the chemical structure
  ! of the mass map. In some cases (advection moment) they can be undefined. 
  ! 
  type Tmass_map
    integer :: quantity                    ! E.g., deposition, mass, moment, etc.
    integer :: nx,ny,n3D,nSrc, nSpecies    !, nSubst  ! Dimensions of %arM
    integer :: nXBorder, nYBorder, n3DBorder, nSrcBorder
    type(silja_grid) :: gridTemplate
    type(silam_vertical) :: vertTemplate
    type(silam_species), dimension(:), pointer :: species
    type(chemical_mapper) :: mapper
    real, dimension(:,:,:,:,:), pointer :: arM  ! the main storage, e.g. (nSpecies,nSrc,nLev,nx,ny)
    logical, dimension(:,:,:), pointer :: ifColumnValid ! .true. if >=1 element is meaningful
    logical, dimension(:,:), pointer :: ifGridValid     ! .true. if >=1 element is meaningful
    type(TPassenger), dimension(:), pointer :: passengers  ! (nSpecies)
    type(silja_logical) :: defined
  end type Tmass_map

  type Tmass_map_ptr
    type(Tmass_map), pointer :: ptrMassMap
  end type Tmass_map_ptr

  type (Tmass_map_ptr), public, parameter :: mass_map_ptr_missing = Tmass_map_ptr(null())

  type (Tmass_map), public, parameter :: mass_map_missing = Tmass_map(&
                        & int_missing, int_missing, int_missing, int_missing, int_missing,&
                        & int_missing, int_missing, int_missing, int_missing, int_missing,&        
                        & grid_missing, vertical_missing, null(), chemical_mapper_missing, &
                        & null(), null(), null(), null(), silja_false )

  !-------------------------------------------------------------------------------------
  !
  ! The mass map can be a moment of general type, which is used via a conservation
  ! equation rather than the transport one.
  ! Here is the complete mapping from the FROM species to TO ones
  !
  type Tmoment_mapping
    integer :: mapping_type           ! advection moment vs centre of mass, number cnc vs mass, etc
    integer :: nMainSpecies           ! the number of related TO species
    real, dimension(:), pointer :: cutoff_thresholds ! in massMapScale
    !type(Tmass_map), pointer :: pMassMapMain, pMassMapScale  ! where the related species come from/to
    type(TspeciesReference), dimension(:), pointer :: refMainSpecies  ! the related species and factors
    type(silja_logical) :: defined
  end type Tmoment_mapping

  !
  ! Types of the moment mapping: 
  ! - between the centre of mass and the advection moment (factor: particle position)
  ! - between the number and mass concentration (factor: mass of a single particle)
  !
  integer, public, parameter :: centre_mass_vs_moment_bulk = 1201
  integer, public, parameter :: volume_vs_number_cnc_flag = 1202
  integer, public, parameter :: centre_mass_vs_moment_species = 1203
  
  integer, public, parameter :: to_moment = 1203    ! direction of conversions
  integer, public, parameter :: from_moment = 1204
  

  !-------------------------------------------------------------------------
  !
  ! Pointers to the meteorological fields, which can be requested
  ! by the corresponding transformation or deposition routines.
  ! There will be a single pool of the meteorological and other input data
  ! that will be delivered to all the transformation and deposition modules
  ! Each of them will have to select their own parts. This is instead of 
  ! prepare_transdep routines in v.4.5.1
  !
  type Tmeteo_input
    integer, dimension(max_quantities) :: quantity, q_type !!! Index in buffer for global nteo input
    integer, dimension(max_quantities) :: idx !!! Index in buffer for global nteo input
    integer :: nQuantities = int_missing                        !!! Index in  global mteo input for local
    type(silja_logical) :: defined = silja_false
  end type Tmeteo_input

  type (Tmeteo_input), parameter :: meteo_input_empty = & 
                & Tmeteo_input(int_missing,int_missing,int_missing,0,silja_true) 

  ! Main reference set for the chemical setup of the run. They key role is to connect different sets of
  ! species, which can be often in the most-general form as many:many
  !
  type TchemicalRunSetup
    type(TspeciesReference), dimension(:), pointer :: &
                               & refEmis2Transp_mass  => null(), & ! from emis to transp mass-species
                               & refEmis2Transp_nbr => null(), &  ! from emis to transp-number species
                               & refsTransp2opt => null(), &      ! from transp to optic species
                               & refsOpt2transp => null()         ! ... and back
    type(Tmoment_mapping), pointer :: mapVolume2NumberCnc => null() ! from aerosol-volume to number cnc
    type(silja_logical) :: defined = silja_false
  end type TchemicalRunSetup

  !
  ! For hybrid simulations with both Lagrangian and Eulerian dynamics and plume-in-grid technology
  ! we have to define the rules for connecting these systems
  !
  type TLagr2Euler_projectionRules
    integer :: criterionType
    real :: LP_age_max
    real :: LP_relative_size_max
  end type TLagr2Euler_projectionRules
  !
  ! For Lagrangian -> Eulerian projection rules, one can imagine several criteria
  !
  integer, parameter, public :: maximum_age_flag = 2301
  integer, parameter, public :: maximum_relative_size_flag = 2302
  integer, parameter, public :: never_project_flag = 2303
  integer, parameter, public :: immediately_project_flag = 2304
  integer, parameter, public :: residence_interval_fraction = 2305

CONTAINS


  !****************************************************************
  !****************************************************************
  !
  !    Routines for the dispersion maps of cocktails
  !
  !****************************************************************
  !****************************************************************

  subroutine dealloc_mass_map(mass_map)
    implicit none
    type(Tmass_map), intent(inout) :: mass_map

    if (fu_true(mass_map%defined)) then
      call set_missing(mass_map%vertTemplate, .true.)
      call set_missing(mass_map%mapper)
      deallocate(mass_map%arm, mass_map%ifColumnValid, mass_map%ifGridValid, mass_map%passengers)
    endif
    mass_map = mass_map_missing
  end subroutine dealloc_mass_map


  !***********************************************************************

  subroutine set_mass_map_from_basic_par(MassMap, quantity,nx,ny,n3D,nSrc,nspecies,nPassengers, &
                                        & nXBorder, nYBorder, n3DBorder, nSrcBorder, val) 
    !
    ! Sets an abstract mass map for the given sizes. Can be used just for 
    ! an arbitrary needs of chemically structured memory
    !
    implicit none
    ! return value
    type(Tmass_map), intent(out) :: MassMap

    ! Imported parameters
    integer, intent(in) :: quantity, nx, ny, n3D, nSrc, &
                         & nXBorder, nYBorder, n3DBorder, nSrcBorder
    integer, dimension(:), intent(in) :: nPassengers  ! (nSpecies)
    integer, intent(in) :: nspecies
    real, intent(in) :: val

    ! Local variables
    integer :: iTmp, iSubst, iMode, iWave

    !
    ! The main map
    !
    MassMap%quantity = quantity
    MassMap%nx = nx
    MassMap%ny = ny
    MassMap%n3D = n3D
    MassMap%nSrc = nSrc
    MassMap%nXBorder = nXBorder
    MassMap%nYBorder = nYBorder
    MassMap%n3DBorder = n3DBorder 
    MassMap%nSrcBorder = nSrcBorder 
    MassMap%gridTemplate = grid_missing
    call set_missing(MassMap%vertTemplate, .true.)
    !
    ! Take care of the fastest-varying dimension: species. 
    ! Once again: numerous allocations lead to trouble due to overhead of each allocation.
    ! Solution: align ALL species (i.e. substance + its phase or size mode) along a single 
    ! dimension. 
    ! Its total size and decoding array are allocated and set here.
    ! At this stage, we have to keep in mind that in many places it is explicitly assumed
    ! that absolute maximum number of species is max_species. Not the most-elegant solution
    ! but have to accept it at least for now.
    !
    if(nSpecies > max_species)then
      call set_error('Too many species:' + fu_str(nSpecies) + &
                   & ', while it cannot exceed max_species=' + fu_str(max_species), &
                   & 'fu_set_mass_map_from_basic_par')
      return
    endif
    allocate(MassMap%species(nspecies), MassMap%passengers(nSpecies), stat=iTmp)
    if(fu_fails(iTmp == 0,'allocate failed for: ' + fu_str(nspecies) + '-species', &
                        & 'fu_set_mass_map_from_basic_par'))return
    MassMap%nspecies = nspecies
    do iTmp = 1, nspecies
      MassMap%species(iTmp) = species_missing
      MassMap%passengers(iTmp)%nPassengers = nPassengers(iTmp)
      if(nPassengers(iTmp) > 0)then
        allocate(MassMap%passengers(iTmp)%iSpeciesGarbage(nPassengers(iTmp)), &
               & MassMap%passengers(iTmp)%iSpeciesPass(nPassengers(iTmp)))
      else
        nullify(MassMap%passengers(iTmp)%iSpeciesGarbage, MassMap%passengers(iTmp)%iSpeciesPass)
      endif
!      ptrMassMap%indTransformType(iTmp) = int_missing
!      ptrMassMap%indDepositionType(iTmp) = int_missing
!      ptrMassMap%indAerDynamicType(iTmp) = int_missing
    end do

    !
    ! Count the total number of coded array iSp3D and allocate the coding tables 
    ! nSubstPerMode and nModes
    !
    if(nx<=0 .or. ny<=0 .or. n3D<=0 .or. nSrc <=0)then
      call msg_warning('Strange parameters of the mass map')
      call msg('nx=', nx)
      call msg('ny=', ny)
      call msg('n3D=', n3D)
      call msg('nSources=', nSrc)
      call set_error('Strange given map size(s)','fu_set_mass_map_from_basic_par')
      return
    endif

    !
    ! Grand memory storage allocation
    ! Allocation takes into account the reserved border lines
    !
    allocate( MassMap%arM(MassMap%nSpecies, &
                          & 1-nSrcBorder:nSrc+nSrcBorder, &
                          & 1-n3DBorder:n3D+n3DBorder, &
                          & 1-nXBorder:nx+nXBorder, &
                          & 1-nYBorder:ny+nYBorder), &
           & MassMap%ifColumnValid(1-nSrcBorder:nSrc+nSrcBorder, &
                                    & 1-nXBorder:nx+nXBorder, 1-nYBorder:ny+nYBorder), &
           &  MassMap%ifGridValid(1-n3DBorder:n3D+n3DBorder, 1-nSrcBorder:nSrc+nSrcBorder), stat=iTmp)
    if(iTmp /= 0)then
      call msg('Failed to allocate memory for map with size:', nx * ny * n3D)
      call set_error('Failed to allocate memory for map','fu_set_mass_map_from_basic_par')
      return
    endif

    MassMap%arM(1:MassMap%nSpecies, &
          & 1-nSrcBorder:nSrc+nSrcBorder, &
          & 1-n3DBorder:n3D+n3DBorder, &
          & 1-nXBorder:nx+nXBorder, &
          & 1-nYBorder:ny+nYBorder) = val
    MassMap%ifColumnValid(1-nSrcBorder:nSrc+nSrcBorder, &
                        & 1-nXBorder:nx+nXBorder, 1-nYBorder:ny+nYBorder) = .not. (val .eps. 0.0)
    MassMap%ifGridValid(1-n3DBorder:n3D+n3DBorder, 1-nSrcBorder:nSrc+nSrcBorder) = &
                                                                          & .not. (val .eps. 0.0)

    MassMap%defined = silja_true

  end subroutine set_mass_map_from_basic_par


  !***********************************************************************

  subroutine set_mass_map(MassMap, quantity, nSrc, nBorderCells, &
                                     & gridTemplate, verticalTemplate, speciesTemplate, &
                                     & nPassengers_, val)
    !
    ! Allocates memory for one complete map of mass vectors. To speed-up
    ! the process we do not call many times the above set_cocktail_mass_vector
    ! but rather repeat it here
    !
    implicit none

    ! return value
    type(Tmass_map), intent(out):: MassMap

    ! Imported parameters
    integer, intent(in) :: quantity, nSrc, nBorderCells
    type(silja_grid), intent(in) :: gridTemplate
    type(silam_vertical), intent(in) :: verticalTemplate
    type(silam_species), dimension(:),  intent(in) :: speciesTemplate
    integer, dimension(:), intent(in), optional :: nPassengers_
    real, intent(in) :: val

    ! Internal variables
    integer :: nx, ny, n3D, nspecies, i
    integer, dimension(:), pointer :: nModes, nWaves, nPassengers

    if(defined(gridTemplate) .and. defined(verticalTemplate) ) then

!     & defined(speciesTemplate(1))) then
!     & defined(speciesTemplate(1)%ptr)) then

      call grid_dimensions(gridTemplate, nx, ny)
      n3D = fu_NbrOfLevels(verticalTemplate)

      ! Number of species is either the number of defined species in the
      ! list, or the size of the list.
      !
     ! nSpecies = fu_index(species_missing, speciesTemplate) - 1
      !if (nspecies == int_missing-1) 
      nspecies = size(speciesTemplate)
      
      if(present(nPassengers_))then
        call set_mass_map_from_basic_par(MassMap, quantity, nx,ny,n3D,nSrc,nspecies, nPassengers_, &
                                                   & nBorderCells, nBorderCells, nBorderCells, 0, &
                                                   & val=val)
      else
        nPassengers => fu_work_int_array(max_species)
        if(error)return
        nPassengers(1:nSpecies) = 0
        call set_mass_map_from_basic_par(MassMap, quantity, nx,ny,n3D,nSrc,nspecies, nPassengers, &
                                                   & nBorderCells, nBorderCells, nBorderCells, 0, &
                                                   & val=val)
        call free_work_array(nPassengers)
      endif

      if(error)return

      MassMap%gridTemplate = gridTemplate
      MassMap%vertTemplate = verticalTemplate
      do i = 1, nspecies
        MassMap%species(i) = speciesTemplate(i) !%ptr
      end do
      call create_chemical_mapper(MassMap%species, MassMap%mapper, nspecies)
      MassMap%defined = silja_true

    else
      call msg('One of the given templates (grid, vertical, first species) is undefined:')
      call report(gridTemplate)
      call report(verticalTemplate)
      call report(speciesTemplate(1))
      call set_error('One of the given templates (grid, vertical, cocktail) is undefined:', &
                   & 'fu_set_mass_map_from_params')
      MassMap%defined = silja_false
    endif

  end subroutine set_mass_map


  !***********************************************************************************

  function fu_quantity_of_mass_map(massmap) result(quantity)
    implicit none
    type(Tmass_map), intent(in) :: massmap
    integer :: quantity
    quantity = massmap%quantity
  end function fu_quantity_of_mass_map


  !**********************************************************************************

  subroutine reset_map_to_val(map_to_reset, val)
    !
    ! Just reset the map to the given value
    !
    implicit none

    ! Imported parameters
    type(Tmass_map), intent(inout) :: map_to_reset
    real, intent(in) :: val

    if(associated(map_to_reset%arM))then
      map_to_reset%arM(:,:,:,:,:) = val
    else
      call set_error('Cannot reset the non-associated mass map','')
    endif

  end subroutine reset_map_to_val


  !**********************************************************************************
  
  subroutine set_map_species(pMassMap, speciesTemplate)
    !
    ! Copies the species to the mass map
    !
    implicit none
    ! Imported parameters
    type(Tmass_map), intent(inout) :: pMassMap
    type(silam_species), dimension(:) :: speciesTemplate

    ! Local variables
    integer :: nSpecies, iSpecies

    ! Number of species is either the number of defined species in the
    ! list, or the size of the list.
    !
    nSpecies = fu_index(species_missing, speciesTemplate) - 1
    if (nspecies == int_missing-1) nspecies = size(speciesTemplate)

    !
    ! Attention! We cannot do anything if the number of species available is not the same 
    ! as the number of coming species
    !
    if(pMassMap%nSpecies /= nSpecies)then
      call msg('Number of species in mass map is not the same as the number of new species:', &
                   & pMassMap%nSpecies, nSpecies)
      call set_error('Number of species in mass map is not the same as the number of new species', &
                   & 'set_map_species')
      return
    endif
    !
    ! Now just copy them
    !
    do iSpecies = 1, nSpecies
      pMassMap%species(iSpecies) = speciesTemplate(iSpecies)
    end do
    !
    ! And create the mapping between species and substance-mode-wavelength. In principle,
    ! this is dangerous operation: the mapper creation involves memory allocation without checking.
    ! However, this function is not supposed to be used regularly and in theory it is applied to 
    ! empty mapper, so no problem
    !
    call create_chemical_mapper(pMassMap%species, pMassMap%mapper, nspecies)

  end subroutine set_map_species


  !**********************************************************************************

  subroutine emis2transpMap_Euler_speciesMmt(mapEmis, mapPx_emis, mapPy_emis, mapPz_emis, &
                                           & mapConc, mapPx_conc, mapPy_conc, mapPz_conc, &
                                           & mapAerosol, &
                                           & ChemRunSetup)
    !
    ! Transfers the mass and momentum from emission mass map to concentration map.
    ! A trick is to use the run chemical setup from materials because emission and
    ! transport species are not necessarily the same.
    !
    ! ATTENTION. 
    ! Emission mapP_emis keeps the momentum (no conversion then) while for concentrations 
    ! we keep centre of masses, so have to be careful here
    ! 
    implicit none

    ! Imported parameters
    type(Tmass_map), intent(in) :: mapEmis, mapPx_emis, mapPy_emis, mapPz_emis
    type(Tmass_map), intent(inout) :: mapConc, mapPx_conc, mapPy_conc, mapPz_conc, &
                              & mapAerosol
    type(TchemicalRunSetup), pointer :: ChemRunSetup

    ! Local variables
    integer :: ix,iy,iz,iSrc,iSpEmis, iref, iSpTransp, iSpTo
    integer :: nCellsEmit
    real :: emission_mass, fraction, mass_cur, mass_new
    type(TspeciesReference), dimension(:), pointer :: references, references_nbr

    !
    ! The only trick is to distribute the emission in accordance with the chemical setup
    !
    !call msg('[NO2] before emission', sum(mapConc%arm(26,1,:,:,:)))
    !call msg('[NO] before emission', sum(mapConc%arm(31,1,:,:,:)))

    references => chemRunSetup%refEmis2Transp_mass
    references_nbr => chemRunSetup%refEmis2Transp_nbr
      
    ! Very rough estimate for emitting cell-species
    ! Do decide if parallel is needed at all... Some 1e5-1e6 is a good number
    nCellsEmit = count(mapEmis%ifGridValid(:,:))* count(mapEmis%ifColumnValid(:,:,:)) / mapEmis%nSrc * mapEmis%nspecies

!!!!    !$OMP PARALLEL IF (nCellsEmit > 100000)  default(shared) &
    !$OMP PARALLEL default(shared) &
    !$OMP  &  private(ix,iy,iz,iSrc,iSpEmis, iref, iSpTransp, iSpTo, emission_mass, &
    !$OMP             & fraction, mass_cur, mass_new)
    !$OMP MASTER
    !$  call msg("Emission to transport in parallel. nthreads="+fu_str(omp_get_num_threads()))
    !$OMP END MASTER
    !$OMP DO COLLAPSE(2)
    do iy = 1, mapEmis%ny
     do ix = 1, mapEmis%nx
       if (any(mapEmis%ifColumnValid(:,ix,iy))) then
          do iz = 1, mapEmis%n3D
           do iSrc = 1,mapEmis%nSrc
             if(.not. (mapEmis%ifGridValid(iz,iSrc) .and. mapEmis%ifColumnValid(iSrc,ix,iy))) cycle ! speedup
               do iSpEmis = 1, mapEmis%nspecies
                 emission_mass = mapEmis%arm(iSpEmis,isrc,iz,ix,iy)
                 
                 if (emission_mass == 0.0) cycle  ! No emission -- true zero 

                 do iref = 1, references(ispEmis)%nRefSpecies
                   iSpTransp = references(ispEmis)%indSpeciesTo(iref)
                   fraction = references(ispEmis)%fract(iref)
                   mass_cur = mapConc%arm(ispTransp,isrc,iz,ix,iy)
                   mass_new = emission_mass*fraction + mass_cur
                   mapConc%arm(ispTransp,isrc,iz,ix,iy) = mass_new
                   mapPx_conc%arm(iSpTransp,isrc,iz,ix,iy) = (mapPx_conc%arm(iSpTransp,isrc,iz,ix,iy) * mass_cur &
                                                            & + mapPx_emis%arm(iSpEmis,isrc,iz,ix,iy) * fraction) &
                                                           & / mass_new
                   mapPy_conc%arm(iSpTransp,isrc,iz,ix,iy) = (mapPy_conc%arm(iSpTransp,isrc,iz,ix,iy) * mass_cur &
                                                            & + mapPy_emis%arm(iSpEmis,isrc,iz,ix,iy) * fraction) &
                                                           & / mass_new
                   mapPz_conc%arm(iSpTransp,isrc,iz,ix,iy) = (mapPz_conc%arm(iSpTransp,isrc,iz,ix,iy) * mass_cur &
                                                            & + mapPz_emis%arm(iSpEmis,isrc,iz,ix,iy) * fraction) &
                                                           & / mass_new

                 end do

                 do iref = 1, references_nbr(iSpEmis)%nRefSpecies
                   iSpTransp = references_nbr(iSpEmis)%indSpeciesTo(iref)
                   mapAerosol%arM(iSpTransp,isrc,iz,ix,iy) = mapAerosol%arM(iSpTransp,isrc,iz,ix,iy) + &
                             & mapEmis%arm(iSpEmis,isrc,iz,ix,iy) * references_nbr(iSpEmis)%fract(iSpTo)
                 end do

               end do  ! iSpEmis
           end do  ! iSrc
          end do  ! iz
        endif
     end do  ! ix
    end do  ! iy
    !$OMP END DO
    !$OMP END PARALLEL


    !call msg('[NO2] after emission', sum(mapConc%arm(26,1,:,:,:)))
    !call msg('[NO] after emission', sum(mapConc%arm(31,1,:,:,:)))

  end subroutine emis2transpMap_Euler_speciesMmt


  !**********************************************************************************

  subroutine emis2transpMap_Euler(mapEmis, mapPx_emis, mapPy_emis, mapPz_emis, &
                                & mapConc, mapPx_conc, mapPy_conc, mapPz_conc, &
                                & mapAerosol, &
                                & ChemRunSetup)
    !
    ! Transfers the mass and momentum from emission mass map to concentration map.
    ! A trick is to use the run chemical setup from materials because emission and
    ! transport substances and species are not necessarily the same.
    !
    ! ATTENTION. 
    ! Emission mapP_emis keeps the moment (no conversion then) while for concentrations 
    ! we keep centre of masses, so have to be careful here
    ! 
    !
    implicit none

    ! Imported parameters
    type(Tmass_map), intent(in) :: mapEmis, mapPx_emis, mapPy_emis, mapPz_emis
    type(Tmass_map), intent(inout) :: mapConc, mapPx_conc, mapPy_conc, mapPz_conc, &
                              & mapAerosol
    type(TchemicalRunSetup), pointer :: ChemRunSetup

    ! Local variables
    integer :: ix,iy,iz,iSrc,iSpEmis, iD, iSpTransp, iSpTo
    real :: fMassEmis !, fMassTransp

    !
    ! The only trick is to distribute the emission in accordance with the chemical setup
    !
    do iy = 1, mapEmis%ny
     do ix = 1, mapEmis%nx
       do iz = 1, mapEmis%n3D
        do iSrc = 1,mapEmis%nSrc
          !
          ! Below, we sum-up the moments of emission and transport. However, the transport map
          ! comes as centre of mass, not moments. Have to make the change now.
          !
          mapPx_conc%arM(1,iSrc,iz,ix,iy) = mapPx_conc%arM(1,iSrc,iz,ix,iy) * &
                                          & sum(mapConc%arM(1:mapConc%nSpecies, isrc, iz, ix, iy))
          mapPy_conc%arM(1,iSrc,iz,ix,iy) = mapPy_conc%arM(1,iSrc,iz,ix,iy) * &
                                          & sum(mapConc%arM(1:mapConc%nSpecies, isrc, iz, ix, iy))
          mapPz_conc%arM(1,iSrc,iz,ix,iy) = mapPz_conc%arM(1,iSrc,iz,ix,iy) * &
                                          & sum(mapConc%arM(1:mapConc%nSpecies, isrc, iz, ix, iy))

          if(.not.(mapEmis%ifGridValid(iz,iSrc) .and. mapEmis%ifColumnValid(iSrc,ix,iy)))cycle  ! speedup


!if(sum(mapEmis%arM(:,iSrc,iz,ix,iy))> 0)then !(1:lstEm%nSpeciesEms,ix,iy,iz,iSrc)) > 0)then
!call msg('Got the mass, ix, real(iy)',ix,real(iy))
!endif
!fMassTransp = sum(mapConc%arM(1:mapConc%nSpecies,iSrc,iz,ix,iy))
!fMassEmis = sum(mapEmis%arM(1:mapEmis%nSpecies,iSrc,iz,ix,iy))

          !
          ! Again speedup and also handling numerics. But here comparison has to be exact, not
          ! approximate (== rathern then .eps.) since the absolute level depends on substance
          !
!          if(fMassEmis == 0.0) cycle

          if(sum(mapEmis%arM(1:mapEmis%nSpecies,iSrc,iz,ix,iy)) == 0.0)cycle

          do iSpEmis = 1, mapEmis%nSpecies
            !
            ! Distribute emitted masses
            !
            do iSpTo = 1, ChemRunSetup%refEmis2Transp_mass(iSpEmis)%nRefSpecies

              iSpTransp = ChemRunSetup%refEmis2Transp_mass(iSpEmis)%indSpeciesTo(iSpTo)

              mapConc%arM(iSpTransp,isrc,iz,ix,iy) = mapConc%arM(iSpTransp,isrc,iz,ix,iy) + &
                                           & mapEmis%arm(iSpEmis,isrc,iz,ix,iy) * &
                                           & ChemRunSetup%refEmis2Transp_mass(iSpEmis)%fract(iSpTo)
            end do
            !
            ! Distribute emitted numbers for aerosols. Note that emission sources
            ! emit both mass and numbers, so no actual conversion is needed here, just distribution
            ! However, the aerosol maps must be in "moment" status (i.e., numbers), not "diameter".
            !
            do iSpTo = 1, ChemRunSetup%refEmis2Transp_nbr(iSpEmis)%nRefSpecies

              iSpTransp = ChemRunSetup%refEmis2Transp_nbr(iSpEmis)%indSpeciesTo(iSpTo)

              mapAerosol%arM(iSpTransp, isrc, iz, ix, iy) = &
                                            & mapAerosol%arM(iSpTransp, isrc, iz, ix, iy) + &
                                            & mapEmis%arm(iSpEmis, isrc, iz, ix, iy) * &
                                            & ChemRunSetup%refEmis2Transp_nbr(iSpEmis)%fract(iSpTo)
            end do
          end do  ! iSpecies
          !
          ! Having masses distributed, sum-up moment and return it on-the-fly back to
          ! center-of-mass coordinate
          !
          mapPx_conc%arM(1,iSrc,iz,ix,iy) = mapPx_conc%arM(1,iSrc,iz,ix,iy) + &
                                          & mapPx_emis%arM(1,iSrc,iz,ix,iy)
          mapPy_conc%arM(1,iSrc,iz,ix,iy) = mapPy_conc%arM(1,iSrc,iz,ix,iy) + &
                                          & mapPy_emis%arM(1,iSrc,iz,ix,iy)
          mapPz_conc%arM(1,iSrc,iz,ix,iy) = mapPz_conc%arM(1,iSrc,iz,ix,iy) + &
                                          & mapPz_emis%arM(1,iSrc,iz,ix,iy)
        end do  ! iSrc
       end do  ! iz
     end do  ! ix
    end do  ! iy
  end subroutine emis2transpMap_Euler


  !**************************************************************************************

  subroutine check_mass_vector(vMass, garbage, species_list, chPlace, fLowCncThresh, nSpecies, ix, iy, iz, print_it)
    !
    ! Scans the mass vector and checks the masses for non-negative values
    !
    implicit none

    ! Imported parameters
    real, dimension(:), intent(inout) :: vMass
    real, dimension(:), intent(inout) :: garbage
    type(silam_species), dimension(:), intent(in) :: species_list
    character(len=*), intent(in) :: chPlace
    real, dimension(:), intent(in) :: fLowCncThresh
    integer, intent(in) :: nSpecies, ix, iy, iz
    logical, intent(out), optional :: print_it

    ! Local variables
    integer :: iSpecies, nSpeciesMax

    do iSpecies = 1, nSpecies
      if(.not.(vMass(iSpecies) .ge.  0.0))then ! Should capture NaNs
        if(-vMass(iSpecies) < fLowCncThresh(iSpecies))then
          garbage(iSpecies) = garbage(iSpecies) + vMass(iSpecies)
          vMass(iSpecies) = 0.
          cycle
        endif
        !$OMP CRITICAL (barkMassVector)
        call msg('Negative mass found at:' + chPlace)
        call report(species_list(iSpecies))
        call msg('Mass itself', vMass(iSpecies))
        call set_error('Negative mass found at:' + chPlace,'check_mass_vector')
        !$OMP END CRITICAL (barkMassVector)
      endif
    end do  ! iSpecies
    !
    ! If negative mass is found in the cell, report it in more details
    ! Note that error is shared OMP variable, so let's set private print_it as well
    !
    if(error)then
      !$OMP CRITICAL (barkMassVector)
      if(present(print_it)) print_it = .true.
      call msg('Place of error (ix,iy,iz):',(/ix,iy,iz/))
      do iSpecies = 1, nSpecies
        call msg('Mass:'+ fu_str(species_list(iSpecies)), vMass(iSpecies))
      end do
      !$OMP END CRITICAL (barkMassVector)
    else
      if(present(print_it)) print_it = .false.
    endif  ! error

  end subroutine check_mass_vector

  !************************************************************************************

  subroutine create_moment_2_cm_mapping(moment_map, mass_map, &
                                      & low_mass_thr, mapping_type, mapping)
    implicit none
    type(Tmass_map), target, intent(in) :: moment_map, mass_map
    real, dimension(:), target, intent(in) :: low_mass_thr
    integer, intent(in) :: mapping_type
    type(Tmoment_mapping), intent(out) :: mapping

    integer :: status

    mapping%nMainSpecies = moment_map%nspecies
    mapping%mapping_type = mapping_type
    mapping%cutoff_thresholds => low_mass_thr

    if (.not. fu_true(moment_map%defined)) then
      call set_error('Moment map not defined', 'create_moment_2_cm_mapping')
      return
    end if
    
    select case(mapping_type) 
    case (centre_mass_vs_moment_species)
      ! 1-to-1
      !
      allocate(mapping%refMainSpecies(mass_map%nspecies),  stat=status)
      if (status /= 0) then
        call set_error('Allocate failed', 'create_moment_2_cm_mapping')
        return
      end if
      call set_speciesReferences(mass_map%species, mass_map%species, mapping%refMainSpecies)
      
    case (centre_mass_vs_moment_bulk)
      ! 1-to-all
      !
      allocate(mapping%refMainSpecies(1), stat=status)
      if (status /= 0) then
        call set_error('Allocate failed', 'create_moment_2_cm_mapping')
        return
      end if     
      call set_speciesReferences((/mass_map%species(1)/), mass_map%species, mapping%refMainSpecies, &
                               & .true., .true.) ! ignore mode & substance

    case default 
      call set_error('Invalid mapping type', 'create_moment_2_cm_mapping')
      
    end select
    
    if (error) return

    mapping%defined = silja_true

  end subroutine create_moment_2_cm_mapping
  

  !****************************************************************************

  subroutine create_mass_2_nbr_mapping(massMapTransport, &    ! transport
                                     & massMapAerosol, &      ! aerosol
                                     & mappingNbr2Volume)  ! mapping of number <-> volume
    !
    ! In case the aerosol pecies are considered carefully, this function will calculate 
    ! the ways to flip from the number particles to their volumes. This corresponds
    ! to switch between the variable, such as particle diameter, and its mass-weighted 
    ! moment, such as the number concentration
    !
    ! The idea is based on the physical meaning of the number concentrations: they are made
    ! for each mode summing-up all the species in that mode.
    !
    implicit none

    ! Imported parameters
    type(Tmass_map), intent(inout) :: massMapTransport, massMapAerosol
    type(Tmoment_mapping), intent(inout) :: mappingNbr2Volume

    ! Local variables
    integer :: iSpAerosol, iSpTransport
    type(silam_material), pointer :: materialPtrTmp
    integer, dimension(max_species) :: indicesTranspSpecies

    if(.not. defined(massMapAerosol))then
      call set_error('Aerosol mass map is not initiated','create_mass_2_nbr_mapping')
      return
    endif
    if(error)return

    !
    ! Basic variables of the mapping
    !
    mappingNbr2Volume%defined = silja_false
    mappingNbr2Volume%mapping_type = volume_vs_number_cnc_flag
    mappingNbr2Volume%nMainSpecies = massMapAerosol%nSpecies
    !mappingNbr2Volume%pMassMapMain => massMapAerosol
    !mappingNbr2Volume%pMassMapScale => massMapTransport
    allocate(mappingNbr2Volume%refMainSpecies(massMapAerosol%nSpecies), stat=iSpAerosol)
    if(iSpAerosol /= 0)then
      call set_error('failed allocation of main species array','create_mass_2_nbr_mapping')
      return
    endif

    !
    ! Go through the aerosol species, get all the transport species they consist of
    ! and set the moment_mapping structure
    !
    do iSpAerosol = 1, massMapAerosol%nSpecies
      !
      ! Find the species in the transport mass map that have the same mode
      !
      call select_species(massMapTransport%species, &     ! transport species
                        & massMapTransport%nSpecies, &     ! nbr of transport species
                        & char_missing, &                 ! substance name: take any
                        & fu_mode(massMapAerosol%species(iSpAerosol)), &  ! mode to find
                        & real_missing, &                 ! wave length: take any
                        & indicesTranspSpecies, &          ! out: species found
                        & mappingNbr2Volume%refMainSpecies(iSpAerosol)%nRefSpecies) ! out: nbr of found species
      if(error)return
      !
      ! Create the space for the found references
      !
      allocate(mappingNbr2Volume%refMainSpecies(iSpAerosol)%indSpeciesTo( &
                                     & mappingNbr2Volume%refMainSpecies(iSpAerosol)%nRefSpecies), &
             & mappingNbr2Volume%refMainSpecies(iSpAerosol)%fract( &
                                     & mappingNbr2Volume%refMainSpecies(iSpAerosol)%nRefSpecies), &
             & stat = iSpTransport)
      if(iSpTransport /= 0)then
        call set_error('Failed allocation of the refernce species','create_mass_2_nbr_mapping')
        return
      endif

      !
      ! Find out the conversion factor between the number concentration and the volume of a 
      ! single particle
      !
      do iSpTransport = 1, mappingNbr2Volume%refMainSpecies(iSpAerosol)%nRefSpecies
        !
        ! Each reference is the index of the transport species and the conversion factor to be used
        ! for this very species.
        !
        materialPtrTmp => massMapTransport%species(indicesTranspSpecies(iSpTransport))%material
        mappingNbr2Volume%refMainSpecies(iSpAerosol)%indSpeciesTo(iSpTransport) = &
                                                          & indicesTranspSpecies(iSpTransport)
        mappingNbr2Volume%refMainSpecies(iSpAerosol)%fract(iSpTransport) = &
                  & fu_conversion_factor(fu_basic_mass_unit(materialPtrTmp),'kg',materialPtrTmp) / &
                  & fu_dry_part_density(materialPtrTmp)

      end do  ! transport species 

    end do ! aerosol species

    mappingNbr2Volume%defined = silja_true

  end subroutine create_mass_2_nbr_mapping


  !****************************************************************************

  subroutine flip_moment_and_basic_var_Euler(moment_mapping, main_map, scale_map, conversion_direction)
    !
    ! This sub flips between the moment and basic variables as described in the given
    ! moment mapping. Note that for the opposite flipping, the other mapping is needed
    !
    implicit none

    ! Imported parameters
    type(Tmoment_mapping), intent(in) :: moment_mapping
    type(Tmass_map), pointer :: main_map
    type(Tmass_map), intent(in) :: scale_map
    integer, intent(in) :: conversion_direction

    ! Local variables
    integer :: iSpMain, iSpScale, ix, iy, iSrc, iLev, isp
    logical :: below_cutoff
    real :: fMTotal
    real, dimension(:), pointer :: cutoff_mass

    !if(.not. associated(moment_mapping))return
    if(.not. fu_true(moment_mapping%defined))return

    select case(moment_mapping%mapping_type)
      case(centre_mass_vs_moment_species, centre_mass_vs_moment_bulk)
        cutoff_mass => moment_mapping%cutoff_thresholds

        if (conversion_direction == to_moment) then

          do iy = 1, main_map%ny
            do ix = 1, main_map%nx
              do iLev = 1, main_map%n3d
                do iSrc = 1, main_map%nSrc
                  do iSpMain = 1, moment_mapping%nMainSpecies
                    fmTotal = 0.0
                    do iSpScale = 1, moment_mapping%refMainSpecies(iSpMain)%nRefSpecies
                      isp = moment_mapping%refMainSpecies(iSpMain)%indSpeciesTo(iSpScale)
                      fMTotal = fMTotal + scale_map%arM(isp, iSrc, iLev, ix, iy)
                    end do
                    main_map%arM(iSpMain,iSrc,iLev,ix,iy) = fMTotal * main_map%arM(iSpMain,iSrc,iLev,ix,iy)
                  end do
                end do
              end do
            end do
          end do

        else if (conversion_direction == from_moment) then
          do iy = 1, main_map%ny
            do ix = 1, main_map%nx
              do iLev = 1, main_map%n3d
                do iSrc = 1, main_map%nSrc
                  do iSpMain = 1, moment_mapping%nMainSpecies
                    below_cutoff = .true.
                    fmTotal = 0.0
                    do iSpScale = 1, moment_mapping%refMainSpecies(iSpMain)%nRefSpecies
                      isp = moment_mapping%refMainSpecies(iSpMain)%indSpeciesTo(iSpScale)
                      fMTotal = fMTotal + scale_map%arM(isp, iSrc, iLev, ix, iy)
                      below_cutoff = below_cutoff &
                                 & .and. scale_map%arM(isp, iSrc, iLev, ix, iy) < cutoff_mass(isp)
                    end do

                    if (below_cutoff) cycle

                    main_map%arM(iSpMain,iSrc,iLev,ix,iy) = main_map%arM(iSpMain,iSrc,iLev,ix,iy) / fMTotal
                    
                  end do
                end do
              end do
            end do
          end do

        else
          call set_error('strange conversion direction', &
                       & 'flip_moment_and_basic_var_Euler')
        end if

      case(volume_vs_number_cnc_flag)
        !
        ! Given: single-particle mass m=ro*d^3 / (6*Pi)
        ! Available: mass concentrations of all species M_i
        ! Obtain: number concentrations N
        ! ... or vice versa: the procedure will be the same: N=sum(M_i)/m or m=sum(M_i)/N
        !
        do iy = 1, main_map%ny
         do ix = 1, main_map%nx
          do iLev = 1, main_map%n3d
           do iSrc = 1, main_map%nSrc
            do iSpMain = 1, moment_mapping%nMainSpecies
             do iSpScale = 1, moment_mapping%refMainSpecies(iSpMain)%nRefSpecies
               fMTotal = fMTotal + &
                       & scale_map%arM( &
                                     & moment_mapping%refMainSpecies(iSpMain)%indSpeciesTo(iSpScale), &
                                     & iSrc, iLev, ix, iy)
             end do
             main_map%arM(iSpMain,iSrc,iLev,ix,iy) = fMTotal / &
                                         & main_map%arM(iSpMain,iSrc,iLev,ix,iy)
            end do
           end do
          end do
         end do
        end do
 
      case default
        call msg('Unknown moment mapping type:',moment_mapping%mapping_type)
        call set_error('Unknown moment mapping type','flip_moment_and_basic_var_Euler')
        return
    end select

  end subroutine flip_moment_and_basic_var_Euler


  !****************************************************************************

  subroutine flip_moment_and_basic_var_Lagr(moment_mapping, main_LPs, scale_LPs, nop, nSpecies, &
                                          & conversion_direction)
    !
    ! This sub flips between the moment and basic variables as described in the given
    ! moment mapping. Note that for the opposite flipping, the other mapping is needed
    !
    implicit none

    ! Imported parameters
    type(Tmoment_mapping), intent(in) :: moment_mapping
    real, dimension(:,:), pointer :: main_LPs
    real, dimension(:,:), intent(in) :: scale_LPs
    integer, intent(in) :: conversion_direction, nop, nSpecies

    ! Local variables
    integer :: iSpMain, iSpScale, iLP, isp
    logical :: below_cutoff
    real :: fMTotal
    real, dimension(:), pointer :: cutoff_mass

    if(.not. fu_true(moment_mapping%defined))return

    select case(moment_mapping%mapping_type)
      case(centre_mass_vs_moment_species, centre_mass_vs_moment_bulk)
        cutoff_mass => moment_mapping%cutoff_thresholds

        if (conversion_direction == to_moment) then

          do iLP = 1, nop
            do iSpMain = 1, moment_mapping%nMainSpecies
              fmTotal = 0.0
              do iSpScale = 1, moment_mapping%refMainSpecies(iSpMain)%nRefSpecies
                isp = moment_mapping%refMainSpecies(iSpMain)%indSpeciesTo(iSpScale)
                fMTotal = fMTotal + scale_LPs(isp, iLP)
              end do
              main_LPs(iSpMain,iLP) = fMTotal * main_LPs(iSpMain,iLP)
            end do
          end do

        else if (conversion_direction == from_moment) then
          do iLP = 1, nop
            do iSpMain = 1, moment_mapping%nMainSpecies
              below_cutoff = .true.
              fmTotal = 0.0
              do iSpScale = 1, moment_mapping%refMainSpecies(iSpMain)%nRefSpecies
                isp = moment_mapping%refMainSpecies(iSpMain)%indSpeciesTo(iSpScale)
                fMTotal = fMTotal + scale_LPs(isp, iLP)
                below_cutoff = below_cutoff .and. scale_LPs(isp, iLP) < cutoff_mass(isp)
              end do

              if (below_cutoff) cycle

              main_LPs(iSpMain,iLP) = main_LPs(iSpMain,iLP) / fMTotal
                    
            end do
          end do

        else
          call set_error('strange conversion direction', 'flip_moment_and_basic_var_Lagr')
        end if

      case(volume_vs_number_cnc_flag)
        !
        ! Given: single-particle mass m=ro*d^3 / (6*Pi)
        ! Available: mass concentrations of all species M_i
        ! Obtain: number concentrations N
        ! ... or vice versa: the procedure will be the same: N=sum(M_i)/m or m=sum(M_i)/N
        !
        do iLP = 1, nop
          do iSpMain = 1, moment_mapping%nMainSpecies
            do iSpScale = 1, moment_mapping%refMainSpecies(iSpMain)%nRefSpecies
               fMTotal = fMTotal + &
                       & scale_LPs(moment_mapping%refMainSpecies(iSpMain)%indSpeciesTo(iSpScale), iLP)
            end do
            main_LPs(iSpMain,iLP) = fMTotal / main_LPs(iSpMain,iLP)
          end do
        end do
 
      case default
        call msg('Unknown moment mapping type:',moment_mapping%mapping_type)
        call set_error('Unknown moment mapping type','flip_moment_and_basic_var_Lagr')
        return
    end select

  end subroutine flip_moment_and_basic_var_Lagr


  !******************************************************************************

  subroutine many_mass_maps_to_grads_file(pMapAr, iSrc, chOutFNm, now, factor)
    !
    ! Receives an array of mass maps and dumps them all into the given grads file
    !
    implicit none
    ! Improted parameters
    type(Tmass_map), dimension(:), intent(in) :: pMapAr
    integer, intent(in) :: iSrc
    character(len=*), intent(in) :: chOutFNm
    type(silja_time), intent(in) :: now
    real, intent(in) :: factor
    
    ! Local variables
    integer :: iMap, iFile, nMaps

    if(size(pMapAr) < 1 .or. size(pMapAr) > 10)then
      call set_error('Strange size of mass map array:'+ fu_str(size(pMapAr)),'many_mass_maps_to_grads_file')
      return
    endif
    !
    ! How many valid mass maps we have?
    !
    nMaps = size(pMapAr)
    do iMap = 1, size(pMapAr)
      if(.not.defined(pMapAr(iMap)))then
        nMaps = iMap -1
        exit
      endif
    end do
    !

    iFile = int_missing  ! the unit of the opened file to keep between writes

    do iMap = 1, nMaps
      call mass_map_to_grads_file(pMapAr(iMap), iSrc, chOutFNm, now, factor, indFile=iFile, ifClose= (iMap == nMaps))
    end do
    
  end subroutine many_mass_maps_to_grads_file


  !******************************************************************************

  subroutine mass_map_to_grads_file(pMap, iSrc, chOutFNm, now, factor, indFile, ifClose)
    !
    ! Sends a content of the mass map to the GrADS file, whose name is to be given
    !
    implicit none

    ! Improted parameters
    type(Tmass_map), intent(in) :: pMap
    integer, intent(in) :: iSrc
    character(len=*), intent(in) :: chOutFNm
    type(silja_time), intent(in) :: now
    real, intent(in) :: factor
    integer, intent(inout), optional :: indFile
    logical, intent(in), optional :: ifClose

    ! Local declarations:
    INTEGER :: igf, iLev, iSpecies, ix, iy, iSrcStart, iSrcEnd, iSrcTmp
    type(silja_field_id) :: idTmp
    real, dimension(:), pointer :: arTmp
    character(len=substNmLen)  :: cocktail_name

    ! Checking
    !
    if(iSrc == int_missing)then
      iSrcStart = 1
      iSrcEnd = pMap%nSrc
    else
      if(iSrc < 1 .or. iSrc > pMap%nSrc)then
        call msg('Strange source index:',iSrc)
        call set_error('Strange source index:','mass_map_to_gards_file')
        return
      endif
      iSrcStart = iSrc
      iSrcEnd = iSrc
    endif
    arTmp => fu_work_array()
    if(error)return

    igf = int_missing
    if(present(indFile)) igf = indFile
    !
    ! Ignore the file name if we have a ready-to-write structure
    !
    if(igf == int_missing)then
      igf = open_gradsfile_o('', &  ! directory not used 
                           & chOutfNm, &
                           & pMap%gridTemplate, ifMPIIO=smpi_use_mpiio_grads, &
                           & ifBuffered = smpi_use_mpiio_grads, &
                           & time_label_position = instant_fields)
      if(present(indFile)) indFile = igf
    endif

    if(igf <= 0)then
      call set_error('Failed map to grads writing: file index is wrong:' + fu_str(igf), &
                   & 'mass_map_to_gards_file')
      return
    endif
    !
    ! Go through the mass map and store field by field
    !
    do iSrcTmp = iSrcStart, iSrcEnd
      do iSpecies = 1, pMap%nSpecies
        do iLev = 1, pMap%n3D
          !
          ! Make the 1D vector of the data to store
          !
          do iy = 1, pMap%ny
            do ix = 1, pMap%nx
              arTmp(ix+(iy-1)*pMap%nx) = pMap%arM(iSpecies,iSrcTmp,iLev,ix,iy) * factor
            end do
          end do
          !
          ! Set the field id
          !
          if(pMap%quantity == emission_intensity_flag .or. pMap%quantity == reaction_rate_flag)then
            if(pMap%quantity == emission_intensity_flag) then
              cocktail_name = fu_substance_name(pMap%species(iSpecies))
            else
              write (unit=cocktail_name, fmt='(A,I3.3)') 'R',iSpecies
            endif

            idTmp = fu_set_field_id(met_src_missing,&
                                  & pMap%quantity, &
                                  & now, &             ! analysis_time
                                  & zero_interval, &   ! forecast_length
                                  & pMap%gridTemplate,&
                                  & fu_level(pMap%vertTemplate,iLev),&
                                  & zero_interval, &     ! & length_of_accumulation
                                  & zero_interval, &     ! & length_of_validity
                                  & forecast_flag, &
!                                  & species = pMap%species(iSpecies), &
                                  & chCocktail = cocktail_name)
          else
            idTmp = fu_set_field_id(met_src_missing,&
                                  & pMap%quantity, &
                                  & now, &             ! analysis_time
                                  & zero_interval, &   ! forecast_length
                                  & pMap%gridTemplate,&
                                  & fu_level(pMap%vertTemplate,iLev),&
                                  & zero_interval, &     ! & length_of_accumulation
                                  & zero_interval, &     ! & length_of_validity
                                  & forecast_flag, &
                                  & species = pMap%species(iSpecies))
          endif
          if(error)return
          !
          ! And write it down
          !
          CALL write_next_field_to_gradsfile(igf, idTmp, arTmp)
          if(error)return
        end do   ! iSpecies
      end do ! iLev
    end do !! iSrcTmp
    !
    ! Done the writing. Close the binary and write the ctl
    ! if this file was open here. In some cases we should close it anyway - if asked by ifClose
    !
    if(present(indFile))then
      if(present(ifClose))then
        if(ifClose) call close_gradsfile_o(igf,"")
      endif
    else
      call close_gradsfile_o(igf,"")
    endif
    
    call free_work_array(arTmp)

  end subroutine mass_map_to_grads_file


  !************************************************************************************
  
  subroutine mass_to_mixing_ratio(mass_map, met_buffer, horiz_interp_struct, if_interpolate)
    implicit none
    type(tmass_map), intent(inout) :: mass_map
    type(Tfield_buffer), intent(in) :: met_buffer
    type(THorizInterpStruct), pointer :: horiz_interp_struct
    logical, intent(in) :: if_interpolate

    real :: surf_pres, area, pres_bottom, pres_top, dx, dy, dp, column_mass, n_air
    integer :: ind_ps, ix, iy, iz
    type(silja_level) :: layer
    if (.not. fu_leveltype(mass_map%vertTemplate) == layer_btw_2_hybrid) then
      call set_error('Works only with hybrid vertical', 'set_const_mixing_ratio')
      return
    end if
    
    ind_ps = fu_index(met_buffer, surface_pressure_flag)
    if (ind_ps < 1) then
      call set_error('Cannot find surface pressure', 'cnc_to_mass_mixing_ratio')
      return
    end if
    
    do iy = 1, mass_map%ny
      do ix = 1, mass_map%nx
        surf_pres = fu_get_value(met_buffer%p2d(ind_ps), nx_meteo, ix, iy, &
                               & met_buffer%weight_past, horiz_interp_struct, if_interpolate)
        do iz = 1, mass_map%n3d
          layer = fu_level(mass_map%vertTemplate, iz)
          pres_bottom = fu_hybrid_level_pressure(fu_lower_boundary_of_layer(layer), surf_pres)
          pres_top = fu_hybrid_level_pressure(fu_upper_boundary_of_layer(layer), surf_pres)
          dx = fu_dx_cell_m(mass_map%gridTemplate, ix, iy)
          dy = fu_dy_cell_m(mass_map%gridTemplate, ix, iy)
          dp = pres_bottom - pres_top
          column_mass = dp / g
          n_air = dx*dy*column_mass / molecular_weight_air ! m.w air is in kg/mol
          mass_map%arm(:,:,iz,ix,iy) =  mass_map%arm(:,:,iz,ix,iy) / n_air
        end do
      end do
    end do

  end subroutine mass_to_mixing_ratio
  
    !************************************************************************************
  
  subroutine mixing_ratio_to_mass(mass_map, met_buffer, horiz_interp_struct, if_interpolate)
    implicit none
    type(tmass_map), intent(inout) :: mass_map
    type(Tfield_buffer), intent(in) :: met_buffer
    type(THorizInterpStruct), intent(in) :: horiz_interp_struct
    logical, intent(in) :: if_interpolate

    real :: surf_pres, area, pres_bottom, pres_top, dx, dy, dp, column_mass, n_air, &
         & weight_past, air_dens
    integer :: ind_ps, ix, iy, iz, ind_dens, i1d
    type(silja_level) :: layer

!!$    if (.not. fu_leveltype(mass_map%vertTemplate) == layer_btw_2_hybrid) then
!!$      call set_error('Works only with hybrid vertical', 'set_const_mixing_ratio')
!!$      return
!!$    end if


    weight_past = met_buffer%weight_past
    
    if (fu_leveltype(dispersion_vertical) == layer_btw_2_height) then
      ind_dens = fu_index(met_buffer, air_density_flag)
      if (ind_dens < 1) then
        call set_error('Cannot find air density', 'cnc_to_mass_mixing_ratio')
        return
      end if

      call vmr_to_mass_hgt_vert()
      return
    end if
    ind_ps = fu_index(met_buffer, surface_pressure_flag)
    if (ind_ps < 1) then
      call set_error('Cannot find surface pressure', 'cnc_to_mass_mixing_ratio')
      return
    end if
    
    do iy = 1, mass_map%ny
      do ix = 1, mass_map%nx
        i1d = (iy-1) * nx_dispersion + ix
        !surf_pres = fu_get_value(met_buffer%p2d(ind_ps), nx_meteo, ix, iy, &
        !                       & met_buffer%weight_past, horiz_interp_struct, if_interpolate)
        surf_pres = met_buffer%p2d(ind_ps)%past%ptr(i1d)*weight_past &
             & +   met_buffer%p2d(ind_ps)%future%ptr(i1d)*(1-weight_past)
          
        do iz = 1, mass_map%n3d
          !air_dens = met_buffer%p4d(ind_dens)%past%p2d(iz)%ptr(i1d)*weight_past &
          !     & +   met_buffer%p4d(ind_dens)%future%p2d(iz)%ptr(i1d)*(1-weight_past)
          layer = fu_level(mass_map%vertTemplate, iz)
          pres_bottom = fu_hybrid_level_pressure(fu_lower_boundary_of_layer(layer), surf_pres)
          pres_top = fu_hybrid_level_pressure(fu_upper_boundary_of_layer(layer), surf_pres)
          dx = fu_dx_cell_m(mass_map%gridTemplate, ix, iy)
          dy = fu_dy_cell_m(mass_map%gridTemplate, ix, iy)
          dp = pres_bottom - pres_top
          column_mass = dp / g
          n_air = dx*dy*column_mass / molecular_weight_air
          !m_air = dx*dy*column_mass / molecular_weight_air
          !n_air = air_dens * 1000.0 / molecular_weight_air
          mass_map%arm(:,:,iz,ix,iy) =  mass_map%arm(:,:,iz,ix,iy) * n_air
          !mass_map%arm(:,:,iz,ix,iy) =  mass_map%arm(:,:,iz,ix,iy) * n_air
        end do
      end do
    end do
    
  contains 
    
    subroutine vmr_to_mass_hgt_vert()
      implicit none
    
      real :: air_dens, dz

      do iy = 1, mass_map%ny
        do ix = 1, mass_map%nx
          i1d = (iy-1) * nx_dispersion + ix
          do iz = 1, mass_map%n3d
            air_dens = met_buffer%p4d(ind_dens)%past%p2d(iz)%ptr(i1d)*weight_past &
                 & +   met_buffer%p4d(ind_dens)%future%p2d(iz)%ptr(i1d)*(1-weight_past)
            layer = fu_level(mass_map%vertTemplate, iz)
            dz = fu_layer_thickness_m(layer)
            dx = fu_dx_cell_m(mass_map%gridTemplate, ix, iy)
            dy = fu_dy_cell_m(mass_map%gridTemplate, ix, iy)
            n_air = dx*dy*dz*air_dens / molecular_weight_air
            !m_air = dx*dy*column_mass / molecular_weight_air
            !n_air = air_dens * 1000.0 / molecular_weight_air
            mass_map%arm(:,:,iz,ix,iy) =  mass_map%arm(:,:,iz,ix,iy) * n_air
            !mass_map%arm(:,:,iz,ix,iy) =  mass_map%arm(:,:,iz,ix,iy) * n_air
          end do
        end do
      end do
      
    end subroutine vmr_to_mass_hgt_vert

  end subroutine mixing_ratio_to_mass

 
  !*****************************************************************

  integer function fu_nbr_of_species_of_massMap(massMap)
    !
    ! Returns the number of the species in the cocktail.
    ! A single species is a substance at some phase - gas phase or mode of some size
    !
    implicit none

    ! Imported parameters
    type(Tmass_map), intent(in) :: massMap

    fu_nbr_of_species_of_massMap = massMap%nSpecies

  end function fu_nbr_of_species_of_massMap


  !************************************************************************

  subroutine report_mass_map(massMap)
    implicit none
    type(Tmass_map), intent(in) :: massMap

    ! Local variables
    integer :: ix,iy,iz,iSrc,iSpecies 
    real, dimension(:), pointer :: fTot
    
    call msg('===================== REPORT of MASS MAP ==============================')
    if(fu_true(massMap%defined))then
      call msg('Quantity:' + fu_quantity_string(massMap%quantity))
      call msg('Grid:')
      call report(massMap%gridTemplate)
      call msg('vertical')
      call report(massMap%vertTemplate)
      call msg('Number of sources:',massMap%nSrc)
      fTot => fu_work_array()
      if(error)return
      fTot(1:massMap%nSpecies) = 0.0
      do iy = 1, massMap%ny
        do ix = 1, massMap%nx
          do iz = 1,massMap%n3d
            do iSrc = 1, massMap%nSrc
              do iSpecies = 1, massMap%nSpecies
                fTot(iSpecies) = fTot(iSpecies) + massMap%arM(iSpecies,iSrc,iz,ix,iy)
              end do  ! species
            end do ! iSrc
          end do ! iz
        end do  ! ix
      end do  ! iy
      do iSpecies = 1, massMap%nSpecies
        call msg('Total for:' + fu_str(massMap%species(iSpecies)), fTot(iSpecies))
      end do
      call free_work_array(fTot)
    else
      call msg('Undefined Mass Map')
    endif
    
    call msg('===================== END OF REPORT of MASS MAP =======================')

  end subroutine report_mass_map
  
  logical function fu_mass_map_defined(mass_map) result(is_defined)
    implicit none
    type(Tmass_map), intent(in) :: mass_map
    is_defined = fu_true(mass_map%defined)

  end function fu_mass_map_defined


  !***************************************************************

  !***********************************************************************************
  !
  !   Input data request and handling for aerosols
  !
  !***********************************************************************************

  subroutine aerosol_processes_input_needs(aerosolRules, meteo_input_local, q_dynamic, q_static)
    !
    ! Returns the required fields for the thermodiffusion computations
    !
    implicit none

    ! Imported parameters
    type(Tmeteo_input), intent(out) :: meteo_input_local
    integer, dimension(:), intent(out) :: q_dynamic, q_static
    type(Taerosol_rules), intent(in) :: aerosolRules

    ! Local variable
    integer :: iTmp

    meteo_input_local = meteo_input_empty
    q_dynamic = int_missing
    q_static = int_missing

  end subroutine aerosol_processes_input_needs


!  !**************************************************************************************
!  !**************************************************************************************
!  !
!  !   Thermodiffusion stuff. MOVED to ADVECTIONS
!  !
!  !**************************************************************************************
!  !**************************************************************************************
!
!
!  real function fu_thermodiff_velocity(vSettling, metdat) result(velocity)
!    !
!    ! Computes the thermodiffusion velocity from the given settling velocity and meteodata
!    !
!    implicit none
!
!    ! Imported parameters
!    real, intent(in) :: vSettling
!    type(Tmeteo_input), intent(in) :: metdat
!
!    !
!    ! Density-fluctuation correction term first: v=-Kz(iLev) * (nabla_T/T - nabla_P/P)
!    ! Only vertical direction is taken as the gradients along the horizontal directions are small
!    !
!    velocity = 0.0  !Kz * (grad_P / pressure - grad_T / temperature)
!
!    if(vSettling .eps. 0.0) return ! this is it for negligible settling
!
!    velocity = velocity - 0.66667 * (metdat%val(ind_pressure) / metdat%val(ind_dpressure_dz)) * &
!                        & log(3.* metdat%val(ind_Kz_1m) / &
!                            & fu_kinematic_viscosity(metdat%val(ind_tempr), &
!                                                   & metdat%val(ind_pressure))) * &
!                        & metdat%val(ind_dTheta_dz) / metdat%val(ind_pot_tempr) * vSettling
!
!  end function fu_thermodiff_velocity


  !****************************************************************************************
  !****************************************************************************************
  !
  ! Tmeteo_input stuff
  !
  !****************************************************************************************
  !****************************************************************************************

  !*************************************************************************************

  subroutine define_meteo_input(meteo_input, MeteoInputLocal, nRequests) !q_dynamic, q_static, q_dispersion)
    !
    ! Defines the meteo_input structure by merging the given arrays of dynamic, static and 
    ! dispersion quantities requested by ALL transformation and deposition modules for the run.
    ! Allocates the necessary memory for the structure.
    ! The sub is called only once, so the efficiency is not too important
    !
    implicit none

    ! Imported parameters
    type(Tmeteo_input), intent(out) :: meteo_input
    type(Tmeteo_input), dimension(:), intent(inout) :: MeteoInputLocal
    integer, intent(in) :: nRequests !! Length of MeteoInputLocal


    ! Local variables
    integer :: iReq, iQLocal, iTmp
    integer, dimension(max_quantities) :: iQTmp, iTypeTmp
    logical :: ifFound

    meteo_input%defined = silja_false

    !
    ! Scan the given local requests checking for repetitions, count
    ! the needed size of the structure
    !
    meteo_input%nQuantities = 0
    !
    ! First, dynamic quantities
    !
    do iReq = 1, nRequests

      if(fu_false( MeteoInputLocal(iReq)%defined))cycle

      do iQLocal = 1, MeteoInputLocal(iReq)%nQuantities

        do iTmp = 1, meteo_input%nQuantities
          if(meteo_input%quantity(iTmp) == MeteoInputLocal(iReq)%Quantity(iQLocal)) exit
        end do

        if (iTmp > meteo_input%nQuantities) then !! quantity not found
          !add to meteo_input
          meteo_input%nQuantities = iTmp 
          meteo_input%quantity(iTmp) = MeteoInputLocal(iReq)%Quantity(iQLocal)
          meteo_input%q_type(iTmp) = MeteoInputLocal(iReq)%q_type(iQLocal)
        endif
        MeteoInputLocal(iReq)%idx(iQLocal) = iTmp

      end do  ! quantities of meteo_input_local

    end do  ! iReq - requests

    !
    ! Having the stuff collected and counted, allocate the memory and fill-in the structure
    !
    if(meteo_input%nQuantities < 1)then
      call msg_warning('No quantities are needed for the transformations','define_meteo_input')
      meteo_input = meteo_input_empty
    else
      meteo_input%defined = silja_true
    endif

    call msg('Meteo input defined. nQuantities:',meteo_input%nQuantities)
    do iTmp = 1, meteo_input%nQuantities
      call msg('Quantity requested:' + fu_quantity_short_string(meteo_input%Quantity(iTmp)))
    end do

  end subroutine define_meteo_input


  !****************************************************************************************

  subroutine pick_meteo_pointers(meteo_input, meteo_buf_ptr, disp_buf_ptr)
    !
    ! Finds out the indices of the corresponding quantities in meteo and dispersion 
    ! buffers and stores them into the meteo_input structure. Called every time step
    ! after the set_buffer_pointers because that sub can change the order of quantities
    !
    implicit none

    ! Imported parameters
    type(Tmeteo_input), intent(inout) :: meteo_input
    type(Tfield_buffer), pointer :: meteo_buf_ptr, disp_buf_ptr

    ! Local variables
    integer :: iQ

    do iQ = 1, meteo_input%nQuantities
      !
      ! Depending on the quantity type, search for the corresponding buffer
      !
      select case(meteo_input%q_type(iQ))
        case(meteo_dynamic_flag, meteo_single_time_flag)
          meteo_input%idx(iQ) = fu_index(meteo_buf_ptr, meteo_input%quantity(iQ))
          
        case (dispersion_dynamic_flag, dispersion_single_time_flag)
          meteo_input%idx(iQ) = fu_index(disp_buf_ptr, meteo_input%quantity(iQ))

        case default
          call msg('Unknown quantity source type:', meteo_input%q_type(iQ))
          call set_error('Unknown quantity source type','pick_meteo_pointers')
      end select
      if(error .or. meteo_input%idx(iQ) == int_missing)then
        call set_error('Problem with quantity:'+fu_quantity_short_string(meteo_input%quantity(iQ)), &
                     & 'pick_meteo_pointers')
        return
      endif
    end do  ! cycle over the meteo_input quantities

  end subroutine pick_meteo_pointers


  !*******************************************************************************************

  subroutine fill_in_meteo_input_column(meteo_input, values, &
                               & meteo_buf, dispersion_buf, &
                               & ix, iy, nz_disp, & 
                               & interpCoefMet2DispHoriz, interpCoefMet2DispVert, &
                               & ifInterpMet2DispHoriz, ifInterpMet2DispVert)
    !
    ! Fills-in the meteo input values for the given grid cell in dispersion or meteorological grid
    ! Note that this function does not know any more about the meteorological buffers.
    ! It works directly with the fields and uses interpolation wherever needed
    ! Idea is to use it for both Lagrangian and Eulerian computations
    !
    implicit none

    ! Imported parameters
    type(Tmeteo_input), intent(in) :: meteo_input
    real, dimension(:,:), intent(out) :: values !!(values, nz)
    type(Tfield_buffer), pointer :: meteo_buf, dispersion_buf
    integer, intent(in) :: ix, iy, nz_disp
    type(THorizInterpStruct), intent(in) :: interpCoefMet2DispHoriz
    type(TVertInterpStruct), intent(in) :: interpCoefMet2DispVert
    logical, intent(in) :: ifInterpMet2DispHoriz, ifInterpMet2DispVert


    ! Local variables
    integer :: iQ, iLev, bufindex
    integer :: quantity
    logical :: if3d

    do iQ = 1, meteo_input%nQuantities
      quantity = meteo_input%quantity(iQ)
      if3d = fu_multi_level_quantity(quantity)
      select case(meteo_input%q_type(iQ))
        case(meteo_dynamic_flag, meteo_single_time_flag)
          do iLev = 1, nz_disp
               values(iQ, iLev) = fu_get_value(meteo_buf, &                 ! data buffer
                                  & meteo_input%idx(iQ), & ! field index in the buffer
                                  & nx_meteo, &                      ! nxFrom        
                                  & ix, iy, iLev, &                  ! 3D position in TO grid
                                  & interpCoefMet2DispHoriz, interpCoefMet2DispVert, & ! space interpolation
                                  & ifInterpMet2DispHoriz, ifInterpMet2DispVert) ! if interp needed
              if (.not. if3d) then
                values(iQ, 2:nz_disp) = values(iQ, 1)
                exit
              endif
          enddo
        case (dispersion_dynamic_flag, dispersion_single_time_flag)
          do iLev = 1, nz_disp
              values(iQ, iLev) = fu_get_value(dispersion_buf, &             ! data buffer
                                  & meteo_input%idx(iQ), &  ! field index in the buffer
                                  & nx_dispersion, &                  ! nxFrom
                                  & ix, iy, iLev, &                   ! 3D position in TO grid
                                  & interpCoefMet2DispHoriz, interpCoefMet2DispVert, & ! void here
                                  & .false., .false.)                 ! no interpolation needed
              if (.not. if3d) then
                values(iQ, 2:nz_disp) = values(iQ, 1)
                exit
              endif
          enddo
        case default
          call msg('Unknown quantity source type:', meteo_input%q_type(iQ))
          call set_error('Unknown quantity source type','fill_in_meteo_input_column')
      end select
      if(error)then
        call msg('Source type:', meteo_input%q_type(iQ))
        call msg('ix,iy:', ix, iy)
        call msg("Ilev", iLev)
        call set_error('Problem with quantity:'+fu_quantity_short_string(meteo_input%quantity(iQ)), &
                     & 'fill_in_meteo_input_column')
        return
      endif
      
    end do !iQ 

  end subroutine fill_in_meteo_input_column

  !*******************************************************************************

  subroutine fill_in_meteo_input(meteo_input, values, &
                               & meteo_buf_ptr, dispersion_buf_ptr, &
                               & ix, iy, iLev, & 
                               & interpCoefMet2DispHoriz, interpCoefMet2DispVert, &
                               & ifInterpMet2DispHoriz, ifInterpMet2DispVert)
    !
    ! Fills-in the meteo input values for the given grid cell in dispersion or meteorological grid
    ! Note that this function does not know any more about the meteorological buffers.
    ! It works directly with the fields and uses interpolation wherever needed
    ! Idea is to use it for both Lagrangian and Eulerian computations
    !
    implicit none

    ! Imported parameters
    type(Tmeteo_input), intent(in) :: meteo_input
    real, dimension(:), intent(inout) :: values
    type(Tfield_buffer), pointer :: meteo_buf_ptr, dispersion_buf_ptr
    integer, intent(in) :: ix, iy, iLev
    type(THorizInterpStruct), pointer :: interpCoefMet2DispHoriz
    type(TVertInterpStruct), pointer :: interpCoefMet2DispVert
    logical, intent(in) :: ifInterpMet2DispHoriz, ifInterpMet2DispVert

    ! Local variables
    integer :: iQ

    do iQ = 1, meteo_input%nQuantities
      
      select case(meteo_input%q_type(iQ))
        case(meteo_dynamic_flag, meteo_single_time_flag)
          values(iQ) = fu_get_value(meteo_buf_ptr, &                 ! data buffer
                                  & meteo_input%idx(iQ), & ! field index in the buffer
                                  & nx_meteo, &                      ! nxFrom        
                                  & ix, iy, iLev, &                  ! 3D position in TO grid
                                  & interpCoefMet2DispHoriz, interpCoefMet2DispVert, & ! space interpolation
                                  & ifInterpMet2DispHoriz, ifInterpMet2DispVert) ! if interp needed
        case (dispersion_dynamic_flag, dispersion_single_time_flag)
          values(iQ) = fu_get_value(dispersion_buf_ptr, &             ! data buffer
                                  & meteo_input%idx(iQ), &  ! field index in the buffer
                                  & nx_dispersion, &                  ! nxFrom
                                  & ix, iy, iLev, &                   ! 3D position in TO grid
                                  & interpCoefMet2DispHoriz, interpCoefMet2DispVert, & ! void here
                                  & .false., .false.)                 ! no interpolation needed
        case default
          call msg('Unknown quantity source type:', meteo_input%q_type(iQ))
          call set_error('Unknown quantity source type','fill_in_meteo_input')
      end select
      if(error)then
        call msg('Source type:', meteo_input%q_type(iQ))
        call msg('ix,iy:', ix, iy)
        call msg("Ilev", iLev)
        call set_error('Problem with quantity:'+fu_quantity_short_string(meteo_input%quantity(iQ)), &
                     & 'fill_in_meteo_input')
        return
      endif
    end do

  end subroutine fill_in_meteo_input


  !*******************************************************************************

  integer function fu_index_in_meteo_input(meteo_input, quantity)
    !
    ! Finds the given quantity and returns its index
    !
    implicit none

    type(Tmeteo_input), intent(in) :: meteo_input
    integer, intent(in) :: quantity

    ! Local variable
    integer :: iQ

    do iQ = 1, meteo_input%nQuantities
      if(meteo_input%quantity(iQ) == quantity)then
        fu_index_in_meteo_input = iQ
        return
      endif
    end do

    call set_error('Quantity:' + fu_quantity_short_string(quantity) + '- is not found in meteo input', &
                 & 'fu_index_in_meteo_input')

  end function fu_index_in_meteo_input


  !************************************************************************************

  real function fu_OH_cnc_forced(fSunZenithAngleCos, fHeight_m)
    !
    ! This function returns the OH value, which is parameterised from the sun zenith angle and 
    ! altitude.
    ! Hypothesis is: near surface we can parameterise the story using brute-force fitting
    ! of the observation data to some functional dependence. With altitude, we will take
    ! that OH concentrations are approaching those in the stratosphere by the tropopause.
    !
    ! Units: SI
    !
    implicit none

    ! Imported variables
    real, intent(in) :: fSunZenithAngleCos, fHeight_m

    ! Tuning:
    ! 29.10.2007: the function from EMEP is taken and the altitude dependence follows the
    !             O3 -> O1d, which roughly doubles by 25km
    ! 30.10.2007: 5-fold night OH, two-fold daytime OH following the Hyytiala observations
    ! 01.11.2007: return back to EMEP model, despite it seems to be somewhat too low. Reason: NO2
    !             degrades much too fast to HNO3. SO2 gas-phase conversion with OH is accompanied with
    !             a surrogate of aqueous-phase conversion.
    !
    ! EMEP:
    if(fSunZenithAngleCos <= 0.01)then  !Should not overflow single precision
      fu_OH_cnc_forced = 1.66e-14  ! 1e4 molecules/cm3 -> moles/m3
    else
      fu_OH_cnc_forced = (1.66e-14 + 6.64e-12 * exp(-0.25 / fSunZenithAngleCos)) * (1.+fHeight_m/25000.)
    endif
!    !
!    ! Adjustment 30.10.2007
!    if(fSunZenithAngleCos <= 0)then
!      fu_OH_cnc_forced = 8.3e-14  ! 5e4 molecules/cm3 -> moles/m3
!    else
!      fu_OH_cnc_forced = (1.66e-13 + 1.33e-11 * exp(-0.25 / fSunZenithAngleCos)) * (1.+fHeight_m/25000.)
!    endif
    
  end function fu_OH_cnc_forced

  !**************************************************************************************
    
  subroutine check_mass_moments_mass_map(MM_cnc, MM, chPlace, chMMName)
    !
    ! Actually checks the mass map
    !    
    implicit none
      
    type(TMass_map), intent(in) :: MM, MM_cnc
    character(len=*), intent(in) :: chPlace, chMMName
      
    integer :: ix, iy, iLev, iSrc, iSpecies, kTmp
    logical :: ifError
    real :: fMass, centre

    !
    ! Eulerian environment
    !
    call msg ("check_mass_moments_mass_map Place:"+chplace+", MMname:"+chMMName)


    !$OMP PARALLEL  default(none) &
    !$OMP & PRIVATE (ix, iy, iLev, iSrc, iSpecies, kTmp, ifError, fMass, centre) &
    !$OMP & SHARED (MM, MM_cnc,  chPlace, chMMName, error)
    ifError = .false.
    !$OMP DO collapse(4)
    do iy=1, MM%ny
      do ix=1, MM%nx
        do iLev = 1, MM%n3D  ! Vertical levels
          do iSrc = 1, MM%nSrc   ! emission sources
            if (error) cycle  ! No further checks
            do iSpecies = 1, MM%nSpecies
              fMass = MM_cnc%arM(iSpecies,iSrc,iLev,ix,iy)
              if (.not. (abs(MM%arM(iSpecies,iSrc,iLev,ix,iy)) <= 0.5*abs(fMass))) then
               !$OMP critical (bark_cm)
                centre = MM%arM(iSpecies,iSrc,iLev,ix,iy) / fMass
                call msg('Strange mass centre found,' + chPlace + ', in:' + chMMName + &
                       & ', ix,iy,Level, source, ind_species:', (/ix,iy, iLev,iSrc, iSpecies/))
                call msg('Mass centre, moment and mass:', &
                        &(/centre, MM%arM(iSpecies,iSrc,iLev,ix,iy), fMass/))
                !$OMP end critical (bark_cm)
                ifError = .true.
              endif
            end do  ! species
            !
            ! If wrong mass centre found, report the cell in more details
            !
            if(ifError)then
              !$OMP critical (bark_cm)
              do iSpecies = 1, MM%nSpecies
                do kTmp = 1, MM%n3D
                  call msg('Level, iSpecies, moment, mass:' + fu_str(kTmp) + ',' + fu_str(iSpecies), &
                            & MM%arM(iSpecies,iSrc,kTmp,ix,iy), MM_cnc%arM(iSpecies,iSrc,kTmp,ix,iy))
                end do
              end do
              call set_error('Wrong mass centre found at:' + chPlace + ', in:' + chMMName, &
                           & 'check_mass_moments_mass_map')
             !$OMP end critical (bark_cm)
            endif  ! error
          end do  ! iSrc
        end do  ! iLev
      end do  ! ix
    end do  ! iy
    !$OMP END DO 
    !$OMP END PARALLEL
  end subroutine check_mass_moments_mass_map

    
  !************************************************************************************
    
    subroutine check_mass_centres_mass_map( MM_cnc, MM, chPlace, chMMName, tolerance_factor)
      !
      ! Actually checks the mass map --- seems to change it as well...
      !    
      implicit none
      
      type(TMass_map), intent(inout) :: MM, MM_cnc
      character(len=*), intent(in) :: chPlace, chMMName
      real :: tolerance_factor
      
      integer :: ix, iy, iLev, iSrc, iSpecies, kTmp
      logical :: ifError
      real :: fTmp

      if (error) return
       call msg ("check_mass_centres_mass_map Place:"+chplace+", MMname:"+chMMName)
      !
      ! Eulerian environment
      !
      !$OMP PARALLEL  default(none) &
      !$OMP & PRIVATE (ix, iy, iLev, iSrc, iSpecies, kTmp, ifError, fTmp) &
      !$OMP & SHARED (MM, MM_cnc,  chPlace, chMMName, error, tolerance_factor)
      ifError = .false.
      !$OMP DO collapse(4)
      do iy=1, MM%ny
        do ix=1, MM%nx
         do iLev = 1, MM%n3D  ! Vertical levels
           do iSrc = 1, MM%nSrc   ! emission sources
            if (error) cycle
            do iSpecies = 1, MM%nSpecies
             if( .not. (abs(MM%arM(iSpecies,iSrc,iLev,ix,iy)) <=  0.5 ))then
               fTmp = MM%arM(iSpecies,iSrc,iLev,ix,iy) * tolerance_factor
               if( .not. (abs(fTmp) <=  0.5 ))then
                 !$OMP critical (bark_cm)
                 call msg('Strange mass centre found,' + chPlace + ', in:' + chMMName + &
                        & ', ix,iy,Level, source, ind_species:', (/ix,iy, iLev,iSrc, iSpecies/))
                 call msg('Mass centre and mass:', MM%arM(iSpecies,iSrc,iLev,ix,iy), &
                                                 & MM_cnc%arM(iSpecies,iSrc,iLev,ix,iy))
                 !$OMP end critical (bark_cm)
                 ifError = .true.
               else
                 MM%arM(iSpecies,iSrc,iLev,ix,iy) = fTmp
               endif
             endif
            end do !iSpecies
            !
            ! If wrong mass centre found, report the cell in more details
            !
            if(ifError)then
              !$OMP critical (bark_cm)
              do iSpecies = 1, MM%nSpecies
                do kTmp = 1, MM%n3D
                  call msg('Level =' + fu_str(kTmp) + ', iSpecies, mass centre:',iSpecies, &
                                                                & MM%arM(iSpecies,iSrc,kTmp,ix,iy))
                end do !iSpecies
              end do
              call set_error('Wrong mass centre found at:' + chPlace + ', in:' + chMMName, &
                           & 'check_mass_centres_mass_map')
              !$OMP end critical (bark_cm)
              cycle
            endif  ! error
           end do  !  iSrc
         end do  ! iLev 
       end do  ! ix
      end do  ! iy
      !$OMP END DO 
      !$OMP END PARALLEL
    end subroutine check_mass_centres_mass_map
  
END MODULE cocktail_basic


