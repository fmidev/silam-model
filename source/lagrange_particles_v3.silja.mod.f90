MODULE lagrange_particles 
  
  ! This module contains the basic type of Lagrange-type simulations
  ! A Lagrange particle is NOT a physical, solid particle, but rather a small volume
  ! of air with chemicals and physical particles inside. The Lagrange particle contains 
  ! info about its position, size, and chemical/physical composition. Additionally, it stores
  ! its own turbulent micro-movement because it depends not only on the
  ! meteorological parameters but also has own "memory" and inertion
  ! connected with the relaxation time of the turbulent motion - see
  ! random_walks_2 for details
  !
  ! Author: Mikhail Sofiev, FMI email Mikhail.Sofiev@fmi.fi
  !
  ! All units: SI
  ! 
  ! Language: ANSI standard Fortran 90
  ! 

  USE cocktail_basic ! field_buffer !chemistry_manager !server    ! The only allowed link to substance-related parts of the model

  IMPLICIT NONE

  ! The public functions and subroutines available in this module:

  PUBLIC init_lagrange_particles_set    ! Initialize all except for particle arrays
  PUBLIC enlarge_lagrange_particles_set    ! Extend particle arrays

  public project_lagr_2_euler_flds    ! Project some particles to Eulerian mass map
  public set_missing
  public report                       ! Report Lagrnagian particle
  
  ! Private subroutines  
  private set_lp_set_missing
  private set_lp_set_ptr_missing
  private report_lagr_particle

  INTERFACE report
    MODULE PROCEDURE report_lagr_particle
  END INTERFACE

  interface set_missing
    module procedure set_lp_set_missing
    module procedure set_lp_set_ptr_missing
  end interface

  !
  ! The concept of lagrangian particle in SILAM v.5
  ! Lagrangian particle is an object with the following properties:
  ! - position in space: x, y, z
  ! - turbulent part of the motion: uT, vT, wT
  ! - masses of the transported species: m_i, i=1,nSpeciesTransport
  ! - size of the puff represented by the particle: dx, dy, dz
  !
  ! That stands for 9+nSpeciesTransport reals made as three arrays in pollution_cloud:
  ! - arDyn(9, nParticles)
  ! - arMass(nSpeciesTransport, nParticles)
  ! - lpStatus (nParticles). Integer = age of LP,[MINUTES] * 100 + iSrc, or int_missing for empty LP
  ! These are addressed via presceribed indices:
  integer, parameter, private :: nDynamicVars = 9
  integer, parameter, public :: lp_x = 1
  integer, parameter, public :: lp_y = 2
  integer, parameter, public :: lp_z = 3
  integer, parameter, public :: lp_uT = 4
  integer, parameter, public :: lp_vT = 5
  integer, parameter, public :: lp_wT = 6
  integer, parameter, public :: lp_dx = 7
  integer, parameter, public :: lp_dy = 8
  integer, parameter, public :: lp_dz = 9
  !
  ! The main structure holding the Lagrangian environment. PUBLIC to speed-up the access
  !
  type Tlagrange_particles_set
    type(silam_species), dimension(:), pointer :: spTransp =>null(), spShortLived=>null(), spAeros=>null()
    INTEGER :: nop=0, nSpeciesTrn=0, nSpeciesSL=0, nSpeciesAer=0, nSrcs=0 ! number of particles, species, and sources
    integer :: iFirstEmptyParticle = 0, nop_active = 0 
    type(silja_grid) :: gridTemplate = grid_missing
    type(silam_vertical) :: verticalTemplate = vertical_missing
    real, dimension(:,:), allocatable :: lpMassTrn, lpMassSL, lpMassAer, lpDyn !(iSp,1:Nop) masses and motion params of LPs
    integer, dimension(:), allocatable :: lpStatus   ! (1:NoP)status of LPs: age in MINUTES * 100 + iSrc
    
    integer, dimension(:), allocatable :: newPartList ! (1:NoP) Indices if just emitted particles
    integer :: nNewPart = 0  ! number of just emitted particles, should be reset before emission
    type(silja_logical) :: defined = silja_false
  end type Tlagrange_particles_set
  public Tlagrange_particles_set

  type Tlagrange_particles_set_ptr
    type(Tlagrange_particles_set), pointer :: ptrLpSet => null()
  end type Tlagrange_particles_set_ptr


  CONTAINS

  !******************************************************************
  
  subroutine set_lp_set_missing(lp_set)
    implicit none
    type (Tlagrange_particles_set), intent(out) :: lp_set
    
    lp_set%spTransp => null()
    lp_set%spShortLived=>null()
    lp_set%spAeros=>null()
    lp_set%nop = 0
    lp_set%nSpeciesTrn = 0
    lp_set%nSpeciesSL = 0
    lp_set%nSpeciesAer = 0
    lp_set%nSrcs = 0 ! number of particles, species, and sources
    lp_set%iFirstEmptyParticle = 0
    lp_set%nop_active = 0 
    lp_set%gridTemplate = grid_missing
    lp_set%verticalTemplate = vertical_missing
    if(allocated(lp_set%lpMassTrn))deallocate(lp_set%lpMassTrn)
    if(allocated(lp_set%lpMassSL))deallocate(lp_set%lpMassSL)
    if(allocated(lp_set%lpMassAer))deallocate(lp_set%lpMassAer)
    if(allocated(lp_set%lpDyn))deallocate(lp_set%lpDyn)
    if(allocated(lp_set%lpStatus))deallocate(lp_set%lpStatus)
    if(allocated(lp_set%newPartList))deallocate(lp_set%newPartList)
    lp_set%nNewPart = 0  ! number of just emitted particles, should be reset before emission
    lp_set%defined = silja_false
  end subroutine set_lp_set_missing
  
  
  !******************************************************************
  
  subroutine set_lp_set_ptr_missing(lp_set_ptr)
    implicit none
    type (Tlagrange_particles_set_ptr), intent(out) :: lp_set_ptr
    
    if(associated(lp_set_ptr%ptrLpSet)) call set_lp_set_missing(lp_set_ptr%ptrLpSet)
  end subroutine set_lp_set_ptr_missing
    
  
  !******************************************************************

  subroutine init_lagrange_particles_set(lpset, spTransp, spShortLived, spAeros, grid, vertical, nSrc)
    
     implicit none
     type (Tlagrange_particles_set), intent(out) :: lpset
     type(silam_species), dimension(:), pointer :: spTransp, spShortLived, spAeros
     type(silam_vertical), intent(in) :: vertical
     type(silja_grid), intent(in) :: grid
     integer, intent(in) :: nSrc

     call  set_lp_set_missing(lpset)

     lpset%nop = 0
     lpset%nop_active = 0
     lpset%nNewPart = 0
     lpset%nSpeciesTrn = 0
     if (associated(spTransp)) then
        lpset%spTransp => spTransp
        lpset%nSpeciesTrn = size(spTransp)
     endif

     lpset%nSpeciesSL = 0
     if (associated(spShortLived)) then
        lpset%spShortLived => spShortLived
        lpset%nSpeciesSL = size(spShortLived)
     endif
   
     lpset%nSpeciesAer = 0
     if (associated(spAeros)) then
        lpset%spAeros => spAeros
        lpset%nSpeciesAer = size(spAeros)
     endif

     lpset%gridTemplate = grid
     lpset%verticalTemplate = vertical
     lpset%nSrcs = nSrc
     lpset%defined = silja_true

  end subroutine init_lagrange_particles_set


  !********************************************************************
  
  subroutine enlarge_lagrange_particles_set(lpset, nParticles) 
    !
    ! Reserves or enlarges arrays for largrangian particles.
    !
    implicit none
    
    ! Imported parameters
    type (Tlagrange_particles_set), intent(inout) :: lpset
    integer, intent(in) :: nParticles

    ! Local variables
    integer :: istat
    type (Tlagrange_particles_set) :: lpset_tmp !Just a container
    character (len=*), parameter :: sub_name="enlarge_lagrange_particles_set"

    if (.not. lpset%defined == silja_true) then
         call set_error("enlarge_lagrange_particles_set called on undefined lpset", sub_name)
         return
    endif

    if (nParticles == 0) then
       if (lpset%nop > 0) then
         call msg("Resetting Number of Particles in lpset to 0")
         lpset%nop = 0
         if (lpset%nSpeciesAer >0) deallocate(lpset_tmp%lpMassAer)
         if (lpset%nSpeciesSL >0)  deallocate(lpset_tmp%lpMassSL)
         deallocate(lpset_tmp%lpMassTrn,lpset_tmp%lpDyn,lpset_tmp%lpStatus,lpset_tmp%newPartList)
      else
         call msg("Resetting empty lpset: nothing to do")
      endif
      return
    endif

    if (nParticles <= lpset%nop) then
       call msg("Attempt to reduce number of particles in lpset old, new",lpset%nop, nParticles)
       call set_error("This is not implemented (yet?)", sub_name)
       return
    endif

    call set_lp_set_missing(lpset_tmp)
    
    if (lpset%nSpeciesAer >0) then 
       allocate(lpset_tmp%lpMassAer(lpset%nSpeciesAer,nParticles), stat=istat)
       if (iStat /=0) call set_error("Failed to allocate lpMassAer", sub_name)
       lpset_tmp%lpMassAer(:,1:lpset%nop) =  lpset%lpMassAer(:,1:lpset%nop)
       call move_alloc(lpset_tmp%lpMassAer, lpset%lpMassAer)
    endif

    if (lpset%nSpeciesSL >0) then 
       allocate(lpset_tmp%lpMassSL(lpset%nSpeciesSL,nParticles), stat=istat)
       if (iStat /=0) call set_error("Failed to allocate lpMassSL", sub_name)
       lpset_tmp%lpMassSL(:,1:lpset%nop) =  lpset%lpMassSL(:,1:lpset%nop)
       call move_alloc(lpset_tmp%lpMassSL, lpset%lpMassSL)
    endif
                
    allocate(lpset_tmp%lpMassTrn(lpset%nSpeciesTrn,nParticles), &
          &  lpset_tmp%lpDyn(nDynamicVars,nParticles), &
          &  lpset_tmp%lpStatus(nParticles),  &
          &  lpset_tmp%newPartList(nParticles) , stat=istat)
    if (iStat /=0) call set_error('Failed to allocate lpMass & co.', sub_name)
    
    if (lpset%nop > 0) then !Copy valuable stuff
       lpset_tmp%lpMassTrn(:,1:lpset%nop) = lpset%lpMassTrn(:,1:lpset%nop)
       lpset_tmp%lpDyn(:,1:lpset%nop)     = lpset%lpDyn(:,1:lpset%nop)
       lpset_tmp%lpStatus(1:lpset%nop)    = lpset%lpStatus(1:lpset%nop)
       lpset_tmp%newPartList(1:lpset%nNewPart) = lpset%newPartList(1:lpset%nNewPart) ! Only new ones
    else
       lpset%iFirstEmptyParticle = 1  !Init 
    endif

    call move_alloc(lpset_tmp%lpMassTrn,   lpset%lpMassTrn)
    call move_alloc(lpset_tmp%lpDyn,       lpset%lpDyn)
    call move_alloc(lpset_tmp%lpStatus,    lpset%lpStatus)
    call move_alloc(lpset_tmp%newPartList, lpset%newPartList)
    

    lpset%lpStatus( lpset%nop+1 : nParticles) = int_missing
    lpset%lpMassTrn(:, lpset%nop+1 : nParticles) = 0.0
    lpset%lpDyn (:, lpset%nop+1 : nParticles) = real_missing

    call msg("Enlarged iNoP in lpset from, to :", lpset%nop, nParticles)
    lpset%nop = nParticles


  end subroutine enlarge_lagrange_particles_set


  !************************************************************************************
  
  subroutine project_lagr_2_euler_flds(lpSet, &
                                     & mapConc, mapAerosol, &
                                     & mapPx_conc, mapPy_conc, mapPz_conc, &
                                     & projectionRules, residence_time_sec)
    !
    ! Projects the lagrangian particles to the Eulerian mass map.
    ! This is the cornerstone of the plume-in-grid technology: partciles that are
    ! of sufficient age are projected to Eulerian grid and eliminated from Lagrangian
    ! dynamics.
    ! There can be several criteria for the projection:
    ! - all particles can be forced to the Eulerian grid
    ! - particles older than smth
    ! - particles with horizontal size bigger than smth
    ! - ...
    !
    implicit none
    
    ! Imported parameters
    type(Tlagrange_particles_set), intent(inout) :: lpSet
    type(Tmass_map), intent(inout) :: mapConc, mapAerosol, &
                              & mapPx_conc, mapPy_conc, mapPz_conc           ! moments of centres of masses
    type(TLagr2Euler_projectionRules), intent(in) :: projectionRules
    real, intent(in) :: residence_time_sec

    ! Local variables
    integer :: iP, ix, iy, iz, iSrc, iSpecies
    real :: fX, fY, fZ, fTmp
    real, save :: dx_grid, dy_grid, dz_grid
    logical, save :: ifFirst = .true.
    !
    ! Just precaution: may be, nothing to project...
    !
    if(projectionRules%criterionType == never_project_flag)return

    call msg('Projecting Lagrnagian particles to Eulerian grid')
    if(ifFirst)then
      ifFirst = .false.
      dx_grid = fu_dx_cell_m(mapConc%gridTemplate, mapConc%nx/2, mapConc%ny/2)
      dy_grid = fu_dy_cell_m(mapConc%gridTemplate, mapConc%nx/2, mapConc%ny/2)
      dz_grid = fu_layer_thickness_m(fu_level(mapConc%vertTemplate,mapConc%n3D/2))
    endif
    !
    ! Cycle over the lagrangian particles. Each should be tested to the specified criterion
    ! and projected to Eulerian grid if needed
    ! Note that the particles fly in the Eulerian grid, so projection is easy.
    !
    do iP = 1, lpSet%nop

      ! Skip particle if not active 
      if(lpSet%lpStatus(iP) == int_missing)cycle  ! lpStatus = status + iSrcNbr*10

      ! ... or not satisfying the crietrion
      if(projectionRules%criterionType == maximum_age_flag)then
        if(int(lpSet%lpStatus(iP)/100) < projectionRules%LP_age_max / 60.)cycle  ! too young

      elseif(projectionRules%criterionType == residence_interval_fraction)then
        if(int(lpSet%lpStatus(iP)/100) < projectionRules%LP_age_max * residence_time_sec / 60. )cycle  ! too young

      elseif(projectionRules%criterionType == maximum_relative_size_flag)then
        if(lpSet%lpDyn(lp_dx,iP)/dx_grid + lpSet%lpDyn(lp_dy,iP)/dy_grid + &
         & lpSet%lpDyn(lp_dz,iP)/dz_grid < projectionRules%LP_relative_size_max) cycle  ! too small 3D perimeter

      else
        call set_error('Strange Lagr->Euler criterion:'+fu_str(projectionRules%criterionType), &
                     & 'project_lagr_2_euler_flds')
        return
      endif
      !
      ! Rest is to be projected and eliminated
      !
      fX = lpSet%lpDyn(lp_x,iP); ix = nint(fX)
      fY = lpSet%lpDyn(lp_y,iP); iy = nint(fY)
      fZ = lpSet%lpDyn(lp_z,iP); iz = nint(fZ)
      iSrc = mod(lpSet%lpStatus(iP),100)
      if(mapPx_conc%nSpecies == 1)then           ! For the case of bulk advection
        fTmp = 0.0
        do iSpecies = 1, mapConc%nSpecies
          fTmp = fTmp + lpSet%lpMassTrn(iSpecies,iP)
          mapConc%arM(iSpecies,iSrc,iz,ix,iy) = mapConc%arM(iSpecies,iSrc,iz,ix,iy) + fTmp
        end do
        mapPx_conc%arM(1,iSrc,iz,ix,iy) = mapPx_conc%arM(1,iSrc,iz,ix,iy) + fTmp * (fX-ix)
        mapPy_conc%arM(1,iSrc,iz,ix,iy) = mapPy_conc%arM(1,iSrc,iz,ix,iy) + fTmp * (fY-iy)
        mapPz_conc%arM(1,iSrc,iz,ix,iy) = mapPz_conc%arM(1,iSrc,iz,ix,iy) + fTmp * (fZ-iz)
      else
        do iSpecies = 1, mapConc%nSpecies        ! Species-moments advection
          fTmp = lpSet%lpMassTrn(iSpecies,iP)
          mapConc%arM(iSpecies,1,iz,ix,iy) = mapConc%arM(iSpecies,1,iz,ix,iy) + fTmp
          mapPx_conc%arM(iSpecies,1,iz,ix,iy) = mapPx_conc%arM(iSpecies,1,iz,ix,iy) + fTmp * (fX-ix)
          mapPy_conc%arM(iSpecies,1,iz,ix,iy) = mapPy_conc%arM(iSpecies,1,iz,ix,iy) + fTmp * (fY-iy)
          mapPz_conc%arM(iSpecies,1,iz,ix,iy) = mapPz_conc%arM(iSpecies,1,iz,ix,iy) + fTmp * (fZ-iz)
        end do
      endif  ! bulk or species momentum
      
      lpSet%lpStatus(iP) = int_missing  ! ready for the next activation
      if(iP < lpSet%iFirstEmptyParticle) lpSet%iFirstEmptyParticle = iP
    end do ! iP

  end subroutine project_lagr_2_euler_flds


  !********************************************************************
  
  SUBROUTINE report_lagr_particle(lpset, iParticle_, ifDump_)
    
    ! Prints the contents of a particle for test purposes. 
    ! 
    IMPLICIT NONE
    !
    ! Imported parameters with intent IN:
    type(Tlagrange_particles_set), intent(inout) :: lpSet
    integer, intent(in), optional :: iParticle_
    logical, intent(in), optional :: ifDump_

    integer :: iP, nFree, nInAir, nSrc, iSrc
    logical :: ifDump
    real :: fAgeHr

    !
    ! Stupidity check
    !
    if( lpset%defined == silja_false)then
      call msg('Lpset undefined')
      return
    endif

    if( lpset%NoP < 1)then
      call msg('Lpset Empty')
      return
    endif
    !
    ! If one particle, report and return
    !
    if(present(iParticle_))then
      call report_single_particle(iParticle_)
      return
    endif
    if(present(ifDump_))then
      ifDump = ifDump_
    else
      ifDump = .false.
    endif
    
    nFree = 0
    nInAir = 0
    nSrc = 0
    fAgeHr = 0.0

    do iP = 1, size(lpset%lpStatus)
      if(lpset%lpStatus(iP) == int_missing)then
        nFree = nFree + 1
      else
        nInAir = nInAir + 1
        iSrc = mod(lpset%lpStatus(iP),100)
        nSrc = max(nSrc,iSrc)
        fAgeHr = fAgeHr + real(lpset%lpStatus(iP)) / 6000.
      endif
      if(ifDump) call report_single_particle(iP)
    end do

    call msg('Particles free and in-air:', nFree, nInAir)
    call msg('Max source nbr and mean age [hr]:', nSrc, fAgeHr / (real(nInAir)+0.001))

    CONTAINS
    
    subroutine report_single_particle(iPart)
      implicit none
      integer, intent(in) :: iPart
      integer :: iTmp
      if(lpset%lpStatus(iPart) == int_missing)then
        call msg('Free:', iPart)
      else
        call msg('In-air:', iPart)
        do iTmp = 1, size(lpset%lpMassTrn,1)
          call msg('Mass of species:', iTmp, lpset%lpMassTrn(iTmp, iPart))
        end do
        call msg('X coord and size:', lpset%lpDyn(iTmp, lp_x), lpset%lpDyn(iTmp,lp_dx))
        call msg('Y coord and size:', lpset%lpDyn(iTmp, lp_y), lpset%lpDyn(iTmp,lp_dy))
        call msg('Z coord and size:', lpset%lpDyn(iTmp, lp_z), lpset%lpDyn(iTmp,lp_dz))
        call msg('X and Y turb velocity:', lpset%lpDyn(iTmp, lp_uT), lpset%lpDyn(iTmp,lp_vT))
        call msg('Z turb velocity:', lpset%lpDyn(iTmp, lp_wT))
        call msg('Source index and age [hr]',mod(lpset%lpStatus(iPart),100), real(lpset%lpStatus(iPart)) / 6000.)
      endif
    end subroutine report_single_particle
  END SUBROUTINE report_lagr_particle

END MODULE lagrange_particles
