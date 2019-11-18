MODULE optical_density
  !
  ! Contains the stuff necessary for for computing optical depth and density.
  ! 
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! The actual computation is done so far in a very stupid and inefficient way - 
  ! maps are computed and stored for all the substances for all wavelenghts and 
  ! only the needed ones are put to output.
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !
  ! Language: FORTRAN-90
  !
  ! Units: NOT NECESSARILY SI. 
  ! Warning is issued in each non-SI case
  !
  ! Author: Marje Prank marje.prank@fmi.fi
  !
  !use cocktail_basic
  use dispersion_server

  implicit none

  !
  ! Public routines
  !
  public set_optical_density_rules
  public defined
  public set_missing

  public fu_if_species_optics_known
  public add_optical_dens_input_needs
  public optical_column_input_needs

  !! Stuff for AOD-dependent photochemistry
  public get_photoatt_aod
  public init_photoatt_lut


  !!! Stuff needed for aod output
  public init_optical_density_data
  public make_optical_dens_map
  public make_optical_column_depth_map
  public compute_optical_dens_coef
  public sum_opt_dns_to_col_depth
  !
  ! Private routines
  !
  private defined_optical_rules
  private set_missing_optical_rules


  private MIECALC
  private Mie_Approx


  interface defined
    module procedure defined_optical_rules
  end interface
  interface set_missing
    module procedure set_missing_optical_rules
  end interface

  ! Type for the internal structure for the optical density computations
  ! Dimensions: - species(substance, mode, wavelength) 
  !             - humidity or temperature
  ! Phase transition for aerosols taken 20% deliquesence RH & 80% cristallisatioin RH

  !
  integer, parameter :: extinction_flag = 80001, scatter_flag = 80002, BackScatter_flag=80003

  integer, private, pointer, save :: ind_tempr => null() ! Used for addressing the meteo data in meteo input
  integer, private, pointer, save :: ind_rh    => null() ! Used for addressing the meteo data in meteo input
  integer, private, pointer, save :: ind_dx_size => null()
  integer, private, pointer, save :: ind_dy_size => null()

  type Textinction_coef
    private
    real, dimension(:,:), allocatable :: coef   ! (nSpecies, nMetOPTlut=21) !!m2/MassUnit for extinction_flag, scatter_flag
                                            ! m2/MassUnit/sr for backscatter
                                            ! Second dimension for dependence on meteo-data. 
    logical, dimension(:), allocatable :: ifRhDepSp, ifTDepSp ! (nSpecies) redundant but unified
    integer :: nSpecies
    integer :: coeff_type
    logical :: defined
  end type Textinction_coef

  public Textinction_coef

  ! Relative humidity for aerosols and temperature for gases,
  ! both divided to 21 intervals.
  ! RH: 0-50, 55, 60, 65, 70, ..., 120
  ! T: 200, 205, 210, ..., 300
  integer, parameter :: nMetOPTlut=21 ! 200, 205, 210, ..., 300K or 0, 50, 53, 56, .., 145%
  real, dimension(nMetOPTlut), parameter :: RHinter = &
      & (/0., 0.50, 0.54, 0.58, 0.62, 0.66, 0.70, 0.73, 0.76, 0.79, 0.82, 0.85, 0.87, 0.89, 0.91, 0.93, 0.95, 0.96, 0.97, 0.98, 0.99/)
  real, dimension(nMetOPTlut), parameter :: Tinter = &
      & (/200., 205., 210., 215., 220., 225., 230., 235., 240., 245., 250., 255., 260., 265., 270., 275., 280., 285., 290., 295., 300./)

  ! .. and the structure itself
  !
  type(Textinction_coef), private, target, save :: Extinction_coef  !!!Structure for optical_density output

  type(Textinction_coef), private, target, save :: Photo_ext, Photo_scat !!! Lookup tables for photolysis attenuation

  !-------------------------------------------------------------
  !
  ! The main set of rules for optical density computations
  !
  type Toptical_density_rules
    private
    logical :: ifRelHumidDep
    logical :: ifTDep
    type(silja_logical) :: defined
  end type Toptical_density_rules

  type (Toptical_density_rules), parameter :: optical_density_rules_missing=Toptical_density_rules(.false.,.false.,silja_false)

  type Toptical_properties
    private
    character(len=clen) :: chOpticStdSubst
    integer :: OpticalType
    real:: min_wl, max_wl
    integer :: nWaves, nRh, nT
    real, dimension(:), allocatable :: waves
    real, dimension(:), allocatable :: RH  !! Allocated only if nRh>1 and OpticalType==1
    real, dimension(:), allocatable :: T  !! Allocated only if nRh>1 and OpticalType==2
    complex, dimension(:,:), allocatable :: refrInd  !!!if OptType=1  !! Refractive index only
    !! (1:nRHs,1:nWavs) 
    real, dimension(:,:), allocatable :: xsect  !! m2/mole !!!if OptType=2    !! xsection only
                                                 !!!(1:nT,nWavs)
    type(silja_logical) :: defined
  end type Toptical_properties

  type (Toptical_properties), dimension(:), allocatable, target, private :: refProperties


CONTAINS

  !*******************************************************************
  !*******************************************************************
  !
  ! Encapsulation for the optical rules
  !
  !*******************************************************************
  !*******************************************************************

  subroutine set_optical_density_rules(rulesOptDns, nlSetup)
    implicit none
    type(Toptical_density_rules), intent(inout) :: rulesOptDns
    type(Tsilam_namelist), pointer :: nlSetup
    rulesOptDns%ifRelHumidDep = &
      & .not.(fu_str_u_case(fu_content(nlSetup, &
                                & 'optical_coefficients_depend_on_relative_humidity')) == 'NO')
    rulesOptDns%ifTDep = &
      & .not.(fu_str_u_case(fu_content(nlSetup,'optical_coefficients_depend_on_temperature')) == 'NO')
    rulesOptDns%defined = silja_true
  end subroutine set_optical_density_rules

  !*********************************************************************
  subroutine set_missing_optical_rules(rules)
    type(Toptical_density_rules), intent(out) :: rules
    rules%defined = silja_false
  end subroutine set_missing_optical_rules

  !*********************************************************************
  logical function defined_optical_rules(rules)
    type(Toptical_density_rules), intent(in) :: rules
    defined_optical_rules = (rules%defined == silja_true)
  end function defined_optical_rules


  !*********************************************************************
  
  logical function fu_if_species_optics_known(species, wavelength)
    !
    ! Just checks whether the optical properties of the given species are in the database
    !
    implicit none
    
    ! Imported parameters
    type(silam_species), intent(in) :: species
    real, intent(in) :: wavelength !! wavelength to use if species has no wavelength

    ! Local variables
    character(len=26), parameter :: subname = 'fu_if_species_optics_known'
    integer :: iTmp, iStatus, nItems, iRefSubst, iVal
    character(len=clen) :: spTmp
    type(Toptical_features), pointer :: pOpticFeatures
    real :: fTmp
    type (Toptical_properties), pointer :: optPrp

    fu_if_species_optics_known = .false.

    !
    ! First of all, does the substance given in the species have the optical parameters?
    !
    nullify(pOpticFeatures)
    pOpticFeatures => fu_optical_features(fu_material(species))
    if(.not. associated(pOpticFeatures))return
    if(pOpticFeatures%nOptSubst < 1 .or. pOpticFeatures%nOptSubst > 100)return
    if(.not. associated(pOpticFeatures%chOpticStdSubst))return
    if(size(pOpticFeatures%chOpticStdSubst) < pOpticFeatures%nOptSubst)return
    if(.not. associated(pOpticFeatures%fractionOptStdSubst))return
    if(size(pOpticFeatures%fractionOptStdSubst) < pOpticFeatures%nOptSubst)return

    if (.not. allocated(refProperties)) then
      call set_error("refProperties not initialized (yet?)", subname)
      return
    endif
    
    do iRefSubst = 1, pOpticFeatures%nOptSubst
      
      ! 
      !Find ref subst
      spTmp = pOpticFeatures%chOpticStdSubst(iRefSubst)
      do iTmp = 1, size(refProperties)
        optPrp => refProperties(iTmp)
        if (trim(optPrp%chOpticStdSubst) == trim(spTmp)) exit
      enddo
      if (iTmp >  size(refProperties)) then
         call set_error("Reference subst '"//trim(spTmp)//"' is missing from metadata!", subname)
         return
      endif

      !Check right mode
      if (fu_mode(species) == in_gas_phase .neqv. (optPrp%OpticalType == 2))then
         call msg("")
         call report(species)
         call msg_warning("Wrong OpticalType for ref subst '"//trim(spTmp)//"'", subname)
         exit
      endif

      !check wavelength range
      fTmp = fu_optical_wave_length(species)

      if (fTmp == real_missing) fTmp = wavelength

      if (.not. ( fTmp >= optPrp%min_wl .and. fTmp <= optPrp%max_wl) ) then
        call msg("")
        call report(species)
        call msg("Wavelength for ref substance '"//trim(spTmp)//"', nm", &
                & nint(optPrp%min_wl*1e9), nint(optPrp%max_wl*1e9))
        call msg_warning("Species wavelength outside the range of optical properties", subname)
        exit
      endif
    enddo

    !Last cycle completed
    if (iRefSubst > pOpticFeatures%nOptSubst) fu_if_species_optics_known = .true.

  end function fu_if_species_optics_known


  !*******************************************************************
  !*******************************************************************
  !
  ! Mapping of substances and waves
  !
  !*******************************************************************
  !*******************************************************************


  !****************************************************************

  subroutine optical_column_input_needs(opticalRules, meteo_input_local)
    !
    ! Fills meteo input for in-transformation 
    !
    implicit none

    ! Imported parameters
    type (Toptical_density_rules), intent(in) :: opticalRules
    type(Tmeteo_input), intent(out), target :: meteo_input_local
    character(len=*), parameter :: subname="optical_column_input_needs"

    ! Local variables
    integer :: iQ, iTmp, nq
    
    meteo_input_local =  meteo_input_empty
    if(.not. opticalRules%defined == silja_true)then
      call set_error('Optical density rules are not defined',subname)
      return
    endif

    nq=2
    meteo_input_local%quantity(1:2) = (/cell_size_x_flag, cell_size_y_flag/)
    meteo_input_local%q_type(1:2) = dispersion_single_time_flag
    ind_dx_size =>  meteo_input_local%idx(1)
    ind_dy_size =>  meteo_input_local%idx(2)


    if(opticalRules%ifRelHumidDep) then
      nq = nq+1
      meteo_input_local%quantity(nq) = relative_humidity_flag
      meteo_input_local%q_type(nq) = meteo_dynamic_flag
      ind_rH => meteo_input_local%idx(nq)
    endif

    if(opticalRules%ifTDep) then
      nq = nq+1
      meteo_input_local%quantity(nq) = temperature_flag    ! For self-degradation
      meteo_input_local%q_type(nq) = meteo_dynamic_flag
      ind_tempr => meteo_input_local%idx(nq)
    endif
    meteo_input_local%nQuantities = nq

  end subroutine optical_column_input_needs

  !***********************************************************************

  subroutine add_optical_dens_input_needs(iOutQuantities, q_met_dyn, q_met_st, &
                               & q_disp_dyn, q_disp_st, opticalRules)
    !
    ! Analyses the required output, finds out the diagnosed variables and 
    ! forms a request for the needed prognostic vars (and other input, i.e. meteo)
    !
    implicit none

    ! Imported parameters
    integer, dimension(:), intent(in) :: iOutQuantities
    integer, dimension(:), intent(inout) :: q_met_dyn, q_met_st, &
                               & q_disp_dyn, q_disp_st
    type (Toptical_density_rules), intent(in) :: opticalRules

    ! Local variables
    integer :: iQ, iTmp

    !
    ! Scan the array of quantities requested for the output, locate the diagnostic
    ! variables and fills-in the needed quantities - separately dynamic and static
    !
    do iQ = 1, size(iOutQuantities)
      if(iOutQuantities(iQ) == int_missing) return   ! all done
      !
      ! Check the quantity
      !
      if((iOutQuantities(iQ) == optical_density_flag) .or. &
     &   (iOutQuantities(iQ) == optical_column_depth_flag))then
        
        if(opticalRules%defined == silja_true)then
          if(opticalRules%ifRelHumidDep) &
                               & iTmp = fu_merge_integer_to_array(relative_humidity_flag, q_met_dyn)
          if(opticalRules%ifTDep) &
                               & iTmp = fu_merge_integer_to_array(temperature_flag, q_met_dyn)
        else
          call set_error('Optical density rules are not defined','add_optical_dens_input_needs')
          return
        endif
        exit
      end if  
    end do ! over the list of output quantities

  end subroutine add_optical_dens_input_needs


  !***********************************************************************
  
  subroutine read_optical_properties(fname, refProp)
    !
    ! Parse optical_properties.dat file and fill refProperties array
    !

    character(len=*), intent(in) :: fname
    type (Toptical_properties), dimension(:), allocatable, intent(inout) :: refProp

    !Locals
    integer :: iStat, iTmp, iRh, iWave, iRefSub, iUnit, iStatus, iTempr
    integer :: nRefSub, nRHs, nTemprs, nVal, nWavs
    integer :: optDataType
    type(Tsilam_namelist_group), pointer :: nlGrpPtr
    type(Tsilam_namelist), pointer :: nlPtr
    type(Tsilam_nl_item_ptr), dimension(:), pointer :: pItems
    character(len=fnlen) :: strTmp
    character(len=clen) :: chrTmp1,chrTmp2
    real :: fTmp1, fTmp2
    character(len=*), parameter :: subname="read_optical_properties"

    iUnit = fu_next_free_unit()

    open(file=fname, unit=iUnit, action='read', status='old', iostat=iStatus)
    
    if(iStatus /= 0)then
      call set_error('Failed to open optical properties list: '// trim(strTmp), subname)
      return
    endif

    nlPtr => fu_read_namelist(iUnit, .false., 'END_LIST = metadata')
    nullify(pItems)
    call get_items(nlPtr, 'reference_substance', pItems, nRefSub)
    if(error .or. nRefSub < 1)then
      call set_error('Failed to get the reference_substance items from the metadata namelist', subname)
    endif

    allocate(refProp(nRefSub), stat=iStat)

    do iRefSub = 1, nRefSub
        refProp(iRefSub)%defined = silja_false
        strTmp = fu_content(pItems(iRefSub))
! # <reference_subst_name>  <data_type>  <MIn_Wavelength_Known> <MIWK_unit>  <MAx_Wavelen_Known> <MAWK_unit>
!LIST = metadata
! reference_substance = seasalt    1    250 nm  40 mkm
        !!read(unit=strTmp, fmt=*, stat=istat) refProp(iRefSub)%chOpticStdSubst, refProp(iRefSub)%OpticalType, &
        read (unit=strTmp, fmt=*, iostat=iStat) refProp(iRefSub)%chOpticStdSubst, refProp(iRefSub)%OpticalType, &
                  & fTmp1, chrTmp1, fTmp2, chrTmp2
        refProp(iRefSub)%min_wl =  fTmp1 * fu_factor_to_basic_unit(chrTmp1, 'm')
        refProp(iRefSub)%max_wl =  fTmp2 * fu_factor_to_basic_unit(chrTmp2, 'm')
        if (error .or. iStat /=0 ) then
          call msg("Failes to parse line")
          call msg(strTmp)
          call set_error("Error parsing metadata", subname)
          return
        endif
    enddo
    call destroy_namelist(nlPtr)

    nlGrpPtr => fu_read_namelist_group(iUnit, .false.,'END_OPTICAL_PROPERTIES')
    close(iUnit)

    do iRefSub = 1, nRefSub
        nlPtr => fu_namelist(nlGrpPtr, refProp(iRefSub)%chOpticStdSubst)  
        !call msg("RefProp"+refProp(iRefSub)%chOpticStdSubst) 

        if(.not. associated(nlPtr))then
          call set_error('Optical properties namelist not associated for '// &
              & trim(refProp(iRefSub)%chOpticStdSubst), subname)
          return
        endif  

        optDataType = fu_content_int(nlPtr, 'dataType')  ! 1 - refractive indexes for aerosols (need Mie computations)
                                                         ! 2 - cross-sections for gases
        if (optDataType /= refProp(iRefSub)%OpticalType) then
          call set_error("OptType mismatch for "//trim(refProp(iRefSub)%chOpticStdSubst), subname)
          return
        endif
        ! Wavelengths
        strTmp = fu_content(nlPtr, 'Wavelength_unit')
        call get_items(nlPtr, 'Lambda', pItems, nWavs)
        if(nWavs < 1)then
           call set_error('No wavelengths in namelist',subname)
           return
        endif
        refProp(iRefSub)%nWaves = nWavs
        allocate(refProp(iRefSub)%waves(nWavs))
        
        fTmp1 = fu_factor_to_basic_unit(strTmp, 'm')
        do iTmp = 1, nWavs
          refProp(iRefSub)%waves(iTmp) = fTmp1 * fu_content_real(pItems(iTmp))
        enddo
        if(error)return

 
        if(optDataType == 1)then  !refractive indices for Mie computations
          if(fu_content(nlPtr,'ifRelativeHumidityDependent')== 'YES')then
            call get_items(nlPtr, 'RH', pItems, nRHs)
            if(fu_fails(nRHs >= 1,'No relative humidities in namelist',subname))return
            allocate(refProp(iRefSub)%RH(nRHs))

            do iTmp = 1, nRHs
            ! Relative humidities in database in %-s, model uses fractions
              refProp(iRefSub)%RH(iTmp) = 0.01*fu_content_real(pItems(iTmp))
            enddo
          else !if not RH dependent
            nRHs = 1
            allocate(refProp(iRefSub)%RH(nRHs))
            refProp(iRefSub)%RH(1) = 0.5
          endif
          refProp(iRefSub)%nRH=nRHs

          call get_items(nlPtr, 'val', pItems, nVal)

          if(fu_fails(nVal == nWavs * nRHs,'Wrong number of refractive indices in namelist', &
                    & subname))return
          allocate(refProp(iRefSub)%refrInd(nRHs,nWavs))

          do iRH = 1, nRHs
            do iWave = 1, nWavs
              iTmp = iWave + (iRH-1)*nWavs
              strTmp = fu_content(pItems(iTmp))
              read(unit=strTmp, iostat=iStatus, fmt=*) fTmp1, fTmp2
              if (iStatus /= 0) then
                call set_error("Failed to parse line '"//trim(strTmp)//"'",subname)
                return
              endif
              refProp(iRefSub)%refrInd(iRH, iWave) = cmplx(fTmp1, fTmp2)
            !  call msg("Table RI", (/refProp(iRefSub)%RH(iRH), refProp(iRefSub)%waves(iWave),fTmp1, fTmp2/))
            enddo
          enddo
        elseif(optDataType == 2)then !cross sections for molecules

          if(fu_content(nlPtr,'ifTemperatureDependent')== 'YES')then

            call get_items(nlPtr, 'Temperature', pItems, nTemprs)
            if(fu_fails(nTemprs >= 1,'No temperatures in namelist',subname))return

            allocate(refProp(iRefSub)%T(nTemprs))
            do iTmp = 1, nTemprs
              refProp(iRefSub)%T(iTmp) = fu_content_real(pItems(iTmp))
            enddo
          else 
            nTemprs = 1
            allocate(refProp(iRefSub)%T(1))
            refProp(iRefSub)%T(1) = 273.
          endif
          refProp(iRefSub)%nT=nTemprs

          !Cross-sections
          !
          strTmp = fu_content(nlPtr, 'Csect_unit')
          fTmp1 =  fu_factor_to_basic_unit(strTmp, 'm')
          
          call get_items(nlPtr, 'val', pItems, nVal)
          if(nVal /= nWavs * nTemprs) then
            call set_error('Wrong number of cross sections in namelist', subname)
            return
          endif
          allocate(refProp(iRefSub)%xsect(nTemprs,nWavs), stat = iTmp)
          if(fu_fails(iTmp == 0,'Failed to allocate temporary CS array',subname))return
          do iTempr = 1, nTemprs
            do iWave = 1, nWavs
              iTmp =  iWave + (iTempr-1)*nWavs 
              refProp(iRefSub)%xsect(iTempr,iWave) = fTmp1 * fu_content_real(pItems(iTmp))
            enddo
          enddo

        else
          call msg("Wrong namelist")
          call  report(nlPtr)
          call set_error("unknown optDataType",subname)
        endif

        refProp(iRefSub)%defined = silja_true
        IF (error) RETURN
    end do  ! reference optic substances

  endsubroutine read_optical_properties



  subroutine init_optical_LUT(XS_lut, XS_type, speciesLst, nspeciesLst, ref_wavelength)  
    !
    ! Prepares the extinction coefficients for the optical depth computations
    ! Receives the cocktail from where the needed species are to be taken, as well as the
    ! aerosol modes
    !
    ! ATTENTION. Units are NOT NECESSARILY SI
    ! However, input data come in SI.
    !
    implicit none

    ! Imported parameters
    type(Textinction_coef), intent(out) :: XS_lut
    integer, intent(in) :: XS_type ! one of extinction_flag, scatter_flag, BackScatter_flag
    type(silam_species), dimension(:), intent(in) :: speciesLst
    integer, intent(in) :: nspeciesLst
    real, intent(in) :: ref_wavelength ! one for soecies that do not have a wavelength


    ! Local variables

    integer :: iSpecies, iMode, iWave, iStatus, iTempr, iRH, iFract, iMet, iRefSubst, iTmp
    real :: fTmp, fTmp1, fTmp2, rad_grth, pDens, wavelength, wavelength_species, wLow
    real(r8k) :: minRad, maxRad, meanRad, dRad, X, qExt, qscat, qback, rad
    complex*16, dimension(nMetOPTlut) :: RImet
    real(r8k), dimension(nMetOPTlut) :: XSmet
    integer :: iRad, nTerms 
    integer, parameter :: ndRads=200  !number of radius intervals for extinction coefficient computation
    type(Toptical_features), pointer :: optFeats

    character(len=fnlen) :: strTmp
    type(TwetParticle) :: wetParticle
    type (Toptical_properties), pointer :: optPrp
    character(len=*), parameter :: subname = "init_optical_LUT"
!#######

    
   XS_lut%defined = .False.

   allocate(XS_lut%coef(nspeciesLst, nMetOPTlut),  &
      & XS_lut%ifRhDepSp(nspeciesLst), XS_lut%ifTDepSp(nspeciesLst), stat = iStatus)

    if(iStatus /= 0)then
      call set_error('Failed to allocate the extinction coefficient structure', subname)
      return
    endif
    XS_lut%coef(:, :) = 0 !(iSp, iMet)
    XS_lut%ifRHDepSp(:) = .false.
    XS_lut%ifTDepSp(:) = .false.

#ifdef DEBUG
      call msg("Making Crossection LUT "//subname)
#endif

    ! Cycle to read the optical data and compute the extinction coefficients
    !
species:  do iSpecies = 1, nspeciesLst
        
      optFeats => Fu_optical_features(fu_material(speciesLst(iSpecies)))
       
      wavelength_species = fu_optical_wave_length(speciesLst(iSpecies))

      if(.not. associated(optFeats))then
#ifdef DEBUG
        call msg('Optical features are not associated for '//fu_str(speciesLst(iSpecies)))
#endif
        !!! If species does not have a wavelength assocoated -- ignore it
        if (wavelength_species == real_missing) cycle  species
        call report(fu_material(speciesLst(iSpecies)))
        call set_error('Optical features are not associated',subname)
        return
      endif 

      wavelength = wavelength_species
      if (wavelength == real_missing) wavelength = ref_wavelength
      if (wavelength == real_missing) then
        call msg("")
        call report(speciesLst(iSpecies))
        call set_error("missing wavelength for species", subname)
        return
      endif

      if(.not. associated(optFeats%chOpticStdSubst))then
        call msg('Optical substance is not associated for ')
        call report(fu_material(speciesLst(iSpecies)))
        call set_error('Optical substance is not associated',subname)
        return
      endif  

      do iFract = 1, optFeats%nOptSubst

        !!!Find ref subst
        do iRefSubst= 1, size(refProperties)
         if (refProperties(iRefSubst)%chOpticStdSubst == optFeats%chOpticStdSubst(iFract) ) exit
        enddo

        if (iRefSubst > size(refProperties)) then 
          call msg("Could not find reference substance: "//trim(optFeats%chOpticStdSubst(iFract)))
          call msg("available substances are:")
          do iRefSubst= 1, size(refProperties)
             call msg(refProperties(iRefSubst)%chOpticStdSubst)
          enddo
          call set_error("Could not find reference substance", subname)
          return
        endif
        optPrp => refProperties(iRefSubst)

        !Find needed wavelength index
        do iWave=1,optPrp%nWaves-1
           if (wavelength <=  optPrp%waves(iWave+1)) then
             if (wavelength >=  optPrp%waves(iWave)) exit
           endif
        enddo
!        call msg("Wave no", iwave, nint( wavelength*1e9))
!        call msg("interpolated between", nint(optPrp%waves(iWave)*1e9), nint(optPrp%waves(iWave+1)*1e9) )
        if (iWave >= optPrp%nWaves) then
          call msg("Wavelength index failed for", nint(wavelength*1e9))
          call msg("RefSubst: "//trim(optPrp%chOpticStdSubst))
          if (wavelength_species == real_missing) cycle species
          call msg("Wl available (nm):", nint(optPrp%waves*1e9))
          call set_error("Failed to find wavelrngth index", subname)
          return
        endif
        wLow = (optPrp%waves(iWave+1)-wavelength)/(optPrp%waves(iWave+1) -  optPrp%waves(iWave))


!######Refractive indices for aerosols ####################
 
        if(optPrp%OpticalType == 1)then  !refractive indices for Mie computations

          fTmp1 = fu_deliquescence_humidity(fu_material(speciesLst(iSpecies)))

          if (optPrp%nRH == 1) then
               RImet(:) = optPrp%refrInd(1,iWave)*wLow + optPrp%refrInd(1,iWave+1)*(1.-wLow)
          else
            if (optPrp%nRh > 1)  XS_lut%ifRhDepSp(iSpecies) = .true.
            iRH = 1
            do iMet = 1,  nMetOPTlut
              do iRh=iRh, optPrp%nRH-1
                if (RHinter(iMet) < optPrp%RH(iRh+1)) exit
              enddo
              if (iRH < optPrp%nRH .and. RHinter(iMet) >= optPrp%RH(iRH)) then  !!!Interpolatable
                if (RHinter(iMet) <= fTmp1) then !! at/below deliquiscence pont
                    RImet(iMet) = optPrp%refrInd(1,iWave)*wLow + optPrp%refrInd(1,iWave+1)*(1.-wLow)
                  !  call msg("MetInterd (/iMet, iRh, iWave/)", (/iMet, 1, iWave/))
                else
                   fTmp2 = (RHinter(iMet) - optPrp%RH(iRH+1))/(optPrp%RH(iRH) - optPrp%RH(iRH+1))
                   RImet(iMet) = (optPrp%refrInd(iRh,  iWave)*wLow + optPrp%refrInd(iRh,  iWave+1)*(1.-wLow)) * fTmp2 + &
                   &             (optPrp%refrInd(iRh+1,iWave)*wLow + optPrp%refrInd(iRh+1,iWave+1)*(1.-wLow)) * (1.-fTmp2)
                  ! call msg("MetInterh (/iMet, iRh, iWave/)", (/iMet, iRh, iWave/))
                endif
              else !!Not interpolatable
                call msg("Could not interpolate to RH for "//trim(optPrp%chOpticStdSubst),RHinter(iMet))
                call msg("optPrp%RH(1:nRH)", optPrp%RH(1:optPrp%nRh))
                call set_error("Failed to interpolate RH", subname)
                return
              endif
            end do  ! iMet
          endif

          !!!Some species like AVB bins cen be gases, but still have optPrp%OpticalType ==1
          !! Ignore them
          !!! 
          if (fu_mode_type(speciesLst(iSpecies)) == gas_phase_flag) then
              XS_lut%coef(iSpecies, :) = 0.

          else 

            ! Conversion from mass extinction to molar extinction, if necessary 
            !
            fTmp = fu_conversion_factor(fu_basic_mass_unit(fu_material(speciesLst(iSpecies))), &
                                      & 'kg', fu_material(speciesLst(iSpecies)))

            pDens = fu_dry_part_density(fu_material(speciesLst(iSpecies)))
            if(error)return

            ! Cycle over relative humidities in the table
            !

            do iMet = 1, nMetOPTlut
              wetParticle = fu_wet_particle_features(fu_material(speciesLst(iSpecies)), RHinter(iMet))
              rad_grth = wetParticle%fGrowthFactor/2    ! includes rad = d/2
              minRad = fu_min_d(speciesLst(iSpecies)) * rad_grth !wetParticle%fGrowthFactor/2 ! RHgrowth(iMet)
              maxRad = fu_max_d(speciesLst(iSpecies)) * rad_grth !wetParticle%fGrowthFactor/2 ! RHgrowth(iMet)

              dRad = (maxRad-minRad)/ndRads

              fTmp1 = 0 !!Accumulator for this fraction over diameters
              !!call msg("for MIECALC",(/real(iMet), real(Real(RImet(iMet))), real(aimag(RImet(iMet)))  /) )
              do iRad = 1, ndRads + 1        !Integration loop over radius of spheres
                 rad = minRad + (iRad-1) * dRad
                 X = 2. * Pi * rad / wavelength
                 nTerms = X + 4.0 * X**0.3334 + 2.

                 call MIECALC (X,  RImet(iMet), nTerms, qExt, qscat, qback)
                !!call msg("X,  nTerms, qExt, qscat, qback",(/real(X), real(nTerms), real(qExt), real(qscat), real(qback)/))

                !Particle extinction efficiency is multiplied by relative humidity dependent particle
                !cross section and divided by dry particle mass to get mass extinction efficiencies
                 if (XS_type == extinction_flag) then
                   fTmp1 = fTmp1 + (6. * qExt* rad_grth**3 / (ndRads * rad * pDens)) 
                 elseif (XS_type == scatter_flag) then
                   fTmp1 = fTmp1 + (6. * qScat* rad_grth**3 / (ndRads * rad * pDens)) 
                 else
                   fTmp1 = fTmp1 + (6. * qBack* rad_grth**3 / (ndRads * rad * pDens)) 
                 endif
              enddo

              
              fTmp1 = fTmp1 * optFeats%fractionOptStdSubst(iFract) * fTmp



              XS_lut%coef(iSpecies,iMet) = XS_lut%coef(iSpecies, iMet) + fTmp1 
              if (iMet == 1) then
                 if ((abs(wetParticle%fGrowthFactor - 1.0) < 1e-5) .and. (optPrp%nRH == 1)) then
                    XS_lut%coef(iSpecies, 2:) = XS_lut%coef(iSpecies, 2:) + fTmp1
                    exit !!iMet loop
                 endif
              endif
            enddo  !!iMet = 1, nMetOPTlut
          endif !!Not gas

!##### Cross-section data for gases ##########

        elseif(optPrp%OpticalType == 2)then !cross sections

          if (optPrp%nT == 1) then
               XSmet(:) = optPrp%xsect(1,iWave)*wLow + optPrp%xsect(1,iWave+1)*(1.-wLow)
          else
            if (optPrp%nT > 1)  XS_lut%ifTDepSp(iSpecies) = .true.
            iTempr = 1
            do iMet = 1,  nMetOPTlut
              !below lower
              if (Tinter(iMet) <= optPrp%T(1)) then
                XSmet(iMet) = optPrp%xsect(1,iWave)*wLow + optPrp%xsect(1,iWave+1)*(1.-wLow)
                cycle
              endif

              !above upper
              if (Tinter(iMet) >= optPrp%T(optPrp%nT)) then
                XSmet(iMet) = optPrp%xsect(optPrp%nT,iWave)*wLow + optPrp%xsect(optPrp%nT,iWave+1)*(1.-wLow)
                cycle
              endif

              do iTempr=iTempr, optPrp%nT-1
                if (Tinter(iMet) < optPrp%T(iTempr+1)) exit
              enddo

              if (iTempr < optPrp%nT .and. Tinter(iMet) >= optPrp%T(iTempr)) then  !!!Interpolatable

                   fTmp2 = (Tinter(iMet) - optPrp%T(iTempr+1))/(optPrp%T(iTempr) - optPrp%T(iTempr+1))
                   XSmet(iMet) = (optPrp%xsect(iTempr,iWave)*wLow + optPrp%xsect(iTempr,iWave+1)*(1.-wLow)) * fTmp2 + &
                     & (optPrp%xsect(iTempr+1,iWave)*wLow + optPrp%xsect(iTempr+1,iWave+1)*(1.-wLow)) * (1.-fTmp2)
              else !!Not interpolatable (Should we extrapolate????)
                call msg("Could not interpolate to T for "//trim(optPrp%chOpticStdSubst),Tinter(iMet))
                call msg("optPrp%T(1:nT)", optPrp%T)
                call set_error("Failed to interpolate T", subname)
                return
              endif
            end do  ! iMet
          endif
          fTmp = avogadro * optFeats%fractionOptStdSubst(iFract)
          XS_lut%coef(iSpecies,:) = XS_lut%coef(iSpecies,:) + XSmet(:) *fTmp

          if (XS_type == backscatter_flag .and. iFract == optFeats%nOptSubst  ) then
            XS_lut%coef(iSpecies,:) = XS_lut%coef(iSpecies,:) * (3./(8.*PI))
          endif

        else
          call set_error('Strange optical data type',subname)
          return
        endif

      end do !optical standard substance
    end do species!subst
   XS_lut%defined = .True.


  end subroutine init_optical_LUT
 

  !************************************************************************
  subroutine init_optical_metadata(nlStdSetup)
    type(Tsilam_namelist), intent(in) :: nlStdSetup

    character(len=fnlen) :: fname
    character(len=*), parameter :: subname="init_optical_metadata"

    if (allocated(refProperties)) then
      call msg_warning("refProperties already initialized", subname)
    else
      fname = fu_content(nlStdSetup,'optical_properties_meta_data_file')
      call read_optical_properties(fname, refProperties)
      if (error) then
        call set_error("after read_optical_properties", subname)
        return
      endif
    endif 
  end subroutine init_optical_metadata


!!!  subroutine optical_column()

  subroutine init_photoatt_lut(speciesTransport, nSpeciesTransport, wavelength)
    !
    ! Prepares the extinction coefficients for the photolysis attenuation
    ! aerosol modes
    !
    !Wrapper tht calls init_optical_LUT for Photo_ext and Photo_scat
    ! that is used bu make_optical_dens_map etc..
    ! Should be removed at some point ()
    !
    implicit none

    type(silam_species), dimension(:), intent(in) :: speciesTransport
    integer, intent(in) :: nSpeciesTransport
    real, intent(in) :: wavelength
    integer :: iTmp


    ! Local variables
    character(len=*), parameter :: subname="init_photoatt_lut"

    if (.not. allocated(refProperties)) then
        call set_error("Optical reference properties not initialized", subname)
        return
    endif 


    call init_optical_LUT(Photo_ext, extinction_flag , speciesTransport, nSpeciesTransport, wavelength)
    if (error) call set_error("after init_optical_LUT Photo_ext", subname)

    call init_optical_LUT(Photo_scat, extinction_flag , speciesTransport, nSpeciesTransport, wavelength)
    if (error) call set_error("after init_optical_LUT Photo_scat,", subname)

    do iTmp=1,nSpeciesTransport
      call report(speciesTransport(iTmp)) 
    enddo
    do iTmp=1,nMetOPTlut
      call msg("Photo_ext, LUT: ", Photo_ext%coef(1:nSpeciesTransport, iTmp))
    enddo
    do iTmp=1,nMetOPTlut
      call msg("Photo_scat LUT: ", Photo_scat%coef(1:nSpeciesTransport, iTmp))
    enddo

  end subroutine init_photoatt_lut

  subroutine get_photoatt_aod(metdat_col,masses, ext_col, scat_col)
    real, dimension(:,:), intent(in)   :: metdat_col !!metind, iLev
    real, dimension(:,:,:), intent(in) :: masses     !nsrc, nsp, nLev  slice of a massmap
    real, dimension(:), intent(out) :: ext_col, scat_col

    integer :: iLev, nLev, imetRh, iMetT, iSp, nSp

    real ::  Rh, T, cellarea, extTmp, scaTmp
    logical :: ifRh, ifT
   ! logical, dimension(max_species) :: mask

    nLev = size(metdat_col,2)
    nSp = size(masses,1)
   !    mask(1:nSp) = .not.(Photo_ext%ifTDepSp(1:nSp))
    imetRh = 1
 !   imetT = 1
    ifRh = associated(ind_rh)
!    ifT = associated(ind_tempr)
     cellarea = metdat_col(ind_dx_size,1) * metdat_col(ind_dy_size,1)

    do iLev = 1,nLev
       if (ifRh) then
         Rh = metdat_col(ind_rh,iLev)
         imetRh = iMet4Rh(Rh)
       endif
!       if (ifT) then
!         T =  metdat_col(ind_tempr,iLev)
!         imetT = iMet4T(T)
!       endif
       !!Sum all that is not gas   
       ext_col(iLev)  = sum(masses(1:nSp,1,iLev)* Photo_ext%coef(1:nSp,imetRh)) / cellarea 
       scat_col(iLev) = sum(masses(1:nSp,1,iLev)*Photo_scat%coef(1:nSp,imetRh)) / cellarea
    enddo

        
  end subroutine get_photoatt_aod



  !***********************************************************************
  !***********************************************************************
  !***********************************************************************
  !***********************************************************************
  !***********************************************************************
  !***********************************************************************
  !***********************************************************************

  subroutine init_optical_density_data(speciesOptical, nSpeciesOptical, nlStdSetup, ChemRunSetup, &
                                     & opticalRules)
    !
    ! Prepares the extinction coefficients for the optical depth computations
    ! Receives the cocktail from where the needed species are to be taken, as well as the
    ! aerosol modes
    !
    ! For now just wrapper that calls init_optical_LUT for in_module Extinction_coef
    ! that is used bu make_optical_dens_map etc..
    ! Should be removed at some point ()
    !
    implicit none

    type(silam_species), dimension(:), intent(in) :: speciesOptical
    integer, intent(in) :: nSpeciesOptical
    type(Tsilam_namelist), pointer :: nlStdSetup
    type(TChemicalRunSetup), pointer :: ChemRunSetup
    type (Toptical_density_rules), intent(inout) :: opticalRules
    integer :: iTmp


    ! Local variables
    character(len=*), parameter :: subname="init_optical_density_data"

    if (.not. allocated(refProperties)) then
        call set_error("Optical reference properties not initialized", subname)
        return
    endif 


    call init_optical_LUT(Extinction_coef, extinction_flag , speciesOptical, nSpeciesOptical, real_missing)
    if (error) call set_error("after init_optical_LUT", subname)

    do iTmp=1,nSpeciesOptical
      call report(speciesOptical(iTmp)) 
    enddo
    do iTmp=1,nMetOPTlut
      call msg("Extinction_coef LUT: ", Extinction_coef%coef(1:nSpeciesOptical, iTmp))
    enddo

  end subroutine init_optical_density_data
 

  !************************************************************************************

  subroutine make_optical_dens_map(pDispFlds, met_buf, &
                                 & interpCoefMet2DispHoriz, interpCoefMet2DispVert, &
                                 & ifInterpMet2DispHoriz, ifInterpMet2DispVert, &
                                 & pOptDns, rulesOptDns, ChemRunSetup)
    !
    ! Actually computes the optical density map from concentration map and the 
    ! pre-computed extinction coefficient matrix.
    ! The concentrations are taken from the pDispFlds map, optical density
    ! stored to the pOptDepth map.

    implicit none

    ! Imported parameters
    type(Tmass_map), intent(in) :: pDispFlds
    type(Tmass_map), intent(inout) :: pOptDns  
    type(THorizInterpStruct), intent(in) :: interpCoefMet2DispHoriz
    type(TVertInterpStruct), intent(in) :: interpCoefMet2DispVert
    logical, intent(in) :: ifInterpMet2DispHoriz, ifInterpMet2DispVert
    TYPE(Tfield_buffer), intent(in), target :: met_buf
    type(TChemicalRunSetup), intent(in) :: ChemRunSetup

    integer, dimension(:), pointer :: mdl_in_q
    integer :: ix, iy, iLev, iSrc, iWave, iQ, iSpeciesOpt, iSpeciesTr, iSpTo
    integer :: nSpeciesOpt
    type(field_4d_data_ptr), pointer :: ptrRelHumid, ptrTemp
    type(Toptical_density_rules), intent(in) :: rulesOptDns
    real, dimension(max_species) :: arrOptDens
    real :: relHumid, tmpr

     mdl_in_q => met_buf%buffer_quantities
     if(rulesOptDns%ifRelHumidDep)then 
       iQ = fu_index(mdl_in_q, relative_humidity_flag)
       if(iQ <= 0)then
         call set_error('No relative humidity','make_optical_column_dens_map')
         return
       endif
       ptrRelHumid => met_buf%p4d(iQ)
     endif
     if (rulesOptDns%ifTDep)then
      iQ = fu_index(mdl_in_q, temperature_flag)
       if(iQ <= 0)then
         call set_error('No temperature','make_optical_column_depth_map')
         return
       endif
       ptrTemp => met_buf%p4d(iQ)
     endif
     !
     ! Nullify the entire optical map
     !
     nSpeciesOpt = pOptDns%nSpecies
     pOptDns%arM(1:pOptDns%nSpecies, &
                 & 1:pOptDns%nSrc, &
                 & 1:pOptDns%n3D, &
                 & 1:pOptDns%nx, &
                 & 1:pOptDns%ny) = 0.0
     !
     ! Now fill it from the transport cocktail
     !
     do iy = 1, pOptDns%ny
      do ix = 1, pOptDns%nx
       do iLev = 1, pOptDns%n3D
        do iSrc  = 1, pOptDns%nSrc

            arrOptDens(1:nSpeciesOpt) = 0

            if(rulesOptDns%ifRelHumidDep)then
                 relHumid = fu_get_value(ptrRelHumid, nx_meteo, ix, iy, iLev, met_buf%weight_past, &  
                                    & interpCoefMet2DispHoriz, interpCoefMet2DispVert, &
                                    & ifInterpMet2DispHoriz, ifInterpMet2DispVert)
            endif
            if (rulesOptDns%ifTDep)then 
                 tmpr = fu_get_value(ptrTemp, nx_meteo, ix, iy, iLev, met_buf%weight_past, &  
                                    & interpCoefMet2DispHoriz, interpCoefMet2DispVert, &
                                    & ifInterpMet2DispHoriz, ifInterpMet2DispVert)
            endif
            do iSpeciesTr = 1, pDispFlds%nSpecies  ! scan the transport cocktail
              !
              ! Contribute to the right optical species pointed by the ChemRunSetup
              !
              do iSpTo = 1, ChemRunSetup%refsTransp2opt(iSpeciesTr)%nRefSpecies
                      
                iSpeciesOpt = ChemRunSetup%refsTransp2opt(iSpeciesTr)%indSpeciesTo(iSpTo)


                arrOptDens(iSpeciesOpt) = arrOptDens(iSpeciesOpt) + &
                                       &  pDispFlds%arM(iSpeciesTr,iSrc,iLev,ix,iy) * &
                                       & fu_Extinction_coef(iSpeciesOpt, relHumid, tmpr )

              end do  ! iSpTo for the given iSpeciesTr
            end do  ! iSpeciesTr - of transport cocktail
          do iSpeciesOpt = 1, pOptDns%nSpecies
            pOptDns%arM(iSpeciesOpt,iSrc,iLev,ix,iy) = arrOptDens(iSpeciesOpt)
          enddo
        end do ! iSrc
       end do ! iLev
      end do ! ix
     end do ! iy

  end subroutine make_optical_dens_map


!*******************************************************************************


  subroutine make_optical_column_depth_map(pDispFlds, met_buf, &
                                         & interpCoefMet2DispHoriz, interpCoefMet2DispVert, &
                                         & ifInterpMet2DispHoriz, ifInterpMet2DispVert, &
                                         & pOpticColDepth, rulesOptDns, ChemRunSetup, pOpticDensity)
    !
    ! Actually computes the optical column depth map from concentration map and the 
    ! pre-computed extinction coefficient matrix.
    ! The concentrations are taken from the pDispFlds map, optical density
    ! stored to the pOptDepth map.

    implicit none

    ! Imported parameters
    type(Tmass_map), intent(in) :: pDispFlds, pOpticColDepth  
    type(THorizInterpStruct), intent(in) :: interpCoefMet2DispHoriz
    type(TVertInterpStruct), intent(in) :: interpCoefMet2DispVert
    logical, intent(in) :: ifInterpMet2DispHoriz, ifInterpMet2DispVert
    TYPE(Tfield_buffer), intent(in), target :: met_buf
    type(TChemicalRunSetup), intent(in) :: ChemRunSetup
    type(Tmass_map), intent(in) :: pOpticDensity !! Can be mass_map_missing

    integer, dimension(:), pointer :: mdl_in_q
    integer :: ix, iy, iSrc, iWave, iQ, nSpeciesOpt, iLev, iSpeciesOpt, iSpeciesTr, iSpTo
    type(field_4d_data_ptr), pointer :: ptrRelHumid, ptrTemp
    real :: RelHumid, tmpr
    type(Toptical_density_rules), intent(in) :: rulesOptDns
    real, dimension(max_species) :: arrOptColDepth



    mdl_in_q => met_buf%buffer_quantities
    if(rulesOptDns%ifRelHumidDep)then 
      iQ = fu_index(mdl_in_q, relative_humidity_flag)
      if(iQ <= 0)then
        call set_error('No relative humidity','make_optical_column_depth_map')
        return
      endif
      ptrRelHumid => met_buf%p4d(iQ)
    endif
    if (rulesOptDns%ifTDep)then
      iQ = fu_index(mdl_in_q, temperature_flag)
      if(iQ <= 0)then
        call set_error('No temperature','make_optical_column_depth_map')
        return
      endif
      ptrTemp => met_buf%p4d(iQ)
    endif


     nSpeciesOpt = pOpticColDepth%nSpecies
     !$OMP PARALLEL default(none) private(ix,iy,isrc,arrOptColDepth, iLev, relHumid, tmpr, iSpeciesTr, iSpeciesOpt, iSpTo) &
     !$OMP & shared(pOpticColDepth, pOpticDensity, pDispFlds, ptrRelHumid, &
     !$OMP & ptrTemp, ifInterpMet2DispHoriz, ifInterpMet2DispVert, &
     !$OMP & interpCoefMet2DispHoriz, interpCoefMet2DispVert, met_buf, &
     !$OMP & ChemRunSetup,rulesOptDns, nSpeciesOpt, nx_meteo)

     ! Now fill it from the transport cocktail
     !

     !
     if(defined(pOpticDensity))then
       !
       ! Just sum-up the optical density 
       !!!!!!!!!!!!Safe now, as optical depth and density are computed both for all species. 
       !!!!!!!!!!!!Later some checking has to be introduced!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       !$OMP DO collapse(2)
       do iy = 1, pOpticColDepth%ny
        do ix = 1, pOpticColDepth%nx
         do iSrc  = 1, pOpticColDepth%nSrc       !
          call sum_opt_dns_to_col_depth(pOpticDensity, arrOptColDepth, ix, iy, iSrc) ! 
          pOpticColDepth%arM(1:nSpeciesOpt, iSrc, 1, ix, iy) = arrOptColDepth(1:nSpeciesOpt)
          end do ! iSrc
        end do ! ix
       end do ! iy   
       !$OMP END DO


     else
       !
       ! Have to compute everything anew
       !
       !$OMP DO collapse(2)
       do iy = 1, pOpticColDepth%ny
        do ix = 1, pOpticColDepth%nx
         do iSrc  = 1, pOpticColDepth%nSrc

             
             arrOptColDepth(1:nSpeciesOpt) = 0

                 
             do iLev = 1, pDispFlds%n3D
               if(rulesOptDns%ifRelHumidDep)then
                 relHumid = fu_get_value(ptrRelHumid, nx_meteo, ix, iy, iLev, met_buf%weight_past, &  
                                       & interpCoefMet2DispHoriz, interpCoefMet2DispVert, &
                                       & ifInterpMet2DispHoriz, ifInterpMet2DispVert)
               endif
               if (rulesOptDns%ifTDep)then
                 tmpr = fu_get_value(ptrTemp, nx_meteo, ix, iy, iLev, met_buf%weight_past, &  
                                   & interpCoefMet2DispHoriz, interpCoefMet2DispVert, &
                                   & ifInterpMet2DispHoriz, ifInterpMet2DispVert)
               endif

               do iSpeciesTr = 1, pDispFlds%nSpecies  ! scan the transport cocktail
                 !
                 ! Contribute to the right optical species pointed by the ChemRunSetup
                 !
                 do iSpTo = 1, chemRunSetup%refsTransp2opt(iSpeciesTr)%nrefSpecies   
                   iSpeciesOpt = chemRunSetup%refsTransp2opt(iSpeciesTr)%indSpeciesTo(iSpTo)
                  
                   arrOptColDepth(iSpeciesOpt) = arrOptColDepth(iSpeciesOpt) + &  
                                              &  pDispFlds%arM(iSpeciesTr,iSrc,iLev,ix,iy) * &
                                              &  fu_Extinction_coef(iSpeciesOpt, relHumid, tmpr )
                 end do  ! iSpTo for the given iSpeciesTr
               end do  ! iSpeciesTr - of transport cocktail
             enddo ! level
             pOpticColDepth%arM(1:nSpeciesOpt, iSrc, 1, ix, iy) = arrOptColDepth(1:nSpeciesOpt)
          end do ! iSrc
        end do ! ix
       end do ! iy
       !$OMP END DO
     endif  ! if optical density is present
     !$OMP END PARALLEL


  end subroutine make_optical_column_depth_map


   !***************************************************************************************

   real function fu_Extinction_coef(iSpeciesOpt, relHumid, tmpr )
     !
     !  Trick: Some species depend on relHumid, others on temperature, but none of BOTH
     !  Extinction_coef in m2/BasicMassUnit
     !
     integer, intent(in) :: iSpeciesOpt
     real, intent(in) :: relHumid, tmpr

     integer :: MetIndex

     if(Extinction_coef%ifRhDepSp(iSpeciesOpt))then
       fu_Extinction_coef = Extinction_coef%coef(iSpeciesOpt, iMet4RH(relHumid))
     elseif(Extinction_coef%ifTDepSp(iSpeciesOpt))then
       fu_Extinction_coef = Extinction_coef%coef(iSpeciesOpt, iMet4T(tmpr))
     else
       fu_Extinction_coef = Extinction_coef%coef(iSpeciesOpt, 1)
     endif

   end function fu_Extinction_coef
   
   !***************************************************************************************

   integer function iMet4T(tmpr) result(MetIndex)
    real, intent(in) :: tmpr
       if (tmpr < 200.) then
        MetIndex = 1
       elseif (tmpr > 300.) then
        MetIndex = 21
       else
        MetIndex = nint(((tmpr - 195.) / 5.))
       endif
   end function iMet4T

   !***************************************************************************************

   integer function iMet4RH(relHumid) result(MetIndex)
    real, intent(in) :: relHumid
       if (relHumid >= 0.99) then
         MetIndex = 21
       elseif (relHumid >= 0.95) then
         MetIndex = int((relHumid - 0.775)*100) 
       elseif (relHumid >= 0.85) then
         MetIndex = int((relHumid - 0.60)*50) 
       elseif (relHumid >= 0.70) then
         MetIndex = int((relHumid - 0.475)*33.333) 
       elseif (relHumid >= 0.50) then
         MetIndex = int((relHumid - 0.40) *25)
       else
         MetIndex = 1
       endif
    end function iMet4RH




   !***************************************************************************************

   subroutine sum_opt_dns_to_col_depth(pOpticDensity, arrOptColDepth, ix, iy, iSrc)
    ! Adds up the optical densities from pOpticDensity to the optical column depth.  
    ! 
    implicit none

    ! Imported parameters
    type(Tmass_map), intent(in) :: pOpticDensity
    integer :: ix, iy, iLev, iSrc, nSp
    real, dimension(:), intent(out) :: arrOptColDepth

     arrOptColDepth(1:pOpticDensity%nSpecies) = 0

     nSp = pOpticDensity%nSpecies
     do iLev = 1, pOpticDensity%n3D
       arrOptColDepth(1:nSp) = arrOptColDepth(1:nSp) + pOpticDensity%arM(1:nSp,iSrc,iLev,ix,iy)
     enddo

  end subroutine sum_opt_dns_to_col_depth


!***************************************************************************************

  subroutine compute_optical_dens_coef(ptrRelHumid, ptrTemp, nOptSpecies, & !pOptCocktail, &
                                & interpCoefMet2DispHoriz, interpCoefMet2DispVert, &
                                & ifInterpMet2DispHoriz, ifInterpMet2DispVert, weight_past, &
                                & rulesOptDns, ix, iy, iLev, arrCoef)
!!
!! Deals with species_optic, should not be used in new developments



    implicit none

    ! Imported parameters

    type(field_4d_data_ptr), pointer :: ptrRelHumid, ptrTemp
    type(THorizInterpStruct), pointer :: interpCoefMet2DispHoriz
    type(TVertInterpStruct), pointer :: interpCoefMet2DispVert
    logical, intent(in) :: ifInterpMet2DispHoriz, ifInterpMet2DispVert
    real, intent(in) :: weight_past
    integer, intent(in) :: nOptSpecies
    integer :: MetIndex, iSpeciesOpt, ix, iy, iLev !, iSubstOpt
    real, dimension(:), intent(out) :: arrCoef
    real :: relHumid, tmpr !, fLayerThickness
!    type(silam_species), dimension(:), intent(in) :: pOptSpecies
!    type(silam_cocktail), intent(in) :: pOptCocktail
    type(Toptical_density_rules), intent(in) :: rulesOptDns


    !fLayerThickness = fu_layer_thickness_m(fu_level(dispersion_vertical,iLev))
    arrCoef(1:nOptSpecies) = 0
    relHumid = 50.
    tmpr = 273.
    if(rulesOptDns%ifRelHumidDep)then
         relHumid = fu_get_value(ptrRelHumid, nx_meteo, ix, iy, iLev, weight_past, &  
                               & interpCoefMet2DispHoriz, interpCoefMet2DispVert, &
                               & ifInterpMet2DispHoriz, ifInterpMet2DispVert)
    endif
    if (rulesOptDns%ifTDep)then
       tmpr = fu_get_value(ptrTemp, nx_meteo, ix, iy, iLev, weight_past, &  
                         & interpCoefMet2DispHoriz, interpCoefMet2DispVert, &
                         & ifInterpMet2DispHoriz, ifInterpMet2DispVert)
    endif

    do iSpeciesOpt = 1, nOptSpecies
      arrCoef(iSpeciesOpt) =  fu_Extinction_coef(iSpeciesOpt, relHumid, tmpr )
    end do 

  end subroutine compute_optical_dens_coef



!****************************************************************************************************'



     SUBROUTINE MIECALC (X, MN, NTERMS, QEXT, qscat, qback)

!   Calculates the volume extinction coefficient Qext, given size parameter X,
!   complex refractive index MN and parameter NTERMS. 
!   (Bohren & Huffman, 1984) 

      IMPLICIT   NONE
      INTEGER, intent(in) :: NTERMS
      real(r8k), intent(in) :: X
      COMPLEX*16, intent(in) :: MN
      real(r8k), intent(out) :: QEXT, qscat, qback
      INTEGER   N, NN
      COMPLEX*16  A(NTERMS), B(NTERMS), D(NTERMS+15)    
      real(r8k)      PSIN, PSIM, CHIN, CHIM, TMP
      COMPLEX*16  M, Y, XIN, XIM, CTMP,  sum3
      real(r8k)     SUM1, sum2



!    Generate the Dn's by down recurrence  D = d(log(PSI(y)))/dy

      M = DCONJG(MN)
      Y = M*X
      NN = NTERMS + 15
      D(NN) = DCMPLX (0.0D0, 0.0D0)
      DO N = NN, 2, -1
          D(N-1) = N/Y - 1.0/ (D(N) + N/Y)
      END DO


!           Generate the PSIn's and XIn'S by upward recurrence
!           and calculate the An's and Bn's from them.
!           (PSIN = PSI(n), PSIM = PSI(n-1), same for CHI)
      PSIM = DCOS(X)
      PSIN = DSIN(X)
      CHIM = -DSIN(X)
      CHIN = DCOS(X)

      DO N = 1, NTERMS

          TMP = PSIN
          PSIN = (2*N-1)/X *PSIN - PSIM
          PSIM = TMP
          TMP = CHIN
          CHIN = (2*N-1)/X *CHIN - CHIM
          CHIM = TMP
          XIN = DCMPLX (PSIN, -CHIN)
          XIM = DCMPLX (PSIM, -CHIM)
          CTMP = D(N)/M + N/X
          A(N) = (CTMP*PSIN - PSIM) / (CTMP*XIN - XIM)
          CTMP = M*D(N) + N/X
          B(N) = (CTMP*PSIN - PSIM) / (CTMP*XIN - XIM)

      END DO

      SUM1 = 0.
      SUM2 = 0.
      SUM3 = 0.

      DO N = 1, NTERMS
          SUM1 = SUM1 + (2*N+1)*( DREAL(A(N)) + DREAL(B(N)) )
          SUM2 = SUM2 + (2*N+1)*( DREAL(A(N)*DCONJG(A(N))) &
     &                          + DREAL(B(N)*DCONJG(B(N))) )
          SUM3 = SUM3 + (2*N+1)*(-1)**N * (A(N) - B(N))
      ENDDO
      QEXT = 2.0D0/X**2 * SUM1
      QSCAT = 2.0D0/X**2 * SUM2
      QBACK = CDABS(SUM3)**2 / X**2


      
     END SUBROUTINE MIECALC

     !*******************************************************************************

     subroutine Mie_Approx (x, n, k, QEXT)

     ! Calculates the volume extinction coefficient Qext, given size parameter x and
     ! components of complex refractive index n - ik (Fournier & Evans, 1991) 

     implicit none
     real ::  X, QEXT, n, k   
     real ::  Qext_HA, Qext_R, F1, F2, fTmp1, fTmp2, fTmp3


     call msg('')
     call msg('x: ', x)
     call msg('n: ', n)
     call msg('k: ', k)

     fTmp1 = (n**2. + k**2.)**2.
     fTmp2 = n**2. + k**2.
     fTmp3 = n*k

     F1 = fTmp1 + 4. * fTmp2 + 4.
     call msg('F1: ', F1)
     F2 = 4. * fTmp1 + 12. * fTmp2 + 9.
     call msg('F2: ', F2)

     Qext_R = 24. * fTmp3 * x / F1 + &
           & (4./15. * fTmp3 + 20./3. * fTmp3 / F2 + 4.8 * fTmp3 * (7. * fTmp1 + 4. * (fTmp2 - 5.)) / F1**2.) * x**3. + &
           &  8./3. * (((fTmp1 + fTmp2 - 2.)**2 - 36.* fTmp3**2.) / F1**2.) * x**4.
     call msg('Qext_R: ', Qext_R)
!  QEXT = Qext_R

     !rho
     F1 = 2. * x * (n - 1.)
     call msg('rho: ', F1)
     !beta
!     F2 = atan(2. * n * k * x / F1)
     F2 = atan(k /(n - 1))
     call msg('beta: ', F2)

     fTmp1 = 4. / F1**2.
     if(k == 0.)then
       call msg_warning('k is 0')
       Qext_HA = 2. - fTmp1 * ( sin(F1) * F1 + 1. - cos(F1))
     else
       fTmp2 = exp(-F1*tan(F2))
       fTmp3 = cos(F2)
       call msg('fTmp1: ', fTmp1)
       call msg('fTmp2: ', fTmp2)
       call msg('fTmp3: ', fTmp3)
       call msg('fTmp1 * fTmp3 * (F1 * fTmp2 * sin(F1 - F2): ', fTmp1 * fTmp3 * (F1 * fTmp2 * sin(F1 - F2)))
       call msg('fTmp2 * fTmp3 * cos(F1 - 2. * F2): ', fTmp2 * fTmp3 * cos(F1 - 2. * F2))
       call msg('fTmp3 * cos(2.*F2): ', fTmp3 * cos(2.*F2))


       Qext_HA = 2. - &
               & fTmp1 * fTmp3 * (F1 * fTmp2 * sin(F1 - F2) - &
                               &  fTmp2 * fTmp3 * cos(F1 - 2. * F2) + &
                               &  fTmp3 * cos(2*F2))
     endif

     call msg('Qext_HA: ', Qext_HA)
!   QEXT = Qext_HA
     if(Qext_HA<0)then
       call msg_warning('Qext_HA<0.')
       Qext_HA = -Qext_HA
       call msg('Qext_HA: ', Qext_HA)
     endif


     ! A
     fTmp1 = 0.5 + n - 1. - 2./3. * sqrt(k) - k/2. + (n - 1. + 2./3. * (sqrt(k) - 5.* k))**2. 
     call msg('A: ', fTmp1)
     ! mu
     fTmp2 = 3./5. - 3./4. * sqrt(n - 1.) + 3. * (n - 1.)**4. + 25. / (6. + 5. * (n - 1.) / k) 
     call msg('mu: ', fTmp2)
     
     ! P
     F1 = fTmp1 + fTmp2/x
     call msg('P: ', F1) 
     ! T
     F2 = 2. - exp(-x**(-2./3.))
     call msg('T: ', F2)

     QEXT = Qext_R  * (1. + (Qext_R / (F2 * Qext_HA))**F1)**(-1./F1)
!     Qext = (1./(Qext_R**F1) + 1./((F2 * Qext_HA)**F1))**(-1./F1)


     call msg('(Qext_R**F1): ', (Qext_R**F1))
     call msg('((F2 * Qext_HA)**F1): ', ((F2 * Qext_HA)**F1))
     call msg('1./(Qext_R**F1):', 1./(Qext_R**F1))
     call msg('1./((F2 * Qext_HA)**F1): ', 1./((F2 * Qext_HA)**F1))
     call msg('QEXT: ', QEXT)
      
if(QEXT .eps. 0.)then
   call msg('QEXT .eps. 0.')
endif

     end subroutine Mie_Approx

END MODULE optical_density

