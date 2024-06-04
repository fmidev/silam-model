MODULE chem_dep_sulphur_dmat
  !
  ! This module contains full description of the DMAT-based SOx chemistry and deposition.
  ! To be exact, there is a description of the SO2 and SO4=, removal computation - both
  ! dry and wet.
  !
  ! One type, which explains how this stuff has to be treated is Tchem_rules_S_DMAT,
  ! which contains all necessary rules for the operation, including the types
  ! of the functions used for the computations. This type is used 
  ! only here and is not seen from outside, except for the command "create".
  !
  ! All units: SI, unless otherwise stated
  !
  ! Code owner Mikhail Sofiev mikhail.sofiev@fmi.fi
  !
  ! Language: ANSI FORTRAN-90 (or close to it)
  !
  use cocktail_basic
  use hydroxyl
  use depositions
  
  implicit none
  private

  !
  ! PUBLIC routines of sulphur DMAT cocktail 
  !
  
  public init_chemicals_dmat
  public sulphur_dmat_input_needs
  public transform_dmat
  
  ! routines of Tchem_rules_DMAT_S 
  public set_chem_rules_sulphur_dmat
  public set_missing
  public inventory_dmat
  public register_species_dmat
  public fu_lifetime_DMAT_S
  public fu_if_specific_deposition

  !
  ! Private routines of the sulphur DMAT cocktail and chemical rules
  !
  private set_missing_chem_rules_DMAT_S
  private fu_if_specific_dep_DMAT_S

  interface set_missing
    module procedure set_missing_chem_rules_DMAT_S
  end interface

  interface fu_if_specific_deposition
    module procedure fu_if_specific_dep_DMAT_S
  end interface

  public fu_tla_size_dmat

  !--------------------------------------------------------------------
  !
  !  Tchem_rules_DMAT_S type definition
  !  Determines the rules of operating with the SOx following the DMAT chemistry
  !
  type Tchem_rules_DMAT_S
    private

!    real :: vdSO2, vdSO4
!    integer :: nPrecFlds
    real, dimension(:,:), pointer :: SO2_OH_2_SO3
    real :: conv_SO2_2_SO4_aqua_basic
    logical :: use_clim_oh = .false.
    logical :: use_mm_oh = .false. ! Use OH from MassMap
    !    real :: massLowThreshold
    logical :: use_detailed_sulphur_solubility = .false.
    type(silja_logical) :: defined
  end type Tchem_rules_DMAT_S
  public Tchem_rules_DMAT_S

  !--------------------------------------------------------------------------
  !
  ! DMAT sulphur chemistry knows only two species SO2 and SO4.
  ! We will point them by iSO2 and iSO4
  ! Transformation and deposition needs meteodata pointed by indices in meteo buffer
  !
  integer, private, save :: iSO2, iH2SO4_sl, iSO4w_sl, iOH = int_missing
  integer, private, pointer, save :: ind_tempr, ind_tempr_2m, ind_RH, ind_hgt, ind_cloud_cvr, &
       & ind_ablh, ind_cwcabove

  !
  ! Stuff needed for computations of transformation and deposition.
  ! Pointers are set during prepare_DMAT_S_transdep
  !
  REAL, DIMENSION(:), private, save, POINTER :: ptrTempr2m, ptrFricVel, ptrVdCorrection
  type(field_4d_data_ptr), private, save, pointer :: ptrTempr, ptrScavStd, ptrRelHumidity

  integer, private, parameter :: nDMAT_S_species_tr = 1
  integer, private, parameter :: nDMAT_S_species_short_lived = 2
  type(silam_species), dimension(nDMAT_S_species_tr), private, target, save :: species_dmat_tr
  type(silam_species), dimension(nDMAT_S_species_short_lived), private, target, save :: &
                                                                          & species_dmat_short_lived

  !
  ! Label of this module
  !
  integer, parameter, public :: transformation_sulphur_dmat = 5005
  !
  ! Size of the lookup chemistry coefficient table
  !
  integer, private, parameter :: nIndZ = 121
  integer, private, parameter :: nIndT = 301


CONTAINS
  

  !************************************************************************************

  subroutine init_chemicals_dmat()
    !
    ! Initialises the chemicals for DMAT_S transformation
    !
    implicit none

    ! Let's first define the DMAT-SOX internal species: SO2, H2SO4 as a result of gas-phase oxidation, 
    ! and SO4 in cloud-processed aersool mode as a result of aquaus-phase oxidation. The corresponding 
    ! lognormal mode is roughly projected to a fixed-diameter nib
    !
    call set_species(species_dmat_tr(1), fu_get_material_ptr('SO2'), in_gas_phase)
    call set_species(species_dmat_short_lived(1), fu_get_material_ptr('H2SO4'), in_gas_phase)
    call set_species(species_dmat_short_lived(2), fu_get_material_ptr('SO4'), in_gas_phase) ! &
                                !  & fu_set_mode(fixed_diameter_flag, 0.5e-6, 2.0e-6, 1.11e-6))

  end subroutine init_chemicals_dmat


  !************************************************************************************

  subroutine inventory_dmat(rules, &
                                   & speciesEmis, speciesTransp, speciesShortlived, speciesAerosol,&
                                   & nSpeciesEmis, nSpeciesTransp, nspeciesShortlived, nspeciesAerosol, &
                                   & iClaimedSpecies, ifActive)
    implicit none
    type(Tchem_rules_DMAT_S), intent(in) :: rules
    type(silam_species), dimension(:), pointer :: speciesEmis, speciesTransp, speciesShortlived, &
                                                & speciesAerosol
    integer, intent(in) :: nSpeciesEmis
    integer, intent(inout) :: nSpeciesTransp, nSpeciesShortlived, nSpeciesAerosol
    integer, dimension(:), intent(inout) :: iClaimedSpecies
    logical, intent(out) :: ifActive

    ! Local variable    
    integer :: iEmis, iCount
    
!    speciesTransp => species_dmat_tr
!    speciesShortLived => species_dmat_short_lived
!
!    nSpeciesTransp = nDMAT_S_species_tr
!    nSpeciesShortlived = nDMAT_S_species_short_lived
!    nullify(speciesAerosol)
!    nSpeciesAerosol = 0

    !
    ! Presence of emission species in the list of transported DMAT-S species means ownership
    ! Presence of emission species in short-living DMAT-S species means error
    !
    iCount = 0
    do iEmis = 1, nSpeciesEmis
      if(fu_index(speciesEmis(iEmis), species_dmat_tr, nDMAT_S_species_tr) > 0)then
        if(iClaimedSpecies(iEmis) < 0)then
          call msg('DMAT_S owns:' + fu_substance_name(speciesEmis(iEmis)))
          iClaimedSpecies(iEmis) = transformation_sulphur_dmat
          iCount = iCount + 1
        else
          call msg('Cannot claim ownership because he claimed it already:',iClaimedSpecies(iEmis))
          call set_error('Cannot claim ownership because someone claimed it already','species_for_transf_dmat')
          return
        endif
      endif  ! emission species belongs to DMAT-S transport list
    end do  ! emission species
    !
    ! If something useful is emitted, we can proceed with transformations
    !
    !if(iCount > 0)then
    if(.true.)then
      call addSpecies(speciesTransp, nSpeciesTransp, species_dmat_tr, nDMAT_S_species_tr, .true.)
      call addSpecies(speciesShortlived, nSpeciesShortlived, species_dmat_short_lived, &
                                                              & nDMAT_S_species_short_lived, .true.)
      ifActive = .not. error
    else
      call set_error('Cannot find anything useful in emission list','species_for_transf_dmat')
      call unset_error('species_for_transf_dmat')
      ifActive = .false.
    endif

  end subroutine inventory_dmat


  !************************************************************************************

  subroutine register_species_dmat(rules, &
                                 & speciesTransp, speciesShortlived, speciesAerosol,&
                                 & nspeciesTransp, nspeciesShortlived, nspeciesAerosol)
    implicit none
    type(Tchem_rules_DMAT_S), intent(in) :: rules
    type(silam_species), dimension(:), pointer :: speciesTransp, speciesShortlived, speciesAerosol
    integer, intent(in) :: nspeciesTransp, nspeciesShortlived, nspeciesAerosol

    ! Local variables
    !
    integer :: iEmis
    type (silam_species) :: SPtmp


    iSO2 = fu_index(species_dmat_tr(1), speciesTransp, nspeciesTransp)
    if (iSO2 < 1) then
      call set_error('No SO2 in transport species', 'register_species_dmat')
      return
    end if

    iH2SO4_sl = fu_index(species_dmat_short_lived(1), speciesShortLived, nspeciesShortLived)
    if (iH2SO4_sl < 1) then
      call set_error('No H2SO4 gas-phase in short-living species', 'register_species_dmat')
      return
    end if

    iSO4w_sl = fu_index(species_dmat_short_lived(2), speciesShortLived, nspeciesShortLived)
    if (iSO4w_sl < 1) then
      call set_error('No SO4 aqueou-phase in short-living species', 'register_species_dmat')
      return
    end if

    if (rules%use_mm_oh) then
       call set_species(SPtmp, fu_get_material_ptr('OH'), in_gas_phase)
       iOH = fu_index(SPtmp, speciesTransp, nspeciesTransp)
       if (fu_fails(iOH > 0, "OH not found", 'register_species_dmat')) return
    endif


  end subroutine register_species_dmat

  
  !************************************************************************************

  subroutine sulphur_dmat_input_needs(rules, meteo_input_local)
    implicit none
    type(Tchem_rules_DMAT_S), intent(in) :: rules
    type(Tmeteo_input), intent(out), target :: meteo_input_local

    integer :: allocstat, nq

    ! Need 3D and 2m temperature for transformation.
    !
    if (.not. rules%defined == silja_true) then
      return
    end if

    meteo_input_local =  meteo_input_empty

    if (rules%use_detailed_sulphur_solubility) then
      meteo_input_local%nquantities = 4

      meteo_input_local%quantity(1) = temperature_flag
      meteo_input_local%quantity(2) = cwcabove_3d_flag
      meteo_input_local%quantity(3) = height_flag
      meteo_input_local%quantity(4) = total_cloud_cover_flag

      meteo_input_local%q_type(1:4) = meteo_dynamic_flag
      
      ind_tempr => meteo_input_local%idx(1)
      ind_cwcabove => meteo_input_local%idx(2)
      ind_hgt => meteo_input_local%idx(3)
      ind_cloud_cvr => meteo_input_local%idx(4)
    else
      meteo_input_local%nquantities = 6

      meteo_input_local%quantity(1) = temperature_flag
      meteo_input_local%quantity(2) = temperature_2m_flag
      meteo_input_local%quantity(3) = relative_humidity_flag
      meteo_input_local%quantity(4) = height_flag
      meteo_input_local%quantity(5) = total_cloud_cover_flag
      meteo_input_local%quantity(6) = abl_height_m_flag

      meteo_input_local%q_type(1:6) = meteo_dynamic_flag
    
      ind_tempr => meteo_input_local%idx(1)
      ind_tempr_2m => meteo_input_local%idx(2)
      ind_RH => meteo_input_local%idx(3)
      ind_hgt => meteo_input_local%idx(4)
      ind_cloud_cvr => meteo_input_local%idx(5)
      ind_ablh => meteo_input_local%idx(6)
    end if

  end subroutine sulphur_dmat_input_needs


  !***********************************************************************


  subroutine transform_dmat(vSp, vSp_SL, rules, metdat_col, i3d, n3d, zenith_cos, now, lat, lon, &
                      & timestep_sec, vtla, cell_volume, cell_area, species, nSpecies)
    !
    ! Implements conversion SO2 -> SO4 as it is done in DMAT
    !
    implicit none

    ! Imported parameters
    real, dimension(:), intent(inout) :: vSp, vSp_SL
    type(Tchem_rules_DMAT_S), intent(in) :: rules
    real, dimension(:,:), intent(in) :: metdat_col
    type(silja_time), intent(in) :: now
    real, intent(in) :: zenith_cos, timestep_sec, lat, lon, cell_volume, cell_area
    real, dimension(:), pointer :: vtla !!!Single value
    type(silam_species), dimension(:), intent(in) :: species
    integer, intent(in) :: nSpecies, i3d, n3d

    ! Local variables
    real(kind=8) :: fTmp, dSO2_air, dSO2_water, rate_air, rate_water, dSO2_sum, alphat
    real :: tempr, cwc, fCloudCover, fHeight
    real :: sun, sunD, cOH_forced, scaling, dSO4_air, dSO4_water, SO4w, SO4a, SO2
    real :: sulphur_in_air, sulphur_in_water
    integer :: indT, indZ
    character (len=*), parameter :: sub_name = 'transform_dmat'

    if (associated(vtla)) then
      if (fu_fails(size(vtla) == 1, "Wrong size TLA", sub_name)) return
    else
      if(vSp(iSO2) == 0. .and. timestep_sec > 0) return
    endif

    tempr = metdat_col(ind_tempr, i3d)
    fHeight = metdat_col(ind_hgt, i3d)
    fCloudCover = min(0.999,max(metdat_col(ind_cloud_cvr, i3d), 1e-2))
    
    if (rules%use_detailed_sulphur_solubility) then
      if (i3d == n3d) then
        cwc = metdat_col(ind_cwcabove, i3d)
      else
        cwc = metdat_col(ind_cwcabove, i3d) - metdat_col(ind_cwcabove, i3d+1)
      end if

      cwc = cwc * cell_area

      if (cwc > 0.5) then !FIXME: small amount, but should depend on the cell height
       call compute_sulphur_in_water(vSp, cwc, species, nSpecies, cell_volume, &
            & tempr, fCloudCover, sulphur_in_air, sulphur_in_water)
      else
        sulphur_in_air = vSp(iSO2) * cell_volume
        sulphur_in_water = 0.0
      end if
      
    else    
      sun = zenith_cos
      if(sun < 0.) sun = 0.  ! night
      sunD = sun * (1. - fCloudCover*0.5)
      if (sunD < 0.)then
        call set_error('Negative sunshine', sub_name)
        call msg('sun, fCloudCover', sun, fCloudCover)
      endif
      if(error)return
    end if

    if (rules%use_clim_oh) then
      cOH_forced = fu_OH_cnc_clim(lat, lon, fHeight, now)
    elseif (rules%use_mm_oh) then
      if (timestep_sec>0) then
        cOH_forced = vSp(iOH)
        if (associated(vtla)) vtla(1) = cOH_forced
      else
        if (associated(vtla)) then
           cOH_forced = vtla(1)
        else
          call set_error("Can't use massmap OH in adjoint..", sub_name)
        endif
      endif
    else
      cOH_forced = fu_OH_cnc_forced(sunD, fHeight)
    end if

    if (error) return
    if ( .not. (cOH_forced >= 0) ) then
      call msg('cOH:', cOH_forced)
      call msg('lat lon hgt', (/lat, lon, fHeight/))
      call set_error('Neg cOH', sub_name)
      return
    end if

    if(vSp(iSO2) == 0. .and. timestep_sec > 0) return

    ! Get temperature from the meteo buffer
    !
    tempr = metdat_col(ind_tempr, i3d)
    if(tempr > 372 .or. tempr < 74)then
      call set_error('Strange temperature :'//trim(fu_str(tempr)),sub_name)
      return
    endif

    indT = min(int(tempr-71.65), nIndT)
    indZ = min(int(fHeight / 500.) + 1, nIndZ)

    rate_air =  rules%SO2_OH_2_SO3(indT,indZ) * cOH_forced

    if (rules%use_detailed_sulphur_solubility) then
      ! Freezing out
      if(tempr < 270.)then
        scaling = 0.
      else   
        scaling = 20.0 ! so that the new method
        ! will result in the same amount of SO2 and PM at European
        ! measurement stations as the old method
      endif

      ! aqueous conversion rate
      rate_water = rules%conv_SO2_2_SO4_aqua_basic * scaling

    else
      scaling = 0.
      ! option 1: linear to rh above 0.5
      if(metdat_col(ind_rh, i3d) > 0.5)then
        scaling = metdat_col(ind_rh, i3d) - 0.5
      endif

      ! in cloud (above abl)
      if(metdat_col(ind_hgt, i3d) > metdat_col(ind_ablh, i3d))then
        scaling = scaling + metdat_col(ind_cloud_cvr, i3d)
      endif

      ! Temperature dependence and freezing out
      if(tempr < 270.)then
        scaling = 0.
      else
        scaling = scaling * (tempr - 270.)
      end if

      ! aqueous conversion rate
      rate_water = rules%conv_SO2_2_SO4_aqua_basic * scaling
    end if

    !! Copy from advection 
    alphat = -(rate_air + rate_water) * abs(timestep_sec)
    if (abs(alphat) > 0.05_r8k) then
      fTmp =   (1.0_r8k - exp(alphat))/(rate_air + rate_water)
    else
      !! ~2e-8 max error
        fTmp =   abs(timestep_sec)*(0.999999993201951d0 - alphat*(0.499999996857818d0 - &
            & alphat* (0.166688704367924d0 - alphat*(0.0416714248986641d0))))
    endif

    if (rules%use_detailed_sulphur_solubility) then
      dSO2_air = sulphur_in_air * rate_air * fTmp / cell_volume
      dSO2_water = sulphur_in_water * rate_water * fTmp / cell_volume
    else
      dSO2_air = vSp(iSO2) * rate_air * fTmp
      dSO2_water = vSp(iSO2) * rate_water * fTmp
    end if 
    !write(15,*)fHeight,cOH_forced,indT,indZ,rules%SO2_OH_2_SO3(indT,indZ),dSO2_air,dSO2_water,vSp(iSO2),vSp_SL(iH2SO4_sl),vSp_SL(iSO4w_sl)

    if (dSO2_air + dSO2_water > vSp(iSO2)) then
      call msg('ERROR with sulphur. vSp(iSO2), dSO2_air, dSO2_water, rate_air, rate_water, fTmp, sulphur_in_air, sulphur_in_water, scaling:', &
            & (/ vSp(iSO2)+0.0_r8k, dSO2_air, dSO2_water, rate_air, rate_water, fTmp+0.0_r8k, sulphur_in_air/cell_volume+0.0_r8k, sulphur_in_water/cell_volume+0.0_r8k, scaling+0.0_r8k /)) 
    end if
    
    if (timestep_sec > 0) then
      !
      ! Forwad mode
      !
!      SO2 = vSp(iSO2) - dSO2_air - dSO2_water
!      SO4a = vSp_SL(iH2SO4_sl) + dSO2_air
!      SO4w = vSp_SL(iSO4w_sl) + dSO2_water
!      if (any( (/SO2, SO4a, SO4w/) < 0)) call ooops("")
!
      vSp(iSO2) = vSp(iSO2) - dSO2_air - dSO2_water
      vSp_SL(iH2SO4_sl) = vSp_SL(iH2SO4_sl) + dSO2_air
      vSp_SL(iSO4w_sl) = vSp_SL(iSO4w_sl) + dSO2_water

    else
      !
      ! adjoint mode
      !
      dSO4_air = vSp_SL(iH2SO4_sl) * rate_air * fTmp 
      dSO4_water = vSp_SL(iSO4w_sl) * rate_water * fTmp 
      vSp(iSO2) = vSp(iSO2) - dSO2_air - dSO2_water + dSO4_air + dSO4_water
    end if
    
  end subroutine transform_dmat



  !************************************************************************
  !
  !  DMAT_S chemistry rules
  !
  !************************************************************************

  subroutine set_chem_rules_sulphur_dmat(nlSetup, nlStdSetup, rules)
    !
    ! Creates the set of rules, in accordance with which the computations are to
    ! be made. Input - SILAM namelist
    ! Taking this chance, let's also set once-and-forever a few chemical coefficients
    !
    implicit none

    ! Imported parameter
    type(Tsilam_namelist), intent(in) :: nlSetup, nlStdSetup
    type(Tchem_rules_DMAT_S) :: rules

    ! Local variables
    integer :: iTmp, jTmp, ilev
    real :: fRelatDensity, fRelatPress, fRelatTempr, fTempr, f300_T
    real, parameter :: vres = 250.0
    character(len=fnlen) :: filename
    type(silam_vertical) :: vert_lut
    !
    ! Stupidity checking
    !
    rules%defined = silja_false
    if(.not.defined(nlSetup))then
      call set_error('Undefined namelist given','set_chem_rules_DMAT_S')
      return
    endif
    !
    ! Dry deposition for SO2 and SO4 - so far we just take them as they are DMAT.
    ! Later one can take something more sophisticated - e.g. link the whole stuff
    ! to resistance laws, etc. So far it is not supported.
    ! Another problem is that SILAM does not have the near-surface concentration profile,
    ! which essentially limits the dry deposition velocity. Therefore, I might need to
    ! cut it manually - and risk to over-state the absolute concentration levels.
    !
!    rules%vdSO2 = 0.002    ! in DMAT: 0.002  ! [m/s]
!    rules%vdSO4 = 0.0005 ! [m/s]

    !
    ! Tuning the scheme:
    ! attempt 2:
    ! vdSO2 = 0.001; vdSO4 = 0.0003
    ! attempt 3:
    ! vdSO2 = 0.001; vdSO4 = 0.0005
    !

!    rules%nPrecFlds = fu_content_int(nlSetup,'number_of_precipitation_fields')

    !
    ! SO2 to SO4 conversion as a fraction of the SO2 mass per second
    ! It is temperature-dependent and in DMAT is set as a tabulated array following 
    ! the formula:
    ! TLA = 4e-4 * 1. * 1.1419**(T-261K) + 3e-3
    ! However, according to paper, it is: TLA = 138500 * exp(-4517/T)
    !
    ! It is not too far for low temperatures but MUCH smaller than that in the code for 
    ! higher temperatures. Therefore, I replace it with a similar-type exponent, which
    ! meets the power-law-based coef at lowest and highest temperatures. In-between it 
    ! is slightly larger (up to a factor of 2).
    ! 
    ! TLA = 4.2e7 * exp(-6000/T)  [1/hour] = 1.17e4 * exp(-6000/T)  [1/sec]
    !
    ! We will set nIndT=301 value from -200 to +100 degrees Celsius, which should be enough
    !

    !
    ! Tuning the scheme:
    ! Attempt 2:
    ! As described above
    ! Attempt 3:
    ! reduced conversion rate: 1e4 *& exp(-6000/T)
    !
    allocate(rules%SO2_OH_2_SO3(nIndT,nIndZ), stat=iTmp)
    if(iTmp /= 0)then
      call set_error('Failed to allocate memory for SO2->SO4 conversion rate', &
                   & 'set_chem_rules_DMAT_S')
      return
    endif
    
    rules%conv_SO2_2_SO4_aqua_basic = 1.e-6  ! see notebook 11, pp.24-26

    !
    ! Gas-phase conversion is a three-body reaction
    !
    do iTmp = 1, 301

      fTempr = real(iTmp) + 72.15
      f300_T = 300. / fTempr

      do jTmp = 1, nIndZ

        call us_standard_atmosphere(real(jTmp-1)*500., fRelatDensity, fRelatPress, fRelatTempr)
        if(error)return

        rules%SO2_OH_2_SO3(iTmp,jTmp) = (7.42e6*fRelatDensity*(f300_T)**3.3) / &          ! n.9,p.32: nom
                   & (1.21 + 6.13*fRelatDensity*(f300_T)**3.3) * &  
                   & 0.45 ** (1./(1.+(log10(5.07*fRelatDensity*(f300_T)**3.3))**2))                      ! denominator
                  ! & 0.45 ** (1./(1.+(log(5.07*fRelatDensity*(f300_T)**3.3)*ln_10_1)**2)) ! broadening
      end do

!      rulesDMAT_S%conversion_rate(iTmp) = 10000. * exp(-6000. / (iTmp + 72.15))
    end do

    !
    ! For Eulerian scheme we will need the minimum threshold of amount of species in the grid cell
    ! The cell will not be processed by ANY Eulerian-pool routine should it have smaller amount 
    ! than that threshold. A data-based setting is possible only when the total emitted amount
    ! is known, as well as a grid size. So far let's put it to something small enough.
    !
!    rules%massLowThreshold = 1.0e-13

    rules%use_clim_oh = fu_str_u_case(fu_content(nlsetup, 'oh_param_method')) == 'CLIMATOLOGY'
    if (rules%use_clim_oh) then
      filename = fu_content(nlStdSetup, 'oh_climatology_file')
      if (fu_fails(filename /= '', 'Missing oh_climatology_file', 'set_chem_rules_sulphur_dmat')) return
      call set_vertical((/(fu_set_level(constant_height, ilev*vres), ilev = 0, int(20e3/vres))/), vert_lut)
      call setup_oh_lut(fu_process_filepath(filename), vert_lut)
      call set_missing(vert_lut, ifNew=.false.)
      if (error) return
    end if

    rules%use_detailed_sulphur_solubility = fu_str_u_case(fu_content(nlsetup, 'use_detailed_sulphur_solubility')) == 'YES'

    if (rules%use_detailed_sulphur_solubility) then
      call msg("sulphur_dmat will use detailed sulphur solubility")
    else
      call msg("sulphur_dmat will use simplified sulphur solubility")
    endif

    rules%use_mm_oh = fu_str_u_case(fu_content(nlsetup, 'oh_param_method')) == 'FROM_MASSMAP'

    if (rules%use_mm_oh) then
      call msg("sulphur_dmat will use OH from MassMap")
    elseif (.not. rules%use_clim_oh) then
      call msg("sulphur_dmat will use Forced OH")
    endif

      
    rules%defined = silja_true

  end subroutine set_chem_rules_sulphur_dmat


  !***********************************************************************

  subroutine set_missing_chem_rules_DMAT_S(rulesDMAT_S)
    implicit none
    type(Tchem_rules_DMAT_S), intent(out) :: rulesDMAT_S

    rulesDMAT_S%defined = silja_false
!    rulesDMAT_S%vdSO2 = real_missing
!    rulesDMAT_S%vdSO4 = real_missing
!    rulesDMAT_S%nPrecFlds = int_missing
!    rulesDMAT_S%massLowThreshold = real_missing
    rulesDMAT_S%conv_SO2_2_SO4_aqua_basic = real_missing
    nullify(rulesDMAT_S%SO2_OH_2_SO3)
  end subroutine set_missing_chem_rules_DMAT_S

  subroutine compute_sulphur_in_water(cnc_in_air, total_water_in_cell, species, nSpecies, volume, &
         & fTemperature, fCloudCover, sulphur_in_air, sulphur_in_water)

      implicit none

      real, dimension(:), intent(in) :: cnc_in_air
      type(silam_species), dimension(:), intent(in) :: species
      real, intent(in) :: total_water_in_cell, volume, fTemperature, fCloudCover
      integer, intent(in) :: nSpecies
      real, intent(out) :: sulphur_in_air, sulphur_in_water

      ! Local variables
      real, dimension(:), allocatable :: mass_in_water
      integer :: iSpTr, iSp, iSpeciesSA, stat
      type(Tgas_deposition_param) :: DepData
      real :: efficiency, acid_mol
      real :: equilibrium_aq_SIVfrac, fB, fC, fD ! stuff for sulphur equilibrium
      real, parameter :: H0_SO2 = 1.2e-5  !True Henry const for SO2 @ 298.15 K  [Mol/kg/Pa]
      real, parameter :: H0_SO2_Tconst = 3000 ! [K] https://webbook.nist.gov/cgi/cbook.cgi?ID=C7446095&Mask=10
      real ::  H_SO2 ! Henry const for SO2
      real :: fKs1_SO2 ! First dissociation of H2SO3
      real :: Henry_const_298K, HenryConst
      real :: total_SIV_in_cell, water_capacitance_factor, dissolved_amt, fTmp
      character(len=*), parameter :: sub_name = 'compute_sulphur_in_water'

      ! Also dissolved amounts that are not needed by the sulphur chemistry are calculated here,
      ! but perhaps this subroutine can be used for deposition later

      allocate(mass_in_water(nSpecies), stat=stat)

      mass_in_water = 0.0

        do iSp = 1, nSpecies
          !
          ! Note that SO2 has to be scavenged last because it is the only one that currently depends
          ! on the rain acidity, i.e. it is affected by other acidic and alcaline species
          !
          iSptr=iSp
          if (iSp == indSO2) iSpTr = nSpecies
          if ((indSO2 > 0) .and. (iSp == nSpecies)) iSpTr = indSO2

          if (species(iSpTr)%mode == in_gas_phase) then
            DepData = fu_gas_deposition_param(species(iSpTr)%material)
             
            Henry_const_298K = max(DepData%Henry_const_298K, 0.0)

            HenryConst = Henry_const_298K &
                    * exp (DepData%Henry_const_T_dep * ((1./fTemperature)-(1./298.)) )
          else
            HenryConst = 1e10
          end if

          water_capacitance_factor =  gas_constant_uni * fTemperature * HenryConst
          
          if(iSpTr == indSO2) then

            !if (fu_fails( timestep_sec > 0, "Negative timestep", sub_name)) return
               
            total_SIV_in_cell = volume * cnc_in_air(iSpTr) * fCloudCover
      
            acid_mol=0.
            do iSpeciesSA = 1, nStrongAcidSpecies
              acid_mol = acid_mol + mass_in_water(indStrongAcidSpecies(iSpeciesSA))  !acid is in moles
            end do
           
            ! Subtract alcalines....
            if (indNH3 > 0) acid_mol = acid_mol - mass_in_water(indNH3)

            acid_mol = acid_mol / total_water_in_cell ! Mole/kg

            if (acid_mol < 2.5e-6) acid_mol = 2.5e-6 ! pH5.6
            H_SO2 = H0_SO2 * exp(H0_SO2_Tconst*(1./fTemperature - 1./298.5)) ![Mol/kg]

            fKs1_SO2 = exp(1964.1/fTemperature - 10.91) ![mol/kg]

            fTmp = (fKs1_SO2 * gas_constant_uni * fTemperature * H_SO2) / (volume * fCloudCover)
            fB =  (acid_mol + fTmp * total_water_in_cell) * total_SIV_in_cell
            fC = - fTmp * total_SIV_in_cell

            if (fB*fB > 10000*abs(fC)  ) then ! No need for quadratic equation, Taylor series sufficient
              equilibrium_aq_SIVfrac = - fC/fB - fC*fC / (fB*fB*fB)
            else
              fD = fB*fB - 4*fC
              equilibrium_aq_SIVfrac =  0.5*(-fB + sqrt(fD))
            endif
            !
            ! In case of no strong acids (acid_mol=0) in the droplet and quadratic equation reduced to
            ! linear approximation, capasitance will be infinity. Should prevent this from happening
            !
            fTmp =  1 -  equilibrium_aq_SIVfrac * total_water_in_cell
                   !! fraction of S_IV that wants to stay in air
            if (abs(fTmp) < 1e-3 ) fTmp = 1e-3

            water_capacitance_factor =  volume * fCloudCover / total_water_in_cell / fTmp

            if (.not. water_capacitance_factor >= 0) then
              call msg('Negative water_capacitance_factor detected')
              water_capacitance_factor = 1e-3 !something small
            endif
          endif ! SO2
            
          if (water_capacitance_factor < 1e-10) cycle

          dissolved_amt = cnc_in_air(iSpTr) * water_capacitance_factor * total_water_in_cell

          dissolved_amt = min(dissolved_amt, volume * cnc_in_air(iSpTr))
          
          mass_in_water(iSpTr) =  mass_in_water(iSpTr) + dissolved_amt

          !call msg(fu_substance_name(species(iSpTr)) +  ' mass_in_air, mass_in_water, water_capacitance_factor, total_water_in_cell', &
          !     & (/ volume * cnc_in_air(iSpTr), mass_in_water(iSpTr), water_capacitance_factor, total_water_in_cell /))

        end do

      sulphur_in_water = mass_in_water(indSO2)
      sulphur_in_air = volume * cnc_in_air(indSO2)  - sulphur_in_water
          
  end subroutine compute_sulphur_in_water
  
  !********************************************************************************

  function fu_lifetime_DMAT_S(rules)result(lifetime)
    !
    ! A typical life time due to degradation and deposition for SOx materials
    ! So far, whatever the species are, we set the "typical" lifetime to be two days
    !
    implicit none

    type(silja_interval) :: lifetime
    type(Tchem_rules_DMAT_S), intent(in) :: rules

    lifetime = one_day * 2.0

  end function fu_lifetime_DMAT_S


  !**********************************************************************************

  logical function fu_if_specific_dep_DMAT_S(rulesDMAT_S)
    implicit none
    type(Tchem_rules_DMAT_S), intent(in) :: rulesDMAT_S
    fu_if_specific_dep_DMAT_S = .false.     ! no specific deposition whatsoever
  end function fu_if_specific_dep_DMAT_S

  !************************************************************************************

  integer function fu_tla_size_dmat(rules) result(n)
    implicit none
    type(Tchem_rules_DMAT_S), intent(in) :: rules
    
    n = 0
    if (rules%use_mm_oh) n  = n + 1
  end function fu_tla_size_dmat



END MODULE chem_dep_sulphur_dmat


