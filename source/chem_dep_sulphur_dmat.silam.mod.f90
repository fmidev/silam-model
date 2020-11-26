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

  public fu_if_tla_required
  private fu_if_tla_required_dmat
  interface fu_if_tla_required
     module procedure fu_if_tla_required_dmat
  end interface

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
  integer, private, pointer, save :: ind_tempr, ind_tempr_2m, ind_RH, ind_hgt, ind_cloud_cvr, ind_ablh

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

  end subroutine sulphur_dmat_input_needs


  !***********************************************************************

  subroutine transform_dmat(vSp, vSp_SL, rules, metdat, zenith_cos, now, lat, lon, timestep_sec, &
                          & low_conc_thresh, garbage, print_it)
    !
    ! Implements conversion SO2 -> SO4 as it is done in DMAT
    !
    implicit none

    ! Imported parameters
    real, dimension(:), intent(inout) :: vSp, vSp_SL
    type(Tchem_rules_DMAT_S), intent(in) :: rules
    real, dimension(:), intent(in) :: metdat
    real, dimension(:), intent(inout) :: low_conc_thresh, garbage
    type(silja_time), intent(in) :: now
    real, intent(in) :: zenith_cos, timestep_sec, lat, lon
    logical, intent(out) :: print_it

    ! Local variables
    real :: tempr, dSO2_air, dSO2_water, fTmp, dSO4_air, dSO4_water, rate_air, rate_water, &
          & fCloudCover, fHeight, sun, sunD, cOH_forced, scaling
   integer :: indT, indZ

    print_it = .false.  ! set to true and chemistry manager will give complete dump for this cell
    if(vSp(iSO2) < low_conc_thresh(iSO2) .and. timestep_sec > 0)then
      garbage(iSO2) = garbage(iSO2) + vSp(iSO2)  ! but do not touch SO4...
      vSp(iSO2) = 0.
      return
    endif
    
    ! Get temperature from the meteo buffer
    !
    tempr = metdat(ind_tempr)
    if(tempr > 372 .or. tempr < 74)then
      call set_error('Strange temperature :'//trim(fu_str(tempr)),'transform_dmat')
      return
    endif

    fCloudCover = metdat(ind_cloud_cvr)
    fHeight = metdat(ind_hgt)

    sun = zenith_cos
    if(sun < 0.) sun = 0.  ! night
    sunD = sun * (1. - fCloudCover*0.5)
    if (sunD < 0.)then
      call set_error('Negative sunshine', 'transform_dmat') 
      call msg('sun, fCloudCover', sun, fCloudCover)
    endif
    if(error)return

    indT = min(int(tempr-71.65), nIndT)
    indZ = min(int(fHeight / 500.) + 1, nIndZ)

    if (rules%use_clim_oh) then
      cOH_forced = fu_OH_cnc_clim(lat, lon, fHeight, now)
    elseif (rules%use_mm_oh) then
      cOH_forced = vSp(iOH)
      if (fu_fails(timestep_sec>0, "Can't use massmap OH in adjoint..", 'transform_dmat')) return
    else
      cOH_forced = fu_OH_cnc_forced(sunD, fHeight)
    end if
    if (error) return
    if ( .not. (cOH_forced >= 0) ) then
      call msg('cOH:', cOH_forced)
      call msg('lat lon hgt', (/lat, lon, fHeight/))
      call set_error('Neg cOH', 'transform_dmat')
      return
    end if
!    !
!    ! Compute the actual transformation rate. Unit = moles, so no problem with mass
!    !
!    fTmp = rules_DMAT_S%conversion_rate(int(tempr-71.65)) * timestep_sec 
!    dSO2 = mass_vector_tr(iSO2) * timestep_sec * &
!         & rules_DMAT_S%conversion_rate(int(tempr-71.65))  ! tempr-72.15+0.5
    !
    ! gas phase oxidation
    !
    rate_air = abs(timestep_sec) * rules%SO2_OH_2_SO3(indT,indZ) * cOH_forced
    !
    ! Aquesous-phase oxidation - see n.11,pp24-26
    !
!    scaling = 0.
    !
    ! This is for the humidity-dependent sub-saturation processing
    !
!    if(metdat(ind_rh) > 0.4 .and. metdat(ind_tempr) > 265.)then ! not too dry, neither too cold
!      if(metdat(ind_rh) > 0.95)then
!!        scaling = 4.32675  !9.0^(2/3)
!!        scaling = 2.08  !9.0^(1/3)
!        scaling = 2.42  !9.0^(0.3) * 1.25
!!        scaling = 1.55184  !9.0^(1/5)
!      else
!!        scaling = ((1.4 - metdat(ind_RH)) / (1. - metdat(ind_RH)))**0.666666666666667
!!        scaling = ((1.4 - metdat(ind_RH)) / (1. - metdat(ind_RH)))**0.3333333333333333
!        scaling = ((1.4 - metdat(ind_RH)) / (1. - metdat(ind_RH)))**0.3 * 1.25
!!        scaling = ((1.4 - metdat(ind_RH)) / (1. - metdat(ind_RH)))**0.2
!      endif
!      if(metdat(ind_tempr) < 275.)then  ! decrease the rate for freezing droplets
!        scaling = scaling * (1.-(275-metdat(ind_tempr)) / 10.0)
!      endif
!!      scaling = scaling * (1.-(285-metdat(ind_tempr)) / 20.0)
!      rate_water = rules%conv_SO2_2_SO4_aqua_basic * scaling * timestep_sec
!    else
!      rate_water = 0.0
!    endif
    !
    ! This is for the aboo-o-ve-ABL cloud-cover dependent in-cloud processing
    ! Possibly, is to be added to the below-ABL mechanism.
    !
!    if(metdat(ind_hgt) > metdat(ind_ablh) .and. metdat(ind_tempr) > 255.)then ! above ABL, not too cold
!!      scaling = 10. * metdat(ind_cloud_cvr)  ! the more clouds the faster conversion
!      scaling = 10. * sqrt(metdat(ind_cloud_cvr))  ! the more clouds the faster conversion
!!      if(metdat(ind_tempr) < 275.)then  ! decrease the rate for freezing droplets
!!        scaling = scaling * (1.-(275-metdat(ind_tempr)) / 20.0)
!!      endif
!      scaling = scaling * (1.-(285-metdat(ind_tempr)) / 30.0)
!      rate_water = rate_water +  rules%conv_SO2_2_SO4_aqua_basic * scaling * timestep_sec
!!    else
!!      rate_water = 0.0
!    endif


  ! Relative humidity, cloud cover and temperature dependent scaling 
  ! for the aqueous conversion rate
    scaling = 0.
   
    ! option 1: linear to rh above 0.5
    if(metdat(ind_rh) > 0.5)then 
      scaling = metdat(ind_rh) - 0.5
    endif
    
!    ! option 2: water in the wet particle volume
!    if(metdat(ind_rh) > 0.95)then
!      scaling = 0.5
!    elseif(metdat(ind_rh) > 0.4)then
!      scaling = (((1.4 - metdat(ind_rh)) / (1. - metdat(ind_rh))) - 1.) * 0.5 / 8.
!    endif
!  
    ! in cloud (above abl)
    if(metdat(ind_hgt) > metdat(ind_ablh))then 
      scaling = scaling + metdat(ind_cloud_cvr)
    endif

    ! Temperature dependence and freezing out
    if(metdat(ind_tempr) < 270.)then
      scaling = 0.
    else   
      scaling = scaling * (metdat(ind_tempr) - 270.)
!    elseif(metdat(ind_tempr) < 275.)then  
!      scaling = scaling * (1.-(275. - metdat(ind_tempr)) / 5.0)
!    else
!      scaling = scaling * (1.-(275. - metdat(ind_tempr)) / 5.0)
    endif     
  
    ! Constant conversion rate
!    scaling = 0.5

!   ! Lotos-Euros parameterisation. Seems to be enormous + illogical, but for a try .. 
!   ! factor 8.3 replaced with 0.5 * rules%conv_SO2_2_SO4_aqua_basic !
!    scaling = 0.5 + metdat(ind_cloud_cvr)
!    if(metdat(ind_rh) >= 0.9)then
!      scaling = scaling * (1. + 10.*(metdat(ind_rh)-0.9))
!    endif
    
    
    ! aqueous conversion rate
    rate_water = rules%conv_SO2_2_SO4_aqua_basic * scaling * abs(timestep_sec)

    !
    ! Take possible fast reduction of the SO2 mass
    !
    if (rate_air + rate_water < 1e-15) then
      ftmp = 0.0
    else
      fTmp = (1.-exp(-rate_air - rate_water)) / (rate_air + rate_water)
    end if

    dSO2_air = vSp(iSO2) * rate_air * fTmp
    dSO2_water = vSp(iSO2) * rate_water * fTmp

!write(15,*)fHeight,cOH_forced,indT,indZ,rules%SO2_OH_2_SO3(indT,indZ),dSO2_air,dSO2_water,vSp(iSO2),vSp_SL(iH2SO4_sl),vSp_SL(iSO4w_sl)

    if (timestep_sec > 0) then
      !
      ! Forwad mode
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
      !if (vsp(iso2) < 0) then
      !  print *, vSp(iSO2), dSO4_air, vSp_SL(iH2SO4_sl), dSO4_water,  vSp_SL(iSO4w_sl), dSO2_air, dSO2_water
      !end if
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
  
  logical function fu_if_tla_required_dmat(rules) result(required)
    ! Collect transformations' requests for tangent linearization. If tangent linear
    ! is not allowed, corresponding subroutine sets error. 
    implicit none
    type(Tchem_rules_dmat_s), intent(in) :: rules
    
    required = .false.

  end function fu_if_tla_required_dmat




END MODULE chem_dep_sulphur_dmat


