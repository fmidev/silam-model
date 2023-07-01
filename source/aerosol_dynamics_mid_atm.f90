
MODULE aer_dyn_middle_atmosph
  !
  ! The module contains a set of routines computing the polar stratospheric clouds largely
  ! following the FinROSE parameterization. 
  ! The module is made as an add-on to whatever chemical and other aerosol dynamic routines
  ! used in the run. All what is done here is:
  ! - computation of PSC parameters, including two tracers in the Aerosol mass map (NAT and STS)
  ! - heterogeneous reactions of chlorine and bromie activation on the surface of PSC particles
  ! Connection with other modules:
  ! - reservoir chlorine and bromine compounds come from gas-phase transformation module and 
  !   are ther active compounds are returned back after activation
  ! - PSC tracers should be sedimented: their size-related Stockes velocity is responsible 
  !   for their removal from the stratosphere: denitrification process. Note that this is 
  !   different from FinROSE, which accepted much cruder approach
  !
  ! ATTENTION
  ! Units: FinROSE uses a big zoo of the units, so SI is NOT GUARANTEED
  !        SILAM modifications are mainly SI, unless otherwise is needed for FinROSE link
  !
  ! Language: ANSI FORTRAN 90 (or something quite close to it)
  !
  ! Authors: FinROSE development team, SILAM extentions and adaptations: M.Sofiev
  !          Recent modifications 2017July->: R.Hanninen
  !
  use cocktail_basic

  implicit none
  
  !
  ! public routines for aer_dyn_simple
  !
  public init_AerDynMidAtm
  public transform_AerDynMidAtm
  public set_rules_AerDynMidAtm
  public full_spec_lst_4_AerDynMidAtm
  public registerSpecies
  public AerDynMidAtm_input_needs
  public fu_if_tla_required

  private registerSpecies_4_AerDynMidAtm
  private polar_stratospheric_clouds    ! The main PSC routine
  private TER              ! Ternary solution calculation: STS PSCs
  private NAT              ! HNO3 x 3*H2O: NAT PSC
  private ICE              ! Ice PSC
  private rates_pseudo_first_order  ! pseudo-first-order rates of heterogeneous Cl and Br reactions

  private fu_if_tla_required_MAAD

  interface fu_if_tla_required
     module procedure fu_if_tla_required_MAAD
  end interface

  interface registerSpecies
    module procedure registerSpecies_4_AerDynMidAtm
  end interface

  !!Without this define the code does not compile 
#define DEBUG_MORE  
  
  !
  ! A simplified and shortened structure of ADB. 
  ! Modifications: 
  ! (i) Three parallel bins NAT-STS-ICE
  ! (ii) only H2O, H2SO4, HNO3 are included, 
  ! (iii) just one computational regime, 
  ! (iv) no relative humidity dependence of particle size for the time being
  !
  integer, parameter, private :: nBins = 2, &             ! background (0.1um) and rock (15um)
                               & nParallelBins = 3, &     ! STS, NAT, ICE
!                               & nAerMAAD = 4, &          ! NAT, H2O, H2SO4, HNO3
                               & nAerMAAD = 7, &          ! NAT, H2O, H2SO4, HNO3, and 2Aug2017 by R.H. added HCl, HOCl, HBr.
                               & iModeFine = 1, iModeCoarse = 2, &  ! indices in masses/naero/core
                               & binSTS = 1, binNAT = 2, binICE = 3, &  
                               & iH2O_loc = 1, iH2SO4_loc = 2, iHNO3_loc = 3, iNAT_loc = 4, &
                               & iHCl_loc = 5, iHOCl_loc = 6, iHBr_loc = 7 !2Aug2017 by R.H. additions for HCl, HOCl, HBr
!RISTO Cl&Br Q: Do we need to include HCl, HOCl, HBr in the above indeces and change nAerMAAD = 7
  character(len=2), dimension(nBins), parameter, private :: chBin = (/'01','15'/)
  character(len=3), dimension(nParallelBins), parameter, private :: chParBin = (/'STS','NAT','ICE'/)
  !character(len=3), dimension(nAerMAAD), parameter, private :: chAer = (/'H2O','_S_','_N_','NAT'/)
  character(len=3), dimension(nAerMAAD), parameter, private :: chAer = (/'H2O','_S_','_N_','NAT','HCl','OCl','HBr'/) !2Aug2017 by R.H. added HCl, HOCl, and HBr
  !
  ! Aerosol modes
  type(Taerosol_mode), dimension(nBins, nParallelBins), private, save :: modesMAAD  ! (nBins)
  !
  ! Aerosol species (mass and number concentrations).
  ! Idea: transported are number concentrations, aerosol masses are just passengers
  type(silam_material), private, pointer, save :: materialMAAD_nbr
  !
  ! Indices for aerosol species in transport and aerosol-mass maps
  ! Note that numbers are transported and thus sit in transported mass map, 
  ! whereas aerosol masses are passenges in aerosol mass map
  integer, dimension(nBins,nParallelBins), private, save :: indMAAD_transp   ! in transport mass map
  integer, dimension(nAerMAAD,nBins,nParallelBins), private, save :: indMAAD_aerosol ! in aerosol mass map
  !
  ! Named indices for gaseous species in transport & shortliving maps
  ! We have (i) gas-phase transported species (further handled by CB4-strato)
  !        (ii) PSC aerosol masses for NAT and STS in transported mass map
  !       (iii) PSC aerosol number for NAT and STS in aerosol mass map
  !
  integer, private, save :: iClONO2, iClNO2, iHCl, iHOCl, iCl2, & ! gas Cl species <-> CBM_strato
                          & iN2O5, iHNO3g, &  ! gas N species <-> CBM_strato
!RISTO Cl&Br: 2Aug2017 by R.H. We included the following Br species (1 line): 
                          & iBrONO2, iBrNO2, iHOBr, iBrCl, iBr2, iHBr, & ! gas Br species  <-> CBM_strato
                          & iH2SO4_gas_sl, iSO4_gas_sl, &            ! Stuff to make STS of; from chemistry
                          & iSO4f_aer, iSO4c_aer, iNH415SO4f_aer, iNH415SO4c_aer, &
                          & iH2Og, &
                          & iNATf, iNATc, iH2Oice, iHNO3ice, iSO4ice   ! 5 transported species
  !
  ! Named indices for meteo input
  integer, private, pointer, save :: ind_tempr, ind_q, ind_pres, ind_ipv, ind_theta

  ! Label of this module
  !
  integer, parameter, public :: aerosol_dynamics_Mid_Atmosph = 5022
  !
  ! Number of temperature intervals used in precomputed coefficients
  !
  integer, private, parameter :: nIndT = 301

  ! Various stuff
  real, private, parameter :: Pi_36 = 113.09733192923255628
  integer, private, parameter :: prescribe_sizes = 6001, &  ! no size changes, nbr and mass concentrations are related 1:1
                               & prescribe_number_concentr = 6002, &  ! prescribed number conc, size and mass cnc are 1:1
                               & dynamic_spectra = 6003

!NOTE by R.H. 2017: The list below is still outdated and should be updated when all the gaseous reactions are fixed!
!                   Also, the list omits, e.g. all the bromine stuff. 
!Species available from CB4_Strato
!  subst_names = (/'CH3   ','HCO   ','CH2O  ','N2O   ','PAN   ','CRO   ', &
!       & 'TOL   ','N2O5  ','XYL   ','XO2N  ','Cl2O2 ','HONO  ', &
!       & 'PNA   ','TO2   ','ClNO2 ','ROR   ','H2O2  ','MGLY  ', &
!       & 'CO    ','CRES  ','HNO3  ','O1D   ','ETH   ','XO2   ', &
!       & 'OPEN  ','PAR   ','ClONO2','HOCl  ','C5H8  ','OLE   ', &
!       & 'C2O3  ','HCHO  ','ALD2  ','OClO  ','HCl   ','Cl2   ', &
!       & 'ClO   ','OH    ','O     ','Cl    ','O3    ','NO3   ', &
!       & 'HO2   ','NO    ','NO2   '

  !--------------------------------------------------------------------
  !
  !  Tchem_rules_AeroDynSimple type definition
  !  Totally, ten reactions are considered:
  !
!c     ( 1)  clono2 +  h2o -> hocl  +  hno3
!c     ( 2)  brono2 +  h2o -> hobr  +  hno3                     
!c     ( 3)  n2o5   +  h2o ->        2 hno3                
!c     ( 4)  clono2 +  hcl -> cl2   +  hno3          
!c     ( 5)  hocl   +  hcl -> cl2   +  h2o          
!c     ( 6)  brono2 +  hcl -> brcl  +  hno3
!c     ( 7)  hobr   +  hcl -> brcl  +  h2o 
!c     ( 8)  n2o5   +  hcl -> clno2 +  hno3
!c     ( 9)  clono2 +  hbr -> brcl  +  hno3
!c     (10)  hocl   +  hbr -> brcl  +  h2o
!c     These three heterogeneous reactions (11)-(13) are now also included:
!c     (11)  brono2 +  hbr -> br2   +  hno3
!c     (12)  hobr   +  hbr -> br2   +  h2o 
!c     (13)  n2o5   +  hbr -> brno2 +  hno3
  !
  ! Can be structured in the following way:
  !
!c                   ----------------------------     
!c                   |  h2o   |  hcl   |  hbr   |
!c     ------------------------------------------
!c     |cl: | clono2 |  h(1)  |  h(4)  |  h(9)  |
!c     |    | hocl   |   -    |  h(5)  |  h(10) |
!c     ------------------------------------------
!c     |br: | brono2 |  h(2)  |  h(6)  |  h(11) |
!c     |    | hobr   |   -    |  h(7)  |  h(12) |
!c     ------------------------------------------
!c     |n:  | n2o5   |  h(3)  |  h(8)  |  h(13) |
!c     ------------------------------------------

  !
  ! Reaction rates are calculated dynamically but probabilities
  ! of reaction in case of collision gamma are partly here.
  ! Each type of the clouds 1a NAT, 1b STS and 2 ICE have their contribution
  !
  type Tchem_rules_AerDynMidAtm
    private
    !
    ! Parameters controlling the execution
    !
    logical :: ifNAT_nucl_satur_dependent = .false. ! if #cnc dynamic, how to form NAT?
    !NOTE: The following selects the particle spectrum:
    integer :: iSpectrumDefinitionSwitch = prescribe_sizes  ! prescribe_sizes / prescribe_number_conc / dynamic_spectra
    integer :: nModesICE = 1, nModesNAT = 1, nModesSTS = 1  ! number of modes in each PSC type. =1 for FinROSE

!    real :: nbrLowThreshold                     ! low-number-concentrations threshold
!    real, dimension(:), pointer :: massLowTrsh  ! low-mass threshold
    real, dimension(max_species) :: mass2vol, mean_diam_foreign_aerosol   ! both have nSpTransport actual size

    real :: ssNAT_crit = 20., ssICE_crit = 1.4 ! supersaturation ratio thresholds
                   ! LB no supersaturation testing for Arctic winter 2004/05 run
                   !    ssnNAT_crit = 1.    ssICE_crit = 1.  ! FinROSE code: ssn, ss2
    ! Molar masses of the chemicals
    real :: mwhcl=0.036461,    mwclono2=0.097457, mwhocl=0.052460, &
          & mwbrono2=0.141908, mwhobr=0.096911,   mwn2o5=0.108009, &
          & mwh2o=0.018015,  mwhbr=0.080912, mwh2so4=0.098072, mwhno3=0.063012, &
          & mwnat=0.117057
    real :: m_hcl_clono2_ratio = 0.6116567441170119561733  ! sqrt(rules%mwhcl/rules%mwclono2)
    real :: nat1 = -2.7836, nat2 = -0.00088, nat3 = 38.9855, &
                     & nat4 = -11397.0, nat5 = 0.009179, & ! parameters for threshold temperatures for NAT
                     & ice1 = -2663.5, ice2 = 12.537, &    ! H2O satur pressure over ice, marti and mauersberger, GRL, 20 (1993) 363
                     & nice1 = 12.298, nice2 = -3968.0 ! parameters for HNO3 threshold temperatures for ICE
    real :: H2O_forcing_pressure = 40000.   ! Pa. Below that, force H2O to humidity (option 1)
    real :: H2O_forcing_ipv = 1.5e-6        ! m2 s-1 K kg-1 Below that, force H2O to humidity (option 2)
    real :: dynamic_tropopause_ipv = 2.0e-6 ! m2 s-1 K kg-1. SHows the dynamic tropopause
    real, dimension(nBins,nParallelBins) :: fixed_number_concentration = reshape((/real_missing, 0., &   ! STS
                                                                       &       1.e6 ,  230.,  & ! NAT 2300?
                                                                       &       0.   , 4.e4 /), &!ICE
                                                              & (/nBins,nParallelBins/)) 
    !NOTE!! These particle diameters should be consistent with the values introduced in init_AerDynMidAtm below:
    !real :: volume_background = 5.236e-22    ! [m3], 0.05 um in radius (0.1 um in diameter)
    !real :: volume_background = 4.1888e-21  ! [m3], 0.10 um in radius (0.2 um in diameter)
    !real :: rock_volume = 1.7671e-15         ! [m3], 7.50 um in radius (15 um in diameter)
    !real :: rock_volume = 5.8e-14           ! [m3],  ~24 um in radius (48 um in diameter)
    !real, dimension(nBins) :: bin_particle_volume = (/5.236e-22, 5.8e-14/)    ![m3] Corresponds to radii of 0.05 and 24 micrometers
    real, dimension(nBins) :: bin_particle_volume = (/5.236e-22, 1.7671e-15/) ![m3] Corresponds to radii of 0.05 and 7.5 micrometers
    real :: rock_formation_rate = 6.7e-4 ! up to 6.7e-3, [#/m3sec]

    !
    ! Parameters and coefficents to be pre-computed for temperature
    !
    real, dimension(:), pointer :: h2oeq_ice, & ! Equilibrium water pressure over ice
                                 & hno3eq_ice   ! Equilibrium HNO3 pressure over ice

    type(silja_logical) :: defined
  end type Tchem_rules_AerDynMidAtm

  !---------------------------------------------------------------------------------------------
  !
  ! A bunch of variables and coefficients is local to the grid cell: computed here and forgotten afterwards
  !
  type TLocalRates
    private
    !
    !
    real :: cnc2vmr
    !
    ! Pseudo-first-order reaction rates - for the same 13 reactions
    !
    real :: h_clono2_h2o, h_brono2_h2o, h_n2o5_h2o, h_clono2_hcl, h_hocl_hcl, h_brono2_hcl, &
          & h_hobr_hcl, h_n2o5_hcl, h_clono2_hbr, h_hocl_hbr, h_brono2_hbr, h_hobr_hbr, h_n2o5_hbr
    !
    ! PSC type 1: STS clouds. Coefficients are dynamic, have to be local
    !
    real :: gamSTS_clono2_h2o, &  ! gamSTS(1) calculated in subroutine rates_pseudo_first_order
          & gamSTS_brono2_h2o, &  ! gamSTS(2) calculated in subroutine rates_pseudo_first_order
          & gamSTS_n2o5_h2o, &    ! gamSTS(3) calculated in subroutine rates_pseudo_first_order
          & gamSTS_clono2_hcl, &  ! gamSTS(4) calculated in subroutine rates_pseudo_first_order
          & gamSTS_hocl_hcl       ! gamSTS(5) calculated in subroutine rates_pseudo_first_order
    real :: gamSTS_brono2_hcl = 0.3, &    ! gamSTS(6) = 0.9       !JLP2015: gamma=0.9   !Sessler et al.,JGR (1996) 28817: gamma=0.008
          & gamSTS_hobr_hcl = 0.2, &      ! gamSTS(7) = 0.2       !jpl1997, see jpl2000, upper limit !NOTE by R.H.: Could be refined e.g. by using the IUPAC recommendation!
          & gamSTS_n2o5_hcl = 0.003, &    ! gamSTS(8) = 0.003     !FinRose
          & gamSTS_clono2_hbr = 0.08, &   ! gamSTS(9) = 0.08      !Sessler et al.,JGR (1996) 28817.
          & gamSTS_hocl_hbr, &    ! gamSTS(10) calculated now in subroutine rates_pseudo_first_order. Used to be gamSTS(10)= 0.07    !iupac ! gamSTS(10)= 0.007
          !2Aug2017 by R.H. enabled the following reactions coefficients:
          & gamSTS_brono2_hbr = 0., &     ! n.a.
          & gamSTS_hobr_hbr = 0.25, &     ! gamSTS(12)= 0.25 !jpl1997,   gamSTS(12)= 0.025 !iupac
          & gamSTS_n2o5_hbr = 0.          ! n.a.
    !
    ! PSC type 2: ICE clouds  [Some gamma values depend on temperature. They are now set in subr rates_pseudo_first_order. OLD (ver5.5) values visible here]
    !                                  old (ver5.5) value
    real :: gamICE_clono2_h2o, &  !    gamICE(1) = 0.30       !JPL2015; NOTE: IUPAC gives equation, which could be used.
          & gamICE_brono2_h2o, &  !    gamICE(2) = 0.30       !JPL2015: gamma>0.2; NOTE: IUPAC gives equation, which could be used.
          & gamICE_n2o5_h2o,   &  !    gamICE(3) = 0.02       !JPL2015
          & gamICE_clono2_hcl, &  !    gamICE(4) = 0.30       !JPL2015; NOTE: IUPAC gives an equation which could be used.
          & gamICE_hocl_hcl,   &  !    gamICE(5) = 0.20       !JPL2015
          & gamICE_brono2_hcl, &  !    gamICE(6) = 0.30       !IUPAC??  !JPL2015: gamma 0.15, 0.34, 0.26, or even very close to 1
          & gamICE_hobr_hcl,   &  !    gamICE(7) = 0.30       !JPL2015
          & gamICE_n2o5_hcl,   &  !    gamICE(8) = 0.03       !JPL2015
          & gamICE_clono2_hbr, &  !    gamICE(9) = 0.30       !JPL2015: gamma>0.3
          & gamICE_hocl_hbr,   &  !    gamICE(10)= 0.05       !JPL2015: gamma from 0.01 to 0.38 depending on temperature and HBr concentration. 
          !2Aug2017 by R.H. enabled the following reactions coefficients:
          & gamICE_brono2_hbr, &  !    gamICE(11)= 0.0(n.a.)  !IUPAC: gamma=6.6e-4*exp(700/T); JPL2015: gamma_0 > 0.1 
          & gamICE_hobr_hbr,   &  !    gamICE(12)= 0.10       !JPL2015: gamma>0.1
          & gamICE_n2o5_hbr       !    gamICE(13)= 0.10       !iupac ; JPL2015 (Seisel et al.) gamma from 3e-3 to 0.1
    !
    ! PSC type 1: NAT clouds [Some gamma values depend on temperature. They are now set in subr rates_pseudo_first_order. OLD (ver5.5) values visible here]
    !                                  old (ver5.5) value
    real :: gamNAT_clono2_h2o, &  !    gamNAT(1) = 0.004      !JPL2015; Sessler et al. JGR1996: 0.006
          & gamNAT_brono2_h2o, &  !    gamNAT(2) = 0.006      !FinROSE, Sessler et al. JGR1996
          & gamNAT_n2o5_h2o,   &  !    gamNAT(3) = 0.0004     !JPL2015; Sessler et al. JGR1996: 0.0006    
          & gamNAT_clono2_hcl, &  !    gamNAT(4) = 0.20       !JPL2015; Sessler et al. JGR1996: 0.30
          & gamNAT_hocl_hcl,   &  !    gamNAT(5) = 0.10       !JPL2015, Sessler et al JGR1996
          & gamNAT_brono2_hcl, &  !    gamNAT(6) = 0.30       !Sessler et al.,JGR (1996) 28817.
          & gamNAT_hobr_hcl,   &  !    gamNAT(7) = 0.25       !Sessler et al.,JGR (1996) 28817.
          & gamNAT_n2o5_hcl,   &  !    gamNAT(8) = 0.003      !JPL2015, Sessler et al JGR1996
          & gamNAT_clono2_hbr, &  !    gamNAT(9) = 0.30       !JPL2015: gamma>0.3; Sessler et al JGR1996: gamma=0.30
          & gamNAT_hocl_hbr,   &  !    gamNAT(10)= 0.10       !Sessler et al.,JGR (1996) 28817.
          !2Aug2017 by R.H. enabled the following reactions coefficients:
          & gamNAT_brono2_hbr, &  !    gamNAT(11)= 0.0(n.a.)  !Sessler et al. JGR1996: 0.30
          & gamNAT_hobr_hbr,   &  !    gamNAT(12)= 0.0(n.a.)  !Sessler et al. JGR1996: 0.12
          & gamNAT_n2o5_hbr       !    gamNAT(13)= 0.005      !JPL2015    
    !
    ! Non-PSC aerosols: Al2O3, soot, etc. Have very few coefs. The values are actually set below in subroutine rates_pseudo_first_order. 
    ! For those, where we have no information, assume sulphate coating and set values for STS or ICE.
    real :: gamOT_clono2_h2o, & ! => gamSTS_clono2_h2o, &
          & gamOT_brono2_h2o, & ! => gamSTS_brono2_h2o, &
          & gamOT_n2o5_h2o, &   ! => gamSTS_n2o5_h2o, &
          & gamOT_clono2_hcl = 0.02, &  ! gamOT(4) = 0.3      !JPL2000: 0.3 -> JPL2015: ClONO+HCl on Alumina gamma=0.02 (likely typo in JPL2000/2002 where gamma=0.3)
          & gamOT_hocl_hcl,   & ! => gamSTS_hocl_hcl, &
          & gamOT_brono2_hcl, & ! => gamSTS_brono2_hcl , &
          & gamOT_hobr_hcl,   & ! => gamSTS_hobr_hcl, &
          & gamOT_n2o5_hcl,   & ! => gamSTS_n2o5_hcl, &
          & gamOT_clono2_hbr, & ! => gamSTS_clono2_hbr, &
          & gamOT_hocl_hbr,   & ! => gamSTS_hocl_hbr
    !2Aug2017 by R.H. enabled the following reactions coefficients:
          & gamOT_brono2_hbr, & ! => gamSTS_brono2_hbr, &
          & gamOT_hobr_hbr,   & ! => gamSTS_hobr_hbr, 
          & gamOT_n2o5_hbr      ! => gamSTS_n2o5_hbr, old value gamOT(13) = 0.005    !JPL2000: 0.005?? Cannot find such value any more??
  end type TLocalRates
  private TLocalRates
  
  
character(len=1024), save :: strTmp, strMetaMasses, strMetaNaero

  CONTAINS


  !***********************************************************************

  subroutine init_AerDynMidAtm(rulesAerDynMidAtm, species_lst_initial)

    ! Set modes, transport- and short-living species arrays (none so far)

    implicit none

    ! Imported parameter
    type(Tchem_rules_AerDynMidAtm), intent(inout), target :: rulesAerDynMidAtm
    type(silam_species), dimension(:), intent(out), allocatable :: species_lst_initial

    ! local variables
    integer :: iBin, iAer, iParBin

    !NOTE: These sizes should be consistent with the bin_particle_volume that is set at the begin of the module!
    modesMAAD(iModeFine,binICE) = fu_set_mode(moving_diameter_flag, 0.09e-6, 0.11e-6, &
                                            & label = 'ICE', solubility = 1)
    modesMAAD(iModeCoarse,binICE) = fu_set_mode(moving_diameter_flag, 14e-6, 16e-6, &
                                              & label = 'ICE', solubility = 1)
    modesMAAD(iModeFine,binNAT) = fu_set_mode(moving_diameter_flag, 0.09e-6, 0.11e-6, &
                                            & label = 'NAT', solubility = 1)
    modesMAAD(iModeCoarse,binNAT) = fu_set_mode(moving_diameter_flag, 14e-6, 16e-6, &
                                              & label = 'NAT', solubility = 1)
    modesMAAD(iModeFine,binSTS) = fu_set_mode(moving_diameter_flag, 0.09e-6, 0.11e-6, &
                                            & label = 'STS', solubility = 1)
    modesMAAD(iModeCoarse,binSTS) = fu_set_mode(moving_diameter_flag, 14e-6, 16e-6, &
                                              & label = 'STS', solubility = 1)

    strMetaMasses = ''
    strMetaNaero = ''
    do iParBin = 1, nParallelBins
      do iBin = 1, nBins
        do iAer = 1, nAerMAAD
          strMetaMasses = strMetaMasses + chAer(iAer) + '_' + chBin(iBin) + '_' + chParBin(iParBin) + '.'
        enddo
        strMetaNaero = strMetaNaero + chBin(iBin) + '_' + chParBin(iParBin) + '....'
      end do
    end do
    
  end subroutine init_AerDynMidAtm


  !***********************************************************************

  subroutine AerDynMidAtm_input_needs(meteo_input_local)
    !
    ! Returns input needs for the transformation routines. 
    !
    implicit none

    ! Imported parameters
    type(Tmeteo_input), intent(out), target :: meteo_input_local

    ! Local variables
    integer :: iTmp

    meteo_input_local = meteo_input_empty

    meteo_input_local%nQuantities = 5
    meteo_input_local%quantity(1) = temperature_flag    
    meteo_input_local%q_type(1) = meteo_dynamic_flag
    ind_tempr => meteo_input_local%idx(1)

    meteo_input_local%quantity(2) = specific_humidity_flag
    meteo_input_local%q_type(2) = meteo_dynamic_flag
    ind_q => meteo_input_local%idx(2)

    meteo_input_local%quantity(3) = pressure_flag
    meteo_input_local%q_type(3) = meteo_dynamic_flag
    ind_pres => meteo_input_local%idx(3)

    meteo_input_local%quantity(4) = ipv_flag
    meteo_input_local%q_type(4) = meteo_dynamic_flag
    ind_ipv => meteo_input_local%idx(4)

    meteo_input_local%quantity(5) = potential_temperature_flag
    meteo_input_local%q_type(5) = meteo_dynamic_flag
    ind_theta => meteo_input_local%idx(5)

    meteo_input_local%defined = silja_true

  end subroutine AerDynMidAtm_input_needs

  
  !***********************************************************************

  subroutine full_spec_lst_4_AerDynMidAtm(rulesAerDynMidAtm, &
                                        & speciesEmis, speciesTransp, speciesSL, speciesAerosol, &
                                        & nSpeciesEmis, nSpeciesTransp, nSpeciesSL, nSpeciesAerosol, &
                                        & iClaimedSpecies)
    !
    ! Here we select the species to deal with.
    ! Note that we are NOT replacing the activity of any other aerosol dynamics. 
    ! Our goal is only to make-up the PSC aerosol and transport species, nothing more
    !
    implicit none

    ! Imported parameters
    type(Tchem_rules_AerDynMidAtm), intent(inout) :: rulesAerDynMidAtm
    type(silam_species), dimension(:), pointer :: speciesTransp, speciesEmis, speciesSL, speciesAerosol
    integer, intent(in) :: nSpeciesEmis
    integer, intent(inout) :: nSpeciesTransp, nSpeciesSL, nSpeciesAerosol
    integer, dimension(:), intent(inout) :: iClaimedSpecies

    ! Local variables
    integer :: nSelected, iTmp
    type(silam_species) :: speciesTmp
    integer, dimension(:), pointer :: indices
    type(silam_material), pointer :: materialPtrTmp

    ! H2SO4 and sulphates
    ! Chemistry puts H2SO4 to short living, we have to use it here. There should also be
    ! sulphates either in short-living or in aerosol. For the time being, however, it is 
    ! enough to ensure that at least one of these folks is around. 
    !
    call set_species(speciesTmp, fu_get_material_ptr('H2SO4'), in_gas_phase)  ! H2SO4
    if(fu_index(speciesTmp, speciesSL, nSpeciesSL) == int_missing)then
      call set_species(speciesTmp, fu_get_material_ptr('SO4'), in_gas_phase)  ! SO4 short-lived
      if(fu_index(speciesTmp, speciesSL, nSpeciesSL) == int_missing)then
        call select_species(speciesTransp, nSpeciesTransp, 'SO4', aerosol_mode_missing, &
                          & real_missing, indices, nSelected)                 ! SO4 aerosol
        if(nSelected < 1)then
          call set_error('Mid-Atmosphere aerosol dynamics requires some sulphates', &
                       & 'full_spec_lst_4_AerDynMidAtm')
          return
        endif
      endif
    endif
    !
    ! STS and NAT require HNO3.
    !
    call set_species(speciesTmp, fu_get_material_ptr('HNO3'), in_gas_phase)  ! HNO3  gas
    if(fu_index(speciesTmp, speciesTransp, nSpeciesTransp) == int_missing)then
      call set_error('Mid-Atmosphere aerosol dynamics requires HNO3 in gas phase', &
                   & 'full_spec_lst_4_AerDynMidAtm')
      return
    endif
    !RISTO Cl&Br: Added HCl, HBr, and HOCl for STS! 
    !2Aug2017 by R.H. Added HCl, HOCl and HBr for STS
    !
    call set_species(speciesTmp, fu_get_material_ptr('HCl'), in_gas_phase)  ! HCl  gas
    if(fu_index(speciesTmp, speciesTransp, nSpeciesTransp) == int_missing)then
      call set_error('Mid-Atmosphere aerosol dynamics requires HCl in gas phase', &
                   & 'full_spec_lst_4_AerDynMidAtm')
      return
   endif
    call set_species(speciesTmp, fu_get_material_ptr('HOCl'), in_gas_phase)  ! HOCl  gas
    if(fu_index(speciesTmp, speciesTransp, nSpeciesTransp) == int_missing)then
      call set_error('Mid-Atmosphere aerosol dynamics requires HOCl in gas phase', &
                   & 'full_spec_lst_4_AerDynMidAtm')
      return
   endif
    call set_species(speciesTmp, fu_get_material_ptr('HBr'), in_gas_phase)  ! HBr  gas
    if(fu_index(speciesTmp, speciesTransp, nSpeciesTransp) == int_missing)then
      call set_error('Mid-Atmosphere aerosol dynamics requires HBr in gas phase', &
                   & 'full_spec_lst_4_AerDynMidAtm')
      return
   endif
    
    !
    ! ICE requires water as a tracer.
    !
    call set_species(speciesTmp, fu_get_material_ptr('H2O'), in_gas_phase)  ! H2O  gas
    if(fu_index(speciesTmp, speciesTransp, nSpeciesTransp) == int_missing)then
      call msg_warning('Mid-Atmosphere aerosol dynamics requires H2O in gas phase. Adding...', &
                     & 'full_spec_lst_4_AerDynMidAtm')
      call addSpecies(speciesTransp, nSpeciesTransp, (/speciesTmp/), 1 ) 
      if(error)return
    endif

    !25Oct2017 by R.H. : The following is now commented since Br2 and BrNO2 are now included in gas reactions (KPP).
    !6Sep2017 by R.H.:
    !Heterogeneous reactions result gaseous Br2 and BrNO2, but these are not in species that
    !are present in normal gaseous reactions (KPP), and therefore need to be added if not present.
    !Whether this is the right place to do this, who knows???
    call set_species(speciesTmp, fu_get_material_ptr('Br2'), in_gas_phase)  ! Br2  gas
    if(fu_index(speciesTmp, speciesTransp, nSpeciesTransp) == int_missing)then
      call msg_warning('Mid-Atmosphere aerosol dynamics requires Br2 in gas phase. Adding...', &
                     & 'full_spec_lst_4_AerDynMidAtm')
      call addSpecies(speciesTransp, nSpeciesTransp, (/speciesTmp/), 1 ) 
      if(error)return
    endif
    call set_species(speciesTmp, fu_get_material_ptr('BrNO2'), in_gas_phase)  ! BrNO2  gas
    if(fu_index(speciesTmp, speciesTransp, nSpeciesTransp) == int_missing)then
      call msg_warning('Mid-Atmosphere aerosol dynamics requires BrNO2 in gas phase. Adding...', &
                     & 'full_spec_lst_4_AerDynMidAtm')
      call addSpecies(speciesTransp, nSpeciesTransp, (/speciesTmp/), 1 ) 
      if(error)return
    endif
    !END 6Sep2017 by R.H.

    !
    ! - number concentrations. Species with missing_material.
    !       (i) STS - small bin      =====> EXCLUDED: no processes, TER sub is not dynamic
    !       (ii) NAT - small and rock bins
    !       (iii) ICE - rock bin
    ! 
    allocate(materialMAAD_nbr)
    call set_material_missing(materialMAAD_nbr)
    materialPtrTmp => materialMAAD_nbr
    
    ! Having all needed species available, can set PSC species.
    ! The list of PSC species depends on the way the spectra are calculated.
    ! We have three options now:
    ! - prescribe_sizes: each bin diameter is fixed and adding-removing mass means changing the number
    !   concentrations. 
    ! - prescribe_number_conc: the number concentration of each particle type is fixed (FinROSE way)
    ! - dynamic spectra.
    ! Evidently, number concentrations as separate advectable species are needed only for dynamic
    ! spectra. If anything is prescribed, advecting the mass concentration is enough.
    !
    ! ATTENTION.
    ! For dynamic spectra, mass concentrations become passengers to the number concentrations.
    ! For any fixed spectra parameters, we do not need number concentrations advected at all.
    !
    select case(rulesAerDynMidAtm%iSpectrumDefinitionSwitch)
      case(prescribe_sizes, prescribe_number_concentr)
        !
        ! Mass concentrations only. Add to transport mass map. Note that the fine-sulphates
        ! are not added here - they are already in the transport mass map, we have just checked it.
        !
        materialPtrTmp => fu_get_material_ptr('H2O')                                          ! ICE
        call set_species(speciesTmp, materialPtrTmp, modesMAAD(iModeCoarse, binICE))
        call addSpecies(speciesTransp, nspeciesTransp, (/speciesTmp/), 1)
!        call set_species(speciesTmp, materialPtrTmp, modesMAAD(iModeFine, binSTS))
!        call addSpecies(speciesTransp, nspeciesTransp, (/speciesTmp/), 1)
        materialPtrTmp => fu_get_material_ptr('SO4')                                          ! ICE
        call set_species(speciesTmp, materialPtrTmp, modesMAAD(iModeCoarse, binICE))
        call addSpecies(speciesTransp, nspeciesTransp, (/speciesTmp/), 1)
        materialPtrTmp => fu_get_material_ptr('HNO3')                                         ! ICE
        call set_species(speciesTmp, materialPtrTmp, modesMAAD(iModeCoarse, binICE))
        call addSpecies(speciesTransp, nspeciesTransp, (/speciesTmp/), 1)
!        call set_species(speciesTmp, materialPtrTmp, modesMAAD(iModeFine, binSTS))
!        call addSpecies(speciesTransp, nspeciesTransp, (/speciesTmp/), 1)
        materialPtrTmp => fu_get_material_ptr('NAT')                                          ! NAT
        call set_species(speciesTmp, materialPtrTmp, modesMAAD(iModeFine, binNAT))
        call addSpecies(speciesTransp, nspeciesTransp, (/speciesTmp/), 1)
        call set_species(speciesTmp, materialPtrTmp, modesMAAD(iModeCoarse, binNAT))          ! NAT
        call addSpecies(speciesTransp, nspeciesTransp, (/speciesTmp/), 1)
!2017 by R.H. related to Cl&Br: No need for HCl, HOCl and HBr because generally binSTS is not needed (where these species are)

      case(dynamic_spectra)
        !
        ! Dynamic spectra will require both mass and number transported. NUmber conc are then transportable
        ! whereas mass sits in aerosol mass map and is transported as passengers.
        !
        call set_species(speciesTmp, materialPtrTmp, modesMAAD(iModeFine,binNAT))    ! small NAT
        call addSpecies(speciesTransp, nSpeciesTransp, (/speciesTmp/), 1) !no duplicates

        call set_species(speciesTmp, materialPtrTmp, modesMAAD(iModeCoarse,binNAT))  ! rock NAT
        call addSpecies(speciesTransp, nSpeciesTransp, (/speciesTmp/), 1) !no duplicates

        call set_species(speciesTmp, materialPtrTmp, modesMAAD(iModeCoarse,binICE))  ! rock ICE
        call addSpecies(speciesTransp, nSpeciesTransp, (/speciesTmp/), 1) !no duplicates

        ! add to aerosol mass map
        !
        ! - mass concentrations
        !       (i) STS - small bin, SO4, H2O, HNO3 ====> EXCLUDED: no processes, TER is not dynamic
        !       (ii) NAT - small and rock bins, NAT = 3xH2O-HNO3
        !       (iii) ICE - rock bin, H2O, HNO3, H2SO4
        !
        materialPtrTmp => fu_get_material_ptr('H2O')
        call set_species(speciesTmp, materialPtrTmp, modesMAAD(iModeCoarse, binICE))
        call addSpecies(speciesAerosol, nSpeciesAerosol, (/speciesTmp/), 1)
        materialPtrTmp => fu_get_material_ptr('SO4')
        call set_species(speciesTmp, materialPtrTmp, modesMAAD(iModeCoarse, binICE))
        call addSpecies(speciesAerosol, nSpeciesAerosol, (/speciesTmp/), 1)
        materialPtrTmp => fu_get_material_ptr('HNO3')
        call set_species(speciesTmp, materialPtrTmp, modesMAAD(iModeCoarse, binICE))
        call addSpecies(speciesAerosol, nSpeciesAerosol, (/speciesTmp/), 1)
        materialPtrTmp => fu_get_material_ptr('NAT')
        call set_species(speciesTmp, materialPtrTmp, modesMAAD(iModeFine, binNAT))
        call addSpecies(speciesAerosol, nSpeciesAerosol, (/speciesTmp/), 1)
        call set_species(speciesTmp, materialPtrTmp, modesMAAD(iModeCoarse, binNAT))
        call addSpecies(speciesAerosol, nSpeciesAerosol, (/speciesTmp/), 1)
      case default
        call msg('Strange iSpectrumDefinitionSwitch:',rulesAerDynMidAtm%iSpectrumDefinitionSwitch)
        call set_error('Strange iSpectrumDefinitionSwitch','full_spec_lst_4_AerDynMidAtm')
        return
    end select
  end subroutine full_spec_lst_4_AerDynMidAtm


  !*******************************************************************

  subroutine registerSpecies_4_AerDynMidAtm(rulesAerDynMidAtm, &
                                          & speciesTrans, speciesSL, speciesAerosol, &
                                          & nSpeciesTrans, nSpeciesSL, nSpeciesAerosol) !, &
!                                          & passenger_tickets)
    !
    ! Set the species indices
    !
    implicit none
    
    ! Imported parameters
    type(Tchem_rules_AerDynMidAtm), intent(inout) :: rulesAerDynMidAtm
    type(silam_species), dimension(:), intent(in) :: speciesTrans, speciesSL, speciesAerosol
    integer, intent(in) :: nSpeciesTrans, nSpeciesSL, nSpeciesAerosol
!    integer, dimension(:), pointer :: passenger_tickets

    ! Local variables
    type(silam_species) :: speciesTmp
    type(silam_material), pointer :: materialPtrTmp
    integer :: nSelected, iTmp
    integer, dimension(:), pointer :: indices

    indices => fu_work_int_array()
    if(error)return
    indMAAD_transp = int_missing
    indMAAD_aerosol = int_missing

    call msg("registerSpecies_4_AerDynMidAtm: nSpeciesTrans, nSpeciesSL, nSpeciesAerosol", &
                &(/nSpeciesTrans, nSpeciesSL, nSpeciesAerosol/))
    !
    ! Components of PSCs. Note many types of sulphates: they all go to STS
    !
    call set_species(speciesTmp, fu_get_material_ptr('H2SO4'), in_gas_phase)
    iH2SO4_gas_sl = fu_index(speciesTmp, speciesSL, nSpeciesSL)

    call set_species(speciesTmp, fu_get_material_ptr('SO4'), in_gas_phase)
    iSO4_gas_sl = fu_index(speciesTmp, speciesSL, nSpeciesSL)
    !
    ! Searching external sulphates is little pleasure: they can have several modes of whatever size.
    ! Note also that we have just added sulphates for ICE. So, there must be 3 sulphates
    !
    call select_species(speciesTrans, nSpeciesTrans, 'SO4', aerosol_mode_missing, &
                      & real_missing, indices, nSelected)
    iSO4f_aer = int_missing
    iSO4c_aer = int_missing
    if(nSelected > 1)then
      call set_species(speciesTmp, fu_get_material_ptr('SO4'), modesMAAD(iModeCoarse,binICE))
      do iTmp = 1, nSelected
        if(.not. (speciesTrans(indices(iTmp)) == speciesTmp))then  ! our own ICE sulphate?
          if(iSO4f_aer == int_missing)then
            iSO4f_aer = indices(iTmp)
          else
            iSO4c_aer = indices(iTmp)
          endif
        endif   ! own ICE sulphate ?
      end do  ! selected sulphates
      if(fu_mode(speciesTrans(iSO4f_aer)) > fu_mode(speciesTrans(iSO4c_aer)))then  ! fine mode must be smaller
        iTmp = iSO4f_aer
        iSO4f_aer = iSO4c_aer
        iSO4c_aer = iTmp
      endif
    else
      call msg('Do not understand the number of transported sulphates:', nSelected)
      do iTmp = 1, min(nSelected,nSpeciesTrans)
        call report(speciesTrans(indices(iTmp)))
      end do
      call set_error('Must be 3 sulphate species','registerSpecies_4_AerDynMidAtm')
      return
    endif
    !
    ! Here we unequivocally link the fine sulphates with STS. 
    !
    indMAAD_aerosol(iH2SO4_loc,iModeFine,binSTS) = iSO4f_aer
    !
    ! Note that nitric acid and water so far were considered only via equilibrium or ignored, 
    ! i.e. theoretically they need extra transported tracers
    ! for the time being, however, TER is fully equilibrium-based sub, i.e. sending HNO3 and H2O to gas
    ! tracers in-between is fully OK. Just do not forget to keep the mass budget.
!    !
!    indMAAD_aerosol(iHNO3_loc,iModeFine,binSTS) = iNO3f_aer
!    indMAAD_aerosol(iH2O_loc,iModeFine,binSTS) = iH2Of_aer
    
    ! Same thing with ammonium sulphate - but there is no way to add it to the main game - 
    ! unless the dimension of masses is increased and ternary equilibrium is computed taking NH3.
    ! Stupid, of course... but one can argue that NH3 is in negligible amounts up there.
    !
    iNH415SO4f_aer = int_missing
    iNH415SO4c_aer = int_missing
    call select_species(speciesTrans, nSpeciesTrans, 'NH415SO4', aerosol_mode_missing, &
                      & real_missing, indices, nSelected)
    select case(nSelected)
      case(0)
      case(1)
        iNH415SO4f_aer = indices(1)  ! may be, no sea salt nitrates
      case(2)
        iNH415SO4f_aer = indices(1)
        iNH415SO4c_aer = indices(2)
      case default
        call msg('Do not understand the number of transported ammonium sulphates:', nSelected)
        do iTmp = 1, min(nSelected,nSpeciesTrans)
          call report(speciesTrans(indices(iTmp)))
        end do
        call set_error('Must be 1 or 2 ammonum sulphates','registerSpecies_4_AerDynMidAtm')
        return
    end select

    call set_species(speciesTmp, fu_get_material_ptr('H2O'), in_gas_phase)
    iH2Og = fu_index(speciesTmp, speciesTrans, nSpeciesTrans) 

    call set_species(speciesTmp, fu_get_material_ptr('HNO3'), in_gas_phase)
    iHNO3g = fu_index(speciesTmp, speciesTrans, nSpeciesTrans) 

    !
    ! Stuff to react: chlorines, bromines, reactive nitrogen
    !
    call set_species(speciesTmp, fu_get_material_ptr('ClONO2'), in_gas_phase)
    iClONO2 = fu_index(speciesTmp, speciesTrans, nSpeciesTrans) 

    call set_species(speciesTmp, fu_get_material_ptr('HCl'), in_gas_phase)
    iHCl = fu_index(speciesTmp, speciesTrans, nSpeciesTrans) 

    call set_species(speciesTmp, fu_get_material_ptr('HOCl'), in_gas_phase)
    iHOCl = fu_index(speciesTmp, speciesTrans, nSpeciesTrans) 

    call set_species(speciesTmp, fu_get_material_ptr('N2O5'), in_gas_phase)
    iN2O5 = fu_index(speciesTmp, speciesTrans, nSpeciesTrans) 

    call set_species(speciesTmp, fu_get_material_ptr('Cl2'), in_gas_phase)
    iCl2 = fu_index(speciesTmp, speciesTrans, nSpeciesTrans) 

    call set_species(speciesTmp, fu_get_material_ptr('ClNO2'), in_gas_phase)
    iClNO2 = fu_index(speciesTmp, speciesTrans, nSpeciesTrans) 

    !3Aug2017 by R.H. Enabled the following gases (BrONO2, BrCl, Br2)
    call set_species(speciesTmp, fu_get_material_ptr('BrONO2'), in_gas_phase)
    iBrONO2 = fu_index(speciesTmp, speciesTrans, nSpeciesTrans) 

    call set_species(speciesTmp, fu_get_material_ptr('BrCl'), in_gas_phase)
    iBrCl = fu_index(speciesTmp, speciesTrans, nSpeciesTrans) 

    call set_species(speciesTmp, fu_get_material_ptr('Br2'), in_gas_phase)
    iBr2 = fu_index(speciesTmp, speciesTrans, nSpeciesTrans) 

    !3Aug2017 by R.H. Added the following gases (BrNO2, HBr, HOBr)
    call set_species(speciesTmp, fu_get_material_ptr('BrNO2'), in_gas_phase)
    iBrNO2 = fu_index(speciesTmp, speciesTrans, nSpeciesTrans) 

    call set_species(speciesTmp, fu_get_material_ptr('HBr'), in_gas_phase)
    iHBr = fu_index(speciesTmp, speciesTrans, nSpeciesTrans) 

    call set_species(speciesTmp, fu_get_material_ptr('HOBr'), in_gas_phase)
    iHOBr = fu_index(speciesTmp, speciesTrans, nSpeciesTrans) 

    !RISTO TEST 6Sep2017
    !call msg('indeces:    iH2Og,  iHNO3g, iClONO2,    iHCl,   iHOCl,   iN2O5,    iCl2,  iClNO2, iBrONO2,   iBrCl,    iBr2,  iBrNO2,    iHBr,   iHOBr')
    !call msg('indeces:', (/iH2Og, iHNO3g, iClONO2, iHCl, iHOCl, iN2O5, iCl2, iClNO2, iBrONO2, iBrCl, iBr2, iBrNO2, iHBr, iHOBr/) )
    !END RISTO TEST 6Sep2017

    
    !
    ! Aerosol and transport indices mapping
    !
    select case(rulesAerDynMidAtm%iSpectrumDefinitionSwitch)

      case(prescribe_sizes, prescribe_number_concentr)
        !
        ! Mass concentrations only - from transport mass map
        !
        call set_species(speciesTmp, fu_get_material_ptr('NAT'), modesMAAD(iModeFine,binNAT))
        iNATf = fu_index(speciesTmp, speciesTrans, nSpeciesTrans) 
        indMAAD_aerosol(iNAT_loc,iModeFine,binNAT) = iNATf

        call set_species(speciesTmp, fu_get_material_ptr('NAT'), modesMAAD(iModeCoarse,binNAT))
        iNATc = fu_index(speciesTmp, speciesTrans, nSpeciesTrans) 
        indMAAD_aerosol(iNAT_loc,iModeCoarse,binNAT) = iNATc
        
        call set_species(speciesTmp, fu_get_material_ptr('H2O'), modesMAAD(iModeCoarse,binICE))
        iH2Oice = fu_index(speciesTmp, speciesTrans, nSpeciesTrans) 
        indMAAD_aerosol(iH2O_loc,iModeCoarse,binICE) = iH2Oice

        call set_species(speciesTmp, fu_get_material_ptr('SO4'), modesMAAD(iModeCoarse,binICE))
        iSO4ice = fu_index(speciesTmp, speciesTrans, nSpeciesTrans)
        indMAAD_aerosol(iH2SO4_loc,iModeCoarse,binICE) = iSO4ice

        call set_species(speciesTmp, fu_get_material_ptr('HNO3'), modesMAAD(iModeCoarse,binICE))
        iHNO3ice = fu_index(speciesTmp, speciesTrans, nSpeciesTrans)
        indMAAD_aerosol(iHNO3_loc,iModeCoarse,binICE) = iHNO3ice
!RISTO Cl&Br: No need for HCl, HOCl, and HBr since STS is not required!
        
call msg('Reporting the transport species and the indMAAD in PRESCRIBE_SIZE mode')
call msg('Transport species')
do iTmp = 1, nSpeciesTrans
  call msg('Transport species:' + fu_str(iTmp) + ','+ fu_str(speciesTrans(iTmp)))
end do
call msg('MAAD indices',(/indMAAD_aerosol(:,:,:)/))
call msg('MAAD index NAT-fine-binNAT:', indMAAD_aerosol(iNAT_loc,iModeFine,binNAT))
call msg('MAAD index NAT-coarse-binNAT:',indMAAD_aerosol(iNAT_loc,iModeCoarse,binNAT))
call msg('MAAD index H2O-coarse-binICE:',indMAAD_aerosol(iH2O_loc,iModeCoarse,binICE))
call msg('MAAD index H2SO4-coarse-binICE:',indMAAD_aerosol(iH2SO4_loc,iModeCoarse,binICE))
call msg('MAAD index HNO3-coarse-binICE:',indMAAD_aerosol(iHNO3_loc,iModeCoarse,binICE))
call msg('MAAD index H2SO4-fine-binSTS:',indMAAD_aerosol(iH2SO4_loc,iModeFine,binSTS))





      case(dynamic_spectra)
        !
        ! Transported number concentrations while aerosol masses are passengers
        !
        materialPtrTmp => materialMAAD_nbr
        call set_species(speciesTmp, materialPtrTmp, modesMAAD(1,binNAT))            ! small NAT
        indMAAD_transp(iModeFine,binNAT) = fu_index(speciesTmp, speciesTrans, nSpeciesTrans)

        call set_species(speciesTmp, materialPtrTmp, modesMAAD(nBins,binNAT))            ! rock NAT
        indMAAD_transp(iModeCoarse,binNAT) = fu_index(speciesTmp, speciesTrans, nSpeciesTrans)

        call set_species(speciesTmp, materialPtrTmp, modesMAAD(nBins,binICE))            ! rock ICE
        indMAAD_transp(iModeCoarse,binICE) = fu_index(speciesTmp, speciesTrans, nSpeciesTrans)
        !
        ! ... and the indices for the aerosol species
        !
        call set_species(speciesTmp, fu_get_material_ptr('NAT'), modesMAAD(iModeFine,binNAT))
        iNATf = fu_index(speciesTmp, speciesAerosol, nSpeciesAerosol) 
        indMAAD_aerosol(iNAT_loc,iModeFine,binNAT) = iNATf

        call set_species(speciesTmp, fu_get_material_ptr('NAT'), modesMAAD(iModeCoarse,binNAT))
        iNATc = fu_index(speciesTmp, speciesAerosol, nSpeciesAerosol) 
        indMAAD_aerosol(iNAT_loc,iModeCoarse,binNAT) = iNATc
        
        call set_species(speciesTmp, fu_get_material_ptr('H2O'), modesMAAD(iModeCoarse,binICE))
        iH2Oice = fu_index(speciesTmp, speciesAerosol, nSpeciesAerosol) 
        indMAAD_aerosol(iH2O_loc,iModeCoarse,binICE) = iH2Oice

        call set_species(speciesTmp, fu_get_material_ptr('SO4'), modesMAAD(iModeCoarse,binICE))
        iSO4ice = fu_index(speciesTmp, speciesAerosol, nSpeciesAerosol)
        indMAAD_aerosol(iH2SO4_loc,iModeCoarse,binICE) = iSO4ice

        call set_species(speciesTmp, fu_get_material_ptr('HNO3'), modesMAAD(iModeCoarse,binICE))
        iHNO3ice = fu_index(speciesTmp, speciesAerosol, nSpeciesAerosol)
        indMAAD_aerosol(iHNO3_loc,iModeCoarse,binICE) = iHNO3ice

call msg('Reporting the transport species and the indMAAD in DYNAMIC mode')
call msg('Transport species')
do iTmp = 1, nSpeciesTrans
  call msg('Transport species:' + fu_str(iTmp) + ','+ fu_str(speciesTrans(iTmp)))
end do
call msg('MAAD indices',(/indMAAD_aerosol(:,:,:)/))
call msg('MAAD index NAT-fine-binNAT:', indMAAD_aerosol(iNAT_loc,iModeFine,binNAT))
call msg('MAAD index NAT-coarse-binNAT:',indMAAD_aerosol(iNAT_loc,iModeCoarse,binNAT))
call msg('MAAD index H2O-coarse-binICE:',indMAAD_aerosol(iH2O_loc,iModeCoarse,binICE))
call msg('MAAD index H2SO4-coarse-binICE:',indMAAD_aerosol(iH2SO4_loc,iModeCoarse,binICE))
call msg('MAAD index HNO3-coarse-binICE:',indMAAD_aerosol(iHNO3_loc,iModeCoarse,binICE))
call msg('MAAD index H2SO4-fine-binSTS:',indMAAD_aerosol(iH2SO4_loc,iModeFine,binSTS))





      case default
        call msg('Strange iSpectrumDefinitionSwitch:',rulesAerDynMidAtm%iSpectrumDefinitionSwitch)
        !call set_error('Strange iSpectrumDefinitionSwitch','full_spec_lst_4_AerDynMidAtm')
        call set_error('Strange iSpectrumDefinitionSwitch','registerSpecies_4_AerDynMidAtm')
        return
    end select

!   FIXME Even more dirty hack!
!    Roux: The trick below does not work!  More than one species is added in
!    global apta runs,
!    Rubbish comes there as mean_diam_foreign_aerosol and mass2vol causing NaNs
!    in sadOTHER and breaking the run
! Thus made them static arrays (not pointers and force safe values)
! Some species will left unaccounted.. Too bad... 


!!    !
!!    ! Finally, scan the transported species and select the yet-uncovered aerosols: their surface area
!!    ! will be needed for heterogeneous reactions
!!    !
!!
!!
!!    ! Dirty hack FIXME  allocate one more for NH4NO3_m_30 added later by SAD. should be done properly
!!    allocate(rulesAerDynMidAtm%mean_diam_foreign_aerosol(nSpeciesTrans + 1), &  
!!           & rulesAerDynMidAtm%mass2vol(nSpeciesTrans + 1 ), stat = iTmp)
!!    if(fu_fails(iTmp==0,'Failed mass2vol and diam allocation','registerSpecies_4_AerDynMidAtm'))return




    rulesAerDynMidAtm%mean_diam_foreign_aerosol(:) = real_missing
    rulesAerDynMidAtm%mass2vol(1:nSpeciesTrans) = real_missing

    do iTmp = 1, nSpeciesTrans
      if(fu_mode(speciesTrans(iTmp)) == in_gas_phase) cycle
      rulesAerDynMidAtm%mass2vol(iTmp) = &
                    & fu_conversion_factor(fu_basic_mass_unit(fu_material(speciesTrans(iTmp))), &
                                         & 'kg', fu_material(speciesTrans(iTmp))) / &
                    & fu_dry_part_density(fu_material(speciesTrans(iTmp)))
      !
      ! The species are either to be unequivocally connected with masses via indMAAD_aerosol
      ! or declared foreign. Among sulphates, for instance, super-um mode and all ammonium-containing
      ! species are foreign.
      !
      if(iTmp /= iSO4f_aer .and. & ! iTmp /= iSO4c_aer .and. &
!       & iTmp /= iNH415SO4f_aer .and. iTmp /= iNH415SO4c_aer .and. &
       & iTmp /= iNATf .and. iTmp /= iNATc .and. &
       & iTmp /= iH2Oice .and. iTmp /= iHNO3ice .and. iTmp /= iSO4ice)then

        rulesAerDynMidAtm%mean_diam_foreign_aerosol(iTmp) = fu_massmean_d(fu_mode(speciesTrans(iTmp)))
      endif
#ifdef DEBUG_MORE
call msg('Species:'+fu_str(iTmp)+','+fu_name(fu_material(speciesTrans(iTmp)))+', mean_diam, mass2vol', &
                             & rulesAerDynMidAtm%mean_diam_foreign_aerosol(iTmp),rulesAerDynMidAtm%mass2vol(iTmp))
call report(speciesTrans(iTmp))
#endif
    end do  ! transport species

    call free_work_array(indices)

#ifdef DEBUG_MORE
call msg('Diameters:',rulesAerDynMidAtm%mean_diam_foreign_aerosol(1:nSpeciesTrans))
call msg('Mass to volumes:',rulesAerDynMidAtm%mass2vol(1:nSpeciesTrans))
#endif

  end subroutine registerSpecies_4_AerDynMidAtm


  !***********************************************************************

  subroutine set_rules_AerDynMidAtm(nlSetup, rules)
    !
    ! A few constants and rates set in this module.
    !
    implicit none

    ! Imported parameter
    type(Tsilam_namelist), intent(in) :: nlSetup 
    type(Tchem_rules_AerDynMidAtm), intent(out) :: rules

    ! Local variables
    integer :: iT, fTempr
    !
    ! Allocate space for precomputed variables
    !
    allocate(rules%h2oeq_ice(nIndT), &          ! H2O equilibrium pressure over ICE
           & rules%hno3eq_ice(nIndT), stat=iT)  ! HNO3 equilibrium pressure over ICE
           
    if(fu_fails(iT==0,'Failed allocation of precomputed coefficients','set_rules_AerDynMidAtm'))return
    
    do iT = 1, nIndT
      fTempr = real(iT) + 72.15
      rules%h2oeq_ice(iT) = 10.**(rules%ice1 / fTempr + rules%ice2) ! p(H2O) above ICE, Pa
!      h2oeq_ice = h2oeq_ice * hnm / ppa         ! convert to # density
      rules%hno3eq_ice(iT) = 10.0**(rules%nice1 + rules%nice2 / fTempr) * 133.322  ! p(HNO3) above ICE [Pa]
    end do   ! fTempr
#ifdef DEBUG_MORE
call msg('Coefficients are set in set_rules_AerDynMidAtm')
#endif

  end subroutine set_rules_AerDynMidAtm


  !*******************************************************************

  subroutine transform_AerDynMidAtm(vCncTrn_mm, vCncSL_mm, vCncAer_mm, garbage, &
                                  & rules, fLowCncThresh, metdat, seconds, nSpeciesTrn, &
                                  & species, ix, iy, iz, print_it) 
    !NOTE: Additional ifFirst argument needed if water is initialized (modify also chemistry_manager.silam.mod.f90)
    !
    !
    ! Calculates the parameters for polar stratospheri clouds, if any, and 
    ! heterogeneous reactions specific for the stratosphere - essentially,
    ! chlorine and bromine activation.
    ! PSC parameters come from polar_stratospheric_clouds subroutine,
    ! which is partly based on hetero subroutine of FinROSE
    ! Aerosol species considered here:
    !            bin1 0.1um               bin2 15um
    !  STS     SO4, H2SO4, HNO3              -
    !  NAT       HNO3+3*H2O              HNO3+3*H2O
    !  ICE           -              SO4, H2SO4, HNO3, H2O
    !
    ! The reactions are grouped into hydrolysis reactions (1-3) and 
    ! reations with the mineral acids HCl (4-8) and HBr (9-10).
    !
    !                   ----------------------------     
    !                   |  h2o   |  hcl   |  hbr   |
    !     ------------------------------------------
    !     |cl: | clono2 |  h(1)  |  h(4)  |  h(9)  |
    !     |    | hocl   |   -    |  h(5)  |  h(10) |
    !     ------------------------------------------
    !     |br: | brono2 |  h(2)  |  h(6)  |  h(11) |
    !     |    | hobr   |   -    |  h(7)  |  h(12) |
    !     ------------------------------------------
    !     |n:  | n2o5   |  h(3)  |  h(8)  |  h(13) |
    !     ------------------------------------------
    !
    implicit none

    ! Imported parameters
    real, dimension(:), intent(inout), target :: vCncTrn_mm, vCncSL_mm, vCncAer_mm
    type(Tchem_rules_AerDynMidAtm), intent(in) :: rules
    real, dimension(:), intent(inout) :: garbage   ! holds transport species, then aerosols
    real, dimension(:), intent(in) :: metdat, fLowCncThresh
    real, intent(in) :: seconds
    type(silam_species), dimension(:), pointer :: species
    integer, intent(in) :: nSpeciesTrn, ix, iy, iz
    logical, intent(out) :: print_it
    !logical, intent(in) :: ifFirst !Needed if one want to initialize water!

    ! Local variables
    real, dimension(:), pointer :: vCncTrn, vCncSL, vCncAer
    integer :: iTmp, iBin, iS, iMat
    real :: hno3eq,ph2o,d
    real :: HNO3r,HClr,HBrr
    real :: wtH2SO4_prc, extra_sulphur_mass_STS, extra_particles_STS, dM, fTmp
    real :: sadSTS, sadNAT, sadICE, sadOTHER
    REAL, dimension(nAerMAAD, nBins, nParallelBins) :: masses, masses_tmp  ! mass conc of aerosols
    real, dimension(nBins, nParallelBins) :: naero, naero_tmp,  &   ! number cnc of aerosols
                                           & core, core_tmp        ! mean volume of single particle
    real, dimension(nParallelBins) :: density      ! densities of the particles for each type
    logical :: ifPSC, ifICE, ifNAT, ifSTS
    type(TLocalRates) :: rulesRates

    !real :: vCncTrn_iHNO3g, vCncSL_iH2SO4_gas_sl, vCncTrn_iH2Og, vCncTrn_iHCl
    !3Aug2017 by R.H. added HOCl and Hbr
    real :: vCncTrn_iHNO3g, vCncSL_iH2SO4_gas_sl, vCncTrn_iH2Og, vCncTrn_iHCl, vCncTrn_iHOCl, vCncTrn_iHBr

    
    
#ifdef DEBUG_MORE
if(print_it)then
call msg('')
call msg('mid-atm aerosol dynamics')
endif
#endif

    print_it = .false.
    !
    ! ATTENTION. 
    ! In case of prescribed spectrum parameters, we do not need number concentrations transported
    ! They will be calculated in thus subroutine from mass and prescribed parameters of the spectrum
    ! (number concentrations or bin diameters).
    ! To save the trouble, we now REDIRECT the local vectors to appropriate input vectors
    !
    select case(rules%iSpectrumDefinitionSwitch)

      case (prescribe_sizes, prescribe_number_concentr)
        vCncTrn => vCncTrn_mm  ! always: gas-phase species are in transport mass map
        vCncSL => vCncSL_mm    ! always
        vCncAer => vCncTrn_mm  ! aerosol masses are in transport mass map

      case(dynamic_spectra)   ! DOES NOT REALLY WORK YET
        vCncTrn => vCncTrn_mm  ! always: gas-phase species are in transport mass map
        vCncSL => vCncSL_mm    ! always
        vCncAer => vCncAer_mm  ! aerosol masses are in aerosol mass map

      case default
        call msg('Funny iSpectrumDefinitionSwitch:',rules%iSpectrumDefinitionSwitch)
        call set_error('Funny rules%iSpectrumDefinitionSwitch','transform_AerDynMidAtm')
        return
    end select

call check_mass_vector(vCncTrn, garbage, species, 'Start of MAAD', fLowCncThresh, nSpeciesTrn, ix, iy, iz, print_it)

    !-----------------------------------------------------------------------
    !
    ! Water tracer can be considered as a tracer only in the upper troposophere or higher
    ! We force it in the lower altitudes from specific humidity.
    ! For separation, will use abs(IPV) value: below 1e-6 is usually troposphere, above 2e-6 is 
    ! usually stratosphere. Let's take 1.5e-6 as a threshold. Note that 2e-6 will trigger PSCs
    !
!    if(metdat(ind_pres) > rules%H2O_forcing_pressure)then  ! about 7km
!    if (abs(metdat(ind_ipv)) < rules%H2O_forcing_ipv .or. ifFirst) then !R.H. 2017: used for testing (always initialize the water at first time-step) 
    if (abs(metdat(ind_ipv)) < rules%H2O_forcing_ipv) then               !R.H. 2017: use if everything works fine!!!!
      vCncTrn(iH2Og) = metdat(ind_q) * molecular_weight_air * metdat(ind_pres) / &  ! kg/kg to mole/m3
                    & (molecular_weight_water * gas_constant_uni * metdat(ind_tempr))
      !if ( ifFirst .and. abs(metdat(ind_ipv)) >= rules%H2O_forcing_ipv ) then
      !   vCncTrn(iH2Og) = 0.2*vCncTrn(iH2Og) !RH TEST: initialize the water at ~stratosphere using 20% water from meteo 
      !end if
      if(.not. vCncTrn(iH2Og) >= 0.)then
        call msg('vCncTrn(iH2Og), q, M_air, pres, M_water, R, tempr:', &
               & (/vCncTrn(iH2Og), metdat(ind_q), molecular_weight_air, metdat(ind_pres), &
                 & molecular_weight_water, gas_constant_uni, metdat(ind_tempr)/))
        call set_error('forced water is negative','transform_AerDynMidAtm')
        return
      endif
    endif

    !-----------------------------------------------------------------------
    !
    ! PSC do not exist if 
    ! - we are too high - mesosphere
    ! - we are too low - troposphere
    ! - it is too warm
    !
    if(metdat(ind_pres) < 500. .or. &          ! low pressure
!     & h2so4 <= 0. .or.  &       ! no sulphuric acid
!     & qn(3)*ppa < 1.e-3 .or. &  ! too little water vapour
!     & qn(3)*ppa > 0.2 .or. &    ! too much water vapour
     ! my idea of using IPV
     !& (metdat(ind_pres) > 6500. .and. abs(metdat(ind_ipv)) < rules%dynamic_tropopause_ipv) .or. &   ! below dynamic tropopause
     ! FinROSE code: if(pmb(k,i,j) > 400. .or. (abs(pv(k,i,j)) <= 3.5 .and. theta(k,i,j) <= 380.)) then
     & (abs(metdat(ind_ipv)) <= 3.5e-6 .and. metdat(ind_theta) <= 380.) .or. &
     & metdat(ind_pres) > 40000. .or. &          ! too close to surface: FinROSE
     & metdat(ind_tempr) > 240.) then            ! high temperature

      !
      ! Conditions are not suitable for PSC. Destroy the whole thing
      ! - H2SO4 -> to sulphates of the same size
      ! - Other species are evaporated to gas phase from NAT, STS, ICE
      ! - number cnc of PSC aerosols is zeroed
      ! - reaction rates on PSC surface are zeroed
      ! - since STS H2SO4 mass is linked to fine sulphates should not touch it
      ! - HNO3 and H2O exist in STS only within PSC sub + SAD calculations. All other time they 
      !   are in gas phase
      !
      ifPSC = .false.
      ifICE = .false.
      ifNAT = .false.
      ifSTS = .false.
      !wtH2SO4_prc = 0.11  ! mass fraction, %, anything that is not too stupid
      wtH2SO4_prc = 60.0   ! mass fraction, %, anything that is not too stupid. This value is now more cosistent with density of 1500 kg/m3
      density = 1500.
      vCncAer(iSO4c_aer) = vCncAer(iSO4c_aer) + vCncAer(indMAAD_aerosol(iH2SO4_loc,iModeCoarse,binICE)) !+ &
!                                              & vCncAer(indMAAD_aerosol(iH2SO4_loc,iModeFine,binSTS))
      vCncTrn(iHNO3g) = vCncTrn(iHNO3g) + vCncAer(indMAAD_aerosol(iHNO3_loc,iModeCoarse,binICE)) + &
                                        & vCncAer(indMAAD_aerosol(iNAT_loc,iModeFine,binNAT)) + &
                                        & vCncAer(indMAAD_aerosol(iNAT_loc,iModeCoarse,binNAT))  ! + &
                                        !& vCncAer(indMAAD_aerosol(iHNO3_loc,iModeFine,binSTS))
      vCncTrn(iH2Og) = vCncTrn(iH2Og) + vCncAer(indMAAD_aerosol(iH2O_loc,iModeCoarse,binICE)) + &
                                      & 3. * (vCncAer(indMAAD_aerosol(iNAT_loc,iModeFine,binNAT)) + &
                                            & vCncAer(indMAAD_aerosol(iNAT_loc,iModeCoarse,binNAT))) !+ &
                                            !& vCncAer(indMAAD_aerosol(iH2O_loc,iModeFine,binSTS))
      vCncAer(indMAAD_aerosol(iH2SO4_loc,iModeCoarse,binICE)) = 0.
      vCncAer(indMAAD_aerosol(iHNO3_loc,iModeCoarse,binICE)) = 0.
      vCncAer(indMAAD_aerosol(iH2O_loc,iModeCoarse,binICE)) = 0.
      vCncAer(indMAAD_aerosol(iNAT_loc,iModeFine,binNAT)) = 0.
      vCncAer(indMAAD_aerosol(iNAT_loc,iModeCoarse,binNAT)) = 0.
!      vCncAer(indMAAD_aerosol(iH2SO4_loc,iModeFine,binSTS)) = 0.
!      vCncAer(indMAAD_aerosol(iHNO3_loc,iModeFine,binSTS)) = 0.
!      vCncAer(indMAAD_aerosol(iH2O_loc,iModeFine,binSTS)) = 0.
      
      if(rules%iSpectrumDefinitionSwitch == dynamic_spectra)then
        vCncTrn(indMAAD_transp(iModeFine,binNAT)) = 0.
        vCncTrn(indMAAD_transp(iModeCoarse,binNAT)) = 0.
        vCncTrn(indMAAD_transp(iModeCoarse,binICE)) = 0.
      endif
      rulesRates%h_clono2_h2o = 0.
      rulesRates%h_brono2_h2o = 0. !3Aug2017 by R.H. enabled the reaction (2)
      rulesRates%h_n2o5_h2o = 0.
      rulesRates%h_clono2_hcl = 0.
      rulesRates%h_hocl_hcl = 0.
      rulesRates%h_brono2_hcl = 0. !3Aug2017 by R.H. enabled the reaction (6)
      rulesRates%h_hobr_hcl = 0.   !3Aug2017 by R.H. enabled the reaction (7)
      rulesRates%h_n2o5_hcl = 0.
      !3Aug2017 by R.H. enabled the reactions (9-13)
      rulesRates%h_clono2_hbr = 0.
      rulesRates%h_hocl_hbr = 0.
      rulesRates%h_brono2_hbr = 0.
      rulesRates%h_hobr_hbr = 0.
      rulesRates%h_n2o5_hbr = 0.

call check_mass_vector(vCncTrn, garbage, species, 'No-PSC conditions', fLowCncThresh, nSpeciesTrn, ix, iy, iz, print_it)

    else
      ifPSC = .true.
    endif

    ! Organize mass concentrations of aerosol compounds into a matrix of volumes.
    ! Also some checking here, that the volume-number ratio is reasonable, as no one outside 
    ! this module should be able to change the particle size to be outside the bin.
    !
    masses = 0.
    core = 0.
    naero = 0.

    select case(rules%iSpectrumDefinitionSwitch)

      case (prescribe_sizes, prescribe_number_concentr)
        !
        ! Need to calculate naero from mass concentrations. ALso fill-in masses and core
        !
        do iBin = 1,nBins 
          do iS = 1, nParallelBins 
            do iMat = 1, nAerMAAD
              if(indMAAD_aerosol(iMat,iBin,iS) /= int_missing)then
                if(vCncAer(indMAAD_aerosol(iMat,iBin,iS)) > &
                                                 & fLowCncThresh(indMAAD_aerosol(iMat,iBin,iS)))then
                  masses(iMat,iBin,iS) = vCncAer(indMAAD_aerosol(iMat,iBin,iS))
#ifdef DEBUG_MORE
if(print_it)then
call msg('Valuable, iBinType, iMat:',iS,iMat)
endif
#endif
                else
                  !
                  ! Send to garbage: transported species have correct indices there
                  !
#ifdef DEBUG_MORE
if(print_it)then
call msg('Garbage, iBinType, iMat:',iS,iMat)
endif
#endif
                  garbage(indMAAD_aerosol(iMat,iBin,iS)) = &
                                          & garbage(indMAAD_aerosol(iMat,iBin,iS)) + &
                                          & vCncAer(indMAAD_aerosol(iMat,iBin,iS))
                  vCncAer(indMAAD_aerosol(iMat,iBin,iS)) = 0.0
                endif  ! low mass threshold
              endif  ! material-bin-type combination exists
            enddo   ! materials
            !
            ! Having masses, get the numbers or sizes
            !
            if(rules%iSpectrumDefinitionSwitch == prescribe_sizes)then
              naero(iBin,iS) = 0.
              do iMat = 1, nAerMAAD
                if (indMAAD_aerosol(iMat,iBin,iS) /= int_missing) then
                   naero(iBin,iS) = naero(iBin,iS) + masses(iMat,iBin,iS) * &
                      & rules%mass2vol(indMAAD_aerosol(iMat,iBin,iS)) / rules%bin_particle_volume(iBin)
                   !RISTO TEST 27OCT2017:
                   !write(*,*) 'For iBin,iS,iMat = ',iBin,iS,iMat,', mass2vol = ', rules%mass2vol(indMAAD_aerosol(iMat,iBin,iS))
                   !END RISTO TEST 27OCT2017:
                end if
             end do
#ifdef DEBUG_MORE
if(print_it)then
call msg('iBin,iType,naero:' + fu_str(iBin) + ',' + fu_str(is), naero(iBin,iS))
endif
#endif
            elseif(rules%iSpectrumDefinitionSwitch == prescribe_number_concentr)then
              if(iBin == 1 .and. iS == binSTS)then  ! STS has tricky number concentration profile
                if (metdat(ind_pres) < 500.)then
                  naero(iBin,iS) = 0.1e6     ! #/m3
                elseif (metdat(ind_pres) > 20000.) then
                  naero(iBin,iS) = 10.e6     ! #/m3
                else
                  naero(iBin,iS) = 1.e6 * amax1(0.1, (((4.6615e-16 * metdat(ind_pres) - 2.3926e-11) * &   ! #/m3
                                                     & metdat(ind_pres) + 3.2760e-7) * &
                                                    & metdat(ind_pres) - 1.1182e-4) * &
                                                   & metdat(ind_pres) - 3.7668e-2)
                endif
              else
                naero(iBin,iS) = rules%fixed_number_concentration(iBin,iS)
              endif  ! if fine STS bin
            else
              call set_error('Unknown spectrum definition:'+fu_str(rules%iSpectrumDefinitionSwitch), &
                           & 'transform_AerDynMidAtm')
              return
            endif  ! id prescribe sizes or number concentrations
          end do  ! iS
        end do  ! iBin
      
      case(dynamic_spectra)
        !
        ! Making naero, masses, and core. Need to check the correct parameters of the districbutions
        ! DOES NOT WORK YET
        !
        do iBin = 1,nBins 
!call msg('d', fu_mean_D(modesMAAD(iBin, 1)))
!call msg('d', fu_mean_D(modesMAAD(iBin, 2)))
!call msg('log_d_bin', log10(fu_max_D(modesMAAD(iBin,1))*1e6) - log10(fu_min_D(modesMAAD(iBin,1))*1e6))
          do iS = 1, nParallelBins 
            if(indMAAD_transp(iBin,iS) == int_missing)cycle
            if(vCncTrn(indMAAD_transp(iBin,iS)) > fLowCncThresh(indMAAD_transp(iBin,iS)))then
              !
              ! Reasonable number cnc. Mass ?
              !
              naero(iBin,iS) = vCncTrn(indMAAD_transp(iBin,iS))
              do iMat = 1, nAerMAAD
                if(indMAAD_aerosol(iMat,iBin,iS) /= int_missing)then
                  if(vCncAer(indMAAD_aerosol(iMat,iBin,iS)) > fLowCncThresh(indMAAD_aerosol(iMat,iBin,iS)))then
                    ! 
                    masses(iMat,iBin,iS) =  vCncAer(indMAAD_aerosol(iMat,iBin,iS)) !* &
                                                      !& rules%mass2vol(indMAAD_aerosol(iMat,iBin,iS))
                  else
                    ! Send to garbage, where the aerosol species are after transport ones
#ifdef DEBUG_MORE
if(print_it)then
call msg('Garbage, iBin, vMass:',iBin,vCncAer(indMAAD_aerosol(iMat,iBin,iS)))
endif
#endif
                    garbage(indMAAD_aerosol(iMat,iBin,iS) + nSpeciesTrn) = &
                                          & garbage(indMAAD_aerosol(iMat,iBin,iS) + nSpeciesTrn) + &
                                          & vCncAer(indMAAD_aerosol(iMat,iBin,iS))
                    vCncAer(indMAAD_aerosol(iMat,iBin,iS)) = 0.0
                  endif
                endif  ! material-bin-type combination exists
              enddo   ! materials
              !
              ! What if all masses are sent to garbage? Have to kill n-concentrations too
              !
              core(iBin,iS) = 0.
              if(sum(masses(1:nAerMAAD,iBin,iS)) == 0.0)then
                call msg_warning('reasonable n-cnc but zero volume-cnc: (iBin,iParalelBin): (' + &
                               & fu_str(iBin) + ',' + fu_str(iS) + ')', &
                               & 'transform_AerDynMidAtm')
                garbage(indMAAD_transp(iBin,iS)) = garbage(indMAAD_transp(iBin,iS)) + &
                                                 & vCncTrn(indMAAD_transp(iBin,iS))
                naero(iBin,iS) = 0.
              else
                !
                ! Mean volume of a single particle in the bin
                !
                do iMat = 1, nAerMAAD
                  core(iBin,iS) = core(iBin,iS) + masses(iMat,iBin,iS) * &
                                                & rules%mass2vol(indMAAD_aerosol(iMat,iBin,iS))
                enddo
                core(iBin,iS) = core(iBin,iS) / naero(iBin, iS)
              endif
              !
              ! Check for its diameter
              !
              d =(6. * core(iBin,iS) / pi)**(1./3.)
              if(d < fu_min_D(modesMAAD(iBin,iS)) .or. d > fu_max_D(modesMAAD(iBin,iS)))then
                call set_error('Particle diameter at the beginning is outside bin limits', &
                             & 'transform_AerDynMidAtm')
                call msg('iBin, iS: ', iBin, iS)
                call msg('d, naero', d, real(naero(iBin, iS)))
                call msg('core, masses:',(/core(iBin, iS),masses(1:nAerMAAD,iBin,iS)/))
                call msg('min, max', fu_min_D(modesMAAD(iBin,iS)), fu_max_D(modesMAAD(iBin,iS)))
              endif
            else
              !
              ! Too low number cnc. Send to garbage but check for mass to be sure in consistency
              !
              do iMat = 1, nAerMAAD
                if(indMAAD_aerosol(iMat,iBin,iS) /= int_missing)then
                  !
                  ! If mass is sufficient, there has to be an error but since the mass and number 
                  ! thresholds are largely independent, can still proceed. Just make some noise.
                  !
                  if(vCncAer(indMAAD_aerosol(iMat,iBin,iS)) > &
                                                 & fLowCncThresh(indMAAD_aerosol(iMat,iBin,iS)))then
                    call set_error('Low number but enough mass cnc at the beginning', 'transform_AerDynMidAtm')
                    call msg('iBin, nbr', iBin, vCncTrn(indMAAD_transp(iBin,iS)))
                    call msg('Material, iS, mass:', iMat, vCncAer(indMAAD_aerosol(iMat,iBin,iS)))
                    call msg('Mass threshold:', fLowCncThresh(indMAAD_aerosol(iMat,iBin,iS)))
                    call unset_error('transform_AerDynMidAtm')
                  endif
                  garbage(indMAAD_aerosol(iMat,iBin,iS) + nSpeciesTrn) = &
                                               & garbage(indMAAD_aerosol(iMat,iBin,iS) + nSpeciesTrn) + &
                                               & vCncAer(indMAAD_aerosol(iMat,iBin,iS))
                endif   ! material-bin-type combination exists
              enddo  ! materials
              garbage(indMAAD_transp(iBin,iS)) = garbage(indMAAD_transp(iBin,iS)) + &
                                               & vCncTrn(indMAAD_transp(iBin,iS))
              vCnCTrn(indMAAD_transp(iBin,iS)) = 0.

            endif  ! sufficient number concentration
          enddo   ! parallel bins
        enddo   ! bins

      case default
        call msg('Funny iSpectrumDefinitionSwitch:',rules%iSpectrumDefinitionSwitch)
        call set_error('Funny rules%iSpectrumDefinitionSwitch','transform_AerDynMidAtm')
        return
    end select  ! spectrum description

    if(error)return

call check_mass_vector(vCncTrn, garbage, species, 'Before PSC', fLowCncThresh, nSpeciesTrn, ix, iy, iz, print_it)

    sadSTS = 0.
    sadNAT = 0.
    sadICE = 0.
    sadOTHER = 0.
    !
    !---------------------------------------------------------------------
    !
    ! Create/evolve the polar clouds if the conditions are suitable
    !
    if(ifPSC)then
      !
      ! Polar stratospheric clouds may exist. Calculate their surface area density, 
      ! mass, and number concentrations: contribution to the heterogeneous reaction rates
      !
#ifdef DEBUG_MORE
if(print_it)then
write(unit=strTmp,fmt=*)masses
call msg('Before polar_stratospheric_clouds: masses:' + strMetaMasses)
call msg('Before polar_stratospheric_clouds: masses:' + strTmp)
write(unit=strTmp,fmt=*)naero
call msg('naero:' + strMetaNaero)
call msg('naero:' + strTmp)
call msg('sad (sadICE, sadNAT, sadSTS, sadOTHER):',(/sadICE, sadNAT, sadSTS, sadOTHER/))
endif
#endif      


do iBin = 1,nBins 
  do iS = 1, nParallelBins 
    if(.not. naero(iBin,iS) >= 0.0)then
      call msg('Error in naero before PSC. naero(1:2_nBins, binSTS=1, binNAT=2, binICE=3_nParallelBins):',(/naero(1:nBins,1:nParallelBins)/))
      call msg('Masses:',(/masses(1:nAerMAAD, 1:nBins, 1:nParallelBins)/))
      call msg('Core:',(/core(1:nBins, 1:nParallelBins)/))
      call set_error('Strange naero before PSC','transform_AerDynMidAtm')
      return
    endif
  end do
end do

vCncTrn_iHNO3g = vCncTrn(iHNO3g)
vCncSL_iH2SO4_gas_sl = vCncSL(iH2SO4_gas_sl)
vCncTrn_iH2Og = vCncTrn(iH2Og)
vCncTrn_iHCl = vCncTrn(iHCl)
!3Aug2017 by R.H. Added HOCl and HBr
vCncTrn_iHOCl = vCncTrn(iHOCl)
vCncTrn_iHBr = vCncTrn(iHBr)
masses_tmp = masses
core_tmp = core
naero_tmp = naero
      if (error) return
      !RISTO Cl&Br: 3Aug2017 by R.H. Added HBr and HOCl concentrations after HCl?
      call polar_stratospheric_clouds(vCncTrn(iHNO3g), vCncSL(iH2SO4_gas_sl), &
                                    !& vCncTrn(iH2Og), vCncTrn(iHCl), & ! gases
                                    & vCncTrn(iH2Og), vCncTrn(iHCl), vCncTrn(iHOCl), vCncTrn(iHBr), & ! gases
                                    & masses, core, naero, density, &      ! aerosols
                                    & ifICE, ifNAT, ifSTS, &  ! what cloud types are present?
                                    & wtH2SO4_prc, &          ! STS sulphur mass fraction (%)
                                    & extra_sulphur_mass_STS, extra_particles_STS, & 
                                    & rules, rulesRates, metdat, fLowCncThresh, seconds, &  ! metadata
                                    & print_it)
      if(error) then
           !$OMP CRITICAL(after_psc_report)               
              call unset_error("Right after PSC")
              call msg("Thread No"&
                !$  &, omp_get_thread_num() &
              &)
              call msg('****Error right after after PSC. naero(1:2_nBins, binSTS=1, binNAT=2, binICE=3_nParallelBins):',(/naero(1:nBins,1:nParallelBins)/))
              call msg('****Masses:',(/masses(1:nAerMAAD, 1:nBins, 1:nParallelBins)/))
              call msg('****Core:',(/core(1:nBins, 1:nParallelBins)/))
              call msg('**Before the PSC call they were:')
              !27July2017 by R.H.: changed naero, masses, and core -> naero_tmp, masses_tmp, and core_tmp (3 lines) in order to show the old values 
              call msg('**naero(1:2_nBins, binSTS=1, binNAT=2, binICE=3_nParallelBins):',(/naero_tmp(1:nBins,1:nParallelBins)/))
              call msg('**Masses:',(/masses_tmp(1:nAerMAAD, 1:nBins, 1:nParallelBins)/))
              call msg('**Core:',(/core_tmp(1:nBins, 1:nParallelBins)/))
              call msg('****Repeating the PSC call verbosely')
              print_it = .true.
              !RISTO Cl&Br: 3Aug2017 by R.H. Added HBr and HOCl concentrations after HCl?
              call polar_stratospheric_clouds(vCncTrn_iHNO3g, vCncSL_iH2SO4_gas_sl, &
                                            !& vCncTrn_iH2Og, vCncTrn_iHCl, & ! gases
                                            & vCncTrn_iH2Og, vCncTrn_iHCl, vCncTrn_iHOCl, vCncTrn_iHBr, & ! gases
                                            & masses_tmp, core_tmp, naero_tmp, density, &      ! aerosols
                                            & ifICE, ifNAT, ifSTS, &  ! what cloud types are present?
                                            & wtH2SO4_prc, &          ! STS sulphur mass fraction (%)
                                            & extra_sulphur_mass_STS, extra_particles_STS, & 
                                            & rules, rulesRates, metdat, fLowCncThresh, seconds, &  ! metadata
                                            & print_it)
              !27July2017 by R.H.: changed naero, masses, and core -> naero_tmp, masses_tmp, and core_tmp (2 lines) in order to show the re-calculated values 
              call msg('****Masses:',(/masses_tmp(1:nAerMAAD, 1:nBins, 1:nParallelBins)/))
              call msg('****Core:',(/core_tmp(1:nBins, 1:nParallelBins)/))
              call set_error('****Setting error back after PSC','transform_AerDynMidAtm')
          !$OMP END CRITICAL(after_psc_report)               
              return
      endif
      
call check_mass_vector(vCncTrn, garbage, species, 'After PSC call', fLowCncThresh, nSpeciesTrn, ix, iy, iz, print_it)

do iBin = 1,nBins 
  do iS = 1, nParallelBins 
    if(.not. naero(iBin,iS) >= 0.0)then
      call msg('****Error in naero after PSC. naero(1:2_nBins, binSTS=1, binNAT=2, binICE=3_nParallelBins):',(/naero(1:nBins,1:nParallelBins)/))
      call msg('****Masses:',(/masses(1:nAerMAAD, 1:nBins, 1:nParallelBins)/))
      call msg('****Core:',(/core(1:nBins, 1:nParallelBins)/))
      call msg('**Before the PSC call they were:')
      !27July2017 by R.H.: changed naero, masses, and core -> naero_tmp, masses_tmp, and core_tmp (3 lines) in order to show the old values 
      call msg('**naero(1:2_nBins, binSTS=1, binNAT=2, binICE=3_nParallelBins):',(/naero_tmp(1:nBins,1:nParallelBins)/))
      call msg('**Masses:',(/masses_tmp(1:nAerMAAD, 1:nBins, 1:nParallelBins)/))
      call msg('**Core:',(/core_tmp(1:nBins, 1:nParallelBins)/))
      call msg('****Repeating the PSC call verbosely')
      print_it = .true.
      !RISTO Cl&Br: 3Aug2017 by R.H. Added HBr and HOCl concentrations after HCl?
      call polar_stratospheric_clouds(vCncTrn_iHNO3g, vCncSL_iH2SO4_gas_sl, &
                                    !& vCncTrn_iH2Og, vCncTrn_iHCl, & ! gases
                                    & vCncTrn_iH2Og, vCncTrn_iHCl, vCncTrn_iHOCl, vCncTrn_iHBr,& ! gases
                                    & masses_tmp, core_tmp, naero_tmp, density, &      ! aerosols
                                    & ifICE, ifNAT, ifSTS, &  ! what cloud types are present?
                                    & wtH2SO4_prc, &          ! STS sulphur mass fraction (%)
                                    & extra_sulphur_mass_STS, extra_particles_STS, & 
                                    & rules, rulesRates, metdat, fLowCncThresh, seconds, &  ! metadata
                                    & print_it)
      !27July2017 by R.H.: changed naero, masses, and core -> naero_tmp, masses_tmp, and core_tmp (2 lines) in order to show the re-calculated values 
      call msg('****Masses:',(/masses_tmp(1:nAerMAAD, 1:nBins, 1:nParallelBins)/))
      call msg('****Core:',(/core_tmp(1:nBins, 1:nParallelBins)/))
      call set_error('****Strange naero after PSC','transform_AerDynMidAtm')
      return
    endif
  end do
end do


if(print_it)then
call msg('ix, iy, iz',(/ix,iy,iz/))
call msg('TRANSFORM, wtH2SO4_prc',wtH2SO4_prc)      
write(unit=strTmp,fmt=*)masses
call msg('After polar_stratospheric_clouds: masses:' + strMetaMasses)
call msg('After polar_stratospheric_clouds: masses:' + strTmp)
write(unit=strTmp,fmt=*)naero
call msg('naero:' + strMetaNaero)
call msg('naero:' + strTmp)
call msg('sad (/sadICE, sadNAT, sadSTS, sadOTHER/):',(/sadICE, sadNAT, sadSTS, sadOTHER/))
endif
      !
      ! Having the PSCs done, calculate their surface area, [m2/m3].
      ! In terms of volume, it is: (Pi * N * 6**2 * v**2) ** (1./3.)
      ! here N is number concentration #/m3, v - volume, m3/m3
      !
      if(ifICE)then                                                                          ! ICE
        do iBin = 1, nBins
          if(naero(iBin,binICE) < 1e-5)cycle
          fTmp = (masses(iH2SO4_loc,iBin,binICE) * rules%mwH2SO4 + &
                & masses(iHNO3_loc,iBin,binICE) * rules%mwHNO3 + &
                & masses(iH2O_loc,nBins,binICE) * rules%mwH2O) / density(binICE) 
          sadICE = sadICE + (Pi_36 * naero(iBin,binICE) * fTmp * fTmp) ** 0.333333333333333333
        end do
      else
        if(ifNAT)then                                                                        ! NAT
          do iBin = 1, nBins
            if(naero(iBin,binNAT) < 1e-5)cycle
!            fTmp = (masses(iH2SO4_loc,iBin,binNAT) * rules%mwH2SO4 + &
!                  & masses(iNAT_loc,iBin,binNAT) * rules%mwNAT) / density(binNAT)
            fTmp = masses(iNAT_loc,iBin,binNAT) * rules%mwNAT / density(binNAT)  ! NAT is only NAT
            sadNAT = sadNAT + (Pi_36 * naero(iBin,binNAT) * fTmp * fTmp) ** 0.333333333333333333
if(.not. sadNAT >= 0.)then
  call msg('sadNAT is not non-negative:',sadNAT)
  call msg('(Pi_36 * naero(iBin,binNAT) * fTmp * fTmp)**0.333. Pi_36,naero,fTmp:', (/Pi_36, naero(iBin,binNAT), fTmp/))
  call msg('fTmp = masses(iNAT_loc,iBin,binNAT) * rules%mwNAT / density(binNAT)', &
                                    & (/masses(iNAT_loc,iBin,binNAT), rules%mwNAT, density(binNAT)/))
  call msg('naero(1:2_nBins, binSTS=1, binNAT=2, binICE=3_nTypes)',(/naero(1:nBins,binSTS), naero(1:nBins,binNAT), naero(1:nBins,binICE)/))
  call msg('masses(iSubstLoc, 1:2_nBins, types_as_above)', (/masses(1:nAerMAAD, 1:nBins, 1:nParallelBins)/))
  call set_error('sadNAT is not positive','transform_AerDynMidAtm')
  return
endif
          end do
        endif   ! ifNAT

        do iBin = 1, nBins                                                                   ! STS
          if(naero(iBin,binSTS) < 1e-5)cycle
          !RISTO Cl&Br: 3Aug2017 by R.H. added the HCl, HOCl, and HBr concentrations.
          !However, note that the density of STS is assumed to be unaffected by these.           
          fTmp = (masses(iH2SO4_loc,iBin,binSTS) * rules%mwH2SO4 + &
                & masses(iHNO3_loc,iBin,binSTS) * rules%mwHNO3 + &
                & masses(iH2O_loc,iBin,binSTS) * rules%mwH2O + &
                & masses(iHCl_loc,iBin,binSTS) * rules%mwHCl + &
                & masses(iHOCl_loc,iBin,binSTS) * rules%mwHOCl + &
                & masses(iHBr_loc,iBin,binSTS) * rules%mwHBr)/ density(binSTS)
          sadSTS = sadSTS + (Pi_36 * naero(iBin,binSTS) * fTmp * fTmp) ** 0.33333333333333333
        end do
        !RISTO TEST 7Sep2017
        !if (.not. (sadSTS > 0.0)) then
        !   call msg('sadSTS not > 0:')
        if (sadSTS < 0.0) then
           call msg('sadSTS < 0:')
           call msg('molar weights: ', (/rules%mwH2SO4,rules%mwHNO3,rules%mwH2O,rules%mwHCl,rules%mwHOCl,rules%mwHBr/) )
           call msg('masses:', (/masses(iH2SO4_loc,1,binSTS), masses(iHNO3_loc,1,binSTS), masses(iH2O_loc,1,binSTS), &
                               & masses(iHCl_loc,1,binSTS),   masses(iHOCl_loc,1,binSTS), masses(iHBr_loc,1,binSTS) /) )
           call msg('masses:', (/masses(iH2SO4_loc,2,binSTS), masses(iHNO3_loc,2,binSTS), masses(iH2O_loc,2,binSTS), &
                               & masses(iHCl_loc,2,binSTS),   masses(iHOCl_loc,2,binSTS), masses(iHBr_loc,2,binSTS) /) )
           call msg('density: ',density(binSTS) )
        end if
        !END RISTO TEST 7Sep2017
      endif  ! ifICE
    else
       extra_sulphur_mass_STS = 0.
       extra_particles_STS = 0.
    endif  ! if PSC

    !
    ! Get the surface area for the non-PSC particles.
    ! For the time being, they are considered entirely outside the story - just usual transported 
    ! species with non-missing aerosol features and reasonable particle diameter.
    !
    do iTmp = 1, nSpeciesTrn
      if(rules%mean_diam_foreign_aerosol(iTmp) > 0.)then
        sadOTHER = sadOTHER + &
                & (6. * vCncTrn(iTmp) * rules%mass2vol(iTmp) / rules%mean_diam_foreign_aerosol(iTmp))
        if (.not. sadOTHER >= 0.) then 
           call msg("Gotcha!"+fu_str(iTmp),(/vCncTrn(iTmp), rules%mass2vol(iTmp), rules%mean_diam_foreign_aerosol(iTmp)/))
        endif
      endif
    end do
    if (.not. sadOTHER >= 0.) then
            call msg("Strange sadOTHER nSpeciesTrans:",nSpeciesTrn)
            call msg( "vCncTrn", vCncTrn(1:nSpeciesTrn))
            call msg( "rules%mass2vol", rules%mass2vol(1:nSpeciesTrn))
            call msg("rules%mean_diam_foreign_aerosol", rules%mean_diam_foreign_aerosol(1:nSpeciesTrn))
    endif


#ifdef DEBUG_MORE
if(print_it)then
write(unit=strTmp,fmt=*)masses
call msg('Before rates: masses:' + strMetamasses)
call msg('Before rates: masses:' + strTmp)
write(unit=strTmp,fmt=*)naero
call msg('naero:' + strMetaNaero)
call msg('naero:' + strTmp)
call msg('sad (/sadICE, sadNAT, sadSTS, sadOTHER/):',(/sadICE, sadNAT, sadSTS, sadOTHER/))
endif
#endif
    !-----------------------------------------------------------------------------------------
    !
    ! heterogeneous reaction rates.
    ! Note that GAMMA gives them as pseudo-first-order rates, [sec-1]
    !
!RISTO Cl&Br: Check that all the required substances are added! Should be ok. Only reactant species are needed, not result species.
    call rates_pseudo_first_order(ifICE, ifNAT, ifSTS, &
                                & seconds, &  ! time step
                                & metdat(ind_tempr), & ! temperature (pressure not needed any more)
                                !& metdat(ind_tempr), metdat(ind_pres), & ! temperature and pressure
!                                & vCncTrn(iH2Og), vCncTrn(iHCl), vCncTrn(iHOCl), vCncTrn(iClONO2), & !gases
                                !3Aug2017 by R.H. added the additional gasses required for Cl and Br reactions
                                & vCncTrn(iH2Og), vCncTrn(iHCl), vCncTrn(iHOCl), vCncTrn(iClONO2), & !gases
                                & vCncTrn(iN2O5), vCncTrn(iHBr), vCncTrn(iHOBr), vCncTrn(iBrONO2), & !gases
                                & vCncTrn(iHNO3g), & !gases (this is needed for some reaction probabilities)
                                & sadSTS,sadNAT,sadICE,sadOTHER, &     ! surface area densities, m2/m3
!                                & sadSTS*2,sadNAT*2,sadICE*2,sadOTHER*2, &     ! surface ara densities, m2/m3
                                & (3.*rules%bin_particle_volume(iModeFine)/(4.*Pi))**(1./3.), & ! STS aerosol radius [m]
                                & wtH2SO4_prc, density(binSTS), &  ! STS sulphur mass fraction (%), aerosol density
                                & rules, rulesRates, &        ! other metadata
                                & print_it)
    if(error)return
    
call check_mass_vector(vCncTrn, garbage, species, 'After rates call', fLowCncThresh, nSpeciesTrn, ix, iy, iz, print_it)

#ifdef DEBUG_MORE
if(print_it)then
write(unit=strTmp,fmt=*)masses
call msg('After rates: messes:' + strTmp)
write(unit=strTmp,fmt=*)naero
call msg('naero:' + strTmp)
call msg('sadICE, sadNAT, sadSTS, sadOTHER:',(/sadICE, sadNAT, sadSTS, sadOTHER/))
if(naero(2,2) > 18.81)then
  call msg('Here')
endif
endif
#endif
    !----------------------------------------------------------------------------------------------
    !
    ! Thermodynamics is over, we need to leave only transported amounts in masses and naero variables
    ! before returning them back. This is somewhat stupid place to put this return but we cannot do it
    ! prior to computing the surface area.
    !
    if(masses(iH2SO4_loc,1,binSTS) >= 0.999 * extra_sulphur_mass_STS)then
      masses(iH2SO4_loc,1,binSTS) = max(0.,masses(iH2SO4_loc,1,binSTS) - extra_sulphur_mass_STS)
      naero(1,binSTS) = max(0., naero(1,binSTS) - extra_particles_STS)
    else
      call msg('Funny H2SO4 mass: less than its artificial part:', &
             & masses(iH2SO4_loc,1,binSTS), extra_sulphur_mass_STS)
      call set_error('Funny H2SO4 volume: less than its artificial part','transform_AerDynMidAtm')
      masses(iH2SO4_loc,1,binSTS) = 0.
    endif
    !
    ! Now return the aerosols: numbers to transported MASsMap, masses to aerosols MM
    !
    do iBin = 1,nBins 
!call msg('d', fu_mean_D(modesMAAD(iBin,1)))
!call msg('d', fu_mean_D(modesMAAD(iBin,1)))
!call msg('log_d_bin', log10(fu_max_D(modesMAAD(iBin,1))*1e6) - log10(fu_min_D(modesMAAD(iBin,1))*1e6))
      do iS = 1, nParallelBins
!        if(iS == binSTS)cycle
        if(indMAAD_transp(iBin,iS) /= int_missing) &
                                     & vCncTrn(indMAAD_transp(iBin,iS)) = naero(iBin,iS)    ! number
        do iMat = 1, nAerMAAD
          if(indMAAD_aerosol(iMat,iBin,iS) /= int_missing)then
#ifdef DEBUG_MORE
if(print_it)then
call msg('Returning the species:' + fu_str(indMAAD_aerosol(iMat,iBin,iS)), vCncAer(indMAAD_aerosol(iMat,iBin,iS)), &
                                                  & masses(iMat,iBin,iS))
endif
#endif
            vCncAer(indMAAD_aerosol(iMat,iBin,iS)) = masses(iMat,iBin,iS) !/  &
                                           !& rules%mass2vol(indMAAD_aerosol(iMat,iBin,iS))      ! mass
          endif  ! material-bin-type combination exists
        enddo   ! materials
      enddo   ! parallel bins
    enddo   ! bins
    !
    ! Additionally, two masses with no correspondence to transport species have to be "returned"
    ! to gas phase: HNO3, H2O, HCl, HOCl and HBr from STS. When needed, they will be picked up again to STS by TER.
    !
    vCncTrn(iHNO3g) = vCncTrn(iHNO3g) + masses(iHNO3_loc,iModeFine,binSTS)
    vCncTrn(iH2Og) = vCncTrn(iH2Og) + masses(iH2O_loc,iModeFine,binSTS)
    !3Aug2017 by R.H. Added the HCl, HOCl, and HBr
    vCncTrn(iHCl) = vCncTrn(iHCl) + masses(iHCl_loc,iModeFine,binSTS)
    vCncTrn(iHOCl) = vCncTrn(iHOCl) + masses(iHOCl_loc,iModeFine,binSTS)
    vCncTrn(iHBr) = vCncTrn(iHBr) + masses(iHBr_loc,iModeFine,binSTS)
    
    
!!$    !BEGIN RISTO TEST
!!$    if (ix==90 .and. iy==44) then
!!$       open(777, FILE='cncX90Y44fine.dat',ACCESS = 'APPEND')
!!$       write(777,*) seconds,iz,wtH2SO4_prc, &
!!$                  & masses(iH2O_loc,iModeFine,binSTS),masses(iH2O_loc,iModeFine,binICE),masses(iH2O_loc,iModeFine,binNAT), &
!!$                  & masses(iH2SO4_loc,iModeFine,binSTS),masses(iH2SO4_loc,iModeFine,binICE),masses(iH2SO4_loc,iModeFine,binNAT), &
!!$                  & masses(iHNO3_loc,iModeFine,binSTS),masses(iHNO3_loc,iModeFine,binICE),masses(iHNO3_loc,iModeFine,binNAT), &
!!$                  & masses(iNAT_loc,iModeFine,binSTS),masses(iNAT_loc,iModeFine,binICE),masses(iNAT_loc,iModeFine,binNAT)
!!$       close(777)
!!$    end if
!!$    if (ix==90 .and. iy==66) then
!!$       open(777, FILE='cncX90Y66fine.dat',ACCESS = 'APPEND')
!!$       write(777,*) seconds,iz,wtH2SO4_prc, &
!!$                  & masses(iH2O_loc,iModeFine,binSTS),masses(iH2O_loc,iModeFine,binICE),masses(iH2O_loc,iModeFine,binNAT), &
!!$                  & masses(iH2SO4_loc,iModeFine,binSTS),masses(iH2SO4_loc,iModeFine,binICE),masses(iH2SO4_loc,iModeFine,binNAT), &
!!$                  & masses(iHNO3_loc,iModeFine,binSTS),masses(iHNO3_loc,iModeFine,binICE),masses(iHNO3_loc,iModeFine,binNAT), &
!!$                  & masses(iNAT_loc,iModeFine,binSTS),masses(iNAT_loc,iModeFine,binICE),masses(iNAT_loc,iModeFine,binNAT)
!!$       close(777)
!!$    end if
!!$    !END RISTO TEST
    
!if(.not. ifPSC) return

!return


    !--------------------------------------------------------------------------------------------
    !
    ! Having the reaction rates known and stored in rules, proceed with the reactions.
    ! The pseudo-first-order coefs mean: dM/dt = -h M, i.e. M=M0*exp(-ht)
    !
    ! ATTENTION.
    ! In the below implementation, I assumed that there is always enough input mass,
    ! i.e. all reactions are independent. However, there may be competition. Then the
    ! reactions are to be grouped together producing the demand, which is then satisfied
    ! proportionally to the reaction rates.
    ! The second problem is over-consumption of the assumed-abudant agent. At least blocking
    ! it is mandatory
    !
!$OMP CRITICAL(Activation)  
  
    call check_mass_vector(vCncTrn, garbage, species, 'Before Cl activation', fLowCncThresh, nSpeciesTrn, ix, iy, iz, print_it)

!!$    !RISTO TEST: disabling the heterogeneous reactions 2017Dec:
!!$    rulesRates%h_clono2_h2o = 0.0 !Reaction (1)
!!$    rulesRates%h_brono2_h2o = 0.0 !Reaction (2) 
!!$    rulesRates%h_n2o5_h2o   = 0.0 !Reaction (3) 
!!$    rulesRates%h_clono2_hcl = 0.0 !Reaction (4) 
!!$    rulesRates%h_hocl_hcl   = 0.0 !Reaction (5) 
!!$    rulesRates%h_brono2_hcl = 0.0 !Reaction (6) 
!!$    rulesRates%h_hobr_hcl   = 0.0 !Reaction (7) 
!!$    rulesRates%h_n2o5_hcl   = 0.0 !Reaction (8) 
!!$    rulesRates%h_clono2_hbr = 0.0 !Reaction (9) 
!!$    rulesRates%h_hocl_hbr   = 0.0 !Reaction (10) 
!!$    rulesRates%h_brono2_hbr = 0.0 !Reaction (11) 
!!$    rulesRates%h_hobr_hbr   = 0.0 !Reaction (12)
!!$    rulesRates%h_n2o5_hbr   = 0.0 !Reaction (13) 
!!$    !END RISTO TEST
!!$    write(*,*) rulesRates%h_clono2_h2o, rulesRates%h_brono2_h2o, rulesRates%h_n2o5_h2o, rulesRates%h_clono2_hcl, rulesRates%h_hocl_hcl, &
!!$             & rulesRates%h_brono2_hcl, rulesRates%h_hobr_hcl, rulesRates%h_n2o5_hcl, rulesRates%h_clono2_hbr, rulesRates%h_hocl_hbr,   &
!!$             & rulesRates%h_brono2_hbr, rulesRates%h_hobr_hbr, rulesRates%h_n2o5_hbr, '% Rates%h_something RISTO rate-test!'
    
    !RISTO Cl&Br:  Not only reactions 2, 6, 7, 9, and 10, but also reactions 11-13 are included!
    !( 1)  clono2 +  h2o -> hocl  +  hno3
    !fTmp = rulesRates%h_clono2_h2o * seconds
    !if(.not. fTmp > 1.e-3)then
    !  if(.not. fTmp >= 0.)call set_error('Problem with rate h_clono2_h2o:'+fu_str(fTmp),'transform_AerDynMidAtm')
    !  dM = min(vCncTrn(iH2Og), vCncTrn(iClONO2)) * fTmp
    !else
    !  dM = min(vCncTrn(iH2Og), vCncTrn(iClONO2)) * (1. - exp(-fTmp))
    !endif
    dM = min(vCncTrn(iH2Og), vCncTrn(iClONO2)) * (1. - exp(-rulesRates%h_clono2_h2o * seconds))
    vCncTrn(iClONO2) = vCncTrn(iClONO2) - dM
    vCncTrn(iH2Og)   = vCncTrn(iH2Og)   - dM
    vCncTrn(iHOCl)   = vCncTrn(iHOCl)   + dM
    vCncTrn(iHNO3g)  = vCncTrn(iHNO3g)  + dM
    
    !call check_mass_vector(vCncTrn, garbage, species, 'After h_clono2_h2o', fLowCncThresh, nSpeciesTrn, ix, iy, iz, print_it)
    !if(print_it)call msg('rules%h_clono2_h2o',rules%h_clono2_h2o)

    !3Aug2017 by R.H. enabled the reaction (2):
    !( 2)  brono2 +  h2o -> hobr  +  hno3                     
    dM = min(vCncTrn(iH2Og), vCncTrn(iBrONO2) * (1. - exp(-rulesRates%h_brono2_h2o * seconds)))
    vCncTrn(iBrONO2) = vCncTrn(iBrONO2) - dM
    vCncTrn(iH2Og)   = vCncTrn(iH2Og)   - dM
    vCncTrn(iHOBr)   = vCncTrn(iHOBr)   + dM
    vCncTrn(iHNO3g)  = vCncTrn(iHNO3g)  + dM

    !( 3)  n2o5   +  h2o ->        2 hno3                
    dM = min(vCncTrn(iH2Og), vCncTrn(iN2O5)) * (1. - exp(-rulesRates%h_n2o5_h2o * seconds))
    vCncTrn(iN2O5)  = vCncTrn(iN2O5)  - dM
    vCncTrn(iH2Og)  = vCncTrn(iH2Og)  - dM
    vCncTrn(iHNO3g) = vCncTrn(iHNO3g) + dM*2.
    !call check_mass_vector(vCncTrn, garbage, species, 'After h_n2o5_h2o', fLowCncThresh, nSpeciesTrn, ix, iy, iz, print_it)
    !if(print_it)call msg('rules%h_n2o5_h2o',rulesRates%h_n2o5_h2o)

    !( 4)  clono2 +  hcl -> cl2   +  hno3          
    dM = min(vCncTrn(iHCl), vCncTrn(iClONO2)) * (1. - exp(-rulesRates%h_clono2_hcl * seconds))
    vCncTrn(iClONO2) = vCncTrn(iClONO2) - dM
    vCncTrn(iHCl)    = vCncTrn(iHCl)    - dM
    vCncTrn(iCl2)    = vCncTrn(iCl2)    + dM
    vCncTrn(iHNO3g)  = vCncTrn(iHNO3g)  + dM
    !call check_mass_vector(vCncTrn, garbage, species, 'After h_clono2_h2o', fLowCncThresh, nSpeciesTrn, ix, iy, iz, print_it)
    !if(print_it) call msg('rules%h_clono2_hcl', rules%h_clono2_hcl)

    !( 5)  hocl   +  hcl -> cl2   +  h2o          
    dM = min(vCncTrn(iHCl), vCncTrn(iHOCl)) * (1. - exp(-rulesRates%h_hocl_hcl * seconds))
    vCncTrn(iHOCl) = vCncTrn(iHOCl) - dM
    vCncTrn(iHCl)  = vCncTrn(iHCl)  - dM
    vCncTrn(iCl2)  = vCncTrn(iCl2)  + dM
    vCncTrn(iH2Og) = vCncTrn(iH2Og) + dM !Corrected the missing water 1Dec2017 by R.H. 

    !call check_mass_vector(vCncTrn, garbage, species, 'After rules%h_hocl_hcl', fLowCncThresh, nSpeciesTrn, ix, iy, iz, print_it)
    !if(print_it)call msg('rules%h_hocl_hcl',rules%h_hocl_hcl)

    !3Aug2017 by R.H. enabled the reaction (6):
    !( 6)  brono2 +  hcl -> brcl  +  hno3
    dM = min(vCncTrn(iHCl), vCncTrn(iBrONO2)) * (1. - exp(-rulesRates%h_brono2_hcl * seconds))
    vCncTrn(iBrONO2) = vCncTrn(iBrONO2) - dM
    vCncTrn(iHCl)    = vCncTrn(iHCl)    - dM
    vCncTrn(iBrCl)   = vCncTrn(iBrCl)   + dM
    vCncTrn(iHNO3g)  = vCncTrn(iHNO3g)  + dM

    !3Aug2017 by R.H. enabled the reaction (7):
    !( 7)  hobr   +  hcl -> brcl  +  h2o 
    dM = min(vCncTrn(iHCl), vCncTrn(iHOBr)) * (1. - exp(-rulesRates%h_hobr_hcl * seconds))
    vCncTrn(iHOBr) = vCncTrn(iHOBr) - dM
    vCncTrn(iHCl)  = vCncTrn(iHCl)  - dM
    vCncTrn(iBrCl) = vCncTrn(iBrCl) + dM
    vCncTrn(iH2Og) = vCncTrn(iH2Og) + dM !Corrected the missing water 1Dec2017 by R.H. 

    !( 8)  n2o5   +  hcl -> clno2 +  hno3
    dM = min(vCncTrn(iHCl), vCncTrn(iN2O5)) * (1. - exp(-rulesRates%h_n2o5_hcl * seconds))
    vCncTrn(iN2O5)  = vCncTrn(iN2O5)  - dM
    vCncTrn(iHCl)   = vCncTrn(iHCl)   - dM
    vCncTrn(iClNO2) = vCncTrn(iClNO2) + dM
    vCncTrn(iHNO3g) = vCncTrn(iHNO3g) + dM

    !call check_mass_vector(vCncTrn, garbage, species, 'After rules%h_n2o5_hcl', fLowCncThresh, nSpeciesTrn, ix, iy, iz, print_it)
    !if(print_it)call msg('rules%h_n2o5_hcl',rulesRates%h_n2o5_hcl)
    
    !3Aug2017 by R.H. enabled the reactions (9-13):
    !( 9)  clono2 +  hbr -> brcl  +  hno3
    dM = min(vCncTrn(iHBr), vCncTrn(iClONO2)) * (1. - exp(-rulesRates%h_clono2_hbr * seconds))
    vCncTrn(iClONO2) = vCncTrn(iClONO2) - dM
    vCncTrn(iHBr)    = vCncTrn(iHBr)    - dM
    vCncTrn(iBrCl)   = vCncTrn(iBrCl)   + dM
    vCncTrn(iHNO3g)  = vCncTrn(iHNO3g)  + dM

    !(10)  hocl   +  hbr -> brcl  +  h2o
    dM = min(vCncTrn(iHBr), vCncTrn(iHOCl)) * (1. - exp(-rulesRates%h_hocl_hbr * seconds))
    vCncTrn(iHOCl) = vCncTrn(iHOCl) - dM
    vCncTrn(iHBr)  = vCncTrn(iHBr)  - dM
    vCncTrn(iBrCl) = vCncTrn(iBrCl) + dM
    vCncTrn(iH2Og) = vCncTrn(iH2Og) + dM !Corrected the missing water 1Dec2017 by R.H. 

    !(11)  brono2 +  hbr -> br2   +  hno3
    dM = min(vCncTrn(iHBr), vCncTrn(iBrONO2)) * (1. - exp(-rulesRates%h_brono2_hbr * seconds))
    vCncTrn(iBrONO2) = vCncTrn(iBrONO2) - dM
    vCncTrn(iHBr)    = vCncTrn(iHBr)    - dM
    vCncTrn(iBr2)    = vCncTrn(iBr2)    + dM
    vCncTrn(iHNO3g)  = vCncTrn(iHNO3g)  + dM

    !(12)  hobr   +  hbr -> br2   +  h2o  !RISTO NOTE: Final species should be Br2 and H2O (not HNO3)
    dM = min(vCncTrn(iHBr), vCncTrn(iHOBr)) * (1. - exp(-rulesRates%h_hobr_hbr * seconds))
    vCncTrn(iHOBr) = vCncTrn(iHOBr) - dM
    vCncTrn(iHBr)  = vCncTrn(iHBr)  - dM
    vCncTrn(iBr2)  = vCncTrn(iBr2)  + dM
    vCncTrn(iH2Og) = vCncTrn(iH2Og) + dM

    !(13)  n2o5   +  hbr -> brno2 +  hno3
    dM = min(vCncTrn(iHBr), vCncTrn(iN2O5)) * (1. - exp(-rulesRates%h_n2o5_hbr * seconds))
    vCncTrn(iN2O5)  = vCncTrn(iN2O5)  - dM
    vCncTrn(iHBr)   = vCncTrn(iHBr)   - dM
    vCncTrn(iBrNO2) = vCncTrn(iBrNO2) + dM
    vCncTrn(iHNO3g) = vCncTrn(iHNO3g) + dM

    !11Oct2017 by R.H. added a check after all heterogeneous reactions:  
    call check_mass_vector(vCncTrn, garbage, species, 'After rules%h_n2o5_hbr', fLowCncThresh, nSpeciesTrn, ix, iy, iz, print_it)
    if(print_it)call msg('rules%h_n2o5_hcl',rulesRates%h_n2o5_hcl)
    
!$OMP END CRITICAL(Activation)  

!!$    !RISTO TEST TEST TEST TEST 5Dec2017
!!$    call msg('fLowCncThres:', fLowCncThresh)
!!$    call msg('indeces:    iH2Og,  iHNO3g, iClONO2,    iHCl,   iHOCl,   iN2O5,    iCl2,  iClNO2, iBrONO2,   iBrCl,    iBr2,  iBrNO2,    iHBr,   iHOBr')
!!$    call msg('indeces:', (/iH2Og, iHNO3g, iClONO2, iHCl, iHOCl, iN2O5, iCl2, iClNO2, iBrONO2, iBrCl, iBr2, iBrNO2, iHBr, iHOBr/) )
!!$    call msg('indeces:', (/fLowCncThresh(iH2Og), fLowCncThresh(iHNO3g), fLowCncThresh(iClONO2), fLowCncThresh(iHCl), fLowCncThresh(iHOCl), fLowCncThresh(iN2O5), &
!!$         & fLowCncThresh(iCl2), fLowCncThresh(iClNO2), fLowCncThresh(iBrONO2), fLowCncThresh(iBrCl), fLowCncThresh(iBr2), fLowCncThresh(iBrNO2), fLowCncThresh(iHBr), &
!!$         & fLowCncThresh(iHOBr)/) )
!!$    STOP
!!$    !RISTO TEST TEST TEST TEST 5Dec2017
    
  end subroutine transform_AerDynMidAtm


  !************************************************************************************
  !************************************************************************************
  !
  ! Below subroutines are copied from FinROSE code and adapted to SILAM environment.
  !
  !************************************************************************************
  !************************************************************************************
  !3Aug2017 by R.H. Included also the HCl, HOCl, and HBr concentrations
  subroutine polar_stratospheric_clouds(cncHNO3, cncH2SO4, cncH2O, cncHCl, cncHOCl, cncHBr, &
                                      & masses, &   ! mass of aerosols in a bin
                                      & core, &   ! mean volume of a single particle in a bin
                                      & naero, &  ! numbers of aerosols in a bin
                                      & density, & ! density of aerosols of each type
                                      & ifICE, ifNAT, ifSTS, & ! what types of clouds are present
                                      & wtH2SO4, & ! weight fraction of H2SO4 (%)
                                      & extra_sulphur_mass_STS, & ! 
                                      & extra_particles_STS, & !
                                      & rules, &  ! rulesAerDynMidAtm
                                      & rulesRates, &  ! local variable: rates, conversion constants, etc
                                      & metdat, fLowCncThresh, & ! meteorology, low-mass threshold
                                      & seconds, &  ! time steps, seconds
                                      & print_it)
    !
    ! Calculates the PSC parameters and heterogeneous reaction rates.
    ! It gets:
    ! - transported species, which include a bunch of gases and number concentration 
    !   of NAT and ICE particles
    ! - short-lived species, from where H2SO4 is taken. #### MUST BE CALLED PRIOR TO SAD ####
    ! - aerosol species, which have no own moments but rather linked to some of the 
    !   transported species, e.g. to the corresponding number concentrations.
    ! It can:
    ! - generate/evaporate NAT, STS, and ICE particles
    ! - crudely describe growth of a small fraction of them, others disappearing
    ! - heterogeneous reactions at PSC surface, mainly chlorine/bromine activation
    ! Output:
    ! - modified values of transported, short-lived, and aerosol species.
    ! By-products:
    !- the reaction rate coefficients for the heterogeneous reactions
    !  in/on sulfate aerosols and polar stratospheric cloud (PSC) particles 
    !- surface area densities (sad) and the composition of sulfate aerosols and PSCs
    !- partitoning between gas phase and condensed phase
    !
    !  There are two basic ways of calcualting PSCs: FinROSE equilibrium and SILAM dynamic.
    !
    ! Within FinROSE algorithm, PSC composition is calculated based on thermodynamic
    ! equilibrium. The PSC types that are considered, in addition to 
    ! aerosols, are STS/TER (PSC 1b), NAT (PSC 1a) and ICE (PSC 2).
    ! A parametrisation for NAT-rocks is also included. The NAT rock
    ! scheme is based on the use of a transported tracer. No activation
    ! barrier is assumed for formation of TER/STS from aerosols. The
    ! composition of PSC 1b is caculated in subroutine TER. 
    ! Equilibrium vapour pressure is used to calculate T(NAT) and T(ICE). 
    !    T(TER) = f(T,p(HNO3),p(H2O),H2SO4)
    !    T(NAT) = f(T,p(HNO3),p(H2O))
    !    T(ICE) = f(T,p(H2O))
    ! 
    ! The formation of NAT and ICE is controlled by the supersaturation
    ! of HNO3 and H2O (ssi = pi/peqi).
    !
    ! Aerosol number-concentrations are calculated for 3 parallel bins:
    ! ICE, STS, and NAT. Each has own nBins bins, identical over the species.
    ! Mass concentrations have additional substance dimension, of course.
    !
    !     This subroutine uses the following subroutines:
    !       - TER (ck/cs, rp, sad, composition)
    !       - NAT (ck/cs, rp, sad)
    !       - ICE (ck/cs, rp, sad)
    !       - GAMMA (heterogeneous reaction rate coefficients)
    !
    ! Authors
    ! The initial code is taken from FinROSE model, modified by M.Sofiev 

    implicit none
    !3Aug2017 by R.H. Added the HCl, HOCl, and HBr concentrations
    real, intent(inout) :: cncHNO3, cncH2SO4, cncH2O, cncHCl, cncHOCl, cncHBr
    real, intent(out) :: extra_sulphur_mass_STS, extra_particles_STS
    REAL, dimension(:,:,:), intent(inout) :: masses  ! masses of aerosols
    real, dimension(:,:), intent(inout) :: naero, core  ! number cnc of aerosols, mean volume of single particle
    real, dimension(:), intent(out) :: density  ! density of aerosols of each type
    logical, intent(out) :: ifICE, ifNAT, ifSTS
    real, intent(out) :: wtH2SO4
    type(Tchem_rules_AerDynMidAtm), intent(in) :: rules
    type(TLocalRates), intent(inout) :: rulesRates
    real, dimension(:), intent(in) :: metdat, fLowCncThresh
    real, intent(in) :: seconds
    logical, intent(inout) :: print_it

    ! Local variables
    integer :: iTmp, indT, iBin
    real :: ssnat,ssice,T, ppa, q, sulnum !, lphno3,lph2o  !, Tnat, Tice
    real :: rho1,patm,ptorr, hnm, hno3eq,ph2o !mlhcl 
    real :: cncH2O_tot, H2Opp, cncH2SO4_tot, cncHNO3_tot, cncHCl_tot, cncHOCl_tot, cncHBr_tot  !, HClpp,HBrpp
    ! psc particle radius and terminal velocity (1-TER, 2-NAT, 3-ICE)      
    real :: rp_HNO3_STS, rp_HNO3_NAT, rp_ICE
    !logical, save :: ifFirst = .true.    
    real :: n_tot, s_tot, fTmp
    
    ifICE = .false.
    ifNAT = .false.
    ifSTS = .true.  ! with forcing the sulnum we always have them
    !wtH2SO4 = 0.01    ! anything not too stupid
    wtH2SO4 = 60.0     ! anything not too stupid. This is now more consistent with density of 1500 kg/m3
    density(:) = 1500. ! anything not too stupid
    extra_sulphur_mass_STS = 0.
    extra_particles_STS = 0.
    
    n_tot = cncHNO3 + sum(masses(iHNO3_loc,1:nBins,1:3)) + & !/ rules%mass2vol(iHNO3ice) + &
                    & sum(masses(iNAT_loc,1:nBins,1:3)) !/ rules%mass2vol(iNATf)
    s_tot = cncH2SO4 + sum(masses(iH2SO4_loc,1:nBins,1:3))

    ! Index for temperature in the pre-defined arrays in rules
    !
    T = metdat(ind_tempr)
    ppa = metdat(ind_pres)
    q = metdat(ind_q)
    indT = int(T - 71.65)

    ! conditions (T,p and m.r.)
    patm = ppa / 101325.        !atm
    ptorr = ppa / 133.322       !torr
    rulesRates%cnc2vmr = gas_constant_uni * T / ppa  ! conversion in current conditions
    hnm = avogadro / rulesRates%cnc2vmr ! air density, #/m3

    !
    ! mixing ratios: water is from specific humidity, to be translated from kg/kg to vmr
    ! the rest is to be converted from concentrations mole/m3 to mixing ratio - and then 
    ! to partial pressure.
    ! Note that for H2O, SO4, HNO3, HCl, and HBr we need total amounts and gas phase separately
    !
    cncH2O_tot = cncH2O + sum(masses(iH2O_loc,1:nBins,1:3)) + 3.0 * sum(masses(iNAT_loc,1:nBins,1:3))
    !RH: TESTcnc2vmrChange
    !H2Opp = cncH2O_tot * rulesRates%cnc2vmr * ppa  ! partial pressure
    H2Opp = cncH2O_tot * gas_constant_uni * T  ! partial pressure
    cncHNO3_tot = cncHNO3 + sum(masses(iHNO3_loc,1:nBins,1:3)) + & 
                          & sum(masses(iNAT_loc,1:nBins,1:3))
!    HClpp = cncHCl * rulesRates%cnc2vmr * ppa
!    HClpp = cncHCl * gas_constant_uni * T
    cncH2SO4_tot = cncH2SO4 + sum(masses(iH2SO4_loc,1:nBins,1:3))
    !3Aug2017 by R.H. Added the HCl, HOCl, HBr concentrations.
    cncHCl_tot = cncHCl + sum(masses(iHCl_loc,1:nBins,binSTS))    !HCl is always only on STS
    cncHOCl_tot = cncHOCl + sum(masses(iHOCl_loc,1:nBins,binSTS)) !HOCl is always only on STS
    cncHBr_tot = cncHBr + sum(masses(iHBr_loc,1:nBins,binSTS))    !HBr is always only on STS
    
!    Tice = rules%ice1 / (alog10(H2Opp) - rules%ice2)
    
!c T(NAT) from hanson and mauersberger, GRL, 15 (1988) 855.
!      lphno3 = alog10(cncHNO3 * rulesRates%cnc2vmr * ptorr)
!      lph2o = alog10(cncH2O * rulesRates%cnc2vmr * ptorr)
!      Tnat = ((lphno3 - rules%nat1 * lph2o - rules%nat3) + &
!            & ((lphno3 - rules%nat1 * lph2o - rules%nat3)**2 - &
!             & 4.*(rules%nat2 * lph2o + rules%nat5) * rules%nat4)**.5) / (2.*(rules%nat2 * lph2o + rules%nat5))
#ifdef DEBUG_MORE
if(print_it)then
!call msg('rules%mass2vol(iH2O):',iH2Oice,rules%mass2vol(iH2Oice))
!call msg('rules%mass2vol(iH2SO4):',iSO4ice,rules%mass2vol(iSO4ice))
!call msg('rules%mass2vol(iHNO3):',iHNO3ice,rules%mass2vol(iHNO3ice))
!call msg('rules%mass2vol:',rules%mass2vol)
!RISTO Cl&Br: Added HCl, HOCl, and HBr!
call msg('cncHNO3, cncH2SO4, cncH2O, cncHCl, cncHOCl, cncHBr',(/cncHNO3, cncH2SO4, cncH2O, cncHCl, cncHOCl, cncHBr/))
call msg('thresholds HNO3, H2SO4, H2O, HCl, HOCl, HBr:',(/fLowCncThresh(iHNO3g), fLowCncThresh(iSO4f_aer), fLowCncThresh(iH2Og), &
                                                        & fLowCncThresh(iHCl), fLowCncThresh(iHOCl), fLowCncThresh(iHBr)/))
write(unit=strTmp,fmt='(1x,1000(E9.3,2x))')masses
call msg('Start PSC: masses:' + strMetaMasses)
call msg('Start PSC: masses:' + strTmp)
write(unit=strTmp,fmt='(1x,1000(E9.3,2x))')naero
call msg('Start PSC naero:' + strMetaNaero)
call msg('Start PSC naero:' + strTmp)
call msg('Start PSC NAT masses:',masses(iNAT_loc,iModeFine,binNAT), masses(iNAT_loc,iModeCoarse,binNAT))
call msg('NAT naero:',naero(iModeFine,binNAT), naero(iModeCoarse,binNAT))
endif
#endif

    if(cncH2O_tot < fLowCncThresh(iH2Og))return     ! no water => no love

!RISTO Cl&Br Q: Why the following lines are here??
!!    cncHCl_tot = cncHCl + sum(masses(iHCl,1:nBins,iSTS)) !/ mass2vol(HCl)
!    HBrpp = cncHBr * cnc2vmr * ppa
!!    cncHBr_tot = cncHBr + sum(masses(iHBr,1:nBins,iSTS)) !/ mass2vol(HBr)
!    HOClpp = cncHOCl * cnc2vmr * ppa
!!    cncHOCl_tot = cncHOCl + sum(masses(iHOCl,1:nBins,iSTS)) !/ mass2vol(HOCl)

    !
    ! threshold temperatures and partial pressure for ICE
    ! do ICE if it is cold or high super-saturation and then skip NAT & TER
    ! FinROSE requires either very cold conditions (T < Tice-3K) or existing ice from previous step with
    ! T<Tice, or high super-saturation. I.e., ice is created either in strong cold or by high SS.
    ! It can therefore stay as long as T<Tice.
    !
    ! So far I put T<Tice and high super-saturation as the only combined criterion. 
    ! That, however, can be wrong if FinROSE formulations are used.
    ! !!!!!!!!!!!!!!!!!!CHECK WITH LEIF.!!!!!!!!!!
    !
    ssice = H2Opp / rules%h2oeq_ice(indT)    ! super-saturation over ice

!    ifICE = (T <= Tice .and. sum(naero(1:nBins,binICE)) >= 1e-5) .or. &
!          & (ssice > rules%ssICE_crit) .or. &
!          & T <= Tice-3.
     
    ifICE = ssice > rules%ssICE_crit .or. &                       ! minimum super-saturation
          & ssice > 0.9 .and. sum(naero(1:nBins,binICE)) >= 1e-5  ! arbitrary numbers 0.9 and 1e-5: cold + exiting ice

#ifdef DEBUG_MORE
if(print_it)then
if(ifICE)then
  call msg('ICE active, ssice, sum(naero):', ssice, sum(naero(1:nBins,binICE)))
else
  call msg('ICE passive, ssice, sum(naero):', ssice, sum(naero(1:nBins,binICE)))
endif
    
write(unit=strTmp,fmt='(1x,1000(E9.3,2x))')masses
call msg('Before ICE: masses:' + strMetamasses)
call msg('Before ICE: masses:' + strTmp)
write(unit=strTmp,fmt='(1x,1000(E9.3,2x))')naero
call msg('naero:' + strMetaNaero)
call msg('naero:' + strTmp)
endif
#endif

    
    if(ifICE)then
      call ICE(naero(:,binICE), &         ! # concentration of ice particles, [m-3]
             !& rules%h2oeq_ice(indT)/(ppa*rulesRates%cnc2vmr), &   ! equilibrium H2O concentration, [mole/m3]
             !& rules%hno3eq_ice(indT)/(ppa*rulesRates%cnc2vmr), &  ! equilibrium HNO3 concentration, [mole/m3]
             & rules%h2oeq_ice(indT)/(gas_constant_uni*T), &   ! equilibrium H2O concentration, [mole/m3]
             & rules%hno3eq_ice(indT)/(gas_constant_uni*T), &  ! equilibrium HNO3 concentration, [mole/m3]
             & cncH2O, &                   ! actual H2O gas concentration, [mole/m3]
             & cncHNO3, &                  ! actual HNO3 gas concentration, [mole/m3]
             & cncH2SO4, &                 ! actual H2SO4 gas concentration, [mole/m3]
             & masses(iH2SO4_loc,iModeFine,binSTS), &  ! Sulphates have memory, have to treat separately
             & masses(:,:,binICE), &         ! masses of all species in particles, [m3/m3]
!             & sadICE, &                  ! surface area density of ice particles, [um2/m3]
!             & cnc2vol, &           ! conversion btw concentratin and volume in particle
             & rules, &                     ! aerosol rules
             & rulesRates, &
             & density(binICE), &    ! density of ICE aerosol
             & print_it)
if(.not. all((/cncH2O, cncHNO3, cncH2SO4/) >= 0.))then
  call msg('Problem with ICE: (/cncH2O, cncHNO3, cncH2SO4/)', (/cncH2O, cncHNO3, cncH2SO4/))
  call set_error('Problem with ICE','polar_stratospheric_clouds')
endif
if(abs(s_tot + extra_sulphur_mass_STS - (cncH2SO4 + sum(masses(iH2SO4_loc,1:nBins,1:3)))) > &
  & 1e-5 * (s_tot + extra_sulphur_mass_STS + cncH2SO4 + sum(masses(iH2SO4_loc,1:nBins,1:3))))then
  call msg('Wrong H2SO4. Place 1. After ICE: s_tot, extra,current-total:', &
                & (/s_tot,extra_sulphur_mass_STS, cncH2SO4 + sum(masses(iH2SO4_loc,1:nBins,1:3))/))
  if(abs(s_tot + extra_sulphur_mass_STS - (cncH2SO4 + sum(masses(iH2SO4_loc,1:nBins,1:3)))) > &
  & 1e-4 * (s_tot + extra_sulphur_mass_STS + cncH2SO4 + sum(masses(iH2SO4_loc,1:nBins,1:3)))) &
  call set_error('Wrong H2SO4. Place 1. After ICE','polar_stratospheric_clouds')
endif
if(abs(n_tot - (cncHNO3 + sum(masses(iHNO3_loc,1:nBins,1:3)) + sum(masses(iNAT_loc,1:nBins,1:3)) )) > &
  & 1e-5 * (n_tot + (cncHNO3 + sum(masses(iHNO3_loc,1:nBins,1:3)) + sum(masses(iNAT_loc,1:nBins,1:3)) )))then
  call msg('Wrong HNO3. Place 1. After ICE:',n_tot, (cncHNO3 + sum(masses(iHNO3_loc,1:nBins,1:3)) + &
                          & sum(masses(iNAT_loc,1:nBins,1:3)) ))
  if(abs(n_tot - (cncHNO3 + sum(masses(iHNO3_loc,1:nBins,1:3)) + sum(masses(iNAT_loc,1:nBins,1:3)) )) > &
  & 1e-4 * (n_tot + (cncHNO3 + sum(masses(iHNO3_loc,1:nBins,1:3)) + sum(masses(iNAT_loc,1:nBins,1:3)) ))) &
  call set_error('Wrong HNO3. Place 1. After ICE','polar_stratospheric_clouds')
endif
fTmp = cncH2O + sum(masses(iH2O_loc,1:nBins,1:3)) + 3.*sum(masses(iNAT_loc,1:nBins,1:3))
if(abs(cncH2O_tot - fTmp) > 1e-4 * (cncH2O_tot + fTmp))then
  call msg("Wrong H2O.  Thread No"&
      !$ & , omp_get_thread_num() &
         &)
  call msg('Wrong H2O. Place 1. After ICE (H2O_tot, H2O_sum, H2Ogas, H2O_part, NAT*3):', &
                                 & (/cncH2O_tot, fTmp, cncH2O, sum(masses(iH2O_loc,1:nBins,1:3)), &
                                                            &  sum(masses(iNAT_loc,1:nBins,1:3))/) )
  call set_error('Wrong H2O. Place 1. After ICE','polar_stratospheric_clouds')
endif
#ifdef DEBUG_MORE
if(print_it)then
call msg('After ICE, NAT masses:',masses(iNAT_loc,iModeFine,binNAT), masses(iNAT_loc,iModeCoarse,binNAT))
call msg('NAT naero:',naero(iModeFine,binNAT), naero(iModeCoarse,binNAT))
endif
#endif
    else
      !
      ! Conditions are not suitable for ice. Evaporate everything
      !
      do iBin = 1, nBins
        ! H2O -> gas
        cncH2O = cncH2O + masses(iH2O_loc,iBin,binICE) !/ rules%mass2vol(iH2Oice)
        masses(iH2O_loc,iBin,binICE) = 0.
        ! HNO3 -> gas
        cncHNO3 = cncHNO3 + masses(iHNO3_loc,iBin,binICE) !/ rules%mass2vol(iHNO3ice)
        masses(iHNO3_loc,iBin,binICE) = 0.
        ! H2SO4 in ice bin -> short-lived gas
        cncH2SO4 = cncH2SO4 + masses(iH2SO4_loc,iBin,binICE) !/ rules%mass2vol(iSO4ice)
        masses(iH2SO4_loc,iBin,binICE) = 0.
        ! Evaporate ice particles
        naero(iBin,binICE) = 0.
      enddo

if(abs(n_tot - (cncHNO3 + sum(masses(iHNO3_loc,1:nBins,1:3)) + sum(masses(iNAT_loc,1:nBins,1:3)))) > &
  & 1e-5 * (n_tot + cncHNO3 + sum(masses(iHNO3_loc,1:nBins,1:3)) + sum(masses(iNAT_loc,1:nBins,1:3))))then
  call msg('Wrong HNO3. Place 2. No-ICE starts:',n_tot, (cncHNO3 + sum(masses(iHNO3_loc,1:nBins,1:3)) + &
                          & sum(masses(iNAT_loc,1:nBins,1:3)) ))
  if(abs(n_tot - (cncHNO3 + sum(masses(iHNO3_loc,1:nBins,1:3)) + sum(masses(iNAT_loc,1:nBins,1:3)))) > &
  & 1e-4 * (n_tot + cncHNO3 + sum(masses(iHNO3_loc,1:nBins,1:3)) + sum(masses(iNAT_loc,1:nBins,1:3)))) &
  call set_error('Wrong HNO3. Place 2. No-ICE starts:','polar_stratospheric_clouds')
endif
if(abs(s_tot + extra_sulphur_mass_STS - (cncH2SO4 + sum(masses(iH2SO4_loc,1:nBins,1:3)))) > &
  & 1e-5 * (s_tot + extra_sulphur_mass_STS + cncH2SO4 + sum(masses(iH2SO4_loc,1:nBins,1:3))))then
  call msg('Wrong H2SO4. Place 2. No-ICE starts: s_tot, extra,current-total:', &
                & (/s_tot,extra_sulphur_mass_STS, cncH2SO4 + sum(masses(iH2SO4_loc,1:nBins,1:3))/))
  if(abs(s_tot + extra_sulphur_mass_STS - (cncH2SO4 + sum(masses(iH2SO4_loc,1:nBins,1:3)))) > &
  & 1e-4 * (s_tot + extra_sulphur_mass_STS + cncH2SO4 + sum(masses(iH2SO4_loc,1:nBins,1:3)))) &
  call set_error('Wrong H2SO4. Place 2. No ICE','polar_stratospheric_clouds')
endif
fTmp = cncH2O + sum(masses(iH2O_loc,1:nBins,1:3)) + 3.*sum(masses(iNAT_loc,1:nBins,1:3))
if(abs(cncH2O_tot - fTmp) > 1e-4 * (cncH2O_tot + fTmp))then
  call msg("Wrong H2O.  Thread No"&
      !$ & , omp_get_thread_num() &
         &)
  call msg('Wrong H2O. Place 2. No ICE (H2O_tot, H2O_sum, H2Ogas, H2O_part, NAT*3):', &
                                 & (/cncH2O_tot, fTmp, cncH2O, sum(masses(iH2O_loc,1:nBins,1:3)), &
                                                            &  sum(masses(iNAT_loc,1:nBins,1:3))/) )
  call set_error('Wrong H2O. Place 2. No ICE','polar_stratospheric_clouds')
endif

      !
      ! In case of no ICE, try NAT and STS particles
      ! The saturation vapour pressure over NAT for HNO3: hanson and mauersberger, GRL, 15 (1988) 855. 
      ! Note torr as pressure unit
      !
      ph2o = amax1(H2Opp / 133.322, 1.e-12)
      hno3eq = 10.**((rules%nat1 + rules%nat2*T)*alog10(ph2o) + rules%nat3 + &
                   & rules%nat4/T + rules%nat5*T) ! part pr, [torr]
      !RH: TESTcnc2vmrChange
      !hno3eq = hno3eq / (ptorr * rulesRates%cnc2vmr)           ! convert to concentration, [mole/m3]
      hno3eq = hno3eq*133.322/(gas_constant_uni*T)           ! convert to concentration, [mole/m3]
      !
      ! Do we have NAT particles already?
      ! From previous time step or due to nucleation in near-saturation conditions (0.9 is arbitrary)
      !
!      ifNAT = T <= Tnat .and. sum(naero(1:nBins,binNAT)) >= 1e-5   ! FinROSE - inverse of Hanson & Maursberger, 1988

      ifNAT = cncHNO3 >= 0.9 * hno3eq .or. sum(naero(1:nBins,binNAT)) >= 1e-5

#ifdef DEBUG_MORE
if(print_it)then
if(ifNAT)then
  call msg('NAT active, ssnat, sum(naero):', cncHNO3 / (0.9 * hno3eq), sum(naero(1:nBins,binNAT)))
else
  call msg('NAT passive, ssnat, sum(naero):', cncHNO3 / (0.9 * hno3eq), sum(naero(1:nBins,binNAT)))
endif
  call msg("cncHNO3, hno3eq",cncHNO3, hno3eq)
endif
#endif
      !
      ! Number concentrations of aerosols. Trouble is: so far, we have no chance to reproduce the bulk
      ! aerosol concentration in the stratosphere. That has to be taken from somewhere, e.g., from 
      ! observations using FinROSE approach.
      ! To make things somewhat more dynamic and responsive, we shall: 
      ! (i) use explicit sulphates wherever we have them higher than the observed background
      ! (ii) add all non-sulphate aerosols on-top of background sulphates
      ! This funny summation is evidently wrong since it assumes sulphate coating of those aerosols
      ! But this is the only way: uptake coefs for halogens are given only for STS/NAT/ICE surfaces
      ! Therefore, we just sum-up these to background sulphates. Care should be taken when considering
      ! long-term eruption history: then we shall actually try to compute the sulphates. Then only
      ! space debris will be the missing component, to be prescribed.
      ! We shall take special care to avoid impact of the forced sulphates on transported masses: use
      ! extra_sulphur_mass
      !
      ! Observed vertical profile for sulnum, a 4th degree polynomial fit to data from
      ! McLinden et al., Observations of Stratospheric Aerosols using CPFM...,
      ! JAS, 56 (1999) 233-240. Fit only to be used between 5 and 200 hPa.
      ! The fit gives # densities between approximately 0.1 and 15, with the
      ! maximum around 150 hPa.
      !
      if (ppa < 500.)then
        sulnum = 1.0e5     ! #/m3
      elseif (ppa > 20000.) then
        sulnum = 1.0e7     ! #/m3
      else
        sulnum = 1.e6 * amax1(0.1, (((4.6615e-16 * ppa - 2.3926e-11) * ppa  + &     ! #/m3
                                                & 3.2760e-7) * ppa - 1.1182e-4) * ppa - 3.7668e-2)
      endif  ! parameters

!!        !
!!        ! Presumably more accurate way: NAT and STS are tracers with actual number concentrations
!!        ! Note that STS are formed on existing sulphates, which concentration is known: just SO4 in 
!!        ! DMAT-S. However, it has to be initialised somewhere.
!!        !
!!        ! STS formed by absorption of HNO3 and H2O on the surface of existing sulphate aerosols.
!!        ! See (Carslaw et al.,1994; Meilinger et al.,1995; Peter,1997; Lowe and MacKenzie,2008, 
!!        ! Khosrawi ea, Atmos. Chem. Phys., 11, 8471-8487, 2011.
!!        !
!!call set_error('Generic STS formation is not ready yet','polar_stratospheric_clouds')
!!return
!!!smth like this:    sulnum = vCncTrn(iSO4) * 0.075 / (rhoSO4 * Pi * diamSO4**3)) ! mole/m3 -> #/m2, 75~96*3/4
!!        sulnum = amax1(1e4, amin1(20.e6, sulnum)) ! who knows...
!!        !
!!        ! NAT particles have two options for nucleation: rate either prescribed or supersat-dependent
!!        !
!!!        if(ifNAT_nucl_satur_dependent)then
!!!          !
!!!          ! A super-saturation-dependent heterogeneous nucleation model of Hoyle ea, ACP, 2013
!!!          !
!!!          !vCncTrn(iNATnbr) = vCncTrn(iNATnbr) + &
!!!          !             & (rules%J_NAT_hetero_ice(indT) + rules%J_NAT_hetero_foreign(indT)) * seconds
!!!          !and so on
!!!call set_error('NAT saturation-dependent nucleation does not work yet','')
!!!return
!!!        else
!!          !
!!          ! Constants used by various authors (funny ones excluded): see Hoyle ea, ACP, 2013
!!          ! Range from 3e-6 # cm-3 air hour-1 up to 2.5e-5 # cm-3 air hour-1
!!          ! we take randomly: 7.2e-6 # cm-3 hr-1 / 3600  = 2e-9 # cm-3 sec-1 = 1e-15 # m-3 sec-1
!!          ! Evidently, the nucleation brings about the smallest size bin
!!          !
!!          naero(1,binNAT) = naero(1,binNAT) + 2.0e-15 * seconds
!!!        endif   ! How to compute NAT # concentration


      !
      ! In case of no existing/new NAT particles, go along the ternary solution branch
      ! It will give equilibrium concentrations in aerosol and gas phases. However, distribution 
      ! should either be made via condensation mechanism as in BAD or simply distributed 
      ! proportionally to surface area density of each bin
      !
      if(.not. ifNAT)then
        !
        ! If too warm (or no previous-time NATs), evaporate all from NAT species
        !
        cncHNO3 = cncHNO3 + sum(masses(iNAT_loc,:,binNAT))
        cncH2O = cncH2O + 3.0 * sum(masses(iNAT_loc,:,binNAT))
        masses(:,:,binNat) = 0.0
        naero(:,binNAT) = 0.0

#ifdef DEBUG_MORE
if(print_it)then
write(unit=strTmp,fmt='(1x,1000(E9.3,2x))')masses
call msg('Before TER')
call msg('masses :' + strMetamasses)
call msg('masses :' + strTmp)
write(unit=strTmp,fmt='(1x,1000(E9.3,2x))')naero
call msg('naero:' + strMetaNaero)
call msg('naero:' + strTmp)
endif
#endif

if(abs(n_tot - (cncHNO3 + sum(masses(iHNO3_loc,1:nBins,1:3)) + sum(masses(iNAT_loc,1:nBins,1:3)) )) > &
  & 1e-5 * (n_tot + (cncHNO3 + sum(masses(iHNO3_loc,1:nBins,1:3)) + &
                          & sum(masses(iNAT_loc,1:nBins,1:3)) )))then
  call msg('Wrong HNO3. Place 3. Before first TER:',n_tot, (cncHNO3 + sum(masses(iHNO3_loc,1:nBins,1:3)) + &
                          & sum(masses(iNAT_loc,1:nBins,1:3)) ))
  if(abs(n_tot - (cncHNO3 + sum(masses(iHNO3_loc,1:nBins,1:3)) + sum(masses(iNAT_loc,1:nBins,1:3)) )) > &
  & 1e-4 * (n_tot + (cncHNO3 + sum(masses(iHNO3_loc,1:nBins,1:3)) + &
                          & sum(masses(iNAT_loc,1:nBins,1:3)) ))) &
  call set_error('Wrong HNO3. Place 3. Before first TER','polar_stratospheric_clouds')
endif
if(abs(s_tot + extra_sulphur_mass_STS - (cncH2SO4 + sum(masses(iH2SO4_loc,1:nBins,1:3)))) > &
  & 1e-5 * (s_tot + extra_sulphur_mass_STS + cncH2SO4 + sum(masses(iH2SO4_loc,1:nBins,1:3))))then
  call msg('Wrong H2SO4. Place 3. Before first TER: s_tot, extra,current-total:', &
                & (/s_tot,extra_sulphur_mass_STS, cncH2SO4 + sum(masses(iH2SO4_loc,1:nBins,1:3))/))
  if(abs(s_tot + extra_sulphur_mass_STS - (cncH2SO4 + sum(masses(iH2SO4_loc,1:nBins,1:3)))) > &
  & 1e-4 * (s_tot + extra_sulphur_mass_STS + cncH2SO4 + sum(masses(iH2SO4_loc,1:nBins,1:3)))) &
  call set_error('Wrong H2SO4. Place 3. Before first TER','polar_stratospheric_clouds')
endif
fTmp = cncH2O + sum(masses(iH2O_loc,1:nBins,1:3)) + 3.*sum(masses(iNAT_loc,1:nBins,1:3))
if(abs(cncH2O_tot - fTmp) > 1e-4 * (cncH2O_tot + fTmp))then
  call msg("Wrong H2O.  Thread No"&
      !$ & , omp_get_thread_num() &
         &)
  call msg('Wrong H2O. Place 3. Before first TER (H2O_tot, H2O_sum, H2Ogas, H2O_part, NAT*3):', &
                                 & (/cncH2O_tot, fTmp, cncH2O, sum(masses(iH2O_loc,1:nBins,1:3)), &
                                                            &  sum(masses(iNAT_loc,1:nBins,1:3))/) )
  call set_error('Wrong H2O. Place 3. Before first TER','polar_stratospheric_clouds')
endif
!11Oct2017 by R.H. Added checkups for HCl, HOCl, and HBr! !RISTO CONSIDER TO REMOVE THESE SINCE HCl,HOCl,and HBr are not touched the earlier routines!
if (abs( cncHCl_tot - (cncHCl + sum(masses(iHCl_loc,1:nBins,binSTS))) ) > &
  & 1e-5*( cncHCl_tot + (cncHCl + sum(masses(iHCl_loc,1:nBins,binSTS))) )) then
  call msg('Wrong HCl. Place 3. Before first TER:',cncHCl_tot, cncHCl + sum(masses(iHCl_loc,1:nBins,binSTS)) )
  call msg('RISTO TEST cncHCl_tot, cncHCl, sum(masses(iHCl_loc,1:nBins,binSTS)):', (/cncHCl_tot, cncHCl, sum(masses(iHCl_loc,1:nBins,binSTS))/) )
  if (abs( cncHCl_tot - (cncHCl + sum(masses(iHCl_loc,1:nBins,binSTS))) ) > &
  & 1e-4*( cncHCl_tot + (cncHCl + sum(masses(iHCl_loc,1:nBins,binSTS))) )) &
  call set_error('Wrong HCl. Place 3. Before first TER:','polar_stratospheric_clouds')
end if
if (abs( cncHOCl_tot - (cncHOCl + sum(masses(iHOCl_loc,1:nBins,binSTS))) ) > &
  & 1e-5*( cncHOCl_tot + (cncHOCl + sum(masses(iHOCl_loc,1:nBins,binSTS))) )) then
  call msg('Wrong HOCl. Place 3. Before first TER:',cncHOCl_tot, cncHOCl + sum(masses(iHOCl_loc,1:nBins,binSTS)) )
  call msg('RISTO TEST cncHOCl_tot, cncHOCl, sum(masses(iHOCl_loc,1:nBins,binSTS)):', (/cncHOCl_tot, cncHOCl, sum(masses(iHOCl_loc,1:nBins,binSTS))/) )
  if (abs( cncHOCl_tot - (cncHOCl + sum(masses(iHOCl_loc,1:nBins,binSTS))) ) > &
  & 1e-4*( cncHOCl_tot + (cncHOCl + sum(masses(iHOCl_loc,1:nBins,binSTS))) )) &
  call set_error('Wrong HOCl. Place 3. Before first TER:','polar_stratospheric_clouds')
end if
if (abs( cncHBr_tot - (cncHBr + sum(masses(iHBr_loc,1:nBins,binSTS))) ) > &
  & 1e-5*( cncHBr_tot + (cncHBr + sum(masses(iHBr_loc,1:nBins,binSTS))) )) then
  call msg('Wrong HBr. Place 3. Befor first TER:',cncHBr_tot, cncHBr + sum(masses(iHBr_loc,1:nBins,binSTS)) )
  call msg('RISTO TEST cncHBr_tot, cncHBr, sum(masses(iHBr_loc,1:nBins,binSTS)):', (/cncHBr_tot, cncHBr, sum(masses(iHBr_loc,1:nBins,binSTS))/) )
  if (abs( cncHBr_tot - (cncHBr + sum(masses(iHBr_loc,1:nBins,binSTS))) ) > &
  & 1e-4*( cncHBr_tot + (cncHBr + sum(masses(iHBr_loc,1:nBins,binSTS))) )) &
  call set_error('Wrong HBr. Place 3. Before first TER:','polar_stratospheric_clouds')
end if

if(print_it)then
call msg('masses before TER-1:',(/masses(:,:,:)/))
call msg('n_tot and cncHNO3 before TER-1',n_tot, cncHNO3) !Changed 3July2017 by R.H.
endif
        !
        ! Note that HNO3 from NAT is not needed: zeroed just above. ICE is also of no interest
        !
        ! RH2018 NOTE: the terms sum(masses(*_loc,1:nBins,binSTS)) are likely zero at this stage and not needed (CHECK!!)
        call TER(T, patm, hnm, &
               & cncH2SO4_tot, sulnum, &   ! H2SO4 mixing ratio, #-density of sulphate aerosol [m-3]
               & cncHNO3 + sum(masses(iHNO3_loc,1:nBins,binSTS)), &    ! semi-totals: gas + STS
               & cncH2O + sum(masses(iH2O_loc,1:nBins,binSTS)), &   ! water in NAT is not included
              !& cncHCl, & !cncHBr_tot, &
               !3Aug2017 by R.H. Added gas + STS semitotals for HCl, HOCl, and HBr plus the gas phase cnc calculated by TER
               & cncHCl + sum(masses(iHCl_loc,1:nBins,binSTS)), &   !gas + STS for HCl
               & cncHOCl + sum(masses(iHOCl_loc,1:nBins,binSTS)), & !gas + STS for HOCl
               & cncHBr + sum(masses(iHBr_loc,1:nBins,binSTS)), &   !gas + STS for HBr
               !& cncHNO3, cncH2O, cncHCl, & !cncHBr, &  ! gas phase, concentration [mole/m3]
               & cncHNO3, cncH2O, cncHCl, cncHOCl, cncHBr, &  ! gas phase, concentration [mole/m3]
               & naero(:,binSTS), &       ! # density of aerosols
               & masses(:,:,binSTS), &      ! species volumes in aerosols
               & extra_sulphur_mass_STS, extra_particles_STS, &  ! sulphate deficit, added to masses, bin1
               & density(binSTS), &        ! density of STS aerosol
               & wtH2SO4, &                ! weight fraction of H2SO4, %
!               & sad1, &                  ! surface area density of aerosol
!RISTO Cl&Br Q: Why the molality of HCl is outputted?? Not used anywhere!
               !& mlhcl, &                   ! molality of HCl in aerosol
               & rules, &
               & rulesRates, &
               & print_it)
        cncH2SO4 = 0.0   ! STS consume all H2SO4 and send it to STS particles

if(print_it)then
call msg('masses after  TER-1:',(/masses(:,:,:)/))
call msg('n_tot and cncHNO3 after TER-1',n_tot, cncHNO3) !Changed 13July2017 by R.H.
endif

if(abs(n_tot - (cncHNO3 + sum(masses(iHNO3_loc,1:nBins,1:3)) + sum(masses(iNAT_loc,1:nBins,1:3)) )) > &
  & 1e-5 * (n_tot + cncHNO3 + sum(masses(iHNO3_loc,1:nBins,1:3)) + sum(masses(iNAT_loc,1:nBins,1:3)) ))then
  call msg('Wrong HNO3. Place 4. After first TER:',n_tot, (cncHNO3 + sum(masses(iHNO3_loc,1:nBins,1:3)) + &
                          & sum(masses(iNAT_loc,1:nBins,1:3)) ))
  if(abs(n_tot - (cncHNO3 + sum(masses(iHNO3_loc,1:nBins,1:3)) + sum(masses(iNAT_loc,1:nBins,1:3)) )) > &
  & 1e-4 * (n_tot + cncHNO3 + sum(masses(iHNO3_loc,1:nBins,1:3)) + sum(masses(iNAT_loc,1:nBins,1:3)) ))&
  call set_error('Wrong HNO3. Place 4. After first TER','polar_stratospheric_clouds')
endif
if(abs(s_tot + extra_sulphur_mass_STS - (cncH2SO4 + sum(masses(iH2SO4_loc,1:nBins,1:3)))) > &
  & 1e-5 * (s_tot + extra_sulphur_mass_STS + cncH2SO4 + sum(masses(iH2SO4_loc,1:nBins,1:3))))then
  call msg('Wrong H2SO4. Place 4. After first TER: s_tot, extra,current-total:', &
                & (/s_tot,extra_sulphur_mass_STS, cncH2SO4 + sum(masses(iH2SO4_loc,1:nBins,1:3))/))
  if(abs(s_tot + extra_sulphur_mass_STS - (cncH2SO4 + sum(masses(iH2SO4_loc,1:nBins,1:3)))) > &
  & 1e-4 * (s_tot + extra_sulphur_mass_STS + cncH2SO4 + sum(masses(iH2SO4_loc,1:nBins,1:3)))) &
  call set_error('Wrong H2SO4. Place 4. After first TER','polar_stratospheric_clouds')
endif
fTmp = cncH2O + sum(masses(iH2O_loc,1:nBins,1:3)) + 3.*sum(masses(iNAT_loc,1:nBins,1:3))
if(abs(cncH2O_tot - fTmp) > 1e-4 * (cncH2O_tot + fTmp))then
  call msg("Wrong H2O.  Thread No"&
      !$ & , omp_get_thread_num() &
         &)
  call msg('Wrong H2O. Place 4. After first TER (H2O_tot, H2O_sum, H2Ogas, H2O_part, NAT*3):', &
                                 & (/cncH2O_tot, fTmp, cncH2O, sum(masses(iH2O_loc,1:nBins,1:3)), &
                                                            &  sum(masses(iNAT_loc,1:nBins,1:3))/) )
  call set_error('Wrong H2O. Place 4. After first TER','polar_stratospheric_clouds')
endif
!3Aug2017 by R.H. Added checkups for HCl, HOCl, and HBr!
if (abs( cncHCl_tot - (cncHCl + sum(masses(iHCl_loc,1:nBins,binSTS))) ) > &
  & 1e-5*( cncHCl_tot + (cncHCl + sum(masses(iHCl_loc,1:nBins,binSTS))) )) then
  call msg('Wrong HCl. Place 4. After first TER:',cncHCl_tot, cncHCl + sum(masses(iHCl_loc,1:nBins,binSTS)) )
  call msg('RISTO TEST cncHCl_tot, cncHCl, sum(masses(iHCl_loc,1:nBins,binSTS)):', (/cncHCl_tot, cncHCl, sum(masses(iHCl_loc,1:nBins,binSTS))/) )
  if (abs( cncHCl_tot - (cncHCl + sum(masses(iHCl_loc,1:nBins,binSTS))) ) > &
  & 1e-4*( cncHCl_tot + (cncHCl + sum(masses(iHCl_loc,1:nBins,binSTS))) )) &
  call set_error('Wrong HCl. Place 4. After first TER:','polar_stratospheric_clouds')
end if
if (abs( cncHOCl_tot - (cncHOCl + sum(masses(iHOCl_loc,1:nBins,binSTS))) ) > &
  & 1e-5*( cncHOCl_tot + (cncHOCl + sum(masses(iHOCl_loc,1:nBins,binSTS))) )) then
  call msg('Wrong HOCl. Place 4. After first TER:',cncHOCl_tot, cncHOCl + sum(masses(iHOCl_loc,1:nBins,binSTS)) )
  call msg('RISTO TEST cncHOCl_tot, cncHOCl, sum(masses(iHOCl_loc,1:nBins,binSTS)):', (/cncHOCl_tot, cncHOCl, sum(masses(iHOCl_loc,1:nBins,binSTS))/) )
  if (abs( cncHOCl_tot - (cncHOCl + sum(masses(iHOCl_loc,1:nBins,binSTS))) ) > &
  & 1e-4*( cncHOCl_tot + (cncHOCl + sum(masses(iHOCl_loc,1:nBins,binSTS))) )) &
  call set_error('Wrong HOCl. Place 4. After first TER:','polar_stratospheric_clouds')
end if
if (abs( cncHBr_tot - (cncHBr + sum(masses(iHBr_loc,1:nBins,binSTS))) ) > &
  & 1e-5*( cncHBr_tot + (cncHBr + sum(masses(iHBr_loc,1:nBins,binSTS))) )) then
  call msg('Wrong HBr. Place 4. After first TER:',cncHBr_tot, cncHBr + sum(masses(iHBr_loc,1:nBins,binSTS)) )
  call msg('RISTO TEST cncHBr_tot, cncHBr, sum(masses(iHBr_loc,1:nBins,binSTS)):', (/cncHBr_tot, cncHBr, sum(masses(iHBr_loc,1:nBins,binSTS))/) )
  write(*,*) 'cncHBr_tot = ', cncHBr_tot
  write(*,*) 'cncHBr     = ', cncHBr
  write(*,*) 'cncHBr_aer = ', sum(masses(iHBr_loc,1:nBins,binSTS))
  write(*,*) 'sum        = ', cncHBr+sum(masses(iHBr_loc,1:nBins,binSTS))
  if (abs( cncHBr_tot - (cncHBr + sum(masses(iHBr_loc,1:nBins,binSTS))) ) > &
  & 1e-4*( cncHBr_tot + (cncHBr + sum(masses(iHBr_loc,1:nBins,binSTS))) )) &
  call set_error('Wrong HBr. Place 4. After first TER:','polar_stratospheric_clouds')
end if

if(.not. all((/cncH2O_tot, cncHNO3, cncH2SO4, cncHCl, cncHOCl, cncHBr/) >= 0.))then
  call msg('Problem with TER-1: (/cncH2O_tot, cncHNO3, cncH2SO4, cncHCl/)', (/cncH2O_tot, cncHNO3, cncH2SO4, cncHCl, cncHOCl, cncHBr/))
  call set_error('Problem with TER-1','polar_stratospheric_clouds')
endif

#ifdef DEBUG_MORE
if(print_it)then
call msg('PSC-1, wtH2SO4:',wtH2SO4)
write(unit=strTmp,fmt='(1x,1000(E9.3,2x))')masses
call msg('After TER')
call msg('masses :' + strMetamasses)
call msg('masses :' + strTmp)
write(unit=strTmp,fmt='(1x,1000(E9.3,2x))')naero
call msg('naero:' + strMetaNaero)
call msg('naero:' + strTmp)
call msg('Extra volume and particles:',extra_sulphur_mass_STS, extra_particles_STS)
call msg('After TER-1, NAT masses:',masses(iNAT_loc,iModeFine,binNAT), masses(iNAT_loc,iModeCoarse,binNAT))
call msg('NAT naero:',naero(iModeFine,binNAT), naero(iModeCoarse,binNAT))
endif
#endif
      endif   ! not ifNAT
      !
      ! do NAT if hno3 ss-ratio reaches threshold, or ifnat 
      !
      if( cncHNO3 / hno3eq > rules%ssNAT_crit .or. ifnat ) then

        call NAT(hno3eq, &
               & cncHNO3, cncH2O, &
               & masses(:,:,binNAT), naero(:,binNAT), &
               & rules, rulesRates, density(binNAT), seconds, print_it) !,sadn)
        ifnat = .true.
if(.not. all((/cncH2O, cncHNO3, cncH2SO4/) >= 0.))then
  call msg('Problem with NAT: (/cncH2O, cncHNO3, cncH2SO4/)', (/cncH2O, cncHNO3, cncH2SO4/))
  call set_error('Problem with NAT','polar_stratospheric_clouds')
endif

#ifdef DEBUG_MORE
if(print_it)then
write(unit=strTmp,fmt='(1x,1000(E9.3,2x))')masses
call msg('After NAT')
call msg(' masses:' + strMetamasses)
call msg(' masses:' + strTmp)
write(unit=strTmp,fmt='(1x,1000(E9.3,2x))')naero
call msg('naero:' + strMetaNaero)
call msg('naero:' + strTmp)
call msg('After NAT, NAT masses:',masses(iNAT_loc,iModeFine,binNAT), masses(iNAT_loc,iModeCoarse,binNAT))
call msg('NAT naero:',naero(iModeFine,binNAT), naero(iModeCoarse,binNAT))
endif
#endif
      endif  ! if NAT

if(abs(n_tot - (cncHNO3 + sum(masses(iHNO3_loc,1:nBins,1:3)) + sum(masses(iNAT_loc,1:nBins,1:3)) )) > &
  & 1e-4 * (n_tot + cncHNO3 + sum(masses(iHNO3_loc,1:nBins,1:3)) + sum(masses(iNAT_loc,1:nBins,1:3)) ))then
  call msg('Wrong HNO3. Place 5. After nat:',n_tot, (cncHNO3 + sum(masses(iHNO3_loc,1:nBins,1:3)) + &
                          & sum(masses(iNAT_loc,1:nBins,1:3)) ))
  if(abs(n_tot - (cncHNO3 + sum(masses(iHNO3_loc,1:nBins,1:3)) + sum(masses(iNAT_loc,1:nBins,1:3)) )) > &
  & 1e-4 * (n_tot + cncHNO3 + sum(masses(iHNO3_loc,1:nBins,1:3)) + sum(masses(iNAT_loc,1:nBins,1:3)) )) &
  call set_error('Wrong HNO3. Place 5. After nat.','polar_stratospheric_clouds')
endif
if(abs(s_tot + extra_sulphur_mass_STS - (cncH2SO4 + sum(masses(iH2SO4_loc,1:nBins,1:3)))) > &
  & 1e-4 * (s_tot + extra_sulphur_mass_STS + cncH2SO4 + sum(masses(iH2SO4_loc,1:nBins,1:3))))then
  call msg('Wrong H2SO4. Place 5. After NAT: s_tot, extra,current-total:', &
                & (/s_tot,extra_sulphur_mass_STS, cncH2SO4 + sum(masses(iH2SO4_loc,1:nBins,1:3))/))
  if(abs(s_tot + extra_sulphur_mass_STS - (cncH2SO4 + sum(masses(iH2SO4_loc,1:nBins,1:3)))) > &
  & 1e-4 * (s_tot + extra_sulphur_mass_STS + cncH2SO4 + sum(masses(iH2SO4_loc,1:nBins,1:3)))) &
  call set_error('Wrong H2SO4. Place 5. After NAT','polar_stratospheric_clouds')
endif
fTmp = cncH2O + sum(masses(iH2O_loc,1:nBins,1:3)) + 3.*sum(masses(iNAT_loc,1:nBins,1:3))
if(abs(cncH2O_tot - fTmp) > 1e-4 * (cncH2O_tot + fTmp))then
  call msg("Wrong H2O.  Thread No"&
      !$ & , omp_get_thread_num() &
         &)
  call msg('Wrong H2O. Place 5. After NAT (H2O_tot, H2O_sum, H2Ogas, H2O_part, NAT*3):', &
                                 & (/cncH2O_tot, fTmp, cncH2O, sum(masses(iH2O_loc,1:nBins,1:3)), &
                                                            &  sum(masses(iNAT_loc,1:nBins,1:3))/) )
  call set_error('Wrong H2O. Place 5. After NAT','polar_stratospheric_clouds')
endif

if(print_it)then
call msg(strMetamasses)
call msg('n_tot and cncHNO3 before TER-2',n_tot, cncHNO3)
call msg(' masses before TER-2:',(/masses(:,:,:)/))
endif
      !
      ! (re)calculate ternary composition if nat 
      ! an error is introduced to the hno3 partition due to NAT sub. Note that NAT should not be touched
      !
      if ( ifnat ) then
        call TER(T, patm, hnm, &
               & cncH2SO4_tot, sulnum, &   ! H2SO4 mixing ratio, #-density of sulphate aerosol [m-3]
               & cncHNO3 + sum(masses(iHNO3_loc,1:nBins,binSTS)), &
               & cncH2O + sum(masses(iH2O_loc,1:nBins,binSTS)), &   ! water in NAT is not included
               !3Aug2017 by R.H. Added gas + STS semitotals for HCl, HOCl, and HBr plus the gas phase cnc calculated by TER
               & cncHCl + sum(masses(iHCl_loc,1:nBins,binSTS)), &   !gas + STS for HCl
               & cncHOCl + sum(masses(iHOCl_loc,1:nBins,binSTS)), & !gas + STS for HOCl
               & cncHBr + sum(masses(iHBr_loc,1:nBins,binSTS)), &   !gas + STS for HBr
               !old v5.5& cncHNO3, cncH2O, cncHCl, & !cncHBr, &  ! gas phase, concentration [mole/m3]
               & cncHNO3, cncH2O, cncHCl, cncHOCl, cncHBr, &  ! gas phase, concentration [mole/m3]
               & naero(:,binSTS), masses(:,:,binSTS), &    ! # density of aerosol and species volumes
               & extra_sulphur_mass_STS, extra_particles_STS, &  ! sulphate deficit, added to masses, bin1
               & density(binSTS), &        ! density of STS aerosol
               & wtH2SO4, &                ! weight fraction of H2SO4, %
!               & sad1, &                  ! surface area density of aerosol
               !& mlhcl, &
               & rules, &
               & rulesRates, &
               & print_it)
        cncH2SO4 = 0.0   ! STS consume all H2SO4 and send it to STS particles
if(print_it)then
call msg('masses after  TER-2:',(/masses(:,:,:)/))
call msg('n_tot and cncHNO3 after TER-2',n_tot, cncHNO3)
endif
if(.not. all((/cncH2O_tot, cncHNO3, cncH2SO4/) >= 0.))then
  call msg('Problem with TER-2: (/cncH2O, cncHNO3, cncH2SO4/)', (/cncH2O, cncHNO3, cncH2SO4/))
  call set_error('Problem with TER-2','polar_stratospheric_clouds')
endif

#ifdef DEBUG_MORE
if(print_it)then
call msg('PSC-2, wtH2SO4:',wtH2SO4)
write(unit=strTmp,fmt='(1x,1000(E9.3,2x))')masses
call msg('After TER 2')
call msg(' masses:' + strMetamasses)
call msg(' masses:' + strTmp)
write(unit=strTmp,fmt='(1x,1000(E9.3,2x))')naero
call msg('naero:' + strMetaNaero)
call msg('naero:' + strTmp)
call msg('After TER-2, NAT masses:',masses(iNAT_loc,iModeFine,binNAT), masses(iNAT_loc,iModeCoarse,binNAT))
call msg('NAT naero:',naero(iModeFine,binNAT), naero(iModeCoarse,binNAT))
endif
#endif
      endif  ! if NAT
if(abs(n_tot - (cncHNO3 + sum(masses(iHNO3_loc,1:nBins,1:3)) + sum(masses(iNAT_loc,1:nBins,1:3)) )) > &
  & 1e-5 * (n_tot + cncHNO3 + sum(masses(iHNO3_loc,1:nBins,1:3)) + sum(masses(iNAT_loc,1:nBins,1:3)) ))then
  call msg('Wrong HNO3. Place 6. After second TER:',n_tot, (cncHNO3 + sum(masses(iHNO3_loc,1:nBins,1:3)) + &
                          & sum(masses(iNAT_loc,1:nBins,1:3)) ))
  if(abs(n_tot - (cncHNO3 + sum(masses(iHNO3_loc,1:nBins,1:3)) + sum(masses(iNAT_loc,1:nBins,1:3)) )) > &
  & 1e-4 * (n_tot + cncHNO3 + sum(masses(iHNO3_loc,1:nBins,1:3)) + sum(masses(iNAT_loc,1:nBins,1:3)) ))&
  call set_error('Wrong HNO3. Place 6. After second TER','polar_stratospheric_clouds')
endif
if(abs(s_tot + extra_sulphur_mass_STS - (cncH2SO4 + sum(masses(iH2SO4_loc,1:nBins,1:3)))) > &
  & 1e-5 * (s_tot + extra_sulphur_mass_STS + cncH2SO4 + sum(masses(iH2SO4_loc,1:nBins,1:3))))then
  call msg('Wrong H2SO4. Place 6. After second TER: s_tot, extra,current-total:', &
                & (/s_tot,extra_sulphur_mass_STS, cncH2SO4 + sum(masses(iH2SO4_loc,1:nBins,1:3))/))
  if(abs(s_tot + extra_sulphur_mass_STS - (cncH2SO4 + sum(masses(iH2SO4_loc,1:nBins,1:3)))) > &
  & 1e-4 * (s_tot + extra_sulphur_mass_STS + cncH2SO4 + sum(masses(iH2SO4_loc,1:nBins,1:3))))&
  call set_error('Wrong H2SO4. Place 6. After second TER','polar_stratospheric_clouds')
endif
fTmp = cncH2O + sum(masses(iH2O_loc,1:nBins,1:3)) + 3.*sum(masses(iNAT_loc,1:nBins,1:3))
if(abs(cncH2O_tot - fTmp) > 1e-4 * (cncH2O_tot + fTmp))then
  call msg("Wrong H2O.  Thread No" &
      !$ & , omp_get_thread_num() &
         &)
  call msg('Wrong H2O. Place 6. After second TER (H2O_tot, H2O_sum, H2Ogas, H2O_part, NAT*3):', &
                                 & (/cncH2O_tot, fTmp, cncH2O, sum(masses(iH2O_loc,1:nBins,1:3)), &
                                                            &  sum(masses(iNAT_loc,1:nBins,1:3))/) )
  call set_error('Wrong H2O. Place 6. After second TER','polar_stratospheric_clouds')
endif
!RISTO Cl&Br: Added checkups for HCl, HOCl, and HBr!
if (abs( cncHCl_tot - (cncHCl + sum(masses(iHCl_loc,1:nBins,binSTS))) ) > &
  & 1e-5*( cncHCl_tot + (cncHCl + sum(masses(iHCl_loc,1:nBins,binSTS))) )) then
  call msg('Wrong HCl. Place 6. After second TER:',cncHCl_tot, cncHCl + sum(masses(iHCl_loc,1:nBins,binSTS)) )
  if (abs( cncHCl_tot - (cncHCl + sum(masses(iHCl_loc,1:nBins,binSTS))) ) > &
  & 1e-4*( cncHCl_tot + (cncHCl + sum(masses(iHCl_loc,1:nBins,binSTS))) )) &
  call set_error('Wrong HCl. Place 6. After second TER:','polar_stratospheric_clouds')
end if
if (abs( cncHOCl_tot - (cncHOCl + sum(masses(iHOCl_loc,1:nBins,binSTS))) ) > &
  & 1e-5*( cncHOCl_tot + (cncHOCl + sum(masses(iHOCl_loc,1:nBins,binSTS))) )) then
  call msg('Wrong HOCl. Place 6. After second TER:',cncHOCl_tot, cncHOCl + sum(masses(iHOCl_loc,1:nBins,binSTS)) )
  !call msg('RISTO TEST cncHOCl_tot, cncHOCl, sum(masses(iHOCl_loc,1:nBins,binSTS)):', (/cncHOCl_tot, cncHOCl, sum(masses(iHOCl_loc,1:nBins,binSTS))/) )
  if (abs( cncHOCl_tot - (cncHOCl + sum(masses(iHOCl_loc,1:nBins,binSTS))) ) > &
  & 1e-4*( cncHOCl_tot + (cncHOCl + sum(masses(iHOCl_loc,1:nBins,binSTS))) )) &
  call set_error('Wrong HOCl. Place 6. After second TER:','polar_stratospheric_clouds')
end if
if (abs( cncHBr_tot - (cncHBr + sum(masses(iHBr_loc,1:nBins,binSTS))) ) > &
  & 1e-5*( cncHBr_tot + (cncHBr + sum(masses(iHBr_loc,1:nBins,binSTS))) )) then
  call msg('Wrong HBr. Place 6. After second TER:',cncHBr_tot, cncHBr + sum(masses(iHBr_loc,1:nBins,binSTS)) )
  !call msg('RISTO TEST cncHBr_tot, cncHBr, sum(masses(iHBr_loc,1:nBins,binSTS)):', (/cncHBr_tot, cncHBr, sum(masses(iHBr_loc,1:nBins,binSTS))/) )
  if (abs( cncHBr_tot - (cncHBr + sum(masses(iHBr_loc,1:nBins,binSTS))) ) > &
  & 1e-4*( cncHBr_tot + (cncHBr + sum(masses(iHBr_loc,1:nBins,binSTS))) )) &
  call set_error('Wrong HBr. Place 6. After second TER:','polar_stratospheric_clouds')
end if

    endif  ! if ICE

#ifdef DEBUG_MORE
if(print_it)then
write(unit=strTmp,fmt='(1x,1000(E9.3,2x))')masses
call msg('End PSC')
call msg(' masses:' + strMetamasses)
call msg(' masses:' + strTmp)
write(unit=strTmp,fmt='(1x,1000(E9.3,2x))')naero
call msg('End PSC naero:' + strMetaNaero)
call msg('End PSC naero:' + strTmp)
call msg('End PSC NAT masses:',masses(iNAT_loc,iModeFine,binNAT), masses(iNAT_loc,iModeCoarse,binNAT))
call msg('NAT naero:',naero(iModeFine,binNAT), naero(iModeCoarse,binNAT))
endif
#endif
if(abs(n_tot - (cncHNO3 + sum(masses(iHNO3_loc,1:nBins,1:3)) + sum(masses(iNAT_loc,1:nBins,1:3)) )) > &
  & 1e-5 * (n_tot + (cncHNO3 + sum(masses(iHNO3_loc,1:nBins,1:3)) + sum(masses(iNAT_loc,1:nBins,1:3)) )))then
  call msg('Wrong HNO3. Place 7. End of sub:',n_tot, (cncHNO3 + sum(masses(iHNO3_loc,1:nBins,1:3)) + &
                          & sum(masses(iNAT_loc,1:nBins,1:3)) ))
  if(abs(n_tot - (cncHNO3 + sum(masses(iHNO3_loc,1:nBins,1:3)) + sum(masses(iNAT_loc,1:nBins,1:3)) )) > &
  & 1e-4 * (n_tot + (cncHNO3 + sum(masses(iHNO3_loc,1:nBins,1:3)) + sum(masses(iNAT_loc,1:nBins,1:3)) )))&
  call set_error('Wrong HNO3. Place 7. End of sub','polar_stratospheric_clouds')
endif
if(abs(s_tot + extra_sulphur_mass_STS - (cncH2SO4 + sum(masses(iH2SO4_loc,1:nBins,1:3)))) > &
  & 1e-5 * (s_tot + extra_sulphur_mass_STS + cncH2SO4 + sum(masses(iH2SO4_loc,1:nBins,1:3))))then
  call msg('Wrong H2SO4. Place 7. End of sub: s_tot, extra,current-total:', &
                & (/s_tot,extra_sulphur_mass_STS, cncH2SO4 + sum(masses(iH2SO4_loc,1:nBins,1:3))/))
  if(abs(s_tot + extra_sulphur_mass_STS - (cncH2SO4 + sum(masses(iH2SO4_loc,1:nBins,1:3)))) > &
  & 1e-4 * (s_tot + extra_sulphur_mass_STS + cncH2SO4 + sum(masses(iH2SO4_loc,1:nBins,1:3)))) &
  call set_error('Wrong H2SO4. Place 7. End of sub','polar_stratospheric_clouds')
endif
fTmp = cncH2O + sum(masses(iH2O_loc,1:nBins,1:3)) + 3.*sum(masses(iNAT_loc,1:nBins,1:3))
if(abs(cncH2O_tot - fTmp) > 1e-4 * (cncH2O_tot + fTmp))then
  call msg("Wrong H2O.  Thread No"&
      !$ & , omp_get_thread_num() &
         &)
  call msg('Wrong H2O. Place 7. End of sub (H2O_tot, H2O_sum, H2Ogas, H2O_part, NAT*3):', &
                                 & (/cncH2O_tot, fTmp, cncH2O, sum(masses(iH2O_loc,1:nBins,1:3)), &
                                                            &  sum(masses(iNAT_loc,1:nBins,1:3))/) )
  call set_error('Wrong H2O. Place 7. End of sub','polar_stratospheric_clouds')
endif
!RISTO Cl&Br: Added checkups for HCl, HOCl, and HBr!
if (abs( cncHCl_tot - (cncHCl + sum(masses(iHCl_loc,1:nBins,binSTS))) ) > &
  & 1e-5*( cncHCl_tot + (cncHCl + sum(masses(iHCl_loc,1:nBins,binSTS))) )) then
  call msg('Wrong HCl. Place 7. End of sub:',cncHCl_tot, cncHCl + sum(masses(iHCl_loc,1:nBins,binSTS)) )
  if (abs( cncHCl_tot - (cncHCl + sum(masses(iHCl_loc,1:nBins,binSTS))) ) > &
  & 1e-4*( cncHCl_tot + (cncHCl + sum(masses(iHCl_loc,1:nBins,binSTS))) )) &
  call set_error('Wrong HCl. Place 7. End of sub,','polar_stratospheric_clouds')
end if
if (abs( cncHOCl_tot - (cncHOCl + sum(masses(iHOCl_loc,1:nBins,binSTS))) ) > &
  & 1e-5*( cncHOCl_tot + (cncHOCl + sum(masses(iHOCl_loc,1:nBins,binSTS))) )) then
  call msg('Wrong HOCl. Place 7. End of sub:',cncHOCl_tot, cncHOCl + sum(masses(iHOCl_loc,1:nBins,binSTS)) )
  call msg('RISTO TEST cncHOCl_tot, cncHOCl, sum(masses(iHOCl_loc,1:nBins,binSTS)):', (/cncHOCl_tot, cncHOCl, sum(masses(iHOCl_loc,1:nBins,binSTS))/) )
  if (abs( cncHOCl_tot - (cncHOCl + sum(masses(iHOCl_loc,1:nBins,binSTS))) ) > &
  & 1e-4*( cncHOCl_tot + (cncHOCl + sum(masses(iHOCl_loc,1:nBins,binSTS))) )) &
  call set_error('Wrong HOCl. Place 7. End of sub,','polar_stratospheric_clouds')
end if
if (abs( cncHBr_tot - (cncHBr + sum(masses(iHBr_loc,1:nBins,binSTS))) ) > &
  & 1e-5*( cncHBr_tot + (cncHBr + sum(masses(iHBr_loc,1:nBins,binSTS))) )) then
  call msg('Wrong HBr. Place 7. End of sub:',cncHBr_tot, cncHBr + sum(masses(iHBr_loc,1:nBins,binSTS)) )
  call msg('RISTO TEST cncHBr_tot, cncHBr, sum(masses(iHBr_loc,1:nBins,binSTS)):', (/cncHBr_tot, cncHBr, sum(masses(iHBr_loc,1:nBins,binSTS))/) )
  if (abs( cncHBr_tot - (cncHBr + sum(masses(iHBr_loc,1:nBins,binSTS))) ) > &
  & 1e-4*( cncHBr_tot + (cncHBr + sum(masses(iHBr_loc,1:nBins,binSTS))) )) &
  call set_error('Wrong HBr. Place 7. End of sub,','polar_stratospheric_clouds')
end if

  end subroutine polar_stratospheric_clouds


  !*******************************************************************

!RISTO Cl&Br: Added HOCl and HBr (totals and gass phase)! Why mcl is outputted?
  subroutine TER(Temp, &   !       (r) = temperature(unchanged), k
               & patm, &   !       (r) = pressure (unchaged), atm
               & hnm, &   !     (r) = air # density (unchanged), m-3
               & h2so4_tot, &   !   (r) = total h2so4 concentration (unchanged), mole/m3
               & sulnum, &   !  (r) = # density of sulfate aerosols (unchanged), m-3
               & hno3_tot, &   !   (r) = total hno3 concentration (unchanged), mole/m3
               & h2o_tot, &   !    (r) = water vapour concentration (unchanged), mole/m3
               & hcl_tot, &   !    (r) = total hcl concentration (unchanged), mole/m3
               !3Aug2017 by R.H. Added the total HOCl and HBr concentrations 
               & hocl_tot, &   !    (r) = total hocl concentration (unchanged), mole/m3
               & hbr_tot, &   !    (r) = total hcbr concentration (unchanged), mole/m3
               !!!!!!!!!!!!!!!!!!!
               & hno3g, &   !   (r) = hno3 concentration in gas phase (changed), mole/m3
               & h2og, &    !   (r) = h2o concentration in gas ohase (changed), mole/m3
               & hclg, &   !    (r) = hcl concentration in gas phase (changed), mole/m3
               !3Aug2017 by R.H. Added the HOCl and HBr gas concentrations 
               & hoclg, &   !    (r) = hcl concentration in gas phase (changed), mole/m3
               & hbrg, &   !    (r) = hbr concentration in gas phase (changed), mole/m3
               !!!!!!!!!!!!!!!!!!!
               & naero, &   !   (r) = aerosol number concentration for each STS bin, m-3
               & masses, &    !   (r) = species in each bin of STS aerosols mole/m3
               & extra_sulphur_mass, extra_particles, &  ! (r) = deficit of sulphur wrt observations
               & rho, &     !   (r) = density of the ternary solution (changed), kg/m3
               & wtH2SO4, &   !    (r) = weight percentage of h2so4 (changed), %
               !& mcl,  &     !     (r) = molality of hcl in aerosol (changed), mol kg-1
               & rules, &
               & rulesRates, &
               & print_it)

!c     Purpose
!c     -------
!c
!c     This subroutine calculates the composition and liquid volume of
!c     stratospheric HNO3-H2SO4-H2O aerosols, i.e. liquid binary aerosols
!c     (LBA, 215 < T <= 240) and supercooled ternary solutions 
!c     (STS/TER, max(Tice-3, 185) < T <= 215).
!c
!c     Results
!c     -------
!c
!c     The results are stored in hno3g, hclg, hbrg, hno3g, hoclg, rho, sad, ... 
!c
!c
!c     Method
!c     ------
!c
!c     This non-iterative method is used for calculation of the 
!c     composition and liquid volume of aqueous HNO3-H2SO4 
!c     stratospheric aerosols is described in detail in 
!c     Carslaw et al. GRL 22 (1995) 1877. The model is valid over
!c     the range of water partial pressure pw (mb)
!c
!c       2e-5 mb < pw < 2e-3 mb
!c
!c     and temperature (T/K)
!c
!c       max(Tice-3, 185) < T < 240.
!        here Tice is ice frost point (~ dew point over ice).
!c
!               ATTENTION
!      That parameterization really fails for low temperatures and high water vapour
!
!c     Liquid binary aerosols (LBA, T > 215K) and supercooled ternary 
!c     solutions (STS/TER, T <= 215K). The partial pressure of water is
!c     assumed not to be affected.
!c
!c     The solubility of a gas phase component in an aerosol can be
!c     calculated from its effective Henry's law constant (H=c/p).
!c     The effective Henry's law coefficient for the ternary solution is
!c     calculated from the ones for the two binary (HNO3-H2O, H2SO4-H2O)
!c     solutions.
!c
!c     The expressions used for calculation of the solubility of HCl and
!c     HBr in stratospheric aerosols have been derived from vapour 
!c     pressure relations given in Luo et al. GRL 22 (1995) 247. 
!c     The Henry's law coefficients are valid over the range of weight
!c     fractions
!c
!c       0 < ms + mn < 0.7
!c
!c     and temperature (T/K)
!c
!c       185 < T < 235.
!c
!c     The expression for calculating the effective Henry's law 
!c     coefficient for HOCl is taken from Huthwelker et al. 
!c     J.Atmos.Chem. 21 (1995) 81. The combined molality of H2SO4
!c     and HNO3 is used.
!c
!c     
!c     Author
!c     ------
!c
!c     FMI developers, 25.11.2002, FMI.
!c
!c-----------------------------------------------------------------------

    implicit none

    ! Imported parameters
    !real, intent(in) :: t_, patm, h2so4_tot, hno3_tot, h2o_tot, hcl_tot, & !hbr_tot, &
    !3Aug2017 by R.H. Added HOCl and HBr total concentrations mol/m3
    real, intent(in) :: Temp, patm, h2so4_tot, hno3_tot, h2o_tot, hcl_tot, hocl_tot, hbr_tot, sulnum, hnm
    real, dimension(:), intent(inout) :: naero   ! aerosol number concentration for each STS bin, m-3
    real, dimension(:,:), intent(inout) :: masses           ! species volumes in each bin of STS aerosols m3/m3
    !3Aug2017 by R.H. Added HOCl and HBr gas concentrations mol/m3
    real, intent(out) :: hno3g, h2og, hclg, hoclg, hbrg !removed the mcl since not needed
    real, intent(out) :: wtH2SO4, extra_sulphur_mass, extra_particles
    real, intent(out) :: rho
    type(Tchem_rules_AerDynMidAtm), intent(in) :: rules
    type(TLocalRates), intent(inout) :: rulesRates
    logical, intent(inout) :: print_it

    ! Local variables
    integer :: i, j, iBin, iTemprCount 
    real :: T, pm, ns, pn0, pcl0, pbr0, pocl0, pw,  ws, wn, v, tt
    real :: pcl,pbr,pocl,wcl,wbr,wocl,xhno3p,xhno3g,delta,fTmp,fLowestTempr
    real, dimension(nBins) :: sad
    real :: lnt,lnpw, tr,pr,lnhnsb,lnhnnb,hnsb,hnnb,xsb,xnb,msb,mnb,ms,mn,wt,rhos,rhon
    real(r8k) :: fA, fB, fC, phi, dTmp, dTmp2, dTmp3

    real :: lnhhcl,lnhhbr,lnhhocl,hhcl,hhbr,hhocl,mcl,mbr,mocl,tm,fcl,fbr,focl  !, aa(2),bb(2)
    real :: h2o_aer, hno3_aer, hcl_aer, hocl_aer, hbr_aer !concentrations of H2O and HNO3 in aerosol
    real :: AA_Cl, BB_Cl, AA_Br, BB_Br, fracFine, fracCoarse
    
    real, parameter ::  RgasAtm = gas_constant_uni/std_pressure_sl ! Gas constant unit R in m3*atm/(K*mol), RgasAtm=8.2057e-5
    ! molar weight of H2O, HCl, HBr, HOCl, H2SO4 and HNO3 (kg mol-1)
    real, parameter :: mwh2o=0.018015,  mwhcl=0.036461,  mwhbr=0.080912, &
                       mwhocl=0.052460, mwh2so4=0.098072, mwhno3=0.063012

    !-----------------------------------------------------------------------
    ! The following data contain the coefficients for calculation of the 
    ! mole fractions and effective Henry's law constants.
    !-----------------------------------------------------------------------
    ! coefficients for mole fracions
    real, dimension(7), parameter :: kn =(/-39.136, 6358.4, 83.29, -17650., 198.53, -11948., -28.469/),&
                                   & ks =(/-21.661, 2724.2, 51.81, -15732., 47.004, -6969.0, -4.6183/)
    ! coefficients for effective Henry's law constants
    real, dimension(10), parameter :: qn = (/ 14.5734, 0.0615994, -1.14895, 0.691693, -0.098863, &
                          & 0.0051579, 0.123472, -0.115574, 0.0110113, 0.0097914/), &
                                    & qs = (/ 14.4700, 0.0638795, -3.29597, 1.778224, -0.223244, &
                          & 0.0086486, 0.536695, -0.335164, 0.0265153, 0.0157550/)
!    !
!    ! coefficients for 'vapour pressure', HCl and HBr
!    !
    real, dimension(8), parameter :: &
                            & a_Cl = (/21.0, 46.61, 4.069, -4.837, 2.186, -63.0, -40.17, -1.571/), &
                            & a_Br = (/17.83, 1.02, -1.08, 3.9, 4.38, -8.87, -17., 3.73/), &
                            & b_Cl = (/-7437.0, -8327.8, 1300.9, 1087.2, -242.71, 18749., 18500., 5632./), &
                            & b_Br = (/-8220.5, -362.76, 658.93, -914., -955.3, 9976.6, 19778.5, 7680./)

#ifdef DEBUG_MORE
if(print_it)then
if(h2so4_tot > 1e-20)then
  call msg('Oxidation started. H2SO4 total and vol(so4,0.1um):', h2so4_tot, masses(iH2SO4_loc,iModeFine))
endif
endif
#endif

    !
    ! If our number concentration is too low, raise it up to the observed value
    !
    naero(iModeFine) = h2so4_tot * rules%mass2vol(iSO4f_aer) / rules%bin_particle_volume(iModeFine)
    if(.not. naero(iModeFine) >= sulnum)then
      extra_sulphur_mass = (sulnum - naero(iModeFine)) * &
                         & rules%bin_particle_volume(iModeFine) / rules%mass2vol(iSO4f_aer)
      extra_particles = sulnum - naero(iModeFine)
#ifdef DEBUG_MORE
if(print_it)then
call msg('Adding volume. H2SO4 total, nAero(0.1um), vol(so4), and extra mass and num:', &
            & (/h2so4_tot, naero(iModeFine), masses(iH2SO4_loc,iModeFine),extra_sulphur_mass, extra_particles/))
endif
#endif
      naero(iModeFine) = sulnum
!      masses(iH2SO4_loc,iModeFine) = masses(iH2SO4_loc,iModeFine) + extra_sulphur_mass
      ns = h2so4_tot + extra_sulphur_mass !/ rules%mass2vol(iSO4f_aer)
    else
#ifdef DEBUG_MORE
if(print_it)then
call msg('Sufficient particles. H2SO4 total, so4, naero(0.1um):', (/h2so4_tot, masses(iH2SO4_loc,iModeFine), naero(iModeFine)/))
endif
#endif
      extra_sulphur_mass = 0.
      extra_particles = 0.
      ns = h2so4_tot                     ! total H2SO4 (mole/m3)
    endif

    !RISTO TEST:
    !write(*,*) 'RHcnc_tot:', h2so4_tot,hno3_tot,h2o_tot,hcl_tot,hocl_tot,hbr_tot
    !write(*,*) 'RHsulnum etc :', sulnum, naero(iModeFine), naero(iModeCoarse), extra_particles, extra_sulphur_mass, ns
    !write(*,*) 'RHvarious :', Temp, patm, rules%mass2vol(iSO4f_aer), rules%bin_particle_volume(iModeFine)
    !END RISTO TEST

    !RH: TESTcnc2vmrChange
    !pn0 = hno3_tot * rulesRates%cnc2vmr * patm ! total HNO3 (atm)
    !pcl0 = hcl_tot * rulesRates%cnc2vmr * patm ! total HCl, HBr and HOCl (atm) 
    pn0 = hno3_tot*RgasAtm*Temp   ! total HNO3 partial pressure (atm) [when nothing is in aerosol]
    pcl0 = hcl_tot*RgasAtm*Temp   ! total HCl partial pressure (atm)
    !3Aug2017 by R.H. Added HOCl and HBr.
    !RH: TESTcnc2vmrChange
    !pbr0 = hbr_tot * rulesRates%cnc2vmr * patm
    !pocl0 = hocl_tot * rulesRates%cnc2vmr * patm
    pbr0 = hbr_tot*RgasAtm*Temp   ! total HBr partial pressure (atm)
    pocl0 = hocl_tot*RgasAtm*Temp ! total HOCl partial pressure (atm)
    !
    ! Are the conditions suitable for the parameterization?
    ! max(Tice-3, 185) < T < 240, Tice is frost point
    ! 2e-5mb < pw < 2e-3mb,       pw is water vapour partial pressure
    !
    !NOTE by R.H. 2017: The water partial pressure (gas phase) is assumed to not to change
    !                   much due to absorption into the aerosol. However, in some situation
    !                   when the air is dry this does not hold and most of the water can be
    !                   taken into the sulfuric-water solution. END NOTE !!!!!!!!!!!!!!!
    !
    pw = h2o_tot*RgasAtm*Temp     ! water partial pressure (atm)

    
!call msg('Forced input concentrations and environment to check against Carslaw')
!call msg('Uncomment below ws and wn, should be 0.262 and 0.226, respectively')
!
!call msg('Before. T,patm,rulesRates%cnc2vmr,hnm,h2o_tot,pw,ns, hno3_tot,pn0:',(/T,patm,rules%cnc2vmr,hnm,h2o_tot,pw,ns, hno3_tot,pn0/))
!
!    T = 192
!    patm = 50.0e-3 / 1.01325        ! 50 mb -> atm
!!    ppa = patm * 101325.            !atm -> Pa
!!    ptorr = ppa / 133.322           !torr
!    rulesRates%cnc2vmr = gas_constant_uni * T / ppa  ! conversion in current conditions
!    hnm = avogadro / rulesRates%cnc2vmr  ! air density, #/m3
!    
!    h2o_tot = 5.e-6 / rulesRates%cnc2vmr
!    pw = 5.e-6 * patm               ! 5 ppmv of H2O -> atm
!    
!    ns = 0.5e-9 / rulesRates%cnc2vmr     ! 0.5 ppbv of H2SO4 -> mole/m3
!    h2so4_tot = ns
!    
!    hno3_tot = 10.0e-9  / rulesRates%cnc2vmr
!    pn0 = 10.e-9 * patm             ! 10 ppbv HNO3 -> atm
!    
!    
!call msg('after. T,patm,rulesRates%cnc2vmr,hnm,h2o_tot,pw,ns, hno3_tot,pn0:',(/T,patm,rulesRates%cnc2vmr,hnm,h2o_tot,pw,ns, hno3_tot,pn0/))
!    


!FIXME iTemprCount and fLowestTempr should be initialised and synchronized somehow....
! Here is just dummi init to make authomatic checks happy
fLowestTempr = 0.
iTemprCount = 0


#ifdef DEBUG_MORE
    if(pw < 2.e-8)then
      if(print_it)call msg('Too low water vapour pressure, extrapolating the Carslaw values, pw in atm:',pw) !Changed 3July2017
      !call msg('Too low water vapour pressure, extrapolating the Carslaw values, pw in atm:',pw) !RISTO TEST 12July2017
      !if(print_it)call msg('Too low water vapour pressure, atm:',pw)
      !pw = 2.e-8 !Commented 3July2017 by R.H.
    elseif(pw > 2.e-6)then
      if(print_it) call msg('Too high water vapour pressure, extrapolating the Carslaw values, pw in atm:',pw) !Changed 3July2017
      !call msg('Too high water vapour pressure, extrapolating the Carslaw values, pw in atm:',pw) !RISTO TEST 12July2017
      !if(print_it) call msg('Too high water vapour pressure, atm:',pw)
      !pw = 2.e-6 !Commented 3July2017 by R.H.
    endif
    if(Temp < 180.0)then
      call msg('Too low temperature, using 180 K for TER parameters. T :',Temp)
      iTemprCount = iTemprCount + 1
      fLowestTempr = min(fLowestTempr, Temp)
      if(mod(iTemprCount,1000) == 0) &
              & call msg('Minimal temperature and number of cases T < 180 K:',fLowestTempr, iTemprCount)
      T = 180.0
    else
      !if (Temp < 185.0) then
      !  call msg('Temperature below 185 K, extrapolating the Carslaw parametrization to temperature:',Temp)
      !end if
      T = Temp
    endif
#else
    !Commented the following 3July2017 by R.H. (extrapolate the Carslaw values)
    !if(pw < 2.e-8)then
    !  !pw = 2.e-8
    !  call msg('Too low water vapour pressure, extrapolating the Carslaw values, pw in atm:',pw) !RISTO TEST 12July2017
    !elseif(pw > 2.e-6)then
    !  !pw = 2.e-6
    !  call msg('Too high water vapour pressure, extrapolating the Carslaw values, pw in atm:',pw) !RISTO TEST 12July2017
    !endif
    if(Temp < 180)then  ! so cold temperature is really rare thing. Should give a warning.
      if(print_it .and. iTemprCount < 100)call msg('Too low temperature:',t_)
      iTemprCount = iTemprCount + 1
      fLowestTempr = min(fLowestTempr, t_)
      if(mod(iTemprCount,1000) == 0) &
              & call msg('Minimal temperature and number of cases T<180K:',fLowestTempr, iTemprCount)
      T = 180.0
    else
      T = Temp
    endif
#endif

    
!patm = 50. / 1035.
!pw = 5.e-6 * patm   ! 5ppmv H2O
!pn0 = 10.e-9 * patm  ! 10ppbv HNO3
!ns = 1.5669e-9  ! 0.5 ppbv H2SO4 -> conc [mol/m3]
!t = 192.

    lnt = log(T)
    lnpw = log(pw)
    tr = 1.e4*(1./T-1./230.)
    pr = lnpw+18.4

    ! mole fraction of H2SO4 in binary solution
    ! Cut out the low temperature with high humidity range: the xsb is just zero over there.
    ! For very warm temperatures, the ratio can be quite high reaching 1 way outside the applicability
    ! of this formula - around 290K, dry atmosphere. Still, let's be on the safe side.
    !
    ! Attention!!
    ! Carslaw et al have overlooked reduction of equation from quadratic to linear. 
    ! For H2SO4 it is 303K, hardly a problem but still...
    !
    if(abs(T + ks(4)/ks(3)) < 0.1)then
      xsb = -(ks(5) + ks(6)/T + ks(7)*lnt - lnpw) / (ks(1) + ks(2)/T)
    else
      fTmp = (ks(1)+ks(2)/T)**2 - 4.*(ks(3)+ks(4)/T)*(ks(5)+ks(6)/T+ks(7)*lnt-lnpw)
      !
      ! Note that in approaching the cold-humid conditions, ratio xsb/xnb decreases from 0.6 
      ! and goes sharply down below zero when xnb flips the sign. So, we shall not let it by min
      !
      if(fTmp > 0.)then
        xsb = min(0.99, max(0.003, (-ks(1)-ks(2)/T-sqrt(fTmp)) / (2.*(ks(3)+ks(4)/T)) ))
      else
        xsb = min(0.99, max(0.003, (-ks(1)-ks(2)/T)/(2.*(ks(3)+ks(4)/T)) ))
      endif
    endif

    msb = 1./mwh2o*xsb/(1.-xsb) ! molality of H2SO4 in binary solution (mol kg(H2O)-1)

    ! molality of H2SO4 and HNO3 (mol kg(H2O)-1) in the aerosol
    ! liquid binary aerosols if T>215K otherwise supercooled ternary solutions
    !
    if ( t >= 215. .or. pn0 <= 0. ) then   ! liquid binary solution (LBI)
       ms = msb             ! mass in LBI is equal to molality (mole H2SO4 / kg H20)
       mn = 0.              ! HNO3 is not there
       xnb = 0.
       mnb = 0.
       tt = real_missing
    else
      !
      ! effective Henry's law coefficients in binary solutions, hnsb and hnnb
      ! (mol kg(H2O)-1 atm-1)
      !
      lnhnsb = qs(1)+qs(2)*tr*tr+(qs(3)+qs(4)*tr+qs(5)*tr*tr &
            &  +qs(6)*tr*tr*tr)*pr+(qs(7)+qs(8)*tr+qs(9)*tr*tr)*pr*pr &
            &  +qs(10)*tr*pr*pr*pr
      hnsb = exp(lnhnsb)
      lnhnnb = qn(1)+qn(2)*tr*tr+(qn(3)+qn(4)*tr+qn(5)*tr*tr &
             & +qn(6)*tr*tr*tr)*pr+(qn(7)+qn(8)*tr+qn(9)*tr*tr)*pr*pr &
             & +qn(10)*tr*pr*pr*pr
      hnnb = exp(lnhnnb)

      ! mole fraction of HNO3 in binary solutions. 
      !
      ! Attention!!
      ! Carslaw et al have overlooked reduction of equation from quadratic to linear. 
      ! It happens at 211.91K, have to take care of it.
      !
      dTmp2=0. ! Reset it to something reasonable
!$OMP CRITICAL(debug_ter)
      if(abs(t + kn(4)/kn(3)) < 0.1)then
        xnb = -(kn(5) + kn(6)/T + kn(7)*lnt - lnpw) / (kn(1) + kn(2)/T)
      else
        fTmp = (kn(1)+kn(2)/T)**2 - 4.*(kn(3)+kn(4)/T)*(kn(5)+kn(6)/T+kn(7)*lnt-lnpw)
        if(fTmp > 0.)then
          xnb = (-kn(1)-kn(2)/T-sqrt(fTmp)) / (2.*(kn(3)+kn(4)/T))
        else
          xnb = -1.
        endif
      endif
      !
      ! Mandatory cut because Carslaw parameterization fails in two corners
      !      
      xnb = min(0.99, max(0.01,xnb))

!      if(.not. (xnb > 0. .and. xnb <= 1.0))then
!        call msg('xnb problem. xnb, fTmp, lnt, lnpw, T, pw',(/xnb, real(fTmp), lnt, lnpw, T, pw/))
!        call msg('(hnnb, hnsb, pn0)',(/hnnb, hnsb, pn0/))
!        call msg('(-kn(1)-kn(2)/t-sqrt(fTmp)) / (2.*(kn(3)+kn(4)/t))', (-kn(1)-kn(2)/t-sqrt(fTmp)) / (2.*(kn(3)+kn(4)/t)))
!        call msg('(-kn(1)-kn(2)/t-sqrt(fTmp)) ', (-kn(1)-kn(2)/t-sqrt(fTmp)) )
!        call msg('(kn(3), kn(4), t, 2.*(kn(3)+kn(4)/t))',(/kn(3), kn(4), t,  (2.*(kn(3)+kn(4)/t)) /) )
!        call msg('kn',kn)
!        print_it = .true.
!      endif

      mnb = 1./mwh2o*xnb/(1.-xnb) ! molality of HNO3 in binary solutions (mol kg(H2O)-1)
      !
      ! Another problem outlooked by Carslaw et al: large coeffiicents and/or mnb ~ msb. In this case, 
      ! fA, fB, fC all become large and their cubes etc become infinitely large and their subtraction 
      ! stops making sense.
      !
      if(abs(mnb - msb) < 1.e-4 * (mnb+msb))then
        ms = 0.0001 * msb                ! not actually sure but seems so
        mn = mnb
      else
        !
        ! functions needed for relating the molalities of H2SO4 and HNO3 in the 
        ! binary solutions to the ternary solution
        !
        tt = RgasAtm*T*ns   ! ns [concentration] -> atm., R = 8.205e-5 m3 atm mol-1 K-1
        fA = (tt*hnnb*mnb-tt*hnsb*msb-2.*mnb*msb+ msb*msb+hnnb*msb*pn0-hnsb*msb/mnb*msb*pn0)/(mnb-msb)
        fB = msb*(-2.*tt*hnnb*mnb+tt*hnsb*msb+mnb*msb-hnnb*msb*pn0) / (mnb-msb)
        fC = (tt*hnnb*mnb*msb*msb) /(mnb-msb)

        dTmp = -2.0_r8k*fA*fA*fA + 9.0_r8k*fA*fB - 27.0_r8k*fC
        dTmp2 = fA*fA-3.0_r8k*fB
        dTmp3 = dTmp2 / dTmp

        if(.not. (4.* dTmp3*dTmp3 * dTmp2 - 1.0_r8k) >= 0.) then
          phi=0.0_r8k
          call msg('phi undefined: sqrt(negative),T:',real(4.* dTmp3*dTmp3 * dTmp2 - 1.0_r8k), t)
          call msg('(hnnb, hnsb, pn0)',(/hnnb, hnsb, pn0/))
          call msg('T,tt,ns,fA,fB,fC,()_for_cube,()_for_square,()_ratio',(/T,tt,ns,real(fA),real(fB),real(fC),real(dTmp),real(dTmp2),real(dTmp3)/))
          call msg('mnb, msb,xnb, lnt, lnpw, T, pw',(/mnb,msb,xnb, lnt, lnpw, T, pw/))
          print_it = .true.
        else
          phi = datan(dsqrt(4.*dTmp3*dTmp3 * dTmp2 - 1.0_r8k)*sign(1.0_r8k,dTmp) )
        endif
        if(phi < 0.) phi = phi + PI
        if (dTmp2 >= 0) then
           ms = -1./3. * (fA+2.*SQRT(dTmp2)*COS((PI+phi)/3.) )
        else
           call msg("Negative ()_for_square. Return to save the run... ")
           !ms= msb !FIXME Dirty hack: attempt to make global suite of v5_5 pass
           ! Helps for a single timestep
           ! ms = 0. !Zero here causes further problems
        endif
        !
        ! Here we have to take care of numerics: third-order equation solved via trigonometrics functions
        ! is never a good thing if someone is after mass budget.
        !
        if(ms <= 0.)then  ! alarm if sulphur fraction goes to zero. Should never happen
          call msg('Negative ms. H2SO4 and ms:',ns,ms)
          ms = 0.   !ms == 0 leads to infinite volume
          mn = mnb
        elseif(ms >= 0.999 * msb)then  ! H2SO4 concentration is higher than or ~ at binary equilibrium without HNO3
          mn = 0.
        else
          mn = mnb/msb*(msb-ms)
        endif
      endif  ! mnb~msb
      if (.not. dTmp2 >=0. ) call msg("Trouble with TER. Will return")
!$OMP END CRITICAL(debug_ter)
      if (.not. dTmp2 >=0. ) return 
      !mn = max(mnb*(1.-ms/msb), 0.)  ! just to be on the safe side...
      ! equilibrium pressure of HNO3 (atm)
      !write(*,*) mn,ms,ns,hnsb,hnnb,(hnnb*mn/(mn+ms)+hnsb*ms/(mn+ms))
    endif  ! Liquid binary solution of upercooled ternary

    ! weight fraction of H2SO4, HNO3
    !wt = 1.+ms*mwh2so4+mn*mwhno3    ! total, kg  !FIXME Nonsense! 1kg???? 
    wt = 1.+ms*mwh2so4+mn*mwhno3    ! total mass of aerosol per 1kg of water 
    ws = ms*mwh2so4/wt              ! sulphuric acid mass fraction
    wn = mn*mwhno3/wt               ! nitric acid mass fraction

    !
    ! call msg('Carlslaw test: should be ws=0.262, wn=0.226. In reality:',ws,wn)
    !    
    
    if(.not. ws > 0. .or. .not. wn >= 0)then
       call msg('tt, tr, T, pw, pr', (/tt, tr, T, pw, pr/))
       call msg('ws, xsb, msb, lnhnsb, hnsb', (/ws, xsb, msb, lnhnsb, hnsb/))
       call msg('wn, xnb, mnb, lnhnnb, hnnb', (/wn, xnb, mnb, lnhnnb, hnnb/))
       call msg('phi, fA, fB, fC, ms, mn,  wt', (/real(phi), real(fA), real(fB), real(fC), ms, mn,  wt/))
       call set_error('Strange sulphur/nitrogen fractions in STS','TER')
       return
    endif
    
    !July2017 by R.H.: In order to avoid negative amount of gaseous H2O/HNO3 we limit to amounts of H2O/HNO3 in aerosol.
    h2o_aer = ns*(max((1.-ws-wn),0.0)/ws)*(mwh2so4/mwh2o) !max to avoid numerical rounding errors
    hno3_aer = ns*(wn/ws)*(mwh2so4/mwhno3)
    if ( h2o_aer > h2o_tot .or. hno3_aer > hno3_tot ) then
       !Commenting the calls to msg (10 x) inside this if structure, since for some runs this produces lots of output! 
       !call msg('Additional Hack by R.H. Likely due to extremely dry air, pw only (in atm):',pw)
       !call msg('Temperature (K), pressure (atm), rulesRates%cnc2vmr: ', (/T, patm, rulesRates%cnc2vmr/))
       !call msg('h2so4_tot, ns: ', (/h2so4_tot, ns/)) 
       !call msg('h2o_tot, h2o_aer, pw (atm): ', (/h2o_tot, h2o_aer, pw/))
       !call msg('hno3_tot, hno3_aer, mn*h2o_aer*mwh2o: ', (/hno3_tot, hno3_aer, mn*h2o_aer*mwh2o/))
       !call msg('before the hack: ws, wn: ', (/ws, wn/))
       if (h2o_aer >= h2o_tot) then
          !call msg('Limiting H2O in aerosol to total amount of H2O. (h2o_tot,h2o_aer):', (/h2o_tot, h2o_aer/))
          h2o_aer = min(h2o_aer,0.99*h2o_tot)   
       end if
       if (hno3_aer > hno3_tot) then
          !call msg('Limiting HNO3 in aerosol to total amount of HNO3. (hno3_tot,hno3_aer):', (/hno3_tot, hno3_aer/))
          hno3_aer = min(hno3_aer,0.99*hno3_tot)
       end if
       if (hno3_aer > mn*h2o_aer*mwh2o) then
          !call msg('Limiting HNO3 in aerosol to value determined by molality. (mn*h2o_aer*mwh2o,hno3_aer):', (/mn*h2o_aer*mwh2o, hno3_aer/))
          hno3_aer = min(hno3_aer,mn*h2o_aer*mwh2o)
       end if
       if (hno3_aer < 1e-30) then
          hno3_aer = 0.0
       end if
       ws = ns*mwh2so4/(h2o_aer*mwh2o+ns*mwh2so4+hno3_aer*mwhno3)
       wn = hno3_aer*mwhno3/(h2o_aer*mwh2o+ns*mwh2so4+hno3_aer*mwhno3)
       !call msg('after the hack: ws, wn: ', (/ws, wn/))
    end if
    !End July2017 by R.H.:

    wtH2SO4 = ws * 100.   ! mass fraction in %, if someone needs it, Changed 3July2017 by R.H.

#ifdef DEBUG_MORE
if(print_it)then
  call msg('TER: OK, wtH2SO4:', wtH2SO4)
endif
#endif
    ! density of the binary solutions (kg m-3)
    rhos = 1000.+(123.64-5.6e-4*T*T)*ms-(29.54-1.81e-4*T*T)*ms**1.5 &
         & +(2.343-1.487e-3*T-1.324e-5*T*T)*ms*ms

    rhon = 1000.+(85.11-5.04e-4*T*T)*mn-(18.96-1.427e-4*T*T)*mn**1.5 &
         & +(1.458-1.198e-3*T-9.703e-6*T*T)*mn*mn

    !NOTE: The HCl, HOCl, and HBr, added below, are assumed not to change the density of the solution!
    rho = 1./((1./rhos)*ms/(mn+ms)+(1./rhon)*mn/(mn+ms)) ! density of the ternary solution (kg m-3)

    !v = ns/(gas_constant_uni*t)*mwh2so4/(ws*rho)  ! volume of aerosol, m3(aer) per m3 of air
    !v = ns * rules%mass2vol(iSO4f_aer) !/ ws      ! volume of sulphur aerosol, m3(aer) per m3 of air
    v = ns*rules%mwh2so4/(ws*rho)        ! H2SO4 cnc mole/m3 -> volume of total aerosol, m3aer / m3air

    ! Next dirty hack....
    if (.not. v>0. .or. .not. ws>0. .or. .not. rho>0.) then
        call msg("OOOps v , ns , rules%mwh2so4, ws, rho ", (/v, ns, rules%mwh2so4, ws,rho/) )
        call msg("rhos rhon t ms mn", (/rhos,rhon,t,ms,mn/))
        write (unit=strTmp, fmt='(20(e10.3,2x))') lnt, lnpw, XNB, XSB, MNB, MSB, tr, pr, HNSB, HNNB, TT, fA, fB, fC, PHI, MS, MN,  WS, WN
        call msg('SILAM code:')
        call msg('logt,       logpw,        XNB,         XSB,         MNB,        MSB,        tr,        pr,        HNSB,        HNNB' + &
                   & ',          TT,        fA,         fB,         fC,        PHI,        MS,        MN,         WS,        WN')
        call msg(strTmp)
        return
    endif
#ifdef DEBUG_MORE
if(print_it)then
write (unit=strTmp, fmt='(20(e10.3,2x))') lnt, lnpw, XNB, XSB, MNB, MSB, tr, pr, HNSB, HNNB, TT, fA, fB, fC, PHI, MS, MN, WS, WN
call msg('SILAM code:')
call msg('logt,       logpw,        XNB,         XSB,         MNB,        MSB,        tr,        pr,        HNSB,        HNNB' + &
           & ',          TT,        fA,         fB,         fC,        PHI,        MS,        MN,         WS,        WN')
call msg(strTmp)
endif
#endif

    !
    ! Now the volume v has to be distributed proportionally to the surface area of all bins.
    ! In principle, here must be condensation process towards this equilibrium.
    ! Also, here we take into account that so far we do not have all sulphur at place, i.e. have to
    ! add some funny mass to get it to the observed values written in sulnum.
    ! Note that HNO3 and water are, hopefully, in right order: HNO3 is produced locally,
    ! water is transported from the troposphere.
    ! Note also that HCl, HOCl, and HBr concentrations are not taken into account when
    ! calculating the SAD's because they fraction is always assumed to be small.
    !
    sad(1:nBins) = 0.    ! calculate total area before taking the new equilibrium into account
    !RISTO Cl&Br Q: Should one update these after HCl, HOCl, and HBr are added? Likely not since their effect is considered to be small??? 
    do iBin = 1, nBins
      fTmp  = (masses(iH2SO4_loc,iBin)*rules%mwH2SO4 + masses(iHNO3_loc,iBin)*rules%mwHNO3 + &
             & masses(iH2O_loc,iBin)*rules%mwH2O) / rho !Total volume of bin
      sad(iBin) = sad(iBin) + (Pi_36*naero(iBin) * fTmp * fTmp ) ** 0.3333333333333333333
    end do
        
    !
    ! 12.06.2023. Faced ws + wn > 1 but very close to it. Since these are fractions, correct.
    ! Potential problem with water amount in aerosol - I put it to zero here. MAS
    !
    if(ws + wn > 1.0)then
      if(ws + wn < 1.0001)then
        fTmp = 1.0 / (ws + wn)
        ws = ws * fTmp
        wn = wn * fTmp
      else
      !$OMP CRITICAL (BARK_TER)
        call msg_warning('ws + wn > 1.0 substantially','TER')
        call msg('v,sad(1:nBins),ws,wn,T,ms,msb,mn,mnb',(/v,sad(1:nBins),ws,wn,t,ms,msb,mn,mnb/))
        write(unit=strTmp,fmt='(1x,1000(E9.3,2x))')masses
        call msg(' Metamasses:' + strMetamasses)
        call msg(' masses    :' + strTmp)
        write(unit=strTmp,fmt='(1x,1000(E9.3,2x))')naero
        call msg('Metanaero:' + strMetaNaero)
        call msg('naero    :' + strTmp)
        call set_error('ws + wn > 1.0','TER')
        !$OMP END CRITICAL (BARK_TER)
      endif
    endif
    !
    ! Distribute H2SO4, HNO3, H2O in the aerosol phase
    !




    if(sum(sad(1:nBins)) < 1.e-20)then
      !
      ! No existing STS, these ones are new, naero forced, formed on existing sulphates and other nucleation cores. 
      ! Put everything into the first bin
      !
#ifdef DEBUG_MORE
if(print_it)then
write(unit=strTmp,fmt='(1x,1000(E9.3,2x))')masses
call msg('Before New STS: Metamasses:' + strMetamasses)
call msg('Before New STS: masses    :' + strTmp)
call msg('v,ws,wn,mwh2so4,mwhno3,mwh2o:',(/v,ws,wn,mwh2so4,mwhno3,mwh2o/))
endif
#endif
      ! Spread the total volume according to molar fractions
      masses(iH2SO4_loc,iModeFine) = ns           ! we say that all H2SO4 is in aerosol phase
      !masses(iHNO3_loc,iModeFine) = ns / ws * wn   !/mwhno3
      !masses(iH2O_loc,iModeFine) = ns / ws * (1.-ws-wn)  !/mwh2o
      !Changed to 3July2017 by R.H.: These are the correct equations
      masses(iHNO3_loc,iModeFine) = ns*(wn/ws)*(mwh2so4/mwhno3)
      masses(iH2O_loc,iModeFine) = ns*(max((1.-ws-wn),0.0)/ws)*(mwh2so4/mwh2o) !max to avoid numerical rounding errors
#ifdef DEBUG_MORE
if(print_it)then
write(unit=strTmp,fmt='(1x,1000(E9.3,2x))')masses
call msg('New STS: Metamasses:' + strMetamasses)
call msg('New STS: masses    :' + strTmp)
endif
#endif
    else
      !
      ! STS already available in the cell. The new mass is ns for H2SO4, which is projected to HNO3 
      ! and H2O via ws and wn. The above sad is the scaling parameter, i.e. relaxation
      ! should be computed via condensation/evaporation - in principle. Here we take it instant.
      ! The excess/deficit distribution goes proportionally to surface area of each bin
      !
      delta = (ns-sum(masses(iH2SO4_loc,1:nBins))) / sum(sad(1:nBins)) ! extra mass of H2SO4, normalised
      delta = max(delta,0.0) !Avoid numerical rounding errors
      do iBin = 1, nBins
        masses(iH2SO4_loc,iBin) = masses(iH2SO4_loc,iBin) + delta * sad(iBin)
        !masses(iHNO3_loc,iBin) = masses(iHNO3_loc,iBin) + delta * sad(iBin) * wn / ws
        !masses(iH2O_loc,iBin) = masses(iH2O_loc,iBin) + delta * sad(iBin) * (1.-ws-wn) / ws
        !Changed to 3July2017 by R.H.: These are the correct equations
        masses(iHNO3_loc,iBin) = masses(iHNO3_loc,iBin) + delta*sad(iBin)*(wn/ws)*(mwh2so4/mwhno3)
        masses(iH2O_loc,iBin) = masses(iH2O_loc,iBin) + delta*sad(iBin)*(max((1.-ws-wn),0.0)/ws)*(mwh2so4/mwh2o) !max to avoid numerical rounding errors
      end do

#ifdef DEBUG_MORE
if(print_it)then
write(unit=strTmp,fmt='(1x,1000(E9.3,2x))')masses
call msg('Existing STS: Metamasses:' + strMetamasses)
call msg('Existing STS: masses    :' + strTmp)
endif
#endif

    endif   ! new or old STS
    
    if(.not. all(masses(:,:) >= 0.))then
      !$OMP CRITICAL (BARK_TER)
      call msg_warning('Negative volume','TER')
      call msg('v,delta,sad(1:nBins),ws,wn,T,ms,msb,mn,mnb',(/v,delta,sad(1:nBins),ws,wn,t,ms,msb,mn,mnb/))
      write(unit=strTmp,fmt='(1x,1000(E9.3,2x))')masses
      call msg(' Metamasses:' + strMetamasses)
      call msg(' masses    :' + strTmp)
      write(unit=strTmp,fmt='(1x,1000(E9.3,2x))')naero
      call msg('Metanaero:' + strMetaNaero)
      call msg('naero    :' + strTmp)
      call set_error('Negative volume','TER')
      !$OMP END CRITICAL (BARK_TER)
      return
    endif

#ifdef DEBUG_MORE
if(print_it)then
write(unit=strTmp,fmt='(1x,1000(E9.3,2x))')naero
call msg('Metanaero:' + strMetaNaero)
call msg('naero    :' + strTmp)
endif
#endif

    !
    ! Partioning of HNO3. Above, the aerosol amount has been calculated. Its connection with
    ! gas phase is mass2vol - at least in principle. THIS MAY BE WRONG.
    ! Enforce it here
    !
    fTmp = sum(masses(iHNO3_loc,1:nBins))  ! amount of HNO3 in aerosols, mole
    if (fTmp > hno3_tot) then
      call msg('Amount in aerosols is bigger than total HNO3:', real(fTmp), real(HNO3_tot))
      call set_error('Amount in aerosols is bigger than HNO3_tot','TER')
      !write(*,*) 'h2so4_tot, ns:', h2so4_tot, ns !RISTO tmp
      !write(*,*) 'h2o_tot: ', h2o_tot !RISTO tmp
      !write(*,*) 'T and patm, pw (atm): ', t_, patm, pw !RISTO tmp 
      return
    else
      hno3g = hno3_tot - fTmp
    endif

    !
    ! Partioning of H2O. Above, the aerosol amount has been calculated. Its connection with
    ! gas phase is just the mass budget. Enforce it here
    !
    fTmp = sum(masses(iH2O_loc,1:nBins))  ! amount of H2O in aerosols, mole
    if (fTmp > h2o_tot) then
      call msg('Amount in aerosols is bigger than total H2O:', real(fTmp), real(H2O_tot))
      !call set_error('Amount in aerosols is bigger than H2O_tot','TER')
      !return
      !Dirty hack for avoiding this 3Jan2019 (appears in puhti)
      if (masses(iH2O_loc,iModeFine) > h2o_tot) then
         masses(iH2O_loc,iModeFine) = h2o_tot
         masses(iH2O_loc,iModeCoarse) = 0.0
      else
         masses(iH2O_loc,iModeCoarse) = h2o_tot-masses(iH2O_loc,iModeFine)
      end if
      h2og = 0.0
    else
      h2og = h2o_tot - fTmp
    endif


    ! ----------------------------------------------------------------------
    !
    ! calculation of solubility of HCl, HBr and HOCl
    !
    !return ! skip partition of hcl, hbr and hocl for now
    
    !Re-calculate the SAD's in order to be able to distribute the HCl, HBr, and HOCl 
    sad(1:nBins) = 0.
    do iBin = 1, nBins
      fTmp  = (masses(iH2SO4_loc,iBin)*rules%mwH2SO4 + masses(iHNO3_loc,iBin)*rules%mwHNO3 + &
             & masses(iH2O_loc,iBin)*rules%mwH2O) / rho
      sad(iBin) = sad(iBin) + (Pi_36*naero(iBin) * fTmp * fTmp ) ** 0.3333333333333333333
    end do

    ! define weight fraction, equilibrium pressure, molality and 
    ! partitioning into liquid phase of HCl, HBr and HOCl
    ! The equations are from Carslaw et al. GRL 22, 1877 (1995) and Luo et al. GRL 22, 247 (1995).
    !
    ! If no aerosol particles then HCl, HBr, and HOCl all remain in gas phase.
    ! Also, if conditions are outside the model validity range, or if ms <= 0 (should not occur) they remain in gas phase.
    if ( (sum(sad(1:nBins)) < 1.0e-20) .or. (ws+wn > 0.7) .or. (t > 235.) .or. (ms <= 0) )then
       !wcl = 0.0; wbr = 0.0; wocl = 0.0; pcl = pcl0; pbr = pbr0; pocl = pocl0; mcl = 0.0; mbr = 0.0; mocl = 0.0
       masses(iHCl_loc,1:nBins) = 0.0
       masses(iHBr_loc,1:nBins) = 0.0
       masses(iHOCl_loc,1:nBins) = 0.0
       hclg  = hcl_tot
       hbrg  = hbr_tot
       hoclg = hocl_tot
    else
       ! Conditions are fine, model can be used. Coefficients for 'vapour pressure'
       AA_Cl=a_Cl(1)+a_Cl(2)*wn+a_Cl(3)*ws+a_Cl(4)*sqrt(wn)+a_Cl(5)*sqrt(ws)+a_Cl(6)*wn*wn+a_Cl(7)*wn*ws+a_Cl(8)*ws*ws
       BB_Cl=b_Cl(1)+b_Cl(2)*wn+b_Cl(3)*ws+b_Cl(4)*sqrt(wn)+b_Cl(5)*sqrt(ws)+b_Cl(6)*wn*wn+b_Cl(7)*wn*ws+b_Cl(8)*ws*ws
       AA_Br=a_Br(1)+a_Br(2)*wn+a_Br(3)*ws+a_Br(4)*sqrt(wn)+a_Br(5)*sqrt(ws)+a_Br(6)*wn*wn+a_Br(7)*wn*ws+a_Br(8)*ws*ws
       BB_Br=b_Br(1)+b_Br(2)*wn+b_Br(3)*ws+b_Br(4)*sqrt(wn)+b_Br(5)*sqrt(ws)+b_Br(6)*wn*wn+b_Br(7)*wn*ws+b_Br(8)*ws*ws
       
       ! Henry's law coefficient for HCl and HBr (mol kg(aerosol)-1 atm-1)
       ! Assuming that the weight fractions of HCl and HBr are small compared to H2NO4 and HNO3 (this implies wcl/mcl = mwhcl etc)
       lnhhcl = -(AA_Cl + BB_Cl/T + log(mwhcl*(wn + 0.61*ws))) !Fine when wcl << wn,ws
       lnhhbr = -(AA_Br + BB_Br/T + log(mwhbr*(wn + 0.41*ws))) !Fine when wbr << wn,ws
       hhcl = 1013.25*exp(lnhhcl) !in mol kg-1, atm-1, also the fractor 1013.25 is due to the fact that Luo et al, GRL1995 give pressure in mbar
       hhbr = 1013.25*exp(lnhhbr)
       
       ! Henry's law coefficient for HOCl (mol kg(aerosol)-1 atm-1) using the combined molalities of H2SO4 and HNO3
       lnhhocl = 6.4946+(mn+ms)*(0.04107-54.56/T)-5862*(1/298.15-1/T)
       hhocl = exp(lnhhocl)
       
       ! molality of HCl, HBr and HOCl (mol kg(aerosol)-1)
       !tm = ns/(ms*(1-ws-wn)) !RISTO: Why this was used earlier? The units are already wrong!!!!
       tm = RgasAtm*T*ns/ms   !R = 8.205e-5 m3 atm mol-1 K-1
       mcl = pcl0/(tm+1/hhcl)
       mbr = pbr0/(tm+1/hhbr)
       mocl = pocl0/(tm+1/hhocl)
       ! weight fraction of HCl, HBr and HOCl (here also assuming that the fractions are small compared with ws and wn)
       wcl = mcl*mwhcl
       wbr = mbr*mwhbr
       wocl = mocl*mwhocl
       ! equilibrium pressure of HCl, HBr and HOCl (atm)
       !pcl = pcl0-tm*mcl
       !pbr = pbr0-tm*mbr
       !pocl = pocl0-tm*mocl
       
       !concentrations of gaseous hcl, hbr, and hocl
       !hclg  = pcl/(patm*rulesRates%cnc2vmr)
       !hbrg  = pbr/(patm*rulesRates%cnc2vmr)
       !hoclg = pocl/(patm*rulesRates%cnc2vmr)

       !Try to avoid using the pressure, since at small concentration when p0=p it might result e.g. hbrg /= hbr_tot since
       !multiplification and division by the scaling factor does not always result the original number if the number is close to zero.
       !hcl_aer  = min(tm*mcl/(RgasAtm*T),hcl_tot)
       !hbr_aer  = min(tm*mbr/(RgasAtm*T),hbr_tot)
       !hocl_aer = min(tm*mocl/(RgasAtm*T),hocl_tot)
       !setting tm in place the above simplifies to
       hcl_aer  = min(ns*mcl/ms,hcl_tot)
       hbr_aer  = min(ns*mbr/ms,hbr_tot)
       hocl_aer = min(ns*mocl/ms,hocl_tot)
       hclg     = max(hcl_tot-hcl_aer,0.0)
       hbrg     = max(hbr_tot-hbr_aer,0.0)
       hoclg     = max(hocl_tot-hocl_aer,0.0)
       fracFine = sad(iModeFine)/sum(sad(1:nBins))
       fracCoarse = sad(iModeCoarse)/sum(sad(1:nBins))
       masses(iHCl_loc,iModeFine)    = hcl_aer*fracFine  !Amount of HCl in aerosols
       masses(iHCl_loc,iModeCoarse)  = hcl_aer*fracCoarse
       masses(iHBr_loc,iModeFine)    = hbr_aer*fracFine  !Amount of HBr in aerosols
       masses(iHBr_loc,iModeCoarse)  = hbr_aer*fracCoarse
       masses(iHOCl_loc,iModeFine)   = hocl_aer*fracFine !Amount of HOCl in aerosols
       masses(iHOCl_loc,iModeCoarse) = hocl_aer*fracCoarse 
       
       if (wcl > 0.10 .or. wcl > 0.3*(wn+0.61*ws)) then
          call msg('*** WARNING: Assumption in TER that HCl amount is small is not valid! wcl, ws, wn: ', (/wcl, ws, wn/))
       end if
       if (wocl > 0.10 .or. wocl > 0.3*(wn+0.41*ws)) then
          call msg('*** WARNING: Assumption in TER that HOCl amount is small is not valid! wcol, ws, wn: ', (/wocl, ws, wn/))
       end if
       !if (wbr > 0.10 .or. wbr > 0.3*(wn+ws)) then
       if (wbr > 0.25 .or. wbr > 0.5*(wn+ws)) then !Since HBr amount can be rather large we only print the warning when assumption is really bad!
          call msg('*** WARNING: Assumption in TER that HBr amount is small is not valid! wbr, ws, wn: ', (/wbr, ws, wn/))
       end if
       
    endif ! partition of hcl, hbr and hocl: if inside the model applicability range and there are aerosol particles present. 
    
    !RISTO TEST 20April2018 outputting the gas phase fraction H2O, HNO3, HCl, HOCl, and HBr, compared with total amount 
    !write(*,'(A,5F10.6)') 'RISTO TEST TERFRAC:', h2og/h2o_tot, hno3g/hno3_tot, hclg/hcl_tot, hoclg/hocl_tot, hbrg/hbr_tot
    !write(*,'(A,5E14.6)') 'RISTO TEST TER GAS:', h2og, hno3g, hclg, hoclg, hbrg
    !write(*,'(A,5E14.6)') 'RISTO TEST TER AER:', h2o_aer, hno3_aer, hcl_aer, hocl_aer, hbr_aer
    !write(*,'(A,5E14.6)') 'RISTO TEST TER TOT:', h2o_tot, hno3_tot, hcl_tot, hocl_tot, hbr_tot
    !END RISTO TEST 20April2018 outputting the gas phase fraction H2O, HNO3, HCl, HOCl, and HBr, compared with total amount 
    
    return
  end subroutine TER


  !*********************************************************************************

  subroutine NAT(hno3eq, &   !  (r) = equil. HNO3 concentration above NAT (unchanged), cm-3
               & hno3g, &   !   (r) = HNO3 concentration in gas phase (changed), mole/m3
               & h2og, &    !   (r) = H2O concentration in gas phase (changed), mole/m3
               & masses, &    !   (r) = # density of hno3 in particle (changed), cm-3
               & natnum, &   !  (r) = # density of NAT particles (changed), cm-3 
!              & sad, &   !     (r) = sad of NAT (changed), um2 cm-3
 !              & cnc2vol, &  !  (r) = scaling factor from concentrations to volume
               & rules, rulesRates, density, seconds, &
               & print_it)
!c
!c     This subroutine calculates
!c       - the partition of HNO3 between gas phase and  NAT
!c       - the radius of the NAT particles
!c       - the sad of NAT
!c
!c     The results are stored in hno3g, hno3c, rpn and sad.     
!c
!c     PSC type 1a, nitric acid trihydrate (NAT).
!c     The [thermodynamic equilibrium of NAT] = f(T,pH2O,pHNO3) 
!c     is used, see Hanson and Mauersberger, GRL, 15 (1988) 855.  
!c     The NAT rock # density is based on Fahey et al., Science, 
!c     291 (2001) 1026.
!c
!c     Author:     FMI developers, 25.11.2002, FMI.
!      Adaptation to SILAM environment & dynamic features: M.Sofiev
!
!      Here we deviate from FinROSE: fixing number of the particles looks wrong if concentration is 
!      not constrained in any sense. Therefore, we fix the size of newly formed NAT particles. This
!      is taken to be 0.1 um (diam), which is the size of background sulphur aerosol in the stratosphere. 
!      An asumption is then that the NAT aerosol is formed on pre-existing sulphate aerosol, which is
!      about that size (also mentioned by Fahey et al).
!c
!c-----------------------------------------------------------------------

    implicit none  

    ! Imported parameters
    real, intent(in) :: hno3eq
    real,intent(inout) :: hno3g, h2og
    real, dimension(:), intent(inout) :: natnum !, cnc2vol
    real, dimension(:,:), intent(inout) :: masses
!    real, intent(out) :: sad
    type(Tchem_rules_AerDynMidAtm), intent(in) :: rules
    type(TLocalRates), intent(inout) :: rulesRates
    real, intent(in) :: seconds
    real, intent(out) :: density
    logical, intent(in) :: print_it

    ! Local variables
    real csNAT,v, mass_tot, n_tot, fTmp
    integer :: iBin
    ! rho kg m(nat)-3, mw kg mol(nat)-3
    real, parameter :: rhonat=1350.  !NOTE: This should be obtained from silam_chemicals.dat


    density = rhonat
    !
    !The number of nat particles (natnum) is currently calculated at the end of the subroutine. 
    !
    ! Is it new NAT or existing for some time already?
    !
    mass_tot = sum(masses(iNAT_loc,1:nBins))
    n_tot = sum(natnum(1:nBins))
    if(n_tot > 1.0e-3)then         ! 1 particle per 1000 m3 is not much
      !
      ! This is old NAT, take care of the distribution. 
      ! Actually, just compute the excess/deficit and attribute it to the smallest bin
      ! Should the deficit be too large and the smallest bin fully depleted, take the next one
      !
      csNAT = hno3g - hno3eq                       ! condensed phase concentration, mole/m3
      if(csNAT > 0.)then            ! Condensation
        if (3.*csNAT > h2og) then   ! Gaseous water is the limiting factor 
          csNAT = h2og / 3.         ! this is Nitric Acid TRIHYDRAT
          h2og = 0.
        else
          h2og = h2og - csNAT * 3.
        endif
        hno3g = hno3g - csNAT                        ! gas phase concentration
        masses(iNAT_loc,iModeFine) = masses(iNAT_loc,iModeFine) + csNAT
        !natnum(iModeFine) = natnum(iModeFine) + masses(iNAT_loc,iModeFine)*rules%mwNAT / &  !RISTO: WRONG, should be csNAT instead of masses, or omit the natnum on the r.h.s
        !                                      & (density*rules%bin_particle_volume(iModeFine))
      else                          ! Evaporation
        do iBin = 1, nBins
          if (masses(iNAT_loc,iBin) >= - csNAT .and. masses(iNAT_loc,iBin) > 0.) then  ! partial evaporation
            !natnum(iBin) = natnum(iBin) * (-csNAT) / masses(iNAT_loc,iBin) !RISTO: WRONG
            masses(iNAT_loc,iBin) = masses(iNAT_loc,iBin) + csNAT
            hno3g = hno3g - csNAT
            h2og  = h2og - 3.*csNAT
            exit
          else
            csNAT = csNAT + masses(iNAT_loc,iBin)
            hno3g = hno3g + masses(iNAT_loc,iBin)
            h2og  = h2og + 3.*masses(iNAT_loc,iBin)
            masses(iNAT_loc,iBin) = 0.
            !natnum(iBin) = 0.
          endif
        end do
      endif
      
    else
      !
      ! This is new NAT, nothing is in masses and vol_tot. Put everything into the first bin
      ! To keep exact budget push the mass_tot (may be small but still >0) to gas phase
      !
      hno3g = hno3g + mass_tot
      h2og = h2og + 3.*mass_tot
      csNAT = amax1(hno3g - hno3eq, 0.)              ! condensed phase concentration, mole/m3
      if (3.*csNAT > h2og) then
        csNAT = h2og/3.  ! check for water availability
        h2og = 0.
      else
        h2og = h2og - 3.*csNAT
      endif
      hno3g = hno3g - csNAT                     ! gas phase concentration
      masses(iNAT_loc,1:nBins) = 0.
      masses(iNAT_loc,iModeFine) = csNAT        ! HNO3 in particles
      !natnum(iModeFine) = masses(iNAT_loc,iModeFine)*rules%mwNAT/(density*rules%bin_particle_volume(iModeFine)) ! 0.1um (diam) particle volume: 4.188790*e-21 m3
    endif
if(.not. hno3g >= 0.)then
  call msg('Negative hno3 after initial partitioning (hno3g and hno3eq):',hno3g,hno3eq)
  call msg('rules%mwHNO3, rules%mwH2O, mass_tot, csNAT, masses(iNAT_loc,1:nBins), natnum(1:nBins)', &
         & (/rules%mwHNO3, rules%mwH2O, mass_tot, csNAT, masses(iNAT_loc,1:nBins), natnum(1:nBins)/))
  call set_error('Non-positive hno3g','NAT')
endif
      
!    hno3c = csNAT
    !
    ! Describe gradual formation of very coarse particles. We know that within 84 hrs,
    ! about 0.0002 - 0.002 #/cm3 = 200-2000 #/m3 of such particles are made of small ones.
    ! The smaller concentration corresponds to the particle size of ~15um, the larger one - to ~3.5um
    ! The trick is that the growth is not homogeneous: many small particles probably disappear
    ! giving out their mass to the growing few. Hense, we shall use the funny description
    ! of growth: the rate of formation of large ready-made particles.
    ! It leads to the increase rate of ~2.4 - 24 #/m3hour ~ 6.7e-4 - 6.7e-3 #/m3sec
    ! These particles are sitting in the coarse bin. Volume of single 24um particle = 5.8e-14 m3
    ! Do not forget: small particles have to provide enough material for the larger ones. Fade-out if not
    !
    csNAT = rules%rock_formation_rate * seconds * &
          & rules%bin_particle_volume(iModeCoarse) * density / rules%mwNAT
    mass_tot = sum(masses(iNAT_loc,1:nBins)) - masses(iNAT_loc,iModeCoarse)  ! all what is available
    csNAT = csNAT * mass_tot / (mass_tot + 100. * csNAT) !Scale. This is the additional amount we try to put into Coarse mode.

    if(mass_tot >= csNAT)then
      !
      ! Small bins have enough mass, part of which is sent to rocks. Eat mid sizes first
      !
      !RH2017note: Below the natnum is calculate incorrectly because csNAT is scaled above.
      !            Also for other bin(s) the natnum should be recalculated! NOW done at the end of subroutine.
      !natnum(iModeCoarse) = natnum(iModeCoarse) + rules%rock_formation_rate * seconds  ! #/m3 of rocks
      
      do iBin = nBins, 1, -1  ! collect the mass from all bins
        if(iBin == iModeCoarse)cycle
        if(masses(iNAT_loc,iBin) > csNAT)then
          ! enough mass in this bin
          masses(iNAT_loc,iModeCoarse) =  masses(iNAT_loc,iModeCoarse) + csNAT
          masses(iNAT_loc,iBin) = masses(iNAT_loc,iBin) - csNAT
          exit    ! all is transferred
        else
          ! not enough mass in this bin
          masses(iNAT_loc,iModeCoarse) =  masses(iNAT_loc,iModeCoarse) + masses(iNAT_loc,iBin)  ! transferred HNO3
          csNAT = csNAT - masses(iNAT_loc,iBin)
          masses(iNAT_loc,iBin) = 0.
        endif
      end do
    else
      !
      ! Not enough mass in smaller bins. Put all to rocks!
      !
      masses(iNAT_loc,iModeCoarse) = masses(iNAT_loc,iModeCoarse) + mass_tot
      !natnum(1:nBins) = 0.
      !natnum(iModeCoarse) = masses(iNAT_loc,iModeCoarse) * rules%mwNAT / &
      !                    & (density * rules%bin_particle_volume(iModeCoarse))
      do iBin = 1, nBins
        if(iBin /= iModeCoarse) masses(iNAT_loc,iBin) = 0.
      end do
    endif  ! enough material for NAT rocks?

    !Re-calcuate the natnum (RH2017update)
    do iBin = 1, nBins
       natnum(iBin) = masses(iNAT_loc,iBin)*rules%mwNAT/(density*rules%bin_particle_volume(iBin))
    end do

!    !
!    ! Having determined the number concentration in each bin, get the sad
!    !
!    sad = 0.
!    do iBin = 1, nBins
!      if(natnum(iBin) < 1.e-5)cycle  ! skip if empty bin
!      sad = sad + (3. * (masses(iH2O_loc,iBin) + masses(iHNO3_loc,nBins)) / natnum(iBin)) ** (2. / 3.)
!    end do
!    sad = sad * 1.e12 * Pi_4_power_1_3    ! remaining scaling and to um2/m3
    
    return
  end subroutine NAT


  !**********************************************************************************

  subroutine ICE(icenum, &  ! (r) = # density of ICE particles (unchanged), m-3 
               & h2oeq, &   ! (r) = equil. H2O concentration above ICE (unchanged), [mole/m3]
               & hno3eq, &  ! (r) = equil. HNO3 concentration above ICE (unchanged),[mole/m3] 
               & h2og, &    ! (r) = actual H2O concentration in gas phase (changed), [mole/m3]
               & hno3g, &   ! (r) = actual HNO3 concentration in gas phase (changed),[mole/m3]
               & h2so4, &  !(r)= actual H2SO4 gas concentration (changed), [mole/m3]
               & mass_SO4, &  ! (r) = concentration of fine sulphates (changed), [mole/m3]
               & masses, &    ! (r) = concentration all species in particles (changed), [mole/m3]
!               & sad, &     ! (r) = sad of ICE (changed), um2 m-3
               & rules, &     ! Aerosol dyanmics rules
               & rulesRates, &
               & density, &
               & print_it)
!c
!c     This subroutine calculates
!c       - the partition of H2O and HNO3 between gas phase and ICE
!c       - the radius of the ICE particles
!c       - the sad of ICE
!c     
!c     The results are stored in H2O_cnc_molec_air, ck(5), cs(2), cs(3), 
!c     rp2 and sad.     
!c
!c     PSC type 2, water ice (ICE).
!c     The equilibrium pressure of H2O and HNO3 above ICE is used,
!c     see Marti and Mauersberger, GRL, 20 (1993) 363., and
!c     Hanson and Mauersberger, GRL, 15 (1988) 855.
!c
!      ATTENTION. FinROSE code does not have ICE particle components transported, i.e. all 
!                 is in gas phase. In our case, however, aerosols are transported, i.e. full
!                 budget needs to be maintained.
   !
    ! Our treatment will include:
    ! - nucleation (heterogeneous and later, may be, homogeneous) of ice
    ! - condensation of water on ice cryctals, evaporation from these
    ! - condensation of HNO3 on ice crystals
    !
!c
!c     Author
!c     ------
!c
!c     FMI developers, 25.11.2002, FMI.
!c
    implicit none  

    ! Imported parameters
    real, intent(in) :: h2oeq,hno3eq
    real, intent(inout) :: h2og,hno3g,h2so4,mass_SO4
    real, dimension(:), intent(inout) :: icenum         ! condensed species in several bins [#/m3]
    real, dimension(:,:), intent(inout) :: masses       ! concentration of condensed species in several bins [mole/m3]
    type(Tchem_rules_AerDynMidAtm), intent(in) :: rules
    type(TLocalRates), intent(inout) :: rulesRates
    real, intent(out) :: density
    logical, intent(in) :: print_it

    ! Local variables
    integer :: iBin
    !
    ! Number concentrations of the ice particles
    ! Get the ice mass, simply requiring that all super-saturation goes into ice
    !
    density = fu_dry_part_density(fu_get_material_ptr('H2O'))
    h2og = h2og + sum(masses(iH2O_loc,1:nBins))
!    masses(iH2O_loc,1:nBins) = 0.
    if(h2og > 1.0001 * h2oeq)then  ! water in air can be very close to the equilibrium state
      masses(iH2O_loc,iModeCoarse) = h2og - h2oeq          ! h2o concentration in condensed phase
      h2og = h2oeq       !g - masses(iH2O_loc,iModeCoarse)  new remaining h2o in gas phase
    else
      masses(iH2O_loc,iModeCoarse) = 0.
    endif
    if(h2og < 0. .or. masses(iH2O_loc,iModeCoarse) < 0.)then
      call msg('h2og or h2oc < 0 in ICE, h2og,h2oc(.),h2oeq:', (/h2og,masses(iH2O_loc,iModeCoarse),h2oeq/))
      call set_error('h2og or h2oc < 0', 'ICE')
      h2og = max(0.,h2og + masses(iH2O_loc,iModeCoarse) )  ! avoid crashes until exit
      masses(iH2O_loc,iModeCoarse) = 0.
    endif
    !
    ! HNO3 in ICE & around. Note that ICE is only coarse bin
    !
    hno3g = hno3g + masses(iHNO3_loc,iModeCoarse)    ! Get the total HNO3 into gas phase
    if(hno3g > hno3eq)then
      masses(iHNO3_loc,iModeCoarse) = hno3g - hno3eq  ! hno3 concentration in condensed phase
      hno3g = hno3eq                                  ! hno3 gas phase reaches equilibrium
    else
      masses(iHNO3_loc,iModeCoarse) = 0.   ! HNO3 pressure < equilibrium => nothing in condensed phase
    endif  ! HNO3 supersaturated over ice
    !
    ! Sulphates and H2SO4 are simply sent to ice fraction
    !
    masses(iH2SO4_loc,iModeCoarse) = masses(iH2SO4_loc,iModeCoarse) + h2so4 + mass_SO4
    h2so4 = 0.
    mass_SO4 = 0.
    !
    ! Two ways to get number density: prescribe it (FinROSE approach) or
    ! calculate it assuming the prescribed particle sizes
    !
    if(rules%iSpectrumDefinitionSwitch == prescribe_number_concentr)then
      !
      ! For reactions we need surface area density, sum over all modes
      ! FinROSE solution: if PSCs exist, ICE will have number concentration 0.04e6 # m-3
      ! 
      !icenum(1:nBins) = 0.
      !icenum(iModeCoarse) = 0.04e6   ! [m-3]
      !R.H. 2017update: Use the values set at the begin of the module and avoid repetitions. 
      icenum(iModeFine) = rules%fixed_number_concentration(iModeFine,binICE)
      icenum(iModeCoarse) = rules%fixed_number_concentration(iModeCoarse,binICE)
    elseif(rules%iSpectrumDefinitionSwitch == prescribe_sizes)then
      !
      ! However, different paradigm can be used: prescribe size of the particles assuming
      ! that ice is about 10-15um, same as NAT bricks. So, we use the same rock bin.
      !
      do iBin = 1, nBins
        icenum(iBin) = masses(iH2O_loc,iBin) * rules%mwH2O / (density * rules%bin_particle_volume(nBins))
      end do
    else
      call set_error('Dyhnamic spectrum repreesntation is not ready yet','ICE')
    endif ! ifPrescribeNumConc
    !!
    !! Having determined the number concentration in each bin, get the sad
    !!
    !sad = 0.
    !do iBin = 1, nBins
    !  if(icenum(iBin) < 1.e-5)cycle  ! skip if empty bin
    !  sad = sad + (3. * masses(iH2O_loc,iBin) / icenum(iBin)) ** (2. / 3.)
    !end do
    !sad = sad * 1.e12 * Pi_4_power_1_3    ! remaining scaling and to um2/m3

  end subroutine ICE


  !**************************************************************************************

!RISTO Cl&Br: Only the concentrations of the input species for reactions (9-13) are needed! 
  subroutine rates_pseudo_first_order(ifICE, ifNAT, ifSTS, &  ! (i) = PSC and ICE flags (unchanged)
                                    & dt, &  !    (r) = chemistry timestep (unchanged), s
                                    & T, &  !     (r) = temperature (unchanged), K 
                                    !& ppa, &  !     (r) = pressure (unchanged), pa
                                    !& cncH2O, cncHCl, cncHOCl, cncClONO2, &  !(r) = gases, mole/m3
                                    !3Aug2017 by R.H. added the missing species
                                    & cncH2O,  cncHCl, cncHOCl, cncClONO2, &  !(r) = gases, mole/m3
                                    & cncN2O5, cncHBr, cncHOBr, cncBrONO2, &  !(r) = gases, mole/m3
                                    & cncHNO3, &
                                    & sadSTS, &  !  (r) = sad of binary and ternary aer. (unchanged), m2/m3
                                    & sadNAT, &  !  (r) = sad of nat (unchanged), m2/m3
                                    & sadICE, &  !  (r) = sad of ice (unchanged), m2/m3
                                    & sadOTH, &  !(r) = sad of non-PSC species (unchanged), m2 cm-3
                                    & rp, &    !  (r) = radius of STS particles, [m]
                                    & wt_h2so4_mass_prc, &    !  (r) = mass wt of h2so4 (unchanged), %
                                    & rho1, &  !  (r) = density of STS aerosol (unchanged), kg m-3
                                    & rules, &  !  (r) = MAD rules. Hold all coefs, (changed), s-1
                                    & rulesRates, &
                                    & print_it)
!c
!c     Subroutine for calculation of reaction rate coefficients (pseudo
!c     first order) for heterogeneous reactions on/in binary and ternary 
!c     aerosols and on nat and ice particles.
!c
!c     The results are stored in rulesRates%h_whatever (h(1)-h(13))
!c
!c     The uptake coefficients are based on recommendations in JPL1997&
!c     2000 (DeMore et al., 1997&2000) and IUPAC (Atkinson et al., 2001)
!c     unless otherwise stated. The effect of hno3 in the ternary
!c     solutions has not been considered. A decrease in the uptake
!c     coefficient is expected with increasing hno3 concentration.
!c     The reaction rates for the heterogeneous reactions are
!c     calculated using the maximum kinetic mass flux (kmax) and
!c     uptake coefficients (ratio between actual mass flux and kmax).
!c
!c       actual mass flux = uptake coefficient * kmax
!c
!c     The uptake coefficient includes all processes controlling
!c     the mass transport (gas and liquid phase diffusion, mass 
!c     accommodation, chemical reactions).
!c     See Seinfeld and Pandis, 2006, p.91.
!c     Collision frequency for 3d random motion Zi, molecules cm-2 s-1
!c         Zi  = 0.25*ck(i)*(mean speed of molecule i)
!c         100. cm/m
!c         Zi  = 0.25*ck(i)*100.*sqrt(8*T*k/(pi*mi)) 
!c         k   = R/Na
!c         mi  = mwi/Na
!c         Zi  = 25.*ck(i)*sqrt(8*T*R/(pi*mwi)) = cf*ck(i)
!c         cf  = 25.*sqrt(8*T*R/(pi*mw))
!c
!c     Reaction rate, molecules cm-3 s-1
!c         rr = gam * (Zi * sad) = gam*cf*sad*ck(i) 
!c
!c     and reaction rate coefficient, s-1
!c         h = gam*cf*sad
!c
!c     The numbering of the heterogeneous reactions:
!c     ( 1)  clono2 +  h2o -> hocl  +  hno3
!c     ( 2)  brono2 +  h2o -> hobr  +  hno3                     
!c     ( 3)  n2o5   +  h2o ->        2 hno3                
!c     ( 4)  clono2 +  hcl -> cl2   +  hno3          
!c     ( 5)  hocl   +  hcl -> cl2   +  h2o          
!c     ( 6)  brono2 +  hcl -> brcl  +  hno3
!c     ( 7)  hobr   +  hcl -> brcl  +  h2o 
!c     ( 8)  n2o5   +  hcl -> clno2 +  hno3
!c     ( 9)  clono2 +  hbr -> brcl  +  hno3
!c     (10)  hocl   +  hbr -> brcl  +  h2o
!c     (11)  brono2 +  hbr -> br2   +  hno3
!c     (12)  hobr   +  hbr -> br2   +  h2o 
!c     (13)  n2o5   +  hbr -> brno2 +  hno3
!c
!c
!    ATTENTION !! UNITS ARE NOT SI
!    ATTENTION !! UNITS ARE NOT SI
!
!c     Author
!c     ------
!c
!c     FMI developers, 25.11.2002, FMI.
!c
!c-----------------------------------------------------------------------
   
    implicit none

    ! Imported parameters
    logical, intent(in) :: ifice, ifNAT, ifSTS
    !real, intent(in) :: dt,T,ppa, cncH2O, cncHCl, cncHOCl, cncClONO2, &
    !3Aug2017 by R.H. added the additional gasses for Cl & Br reactions
    real, intent(in) :: dt,T,cncH2O, cncHCl, cncHOCl, cncClONO2, cncHBr, cncHOBr, cncBrONO2, cncN2O5, cncHNO3, &         
                      & sadSTS,sadNAT,sadICE,sadOTH, rp,wt_h2so4_mass_prc,rho1 !,ppa
    type(Tchem_rules_AerDynMidAtm), intent(in) :: rules
    type(TLocalRates), intent(inout) :: rulesRates
    logical, intent(in) :: print_it

    ! Local variables
    real :: cf_clono2, cf_brono2, cf_n2o5, cf_hocl, cf_hobr, wt_mole_h2so4
    real :: pclono2,phcl,ph2o0,ph2o,pHBr,Mh2so4,A,t0,T_over_nya,ah,x,aw,b0,b1,b2
    real :: Sclono2,Hclono2,Dclono2,Hhcl,Mhcl,kh2o,kh,khydr,khcl,Gbh2o,RgasAtm,Hhbr
    real :: lclono2_over_rp, lhocl_over_rp
    real :: fclono2,Grxn,Gbhcl,Gs,Fhcl,Gsp,Gbphcl,Gb,gam,thetaHBr,thetaHCl,KlangHCl,KlangHNO3
    real(r8k) :: RgasScaled,Shocl,Hhocl,Dhocl,k5,lhocl,fhocl,Ghocl
    logical :: ifH2O, ifHCl

    !
    ! Mean velocity of the molecules, [m/s]. Note molar mass and R instead
    ! of molecular mass amd k. Ratio is the same.
    ! A quarter comes from 1/4*v_gas*cnc_gas - the collision frequency
    !
    cf_clono2 = 0.25*sqrt(8*T*gas_constant_uni/(pi*rules%mwclono2))
    cf_brono2 = 0.25*sqrt(8*T*gas_constant_uni/(pi*rules%mwbrono2)) !Enabled 2017 by R.H. to include remaining equations.
    cf_n2o5   = 0.25*sqrt(8*T*gas_constant_uni/(pi*rules%mwn2o5))
    cf_hocl   = 0.25*sqrt(8*T*gas_constant_uni/(pi*rules%mwhocl))
    cf_hobr   = 0.25*sqrt(8*T*gas_constant_uni/(pi*rules%mwhobr)) !Enabled 2017 by R.H. to include remaining equations.

    if(ifice)then
      !
      ! reaction rate coefficients, pseudo first order for ice
      !
      ! The rate is the product of:
      ! - gamICE probabilities of reaction to happen in case of collision
      ! - cf* collision frequency without gas concentration: v_gas/4, [m/s]
      ! - sad surface area density, [m2/m3]
      !
#ifdef DEBUG_MORE
if(print_it)then
call msg('Rates. ICE')
endif
#endif
      !
      ! Reaction rate gamma values for ICE particles. The old (ver5.5) values are stored in gamICE(?) = old_value. 
      !
      !IUPAC for ClONO2+H2O on ice: 1/gamma = 1/alpha_s+c/(4*k_s*K_LinC*[H2O]_s)
      !      with alpha_s=0.5,    k_s*K_LinC=5.2e-17*exp(2032/T) (in cm3/molecules*s)   [H2O]_s = 1e15*(1-0.81*theta_HNO3)  (in molecule/cm2)
      !      Assuming steady state-> theta_HNO3=1. Note also that c/4 corresponds to our cf_clono2 and needs to multiply by 100 to get in cm/s.
      rulesRates%gamICE_clono2_h2o = 1.0/(2.0+100*cf_clono2/(5.2e-2*exp(2032./T)*(1.0-0.81))) !gamICE(1) !IUPAC
                                            !    gamICE(1) = 0.30       !JPL2015: gamma=0.3;
      rulesRates%gamICE_brono2_h2o = 5.3e-4*exp(1100.0/T) !gamICE(2)    !IUPAC
                                            !    gamICE(2) = 0.30       !IUPAC: gamma=5.3e-4*exp(1100.0/T); JPL2015: gamma>0.2;
      rulesRates%gamICE_n2o5_h2o   = 0.02   !    gamICE(3) = 0.02       !JPL2015: gamma=0.02; IUPAC: gamma=0.02
      !rulesRates%gamICE_clono2_hcl = 0.24  !    gamICE(4) = 0.30       !IUPAC: gamma_max=0.24; JPL2015: gamma=0.3;
      !IUPAC: gamma_clono2_hcl = 0.24*theta_HCl where theta_HCl= K_LangC(HCl)[HCl]/(1+K_LangC(HCl)[HCl]+K_LangC(HNO3)[HNO3])
      !       Here [HCl] and [HNO3] are concentrations in molecule/cm3 for HCl and HNO3 gas, respectively.
      !       NOTE that IUPAC wrongly uses the same K_LangC also for HNO3. Temperature dependence of gamma_ER=0.24 is omitted due to scatter. 
      !       The Langmuir equilibrium constants K_LangC=K_LinC/N_max are
      !       for HCl : N_max = 3.0e14 molecule/cm2 and K_LinC = 1.3e-5*exp(4600/T) [These IUPAC values produce the Fig.1 but differ from newest HCl+ice_V.A1.27.pdf]
      !       for HNO3: N_max = 2.7e14 molecule/cm2 and K_LinC = 7.5e-5*exp(4585/T) [See IUPAC Sheet HNO3+ice_V.A1.12.pdf for HNO3+ice.]
      KlangHCl  = 1.3e-5*exp(4600.0/T)/3.0e14;
      KlangHNO3 = 7.5e-5*exp(4585.0/T)/2.7e14;
      rulesRates%gamICE_clono2_hcl = 0.24*KlangHCl*(cncHCl*1e-6*avogadro)/(1+KlangHCl*(cncHCl*1e-6*avogadro)+KlangHNO3*(cncHNO3*1e-6*avogadro)) !gamICE(4)   !IUPAC
      !rulesRates%gamICE_hocl_hcl   = 0.20  !    gamICE(5) = 0.20       !JPL2015: gamma=0.2; 
      !IUPAC: Equation for HOCl+HCl should be read with care and see the comments. The equation misses 1/N_max.
      !       We also use gamma_ER=0.15 to better reproduce the figure. One might want to add the HNO3 surface coverage limitation as for ClONO2+HCl.
      KlangHCl=19200/2e14  !This is taken from the new 2016 June evaluation for HCl+ice_V.A1.27.pdf
      rulesRates%gamICE_hocl_hcl = 0.15*KlangHCL*(cncHCl*1e-6*avogadro)/(1+(cncHCl*1e-6*avogadro));   ! gammaICE(5)  !IUPAC
      rulesRates%gamICE_brono2_hcl = 0.30   !    gamICE(6) = 0.30       !IUPAC: gamma=0.30;  JPL2015: gamma 0.15, 0.34, 0.26, or even very close to 1
      rulesRates%gamICE_hobr_hcl   = 0.25   !    gamICE(7) = 0.30       !IUPAC: gamma=0.25; JPL2015: gamma=0.3
      rulesRates%gamICE_n2o5_hcl   = 0.03   !    gamICE(8) = 0.03       !JPL2015: gamma=0.03; IUPAC: no recommendation, refers e.g. to Leu 1988, where gamma=0.028, 0.037, 0.064.
      !IUPAC for ClONO2+HBr on ice: gamma = 0.56*theta_HBr and surface coverage theta_HBr = 4.14e-10*[HBr]**0.88 (valid at 188 K) where [HBr] is concentration in molecules/cm3.
      thetaHBr = min(1.0,4.14e-10*(cncHBr*1e-6*avogadro)**0.88)
      rulesRates%gamICE_clono2_hbr = 0.56*thetaHBr ! gamICE(9) = 0.30   !IUPAC: NOTE that thetaHBr is taken at 188K; JPL2015: gamma>0.3
      !For gamICE_hocl_hbr IUPAC gives an equation that is valid at 188 K, but the form of equation is unclear (buggy or wrong units).
      rulesRates%gamICE_hocl_hbr   = 0.05   !    gamICE(10)= 0.05       !IUPAC; JPL2015: gamma = 0.01..0.07 (220 K) and to 0.06..0.38 (189 K) depending on HBr concentration. 
      !2Aug2017 by R.H. enabled the following reactions coefficients:
      rulesRates%gamICE_brono2_hbr = 6.6e-4*exp(700.0/T)    !gamICE(11)= 0.0(n.a.)  !IUPAC: gamma=6.6e-4*exp(700/T); JPL2015: gamma_0 > 0.1 
      rulesRates%gamICE_hobr_hbr   = 4.8e-4*exp(1240.0/T)   !gamICE(12)= 0.10       !IUPAC: gamma=4.8e-4*exp(1240/T), provided that [HBr]>[HOBr]; JPL2015: gamma>0.1
      rulesRates%gamICE_n2o5_hbr   = 0.08   !    gamICE(13)= 0.10       !IUPAC: no recommendation (gamma=0.02...0.15); JPL2015 (Seisel et al.) gamma from 3e-3 to 0.1
    
      !
      ! Non-PSC aerosols: Al2O3, soot, etc. Have very few coefs
      ! Here, we do not have gamSTS for sulphates, have to use values of ice clouds
      ! for foreign aerosol if other value does not exist
      !
      rulesRates%gamOT_clono2_h2o = rulesRates%gamICE_clono2_h2o    ! gamOT(1)
      rulesRates%gamOT_brono2_h2o = rulesRates%gamICE_brono2_h2o    ! gamOT(2)
      rulesRates%gamOT_n2o5_h2o = rulesRates%gamICE_n2o5_h2o        ! gamOT(3)
      !rate for reaction (4): clono2_hcl is set at the begin of the module in type "TLocalRates"
      rulesRates%gamOT_hocl_hcl = rulesRates%gamICE_hocl_hcl        ! gamOT(5)
      rulesRates%gamOT_brono2_hcl = rulesRates%gamICE_brono2_hcl    ! gamOT(6)
      rulesRates%gamOT_hobr_hcl = rulesRates%gamICE_hobr_hcl        ! gamOT(7)
      rulesRates%gamOT_n2o5_hcl = rulesRates%gamICE_n2o5_hcl        ! gamOT(8)
      rulesRates%gamOT_clono2_hbr = rulesRates%gamICE_clono2_hbr    ! gamOT(9)
      rulesRates%gamOT_hocl_hbr = rulesRates%gamICE_hocl_hbr        ! gamOT(10)
      !2Aug2017 by R.H. enabled the following 2 coefficients for reactions (11-12)
      rulesRates%gamOT_brono2_hbr = rulesRates%gamICE_brono2_hbr    ! gamOT(11)
      rulesRates%gamOT_hobr_hbr = rulesRates%gamICE_hobr_hbr        ! gamOT(12)
      rulesRates%gamOT_n2o5_hbr = rulesRates%gamICE_n2o5_hbr        ! gamOT(13) = 0.005 <= This JPL2000 value does not exist in their tables!!??

      rulesRates%h_clono2_h2o =(rulesRates%gamICE_clono2_h2o * sadICE + rulesRates%gamOT_clono2_h2o * sadOTH)* cf_clono2
      !2Aug2017 by R.H. enabled the following coefficient for reaction (2)
      rulesRates%h_brono2_h2o =(rulesRates%gamICE_brono2_h2o * sadICE + rulesRates%gamOT_brono2_h2o * sadOTH)* cf_brono2
      rulesRates%h_n2o5_h2o =  (rulesRates%gamICE_n2o5_h2o * sadICE + rulesRates%gamOT_n2o5_h2o * sadOTH)* cf_n2o5
      rulesRates%h_clono2_hcl =(rulesRates%gamICE_clono2_hcl * sadICE + rulesRates%gamOT_clono2_hcl * sadOTH)* cf_clono2
      rulesRates%h_hocl_hcl =  (rulesRates%gamICE_hocl_hcl * sadICE + rulesRates%gamOT_hocl_hcl * sadOTH)* cf_hocl
      !2Aug2017 by R.H. enabled the following 2 coefficients for reactions (6-7)
      rulesRates%h_brono2_hcl =(rulesRates%gamICE_brono2_hcl * sadICE + rulesRates%gamOT_brono2_hcl * sadOTH)* cf_brono2
      rulesRates%h_hobr_hcl =  (rulesRates%gamICE_hobr_hcl * sadICE + rulesRates%gamOT_hobr_hcl * sadOTH)* cf_hobr
      rulesRates%h_n2o5_hcl =  (rulesRates%gamICE_n2o5_hcl * sadICE + rulesRates%gamOT_n2o5_hcl * sadOTH)* cf_n2o5
      !2Aug2017 by R.H. enabled the following 5 coefficients for reactions (9-13)
      rulesRates%h_clono2_hbr =(rulesRates%gamICE_clono2_hbr * sadICE + rulesRates%gamOT_clono2_hbr * sadOTH)* cf_clono2
      rulesRates%h_hocl_hbr =  (rulesRates%gamICE_hocl_hbr * sadICE + rulesRates%gamOT_hocl_hbr * sadOTH)* cf_hocl
      rulesRates%h_brono2_hbr =(rulesRates%gamICE_brono2_hbr * sadICE + rulesRates%gamOT_brono2_hbr * sadOTH)* cf_brono2
      rulesRates%h_hobr_hbr =  (rulesRates%gamICE_hobr_hbr * sadICE + rulesRates%gamOT_hobr_hbr * sadOTH)* cf_hobr
      rulesRates%h_n2o5_hbr =  (rulesRates%gamICE_n2o5_hbr * sadICE + rulesRates%gamOT_n2o5_hbr * sadOTH)* cf_n2o5
      if (.not. all((/rulesRates%h_clono2_h2o, rulesRates%h_n2o5_h2o, rulesRates%h_clono2_hcl, rulesRates%h_hocl_hcl, rulesRates%h_n2o5_hcl, & !Cl-reactions (1,3-5,8)
                    & rulesRates%h_brono2_h2o, rulesRates%h_brono2_hcl, rulesRates%h_hobr_hcl, rulesRates%h_clono2_hbr, rulesRates%h_hocl_hbr, & !Br-reactions (2,6,7,9,10)
                    & rulesRates%h_brono2_hbr, rulesRates%h_hobr_hbr, rulesRates%h_n2o5_hbr /) >= 0.)) then !new Br-reactions (11-13)
  call set_error('Trouble with heterogeneous rates in ICE','rates_pseudo_first_order')
  call msg('Rates with ICE: h_clono2_h2o  h_n2o5_h2o  h_clono2_hcl  h_hocl_hcl  h_n2o5_hcl')
  call msg('Rates with ICE:',(/rulesRates%h_clono2_h2o, rulesRates%h_n2o5_h2o, rulesRates%h_clono2_hcl, rulesRates%h_hocl_hcl, rulesRates%h_n2o5_hcl/))
  call msg('Rates with ICE: h_brono2_h2o h_brono2_hcl  h_hobr_hcl  h_clono2_hbr h_hocl_hbr')
  call msg('Rates with ICE:',(/rulesRates%h_brono2_h2o, rulesRates%h_brono2_hcl, rulesRates%h_hobr_hcl, rulesRates%h_clono2_hbr, rulesRates%h_hocl_hbr/))
  call msg('Rates with ICE: h_brono2_hbr h_hobr_hbr    h_n2o5_hbr')
  call msg('Rates with ICE:',(/rulesRates%h_brono2_hbr, rulesRates%h_hobr_hbr, rulesRates%h_n2o5_hbr/))
  
  call msg("rulesRates%gamICE_clono2_h2o, sadICE, rulesRates%gamOT_clono2_h2o, sadOTH, cf_clono2")
  call msg("",(/rulesRates%gamICE_clono2_h2o, sadICE,  rulesRates%gamOT_clono2_h2o , sadOTH , cf_clono2/))
      end if

    endif  ! ifICE
    
    if(ifSTS .or. ifNAT)then
      !
      ! aerosol properties in non-ice conditions
      !
#ifdef DEBUG_MORE
if(print_it)then
call msg('Rates. STS/NAT')
endif
#endif

      RgasAtm = gas_constant_uni/std_pressure_sl ! Gas constant unit R in m3*atm/(K*mol)
      !r2 = 0.082057        !Gas constant R in atm/(K*M) (M=mol/liter)
      !Below some equations require the gas constant R in units of atm*L/(K*mol)=atm*dm3/(K*mol) and cf_clono2 and cf_hocl in cm/s.
      !Because above the cf_* are in m/s we scale the gas unit additionlly to get things right. Also, note that cf_* already contains a factor 0.25.  
      RgasScaled = 10*RgasAtm    !0.01 times the gas constant in units of L*atm/(K*mol)

      Mh2so4 = 1.e-3*rho1*wt_h2so4_mass_prc/9.8 ! molarity, mol dm-3
      A =  169.5+5.18*wt_h2so4_mass_prc-0.0825*wt_h2so4_mass_prc**2+3.27e-3*wt_h2so4_mass_prc**3.
      t0 = 144.11+0.166*wt_h2so4_mass_prc-0.015*wt_h2so4_mass_prc**2+2.18e-4*wt_h2so4_mass_prc**3.
      ! nya -- viscosity
      T_over_nya = T**(2.43)*exp(-448./(max(T-t0,1.))) / A   ! T / nya 

      ah = exp(60.51-0.095*wt_h2so4_mass_prc+0.0077*wt_h2so4_mass_prc**2-1.61e-5*wt_h2so4_mass_prc**3. &  !   acid activity (in molarity)
         & -(1.76+2.52e-4*wt_h2so4_mass_prc**2)*sqrt(T)+(-805.89+253.05*wt_h2so4_mass_prc**0.076)/sqrt(T))

      wt_mole_h2so4 = wt_h2so4_mass_prc/(wt_h2so4_mass_prc+(100.-wt_h2so4_mass_prc)*98./18.) ! mole fraction

      ! saturation water vapour, Pa
      ph2o0 = 100.0*exp(18.452406985-3505.1578807/T-330918.55082/(T*T)+ 12725068.262/T**3.) 
      !ph2o = ppa * cncH2O * rulesRates%cnc2vmr    ! qn(3) in Pa
      ph2o = cncH2O*gas_constant_uni*T             ! qn(3) in Pa
      ifH2O = ph2o > 1e-10
      
      ! water activity
      !aw = exp((-69.775*wt_mole_h2so4-18253.7*wt_mole_h2so4*wt_mole_h2so4+31072.2*wt_mole_h2so4**3.- &
      !         & 25668.8*wt_mole_h2so4**4.)*(1/T-26.9033/(T*T)))
      aw = ph2o/ph2o0 !Use this since the water partial pressure is known.

      !phcl = ppa / 101325 * cncHCl * rulesRates%cnc2vmr     ! qn(13)     !atm
      phcl = cncHCl*RgasAtm*T                                ! qn(13)     !atm
      Hhcl = (0.094-0.61*wt_mole_h2so4+1.2*wt_mole_h2so4**2)*exp(-8.68+(8515.-10718.*wt_mole_h2so4**0.7)/T)
      Mhcl = Hhcl*phcl
      ifHCl = phcl > 1e-15 !5Oct2017 by R.H.: Changed 1e-10 -> 1e-15
      !
      ! binary and ternary aerosols with water and HCl reactions
      !
      !pclono2 = (ppa/101325.)*cncClONO2*rulesRates%cnc2vmr   ! *qn(14)   !atm
      pclono2 = cncClONO2*RgasAtm*T                      ! *qn(14)   !atm
      Sclono2 = 0.306+24.0/T                             ! Setchenow coefficient, M-1
      Hclono2 = 1.6e-6*exp(4710./T)*exp(-Sclono2*Mh2so4) ! Henry's law constant
      Dclono2 = 5.e-8*T_over_nya                              ! liquid phase diffusion constant
      kh = 1.22e12*exp(-6200./T)

      if(ifHCl .and. ifH2O)then
        !
        ! Both water and HCl are present
        !
#ifdef DEBUG_MORE
if(print_it)then
call msg('Rates. HCl and H2O')
endif
#endif
        kh2o = 1.95e10*exp(-2800./T)
        khydr = aw*(kh2o+kh*ah)             ! reaction (1) liquid phase reaction rate, 2nd order
        khcl = 7.9e11*ah*Dclono2*Mhcl       ! reaction (4) liquid phase reaction rate, 2nd order
        Gbh2o = Hclono2*RgasScaled*T*sqrt((Dclono2*khydr))/cf_clono2
        !lclono2 = sqrt(Dclono2/(khydr+khcl))         ! reaction diffusive length, cm (NOTE: should be in meters as rp)
        lclono2_over_rp = 0.01*sqrt(Dclono2/(khydr+khcl))/rp    ! reaction diffusive length, meters (since particle radius rp is in meters)

        fclono2 = fclono(lclono2_over_rp)
        Grxn = fclono2*Gbh2o*sqrt(1+khcl/(khydr + 1e-10))
        Gbhcl = Grxn*khcl/(khcl+khydr)
        Gs = 66.12*Hclono2*Mhcl*exp(-1374./T)
        Fhcl = phcl / (phcl+rules%m_hcl_clono2_ratio*(Gs+Gbhcl)*pclono2)
        Gsp = Fhcl*Gs 
        Gbphcl = Fhcl*Gbhcl
        Gb = Gbphcl+Grxn*khydr/(khcl+khydr)
        gam = (Gsp+Gb) / (Gsp + Gb + 1.0)                                   !=1./(1.+1./(Gsp+Gb))
        rulesRates%gamSTS_clono2_hcl = gam*(Gsp+Gbphcl)/(Gsp+Gb)              !gamSTS(4) !JPL2000 
        !rulesRates%gamSTS_clono2_h2o = gam-rulesRates%gamSTS_clono2_hcl      !gamSTS(1) !JPL2000   !!DANGEROUS!!!
        rulesRates%gamSTS_clono2_h2o = gam*Grxn*khydr/((khcl+khydr)*(Gsp+Gb)) !gamSTS(1)
        
      elseif(ifHCl)then
        !
        ! No water
        !
#ifdef DEBUG_MORE
if(print_it)then
call msg('Rates. HCl')
endif
#endif
        !khcl = 7.9e11*ah*Dclono2*Mhcl      ! reaction (4) liquid phase reaction rate, 2nd order
        !lclono2 = sqrt(Dclono2/khcl)       ! reaction diffusive length, cm (NOTE: should be in meters as rp)
        lclono2_over_rp = 0.01/sqrt(7.9e11*ah*Mhcl)/rp ! reaction diffusive length, m  (Here we avoid division by zero when Dclono2=0, when nya is large/inf)
        fclono2 = fclono(lclono2_over_rp)
        Gs = 66.12*Hclono2*Mhcl*exp(-1374./T)
        Fhcl = phcl / (phcl+rules%m_hcl_clono2_ratio*Gs*pclono2)
        Gsp = Fhcl*Gs 
        rulesRates%gamSTS_clono2_hcl = Gsp/(1.+Gsp) ! gamSTS(4)
        rulesRates%gamSTS_clono2_h2o = 0.           ! gamSTS(1)
#ifdef DEBUG
if(isNaN(rulesRates%gamSTS_clono2_hcl))then
  call msg('cnc2vmr, rulesRates%gamSTS_clono2_hcl',rulesRates%cnc2vmr, rulesRates%gamSTS_clono2_hcl)
  !call msg('ppa,cnc2vmr',ppa,rulesRates%cnc2vmr)
  call msg('T,mh2so4',T,mh2so4)
  call msg('rulesRates%gamSTS_clono2_hcl is ?. Gsp, Gs:',GsP,Gs)
  call msg('pclono2,Sclono2 ',pclono2,Sclono2)
  call msg('Hclono2,Dclono2',Hclono2,Dclono2)
  call msg('Mhcl',Mhcl)
  call msg('wt_h2so4_mass_prc,wt_mole_h2so4',wt_h2so4_mass_prc) !,wt_mole_h2so4)
  call msg('kh, khcl',kh,khcl)
  call msg('lclono2_over_rp,fclono2',lclono2_over_rp,fclono2)
  call msg('FHcl,cncClONO2',fhcl,cncClONO2)
endif
#endif
      elseif(ifH2O)then
        !
        ! No HCl
        !
#ifdef DEBUG_MORE
if(print_it)then
call msg('Rates. H2O')
endif
#endif
        kh2o = 1.95e10*exp(-2800./T)
        khydr = aw*(kh2o+kh*ah)            ! reaction (1) liquid phase reaction rate, 2nd order
        Gbh2o = Hclono2*RgasScaled*T*sqrt((Dclono2*khydr))/cf_clono2
        !lclono2 = sqrt(Dclono2/khydr)     ! reaction diffusive length, cm (NOTE: should be in meters as rp)
        lclono2_over_rp = 0.01*sqrt(Dclono2/khydr)/rp ! reaction diffusive length, meters (since particle radius rp is in meters)
        fclono2 = fclono(lclono2_over_rp)
        Grxn = fclono2*Gbh2o
        rulesRates%gamSTS_clono2_hcl = 0.               ! gamSTS(4)
        !rulesRates%gamSTS_clono2_h2o = 1./(1.+1./Grxn) ! gamSTS(1)
        rulesRates%gamSTS_clono2_h2o = Grxn/(1.+Grxn)   ! gamSTS(1) !avoid division by zero if Grxn=0
      else
#ifdef DEBUG_MORE
if(print_it)then
call msg('Rates. Clean air')
endif
#endif
        ! Nothing is in air
        rulesRates%gamSTS_clono2_hcl = 0.  ! gamSTS(4)
        rulesRates%gamSTS_clono2_h2o = 0.  ! gamSTS(1)
      endif
! option:     rulesRates%gamSTS_clono2_hcl = 10.**(1.84-0.075*wt)   !iupac, gamSTS(1)
!             rulesRates%gamSTS_clono2_hcl = 0.2                    !FinRose, or see above,  gamSTS(4)


      !The following SILAM value is from Hanson JGR 108, 4239 (2003), but somebody has invented an extra decimal in coefficients :O
      rulesRates%gamSTS_brono2_h2o = 1./(1./0.805+ 1./(exp(29.24-0.396*wt_h2so4_mass_prc)+0.114)) !jpl2000, gamSTS(2) !13Sep2017 enabled by R.H.

      !The following gamSTS(3) rate is consistent with JPL 2015 formula:
      b0 = -25.5265-0.133188*wt_h2so4_mass_prc+0.00930846*wt_h2so4_mass_prc**2-9.0194e-5*wt_h2so4_mass_prc**3.
      b1 = 9283.76+115.354*wt_h2so4_mass_prc-5.19258*wt_h2so4_mass_prc**2+0.0483464*wt_h2so4_mass_prc**3.
      b2 = -851801.-22191.2*wt_h2so4_mass_prc+766.916*wt_h2so4_mass_prc**2-6.85427*wt_h2so4_mass_prc**3.
      rulesRates%gamSTS_n2o5_h2o = exp(b0+b1/T+b2/(T*T)) !jpl2000  gamSTS(3)           

      Dhocl = 6.4e-8*T_over_nya !moved here (to calculate Dhocl and Hhocl) to be used for calculation of rulesRates%gamSTS_hocl_hbr
      Shocl = 0.0776+59.18/T
      Hhocl = 1.91e-6*exp(5862.4/T)*exp(-Shocl*Mh2so4) ! Hhocl could be obtained from subroutine TER
      !if(ifHCl .and. cncHOCl > 1e-20)then
      if (ifHCl) then !Use ifHCl since otherwise some things (above) may not be calculated that are needed for gamSTS_hocl_hcl!
        k5 = 1.25e9*ah*Dhocl*Mhcl
        !lhocl = sqrt(Dhocl/k5)           ! reaction diffusive length, cm (NOTE: should be in meters as rp)
        lhocl = 0.01/sqrt(1.25e9*ah*Mhcl) ! reaction diffusive length, meters. Avoiding division by zero if Dhocl is zero (if viscosity nya is infinity/large)
        !fhocl = 1./tanh(rp/lhocl)-lhocl/rp
        lhocl_over_rp  =  lhocl/rp
        fhocl = fclono(lhocl_over_rp)
        !!!     fhocl = max(min(1./tanh(rp/lhocl)-lhocl/rp,1.0),0.0) !avoid numerics to go beyond the actual limits of 0...1

        Ghocl = Hhocl*RgasScaled*T*sqrt(Dhocl*k5)/cf_hocl !19Sep2017 by R.H. Corrected the equation. 
        if(fhocl*Ghocl*Fhcl <= 0.)then
          ! fhocl can become negative if lhocl is too big ?
          rulesRates%gamSTS_hocl_hcl = 0.      !gamSTS(5)
        else
          rulesRates%gamSTS_hocl_hcl = 1./(1.+1./(fhocl*Ghocl*Fhcl)) !jpl2000,  gamSTS(5) NOTE: Bug in JPL equation (F_HOCl in stead of F_HCl)
        endif
      else
        rulesRates%gamSTS_hocl_hcl = 0.      !gamSTS(5)
      endif
      !Uncomment if later check (below) gives NaN for rulesRates%gamSTS_hocl_hcl.
      !if (.not.( rulesRates%gamSTS_hocl_hcl >= 0.)) then
      !  call msg('gamSTS_hocl_hcl, fhocl,lhocl,rp,cncH2O,T:', (/rulesRates%gamSTS_hocl_hcl, real(fhocl),real(lhocl),rp,cncH2O,T/))
      !  call msg('gamSTS_hocl_hcl, Ghocl,Fhcl,Dhocl,k5:', (/rulesRates%gamSTS_hocl_hcl, real(Ghocl),real(Fhcl),real(Dhocl),real(k5)/))
      !  call msg('phcl, m_hcl_clono2_ratio, Gs, Gbhcl, pclono2:', (/phcl, rules%m_hcl_clono2_ratio, Gs, Gbhcl, pclono2/))
      !  write(*,*) ifHCl, ifH2O
      !  call set_error('rulesRates%gamSTS_hocl_hcl','rates_pseudo_first_order')
      !endif

      !NEW RATES:
      !gamSTS_hocl_hbr should now work, produces the same figure than IUPAC for 228 K, pH2O = 50 ppm (in order to get 69% sulfuric acid)
      !NOTE: Henry law constant for HBr could be abtained from TER routine, but we use the formula here to get the recommended gamma value. 
      Hhbr = 10.0**(1000.0*(-1.977e-4*wt_h2so4_mass_prc**2-2.096e-2*wt_h2so4_mass_prc+4.445)/T &
           & -8.979e-5*wt_h2so4_mass_prc**2+2.141e-2*wt_h2so4_mass_prc-6.067)
      pHBr = cncHBr*RgasAtm*T !HBr partial pressure in atm
      !Gb = Hhocl*RgasScaled*T*sqrt(Dhocl*1e-7*Hhbr*pHBr)/cf_hocl !IUPAC states k_{HOCl+HBr}=1e-7 but this must contain a sign mistake in the exponent.
      Gb = Hhocl*RgasScaled*T*sqrt(Dhocl*1e7*Hhbr*pHBr)/cf_hocl  !Here we use k_{HOCl+HBr}=1e7 1/(M*s) which is close to Abbatt & Nowak value
      rulesRates%gamSTS_hocl_hbr = Gb/(1.0+Gb) !gamSTS(10) !use this, instead of 1.0/(1.0+1.0/Gb), in order to avoid division by zero if Gb=0.
      !Uncomment if later check (below) gives NaN for rulesRates%gamSTS_hocl_hbr.
      !if (.not.( rulesRates%gamSTS_hocl_hbr >= 0.)) then
      !  call msg('gamSTS_hocl_hbr, Hhbr,pHBr,Gb,Hhocl,Dhocl:', (/rulesRates%gamSTS_hocl_hcl, Hhbr,pHBr,Gb,real(Hhocl),real(Dhocl)/))
      !  call msg('gamSTS_hocl_hbr, wt_h2so4_mass_prc, T:', (/rulesRates%gamSTS_hocl_hcl, wt_h2so4_mass_prc, T/))
      !  write(*,*) ifHCl, ifH2O
      !  call set_error('rulesRates%gamSTS_hocl_hbr','rates_pseudo_first_order')
      !endif

      
      ! reaction rate coefficients, pseudo first order in no-ice conditions
      !
      ! ! ATTENTION. FinROSE has 5*gamSTS
      ! and commented-out option for gamSTS=min(0.3,gamSTS) to reflect the mass transfer resistance
      ! The rate is the product of:
      ! - gam* probabilities of reaction to happen in case of collision
      ! - cf* collision frequency without gas concentration, [m/sec]
      ! - sad surface area density,  [m2/m3]
      !
      ! Non-PSC aerosols: Al2O3, soot, etc. Have very few coefs
      ! Here, for foreign aerosol we use gamSTS for sulphates if not available other value:
      !
      rulesRates%gamOT_clono2_h2o = rulesRates%gamSTS_clono2_h2o    ! gamOT(1)
      rulesRates%gamOT_brono2_h2o = rulesRates%gamSTS_brono2_h2o    ! gamOT(2)
      rulesRates%gamOT_n2o5_h2o = rulesRates%gamSTS_n2o5_h2o        ! gamOT(3)
      rulesRates%gamOT_clono2_hcl = 0.02 !rate for reaction (4): clono2_hcl is set at the begin of the module in type "TLocalRates"
      rulesRates%gamOT_hocl_hcl = rulesRates%gamSTS_hocl_hcl        ! gamOT(5)
      rulesRates%gamOT_brono2_hcl = rulesRates%gamSTS_brono2_hcl    ! gamOT(6)
      rulesRates%gamOT_hobr_hcl = rulesRates%gamSTS_hobr_hcl        ! gamOT(7)
      rulesRates%gamOT_n2o5_hcl = rulesRates%gamSTS_n2o5_hcl        ! gamOT(8)
      rulesRates%gamOT_clono2_hbr = rulesRates%gamSTS_clono2_hbr    ! gamOT(9)
      rulesRates%gamOT_hocl_hbr = rulesRates%gamSTS_hocl_hbr        ! gamOT(10)
      !2Aug2017 by R.H. enabled the following 2 coefficients for reactions (11-13). 
      rulesRates%gamOT_brono2_hbr = rulesRates%gamSTS_brono2_hbr    ! gamOT(11)
      rulesRates%gamOT_hobr_hbr = rulesRates%gamSTS_hobr_hbr        ! gamOT(12)
      rulesRates%gamOT_n2o5_hbr = rulesRates%gamSTS_n2o5_hbr        ! gamOT(13) = 0.005 <= This JPL2000 value does not exist in their tables!!??

      !!
      !! Reaction rate gamma values for NAT-particles. The old (ver5.5) values are stored in gamNAT(?) = value 
      !!
      rulesRates%gamNAT_clono2_h2o = 0.004      ! gamNAT(1) = 0.004     !JPL2015; Sessler et al. JGR1996: 0.006
      rulesRates%gamNAT_brono2_h2o = 0.006      ! gamNAT(2) = 0.006     !FinROSE, Sessler et al. JGR1996
      rulesRates%gamNAT_n2o5_h2o   = 0.0004     ! gamNAT(3) = 0.0004    !JPL2015; Sessler et al. JGR1996: 0.0006
      !rulesRates%gamNAT_clono2_hcl = 0.20      ! gamNAT(4) = 0.20      !JPL2015; Sessler et al. JGR1996: 0.30
      !rulesRates%gamNAT_clono2_hcl = 1./(4.34 +1./1.4*exp(-0.518*amax1(t-190.,0.))) !Hanson et al., JGR 98, 22931 (1993) !Original buggy value for p(HCl) = 6.7e-8 mbar
      !IUPAC: For ClONO2+HCl on NAT a format 1/gamma = 1/gamma_max + 1/(A*exp(B*(T-190)) where parameters depend on HCl partial pressure.
      !       p(HCl) = 6.7e-8 mbar: gamma_max = 0.23 (or 0.25), A = 0.7022, and B = -0.518
      !       p(HCl) = 4.6e-7 mbar: gamma_max = 0.20, A = 2.2543, and B = -0.558
      !rulesRates%gamNAT_clono2_hcl = 1./( 1.0/0.20 +1./(2.2543*exp(-0.558*amax1(T-190.,0.))) ) ! gamNAT(4)   IUPAC for  p(HCl) = 4.6e-7 mbar      
      rulesRates%gamNAT_clono2_hcl = 1./( 1.0/0.23 +1./(0.7022*exp(-0.518*amax1(T-190.,0.))) )  ! gamNAT(4)   IUPAC for  p(HCl) = 6.7e-8 mbar
      rulesRates%gamNAT_hocl_hcl   = 0.10       ! gamNAT(5) = 0.10      !JPL2015, Sessler et al JGR1996
      rulesRates%gamNAT_brono2_hcl = 0.30       ! gamNAT(6) = 0.30      !Sessler et al.,JGR (1996) 28817.
      rulesRates%gamNAT_hobr_hcl   = 0.25       ! gamNAT(7) = 0.25      !Sessler et al.,JGR (1996) 28817.
      !IUPAC: For N2O5+HCl the HCl coverage again takes a somewhat different value than for the recent June 2016 evaluation for HCl+ice. 
      thetaHCl=7.3e-17*exp(2858/T)*(cncHCl*1e-6*avogadro)/(1.0+7.3e-17*exp(2858/T)*(cncHCl*1e-6*avogadro));
      rulesRates%gamNAT_n2o5_hcl = 4e-3*thetaHCl+6e-4;       ! gamNAT(8) = 0.003   !IUPAC;  JPL: gam=0.003 
      !IUPAC: For ClONO2+HBr on NAT: gamma = 0.56*theta_HBr. Surface coverage theta_HBr = 4.14e-10*[HBr]**0.88 (188 K) where [HBr] is concentration in molecule/cm3.
      thetaHBr = 4.14e-10*(cncHBr*1e-6*avogadro)**0.88
      rulesRates%gamNAT_clono2_hbr = 0.56*min(thetaHBr,1.0)  ! gamNAT(9) = 0.3     !IUPAC: Same value than for ICE. Also, thetaHBr evaluated at 188 K;  JPL2015: gamma>0.3
      rulesRates%gamNAT_hocl_hbr   = 0.10       ! gamNAT(10)= 0.10      !Sessler et al.,JGR (1996) 28817.
      !2Aug2017 by R.H. enabled the following reactions coefficients:
      rulesRates%gamNAT_brono2_hbr = 0.30       ! gamNAT(11)= 0.0(n.a.) !Sessler et al. JGR1996: 0.30
      rulesRates%gamNAT_hobr_hbr   = 0.12       ! gamNAT(12)= 0.0(n.a.) !Sessler et al. JGR1996: 0.12
      !IUPAC: For N2O5+HBr on NAT: The formula seems to allow values larger than one for surface coverage thetaHBr (gamma > 0.02 when [HBr] is large)  
      rulesRates%gamNAT_n2o5_hbr = min(0.02*thetaHBr,0.1)    ! gamNAT(13) = 0.005  !IUPAC but additionally limiting max to 0.1;  JPL: gamma = 0.005
      !!
      !! END of NAT-gamma
      !!
      
      rulesRates%h_clono2_h2o = (rulesRates%gamSTS_clono2_h2o * sadSTS + rulesRates%gamNAT_clono2_h2o * sadNAT + &
                               & rulesRates%gamOT_clono2_h2o * sadOTH)* cf_clono2
      !2Aug2017 by R.H. enabled the following coefficient for reaction (2)
      rulesRates%h_brono2_h2o = (rulesRates%gamSTS_brono2_h2o * sadSTS + rulesRates%gamNAT_brono2_h2o * sadNAT + &
                               & rulesRates%gamOT_brono2_h2o * sadOTH)* cf_brono2
      rulesRates%h_n2o5_h2o   =  (rulesRates%gamSTS_n2o5_h2o * sadSTS  + rulesRates%gamNAT_n2o5_h2o * sadNAT   + &
                               & rulesRates%gamOT_n2o5_h2o * sadOTH)* cf_n2o5
      rulesRates%h_clono2_hcl = (rulesRates%gamSTS_clono2_hcl * sadSTS + rulesRates%gamNAT_clono2_hcl * sadNAT + &
                               & rulesRates%gamOT_clono2_hcl * sadOTH)* cf_clono2
      rulesRates%h_hocl_hcl   = (rulesRates%gamSTS_hocl_hcl * sadSTS   + rulesRates%gamNAT_hocl_hcl * sadNAT + & !BUGGY corrected 14Nov2017: gamNAT_clono2_hcl->gamNAT_hocl_hcl
                               & rulesRates%gamOT_hocl_hcl * sadOTH)* cf_hocl
      !2Aug2017 by R.H. enabled the following 2 coefficients for reactions (6-7)
      rulesRates%h_brono2_hcl = (rulesRates%gamSTS_brono2_hcl * sadSTS+ rulesRates%gamNAT_brono2_hcl * sadNAT + &
                               & rulesRates%gamOT_brono2_hcl * sadOTH)* cf_brono2
      rulesRates%h_hobr_hcl   = (rulesRates%gamSTS_hobr_hcl * sadSTS  + rulesRates%gamNAT_hobr_hcl * sadNAT + &
                               & rulesRates%gamOT_hobr_hcl * sadOTH)* cf_hobr
      rulesRates%h_n2o5_hcl   = (rulesRates%gamSTS_n2o5_hcl * sadSTS  + rulesRates%gamNAT_n2o5_hcl * sadNAT + &
                               & rulesRates%gamOT_n2o5_hcl * sadOTH)* cf_n2o5
      !2Aug2017 by R.H. enabled the following 5 coefficients for reactions (9-13)
      rulesRates%h_clono2_hbr = (rulesRates%gamSTS_clono2_hbr * sadSTS+ rulesRates%gamNAT_clono2_hbr * sadNAT + &
                               & rulesRates%gamOT_clono2_hbr * sadOTH)* cf_clono2
      rulesRates%h_hocl_hbr   = (rulesRates%gamSTS_hocl_hbr * sadSTS  + rulesRates%gamNAT_hocl_hbr * sadNAT + &
                               & rulesRates%gamOT_hocl_hbr * sadOTH)* cf_hocl
      !R.H. 1Aug2017: corrected the typo in the following h_brono2_nbr -> h_brono2_hbr (PLUS SEVERAL OTHER PLACES, TOTALLY 19)
      rulesRates%h_brono2_hbr = (rulesRates%gamSTS_brono2_hbr * sadSTS + rulesRates%gamNAT_brono2_hbr * sadNAT + &
                               & rulesRates%gamOT_brono2_hbr * sadOTH)* cf_brono2
      rulesRates%h_hobr_hbr   = (rulesRates%gamSTS_hobr_hbr * sadSTS + rulesRates%gamNAT_hobr_hbr * sadNAT + &
                               & rulesRates%gamOT_hobr_hbr * sadOTH)* cf_hobr
      !write(*,*) cf_hobr,rulesRates%gamSTS_hobr_hbr,sadSTS,rulesRates%gamNAT_hobr_hbr,sadNAT, '% HOBr+HBr: cf_hobr,gamSTS,sadSTS,gamNAT,sadNAT' !TEST TEST RISTO
      rulesRates%h_n2o5_hbr   = (rulesRates%gamSTS_n2o5_hbr * sadSTS + rulesRates%gamNAT_n2o5_hbr * sadNAT + &
                               & rulesRates%gamOT_n2o5_hbr * sadOTH)* cf_n2o5 
!RISTO Cl&Br: NOTE: Added wider checks for all rates!!!
      if (.not. all((/rulesRates%h_clono2_h2o, rulesRates%h_n2o5_h2o, rulesRates%h_clono2_hcl, rulesRates%h_hocl_hcl, rulesRates%h_n2o5_hcl, & !Cl-reactions (1,3-5,8)
                    & rulesRates%h_brono2_h2o, rulesRates%h_brono2_hcl, rulesRates%h_hobr_hcl, rulesRates%h_clono2_hbr, rulesRates%h_hocl_hbr, & !Br-reactions (2,6,7,9,10)
                    & rulesRates%h_brono2_hbr, rulesRates%h_hobr_hbr, rulesRates%h_n2o5_hbr /) >= 0.)) then !new Br-reactions (11-13)
  call set_error('Trouble with heterogeneous rates in STS or NAT','rates_pseudo_first_order')
  call msg('Rates with STS or NAT1: h_clono2_h2o  h_n2o5_h2o  h_clono2_hcl  h_hocl_hcl  h_n2o5_hcl')
  call msg('Rates with STS or NAT1:',(/rulesRates%h_clono2_h2o, rulesRates%h_n2o5_h2o, rulesRates%h_clono2_hcl, rulesRates%h_hocl_hcl, rulesRates%h_n2o5_hcl/))
  call msg('Rates with STS or NAT1: h_brono2_h2o h_brono2_hcl  h_hobr_hcl  h_clono2_hbr h_hocl_hbr')
  call msg('Rates with STS or NAT1:',(/rulesRates%h_brono2_h2o, rulesRates%h_brono2_hcl, rulesRates%h_hobr_hcl, rulesRates%h_clono2_hbr, rulesRates%h_hocl_hbr/))
  call msg('Rates with STS or NAT1: h_brono2_hbr h_hobr_hbr    h_n2o5_hbr')
  call msg('Rates with STS or NAT1:',(/rulesRates%h_brono2_hbr, rulesRates%h_hobr_hbr, rulesRates%h_n2o5_hbr/))
  
  call msg("rulesRates%gamSTS_clono2_h2o, sadSTS, rulesRates%gamNAT_clono2_h2o, sadNAT, rulesRates%gamOT_clono2_h2o, sadOTH,  cf_clono2")
  call msg("",(/rulesRates%gamSTS_clono2_h2o, sadSTS,  rulesRates%gamNAT_clono2_h2o , sadNAT , rulesRates%gamOT_clono2_h2o , sadOTH,  cf_clono2/))

!!$!RISTO: DEBUG
!!$  write(*,*) 'ifH2O = ', ifH2O
!!$  write(*,*) 'ifHCl = ', ifHCl
!!$  call msg('T, rho1, wt_h2so4_mass_prc')
!!$  call msg("",(/T, rho1, wt_h2so4_mass_prc/))
!!$  call msg('Mh2so4, A, t0, nya, ah, wt_mole_h2so4, ph2o0, ph2o, aw, phcl, Hhcl, Mhcl')
!!$  call msg("",(/Mh2so4, A, t0, nya, ah, wt_mole_h2so4, ph2o0, ph2o, aw, phcl, Hhcl, Mhcl/))
!!$  if (ifH2O .and. ifHCl) then
!!$     call msg('pclono2, Sclono2, Dclono2, kh2o, kh, khydr, khcl, Gbh2o, lclono2, fclono2, Grxn, Gbhcl, Gsp, Gbphcl, Gb, gam')
!!$     call msg("",(/pclono2, Sclono2, Dclono2, kh2o, kh, khydr, khcl, Gbh2o, lclono2, fclono2, Grxn, Gbhcl, Gsp, Gbphcl, Gb, gam/))     
!!$  elseif (ifHCl) then
!!$     !call msg('')
!!$     !call msg("",(/ /))     
!!$  elseif(ifH2O)then
!!$     call msg('pclono2, Sclono2, Hclono2, Dclono2, kh2o, kh, khydr, Gbh2o, lclono2, fclono2, Grxn')
!!$     call msg("",(/pclono2, Sclono2, Hclono2, Dclono2, kh2o, kh, khydr, Gbh2o, lclono2, fclono2, Grxn/))     
!!$  end if
!!$!RISTO: END DEBUG

      endif
#ifdef DEBUG_MORE
if(print_it)then
  call msg('Rates with STS or NAT: h_clono2_h2o  h_n2o5_h2o  h_clono2_hcl  h_hocl_hcl  h_n2o5_hcl')
  call msg('Rates with STS or NAT:',(/rulesRates%h_clono2_h2o, rulesRates%h_n2o5_h2o, rulesRates%h_clono2_hcl, rulesRates%h_hocl_hcl, rulesRates%h_n2o5_hcl/))
  call msg('Rates with STS or NAT: h_brono2_h2o h_brono2_hcl  h_hobr_hcl  h_clono2_hbr h_hocl_hbr')
  call msg('Rates with STS or NAT:',(/rulesRates%h_brono2_h2o, rulesRates%h_brono2_hcl, rulesRates%h_hobr_hcl, rulesRates%h_clono2_hbr, rulesRates%h_hocl_hbr/))
  call msg('Rates with STS or NAT: h_brono2_hbr h_hobr_hbr    h_n2o5_hbr')
  call msg('Rates with STS or NAT:',(/rulesRates%h_brono2_hbr, rulesRates%h_hobr_hbr, rulesRates%h_n2o5_hbr/))
endif
#endif
    
    endif ! ifSTS

    return
  end subroutine rates_pseudo_first_order


   !*******************************************************
  
  function fclono(x) result(fclono2)
      ! Calculates  1./tanh(1/x) - x with correcte assymptotes for x \in [0,\infty)
      implicit none
      real, intent(in) :: x
      real :: fclono2
      character(len = *), parameter :: sub_name = 'fclono'

        if (x > 10.) then
          fclono2 = 1./(3.*x) - 1./(45.*x*x*x)
        elseif (x > 0.1) then
          fclono2 = 1./tanh(1./x) - x
        elseif (x >= 0.) then
          fclono2 =  1. - x
        else
          fclono2 = F_NAN
        endif
  end function fclono
  !************************************************************************************

  logical function fu_if_tla_required_MAAD(rules) result(required)
    ! Collect transformations' requests for tangent linearization. If tangent linear
    ! is not allowed, corresponding subroutine sets error. 
    implicit none
    type(Tchem_rules_AerDynMidAtm), intent(out) :: rules
    
    call set_error('Aerosol dynamics mid-atmosphere does not have TLA/ADJOINT','fu_if_tla_required_MAAD')
    required = .false.
    
  end function fu_if_tla_required_MAAD

END MODULE aer_dyn_middle_atmosph

