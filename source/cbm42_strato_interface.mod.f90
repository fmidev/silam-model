module cbm42_strato_interface
  use cbm42_strato_Parameters
  use cbm42_strato_Integrator, only : rosenbrock, init_solver
  use cbm42_strato_Global
  use cbm42_strato_Precision
  use cbm42_strato_adj_Integrator, only : rosenbrockAdj, init_solver_adj
  use cocktail_basic
  use photolysis
  !$use omp_lib
  implicit none

  private
  
  ! public subroutines - but exported only one level up!
  public rosenbrock
  public set_rates_cbm42_strato
  public set_const_rates_cbm42_strato
  public set_fixed_cbm42_strato
  public register_species_int
  public inventory_int
  public init_solver
  public report_rates_cbm42_strato

  public init_solver_adj
  public rosenbrockAdj

  ! only this will be exported globally
  public init_chemicals_cbm42_strato

  
  integer, parameter, private :: num_species = nvar ! from _Parameters module
  integer, parameter, public :: num_species_cbm42_strato = num_species ! just for laziness
  integer, parameter, public :: num_react_cbm42_strato = nreact
  integer, parameter, public :: num_fixed_cbm42_strato = nfix
  integer, parameter, public :: precision_cbm42_strato = sp ! from Precision module
  character(len=substNmLen), dimension(num_species), parameter, private :: &
       & subst_names = (/ &
                      'Br2    ','BrNO2  ','CCl3F  ','CCl2F2 ','CCl4   ','CH3CCl3', &
& 'PAN    ','CRO    ','BrONO2 ','BrCl   ','TOL    ','CH3Cl  ', &
& 'CH3Br  ','Cl2O2  ','N2O5   ','XYL    ','XO2N   ','N      ', &
& 'N2O    ','HONO   ','TO2    ','PNA    ','ClNO2  ','HOBr   ', &
& 'ROR    ','MGLY   ','H2O2   ','CRES   ','HNO3   ','CO     ', &
& 'Cl2    ','XO2    ','PAR    ','OPEN   ','ETH    ','ClONO2 ', &
& 'HOCl   ','HCl    ','HBr    ','C2O3   ','C5H8   ','OLE    ', &
& 'ALD2   ','H      ','OClO   ','NO3    ','NO2    ','O      ', &
& 'OH     ','NO     ','O3     ','HO2    ','HCHO   ','ClO    ', &
& 'O1D    ','BrO    ','Br     ','Cl     ' &
                      & /)
  type(silam_species), dimension(num_species), public, target, save :: species_cbm42_strato

  
  
contains

 
  subroutine set_rates_cbm42_strato(photoarr, temp, cAir, press, cH2O)
  !subroutine set_rates_cbm42_strato(photoarr, temp, cAir, press, cH2O, &
  !                             & SOA_b, ASOA_aging, BSOA_aging, IVOC_aging, ifMTprods)

    implicit none
    real, dimension(:), intent(in) :: photoarr
    real, intent(in) :: temp, cH2O, press, cAir
    !real, intent(in) :: SOA_b, ASOA_aging, BSOA_aging, IVOC_aging, ifMTprods
    ! cH2O is abs. humidity in concentration units
    real, parameter :: mol2molec = avogadro*1e-6, press_trp = 200.0e2 !Pa
      
    real(r8k) :: temp_dp

    temp_dp = real(temp, r8k)
        RCONST(1) = (photoarr(pd_no2))
    RCONST(2) = (k0_pow(6.0d-34,2.6d0))
    RCONST(3) = (ARRZ(2.07d-12,1400.0d0))
    RCONST(4) = (ARRZ(5.1d-12,-198.0d0))
    RCONST(5) = (k3rd_iup_pow(1.3d-31,1.5d0,2.3d-11,-0.24d0,0.6d0))
    RCONST(6) = (k3rd_iup_pow(1.0d-31,1.6d0,5.0d-11,0.3d0,0.85d0))
    RCONST(7) = (ARRZ(1.4d-13,2470.0d0))
    RCONST(8) = (photoarr(pd_o3))
    RCONST(9) = (photoarr(pd_o3_o1d))
    RCONST(10) = (arrz(0.79*2.15d-11,-110.0d0)+arrz(0.21*3.2d-11,-67.0d0))
    ! RCONST(11) = constant rate coefficient
    RCONST(12) = (ARRZ(1.7d-12,940.0d0))
    RCONST(13) = (2.03d-16*(temp_dp/300.0d0)**4.57*exp(693.0d0/temp_dp))
    RCONST(14) = (photoarr(pd_no3_no2_o))
    RCONST(15) = (photoarr(pd_no3_no_o2))
    RCONST(16) = (ARRZ(1.8d-11,-110.0d0))
    RCONST(17) = (ARRZ(4.5d-14,1260.0d0))
    RCONST(18) = (k3rd_iup_pow(3.6d-30,4.1d0,1.9d-12,-0.2d0,0.35d0))
    ! RCONST(19) = constant rate coefficient
    RCONST(20) = (k3rd_iup_exp(1.3d-3,3.5d0,11000.0d0,9.7d14,-0.1d0,11080.0d0,0.35d0))
    RCONST(21) = (ARRZ(4.25d-39,-663.5d0))
    ! RCONST(22) = constant rate coefficient
    RCONST(23) = (k3rd_iup_pow(7.4d-31,2.4d0,3.3d-11,0.3d0,0.81d0))
    RCONST(24) = (photoarr(pd_hono))
    RCONST(25) = (ARRZ(2.5d-12,-260.0d0))
    ! RCONST(26) = constant rate coefficient
    RCONST(27) = (k3rd_iup_pow(3.3d-30,4.5d0,3.0d-11,0.0d0,0.41d0))
    RCONST(28) = (ksum(2.4d-14,-460.0d0,6.5d-34,-1335.0d0,2.7d-17,-2199.0d0))
    RCONST(29) = (ARRZ(3.45d-12,-270.0d0))
    RCONST(30) = (k3rd_iup_pow(1.4d-31,3.1d0,4.0d-12,0.0d0,0.4d0))
    RCONST(31) = (k3rd_iup_exp(4.1d-5,0.0d0,10650.0d0,6.0d15,0.0d0,11170.0d0,0.4d0))
    RCONST(32) = (ARRZ(3.2d-13,-690.0d0))
    RCONST(33) = (kwtf(2.2d-13,-600.0d0,1.4d-21,-2200.0d0))
    RCONST(34) = (kwtf(1.9d-33,-980.0d0,1.4d-21,-2200.0d0))
    RCONST(35) = (photoarr(pd_h2o2))
    RCONST(36) = (ARRZ(2.9d-12,160.0d0))
    RCONST(37) = (1.44d-13*(1.0d0+real(cAir,dp)/4.2d19))
    RCONST(38) = (ARRZ(5.4d-12,-135.0d0))
    RCONST(39) = (photoarr(pd_hcho_2h))
    RCONST(40) = (photoarr(pd_hcho_h2))
    RCONST(41) = (ARRZ(3.4d-11,1600.0d0))
    ! RCONST(42) = constant rate coefficient
    RCONST(43) = (ARRZ(1.8d-11,1100.0d0))
    RCONST(44) = (ARRZ(4.7d-12,-345.0d0))
    RCONST(45) = (ARRZ(1.4d-12,1860.0d0))
    RCONST(46) = (photoarr(pd_ald2))
    RCONST(47) = (ARRZ(7.5d-12,-290.0d0))
    RCONST(48) = (k3rd_iup_pow(2.7d-28,7.1d0,1.2d-11,0.9d0,0.3d0))
    RCONST(49) = (k3rd_iup_exp(4.9d-3,0.0d0,12100.0d0,5.4d16,0.0d0,13830.0d0,0.3d0))
    RCONST(50) = (ARRZ(2.9d-12,-500.0d0))
    RCONST(51) = (ARRZ(5.2d-13,-980.0d0))
    RCONST(52) = (ARRZ(1.85d-12,1690.0d0))
    RCONST(53) = (ARRZ(6.9d-12,1000.0d0))
    RCONST(54) = (ARRZ(1.0d+15,8000.0d0))
    ! RCONST(55) = constant rate coefficient
    ! RCONST(56) = constant rate coefficient
    RCONST(57) = (ARRZ(1.0d-11,280.0d0))
    RCONST(58) = (ARRZ(5.2d-12,-504.0d0))
    RCONST(59) = (ARRZ(1.4d-14,2105.0d0))
    ! RCONST(60) = constant rate coefficient
    RCONST(61) = (ARRZ(1.0d-11,792.0d0))
    RCONST(62) = (k3rd_iup_pow(8.6d-29,3.1d0,9.0d-12,0.85d0,0.48d0))
    RCONST(63) = (ARRZ(9.1d-15,2580.0d0))
    RCONST(64) = (ARRZ(1.8d-12,-340.0d0))
    ! RCONST(65) = constant rate coefficient
    ! RCONST(66) = constant rate coefficient
    ! RCONST(67) = constant rate coefficient
    ! RCONST(68) = constant rate coefficient
    ! RCONST(69) = constant rate coefficient
    RCONST(70) = (ARRZ(1.7d-11,-116.0d0))
    ! RCONST(71) = constant rate coefficient
    RCONST(72) = (photoarr(pd_open))
    RCONST(73) = (ARRZ(5.4d-17,500.0d0))
    RCONST(74) = (ARRZ(1.9d-12,-575.0d0))
    RCONST(75) = (photoarr(pd_mgly))
    ! RCONST(76) = constant rate coefficient
    RCONST(77) = (ARRZ(2.7d-11,-390.0d0))
    RCONST(78) = (ARRZ(1.03d-14,1995.0d0))
    RCONST(79) = (ARRZ(3.15d-12,450.0d0))
    RCONST(80) = (ARRZ(2.6d-12,-365.0d0))
    ! RCONST(81) = constant rate coefficient
    RCONST(82) = (ARRZ(2.6d-12,-365.0d0))
    RCONST(83) = (ARRZ(2.8d-11,250d0))
    RCONST(84) = (ARRZ(2.5d-11,-110d0))
    RCONST(85) = (ARRZ(6.2d-12,-295d0))
    RCONST(86) = (ARRZ(6.6d-12,1240d0))
    RCONST(87) = (ARRZ(3.9d-11,2310d0))
    RCONST(88) = (4.4E-11-ARRZ(7.48d-11,620d0))
    RCONST(89) = (ARRZ(6.132d-12,-300d0))
    RCONST(90) = (ARRZ(8.1d-11,30d0))
    RCONST(91) = (ARRZ(1.7d-12,230d0))
    RCONST(92) = (k3rd_iup_pow(1.6d-31,3.4d0,7.0d-11,0.0d0,0.4d0))
    RCONST(93) = (ARRZ(4.5d-12,900d0))
    RCONST(94) = (ARRZ(2.2d-12,-340d0))
    ! RCONST(95) = constant rate coefficient
    ! RCONST(96) = constant rate coefficient
    RCONST(97) = (k3rd_jpl_pow(1.8d-31,2.0d0,1.0d-10,1.0d0))
    RCONST(98) = (ARRZ(3.4d-12,130d0))
    RCONST(99) = (ARRZ(4.38d-13,-300d0))
    RCONST(100) = (ARRZ(3.5d-13,1370d0))
    RCONST(101) = (ARRZ(1.0d-12,1590d0))
    RCONST(102) = (ARRZ(3.0d-11,2450d0))
    RCONST(103) = (k3rd_iup_pow(2.0d-32,4.0d0,1.0d-11,0.0d0,0.45d0))
    RCONST(104) = (k3rd_iup_exp(3.7d-7,0.0d0,7690d0,1.8d14,0.0d0,7690d0,0.45d0))
    RCONST(105) = (ARRZ(1.4d-12,-600d0))
    RCONST(106) = (ARRZ(3.2d-11,-170d0))
    RCONST(107) = (ARRZ(2.4d-12,960d0))
    RCONST(108) = (ARRZ(2.1d-12,4700d0))
    RCONST(109) = (ARRZ(1.1d-13,350d0))
    ! RCONST(110) = constant rate coefficient
    RCONST(111) = (ARRZ(7.6d-11,-65.0d0))
    ! RCONST(112) = constant rate coefficient
    ! RCONST(113) = constant rate coefficient
    ! RCONST(114) = constant rate coefficient
    RCONST(115) = (ARRZ(3.6d-12,1200d0))
    RCONST(116) = (ARRZ(6.2d-12,-145d0))
    RCONST(117) = (ARRZ(7.48d-11,620d0))
    RCONST(118) = (ARRZ(1.1d-11,980d0))
    RCONST(119) = (ARRZ(1.0d-11,3300d0))
    RCONST(120) = (ARRZ(1.2d-12,330d0))
    RCONST(121) = (ARRZ(2.4d-12,1250d0))
    RCONST(122) = (ARRZ(8.0d-12,2060d0))
    ! RCONST(123) = constant rate coefficient
    RCONST(124) = (0.78*k0_pow(2.8d-36,0.9d0))
    ! RCONST(125) = constant rate coefficient
    RCONST(126) = (ARRZ(2.7d-11,-224d0))
    RCONST(127) = (ARRZ(4.8d-11,-250d0))
    RCONST(128) = (ARRZ(1.4d-12,2000.0d0))
    RCONST(129) = (6.2d-14*((temp_dp/298.0d0)**2.6)*exp(945.0d0/temp_dp))
    RCONST(130) = (k3rd_iup_pow(9.0d-31,3.2d0,3.9d-11,0.47d0,0.43d0))
    RCONST(131) = (1.25d-10)
    RCONST(132) = (1.7d-11)
    RCONST(133) = (2.0d-11)
    ! RCONST(134) = constant rate coefficient
    RCONST(135) = (ARRZ(8.5d-13,2450d0))
    RCONST(136) = (photoarr(pd_o2))
    RCONST(137) = (photoarr(pd_n2o))
    RCONST(138) = (photoarr(pd_n2o5))
    RCONST(139) = (photoarr(pd_hno3))
    RCONST(140) = (photoarr(pd_hocl))
    RCONST(141) = (photoarr(pd_clono2_cl_no3))
    RCONST(142) = (photoarr(pd_clono2))
    RCONST(143) = (photoarr(pd_oclo))
    RCONST(144) = (photoarr(pd_cl2o2))
    RCONST(145) = (photoarr(pd_hcl))
    RCONST(146) = (photoarr(pd_cl2))
    RCONST(147) = (photoarr(pd_clno2))
    RCONST(148) = (ARRZ(1.7d-11,800d0))
    RCONST(149) = (ARRZ(1.9d-11,-230d0))
    RCONST(150) = (ARRZ(8.7d-12,-260d0))
    RCONST(151) = (ARRZ(1.6d-12,-430d0))
    RCONST(152) = (ARRZ(2.9d-12,-220d0))
    RCONST(153) = (ARRZ(5.8d-13,-170d0))
    RCONST(154) = (ARRZ(1.6d-12,-210d0))
    RCONST(155) = (ARRZ(7.7d-12,450d0))
    RCONST(156) = (ARRZ(2.7d-11,1300d0))
    RCONST(157) = (ARRZ(1.7d-11,800d0))
    RCONST(158) = (ARRZ(6.7d-12,-155d0))
    RCONST(159) = (k3rd_iup_pow(4.7d-31,3.1d0,1.8d-11,0.0d0,0.40d0))
    RCONST(160) = (ARRZ(4.5d-12,-500d0))
    ! RCONST(161) = constant rate coefficient
    RCONST(162) = (ARRZ(1.7d-11,-250d0))
    RCONST(163) = (ARRZ(5.8d-12,1500d0))
    RCONST(164) = (ARRZ(1.2d-10,430d0))
    RCONST(165) = (photoarr(pd_brono2))
    RCONST(166) = (photoarr(pd_brcl))
    RCONST(167) = (photoarr(pd_hobr))
    RCONST(168) = (photoarr(pd_bro))
    RCONST(169) = (photoarr(pd_br2))
    RCONST(170) = (photoarr(pd_brno2))
    RCONST(171) = (photoarr(pd_ho2no2))
    RCONST(172) = (photoarr(pd_ho2no2_oh_no3))
    RCONST(173) = (photoarr(pd_h2o))
    RCONST(174) = (photoarr(pd_ch4))
    RCONST(175) = (photoarr(pd_co2))
    RCONST(176) = (photoarr(pd_no))
    RCONST(177) = (ARRZ(2.4d-11,-110d0))
    RCONST(178) = (ARRZ(7.7d-12,2100d0))
    RCONST(179) = (1.5d-10)
    RCONST(180) = (ARRZ(2.1d-11,-100d0))
    RCONST(181) = (ARRZ(1.5d-11,3600d0))
    RCONST(182) = (ARRZ(5.8d-12,-220d0))
    RCONST(183) = (k3rd_iup_pow(4.3d-32,1.2d0,9.5d-11,-0.44d0,0.5d0))
    RCONST(184) = (ARRZ(1.4d-10,470d0))
    RCONST(185) = (8.0d-11)
    RCONST(186) = (ARRZ(4.0d-10,340d0))
    RCONST(187) = (photoarr(pd_ch3cl))
    RCONST(188) = (ARRZ(2.1d-12,1210d0))
    RCONST(189) = (photoarr(pd_ch3br))
    RCONST(190) = (ARRZ(1.7d-12,1215d0))
    RCONST(191) = (photoarr(pd_cfc11))
    RCONST(192) = (photoarr(pd_cfc12))
    RCONST(193) = (photoarr(pd_ccl4))
    RCONST(194) = (photoarr(pd_ch3ccl3))
    RCONST(195) = (ARRZ(1.2d-12,1440d0))

  contains

    REAL(kind=sp) FUNCTION ARR2( A0,B0 )
      REAL A0,B0           
      ARR2 =  DBLE(A0) * EXP( DBLE(B0)/TEMP )              
    END FUNCTION ARR2

    real function photo(l, m, n, cos_theta)
      implicit none
      real, intent(in) :: l, n, m, cos_theta
      if (cos_theta < 1e-5) then
        photo = 0.0
        return
      end if
      photo = l * cos_theta**(m) * exp(-n / cos_theta)
    end function photo
    
!!$    real function phux(x, y, z) result(rate)
!!$      implicit none
!!$      real, intent(in) :: x, y, z
!!$      real, parameter :: minyz = -30.0, eminyz=9.357623e-14
!!$      
!!$      real :: chiz, ychiz, eychiz, chi
!!$      
!!$      chi = acos(cos_theta)
!!$      if (cos_theta > 0) then
!!$        chiz = chi*z
!!$      else
!!$        chiz = chi
!!$      end if
!!$      
!!$      if (chiz < pi/2-1e-5) then
!!$        ychiz = y * (1 - (1/cos(chiz)))
!!$        if (ychiz > minyz) then
!!$          eychiz = exp(ychiz)
!!$        else
!!$          eychiz = eminyz
!!$        end if
!!$      else
!!$        eychiz = eminyz
!!$      end if
!!$      rate = x * eychiz * cloud_att
!!$      !print *, cos_theta, chi, chiz, rate, x, eychiz, ychiz
!!$
!!$    end function phux
    	
    real(kind=sp) function arrz(a, b) result(coef)
      implicit none
      real(r8k), intent(in) :: a, b

      coef = a * exp(-b / temp_dp)
      
    end function arrz

    REAL(kind=sp) FUNCTION k3rd_jpl_pow(k0_300K, n, kinf_300K, m) result(rate)
    ! rate constants in form of JPL (fc = 0.6, N = 1 in IUPAC fit), with 
    ! asymptotic values given as k0 = k0_300 * (T/300)**-n [M] and similar for kinf.
      !REAL(kind=sp), INTENT(IN) :: temp      ! temperature [K]
      !REAL(kind=sp), INTENT(IN) :: cair      ! air concentration [molecules/cm3]
      REAL(r8k), INTENT(IN) :: k0_300K   ! low pressure limit at 300 K
      REAL(r8k), INTENT(IN) :: n         ! exponent for low pressure limit
      REAL(r8k), INTENT(IN) :: kinf_300K ! high pressure limit at 300 K
      REAL(r8k), INTENT(IN) :: m         ! exponent for high pressure limit

      REAL(kind=r8k) :: zt_help, k0_T, kinf_T, k_ratio
      
      real, parameter :: fc = 0.6 
      
      zt_help = 300.0_dp/temp_dp
      k0_T    = k0_300K   * zt_help**(n) * cair ! k_0   at current T
      kinf_T  = kinf_300K * zt_help**(m)        ! k_inf at current T
      k_ratio = k0_T/kinf_T
      rate   = k0_T/(1.0_dp+k_ratio)*fc**(1.0_dp/(1.0_dp+LOG10(k_ratio)**2))
      
    END FUNCTION k3rd_jpl_pow

    REAL(kind=sp) FUNCTION k3rd_iup_pow(k0_300K, n, kinf_300K, m, fc) result(rate)
    ! rate constants as defined by IUPAC. Variable Fc, N defined as in IUPAC. 
    ! Limiting cases with power expressions, as above.
      !REAL(kind=sp), INTENT(IN) :: temp      ! temperature [K]
      !REAL(kind=sp), INTENT(IN) :: cair      ! air concentration [molecules/cm3]
      REAL(r8k), INTENT(IN) :: k0_300K   ! low pressure limit at 300 K
      REAL(r8k), INTENT(IN) :: n         ! exponent for low pressure limit
      REAL(r8k), INTENT(IN) :: kinf_300K ! high pressure limit at 300 K
      REAL(r8k), INTENT(IN) :: m         ! exponent for high pressure limit
      REAL(r8k), INTENT(IN) :: fc        ! broadening factor (usually fc=0.6)
      REAL(kind=r8k) :: zt_help, k0_T, kinf_T, k_ratio
      real(r8k) :: z, Nbig
      
      zt_help = 300.0_dp/temp_dp
      k0_T    = k0_300K   * zt_help**(n) * cair ! k_0   at current T
      kinf_T  = kinf_300K * zt_help**(m)        ! k_inf at current T
      k_ratio = k0_T/kinf_T
      Nbig = 0.75_r8k -1.27_r8k*log10(fc)
      z = 1_r8k + (log10(k_ratio)/Nbig)**2
      rate  = k0_T/(1._dp+k_ratio)*fc**(1_dp/z)
      
    END FUNCTION k3rd_iup_pow

    REAL(kind=sp) FUNCTION k3rd_iup_exp(k0_300K, n, b1, kinf_300K, m, b2, fc) result(rate)
      !
      ! Reaction rates in form 
      ! k0 = k0_300k * [M] * (T/300)^-n * exp(-b1/T)
      ! kinf = kinf_300k * (T/300)^-m *exp(-b2/T)
      ! 
      !REAL(kind=sp), INTENT(IN) :: temp      ! temperature [K]
      !REAL(kind=sp), INTENT(IN) :: cair      ! air concentration [molecules/cm3]
      REAL(r8k), INTENT(IN) :: k0_300K   ! low pressure limit at 300 K
      REAL(r8k), INTENT(IN) :: n         ! exponent for low pressure limit
      REAL(r8k), INTENT(IN) :: kinf_300K ! high pressure limit at 300 K
      REAL(r8k), INTENT(IN) :: m         ! exponent for high pressure limit
      REAL(r8k), INTENT(IN) :: fc        ! broadening factor (usually fc=0.6)
      real(r8k), intent(in) :: b1, b2    ! the factors inside exponents
      REAL(kind=r8k) :: zt_help, k0_T, kinf_T, k_ratio
      real(r8k) :: z, Nbig

      zt_help = 300._dp/temp_dp
      k0_T    = k0_300K   * zt_help**(n) * cair * exp(-b1/temp_dp)! k_0   at current T
      kinf_T  = kinf_300K * zt_help**(m) * exp(-b2/temp_dp)       ! k_inf at current T
      k_ratio = k0_T/kinf_T
      Nbig = 0.75_r8k -1.27_r8k*log10(fc)
      z = 1_r8k + (log10(k_ratio)/Nbig)**2
      rate   = k0_T/(1._dp+k_ratio)*fc**(1/z)
      
    END FUNCTION k3rd_iup_exp
    

    real(kind=sp) function k0_pow(k0_300K, n) result(k)
      implicit none
      real(r8k), intent(in) :: k0_300K, n
      
      real(r8k) :: zt_help
      
      zt_help = 300._dp/temp_dp
      k = k0_300K   * zt_help**(n) * cAir ! k_0   at current T
      
    end function k0_pow

    real(kind=sp) function ksum(a1, b1, a3, b3, a4, b4) result(k)
      ! Expression used by IUPAC for the CB4 reaction 27 (HNO3+OH+M->NO3), with
      ! k = k1(T) + k2(M, T) where
      ! k1, k3, k4 = arrhenius-like expressions and
      ! k2(M,T) = k3 * [M] / (1+k3*[M]/k4)
      implicit none
      real(kind=r8k), intent(in) :: a1, b1, a3, b3, a4, b4
      
      real(kind=r8k) :: k1, k2, k3, k4
  
      k1 = arrz(a1, b1)
      k3 = arrz(a3, b3) ! note: k3 has different unit
      k4 = arrz(a4, b4)
      
      k2 = k3*cAir / (1.0_r8k + k3*cAir / k4)
      
      k = k1 + k2
    end function ksum

    real(kind=sp) function kwtf(a, b, a_wet, b_wet) result(k)
      implicit none
      ! Rate expression for the cb4 reaction 32 (2 HO2 -> H2O2).
      ! k = arrz(a, b) * (1 + [H2O]*arrz(a_wet, b_wet))
      real(kind=r8k), intent(in) :: a, b, a_wet, b_wet
      
      real :: k_wet
      
      k_wet = arrz(a_wet, b_wet)
      k = arrz(a, b) * (1.0d0 + cH2O*k_wet)
      
    end function kwtf

!!$    real(kind=sp) function photo_var(j_0, j_25) result(rate)
!!$      ! Compute photolysis rate based on the values in Jacobson. Assume:
!!$      ! constant rate j_0 until tropopause
!!$      ! then pressure-linear increase up to 25 hPa
!!$      ! then constant rate j_25
!!$      ! the rate is (incorrectly) scaled with cosine of zenith angle,
!!$      ! cloud_att is applied up to 700 hPa.
!!$      real(kind=sp), intent(in) :: j_0, j_25
!!$
!!$      real :: weight_srf
!!$      
!!$      if (press > press_trp) then
!!$        rate = j_0
!!$      else if (press > 2500.0) then
!!$        weight_srf = (press - 2500) / (press_trp - 2500)
!!$        rate = j_0 * weight_srf + j_25 * (1 - weight_srf)
!!$      else
!!$        rate = j_25
!!$      end if
!!$      rate = rate * cos_theta
!!$      if (press > 70000) rate = rate*cloud_att
!!$      
!!$    end function photo_var

  end subroutine set_rates_cbm42_strato

  subroutine report_rates_cbm42_strato()
    implicit none
    
    integer :: ind_react

    do ind_react = 1, nreact
      call msg('Reaction #, rate:', ind_react, rconst(ind_react))
    end do

  end subroutine report_rates_cbm42_strato

  subroutine set_const_rates_cbm42_strato()
    implicit none
        RCONST(11) = 2.14e-10
    RCONST(19) = 1e-22
    RCONST(22) = 4.41e-40
    RCONST(26) = 1e-20
    RCONST(42) = 5.5e-16
    RCONST(55) = 1600
    RCONST(56) = 1.5e-11
    RCONST(60) = 7.7e-15
    RCONST(65) = 8.1e-12
    RCONST(66) = 4.2
    RCONST(67) = 4.1e-11
    RCONST(68) = 2.2e-11
    RCONST(69) = 1.4e-11
    RCONST(71) = 3e-11
    RCONST(76) = 3.6e-11
    RCONST(81) = 6.8e-14
    RCONST(95) = 5e-13
    RCONST(96) = 1.7e-13
    RCONST(110) = 2.7e-10
    RCONST(112) = 2.4e-11
    RCONST(113) = 3.4e-13
    RCONST(114) = 1.5e-10
    RCONST(123) = 2.4e-10
    RCONST(125) = 1.2e-10
    RCONST(134) = 4e-12
    RCONST(161) = 1.5e-10
  end subroutine set_const_rates_cbm42_strato

  !subroutine set_fixed_cbm42_strato(cnc_m, cnc_H2O,cnc_H2,cnc_O2,cnc_CH4,cnc_CO2,cnc_M)
  subroutine set_fixed_cbm42_strato(cnc_H2O,cnc_H2,cnc_O2,cnc_CH4,cnc_CO2,cnc_M)
    implicit none

    !real :: cnc_m
    real, intent(in) :: cnc_H2O,cnc_H2,cnc_O2,cnc_CH4,cnc_CO2,cnc_M
    
    fix(1) = cnc_H2O
fix(2) = cnc_H2
fix(3) = cnc_O2
fix(4) = cnc_CH4
fix(5) = cnc_CO2
fix(6) = cnc_M

  end subroutine set_fixed_cbm42_strato

  subroutine init_chemicals_cbm42_strato()
    implicit none

    integer :: i

    do i = 1, num_species_cbm42_strato
      call set_species(species_cbm42_strato(i), fu_get_material_ptr(subst_names(i)), in_gas_phase)
      if (error) return
    end do

  end subroutine init_chemicals_cbm42_strato
  
  subroutine inventory_int(flag_transf, &
                         & speciesEmis, speciesTransp, speciesShortlived, speciesAerosol,&
                         & nSpeciesEmis, nSpeciesTransp, nSpeciesShortlived, nSpeciesAerosol, &
                         & iClaimedSpecies)
    implicit none
    integer, intent(in) :: flag_transf ! the transformation id defined in an upper module
    type(silam_species), dimension(:), pointer :: speciesEmis, speciesTransp, speciesShortlived, &
                                                & speciesAerosol
    integer, intent(in) :: nSpeciesEmis
    integer, intent(inout) :: nspeciesTransp, nspeciesShortlived, nspeciesAerosol
    integer, dimension(:), intent(inout) :: iClaimedSpecies

    ! Local variables
    integer :: iEmis

    call addSpecies(speciesTransp, nSpeciesTransp, species_cbm42_strato, num_species, .true.)
    if(error)return

    !
    ! Presence of emission species in the list of transported cbm species means ownership
    !
    do iEmis = 1, nSpeciesEmis
      if(fu_index(speciesEmis(iEmis), species_cbm42_strato, num_species) > 0)then
        if(iClaimedSpecies(iEmis) < 0)then
          call msg('cbm42_strato owns:' + fu_substance_name(speciesEmis(iEmis)))
          iClaimedSpecies(iEmis) = flag_transf
        else
          call msg('Cannot claim ownership because he claimed it already:',iClaimedSpecies(iEmis))
          call set_error('Cannot claim ownership because someone claimed it already',&
                       & 'inventory_int')
          return
        endif
      endif  ! emission species belongs to my transport
    end do  ! emission species

  end subroutine inventory_int
  
    !************************************************************************************

  subroutine register_species_int(speciesTransp, speciesShortlived, speciesAerosol,&
                              & nspeciesTransp, nspeciesShortlived, nspeciesAerosol)                     !Without SOA
  !                            & nspeciesTransp, nspeciesShortlived, nspeciesAerosol, iNO, iXO2, iHO2)     !With SOA
    !
    ! The integrator requires fixed indices for the species. At this
    ! point, we just check that the transport species follow this ordering.
    implicit none
    type(silam_species), dimension(:), pointer :: speciesTransp, speciesShortlived, &
                                                & speciesAerosol
    integer, intent(in) :: nspeciesTransp, nspeciesShortlived, nspeciesAerosol
    !integer, intent(out) :: iNO, iXO2, iHO2
    
    ! Local variable
    integer :: i
    !iNO = int_missing
    !iXO2 = int_missing
    !iHO2 = int_missing

    do i = 1, num_species
      if (.not. fu_index(species_cbm42_strato(i), speciesTransp) == i) then
        call msg('Looking for:')
        call report(species_cbm42_strato(i))
        call msg('Index in cbm42_strato, transport:', i, fu_index(species_cbm42_strato(i), speciesTransp))
        call set_error('Incompatible ordering of transport species', &
                     & 'register_species_int')
        return
      end if
      !if (fu_name(fu_material(species_cbm42_strato(i)))=='NO') iNO = i
      !if (fu_name(fu_material(species_cbm42_strato(i)))=='XO2') iXO2 = i
      !if (fu_name(fu_material(species_cbm42_strato(i)))=='HO2') iHO2 = i
    end do  ! species
    !if(iNO<0 .or. iXO2<0 .or. iHO2<0)then
    !    call set_error('Failed to find something', 'register_species_int')
    !endif
    
  end subroutine register_species_int

end module cbm42_strato_interface
