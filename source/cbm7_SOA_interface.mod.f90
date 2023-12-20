module cbm7_SOA_interface
  use cbm7_SOA_Parameters
  use cbm7_SOA_Integrator, only : rosenbrock, init_solver
  use cbm7_SOA_Global
  use cbm7_SOA_Precision
  use cbm7_SOA_adj_Integrator, only : rosenbrockAdj, init_solver_adj
  use cocktail_basic
  use photolysis
  !$use omp_lib
  implicit none

  private
  
  ! public subroutines - but exported only one level up!
  public rosenbrock
  public set_rates_cbm7_SOA
  public set_const_rates_cbm7_SOA
  public set_fixed_cbm7_SOA
  public register_species_int
  public inventory_int
  public init_solver
  public report_rates_cbm7_SOA

  public init_solver_adj
  public rosenbrockAdj

  ! only this will be exported globally
  public init_chemicals_cbm7_SOA

  
  integer, parameter, private :: num_species = nvar ! from _Parameters module
  integer, parameter, public :: num_species_cbm7_SOA = num_species ! just for laziness
  integer, parameter, public :: num_react_cbm7_SOA = nreact
  integer, parameter, public :: num_fixed_cbm7_SOA = nfix
  integer, parameter, public :: precision_cbm7_SOA = sp ! from Precision module
  character(len=substNmLen), dimension(num_species), parameter, private :: &
       & subst_names = (/ &
                      'AVB0   ','BVB0   ','O1D    ','PRPA   ','AVB1e6 ','AVB1e5 ', &
& 'AVB1e4 ','ECH4   ','ETHA   ','ETHY   ','AVB1e0 ','AVB1e1 ', &
& 'AVB1e2 ','AVB1e3 ','BENZENE','N2O5   ','XPRP   ','TOL    ', &
& 'XYL    ','ETOH   ','BVB1e0 ','BVB1e1 ','BVB1e2 ','BVB1e3 ', &
& 'HACT   ','HONO   ','PACD   ','HPLD   ','PAN    ','MEOH   ', &
& 'MEPX   ','FACD   ','ROR7   ','PANX   ','PNA    ','NTR1   ', &
& 'PAR7   ','EPOX   ','CAT1   ','OPAN   ','NTR2   ','INTR   ', &
& 'KET7   ','CRO    ','ISPX   ','ARPX   ','H2O2   ','AACD   ', &
& 'ETH    ','SESQ   ','APIN   ','CRON   ','CO     ','C5H8   ', &
& 'OLE7   ','IOLE   ','XPAR7  ','ALDX   ','XO2H   ','XO2    ', &
& 'SQO2   ','EPX2   ','TPO2   ','CRES   ','BZO2   ','XLO2   ', &
& 'TO2    ','HNO3   ','C5H8_2 ','HCHO   ','APO2   ','XO2N   ', &
& 'ACET   ','ISPD   ','GLYD   ','ALD2   ','ISO2   ','GLY    ', &
& 'O      ','MGLY   ','ROOH   ','XOPN   ','OPEN   ','MEO2   ', &
& 'C2O3   ','CXO3   ','RO2    ','HO2    ','OH     ','NO3    ', &
& 'OPO3   ','O3     ','NO     ','NO2    ' &
                      & /)
  type(silam_species), dimension(num_species), public, target, save :: species_cbm7_SOA

  
  
contains

 
   subroutine set_rates_cbm7_SOA(photoarr, temp, cAir, press, cH2O&
                               &, SOA_b, ASOA_aging, BSOA_aging, IVOC_aging, ifMTprods&
  !                             &, MCF_factor&
                               &)

    implicit none
    real, dimension(:), intent(in) :: photoarr
    real, intent(in) :: temp, cH2O, press, cAir
    real, intent(in) :: SOA_b, ASOA_aging, BSOA_aging, IVOC_aging, ifMTprods
    !real, intent(in) :: MCF_factor
    ! cH2O is abs. humidity in concentration units
    real, parameter :: mol2molec = avogadro*1e-6, press_trp = 200.0e2 !Pa
      
    real(r8k) :: temp_dp

    temp_dp = real(temp, r8k)
        RCONST(1) = (photoarr(pd_no2))
    RCONST(2) = (k0_pow(6.0d-34,2.6d0))
    RCONST(3) = (ARRZ(2.07d-12,1400.0d0))
    RCONST(4) = (k3rd_iup_pow(1.0d-31,1.6d0,5.0d-11,0.3d0,0.85d0))
    RCONST(5) = (ARRZ(5.1d-12,-198.0d0))
    RCONST(6) = (k3rd_iup_pow(1.3d-31,1.5d0,2.3d-11,-0.24d0,0.6d0))
    RCONST(7) = (ARRZ(8.0d-12,2060d0))
    RCONST(8) = (photoarr(pd_o3))
    RCONST(9) = (photoarr(pd_o3_o1d))
    RCONST(10) = (arrz(0.79*2.15d-11,-110.0d0)+arrz(0.21*3.2d-11,-67.0d0))
    ! RCONST(11) = constant rate coefficient
    RCONST(12) = (ARRZ(1.7d-12,940.0d0))
    RCONST(13) = (2.03d-16*(temp_dp/300.0d0)**4.57*exp(693.0d0/temp_dp))
    RCONST(14) = (ARRZ(2.4d-11,-110d0))
    RCONST(15) = (ARRZ(3.0d-11,-200d0))
    RCONST(16) = (6.2d-14*((temp_dp/298.0d0)**2.6)*exp(945.0d0/temp_dp))
    RCONST(17) = (k3rd_iup_pow(9.0d-31,3.2d0,3.9d-11,0.47d0,0.43d0))
    RCONST(18) = (ARRZ(4.8d-11,-250d0))
    RCONST(19) = (kwtf(2.2d-13,-600.0d0,1.4d-21,-2200.0d0))
    RCONST(20) = (kwtf(1.9d-33,-980.0d0,1.4d-21,-2200.0d0))
    RCONST(21) = (photoarr(pd_h2o2))
    ! RCONST(22) = constant rate coefficient
    RCONST(23) = (ARRZ(1.4d-12,2000.0d0))
    RCONST(24) = (ARRZ(4.25d-39,-663.5d0))
    RCONST(25) = (ARRZ(3.45d-12,-270.0d0))
    RCONST(26) = (ARRZ(1.40d-13,2470.0d0))
    RCONST(27) = (photoarr(pd_no3_no2_o))
    RCONST(28) = (photoarr(pd_no3_no_o2))
    RCONST(29) = (ARRZ(1.8d-11,-110.0d0))
    RCONST(30) = (ARRZ(4.5d-14,1260.0d0))
    RCONST(31) = (2.0d-11)
    ! RCONST(32) = constant rate coefficient
    RCONST(33) = (ARRZ(8.5d-13,2450d0))
    RCONST(34) = (k3rd_iup_pow(3.6d-30,4.1d0,1.9d-12,-0.2d0,0.35d0))
    RCONST(35) = (k3rd_iup_exp(1.3d-3,3.5d0,11000.0d0,9.7d14,-0.1d0,11080.0d0,0.35d0))
    RCONST(36) = (photoarr(pd_n2o5))
    ! RCONST(37) = constant rate coefficient
    RCONST(38) = (k3rd_iup_pow(7.4d-31,2.4d0,3.3d-11,0.3d0,0.81d0))
    RCONST(39) = (photoarr(pd_hono))
    RCONST(40) = (ARRZ(2.5d-12,-260.0d0))
    RCONST(41) = (k3rd_iup_pow(3.2d-30,4.5d0,3.0d-11,0.0d0,0.41d0))
    ! RCONST(42) = constant rate coefficient
    RCONST(43) = (ksum(2.4d-14,-460.0d0,6.5d-34,-1335.0d0,2.7d-17,-2199.0d0))
    RCONST(44) = (photoarr(pd_hno3))
    RCONST(45) = (k3rd_iup_pow(1.4d-31,3.1d0,4.0d-12,0.0d0,0.4d0))
    RCONST(46) = (k3rd_iup_exp(4.1d-5,0.0d0,10650.0d0,6.0d15,0.0d0,11170.0d0,0.4d0))
    RCONST(47) = (photoarr(pd_ho2no2))
    RCONST(48) = (photoarr(pd_ho2no2_oh_no3))
    RCONST(49) = (ARRZ(3.2d-13,-690.0d0))
    RCONST(50) = (ARRZ(7.7d-12,2100d0))
    RCONST(51) = (1.5E-13+k3rd_jpl_pow(5.9d-33,1.0d0,1.1d-12,-1.3d0))
    RCONST(52) = (ARRZ(7.5d-12,-290.0d0))
    RCONST(53) = (k3rd_iup_pow(3.28d-28,6.87d0,1.125d-11,1.105d0,0.3d0))
    RCONST(54) = (k3rd_iup_exp(1.1d-5,0.0d0,10100.0d0,1.9d17,0.0d0,14100.0d0,0.3d0))
    RCONST(55) = (photoarr(pd_pan))
    RCONST(56) = (ARRZ(3.14d-12,-580.0d0))
    RCONST(57) = (ARRZ(4.4d-13,-1070.0d0))
    RCONST(58) = (ARRZ(2.9d-12,-500.0d0))
    RCONST(59) = (ARRZ(6.7d-12,-340.0d0))
    RCONST(60) = (1.19d0*k3rd_iup_pow(3.28d-28,6.87d0,1.125d-11,1.105d0,0.3d0))
    RCONST(61) = (1.19d0*k3rd_iup_exp(1.1d-5,0.0d0,10100.0d0,1.9d17,0.0d0,14100.0d0,0.3d0))
    RCONST(62) = (3.00d-12)
    RCONST(63) = (ARRZ(3.14d-12,-580.0d0))
    RCONST(64) = (ARRZ(4.4d-13,-1070.0d0))
    RCONST(65) = (ARRZ(6.7d-12,-340.0d0))
    RCONST(66) = (1.19d0*k3rd_iup_pow(3.28d-28,6.87d0,1.125d-11,1.105d0,0.3d0))
    RCONST(67) = (1.19d0*k3rd_iup_exp(1.1d-5,0.0d0,10100.0d0,1.9d17,0.0d0,14100.0d0,0.3d0))
    RCONST(68) = (3.60d-11)
    RCONST(69) = (ARRZ(3.14d-12,-580.0d0))
    RCONST(70) = (ARRZ(4.4d-13,-1070.0d0))
    RCONST(71) = (ARRZ(2.7d-12,-360.0d0))
    RCONST(72) = (ARRZ(2.12d-13,-1300.0d0))
    RCONST(73) = (ARRZ(1.55d-13,-350.0d0))
    RCONST(74) = (ARRZ(2.3d-12,-360.0d0))
    RCONST(75) = (ARRZ(3.8d-13,-780.0d0))
    RCONST(76) = (ARRZ(2.0d-12,-500.0d0))
    RCONST(77) = (ARRZ(1.55d-13,-350.0d0))
    RCONST(78) = (ARRZ(5.3d-12,-190.0d0))
    RCONST(79) = (photoarr(pd_ch3ooh))
    RCONST(80) = (ARRZ(2.7d-12,-360.0d0))
    RCONST(81) = (ARRZ(2.12d-13,-1300.0d0))
    RCONST(82) = (ARRZ(3.1d-13,-350.0d0))
    RCONST(83) = (ARRZ(2.7d-12,-360.0d0))
    RCONST(84) = (ARRZ(2.12d-13,-1300.0d0))
    RCONST(85) = (ARRZ(3.1d-13,-350.0d0))
    RCONST(86) = (ARRZ(2.7d-12,-360.0d0))
    RCONST(87) = (ARRZ(2.12d-13,-1300.0d0))
    RCONST(88) = (ARRZ(3.1d-13,-350.0d0))
    RCONST(89) = (ARRZ(5.3d-12,-190.0d0))
    RCONST(90) = (photoarr(pd_ch3ooh))
    RCONST(91) = (2.00d-12)
    RCONST(92) = (photoarr(pd_c3h7ono2))
    RCONST(93) = (2.30d-5)
    RCONST(94) = (ARRZ(2.85d-12,345.0d0))
    RCONST(95) = (ARRZ(3.0d-12,-20.0d0))
    RCONST(96) = (ARRZ(5.4d-12,-135.0d0))
    RCONST(97) = (photoarr(pd_hcho_2h))
    RCONST(98) = (photoarr(pd_hcho_h2))
    ! RCONST(99) = constant rate coefficient
    RCONST(100) = (ARRZ(4.7d-12,-345.0d0))
    RCONST(101) = (ARRZ(1.4d-12,1860.0d0))
    RCONST(102) = (photoarr(pd_ald2))
    RCONST(103) = (ARRZ(4.9d-12,-405.0d0))
    RCONST(104) = (6.3d-15)
    RCONST(105) = (photoarr(pd_c2h5cho))
    RCONST(106) = (8.00d-12)
    RCONST(107) = (photoarr(pd_glyd))
    RCONST(108) = (ARRZ(1.4d-12,1860.0d0))
    RCONST(109) = (ARRZ(3.1d-12,-340.0d0))
    RCONST(110) = (photoarr(pd_gly))
    RCONST(111) = (4.0d-16)
    RCONST(112) = (photoarr(pd_mgly))
    RCONST(113) = (5.0d-16)
    RCONST(114) = (ARRZ(1.9d-12,-575.0d0))
    RCONST(115) = (photoarr(pd_acet))
    RCONST(116) = (ARRZ(1.41d-12,620.6d0))
    RCONST(117) = (photoarr(pd_mek))
    RCONST(118) = (1.0d-12)
    RCONST(119) = (ARRZ(2.0d-12,-320.0d0))
    RCONST(120) = (4.5d-13)
    RCONST(121) = (ARRZ(4.0d-14,-850.0d0))
    RCONST(122) = (ARRZ(5.3d-12,-190.0d0))
    RCONST(123) = (ARRZ(1.85d-12,1690.0d0))
    RCONST(124) = (ARRZ(1.85d-12,1690.0d0))
    RCONST(125) = (ARRZ(6.9d-12,1000.0d0))
    RCONST(126) = (ARRZ(7.6d-12,585.0d0))
    RCONST(127) = (k3rd_iup_pow(2.37d-21,0d0,4.3d-1,8.0d0,0.41d0))
    RCONST(128) = (1.0d0)
    RCONST(129) = (3.09d-13*(temp_dp/300.0d0)**2*exp(300.0d0/temp_dp))
    RCONST(130) = (k3rd_iup_pow(4.81d-20,0d0,4.3d-1,8.0d0,0.41d0))
    RCONST(131) = (1.0d0)
    RCONST(132) = (k3rd_iup_exp(5.0d-15,0.0d0,250.0d0,2.4d12,0.0d0,5000.0d0,0.6d0))
    RCONST(133) = (ARRZ(2.00d-14,250.0d0))
    RCONST(134) = (k3rd_iup_pow(5.00d-30,1.5d0,1.0d-12,0.0d0,0.37d0))
    RCONST(135) = (k3rd_iup_pow(8.6d-29,3.1d0,9.0d-12,0.85d0,0.48d0))
    RCONST(136) = (ARRZ(6.82d-15,2500.0d0))
    RCONST(137) = (ARRZ(3.3d-12,2880.0d0))
    RCONST(138) = (k3rd_iup_pow(8.0d-27,3.5d0,3.0d-11,1.0d0,0.5d0))
    RCONST(139) = (ARRZ(5.77d-15,1880.0d0))
    RCONST(140) = (ARRZ(4.6d-13,1155.0d0))
    RCONST(141) = (ARRZ(1.05d-11,-519.0d0))
    RCONST(142) = (ARRZ(4.7d-15,1013.0d0))
    RCONST(143) = (3.7d-13)
    RCONST(144) = (ARRZ(2.30d-12,190.0d0))
    RCONST(145) = (ARRZ(2.7d-12,-360.0d0))
    RCONST(146) = (ARRZ(2.24d-13,-1300.0d0))
    RCONST(147) = (ARRZ(3.1d-13,-350.0d0))
    RCONST(148) = (ARRZ(1.8d-12,-340.0d0))
    RCONST(149) = (ARRZ(2.7d-12,-360.0d0))
    RCONST(150) = (ARRZ(2.39d-13,-1300.0d0))
    RCONST(151) = (ARRZ(3.1d-13,-350.0d0))
    RCONST(152) = (1.51d-11)
    RCONST(153) = (ARRZ(2.7d-12,-360.0d0))
    RCONST(154) = (ARRZ(2.50d-13,-1300.0d0))
    RCONST(155) = (ARRZ(3.1d-13,-350.0d0))
    RCONST(156) = (0.028*photoarr(pd_no2))
    RCONST(157) = (4.4d-11)
    RCONST(158) = (ARRZ(5.4d-17,500.0d0))
    RCONST(159) = (3.80d-12)
    RCONST(160) = (photoarr(pd_c5h6o2))
    RCONST(161) = (9.00d-11)
    RCONST(162) = (ARRZ(1.08d-16,500.0d0))
    RCONST(163) = (3.00d-12)
    RCONST(164) = (ARRZ(1.70d-12,-950.0d0))
    RCONST(165) = (1.40d-11)
    RCONST(166) = (2.10d-12)
    RCONST(167) = (5.50d-12)
    RCONST(168) = (1.53d-12)
    RCONST(169) = (3.80d-12)
    RCONST(170) = (photoarr(pd_cron))
    RCONST(171) = (5.00d-11)
    RCONST(172) = (1.70d-10)
    RCONST(173) = (8.00d-11)
    RCONST(174) = (ARRZ(2.70d-11,-390.0d0))
    RCONST(175) = (ARRZ(2.7d-12,-360.0d0))
    RCONST(176) = (ARRZ(2.12d-13,-1300.0d0))
    RCONST(177) = (ARRZ(3.1d-13,-350.0d0))
    RCONST(178) = (ARRZ(3.3d9,8300.0d0))
    RCONST(179) = (ARRZ(1.05d-14,2000.0d0))
    RCONST(180) = (ARRZ(2.95d-12,450.0d0))
    RCONST(181) = (ARRZ(7.00d-12,-430.0d0))
    RCONST(182) = (ARRZ(3.94d-14,-475.0d0))
    RCONST(183) = (0.0036*photoarr(pd_ch2chcho))
    RCONST(184) = (ARRZ(2.80d-11,-370.0d0))
    RCONST(185) = (0.0036*photoarr(pd_ch2chcho))
    RCONST(186) = (ARRZ(1.17d-11,-450.0d0))
    RCONST(187) = (ARRZ(5.43d-11,450.0d0))
    RCONST(188) = (ARRZ(2.7d-12,-360.0d0))
    RCONST(189) = (ARRZ(2.38d-13,-1300.0d0))
    RCONST(190) = (ARRZ(3.1d-13,-350.0d0))
    RCONST(191) = (ARRZ(1.34d-11,-410.0d0))
    RCONST(192) = (ARRZ(1.00d-11,-300.0d0))
    RCONST(193) = (ARRZ(2.7d-12,-360.0d0))
    RCONST(194) = (ARRZ(2.60d-13,-1300.0d0))
    RCONST(195) = (ARRZ(3.1d-13,-350.0d0))
    RCONST(196) = (ARRZ(8.22d-16,640.0d0))
    RCONST(197) = (ARRZ(1.20d-12,-490.0d0))
    RCONST(198) = (ARRZ(4.07d-11,-350.0d0))
    RCONST(199) = (ARRZ(2.7d-12,-360.0d0))
    RCONST(200) = (ARRZ(2.60d-13,-1300.0d0))
    RCONST(201) = (ARRZ(3.1d-13,-350.0d0))
    RCONST(202) = (ARRZ(4.46d-15,960.0d0))
    RCONST(203) = (7.00d-12)
    RCONST(204) = (2.00d-10)
    RCONST(205) = (ARRZ(2.7d-12,-360.0d0))
    RCONST(206) = (ARRZ(2.60d-13,-1300.0d0))
    RCONST(207) = (ARRZ(3.1d-13,-350.0d0))
    RCONST(208) = (1.20d-14)
    RCONST(209) = (1.90d-11)
    RCONST(210) = (photoarr(pd_o2))
    RCONST(211) = (ARRZ(1.8d-12,-340.0d0)*SOA_b)
    RCONST(212) = (1.51e-11*SOA_b)
    RCONST(213) = (ARRZ(2.3d-12,190.0d0)*SOA_b)
    RCONST(214) = (3.60E-11*SOA_b)
    RCONST(215) = (ARRZ(2.7d-11,-390.0d0)*SOA_b)
    RCONST(216) = (ARRZ(1.05d-14,2000.0d0)*SOA_b)
    RCONST(217) = (ARRZ(2.95d-12,450.0d0)*SOA_b)
    RCONST(218) = (3.60E-11*SOA_b)
    RCONST(219) = (ARRZ(4.07d-11,-350.0d0)*SOA_b)
    RCONST(220) = (ARRZ(4.46d-15,960.0d0)*SOA_b)
    RCONST(221) = (7.00d-12*SOA_b)
    RCONST(222) = (ARRZ(1.8d-12,-340.0d0)*(1.0-SOA_b))
    RCONST(223) = (1.51e-11*(1.0-SOA_b))
    RCONST(224) = (ARRZ(2.3d-12,190.0d0)*(1.0-SOA_b))
    RCONST(225) = (3.60E-11*(1.0-SOA_b))
    RCONST(226) = (ARRZ(2.7d-11,-390.0d0)*(1.0-SOA_b))
    RCONST(227) = (ARRZ(1.05d-14,2000.0d0)*(1.0-SOA_b))
    RCONST(228) = (ARRZ(2.95d-12,450.0d0)*(1.0-SOA_b))
    RCONST(229) = (3.60E-11*(1.0-SOA_b))
    RCONST(230) = (ARRZ(4.07d-11,-350.0d0)*(1.0-SOA_b))
    RCONST(231) = (ARRZ(4.46d-15,960.0d0)*(1.0-SOA_b))
    RCONST(232) = (7.00d-12*(1.0-SOA_b))
    ! RCONST(233) = constant rate coefficient
    ! RCONST(234) = constant rate coefficient
    ! RCONST(235) = constant rate coefficient
    RCONST(236) = (ASOA_aging)
    RCONST(237) = (ASOA_aging)
    RCONST(238) = (ASOA_aging)
    RCONST(239) = (ASOA_aging)
    RCONST(240) = (IVOC_aging)
    RCONST(241) = (IVOC_aging)
    RCONST(242) = (IVOC_aging)
    RCONST(243) = (BSOA_aging)
    RCONST(244) = (BSOA_aging)
    RCONST(245) = (BSOA_aging)
    RCONST(246) = (BSOA_aging)

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

  end subroutine set_rates_cbm7_SOA

  subroutine report_rates_cbm7_SOA()
    implicit none
    
    integer :: ind_react

    do ind_react = 1, nreact
      call msg('Reaction #, rate:', ind_react, rconst(ind_react))
    end do

  end subroutine report_rates_cbm7_SOA

  subroutine set_const_rates_cbm7_SOA()
    implicit none
        RCONST(11) = 2.14e-10
    RCONST(22) = 1.8e-12
    RCONST(32) = 4e-12
    RCONST(37) = 1e-22
    RCONST(42) = 1.1e-30
    RCONST(99) = 5.5e-16
    RCONST(233) = 2e-10
    RCONST(234) = 1.2e-14
    RCONST(235) = 1.9e-11
  end subroutine set_const_rates_cbm7_SOA

  !subroutine set_fixed_cbm7_SOA(cnc_m, cnc_H2O,cnc_H2,cnc_O2,cnc_CH4,cnc_M)
  subroutine set_fixed_cbm7_SOA(cnc_H2O,cnc_H2,cnc_O2,cnc_CH4,cnc_M)
    implicit none

    !real :: cnc_m
    real, intent(in) :: cnc_H2O,cnc_H2,cnc_O2,cnc_CH4,cnc_M
    
    fix(1) = cnc_H2O
fix(2) = cnc_H2
fix(3) = cnc_O2
fix(4) = cnc_CH4
fix(5) = cnc_M

  end subroutine set_fixed_cbm7_SOA

  subroutine init_chemicals_cbm7_SOA()
    implicit none

    integer :: i

    do i = 1, num_species_cbm7_SOA
      call set_species(species_cbm7_SOA(i), fu_get_material_ptr(subst_names(i)), in_gas_phase)
      if (error) return
    end do

  end subroutine init_chemicals_cbm7_SOA
  
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

    call addSpecies(speciesTransp, nSpeciesTransp, species_cbm7_SOA, num_species, .true.)
    if(error)return

    !
    ! Presence of emission species in the list of transported cbm species means ownership
    !
    do iEmis = 1, nSpeciesEmis
      if(fu_index(speciesEmis(iEmis), species_cbm7_SOA, num_species) > 0)then
        if(iClaimedSpecies(iEmis) < 0)then
          call msg('cbm7_SOA owns:' + fu_substance_name(speciesEmis(iEmis)))
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
  !                            & nspeciesTransp, nspeciesShortlived, nspeciesAerosol)                     !Without SOA
                              & nspeciesTransp, nspeciesShortlived, nspeciesAerosol, iNO, iXO2, iHO2)     !With SOA
    !
    ! The integrator requires fixed indices for the species. At this
    ! point, we just check that the transport species follow this ordering.
    implicit none
    type(silam_species), dimension(:), pointer :: speciesTransp, speciesShortlived, &
                                                & speciesAerosol
    integer, intent(in) :: nspeciesTransp, nspeciesShortlived, nspeciesAerosol
    integer, intent(out) :: iNO, iXO2, iHO2
    
    ! Local variable
    integer :: i
    iNO = int_missing
    iXO2 = int_missing
    iHO2 = int_missing

    do i = 1, num_species
      if (.not. fu_index(species_cbm7_SOA(i), speciesTransp) == i) then
        call msg('Looking for:')
        call report(species_cbm7_SOA(i))
        call msg('Index in cbm7_SOA, transport:', i, fu_index(species_cbm7_SOA(i), speciesTransp))
        call set_error('Incompatible ordering of transport species', &
                     & 'register_species_int')
        return
      end if
      if (fu_name(fu_material(species_cbm7_SOA(i)))=='NO') iNO = i
      if (fu_name(fu_material(species_cbm7_SOA(i)))=='XO2') iXO2 = i
      if (fu_name(fu_material(species_cbm7_SOA(i)))=='HO2') iHO2 = i
    end do  ! species
    if(iNO<0 .or. iXO2<0 .or. iHO2<0)then
        call set_error('Failed to find something', 'register_species_int')
    endif
    
  end subroutine register_species_int

end module cbm7_SOA_interface
