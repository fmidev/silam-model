module cbm7_strato_SOA_interface
  use cbm7_strato_SOA_Parameters
  use cbm7_strato_SOA_Integrator, only : rosenbrock, init_solver
  use cbm7_strato_SOA_Global
  use cbm7_strato_SOA_Precision
  use cbm7_strato_SOA_adj_Integrator, only : rosenbrockAdj, init_solver_adj
  use cocktail_basic
  use photolysis
  !$use omp_lib
  implicit none

  private
  
  ! public subroutines - but exported only one level up!
  public rosenbrock
  public set_rates_cbm7_strato_SOA
  public set_const_rates_cbm7_strato_SOA
  public set_fixed_cbm7_strato_SOA
  public register_species_int
  public inventory_int
  public init_solver
  public report_rates_cbm7_strato_SOA

  public init_solver_adj
  public rosenbrockAdj

  ! only this will be exported globally
  public init_chemicals_cbm7_strato_SOA

  
  integer, parameter, private :: num_species = nvar ! from _Parameters module
  integer, parameter, public :: num_species_cbm7_strato_SOA = num_species ! just for laziness
  integer, parameter, public :: num_react_cbm7_strato_SOA = nreact
  integer, parameter, public :: num_fixed_cbm7_strato_SOA = nfix
  integer, parameter, public :: precision_cbm7_strato_SOA = sp ! from Precision module
  character(len=substNmLen), dimension(num_species), parameter, private :: &
       & subst_names = (/ &
                      'CCl3F  ','CCl2F2 ','CCl4   ','AVB0   ','BVB0   ','PRPA   ', &
& 'CH3CCl3','AVB1e6 ','AVB1e5 ','AVB1e4 ','ECH4   ','BrNO2  ', &
& 'ETHY   ','ClO2   ','CH3Br  ','AVB1e0 ','AVB1e1 ','AVB1e2 ', &
& 'AVB1e3 ','BENZENE','XPRP   ','Cl2O2  ','N2O5   ','XYL    ', &
& 'TOL    ','N      ','BVB1e0 ','BVB1e1 ','BVB1e2 ','BVB1e3 ', &
& 'HACT   ','N2O    ','HONO   ','PACD   ','HPLD   ','PAN    ', &
& 'MEPX   ','FACD   ','PANX   ','XPAR7  ','ClNO2  ','BrCl   ', &
& 'Br2    ','PNA    ','NTR1   ','CH3Cl  ','EPOX   ','CAT1   ', &
& 'ETOH   ','OPAN   ','ETHA   ','BrONO2 ','NTR2   ','INTR   ', &
& 'KET7   ','ISPX   ','CRO    ','FMCL   ','ARPX   ','MEOH   ', &
& 'HOBr   ','Cl2    ','AACD   ','SESQ   ','ClONO2 ','HOCl   ', &
& 'ROR7   ','APIN   ','H2O2   ','CRON   ','ETH    ','CO     ', &
& 'HBr    ','SQO2   ','EPX2   ','TPO2   ','ROOH   ','BZO2   ', &
& 'XLO2   ','TO2    ','CRES   ','ISO2   ','HNO3   ','APO2   ', &
& 'ACET   ','IOLE   ','PAR7   ','C5H8_2 ','ISPD   ','XOPN   ', &
& 'GLY    ','OLE7   ','MGLY   ','XO2N   ','OPEN   ','XO2    ', &
& 'C5H8   ','XO2H   ','GLYD   ','HCl    ','OClO   ','OPO3   ', &
& 'H      ','CXO3   ','C2O3   ','ALDX   ','ALD2   ','HCHO   ', &
& 'O1D    ','ClO    ','NO2    ','O3     ','MEO2   ','OH     ', &
& 'O      ','Br     ','NO3    ','NO     ','RO2    ','BrO    ', &
& 'HO2    ','Cl     ' &
                      & /)
  type(silam_species), dimension(num_species), public, target, save :: species_cbm7_strato_SOA

  
  
contains

 
   subroutine set_rates_cbm7_strato_SOA(photoarr, temp, cAir, press, cH2O&
                               &, SOA_b, ASOA_aging, BSOA_aging, IVOC_aging, ifMTprods&
                               &, MCF_factor&
                               &)

    implicit none
    real, dimension(:), intent(in) :: photoarr
    real, intent(in) :: temp, cH2O, press, cAir
    real, intent(in) :: SOA_b, ASOA_aging, BSOA_aging, IVOC_aging, ifMTprods
    real, intent(in) :: MCF_factor
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
    RCONST(211) = (photoarr(pd_cl2))
    RCONST(212) = (photoarr(pd_hocl))
    RCONST(213) = (ARRZ(2.8d-11,250d0))
    RCONST(214) = (ARRZ(1.0d-12,1590d0))
    RCONST(215) = (ARRZ(3.0d-11,2450d0))
    RCONST(216) = (ARRZ(6.2d-12,-295d0))
    RCONST(217) = (ARRZ(2.2d-12,-340d0))
    ! RCONST(218) = constant rate coefficient
    RCONST(219) = (photoarr(pd_fmcl))
    RCONST(220) = (ARRZ(6.6d-12,1240d0))
    ! RCONST(221) = constant rate coefficient
    RCONST(222) = (ARRZ(7.05d-11,60d0))
    RCONST(223) = (k3rd_iup_pow(1.85d-29,3.3d0,6.0d-10,0.0d0,0.4d0))
    RCONST(224) = (k3rd_iup_pow(4.0d-28,0.0d0,2.7d-10,0.0d0,0.6d0))
    ! RCONST(225) = constant rate coefficient
    ! RCONST(226) = constant rate coefficient
    RCONST(227) = (ARRZ(8.1d-11,30d0))
    ! RCONST(228) = constant rate coefficient
    ! RCONST(229) = constant rate coefficient
    RCONST(230) = (ARRZ(7.1d-11,75d0))
    RCONST(231) = (ARRZ(6.0d-11,-155d0))
    RCONST(232) = (ARRZ(1.7d-12,230d0))
    RCONST(233) = (ARRZ(1.3d-11,360d0))
    ! RCONST(234) = constant rate coefficient
    ! RCONST(235) = constant rate coefficient
    RCONST(236) = (ARRZ(2.1d-11,-240d0))
    ! RCONST(237) = constant rate coefficient
    RCONST(238) = (ARRZ(2.65d-12,-1600.0d0))
    ! RCONST(239) = constant rate coefficient
    RCONST(240) = (ARRZ(2.5d-11,-110d0))
    RCONST(241) = (ARRZ(3.9d-11,2310d0))
    RCONST(242) = (4.4E-11-ARRZ(7.48d-11,620d0))
    RCONST(243) = (ARRZ(7.3d-12,-300d0))
    RCONST(244) = (k3rd_iup_pow(1.6d-31,3.4d0,7.0d-11,0.0d0,0.4d0))
    RCONST(245) = (ARRZ(4.5d-12,900d0))
    ! RCONST(246) = constant rate coefficient
    ! RCONST(247) = constant rate coefficient
    RCONST(248) = (k3rd_jpl_pow(1.8d-31,2.0d0,1.0d-10,1.0d0))
    RCONST(249) = (ARRZ(3.4d-12,130d0))
    ! RCONST(250) = constant rate coefficient
    RCONST(251) = (ARRZ(3.5d-13,1370d0))
    RCONST(252) = (k3rd_iup_pow(2.0d-32,4.0d0,1.0d-11,0.0d0,0.45d0))
    RCONST(253) = (k3rd_iup_exp(3.7d-7,0.0d0,7690d0,1.8d14,0.0d0,7690d0,0.45d0))
    RCONST(254) = (ARRZ(1.4d-12,-600d0))
    RCONST(255) = (ARRZ(3.2d-11,-170d0))
    RCONST(256) = (ARRZ(2.4d-12,960d0))
    RCONST(257) = (ARRZ(2.1d-12,4700d0))
    RCONST(258) = (ARRZ(1.1d-13,350d0))
    ! RCONST(259) = constant rate coefficient
    RCONST(260) = (ARRZ(7.6d-11,-65.0d0))
    ! RCONST(261) = constant rate coefficient
    ! RCONST(262) = constant rate coefficient
    ! RCONST(263) = constant rate coefficient
    RCONST(264) = (ARRZ(3.6d-12,1200d0))
    RCONST(265) = (ARRZ(6.2d-12,-145d0))
    RCONST(266) = (ARRZ(7.48d-11,620d0))
    RCONST(267) = (ARRZ(1.1d-11,980d0))
    RCONST(268) = (ARRZ(1.0d-11,3300d0))
    RCONST(269) = (ARRZ(1.2d-12,330d0))
    RCONST(270) = (ARRZ(2.4d-12,1250d0))
    ! RCONST(271) = constant rate coefficient
    RCONST(272) = (0.78*k0_pow(2.8d-36,0.9d0))
    RCONST(273) = (1.25d-10)
    RCONST(274) = (ARRZ(2.8d-10,1820d0))
    RCONST(275) = (0.79*1.4d-33*(300.0d0/temp_dp)**3.9+0.21*1.6d-33*(300.0d0/temp_dp)**2.9)
    RCONST(276) = (photoarr(pd_n2o))
    RCONST(277) = (photoarr(pd_clono2_cl_no3))
    RCONST(278) = (photoarr(pd_clono2))
    RCONST(279) = (photoarr(pd_oclo))
    RCONST(280) = (photoarr(pd_cl2o2))
    RCONST(281) = (photoarr(pd_hcl))
    RCONST(282) = (photoarr(pd_clno2))
    RCONST(283) = (photoarr(pd_cloo))
    RCONST(284) = (ARRZ(1.7d-11,800d0))
    RCONST(285) = (ARRZ(1.9d-11,-230d0))
    RCONST(286) = (ARRZ(8.7d-12,-260d0))
    RCONST(287) = (ARRZ(1.6d-12,-430d0))
    RCONST(288) = (ARRZ(2.9d-12,-220d0))
    RCONST(289) = (ARRZ(5.8d-13,-170d0))
    ! RCONST(290) = constant rate coefficient
    RCONST(291) = (ARRZ(2.9d-14,-840d0))
    RCONST(292) = (ARRZ(7.7d-12,450d0))
    RCONST(293) = (ARRZ(2.7d-11,1300d0))
    RCONST(294) = (ARRZ(1.7d-11,800d0))
    RCONST(295) = (ARRZ(6.7d-12,-155d0))
    RCONST(296) = (k3rd_iup_pow(4.7d-31,3.1d0,1.8d-11,0.0d0,0.40d0))
    RCONST(297) = (k3rd_iup_pow(4.2d-31,2.4d0,2.7d-11,0.0d0,0.55d0))
    RCONST(298) = (ARRZ(4.5d-12,-500d0))
    ! RCONST(299) = constant rate coefficient
    RCONST(300) = (ARRZ(1.8d-11,-250d0))
    RCONST(301) = (ARRZ(5.8d-12,1500d0))
    RCONST(302) = (ARRZ(1.2d-10,430d0))
    RCONST(303) = (photoarr(pd_brono2))
    RCONST(304) = (photoarr(pd_brcl))
    RCONST(305) = (photoarr(pd_hobr))
    RCONST(306) = (photoarr(pd_bro))
    RCONST(307) = (photoarr(pd_br2))
    RCONST(308) = (photoarr(pd_brno2))
    RCONST(309) = (photoarr(pd_clo))
    RCONST(310) = (photoarr(pd_h2o))
    RCONST(311) = (photoarr(pd_ch4))
    RCONST(312) = (photoarr(pd_co2))
    RCONST(313) = (photoarr(pd_no))
    RCONST(314) = (1.5d-10)
    RCONST(315) = (ARRZ(2.1d-11,-100d0))
    RCONST(316) = (ARRZ(1.5d-11,3600d0))
    RCONST(317) = (ARRZ(5.8d-12,-220d0))
    RCONST(318) = (k3rd_iup_pow(4.3d-32,1.2d0,9.5d-11,-0.44d0,0.5d0))
    RCONST(319) = (ARRZ(1.4d-10,470d0))
    RCONST(320) = (8.0d-11)
    RCONST(321) = (ARRZ(4.0d-10,340d0))
    RCONST(322) = (photoarr(pd_ch3cl))
    RCONST(323) = (ARRZ(2.1d-12,1210d0))
    RCONST(324) = (ARRZ(1.8d-11,1081d0))
    RCONST(325) = (photoarr(pd_ch3br))
    RCONST(326) = (ARRZ(1.7d-12,1215d0))
    RCONST(327) = (photoarr(pd_cfc11))
    RCONST(328) = (photoarr(pd_cfc12))
    RCONST(329) = (photoarr(pd_ccl4))
    RCONST(330) = (photoarr(pd_ch3ccl3))
    RCONST(331) = (ARRZ(1.2d-12,1440d0)*MCF_factor)
    RCONST(332) = (ARRZ(1.8d-12,-340.0d0)*SOA_b)
    RCONST(333) = (1.51e-11*SOA_b)
    RCONST(334) = (ARRZ(2.3d-12,190.0d0)*SOA_b)
    RCONST(335) = (3.60E-11*SOA_b)
    RCONST(336) = (ARRZ(2.7d-11,-390.0d0)*SOA_b)
    RCONST(337) = (ARRZ(1.05d-14,2000.0d0)*SOA_b)
    RCONST(338) = (ARRZ(2.95d-12,450.0d0)*SOA_b)
    RCONST(339) = (3.60E-11*SOA_b)
    RCONST(340) = (ARRZ(4.07d-11,-350.0d0)*SOA_b)
    RCONST(341) = (ARRZ(4.46d-15,960.0d0)*SOA_b)
    RCONST(342) = (7.00d-12*SOA_b)
    RCONST(343) = (ARRZ(1.8d-12,-340.0d0)*(1.0-SOA_b))
    RCONST(344) = (1.51e-11*(1.0-SOA_b))
    RCONST(345) = (ARRZ(2.3d-12,190.0d0)*(1.0-SOA_b))
    RCONST(346) = (3.60E-11*(1.0-SOA_b))
    RCONST(347) = (ARRZ(2.7d-11,-390.0d0)*(1.0-SOA_b))
    RCONST(348) = (ARRZ(1.05d-14,2000.0d0)*(1.0-SOA_b))
    RCONST(349) = (ARRZ(2.95d-12,450.0d0)*(1.0-SOA_b))
    RCONST(350) = (3.60E-11*(1.0-SOA_b))
    RCONST(351) = (ARRZ(4.07d-11,-350.0d0)*(1.0-SOA_b))
    RCONST(352) = (ARRZ(4.46d-15,960.0d0)*(1.0-SOA_b))
    RCONST(353) = (7.00d-12*(1.0-SOA_b))
    ! RCONST(354) = constant rate coefficient
    ! RCONST(355) = constant rate coefficient
    ! RCONST(356) = constant rate coefficient
    RCONST(357) = (ASOA_aging)
    RCONST(358) = (ASOA_aging)
    RCONST(359) = (ASOA_aging)
    RCONST(360) = (ASOA_aging)
    RCONST(361) = (IVOC_aging)
    RCONST(362) = (IVOC_aging)
    RCONST(363) = (IVOC_aging)
    RCONST(364) = (BSOA_aging)
    RCONST(365) = (BSOA_aging)
    RCONST(366) = (BSOA_aging)
    RCONST(367) = (BSOA_aging)

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

  end subroutine set_rates_cbm7_strato_SOA

  subroutine report_rates_cbm7_strato_SOA()
    implicit none
    
    integer :: ind_react

    do ind_react = 1, nreact
      call msg('Reaction #, rate:', ind_react, rconst(ind_react))
    end do

  end subroutine report_rates_cbm7_strato_SOA

  subroutine set_const_rates_cbm7_strato_SOA()
    implicit none
        RCONST(11) = 2.14e-10
    RCONST(22) = 1.8e-12
    RCONST(32) = 4e-12
    RCONST(37) = 1e-22
    RCONST(42) = 1.1e-30
    RCONST(99) = 5.5e-16
    RCONST(218) = 5e-13
    RCONST(221) = 5e-11
    RCONST(225) = 3.5e-10
    RCONST(226) = 4.3e-10
    RCONST(228) = 8e-11
    RCONST(229) = 1.3e-10
    RCONST(234) = 4.9e-11
    RCONST(235) = 1.6e-11
    RCONST(237) = 9.7e-12
    RCONST(239) = 1.7e-12
    RCONST(246) = 5e-13
    RCONST(247) = 1.7e-13
    RCONST(250) = 2.1e-11
    RCONST(259) = 2.7e-10
    RCONST(261) = 2.4e-11
    RCONST(262) = 4.6e-13
    RCONST(263) = 1.5e-10
    RCONST(271) = 2.4e-10
    RCONST(290) = 2.7e-12
    RCONST(299) = 1.5e-10
    RCONST(354) = 2e-10
    RCONST(355) = 1.2e-14
    RCONST(356) = 1.9e-11
  end subroutine set_const_rates_cbm7_strato_SOA

  !subroutine set_fixed_cbm7_strato_SOA(cnc_m, cnc_H2O,cnc_H2,cnc_O2,cnc_CH4,cnc_CO2,cnc_M)
  subroutine set_fixed_cbm7_strato_SOA(cnc_H2O,cnc_H2,cnc_O2,cnc_CH4,cnc_CO2,cnc_M)
    implicit none

    !real :: cnc_m
    real, intent(in) :: cnc_H2O,cnc_H2,cnc_O2,cnc_CH4,cnc_CO2,cnc_M
    
    fix(1) = cnc_H2O
fix(2) = cnc_H2
fix(3) = cnc_O2
fix(4) = cnc_CH4
fix(5) = cnc_CO2
fix(6) = cnc_M

  end subroutine set_fixed_cbm7_strato_SOA

  subroutine init_chemicals_cbm7_strato_SOA()
    implicit none

    integer :: i

    do i = 1, num_species_cbm7_strato_SOA
      call set_species(species_cbm7_strato_SOA(i), fu_get_material_ptr(subst_names(i)), in_gas_phase)
      if (error) return
    end do

  end subroutine init_chemicals_cbm7_strato_SOA
  
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

    call addSpecies(speciesTransp, nSpeciesTransp, species_cbm7_strato_SOA, num_species, .true.)
    if(error)return

    !
    ! Presence of emission species in the list of transported cbm species means ownership
    !
    do iEmis = 1, nSpeciesEmis
      if(fu_index(speciesEmis(iEmis), species_cbm7_strato_SOA, num_species) > 0)then
        if(iClaimedSpecies(iEmis) < 0)then
          call msg('cbm7_strato_SOA owns:' + fu_substance_name(speciesEmis(iEmis)))
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
      if (.not. fu_index(species_cbm7_strato_SOA(i), speciesTransp) == i) then
        call msg('Looking for:')
        call report(species_cbm7_strato_SOA(i))
        call msg('Index in cbm7_strato_SOA, transport:', i, fu_index(species_cbm7_strato_SOA(i), speciesTransp))
        call set_error('Incompatible ordering of transport species', &
                     & 'register_species_int')
        return
      end if
      if (fu_name(fu_material(species_cbm7_strato_SOA(i)))=='NO') iNO = i
      if (fu_name(fu_material(species_cbm7_strato_SOA(i)))=='XO2') iXO2 = i
      if (fu_name(fu_material(species_cbm7_strato_SOA(i)))=='HO2') iHO2 = i
    end do  ! species
    if(iNO<0 .or. iXO2<0 .or. iHO2<0)then
        call set_error('Failed to find something', 'register_species_int')
    endif
    
  end subroutine register_species_int

end module cbm7_strato_SOA_interface
