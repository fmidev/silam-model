module cbm4_SOA_interface
  use cbm4_SOA_Parameters
  use cbm4_SOA_Integrator, only : rosenbrock, init_solver
  use cbm4_SOA_Global
  use cbm4_SOA_Precision
  use cbm4_SOA_adj_Integrator, only : rosenbrockAdj, init_solver_adj
  use cocktail_basic
  use photolysis
  !$use omp_lib
  implicit none

  private
  
  ! public subroutines - but exported only one level up!
  public rosenbrock
  public set_rates_cbm4_SOA
  public set_const_rates_cbm4_SOA
  public set_fixed_cbm4_SOA
  public register_species_int
  public inventory_int
  public init_solver
  public report_rates_cbm4_SOA

  public init_solver_adj
  public rosenbrockAdj

  ! only this will be exported globally
  public init_chemicals_cbm4_SOA

  
  integer, parameter, private :: num_species = nvar ! from _Parameters module
  integer, parameter, public :: num_species_cbm4_SOA = num_species ! just for laziness
  integer, parameter, public :: num_react_cbm4_SOA = nreact
  integer, parameter, public :: num_fixed_cbm4_SOA = nfix
  integer, parameter, public :: precision_cbm4_SOA = sp ! from Precision module
  character(len=substNmLen), dimension(num_species), parameter, private :: &
       & subst_names = (/ &
                      'NTR   ','AVB0  ','AVB1e0','AVB1e1','AVB1e2','AVB1e3', &
& 'AVB1e4','AVB1e5','AVB1e6','BVB0  ','BVB1e0','BVB1e1', &
& 'BVB1e2','BVB1e3','O1D   ','MEOH  ','ETOH  ','CRO   ', &
& 'PAN   ','H2O2  ','TOL   ','N2O5  ','XYL   ','HONO  ', &
& 'PNA   ','TO2   ','HNO3  ','ROR4  ','CRES  ','MGLY  ', &
& 'CO    ','ETH   ','OPEN  ','XO2N  ','XO2   ','C5H8_2', &
& 'PAR4  ','OLE4  ','HCHO  ','C5H8  ','ISPD  ','ALD2  ', &
& 'O3    ','O     ','HO2   ','NO3   ','C2O3  ','NO    ', &
& 'NO2   ','OH    ' &
                      & /)
  type(silam_species), dimension(num_species), public, target, save :: species_cbm4_SOA

  
  
contains

 
   subroutine set_rates_cbm4_SOA(photoarr, temp, cAir, press, cH2O&
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
    RCONST(2) = (k0_pow(5.7d-34,2.6d0))
    RCONST(3) = (ARRZ(1.4d-12,1310.0d0))
    RCONST(4) = (ARRZ(5.5d-12,-188.0d0))
    RCONST(5) = (k3rd_iup_pow(1.3d-31,1.5d0,2.3d-11,-0.24d0,0.6d0))
    RCONST(6) = (k0_pow(1.0d-31,1.6d0))
    RCONST(7) = (ARRZ(1.4d-13,2470.0d0))
    RCONST(8) = (photoarr(pd_o3))
    RCONST(9) = (photoarr(pd_o3_o1d))
    RCONST(10) = (arrz(0.79*2.0d-11,-130.0d0)+arrz(0.21*3.2d-11,-67.0d0))
    ! RCONST(11) = constant rate coefficient
    RCONST(12) = (ARRZ(1.7d-12,940.0d0))
    RCONST(13) = (2.03d-16*(temp_dp/300.0d0)**4.57*exp(693.0d0/temp_dp))
    RCONST(14) = (photoarr(pd_no3_no2_o))
    RCONST(15) = (photoarr(pd_no3_no_o2))
    RCONST(16) = (ARRZ(1.8d-11,-110.0d0))
    RCONST(17) = (ARRZ(4.5d-14,1260.0d0))
    RCONST(18) = (k3rd_iup_pow(3.6d-30,4.1d0,1.9d-12,0.2d0,0.35d0))
    ! RCONST(19) = constant rate coefficient
    RCONST(20) = (k3rd_iup_exp(1.3d-3,3.5d0,11000.0d0,9.7d14,-0.1d0,11080.0d0,0.35d0))
    RCONST(21) = (ARRZ(3.3d-39,-530.0d0))
    ! RCONST(22) = constant rate coefficient
    RCONST(23) = (k3rd_iup_pow(7.4d-31,2.4d0,3.3d-11,0.3d0,0.81d0))
    RCONST(24) = (photoarr(pd_hono))
    RCONST(25) = (ARRZ(2.5d-12,-260.0d0))
    ! RCONST(26) = constant rate coefficient
    RCONST(27) = (k3rd_iup_pow(3.3d-30,3.0d0,6.0d-11,0.0d0,0.4d0))
    RCONST(28) = (ksum(2.4d-14,-460.0d0,6.5d-34,-1335.0d0,2.7d-17,-2199.0d0))
    RCONST(29) = (ARRZ(3.45d-12,-270.0d0))
    RCONST(30) = (k3rd_iup_pow(1.3d-31,3.1d0,4.0d-12,0.0d0,0.6d0))
    RCONST(31) = (k3rd_iup_exp(4.1d-5,0.0d0,10650.0d0,4.8d15,0.0d0,11170.0d0,0.5d0))
    RCONST(32) = (ARRZ(3.2d-13,-690.0d0))
    RCONST(33) = (kwtf(2.2d-13,-600.0d0,1.4d-21,-2200.0d0))
    RCONST(34) = (kwtf(1.9d-33,-980.0d0,1.4d-21,-2200.0d0))
    RCONST(35) = (photoarr(pd_h2o2))
    RCONST(36) = (ARRZ(2.9d-12,160.0d0))
    RCONST(37) = (1.5E-13+k3rd_jpl_pow(5.9d-33,1.0d0,1.1d-12,-1.3d0))
    RCONST(38) = (ARRZ(5.5d-12,-125.0d0))
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
    RCONST(52) = (ARRZ(2.45d-12,1775.0d0))
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
    RCONST(78) = (ARRZ(1.05d-14,2000.0d0))
    RCONST(79) = (ARRZ(2.95d-12,450.0d0))
    RCONST(80) = (ARRZ(2.6d-12,-365.0d0))
    ! RCONST(81) = constant rate coefficient
    RCONST(82) = (ARRZ(2.6d-12,-365.0d0))
    ! RCONST(83) = constant rate coefficient
    ! RCONST(84) = constant rate coefficient
    ! RCONST(85) = constant rate coefficient
    RCONST(86) = (0.0036*photoarr(pd_ch2chcho))
    ! RCONST(87) = constant rate coefficient
    RCONST(88) = (ARRZ(7.5d-13,-700.0d0))
    RCONST(89) = (ARRZ(7.5d-13,-700.0d0))
    ! RCONST(90) = constant rate coefficient
    ! RCONST(91) = constant rate coefficient
    RCONST(92) = (ARRZ(7.3d-12,620.0d0))
    RCONST(93) = (ARRZ(3.0d-12,-20.0d0))
    RCONST(94) = (ARRZ(4.8d-11,-250d0))
    RCONST(95) = (3.60E-11*ifMTprods)
    RCONST(96) = (ARRZ(1.5d-11,-449.0d0)*ifMTprods)
    RCONST(97) = (ARRZ(1.2d-15,821.0d0)*ifMTprods)
    RCONST(98) = (ARRZ(3.7d-12,-175.0d0)*ifMTprods)
    RCONST(99) = (3.60E-11*(1.0-ifMTprods))
    RCONST(100) = (ARRZ(1.5d-11,-449.0d0)*(1.0-ifMTprods))
    RCONST(101) = (ARRZ(1.2d-15,821.0d0)*(1.0-ifMTprods))
    RCONST(102) = (ARRZ(3.7d-12,-175.0d0)*(1.0-ifMTprods))
    RCONST(103) = (ARRZ(1.8d-12,-340.0d0)*SOA_b)
    RCONST(104) = (ARRZ(1.7d-11,-116.0d0)*SOA_b)
    RCONST(105) = (3.60E-11*SOA_b)
    RCONST(106) = (ARRZ(2.7d-11,-390.0d0)*SOA_b)
    RCONST(107) = (ARRZ(1.03d-14,1995.0d0)*SOA_b)
    RCONST(108) = (ARRZ(3.15d-12,450.0d0)*SOA_b)
    RCONST(109) = (3.60E-11*SOA_b)
    RCONST(110) = (ARRZ(1.5d-11,-449.0d0)*SOA_b)
    RCONST(111) = (ARRZ(1.2d-15,821.0d0)*SOA_b)
    RCONST(112) = (ARRZ(3.7d-12,-175.0d0)*SOA_b)
    RCONST(113) = (ARRZ(1.8d-12,-340.0d0)*(1.0-SOA_b))
    RCONST(114) = (ARRZ(1.7d-11,-116.0d0)*(1.0-SOA_b))
    RCONST(115) = (3.60E-11*(1.0-SOA_b))
    RCONST(116) = (ARRZ(2.7d-11,-390.0d0)*(1.0-SOA_b))
    RCONST(117) = (ARRZ(1.03d-14,1995.0d0)*(1.0-SOA_b))
    RCONST(118) = (ARRZ(3.15d-12,450.0d0)*(1.0-SOA_b))
    RCONST(119) = (3.60E-11*(1.0-SOA_b))
    RCONST(120) = (ARRZ(1.5d-11,-449.0d0)*(1.0-SOA_b))
    RCONST(121) = (ARRZ(1.2d-15,821.0d0)*(1.0-SOA_b))
    RCONST(122) = (ARRZ(3.7d-12,-175.0d0)*(1.0-SOA_b))
    RCONST(123) = (ASOA_aging)
    RCONST(124) = (ASOA_aging)
    RCONST(125) = (ASOA_aging)
    RCONST(126) = (ASOA_aging)
    RCONST(127) = (IVOC_aging)
    RCONST(128) = (IVOC_aging)
    RCONST(129) = (IVOC_aging)
    RCONST(130) = (BSOA_aging)
    RCONST(131) = (BSOA_aging)
    RCONST(132) = (BSOA_aging)
    RCONST(133) = (BSOA_aging)

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

  end subroutine set_rates_cbm4_SOA

  subroutine report_rates_cbm4_SOA()
    implicit none
    
    integer :: ind_react

    do ind_react = 1, nreact
      call msg('Reaction #, rate:', ind_react, rconst(ind_react))
    end do

  end subroutine report_rates_cbm4_SOA

  subroutine set_const_rates_cbm4_SOA()
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
    RCONST(83) = 3.36e-11
    RCONST(84) = 7.1e-18
    RCONST(85) = 1e-15
    RCONST(87) = 1.5e-19
    RCONST(90) = 6.8e-14
    RCONST(91) = 6.8e-14
  end subroutine set_const_rates_cbm4_SOA

  !subroutine set_fixed_cbm4_SOA(cnc_m, cnc_H2O,cnc_O2,cnc_CH4,cnc_M)
  subroutine set_fixed_cbm4_SOA(cnc_H2O,cnc_O2,cnc_CH4,cnc_M)
    implicit none

    !real :: cnc_m
    real, intent(in) :: cnc_H2O,cnc_O2,cnc_CH4,cnc_M
    
    fix(1) = cnc_H2O
fix(2) = cnc_O2
fix(3) = cnc_CH4
fix(4) = cnc_M

  end subroutine set_fixed_cbm4_SOA

  subroutine init_chemicals_cbm4_SOA()
    implicit none

    integer :: i

    do i = 1, num_species_cbm4_SOA
      call set_species(species_cbm4_SOA(i), fu_get_material_ptr(subst_names(i)), in_gas_phase)
      if (error) return
    end do

  end subroutine init_chemicals_cbm4_SOA
  
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

    call addSpecies(speciesTransp, nSpeciesTransp, species_cbm4_SOA, num_species, .true.)
    if(error)return

    !
    ! Presence of emission species in the list of transported cbm species means ownership
    !
    do iEmis = 1, nSpeciesEmis
      if(fu_index(speciesEmis(iEmis), species_cbm4_SOA, num_species) > 0)then
        if(iClaimedSpecies(iEmis) < 0)then
          call msg('cbm4_SOA owns:' + fu_substance_name(speciesEmis(iEmis)))
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
      if (.not. fu_index(species_cbm4_SOA(i), speciesTransp) == i) then
        call msg('Looking for:')
        call report(species_cbm4_SOA(i))
        call msg('Index in cbm4_SOA, transport:', i, fu_index(species_cbm4_SOA(i), speciesTransp))
        call set_error('Incompatible ordering of transport species', &
                     & 'register_species_int')
        return
      end if
      if (fu_name(fu_material(species_cbm4_SOA(i)))=='NO') iNO = i
      if (fu_name(fu_material(species_cbm4_SOA(i)))=='XO2') iXO2 = i
      if (fu_name(fu_material(species_cbm4_SOA(i)))=='HO2') iHO2 = i
    end do  ! species
    if(iNO<0 .or. iXO2<0 .or. iHO2<0)then
        call set_error('Failed to find something', 'register_species_int')
    endif
    
  end subroutine register_species_int

end module cbm4_SOA_interface
