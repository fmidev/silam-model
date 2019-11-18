! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! 
! Sparse Hessian Data Structures File
! 
! Generated by KPP-2.2.3 symbolic chemistry Kinetics PreProcessor
!       (http://www.cs.vt.edu/~asandu/Software/KPP)
! KPP is distributed under GPL, the general public licence
!       (http://www.gnu.org/copyleft/gpl.html)
! (C) 1995-1997, V. Damian & A. Sandu, CGRER, Univ. Iowa
! (C) 1997-2005, A. Sandu, Michigan Tech, Virginia Tech
!     With important contributions from:
!        M. Damian, Villanova University, USA
!        R. Sander, Max-Planck Institute for Chemistry, Mainz, Germany
! 
! File                 : cbm4_SOA_HessianSP.f90
! Time                 : Fri Sep 21 17:33:56 2018
! Working directory    : /home/hanniner/tmp/silam_v5_6-test/kpp/cbm4_SOA
! Equation file        : cbm4_SOA.kpp
! Output root filename : cbm4_SOA
! 
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



MODULE cbm4_SOA_HessianSP

  PUBLIC
  SAVE


! Hessian Sparse Data
! 

  INTEGER, PARAMETER, DIMENSION(345) :: IHESS_I = (/ &
       1,  2,  2,  2,  2,  3,  3,  3,  3,  4,  4,  4, &
       4,  5,  5,  5,  5,  6,  6,  7,  7,  8,  9, 10, &
      10, 10, 10, 10, 10, 10, 10, 10, 10, 11, 11, 11, &
      11, 11, 11, 11, 11, 11, 11, 12, 12, 12, 12, 12, &
      12, 12, 12, 12, 12, 13, 13, 13, 13, 13, 15, 16, &
      16, 16, 17, 17, 18, 19, 20, 21, 21, 21, 21, 22, &
      22, 23, 23, 23, 24, 24, 24, 24, 24, 24, 24, 24, &
      24, 25, 25, 25, 25, 25, 26, 26, 26, 26, 27, 27, &
      28, 28, 28, 28, 28, 29, 29, 29, 29, 29, 29, 29, &
      29, 29, 29, 29, 29, 29, 30, 30, 30, 30, 31, 31, &
      31, 31, 31, 31, 31, 31, 31, 31, 31, 31, 31, 31, &
      31, 31, 31, 31, 31, 31, 31, 31, 31, 32, 32, 32, &
      32, 32, 32, 33, 33, 33, 33, 34, 34, 34, 34, 34, &
      34, 34, 34, 34, 34, 34, 35, 35, 35, 35, 35, 35, &
      35, 35, 35, 35, 35, 35, 35, 35, 35, 35, 35, 35, &
      35, 36, 36, 36, 36, 36, 37, 37, 37, 37, 38, 38, &
      38, 38, 38, 38, 38, 38, 38, 38, 38, 38, 38, 38, &
      38, 38, 38, 39, 39, 39, 39, 39, 39, 39, 39, 39, &
      40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 41, &
      41, 41, 41, 41, 41, 41, 41, 41, 41, 41, 41, 41, &
      41, 41, 41, 41, 41, 41, 41, 41, 41, 41, 41, 41, &
      41, 41, 41, 41, 41, 41, 42, 42, 42, 42, 42, 42, &
      42, 42, 42, 42, 42, 42, 42, 42, 42, 42, 42, 42, &
      42, 42, 42, 42, 43, 43, 43, 43, 43, 43, 43, 43, &
      43, 43, 43, 43, 43, 43, 43, 43, 43, 43, 43, 43, &
      43, 43, 43, 43, 43, 43, 43, 43, 43, 43, 43, 43, &
      44, 44, 44, 44, 44, 44, 44, 44, 44, 44, 44, 45, &
      45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, &
      45, 46, 46, 46, 46, 46, 46, 46, 46 /)

  INTEGER, PARAMETER, DIMENSION(345) :: IHESS_J = (/ &
       2,  2,  3, 18, 20,  3,  4, 18, 20,  4,  5, 18, &
      20,  5,  6, 18, 20,  6,  7,  7,  8,  8, 10, 10, &
      11, 33, 33, 33, 33, 37, 37, 37, 37, 11, 12, 33, &
      33, 33, 33, 37, 37, 37, 37, 12, 13, 33, 33, 33, &
      33, 37, 37, 37, 37, 13, 33, 33, 33, 33, 42, 16, &
      26, 26, 17, 43, 18, 40, 20, 21, 21, 41, 42, 22, &
      42, 18, 20, 23, 24, 33, 33, 33, 34, 36, 36, 37, &
      37, 25, 26, 35, 38, 41, 18, 20, 26, 26, 27, 34, &
      20, 28, 30, 37, 37, 29, 30, 30, 32, 32, 33, 35, &
      35, 35, 36, 36, 37, 37, 23, 26, 30, 30, 18, 20, &
      26, 28, 30, 30, 31, 31, 32, 32, 33, 33, 33, 34, &
      36, 36, 36, 36, 37, 37, 43, 44, 44, 32, 32, 32, &
      37, 37, 37, 33, 33, 33, 33, 20, 33, 33, 33, 34, &
      36, 36, 36, 36, 37, 37, 30, 30, 32, 32, 32, 33, &
      33, 35, 35, 35, 36, 36, 36, 36, 37, 37, 43, 44, &
      44, 36, 36, 36, 36, 37, 37, 37, 37, 37, 30, 32, &
      33, 33, 33, 33, 34, 36, 36, 36, 36, 37, 37, 37, &
      38, 38, 38, 30, 32, 33, 36, 37, 39, 39, 39, 39, &
      25, 26, 33, 35, 36, 37, 38, 39, 40, 40, 42, 17, &
      18, 20, 21, 22, 25, 26, 28, 29, 30, 30, 32, 32, &
      33, 33, 34, 35, 35, 36, 36, 36, 37, 37, 38, 38, &
      39, 39, 41, 41, 43, 43, 16, 21, 21, 22, 23, 27, &
      31, 33, 36, 39, 39, 40, 40, 41, 42, 42, 42, 42, &
      43, 44, 45, 45, 17, 18, 20, 23, 26, 29, 30, 30, &
      32, 32, 32, 33, 33, 33, 34, 35, 35, 35, 36, 36, &
      36, 37, 37, 37, 39, 39, 42, 43, 43, 43, 44, 44, &
      28, 30, 30, 37, 38, 38, 38, 42, 43, 44, 44, 21, &
      23, 24, 31, 39, 40, 40, 41, 42, 42, 43, 44, 45, &
      45, 32, 33, 35, 36, 37, 38, 42, 45 /)

  INTEGER, PARAMETER, DIMENSION(345) :: IHESS_K = (/ &
      41, 41, 41, 41, 41, 41, 41, 41, 41, 41, 41, 41, &
      41, 41, 41, 41, 41, 41, 41, 41, 41, 41, 41, 41, &
      41, 39, 40, 41, 46, 39, 40, 41, 46, 41, 41, 39, &
      40, 41, 46, 39, 40, 41, 46, 41, 41, 39, 40, 41, &
      46, 39, 40, 41, 46, 41, 39, 40, 41, 46, 44, 42, &
      40, 41, 41, 43, 41, 42, 41, 21, 41, 45, 45, 41, &
      43, 41, 41, 45, 45, 39, 40, 41, 41, 40, 46, 40, &
      41, 41, 40, 40, 40, 42, 41, 41, 40, 41, 42, 41, &
      41, 41, 39, 39, 41, 41, 39, 41, 39, 46, 39, 40, &
      41, 46, 39, 46, 39, 46, 45, 41, 39, 41, 41, 41, &
      41, 41, 39, 41, 31, 45, 41, 46, 39, 40, 41, 41, &
      39, 40, 41, 46, 41, 46, 44, 44, 45, 39, 41, 46, &
      39, 41, 46, 39, 40, 41, 46, 41, 39, 41, 46, 41, &
      39, 40, 41, 46, 39, 46, 39, 41, 39, 41, 46, 39, &
      41, 40, 41, 46, 39, 40, 41, 46, 39, 41, 44, 44, &
      45, 39, 40, 41, 46, 46, 39, 40, 41, 46, 39, 41, &
      39, 40, 41, 46, 41, 39, 40, 41, 46, 39, 41, 46, &
      40, 41, 46, 39, 39, 39, 39, 39, 41, 42, 43, 45, &
      41, 40, 40, 40, 40, 40, 40, 42, 42, 45, 46, 41, &
      41, 41, 41, 41, 41, 41, 41, 41, 39, 41, 41, 46, &
      39, 41, 41, 41, 46, 39, 41, 46, 39, 41, 41, 46, &
      41, 43, 42, 45, 44, 45, 42, 21, 41, 41, 45, 42, &
      45, 40, 40, 42, 45, 42, 45, 42, 43, 44, 45, 46, &
      45, 45, 45, 46, 41, 41, 41, 45, 41, 41, 39, 41, &
      39, 41, 46, 39, 40, 41, 41, 40, 41, 46, 39, 41, &
      46, 39, 41, 46, 41, 43, 43, 43, 44, 45, 44, 45, &
      41, 39, 41, 41, 40, 41, 46, 44, 44, 44, 45, 21, &
      45, 45, 45, 45, 42, 45, 45, 45, 46, 45, 45, 45, &
      46, 46, 46, 46, 46, 46, 46, 46, 46 /)


END MODULE cbm4_SOA_HessianSP

