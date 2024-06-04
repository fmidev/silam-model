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
! File                 : cbm5_SOA_HessianSP.f90
! Time                 : Thu Apr  4 16:17:23 2024
! Working directory    : /mnt/d/kpp/kpp/cbm5_SOA
! Equation file        : cbm5_SOA.kpp
! Output root filename : cbm5_SOA
! 
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



MODULE cbm5_SOA_HessianSP

  PUBLIC
  SAVE


! Hessian Sparse Data
! 

  INTEGER, PARAMETER, DIMENSION(360) :: IHESS_I_0 = (/ &
       1,  2,  3,  4,  4,  5,  5,  7,  8,  8,  8,  8, &
       8,  9,  9,  9,  9,  9, 10, 10, 10, 10, 10, 11, &
      11, 11, 11, 11, 12, 13, 14, 14, 15, 16, 17, 18, &
      19, 19, 19, 19, 19, 19, 19, 19, 19, 19, 19, 19, &
      19, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, &
      20, 20, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, &
      21, 21, 21, 22, 22, 22, 23, 23, 23, 23, 23, 23, &
      23, 23, 24, 24, 24, 24, 25, 25, 25, 25, 26, 26, &
      26, 27, 27, 27, 27, 27, 27, 27, 28, 28, 28, 29, &
      29, 30, 30, 30, 31, 31, 31, 32, 32, 32, 33, 33, &
      33, 33, 33, 34, 34, 34, 34, 35, 35, 36, 36, 37, &
      37, 37, 38, 38, 38, 38, 38, 38, 38, 38, 38, 38, &
      38, 38, 38, 38, 38, 38, 38, 39, 39, 39, 39, 40, &
      40, 40, 40, 40, 40, 40, 40, 40, 41, 41, 41, 41, &
      42, 42, 42, 42, 43, 43, 43, 43, 44, 44, 44, 44, &
      45, 45, 45, 45, 45, 46, 46, 46, 46, 46, 47, 47, &
      47, 47, 47, 47, 47, 47, 48, 48, 48, 48, 48, 48, &
      48, 48, 48, 49, 49, 49, 49, 49, 49, 49, 49, 49, &
      49, 49, 49, 49, 49, 49, 49, 49, 49, 49, 49, 49, &
      49, 49, 49, 49, 49, 50, 50, 50, 50, 50, 50, 50, &
      50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, &
      50, 50, 50, 50, 50, 51, 51, 51, 51, 51, 51, 51, &
      51, 51, 51, 51, 51, 51, 51, 51, 51, 51, 51, 51, &
      51, 51, 51, 51, 51, 51, 51, 51, 51, 51, 51, 52, &
      52, 52, 52, 52, 52, 52, 52, 52, 52, 52, 52, 52, &
      52, 53, 53, 53, 53, 53, 53, 53, 53, 53, 53, 53, &
      53, 54, 54, 54, 54, 54, 54, 54, 54, 54, 54, 54, &
      54, 54, 54, 54, 54, 54, 54, 54, 54, 54, 54, 54, &
      54, 54, 54, 54, 54, 54, 54, 54, 54, 54, 54, 54, &
      54, 54, 54, 54, 54, 54, 54, 54, 54, 54, 54, 54 /)
  INTEGER, PARAMETER, DIMENSION(287) :: IHESS_I_1 = (/ &
      54, 54, 54, 54, 54, 54, 54, 54, 54, 54, 54, 54, &
      54, 54, 54, 54, 54, 54, 54, 54, 55, 55, 55, 55, &
      55, 55, 55, 55, 55, 55, 55, 55, 56, 56, 56, 56, &
      56, 56, 56, 56, 56, 56, 56, 56, 56, 56, 56, 56, &
      57, 57, 57, 57, 57, 57, 57, 57, 57, 57, 57, 57, &
      57, 57, 57, 57, 58, 58, 58, 58, 58, 58, 58, 58, &
      58, 58, 58, 58, 58, 58, 58, 58, 58, 58, 58, 58, &
      59, 59, 59, 59, 59, 59, 59, 59, 59, 59, 59, 59, &
      59, 59, 59, 59, 59, 59, 59, 59, 59, 59, 59, 59, &
      59, 59, 59, 59, 59, 59, 59, 59, 59, 59, 59, 59, &
      59, 59, 59, 59, 59, 59, 59, 59, 59, 59, 59, 59, &
      59, 59, 59, 59, 59, 59, 59, 59, 59, 59, 59, 59, &
      59, 59, 59, 59, 59, 59, 60, 60, 60, 60, 60, 60, &
      60, 60, 60, 60, 60, 60, 60, 60, 60, 60, 60, 60, &
      60, 60, 60, 60, 60, 60, 60, 60, 60, 60, 60, 60, &
      60, 60, 60, 60, 60, 60, 60, 60, 60, 60, 60, 60, &
      61, 61, 61, 61, 61, 61, 61, 61, 61, 61, 61, 61, &
      61, 61, 61, 61, 61, 61, 61, 61, 61, 61, 61, 61, &
      61, 61, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, &
      62, 62, 62, 62, 62, 63, 63, 63, 63, 63, 63, 63, &
      63, 63, 63, 63, 63, 63, 63, 63, 63, 63, 63, 64, &
      64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, &
      64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, &
      64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64 /)
  INTEGER, PARAMETER, DIMENSION(647) :: IHESS_I = (/&
    IHESS_I_0, IHESS_I_1 /)

  INTEGER, PARAMETER, DIMENSION(360) :: IHESS_J_0 = (/ &
       8, 19,  3,  3,  4,  4,  5,  7,  7,  8,  9, 12, &
      17,  7,  9, 10, 12, 17,  7, 10, 11, 12, 17,  5, &
       7, 11, 12, 17, 12, 13, 14, 53, 15, 58, 17, 57, &
      19, 20, 22, 22, 22, 44, 44, 44, 44, 46, 46, 46, &
      46, 20, 21, 22, 22, 22, 44, 44, 44, 44, 46, 46, &
      46, 46, 21, 22, 22, 22, 23, 44, 44, 44, 44, 46, &
      46, 46, 46, 22, 22, 22, 22, 22, 22, 23, 44, 44, &
      44, 44, 24, 24, 54, 59, 25, 25, 54, 63, 26, 31, &
      41, 27, 53, 53, 57, 57, 59, 60, 28, 57, 59, 29, &
      59, 12, 17, 30, 31, 31, 51, 32, 55, 59, 17, 33, &
      39, 47, 47, 34, 34, 45, 45, 35, 62, 36, 54, 31, &
      37, 53, 38, 39, 39, 41, 41, 42, 42, 43, 43, 44, &
      46, 47, 47, 47, 51, 51, 51, 30, 39, 39, 45, 40, &
      45, 47, 48, 49, 50, 51, 54, 58, 41, 41, 41, 41, &
      42, 42, 42, 42, 43, 43, 43, 43, 44, 44, 44, 44, &
      12, 17, 34, 45, 45, 46, 46, 46, 46, 46, 46, 46, &
      46, 46, 46, 47, 47, 47, 30, 34, 36, 44, 46, 46, &
      47, 48, 55, 15, 32, 39, 41, 42, 42, 42, 42, 43, &
      43, 43, 43, 44, 44, 44, 44, 46, 46, 46, 47, 47, &
      48, 49, 49, 49, 54, 13, 15, 32, 35, 42, 42, 42, &
      42, 43, 43, 43, 43, 47, 47, 48, 50, 50, 50, 53, &
      54, 57, 60, 62, 62, 14, 15, 39, 39, 41, 41, 41, &
      41, 42, 43, 43, 43, 43, 44, 44, 46, 46, 46, 47, &
      47, 47, 48, 51, 51, 51, 51, 53, 53, 53, 53, 39, &
      41, 42, 43, 44, 46, 47, 52, 52, 52, 52, 52, 57, &
      59, 27, 37, 53, 53, 53, 53, 53, 57, 57, 57, 57, &
      57,  3,  4,  5,  8,  9, 10, 11, 12, 13, 14, 15, &
      17, 19, 20, 21, 23, 24, 24, 25, 26, 27, 28, 29, &
      32, 33, 35, 37, 38, 39, 39, 40, 41, 41, 41, 42, &
      42, 43, 43, 43, 44, 44, 45, 46, 46, 47, 47, 48 /)
  INTEGER, PARAMETER, DIMENSION(287) :: IHESS_J_1 = (/ &
      49, 49, 50, 50, 51, 51, 52, 52, 54, 54, 54, 54, &
      54, 54, 54, 56, 57, 58, 59, 59, 13, 43, 43, 44, &
      44, 44, 46, 54, 55, 55, 55, 55, 24, 41, 42, 42, &
      43, 44, 46, 49, 50, 51, 54, 54, 56, 56, 56, 56, &
      28, 33, 39, 39, 47, 47, 50, 50, 50, 53, 57, 57, &
      57, 57, 57, 57, 40, 41, 42, 43, 44, 45, 46, 47, &
      49, 50, 51, 52, 52, 54, 56, 56, 58, 58, 58, 58, &
      12, 13, 14, 15, 17, 24, 24, 26, 30, 31, 31, 34, &
      37, 38, 39, 39, 41, 41, 41, 42, 42, 42, 42, 43, &
      43, 43, 44, 44, 44, 45, 46, 46, 46, 46, 46, 47, &
      47, 47, 48, 51, 51, 51, 51, 52, 52, 53, 53, 53, &
      53, 53, 54, 54, 54, 54, 55, 56, 57, 57, 58, 59, &
      59, 59, 59, 59, 62, 62, 12, 13, 15, 17, 32, 33, &
      37, 39, 39, 41, 41, 41, 42, 42, 43, 43, 43, 43, &
      44, 44, 44, 45, 46, 46, 46, 46, 46, 47, 47, 47, &
      53, 54, 55, 57, 57, 59, 59, 60, 60, 60, 62, 62, &
      17, 30, 36, 42, 42, 42, 42, 43, 43, 43, 44, 44, &
      44, 44, 45, 46, 46, 46, 46, 47, 47, 47, 54, 55, &
      55, 55, 44, 46, 46, 47, 47, 49, 49, 49, 53, 57, &
      59, 60, 62, 62, 62, 25, 30, 31, 46, 52, 53, 54, &
      55, 56, 56, 57, 58, 58, 59, 60, 62, 63, 63, 25, &
      25, 29, 30, 31, 34, 35, 36, 41, 42, 43, 44, 46, &
      46, 52, 52, 52, 53, 54, 54, 56, 56, 56, 57, 57, &
      58, 58, 58, 58, 59, 59, 60, 62, 62, 63, 63 /)
  INTEGER, PARAMETER, DIMENSION(647) :: IHESS_J = (/&
    IHESS_J_0, IHESS_J_1 /)

  INTEGER, PARAMETER, DIMENSION(360) :: IHESS_K_0 = (/ &
      54, 54, 54, 54, 54, 54, 54, 54, 54, 54, 54, 54, &
      54, 54, 54, 54, 54, 54, 54, 54, 54, 54, 54, 54, &
      54, 54, 54, 54, 54, 54, 54, 53, 54, 64, 54, 64, &
      54, 54, 52, 54, 58, 52, 54, 56, 58, 52, 54, 56, &
      58, 54, 54, 52, 54, 58, 52, 54, 56, 58, 52, 54, &
      56, 58, 54, 52, 54, 58, 54, 52, 54, 56, 58, 52, &
      54, 56, 58, 52, 54, 58, 52, 54, 58, 54, 52, 54, &
      56, 58, 54, 56, 54, 59, 25, 54, 63, 64, 54, 63, &
      52, 54, 57, 62, 59, 60, 62, 62, 54, 59, 62, 54, &
      64, 54, 54, 63, 59, 63, 59, 54, 59, 60, 54, 54, &
      52, 52, 54, 59, 64, 54, 58, 54, 64, 64, 61, 59, &
      54, 59, 54, 52, 54, 52, 56, 52, 56, 52, 56, 52, &
      52, 52, 54, 58, 54, 56, 58, 63, 52, 54, 54, 54, &
      58, 58, 54, 58, 58, 58, 64, 59, 52, 54, 56, 58, &
      52, 54, 56, 58, 52, 54, 56, 58, 52, 54, 56, 58, &
      54, 54, 59, 54, 58, 52, 54, 56, 58, 64, 52, 54, &
      56, 58, 64, 52, 54, 58, 63, 64, 64, 58, 58, 64, &
      58, 54, 63, 54, 54, 52, 54, 52, 54, 56, 58, 52, &
      54, 56, 58, 52, 54, 56, 58, 52, 58, 64, 54, 58, &
      54, 54, 56, 58, 61, 54, 54, 54, 54, 52, 54, 56, &
      58, 52, 54, 56, 58, 52, 54, 54, 54, 56, 58, 62, &
      61, 62, 62, 62, 63, 54, 54, 52, 54, 52, 54, 56, &
      58, 52, 52, 54, 56, 58, 52, 54, 52, 54, 56, 52, &
      54, 58, 54, 54, 56, 58, 59, 53, 57, 62, 63, 52, &
      52, 52, 52, 52, 52, 52, 54, 58, 59, 63, 64, 59, &
      62, 54, 54, 53, 57, 59, 62, 63, 57, 59, 60, 62, &
      63, 54, 54, 54, 54, 54, 54, 54, 54, 54, 54, 54, &
      54, 54, 54, 54, 54, 54, 56, 54, 54, 54, 54, 54, &
      54, 54, 54, 54, 54, 52, 54, 54, 52, 54, 56, 52, &
      54, 52, 54, 56, 52, 54, 54, 52, 54, 52, 54, 54 /)
  INTEGER, PARAMETER, DIMENSION(287) :: IHESS_K_1 = (/ &
      54, 56, 54, 56, 54, 56, 54, 59, 54, 56, 58, 59, &
      61, 63, 64, 59, 59, 59, 62, 63, 54, 56, 58, 52, &
      54, 58, 54, 61, 55, 59, 60, 63, 56, 56, 52, 56, &
      56, 56, 56, 56, 56, 56, 54, 56, 58, 59, 63, 64, &
      54, 54, 52, 54, 52, 54, 54, 56, 58, 57, 57, 59, &
      60, 62, 63, 64, 54, 58, 58, 58, 58, 58, 58, 58, &
      58, 58, 58, 58, 64, 58, 58, 64, 58, 59, 63, 64, &
      54, 54, 54, 54, 54, 54, 56, 54, 63, 59, 63, 59, &
      54, 54, 52, 54, 52, 54, 56, 52, 54, 56, 58, 52, &
      54, 56, 52, 54, 58, 54, 52, 54, 56, 58, 64, 52, &
      54, 58, 54, 54, 56, 58, 59, 54, 59, 53, 57, 59, &
      62, 63, 56, 58, 59, 61, 59, 59, 59, 62, 59, 59, &
      60, 62, 63, 64, 62, 63, 54, 54, 54, 54, 54, 54, &
      54, 52, 54, 54, 56, 58, 54, 56, 52, 54, 56, 58, &
      52, 54, 58, 54, 52, 54, 56, 58, 64, 52, 54, 58, &
      62, 61, 60, 60, 62, 60, 62, 60, 62, 63, 62, 63, &
      54, 63, 64, 52, 54, 56, 58, 52, 54, 56, 52, 54, &
      56, 58, 54, 52, 56, 58, 64, 52, 54, 58, 61, 55, &
      59, 60, 52, 52, 56, 54, 58, 54, 56, 58, 62, 62, &
      62, 62, 62, 63, 64, 25, 63, 63, 64, 63, 63, 63, &
      63, 63, 64, 63, 63, 64, 63, 63, 63, 63, 64, 25, &
      54, 54, 63, 63, 64, 54, 64, 58, 58, 58, 58, 58, &
      64, 58, 63, 64, 63, 58, 64, 58, 63, 64, 63, 64, &
      58, 59, 63, 64, 63, 64, 63, 63, 64, 63, 64 /)
  INTEGER, PARAMETER, DIMENSION(647) :: IHESS_K = (/&
    IHESS_K_0, IHESS_K_1 /)


END MODULE cbm5_SOA_HessianSP

