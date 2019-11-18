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
! File                 : cbm42_strato_SOA_HessianSP.f90
! Time                 : Fri Sep 21 17:34:02 2018
! Working directory    : /home/hanniner/tmp/silam_v5_6-test/kpp/cbm42_strato_SOA
! Equation file        : cbm42_strato_SOA.kpp
! Output root filename : cbm42_strato_SOA
! 
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



MODULE cbm42_strato_SOA_HessianSP

  PUBLIC
  SAVE


! Hessian Sparse Data
! 

  INTEGER, PARAMETER, DIMENSION(360) :: IHESS_I_0 = (/ &
       6,  7,  7,  7,  7,  8,  8,  8,  8,  9,  9,  9, &
       9, 10, 10, 10, 10, 11, 11, 12, 12, 13, 14, 15, &
      15, 15, 15, 15, 15, 15, 15, 15, 15, 16, 16, 16, &
      16, 16, 16, 16, 16, 16, 16, 17, 17, 17, 17, 17, &
      17, 17, 17, 17, 17, 18, 18, 18, 18, 18, 19, 20, &
      20, 20, 21, 22, 23, 24, 25, 26, 27, 27, 28, 29, &
      30, 30, 31, 31, 32, 32, 32, 32, 33, 33, 33, 33, &
      33, 33, 33, 33, 33, 34, 34, 34, 35, 35, 36, 36, &
      37, 37, 38, 38, 39, 39, 39, 39, 39, 40, 40, 40, &
      40, 40, 41, 41, 41, 41, 42, 42, 42, 42, 42, 42, &
      43, 43, 43, 43, 43, 43, 43, 43, 43, 43, 43, 43, &
      43, 43, 43, 44, 44, 44, 44, 44, 44, 45, 45, 45, &
      45, 45, 45, 46, 46, 46, 46, 47, 47, 47, 47, 47, &
      47, 47, 47, 47, 47, 47, 47, 47, 47, 47, 47, 47, &
      47, 47, 47, 47, 47, 47, 48, 48, 48, 48, 49, 49, &
      49, 49, 49, 49, 49, 49, 50, 50, 50, 50, 50, 50, &
      50, 50, 50, 50, 50, 51, 51, 51, 51, 52, 52, 52, &
      52, 52, 52, 53, 53, 53, 53, 53, 53, 53, 54, 54, &
      54, 54, 54, 54, 54, 54, 54, 54, 54, 55, 55, 55, &
      55, 56, 56, 56, 56, 56, 57, 57, 57, 57, 57, 57, &
      57, 57, 57, 57, 57, 57, 57, 57, 57, 57, 57, 58, &
      58, 58, 58, 58, 58, 58, 59, 59, 59, 59, 59, 59, &
      59, 59, 60, 60, 60, 60, 60, 60, 60, 60, 60, 60, &
      60, 60, 60, 60, 60, 60, 60, 60, 60, 60, 61, 61, &
      61, 61, 61, 61, 61, 61, 61, 61, 61, 61, 61, 62, &
      62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, &
      62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, 62, &
      62, 62, 62, 62, 63, 63, 63, 63, 63, 63, 63, 63, &
      63, 63, 63, 63, 63, 63, 63, 63, 63, 63, 63, 63, &
      63, 63, 63, 63, 63, 63, 63, 63, 63, 63, 63, 63 /)
  INTEGER, PARAMETER, DIMENSION(227) :: IHESS_I_1 = (/ &
      63, 63, 63, 63, 63, 63, 63, 63, 63, 63, 63, 63, &
      63, 63, 63, 63, 63, 63, 63, 63, 63, 63, 63, 63, &
      63, 63, 63, 63, 63, 64, 64, 64, 64, 64, 64, 64, &
      64, 64, 64, 64, 64, 64, 64, 64, 65, 65, 65, 65, &
      65, 65, 65, 65, 65, 65, 65, 65, 65, 65, 65, 65, &
      65, 65, 65, 65, 65, 65, 65, 65, 65, 65, 65, 65, &
      65, 65, 65, 65, 65, 65, 65, 65, 65, 65, 65, 65, &
      65, 65, 65, 65, 65, 65, 66, 66, 66, 66, 66, 67, &
      67, 67, 67, 67, 67, 67, 67, 67, 67, 67, 67, 67, &
      67, 67, 67, 67, 67, 67, 67, 67, 67, 68, 68, 68, &
      68, 68, 68, 68, 68, 68, 68, 68, 68, 68, 68, 68, &
      68, 68, 68, 68, 68, 68, 69, 69, 69, 69, 69, 69, &
      69, 69, 69, 69, 69, 70, 70, 70, 70, 70, 70, 70, &
      70, 70, 70, 70, 70, 70, 70, 70, 70, 70, 70, 70, &
      70, 71, 71, 71, 71, 71, 71, 71, 71, 71, 71, 71, &
      71, 71, 71, 71, 71, 71, 71, 71, 71, 71, 71, 71, &
      71, 71, 71, 71, 71, 71, 71, 71, 71, 71, 71, 71, &
      71, 71, 72, 72, 72, 72, 72, 72, 72, 72, 72, 72, &
      72, 72, 72, 72, 72, 72, 72, 72, 72, 72, 72 /)
  INTEGER, PARAMETER, DIMENSION(587) :: IHESS_I = (/&
    IHESS_I_0, IHESS_I_1 /)

  INTEGER, PARAMETER, DIMENSION(360) :: IHESS_J_0 = (/ &
       7,  7,  8, 24, 29,  8,  9, 24, 29,  9, 10, 24, &
      29, 10, 11, 24, 29, 11, 12, 12, 13, 13, 15, 15, &
      16, 51, 51, 51, 51, 55, 55, 55, 55, 16, 17, 51, &
      51, 51, 51, 55, 55, 55, 55, 17, 18, 51, 51, 51, &
      51, 55, 55, 55, 55, 18, 51, 51, 51, 51, 19, 20, &
      41, 41, 69, 68, 54, 24, 25, 26, 27, 68, 60, 29, &
      30, 30, 30, 31, 32, 32, 63, 70, 33, 50, 51, 51, &
      51, 55, 55, 56, 56, 24, 29, 34, 35, 67, 36, 65, &
      37, 65, 38, 50, 29, 39, 46, 55, 55, 40, 40, 40, &
      63, 65, 24, 29, 41, 41, 41, 42, 57, 60, 60, 63, &
      43, 45, 45, 46, 46, 51, 55, 55, 56, 56, 60, 61, &
      62, 63, 67, 27, 44, 44, 48, 49, 68, 45, 45, 45, &
      55, 55, 55, 34, 41, 46, 46, 24, 29, 39, 41, 45, &
      45, 46, 46, 47, 47, 50, 51, 51, 51, 54, 54, 54, &
      55, 55, 56, 56, 56, 56, 48, 48, 48, 68, 35, 44, &
      48, 49, 49, 49, 59, 65, 29, 50, 51, 51, 51, 55, &
      55, 56, 56, 56, 56, 51, 51, 51, 51, 52, 52, 52, &
      61, 61, 63, 40, 53, 53, 53, 63, 65, 67, 39, 46, &
      46, 54, 54, 54, 54, 55, 57, 57, 57, 55, 55, 55, &
      55, 55, 56, 56, 56, 56, 45, 46, 50, 51, 51, 51, &
      51, 55, 55, 55, 56, 56, 56, 56, 57, 57, 57, 43, &
      52, 53, 58, 58, 58, 62, 59, 59, 59, 59, 59, 59, &
      68, 68, 41, 42, 48, 48, 48, 51, 55, 56, 57, 60, &
      60, 60, 60, 60, 60, 60, 60, 60, 62, 64, 26, 52, &
      52, 52, 59, 61, 61, 61, 62, 63, 68, 69, 69, 30, &
      30, 31, 36, 40, 45, 48, 49, 51, 52, 52, 53, 53, &
      55, 56, 57, 58, 59, 60, 62, 62, 62, 62, 62, 62, &
      62, 62, 63, 64, 19, 24, 25, 26, 29, 32, 35, 36, &
      37, 39, 40, 40, 41, 42, 43, 44, 45, 45, 46, 46, &
      48, 49, 49, 49, 50, 51, 51, 52, 52, 52, 53, 53 /)
  INTEGER, PARAMETER, DIMENSION(227) :: IHESS_J_1 = (/ &
      53, 54, 55, 55, 56, 56, 56, 57, 57, 58, 58, 58, &
      59, 60, 60, 62, 62, 62, 63, 63, 63, 63, 63, 63, &
      63, 63, 64, 65, 65, 45, 46, 51, 55, 56, 58, 59, &
      61, 62, 63, 64, 64, 64, 64, 64, 24, 29, 34, 40, &
      40, 40, 41, 45, 45, 45, 46, 46, 50, 51, 51, 51, &
      54, 54, 54, 55, 55, 55, 56, 56, 56, 58, 60, 60, &
      60, 61, 61, 62, 62, 63, 63, 63, 63, 63, 64, 65, &
      65, 65, 65, 65, 65, 67, 31, 44, 52, 53, 64, 19, &
      25, 40, 44, 44, 48, 49, 53, 53, 53, 59, 60, 60, &
      62, 63, 64, 65, 67, 67, 68, 68, 68, 44, 48, 49, &
      49, 53, 59, 59, 59, 59, 59, 60, 60, 62, 63, 64, &
      65, 65, 68, 68, 68, 68, 36, 52, 59, 61, 62, 63, &
      65, 68, 69, 69, 69, 30, 31, 32, 33, 34, 47, 54, &
      58, 59, 60, 60, 62, 62, 63, 64, 65, 68, 69, 70, &
      70, 20, 30, 32, 32, 34, 35, 37, 38, 47, 51, 54, &
      54, 56, 58, 59, 60, 60, 60, 60, 60, 60, 60, 60, &
      62, 62, 63, 64, 64, 65, 65, 67, 68, 68, 69, 69, &
      70, 70, 45, 45, 45, 46, 46, 51, 51, 54, 54, 54, &
      55, 55, 56, 56, 56, 56, 60, 61, 62, 63, 67 /)
  INTEGER, PARAMETER, DIMENSION(587) :: IHESS_J = (/&
    IHESS_J_0, IHESS_J_1 /)

  INTEGER, PARAMETER, DIMENSION(360) :: IHESS_K_0 = (/ &
      63, 63, 63, 63, 63, 63, 63, 63, 63, 63, 63, 63, &
      63, 63, 63, 63, 63, 63, 63, 63, 63, 63, 63, 63, &
      63, 60, 62, 63, 64, 60, 62, 63, 64, 63, 63, 60, &
      62, 63, 64, 60, 62, 63, 64, 63, 63, 60, 62, 63, &
      64, 60, 62, 63, 64, 63, 60, 62, 63, 64, 63, 71, &
      60, 63, 71, 69, 71, 63, 63, 63, 67, 68, 71, 63, &
      70, 71, 71, 66, 32, 63, 70, 71, 70, 63, 60, 63, &
      64, 60, 63, 60, 62, 63, 63, 70, 63, 71, 62, 69, &
      63, 71, 71, 63, 63, 63, 64, 63, 64, 62, 63, 67, &
      63, 65, 63, 63, 60, 63, 60, 63, 60, 65, 72, 71, &
      63, 62, 64, 63, 64, 64, 62, 64, 62, 64, 72, 72, &
      72, 72, 72, 67, 63, 66, 67, 67, 68, 62, 63, 64, &
      62, 63, 64, 70, 63, 63, 64, 63, 63, 63, 63, 62, &
      63, 63, 64, 47, 70, 63, 60, 63, 64, 54, 65, 70, &
      62, 63, 60, 62, 63, 64, 62, 63, 67, 71, 63, 63, &
      63, 62, 63, 67, 63, 68, 63, 63, 62, 63, 64, 62, &
      64, 60, 62, 63, 64, 60, 62, 63, 64, 62, 63, 66, &
      65, 72, 69, 67, 62, 63, 66, 68, 67, 72, 63, 63, &
      64, 54, 65, 70, 71, 63, 60, 62, 63, 60, 62, 63, &
      64, 62, 60, 62, 63, 64, 63, 64, 63, 60, 62, 63, &
      64, 62, 63, 64, 60, 62, 63, 64, 60, 62, 63, 63, &
      66, 66, 64, 65, 71, 63, 61, 62, 63, 64, 67, 70, &
      68, 69, 60, 63, 62, 63, 67, 60, 60, 60, 60, 60, &
      62, 63, 65, 67, 68, 70, 71, 72, 71, 71, 63, 62, &
      63, 66, 61, 64, 65, 72, 69, 69, 69, 69, 70, 70, &
      71, 66, 62, 62, 62, 62, 62, 62, 62, 66, 62, 66, &
      62, 62, 62, 65, 62, 62, 63, 64, 65, 68, 69, 70, &
      71, 72, 63, 66, 63, 63, 63, 63, 63, 63, 63, 62, &
      63, 63, 62, 63, 63, 63, 63, 63, 62, 63, 63, 64, &
      63, 62, 63, 67, 63, 63, 64, 62, 63, 66, 62, 63 /)
  INTEGER, PARAMETER, DIMENSION(227) :: IHESS_K_1 = (/ &
      66, 65, 63, 64, 62, 63, 64, 62, 63, 64, 65, 71, &
      63, 63, 65, 63, 65, 72, 63, 64, 65, 68, 69, 70, &
      71, 72, 65, 67, 70, 64, 64, 64, 64, 64, 64, 64, &
      64, 64, 64, 65, 66, 67, 70, 71, 63, 63, 70, 62, &
      63, 67, 63, 62, 63, 64, 63, 64, 63, 60, 63, 64, &
      54, 65, 70, 62, 63, 64, 62, 63, 64, 65, 63, 65, &
      72, 65, 72, 65, 72, 64, 65, 68, 69, 72, 65, 65, &
      67, 68, 69, 70, 71, 72, 66, 66, 66, 66, 66, 63, &
      63, 67, 63, 66, 67, 67, 62, 63, 66, 67, 67, 68, &
      68, 68, 67, 67, 71, 72, 68, 69, 70, 66, 62, 62, &
      63, 66, 61, 62, 64, 67, 70, 67, 68, 68, 68, 67, &
      67, 68, 68, 69, 70, 71, 62, 66, 61, 64, 69, 69, &
      69, 69, 69, 70, 71, 70, 66, 32, 70, 70, 70, 70, &
      71, 70, 70, 71, 70, 71, 70, 70, 70, 70, 70, 70, &
      71, 71, 71, 32, 63, 70, 63, 63, 71, 70, 60, 70, &
      71, 60, 71, 70, 60, 62, 63, 65, 67, 68, 70, 71, &
      70, 71, 71, 70, 71, 70, 71, 71, 70, 71, 70, 71, &
      70, 71, 62, 63, 64, 63, 64, 63, 64, 54, 65, 70, &
      63, 64, 60, 62, 63, 64, 72, 72, 72, 72, 72 /)
  INTEGER, PARAMETER, DIMENSION(587) :: IHESS_K = (/&
    IHESS_K_0, IHESS_K_1 /)


END MODULE cbm42_strato_SOA_HessianSP

