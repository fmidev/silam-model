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
! Time                 : Wed Nov  6 16:46:40 2019
! Working directory    : /home/hanniner/Silam/silam_v5_7-risto/kpp/cbm42_strato_SOA
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
       5,  5,  5,  5,  5,  5,  5,  6,  7,  7,  7,  7, &
       8,  8,  8,  8,  9,  9,  9,  9, 10, 10, 10, 10, &
      11, 11, 12, 12, 13, 14, 15, 15, 15, 15, 15, 15, &
      15, 15, 15, 15, 16, 16, 16, 16, 16, 16, 16, 16, &
      16, 16, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, &
      18, 18, 18, 18, 18, 19, 20, 21, 22, 23, 23, 23, &
      24, 25, 26, 27, 28, 29, 29, 30, 31, 32, 32, 33, &
      33, 34, 34, 34, 34, 35, 35, 35, 36, 36, 37, 37, &
      38, 38, 39, 39, 40, 40, 41, 41, 41, 41, 41, 42, &
      42, 42, 42, 42, 43, 43, 43, 43, 44, 44, 44, 45, &
      45, 45, 45, 45, 45, 45, 46, 46, 46, 46, 46, 46, &
      46, 46, 46, 46, 46, 46, 46, 46, 46, 46, 46, 46, &
      46, 47, 47, 47, 47, 48, 48, 48, 48, 48, 48, 49, &
      49, 49, 49, 49, 49, 49, 49, 49, 49, 49, 50, 50, &
      50, 50, 51, 51, 51, 51, 51, 51, 51, 51, 52, 52, &
      52, 52, 52, 52, 53, 53, 53, 53, 53, 53, 53, 53, &
      54, 54, 54, 54, 55, 55, 55, 55, 56, 56, 56, 56, &
      56, 56, 56, 56, 56, 56, 56, 56, 56, 56, 56, 56, &
      56, 56, 56, 56, 56, 56, 56, 56, 56, 56, 56, 56, &
      56, 56, 56, 56, 56, 56, 57, 57, 57, 57, 57, 57, &
      57, 57, 57, 57, 57, 57, 57, 57, 58, 58, 58, 58, &
      58, 58, 58, 58, 58, 58, 58, 58, 58, 58, 58, 59, &
      59, 59, 59, 59, 59, 59, 59, 60, 60, 60, 60, 60, &
      61, 61, 61, 61, 61, 61, 61, 61, 61, 61, 61, 61, &
      61, 61, 61, 61, 61, 61, 61, 61, 61, 62, 62, 62, &
      62, 62, 62, 63, 63, 63, 63, 63, 63, 63, 63, 64, &
      64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, &
      64, 64, 64, 64, 64, 64, 64, 64, 65, 65, 65, 65, &
      65, 65, 65, 65, 65, 65, 65, 65, 65, 66, 66, 66, &
      66, 66, 66, 66, 66, 66, 66, 66, 66, 66, 66, 66 /)
  INTEGER, PARAMETER, DIMENSION(299) :: IHESS_I_1 = (/ &
      66, 66, 66, 66, 66, 66, 66, 66, 66, 66, 66, 66, &
      66, 66, 66, 66, 66, 66, 66, 66, 66, 66, 66, 66, &
      66, 66, 66, 66, 66, 66, 66, 66, 66, 66, 66, 66, &
      66, 66, 66, 66, 66, 66, 66, 66, 66, 66, 66, 66, &
      66, 66, 67, 67, 67, 67, 67, 67, 67, 67, 67, 67, &
      67, 67, 67, 67, 67, 67, 68, 68, 68, 68, 68, 68, &
      68, 68, 68, 68, 68, 68, 68, 68, 68, 68, 68, 68, &
      68, 68, 68, 68, 68, 68, 68, 68, 68, 68, 68, 68, &
      68, 68, 68, 68, 68, 68, 68, 68, 68, 68, 68, 68, &
      68, 68, 68, 68, 68, 68, 68, 68, 68, 68, 68, 68, &
      68, 68, 68, 68, 69, 69, 69, 69, 69, 69, 69, 69, &
      69, 69, 69, 69, 69, 69, 69, 69, 69, 69, 69, 69, &
      69, 69, 70, 70, 70, 70, 70, 70, 70, 70, 70, 70, &
      70, 70, 70, 70, 70, 70, 70, 70, 70, 70, 70, 70, &
      70, 70, 70, 70, 70, 70, 70, 71, 71, 71, 71, 71, &
      72, 72, 72, 72, 72, 72, 72, 72, 72, 72, 72, 72, &
      72, 72, 72, 72, 72, 72, 72, 72, 72, 73, 73, 73, &
      73, 73, 73, 73, 73, 73, 73, 73, 73, 73, 73, 73, &
      73, 73, 73, 73, 73, 73, 73, 73, 73, 73, 73, 73, &
      73, 73, 73, 73, 73, 73, 73, 73, 73, 73, 73, 73, &
      74, 74, 74, 74, 74, 74, 74, 74, 74, 74, 74, 74, &
      74, 74, 74, 74, 74, 74, 74, 74, 74, 75, 75, 75, &
      75, 75, 75, 75, 75, 75, 75, 75, 76, 76, 76, 76, &
      76, 76, 76, 76, 76, 76, 76, 76, 76, 76, 76, 76, &
      76, 76, 76, 76, 76, 76, 76, 76, 76, 76, 76 /)
  INTEGER, PARAMETER, DIMENSION(659) :: IHESS_I = (/&
    IHESS_I_0, IHESS_I_1 /)

  INTEGER, PARAMETER, DIMENSION(360) :: IHESS_J_0 = (/ &
      23, 35, 40, 49, 59, 60, 60,  7,  7,  8, 28, 31, &
       8,  9, 28, 31,  9, 10, 28, 31, 10, 11, 28, 31, &
      11, 12, 12, 13, 13, 15, 15, 16, 55, 55, 55, 55, &
      60, 60, 60, 60, 16, 17, 55, 55, 55, 55, 60, 60, &
      60, 60, 17, 18, 55, 55, 55, 55, 60, 60, 60, 60, &
      18, 55, 55, 55, 55, 75, 20, 21, 22, 23, 43, 43, &
      73, 74, 58, 27, 28, 29, 74, 64, 31, 32, 32, 32, &
      33, 34, 34, 66, 72, 28, 31, 35, 36, 69, 37, 68, &
      38, 68, 39, 39, 40, 57, 31, 41, 47, 59, 59, 42, &
      42, 42, 66, 68, 28, 31, 43, 43, 44, 44, 44, 43, &
      45, 59, 61, 64, 64, 66, 39, 39, 44, 44, 46, 47, &
      47, 54, 54, 55, 59, 59, 59, 60, 64, 65, 66, 69, &
      70, 35, 43, 47, 47, 29, 48, 48, 50, 51, 74, 49, &
      49, 49, 49, 54, 54, 55, 55, 55, 57, 60, 50, 50, &
      50, 73, 36, 48, 50, 51, 51, 51, 63, 68, 52, 52, &
      52, 65, 65, 66, 39, 42, 53, 53, 53, 66, 68, 69, &
      54, 54, 54, 54, 55, 55, 55, 55, 27, 28, 31, 39, &
      39, 41, 43, 44, 44, 47, 47, 49, 54, 54, 54, 54, &
      55, 55, 55, 56, 56, 56, 57, 58, 58, 58, 59, 59, &
      59, 60, 60, 60, 60, 60, 31, 40, 54, 55, 55, 55, &
      57, 59, 59, 59, 60, 60, 60, 60, 41, 47, 47, 58, &
      58, 58, 58, 59, 59, 59, 60, 60, 61, 61, 61, 59, &
      59, 59, 60, 60, 60, 60, 60, 60, 60, 60, 60, 60, &
      22, 44, 47, 54, 54, 54, 54, 55, 55, 55, 55, 57, &
      59, 59, 59, 60, 60, 60, 61, 61, 61, 52, 53, 62, &
      62, 62, 66, 63, 63, 63, 63, 63, 63, 74, 74, 43, &
      45, 50, 50, 50, 54, 55, 59, 60, 61, 64, 64, 64, &
      64, 64, 64, 64, 64, 64, 67, 70, 27, 52, 52, 52, &
      63, 65, 65, 65, 66, 70, 72, 74, 75, 20, 21, 22, &
      27, 28, 31, 34, 36, 37, 38, 39, 41, 42, 42, 43 /)
  INTEGER, PARAMETER, DIMENSION(299) :: IHESS_J_1 = (/ &
      44, 44, 45, 46, 47, 47, 48, 50, 51, 51, 51, 52, &
      52, 52, 53, 53, 53, 54, 54, 54, 55, 55, 57, 58, &
      59, 59, 60, 60, 61, 61, 62, 62, 62, 63, 64, 64, &
      66, 66, 66, 66, 66, 66, 66, 66, 66, 67, 68, 68, &
      68, 70, 44, 47, 54, 55, 59, 60, 62, 63, 65, 66, &
      67, 67, 67, 67, 67, 67, 21, 22, 28, 31, 35, 39, &
      39, 42, 42, 42, 43, 44, 44, 44, 46, 47, 47, 49, &
      54, 54, 54, 55, 55, 55, 56, 57, 58, 58, 58, 59, &
      59, 59, 60, 60, 60, 60, 60, 62, 64, 64, 64, 65, &
      65, 66, 66, 66, 66, 66, 67, 68, 68, 68, 68, 68, &
      68, 68, 69, 70, 20, 39, 42, 48, 48, 50, 51, 53, &
      53, 53, 63, 64, 64, 66, 67, 68, 69, 69, 70, 72, &
      74, 74, 32, 32, 33, 37, 42, 44, 50, 51, 52, 52, &
      53, 53, 54, 55, 60, 61, 62, 63, 64, 66, 66, 67, &
      67, 68, 70, 70, 70, 70, 70, 33, 48, 52, 53, 67, &
      32, 33, 34, 35, 49, 56, 58, 60, 62, 63, 64, 64, &
      66, 67, 68, 70, 70, 72, 72, 72, 72, 23, 32, 34, &
      34, 35, 36, 38, 40, 54, 55, 56, 58, 58, 60, 60, &
      62, 63, 64, 64, 64, 64, 64, 64, 64, 64, 66, 67, &
      67, 68, 68, 69, 70, 70, 72, 72, 72, 72, 73, 73, &
      48, 50, 51, 51, 53, 63, 63, 63, 63, 63, 64, 64, &
      66, 67, 68, 68, 70, 72, 73, 74, 74, 37, 52, 63, &
      65, 66, 68, 70, 72, 73, 74, 75, 21, 27, 44, 44, &
      44, 47, 47, 54, 54, 54, 54, 55, 55, 58, 58, 58, &
      59, 59, 59, 60, 60, 60, 64, 65, 66, 69, 70 /)
  INTEGER, PARAMETER, DIMENSION(659) :: IHESS_J = (/&
    IHESS_J_0, IHESS_J_1 /)

  INTEGER, PARAMETER, DIMENSION(360) :: IHESS_K_0 = (/ &
      73, 72, 73, 72, 64, 64, 73, 66, 66, 66, 66, 66, &
      66, 66, 66, 66, 66, 66, 66, 66, 66, 66, 66, 66, &
      66, 66, 66, 66, 66, 66, 66, 66, 64, 66, 67, 70, &
      64, 66, 67, 70, 66, 66, 64, 66, 67, 70, 64, 66, &
      67, 70, 66, 66, 64, 66, 67, 70, 64, 66, 67, 70, &
      66, 64, 66, 67, 70, 75, 66, 66, 66, 73, 64, 66, &
      75, 75, 73, 66, 66, 69, 74, 73, 66, 72, 73, 73, &
      71, 34, 66, 72, 73, 66, 66, 72, 66, 73, 70, 75, &
      66, 73, 66, 69, 73, 66, 66, 66, 67, 66, 67, 66, &
      69, 70, 66, 68, 66, 66, 64, 66, 66, 67, 70, 64, &
      66, 64, 64, 68, 76, 73, 66, 69, 67, 70, 66, 66, &
      67, 67, 70, 67, 64, 66, 67, 67, 76, 76, 76, 76, &
      76, 72, 66, 66, 67, 69, 66, 71, 69, 69, 74, 49, &
      56, 68, 72, 64, 70, 64, 66, 67, 66, 66, 66, 69, &
      70, 74, 66, 66, 66, 66, 69, 70, 66, 74, 66, 70, &
      71, 68, 76, 75, 69, 69, 66, 70, 71, 74, 69, 76, &
      64, 66, 67, 70, 64, 66, 67, 70, 66, 66, 66, 66, &
      69, 66, 66, 66, 70, 66, 67, 56, 64, 66, 67, 70, &
      64, 66, 67, 56, 68, 72, 66, 58, 68, 72, 64, 66, &
      67, 64, 66, 67, 70, 73, 66, 73, 70, 66, 67, 70, &
      66, 64, 66, 67, 64, 67, 70, 73, 66, 66, 67, 58, &
      68, 72, 73, 64, 66, 67, 67, 70, 64, 66, 70, 64, &
      66, 67, 64, 66, 67, 70, 73, 64, 66, 67, 70, 73, &
      66, 66, 67, 64, 66, 67, 70, 64, 66, 67, 70, 66, &
      64, 66, 67, 64, 67, 73, 64, 66, 70, 71, 71, 67, &
      68, 73, 70, 65, 66, 67, 69, 70, 72, 74, 75, 64, &
      66, 66, 69, 70, 64, 64, 64, 64, 64, 64, 66, 68, &
      69, 70, 72, 73, 74, 76, 73, 73, 66, 66, 70, 71, &
      65, 67, 68, 76, 75, 75, 75, 75, 75, 66, 66, 66, &
      66, 66, 66, 66, 66, 70, 66, 66, 66, 66, 70, 66 /)
  INTEGER, PARAMETER, DIMENSION(299) :: IHESS_K_1 = (/ &
      66, 70, 66, 66, 66, 67, 66, 66, 66, 69, 70, 66, &
      70, 71, 66, 70, 71, 66, 67, 70, 66, 67, 66, 68, &
      66, 67, 66, 67, 66, 70, 67, 68, 73, 66, 66, 68, &
      66, 67, 68, 70, 72, 73, 74, 75, 76, 68, 69, 70, &
      72, 76, 67, 67, 67, 67, 67, 67, 67, 67, 67, 67, &
      68, 69, 70, 71, 72, 73, 66, 66, 66, 66, 72, 66, &
      69, 66, 69, 70, 66, 66, 67, 70, 66, 66, 67, 68, &
      66, 67, 70, 64, 66, 67, 68, 66, 58, 68, 72, 64, &
      66, 67, 64, 66, 67, 70, 73, 68, 66, 68, 76, 68, &
      76, 67, 68, 74, 75, 76, 68, 68, 69, 70, 72, 73, &
      74, 75, 76, 76, 66, 66, 69, 66, 71, 69, 69, 66, &
      70, 71, 69, 69, 74, 74, 69, 69, 73, 76, 74, 74, &
      74, 75, 72, 73, 71, 70, 70, 70, 70, 70, 70, 71, &
      70, 71, 70, 70, 70, 70, 68, 70, 70, 66, 70, 70, &
      71, 70, 72, 73, 74, 75, 76, 71, 71, 71, 71, 71, &
      72, 71, 34, 72, 72, 72, 72, 73, 73, 72, 72, 73, &
      72, 72, 72, 72, 73, 72, 73, 74, 75, 73, 73, 34, &
      66, 72, 66, 66, 73, 64, 64, 72, 72, 73, 64, 73, &
      73, 72, 64, 66, 68, 69, 70, 72, 73, 74, 73, 72, &
      73, 72, 73, 73, 72, 73, 72, 73, 74, 75, 74, 75, &
      71, 70, 66, 70, 71, 65, 67, 69, 70, 72, 69, 74, &
      74, 69, 69, 74, 74, 74, 74, 74, 75, 70, 71, 65, &
      67, 75, 75, 75, 75, 75, 75, 75, 66, 66, 66, 67, &
      70, 66, 67, 64, 66, 67, 70, 66, 67, 58, 68, 72, &
      64, 66, 67, 66, 67, 70, 76, 76, 76, 76, 76 /)
  INTEGER, PARAMETER, DIMENSION(659) :: IHESS_K = (/&
    IHESS_K_0, IHESS_K_1 /)


END MODULE cbm42_strato_SOA_HessianSP

