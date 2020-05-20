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
! File                 : cbm5_strato_SOA_HessianSP.f90
! Time                 : Wed Nov  6 16:47:03 2019
! Working directory    : /home/hanniner/Silam/silam_v5_7-risto/kpp/cbm5_strato_SOA
! Equation file        : cbm5_strato_SOA.kpp
! Output root filename : cbm5_strato_SOA
! 
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



MODULE cbm5_strato_SOA_HessianSP

  PUBLIC
  SAVE


! Hessian Sparse Data
! 

  INTEGER, PARAMETER, DIMENSION(360) :: IHESS_I_0 = (/ &
       5,  6,  6,  6,  6,  6,  7,  7,  7,  7,  7,  8, &
       8,  8,  8,  8,  9,  9,  9,  9,  9, 10, 11, 11, &
      12, 12, 13, 14, 15, 15, 15, 15, 15, 15, 15, 15, &
      15, 15, 15, 15, 15, 16, 16, 16, 16, 16, 16, 16, &
      16, 16, 16, 16, 16, 16, 17, 17, 17, 17, 17, 17, &
      17, 17, 17, 17, 17, 17, 17, 18, 18, 18, 18, 18, &
      18, 18, 18, 19, 19, 19, 20, 21, 22, 23, 24, 24, &
      25, 26, 27, 27, 28, 29, 29, 30, 30, 30, 30, 31, &
      31, 31, 32, 32, 33, 33, 33, 34, 34, 34, 34, 34, &
      34, 34, 35, 35, 35, 36, 36, 37, 37, 37, 38, 38, &
      39, 39, 40, 40, 40, 41, 41, 42, 42, 43, 43, 43, &
      43, 44, 44, 44, 45, 45, 45, 46, 46, 46, 46, 46, &
      47, 47, 47, 47, 47, 48, 48, 49, 49, 49, 50, 50, &
      50, 50, 50, 50, 50, 51, 51, 51, 51, 52, 52, 52, &
      52, 52, 52, 53, 53, 53, 53, 53, 53, 53, 53, 53, &
      54, 54, 54, 54, 55, 55, 55, 56, 56, 56, 56, 56, &
      56, 56, 56, 56, 56, 56, 56, 56, 56, 56, 56, 56, &
      56, 56, 56, 57, 57, 57, 57, 58, 58, 58, 58, 58, &
      58, 58, 58, 59, 59, 59, 59, 60, 60, 60, 60, 60, &
      60, 60, 60, 61, 61, 61, 61, 61, 62, 62, 62, 62, &
      62, 63, 63, 63, 63, 63, 64, 64, 64, 64, 64, 64, &
      65, 65, 65, 65, 65, 65, 65, 65, 66, 66, 66, 66, &
      66, 66, 66, 66, 66, 67, 67, 67, 67, 67, 67, 67, &
      67, 67, 68, 68, 68, 68, 68, 68, 68, 68, 68, 68, &
      68, 68, 68, 68, 68, 68, 69, 69, 69, 69, 69, 69, &
      69, 70, 70, 70, 70, 70, 70, 71, 71, 71, 71, 71, &
      71, 71, 71, 71, 71, 71, 71, 71, 71, 72, 72, 72, &
      72, 72, 72, 72, 72, 72, 72, 72, 72, 72, 72, 72, &
      72, 72, 72, 72, 72, 72, 72, 72, 72, 72, 72, 72, &
      72, 72, 73, 73, 73, 73, 73, 73, 73, 73, 73, 73 /)
  INTEGER, PARAMETER, DIMENSION(360) :: IHESS_I_1 = (/ &
      73, 73, 73, 73, 73, 73, 73, 73, 73, 73, 73, 73, &
      73, 73, 73, 73, 73, 73, 73, 73, 73, 74, 74, 74, &
      74, 74, 74, 74, 74, 74, 74, 74, 74, 74, 74, 74, &
      74, 74, 75, 75, 75, 75, 75, 75, 75, 75, 75, 75, &
      75, 75, 75, 75, 75, 75, 75, 75, 75, 76, 76, 76, &
      76, 76, 76, 76, 76, 76, 76, 76, 76, 76, 76, 76, &
      76, 76, 76, 76, 76, 76, 77, 77, 77, 77, 77, 78, &
      78, 78, 78, 78, 78, 78, 78, 78, 78, 78, 78, 78, &
      78, 78, 78, 78, 78, 78, 78, 78, 78, 78, 78, 78, &
      78, 78, 78, 78, 78, 78, 79, 79, 79, 79, 79, 79, &
      79, 79, 79, 79, 79, 79, 79, 79, 79, 79, 79, 79, &
      79, 79, 80, 80, 80, 80, 80, 80, 80, 80, 80, 80, &
      80, 80, 80, 80, 80, 80, 80, 80, 80, 80, 80, 80, &
      80, 80, 80, 80, 80, 80, 80, 80, 80, 80, 80, 80, &
      80, 80, 80, 80, 80, 80, 80, 80, 80, 80, 80, 80, &
      80, 80, 80, 80, 80, 80, 80, 80, 80, 80, 80, 80, &
      80, 80, 80, 80, 80, 80, 80, 80, 80, 80, 80, 80, &
      80, 80, 80, 80, 80, 80, 80, 80, 80, 80, 80, 81, &
      81, 81, 81, 81, 81, 81, 81, 81, 81, 81, 81, 81, &
      81, 81, 81, 81, 81, 81, 81, 81, 81, 81, 81, 81, &
      81, 81, 81, 81, 81, 81, 81, 81, 82, 82, 82, 82, &
      82, 82, 82, 82, 82, 82, 82, 82, 82, 82, 82, 82, &
      82, 82, 82, 82, 82, 82, 82, 82, 82, 82, 82, 82, &
      82, 82, 82, 82, 83, 83, 83, 83, 83, 83, 83, 83, &
      83, 83, 83, 83, 83, 83, 83, 83, 83, 83, 83, 83, &
      84, 84, 84, 84, 84, 84, 84, 84, 84, 84, 84, 84, &
      84, 84, 84, 84, 84, 84, 84, 84, 84, 84, 84, 84, &
      84, 84, 84, 84, 84, 84, 84, 84, 84, 84, 84, 84, &
      84, 84, 84, 84, 84, 84, 84, 84, 84, 84, 84, 84, &
      84, 84, 84, 84, 84, 84, 84, 84, 84, 84, 84, 84 /)
  INTEGER, PARAMETER, DIMENSION(239) :: IHESS_I_2 = (/ &
      84, 84, 84, 84, 84, 84, 84, 84, 84, 84, 84, 84, &
      84, 84, 84, 84, 84, 84, 84, 84, 84, 84, 84, 84, &
      84, 85, 85, 85, 85, 85, 85, 85, 85, 85, 85, 85, &
      85, 85, 85, 86, 86, 86, 86, 86, 86, 86, 86, 86, &
      86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, &
      86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, &
      86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, &
      86, 86, 86, 87, 87, 87, 87, 87, 87, 87, 87, 87, &
      87, 87, 87, 87, 87, 87, 87, 87, 87, 87, 87, 87, &
      87, 87, 87, 87, 87, 87, 88, 88, 88, 88, 88, 88, &
      88, 88, 88, 88, 88, 88, 88, 88, 88, 88, 88, 88, &
      88, 88, 88, 88, 88, 88, 89, 89, 89, 89, 89, 89, &
      89, 89, 89, 89, 89, 89, 89, 89, 90, 90, 90, 90, &
      90, 90, 90, 90, 90, 90, 90, 90, 90, 90, 90, 90, &
      90, 90, 90, 90, 90, 90, 90, 90, 90, 90, 90, 90, &
      90, 90, 90, 90, 90, 90, 90, 90, 90, 90, 90, 90, &
      90, 90, 90, 90, 90, 90, 90, 90, 90, 90, 90, 91, &
      91, 91, 91, 91, 91, 91, 91, 91, 91, 91, 91, 91, &
      91, 91, 91, 91, 91, 91, 91, 91, 91, 91, 91, 91, &
      91, 91, 91, 91, 91, 91, 91, 91, 91, 91, 91 /)
  INTEGER, PARAMETER, DIMENSION(959) :: IHESS_I = (/&
    IHESS_I_0, IHESS_I_1, IHESS_I_2 /)

  INTEGER, PARAMETER, DIMENSION(360) :: IHESS_J_0 = (/ &
       6,  6,  7, 10, 22, 26,  7,  8, 10, 22, 26,  8, &
       9, 10, 22, 26,  9, 10, 11, 22, 26, 10, 11, 12, &
      12, 13, 13, 15, 15, 16, 19, 19, 19, 59, 59, 59, &
      59, 70, 70, 70, 70, 16, 17, 19, 19, 19, 59, 59, &
      59, 59, 70, 70, 70, 70, 17, 18, 19, 19, 19, 59, &
      59, 59, 59, 70, 70, 70, 70, 18, 19, 19, 19, 59, &
      59, 59, 59, 19, 19, 19, 20, 76, 22, 23, 24, 76, &
      86, 26, 27, 27, 75, 27, 29, 30, 30, 80, 86, 31, &
      49, 61, 32, 81, 33, 42, 89, 34, 74, 74, 74, 75, &
      75, 75, 35, 74, 75, 36, 84, 22, 26, 37, 38, 38, &
      39, 39, 40, 40, 85, 41, 41, 42, 86, 43, 43, 62, &
      62, 44, 49, 84, 45, 71, 84, 46, 46, 46, 80, 84, &
      26, 47, 54, 67, 67, 48, 74, 49, 49, 84, 38, 38, &
      50, 61, 63, 64, 70, 33, 51, 84, 85, 24, 52, 52, &
      57, 58, 76, 53, 62, 66, 67, 73, 78, 80, 84, 87, &
      37, 54, 54, 62, 55, 72, 72, 50, 54, 54, 56, 59, &
      61, 61, 63, 63, 64, 64, 67, 67, 67, 70, 79, 80, &
      81, 82, 87, 57, 57, 57, 76, 32, 52, 57, 58, 58, &
      58, 65, 76, 59, 59, 59, 59, 60, 60, 60, 73, 78, &
      79, 79, 80, 61, 61, 61, 61, 61, 22, 26, 43, 62, &
      62, 63, 63, 63, 63, 63, 63, 64, 64, 64, 64, 64, &
      65, 65, 65, 65, 65, 65, 76, 76, 37, 43, 55, 59, &
      66, 67, 70, 70, 71, 67, 67, 67, 70, 70, 70, 70, &
      70, 70, 38, 39, 40, 41, 46, 63, 68, 68, 68, 70, &
      72, 73, 76, 78, 81, 81, 31, 60, 68, 69, 69, 69, &
      80, 70, 70, 70, 70, 70, 70, 39, 39, 59, 59, 59, &
      64, 64, 70, 71, 71, 71, 71, 72, 72, 26, 37, 55, &
      59, 59, 59, 59, 62, 63, 63, 63, 63, 64, 64, 64, &
      64, 67, 67, 67, 70, 70, 70, 70, 70, 71, 71, 71, &
      72, 72, 41, 45, 54, 59, 59, 59, 59, 61, 63, 63 /)
  INTEGER, PARAMETER, DIMENSION(360) :: IHESS_J_1 = (/ &
      63, 63, 63, 64, 64, 64, 64, 64, 66, 67, 67, 70, &
      70, 70, 72, 72, 73, 73, 73, 73, 73, 59, 67, 67, &
      70, 70, 73, 73, 73, 73, 73, 74, 74, 74, 74, 74, &
      74, 74, 35, 47, 54, 54, 67, 67, 74, 75, 75, 75, &
      75, 75, 75, 75, 78, 78, 78, 78, 78, 52, 57, 58, &
      58, 65, 65, 65, 65, 65, 68, 76, 76, 76, 76, 76, &
      76, 76, 76, 81, 81, 81, 29, 52, 60, 68, 77, 39, &
      39, 41, 41, 45, 48, 63, 63, 63, 63, 63, 64, 64, &
      64, 64, 64, 66, 67, 67, 72, 72, 74, 74, 74, 74, &
      74, 78, 78, 78, 78, 78, 23, 33, 42, 60, 60, 60, &
      65, 73, 75, 76, 78, 79, 79, 79, 79, 80, 82, 85, &
      88, 89, 20, 22, 23, 26, 30, 31, 32, 33, 34, 35, &
      36, 38, 39, 40, 41, 44, 45, 46, 46, 47, 48, 50, &
      51, 52, 53, 54, 54, 56, 57, 58, 58, 58, 59, 59, &
      60, 60, 60, 61, 61, 61, 62, 63, 63, 64, 64, 64, &
      65, 66, 67, 67, 68, 68, 68, 69, 69, 69, 70, 70, &
      72, 73, 73, 74, 75, 76, 78, 78, 80, 80, 80, 80, &
      80, 80, 80, 80, 80, 81, 82, 82, 83, 84, 84, 20, &
      38, 39, 40, 41, 46, 50, 52, 52, 57, 58, 61, 63, &
      64, 65, 68, 68, 68, 70, 72, 73, 76, 76, 76, 76, &
      76, 76, 78, 81, 81, 81, 81, 81, 27, 27, 29, 46, &
      51, 57, 58, 59, 60, 60, 61, 63, 63, 64, 65, 68, &
      68, 69, 70, 73, 76, 77, 78, 80, 80, 82, 82, 82, &
      82, 82, 82, 82, 54, 59, 61, 63, 64, 65, 67, 69, &
      70, 74, 75, 77, 79, 80, 81, 82, 83, 83, 83, 83, &
      22, 26, 37, 38, 38, 39, 39, 40, 40, 41, 41, 43, &
      44, 46, 46, 46, 49, 49, 54, 54, 56, 59, 59, 59, &
      61, 61, 61, 61, 62, 63, 63, 63, 63, 63, 64, 64, &
      64, 64, 66, 67, 67, 67, 69, 70, 70, 70, 70, 70, &
      70, 71, 72, 72, 74, 74, 74, 74, 74, 75, 75, 76 /)
  INTEGER, PARAMETER, DIMENSION(239) :: IHESS_J_2 = (/ &
      76, 79, 79, 80, 80, 80, 80, 80, 81, 81, 82, 82, &
      83, 84, 84, 84, 84, 84, 84, 84, 84, 85, 85, 85, &
      87, 34, 44, 74, 74, 75, 75, 75, 75, 75, 75, 84, &
      85, 85, 85, 27, 30, 30, 32, 36, 37, 43, 48, 49, &
      55, 59, 61, 63, 64, 65, 69, 70, 70, 74, 74, 75, &
      75, 76, 76, 76, 79, 80, 80, 81, 81, 82, 82, 82, &
      83, 83, 83, 84, 84, 84, 85, 86, 86, 86, 87, 87, &
      88, 88, 88, 42, 53, 57, 57, 57, 59, 61, 62, 63, &
      64, 67, 70, 73, 76, 78, 79, 80, 81, 82, 82, 83, &
      83, 84, 86, 87, 87, 87, 27, 29, 30, 37, 49, 65, &
      69, 70, 71, 74, 75, 76, 80, 82, 82, 83, 84, 85, &
      86, 86, 87, 88, 88, 88, 51, 60, 65, 75, 76, 79, &
      79, 80, 82, 84, 85, 86, 88, 89, 22, 23, 26, 38, &
      38, 39, 39, 41, 44, 45, 47, 54, 54, 59, 59, 59, &
      61, 61, 61, 61, 62, 63, 63, 63, 64, 64, 64, 64, &
      64, 67, 67, 67, 70, 70, 70, 70, 70, 70, 71, 72, &
      72, 74, 74, 74, 74, 74, 74, 75, 84, 88, 90, 23, &
      40, 40, 41, 54, 54, 59, 59, 61, 61, 61, 61, 61, &
      63, 64, 64, 64, 64, 66, 67, 67, 67, 70, 70, 70, &
      74, 75, 79, 80, 81, 82, 84, 85, 85, 85, 87 /)
  INTEGER, PARAMETER, DIMENSION(959) :: IHESS_J = (/&
    IHESS_J_0, IHESS_J_1, IHESS_J_2 /)

  INTEGER, PARAMETER, DIMENSION(360) :: IHESS_K_0 = (/ &
      80, 80, 80, 80, 80, 80, 80, 80, 80, 80, 80, 80, &
      80, 80, 80, 80, 80, 80, 80, 80, 80, 80, 80, 80, &
      80, 80, 80, 80, 80, 80, 80, 83, 87, 80, 82, 83, &
      87, 80, 82, 83, 87, 80, 80, 80, 83, 87, 80, 82, &
      83, 87, 80, 82, 83, 87, 80, 80, 80, 83, 87, 80, &
      82, 83, 87, 80, 82, 83, 87, 80, 80, 83, 87, 80, &
      82, 83, 87, 80, 83, 87, 80, 89, 80, 80, 81, 76, &
      87, 80, 86, 88, 86, 86, 77, 30, 80, 88, 88, 80, &
      88, 83, 80, 86, 80, 79, 89, 80, 84, 85, 90, 84, &
      85, 90, 80, 84, 84, 80, 86, 80, 80, 88, 80, 81, &
      80, 81, 80, 81, 85, 80, 81, 79, 89, 84, 86, 80, &
      87, 80, 84, 85, 80, 84, 90, 80, 81, 82, 80, 84, &
      80, 80, 83, 80, 83, 80, 86, 84, 88, 91, 80, 81, &
      80, 81, 81, 81, 81, 80, 82, 89, 89, 81, 77, 80, &
      81, 81, 76, 80, 87, 80, 87, 87, 87, 86, 87, 91, &
      88, 80, 83, 80, 86, 80, 81, 80, 80, 83, 80, 83, &
      82, 83, 82, 83, 82, 83, 80, 83, 87, 83, 91, 91, &
      91, 91, 91, 80, 81, 82, 86, 80, 80, 80, 80, 81, &
      82, 80, 84, 80, 82, 83, 87, 77, 80, 82, 79, 79, &
      84, 91, 89, 80, 81, 82, 83, 87, 80, 80, 84, 80, &
      87, 80, 81, 82, 83, 87, 81, 80, 81, 82, 83, 87, &
      79, 80, 81, 82, 83, 88, 76, 89, 88, 86, 86, 87, &
      80, 87, 86, 87, 88, 80, 83, 87, 80, 81, 82, 83, &
      86, 87, 81, 81, 81, 81, 81, 81, 77, 80, 82, 81, &
      81, 81, 80, 81, 84, 91, 80, 77, 77, 83, 84, 86, &
      82, 80, 81, 82, 83, 86, 87, 80, 81, 80, 83, 87, &
      82, 87, 80, 71, 84, 88, 90, 80, 81, 80, 88, 86, &
      80, 82, 83, 87, 80, 80, 82, 83, 87, 80, 81, 82, &
      83, 80, 83, 87, 81, 82, 83, 86, 87, 71, 84, 90, &
      80, 81, 80, 80, 83, 80, 82, 83, 87, 80, 80, 81 /)
  INTEGER, PARAMETER, DIMENSION(360) :: IHESS_K_1 = (/ &
      82, 83, 87, 80, 81, 82, 83, 87, 80, 80, 87, 83, &
      86, 87, 80, 81, 79, 80, 81, 82, 87, 83, 80, 87, &
      82, 83, 79, 80, 81, 82, 87, 74, 75, 84, 85, 86, &
      88, 90, 80, 80, 80, 83, 80, 83, 75, 75, 84, 85, &
      86, 88, 89, 90, 79, 80, 81, 82, 87, 77, 82, 80, &
      82, 79, 81, 82, 83, 88, 77, 76, 80, 82, 84, 86, &
      87, 88, 89, 83, 84, 87, 77, 77, 77, 77, 83, 80, &
      81, 80, 81, 80, 80, 80, 81, 82, 83, 87, 80, 81, &
      82, 83, 87, 80, 80, 83, 80, 81, 74, 75, 85, 88, &
      90, 79, 80, 81, 82, 87, 80, 80, 79, 77, 80, 82, &
      79, 79, 89, 89, 79, 83, 84, 87, 91, 89, 89, 89, &
      89, 89, 80, 80, 80, 80, 80, 80, 80, 80, 80, 80, &
      80, 80, 80, 80, 80, 80, 80, 80, 82, 80, 80, 80, &
      82, 80, 80, 80, 83, 80, 80, 80, 81, 82, 80, 83, &
      77, 80, 82, 80, 82, 83, 80, 80, 83, 80, 82, 83, &
      80, 80, 80, 83, 77, 80, 82, 83, 84, 86, 80, 83, &
      80, 80, 82, 84, 84, 80, 80, 82, 80, 82, 83, 84, &
      86, 87, 88, 89, 91, 84, 84, 91, 84, 87, 88, 80, &
      81, 81, 81, 81, 81, 80, 77, 80, 81, 81, 81, 81, &
      81, 81, 77, 80, 82, 81, 81, 81, 76, 80, 82, 87, &
      88, 89, 81, 83, 84, 86, 87, 91, 86, 88, 77, 82, &
      82, 82, 82, 82, 77, 82, 82, 82, 83, 82, 82, 77, &
      82, 84, 82, 82, 82, 83, 82, 80, 82, 83, 84, 86, &
      87, 88, 89, 91, 83, 83, 83, 83, 83, 83, 83, 83, &
      83, 84, 84, 83, 83, 83, 83, 83, 84, 86, 87, 88, &
      80, 80, 88, 80, 81, 80, 81, 80, 81, 80, 81, 84, &
      80, 80, 81, 82, 84, 88, 80, 83, 80, 80, 83, 87, &
      80, 81, 82, 83, 80, 80, 81, 82, 83, 87, 80, 81, &
      82, 83, 80, 80, 83, 87, 84, 80, 81, 82, 83, 86, &
      87, 84, 80, 81, 74, 75, 84, 85, 88, 84, 85, 80 /)
  INTEGER, PARAMETER, DIMENSION(239) :: IHESS_K_2 = (/ &
      84, 84, 91, 83, 84, 87, 89, 91, 84, 91, 84, 91, &
      84, 84, 85, 86, 87, 88, 89, 90, 91, 85, 88, 89, &
      91, 80, 80, 75, 85, 75, 84, 85, 88, 89, 90, 85, &
      85, 88, 89, 86, 30, 80, 80, 80, 88, 86, 80, 88, &
      86, 87, 87, 87, 87, 88, 86, 86, 87, 86, 88, 86, &
      88, 86, 87, 88, 87, 86, 87, 86, 87, 86, 87, 88, &
      86, 87, 88, 86, 87, 88, 88, 87, 88, 89, 87, 88, &
      88, 89, 90, 79, 80, 80, 81, 82, 87, 87, 87, 87, &
      87, 87, 87, 87, 87, 87, 87, 87, 87, 86, 87, 86, &
      87, 87, 87, 87, 88, 91, 88, 77, 30, 88, 88, 88, &
      86, 86, 88, 88, 88, 88, 88, 86, 88, 88, 88, 88, &
      87, 88, 88, 88, 89, 90, 82, 77, 79, 89, 89, 83, &
      87, 89, 89, 89, 89, 89, 89, 89, 80, 80, 80, 80, &
      81, 80, 81, 80, 80, 80, 80, 80, 83, 80, 83, 87, &
      80, 81, 82, 87, 80, 80, 81, 82, 80, 81, 82, 83, &
      87, 80, 83, 87, 80, 81, 82, 83, 86, 87, 90, 80, &
      81, 74, 75, 84, 85, 88, 90, 90, 90, 90, 90, 80, &
      80, 81, 80, 80, 83, 80, 83, 80, 81, 82, 83, 87, &
      83, 80, 82, 83, 87, 80, 80, 83, 87, 80, 82, 83, &
      85, 85, 91, 91, 91, 91, 91, 85, 88, 89, 91 /)
  INTEGER, PARAMETER, DIMENSION(959) :: IHESS_K = (/&
    IHESS_K_0, IHESS_K_1, IHESS_K_2 /)


END MODULE cbm5_strato_SOA_HessianSP
