! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! 
! Sparse Jacobian Data Structures File
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
! File                 : cbm5_strato_SOA_JacobianSP.f90
! Time                 : Wed Nov  6 16:47:03 2019
! Working directory    : /home/hanniner/Silam/silam_v5_7-risto/kpp/cbm5_strato_SOA
! Equation file        : cbm5_strato_SOA.kpp
! Output root filename : cbm5_strato_SOA
! 
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



MODULE cbm5_strato_SOA_JacobianSP

  PUBLIC
  SAVE


! Sparse Jacobian Data


  INTEGER, PARAMETER, DIMENSION(360) :: LU_IROW_0 = (/ &
       1,  2,  3,  4,  5,  5,  5,  6,  6,  6,  6,  6, &
       6,  7,  7,  7,  7,  7,  7,  8,  8,  8,  8,  8, &
       8,  9,  9,  9,  9,  9,  9, 10, 10, 11, 11, 11, &
      12, 12, 12, 13, 13, 14, 14, 14, 15, 15, 15, 15, &
      15, 15, 15, 15, 15, 16, 16, 16, 16, 16, 16, 16, &
      16, 16, 17, 17, 17, 17, 17, 17, 17, 17, 17, 18, &
      18, 18, 18, 18, 18, 18, 19, 19, 19, 19, 20, 20, &
      21, 21, 21, 22, 22, 23, 23, 24, 24, 24, 25, 25, &
      25, 26, 26, 27, 27, 27, 28, 28, 28, 29, 29, 29, &
      29, 29, 30, 30, 30, 30, 31, 31, 31, 31, 31, 31, &
      32, 32, 32, 32, 33, 33, 33, 33, 33, 34, 34, 34, &
      34, 34, 34, 34, 35, 35, 35, 35, 35, 36, 36, 36, &
      36, 37, 37, 37, 37, 37, 38, 38, 38, 39, 39, 39, &
      40, 40, 40, 40, 41, 41, 41, 42, 42, 42, 42, 43, &
      43, 43, 43, 43, 43, 44, 44, 44, 44, 44, 45, 45, &
      45, 45, 45, 46, 46, 46, 46, 46, 47, 47, 47, 47, &
      47, 47, 48, 48, 48, 48, 49, 49, 49, 49, 50, 50, &
      50, 50, 50, 50, 50, 50, 51, 51, 51, 51, 51, 51, &
      51, 51, 51, 51, 52, 52, 52, 52, 52, 52, 52, 52, &
      53, 53, 53, 53, 53, 53, 53, 53, 53, 53, 53, 53, &
      54, 54, 54, 54, 54, 54, 55, 55, 55, 55, 55, 56, &
      56, 56, 56, 56, 56, 56, 56, 56, 56, 56, 56, 56, &
      56, 56, 56, 56, 56, 56, 56, 56, 57, 57, 57, 57, &
      57, 57, 58, 58, 58, 58, 58, 58, 58, 58, 58, 58, &
      58, 58, 59, 59, 59, 59, 59, 60, 60, 60, 60, 60, &
      60, 60, 60, 60, 60, 61, 61, 61, 61, 61, 61, 62, &
      62, 62, 62, 62, 62, 62, 62, 62, 62, 63, 63, 63, &
      63, 63, 63, 64, 64, 64, 64, 64, 64, 64, 65, 65, &
      65, 65, 65, 65, 65, 65, 65, 66, 66, 66, 66, 66, &
      66, 66, 66, 66, 66, 66, 66, 66, 66, 66, 66, 66 /)
  INTEGER, PARAMETER, DIMENSION(360) :: LU_IROW_1 = (/ &
      66, 67, 67, 67, 67, 67, 67, 67, 67, 68, 68, 68, &
      68, 68, 68, 68, 68, 68, 68, 68, 68, 68, 68, 68, &
      68, 68, 68, 68, 68, 68, 69, 69, 69, 69, 69, 69, &
      69, 69, 69, 69, 69, 69, 69, 69, 69, 69, 69, 69, &
      69, 69, 69, 69, 69, 69, 70, 70, 70, 70, 70, 70, &
      70, 71, 71, 71, 71, 71, 71, 71, 71, 71, 71, 71, &
      71, 71, 71, 71, 71, 72, 72, 72, 72, 72, 72, 72, &
      72, 72, 72, 72, 72, 72, 72, 72, 72, 72, 72, 72, &
      72, 73, 73, 73, 73, 73, 73, 73, 73, 73, 73, 73, &
      73, 73, 73, 73, 73, 73, 73, 73, 73, 73, 73, 73, &
      73, 73, 74, 74, 74, 74, 74, 74, 74, 74, 74, 74, &
      74, 74, 74, 74, 74, 74, 74, 74, 75, 75, 75, 75, &
      75, 75, 75, 75, 75, 75, 75, 75, 75, 75, 75, 75, &
      75, 75, 75, 75, 75, 75, 76, 76, 76, 76, 76, 76, &
      76, 76, 76, 76, 76, 76, 76, 76, 76, 76, 76, 76, &
      76, 76, 76, 76, 76, 76, 76, 77, 77, 77, 77, 77, &
      77, 77, 77, 77, 77, 77, 77, 77, 77, 77, 77, 77, &
      77, 77, 77, 77, 77, 77, 77, 77, 77, 78, 78, 78, &
      78, 78, 78, 78, 78, 78, 78, 78, 78, 78, 78, 78, &
      78, 78, 78, 78, 78, 78, 78, 78, 78, 78, 78, 78, &
      79, 79, 79, 79, 79, 79, 79, 79, 79, 79, 79, 79, &
      79, 79, 79, 79, 79, 79, 79, 79, 79, 79, 79, 79, &
      79, 79, 80, 80, 80, 80, 80, 80, 80, 80, 80, 80, &
      80, 80, 80, 80, 80, 80, 80, 80, 80, 80, 80, 80, &
      80, 80, 80, 80, 80, 80, 80, 80, 80, 80, 80, 80, &
      80, 80, 80, 80, 80, 80, 80, 80, 80, 80, 80, 80, &
      80, 80, 80, 80, 80, 80, 80, 80, 80, 80, 80, 80, &
      80, 80, 80, 80, 80, 81, 81, 81, 81, 81, 81, 81, &
      81, 81, 81, 81, 81, 81, 81, 81, 81, 81, 81, 81, &
      81, 81, 81, 81, 81, 81, 81, 81, 81, 81, 81, 81 /)
  INTEGER, PARAMETER, DIMENSION(349) :: LU_IROW_2 = (/ &
      81, 81, 81, 81, 81, 81, 81, 81, 81, 82, 82, 82, &
      82, 82, 82, 82, 82, 82, 82, 82, 82, 82, 82, 82, &
      82, 82, 82, 82, 82, 82, 82, 82, 82, 82, 82, 82, &
      82, 82, 82, 82, 82, 82, 83, 83, 83, 83, 83, 83, &
      83, 83, 83, 83, 83, 83, 83, 83, 83, 83, 83, 83, &
      83, 83, 83, 83, 83, 83, 83, 83, 83, 83, 83, 83, &
      84, 84, 84, 84, 84, 84, 84, 84, 84, 84, 84, 84, &
      84, 84, 84, 84, 84, 84, 84, 84, 84, 84, 84, 84, &
      84, 84, 84, 84, 84, 84, 84, 84, 84, 84, 84, 84, &
      84, 84, 84, 84, 84, 84, 84, 84, 84, 84, 84, 84, &
      84, 85, 85, 85, 85, 85, 85, 85, 85, 85, 85, 85, &
      85, 85, 85, 85, 85, 85, 85, 85, 85, 85, 85, 85, &
      85, 85, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, &
      86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, &
      86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, &
      86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, &
      87, 87, 87, 87, 87, 87, 87, 87, 87, 87, 87, 87, &
      87, 87, 87, 87, 87, 87, 87, 87, 87, 87, 87, 87, &
      87, 87, 87, 87, 87, 87, 87, 87, 87, 87, 87, 87, &
      88, 88, 88, 88, 88, 88, 88, 88, 88, 88, 88, 88, &
      88, 88, 88, 88, 88, 88, 88, 88, 88, 88, 88, 88, &
      88, 88, 88, 88, 88, 89, 89, 89, 89, 89, 89, 89, &
      89, 89, 89, 89, 89, 89, 89, 89, 89, 89, 89, 89, &
      89, 89, 89, 90, 90, 90, 90, 90, 90, 90, 90, 90, &
      90, 90, 90, 90, 90, 90, 90, 90, 90, 90, 90, 90, &
      90, 90, 90, 90, 90, 90, 90, 90, 90, 90, 90, 90, &
      90, 90, 90, 90, 90, 91, 91, 91, 91, 91, 91, 91, &
      91, 91, 91, 91, 91, 91, 91, 91, 91, 91, 91, 91, &
      91, 91, 91, 91, 91, 91, 91, 91, 91, 91, 91, 91, &
      91 /)
  INTEGER, PARAMETER, DIMENSION(1069) :: LU_IROW = (/&
    LU_IROW_0, LU_IROW_1, LU_IROW_2 /)

  INTEGER, PARAMETER, DIMENSION(360) :: LU_ICOL_0 = (/ &
       1,  2,  3,  4,  5,  6, 80,  6,  7, 10, 22, 26, &
      80,  7,  8, 10, 22, 26, 80,  8,  9, 10, 22, 26, &
      80,  9, 10, 11, 22, 26, 80, 10, 80, 11, 12, 80, &
      12, 13, 80, 13, 80, 14, 15, 80, 15, 16, 19, 59, &
      70, 80, 82, 83, 87, 16, 17, 19, 59, 70, 80, 82, &
      83, 87, 17, 18, 19, 59, 70, 80, 82, 83, 87, 18, &
      19, 59, 80, 82, 83, 87, 19, 80, 83, 87, 20, 80, &
      21, 76, 89, 22, 80, 23, 80, 24, 76, 81, 25, 86, &
      87, 26, 80, 27, 86, 88, 28, 75, 86, 27, 29, 77, &
      86, 88, 30, 80, 86, 88, 31, 49, 61, 80, 83, 88, &
      32, 80, 81, 86, 33, 42, 79, 80, 89, 34, 74, 75, &
      80, 84, 85, 90, 35, 74, 75, 80, 84, 36, 80, 84, &
      86, 22, 26, 37, 80, 88, 38, 80, 81, 39, 80, 81, &
      40, 80, 81, 85, 41, 80, 81, 42, 79, 86, 89, 43, &
      62, 80, 84, 86, 87, 44, 49, 80, 84, 85, 45, 71, &
      80, 84, 90, 46, 80, 81, 82, 84, 26, 47, 54, 67, &
      80, 83, 48, 74, 80, 86, 49, 84, 88, 91, 38, 50, &
      61, 63, 64, 70, 80, 81, 33, 42, 51, 79, 80, 82, &
      84, 85, 86, 89, 24, 52, 57, 58, 76, 77, 80, 81, &
      25, 53, 62, 66, 67, 73, 78, 80, 84, 86, 87, 91, &
      37, 54, 62, 80, 83, 88, 55, 72, 80, 81, 86, 47, &
      50, 54, 56, 59, 61, 62, 63, 64, 67, 70, 73, 78, &
      79, 80, 81, 82, 83, 87, 88, 91, 57, 76, 80, 81, &
      82, 86, 32, 52, 57, 58, 65, 76, 77, 80, 81, 82, &
      84, 86, 59, 80, 82, 83, 87, 60, 73, 77, 78, 79, &
      80, 82, 84, 89, 91, 61, 80, 81, 82, 83, 87, 22, &
      26, 37, 43, 62, 80, 84, 86, 87, 88, 63, 80, 81, &
      82, 83, 87, 63, 64, 80, 81, 82, 83, 87, 65, 76, &
      79, 80, 81, 82, 83, 88, 89, 37, 43, 55, 59, 62, &
      66, 67, 70, 71, 72, 80, 81, 82, 83, 84, 86, 87 /)
  INTEGER, PARAMETER, DIMENSION(360) :: LU_ICOL_1 = (/ &
      88, 67, 70, 80, 81, 82, 83, 86, 87, 38, 39, 40, &
      41, 46, 63, 68, 70, 72, 73, 76, 77, 78, 80, 81, &
      82, 83, 84, 85, 87, 91, 31, 49, 60, 61, 68, 69, &
      70, 72, 73, 76, 77, 78, 79, 80, 81, 82, 83, 84, &
      85, 86, 87, 88, 89, 91, 70, 80, 81, 82, 83, 86, &
      87, 39, 55, 59, 64, 70, 71, 72, 80, 81, 82, 83, &
      84, 86, 87, 88, 90, 26, 37, 55, 59, 62, 63, 64, &
      67, 70, 71, 72, 80, 81, 82, 83, 84, 86, 87, 88, &
      90, 41, 45, 54, 55, 59, 61, 62, 63, 64, 66, 67, &
      70, 71, 72, 73, 79, 80, 81, 82, 83, 84, 86, 87, &
      88, 90, 48, 59, 67, 70, 73, 74, 75, 79, 80, 81, &
      82, 83, 84, 85, 86, 87, 88, 90, 28, 35, 47, 54, &
      62, 67, 70, 74, 75, 78, 79, 80, 81, 82, 83, 84, &
      85, 86, 87, 88, 89, 90, 24, 52, 57, 58, 65, 68, &
      70, 72, 73, 76, 77, 78, 79, 80, 81, 82, 83, 84, &
      85, 86, 87, 88, 89, 90, 91, 29, 52, 57, 58, 60, &
      65, 68, 70, 72, 73, 76, 77, 78, 79, 80, 81, 82, &
      83, 84, 85, 86, 87, 88, 89, 90, 91, 39, 41, 45, &
      48, 55, 63, 64, 66, 67, 70, 71, 72, 74, 75, 78, &
      79, 80, 81, 82, 83, 84, 85, 86, 87, 88, 89, 90, &
       1, 21, 23, 33, 42, 51, 60, 65, 73, 75, 76, 77, &
      78, 79, 80, 81, 82, 83, 84, 85, 86, 87, 88, 89, &
      90, 91, 20, 22, 23, 26, 30, 31, 32, 33, 34, 35, &
      36, 38, 39, 40, 41, 42, 44, 45, 46, 47, 48, 49, &
      50, 51, 52, 53, 54, 56, 57, 58, 59, 60, 61, 62, &
      63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, &
      75, 76, 77, 78, 79, 80, 81, 82, 83, 84, 85, 86, &
      87, 88, 89, 90, 91,  2,  3,  4, 20, 21, 24, 32, &
      38, 39, 40, 41, 46, 50, 52, 57, 58, 61, 63, 64, &
      65, 68, 70, 72, 73, 76, 77, 78, 79, 80, 81, 82 /)
  INTEGER, PARAMETER, DIMENSION(349) :: LU_ICOL_2 = (/ &
      83, 84, 85, 86, 87, 88, 89, 90, 91, 27, 29, 46, &
      51, 57, 58, 59, 60, 61, 63, 64, 65, 68, 69, 70, &
      72, 73, 76, 77, 78, 79, 80, 81, 82, 83, 84, 85, &
      86, 87, 88, 89, 90, 91, 54, 59, 61, 62, 63, 64, &
      65, 67, 69, 70, 72, 73, 74, 75, 76, 77, 78, 79, &
      80, 81, 82, 83, 84, 85, 86, 87, 88, 89, 90, 91, &
      22, 26, 36, 37, 38, 39, 40, 41, 43, 44, 45, 46, &
      47, 48, 49, 50, 54, 55, 56, 59, 61, 62, 63, 64, &
      66, 67, 69, 70, 71, 72, 73, 74, 75, 76, 77, 78, &
      79, 80, 81, 82, 83, 84, 85, 86, 87, 88, 89, 90, &
      91, 23, 28, 34, 35, 38, 44, 49, 73, 74, 75, 77, &
      78, 79, 80, 81, 82, 83, 84, 85, 86, 87, 88, 89, &
      90, 91,  1, 25, 27, 28, 30, 32, 36, 37, 42, 43, &
      48, 49, 53, 55, 57, 59, 61, 62, 63, 64, 65, 66, &
      67, 69, 70, 71, 72, 73, 74, 75, 76, 77, 78, 79, &
      80, 81, 82, 83, 84, 85, 86, 87, 88, 89, 90, 91, &
      25, 28, 36, 42, 48, 53, 57, 59, 61, 62, 63, 64, &
      66, 67, 70, 71, 72, 73, 74, 75, 76, 77, 78, 79, &
      80, 81, 82, 83, 84, 85, 86, 87, 88, 89, 90, 91, &
      27, 29, 30, 37, 49, 65, 69, 70, 71, 72, 73, 74, &
      75, 76, 77, 78, 79, 80, 81, 82, 83, 84, 85, 86, &
      87, 88, 89, 90, 91, 42, 51, 60, 65, 73, 75, 76, &
      77, 78, 79, 80, 81, 82, 83, 84, 85, 86, 87, 88, &
      89, 90, 91, 22, 23, 26, 38, 39, 41, 44, 45, 47, &
      48, 49, 54, 55, 59, 61, 62, 63, 64, 67, 70, 71, &
      72, 74, 75, 78, 79, 80, 81, 82, 83, 84, 85, 86, &
      87, 88, 89, 90, 91, 23, 40, 41, 44, 49, 54, 59, &
      61, 62, 63, 64, 66, 67, 70, 71, 72, 74, 75, 78, &
      79, 80, 81, 82, 83, 84, 85, 86, 87, 88, 89, 90, &
      91 /)
  INTEGER, PARAMETER, DIMENSION(1069) :: LU_ICOL = (/&
    LU_ICOL_0, LU_ICOL_1, LU_ICOL_2 /)

  INTEGER, PARAMETER, DIMENSION(92) :: LU_CROW = (/ &
       1,  2,  3,  4,  5,  8, 14, 20, 26, 32, 34, 37, &
      40, 42, 45, 54, 63, 72, 79, 83, 85, 88, 90, 92, &
      95, 98,100,103,106,111,115,121,125,130,137,142, &
     146,151,154,157,161,164,168,174,179,184,189,195, &
     199,203,211,221,229,241,247,252,273,279,291,296, &
     306,312,322,328,335,344,362,370,391,415,422,438, &
     458,483,501,523,548,574,601,627,690,730,763,793, &
     842,867,913,949,978,1000,1038,1070 /)

  INTEGER, PARAMETER, DIMENSION(92) :: LU_DIAG = (/ &
       1,  2,  3,  4,  5,  8, 14, 20, 26, 32, 34, 37, &
      40, 42, 45, 54, 63, 72, 79, 83, 85, 88, 90, 92, &
      95, 98,100,103,107,111,115,121,125,130,137,142, &
     148,151,154,157,161,164,168,174,179,184,190,195, &
     199,204,213,222,230,242,247,255,273,282,291,296, &
     306,316,322,329,335,349,362,376,396,415,427,448, &
     472,488,509,532,559,588,614,678,719,753,784,834, &
     860,907,944,974,997,1036,1069,1070 /)


END MODULE cbm5_strato_SOA_JacobianSP
