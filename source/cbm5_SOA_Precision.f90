
MODULE cbm5_SOA_Precision

!
! Definition of different levels of accuracy
! for REAL variables using KIND parameterization
!
! KPP SP - Single precision kind
#ifdef DOUBLE_PRECISION
  INTEGER, PARAMETER :: sp = SELECTED_REAL_KIND(14,300)
#else
  INTEGER, PARAMETER :: sp = SELECTED_REAL_KIND(6,30)
#endif
! KPP DP - Double precision kind
  INTEGER, PARAMETER :: dp = SELECTED_REAL_KIND(14,300)
! KPP QP - Quadruple precision kind
  INTEGER, PARAMETER :: qp = SELECTED_REAL_KIND(18,400)

END MODULE cbm5_SOA_Precision


