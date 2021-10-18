! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! 
! Parameter Module File
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
! File                 : cbm4_SOA_Parameters.f90
! Time                 : Fri Apr 23 13:17:27 2021
! Working directory    : /fmi/scratch/project_2001411/risto/silam_v5_7/kpp/cbm4_SOA
! Equation file        : cbm4_SOA.kpp
! Output root filename : cbm4_SOA
! 
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



MODULE cbm4_SOA_Parameters

  USE cbm4_SOA_Precision
  PUBLIC
  SAVE


! NSPEC - Number of chemical species
  INTEGER, PARAMETER :: NSPEC = 54 
! NVAR - Number of Variable species
  INTEGER, PARAMETER :: NVAR = 50 
! NVARACT - Number of Active species
  INTEGER, PARAMETER :: NVARACT = 47 
! NFIX - Number of Fixed species
  INTEGER, PARAMETER :: NFIX = 4 
! NREACT - Number of reactions
  INTEGER, PARAMETER :: NREACT = 133 
! NVARST - Starting of variables in conc. vect.
  INTEGER, PARAMETER :: NVARST = 1 
! NFIXST - Starting of fixed in conc. vect.
  INTEGER, PARAMETER :: NFIXST = 51 
! NONZERO - Number of nonzero entries in Jacobian
  INTEGER, PARAMETER :: NONZERO = 405 
! LU_NONZERO - Number of nonzero entries in LU factoriz. of Jacobian
  INTEGER, PARAMETER :: LU_NONZERO = 426 
! CNVAR - (NVAR+1) Number of elements in compressed row format
  INTEGER, PARAMETER :: CNVAR = 51 
! NHESS - Length of Sparse Hessian
  INTEGER, PARAMETER :: NHESS = 408 
! NLOOKAT - Number of species to look at
  INTEGER, PARAMETER :: NLOOKAT = 54 
! NMONITOR - Number of species to monitor
  INTEGER, PARAMETER :: NMONITOR = 1 
! NMASS - Number of atoms to check mass balance
  INTEGER, PARAMETER :: NMASS = 1 

! Index declaration for variable species in C and VAR
!   VAR(ind_spc) = C(ind_spc)

  INTEGER, PARAMETER :: ind_NTR = 1 
  INTEGER, PARAMETER :: ind_AVB0 = 2 
  INTEGER, PARAMETER :: ind_AVB1e0 = 3 
  INTEGER, PARAMETER :: ind_AVB1e1 = 4 
  INTEGER, PARAMETER :: ind_AVB1e2 = 5 
  INTEGER, PARAMETER :: ind_AVB1e3 = 6 
  INTEGER, PARAMETER :: ind_AVB1e4 = 7 
  INTEGER, PARAMETER :: ind_AVB1e5 = 8 
  INTEGER, PARAMETER :: ind_AVB1e6 = 9 
  INTEGER, PARAMETER :: ind_BVB0 = 10 
  INTEGER, PARAMETER :: ind_BVB1e0 = 11 
  INTEGER, PARAMETER :: ind_BVB1e1 = 12 
  INTEGER, PARAMETER :: ind_BVB1e2 = 13 
  INTEGER, PARAMETER :: ind_BVB1e3 = 14 
  INTEGER, PARAMETER :: ind_O1D = 15 
  INTEGER, PARAMETER :: ind_MEOH = 16 
  INTEGER, PARAMETER :: ind_ETOH = 17 
  INTEGER, PARAMETER :: ind_CRO = 18 
  INTEGER, PARAMETER :: ind_PAN = 19 
  INTEGER, PARAMETER :: ind_H2O2 = 20 
  INTEGER, PARAMETER :: ind_TOL = 21 
  INTEGER, PARAMETER :: ind_N2O5 = 22 
  INTEGER, PARAMETER :: ind_XYL = 23 
  INTEGER, PARAMETER :: ind_HONO = 24 
  INTEGER, PARAMETER :: ind_PNA = 25 
  INTEGER, PARAMETER :: ind_TO2 = 26 
  INTEGER, PARAMETER :: ind_HNO3 = 27 
  INTEGER, PARAMETER :: ind_ROR4 = 28 
  INTEGER, PARAMETER :: ind_CRES = 29 
  INTEGER, PARAMETER :: ind_MGLY = 30 
  INTEGER, PARAMETER :: ind_CO = 31 
  INTEGER, PARAMETER :: ind_ETH = 32 
  INTEGER, PARAMETER :: ind_OPEN = 33 
  INTEGER, PARAMETER :: ind_XO2N = 34 
  INTEGER, PARAMETER :: ind_XO2 = 35 
  INTEGER, PARAMETER :: ind_C5H8_2 = 36 
  INTEGER, PARAMETER :: ind_PAR4 = 37 
  INTEGER, PARAMETER :: ind_OLE4 = 38 
  INTEGER, PARAMETER :: ind_HCHO = 39 
  INTEGER, PARAMETER :: ind_C5H8 = 40 
  INTEGER, PARAMETER :: ind_ISPD = 41 
  INTEGER, PARAMETER :: ind_ALD2 = 42 
  INTEGER, PARAMETER :: ind_O3 = 43 
  INTEGER, PARAMETER :: ind_O = 44 
  INTEGER, PARAMETER :: ind_HO2 = 45 
  INTEGER, PARAMETER :: ind_NO3 = 46 
  INTEGER, PARAMETER :: ind_C2O3 = 47 
  INTEGER, PARAMETER :: ind_NO = 48 
  INTEGER, PARAMETER :: ind_NO2 = 49 
  INTEGER, PARAMETER :: ind_OH = 50 

! Index declaration for fixed species in C
!   C(ind_spc)

  INTEGER, PARAMETER :: ind_H2O = 51 
  INTEGER, PARAMETER :: ind_O2 = 52 
  INTEGER, PARAMETER :: ind_CH4 = 53 
  INTEGER, PARAMETER :: ind_M = 54 

! Index declaration for fixed species in FIX
!    FIX(indf_spc) = C(ind_spc) = C(NVAR+indf_spc)

  INTEGER, PARAMETER :: indf_H2O = 1 
  INTEGER, PARAMETER :: indf_O2 = 2 
  INTEGER, PARAMETER :: indf_CH4 = 3 
  INTEGER, PARAMETER :: indf_M = 4 

END MODULE cbm4_SOA_Parameters

