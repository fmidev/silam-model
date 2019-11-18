MODULE cbm4_SOA_Global
  
  !$use omp_lib
  USE cbm4_SOA_Parameters, ONLY: sp, NSPEC, NVAR, NFIX, NREACT
  PUBLIC
  SAVE

  REAL(kind=sp) :: FIX(NFIX)
  !$OMP THREADPRIVATE(FIX)
  REAL(kind=sp) :: RCONST(NREACT)
  !$OMP THREADPRIVATE(RCONST)
  
END MODULE cbm4_SOA_Global

