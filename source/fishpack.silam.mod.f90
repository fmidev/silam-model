!Copyright © 2004 the University Corporation for Atmospheric Research ("UCAR"). 
!  All rights reserved. 
!  Developed by NCAR's Computational and Information Systems Laboratory, UCAR.
!
!Redistribution and use of the Software in source and binary forms, with or without 
!  modification, is permitted provided that the following conditions are met:
!
!    Neither the names of NCAR's Computational and Information Systems Laboratory, the 
!    University Corporation for Atmospheric Research, nor the names of its sponsors or 
!    contributors may be used to endorse or promote products derived from this Software 
!    without specific prior written permission.
!    Redistributions of source code must retain the above copyright notices, this list 
!    of conditions, and the disclaimer below.
!    Redistributions in binary form must reproduce the above copyright notice, this 
!    list of conditions, and the disclaimer below in the documentation and/or other 
!    materials provided with the distribution.
!
!THIS SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, 
!INCLUDING, BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR 
!PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE CONTRIBUTORS OR COPYRIGHT HOLDERS BE 
!LIABLE FOR ANY CLAIM, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES 
!OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, 
!OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS WITH THE SOFTWARE.
!
!FISHPACK is a product of the Computational & Information Systems Laboratory at the National 
!Center for Atmospheric Research (NCAR)


module fishpack_silam
  
   ! Cut-out of fishpack90 library for blktri solver
   ! Basically join of fish.f comf.f blktri.f pois3d.f
   ! from the original library with some syntax adjustments and bugfixes:
   ! Removed common blocks
   ! Made explicit array boundaries
   ! Worked around addressing zero-th element in scratch arrays

use globals
  
private
  
  public blktri
  public tblktri
  public fishfin
  public fish_kind
  public fishworkspace
  public pois3d 

  integer, parameter :: fish_kind = 8


!
!     file fish.f
!
!
!     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
!     *                                                               *
!     *                  copyright (c) 2005 by UCAR                   *
!     *                                                               *
!     *       University Corporation for Atmospheric Research         *
!     *                                                               *
!     *                      all rights reserved                      *
!     *                                                               *
!     *                    FISHPACK90  version 1.1                    *
!     *                                                               *
!     *                 A Package of Fortran 77 and 90                *
!     *                                                               *
!     *                Subroutines and Example Programs               *
!     *                                                               *
!     *               for Modeling Geophysical Processes              *
!     *                                                               *
!     *                             by                                *
!     *                                                               *
!     *        John Adams, Paul Swarztrauber and Roland Sweet         *
!     *                                                               *
!     *                             of                                *
!     *                                                               *
!     *         the National Center for Atmospheric Research          *
!     *                                                               *
!     *                Boulder, Colorado  (80307)  U.S.A.             *
!     *                                                               *
!     *                   which is sponsored by                       *
!     *                                                               *
!     *              the National Science Foundation                  *
!     *                                                               *
!     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
!
!

!     this module is used by all fishpack solvers to allocate
!     real and complex work space
!      MODULE fish
        TYPE fishworkspace
          private
          REAL(fish_kind),POINTER,DIMENSION(:) :: rew => NULL()
          COMPLEX,POINTER,DIMENSION(:) :: cxw => NULL()
          INTEGER :: nrew
          INTEGER :: ncxw
        END TYPE fishworkspace


  ! Used to be a common block
  !    COMMON /CBLKT/ NPP, K, EPS, CNV, NM, NCMPLX, IK
  INTEGER  :: NPP, K, NM, NCMPLX, IK
  REAL(fish_kind) ::   EPS, CNV

 CONTAINS
    
        SUBROUTINE allocatfish(irwk,icwk,wsave,ierror)
        IMPLICIT NONE
        TYPE (fishworkspace) :: wsave
!       irwk is the required real work space length
!       icwk is the required integer work space length
        INTEGER, INTENT(IN) :: irwk,icwk
!       ierror is set to 20 if the dynamic allocation is unsuccessful
!       (e.g., this would happen if m,n are too large for the computers memory
        INTEGER, INTENT(INOUT) :: ierror
        INTEGER :: istatus
!       first deallocate to avoid memory leakage
#ifndef G95
#ifndef GFORTRAN
        if(associated(wsave%rew)) then 
         DEALLOCATE(wsave%rew)
         wsave%rew => NULL()
         wsave%nrew = 0
   endif
        if(associated(wsave%cxw)) then 
      DEALLOCATE(wsave%cxw)
      wsave%cxw => NULL()
      wsave%ncxw = 0
   endif
#endif
#endif
!       allocate irwk words of real work space
        if (irwk > 0) then
             ALLOCATE(wsave%rew(irwk),STAT = istatus)
             wsave%nrew = irwk
             wsave%rew(1:irwk) = 0
        end if
!       allocate icwk words of complex work space
        if (icwk > 0) then
             ALLOCATE(wsave%cxw(icwk),STAT = istatus)
             wsave%ncxw = icwk
             wsave%cxw(1:icwk) = 0
        end if
        ierror = 0
!       flag fatal error if allocation fails
!       IF (istatus /= 0) THEN
        if (istatus .ne. 0 ) then
          ierror = 20
        END IF
        RETURN
        END SUBROUTINE allocatfish

        SUBROUTINE BLK_space(N,M,irwk,icwk)
!       this subroutine computes the real and complex work space
!       requirements (generous estimate) of blktri for N,M values
        IMPLICIT NONE
        INTEGER,INTENT(IN) :: N,M
        INTEGER,INTENT(OUT) :: irwk,icwk
        INTEGER :: L,log2n
!       compute nearest integer greater than or equal to
!       log base 2 of n+1, i.e., log2n is smallest integer
!       such that 2**log2n >= n+1
        log2n = 1
        do
           log2n = log2n+1
           if (n+1 <= 2**log2n) EXIT
        end do
        L = 2**(log2n+1)
!                    !From  BLKTRI manual for fishpack 4.1
!                          IF NP=1 DEFINE K=INT(LOG2(N))+1 AND          
!                          SET L=2**(K+1) THEN W MUST HAVE DIMENSION    
!                          (K-2)*L+K+5+MAX(2N,6M)                       
!                                                                       
!                          IF NP=0 DEFINE K=INT(LOG2(N-1))+1 AND        
!                          SET L=2**(K+1) THEN W MUST HAVE DIMENSION    
!                          (K-2)*L+K+5+2N+MAX(2N,6M)                    

            
        irwk = (log2n-2)*L+5+MAX0(2*N,6*M)+log2n+2*n
        icwk = ((log2n-2)*L+5+log2n)/2+3*M+N
!        write(*,*) "N, M, irwk, icwk", N, M, irwk, icwk
        icwk = irwk  ! For some reason irwk is not enough....
        RETURN
        END SUBROUTINE BLK_space

        SUBROUTINE GEN_space(N,M,irwk)
!       this subroutine computes the real work space
!       requirement (generously) of genbun for the current N,M
        IMPLICIT NONE
        INTEGER,INTENT(IN) :: N,M
        INTEGER,INTENT(OUT) :: irwk
        INTEGER :: log2n
!       compute nearest integer greater than or equal to
!       log base 2 of n+1, i.e., log2n is smallest integer
!       such that 2**log2n >= n+1
        log2n = 1
        do
           log2n = log2n+1
           if (n+1 <= 2**log2n) EXIT
        end do
        irwk = 4*N + (10 + log2n)*M
        RETURN
        END SUBROUTINE GEN_space

        SUBROUTINE fishfin(wsave)
!       this subroutine releases allocated work space
!       fishfin should be called after a fishpack solver has finished
!       TYPE (fishworkspace) variable wsave.
        IMPLICIT NONE
        TYPE (fishworkspace) :: wsave
        INTEGER :: istatus
#ifndef G95
#ifndef GFORTRAN
        if(associated(wsave%rew))DEALLOCATE(wsave%rew)
        if(associated(wsave%cxw))DEALLOCATE(wsave%cxw)
#endif
#endif
        RETURN
        END SUBROUTINE fishfin

!      END MODULE fish
!
!     file comf.f
!
!
!     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
!     *                                                               *
!     *                  copyright (c) 2005 by UCAR                   *
!     *                                                               *
!     *       University Corporation for Atmospheric Research         *
!     *                                                               *
!     *                      all rights reserved                      *
!     *                                                               *
!     *                    FISHPACK90  version 1.1                    *
!     *                                                               *
!     *                 A Package of Fortran 77 and 90                *
!     *                                                               *
!     *                Subroutines and Example Programs               *
!     *                                                               *
!     *               for Modeling Geophysical Processes              *
!     *                                                               *
!     *                             by                                *
!     *                                                               *
!     *        John Adams, Paul Swarztrauber and Roland Sweet         *
!     *                                                               *
!     *                             of                                *
!     *                                                               *
!     *         the National Center for Atmospheric Research          *
!     *                                                               *
!     *                Boulder, Colorado  (80307)  U.S.A.             *
!     *                                                               *
!     *                   which is sponsored by                       *
!     *                                                               *
!     *              the National Science Foundation                  *
!     *                                                               *
!     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
!
!
! PACKAGE COMF           THE ENTRIES IN THIS PACKAGE ARE LOWLEVEL
!                        ENTRIES, SUPPORTING FISHPACK ENTRIES BLKTRI
!                        AND CBLKTRI. THAT IS, THESE ROUTINES ARE
!                        NOT CALLED DIRECTLY BY USERS, BUT RATHER
!                        BY ENTRIES WITHIN BLKTRI AND CBLKTRI.
!                        DESCRIPTION OF ENTRIES EPMACH
!                        FOLLOW BELOW.
!
! LATEST REVISION        JUNE 2004
!
! SPECIAL CONDITIONS     NONE
!
! I/O                    NONE
!
! PRECISION              SINGLE
!
! REQUIRED LIBRARY       NONE
! FILES
!
! LANGUAGE               FORTRAN 90
! ********************************************************************
!
! FUNCTION EPMACH (DUM)
!
! PURPOSE                TO COMPUTE AN APPROXIMATE MACHINE ACCURACY
!                        EPSILON ACCORDING TO THE FOLLOWING DEFINITION:
!                        EPSILON IS THE SMALLEST NUMBER SUCH THAT
!                        (1.+EPSILON).GT.1.)
!
! USAGE                  EPS = EPMACH (DUM)
!
! ARGUMENTS
! ON INPUT               DUM
!                          DUMMY VALUE
!
! ARGUMENTS
! ON OUTPUT              NONE
!
! HISTORY                THE ORIGINAL VERSION, WRITTEN WHEN THE
!                        BLKTRI PACKAGE WAS CONVERTED FROM THE
!                        CDC 7600 TO RUN ON THE CRAY-1, CALCULATED
!                        MACHINE ACCURACY BY SUCCESSIVE DIVISIONS
!                        BY 10.  USE OF THIS CONSTANT CAUSED BLKTRI
!                        TO COMPUTE SOLUTIONS ON THE CRAY-1 WITH FOUR
!                        FEWER PLACES OF ACCURACY THAN THE VERSION
!                        ON THE 7600.  IT WAS FOUND THAT COMPUTING
!                        MACHINE ACCURACY BY SUCCESSIVE DIVISIONS
!                        OF 2 PRODUCED A MACHINE ACCURACY 29% LESS
!                        THAN THE VALUE GENERATED BY SUCCESSIVE
!                        DIVISIONS BY 10, AND THAT USE OF THIS
!                        MACHINE CONSTANT IN THE BLKTRI PACKAGE
!                        RECOVERED THE ACCURACY THAT APPEARED TO
!                        BE LOST ON CONVERSION.
!
! ALGORITHM              COMPUTES MACHINE ACCURACY BY SUCCESSIVE
!                        DIVISIONS OF TWO.
!
! PORTABILITY            THIS CODE WILL EXECUTE ON MACHINES OTHER
!                        THAN THE CRAY1, BUT THE RETURNED VALUE MAY
!                        BE UNSATISFACTORY.  SEE HISTORY ABOVE.
! Note from M.Sofiev, 2015.07.08.
! The original algorithm bluntly fails in an optimising compiler,
! which gives EPS=0.0 by optimising V out and rewriting the while condition
! as EPS > 0. A simpler and robust approach is not to use V at all, to take the 
! constants of the correct kind, and to declare EPS as a result value.
!
!
! ********************************************************************
!
      REAL(fish_kind) FUNCTION EPMACH () 
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      REAL(fish_kind)  :: DUM 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
!      REAL(fish_kind) :: EPMACH,V
      real(fish_kind), parameter :: one = 1.0, two = 2.0
      
      EPMACH = one
      do while(EPMACH + one > one)
        EPMACH = EPMACH / two
      end do
!      EPMACH = epsilon(one)  ??????????????
!      V = EPS + 1. 
!      DO WHILE(V - 1. > 0.) 
!         EPS = EPS/2. 
!         V =  EPS + 1. 
!      END DO 

      EPMACH = 100.*EPMACH 
      RETURN  

      END FUNCTION EPMACH 



      REAL(fish_kind) FUNCTION PPSGF (X, IZ, C, A, BH) 
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER , INTENT(IN) :: IZ 
      REAL(fish_kind) , INTENT(IN) :: X 
      REAL(fish_kind)  :: C(*) 
      REAL(fish_kind)  :: A(*) 
      REAL(fish_kind) , INTENT(IN) :: BH(*) 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: J 
      REAL(fish_kind) :: ALL, SUM 

      SAVE ALL 
!-----------------------------------------------
      SUM = 0. 
      DO J = 1, IZ 
         SUM = SUM - 1./(X - BH(J))**2 
      END DO 
      PPSGF = SUM 
      RETURN  
      END FUNCTION PPSGF 


      REAL(fish_kind) FUNCTION PPSPF (X, IZ, C, A, BH) 
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER , INTENT(IN) :: IZ 
      REAL(fish_kind) , INTENT(IN) :: X 
      REAL(fish_kind)  :: C(*) 
      REAL(fish_kind)  :: A(*) 
      REAL(fish_kind) , INTENT(IN) :: BH(*) 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: J 
      REAL(fish_kind) :: ALL, SUM 

      SAVE ALL 
!-----------------------------------------------
      SUM = 0. 
      DO J = 1, IZ 
         SUM = SUM + 1./(X - BH(J)) 
      END DO 
      PPSPF = SUM 
      RETURN  
      END FUNCTION PPSPF 


      REAL(fish_kind) FUNCTION PSGF (X, IZ, C, A, BH) 
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER , INTENT(IN) :: IZ 
      REAL(fish_kind) , INTENT(IN) :: X 
      REAL(fish_kind) , INTENT(IN) :: C(*) 
      REAL(fish_kind) , INTENT(IN) :: A(*) 
      REAL(fish_kind) , INTENT(IN) :: BH(*) 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: J 
      REAL(fish_kind) :: ALL, FSG, HSG, DD 

      SAVE ALL 
!-----------------------------------------------
      FSG = 1. 
      HSG = 1. 
      DO J = 1, IZ 
         DD = 1./(X - BH(J)) 
         FSG = FSG*A(J)*DD 
         HSG = HSG*C(J)*DD 
      END DO 
      IF (MOD(IZ,2) == 0) THEN 
         PSGF = 1. - FSG - HSG 
         RETURN  
      ENDIF 
      PSGF = 1. + FSG + HSG 
      RETURN  
!
! REVISION HISTORY---
!
! SEPTEMBER 1973    VERSION 1
! APRIL     1976    VERSION 2
! JANUARY   1978    VERSION 3
! DECEMBER  1979    VERSION 3.1
! FEBRUARY  1985    DOCUMENTATION UPGRADE
! NOVEMBER  1988    VERSION 3.2, FORTRAN 77 CHANGES
! June      2004    Version 5.0, Fortran 90 changes
!-----------------------------------------------------------------------
      END FUNCTION PSGF 
!
!     file blktri.f
!
!     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
!     *                                                               *
!     *                  copyright (c) 2005 by UCAR                   *
!     *                                                               *
!     *       University Corporation for Atmospheric Research         *
!     *                                                               *
!     *                      all rights reserved                      *
!     *                                                               *
!     *                    FISHPACK90  version 1.1                    *
!     *                                                               *
!     *                 A Package of Fortran 77 and 90                *
!     *                                                               *
!     *                Subroutines and Example Programs               *
!     *                                                               *
!     *               for Modeling Geophysical Processes              *
!     *                                                               *
!     *                             by                                *
!     *                                                               *
!     *        John Adams, Paul Swarztrauber and Roland Sweet         *
!     *                                                               *
!     *                             of                                *
!     *                                                               *
!     *         the National Center for Atmospheric Research          *
!     *                                                               *
!     *                Boulder, Colorado  (80307)  U.S.A.             *
!     *                                                               *
!     *                   which is sponsored by                       *
!     *                                                               *
!     *              the National Science Foundation                  *
!     *                                                               *
!     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
!
!
!     SUBROUTINE BLKTRI (IFLG,NP,N,AN,BN,CN,MP,M,AM,BM,CM,IDIMY,Y,
!    +                   IERROR,W)
!
!
!                                                                       
! DIMENSION OF           AN(N),BN(N),CN(N),AM(M),BM(M),CM(M),Y(IDIMY,N),
! ARGUMENTS
!                                                                       
! LATEST REVISION        JUNE 2004
!                                                                       
! USAGE                  CALL BLKTRI (IFLG,NP,N,AN,BN,CN,MP,M,AM,BM,    
!                                     CM,IDIMY,Y,IERROR,W)              
!                                                                       
! PURPOSE                BLKTRI SOLVES A SYSTEM OF LINEAR EQUATIONS     
!                        OF THE FORM                                    
!                                                                       
!                        AN(J)*X(I,J-1) + AM(I)*X(I-1,J) +              
!                        (BN(J)+BM(I))*X(I,J) + CN(J)*X(I,J+1) +        
!                        CM(I)*X(I+1,J) = Y(I,J)                        
!                                                                       
!                        FOR I = 1,2,...,M  AND  J = 1,2,...,N.         
!                                                                       
!                        I+1 AND I-1 ARE EVALUATED MODULO M AND         
!                        J+1 AND J-1 MODULO N, I.E.,                    
!                                                                       
!                        X(I,0) = X(I,N),  X(I,N+1) = X(I,1),           
!                        X(0,J) = X(M,J),  X(M+1,J) = X(1,J).           
!                                                                       
!                        THESE EQUATIONS USUALLY RESULT FROM THE        
!                        DISCRETIZATION OF SEPARABLE ELLIPTIC           
!                        EQUATIONS.  BOUNDARY CONDITIONS MAY BE         
!                        DIRICHLET, NEUMANN, OR PERIODIC.               
!                                                                       
! ARGUMENTS                                                             
!                                                                       
! ON INPUT               IFLG                                           
!                                                                       
!                          = 0  INITIALIZATION ONLY.                    
!                               CERTAIN QUANTITIES THAT DEPEND ON NP,   
!                               N, AN, BN, AND CN ARE COMPUTED AND      
!                               STORED IN DERIVED data type w (see
!                               description of w below)
!                                                                       
!                          = 1  THE QUANTITIES THAT WERE COMPUTED       
!                               IN THE INITIALIZATION ARE USED          
!                               TO OBTAIN THE SOLUTION X(I,J).          
!                                                                       
!                               NOTE:                                   
!                               A CALL WITH IFLG=0 TAKES                
!                               APPROXIMATELY ONE HALF THE TIME         
!                               AS A CALL WITH IFLG = 1.                
!                               HOWEVER, THE INITIALIZATION DOES        
!                               NOT HAVE TO BE REPEATED UNLESS NP,      
!                               N, AN, BN, OR CN CHANGE.                
!                                                                       
!                        NP                                             
!                          = 0  IF AN(1) AND CN(N) ARE NOT ZERO,        
!                               WHICH CORRESPONDS TO PERIODIC           
!                               BOUNARY CONDITIONS.                     
!                                                                       
!                          = 1  IF AN(1) AND CN(N) ARE ZERO.            
!                                                                       
!                        N                                              
!                          THE NUMBER OF UNKNOWNS IN THE J-DIRECTION.   
!                          N MUST BE GREATER THAN 4.                    
!                          THE OPERATION COUNT IS PROPORTIONAL TO       
!                          MNLOG2(N), HENCE N SHOULD BE SELECTED        
!                          LESS THAN OR EQUAL TO M.                     
!                                                                       
!                        AN,BN,CN                                       
!                          ONE-DIMENSIONAL ARRAYS OF LENGTH N           
!                          THAT SPECIFY THE COEFFICIENTS IN THE         
!                          LINEAR EQUATIONS GIVEN ABOVE.                
!                                                                       
!                        MP                                             
!                          = 0  IF AM(1) AND CM(M) ARE NOT ZERO,        
!                               WHICH CORRESPONDS TO PERIODIC           
!                               BOUNDARY CONDITIONS.                    
!                                                                       
!                          = 1  IF AM(1) = CM(M) = 0  .                 
!                                                                       
!                        M                                              
!                          THE NUMBER OF UNKNOWNS IN THE I-DIRECTION.   
!                           M MUST BE GREATER THAN 4.                   
!                                                                       
!                        AM,BM,CM                                       
!                          ONE-DIMENSIONAL ARRAYS OF LENGTH M THAT      
!                          SPECIFY THE COEFFICIENTS IN THE LINEAR       
!                          EQUATIONS GIVEN ABOVE.                       
!                                                                       
!                        IDIMY                                          
!                          THE ROW (OR FIRST) DIMENSION OF THE          
!                          TWO-DIMENSIONAL ARRAY Y AS IT APPEARS        
!                          IN THE PROGRAM CALLING BLKTRI.               
!                          THIS PARAMETER IS USED TO SPECIFY THE        
!                          VARIABLE DIMENSION OF Y.                     
!                          IDIMY MUST BE AT LEAST M.                    
!                                                                       
!                        Y                                              
!                          A TWO-DIMENSIONAL ARRAY THAT SPECIFIES       
!                          THE VALUES OF THE RIGHT SIDE OF THE LINEAR   
!                          SYSTEM OF EQUATIONS GIVEN ABOVE.             
!                          Y MUST BE DIMENSIONED AT LEAST M*N.          
!                                                                       
!                        W
!                          A fortran 90 derived TYPE (fishworkspace) variable
!                          that must be declared by the user.  The first
!                          two declarative statements in the user program
!                          calling BLKTRI must be:
!
!                               USE fish
!                               TYPE (fishworkspace) :: W
!
!                          The first statement makes the fishpack module
!                          defined in the file "fish.f" available to the
!                          user program calling BLKTRI.  The second statement
!                          declares a derived type variable (defined in
!                          the module "fish.f") which is used internally
!                          in BLKTRI to dynamically allocate real and complex
!                          work space used in solution.  An error flag
!                          (IERROR = 20) is set if the required work space
!                          allocation fails (for example if N,M are too large)
!                          Real and complex values are set in the components
!                          of W on a initial (IFLG=0) call to BLKTRI.  These
!                          must be preserved on non-initial calls (IFLG=1)
!                          to BLKTRI.  This eliminates redundant calculations
!                          and saves compute time.
!               ****       IMPORTANT!  The user program calling BLKTRI should
!                          include the statement:
!
!                               CALL FISHFIN(W)
!
!                          after the final approximation is generated by
!                          BLKTRI.  The will deallocate the real and complex
!                          work space of W.  Failure to include this statement
!                          could result in serious memory leakage.
!
!                                                                       
! ARGUMENTS                                                             
!                                                                       
! ON OUTPUT              Y                                              
!                          CONTAINS THE SOLUTION X.                     
!                                                                       
!                        IERROR                                         
!                          AN ERROR FLAG THAT INDICATES INVALID         
!                          INPUT PARAMETERS.  EXCEPT FOR NUMBER ZER0,   
!                          A SOLUTION IS NOT ATTEMPTED.                 
!                                                                       
!                        = 0  NO ERROR.                                 
!                        = 1  M IS LESS THAN 5                          
!                        = 2  N IS LESS THAN 5                          
!                        = 3  IDIMY IS LESS THAN M.                     
!                        = 4  BLKTRI FAILED WHILE COMPUTING RESULTS     
!                             THAT DEPEND ON THE COEFFICIENT ARRAYS     
!                             AN, BN, CN.  CHECK THESE ARRAYS.          
!                        = 5  AN(J)*CN(J-1) IS LESS THAN 0 FOR SOME J.  
!                                                                       
!                             POSSIBLE REASONS FOR THIS CONDITION ARE   
!                             1. THE ARRAYS AN AND CN ARE NOT CORRECT   
!                             2. TOO LARGE A GRID SPACING WAS USED      
!                                IN THE DISCRETIZATION OF THE ELLIPTIC  
!                                EQUATION.                              
!                             3. THE LINEAR EQUATIONS RESULTED FROM A   
!                                PARTIAL DIFFERENTIAL EQUATION WHICH    
!                                WAS NOT ELLIPTIC.                      
!                                                                       
!                        = 20 If the dynamic allocation of real and
!                             complex work space in the derived type
!                             (fishworkspace) variable W fails (e.g.,
!                             if N,M are too large for the platform used)
!
!                                                                       
!                        W
!                             The derived type (fishworkspace) variable W
!                             contains real and complex values that must not
!                             be destroyed if BLKTRI is called again with
!                             IFLG=1.
!                                                                       
!                                                                       
! SPECIAL CONDITIONS     THE ALGORITHM MAY FAIL IF ABS(BM(I)+BN(J))     
!                        IS LESS THAN ABS(AM(I))+ABS(AN(J))+            
!                        ABS(CM(I))+ABS(CN(J))                          
!                        FOR SOME I AND J. THE ALGORITHM WILL ALSO      
!                        FAIL IF AN(J)*CN(J-1) IS LESS THAN ZERO FOR    
!                        SOME J.                                        
!                        SEE THE DESCRIPTION OF THE OUTPUT PARAMETER    
!                        IERROR.                                        
!                                                                       
! I/O                    NONE                                           
!                                                                       
! PRECISION              SINGLE                                         
!                                                                       
! REQUIRED FILES         fish.f,comf.f
!                                                                       
! LANGUAGE               FORTRAN 90
!                                                                       
! HISTORY                WRITTEN BY PAUL SWARZTRAUBER AT NCAR IN THE    
!                        EARLY 1970'S.  REWRITTEN AND RELEASED IN       
!                        LIBRARIES IN JANUARY 1980. Revised in June
!                        2004 using Fortan 90 dynamically allocated work
!                        space and derived data types to eliminate mixed
!                        mode conflicts in the earlier versions.
!                                                                       
! ALGORITHM              GENERALIZED CYCLIC REDUCTION                   
!                                                                       
! PORTABILITY            FORTRAN 90.  APPROXIMATE MACHINE ACCURACY
!                        IS COMPUTED IN FUNCTION EPMACH.                
!                                                                       
! REFERENCES             SWARZTRAUBER,P. AND R. SWEET, 'EFFICIENT       
!                        FORTRAN SUBPROGRAMS FOR THE SOLUTION OF        
!                        ELLIPTIC EQUATIONS'                            
!                        NCAR TN/IA-109, JULY, 1975, 138 PP.            
!                                                                       
!                        SWARZTRAUBER P. N.,A DIRECT METHOD FOR         
!                        THE DISCRETE SOLUTION OF SEPARABLE             
!                        ELLIPTIC EQUATIONS, S.I.A.M.                   
!                        J. NUMER. ANAL.,11(1974) PP. 1136-1150.        
!***********************************************************************
      SUBROUTINE BLKTRI(IFLG, NP, N, AN, BN, CN, MP, M, AM, BM, CM, IDIMY, Y, IERROR, W)
!      USE fish
      Implicit none
      TYPE (fishworkspace) :: w
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER  :: IFLG
      INTEGER  :: NP
      INTEGER  :: N
      INTEGER  :: MP
      INTEGER  :: M
      INTEGER  :: IDIMY
      INTEGER  :: IERROR
      REAL(fish_kind)  :: AN(*)
      REAL(fish_kind)  :: BN(*)
      REAL(fish_kind)  :: CN(*)
      REAL(fish_kind)  :: AM(*)
      REAL(fish_kind)  :: BM(*)
      REAL(fish_kind)  :: CM(*)
      REAL(fish_kind)  :: Y(IDIMY,*)
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: IRWK, ICWK
!-----------------------------------------------
      IF (M < 5) THEN
         IERROR = 1
         RETURN 
      ENDIF
      IF (N < 3) THEN
         IERROR = 2
         RETURN 
      ENDIF
      IF (IDIMY < M) THEN
         IERROR = 3
         RETURN 
      ENDIF
      IF (IFLG == 0) THEN
!     compute and allocate real and complex required work space
         CALL BLK_SPACE (N, M, IRWK, ICWK)
         CALL ALLOCATFISH (IRWK, ICWK, W, IERROR)
         IF (IERROR == 20) RETURN 
      ENDIF
      call blktrii (IFLG,NP,N,AN,BN,CN,MP,M,AM,BM,CM,IDIMY,Y, IERROR,  w%nrew,  w%rew,  w%ncxw,  w%cxw)
      RETURN 
      END SUBROUTINE BLKTRI


 
      SUBROUTINE BLKTRII(IFLG, NP, N, AN, BN, CN, MP, M, AM, BM, CM, IDIMY, Y, IERROR, NW, W, NWC, WC)
      Implicit none
!      external prod,prodp,cprod,cprodp
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER , INTENT(IN) :: IFLG
      INTEGER , INTENT(IN) :: NP
      INTEGER , INTENT(IN) :: N
      INTEGER , INTENT(IN) :: MP
      INTEGER , INTENT(IN) :: M
      INTEGER , INTENT(IN) :: IDIMY
      INTEGER , INTENT(IN) :: NW
      INTEGER , INTENT(IN) :: NWC
      INTEGER  :: IERROR
      REAL(fish_kind)  :: AN(N)
      REAL(fish_kind)  :: BN(N)
      REAL(fish_kind)  :: CN(N)
      REAL(fish_kind)  :: AM(M)
      REAL(fish_kind)  :: BM(M)
      REAL(fish_kind)  :: CM(M)
      REAL(fish_kind)  :: Y(IDIMY,N)
      REAL(fish_kind)  :: W(NW)
      COMPLEX  :: WC(NWC)
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: IW1, IW2, IW3, IWW, IWU, IWD, NL, NH, IWAH, IWBH

      SAVE IW1, IW2, IW3, IWW, IWU, IWD, NL
!-----------------------------------------------
 
!
! TEST M AND N FOR THE PROPER FORM
!
      NM = N
!     check again for solvers which call blktrii directly
      IF (M < 5) THEN
         IERROR = 1
         RETURN 
      ENDIF
      IF (NM < 3) THEN
         IERROR = 2
         RETURN 
      ENDIF
      IF (IDIMY < M) THEN
         IERROR = 3
         RETURN 
      ENDIF
 
      IF (IFLG == 0) THEN
         NH = N
         NPP = NP
         IF (NPP /= 0) THEN
            NH = NH + 1
         ENDIF
         IK = 2
         K = 1
         IK = IK + IK
         K = K + 1
         DO WHILE(NH - IK > 0)
            IK = IK + IK
            K = K + 1
         END DO
         NL = IK
         IK = IK + IK
         NL = NL - 1
         IWAH = (K - 2)*IK + K + 5
         IF (NPP == 0) THEN
            IWBH = IWAH + NM + NM
            IW1 = IWBH
            NM = N - 1
         ELSE
            IW1 = IWAH
            IWBH = IW1 + NM
         ENDIF
!     set pointers in real,complex work space
         IW2 = IW1 + M
         IW3 = IW2 + M
         IWD = IW3 + M
         IWW = IWD + M
         IWU = IWW + M
         CALL COMPB (NL, IERROR, AN, BN, CN, W, WC, W(IWAH), W(IWBH))
         RETURN 
      ENDIF
! *** important to reset nm for np = 0
      IF (NPP == 0) NM = N - 1

!      write (*,*) "Before BLKTR1: Common block content:", NPP, K, 
!     1 EPS, CNV, NM, NCMPLX, IK
!      write (*,*) "BLKTR1 called with N", NL
!      write (*,*) "BLKTR1 called with AN", AN
!      write (*,*) "BLKTR1 called with BN", BN
!      write (*,*) "BLKTR1 called with CN", CN 
!      write (*,*) "BLKTR1 called with M",  M 
!      write (*,*) "BLKTR1 called with AM", AM
!      write (*,*) "BLKTR1 called with BM", BM
!      write (*,*) "BLKTR1 called with CM", CM
!      write (*,*) "BLKTR1 called with IDIMY", IDIMY
!      write (*,*) "BLKTR1 called with NB", NW
!      write (*,*) "BLKTR1 called with B", W
!      write (*,*) "BLKTR1 called with NBC", NWC
!      write (*,*) "BLKTR1 called with BC", WC
!      write (*,*) "BLKTR1 called with  IW1,IW2, IW3",  IW1,IW2, IW3
!      write (*,*) "BLKTR1 called with  IWD,IWW, IWU",  IWD,IWW, IWU
      
 
      IF (MP /= 0) THEN
         CALL BLKTR1 (NL, AN, BN, CN, M, AM, BM, CM, IDIMY, Y, NW, W, &
        &      NWC, WC, W(IW1), W(IW2), W(IW3), W(IWD), W(IWW), W(IWU),& 
        &      WC(IW1), WC(IW2), WC(IW3), PROD, CPROD)
      ELSE
         CALL BLKTR1 (NL, AN, BN, CN, M, AM, BM, CM, IDIMY, Y, NW, W, &
        &      NWC, WC, W(IW1), W(IW2), W(IW3), W(IWD), W(IWW), W(IWU),& 
        &      WC(IW1), WC(IW2), WC(IW3), PRODP, CPRODP)
      ENDIF
      RETURN 
      END SUBROUTINE BLKTRII


      SUBROUTINE BLKTR1(N, AN, BN, CN, M, AM, BM, CM, IDIMY, Y, NB, B, &
        &   NBC, BC, W1, W2, W3, WD, WW, WU, CW1, CW2, CW3, PRDCT, CPRDCT)
      Implicit none
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER , INTENT(IN) :: N
      INTEGER , INTENT(IN) :: M
      INTEGER , INTENT(IN) :: IDIMY
      INTEGER , INTENT(IN) :: NB
      INTEGER , INTENT(IN) :: NBC
      REAL(fish_kind)  :: AN(N)
      REAL(fish_kind)  :: BN(N)
      REAL(fish_kind)  :: CN(N)
      REAL(fish_kind)  :: AM(M)
      REAL(fish_kind)  :: BM(M)
      REAL(fish_kind)  :: CM(M)
      REAL(fish_kind)  :: Y(IDIMY,N)
      REAL(fish_kind)  :: B(NB)
      REAL(fish_kind)  :: W1(*)
      REAL(fish_kind)  :: W2(*)
      REAL(fish_kind)  :: W3(*)
      REAL(fish_kind)  :: WD(*)
      REAL(fish_kind)  :: WW(*)
      REAL(fish_kind)  :: WU(*)
      COMPLEX  :: BC(NBC)
      COMPLEX  :: CW1(*)
      COMPLEX  :: CW2(*)
      COMPLEX  :: CW3(*)
      external PRDCT, CPRDCT
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: KDO, L, IR, I2, I1, I3, I4, IRM1, IM2, NM2, IM3, NM3 
      INTEGER :: IM1, NM1, I0, IF, I, IPI1, IPI2, IPI3, IDXC, NC, IDXA, NA, IP2
      INTEGER :: NP2, IP1, NP1, IP3, NP3, J, IZ, NZ, IZR, LL, IFD, IP, NP
      INTEGER :: IMI1, IMI2
      REAL(fish_kind) :: DUM
!-----------------------------------------------
!
! BLKTR1 SOLVES THE LINEAR SYSTEM
!
! B  CONTAINS THE ROOTS OF ALL THE B POLYNOMIALS
! W1,W2,W3,WD,WW,WU  ARE ALL WORKING ARRAYS
! PRDCT  IS EITHER PRODP OR PROD DEPENDING ON WHETHER THE BOUNDARY
! CONDITIONS IN THE M DIRECTION ARE PERIODIC OR NOT
! CPRDCT IS EITHER CPRODP OR CPROD WHICH ARE THE COMPLEX VERSIONS
! OF PRODP AND PROD. THESE ARE CALLED IN THE EVENT THAT SOME
! OF THE ROOTS OF THE B SUB P POLYNOMIAL ARE COMPLEX
!
!
!
! BEGIN REDUCTION PHASE
!
!      write (*,*) "Entering BLKTR1: Common block content:", NPP, K, 
!     1 EPS, CNV, NM, NCMPLX, IK
!      write (*,*) "BLKTR1 called with N", N
!      write (*,*) "BLKTR1 called with AN", AN
!      write (*,*) "BLKTR1 called with BN", BN
!      write (*,*) "BLKTR1 called with CN", CN 
!      write (*,*) "BLKTR1 called with M",  M 
!      write (*,*) "BLKTR1 called with AM", AM
!      write (*,*) "BLKTR1 called with BM", BM
!      write (*,*) "BLKTR1 called with CM", CM
!      write (*,*) "BLKTR1 called with IDIMY", IDIMY
!      write (*,*) "BLKTR1 called with NB",  NB
!      write (*,*) "BLKTR1 called with B",  B
!      write (*,*) "BLKTR1 called with NBC", NBC
!      write (*,*) "BLKTR1 called with W1", W1
!      write (*,*) "BLKTR1 called with W2", W2
!      write (*,*) "BLKTR1 called with W3", W3
!      write (*,*) "BLKTR1 called with WD", WD
!      write (*,*) "BLKTR1 called with WW", WW
!      write (*,*) "BLKTR1 called with WU", WU
!      write (*,*) "BLKTR1 called with CW1", CW1
!      write (*,*) "BLKTR1 called with CW2", CW2
!      write (*,*) "BLKTR1 called with CW3", CW3
      KDO = K - 1
      DO L = 1, KDO
         IR = L - 1
         I2 = 2**IR
         I1 = I2/2
         I3 = I2 + I1
         I4 = I2 + I2
         IRM1 = IR - 1
!         write (*,*) "Before 1 INDXB: Common block content:", NPP, K, 
!     1       EPS, CNV, NM, NCMPLX, IK
!          write (*,*) "Calling idxb I2, IR, IM2, NM2",I2, IR, IM2, NM2
         CALL INDXB (I2, IR, IM2, NM2)
!         write (*,*) "Before 2 INDXB: Common block content:", NPP, K, 
!     1       EPS, CNV, NM, NCMPLX, IK
!         write (*,*) "Calling idxb I1, IRM1, IM3, NM3",I1, IRM1,IM3,NM3
         CALL INDXB (I1, IRM1, IM3, NM3)
!         write (*,*) "Before 3 INDXB: Common block content:", NPP, K, 
!     1       EPS, CNV, NM, NCMPLX, IK
!         write (*,*) "Calling idxb I3, IRM1,IM1,NM1",I3, IRM1, IM1, NM1
         CALL INDXB (I3, IRM1, IM1, NM1)
!         write (*,*) "After 3 INDXB: Common block content:", NPP, K, 
!     1       EPS, CNV, NM, NCMPLX, IK
         I0 = 0
         if (im1==0) IM1=1
         if (im2==0) IM2=1
         if (im3==0) IM3=1
!         write (*,*), "BLKTR1  IR IRM1:", L,IR,IRM1
!         write (*,*), "BLKTR1  I123:", I1,   I2,  I3
!         write (*,*), "BLKTR1 IM123:", IM1, IM2, IM3
!         write (*,*), "BLKTR1    B:", B(IM1), B(IM2), B(IM3)

         CALL PRDCT (NM2, B(IM2), NM3, B(IM3), NM1, B(IM1), I0, DUM, Y(1,I2), W3, M, AM, BM, CM, WD, WW, WU)
         IF = 2**K
         DO I = I4, IF, I4
            IF (I - NM > 0) CYCLE 
            IPI1 = I + I1
            IPI2 = I + I2
            IPI3 = I + I3
            CALL INDXC (I, IR, IDXC, NC)
            IF (I - IF >= 0) CYCLE 
            CALL INDXA (I, IR, IDXA, NA)
            CALL INDXB (I - I1, IRM1, IM1, NM1)
            CALL INDXB (IPI2, IR, IP2, NP2)
            CALL INDXB (IPI1, IRM1, IP1, NP1)
            CALL INDXB (IPI3, IRM1, IP3, NP3)
            if (im1==0) IM1=1
            CALL PRDCT (NM1, B(IM1), 0, DUM, 0, DUM, NA, AN(IDXA), W3, W1, M, AM, BM, CM, WD, WW, WU)
            IF (IPI2 - NM > 0) THEN
               W3(:M) = 0.
               W2(:M) = 0.
            ELSE
!               write (*,*) IP1, IP2, IP3
               if (ip1==0) Ip1=1
               if (ip2==0) Ip2=1
               if (ip3==0) Ip3=1

               CALL PRDCT (NP2, B(IP2), NP1, B(IP1), NP3, B(IP3), 0, DUM, Y(1,IPI2), W3, M, AM, BM, CM, WD, WW, WU)
               CALL PRDCT (NP1, B(IP1), 0, DUM, 0, DUM, NC, CN(IDXC), W3, W2, M, AM, BM, CM, WD, WW, WU)
            ENDIF
            Y(:M,I) = W1(:M) + W2(:M) + Y(:M,I)
         END DO
      END DO
      IF (NPP == 0) THEN
         IF = 2**K
         I = IF/2
         I1 = I/2
         CALL INDXB (I - I1, K - 2, IM1, NM1)
         CALL INDXB (I + I1, K - 2, IP1, NP1)
         CALL INDXB (I, K - 1, IZ, NZ)
         CALL PRDCT (NZ, B(IZ), NM1, B(IM1), NP1, B(IP1), 0, DUM, Y(1,I), W1, M, AM, BM, CM, WD, WW, WU)
         IZR = I
         W2(:M) = W1(:M)
         DO LL = 2, K
            L = K - LL + 1
            IR = L - 1
            I2 = 2**IR
            I1 = I2/2
            I = I2
            CALL INDXC (I, IR, IDXC, NC)
            CALL INDXB (I, IR, IZ, NZ)
            CALL INDXB (I - I1, IR - 1, IM1, NM1)
            CALL INDXB (I + I1, IR - 1, IP1, NP1)
            CALL PRDCT (NP1, B(IP1), 0, DUM, 0, DUM, NC, CN(IDXC), W1, W1, M, AM, BM, CM, WD, WW, WU)
            W1(:M) = Y(:M,I) + W1(:M)
            CALL PRDCT (NZ, B(IZ), NM1, B(IM1), NP1, B(IP1), 0, DUM, W1, W1, M, AM, BM, CM, WD, WW, WU)
         END DO
         L118: DO LL = 2, K
            L = K - LL + 1
            IR = L - 1
            I2 = 2**IR
            I1 = I2/2
            I4 = I2 + I2
            IFD = IF - I2
            DO I = I2, IFD, I4
               IF (I - I2 - IZR /= 0) CYCLE 
               IF (I - NM > 0) CYCLE  L118
               CALL INDXA (I, IR, IDXA, NA)
               CALL INDXB (I, IR, IZ, NZ)
               CALL INDXB (I - I1, IR - 1, IM1, NM1)
               CALL INDXB (I + I1, IR - 1, IP1, NP1)
               CALL PRDCT (NM1, B(IM1), 0, DUM, 0, DUM, NA, AN(IDXA), W2, W2, M, AM, BM, CM, WD, WW, WU)
               W2(:M) = Y(:M,I) + W2(:M)
               CALL PRDCT (NZ, B(IZ), NM1, B(IM1), NP1, B(IP1), 0, DUM, W2, W2, M, AM, BM, CM, WD, WW, WU)
               IZR = I
               IF (I - NM == 0) EXIT  L118
            END DO
         END DO L118
  119    CONTINUE
         Y(:M,NM+1) = Y(:M,NM+1) - CN(NM+1)*W1(:M) - AN(NM+1)*W2(:M)
         CALL INDXB (IF/2, K - 1, IM1, NM1)
         CALL INDXB (IF, K - 1, IP, NP)
         IF (NCMPLX /= 0) THEN
            CALL CPRDCT (NM + 1, BC(IP), NM1, B(IM1), 0, DUM, 0, DUM, Y(1,NM+1), Y(1,NM+1), M, AM, BM, CM, CW1, CW2, CW3)
         ELSE
            CALL PRDCT (NM + 1, B(IP), NM1, B(IM1), 0, DUM, 0, DUM, Y(1,NM+1), Y(1,NM+1), M, AM, BM, CM, WD, WW, WU)
         ENDIF
         W1(:M) = AN(1)*Y(:M,NM+1)
         W2(:M) = CN(NM)*Y(:M,NM+1)
         Y(:M,1) = Y(:M,1) - W1(:M)
         Y(:M,NM) = Y(:M,NM) - W2(:M)
         DO L = 1, KDO
            IR = L - 1
            I2 = 2**IR
            I4 = I2 + I2
            I1 = I2/2
            I = I4
            CALL INDXA (I, IR, IDXA, NA)
            CALL INDXB (I - I2, IR, IM2, NM2)
            CALL INDXB (I - I2 - I1, IR - 1, IM3, NM3)
            CALL INDXB (I - I1, IR - 1, IM1, NM1)
            CALL PRDCT (NM2, B(IM2), NM3, B(IM3), NM1, B(IM1), 0, DUM, W1, W1, M, AM, BM, CM, WD, WW, WU)
            CALL PRDCT (NM1, B(IM1), 0, DUM, 0, DUM, NA, AN(IDXA), W1, W1, M, AM, BM, CM, WD, WW, WU)
            Y(:M,I) = Y(:M,I) - W1(:M)
         END DO
!
         IZR = NM
         L131: DO L = 1, KDO
            IR = L - 1
            I2 = 2**IR
            I1 = I2/2
            I3 = I2 + I1
            I4 = I2 + I2
            IRM1 = IR - 1
            DO I = I4, IF, I4
               IPI1 = I + I1
               IPI2 = I + I2
               IPI3 = I + I3
               IF (IPI2 - IZR /= 0) THEN
                  IF (I - IZR /= 0) CYCLE 
                  CYCLE  L131
               ENDIF
               CALL INDXC (I, IR, IDXC, NC)
               CALL INDXB (IPI2, IR, IP2, NP2)
               CALL INDXB (IPI1, IRM1, IP1, NP1)
               CALL INDXB (IPI3, IRM1, IP3, NP3)
               CALL PRDCT (NP2, B(IP2), NP1, B(IP1), NP3, B(IP3), 0, DUM, W2, W2, M, AM, BM, CM, WD, WW, WU)
               CALL PRDCT (NP1, B(IP1), 0, DUM, 0, DUM, NC, CN(IDXC), W2, W2, M, AM, BM, CM, WD, WW, WU)
               Y(:M,I) = Y(:M,I) - W2(:M)
               IZR = I
               CYCLE  L131
            END DO
         END DO L131
      ENDIF
!
! BEGIN BACK SUBSTITUTION PHASE
!
      DO LL = 1, K
         L = K - LL + 1
         IR = L - 1
         IRM1 = IR - 1
         I2 = 2**IR
         I1 = I2/2
         I4 = I2 + I2
         IFD = IF - I2
         DO I = I2, IFD, I4
            IF (I - NM > 0) CYCLE 
            IMI1 = I - I1
            IMI2 = I - I2
            IPI1 = I + I1
            IPI2 = I + I2
            CALL INDXA (I, IR, IDXA, NA)
            CALL INDXC (I, IR, IDXC, NC)
            CALL INDXB (I, IR, IZ, NZ)
            CALL INDXB (IMI1, IRM1, IM1, NM1)
            CALL INDXB (IPI1, IRM1, IP1, NP1)
            if (im1==0) Im1=1
            if (ip1==0) Ip1=1
            IF (I - I2 <= 0) THEN
               W1(:M) = 0.
            ELSE
               CALL PRDCT (NM1, B(IM1), 0, DUM, 0, DUM, NA, AN(IDXA), Y(1,IMI2), W1, M, AM, BM, CM, WD, WW, WU)
            ENDIF
            IF (IPI2 - NM > 0) THEN
               W2(:M) = 0.
            ELSE
               CALL PRDCT (NP1, B(IP1), 0, DUM, 0, DUM, NC, CN(IDXC), Y(1,IPI2), W2, M, AM, BM, CM, WD, WW, WU)
            ENDIF
            W1(:M) = Y(:M,I) + W1(:M) + W2(:M)
            CALL PRDCT (NZ, B(IZ), NM1, B(IM1), NP1, B(IP1), 0, DUM, W1, Y(1,I), M, AM, BM, CM, WD, WW, WU)
         END DO
      END DO
      RETURN 
      END SUBROUTINE BLKTR1


      REAL(fish_kind) FUNCTION BSRH (XLL, XRR, IZ, C, A, BH, F, SGN)
      implicit none
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER  :: IZ
      REAL(fish_kind) , INTENT(IN) :: XLL
      REAL(fish_kind) , INTENT(IN) :: XRR
      REAL(fish_kind)  :: F
      REAL(fish_kind) , INTENT(IN) :: SGN
      REAL(fish_kind)  :: C(*)
      REAL(fish_kind)  :: A(*)
      REAL(fish_kind)  :: BH(*)
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      REAL(fish_kind) :: R1, XL, XR, DX, X
!-----------------------------------------------
!   E x t e r n a l   F u n c t i o n s
!-----------------------------------------------
!-----------------------------------------------
      XL = XLL
      XR = XRR
      DX = 0.5*ABS(XR - XL)
      X = 0.5*(XL + XR)
      R1 = SGN*F(X,IZ,C,A,BH)
      IF (R1 >= 0.) THEN
         IF (R1 == 0.) GO TO 105
         XR = X
      ELSE
         XL = X
      ENDIF
      DX = 0.5*DX
      DO WHILE(DX - CNV > 0.)
         X = 0.5*(XL + XR)
         R1 = SGN*F(X,IZ,C,A,BH)
         IF (R1 >= 0.) THEN
            IF (R1 == 0.) GO TO 105
            XR = X
         ELSE
            XL = X
         ENDIF
         DX = 0.5*DX
      END DO
  105 CONTINUE
      BSRH = 0.5*(XL + XR)
      RETURN 
      END FUNCTION BSRH


      SUBROUTINE COMPB(N, IERROR, AN, BN, CN, B, BC, AH, BH)
      implicit none
!      real :: epmach
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER  :: N
      INTEGER  :: IERROR
      REAL(fish_kind)  :: AN(*)
      REAL(fish_kind) , INTENT(IN) :: BN(*)
      REAL(fish_kind)  :: CN(*)
      REAL(fish_kind)  :: B(*)
      REAL(fish_kind)  :: AH(*)
      REAL(fish_kind)  :: BH(*)
      COMPLEX  :: BC(*)
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: J, IF, KDO, L, IR, I2, I4, IPL, IFD, I, IB, NB, JS, JF
      INTEGER :: LS, LH, NMP, L1, L2, J2, J1, N2M2
      REAL(fish_kind) :: DUM, BNORM, ARG, D1, D2, D3
!-----------------------------------------------
!
!     COMPB COMPUTES THE ROOTS OF THE B POLYNOMIALS USING SUBROUTINE
!     TEVLS WHICH IS A MODIFICATION THE EISPACK PROGRAM TQLRAT.
!     IERROR IS SET TO 4 IF EITHER TEVLS FAILS OR IF A(J+1)*C(J) IS
!     LESS THAN ZERO FOR SOME J.  AH,BH ARE TEMPORARY WORK ARRAYS.
!
      EPS = EPMACH()
call msg('EPS from COMPB',EPS)
      BNORM = ABS(BN(1))
      DO J = 2, NM
         BNORM = AMAX1(BNORM,ABS(BN(J)))
         ARG = AN(J)*CN(J-1)
         IF (ARG < 0.) GO TO 119
         B(J) = SIGN(SQRT(ARG),AN(J))
      END DO
      CNV = EPS*BNORM
      IF = 2**K
      KDO = K - 1
      L108: DO L = 1, KDO
         IR = L - 1
         I2 = 2**IR
         I4 = I2 + I2
         IPL = I4 - 1
         IFD = IF - I4
         DO I = I4, IFD, I4
            CALL INDXB (I, L, IB, NB)
            IF (NB <= 0) CYCLE  L108
            JS = I - IPL
            JF = JS + NB - 1
            LS = 0
            BH(:JF-JS+1) = BN(JS:JF)
            AH(:JF-JS+1) = B(JS:JF)
            CALL TEVLS (NB, BH, AH, IERROR)
            IF (IERROR /= 0) GO TO 118
            LH = IB - 1
            IF (NB > 0) THEN
               B(LH+1:NB+LH) = -BH(:NB)
               LH = NB + LH
            ENDIF
         END DO
      END DO L108
      B(:NM) = -BN(:NM)
      IF (NPP == 0) THEN
         NMP = NM + 1
         NB = NM + NMP
         DO J = 1, NB
            L1 = MOD(J - 1,NMP) + 1
            L2 = MOD(J + NM - 1,NMP) + 1
            ARG = AN(L1)*CN(L2)
            IF (ARG < 0.) GO TO 119
            BH(J) = SIGN(SQRT(ARG),(-AN(L1)))
            AH(J) = -BN(L1)
         END DO
         CALL TEVLS (NB, AH, BH, IERROR)
         IF (IERROR /= 0) GO TO 118
         CALL INDXB (IF, K - 1, J2, LH)
         CALL INDXB (IF/2, K - 1, J1, LH)
         J2 = J2 + 1
         LH = J2
         N2M2 = J2 + NM + NM - 2
  114    CONTINUE
         D1 = ABS(B(J1)-B(J2-1))
         D2 = ABS(B(J1)-B(J2))
         D3 = ABS(B(J1)-B(J2+1))
         IF (D2>=D1 .OR. D2>=D3) THEN
            B(LH) = B(J2)
            J2 = J2 + 1
            LH = LH + 1
            IF (J2 - N2M2 <= 0) GO TO 114
         ELSE
            J2 = J2 + 1
            J1 = J1 + 1
            IF (J2 - N2M2 <= 0) GO TO 114
         ENDIF
         B(LH) = B(N2M2+1)
         CALL INDXB (IF, K - 1, J1, J2)
         J2 = J1 + NMP + NMP
         CALL PPADD (NM + 1, IERROR, AN, CN, BC(J1), B(J1), B(J2))
      ENDIF
      RETURN 
  118 CONTINUE
      IERROR = 4
      RETURN 
  119 CONTINUE
      IERROR = 5
      RETURN 
      END SUBROUTINE COMPB


      SUBROUTINE CPROD(ND,BD,NM1,BM1,NM2,BM2,NA,AA,X,YY,M,A,B,C,D,W,Y)
      implicit none
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER , INTENT(IN) :: ND
      INTEGER , INTENT(IN) :: NM1
      INTEGER , INTENT(IN) :: NM2
      INTEGER , INTENT(IN) :: NA
      INTEGER , INTENT(IN) :: M
      REAL(fish_kind) , INTENT(IN) :: BM1(*)
      REAL(fish_kind) , INTENT(IN) :: BM2(*)
      REAL(fish_kind) , INTENT(IN) :: AA(*)
      REAL(fish_kind) , INTENT(IN) :: X(*)
      REAL(fish_kind) , INTENT(OUT) :: YY(*)
      REAL(fish_kind) , INTENT(IN) :: A(*)
      REAL(fish_kind) , INTENT(IN) :: B(*)
      REAL(fish_kind) , INTENT(IN) :: C(*)
      COMPLEX , INTENT(IN) :: BD(*)
      COMPLEX , INTENT(INOUT) :: D(*)
      COMPLEX , INTENT(INOUT) :: W(*)
      COMPLEX , INTENT(INOUT) :: Y(*)
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: J, MM, ID, M1, M2, IA, IFLG, K
      REAL(fish_kind) :: RT
      COMPLEX :: CRT, DEN, Y1, Y2
!-----------------------------------------------
!
! PROD APPLIES A SEQUENCE OF MATRIX OPERATIONS TO THE VECTOR X AND
! STORES THE RESULT IN YY           (COMPLEX CASE)
! AA   ARRAY CONTAINING SCALAR MULTIPLIERS OF THE VECTOR X
! ND,NM1,NM2 ARE THE LENGTHS OF THE ARRAYS BD,BM1,BM2 RESPECTIVELY
! BD,BM1,BM2 ARE ARRAYS CONTAINING ROOTS OF CERTIAN B POLYNOMIALS
! NA IS THE LENGTH OF THE ARRAY AA
! X,YY THE MATRIX OPERATIONS ARE APPLIED TO X AND THE RESULT IS YY
! A,B,C  ARE ARRAYS WHICH CONTAIN THE TRIDIAGONAL MATRIX
! M  IS THE ORDER OF THE MATRIX
! D,W,Y ARE WORKING ARRAYS
! ISGN  DETERMINES WHETHER OR NOT A CHANGE IN SIGN IS MADE
!
      DO J = 1, M
         Y(J) = CMPLX(X(J),0.)
      END DO
      MM = M - 1
      ID = ND
      M1 = NM1
      M2 = NM2
      IA = NA
  102 CONTINUE
      IFLG = 0
      IF (ID > 0) THEN
         CRT = BD(ID)
         ID = ID - 1
!
! BEGIN SOLUTION TO SYSTEM
!
         D(M) = A(M)/(B(M)-CRT)
         W(M) = Y(M)/(B(M)-CRT)
         DO J = 2, MM
            K = M - J
            DEN = B(K+1) - CRT - C(K+1)*D(K+2)
            D(K+1) = A(K+1)/DEN
            W(K+1) = (Y(K+1)-C(K+1)*W(K+2))/DEN
         END DO
         DEN = B(1) - CRT - C(1)*D(2)
         IF (CABS(DEN) /= 0.) THEN
            Y(1) = (Y(1)-C(1)*W(2))/DEN
         ELSE
            Y(1) = (1.,0.)
         ENDIF
         DO J = 2, M
            Y(J) = W(J) - D(J)*Y(J-1)
         END DO
      ENDIF
      IF (M1 <= 0) THEN
         IF (M2 <= 0) GO TO 121
         RT = BM2(M2)
         M2 = M2 - 1
      ELSE
         IF (M2 <= 0) THEN
            RT = BM1(M1)
            M1 = M1 - 1
         ELSE
            IF (ABS(BM1(M1)) - ABS(BM2(M2)) > 0.) THEN
               RT = BM1(M1)
               M1 = M1 - 1
            ELSE
               RT = BM2(M2)
               M2 = M2 - 1
            ENDIF
         ENDIF
      ENDIF
      Y1 = (B(1)-RT)*Y(1) + C(1)*Y(2)
      IF (MM - 2 >= 0) THEN
         DO J = 2, MM
            Y2 = A(J)*Y(J-1) + (B(J)-RT)*Y(J) + C(J)*Y(J+1)
            Y(J-1) = Y1
            Y1 = Y2
         END DO
      ENDIF
      Y(M) = A(M)*Y(M-1) + (B(M)-RT)*Y(M)
      Y(M-1) = Y1
      IFLG = 1
      GO TO 102
  121 CONTINUE
      IF (IA > 0) THEN
         RT = AA(IA)
         IA = IA - 1
         IFLG = 1
!
! SCALAR MULTIPLICATION
!
         Y(:M) = RT*Y(:M)
      ENDIF
      IF (IFLG > 0) GO TO 102
      DO J = 1, M
         YY(J) = REAL(Y(J))
      END DO
      RETURN 
      END SUBROUTINE CPROD


      SUBROUTINE CPRODP(ND, BD, NM1, BM1, NM2, BM2, NA, AA, X, YY, M, A , B, C, D, U, Y)
      implicit none
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER , INTENT(IN) :: ND
      INTEGER , INTENT(IN) :: NM1
      INTEGER , INTENT(IN) :: NM2
      INTEGER , INTENT(IN) :: NA
      INTEGER , INTENT(IN) :: M
      REAL(fish_kind) , INTENT(IN) :: BM1(*)
      REAL(fish_kind) , INTENT(IN) :: BM2(*)
      REAL(fish_kind) , INTENT(IN) :: AA(*)
      REAL(fish_kind) , INTENT(IN) :: X(*)
      REAL(fish_kind) , INTENT(OUT) :: YY(*)
      REAL(fish_kind) , INTENT(IN) :: A(*)
      REAL(fish_kind) , INTENT(IN) :: B(*)
      REAL(fish_kind) , INTENT(IN) :: C(*)
      COMPLEX , INTENT(IN) :: BD(*)
      COMPLEX , INTENT(INOUT) :: D(*)
      COMPLEX , INTENT(INOUT) :: U(*)
      COMPLEX , INTENT(INOUT) :: Y(*)
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: J, MM, MM2, ID, M1, M2, IA, IFLG, K
      REAL(fish_kind) :: RT
      COMPLEX :: V, DEN, BH, YM, AM, Y1, Y2, YH, CRT
!-----------------------------------------------
!
! PRODP APPLIES A SEQUENCE OF MATRIX OPERATIONS TO THE VECTOR X AND
! STORES THE RESULT IN YY       PERIODIC BOUNDARY CONDITIONS
! AND  COMPLEX  CASE
!
! BD,BM1,BM2 ARE ARRAYS CONTAINING ROOTS OF CERTIAN B POLYNOMIALS
! ND,NM1,NM2 ARE THE LENGTHS OF THE ARRAYS BD,BM1,BM2 RESPECTIVELY
! AA   ARRAY CONTAINING SCALAR MULTIPLIERS OF THE VECTOR X
! NA IS THE LENGTH OF THE ARRAY AA
! X,YY THE MATRIX OPERATIONS ARE APPLIED TO X AND THE RESULT IS YY
! A,B,C  ARE ARRAYS WHICH CONTAIN THE TRIDIAGONAL MATRIX
! M  IS THE ORDER OF THE MATRIX
! D,U,Y ARE WORKING ARRAYS
! ISGN  DETERMINES WHETHER OR NOT A CHANGE IN SIGN IS MADE
!
      DO J = 1, M
         Y(J) = CMPLX(X(J),0.)
      END DO
      MM = M - 1
      MM2 = M - 2
      ID = ND
      M1 = NM1
      M2 = NM2
      IA = NA
  102 CONTINUE
      IFLG = 0
      IF (ID > 0) THEN
         CRT = BD(ID)
         ID = ID - 1
         IFLG = 1
!
! BEGIN SOLUTION TO SYSTEM
!
         BH = B(M) - CRT
         YM = Y(M)
         DEN = B(1) - CRT
         D(1) = C(1)/DEN
         U(1) = A(1)/DEN
         Y(1) = Y(1)/DEN
         V = CMPLX(C(M),0.)
         IF (MM2 - 2 >= 0) THEN
            DO J = 2, MM2
               DEN = B(J) - CRT - A(J)*D(J-1)
               D(J) = C(J)/DEN
               U(J) = -A(J)*U(J-1)/DEN
               Y(J) = (Y(J)-A(J)*Y(J-1))/DEN
               BH = BH - V*U(J-1)
               YM = YM - V*Y(J-1)
               V = -V*D(J-1)
            END DO
         ENDIF
         DEN = B(M-1) - CRT - A(M-1)*D(M-2)
         D(M-1) = (C(M-1)-A(M-1)*U(M-2))/DEN
         Y(M-1) = (Y(M-1)-A(M-1)*Y(M-2))/DEN
         AM = A(M) - V*D(M-2)
         BH = BH - V*U(M-2)
         YM = YM - V*Y(M-2)
         DEN = BH - AM*D(M-1)
         IF (CABS(DEN) /= 0.) THEN
            Y(M) = (YM - AM*Y(M-1))/DEN
         ELSE
            Y(M) = (1.,0.)
         ENDIF
         Y(M-1) = Y(M-1) - D(M-1)*Y(M)
         DO J = 2, MM
            K = M - J
            Y(K) = Y(K) - D(K)*Y(K+1) - U(K)*Y(M)
         END DO
      ENDIF
      IF (M1 <= 0) THEN
         IF (M2 <= 0) GO TO 123
         RT = BM2(M2)
         M2 = M2 - 1
      ELSE
         IF (M2 <= 0) THEN
            RT = BM1(M1)
            M1 = M1 - 1
         ELSE
            IF (ABS(BM1(M1)) - ABS(BM2(M2)) > 0.) THEN
               RT = BM1(M1)
               M1 = M1 - 1
            ELSE
               RT = BM2(M2)
               M2 = M2 - 1
!
! MATRIX MULTIPLICATION
!
            ENDIF
         ENDIF
      ENDIF
      YH = Y(1)
      Y1 = (B(1)-RT)*Y(1) + C(1)*Y(2) + A(1)*Y(M)
      IF (MM - 2 >= 0) THEN
         DO J = 2, MM
            Y2 = A(J)*Y(J-1) + (B(J)-RT)*Y(J) + C(J)*Y(J+1)
            Y(J-1) = Y1
            Y1 = Y2
         END DO
      ENDIF
      Y(M) = A(M)*Y(M-1) + (B(M)-RT)*Y(M) + C(M)*YH
      Y(M-1) = Y1
      IFLG = 1
      GO TO 102
  123 CONTINUE
      IF (IA > 0) THEN
         RT = AA(IA)
         IA = IA - 1
         IFLG = 1
!
! SCALAR MULTIPLICATION
!
         Y(:M) = RT*Y(:M)
      ENDIF
      IF (IFLG > 0) GO TO 102
      DO J = 1, M
         YY(J) = REAL(Y(J))
      END DO
      RETURN 
      END SUBROUTINE CPRODP


      SUBROUTINE INDXA(I, IR, IDXA, NA)
      implicit none
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER , INTENT(IN) :: I
      INTEGER , INTENT(IN) :: IR
      INTEGER , INTENT(OUT) :: IDXA
      INTEGER , INTENT(OUT) :: NA
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
!-----------------------------------------------
      NA = 2**IR
      IDXA = I - NA + 1
      IF (I - NM > 0) THEN
         NA = 0
      ENDIF
      RETURN 
      END SUBROUTINE INDXA


      SUBROUTINE INDXB(I, IR, IDX, IDP)
      implicit none
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER , INTENT(IN) :: I
      INTEGER , INTENT(IN) :: IR
      INTEGER , INTENT(OUT) :: IDX
      INTEGER , INTENT(OUT) :: IDP
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: IZH, ID, IPL
!-----------------------------------------------
!
! B(IDX) IS THE LOCATION OF THE FIRST ROOT OF THE B(I,IR) POLYNOMIAL
!
      IDP = 0
      IDX = 0   !! SHOULD BE INITIALIZED!!!!
      IF (IR >= 0) THEN
         IF (IR <= 0) THEN
            IF (I - NM > 0) GO TO 107
            IDX = I
            IDP = 1
            RETURN 
         ENDIF
         IZH = 2**IR
         ID = I - IZH - IZH
         IDX = ID + ID + (IR - 1)*IK + IR + (IK - I)/IZH + 4
         IPL = IZH - 1
         IDP = IZH + IZH - 1
         IF (I - IPL - NM > 0) THEN
            IDP = 0
            RETURN 
         ENDIF
         IF (I + IPL - NM > 0) THEN
            IDP = NM + IPL - I + 1
         ENDIF
      ENDIF
  107 CONTINUE
      RETURN 
      END SUBROUTINE INDXB

!      SUBROUTINE INDXB (I,IR,IDX,IDP)                                   
!C                                                                       
!C B(IDX) IS THE LOCATION OF THE FIRST ROOT OF THE B(I,IR) POLYNOMIAL    
!C                                                                       
!      COMMON /CBLKT/  NPP        ,K          ,EPS        ,CNV        ,  
!     1                NM         ,NCMPLX     ,IK                        
!      IDP = 0                                                           
!      IF (IR) 107,101,103                                               
!  101 IF (I-NM) 102,102,107                                             
!  102 IDX = I                                                           
!      IDP = 1                                                           
!      RETURN                                                            
!  103 IZH = 2**IR                                                       
!      ID = I-IZH-IZH                                                    
!      IDX = ID+ID+(IR-1)*IK+IR+(IK-I)/IZH+4                             
!      IPL = IZH-1                                                       
!      IDP = IZH+IZH-1                                                   
!      IF (I-IPL-NM) 105,105,104                                         
!  104 IDP = 0                                                           
!      RETURN                                                            
!  105 IF (I+IPL-NM) 107,107,106                                         
!  106 IDP = NM+IPL-I+1                                                  
!  107 RETURN                                                            
!      END                                                               

      SUBROUTINE INDXC(I, IR, IDXC, NC)
      implicit none
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER , INTENT(IN) :: I
      INTEGER , INTENT(IN) :: IR
      INTEGER , INTENT(OUT) :: IDXC
      INTEGER , INTENT(OUT) :: NC
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
!-----------------------------------------------
      NC = 2**IR
      IDXC = I
      IF (IDXC + NC - 1 - NM > 0) THEN
         NC = 0
      ENDIF
      RETURN 
      END SUBROUTINE INDXC


      SUBROUTINE PPADD(N, IERROR, A, C, CBP, BP, BH)
      implicit none
!      external psgf,ppspf,ppsgf,bsrh
!   real psgf,ppspf,ppsgf,bsrh
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER , INTENT(IN) :: N
      INTEGER , INTENT(OUT) :: IERROR
      REAL(fish_kind)  :: A(*)
      REAL(fish_kind)  :: C(*)
      REAL(fish_kind) , INTENT(INOUT) :: BP(*)
      REAL(fish_kind)  :: BH(*)
      COMPLEX , INTENT(INOUT) :: CBP(*)
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER::IZ,IZM,IZM2,J,NT,MODIZ,IS,IF,IG,IT,ICV,I3,I2,NHALF
      REAL(fish_kind) :: R4, R5, R6, SCNV, XL, DB, SGN, XR, XM, PSG
      COMPLEX :: CF, CX, FSG, HSG, DD, F, FP, FPP, CDIS, R1, R2, R3
!-----------------------------------------------
!
!     PPADD COMPUTES THE EIGENVALUES OF THE PERIODIC TRIDIAGONAL MATRIX
!     WITH COEFFICIENTS AN,BN,CN
!
! N IS THE ORDER OF THE BH AND BP POLYNOMIALS
!     ON OUTPUT BP CONTIANS THE EIGENVALUES
! CBP IS THE SAME AS BP EXCEPT TYPE COMPLEX
! BH IS USED TO TEMPORARILY STORE THE ROOTS OF THE B HAT POLYNOMIAL
! WHICH ENTERS THROUGH BP
!
      SCNV = SQRT(CNV)
      IZ = N
      IZM = IZ - 1
      IZM2 = IZ - 2
      IF (BP(N) - BP(1) <= 0.) THEN
         IF (BP(N) - BP(1) == 0.) GO TO 142
         BH(:N) = BP(N:1:(-1))
      ELSE
         BH(:N) = BP(:N)
      ENDIF
      NCMPLX = 0
      MODIZ = MOD(IZ,2)
      IS = 1
      IF (MODIZ /= 0) THEN
         IF (A(1) < 0.) GO TO 110
         IF (A(1) == 0.) GO TO 142
      ENDIF
      XL = BH(1)
      DB = BH(3) - BH(1)
      XL = XL - DB
      R4 = PSGF(XL,IZ,C,A,BH)
      DO WHILE(R4 <= 0.)
         XL = XL - DB
         R4 = PSGF(XL,IZ,C,A,BH)
      END DO
      SGN = -1.
      CBP(1) = CMPLX(BSRH(XL,BH(1),IZ,C,A,BH,PSGF,SGN),0.)
      BP(1) = REAL(CBP(1))
      IS = 2
  110 CONTINUE
      IF = IZ - 1
      IF (MODIZ /= 0) THEN
         IF (A(1) > 0.) GO TO 115
         IF (A(1) == 0.) GO TO 142
      ENDIF
      XR = BH(IZ)
      DB = BH(IZ) - BH(IZ-2)
      XR = XR + DB
      R5 = PSGF(XR,IZ,C,A,BH)
      DO WHILE(R5 < 0.)
         XR = XR + DB
         R5 = PSGF(XR,IZ,C,A,BH)
      END DO
      SGN = 1.
      CBP(IZ) = CMPLX(BSRH(BH(IZ),XR,IZ,C,A,BH,PSGF,SGN),0.)
      IF = IZ - 2
  115 CONTINUE
      DO IG = IS, IF, 2
         XL = BH(IG)
         XR = BH(IG+1)
         SGN = -1.
         XM = BSRH(XL,XR,IZ,C,A,BH,PPSPF,SGN)
         PSG = PSGF(XM,IZ,C,A,BH)
         IF (ABS(PSG) - EPS <= 0.) GO TO 118
         R6 = PSG*PPSGF(XM,IZ,C,A,BH)
         IF (R6 > 0.) GO TO 119
         IF (R6 == 0.) GO TO 118
         SGN = 1.
         CBP(IG) = CMPLX(BSRH(BH(IG),XM,IZ,C,A,BH,PSGF,SGN),0.)
!        bp(ig) = real(cbp(ig))
         SGN = -1.
         CBP(IG+1) = CMPLX(BSRH(XM,BH(IG+1),IZ,C,A,BH,PSGF,SGN),0.)
!        bp(ig) = real(cbp(ig))
!        bp(ig+1) = real(cbp(ig+1))
         CYCLE 
!
!     CASE OF A MULTIPLE ZERO
!
  118    CONTINUE
         CBP(IG) = CMPLX(XM,0.)
         CBP(IG+1) = CMPLX(XM,0.)
!        bp(ig) = real(cbp(ig))
!        bp(ig+1) = real(cbp(ig+1))
         CYCLE 
!
!     CASE OF A COMPLEX ZERO
!
  119    CONTINUE
         IT = 0
         ICV = 0
         CX = CMPLX(XM,0.)
  120    CONTINUE
         FSG = (1.,0.)
         HSG = (1.,0.)
         FP = (0.,0.)
         FPP = (0.,0.)
         DO J = 1, IZ
            DD = 1./(CX - BH(J))
            FSG = FSG*A(J)*DD
            HSG = HSG*C(J)*DD
            FP = FP + DD
            FPP = FPP - DD*DD
         END DO
         IF (MODIZ == 0) THEN
            F = (1.,0.) - FSG - HSG
         ELSE
            F = (1.,0.) + FSG + HSG
         ENDIF
         I3 = 0
         IF (CABS(FP) > 0.) THEN
            I3 = 1
            R3 = -F/FP
         ENDIF
         I2 = 0
         IF (CABS(FPP) > 0.) THEN
            I2 = 1
            CDIS = CSQRT(FP**2 - 2.*F*FPP)
            R1 = CDIS - FP
            R2 = (-FP) - CDIS
            IF (CABS(R1) - CABS(R2) > 0.) THEN
               R1 = R1/FPP
            ELSE
               R1 = R2/FPP
            ENDIF
            R2 = 2.*F/FPP/R1
            IF (CABS(R2) < CABS(R1)) R1 = R2
            IF (I3 <= 0) GO TO 133
            IF (CABS(R3) < CABS(R1)) R1 = R3
            GO TO 133
         ENDIF
         R1 = R3
  133    CONTINUE
         CX = CX + R1
         IT = IT + 1
         IF (IT > 50) GO TO 142
         IF (CABS(R1) > SCNV) GO TO 120
         IF (ICV > 0) GO TO 135
         ICV = 1
         GO TO 120
  135    CONTINUE
         CBP(IG) = CX
         CBP(IG+1) = CONJG(CX)
      END DO
      IF (CABS(CBP(N)) - CABS(CBP(1)) <= 0.) THEN
         IF (CABS(CBP(N)) - CABS(CBP(1)) == 0.) GO TO 142
         NHALF = N/2
         DO J = 1, NHALF
            NT = N - J
            CX = CBP(J)
            CBP(J) = CBP(NT+1)
            CBP(NT+1) = CX
         END DO
      ENDIF
      NCMPLX = 1
      DO J = 2, IZ
         IF (AIMAG(CBP(J)) /= 0.) GO TO 143
      END DO
      NCMPLX = 0
      BP(1) = REAL(CBP(1))
      DO J = 2, IZ
         BP(J) = REAL(CBP(J))
      END DO
      GO TO 143
  142 CONTINUE
      IERROR = 4
  143 CONTINUE
      RETURN 
      END SUBROUTINE PPADD


      SUBROUTINE PROD(ND,BD,NM1,BM1,NM2,BM2,NA,AA,X,Y,M,A,B,C,D,W,U)
      implicit none
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER , INTENT(IN) :: ND
      INTEGER , INTENT(IN) :: NM1
      INTEGER , INTENT(IN) :: NM2
      INTEGER , INTENT(IN) :: NA
      INTEGER , INTENT(IN) :: M
      REAL(fish_kind) , INTENT(IN) :: BD(*)
      REAL(fish_kind) , INTENT(IN) :: BM1(*)
      REAL(fish_kind) , INTENT(IN) :: BM2(*)
      REAL(fish_kind) , INTENT(IN) :: AA(*)
      REAL(fish_kind) , INTENT(IN) :: X(*)
      REAL(fish_kind) , INTENT(INOUT) :: Y(*)
      REAL(fish_kind) , INTENT(IN) :: A(*)
      REAL(fish_kind) , INTENT(IN) :: B(*)
      REAL(fish_kind) , INTENT(IN) :: C(*)
      REAL(fish_kind) , INTENT(INOUT) :: D(*)
      REAL(fish_kind) , INTENT(INOUT) :: W(*)
      REAL(fish_kind)  :: U(*)
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: J, MM, ID, IBR, M1, M2, IA, K
      REAL(fish_kind) :: RT, DEN
!-----------------------------------------------
!
! PROD APPLIES A SEQUENCE OF MATRIX OPERATIONS TO THE VECTOR X AND
! STORES THE RESULT IN Y
! BD,BM1,BM2 ARE ARRAYS CONTAINING ROOTS OF CERTIAN B POLYNOMIALS
! ND,NM1,NM2 ARE THE LENGTHS OF THE ARRAYS BD,BM1,BM2 RESPECTIVELY
! AA   ARRAY CONTAINING SCALAR MULTIPLIERS OF THE VECTOR X
! NA IS THE LENGTH OF THE ARRAY AA
! X,Y  THE MATRIX OPERATIONS ARE APPLIED TO X AND THE RESULT IS Y
! A,B,C  ARE ARRAYS WHICH CONTAIN THE TRIDIAGONAL MATRIX
! M  IS THE ORDER OF THE MATRIX
! D,W,U ARE WORKING ARRAYS
! IS  DETERMINES WHETHER OR NOT A CHANGE IN SIGN IS MADE
!
      W(:M) = X(:M)
      Y(:M) = W(:M)
      MM = M - 1
      ID = ND
      IBR = 0
      M1 = NM1
      M2 = NM2
      IA = NA
  102 CONTINUE
      IF (IA > 0) THEN
         RT = AA(IA)
         IF (ND == 0) RT = -RT
         IA = IA - 1
!
! SCALAR MULTIPLICATION
!
         Y(:M) = RT*W(:M)
      ENDIF
      IF (ID <= 0) GO TO 125
      RT = BD(ID)
      ID = ID - 1
      IF (ID == 0) IBR = 1
!
! BEGIN SOLUTION TO SYSTEM
!
      D(M) = A(M)/(B(M)-RT)
      W(M) = Y(M)/(B(M)-RT)
      DO J = 2, MM
         K = M - J
         DEN = B(K+1) - RT - C(K+1)*D(K+2)
         D(K+1) = A(K+1)/DEN
         W(K+1) = (Y(K+1)-C(K+1)*W(K+2))/DEN
      END DO
      DEN = B(1) - RT - C(1)*D(2)
      W(1) = 1.
      IF (DEN /= 0.) THEN
         W(1) = (Y(1)-C(1)*W(2))/DEN
      ENDIF
      DO J = 2, M
         W(J) = W(J) - D(J)*W(J-1)
      END DO
      IF (NA > 0) GO TO 102
      GO TO 113
  111 CONTINUE
      Y(:M) = W(:M)
      IBR = 1
      GO TO 102
  113 CONTINUE
      IF (M1 <= 0) THEN
         IF (M2 <= 0) GO TO 111
      ELSE
         IF (M2 > 0) THEN
            IF (ABS(BM1(M1)) - ABS(BM2(M2)) <= 0.) GO TO 120
         ENDIF
         IF (IBR <= 0) THEN
            IF (ABS(BM1(M1)-BD(ID)) - ABS(BM1(M1)-RT) < 0.) GO TO 111
         ENDIF
         RT = RT - BM1(M1)
         M1 = M1 - 1
         GO TO 123
      ENDIF
  120 CONTINUE
      IF (IBR <= 0) THEN
         IF (ABS(BM2(M2)-BD(ID)) - ABS(BM2(M2)-RT) < 0.) GO TO 111
      ENDIF
      RT = RT - BM2(M2)
      M2 = M2 - 1
  123 CONTINUE
      Y(:M) = Y(:M) + RT*W(:M)
      GO TO 102
  125 CONTINUE
      RETURN 
      END SUBROUTINE PROD


      SUBROUTINE PRODP(ND,BD,NM1,BM1,NM2,BM2,NA,AA,X,Y,M,A,B,C,D,U,W)
      implicit none
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER , INTENT(IN) :: ND
      INTEGER , INTENT(IN) :: NM1
      INTEGER , INTENT(IN) :: NM2
      INTEGER , INTENT(IN) :: NA
      INTEGER , INTENT(IN) :: M
      REAL(fish_kind) , INTENT(IN) :: BD(*)
      REAL(fish_kind) , INTENT(IN) :: BM1(*)
      REAL(fish_kind) , INTENT(IN) :: BM2(*)
      REAL(fish_kind) , INTENT(IN) :: AA(*)
      REAL(fish_kind) , INTENT(IN) :: X(*)
      REAL(fish_kind) , INTENT(INOUT) :: Y(*)
      REAL(fish_kind) , INTENT(IN) :: A(*)
      REAL(fish_kind) , INTENT(IN) :: B(*)
      REAL(fish_kind) , INTENT(IN) :: C(*)
      REAL(fish_kind) , INTENT(INOUT) :: D(*)
      REAL(fish_kind) , INTENT(INOUT) :: U(*)
      REAL(fish_kind) , INTENT(INOUT) :: W(*)
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: J, MM, MM2, ID, IBR, M1, M2, IA, K
      REAL(fish_kind) :: RT, BH, YM, DEN, V, AM
!-----------------------------------------------
!
! PRODP APPLIES A SEQUENCE OF MATRIX OPERATIONS TO THE VECTOR X AND
! STORES THE RESULT IN Y        PERIODIC BOUNDARY CONDITIONS
!
! BD,BM1,BM2 ARE ARRAYS CONTAINING ROOTS OF CERTIAN B POLYNOMIALS
! ND,NM1,NM2 ARE THE LENGTHS OF THE ARRAYS BD,BM1,BM2 RESPECTIVELY
! AA   ARRAY CONTAINING SCALAR MULTIPLIERS OF THE VECTOR X
! NA IS THE LENGTH OF THE ARRAY AA
! X,Y  THE MATRIX OPERATIONS ARE APPLIED TO X AND THE RESULT IS Y
! A,B,C  ARE ARRAYS WHICH CONTAIN THE TRIDIAGONAL MATRIX
! M  IS THE ORDER OF THE MATRIX
! D,U,W ARE WORKING ARRAYS
! IS  DETERMINES WHETHER OR NOT A CHANGE IN SIGN IS MADE
!
 
      Y(:M) = X(:M)
      W(:M) = Y(:M)
      MM = M - 1
      MM2 = M - 2
      ID = ND
      IBR = 0
      M1 = NM1
      M2 = NM2
      IA = NA
  102 CONTINUE
      IF (IA > 0) THEN
         RT = AA(IA)
         IF (ND == 0) RT = -RT
         IA = IA - 1
         Y(:M) = RT*W(:M)
      ENDIF
      IF (ID <= 0) GO TO 128
      RT = BD(ID)
      ID = ID - 1
      IF (ID == 0) IBR = 1
!
! BEGIN SOLUTION TO SYSTEM
!
      BH = B(M) - RT
      YM = Y(M)
      DEN = B(1) - RT
      D(1) = C(1)/DEN
      U(1) = A(1)/DEN
      W(1) = Y(1)/DEN
      V = C(M)
      IF (MM2 - 2 >= 0) THEN
         DO J = 2, MM2
            DEN = B(J) - RT - A(J)*D(J-1)
            D(J) = C(J)/DEN
            U(J) = -A(J)*U(J-1)/DEN
            W(J) = (Y(J)-A(J)*W(J-1))/DEN
            BH = BH - V*U(J-1)
            YM = YM - V*W(J-1)
            V = -V*D(J-1)
         END DO
      ENDIF
      DEN = B(M-1) - RT - A(M-1)*D(M-2)
      D(M-1) = (C(M-1)-A(M-1)*U(M-2))/DEN
      W(M-1) = (Y(M-1)-A(M-1)*W(M-2))/DEN
      AM = A(M) - V*D(M-2)
      BH = BH - V*U(M-2)
      YM = YM - V*W(M-2)
      DEN = BH - AM*D(M-1)
      IF (DEN /= 0.) THEN
         W(M) = (YM - AM*W(M-1))/DEN
      ELSE
         W(M) = 1.
      ENDIF
      W(M-1) = W(M-1) - D(M-1)*W(M)
      DO J = 2, MM
         K = M - J
         W(K) = W(K) - D(K)*W(K+1) - U(K)*W(M)
      END DO
      IF (NA > 0) GO TO 102
      GO TO 116
  114 CONTINUE
      Y(:M) = W(:M)
      IBR = 1
      GO TO 102
  116 CONTINUE
      IF (M1 <= 0) THEN
         IF (M2 <= 0) GO TO 114
      ELSE
         IF (M2 > 0) THEN
            IF (ABS(BM1(M1)) - ABS(BM2(M2)) <= 0.) GO TO 123
         ENDIF
         IF (IBR <= 0) THEN
            IF (ABS(BM1(M1)-BD(ID)) - ABS(BM1(M1)-RT) < 0.) GO TO 114
         ENDIF
         RT = RT - BM1(M1)
         M1 = M1 - 1
         GO TO 126
      ENDIF
  123 CONTINUE
      IF (IBR <= 0) THEN
         IF (ABS(BM2(M2)-BD(ID)) - ABS(BM2(M2)-RT) < 0.) GO TO 114
      ENDIF
      RT = RT - BM2(M2)
      M2 = M2 - 1
  126 CONTINUE
      Y(:M) = Y(:M) + RT*W(:M)
      GO TO 102
  128 CONTINUE
      RETURN 
      END SUBROUTINE PRODP

      !*********************************************************************************

      SUBROUTINE TEVLS(N, D, E2, IERR)
      implicit none
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER , INTENT(IN) :: N
      INTEGER , INTENT(OUT) :: IERR
      REAL(fish_kind) , INTENT(INOUT) :: D(N)
      REAL(fish_kind) , INTENT(INOUT) :: E2(N)
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: I, J, L, M, L1, NHALF, NTOP
      REAL(fish_kind) :: B, C, F, G, H, P, R, S, DHOLD
!-----------------------------------------------
!
!
!     REAL SQRT,ABS,SIGN
!
!
!     THIS SUBROUTINE IS A MODIFICATION OF THE EISPACK SUBROUTINE TQLRAT
!     ALGORITHM 464, COMM. ACM 16, 689(1973) BY REINSCH.
!
!     THIS SUBROUTINE FINDS THE EIGENVALUES OF A SYMMETRIC
!     TRIDIAGONAL MATRIX BY THE RATIONAL QL METHOD.
!
!     ON INPUT-
!
!        N IS THE ORDER OF THE MATRIX,
!
!        D CONTAINS THE DIAGONAL ELEMENTS OF THE INPUT MATRIX,
!
!        E2 CONTAINS THE                SUBDIAGONAL ELEMENTS OF THE
!          INPUT MATRIX IN ITS LAST N-1 POSITIONS.  E2(1) IS ARBITRARY.
!
!      ON OUTPUT-
!
!        D CONTAINS THE EIGENVALUES IN ASCENDING ORDER.  IF AN
!          ERROR EXIT IS MADE, THE EIGENVALUES ARE CORRECT AND
!          ORDERED FOR INDICES 1,2,...IERR-1, BUT MAY NOT BE
!          THE SMALLEST EIGENVALUES,
!
!        E2 HAS BEEN DESTROYED,
!
!        IERR IS SET TO
!          ZERO       FOR NORMAL RETURN,
!          J          IF THE J-TH EIGENVALUE HAS NOT BEEN
!                     DETERMINED AFTER 30 ITERATIONS.
!
!     QUESTIONS AND COMMENTS SHOULD BE DIRECTED TO B. S. GARBOW,
!     APPLIED MATHEMATICS DIVISION, ARGONNE NATIONAL LABORATORY
!
!
!     ********** MACHEP IS A MACHINE DEPENDENT PARAMETER SPECIFYING
!                THE RELATIVE PRECISION OF FLOATING POINT ARITHMETIC.
!
!                **********
!      
      IERR = 0
      IF (N /= 1) THEN
!
         E2(:N-1) = E2(2:N)*E2(2:N)
!
         F = 0.0
         B = 0.0
         E2(N) = 0.0
!
         DO L = 1, N
            J = 0
            H = EPS*(ABS(D(L))+SQRT(E2(L)))
            IF (B < H) B = H
            C = B*B
!
!     ********** LOOK FOR SMALL SQUARED SUB-DIAGONAL ELEMENT **********
!
            DO M = L, N
               IF (E2(M) > C) CYCLE 
               EXIT 
!
!     ********** E2(N) IS ALWAYS ZERO, SO THERE IS NO EXIT
!                THROUGH THE BOTTOM OF THE LOOP **********
!
            END DO
!
            IF (M /= L) THEN
  105          CONTINUE
               IF (J == 30) GO TO 114
               J = J + 1
!
!     ********** FORM SHIFT **********
!
               L1 = L + 1
               S = SQRT(E2(L))
               G = D(L)
               P = (D(L1)-G)/(2.0*S)
               R = SQRT(P*P + 1.0)
               D(L) = S/(P + SIGN(R,P))
               H = G - D(L)
!
               D(L1:N) = D(L1:N) - H
!
               F = F + H
!
!     ********** RATIONAL QL TRANSFORMATION **********
!
               G = D(M)
               IF (abs(G) < 1.0D-100) G = B
               H = G
               S = 0.0
!
!     ********** FOR I=M-1 STEP -1 UNTIL L DO -- **********
!
               DO I = M-1, L, -1
                  P = G*H
                  R = P + E2(I)
                  E2(I+1) = S*R
                  S = E2(I)/R
                  D(I+1) = H + S*(H + D(I))
                  G = D(I) - E2(I)/G
                  IF (abs(G) < 1.0D-100) G = B
                  H = G*P/R
               END DO
!
               E2(L) = S*G
               D(L) = H
!
!     ********** GUARD AGAINST UNDERFLOWED H **********
!
               IF (abs(H) < 1.0D-100) GO TO 108
               IF (ABS(E2(L)) <= ABS(C/H)) GO TO 108
               E2(L) = H*E2(L)
               IF (abs(E2(L)) >= 1.0D-100)then
                 GO TO 105
               else
                 E2(L) = 0.0
               endif
            ENDIF
  108       CONTINUE
            P = D(L) + F
!
!     ********** ORDER EIGENVALUES **********
!
            IF (L /= 1) THEN
!
!     ********** FOR I=L STEP -1 UNTIL 2 DO -- **********
!
               DO I = L, 2, -1
                  IF (P >= D(I-1)) GO TO 111
                  D(I) = D(I-1)
               END DO
            ENDIF
!
            I = 1
  111       CONTINUE
            D(I) = P
         END DO    ! L = 1..N eigen values
!
         IF (ABS(D(N)) >= ABS(D(1))) GO TO 115
         NHALF = N/2
         DO I = 1, NHALF
            NTOP = N - I
            DHOLD = D(I)
            D(I) = D(NTOP+1)
            D(NTOP+1) = DHOLD
         END DO
         GO TO 115
!
!     ********** SET ERROR -- NO CONVERGENCE TO AN
!                EIGENVALUE AFTER 30 ITERATIONS **********
!
  114 CONTINUE
      IERR = L
      call msg('The L-th eigen value has not beein found after 30 iterations:',L)
      call set_error('The L-th eigen value has not beein found after 30 iterations','TEVLS')
      ENDIF
  115 CONTINUE
      RETURN 
!
!     ********** LAST CARD OF TQLRAT **********
!
!
! REVISION HISTORY---
!
! SEPTEMBER 1973    VERSION 1
! APRIL     1976    VERSION 2
! JANUARY   1978    VERSION 3
! DECEMBER  1979    VERSION 3.1
! FEBRUARY  1985    DOCUMENTATION UPGRADE
! NOVEMBER  1988    VERSION 3.2, FORTRAN 77 CHANGES
! June      2004    Version 5.0, Fortran 90 changes
!-----------------------------------------------------------------------
      END SUBROUTINE TEVLS

!
!     file fftpack.f
!
!
!     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
!     *                                                               *
!     *                  copyright (c) 2005 by UCAR                   *
!     *                                                               *
!     *       University Corporation for Atmospheric Research         *
!     *                                                               *
!     *                      all rights reserved                      *
!     *                                                               *
!     *                    FISHPACK90  version 1.1                    *
!     *                                                               *
!     *                 A Package of Fortran 77 and 90                *
!     *                                                               *
!     *                Subroutines and Example Programs               *
!     *                                                               *
!     *               for Modeling Geophysical Processes              *
!     *                                                               *
!     *                             by                                *
!     *                                                               *
!     *        John Adams, Paul Swarztrauber and Roland Sweet         *
!     *                                                               *
!     *                             of                                *
!     *                                                               *
!     *         the National Center for Atmospheric Research          *
!     *                                                               *
!     *                Boulder, Colorado  (80307)  U.S.A.             *
!     *                                                               *
!     *                   which is sponsored by                       *
!     *                                                               *
!     *              the National Science Foundation                  *
!     *                                                               *
!     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
!
!
! LATEST REVISION
! ---------------
!     June 2004      (VERSION 5.0) FORTRAN 90 CHANGES
!
! PURPOSE
! -------
!     THIS PACKAGE CONSISTS OF PROGRAMS WHICH PERFORM FAST FOURIER
!     TRANSFORMS FOR BOTH COMPLEX AND REAL PERIODIC SEQUENCES AND
!     CERTAIN OTHER SYMMETRIC SEQUENCES THAT ARE LISTED BELOW.
!
! USAGE
! -----
!     1.   RFFTI     INITIALIZE  RFFTF AND RFFTB
!     2.   RFFTF     FORWARD TRANSFORM OF A REAL PERIODIC SEQUENCE
!     3.   RFFTB     BACKWARD TRANSFORM OF A REAL COEFFICIENT ARRAY
!
!     4.   EZFFTI    INITIALIZE EZFFTF AND EZFFTB
!     5.   EZFFTF    A SIMPLIFIED REAL PERIODIC FORWARD TRANSFORM
!     6.   EZFFTB    A SIMPLIFIED REAL PERIODIC BACKWARD TRANSFORM
!
!     7.   SINTI     INITIALIZE SINT
!     8.   SINT      SINE TRANSFORM OF A REAL ODD SEQUENCE
!
!     9.   COSTI     INITIALIZE COST
!     10.  COST      COSINE TRANSFORM OF A REAL EVEN SEQUENCE
!
!     11.  SINQI     INITIALIZE SINQF AND SINQB
!     12.  SINQF     FORWARD SINE TRANSFORM WITH ODD WAVE NUMBERS
!     13.  SINQB     UNNORMALIZED INVERSE OF SINQF
!
!     14.  COSQI     INITIALIZE COSQF AND COSQB
!     15.  COSQF     FORWARD COSINE TRANSFORM WITH ODD WAVE NUMBERS
!     16.  COSQB     UNNORMALIZED INVERSE OF COSQF
!
!     17.  CFFTI     INITIALIZE CFFTF AND CFFTB
!     18.  CFFTF     FORWARD TRANSFORM OF A COMPLEX PERIODIC SEQUENCE
!     19.  CFFTB     UNNORMALIZED INVERSE OF CFFTF
!
! SPECIAL CONDITIONS
! ------------------
!     BEFORE CALLING ROUTINES EZFFTB AND EZFFTF FOR THE FIRST TIME,
!     OR BEFORE CALLING EZFFTB AND EZFFTF WITH A DIFFERENT LENGTH,
!     USERS MUST INITIALIZE BY CALLING ROUTINE EZFFTI.
!
! I/O
! ---
!     NONE
!
! PRECISION
! ---------
!     NONE
!
! REQUIRED LIBRARY FILES
! ----------------------
!     NONE
!
! LANGUAGE
! --------
!     FORTRAN
!
! HISTORY
! -------
!     DEVELOPED AT NCAR IN BOULDER, COLORADO BY PAUL N. SWARZTRAUBER
!     OF THE SCIENTIFIC COMPUTING DIVISION.  RELEASED ON NCAR'S PUBLIC
!     SOFTWARE LIBRARIES IN JANUARY 1980.  MODIFIED MAY 29, 1985 TO
!     INCREASE EFFICIENCY.
!
! PORTABILITY
! -----------
!     FORTRAN 77
!
! **********************************************************************
!
!     SUBROUTINE RFFTI(N,WSAVE)
!
!     SUBROUTINE RFFTI INITIALIZES THE ARRAY WSAVE WHICH IS USED IN
!     BOTH RFFTF AND RFFTB. THE PRIME FACTORIZATION OF N TOGETHER WITH
!     A TABULATION OF THE TRIGONOMETRIC FUNCTIONS ARE COMPUTED AND
!     STORED IN WSAVE.
!
!     INPUT PARAMETER
!
!     N       THE LENGTH OF THE SEQUENCE TO BE TRANSFORMED.
!
!     OUTPUT PARAMETER
!
!     WSAVE   A WORK ARRAY WHICH MUST BE DIMENSIONED AT LEAST 2*N+15.
!             THE SAME WORK ARRAY CAN BE USED FOR BOTH RFFTF AND RFFTB
!             AS LONG AS N REMAINS UNCHANGED. DIFFERENT WSAVE ARRAYS
!             ARE REQUIRED FOR DIFFERENT VALUES OF N. THE CONTENTS OF
!             WSAVE MUST NOT BE CHANGED BETWEEN CALLS OF RFFTF OR RFFTB.
!
! **********************************************************************
!
!     SUBROUTINE RFFTF(N,R,WSAVE)
!
!     SUBROUTINE RFFTF COMPUTES THE FOURIER COEFFICIENTS OF A REAL
!     PERODIC SEQUENCE (FOURIER ANALYSIS). THE TRANSFORM IS DEFINED
!     BELOW AT OUTPUT PARAMETER R.
!
!     INPUT PARAMETERS
!
!     N       THE LENGTH OF THE ARRAY R TO BE TRANSFORMED.  THE METHOD
!             IS MOST EFFICIENT WHEN N IS A PRODUCT OF SMALL PRIMES.
!             N MAY CHANGE SO LONG AS DIFFERENT WORK ARRAYS ARE PROVIDED
!
!     R       A REAL ARRAY OF LENGTH N WHICH CONTAINS THE SEQUENCE
!             TO BE TRANSFORMED
!
!     WSAVE   A WORK ARRAY WHICH MUST BE DIMENSIONED AT LEAST 2*N+15.
!             IN THE PROGRAM THAT CALLS RFFTF. THE WSAVE ARRAY MUST BE
!             INITIALIZED BY CALLING SUBROUTINE RFFTI(N,WSAVE) AND A
!             DIFFERENT WSAVE ARRAY MUST BE USED FOR EACH DIFFERENT
!             VALUE OF N. THIS INITIALIZATION DOES NOT HAVE TO BE
!             REPEATED SO LONG AS N REMAINS UNCHANGED THUS SUBSEQUENT
!             TRANSFORMS CAN BE OBTAINED FASTER THAN THE FIRST.
!             THE SAME WSAVE ARRAY CAN BE USED BY RFFTF AND RFFTB.
!
!
!     OUTPUT PARAMETERS
!
!     R       R(1) = THE SUM FROM I=1 TO I=N OF R(I)
!
!             IF N IS EVEN SET L =N/2   , IF N IS ODD SET L = (N+1)/2
!
!               THEN FOR K = 2,...,L
!
!                  R(2*K-2) = THE SUM FROM I = 1 TO I = N OF
!
!                       R(I)*COS((K-1)*(I-1)*2*PI/N)
!
!                  R(2*K-1) = THE SUM FROM I = 1 TO I = N OF
!
!                      -R(I)*SIN((K-1)*(I-1)*2*PI/N)
!
!             IF N IS EVEN
!
!                  R(N) = THE SUM FROM I = 1 TO I = N OF
!
!                       (-1)**(I-1)*R(I)
!
!      *****  NOTE
!                  THIS TRANSFORM IS UNNORMALIZED SINCE A CALL OF RFFTF
!                  FOLLOWED BY A CALL OF RFFTB WILL MULTIPLY THE INPUT
!                  SEQUENCE BY N.
!
!     WSAVE   CONTAINS RESULTS WHICH MUST NOT BE DESTROYED BETWEEN
!             CALLS OF RFFTF OR RFFTB.
!
!
! **********************************************************************
!
!     SUBROUTINE RFFTB(N,R,WSAVE)
!
!     SUBROUTINE RFFTB COMPUTES THE REAL PERODIC SEQUENCE FROM ITS
!     FOURIER COEFFICIENTS (FOURIER SYNTHESIS). THE TRANSFORM IS DEFINED
!     BELOW AT OUTPUT PARAMETER R.
!
!     INPUT PARAMETERS
!
!     N       THE LENGTH OF THE ARRAY R TO BE TRANSFORMED.  THE METHOD
!             IS MOST EFFICIENT WHEN N IS A PRODUCT OF SMALL PRIMES.
!             N MAY CHANGE SO LONG AS DIFFERENT WORK ARRAYS ARE PROVIDED
!
!     R       A REAL ARRAY OF LENGTH N WHICH CONTAINS THE SEQUENCE
!             TO BE TRANSFORMED
!
!     WSAVE   A WORK ARRAY WHICH MUST BE DIMENSIONED AT LEAST 2*N+15.
!             IN THE PROGRAM THAT CALLS RFFTB. THE WSAVE ARRAY MUST BE
!             INITIALIZED BY CALLING SUBROUTINE RFFTI(N,WSAVE) AND A
!             DIFFERENT WSAVE ARRAY MUST BE USED FOR EACH DIFFERENT
!             VALUE OF N. THIS INITIALIZATION DOES NOT HAVE TO BE
!             REPEATED SO LONG AS N REMAINS UNCHANGED THUS SUBSEQUENT
!             TRANSFORMS CAN BE OBTAINED FASTER THAN THE FIRST.
!             THE SAME WSAVE ARRAY CAN BE USED BY RFFTF AND RFFTB.
!
!
!     OUTPUT PARAMETERS
!
!     R       FOR N EVEN AND FOR I = 1,...,N
!
!                  R(I) = R(1)+(-1)**(I-1)*R(N)
!
!                       PLUS THE SUM FROM K=2 TO K=N/2 OF
!
!                        2.*R(2*K-2)*COS((K-1)*(I-1)*2*PI/N)
!
!                       -2.*R(2*K-1)*SIN((K-1)*(I-1)*2*PI/N)
!
!             FOR N ODD AND FOR I = 1,...,N
!
!                  R(I) = R(1) PLUS THE SUM FROM K=2 TO K=(N+1)/2 OF
!
!                       2.*R(2*K-2)*COS((K-1)*(I-1)*2*PI/N)
!
!                      -2.*R(2*K-1)*SIN((K-1)*(I-1)*2*PI/N)
!
!      *****  NOTE
!                  THIS TRANSFORM IS UNNORMALIZED SINCE A CALL OF RFFTF
!                  FOLLOWED BY A CALL OF RFFTB WILL MULTIPLY THE INPUT
!                  SEQUENCE BY N.
!
!     WSAVE   CONTAINS RESULTS WHICH MUST NOT BE DESTROYED BETWEEN
!             CALLS OF RFFTB OR RFFTF.
!
!
! **********************************************************************
!
!     SUBROUTINE EZFFTI(N,WSAVE)
!
!     SUBROUTINE EZFFTI INITIALIZES THE ARRAY WSAVE WHICH IS USED IN
!     BOTH EZFFTF AND EZFFTB. THE PRIME FACTORIZATION OF N TOGETHER WITH
!     A TABULATION OF THE TRIGONOMETRIC FUNCTIONS ARE COMPUTED AND
!     STORED IN WSAVE.
!
!     INPUT PARAMETER
!
!     N       THE LENGTH OF THE SEQUENCE TO BE TRANSFORMED.
!
!     OUTPUT PARAMETER
!
!     WSAVE   A WORK ARRAY WHICH MUST BE DIMENSIONED AT LEAST 3*N+15.
!             THE SAME WORK ARRAY CAN BE USED FOR BOTH EZFFTF AND EZFFTB
!             AS LONG AS N REMAINS UNCHANGED. DIFFERENT WSAVE ARRAYS
!             ARE REQUIRED FOR DIFFERENT VALUES OF N.
!
!
! **********************************************************************
!
!     SUBROUTINE EZFFTF(N,R,AZERO,A,B,WSAVE)
!
!     SUBROUTINE EZFFTF COMPUTES THE FOURIER COEFFICIENTS OF A REAL
!     PERODIC SEQUENCE (FOURIER ANALYSIS). THE TRANSFORM IS DEFINED
!     BELOW AT OUTPUT PARAMETERS AZERO,A AND B. EZFFTF IS A SIMPLIFIED
!     BUT SLOWER VERSION OF RFFTF.
!
!     INPUT PARAMETERS
!
!     N       THE LENGTH OF THE ARRAY R TO BE TRANSFORMED.  THE METHOD
!             IS MUST EFFICIENT WHEN N IS THE PRODUCT OF SMALL PRIMES.
!
!     R       A REAL ARRAY OF LENGTH N WHICH CONTAINS THE SEQUENCE
!             TO BE TRANSFORMED. R IS NOT DESTROYED.
!
!
!     WSAVE   A WORK ARRAY WHICH MUST BE DIMENSIONED AT LEAST 3*N+15.
!             IN THE PROGRAM THAT CALLS EZFFTF. THE WSAVE ARRAY MUST BE
!             INITIALIZED BY CALLING SUBROUTINE EZFFTI(N,WSAVE) AND A
!             DIFFERENT WSAVE ARRAY MUST BE USED FOR EACH DIFFERENT
!             VALUE OF N. THIS INITIALIZATION DOES NOT HAVE TO BE
!             REPEATED SO LONG AS N REMAINS UNCHANGED THUS SUBSEQUENT
!             TRANSFORMS CAN BE OBTAINED FASTER THAN THE FIRST.
!             THE SAME WSAVE ARRAY CAN BE USED BY EZFFTF AND EZFFTB.
!
!     OUTPUT PARAMETERS
!
!     AZERO   THE SUM FROM I=1 TO I=N OF R(I)/N
!
!     A,B     FOR N EVEN B(N/2)=0. AND A(N/2) IS THE SUM FROM I=1 TO
!             I=N OF (-1)**(I-1)*R(I)/N
!
!             FOR N EVEN DEFINE KMAX=N/2-1
!             FOR N ODD  DEFINE KMAX=(N-1)/2
!
!             THEN FOR  K=1,...,KMAX
!
!                  A(K) EQUALS THE SUM FROM I=1 TO I=N OF
!
!                       2./N*R(I)*COS(K*(I-1)*2*PI/N)
!
!                  B(K) EQUALS THE SUM FROM I=1 TO I=N OF
!
!                       2./N*R(I)*SIN(K*(I-1)*2*PI/N)
!
!
! **********************************************************************
!
!     SUBROUTINE EZFFTB(N,R,AZERO,A,B,WSAVE)
!
!     SUBROUTINE EZFFTB COMPUTES A REAL PERODIC SEQUENCE FROM ITS
!     FOURIER COEFFICIENTS (FOURIER SYNTHESIS). THE TRANSFORM IS
!     DEFINED BELOW AT OUTPUT PARAMETER R. EZFFTB IS A SIMPLIFIED
!     BUT SLOWER VERSION OF RFFTB.
!
!     INPUT PARAMETERS
!
!     N       THE LENGTH OF THE OUTPUT ARRAY R.  THE METHOD IS MOST
!             EFFICIENT WHEN N IS THE PRODUCT OF SMALL PRIMES.
!
!     AZERO   THE CONSTANT FOURIER COEFFICIENT
!
!     A,B     ARRAYS WHICH CONTAIN THE REMAINING FOURIER COEFFICIENTS
!             THESE ARRAYS ARE NOT DESTROYED.
!
!             THE LENGTH OF THESE ARRAYS DEPENDS ON WHETHER N IS EVEN OR
!             ODD.
!
!             IF N IS EVEN N/2    LOCATIONS ARE REQUIRED
!             IF N IS ODD (N-1)/2 LOCATIONS ARE REQUIRED
!
!     WSAVE   A WORK ARRAY WHICH MUST BE DIMENSIONED AT LEAST 3*N+15.
!             IN THE PROGRAM THAT CALLS EZFFTB. THE WSAVE ARRAY MUST BE
!             INITIALIZED BY CALLING SUBROUTINE EZFFTI(N,WSAVE) AND A
!             DIFFERENT WSAVE ARRAY MUST BE USED FOR EACH DIFFERENT
!             VALUE OF N. THIS INITIALIZATION DOES NOT HAVE TO BE
!             REPEATED SO LONG AS N REMAINS UNCHANGED THUS SUBSEQUENT
!             TRANSFORMS CAN BE OBTAINED FASTER THAN THE FIRST.
!             THE SAME WSAVE ARRAY CAN BE USED BY EZFFTF AND EZFFTB.
!
!
!     OUTPUT PARAMETERS
!
!     R       IF N IS EVEN DEFINE KMAX=N/2
!             IF N IS ODD  DEFINE KMAX=(N-1)/2
!
!             THEN FOR I=1,...,N
!
!                  R(I)=AZERO PLUS THE SUM FROM K=1 TO K=KMAX OF
!
!                  A(K)*COS(K*(I-1)*2*PI/N)+B(K)*SIN(K*(I-1)*2*PI/N)
!
!     ********************* COMPLEX NOTATION **************************
!
!             FOR J=1,...,N
!
!             R(J) EQUALS THE SUM FROM K=-KMAX TO K=KMAX OF
!
!                  C(K)*EXP(I*K*(J-1)*2*PI/N)
!
!             WHERE
!
!                  C(K) = .5*CMPLX(A(K),-B(K))   FOR K=1,...,KMAX
!
!                  C(-K) = CONJG(C(K))
!
!                  C(0) = AZERO
!
!                       AND I=SQRT(-1)
!
!     *************** AMPLITUDE - PHASE NOTATION ***********************
!
!             FOR I=1,...,N
!
!             R(I) EQUALS AZERO PLUS THE SUM FROM K=1 TO K=KMAX OF
!
!                  ALPHA(K)*COS(K*(I-1)*2*PI/N+BETA(K))
!
!             WHERE
!
!                  ALPHA(K) = SQRT(A(K)*A(K)+B(K)*B(K))
!
!                  COS(BETA(K))=A(K)/ALPHA(K)
!
!                  SIN(BETA(K))=-B(K)/ALPHA(K)
!
! **********************************************************************
!
!     SUBROUTINE SINTI(N,WSAVE)
!
!     SUBROUTINE SINTI INITIALIZES THE ARRAY WSAVE WHICH IS USED IN
!     SUBROUTINE SINT. THE PRIME FACTORIZATION OF N TOGETHER WITH
!     A TABULATION OF THE TRIGONOMETRIC FUNCTIONS ARE COMPUTED AND
!     STORED IN WSAVE.
!
!     INPUT PARAMETER
!
!     N       THE LENGTH OF THE SEQUENCE TO BE TRANSFORMED.  THE METHOD
!             IS MOST EFFICIENT WHEN N+1 IS A PRODUCT OF SMALL PRIMES.
!
!     OUTPUT PARAMETER
!
!     WSAVE   A WORK ARRAY WITH AT LEAST INT(2.5*N+15) LOCATIONS.
!             DIFFERENT WSAVE ARRAYS ARE REQUIRED FOR DIFFERENT VALUES
!             OF N. THE CONTENTS OF WSAVE MUST NOT BE CHANGED BETWEEN
!             CALLS OF SINT.
!
! **********************************************************************
!
!     SUBROUTINE SINT(N,X,WSAVE)
!
!     SUBROUTINE SINT COMPUTES THE DISCRETE FOURIER SINE TRANSFORM
!     OF AN ODD SEQUENCE X(I). THE TRANSFORM IS DEFINED BELOW AT
!     OUTPUT PARAMETER X.
!
!     SINT IS THE UNNORMALIZED INVERSE OF ITSELF SINCE A CALL OF SINT
!     FOLLOWED BY ANOTHER CALL OF SINT WILL MULTIPLY THE INPUT SEQUENCE
!     X BY 2*(N+1).
!
!     THE ARRAY WSAVE WHICH IS USED BY SUBROUTINE SINT MUST BE
!     INITIALIZED BY CALLING SUBROUTINE SINTI(N,WSAVE).
!
!     INPUT PARAMETERS
!
!     N       THE LENGTH OF THE SEQUENCE TO BE TRANSFORMED.  THE METHOD
!             IS MOST EFFICIENT WHEN N+1 IS THE PRODUCT OF SMALL PRIMES.
!
!     X       AN ARRAY WHICH CONTAINS THE SEQUENCE TO BE TRANSFORMED
!
!
!     WSAVE   A WORK ARRAY WITH DIMENSION AT LEAST INT(2.5*N+15)
!             IN THE PROGRAM THAT CALLS SINT. THE WSAVE ARRAY MUST BE
!             INITIALIZED BY CALLING SUBROUTINE SINTI(N,WSAVE) AND A
!             DIFFERENT WSAVE ARRAY MUST BE USED FOR EACH DIFFERENT
!             VALUE OF N. THIS INITIALIZATION DOES NOT HAVE TO BE
!             REPEATED SO LONG AS N REMAINS UNCHANGED THUS SUBSEQUENT
!             TRANSFORMS CAN BE OBTAINED FASTER THAN THE FIRST.
!
!     OUTPUT PARAMETERS
!
!     X       FOR I=1,...,N
!
!                  X(I)= THE SUM FROM K=1 TO K=N
!
!                       2*X(K)*SIN(K*I*PI/(N+1))
!
!                  A CALL OF SINT FOLLOWED BY ANOTHER CALL OF
!                  SINT WILL MULTIPLY THE SEQUENCE X BY 2*(N+1).
!                  HENCE SINT IS THE UNNORMALIZED INVERSE
!                  OF ITSELF.
!
!     WSAVE   CONTAINS INITIALIZATION CALCULATIONS WHICH MUST NOT BE
!             DESTROYED BETWEEN CALLS OF SINT.
!
! **********************************************************************
!
!     SUBROUTINE COSTI(N,WSAVE)
!
!     SUBROUTINE COSTI INITIALIZES THE ARRAY WSAVE WHICH IS USED IN
!     SUBROUTINE COST. THE PRIME FACTORIZATION OF N TOGETHER WITH
!     A TABULATION OF THE TRIGONOMETRIC FUNCTIONS ARE COMPUTED AND
!     STORED IN WSAVE.
!
!     INPUT PARAMETER
!
!     N       THE LENGTH OF THE SEQUENCE TO BE TRANSFORMED.  THE METHOD
!             IS MOST EFFICIENT WHEN N-1 IS A PRODUCT OF SMALL PRIMES.
!
!     OUTPUT PARAMETER
!
!     WSAVE   A WORK ARRAY WHICH MUST BE DIMENSIONED AT LEAST 3*N+15.
!             DIFFERENT WSAVE ARRAYS ARE REQUIRED FOR DIFFERENT VALUES
!             OF N. THE CONTENTS OF WSAVE MUST NOT BE CHANGED BETWEEN
!             CALLS OF COST.
!
! **********************************************************************
!
!     SUBROUTINE COST(N,X,WSAVE)
!
!     SUBROUTINE COST COMPUTES THE DISCRETE FOURIER COSINE TRANSFORM
!     OF AN EVEN SEQUENCE X(I). THE TRANSFORM IS DEFINED BELOW AT OUTPUT
!     PARAMETER X.
!
!     COST IS THE UNNORMALIZED INVERSE OF ITSELF SINCE A CALL OF COST
!     FOLLOWED BY ANOTHER CALL OF COST WILL MULTIPLY THE INPUT SEQUENCE
!     X BY 2*(N-1). THE TRANSFORM IS DEFINED BELOW AT OUTPUT PARAMETER X
!
!     THE ARRAY WSAVE WHICH IS USED BY SUBROUTINE COST MUST BE
!     INITIALIZED BY CALLING SUBROUTINE COSTI(N,WSAVE).
!
!     INPUT PARAMETERS
!
!     N       THE LENGTH OF THE SEQUENCE X. N MUST BE GREATER THAN 1.
!             THE METHOD IS MOST EFFICIENT WHEN N-1 IS A PRODUCT OF
!             SMALL PRIMES.
!
!     X       AN ARRAY WHICH CONTAINS THE SEQUENCE TO BE TRANSFORMED
!
!     WSAVE   A WORK ARRAY WHICH MUST BE DIMENSIONED AT LEAST 3*N+15
!             IN THE PROGRAM THAT CALLS COST. THE WSAVE ARRAY MUST BE
!             INITIALIZED BY CALLING SUBROUTINE COSTI(N,WSAVE) AND A
!             DIFFERENT WSAVE ARRAY MUST BE USED FOR EACH DIFFERENT
!             VALUE OF N. THIS INITIALIZATION DOES NOT HAVE TO BE
!             REPEATED SO LONG AS N REMAINS UNCHANGED THUS SUBSEQUENT
!             TRANSFORMS CAN BE OBTAINED FASTER THAN THE FIRST.
!
!     OUTPUT PARAMETERS
!
!     X       FOR I=1,...,N
!
!                 X(I) = X(1)+(-1)**(I-1)*X(N)
!
!                  + THE SUM FROM K=2 TO K=N-1
!
!                      2*X(K)*COS((K-1)*(I-1)*PI/(N-1))
!
!                  A CALL OF COST FOLLOWED BY ANOTHER CALL OF
!                  COST WILL MULTIPLY THE SEQUENCE X BY 2*(N-1)
!                  HENCE COST IS THE UNNORMALIZED INVERSE
!                  OF ITSELF.
!
!     WSAVE   CONTAINS INITIALIZATION CALCULATIONS WHICH MUST NOT BE
!             DESTROYED BETWEEN CALLS OF COST.
!
! **********************************************************************
!
!     SUBROUTINE SINQI(N,WSAVE)
!
!     SUBROUTINE SINQI INITIALIZES THE ARRAY WSAVE WHICH IS USED IN
!     BOTH SINQF AND SINQB. THE PRIME FACTORIZATION OF N TOGETHER WITH
!     A TABULATION OF THE TRIGONOMETRIC FUNCTIONS ARE COMPUTED AND
!     STORED IN WSAVE.
!
!     INPUT PARAMETER
!
!     N       THE LENGTH OF THE SEQUENCE TO BE TRANSFORMED. THE METHOD
!             IS MOST EFFICIENT WHEN N IS A PRODUCT OF SMALL PRIMES.
!
!     OUTPUT PARAMETER
!
!     WSAVE   A WORK ARRAY WHICH MUST BE DIMENSIONED AT LEAST 3*N+15.
!             THE SAME WORK ARRAY CAN BE USED FOR BOTH SINQF AND SINQB
!             AS LONG AS N REMAINS UNCHANGED. DIFFERENT WSAVE ARRAYS
!             ARE REQUIRED FOR DIFFERENT VALUES OF N. THE CONTENTS OF
!             WSAVE MUST NOT BE CHANGED BETWEEN CALLS OF SINQF OR SINQB.
!
! **********************************************************************
!
!     SUBROUTINE SINQF(N,X,WSAVE)
!
!     SUBROUTINE SINQF COMPUTES THE FAST FOURIER TRANSFORM OF QUARTER
!     WAVE DATA. THAT IS , SINQF COMPUTES THE COEFFICIENTS IN A SINE
!     SERIES REPRESENTATION WITH ONLY ODD WAVE NUMBERS. THE TRANSFORM
!     IS DEFINED BELOW AT OUTPUT PARAMETER X.
!
!     SINQB IS THE UNNORMALIZED INVERSE OF SINQF SINCE A CALL OF SINQF
!     FOLLOWED BY A CALL OF SINQB WILL MULTIPLY THE INPUT SEQUENCE X
!     BY 4*N.
!
!     THE ARRAY WSAVE WHICH IS USED BY SUBROUTINE SINQF MUST BE
!     INITIALIZED BY CALLING SUBROUTINE SINQI(N,WSAVE).
!
!
!     INPUT PARAMETERS
!
!     N       THE LENGTH OF THE ARRAY X TO BE TRANSFORMED.  THE METHOD
!             IS MOST EFFICIENT WHEN N IS A PRODUCT OF SMALL PRIMES.
!
!     X       AN ARRAY WHICH CONTAINS THE SEQUENCE TO BE TRANSFORMED
!
!     WSAVE   A WORK ARRAY WHICH MUST BE DIMENSIONED AT LEAST 3*N+15.
!             IN THE PROGRAM THAT CALLS SINQF. THE WSAVE ARRAY MUST BE
!             INITIALIZED BY CALLING SUBROUTINE SINQI(N,WSAVE) AND A
!             DIFFERENT WSAVE ARRAY MUST BE USED FOR EACH DIFFERENT
!             VALUE OF N. THIS INITIALIZATION DOES NOT HAVE TO BE
!             REPEATED SO LONG AS N REMAINS UNCHANGED THUS SUBSEQUENT
!             TRANSFORMS CAN BE OBTAINED FASTER THAN THE FIRST.
!
!     OUTPUT PARAMETERS
!
!     X       FOR I=1,...,N
!
!                  X(I) = (-1)**(I-1)*X(N)
!
!                     + THE SUM FROM K=1 TO K=N-1 OF
!
!                     2*X(K)*SIN((2*I-1)*K*PI/(2*N))
!
!                  A CALL OF SINQF FOLLOWED BY A CALL OF
!                  SINQB WILL MULTIPLY THE SEQUENCE X BY 4*N.
!                  THEREFORE SINQB IS THE UNNORMALIZED INVERSE
!                  OF SINQF.
!
!     WSAVE   CONTAINS INITIALIZATION CALCULATIONS WHICH MUST NOT
!             BE DESTROYED BETWEEN CALLS OF SINQF OR SINQB.
!
! **********************************************************************
!
!     SUBROUTINE SINQB(N,X,WSAVE)
!
!     SUBROUTINE SINQB COMPUTES THE FAST FOURIER TRANSFORM OF QUARTER
!     WAVE DATA. THAT IS , SINQB COMPUTES A SEQUENCE FROM ITS
!     REPRESENTATION IN TERMS OF A SINE SERIES WITH ODD WAVE NUMBERS.
!     THE TRANSFORM IS DEFINED BELOW AT OUTPUT PARAMETER X.
!
!     SINQF IS THE UNNORMALIZED INVERSE OF SINQB SINCE A CALL OF SINQB
!     FOLLOWED BY A CALL OF SINQF WILL MULTIPLY THE INPUT SEQUENCE X
!     BY 4*N.
!
!     THE ARRAY WSAVE WHICH IS USED BY SUBROUTINE SINQB MUST BE
!     INITIALIZED BY CALLING SUBROUTINE SINQI(N,WSAVE).
!
!
!     INPUT PARAMETERS
!
!     N       THE LENGTH OF THE ARRAY X TO BE TRANSFORMED.  THE METHOD
!             IS MOST EFFICIENT WHEN N IS A PRODUCT OF SMALL PRIMES.
!
!     X       AN ARRAY WHICH CONTAINS THE SEQUENCE TO BE TRANSFORMED
!
!     WSAVE   A WORK ARRAY WHICH MUST BE DIMENSIONED AT LEAST 3*N+15.
!             IN THE PROGRAM THAT CALLS SINQB. THE WSAVE ARRAY MUST BE
!             INITIALIZED BY CALLING SUBROUTINE SINQI(N,WSAVE) AND A
!             DIFFERENT WSAVE ARRAY MUST BE USED FOR EACH DIFFERENT
!             VALUE OF N. THIS INITIALIZATION DOES NOT HAVE TO BE
!             REPEATED SO LONG AS N REMAINS UNCHANGED THUS SUBSEQUENT
!             TRANSFORMS CAN BE OBTAINED FASTER THAN THE FIRST.
!
!     OUTPUT PARAMETERS
!
!     X       FOR I=1,...,N
!
!                  X(I)= THE SUM FROM K=1 TO K=N OF
!
!                    4*X(K)*SIN((2K-1)*I*PI/(2*N))
!
!                  A CALL OF SINQB FOLLOWED BY A CALL OF
!                  SINQF WILL MULTIPLY THE SEQUENCE X BY 4*N.
!                  THEREFORE SINQF IS THE UNNORMALIZED INVERSE
!                  OF SINQB.
!
!     WSAVE   CONTAINS INITIALIZATION CALCULATIONS WHICH MUST NOT
!             BE DESTROYED BETWEEN CALLS OF SINQB OR SINQF.
!
! **********************************************************************
!
!     SUBROUTINE COSQI(N,WSAVE)
!
!     SUBROUTINE COSQI INITIALIZES THE ARRAY WSAVE WHICH IS USED IN
!     BOTH COSQF AND COSQB. THE PRIME FACTORIZATION OF N TOGETHER WITH
!     A TABULATION OF THE TRIGONOMETRIC FUNCTIONS ARE COMPUTED AND
!     STORED IN WSAVE.
!
!     INPUT PARAMETER
!
!     N       THE LENGTH OF THE ARRAY TO BE TRANSFORMED.  THE METHOD
!             IS MOST EFFICIENT WHEN N IS A PRODUCT OF SMALL PRIMES.
!
!     OUTPUT PARAMETER
!
!     WSAVE   A WORK ARRAY WHICH MUST BE DIMENSIONED AT LEAST 3*N+15.
!             THE SAME WORK ARRAY CAN BE USED FOR BOTH COSQF AND COSQB
!             AS LONG AS N REMAINS UNCHANGED. DIFFERENT WSAVE ARRAYS
!             ARE REQUIRED FOR DIFFERENT VALUES OF N. THE CONTENTS OF
!             WSAVE MUST NOT BE CHANGED BETWEEN CALLS OF COSQF OR COSQB.
!
! **********************************************************************
!
!     SUBROUTINE COSQF(N,X,WSAVE)
!
!     SUBROUTINE COSQF COMPUTES THE FAST FOURIER TRANSFORM OF QUARTER
!     WAVE DATA. THAT IS , COSQF COMPUTES THE COEFFICIENTS IN A COSINE
!     SERIES REPRESENTATION WITH ONLY ODD WAVE NUMBERS. THE TRANSFORM
!     IS DEFINED BELOW AT OUTPUT PARAMETER X
!
!     COSQF IS THE UNNORMALIZED INVERSE OF COSQB SINCE A CALL OF COSQF
!     FOLLOWED BY A CALL OF COSQB WILL MULTIPLY THE INPUT SEQUENCE X
!     BY 4*N.
!
!     THE ARRAY WSAVE WHICH IS USED BY SUBROUTINE COSQF MUST BE
!     INITIALIZED BY CALLING SUBROUTINE COSQI(N,WSAVE).
!
!
!     INPUT PARAMETERS
!
!     N       THE LENGTH OF THE ARRAY X TO BE TRANSFORMED.  THE METHOD
!             IS MOST EFFICIENT WHEN N IS A PRODUCT OF SMALL PRIMES.
!
!     X       AN ARRAY WHICH CONTAINS THE SEQUENCE TO BE TRANSFORMED
!
!     WSAVE   A WORK ARRAY WHICH MUST BE DIMENSIONED AT LEAST 3*N+15
!             IN THE PROGRAM THAT CALLS COSQF. THE WSAVE ARRAY MUST BE
!             INITIALIZED BY CALLING SUBROUTINE COSQI(N,WSAVE) AND A
!             DIFFERENT WSAVE ARRAY MUST BE USED FOR EACH DIFFERENT
!             VALUE OF N. THIS INITIALIZATION DOES NOT HAVE TO BE
!             REPEATED SO LONG AS N REMAINS UNCHANGED THUS SUBSEQUENT
!             TRANSFORMS CAN BE OBTAINED FASTER THAN THE FIRST.
!
!     OUTPUT PARAMETERS
!
!     X       FOR I=1,...,N
!
!                  X(I) = X(1) PLUS THE SUM FROM K=2 TO K=N OF
!
!                     2*X(K)*COS((2*I-1)*(K-1)*PI/(2*N))
!
!                  A CALL OF COSQF FOLLOWED BY A CALL OF
!                  COSQB WILL MULTIPLY THE SEQUENCE X BY 4*N.
!                  THEREFORE COSQB IS THE UNNORMALIZED INVERSE
!                  OF COSQF.
!
!     WSAVE   CONTAINS INITIALIZATION CALCULATIONS WHICH MUST NOT
!             BE DESTROYED BETWEEN CALLS OF COSQF OR COSQB.
!
! **********************************************************************
!
!     SUBROUTINE COSQB(N,X,WSAVE)
!
!     SUBROUTINE COSQB COMPUTES THE FAST FOURIER TRANSFORM OF QUARTER
!     WAVE DATA. THAT IS , COSQB COMPUTES A SEQUENCE FROM ITS
!     REPRESENTATION IN TERMS OF A COSINE SERIES WITH ODD WAVE NUMBERS.
!     THE TRANSFORM IS DEFINED BELOW AT OUTPUT PARAMETER X.
!
!     COSQB IS THE UNNORMALIZED INVERSE OF COSQF SINCE A CALL OF COSQB
!     FOLLOWED BY A CALL OF COSQF WILL MULTIPLY THE INPUT SEQUENCE X
!     BY 4*N.
!
!     THE ARRAY WSAVE WHICH IS USED BY SUBROUTINE COSQB MUST BE
!     INITIALIZED BY CALLING SUBROUTINE COSQI(N,WSAVE).
!
!
!     INPUT PARAMETERS
!
!     N       THE LENGTH OF THE ARRAY X TO BE TRANSFORMED.  THE METHOD
!             IS MOST EFFICIENT WHEN N IS A PRODUCT OF SMALL PRIMES.
!
!     X       AN ARRAY WHICH CONTAINS THE SEQUENCE TO BE TRANSFORMED
!
!     WSAVE   A WORK ARRAY THAT MUST BE DIMENSIONED AT LEAST 3*N+15
!             IN THE PROGRAM THAT CALLS COSQB. THE WSAVE ARRAY MUST BE
!             INITIALIZED BY CALLING SUBROUTINE COSQI(N,WSAVE) AND A
!             DIFFERENT WSAVE ARRAY MUST BE USED FOR EACH DIFFERENT
!             VALUE OF N. THIS INITIALIZATION DOES NOT HAVE TO BE
!             REPEATED SO LONG AS N REMAINS UNCHANGED THUS SUBSEQUENT
!             TRANSFORMS CAN BE OBTAINED FASTER THAN THE FIRST.
!
!     OUTPUT PARAMETERS
!
!     X       FOR I=1,...,N
!
!                  X(I)= THE SUM FROM K=1 TO K=N OF
!
!                    4*X(K)*COS((2*K-1)*(I-1)*PI/(2*N))
!
!                  A CALL OF COSQB FOLLOWED BY A CALL OF
!                  COSQF WILL MULTIPLY THE SEQUENCE X BY 4*N.
!                  THEREFORE COSQF IS THE UNNORMALIZED INVERSE
!                  OF COSQB.
!
!     WSAVE   CONTAINS INITIALIZATION CALCULATIONS WHICH MUST NOT
!             BE DESTROYED BETWEEN CALLS OF COSQB OR COSQF.
!
! **********************************************************************
!
!     SUBROUTINE CFFTI(N,WSAVE)
!
!     SUBROUTINE CFFTI INITIALIZES THE ARRAY WSAVE WHICH IS USED IN
!     BOTH CFFTF AND CFFTB. THE PRIME FACTORIZATION OF N TOGETHER WITH
!     A TABULATION OF THE TRIGONOMETRIC FUNCTIONS ARE COMPUTED AND
!     STORED IN WSAVE.
!
!     INPUT PARAMETER
!
!     N       THE LENGTH OF THE SEQUENCE TO BE TRANSFORMED
!
!     OUTPUT PARAMETER
!
!     WSAVE   A WORK ARRAY WHICH MUST BE DIMENSIONED AT LEAST 4*N+15
!             THE SAME WORK ARRAY CAN BE USED FOR BOTH CFFTF AND CFFTB
!             AS LONG AS N REMAINS UNCHANGED. DIFFERENT WSAVE ARRAYS
!             ARE REQUIRED FOR DIFFERENT VALUES OF N. THE CONTENTS OF
!             WSAVE MUST NOT BE CHANGED BETWEEN CALLS OF CFFTF OR CFFTB.
!
! **********************************************************************
!
!     SUBROUTINE CFFTF(N,C,WSAVE)
!
!     SUBROUTINE CFFTF COMPUTES THE FORWARD COMPLEX DISCRETE FOURIER
!     TRANSFORM (THE FOURIER ANALYSIS). EQUIVALENTLY , CFFTF COMPUTES
!     THE FOURIER COEFFICIENTS OF A COMPLEX PERIODIC SEQUENCE.
!     THE TRANSFORM IS DEFINED BELOW AT OUTPUT PARAMETER C.
!
!     THE TRANSFORM IS NOT NORMALIZED. TO OBTAIN A NORMALIZED TRANSFORM
!     THE OUTPUT MUST BE DIVIDED BY N. OTHERWISE A CALL OF CFFTF
!     FOLLOWED BY A CALL OF CFFTB WILL MULTIPLY THE SEQUENCE BY N.
!
!     THE ARRAY WSAVE WHICH IS USED BY SUBROUTINE CFFTF MUST BE
!     INITIALIZED BY CALLING SUBROUTINE CFFTI(N,WSAVE).
!
!     INPUT PARAMETERS
!
!
!     N      THE LENGTH OF THE COMPLEX SEQUENCE C. THE METHOD IS
!            MORE EFFICIENT WHEN N IS THE PRODUCT OF SMALL PRIMES. N
!
!     C      A COMPLEX ARRAY OF LENGTH N WHICH CONTAINS THE SEQUENCE
!
!     WSAVE   A REAL WORK ARRAY WHICH MUST BE DIMENSIONED AT LEAST 4N+15
!             IN THE PROGRAM THAT CALLS CFFTF. THE WSAVE ARRAY MUST BE
!             INITIALIZED BY CALLING SUBROUTINE CFFTI(N,WSAVE) AND A
!             DIFFERENT WSAVE ARRAY MUST BE USED FOR EACH DIFFERENT
!             VALUE OF N. THIS INITIALIZATION DOES NOT HAVE TO BE
!             REPEATED SO LONG AS N REMAINS UNCHANGED THUS SUBSEQUENT
!             TRANSFORMS CAN BE OBTAINED FASTER THAN THE FIRST.
!             THE SAME WSAVE ARRAY CAN BE USED BY CFFTF AND CFFTB.
!
!     OUTPUT PARAMETERS
!
!     C      FOR J=1,...,N
!
!                C(J)=THE SUM FROM K=1,...,N OF
!
!                      C(K)*EXP(-I*(J-1)*(K-1)*2*PI/N)
!
!                            WHERE I=SQRT(-1)
!
!     WSAVE   CONTAINS INITIALIZATION CALCULATIONS WHICH MUST NOT BE
!             DESTROYED BETWEEN CALLS OF SUBROUTINE CFFTF OR CFFTB
!
! **********************************************************************
!
!     SUBROUTINE CFFTB(N,C,WSAVE)
!
!     SUBROUTINE CFFTB COMPUTES THE BACKWARD COMPLEX DISCRETE FOURIER
!     TRANSFORM (THE FOURIER SYNTHESIS). EQUIVALENTLY , CFFTB COMPUTES
!     A COMPLEX PERIODIC SEQUENCE FROM ITS FOURIER COEFFICIENTS.
!     THE TRANSFORM IS DEFINED BELOW AT OUTPUT PARAMETER C.
!
!     A CALL OF CFFTF FOLLOWED BY A CALL OF CFFTB WILL MULTIPLY THE
!     SEQUENCE BY N.
!
!     THE ARRAY WSAVE WHICH IS USED BY SUBROUTINE CFFTB MUST BE
!     INITIALIZED BY CALLING SUBROUTINE CFFTI(N,WSAVE).
!
!     INPUT PARAMETERS
!
!
!     N      THE LENGTH OF THE COMPLEX SEQUENCE C. THE METHOD IS
!            MORE EFFICIENT WHEN N IS THE PRODUCT OF SMALL PRIMES.
!
!     C      A COMPLEX ARRAY OF LENGTH N WHICH CONTAINS THE SEQUENCE
!
!     WSAVE   A REAL WORK ARRAY WHICH MUST BE DIMENSIONED AT LEAST 4N+15
!             IN THE PROGRAM THAT CALLS CFFTB. THE WSAVE ARRAY MUST BE
!             INITIALIZED BY CALLING SUBROUTINE CFFTI(N,WSAVE) AND A
!             DIFFERENT WSAVE ARRAY MUST BE USED FOR EACH DIFFERENT
!             VALUE OF N. THIS INITIALIZATION DOES NOT HAVE TO BE
!             REPEATED SO LONG AS N REMAINS UNCHANGED THUS SUBSEQUENT
!             TRANSFORMS CAN BE OBTAINED FASTER THAN THE FIRST.
!             THE SAME WSAVE ARRAY CAN BE USED BY CFFTF AND CFFTB.
!
!     OUTPUT PARAMETERS
!
!     C      FOR J=1,...,N
!
!                C(J)=THE SUM FROM K=1,...,N OF
!
!                      C(K)*EXP(I*(J-1)*(K-1)*2*PI/N)
!
!                            WHERE I=SQRT(-1)
!
!     WSAVE   CONTAINS INITIALIZATION CALCULATIONS WHICH MUST NOT BE
!             DESTROYED BETWEEN CALLS OF SUBROUTINE CFFTF OR CFFTB
! **********************************************************************
      SUBROUTINE EZFFTF(N, R, AZERO, A, B, WSAVE)
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER  :: N
      REAL(fish_kind) , INTENT(OUT) :: AZERO
      REAL(fish_kind) , INTENT(IN) :: R(*)
      REAL(fish_kind) , INTENT(OUT) :: A(*)
      REAL(fish_kind) , INTENT(OUT) :: B(*)
      REAL(fish_kind)  :: WSAVE(*)
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: I, NS2, NS2M
      REAL :: CF, CFM
!-----------------------------------------------
!
      IF (N - 2 <= 0) THEN
         IF (N - 2 /= 0) THEN
            AZERO = R(1)
            RETURN 
         ENDIF
         AZERO = 0.5*(R(1)+R(2))
         A(1) = 0.5*(R(1)-R(2))
         RETURN 
      ENDIF
      WSAVE(:N) = R(:N)
      CALL RFFTF (N, WSAVE, WSAVE(N+1))
      CF = 2./FLOAT(N)
      CFM = -CF
      AZERO = 0.5*CF*WSAVE(1)
      NS2 = (N + 1)/2
      NS2M = NS2 - 1
      A(:NS2M) = CF*WSAVE(2:NS2M*2:2)
      B(:NS2M) = CFM*WSAVE(3:NS2M*2+1:2)
      IF (MOD(N,2) == 1) RETURN 
      A(NS2) = 0.5*CF*WSAVE(N)
      B(NS2) = 0.
      RETURN 
      END SUBROUTINE EZFFTF


      SUBROUTINE EZFFTB(N, R, AZERO, A, B, WSAVE)
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER  :: N
      REAL(fish_kind) , INTENT(IN) :: AZERO
      REAL(fish_kind)  :: R(*)
      REAL(fish_kind) , INTENT(IN) :: A(*)
      REAL(fish_kind) , INTENT(IN) :: B(*)
      REAL(fish_kind)  :: WSAVE(*)
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: NS2, I
!-----------------------------------------------
!
      IF (N - 2 <= 0) THEN
         IF (N - 2 /= 0) THEN
            R(1) = AZERO
            RETURN 
         ENDIF
         R(1) = AZERO + A(1)
         R(2) = AZERO - A(1)
         RETURN 
      ENDIF
      NS2 = (N - 1)/2
      R(2:NS2*2:2) = 0.5*A(:NS2)
      R(3:NS2*2+1:2) = -0.5*B(:NS2)
      R(1) = AZERO
      IF (MOD(N,2) == 0) R(N) = A(NS2+1)
      CALL RFFTB (N, R, WSAVE(N+1))
      RETURN 
      END SUBROUTINE EZFFTB


      SUBROUTINE EZFFTI(N, WSAVE)
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER  :: N
      REAL(fish_kind)  :: WSAVE(*)
!-----------------------------------------------
!
      IF (N == 1) RETURN 
      CALL EZFFT1 (N, WSAVE(2*N+1), WSAVE(3*N+1))
      RETURN 
      END SUBROUTINE EZFFTI


      SUBROUTINE EZFFT1(N, WA, IFAC)
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER , INTENT(IN) :: N
      REAL(fish_kind) , INTENT(INOUT) :: IFAC(*)
      REAL(fish_kind) , INTENT(INOUT) :: WA(*)
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER , DIMENSION(4) :: NTRYH
      INTEGER::NL,NF,J,NTRY,NQ,NR,I,IB,IS,NFM1,L1,K1,IP,L2,IDO,IPM,II
      REAL :: TPI, DUM, ARGH, ARG1, CH1, SH1, DCH1, DSH1, CH1H
!-----------------------------------------------
      DATA NTRYH(1), NTRYH(2), NTRYH(3), NTRYH(4)/ 4, 2, 3, 5/ 
      TPI = 8.0*ATAN(1.0)
      NL = N
      NF = 0
      J = 0
  101 CONTINUE
      J = J + 1
      IF (J - 4 <= 0) THEN
         NTRY = NTRYH(J)
      ELSE
         NTRY = NTRY + 2
      ENDIF
  104 CONTINUE
      NQ = NL/NTRY
      NR = NL - NTRY*NQ
      IF (NR /= 0) GO TO 101
      NF = NF + 1
      IFAC(NF+2) = NTRY
      NL = NQ
      IF (NTRY == 2) THEN
         IF (NF /= 1) THEN
            IFAC(NF+2:4:(-1)) = IFAC(NF+1:3:(-1))
            IFAC(3) = 2
         ENDIF
      ENDIF
      IF (NL /= 1) GO TO 104
      IFAC(1) = N
      IFAC(2) = NF
      ARGH = TPI/FLOAT(N)
      IS = 0
      NFM1 = NF - 1
      L1 = 1
      IF (NFM1 == 0) RETURN 
      DO K1 = 1, NFM1
         IP = IFAC(K1+2)
         L2 = L1*IP
         IDO = N/L2
         IPM = IP - 1
         ARG1 = FLOAT(L1)*ARGH
         CH1 = 1.
         SH1 = 0.
         DCH1 = COS(ARG1)
         DSH1 = SIN(ARG1)
         DO J = 1, IPM
            CH1H = DCH1*CH1 - DSH1*SH1
            SH1 = DCH1*SH1 + DSH1*CH1
            CH1 = CH1H
            I = IS + 2
            WA(I-1) = CH1
            WA(I) = SH1
            IF (IDO >= 5) THEN
               DO II = 5, IDO, 2
                  I = I + 2
                  WA(I-1) = CH1*WA(I-3) - SH1*WA(I-2)
                  WA(I) = CH1*WA(I-2) + SH1*WA(I-3)
               END DO
            ENDIF
            IS = IS + IDO
         END DO
         L1 = L2
      END DO
      RETURN 
      END SUBROUTINE EZFFT1


      SUBROUTINE COSTI(N, WSAVE)
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER , INTENT(IN) :: N
      REAL(fish_kind)  :: WSAVE(*)
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: NM1, NP1, NS2, K, KC
      REAL :: PI, DUM, DT, FK
!-----------------------------------------------
!
      PI = 4.0*ATAN(1.0)
      IF (N <= 3) RETURN 
      NM1 = N - 1
      NP1 = N + 1
      NS2 = N/2
      DT = PI/FLOAT(NM1)
      FK = 0.
      DO K = 2, NS2
         KC = NP1 - K
         FK = FK + 1.
         WSAVE(K) = 2.*SIN(FK*DT)
         WSAVE(KC) = 2.*COS(FK*DT)
      END DO
      CALL RFFTI (NM1, WSAVE(N+1))
      RETURN 
      END SUBROUTINE COSTI


      SUBROUTINE COST(N, X, WSAVE)
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER , INTENT(IN) :: N
      REAL(fish_kind)  :: X(*)
      REAL(fish_kind)  :: WSAVE(*)
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: NM1, NP1, NS2, K, KC, MODN, I
      REAL :: X1H, X1P3, TX2, C1, T1, T2, XIM2, XI
!-----------------------------------------------
!
      NM1 = N - 1
      NP1 = N + 1
      NS2 = N/2
      IF (N - 2 >= 0) THEN
         IF (N - 2 <= 0) THEN
            X1H = X(1) + X(2)
            X(2) = X(1) - X(2)
            X(1) = X1H
            RETURN 
         ENDIF
         IF (N <= 3) THEN
            X1P3 = X(1) + X(3)
            TX2 = X(2) + X(2)
            X(2) = X(1) - X(3)
            X(1) = X1P3 + TX2
            X(3) = X1P3 - TX2
            RETURN 
         ENDIF
         C1 = X(1) - X(N)
         X(1) = X(1) + X(N)
         DO K = 2, NS2
            KC = NP1 - K
            T1 = X(K) + X(KC)
            T2 = X(K) - X(KC)
            C1 = C1 + WSAVE(KC)*T2
            T2 = WSAVE(K)*T2
            X(K) = T1 - T2
            X(KC) = T1 + T2
         END DO
         MODN = MOD(N,2)
         IF (MODN /= 0) X(NS2+1) = X(NS2+1) + X(NS2+1)
         CALL RFFTF (NM1, X, WSAVE(N+1))
         XIM2 = X(2)
         X(2) = C1
         DO I = 4, N, 2
            XI = X(I)
            X(I) = X(I-2) - X(I-1)
            X(I-1) = XIM2
            XIM2 = XI
         END DO
         IF (MODN /= 0) X(N) = XIM2
      ENDIF
      RETURN 
      END SUBROUTINE COST


      SUBROUTINE SINTI(N, WSAVE)
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER , INTENT(IN) :: N
      REAL(fish_kind)  :: WSAVE(*)
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: NS2, NP1, K
      REAL :: PI, DUM, DT
!-----------------------------------------------
!
      PI = 4.0*ATAN(1.0)
      IF (N <= 1) RETURN 
      NS2 = N/2
      NP1 = N + 1
      DT = PI/FLOAT(NP1)
      DO K = 1, NS2
         WSAVE(K) = 2.*SIN(K*DT)
      END DO
      CALL RFFTI (NP1, WSAVE(NS2+1))
      RETURN 
      END SUBROUTINE SINTI


      SUBROUTINE SINT(N, X, WSAVE)
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER  :: N
      REAL(fish_kind)  :: X(*)
      REAL(fish_kind)  :: WSAVE(*)
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: NP1, IW1, IW2, IW3
!-----------------------------------------------
!
      NP1 = N + 1
      IW1 = N/2 + 1
      IW2 = IW1 + NP1
      IW3 = IW2 + NP1
      CALL SINT1 (N, X, WSAVE, WSAVE(IW1), WSAVE(IW2), WSAVE(IW3))
      RETURN 
      END SUBROUTINE SINT


      SUBROUTINE SINT1(N, WAR, WAS, XH, X, IFAC)
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER , INTENT(IN) :: N
!      INTEGER  :: IFAC(*)
      REAL(fish_kind)   :: IFAC(*)
      REAL(fish_kind)  :: WAR(*)
      REAL(fish_kind) , INTENT(IN) :: WAS(*)
      REAL(fish_kind)  :: XH(*)
      REAL(fish_kind)  :: X(*)
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: I, NP1, NS2, K, KC, MODN
      REAL :: SQRT3, XHOLD, T1, T2
!-----------------------------------------------
      DATA SQRT3/ 1.73205080756888/ 
      XH(:N) = WAR(:N)
      WAR(:N) = X(:N)
      IF (N - 2 <= 0) THEN
         IF (N - 2 /= 0) THEN
            XH(1) = XH(1) + XH(1)
            GO TO 106
         ENDIF
         XHOLD = SQRT3*(XH(1)+XH(2))
         XH(2) = SQRT3*(XH(1)-XH(2))
         XH(1) = XHOLD
         GO TO 106
      ENDIF
      NP1 = N + 1
      NS2 = N/2
      X(1) = 0.
      DO K = 1, NS2
         KC = NP1 - K
         T1 = XH(K) - XH(KC)
         T2 = WAS(K)*(XH(K)+XH(KC))
         X(K+1) = T1 + T2
         X(KC+1) = T2 - T1
      END DO
      MODN = MOD(N,2)
      IF (MODN /= 0) X(NS2+2) = 4.*XH(NS2+1)
      CALL RFFTF1 (NP1, X, XH, WAR, IFAC)
      XH(1) = 0.5*X(1)
      DO I = 3, N, 2
         XH(I-1) = -X(I)
         XH(I) = XH(I-2) + X(I-1)
      END DO
      IF (MODN == 0) XH(N) = -X(N+1)
  106 CONTINUE
      X(:N) = WAR(:N)
      WAR(:N) = XH(:N)
      RETURN 
      END SUBROUTINE SINT1


      SUBROUTINE COSQI(N, WSAVE)
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER  :: N
      REAL(fish_kind)  :: WSAVE(*)
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: K
      REAL :: PIH, DUM, DT, FK
!-----------------------------------------------
!
      PIH = 2.0*ATAN(1.0)
      DT = PIH/FLOAT(N)
      FK = 0.
      DO K = 1, N
         FK = FK + 1.
         WSAVE(K) = COS(FK*DT)
      END DO
      CALL RFFTI (N, WSAVE(N+1))
      RETURN 
      END SUBROUTINE COSQI


      SUBROUTINE COSQF(N, X, WSAVE)
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER  :: N
      REAL(fish_kind)  :: X(*)
      REAL(fish_kind)  :: WSAVE(*)
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      REAL :: SQRT2, TSQX
!-----------------------------------------------
      DATA SQRT2/ 1.4142135623731/ 
!
      IF (N - 2 >= 0) THEN
         IF (N - 2 > 0) GO TO 103
         TSQX = SQRT2*X(2)
         X(2) = X(1) - TSQX
         X(1) = X(1) + TSQX
      ENDIF
      RETURN 
  103 CONTINUE
      CALL COSQF1 (N, X, WSAVE, WSAVE(N+1))
      RETURN 
      END SUBROUTINE COSQF


      SUBROUTINE COSQF1(N, X, W, XH)
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER  :: N
      REAL(fish_kind)  :: X(*)
      REAL(fish_kind) , INTENT(IN) :: W(*)
      REAL(fish_kind)  :: XH(*)
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: NS2, NP2, K, KC, MODN, I
      REAL :: XIM1
!-----------------------------------------------
      NS2 = (N + 1)/2
      NP2 = N + 2
      DO K = 2, NS2
         KC = NP2 - K
         XH(K) = X(K) + X(KC)
         XH(KC) = X(K) - X(KC)
      END DO
      MODN = MOD(N,2)
      IF (MODN == 0) XH(NS2+1) = X(NS2+1) + X(NS2+1)
      DO K = 2, NS2
         KC = NP2 - K
         X(K) = W(K-1)*XH(KC) + W(KC-1)*XH(K)
         X(KC) = W(K-1)*XH(K) - W(KC-1)*XH(KC)
      END DO
      IF (MODN == 0) X(NS2+1) = W(NS2)*XH(NS2+1)
      CALL RFFTF (N, X, XH)
      DO I = 3, N, 2
         XIM1 = X(I-1) - X(I)
         X(I) = X(I-1) + X(I)
         X(I-1) = XIM1
      END DO
      RETURN 
      END SUBROUTINE COSQF1


      SUBROUTINE COSQB(N, X, WSAVE)
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER  :: N
      REAL(fish_kind)  :: X(*)
      REAL(fish_kind)  :: WSAVE(*)
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      REAL :: TSQRT2, X1
!-----------------------------------------------
      DATA TSQRT2/ 2.82842712474619/ 
!
      IF (N - 2 <= 0) THEN
         IF (N - 2 /= 0) THEN
            X(1) = 4.*X(1)
            RETURN 
         ENDIF
         X1 = 4.*(X(1)+X(2))
         X(2) = TSQRT2*(X(1)-X(2))
         X(1) = X1
         RETURN 
      ENDIF
      CALL COSQB1 (N, X, WSAVE, WSAVE(N+1))
      RETURN 
      END SUBROUTINE COSQB


      SUBROUTINE COSQB1(N, X, W, XH)
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER  :: N
      REAL(fish_kind)  :: X(*)
      REAL(fish_kind) , INTENT(IN) :: W(*)
      REAL(fish_kind)  :: XH(*)
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: NS2, NP2, I, MODN, K, KC
      REAL :: XIM1
!-----------------------------------------------
      NS2 = (N + 1)/2
      NP2 = N + 2
      DO I = 3, N, 2
         XIM1 = X(I-1) + X(I)
         X(I) = X(I) - X(I-1)
         X(I-1) = XIM1
      END DO
      X(1) = X(1) + X(1)
      MODN = MOD(N,2)
      IF (MODN == 0) X(N) = X(N) + X(N)
      CALL RFFTB (N, X, XH)
      DO K = 2, NS2
         KC = NP2 - K
         XH(K) = W(K-1)*X(KC) + W(KC-1)*X(K)
         XH(KC) = W(K-1)*X(K) - W(KC-1)*X(KC)
      END DO
      IF (MODN == 0) X(NS2+1) = W(NS2)*(X(NS2+1)+X(NS2+1))
      DO K = 2, NS2
         KC = NP2 - K
         X(K) = XH(K) + XH(KC)
         X(KC) = XH(K) - XH(KC)
      END DO
      X(1) = X(1) + X(1)
      RETURN 
      END SUBROUTINE COSQB1


      SUBROUTINE SINQI(N, WSAVE)
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER  :: N
      REAL(fish_kind)  :: WSAVE(*)
!-----------------------------------------------
!
      CALL COSQI (N, WSAVE)
      RETURN 
      END SUBROUTINE SINQI


      SUBROUTINE SINQF(N, X, WSAVE)
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER  :: N
      REAL(fish_kind)  :: X(*)
      REAL(fish_kind)  :: WSAVE(*)
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: NS2, K, KC
      REAL :: XHOLD
!-----------------------------------------------
!
      IF (N == 1) RETURN 
      NS2 = N/2
      DO K = 1, NS2
         KC = N - K
         XHOLD = X(K)
         X(K) = X(KC+1)
         X(KC+1) = XHOLD
      END DO
      CALL COSQF (N, X, WSAVE)
      X(2:N:2) = -X(2:N:2)
      RETURN 
      END SUBROUTINE SINQF


      SUBROUTINE SINQB(N, X, WSAVE)
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER  :: N
      REAL(fish_kind)  :: X(*)
      REAL(fish_kind)  :: WSAVE(*)
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: NS2, K, KC
      REAL :: XHOLD
!-----------------------------------------------
!
      IF (N <= 1) THEN
         X(1) = 4.*X(1)
         RETURN 
      ENDIF
      NS2 = N/2
      X(2:N:2) = -X(2:N:2)
      CALL COSQB (N, X, WSAVE)
      DO K = 1, NS2
         KC = N - K
         XHOLD = X(K)
         X(K) = X(KC+1)
         X(KC+1) = XHOLD
      END DO
      RETURN 
      END SUBROUTINE SINQB


      SUBROUTINE CFFTI(N, WSAVE)
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER  :: N
      REAL(fish_kind)  :: WSAVE(*)
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: IW1, IW2
!-----------------------------------------------
!
      IF (N == 1) RETURN 
      IW1 = N + N + 1
      IW2 = IW1 + N + N
      CALL CFFTI1 (N, WSAVE(IW1), WSAVE(IW2))
      RETURN 
      END SUBROUTINE CFFTI


      SUBROUTINE CFFTI1(N, WA, IFAC)
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER , INTENT(IN) :: N
      REAL(fish_kind) , INTENT(INOUT) :: IFAC(*)
      REAL(fish_kind) , INTENT(INOUT) :: WA(*)
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER , DIMENSION(4) :: NTRYH
      INTEGER :: NL, NF, J, NTRY, NQ, NR, I, IB, L1, K1, IP, LD, L2, IDO , IDOT, IPM, I1, II
      REAL :: TPI, DUM, ARGH, FI, ARGLD, ARG
!-----------------------------------------------
      DATA NTRYH(1), NTRYH(2), NTRYH(3), NTRYH(4)/ 3, 4, 2, 5/ 
      NL = N
      NF = 0
      J = 0
  101 CONTINUE
      J = J + 1
      IF (J - 4 <= 0) THEN
         NTRY = NTRYH(J)
      ELSE
         NTRY = NTRY + 2
      ENDIF
  104 CONTINUE
      NQ = NL/NTRY
      NR = NL - NTRY*NQ
      IF (NR /= 0) GO TO 101
      NF = NF + 1
      IFAC(NF+2) = NTRY
      NL = NQ
      IF (NTRY == 2) THEN
         IF (NF /= 1) THEN
            IFAC(NF+2:4:(-1)) = IFAC(NF+1:3:(-1))
            IFAC(3) = 2
         ENDIF
      ENDIF
      IF (NL /= 1) GO TO 104
      IFAC(1) = N
      IFAC(2) = NF
      TPI = 8.0*ATAN(1.0)
      ARGH = TPI/FLOAT(N)
      I = 2
      L1 = 1
      DO K1 = 1, NF
         IP = IFAC(K1+2)
         LD = 0
         L2 = L1*IP
         IDO = N/L2
         IDOT = IDO + IDO + 2
         IPM = IP - 1
         DO J = 1, IPM
            I1 = I
            WA(I-1) = 1.
            WA(I) = 0.
            LD = LD + L1
            FI = 0.
            ARGLD = FLOAT(LD)*ARGH
            DO II = 4, IDOT, 2
               I = I + 2
               FI = FI + 1.
               ARG = FI*ARGLD
               WA(I-1) = COS(ARG)
               WA(I) = SIN(ARG)
            END DO
            IF (IP <= 5) CYCLE 
            WA(I1-1) = WA(I-1)
            WA(I1) = WA(I)
         END DO
         L1 = L2
      END DO
      RETURN 
      END SUBROUTINE CFFTI1


      SUBROUTINE CFFTB(N, C, WSAVE)
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER  :: N
      REAL(fish_kind)  :: C(*)
      REAL(fish_kind)  :: WSAVE(*)
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: IW1, IW2
!-----------------------------------------------
!
      IF (N == 1) RETURN 
      IW1 = N + N + 1
      IW2 = IW1 + N + N
      CALL CFFTB1 (N, C, WSAVE, WSAVE(IW1), WSAVE(IW2))
      RETURN 
      END SUBROUTINE CFFTB


      SUBROUTINE CFFTB1(N, C, CH, WA, IFAC)
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER , INTENT(IN) :: N
      REAL(fish_kind) , INTENT(IN) :: IFAC(*)
      REAL(fish_kind)  :: C(*)
      REAL(fish_kind)  :: CH(*)
      REAL(fish_kind)  :: WA(*)
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER::NF,NA,L1,IW,K1,IP,L2,IDO,IDOT,IDL1,IX2,IX3,IX4,NAC,N2,I
!-----------------------------------------------
      NF = IFAC(2)
      NA = 0
      L1 = 1
      IW = 1
      DO K1 = 1, NF
         IP = IFAC(K1+2)
         L2 = IP*L1
         IDO = N/L2
         IDOT = IDO + IDO
         IDL1 = IDOT*L1
         IF (IP == 4) THEN
            IX2 = IW + IDOT
            IX3 = IX2 + IDOT
            IF (NA == 0) THEN
               CALL PASSB4 (IDOT, L1, C, CH, WA(IW), WA(IX2), WA(IX3))
            ELSE
               CALL PASSB4 (IDOT, L1, CH, C, WA(IW), WA(IX2), WA(IX3))
            ENDIF
            NA = 1 - NA
         ELSE
            IF (IP == 2) THEN
               IF (NA == 0) THEN
                  CALL PASSB2 (IDOT, L1, C, CH, WA(IW))
               ELSE
                  CALL PASSB2 (IDOT, L1, CH, C, WA(IW))
               ENDIF
               NA = 1 - NA
            ELSE
               IF (IP == 3) THEN
                  IX2 = IW + IDOT
                  IF (NA == 0) THEN
                     CALL PASSB3 (IDOT, L1, C, CH, WA(IW), WA(IX2))
                  ELSE
                     CALL PASSB3 (IDOT, L1, CH, C, WA(IW), WA(IX2))
                  ENDIF
                  NA = 1 - NA
               ELSE
                  IF (IP == 5) THEN
                     IX2 = IW + IDOT
                     IX3 = IX2 + IDOT
                     IX4 = IX3 + IDOT
                     IF (NA == 0) THEN
                        CALL PASSB5 (IDOT, L1, C, CH, WA(IW), WA(IX2), WA(IX3), WA(IX4))
                     ELSE
                        CALL PASSB5 (IDOT, L1, CH, C, WA(IW), WA(IX2), WA(IX3), WA(IX4))
                     ENDIF
                     NA = 1 - NA
                  ELSE
                     IF (NA == 0) THEN
                        CALL PASSB (NAC, IDOT, IP, L1, IDL1, C, C, C, CH, CH, WA(IW))
                     ELSE
                        CALL PASSB (NAC, IDOT, IP, L1, IDL1, CH, CH, CH, C, C, WA(IW))
                     ENDIF
                     IF (NAC /= 0) NA = 1 - NA
                  ENDIF
               ENDIF
            ENDIF
         ENDIF
         L1 = L2
         IW = IW + (IP - 1)*IDOT
      END DO
      IF (NA == 0) RETURN 
      N2 = N + N
      C(:N2) = CH(:N2)
      RETURN 
      END SUBROUTINE CFFTB1


      SUBROUTINE PASSB2(IDO, L1, CC, CH, WA1)
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER , INTENT(IN) :: IDO
      INTEGER , INTENT(IN) :: L1
      REAL(fish_kind) , INTENT(IN) :: CC(IDO,2,L1)
      REAL(fish_kind) , INTENT(OUT) :: CH(IDO,L1,2)
      REAL(fish_kind) , INTENT(IN) :: WA1(1)
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: K, I
      REAL :: TR2, TI2
!-----------------------------------------------
      IF (IDO <= 2) THEN
         CH(1,:,1) = CC(1,1,:) + CC(1,2,:)
         CH(1,:,2) = CC(1,1,:) - CC(1,2,:)
         CH(2,:,1) = CC(2,1,:) + CC(2,2,:)
         CH(2,:,2) = CC(2,1,:) - CC(2,2,:)
         RETURN 
      ENDIF
      DO K = 1, L1
         DO I = 2, IDO, 2
            CH(I-1,K,1) = CC(I-1,1,K) + CC(I-1,2,K)
            TR2 = CC(I-1,1,K) - CC(I-1,2,K)
            CH(I,K,1) = CC(I,1,K) + CC(I,2,K)
            TI2 = CC(I,1,K) - CC(I,2,K)
            CH(I,K,2) = WA1(I-1)*TI2 + WA1(I)*TR2
            CH(I-1,K,2) = WA1(I-1)*TR2 - WA1(I)*TI2
         END DO
      END DO
      RETURN 
      END SUBROUTINE PASSB2


      SUBROUTINE PASSB3(IDO, L1, CC, CH, WA1, WA2)
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER , INTENT(IN) :: IDO
      INTEGER , INTENT(IN) :: L1
      REAL(fish_kind) , INTENT(IN) :: CC(IDO,3,L1)
      REAL(fish_kind) , INTENT(OUT) :: CH(IDO,L1,3)
      REAL(fish_kind) , INTENT(IN) :: WA1(*)
      REAL(fish_kind) , INTENT(IN) :: WA2(*)
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: K, I
      REAL::TAUR,TAUI,TR2,CR2,TI2,CI2,CR3,CI3,DR2,DR3,DI2,DI3
!-----------------------------------------------
      DATA TAUR, TAUI/ -.5, 0.866025403784439/ 
      IF (IDO == 2) THEN
         DO K = 1, L1
            TR2 = CC(1,2,K) + CC(1,3,K)
            CR2 = CC(1,1,K) + TAUR*TR2
            CH(1,K,1) = CC(1,1,K) + TR2
            TI2 = CC(2,2,K) + CC(2,3,K)
            CI2 = CC(2,1,K) + TAUR*TI2
            CH(2,K,1) = CC(2,1,K) + TI2
            CR3 = TAUI*(CC(1,2,K)-CC(1,3,K))
            CI3 = TAUI*(CC(2,2,K)-CC(2,3,K))
            CH(1,K,2) = CR2 - CI3
            CH(1,K,3) = CR2 + CI3
            CH(2,K,2) = CI2 + CR3
            CH(2,K,3) = CI2 - CR3
         END DO
         RETURN 
      ENDIF
      DO K = 1, L1
         DO I = 2, IDO, 2
            TR2 = CC(I-1,2,K) + CC(I-1,3,K)
            CR2 = CC(I-1,1,K) + TAUR*TR2
            CH(I-1,K,1) = CC(I-1,1,K) + TR2
            TI2 = CC(I,2,K) + CC(I,3,K)
            CI2 = CC(I,1,K) + TAUR*TI2
            CH(I,K,1) = CC(I,1,K) + TI2
            CR3 = TAUI*(CC(I-1,2,K)-CC(I-1,3,K))
            CI3 = TAUI*(CC(I,2,K)-CC(I,3,K))
            DR2 = CR2 - CI3
            DR3 = CR2 + CI3
            DI2 = CI2 + CR3
            DI3 = CI2 - CR3
            CH(I,K,2) = WA1(I-1)*DI2 + WA1(I)*DR2
            CH(I-1,K,2) = WA1(I-1)*DR2 - WA1(I)*DI2
            CH(I,K,3) = WA2(I-1)*DI3 + WA2(I)*DR3
            CH(I-1,K,3) = WA2(I-1)*DR3 - WA2(I)*DI3
         END DO
      END DO
      RETURN 
      END SUBROUTINE PASSB3


      SUBROUTINE PASSB4(IDO, L1, CC, CH, WA1, WA2, WA3)
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER , INTENT(IN) :: IDO
      INTEGER , INTENT(IN) :: L1
      REAL(fish_kind) , INTENT(IN) :: CC(IDO,4,L1)
      REAL(fish_kind) , INTENT(OUT) :: CH(IDO,L1,4)
      REAL(fish_kind) , INTENT(IN) :: WA1(*)
      REAL(fish_kind) , INTENT(IN) :: WA2(*)
      REAL(fish_kind) , INTENT(IN) :: WA3(*)
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: K, I
      REAL::TI1,TI2,TR4,TI3,TR1,TR2,TI4,TR3,CR3,CI3,CR2,CR4,CI2,CI4
!-----------------------------------------------
      IF (IDO == 2) THEN
         DO K = 1, L1
            TI1 = CC(2,1,K) - CC(2,3,K)
            TI2 = CC(2,1,K) + CC(2,3,K)
            TR4 = CC(2,4,K) - CC(2,2,K)
            TI3 = CC(2,2,K) + CC(2,4,K)
            TR1 = CC(1,1,K) - CC(1,3,K)
            TR2 = CC(1,1,K) + CC(1,3,K)
            TI4 = CC(1,2,K) - CC(1,4,K)
            TR3 = CC(1,2,K) + CC(1,4,K)
            CH(1,K,1) = TR2 + TR3
            CH(1,K,3) = TR2 - TR3
            CH(2,K,1) = TI2 + TI3
            CH(2,K,3) = TI2 - TI3
            CH(1,K,2) = TR1 + TR4
            CH(1,K,4) = TR1 - TR4
            CH(2,K,2) = TI1 + TI4
            CH(2,K,4) = TI1 - TI4
         END DO
         RETURN 
      ENDIF
      DO K = 1, L1
         DO I = 2, IDO, 2
            TI1 = CC(I,1,K) - CC(I,3,K)
            TI2 = CC(I,1,K) + CC(I,3,K)
            TI3 = CC(I,2,K) + CC(I,4,K)
            TR4 = CC(I,4,K) - CC(I,2,K)
            TR1 = CC(I-1,1,K) - CC(I-1,3,K)
            TR2 = CC(I-1,1,K) + CC(I-1,3,K)
            TI4 = CC(I-1,2,K) - CC(I-1,4,K)
            TR3 = CC(I-1,2,K) + CC(I-1,4,K)
            CH(I-1,K,1) = TR2 + TR3
            CR3 = TR2 - TR3
            CH(I,K,1) = TI2 + TI3
            CI3 = TI2 - TI3
            CR2 = TR1 + TR4
            CR4 = TR1 - TR4
            CI2 = TI1 + TI4
            CI4 = TI1 - TI4
            CH(I-1,K,2) = WA1(I-1)*CR2 - WA1(I)*CI2
            CH(I,K,2) = WA1(I-1)*CI2 + WA1(I)*CR2
            CH(I-1,K,3) = WA2(I-1)*CR3 - WA2(I)*CI3
            CH(I,K,3) = WA2(I-1)*CI3 + WA2(I)*CR3
            CH(I-1,K,4) = WA3(I-1)*CR4 - WA3(I)*CI4
            CH(I,K,4) = WA3(I-1)*CI4 + WA3(I)*CR4
         END DO
      END DO
      RETURN 
      END SUBROUTINE PASSB4


      SUBROUTINE PASSB5(IDO, L1, CC, CH, WA1, WA2, WA3, WA4)
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER , INTENT(IN) :: IDO
      INTEGER , INTENT(IN) :: L1
      REAL(fish_kind) , INTENT(IN) :: CC(IDO,5,L1)
      REAL(fish_kind) , INTENT(OUT) :: CH(IDO,L1,5)
      REAL(fish_kind) , INTENT(IN) :: WA1(*)
      REAL(fish_kind) , INTENT(IN) :: WA2(*)
      REAL(fish_kind) , INTENT(IN) :: WA3(*)
      REAL(fish_kind) , INTENT(IN) :: WA4(*)
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: K, I
      REAL :: TR11, TI11, TR12, TI12, TI5, TI2, TI4, TI3, TR5, TR2, TR4
      REAL :: TR3, CR2, CI2, CR3, CI3, CR5, CI5, CR4, CI4, DR3, DR4, DI3
      REAL :: DI4, DR5, DR2, DI5, DI2
!-----------------------------------------------
      DATA TR11, TI11, TR12, TI12/ 0.309016994374947, 0.951056516295154, -.809016994374947, 0.587785252292473/ 
      IF (IDO == 2) THEN
         DO K = 1, L1
            TI5 = CC(2,2,K) - CC(2,5,K)
            TI2 = CC(2,2,K) + CC(2,5,K)
            TI4 = CC(2,3,K) - CC(2,4,K)
            TI3 = CC(2,3,K) + CC(2,4,K)
            TR5 = CC(1,2,K) - CC(1,5,K)
            TR2 = CC(1,2,K) + CC(1,5,K)
            TR4 = CC(1,3,K) - CC(1,4,K)
            TR3 = CC(1,3,K) + CC(1,4,K)
            CH(1,K,1) = CC(1,1,K) + TR2 + TR3
            CH(2,K,1) = CC(2,1,K) + TI2 + TI3
            CR2 = CC(1,1,K) + TR11*TR2 + TR12*TR3
            CI2 = CC(2,1,K) + TR11*TI2 + TR12*TI3
            CR3 = CC(1,1,K) + TR12*TR2 + TR11*TR3
            CI3 = CC(2,1,K) + TR12*TI2 + TR11*TI3
            CR5 = TI11*TR5 + TI12*TR4
            CI5 = TI11*TI5 + TI12*TI4
            CR4 = TI12*TR5 - TI11*TR4
            CI4 = TI12*TI5 - TI11*TI4
            CH(1,K,2) = CR2 - CI5
            CH(1,K,5) = CR2 + CI5
            CH(2,K,2) = CI2 + CR5
            CH(2,K,3) = CI3 + CR4
            CH(1,K,3) = CR3 - CI4
            CH(1,K,4) = CR3 + CI4
            CH(2,K,4) = CI3 - CR4
            CH(2,K,5) = CI2 - CR5
         END DO
         RETURN 
      ENDIF
      DO K = 1, L1
         DO I = 2, IDO, 2
            TI5 = CC(I,2,K) - CC(I,5,K)
            TI2 = CC(I,2,K) + CC(I,5,K)
            TI4 = CC(I,3,K) - CC(I,4,K)
            TI3 = CC(I,3,K) + CC(I,4,K)
            TR5 = CC(I-1,2,K) - CC(I-1,5,K)
            TR2 = CC(I-1,2,K) + CC(I-1,5,K)
            TR4 = CC(I-1,3,K) - CC(I-1,4,K)
            TR3 = CC(I-1,3,K) + CC(I-1,4,K)
            CH(I-1,K,1) = CC(I-1,1,K) + TR2 + TR3
            CH(I,K,1) = CC(I,1,K) + TI2 + TI3
            CR2 = CC(I-1,1,K) + TR11*TR2 + TR12*TR3
            CI2 = CC(I,1,K) + TR11*TI2 + TR12*TI3
            CR3 = CC(I-1,1,K) + TR12*TR2 + TR11*TR3
            CI3 = CC(I,1,K) + TR12*TI2 + TR11*TI3
            CR5 = TI11*TR5 + TI12*TR4
            CI5 = TI11*TI5 + TI12*TI4
            CR4 = TI12*TR5 - TI11*TR4
            CI4 = TI12*TI5 - TI11*TI4
            DR3 = CR3 - CI4
            DR4 = CR3 + CI4
            DI3 = CI3 + CR4
            DI4 = CI3 - CR4
            DR5 = CR2 + CI5
            DR2 = CR2 - CI5
            DI5 = CI2 - CR5
            DI2 = CI2 + CR5
            CH(I-1,K,2) = WA1(I-1)*DR2 - WA1(I)*DI2
            CH(I,K,2) = WA1(I-1)*DI2 + WA1(I)*DR2
            CH(I-1,K,3) = WA2(I-1)*DR3 - WA2(I)*DI3
            CH(I,K,3) = WA2(I-1)*DI3 + WA2(I)*DR3
            CH(I-1,K,4) = WA3(I-1)*DR4 - WA3(I)*DI4
            CH(I,K,4) = WA3(I-1)*DI4 + WA3(I)*DR4
            CH(I-1,K,5) = WA4(I-1)*DR5 - WA4(I)*DI5
            CH(I,K,5) = WA4(I-1)*DI5 + WA4(I)*DR5
         END DO
      END DO
      RETURN 
      END SUBROUTINE PASSB5


      SUBROUTINE PASSB(NAC, IDO, IP, L1, IDL1, CC, C1, C2, CH, CH2, WA)
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER , INTENT(OUT) :: NAC
      INTEGER , INTENT(IN) :: IDO
      INTEGER , INTENT(IN) :: IP
      INTEGER , INTENT(IN) :: L1
      INTEGER , INTENT(IN) :: IDL1
      REAL(fish_kind) , INTENT(IN) :: CC(IDO,IP,L1)
      REAL(fish_kind) , INTENT(OUT) :: C1(IDO,L1,IP)
      REAL(fish_kind) , INTENT(INOUT) :: C2(IDL1,IP)
      REAL(fish_kind) , INTENT(INOUT) :: CH(IDO,L1,IP)
      REAL(fish_kind) , INTENT(INOUT) :: CH2(IDL1,IP)
      REAL(fish_kind) , INTENT(IN) :: WA(*)
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: IDOT, NT, IPP2, IPPH, IDP, J, JC, K, I, IDL, INC, L, LC
      INTEGER :: IK, IDLJ, IDIJ, IDJ
      REAL :: WAR, WAI
!-----------------------------------------------
      IDOT = IDO/2
      NT = IP*IDL1
      IPP2 = IP + 2
      IPPH = (IP + 1)/2
      IDP = IP*IDO
!
      IF (IDO >= L1) THEN
         DO J = 2, IPPH
            JC = IPP2 - J
            CH(:,:,J) = CC(:,J,:) + CC(:,JC,:)
            CH(:,:,JC) = CC(:,J,:) - CC(:,JC,:)
         END DO
         CH(:,:,1) = CC(:,1,:)
      ELSE
         DO J = 2, IPPH
            JC = IPP2 - J
            CH(:,:,J) = CC(:,J,:) + CC(:,JC,:)
            CH(:,:,JC) = CC(:,J,:) - CC(:,JC,:)
         END DO
         CH(:,:,1) = CC(:,1,:)
      ENDIF
      IDL = 2 - IDO
      INC = 0
      DO L = 2, IPPH
         LC = IPP2 - L
         IDL = IDL + IDO
         C2(:,L) = CH2(:,1) + WA(IDL-1)*CH2(:,2)
         C2(:,LC) = WA(IDL)*CH2(:,IP)
         IDLJ = IDL
         INC = INC + IDO
         DO J = 3, IPPH
            JC = IPP2 - J
            IDLJ = IDLJ + INC
            IF (IDLJ > IDP) IDLJ = IDLJ - IDP
            WAR = WA(IDLJ-1)
            WAI = WA(IDLJ)
            C2(:,L) = C2(:,L) + WAR*CH2(:,J)
            C2(:,LC) = C2(:,LC) + WAI*CH2(:,JC)
         END DO
      END DO
      DO J = 2, IPPH
         CH2(:,1) = CH2(:,1) + CH2(:,J)
      END DO
      DO J = 2, IPPH
         JC = IPP2 - J
         CH2(:IDL1-1:2,J) = C2(:IDL1-1:2,J) - C2(2:IDL1:2,JC)
         CH2(:IDL1-1:2,JC) = C2(:IDL1-1:2,J) + C2(2:IDL1:2,JC)
         CH2(2:IDL1:2,J) = C2(2:IDL1:2,J) + C2(:IDL1-1:2,JC)
         CH2(2:IDL1:2,JC) = C2(2:IDL1:2,J) - C2(:IDL1-1:2,JC)
      END DO
      NAC = 1
      IF (IDO == 2) RETURN 
      NAC = 0
      C2(:,1) = CH2(:,1)
      C1(1,:,2:IP) = CH(1,:,2:IP)
      C1(2,:,2:IP) = CH(2,:,2:IP)
      IF (IDOT <= L1) THEN
         IDIJ = 0
         DO J = 2, IP
            IDIJ = IDIJ + 2
            DO I = 4, IDO, 2
               IDIJ = IDIJ + 2
               C1(I-1,:,J) = WA(IDIJ-1)*CH(I-1,:,J) - WA(IDIJ)*CH(I,:,J)
               C1(I,:,J) = WA(IDIJ-1)*CH(I,:,J) + WA(IDIJ)*CH(I-1,:,J)
            END DO
         END DO
         RETURN 
      ENDIF
      IDJ = 2 - IDO
      DO J = 2, IP
         IDJ = IDJ + IDO
         DO K = 1, L1
            IDIJ = IDJ
            C1(3:IDO-1:2,K,J) = WA(IDIJ+1:IDO-3+IDIJ:2)*CH(3:IDO-1:2,K,J) - WA(IDIJ+2:IDO-2+IDIJ:2)*CH(4:IDO:2,K,J)
            C1(4:IDO:2,K,J) = WA(IDIJ+1:IDO-3+IDIJ:2)*CH(4:IDO:2,K,J) +    WA(IDIJ+2:IDO-2+IDIJ:2)*CH(3:IDO-1:2,K,J)
         END DO
      END DO
      RETURN 
      END SUBROUTINE PASSB


      SUBROUTINE CFFTF(N, C, WSAVE)
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER  :: N
      REAL(fish_kind)  :: C(*)
      REAL(fish_kind)  :: WSAVE(*)
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: IW1, IW2
!-----------------------------------------------
!
      IF (N == 1) RETURN 
      IW1 = N + N + 1
      IW2 = IW1 + N + N
      CALL CFFTF1 (N, C, WSAVE, WSAVE(IW1), WSAVE(IW2))
      RETURN 
      END SUBROUTINE CFFTF


      SUBROUTINE CFFTF1(N, C, CH, WA, IFAC)
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER , INTENT(IN) :: N
      REAL(fish_kind) , INTENT(IN) :: IFAC(*)
      REAL(fish_kind)  :: C(*)
      REAL(fish_kind)  :: CH(*)
      REAL(fish_kind)  :: WA(*)
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER::NF,NA,L1,IW,K1,IP,L2,IDO,IDOT,IDL1,IX2,IX3,IX4,NAC,N2,I
!-----------------------------------------------
      NF = IFAC(2)
      NA = 0
      L1 = 1
      IW = 1
      DO K1 = 1, NF
         IP = IFAC(K1+2)
         L2 = IP*L1
         IDO = N/L2
         IDOT = IDO + IDO
         IDL1 = IDOT*L1
         IF (IP == 4) THEN
            IX2 = IW + IDOT
            IX3 = IX2 + IDOT
            IF (NA == 0) THEN
               CALL PASSF4 (IDOT, L1, C, CH, WA(IW), WA(IX2), WA(IX3))
            ELSE
               CALL PASSF4 (IDOT, L1, CH, C, WA(IW), WA(IX2), WA(IX3))
            ENDIF
            NA = 1 - NA
         ELSE
            IF (IP == 2) THEN
               IF (NA == 0) THEN
                  CALL PASSF2 (IDOT, L1, C, CH, WA(IW))
               ELSE
                  CALL PASSF2 (IDOT, L1, CH, C, WA(IW))
               ENDIF
               NA = 1 - NA
            ELSE
               IF (IP == 3) THEN
                  IX2 = IW + IDOT
                  IF (NA == 0) THEN
                     CALL PASSF3 (IDOT, L1, C, CH, WA(IW), WA(IX2))
                  ELSE
                     CALL PASSF3 (IDOT, L1, CH, C, WA(IW), WA(IX2))
                  ENDIF
                  NA = 1 - NA
               ELSE
                  IF (IP == 5) THEN
                     IX2 = IW + IDOT
                     IX3 = IX2 + IDOT
                     IX4 = IX3 + IDOT
                     IF (NA == 0) THEN
                        CALL PASSF5 (IDOT, L1, C, CH, WA(IW), WA(IX2), WA(IX3), WA(IX4))
                     ELSE
                        CALL PASSF5 (IDOT, L1, CH, C, WA(IW), WA(IX2), WA(IX3), WA(IX4))
                     ENDIF
                     NA = 1 - NA
                  ELSE
                     IF (NA == 0) THEN
                        CALL PASSF (NAC, IDOT, IP, L1, IDL1, C, C, C, CH, CH, WA(IW))
                     ELSE
                        CALL PASSF (NAC, IDOT, IP, L1, IDL1, CH, CH, CH, C, C, WA(IW))
                     ENDIF
                     IF (NAC /= 0) NA = 1 - NA
                  ENDIF
               ENDIF
            ENDIF
         ENDIF
         L1 = L2
         IW = IW + (IP - 1)*IDOT
      END DO
      IF (NA == 0) RETURN 
      N2 = N + N
      C(:N2) = CH(:N2)
      RETURN 
      END SUBROUTINE CFFTF1


      SUBROUTINE PASSF2(IDO, L1, CC, CH, WA1)
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER , INTENT(IN) :: IDO
      INTEGER , INTENT(IN) :: L1
      REAL(fish_kind) , INTENT(IN) :: CC(IDO,2,L1)
      REAL(fish_kind) , INTENT(OUT) :: CH(IDO,L1,2)
      REAL(fish_kind) , INTENT(IN) :: WA1(*)
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: K, I
      REAL :: TR2, TI2
!-----------------------------------------------
      IF (IDO <= 2) THEN
         CH(1,:,1) = CC(1,1,:) + CC(1,2,:)
         CH(1,:,2) = CC(1,1,:) - CC(1,2,:)
         CH(2,:,1) = CC(2,1,:) + CC(2,2,:)
         CH(2,:,2) = CC(2,1,:) - CC(2,2,:)
         RETURN 
      ENDIF
      DO K = 1, L1
         DO I = 2, IDO, 2
            CH(I-1,K,1) = CC(I-1,1,K) + CC(I-1,2,K)
            TR2 = CC(I-1,1,K) - CC(I-1,2,K)
            CH(I,K,1) = CC(I,1,K) + CC(I,2,K)
            TI2 = CC(I,1,K) - CC(I,2,K)
            CH(I,K,2) = WA1(I-1)*TI2 - WA1(I)*TR2
            CH(I-1,K,2) = WA1(I-1)*TR2 + WA1(I)*TI2
         END DO
      END DO
      RETURN 
      END SUBROUTINE PASSF2


      SUBROUTINE PASSF3(IDO, L1, CC, CH, WA1, WA2)
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER , INTENT(IN) :: IDO
      INTEGER , INTENT(IN) :: L1
      REAL(fish_kind) , INTENT(IN) :: CC(IDO,3,L1)
      REAL(fish_kind) , INTENT(OUT) :: CH(IDO,L1,3)
      REAL(fish_kind) , INTENT(IN) :: WA1(*)
      REAL(fish_kind) , INTENT(IN) :: WA2(*)
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: K, I
      REAL::TAUR,TAUI,TR2,CR2,TI2,CI2,CR3,CI3,DR2,DR3,DI2,DI3
!-----------------------------------------------
      DATA TAUR, TAUI/ -.5,  - 0.866025403784439/ 
      IF (IDO == 2) THEN
         DO K = 1, L1
            TR2 = CC(1,2,K) + CC(1,3,K)
            CR2 = CC(1,1,K) + TAUR*TR2
            CH(1,K,1) = CC(1,1,K) + TR2
            TI2 = CC(2,2,K) + CC(2,3,K)
            CI2 = CC(2,1,K) + TAUR*TI2
            CH(2,K,1) = CC(2,1,K) + TI2
            CR3 = TAUI*(CC(1,2,K)-CC(1,3,K))
            CI3 = TAUI*(CC(2,2,K)-CC(2,3,K))
            CH(1,K,2) = CR2 - CI3
            CH(1,K,3) = CR2 + CI3
            CH(2,K,2) = CI2 + CR3
            CH(2,K,3) = CI2 - CR3
         END DO
         RETURN 
      ENDIF
      DO K = 1, L1
         DO I = 2, IDO, 2
            TR2 = CC(I-1,2,K) + CC(I-1,3,K)
            CR2 = CC(I-1,1,K) + TAUR*TR2
            CH(I-1,K,1) = CC(I-1,1,K) + TR2
            TI2 = CC(I,2,K) + CC(I,3,K)
            CI2 = CC(I,1,K) + TAUR*TI2
            CH(I,K,1) = CC(I,1,K) + TI2
            CR3 = TAUI*(CC(I-1,2,K)-CC(I-1,3,K))
            CI3 = TAUI*(CC(I,2,K)-CC(I,3,K))
            DR2 = CR2 - CI3
            DR3 = CR2 + CI3
            DI2 = CI2 + CR3
            DI3 = CI2 - CR3
            CH(I,K,2) = WA1(I-1)*DI2 - WA1(I)*DR2
            CH(I-1,K,2) = WA1(I-1)*DR2 + WA1(I)*DI2
            CH(I,K,3) = WA2(I-1)*DI3 - WA2(I)*DR3
            CH(I-1,K,3) = WA2(I-1)*DR3 + WA2(I)*DI3
         END DO
      END DO
      RETURN 
      END SUBROUTINE PASSF3


      SUBROUTINE PASSF4(IDO, L1, CC, CH, WA1, WA2, WA3)
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER , INTENT(IN) :: IDO
      INTEGER , INTENT(IN) :: L1
      REAL(fish_kind) , INTENT(IN) :: CC(IDO,4,L1)
      REAL(fish_kind) , INTENT(OUT) :: CH(IDO,L1,4)
      REAL(fish_kind) , INTENT(IN) :: WA1(*)
      REAL(fish_kind) , INTENT(IN) :: WA2(*)
      REAL(fish_kind) , INTENT(IN) :: WA3(*)
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: K, I
      REAL::TI1,TI2,TR4,TI3,TR1,TR2,TI4,TR3,CR3,CI3,CR2,CR4,CI2,CI4
!-----------------------------------------------
      IF (IDO == 2) THEN
         DO K = 1, L1
            TI1 = CC(2,1,K) - CC(2,3,K)
            TI2 = CC(2,1,K) + CC(2,3,K)
            TR4 = CC(2,2,K) - CC(2,4,K)
            TI3 = CC(2,2,K) + CC(2,4,K)
            TR1 = CC(1,1,K) - CC(1,3,K)
            TR2 = CC(1,1,K) + CC(1,3,K)
            TI4 = CC(1,4,K) - CC(1,2,K)
            TR3 = CC(1,2,K) + CC(1,4,K)
            CH(1,K,1) = TR2 + TR3
            CH(1,K,3) = TR2 - TR3
            CH(2,K,1) = TI2 + TI3
            CH(2,K,3) = TI2 - TI3
            CH(1,K,2) = TR1 + TR4
            CH(1,K,4) = TR1 - TR4
            CH(2,K,2) = TI1 + TI4
            CH(2,K,4) = TI1 - TI4
         END DO
         RETURN 
      ENDIF
      DO K = 1, L1
         DO I = 2, IDO, 2
            TI1 = CC(I,1,K) - CC(I,3,K)
            TI2 = CC(I,1,K) + CC(I,3,K)
            TI3 = CC(I,2,K) + CC(I,4,K)
            TR4 = CC(I,2,K) - CC(I,4,K)
            TR1 = CC(I-1,1,K) - CC(I-1,3,K)
            TR2 = CC(I-1,1,K) + CC(I-1,3,K)
            TI4 = CC(I-1,4,K) - CC(I-1,2,K)
            TR3 = CC(I-1,2,K) + CC(I-1,4,K)
            CH(I-1,K,1) = TR2 + TR3
            CR3 = TR2 - TR3
            CH(I,K,1) = TI2 + TI3
            CI3 = TI2 - TI3
            CR2 = TR1 + TR4
            CR4 = TR1 - TR4
            CI2 = TI1 + TI4
            CI4 = TI1 - TI4
            CH(I-1,K,2) = WA1(I-1)*CR2 + WA1(I)*CI2
            CH(I,K,2) = WA1(I-1)*CI2 - WA1(I)*CR2
            CH(I-1,K,3) = WA2(I-1)*CR3 + WA2(I)*CI3
            CH(I,K,3) = WA2(I-1)*CI3 - WA2(I)*CR3
            CH(I-1,K,4) = WA3(I-1)*CR4 + WA3(I)*CI4
            CH(I,K,4) = WA3(I-1)*CI4 - WA3(I)*CR4
         END DO
      END DO
      RETURN 
      END SUBROUTINE PASSF4


      SUBROUTINE PASSF5(IDO, L1, CC, CH, WA1, WA2, WA3, WA4)
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER , INTENT(IN) :: IDO
      INTEGER , INTENT(IN) :: L1
      REAL(fish_kind) , INTENT(IN) :: CC(IDO,5,L1)
      REAL(fish_kind) , INTENT(OUT) :: CH(IDO,L1,5)
      REAL(fish_kind) , INTENT(IN) :: WA1(*)
      REAL(fish_kind) , INTENT(IN) :: WA2(*)
      REAL(fish_kind) , INTENT(IN) :: WA3(*)
      REAL(fish_kind) , INTENT(IN) :: WA4(*)
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: K, I
      REAL :: TR11, TI11, TR12, TI12, TI5, TI2, TI4, TI3, TR5, TR2, TR4
      REAL :: TR3, CR2, CI2, CR3, CI3, CR5, CI5, CR4, CI4, DR3, DR4, DI3
      REAL :: DI4, DR5, DR2, DI5, DI2
!-----------------------------------------------
      DATA TR11, TI11, TR12, TI12/ 0.309016994374947, -.951056516295154, -.809016994374947,  - 0.587785252292473/ 
      IF (IDO == 2) THEN
         DO K = 1, L1
            TI5 = CC(2,2,K) - CC(2,5,K)
            TI2 = CC(2,2,K) + CC(2,5,K)
            TI4 = CC(2,3,K) - CC(2,4,K)
            TI3 = CC(2,3,K) + CC(2,4,K)
            TR5 = CC(1,2,K) - CC(1,5,K)
            TR2 = CC(1,2,K) + CC(1,5,K)
            TR4 = CC(1,3,K) - CC(1,4,K)
            TR3 = CC(1,3,K) + CC(1,4,K)
            CH(1,K,1) = CC(1,1,K) + TR2 + TR3
            CH(2,K,1) = CC(2,1,K) + TI2 + TI3
            CR2 = CC(1,1,K) + TR11*TR2 + TR12*TR3
            CI2 = CC(2,1,K) + TR11*TI2 + TR12*TI3
            CR3 = CC(1,1,K) + TR12*TR2 + TR11*TR3
            CI3 = CC(2,1,K) + TR12*TI2 + TR11*TI3
            CR5 = TI11*TR5 + TI12*TR4
            CI5 = TI11*TI5 + TI12*TI4
            CR4 = TI12*TR5 - TI11*TR4
            CI4 = TI12*TI5 - TI11*TI4
            CH(1,K,2) = CR2 - CI5
            CH(1,K,5) = CR2 + CI5
            CH(2,K,2) = CI2 + CR5
            CH(2,K,3) = CI3 + CR4
            CH(1,K,3) = CR3 - CI4
            CH(1,K,4) = CR3 + CI4
            CH(2,K,4) = CI3 - CR4
            CH(2,K,5) = CI2 - CR5
         END DO
         RETURN 
      ENDIF
      DO K = 1, L1
         DO I = 2, IDO, 2
            TI5 = CC(I,2,K) - CC(I,5,K)
            TI2 = CC(I,2,K) + CC(I,5,K)
            TI4 = CC(I,3,K) - CC(I,4,K)
            TI3 = CC(I,3,K) + CC(I,4,K)
            TR5 = CC(I-1,2,K) - CC(I-1,5,K)
            TR2 = CC(I-1,2,K) + CC(I-1,5,K)
            TR4 = CC(I-1,3,K) - CC(I-1,4,K)
            TR3 = CC(I-1,3,K) + CC(I-1,4,K)
            CH(I-1,K,1) = CC(I-1,1,K) + TR2 + TR3
            CH(I,K,1) = CC(I,1,K) + TI2 + TI3
            CR2 = CC(I-1,1,K) + TR11*TR2 + TR12*TR3
            CI2 = CC(I,1,K) + TR11*TI2 + TR12*TI3
            CR3 = CC(I-1,1,K) + TR12*TR2 + TR11*TR3
            CI3 = CC(I,1,K) + TR12*TI2 + TR11*TI3
            CR5 = TI11*TR5 + TI12*TR4
            CI5 = TI11*TI5 + TI12*TI4
            CR4 = TI12*TR5 - TI11*TR4
            CI4 = TI12*TI5 - TI11*TI4
            DR3 = CR3 - CI4
            DR4 = CR3 + CI4
            DI3 = CI3 + CR4
            DI4 = CI3 - CR4
            DR5 = CR2 + CI5
            DR2 = CR2 - CI5
            DI5 = CI2 - CR5
            DI2 = CI2 + CR5
            CH(I-1,K,2) = WA1(I-1)*DR2 + WA1(I)*DI2
            CH(I,K,2) = WA1(I-1)*DI2 - WA1(I)*DR2
            CH(I-1,K,3) = WA2(I-1)*DR3 + WA2(I)*DI3
            CH(I,K,3) = WA2(I-1)*DI3 - WA2(I)*DR3
            CH(I-1,K,4) = WA3(I-1)*DR4 + WA3(I)*DI4
            CH(I,K,4) = WA3(I-1)*DI4 - WA3(I)*DR4
            CH(I-1,K,5) = WA4(I-1)*DR5 + WA4(I)*DI5
            CH(I,K,5) = WA4(I-1)*DI5 - WA4(I)*DR5
         END DO
      END DO
      RETURN 
      END SUBROUTINE PASSF5


      SUBROUTINE PASSF(NAC, IDO, IP, L1, IDL1, CC, C1, C2, CH, CH2, WA)
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER , INTENT(OUT) :: NAC
      INTEGER , INTENT(IN) :: IDO
      INTEGER , INTENT(IN) :: IP
      INTEGER , INTENT(IN) :: L1
      INTEGER , INTENT(IN) :: IDL1
      REAL(fish_kind) , INTENT(IN) :: CC(IDO,IP,L1)
      REAL(fish_kind) , INTENT(OUT) :: C1(IDO,L1,IP)
      REAL(fish_kind) , INTENT(INOUT) :: C2(IDL1,IP)
      REAL(fish_kind) , INTENT(INOUT) :: CH(IDO,L1,IP)
      REAL(fish_kind) , INTENT(INOUT) :: CH2(IDL1,IP)
      REAL(fish_kind) , INTENT(IN) :: WA(*)
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: IDOT, NT, IPP2, IPPH, IDP, J, JC, K, I, IDL, INC, L, LC
      INTEGER :: IK, IDLJ, IDIJ, IDJ
      REAL :: WAR, WAI
!-----------------------------------------------
      IDOT = IDO/2
      NT = IP*IDL1
      IPP2 = IP + 2
      IPPH = (IP + 1)/2
      IDP = IP*IDO
!
      IF (IDO >= L1) THEN
         DO J = 2, IPPH
            JC = IPP2 - J
            CH(:,:,J) = CC(:,J,:) + CC(:,JC,:)
            CH(:,:,JC) = CC(:,J,:) - CC(:,JC,:)
         END DO
         CH(:,:,1) = CC(:,1,:)
      ELSE
         DO J = 2, IPPH
            JC = IPP2 - J
            CH(:,:,J) = CC(:,J,:) + CC(:,JC,:)
            CH(:,:,JC) = CC(:,J,:) - CC(:,JC,:)
         END DO
         CH(:,:,1) = CC(:,1,:)
      ENDIF
      IDL = 2 - IDO
      INC = 0
      DO L = 2, IPPH
         LC = IPP2 - L
         IDL = IDL + IDO
         C2(:,L) = CH2(:,1) + WA(IDL-1)*CH2(:,2)
         C2(:,LC) = -WA(IDL)*CH2(:,IP)
         IDLJ = IDL
         INC = INC + IDO
         DO J = 3, IPPH
            JC = IPP2 - J
            IDLJ = IDLJ + INC
            IF (IDLJ > IDP) IDLJ = IDLJ - IDP
            WAR = WA(IDLJ-1)
            WAI = WA(IDLJ)
            C2(:,L) = C2(:,L) + WAR*CH2(:,J)
            C2(:,LC) = C2(:,LC) - WAI*CH2(:,JC)
         END DO
      END DO
      DO J = 2, IPPH
         CH2(:,1) = CH2(:,1) + CH2(:,J)
      END DO
      DO J = 2, IPPH
         JC = IPP2 - J
         CH2(:IDL1-1:2,J) = C2(:IDL1-1:2,J) - C2(2:IDL1:2,JC)
         CH2(:IDL1-1:2,JC) = C2(:IDL1-1:2,J) + C2(2:IDL1:2,JC)
         CH2(2:IDL1:2,J) = C2(2:IDL1:2,J) + C2(:IDL1-1:2,JC)
         CH2(2:IDL1:2,JC) = C2(2:IDL1:2,J) - C2(:IDL1-1:2,JC)
      END DO
      NAC = 1
      IF (IDO == 2) RETURN 
      NAC = 0
      C2(:,1) = CH2(:,1)
      C1(1,:,2:IP) = CH(1,:,2:IP)
      C1(2,:,2:IP) = CH(2,:,2:IP)
      IF (IDOT <= L1) THEN
         IDIJ = 0
         DO J = 2, IP
            IDIJ = IDIJ + 2
            DO I = 4, IDO, 2
               IDIJ = IDIJ + 2
               C1(I-1,:,J) = WA(IDIJ-1)*CH(I-1,:,J) + WA(IDIJ)*CH(I,:,J)
               C1(I,:,J) = WA(IDIJ-1)*CH(I,:,J) - WA(IDIJ)*CH(I-1,:,J)
            END DO
         END DO
         RETURN 
      ENDIF
      IDJ = 2 - IDO
      DO J = 2, IP
         IDJ = IDJ + IDO
         DO K = 1, L1
            IDIJ = IDJ
            C1(3:IDO-1:2,K,J) = WA(IDIJ+1:IDO-3+IDIJ:2)*CH(3:IDO-1:2,K,J) + WA(IDIJ+2:IDO-2+IDIJ:2)*CH(4:IDO:2,K,J)
            C1(4:IDO:2,K,J) = WA(IDIJ+1:IDO-3+IDIJ:2)*CH(4:IDO:2,K,J) -  WA(IDIJ+2:IDO-2+IDIJ:2)*CH(3:IDO-1:2,K,J)
         END DO
      END DO
      RETURN 
      END SUBROUTINE PASSF


      SUBROUTINE RFFTI(N, WSAVE)
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER  :: N
      REAL(fish_kind)  :: WSAVE(*)
!-----------------------------------------------
!
      IF (N == 1) RETURN 
      CALL RFFTI1 (N, WSAVE(N+1), WSAVE(2*N+1))
      RETURN 
      END SUBROUTINE RFFTI


      SUBROUTINE RFFTI1(N, WA, IFAC)
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER , INTENT(IN) :: N
      REAL(fish_kind) , INTENT(INOUT) :: IFAC(*)
      REAL(fish_kind) , INTENT(OUT) :: WA(*)
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER , DIMENSION(4) :: NTRYH
      INTEGER :: NL, NF, J, NTRY, NQ, NR, I, IB, IS, NFM1, L1, K1, IP
      INTEGER :: LD, L2, IDO, IPM, II
      REAL :: TPI, DUM, ARGH, ARGLD, FI, ARG
!-----------------------------------------------
      DATA NTRYH(1), NTRYH(2), NTRYH(3), NTRYH(4)/ 4, 2, 3, 5/ 
      NL = N
      NF = 0
      J = 0
  101 CONTINUE
      J = J + 1
      IF (J - 4 <= 0) THEN
         NTRY = NTRYH(J)
      ELSE
         NTRY = NTRY + 2
      ENDIF
  104 CONTINUE
      NQ = NL/NTRY
      NR = NL - NTRY*NQ
      IF (NR /= 0) GO TO 101
      NF = NF + 1
      IFAC(NF+2) = NTRY
      NL = NQ
      IF (NTRY == 2) THEN
         IF (NF /= 1) THEN
            IFAC(NF+2:4:(-1)) = IFAC(NF+1:3:(-1))
            IFAC(3) = 2
         ENDIF
      ENDIF
      IF (NL /= 1) GO TO 104
      IFAC(1) = N
      IFAC(2) = NF
      TPI = 8.0*ATAN(1.0)
      ARGH = TPI/FLOAT(N)
      IS = 0
      NFM1 = NF - 1
      L1 = 1
      IF (NFM1 == 0) RETURN 
      DO K1 = 1, NFM1
         IP = IFAC(K1+2)
         LD = 0
         L2 = L1*IP
         IDO = N/L2
         IPM = IP - 1
         DO J = 1, IPM
            LD = LD + L1
            I = IS
            ARGLD = FLOAT(LD)*ARGH
            FI = 0.
            DO II = 3, IDO, 2
               I = I + 2
               FI = FI + 1.
               ARG = FI*ARGLD
               WA(I-1) = COS(ARG)
               WA(I) = SIN(ARG)
            END DO
            IS = IS + IDO
         END DO
         L1 = L2
      END DO
      RETURN 
      END SUBROUTINE RFFTI1


      SUBROUTINE RFFTB(N, R, WSAVE)
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER  :: N
      REAL(fish_kind)  :: R(*)
      REAL(fish_kind)  :: WSAVE(*)
!-----------------------------------------------
!
      IF (N == 1) RETURN 
      CALL RFFTB1 (N, R, WSAVE, WSAVE(N+1), WSAVE(2*N+1))
      RETURN 
      END SUBROUTINE RFFTB


      SUBROUTINE RFFTB1(N, C, CH, WA, IFAC)
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER , INTENT(IN) :: N
      REAL(fish_kind) , INTENT(IN) :: IFAC(*)
      REAL(fish_kind)  :: C(*)
      REAL(fish_kind)  :: CH(*)
      REAL(fish_kind)  :: WA(*)
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: NF, NA, L1, IW, K1, IP, L2, IDO, IDL1, IX2, IX3, IX4, I
!-----------------------------------------------
      NF = IFAC(2)
      NA = 0
      L1 = 1
      IW = 1
      DO K1 = 1, NF
         IP = IFAC(K1+2)
         L2 = IP*L1
         IDO = N/L2
         IDL1 = IDO*L1
         IF (IP == 4) THEN
            IX2 = IW + IDO
            IX3 = IX2 + IDO
            IF (NA == 0) THEN
               CALL RADB4 (IDO, L1, C, CH, WA(IW), WA(IX2), WA(IX3))
            ELSE
               CALL RADB4 (IDO, L1, CH, C, WA(IW), WA(IX2), WA(IX3))
            ENDIF
            NA = 1 - NA
         ELSE
            IF (IP == 2) THEN
               IF (NA == 0) THEN
                  CALL RADB2 (IDO, L1, C, CH, WA(IW))
               ELSE
                  CALL RADB2 (IDO, L1, CH, C, WA(IW))
               ENDIF
               NA = 1 - NA
            ELSE
               IF (IP == 3) THEN
                  IX2 = IW + IDO
                  IF (NA == 0) THEN
                     CALL RADB3 (IDO, L1, C, CH, WA(IW), WA(IX2))
                  ELSE
                     CALL RADB3 (IDO, L1, CH, C, WA(IW), WA(IX2))
                  ENDIF
                  NA = 1 - NA
               ELSE
                  IF (IP == 5) THEN
                     IX2 = IW + IDO
                     IX3 = IX2 + IDO
                     IX4 = IX3 + IDO
                     IF (NA == 0) THEN
                        CALL RADB5 (IDO, L1, C, CH, WA(IW), WA(IX2), WA(IX3), WA(IX4))
                     ELSE
                        CALL RADB5 (IDO, L1, CH, C, WA(IW), WA(IX2), WA(IX3), WA(IX4))
                     ENDIF
                     NA = 1 - NA
                  ELSE
                     IF (NA == 0) THEN
                        CALL RADBG(IDO,IP,L1,IDL1,C,C,C,CH,CH,WA(IW))
                     ELSE
                        CALL RADBG(IDO,IP,L1,IDL1,CH,CH,CH,C,C,WA(IW))
                     ENDIF
                     IF (IDO == 1) NA = 1 - NA
                  ENDIF
               ENDIF
            ENDIF
         ENDIF
         L1 = L2
         IW = IW + (IP - 1)*IDO
      END DO
      IF (NA == 0) RETURN 
      C(:N) = CH(:N)
      RETURN 
      END SUBROUTINE RFFTB1


      SUBROUTINE RADB2(IDO, L1, CC, CH, WA1)
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER , INTENT(IN) :: IDO
      INTEGER , INTENT(IN) :: L1
      REAL(fish_kind) , INTENT(IN) :: CC(IDO,2,L1)
      REAL(fish_kind) , INTENT(OUT) :: CH(IDO,L1,2)
      REAL(fish_kind) , INTENT(IN) :: WA1(*)
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: K, IDP2, I, IC
      REAL :: TR2, TI2
!-----------------------------------------------
      CH(1,:,1) = CC(1,1,:) + CC(IDO,2,:)
      CH(1,:,2) = CC(1,1,:) - CC(IDO,2,:)
      IF (IDO - 2 >= 0) THEN
         IF (IDO - 2 /= 0) THEN
            IDP2 = IDO + 2
            DO K = 1, L1
               DO I = 3, IDO, 2
                  IC = IDP2 - I
                  CH(I-1,K,1) = CC(I-1,1,K) + CC(IC-1,2,K)
                  TR2 = CC(I-1,1,K) - CC(IC-1,2,K)
                  CH(I,K,1) = CC(I,1,K) - CC(IC,2,K)
                  TI2 = CC(I,1,K) + CC(IC,2,K)
                  CH(I-1,K,2) = WA1(I-2)*TR2 - WA1(I-1)*TI2
                  CH(I,K,2) = WA1(I-2)*TI2 + WA1(I-1)*TR2
               END DO
            END DO
            IF (MOD(IDO,2) == 1) RETURN 
         ENDIF
         CH(IDO,:,1) = CC(IDO,1,:) + CC(IDO,1,:)
         CH(IDO,:,2) = -(CC(1,2,:)+CC(1,2,:))
      ENDIF
      RETURN 
      END SUBROUTINE RADB2


      SUBROUTINE RADB3(IDO, L1, CC, CH, WA1, WA2)
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER , INTENT(IN) :: IDO
      INTEGER , INTENT(IN) :: L1
      REAL(fish_kind) , INTENT(IN) :: CC(IDO,3,L1)
      REAL(fish_kind) , INTENT(OUT) :: CH(IDO,L1,3)
      REAL(fish_kind) , INTENT(IN) :: WA1(*)
      REAL(fish_kind) , INTENT(IN) :: WA2(*)
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: K, IDP2, I, IC
      REAL::TAUR,TAUI,TR2,CR2,CI3,TI2,CI2,CR3,DR2,DR3,DI2,DI3
!-----------------------------------------------
      DATA TAUR, TAUI/ -.5, 0.866025403784439/ 
      DO K = 1, L1
         TR2 = CC(IDO,2,K) + CC(IDO,2,K)
         CR2 = CC(1,1,K) + TAUR*TR2
         CH(1,K,1) = CC(1,1,K) + TR2
         CI3 = TAUI*(CC(1,3,K)+CC(1,3,K))
         CH(1,K,2) = CR2 - CI3
         CH(1,K,3) = CR2 + CI3
      END DO
      IF (IDO == 1) RETURN 
      IDP2 = IDO + 2
      DO K = 1, L1
         DO I = 3, IDO, 2
            IC = IDP2 - I
            TR2 = CC(I-1,3,K) + CC(IC-1,2,K)
            CR2 = CC(I-1,1,K) + TAUR*TR2
            CH(I-1,K,1) = CC(I-1,1,K) + TR2
            TI2 = CC(I,3,K) - CC(IC,2,K)
            CI2 = CC(I,1,K) + TAUR*TI2
            CH(I,K,1) = CC(I,1,K) + TI2
            CR3 = TAUI*(CC(I-1,3,K)-CC(IC-1,2,K))
            CI3 = TAUI*(CC(I,3,K)+CC(IC,2,K))
            DR2 = CR2 - CI3
            DR3 = CR2 + CI3
            DI2 = CI2 + CR3
            DI3 = CI2 - CR3
            CH(I-1,K,2) = WA1(I-2)*DR2 - WA1(I-1)*DI2
            CH(I,K,2) = WA1(I-2)*DI2 + WA1(I-1)*DR2
            CH(I-1,K,3) = WA2(I-2)*DR3 - WA2(I-1)*DI3
            CH(I,K,3) = WA2(I-2)*DI3 + WA2(I-1)*DR3
         END DO
      END DO
      RETURN 
      END SUBROUTINE RADB3


      SUBROUTINE RADB4(IDO, L1, CC, CH, WA1, WA2, WA3)
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER , INTENT(IN) :: IDO
      INTEGER , INTENT(IN) :: L1
      REAL(fish_kind) , INTENT(IN) :: CC(IDO,4,L1)
      REAL(fish_kind) , INTENT(OUT) :: CH(IDO,L1,4)
      REAL(fish_kind) , INTENT(IN) :: WA1(*)
      REAL(fish_kind) , INTENT(IN) :: WA2(*)
      REAL(fish_kind) , INTENT(IN) :: WA3(*)
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: K, IDP2, I, IC
      REAL :: SQRT2, TR1, TR2, TR3, TR4, TI1, TI2, TI3, TI4, CR3, CI3 
      REAL :: CR2, CR4, CI2, CI4
!-----------------------------------------------
      DATA SQRT2/ 1.414213562373095/ 
      DO K = 1, L1
         TR1 = CC(1,1,K) - CC(IDO,4,K)
         TR2 = CC(1,1,K) + CC(IDO,4,K)
         TR3 = CC(IDO,2,K) + CC(IDO,2,K)
         TR4 = CC(1,3,K) + CC(1,3,K)
         CH(1,K,1) = TR2 + TR3
         CH(1,K,2) = TR1 - TR4
         CH(1,K,3) = TR2 - TR3
         CH(1,K,4) = TR1 + TR4
      END DO
      IF (IDO - 2 >= 0) THEN
         IF (IDO - 2 /= 0) THEN
            IDP2 = IDO + 2
            DO K = 1, L1
               DO I = 3, IDO, 2
                  IC = IDP2 - I
                  TI1 = CC(I,1,K) + CC(IC,4,K)
                  TI2 = CC(I,1,K) - CC(IC,4,K)
                  TI3 = CC(I,3,K) - CC(IC,2,K)
                  TR4 = CC(I,3,K) + CC(IC,2,K)
                  TR1 = CC(I-1,1,K) - CC(IC-1,4,K)
                  TR2 = CC(I-1,1,K) + CC(IC-1,4,K)
                  TI4 = CC(I-1,3,K) - CC(IC-1,2,K)
                  TR3 = CC(I-1,3,K) + CC(IC-1,2,K)
                  CH(I-1,K,1) = TR2 + TR3
                  CR3 = TR2 - TR3
                  CH(I,K,1) = TI2 + TI3
                  CI3 = TI2 - TI3
                  CR2 = TR1 - TR4
                  CR4 = TR1 + TR4
                  CI2 = TI1 + TI4
                  CI4 = TI1 - TI4
                  CH(I-1,K,2) = WA1(I-2)*CR2 - WA1(I-1)*CI2
                  CH(I,K,2) = WA1(I-2)*CI2 + WA1(I-1)*CR2
                  CH(I-1,K,3) = WA2(I-2)*CR3 - WA2(I-1)*CI3
                  CH(I,K,3) = WA2(I-2)*CI3 + WA2(I-1)*CR3
                  CH(I-1,K,4) = WA3(I-2)*CR4 - WA3(I-1)*CI4
                  CH(I,K,4) = WA3(I-2)*CI4 + WA3(I-1)*CR4
               END DO
            END DO
            IF (MOD(IDO,2) == 1) RETURN 
         ENDIF
         DO K = 1, L1
            TI1 = CC(1,2,K) + CC(1,4,K)
            TI2 = CC(1,4,K) - CC(1,2,K)
            TR1 = CC(IDO,1,K) - CC(IDO,3,K)
            TR2 = CC(IDO,1,K) + CC(IDO,3,K)
            CH(IDO,K,1) = TR2 + TR2
            CH(IDO,K,2) = SQRT2*(TR1 - TI1)
            CH(IDO,K,3) = TI2 + TI2
            CH(IDO,K,4) = -SQRT2*(TR1 + TI1)
         END DO
      ENDIF
      RETURN 
      END SUBROUTINE RADB4


      SUBROUTINE RADB5(IDO, L1, CC, CH, WA1, WA2, WA3, WA4)
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER , INTENT(IN) :: IDO
      INTEGER , INTENT(IN) :: L1
      REAL(fish_kind) , INTENT(IN) :: CC(IDO,5,L1)
      REAL(fish_kind) , INTENT(OUT) :: CH(IDO,L1,5)
      REAL(fish_kind) , INTENT(IN) :: WA1(*)
      REAL(fish_kind) , INTENT(IN) :: WA2(*)
      REAL(fish_kind) , INTENT(IN) :: WA3(*)
      REAL(fish_kind) , INTENT(IN) :: WA4(*)
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: K, IDP2, I, IC
      REAL :: TR11, TI11, TR12, TI12, TI5, TI4, TR2, TR3, CR2, CR3, CI5
      REAL :: CI4, TI2, TI3, TR5, TR4, CI2, CI3, CR5, CR4, DR3, DR4, DI3
      REAL :: DI4, DR5, DR2, DI5, DI2
!-----------------------------------------------
      DATA TR11, TI11, TR12, TI12/ 0.309016994374947, 0.951056516295154, -.809016994374947, 0.587785252292473/ 
      DO K = 1, L1
         TI5 = CC(1,3,K) + CC(1,3,K)
         TI4 = CC(1,5,K) + CC(1,5,K)
         TR2 = CC(IDO,2,K) + CC(IDO,2,K)
         TR3 = CC(IDO,4,K) + CC(IDO,4,K)
         CH(1,K,1) = CC(1,1,K) + TR2 + TR3
         CR2 = CC(1,1,K) + TR11*TR2 + TR12*TR3
         CR3 = CC(1,1,K) + TR12*TR2 + TR11*TR3
         CI5 = TI11*TI5 + TI12*TI4
         CI4 = TI12*TI5 - TI11*TI4
         CH(1,K,2) = CR2 - CI5
         CH(1,K,3) = CR3 - CI4
         CH(1,K,4) = CR3 + CI4
         CH(1,K,5) = CR2 + CI5
      END DO
      IF (IDO == 1) RETURN 
      IDP2 = IDO + 2
      DO K = 1, L1
         DO I = 3, IDO, 2
            IC = IDP2 - I
            TI5 = CC(I,3,K) + CC(IC,2,K)
            TI2 = CC(I,3,K) - CC(IC,2,K)
            TI4 = CC(I,5,K) + CC(IC,4,K)
            TI3 = CC(I,5,K) - CC(IC,4,K)
            TR5 = CC(I-1,3,K) - CC(IC-1,2,K)
            TR2 = CC(I-1,3,K) + CC(IC-1,2,K)
            TR4 = CC(I-1,5,K) - CC(IC-1,4,K)
            TR3 = CC(I-1,5,K) + CC(IC-1,4,K)
            CH(I-1,K,1) = CC(I-1,1,K) + TR2 + TR3
            CH(I,K,1) = CC(I,1,K) + TI2 + TI3
            CR2 = CC(I-1,1,K) + TR11*TR2 + TR12*TR3
            CI2 = CC(I,1,K) + TR11*TI2 + TR12*TI3
            CR3 = CC(I-1,1,K) + TR12*TR2 + TR11*TR3
            CI3 = CC(I,1,K) + TR12*TI2 + TR11*TI3
            CR5 = TI11*TR5 + TI12*TR4
            CI5 = TI11*TI5 + TI12*TI4
            CR4 = TI12*TR5 - TI11*TR4
            CI4 = TI12*TI5 - TI11*TI4
            DR3 = CR3 - CI4
            DR4 = CR3 + CI4
            DI3 = CI3 + CR4
            DI4 = CI3 - CR4
            DR5 = CR2 + CI5
            DR2 = CR2 - CI5
            DI5 = CI2 - CR5
            DI2 = CI2 + CR5
            CH(I-1,K,2) = WA1(I-2)*DR2 - WA1(I-1)*DI2
            CH(I,K,2) = WA1(I-2)*DI2 + WA1(I-1)*DR2
            CH(I-1,K,3) = WA2(I-2)*DR3 - WA2(I-1)*DI3
            CH(I,K,3) = WA2(I-2)*DI3 + WA2(I-1)*DR3
            CH(I-1,K,4) = WA3(I-2)*DR4 - WA3(I-1)*DI4
            CH(I,K,4) = WA3(I-2)*DI4 + WA3(I-1)*DR4
            CH(I-1,K,5) = WA4(I-2)*DR5 - WA4(I-1)*DI5
            CH(I,K,5) = WA4(I-2)*DI5 + WA4(I-1)*DR5
         END DO
      END DO
      RETURN 
      END SUBROUTINE RADB5


      SUBROUTINE RADBG(IDO, IP, L1, IDL1, CC, C1, C2, CH, CH2, WA)
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER , INTENT(IN) :: IDO
      INTEGER , INTENT(IN) :: IP
      INTEGER , INTENT(IN) :: L1
      INTEGER , INTENT(IN) :: IDL1
      REAL(fish_kind) , INTENT(IN) :: CC(IDO,IP,L1)
      REAL(fish_kind) , INTENT(INOUT) :: C1(IDO,L1,IP)
      REAL(fish_kind) , INTENT(INOUT) :: C2(IDL1,IP)
      REAL(fish_kind) , INTENT(INOUT) :: CH(IDO,L1,IP)
      REAL(fish_kind) , INTENT(INOUT) :: CH2(IDL1,IP)
      REAL(fish_kind) , INTENT(IN) :: WA(*)
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER::IDP2,NBD,IPP2,IPPH,K,I,J,JC,J2,IC,L,LC,IK,IS,IDIJ
      REAL::TPI,DUM,ARG,DCP,DSP,AR1,AI1,AR1H,DC2,DS2,AR2,AI2,AR2H
!-----------------------------------------------
      TPI = 8.0*ATAN(1.0)
      ARG = TPI/FLOAT(IP)
      DCP = COS(ARG)
      DSP = SIN(ARG)
      IDP2 = IDO + 2
      NBD = (IDO - 1)/2
      IPP2 = IP + 2
      IPPH = (IP + 1)/2
      IF (IDO >= L1) THEN
         CH(:,:,1) = CC(:,1,:)
      ELSE
         CH(:,:,1) = CC(:,1,:)
      ENDIF
      DO J = 2, IPPH
         JC = IPP2 - J
         J2 = J + J
         CH(1,:,J) = CC(IDO,J2-2,:) + CC(IDO,J2-2,:)
         CH(1,:,JC) = CC(1,J2-1,:) + CC(1,J2-1,:)
      END DO
      IF (IDO /= 1) THEN
         IF (NBD >= L1) THEN
            DO J = 2, IPPH
               JC = IPP2 - J
               CH(2:IDO-1:2,:,J) = CC(2:IDO-1:2,2*J-1,:) + CC(IDP2-4:IDP2-1-IDO:(-2),2*J-2,:)
               CH(2:IDO-1:2,:,JC) = CC(2:IDO-1:2,2*J-1,:) - CC(IDP2-4:IDP2-1-IDO:(-2),2*J-2,:)
               CH(3:IDO:2,:,J) = CC(3:IDO:2,2*J-1,:) - CC(IDP2-3:IDP2-IDO:(-2),2*J-2,:)
               CH(3:IDO:2,:,JC) = CC(3:IDO:2,2*J-1,:) + CC(IDP2-3:IDP2-IDO:(-2),2*J-2,:)
            END DO
         ELSE
            DO J = 2, IPPH
               JC = IPP2 - J
               CH(2:IDO-1:2,:,J) = CC(2:IDO-1:2,2*J-1,:) + CC(IDP2-4:IDP2-1-IDO:(-2),2*J-2,:)
               CH(2:IDO-1:2,:,JC) = CC(2:IDO-1:2,2*J-1,:) - CC(IDP2-4:IDP2-1-IDO:(-2),2*J-2,:)
               CH(3:IDO:2,:,J) = CC(3:IDO:2,2*J-1,:) - CC(IDP2-3:IDP2-IDO:(-2),2*J-2,:)
               CH(3:IDO:2,:,JC) = CC(3:IDO:2,2*J-1,:) + CC(IDP2-3:IDP2-IDO:(-2),2*J-2,:)
            END DO
         ENDIF
      ENDIF
      AR1 = 1.
      AI1 = 0.
      DO L = 2, IPPH
         LC = IPP2 - L
         AR1H = DCP*AR1 - DSP*AI1
         AI1 = DCP*AI1 + DSP*AR1
         AR1 = AR1H
         C2(:,L) = CH2(:,1) + AR1*CH2(:,2)
         C2(:,LC) = AI1*CH2(:,IP)
         DC2 = AR1
         DS2 = AI1
         AR2 = AR1
         AI2 = AI1
         DO J = 3, IPPH
            JC = IPP2 - J
            AR2H = DC2*AR2 - DS2*AI2
            AI2 = DC2*AI2 + DS2*AR2
            AR2 = AR2H
            C2(:,L) = C2(:,L) + AR2*CH2(:,J)
            C2(:,LC) = C2(:,LC) + AI2*CH2(:,JC)
         END DO
      END DO
      DO J = 2, IPPH
         CH2(:,1) = CH2(:,1) + CH2(:,J)
      END DO
      DO J = 2, IPPH
         JC = IPP2 - J
         CH(1,:,J) = C1(1,:,J) - C1(1,:,JC)
         CH(1,:,JC) = C1(1,:,J) + C1(1,:,JC)
      END DO
      IF (IDO /= 1) THEN
         IF (NBD >= L1) THEN
            DO J = 2, IPPH
               JC = IPP2 - J
               CH(2:IDO-1:2,:,J) = C1(2:IDO-1:2,:,J) - C1(3:IDO:2,:,JC)
               CH(2:IDO-1:2,:,JC) = C1(2:IDO-1:2,:,J) + C1(3:IDO:2,:,JC)
               CH(3:IDO:2,:,J) = C1(3:IDO:2,:,J) + C1(2:IDO-1:2,:,JC)
               CH(3:IDO:2,:,JC) = C1(3:IDO:2,:,J) - C1(2:IDO-1:2,:,JC)
            END DO
         ELSE
            DO J = 2, IPPH
               JC = IPP2 - J
               CH(2:IDO-1:2,:,J) = C1(2:IDO-1:2,:,J) - C1(3:IDO:2,:,JC)
               CH(2:IDO-1:2,:,JC) = C1(2:IDO-1:2,:,J) + C1(3:IDO:2,:,JC)
               CH(3:IDO:2,:,J) = C1(3:IDO:2,:,J) + C1(2:IDO-1:2,:,JC)
               CH(3:IDO:2,:,JC) = C1(3:IDO:2,:,J) - C1(2:IDO-1:2,:,JC)
            END DO
         ENDIF
      ENDIF
      IF (IDO == 1) RETURN 
      C2(:,1) = CH2(:,1)
      C1(1,:,2:IP) = CH(1,:,2:IP)
      IF (NBD <= L1) THEN
         IS = -IDO
         DO J = 2, IP
            IS = IS + IDO
            IDIJ = IS
            DO I = 3, IDO, 2
               IDIJ = IDIJ + 2
               C1(I-1,:,J) = WA(IDIJ-1)*CH(I-1,:,J) - WA(IDIJ)*CH(I,:,J)
               C1(I,:,J) = WA(IDIJ-1)*CH(I,:,J) + WA(IDIJ)*CH(I-1,:,J)
            END DO
         END DO
      ELSE
         IS = -IDO
         DO J = 2, IP
            IS = IS + IDO
            DO K = 1, L1
               IDIJ = IS
               C1(2:IDO-1:2,K,J) = WA(IDIJ+1:IDO-2+IDIJ:2)*CH(2:IDO-1:2, K,J) - WA(IDIJ+2:IDO-1+IDIJ:2)*CH(3:IDO:2,K,J)
               C1(3:IDO:2,K,J) = WA(IDIJ+1:IDO-2+IDIJ:2)*CH(3:IDO:2,K,J)      + WA(IDIJ+2:IDO-1+IDIJ:2)*CH(2:IDO-1:2,K,J)
            END DO
         END DO
      ENDIF
      RETURN 
      END SUBROUTINE RADBG


      SUBROUTINE RFFTF(N, R, WSAVE)
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER  :: N
      REAL(fish_kind)  :: R(*)
      REAL(fish_kind)  :: WSAVE(*)
!-----------------------------------------------
!
      IF (N == 1) RETURN 
      CALL RFFTF1 (N, R, WSAVE, WSAVE(N+1), WSAVE(2*N+1))
      RETURN 
      END SUBROUTINE RFFTF


      SUBROUTINE RFFTF1(N, C, CH, WA, IFAC)
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER , INTENT(IN) :: N
      REAL(fish_kind) , INTENT(IN) :: IFAC(*)
      REAL(fish_kind)  :: C(*)
      REAL(fish_kind)  :: CH(*)
      REAL(fish_kind)  :: WA(*)
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER::NF,NA,L2,IW,K1,KH,IP,L1,IDO,IDL1,IX2,IX3,IX4,I
!-----------------------------------------------
      NF = IFAC(2)
      NA = 1
      L2 = N
      IW = N
      DO K1 = 1, NF
         KH = NF - K1
         IP = IFAC(KH+3)
         L1 = L2/IP
         IDO = N/L2
         IDL1 = IDO*L1
         IW = IW - (IP - 1)*IDO
         NA = 1 - NA
         IF (IP == 4) THEN
            IX2 = IW + IDO
            IX3 = IX2 + IDO
            IF (NA == 0) THEN
               CALL RADF4 (IDO, L1, C, CH, WA(IW), WA(IX2), WA(IX3))
               GO TO 110
            ENDIF
            CALL RADF4 (IDO, L1, CH, C, WA(IW), WA(IX2), WA(IX3))
            GO TO 110
         ENDIF
         IF (IP == 2) THEN
            IF (NA == 0) THEN
               CALL RADF2 (IDO, L1, C, CH, WA(IW))
               GO TO 110
            ENDIF
            CALL RADF2 (IDO, L1, CH, C, WA(IW))
            GO TO 110
         ENDIF
  104    CONTINUE
         IF (IP == 3) THEN
            IX2 = IW + IDO
            IF (NA == 0) THEN
               CALL RADF3 (IDO, L1, C, CH, WA(IW), WA(IX2))
               GO TO 110
            ENDIF
            CALL RADF3 (IDO, L1, CH, C, WA(IW), WA(IX2))
            GO TO 110
         ENDIF
  106    CONTINUE
         IF (IP == 5) THEN
            IX2 = IW + IDO
            IX3 = IX2 + IDO
            IX4 = IX3 + IDO
            IF (NA == 0) THEN
               CALL RADF5(IDO,L1,C,CH,WA(IW),WA(IX2),WA(IX3),WA(IX4))
               GO TO 110
            ENDIF
            CALL RADF5(IDO,L1,CH,C,WA(IW),WA(IX2),WA(IX3),WA(IX4))
            GO TO 110
         ENDIF
  108    CONTINUE
         IF (IDO == 1) NA = 1 - NA
         IF (NA == 0) THEN
            CALL RADFG (IDO, IP, L1, IDL1, C, C, C, CH, CH, WA(IW))
            NA = 1
         ELSE
            CALL RADFG (IDO, IP, L1, IDL1, CH, CH, CH, C, C, WA(IW))
            NA = 0
         ENDIF
  110    CONTINUE
         L2 = L1
      END DO
      IF (NA == 1) RETURN 
      C(:N) = CH(:N)
      RETURN 
      END SUBROUTINE RFFTF1


      SUBROUTINE RADF2(IDO, L1, CC, CH, WA1)
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER , INTENT(IN) :: IDO
      INTEGER , INTENT(IN) :: L1
      REAL(fish_kind) , INTENT(IN) :: CC(IDO,L1,2)
      REAL(fish_kind) , INTENT(OUT) :: CH(IDO,2,L1)
      REAL(fish_kind) , INTENT(IN) :: WA1(*)
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: K, IDP2, I, IC
      REAL :: TR2, TI2
!-----------------------------------------------
      CH(1,1,:) = CC(1,:,1) + CC(1,:,2)
      CH(IDO,2,:) = CC(1,:,1) - CC(1,:,2)
      IF (IDO - 2 >= 0) THEN
         IF (IDO - 2 /= 0) THEN
            IDP2 = IDO + 2
            DO K = 1, L1
               DO I = 3, IDO, 2
                  IC = IDP2 - I
                  TR2 = WA1(I-2)*CC(I-1,K,2) + WA1(I-1)*CC(I,K,2)
                  TI2 = WA1(I-2)*CC(I,K,2) - WA1(I-1)*CC(I-1,K,2)
                  CH(I,1,K) = CC(I,K,1) + TI2
                  CH(IC,2,K) = TI2 - CC(I,K,1)
                  CH(I-1,1,K) = CC(I-1,K,1) + TR2
                  CH(IC-1,2,K) = CC(I-1,K,1) - TR2
               END DO
            END DO
            IF (MOD(IDO,2) == 1) RETURN 
         ENDIF
         CH(1,2,:) = -CC(IDO,:,2)
         CH(IDO,1,:) = CC(IDO,:,1)
      ENDIF
      RETURN 
      END SUBROUTINE RADF2


      SUBROUTINE RADF3(IDO, L1, CC, CH, WA1, WA2)
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER , INTENT(IN) :: IDO
      INTEGER , INTENT(IN) :: L1
      REAL(fish_kind) , INTENT(IN) :: CC(IDO,L1,3)
      REAL(fish_kind) , INTENT(OUT) :: CH(IDO,3,L1)
      REAL(fish_kind) , INTENT(IN) :: WA1(*)
      REAL(fish_kind) , INTENT(IN) :: WA2(*)
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: K, IDP2, I, IC
      REAL::TAUR,TAUI,CR2,DR2,DI2,DR3,DI3,CI2,TR2,TI2,TR3,TI3
!-----------------------------------------------
      DATA TAUR, TAUI/ -.5, 0.866025403784439/ 
      DO K = 1, L1
         CR2 = CC(1,K,2) + CC(1,K,3)
         CH(1,1,K) = CC(1,K,1) + CR2
         CH(1,3,K) = TAUI*(CC(1,K,3)-CC(1,K,2))
         CH(IDO,2,K) = CC(1,K,1) + TAUR*CR2
      END DO
      IF (IDO == 1) RETURN 
      IDP2 = IDO + 2
      DO K = 1, L1
         DO I = 3, IDO, 2
            IC = IDP2 - I
            DR2 = WA1(I-2)*CC(I-1,K,2) + WA1(I-1)*CC(I,K,2)
            DI2 = WA1(I-2)*CC(I,K,2) - WA1(I-1)*CC(I-1,K,2)
            DR3 = WA2(I-2)*CC(I-1,K,3) + WA2(I-1)*CC(I,K,3)
            DI3 = WA2(I-2)*CC(I,K,3) - WA2(I-1)*CC(I-1,K,3)
            CR2 = DR2 + DR3
            CI2 = DI2 + DI3
            CH(I-1,1,K) = CC(I-1,K,1) + CR2
            CH(I,1,K) = CC(I,K,1) + CI2
            TR2 = CC(I-1,K,1) + TAUR*CR2
            TI2 = CC(I,K,1) + TAUR*CI2
            TR3 = TAUI*(DI2 - DI3)
            TI3 = TAUI*(DR3 - DR2)
            CH(I-1,3,K) = TR2 + TR3
            CH(IC-1,2,K) = TR2 - TR3
            CH(I,3,K) = TI2 + TI3
            CH(IC,2,K) = TI3 - TI2
         END DO
      END DO
      RETURN 
      END SUBROUTINE RADF3


      SUBROUTINE RADF4(IDO, L1, CC, CH, WA1, WA2, WA3)
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER , INTENT(IN) :: IDO
      INTEGER , INTENT(IN) :: L1
      REAL(fish_kind) , INTENT(IN) :: CC(IDO,L1,4)
      REAL(fish_kind) , INTENT(OUT) :: CH(IDO,4,L1)
      REAL(fish_kind) , INTENT(IN) :: WA1(*)
      REAL(fish_kind) , INTENT(IN) :: WA2(*)
      REAL(fish_kind) , INTENT(IN) :: WA3(*)
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: K, IDP2, I, IC
      REAL :: HSQT2, TR1, TR2, CR2, CI2, CR3, CI3, CR4, CI4, TR4, TI1 
      REAL ::  TI4, TI2, TI3, TR3
!-----------------------------------------------
      DATA HSQT2/ 0.7071067811865475/ 
      DO K = 1, L1
         TR1 = CC(1,K,2) + CC(1,K,4)
         TR2 = CC(1,K,1) + CC(1,K,3)
         CH(1,1,K) = TR1 + TR2
         CH(IDO,4,K) = TR2 - TR1
         CH(IDO,2,K) = CC(1,K,1) - CC(1,K,3)
         CH(1,3,K) = CC(1,K,4) - CC(1,K,2)
      END DO
      IF (IDO - 2 >= 0) THEN
         IF (IDO - 2 /= 0) THEN
            IDP2 = IDO + 2
            DO K = 1, L1
               DO I = 3, IDO, 2
                  IC = IDP2 - I
                  CR2 = WA1(I-2)*CC(I-1,K,2) + WA1(I-1)*CC(I,K,2)
                  CI2 = WA1(I-2)*CC(I,K,2) - WA1(I-1)*CC(I-1,K,2)
                  CR3 = WA2(I-2)*CC(I-1,K,3) + WA2(I-1)*CC(I,K,3)
                  CI3 = WA2(I-2)*CC(I,K,3) - WA2(I-1)*CC(I-1,K,3)
                  CR4 = WA3(I-2)*CC(I-1,K,4) + WA3(I-1)*CC(I,K,4)
                  CI4 = WA3(I-2)*CC(I,K,4) - WA3(I-1)*CC(I-1,K,4)
                  TR1 = CR2 + CR4
                  TR4 = CR4 - CR2
                  TI1 = CI2 + CI4
                  TI4 = CI2 - CI4
                  TI2 = CC(I,K,1) + CI3
                  TI3 = CC(I,K,1) - CI3
                  TR2 = CC(I-1,K,1) + CR3
                  TR3 = CC(I-1,K,1) - CR3
                  CH(I-1,1,K) = TR1 + TR2
                  CH(IC-1,4,K) = TR2 - TR1
                  CH(I,1,K) = TI1 + TI2
                  CH(IC,4,K) = TI1 - TI2
                  CH(I-1,3,K) = TI4 + TR3
                  CH(IC-1,2,K) = TR3 - TI4
                  CH(I,3,K) = TR4 + TI3
                  CH(IC,2,K) = TR4 - TI3
               END DO
            END DO
            IF (MOD(IDO,2) == 1) RETURN 
         ENDIF
         DO K = 1, L1
            TI1 = -HSQT2*(CC(IDO,K,2)+CC(IDO,K,4))
            TR1 = HSQT2*(CC(IDO,K,2)-CC(IDO,K,4))
            CH(IDO,1,K) = TR1 + CC(IDO,K,1)
            CH(IDO,3,K) = CC(IDO,K,1) - TR1
            CH(1,2,K) = TI1 - CC(IDO,K,3)
            CH(1,4,K) = TI1 + CC(IDO,K,3)
         END DO
      ENDIF
      RETURN 
      END SUBROUTINE RADF4


      SUBROUTINE RADF5(IDO, L1, CC, CH, WA1, WA2, WA3, WA4)
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER , INTENT(IN) :: IDO
      INTEGER , INTENT(IN) :: L1
      REAL(fish_kind) , INTENT(IN) :: CC(IDO,L1,5)
      REAL(fish_kind) , INTENT(OUT) :: CH(IDO,5,L1)
      REAL(fish_kind) , INTENT(IN) :: WA1(*)
      REAL(fish_kind) , INTENT(IN) :: WA2(*)
      REAL(fish_kind) , INTENT(IN) :: WA3(*)
      REAL(fish_kind) , INTENT(IN) :: WA4(*)
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: K, IDP2, I, IC
      REAL :: TR11, TI11, TR12, TI12, CR2, CI5, CR3, CI4, DR2, DI2, DR3
      REAL :: DI3, DR4, DI4, DR5, DI5, CR5, CI2, CR4, CI3, TR2, TI2, TR3
      REAL :: TI3, TR5, TI5, TR4, TI4
!-----------------------------------------------
      DATA TR11, TI11, TR12, TI12/ 0.309016994374947, 0.951056516295154, -.809016994374947, 0.587785252292473/ 
      DO K = 1, L1
         CR2 = CC(1,K,5) + CC(1,K,2)
         CI5 = CC(1,K,5) - CC(1,K,2)
         CR3 = CC(1,K,4) + CC(1,K,3)
         CI4 = CC(1,K,4) - CC(1,K,3)
         CH(1,1,K) = CC(1,K,1) + CR2 + CR3
         CH(IDO,2,K) = CC(1,K,1) + TR11*CR2 + TR12*CR3
         CH(1,3,K) = TI11*CI5 + TI12*CI4
         CH(IDO,4,K) = CC(1,K,1) + TR12*CR2 + TR11*CR3
         CH(1,5,K) = TI12*CI5 - TI11*CI4
      END DO
      IF (IDO == 1) RETURN 
      IDP2 = IDO + 2
      DO K = 1, L1
         DO I = 3, IDO, 2
            IC = IDP2 - I
            DR2 = WA1(I-2)*CC(I-1,K,2) + WA1(I-1)*CC(I,K,2)
            DI2 = WA1(I-2)*CC(I,K,2) - WA1(I-1)*CC(I-1,K,2)
            DR3 = WA2(I-2)*CC(I-1,K,3) + WA2(I-1)*CC(I,K,3)
            DI3 = WA2(I-2)*CC(I,K,3) - WA2(I-1)*CC(I-1,K,3)
            DR4 = WA3(I-2)*CC(I-1,K,4) + WA3(I-1)*CC(I,K,4)
            DI4 = WA3(I-2)*CC(I,K,4) - WA3(I-1)*CC(I-1,K,4)
            DR5 = WA4(I-2)*CC(I-1,K,5) + WA4(I-1)*CC(I,K,5)
            DI5 = WA4(I-2)*CC(I,K,5) - WA4(I-1)*CC(I-1,K,5)
            CR2 = DR2 + DR5
            CI5 = DR5 - DR2
            CR5 = DI2 - DI5
            CI2 = DI2 + DI5
            CR3 = DR3 + DR4
            CI4 = DR4 - DR3
            CR4 = DI3 - DI4
            CI3 = DI3 + DI4
            CH(I-1,1,K) = CC(I-1,K,1) + CR2 + CR3
            CH(I,1,K) = CC(I,K,1) + CI2 + CI3
            TR2 = CC(I-1,K,1) + TR11*CR2 + TR12*CR3
            TI2 = CC(I,K,1) + TR11*CI2 + TR12*CI3
            TR3 = CC(I-1,K,1) + TR12*CR2 + TR11*CR3
            TI3 = CC(I,K,1) + TR12*CI2 + TR11*CI3
            TR5 = TI11*CR5 + TI12*CR4
            TI5 = TI11*CI5 + TI12*CI4
            TR4 = TI12*CR5 - TI11*CR4
            TI4 = TI12*CI5 - TI11*CI4
            CH(I-1,3,K) = TR2 + TR5
            CH(IC-1,2,K) = TR2 - TR5
            CH(I,3,K) = TI2 + TI5
            CH(IC,2,K) = TI5 - TI2
            CH(I-1,5,K) = TR3 + TR4
            CH(IC-1,4,K) = TR3 - TR4
            CH(I,5,K) = TI3 + TI4
            CH(IC,4,K) = TI4 - TI3
         END DO
      END DO
      RETURN 
      END SUBROUTINE RADF5


      SUBROUTINE RADFG(IDO, IP, L1, IDL1, CC, C1, C2, CH, CH2, WA)
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER , INTENT(IN) :: IDO
      INTEGER , INTENT(IN) :: IP
      INTEGER , INTENT(IN) :: L1
      INTEGER , INTENT(IN) :: IDL1
      REAL(fish_kind) , INTENT(OUT) :: CC(IDO,IP,L1)
      REAL(fish_kind) , INTENT(INOUT) :: C1(IDO,L1,IP)
      REAL(fish_kind) , INTENT(INOUT) :: C2(IDL1,IP)
      REAL(fish_kind) , INTENT(INOUT) :: CH(IDO,L1,IP)
      REAL(fish_kind) , INTENT(INOUT) :: CH2(IDL1,IP)
      REAL(fish_kind) , INTENT(IN) :: WA(*)
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER::IPPH,IPP2,IDP2,NBD,IK,J,K,IS,IDIJ,I,JC,L,LC,J2,IC
      REAL::TPI,DUM,ARG,DCP,DSP,AR1,AI1,AR1H,DC2,DS2,AR2,AI2,AR2H
!-----------------------------------------------
      TPI = 8.0*ATAN(1.0)
      ARG = TPI/FLOAT(IP)
      DCP = COS(ARG)
      DSP = SIN(ARG)
      IPPH = (IP + 1)/2
      IPP2 = IP + 2
      IDP2 = IDO + 2
      NBD = (IDO - 1)/2
      IF (IDO /= 1) THEN
         CH2(:,1) = C2(:,1)
         CH(1,:,2:IP) = C1(1,:,2:IP)
         IF (NBD <= L1) THEN
            IS = -IDO
            DO J = 2, IP
               IS = IS + IDO
               IDIJ = IS
               DO I = 3, IDO, 2
                  IDIJ = IDIJ + 2
                  CH(I-1,:,J)=WA(IDIJ-1)*C1(I-1,:,J)+WA(IDIJ)*C1(I,:,J)
                  CH(I,:,J)=WA(IDIJ-1)*C1(I,:,J)-WA(IDIJ)*C1(I-1,:,J)
               END DO
            END DO
         ELSE
            IS = -IDO
            DO J = 2, IP
               IS = IS + IDO
               DO K = 1, L1
                  IDIJ = IS
                  CH(2:IDO-1:2,K,J) = WA(IDIJ+1:IDO-2+IDIJ:2)*C1(2:IDO-1:2,K,J) + WA(IDIJ+2:IDO-1+IDIJ:2)*C1(3:IDO:2,K,J)
                  CH(3:IDO:2,K,J) = WA(IDIJ+1:IDO-2+IDIJ:2)*C1(3:IDO:2,K,J) - WA(IDIJ+2:IDO-1+IDIJ:2)*C1(2:IDO-1:2,K,J)
               END DO
            END DO
         ENDIF
         IF (NBD >= L1) THEN
            DO J = 2, IPPH
               JC = IPP2 - J
               C1(2:IDO-1:2,:,J)=CH(2:IDO-1:2,:,J)+CH(2:IDO-1:2,:,JC)
               C1(2:IDO-1:2,:,JC) = CH(3:IDO:2,:,J) - CH(3:IDO:2,:,JC)
               C1(3:IDO:2,:,J) = CH(3:IDO:2,:,J) + CH(3:IDO:2,:,JC)
               C1(3:IDO:2,:,JC) = CH(2:IDO-1:2,:,JC) - CH(2:IDO-1:2,:,J)
            END DO
            GO TO 121
         ENDIF
         DO J = 2, IPPH
            JC = IPP2 - J
            C1(2:IDO-1:2,:,J) = CH(2:IDO-1:2,:,J) + CH(2:IDO-1:2,:,JC)
            C1(2:IDO-1:2,:,JC) = CH(3:IDO:2,:,J) - CH(3:IDO:2,:,JC)
            C1(3:IDO:2,:,J) = CH(3:IDO:2,:,J) + CH(3:IDO:2,:,JC)
            C1(3:IDO:2,:,JC) = CH(2:IDO-1:2,:,JC) - CH(2:IDO-1:2,:,J)
         END DO
         GO TO 121
      ENDIF
      C2(:,1) = CH2(:,1)
  121 CONTINUE
      DO J = 2, IPPH
         JC = IPP2 - J
         C1(1,:,J) = CH(1,:,J) + CH(1,:,JC)
         C1(1,:,JC) = CH(1,:,JC) - CH(1,:,J)
      END DO
!
      AR1 = 1.
      AI1 = 0.
      DO L = 2, IPPH
         LC = IPP2 - L
         AR1H = DCP*AR1 - DSP*AI1
         AI1 = DCP*AI1 + DSP*AR1
         AR1 = AR1H
         CH2(:,L) = C2(:,1) + AR1*C2(:,2)
         CH2(:,LC) = AI1*C2(:,IP)
         DC2 = AR1
         DS2 = AI1
         AR2 = AR1
         AI2 = AI1
         DO J = 3, IPPH
            JC = IPP2 - J
            AR2H = DC2*AR2 - DS2*AI2
            AI2 = DC2*AI2 + DS2*AR2
            AR2 = AR2H
            CH2(:,L) = CH2(:,L) + AR2*C2(:,J)
            CH2(:,LC) = CH2(:,LC) + AI2*C2(:,JC)
         END DO
      END DO
      DO J = 2, IPPH
         CH2(:,1) = CH2(:,1) + C2(:,J)
      END DO
!
      IF (IDO >= L1) THEN
         CC(:,1,:) = CH(:,:,1)
      ELSE
         CC(:,1,:) = CH(:,:,1)
      ENDIF
      CC(IDO,2:(IPPH-1)*2:2,:) = TRANSPOSE(CH(1,:,2:IPPH))
      CC(1,3:IPPH*2-1:2,:) = TRANSPOSE(CH(1,:,IPP2-2:IPP2-IPPH:(-1)))
      IF (IDO == 1) RETURN 
      IF (NBD >= L1) THEN
         CC(2:IDO-1:2,3:IPPH*2-1:2,:) = RESHAPE(SOURCE = CH(2:IDO-1:2,:, 2:IPPH)+CH(2:IDO-1:2,:,IPP2-2:IPP2-IPPH:(-1)),SHAPE = (/(IDO -1)/2,IPPH-1,L1/),ORDER = (/1,3,2/))
         CC(IDP2-4:IDP2-1-IDO:(-2),2:(IPPH-1)*2:2,:) = RESHAPE(SOURCE = CH(2:IDO-1:2,:,2:IPPH)-CH(2:IDO-1:2,:,IPP2-2:IPP2-IPPH:(-1)) ,SHAPE = (/(IDO-1)/2,IPPH-1,L1/),ORDER = (/1,3,2/))
         CC(3:IDO:2,3:IPPH*2-1:2,:) = RESHAPE(SOURCE = CH(3:IDO:2,:,2: IPPH)+CH(3:IDO:2,:,IPP2-2:IPP2-IPPH:(-1)),SHAPE = (/(IDO-1)/ 2,IPPH-1,L1/),ORDER = (/1,3,2/))
         CC(IDP2-3:IDP2-IDO:(-2),2:(IPPH-1)*2:2,:) = RESHAPE(SOURCE = CH (3:IDO:2,:,IPP2-2:IPP2-IPPH:(-1))-CH(3:IDO:2,:,2:IPPH),SHAPE = (/(IDO-1)/2,IPPH-1,L1/),ORDER = (/1,3,2/))
         RETURN 
      ENDIF
      CC(2:IDO-1:2,3:IPPH*2-1:2,:) = RESHAPE(SOURCE = CH(2:IDO-1:2,:,2: IPPH)+CH(2:IDO-1:2,:,IPP2-2:IPP2-IPPH:(-1)),SHAPE = (/(IDO-1)/2 ,IPPH-1,L1/),ORDER = (/1,3,2/))
      CC(IDP2-4:IDP2-1-IDO:(-2),2:(IPPH-1)*2:2,:) = RESHAPE(SOURCE = CH( 2:IDO-1:2,:,2:IPPH)-CH(2:IDO-1:2,:,IPP2-2:IPP2-IPPH:(-1)),SHAPE = (/(IDO-1)/2,IPPH-1,L1/),ORDER = (/1,3,2/))
      CC(3:IDO:2,3:IPPH*2-1:2,:) = RESHAPE(SOURCE = CH(3:IDO:2,:,2:IPPH) +CH(3:IDO:2,:,IPP2-2:IPP2-IPPH:(-1)),SHAPE = (/(IDO-1)/2,IPPH-1 ,L1/),ORDER = (/1,3,2/))
      CC(IDP2-3:IDP2-IDO:(-2),2:(IPPH-1)*2:2,:) = RESHAPE(SOURCE = CH(3: IDO:2,:,IPP2-2:IPP2-IPPH:(-1))-CH(3:IDO:2,:,2:IPPH),SHAPE = (/( IDO-1)/2,IPPH-1,L1/),ORDER = (/1,3,2/))
      RETURN 
      END SUBROUTINE RADFG
!     this function is define in the file comf.f
! SEPTEMBER 1973    VERSION 1
! APRIL     1976    VERSION 2
! JANUARY   1978    VERSION 3
! DECEMBER  1979    VERSION 3.1
! FEBRUARY  1985    DOCUMENTATION UPGRADE
! NOVEMBER  1988    VERSION 3.2, FORTRAN 77 CHANGES
! June 2004 2004    fortran 90 updates
!-----------------------------------------------------------------------
!     END



!
!     file tblktri.f
!
!     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
!     *                                                               *
!     *                  copyright (c) 2005 by UCAR                   *
!     *                                                               *
!     *       University Corporation for Atmospheric Research         *
!     *                                                               *
!     *                      all rights reserved                      *
!     *                                                               *
!     *                    FISHPACK90  version 1.1                    *
!     *                                                               *
!     *                 A Package of Fortran 77 and 90                *
!     *                                                               *
!     *                Subroutines and Example Programs               *
!     *                                                               *
!     *               for Modeling Geophysical Processes              *
!     *                                                               *
!     *                             by                                *
!     *                                                               *
!     *        John Adams, Paul Swarztrauber and Roland Sweet         *
!     *                                                               *
!     *                             of                                *
!     *                                                               *
!     *         the National Center for Atmospheric Research          *
!     *                                                               *
!     *                Boulder, Colorado  (80307)  U.S.A.             *
!     *                                                               *
!     *                   which is sponsored by                       *
!     *                                                               *
!     *              the National Science Foundation                  *
!     *                                                               *
!     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
!
!     PROGRAM TO ILLUSTRATE THE USE OF SUBROUTINE BLKTRI TO
!     SOLVE THE EQUATION
!
!     .5/S*(D/DS)(.5/S*DU/DS)+.5/T*(D/DT)(.5/T*DU/DT)
!                                                          (1)
!                   = 15/4*S*T*(S**4+T**4)
!
!     ON THE RECTANGLE 0 .LT. S .LT. 1 AND 0 .LT. T .LT. 1
!     WITH THE BOUNDARY CONDITIONS
!
!     U(0,T) = 0
!                            0 .LE. T .LE. 1
!     U(1,T) = T**5
!
!     AND
!
!     U(S,0) = 0
!                            0 .LE. S .LE. 1
!     U(S,1) = S**5
!
!     THE EXACT SOLUTION OF THIS PROBLEM IS U(S,T) = (S*T)**5
!
!     DEFINE THE INTEGERS M = 50 AND N = 63. THEN DEFINE THE
!     GRID INCREMENTS DELTAS = 1/(M+1) AND DELTAT = 1/(N+1).
!
!     THE GRID IS THEN GIVEN BY S(I) = I*DELTAS FOR I = 1,...,M
!     AND T(J) = J*DELTAT FOR J = 1,...,N.
!
!     THE APPROXIMATE SOLUTION IS GIVEN AS THE SOLUTION TO
!     THE FOLLOWING FINITE DIFFERENCE APPROXIMATION OF EQUATION (1).
!
!     .5/(S(I)*DELTAS)*((U(I+1,J)-U(I,J))/(2*S(I+.5)*DELTAS)
!                     -(U(I,J)-U(I-1,J))/(2*S(I-.5)*DELTAS))
!     +.5/(T(I)*DELTAT)*((U(I,J+1)-U(I,J))/(2*T(I+.5)*DELTAT) (2)
!                     -(U(I,J)-U(I,J-1))/(2*T(I-.5)*DELTAT))
!               = 15/4*S(I)*T(J)*(S(I)**4+T(J)**4)
!
!             WHERE S(I+.5) = .5*(S(I+1)+S(I))
!                   S(I-.5) = .5*(S(I)+S(I-1))
!                   T(I+.5) = .5*(T(I+1)+T(I))
!                   T(I-.5) = .5*(T(I)+T(I-1))
!
!     THE APPROACH IS TO WRITE EQUATION (2) IN THE FORM
!
!     AM(I)*U(I-1,J)+BM(I,J)*U(I,J)+CM(I)*U(I+1,J)
!       +AN(J)*U(I,J-1)+BN(J)*U(I,J)+CN(J)*U(I,J+1)      (3)
!           = Y(I,J)
!
!     AND THEN CALL SUBROUTINE BLKTRI TO DETERMINE U(I,J)
!
!
!
   logical function  TBLKTRI(DUMMY)


      TYPE ( fishworkspace) :: w
      integer, intent(in) :: DUMMY ! Does not mean anything
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: IFLG, NP, N, MP, M, IDIMY, I, J, IERROR
      REAL(fish_kind) , DIMENSION(75,105) :: Y
      REAL(fish_kind) , DIMENSION(75) :: AM, BM, CM
      REAL(fish_kind) , DIMENSION(105) :: AN, BN, CN
      REAL(fish_kind) , DIMENSION(75) :: S
      REAL(fish_kind) , DIMENSION(105) :: T
      REAL(fish_kind) ::DELTAS,DELTAT,HDS,TDS,TEMP1,TEMP2,TEMP3,HDT,TDT,ERR,Z
!-----------------------------------------------
      IFLG = 0
      NP = 1
      N = 63
      MP = 1
      M = 50
      IDIMY = 75
!
!     GENERATE AND STORE GRID POINTS FOR THE PURPOSE OF COMPUTING THE
!     COEFFICIENTS AND THE ARRAY Y.
!
      DELTAS = 1./FLOAT(M + 1)
      DO I = 1, M
         S(I) = FLOAT(I)*DELTAS
      END DO
      DELTAT = 1./FLOAT(N + 1)
      DO J = 1, N
         T(J) = FLOAT(J)*DELTAT
      END DO
!
!     COMPUTE THE COEFFICIENTS AM,BM,CM CORRESPONDING TO THE S DIRECTION
!
      HDS = DELTAS/2.
      TDS = DELTAS + DELTAS
      DO I = 1, M
         TEMP1 = 1./(S(I)*TDS)
         TEMP2 = 1./((S(I)-HDS)*TDS)
         TEMP3 = 1./((S(I)+HDS)*TDS)
         AM(I) = TEMP1*TEMP2
         CM(I) = TEMP1*TEMP3
         BM(I) = -(AM(I)+CM(I))
      END DO
!
!     COMPUTE THE COEFFICIENTS AN,BN,CN CORRESPONDING TO THE T DIRECTION
!
      HDT = DELTAT/2.
      TDT = DELTAT + DELTAT
      DO J = 1, N
         TEMP1 = 1./(T(J)*TDT)
         TEMP2 = 1./((T(J)-HDT)*TDT)
         TEMP3 = 1./((T(J)+HDT)*TDT)
         AN(J) = TEMP1*TEMP2
         CN(J) = TEMP1*TEMP3
         BN(J) = -(AN(J)+CN(J))
      END DO
!
!     COMPUTE RIGHT SIDE OF EQUATION
!
      DO J = 1, N
         Y(:M,J) = 3.75*S(:M)*T(J)*(S(:M)**4+T(J)**4)
      END DO
!
!     THE NONZERO BOUNDARY CONDITIONS ENTER THE LINEAR SYSTEM VIA
!     THE RIGHT SIDE Y(I,J). IF THE EQUATIONS (3) GIVEN ABOVE
!     ARE EVALUATED AT I=M AND J=1,...,N THEN THE TERM CM(M)*U(M+1,J)
!     IS KNOWN FROM THE BOUNDARY CONDITION TO BE CM(M)*T(J)**5.
!     THEREFORE THIS TERM CAN BE INCLUDED IN THE RIGHT SIDE Y(M,J).
!     THE SAME ANALYSIS APPLIES AT J=N AND I=1,..,M. NOTE THAT THE
!     CORNER AT J=N,I=M INCLUDES CONTRIBUTIONS FROM BOTH BOUNDARIES.
!
      Y(M,:N) = Y(M,:N) - CM(M)*T(:N)**5
      Y(:M,N) = Y(:M,N) - CN(N)*S(:M)**5
!
!     DETERMINE THE APPROXIMATE SOLUTION U(I,J)
!
      CALL BLKTRI(IFLG,NP,N,AN,BN,CN,MP,M,AM,BM,CM,IDIMY,Y,IERROR,W)
      IFLG = IFLG + 1
      DO WHILE(IFLG - 1 <= 0)
         CALL BLKTRI (IFLG, NP, N, AN, BN, CN, MP, M, AM, BM, CM, IDIMY, Y, IERROR, W)
         IFLG = IFLG + 1
      END DO
      ERR = 0.
      DO J = 1, N
         DO I = 1, M
            Z = ABS(Y(I,J)-(S(I)*T(J))**5)
            ERR = AMAX1(Z,ERR)
         END DO
      END DO
!     Print earlier output from platforms with 32 and 64 bit floating point
!     arithemtic followed by the output from this computer
      call msg('    BLKTRI TEST RUN *** ')
      call msg('    Previous 64 bit floating point arithmetic result ')
      call msg('    IERROR = 0,  Discretization Error = 1.6478E-05')
      call msg('    Previous 32 bit floating point arithmetic result ')
      call msg('    IERROR = 0,  Discretization Error = 1.2737E-02')
      call msg('    The output from your computer is: ')
      call msg('    IERROR ='+fu_str( IERROR) +', Discretization Error = ',ERR)
!     release dynamically allocated work space
      CALL FISHFIN (W)
      if (fish_kind == 4) then 
            TBLKTRI = (ERR<1.3E-02)
      elseif (fish_kind == 8) then 
            TBLKTRI = (ERR<1.7E-05)
      else
           call msg("Unknown fish_kind",fish_kind)
           TBLKTRI = .false.
      endif
      return 
   END function TBLKTRI



!
!     file pois3d.f
!
!     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
!     *                                                               *
!     *                  copyright (c) 2005 by UCAR                   *
!     *                                                               *
!     *       University Corporation for Atmospheric Research         *
!     *                                                               *
!     *                      all rights reserved                      *
!     *                                                               *
!     *                    FISHPACK90  version 1.1                    *
!     *                                                               *
!     *                 A Package of Fortran 77 and 90                *
!     *                                                               *
!     *                Subroutines and Example Programs               *
!     *                                                               *
!     *               for Modeling Geophysical Processes              *
!     *                                                               *
!     *                             by                                *
!     *                                                               *
!     *        John Adams, Paul Swarztrauber and Roland Sweet         *
!     *                                                               *
!     *                             of                                *
!     *                                                               *
!     *         the National Center for Atmospheric Research          *
!     *                                                               *
!     *                Boulder, Colorado  (80307)  U.S.A.             *
!     *                                                               *
!     *                   which is sponsored by                       *
!     *                                                               *
!     *              the National Science Foundation                  *
!     *                                                               *
!     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
!
!     SUBROUTINE POIS3D (LPEROD,L,C1,MPEROD,M,C2,NPEROD,N,A,B,C,LDIMF,
!    +                   MDIMF,F,IERROR)
!
!
! DIMENSION OF           A(N), B(N), C(N), F(LDIMF,MDIMF,N)
! ARGUMENTS
!
! LATEST REVISION        June 2004
!
! PURPOSE                SOLVES THE LINEAR SYSTEM OF EQUATIONS
!                        FOR UNKNOWN X VALUES, WHERE I=1,2,...,L,
!                        J=1,2,...,M, AND K=1,2,...,N
!
!                        C1*(X(I-1,J,K) -2.*X(I,J,K) +X(I+1,J,K)) +
!                        C2*(X(I,J-1,K) -2.*X(I,J,K) +X(I,J+1,K)) +
!                        A(K)*X(I,J,K-1) +B(K)*X(I,J,K)+ C(K)*X(I,J,K+1)
!                        = F(I,J,K)
!
!                        THE INDICES K-1 AND K+1 ARE EVALUATED MODULO N,
!                        I.E. X(I,J,0)=X(I,J,N) AND X(I,J,N+1)=X(I,J,1).
!                        THE UNKNOWNS
!                        X(0,J,K), X(L+1,J,K), X(I,0,K), AND X(I,M+1,K)
!                        ARE ASSUMED TO TAKE ON CERTAIN PRESCRIBED
!                        VALUES DESCRIBED BELOW.
!
! USAGE                  CALL POIS3D (LPEROD,L,C1,MPEROD,M,C2,NPEROD,
!                        N,A,B,C,LDIMF,MDIMF,F,IERROR)
!
! ARGUMENTS
!
! ON INPUT
!                        LPEROD
!                          INDICATES THE VALUES THAT X(0,J,K) AND
!                          X(L+1,J,K) ARE ASSUMED TO HAVE.
!                          = 0  X(0,J,K)=X(L,J,K), X(L+1,J,K)=X(1,J,K)
!                          = 1  X(0,J,K) = 0,      X(L+1,J,K) = 0
!                          = 2  X(0,J,K)=0,        X(L+1,J,K)=X(L-1,J,K)
!                          = 3  X(0,J,K)=X(2,J,K), X(L+1,J,K)=X(L-1,J,K)
!                          = 4  X(0,J,K)=X(2,J,K), X(L+1,J,K) = 0.
!
!                        L
!                          THE NUMBER OF UNKNOWNS IN THE I-DIRECTION.
!                          L MUST BE AT LEAST 3.
!
!                        C1
!                          REAL CONSTANT IN THE ABOVE LINEAR SYSTEM
!                          OF EQUATIONS TO BE SOLVED.
!
!                        MPEROD
!                          INDICATES THE VALUES THAT X(I,0,K) AND
!                          X(I,M+1,K) ARE ASSUMED TO HAVE.
!                          = 0  X(I,0,K)=X(I,M,K), X(I,M+1,K)=X(I,1,K)
!                          = 1  X(I,0,K)=0,        X(I,M+1,K)=0
!                          = 2  X(I,0,K)=0,        X(I,M+1,K)=X(I,M-1,K)
!                          = 3  X(I,0,K)=X(I,2,K)  X(I,M+1,K)=X(I,M-1,K)
!                          = 4  X(I,0,K)=X(I,2,K)  X(I,M+1,K)=0
!
!                        M
!                          THE NUMBER OF UNKNOWNS IN THE J-DIRECTION.
!                          M MUST BE AT LEAST 3.
!
!                        C2
!                          REAL CONSTANT IN THE ABOVE LINEAR SYSTEM
!                          OF EQUATIONS TO BE SOLVED.
!
!                        NPEROD
!                          = 0  IF A(1) AND C(N) ARE NOT ZERO.
!                          = 1  IF A(1) = C(N) = 0.
!
!                        N
!                          THE NUMBER OF UNKNOWNS IN THE K-DIRECTION.
!                          N MUST BE AT LEAST 3.
!
!                        A, B, C
!                          ONE-DIMENSIONAL ARRAYS OF LENGTH N THAT
!                          SPECIFY THE COEFFICIENTS IN THE LINEAR
!                          EQUATIONS GIVEN ABOVE.
!
!                          IF NPEROD = 0 THE ARRAY ELEMENTS MUST NOT
!                          DEPEND UPON INDEX K, BUT MUST BE CONSTANT.
!                          SPECIFICALLY,THE SUBROUTINE CHECKS THE
!                          FOLLOWING CONDITION
!                            A(K) = C(1)
!                            C(K) = C(1)
!                            B(K) = B(1)
!                          FOR K=1,2,...,N.
!
!                        LDIMF
!                          THE ROW (OR FIRST) DIMENSION OF THE THREE-
!                          DIMENSIONAL ARRAY F AS IT APPEARS IN THE
!                          PROGRAM CALLING POIS3D.  THIS PARAMETER IS
!                          USED TO SPECIFY THE VARIABLE DIMENSION
!                          OF F.  LDIMF MUST BE AT LEAST L.
!
!                        MDIMF
!                          THE COLUMN (OR SECOND) DIMENSION OF THE THREE
!                          DIMENSIONAL ARRAY F AS IT APPEARS IN THE
!                          PROGRAM CALLING POIS3D.  THIS PARAMETER IS
!                          USED TO SPECIFY THE VARIABLE DIMENSION
!                          OF F.  MDIMF MUST BE AT LEAST M.
!
!                        F
!                          A THREE-DIMENSIONAL ARRAY THAT SPECIFIES THE
!                          VALUES OF THE RIGHT SIDE OF THE LINEAR SYSTEM
!                          OF EQUATIONS GIVEN ABOVE.  F MUST BE
!                          DIMENSIONED AT LEAST L X M X N.
!
! ON OUTPUT
!
!                        F
!                          CONTAINS THE SOLUTION X.
!
!                        IERROR
!                          AN ERROR FLAG THAT INDICATES INVALID INPUT
!                          PARAMETERS.  EXCEPT FOR NUMBER ZERO, A
!                          SOLUTION IS NOT ATTEMPTED.
!                          = 0  NO ERROR
!                          = 1  IF LPEROD .LT. 0 OR .GT. 4
!                          = 2  IF L .LT. 3
!                          = 3  IF MPEROD .LT. 0 OR .GT. 4
!                          = 4  IF M .LT. 3
!                          = 5  IF NPEROD .LT. 0 OR .GT. 1
!                          = 6  IF N .LT. 3
!                          = 7  IF LDIMF .LT. L
!                          = 8  IF MDIMF .LT. M
!                          = 9  IF A(K) .NE. C(1) OR C(K) .NE. C(1)
!                               OR B(I) .NE.B(1) FOR SOME K=1,2,...,N.
!                          = 10 IF NPEROD = 1 AND A(1) .NE. 0
!                               OR C(N) .NE. 0
!                          = 20 If the dynamic allocation of real and
!                               complex work space required for solution
!                               fails (for example if N,M are too large
!                               for your computer)
!
!                          SINCE THIS IS THE ONLY MEANS OF INDICATING A
!                          POSSIBLY INCORRECT CALL TO POIS3D, THE USER
!                          SHOULD TEST IERROR AFTER THE CALL.
!
! SPECIAL CONDITIONS     NONE
!
! I/O                    NONE
!
! PRECISION              SINGLE
!
! REQUIRED files         fish.f,comf.f,fftpack.f
!
! LANGUAGE               FORTRAN 90
!
! HISTORY                WRITTEN BY ROLAND SWEET AT NCAR IN THE LATE
!                        1970'S.  RELEASED ON NCAR'S PUBLIC SOFTWARE
!                        LIBRARIES IN JANUARY, 1980.
!                        Revised in June 2004 by John Adams using
!                        Fortran 90 dynamically allocated work space.
!
! PORTABILITY            FORTRAN 90
!
! ALGORITHM              THIS SUBROUTINE SOLVES THREE-DIMENSIONAL BLOCK
!                        TRIDIAGONAL LINEAR SYSTEMS ARISING FROM FINITE
!                        DIFFERENCE APPROXIMATIONS TO THREE-DIMENSIONAL
!                        POISSON EQUATIONS USING THE FFT PACKAGE
!                        FFTPACK WRITTEN BY PAUL SWARZTRAUBER.
!
! TIMING                 FOR LARGE L, M AND N, THE OPERATION COUNT
!                        IS ROUGHLY PROPORTIONAL TO
!                          L*M*N*(LOG2(L)+LOG2(M)+5)
!                        BUT ALSO DEPENDS ON INPUT PARAMETERS LPEROD
!                        AND MPEROD.
!
! ACCURACY               TO MEASURE THE ACCURACY OF THE ALGORITHM A
!                        UNIFORM RANDOM NUMBER GENERATOR WAS USED TO
!                        CREATE A SOLUTION ARRAY X FOR THE SYSTEM GIVEN
!                        IN THE 'PURPOSE' SECTION WITH
!                          A(K) = C(K) = -0.5*B(K) = 1,  K=1,2,...,N
!                        AND, WHEN NPEROD = 1
!                          A(1) = C(N) = 0
!                          A(N) = C(1) = 2.
!
!                        THE SOLUTION X WAS SUBSTITUTED INTO THE GIVEN
!                        SYSTEM AND, USING DOUBLE PRECISION, A RIGHT
!                        SIDE Y WAS COMPUTED.  USING THIS ARRAY Y
!                        SUBROUTINE POIS3D WAS CALLED TO PRODUCE AN
!                        APPROXIMATE SOLUTION Z.  RELATIVE ERROR
!
!                        E = MAX(ABS(Z(I,J,K)-X(I,J,K)))/MAX(ABS(X(I,J,K
!
!                        WAS COMPUTED, WHERE THE TWO MAXIMA ARE TAKEN
!                        OVER I=1,2,...,L, J=1,2,...,M AND K=1,2,...,N.
!                        VALUES OF E ARE GIVEN IN THE TABLE BELOW FOR
!                        SOME TYPICAL VALUES OF L,M AND N.
!
!                        L(=M=N)   LPEROD    MPEROD       E
!                        ------    ------    ------     ------
!
!                          16        0         0        1.E-13
!                          15        1         1        4.E-13
!                          17        3         3        2.E-13
!                          32        0         0        2.E-13
!                          31        1         1        2.E-12
!                          33        3         3        7.E-13
!
! REFERENCES              NONE
! ********************************************************************
      SUBROUTINE POIS3D(LPEROD, L, C1, MPEROD, M, C2, NPEROD, N, A, B, C, LDIMF, MDIMF, F, IERROR)
!               USE fish
      implicit none
      TYPE (fishworkspace) :: w
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER  :: LPEROD
      INTEGER  :: L
      INTEGER  :: MPEROD
      INTEGER  :: M
      INTEGER  :: NPEROD
      INTEGER  :: N
      INTEGER  :: LDIMF
      INTEGER  :: MDIMF
      INTEGER  :: IERROR
      REAL(fish_kind)  :: C1
      REAL(fish_kind)  :: C2
      REAL(fish_kind)  :: A(*)
      REAL(fish_kind)  :: B(*)
      REAL(fish_kind)  :: C(*)
      REAL(fish_kind)  :: F(LDIMF,MDIMF,*)
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: LP, MP, NP, K, IRWK, ICWK
      REAL, DIMENSION(6) :: SAVE
!-----------------------------------------------
      LP = LPEROD + 1
      MP = MPEROD + 1
      NP = NPEROD + 1
!
!     CHECK FOR INVALID INPUT.
!
      IERROR = 0
      IF (LP<1 .OR. LP>5) IERROR = 1
      IF (L < 3) IERROR = 2
      IF (MP<1 .OR. MP>5) IERROR = 3
      IF (M < 3) IERROR = 4
      IF (NP<1 .OR. NP>2) IERROR = 5
      IF (N < 3) IERROR = 6
      IF (LDIMF < L) IERROR = 7
      IF (MDIMF < M) IERROR = 8
      IF (NP == 1) THEN
         DO K = 1, N
            IF (A(K) /= C(1)) GO TO 102
            IF (C(K) /= C(1)) GO TO 102
            IF (B(K) /= B(1)) GO TO 102
         END DO
         GO TO 104
  102    CONTINUE
         IERROR = 9
      ENDIF
      IF (NPEROD==1 .AND. (A(1)/=0. .OR. C(N)/=0.)) IERROR = 10
! 104 IF (IERROR .NE. 0) GO TO 122
  104 CONTINUE
      IF (IERROR /= 0) RETURN 
!     allocate required work space length (generous estimate)
      IRWK=30+L+M+2*N+MAX0(L,M,N)+7*(INT((L+1)/2)+INT((M+1)/2))
      ICWK = 0
      CALL ALLOCATFISH (IRWK, ICWK, W, IERROR)
!     check that allocation was successful
      IF (IERROR == 20) RETURN 
      call pois3dd(LPEROD,L,C1,MPEROD,M,C2,NPEROD,N,A,B,C,LDIMF, MDIMF,F,IERROR,w%rew)
!     release work space
      CALL FISHFIN (W)
      RETURN 
      END SUBROUTINE POIS3D


 
      SUBROUTINE POIS3DD(LPEROD, L, C1, MPEROD, M, C2, NPEROD, N, A, B,  C, LDIMF, MDIMF, F, IERROR, W)
      implicit none
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER , INTENT(IN) :: LPEROD
      INTEGER  :: L
      INTEGER , INTENT(IN) :: MPEROD
      INTEGER  :: M
      INTEGER , INTENT(IN) :: NPEROD
      INTEGER  :: N
      INTEGER  :: LDIMF
      INTEGER  :: MDIMF
      INTEGER  :: IERROR
      REAL(fish_kind)  :: C1
      REAL(fish_kind)  :: C2
      REAL(fish_kind)  :: A(*)
      REAL(fish_kind)  :: B(*)
      REAL(fish_kind)  :: C(*)
      REAL(fish_kind)  :: F(LDIMF,MDIMF,*)
      REAL (fish_kind)  :: W(*)
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: LP, MP, NP, IWYRT, IWT, IWD, IWBB, IWX, IWY, NH, NHM1, & 
                     & NODD, I, J, K, NHPK, NHMK
      REAL, DIMENSION(6) :: SAVE
!-----------------------------------------------
      LP = LPEROD + 1
      MP = MPEROD + 1
      NP = NPEROD + 1
      IWYRT = L + 1
      IWT = IWYRT + M
      IWD = IWT + MAX0(L,M,N) + 1
      IWBB = IWD + N
      IWX = IWBB + N
      IWY = IWX + 7*((L + 1)/2) + 15
      GO TO (105,114) NP
!
!     REORDER UNKNOWNS WHEN NPEROD = 0.
!
  105 CONTINUE
      NH = (N + 1)/2
      NHM1 = NH - 1
      NODD = 1
      IF (2*NH == N) NODD = 2
      DO I = 1, L
         DO J = 1, M
            DO K = 1, NHM1
               W(K) = F(I,J,NH-K) - F(I,J,K+NH)
               W(K+NH) = F(I,J,NH-K) + F(I,J,K+NH)
            END DO
            W(NH) = 2.*F(I,J,NH)
            GO TO (108,107) NODD
  107       CONTINUE
            W(N) = 2.*F(I,J,N)
  108       CONTINUE
            F(I,J,:N) = W(:N)
         END DO
      END DO
      SAVE(1) = C(NHM1)
      SAVE(2) = A(NH)
      SAVE(3) = C(NH)
      SAVE(4) = B(NHM1)
      SAVE(5) = B(N)
      SAVE(6) = A(N)
      C(NHM1) = 0.
      A(NH) = 0.
      C(NH) = 2.*C(NH)
      SELECT CASE (NODD) 
      CASE DEFAULT
         B(NHM1) = B(NHM1) - A(NH-1)
         B(N) = B(N) + A(N)
      CASE (2) 
         A(N) = C(NH)
      END SELECT
  114 CONTINUE
      CALL POS3D1 (LP, L, MP, M, N, A, B, C, LDIMF, MDIMF, F, W, W(IWYRT), W(IWT), W(IWD), W(IWX), W(IWY), C1, C2, W(IWBB))
      GO TO (115,122) NP
  115 CONTINUE
      DO I = 1, L
         DO J = 1, M
            W(NH-1:NH-NHM1:(-1))=0.5*(F(I,J,NH+1:NHM1+NH)+F(I,J,:NHM1))
            W(NH+1:NHM1+NH) = 0.5*(F(I,J,NH+1:NHM1+NH)-F(I,J,:NHM1))
            W(NH) = 0.5*F(I,J,NH)
            GO TO (118,117) NODD
  117       CONTINUE
            W(N) = 0.5*F(I,J,N)
  118       CONTINUE
            F(I,J,:N) = W(:N)
         END DO
      END DO
      C(NHM1) = SAVE(1)
      A(NH) = SAVE(2)
      C(NH) = SAVE(3)
      B(NHM1) = SAVE(4)
      B(N) = SAVE(5)
      A(N) = SAVE(6)
  122 CONTINUE
      RETURN 
      END SUBROUTINE POIS3DD


      SUBROUTINE POS3D1(LP, L, MP, M, N, A, B, C, LDIMF, MDIMF, F, XRT, YRT, T, D, WX, WY, C1, C2, BB)
      implicit none
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER , INTENT(IN) :: LP
      INTEGER , INTENT(IN) :: L
      INTEGER , INTENT(IN) :: MP
      INTEGER , INTENT(IN) :: M
      INTEGER , INTENT(IN) :: N
      INTEGER , INTENT(IN) :: LDIMF
      INTEGER , INTENT(IN) :: MDIMF
      REAL(fish_kind) , INTENT(IN) :: C1
      REAL(fish_kind) , INTENT(IN) :: C2
      REAL(fish_kind)  :: A(*)
      REAL(fish_kind) , INTENT(IN) :: B(*)
      REAL(fish_kind)  :: C(*)
      REAL(fish_kind) , INTENT(INOUT) :: F(LDIMF,MDIMF,1)
      REAL(fish_kind) , INTENT(INOUT) :: XRT(*)
      REAL(fish_kind) , INTENT(INOUT) :: YRT(*)
      REAL(fish_kind)  :: T(*)
      REAL(fish_kind)  :: D(*)
      REAL(fish_kind)  :: WX(*)
      REAL(fish_kind)  :: WY(*)
      REAL(fish_kind)  :: BB(*)
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: LR, MR, NR, LRDEL, I, MRDEL, J, IFWRD, IS, K
      REAL :: PI, DUM, SCALX, DX, DI, SCALY, DY, DJ
!-----------------------------------------------
      PI = 4.0*ATAN(1.0)
      LR = L
      MR = M
      NR = N
!
!     GENERATE TRANSFORM ROOTS
!
      LRDEL = ((LP - 1)*(LP - 3)*(LP - 5))/3
      SCALX = LR + LRDEL
      DX = PI/(2.*SCALX)
      GO TO (108,103,101,102,101) LP
  101 CONTINUE
      DI = 0.5
      SCALX = 2.*SCALX
      GO TO 104
  102 CONTINUE
      DI = 1.0
      GO TO 104
  103 CONTINUE
      DI = 0.0
  104 CONTINUE
      DO I = 1, LR
         XRT(I) = -4.*C1*SIN((FLOAT(I) - DI)*DX)**2
      END DO
      SCALX = 2.*SCALX
      GO TO (112,106,110,107,111) LP
  106 CONTINUE
      CALL SINTI (LR, WX)
      GO TO 112
  107 CONTINUE
      CALL COSTI (LR, WX)
      GO TO 112
  108 CONTINUE
      XRT(1) = 0.
      XRT(LR) = -4.*C1
      DO I = 3, LR, 2
         XRT(I-1) = -4.*C1*SIN(FLOAT(I - 1)*DX)**2
         XRT(I) = XRT(I-1)
      END DO
      CALL RFFTI (LR, WX)
      GO TO 112
  110 CONTINUE
      CALL SINQI (LR, WX)
      GO TO 112
  111 CONTINUE
      CALL COSQI (LR, WX)
  112 CONTINUE
      MRDEL = ((MP - 1)*(MP - 3)*(MP - 5))/3
      SCALY = MR + MRDEL
      DY = PI/(2.*SCALY)
      GO TO (120,115,113,114,113) MP
  113 CONTINUE
      DJ = 0.5
      SCALY = 2.*SCALY
      GO TO 116
  114 CONTINUE
      DJ = 1.0
      GO TO 116
  115 CONTINUE
      DJ = 0.0
  116 CONTINUE
      DO J = 1, MR
         YRT(J) = -4.*C2*SIN((FLOAT(J) - DJ)*DY)**2
      END DO
      SCALY = 2.*SCALY
      GO TO (124,118,122,119,123) MP
  118 CONTINUE
      CALL SINTI (MR, WY)
      GO TO 124
  119 CONTINUE
      CALL COSTI (MR, WY)
      GO TO 124
  120 CONTINUE
      YRT(1) = 0.
      YRT(MR) = -4.*C2
      DO J = 3, MR, 2
         YRT(J-1) = -4.*C2*SIN(FLOAT(J - 1)*DY)**2
         YRT(J) = YRT(J-1)
      END DO
      CALL RFFTI (MR, WY)
      GO TO 124
  122 CONTINUE
      CALL SINQI (MR, WY)
      GO TO 124
  123 CONTINUE
      CALL COSQI (MR, WY)
  124 CONTINUE
      IFWRD = 1
      IS = 1
  125 CONTINUE
!
!     TRANSFORM X
!
      DO J=1,MR
         DO K=1,NR
            DO I=1,LR
                  T(I) = F(I,J,K)
            END DO
            GO TO (127,130,131,134,135),LP
  127       GO TO (128,129),IFWRD
  128       CALL RFFTF (LR,T,WX)
            GO TO 138
  129       CALL RFFTB (LR,T,WX)
            GO TO 138
  130       CALL SINT (LR,T,WX)
            GO TO 138
  131       GO TO (132,133),IFWRD
  132       CALL SINQF (LR,T,WX)
            GO TO 138
  133       CALL SINQB (LR,T,WX)
            GO TO 138
  134       CALL COST (LR,T,WX)
            GO TO 138
  135       GO TO (136,137),IFWRD
  136       CALL COSQF (LR,T,WX)
            GO TO 138
  137       CALL COSQB (LR,T,WX)
  138       CONTINUE
             DO I=1,LR
               F(I,J,K) = T(I)
             END DO
         END DO
      END DO
      GO TO (142,164) IFWRD
!
!     TRANSFORM Y
!
  142 CONTINUE
      DO I=1,LR
	 DO K=1,NR
	    DO J=1,MR
               T(J) = F(I,J,K)
	    END DO
            GO TO (144,147,148,151,152),MP
  144       GO TO (145,146),IFWRD
  145       CALL RFFTF (MR,T,WY)
            GO TO 155
  146       CALL RFFTB (MR,T,WY)
            GO TO 155
  147       CALL SINT (MR,T,WY)
            GO TO 155
  148       GO TO (149,150),IFWRD
  149       CALL SINQF (MR,T,WY)
            GO TO 155
  150       CALL SINQB (MR,T,WY)
            GO TO 155
  151       CALL COST (MR,T,WY)
            GO TO 155
  152       GO TO (153,154),IFWRD
  153       CALL COSQF (MR,T,WY)
            GO TO 155
  154       CALL COSQB (MR,T,WY)
  155       CONTINUE
	    DO J=1,MR
               F(I,J,K) = T(J)
	    END DO
	 END DO
      END DO
      GO TO (159,125) IFWRD
  159 CONTINUE
      DO I = 1, LR
         DO J = 1, MR
            BB(:NR) = B(:NR) + XRT(I) + YRT(J)
            T(:NR) = F(I,J,:NR)
            CALL TRID (NR, A, BB, C, T, D)
            F(I,J,:NR) = T(:NR)
         END DO
      END DO
      IFWRD = 2
      IS = -1
      GO TO 142
  164 CONTINUE
      F(:LR,:MR,:NR) = F(:LR,:MR,:NR)/(SCALX*SCALY)
      RETURN 
      END SUBROUTINE POS3D1


      SUBROUTINE TRID(MR, A, B, C, Y, D)
      implicit none
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER , INTENT(IN) :: MR
      REAL(fish_kind) , INTENT(IN) :: A(*)
      REAL(fish_kind) , INTENT(IN) :: B(*)
      REAL(fish_kind) , INTENT(IN) :: C(*)
      REAL(fish_kind) , INTENT(INOUT) :: Y(*)
      REAL(fish_kind) , INTENT(INOUT) :: D(*)
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: M, MM1, I, IP
      REAL :: Z
!-----------------------------------------------
      M = MR
      MM1 = M - 1
      Z = 1./B(1)
      D(1) = C(1)*Z
      Y(1) = Y(1)*Z
      DO I = 2, MM1
         Z = 1./(B(I)-A(I)*D(I-1))
         D(I) = C(I)*Z
         Y(I) = (Y(I)-A(I)*Y(I-1))*Z
      END DO
      Z = B(M) - A(M)*D(MM1)
      IF (Z == 0.) THEN
         Y(M) = 0.
      ELSE
         Y(M) = (Y(M)-A(M)*Y(MM1))/Z
      ENDIF
      DO IP = 1, MM1
         I = M - IP
         Y(I) = Y(I) - D(I)*Y(I+1)
      END DO
      RETURN 
!
! REVISION HISTORY---
!
! SEPTEMBER 1973    VERSION 1
! APRIL     1976    VERSION 2
! JANUARY   1978    VERSION 3
! DECEMBER  1979    VERSION 3.1
! FEBRUARY  1985    DOCUMENTATION UPGRADE
! NOVEMBER  1988    VERSION 3.2, FORTRAN 77 CHANGES
! June      2004    Version 5.0, Fortran 90 Changes
!-----------------------------------------------------------------------
      END SUBROUTINE TRID

   end module fishpack_silam
