!+
SUBROUTINE CUBIC (T,P,ROOTS,NROOTS)
! ---------------------------------------------------------------------------
! PURPOSE - SOLVES FOR DENSITY IN REDUCED VAN DER WAALS GAS:
!        D**3 -3*D**2 +((8*T+P)/3)*D -P = 0                             
!     IF MORE THAN ONE DISTINCT ROOT EXISTS, ONLY THE MINIMUM AND       
!     MAXIMUM ONES ARE RETURNED AND NROOTS=2.                           

!      WRITTEN BY: THEODORE E. FESSLER, NASA TM X-3572, AUGUST 1977
!      REWRITTEN BY: LT. MARK D. KLEM & MARGARET P. PROCTOR, JANUARY 198
!                                                                       
  IMPLICIT NONE

  REAL(KIND=4),INTENT(IN):: t,p
  REAL(KIND=4),DIMENSION(:),INTENT(OUT):: roots   ! dimension=2
  INTEGER,INTENT(OUT):: nroots
                                                                        

                                                                        
                                                                        
!!      IMPLICIT REAL*8 (A-H,O-Z) 
!!      REAL*8 K1,K2 
!!      REAL*4 T,P,ROOTS(2) 
  REAL(KIND=8),PARAMETER:: K1=2.094395102393195D+0
  REAL(KIND=8),PARAMETER:: K2=4.188790204786391D+0
  REAL(KIND=8),PARAMETER:: ZERO=0.0, ONE=1.0, TWO=2.0, THREE=3.0, EIGHT=8.0

!!!      DATA THIRD / 3.333333333333333D-1 / 
  REAL(KIND=8):: a,a3,aa,bb,b,b2,cr
  REAL(KIND=8):: rad,phi
  REAL(KIND=8):: q
  REAL(KIND=8):: r1,r2,r3
                                                                        
!.....MAIN PROGRAM FLOW.                                                
!     PROCEDURE (REDUCE TO NORMAL FORM)                                 
      Q=(EIGHT*T+P)/THREE
      A=Q/THREE-ONE
      A3=A*A*A 
      B=(Q-P)/TWO-ONE
      B2=B*B 
      RAD=A3+B2 
      IF (RAD < ZERO) THEN 
!        PROCEDURE (CASE FOR 3 DISTINCT REAL ROOTS)                     
         PHI=ACOS(SIGN(SQRT(-B2/A3),-B))/THREE
         CR=TWO*SQRT(-A) 
         R1=CR*COS(PHI) 
         R2=CR*COS(PHI+K1) 
         R3=CR*COS(PHI+K2) 
         ROOTS(1)=MIN(R1,R2,R3)+ONE
         ROOTS(2)=MAX(R1,R2,R3)+ONE
         NROOTS=2 
      ELSE 
         IF (RAD .EQ. ZERO) THEN 
            IF (B .EQ. ZERO) THEN 
!              PROCEDURE (CASE FOR ONLY ONE REAL ROOT)                  
               IF (RAD .NE. ZERO) RAD=SQRT(RAD) 
               AA=-B+RAD 
               BB=-B-RAD 
               IF (AA .NE. ZERO) AA=QRT(AA) 
               IF (BB .NE. ZERO) BB=QRT(BB) 
               ROOTS(1)=AA+BB+ONE
               NROOTS=1 
            ELSE 
!              PROCEDURE (CASE FOR 3 REAL ROOTS BUT TWO ARE EQUAL)      
               R1=QRT(B) 
               R2=-(R1+R1)
               ROOTS(1)=MIN(R1,R2)+ONE
               ROOTS(2)=MAX(R1,R2)+ONE
               NROOTS=2 
            END IF 
         ELSE 
!           PROCEDURE (CASE FOR ONLY ONE REAL ROOT)                     
            IF (RAD.NE.0.0D0) RAD=SQRT(RAD) 
            AA=-B+RAD 
            BB=-B-RAD 
            IF (AA .NE. ZERO) AA=QRT(AA) 
            IF (BB .NE. ZERO) BB=QRT(BB) 
            ROOTS(1)=AA+BB+ONE
            NROOTS=1 
         END IF 
      END IF 

  RETURN

CONTAINS
!.....FUNCTION DEFINITION.
!+
FUNCTION Qrt(arg) RESULT(f)
! ---------------------------------------------------------------------------
  REAL(KIND=8),INTENT(IN):: arg
  REAL(KIND=8):: f
  REAL(KIND=8),PARAMETER:: THIRD=1.0D0/3.0D0
!----------------------------------------------------------------------------
  f=SIGN(ABS(ARG)**THIRD, ARG) 
  RETURN
END Function Qrt   ! ========================================================

END Subroutine Cubic   ! ====================================================
