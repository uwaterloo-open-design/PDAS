!+
SUBROUTINE NTRP (A,X,N0,IA,OUT,IGOTO,F,FF,VALUE)
! ---------------------------------------------------------------------------
! PURPOSE - DOUBLE 3-POINT INTERPOLATION ROUTINES FOR SUBPROGRAM FLUID.
!      WRITTEN BY: THEODORE E. FESSLER, NASA TM X-3572, AUGUST 1977
!      REWRITTEN BY: LT. MARK D. KLEM & MARGARET P. PROCTOR, JANUARY 198
!                                                                       


  IMPLICIT NONE
  REAL,INTENT(IN),DIMENSION(:):: a   ! TABLE OF INDEPENDENT-VARIABLE VALUES
  REAL,INTENT(IN):: x   ! ENTRY VALUE OF INDEPENDENT VARIABLE
  INTEGER,INTENT(IN):: n0   !  (OR N1/N2) = LENGTH OF A TABLE
  INTEGER,INTENT(OUT):: ia   !  CURRENT INDEX VALUE ASSOCIATED WITH A TABLE
  INTEGER,INTENT(OUT):: out  !  = 0,1 IF X IS INSIDE/NOT INSIDE LIMITS OF A
  INTEGER,INTENT(IN):: igoto ! ???
  REAL,INTENT(IN),DIMENSION(:):: f
  REAL,INTENT(IN),DIMENSION(:,:):: ff
  REAL,INTENT(OUT):: value

!      THE ORIGINAL PROGRAM WAS WRITTEN IN SFTRAN, BUT THIS PROGRAM IS  
!      REWRITTEN IN STANDARD FORTRAN 77 TO RUN PRIMARILY ON PERSONAL    
!      COMPUTERS AND CAN ALSO RUN ON MAINFRAMES. MUCH OF THE OLD CODE WA
!      COMMENTED OUT BUT LEFT INTACT FOR DOCUMENTATION AND REFERENCE    
!      PURPOSES.

! NOTES - This routine is an example of a situation that requires the
!   SAVE attribute of Fortran 90. The routine must first be called with
!   IGOTO = 1,2, or 3. Then, a later call with IGOTO=4 or 5 actually
!   produces the result. BUT, certain values must be saved between
!   these calls in order for this to work. Most Fortran IV compilers of
!   the early days saved all the variables and many programs were written
!   with this as a basic assumption. Since modern Fortran uses variables
!   with the AUTOMATIC attribute, routines such as this will fail
!   unless the appropriate variables are declared SAVE.

  REAL,SAVE,DIMENSION(4,3):: c
  REAL,SAVE,DIMENSION(4):: c0,c1,c2     ! SAVE is redundant (EQUIVALENCE)
  REAL,SAVE,DIMENSION(16):: cc 
  INTEGER,SAVE,DIMENSION(3):: i1,i2,nn
  LOGICAL,SAVE:: NEW12

  INTEGER:: i,j,ij,k,l,m   ! various counters
  INTEGER:: l1,l2, m1,m2, mm1,mm2
  INTEGER:: na
  INTEGER:: ni,nj

!!!      REAL F(40) 
!!!      REAL FF(27,30)


      EQUIVALENCE (c0(1),c(1,1)), (c1(1),c(1,2)), (c2,c(1,3))
      ! so, c0 is col. 1 of c, c1 is col. 2, c2 is col.3

      EQUIVALENCE (L1,i1(1)), (m1,i1(2)), (mm1,i1(3))
      EQUIVALENCE (L2,i2(1)), (m2,i2(2)), (mm2,i2(3))
      EQUIVALENCE             (ni,nn(2)), (nj,nn(3))

  REAL:: p1,p2,p3
  REAL:: a11,a12,a13,a14, a22,a23,a24, a33,a34, a44
!----------------------------------------------------------------------------                                                                        

  IF (IGOTO .EQ. 4) THEN      ! function of one variable
    VALUE=0.0 
    K=0 
    DO I=L1,L2 
      K=K+1 
      VALUE=VALUE+C0(K)*F(I)
    END DO
    RETURN 
  END IF                                                                      

  IF (IGOTO .EQ. 5) THEN 
    IF (NEW12) THEN 
      IJ=0 
      DO I=1,NI 
        DO J=1,NJ 
          IJ=IJ+1 
          CC(IJ)=C1(I)*C2(J) 
        END DO   
      END DO   
      NEW12=.FALSE. 
    END IF

    VALUE=0.0
    IJ=0 
    DO I=M1,M2 
      DO J=MM1,MM2 
        IJ=IJ+1 
        VALUE=VALUE+CC(IJ)*FF(I,J) 
      END DO   
    END DO   
    RETURN 
  END IF 

  IF (IGOTO .EQ. 1) THEN 
    K=1 
    NA=N0 
  ELSE IF (IGOTO .EQ. 2) THEN 
    K=2 
    NA=N0 
    NEW12=.TRUE. 
  ELSE IF (IGOTO .EQ. 3) THEN 
    K=3 
    NA=N0 
    NEW12=.TRUE. 
  END IF 


  L=Limit(1,NA,3) 
  M=NA-2 
  IF (X .LT. A(1)) THEN     !     PROCEDURE (TABLE LOOK-UP)                                         
    OUT=1 
    IA=1 
  ELSE IF (X .GT. A(NA)) THEN 
    OUT=1 
    IA=NA 
  ELSE 
    OUT=0 
    IA=LIMIT(1,IA,NA) 
    IF (X .LT. A(IA)) THEN 
10    IA=IA-1                         ! look backwards
      IF (X .LT. A(IA)) GO TO 10 
    ELSE 
20    IF (X .GT. A(IA+1)) THEN 
        IA=IA+1                        ! look forwards
        GO TO 20 
      END IF 
    END IF 
  END IF

  IF (IA.GE.2 .AND. IA.LE.M) L=4 
  I=LIMIT(1,IA-1,M) 
  I1(K)=I 
  I2(K)=I+L-1 
  NN(K)=L 

  IF (L.GT.1) THEN      !     PROCEDURE (CALCULATE AIJ VALUES)                                  
    IF (L.GT.2) THEN 
      IF (L.GT.3) THEN 
        A44=A(I+3) 
        A14=A(I)-A44 
        A24=A(I+1)-A44 
        A34=A(I+2)-A44 
        A44=X-A44 
      END IF 
      A33=A(I+2) 
      A13=A(I)-A33 
      A23=A(I+1)-A33 
      A33=X-A33 
    END IF 
    A22=A(I+1) 
    A12=A(I)-A22 
    A22=X-A22 
  END IF 
  A11=X-A(I) 


  SELECT CASE(L)             !     PROCEDURE (CALCULATE COEFFICIENTS)
    CASE(1)
      C(1,K)=1.0 
    CASE(2)
      C(1,K)=+A22/A12 
      C(2,K)=-A11/A12 
    CASE(3)
      C(1,K)=+A22/A12*A33/A13 
      C(2,K)=-A11/A12*A33/A23 
      C(3,K)=+A11/A13*A22/A23 
    CASE(4)
      P1=A22/A23*A33 
      C(1,K)=+P1/A12*A33/A13 
      C(4,K)=-P1/A34*A22/A24 
      P2=A33/A23*A11/A23 
      P3=A22/A23*A44/A23 
      C(2,K)=-A33*(P2/A12+P3/A24) 
      C(3,K)=+A22*(P3/A34+P2/A13) 
  END SELECT
  RETURN 
                                                                        

CONTAINS
FUNCTION Limit(imin, i, imax) RESULT(k)
  INTEGER,INTENT(IN):: imin,i,imax
  INTEGER:: k
  k=MAX(IMIN,MIN(I,IMAX))
  RETURN
END Function Limit   ! ======================================================
END Subroutine Ntrp   ! =====================================================
