!+
! PROGRAM MassProperties
! ------------------------------------------------------------------------------
! PURPOSE - Calculate the mass properties of complex rigid structural systems.

! AUTHORS - Reid A. Hull, John L. Gilbert, and Phillip J. Klitch, NASA Langley
!         - R. C. WARD / CSC
!         - B.S. Garbow, Argonne National Laboratory
!         - Bowdler, Martin, Reinsch, and Wilkinson, 
!             Handbook of Automatic Computation
!         - Ralph L. Carmichael, Public Domain Aeronautical Software

! REVISION HISTORY
!   DATE  VERS PERSON  STATEMENT OF CHANGES
!   1975   1.0    "    Publication of NASA Technical Memorandum (March 1978)
!   1996   1.1   RLC   Acquisition of COSMIC Program LAR-11629
!   2005   1.2   RLC   Converted all source to free-format Fortran 90
! 29Aug08  1.3   RLC   Put procedures in modules
! 22Dec08  1.4   RLC   Changed input and output file numbers
! 10Feb09  1.5   RLC   Added WP to allow easy change from single to double precision


! The original program description from COSMIC...
! The computer program MASPROP was developed to rapidly calculate the mass 
! properties of complex rigid structural systems. This program's basic premise
! is that complex systems can be adequately described by a combination of basic 
! elementary structural shapes. Thirteen widely used basic structural shapes 
! are available in this program. They are as follows: Discrete Mass, Cylinder, 
! Truncated Cone, Torus, Beam (arbitrary cross section), Circular Rod 
! (arbitrary cross section), Spherical Segment, Sphere, Hemisphere, 
! Parallelepiped, Swept Trapezoidal Panel, Symmetric Trapezoidal Panels, and a 
! Curved Rectangular Panel. MASPROP provides a designer with a simple technique 
! that requires minimal input to calculate the mass properties of a complex 
! rigid structure and should be useful in any situation where one needs to 
! calculate the center of gravity and moments of inertia of a complex structure.

! Rigid body analysis is used to calculate mass properties. Mass properties are 
! calculated about component axes that have been rotated to be parallel to the 
! system coordinate axes. Then the system center of gravity is calculated and 
! the mass properties are transferred to axes through the system center of 
! gravity by using the parallel axis theorem. System weight, moments of inertia 
! about the system origin, and the products of inertia about the system center 
! of mass are calculated and printed. From the information about the system 
! center of mass the principal axes of the system and the moments of inertia 
! about them are calculated and printed.

! The only input required is simple geometric data describing the size and 
! location of each element and the respective material density or weight of 
! each element.





!+
MODULE MassPropertiesProcedures
! ------------------------------------------------------------------------------

IMPLICIT NONE

  CHARACTER(LEN=*),PARAMETER:: VERSION = "1.5 (10 February 2009)"
  INTEGER,PARAMETER,PRIVATE:: DBG = 3
  INTEGER,PARAMETER:: SP = SELECTED_REAL_KIND(6)    ! single precision
  INTEGER,PARAMETER:: DP = SELECTED_REAL_KIND(14)   ! double precision
  INTEGER,PARAMETER:: WP = DP   ! working precision

CONTAINS

!+
SUBROUTINE ReadAndScanInputData(efu1,efu2, nitems)
! ------------------------------------------------------------------------------
! PURPOSE - Read the records in the input file. Each object should have three
!  records. The first is a description of the object; the second has seven
! numbers; and the third has nine numbers. None are changed by this routine. 
! Returns the number of valid objects.

  INTEGER,INTENT(IN):: efu1   ! unit number of the input file  
  INTEGER,INTENT(IN):: efu2   ! unit number of the output file
  INTEGER,INTENT(OUT):: nItems   ! number of valid items found on input file

  REAL(WP):: a,b,c,d,f
  CHARACTER(LEN=80):: description, buffer
  INTEGER:: errCode
  INTEGER:: kgood=0
  REAL(WP):: rho
  INTEGER:: shapeCode
  REAL(WP):: xi,yi,zi, xj,yj,zj, xk,yk,zk
!-------------------------------------------------------------------------------
  WRITE(efu2,*) 'INPUT DATA'
  DO
    READ(efu1,'(A)',IOSTAT=errCode) description
    IF (errCode /= 0) EXIT
    WRITE(efu2,*) TRIM(description)

    READ(efu1,'(A)',IOSTAT=errCode) buffer
    IF (errCode /= 0) EXIT
    WRITE(efu2,*) TRIM(buffer)
    READ(buffer,*,IOSTAT=errCode) shapeCode, rho,a,b,c,d,f
    IF (errCode /= 0) EXIT


    READ(efu1,'(A)',IOSTAT=errCode) buffer
    IF (errCode /= 0) EXIT
    WRITE(efu2,*) TRIM(buffer)
    READ(buffer,*,IOSTAT=errCode) xi,yi,zi, xj,yj,zj, xk,yk,zk
    IF (errCode /= 0) EXIT
    
    kgood=kgood+1   ! has to read all three successfully to get here
    
  END DO

  WRITE(*,*)    kgood, " objects read from input file"
  WRITE(efu2,*) kgood, " objects read from input file"
  nitems=kgood
  RETURN
END Subroutine ReadAndScanInputData   ! ----------------------------------------

!+
SUBROUTINE SymQL(a, e, wk, ierr)
! ------------------------------------------------------------------------------
! PURPOSE: computes the eigenvalues and eigenvectors of a symmetric matrix
! NOTES - on input, a contains a symmetric matrix (only full lower triangle 
!  need be supplied. On output, a contains the contains orthonormal eigenvectors
! NOTES - Eigenvalues are stored in ascending order.
!   Upon error exit, eigenvalues are correct but unordered for indices 1,2,...,ierr-1.
!   The vector associated with the i-th eigenvalue is found in the i-th column of a.
!   Upon error exit, a contains eigenvectors associated with the stored eigenvalues.

!   REQUIRED ROUTINES    - TRED2,TQL2
!   AUTHOR/IMPLEMENTOR   - R. C. WARD / CSC
!   DATE RELEASED        - OCT. 19, 1972
!   LATEST REVISION      - DECEMBER 15, 1976

  REAL(WP),INTENT(IN OUT),DIMENSION(:,:):: a   
  REAL(WP),INTENT(OUT),DIMENSION(:):: e   ! contains eigenvalues in ascending order
  REAL(WP),INTENT(IN OUT),DIMENSION(:):: wk   ! working storage of dimension n
  INTEGER,INTENT(OUT):: ierr   ! error code
                               ! = 0   normal return
                               ! = j   j-th eigenvalue has not been determined 
                               !            after 30 iterations

  INTEGER:: n
!-------------------------------------------------------------------------------
  n=SIZE(a,1)
  write(dbg,*) 'Entering SymQl'
  write(dbg,*) a(1,:)
  write(dbg,*) a(2,:)
  write(dbg,*) a(3,:)
  CALL Tred2(a(1:n,1:n), e, wk, a(1:n,1:n))
  CALL Tql2(e, wk, a(1:n,1:n), ierr)
  RETURN
END Subroutine Symql   ! -------------------------------------------------------

!+
SUBROUTINE Tred2(a, d, e, z)
! ------------------------------------------------------------------------------
! PURPOSE: Reduce a real symmetric matrix to a symmetric tridiagonal 
!  matrix (accumulating transformations)
IMPLICIT NONE
!!!  INTEGER,INTENT(IN):: NM    !  - MAXIMUM ROW DIMENSION OF A
!!!  INTEGER,INTENT(IN):: N     !  - ORDER OF A
  REAL(WP),INTENT(IN OUT),DIMENSION(:,:):: a   ! symmetric matrix 
                       !   (only full lower triangle of matrix need be supplied)
  REAL(WP),INTENT(OUT),DIMENSION(:):: d   ! output array containing diagonal 
                                      !   elements of the tridiagonal matrix
  REAL(WP),INTENT(OUT),DIMENSION(:):: e   ! output array containing subdiagonal 
    !   elements of tridiagonal matrix in its last n-1 positions. e(1) = 0
  REAL(WP),INTENT(OUT),DIMENSION(:,:):: z   ! contains the orthogonal transformation 
    ! matrix produced in the reduction  (may occupy same locations as a)

!   REQUIRED ROUTINES   - NONE
!   AUTHOR/IMPLEMENTER   - R.C. WARD / R.C. WARD
!   LANGUAGE            - FORTRAN
!   DATE RELEASED       - OCT. 19, 1972
!   LATEST REVISION     - FEB. 28, 1973

  INTEGER:: i,j,k,l   !!! ii,jp1  (never referenced)
  INTEGER:: n 
  REAL(WP):: f,g,h,hh,scale
!----------------------------------------------------------------------------
!!!  f=0.0   ! added by RLC
  n=SIZE(a,1)
  write(dbg,*) 'Entering Tred2'
  write(dbg,*) a(1,:)
  write(dbg,*) a(2,:)
  write(dbg,*) a(3,:)

  DO i = 1, n
    DO j = 1, i
      z(i,j) = a(i,j)
    END DO
  END DO
  write(dbg,*) 'z'
  write(dbg,*) z(1,:)
  write(dbg,*) z(2,:)
  write(dbg,*) z(3,:)
   
  IF (n .EQ. 1) GO TO 320


  DO i=n,2,-1
!     ********** FOR I=N STEP -1 UNTIL 2 DO -- **********
!!!      DO  II = 2, N   ! big loop
!!!         I = N + 2 - II
    L = i - 2
    h = 0.0


    IF (L .LT. 1) GO TO 130

    scale = 0.0   !     ********** SCALE ROW (ALGOL TOL THEN NOT NEEDED) 
    DO k = 1, l
      scale = scale + ABS(z(i,k))
    END DO

    IF (scale .NE. 0.0) GO TO 140

130 e(i) = z(i,L+1)   ! was 130
    GO TO 290

  140    L = L + 1
    scale = scale + ABS(z(i,l))


    DO K = 1, L
      Z(I,K) = Z(I,K) / SCALE
      H = H + Z(I,K) * Z(I,K)
    END DO

    F = Z(I,L)
    G = -SIGN(SQRT(H),F)
    E(I) = SCALE * G
    H = H - F * G
    Z(I,L) = F - G
    F = 0.0

    DO J = 1, L
      Z(J,I) = Z(I,J) / (SCALE * H)
      G = 0.0   !     ********** FORM ELEMENT OF A*U **********
      DO K = 1, J
        G = G + Z(J,K) * Z(I,K)
      END DO
!!!            JP1 = J + 1
!!!            IF (L .LT. JP1) GO TO 220

      DO K = J+1, L
        G = G + Z(K,J)*Z(I,K)
      END DO

      E(J) = G / H   !  was 220   ********** FORM ELEMENT OF P **********
  write(dbg,*) 'line 229', j,e 
     F = F + E(J) * Z(I,J)

    END DO

    HH = F / (H + H)
!     ********** FORM REDUCED A **********
    DO J = 1, L
      F = Z(I,J)
      G = E(J) - HH * F
      E(J) = G
      DO K = 1, J
        Z(J,K) = Z(J,K) - F * E(K) - G * Z(I,K)
      END DO
    END DO

    DO K = 1, L
      Z(I,K) = SCALE * Z(I,K)
    END DO

  290    D(I) = H
  END DO   ! end of big loop

320 D(1) = 0.0
  E(1) = 0.0
!     ********** ACCUMULATION OF TRANSFORMATION MATRICES **********
  DO I = 1, N
    L = I - 1
    IF (D(I) .NE. 0.0) THEN
      DO  J = 1, L
        G = 0.0
        DO  K = 1, L
          G = G + Z(I,K) * Z(K,J)
        END DO
        DO K = 1, L
          Z(K,J) = Z(K,J) - G * Z(K,I)
        END DO
      END DO
    END IF

    d(i) = z(i,i)
    z(i,i) = 1.0
    DO J = 1, l
      z(i,j) = 0.0
      z(j,i) = 0.0
    END DO
  END DO

  RETURN
END Subroutine Tred2   ! -------------------------------------------------------


!+
SUBROUTINE TQL2(D,E,Z,IERR)
! PURPOSE - Find the eigenvalues and eigenvectors of a symmetric tridiagonal 
!  matrix by the QL method. The eigenvectors of a full symmetric matrix can 
!  also be found if  Tred2  has been used to reduce this full matrix to 
!  tridiagonal form.
!  This subroutine is a translation of the algol procedure TQL2,
!  Num. Math. 11, 293-306(1968) by Bowdler, Martin, Reinsch, and Wilkinson.
!  Handbook for Auto. Comp., vol.II-Linear Algebra, 227-240(1971).
IMPLICIT NONE
  REAL(WP),INTENT(IN OUT),DIMENSION(:):: d   ! On input, d contains the diagonal
                                             ! elements of the input matrix.
                                             ! On output, d contains the 
                                             ! eigenvalues in ascending order. 
                                             ! If an error exit is made, the 
                                             ! eigenvalues are correct but
                                             ! unordered for indices 1,2,...,ierr-1,

  REAL(WP),INTENT(IN OUT),DIMENSION(:):: e   ! On input, e contains the subdiagonal 
                                             ! elements of the input matrix in its 
                                             ! last n-1 positions;  e(1) is arbitrary.
                                             ! On output, e has been destroyed.

  REAL(WP),INTENT(IN OUT),DIMENSION(:,:):: z   ! On input, z contains the 
                                               ! transformation matrix produced in the
                                               ! reduction by Tred2, if performed.  
                                               ! If the eigenvectors of the tridiagonal 
                                               ! matrix are desired, z must contain
                                               ! the identity matrix.
                                               ! On output, z contains orthonormal 
                                               ! eigenvectors of the symmetric 
                                               ! tridiagonal (or full) matrix.  
                                               ! If an error exit is made, z contains 
                                               ! the eigenvectors associated with the 
                                               ! stored eigenvalues.

  INTEGER,INTENT(OUT):: ierr   !   =0 for normal return
                               !   =j if the j-th eigenvalue has not been determined 
                               !        after 30 iterations.

  INTEGER:: I,J,K,L,M,II,L1,MML
  INTEGER:: n
!      REAL D(N),E(N),Z(NM,N)
  REAL(WP):: B,C,F,G,H,P,R,S,MACHEP

  INTRINSIC:: SQRT,ABS,SIGN
!

!     QUESTIONS AND COMMENTS SHOULD BE DIRECTED TO B. S. GARBOW,
!     APPLIED MATHEMATICS DIVISION, ARGONNE NATIONAL LABORATORY
!
!     ********** MACHEP IS A MACHINE DEPENDENT PARAMETER SPECIFYING
!                THE RELATIVE PRECISION OF FLOATING POINT ARITHMETIC.
!

      MACHEP = 2.**(-47)
!----------------------------------------------------------------------------
      IERR = 0
  n=SIZE(z,1)

      IF (N .EQ. 1) GO TO 1001

      DO  I = 2, N
        E(I-1) = E(I)
      END DO

      F = 0.0
      B = 0.0
      E(N) = 0.0

      DO L = 1, N   ! start of big loop  (was do 240)
         J = 0
         H = MACHEP * (ABS(D(L)) + ABS(E(L)))
         IF (B .LT. H) B = H

!     ********** LOOK FOR SMALL SUB-DIAGONAL ELEMENT **********
!     ********** E(N) IS ALWAYS ZERO, SO THERE IS NO EXIT
!                THROUGH THE BOTTOM OF THE LOOP **********
         DO  M = L, N
            IF (ABS(E(M)) .LE. B) GO TO 120
         END DO
!
  120    IF (M .EQ. L) GO TO 220
  130    IF (J .EQ. 30) GO TO 1000
         J = J + 1
!     ********** FORM SHIFT **********
         L1 = L + 1
         G = D(L)
         P = (D(L1) - G) / (2.0 * E(L))
         R = SQRT(P*P+1.0)
         D(L) = E(L) / (P + SIGN(R,P))
         H = G - D(L)
!
         DO  I = L1, N
           D(I) = D(I) - H
         END DO
!
         F = F + H
!     ********** QL TRANSFORMATION **********
         P = D(M)
         C = 1.0
         S = 0.0
         MML = M - L
!     ********** FOR I=M-1 STEP -1 UNTIL L DO -- **********
         DO II = 1, MML   ! was do 200
            I = M - II
            G = C * E(I)
            H = C * P
            IF (ABS(P) .LT. ABS(E(I))) GO TO 150
            C = E(I) / P
            R = SQRT(C*C+1.0)
            E(I+1) = S * P * R
            S = C / R
            C = 1.0 / R
            GO TO 160
  150       C = P / E(I)
            R = SQRT(C*C+1.0)
            E(I+1) = S * E(I) * R
            S = 1.0 / R
            C = C * S
  160       P = C * D(I) - S * G
            D(I+1) = H + S * (C * G + S * D(I))
!     ********** FORM VECTOR **********
            DO K = 1, N
               H = Z(K,I+1)
               Z(K,I+1) = S * Z(K,I) + C * H
               Z(K,I) = C * Z(K,I) - S * H
            END DO
!
         END DO   ! was 200
!
         E(L) = S * P
         D(L) = C * P
         IF (ABS(E(L)) .GT. B) GO TO 130
  220    D(L) = D(L) + F
      END DO   ! was 240




!     ********** ORDER EIGENVALUES AND EIGENVECTORS **********
      DO II = 2, N
         I = II - 1
         K = I
         P = D(I)
!
         DO J = II, N
            IF (D(J) .GE. P) CYCLE
            K = J
            P = D(J)
         END DO
!
         IF (K .EQ. I) CYCLE
         D(K) = D(I)
         D(I) = P
!
         DO J = 1, N
            P = Z(J,I)
            Z(J,I) = Z(J,K)
            Z(J,K) = P
         END DO
!
      END DO

      GO TO 1001
!     ********** SET ERROR -- NO CONVERGENCE TO AN
!                EIGENVALUE AFTER 30 ITERATIONS **********
 1000 IERR = L
 1001 RETURN
END Subroutine Tql2   ! -----------------------------------------------------

END Module MassPropertiesProcedures   ! =====================================

!+
PROGRAM MassProperties
! ---------------------------------------------------------------------------
USE MassPropertiesProcedures
IMPLICIT NONE

!     LATEST DIRCOS IN MAIN DECK AND 13 SHAPES 12/20/72
  INTEGER,PARAMETER:: MAXITEMS = 200
  CHARACTER(LEN=80),DIMENSION(MAXITEMS):: description
  INTEGER,DIMENSION(MAXITEMS):: shapeCode
  REAL(WP),DIMENSION(MAXITEMS):: rho, a,b,c,d,f
  REAL(WP),DIMENSION(MAXITEMS):: xi,yi,zi, xj,yj,zj, xk,yk,zk
  REAL(WP),DIMENSION(MAXITEMS):: ixxcg,iyycg,izzcg, ixycg,iyzcg,ixzcg
  REAL(WP),DIMENSION(MAXITEMS):: xcg,ycg,zcg, xl,yl,zl
  REAL(WP),DIMENSION(MAXITEMS):: ixxco,iyyco,izzco, ixyco,ixzco,iyzco
  REAL(WP),DIMENSION(MAXITEMS):: w     !!! ,o,p

  REAL(WP),DIMENSION(3,3):: arr
  REAL(WP),DIMENSION(3):: crr,e

!      DIMENSION  ITEM(200),SHAPE(200),RHO(200),A(200),B(200),C(200),    &
!     & D(200),F(200),       XI(200),YI(200),ZI(200),XJ(200),YJ(200),ZJ(2&
!     &00),XK(200),YK(200),ZK(200),IXXCG(200),IYYCG(200),IZZCG(200),IXYCG&
!     &(200),IYZCG(200),IXZCG(200),XCG(200),YCG(200),ZCG(200),XL(200),YL(&
!     &200),ZL(200),IXXCO(200),IYYCO(200),IZZCO(200),IXYCO(200),IXZCO(200&
!     &),IYZCO(200), W(200),O(200),P(200),ARR(3,3),E(3),CRR(3)

!!!      DIMENSION DES(3,200)
!!!      CHARACTER(LEN=18),DIMENSION(200):: des
  

!      REAL IXX,IYY,IZZ,IXXCG,IYYCG,IZZCG,IXYCG,IXZCG,                   &
!     &IYZCG,LX,MX,NX,LY,MY,NY,LZ,MZ,NZ,IXY,IXZ,IYZ,LGTH,IXXCO,IYYCO,    &
!     &IZZCO,IXYCO,IXZCO,IYZCO,IXXO,IYYO,IZZO,KI,KII,KIII,IXYK,IXYCGP
  REAL(WP):: ixx,ixy,ixz, iyy,iyz, izz
  REAL(WP):: ixycgp
  REAL(WP):: ixyk
  REAL(WP):: ixxo,iyyo,izzo
  REAL(WP):: lx,mx,nx, ly,my,ny, lz,mz,nz
  REAL(WP):: lgth

  REAL(WP):: actemp,acvol
  REAL(WP):: aftan
  REAL(WP):: apxcg
  REAL(WP):: bcg
  REAL(WP):: cgk
  REAL(WP),PARAMETER:: CONS=4632.0 ! = 32.16 * 12 * 12
  INTEGER,PARAMETER:: DBG = 3
  REAL(WP):: delx,dely,delz
  INTEGER:: errCode
  REAL(WP):: fetan
  CHARACTER(LEN=80):: fileName
  CHARACTER(LEN=*),PARAMETER:: FMT1 = '(" eigenvalue  ", I3, ES15.5)'
  CHARACTER(LEN=*),PARAMETER:: FMT2 = '(" eigenvector ", I3, 3ES15.5)'
  CHARACTER(LEN=*),PARAMETER:: FMTD = '(A,4ES13.5)'
  REAL(WP):: g,h
  INTEGER:: ierr
  INTEGER:: i,j
  INTEGER,PARAMETER:: IN = 1
  INTEGER:: ki,kii,kiii
!!!  INTEGER:: n    ! never referenced
  INTEGER:: nGood 
  INTEGER,PARAMETER:: OUT=2
  REAL(WP):: paxrc
  REAL(WP),PARAMETER:: PI = 3.1415265
  REAL(WP):: prb
  REAL(WP):: prixyk
  REAL(WP):: rscg
  REAL(WP):: shfti
  REAL(WP):: smh
  REAL(WP):: t1,t2,t3,t4, temp1,temp2,temp3,temp4
  REAL(WP):: thrb
  REAL(WP):: tm
  REAL(WP):: tmpixx
  REAL(WP):: tntau
  REAL(WP):: tw
  REAL(WP):: vol,vol1,vol2, volcyl
  REAL(WP):: xbar,ybar,zbar, xbar1,xbar2
  REAL(WP):: xfrxx
  REAL(WP):: xi1,xi2,xi3,xi4,xi5,xi6
  REAL(WP):: xlt
  REAL(WP):: xm1,xm2
  REAL(WP):: xmom,ymom,zmom
  REAL(WP):: xpan
  REAL(WP):: xt
  REAL(WP):: xtlgth
  REAL(WP):: yi1,yi2
!----------------------------------------------------------------------------
  WRITE(*,*) 'Mass properties of a rigid structure'
  WRITE(*,*) 'Version '//VERSION
!!!  WRITE(*,*) 'WP', WP
  DO
    WRITE(*,*) 'Enter the name of the input file: '
    READ(*,'(A)') fileName
    IF (LEN_TRIM(fileName)==0) STOP
    OPEN(UNIT=IN,FILE=fileName,STATUS='OLD',IOSTAT=errCode,ACTION='READ')
    IF (errCode==0) EXIT
    OPEN(UNIT=IN,FILE=Trim(fileName)//'.inp',STATUS='OLD',IOSTAT=errCode,ACTION='READ')
    IF (errCode==0) EXIT
    WRITE(*,*) 'Unable to open this file. Try again.'
  END DO
  OPEN(UNIT=OUT,FILE='massprop.out',STATUS='REPLACE',ACTION='WRITE')
  OPEN(UNIT=DBG,FILE='massprop.dbg',STATUS='REPLACE',ACTION='WRITE')

! First, use ReadAndScanInputData to see how many good items are on the input file.
! Then, read them into memory. No need to do error checking here, as it was
! done in ReadAndScanInputData.
  CALL ReadAndScanInputData(IN,OUT,nGood)
  REWIND(UNIT=IN)
  DO i=1,nGood  
    READ(IN,'(A)') description(i)
    READ(IN,*) shapeCode(i),rho(i),a(i),b(i),c(i),d(i),f(i)
    READ(IN,*) xi(i),yi(i),zi(i), xj(i),yj(i),zj(i), xk(i),yk(i),zk(i)
  END DO

!!!      IF(EOF(5))1009,1010
!!! 1009 IMAX=I-1
!!!      WRITE(6,271)
!!!  271 FORMAT(1H0' INPUT DATA LISTED BELOW'//)
!!!      WRITE(6,103)(ITEM(I),DES(I),SHAPE(I),RHO(I),A(I),B(I), &
!!!     &C(I),D(I),F(I),     XI(I),YI(I),ZI(I),XJ(I),YJ(I),ZJ(I),XK(I),    &
!!!     &YK(I),ZK(I),I=1,IMAX)
!!!  103 FORMAT (2X'ITEM       DESCRIPTION       SHAPE       RHO           &
!!!     &     A              B              C              D              F&
!!!     &  '/I5,6X,2A9,I5,1X,E16.8,5F15.3/'        XI            YI        &
!!!     &    ZI            XJ            YJ            ZJ            XK    &
!!!     &        YK            ZK'/9F14.3//)


  WRITE(OUT,*) '   COMPONENT  DATA  LISTED  BELOW'
!!!  WRITE(OUT,*) ' '
!!!  WRITE(OUT,*) ' ITEM       DESCRIPTION          WT         IXXCO      '// &
!!!     ' IYYCO     IZZCO     XCGCO      YCGCO    ZCGCO    IXYCO     IYZCO   IXZCO'
!!!  WRITE(OUT,*) ' '

  tw=0.
  xmom=0.
  ymom=0.
  zmom=0.
  ixx=0.
  iyy=0.
  izz=0.
  ixy=0.
  ixz=0.
  iyz=0.
  ixxo=0.
  iyyo=0.
  izzo=0.
  
  DO  i=1,nGood
    xl(i)=0.
    yl(i)=0.
    zl(i)=0.
    ixycg(i)=0.
    ixzcg(i)=0.
    iyzcg(i)=0.
    lgth=SQRT((xj(i)-xi(i))**2 + (yj(i)-yi(i))**2 + (zj(i)-zi(i))**2)
!    itemp=shape(i)

    WRITE(DBG,*) ' '
    WRITE(DBG,*) 'Evaluating element with shape code ', shapeCode(i)
    WRITE(DBG,*) TRIM(description(i))
    WRITE(DBG,FMTD) 'Values used: a,b,c,rho=',a(i),b(i),c(i),rho(i)
    SELECT CASE(shapeCode(i))
!!!      GO TO (1,2,3,4,5,6,7,8,9,10,11,12,13),ITEMP
      CASE(1)   ! discrete mass   *shape 1*
                !     point  j must be selected anywhere on the x axis for dir. cosines
        W(I)=RHO(I)
        ixxcg(i)=a(i)*CONS
        iyycg(i)=b(i)*CONS
        izzcg(i)=c(i)*CONS
        IF((ABS(xj(i)).NE. 0.).OR.(ABS(yj(i)).NE.0.).OR. (ABS(zj(i)).NE.0.0)) GO TO 704
        lgth=10.
        xj(i)=1.1*xi(i)+10.
        IF (xj(i).EQ.xi(i))xj(i)=xi(i)+10.
  704 xl(i)=0.0
!      WRITE(DBG,*) 
      
      CASE(2)   !    CYLINDER   *SHAPE 2*
        IF (B(I) .EQ. A(I)) THEN
          W(I)=RHO(I)   !     CYLINDER (THIN WALL)   *SHAPE 2*
          XL(I)=LGTH/2.
          IXXCG(I)=RHO(I)*A(I)**2
          IYYCG(I)=0.5*RHO(I)*A(I)**2+RHO(I)*(LGTH**3-C(I)**3)/              &
         &((LGTH-C(I))*12.0)
          IZZCG(I)=IYYCG(I)
        ELSE
          IF (B(I).LT.0.) B(I)=0.
          IF (C(I).LT.0.) C(I)=0.
          VOLCYL = PI*(A(I)**2-B(I)**2)*(LGTH-C(I))
          W(I)=RHO(I)*VOLCYL
          IF(RHO(I).GT.0.4) W(I)=RHO(I)
          rho(i)=w(i)/volcyl
          ixxcg(i)=0.5*w(i)*(a(i)**2+b(i)**2)
          iyycg(i)=.25*w(i)*(a(i)**2+b(i)**2) + &
                                      w(i)*(lgth**3-c(i)**3)/(12.*(lgth-c(i)))
          izzcg(i)=iyycg(i)
          xl(i)=lgth/2.
          WRITE(DBG,*) 'volcyl', volcyl
        END IF

      CASE(3)   !     TRUNCATED CONE   *SHAPE 3* (REF--R. HULL)
        IF(B(I) .EQ. A(I)) THEN
          VOL=F(I)*SQRT((A(I)-C(I))**2 + LGTH**2)*3.14159*(A(I)+C(I)) !     TRUNCATED CONE (THIN WALL)   *SHAPE 3*
          W(I) = VOL*RHO(I)
          XL(I)=LGTH/3.*((2.*C(I)+A(I))/(C(I)+A(I)))
          IYYCG(I)=(RHO(I)/4.*(A(I)**2+C(I)**2)+RHO(I)*LGTH**2/18.*(1.+2.*A(&
         &I)*C(I)/(A(I)+C(I))**2))
          IZZCG(I)=IYYCG(I)
          IXXCG(I)=RHO(I)/2.*(A(I)**2+C(I)**2)
        ELSE
          VOL = 1.0472*LGTH*(A(I)**2+A(I)*C(I)+C(I)**2-B(I)**2-B(I)*D(I)-D(I)**2)
          W(I) = RHO(I)*VOL
          IF(RHO(I).GT. 0.4) W(I)=RHO(I)
          RHO(I) = W(I)/VOL
          XL(I)= (LGTH**2)*(A(I)**2+2.*A(I)*C(I)+3.*C(I)**2-B(I)**2-2.*B(I)*&
         &D(I)-3.*D(I)**2)*.2618/VOL
          KI =.05*3.141593*RHO(I)*((C(I)**5-A(I)**5)*LGTH/(C(I)-A(I))-      &
         &(D(I)**5-B(I)**5)*LGTH/(D(I)-B(I)))
          KII=3.1416*RHO(I)*LGTH**3*((A(I)**2-B(I)**2)/3.0-.5*((A(I)-       &
         &C(I))*A(I)-(B(I)-D(I))*B(I))+.2*((A(I)-C(I))**2-(B(I)-D(I))**2))
          KIII=1.0472*(LGTH)*(A(I)**2+A(I)*C(I)+C(I)**2-B(I)**2-B(I)*D(I)-  &
         &D(I)**2)*RHO(I)*XL(I)**2
          IXXCG(I)=2.0*KI
          IYYCG(I)=KI+KII-KIII
          IZZCG(I)=IYYCG(I)
          WRITE(DBG,*) 'ki,kii,kiii=', ki,kii,kiii
          WRITE(DBG,*) 'vol', vol
        END IF

      CASE(4)   !     TORUS   *SHAPE 4*
        A(I)=A(I)-LGTH
        IF(D(I).LT.0.)D(I)=0.
        VOL1=2.*PI**2*LGTH**2*A(I)
        VOL2=2.*PI**2*D(I)**2*A(I)
        ACVOL=VOL1-VOL2
        W(I)=RHO(I)*ACVOL
        XM1=(RHO(I)*VOL1)
        XI1=.125*XM1*(4.*A(I)**2+5.*LGTH**2)
        XM2=(RHO(I)*VOL2)
        XI2=.125*XM2*(4.*A(I)**2+5.*D(I)**2)
        IYYCG(I)=(XI1-XI2)
        IZZCG(I)=IYYCG(I)
        YI1=.25*XM1*(4.*A(I)**2+3.*LGTH**2)
        YI2=.25*XM2*(4.*A(I)**2+3.*D(I)**2)
        IXXCG(I)=YI1-YI2
        WRITE(DBG,FMTD) 'Input values used: d', d(i)
        WRITE(DBG,FMTD) 'vol1,vol2,acvol', vol1,vol2,acvol
        WRITE(DBG,FMTD) 'xm1,xm2,xi1,xi2', xm1,xm2,xi1,xi2
      
      CASE(5)   !     BEAM (ARBITRARY CROSS SECTION)   *SHAPE 5*
        XL(I)=LGTH/2.
        IXXCG(I)=(RHO(I)*LGTH*(B(I)+C(I)))
        IYYCG(I)=(RHO(I)*(B(I)*LGTH+.0833*A(I)*LGTH**3))
        IZZCG(I)=(RHO(I)*(C(I)*LGTH+.0833*A(I)*LGTH**3))
        VOL=A(I)*LGTH
        W(I)=RHO(I)*VOL
        WRITE(DBG,FMTD) 'vol', vol

      CASE(6)   !     CIRCULAR ROD (ARBITRARY CROSS SECTION)   *SHAPE 6*
        W(I)=RHO(I)
        RSCG=(A(I)**2*C(I)+B(I))/(A(I)*C(I))
        IXXCG(I)=RSCG**2*W(I)
        IYYCG(I)=.5*RSCG**2*W(I)
        IZZCG(I)=IYYCG(I)
        WRITE(DBG,FMTD) 'rscg', rscg
    
      CASE(7)   !     SPHERICAL SEGMENT   *SHAPE 7*
        F(I)=C(I)-B(I)
        D(I)=C(I)-LGTH
        G   =LGTH-B(I)
        VOL1=1.0472*LGTH**2*(3.*C(I)-LGTH)
        VOL2=1.0472*G   **2*(3.*F(I)-G   )
        ACVOL=VOL1-VOL2
        W(I)=ACVOL*RHO(I)
        IF(RHO(I).GT.0.4) W(I)=RHO(I)
        RHO(I) = W(I)/ACVOL
        XBAR1=.75*(2.*C(I)-LGTH)**2/(3.*C(I)-LGTH)
        XBAR2=.75*(2.*F(I)-G   )**2/(3.*F(I)-G   )
        XLT  =(XBAR1*VOL1-XBAR2*VOL2)/ACVOL
        XL(I)=XLT  -D(I)
        XM1=VOL1*RHO(I)
        XM2=VOL2*RHO(I)
        TM=XM1-XM2
        XI1=(2.*LGTH*XM1/(3.*C(I)-LGTH))*(C(I)**2-.75*C(I)*LGTH+.15*LGTH**2)
        XI2=(2.*G   *XM2/(3.*F(I)-G   ))*(F(I)**2-.75*F(I)*G   +.15*G**2)
        IXXCG(I)=(XI1-XI2)
        TEMP1=.05236*RHO(I)*(15.*C(I)**4*C(I)-10.*C(I)**2*C(I)**3+3.*C(I) &
       &**5)+.20944*RHO(I)*(5.*C(I)**2*C(I)**3-3.*C(I)**5)
        TEMP2=.05236*RHO(I)*(15.*C(I)**4*D(I)-10.*C(I)**2*D(I)**3+3.*D(I) &
       &**5)+.20944*RHO(I)*(5.*C(I)**2*D(I)**3-3.*D(I)**5)
        TEMP3=.05236*RHO(I)*(15.*F(I)**4*F(I)-10.*F(I)**2*F(I)**3+3.*F(I) &
       &**5)+.20944*RHO(I)*(5.*F(I)**2*F(I)**3-3.*F(I)**5)
        TEMP4=.05236*RHO(I)*(15.*F(I)**4*D(I)-10.*F(I)**2*D(I)**3+3.*D(I) &
       &**5)+.20944*RHO(I)*(5.*F(I)**2*D(I)**3-3.*D(I)**5)
        ACTEMP=((TEMP1-TEMP2)-(TEMP3-TEMP4))-(TM*XLT  **2)
        IYYCG(I)=ACTEMP
        IZZCG(I)=IYYCG(I)
        WRITE(DBG,FMTD) 'vol1,vol2,acvol', vol1,vol2,acvol
        WRITE(DBG,FMTD) 'xbar1,xbar2,xlt', xbar1,xbar2,xlt
        WRITE(DBG,FMTD) 'xm1,xm2,tm', xm1,xm2,tm
        WRITE(DBG,FMTD) 'xi1,xi2', xi1,xi2
        WRITE(DBG,FMTD) 'temp1,temp2,temp3,temp4', temp1,temp2,temp3,temp4
        WRITE(DBG,FMTD) 'actemp', actemp
    
      CASE(8)   !     SPHERE   *SHAPE 8*
        LGTH=A(I)
        XJ(I)=XI(I)+LGTH
        YJ(I)=YI(I)
        ZJ(I)=ZI(I)
        IF (B(I).LT. 0.) B(I)=0.
        IF (B(I).EQ. A(I)) THEN
          W(I)=RHO(I)   !     SPHERE (THIN WALL)   *SHAPE 8*
          IXXCG(I)=.6667*RHO(I)*A(I)**2
          IYYCG(I)=IXXCG(I)
          IZZCG(I)=IXXCG(I)
          XL(I)=0.0
        ELSE
          VOL1=4.188791*A(I)**3
          VOL2=4.188791*B(I)**3
          ACVOL=VOL1-VOL2
          W(I)=RHO(I)*ACVOL
          XM1=RHO(I)*VOL1
          XM2=RHO(I)*VOL2
          XI1=0.4*XM1*LGTH**2
          XI2=0.4*XM2*B(I)**2
          IXXCG(I)=XI1-XI2
          IYYCG(I)=IXXCG(I)
          IZZCG(I)=IXXCG(I)
          WRITE(DBG,FMTD) 'vol1,vol2,acvol', vol1,vol2,acvol
          WRITE(DBG,FMTD) 'xm1,xm2,xi1,xi2', xm1,xm2, xi1,xi2
        END IF


      CASE(9)   !     HEMISPHERE  *SHAPE 9*
        IF (LGTH.EQ.A(I)) THEN
          W(I)=RHO(I)   !     HEMISPHERE (THIN WALL)   *SHAPE 9*
          XL(I)=LGTH/2.
          IXXCG(I)=.666*RHO(I)*LGTH**2
          IYYCG(I)=.4166*RHO(I)*LGTH**2
          IZZCG(I)=IYYCG(I)
        ELSE
          IF(A(I).LT.0.)A(I)=0.
          VOL1=2.09439*LGTH**3
          VOL2=2.09439*A(I)**3
          ACVOL=VOL1-VOL2
          W(I)=RHO(I)*ACVOL
          XBAR1=.375*LGTH
          XBAR2=.375*A(I)
          XL(I)=(XBAR1*VOL1-XBAR2*VOL2)/ACVOL
          XM1=(RHO(I)*VOL1)
          XM2=(RHO(I)*VOL2)
          XI1=(.4*XM1*LGTH**2)
          XI2=(.4*XM2*A(I)**2)
          IXXCG(I)=(XI1-XI2)
          XI3=(.26*XM1*LGTH**2)
          XI4=(.26*XM2*A(I)**2)
          IYYCG(I)=(XI3-XI4)
          IZZCG(I)=IYYCG(I)
        END IF

      CASE(10)   !     PARALLELEPIPED  *SHAPE 10*

        IF (D(I).EQ.A(I)) THEN
          XL(I)=LGTH/2.   !     PARALLELEPIPED (THIN WALL)  *SHAPE 10*
          W(I)=RHO(I)
          TEMP1=(LGTH*B(I)*A(I))
          TEMP2=(LGTH*B(I)+B(I)*A(I)+LGTH*A(I))
          IXXCG(I)=(.083333*RHO(I)*(B(I)**2+A(I)**2)+(RHO(I)/6.)*(TEMP1*    &
         &(B(I)+A(I))/TEMP2))
          IYYCG(I)=(.083333*RHO(I)*(LGTH**2+A(I)**2)+(RHO(I)/6.)*(TEMP1*    &
         &(LGTH+A(I))/TEMP2))
          IZZCG(I)=(.083333*RHO(I)*(LGTH**2+B(I)**2)+(RHO(I)/6.)*(TEMP1*    &
         &(LGTH+B(I))/TEMP2))
        ELSE
          IF(D(I).LT.0.)D(I)=0.
          VOL1=LGTH*B(I)*A(I)
          VOL2=C(I)*F(I)*D(I)
          ACVOL=VOL1-VOL2
          W(I)=RHO(I)*ACVOL
          IF(RHO(I).GT.0.4) W(I)=RHO(I)
          RHO(I) = W(I)/ACVOL
          XL(I)=LGTH/2.
          XM1=VOL1*RHO(I)
          XM2=VOL2*RHO(I)
          XI1=(.083333*XM1*(B(I)**2+A(I)**2))
          XI2=(.083333*XM2*(F(I)**2+D(I)**2))
          IXXCG(I)=(XI1-XI2)
          XI3=(.083333*XM1*(LGTH**2+A(I)**2))
          XI4=(.083333*XM2*(C(I)**2+D(I)**2))
          IYYCG(I)=(XI3-XI4)
          XI5=(.083333*XM1*(LGTH**2+B(I)**2))
          XI6=(.083333*XM2*(C(I)**2+F(I)**2))
          IZZCG(I)=(XI5-XI6)
        END IF
!!!      WRITE(6,828) XI1,XI2,XI3,XI4,XI5,XI6
!!!  828 FORMAT ('BXI1'6E16.8)
   

      CASE(11)   !     SWEPT TRAPEZOIDAL PANEL (THICK WALL) *SHAPE 11*  (REF--R. HULL)
        W(I)=A(I)*RHO(I)*(LGTH*(B(I)+F(I))/2.)
        IF(RHO(I).GT.0.4) W(I)=RHO(I)
        RHO(I)=(W(I)) /(((B(I)+F(I))/2.0)*LGTH*A(I))
        XL(I)=LGTH*(B(I)+2.0*F(I))/(3.0*(B(I)+F(I)))
        XT = F(I)*LGTH/(B(I)-F(I))
        TNTAU = C(I)/LGTH
        FETAN = (F(I)/2.- B(I)/2.+ C(I))/LGTH
        AFTAN = (C(I)-F(I)/2.0+B(I)/2.0)/LGTH
        THRB=LGTH+XT
        BCG=((FETAN+AFTAN)/2.)*(THRB-XL(I))
        XPAN=AFTAN-FETAN
        IYYCG(I)= RHO(I)*((XPAN*A(I)**3*(THRB**2-XT**2))/24.+  XPAN*A(I)* &
       &(THRB**4-XT**4)/4.)-((THRB-XL(I))**2)*W(I)
        XFRXX=   XPAN*A(I)* (BCG**2*(THRB**2-XT**2)/2.-BCG*(FETAN+AFTAN)* &
       &(THRB**3-XT**3)/3.+(FETAN+AFTAN)**2*(THRB**4-XT**4)/16.)
        IXXCG(I)= RHO(I)* (XPAN**3 *A(I)*(THRB**4-XT**4)/48. +XPAN*A(I)**3&
       &* (THRB**2-XT**2)/24. +XFRXX)
        IZZCG(I)=A(I)*RHO(I)*(XPAN**3*(THRB**4-XT**4)/48.+XPAN*(((THRB-XL(&
       &I))**2+BCG**2)*(THRB**2-XT**2)/2.-(2.*(THRB-XL(I))+(FETAN+AFTAN)*B&
       &CG)*(THRB**3-XT**3)/3.+(4.+(FETAN+AFTAN)**2)*(THRB**4-XT**4)/16.))
        IF (ABS(C(I)).LT. .001) GO TO 60
        APXCG=B(I)/2.+ABS(TNTAU*XL(I))
        H=APXCG+SQRT((XK(I)-XI(I))**2 + (YK(I)-YI(I))**2 + (ZK(I)-ZI(I))**2)
        CGK=H-APXCG
        PRB=THRB-LGTH
        SMH=H*PRB/THRB
        PRIXYK=.5*F(I)*(2.*SMH-F(I)) *(LGTH*PRB/3.+PRB**2/12.)
        IXYK=THRB**2*(H**2-(H-B(I))**2)/24.-PRIXYK
        IXYCGP = ABS( IXYK*A(I)*RHO(I) - W(I)*CGK*XL(I) )*C(I)/ABS(C(I))
        PAXRC=.7854*ABS(C(I))/C(I)
        IF (ABS(IXXCG(I)-IYYCG(I)) .LT. .001) GO TO 20
        PAXRC= ATAN((2.0*IXYCGP)/(IYYCG(I)-IXXCG(I)))/2.0
        IF(IXXCG(I).GT.IYYCG(I)) PAXRC =.5*(ABS(2.*PAXRC)-3.1416)*       &
       & PAXRC/ABS(PAXRC)
     20  TMPIXX = IXXCG(I)*COS(PAXRC)**2 + IYYCG(I)*SIN(PAXRC)**2-           &
       &  2.0*IXYCGP*SIN(PAXRC)*COS(PAXRC)
        IYYCG(I)= IXXCG(I)*SIN(PAXRC)**2 + IYYCG(I)*COS(PAXRC)**2 +        &
       &  2.0*IXYCGP*SIN(PAXRC)*COS(PAXRC) 
        IXXCG(I)=TMPIXX
        SHFTI=ABS(TAN(PAXRC))*XL(I)   ! was labelled 21
        XJ(I)=(XI(I)*(LGTH-XL(I)) + XJ(I)*XL(I))/LGTH
        YJ(I)=(YI(I)*(LGTH-XL(I)) + YJ(I)*XL(I))/LGTH
        ZJ(I)=(ZI(I)*(LGTH-XL(I)) + ZJ(I)*XL(I))/LGTH
        XI(I)=(XI(I)*(SHFTI+CGK) - XK(I)*SHFTI)/CGK
        YI(I)=(YI(I)*(SHFTI+CGK) - YK(I)*SHFTI)/CGK
        ZI(I)=(ZI(I)*(SHFTI+CGK) - ZK(I)*SHFTI)/CGK
        LGTH=SQRT( (XJ(I)-XI(I))**2 + (YJ(I)-YI(I))**2 + (ZJ(I)-ZI(I))**2 )
        XL(I) = LGTH
      
      CASE(12)   !     SYMMETRIC TRAPEZOIDAL PANELS (THICK WALL)*SHAPE 12* (REF--R. HULL)
        LGTH = (LGTH - 2.0* D(I))/2.0
        VOL=(LGTH*(B(I)+F(I))   )*A(I)
        W(I)= RHO(I)*VOL
        IF(RHO(I).GT.0.4) W(I)=RHO(I)
        RHO(I)=(W(I)) /(((B(I)+F(I))    )*LGTH*A(I))
        XL(I)=LGTH*(B(I)+2.0*F(I))/(3.0*(B(I)+F(I)))
        IYYCG(I) = (LGTH**3)*  (F(I)**2+4.*F(I)*B(I)+B(I)**2)*A(I)*RHO(I)/    &
       &((36.*(F(I)+B(I)))*.5)+(A(I)**2/12.+(D(I)+XL(I))**2)*W(I)
        XT = F(I)*LGTH/(B(I)-F(I))
        IF(B(I).LT.F(I)) XT=B(I)*LGTH/(B(I)-F(I))
        TNTAU = C(I)/LGTH
        FETAN = (F(I)/2.- B(I)/2.+ C(I))/LGTH
        AFTAN = (C(I)-F(I)/2.0+B(I)/2.0)/LGTH
        IF (B(I).LT.F(I)) LGTH = 0.0
        KI=ABS(RHO(I)*A(I)*(FETAN**3-AFTAN**3)/3.0)
        KII = (B(I)*A(I)*(RHO(I))/(LGTH+XT))*(TNTAU*(LGTH+XT-XL(I)))**2
        KII=ABS(KII)
        IF(B(I).LT.F(I)) LGTH = (1.0-F(I)/B(I))*XT
        XTLGTH = LGTH + ABS(XT)
        IXXCG(I)=(KI*(XTLGTH**4-XT**4)/2.0)-KII*(XTLGTH**2-XT**2)+ A(I)**2*W(I)/12.
        IZZCG(I) = IXXCG(I) + IYYCG(I)-A(I)**2*W(I)/6.
        XL(I) = LGTH + D(I)
        LGTH=SQRT( (XJ(I)-XI(I))**2 + (YJ(I)-YI(I))**2 + (ZJ(I)-ZI(I))**2 )
      
      CASE(13)   ! curved thin wall panel *shape 13* (Ref -- R Hull)
        w(i)= 2.*lgth*a(i)*b(i) *c(i) *rho(i)
        IF (rho(i) .GT. 0.4) w(i)=rho(i)
        rho(i) = w(i)/(2.*lgth*a(i)*b(i)*c(i))
        ki=2.*c(i)*a(i)*b(i)*lgth*rho(i)
        kii=lgth**2/12.
        xl(i)=lgth/2.
        ixxcg(i)=lgth*rho(i)*a(i)**3*b(i)*(2.*c(i)-(2.*SIN(c(i))**2)/c(i))
        iyycg(i)=ki*(a(i)**2*(c(i)-SIN(c(i))*COS(c(i)))/(2.*c(i)) +kii)
        izzcg(i)=ki*(a(i)**2*(c(i)+SIN(c(i))*COS(c(i))-2.*SIN(c(i))**2/c(i))/(2.*c(i))+kii)

    END SELECT
    
    WRITE(DBG,FMTD) 'Computed IxxCG,IyyCG,IzzCG:', ixxcg(i),iyycg(i),izzcg(i)
    WRITE(DBG,FMTD) 'mass', w(i)
    WRITE(DBG,FMTD) 'xl', xl(i)

!     BEGIN DIRCOS

60  IF ( ABS(XK(I))==0.0 .AND. ABS(YK(I))==0.0 .AND. ABS(ZK(I))==0.0 ) THEN
      XK(I)=XI(I)-(YJ(I)-YI(I))
      YK(I)=YI(I)+(XJ(I)-XI(I))
      ZK(I)=ZI(I)
      IF((YJ(I).EQ.YI(I)).AND.(XJ(I).EQ.XI(I))) YK(I)=LGTH+YI(I)
    END IF
    lx=(xj(i)-xi(i))/lgth
    mx=(yj(i)-yi(i))/lgth
    nx=(zj(i)-zi(i))/lgth
    t1= xk(i)-xi(i)
    t2= yk(i)-yi(i)
    t3= zk(i)-zi(i)
    lz= mx*t3-t2*nx
    mz=nx*t1-t3*lx
    nz=t2*lx-t1*mx
    t4= SQRT(lz**2+ mz**2 +nz**2)
    lz=lz/t4
    mz=mz/t4
    nz=nz/t4
    ly= mz*nx-mx*nz
    my=nz*lx-nx*lz
    ny=mx*lz-lx*mz
    WRITE(DBG,FMTD) 'lgth', lgth
    WRITE(DBG,FMTD) 'lx,ly,lz', lx,ly,lz
    WRITE(DBG,FMTD) 'mx,my,mz', mx,my,mz
    WRITE(DBG,FMTD) 'nx,ny,nz', nx,ny,nz


! Rotate component moment of inertia to system coordinates
    ixxco(i)=(ixxcg(i)*lx**2 + iyycg(i)*ly**2 + izzcg(i)*lz**2)/CONS
    iyyco(i)=(ixxcg(i)*mx**2 + iyycg(i)*my**2 + izzcg(i)*mz**2)/CONS
    izzco(i)=(ixxcg(i)*nx**2 + iyycg(i)*ny**2 + izzcg(i)*nz**2)/CONS
    ixyco(i)=(ixxcg(i)*lx*mx + iyycg(i)*ly*my + izzcg(i)*lz*mz)/CONS
    ixzco(i)=(ixxcg(i)*lx*nx + iyycg(i)*ly*ny + izzcg(i)*lz*nz)/CONS
    iyzco(i)=(ixxcg(i)*mx*nx + iyycg(i)*my*ny + izzcg(i)*mz*nz)/CONS

! Calculate component center of mass coordinates and write out
    xcg(i)=xi(i)+xl(i)*lx
    ycg(i)=yi(i)+xl(i)*mx
    zcg(i)=zi(i)+xl(i)*nx
!!!      WRITE(6,300)ITEM(I),DES(I),W(I),IXXCO(I),IYYCO(I),IZZCO(I),        &
!!!     & XCG(I),YCG(I),ZCG(I),IXYCO(I),IYZCO(I),IXZCO(I)
!!!  300 FORMAT(I5,5X,A18,6F11.5,4F9.5/)
    WRITE(DBG,*) 'Component ', i, "  weight=", w(i)
    WRITE(DBG,FMTD) 'IxxCO,IyyCO,IzzCO=', ixxco(i),iyyco(i),izzco(i)
    WRITE(DBG,FMTD) 'xcg,ycg,zcg=', xcg(i),ycg(i),zcg(i)
    WRITE(DBG,FMTD) 'IxyCO,IyzCO,IxzCO=', ixyco(i),iyzco(i),ixzco(i)

    WRITE(OUT,*) 'Component ', i, '   '//TRIM(description(i))
    WRITE(OUT,FMTD) 'Weight, lbs.      ', w(i)
    WRITE(OUT,FMTD) 'IxxCO,IyyCO,IzzCO=', ixxco(i),iyyco(i),izzco(i)
    WRITE(OUT,FMTD) 'xcg,ycg,zcg=      ', xcg(i),ycg(i),zcg(i)
    WRITE(OUT,FMTD) 'IxyCO,IyzCO,IxzCO=', ixyco(i),iyzco(i),ixzco(i)

! Calculate system weight, second moment at orgin and c.g.
      tw=tw+w(i)
      ixxo=ixxo+ixxco(i)+w(i)*(zcg(i)**2+ycg(i)**2)/cons
      iyyo=iyyo+iyyco(i)+w(i)*(zcg(i)**2+xcg(i)**2)/cons
      izzo=izzo+izzco(i)+w(i)*(xcg(i)**2+ycg(i)**2)/cons
      xmom=xmom+w(i)*xcg(i)
      ymom=ymom+w(i)*ycg(i)
      zmom=zmom+w(i)*zcg(i)
  END DO    ! was 50    **** End of big loop on index i


! Compute system c.g. coordinates
  xbar=xmom/tw
  ybar=ymom/tw
  zbar=zmom/tw

  WRITE(OUT,*) ' '
  WRITE(OUT,*) ' SYSTEM DATA LISTED BELOW'
  WRITE(OUT,*) ' WT=LBS, INERTIAS=SLUGS FT SQUARED'
  WRITE(OUT,*) ' C.G.=INS, SECOND MOMENT=SLUG FT SQUARED)'
  WRITE(OUT,*) ' '
!  WRITE(OUT,*) '   SYS WT     IXXO     IYYO     IZZO     XBAR     YBAR     ZBAR'
!  WRITE(OUT,'(7F12.3)' ) TW,IXXO,IYYO,IZZO,XBAR,YBAR,ZBAR
  WRITE(OUT,FMTD) " System weight   ", tw
  WRITE(OUT,FMTD) " Ixx0,Iyy0,Izz0  ", ixxo,iyyo,izzo
  WRITE(OUT,FMTD) " xbap,ybar,zbar  ", xbar,ybar,zbar

!     TRANSFER MASS PROPERTIES TO SYSTEM C.G., SUM AND WRITE OUT
  DO i=1,nGood
    delx=xcg(i)-xbar
    dely=ycg(i)-ybar
    delz=zcg(i)-zbar
    ixx=ixx + ixxco(i) + w(i)*(dely**2+delz**2)/cons
    iyy=iyy + iyyco(i) + w(i)*(delx**2+delz**2)/cons
    izz=izz + izzco(i) + w(i)*(dely**2+delx**2)/cons
    ixy=ixy + ixyco(i) - w(i)*delx*dely/cons
    ixz=ixz + ixzco(i) - w(i)*delx*delz/cons
    iyz=iyz + iyzco(i) - w(i)*dely*delz/cons
  END DO
  WRITE(OUT,*) ' '
!  WRITE(OUT,*) '   IXX     IYY     IZZ     IXY     IXZ    IYZ'
!  WRITE(OUT,'(6F12.5)') IXX,IYY,IZZ,IXY,IXZ,IYZ
  WRITE(OUT,FMTD) " Ixx,Iyy,Izz ", ixx,iyy,izz
  WRITE(OUT,FMTD) " Ixy,Ixz,Iyz ", ixy,ixz,iyz

!     COMPUTE INERTIAS (EIGENVALUES) ABOUT PRINCIPAL AXES AND EACH AXIS
!     DIRECTION COSINES (EIGENVECTORS) AND WRITE OUT
  WRITE(OUT,*) ' '
  WRITE(OUT,*) ' inertias (eigenvalues) about system principal axes with axis'
  WRITE(OUT,*) ' direction cosines (eigenvectors) relating the principal axes'
  WRITE(OUT,*) ' to the x, y, and z system axes in that sequence'
  WRITE(OUT,*) ' '
  WRITE(OUT,*) ' '

!!!  max=3
!!!  n=3
  arr(1,1)=ixx
  arr(1,2)=ixy
  arr(2,1)=arr(1,2)
  arr(3,1)=ixz
  arr(1,3)=arr(3,1)
  arr(2,2)=iyy
  arr(2,3)=iyz
  arr(3,2)=arr(2,3)
  arr(3,3)=izz

  WRITE(DBG,*) 'arr matrix'
  DO i=1,3
    WRITE(DBG,'(3ES15.5)') arr(i,1:3)
  END DO

  CALL Symql(arr,e(1:3),crr(1:3),ierr)
  IF (ierr/=0) THEN
    WRITE(*,*) 'Error!  ierr=', ierr
    WRITE(OUT,*) 'Error!  ierr=', ierr
  ELSE  
    DO j=1,3
      WRITE(OUT,FMT1) j, e(j)
      WRITE(OUT,FMT2) j, arr(1:3,j)
    END DO
  END IF

  WRITE(*,*) 'Normal stop. File massprop.out added to your directory.'
  STOP
END Program MassProperties   ! =================================================
