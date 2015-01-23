!+
MODULE VmacoProcedures
! ------------------------------------------------------------------------------

! AUTHORS  - J. D. Frick, MDAC/HOUSTON
!            Ralph L. Carmichael, Public Domain Aeronautical Software

! REVISION HISTORY                                                  
!   DATE  VERS PERSON  STATEMENT OF CHANGES   
! 1986-1987      JDF   Original coding
! 1998           RLC   Converted to free format
! 25Oct05  0.1   RLC   Collected procedures into a module
! 30Oct05  0.2   RLC   Added INTENT to all procedure arguments
! 01Nov05  0.3   RLC   Replaced all calls of DotPrd with intrinsic DOT_PRODUCT
! 01Nov05  0.3   RLC   Replaced all calls of Abig with intrinsic MAX
! 14Nov08  0.31  RLC   Changed coding of depend in Chkdep
! 23Nov08  0.5   RLC   Removed nsize from argument list of Vmaco and NX from Matin


! REFERENCES

!  Equations And Algorithm Cited From "The Calculation Of Feasible
!  Points For Linearly Constrained Optimization Problems" 
!  Aere-R.6354 By R. Fletcher

!  Algorithm And Equations Cited Based On The Paper "A Fast
!  Algorithm For Nonlinearly Constrained Optimization
!  Calculations" By M.J.D. Powell At The 1977 Dundee
!  Conference On Numerical Analysis

!  Equations And Algorithm Cited From "A General Quadratic
!  Algorithm" Aere-Tp.401 By R. Fletcher

IMPLICIT NONE



  INTEGER,PARAMETER:: MAXC=5   ! max size of arrays dealing with constraints
  INTEGER,PARAMETER:: MAXX=5   ! max size of arrays dealing with controls
  INTEGER,PARAMETER:: NCOM=15  ! size of array comvma
  INTEGER,PARAMETER:: NKOM=15  ! size of array komvma
  INTEGER,PARAMETER:: MAXACT = MAXX + MAXX + MAXC ! maximum size of array dealing with iactiv array
  INTEGER,PARAMETER:: MAXN   = MAXX - 1 ! maximum size of arrays dealing with control variables

! Variables from old COMMON /COMVMA/
  REAL,DIMENSION(NCOM):: comvma
  REAL:: amdamx ! the maximimum Lagrange multiplier
  REAL:: cac ! curvature of the constraint vector ct (= ct*a*c)
  REAL:: cc  ! particular value of cstarc which corresponds to
             ! the constraint being removed
  REAL:: ccgcac ! required direction of search
  REAL:: chc ! product of the normal of the constraint and hdotc
  REAL:: cnear ! distance to the nearest currently inactive constraint
  REAL:: gdelta ! scalar product of gradient of the objective function and delta.
  LOGICAL:: passiv ! indicates a constraint has been removed to get a
                   ! search direction
  REAL:: y   ! product of chc*cac + cc**2
  REAL:: fdelta ! absolute change in the objective function f. Checked
        ! against accuracy desired for convergence criteria.
  REAL:: dflsav ! saved value of dflsa
  REAL:: flsa ! initial value of fls. corresponds to phi(0) in eqn(4.7)

! Variables from old COMMON /CVMAC/
  REAL,DIMENSION(MAXC):: d
  REAL,DIMENSION(MAXC):: vlam
  REAL,DIMENSION(MAXC):: vmu

! Variables from old COMMON /CVMACT/
  REAL,DIMENSION(MAXACT):: resid  ! product of the ith residual in the basis and recipr

! Variables from old COMMON /CVMAN/
  REAL,DIMENSION(MAXN):: glf ! gradient of the lagrange function. See eqns (3.2) and (3.3).
  REAL,DIMENSION(MAXN):: glfa ! previous gradient of the lagrange function
  REAL,DIMENSION(MAXN):: oxa ! saved value of the variable xa
  REAL,DIMENSION(MAXN):: xa ! initial guess at control varible values

! Variables from old COMMON /CVMAX1/
  REAL,DIMENSION(MAXX,MAXX):: a,b
  REAL,DIMENSION(MAXX):: bdl,bdu ! lower and upper bounds to the controls and 
                                 ! artificial variable. See eqn (1a)
  REAL,DIMENSION(MAXX):: delta ! change in the control variable. Corresponds to
                               ! search direction times the step-length.
  REAL,DIMENSION(MAXX):: gm ! negative of the gradient vector (= -g). Corresponds
                            ! to b in eqns (4),(5),(6), and (8).
  REAL,DIMENSION(MAXX,MAXC):: cinvrs ! inverse matrix corresponding to the normals
                                     ! of the constraints in the basis. see eqn (2)
  REAL,DIMENSION(MAXX,MAXC):: cn ! matrix of contstraint normals
  REAL,DIMENSION(MAXX,MAXC):: cstar ! inverse matrix corresponding to the normals 
                                    ! of the constraints in the basis. see eqn (11)

  EQUIVALENCE(a,b)   ! not sure why?  RLC
  EQUIVALENCE(cstar(1,1), cinvrs(1,1) )  ! not sure why?  RLC

! Variables from old COMMON /CVMAX2/
  REAL,DIMENSION(MAXX):: ghat ! current gradient vector of the objective function.
                              ! see eqn (9)
  REAL,DIMENSION(MAXX,MAXX):: h  ! Can be considered as a reduced inverse 
                                 ! Hessian operator appropriate to the 
                                 ! manifold formed by the intersection of the
                                 ! constraints. See eqn (12)

! Variables from old COMMON /CVMAX3/
  REAL,DIMENSION(MAXX):: ac !     product of the a matrix and ct vector
  REAL,DIMENSION(MAXX):: ckac ! product of the kth row of cstar and ac matrix
  REAL,DIMENSION(MAXX):: cstarc  ! product of the cstar operator and normal of the constraint
  REAL,DIMENSION(MAXX):: cstep ! component of the gradient along the search direction
  REAL,DIMENSION(MAXX):: ct ! constraint vector, corresponding to the largest
                            ! lagrange multiplier, which is to be removed from the basis.
  REAL,DIMENSION(MAXX):: hdotc  ! product of the h operator and normal of the constraint

! Variables from old COMMON /KOMVMA/
  INTEGER,DIMENSION(NKOM):: komvma
  INTEGER:: icnear  ! index of the nearest currently inactive constraint
  INTEGER:: iter    ! counter for number of iterations completed
  INTEGER:: lamrow  ! index of the constraint, corresponding to the largest
               ! Lagrange multiplier, which is to be removed from the basis
  INTEGER:: mtotal  ! total number of constraints (= m+2*nvar). Includes the
                    ! input constraints and the upper and lower boundaries
  INTEGER:: nc      ! number of constraints in the basis. Initially
                    ! represents the number of designated constraints in
                    ! the basis.
  INTEGER:: nvar    ! (=n+1) number of controls plus the artificial variable.
  INTEGER:: nvar2   ! twice the value of nvar. (= 2*(n+1))
  INTEGER:: nvma    ! number of calls to vmaco
  INTEGER:: nvmitr  ! value of nvma at the start of an iteration
  INTEGER:: icrow   ! index of the row of cinvrs which is to be removed from the basis
  INTEGER:: ibugfg  ! flag to indicate extra printout for debugging purposes
                    !  (= 0) no extra printout; (= 1) extra printout

! Variables from old COMMON /KVMACT/
  INTEGER,DIMENSION(MAXACT):: iactiv ! status of each of the mtotal constraints
                                     !  =-1 equality
                                     !  = 0 active
                                     !  = 1 inactive
                                     !  = 2 violated

! Variables from old COMMON /KVMAX/
  INTEGER,DIMENSION(MAXX):: ibasis  ! index numbers of the constraints in the basis

! other module variables
!  INTEGER,PRIVATE:: IECHO = 6   ! file unit number to which printout is written
  REAL,PARAMETER,PRIVATE:: ZERO = 0.0
  REAL,PARAMETER,PRIVATE:: FP1  = 1.0
  REAL,PARAMETER,PRIVATE:: HALF = 0.5

  PRIVATE:: Add
  PRIVATE:: Bugaid
  PRIVATE:: Changp
  PRIVATE:: Chkdep
  PRIVATE:: Ckvrtx
  PRIVATE:: Colum3
  PRIVATE:: Defer
  PRIVATE:: Desig
  PRIVATE:: Exchng
  PRIVATE:: Gradlf
  PRIVATE:: Initvm
  PRIVATE:: Lgrang
  PRIVATE:: Lnsrch
  PRIVATE:: Matin
  PRIVATE:: Move
  PRIVATE:: Newbas
  PRIVATE:: Newdir
  PRIVATE:: Nodesg
  PRIVATE:: Quad
  PRIVATE:: Remove
  PRIVATE:: Reset
  PRIVATE:: Serch
  PRIVATE:: Slnear
  PRIVATE:: Stepc
  PRIVATE:: Stinit
  PRIVATE:: Stovmu
  PRIVATE:: Updatb
  PRIVATE:: Vertex
  PUBLIC::  Vmaco  ! in normal usage, the only procedure to be called
  PRIVATE:: Vmamsg
  PRIVATE:: Vmaout
  PRIVATE:: Vmaset

CONTAINS

!+
SUBROUTINE Add()
! ------------------------------------------------------------------------------
! PURPOSE - Adds a constraint to the basis. see eqns (13a) and (13b)
! NOTE - Equations and algorithm cited from "A general quadratic algorithm" 
!  aere-tp.401 by r. fletcher
! called by Serch

! Global variables used but not changed:
!USE GlobalVmaco,ONLY: chc,cstarc,hdotc,icnear,nvar
! Global variables modified in Subroutine Addd
!USE GlobalVmaco,ONLY: cstar,h,ibasis,nc

  REAL:: ccdchc  ! value of ck'c/chc
  REAL:: hcdchc  ! value of hdotc/chc
  INTEGER:: i,j
!-------------------------------------------------------------------------------
      IF (NC .GT. 0) THEN
        DO 20 I=1,NC
          CCDCHC = CSTARC(I)/CHC
          DO 10 J=1,NVAR
            CSTAR(I,J) = CSTAR(I,J) - CCDCHC*HDOTC(J)
   10     END DO
   20   END DO
      END IF

      nc = nc + 1
      ibasis(nc) = icnear

!...the new constraint is add to the cstar operator

      DO 30 i=1,nvar
        cstar(nc,i) = hdotc(i)/chc
   30 END DO

!...the h operator is updated

      IF (NC .GE. NVAR) THEN
        DO 50 I=1,NVAR
          DO 40 J=1,NVAR
            H(I,J) = ZERO
   40     END DO
   50   CONTINUE
      ELSE
        DO 70 I=1,NVAR
          HCDCHC = HDOTC(I)/CHC
          DO 60 J=1,I
            H(I,J) = H(I,J) - HCDCHC*HDOTC(J)
            H(J,I) = H(I,J)
   60     END DO
   70   END DO
      END IF

  RETURN
END Subroutine Add   ! ---------------------------------------------------------

!+
SUBROUTINE BUGAID (IBUG,NUM,XNUM)
! ------------------------------------------------------------------------------
! PURPOSE - Provides extra printout for clarification and debugging purposes.
  INTEGER,INTENT(IN):: IBUG   ! flag indicating which variable values are to be output
  INTEGER,INTENT(IN):: NUM    ! number of data array points to be printed
  REAL,INTENT(IN):: XNUM   ! values to be printed

  INTEGER:: i
  INTEGER:: ndummy
!!!      PARAMETER ( IECHO  =    6 )
  REAL,DIMENSION(1):: xdummy   ! added by RLC, 30Oct2005
!-------------------------------------------------------------------------------

      IF (IBUG .EQ. 1) THEN
!
!...    THE DERIVATIVE OF THE LINE SEARCH FUNCTION IS OUTPUT
!
        WRITE(*,510) XNUM
!
      ELSE IF (IBUG .EQ. 2) THEN
!
!...    THE B MATRIX IS OUTPUT
!
        DO 10 I=1,NUM
          CALL COLUM3 (1,NUM,XDUMMY,I)
   10   CONTINUE
!
      ELSE IF (IBUG .EQ. 3) THEN
!
!...    LINE SEARCH FUNCTION VALUE IS OUTPUT
!
        WRITE(*,530) XNUM
!
      ELSE IF (IBUG .EQ. 4) THEN
!
!...    STEP LENGTH IS OUTPUT
!
        WRITE(*,540) XNUM
!
      ELSE IF (IBUG .EQ. 5) THEN
!
!...    CHANGE IN THE CONTROL VARIABLE (DELTA) IS OUTPUT
!
        CALL COLUM3 (2,NUM,XDUMMY,NDUMMY)
!
      ELSEIF (IBUG .EQ. 6) THEN
!
!...    THE MAXIMUM LAGRANGE MULTIPLIER IS OUTPUT
!
         WRITE(*,560) XNUM
!
       ELSE IF (IBUG .EQ. 7) THEN
!
!...     CURRENT GRADIENT OF THE OBJECTIVE FUNCION (GHAT) IS OUTPUT
!
         CALL COLUM3 (3,NUM,XDUMMY,NDUMMY)
!
      ELSE IF (IBUG .EQ. 8) THEN
!
!...    GRADIENT OF THE LINE SEARCH FUNCTION (GLF) IS OUTPUT
        CALL COLUM3 (4,NUM,XDUMMY,NDUMMY)
!
      ELSE IF (IBUG .EQ. 9) THEN

!...    VECTOR OF LAGRANGE MULTIPLIERS (VLAM) ARE OUTPUT
        CALL COLUM3 (5,NUM,XDUMMY,NDUMMY)

      ELSE IF (IBUG .EQ. 10) THEN

!...    RMS VALUE OF THE GRADIENT OF THE LINE SEARCH
!...    IS OUTPUT
        WRITE(*,570) XNUM

      END IF

  RETURN

  510 FORMAT (' DFLSA =',1PE16.7)
  530 FORMAT (' FLS =',1PE16.7 )
  540 FORMAT (' STEPL =',1PE16.7 )
  560 FORMAT (' LAMBDA MAX =',1PE16.7 )
  570 FORMAT (' GLFMAG =',1PE16.7 )
END Subroutine Bugaid   ! ------------------------------------------------------

!+
SUBROUTINE CHANGP()
! ------------------------------------------------------------------------------
! PURPOSE - Exchanges the new constraint with the passive constraint.
!D7
!D7   LIMITATIONS/ASSUMPTIONS
!D7
!D7      EQUATIONS AND ALGORITH CITED FROM "A GENERAL QUADRATIC
!D7      ALGORITHM" AERE-TP.401 BY R. FLETCHER
! called by Serch
! calls DotPrd


! Global variables used but not changed
!USE GlobalVmaco,ONLY: ac,cac,cc,chc,cstar,ct,hdotc,icnear,lamrow,MAXX,nc,nvar,y
! Global variables modified by Subroutine ChangP
!USE GlobalVmaco,ONLY: ckac,cstarc,h,ibasis
IMPLICIT NONE
  REAL,DIMENSION(MAXX):: udivy !  value of u divided by y where u = c*(c'hkc)-hkc(c'c*)
  REAL,DIMENSION(MAXX):: wdivy !  value w divided by y where  w=hkc(c*'ac*)+c*(c'c*)
  INTEGER:: i,j
  INTRINSIC:: DOT_PRODUCT
!-------------------------------------------------------------------------------

!      LSIZE = 1 + MAXX*(NVAR - 1)
      DO I=1,NC
!!!        CKAC(I) = DOTPRD (NVAR,CSTAR(I:,1),LSIZE,MAXX,AC(1:),NVAR,1)
        ckac(i)=DOT_PRODUCT(cstar(i,1:nvar), ac(1:nvar) )
      END DO

      DO 20 I=1,NVAR
        UDIVY(I) = (CHC*CT(I) - CC*HDOTC(I))/Y
        WDIVY(I) = (CAC*HDOTC(I) + CC*CT(I))/Y
   20 END DO

!...update the h operator as per formula eqn (15b)

      DO 40 I=1,NVAR
        DO 30 J=1,I
          H(I,J) = H(I,J) + UDIVY(I)*CT(J) - WDIVY(I)*HDOTC(J)
          H(J,I) = H(I,J)
   30   END DO
   40 END DO

!...update the cstar operator as per formula (15a)

      CSTARC (LAMROW) = CSTARC (LAMROW) - FP1

      DO 60 I=1,NC
        DO 50 J=1,NVAR
          CSTAR(I,J) = CSTAR(I,J)-CSTARC(I)*WDIVY(J)-CKAC(I)*UDIVY(J)
   50   END DO
   60 END DO

      IBASIS(LAMROW) = ICNEAR
  RETURN
END Subroutine Changp   ! ------------------------------------------------------

!+
SUBROUTINE Chkdep (depend)
! ------------------------------------------------------------------------------
! PURPOSE - Checks dependency of the new constraint.
!D7      EQUATIONS AND ALGORITHM CITED FROM  "A GENERAL QUADRATIC
!D7      ALGORITHM" AERE-TP.401 BY R. FLETCHER
! called by Serch
! calls DotPrd

! Global variables used but not changed
!USE GlobalVmaco,ONLY: amdamx,cac,ccgcac,chc,cstarc,ghat,hdotc,lamrow,nvar
! Global variables modified by Subroutine Chkdep
!USE GlobalVmaco,ONLY: cc,y
IMPLICIT NONE
  LOGICAL,INTENT(OUT):: depend   ! flag to indicate the dependency of the new constraint

  REAL:: ghc   ! product of ghat and hdotc
  REAL:: t1,t2 ! temporary products to check dependency
  INTRINSIC:: DOT_PRODUCT
!-------------------------------------------------------------------------------
      cc = cstarc (lamrow)
      y = chc*cac + cc*cc
!!!      GHC = DOTPRD ( NVAR,GHAT(1:),NVAR,1,HDOTC(1:),NVAR,1)
      ghc=DOT_PRODUCT(ghat(1:nvar), hdotc(1:nvar) )
      t1 = ccgcac*y
      t2 = chc*(amdamx - ccgcac*cac) + ghc*cc
!
!      IF (T1 .LT. T2) THEN
!        DEPEND = .FALSE.
!      ENDIF
  depend = t1 >= t2   ! changed by RLC 14Nov2008. depend might be undefined
  RETURN
END Subroutine ChkDep   ! ------------------------------------------------------

!+
SUBROUTINE ckvrtx (nviol, snorm)
! ------------------------------------------------------------------------------
! PURPOSE - Checks the vertex against the inactive constraints to make sure it
! is a feasible point.

! Global variables used but not changed
!USE GlobalVmaco,ONLY: bdl,bdu,cn,d,delta,iactiv,MAXC,mtotal,nvar,nvar2
! Global variables modified by Subroutine CkVrtx
!USE GlobalVmaco,ONLY: resid
IMPLICIT NONE
  INTEGER,INTENT(OUT):: nviol   ! number of violated constraints
  REAL,INTENT(OUT),DIMENSION(MAXC):: snorm   ! sum of the violated constraint normals
  INTEGER:: i,j

  INTRINSIC:: DOT_PRODUCT
! called by Vertex
!-------------------------------------------------------------------------------
      NVIOL = 0
      DO 10 I=1,NVAR
        SNORM(I) = ZERO
   10 END DO

!...the constraint residuals, the number of violated constraints, and
!...the sum of their normals are calculated

      DO 30 I=1,MTOTAL
        IF (IACTIV(I) .GE. 1) THEN

!...      the constraint is inactive and checked for violation

          IF (I .LE. NVAR) THEN
            RESID(I) = DELTA(I) - BDL(I)
          ELSE IF (I .LE. NVAR2) THEN
            J = I-NVAR
            RESID(I) = BDU(J) - DELTA(J)
          ELSE
            J = I-NVAR2
!!!            RESID(I) = DOTPRD (NVAR,CN(1:,J),NVAR,1,DELTA(1:),NVAR,1)
            resid(i)=DOT_PRODUCT(cn(1:nvar,j), delta(1:nvar) ) 
            RESID(I) = RESID(I) -D(J)
          END IF

          IF (RESID(I) .LT. 0.) THEN
!
!...        this constraint has been violated
!
            NVIOL = NVIOL + 1
            IACTIV(I) = 2

            IF (I .LE. NVAR) THEN
              SNORM(I) = SNORM(I) + FP1
            ELSE IF (I .LE. NVAR2) THEN
              J = I - NVAR
              SNORM(J) = SNORM(J) - FP1
            ELSE
              DO 20 J=1,NVAR
                SNORM(J) = SNORM(J) + CN(J,I-NVAR2)
   20         END DO
            END IF
          END IF
        END IF
   30 END DO
  RETURN
END Subroutine CKVRTX   ! ------------------------------------------------------

!+
SUBROUTINE Colum3(iform,num,xnum,numx)
! ------------------------------------------------------------------------------
! PURPOSE - Provides printout data in columns of 3.

! Global variables used but not changed
!USE GlobalVmaco,ONLY: b=>a,cn,delta,ghat,glf,vlam
! No global variables modified
IMPLICIT NONE
  INTEGER,INTENT(IN):: IFORM   ! FLAG INDICATING WHICH DATA ARRAY TO PRINT
  INTEGER,INTENT(IN):: NUM     ! NUMBER OF DATA ARRAY POINTS TO BE PRINTED
  INTEGER,INTENT(IN):: NUMX    ! EXTRA INTEGER NUMBER USED FOR MATRIX COLUMN IDENTIFICATION
  REAL,INTENT(IN),DIMENSION(:):: XNUM   ! A DATA ARRAY TO BE PRINTED


  INTEGER:: i
  INTEGER:: iend  !  last point of data array print in a row
  INTEGER:: irpt  !  number of pieces of data to print in a row
  INTEGER:: istart ! start point of data array print in a row
!D7      LIMITED TO 72 COLUMNS OF OUTPUT
! called by BugAid and VmaOut

      INTEGER,PARAMETER:: NFORM=9, NREPEAT=3

      CHARACTER(LEN=1),DIMENSION(NREPEAT):: REPEAT
      DATA REPEAT / '1','2','3' /
  CHARACTER(LEN=45),DIMENSION(NFORM):: FORM = (/ &
    "(Z(1X, 'B(',     I2, ',', I2, ')=', ES14.7)) ", &
    "(Z(1X, 'DELTA(', I2, ')=', ES13.6))          ", &
    "(Z(1X, 'GHAT(',  I2, ')=', ES14.7))          ", &
    "(Z(1X, 'GLF(',   I2, ')=', ES14.7))          ", &
    "(Z(1X, 'VLAM(',  I2, ')=', ES14.7))          ", &
    "(Z(1X, 'X(',     I2, ')=',  ES16.7))         ", &
    "(Z(1X, 'C(',     I2, ')=', ES16.7))          ", &
    "(Z(1X, 'G(',     I2, ')=', ES16.7))          ", &
    "(Z(1X, 'CN(',    I2, ',', I2, ')=', ES14.7)) " /)



!-------------------------------------------------------------------------------
if (num==0) write(*,*) "num=0 in Colum3"
      ISTART = 1
      IF (NUM .LT. 3) THEN
        IRPT = NUM
        IEND = NUM
      ELSE
        IRPT = 3
        IEND = 3
      ENDIF
!
   10 CONTINUE

      FORM (IFORM)(2:2) = REPEAT (IRPT)
!  write(*,*) num
!  write(*,*) form(iform)

      IF (IFORM .EQ. 1) THEN
        WRITE(*,FORM(1))(I,NUMX,B(I,NUMX),I=ISTART,IEND)
      ELSE IF (IFORM .EQ. 2) THEN
        WRITE(*,FORM(2))(I,DELTA(I),I=ISTART,IEND)
      ELSE IF (IFORM .EQ. 3) THEN
        WRITE(*,FORM(3))(I,GHAT(I),I=ISTART,IEND)
      ELSE IF (IFORM .EQ. 4) THEN
        WRITE(*,FORM(4))(I,GLF(I),I=ISTART,IEND)
      ELSE IF (IFORM .EQ. 5) THEN
        WRITE(*,FORM(5))(I,VLAM(I),I=ISTART,IEND)
      ELSE IF (IFORM .EQ. 9) THEN
        WRITE(*,FORM(9))(I,NUMX,CN(I,NUMX),I=ISTART,IEND)
      ELSE
        WRITE(*,FORM(IFORM))(I,XNUM(I),I=ISTART,IEND)
      ENDIF
!
      IF (IEND .LT. NUM) THEN
        ISTART = IEND + 1
        IF (NUM .LE. IEND+3) THEN
          IEND = NUM
          IRPT = NUM - ISTART + 1
        ELSE
          IEND = IEND + 3
        ENDIF
        GO TO 10
      ENDIF

  RETURN
END Subroutine Colum3   ! ------------------------------------------------------

!+
SUBROUTINE DEFER()
! ------------------------------------------------------------------------------
! PURPOSE - Calculate the maximum Lagrange multiplier and the gradient vector as
!  though the constraint is being removed from the augmented basis.

! Global variables used but not changed
!USE GlobalVmaco,ONLY: a,delta,gm,nvar
! Global variables modified by Subroutine Defer
!USE GlobalVmaco,ONLY: amdamx,ct,ghat
IMPLICIT NONE
  INTEGER:: i
!D7
!D7      EQUATIONS AND ALGORITHM CITED FROM "A GENERAL QUADRATIC
!D7      ALGORITHM" AERE-TP.401 BY R. FLETCHER


! called by Serch
! calls DotPrd
  INTRINSIC:: DOT_PRODUCT
!-------------------------------------------------------------------------------
!      LSIZE = 1 + MAXX*(NVAR - 1)
      DO I=1,NVAR
!!!        GHAT(I) = DOTPRD(NVAR,A(I:,1),LSIZE,MAXX,DELTA(1:),NVAR,1) - GM(I)
        ghat(i)=DOT_PRODUCT(a(i,1:nvar), delta(1:nvar)) - gm(i)
      END DO
!
!!!      AMDAMX = -DOTPRD (NVAR,GHAT(1:),NVAR,1,CT(1:),NVAR,1)
      amdamx=-DOT_PRODUCT(ghat(1:nvar), ct(1:nvar) )

  RETURN
END Subroutine Defer   ! -------------------------------------------------------

!+
SUBROUTINE DESIG (MEQ)
! ------------------------------------------------------------------------------
! PURPOSE - There are designated constraints in the basis. A feasible vertex 
!  must be found with the designated constraints.

! Global variables used but not changed
!USE GlobalVmaco,ONLY: bdl,bdu,cn,delta,ibasis,MAXX,nc,nvar,nvar2

! Global variables modified by Subroutine 
!USE GlobalVmaco,ONLY: cinvrs,iactiv,ibasis

IMPLICIT NONE
  INTEGER,INTENT(IN):: MEQ   !  TOTAL NUMBER OF EQUALITY CONSTRAINTS

  REAL,DIMENSION(MAXX):: ctemp !  temporary array to store rows of cinvrs data
  REAL:: determ
  REAL:: dummy
  INTEGER:: i,j,jx
  INTEGER:: ip !     flag to indentify which element of the projection is the smallest
  INTEGER,DIMENSION(MAXX):: ipaflg ! flag to avoid recalculating p ab-initio at each stage
                                   !  = 0 vector i not in basis update
                                   !  = 1 vector i already in basis update
  INTEGER:: key ! rank of the inverted matrix
  INTEGER:: lnew !   index number of the constraint being brought into the basis
  INTEGER:: ndesx !  counter for adding remaining constraints to the cinvrs matrix
  REAL,DIMENSION(MAXX):: p !      diagonal elements of the projection matrix
  REAL:: pmult ! multiplier determined by whether the vertex is closer
               ! to the lower or upper boundary
  REAL:: psmall ! smallest diagonal element of the projection matrix
  REAL:: u !      value of e(i) - v*vplus*e(i)
  REAL:: uuei !   value of u transpose/u transpose*e(i)
  REAL,DIMENSION(MAXX,MAXX):: v !      normals of the designated constraints
  REAL,DIMENSION(1):: vecdum   ! addedby RLC, 30Oct2005 dummy value for routine matin
  REAL,DIMENSION(MAXX):: vpluse ! vector v+ * e(i)
  REAL,DIMENSION(MAXX):: vvplse ! vector v * v+ * e(i)

! called by Vertex
! calls Matin, Vmamsg and DotPrd
  INTRINSIC:: DOT_PRODUCT
!-------------------------------------------------------------------------------
!  write(*,*) "Entering Desig, nc,nvar=", nc,nvar
!...The normals v of the nc designated constraints in the basis are set up.

  v=ZERO
!  v=-999.9


      DO 30 I=1,NC
        IF (I .LE. MEQ) THEN
          J = IBASIS(I)
          IACTIV(J) = -1
        ELSE
          J = IBASIS(I)
          IACTIV(J) = 0
        ENDIF

        IF (IBASIS(I) .LE. NVAR2) THEN
          DO 10 J=1,NVAR
            V(J,I) = ZERO
   10     CONTINUE

          IF (IBASIS(I) .LE. NVAR) THEN
            V(IBASIS(I),I) = FP1 ! lower boundary constraint is in the basis
          ELSE
            JX = IBASIS(I) - NVAR ! upper boundary constraint is in the basis
            V(JX,I) = -FP1
          END IF
        ELSE
          JX = IBASIS(I) - NVAR2 ! equality or inequality constraint
                                 ! of type eqn (1b) is in the basis
          DO 20 J=1,NVAR
            V(J,I) = CN(J,JX)
   20     CONTINUE
        END IF
   30 END DO

!  write(*,*) 'After 30', v(1:nvar,1:nc)


      IF (NC .EQ. NVAR) THEN
        DO 50 I=1,NVAR
          DO 40 J=1,NVAR
            CINVRS(I,J) = V(I,J)
   40     CONTINUE
   50   CONTINUE
!        write(*,*) 'Why am I here?, nc,nvar=', nc,nvar
        CALL MATIN(CINVRS,NVAR,VECDUM,0,KEY,DETERM)
        IF (KEY .LT. NVAR) THEN
          CALL VMAMSG (15,1,DUMMY)
        ENDIF
        GO TO 999
      ENDIF
!
!...calculate the inverse of (vtranpose * v) and store it temporarily in cinvrs

!  write(*,*) 'Before 60', v(1:nvar,1:nc)
      DO 60 I=1,NC
        DO 55 J=I,NC
!  write(*,*) 'In loop 55',i,j,v(1:nvar,i),v(1:nvar,j)
!!!          CINVRS(I,J) = DOTPRD (NVAR,V(1:,I),NVAR,1,V(1:,J),NVAR,1)
          cinvrs(i,j)=DOT_PRODUCT(v(1:nvar,i), v(1:nvar,j) )
          CINVRS(J,I) = CINVRS (I,J)
   55   CONTINUE
   60 END DO

      IF (NC .EQ. 1) THEN
        CINVRS(1,1) = FP1/CINVRS(1,1)
      ELSE
        CALL MATIN (CINVRS,NC,VECDUM,0,KEY,DETERM)
        IF (KEY .LT. NC) THEN
          CALL VMAMSG (15,1,DUMMY)
        ENDIF
      ENDIF
!
!...calculate the generalized inverse of v (i.e. v+) and store it in cinvrs. 
!   v+ = (1/(vtranspose*v)) * vtranspose

!!!      LSIZE = 1 + MAXX*(NC - 1)
!  write(*,*) 'Before 80', v(1:nvar,1:nc)
      DO 80 I=1,NC
        DO 65 J=1,NC
          CTEMP(J) = CINVRS(I,J)
   65   CONTINUE
        DO 70 J=1,NVAR
!  write(*,*) 'In 70 loop ',j,v(j,1:nc)
!!!          CINVRS(I,J) = DOTPRD(NC,CTEMP(:),NC,1,V(J:,1),LSIZE,MAXX)  ! modified RLC
          cinvrs(i,j)=DOT_PRODUCT(ctemp(1:nc), v(j,1:nc) )
   70   CONTINUE
   80 END DO
!
!...calculate the diagonal elements of the projection matrix
!...(i.e p(i) = diagonal of v*v+) and initialize bound variables.
!

!  write(*,*) 'First ',v(1,1:2)
      DO I=1,NVAR
!!!        P(I) = DOTPRD (NC,V(I:,1),LSIZE,MAXX,CINVRS(1:,I),NC,1) ! modified RLC
        p(i)=DOT_PRODUCT(v(i,1:nc), cinvrs(1:nc,i) )
        IPAFLG(I) = 0
      END DO
      NDESX = NC + 1
!
!...the remaining nvar-nc constraints are added to the cinvrs matrix
!
      DO 160 LNEW=NDESX,NVAR
        PSMALL=FP1
!
!...    calculate the smallest bound e(i) from the diagonal of the
!...    projection matrix and store it in psmall
!
        DO 100 I=1,NVAR
          IF (IPAFLG(I) .EQ. 0) THEN
            IF (P(I) .LT. PSMALL) THEN
              PSMALL = P(I)
              IP = I
            ENDIF
          ENDIF
  100   CONTINUE
!
        IF (DELTA(IP) - BDL(IP) .GT. BDU(IP) - DELTA(IP)) THEN
          PMULT = -FP1
        ELSE
          PMULT = FP1
        ENDIF
!
!...    calculate the vectors v+ * e(i)
!
        DO 110 I=1,NC
          VPLUSE(I) = PMULT*CINVRS(I,IP)
  110   CONTINUE
!
!...calculate vectors v * v+ * e(i)
!

!  write(*,*) 'Before 120', v(1:nvar,1:nc)
!  write(*,*) 'nvar,nc', nvar,nc

        DO 120 I=1,NVAR
!!!  write(*,*) 'In loop ',v(1,1:nc)
          IF (IPAFLG(I) .EQ. 0) THEN
!!!            VVPLSE(I) = -DOTPRD (NC,V(I:,1),LSIZE,MAXX,VPLUSE(1:),NC,1)
            vvplse(i)=-DOT_PRODUCT(v(i,1:nc), vpluse(1:nc) )
          END IF
          CINVRS(I,IP) = ZERO
  120   CONTINUE
!
!...    calculate u = e(i) - v * v+ * e(i)
!
        U = FP1 + VVPLSE(IP)*PMULT
        IPAFLG(IP) = 1
!
!...    update v+
!
        DO 140 I=1,NVAR
          IF (IPAFLG(I) .EQ. 0) THEN
            UUEI = VVPLSE(I)/U
            CINVRS (LNEW,I) = UUEI
            DO 130 J=1,NC
              CINVRS(J,I) = CINVRS(J,I) - VPLUSE(J)*UUEI
  130       CONTINUE
          ENDIF
  140   CONTINUE
!
!...    update the diagonal elements of the projection matrix
!
        DO 150 I=1,NVAR
          IF (IPAFLG(I) .EQ. 0) THEN
            P(I) = P(I) + (VVPLSE(I) * VVPLSE(I)) / U
          ENDIF
  150   CONTINUE
!
!...    the constraint brought into the basis is identified and indexed
!
        CINVRS(LNEW,IP) = PMULT
        IF (PMULT .LT. ZERO) THEN
          IACTIV(IP+NVAR) = 0
          IBASIS(LNEW) = IP + NVAR
        ELSE
          IACTIV(IP) = 0
          IBASIS(LNEW) = IP
        ENDIF
!
        NC = NC + 1
!
  160 END DO
!
  999 CONTINUE

  RETURN
END Subroutine Desig   ! -------------------------------------------------------


!+
SUBROUTINE EXCHNG()
! ------------------------------------------------------------------------------
! PURPOSE - Utilize the simplex method to exchange constraints within the basis.

! Global variables used but not changed
!USE GlobalVmaco,ONLY: cstarc,ct,icnear,lamrow,nvar
! Global variables modified by Subroutine Exchng
!USE GlobalVmaco,ONLY: cstar,ibasis
IMPLICIT NONE
  INTEGER:: i,j
  REAL:: recip
!-------------------------------------------------------------------------------
      recip = FP1/cstarc(lamrow)

      DO 30 i=1,nvar
        IF (i .EQ. lamrow) THEN
          DO 10 j=1,nvar
            cstar(i,j) = cstar(i,j)*recip
   10     END DO
        ELSE
          DO 20 j=1,nvar
            cstar(i,j) = cstar(i,j) - recip*cstarc(i)*ct(j)
   20     END DO
        END IF
   30 END DO

  ibasis(lamrow) = icnear

  RETURN
END Subroutine EXCHNG   ! ------------------------------------------------------

!+
SUBROUTINE Gradlf(g,m,n,iprint)
! ------------------------------------------------------------------------------
! PURPOSE - Calculate the gradient of the Lagrange function.

! Global variables used but not changed
!USE GlobalVmaco,ONLY: cn,vlam

! Global variables modified by Subroutine GradLf
!USE GlobalVmaco,ONLY: glf

IMPLICIT NONE
  REAL,INTENT(IN),DIMENSION(:):: G   !    GRADIENT OF THE OBJECTIVE FUNCTION
  INTEGER,INTENT(IN):: IPRINT   ! FLAG TO INDICATE AMOUNT OF PRINTOUT
!D2              = 0 NO PRINTOUT
!D2              = 1 VALUES OF X, C, AND F ARE PRINTED
!D2              = 2 DEBUGGING PRINTOUT
  INTEGER,INTENT(IN):: M   !    TOTAL NUMBER OF CONSTRAINTS (I.E. CONSTRAINT EQUATIONS)
  INTEGER,INTENT(IN):: N   !    NUMBER OF CONTROL VARIABLES

  REAL:: glfmag ! rms value of the glf vector
  INTEGER:: i,j
  INTEGER:: ibug
  REAL,PARAMETER:: TOL = 1E-38 ! TOLERANCE FOR EQUALITY ON 0.0
  REAL,DIMENSION(1):: xdummy ! ???
!D7
!D7   LIMITATIONS/ASSUMPTIONS
!D7
!D7      ALGORITHM AND EQUATIONS CITED BASED ON THE PAPER "A FAST
!D7      ALGORITHM FOR NONLINEARLY CONSTRAINED OPTIMIZATION
!D7      CALCULATIONS" BY M.J.D. POWELL AT THE 1977 DUNDEE
!D7      CONFERENCE ON NUMERICAL ANALYSIS

  INTRINSIC:: DOT_PRODUCT
!-------------------------------------------------------------------------------
      DO 10 I=1,N
        GLF(I) = G(I)
   10 END DO

      IF (M .GE. 1) THEN
        DO 30 I=1,M
          IF (ABS(VLAM(I)) .GT. TOL ) THEN
            DO 20 J=1,N
              GLF(J) =GLF(J) - CN(J,I)*VLAM(I)
   20       CONTINUE
          ENDIF
   30   CONTINUE
      ENDIF

      IF (IPRINT .EQ. 2) THEN
!!!        GLFMAG = SQRT ( DOTPRD(N,GLF(:),MAXN,1,GLF(:),MAXN,1) )
        glfmag=SQRT(DOT_PRODUCT(glf(1:n), glf(1:n) ))
        IBUG = 10
        CALL BUGAID (IBUG,1,GLFMAG)
        IBUG = 8
        CALL BUGAID (IBUG,N,XDUMMY(1))
        IF (M .GE. 1) THEN
          IBUG = 9
          CALL BUGAID (IBUG,M,XDUMMY(1))
        ENDIF
      ENDIF
!
  RETURN
END Subroutine GradLF   ! ------------------------------------------------------

!+
SUBROUTINE INITVM (ISTAT,N,M)
! ------------------------------------------------------------------------------
! PURPOSE - Initialize variable values needed by Vmaco.

! Global variables modified by Subroutine 
!USE GlobalVmaco,ONLY: b=>a,ibugfg,iter,nvar,nvma,vmu

IMPLICIT NONE
  INTEGER,INTENT(IN OUT):: ISTAT   ! indicate the status of the optimization algorithm.
!              =-2 initial pass through vmaco with the b matrix
!                  already initialized
!              =-1 initial pass through vmaco
!              = 0 continue algorithm calculations
!              = 1 converge within required acurracy
!              > 1 error encountered
  INTEGER,INTENT(IN):: M   !    total number of constraints (i.e. constraint equations)
  INTEGER,INTENT(IN):: N   !    number of control variables

  INTEGER:: i,j


!D7      EQUATIONS CITED FROM "A FAST ALGORITHM FOR NONLINEARLY
!D7      CONSTRAINED OPTIMIZATION CALCULATIONS" BY M.J.D. POWELL
! called by Vmaco

!-------------------------------------------------------------------------------
!write(*,*) "Entering Initvm, istat,n,m=", istat,n,m
      IBUGFG = 0
      ITER = 0
      NVMA = 1
      NVAR = N + 1
!
      IF (ISTAT .EQ. -1) THEN
!
!...    initialize the b matrix
!
        DO 20 I=1,N
          DO 10 J=1,N
            B(I,J)= ZERO
   10     CONTINUE
          B(I,I)= FP1
   20   CONTINUE
      ENDIF
!
      IF (M .GE. 1) THEN
        DO 30 I=1,M
          VMU(I) = ZERO
   30   CONTINUE
      ENDIF
!write(*,*) "Leaving InitVm, istat,n,m,nvar=", istat,n,m,nvar
  RETURN
END Subroutine InitVM   ! ------------------------------------------------------

!+
SUBROUTINE LGRANG()
! ------------------------------------------------------------------------------
! PURPOSE - Set up the Lagrangian function for quadratic programming iteration.

! Global variables used but not changed
!USE GlobalVmaco,ONLY: a,cstar,delta,gm,iactiv,ibasis,nc,nvar
! Global variables modified by Subroutine Lgrang
!USE GlobalVmaco,ONLY: amdamx,ghat,lamrow

IMPLICIT NONE
  REAL:: ambda !  lagrange multiplier (lambda) (= -cstar*ghat)
  INTEGER:: i
  REAL,PARAMETER:: XINF = 1E38 ! number for infinity

!D7   LIMITATIONS/ASSUMPTIONS
!D7
!D7      EQUATIONS AND ALGORITHM CITED FROM "A GENERAL QUADRATIC
!D7      ALGORITHM" AERE-TP.401 BY R. FLETCHER
!
! called by Quad
  INTRINSIC:: DOT_PRODUCT 
!-------------------------------------------------------------------------------
!...the gradient vector of the objective function is calculated

!      LSIZE = 1 + MAXX*(NVAR-1)
      DO I=1,NVAR
!!!        GHAT(I) = DOTPRD(NVAR,A(I:,1),LSIZE,MAXX,DELTA(1:),NVAR,1) - GM(I)   ! modified RLC
        ghat(i)=DOT_PRODUCT(a(i,1:nvar), delta(1:nvar) ) - gm(i)
      END DO
!
!...the maximum lagrange multiplier is calculated
!
      AMDAMX = -XINF
!
      DO 20 I=1,NC
        IF (iactiv(ibasis(i)) .NE. -1) THEN
!!!          AMBDA = -DOTPRD (NVAR,CSTAR(I:,1),LSIZE,MAXX,GHAT(1:),NVAR,1)   ! modified RLC
          ambda=-DOT_PRODUCT(cstar(i,1:nvar), ghat(1:nvar) )
          IF (ambda .GT. amdamx) THEN
            amdamx = ambda
            lamrow = i
          END IF
        END IF
   20 END DO

  RETURN
END Subroutine Lgrang   ! ------------------------------------------------------

!+
SUBROUTINE LNSRCH (c,f,g,istat,m,maxfun,meq,n,nvmitr,nvmlin,x,bflag,iprint)
! PURPOSE - Calculate the necessary information for the line search objective function.


! Global variables used but not changed
!USE GlobalVmaco,ONLY: delta,gdelta,nvma,vmu
! Global variables modified by Subroutine Lnsrch
!USE GlobalVmaco,ONLY: dflsav,flsa,oxa,xa

IMPLICIT NONE
  REAL,INTENT(IN),DIMENSION(:):: c   !     vector of constraint variables
  REAL,INTENT(IN):: f   !    value of the objective function (optimized variable)
  REAL,INTENT(IN),DIMENSION(:):: g   !    gradient of the objective function
  INTEGER,INTENT(IN OUT):: istat   ! status of the optimization variable
    ! =-1 initial pass through vmaco
    ! = 0 continue algorithm calculations
    ! = 1 converge within required accuracy
    ! > 1 error encountered
  INTEGER,INTENT(IN):: m   !    total number of constraints (i.e. constraint equations)
  INTEGER,INTENT(IN):: maxfun   ! maximum number of times vmaco can be called
  INTEGER,INTENT(IN):: meq   !  total number of equality constraints
  INTEGER,INTENT(IN):: n   !    number of control variables
  INTEGER,INTENT(IN):: nvmitr   ! value of nvma at the start of an iteration
  INTEGER,INTENT(IN):: nvmlin   ! controls the error return from the line search 
!                                  objective function. eqn(4.2)
  REAL,INTENT(IN OUT),DIMENSION(:):: x   !    vector of control variables
  LOGICAL,INTENT(IN OUT):: bflag   ! flag for whether the b matrix should be updated
  INTEGER,INTENT(IN):: iprint   ! flag to indicate amount of printout
!              = 0 no printout
!              = 1 values of x, c, and f are printed
!              = 2 debugging printout

  REAL,PARAMETER:: CHORD = 0.1 !  slope of the armijo chord
  REAL:: dflsa ! derivative of fls. corresponds to capital delta or
               !   phi'(0) in eqn (4.7)
  REAL:: diffls ! difference between the current line search objective
                ! function and the previous value. eqn (4.8)
  REAL:: fls !    value of the line search objective function. eqn (4.2)
  INTEGER:: i
  INTEGER:: ibug
  INTEGER:: idummy
  REAL:: stepl !  step-length for the search direction. corresponds to
               !   alpha in eqn (2.2)
  REAL,PARAMETER:: STPMAX = 0.1 ! limits the reduction in the line search step-length
  REAL:: suminf ! weighted sum of infeasibilities. eqn (4.2)
  REAL:: xdummy
!D7   LIMITATIONS/ASSUMPTIONS
!D7
!D7      EQUATIONS CITED FROM "A FAST ALGORITHM FOR NONLINEARLY
!D7      CONSTRAINED OPTIMIZATION CALCULATIONS" BY M.J.D. POWELL
! called by Vmaco
! calls Abig, BugAid, GradLf, VmaMsg

  INTRINSIC:: MAX
!-------------------------------------------------------------------------------
      suminf = ZERO

      IF (m .GT. 0) THEN
        DO i=1,m
          IF (i .LE. meq) THEN
            suminf = suminf + vmu(i)*ABS(c(i))
          ELSE
!!!            SUMINF = SUMINF + VMU(I)*ABIG(ZERO,-C(I))
            suminf=suminf+vmu(i)*MAX(zero, -c(i))
          END IF
        END DO
      ENDIF
!
      fls = f + suminf
!
      IF (nvma .EQ. nvmitr) THEN
!
!...    set the initial conditions for the line search
!
        flsa = fls
        dflsa = gdelta - delta(n+1)*suminf
!
        IF (dflsa .GE. ZERO) THEN
          istat = 7
          CALL Vmamsg(istat,0,ZERO)
          GO TO 999
        ELSE
          stepl = FP1
!
!...      calculate the new control variable guesses and save the old
!
          DO 20 i=1,n
            xa(i) = oxa(i)
            x(i) = xa(i) + delta(i)
   20     END DO

          dflsa = stepl*dflsa
          dflsav = dflsa
        END IF

        IF (nvma .GE. MAXFUN) THEN
          DO 30 i=1,n
            x(i) = xa(i)
   30     END DO

          istat = 5
          CALL Vmamsg (istat,maxfun,ZERO)
          GO TO 999
        END IF

        nvma = nvma + 1

      ELSE
!
!...    continue to calculate necessary information for the line search
!
        diffls = fls - flsa
        dflsa = dflsav
!
        IF (diffls .LE. CHORD*dflsa) THEN
          CALL Gradlf(g,m,n,iprint)
          bflag = .TRUE.
        ELSE
!
          IF (nvma .LT. nvmlin + nvmitr) THEN
!!!            STEPL = ABIG(STPMAX,HALF*DFLSA/(DFLSA-DIFFLS))
            stepl = MAX(stpmax, HALF*dflsa/(dflsa-diffls))
!
!...        calculate the new control variable guesses and save the old
!
            DO 40 I=1,N
              delta(i) = stepl*delta(i)
              xa(i) = oxa(i)
              x(i) = xa(i) + delta(i)
   40       END DO

            dflsa = stepl*dflsa
            dflsav = dflsa

            IF (nvma .LT. maxfun) THEN
              nvma = nvma + 1
              GO TO 999
            END IF
          END IF
!
          DO 50 i=1,n
            x(i) = xa(i)
   50     CONTINUE

          nvma = nvma + 1
        END IF
      END IF

  999 CONTINUE

      IF (iprint .EQ. 2) THEN
        IBUG = 3
        CALL Bugaid (ibug,idummy,fls)
        ibug = 1
        CALL Bugaid (ibug,idummy,dflsa)
        IF (.NOT. bflag) THEN
          ibug = 4
          CALL Bugaid (ibug,idummy,stepl)
          ibug = 5
!          write(*,*) 'calling BugAid in Lnsrch with n=', n
          CALL Bugaid (ibug,n,xdummy)
        END IF
      END IF

  RETURN
END Subroutine Lnsrch   ! ------------------------------------------------------

!+
SUBROUTINE MATIN(C,N,D,M,KEY,DETRMX)
! ------------------------------------------------------------------------------
! PURPOSE - Derive the inverse rank and determinant of an input matrix (c). 
! If desired, the inverted matrix can be multiplied by a vector (d) by turning 
! on the switch m=1. This version converts single precision inputs to double 
! precision accuracy internally. When matrix inversion is completed the
! outputs are converted back to single precision.


! Global variables used but not changed
!USE GlobalVmaco,ONLY: MAXX
IMPLICIT NONE
  REAL,INTENT(IN OUT),DIMENSION(:,:):: c   !    the matrix to be inverted
  REAL,INTENT(IN OUT),DIMENSION(:):: d   !    vector to multiply times the inverted matrix
  REAL,INTENT(IN OUT):: detrmx  ! determinant of the matrix
  INTEGER,INTENT(OUT):: key   !  rank of the inverted matrix
  INTEGER,INTENT(IN):: m   !    switch for the vector multiplication  
  INTEGER,INTENT(IN):: n   !    dimension of the matrix a and vector b
!!!  INTEGER,INTENT(IN):: nx  !    total size of the matrix coming in    

  INTEGER:: i,j,k
  INTEGER:: icolum,jcolum
  INTEGER:: index
  INTEGER:: ipivot
  INTEGER:: irow,jrow
  INTEGER:: l,l1
!
      DOUBLE PRECISION A
      DOUBLE PRECISION AMAX
      DOUBLE PRECISION B
      DOUBLE PRECISION DETERM
      DOUBLE PRECISION PIVOT
      DOUBLE PRECISION SWAP
      DOUBLE PRECISION T
      DOUBLE PRECISION ZERO_TOLERANCE   ! changed from ZERO by RLC 30Oct2005
!
!      PARAMETER ( MAXX   =    5 )
      INTEGER,PARAMETER:: MAXA   = MAXX * 2 
!
      DIMENSION A      (MAXA  ,MAXA  )
      DIMENSION B      (     1)
!!!      DIMENSION C      (NX    ,NX    )
!!!      DIMENSION D      (     1)
      DIMENSION IPIVOT (MAXA  )
      DIMENSION INDEX  (MAXA  ,     2)
      DIMENSION PIVOT  (MAXA  )
!-------------------------------------------------------------------------------

!...convert the c matrix and b vector to double precison (a and b)

      B(1) =  D(1)
!!!      DETERM =  DETRMX   ! removed by RLC  12Nov2008

      DO 15 I = 1,N
        DO 10 J=1,N
          A(I,J) =  C(I,J)
   10   CONTINUE
   15 END DO
!
!...initialization
!
      ZERO_TOLERANCE = 1.D-16   ! changed from ZERO by RLC 30Oct2005 not used???
      KEY = N
      DETERM = 1.D+0
      DO 20 J=1,N
        IPIVOT(J) = 0
   20 END DO
!
      DO 550 I=1,N
!
!...    search for pivot element
!
        AMAX = 0.D+0
!
        DO 105 J=1,N
          IF (IPIVOT(J) - 1 .NE. 0) THEN
            DO 100 K=1,N
              IF (IPIVOT(K) - 1 .LT. 0) THEN
                IF (DABS(AMAX) - DABS(A(J,K)) .LT. 0.) THEN
                  IROW = J
                  ICOLUM = K
                  AMAX = A(J,K)
                END IF
              ELSEIF (IPIVOT(J) - 1 .GT. 0) THEN
                GO TO 740
              ENDIF
  100       CONTINUE
          ENDIF
  105   CONTINUE
!
        IF (DABS(AMAX) .LE. ZERO_TOLERANCE ) THEN
          DETERM = 0.D0
          KEY = I-1
          GO TO 740
        ENDIF
        IPIVOT(ICOLUM) = IPIVOT(ICOLUM) + 1
!
!...    INTERCHANGE ROWS TO PUT PIVOT ELEMENT ON THE DIAGONAL
!
        IF (IROW - ICOLUM .NE. 0) THEN
          DETERM = -DETERM
          DO 200 L=1,N
            SWAP = A(IROW,L)
            A(IROW,L) = A(ICOLUM,L)
            A(ICOLUM,L) = SWAP
  200     CONTINUE
!
          IF (M .NE. 0) THEN
            SWAP = B(IROW)
            B(IROW) = B(ICOLUM)
            B(ICOLUM) = SWAP
          ENDIF
        ENDIF
!
        INDEX(I,1) = IROW
        INDEX(I,2) = ICOLUM
        PIVOT(I) = A(ICOLUM,ICOLUM)
!
!...    DIVIDE PIVOT ROW BY PIVOT ELEMENT
!
        A(ICOLUM,ICOLUM) = 1.D+0
!
        DO 350 L=1,N
          A(ICOLUM,L) = A(ICOLUM,L)/PIVOT(I)
  350   CONTINUE
!
        IF (M .NE. 0) THEN
          B(ICOLUM) = B(ICOLUM)/PIVOT(I)
        ENDIF
!
!...    REDUCE NON-PIVOT ROWS
!
        DO 540 L1=1,N
          IF (L1 - ICOLUM .NE. 0) THEN
            T = A(L1,ICOLUM)
            A(L1,ICOLUM) = 0.D0
            DO 450 L=1,N
              A(L1,L) = A(L1,L) - A(ICOLUM,L)*T
  450       CONTINUE
!
            IF (M .NE. 0) THEN
              B(L1) = B(L1) - B(ICOLUM)*T
            ENDIF
          ENDIF
  540   CONTINUE
  550 END DO
!
!...interchange columns
!
      DO 710 I=1,N
        L= N+ 1 - I
        IF (INDEX(L,1) - INDEX(L,2) .NE. 0) THEN
          JROW = INDEX (L,1)
          JCOLUM = INDEX(L,2)
          DO 705 K=1,N
            SWAP = A(K,JROW)
            A(K,JROW) = A(K,JCOLUM)
            A(K,JCOLUM) = SWAP
  705     CONTINUE
        ENDIF
  710 END DO
!
      DO 720 I=1,N
        J = N + 1 - I
        DETERM = DETERM*PIVOT(J)
  720 END DO
!
!...convert back to single precision accuracy
!
  740 CONTINUE
      D(1) = (B(1))
      DETRMX = (DETERM)
      DO 770 I=1,N
        DO 760 J=1,N
          C(I,J) = (A(I,J))
  760   CONTINUE
  770 END DO

  RETURN
END Subroutine Matin   ! -------------------------------------------------------

!+
SUBROUTINE Move()
! ------------------------------------------------------------------------------
! PURPOSE - Step along the search direction and distance calculated to get a new 
!  vertex. Also calculates necessary information for adding or exchanging a 
!  constraint in the basis.


! Global variables used but not changed
!USE GlobalVmaco,ONLY: amdamx,cac,cn,cnear,cstar,cstep,h,icnear,nc,nvar,nvar2

! Global variables modified by Subroutine Move
!USE GlobalVmaco,ONLY: ccgcac,chc,cstarc,delta,hdotc,iactiv

IMPLICIT NONE
  INTEGER:: i,j
!D7      EQUATIONS AND ALGORITHM CITED FROM "A GENERAL QUADRATIC
!D7      ALGORITHM" AERE-TP.401 BY R. FLETCHER

! called by Serch

  INTRINSIC:: DOT_PRODUCT
!-------------------------------------------------------------------------------
      DO 10 I=1,NVAR
        DELTA(I) = DELTA(I) + CNEAR*CSTEP(I)
   10 END DO
!
      CCGCAC = CNEAR*AMDAMX/CAC
!
      IF (ICNEAR .LE. NVAR) THEN
        IF (NVAR .NE. NC) THEN
          DO 20 I=1,NVAR
            HDOTC(I) = H(I,ICNEAR)
   20     CONTINUE
        ENDIF
!
        DO 30 I=1,NC
          CSTARC(I) = CSTAR(I,ICNEAR)
   30   CONTINUE
        CHC = HDOTC(ICNEAR)
!
      ELSEIF (ICNEAR .LE. NVAR2) THEN
        J = ICNEAR - NVAR
        IF (NVAR .NE. NC) THEN
          DO 40 I=1,NVAR
            HDOTC(I) = -H(I,J)
   40     CONTINUE
        ENDIF
!
        DO 50 I=1,NC
          CSTARC(I) = -CSTAR(I,J)
   50   CONTINUE
        CHC = -HDOTC(J)
!
      ELSE
        J = ICNEAR - NVAR2
!        LSIZE = 1 + MAXX*(NVAR - 1)
        IF (NVAR .NE. NC) THEN
          DO I=1,NVAR
!!!            HDOTC(I) = DOTPRD (NVAR, H(I:,1),LSIZE,MAXX,CN(1:,J),NVAR,1) ! modified RLC
            hdotc(i)=DOT_PRODUCT(h(i,1:nvar), cn(1:nvar,j) )
          END DO
!!!          CHC = DOTPRD (NVAR,CN(1:,J),NVAR,1,HDOTC(1:),NVAR,1)
          chc=DOT_PRODUCT(cn(1:nvar,j), hdotc(1:nvar) )
        ENDIF
!
        DO I=1,NC
!!!          CSTARC(I) = DOTPRD (NVAR,CSTAR(I:,1),LSIZE,MAXX,CN(1:,J),NVAR,1) ! modified RLC
          cstarc(i)=DOT_PRODUCT(cstar(i,1:nvar), cn(1:nvar,j) )
        END DO
      END IF
!
      IACTIV(ICNEAR) = 0

  RETURN
END Subroutine Move   ! --------------------------------------------------------

!+
SUBROUTINE NEWBAS (ICROW)
! ------------------------------------------------------------------------------
! PURPOSE - A new basis is obtained by exchanging constraints.

! Global variables used but not changed
!USE GlobalVmaco,ONLY: cinvrs,cn,MAXX,mtotal,nvar,nvar2,resid

! Global variables modified by Subroutine NewBas
!USE GlobalVmaco,ONLY: delta,iactiv,ibasis

IMPLICIT NONE
  INTEGER,INTENT(IN):: ICROW   ! index of the row of cinvrs which is to be 
                               ! removed from the basis.
  REAL:: alpha !  nearest nonviolated nonbasic constraint distance
  REAL:: beta  !  furthest violated constraint distance and then the
!d6              minimum between alpha and beta
  REAL,DIMENSION(MAXX):: crow !   constraint vector of cinvrs which is being removed
!d6              from the basis
  REAL:: distnc ! distance moved away from the vertex in the direction
!d6              of ui

  INTEGER:: i,j
  INTEGER:: iaindx ! index of the constraint corresponding to alpha
  INTEGER:: ibindx ! index of the constraint corresponding to beta
  REAL:: recipr ! reciprical of the residual constraint to be removed
!d6              from the basis
  REAL:: rrecip ! product of the ith residual in the basis and recipr
  REAL,DIMENSION(MAXX):: ui !     ith row of cinvrs. direction of search from the
!d6              current vertex
  REAL:: uicj !   product of the ith row of cinvrs and the normal vectors
!d6              of constraints in the basis
  REAL,PARAMETER:: XINF = 1E38 ! number for infinity
!D7
!D7      EQUATIONS AND ALGORITHM CITED FROM "THE CALCULATION OF FEASIBLE
!D7      POINTS FOR LINEARLY CONSTRAINED OPTIMIZATION PROBLEMS" AERE-
!D7      R.6354 BY R. FLETCHER

! called by Vertex
! calls DotPrd
  INTRINSIC:: DOT_PRODUCT
!-------------------------------------------------------------------------------
!...move away from the vertex in the direction of ui

      DO 10 I=1,NVAR
        UI(I) = CINVRS(ICROW,I)
   10 END DO

      ALPHA = XINF
      BETA = ZERO

      DO 20 I=1,MTOTAL
        IF (IACTIV(I) .GE. 1) THEN
          IF (I .LE. NVAR) THEN
            UICJ = -UI(I)
          ELSEIF (I .LE. NVAR2) THEN
            J = I - NVAR
            UICJ = UI(J)
          ELSE
            J = I - NVAR2
!!!            UICJ = -DOTPRD (NVAR,UI(1:),NVAR,1,CN(1:,J),NVAR,1)  ! modified RLC
            uicj=-DOT_PRODUCT(ui(1:nvar), cn(1:nvar,j) )
          END IF

          IF (IACTIV(I) .EQ. 1) THEN

!...        find the nearest nonviolated nonbasis constraint
            IF (UICJ .GT. 0.) THEN
              DISTNC = RESID(I)/UICJ
              IF (DISTNC .LT. ALPHA) THEN
                ALPHA = DISTNC
                IAINDX = I
              END IF
            END IF
          ELSE

!...        find the furthest violated constraint
            IACTIV(I) = 1
            IF (UICJ .LT. 0.) THEN
              DISTNC = RESID(I)/UICJ
              IF (DISTNC .GT. BETA) THEN
                BETA = DISTNC
                IBINDX = I
              END IF
            END IF
          END IF
        END IF
   20 END DO

      IF (ALPHA .LE. BETA) THEN

!...    new vertex at x + alpha*ui
        IBINDX = IAINDX
        BETA = ALPHA
      ENDIF

!...exchange the constraints within the basis
      IACTIV(IBASIS(ICROW)) = 1
      IACTIV(IBINDX) = 0
      IBASIS(ICROW) = IBINDX

!...utilize the simplex method to exchange within cinvrs
      IF (IBINDX .LE. NVAR) THEN
        DO 30 I=1,NVAR
          CROW(I) = CINVRS(I,IBINDX)
   30   CONTINUE
      ELSEIF (IBINDX .LE. NVAR2) THEN
        IBINDX = IBINDX - NVAR
        DO 40 I=1,NVAR
          CROW(I) = -CINVRS(I,IBINDX)
   40   CONTINUE
      ELSE
        IBINDX = IBINDX - NVAR2
!!!        LSIZE = 1 + MAXX*(NVAR - 1)
        DO I=1,NVAR
!!!          CROW(I) = DOTPRD (NVAR,CINVRS(I:,1),LSIZE,MAXX,CN(1:,IBINDX),   &  ! modified RLC
!!!     &                      NVAR,1)
          crow(i)=DOT_PRODUCT(cinvrs(i,1:nvar), cn(1:nvar,ibindx) )
        END DO
      END IF
!
      recipr = FP1/crow(icrow)
      DO 80 i=1,nvar
        delta(i) = delta(i) + beta*ui(i)
        IF (i .EQ. icrow) THEN
          DO 60 j=1,nvar
            cinvrs(i,j) = cinvrs(i,j)*recipr
   60     END DO
        ELSE
          rrecip = recipr*crow(i)
          DO 70 j=1,nvar
            cinvrs(i,j) = cinvrs(i,j) - rrecip*ui(j)
   70     END DO
        END IF
   80 END DO

  RETURN
END Subroutine NewBas   ! ------------------------------------------------------

!+
SUBROUTINE NEWDIR()
! ------------------------------------------------------------------------------
! PURPOSE - Calculate the new maximum Lagrange multiplier, gradient direction, 
!  and search direction.

! Global variables used but not changed
!USE GlobalVmaco,ONLY: a,cstep,delta,gm,h,nvar
! Global variables modified by Subroutine NewDir
!USE GlobalVmaco,ONLY: amdamx,ghat,passiv

IMPLICIT NONE
  INTEGER:: i
  REAL,PARAMETER:: TOL = 1E-38 !   tolerance for equality of zero

!D7   LIMITATIONS/ASSUMPTIONS
!D7
!D7      EQUATIONS AND ALGORITHM CITED FROM "A GENERAL QUADRATIC
!D7      ALGORITHM" AERE-TP.401 BY R. FLETCHER
! called by Serch

  INTRINSIC:: DOT_PRODUCT
!-------------------------------------------------------------------------------

!      LSIZE = 1 + MAXX*(NVAR - 1)
      DO I=1,NVAR
!!!        GHAT(I) = DOTPRD(NVAR,A(I:,1),LSIZE,MAXX,DELTA(1:),NVAR,1) - GM(I) ! modified RLC
        ghat(i)=DOT_PRODUCT(a(i,1:nvar), delta(1:nvar) )
      END DO

      AMDAMX = ZERO

      DO I=1,NVAR
!!!        CSTEP(I) = -DOTPRD (NVAR,H(I:,1),LSIZE,MAXX,GHAT(I:),NVAR,1) ! modified RLC
        cstep(i)=-DOT_PRODUCT(h(i,1:nvar), ghat(1:nvar) )
        IF (ABS(CSTEP(I)) .GT. TOL) THEN
          AMDAMX = FP1
        ENDIF
      END DO

      PASSIV = .FALSE.
  RETURN
END Subroutine NewDir   ! ------------------------------------------------------

!+
SUBROUTINE NODESG()
! ------------------------------------------------------------------------------
! PURPOSE - No designated constraints, vertex chosen from the lower and upper 
! bounds. (whichever is closer)

! Global variables used but not changed
!USE GlobalVmaco,ONLY: bdl,bdu,delta,nvar
! Global variables modified by Subroutine NoDesg
!USE GlobalVmaco,ONLY: cinvrs,iactiv,ibasis,nc

IMPLICIT NONE
  INTEGER:: i,j

!D7   LIMITATIONS/ASSUMPTIONS
!D7
!D7      EQUATIONS AND ALGORITHM CITED FROM "THE CALCULATION OF FEASIBLE
!D7      POINTS FOR LINEARLY CONSTRAINED OPTIMIZATION PROBLEMS" AERE-
!D7      R.6354 BY R. FLETCHER
! called by Vertex
!-------------------------------------------------------------------------------
      DO 20 I=1,NVAR
        DO 10 J=1,NVAR
          CINVRS(I,J) = ZERO
   10   CONTINUE

        IF (DELTA(I) - BDL(I) .LE. BDU(I) - DELTA(I)) THEN
          IBASIS(I) = I
          CINVRS(I,I) = FP1
        ELSE
          IBASIS(I) = NVAR + I
          CINVRS(I,I) = -FP1
        END IF

        J = IBASIS(I)
        IACTIV(J) = 0
   20 END DO

      NC = NVAR

  RETURN
END Subroutine NoDesg   ! ------------------------------------------------------

!+
SUBROUTINE Quad (istat,meq,iprint)
! ------------------------------------------------------------------------------
! PURPOSE - Minimize a quadratic function of n variables subject to m linear 
!  equality and inequality constraints.

! Global variables used but not changed
!USE GlobalVmaco,ONLY: amdamx,MAXX,nc,nvar
! Global variables modified by Subroutine Quad
!USE GlobalVmaco,ONLY: cstar,h

IMPLICIT NONE
  INTEGER,INTENT(IN):: iprint  ! flag to indicate amount of printout
                               ! = 0 no printout
                               ! = 1 values of x, c, and f are printed
                               ! = 2 debugging printout
  INTEGER,INTENT(IN OUT):: istat   ! status of the optimization algorithm
    !  =-2 initial pass through vmaco with the b matrix already initialized
    !  =-1 initial pass through vmaco
    !  = 0 continue algorithm calculations
    !  = 1 converge within required accuracy
    !  > 1 error encountered
  INTEGER,INTENT(IN):: meq   !  total number of equality constraints
  INTEGER,PARAMETER:: MAXHC  = MAXX*2 !   maximum size of square array dealing with hcmat
  REAL,DIMENSION(MAXHC,MAXHC):: hcmat ! general matrix containing equations
    ! governing equality problem eqn (7). Partitions of this matrix are the 
    ! h and cstar operators.
  INTEGER:: i,j, ix
  INTEGER:: ibug
  INTEGER:: idummy
  INTEGER:: nserch ! counter on the iteration number curently on in the search
                   ! for the minimum of the quadratic function
  INTEGER,PARAMETER:: NSCHMX = 100 ! maximum value allowed for nserch
  REAL:: residu
  LOGICAL:: retest ! flag to indicate whether the final solution has been verified
  REAL,PARAMETER:: TOL = 1E-38 !    tolerance for equality of zero
  REAL:: xdummy


!D7      EQUATIONS AND ALGORITHM CITED FROM "A GENERAL QUADRATIC
!D7      ALGORITHM" AERE-TP.401 BY R. FLETCHER
!-------------------------------------------------------------------------------
!  write(*,*) 'Entering Quad with nvar=', nvar
      RETEST = .FALSE.
!

   10 CONTINUE
      CALL Vertex(meq) ! call the feasible vertex routine

      IF (ISTAT .NE. 0) THEN
        NC = 0
        GO TO 999
      END IF

      NC = NVAR
      DO I=1,nvar  ! initialize the operator h=0.
        DO J=1,nvar
          h(i,j) = ZERO
        END DO
      END DO

      NSERCH = 0
   40 CONTINUE
      NSERCH = NSERCH + 1

!...protect against an infinite loop occurring
      IF (nserch .GE. NSCHMX) THEN
        ISTAT = 14
        CALL Vmamsg(istat,NSCHMX,XDUMMY)
        IBUG = 5
        CALL Bugaid(ibug,nvar,xdummy)
        IBUG = 7
        CALL Bugaid(ibug,nvar,xdummy)
        IBUG = 6
        CALL Bugaid(ibug,idummy,amdamx)
        GO TO 999
      END IF

      CALL LGRANG()   !...initialize the Lagrangian function

      IF (iprint .EQ. 2) THEN
        IBUG = 5
        CALL Bugaid(ibug,nvar,xdummy)
        IBUG = 7
        CALL Bugaid(ibug,nvar,xdummy)
        IBUG = 6
        CALL Bugaid(ibug,idummy,amdamx)
      END IF
!
  IF (amdamx .GT. ZERO) THEN
    CALL SERCH()  ! downhill progress can be made by dropping a constraint from the basis
    GO TO 40
  END IF

  IF (retest) THEN   !...    a verified solution has been found
    IF (ABS(amdamx) .LT. TOL) THEN
      CALL Vmamsg(10,0,ZERO)   !...      solution may be a degenerate local minimum
    END IF
  ELSE   
    retest = .TRUE.   ! retest the Lagrange function because the successive updating
                      ! of the operators (cstar and h) might cause an accumulation error
    CALL RESET (HCMAT,RESIDU)

    IF (residu .LT. ZERO) THEN
      GO TO 10
    ELSE
      DO i=1,nvar
        ix = i+nvar
        DO j=1,nvar
          h(i,j)=hcmat(i,j)
          cstar(i,j)=hcmat(ix,j)
        END DO
      END DO
      GO TO 40
    END IF
  END IF

  999 CONTINUE

  RETURN
END Subroutine Quad   ! --------------------------------------------------------

!+
SUBROUTINE REMOVE()
! PURPOSE - Constraint is removed from the basis and the corresponding
!  h and cstar operators are updated.


! Global variables used but not changed
!USE GlobalVmaco,ONLY: ac,cac,ct,lamrow,nc,nvar
! Global variables modified by Subroutine Remove
!USE GlobalVmaco,ONLY: ckac,cstar,h,ibasis

IMPLICIT NONE
  REAL:: ckacdc ! value of ckac divided by cac
  REAL:: ctdcac ! value of ct divided by cac
  INTEGER:: i,j
!D7      EQUATIONS AND ALGORITHM CITED FROM "A GENERAL QUADRATIC
!D7      ALGORITHM" AERE-TP.401 BY R. FLETCHER
!
  INTRINSIC:: DOT_PRODUCT
!-------------------------------------------------------------------------------
!...first the h operator is updated
      DO 20 I=1,NVAR
        CTDCAC = CT(I)/CAC
        DO 10 J=1,I
          H(I,J) = H(I,J) + CTDCAC*CT(J)
          H(J,I) = H(I,J)
   10   CONTINUE
   20 END DO

      IF (NC .GT. 1) THEN
        IF (LAMROW .NE. NC) THEN
!
!...      replace the delegated constraint by moving the last constraint
!...      in the basis into its place
          DO 30 I=1,NVAR
            CSTAR(LAMROW,I) = CSTAR(NC,I)
   30     CONTINUE
          IBASIS(LAMROW) = IBASIS(NC)
        END IF
!
!...    update the cstar operator matrix
        NC = NC - 1

!!!        LSIZE = 1 + MAXX*(NVAR - 1)
        DO 40 I=1,NC
!!!          CKAC(I) = DOTPRD (NVAR,CSTAR(I:,1),LSIZE,MAXX,AC(1:),NVAR,1)
          ckac(i)=DOT_PRODUCT(cstar(i,1:nvar), ac(1:nvar) )
   40   CONTINUE

        DO 60 I=1,NC
          CKACDC = CKAC(I)/CAC
          DO 50 J=1,NVAR
            CSTAR(I,J) =CSTAR(I,J) - CKACDC*CT(J)
   50     CONTINUE
   60   CONTINUE

      ELSE
        NC = 0
      ENDIF

  RETURN
END Subroutine Remove   ! ------------------------------------------------------

!+
SUBROUTINE RESET(HCMAT,RESIDU)
! ------------------------------------------------------------------------------
! PURPOSE - Set up the matrix and right hand side of equations governing the 
!  quadratic (and Lagrangian function) equality problem.


! Global variables used but not changed
!USE GlobalVmaco,ONLY: a,bdl,bdu,cn,d,delta,gm,iactiv,ibasis,MAXX,mtotal,nc,nvar,nvar2
! No global variables modified by Subroutine Reset

IMPLICIT NONE
  REAL,INTENT(OUT),DIMENSION(:,:):: HCMAT   ! general matrix containing 
!             equations governing equality problem. See eqn (7). Partitions of 
!             this matrix are h and cstar operators.
  REAL,INTENT(OUT):: RESIDU  ! constraint residual for the quadratic equality problem

  REAL:: determ
  REAL:: dummy
  INTEGER:: i,j, ix,jx
  INTEGER:: key ! rank of the inverted matrix
  INTEGER,PARAMETER:: MAXHC = MAXX*2 ! maximum size of arrays dealing with the hcmat matrix
  INTEGER:: nnc !    sum of nvar + nc
  REAL,DIMENSION(MAXHC):: rhsep ! right hand side of the euation governing the equality problem.
!D7   LIMITATIONS/ASSUMPTIONS
!D7
!D7      EQUATIONS AND ALGORITHM CITED FROM "A GENERAL QUADRATIC
!D7      ALGORITHM" AERE-TP.401 BY R. FLETCHER
! called by Quad
! calls Matin, VmaMsg and DotPrd
  REAL,DIMENSION(1):: vecdum    ! added by RLC, 30 Oct 2005

  INTRINSIC:: DOT_PRODUCT
!-------------------------------------------------------------------------------

  hcmat=ZERO

      DO 20 I=1,NVAR
        RHSEP(I) = GM(I)
        DO 10 J=1,NVAR
          HCMAT(I,J) = A(I,J)
   10   CONTINUE
   20 END DO
!
      IF (NC .NE. 0) THEN
        DO 60 I=1,NC
          IF (IBASIS(I) .LE. NVAR2) THEN
            DO 30 J=1,NVAR
              HCMAT(J,NVAR+I) = ZERO
              HCMAT(NVAR+I,J) = ZERO
   30       END DO
!
            IF (IBASIS(I) .LE. NVAR) THEN
              HCMAT(NVAR+I,IBASIS(I)) = FP1
              HCMAT(IBASIS(I),NVAR+I) = FP1
              RHSEP(NVAR+I) = BDL(IBASIS(I))
            ELSE
              J = IBASIS(I) - NVAR
              HCMAT(NVAR+I,J) = -FP1
              HCMAT(J,NVAR+I) = -FP1
              RHSEP(NVAR+I) = -BDU(J)
            ENDIF
          ELSE
            J = IBASIS(I) - NVAR2
            DO 40 JX = 1,NVAR
              HCMAT(NVAR+I,JX) = CN(JX,J)
              HCMAT(JX,NVAR+I) = CN(JX,J)
   40       CONTINUE
            RHSEP(NVAR+I) = D(J)
          ENDIF
!
          DO 50 J=1,NC
            IX = NVAR + I
            JX = NVAR + J
            HCMAT(IX,JX) = ZERO
   50     CONTINUE
   60   CONTINUE
      ENDIF
!
      NNC = NVAR + NC

      CALL MATIN (HCMAT,NNC,VECDUM,0,KEY,DETERM)   !...invert the general matrix
      IF (KEY .LT. NNC) THEN
        CALL VMAMSG (15,2,DUMMY)
      END IF

      DO I=1,NVAR
!!!        DELTA(I) = DOTPRD (NNC,HCMAT(1:,I),NNC,1,RHSEP(1:),NNC,1) ! modified RLC
        delta(i)=DOT_PRODUCT(hcmat(1:nnc,i), rhsep(1:nnc) )
      END DO
!
!...check the feasibility of the new point
!
      DO I=1,MTOTAL
        IF (IACTIV(I) .GE. 1 ) THEN
          IF (I .LE. NVAR) THEN
            RESIDU = DELTA(I) - BDL(I)
          ELSEIF (I .LE. NVAR2) THEN
            J = I-NVAR
            RESIDU = BDU(J) - DELTA(J)
          ELSE
            J = I-NVAR2
!!!            RESIDU = DOTPRD (NVAR,CN(1:,J),NVAR,1,DELTA(1:),NVAR,1) - D(J) ! modified RLC
            residu=DOT_PRODUCT(cn(1:nvar,j), delta(1:nvar) )
          END IF
        END IF
      END DO

  RETURN
END Subroutine Reset   ! -------------------------------------------------------

!+
SUBROUTINE SERCH()
! ------------------------------------------------------------------------------
! PURPOSE - Search for an optimal solution by dropping a constraint from the basis

! Global variables used but not changed
!USE GlobalVmaco,ONLY: amdamx,cnear,cstar,cstep,nc,nvar,passiv
! Global variables modified by Subroutine Serch
!USE GlobalVmaco,ONLY: ct,delta,lamrow

IMPLICIT NONE
  LOGICAL:: depend
  LOGICAL:: postiv
  INTEGER:: i
  REAL,PARAMETER:: TOL = 1E-38 !    tolerance equality for zero
!D7
!D7      EQUATIONS AND ALGORITHM CITED FROM "A GENERAL QUADRATIC
!D7      ALGORITHM" AERE-TP.401 BY R. FLETCHER
! called by Quad
! calls Add,Changp,Chkdep,Defer,Exchng,Move,NewDir,Remove,Slnear,StepC
!-------------------------------------------------------------------------------

!...set the direction of search as corresponding to the row of cstar
!...with the largest positive lagrange multiplier
!
      DO 10 I=1,NVAR
        CT(I) = CSTAR(LAMROW,I)
   10 END DO
!
   20 CONTINUE
      CALL STEPC (POSTIV)
!
   30 CONTINUE
      CALL SLNEAR
!
      IF (POSTIV .AND. CNEAR .GE. FP1) THEN
!
!...    a minimum has been found
!
        DO 40 I=1,NVAR
          DELTA(I) = DELTA(I) + CSTEP(I)
   40   CONTINUE
        IF (PASSIV) THEN
          CALL REMOVE
        ENDIF
        GO TO 999
      ENDIF
!
!...move to a new vertex and calculate pertinent information
!
      CALL MOVE
!
      IF (NC .EQ. NVAR) THEN
!
!...    apply the simplex formula to exchange constraints
!
        CALL EXCHNG
        GO TO 999
      ENDIF
!
      DEPEND = .TRUE.
      IF (PASSIV) THEN
!
!...    check the dependency of the new constraint
!
        CALL CHKDEP (DEPEND)
      ENDIF
!
      IF (DEPEND) THEN
!
!...    apply formula for exchanging new constraint with a passive one
!
        CALL CHANGP
      ELSE
!
!...    apply formula for adding a constraint to the basis
!
        CALL ADD
        IF (PASSIV) THEN
!
!...      removal of a constraint has been deferred. set up as if the
!...      constraint is being removed from the augmented basis
!
          CALL DEFER
!
          IF (ABS(AMDAMX) .LT. TOL) THEN
            CALL REMOVE
            GO TO 999
          ELSE
            GO TO 20
          ENDIF
        ENDIF
      ENDIF
!
      IF (NC .NE. NVAR) THEN
!
!...    calculate the new search direction
!
        CALL NEWDIR
        IF (ABS(AMDAMX) .GT. TOL) THEN
          POSTIV = .TRUE.
          GO TO 30
        ENDIF
      ENDIF
!
  999 CONTINUE

  RETURN
END Subroutine Serch   ! -------------------------------------------------------

!+
SUBROUTINE SLNEAR()
! ------------------------------------------------------------------------------
! PURPOSE - A linear search along the direction of search is conducted.


! Global variables used but not changed
!USE GlobalVmaco,ONLY: bdl,bdu,cn,cstep,d,delta,ibasis,lamrow,mtotal,nvar,nvar2,passiv
! Global variables modified by Subroutine Slnear
!USE GlobalVmaco,ONLY: cnear,iactiv,icnear

IMPLICIT NONE
  REAL:: ccstep
  REAL:: distc
  INTEGER:: i,j
  REAL,PARAMETER:: XINF =1E38 !   number for infinity
!D7   LIMITATIONS/ASSUMPTIONS
!D7
!D7      EQUATIONS AND ALGORITHM CITED FROM "A GENERAL QUADRATIC
!D7      ALGORITHM" AERE-TP.401 BY R. FLETCHER
! called by Serch

  INTRINSIC:: DOT_PRODUCT
!-------------------------------------------------------------------------------

      CNEAR = XINF

      DO 10 I=1,MTOTAL
        IF (IACTIV(I) .EQ. 1) THEN
!
!...      the constraint being considered is inactive
!
          IF (I .LE. NVAR) THEN
            IF (CSTEP(I) .LT. ZERO) THEN
              DISTC = (BDL(I)-DELTA(I))/CSTEP(I)
              IF (DISTC .LT. CNEAR) THEN
                CNEAR = DISTC
                ICNEAR = I
              ENDIF
            ENDIF
          ELSEIF (I .LE. NVAR2) THEN
            J = I - NVAR
            IF (CSTEP(J) .GT. ZERO) THEN
              DISTC = (BDU(J) - DELTA(J))/CSTEP(J)
              IF (DISTC .LT. CNEAR) THEN
                CNEAR = DISTC
                ICNEAR = I
              ENDIF
            ENDIF
          ELSE
            J = I - NVAR2
!!!            CCSTEP = DOTPRD (NVAR,CN(1:,J),NVAR,1,CSTEP(1:),NVAR,1) ! modified RLC
            ccstep=DOT_PRODUCT(cn(1:nvar,j), cstep(1:nvar) )
            IF (CCSTEP .LT. ZERO) THEN
!!!              DISTC = D(J) -DOTPRD (NVAR,CN(1:,J),NVAR,1,DELTA(1:),NVAR,1) ! modified RLC
              distc=d(j)-DOT_PRODUCT(cn(1:nvar,j), delta(1:nvar) )
              DISTC = DISTC/CCSTEP
              IF (DISTC .LT. CNEAR) THEN
                CNEAR = DISTC
                ICNEAR = I
              END IF
            END IF
          END IF
        END IF
   10 END DO
!
!...the current constraint (direction of search) is removed from the
!...basis
!
      IF (PASSIV) THEN
        IACTIV(IBASIS(LAMROW)) = 1
      END IF

  RETURN
END Subroutine Slnear   ! ------------------------------------------------------

!+
SUBROUTINE Stepc(postiv)
! ------------------------------------------------------------------------------
! PURPOSE - Calculate the step to the stationary point in the new basis.

! Global variables used but not changed
!USE GlobalVmaco,ONLY: a,amdamx,ct,nvar
! Global variables modified by Subroutine Stepc
!USE GlobalVmaco,ONLY: ac,cac,cstep,passiv

IMPLICIT NONE
  LOGICAL,INTENT(OUT):: postiv  ! curvature of the constraint vector being 
                                ! removed from the basis.
                                !  = .true. positive curvature (optimal)
                                !  = .false. negative curvature

  INTEGER:: i
!D7      EQUATIONS AND ALGORITH CITED FROM "A GENERAL QUADRATIC
!D7      ALGORITHM" AERE-TP.401 BY R. FLETCHER
! called by Serch
  INTRINSIC:: DOT_PRODUCT
!-------------------------------------------------------------------------------
  DO I=1,NVAR
    ac(i)=DOT_PRODUCT(a(i,1:nvar), ct(1:nvar) )
  END DO
  cac=DOT_PRODUCT(ct(1:nvar), ac(1:nvar) )

  IF (cac .LE. ZERO) THEN
!
!...    curvature of the constraint vector being remove from the basis
!...    is negative and therefore the step to a stationary point is
!...    not a minimum
!
        postiv = .false.
        DO 20 i=1,nvar
          cstep(i) = ct(i)
   20   END DO
      ELSE
!
!...    curvature of constraint vector being removed from the basis
!...    is positive and the minimum can be found at cstep
!
        POSTIV = .TRUE.
        DO 30 i=1,nvar
          cstep(i) = ct(i)*amdamx/cac
   30   END DO
      ENDIF

      passiv = .TRUE.

  RETURN
END Subroutine StepC   ! -------------------------------------------------------

!+
SUBROUTINE Stinit (istat,m,meq,n,vlarge)
! ------------------------------------------------------------------------------
! PURPOSE - Currently on the initial call to vmaco (istat=-1) and therefore 
!  variables and elements are set up accordingly for the 
!  quadratic programming function

! Global variables used but not changed
!USE GlobalVmaco,ONLY: nvar
! Global variables modified by Subroutine Stinit
!USE GlobalVmaco,ONLY: b=>a,bdl,bdu,delta,gm,ibasis,mtotal,nc,nvar2

IMPLICIT NONE
  INTEGER,INTENT(IN OUT):: istat   ! status of the optimization algorithm
    ! =-2 initial pass through vmaco with the b matrix already initialized
    ! =-1 initial pass through vmaco
    ! = 0 continue algorithm calculations
    ! = 1 converge within required accuracy
    ! > 1 error encountered
  INTEGER,INTENT(IN):: m   ! total number of constraints (i.e. constraint equations)
  INTEGER,INTENT(IN):: meq ! total number of equality constraints                   
  INTEGER,INTENT(IN):: n   ! number of control variables                            
  REAL,INTENT(IN):: vlarge ! a very large number 

  INTEGER:: i

!D7      EQUATIONS AND ALGORITHM CITED FROM "A GENERAL QUADRATIC
!D7      ALGORITHM" AERE-TP.401 BY R. FLETCHER
! called by VmaSet
!-------------------------------------------------------------------------------
!...  set initial variables and elements
      nvar2 = nvar*2
      mtotal = nvar2 + m
      nc = meq + 1

      DO 10 i=1,n
        bdl(i) = -vlarge
        bdu(i) = vlarge
        delta(i) = ZERO
   10 END DO

      bdl(nvar) = ZERO
      delta(nvar) = FP1

      IF (meq .GE. 1) THEN
        DO 20 i=1,meq
          ibasis(i) = i + nvar2
   20   END DO
      END IF

      ibasis(nc) = nvar2

!...set the extra elements for gm and b to allow for infeasibility
      gm(nvar) = vlarge
      IF (istat .EQ. -1) THEN
        DO 30 I=1,nvar
          b(i,nvar) = ZERO
          b(nvar,i) = ZERO
   30   END DO
      ENDIF

      istat = 0

  RETURN
END Subroutine StInit   ! ------------------------------------------------------

!+
SUBROUTINE STOVMU (G,M,N,X)
! ------------------------------------------------------------------------------
! PURPOSE - store the necessary information for the gradients of the Lagrange 
!  function and also the multipliers for the line search function.


! Global variables used but not changed
!USE GlobalVmaco,ONLY: delta,glf,vlam
! Global variables modified by Subroutine Stovmu
!USE GlobalVmaco,ONLY: gdelta,glfa,oxa,vmu,xa

IMPLICIT NONE
  REAL,INTENT(IN),DIMENSION(:):: g   !    gradient of the objective function
  INTEGER,INTENT(IN):: m   !    total number of constraints (i.e. constraint equations)
  INTEGER,INTENT(IN):: n   !    number of control variables
  REAL,INTENT(IN),DIMENSION(:):: x   !    vector of control variables

  REAL:: abslam ! absolute value of Lagrange multipliers (vlam)
  INTEGER:: i
!D7      EQUATIONS CITED FROM "A FAST ALGORITHM FOR NONLINEARLY
!D7      CONSTRAINED OPTIMIZATION CALCULATIONS" BY M.J.D. POWELL
! called by Vmaco

  INTRINSIC:: MAX
!-------------------------------------------------------------------------------
!
!...first the scalar product of g and delta is calculated. The control
!...value guesses and the gradient are also saved

      gdelta = ZERO
      DO 10 i=1,n
        gdelta = gdelta + g(i)*delta(i)
        glfa(i) = glf(i)
        xa(i) = x(i)
        oxa(i) = xa(i)
   10 END DO

      IF (m .GT. 0) THEN  ! revise the vmu vector
        DO i=1,m
          abslam = ABS(vlam(i))
!!!          VMU(I) = ABIG(ABSLAM,HALF*(ABSLAM + VMU(I)))
          vmu(i)=MAX(abslam, HALF*(abslam + vmu(i)))
        END DO
      END IF

  RETURN
END Subroutine STOVMU   ! ------------------------------------------------------

!+
SUBROUTINE UpdatB(N)
! ------------------------------------------------------------------------------
! PURPOSE - Update the variable metric matrix b

! Global variables used but not changed
!USE GlobalVmaco,ONLY: delta,glf,glfa,MAXN
! Global variables modified by Subroutine UpdatB
!USE GlobalVmaco,ONLY: b=>a

IMPLICIT NONE
  INTEGER,INTENT(IN):: n   !    number of control variables
  REAL:: bddbd !  bdelta divided by dbd
  REAL,DIMENSION(MAXN):: bdelta ! product of the variable metric matrix b and delta       
  REAL,PARAMETER:: BFACTR = 0.2 ! b matrix factor. chosen empirically to use in the       
                                !   b-f-g-s formula. see eqn (3.6)                          
  REAL:: dbd !    scalar product of delta*b*delta. see eqn (3.6)          
  REAL:: dbdfac ! product of dbd and bfactr. see eqn (3.6)
  REAL:: deleta ! scalar product of delta and eta. see eqn (3.8)
  REAL:: delgam ! scalar product of delta and gamma where gamma           
                !  corresponds to the difference in gradients. eqn (3.3)   
  REAL,DIMENSION(MAXN):: difgrd ! difference between the current and previous gradient of 
                                ! the Lagrange function. corresponds to gamma in eqn(3.3) 
  REAL,DIMENSION(MAXN):: eta    ! replaces difgrd (i.e. gamma) in the b-f-g-s formula.    
                                !  see eqn (3.5)                                           
  REAL:: etade !  eta divided by deleta
  INTEGER:: i,j
  REAL:: theta !  the theta of eqn (3.6) and (3.7)
  REAL:: thcomp ! compliment of theta (i.e. 1-theta)                      


!D7      EQUATIONS CITED FROM "A FAST ALGORITHM FOR NONLINEARLY
!D7      CONSTRAINED OPTIMIZATION CALCULATIONS" BY M.J.D. POWELL
!-------------------------------------------------------------------------------
!
!...calculate the necessary matrices and information to revise the B matrix
  delgam = ZERO
  dbd = ZERO

      DO 20 I=1,N
        difgrd(i) = glf(i)-glfa(i)
        bdelta(i) = ZERO
        DO 10 J=1,N
          bdelta(i) = bdelta(i) + b(i,j)*delta(j)
   10   END DO

        delgam = delgam + delta(i)*difgrd(i)
        dbd = dbd + delta(i)*bdelta(i)
   20 END DO
!
!...  calculate the vector eta for the b-f-g-s formula
!
      dbdfac = bfactr*dbd
!
      IF (delgam .LT. dbdfac) THEN
        theta = (dbd - dbdfac)/(dbd-delgam)
        thcomp = FP1-theta
        DO 30 i=1,n
          eta(i) = theta*difgrd(i) + thcomp*bdelta(i)
   30   END DO

        deleta = dbdfac
      ELSE
        DO 40 i=1,n
          eta(i) = difgrd(i)
   40   END DO
        deleta = delgam
      END IF
!
!...revise the b matrix as per eqn (3.8)
      DO 60 i=1,n
        bddbd = bdelta(i)/dbd
        etade = eta(i)/deleta
!
        DO 50 j=i,n
          b(i,j) = b(i,j) - bddbd*bdelta(j) + etade*eta(j)
          b(j,i) = b(i,j)
   50   END DO
   60 END DO

  RETURN
END Subroutine UpdatB   ! ------------------------------------------------------

!+
SUBROUTINE VERTEX(MEQ)
! ------------------------------------------------------------------------------
! PURPOSE - Obtain an initial feasible vertex (a point formed by the 
! intersection of nvar inequalities) and works by exchanging constraints 
! in a trial vertex.

! Global variables used but not changed
!USE GlobalVmaco,ONLY: bdl,bdu,cinvrs,d,ibasis,MAXC,MAXX,mtotal,nc,nvar,nvar2
! Global variables modified by Subroutine Vertex
!USE GlobalVmaco,ONLY: delta,iactiv

IMPLICIT NONE
  INTEGER,INTENT(IN):: meq   !  total number of equality constraints

  INTEGER:: i
  INTEGER:: icrow
  INTEGER:: nviol !  number of violated constraints
  REAL,DIMENSION(MAXX):: rhs !    right hand side of the constaint eqns in the basis
  REAL,DIMENSION(MAXC):: snorm !  sum of the violated constraint normals
  REAL:: uicj !   product of the ith row of cinvrs and the normal vectors of constraints in the basis
  REAL:: uicmax ! maximum value of uicj

!D7      EQUATIONS AND ALGORITHM CITED FROM "THE CALCULATION OF FEASIBLE
!D7      POINTS FOR LINEARLY CONSTRAINED OPTIMIZATION PROBLEMS" AERE-
!D7      R.6354 BY R. FLETCHER
!
! called by Quad
! calls Ckvrtx, Desig, DotPrd, NewBas, NoDesg, Vmamsg
  INTRINSIC:: DOT_PRODUCT
!-------------------------------------------------------------------------------
  DO i=1,mtotal ! initialize the status of the constraints as inactive
        iactiv(i) = 1
  END DO

!...get the inverse matrix corresponding to the constraints in the basis
  IF (nc .EQ. 0) THEN
    CALL Nodesg
  ELSE
    CALL Desig(meq)
  END IF

!...get the right hand side of the constraints in the basis
 DO i=1,nvar
    IF (ibasis(i) .GT. nvar2) THEN
      rhs(i) = d(ibasis(i)-nvar2)   !...      equality or inequality constraint in the basis
    ELSE IF (ibasis(i) .GT. nvar) THEN
      rhs(i) = -bdu(ibasis(i)-nvar)   !...      upper boundary constraint in the basis
    ELSE
      RHS(I) = BDL(IBASIS(I))   !...      lower boundary constraint in the basis
    END IF
  END DO

!...calculate the position of the vertex
  DO I=1,NVAR
!!!        DELTA(I) = DOTPRD(NVAR,CINVRS(1:,I),NVAR,1,RHS(1:),NVAR,1)   ! modified by RLC
    delta(i)=DOT_PRODUCT(cinvrs(1:nvar,i), rhs(1:nvar) )
  END DO

!...A trial vertex has been found. It must be determined if it is a feasible point.
  40 CONTINUE
  CALL Ckvrtx(nviol,snorm)

!... Possible directions of search are rows of cinvrs obtained by
!...  removing a constraint. Calculate the optimum direction.
  IF (nviol .NE. 0) THEN
    uicmax = ZERO
    DO i=1,nvar
      IF (iactiv(ibasis(i)) .NE. -1) THEN
        uicj=DOT_PRODUCT(cinvrs(i,1:nvar), snorm(1:nvar) )
        IF (uicj .GT. uicmax) THEN
          uicmax = uicj
          icrow = i
        END IF
      END IF
    END DO

    IF (uicmax .LE. ZERO) THEN
      CALL Vmamsg(11,0,ZERO)
    ELSE
      CALL Newbas (icrow)
      GO TO 40
    END IF
  END IF

  RETURN
END Subroutine Vertex   ! ------------------------------------------------------

!+
SUBROUTINE Vmaco(n,m,meq,x,f,g,c,cnorml,maxfun,acc,istat,iprint)
! ------------------------------------------------------------------------------
! PURPOSE - Variable metric algorithm for constrained optimization.
!  Calculate the least value of a function of several variables subject to 
!  equality and inequality constraints.

! Global variables used but not changed
!USE GlobalVmaco,ONLY: cn,gdelta,MAXC,MAXX,nvma,vlam
! Global variables modified by Subroutine Vmaco
!USE GlobalVmaco,ONLY: iter,ibugfg

IMPLICIT NONE
  REAL,INTENT(IN):: ACC   !  controls the final accuracy for convergence. 
    ! Convergence occurs when the value for the change in objective function 
    ! plus suitably weighted multipliers of the constraint function differs 
    ! from its previous value by at most acc.
  REAL,INTENT(IN),DIMENSION(:):: c   !    vector of constraint variables.
  REAL,INTENT(IN),DIMENSION(:,:):: cnorml   ! matrix of constraint normals
  REAL,INTENT(IN):: f    !   value of the objective function (optimized variable)
  REAL,INTENT(IN),DIMENSION(:):: g    !   gradient of the objective function
  INTEGER,INTENT(IN):: iprint   ! flag to indicate amount of printout
                                    ! = 0 no printout
                                    ! = 1 values of x, c, nd f are printed
                                    ! = 2 debugging printout
  INTEGER,INTENT(IN OUT):: istat   ! status of the optimization algorithm.
    ! =-2 initial pass through vmaco with the b matrix already initialized
    ! =-1 initial pass through vmaco
    ! = 0 continue algorithm calculations
    ! = 1 converge within required acurracy
    ! > 1 error encountered
  INTEGER,INTENT(IN):: m    !   total number of constraints (i.e. constraint equations)
  INTEGER,INTENT(IN):: maxfun  ! maximum number of times vmaco can be called
  INTEGER,INTENT(IN):: meq     ! total number of equality constraints
  INTEGER,INTENT(IN):: n       ! number of control variables
!  INTEGER,INTENT(IN):: nsize   ! corresponds to the maximum dimension of the innermost
    ! array size from the cnorml matrix (i.e. max control variable size). This allows 
    ! correct communication between the cnorml and cn matrix.
!!! NOTE - PDAS version makes cnorml an assumed shape array and 
!!!   therefore NSIZE is no longer needed.

  REAL,INTENT(IN OUT),DIMENSION(:):: x   !    vector of control variables


  LOGICAL::BFLAG   ! flag for whether the b matrix should be updated
  REAL:: fdelta ! absolute change in the objective function f. Checked against
                ! accuracy desired for convergence criteria.
  INTEGER,PARAMETER:: MAXCM1 = MAXC-1  ! value of maxc subtract 1
  INTEGER,PARAMETER:: MAXXM1 = MAXX-1 ! value of maxx subtract 1
  INTEGER:: i,j
  INTEGER:: ibug
  INTEGER:: nvmitr ! value of nvma at the start of an iteration
  INTEGER,PARAMETER:: NVMLIN =5 ! controls the error return from the line search
                   !  objective function. see eqn (4.2)
  REAL:: xdummy
!D7      ALGORITHM AND EQUATIONS CITED BASED ON THE PAPER "A FAST
!D7      ALGORITHM FOR NONLINEARLY CONSTRAINED OPTIMIZATION
!D7      CALCULATIONS" BY M.J.D POWELL PRESENTED AT THE 1977 DUNDEE
!D7      CONFERENCE ON NUMERICAL ANALYSIS

!-------------------------------------------------------------------------------
!...check n, m, and meq are logical values (all istat > 1 are error exits)
!
      IF (n .LE. 0) THEN
        istat=2
        CALL VMAMSG (istat,0,ZERO)
        GO TO 999
      ELSEIF (m .LT. meq) THEN
        istat=3
        CALL VMAMSG (istat,0,ZERO)
        GO TO 999
      ELSEIF (meq .LT. 0) THEN
        istat = 4
        CALL VMAMSG (istat,0,ZERO)
        GO TO 999
      ELSEIF (n .GT. MAXXM1) THEN
        istat = 12
        CALL VMAMSG (istat, MAXXM1, ZERO)
        GO TO 999
      ELSEIF (m .GT. MAXCM1) THEN
        istat = 13
        CALL VMAMSG (istat, MAXCM1, ZERO)
        GO TO 999
      ENDIF
!
!...equate cnorml to cn to allow for adding an extra
!...variable for infeasibility
!
      IF (M .GE. 1) THEN
        DO 5 I=1,N
          DO 4 J=1,M
            CN(I,J)= CNORML(I,J)
    4     CONTINUE
    5   CONTINUE
      ENDIF

! write(*,*) "Before initial call to Initvm, istat,n,m=", istat, n, m
  IF (istat.LE.-1) THEN   !  currently on the initial iteration pass
    CALL Initvm (istat,n,m)
!	write(*,*) 'After call to Initvm, istat,nvar=', istat,nvar
    GO TO 20
  END IF

  IF (iprint .EQ. 2) THEN
    IF (ibugfg .EQ. 1) THEN
      CALL Vmaout (c,f,g,m,n,x,iprint)
      ibugfg = 0
    END IF
  END IF

  IF (nvma .GE. MAXFUN) THEN
    istat = 5
    CALL Vmamsg (istat,MAXFUN,ZERO)
    GO TO 999
  ELSEIF (nvma .GT. NVMITR+NVMLIN) THEN
    istat = 6
    CALL Vmamsg (istat,NVMLIN,ZERO)
    GO TO 999
  END IF

10 CONTINUE   !...calculate the line search objective function
  bflag = .FALSE.

  CALL Lnsrch(c,f,g,istat,m,maxfun,meq,n,nvmitr,nvmlin,x,bflag,iprint)
  IF (bflag .AND. (istat .LE. 0) ) THEN
    CALL UPDATB (N)   !...    update the b matrix
  ELSE
    ibugfg = 1
    GO TO 999
  END IF

20 CONTINUE
  iter = iter + 1
! write(*,*) 'At 20 in Vmaco, iter,nvar=', iter,nvar  

  IF (iprint .GE. 1) THEN
    CALL Vmaout (c,f,g,m,n,x,iprint)
    IF (iprint .EQ. 2) THEN
      ibug = 2
!!!      write(*,*) 'In Vmaco, calling Bugaid with n=',n
      CALL Bugaid (ibug,n,xdummy)
    END IF
  END IF

  CALL Vmaset(c,g,istat,m,meq,n,iprint)

  IF (istat .GE. 2) THEN
    GO TO 999
  END IF

  nvmitr = nvma
  CALL Gradlf (g,m,n,iprint)
  CALL Stovmu (g,m,n,x)

!...calculate whether the iteration has converged on the optimal solution
  fdelta = ABS(gdelta)
  DO i=1,m
   fdelta = fdelta + ABS(vlam(i)*c(i))
  END DO

  IF (fdelta .LE. acc) THEN
    istat = 1   !  convergence has ocurred
    IF (iprint .GE. 1) CALL Vmamsg (istat,0,ZERO)
  ELSE
    GO TO 10   !  do a line search for new values
  END IF

999 CONTINUE

  RETURN
END Subroutine Vmaco   ! -------------------------------------------------------

!+
SUBROUTINE VmaMsg(icode,num,xnum)
! ---------------------------------------------------------------------------
! PURPOSE - Write diagnostic messages
! NOTES - Originally written by J.D. Frick of MDAC/Houston, Jan. 1986
!   Rewritten Dec 2005 by Ralph Carmichael, as Absoft Fortran 95 would not
!   accept the DATA statements in the original code.
IMPLICIT NONE
  INTEGER,INTENT(IN):: icode   ! message code number
  INTEGER,INTENT(IN):: num     ! integer number
  REAL,INTENT(IN):: xnum   ! real number

  CHARACTER(LEN=*),PARAMETER:: BORDER = &
    '------------------------------------------------------------------------'

!----------------------------------------------------------------------------
  WRITE(*,*) BORDER
  SELECT CASE(icode)
    CASE(1)
      WRITE(*,*) 'The variable metric algorithm has determined an optimal'
      WRITE(*,*) 'solution within the accuracy desired.'

    CASE(2)
      WRITE(*,*) '*** ERROR *** the input number of control variables is less than one.'

    CASE(3)
      WRITE(*,*) '*** ERROR *** the input number of constraints is less than'
      WRITE(*,*) ' the number of equality constraints.'

    CASE(4)
      WRITE(*,*) '*** ERROR *** the input number of equality constraints is less then zero'

    CASE(5)
      WRITE(*,*) '*** WARNING *** maximum number of calls (=', num, ')'
      WRITE(*,*) '  to Vmaco has been exceeded.'

    CASE(6)
      WRITE(*,*) '*** ERROR *** line search requires ', num, ' calls to Vmaco.'

    CASE(7)
      WRITE(*,*) '*** ERROR *** an uphill search direction has been encountered.'
      WRITE(*,*) '  Usually this is due to loss of accuracy in the quadratic'
      WRITE(*,*) '  programming calculation.'

    CASE(8)
      WRITE(*,*) '*** ERROR *** the given constraint seems to be inconsistent.'

    CASE(9)
      WRITE(*,*) '*** ERROR *** an artificial bound is active.'
      WRITE(*,*) '  The predicted change in the variable exceeds ', xnum

    CASE(10)
      WRITE(*,*) '*** WARNING *** the solution may degenerate into a local minimum.'

    CASE(11)
      WRITE(*,*) '*** WARNING *** no feasible vertex found.'

    CASE(12)
      WRITE(*,*) '*** ERROR *** the maximum size of the control'
      WRITE(*,*) '  array (nsize) is greater than ', num

    CASE(13)
      WRITE(*,*) '*** ERROR *** the maximum size of the constraint'
      write(*,*) '  array ( m ) is greater than ', num

    CASE(14)
      WRITE(*,*) '*** ERROR *** quadratic search requires more than ', num
      WRITE(*,*) '  loops. Limit or adjust the appropriate controls.'

    CASE(15)
      WRITE(*,*) '*** WARNING *** rank of an inverted matrix in'
      SELECT CASE(num)
        CASE(1)
          WRITE(*,*)  '  subroutine Desig has decreased.'
        CASE(2)
          WRITE(*,*)  '  subroutine Reset has decreased.'
      END SELECT

  END SELECT
  WRITE(*,*) BORDER
  
  RETURN
END Subroutine VmaMsg   !----------------------------------------------------

!+
SUBROUTINE Vmaout(c,f,g,m,n,x,iprint)
! ------------------------------------------------------------------------------
! PURPOSE - Provide the printout of the objective function, control variable, 
!  and constraint variable values upon each return from Vmaco.

! Global variables used but not changed
!USE GlobalVmaco,ONLY: iter,nvma
! No global variables modified by Subroutine VmaOut

IMPLICIT NONE
  REAL,INTENT(IN),DIMENSION(:):: c   !    vector of constraint variables
  REAL,INTENT(IN):: f      !    value of the objective function (optimized variable)
  REAL,INTENT(IN),DIMENSION(:):: g   !    gradient of the objective function
  INTEGER,INTENT(IN):: m   !    total number of constraints (i.e. constraint equations)
  INTEGER,INTENT(IN):: n   !    number of control variables
  REAL,INTENT(IN),DIMENSION(:):: x   !    vector of control variables
  INTEGER,INTENT(IN):: iprint  ! flag to indicate amount of printout
                               !  =0 no printout
                               !  =1 values of x, c, and f are printed
                               !  =2 debugging printout

  CHARACTER(LEN=*),PARAMETER:: FMT60 = '(/" OBJECTIVE FUNCTION VALUE IS: ", ES16.7)'

  INTEGER:: i
  INTEGER:: ndummy
  REAL,DIMENSION(1):: xdummy
! called by Vmaco
!-------------------------------------------------------------------------------


      WRITE(*,  5)
      WRITE(*, 10) iter, nvma
      WRITE(*, 20)
      CALL Colum3 (6,n,x,ndummy)

      IF (m .GE. 1) THEN
        WRITE(*, 40)
        CALL Colum3 (7,m,c,ndummy)
      END IF

      WRITE(*, FMT60) f
      IF (iprint .EQ. 2) THEN
        WRITE(*,70)
        CALL Colum3 (8,n,g,ndummy)
        WRITE(*,80)
        DO 3 I=1,m
          CALL Colum3 (9,n,xdummy,i)
    3   END DO
      END IF
      WRITE(*,5)


  RETURN

    5 FORMAT (1X,72('*'))
   10 FORMAT (/5X, 'ITERATIONS =', I5, 5X,                             &
     &        'NUMBER OF CALLS TO VMACO =', I5)
   20 FORMAT (/1X,' THE VECTOR OF CONTROL VARIABLE VALUES ARE:'/)
   40 FORMAT (/1X,' THE VECTOR OF CONSTRAINT VALUES ARE:'/)
!   60 FORMAT (/1X,' THE OBJECTIVE FUNCTION VALUE IS:'//5X, 3HF= ,      &
!     &        ES16.7)
   70 FORMAT (/2X,'  THE VECTOR OF GRADIENT VALUES ARE:'/)
   80 FORMAT (/2X,'  THE JACOBIAN MATRIX IS:'/)

END Subroutine VmaOut   ! ------------------------------------------------------

!+
SUBROUTINE VmaSet(c,g,istat,m,meq,n,iprint)
! ------------------------------------------------------------------------------
! PURPOSE - Set up and calculate variables needed by Vmaco and the
!  quadratic programming routines

! Global variables used but not changed
!USE GlobalVmaco,ONLY: cstar,delta,ghat,ibasis,nc,nvar,nvar2,vlam
! Global variables modified by Subroutine VmaSet
!USE GlobalVmaco,ONLY: bdu,cn,d,gm,vlam

IMPLICIT NONE
  REAL,INTENT(IN),DIMENSION(:):: c   !    vector of constraint variables
  REAL,INTENT(IN),DIMENSION(:):: g   !    gradient of the objective function
  INTEGER,INTENT(IN OUT):: istat   ! status of the optimization algorithm
    ! =-2 initial pass through vmaco with the b matrix already initialized
    ! =-1 initial pass throught vmaco
    ! = 0 continue algorithm calculations
    ! = 1 converge within required accuracy
    ! > 1 error encountered
  INTEGER,INTENT(IN):: m   !    total number of constraints (i.e. constraint equations)
  INTEGER,INTENT(IN):: meq !    total number of equality constraints
  INTEGER,INTENT(IN):: n   !    number of control variables
  INTEGER,INTENT(IN):: iprint   ! flag to indicate amount of printout
                                ! = 0 no printout
                                ! = 1 values of x, c, and f are printed
                                ! = 2 debugging printout

  REAL,PARAMETER:: FEASP = 0.9 ! scaling factor used to achieve feasibility. 
                               ! Corresponds to greek letter xi, 
                               ! eqn (2.7) in Powell's paper.
  INTEGER:: ifeas   ! flag used to mark feasibility conditions
  REAL,PARAMETER:: VLARGE = 1.0E+6   ! a very large number
  REAL,PARAMETER:: VSMALL = 1.0E-6   ! a very small positive number
  INTEGER:: i,j,k

!D7      EQUATIONS AND ALGORITHM CITED FROM "A GENERAL QUADRATIC
!D7      ALGORITHM" AERE-TP.401 BY R. FLETCHER
! called by Vmaco
! calls Quad, StInit, VmaMsg
!-------------------------------------------------------------------------------
!write(*,*) "Entering VmaSet with nvar=", nvar
  IF (istat .LE. -1) THEN   !...initialize necessary values
    CALL Stinit(istat,m,meq,n,vlarge)
  END IF

  DO i=1,N
    gm(i) = -g(i)
  END DO

  IF (m .GE. 1) THEN
    DO i=1,m
      IF (i .LE. meq .OR. c(i) .LT. ZERO) THEN
        d(i) = ZERO
        cn(nvar,i) = c(i)
      ELSE
        d(i) = -c(i)
        cn(nvar,i) = ZERO
      END IF
      vlam(i) = ZERO
    END DO
  END IF

  bdu(nvar) = FP1
  ifeas = -1
30 CONTINUE


  CALL QUAD(istat,meq,iprint)      ! the quadratic programming routine is called

!...determine if the required feasibility conditions hold
  IF (delta(nvar) .LE. VSMALL) THEN
    istat = 8
    CALL Vmamsg(istat,0,ZERO)
    GO TO 999
  END IF

  DO i=1,nc
    IF (ibasis(i) .LT. nvar2) THEN
      ISTAT = 9
      CALL Vmamsg(istat,0,VLARGE)
      GO TO 999
    ELSEIF (ibasis(i) .EQ. nvar2) THEN
      ifeas = 1
    END IF
  END DO

  IF (ifeas .EQ. 0) THEN
    istat = 8
    CALL Vmamsg (istat,0,ZERO)
    GO TO 999
  ELSEIF (ifeas .LT. 0) THEN
    bdu(nvar) = feasp*delta(nvar)
    ifeas = 0
    GO TO 30
  END IF


  DO i=1,nc    ! calculate the Lagrange multipliers...
    k = ibasis(i) - nvar2
    IF (k .GT. 0) THEN
      DO J=1,N
        vlam(k) = vlam(k) + cstar(i,j)*ghat(j)
      END DO
    END IF
  END DO

  999 CONTINUE

  RETURN
END Subroutine VmaSet   ! ------------------------------------------------------

END Module VmacoProcedures   ! =================================================
