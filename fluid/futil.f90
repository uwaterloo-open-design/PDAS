
!+
MODULE CommonFludpc
! ------------------------------------------------------------------------------
! PURPOSE - Emulate the old common block /FLUDPC/

!      DIMENSION G(29,4),TT(27),DD(30),GG(27,30,8),SAT(40,8),C(8),M(4)
!       INTEGER ENTRY,ERROR
!       DIMENSION PROPS(8)
!       LOGICAL VAPOR
!!!       COMMON /FLUDPC/ G,TT,DD,GG,SAT,C,M,NG,N1,N2,N3,NS
  REAL,DIMENSION(29,4):: g
  REAL,DIMENSION(27):: tt
  REAL,DIMENSION(30):: dd
  REAL,DIMENSION(27,30,8):: gg
  REAL,DIMENSION(40,8):: sat
  REAL,DIMENSION(8):: c
  INTEGER,DIMENSION(4):: mIndexValue
  INTEGER:: ng
  INTEGER:: n1
  INTEGER:: n2
  INTEGER:: n3
  INTEGER:: ns

!-------------------------------------------------------------------------------

END Module CommonFludpc   !=====================================================

!+
MODULE FluidProperties
! ------------------------------------------------------------------------------
! PURPOSE - Hold all of the subroutines that define the tables
!!! USE FluidProcedures
!-------------------------------------------------------------------------------

CONTAINS
INCLUDE 'newar.f90'
INCLUDE 'newch4.f90'
INCLUDE 'newco2.f90'
INCLUDE 'newdryair.f90'
INCLUDE 'newf2.f90'
INCLUDE 'newn2.f90'
INCLUDE 'newO2.f90'
INCLUDE 'newph2.f90'
INCLUDE 'newsteam.f90'
!INCLUDE 'others.f90'

END Module FluidProperties
!+
MODULE FluidUtilities
! ---------------------------------------------------------------------------

!      WRITTEN BY: THEODORE E. FESSLER, NASA TM X-3572, AUGUST 1977
!      REWRITTEN BY: LT. MARK D. KLEM & MARGARET P. PROCTOR, JANUARY 198

  IMPLICIT NONE

  INTEGER,PARAMETER,PRIVATE:: SP=KIND(1.0), DP=KIND(1.0D0)

CONTAINS

!+
PURE FUNCTION Limit(imin, i, imax) RESULT(k)
! ------------------------------------------------------------------------------
! PURPOSE - If integer i is in [imin,imax], return that value; otherwise, move
!  it to nearest boundary.
  INTEGER,INTENT(IN):: imin,i,imax
  INTEGER:: k
  INTRINSIC:: MAX,MIN
!-------------------------------------------------------------------------------
  k=MAX(imin,MIN(i,imax))
  RETURN
END Function Limit   ! =========================================================

!+
PURE FUNCTION RLimit(amin, a, amax) RESULT(b)
! ------------------------------------------------------------------------------
! PURPOSE - If real number a is in [amin,amax], return that value; otherwise, move
!  it to nearest boundary.

  REAL,INTENT(IN):: amin,a,amax
  REAL:: b
!-------------------------------------------------------------------------------
  b=MAX(amin,MIN(a,amax))
  RETURN
END Function RLimit   ! ========================================================

!+
FUNCTION Lookup(h,htab) RESULT(i)
! ------------------------------------------------------------------------------
! PURPOSE - Find the unique value of i such that htab(i)<=h<htab(i+1).
!  Handle the special cases: i=1 if h<h(1) and i=SIZE(htab) if h>=htab(SIZE(htab))
!  Assumes that htab is a strictly increasing array. No check is made of this.
  REAL,INTENT(IN):: h
  REAL,INTENT(IN),DIMENSION(:):: htab
  INTEGER:: i,j,k
!----------------------------------------------------------------------------
  i=1
  j=SIZE(htab)                               ! setting up for binary search
  DO
    k=(i+j)/2                                              ! integer division
    IF (h < htab(k)) THEN
      j=k
    ELSE
      i=k
    END IF
    IF (j <= i+1) EXIT
  END DO
  RETURN
END Function LookUp   ! =====================================================



!+
PURE FUNCTION TF(P,D,Z) RESULT(f)
! ---------------------------------------------------------------------------
! PURPOSE - Solve the van der Waals equation with correction factor for
!  temperature, given pressure and density. See p. 4 of TM X-3572
  REAL,INTENT(IN):: p,d,z  ! reduced pressure, reduced density, correction
  REAL:: f
!----------------------------------------------------------------------------
  f=(P+3.0*D**2)*(3.0-D)/(8.0*Z*D)
  RETURN
END Function TF   ! =========================================================

!+
PURE FUNCTION PF(T,D,Z) RESULT(f)
! ---------------------------------------------------------------------------
! PURPOSE - Solve the van der Waals equation with correction factor for
!  pressure, given temperature and density. See p. 4 of TM X-3572

  REAL,INTENT(IN):: t,d,z  ! reduced temperature, density, correction
  REAL:: f
!----------------------------------------------------------------------------
  f=8.0*Z*D*T/(3.0-D)-3.0*D**2
  RETURN
END Function PF   ! =========================================================

!+
PURE FUNCTION HF(P,D) RESULT(f)
! ------------------------------------------------------------------------------
! PURPOSE - Solve for the correction to the ideal gas enthalpy using van der 
!  Waals equation with correction factor, given pressure and density. 
!  See p. 7 of TM X-3572. You also need the internal energy and the correction
!  factor to compute the enthalpy.

  REAL,INTENT(IN):: p,d  ! reduced pressure and density
  REAL:: f
!----------------------------------------------------------------------------
  f=1.125*(P/(3.0*D)-D)
  RETURN
END Function HF   ! =========================================================

!+
PURE FUNCTION SF(D) RESULT(f)
! ------------------------------------------------------------------------------
! PURPOSE - Solve for the correction to the ideal gas entropy using van der 
!  Waals equation with correction factor, given pressure and density. 
!  See p. 7 of TM X-3572. You also need the ideal gas entropy and the correction
!  factor to compute the entropy. 

  REAL,INTENT(IN):: d   ! reduced density
  REAL:: f
!----------------------------------------------------------------------------
  f=LOG((3.0-D)/D)
  RETURN
END Function SF   ! =========================================================

!+
FUNCTION CPF(t,d) RESULT(f)
! ------------------------------------------------------------------------------
! PURPOSE - Solve for the final term in the equation for specific heat at 
!  constant pressure.
!  See p. 7 of TM X-3572. 

  REAL,INTENT(IN):: t,d
  REAL:: f
!-------------------------------------------------------------------------------
  f=1.0/(1.0-d*(3.0-d)**2/(4.0*t))
  RETURN
END Function CPF   ! ===========================================================

!+
PURE FUNCTION A2F(t,d,grt) RESULT(f)
! ------------------------------------------------------------------------------
! PURPOSE - Solve for the term under the root in the equation for the speed
!  of sound.
!  See p. 7 of TM X-3572. 

  REAL,INTENT(IN):: t,d   ! reduced temperature and density
  REAL,INTENT(IN):: grt  ! grt is gamma*R*tc
  REAL:: f
!-------------------------------------------------------------------------------
  f=0.375*grt*(8.0*t/(3.0-d)*(1.0+d/(3.0-d))-6.0*d)
  RETURN
END Function A2F   ! ===========================================================



!+
SUBROUTINE NTRP (A,X,N0,IA,OUT,IGOTO,F,FF,VALUE)
! ---------------------------------------------------------------------------
! PURPOSE - DOUBLE 3-POINT INTERPOLATION ROUTINES FOR SUBPROGRAM FLUID.

  REAL,INTENT(IN),DIMENSION(:):: a   ! TABLE OF INDEPENDENT-VARIABLE VALUES
  REAL,INTENT(IN):: x   ! ENTRY VALUE OF INDEPENDENT VARIABLE
  INTEGER,INTENT(IN):: n0   !  (OR N1/N2) = LENGTH OF A TABLE
  INTEGER,INTENT(IN OUT):: ia   !  CURRENT INDEX VALUE ASSOCIATED WITH A TABLE
  INTEGER,INTENT(OUT):: out  !  = 0,1 IF X IS INSIDE/NOT INSIDE LIMITS OF A
  INTEGER,INTENT(IN):: igoto ! ???
  REAL,INTENT(IN),DIMENSION(:):: f
  REAL,INTENT(IN),DIMENSION(:,:):: ff
  REAL,INTENT(OUT):: value

  REAL,SAVE,DIMENSION(4,3):: c
  REAL,SAVE,DIMENSION(4):: c0,c1,c2
  INTEGER,SAVE,DIMENSION(3):: i1,i2,nn
  LOGICAL,SAVE:: NEW12

  INTEGER:: i,j,ij,k,l,m   ! various counters
  INTEGER:: l1,l2, m1,m2, mm1,mm2
  INTEGER:: na
  INTEGER:: ni,nj

  REAL,SAVE,DIMENSION(16):: cc
!!!      EQUIVALENCE (C0,C(1,1)), (C1,C(1,2)), (C2,C(1,3)), (L1,I1(1)),    &
!!!     &            (L2,I2(1)),  (M1,I1(2)),  (M2,I2(2)),  (MM1,I1(3)),   &
!!!     &            (MM2,I2(3)), (NI,NN(2)),  (NJ,NN(3))


      EQUIVALENCE (c0(1),c(1,1)), (c1(1),c(1,2)), (c2,c(1,3))
      ! so, c0 is col. 1 of c, c1 is col. 2, c2 is col.3

      EQUIVALENCE (L1,i1(1)), (m1,i1(2)), (mm1,i1(3))
      EQUIVALENCE (L2,i2(1)), (m2,i2(2)), (mm2,i2(3))
      EQUIVALENCE             (ni,nn(2)), (nj,nn(3))

  REAL:: p1,p2,p3
  REAL:: a11,a12,a13,a14, a21,a22,a23,a24, a31,a32,a33,a34, a41,a42,a43,a44

!----------------------------------------------------------------------------
  IF (IGOTO==4) THEN
    VALUE=0.0    !.....EVALUATION OF FUNCTION IN ONE VARIABLE.
    K=0
    DO I=L1,L2
      K=K+1
      VALUE=VALUE+C0(K)*F(I)   !.....ARGUMENT F IS THE TABLE OF DEPENDENT-VARIABLE VALUES.
    END DO
    RETURN
  END IF

  IF (IGOTO==5) THEN     !.....ENTRY FOR EVALUATION OF FUNCTIONS IN TWO VARIABLES.
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
        VALUE=VALUE+CC(IJ)*FF(I,J) !.....ARGUMENT FF IS THE (2-D) TABLE OF DEPENDENT-VARIABLE VALUES.
      END DO
    END DO
    RETURN
  END IF

  SELECT CASE(IGOTO)
    CASE (1)      !.....ENTER HERE TO CALCULATE COEFFICIENTS FOR GNTRP.
      K=1
      NA=N0
    CASE (2)   !.....ENTER HERE TO CALCULATE FIRST COEFFICIENTS FOR GGNTRP.
      K=2
      NA=N0
      NEW12=.TRUE.
    CASE (3)  !.....ENTER HERE TO CALCULATE SECOND COEFFICIENTS FOR GGNTRP.
      K=3
      NA=N0
      NEW12=.TRUE.
  END SELECT
                  !     PROCEDURE (CALCULATE COEFFICIENTS FOR ONE DIMENSION)
  L=LIMIT(1,NA,3) !     PROCEDURE (CALCULATE INDEX VALUES)
  M=NA-2

  IF (X < A(1)) THEN      !     PROCEDURE (TABLE LOOK-UP)
    OUT=1
    IA=1
  ELSE
    IF (X > A(NA)) THEN
      OUT=1
      IA=NA
    ELSE
      OUT=0
      IA=LIMIT(1,IA,NA)              ! assumes there is a previous value
      IF (X.LT.A(IA)) THEN
   10   IA=IA-1
        IF (X.LT.A(IA)) GO TO 10
      ELSE
   20   IF (X.GT.A(IA+1)) THEN
          IA=IA+1
          GO TO 20
        END IF
      END IF
    END IF
  END IF

  IF (IA.GE.2 .AND. IA.LE.M) L=4
  I=LIMIT(1,IA-1,M)
  I1(K)=I
  I2(K)=I+L-1
  NN(K)=L

  IF (L.GT.1) THEN    !     PROCEDURE (CALCULATE AIJ VALUES)
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

END Subroutine Ntrp   ! =====================================================


!+
SUBROUTINE CUBIC (T,P,ROOTS,NROOTS)
! ---------------------------------------------------------------------------
! PURPOSE - SOLVES FOR DENSITY IN REDUCED VAN DER WAALS GAS:
!        D**3 -3*D**2 +((8*T+P)/3)*D -P = 0
!     IF MORE THAN ONE DISTINCT ROOT EXISTS, ONLY THE MINIMUM AND
!     MAXIMUM ONES ARE RETURNED AND NROOTS=2.

  REAL(KIND=4),INTENT(IN):: t,p
  REAL(KIND=4),DIMENSION(:),INTENT(OUT):: roots   ! dimension=2
  INTEGER,INTENT(OUT):: nroots

  REAL(DP),PARAMETER:: ZERO=0.0,ONE=1.0,TWO=2.0,THREE=3.0,EIGHT=8.0
  REAL(DP),PARAMETER:: PI=3.14159265
  REAL(DP),PARAMETER:: K1=TWO*PI/THREE, K2=K1+K1

  REAL(KIND=8):: a,a3,aa,bb,b,b2,cr
  REAL(KIND=8):: rad,phi
  REAL(KIND=8):: q
  REAL(KIND=8):: r1,r2,r3


  Q=(EIGHT*T+P)/THREE   !     PROCEDURE (REDUCE TO NORMAL FORM)
  A=Q/THREE - ONE
  A3=A*A*A
  B=(Q-P)/TWO-ONE
  B2=B*B
  RAD=A3+B2
  IF (RAD < ZERO) THEN
    PHI=ACOS(SIGN(SQRT(-B2/A3),-B))/THREE !        PROCEDURE (CASE FOR 3 DISTINCT REAL ROOTS)
    CR=TWO*SQRT(-A)
    R1=CR*COS(PHI)
    R2=CR*COS(PHI+K1)
    R3=CR*COS(PHI+K2)
    ROOTS(1)=MIN(R1,R2,R3)+ONE
    ROOTS(2)=MAX(R1,R2,R3)+ONE
    NROOTS=2
  ELSE
    IF (RAD==ZERO) THEN
      IF (B==ZERO) THEN
        IF (RAD.NE.ZERO) RAD=SQRT(RAD) !              PROCEDURE (CASE FOR ONLY ONE REAL ROOT)
        AA=-B+RAD
        BB=-B-RAD
        IF (AA.NE.ZERO) AA=QRT(AA)
        IF (BB.NE.ZERO) BB=QRT(BB)
        ROOTS(1)=AA+BB+ONE
        NROOTS=1
      ELSE
        R1=QRT(B) !              PROCEDURE (CASE FOR 3 REAL ROOTS BUT TWO ARE EQUAL)
        R2=-TWO*R1
        ROOTS(1)=MIN(R1,R2)+ONE
        ROOTS(2)=MAX(R1,R2)+ONE
        NROOTS=2
      END IF
    ELSE
      IF (RAD.NE.ZERO) RAD=SQRT(RAD) !           PROCEDURE (CASE FOR ONLY ONE REAL ROOT)
      AA=-B+RAD
      BB=-B-RAD
      IF (AA.NE.ZERO) AA=QRT(AA)
      IF (BB.NE.ZERO) BB=QRT(BB)
      ROOTS(1)=AA+BB+ONE
      NROOTS=1
    END IF
  END IF

  RETURN

CONTAINS
!+
FUNCTION Qrt(arg) RESULT(f)
! ------------------------------------------------------------------------------
! PURPOSE - Compute cube root
  REAL(KIND=8),INTENT(IN):: arg
  REAL(KIND=8):: f
  REAL(KIND=8),PARAMETER:: THIRD=1.0D0/3.0D0
!-------------------------------------------------------------------------------
  f=SIGN(ABS(ARG)**THIRD, ARG)
  RETURN
END Function Qrt   ! -----------------------------------------------------------

END Subroutine Cubic   ! -------------------------------------------------------



END Module FluidUtilities   ! ===============================================

!+
MODULE FluidProcedures
! ------------------------------------------------------------------------------
! PURPOSE - Encapsulate the FLUID procedure and associated routines

USE CommonFludPC
USE FluidUtilities
IMPLICIT NONE


  INTEGER,PARAMETER:: DBG = 3
!-------------------------------------------------------------------------------

CONTAINS

!+
SUBROUTINE FLUID (TEMP,PRES,DENS,PROPS,NP,ENTRY,VAPOR,ERROR,IGOTO)
! ---------------------------------------------------------------------------
! PURPOSE - A Numerical Interpolation Procedure for Obtaining
!    Thermodynamic and Transport Properties of Fluids
!
! AUTHORS - THEODORE E. FESSLER, NASA Lewis Research Center
!           LT. MARK D. KLEM & MARGARET P. PROCTOR, NASA?
!           Ralph L. Carmichael, Public Domain Aeronautical Software

! REVISION HISTORY
!   DATE  VERS PERSON  STATEMENT OF CHANGES
!   1977   1.0   TEF   Publication of NASA TM X-3572 (August 1977)
! Jan1986  1.1 MDK&MPP Rewritten in Fortran 77
! 11Aug97  1.2   RLC   Acquisition of COSMIC Program LEW-14418
! 20Aug97  1.21  RLC   Removed tabs from source code
! 31Aug97  1.3   RLC   Converted all source to free-format Fortran 90
! 15Dec97  1.31  RLC   Added + in all fields with E 0 as exponent
! 15Dec97  1.32  RLC   DUMOUT (numerous places) is INTEGER; RDUMOUT is REAL
! 18Dec97  1.4   RLC   Added INTERFACE where needed for proper F90 usage
! 02Jan05  1.5   RLC   Changed INTENT of props to IN OUT
! 11Jun10  1.6   RLC   Re-edited the indentation of the text

! The original program description from COSMIC...
! The accurate computation of the thermodynamic and transport properties of
! fluids is a necessity for many engineering calculations. The FLUID program
! was developed to calculate the thermodynamic and transport properties of
! pure fluids in both the liquid and gas phases. Fluid properties are
! calculated using a simple gas model, empirical corrections, and an efficient
! numerical interpolation scheme. FLUID produces results that are in very good
! agreement with measured values, while being much faster than older more
! complex programs developed for the same purpose.
! A Van der Waals equation of state model is used to obtain approximate state
! values. These values are corrected for real gas effects by model correction
! factors obtained from tables based on experimental data. These tables also
! accurately compensate for the special circumstances which arise whenever
! phase conditions occur. Viscosity and thermal conductivity values are
! computed directly from tables. Interpolation within tables is based on
! Lagrange's three point formula. A set of tables must be generated for each
! fluid implemented.
! FLUID currently contains tables for nine fluids including dry air and steam.
! The user can add tables for any fluid for which adequate thermal property
! data is available. The FLUID routine is structured so that it may easily be
! incorporated into engineering programs.


USE FluidUtilities
  IMPLICIT NONE

  REAL,INTENT(IN OUT):: temp    ! fluid temp, Kelvins
  REAL,INTENT(IN OUT):: pres    ! fluid pressure, megaPascals
  REAL,INTENT(IN OUT):: dens    ! density, grams per cubic centimeter
  REAL,INTENT(IN OUT),DIMENSION(:):: props ! the output quantities
                                 ! except that sometimes they are input quantities
                                 !   props(1) = compressibility, pres/(r*dens*temp).
                                 !   props(2) = entropy,  same units as r.
                                 !   props(3) = enthalpy,  same units as r*tc.
                                 !   props(4) = specific heat, cv,  same units as r.
                                 !   props(5) = specific heat, cp,  same units as r.
                                 !   props(6) = sonic velocity,  same units as sqrt(r*tc).
                                 !   props(7) = viscosity, grams/cm-sec
                                 !   props(8) = thermal conductivity, watts/cm-k

  INTEGER,INTENT(IN):: np        ! how many output quantities are wanted
  INTEGER,INTENT(IN):: entry  ! integer that specifies which variables are input:
                              !  = 1 if temperature and density are given.
                              !  = 2 if pressure and density are given.
                              !  = 3 if temperature and pressure are given.
                              !  = 4 if pressure and entropy are given.
                              !  = 5 if pressure and enthalpy are given.
  LOGICAL,INTENT(OUT):: vapor
  INTEGER,INTENT(OUT):: error
  INTEGER,INTENT(IN):: igoto  


!.....CALLING ARGUMENTS:
!       TEMP = FLUID TEMPERATURE,  SAME UNITS AS TC.
!       PRES = FLUID PRESSURE,  SAME UNITS AS PC.
!       DENS = FLUID DENSITY,  SAME UNITS AS DC.
!       PROPS = OTHER FLUID PROPERTIES:
!         PROPS(1) = COMPRESSIBILITY, PRES/(R*DENS*TEMP).
!         PROPS(2) = ENTROPY,  SAME UNITS AS R.
!         PROPS(3) = ENTHALPY,  SAME UNITS AS R*TC.
!         PROPS(4) = SPECIFIC HEAT, CV,  SAME UNITS AS R.
!         PROPS(5) = SPECIFIC HEAT, CP,  SAME UNITS AS R.
!         PROPS(6) = SONIC VELOCITY,  SAME UNITS AS SQRT(R*TC).
!         PROPS(7) = VISCOSITY, GRAMS/CM-SEC
!         PROPS(8) = THERMAL CONDUCTIVITY, WATTS/CM-K
!       NP = NUMBER OF PROPERTIES TO BE CALCULATED.
!       ENTRY = INTEGER THAT SPECIFIES WHICH VARIABLES ARE INPUT:
!             = 1 IF TEMPERATURE AND DENSITY ARE GIVEN.
!             = 2 IF PRESSURE AND DENSITY ARE GIVEN.
!             = 3 IF TEMPERATURE AND PRESSURE ARE GIVEN.
!             = 4 IF PRESSURE AND ENTROPY ARE GIVEN.
!             = 5 IF PRESSURE AND ENTHALPY ARE GIVEN.
!       VAPOR = .TRUE. IF THE FLUID IS SATURATED.  IN THAT CASE,
!               VALUES OF LIQUID AND GAS PHASES ARE GIVEN IN THE
!               COMMON BLOCK /FLUIDC/.
!       ERROR = ERROR FLAGS (BITS -- LEAST SIGNIFICANT = 1):
!             ** IF ERROR = 0,  ALL IS WELL **
!         BIT 1 = OUT OF RANGE IN SAT TABLE.
!         BIT 2 = OUT OF RANGE IN G TABLE.
!         BIT 3 = OUT OF RANGE IN TT TABLE.
!         BIT 4 = OUT OF RANGE IN DD TABLE.
!         BIT 5 = CONVERGENCE NOT ACHIEVED IN SOLVING INVERSE FUNCTION.


!..... Quantities expected to be found in CommonFludPC
!       G = (NON-DIMENSIONAL) TABLES IN TEMPERATURE ONLY.
!         G(1,1) = ARGUMENT VALUES (TEMPERATURE).
!         G(1,2) = ENTROPY, S0.
!         G(1,3) = ENERGY, U0.
!         G(1,4) = SPECIFIC HEAT AT CONSTANT VOLUME, CV0.
!       NG = LENGTH OF EACH G TABLE.
!       GG = TABLES IN TEMPERATURE AND DENSITY.
!         GG(1,1,1) = Z1, THE MODEL-DEPARTURE FACTOR.
!         GG(1,1,2) = Z2 FACTOR FOR ENTROPY (S).
!         GG(1,1,3) = Z3 FACTOR FOR ENTHALPY (H).
!         GG(1,1,4) = Z4 FACTOR FOR SPECIFIC HEAT (CV).
!         GG(1,1,5) = Z5 FACTOR FOR SPECIFIC HEAT (CP).
!         GG(1,1,6) = Z6 FACTOR FOR VELOCITY OF SOUND.
!         GG(1,1,N) = OTHER FACTORS (EG. TRANSPORT PROPERTIES).
!       TT = TEMPERATURE ARGUMENTS (REDUCED) FOR GG TABLES.
!       DD = DENSITY ARGUMENTS (REDUCED) FOR GG TABLES.
!       N1 = LENGTH OF TT TABLE AND FIRST DIMENSION OF GG TABLE.
!       N2 = LENGTH OF DD TABLE AND SECOND DIMENSION OF GG TABLE.
!       N3 = THIRD DIMENSION (NUMBER OF) GG TABLES.
!       SAT = (NON-DIMENSIONAL) SATURATION TABLES.
!         SAT(1,1) = SATURATION TEMPERATURE.
!         SAT(1,2) = SATURATION PRESSURE.
!         SAT(1,3) = SATURATED LIQUID DENSITY.
!         SAT(1,4) = SATURATED GAS DENSITY.
!         SAT(1,5) = SATURATED LIQUID ENTROPY.
!         SAT(1,6) = SATURATED GAS ENTROPY.
!         SAT(1,7) = SATURATED LIQUID ENTHALPY.
!         SAT(1,8) = SATURATED GAS ENTHALPY.
!       NS = LENGTH OF EACH SATURATION TABLE.
!       C = MISCELLANEOUS CONSTANTS:
!         C(1) = SPECIFIC GAS CONSTANT, R.
!         C(2) = CRITICAL TEMPERATURE, TC.
!         C(3) = CRITICAL PRESSURE, PC.
!         C(4) = CRITICAL DENSITY, DC (OR PSEUDO VALUE FOR BEST RANGE).
!         C(5) = CONVERGENCE REQUIREMENT FOR TEMPERATURE.
!         C(6) = CONVERGENCE REQUIREMENT FOR DENSITY.
!         C(7) = CONVERGENCE REQUIREMENT FOR ENTROPY.
!         C(8) = CONVERGENCE REQUIREMENT FOR ENTHALPY.
!       L = MEMORY (4 WORDS) FOR TABLE INDICES.

!.....VARIABLES IN COMMON BLOCK /FLUIDC/:
!       GAMMA = RATIO OF SPECIFIC HEATS, CP/CV.
!       WL = MASS-FRACTION IN THE LIQUID PHASE.
!       WG = MASS-FRACTION IN THE GAS PHASE.
!       DENSL = DENSITY OF SATURATED LIQUID,  SAME UNITS AS DC.
!       DENSG = DENSITY OF SATURATED GAS,     SAME UNITS AS DC.
!       ENTL = ENTROPY OF SATURATED LIQUID,   SAME UNITS AS R.
!       ENTG = ENTROPY OF SATURATED GAS,      SAME UNITS AS R.
!       ENTHL = ENTHALPY OF SATURATED LIQUID, SAME UNITS AS R*TC.
!       ENTHG = ENTHALPY OF SATURATED GAS,    SAME UNITS AS R*TC.

  INTEGER,PARAMETER:: BIT1=1,BIT2=2,BIT3=4,BIT4=8,BIT5=16

  REAL:: a,a2
  REAL:: cv,cp
  REAL:: cv0
  INTEGER:: ERROUT    ! never referenced???
  INTEGER:: DUMN0,DUMIA
  REAL:: DUMF(40),DUMFF(27,30),DUMA(40)
  REAL:: dumval,dumx,dg
  REAL:: d1,dl
  REAL:: deld,delx
  REAL:: dlim
  REAL:: dxds
  REAL:: dxmax = 3.0   ! might be a constant (PARAMETER)
  REAL:: grt
  REAL:: h1
  REAL:: hin,hl,hg,h
  REAL:: delh
  REAL:: p1,p2, t1,t2, delp,dels,delt
  REAL:: small= 0.001
  REAL:: tsat
  REAL:: s,s0,s1
  REAL:: sin,sl,sg
  REAL:: u0
  REAL:: huge=1.0E37
  REAL:: big = 1000.0
  REAL:: r,tc,pc,dc
  REAL:: t,d,p,d2
  REAL:: arg0,arg1,arg2
  INTEGER:: dumout           ! not sure of this ***************
  REAL:: rdumout
  REAL:: psat
  REAL:: dxdh
  REAL:: x,x1
  REAL:: z,z2,z3,z4,z5,z6

  INTEGER:: iprop
  INTEGER:: i2,i3,i4,i5
  INTEGER:: k0,k1,k2,ks
  INTEGER,PARAMETER:: MAX2=12, MAX3=12, MAX4=12, MAX5=12
  INTEGER:: nr
  INTEGER:: nprops
  INTEGER:: nextp
  LOGICAL:: GAS

  REAL:: ROOTS(2)


!!!  REAL:: G(29,4),GG(27,30,8),TT(27),DD(30),SAT(40,8),C(8)
!!!  INTEGER:: L(4)
!!!  INTEGER:: n1,n2,n3,ng,ns
!!!  COMMON /FLUDPC/ G,TT,DD,GG,SAT,C,L,NG,N1,N2,N3,NS

  REAL:: GAMMA,WL,WG,DENSL,DENSG,ENTL,ENTG,ENTHL,ENTHG
  COMMON /FLUIDC/ GAMMA,WL,WG,DENSL,DENSG,ENTL,ENTG,ENTHL,ENTHG


! Arithmetic Statement Functions

!      real tf,pf,hf,sf,cpf,a2f

!      TF(P,D)=(P+3.0*D**2)*(3.0-D)/(8.0*Z*D)
!      PF(T,D)=8.0*Z*D*T/(3.0-D)-3.0*D**2
!      HF(P,D)=1.125*(P/(3.0*D)-D)
!      SF(D)=LOG((3.0-D)/D)
!      CPF(T,D)=1.0/(1.0-D*(3.0-D)**2/(4.0*T))
!      A2F(T,D)=0.375*GAMMA*R*TC*(8.0*T/(3.0-D)*(1.0+D/(3.0-D))-6.0*D)

!----------------------------------------------------------------------------
  write(dbg,*) 'Entering Fluid, entry=', entry
  write(dbg,*) 'Entering Fluid, igoto,np=', igoto,np
  write(dbg,*) 'Entering Fluid, temp,dens,pres=', temp,dens,pres


!.....initialize normalizing values. Do this every time. RLC 12Jun2010
  R=C(1)
  TC=C(2)
  PC=C(3)
  DC=C(4)

  IF (IGOTO.EQ.1) THEN
    IF (NP.GT.N3) THEN       !     PROCEDURE (MAIN ENTRY INITIALIZATION)
      WRITE(*,*) "Fatal Error in Fluid, np > n3", np,n3
      STOP "Fatal Error #1 in Fluid"
    END IF
    IF (ENTRY < 1 .OR. ENTRY > 5) THEN
      WRITE(*,*) "Fatal error in Fluid. Bad Entry"
      STOP "Fatal Error #2 in Fluid"
    END IF
    IF (n2 > SIZE(DD) ) THEN
      WRITE(*,*) "Fatal error in Fluid. n2 > SIZE(DD)"
      STOP "Fatal error #5 in Fluid"
    END IF

    KS=0
    K0=0
    K1=0
    K2=0
    ERROR=0
    ARG0=HUGE
    ARG1=HUGE
    ARG2=HUGE
    T=TEMP/TC
    D=DENS/DC
    P=PRES/PC
    IF (ENTRY.EQ.1) THEN     !     Case 1  given temperature and density.
      write(dbg,*) 'entry=1, reduced t,d,p=', t,d,p
      IF (T.LT.1.0) THEN
        CALL NTRP(SAT(:,1),T,NS,mIndexValue(1),KS,1,DUMF,DUMFF,DUMVAL) ! test for vapor condition
        CALL NTRP(DUMA,DUMX,DUMN0,DUMIA,DUMOUT,4,SAT(:,3), DUMFF,DL) ! saturated liquid density
        CALL NTRP(DUMA,DUMX,DUMN0,DUMIA,DUMOUT,4,SAT(:,4), DUMFF,DG) ! saturated gas density
        VAPOR=D.GT.DG .AND. D.LT.DL
      ELSE
        VAPOR=.FALSE.
      END IF
      write(dbg,*) 'vapor test', vapor

      IF (VAPOR) THEN  !     PROCEDURE (GET SATURATION PRESSURE, PSAT)
                           !     PROCEDURE (GET Z BY DOUBLE INTERPOLATION AT (T,D))
                           !     PROCEDURE (SET UP FOR DOUBLE INTERPOLATION AT (T,D))
        CALL NTRP(DUMA,DUMX,DUMN0,DUMIA,DUMOUT,4, SAT(:,2),DUMFF,PSAT)
        P=PSAT
      ELSE
        IF (ARG1.NE.T) THEN
          CALL NTRP(TT,T,N1,mIndexValue(3),K1,2,DUMF,DUMFF,DUMVAL)
          ARG1=T
        END IF
        IF (ARG2.NE.D) THEN
          CALL NTRP(DD,D,N2,mIndexValue(4),K2,3,DUMF,DUMFF,DUMVAL)
          ARG2=D
        END IF
        CALL NTRP(DUMA,DUMX,DUMN0,DUMIA,DUMOUT,5,DUMF, GG(:,:,1),Z)
        write(dbg,*) 'computed Z', z
        Z=RLIMIT(SMALL,Z,BIG)
        write(dbg,*) 'calling Pf with t,d,z=', t,d,z
        P=PF(T,D,Z)
      END IF
      PRES=P*PC
      write(dbg,*) 'reduced and final temp', p,pres
      write(dbg,*) 'end of entry==1'

    ELSE IF (ENTRY.EQ.2) THEN    ! CASE 2  .GIVEN PRESSURE AND DENSITY.
                                 !     PROCEDURE (TEST FOR VAPOR CONDITION AT (P,D))
                                 !     PROCEDURE (SET UP FOR INTERPOLATION OF SATURATION DATA AT P)
                                 !     PROCEDURE (GET SATURATED LIQUID (DL) AND GAS (DG) DENSITIES)
      write(dbg,*) 'entry=2, reduced t,d,p=', t,d,p
      IF (P.LT.1.0) THEN
        CALL NTRP(SAT(:,2),MAX(P,SAT(1,2)),NS,mIndexValue(1),KS,1,DUMF,DUMFF,DUMVAL)
        CALL NTRP(DUMA,DUMX,DUMN0,DUMIA,DUMOUT,4, SAT(:,3),DUMFF,DL)
        CALL NTRP(DUMA,DUMX,DUMN0,DUMIA,DUMOUT,4, SAT(:,4),DUMFF,DG)
        VAPOR=D.GT.DG .AND. D.LT.DL
      ELSE
        VAPOR=.FALSE.
      END IF



      IF (VAPOR) THEN   !     PROCEDURE (GET SATURATION TEMPERATURE, TSAT)
                        !     PROCEDURE (SOLVE FOR T AT (P,D) BY ITERATION)
        CALL NTRP(DUMA,DUMX,DUMN0,DUMIA,DUMOUT,4, SAT(:,1),DUMFF,TSAT)
        T=TSAT
      ELSE
        Z=1.0
        T=TF(P,D,Z)
        DELT=HUGE
        DO 10 I2=1,MAX2 !     PROCEDURE (GET Z BY DOUBLE INTERPOLATION AT (T,D))
          IF (ARG1.NE.T) THEN   !     PROCEDURE (SET UP FOR DOUBLE INTERPOLATION AT (T,D))
            CALL NTRP(TT,T,N1,mIndexValue(3),K1,2,DUMF,DUMFF,DUMVAL)
            ARG1=T
          END IF
          IF (ARG2.NE.D) THEN
            CALL NTRP(DD,D,N2,mIndexValue(4),K2,3,DUMF,DUMFF,DUMVAL)
            ARG2=D
          END IF
          CALL NTRP(DUMA,DUMX,DUMN0,DUMIA,DUMOUT,5,DUMF, GG(:,:,1),Z)
          Z=RLIMIT(SMALL,Z,BIG)
          IF (ABS(DELT).LE.C(5)*T) GO TO 20
          T1=T
          P1=PF(T1,D,Z)
          T=RLIMIT(TT(1),TF(P,D,Z),TT(N1))
          IF (ARG1.NE.T) THEN    !     PROCEDURE (GET Z BY DOUBLE INTERPOLATION AT (T,D))
            CALL NTRP(TT,T,N1,mIndexValue(3),K1,2,DUMF,DUMFF,DUMVAL)
            ARG1=T
          END IF
          IF (ARG2.NE.D) THEN
            CALL NTRP(DD,D,N2,mIndexValue(4),K2,3,DUMF,DUMFF,DUMVAL)
            ARG2=D
          END IF
          CALL NTRP(DUMA,DUMX,DUMN0,DUMIA,DUMOUT,5,DUMF, GG(:,:,1),Z)
          Z=RLIMIT(SMALL,Z,BIG)
          T2=T
          P2=PF(T2,D,Z)
          DELP=P2-P1
          IF (DELP.EQ.0.0) GO TO 20
          DELT=(P-P2)/DELP*(T2-T1)
          T=T2+DELT
   10   END DO
        ERROR=BIT5
   20   CONTINUE
      END IF
      TEMP=T*TC

    ELSE IF (ENTRY.EQ.3) THEN !     CASE 3  ..GIVEN TEMPERATURE AND PRESSURE.
      write(dbg,*) 'entry=3, reduced t,d,p=', t,d,p
      VAPOR=.FALSE.
      IF (T.LT.1.0) THEN  !  get first guess for d from saturated value
                          !  set up for interpolation of saturation data at t
        CALL NTRP(SAT(:,1),T,NS,mIndexValue(1),KS,1,DUMF,DUMFF,DUMVAL)
        CALL NTRP(DUMA,DUMX,DUMN0,DUMIA,DUMOUT,4, SAT(:,2),DUMFF,PSAT)  !  PSAT
        GAS=P.LT.PSAT
        CALL NTRP(DUMA,DUMX,DUMN0,DUMIA,DUMOUT,4, SAT(:,3),DUMFF,DL) ! DL (liquid)
        CALL NTRP(DUMA,DUMX,DUMN0,DUMIA,DUMOUT,4, SAT(:,4),DUMFF,DG) ! DG (gas)
        IF (GAS) THEN
          D=DG
        ELSE
          D=DL
        END IF
        write(dbg,*) 'T<1', dl,dg,d
      ELSE
        GAS=.TRUE.
        Z=1.0
        CALL CUBIC (Z*T,P,ROOTS,NR) !     PROCEDURE (SOLVE CUBIC FOR D AT (T,P))
        write(dbg,*) 'cubic roots', roots
        IF (GAS) NR=1
        D=ROOTS(NR)
        write(dbg,*) d, ' chosen for D'
      END IF

      DELD=HUGE
      DO 30 I3=1,MAX3 !  solve for d at (t,p) by iteration)
        IF (ARG1.NE.T) THEN ! get z by double interpolation at (t,d))
          CALL NTRP(TT,T,N1,mIndexValue(3),K1,2,DUMF,DUMFF,DUMVAL)
          ARG1=T
        END IF
        IF (ARG2.NE.D) THEN
          CALL NTRP(DD,D,N2,mIndexValue(4),K2,3,DUMF,DUMFF,DUMVAL)
          ARG2=D
        END IF
        CALL NTRP(DUMA,DUMX,DUMN0,DUMIA,DUMOUT,5,DUMF, GG(:,:,1),Z)  ! gives Z
        Z=RLIMIT(SMALL,Z,BIG)
        IF (ABS(DELD).LE.C(6)*D) GO TO 40    !   exit the DO loop
        D1=D
        P1=PF(T,D1,Z)
        CALL CUBIC (Z*T,P,ROOTS,NR)   !  solve cubic for d at (t,p))
        IF (GAS) NR=1
        D=ROOTS(NR)
        D=RLIMIT(DD(1),D,DD(N2))   
        IF (ARG1.NE.T) THEN    !  get z by double interpolation at (t,d))
          CALL NTRP(TT,T,N1,mIndexValue(3),K1,2,DUMF,DUMFF,DUMVAL)
          ARG1=T
        END IF
        IF (ARG2.NE.D) THEN
          CALL NTRP(DD,D,N2,mIndexValue(4),K2,3,DUMF,DUMFF,DUMVAL)
          ARG2=D
        END IF
        CALL NTRP(DUMA,DUMX,DUMN0,DUMIA,DUMOUT,5,DUMF, GG(:,:,1),Z)
        Z=RLIMIT(SMALL,Z,BIG)
        D2=D
        P2=PF(T,D2,Z)
        DELP=P2-P1
        IF (DELP.EQ.0.0) GO TO 40 !              IF (DELP.EQ.0.0) EXIT (I3)
        DELD=(P-P2)/DELP*(D2-D1)
        D=D2+DELD
        write(dbg,*) 'end of iteration',i3, 'p1,p2', p1,p2
        write(dbg,*) '                       d1,d2', d1,d2
   30 END DO
      ERROR=BIT5   !     OTHERWISE
   40 CONTINUE
      DENS=D*DC
      write(dbg,*) 'final. after', i3, ' iterations, dens,d=', dens,d

    ELSE IF (ENTRY.EQ.4) THEN !     CASE 4   !.....GIVEN PRESSURE AND ENTROPY.
      write(dbg,*) 'entry=4, reduced t,d,p=', t,d,p
      SIN=PROPS(2)/R
      IF (P.LT.1.0) THEN    !     PROCEDURE (TEST FOR VAPOR CONDITION AT (P,SIN))
        CALL NTRP(SAT(:,2),MAX(P,SAT(1,2)),NS,mIndexValue(1),KS,1,DUMF,DUMFF,DUMVAL)
        CALL NTRP(DUMA,DUMX,DUMN0,DUMIA,DUMOUT,4, SAT(:,5),DUMFF,SL)  !     PROCEDURE (GET SATURATED LIQUID (SL) AND GAS (SG) ENTROPIES)
        CALL NTRP(DUMA,DUMX,DUMN0,DUMIA,DUMOUT,4, SAT(:,6),DUMFF,SG)
        VAPOR=SIN.GT.SL .AND. SIN.LT.SG
      ELSE
        VAPOR=.FALSE.
      END IF
      IF (VAPOR) THEN   !     PROCEDURE (GET (TWO-PHASE) T, WL, WG AND D AT (P,SIN))
        CALL NTRP(DUMA,DUMX,DUMN0,DUMIA,DUMOUT,4, SAT(:,1),DUMFF,TSAT)
        CALL NTRP(DUMA,DUMX,DUMN0,DUMIA,DUMOUT,4, SAT(:,3),DUMFF,DL)  !     PROCEDURE (GET SATURATED LIQUID (DL) AND GAS (DG) DENSITIES)
        CALL NTRP(DUMA,DUMX,DUMN0,DUMIA,DUMOUT,4, SAT(:,4),DUMFF,DG)
        T=TSAT
        WL=(SIN-SG)/(SL-SG)
        WG=1.0-WL
        D=1.0/(WL/DL+WG/DG)
      ELSE
        GAS=SIN.GT.SAT(NS,5)   !     PROCEDURE (GET FIRST GUESS FOR D AT (P,SIN))
        Z=1.0
        IF (P.LT.1.0) THEN
          CALL NTRP(DUMA,DUMX,DUMN0,DUMIA,DUMOUT,4, SAT(:,3),DUMFF,DL)  !     PROCEDURE (GET SATURATED LIQUID (DL) AND GAS (DG) DENSITIES)
          CALL NTRP(DUMA,DUMX,DUMN0,DUMIA,DUMOUT,4, SAT(:,4),DUMFF,DG)
          IF (GAS) THEN
            D=DG
            S1=SG
          ELSE
             D=DL
             S1=SL
          END IF
        ELSE
          T=TT(N1) !     PROCEDURE (ASSUME WORST-CASE STARTING VALUE FOR DENSITY)
          CALL CUBIC (Z*T,P,ROOTS,NR) !     PROCEDURE (SOLVE CUBIC FOR D AT (T,P))
          IF (GAS) NR=1
          D=ROOTS(NR)
          D=MAX(1.0,ROOTS(1))
          T=TF(P,D,Z) !     PROCEDURE (SOLVE FOR T AT (P,D) BY ITERATION)
          DELT=HUGE
          DO 50 I2=1,MAX2
            IF (ARG1.NE.T) THEN !     PROCEDURE (GET Z BY DOUBLE INTERPOLATION AT (T,D))
              CALL NTRP(TT,T,N1,mIndexValue(3),K1,2,DUMF,DUMFF,DUMVAL)
              ARG1=T
            END IF
            IF (ARG2.NE.D) THEN
              CALL NTRP(DD,D,N2,mIndexValue(4),K2,3,DUMF,DUMFF,DUMVAL)
              ARG2=D
            END IF
            CALL NTRP(DUMA,DUMX,DUMN0,DUMIA,DUMOUT,5,DUMF, GG(:,:,1),Z)
            Z=RLIMIT(SMALL,Z,BIG)
            IF (ABS(DELT).LE.C(5)*T) GO TO 60    !                    IF (ABS(DELT).LE.C(5)*T) EXIT (I2)
            T1=T
            P1=PF(T1,D,Z)
            T=RLIMIT(TT(1),TF(P,D,Z),TT(N1))
            IF (ARG1.NE.T) THEN !     PROCEDURE (GET Z BY DOUBLE INTERPOLATION AT (T,D))
!                       CALL NTRP1(TT,T,N1,L(3),K1,2)
              CALL NTRP(TT,T,N1,mIndexValue(3),K1,2,DUMF,DUMFF,DUMVAL)
              ARG1=T
            END IF
            IF (ARG2.NE.D) THEN
              CALL NTRP(DD,D,N2,mIndexValue(4),K2,3,DUMF,DUMFF,DUMVAL)
              ARG2=D
            END IF
            CALL NTRP(DUMA,DUMX,DUMN0,DUMIA,DUMOUT,5,DUMF, GG(:,:,1),Z)
            Z=RLIMIT(SMALL,Z,BIG)
            T2=T
            P2=PF(T2,D,Z)
            DELP=P2-P1
            IF (DELP.EQ.0.0) GO TO 60    !                    IF (DELP.EQ.0.0) EXIT (I2)
            DELT=(P-P2)/DELP*(T2-T1)
            T=T2+DELT
   50     END DO
          ERROR=BIT5    !                 OTHERWISE
   60     CONTINUE
          IF (ARG0.NE.T) THEN !     PROCEDURE (GET S0 BY INTERPOLATION OF ZERO-DENSITY DATA)
            CALL NTRP(G(:,1),T,NG,mIndexValue(2),K0,1,DUMF,DUMFF,DUMVAL)
            ARG0=T
          END IF
          CALL NTRP(DUMA,DUMX,DUMN0,DUMIA,DUMOUT,4, G(:,2),DUMFF,S0)
          IF (ARG1.NE.T) THEN !     PROCEDURE (GET S BY DOUBLE INTERPOLATION FOR Z2)
            CALL NTRP(TT,T,N1,mIndexValue(3),K1,2,DUMF,DUMFF,DUMVAL)
            ARG1=T
          END IF
          IF (ARG2.NE.D) THEN
            CALL NTRP(DD,D,N2,mIndexValue(4),K2,3,DUMF,DUMFF,DUMVAL)
            ARG2=D
          END IF
          CALL NTRP(DUMA,DUMX,DUMN0,DUMIA,DUMOUT,5,DUMF, GG(:,:,2),Z2)
          S=Z2*(S0+SF(D))
          S1=S
        END IF

        IF (GAS) THEN !     PROCEDURE (SOLVE FOR T AND D AT (P,SIN) BY ITERATION)
          X1=LOG(D)
          D=D/1.2
          X=LOG(D)
        ELSE
          X1=LOG(3.0-D)
          D=1.01*D
          X=LOG(3.0-D)
        END IF
        DELS=HUGE
        DO 90  I4=1,MAX4
          T=TF(P,D,Z) !     PROCEDURE (SOLVE FOR T AT (P,D) BY ITERATION)
          DELT=HUGE
          DO 70 I2=1,MAX2
            IF (ARG1.NE.T) THEN !     PROCEDURE (GET Z BY DOUBLE INTERPOLATION AT (T,D))
              CALL NTRP(TT,T,N1,mIndexValue(3),K1,2,DUMF,DUMFF,DUMVAL)
              ARG1=T
            END IF
            IF (ARG2.NE.D) THEN
              CALL NTRP(DD,D,N2,mIndexValue(4),K2,3,DUMF,DUMFF,DUMVAL)
              ARG2=D
            END IF
            CALL NTRP(DUMA,DUMX,DUMN0,DUMIA,DUMOUT,5,DUMF, GG(:,:,1),Z)
            Z=RLIMIT(SMALL,Z,BIG)
            IF (ABS(DELT).LE.C(5)*T) GO TO 80 !                    IF (ABS(DELT).LE.C(5)*T) EXIT (I2)
            T1=T
            P1=PF(T1,D,z)
            T=RLIMIT(TT(1),TF(P,D,z),TT(N1))
            IF (ARG1.NE.T) THEN   !     PROCEDURE (GET Z BY DOUBLE INTERPOLATION AT (T,D))
              CALL NTRP(TT,T,N1,mIndexValue(3),K1,2,DUMF,DUMFF,DUMVAL)
              ARG1=T
            END IF
            IF (ARG2.NE.D) THEN
              CALL NTRP(DD,D,N2,mIndexValue(4),K2,3,DUMF,DUMFF,DUMVAL)
              ARG2=D
            END IF
            CALL NTRP(DUMA,DUMX,DUMN0,DUMIA,DUMOUT,5,DUMF, GG(:,:,1),Z)
            Z=RLIMIT(SMALL,Z,BIG)
            T2=T
            P2=PF(T2,D,z)
            DELP=P2-P1
            IF (DELP.EQ.0.0) GO TO 80 !                    IF (DELP.EQ.0.0) EXIT (I2)
            DELT=(P-P2)/DELP*(T2-T1)
            T=T2+DELT
   70     END DO
          ERROR=BIT5 !                 OTHERWISE
   80     CONTINUE
          IF (ABS(DELS).LE.C(7)) GO TO 100 !                 IF (ABS(DELS).LE.C(7)) EXIT (I4)
          IF (X.NE.X1) THEN
          ERROR=0
          IF (ARG0.NE.T) THEN !     PROCEDURE (GET S0 BY INTERPOLATION OF ZERO-DENSITY DATA)
            CALL NTRP(G(:,1),T,NG,mIndexValue(2),K0,1,DUMF,DUMFF, DUMVAL)
            ARG0=T
          END IF
          CALL NTRP(DUMA,DUMX,DUMN0,DUMIA,DUMOUT,4, G(:,2),DUMFF,S0)
          IF (ARG1.NE.T) THEN !     PROCEDURE (GET S BY DOUBLE INTERPOLATION FOR Z2)
            CALL NTRP(TT,T,N1,mIndexValue(3),K1,2,DUMF,DUMFF,DUMVAL)
            ARG1=T
          END IF
          IF (ARG2.NE.D) THEN
            CALL NTRP(DD,D,N2,mIndexValue(4),K2,3,DUMF,DUMFF,DUMVAL)
            ARG2=D
          END IF
          CALL NTRP(DUMA,DUMX,DUMN0,DUMIA,DUMOUT,5,DUMF, GG(:,:,2),Z2)
          S=Z2*(S0+SF(D))
          DXDS=(X1-X)/(S1-S)
          DELS=SIN-S
          DELX=RLIMIT(-DXMAX,DXDS*DELS,DXMAX)
          X1=X
          S1=S
          X=X+DELX
          IF (GAS) THEN
            D=EXP(X)
          ELSE
            D=3.0-EXP(X)
          END IF
          DLIM=RLIMIT(DD(1),D,DD(N2))
          IF (DLIM.NE.D) THEN
            D=DLIM
            IF (GAS) THEN
              X=ALOG(D)
            ELSE
              X=ALOG(3.0-D)
            END IF
          END IF
        END IF
   90   END DO
        ERROR=BIT5 !              OTHERWISE
  100          CONTINUE
      END IF
      DENS=D*DC
      TEMP=T*TC

    ELSE IF (ENTRY.EQ.5) THEN !     CASE 5      !.....GIVEN PRESSURE AND ENTHALPY. ***************
      write(dbg,*) 'entry=5, reduced t,d,p=', t,d,p
      HIN=PROPS(3)/R/TC
      IF (P.LT.1.0) THEN    !     PROCEDURE (TEST FOR VAPOR CONDITION AT (P,HIN))
                            !     PROCEDURE (SET UP FOR INTERPOLATION OF SATURATION DATA AT P)
        CALL NTRP(SAT(:,2),AMAX1(P,SAT(1,2)),NS,mIndexValue(1),KS,1, DUMF,DUMFF,DUMVAL)
        CALL NTRP(DUMA,DUMX,DUMN0,DUMIA,DUMOUT,4, SAT(:,7),DUMFF,HL)  !     PROCEDURE (GET SATURATED LIQUID (HL) AND GAS (HG) ENTHALPIES)
        CALL NTRP(DUMA,DUMX,DUMN0,DUMIA,DUMOUT,4, SAT(:,8),DUMFF,HG)
        VAPOR=HIN.GT.HL .AND. HIN.LT.HG
      ELSE
        VAPOR=.FALSE.
      END IF
      IF (VAPOR) THEN    !     PROCEDURE (GET (TWO-PHASE) T, WL, WG AND D AT (P,HIN))
                         ! get saturation temperature tsat, 
                         ! and then dl and dg, saturated liguid and gas densities
        CALL NTRP(DUMA,DUMX,DUMN0,DUMIA,DUMOUT,4, SAT(:,1),DUMFF,TSAT)
        CALL NTRP(DUMA,DUMX,DUMN0,DUMIA,DUMOUT,4, SAT(:,3),DUMFF,DL)  
        CALL NTRP(DUMA,DUMX,DUMN0,DUMIA,DUMOUT,4, SAT(:,4),DUMFF,DG)  
        T=TSAT
        WL=(HIN-HG)/(HL-HG)
        WG=1.0-WL
        D=1.0/(WL/DL+WG/DG)
      ELSE
        Z=1.0    !     PROCEDURE (GET FIRST GUESS FOR D AT (P,HIN))
        IF (P.LT.1.0) THEN
          GAS=HG.LT.HIN
          CALL NTRP(DUMA,DUMX,DUMN0,DUMIA,DUMOUT,4, SAT(:,3),DUMFF,DL)  ! dl
          CALL NTRP(DUMA,DUMX,DUMN0,DUMIA,DUMOUT,4, SAT(:,4),DUMFF,DG)  ! dg
          IF (GAS) THEN
            D=DG
            H1=HG
          ELSE
            D=DL
            H1=HL
          END IF
        ELSE
          T=TT(N1) ! assume worst-case starting value for density
          CALL CUBIC (Z*T,P,ROOTS,NR) !  solve cubic for d at (t,p)
          IF (GAS) NR=1
          D=ROOTS(NR)
          D=AMAX1(1.0,ROOTS(1))
          T=TF(P,D,z) !   solve for t at (p,d) by iteration
          DELT=HUGE
          DO 110 I2=1,MAX2
            IF (ARG1.NE.T) THEN !  get z by double interpolation at (t,d)
              CALL NTRP(TT,T,N1,mIndexValue(3),K1,2,DUMF,DUMFF,rDUMOUT)
              ARG1=T
            END IF
            IF (ARG2.NE.D) THEN
              CALL NTRP(DD,D,N2,mIndexValue(4),K2,3,DUMF,DUMFF,rDUMOUT)
              ARG2=D
            END IF
            CALL NTRP(DUMA,DUMX,DUMN0,DUMIA,DUMOUT,5,DUMF, GG(:,:,1),Z)
            Z=RLIMIT(SMALL,Z,BIG)
            IF (ABS(DELT).LE.C(5)*T) GO TO 120 !   IF (ABS(DELT).LE.C(5)*T) EXIT (I2)
            T1=T
            P1=PF(T1,D,z)
            T=RLIMIT(TT(1),TF(P,D,z),TT(N1))
            IF (ARG1.NE.T) THEN ! get z by double interpolation at (t,d))
              CALL NTRP(TT,T,N1,mIndexValue(3),K1,2,DUMF,DUMFF,DUMVAL)
              ARG1=T
            END IF
            IF (ARG2.NE.D) THEN
              CALL NTRP(DD,D,N2,mIndexValue(4),K2,3,DUMF,DUMFF,DUMVAL)
              ARG2=D
            END IF
            CALL NTRP(DUMA,DUMX,DUMN0,DUMIA,DUMOUT,5,DUMF, GG(:,:,1),Z)
            Z=RLIMIT(SMALL,Z,BIG)
            T2=T
            P2=PF(T2,D,z)
            DELP=P2-P1
            IF (DELP.EQ.0.0) GO TO 120 !     IF (DELP.EQ.0.0) EXIT (I2)
            DELT=(P-P2)/DELP*(T2-T1)
            T=T2+DELT
  110     END DO
          ERROR=BIT5 !                 OTHERWISE
  120     CONTINUE
          IF (ARG0.NE.T) THEN ! get u0 by interpolation of zero-density data)
            CALL NTRP(G(:,1),T,NG,mIndexValue(2),K0,1,DUMF,DUMFF,DUMVAL)
            ARG0=T
          END IF
          CALL NTRP(DUMA,DUMX,DUMN0,DUMIA,DUMOUT,4, G(:,3),DUMFF,U0)
          IF (ARG1.NE.T) THEN ! get h by double interpolation for z3)
            CALL NTRP(TT,T,N1,mIndexValue(3),K1,2,DUMF,DUMFF,DUMVAL)
            ARG1=T
          END IF
          IF (ARG2.NE.D) THEN
            CALL NTRP(DD,D,N2,mIndexValue(4),K2,3,DUMF,DUMFF,DUMVAL)
            ARG2=D
          END IF
          CALL NTRP(DUMA,DUMX,DUMN0,DUMIA,DUMOUT,5,DUMF, GG(:,:,3),Z3)
          H=Z3*(U0+HF(P,D))
          H1=H
          GAS=H1.LT.HIN
        END IF

        IF (GAS) THEN ! solve for t and d at (p,hin) by iteration)
          X1=1.0/D
          D=D/1.2
          X=1.0/D
        ELSE
          X1=D
          D=1.01*D
          X=D
        END IF
        DELH=HUGE

        DO 150 I5=1,MAX5    !  solve for t at (p,d) by iteration)
          T=TF(P,D,z)
          DELT=HUGE
          DO 130 I2=1,MAX2
            IF (ARG1.NE.T) THEN    !  get z by double interpolation at (t,d))
              CALL NTRP(TT,T,N1,mIndexValue(3),K1,2,DUMF,DUMFF,DUMVAL)
              ARG1=T
            END IF
            IF (ARG2.NE.D) THEN
              CALL NTRP(DD,D,N2,mIndexValue(4),K2,3,DUMF,DUMFF,DUMVAL)
              ARG2=D
            END IF
            CALL NTRP(DUMA,DUMX,DUMN0,DUMIA,DUMOUT,5,DUMF, GG(:,:,1),Z)
            Z=RLIMIT(SMALL,Z,BIG)
            IF (ABS(DELT).LE.C(5)*T) GO TO 140    !    IF (ABS(DELT).LE.C(5)*T) EXIT (I2)
            T1=T
            P1=PF(T1,D,z)
            T=RLIMIT(TT(1),TF(P,D,z),TT(N1))
            IF (ARG1.NE.T) THEN   !     PROCEDURE (GET Z BY DOUBLE INTERPOLATION AT (T,D))
              CALL NTRP(TT,T,N1,mIndexValue(3),K1,2,DUMF,DUMFF,rDUMOUT)
              ARG1=T
            END IF
            IF (ARG2.NE.D) THEN
              CALL NTRP(DD,D,N2,mIndexValue(4),K2,3,DUMF,DUMFF,rDUMOUT)
              ARG2=D
            END IF
            CALL NTRP(DUMA,DUMX,DUMN0,DUMIA,DUMOUT,5,DUMF, GG(:,:,1),Z)
            Z=RLIMIT(SMALL,Z,BIG)
            T2=T
            P2=PF(T2,D,z)
            DELP=P2-P1
            IF (DELP.EQ.0.0) GO TO 140    !  IF (DELP.EQ.0.0) EXIT (I2)
            DELT=(P-P2)/DELP*(T2-T1)
            T=T2+DELT
  130     END DO
          ERROR=BIT5   !                 OTHERWISE
  140     CONTINUE

          IF (ABS(DELH).LE.C(8)) GO TO 160   !                 IF (ABS(DELH).LE.C(8)) EXIT (I5)
          IF (X.NE.X1) THEN
            ERROR=0
            IF (ARG0.NE.T) THEN   ! set up for interpolation of zero-density data at t
              CALL NTRP(G(:,1),T,NG,mIndexValue(2),K0,1,DUMF,DUMFF, DUMVAL)
              ARG0=T
            END IF
            CALL NTRP(DUMA,DUMX,DUMN0,DUMIA,DUMOUT,4, G(:,3),DUMFF,U0)
            IF (ARG1.NE.T) THEN   !     PROCEDURE (GET H BY DOUBLE INTERPOLATION FOR Z3)
              CALL NTRP(TT,T,N1,mIndexValue(3),K1,2,DUMF,DUMFF,DUMVAL)
              ARG1=T
            END IF
            IF (ARG2.NE.D) THEN
              CALL NTRP(DD,D,N2,mIndexValue(4),K2,3,DUMF,DUMFF,DUMVAL)
              ARG2=D
            END IF
            CALL NTRP(DUMA,DUMX,DUMN0,DUMIA,DUMOUT,5,DUMF, GG(:,:,3),Z3)
            H=Z3*(U0+HF(P,D))
            DXDH=(X1-X)/(H1-H)
            DELH=HIN-H
            DELX=DXDH*DELH
            X1=X
            H1=H
            X=X+DELX
            IF (GAS) THEN
              D=1.0/X
            ELSE
              D=X
            END IF
            DLIM=RLIMIT(DD(1),D,DD(N2))
            IF (D.NE.DLIM) THEN
              D=DLIM
              IF (GAS) THEN
                X=1.0/D
              ELSE
                X=D
              END IF
            END IF
          END IF

  150   END DO
        ERROR=BIT5   !              OTHERWISE
  160   CONTINUE
            END IF
            DENS=D*DC
            TEMP=T*TC

      END IF

! This concludes the initial calculations, solving for the thermodynamic state.
! Next, we compute the remaining properties at this condition.

      IF (T.GT.0.0 .AND. D.GT.0.0) THEN   !.....calculate remaining properties at (t,d).
        IF (VAPOR) THEN
          NPROPS=MIN0(3,NP)   !  calculate remaining properties for 2-phase fluid)
          DO 170 IPROP=1,NPROPS
            IF (IPROP.EQ.1) THEN    !     CASE 1  ***
              IF (ENTRY.LE.3) THEN
                WL=DL/D*(D-DG)/(DL-DG)
                WG=1.0-WL
              END IF
              DENSL=DL*DC
              DENSG=DG*DC
              PROPS(1)=PRES/(R*DENS*TEMP)
            ELSE IF (IPROP.EQ.2) THEN   !     CASE 2  ***
              IF (ENTRY.NE.4) THEN   !  get saturated liquid (sl) and gas (sg) entropies
                CALL NTRP(DUMA,DUMX,DUMN0,DUMIA,DUMOUT,4, SAT(:,5),DUMFF,SL)
                CALL NTRP(DUMA,DUMX,DUMN0,DUMIA,DUMOUT,4, SAT(:,6),DUMFF,SG)
                PROPS(2)=(WL*SL+WG*SG)*R
              END IF
              ENTL=SL*R
              ENTG=SG*R
            ELSE IF (IPROP.EQ.3) THEN !     CASE 3  ***
              IF (ENTRY.NE.5) THEN   ! get saturated liquid (hl) and gas (hg) enthalpies
                CALL NTRP(DUMA,DUMX,DUMN0,DUMIA,DUMOUT,4, SAT(:,7),DUMFF,HL)
                CALL NTRP(DUMA,DUMX,DUMN0,DUMIA,DUMOUT,4, SAT(:,8),DUMFF,HG)
                PROPS(3)=(WL*HL+WG*HG)*R*TC
              END IF
              ENTHL=HL*R*TC
              ENTHG=HG*R*TC
            END IF
  170     END DO
        ELSE
          NEXTP=1 !   calculate remaining properties for homogeneous fluid)  (vapor=.FALSE.)
          DO 180 IPROP=NEXTP,NP
            IF (IPROP.LE.6) THEN
              IF (IPROP.EQ.1) THEN   !     CASE 1  ***
                PROPS(1)=PRES/(R*DENS*TEMP)
              ELSE IF (IPROP.EQ.2) THEN !     CASE 2 ***
                IF (ENTRY.NE.4) THEN  !  set up for interpolation of zero-density data at t
                  IF (ARG0.NE.T) THEN
                    CALL NTRP(G(:,1),T,NG,mIndexValue(2),K0,1,DUMF,DUMFF,DUMVAL)
                    ARG0=T
                  END IF
                  CALL NTRP(DUMA,DUMX,DUMN0,DUMIA,DUMOUT,4,G(:,2),DUMFF,S0)
                  IF (ARG1.NE.T) THEN  ! get s by double interpolation for z2
                    CALL NTRP(TT,T,N1,mIndexValue(3),K1,2,DUMF, DUMFF,DUMVAL)
                    ARG1=T
                  END IF
                  IF (ARG2.NE.D) THEN
                    CALL NTRP(DD,D,N2,mIndexValue(4),K2,3,DUMF, DUMFF,DUMVAL)
                    ARG2=D
                  END IF
                  CALL NTRP(DUMA,DUMX,DUMN0,DUMIA,DUMOUT,5,DUMF,GG(:,:,2),Z2)
                  S=Z2*(S0+SF(D))
                  PROPS(2)=S*R
                END IF
              ELSE IF (IPROP.EQ.3) THEN   !     CASE 3  ***
                IF (ENTRY.NE.5) THEN  ! set up for interpolation of zero-density data at t
                  IF (ARG0.NE.T) THEN
                    CALL NTRP(G(:,1),T,NG,mIndexValue(2),K0,1, DUMF,DUMFF,DUMVAL)
                    ARG0=T
                  END IF
                  CALL NTRP(DUMA,DUMX,DUMN0,DUMIA,DUMOUT,4, G(:,3),DUMFF,U0)
                  IF (ARG1.NE.T) THEN !  get h by double interpolation for z3
                    CALL NTRP(TT,T,N1,mIndexValue(3),K1,2,  DUMF,DUMFF,DUMVAL)
                    ARG1=T
                  END IF
                  IF (ARG2.NE.D) THEN
                    CALL NTRP(DD,D,N2,mIndexValue(4),K2,3, DUMF,DUMFF,DUMVAL)
                    ARG2=D
                  END IF
                  CALL NTRP(DUMA,DUMX,DUMN0,DUMIA,DUMOUT,5, DUMF,GG(:,:,3),Z3)
                  H=Z3*(U0+HF(P,D))
                  PROPS(3)=H*R*TC
                END IF
              ELSE IF (IPROP.EQ.4) THEN   !     CASE 4   ***
                CALL NTRP(DUMA,DUMX,DUMN0,DUMIA,DUMOUT,4, G(:,4),DUMFF,CV0) !  get cv0
                CALL NTRP(DUMA,DUMX,DUMN0,DUMIA,DUMOUT,5, DUMF,GG(:,:,4),Z4) ! Z4
                CV=Z4*CV0
                PROPS(4)=CV*R
              ELSE IF (IPROP.EQ.5) THEN   !     CASE 5  ***
                CALL NTRP(DUMA,DUMX,DUMN0,DUMIA,DUMOUT,5, DUMF,GG(:,:,5),Z5) !  Z5
                CP=Z5*(CV+CPF(T,D))
                PROPS(5)=CP*R
                GAMMA=CP/CV
              ELSE IF (IPROP.EQ.6) THEN !     CASE 6  ***
                CALL NTRP(DUMA,DUMX,DUMN0,DUMIA,DUMOUT,5,DUMF,GG(:,:,6),Z6) !  Z6
                grt=gamma*R*tc
                A2=A2F(T,D,grt)
                IF (A2.GT.0.0) THEN
                  A=Z6*SQRT(A2)
                ELSE
                  A=0.0
                END IF
                PROPS(6)=A
              END IF
            ELSE    !  DOUBLE INTERPOLATION FOR OTHER PROPERTIES)
              CALL NTRP(DUMA,DUMX,DUMN0,DUMIA,DUMOUT,5,DUMF,GG(:,:,IPROP),PROPS(IPROP))
            END IF
  180     END DO
        END IF
      END IF
      ERROR=KS*BIT1+K0*BIT2+K1*BIT3+K2*BIT4+ERROR
!     RETURN

    ELSE IF (IGOTO.EQ.2.OR.IGOTO.EQ.3) THEN
      IF (IGOTO.EQ.2) THEN   !.....ENTRY FOR GETTING LIQUID-PHASE PROPERTIES NOT IN COMMON.
        D=DL
!     RETURN
      ELSE IF (IGOTO.EQ.3) THEN   !.....ENTRY FOR GETTING GAS-PHASE PROPERTIES NOT IN COMMON.
        D=DG
      END IF

      IF (.NOT.VAPOR) THEN   !   calculate other saturation-locus properties)
        WRITE(*,*) "Fatal Error in Fluid. Vapor is FALSE"
        STOP "Fatal Error #3 in Fluid"
      END IF
      IF (NP.GT.N3) THEN
        WRITE(*,*) "Fatal error in Fluid. np > n3 ", NP,N3
        STOP "Fatal error #4 in Fluid"
      END IF
      IF (ARG0.NE.T) THEN   !  set up for interpolation of zero-density data at t
        CALL NTRP(G(:,1),T,NG,mIndexValue(2),K0,1,DUMF,DUMFF,DUMVAL)
        ARG0=T
      END IF
      IF (ARG1.NE.T) THEN   ! set up for double interpolation at (t,d))
        CALL NTRP(TT,T,N1,mIndexValue(3),K1,2,DUMF,DUMFF,DUMVAL)
        ARG1=T
      END IF
      IF (ARG2.NE.D) THEN
        CALL NTRP(DD,D,N2,mIndexValue(4),K2,3,DUMF,DUMFF,DUMVAL)
        ARG2=D
      END IF
      IF (T.GT.0.0 .AND. NP.GE.4) THEN
        NEXTP=4
        DO 190 IPROP=NEXTP,NP   ! calculate remaining properties for homogeneous fluid
          IF (IPROP.LE.6) THEN
            IF (IPROP.EQ.1) THEN !     CASE 1  ***
              PROPS(1)=PRES/(R*DENS*TEMP)
            ELSE IF (IPROP.EQ.2) THEN !     CASE 2 ***
              IF (ENTRY.NE.4) THEN
                IF (ARG0.NE.T) THEN ! get s0 by interpolation of zero-density data at t
                  CALL NTRP(G(:,1),T,NG,mIndexValue(2),K0,1,DUMF,DUMFF,DUMVAL)
                  ARG0=T
                END IF
                CALL NTRP(DUMA,DUMX,DUMN0,DUMIA,DUMOUT,4, G(:,2),DUMFF,S0)
                IF (ARG1.NE.T) THEN  !  get s by double interpolation for z2 at (t,d)
                  CALL NTRP(TT,T,N1,mIndexValue(3),K1,2,DUMF,DUMFF, DUMVAL)
                  ARG1=T
                END IF
                IF (ARG2.NE.D) THEN
                  CALL NTRP(DD,D,N2,mIndexValue(4),K2,3,DUMF,DUMFF, DUMVAL)
                  ARG2=D
                END IF
                CALL NTRP(DUMA,DUMX,DUMN0,DUMIA,DUMOUT,5, DUMF,GG(:,:,2),Z2)
                S=Z2*(S0+SF(D))
                PROPS(2)=S*R
              END IF
            ELSE IF (IPROP.EQ.3) THEN  !     CASE 3 ***
              IF (ENTRY.NE.5) THEN
                IF (ARG0.NE.T) THEN  ! get u0 by interpolation of zero-density data at t
                  CALL NTRP(G(:,1),T,NG,mIndexValue(2),K0,1,DUMF, DUMFF,DUMVAL)
                  ARG0=T
                END IF
                CALL NTRP(DUMA,DUMX,DUMN0,DUMIA,DUMOUT,4,G(:,3),DUMFF,U0)
                IF (ARG1.NE.T) THEN  ! get h by double interpolation for z3 at (t,d)
                  CALL NTRP(TT,T,N1,mIndexValue(3),K1,2,DUMF,DUMFF, DUMVAL)
                  ARG1=T
                END IF
                IF (ARG2.NE.D) THEN
                  CALL NTRP(DD,D,N2,mIndexValue(4),K2,3,DUMF,DUMFF, DUMVAL)
                  ARG2=D
                END IF
                CALL NTRP(DUMA,DUMX,DUMN0,DUMIA,DUMOUT,5,DUMF,GG(:,:,3),Z3)
                H=Z3*(U0+HF(P,D))
                PROPS(3)=H*R*TC
              END IF
            ELSE IF (IPROP.EQ.4) THEN  !     CASE 4 ***
              CALL NTRP(DUMA,DUMX,DUMN0,DUMIA,DUMOUT,4, G(:,4),DUMFF,CV0)  !   cv0
              CALL NTRP(DUMA,DUMX,DUMN0,DUMIA,DUMOUT,5, DUMF,GG(:,:,4),Z4)  !  Z4
              CV=Z4*CV0
              PROPS(4)=CV*R
            ELSE IF (IPROP.EQ.5) THEN   !     CASE 5 **
              CALL NTRP(DUMA,DUMX,DUMN0,DUMIA,DUMOUT,5, DUMF,GG(:,:,5),Z5) !  Z5
              CP=Z5*(CV+CPF(T,D))
              PROPS(5)=CP*R
              GAMMA=CP/CV
            ELSE IF (IPROP.EQ.6) THEN   !     CASE 6 **
              CALL NTRP(DUMA,DUMX,DUMN0,DUMIA,DUMOUT,5,DUMF,GG(:,:,6),Z6)   !  Z6
              grt=gamma*R*tc
              A2=A2F(T,D,grt)
              IF (A2.GT.0.0) THEN
                A=Z6*SQRT(A2)
              ELSE
                A=0.0
              END IF
              PROPS(6)=A
            END IF
          ELSE
            CALL NTRP(DUMA,DUMX,DUMN0,DUMIA,DUMOUT,5, DUMF,GG(:,:,IPROP),PROPS(IPROP))         !     PROCEDURE (DOUBLE INTERPOLATION FOR OTHER PROPERTIES)
          END IF
  190   END DO
      END IF
      ERROR=KS*BIT1+K0*BIT2+K1*BIT3+K2*BIT4
!     RETURN

    ELSE
      WRITE(*,*) "FATAL ERROR IN FLUID: BAD IGOTO ", IGOTO
    END IF
!!!    1 FORMAT (' FATAL ERROR IN FLUID: NP > N3',' NP =',I4,' N3=',I4)
!!!    2 FORMAT (' FATAL ERROR IN FLUID: BAD ENTRY')
!!!    3 FORMAT (' FATAL ERROR IN FLUID: VAPOR=.FALSE.')

  RETURN
END Subroutine Fluid   ! ====================================================




END Module FluidProcedures   ! =================================================


