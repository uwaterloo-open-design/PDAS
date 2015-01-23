!+                                                                      
! MODULE Rational-Spline Approximation with Automatic Tension Adjustment
! ------------------------------------------------------------------------------
! PURPOSE - Coding of the algorithms described in NASA Technical Paper 2366
!  by James R. Schiess and Patricia A. Kerr


! AUTHORS - James R. Schiess and Patricia A. Kerr, NASA Langley
!         - Ralph L. Carmichael, Public Domain Aeronautical Software

!     REVISION HISTORY
!   DATE  VERS PERSON  STATEMENT OF CHANGES
!   1984   0.1 JRS&PAK Original coding of program RASPFIT at NASA Langley
!   1988   0.2 JRS&PAK Release of program LAR-13694 to COSMIC
!   1996   0.3   RLC   Acquisition of program LAR-13694 from COSMIC
! 30Dec01  0.4   RLC   Converted to Fortran 90 free-format
! 22Jul08  0.5   RLC   Moved procedures into module.

! Scientific data often contains random errors that make plotting and
! curve-fitting difficult. The Rational-Spline Approximation with 
! Automatic Tension Adjustment algorithm lead to a flexible, smooth 
! representation of experimental data.

! The user sets the conditions for each consecutive pair of knots: (knots
! are user-defined divisions in the data set) to apply no tension; to 
! apply fixed tension; or to determine tension with a tension adjustment 
! algorithm. The user also selects the number of knots, the knot 
! abscissas, and the allowed maximum deviations from line segments. The 
! selection of these quantities depends on the actual data and on the 
! requirements of a particular application. This program differs from 
! the usual spline under tension in that it allows the user to specify 
! different tension values between each adjacent pair of knots rather 
! than a constant tension over the entire data range. The subroutines use 
! an automatic adjustment scheme that varies the tension parameter for 
! each interval until the maximum deviation of the spline from the line 
! joining the knots is less than or equal to a user-specified amount. 
! This procedure frees the user from the drudgery of adjusting individual 
! tension parameters while still giving control over the local behavior 
! of the spline.






!+
      MODULE RationalSpline
! ------------------------------------------------------------------------------
! PURPOSE - RATIONAL SPLINE SUBROUTINES
!  ( NASA Langley Research Center )



      PRIVATE:: Qxz019
      PRIVATE:: Qxz063
      PRIVATE:: Qxz255
      PRIVATE:: Qxz275
      PRIVATE:: Qxz276
      PRIVATE:: Qxz277
      PRIVATE:: Qxz278
      PRIVATE:: Qxz279
      PRIVATE:: Qxz280
      PRIVATE:: Qxz281
      PRIVATE:: qxz282
      PUBLIC:: Rsdsdr
      PUBLIC:: Rsutds
      PRIVATE:: Secbi
      PRIVATE:: Sympds


      CONTAINS


      SUBROUTINE RSUTDS(N,NK,X,Y,DF,XK,IPK,PK,DEV,NN,IW,                &
     & NKMAX,ITMAX,XR,R,C,WK,IERR)
!***********************************************************************
!*                                                                     *
!*    PURPOSE:                                                         *
!*                 SUBROUTINE RSUTDS FITS A SMOOTH RATIONAL SPLINE TO A*
!*                 UNIVARIATE FUNCTION.  DATA MAY BE UNEQUALLY SPACED. *
!*                                                                     *
!     E3.4                                                             *
!*                                                                     *
!*        USE:                                                         *
!*                 CALL RSUTDS(N,NK,X,Y,DF,XK,IPK,PK,DEV,NN,           *
!*                             IW,NKMAX,ITMAX,XR,R,C,WK,IERR)          *
!*                                                                     *
!*         N       AN INPUT INTEGER SPECIFYING THE NUMBER OF DATA      *
!*                 POINTS FOR THE INDEPENDENT VARIABLE.                *
!*                                                                     *
!*         NK      AN INPUT INTEGER SPECIFYING THE NUMBER OF KNOTS.    *
!*                 NK>=3.                                              *
!*                                                                     *
!*         X       AN INPUT ONE-DIMENSIONAL REAL ARRAY DIMENSIONED AT  *
!*                 AT LEAST N IN THE CALLING PROGRAM.  X CONTAINS THE  *
!*                 X-COORDINATES OF THE DATA POINTS.                   *
!*                                                                     *
!*         Y       AN INPUT ONE-DIMENSIONAL REAL ARRAY WHICH           *
!*                 CONTAINS THE FUNCTION VALUES OF THE DATA POINTS.    *
!*                 Y MUST BE DIMENSIONED AT LEAST N.                   *
!*                                                                     *
!*         DF      AN INPUT ONE-DIMENSIONAL REAL ARRAY WHICH CONTAINS  *
!*                 A SET OF POSITIVE WEIGHTS INDICATING THE STANDARD   *
!*                 DEVIATION OF THE ERROR OF Y(I).  IF THE STANDARD    *
!*                 DEVIATION IS NOT KNOWN, A SUGGESTED VALUE IS        *
!*                 DF(I)=1.0E-4, I=1,...,N.  DF MUST BE DIMENSIONED    *
!*                 AT LEAST N.                                         *
!*                                                                     *
!*         XK      AN INPUT ONE-DIMENSIONAL REAL ARRAY WHICH           *
!*                 CONTAINS THE X-COORDINATES OF THE SELECTED KNOTS.   *
!*                 XK MUST BE DIMENSIONED AT LEAST NK.  THE FIRST      *
!*                 AND LAST KNOTS MUST COINCIDE WITH X(1) AND X(N),    *
!*                 RESPECTIVELY, WHILE THE INTERVENING KNOTS MUST BE   *
!*                 CHOSEN TO LIE BETWEEN X(1) AND X(N),                *
!*                 IE, X(1)=XK(1)  XK(2)...XK(NK-1) XK(NK)=X(N).       *
!*                 KNOTS MUST BE CHOSEN SO THAT THERE ARE AT LEAST     *
!*                 THREE DATA POINTS IN EACH INTERVAL DEFINED BY       *
!*                 CONSECUTIVE PAIRS OF KNOTS.                         *
!*                                                                     *
!*         IPK     AN INPUT ONE-DIMENSIONAL INTEGER ARRAY WHICH        *
!*                 CONTAINS THE TENSION FACTOR ADJUSTMENT FLAGS.       *
!*                 IPK MUST BE DIMENSIONED AT LEAST NK-1.              *
!*                                                                     *
!*                 IPK(K)=0  FACTOR K MAY BE ADJUSTED.                 *
!*                                                                     *
!*                 IPK(K)=1  FACTOR K IS HELD CONSTANT.                *
!*                                                                     *
!*         PK      AN INPUT/OUTPUT ONE-DIMENSIONAL REAL ARRAY WHICH    *
!*                 CONTAINS THE TENSION FACTORS APPLIED ON THE         *
!*                 INTERVALS BETWEEN KNOTS.  PK MUST BE DIMENSIONED AT *
!*                 LEAST NK-1.  PK(K) MUST BE GREATER THAN -1. ON ENTRY*
!*                 IF IPK(K)=1, PK(K) WILL CONTAIN THE FINAL TENSION   *
!*                 FACTOR TO BE USED IN INTERVAL K.  IF IPK(K)=0, PK(K)*
!*                 MUST CONTAIN AN INITIAL VALUE OF TENSION FACTOR FOR *
!*                 INTERVAL K.  A VALUE OF 0.0 IS SUGGESTED.  ON       *
!*                 RETURN, PK(K) WILL CONTAIN THE FINAL TENSION FACTOR *
!*                 FOR INTERVAL K.                                     *
!*                                                                     *
!*         DEV     AN INPUT ONE-DIMENSIONAL REAL ARRAY WHICH           *
!*                 CONTAINS THE ALLOWED PERCENTAGE OF DEVIATION OF     *
!*                 THE APPROXIMATED SOLUTION FROM THE LINEAR SOLUTION  *
!*                 IN THE INTERVAL BETWEEN KNOT I AND KNOT I+1, FOR    *
!*                 I=1,...,NK-1. DEV MUST BE DIMENSIONED AT LEAST NK-1.*
!*                                                                     *
!*         NN      AN INPUT INTEGER WHOSE ABSOLUTE VALUE SPECIFIES     *
!*                 THE NUMBER OF POINTS AT WHICH INTERPOLATION         *
!*                 IS DESIRED.  IF NN IS POSITIVE, THE POINTS AT WHICH *
!*                 INTERPOLATION IS DESIRED MAY BE UNEQUALLY SPACED.   *
!*                 IF NN IS NEGATIVE, THE POINTS WILL BE EQUALLY SPACED*
!*                 BETWEEN XR(1) AND XR(ABS(NN)).  NN CANNOT BE -1     *
!*                 OR ZERO.                                            *
!*                                                                     *
!*                                                                     *
!*         IW      AN INPUT/OUTPUT INTEGER.                            *
!*                                                                     *
!*                 INPUT   ON THE FIRST CALL TO RSUTDS , IW MUST BE    *
!*                         SET TO 0.  ON SUBSEQUENT CALLS, IW SHOULD BE*
!*                         SET TO 1 TO PROVIDE ADDITIONAL INTERPOLATION*
!*                         ON THE SAME FITTING CALCULATIONS.           *
!*                                                                     *
!*                 OUTPUT  IW IS AUTOMATICALLY SET TO 1.               *
!*                                                                     *
!*                                                                     *
!*                                                                     *
!*         NKMAX   AN INPUT INTEGER WHICH AGREES EXACTLY WITH THE      *
!*                 FIRST DIMENSION OF ARRAY C. NKMAX MUST BE GREATER   *
!*                 THAN OR EQUAL TO 2*NK.                              *
!*                                                                     *
!*         ITMAX   AN INPUT INTEGER SPECIFYING THE MAXIMUM NUMBER      *
!*                 OF ITERATIONS PERMITTED IN ACHIEVING THE ALLOWED    *
!*                 DEVIATIONS.   ITMAX<=35.                            *
!*                                                                     *
!*         XR      AN INPUT/OUTPUT  ONE-DIMENSIONAL REAL ARRAY         *
!*                 DIMENSIONED AT LEAST ABS(NN), WHICH CONTAINS        *
!*                 POINTS TO BE INTERPOLATED.                          *
!*
!*                 INPUT    IF NN IS POSITIVE, THE FIRST NN ELEMENTS   *
!*                          OF XR MUST ALL BE DEFINED.  IF NN IS       *
!*                          NEGATIVE, ONLY ELEMENTS XR(1) AND          *
!*                          XR(ABS(NN)) NEED BE DEFINED.               *
!*                                                                     *
!*                 OUTPUT   IF NN IS POSITIVE, XR IS UNCHANGED.  IF NN *
!*                          IS NEGATIVE, THE FIRST ABS(NN) ELEMENTS OF *
!*                          XR CONTAIN EQUALLY SPACED INTERPOLATION    *
!*                          POINTS.                                    *
!*                                                                     *
!*         R       AN OUTPUT ONE-DIMENSIONAL REAL ARRAY WHICH CONTAINS *
!*                 THE VALUES OF THE INTERPOLATED DEPENDENT VARIABLE   *
!*                 CORRESPONDING TO XR.  R MUST BE DIMENSIONED AT LEAST*
!*                 ABS(NN).                                            *
!*                                                                     *
!*         C       AN OUTPUT TWO-DIMENSIONAL REAL ARRAY WHICH          *
!*                 CONTAINS THE CORRELATION MATRIX.  C MUST  HAVE      *
!*                 A FIRST DIMENSION OF NKMAX AND A SECOND DIMENSION   *
!*                 OF AT LEAST 2*NK.                                   *
!*                                                                     *
!*         WK      AN INPUT/OUTPUT ONE-DIMENSIONAL REAL ARRAY WHICH    *
!*                 MUST BE DIMENSIONED AT LEAST 16*NK**2-6*NK+5.       *
!*                                                                     *
!*                 INPUT   IF IW = 0, WK NEED NOT CONTAIN VALUES.      *
!*                                                                     *
!*                         IF IW = 1, WK SHOULD CONTAIN THE VALUES     *
!*                         WHICH WERE STORED IN WK BY A PREVIOUS CALL  *
!*                         TO RSUTDS.                                  *
!*                                                                     *
!*                 OUTPUT  WK CONTAINS INFORMATION NECESSARY FOR       *
!*                         SUCCEEDING CALLS TO RSUTDS USING THE SAME   *
!*                         INPUT DATA WITHOUT REINITIALIZATION.        *
!*                                                                     *
!*        IERR     OUTPUT ERROR PARAMETER:                             *
!*                 = 0  NORMAL RETURN.  NO ERROR DETECTED.             *
!*                 =-1  AN INTERVAL CONTAINS LESS THAN THREE DATA      *
!*                      POINTS.                                        *
!*                 =-2  DID NOT CONVERGE IN ALLOWED NUMBER OF ITERAT-  *
!*                      IONS, INDICATING THAT TENSION MAY NOT HAVE     *
!*                      BEEN ADJUSTED IN ONE OR MORE INTERVALS.        *
!*                 =-3  NO SOLUTION DETERMINED BECAUSE SECBI FAILED.   *
!*                 =-4  ONE OR MORE VALUES FOR WHICH INTERPOLATION     *
!*                      HAS BEEN REQUESTED REQUIRES EXTRAPOLATION      *
!*                      WHICH IS NOT PERFORMED.                        *
!*                 =-5  NO SOLUTION DETERMINED AS A MATRIX USED IN     *
!*                      CALCULATING THE SPLINE WAS NOT SYMMETRIC       *
!*                      POSITIVE DEFINITE AS REQUIRED.  IT MAY BE      *
!*                      POSSIBLE TO ELIMINATE THE PROBLEM BY SPACING   *
!*                      KNOTS FURTHER APART.                           *
!*                 =-6  NN IS -1 OR ZERO.                              *
!*                 =-7  XK(1) OR XK(ABS(NK)) DOES NOT COINCIDE WITH    *
!*                      X(1) OR X(N), RESPECTIVELY.                    *
!*                 =-8  NK IS LESS THAN 3.                             *
!*                 =-9  PK(I) <= -1.0                                  *
!*                 = J  THE J-TH ELEMENT OF THE X ARRAY IS NOT IN      *
!*                      STRICTLY INCREASING ORDER.                     *
!*                                                                     *
!*                                                                     *
!*                 UPON RETURN FROM RSUTDS, THIS PARAMETER SHOULD BE   *
!*                 TESTED IN THE CALLING PROGRAM.                      *
!*                                                                     *
!*                                                                     *
!*                                                                     *
!*    REQUIRED ROUTINES                QXZ063,QXZ275,QXZ276,QXZ277,    *
!*                                     QXZ278,QXZ279,QXZ280,SECBI(C2.5)*
!*                                     SYMPDS(F4.3)                    *
!*                                                                     *
!*    SOURCE                           PROGRAM RASPFIT BY JAMES R.     *
!*                                     SCHIESS AND PATRICIA A. KERR    *
!*                                     MODIFIED BY COMPUTER SCIENCES   *
!*                                     CORPORATION                     *
!*                                                                     *
!*    LANGUAGE                         FORTRAN                         *
!*                                                                     *
!*    DATE RELEASED                    FEBRUARY, 1986                  *
!*                                                                     *
!*    LATEST REVISION                  FEBRUARY, 1986                  *
!***********************************************************************
!
!
      DIMENSION X(*),Y(*),XK(*),IPK(*),PK(*),DEV(*),                    &
     & R(*),WK(*),C(NKMAX,*),DF(*),XR(*)
      IERR=0
      NN1=IABS(NN)
      NK1=NK-1
      IZ=N-1

      DO 10 I=1,IZ
      IF (X(I+1)-X(I).GT.0) CYCLE   !!! GO TO 10
      IERR=I+1
      RETURN
   10 END DO

      IF(NK.GE.3)GO TO 12
      IERR=-8
      GO TO 120
   12 IF((NN.NE.-1).AND.(NN.NE.0))GO TO 15
      IERR=-6
      GO TO 120

   15 IF(NN.GT.0)GO TO 30
      XDIF=XR(NN1)-XR(1)
      XDIF=XDIF/(NN1-1)
      NN2=NN1-2

      DO 20 I=1,NN2
      XR(I+1)=XR(I)+XDIF
   20 END DO
   30 CONTINUE

      DO 40 I=1,NN1
      IF((XR(I).GE.XK(1)).AND.(XR(I).LE.XK(NK))) CYCLE   !!! GO TO 40
      IERR=-4
      GO TO 120
   40 END DO

      IF(XK(1).EQ.X(1).AND.XK(NK).EQ.X(N))GO TO 45
      IERR=-7
      GO TO 120

   45 DO 47 I=1,NK1
      IF (PK(I).GT.-1.0) CYCLE !!! GO TO 47
      IERR=-9
      GO TO 120
   47 END DO

   60 DO 77 I=1,N
      DF(I)=1.0/(DF(I)+1.0)
   77 END DO

      L=NK
      LM1=L-1
      LM2=L-2
      L2=L*2
      LL=L2+1
      INDA=1
      INDH=2*L+1
      INDE=3*L
      INDD=4*L**2+5*L
      INDGG=8*L**2+5*L
      INDG=9*L**2+L+4
      INDB=10*L**2-3*L+8
      INDXLAM=12*L**2-7*L+8
      INDT=12*L**2-6*L+6
      IF(IW.EQ.1)GO TO 500
      CALL QXZ275(X,Y,DF,XK,PK,IPK,DEV,LM1,N,L,ITMAX,XR,R,NN,NN1,       &
     & C,WK(INDA),L2,LM2,LL,WK(INDE),WK(INDD),WK(INDGG),                 &
     & WK(INDG),WK(INDB),WK(INDXLAM),WK(INDT),IW,IERR,NKMAX)
      GO TO 130
  500 CALL QXZ275(X,Y,DF,XK,PK,IPK,DEV,LM1,N,L,ITMAX,XR,R,NN,NN1,       &
     & C,WK(INDA),L2,LM2,LL,WK(INDE),WK(INDD),WK(INDGG),                 &
     & WK(INDG),WK(INDB),WK(INDXLAM),WK(INDT),IW,IERR,NKMAX)
  130 DO 88 I=1,N
      DF(I)=1.0/DF(I)-1.0
   88 END DO

      NK1=NK-1
      DO 133 I=1,NK1
      DXK=(XK(I+1)-XK(I))**2
      WK(INDH+I-1)=DXK/(2.*PK(I)**2+3.*PK(I)+3.)
  133 END DO
  120 RETURN
      END Subroutine RSUTDS

      SUBROUTINE RSDSDR(N,NK,X,Y,DF,XK,IPK,PK,DEV,NN,IW,                &
     &NKMAX,ITMAX,XR,R,DER1,DER2,C,WK,IERR)
!***********************************************************************
!*                                                                     *
!*    PURPOSE:                                                         *
!*                 SUBROUTINE RSDSDR FITS A SMOOTH RATIONAL SPLINE TO A*
!*                 UNIVARIATE FUNCTION AND PROVIDES INTERPOLATORY      *
!*                 CAPABILITIES FOR THE DEPENDENT VARIABLE AND THE     *
!*                 FIRST AND SECOND DERIVATIVES.  DATA MAY BE UNEQUALLY*
!*                 SPACED.                                             *
!*                                                                     *
!     D4.7                                                             *
!*                                                                     *
!*        USE:                                                         *
!*                 CALL RSDSDR(N,NK,X,Y,DF,XK,IPK,PK,DEV,NN,IW,NKMAX,  *
!*                             ITMAX,XR,R,DER1,DER2,C,WK,IERR)         *
!*                                                                     *
!*         N       AN INPUT INTEGER SPECIFYING THE NUMBER OF DATA      *
!*                 POINTS FOR THE INDEPENDENT VARIABLE.                *
!*                                                                     *
!*         NK      AN INPUT INTEGER SPECIFYING THE NUMBER OF KNOTS.    *
!*                 NK>=3.                                              *
!*                                                                     *
!*         X       AN INPUT ONE-DIMENSIONAL REAL ARRAY DIMENSIONED AT  *
!*                 AT LEAST N IN THE CALLING PROGRAM.  X CONTAINS THE  *
!*                 X-COORDINATES OF THE DATA POINTS.                   *
!*                                                                     *
!*         Y       AN INPUT ONE-DIMENSIONAL REAL ARRAY WHICH           *
!*                 CONTAINS THE FUNCTION VALUES OF THE DATA POINTS.    *
!*                 Y MUST BE DIMENSIONED AT LEAST N.                   *
!*                                                                     *
!*         DF      AN INPUT ONE-DIMENSIONAL REAL ARRAY WHICH CONTAINS  *
!*                 A SET OF POSITIVE WEIGHTS INDICATING THE STANDARD   *
!*                 DEVIATION OF THE ERROR OF Y(I).  IF THE STANDARD    *
!*                 DEVIATION IS NOT KNOWN, A SUGGESTED VALUE IS        *
!*                 DF(I)=1.0E-4, I=1,...,N.  DF MUST BE DIMENSIONED    *
!*                 AT LEAST N.                                         *
!*                                                                     *
!*         XK      AN INPUT ONE-DIMENSIONAL REAL ARRAY WHICH           *
!*                 CONTAINS THE X-COORDINATES OF THE SELECTED KNOTS.   *
!*                 XK MUST BE DIMENSIONED AT LEAST NK.  THE FIRST      *
!*                 AND LAST KNOTS MUST COINCIDE WITH X(1) AND X(N),    *
!*                 RESPECTIVELY, WHILE THE INTERVENING KNOTS MUST BE   *
!*                 CHOSEN TO LIE BETWEEN X(1) AND X(N),                *
!*                 IE, X(1)=XK(1)  XK(2)...XK(NK-1) XK(NK)=X(N).       *
!*                 KNOTS MUST BE CHOSEN SO THAT THERE ARE AT LEAST     *
!*                 THREE DATA POINTS IN EACH INTERVAL DEFINED BY       *
!*                 CONSECUTIVE PAIRS OF KNOTS.                         *
!*                                                                     *
!*         IPK     AN INPUT ONE-DIMENSIONAL INTEGER ARRAY WHICH        *
!*                 CONTAINS THE TENSION FACTOR ADJUSTMENT FLAGS.       *
!*                 IPK MUST BE DIMENSIONED AT LEAST NK-1.              *
!*                                                                     *
!*                 IPK(K)=0  FACTOR K MAY BE ADJUSTED.                 *
!*                                                                     *
!*                 IPK(K)=1  FACTOR K IS HELD CONSTANT.                *
!*                                                                     *
!*         PK      AN INPUT/OUTPUT ONE-DIMENSIONAL REAL ARRAY WHICH    *
!*                 CONTAINS THE TENSION FACTORS APPLIED ON THE         *
!*                 INTERVALS BETWEEN KNOTS.  PK MUST BE DIMENSIONED AT *
!*                 LEAST NK-1.  PK(K) MUST BE GREATER THAN -1. ON ENTRY*
!*                 IF IPK(K)=1, PK(K) WILL CONTAIN THE FINAL TENSION   *
!*                 FACTOR TO BE USED IN INTERVAL K.  IF IPK(K)=0, PK(K)*
!*                 MUST CONTAIN AN INITIAL VALUE OF TENSION FACTOR FOR *
!*                 INTERVAL K.  A VALUE OF 0.0 IS SUGGESTED.  ON       *
!*                 RETURN, PK(K) WILL CONTAIN THE FINAL TENSION FACTOR *
!*                 FOR INTERVAL K.                                     *
!*                                                                     *
!*         DEV     AN INPUT ONE-DIMENSIONAL REAL ARRAY WHICH           *
!*                 CONTAINS THE ALLOWED PERCENTAGE OF DEVIATION OF     *
!*                 THE APPROXIMATED SOLUTION FROM THE LINEAR SOLUTION. *
!*                 DEV MUST BE DIMENSIONED AT LEAST NK-1.              *
!*                 IN THE INTERVAL BETWEEN KNOT I AND KNOT I+1, FOR    *
!*                 I=1,...,NK-1. DEV MUST BE DIMENSIONED AT LEAST NK-1.*
!*                                                                     *
!*         NN      AN INPUT INTEGER WHOSE ABSOLUTE VALUE SPECIFIES     *
!*                 THE NUMBER OF POINTS AT WHICH INTERPOLATION         *
!*                 IS DESIRED.  IF NN IS POSITIVE, THE POINTS AT WHICH *
!*                 INTERPOLATION IS DESIRED MAY BE UNEQUALLY SPACED.   *
!*                 IF NN IS NEGATIVE, THE POINTS WILL BE EQUALLY SPACED*
!*                 BETWEEN XR(1) AND XR(ABS(NN)).  NN CANNOT BE -1     *
!*                 OR ZERO.                                            *
!*                                                                     *
!*                                                                     *
!*         IW      AN INPUT/OUTPUT INTEGER.                            *
!*                                                                     *
!*                 INPUT   ON THE FIRST CALL TO RSDSDR , IW MUST BE    *
!*                         SET TO 0.  ON SUBSEQUENT CALLS, IW SHOULD BE*
!*                         SET TO 1 TO PROVIDE ADDITIONAL INTERPOLATION*
!*                         ON THE SAME FITTING CALCULATIONS.           *
!*                                                                     *
!*                 OUTPUT  IW IS AUTOMATICALLY SET TO 1.               *
!*                                                                     *
!*                                                                     *
!*                                                                     *
!*         NKMAX   AN INPUT INTEGER WHICH AGREES EXACTLY WITH THE      *
!*                 FIRST DIMENSION OF ARRAY C.  NKMAX MUST BE GREATER  *
!*                 THAN OR EQUAL TO 2*NK.
!*                                                                     *
!*         ITMAX   AN INPUT INTEGER SPECIFYING THE MAXIMUM NUMBER      *
!*                 OF ITERATIONS PERMITTED IN ACHIEVING THE ALLOWED    *
!*                 DEVIATIONS.   ITMAX <=35.
!*                                                                     *
!*         XR      AN INPUT/OUTPUT  ONE-DIMENSIONAL REAL ARRAY         *
!*                 DIMENSIONED AT LEAST ABS(NN), WHICH CONTAINS        *
!*                 POINTS TO BE INTERPOLATED.                          *
!*
!*                 INPUT    IF NN IS POSITIVE, THE FIRST NN ELEMENTS   *
!*                          OF XR MUST ALL BE DEFINED.  IF NN IS       *
!*                          NEGATIVE, ONLY ELEMENTS XR(1) AND
!*                          XR(ABS(NN)) NEED BE DEFINED.
!*
!*                 OUTPUT   IF NN IS POSITIVE, XR IS UNCHANGED.  IF NN *
!*                          IS NEGATIVE, THE FIRST ABS(NN) ELEMENTS OF *
!*                          XR CONTAIN EQUALLY SPACED INTERPOLATION    *
!*                          POINTS.                                    *
!*                                                                     *
!*                                                                     *
!*         R       AN OUTPUT ONE-DIMENSIONAL REAL ARRAY WHICH CONTAINS *
!*                 THE VALUES OF THE INTERPOLATED DEPENDENT VARIABLE   *
!*                 CORRESPONDING TO XR.  R MUST BE DIMENSIONED AT LEAST*
!*                 ABS( NN).                                           *
!*                                                                     *
!*         DER1     AN OUTPUT ONE-DIMENSIONAL REAL ARRAY WHICH CONTAINS*
!*                 THE VALUES OF THE FIRST DERIVATIVE OF THE INTER-    *
!*                 POLATED DEPENDENT VARIABLE CORRESPONDING TO XR.     *
!*                 DER1 MUST BE DIMENSIONED AT LEAST ABS(NN).          *
!*                                                                     *
!*         DER2     AN OUTPUT ONE-DIMENSIONAL REAL ARRAY WHICH CONTAINS*
!*                 THE VALUES OF THE SECOND DERIVATIVE OF THE INTER-   *
!*                 POLATED DEPENDENT VARIABLE CORRESPONDING TO XR.     *
!*                 DER2 MUST BE DIMENSIONED AT LEAST ABS(NN).          *
!*                                                                     *
!*         C       AN OUTPUT TWO-DIMENSIONAL REAL ARRAY WHICH          *
!*                 CONTAINS THE CORRELATION MATRIX.  C MUST  HAVE      *
!*                 A FIRST DIMENSION OF NKMAX AND A SECOND DIMENSION   *
!*                 OF AT LEAST 2*NK.                                   *
!*                                                                     *
!*         WK      AN INPUT/OUTPUT ONE-DIMENSIONAL REAL ARRAY WHICH    *
!*                 MUST BE DIMENSIONED AT LEAST 16*NK**2-6*NK+5.       *
!*                                                                     *
!*                 INPUT   IF IW = 0, WK NEED NOT CONTAIN VALUES.      *
!*                                                                     *
!*                         IF IW = 1, WK SHOULD CONTAIN THE VALUES     *
!*                         WHICH WERE STORED IN WK BY A PREVIOUS CALL  *
!*                         TO RSDSDR.                                  *
!*                                                                     *
!' *               OUTPUT  WK CONTAINS INFORMATION NECESSARY FOR       *
!*                         SUCCEEDING CALLS TO RSDSDR USING THE SAME   *
!*                         INPUT DATA WITHOUT REINITIALIZATION.        *
!*                                                                     *
!*        IERR     OUTPUT ERROR PARAMETER:                             *
!*                 = 0  NORMAL RETURN.  NO ERROR DETECTED.             *
!*                 =-1  AN INTERVAL CONTAINS LESS THAN THREE DATA      *
!*                      POINTS.                                        *
!*                 =-2  DID NOT CONVERGE IN ALLOWED NUMBER OF ITERAT-  *
!*                      IONS, INDICATING THAT TENSION MAY NOT HAVE     *
!*                      BEEN ADJUSTED IN ONE OR MORE INTERVALS.        *
!*                 =-3  NO SOLUTION DETERMINED BECAUSE SECBI FAILED.   *
!*                 =-4  ONE OR MORE VALUES FOR WHICH INTERPOLATION     *
!*                      HAS BEEN REQUESTED REQUIRES EXTRAPOLATION      *
!*                      WHICH IS NOT PERFORMED.                        *
!*                 =-5  NO SOLUTION DETERMINED AS A MATRIX USED IN     *
!*                      CALCULATING THE SPLINE WAS NOT SYMMETRIC       *
!*                      POSITIVE DEFINITE AS REQUIRED.  IT MAY BE      *
!*                      POSSIBLE TO ELIMINATE THE PROBLEM BY SPACING   *
!*                      KNOTS FURTHER APART.                           *
!*                 =-6  NN IS -1 OR ZERO.                              *
!*                 =-7  XK(1) OR XK(ABS(NK)) DOES NOT COINCIDE WITH    *
!*                      X(1) OR X(N), RESPECTIVELY.                    *
!*                 =-8  NK IS LESS THAN 3.                             *
!*                 = J  THE J-TH ELEMENT OF THE X ARRAY IS NOT IN      *
!*                      STRICTLY INCREASING ORDER.                     *
!*                                                                     *
!*                                                                     *
!*                 UPON RETURN FROM RSDSDR, THIS PARAMETER SHOULD BE   *
!*                 TESTED IN THE CALLING PROGRAM.                      *
!*                                                                     *
!*                                                                     *
!*                                                                     *
!*    REQUIRED ROUTINES                QXZ063,QXZ276,QXZ277,QXZ278,    *
!*                                     QXZ279,QXZ280,QXZ281,QXZ282,    *
!*                                     SECBI(C2.5),SYMPDS(F4.3)        *
!*                                                                     *
!*    SOURCE                           PROGRAM RASPFIT BY JAMES R.     *
!*                                     SCHIESS AND PATRICIA A. KERR    *
!*                                     MODIFIED BY COMPUTER SCIENCES   *
!*                                     CORPORATION                     *
!*                                                                     *
!*    LANGUAGE                         FORTRAN                         *
!*                                                                     *
!*    DATE RELEASED                    FEBRUARY, 1986                  *
!*                                                                     *
!*    LATEST REVISION                  FEBRUARY, 1986                  *
!***********************************************************************
!
!
      DIMENSION X(1),Y(1),XK(1),IPK(1),PK(1),DEV(1),                    &
     &R(1),WK(1),C(NKMAX,1),DF(1),XR(1),DER1(1),DER2(1)
      IERR=0
      NN1=IABS(NN)
      IZ=N-1

      DO 10 I=1,IZ
      IF (X(I+1)-X(I).GT.0) CYCLE   !!!GO TO 10
      IERR=I+1
      RETURN
   10 END DO

      IF(NK.GE.3)GO TO 12
      IERR=-8
      GO TO 120
   12 IF((NN.NE.-1).AND.(NN.NE.0))GO TO 15
      IERR=-6
      GO TO 120
   15 IF(NN.GT.0)GO TO 30
      XDIF=XR(NN1)-XR(1)
      XDIF=XDIF/(NN1-1)
      NN2=NN1-2

      DO 20 I=1,NN2
      XR(I+1)=XR(I)+XDIF
   20 END DO

   30 CONTINUE
      DO 40 I=1,NN1
      IF ((XR(I).GE.XK(1)).AND.(XR(I).LE.XK(NK))) CYCLE   !!! GO TO 40
      IERR=-4
      GO TO 120
   40 END DO

      IF(XK(1).EQ.X(1).AND.XK(NK).EQ.X(N))GO TO 60
      IERR=-7
      GO TO 120

   60 DO 77 I=1,N
      DF(I)=1.0/(DF(I)+1.0)
   77 END DO

      L=NK
      LM1=L-1
      LM2=L-2
      L2=L*2
      LL=L2+1
      INDA=1
      INDH=2*L+1
      INDE=3*L
      INDD=4*L**2+5*L
      INDGG=8*L**2+5*L
      INDG=9*L**2+L+4
      INDB=10*L**2-3*L+8
      INDXLAM=12*L**2-7*L+8
      INDT=12*L**2-6*L+6
      IF(IW.EQ.1)GO TO 500
      CALL QXZ281(X,Y,DF,XK,PK,IPK,DEV,LM1,N,L,ITMAX,XR,R,NN,NN1,DER1,  &
     &DER2,C,WK(INDA),L2,LM2,LL,WK(INDE),WK(INDD),WK(INDGG),            &
     &WK(INDG),WK(INDB),WK(INDXLAM),WK(INDT),IW,IERR,NKMAX)
      GO TO 130
  500 CALL QXZ281(X,Y,DF,XK,PK,IPK,DEV,LM1,N,L,ITMAX,XR,R,NN,NN1,DER1,  &
     &DER2,C,WK(INDA),L2,LM2,LL,WK(INDE),WK(INDD),WK(INDGG),            &
     &WK(INDG),WK(INDB),WK(INDXLAM),WK(INDT),IW,IERR,NKMAX)
  130 DO 88 I=1,N
      DF(I)=1.0/DF(I)-1.0
   88 END DO
      NK1=NK-1
      DO 133 I=1,NK1
      DXK=(XK(I+1)-XK(I))**2
      WK(INDH+I-1)=DXK/(2.*PK(I)**2+3.*PK(I)+3.)
  133 END DO
  120 RETURN
      END Subroutine RSDSDR

      SUBROUTINE SECBI (X1,X2,DX,FX,EPS,ROOT,IERR)
!     C2.5
!***********************************************************************
!*
!*            PURPOSE
!*               TO DETERMINE A ROOT OF THE REAL VALUED FUNCTION  F(X)=0
!*               GIVEN A SPECIFIED INTERVAL BY EMPLOYING A FRONT END
!*               SEEKER AND A COMBINATION BISECTION/LINEAR INTERPOLATION
!*               INVERSE QUADRATIC INTERPOLATION ITERATION TECHNIQUE
!*
!*            USAGE
!*               CALL SECBI (X1,X2,DX,FX,EPS,ROOT,IERR)
!*
!*            PARAMETERS
!*               X1      - INPUT  VARIABLE  CONTAINING THE INITIAL VALUE
!*                         OF THE INDEPENDENT VARIABLE TO START THE
!*                         ITERATION
!*               X2      - INPUT VARIABLE CONTAINING THE FINAL VALUE OF
!*                         THE INDEPENDENT VARIABLE TO END THE ITERATION
!*                         FX MUST BE NAMED IN AN EXTERAL STATEMENT IN
!*                         THE CALLING PROGRAM
!*               DX      - AN INPUT PARAMETER CONTAINING THE INCREMENT
!*                         STEP SIZE
!*               FX      - A VALUE OF THE USER SUPPLIED EXTERAL FUNCTION
!*                         FX MUST BE IN AN EXTERNAL STATEMENT IN
!*                         THE CALLING PROGRAM
!*               EPS     - A ONE DIMENSIONAL ERROR ARRAY DIMESIONED 3
!*                         IN THE CALLING PROGRAM
!*               EPS(1)  - MAXIMUM RELATIVE ERROR
!*               EPS(2)  - MAXIMUM ABSOLUTE ERROR
!*               EPS(3)  - THRESHOLD ERROR
!*                         SEE SUBROUTINE WRITE-UP FOR DETAILED DESCRIP-
!*                         TION OF EPS ARRAY
!*               ROOT    - AN OUTPUT PARAMETER CONTAINING THE ROOT OF
!*                         THE FUNCTION FX
!*               IERR    - AN INTEGER RETURNED BY SECBI AS AN ERROR CODE
!*                         IT SHOULD BE TESTED BY THE USER UPON RETURN
!*                         TO THE CALLING PROGRAM
!*                            - 0 NORMAL RETURN
!*                            - 3 MAXIMUM OF 50 ITERATIONS IN SUBROUTINE
!*                                QXZ019 IS EXCEEDED
!*                            - 8 IMPROPER INITIAL CONDITIONS WERE
!*                                SPECIFIED (SEE SUBROUTINE WRITE-UP)
!*                            - 9 A ROOT CANNOT BE FOUND IN THE INTER-
!*                                VAL BOUNDED BY X1 AND X2 WITH THE
!*                                GIVEN DX. (SEE SUBROUTINE WRITE-UP)
!*
!*            REQUIRED ROUTINES
!                 SUBROUTINE QXZ019(F,EPS,A,B,ITMAX,IERR)
!*               SUBROUTINE QXZ255 (EPSILON,VAL,IFLAG)
!*               FUNCTION FX(X)
!*
!*            AUTHOR / IMPLEMENTOR
!*               COMPUTER SCIECES CORPORATION / ACD PROGRAMMER SUPPORT
!*                                              GROUP EXT. 3548
!*
!*            LANGUAGE                 FORTRAN
!*
!*            DATE RELEASED            MARCH 1973
!*
!*            LATEST REVISION          AUGUST 15, 1980
!*
!***********************************************************************
!
      DIMENSION EPS(3),X(2),F(2)
      EXTERNAL FX
      LOGICAL IFLAG
!
!             DETERMINE IF INPUT MEETS INITIAL CONDITIONS AND RESET THE
!             THE INCREMENT STEP SIZE WHEN NECESSARY
      IFLAG = .FALSE.
      IF(DX.EQ.0)GO TO 100
      DX=ABS(DX)
      IF (X2-X1) 10,100,20
   10 DX = -DX
   20 MAX=ABS((X2-X1)/DX)
      IF(MAX.GT.100)GO TO 100
!
!             SET UP FIRST LOWER BOUND ON SEQUENCE OF SUBINTERVALS AND
!             COMPUTE FUNCTION VALUE
      X(1)=X1
      F(1)=FX(X(1))
      IF (F(1)) 40,30,40
!
   30 ROOT = X(1)
      GO TO 80
!             SECTION TO SEARCH FOR SUBINTERVAL (X1,X2) CONTAINING A
!
!             POSSIBLE ROOT
   40 X(2) = X(1)+DX
      IF (DX*X(2) .LE. DX*X2) GO TO 45
      IF (IFLAG) GO TO 90
      X(2) = X2
      IFLAG = .TRUE.
   45 F(2) = FX(X(2))
      IF ( F(1)*F(2) ) 70,50,60
!
   50 ROOT = X(2)
      GO TO 80
   60 X(1) = X(2)
      F(1)=F(2)
      GO TO 40
!
!            CALL SUBROUTINE QXZ019 TO APPLY BISECTION/SECANT METHOD TO
!             DETERMINE ROOT (ACCURACY DETERMINED BY SUBROUTINE QXZ255)
   70 ITMAX = 50
      CALL QXZ019(FX,EPS,X(1),X(2),ITMAX,IERR)
      ROOT = X(2)
      RETURN
!
   80 IERR = 0
      RETURN
   90 IERR = 9
      RETURN
  100 IERR = 8
      RETURN
      END Subroutine SECBI

      SUBROUTINE SYMPDS (MAXN,N,A,NRHS,B,IOPT,IFAC,DETERM,ISCALE,P,IERR)
!     F4.3
!***********************************************************************
!*                                                                     *
!*            PURPOSE                                                  *
!*               TO SOLVE THE MATRIX EQUATION AX=B WHERE A IS SYMMETRIC*
!*               POSITIVE DEFINITE MATRIX AND B IS A MATRIX OF CONSTANT*
!*               VECTORS. THE DETERMINANT OF MATRIX A MAY ALSO BE      *
!*               EVALUATED.                                            *
!*                                                                     *
!*            USAGE                                                    *
!*               CALL SYMPDS (MAXN,N,A,NRHS,B,IOPT,IFAC,DETERM,ISCALE,P*
!*                            ,IERR)                                   *
!*                                                                     *
!*            PARAMETERS                                               *
!*               MAXN    - AN INPUT INTEGER SPECIFYING THE MAXIMUN     *
!*                         FIRST DIMENSION OF THE COEFFICIENT MATRIX A *
!*                         AS GIVEN IN THE DIMENSION STATEMENT OF      *
!*                         CALLING PROGRAM                             *
!*               N       - AN INPUT INTEGER SPECIFYING THE ORDER OF    *
!*                         MATRIX A (THE ACTUAL FIRST DIMENSION OF A)  *
!*                         N@MAXN                                      *
!*               A       - A TWO DIMENSIONAL INPUT-OUTPUT ARRAY THAT   *
!*                         MUST HAVE FIRST DIMENSION MAXN AND SECOND   *
!*                         DIMENSION AT LEAST N IN THE CALLING PROGRAM.*
!*                              AS AN INPUT ARRAY, A CONTAINS THE ELE- *
!*                         MENTS OF THE COEFFICIENT MATRIX             *
!*                              UPON RETURN FROM THE SUBROUTINE A      *
!*                         CONTAINS THE UPPER TRIANGULAR ELEMENTS OF   *
!*                         THE COEFFICIENT MATRIX AND THE UNIT LOWER   *
!*                         TRIANGULAR MATRIX L.                        *
!*               NRHS    - AN INPUT INTEGER CONTAINING THE NUMBER OF   *
!*                         COLUMN VECTORS OF MATRIX B                  *
!*               B       - A TWO DIMENSIONAL INPUT-OUTPUT ARRAY THAT   *
!*                         MUST HAVE FIRST DIMENSION MAXN AND SECOND   *
!*                         DIMENSION AT LEAST NRHS IN CALLING PROGRAM  *
!*                              AS INPUT B CONTAINS THE ELEMENTS OF THE*
!*                         CONSTANT VECTORS                            *
!*                              UPON RETURN FROM THE SUBROUTINE, EACH  *
!*                         SOLUTION VECTOR IS STORED OVER THE CORRES-  *
!*                         PONDING CONSTANT COLUMN VECTOR OF MATRIX B  *
!*               IOPT    - AN INPUT INTEGER CONTAINING THE DETERMINANT *
!*                         EVALUATION CODE                             *
!*                         = 0 DETERMINANT IS NOT EVALUATED            *
!*                         = 1 DETERMINANT IS EVALUATED                *
!*               IFAC    - INPUT INTEGER PARAMETER SPECIFYING WHETHER  *
!*                         OR NOT CHOLESKY DECOMPOSITION OF THE COEF-  *
!*                         FICIENT MATRIX IS TO BE COMPUTED            *
!*                         = 0 CHOLESKY DECOMPOSITION FOR THE MATRIX A *
!*                             IS TO BE COMPUTED. IFAC IS RESET TO 1 BY*
!*                             THE ROUTINE.                            *
!*                         = 1 THE CHOLESKY DECOMPOSED FORM OF MATRIX  *
!*                             A IS INPUT SO NO DECOMPOSITION IS COM-  *
!*                             PUTED.                                  *
!*               DETERM  - DETERMINANT EVALUATION PARAMETERS RETURNED  *
!*               ISCALE    BY THE ROUTINE TO EVALUTE THE DETERMINANT   *
!*                         VALUE OF DETERMINANT                        *
!*                              = DETERM*(10.**(100*ISCALE))           *
!*                         (CALCULATION SHOULD BE DONE OUTSIDE CALLING *
!*                         PROGRAM TO PREVENT MACHINE OVER(UNDER) FLOW)*
!*               P       - A ONE DIMENSIONAL OUTPUT ARRAY CONTAINING   *
!*                         THE RECIPROCALS OF THE ELEMENTS IN THE      *
!*                         DIAGONAL MATRIX FORMED DURING THE DECOMPO-  *
!*                         SITION OF THE COEFFICIENT MATRIX A.  P MUST *
!*                         BE DIMENSIONED AT LEAST N IN THE CALLING    *
!*                         PROGRAM                                     *
!*               IERR    - AN INTEGER RETURNED BY THE SUBROUTINE AS AN *
!*                         ERROR CODE                                  *
!*                         = 0 NORMAL RETURN (SOLUTION VECTORS ARE     *
!*                             COMPUTED)                               *
!*                         = 1 COEFFICIENT MATRIX IS NOT SYMMETRIC     *
!*                             POSITIVE DEFINITE                       *
!*                                                                     *
!*            REQUIRED ROUTINES                                        *
!*               NONE                                                  *
!*                                                                     *
!*            AUTHOR / IMPLEMENTOR                                     *
!*               COMPUTER SCIENCES CORPORATION / ACD PROGRAMMER SUPPORT*
!*                                               GROUP EXT. 3548       *
!*                                                                     *
!*            LANGUAGE                  FORTRAN                        *
!*                                                                     *
!*            DATE RELEASED             AUGUST  1972                   *
!*                                                                     *
!*            LATEST REVISION           MAY 1985
!*                                                                     *
!***********************************************************************
!
!
      DIMENSION A(MAXN,N),B(MAXN,NRHS),P(N)
!
!!!      DATA R1,R2/1.0E+100,1.0E-100/
      REAL,PARAMETER:: R1 = HUGE(R1)
      REAL,PARAMETER:: R2 = TINY(R1)
!
!             TEST FOR A SCALAR MATRIX (IF COEFFICIENT MATRIX IS A
!             SCALAR MATRIX, SOLVE AND COMPUTE DETERMINANT IF DESIRED)
      NM1 = N-1
      IF (NM1 .GT. 0) GO TO 20
!
      IF (A(1,1) .LE. 0.0) GO TO 333
      ISCALE = 0
      DETERM = A(1,1)
      P(1) = 1.0/A(1,1)
      DO 10 J=1,NRHS
      B(1,J) = B(1,J)/DETERM
   10 END DO
      IERR = 0
      RETURN
!
!             TEST TO DETERMINE IF CHOLESKY DECOMPOSITION OF COEFFICIENT
!             MATRIX IS DESIRED
   20 IF (IFAC .EQ. 1) GO TO 1050
      IFAC = 1
!
!
!              "LOOP" TO PERFORM CHOLESKY DECOMPOSTION ON THE COEF-
!             FICIENT MATRIX A (I.E. MATRIX A WILL BE DECOMPOSED INTO
!             THE PRODUCT OF A UNIT LOWER TRIANGULAR MATRIX (L), A
!             DIAGONAL MATRIX (D), AND THE TRANSPOSE OF L (LTRANSPOSE).)
!
   30 DO 150 I=1,N
      IM1 = I-1
!
      DO 150 J=1,I
      X = A(J,I)
!
!             DETERMINE IF ELEMENTS ARE ABOVE OR BELOW THE DIAGONAL
      IF (I .GT. J) GO TO 110
!
!             USING THE DIAGONAL ELEMENTS OF MATRIX A, THIS SECTION
!             COMPUTES DIAGONAL MATRIX AND DETERMINES IF MATRIX A IS
!             SYMMETRIC POSITIVE DEFINITE
      IF (IM1 .EQ. 0) GO TO 50
      DO 40 K=1,IM1
      Y = A(I,K)
      A(I,K) = Y*P(K)
      X = X - Y*A(I,K)
   40 END DO
!
!             TEST IF COEFFICIENT MATRIX IS POSITIVE DEFINITE
   50 IF (X .LE. 0) GO TO 333
!
!             COMPUTE INVERSE OF DIAGONAL MATRIX D**-1 = 1/P
      P(I) = 1.0 / X
!
!             TEST TO SEE IF DETERMINANT IS TO BE EVALUATED
      GO TO 150
!
!
!             USING THE LOWER TRIANGULAR ELEMENTS OF MATRIX A, THIS
!             SECTION COMPUTES THE UNIT LOWER TRIANGULAR MATRIX
  110 JM1 = J-1
      IF (JM1 .EQ. 0) GO TO 140
      DO 120 K=1,JM1
      X = X - A(I,K)*A(J,K)
  120 END DO
!
  140 A(I,J) = X
!
  150 CONTINUE
 1050 CONTINUE
!
!             SECTION TO APPLY BACK SUBSTITUTION TO SOLVE L*Y = B FOR
!             UNIT LOWER TRIANGULAR MATRIX AND CONSTANT COLUMN VECTOR B
!
!      TEST TO SEE IF DETERMINANT IS TO BE EVALUATED
      IF(IOPT.EQ.0)GO TO 160
!     INITIALIZE DETERMINANT PARAMETERS
      DETERM = 1.0
      ISCALE = 0
!     COMPUTE DETERMINANT SCALING AS NECESSARY
      DO 1110 I=1,N
      PIVOTI=1.0/P(I)
 1060 IF(ABS(DETERM).LT.R1) GO TO 1070
      DETERM = DETERM*R2
      ISCALE = ISCALE+1
      GO TO 1060
 1070 IF(ABS(DETERM).GT.R2) GO TO 1080
      DETERM = DETERM*R1
      ISCALE = ISCALE-1
      GO TO 1070
 1080 IF(ABS(PIVOTI).LT.R1) GO TO 1090
      PIVOTI = PIVOTI*R2
      ISCALE = ISCALE+1
      GO TO 1080
 1090 IF(ABS(PIVOTI).GT.R2) GO TO 1100
      PIVOTI = PIVOTI*R1
      ISCALE = ISCALE-1
      GO TO 1090
 1100 DETERM = DETERM*PIVOTI
 1110 END DO
  160 DO 180 I=2,N
      IM1=I-1
!
      DO 180 J=1,NRHS
      X = B(I,J)
!
      DO 170 K=1,IM1
      X = X - A(I,K)*B(K,J)
  170 END DO
!
      B(I,J) = X
  180 CONTINUE
!
!             SECTION TO SOLVE (LTRANSPOSE)*X = (D**-1)*Y FOR TRANSPOSE
!             OF UNIT LOWER TRIANGULAR MATRIX AND INVERSE OF DIAGONAL
!             MATRIX
!
      Y = P(N)
      DO 190 J=1,NRHS
      B(N,J) = B(N,J)*Y
  190 END DO
!
  200 I = NM1+1
      Y = P(NM1)
!
      DO 220 J=1,NRHS
      X = B(NM1,J)*Y
!
      DO 210 K=I,N
      X = X - A(K,NM1)*B(K,J)
  210 END DO
!
      B(NM1,J) = X
  220 END DO
!
!
!             TEST TO DETERMINE IF SOLUTIONS HAVE BEEN DETERMINED FOR
!             ALL COLUMN VECTORS
      NM1 = NM1-1
      IF (NM1 .GT. 0) GO TO 200
!
      IERR = 0
      RETURN
!
  333 IERR = 1
      RETURN
      END Subroutine SYMPDS

! --- SUBPROGRAM QXZ019 --- FORMERLY KNOWN AS ROUTINE  ZBRENTL ---
       SUBROUTINE QXZ019(F,EPS,A,B,ITMAX,IERR)
!
!-ZBRENT--------S-------LIBRARY 3---------------------------------------
!
!   FUNCTION            - TO FIND A ZERO OF A FUNCTION WHICH CHANGES
!                           SIGN IN A GIVEN INTERVAL
!   USAGE                - CALL QXZ019(F,EPS,A,B,ITMAX,IERR)
!   PARAMETERS   F      - AN EXTERNAL FUNCTION SUBPROGRAM F(X)
!                           PROVIDED BY THE USER WHICH COMPUTES F FOR
!                           ANY X IN THE INTERVAL (A,B).
!                 EPS   - A ONE DIMENSIONAL ERROR ARRARY SUPPLIED BY THE
!                          USER. THIS ARRAY MUST BE DIMENSIONED 3 IN THE
!                          CALLING PROGRAM
!                 EPS(1) - MAXIMUM RELATIVE ERROR
!                 EPS(2) - MAXIMUM ABSOLUTE ERROR
!                 EPS(3) - THRESHOLD ERROR
!                          SEE SUBROUTINE SECBI WRITE-UP FOR A MORE
!                          DETAILED DESCRIPTION
!                A,B    - ON INPUT, THE USER MUST SUPPLY TWO POINTS, A
!                           AND B, SUCH THAT F(A) AND F(B) ARE OPPOSITE
!                           IN SIGN.
!                           ON OUTPUT, BOTH A AND B ARE ALTERED.  B
!                           WILL CONTAIN THE BEST APPROXIMATION TO THE
!                           ROOT OF F.
!                ITMAX  - ON INPUT, ITMAX SHOULD CONTAIN AN UPPER BOUND
!                           ON THE NUMBER OF ITERATIONS REQUIRED FOR
!                           CONVERGENCE. ON OUTPUT, ITMAX WILL CONTAIN
!                           THE ACTUAL ITERATION COUNT.
!               IERR   - ERROR PARAMETER
!                        TERMINAL ERROR
!                          3 - INDICATES THE ALGORITHM FAILED TO
!                           CONVERGE IN ITMAX ITERATIONS.
!                          8 - INDICATES F(A) AND F(B) HAVE
!                           THE SAME SIGN
!   PRECISION           - SINGLE
!   REQD. ROUTINES       - NONE
!   LANGUAGE            - FORTRAN
!   LATEST REVISION      - AUGUST 4,1975
!-----------------------------------------------------------------------
!
      DIMENSION EPS(3),V(2)
      DATA               ZERO,HALF,ONE,THREE,TEN/0.0,.5,1.0,3.0,10.0/
      IC = 0
      IERR = 0
      FA = F(A)
      FB = F(B)
      IF (FA*FB) 5,90,95
    5 C = A
      FC = FA
      D = B - A
      E = D
   10 IF (ABS(FC) .GE. ABS(FB)) GO TO 15
      A = B
      B = C
      C = A
      FA = FB
      FB = FC
      FC = FA
   15 TOL = EPS(1)*ABS(B)
      IF (FB) 16,111,16
  111 IERR = 0
      GO TO 80
   16 RM = (C-B)*HALF
!                       TEST FOR B AS A ROOT. A IS
!                      PREVIOUS BEST APPROXIMATION.
      V(1) = A
      V(2) = B
      CALL QXZ255(EPS,V,IERR)
      IF (IERR) 20,80,20
   20 IF (ABS(E) .LT. TOL .OR. ABS(FA) .LE. ABS(FB)) GO TO 60
      S = FB/FA
      IF (A-C) 30,25,30
!                                  LINEAR INTERPOLATION
   25 P = (RM+RM) * S
      Q = ONE-S
      GO TO 35
!                                  INVERSE QUADRATIC INTERPOLATION
   30 Q = FA/FC
      R = FB/FC
      RONE = R-ONE
      P = S * ((RM+RM) * Q * (Q-R) - (B-A) * RONE)
      Q = (Q-ONE) * RONE * (S-ONE)
   35 IF (P) 40,40,45
   40 P = -P
      GO TO 50
   45 Q = -Q
   50 S = E
      E = D
      IF ((P+P).GE.THREE*RM*Q-ABS(TOL*Q).AND.P.GE.ABS(HALF*S*Q))GO TO 55
      D = P/Q
      GO TO 65
   55 E = RM
      D = E
      GO TO 65
   60 E = RM
      D = E
   65 A = B
      FA = FB
   70 TEMP = D
!                                  INCREMENT ITERATION COUNTER
   75 IC = IC + 1
!                                  TERMINAL ERROR--MAXIMUM NUMBER OF
!                                                  ITERATIONS EXCEDED
      IF (IC .LE. ITMAX) GO TO 85
      IERR = 3
   80 ITMAX = IC
      GO TO 9005
   85 B = B + TEMP
      FB = F(B)
      FBC = FB*FC
      IF (FBC) 10,10,5
!                                  EITHER F(A) OR F(B) IS ZERO
   90 IF (FA .EQ. ZERO) B = A
      GO TO 80
!                                  TERMINAL ERROR - F(A) AND F(B) HAVE
!                                  THE SAME SIGN
   95 IERR = 8
      IC = 1
      GO TO 80
 9005 RETURN
      END Subroutine QXZ019

! --- SUBPROGRAM QXZ063 --- FORMERLY KNOWN AS ROUTINE  QINTRVL ---
      INTEGER FUNCTION QXZ063(T,X,N)
!
!***********************************************************************
!*                                                                     *
!*                                                                     *
!*  PURPOSE             - DETERMINE THE INDEX OF THE INTERVAL          *
!*                        (DETERMINED BY A GIVEN INCREASING SEQUENCE)  *
!*                        IN WHICH A GIVEN VALUE LIES.                 *
!*                                                                     *
!*  USE                 - I = QXZ063(T,X,N)                           *
!*                                                                     *
!*  PARAMETERS   T      - AN INPUT REAL NUMBER SPECIFYING THE          *
!*                        GIVEN VALUE.                                 *
!*               X      - AN INPUT ONE-DIMENSIONAL REAL ARRAY OF LENGTH*
!*                        N  SPECIFYING THE INCREASING SEQUENCE.  X    *
!*                        MUST BE STRICTLY INCREASING.                 *
!*               N      - AN INPUT INTEGER SPECIFYING THE LENGTH OF    *
!*                        X .  N  > 1.                                 *
!*                                                                     *
!*  OUTPUT       I      - IF  T  .LE.  X(2) ,  I  = 1.                 *
!*                        IF  T  .GE.  X(N-1) ,  I  =  N - 1 .         *
!*                        OTHERWISE,  X(I)  .LE.  T  .LE.  X(I+1)      *
!*                                                                     *
!*  PRECISION           - SINGLE.                                      *
!*                                                                     *
!*  REQUIRED ROUTINES   - NONE.                                        *
!*                                                                     *
!*  DATE RELEASED       - MARCH 1, 1979.                               *
!*                                                                     *
!*  SOURCE              - A. K. CLINE AND R. J. RENKA                  *
!*                        UNIVERSITY OF TEXAS AT AUSTIN                *
!*                                                                     *
!*  LATEST REVISION     - NONE.                                        *
!*                                                                     *
!*                                                                     *
!***********************************************************************
!
!
!             FORMAL PARAMETERS
!
      INTEGER N
!
      REAL T,X(N)
!
!             INTERNAL VARIABLES
!
      INTEGER IH,IL
!
      REAL TT
!
      TT = T
      IF (TT .LE. X(2)) GO TO 4
         IF (TT .GE. X(N-1)) GO TO 5
            IL = 2
            IH = N-1
!
!             INTERPOLATE LINEARLY FOR I
!
    1       I = IL+IFIX(FLOAT(IH-IL)*(TT-X(IL))/(X(IH)-X(IL)))
            IF (I .EQ. IH) I = I - 1
            IF (TT .LT. X(I)) GO TO 2
               IF (TT .LE. X(I+1)) GO TO 3
!
!             I  IS TOO SMALL - ADJUST AND TRY AGAIN
!
                  IL = I+1
                  GO TO 1
!
!             I  IS TOO LARGE - ADJUST AND TRY AGAIN
!
    2       IH = I
            GO TO 1
!
    3          QXZ063 = I
!             I  IS JUST RIGHT - RETURN
!
               RETURN
!
!             LEFT END
!
    4 QXZ063 = 1
      RETURN
!
!             RIGHT END
!
    5    QXZ063 = N-1
         RETURN
!
      END Function QXZ063


! --- SUBPROGRAM QXZ255 --- FORMERLY KNOWN AS ROUTINE  ERROR   ---
      SUBROUTINE QXZ255 (EPSILON,VAL,IFLAG)
!***********************************************************************
!*
!*            PURPOSE
!*              TO TEST THE CONVERGENGE OF A COMPUTED RESULT BASED ON A
!*              RELATIVE CONVERGENCE CRITERION OR AN ABSOLUTE
!*              CONVERGENCE CRITERION
!*
!*            USAGE
!*              CALL ERROR (EPSILON,VAL,IFLAG)
!*            PARAMETERS
!*               EPSILON - A ONE DIMENSION INPUT ARRAY DIMENSIONED 3 IN
!*                         THE CALLING PROGRAM
!*               EPSILON(1) - MAXIMUM RELATIVE ERROR
!*               EPSILON(2) - MAXIMUM ABSOLUTE ERROR
!*               EPSILON(3) - THRESHOLD ERROR
!*                         SEE ACCURACY SECTION OF WRITE-UP FOR SUB-
!*                         ROUTINES SECBI OR RMULL FOR A MORE DETAILED
!*                         DESCRIPTION
!*               VAL     - A ONE DIMENSIONAL INPUT ARRAY CONTAINING THE
!*                         COMPUTED VALUES FROM SUCCEEDING ITERATIONS
!*                         WHOSE CONVERGENGE IS TO BE DETERMINED
!*                            VAL(1) - COMPUTED VALUE FROM ITH ITERATION
!*                            VAL(2) - COMPUTED VALUE FROM I+1TH
!*                                     ITERATION
!*               IFLAG   - AN OUTPUT ERROR INTEGER
!*                         = 0 CONVERGENGE
!*                         = 1 NON CONVERGENCE
!*
!*            REQUIRED ROUTINES
!*               NONE
!*
!***********************************************************************
!
      DIMENSION EPSILON (*),VAL(*)
!
!
      IFLAG=1
!
!             SET THRESHOLD ERROR
      TE=EPSILON(3)
!
!             SET THRESHOLD ERROR TO DEFAULT VALUE IF 0.0
      IF(TE.EQ.0.)TE=EPSILON(1)
!
      TEMP=ABS(VAL(2)-VAL(1))
!
!             DETERMINE WHETHER RELATIVE ERROR OR ABSOLUTE ERROR WILL
!             BE USED TO DETERMINE CONVERGENCE
      IF(ABS(VAL(2)).LT.TE) GOTO 2
      GOTO 3
!
!             ABSOLUTE CONVERGENCE CRITERION
    2 IF(TEMP.LT.EPSILON(2))IFLAG=0
      GO TO 999
!
!             RELATIVE CONVERGENCE CRITERION
    3 IF(TEMP/ABS(VAL(2)).LT.EPSILON(1))IFLAG=0
!
  999 RETURN
      END Subroutine QXZ255

      SUBROUTINE QXZ275(X,Y,W,XK,PK,IPK,DEV,LM1,N,L,ITMAX,XR,R,NN,NN1,  &
     &         C,ANS,L2,LM2,LL,E,D,GG,G,B,XLAM,TEMP3,IW,IERR,NKMAX)
!***********************************************************************
!*                                                                     *
!*    PURPOSE:                                                         *
!*                 SUBROUTINE QXZ275 (RASP) FITS A SMOOTH CUBIC SPLINE *
!*                 TO A UNIVARIATE FUNCTION.  DATA MAY BE UNEQUALLY    *
!*                 SPACED.                                             *
!*                                                                     *
!*                                                                     *
!*        USE:                                                         *
!*                 CALL QXZ275(X,Y,W,XK,PK,IPK,DEV,LM1,N,L,ITMAX,XR,R, *
!*                             NN,NN1,C,ANS,L2,LM2,LL,E,D,GG,G,B,XLAM, *
!*                             TEMP3,IW,IERR,NKMAX)                    *
!*                                                                     *
!*         X       AN INPUT ONE-DIMENSIONAL REAL ARRAY DIMENSIONED AT  *
!*                 LEAST N IN THE CALLING PROGRAM.  X CONTAINS THE     *
!*                 X-COORDINATES OF THE DATA POINTS.                   *
!*                                                                     *
!*         Y       AN INPUT ONE-DIMENSIONAL REAL ARRAY WHICH           *
!*                 CONTAINS THE FUNCTION VALUES OF THE DATA POINTS.    *
!*                 Y MUST BE DIMENSIONED AT LEAST N.                   *
!*                                                                     *
!*         W       AN INPUT ONE-DIMENSIONAL REAL ARRAY WHICH CONTAINS  *
!*                 THE WEIGHTING VALUES.  W MUST BE DIMENSIONED AT     *
!*                 LEAST N.                                            *
!*                                                                     *
!*         XK      AN INPUT ONE-DIMENSIONAL REAL ARRAY WHICH           *
!*                 CONTAINS THE X-COORDINATES OF THE SELECTED KNOTS.   *
!*                 XK MUST BE DIMENSIONED AT LEAST L.                  *
!*                                                                     *
!*         PK      AN INPUT ONE-DIMENSIONAL REAL ARRAY WHICH CONTAINS  *
!*                 THE TENSION FACTORS APPLIED ON THE INTERVALS        *
!*                 BETWEEN KNOTS.  PK MUST BE DIMENSIONED AT           *
!*                 LEAST L-1.  PK(K) MUST BE GREATER THAN -1.          *
!*                                                                     *
!*         IPK     AN INPUT ONE-DIMENSIONAL INTEGER ARRAY WHICH        *
!*                 CONTAINS THE TENSION FACTOR ADJUSTMENT FLAGS.       *
!*                 IPK MUST BE DIMENSIONED AT LEAST L-1.               *
!*                                                                     *
!*                 IPK(K)=0  FACTOR K MAY BE ADJUSTED.                 *
!*                                                                     *
!*                 IPK(K)=1  FACTOR K IS HELD CONSTANT.                *
!*                                                                     *
!*         DEV     AN INPUT ONE-DIMENSIONAL REAL ARRAY WHICH           *
!*                 CONTAINS THE ALLOWED PERCENTAGE OF DEVIATION OF     *
!*                 THE APPROXIMATED SOLUTION FROM THE LINEAR SOLUTION. *
!*                 DEV MUST BE DIMENSIONED AT LEAST L-1.               *
!*                                                                     *
!*         LM1     AN INPUT INTEGER WHICH IS L-1                       *
!*                                                                     *
!*         N       AN INPUT INTEGER SPECIFYING THE NUMBER OF DATA      *
!*                 POINTS FOR THE INDEPENDENT VARIABLE.                *
!*                                                                     *
!*         L       AN INPUT INTEGER SPECIFYING THE NUMBER OF KNOTS.    *
!*                 KNOTS MUST BE CHOSEN SO THAT THERE ARE AT LEAST     *
!*                 THREE DATA POINTS IN EACH INTERVAL DEFINED BY       *
!*                 CONSECUTIVE PAIRS OF KNOTS.                         *
!*                                                                     *
!*         ITMAX   AN INPUT INTEGER SPECIFYING THE MAXIMUM NUMBER      *
!*                 OF ITERATIONS PERMITTED IN ACHIEVING THE ALLOWED    *
!*                 DEVIATIONS.   ITMAX  35.                            *
!*                                                                     *
!*         XR      AN INPUT ONE-DIMENSIONAL REAL ARRAY WHICH CONTAINS  *
!*                 THE X-COORDINATES OF THE POINTS AT WHICH INTER-     *
!*                 POLATION IS DESIRED.  XR MUST BE DIMENSIONED AT     *
!*                 LEAST ABS(NN).  IF NN IS POSITIVE, ALL ELEMENTS OF  *
!*                 XR MUST CONTAIN VALUES.  IF NN IS NEGATIVE, ONLY    *
!*                 XR(1) AND XR(ABS(NN)) NEED CONTAIN VALUES.          *
!*                                                                     *
!*         R       AN OUTPUT REAL ARRAY WHICH CONTAINS THE             *
!*                 RATIONAL SPLINE APPROXIMATIONS TO THE DATA.         *
!*                 R MUST BE DIMENSIONED AT LEAST N.                   *
!*                                                                     *
!*         NN      AN INPUT INTEGER. ABS(NN) IS  THE NUMBER            *
!*                 OF INTERPOLATED VALUES TO BE CALCULATED.            *
!*                                                                     *
!*         NN1     AN INPUT INTEGER WHICH IS ABS(NN).                  *
!*
!*         C       AN OUTPUT TWO-DIMENSIONAL REAL ARRAY WHICH          *
!*                 CONTAINS THE CORRELATION MATRIX.  C SHOULD HAVE     *
!*                 A FIRST DIMENSION OF L2 AND A SECOND DIMENSION      *
!*                 OF AT LEAST 2*L.                                    *
!*                                                                     *
!*         ANS     A ONE-DIMENSIONAL REAL ARRAY WHICH CONTAINS THE     *
!*                 Y-COORDINATES AND SECOND DERIVATIVES OF THE KNOTS.  *
!*                                                                     *
!*         L2      AN INTEGER WHICH IS L*2                             *
!*                                                                     *
!*         LM2     AN INTEGER WHICH IS L-2                             *
!*                                                                     *
!*         LL      AN INTEGER WHICH IS 2*L +1                          *
!*                                                                     *
!*         E       AN OUTPUT TWO-DIMENSIONAL REAL ARRAY WHICH CONTAINS *
!*                 ((ETWE)INVERSE)ETWY ON OUTPUT (ET=E TRANSPOSE)      *
!*                                                                     *
!*         D       AN OUTPUT TWO-DIMENSIONAL REAL ARRAY WHICH CONTAINS *
!*                 (ETWE)INVERSE ON OUTPUT.                            *
!*                                                                     *
!*         GG      AN OUTPUT TWO-DIMENSIONAL REAL ARRAY WHICH CONTAINS *
!*                 ((S(ETWE)INVERSE)ST)INVERSE ON OUTPUT.              *
!*                                                                     *
!*         G       AN OUTPUT TWO-DIMENSIONAL REAL ARRAY WHICH CONTAINS *
!*                 (S(ETWE)INVERSE)ST ON OUTPUT                        *
!*                                                                     *
!*         B       ELEMENTS OF MATRIX S, COEFFICIENTS OF CONSTRAINT    *
!*                 EQUATION 15, PAGE 2.                                *
!*                                                                     *
!*                                                                     *
!*         XLAM    AN OUTPUT ONE-DIMENSIONAL ARRAY WHICH CONTAINS      *
!*                 THE LAGRANGE MULTIPLIERS.                           *
!*                                                                     *
!*         TEMP3   A TWO-DIMENSIONAL REAL ARRAY USED AS TEMPORARY      *
!*                 STORAGE IN COVARIANCE MATRIX CALCULATIONS.          *
!*                 TEMP3= -B*BB*B(INVERSE)*E                           *
!*                                                                     *
!*         IW      AN INPUT/OUTPUT INTEGER.                            *
!*                                                                     *
!*                 INPUT   ON THE FIRST CALL TO RSUTDS ,               *
!*                         IW SHOULD BE SET TO 0.  ON SUBSEQUENT CALLS,*
!*                         SHOULD BE SET TO 1 TO PROVIDE ADDITIONAL    *
!*                         INTERPOLATION OR DIFFERENTIATION ON THE SAME*
!*                         FITTING CALCULATIONS.                       *
!*                                                                     *
!*                 OUTPUT  IW IS AUTOMATICALLY SET TO 1.               *
!*                                                                     *
!*         IERR    OUTPUT ERROR PARAMETER:                             *
!*                                                                     *
!*                 = 0  NORMAL RETURN.  NO ERROR DETECTED.             *
!*                 =-2  DID NOT CONVERGE IN ALLOWED NUMBER OF ITERAT-  *
!*                      IONS.                                          *
!*                 =-3  NO SOLUTION DETERMINED BECAUSE THE MATRIX      *
!*                      SENT TO SYMPDS WAS NOT SPD OR SECBI FAILED.    *
!*                                                                     *
!*                                                                     *
!*         NKMAX   AN INPUT INTEGER WHICH AGREES EXACTLY WITH THE      *
!*                 FIRST DIMENSION OF ARRAY C. NKMAX MUST BE GREATER   *
!*                 THAN OR EQUAL TO 2*NK.                              *
!*                                                                     *
!*                                                                     *
!*                                                                     *
!*                                                                     *
!*    REQUIRED ROUTINES                QXZ063,QXZ276,QXZ277,QXZ278,    *
!*                                     QXZ279,QXZ280,SYMPDS(F4.3),     *
!*                                     SECBI(C2.5)                     *
!*                                                                     *
!*    SOURCE                                                           *
!*                                     PROGRAM RASPFIT BY JAMES R.     *
!*                                     SCHIESS AND PATRICIA A. KERR    *
!*                                     MODIFIED BY COMPUTER SCIENCES   *
!*                                     CORPORATION                     *
!*                                                                     *
!*    LANGUAGE                         -FORTRAN                        *
!*                                                                     *
!*    DATE RELEASED                    FEBRUARY, 1986                  *
!*                                                                     *
!*    LATEST REVISION                  FEBRUARY, 1986                  *
!***********************************************************************
!
!
!
!     PARAMETERS USED
!
!
      DIMENSION X(N),Y(N),W(N),DEV(LM1),IPK(LM1),PK(LM1)
      DIMENSION XK(L),R(NN1),XR(NN1)
      DIMENSION C(NKMAX,L2),ANS(L2),E(L2,LL),B(LM2,L2),GG(LM2,LM2)
!
!     THE FIRST DIMENSION OF C WAS CHANGED FROM L2 TO NKMAX
!
      DIMENSION D(L2,L2),G(LM2,LM2),XLAM(LM2)
      DIMENSION TEMP3(L2,L2)
!
!!!      INTEGER QXZ063
      IF(IW.EQ.1)GO TO 70
!
!
      CALL QXZ278(X,N,XK,L,ISTOP)
      IERR=-ISTOP
      IF(IERR .NE. 0)  GO TO 120
      IPKZERO=LM1
!
      DO 30 I=1,LM1
   30 IPKZERO=IPKZERO-IPK(I)
!
      DO 50 J=1,ITMAX
      CALL QXZ277(X,Y,W,N,XK,L,PK,LM1,ANS,L2,E,GG,B,LM2,LL,D,XLAM,G,    &
     &IERR)
      IF(IERR.NE.0)GO TO 120
      CALL QXZ279(IPK,XK,PK,DEV,LM1,ANS,L,L2,IPKZERO,IERR)
      IF(IERR.NE.0)GO TO 120
      IF(IPKZERO .EQ. 0)GO TO 70
   50 END DO
!
!              CALCULATE STANDARD ERROR
!
      IERR=-2
   70 STDERR=0.
      FLGRP=0.
      FLGRPP=0.
      K = 1
      IF(NN.GT.0)GO TO 960
      I=QXZ063(XR(1),XK,L)
      K=I
      DXK = XK(K+1) - XK(K)
      CALL QXZ276(DXK,PK(K),ANS(2*K-1),ANS(2*K),ANS(2*K+1),ANS(2*K+2),  &
     &     XK(K),XR(1),R(1),FLGRP,FLGRPP)
!
      IF(XR(1).GT.XR(NN1))GO TO 900

      DO 940 J=2,NN1
   40 IF(XR(J).GT.XK(I+1))GO TO 200
      K=I
      DXK = XK(K+1) - XK(K)
      CALL QXZ276(DXK,PK(K),ANS(2*K-1),ANS(2*K),ANS(2*K+1),ANS(2*K+2),  &
     &     XK(K),XR(J),R(J),FLGRP,FLGRPP)
      CYCLE   !!! GO TO 940
  200 I=I+1
      GO TO 40
  940 END DO
      GO TO 980

  900 DO 920 J=2,NN1
  930 IF(XR(J).LT.XK(I))GO TO 910
      K=I
      DXK = XK(K+1) - XK(K)
      CALL QXZ276(DXK,PK(K),ANS(2*K-1),ANS(2*K),ANS(2*K+1),ANS(2*K+2),  &
     &     XK(K),XR(J),R(J),FLGRP,FLGRPP)
      CYCLE   !!! GO TO 920
  910 I=I-1
      GO TO 930
  920 END DO

      GO TO 980
  960 DO 280 I=1,NN
      K=QXZ063(XR(I),XK,L)
      DXK = XK(K+1) - XK(K)
      CALL QXZ276(DXK,PK(K),ANS(2*K-1),ANS(2*K),ANS(2*K+1),ANS(2*K+2),  &
     &     XK(K),XR(I),R(I),FLGRP,FLGRPP)
  280  CONTINUE
!
  980 IF(IW.EQ.1)GO TO 120
      IW=1
      K=1
      DO 80 I=1,N
      IF(X(I).GT.XK(K+1))K=K+1
      DXK = XK(K+1) - XK(K)
      CALL QXZ276(DXK,PK(K),ANS(2*K-1),ANS(2*K),ANS(2*K+1),ANS(2*K+2),  &
     &     XK(K),X(I),F,FLGRP,FLGRPP)
      STDERR= STDERR +(Y(I)-F)**2
   80 END DO
      STDERR= SQRT(STDERR/(N-2*L))
!
!
!      CALCULATE COVARIANCE MATRIX
!
      DO 230 KK=1,L2
      DO 230 II=1,L2
      TEMP3(II,KK) = 0.0
      DO 230 JJ=1,L2
      DO 230  J=1,LM2
      DO 230  M=1,LM2
  230 TEMP3(II,KK)=TEMP3(II,KK)-B(M,II)*GG(M,J)*B(J,JJ)*E(JJ,KK+1)
      DO 260 III=1,L2
      DO 260 KKK=1,L2
      C(III,KKK) = 0.0
      IF(III .EQ. KKK)  C(III,KKK) = C(III,KKK)+E(III,III+1)
      DO 260 JJJ=1,L2
  260 C(III,KKK) = C(III,KKK) + E(III,JJJ+1)*TEMP3(JJJ,KKK)
!
!     CALCULATE CORRELATION MATRIX WITH STANDARD DEVIATIONS ON DIAGONAL
!
      DO 100 I=1,L2
      TEMP = SQRT(STDERR*C(I,I))
      DO 110 J=1,L2
      C(I,J) = STDERR*C(I,J)/TEMP
      C(J,I) = STDERR*C(J,I)/TEMP
  110 END DO
      C(I,I) = TEMP
  100 END DO
  120 RETURN
      END Subroutine QXZ275

      SUBROUTINE QXZ276(DX,PK,YBAR,YBARPP,YBAR1,YBAR1PP,XK,X,F,FP,FPP)
!***********************************************************************
!*                                                                     *
!*    PURPOSE:                                                         *
!*                 SUBROUTINE QXZ276 (FRATS) CALCULATES THE            *
!                  RATIONAL SPLINE APPROXIMATIONS OF F,F',F'',         *
!*                 BUT RETURNS ONLY F IN RSUTDS.                       *
!*                                                                     *
!*         DX      LENGTH OF SUBINTERVAL DEFINED BY TWO CONSECUTIVE    *
!*                 KNOT ABSCISSAS (XK(K+1) - X(K))                     *
!*                                                                     *
!*         PK      A ONE-DIMENSIONAL REAL ARRAY WHICH CONTAINS THE     *
!*                 TENSION FACTORS                                     *
!*                                                                     *
!*         YBAR    Y-COORDINATE OF KNOT XK(J)                          *
!*                                                                     *
!*         YBARPP  SECOND DERIVATIVE OF KNOT XK(J)                     *
!*                                                                     *
!*         YBAR1   Y-COORDINATE OF KNOT XK(J+1)                        *
!*                                                                     *
!*         YBAR1PP SECOND DERIVATIVE AT KNOT XK(J+1)                   *
!*                                                                     *
!*         XK      X-COORDINATE OF KNOT XK(J)                          *
!*                                                                     *
!*         X       X-COORDINATE OF DATA POINT                          *
!*                                                                     *
!*         F       RATIONAL SPLINE APPROXIMATION OF FUNCTION VALUE     *
!*                                                                     *
!*         FPP     RATIONAL SPLINE APPROXIMATION OF FIRST DERIVATIVE   *
!*                                                                     *
!*         FPP     RATIONAL SPLINE APPROXIMATION OF SECOND DERIVATIVE  *
!                                                                      *
!***********************************************************************
      T=(X-XK)/DX
      U=1.-T
      PT=PK*T+1.
      PU=PK*U+1.
      HK = DX**2/(2.*(PK**2+3.*PK+3.))
!
      F=U*YBAR + HK*(U**3/PT-U)*YBARPP + T*YBAR1 +HK*(T**3/PU-T)*YBAR1PP
      GO TO 30
!
   10 IF(FP .LT. 1.0)  GO TO 20
      FP=(-YBAR-HK*(((3.*U**2*PT+U**3*PK)/PT**2)-1.)*YBARPP+YBAR1+HK*   &
     &  (((3.*T**2*PU+T**3*PK)/PU**2)-1.)*YBAR1PP)/DX
!
   20 IF(FPP .LT. 1.0)  GO TO 30
      FPP=HK*YBARPP*(2.*PK**2*U**3+6.*PK*PT*U**2+6.*PT**2*U)/PT**3      &
     &   +HK*YBAR1PP*(2.*PK**2*T**3+6.*PK*PU*T**2+6.*PU**2*T)/PU**3
      FPP=FPP/DX**2
!
   30 RETURN
      END Subroutine QXZ276

      SUBROUTINE QXZ277(X,Y,W,N,XK,L,PK,LM1,A,L2,E,GG,B,LM2,LLL,D,      &
     &      XLAM,G,IERR)
!***********************************************************************
!*                                                                     *
!*    PURPOSE:                                                         *
!*                 SUBROUTINE QXZ277 (LSFIT) IS USED                   *
!*                 TO FIND THE LEAST SQUARES FIT FOR RSUTDS            *
!*                                                                     *
!*                                                                     *
!*        USE:                                                         *
!*                 CALL QXZ277(X,Y,W,N,XK,L,PK,LM1,A,L2,E,GG,B,LM2,    *
!*                            LLL,D,XLAM,G,IERR)                       *
!*                                                                     *
!*         X       AN INPUT ONE-DIMENSIONAL ARRAY DIMENSIONED AT LEAST *
!*                 N IN THE CALLING PROGRAM.  X CONTAINS THE           *
!*                 X-COORDINATES OF THE DATA POINTS.                   *
!*                                                                     *
!*         Y       AN INPUT ONE-DIMENSIONAL REAL ARRAY WHICH           *
!*                 CONTAINS THE FUNCTION VALUES OF THE DATA POINTS.    *
!*                 Y MUST BE DIMENSIONED AT LEAST N.                   *
!*                                                                     *
!*         W       AN INPUT ONE-DIMENSIONAL REAL ARRAY WHICH CONTAINS  *
!*                 THE WEIGHTING VALUES.  W MUST BE DIMENSIONED AT     *
!*                 LEAST N.                                            *
!*                                                                     *
!*         N       AN INPUT INTEGER SPECIFYING THE NUMBER OF DATA      *
!*                 POINTS FOR THE INDEPENDENT VARIABLE.                *
!*                                                                     *
!*         XK      AN INPUT ONE-DIMENSIONAL REAL ARRAY WHICH           *
!*                 CONTAINS THE X-COORDINATES OF THE SELECTED KNOTS.   *
!*                 XK MUST BE DIMENSIONED AT LEAST L.                  *
!*                                                                     *
!*         L       AN INPUT INTEGER SPECIFYING THE NUMBER OF KNOTS.    *
!*                 KNOTS MUST BE CHOSEN SO THAT THERE ARE AT LEAST     *
!*                 THREE DATA POINTS IN EACH INTERVAL DEFINED BY       *
!*                 CONSECUTIVE PAIRS OF KNOTS.                         *
!*                                                                     *
!*                                                                     *
!*         PK      AN INPUT ONE-DIMENSIONAL REAL ARRAY WHICH CONTAINS  *
!*                 THE TENSION FACTORS APPLIED ON THE INTERVALS        *
!*                 BETWEEN KNOTS.  PK MUST BE DIMENSIONED AT           *
!*                 LEAST L-1.  PK(K) MUST BE GREATER THAN -1.          *
!*                 PK(K) MUST BE GREATER THAN -1.                      *
!*                                                                     *
!*         LM1     AN INTEGER WHICH IS L-1                             *
!*                                                                     *
!*         A       AN OUTPUT ONE-DIMENSIONAL ARRAY WHICH               *
!*                 CONTAINS THE FUNCTION VALUES AND SECOND DERIVATIVES *
!*                 OF THE DATA.                                        *
!*                                                                     *
!*         L2      AN INTEGER WHICH IS L*2                             *
!*                                                                     *
!*         E       AN OUTPUT TWO-DIMENSIONAL REAL ARRAY WHICH CONTAINS *
!*                 ((ETWE)INVERSE)ETWY ON OUTPUT (ET=E TRANSPOSE)      *
!*                                                                     *
!*                                                                     *
!*         GG      AN OUTPUT TWO-DIMENSIONAL REAL ARRAY WHICH CONTAINS *
!*                 ((S(ETWE)INVERSE)ST)INVERSE ON OUTPUT.              *
!*                                                                     *
!*         B       ELEMENTS OF MATRIX S, COEFFICIENTS OF CONSTRAINT    *
!*                 EQUATION 15, PAGE 2.                                *
!*                                                                     *
!*         LM2     AN INTEGER WHICH IS L-2                             *
!*                                                                     *
!*         LLL     AN INTEGER WHICH IS 2*L +1                          *
!*                                                                     *
!*         D       AN OUTPUT TWO-DIMENSIONAL REAL ARRAY WHICH CONTAINS *
!*                 (ETWE)INVERSE ON OUTPUT.                            *
!*                                                                     *
!*                                                                     *
!*         XLAM    AN OUTPUT ONE-DIMENSIONAL ARRAY WHICH CONTAINS      *
!*                 THE LAGRANGE MULTIPLIERS.                           *
!*                                                                     *
!*         G       AN OUTPUT TWO-DIMENSIONAL REAL ARRAY WHICH CONTAINS *
!*                 (S(ETWE)INVERSE)ST ON OUTPUT                        *
!*                                                                     *
!*         IERR    OUTPUT ERROR PARAMETER:                             *
!*                                                                     *
!*                 = 0  NORMAL RETURN.  NO ERROR DETECTED.             *
!*                 =-2  DID NOT CONVERGE IN ALLOWED NUMBER OF ITERAT-  *
!*                      IONS.                                          *
!*                 =-3  NO SOLUTION DETERMINED BECAUSE THE MATRIX      *
!*                      SENT TO SYMPDS WAS NOT SPD OR SECBI FAILED.    *
!*                                                                     *
!*                                                                     *
!***********************************************************************
!
      DIMENSION GG(LM2,LM2),XLAM(LM2),B(LM2,L2)
      DIMENSION X(N),Y(N),W(N),PK(LM1),XK(L)
      DIMENSION  A(L2),D(L2,L2),G(LM2,LM2),E(L2,LLL)
!
!     ZERO  ARRAYS
!
      DO 5 I=1,L
      A(I) = 0.0
      A(L+I) = 0.0
      DO 6 J=1,L
      D(I,J) = 0.0
      D(L+I,J) = 0.0
      D(I,L+J) = 0.0
      D(L+I,L+J) = 0.0
      E(I,J) = 0.0
      E(L+I,J) = 0.0
      E(I,L+J) = 0.0
      E(L+I,L+J) = 0.0
    6 END DO
      E(I,2*L+1) = 0.0
      E(L+I,2*L+1) = 0.0
      E(I,I+1) = 1.
      E(L+I,L+I+1) = 1.
    5 END DO
!
!     LEAST SQUARES SET-UP
!
      K = 1
      DO 40 I=1,N
      IF(X(I) .GT. XK(K+1)) K = K + 1
      DXK = XK(K+1) - XK(K)
      HK=(DXK**2)/(2.*(PK(K)**2 +3.*PK(K)+3.))
      KU=2*(K+1)
      KL=KU-3
      T=(X(I)-XK(K))/DXK
      U=1.-T
      A(KL)=U
      A(KL+1)= (U**3/(PK(K)*T+1.)-U)*HK
      A(KL+2)= T
      A(KU)= (T**3/(PK(K)*U+1.)-T)*HK

      DO 30 KK=KL,KU
      E(KK,1)= E(KK,1)+A(KK)*Y(I)*W(I)
      DO 31 LL=KL,KU
        D(KK,LL)=D(KK,LL) +A(KK)*W(I)*A(LL)
   31 END DO
   30 END DO
   40 END DO
!
!     LEAST SQUARES SOLUTION
!
      IOPT=1
      IFAC=0
      CALL SYMPDS(L2,L2,D,LLL,E,IOPT,IFAC,DET,ISC,A,IERR)
      IF(IERR.NE.0)IERR=-5
      IF(IERR.NE.0)GO TO 140
!
!
!     SET UP CONSTRAINTS
!
      DO 60 I=1,L2
      DO 60 K=1,LM2
   60 B(K,I)=0.
      HK=((XK(2)-XK(1))**2)*0.5/(PK(1)**2 +3.*PK(1)+3.)
!
      DO 70 K=1,LM2
      DXK = XK(K+1) - XK(K)
      DXK1 = XK(K+2) - XK(K+1)
      KL=2*K-1
      HK1=0.5*(DXK1**2)/(PK(K+1)**2 +3.*PK(K+1)+3.)
      B(K,KL)=-1./DXK
      B(K,KL+1)=HK/DXK
      B(K,KL+2)=1./DXK +1./DXK1
      B(K,KL+3)=(PK(K)+2.)*HK/DXK +(PK(K+1)+2.)*HK1/DXK1
      B(K,KL+4)= -1./DXK1
      B(K,KL+5)= HK1/DXK1
      HK=HK1
   70 END DO
!
      DO 100 K=1,LM2
      DO 90 I=1,LM2
      G(K,I)=0.
      GG(I,K)=0.
      DO 80 KK=1,L2
      DO 80 LL=1,L2
   80 G(K,I) = G(K,I)+B(K,KK)*E(KK,LL+1)*B(I,LL)
   90 END DO
      GG(K,K)=1.
  100 END DO
!
!     CONSTRAINED  SOLUTION
!
      IOPT=1
      IFAC=0
      CALL SYMPDS(LM2,LM2,G,LM2,GG,IOPT,IFAC,DET,ISC,A,IERR)
      IF(IERR.NE.0)IERR=-5
      IF(IERR .NE. 0) GO TO 140
!
!     CALCULATE XLAM FOR USE WITH LAGRANGE MULTIPLIER METHOD
!
      DO 120 I=1,LM2
      XLAM(I)=0.
      DO 120 J=1,LM2
      DO 120 KK=1,L2
  120 XLAM(I)= XLAM(I) +GG(I,J)*B(J,KK)*E(KK,1)
!
!     UPDATE SOLUTION WITH CONSTRAINTS
!
      DO 130 I=1,L2
      A(I) = E(I,1)
      DO 130 J=1,L2
      DO 130 KK=1,LM2
  130 A(I)=A(I)-E(I,J+1)*B(KK,J)*XLAM(KK)
!
  140 RETURN
      END Subroutine QXZ277

      SUBROUTINE QXZ278(X,N,XK,L,ISTOP)
!***********************************************************************
!*                                                                     *
!*    PURPOSE:                                                         *
!*                 SUBROUTINE QXZ278 (INTERVL) IS USED                 *
!*                 TO DETERMINE WHETHER THERE ARE AT LEAST 3 POINTS IN *
!*                 EACH INTERVAL.                                      *
!*                                                                     *
!*                                                                     *
!*        USE:                                                         *
!*                 CALL QXZ278(X,N,XK,L,ISTOP)                         *
!*                                                                     *
!*         X       AN INPUT ONE-DIMENSIONAL ARRAY DIMENSIONED AT LEAST *
!*                 N IN THE CALLING PROGRAM.  X CONTAINS THE           *
!*                 X-COORDINATES OF THE DATA POINTS.                   *
!*                                                                     *
!*         N       AN INPUT INTEGER SPECIFYING THE NUMBER OF DATA      *
!*                 POINTS FOR THE INDEPENDENT VARIABLE.                *
!*                                                                     *
!*         XK      AN INPUT ONE-DIMENSIONAL REAL ARRAY WHICH           *
!*                 CONTAINS THE X-COORDINATES OF THE SELECTED KNOTS.   *
!*                 XK MUST BE DIMENSIONED AT LEAST L.                  *
!*                                                                     *
!*         L       AN INPUT INTEGER SPECIFYING THE NUMBER OF KNOTS.    *
!*                 KNOTS MUST BE CHOSEN SO THAT THERE ARE AT LEAST     *
!*                 THREE DATA POINTS IN EACH INTERVAL DEFINED BY       *
!*                 CONSECUTIVE PAIRS OF KNOTS.                         *
!*                                                                     *
!*         ISTOP   INTEGER FLAG SET TO 1 IF THERE ARE LESS THAN        *
!*                 3 DATA POINTS IN ANY SUBINTERVAL.                   *
!*                                                                     *
!***********************************************************************
!
      DIMENSION X(N),XK(L)
!
      K=1
      ICOUNT = 0
      ISTOP = 0
      DO 10 J=2,N-1
      ICOUNT = ICOUNT + 1
      IF(X(J) .LT. XK(K+1)) CYCLE   !!!  GO TO 10
      IF(ICOUNT .LT. 3)  ISTOP=1
      ICOUNT = 1
      K=K+1
   10 END DO
      IF(ICOUNT .LT.3)  ISTOP=1
!
      RETURN
      END Subroutine QXZ278

      SUBROUTINE QXZ279(IPK,XK,PK,DEV,LM1,ANS,L,L2,IPKZERO,IERR)
!
!***********************************************************************
!*                                                                     *
!*    PURPOSE:                                                         *
!*                 SUBROUTINE QXZ279 (PKADJ) IS USED                   *
!*                 TO INCREMENT EACH PK IF ERROR EXCEEDS               *
!*                 DESIRED DEVIATION IN INTERVAL                       *
!*                                                                     *
!*                                                                     *
!*        USE:                                                         *
!*                 CALL QXZ279(IPK,XK,PK,DEV,LM1,ANS,L,L2,IPKZERO,IERR)*
!*                                                                     *
!*         IPK     AN INPUT ONE-DIMENSIONAL INTEGER ARRAY WHICH        *
!*                 CONTAINS THE TENSION FACTOR ADJUSTMENT FLAGS.       *
!*                                                                     *
!*                 IPK(K)=0  FACTOR K MAY BE ADJUSTED.                 *
!*                                                                     *
!*                 IPK(K)=1  FACTOR K IS HELD CONSTANT.                *
!*                                                                     *
!*         XK      AN INPUT ONE-DIMENSIONAL REAL ARRAY WHICH           *
!*                 CONTAINS THE X-COORDINATES OF THE SELECTED KNOTS.   *
!*                 XK MUST BE DIMENSIONED AT LEAST L.                  *
!*                                                                     *
!*         PK      AN INPUT ONE-DIMENSIONAL REAL ARRAY WHICH CONTAINS  *
!*                 THE TENSION FACTORS APPLIED ON THE INTERVALS        *
!*                 BETWEEN KNOTS.  PK MUST BE DIMENSIONED AT           *
!*                 LEAST L-1.  PK(K) MAY BE UNDEFINED IF IPK(K)=0.     *
!*                 PK(K) MUST BE GREATER THAN -1.                      *
!*                                                                     *
!*         DEV     AN INPUT ONE-DIMENSIONAL REAL ARRAY WHICH           *
!*                 CONTAINS THE ALLOWED PERCENTAGE OF DEVIATION OF     *
!*                 THE APPROXIMATED SOLUTION FROM THE LINEAR SOLUTION. *
!*                 DEV MUST BE DIMENSIONED AT LEAST L-1.               *
!*                                                                     *
!*         LM1     AN INTEGER WHICH IS L-1                             *
!*                                                                     *
!*         ANS     AN INPUT ONE-DIMENSIONAL REAL ARRAY WHICH HOLDS     *
!*                 Y-COORDINATES AND SECOND DERIVATIVES AT THE KNOTS.  *
!*                 ANS MUST BE DIMENSIONED L2.                         *
!*                                                                     *
!*         L       AN INPUT INTEGER SPECIFYING THE NUMBER OF KNOTS.    *
!*                 KNOTS MUST BE CHOSEN SO THAT THERE ARE AT LEAST     *
!*                 THREE DATA POINTS IN EACH INTERVAL DEFINED BY       *
!*                 CONSECUTIVE PAIRS OF KNOTS.                         *
!*                                                                     *
!*         L2      AN INTEGER WHICH IS L*2                             *
!*                                                                     *
!*         IPKZERO NUMBER OF TENSION FACTORS TO BE ADJUSTED            *
!*                                                                     *
!*         IERR    OUTPUT ERROR PARAMETER:                             *
!*                                                                     *
!*                 = 0  NORMAL RETURN.  NO ERROR DETECTED.             *
!*                 =-2  DID NOT CONVERGE IN ALLOWED NUMBER OF ITERAT-  *
!*                      IONS.                                          *
!*                 =-3  NO SOLUTION DETERMINED BECAUSE THE MATRIX      *
!*                      SENT TO SYMPDS WAS NOT SPD OR SECBI FAILED.    *
!*                                                                     *
!***********************************************************************
      COMMON /INTRVAL/KK,PKK,ANS1,ANS2,ANS3,ANS4,DXK,DYK
!
      DIMENSION EPS(3)
      DIMENSION XK(L),IPK(LM1),PK(LM1),DEV(LM1),ANS(L2)
!
!!!      EXTERNAL QXZ280
!
!     INCREMENTS EACH PK IF ERROR EXCEEDS DESIRED DEVIATION IN INTERVAL
!
      DATA T1,T2,DT/0.,1.,0.125/
      DATA (EPS(J),J=1,3)/1.0E-8,1.0E-12,1.0E-8/
!
      FLGRP = 0.0
      FLGRPP = 0.0

      DO 100 K=1,LM1
      IF(IPK(K) .EQ.1) CYCLE   !!! GO TO 100
      KK=K
      PKK = PK(K)
      DXK = XK(K+1) - XK(K)
      DYK = ANS(2*K+1) - ANS(2*K-1)
      ANS1 = ANS(2*K-1)
      ANS2 = ANS(2*K)
      ANS3 = ANS(2*K+1)
      ANS4 = ANS(2*K+2)
      IF(ANS2.EQ.0.0.AND.ANS4.EQ.0.0)GO TO 70
      CALL SECBI(T1,T2,DT,QXZ280,EPS,T,IERR)
      IF(IERR.NE.0.OR.T.EQ.0.0.OR.DXK.EQ.0.0)IERR=-3
      IF(IERR .NE. 0) CYCLE   !!! GO TO 100
   10 XX = DXK*T + XK(K)
      YY=2.
      CALL QXZ276(DXK,PKK,ANS1,ANS2,ANS3,ANS4,XK(K),XX,YY,FLGRP,FLGRPP)
      XM1=DYK/DXK
      XM2=(YY-ANS(2*K-1))/(XX-XK(K))
      THET= ATAN2( (XM2-XM1),1.+XM1*XM2)
      EL2=SQRT((XX-XK(K))**2 +(YY-ANS(2*K-1))**2)
      EL1=SQRT(DXK**2 +DYK**2)
      DIS=EL2*100.*ABS(SIN(THET))/EL1
!
      IF(DIS .LE. DEV(K)) GO TO 70
      PK(K) = PK(K)+1.0
      CYCLE   !!! GO TO 100
   70 CONTINUE
      IPK(K)=1
      IPKZERO = IPKZERO -1
  100 END DO
!
      RETURN
      END Subroutine QXZ279

      FUNCTION QXZ280(T)
!***********************************************************************
!*                                                                     *
!     FUNCTION (FUNC) WHICH MEASURES  DEVIATION OF SPLINE FROM LINE    *
!*                                                                     *
!***********************************************************************
      COMMON /INTRVAL/K,PKK,YK,YKPP,YKP1,YKP1PP,DXK,DYK
      U=1.-T
      PT=PKK*T+1.
      PU=PKK*U+1.
      HK = DXK**2/(2.*(PKK**2+3.*PKK+3.))
      QXZ280=(-YK-HK*((3.*U**2*PT+U**3*PKK)/PT**2-1.)*YKPP              &
     &      +YKP1+HK*((3.*T**2*PU+T**3*PKK)/PU**2-1.)*YKP1PP)           &
     &    -DYK
!
      RETURN
      END Function QXZ280

      SUBROUTINE QXZ281(X,Y,W,XK,PK,IPK,DEV,LM1,N,L,ITMAX,XR,R,NN,NN1,  &
     & DER1,DER2,C,ANS,L2,LM2,LL,E,D,GG,G,B,XLAM,TEMP3,IW,IERR,NKMAX)
!***********************************************************************
!*                                                                     *
!*    PURPOSE:                                                         *
!*                 SUBROUTINE QXZ281 (RASP) FITS A SMOOTH CUBIC SPLINE *
!*                 TO A UNIVARIATE FUNCTION.  DATA MAY BE UNEQUALLY    *
!*                 SPACED.  FIRST AND SECOND DERIVATIVES ARE ALSO      *
!*                 CALCULATED.                                         *
!*                                                                     *
!*                                                                     *
!*        USE:                                                         *
!*                 CALL QXZ281(X,Y,W,XK,PK,IPK,DEV,LM1,N,L,ITMAX,XR,R, *
!*                             NN,NN1,DER1,DER2,C,ANS,L2,LM2,LL,E,D,GG,*
!*                             G,B,XLAM,TEMP3,IW,IERR,NKMAX)           *
!*                                                                     *
!*         X       AN INPUT ONE-DIMENSIONAL REAL ARRAY DIMENSIONED AT  *
!*                 LEAST N IN THE CALLING PROGRAM.  X CONTAINS THE     *
!*                 X-COORDINATES OF THE DATA POINTS.                   *
!*                                                                     *
!*         Y       AN INPUT ONE-DIMENSIONAL REAL ARRAY WHICH           *
!*                 CONTAINS THE FUNCTION VALUES OF THE DATA POINTS.    *
!*                 Y MUST BE DIMENSIONED AT LEAST N.                   *
!*                                                                     *
!*         W       AN INPUT ONE-DIMENSIONAL REAL ARRAY WHICH CONTAINS  *
!*                 THE WEIGHTING VALUES.  W MUST BE DIMENSIONED AT     *
!*                 LEAST N.                                            *
!*                                                                     *
!*         XK      AN INPUT ONE-DIMENSIONAL REAL ARRAY WHICH           *
!*                 CONTAINS THE X-COORDINATES OF THE SELECTED KNOTS.   *
!*                 XK MUST BE DIMENSIONED AT LEAST L.                  *
!*                                                                     *
!*         PK      AN INPUT ONE-DIMENSIONAL REAL ARRAY WHICH CONTAINS  *
!*                 THE TENSION FACTORS APPLIED ON THE INTERVALS        *
!*                 BETWEEN KNOTS.  PK MUST BE DIMENSIONED AT           *
!*                 LEAST L-1.  PK(K) MUST BE GREATER THAN -1.          *
!*                                                                     *
!*         IPK     AN INPUT ONE-DIMENSIONAL INTEGER ARRAY WHICH        *
!*                 CONTAINS THE TENSION FACTOR ADJUSTMENT FLAGS.       *
!*                 IPK MUST BE DIMENSIONED AT LEAST L-1.               *
!*                                                                     *
!*                 IPK(K)=0  FACTOR K MAY BE ADJUSTED.                 *
!*                                                                     *
!*                 IPK(K)=1  FACTOR K IS HELD CONSTANT.                *
!*                                                                     *
!*         DEV     AN INPUT ONE-DIMENSIONAL REAL ARRAY WHICH           *
!*                 CONTAINS THE ALLOWED PERCENTAGE OF DEVIATION OF     *
!*                 THE APPROXIMATED SOLUTION FROM THE LINEAR SOLUTION. *
!*                 DEV MUST BE DIMENSIONED AT LEAST L-1.               *
!*                                                                     *
!*         LM1     AN INPUT INTEGER WHICH IS L-1                       *
!*                                                                     *
!*         N       AN INPUT INTEGER SPECIFYING THE NUMBER OF DATA      *
!*                 POINTS FOR THE INDEPENDENT VARIABLE.                *
!*                                                                     *
!*         L       AN INPUT INTEGER SPECIFYING THE NUMBER OF KNOTS.    *
!*                 KNOTS MUST BE CHOSEN SO THAT THERE ARE AT LEAST     *
!*                 THREE DATA POINTS IN EACH INTERVAL DEFINED BY       *
!*                 CONSECUTIVE PAIRS OF KNOTS.                         *
!*                                                                     *
!*         ITMAX   AN INPUT INTEGER SPECIFYING THE MAXIMUM NUMBER      *
!*                 OF ITERATIONS PERMITTED IN ACHIEVING THE ALLOWED    *
!*                 DEVIATIONS.   ITMAX  35.                            *
!*                                                                     *
!*         XR      AN INPUT ONE-DIMENSIONAL REAL ARRAY WHICH CONTAINS  *
!*                 THE X-COORDINATES OF THE POINTS AT WHICH INTER-     *
!*                 POLATION IS DESIRED.  XR MUST BE DIMENSIONED AT     *
!*                 LEAST ABS(NN).  IF NN IS POSITIVE, ALL ELEMENTS OF  *
!*                 XR MUST CONTAIN VALUES.  IF NN IS NEGATIVE, ONLY    *
!*                 XR(1) AND XR(ABS(NN)) NEED CONTAIN VALUES.          *
!*                                                                     *
!*         R       AN OUTPUT REAL ARRAY WHICH CONTAINS THE             *
!*                 RATIONAL SPLINE APPROXIMATIONS TO THE DATA.         *
!*                 R MUST BE DIMENSIONED AT LEAST N.                   *
!*                                                                     *
!*         NN      AN INPUT INTEGER. ABS(NN) IS  THE NUMBER            *
!*                 OF INTERPOLATED VALUES TO BE CALCULATED.            *
!*                                                                     *
!*         NN1     AN INPUT INTEGER WHICH IS ABS(NN).                  *
!*
!*         DER1     AN OUTPUT ONE-DIMENSIONAL REAL ARRAY WHICH CONTAINS*
!*                 THE VALUES OF THE FIRST DERIVATIVE OF THE INTER-    *
!*                 POLATED DEPENDENT VARIABLE CORRESPONDING TO XR.     *
!*                 DER1 MUST BE DIMENSIONED AT LEAST ABS(NN).          *
!*                                                                     *
!*         DER2     AN OUTPUT ONE-DIMENSIONAL REAL ARRAY WHICH CONTAINS*
!*                 THE VALUES OF THE SECOND DERIVATIVE OF THE INTER-   *
!*                 POLATED DEPENDENT VARIABLE CORRESPONDING TO XR.     *
!*                 DER1 MUST BE DIMENSIONED AT LEAST ABS(NN).          *
!*                                                                     *
!*         C       AN OUTPUT TWO-DIMENSIONAL REAL ARRAY WHICH          *
!*                 CONTAINS THE CORRELATION MATRIX.  C SHOULD HAVE     *
!*                 A FIRST DIMENSION OF L2 AND A SECOND DIMENSION      *
!*                 OF AT LEAST 2*L.                                    *
!*                                                                     *
!*         ANS     A ONE-DIMENSIONAL REAL ARRAY WHICH CONTAINS THE     *
!*                 Y-COORDINATES AND SECOND DERIVATIVES OF THE KNOTS.  *
!*                                                                     *
!*         L2      AN INTEGER WHICH IS L*2                             *
!*                                                                     *
!*         LM2     AN INTEGER WHICH IS L-2                             *
!*                                                                     *
!*         LL      AN INTEGER WHICH IS 2*L +1                          *
!*                                                                     *
!*         E       AN OUTPUT TWO-DIMENSIONAL REAL ARRAY WHICH CONTAINS *
!*                 ((ETWE)INVERSE)ETWY ON OUTPUT (ET=E TRANSPOSE)      *
!*                                                                     *
!*         D       AN OUTPUT TWO-DIMENSIONAL REAL ARRAY WHICH CONTAINS *
!*                 (ETWE)INVERSE ON OUTPUT.                            *
!*                                                                     *
!*         GG      AN OUTPUT TWO-DIMENSIONAL REAL ARRAY WHICH CONTAINS *
!*                 ((S(ETWE)INVERSE)ST)INVERSE ON OUTPUT.              *
!*                                                                     *
!*         G       AN OUTPUT TWO-DIMENSIONAL REAL ARRAY WHICH CONTAINS *
!*                 (S(ETWE)INVERSE)ST ON OUTPUT                        *
!*                                                                     *
!*         B       ELEMENTS OF MATRIX S, COEFFICIENTS OF CONSTRAINT    *
!*                 EQUATION 15, PAGE 2.                                *
!*                                                                     *
!*                                                                     *
!*         XLAM    AN OUTPUT ONE-DIMENSIONAL ARRAY WHICH CONTAINS      *
!*                 THE LAGRANGE MULTIPLIERS.                           *
!*                                                                     *
!*         TEMP3   A TWO-DIMENSIONAL REAL ARRAY USED AS TEMPORARY      *
!*                 STORAGE IN COVARIANCE MATRIX CALCULATIONS.          *
!*                 TEMP3= -B*BB*B(INVERSE)*E                           *
!*                                                                     *
!*         IW      AN INPUT/OUTPUT INTEGER.                            *
!*                                                                     *
!*                 INPUT   ON THE FIRST CALL TO RSDSDR ,               *
!*                         IW SHOULD BE SET TO 0.  ON SUBSEQUENT CALLS,*
!*                         SHOULD BE SET TO 1 TO PROVIDE ADDITIONAL    *
!*                         INTERPOLATION OR DIFFERENTIATION ON THE SAME*
!*                         FITTING CALCULATIONS.                       *
!*                                                                     *
!*                 OUTPUT  IW IS AUTOMATICALLY SET TO 1.               *
!*                                                                     *
!*         IERR    OUTPUT ERROR PARAMETER:                             *
!*                                                                     *
!*                 = 0  NORMAL RETURN.  NO ERROR DETECTED.             *
!*                 =-2  DID NOT CONVERGE IN ALLOWED NUMBER OF ITERAT-  *
!*                      IONS.                                          *
!*                 =-3  NO SOLUTION DETERMINED BECAUSE THE MATRIX      *
!*                      SENT TO SYMPDS WAS NOT SPD OR SECBI FAILED.    *
!*                                                                     *
!*         NKMAX   AN INPUT INTEGER WHICH AGREES EXACTLY WITH THE      *
!*                 FIRST DIMENSION OF ARRAY C. NKMAX MUST BE GREATER   *
!*                 THAN OR EQUAL TO 2*NK.                              *
!*                                                                     *
!*                                                                     *
!*                                                                     *
!*                                                                     *
!*                                                                     *
!*    REQUIRED ROUTINES                QXZ063,QXZ277,QXZ278,QXZ279,    *
!*                                     QXZ280,QXZ282,SYMPDS(F4.3),     *
!*                                     SECBI(C2.5)                     *
!*                                                                     *
!*    SOURCE                                                           *
!*                                     PROGRAM RASPFIT BY JAMES R.     *
!*                                     SCHIESS AND PATRICIA A. KERR    *
!*                                     MODIFIED BY COMPUTER SCIENCES   *
!*                                     CORPORATION                     *
!*                                                                     *
!*    LANGUAGE                         -FORTRAN                        *
!*                                                                     *
!*    DATE RELEASED                    FEBRUARY, 1986                  *
!*                                                                     *
!*    LATEST REVISION                  FEBRUARY, 1986                  *
!***********************************************************************
!
!
!
!     PARAMETERS USED
!
!
      DIMENSION X(N),Y(N),W(N),DEV(LM1),IPK(LM1),PK(LM1)
      DIMENSION XK(L),R(NN1),XR(NN1),DER1(NN1),DER2(NN1)
      DIMENSION C(NKMAX,L2),ANS(L2),E(L2,LL),B(LM2,L2),GG(LM2,LM2)
!
!     THE FIRST DIMENSION OF C WAS CHANGED FROM L2 TO NKMAX
!
      DIMENSION D(L2,L2),G(LM2,LM2),XLAM(LM2)
      DIMENSION TEMP3(L2,L2)
!
!!!      INTEGER QXZ063
      IF(IW.EQ.1)GO TO 70
!
!
      CALL QXZ278(X,N,XK,L,ISTOP)
      IERR=-ISTOP
      IF(IERR .NE. 0)  GO TO 120
      IPKZERO=LM1
!
      DO I=1,LM1
        IPKZERO=IPKZERO-IPK(I)
      END DO
!
      DO 50 J=1,ITMAX
      CALL QXZ277(X,Y,W,N,XK,L,PK,LM1,ANS,L2,E,GG,B,LM2,LL,D,XLAM,G,IERR)
      IF(IERR.NE.0)GO TO 120
      CALL QXZ279(IPK,XK,PK,DEV,LM1,ANS,L,L2,IPKZERO,IERR)
      IF(IERR.NE.0)GO TO 120
      IF(IPKZERO .EQ. 0)GO TO 70
   50 END DO
!
!              CALCULATE STANDARD ERROR
!
      IERR=-2
   70 STDERR=0.
      FLGRP=0.
      FLGRPP=0.
      K = 1
      IF(NN.GT.0)GO TO 960
      I=QXZ063(XR(1),XK,L)
      K=I
      DXK = XK(K+1) - XK(K)
      CALL QXZ282(DXK,PK(K),ANS(2*K-1),ANS(2*K),ANS(2*K+1),ANS(2*K+2),  &
     &     XK(K),XR(1),R(1),DER1(1),DER2(1))
!
      IF(XR(1).GT.XR(NN1))GO TO 900

      DO 940 J=2,NN1
   40 IF(XR(J).GT.XK(I+1))GO TO 200
      K=I
      DXK = XK(K+1) - XK(K)
      CALL QXZ282(DXK,PK(K),ANS(2*K-1),ANS(2*K),ANS(2*K+1),ANS(2*K+2),  &
     &     XK(K),XR(J),R(J),DER1(J),DER2(J))
      CYCLE   !!! GO TO 940
  200 I=I+1
      GO TO 40
  940 END DO
      GO TO 980
  900 DO 920 J=2,NN1
  930 IF(XR(J).LT.XK(I))GO TO 910
      K=I
      DXK = XK(K+1) - XK(K)
      CALL QXZ282(DXK,PK(K),ANS(2*K-1),ANS(2*K),ANS(2*K+1),ANS(2*K+2),  &
     &     XK(K),XR(J),R(J),DER1(J),DER2(J))
      CYCLE   !!! GO TO 920
  910 I=I-1
      GO TO 930
  920 END DO
      GO TO 980
  960 DO 280 I=1,NN
      K=QXZ063(XR(I),XK,L)
      DXK = XK(K+1) - XK(K)
      CALL QXZ282(DXK,PK(K),ANS(2*K-1),ANS(2*K),ANS(2*K+1),ANS(2*K+2),  &
     &     XK(K),XR(I),R(I),DER1(I),DER2(I))
  280  CONTINUE
!
  980 IF(IW.EQ.1)GO TO 120
      IW=1
      K=1
      DO 80 I=1,N
      IF (X(I).GT.XK(K+1)) K=K+1
      DXK = XK(K+1) - XK(K)
      CALL QXZ282(DXK,PK(K),ANS(2*K-1),ANS(2*K),ANS(2*K+1),ANS(2*K+2),  &
     &     XK(K),X(I),F,FLGRP,FLGRPP)
      STDERR= STDERR +(Y(I)-F)**2
   80 END DO
      STDERR= SQRT(STDERR/(N-2*L))
!
!
!      CALCULATE COVARIANCE MATRIX
!
      DO 230 KK=1,L2
      DO 231 II=1,L2
      TEMP3(II,KK) = 0.0
      DO 232 JJ=1,L2
      DO 233  J=1,LM2
      DO 234  M=1,LM2
        TEMP3(II,KK)=TEMP3(II,KK)-B(M,II)*GG(M,J)*B(J,JJ)*E(JJ,KK+1)
  234 END DO
  233 END DO
  232 END DO
  231 END DO
  230 END DO

      DO 260 III=1,L2
      DO 261 KKK=1,L2
      C(III,KKK) = 0.0
      IF(III .EQ. KKK)  C(III,KKK) = C(III,KKK)+E(III,III+1)
      DO 262 JJJ=1,L2
        C(III,KKK) = C(III,KKK) + E(III,JJJ+1)*TEMP3(JJJ,KKK)
  262 END DO
  261 END DO
  260 END DO
!
!     CALCULATE CORRELATION MATRIX WITH STANDARD DEVIATIONS ON DIAGONAL
!
      DO 100 I=1,L2
      TEMP = SQRT(STDERR*C(I,I))
      DO 110 J=1,L2
      C(I,J) = STDERR*C(I,J)/TEMP
      C(J,I) = STDERR*C(J,I)/TEMP
  110 END DO
      C(I,I) = TEMP
  100 END DO

  120 RETURN
      END Subroutine QXZ281

!+
SUBROUTINE QXZ282(dx,pk,ybar,ybarpp,ybar1,ybar1pp,xk,x,f,fp,fpp)
! ------------------------------------------------------------------------------
! once called FRATS
! PURPOSE - Calculates the rational spline approximations of f,f',f''
IMPLICIT NONE
  REAL,INTENT(IN):: dx  ! length of subinterval defined by two consecutive
                        ! knot abscissas (xk(k+1) - x(k))

  REAL,INTENT(IN):: pk      ! tension factor
  REAL,INTENT(IN):: ybar    ! y-coordinate of knot xk(j)
  REAL,INTENT(IN):: ybarpp  ! second derivative of knot xk(j)
  REAL,INTENT(IN):: ybar1   ! y-coordinate of knot xk(j+1)
  REAL,INTENT(IN):: ybar1pp ! second derivative at knot xk(j+1)
  REAL,INTENT(IN):: xk      ! x-coordinate of knot xk(j)
  REAL,INTENT(IN):: x       ! x-coordinate of data point
                   
  REAL,INTENT(OUT):: f      ! approximation of function value at x
  REAL,INTENT(OUT):: fp     ! approximation of first derivative
  REAL,INTENT(OUT):: fpp    ! approximation of second derivative

  REAL:: t,u,pt,pu,hk
!-------------------------------------------------------------------------------
  t=(x-xk)/dx
  u=1.-t
  pt=pk*t+1.
  pu=pk*u+1.
  hk = dx**2/(2.*(pk**2+3.*pk+3.))

  f=u*ybar + hk*(u**3/pt-u)*ybarpp + t*ybar1 +hk*(t**3/pu-t)*ybar1pp

  fp=(-ybar-hk*(((3.*u**2*pt+u**3*pk)/pt**2)-1.)*ybarpp+ybar1+hk*   &
       (((3.*t**2*pu+t**3*pk)/pu**2)-1.)*ybar1pp)/dx

  fpp=hk*ybarpp*(2.*pk**2*u**3+6.*pk*pt*u**2+6.*pt**2*u)/pt**3      &
        +hk*ybar1pp*(2.*pk**2*t**3+6.*pk*pu*t**2+6.*pu**2*t)/pu**3
  fpp=fpp/dx**2

  RETURN
END Subroutine QXZ282   ! ------------------------------------------------------


END Module RationalSpline   ! ==================================================
