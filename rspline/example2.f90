INCLUDE 'ratspl.f90'

!+
PROGRAM Example2
! ------------------------------------------------------------------------------
! PURPOSE - An example of the use of the RationalSpline Module

! AUTHORS - James R. Schiess and Patricia A. Kerr, NASA Langley
!         - Ralph L. Carmichael, Public Domain Aeronautical Software

!     REVISION HISTORY
!   DATE  VERS PERSON  STATEMENT OF CHANGES
!   1984   0.1 JRS&PAK Original coding of program RASPFIT at NASA Langley
!   1988   0.2 JRS&PAK Release of program LAR-13694 to COSMIC
!   1996   0.3   RLC   Acquisition of program LAR-13694 from COSMIC
! 30Dec01  0.4   RLC   Converted to Fortran 90 free-format
! 29Dec08  0.5   RLC   Made some coding more in line with Fortran 2003


USE RationalSpline

  REAL,DIMENSION(84):: x,df
  REAL,DIMENSION(8):: dev = (/.8E+2,.8E+2,.8E+2,.15E+2,.15E+2,.1E+2, .8E+2,.2E+2/)
  REAL,DIMENSION(18,18):: c
  REAL,DIMENSION(1247):: wk
  REAL,DIMENSION(7):: r, der1, der2 
  REAL,PARAMETER,DIMENSION(9):: XK = (/ &
    0.1E+01,.4E+1,.7E+1,.12E+02,.25E+02,0.49E+02,.68E+02,.74E+02,.84E+02/)
  REAL,DIMENSION(8):: PK = (/ &
    0.1E+2,.0E+0,.0E+0,0.E+0,0.E+0,.0E+00, 0.E+0,0.E+0/)
  INTEGER,PARAMETER,DIMENSION(8):: IPK = (/1,1,1,1,0,0,1,0 /)
       REAL,PARAMETER,DIMENSION(84):: Y = (/                              &
     & 13.0,     12.5,     12.0,     11.5,     11.0,                     &
     & 13.7375,  16.475,   19.2125,  21.95,    23.5875,                  &
     & 25.225,   26.8625,  28.5,     29.3375,  30.175,                   &
     & 31.0125,  31.85,    32.10,    32.35,    32.60,                    &
     & 32.85,    32.75625, 32.6625,  32.56875, 32.475,                   &
     & 32.81875, 33.1625,  33.50625, 33.85,    34.00625,                 &
     & 34.1625,  34.31875, 34.475,   34.725,   34.975,                   &
     & 35.225,   35.475,   35.725,   35.975,   36.225,                   &
     & 36.475,   36.725,   36.975,   37.225,   37.475,                   &
     & 37.35,    37.225,   37.1,     36.975,   37.1,                     &
     & 37.225,   37.35,    37.475,   37.225,   36.975,                   &
     & 36.725,   36.475,   36.3,     36.125,   35.95,                    &
     & 35.775,   35.45,    35.125,   34.8,     34.475,                   &
     & 34.225,   33.975,   33.725,   33.475,   29.16875,                 &
     & 24.8625,  20.55625, 16.25,    15.7625,  15.275,                   &
     & 14.7875,  14.3,     14.475,   14.650,   14.825,                   &
     & 15.0,     15.175,   15.35,    15.525 /)

  REAL,DIMENSION(7):: xr = (/ 3.0, 0.0, 0.0, 0.0, 0.0, 0.0, 53.0 /)

  INTEGER:: i
  INTEGER,PARAMETER:: ITMAX = 35
  INTEGER,PARAMETER:: N=84
  INTEGER,PARAMETER:: NK=9
  INTEGER,PARAMETER:: NKMAX = 18
!-------------------------------------------------------------------------------
  DO I=1,84
    X(I)=REAL(I)
  END DO
  df(:)=0.0

      NN=-7
      IW=1
!!!      NKMAX=18
!!!      ITMAX=35
      READ(7,110) (WK(I),I=1,1247)
      READ(7,110) (PK(I),I=1,8)
      CALL RSDSDR(N,NK,X,Y,DF,XK,IPK,PK,DEV,NN,IW,                      &
     &NKMAX,ITMAX,XR,R,DER1,DER2,C,WK,IERR)
      WRITE(6,*) " IERR=", ierr

  IF (IERR.EQ.0) THEN
    WRITE(6,'(/A)') "       XR  INTERPOLATED    FIRST      SECOND"
    WRITE(6,'(A/)') "             VALUES      DERIVATIVE DERIVATIVE"

    DO i=1,7
      WRITE(6,'(F10.2,3F11.3)') XR(I),R(I),DER1(I),DER2(I)
    END DO
  END IF
  STOP
!!!   30 FORMAT(/5X,27H  XR  INTERPOLATED    FIRST,                        &
!!!     &7X,6HSECOND,/10X,26H   VALUES     DERIVATIVE  ,                   &
!!!     &10HDERIVATIVE/)
!!!   40 FORMAT(6X,F5.2,3X,F7.3,4X,F7.3,5X,F7.3)
!!!   90 FORMAT(12H    IERR IS ,I3)
  110 FORMAT(3E24.14)
      END
