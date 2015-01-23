INCLUDE 'ga.f90'
!+
PROGRAM Problem3
! ---------------------------------------------------------------------------
! PURPOSE - Solve Problem 3 in NASA TN D-7808, p. 159

! Determine a temperature-entropy (T-S) diagram for fluorine from 55 to
!    200 K with pressure as a parameter. Also determine the isoquality
!    lines from 0 to 1 in increments of one tenth. Let the pressures be
!    as given in the following table:
!    0.1, 1.0, 2.0, 3.0, 4.0, 5.215, 6.0, 7.0, 8.0, 10.0

! The problem can be solved by incrementing temperature (KS=1) or
!    entropy (KS=5).  However, it is faster to increment temperatures.
! Therefore, select the temperature increment to be 5 K.
USE GaspPdas
IMPLICIT NONE


!      COMMON /PROPTY/ ku,dl,dv,hl,hv,s,sl,sv,cv,cvl,cvv,cp,cpl,cpv,     &
!     & gamma,gammal,gammav,c,ll,cvp, mu,mul,muv,k, kl, kv, sigma,        &
!     & excesk,excl,excv
!       REAL ML,MUL,MUV,K,KL,KV

  REAL:: d,h
  INTEGER:: i,j,ij, kp,kr,ks, m
  INTEGER:: n = 14
  INTEGER:: nn

        REAL p(10), x(2000),y(2000)
!!!        DIMENSION TITL(7)
!!        DATA TL/ 42H TEMPERATURE ENTROPY DIAGRAM FOR FLUORINE/
        DATA p/0.1, 1.0, 2.0, 3.0, 4.0, 5.215, 6.0, 7.0, 8.0, 10.0/
!!!        DATA TSTART,PC,N/ 55., 5.215, 40/
!!!        DATA tstart,pc,n/80.,5.215,14/
  REAL:: pc = 5.215
  REAL:: t
  REAL:: tsat
  REAL:: tstart = 80.0
        CHARACTER(LEN=3),PARAMETER:: NAM = "F2 "
  REAL:: xx
!----------------------------------------------------------------------------
      CALL Setup(nam)
      WRITE(6,*) "TEMPERATURE ENTROPY DIAGRAM FOR FLUORINE"
      ku=1
      kp=2
      ks=1
      m=0
      DO 3 i=1,10 
        t=tstart
        IF (p(i) .GT. pc) GO TO 5
        kr=1
        tsat=0.

        CALL GasProperties(ks,kp,tsat,p(i),d,h,kr)

        x(m+1)=sl
        y(m+1)=tsat
        x(m+11)=sv
        y(m+11)=tsat
        DO 2 ij=1,9
          nn=m+ij + i
          xx=0.1*REAL(IJ)
          x(nn)=xx*sv+(1.-xx)*sl
          y(nn)=tsat
    2   END DO
        m=m+11
    5   DO 33 j=1,n
          IF (t.EQ.tsat) CYCLE
          kr=0
          CALL GasProperties(ks,kp,t,p(i),d,h,kr)
          m=m+1
          x(m)=s
          y(m)=t
          write(*,'(2F8.3)') x(m),y(m)
          t=t+5.0
   33   END DO
    3 END DO
!!!      WRITE(6,1000) (X(J),Y(J), J=1,M)
    1000 FORMAT(5(F6.3,F8.3))
!!!      CALL  LRPLOT (X,Y,M)
  WRITE(*,*) "End of Problem 3"
  STOP 
END Program Problem3   ! ==============================================

