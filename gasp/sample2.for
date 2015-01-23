      INCLUDE 'g.for'
!+
       PROGRAM ShortTransport
! ---------------------------------------------------------------------------
! PURPOSE - Prepare a table of transport properties for each gas
!

!     THIS TEST PROGRAM FOR -GASP- IS A SHORTENED VERSION.......        
!     THE LIMITS OF THE TABLE ARE NOT NECESSARILY THE MAXIMUM RANGE OF G
!      T, P, OR D.   THE RANGE IS REPRESENTATIVE AND IS MEANT TO BE OF
!     HELP IN CHECKING OUT THE PROGRAM ON VARIOUS COMPUTERS.                    
!
! from Appendix I in TN D-7808
      IMPLICIT NONE
      INTEGER:: ku
      REAL:: dl,dv,hl,hv,s,sl,sv,cv,cvl,cvv,cp,cpl,cpv
      REAL:: gamma,gammal,gammav,c,cl,cvp,mu,mul,muv,k,kl,kv,sigma
      REAL:: excl,excv,excesk    
      COMMON /PROPTY/ ku,dl,dv,hl,hv,s,sl,sv,cv,cvl,cvv,cp,cpl,cpv,
     1 gamma,gammal,gammav,c,cl,cvp,mu,mul,muv,k,kl,kv,sigma,excl,   
     2 excv,excesk                                      

      REAL:: dch1,dch2,pch1,pch2,pch3,tch1,tch2,tch3
      REAL:: dst,tst,hsch1,hsch2 
      COMMON/CHECKS/dch1,dch2,pch1,pch2,pch3,tch1,tch2,tch3,dst,tst,
     1 hsch1,hsch2  
      REAL:: pdt,ptv
      COMMON/PARTLS/pdt,ptv    

      REAL:: d,dt
      REAL:: h
      INTEGER:: i,j,kgas,kp,kpt,kr,ks,n
      INTEGER,PARAMETER:: NGAS = 10   ! the number of gases 

      REAL,PARAMETER,DIMENSION(3,NGAS):: tstart = RESHAPE(              &
     & (/ 95.0,295.0,400.0, 65.0,200.0,400.0, 60.0,225.0,400.0,         &
     &    85.0,225.0,400.0, 220.0,420.0,500.0, 25.0,80.0,400.0,         &
     &    70.0,200.0,400.0, 3.0,10.0,100.0, 70.0,200.0,400.0,           &
     &    16.0,40.0,200.0 /), (/ 3, NGAS /) )
      REAL,PARAMETER,DIMENSION(3,NGAS):: dta = RESHAPE(                 &
     &  (/ 10.0,20.0,100.0, 10.0,20.0,100.0, 10.0,20.0,100.0,           &
     &     10.0,20.0,100.0, 10.0,20.0,100.0, 10.0,20.0,100.0,           &
     &     10.0,20.0,100.0,  0.5,10.0,100.0,  5.0,20.0,100.0,           &
     &      2.0,10.0,100.0 /), (/ 3,NGAS /) )
      CHARACTER(LEN=3),PARAMETER,DIMENSION(10):: nam =                  &
     &   (/ 'CH4', 'N2 ', 'O2 ', 'AR ', 'CO2',                          &
     &      'NE ', 'CO ', 'HE ', 'F2 ', 'H2 ' /)
      CHARACTER(LEN=*),PARAMETER:: STRING1 =                            &
     &"  T     viscosity conductivity surf.tension"
      CHARACTER(LEN=*),PARAMETER:: STRING2 =                            &
     &  " degK   gm/cc     g/cm-s J/cm-s-K   dyne/cm"

      REAL,DIMENSION(3):: p
      REAL,DIMENSION(250):: t
      REAL,DIMENSION(250,3,3):: ypl

      REAL:: tj,ts
      REAL:: z
!----------------------------------------------------------------------------
      WRITE(*,*) "Starting Sample2"
      OPEN(UNIT=6, FILE='sample2.out',STATUS='REPLACE',ACTION='WRITE')
      KS=1
      KU=1                                                              

      DO 1000 kGas=1,NGAS    ! begin big loop for each gas
        CALL SETUP(NAM(kGas))

        p(1)=1.0
        p(2)=PCH2/.101325    ! critical pressure for the gas
        p(3)=100.0

        t(1)=TSTART(1,kGas)
        DT=DTA(1,kGas)              ! build the table t
        ts=TSTART(2,kGas)
        n=2
        KPT=0
        DO J=1,50  
          KPT=J
          t(j+1)=t(j)+dt
          tj=t(j+1)
          IF (tj >= TCH3) EXIT
          IF (kGas==10 .AND. tj >= 1999.5) EXIT
          IF (tj < ts-0.05) CYCLE
          t(j+1)=ts
          dt=dta(n,kGas)
          n=n+1
          ts=tch3
          IF (n < 4) ts=tstart(n,kGas)
        END DO
        kpt=kpt+1

        DO 100 I=1,3   ! three different pressure isobars
          Z=P(I)*0.101325
          DO 50  J=1,KPT
            KR=0
            KP=56
            CALL GASP(KS,KP,T(J),Z,D,H,KR)
            ypl(j,i,3)=0.0
            IF (I.GT.1.OR.T(J).GT.TCH2) GO TO 30
            YPL(j,i,3)=sigma
   30       IF(KR.EQ.1) GO TO 35
            YPL(J,I,1)= MU
            YPL(J,I,2)= K
            GO TO 50
   35       YPL(J,I,1) = MUL
            YPL(J,I,2)= KL
   50     CONTINUE
          WRITE(6,51) P(I)


          WRITE(6,*) "THERMODYNAMIC PROPERTY TABLE --- " // nam(kGas)
          WRITE(6,801)  P(I)     
          WRITE(6,*) STRING1
          WRITE(6,*) STRING2
          WRITE(6,803) (T(J),(YPL(J,I,N),N=1,3),J=1,KPT)
  100   CONTINUE                                      
 1000 CONTINUE 
      WRITE(*,*)  "Normal termination of ShortTransport"
      STOP 

   51 FORMAT('  ISOBARS  =',  5F10.3 )
!!  800 FORMAT(' THERMODYNAMIC PROPERTY TABLE --  ',A4 ) 
  801 FORMAT(F15.2, ' ATMOSPHERE ISOBAR')                  
  803 FORMAT(F6.1,3ES12.4)      
      END Program ShortTransport   ! ========================================

