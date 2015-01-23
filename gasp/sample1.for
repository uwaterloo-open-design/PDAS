      INCLUDE 'g.for'
!+
       PROGRAM ShortIsobar
! ---------------------------------------------------------------------------
! PURPOSE - Prepare a table of density, entropy, enthalpy, Cp, Cv, gamma,
!   and sonic velocity for each gas over a range of PVT.
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

      INTEGER,PARAMETER:: NGAS = 10   ! the number of gases 

!!!      REAL,DIMENSION(NGAS):: tstart,tchg,dt1,dt2
      REAL,PARAMETER,DIMENSION(3,NGAS):: TSTART = RESHAPE(              &
     & (/ 95.0,295.0,400.0, 65.0,200.0,400.0, 60.0,225.0,400.0,         &
     &    85.0,225.0,400.0, 220.0,420.0,500.0, 25.0,80.0,400.0,         &
     &    70.0,200.0,400.0, 3.0,10.0,100.0, 70.0,200.0,400.0,           &
     &    16.0,40.0,200.0 /), (/ 3, NGAS /) )
      REAL,PARAMETER,DIMENSION(3,NGAS):: DTA = RESHAPE(                 &
     &  (/ 10.0,20.0,100.0, 10.0,20.0,100.0, 10.0,20.0,100.0,           &
     &     10.0,20.0,100.0, 10.0,20.0,100.0, 10.0,20.0,100.0,           &
     &     10.0,20.0,100.0,  0.5,10.0,100.0,  5.0,20.0,100.0,           &
     &      2.0,10.0,100.0 /), (/ 3,NGAS /) )
      REAL p(5),t(250),ypl(250,5,7)
      CHARACTER*3 nam(10)
      CHARACTER(LEN=*),PARAMETER:: STRING1 =                            &
     &"  T     density  entropy  enthalpy    Cp     Cv       c    gamma"
      CHARACTER(LEN=*),PARAMETER:: STRING2 =                            &
     &  "   K    gm/cc     J/G-K     J/gm    J/gm-K  J/gm-K  cm/sec"

      REAL:: d,dt
      REAL:: h
      INTEGER:: i,j
      INTEGER:: kgas
      INTEGER:: kp,kpt,kr,ks
      INTEGER:: n
!!      DATA DT1/7*10.,.2,5.,2./                                                
!!      DATA DT2/10*50./                                                        
!!      DATA TCHG/295.,200.,2*225.,420.,80.,200.,10.,200.,80./                  
      DATA NAM/ 'CH4', 'N2 ', 'O2 ', 'AR ', 'CO2',                      &
     &          'NE ', 'CO ', 'HE ', 'F2 ', 'H2 '/
!!      DATA TSTART/95.,65.,60.,85.,220.,25.,70.,3.,70.,16./       
      REAL:: tj,ts
      REAL:: z    
!----------------------------------------------------------------------------
      WRITE(*,*) "Starting Sample1"
      OPEN(UNIT=6, FILE='sample1.out',STATUS='REPLACE',ACTION='WRITE')
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
          DO 50 J=1,KPT                                                          
            KR=0                                                                    
            KP=7                                                                    
            CALL Gasp(KS,KP,T(J),Z,D,H,KR)                                         
   30       IF (KR.EQ.1) GO TO 35                                                    
            YPL(J,I,1)  = D                                                         
            YPL(J,I,2)  = S                                                         
            YPL(J,I,3)  = H                                                         
            YPL(J,I,4)  = CP                                                        
            YPL(J,I,5)  = CV                                                        
            YPL(J,I,6)  = C/100.                                                    
            YPL(J,I,7)  = CP/CV
            GO TO 45         
   35       YPL(J,I,1)  = DL 
            YPL(J,I,2)  = SL 
            YPL(J,I,3)  = HL 
            YPL(J,I,4)  = CPL     
            YPL(J,I,5)  = CVL     
            YPL(J,I,6)  = CL /100.
            YPL(J,I,7)  = CPL/CVL
   45       CONTINUE
!!!   45 ZSAV=Z           ! never used ???
   46 FORMAT(I7,6E14.6)
   50     CONTINUE

          WRITE(6,*) "THERMODYNAMIC PROPERTY TABLE --- " // nam(kGas)
          WRITE(6,801)  P(I)     
          WRITE(6,*) STRING1
          WRITE(6,*) STRING2
          WRITE(6,803) (t(j),(ypl(j,i,n),n=1,7),j=1,kpt)
  100   END DO                                      
 1000 END DO
      WRITE(*,*) "Normal termination of ShortIsobar"
      STOP 

   51 FORMAT('  ISOBARS  =',  5F10.3 )
!!  800 FORMAT(' THERMODYNAMIC PROPERTY TABLE --  ',A4 ) 
  801 FORMAT(F15.2, ' ATMOSPHERE ISOBAR')                  
  802 FORMAT(' T-K   DENSITY (GM/CC)  ENTROPY(J/G-K) ',
     1 'ENTHALPY(J/GM)  CP(J/GM-K)  CV(J/GM-K)   C( M/SEC)    GAMMA')
  803 FORMAT(F6.1,ES12.4,F7.2,F9.2,F9.4,F8.4,F9.3,F6.3)      
      END Program ShortIsobar   ! ===========================================

