      INCLUDE 'g.for'
      PROGRAM SaturationProperties
! ---------------------------------------------------------------------------
! PURPOSE - Print a table of saturation properties for each gas

      IMPLICIT NONE
      INTEGER:: ku
      REAL:: dl,dv,hl,hv,s,sl,sv,cv,cvl,cvv,cp,cpl,cpv   
      REAL:: gamma,gammal,gammav,c,cl,cvp,mu,mul,muv,k,kl,kv,sigma
      REAL:: excl,excv,excesk  
      COMMON/PROPTY/ku,dl,dv,hl,hv,s,sl,sv,cv,cvl,cvv,cp,cpl,cpv   
     1,gamma,gammal,gammav,c,cl,cvp,mu,mul,muv,k,kl,kv,sigma,excl,  
     2 excv,excesk                      
      REAL:: dch1,dch2,pch1,pch2,pch3,tch1,tch2,tch3
      REAL:: dst,tst,hsch1,hsch2                      
     
      COMMON/CHECKS/dch1,dch2,pch1,pch2,pch3,tch1,tch2,tch3,dst,tst,     
     1 hsch1,hsch2                      

      REAL:: d,h
      INTEGER:: i,j
      INTEGER:: kp,kr,ks 
      INTEGER,PARAMETER:: NGAS=10

      REAL:: t
      REAL,DIMENSION(NGAS):: tstart

      CHARACTER(LEN=3),PARAMETER,DIMENSION(NGAS):: nam =                &
     & (/'CH4', 'N2 ','O2 ','AR ','CO2','NE ', 'CO ','HE ','F2 ','H2 '/)
      CHARACTER(LEN=*),PARAMETER:: STRING1 =                            &
     &"  T     density  entropy  enthalpy    Cp     Cv       c    gamma"
      CHARACTER(LEN=*),PARAMETER:: STRING2 =                            &
     &  " degK   gm/cc     J/G-K     J/gm    J/gm-K  J/gm-K  cm/sec"



      DATA TSTART/92.,65.,70.,84.,217.,25.,70.,3.,70.,15./  
      REAL:: z
!----------------------------------------------------------------------------
      WRITE(*,*) "Starting sample3"
      OPEN(UNIT=6,FILE='sample3.out',STATUS='REPLACE',ACTION='WRITE')
      DO 1000 I=1,10
        WRITE(6,169)      
        CALL SETUP(NAM(I))
!        WRITE(6,101)
        WRITE(6,*) "SATURATION PROPERTIES"
        WRITE(6,*) string1
        WRITE(6,*) string2
C                               
C    TEST SATURATION PROPERTIES 
C                   
      KU=1          
      T=TSTART(I)   

      DO 100 J=1,50  
      Z=0.0          
      KR=1           
      KS=1           
      KP=7           
      CALL Gasp(ks,kp,t,z,d,h,kr) 
!!!      WRITE(6,168)T,Z             
      WRITE(6,170) t,z,DL,HL,SL,CVL,CPL,GAMMAL,CL
      WRITE(6,170) t,z,DV,HV,SV,CVV,CPV,GAMMAV,CVP                         

      T=T+5.0
      IF (T.GT.TCH2) EXIT
  100 END DO   

 1000 END DO   

  169 FORMAT   ('1')
  101 FORMAT (' SATURATION PROPERTIES'/ ' DENSITY-G/CC     ',
     1 ' ENTHALPY-J/GM        ENTROPY-J/GM-K        CV-J/GM-K',
     2 '    CP-J/GM-K         GAMMA          SONIC VEL-CM/SEC'    )                                                               36
  168 FORMAT ( ' T =', F9.3,'K     PSAT =',  F10.4, 'MN/M2')
!  170 FORMAT('0',3(E16.8,4X),3(E14.6,4X),E14.6 /
!     1        1X,3(E16.8,4X),3(E14.6,4X),E14.6 )
  170 FORMAT(F6.1,8ES11.3)
      WRITE(*,*) "Sample 3"
      STOP 
      END Program SaturationProperties   ! =============================

