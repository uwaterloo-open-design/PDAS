INCLUDE 'fluid.f90'
INCLUDE 'ar.f90'
INCLUDE 'ch4.f90'
INCLUDE 'co2.f90'
INCLUDE 'dryair.f90'
INCLUDE 'f2.f90'
INCLUDE 'n2.f90'
INCLUDE 'o2.f90'
INCLUDE 'ph2.f90'
INCLUDE 'steam.f90'

!+
PROGRAM SampleCase
! ---------------------------------------------------------------------------
! PURPOSE - Test for the proper operation of the FLUID subroutine
!    package and enable the user to compute thermodynamic properties.
!    User can determine properties of any of the gases in Fluid
!    This version uses US units for input and output, although Fluid
!    does all its calculations in SI units

! AUTHOR - Ralph L. Carmichael, Public Domain Aeronautical Software
!   (based on the sample program supplied with COSMIC program LEW-14418)

! REVISION HISTORY
!   DATE  VERS PERSON  STATEMENT OF CHANGES
! 15Dec97  0.1   RLC   Composed alternate to supplied sample case
!  2Jan98  0.5   RLC   Added units to dialog

  IMPLICIT NONE

  INTEGER,PARAMETER:: NP=8

  REAL:: MU,K,MUL,MUV,KL,KV 
  REAL G(29,4),TT(27),DD(30),GG(27,30,8),SAT(40,8),C(8)
  INTEGER:: M(4)
  INTEGER:: ENTRY,ERROR
  INTEGER:: kprop
  LOGICAL:: vapor
  REAL,DIMENSION(8):: props

  REAL:: T,P,RHO,z,s,h,cv,cp,sonic
  REAL:: ZL,ZV,SL,SV,HL,HV,RHOL,RHOV,CVL,CVV
  REAL:: CPL,CPV,CL,GAMMAL,GAMMAV


      
  INTEGER:: NG,N1,N2,N3,NS
  COMMON /FLUDPC/ G,TT,DD,GG,SAT,C,M,NG,N1,N2,N3,NS

  REAL:: GAMMA,WL,WG,DENSL,DENSG,ENTL,ENTG,ENTHL,ENTHG 
  COMMON /FLUIDC/ GAMMA,WL,WG,DENSL,DENSG,ENTL,ENTG,ENTHL,ENTHG 
!                                                                       
!            KPROP = FLUID OPTION                                       
!              1   = HYDROGEN                                           
!              2   = OXYGEN                                             
!              3   = NITROGEN                                           
!              4   = ARGON                                              
!              5   = FLUORINE                                           
!              6   = STEAM                                              
!              7   = DRYAIR                                             
!              8   = CO2                                                
!              9   = CH4                                                
!     UNITS T = RANKINE (R)                                             
!           P = PSIA                                                    
!           D = LBM/CU FT                                               
!           H = BTU/LBM                                                 
!           S = BTU/LBM-R                                               
!           CV = BTU/LBM-R                                              
!           CP = BTU/LBM-R                                              
!           SONIC = FT/SEC                                              
!           MU = LBM/FT-SEC                                             
!           K = BTU/FT-SEC-R                                            
!
  INTERFACE
    SUBROUTINE Fluid(T,P,RHO,PROPS,NP,ENTRY,VAPOR,ERROR,igoto)
      REAL,INTENT(IN OUT):: t
      REAL,INTENT(IN OUT):: p
      REAL,INTENT(IN OUT):: rho
      REAL,INTENT(OUT),DIMENSION(:):: props
      INTEGER,INTENT(IN):: np
      INTEGER,INTENT(IN):: entry
      LOGICAL,INTENT(OUT):: vapor
      INTEGER,INTENT(OUT):: error
      INTEGER,INTENT(IN):: igoto
    END Subroutine Fluid
  END INTERFACE
!----------------------------------------------------------------------------
  DO
    WRITE(*,*) "Select entry:"
    WRITE(*,*) "  1 to enter temperature and density   "
    WRITE(*,*) "  2 to enter pressure and density   "
    WRITE(*,*) "  3 to enter temperature and pressure   "
    WRITE(*,*) "  4 to enter pressure and entropy   "
    WRITE(*,*) "  5 to enter pressure and enthalpy    "
    READ(*,*) entry
    IF (entry > 0 .AND. entry < 6) EXIT
  END DO

  SELECT CASE (entry)
    CASE (1)
      WRITE(*,*) "Enter temperature in degrees Rankine   "
      READ(*,*) t
      t=t/1.8
      WRITE(*,*) "Enter density in pounds per cubic foot   "
      READ(*,*) rho
      rho=rho/62.428
    CASE (2)
      WRITE(*,*) "Enter pressure in pounds per square inch   "
      READ(*,*) p
      p=p/145.04
      WRITE(*,*) "Enter density in pounds per cubic foot   "
      READ(*,*) rho
      rho=rho/62.428
    CASE (3)
      WRITE(*,*) "Enter temperature in degrees Rankine   "
      READ(*,*) t
      t=t/1.8
      WRITE(*,*) "Enter pressure in pounds per cubic foot   "
      READ(*,*) p
      p=p/145.04
    CASE (4)
      WRITE(*,*) "Enter pressure in pounds per square inch   "
      READ(*,*) p
      p=p/145.04
      WRITE(*,*) "Enter entropy in BTU per pound per degree Rankine   "
      READ(*,*) s
      props(2)=s/0.23901
    CASE (5)
      WRITE(*,*) "Enter pressure in pounds per cubic foot   "
      READ(*,*) p
      p=p/145.04
      WRITE(*,*) "Enter enthalpy in BTU per pound   "
      READ(*,*) h
      props(3)=h/0.43022
  END SELECT

  DO
    WRITE(*,*) "Select your fluid:"
    WRITE(*,*) "  1 is for Hydrogen  "
    WRITE(*,*) "  2 is for Oxygen  "
    WRITE(*,*) "  3 is for Nitrogen  "
    WRITE(*,*) "  4 is for Argon  "
    WRITE(*,*) "  5 is for Fluorine  "
    WRITE(*,*) "  6 is for Steam  "
    WRITE(*,*) "  7 is for Dry Air  "
    WRITE(*,*) "  8 is for Carbon Dioxide  "
    WRITE(*,*) "  9 is for Methane  "

    READ(*,*) kprop
    IF (kprop > 0 .AND. kprop < 10) EXIT
  END DO

!       WRITE (*,*) 'ENTER:P(PSIA),T(R),S(BTU/LBM-R),H(BTU/LBM),RHO(LBM/C&
!     & FT),KPROP(1-9),ENTRY(1-5)'                                       
!      READ (*,*)P,T,S,H,RHO,KPROP,ENTRY 
!       WRITE(6,'(T35,A)')' INPUT DATA:' 
!       WRITE(6,*)' P=',P,' T  =',T,  ' S    =',S 
!      WRITE(6,*)' H=',H,' RHO=',RHO,' ENTRY=',ENTRY 
!      INITIALIZE FLUID VARIABLES FOR THE DESIRED FLUID.                
  SELECT CASE(kprop)
    CASE (1)
      WRITE(*,*) 'HYDROGEN' 
      CALL PH2
      WRITE(*,*) "Setup for hydrogen completed"
    CASE (2)
      WRITE(*,*)'OXYGEN' 
      CALL O2 
      WRITE(*,*) "Setup for oxygen completed"
    CASE (3)
      WRITE(*,*)'NITROGEN' 
      CALL NITRO 
      WRITE(*,*) "Setup for nitrogen completed"
    CASE (4)
      WRITE(*,*)'ARGON' 
      CALL AR 
      WRITE(*,*) "Setup for argon completed"
    CASE (5)
      WRITE(*,*)'FLUORINE' 
      CALL F2 
      WRITE(*,*) "Setup for fluorine completed"
    CASE (6)
      WRITE(*,*)'STEAM' 
      CALL STEAM 
      WRITE(*,*) "Setup for steam completed"
    CASE (7)
      WRITE(*,*)'DRYAIR' 
      CALL DRYAIR 
      WRITE(*,*) "Setup for dry air completed"
    CASE (8)
      WRITE(*,*)'CARBON DIOXIDE' 
      CALL CO2 
      WRITE(*,*) "Setup for carbon dioxide completed"
    CASE (9)
      WRITE(*,*)'METHANE' 
      CALL CH4 
      WRITE(*,*) "Setup for methane completed"
    END SELECT


  WRITE(*,*) "OUTPUT DATA" 
!      CONVERSIONS TO INPUT DATA FOR FLUID                              
!      CONVERT PRESSURE IN PSI TO MPA AND TEMPERATURE IN R TO K         
!      CONVERT DENSITY IN LB/FT3 TO G/CM3                               
!      CONVERT ENTROPY IN BTU/LB-R TO J/GM-K                            
!      CONVERT ENTHALPY IN BTU/LBM TO J/GM                              
!***      P = P/145.04          ! convert psia to MPa
!***      T = T/1.8             ! convert deg R to deg K
!***      RHO = RHO/62.428      ! convert lb/cu ft to g/cu cm
!***      PROPS(2) = S/.23901   ! convert BTU/lb-R to joules/gm-K
!***      PROPS(3) = H/.43022   ! convert BTU/lbm to J/gm

                                             
  CALL Fluid(T,P,RHO,PROPS,NP,ENTRY,VAPOR,ERROR,1) 
  IF (error==0) THEN
    WRITE(*,*) "Normal return from Fluid"
  ELSE
    WRITE(*,*) "Error return from fluid, error= ", error
  END IF

!***  T = T*1.8       ???
!***  P = P*145.04    ???

  IF (VAPOR) THEN          ! see if fluid is saturated
    SL = ENTL*.23901 
    SV = ENTG*.23901 
    HL = ENTHL*.43022 
    HV = ENTHG*.43022 
    RHOL = DENSL*62.428 
    RHOV = DENSG*62.428
!         GET THE TRANSPORT PROPERTIES FOR SATURATED LIQUID             
    CALL FLUID(T,P,RHO,PROPS,NP,ENTRY,VAPOR,ERROR,2) 
    ZL=PROPS(1) 
    CVL = PROPS(4)*.23901 
    CPL = PROPS(5)*.23901 
    CL = PROPS(6)*103.83 
    MUL = PROPS(7)*.067198 
    KL = PROPS(8)*.016061 
    GAMMAL = CPL/CVL 
!         GET THE TRANSPORT PROPERTIES FOR SATURATED GAS                
    CALL FLUID(T,P,RHO,PROPS,8,ENTRY,VAPOR,ERROR,3) 
    ZV=PROPS(1) 
    CVV = PROPS(4)*0.23901 
    CPV = PROPS(5)*0.23901 
    CV = PROPS(6)*103.83 
    MUV = PROPS(7)*0.067198 
    KV = PROPS(8)*0.016061 
    GAMMAV = CPV/CVV 

    p=p*145.04
    t=t*1.8
    WRITE(*,*) "Pressure=", p, " psia"
    WRITE(*,*) "Temperature=", t, " deg. R"
    WRITE(*,*) "                  SATURATED LIQUID    SATURATED VAPOR"
    WRITE(*,*) "Compressibility  ", zl,zv
    WRITE(*,*) "Entropy", sl,sv, " BTU/lb/degR"
    WRITE(*,*) "Enthalpy", hl,hv, " BTU/lb"
    WRITE(*,*) "Density", rhol,rhov, " lb/cubic foot"
    WRITE(*,*) "Cv", cvl,cvv, " BTU/lb-degR"
    WRITE(*,*) "Cp", cpl,cpv, " BTU/lb-degR"
    WRITE(*,*) "Sonic Velocity", cl,cv, " ft/sec"
    WRITE(*,*) "Cp/Cv", gammal,gammav
    WRITE(*,*) "Viscosity", mul,muv, "  lbm/ft-s"
    WRITE(*,*) "Thermal conductivity", kl,kv, "  BTU/ft-sec-degR"

!    WRITE(6,1000) P,T,ZL,ZV,SL,SV,HL,HV,RHOL,RHOV,CVL,CVV,         &
!     &              CPL,CPV,CL,CV,GAMMAL,GAMMAV,MUL,MUV,KL,KV           
! 1000  FORMAT (' P=',G10.3,' PSI',T40,'T=',G10.3,' R'/                  &
!     &T20,'SATURATED LIQUID',T50,'SATURATED VAPOR'/                     &
!     &  ' COMPRESSIBILITY',T25,G10.3,T55,G10.3/                         &
!     &  ' ENTROPY(BTU/LBM-R)',T25,G10.3,T55,G10.3/                      &
!     &  ' ENTHALPY(BTU/LBM)',T25,G10.3,T55,G10.3/                       &
!     &  ' DENSITY(LBM/FT3)',T25,G10.3,T55,G10.3/                        &
!     &  ' CV(BTU/LBM-R)',T25,G10.3,T55,G10.3/                           &
!     &  ' CP(BTU/LBM-R)',T25,G10.3,T55,G10.3/                           &
!     &  ' SONIC VELOCITY(FT/S)',T25,G10.3,T55,G10.3/                    &
!     &  ' CP/CV',T25,G10.3,T55,G10.3/                                   &
!     &  ' VISCOSITY(LBM/FT-S)',T25,G10.3,T55,G10.3/                     &
!     &  ' THERM. COND.(BTU/FT-S-R)',T25,G10.3,T55,G10.3)                
                                                                        
  ELSE 
    p=p*145.04
    t=t*1.8
    RHO = RHO*62.428
    Z=PROPS(1) 
    S = PROPS(2)*0.23901 
    H = PROPS(3)*0.43022 
    CV = PROPS(4)*0.23901 
    CP = PROPS(5)*0.23901 
    SONIC = PROPS(6)*103.83 
    MU = PROPS(7)*0.067198 
    K = PROPS(8)*0.016061 

    WRITE(*,*) "P=", p, "   T=",t
    WRITE(*,*) "Compressibility  ", z
    WRITE(*,*) "Entropy", s
    WRITE(*,*) "Enthalpy", h
    WRITE(*,*) "Density", rho
    WRITE(*,*) "Cv", cv
    WRITE(*,*) "Cp", cp
    WRITE(*,*) "Sonic Velocity", sonic
    WRITE(*,*) "Cp/Cv", gamma
    WRITE(*,*) "Viscosity", mu
    WRITE(*,*) "Thermal conductivity", k
  END IF

!     WRITE(6,1100) P,T,Z,S,H,RHO,CV,CP,GAMMA,SONIC,MU,K
!      ENDIF 
!      STOP 
! 1100  FORMAT (' P=',G10.3,' PSI',T40,'T=',G10.3,' R'/                  &
!     &  ' COMPRESSIBILITY',T25,G10.3/                                   &
!     &  ' ENTROPY(BTU/LBM-R)',T25,G10.3/                                &
!     &  ' ENTHALPY(BTU/LBM)',T25,G10.3/                                 &
!     &  ' DENSITY(LBM/FT3)',T25,G10.3/                                  &
!     &  ' CV(BTU/LBM-R)',T25,G10.3/                                     &
!     &  ' CP(BTU/LBM-R)',T25,G10.3/                                     &
!     &  ' CP/CV',T25,G10.3/                                             &
!     &  ' SONIC VELOCITY(FT/S)',T25,G10.3/                              &
!     &  ' VISCOSITY(LBM/FT-S)',T25,G10.3/                               &
!     &  ' THERM.COND.(BTU/FT-S-R)',T25,G10.3)                           

  STOP "End of Sample Case"
END Program SampleCase   ! ==================================================

!  INCLUDE 'fluid.f90'
!  INCLUDE 'ar.f90'
!  INCLUDE 'ch4.f90'
!  INCLUDE 'co2.f90'
!  INCLUDE 'f2.f90'
!  INCLUDE 'ph2.f90'
!  INCLUDE 'dryair.f90'
!  INCLUDE 'n2.f90'
!  INCLUDE 'o2.f90'
!  INCLUDE 'steam.f90'
!!!  INCLUDE 'ntrp.f90'
!!!  INCLUDE 'cubic.f90'

