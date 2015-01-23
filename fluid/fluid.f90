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
  INTEGER,INTENT(IN):: entry
  LOGICAL,INTENT(OUT):: vapor
  INTEGER,INTENT(OUT):: error
  INTEGER,INTENT(IN):: igoto  ! integer that specifies which variables are input:       
                              !  = 1 if temperature and density are given.                 
                              !  = 2 if pressure and density are given.                    
                              !  = 3 if temperature and pressure are given.                
                              !  = 4 if pressure and entropy are given.                    
                              !  = 5 if pressure and enthalpy are given. 

                                                                        
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


!..... Quantities expected to be found in /FLUDPC/
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
  REAL:: dxds, dxmax
  REAL:: h1
  REAL:: hin,hl,hg,h
  REAL:: delh
  REAL:: p1,p2, t1,t2, delp,dels,delt
  REAL:: small,tsat
  REAL:: s,s0,s1
  REAL:: sin,sl,sg
  REAL:: u0
  REAL:: huge,big   ! WATCH OUT
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
  INTEGER:: max2,max3,max4,max5
  INTEGER:: nr
  INTEGER:: nprops
  INTEGER:: nextp
  LOGICAL:: GAS 
     
  REAL:: ROOTS(2)


  REAL:: G(29,4),GG(27,30,8),TT(27),DD(30),SAT(40,8),C(8)
  INTEGER:: L(4) 
  INTEGER:: n1,n2,n3,ng,ns
  COMMON /FLUDPC/ G,TT,DD,GG,SAT,C,L,NG,N1,N2,N3,NS

  REAL:: GAMMA,WL,WG,DENSL,DENSG,ENTL,ENTG,ENTHL,ENTHG
  COMMON /FLUIDC/ GAMMA,WL,WG,DENSL,DENSG,ENTL,ENTG,ENTHL,ENTHG


      DATA MAX2,MAX3,MAX4,MAX5/12,12,12,12/
      DATA DXMAX/3.0/
      DATA HUGE,SMALL,BIG/1.0E37, 0.001, 1000.0/
     
! Arithmetic Statement Functions

      real tf,pf,hf,sf,cpf,a2f

      TF(P,D)=(P+3.0*D**2)*(3.0-D)/(8.0*Z*D) 
      PF(T,D)=8.0*Z*D*T/(3.0-D)-3.0*D**2 
      HF(P,D)=1.125*(P/(3.0*D)-D) 
      SF(D)=LOG((3.0-D)/D) 
      CPF(T,D)=1.0/(1.0-D*(3.0-D)**2/(4.0*T)) 
      A2F(T,D)=0.375*GAMMA*R*TC*(8.0*T/(3.0-D)*(1.0+D/(3.0-D))-6.0*D) 

!----------------------------------------------------------------------------

       IF (IGOTO.EQ.4) THEN            !.....INITIALIZE NORMALIZING VALUES.
         R=C(1) 
         TC=C(2) 
         PC=C(3) 
         DC=C(4) 
       ELSE IF (IGOTO.EQ.1) THEN 

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

         IF (ENTRY.EQ.1) THEN     !     CASE 1  GIVEN TEMPERATURE AND DENSITY. 

!     PROCEDURE (TEST FOR VAPOR CONDITION AT (T,D))                     
!     PROCEDURE (SET UP FOR INTERPOLATION OF SATURATION DATA AT T)      
!     PROCEDURE (GET SATURATED LIQUID (DL) AND GAS (DG) DENSITIES)      
          IF (T.LT.1.0) THEN 
            CALL NTRP(SAT(:,1),T,NS,L(1),KS,1,DUMF,DUMFF,DUMVAL) 
            CALL NTRP(DUMA,DUMX,DUMN0,DUMIA,DUMOUT,4,SAT(:,3), DUMFF,DL)                                             
            CALL NTRP(DUMA,DUMX,DUMN0,DUMIA,DUMOUT,4,SAT(:,4), DUMFF,DG)                                             
            VAPOR=D.GT.DG .AND. D.LT.DL 
          ELSE 
            VAPOR=.FALSE. 
          END IF

!     PROCEDURE (GET SATURATION PRESSURE, PSAT)
!     PROCEDURE (GET Z BY DOUBLE INTERPOLATION AT (T,D))                
!     PROCEDURE (SET UP FOR DOUBLE INTERPOLATION AT (T,D))              
          IF (VAPOR) THEN
            CALL NTRP(DUMA,DUMX,DUMN0,DUMIA,DUMOUT,4, SAT(:,2),DUMFF,PSAT)                                  
            P=PSAT 
          ELSE 
            IF (ARG1.NE.T) THEN 
              CALL NTRP(TT,T,N1,L(3),K1,2,DUMF,DUMFF,DUMVAL) 
              ARG1=T 
            END IF 
            IF (ARG2.NE.D) THEN 
              CALL NTRP(DD,D,N2,L(4),K2,3,DUMF,DUMFF,DUMVAL) 
              ARG2=D 
            END IF 
            CALL NTRP(DUMA,DUMX,DUMN0,DUMIA,DUMOUT,5,DUMF, GG(:,:,1),Z)  
            Z=RLIMIT(SMALL,Z,BIG) 
            P=PF(T,D) 
          END IF 
          PRES=P*PC 

         ELSE IF (ENTRY.EQ.2) THEN    ! CASE 2  .GIVEN PRESSURE AND DENSITY.
!     PROCEDURE (TEST FOR VAPOR CONDITION AT (P,D))                     
!     PROCEDURE (SET UP FOR INTERPOLATION OF SATURATION DATA AT P)      
!     PROCEDURE (GET SATURATED LIQUID (DL) AND GAS (DG) DENSITIES)      

        IF (P.LT.1.0) THEN 
          CALL NTRP(SAT(:,2),MAX(P,SAT(1,2)),NS,L(1),KS,1,DUMF,DUMFF,DUMVAL)                                       
          CALL NTRP(DUMA,DUMX,DUMN0,DUMIA,DUMOUT,4, SAT(:,3),DUMFF,DL)                                    
          CALL NTRP(DUMA,DUMX,DUMN0,DUMIA,DUMOUT,4, SAT(:,4),DUMFF,DG)                                    
          VAPOR=D.GT.DG .AND. D.LT.DL 
        ELSE 
          VAPOR=.FALSE. 
        END IF 

!     PROCEDURE (GET SATURATION TEMPERATURE, TSAT)
!     PROCEDURE (SOLVE FOR T AT (P,D) BY ITERATION)                     
        IF (VAPOR) THEN
          CALL NTRP(DUMA,DUMX,DUMN0,DUMIA,DUMOUT,4, SAT(:,1),DUMFF,TSAT)                                  
          T=TSAT 
        ELSE 
          Z=1.0 
          T=TF(P,D) 
          DELT=HUGE 
          DO 10 I2=1,MAX2 !     PROCEDURE (GET Z BY DOUBLE INTERPOLATION AT (T,D))                
            IF (ARG1.NE.T) THEN   !     PROCEDURE (SET UP FOR DOUBLE INTERPOLATION AT (T,D))              
              CALL NTRP(TT,T,N1,L(3),K1,2,DUMF,DUMFF,DUMVAL) 
              ARG1=T 
            END IF 
            IF (ARG2.NE.D) THEN 
              CALL NTRP(DD,D,N2,L(4),K2,3,DUMF,DUMFF,DUMVAL) 
              ARG2=D 
            END IF 
            CALL NTRP(DUMA,DUMX,DUMN0,DUMIA,DUMOUT,5,DUMF, GG(:,:,1),Z)                                       
            Z=RLIMIT(SMALL,Z,BIG) 
            IF (ABS(DELT).LE.C(5)*T) GO TO 20 
            T1=T 
            P1=PF(T1,D) 
            T=RLIMIT(TT(1),TF(P,D),TT(N1)) 
!     PROCEDURE (GET Z BY DOUBLE INTERPOLATION AT (T,D))                
!     PROCEDURE (SET UP FOR DOUBLE INTERPOLATION AT (T,D))              
            IF (ARG1.NE.T) THEN 
              CALL NTRP(TT,T,N1,L(3),K1,2,DUMF,DUMFF,DUMVAL) 
              ARG1=T 
            END IF 
            IF (ARG2.NE.D) THEN 
              CALL NTRP(DD,D,N2,L(4),K2,3,DUMF,DUMFF,DUMVAL) 
              ARG2=D 
            END IF 
            CALL NTRP(DUMA,DUMX,DUMN0,DUMIA,DUMOUT,5,DUMF, GG(:,:,1),Z)                                       
            Z=RLIMIT(SMALL,Z,BIG) 
            T2=T 
            P2=PF(T2,D) 
            DELP=P2-P1 
            IF (DELP.EQ.0.0) GO TO 20 
            DELT=(P-P2)/DELP*(T2-T1) 
            T=T2+DELT 
   10       CONTINUE 
            ERROR=BIT5 
   20       CONTINUE 
            END IF 
            TEMP=T*TC 
                                                                        
         ELSE IF (ENTRY.EQ.3) THEN !     CASE 3  ..GIVEN TEMPERATURE AND PRESSURE.

            VAPOR=.FALSE. 
            IF (T.LT.1.0) THEN 
!     PROCEDURE (GET FIRST GUESS FOR D FROM SATURATED VALUE)            
!     PROCEDURE (SET UP FOR INTERPOLATION OF SATURATION DATA AT T)      
        CALL NTRP(SAT(:,1),T,NS,L(1),KS,1,DUMF,DUMFF,DUMVAL) 
!     PROCEDURE (GET SATURATION PRESSURE, PSAT)                         
        CALL NTRP(DUMA,DUMX,DUMN0,DUMIA,DUMOUT,4, SAT(:,2),DUMFF,PSAT)                                  
               GAS=P.LT.PSAT 
!     PROCEDURE (GET SATURATED LIQUID (DL) AND GAS (DG) DENSITIES)      
        CALL NTRP(DUMA,DUMX,DUMN0,DUMIA,DUMOUT,4, SAT(:,3),DUMFF,DL)                                    
        CALL NTRP(DUMA,DUMX,DUMN0,DUMIA,DUMOUT,4, SAT(:,4),DUMFF,DG)                                    
               IF (GAS) THEN 
                  D=DG 
               ELSE 
                  D=DL 
               END IF 
            ELSE 
               GAS=.TRUE. 
               Z=1.0 
!     PROCEDURE (SOLVE CUBIC FOR D AT (T,P))                            
               CALL CUBIC (Z*T,P,ROOTS,NR) 
               IF (GAS) NR=1 
               D=ROOTS(NR) 
            END IF 
!     PROCEDURE (SOLVE FOR D AT (T,P) BY ITERATION)                     
            DELD=HUGE 
            DO 30 I3=1,MAX3 
!     PROCEDURE (GET Z BY DOUBLE INTERPOLATION AT (T,D))                
!     PROCEDURE (SET UP FOR DOUBLE INTERPOLATION AT (T,D))              
               IF (ARG1.NE.T) THEN 
!                 CALL NTRP1 (TT,T,N1,L(3),K1,2)                        
                  CALL NTRP(TT,T,N1,L(3),K1,2,DUMF,DUMFF,DUMVAL) 
                  ARG1=T 
               END IF 
               IF (ARG2.NE.D) THEN 
!                 CALL NTRP2 (DD,D,N2,L(4),K2,3)                        
                  CALL NTRP(DD,D,N2,L(4),K2,3,DUMF,DUMFF,DUMVAL) 
                  ARG2=D 
               END IF 
    CALL NTRP(DUMA,DUMX,DUMN0,DUMIA,DUMOUT,5,DUMF, GG(:,:,1),Z)                                          
               Z=RLIMIT(SMALL,Z,BIG) 
               IF (ABS(DELD).LE.C(6)*D) GO TO 40 
!              IF (ABS(DELD).LE.C(6)*D) EXIT(I3)                        
               D1=D 
               P1=PF(T,D1) 
!     PROCEDURE (SOLVE CUBIC FOR D AT (T,P))                            
               CALL CUBIC (Z*T,P,ROOTS,NR) 
               IF (GAS) NR=1 
               D=ROOTS(NR) 
               D=RLIMIT(DD(1),D,DD(N2))    ! ******* over/under ???
!     PROCEDURE (GET Z BY DOUBLE INTERPOLATION AT (T,D))                
!     PROCEDURE (SET UP FOR DOUBLE INTERPOLATION AT (T,D))              
               IF (ARG1.NE.T) THEN 
!                 CALL NTRP1 (TT,T,N1,L(3),K1,2)                        
                  CALL NTRP(TT,T,N1,L(3),K1,2,DUMF,DUMFF,DUMVAL) 
                  ARG1=T 
               END IF 
               IF (ARG2.NE.D) THEN 
!                 CALL NTRP2 (DD,D,N2,L(4),K2,3)                        
                  CALL NTRP(DD,D,N2,L(4),K2,3,DUMF,DUMFF,DUMVAL) 
                  ARG2=D 
               END IF 
      CALL NTRP(DUMA,DUMX,DUMN0,DUMIA,DUMOUT,5,DUMF, GG(:,:,1),Z)                                          
               Z=RLIMIT(SMALL,Z,BIG) 
               D2=D 
               P2=PF(T,D2) 
               DELP=P2-P1 
               IF (DELP.EQ.0.0) GO TO 40 
!              IF (DELP.EQ.0.0) EXIT (I3)                               
               DELD=(P-P2)/DELP*(D2-D1) 
               D=D2+DELD 
   30        CONTINUE 
!     OTHERWISE                                                         
            ERROR=BIT5 
   40       CONTINUE 
            DENS=D*DC 
                                                                        
         ELSE IF (ENTRY.EQ.4) THEN 
!     CASE 4                                                            
!.....GIVEN PRESSURE AND ENTROPY.                                       
            SIN=PROPS(2)/R 
            IF (P.LT.1.0) THEN 
!     PROCEDURE (TEST FOR VAPOR CONDITION AT (P,SIN))                   
!     PROCEDURE (SET UP FOR INTERPOLATION OF SATURATION DATA AT P)      
       CALL NTRP(SAT(:,2),MAX(P,SAT(1,2)),NS,L(1),KS,1,DUMF,DUMFF,DUMVAL)                                    
!     PROCEDURE (GET SATURATED LIQUID (SL) AND GAS (SG) ENTROPIES)      
     CALL NTRP(DUMA,DUMX,DUMN0,DUMIA,DUMOUT,4, SAT(:,5),DUMFF,SL)                                    
     CALL NTRP(DUMA,DUMX,DUMN0,DUMIA,DUMOUT,4, SAT(:,6),DUMFF,SG)                                    
               VAPOR=SIN.GT.SL .AND. SIN.LT.SG 
            ELSE 
               VAPOR=.FALSE. 
            END IF 
            IF (VAPOR) THEN 
!     PROCEDURE (GET (TWO-PHASE) T, WL, WG AND D AT (P,SIN))            
!     PROCEDURE (GET SATURATION TEMPERATURE, TSAT)                      
    CALL NTRP(DUMA,DUMX,DUMN0,DUMIA,DUMOUT,4, SAT(:,1),DUMFF,TSAT)                                  
!     PROCEDURE (GET SATURATED LIQUID (DL) AND GAS (DG) DENSITIES)      
!              CALL GNTRP (4,SAT(1,3),DL)                               
      CALL NTRP(DUMA,DUMX,DUMN0,DUMIA,DUMOUT,4, SAT(:,3),DUMFF,DL)                                    
      CALL NTRP(DUMA,DUMX,DUMN0,DUMIA,DUMOUT,4, SAT(:,4),DUMFF,DG)                                    
               T=TSAT 
               WL=(SIN-SG)/(SL-SG) 
               WG=1.0-WL 
               D=1.0/(WL/DL+WG/DG) 
            ELSE 
!     PROCEDURE (GET FIRST GUESS FOR D AT (P,SIN))                      
               GAS=SIN.GT.SAT(NS,5) 
               Z=1.0 
               IF (P.LT.1.0) THEN 
!     PROCEDURE (GET SATURATED LIQUID (DL) AND GAS (DG) DENSITIES)      
!                 CALL GNTRP (4,SAT(1,3),DL)                            
         CALL NTRP(DUMA,DUMX,DUMN0,DUMIA,DUMOUT,4, SAT(:,3),DUMFF,DL)                                 
         CALL NTRP(DUMA,DUMX,DUMN0,DUMIA,DUMOUT,4, SAT(:,4),DUMFF,DG)                                 
                  IF (GAS) THEN 
                     D=DG 
                     S1=SG 
                  ELSE 
                     D=DL 
                     S1=SL 
                  END IF 
               ELSE 
!     PROCEDURE (ASSUME WORST-CASE STARTING VALUE FOR DENSITY)          
                  T=TT(N1) 
!     PROCEDURE (SOLVE CUBIC FOR D AT (T,P))                            
                  CALL CUBIC (Z*T,P,ROOTS,NR) 
                  IF (GAS) NR=1 
                  D=ROOTS(NR) 
                  D=MAX(1.0,ROOTS(1)) 
!     PROCEDURE (SOLVE FOR T AT (P,D) BY ITERATION)                     
                  T=TF(P,D) 
                  DELT=HUGE 
                  DO 50 I2=1,MAX2 
!     PROCEDURE (GET Z BY DOUBLE INTERPOLATION AT (T,D))                
!     PROCEDURE (SET UP FOR DOUBLE INTERPOLATION AT (T,D))              
                     IF (ARG1.NE.T) THEN 
                     CALL NTRP(TT,T,N1,L(3),K1,2,DUMF,DUMFF,DUMVAL) 
                        ARG1=T 
                     END IF 
                     IF (ARG2.NE.D) THEN 
                        CALL NTRP(DD,D,N2,L(4),K2,3,DUMF,DUMFF,DUMVAL) 
                        ARG2=D 
                     END IF 
!                    CALL GGNTRP (5,GG(1,1,1),Z)                        
           CALL NTRP(DUMA,DUMX,DUMN0,DUMIA,DUMOUT,5,DUMF, GG(:,:,1),Z)                                    
                     Z=RLIMIT(SMALL,Z,BIG) 
                     IF (ABS(DELT).LE.C(5)*T) GO TO 60 
!                    IF (ABS(DELT).LE.C(5)*T) EXIT (I2)                 
                     T1=T 
                     P1=PF(T1,D) 
                     T=RLIMIT(TT(1),TF(P,D),TT(N1)) 
!     PROCEDURE (GET Z BY DOUBLE INTERPOLATION AT (T,D))                
!     PROCEDURE (SET UP FOR DOUBLE INTERPOLATION AT (T,D))              
                     IF (ARG1.NE.T) THEN 
!                       CALL NTRP1(TT,T,N1,L(3),K1,2)                   
                        CALL NTRP(TT,T,N1,L(3),K1,2,DUMF,DUMFF,DUMVAL) 
                        ARG1=T 
                     END IF 
                     IF (ARG2.NE.D) THEN 
!                       CALL NTRP2(DD,D,N2,L(4),K2,3)                   
                        CALL NTRP(DD,D,N2,L(4),K2,3,DUMF,DUMFF,DUMVAL) 
                        ARG2=D 
                     END IF 
!                    CALL GGNTRP (5,GG(1,1,1),Z)                        
      CALL NTRP(DUMA,DUMX,DUMN0,DUMIA,DUMOUT,5,DUMF, GG(:,:,1),Z)                                    
                     Z=RLIMIT(SMALL,Z,BIG) 
                     T2=T 
                     P2=PF(T2,D) 
                     DELP=P2-P1 
                     IF (DELP.EQ.0.0) GO TO 60 
!                    IF (DELP.EQ.0.0) EXIT (I2)                         
                     DELT=(P-P2)/DELP*(T2-T1) 
                     T=T2+DELT 
   50             CONTINUE 
!                 OTHERWISE                                             
      ERROR=BIT5 
   60             CONTINUE 
!     PROCEDURE (GET S0 BY INTERPOLATION OF ZERO-DENSITY DATA)          
!     PROCEDURE (SET UP FOR INTERPOLATION OF ZERO-DENSITY DATA AT T)    
                  IF (ARG0.NE.T) THEN 
          CALL NTRP(G(:,1),T,NG,L(2),K0,1,DUMF,DUMFF,DUMVAL) 
                     ARG0=T 
                  END IF 
!                 CALL GNTRP (4,G(1,2),S0)                              
         CALL NTRP(DUMA,DUMX,DUMN0,DUMIA,DUMOUT,4, G(:,2),DUMFF,S0)                                   
!     PROCEDURE (GET S BY DOUBLE INTERPOLATION FOR Z2)                  
!     PROCEDURE (SET UP FOR DOUBLE INTERPOLATION AT (T,D))              
                  IF (ARG1.NE.T) THEN 
!                    CALL NTRP1 (TT,T,N1,L(3),K1,2)                     
               CALL NTRP(TT,T,N1,L(3),K1,2,DUMF,DUMFF,DUMVAL) 
                     ARG1=T 
                  END IF 
                  IF (ARG2.NE.D) THEN 
!                    CALL NTRP2 (DD,D,N2,L(4),K2,3)                     
                   CALL NTRP(DD,D,N2,L(4),K2,3,DUMF,DUMFF,DUMVAL) 
                     ARG2=D 
                  END IF 
!                 CALL GGNTRP (5,GG(1,1,2),Z2)                          
      CALL NTRP(DUMA,DUMX,DUMN0,DUMIA,DUMOUT,5,DUMF, GG(:,:,2),Z2)                                      
                  S=Z2*(S0+SF(D)) 
                  S1=S 
               END IF 
!     PROCEDURE (SOLVE FOR T AND D AT (P,SIN) BY ITERATION)             
               IF (GAS) THEN 
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
!     PROCEDURE (SOLVE FOR T AT (P,D) BY ITERATION)                     
                  T=TF(P,D) 
                  DELT=HUGE 
                  DO 70 I2=1,MAX2 
!     PROCEDURE (GET Z BY DOUBLE INTERPOLATION AT (T,D))                
!     PROCEDURE (SET UP FOR DOUBLE INTERPOLATION AT (T,D))              
                     IF (ARG1.NE.T) THEN 
!                       CALL NTRP1(TT,T,N1,L(3),K1,2)                   
                        CALL NTRP(TT,T,N1,L(3),K1,2,DUMF,DUMFF,DUMVAL) 
                        ARG1=T 
                     END IF 
                     IF (ARG2.NE.D) THEN 
!                       CALL NTRP2(DD,D,N2,L(4),K2,3)                   
                        CALL NTRP(DD,D,N2,L(4),K2,3,DUMF,DUMFF,DUMVAL) 
                        ARG2=D 
                     END IF 
!                    CALL GGNTRP (5,GG(1,1,1),Z)                        
        CALL NTRP(DUMA,DUMX,DUMN0,DUMIA,DUMOUT,5,DUMF, GG(:,:,1),Z)                                    
                     Z=RLIMIT(SMALL,Z,BIG) 
                     IF (ABS(DELT).LE.C(5)*T) GO TO 80 
!                    IF (ABS(DELT).LE.C(5)*T) EXIT (I2)                 
                     T1=T 
                     P1=PF(T1,D) 
                     T=RLIMIT(TT(1),TF(P,D),TT(N1)) 
!     PROCEDURE (GET Z BY DOUBLE INTERPOLATION AT (T,D))                
!     PROCEDURE (SET UP FOR DOUBLE INTERPOLATION AT (T,D))              
                     IF (ARG1.NE.T) THEN 
!                       CALL NTRP1(TT,T,N1,L(3),K1,2)                   
                        CALL NTRP(TT,T,N1,L(3),K1,2,DUMF,DUMFF,DUMVAL) 
                        ARG1=T 
                     END IF 
                     IF (ARG2.NE.D) THEN 
!                       CALL NTRP2(DD,D,N2,L(4),K2,3)                   
                        CALL NTRP(DD,D,N2,L(4),K2,3,DUMF,DUMFF,DUMVAL) 
                        ARG2=D 
                     END IF 
!                    CALL GGNTRP(5,GG(1,1,1),Z)                        
        CALL NTRP(DUMA,DUMX,DUMN0,DUMIA,DUMOUT,5,DUMF, GG(:,:,1),Z)                                    
                     Z=RLIMIT(SMALL,Z,BIG) 
                     T2=T 
                     P2=PF(T2,D) 
                     DELP=P2-P1 
                     IF (DELP.EQ.0.0) GO TO 80 
!                    IF (DELP.EQ.0.0) EXIT (I2)                         
                     DELT=(P-P2)/DELP*(T2-T1) 
                     T=T2+DELT 
   70             CONTINUE 
!                 OTHERWISE                                             
                  ERROR=BIT5 
   80             CONTINUE 
                  IF (ABS(DELS).LE.C(7)) GO TO 100 
!                 IF (ABS(DELS).LE.C(7)) EXIT (I4)                      
                  IF (X.NE.X1) THEN 
                     ERROR=0 
!     PROCEDURE (GET S0 BY INTERPOLATION OF ZERO-DENSITY DATA)          
!     PROCEDURE (SET UP FOR INTERPOLATION OF ZERO-DENSITY DATA AT T)    
                     IF (ARG0.NE.T) THEN 
          CALL NTRP(G(:,1),T,NG,L(2),K0,1,DUMF,DUMFF, DUMVAL)                                      
                        ARG0=T 
                     END IF 
!                    CALL GNTRP (4,G(1,2),S0)                           
          CALL NTRP(DUMA,DUMX,DUMN0,DUMIA,DUMOUT,4, G(:,2),DUMFF,S0)                                
!     PROCEDURE (GET S BY DOUBLE INTERPOLATION FOR Z2)                  
!     PROCEDURE (SET UP FOR DOUBLE INTERPOLATION AT (T,D))              
                     IF (ARG1.NE.T) THEN 
!                       CALL NTRP1(TT,T,N1,L(3),K1,2)                   
                        CALL NTRP(TT,T,N1,L(3),K1,2,DUMF,DUMFF,DUMVAL) 
                        ARG1=T 
                     END IF 
                     IF (ARG2.NE.D) THEN 
!                       CALL NTRP2(DD,D,N2,L(4),K2,3)                   
                        CALL NTRP(DD,D,N2,L(4),K2,3,DUMF,DUMFF,DUMVAL) 
                        ARG2=D 
                     END IF 
!                    CALL GGNTRP(5,GG(1,1,2),Z2)                       
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
   90          CONTINUE 
!              OTHERWISE                                                
               ERROR=BIT5 
  100          CONTINUE 
            END IF 
            DENS=D*DC 
            TEMP=T*TC 
                                                                        
         ELSE IF (ENTRY.EQ.5) THEN 
!     CASE 5                                                            
!.....GIVEN PRESSURE AND ENTHALPY.                                      
            HIN=PROPS(3)/R/TC 
            IF (P.LT.1.0) THEN 
!     PROCEDURE (TEST FOR VAPOR CONDITION AT (P,HIN))                   
!     PROCEDURE (SET UP FOR INTERPOLATION OF SATURATION DATA AT P)      
      CALL NTRP(SAT(:,2),AMAX1(P,SAT(1,2)),NS,L(1),KS,1, DUMF,DUMFF,DUMVAL)                                    
!     PROCEDURE (GET SATURATED LIQUID (HL) AND GAS (HG) ENTHALPIES)     
!                             CALL GNTRP (4,SAT(1,7),HL)                
        CALL NTRP(DUMA,DUMX,DUMN0,DUMIA,DUMOUT,4, SAT(:,7),DUMFF,HL)                                    
        CALL NTRP(DUMA,DUMX,DUMN0,DUMIA,DUMOUT,4, SAT(:,8),DUMFF,HG)                                    
               VAPOR=HIN.GT.HL .AND. HIN.LT.HG 
            ELSE 
               VAPOR=.FALSE. 
            END IF 
            IF (VAPOR) THEN 
!     PROCEDURE (GET (TWO-PHASE) T, WL, WG AND D AT (P,HIN))            
!     PROCEDURE (GET SATURATION TEMPERATURE, TSAT)                      
!                             CALL GNTRP (4,SAT(1,1),TSAT)              
        CALL NTRP(DUMA,DUMX,DUMN0,DUMIA,DUMOUT,4, SAT(:,1),DUMFF,TSAT)                                  
!     PROCEDURE (GET SATURATED LIQUID (DL) AND GAS (DG) DENSITIES)      
!                              CALL GNTRP (4,SAT(1,3),DL)               
        CALL NTRP(DUMA,DUMX,DUMN0,DUMIA,DUMOUT,4, SAT(:,3),DUMFF,DL)                                    
        CALL NTRP(DUMA,DUMX,DUMN0,DUMIA,DUMOUT,4, SAT(:,4),DUMFF,DG)                                    
               T=TSAT 
               WL=(HIN-HG)/(HL-HG) 
               WG=1.0-WL 
               D=1.0/(WL/DL+WG/DG) 
            ELSE 
!     PROCEDURE (GET FIRST GUESS FOR D AT (P,HIN))                      
               Z=1.0 
               IF (P.LT.1.0) THEN 
                  GAS=HG.LT.HIN 
!     PROCEDURE (GET SATURATED LIQUID (DL) AND GAS (DG) DENSITIES)      
!                              CALL GNTRP (4,SAT(1,3),DL)               
   CALL NTRP(DUMA,DUMX,DUMN0,DUMIA,DUMOUT,4, SAT(:,3),DUMFF,DL)                                 
   CALL NTRP(DUMA,DUMX,DUMN0,DUMIA,DUMOUT,4, SAT(:,4),DUMFF,DG)                                 
                  IF (GAS) THEN 
                     D=DG 
                     H1=HG 
                  ELSE 
                     D=DL 
                     H1=HL 
                  END IF 
               ELSE 
!     PROCEDURE (ASSUME WORST-CASE STARTING VALUE FOR DENSITY)          
                  T=TT(N1) 
!     PROCEDURE (SOLVE CUBIC FOR D AT (T,P))                            
                  CALL CUBIC (Z*T,P,ROOTS,NR) 
                  IF (GAS) NR=1 
                  D=ROOTS(NR) 
                  D=AMAX1(1.0,ROOTS(1)) 
!     PROCEDURE (SOLVE FOR T AT (P,D) BY ITERATION)                     
                  T=TF(P,D) 
                  DELT=HUGE 
                  DO 110 I2=1,MAX2 
!     PROCEDURE (GET Z BY DOUBLE INTERPOLATION AT (T,D))                
!     PROCEDURE (SET UP FOR DOUBLE INTERPOLATION AT (T,D))              
                     IF (ARG1.NE.T) THEN 
!                       CALL NTRP1(TT,T,N1,L(3),K1,2)                   
                        CALL NTRP(TT,T,N1,L(3),K1,2,DUMF,DUMFF,rDUMOUT) 
                        ARG1=T 
                     END IF 
                     IF (ARG2.NE.D) THEN 
!                       CALL NTRP2(DD,D,N2,L(4),K2,3)                   
                        CALL NTRP(DD,D,N2,L(4),K2,3,DUMF,DUMFF,rDUMOUT) 
                        ARG2=D 
                     END IF 
!                    CALL GGNTRP (5,GG(1,1,1),Z)                        
           CALL NTRP(DUMA,DUMX,DUMN0,DUMIA,DUMOUT,5,DUMF, GG(:,:,1),Z)                                    
                     Z=RLIMIT(SMALL,Z,BIG) 
                     IF (ABS(DELT).LE.C(5)*T) GO TO 120 
!                    IF (ABS(DELT).LE.C(5)*T) EXIT (I2)                 
                     T1=T 
                     P1=PF(T1,D) 
                     T=RLIMIT(TT(1),TF(P,D),TT(N1)) 
!     PROCEDURE (GET Z BY DOUBLE INTERPOLATION AT (T,D))                
!     PROCEDURE (SET UP FOR DOUBLE INTERPOLATION AT (T,D))              
                     IF (ARG1.NE.T) THEN 
                        CALL NTRP(TT,T,N1,L(3),K1,2,DUMF,DUMFF,DUMVAL) 
                        ARG1=T 
                     END IF 
                     IF (ARG2.NE.D) THEN 
!                       CALL NTRP2(DD,D,N2,L(4),K2,3)                   
                        CALL NTRP(DD,D,N2,L(4),K2,3,DUMF,DUMFF,DUMVAL) 
                        ARG2=D 
                     END IF 
!                    CALL GGNTRP (5,GG(1,1,1),Z)                        
      CALL NTRP(DUMA,DUMX,DUMN0,DUMIA,DUMOUT,5,DUMF, GG(:,:,1),Z)                                    
                     Z=RLIMIT(SMALL,Z,BIG) 
                     T2=T 
                     P2=PF(T2,D) 
                     DELP=P2-P1 
                     IF (DELP.EQ.0.0) GO TO 120 
!                    IF (DELP.EQ.0.0) EXIT (I2)                         
                     DELT=(P-P2)/DELP*(T2-T1) 
                     T=T2+DELT 
  110             CONTINUE 
!                 OTHERWISE                                             
                  ERROR=BIT5 
  120             CONTINUE 
!     PROCEDURE (GET U0 BY INTERPOLATION OF ZERO-DENSITY DATA)          
!     PROCEDURE (SET UP FOR INTERPOLATION OF ZERO-DENSITY DATA AT T)    
                  IF (ARG0.NE.T) THEN 
                     CALL NTRP(G(:,1),T,NG,L(2),K0,1,DUMF,DUMFF,DUMVAL) 
                     ARG0=T 
                  END IF 
!                 CALL GNTRP (4,G(1,3),U0)                              
          CALL NTRP(DUMA,DUMX,DUMN0,DUMIA,DUMOUT,4, G(:,3),DUMFF,U0)                                   
!     PROCEDURE (GET H BY DOUBLE INTERPOLATION FOR Z3)                  
!     PROCEDURE (SET UP FOR DOUBLE INTERPOLATION AT (T,D))              
                  IF (ARG1.NE.T) THEN 
!                    CALL NTRP1 (TT,T,N1,L(3),K1,2)                     
                     CALL NTRP(TT,T,N1,L(3),K1,2,DUMF,DUMFF,DUMVAL) 
                     ARG1=T 
                  END IF 
                  IF (ARG2.NE.D) THEN 
!                    CALL NTRP2 (DD,D,N2,L(4),K2,3)                     
                     CALL NTRP(DD,D,N2,L(4),K2,3,DUMF,DUMFF,DUMVAL) 
                     ARG2=D 
                  END IF 
!                 CALL GGNTRP (5,GG(1,1,1),Z3)                          
        CALL NTRP(DUMA,DUMX,DUMN0,DUMIA,DUMOUT,5,DUMF, GG(:,:,3),Z3)                                      
                  H=Z3*(U0+HF(P,D)) 
                  H1=H 
                  GAS=H1.LT.HIN 
               END IF 
!     PROCEDURE (SOLVE FOR T AND D AT (P,HIN) BY ITERATION)             
               IF (GAS) THEN 
                  X1=1.0/D 
                  D=D/1.2 
                  X=1.0/D 
               ELSE 
                  X1=D 
                  D=1.01*D 
                  X=D 
               END IF 
               DELH=HUGE 
               DO 150 I5=1,MAX5 
!     PROCEDURE (SOLVE FOR T AT (P,D) BY ITERATION)                     
                  T=TF(P,D) 
                  DELT=HUGE 
                  DO 130 I2=1,MAX2 
!     PROCEDURE (GET Z BY DOUBLE INTERPOLATION AT (T,D))                
!     PROCEDURE (SET UP FOR DOUBLE INTERPOLATION AT (T,D))              
                     IF (ARG1.NE.T) THEN 
!                       CALL NTRP1(TT,T,N1,L(3),K1,2)                   
                        CALL NTRP(TT,T,N1,L(3),K1,2,DUMF,DUMFF,DUMVAL) 
                        ARG1=T 
                     END IF 
                     IF (ARG2.NE.D) THEN 
!                       CALL NTRP2(DD,D,N2,L(4),K2,3)                   
                        CALL NTRP(DD,D,N2,L(4),K2,3,DUMF,DUMFF,DUMVAL) 
                        ARG2=D 
                     END IF 
!                    CALL GGNTRP (5,GG(1,1,1),Z)                        
      CALL NTRP(DUMA,DUMX,DUMN0,DUMIA,DUMOUT,5,DUMF, GG(:,:,1),Z)                                    
                     Z=RLIMIT(SMALL,Z,BIG) 
                     IF (ABS(DELT).LE.C(5)*T) GO TO 140 
!                    IF (ABS(DELT).LE.C(5)*T) EXIT (I2)                 
                     T1=T 
                     P1=PF(T1,D) 
                     T=RLIMIT(TT(1),TF(P,D),TT(N1)) 
!     PROCEDURE (GET Z BY DOUBLE INTERPOLATION AT (T,D))                
!     PROCEDURE (SET UP FOR DOUBLE INTERPOLATION AT (T,D))              
                     IF (ARG1.NE.T) THEN 
!                       CALL NTRP1(TT,T,N1,L(3),K1,2)                   
                        CALL NTRP(TT,T,N1,L(3),K1,2,DUMF,DUMFF,rDUMOUT) 
                        ARG1=T 
                     END IF 
                     IF (ARG2.NE.D) THEN 
!                       CALL NTRP2(DD,D,N2,L(4),K2,3)                   
                        CALL NTRP(DD,D,N2,L(4),K2,3,DUMF,DUMFF,rDUMOUT) 
                        ARG2=D 
                     END IF 
!                    CALL GGNTRP (5,GG(1,1,1),Z)                        
         CALL NTRP(DUMA,DUMX,DUMN0,DUMIA,DUMOUT,5,DUMF, GG(:,:,1),Z)                                    
                     Z=RLIMIT(SMALL,Z,BIG) 
                     T2=T 
                     P2=PF(T2,D) 
                     DELP=P2-P1 
                     IF (DELP.EQ.0.0) GO TO 140 
!                    IF (DELP.EQ.0.0) EXIT (I2)                         
                     DELT=(P-P2)/DELP*(T2-T1) 
                     T=T2+DELT 
  130             CONTINUE 
!                 OTHERWISE                                             
                  ERROR=BIT5 
  140             CONTINUE 
                  IF (ABS(DELH).LE.C(8)) GO TO 160 
!                 IF (ABS(DELH).LE.C(8)) EXIT (I5)                      
                  IF (X.NE.X1) THEN 
                     ERROR=0 
!     PROCEDURE (GET U0 BY INTERPOLATION OF ZERO-DENSITY DATA)          
!     PROCEDURE (SET UP FOR INTERPOLATION OF ZERO-DENSITY DATA AT T)    
                     IF (ARG0.NE.T) THEN 
                CALL NTRP(G(:,1),T,NG,L(2),K0,1,DUMF,DUMFF, DUMVAL)                                      
                        ARG0=T 
                     END IF 
!                    CALL GNTRP (4,G(1,3),U0)                           
           CALL NTRP(DUMA,DUMX,DUMN0,DUMIA,DUMOUT,4, G(:,3),DUMFF,U0)                                
!     PROCEDURE (GET H BY DOUBLE INTERPOLATION FOR Z3)                  
!     PROCEDURE (SET UP FOR DOUBLE INTERPOLATION AT (T,D))              
                     IF (ARG1.NE.T) THEN 
!                       CALL NTRP1(TT,T,N1,L(3),K1,2)                   
                        CALL NTRP(TT,T,N1,L(3),K1,2,DUMF,DUMFF,DUMVAL) 
                        ARG1=T 
                     END IF 
                     IF (ARG2.NE.D) THEN 
!                       CALL NTRP2(DD,D,N2,L(4),K2,3)                   
                        CALL NTRP(DD,D,N2,L(4),K2,3,DUMF,DUMFF,DUMVAL) 
                        ARG2=D 
                     END IF 
!                    CALL GGNTRP (5,GG(1,1,1),Z3)                       
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
  150          CONTINUE 
!              OTHERWISE                                                
               ERROR=BIT5 
  160          CONTINUE 
            END IF 
            DENS=D*DC 
            TEMP=T*TC 
                                                                        
         END IF 
                                                                        
!.....CALCULATE REMAINING PROPERTIES AT (T,D).                          
         IF (T.GT.0.0 .AND. D.GT.0.0) THEN 
            IF (VAPOR) THEN 
!     PROCEDURE (CALCULATE REMAINING PROPERTIES FOR 2-PHASE FLUID)      
               NPROPS=MIN0(3,NP) 
               DO 170 IPROP=1,NPROPS 
!      DO CASE (IPROP,3)                                                
                                       IF (IPROP.EQ.1) THEN 
!     CASE 1                                                            
                     IF (ENTRY.LE.3) THEN 
                        WL=DL/D*(D-DG)/(DL-DG) 
                        WG=1.0-WL 
                     END IF 
                     DENSL=DL*DC 
                     DENSG=DG*DC 
                     PROPS(1)=PRES/(R*DENS*TEMP) 
                                       ELSE IF (IPROP.EQ.2) THEN 
!     CASE 2                                                            
                     IF (ENTRY.NE.4) THEN 
!     PROCEDURE (GET SATURATED LIQUID (SL) AND GAS (SG) ENTROPIES)      
!                         CALL GNTRP (4,SAT(1,5),SL)                    
        CALL NTRP(DUMA,DUMX,DUMN0,DUMIA,DUMOUT,4, SAT(:,5),DUMFF,SL)                           
        CALL NTRP(DUMA,DUMX,DUMN0,DUMIA,DUMOUT,4, SAT(:,6),DUMFF,SG)                           
                        PROPS(2)=(WL*SL+WG*SG)*R 
                     END IF 
                     ENTL=SL*R 
                     ENTG=SG*R 
                                       ELSE IF (IPROP.EQ.3) THEN 
!     CASE 3                                                            
                     IF (ENTRY.NE.5) THEN 
!     PROCEDURE (GET SATURATED LIQUID (HL) AND GAS (HG) ENTHALPIES)     
   CALL NTRP(DUMA,DUMX,DUMN0,DUMIA,DUMOUT,4, SAT(:,7),DUMFF,HL)                           
   CALL NTRP(DUMA,DUMX,DUMN0,DUMIA,DUMOUT,4, SAT(:,8),DUMFF,HG)                           
                        PROPS(3)=(WL*HL+WG*HG)*R*TC 
                     END IF 
                     ENTHL=HL*R*TC 
                     ENTHG=HG*R*TC 
                  END IF 
  170          CONTINUE 
            ELSE 
               NEXTP=1 
!     PROCEDURE (CALCULATE REMAINING PROPERTIES FOR HOMOGENEOUS FLUID)  
               DO 180 IPROP=NEXTP,NP 
                  IF (IPROP.LE.6) THEN 
!     DO CASE (IPROP,6)                                                 
                               IF (IPROP.EQ.1) THEN 
!     CASE 1                                                            
                        PROPS(1)=PRES/(R*DENS*TEMP) 
                               ELSE IF (IPROP.EQ.2) THEN 
!     CASE 2                                                            
                        IF (ENTRY.NE.4) THEN 
!     PROCEDURE (GET S0 BY INTERPOLATION OF ZERO-DENSITY DATA)          
!     PROCEDURE (SET UP FOR INTERPOLATION OF ZERO-DENSITY DATA AT T)    
                           IF (ARG0.NE.T) THEN 
        CALL NTRP(G(:,1),T,NG,L(2),K0,1,DUMF,DUMFF,DUMVAL)                      
                              ARG0=T 
                           END IF 
!               CALL GNTRP (4,G(1,2),S0)                                
       CALL NTRP(DUMA,DUMX,DUMN0,DUMIA,DUMOUT,4,G(:,2),DUMFF,S0)                             
!     PROCEDURE (GET S BY DOUBLE INTERPOLATION FOR Z2)                  
!     PROCEDURE (SET UP FOR DOUBLE INTERPOLATION AT (T,D))              
           IF (ARG1.NE.T) THEN 
!                             CALL NTRP1 (TT,T,N1,L(3),K1,2)            
           CALL NTRP(TT,T,N1,L(3),K1,2,DUMF, DUMFF,DUMVAL)                             
                              ARG1=T 
                           END IF 
                           IF (ARG2.NE.D) THEN 
!                             CALL NTRP2 (DD,D,N2,L(4),K2,3)            
        CALL NTRP(DD,D,N2,L(4),K2,3,DUMF, DUMFF,DUMVAL)                             
                              ARG2=D 
                           END IF 
!                          CALL GGNTRP (5,GG(1,1,2),Z2)                 
    CALL NTRP(DUMA,DUMX,DUMN0,DUMIA,DUMOUT,5,DUMF,GG(:,:,2),Z2)                           
                           S=Z2*(S0+SF(D)) 
                           PROPS(2)=S*R 
                        END IF 
                               ELSE IF (IPROP.EQ.3) THEN 
!     CASE 3                                                            
                        IF (ENTRY.NE.5) THEN 
!     PROCEDURE (GET U0 BY INTERPOLATION OF ZERO-DENSITY DATA)          
!     PROCEDURE (SET UP FOR INTERPOLATION OF ZERO-DENSITY DATA AT T)    
                           IF (ARG0.NE.T) THEN 
         CALL NTRP(G(:,1),T,NG,L(2),K0,1, DUMF,DUMFF,DUMVAL)                        
                              ARG0=T 
                           END IF 
!                          CALL GNTRP (4,G(1,3),U0)                     
          CALL NTRP(DUMA,DUMX,DUMN0,DUMIA,DUMOUT,4, G(:,3),DUMFF,U0)                             
!     PROCEDURE (GET H BY DOUBLE INTERPOLATION FOR Z3)                  
!     PROCEDURE (SET UP FOR DOUBLE INTERPOLATION AT (T,D))              
                           IF (ARG1.NE.T) THEN 
!                             CALL NTRP1 (TT,T,N1,L(3),K1,2)            
      CALL NTRP(TT,T,N1,L(3),K1,2,  DUMF,DUMFF,DUMVAL)                        
                              ARG1=T 
                           END IF 
                           IF (ARG2.NE.D) THEN 
!                             CALL NTRP2 (DD,D,N2,L(4),K2,3)            
       CALL NTRP(DD,D,N2,L(4),K2,3, DUMF,DUMFF,DUMVAL)                        
                              ARG2=D 
                           END IF 
!                          CALL GGNTRP (5,GG(1,1,3),Z3)                 
      CALL NTRP(DUMA,DUMX,DUMN0,DUMIA,DUMOUT,5, DUMF,GG(:,:,3),Z3)                           
                           H=Z3*(U0+HF(P,D)) 
                           PROPS(3)=H*R*TC 
                        END IF 
                                       ELSE IF (IPROP.EQ.4) THEN 
!     CASE 4                                                            
!     PROCEDURE (GET CV0 BY INTERPOLATION OF ZERO-DENSITY DATA)         
!                       CALL GNTRP (4,G(1,4),CV0)                       
        CALL NTRP(DUMA,DUMX,DUMN0,DUMIA,DUMOUT,4, G(:,4),DUMFF,CV0)                               
!     PROCEDURE (GET CV BY DOUBLE INTERPOLATION FOR Z4)                 
!                       CALL GGNTRP (5,GG(1,1,4),Z4)                    
         CALL NTRP(DUMA,DUMX,DUMN0,DUMIA,DUMOUT,5, DUMF,GG(:,:,4),Z4)                              
                        CV=Z4*CV0 
                        PROPS(4)=CV*R 
                               ELSE IF (IPROP.EQ.5) THEN 
!     CASE 5                                                            
!     PROCEDURE (GET CP BY DOUBLE INTERPOLATION FOR Z5)                 
!                       CALL GGNTRP (5,GG(1,1,5),Z5)                    
      CALL NTRP(DUMA,DUMX,DUMN0,DUMIA,DUMOUT,5, DUMF,GG(:,:,5),Z5)                              
                        CP=Z5*(CV+CPF(T,D)) 
                        PROPS(5)=CP*R 
                        GAMMA=CP/CV 
                                       ELSE IF (IPROP.EQ.6) THEN 
!     CASE 6                                                            
!     PROCEDURE (GET SONIC VELOCITY BY DOUBLE INTERPOLATION FOR Z6)     
!                       CALL GGNTRP (5,GG(1,1,6),Z6)                    
       CALL NTRP(DUMA,DUMX,DUMN0,DUMIA,DUMOUT,5,DUMF,GG(:,:,6),Z6)                              
                        A2=A2F(T,D) 
                        IF (A2.GT.0.0) THEN 
                           A=Z6*SQRT(A2) 
                        ELSE 
                           A=0.0 
                        END IF 
                        PROPS(6)=A 
                     END IF 
                  ELSE 
!     PROCEDURE (DOUBLE INTERPOLATION FOR OTHER PROPERTIES)             
!                    CALL GGNTRP (5,GG(1,1,IPROP),PROPS(IPROP))         
      CALL NTRP(DUMA,DUMX,DUMN0,DUMIA,DUMOUT,5,DUMF,GG(:,:,IPROP),PROPS(IPROP))                   
                  END IF 
  180          CONTINUE 
            END IF 
         END IF 
         ERROR=KS*BIT1+K0*BIT2+K1*BIT3+K2*BIT4+ERROR 
!     RETURN                                                            
                                                                        
      ELSE IF (IGOTO.EQ.2.OR.IGOTO.EQ.3) THEN 
               IF (IGOTO.EQ.2) THEN 
!     ENTRY FLUIDL (PROPS,NP,ERROR)                                     
!.....ENTRY FOR GETTING LIQUID-PHASE PROPERTIES NOT IN COMMON.          
            D=DL 
!     RETURN                                                            
                                                                        
               ELSE IF (IGOTO.EQ.3) THEN 
!     ENTRY FLUIDG (PROPS,NP,ERROR)                                     
!.....ENTRY FOR GETTING GAS-PHASE PROPERTIES NOT IN COMMON.             
            D=DG 
               END IF 
!     PROCEDURE (CALCULATE OTHER SATURATION-LOCUS PROPERTIES)           
         IF (.NOT.VAPOR) THEN 
            WRITE(*,*) "Fatal Error in Fluid. Vapor is FALSE"
            STOP "Fatal Error #3 in Fluid"
         END IF 
         IF (NP.GT.N3) THEN 
            WRITE(*,*) "Fatal error in Fluid. np > n3 ", NP,N3
            STOP "Fatal error #4 in Fluid"
         END IF 
!     PROCEDURE (SET UP FOR INTERPOLATION OF ZERO-DENSITY DATA AT T)    
         IF (ARG0.NE.T) THEN 
            CALL NTRP(G(:,1),T,NG,L(2),K0,1,DUMF,DUMFF,DUMVAL) 
            ARG0=T 
         END IF 
!     PROCEDURE (SET UP FOR DOUBLE INTERPOLATION AT (T,D))              
         IF (ARG1.NE.T) THEN 
!           CALL NTRP1 (TT,T,N1,L(3),K1,2)                              
            CALL NTRP(TT,T,N1,L(3),K1,2,DUMF,DUMFF,DUMVAL) 
            ARG1=T 
         END IF 
         IF (ARG2.NE.D) THEN 
!           CALL NTRP2 (DD,D,N2,L(4),K2,3)                              
            CALL NTRP(DD,D,N2,L(4),K2,3,DUMF,DUMFF,DUMVAL) 
            ARG2=D 
         END IF 
         IF (T.GT.0.0 .AND. NP.GE.4) THEN 
            NEXTP=4 
!     PROCEDURE (CALCULATE REMAINING PROPERTIES FOR HOMOGENEOUS FLUID)  
            DO 190 IPROP=NEXTP,NP 
               IF (IPROP.LE.6) THEN 
!     DO CASE (IPROP,6)                      ! SELECT CASE would be better coding..
                                  IF (IPROP.EQ.1) THEN 
!     CASE 1                                                            
                     PROPS(1)=PRES/(R*DENS*TEMP) 
                               ELSE IF (IPROP.EQ.2) THEN 
!     CASE 2                                                            
                     IF (ENTRY.NE.4) THEN 
!     PROCEDURE (GET S0 BY INTERPOLATION OF ZERO-DENSITY DATA)          
!     PROCEDURE (SET UP FOR INTERPOLATION OF ZERO-DENSITY DATA AT T)    
                        IF (ARG0.NE.T) THEN 
         CALL NTRP(G(:,1),T,NG,L(2),K0,1,DUMF,DUMFF,DUMVAL)                                
                           ARG0=T 
                        END IF 
!                       CALL GNTRP (4,G(1,2),S0)                        
         CALL NTRP(DUMA,DUMX,DUMN0,DUMIA,DUMOUT,4, G(:,2),DUMFF,S0)                                
!     PROCEDURE (GET S BY DOUBLE INTERPOLATION FOR Z2)                  
!     PROCEDURE (SET UP FOR DOUBLE INTERPOLATION AT (T,D))              
                        IF (ARG1.NE.T) THEN 
!                          CALL NTRP1 (TT,T,N1,L(3),K1,2)               
         CALL NTRP(TT,T,N1,L(3),K1,2,DUMF,DUMFF, DUMVAL)                                      
                           ARG1=T 
                        END IF 
                        IF (ARG2.NE.D) THEN 
!                          CALL NTRP2 (DD,D,N2,L(4),K2,3)               
        CALL NTRP(DD,D,N2,L(4),K2,3,DUMF,DUMFF, DUMVAL)                                      
                           ARG2=D 
                        END IF 
!                       CALL GGNTRP (5,GG(1,1,2),Z2)                    
        CALL NTRP(DUMA,DUMX,DUMN0,DUMIA,DUMOUT,5, DUMF,GG(:,:,2),Z2)                              
                        S=Z2*(S0+SF(D)) 
                        PROPS(2)=S*R 
                     END IF 
                               ELSE IF (IPROP.EQ.3) THEN 
!     CASE 3                                                            
                     IF (ENTRY.NE.5) THEN 
!     PROCEDURE (GET U0 BY INTERPOLATION OF ZERO-DENSITY DATA)          
!     PROCEDURE (SET UP FOR INTERPOLATION OF ZERO-DENSITY DATA AT T)    
                        IF (ARG0.NE.T) THEN 
          CALL NTRP(G(:,1),T,NG,L(2),K0,1,DUMF, DUMFF,DUMVAL)                                
                           ARG0=T 
                        END IF 
!                       CALL GNTRP (4,G(1,3),U0)                        
           CALL NTRP(DUMA,DUMX,DUMN0,DUMIA,DUMOUT,4,G(:,3),DUMFF,U0)                                
!     PROCEDURE (GET H BY DOUBLE INTERPOLATION FOR Z3)                  
!     PROCEDURE (SET UP FOR DOUBLE INTERPOLATION AT (T,D))              
                        IF (ARG1.NE.T) THEN 
!                          CALL NTRP1 (TT,T,N1,L(3),K1,2)               
        CALL NTRP(TT,T,N1,L(3),K1,2,DUMF,DUMFF, DUMVAL)                                      
                           ARG1=T 
                        END IF 
                        IF (ARG2.NE.D) THEN 
!                          CALL NTRP2 (DD,D,N2,L(4),K2,3)               
                  CALL NTRP(DD,D,N2,L(4),K2,3,DUMF,DUMFF, DUMVAL)                                      
                           ARG2=D 
                        END IF 
!                       CALL GGNTRP (5,GG(1,1,3),Z3)                    
         CALL NTRP(DUMA,DUMX,DUMN0,DUMIA,DUMOUT,5,DUMF,GG(:,:,3),Z3)                              
                        H=Z3*(U0+HF(P,D)) 
                        PROPS(3)=H*R*TC 
                     END IF 
                                  ELSE IF (IPROP.EQ.4) THEN 
!     CASE 4                                                            
!     PROCEDURE (GET CV0 BY INTERPOLATION OF ZERO-DENSITY DATA)         
!                    CALL GNTRP(4,G(1,4),CV0)                           
      CALL NTRP(DUMA,DUMX,DUMN0,DUMIA,DUMOUT,4, G(:,4),DUMFF,CV0)                                  
!     PROCEDURE (GET CV BY DOUBLE INTERPOLATION FOR Z4)                 
!                    CALL GGNTRP (5,GG(1,1,4),Z4)                       
 CALL NTRP(DUMA,DUMX,DUMN0,DUMIA,DUMOUT,5, DUMF,GG(:,:,4),Z4)                                 
                     CV=Z4*CV0 
                     PROPS(4)=CV*R 
                               ELSE IF (IPROP.EQ.5) THEN 
!     CASE 5                                                            
!     PROCEDURE (GET CP BY DOUBLE INTERPOLATION FOR Z5)                 
!                    CALL GGNTRP (5,GG(1,1,5),Z5)                       
        CALL NTRP(DUMA,DUMX,DUMN0,DUMIA,DUMOUT,5, DUMF,GG(:,:,5),Z5)                                 
                     CP=Z5*(CV+CPF(T,D)) 
                     PROPS(5)=CP*R 
                     GAMMA=CP/CV 
                               ELSE IF (IPROP.EQ.6) THEN 
!     CASE 6                                                            
!     PROCEDURE (GET SONIC VELOCITY BY DOUBLE INTERPOLATION FOR Z6)     
!                    CALL GGNTRP (5,GG(1,1,6),Z6)                       
          CALL NTRP(DUMA,DUMX,DUMN0,DUMIA,DUMOUT,5,DUMF,GG(:,:,6),Z6)                                 
                     A2=A2F(T,D) 
                     IF (A2.GT.0.0) THEN 
                        A=Z6*SQRT(A2) 
                     ELSE 
                        A=0.0 
                     END IF 
                     PROPS(6)=A 
                  END IF 
               ELSE 
!     PROCEDURE (DOUBLE INTERPOLATION FOR OTHER PROPERTIES)             
!                 CALL GGNTRP (5,GG(1,1,IPROP),PROPS(IPROP))            
      CALL NTRP(DUMA,DUMX,DUMN0,DUMIA,DUMOUT,5, DUMF,GG(:,:,IPROP),PROPS(IPROP))                      
               END IF 
  190       CONTINUE 
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

CONTAINS

!+
FUNCTION Limit(imin, i, imax) RESULT(k)
  INTEGER,INTENT(IN):: imin,i,imax
  INTEGER:: k
  k=MAX(IMIN,MIN(I,IMAX))
  RETURN
END Function Limit   ! ======================================================



!+
FUNCTION RLIMIT(XMIN,X,XMAX) RESULT(f)
  REAL,INTENT(IN):: xmin,x,xmax
  REAL:: f
!----------------------------------------------------------------------------
  f=MAX(XMIN,MIN(X,XMAX))
  RETURN
END Function Rlimit   ! =====================================================


!+
SUBROUTINE NTRP (A,X,N0,IA,OUT,IGOTO,F,FF,VALUE)
! ---------------------------------------------------------------------------
! PURPOSE - DOUBLE 3-POINT INTERPOLATION ROUTINES FOR SUBPROGRAM FLUID.
!      WRITTEN BY: THEODORE E. FESSLER, NASA TM X-3572, AUGUST 1977
!      REWRITTEN BY: LT. MARK D. KLEM & MARGARET P. PROCTOR, JANUARY 198
!                                                                       


  IMPLICIT NONE
  REAL,INTENT(IN),DIMENSION(:):: a   ! TABLE OF INDEPENDENT-VARIABLE VALUES
  REAL,INTENT(IN):: x   ! ENTRY VALUE OF INDEPENDENT VARIABLE
  INTEGER,INTENT(IN):: n0   !  (OR N1/N2) = LENGTH OF A TABLE
  INTEGER,INTENT(OUT):: ia   !  CURRENT INDEX VALUE ASSOCIATED WITH A TABLE
  INTEGER,INTENT(OUT):: out  !  = 0,1 IF X IS INSIDE/NOT INSIDE LIMITS OF A
  INTEGER,INTENT(IN):: igoto ! ???
  REAL,INTENT(IN),DIMENSION(:):: f
  REAL,INTENT(IN),DIMENSION(:,:):: ff
  REAL,INTENT(OUT):: value

!      THE ORIGINAL PROGRAM WAS WRITTEN IN SFTRAN, BUT THIS PROGRAM IS  
!      REWRITTEN IN STANDARD FORTRAN 77 TO RUN PRIMARILY ON PERSONAL    
!      COMPUTERS AND CAN ALSO RUN ON MAINFRAMES. MUCH OF THE OLD CODE WA
!      COMMENTED OUT BUT LEFT INTACT FOR DOCUMENTATION AND REFERENCE    
!      PURPOSES.

! NOTES - This routine is an example of a situation that requires the
!   SAVE attribute of Fortran 90. The routine must first be called with
!   IGOTO = 1,2, or 3. Then, a later call with IGOTO=4 or 5 actually
!   produces the result. BUT, certain values must be saved between
!   these calls in order for this to work. Most Fortran IV compilers of
!   the early days saved all the variables and many programs were written
!   with this as a basic assumption. Since modern Fortran uses variables
!   with the AUTOMATIC attribute, routines such as this will fail
!   unless the appropriate variables are declared SAVE.

  REAL,SAVE,DIMENSION(4,3):: c
  REAL,SAVE,DIMENSION(4):: c0,c1,c2     ! SAVE is redundant (EQUIVALENCE)
  REAL,SAVE,DIMENSION(16):: cc 
  INTEGER,SAVE,DIMENSION(3):: i1,i2,nn
  LOGICAL,SAVE:: NEW12

  INTEGER:: i,j,ij,k,l,m   ! various counters
  INTEGER:: l1,l2, m1,m2, mm1,mm2
  INTEGER:: na
  INTEGER:: ni,nj


  EQUIVALENCE (c0(1),c(1,1)), (c1(1),c(1,2)), (c2,c(1,3))
      ! so, c0 is col. 1 of c, c1 is col. 2, c2 is col.3

  EQUIVALENCE (L1,i1(1)), (m1,i1(2)), (mm1,i1(3))
  EQUIVALENCE (L2,i2(1)), (m2,i2(2)), (mm2,i2(3))
  EQUIVALENCE             (ni,nn(2)), (nj,nn(3))

  REAL:: p1,p2,p3
  REAL:: a11,a12,a13,a14, a22,a23,a24, a33,a34, a44
!----------------------------------------------------------------------------                                                                        

  IF (IGOTO .EQ. 4) THEN      ! function of one variable
    VALUE=0.0 
    K=0 
    DO I=L1,L2 
      K=K+1 
      VALUE=VALUE+C0(K)*F(I)
    END DO
    RETURN 
  END IF                                                                      

  IF (IGOTO .EQ. 5) THEN 
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
        VALUE=VALUE+CC(IJ)*FF(I,J) 
      END DO   
    END DO   
    RETURN 
  END IF 

  IF (IGOTO .EQ. 1) THEN 
    K=1 
    NA=N0 
  ELSE IF (IGOTO .EQ. 2) THEN 
    K=2 
    NA=N0 
    NEW12=.TRUE. 
  ELSE IF (IGOTO .EQ. 3) THEN 
    K=3 
    NA=N0 
    NEW12=.TRUE. 
  END IF 

  IF (na > SIZE(a)) WRITE(*,*) "NTRP called with inconsistent data"
  IF (na < 1 ) THEN
    WRITE(*,*) "NTRP called with n0 <= 0;   Set to 1 for safety"
    na=1
  END IF


  L=Limit(1,NA,3) 
  M=NA-2 
  IF (X .LT. A(1)) THEN     !     PROCEDURE (TABLE LOOK-UP)                                         
    OUT=1 
    IA=1 
  ELSE IF (X .GT. A(NA)) THEN 
    OUT=1 
    IA=NA 
  ELSE 
    OUT=0 
    IA=LIMIT(1,IA,NA) 
    IF (X .LT. A(IA)) THEN 
10    IA=IA-1                         ! look backwards
      IF (X .LT. A(IA)) GO TO 10 
    ELSE 
20    IF (X .GT. A(IA+1)) THEN 
        IA=IA+1                        ! look forwards
        GO TO 20 
      END IF 
    END IF 
  END IF

  IF (IA.GE.2 .AND. IA.LE.M) L=4 
  I=LIMIT(1,IA-1,M) 
  I1(K)=I 
  I2(K)=I+L-1 
  NN(K)=L 

  IF (L.GT.1) THEN      !     PROCEDURE (CALCULATE AIJ VALUES)                                  
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

!      WRITTEN BY: THEODORE E. FESSLER, NASA TM X-3572, AUGUST 1977
!      REWRITTEN BY: LT. MARK D. KLEM & MARGARET P. PROCTOR, JANUARY 198
!                                                                       
  IMPLICIT NONE

  REAL(KIND=4),INTENT(IN):: t,p
  REAL(KIND=4),DIMENSION(:),INTENT(OUT):: roots   ! dimension=2
  INTEGER,INTENT(OUT):: nroots
                                                                        

                                                                        
                                                                        
!!      IMPLICIT REAL*8 (A-H,O-Z) 
!!      REAL*8 K1,K2 
!!      REAL*4 T,P,ROOTS(2) 
  REAL(KIND=8),PARAMETER:: K1=2.094395102393195D+0
  REAL(KIND=8),PARAMETER:: K2=4.188790204786391D+0
  REAL(KIND=8),PARAMETER:: ZERO=0.0, ONE=1.0, TWO=2.0, THREE=3.0, EIGHT=8.0

!!!      DATA THIRD / 3.333333333333333D-1 / 
  REAL(KIND=8):: a,a3,aa,bb,b,b2,cr
  REAL(KIND=8):: rad,phi
  REAL(KIND=8):: q
  REAL(KIND=8):: r1,r2,r3
                                                                        
!.....MAIN PROGRAM FLOW.                                                
!     PROCEDURE (REDUCE TO NORMAL FORM)                                 
      Q=(EIGHT*T+P)/THREE
      A=Q/THREE-ONE
      A3=A*A*A 
      B=(Q-P)/TWO-ONE
      B2=B*B 
      RAD=A3+B2 
      IF (RAD < ZERO) THEN 
!        PROCEDURE (CASE FOR 3 DISTINCT REAL ROOTS)                     
         PHI=ACOS(SIGN(SQRT(-B2/A3),-B))/THREE
         CR=TWO*SQRT(-A) 
         R1=CR*COS(PHI) 
         R2=CR*COS(PHI+K1) 
         R3=CR*COS(PHI+K2) 
         ROOTS(1)=MIN(R1,R2,R3)+ONE
         ROOTS(2)=MAX(R1,R2,R3)+ONE
         NROOTS=2 
      ELSE 
         IF (RAD .EQ. ZERO) THEN 
            IF (B .EQ. ZERO) THEN 
!              PROCEDURE (CASE FOR ONLY ONE REAL ROOT)                  
               IF (RAD .NE. ZERO) RAD=SQRT(RAD) 
               AA=-B+RAD 
               BB=-B-RAD 
               IF (AA .NE. ZERO) AA=QRT(AA) 
               IF (BB .NE. ZERO) BB=QRT(BB) 
               ROOTS(1)=AA+BB+ONE
               NROOTS=1 
            ELSE 
!              PROCEDURE (CASE FOR 3 REAL ROOTS BUT TWO ARE EQUAL)      
               R1=QRT(B) 
               R2=-(R1+R1)
               ROOTS(1)=MIN(R1,R2)+ONE
               ROOTS(2)=MAX(R1,R2)+ONE
               NROOTS=2 
            END IF 
         ELSE 
!           PROCEDURE (CASE FOR ONLY ONE REAL ROOT)                     
            IF (RAD.NE.0.0D0) RAD=SQRT(RAD) 
            AA=-B+RAD 
            BB=-B-RAD 
            IF (AA .NE. ZERO) AA=QRT(AA) 
            IF (BB .NE. ZERO) BB=QRT(BB) 
            ROOTS(1)=AA+BB+ONE
            NROOTS=1 
         END IF 
      END IF 

  RETURN
END Subroutine Cubic   ! ====================================================

!+
FUNCTION Qrt(arg) RESULT(f)
! ---------------------------------------------------------------------------
  REAL(KIND=8),INTENT(IN):: arg
  REAL(KIND=8):: f
  REAL(KIND=8),PARAMETER:: THIRD=1.0D0/3.0D0
!----------------------------------------------------------------------------
  f=SIGN(ABS(ARG)**THIRD, ARG) 
  RETURN
END Function Qrt   ! ========================================================


END Subroutine Fluid   ! ====================================================
