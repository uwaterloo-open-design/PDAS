!+
! PROGRAM GASP
! ------------------------------------------------------------------------------
! PURPOSE - Compute the thermodynamic and transport properties of helium, 
!  methane, neon, nitrogen, carbon monoxide, carbon dioxide, oxygen, and argon.


! AUTHORS - Robert C. Hendricks, Anne K. Baron, and 
!           Ildiko C. Peller, NASA Glenn (Lewis) Research Center
!           Ralph L. Carmichael, Public Domain Aeronautical Software

! REVISION HISTORY
!   DATE  VERS PERSON  STATEMENT OF CHANGES
!   1975   1.0    "    Publication of NASA TN D-7808 (February 1975)
! 11Aug97  1.1   RLC   Acquisition of COSMIC Program LEW-11629
! 02Sep08  1.2   RLC   Replaced ASSIGN statements to conform to Fortran 95
!
!******************  SUBMISSION  LEW-11629  ****************************
!C$IBFTC GASPY
!
!     ------------------NASA  LEWIS  RESEARCH---------------------------
!     ---------------------------VERSION 2/1/72-------------------------
!     -----------FOR  INFORMATION --SEE  A.BARON, R.C.HENDRICKS,I.PELLER
      SUBROUTINE GASP(KS,KP,T,P,D,H,KR)
!      BENDER-S EQUATION OF STATE FOR THE FIVE GASES--N2,CH4,AR,O2,CO2
!       ALSO EQUATIONS OF STATE FOR CO, NE, AND HELIUM...
! -------------------------- PRYDZ EQUATION FOR FLUORINE ---------------
!
!        COMPUTE THE STATE RELATIONS AND THERMODYNAMIC AND
!        TRANSPORT PROPERITES OF SPECIFIED FLUID GIVEN TEMPERATURE T,
!        PRESSURE P, DENSITY D, OR ENTHALPY H.  STATE RELATIONS ARE
!        SPECIFIED BY KS.  THERMODYNAMIC AND TRANSPORT PROPERTIES
!        ARE SPECIFIED BY KP.  IF KR IS RETURNED OR SPECIFIED AS 1,
!        PROPERTIES ARE COMPUTED AT SATURATION.
!
      DIMENSION KPC1(32), KPC2(32), KPC3(32),KPC4(32)
      COMMON/PROPTY/KU,DL,DV,HL,HV,S,SL,SV,CV,CVL,CVV,CP,CPL,CPV,GAMMA, &
     & GAMMAL,GAMMAV,C,CL,CVP,MU,MUL,MUV,K,KL,KV,SIGMA,EXCL,EXCV,EXCESK
      REAL MU,MUL,MUV,K,KL,KV
      COMMON/DERIV/PDT,PTV,PDTL,PDTV,PTVL,PTVV
      COMMON/PARTLS/PTV1,PDT1
      COMMON/CHECKS/DCH1,DCH2,PCH1,PCH2,PCH3,TCH1,TCH2,TCH3,DST,TST,    &
     & HSCH1,HSCH2
      COMMON/HELFLU/ IHE,IFY,IHY
      DATA KPC1 /2,3,6,7,10,11,14,15,18,19,22,23,26,27,30,31,34,35,38,  &
     &39,42,43,46,47,50,51,54,55,58,59,62,63/
      DATA KPC2 /4,5,6,7,12,13,14,15,20,21,22,23,28,29,30,31,36,37,38,  &
     &39,44,45,46,47,52,53,54,55,60,61,62,63/
      DATA KPC3 /8,9,10,11,12,13,14,15,24,25,26,27,28,29,30,31,40,41,42,&
     &43,44,45,46,55,56,57,58,59,60,61,62,63/
      DATA KPC4 /16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,       &
     &48,49,50,51,52,53,54,55,56,57,58,59,60,61,62,63/
      GO TO (10,20,30,40,45),KS
!
                                    ! compute density
   10 CALL DENS(KU,T,P,D,DL,DV,KR)
      GO TO 50
!
                                 ! compute pressure
   20 CALL PRESS(KU,T,D,P,KR)
      GO TO 50
!
                               ! compute temperature
   30 CALL TEMP(KU,P,D,T,KR)
      GO TO 50
!
                                        ! compute temperature & density
   40 CALL TEMPPH(KU,P,H,T,D,DL,DV,KR)
                                        ! given pressure & enthalpy
      GO TO 50
!
                                          !compute temperature & density
   45 CALL TEMPPS ( KU,P,S,T,D,DL,DV,KR )
                                          ! given pressure & entropy
   50 IF (MOD(KP,2)) 60,70,60
!
                                       ! compute enthalpy
   60 CALL ENTH(KU,KR,T,P,D,H,HL,HV)
   70 DO 81 I=1,32
      IF (KP-KPC1(I)) 110,100,80
   80 CONTINUE
   81 END DO
      GO TO 110
!
!        COMPUTE ENTROPY
!
  100 CALL ENT(KU,KR,T,P,D,S,SL,SV)
  110 DO 121 I=1,32
      IF(KP-KPC2(I)) 140,130,120
  120 CONTINUE
  121 END DO
      GO TO 140
  130 KCP=0
!
!        COMPUTE  SPECIFIC  HEATS  AND  SONIC  VELOCITY
!
      IF (KR.NE.1) GO TO 200
      CALL DENS(KU,T,P,D,DL,DV,1)
      CALL CPPRL(P,DL,T,CPL,CVL,KU,KR,KCP)
      PTVL=PTV1
      PDTL=PDT1
      GAMMAL=CPL/CVL
      CALL SONIC(KU,KR,T,DL,GAMMAL,CL)
      CALL CPPRL(P,DV,T,CPV,CVV,KU,KR,KCP)
      PDTV=PDT1
      PTVV=PTV1
      GAMMAV=CPV/CVV
      CALL SONIC(KU,KR,T,DV,GAMMAV,CVP)
      GO TO 140
                                         ! ????? in right spot????
  200 CALL CPPRL(P,D,T,CP,CV,KU,KR,KCP)
      PDT=PDT1
      PTV=PTV1
      GAMMA=CP/CV
      CALL SONIC(KU,KR,T,D,GAMMA,C)
  140 DO 151 I=1,32
      IF (KP-KPC3(I)) 170,160,150
  150 CONTINUE
  151 END DO
      GO TO 170
!
  160 IF (KR.NE.1) GO TO 165
                                        ! compute viscosity
      CALL DENS(KU,T,P,D,DL,DV,1)
      CALL VISC(KU,KR,T,DL,MUL)
      CALL VISC(KU,KR,T,DV,MUV)
      GO TO 170
  165 CALL VISC(KU,KR,T,D,MU)
  170 DO 176 I=1,32
      IF(KP-KPC4(I)) 190,180,175
  175 CONTINUE
  176 END DO
      GO TO 190
!
!        COMPUTE THERMAL CONDUCTIVITY
!        NOTE-- FROZEN  VALUE  AVAILABLE  IN K,KL,KY
!           --  REACTING  CONDUCTIVITY  RETURNED  IN  EXCESK, EXCL,EXCV
!
  180 IF (KR.NE.1) GO TO 220
      CALL DENS(KU,T,P,D,DL,DV,1)
      CALL THERM (KU,KR,P,T,DL,EXCL,KL)
      CALL THERM (KU,KR,P,T,DV,EXCV,KV)
       GO TO 190
  220 CALL THERM (KU,KR,P,T,D,EXCESK,K)
  190 IF(KP-32) 230,240,240
!
                                ! compute surface tension
  240 CALL SURF (KU,KR,T,SIGMA)
  230 RETURN
      END
!$IBFTC BLOCD
!     ---------------------------VERSION 2/1/72-------------------------
!
!   STORES THE COEFFICIENTS FOR ALL FLUIDS FOR THE EQUATION OF STATE
!    AND THE TRANSPORT EQUATIONS.
!        STORES  CONVERSION  CONSTANTS  NEEDED  BY  ALL  FLUIDS
!
      BLOCK DATA
      CHARACTER*3 MATCH(10)
!      CHARACTER*90 MESSAG(10)
      CHARACTER(LEN=90),DIMENSION(10):: MESSAG = (/                     &
     & ' THERMODYNAMIC AND TRANSPORT PROPERTIES FOR'//                  &
     & ' METHANE PC=45.66 ATM,TC=190.77 K,ROC=.162 G/CC',               &
     & ' THERMODYNAMIC AND TRANSPORT PROPERTIES FOR'//                  &
     & ' NITROGEN PC=33.72ATM,TC=126.3 K,ROC=.3105 G/CC',               &
     & ' THERMODYNAMIC AND TRANSPORT PROPERTIES FOR'//                  &
     & ' OXYGEN PC=50.16ATM, TC=154.78 K,ROC=.4325 G/CC',               &
     & ' THERMODYNAMIC AND TRANSPORT PROPERTIES FOR'//                  &
     & ' ARGON PC= 48.014ATM,  TC=150.7 K,ROC=.531 G/CC',               &
     & ' THERMODYNAMIC AND TRANSPORT PROPERTIES FOR'//                  &
     & ' CO2 PC=72.869ATM,  TC=304.21 K,  ROC=.464 G/CC',               &
     & ' THERMODYNAMIC AND TRANSPORT PROPERTIES FOR'//                  &
     & ' NEON PC= 26.19ATM , TC=44.4 K,  RHOC=.483 G/CC',               &
     & ' THERMODYNAMIC AND TRANSPORT PROPERTIES FOR'//                  &
     & ' CO  PC=34.529ATM, TC=132.91 K, RHOC=.2997 G/CC',               &
     & ' THERMODYNAMIC AND TRANSPORT PROPERTIES FOR'//                  &
     & ' HELIUM   PC = 2ATM,  TC=5 K,     RHOC=1/25G/CC',               &
     & ' THERMODYNAMIC AND TRANSPORT PROPERTIES FOR'//                  &
     & ' FLUORINE  PC=54ATM, TC=144.31K, RHOC=.3974G/CC',               &
     & ' THERMODYNAMIC AND TRANSPORT PROPERTIES FOR'//                  &
     & ' H2     PC=12.759ATM,TC=32.976K,RHOC=.03143G/CC'  /)

      COMMON/GASES/match,messag
!!!   COMMON/GASES/MATCH(10),MESSAG(15,10)
      COMMON/CON123/ DCONV(5),TCONV(5),PCONV(5)
      COMMON/CONV4/SCONV(5)
      COMMON/CONV5/CCONV(5)
      COMMON/CONV6/HCONV(5)
      COMMON/ALLCOF/COF(57,10),TCCOF(15,10),AKSTCO(18,4),DIST(10),      &
     & WTMOL(10),EPSOK(10),ZETA(10),FF(10),SWT(10),KSWT(10),DIFTT(10),  &
     & RHOSWT(10),DTRIPL(10)
!!!      DIMENSION MCH4(15),MN2(15),MO2(15),MAR(15),MCO2(15),MNE(15),MCO
!!!     1),MHE(15),MF2(15),MH2(15)
!!!      CHARACTER*90 MCH4,MN2,MO2,MAR,MCO2,MNE,MCO,MHE,MF2,MH2
!!!      EQUIVALENCE(MESSAG(1),MCH4),(MESSAG(2),MN2), (MESSAG(3),MO2)
!!!      EQUIVALENCE(MESSAG(4),MAR), (MESSAG(5),MCO2),(MESSAG(6),MNE)
!!!      EQUIVALENCE(MESSAG(7),MCO), (MESSAG(8),MHE), (MESSAG(9),MF2)
!!!      EQUIVALENCE(MESSAG(10),MH2)
      DATA match/'CH4', 'N2 ', 'O2 ', 'AR ', 'CO2',                     &
     & 'NE ', 'CO ', 'HE ', 'F2 ', 'H2 ' /
!!!      DATA MCH4       /' THERMODYNAMIC AND TRANSPORT PROPERTIES FOR'//  &
!!!     & ' METHANE PC=45.66 ATM,TC=190.77 K,ROC=.162 G/CC'/
!!!      DATA MN2         /90H THERMODYNAMIC AND TRANSPORT PROPERTIES FOR N&
!!!     &ITROGEN PC=33.72ATM,TC=126.3 K,ROC=.3105G/CC  /
!!!      DATA MO2        /90H THERMODYNAMIC AND TRANSPORT PROPERTIES FOR OX&
!!!     &YGEN PC=50.16ATM,TC=154.78 K,ROC=.4325G/CC   /
!!!      DATA MAR        /90H THERMODYNAMIC AND TRANSPORT PROPERTIES FOR AR&
!!!     &GON PC= 48.014ATM, TC=150.7 K,ROC=.531 G/CC  /
!!!      DATA MCO2       /90H THERMODYNAMIC AND TRANSPORT PROPERTIES FOR CO&
!!!     & 2 PC=72.869ATM, TC=304.21 K, ROC=.464 G/CC  /
!!!      DATA MNE        /90H THERMODYNAMIC AND TRANSPORT PROPERTIES FOR NE&
!!!     &ON PC= 26.19ATM ,TC=44.4 K, RHOC=.483 G/CC   /
!!!      DATA MCO        /90H THERMODYNAMIC AND TRANSPORT PROPERTIES FOR CO&
!!!     &  PC=34.529ATM, TC=132.91 K, RHOC=.2997 G/CC /
!!!      DATA MHE        /90H THERMODYNAMIC AND TRANSPORT PROPERTIES FOR HE&
!!!     & PC = 2ATM, TC=5 K, RHOC=1/25G/CC                   /
!!!      DATA MF2        /90HTHERMODYNAMIC AND TRANSPORT PROPERTIES FOR FLU&
!!!     &ORINE PC=54ATM,TC=144.31K,RHOC=.3974G/CC               /
!!!      DATA MH2   /   90HTHERMODYNAMIC + TRANSPORT PROPERTIES FOR HYDROGE&
!!!     &N PC=12.759ATM,TC=32.976K,RHOC=.03143G/CC      /

      DATA TCONV/1., 1., 1.8, 2*1./
      DATA PCONV/1., 9.8692327, 145.03774 ,2*1./
      DATA DCONV/2*1., 62.42796, 2*1./
      DATA SCONV/2*1., 0.23900574, 2*1./
      DATA CCONV/2*1., 0.0328084, 2*1./
      DATA HCONV/2*1., 0.4302103, 2*1./


!        METHANE-------------------BENDER 1971
      DATA (COF(I,1),I=1,57)/       .518251          , .17191020E1     ,&
     &-.86366402E3   , .25005236E5  ,-.12533848E8    , .34169547E9     ,&
     & .75523689     ,-.12111233E3  , .20547188E6    , .32337540E2     ,&
     &-.61948317E4   ,-.25603803E2  , .11556713E4    , .27425297E5     ,&
     & .49499630E8   ,-.11956135E11 , .88302981E12   ,-.31713486E10    ,&
     & .12302028E13  ,-.64499295E14 , .37E2  ,8*0.  , -.32097996E2     ,&
     & .29976933E3, .70773219      , -.72349122E-2    , .40851153E-4   ,&
     &-.12130439E-6, .148489295E-9   , .11167E3      ,7.96986252       ,&
     & .39225023E-2  ,-.56781803E-4  , .23442607E-6  ,-.22395007E-9    ,&
     &8.35517346     , .803207E3     , .260815       , .1E-4           ,&
     & .57           , .01167        ,4.627          , .5067E2         ,&
     & .9066E2       , .19077E3      , .6E3          , .45             ,&
     & .25E3         , .0            , .2E4          /
!        NITROGEN------------------BENDER 1971
      DATA (COF(I,2),I=1,57)/       .296797          , .48058175       ,&
     &-.15047010E3   , .26071365E4  ,-.12792742E7    , .29436228E8     ,&
     & .37500265     ,-.50738465E2  , .14499236E5    , .1440998E1      ,&
     &-.24136776E3   , .28954771    ,-.27613799E2    , .36048264E3     ,&
     &-.20083357E7   , .43265184E9  ,-.16513521E11   ,-.10141439E8     ,&
     & .47216396E10  ,-.12016418E12 , .1E2 , 8*0. ,   -.66869126E1     ,&
     &-.20086668E3 , .37077481       ,-.62602129E-2  , .56187159E-4   , &
     &-.25984366E-6 , .49025967E-9  , .77364E2 ,  .29090035E+2      ,   &
     & .89930170E-3, -.92394038E-5, .26008592E-7, -.14102926E-10,       &
     & .23351167E1  , .231188578E3 , .35693888E-1  , .1449E-3  ,        &
     &1.1208         ,  .010132     ,3.417          , .5067E2   ,       &
     & .64E2   , .1263E3  , .100E4   , .92           ,                  &
     & .18E3   , .0           , .8E3  /
!        OXYGEN--------------------BENDER 1970
      DATA(COF(I,3),I=1,57) /          .259832       , .34811077       ,&
     &-.14070678E3   , .25061744E4   ,-.10081345E7   , .19074164E8     ,&
     &-.40134966E-2  , .65172112E2   , .10962206E5   , .6972158        ,&
     &-.26242449E3   , .191378       , .29416771E2  , .78932076E2     , &
     &-.1923158E7    , .4610824E9    ,-.39936263E11  ,-.5668995E7      ,&
     & .13644286E10  , .91977198E11  ,5.4 , 8*0.  ,  -.51504418E1     , &
     & -.25626822E3, .26280757   , -.36222681E-2    , .26516635E-4     ,&
     &-.10004384E-6  ,.15423155E-9   , .90E2         , .29145189E2     ,&
     &-.5691986E-3   , .68032634E-5  ,-.48456604E-7  , .13822647E-9    ,&
     &4.72426938     , .35642E3      , .03125117     , .15E-3          ,&
     &  1.50         , .0101325      ,5.083          ,1.0000E2         ,&
     & .543507E2     , .15478E3      , .5E+3         ,1.3              ,&
     & .23E3         , .0            , .9E3          /
!        ARGON---------------------BENDER 1970
      DATA (COF(I,4),I=1,57)/         .208128    ,     .19825921       ,&
     &-.81733119E2   , .1777747E4    ,-.82406544E6   , .31666098E8     ,&
     &-.44202671E-1  , .6216142E2    , .11443248E4   , .4779752        ,&
     &-.19645227E3   ,-.21572754     , .16544141E3   ,-.28142112E2     ,&
     & .82532059E5   ,-.91538377E7   ,-.18340752E10  ,-.33858136E7     ,&
     & .15532886E10  ,-.67479568E11 ,3.5  , 8*0.    , .13814448E+2    , &
     &-.57396331E+3  , -.18072036E+0 , .18350292E-2  ,-.10957549E-4    ,&
     & .36123787E-7  , -.50097228E-10, .8728E2       ,2.5              ,&
     & .0            , .0            , .0            , .0              ,&
     & 2.757177  , .237932E3         , .208152       , .1E-3           ,&
     &1.48           , .0101325      ,4.865          , .5066E2         ,&
     & .8378E2       , .1507E3       , .1E4         ,1.3               ,&
     & .25E3         , .0            , .25E4         /
!        CARBON DIOXIDE------------BENDER 1970
      DATA(COF (I,5)  ,I=1,57)/          .188918       , .22488558     ,&
     &-.13717965E3    ,-.14430214E5    ,-.29630491E7   ,-.20606039E9   ,&
     & .45554393E-1   , .77042840E2    , .40602371E5   , .40029509     ,&
     &-.39436077E3    , .12115286      , .10783386E3   , .43962336E2   ,&
     &-.36505545E8    , .19490511E11   ,-.29186718E13  , .24358627E8   ,&
     &-.37546530E11   , .11898141E14   ,  5. ,8*0.     ,-.95913112E+2  ,&
     & .26738653E+4  , .11405373E+1  , -.67026553E-2 , .21745017E-4    ,&
     &-.37153665E-7  , .26225608E-10   , .19471E3     ,5.152          , &
     & .015224        ,-.9681E-5       , .2313E-8      ,0.             ,&
     &4.0851019       , .72665417E3    , .095067       , .83E-5        ,&
     &1.25            ,0.010132        ,7.3835         ,.50669E2       ,&
     & .21656E3       , .30421E3       , .100E4        ,1.1            ,&
     & .37E3        ,0.              , .25E4         /
!        NEON----------------------MCCARTY AND STEWART 1965
      DATA (COF(I,6),I=1,57)/        .41185435       , .42879277       ,&
     &-.41617333E2 , -.15368388E4    , .21705583E5   , .0              ,&
     &.5237778E-1    , .31832376E2  ,.15304725E4     , .86991276       ,&
     &-.86565802E2   , .0            , .0            , .41994281E2     ,&
     & .23856729E6   ,-.10444366E8   , .0           ,-.19828098E7      ,&
     & .70965701E8   , .0             ,+.13590153E2  ,0.0              ,&
     &-.12553439E-2  , -.76888021E3  ,-.84938358E3,4*0., .67422985E1   ,&
     &-.11786837E3   ,-.19946844     , .69680979E-2  ,-.14623832E-3    ,&
     & .17229114E-5  ,-.83773053E-8  , 27.09         , 1.02986         ,&
     & .0            , .0            , .0            ,.0               ,&
     &2.4540781      , .9234824E2    , 1.            ,.00015           ,&
     & 1.350         ,  .0101325     , 2.6537        ,.2029E2          ,&
     & 24.54         , 44.4         ,600.            , 1.2             ,&
     & 100.          , 0.            , 2000.         /
!        CARBON MONOXIDE-----------HUST AND STEWART    1963
      DATA (COF(I,7),I=1,57)/        .29692807       , .36547544       ,&
     &-.80241533E2  ,-.16713841E5  , .13129095E6     , .0              ,&
     & .65308710    ,-.80131497E2  , .0              , .73365157       ,&
     & .0           , .0           , .0              , .41363247E3     ,&
     & .18061147E7  ,+.63943277E9  ,-.66472139E11    ,-.24767075E8     ,&
     & .96546154E10 ,-.47286423E12 , .74632444E1     , .58725319E9     ,&
     & .0           , .0           , 5*0.            , .13172848E2     ,&
     &-.52482838E3   ,-.18207888     , .20385068E-2  ,-.13699320E-4    ,&
     & .51448459E-7  ,-.81676932E-10 ,68.14          , 1.0392602       ,&
     &-.5033322E-5   , .26032523E-7  , .0            , .0              ,&
     &4.7928763      , .35387E3      , 1.            , .5E-4           ,&
     &.86            ,.0101325       ,3.4986         ,.2027E2          ,&
     & 68.14         , 132.91        , 600.          , .81             ,&
     & 180.          , .0            , 2000.         /
!        HELIUM--------------------MANN 1962
      DATA (COF(I,8),I=1,57)/         2.07722241    ,  8.44702446      ,&
     &-.234056435E3  , .478576622E2  ,-.119375479E3  , .0              ,&
     & .202301846E2  , .146535214E4  ,  .0           , .0              ,&
     &-.11186866E5   , .0            , .763244863E5  , -.70728125E5    ,&
     &-.790460718E4  ,.545307676E5   ,  .159411473E6 ,-.182597797E7    ,&
     &.414625688E7   , -.168878287E8 , .405704262E3  , .303082504E3    ,&
     & .0           , .0             , 5*0.          ,  -.21700826E1   ,&
     &-.32255748E1  , .23561801E1    ,-.95999492     ,   .22307327     ,&
     &-.26190460E-1 , .12190675E-2   , .42144E1      ,   .519306235E1  ,&
     & .0           , .0             , .0            , .0              ,&
     &4.5053       , 36.5582         , .1E1          , .4E-4           ,&
     &  .2100       ,  .010132       , .2274644      ,   .101325E2     ,&
     & .3E1         , .52E1          , .600E3        , .18             ,&
     & .35E2        , .0             , .25E4           /
!     FLUORINE  BARON (BENDER) 1/72 NBS DATA USED
      DATA ( COF(I,9),I=1,57) / .2188159E+0    ,     .21051539         ,&
     &-.81175243E+02   , .23064630E+04 ,-.63183150E+06,  .96908966E+07 ,&
     & .25933660E+00   ,-.63689419E+02 , .68240310E+04, -.13779899E+00 ,&
     & .74842154E+02   , .30516478E+00 ,-.13304083E+03,  .50315749E+02 ,&
     & .40392186E+06   ,-.27610459E+08 ,-.42709623E+10, -.33659737E+07 ,&
     & .65832302E+09   ,-.61907880E+09 ,3.06E+00      ,   8*0.         ,&
     & .51121868E+00   ,-.35538649E+03 , .16184765E+00, -.28441699E-02 ,&
     & .25041457E-04   ,-.11104809E-06 , .1987869E-09 ,   84.953       ,&
     &7.0981845E+00    ,-.37629101E-02 , .26888016E-04, -.38658913E-07 ,&
     & .17515391E-10   ,3.8416E+00     , 63.23E+00    ,  .1101145E+00  ,&
     &1.0E-06          ,2.25           ,  .011E+00   , 5.215E+00      , &
     & 20.3E+00        , 53.48E+00   , 144.31E+00   ,  500.E+00       , &
     & 1.9E+0          ,  170.E+00     ,  0.E+0       ,  450.E+01 /
!        PARAHYDROGEN -------------BARON (BENDER) 1971
      DATA (COF(I,10),I=1,57)/4.1260486, .44446150E+2,                  &
     &-.38659604E+4,-.10966550E+6, .12080022E+7, -.54747655E+7,         &
     &-.33278647E+3, .81345734E+5, .26294257E+6, .30063983E+5,          &
     &-.33024955E+7,-.24686707E+6, .47555234E+8, -.12064332E+9,         &
     &-.49289827E+8, .15925894E+10, -.87182365E+10, .66330266E+11,      &
     &.16366622E+11, -.12954419E+14, .1050E+4, 8*0.,                    &
     &-.10593817E+1,  -.35249570E+2,  .37870039E+0,  -.23601614E-1,     &
     &.85529568E-3,  -.16180625E-4,  .12500233E-6, 20.268,              &
     &14.7599360,-.21977388,.32100769E-2,-.12061502E-4,.57121808E-8,    &
     &21.17642,209.9406,1.0,                                            &
     &.1E-4,  .11,  .101325E-3,  1.2928,  100.,  13.8,  32.976,         &
     &3000., .09,  60.,  0.,  1000./
      DATA (AKSTCO(I,1),I=1,18)/3.,3.,-4.,0.,0.,.983,.3891166E-2,       &
     & -.96910013E-1,-.16491245E-1,-.84178805E-1,.85181928,.17517383E-1 &
     & ,.903089987 ,0.,0., .631,.136720567, 1.778 /
!
!        TRANSPORT PROPERTIES ARE LIMITED BY THE STATE OF THE ART
!           READ REPORT ON GASP FOR COMPLETE DETAILS
!
!
      DATA (AKSTCO(I,2),I=1,18)/3.,3.,-4.,0.,0., 1.0, 0.,               &
     &0.0,-1.12398985,.79352374,.80465682,-.21370471E-2,.39794,         &
     &0.,0., .729, .081707270, 1.0   /
      DATA (AKSTCO(I,3),I=1,8 )/1., 3., -.22,-.10587949, -.21522105,    &
     &1.31584644, .16196636E-2  ,1.0/
!   CH4
      DATA(AKSTCO(I,4),I=1,15)/ 2.,4.,-.280,0.0, .86187029,.26550054,   &
     & 1.0243779, -.3075470, .238, 1.6717516,-3.6259505,2.7281748,      &
     & .48625031, -.26386826, .71948992/
!     CH4
      DATA (TCCOF(I,1) , I=1,15 )  /.75226444,.20119504E-2, 190.77  ,   &
     &.1620  , 1. ,  4. ,-.9586073  ,.92583319, 1.55004971, .87042136 , &
     &1.29907897  ,-7.01464957, .617  , 4.,   1.E-4/
!     N2
      DATA (TCCOF(I,2) , I=1,15)   /.11425981E1,.20865873E-2, 126.3  ,  &
     &.3105   , 1.,  4., -1.39794 ,.51372790E-1,.17155947,.39312360  ,  &
     &1.56022809  ,-6.85395932, .39967,2., 2.77E-5   /
!     O2
      DATA (TCCOF(I,3) ,I=1,15)    /.96364877,.21320779E-2,  154.78    ,&
     &.4325  , 1.  , 4.  ,-2.6383, 0.,.58701320E-1,.23887321,1.39461538,&
     &-6.86562723 , .5132176,2.,3.40E-5/
!     AR
      DATA (TCCOF(I,4), I=1,15)    /.11037087E1,.21277034E-2,150.70,    &
     &.5310   ,1.,4., -1.744727, .47325357E-1, .29305409, .72199889,    &
     &1.56855570, -6.94782507, .45108,1.,2.24E-5/
!     CO2
      DATA (TCCOF(I,5),I=1,15)     /.98614941,.16279794E-2, 304.2,      &
     &.4640   ,1., 4., -4.699   , .25170328E-1, .28823104, .59310573,   &
     &1.5385090,-6.98446882,.4226,     3.,4.03E-5/
!      NE        USING K-K* FOR ARGON WITH NE XLAM AND ZC5
      DATA (TCCOF(I,6),I=1,15)/.95859  , .244676E-2  , 44.4, .483,      &
     &1.,4.,-1.744727,.47325357E-1,.29305409,.72199889,1.56855570,      &
     &-6.94782507,.45108,1.,2.6E-5/
!      CO      USE N2 WITH CO CONSTANTS
      DATA (TCCOF(I,7),I=1,15)/1.127527,  .226934E-2 ,132.91,           &
     &.2997,1.,4.,-1.39794,.51372790E-1,.17155947,.39312360,1.56022809, &
     &-6.85395932,.39967,2.,2.79E-5/
!      HE
      DATA (TCCOF(I,8),I=1,15)/1.536,.2579E-2,5.2,.0693,1.,4.,-1.744727,&
     &.47325357E-1,.29305409,.72199889,1.56855570,-6.94782507,.45108,1.,&
     &5.E-5/
!
!    F2     SETUP TO LOAD DIATOMIC CURVE FOR K* EVEN THOUGH IS NOT USED
      DATA (TCCOF(I,9),I=1,15)/ 1.02031, .19754E-2, 144.31, .57375,     &
     & 1.,4.,-1.301, .8474771,1.3197988, .47309461, 1.3532272,          &
     &-6.3527256, .4625, 2.,1./
!     DATA   NEEDED    FOR   H2  TRANSPORT  PROPERTIES
!   LOAD CRITICAL RHO AND T FOR THERM AND H2
      DATA (TCCOF(I,10),I=1,15)/2*0.,32.976,.03143,9*0.,1.0,0.0/
       ! item 14 changed tto 1.o by RLC   16Sept1997
      DATA DIST/3.822,3.681,3.499,3.421,3.952,2.82,3.59,2.551,3.357,    &
     &0.0/
      DATA EPSOK/.73E-2,1.093E-2,1.E-2,.837E-2,.5E-2,3.0488E-2,.9066E-2,&
     &9.7847E-2,.8881E-2,0.0/
      DATA WTMOL/16.04,28.016,31.9988,39.94,44.01,20.183,28.01,4.003,   &
     &39.98,0.0/
      DATA ZETA           /.46890513E-1,.40786245E-1,.30115154E-1,      &
     &.27628636E-1,.22407394E-1,.47495E-1,.402544E-1,.3843, .268525E-1, &
     & .40786245E-1/
      DATA FF/ 1.015, 3*1.0, 1.02, 1.1, 1., 2.27, 1.0, 1.0/
      DATA DTRIPL/.57,1.121,1.400,1.415,1.17,1.247,.836,.21,1.71,.11/
      DATA SWT/113.,2*68.,86.,216.,25.,70.,6.3,60.,68. /
      DATA KSWT /5*1,3*2,2*1/
      DATA DIFTT/7*1.,.1,2*1./
      DATA RHOSWT/2.5,2.2,6*2.4,2.2,2.2/

      END
!$IBFTC SETUPP
      SUBROUTINE  SETUP (NAMGAS)
!
!     ---------------------------VERSION 2/1/72-------------------------
!     ------THIS ROUTINE OVERLAYS COFFICIENTS FOR SPECIFIED FLUID-NAMGAS
!     ------ IN PROGRAM COMMON BLOCKS
!     ------COEFFICIENTS FOR ALL FLUIDS ARE PRESTORED IN BLOCK DATA
!
      CHARACTER*3 NAMGAS
      COMMON/ALLCOF/COF(57,10),TCCOF(15,10),AKSTCO(18,4),DIST(10),      &
     &WTMOL(10),EPSOK(10),ZETA(10),FF(10),SWT(10),KSWT(10),DIFTT(10),   &
     &RHOSWT(10),DTRIPL(10)
!       COMMON TO DETERMINE FLUID AND INFORM USER OF HIS CHOICE
      CHARACTER*3 MATCH(10)
      CHARACTER*90 MESSAG(10)
      COMMON/GASES/match,messag
!!!      COMMON/GASES/MATCH(10),MESSAG(15,10)
      COMMON/WHAT/KGAS
!          SPECIAL SWITCHES FOR HELIUM,FLUORINE,HYDROGEN
      COMMON/HELFLU/ IHE, IFL, IHY
!     COMMON FOR EQUATION OF STATE COEFFICIENTS
      COMMON/BEND/ R(1),CP1,CP2,CP3,CP4,CP5,CP6,CP7,CP8,CP9,CP10,       &
     & CP11,CP12,CP13,CP14,CP15,CP16,CP17,CP18,CP19,CP20,               &
     & CP21,CP22,CP23,CP24 ,CP25,CP26,CP27,CP28
!      SATURATION CURVE COEFFICIENTS
      COMMON/COSAT/ CPS1,CPS2,CPS3,CPS4,CPS5,CPS6,CPS7
!      CPO CURVE COEFFICIENTS
      COMMON/COCPO/TO(1),CCOP1,CCOP2,CCOP3,CCOP4,CCOP5
!      REFERENCE ENTROPY, ENTHALPY , CP COERRECTION FACTOR
      COMMON/REFNO/SOTO,HOTO,CPOCOR,HTERM,STERM
!      PARAMETERS FOR CHECKING REGION AND RANGE LIMITS ON DENSITY
!      PRESSURE, TEMPERATURE, ENTHALPY.  DENSITY AND TEMPERATURE ESTIMAT
!       NEWTON-RAPHSON ITERATIONS
      COMMON/CHECKS/DCH1(1),DCH2,PCH1,PCH2,PCH3,TCH1,TCH2,TCH3,DST,TST, &
     & HSCH1,HSCH2
!      CONSTANTS FOR THERMAL CONDUCTIVITY CALCULATION.
      COMMON/TCOND/CKMKST(9), XLAMB,ZC5, RHOCR, TCR, TCSTAR, CKSTAR(18)
!      CONSTANTS FOR VISCOSITY CALCULATION.
      COMMON/COFMU/EPSK,WM,DIS,RHOCRT,ZETAA
!      CONSTANTS FOR SURFACE TENSION CALCULATION.
      COMMON/SURCON/PCTC,TCRIT,FIXIT,ZET
!      SWITCH FOR CP-CV CALCULATION.  VALUES IN CERTAIN REGIONS ARE
!       CALCULATED BY NUMERICAL DIFFERENTIATION ...WHENEVER THE ANALYTIC
!       DERIVATIVES ARE GOOD, THEY ARE USED.........
      COMMON/SWITS/ KSWIT,TSWIT,DIFT,RSWIT
      COMMON/SOLID/DTRIP
      KGAS=0
      DO 10 I=1,10
      IF(NAMGAS.EQ.MATCH(I)) KGAS=I
   10 END DO
      IF(KGAS.EQ.0) GO TO 70
      IHE=0
      IF (KGAS.EQ.8) IHE=1
      IFL=0
      IF (KGAS.EQ.9) IFL=1
      IF (IFL.EQ.1) WRITE(6,72)
      IHY=0
      IF(KGAS.EQ.10) IHY=1
      WRITE(6,101) MESSAG(KGAS)
  101 FORMAT ('0',A)
!           STORE  CONSTANTS  FOR  BENDER-S EQ. OF  STATE
      DO 20 I=1,29
        R(I)= COF(I,KGAS)
   20 END DO
!           STORE  SATURATION  CURVE  FOR  GAS
!!!      DO 25 I=30,36              ! this is cute, but not good Fortran
!!!   25 CPS1(I-29) = COF(I,KGAS)
      CPS1=COF(30,KGAS)
      CPS2=COF(31,KGAS)
      CPS3=COF(32,KGAS)
      CPS4=COF(33,KGAS)
      CPS5=COF(34,KGAS)
      CPS6=COF(35,KGAS)
      CPS7=COF(36,KGAS)
!           STORE  CPO  COEFFICIENTS  AND  REFERENCE  TO
!!!      DO 30 I=37,42
!!!        TO(I-36) = COF(I,KGAS)
!!!   30 END DO
      TO(1)=COF(37,KGAS)
      CCOP1=COF(38,KGAS)
      CCOP2=COF(39,KGAS)
      CCOP3=COF(40,KGAS)
      CCOP4=COF(41,KGAS)
      CCOP5=COF(42,KGAS)
!           STORE  REFERENCE  ENTROPY,  ENTHALPY,  CPO  CORRECTION FACTO
      TX=TO(1)
!    COMPUTE LOWER BOUND FOR THE ENTHALPY AND ENTROPY CALCULATIONS IN HS
      A1=CCOP5/4.
      A2=CCOP4/3.
      A3=CCOP3/2.
      STERM =(((A1*TX+A2)*TX+A3)*TX+CCOP2)*TX+CCOP1*ALOG(TX)
      A1=CCOP5/5.
      A2=CCOP4/4.
      A3=CCOP3/3.
      A4=CCOP2/2.
      HTERM     =((((A1*TX+A2)*TX+A3)*TX+A4)*TX+CCOP1)*TX
      IF(IHY.EQ.1) STERM=4.968*4.184/2.01572*ALOG(TX)
      IF (IHY.EQ.1) HTERM=4.968*4.184/2.01572*TX
      SOTO= COF(43,KGAS)
      HOTO= COF(44,KGAS)
      CPOCOR=COF(45,KGAS)
!           STORE CRITICAL VALUES AND REGION BOUNDARY CONSTANTS
!!!      DO 35 I=46,57
!!!        DCH1(I-45) = COF(I,KGAS)
!!!   35 END DO
      DCH1(1)=COF(46,KGAS)
      DCH2=COF(47,KGAS)
      PCH1=COF(48,KGAS)
      PCH2=COF(49,KGAS)
      PCH3=COF(50,KGAS)
      TCH1=COF(51,KGAS)
      TCH2=COF(52,KGAS)
      TCH3=COF(53,KGAS)
      DST=COF(54,KGAS)
      TST=COF(55,KGAS)
      HSCH1=COF(56,KGAS)
      HSCH2=COF(57,KGAS)
      XLAMB= TCCOF(1,KGAS)
      ZC5= TCCOF(2,KGAS)
      TCR= TCCOF(3,KGAS)
      RHOCR = TCCOF(4,KGAS)
!           LOAD K-K*  CURVE COEFFICIENTS
      DO 45  I=5,13
        CKMKST(I-4)= TCCOF(I,KGAS)
   45 END DO
      TCSTAR= TCCOF(15,KGAS)
      I= TCCOF(14,KGAS)+.1
!                     LOAD  KSTAR/KSTARTC   CURVE  COEFFICIENTS
      DO 50  J=1,18
        CKSTAR(J)= AKSTCO(J,I)
   50 END DO
!      STORE CONSTANTS FOR U-U* AND VISCOSITY CALCULATION
       WM=WTMOL(KGAS)
       ZETAA=ZETA(KGAS)
       RHOCRT=TCCOF(4,KGAS)
       EPSK =EPSOK(KGAS)
       DIS=DIST(KGAS)
!       STORE CONSTANTS FOR SURFACE TENSION (  RF IS RIEDEL FACTOR)
      TCRIT=TCH2
      PCTC=(PCH2/.101325)**(2./3.)*TCH2**(1./3.)
      FIXIT=FF(KGAS)
      ZET=ZC5**(1./5.)
!  ENTER SWITCHING PARAMETER FOR CP,CV CALCULATION
      KSWIT=KSWT(KGAS)
      TSWIT=SWT(KGAS)
      IF (KSWIT.EQ.2) TSWIT=TCH2
      DIFT=DIFTT(KGAS)
      RSWIT=RHOSWT(KGAS) *RHOCRT
      DTRIP=DTRIPL(KGAS)
      RETURN
   70 WRITE(6,71)
   71 FORMAT(' ERROR IN CODE FOR NAMGAS - NO CONSTANTS STORED'/         &
     & ' PROGRAM  STOP.')
   72 FORMAT(' THE REGION 125-145K FOR T AND P GREATER THAN 10 ATM'/    &
     & ' YIELDS POOR RESULTS FOR THE DERIVED PROEPRTIES--BEWARE. ' )
      STOP
      END
!$IBFTC SCHEC
      FUNCTION   CHECK(KU,KR,T)
!
!     ---------------------------VERSION 2/1/72-------------------------
!
      COMMON/CON123/ DCONV(5),TCONV(5),PCONV(5)
      COMMON/IERROR/ IROUT
      COMMON/CHECKS/DCH1,DCH2,PCH1,PCH2,PCH3,TCH1,TCH2,TCH3,DST,TST,    &
     & HSCH1,HSCH2
!!!      DIMENSION   FM1(9), FM2(9), FM3(9), FMT(9), ROUT(13)
      CHARACTER*54 fm1,fm2,fm3,fmt
      CHARACTER*6 rout(13)
      DATA  FM1 /'(G12.4 ,''IS OUT OF RANGE FOR T  IN SUB.-'',A6 )'/
      DATA  FM2 /'(G12.4 ,''IS OUT OF RANGE FOR P  IN SUB.-'',A6 )'/
      DATA  FM3 /'(G12.4 ,''IS OUT OF RANGE FOR D  IN SUB.-'',A6 )'/
      DATA  ROUT /'DENS  ', 'PRESS ', 'TEMP  ', 'ENTH  ', 'ENT   ',     &
     &  'TEMPPH', 'TEMPPS', 'SONIC ', 'SPCHP ', 'SPCHV ', 'THERM ',     &
     &  'VISC  ', 'SURF  '   /
!
!        CONVERT TEMPERATURE T TO DEGREES KELVIN AND CHECK
!        FOR OUT OF RANGE.  UNITS ARE SPECIFIED BY KU.  IF KR
!        IS SPECIFIED AS 1, T IS CHECKED FOR OUT OF SATURATION
!        RANGE.
!
      ENTRY   TCHECK (KU,KR,T)
      CHECK=T/TCONV(KU)
      CH1=TCH1
      CH2=TCH2
      CH3=TCH3
!!!      DO 1 J=1,9
!!!    1 FMT(J)=FM1(J)
      fmt=fm1
      GO TO 10
!
!        CONVERT PRESSURE P TO ATMOSPHERES AND CHECK
!        FOR OUT OF RANGE.  UNITS ARE SPECIFIED BY KU. IF KR IS
!        SPECIFIED AS 1, P IS CHECKED FOR OUT OF SATURATION
      ENTRY  PCHECK(KU,KR,P)
      CHECK=P/PCONV(KU)
      CH1= PCH1
      CH2= PCH2
      CH3= PCH3
!!!      DO 2 J=1,9
!!!    2 FMT(J)=FM2(J)
      fmt=fm2
      GO TO 10
!
!        CONVERT DENSITY D TO GRAM-MOLES/LITER AND CHECK
!        FOR OUT OF RANGE.  UNITS ARE SPECIFIED BY KU.
!
      ENTRY DCHECK(KU,D)
      CHECK =D/DCONV(KU)
      CH1=DCH1
      CH3=DCH2
!!!      DO 3 J=1,9
      fmt=fm3
      GO TO 20
   10 IF(KR.EQ.1) GO TO 30
   20 IF(CHECK.LT.CH1) GO TO 40
      IF(CHECK.GT.CH3) GO TO 40
   25 RETURN
   30 IF(CHECK.LT.CH1) GO TO 40
      IF(CHECK.LE.CH2) GO TO 25
   40 WRITE(6,FMT) CHECK,ROUT(IROUT)
      RETURN
      END
!$IBFTC ROOT
      SUBROUTINE ROOT (X0,X2,FOFX,FUNC,X1)
!
!     ---------------------------VERSION 2/1/72-------------------------
!        SAME AS ROOTX - NEEDED  TO  PREVENT  RECURSION
!        SOLVE FOR X1 SUCH THAT FUNC(X1) = FOFX, WHERE X1 LIES
!        BETWEEN X0 AND X2
!
      real func
      external func
      INTEGER:: jump
      COMMON /CHECK1/KOUNT
      TOL=1.E-5
      XX0 = X0
      XX2 = X2
      F0 = FUNC(XX0)
      F2 = FUNC(XX2)
      A=(FOFX-F0)/(F2-F0)
      IF (A) 1007,120,120
  120 IF (A-1.) 130,130,1008
  130 IF (FOFX-0.) 80,70,80
   70 jump=1   !!! ASSIGN 100 TO JUMP
      GO TO 90
   80 jump=2   !!! ASSIGN 110 TO JUMP
   90 X = (XX0+XX2)/2.
      KOUNT = 0
  150 X1 = X
      KOUNT = KOUNT + 1
      A = FOFX - F2
      FX = FUNC(X)
      FXL=F0+(X-XX0)*(F2-F0)/(XX2-XX0)
      B=ABS((FX-FXL)/(F2-F0))
      IF (A*(FX-FOFX) .LT. 0.) GO TO 1001
      XX0 = X
      F0=FX
      IF (B-.3) 10,20,20
   20 X = (X+XX2)/2.
      GO TO 40
 1001 XX2 = X
      F2 = FX
      IF (B-.3) 10,30,30
   30 X = (XX0+X)/2.
      GO TO 40
   10 X=XX0+(FOFX-F0)*(XX2-XX0)/(F2-F0)
   40 IF (ABS((X-X1)/X)-TOL  ) 50,1000,1000
   50 CONTINUE   !!! GO TO JUMP,(100,110)
      IF (jump==1) GO TO 100
      IF (jump==2) GO TO 110

  100 IF (ABS(FUNC(X))-TOL*10. )60,1000,1000
  110 IF (ABS((FOFX-FUNC(X))/FOFX)-TOL  ) 60,1000,1000
 1000 IF (KOUNT.GT.40) TOL=TOL*10.
      IF (KOUNT.GT.60) TOL=TOL*10.
      IF (KOUNT.GT.80) TOL=TOL*10.
      IF (KOUNT.LT.100) GO TO 150
  160 WRITE (6,170) X1,X
  170 FORMAT ('AN ITERATION HAS BEEN TERMINATED AT 100 ITERATIONS.'/    &
     & ' THE LAST TWO VALUES WERE ', 3G15.5)
   60 X1=X
      RETURN
 1007 X1 = X0
      GO TO 140
 1008 X1 = X2
  140 WRITE(6,141)
  141 FORMAT(' SOLUTION OUT OF RANGE')
      RETURN
      END
!$IBFTC ROOTX
      SUBROUTINE ROOTX(X0,X2,FOFX,FUNC,X1)
!     ---------------------------VERSION 2/1/72-------------------------
!
!        SOLVE FOR X1 SUCH THAT FUNC(X1) = FOFX, WHERE X1 LIES
!        BETWEEN X0 AND X2
!
      real func
      INTEGER:: jump
      external func
      COMMON /CHECK2/KOUNT
      TOL=1.E-5
      XX0 = X0
      XX2 = X2
      F0 = FUNC(XX0)
      F2 = FUNC(XX2)
      A=(FOFX-F0)/(F2-F0)
      IF (A) 1007,120,120
  120 IF (A-1.) 130,130,1008
  130 IF (FOFX-0.) 80,70,80
   70 jump=1   !!! ASSIGN 100 TO JUMP
      GO TO 90
   80 jump=2   !!! ASSIGN 110 TO JUMP
   90 X = (XX0+XX2)/2.
      KOUNT = 0
  150 X1 = X
      KOUNT = KOUNT + 1
      A = FOFX - F2
      FX = FUNC(X)
      FXL=F0+(X-XX0)*(F2-F0)/(XX2-XX0)
      B=ABS((FX-FXL)/(F2-F0))
      IF (A*(FX-FOFX) .LT. 0.) GO TO 1001
      XX0 = X
      F0=FX
      IF (B-.3) 10,20,20
   20 X = (X+XX2)/2.
      GO TO 40
 1001 XX2 = X
      F2 = FX
      IF (B-.3) 10,30,30
   30 X = (XX0+X)/2.
      GO TO 40
   10 X=XX0+(FOFX-F0)*(XX2-XX0)/(F2-F0)
   40 IF (ABS((X-X1)/X)-TOL  ) 50,1000,1000
   50 CONTINUE   !!! GO TO JUMP,(100,110)
      IF (jump==1) GO TO 100
      IF (jump==2) GO TO 110
  100 IF (ABS(FUNC(X))-TOL*10. )60,1000,1000
  110 IF (ABS((FOFX-FUNC(X))/FOFX)-TOL  ) 60,1000,1000
 1000 IF (KOUNT.GT.40) TOL=TOL*10.
      IF (KOUNT.GT.60) TOL=TOL*10.
      IF (KOUNT.GT.80) TOL=TOL*10.
      IF (KOUNT.LT.100) GO TO 150
  160 WRITE (6,170) X1,X
  170 FORMAT (' AN ITERATION HAS BEEN TERMINATED AT 100 ITERATIONS.'/   &
     & ' THE LAST TWO VALUES WERE ', 3G15.5)
   60 X1=X
      RETURN
 1007 X1 = X0
      GO TO 140
 1008 X1 = X2
  140 WRITE(6,141)
  141 FORMAT(' SOLUTION OUT OF RANGE')
      RETURN
      END
!$IBFTC SPLIN
!     ---------------------------VERSION 2/1/72-------------------------
      SUBROUTINE SPLINA(X,Y,NX,T,NT,YINT,KFD,KERROR)
      DIMENSION X(*),Y(*),T(*),YINT(*),U(50),XM(50)
      EQUIVALENCE (XM(1), U(1))
      NXM1= NX-1
      NXP1 = NX+1
      H= X(2)-X(1)
      KERROR=0
      H2 = H*H
      YINT(1)= .5
      U(1)= 0.
      DO 10 I=2,NXM1
        P= .5*YINT(I-1)+2.
        YINT(I)= -.5/P
        D = 3.*((Y(I+1)-Y(I))-(Y(I)-Y(I-1)))/H2
        U(I)= (D-.5*U(I-1))/P
   10 END DO
      P =-.5*YINT(NXM1)+1.
      YINT(NX)= 0.
      U(NX)= .5*U(NXM1)/P
      XM(NX) = U(NX)
      DO 20 I=2,NX
        K= NXP1 -I
        XM(K)= YINT(K)*XM(K+1)+U(K)
   20 END DO
      DO 90 I=1,NT
      K=2
      IF(T(I)-X(1))60,85,70
   60 WRITE(6,600) T(1), X(1), X(NX)
  600 FORMAT(' AN ERROR OCCURRED DURING INTERPOLATION' /                &
     &  ' THE X WHERE INTERPOLATION DESIRED (XINT) IS OUT OF RANGE.'/   &
     &  ' THE ROUTINE WILL TRY TO EXTRAPOLATE.'/                        &
     &  ' XINT=', G15.7, ' X(1)=', G15.7, ' X(N)=', G15.7  )
      KERROR = KERROR+1
      GO TO 85
   70 IF(T(I)-X(K)) 85,85,80
   80 K= K+1
      IF(K-NX) 70,70,81
   81 KERROR = KERROR+1
      WRITE(6,600) T(I), X(1), X(NX)
      K=NX
   85 IF (KFD) 40,40,50
   40 YINT(I)=(XM(K-1)*(X(K)-T(I))**3+XM(K)*(T(I)-X(K-1))**3+(6.*Y(K-1)-&
     & XM(K-1)*H2)*(X(K)-T(I))+(6.*Y(K)-XM(K)*H2)*(T(I)-X(K-1)))/(6.*H)
      GO TO 91
   50 YINT(I)=(-XM(K-1)*(X(K)-T(I))**2/2.+XM(K)*(T(I)-X(K-1))**2/2.+    &
     & Y(K)-Y(K-1)-(XM(K)-XM(K-1))*H2/6.)/H
   91 CONTINUE
   90 END DO
      RETURN
      END
!$IBFTC POLY
      FUNCTION POLY(X,COEF)
!     ---------------------------VERSION 2/1/72-------------------------
!
!        EVALUATE THE POLYNOMIAL IN X DESCRIBED BY COEF
!
      DIMENSION COEF(*)
      XS=X
      NRANGE=COEF(1)+.1
      NDEG=COEF(2)+.1
      ISTEP=NDEG+2
      ITEST=3+NRANGE*ISTEP
      IF (XS-COEF(3)) 10,20,20
   20 IF (XS-COEF(ITEST)) 30,30,40
   30 IBEG=3+ISTEP
      IEND=ITEST-ISTEP
      DO 50 I=IBEG,IEND,ISTEP
        IF (xs <= coef(i)) GO TO 60   !!! IF (XS-COEF(I)) 60,60,50
   50 END DO
      I=ITEST
      GO TO 60
   10 I=3+ISTEP
      GO TO 70
   40 I=ITEST
   70 WRITE(6,71) XS
   71 FORMAT(E15.5,                                                     &
     & ' IS OUT OF RANGE OF CURVE FIT  VALUE IS EXTRAPOLATED.')
   60 IBEG=I-ISTEP+2
      IEND=IBEG+NDEG-1
      POLY=COEF(IBEG-1)
      DO 80 I=IBEG,IEND
        POLY=POLY*XS+COEF(I)
   80 END DO
      RETURN
      END
!$IBFTC SOLVE
      FUNCTION SOLVE(XI,F,DF)
!     ---------------------------VERSION 2/1/72-------------------------
!
!        NEWTON-RAPHSON ITERATION GIVEN AN INITIAL ESTIMATE XI
!        AND THE FUNCTIONS F AND DF
      real f,df
      external f,df
      COMMON/CHECK1/NI
      TOL=1.E-5
      NI=0
      XO=XI
      XN=XI
   10 XOO=XO
      XO=XN
      XN=XO-F(XO)/DF(XO)
   12 NI=NI+1
      IF (ABS((XN-XO)/XN)-TOL  ) 70,20,20
   20 IF (NI.GT.40) TOL=TOL*10.
      IF (NI.GT.60) TOL=TOL*10.
      IF (NI.GT.80) TOL=TOL*10.
      IF (NI-100) 30,50,50
   30 IF (ABS((XN-XOO)/XN)-TOL  ) 40,10,10
   40 XN=(XO+XN)/2.
      GO TO 10
   50 WRITE (6,60) XOO,XO,XN
   60 FORMAT (' AN ITERATION HAS BEEN TERMINATED AT 100 ITERATIONS.'/   &
     & ' THE LAST THREE VALUES WERE ', 3G15.5)
   70 SOLVE=XN
      IF (XN.GE.0.) RETURN
      SOLVE=XI
      WRITE(6,75) XI
   75 FORMAT(' N-R ITERATION DID NOT FIND A VALID ANSWER.'/             &
     & '  INITIAL ESTIMATE RETURNED AS ANSWER. ', E12.6)
      RETURN
      END
!$IBFTC GUES
      SUBROUTINE DGUESS(TS,TCR,DST)
!    ---------------------------VERSION 2/1/72
      COMMON/COFMU/EPSK,WM,DIS,RHOCRT,ZETAA
      IF (TS/TCR.GT..9) DST=RHOCRT*1.75
      IF (TS/TCR.GT..96)DST=RHOCRT*1.60
      IF (TS/TCR.GT..98)DST=RHOCRT*1.50
      IF (TS/TCR.GT..993)DST = RHOCRT*1.4
      IF (TS/TCR.GT..998)DST = RHOCRT*1.30
      IF (TS/TCR.GT..999)  DST=RHOCRT*1.20
      IF (TS/TCR.GT..9995) DST=RHOCRT*1.10
      IF (TS/TCR.GT..9999) DST=RHOCRT*1.05
      RETURN
      END
!$IBFTC DENS1
      SUBROUTINE DENS(KU,T,P,D,DL,DV,KR)
!     ---------------------------VERSION 2/1/72-------------------------
!        COMPUTE DENSITY D GIVEN TEMPERATURE T AND PRESSURE P.
!        UNITS ARE SPECIFIED BY KU.  IF KR IS RETURNED OR
!        SPECIFIED AS 1, THE SATURATED LIQUID AND VAPOR DENSITIES,
!        DL AND DV RESPECTIVELY, ARE COMPUTED AS A FUNCTION
!        OF T OR P.  THE OTHER VALUE MUST BE INPUT AS 0.0 .
!
      COMMON/ CON123/ DCONV(5),TCONV(5),PCONV(5)
      COMMON /CHECK1/NI
      COMMON/BEND/ R,CP1,CP2,CP3,CP4,CP5,CP6,CP7,CP8,CP9,CP10,          &
     & CP11,CP12,CP13,CP14,CP15,CP16,CP17,CP18,CP19,CP20,               &
     & CP21,CP22,CP23,CP24,CP25,CP26,CP27,CP28
      COMMON/BEND13/A1,A4,A5,A6,A7,A8,A9,A99,PS
      COMMON/BEND15/A10,A11,A12,A13,A14
      COMMON/BEND17/TS
      COMMON/CHECKS/DCH1,DCH2,PCH1,PCH2,PCH3,TCH1,TCH2,TCH3,DST,TST,    &
     & HSCH1,HSCH2
      COMMON/COFMU/EPSK,WM,DIS,RHOCRT,ZETAA
      COMMON/HELFLU/ IHE, IFL, IHY
      COMMON/SOLID/DTRIP
      COMMON /IERROR/IROUT
      REAL dsf,ddsf
      EXTERNAL DSF,DDSF
      IROUT=1
      IF (KR.EQ.1) GO TO 70
      TS=TCHECK(KU,KR,T)
      GO TO 5
   70 IF (T.GT.0.0) GO TO 75
      PS=PCHECK(KU,KR,P)
      TS=TSS(PS)
      IF (T.LE.0.) T=TS*TCONV(KU)
      GO TO  5
   75 TS=TCHECK(KU,KR,T)
      CALL  PSSS(PS)
      IF (P.LE.0.) P=PS*PCONV(KU)
    5 T2=1.0
      IF (IHE.EQ.1) T2=TS
!
!   COMPUTATIONS COMMON TO ALL REGIONS
!
      A1= R*TS
      A2= 1./TS
      A3= A2*A2
      A4=((((CP25*A3 + CP21)*A2+CP5)*A2+CP4)*A2+CP3)*A2+CP2+CP1*TS
      A5=(CP23*A2+CP8)*A2+CP7+(CP22*TS+CP6)*TS
      A6= CP9*TS+CP10
      A7= CP11*TS+CP12 +CP24*A2
      A8=((CP16*A2+CP15)*A2+CP14)*A3
      A9=((CP19*A2+CP18)*A2+CP17)*A3
      A99=((CP28*A2+CP27)*A2+CP26)*A3
      A10=-2.*CP20/T2
      A11=6.*CP13
      A12=5.*A7
      A13= 4.*A6
      A14= 2.*A4
!
!        DETERMINE REGION
!
      IF (KR-1) 10,80,10
   10 PS=PCHECK(KU,KR,P)
      IF (PS-PCH2)110,110,100
  100 IF (TS-TCH2)120,120,130
  120 KR=2
      EST=RHOCRT*.90
      CALL ROOT(DCH2,EST,0.,DSF,DS)
      GO TO 150
  130 KR=3
      EST=RHOCRT*3.
      CALL ROOT(EST,DCH1,0.,DSF,DS)
      GO TO 150
  110 IF (TS-TCH2) 20,20,50
   20 CALL PSSS(PSS)
      IF (ABS((PSS-PS)/PSS)-1.E-4) 60,30,30
   30 IF (PS-PSS) 50,60,40
   40 KR=2
      DS=DST
      IF (TS/TCH2.LT..90.OR.PS/PSS.LT..95*PSS) GO TO 90
      CALL DGUESS(TS,TCH2,DS)
      GO TO 90
   50 KR=3
      DS=PS/(R*TS)/2.
      GO TO 90
!
!        REGION 1
!
   60 KR=1
   80 CONTINUE
      DS=DST
      IF (TCH2-TS.LT.1.E-4) GO TO 84
      IF (TCH2-TS.LT..10) GO TO 81
      CALL DGUESS(TS,TCH2,DS)
      DSL=SOLVE(DS, DSF,DDSF)
      DS=PS/(R*TS)/2.
      IF (TS/TCH2.GT..985)DS=.65*RHOCRT
      IF (TS/TCH2.GT..995) DS=.75*RHOCRT
      IF (TS/TCH2.GT..999) DS=.85*RHOCRT
      IF (TS/TCH2.GT..9995) DS=.90*RHOCRT
      DSV=SOLVE(DS,DSF,DDSF)
      DL=DSL*DCONV(KU)
      DV=DSV*DCONV(KU)
      RETURN
   81 DL=RHOCRT*DCONV(KU)*1.001
      DV=RHOCRT*DCONV(KU)*.999
!
!      THIS FIX IS APPLIED IN THE NEAR SUBCRITICAL AREA TO COMPENSATE FO
!     THE FAILURE OF NEWTON-RAPHSON TO CONVERGE ON A DENSITY WITHOUT A V
!     ACCURATE ESTIMATE.  PERHAPS ROOTX SHOULD BE USED IN THE REGION
!       T/TCRT .GT..99 TO PREVENT THIS PROBLEM
!      IF YOU WANT AN ERROR MESSAGE WHEN THIS HAPPENS REMOVE THE C-S FRO
!      FOLLOWING THREE CARDS.
!     WRITE(6,82) TS
!  82 FORMAT(' SATURATION CALL FOR TS= ', F10.3,
!     1 ' IS TOO NEAR CRITICAL.  DL=1.001*RHOCRT  DV=.999*RHOCRT')
      RETURN
   84 DL=RHOCRT*DCONV(KU)
      DV=DL
      RETURN
!
!        REGIONS 2 AND 3
!
   90 DS=SOLVE(DS,DSF,DDSF)
  150 D=DS*DCONV(KU)
      IF (DS.GT.DTRIP) WRITE(6,152) T,P
  152 FORMAT(' DENSITY SOLUTION IS BEYOND THE GRASP OF GASP.'/          &
     & ' ANSWER IS ON SOLID SIDE OF MELTING LOCUS FOR T =', F9.0,       &
     & ' P = ', F10.3, ' MN.')
      RETURN
      END
!$IBFTC DSUB1
      FUNCTION DSF(DS)
!
!     ---------------------------VERSION 2/1/72-------------------------
!        FUNCTION USED TO SOLVE FOR DENSITY DS GIVEN TEMPERATURE
!        AND PRESSURE
!
      COMMON/BEND/ R,CP1,CP2,CP3,CP4,CP5,CP6,CP7,CP8,CP9,CP10,          &
     & CP11,CP12,CP13,CP14,CP15,CP16,CP17,CP18,CP19,CP20,               &
     & CP21,CP22,CP23,CP24 ,CP25,CP26,CP27,CP28
      COMMON/BEND13/A1,A4,A5,A6,A7,A8,A9,A99,PS
      COMMON/BEND15/A10,A11,A12,A13,A14
      COMMON/HELFLU/ IHE, IFL, IHY
      COMMON/BEND17/TS
      T2=1.0
      IF (IHE.EQ.1) T2=TS
      B1=EXP(-CP20*DS*DS/T2)
      B2= A8*B1
      B3= A9*B1
      B4=A99*B1
      DSF=((((((B4*DS+CP13)*DS+A7+B3)*DS+A6)*DS+A5+B2)*DS+A4)*DS        &
     & +A1)*DS-PS
      RETURN
      ENTRY DDSF(DS)
!
!        DERIVATIVE OF FUNCTION USED TO SOLVE FOR DENSITY DS GIVEN
!        TEMPERATURE AND PRESSURE
!
      DDSF=((((((A10*(B4*DS*DS+B3)+7.*B4)*DS+A11)*DS+A12+(5.*A9+A10*A8)*&
     &B1)*DS+A13)*DS+3.*(A5+A8*B1))*DS+A14)*DS+A1
      DSF=DDSF
      RETURN
      END
!$IBFTC PRES1
      SUBROUTINE PRESS(KU,T,D,P,KR)
!
!     ---------------------------VERSION 2/1/72-------------------------
!        COMPUTE PRESSURE P GIVEN TEMPERATURE T AND DENSITY D.
!        UNITS ARE SPECIFIED BY KU.  IF KR IS RETURNED OR
!        SPECIFIED AS 1, P IS COMPUTED AT SATURATION AS A
!        FUNCTION OF T ONLY.
!
      COMMON/ CON123/ DCONV(5),TCONV(5),PCONV(5)
      COMMON/BEND/ R,CP1,CP2,CP3,CP4,CP5,CP6,CP7,CP8,CP9,CP10,          &
     & CP11,CP12,CP13,CP14,CP15,CP16,CP17,CP18,CP19,CP20,               &
     & CP21,CP22,CP23,CP24 ,CP25,CP26,CP27,CP28
      COMMON/CHECKS/DCH1,DCH2,PCH1,PCH2,PCH3,TCH1,TCH2,TCH3,DST,TST,    &
     &HSCH1,HSCH2
      COMMON/BEND17/TS
      COMMON/HELFLU/ IHE, IFL, IHY
      COMMON /IERROR/IROUT
      IROUT=2
!
!        DETERMINE REGION
!
      TS=TCHECK(KU,KR,T)
      T2=1.0
      IF (IHE.EQ.1) T2=TS
      IF (KR-1) 10,70,10
   10 DS=DCHECK(KU,D)
      IF (TS-TCH2) 20,20,50
   20 CALL DENS(1,TS,ZE,ZE,DSL,DSV,1)
      IF (DS-DSL) 30,60,40
   30 IF (DS-DSV) 50,60,60
   40 KR=2
      GO TO 80
   50 KR=3
      GO TO 80
!
!        REGION 1
!
   60 KR=1
   70 CALL PSSS(PS)
      GO TO 90
!
!        REGIONS 2 AND 3
!
   80 A1= R*TS
      A2 = 1./TS
      A3 = A2*A2
      A4=((((CP25*A3+CP21)*A2+CP5)*A2+CP4)*A2+CP3)*A2+CP2+CP1*TS
      A5=(CP23*A2+CP8)*A2+CP7+(CP22*TS+CP6)*TS
      A6 = CP9*TS+CP10
      A7= CP11*TS+CP12 +CP24*A2
      A8 = ((CP16*A2+CP15)*A2+CP14)*A3
      A9 = ((CP19*A2+CP18)*A2+CP17)*A3
      B1=EXP(-CP20*DS*DS/T2)
      A10=((CP28*A2+CP27)*A2+CP26)*A3*B1
      PS=((((((A10*DS+CP13)*DS+A7+A9*B1)*DS+A6)*DS+A5+A8*B1)*DS+A4)*DS+ &
     &A1)*DS
   90 P=PS*PCONV(KU)
      RETURN
      END
!$IBFTC TEMP1
      SUBROUTINE TEMP(KU,P,D,T,KR)
!
!     ---------------------------VERSION 2/1/72-------------------------
!        COMPUTE TEMPERATURE T GIVEN PRESSURE P AND DENSITY D.
!        UNITS ARE SPECIFIED BY KU.  IF KR IS RETURNED OR
!
      COMMON/ CON123/ DCONV(5),TCONV(5),PCONV(5)
      COMMON/BEND/ R,CP1,CP2,CP3,CP4,CP5,CP6,CP7,CP8,CP9,CP10,          &
     & CP11,CP12,CP13,CP14,CP15,CP16,CP17,CP18,CP19,CP20,               &
     & CP21,CP22,CP23,CP24 ,CP25,CP26,CP27,CP28
      COMMON/HELFLU/ IHE, IFL, IHY
      COMMON/CHECKS/DCH1,DCH2,PCH1,PCH2,PCH3,TCH1,TCH2,TCH3,DST,TST,    &
     &HSCH1,HSCH2
      COMMON/BEND3/DS,A1,A3,A4,A5,A6B,A7B,A8B,A12,G1
      COMMON /IERROR/IROUT
      REAL tsf,dtsf
      EXTERNAL TSF,DTSF
      IROUT=3
      PS=PCHECK(KU,KR,P)
!
!        DETERMINE REGION
!
      IF (KR-1) 10,70,10
   10 DS=DCHECK(KU,D)
      IF (PS-PCH2) 20,20,50
   20 TS=TSS(PS)
      CALL DENS(1,TS,ZE,ZE,DSL,DSV,1)
      IF (DS-DSL) 30,60,40
   30 IF (DS-DSV) 50,60,60
   40 KR=2
      TS=TS-10.
      IF (IHE.EQ.1) TS=3.0
      GO TO 80
   50 KR=3
      TS=TST
      GO TO 80
!
!        REGION 1
!
   60 KR=1
      GO TO 110
   70 TS=TSS(PS)
      GO TO 110
!
!        REGIONS 2 AND 3
!
   80 A1= DS*DS
      A3 =((((CP13*DS+CP12)*DS+CP10)*DS+CP7)*DS+CP2)*A1 -PS
      A4 = ((((CP11*DS+CP9)*DS+CP6)*DS+CP1)*DS+R)*DS
      A5=((CP24*A1+CP8)*DS+CP3)*A1
      A6B= (CP23*DS+CP4)*A1
      A7B=CP5*A1
      A8B=CP21*A1
      A12 =2.*CP22*A1*DS
      G1=CP25*A1
      TS=SOLVE(TS,TSF,DTSF)
!
!        VERIFY REGION
!
      IF (PS-PCH2)110,110,90
   90 IF (TS-TCH2) 100,100,110
  100 KR=2
  110 T=TS*TCONV(KU)
      RETURN
      END
!$IBFTC TSUB1
      FUNCTION TSF(TS)
!
!     ---------------------------VERSION 2/1/72-------------------------
!        FUNCTION USED TO SOLVE FOR TEMPERATURE TS GIVEN PRESSURE
!        AND DENSITY
!
      COMMON/BEND/ R,CP1,CP2,CP3,CP4,CP5,CP6,CP7,CP8,CP9,CP10,          &
     & CP11,CP12,CP13,CP14,CP15,CP16,CP17,CP18,CP19,CP20,               &
     & CP21,CP22,CP23,CP24,CP25,CP26,CP27,CP28
      COMMON/HELFLU/ IHE, IFL, IHY
      COMMON/BEND3/DS,A1,A3,A4,A5,A6B,A7B,A8B,A12,G1
      T2=1.0
      IF (IHE.EQ.1) T2=TS
      A2=DS*EXP(-CP20*DS*DS/T2)
      A6A=((CP26*A1+CP17)*A1+CP14)*A1*A2
      A7A=((CP27*A1+CP18)*A1+CP15)*A1*A2
      A8A=((CP28*A1+CP19)*A1+CP16)*A1*A2
      A6=A6A+A6B
      A7=A7A+A7B
      A8=A8A+A8B
      A9 = 2.*A6
      A10= 3.*A7
      A11= 4.*A8
      B1=1./TS
      B2=B1*B1
      A13=0.0
      IF (IHE.EQ.1) A13=(A6A+A7A*B1+A8A*B2)*CP20*A1*B2*B2
      TSF=((((G1*B2+A8)*B1+A7)*B1+A6)*B1+A5)*B1+(CP22*TS*A1*DS+A4)*TS+A3
      RETURN
      ENTRY DTSF(TS)
!
!        DERIVATIVE OF FUNCTION USED TO SOLVE FOR TEMPERATURE TS
!        GIVEN PRESSURE AND DENSITY
!
      DTSF=A4+A12*TS-((((6.*G1*B1*B1+A11)*B1+A10)*B1+A9)*B1+A5)*B1*B1   &
     &+A13
      TSF=DTSF
      RETURN
      END
!$IBFTC TSSS1
      FUNCTION TSSF(TSS)
!
!        FUNCTION USED TO SOLVE FOR SATURATION TEMPERATURE TSS
!        GIVEN PRESSURE
!
      COMMON/COSAT/ CPS1,CPS2,CPS3,CPS4,CPS5,CPS6,CPS7
      COMMON/BEND9/A1,A2,A3,A4,A5
      TSSF=((((CPS7*TSS+CPS6)*TSS+CPS5)*TSS+CPS4)*TSS+CPS3)*TSS+CPS2/   &
     &  TSS+A1
      RETURN
      ENTRY  DTSSF(TSS)
!
!        DERIVATIVE OF FUNCTION USED TO SOLVE FOR SATURATION
!        TEMPERATURE TSS GIVEN PRESSURE
!
      DTSSF=(((A2*TSS+A3)*TSS+A4)*TSS+A5)*TSS+CPS3-CPS2/(TSS*TSS)
      TSSF=DTSSF
      RETURN
      END
!$IBFTC TSST
      FUNCTION TSS(PS)
!
!     ---------------------------VERSION 2/1/72-------------------------
!        COMPUTE SATURATION TEMPERATURE GIVEN PRESSURE PS
!
      COMMON/CHECKS/DCH1(1),DCH2,PCH1,PCH2,PCH3,TCH1,TCH2,TCH3,DST,TST, &
     & HSCH1,HSCH2
      COMMON/COSAT/ CPS1   ,CPS2,CPS3,CPS4,CPS5,CPS6,CPS7
      COMMON/BEND9/A1,A2,A3,A4,A5
      EXTERNAL TSSF,DTSSF
      PS1=PS/.101325
      A1=CPS1-ALOG10(PS1)
      A2=5.*CPS7
      A3=4.*CPS6
      A4=3.*CPS5
      A5=2.*CPS4
      TESTM=PS/(PCH1-PCH2)*(TCH2-TCH1)+TCH1 -10.
      IF (TESTM.LT.TCH1) TESTM=TCH1
      TSS=SOLVE(TESTM,TSSF,DTSSF)
      RETURN
      END
!$IBFTC PSSS
      SUBROUTINE PSSS(PSS)
!     ---------------------------VERSION 2/1/72-------------------------
!
!        COMPUTE SATURATION PRESSURE PSS GIVEN TEMPERATURE
!
      COMMON/COSAT/ CPS1,CPS2,CPS3,CPS4,CPS5,CPS6,CPS7
      COMMON/BEND17/TS
      PSS=10.**(((((CPS7*TS+CPS6)*TS+CPS5)*TS+CPS4)*TS+CPS3)*TS+CPS2/TS+&
     & CPS1)
      PSS=PSS*.101325
      RETURN
      END
!$IBFTC ENTH1
      SUBROUTINE ENTH(KU,KR,T,P,D,H,HL,HV)
!
!     ---------------------------VERSION 2/1/72-------------------------
!        COMPUTE ENTHALPY H GIVEN TEMPERATURE T, PRESSURE P, AND
!        DENSITY D.  UNITS ARE SPECIFIED BY KU.  REGION IS
!        SPECIFIED BY KR.  IF KR IS SPECIFIED AS 1, THE SATURATED
!        LIQUID AND VAPOR ENTHALPIES, HL AND HV RESPECTIVELY,
!        ARE COMPUTED AS A FUNCTION OF T ONLY.
!
      COMMON/CONV6/HCONV(5)
      COMMON/BEND17/ TS
      COMMON/BEND23/ DSL,DSV
      COMMON/BEND29/ HSL,HSV
      COMMON/IERROR/ IROUT
      IROUT=4
      TS=TCHECK(KU,KR,T)
      PS=PCHECK(KU,KR,P)
      IF (KR.NE.1) DS=DCHECK(KU,D)
      GO TO (10,20,30),KR
!
!        REGION 1
!
   10 CALL HSLV(PS)
      HL=HSL*HCONV(KU)
      HV=HSV*HCONV(KU)
      RETURN
!
!        REGION 2
!
   20 CALL PSSS(PSS)
      CALL HSLV(PSS)
      HS = HSL+ PS/DS-PSS/DSL+HDINT(DS,DSL)
      GO TO 40
!
!        REGION 3
!
   30 HS=HSS(PS,DS)
   40 H=HS*HCONV(KU)
      RETURN
      END
!$IBFTC HSSLVF
      FUNCTION HSSLVF(PS)
!
!     ---------------------------VERSION 2/1/72-------------------------
!        FUNCTION USED TO COMPUTE SATURATED LIQUID ENTHALPY
!        FROM SATURATED VAPOR ENTHALPY OR SATURATED LIQUID ENTROPY
!        FROM SATURATED VAPOR ENTROPY GIVEN PRESSUE PS AND
!        TEMPERATURE
!
      COMMON/COSAT/ CPS1,CPS2,CPS3,CPS4,CPS5,CPS6,CPS7
      COMMON/BEND17/TS
      COMMON/BEND23/DSL,DSV
      DPDT=2.30258509*PS*((((5.*CPS7*TS+4.*CPS6)*TS+3.*CPS5)*TS+        &
     &2.*CPS4)*TS+CPS3-CPS2/(TS*TS))
      HSSLVF=DPDT*(1./DSV-1./DSL)
      RETURN
      END
!$IBFTC ENT1
      SUBROUTINE ENT(KU,KR,T,P,D,S,SL,SV)
!
!     ---------------------------VERSION 2/1/72-------------------------
!        COMPUTE ENTROPY S GIVEN TEMPERATURE T, PRESSURE P, AND
!        DENSITY D.  UNITS ARE SPECIFIED BY KU.  REGION IS
!        SPECIFIED BY KR.  IF KR IS SPECIFIED AS 1, THE SATURATED
!        LIQUID AND VAPOR ENTROPIES, SL AND SV RESPECTIVELY,
!        ARE COMPUTED AS A FUNCTION OF T ONLY.
!
      COMMON/BEND/ R,CP1,CP2,CP3,CP4,CP5,CP6,CP7,CP8,CP9,CP10,          &
     & CP11,CP12,CP13,CP14,CP15,CP16,CP17,CP18,CP19,CP20,               &
     & CP21,CP22,CP23,CP24,CP25,CP26,CP27,CP28
      COMMON/BEND21/SSL,SSV
      COMMON/BEND23/DSL,DSV
      COMMON/BEND17/TS
      COMMON/CONV4/SCONV(5)
      COMMON/IERROR/IROUT
      IROUT=5
      TS=TCHECK(KU,KR,T)
      PS=PCHECK(KU,KR,P)
      IF (KR.NE.1) DS=DCHECK(KU,D)
      GO TO (10,20,30),KR
!
!        REGION 1
!
   10 CALL SSLV(PS)
      SL=SSL*SCONV(KU)
      SV=SSV*SCONV(KU)
      RETURN
!
!        REGION 2
!
   20 CALL PSSS(PSS)
      CALL SSLV(PSS)
      SS = SSL + R*(ALOG(DSL)-ALOG(DS)) + SDINT(DS,DSL)
      GO TO 40
!
!        REGION 3
!
   30 SS=SSS(PS,DS)
   40 S=SS*SCONV(KU)
      RETURN
      END
!$IBFTC HDINT1
      FUNCTION HDINT(DS,DSL)
!
!     ---------------------------VERSION 2/1/72-------------------------
!        COMPUTE THE INTEGRAL IN THE ENTHALPY COMPUTATION FROM
!        DENSITY DSL TO DENSITY DS
!
      COMMON/BEND/ R,CP1,CP2,CP3,CP4,CP5,CP6,CP7,CP8,CP9,CP10,          &
     & CP11,CP12,CP13,CP14,CP15,CP16,CP17,CP18,CP19,CP20,               &
     & CP21,CP22,CP23,CP24,CP25,CP26,CP27,CP28
      COMMON/BEND17/TS
      COMMON/BEND27/A3,A4,A5,A6,A7,A8,A9,A10,A11,A12,A13,A66
      COMMON/HELFLU/ IHE, IFL, IHY
      T2=1.0
      IF (IHE.EQ.1) T2=TS
      A1= 1./TS
      A2= A1*A1
      A3=-A2*T2/(2.*CP20)
      A4=((((7.*CP25*A2+5.*CP21)*A1+4.*CP5)*A1+3.*CP4)*A1+2.*CP3)*A1+CP2
      A5=(-CP22/A2+CP7+2.*CP8*A1+3.*CP23*A2)/2.
      A6 = ((5.*CP16*A1+ 4.*CP15) * A1 + 3.*CP14)
      A66=((5.*CP28*A1+4.*CP27)*A1+3.*CP26)
      A7 = ((5.*CP19*A1 + 4.*CP18)*A1 + 3.*CP17)
      A8 = CP10/3.
      A9= (CP12+2.*CP24*A1)/4.
      A10 = CP13/5.
      A11=T2/CP20
      A12= ((CP16*A1+CP15)*A1+CP14)*A2
      A13= ((CP19*A1+CP18)*A1+CP17)*A2
      HDINT=HDINTF(DS)-HDINTF(DSL)
      RETURN
      ENTRY    SDINT(DS,DSL)
!
!        COMPUTE THE INTEGRAL IN THE ENTROPY COMPUTATION FROM
!        DENSITY DSL TO DENSITY DS
!
      T2=1.0
      IF (IHE.EQ.1) T2=TS
      A1=1./TS
      A2=A1*A1
      A3=((((6.*CP25*A2+4.*CP21)*A1+3.*CP5)*A1+2.*CP4)*A1+CP3)*A2-CP1
      A4=(-2.*CP22*TS-CP6+(2.*CP23*A1+CP8)*A2)/2.
      A5 = ((4.*CP16*A1 + 3.*CP15)*A1 + 2.*CP14)
      A6 = ((4.0*CP19*A1 + 3.0*CP18)*A1 + 2.0*CP17)
      A66=((4.0*CP28*A1+3.0*CP27)*A1+2.0*CP26)
      A7 = -A2/(2.*CP20)*A1*T2
      A8 =-CP9/3.
      A9=(-CP11+CP24*A2)/4.
      A10=T2/CP20
      A12= ((CP16*A1+CP15)*A1+CP14)*A2
      A13= ((CP19*A1+CP18)*A1+CP17)*A2
      SDINT=SDINTF(DS)-SDINTF(DSL)
      HDINT=SDINT
      RETURN
      END
!$IBFTC HDINT2
      FUNCTION HDINTF(DS)
!
!     ---------------------------VERSION 2/1/72-------------------------
!        FUNCTION USED TO COMPUTE THE INTEGRAL IN THE ENTHALPY
!        COMPUTATION BASED ON DENSITY DS
!
      COMMON/BEND/ R,CP1,CP2,CP3,CP4,CP5,CP6,CP7,CP8,CP9,CP10,          &
     & CP11,CP12,CP13,CP14,CP15,CP16,CP17,CP18,CP19,CP20,               &
     & CP21,CP22,CP23,CP24,CP25,CP26,CP27,CP28
      COMMON/HELFLU/ IHE, IFL, IHY
      COMMON/BEND17/TS
      COMMON/BEND27/A3,A4,A5,A6,A7,A8,A9,A10,A11,A12,A13,A66
      T2=1.0
      IF (IHE.EQ.1) T2=TS
      B1=DS*DS
      B2 = EXP(-CP20*B1/T2)
      A15=0.0
      IF (IHE.EQ.1) A15=(A12*(B1+A11)+A13*((B1+2.*A11)*B1+2.*A11*A11))  &
     &*B2/2.
      HDINTF= ((((A10*DS + A9 )*DS+A8)*DS+A5)*DS+A4)*DS+                &
     & A3*B2*(A6 +(B1+A11)*A7+A66*(B1*B1+2.*B1/CP20+2./(CP20*CP20)))+A15
      RETURN
      ENTRY SDINTF(DS)
!
!        FUNCTION USED TO COMPUTE THE INTEGRAL IN THE ENTROPY
!        COMPUTATION BASED ON DENSITY DS
!
      T2=1.0
      IF (IHE.EQ.1) T2=TS
      B1=DS*DS
      B2= EXP(-CP20*B1/T2)
      A15=0.0
      IF (IHE.EQ.1) A15=A12*B2*(B1/T2+1./CP20)/2.+A13*B2*(B1*B1/(2.*T2) &
     & +B1/CP20+A10/CP20)
      SDINTF= (((A9*DS+A8)*DS + A4)*DS+A3)*DS+B2*A7*(A5+A6*(B1+A10)+A66*&
     & (B1*B1+2.*B1/CP20+2./(CP20*CP20)))+A15
      HDINTF=SDINTF
      RETURN
      END
!$IBFTC STEMPS
      SUBROUTINE TEMPPS (KU,P,S,T,D,DL,DV,KR )
!     ---------------------------VERSION 2/1/72-------------------------
      COMMON/CON123/ DCONV(5),TCONV(5),PCONV(5)
      COMMON/CONV4/ SCONV(5)
      COMMON/BEND31/PS
      COMMON/BEND55/SS
      COMMON/CHECKS/ DCH1,DCH2,PCH1,PCH2,PCH3,TCH1,TCH2,TCH3,DST,TST,   &
     & HSCH1,HSCH2
      COMMON /IERROR/IROUT
      REAL tpsf
      EXTERNAL TPSF
      IROUT=7
      PS=PCHECK(KU,KR,P)
      SS = S/SCONV(KU)
   40 IF (PS-PCH2) 140,140,130
  130 TS1=TCH1
      TS2=TCH3
      GO TO 110
  140 CALL TEMP(1,PS,ZE,TS,1)
      IF (KR-1) 50,70,50
   50 CALL ENT ( 1,1,TS,PS,ZE ,ZE ,SSL,SSV )
      IF ( SS-SSL )  90,70,60
   60 IF ( SS-SSV) 70,70,100
!
!        REGION 1
!
   70 KR=1
   80 CALL DENS(1,TS,ZE,ZE,DSL,DSV,1)
      DL=DSL*DCONV(KU)
      DV=DSV*DCONV(KU)
      GO TO 120
!
!        REGION 2
!
   90 KR=2
      TS1=TCH1
      PS=PS*1.00011
      TS2= TS-1.E-5
      GO TO 110
!
!        REGION 3
!
  100 KR=3
      TS1= TS+1.E-5
      PS= PS*.99988
      TS2=TCH3
!
!        REGIONS 2 AND 3
!
  110 CALL ROOTX ( TS1,TS2,SS,TPSF,TS )
      CALL DENS(1,TS,PS,DS,ZE,ZE,KR)
      D=DS*DCONV(KU)
!
!        VERIFY REGION
!
      IF (PS-PCH2) 120,120,150
  150 IF (TS-TCH2) 160,160,170
  160 KR=2
      GO TO 120
  170 KR=3
  120 T=TS*TCONV(KU)
      RETURN
      END
!$IBFTC TSHF1
      FUNCTION TSHF(TS)
!     ---------------------------VERSION 2/1/72-------------------------
      COMMON/BEND31/ PS
      COMMON/BEND33/ HS
      COMMON/BEND55/ SS
      KR=0
      CALL DENS(1,TS,PS,DS,ZE,ZE,KR)
      CALL ENTH(1,KR,TS,PS,DS,HSC,ZE,ZE)
      TSHF=HSC
      RETURN
      ENTRY TPSF(TS)
      KR = 0
      CALL DENS(1,TS,PS,DS,ZE,ZE,KR)
      CALL ENT(1,KR,TS,PS,DS,SSC,ZE,ZE)
      TSHF = SSC
      RETURN
      END
!$IBFTC HSLV1
      SUBROUTINE HSLV(PS)
!
!        COMPUTE SATURATED LIQUID AND VAPOR ENTHALPIES GIVEN
!        PRESSURE PS AND TEMPERATURE
!
!     ---------------------------VERSION 2/1/72-------------------------
      COMMON/BEND17/TS
      COMMON/BEND21/ SSL,SSV
      COMMON/BEND23/ DSL,DSV
      COMMON/BEND29/ HSL,HSV
      CALL DENS(1,TS,ZE,ZE,DSL,DSV,1)
      HSV=HSS(PS,DSV)
      HSL=HSV-TS*HSSLVF(PS)
      RETURN
      ENTRY SSLV(PS)
!
!     ------------------VERSION  2/23/71--------------------------------
!        COMPUTE SATURATED LIQUID AND VAPOR ENTROPIES GIVEN
!        PRESSURE PS AND TEMPERATURE
!
      CALL DENS(1,TS,ZE,ZE,DSL,DSV,1)
      SSV=SSS(PS,DSV)
      SSL=SSV-HSSLVF(PS)
      RETURN
      END
!$IBFTC TEMPP1
      SUBROUTINE TEMPPH(KU,P,H,T,D,DL,DV,KR)
!     ---------------------------VERSION 2/1/72-------------------------
      COMMON/CON123/DCONV(5),TCONV(5),PCONV(5)
      COMMON/CONV6/HCONV(5)
      COMMON/BEND31/PS
      COMMON/BEND33/HS
      COMMON/CHECKS/DCH1,DCH2,PCH1,PCH2,PCH3,TCH1,TCH2,TCH3,DST,TST,    &
     & HSCH1,HSCH2
      COMMON /IERROR/IROUT
      REAL tshf
      EXTERNAL TSHF
      IROUT=6
      PS=PCHECK(KU,KR,P)
      HS=H/HCONV(KU)
      IF (HS-HSCH1) 20,10,10
   10 IF (HS-HSCH2) 40,40,20
   20 WRITE(6,21) HS
   21 FORMAT(E15.6,' INPUT H IS OUT OF RANGE--CALC.CONT.' )
   40 IF (PS-PCH2) 140,140,130
  130 TS1=TCH1
      TS2=TCH3
      GO TO 110
  140 CALL TEMP(1,PS,ZE,TS,1)
      IF (KR-1) 50,70,50
   50 CALL ENTH(1,1,TS,PS,ZE,ZE,HSL,HSV)
      IF (HS-HSL) 90,70,60
   60 IF (HS-HSV) 70,70,100
!
!        REGION 1
!
   70 KR=1
   80 CALL DENS(1,TS,ZE,ZE,DSL,DSV,1)
      DL=DSL*DCONV(KU)
      DV=DSV*DCONV(KU)
      GO TO 120
!
!        REGION 2
!
   90 KR=2
      TS1=TCH1
      PS=PS*1.00011
      TS2= TS-1.E-5
      GO TO 110
!
!        REGION 3
!
  100 KR=3
      TS1= TS+1.E-5
      PS= PS*.99988
      TS2=TCH3
!
!        REGIONS 2 AND 3
!
  110 CALL ROOTX(TS1,TS2,HS,TSHF,TS)
      CALL DENS(1,TS,PS,DS,ZE,ZE,KR)
      D=DS*DCONV(KU)
!
!        VERIFY REGION
!
      IF (PS-PCH2) 120,120,150
  150 IF (TS-TCH2) 160,160,170
  160 KR=2
      GO TO 120
  170 KR=3
  120 T=TS*TCONV(KU)
      RETURN
      END
!$IBFTC SPCHTV
      SUBROUTINE SPCHV(KU,KR,KCP,T ,P ,D ,CV)
!
!     ---------------------------VERSION 2/1/72-------------------------
!
!        COMPUTE SPECIFIC HEAT CV GIVEN TEMPERATURE T, PRESSURE P,
!        AN DENSITY D.  UNITS ARE SPECIFIED BY KU.  REGION IS
!        SPECIFIED BY KR.  IF KR IS SPECIFIED AS 1, THE SATURATED
!        LIQUID AND VAPOR SPECIFIC HEATS, CVL AND CVV RESPECTIVELY,
!        ARE COMPUTED.
!
      DIMENSION CVS(5)
      COMMON/CONV4/ SCONV(5)
      COMMON/BEND38/TSH(6),TS
      COMMON/BEND39/DS
      COMMON/CHECKS/DCH1(1),DCH2,PCH1,PCH2,PCH3,TCH1,TCH2,TCH3,DST,TST, &
     & HSCH1,HSCH2
      COMMON/SWITS/ KSWIT,TSWIT,DIFT,RSWIT
      COMMON /IERROR/IROUT
      IROUT=10
      TS=T
      PS=P
      DS=D
      IF (KR-1) 100,10,100
!
!        REGIONS 2 AND 3
!
  100 IF(PS-PCH2) 40,40,150
  150 TSS = TCH2
      GO TO 160
   40 CALL TEMP(1,PS,ZE,TSS,1)
  160 GO TO (10,20,30),KR
!
!        REGION 1
   10 GO TO (12,14),KCP
   12 TSH(1)=TS
      CALL CVPS(1,2,CVS)
      CV=CVS(1)
      RETURN
   14 TSH(1)=TS
!
      CALL CVPS(1,3,CVS)
      CV=CVS(1 )
      RETURN
!
!        REGION 2
!
   20 TSH(1)=TS
      IF (TSS.LE.TS-DIFT) TSH(1)=TS-DIFT
      IF (TSS.LE.TS-2.*DIFT) TSH(1)=TS-2.*DIFT
      GO TO 90
!
!        REGION 3
!
   30 TSH(1)=TS-2.*DIFT
      IF (TSS.GE.TS-2.*DIFT) TSH(1)=TS-DIFT
      IF (TSS.GE.TS-DIFT) TSH(1)=TS
   90 CONTINUE
      CALL CVPS(1,KR,CVS)
      CV=CVS(1)
      RETURN
      END
!$IBFTC CVPSS
      SUBROUTINE CVPS(KVP,KR,CVS)
!
!
!     ---------------------------VERSION 2/1/72-------------------------
!        ROUTINE USED TO COMPUTE SPECIFIC HEAT CVS GIVEN
!        TEMPERATURE, PRESSURE, AND DENSITY.  WHETHER DENSITY OR
!        PRESSURE IS CONSTANT IS SPECIFIED BY KVP.  REGION IS
!        SPECIFIED BY KR.
!
      DIMENSION HS(5),CVS(5)
      COMMON/SWITS/ KSWIT,TSWIT,DIFT,RSWIT
      COMMON/BEND37/PS
      COMMON/BEND38/TSH(6),TS
      COMMON/BEND39/DS
      INTEGER:: jump
      REAL TSdummy(1)
      IF (KVP-1) 60,70,60
   70 jump=1   !!! ASSIGN 80 TO JUMP
      GO TO 90
   60 jump=2   !!! ASSIGN 100 TO JUMP
   90 DO 50 I=1,5
!!!      GO TO JUMP,(80,100)
      IF (jump==1) GO TO 80
      IF (jump==2) GO TO 100
   80 CALL PRESS(1,TSH(I),DS,PSH,KRH)
      CALL ENTH(1,KRH,TSH(I),PSH,DS,HS(I),HSL,HSV)
      IF (KRH.EQ.1) GO TO 95
      HS(I)=HS(I)-PSH/DS
      UHS = HS(I)
      GO TO 50
   95 IF (KR.EQ.2) HS(I)=HSL-PSH/DS
      IF (KR.EQ.3) HS(I)=HSV-PSH/DS
      KRH=0
      GO TO 50
  100 CALL DENS(1,TSH(I),PS,DSH,DLS,DVS,KRH)
      CALL ENTH(1,KRH,TSH(I),PS,DSH,HS(I),HSL,HSV)
  110 IF (KRH-1) 50,10,50
   10 IF (KR-2) 30,20,30
   20 HS(I)=HSL
      GO TO 40
   30 HS(I)=HSV
   40 KRH=0
   50 TSH(I+1)=TSH(I)+DIFT
                      ! it must be an array to call SPLINA
      TSdummy(1)=ts
      CALL SPLINA(TSH,HS,5,TSdummy,1,CVS,1,KERROR)
      RETURN
      END
!$IBFTC PTRHO1
      SUBROUTINE PTRHO(DS,TS)
!
!     ---------------------------VERSION 2/1/72-------------------------
      COMMON/BEND/ R,CP1,CP2,CP3,CP4,CP5,CP6,CP7,CP8,CP9,CP10,          &
     & CP11,CP12,CP13,CP14,CP15,CP16,CP17,CP18,CP19,CP20,               &
     & CP21,CP22,CP23,CP24,CP25,CP26,CP27,CP28
      COMMON/BEND35/A2,A3,A4,A5,A6A,A6B,A7A,A7B,A8A,A8B,A11,A12,A13,A14,&
     &A15,A16,A16A
      COMMON/HELFLU/ IHE, IFL, IHY
      COMMON/PARTLS/PTV,PDT
      T2=1.
      IF (IHE.EQ.1) T2=TS
      ADDFL=A16A*A2*(7.-2.*CP20*DS*DS)
      PDT= (((((ADDFL   +6.*CP13-2.*CP20*A16*A2/T2)*DS+5.*A14)*DS+(5.*A1&
     &6-2.*CP20/T2*A15)*A2+4.*A13)*DS+3.*A12)*DS+(2.*A11+3.*A15*A2))*DS+&
     &R*TS
      A9=1./TS
      A10=A9*A9
      A6=A6A+A6B
      A7=A7A+A7B
      A8=A8A+A8B
      A17=0.0
      IF (IHE.EQ.1) A17=(A6A+A7A*A9+A8A*A10)*CP20*DS*DS*(A10*A10)
      PTV=A4-A10*((((6.*CP25*A10*DS*DS+4.*A8)*A9+3.*A7)*A9+2.*A6)*A9+A5)&
     &+A17
      RETURN
      END
!$IBFTC SONIC1
      SUBROUTINE SONIC(KU,KR,T,D,GAMMA,C)
!
!     ---------------------------VERSION 2/1/72-------------------------
!        COMPUTE SONIC VELOCITY C GIVEN TEMPERATURE T, DENSITY D,
!        AND SPECIFIC HEAT RATIO GAMMA.  UNITS ARE SPECIFIED BY KU.
!        REGION IS SPECIFIED BY KR.
!
      COMMON/BEND/ R,CP1,CP2,CP3,CP4,CP5,CP6,CP7,CP8,CP9,CP10,          &
     & CP11,CP12,CP13,CP14,CP15,CP16,CP17,CP18,CP19,CP20,               &
     & CP21,CP22,CP23,CP24,CP25,CP26,CP27,CP28
      COMMON/CONV5/CCONV(5)
      COMMON/HELFLU/ IHE, IFL, IHY
      COMMON/PARTLS/PTV,PDT
      COMMON /IERROR/IROUT
      IROUT=8
      TS=TCHECK(KU,KR,T)
      T2=1.0
      IF (IHE.EQ.1) T2=TS
      DS=DCHECK(KU,D)
!
!        REGIONS 1, 2, AND 3
!
      GAMPR=GAMMA*10.*PDT
      CS=0.
      IF (GAMPR.GT.0.)CS=1000.*SQRT(GAMPR)
      C=CS*CCONV(KU)
      RETURN
      END
!$IBFTC CPPRLF
      FUNCTION CPPRLF(DS)
!
!     ---------------------------VERSION 2/1/72-------------------------
      COMMON/BEND/ R,CP1,CP2,CP3,CP4,CP5,CP6,CP7,CP8,CP9,CP10,          &
     & CP11,CP12,CP13,CP14,CP15,CP16,CP17,CP18,CP19,CP20,               &
     & CP21,CP22,CP23,CP24,CP25,CP26,CP27,CP28
      COMMON/BEND17/TS
      COMMON/HELFLU/ IHE, IFL, IHY
!
!          INTEGRAL OF ((T /RHO**2)(D2P/DT2) AT CONSTANT RHO)
!
      COMMON/BEND35/A2,A3,A4,A5,A6A,A6B,A7A,A7B,A8A,A8B,A11,A12,A13,A14,&
     &A15,A16,A16A
      T2=1.0
      IF (IHE.EQ.1) T2=TS
      B2 = DS*DS
      B1=-EXP(-CP20*B2/T2)/(2.*CP20*TS/T2)
      A1 = 1./TS
      D2=B2+ T2/CP20
      D3=B2*B2+2.*B2/CP20+2./(CP20*CP20)
      D4 = A1*A1
      A17=(10.*CP16*A1+6.*CP15)*A1+3.*CP14
      A18=(10.*CP19*A1+6.*CP18)*A1+3.*CP17
      A18A=(10.*CP28*A1+6.*CP27)*A1+3.*CP26
      A19=T2/CP20
      A20=A19*A19
      A21=(10.*CP16*A1+8.*CP15)*A1+6.*CP14
      A22=(10.*CP19*A1+8.*CP18)*A1+6.*CP17
      D5=1./CP20+B2/T2+B2*B2/(2.*T2*A19)
      D52=D5/CP20
      D7=((B2/A19+3.)*B2+6.*A19)*B2+6.*A20
      B3= EXP(-CP20*B2/T2)
      TERM=0
      IF (IHE.NE.1) GO TO 10
      TERM=B3*(D2*A21/(2.*T2*T2*T2)-A15*D5+A22*D52/T2-A16*D7/(2.*T2))
   10 CONTINUE
      CPPRLF=2.*D4*(((((21.*CP25*D4+10.*CP21)*A1+6.*CP5)*A1+3.*CP4)*A1+C&
     &P3)*DS+(2.*CP22*TS+CP8+3.*CP23*A1+CP24*B2/2.)*B2/2.+(A17+D2*A18+D3&
     &*A18A)*B1)+TERM
      RETURN
      END
!$IBFTC HEK
!     ---------------------------VERSION 2/1/72-------------------------
      FUNCTION CONZ(TEMP)
!     KZERO TABLE IN MW/CM-K, T IN KELVIN, DECK OF 18 SEPT 70
!           COURTESY OF HANS RODER     NBS LAB  BOULDER COLORADO
      DIMENSION A(55),T(55)
      DATA A/0.0376,0.0496,0.0608,0.0720,0.0829,0.0932,0.1028,0.1200,   &
     & 0.1350,0.1486,0.1610,0.1725,0.1933,0.2121,0.2292,0.2448,0.2594,  &
     & 0.2983,0.3357,0.3709,0.4042,0.4360,0.4666,0.5259,0.5821,0.6359,  &
     & 0.6877,0.7381,0.7862,0.8333,0.8791,0.9241,0.9680,1.0113,1.054,   &
     & 1.096,1.136,1.176,1.216,1.255,1.293,1.331,1.369,1.406,1.443,     &
     & 1.479,1.516,1.551,                 1.912,2.208,2.4836,2.743,2.99,&
     &3.226,3.453/
      DATA T/1.9,2.5,3.0,3.5,4.0,4.5,5.0,6.,7.,8.,9.,10.,12.,14.,16.,   &
     & 18.,20.,25.,30.,35.,40.,45.,50.,60.,70.,80.,90.,100.,110.,120.,  &
     & 130.,140.,150.,160.,170.,180.,190.,200.,210.,220.,230.,240.,250.,&
     & 260.,270.,280.,290.,300.,400.,500.,600.,700.,800.,900.,1000./
      TT= TEMP
    9 DO 11  J=1,58
   10 IF (tt <= t(j)) GO TO 12   !!! IF(TT-T(J))12,12,11
   11 END DO
   12 CONTINUE
   13 CONZ=A(J-1)+(TT-T(J-1))*(A(J)-A(J-1))/(T(J)-T(J-1))
      RETURN
      END
!$IBFTC CPPRL1
      SUBROUTINE CPPRL(P,D,T,CPPART,CVPART,KU,KR,KCP)
!
!     ---------------------------VERSION 2/1/72-------------------------
      COMMON/BEND/ R,CP1,CP2,CP3,CP4,CP5,CP6,CP7,CP8,CP9,CP10,          &
     & CP11,CP12,CP13,CP14,CP15,CP16,CP17,CP18,CP19,CP20,               &
     & CP21,CP22,CP23,CP24,CP25,CP26,CP27,CP28
      COMMON/COCPO/TO,CCOP1,CCOP2,CCOP3,CCOP4,CCOP5
      COMMON/BEND17/TT
      COMMON/BEND35/A2,A3,A4,A5,A6A,A6B,A7A,A7B,A8A,A8B,A11,A12,A13,A14,&
     &A15,A16,A16A
      COMMON/HELFLU/ IHE, IFL, IHY
      COMMON/BEND36/CPO
      COMMON/PARTLS/PTV,PDT
      COMMON/REFNO/SOTO,HOTO,CPOCOR,HTERM,STERM
      COMMON/SWITS/ KSWIT,TSWIT,DIFT,RSWIT
      COMMON/CONV4/SCONV(5)
      INTEGER:: jump

      TS=TCHECK(KU,KR,T)
      PS=PCHECK(KU,KR,P)
      DS=DCHECK(KU,D)
      IF (IHY.EQ.1) CALL SETCPO(TS)
        TT=TS
      IF (KR.EQ.1) KCP=KCP+1
      jump=1   !!! ASSIGN 500 TO JUMP
      IF (DS.GT.RSWIT) jump=2   !!! ASSIGN 600 TO JUMP
      IF (KSWIT.EQ.2.AND.(TS.LE.TSWIT.OR. DS.GT.RSWIT)) jump=2   !!! ASSIGN 600 TO JUMP
!!!      GO TO JUMP,(500,600)
      IF (jump==1) GO TO 500
      IF (jump==2) GO TO 600
  500 A9 = 1./TS
      A10 = A9*A9
      A15 = ((CP16*A9 + CP15) * A9 + CP14) * A10
      A16 = ((CP19*A9 + CP18)*A9 +CP17)*A10
!
!           TEST CP COMPUTED BY THE INTEGRAL OF THE SECOND
!           PARTIAL OF P BY T
!
      CPO = (((CCOP5*TS+CCOP4) *TS+CCOP3)*TS+CCOP2)*TS+CCOP1
      CVO= CPO*CPOCOR-R
      CVPART=    -(CPPRLF(DS) - CPPRLF(0.)) + CVO
      GO TO 1000
  600 CALL SPCHV(KU,KR,KCP,TS,PS,DS,CVPART)
 1000 A1= DS*DS
      T2=1.0
      IF (IHE.EQ.1) T2=TS
      A2= DS*EXP(-CP20*DS*DS/T2)
      A3=((((CP13*DS+CP12)*DS+CP10)*DS+CP7)*DS+CP2)*A1
      A4=((((CP11*DS+CP9)*DS+CP6)*DS+CP1)*DS+R) *DS
      A4=A4+2.*CP22*A1*DS*TS
      A5=((CP24*A1+CP8)*DS+CP3)*A1
      A6A=((CP26*A1+CP17)*A1+CP14)*A1*A2
      A6B= (CP23*DS+CP4)*A1
      A7A=((CP27*A1+CP18)*A1+CP15)*A1*A2
      A7B=CP5*A1
      A8A=((CP28*A1+CP19)*A1+CP16)*A1*A2
      A8B=CP21*A1
      A9 = 1./TS
      A10 = A9*A9
      A11=((((CP25*A10+CP21)*A9+CP5)*A9+CP4)*A9+CP3)*A9+CP2+CP1*TS
      A12=(CP23*A9+CP8)*A9+CP7+(CP22*TS+CP6)*TS
      A13 = CP10 + CP9 * TS
      A14= CP11*TS+CP12 +CP24*A9
      A15 = ((CP16*A9 + CP15) * A9 + CP14) * A10
      A16 = ((CP19*A9 + CP18)*A9 +CP17)*A10
      A16A=((CP28*A9+CP27)*A9+CP26)*A10
      CALL PTRHO(DS,TS)
      CPPART = CVPART+ TS*PTV*PTV/(PDT*DS*DS)
      IF (KCP.EQ.2) KCP=0
      CPPART=CPPART*SCONV(KU)
      CVPART=CVPART*SCONV(KU)
      RETURN
      END
!C$IBFTC SCPO
!     ---------------------------VERSION 2/1/72-------------------------
      SUBROUTINE  SETCPO(TIN)
!         NEEDED FOR STEPWISE INTEGRATION
      DIMENSION ADDH(4),ADDS(4),HCPO(5,3)
      DATA ADDH/0.,-62.9775,-147.6,61.7/
      DATA ADDS/0.,-9.9313,5.361,-16.707/
      COMMON/COCPO/TO(1),CCOP1,CCOP2,CCOP3,CCOP4,CCOP5
      COMMON/HYDRO/DELH,DELS
      DATA HCPO/                                                        &
     &14.7599360,-.21977388,.32100769E-2,-.12061502E-4,.57121808E-8,    &
     &6.6557899 ,.15621077,-.86913643E-03,.18972274E-05,-.14418461E-08, &
     &14.4114861,-.71767870E-03,.18638538E-05,-.53065470E-09,.46649305E-&
     &13/
      TEST=TIN
!
!        MULTIPLE FIT REQUIRES STEPWISE INTEGRATION FOR H2
!
      IF (TEST.LT.40.) GO TO 40
   21 FORMAT(' ',2E16.8)
      K=1
      IF (TEST.LE.150.) GO TO 5
      K=2
      IF (TEST.LE.500.) GO TO 5
      K=3
    5 DO 10 I=2,6
   10 TO(I) = HCPO(I-1,K)
      GO TO 50
   40 DO 41 I=3,6
   41 TO(I)=0.0
      TO(2)=4.968*4.184/2.01572
      RETURN
   50 DELH=ADDH(K+1)
      DELS=ADDS(K+1)
      RETURN
      END
!C$IBFTC HSS1
      FUNCTION HSS(PS,DS)
!
!     ---------------------------VERSION 2/1/72-------------------------
!        COMPUTE ENTHALPY IN REGION 3 OR SATURATED VAPOR ENTHALPY
!        GIVEN PRESSURE PS, DENSITY DS, AND TEMPERATURE
!
      COMMON/BEND/ R,CP1,CP2,CP3,CP4,CP5,CP6,CP7,CP8,CP9,CP10,          &
     & CP11,CP12,CP13,CP14,CP15,CP16,CP17,CP18,CP19,CP20,               &
     & CP21,CP22,CP23,CP24,CP25,CP26,CP27,CP28
      COMMON/HELFLU/ IHE, IFL, IHY
      COMMON/COCPO/TO   ,CCOP1,CCOP2,CCOP3,CCOP4,CCOP5
      COMMON/BEND17/TS
      COMMON/BEND25/ A1,A2,A3,A4
      COMMON/HYDRO/DELH,DELS
      COMMON/REFNO/SOTO,HOTO,CPOCOR,HTERM,STERM
      HTINTF(TX)=((((A1*TX+A2)*TX+A3)*TX+A4)*TX+CCOP1)*TX
      STINTF(TX)=(((A1*TX+A2)*TX+A3)*TX+CCOP2)*TX+CCOP1*ALOG(TX)
      DELH=0.0
      IF (IHY.EQ.1) CALL SETCPO(TS)
      A1=CCOP5/5.
      A3=CCOP3/3.
      A4=CCOP2/2.
      A2=CCOP4/4.
      HTERM1=HTINTF(TS)
      HSS = HOTO + (HTERM1    +DELH-HTERM )*CPOCOR + PS/DS              &
     & -R*TS+HDINT(DS, 0.0)
      RETURN
      ENTRY  SSS(PS,DS)
!
!        COMPUTE ENTROPY IN REGION 3 OR SATURATED VAPOR ENTROPY
!        GIVEN PRESSURE PS, DENSITY DS, AND TEMPERATURE
!
!  THE TERM R*LN(RTD) IS NOT DIMENSIONLESS   CONSEQUENTLY  THE TERM
!  SHOULD BE       R*LN(RTD*CONV) = R*LN(CONV) + R*LN(RTD)
!     THEREFORE ADD THE TERM  R*LN(CONV)
!          FOR   R  ATM-LITER/(GM-MOLE K)  DINNS  VALUE IMPLIES CONV=1
!          FOR   R  PSIA-CC/(GM  K)                   CONV= 14.696
!          FOR   R  J/(GM K)                          CONV= .101325
!    FOR  NITROGEN   THE CORRECTION BECOMES
!   BENDER   SCORR=R*ALOG(.101325)  =-.6794936
      DELS=0.0
      IF (IHY.EQ.1) CALL SETCPO(TS)
      A1=CCOP5/4.
      A2=CCOP4/3.
      A3=CCOP3/2.
      STERM1=STINTF(TS)
      SSS= SOTO + (STERM1 +   DELS-STERM   )*CPOCOR-R*ALOG(R*TS*DS)     &
     & + SDINT(DS,0.)
      HSS=SSS
      RETURN
      END
!C$IBFTC SURFAC
      SUBROUTINE SURF(KU,KR,T,SIGMA)
      DIMENSION STCONV(5)
      DATA STCONV/2*1.0, 6.8521766E-5, 2*1.0/
      COMMON/SURCON/PCTC,TCRIT,FIXIT,ZET
      COMMON/HELFLU/IHE,IFL,IHY
!
!     ---------------------------VERSION 2/1/72-------------------------
!      SURFACE TENSION  DYNE/CM
!
      COMMON /IERROR/IROUT
      IROUT=13
      SIGMA=0.
      TS=TCHECK(KU,1,T)
      IF ( TS .GT. TCRIT) RETURN
      IF (IHY.EQ.1) GO TO 10
      SIG=PCTC*(1.-TS/TCRIT)**(11./9.)*(.432/ZET-.951)*FIXIT
      SIGMA=SIG  *STCONV(KU)
      RETURN
   10 TR=TS/TCRIT
      SIGMA=5.369*(1.0-TR)**1.065*STCONV(KU)
      RETURN
      END
!C$IBFTC VISC
      SUBROUTINE VISC(KU,KR,T,D,MU)
!
!     ---------------------------VERSION 2/1/72-------------------------
!        COMPUTE VISCOSITY MU GIVEN TEMPERATURE T AND DENSITY D.
!        UNITS ARE SPECIFIED BY KU.  REGION IS SPECIFIED BY KR.
!
!      (U-U*)Z PARAMETER IS CALCULATED FROM JOSSI,STIEL,AND THODOS
!          VOL.8 NO.1   A.I.CH.E.JOURNAL  PAGE 60
!
      DIMENSION COMEGA(15),MCONV(5)
      REAL MUS,MCONV,MU,MUMZ
      DATA (COMEGA(I),I=1,15)/2.,4.,-.52288,-1.33115022,.56160244,      &
     &1.49894829,-1.84681526,1.58946684,.60206,.32481334E-1,            &
     &-.23524079,.65223929,-1.03095995,1.40383078,2.60206/
      COMMON/COFMU/EPSK,WM,DIS,RHOCRT,ZETA
      DATA MCONV/2*1.,.67196899E-1,2*1./
      COMMON/HELFLU/IHE,IFL,IHY
      COMMON /IERROR/IROUT
      IROUT=12
      TS=TCHECK(KU,KR,T)
      DS=DCHECK(KU,D)
      IF (IHY.EQ.1) GO TO 20
!
!        REGIONS 1, 2, AND 3
!
      XMS=.26693E-4*SQRT(WM*TS)/(DIS*DIS*POLY(ALOG10(TS*EPSK),COMEGA))
      RHOR=DS/RHOCRT
      MUMZ=((RHOR*(RHOR*(RHOR*(RHOR*.0093324-.040758)+.058533)+.023364) &
     &+.10230)**4-1.E-4)/ZETA/10.**2
!!!      DEBUG XMS,MUMS,DIS,RHOR,EPSK,DIST,ZETA
!
!             CORRECTION APPLIED TO FLUORINE DERIVED BY RCH FROM HANLEYS
!
      IF (IFL.EQ.1) MUMZ=MUMZ*(1.-.25*SIN(2.3808*ALOG(RHOR)))
      MUS=MUMZ+XMS
      GO TO 30
   20 XMS=8.5558*(TS**1.5/(TS+19.55)*(TS+650.39)/(TS+1175.9))*1.E-6
      AFUNC=EXP(5.7694+ALOG(DS)+65.*DS**1.5-6.E-6*EXP(127.2*DS))
      DSTERM=DS/.07
      RANGE=-58.75*(DS/.07)**3
      IF (RANGE.LT.-80.) RANGE=-80.
      BFUNC=10.0+7.2*(DSTERM**6-DSTERM**1.5)-17.63*EXP(RANGE)
      MUS=XMS+AFUNC*EXP(BFUNC/TS)*1.E-6
   30 MU=MUS*MCONV(KU)
      RETURN
      END
!C$IBFTC THER1
      SUBROUTINE THERM (KU,KR,P,T,D,EXCESK, K)
!
!     ---------------------------VERSION 2/1/72-------------------------
!        COMPUTE THERMAL CONDUCTIVITY K GIVEN TEMPERATURE T AND
!        DENSITY D.  UNITS ARE SPECIFIED BY KU.  REGION IS SPECIFIED
!        BY KR.
!
      COMMON/HELFLU/ IHE, IFL, IHY
      COMMON/TCOND/CKMKST(9), XLAMB,ZC5, RHOCR, TCR, TCSTAR, CKSTAR(18)
      COMMON/BEND/R,DUMY(28)
      COMMON/COFMU/EPSK,WM,DIS,RHOCRT,ZETA
      COMMON/REFNO/SOTO,HOTO,CPOCOR,HTERM,STERM
      COMMON/COCPO/TO,CCOP1,CCOP2,CCOP3,CCOP4,CCOP5
      DIMENSION OMEGA1(9),OMEGA2(15)
      REAL KMKS, KSRAT, K, KCONV ,KZERO
      DIMENSION KCONV(5)
      DATA (OMEGA1(I),I=1,9)/1.,4.,-.5228,.82359030E-1,-.53982316,      &
     &1.28552593,-1.53509952,1.44223760,2.602 /
      DATA (OMEGA2(I),I=1,15)/2.,4.,-.52288,-1.33115022,.56160244,      &
     &1.49894829,-1.84681526,1.58946684,.60206,.32481334E-1,            &
     &-.23524079,.65223929,-1.03095995,1.40383078,2.60206/
      DATA KCONV/2*1.,.01606044,2*1./
      COMMON /IERROR/IROUT
      IROUT=11
      TS=TCHECK(KU,KR,T)
      DS=DCHECK(KU,D)
!
!        REGIONS 1, 2, AND 3
!
      IF (IHE.EQ.1) GO TO 30
      IF (IHY.EQ.1) GO TO 38
      RHO = DS/RHOCR
      ARHO=ALOG10(RHO)
      IF (ARHO.GT.CKMKST(3)) GO TO 10
      KMKS=10.**(ARHO-7.0)
      GO TO 20
   10 KMKS = 10.**(POLY(ALOG10(RHO),CKMKST)) / (ZC5*XLAMB)
   20 IF (IFL.EQ.1) GO TO 35
      TR = TS/TCR
      KSRAT = 10.** (POLY (ALOG10(TR),CKSTAR))*TCSTAR
      K= (KMKS + KSRAT) * KCONV(KU) * 4.184
      GO TO 40
   30 SLOPE=10.**(DS*(DS*(DS*(-621.369011)+224.2564)-29.48514)+2.0941961&
     &)
      KMKS= ALOG10(TS)/SLOPE*1.E-3
!       USE HANS RODERS  KO
      KZERO=CONZ(TS) /1000.
      K=(KMKS+KZERO)*KCONV(KU)
      GO TO 40
!
!    FLUORINE KZERO CURVE BASED ON MASON-MONCHICK ANALYSIS
!          FROM REID AND SHERWOOD PAGE 461 EQ 10-12
!
   35 OM1=POLY(ALOG10(TS*EPSK),OMEGA1)
      OM2=POLY(ALOG10(TS*EPSK),OMEGA2)
      OM21=6.*OM2/OM1/5.
      XMO=.26693E-4*SQRT(WM*TS)/(DIS*DIS*OM2)
      CVO=CPOCOR*((((CCOP5*TS+CCOP4)*TS+CCOP3)*TS+CCOP2)*TS+CCOP1)
      KZERO=(1.77*R+OM21*CVO-2./3.14159*(2.5-OM21)**2*(CVO-1.5*R)/      &
     &(2.+.02*TS**.83333))*XMO
      K=(KMKS+KZERO)*KCONV(KU)
      GO TO 40
!    USES RODERS HYDROGEN CALCULATION BUT NOT HIS CRITICAL SCALING FOR
!       REACTING CONDUCTIVITY NEAR THE CRITICAL REGION.
!
   38 K=CONC(TS,DS)/1.E3*KCONV(KU)
!
!          REACTING CONDUCTIVITY CALCULATED BY THE SENGERS-KEYES METHOD
!  GENERAL FLUID
!
   40 DRHOC=ABS(DS-RHOCR)/RHOCR
      DELAMB=0.0
      IF (DRHOC.GT..6) GO TO 50
      DELTC=ABS((TS-TCR)/TCR)
      RAT=DS/RHOCR
      IF (DRHOC.LT..00001) GO TO 102
      IF (DELTC.LT.1.E-7) GO TO 101
      XBETA= DELTC**.35/DRHOC
      IF(XBETA.GT..4)  GO TO 104
  101 DELAMB= 3.05E-5/(SQRT(RAT)*DRHOC**1.71)
      GO TO 50
  104 IF(XBETA.GT.3.) GO TO 102
       DELAMB=3.05E-5/(SQRT(RAT)*DELTC**.6)/(1.+.90/XBETA**(1./.35))    &
     &**(3./5.)
      GO TO 50
  102 IF (DELTC.LT.1.E-7) GO TO 105
      DELAMB=3.05E-5/(SQRT(RAT)*DELTC**.6)
   50 EXCESK=DELAMB*KCONV(KU)
      RETURN
  105 EXCESK=1.E30
      RETURN
      END
!C$IBFTC HYDROK
      FUNCTION CONC(TMEAN,DEN)
!
!    THIS ROUTINE BY HANS RODER HAS BEEN CHANGED TO COMMENT OUT THE
!     REACTING CONDUCTIVITY CALCULATIONS WHICH ARE NOT USED AT PRESENT
!
!     T KELVIN, DEN MUST BE IN GR/CC, CP IN CAL/MOLE
!     THIS VERSION CHANGED TO INTERPOLATE LINEARLY IN K-0
      DIMENSION AAA(12),TTT(12),D(3)
      DATA TTT /13.000,17.020,19.587,25.100,30.010,33.063,              &
     & 40.145,59.187,79.845,99.852,122.909,153.0/
      DATA AAA /0.1765,0.3203,0.4121,0.6679,0.8491,0.9363,              &
     & 1.1153,1.4781,1.7895,2.0757,2.3357,2.5450/
      DATA D /9.88531118E-01, 3.20886940E+01, -9.10140989E+02/
      CONC=0.
!       THIS CALCULATION BLOWS UP IF EXTRAPOLATED TOO FAR
      IF (TMEAN .GT. 248.) RETURN
      DO 9 J=1,12
      IF (TMEAN .LE. TTT(J)) GO TO 10
    9 END DO
                    ! added by RLC 16Sept1997
   10 j=MIN(j,12)
      CONZ1=0.1*EXP(AAA(J-1))
      CONZJ=0.1*EXP(AAA(J))
      AAVE=CONZ1+(TMEAN-TTT(J-1))*(CONZJ-CONZ1)/(TTT(J)-TTT(J-1))
      AAVE=LOG(10.0*AAVE)
! 203 CADJ=TCADJ(TMEAN,DEN,CP)
      BB= 39.6-2.0*SQRT(248.0-(TMEAN/10.0-17.)**2)
      FUNC=AAVE+D(1)*BB*DEN+D(2)*DEN*DEN+D(3)/(TMEAN-9.0)*DEN*DEN
      CONC=0.1*EXP(FUNC)
      RETURN
      END
