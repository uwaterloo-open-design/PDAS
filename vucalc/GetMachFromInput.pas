UNIT GetMachFromInput;
(* PURPOSE - For the isentropic flow, normal shock, Rayleigh and Fanno     *)
(*   pages, compute theupstream Mach number from whatever input value has  *) 
(*   been specified. For the oblique shock page, compute the shock angle.  *)

(* AUTHORS - Tom Benson, NASA Glenn Research Center (formerly Lewis)       *)
(*           Ralph L. Carmichael, Public Domain Aeronautical Software      *)

(* REVISION HISTORY                                                        *)
(*   DATE  VERS PERSON  STATEMENT OF CHANGES                               *)
(*    ?     1.0    TB   Original coding (in C for Unix, X-windows, Forms)  *)
(*   1997   2.0   RLC   Recoded in Delphi for Microsoft Windows            *)
(* 21Nov00  5.5   RLC   Added the Fanno and Rayleigh calculations          *)


INTERFACE

FUNCTION GetMachFromIsenInput(const inputCase : INTEGER;
                              const inputValue : EXTENDED;
                              const gamma : EXTENDED) : EXTENDED;

FUNCTION GetMachFromNormalInput(const inputCase : INTEGER;
                                const inputValue : EXTENDED;
                                const gamma : EXTENDED) : EXTENDED;

FUNCTION GetThetaFromObliqueInput(const inputCase : INTEGER;
                                  const inputValue : EXTENDED;
                                  const mach1 : EXTENDED;
                                  const gamma : EXTENDED) : EXTENDED;
FUNCTION GetMachFromRayleighInput(const gamma : EXTENDED;
                                  const input : EXTENDED;
                                  const itemIndex : INTEGER) : EXTENDED;

FUNCTION GetMachFromFannoInput(const gamma : EXTENDED;
                               const input : EXTENDED;
                               const itemIndex : INTEGER) : EXTENDED;

IMPLEMENTATION

USES CalcSubs, Dialogs, Math, NACA1135, RayFanno;

PROCEDURE ErrorMsg(const s : STRING);
BEGIN
  MessageDlg(s, mtError, [mbOk], 0)
END;   (* -------------------------------------- End of Procedure ErrorMsg *)


FUNCTION GetMachFromIsenInput(const inputCase : INTEGER;
                              const inputValue : EXTENDED;
                              const gamma : EXTENDED) : EXTENDED;
  VAR
    gm1,gp1 : EXTENDED;
    ma : EXTENDED;
    super : BOOLEAN;
BEGIN
  gm1:=gamma-1.0;
  gp1:=gamma+1.0;
  ma:=-1.0;
  CASE inputCase OF
    0: IF inputValue < 0.0 THEN
         ErrorMsg('Mach number cannot be negative')
       ELSE
         ma:=inputValue;
    1: BEGIN                                                   {q/p}
         IF inputValue < 0.0 THEN
           ErrorMsg('q/p cannot be negative')
         ELSE
           ma:=Sqrt(2.0*inputValue/gamma); { inverse of Eq47 }
       END;
    2: BEGIN                                                        {p/pt}
         IF (inputValue > 1.0) OR (inputValue <=0.0) THEN
           ErrorMsg('p/pt cannot be greater than 1 or <= 0')
         ELSE
           ma:=Sqrt(2.0*(Power(inputValue, -gm1/gamma)-1.0)/gm1)
       END;
    3: BEGIN                                                 { A/A*  for M <1}
         IF inputValue < 1.0 THEN
           ErrorMsg('Area ratio cannot be less than 1')
         ELSE
           BEGIN
             super:=FALSE;         {    }
             ma:=GetMachFromAreaRatio(inputValue, gamma, super)
             END
       END;
    4: BEGIN                                                  { A/A* for M>1}
         IF inputValue < 1.0 THEN
           ErrorMsg('Area atio cannot be less than 1')
         ELSE
           BEGIN
             super:=TRUE;
             ma:=GetMachFromAreaRatio(inputValue, gamma, super)
           END
       END;
    5: BEGIN                                                   {rho/rhot}
         IF (inputValue > 1.0) OR (inputValue <= 0.0) THEN
           ErrorMsg('rho/rhot cannot be greater than 1')
         ELSE
           ma:=Sqrt(2.0*(Power(inputValue, -gm1)-1.0)/gm1);
       END;
    6: BEGIN                                               {T/Tt}
         IF (inputValue > 1.0) OR (inputValue <= 0.0) THEN
           ErrorMsg('T/Tt cannot be greater than 1 or less than 0')
         ELSE
           ma:=Sqrt(2.0*(1.0/inputValue-1.0)/gm1)
       END;
    7: BEGIN                                                   {V/a*}
         IF (inputValue < 0.0) OR (Sqr(inputValue) > 6.0)  THEN
           ErrorMsg('V/a* cannot be < 0 or > Sqrt(6)')
         ELSE
           ma:=Sqrt(2.0*Sqr(inputValue)/(gp1*(1.0-gm1*Sqr(inputValue)/gp1)));
       END;
    8: BEGIN                                                {nu}
         IF (inputValue < 0.0) OR (inputValue > 130.45) THEN
           ErrorMsg('Angle cannot be < 0 or > 130.45')
         ELSE
           ma:=GetMachFromPM(inputValue, gamma);
       END;
    9: BEGIN                                         {mu}
         IF (inputValue <=0.0) OR (inputValue >= 90.0) THEN
           ErrorMsg('Mach angle must be between 0 and 90 degrees')
         ELSE
           ma:=1.0/Sin(DegToRad(inputValue))
       END;
  END;  { END of CASE }

  result:=ma  { return -1.0 if failure }
END;   (* ------------------------------- End of Function GetMachFromInput *)


FUNCTION GetMachFromNormalInput(const inputCase : INTEGER;
                                const inputValue : EXTENDED;
                                const gamma : EXTENDED) : EXTENDED;
  VAR
    ma : EXTENDED;
BEGIN
  ma:=-1.0;
  CASE inputCase OF
    0: BEGIN    { upstream Mach }
         IF inputValue < 1.0 THEN
           ErrorMsg('Upstream Mach must not be < 1')
         ELSE
           ma:=inputValue;
       END;
    1: BEGIN                    {downstream Mach }
         IF (inputValue > 1.0) OR (inputValue <= 0.0) THEN
           ErrorMsg('Downstream Mach must not be > 1 or < 0')
         ELSE
           ma:=Eq96inverse(inputValue,gamma)
       END;
    2: BEGIN                         { pt2/pt1 }
         IF (inputValue > 1.0) OR (inputValue <= 0.0) THEN
           ErrorMsg('Ratio must not be > 1 or < 0')
         ELSE
           ma:=GetMachFromPtRatio(inputValue,gamma)
       END;
    3: BEGIN                               {p2/p1}
         IF inputValue < 1.0  THEN
           ErrorMsg('Ratio must not be < 1')
         ELSE
           ma:=Eq93inverse(inputValue,gamma)
       END;
    4: BEGIN                      { T2/T1}
         IF inputValue < 1.0 THEN
           ErrorMsg('T2/T1 must not be < 1')
         ELSE
           ma:=GetMachFromTempRatio(inputValue,gamma)
       END;
    5: BEGIN         { rho2/rho1 }
         IF inputValue < 1.0 THEN
           ErrorMsg('Ratio must not be < 1' )
         ELSE
           ma:=Eq94inverse(inputValue,gamma)
       END;
  END;  { END of CASE }

  result:=ma
END;   (* ------------------------------- End of Function GetMachFromNormalInput *)

FUNCTION GetThetaFromObliqueInput(const inputCase : INTEGER;
                           const inputValue : EXTENDED;
                           const mach1 : EXTENDED;
                           const gamma : EXTENDED) : EXTENDED;
  VAR
    deltaMax : EXTENDED;  { max delta for a given Mach1 }
    theta : EXTENDED;
    thetaMax : EXTENDED;  { the corresponding theta to deltaMax }
    rhoMax : EXTENDED;
    tempMax : EXTENDED;
    ptMin : EXTENDED;
    pmax : EXTENDED;
    m2Min : EXTENDED;
BEGIN
  { compute deltaMax and thetaMax for this upstream Mach number }
  theta:=-1.0;
  MaxWedgeAngleForAttachedShock(mach1,gamma,deltaMax,thetaMax);

  CASE inputCase OF
    0: BEGIN                                 {ramp}
         IF inputValue > RadToDeg(deltaMax) THEN
           ErrorMsg('Ramp angle too large')
         ELSE
           theta:=GetShockAngleFromRampAndMach(DegToRad(inputValue),
                                                         mach1,gamma);
       END;
    1: BEGIN                                            {wave}
         theta:=DegToRad(inputValue);
         IF mach1*Sin(theta) < 1.0 THEN
           ErrorMsg('wave angle too small')  { then what?? }
       END;
    2: BEGIN                                 {pt2/pt1}
         IF (inputValue > 1.0) OR (inputValue <= 0.0) THEN
           ErrorMsg('Total pressure cannot increase or be < 0')
         ELSE
           BEGIN
             ptMin:=Eq99(mach1,gamma);
             IF inputValue < ptMin THEN
               ErrorMsg('Total pressure ratio too small')
             ELSE
               BEGIN
                 theta:=GetMachFromPtRatio(inputValue,gamma);
                 theta:=ArcSin(theta/mach1);
               END
           END
       END;
    3: BEGIN                                        {p2/p1}
         IF inputValue < 1.0 THEN
           ErrorMsg('Static pressure cannot decrease')
         ELSE
           BEGIN
             pmax:=Eq93(mach1,gamma);  { normal shock }
             IF inputValue > pmax THEN
               ErrorMsg('Static pressure ratio too large')
             ELSE
               theta:=ArcSin(Eq93inverse(inputValue,gamma)/mach1);
           END
       END;
    4: BEGIN                                           {T2/T1}
         IF inputValue < 1.0 THEN
           ErrorMsg('Static temperature cannot decrease')
         ELSE
           BEGIN
             tempMax:=Eq95(mach1,gamma);
             IF inputValue > tempMax THEN
               ErrorMsg('Static pressure ratio too large')
             ELSE
               BEGIN
                 theta:=GetMachFromTempRatio(inputValue,gamma);
                 theta:=ArcSin(theta/mach1);
               END
           END
       END;
    5: BEGIN                                           {rho2/rho1}
         IF inputValue < 1.0 THEN
           ErrorMsg('Static density cannot decrease')
         ELSE
           BEGIN
             rhoMax:=Eq94(mach1,gamma);
             IF inputValue > rhoMax THEN
               ErrorMsg('Static density ratio too large')
             ELSE
               BEGIN
                 theta:=Eq94inverse(inputValue,gamma);
                 theta:=ArcSin(theta/mach1);
               END
           END
       END;
    6: BEGIN           { M2 }
         IF inputValue > mach1 THEN
           ErrorMsg('Mach cannot increase through a shock')
         ELSE
           BEGIN
             m2Min:=Sqrt(Eq96(mach1,gamma));  { normal shock }
             IF inputValue < m2Min THEN
               ErrorMsg('Downstream Mach too small')
             ELSE
               theta:=GetThetaFromM1andM2(mach1,inputValue,gamma)
            END
       END;
  END;  { END of CASE }

  IF theta >= thetaMax THEN
        MessageDlg('Strong shock solution', mtInformation, [mbOk], 0);


  result:=theta
END;   (* ----------------------- End of Function GetThetaFromObliqueInput *)



(* ---------- Rayleigh Line ---------------------- *)

FUNCTION GetMachFromRayleighInput(const gamma : EXTENDED;
                         const input : EXTENDED;
                         const itemIndex : INTEGER) : EXTENDED;
  VAR
    gp : EXTENDED;
BEGIN
  gp:=gamma+1.0;
  result:=-1.0;
  CASE itemIndex OF
    0: IF input < 0.0 THEN    { Mach number }
         ErrorMsg('Mach number cannot be negative')
       ELSE
         result:=input;
    1: IF input < 0.0 THEN                 { Tt / Tt*  M < 1  }
         ErrorMsg('Tt / Tt* cannot be negative')
       ELSE
         IF input > 1.0 THEN
           ErrorMsg('Tt / Tt* cannot exceed 1.0')
         ELSE
           result:=RayleighMachFromTtSubsonic(input,gamma);
    2: IF input < 0.0 THEN                 { Tt / Tt*  M > 1 }
         ErrorMsg('Tt / Tt* cannot be negative')
       ELSE
         IF input > 1.0 THEN
           ErrorMsg('Tt / Tt* cannot exceed 1.0')
         ELSE
           result:=RayleighMachFromTtSupersonic(input,gamma);
    3: IF input < 0.0 THEN                      { T / T*   M < 1 }
         ErrorMsg('T/T* cannot be negative')
       ELSE
         IF input > 1.0 THEN
           ErrorMsg('T/T* cannot exceed 1')
         ELSE
           result:=RayleighMachFromTsubsonic(input, gamma);
    4: IF input < 0.0 THEN                      { T / T*   M > 1 }
         ErrorMsg('T/T* cannot be negative')
       ELSE
         IF input > 1.0 THEN
           ErrorMsg('T/T* cannot exceed 1')
         ELSE
           result:=RayleighMachFromTsupersonic(input, gamma);
    5: IF input < 0.0 THEN                    { P/ P* }
         ErrorMsg('p/p* cannot be negative')
       ELSE
         IF input > gp THEN
           ErrorMsg('p/p* cannot exceed gamma+1')
         ELSE
           result:=RayleighMachFromP(input,gamma);
    6: IF input < 1.0 THEN                    {Pt/Pt* M < 1 }
         ErrorMsg('Pt/Pt* cannot be less than 1')
       ELSE
         IF input > 1.2678 THEN
           ErrorMsg('Pt/Pt* is too large')
         ELSE
           result:=RayleighMachFromPtSubsonic(input,gamma);
    7: IF input < 0.0 THEN                    {Pt/Pt* M > 1 }
         ErrorMsg('Pt/Pt* cannot be negative')
       ELSE
         IF input > 1.0 THEN
           ErrorMsg('Pt/Pt* cannot exceed 1')
         ELSE
           result:=RayleighMachFromPtSupersonic(input,gamma);
    8: IF input < 0.0 THEN                     { V / V* }
         ErrorMsg('V/V* cannot be negative')
       ELSE
         IF input > gp/gamma THEN
           ErrorMsg('V/V* vannot exceed (gamma+1)/gamma')
         ELSE
           result:=RayleighMachFromV(input,gamma);
    9: IF input < gamma/gp THEN                     { rho / rho* }
         ErrorMsg('rho/rho* cannot be less than gamma/(gamma=1)')
       ELSE
         result:=RayleighMachFromRho(input,gamma);
  END
END;   (* ------------------- End of Function GetRayleighMach *)


(* ---------- Fanno Line ---------------------- *)

FUNCTION GetMachFromFannoInput(const gamma : EXTENDED;
                      const input : EXTENDED;
                      const itemIndex : INTEGER) : EXTENDED;
  VAR
    gp,gm : EXTENDED;     { gamma+1  and gamma-1 }
    gp2,gm2 : EXTENDED;   { one-half the above quantities }
    msq : EXTENDED;       { mach squared }
BEGIN
  gp:=gamma+1.0; gp2:=0.5*gp;
  gm:=gamma-1.0; gm2:=0.5*gm;
  result:=-1.0;
  CASE itemIndex OF
    0 : IF input < 0.0 THEN    { Mach number }
          ErrorMsg('Mach number cannot be negative')
        ELSE
          result:=input;

    1 : IF input < 0.0 THEN          { T/T* }
          ErrorMsg('T/T* cannot be negative')
        ELSE
          IF 2.0*input > gp THEN
            ErrorMsg('T/T* cannot exceed 0.5*(gamma+1)')
          ELSE
            BEGIN
              msq:=(gp/input - 2.0)/gm;
              result:=Sqrt(msq)
            END;


    2 : IF (input < 0.0) THEN                        { p/p* }
          ErrorMsg('p/p* cannot be negative')
        ELSE
          BEGIN
            msq:=(-1.0+Sqrt(1.0+gp*gm/Sqr(input)))/gm;
            result:=Sqrt(msq)
          END;

    3 : IF input < 1.0 THEN               { pt/pt*  M < 1 }
          ErrorMsg('pt/pt* cannot be less than 1')
        ELSE
          BEGIN
            result:=FannoMachFromPtPtStarSubsonic(input,gamma)
          END;

    4 : IF input < 1.0 THEN     { pt/pt*  M > 1 }
          ErrorMsg('pt/pt* cannot be less than 1')
        ELSE
          BEGIN
            result:=FannoMachFromPtPtStarSupersonic(input,gamma)
          END;

    5 : IF (input < 0.0) THEN      { V / V* }
          ErrorMsg('V/V* cannot be negative')
        ELSE
          IF Sqr(input) > gp/gm THEN
            ErrorMsg('V/V* cannot exceed Sqrt((gamma+1)/gamma-1))')
          ELSE
            result:=input/Sqrt(gp2-gm2*Sqr(input));

    6 : IF (input < 0.0) THEN     { f L / Dmax   M < 1 }
          ErrorMsg('Invalid input for fL/Dmax')
        ELSE
          BEGIN
            result:=FannoMachFromFLDsubsonic(input,gamma)
          END;

    7 : IF input < 0.0 THEN     { f L / Dmax   M > 1 }
          ErrorMsg('fL/Dmax cannot be negative')
        ELSE
          IF input > 0.8215 THEN
            ErrorMsg('fL/Dmax cannot exceed 0.8215')
          ELSE
            BEGIN
              result:=FannoMachFromFLDsupersonic(input,gamma)
            END;
  END;
END;   (* ------ End of Function GetMachFromFannoInput *)

END.
