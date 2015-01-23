UNIT NACA1135;
(* ----------------------------------------------------------------------- *)
(* PURPOSE - Compute the functions for compressible flow from the          *)
(*    NACA Report 1135                                                     *)
(* NOTES - Note that Equations 49,50,51 return the square of the velocity  *)
(*         ratio.

(* AUTHOR - Ralph L. Carmichael, Public Domain Aeronautical Software       *)

(* REVISION HISTORY                                                        *)
(*   DATE  VERS PERSON  STATEMENT OF CHANGES                               *)
(*   1997   2.0   RLC   Coded in Delphi for VuCalc                         *)
(* 30Nov00  5.5   RLC   Added Fanno and Rayliegh functions (not in 1135)   *)
(* 31Dec00  5.6   RLC   Moved Fanno and Rayleigh to their own unit         *)
(* 04May05  10.1  RLC   Added Eq80prime                                    *)

INTERFACE

CONST
  NACA1135_UNIT_VERSION = '10.1 (4May2005)';

{ Isentropic Flow - parameters: mach and gamma }
FUNCTION Eq43(const mach,gamma : EXTENDED): EXTENDED;                  {T/Tt}
FUNCTION Eq44(const mach,gamma : EXTENDED): EXTENDED;                  {p/pt}
FUNCTION Eq45(const mach,gamma : EXTENDED): EXTENDED;            { rho/rhot }
FUNCTION Eq46(const mach,gamma : EXTENDED): EXTENDED;                { a/at }
FUNCTION Eq47(const mach,gamma : EXTENDED): EXTENDED;                 { q/p }
FUNCTION Eq48(const mach,gamma : EXTENDED): EXTENDED;                { q/pt }
FUNCTION Eq49(const mach,gamma : EXTENDED): EXTENDED;           { (V/at)**2 }
FUNCTION Eq50(const mach,gamma : EXTENDED): EXTENDED;           { (V/a*)**2 }
FUNCTION Eq51(const mach,gamma : EXTENDED): EXTENDED;           { (V/Vm)**2 }
FUNCTION Eq80(const mach,gamma : EXTENDED): EXTENDED;                { A*/A }
FUNCTION Eq80prime(const mach,gamma : EXTENDED): EXTENDED;     // d(A*/A)/dM

{ Normal Shock - parameters: upstream Mach and gamma }
FUNCTION Eq93(const mach,gamma : EXTENDED): EXTENDED;                 {p2/p1}
FUNCTION Eq93inverse(const pratio,gamma : EXTENDED): EXTENDED;        {p2/p1}
FUNCTION Eq94(const mach,gamma : EXTENDED): EXTENDED;             {rho2/rho1}
FUNCTION Eq94inverse(const rhoRatio,gamma : EXTENDED): EXTENDED;      {p2/p1}
FUNCTION Eq95(const mach,gamma : EXTENDED): EXTENDED;                 {T2/T1}
FUNCTION Eq96(const mach,gamma : EXTENDED): EXTENDED;                 {M2}
FUNCTION Eq96inverse(const mach2,gamma : EXTENDED) : EXTENDED;
FUNCTION Eq97(const mach,gamma : EXTENDED): EXTENDED;                {p2/pt1}
FUNCTION Eq98(const mach,gamma : EXTENDED): EXTENDED;                {p2/pt2}
FUNCTION Eq99(const mach,gamma : EXTENDED): EXTENDED;               {pt2/pt1}
FUNCTION Eq100(const mach,gamma : EXTENDED): EXTENDED;               {pt2/p1}
FUNCTION Eq101(const mach,gamma : EXTENDED): EXTENDED;            {deltaS/cv}
FUNCTION Eq102(const mach,gamma : EXTENDED): EXTENDED;            {(p2-p1)/q}

{ Oblique Shock - Parameters: upstream Mach, wave angle, gamma }
FUNCTION Eq128(const mach1,theta,gamma : EXTENDED): EXTENDED;       { p2/p1 }
FUNCTION Eq129(const mach1,theta,gamma : EXTENDED): EXTENDED;   { rho2/rho1 }
FUNCTION Eq130(const mach1,theta,gamma : EXTENDED): EXTENDED;       { T2/T1 }
FUNCTION Eq132(const mach1,theta,gamma : EXTENDED): EXTENDED;          { M2 }
FUNCTION Eq138(const mach1,theta,gamma : EXTENDED): EXTENDED;  { cot(delta) }
FUNCTION Eq142(const mach1,theta,gamma : EXTENDED): EXTENDED;     { pt2/pt1 }
FUNCTION Eq144(const mach1,theta,gamma : EXTENDED): EXTENDED;   { deltaS/cv }

{ Prandtl-Meyer Flow }
FUNCTION Eq171c(const mach,gamma : EXTENDED) : EXTENDED;  { Prandtl-Meyer }


//FUNCTION RayleighP(const mach,gamma : EXTENDED): EXTENDED;             {P/P*}
//FUNCTION RayleighT(const mach,gamma : EXTENDED): EXTENDED;             {T/T*}
//FUNCTION RayleighTt(const mach,gamma : EXTENDED): EXTENDED;          {Tt/Tt*}
//FUNCTION RayleighRho(const mach,gamma : EXTENDED): EXTENDED;       {rho/rho*}
//FUNCTION RayleighV(const mach,gamma : EXTENDED): EXTENDED;             {V/V*}
//FUNCTION RayleighPt(const mach,gamma : EXTENDED): EXTENDED;          {Pt/Pt*}
//FUNCTION RayleighMachFromP(const p,gamma : EXTENDED): EXTENDED;
//FUNCTION RayleighMachFromRho(const rho,gamma : EXTENDED): EXTENDED;
//FUNCTION RayleighMachFromV(const v,gamma : EXTENDED): EXTENDED;
//FUNCTION RayleighMachFromTSubsonic(const t,gamma : EXTENDED): EXTENDED;
//FUNCTION RayleighMachFromTSupersonic(const t,gamma : EXTENDED): EXTENDED;
//FUNCTION RayleighMachFromTtSubsonic(const tt,gamma : EXTENDED): EXTENDED;
//FUNCTION RayleighMachFromTtSupersonic(const tt,gamma : EXTENDED): EXTENDED;
//FUNCTION RayleighMachFromPtSubsonic(const pt,gamma : EXTENDED): EXTENDED;
//FUNCTION RayleighMachFromPtSupersonic(const pt,gamma : EXTENDED): EXTENDED;

IMPLEMENTATION

USES Math;


FUNCTION Eq43(const mach,gamma : EXTENDED): EXTENDED;                  {T/Tt}
  VAR compFactor : EXTENDED;
BEGIN
  compFactor:=1.0+0.5*(gamma-1.0)*mach*mach;
  Eq43:= Power(compFactor, -1.0)
END;   (* ----------------------------------- End of NACA 1135 Function 43 *)

FUNCTION Eq44(const mach,gamma : EXTENDED): EXTENDED;                  {p/pt}
  VAR compFactor : EXTENDED;
BEGIN
  compFactor:=1.0+0.5*(gamma-1.0)*mach*mach;
  Eq44:= Power(compFactor, -gamma/(gamma-1.0))
END;   (* ----------------------------------- End of NACA 1135 Function 44 *)

FUNCTION Eq45(const mach,gamma : EXTENDED): EXTENDED;            { rho/rhot }
  VAR compFactor : EXTENDED;
BEGIN
  compFactor:=1.0+0.5*(gamma-1.0)*mach*mach;
  Eq45:= Power(compFactor, -1.0/(gamma-1.0))
END;   (* ----------------------------------- End of NACA 1135 Function 45 *)

FUNCTION Eq46(const mach,gamma : EXTENDED): EXTENDED;                { a/at }
  VAR compFactor : EXTENDED;
BEGIN
  compFactor:=1.0+0.5*(gamma-1.0)*mach*mach;
  Eq46:= Power(compFactor, -0.5)
END;   (* ----------------------------------- End of NACA 1135 Function 46 *)

FUNCTION Eq47(const mach,gamma : EXTENDED): EXTENDED;                 { q/p }
BEGIN
  Eq47:= 0.5*gamma*mach*mach
END;   (* ----------------------------------- End of NACA 1135 Function 47 *)

FUNCTION Eq48(const mach,gamma : EXTENDED): EXTENDED;                { q/pt }
  VAR compFactor : EXTENDED;
BEGIN
  compFactor:=1.0+0.5*(gamma-1.0)*mach*mach;
  Eq48:= 0.5*gamma*mach*mach*Power(compFactor, -gamma/(gamma-1.0))
END;   (* ----------------------------------- End of NACA 1135 Function 48 *)

FUNCTION Eq49(const mach,gamma : EXTENDED): EXTENDED;           { (V/at)**2 }
  VAR compFactor : EXTENDED;
BEGIN
  compFactor:=1.0+0.5*(gamma-1.0)*mach*mach;
  Eq49:= mach*mach*Power(compFactor, -1.0)
END;   (* ----------------------------------- End of NACA 1135 Function 49 *)

FUNCTION Eq50(const mach,gamma : EXTENDED): EXTENDED;           { (V/a*)**2 }
BEGIN
  Eq50:= 0.5*(gamma+1.0)*Eq49(mach,gamma)
END;   (* ----------------------------------- End of NACA 1135 Function 50 *)

FUNCTION Eq51(const mach,gamma : EXTENDED): EXTENDED;           { (V/Vm)**2 }
BEGIN
  Eq51:= 0.5*(gamma-1.0)*Eq49(mach,gamma)
END;   (* ----------------------------------- End of NACA 1135 Function 51 *)

FUNCTION Eq80(const mach,gamma : EXTENDED): EXTENDED;         { A*/A }
  VAR
    compFactor, ex : EXTENDED;
BEGIN
  ex:=0.5*(gamma+1.0)/(gamma-1.0);
  compFactor:=1.0+0.5*(gamma-1.0)*mach*mach;
  Eq80:=Power(0.5*(gamma+1.0), ex)*mach*Power(compFactor, -ex)
END;   (* ----------------------------------- End of NACA 1135 Function 80 *)

FUNCTION Eq80prime(const mach,gamma : EXTENDED): EXTENDED;      // d(A*/A)/dM
  VAR
    b,c,compFactor,d,ex,fb : EXTENDED;
BEGIN
  ex:=0.5*(gamma+1.0)/(gamma-1.0);
  b:=Power(0.5*(gamma+1.0),ex);
  c:=0.5*(gamma-1.0);
  compFactor:=1.0+c*mach*mach;
  d:=-ex;
  fb:=2.0*c*d*mach*mach*Power(compFactor,d-1.0) + Power(compFactor,d);
  Eq80prime:=b*fb
END;  // --------------------------------- End of NACA 1135 Function 80-prime

FUNCTION Eq93(const mach,gamma : EXTENDED): EXTENDED;                 {p2/p1}
BEGIN
  Eq93:=(2.0*gamma*Sqr(mach)-(gamma-1.0))/(gamma+1.0)
END;   (* ----------------------------------- End of NACA 1135 Function 93 *)

FUNCTION Eq93inverse(const pRatio,gamma : EXTENDED) : EXTENDED;
BEGIN
  Eq93inverse:=Sqrt(((gamma+1.0)*pRatio+gamma-1.0)/(2.0*gamma))
END;   (* ------------------------- End of NACA 1135 Function 93 (inverse) *)

FUNCTION Eq94(const mach,gamma : EXTENDED): EXTENDED;             {rho2/rho1}
BEGIN
  Eq94:=(gamma+1.0)*Sqr(mach)/((gamma-1.0)*Sqr(mach)+2.0)
END;   (* ----------------------------------- End of NACA 1135 Function 94 *)

FUNCTION Eq94inverse(const rhoRatio,gamma : EXTENDED) : EXTENDED;
BEGIN
  Eq94inverse:=Sqrt(2.0*rhoRatio/(gamma+1.0-(gamma-1.0)*rhoRatio))
END;   (* ------------------------- End of NACA 1135 Function 94 (inverse) *)

FUNCTION Eq95(const mach,gamma : EXTENDED): EXTENDED;                 {T2/T1}
BEGIN
  Eq95:=Eq93(mach,gamma)/Eq94(mach,gamma)
END;   (* ----------------------------------- End of NACA 1135 Function 95 *)

FUNCTION Eq96(const mach,gamma : EXTENDED): EXTENDED;                    {M2}
  VAR gm1 : EXTENDED;
BEGIN
  gm1:=gamma-1.0;
  Eq96:=(gm1*mach*mach+2.0)/(2.0*gamma*mach*mach-gm1);
END;   (* ----------------------------------- End of NACA 1135 Function 96 *)

FUNCTION Eq96inverse(const mach2,gamma : EXTENDED) : EXTENDED;
  VAR gm1,m2sq : EXTENDED;
BEGIN
  gm1:=gamma-1.0;
  m2sq:=mach2*mach2;
  Eq96inverse:=Sqrt((2.0+gm1*m2sq)/(2.0*gamma*m2sq-gm1))
END;   (* --------------------------- End of NACA 1135 Function 96 inverse *)

FUNCTION Eq97(const mach,gamma : EXTENDED): EXTENDED;                {p2/pt1}
  VAR
    gm1,gp1,ex : EXTENDED;
BEGIN
  gm1:=gamma-1.0;
  gp1:=gamma+1.0;
  ex:=gamma/gm1;
  Eq97:=((2.0*gamma*mach*mach-gm1)/gp1)*Power(2.0/(gm1*mach*mach+2.0), ex)
END;   (* ----------------------------------- End of NACA 1135 Function 97 *)

FUNCTION Eq98(const mach,gamma : EXTENDED): EXTENDED;                {p2/pt2}
  VAR
    gm1,gp1,ex,msq : EXTENDED;
BEGIN
  gm1:=gamma-1.0;
  gp1:=gamma+1.0;
  ex:=gamma/gm1;
  msq:=mach*mach;
  Eq98:=Power((4.0*gamma*msq-2.0*gm1)/(gp1*gp1*msq), ex)
END;   (* ----------------------------------- End of NACA 1135 Function 98 *)

FUNCTION Eq99(const mach,gamma : EXTENDED): EXTENDED;               {pt2/pt1}
BEGIN
  Eq99:=Power(Eq94(mach,gamma), gamma/(gamma-1.0))*
        Power(Eq93(mach,gamma), -1.0/(gamma-1.0));
END;   (* ----------------------------------- End of NACA 1135 Function 99 *)

FUNCTION Eq100(const mach,gamma : EXTENDED): EXTENDED;               {pt2/p1}
BEGIN
  Eq100:=1.0;  { not coded yet }
END;   (* ---------------------------------- End of NACA 1135 Function 100 *)


FUNCTION Eq101(const mach,gamma : EXTENDED): EXTENDED;            {deltaS/cv}
BEGIN
  Eq101:= -(gamma-1.0)*Eq99(mach,gamma)
END;   (* ---------------------------------- End of NACA 1135 Function 101 *)

FUNCTION Eq102(const mach,gamma : EXTENDED): EXTENDED;            {(p2-p1)/q}
BEGIN
  Eq102:= 4.0*(mach*mach-1.0)/((gamma+1.0)*mach*mach)
END;   (* ---------------------------------- End of NACA 1135 Function 102 *)





{ Oblique Shock - Parameters: upstream Mach, wave angle, and gamma }
{  wave angle, theta and wedge angle, delta, are always in radians }

FUNCTION Eq128(const mach1,theta,gamma : EXTENDED): EXTENDED;     { pt2/pt1 }
BEGIN
  Eq128:=Eq93(mach1*Sin(theta), gamma)
END;   (* ---------------------------------- End of NACA 1135 Function 128 *)

FUNCTION Eq129(const mach1,theta,gamma : EXTENDED): EXTENDED;   { rho2/rho1 }
BEGIN
  Eq129:= Eq94(mach1*Sin(theta), gamma)
END;   (* ---------------------------------- End of NACA 1135 Function 129 *)

FUNCTION Eq130(const mach1,theta,gamma : EXTENDED): EXTENDED;       { T2/T1 }
BEGIN
  Eq130:=Eq95(mach1*Sin(theta), gamma)
END;   (* ---------------------------------- End of NACA 1135 Function 130 *)

FUNCTION Eq131(const mach1,theta,gamma : EXTENDED): EXTENDED;
BEGIN
  Eq131:=Eq96(mach1*Sin(theta), gamma)
END;   (* ---------------------------------- End of NACA 1135 Function 131 *)

FUNCTION Eq132(const mach1,theta,gamma : EXTENDED): EXTENDED;          { M2 }
  VAR
    num,denom : EXTENDED;
    gm1,ms : EXTENDED;
BEGIN
  gm1:=gamma-1.0;
  ms:=Sqr(mach1*Sin(theta));
  num:=Sqr(gamma+1.0)*ms*Sqr(mach1) - 4.0*(ms-1.0)*(gamma*ms+1.0);
  denom:=(2.0*gamma*ms-gm1)*(gm1*ms+2.0);
  Eq132:=Sqrt(num/denom)
END;   (* ---------------------------------- End of NACA 1135 Function 132 *)

FUNCTION Eq132A(const mach1,sinth,gamma : EXTENDED): EXTENDED;          { M2 }
  VAR
    num,denom : EXTENDED;
    gm1,ms : EXTENDED;
BEGIN
  gm1:=gamma-1.0;
  ms:=Sqr(mach1*sinth);
  num:=Sqr(gamma+1.0)*ms*Sqr(mach1) - 4.0*(ms-1.0)*(gamma*ms+1.0);
  denom:=(2.0*gamma*ms-gm1)*(gm1*ms+2.0);
  Eq132A:=Sqrt(num/denom)
END;   (* --------------------------------- End of NACA 1135 Function 132A *)

FUNCTION Eq138(const mach1,theta,gamma : EXTENDED): EXTENDED;  { cot(delta) }
  VAR
//    num,denom : EXTENDED;
    gp1,ms : EXTENDED;
BEGIN
{   cotd = tan(shkang*convdr)*((gp1*mach1s)/                    /* EQ. 138 */
            (2.0*(mach1s * sints - 1.0))-1.0);}

//  gm1:=gamma-1.0;
  gp1:=gamma+1.0;
  ms:=Sqr(mach1*Sin(theta));
  IF ms <= 1.0 THEN
    Eq138:=1E38    { something really huge }
  ELSE
    BEGIN
//      num:=Sqr(gamma+1.0)*ms*Sqr(mach1) - 4.0*(ms-1.0)*(gamma*ms+1.0);
//      denom:=(2.0*gamma*ms-gm1)*(gm1*ms+2.0);
      Eq138:=Tan(theta)*(gp1*Sqr(mach1)/(2.0*(ms-1.0))-1.0);
    END;
END;   (* ---------------------------------- End of NACA 1135 Function 138 *)

FUNCTION Eq140(const mach1,theta,gamma : EXTENDED) : EXTENDED;
BEGIN
  Eq140:=Eq97(mach1*Sin(theta), gamma)
END;   (* ---------------------------------- End of NACA 1135 Function 140 *)

FUNCTION Eq142(const mach1,theta,gamma : EXTENDED): EXTENDED;     { pt2/pt1 }
  VAR
//    num,denom : EXTENDED;
    factor1,factor2 : EXTENDED;
    gm1,gp1,ms : EXTENDED;
BEGIN

{   ptrat = (pow(((gp1*mach1s*sints)/(gm1*mach1s*sints + 2.0)),(gama/gm1)))
       * pow((gp1/(2.0*gama*mach1s*sints - gm1)),(1.0/gm1)) ;   /* EQ. 142 */}

  gm1:=gamma-1.0;
  gp1:=gamma+1.0;
  ms:=Sqr(mach1*Sin(theta));
  factor1:=Power((gp1*ms)/(gm1*ms+2.0), (gamma/gm1) );
  factor2:=Power(gp1/(2.0*gamma*ms-gm1), 1.0/gm1);
  Eq142:=factor1*factor2
END;   (* ---------------------------------- End of NACA 1135 Function 142 *)

FUNCTION Eq144(const mach1,theta,gamma : EXTENDED): EXTENDED;   { deltaS/cv }
BEGIN
  Eq144:= -(gamma-1.0)*Eq142(mach1,theta,gamma)
END;   (* ---------------------------------- End of NACA 1135 Function 144 *)

FUNCTION Eq171c(const mach,gamma : EXTENDED) : EXTENDED;    { Prandtl-Meyer }
  VAR
    g : EXTENDED;
BEGIN
  IF mach > 1.0 THEN
    BEGIN
      g:=(gamma+1.0)/(gamma-1.0);
      Eq171c:=Sqrt(g)*ArcTan(Sqrt((mach*mach-1.0)/g)) -
                      ArcTan(Sqrt(mach*mach-1.0))
    END
  ELSE
    Eq171c:=0.0
END;   (* --------------------------------- End of NACA 1135 Function 171c *)


{
FUNCTION Eqxx(const mach,gamma : EXTENDED): EXTENDED;
  VAR
    compFactor : EXTENDED;
BEGIN
  compFactor:=1.0+0.5*(gamma-1.0)*mach*mach;
  eqxx:=1.0
END; }  (* --------------------------------- End of NACA 1135 Function xxxx *)

    (* -------------------------------------------- End of Unit NACA1135 *)
END.
