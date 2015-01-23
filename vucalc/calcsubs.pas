UNIT CalcSubs;
(* PURPOSE - Compute various flow quantities from inputs *)
(* AUTHORS - Tom Benson, NASA Glenn Research Center (formerly Lewis)       *)
(*           Ralph L. Carmichael, Public Domain Aeronautical Software      *)
(* REVISION HISTORY                                                        *)
(*   DATE  VERS PERSON  STATEMENT OF CHANGES                               *)
(*    ?     1.0    TB   Original coding (in C for Unix, X-windows, Forms)  *)
(*   1997   2.0   RLC   Recoded in Delphi for Microsoft Windows            *)
(* 21Nov00  6.0   RLC   Gathered all the procedures into a unit            *)
(* 04May05 10.1   RLC   Revised GetMachFromAreaRatio (use Newton's method) *)

(* NOTE - should probably fold this into NACA1135.pas soon *)

INTERFACE

CONST
  CALCSUBS_UNIT_VERSION = '0.5 (22Nov00)';

(* Get the Mach number from the area ratio A/A*  *)
FUNCTION GetMachFromAreaRatio(const aRatio,gamma : EXTENDED;
                              const super : BOOLEAN) : EXTENDED;
(* get the Mach number, given the Prandtl-Meyer angle *)
FUNCTION GetMachFromPM(const nu,gamma : EXTENDED) : EXTENDED;

(* ------ Inverse Functions for Normal Shock ----------------- *)
FUNCTION GetMachFromPtRatio(const inputValue,gamma : EXTENDED) : EXTENDED;

FUNCTION GetMachFromTempRatio(const inputValue,gamma : EXTENDED) : EXTENDED;

(* ------ Inverse Functions for Oblique Shock -------------------- *)
FUNCTION GetShockAngleFromRampAndMach(const delta,mach,gamma : EXTENDED)
   : EXTENDED;

PROCEDURE MaxWedgeAngleForAttachedShock(const mach1,gamma : EXTENDED;
                                        var delta,theta : EXTENDED);
FUNCTION GetThetaFromM1AndM2(const mach1,mach2,gamma : EXTENDED) : EXTENDED;

(* ------ Inverse Functions for Fanno Flow ---------- *)
FUNCTION GetMachFromFLDsubsonic(
    CONST fld, gamma : EXTENDED): EXTENDED;

FUNCTION GetMachFromFLDsupersonic(
    CONST fld, gamma : EXTENDED): EXTENDED;

FUNCTION GetMachFromPtPtStarSubsonic(
    CONST ptpts, gamma : EXTENDED): EXTENDED;

FUNCTION GetMachFromPtPtStarSupersonic(
    CONST ptpts, gamma : EXTENDED): EXTENDED;

IMPLEMENTATION

Uses Math,NACA1135,Dialogs;


PROCEDURE ErrorMsg(CONST s : STRING);
BEGIN
  MessageDlg(s, mtError, [mbOK], 0)
END;   (* ------------------------ End of Procedure ErrorMsg *)

(* Get the Mach number from the area ratio A/A*  *)
FUNCTION GetMachFromAreaRatio(const aRatio,gamma : EXTENDED;
                              const super : BOOLEAN) : EXTENDED;
  CONST
    EPS = 1E-8;
  VAR
//    factor1,compFactor : EXTENDED;
//    machOld,machNew, arOld,arNew, deriv : EXTENDED;
    a,aprime,ar,mach : EXTENDED;
BEGIN
  IF aRatio < 1.0 THEN BEGIN GetMachFromAreaRatio:=0.0; Exit END;
//  factor1:=0.5*(gamma+1.0)/(gamma-1.0);
  ar:=1.0/aratio;  // ar is A*/A, but aratio is A/A*
// Use Newton's method to find the mach such that Eq80(mach,gamma)=ar

  IF super THEN mach:=2 ELSE mach:=ar/1.728;

  REPEAT
    a:=Eq80(mach,gamma);
    aprime:=Eq80prime(mach,gamma);
    mach:=mach-(a-ar)/aprime;
  UNTIL Abs(a-ar) < EPS;   { test should be on machs }
  GetMachFromAreaRatio:=mach
END;   // ------------------------------ End of Function GetMachFromAreaRatio


(* get the Mach number, given the Prandtl-Meyer angle, in degrees *)
FUNCTION GetMachFromPM(const nu,gamma : EXTENDED) : EXTENDED;
  CONST
    EPS = 1E-6;
  VAR
    nurad : EXTENDED;
    machOld,machNew, nuOld,nuNew, deriv: EXTENDED;
BEGIN
  nurad:=DegToRad(nu);
{  WriteLn('GetMachFromPM for nu=', nurad:10:6, nu:10:6); }
 { g1:=(gamma+1.0)/(gamma-1.0);}
  machOld:=2.0; machNew:=2.1;
  nuOld:=Eq171c(machOld,gamma);
  REPEAT
    nuNew:=Eq171c(machNew,gamma);
    deriv:=(nuNew-nuOld)/(machNew-machOld);
    nuOld:=nuNew;
    machOld:=machNew;
    machNew:=machOld + (nurad-nuOld)/deriv;
    IF machNew < 1.0 THEN machNew:=1.0;
{    WriteLn(betaOld:8:4, betaNew:8:4, nuOld:8:4, nuNew:8:4, deriv:8:4);}
  UNTIL Abs(nurad-nuOld) < EPS;
  GetMachFromPM:=machOld
END;   (* ---------------------------------- End of Function GetMachFromPM *)


(* ------ Inverse Functions for Normal Shock ----------------- *)

(*   --- Inverse of Eq 99 ---------- *)
FUNCTION GetMachFromPtRatio(const inputValue,gamma : EXTENDED) : EXTENDED;
  CONST
    EPS = 1E-4;
  VAR
    msOld,msNew : EXTENDED;  { iterated values of mach-squared }
    pOld,pNew : EXTENDED; { iterated values of total pressure ratio }
    deriv : EXTENDED;
BEGIN
  msOld:=2.25;
  pOld:=Eq99(msOld,gamma);
  msNew:=msOld+0.2;
  REPEAT
    pNew:=Eq99(msNew,gamma);
    deriv:=(pNew-pOld)/(msNew-msOld);
    pOld:=pNew;
    msOld:=msNew;
    msNew:=msOld+(inputValue-pOld)/deriv;
  UNTIL Abs(pOld-inputValue) < EPS;
  GetMachFromPtRatio:=msOld
END;   (* ----------------------------- End of Function GetMachFromPtRatio *)

(* ---------- Inverse of Eq 95 -------------------------- *)
FUNCTION GetMachFromTempRatio(const inputValue,gamma : EXTENDED) : EXTENDED;
  CONST
    EPS = 1E-4;
  VAR
    msOld,msNew : EXTENDED;  { iterated values of mach-squared }
    tOld,tNew : EXTENDED; { iterated values of total pressure ratio }
    deriv : EXTENDED;
BEGIN
  msOld:=2.25;
  tOld:=Eq95(msOld,gamma);
  msNew:=msOld+0.1;
  REPEAT
    tNew:=Eq95(msNew,gamma);
    deriv:=(tNew-tOld)/(msNew-msOld);
    tOld:=tNew;
    msOld:=msNew;
    msNew:=msOld+(inputValue-tOld)/deriv;
  UNTIL Abs(tOld-inputValue) < EPS;
  GetMachFromTempRatio:=msOld
END;   (* --------------------------- End of Function GetMachFromTempRatio *)

FUNCTION GetShockAngleFromRampAndMach(const delta,mach,gamma : EXTENDED)
   : EXTENDED;
  CONST
    EPS = 1E-6;
  VAR
    thetaOld,thetaNew, deltaOld,deltaNew, deriv : EXTENDED;
    cotd : EXTENDED;
BEGIN
  thetaOld:=ArcSin(1.0/mach);
  deltaOld:=0.0;
  thetaNew:=thetaOld + 0.05;
  REPEAT
    cotd:=Eq138(mach,thetaNew,gamma);
    deltaNew:=ArcTan(1.0/cotd);
    deriv:=(deltaNew-deltaOld)/(thetaNew-thetaOld);
    deltaOld:=deltaNew;
    thetaOld:=thetaNew;
    thetaNew:=thetaOld + (delta-deltaOld)/deriv;
  UNTIL Abs(delta-deltaOld) < EPS;
  GetShockAngleFromRampAndMach:=thetaOld
END;


FUNCTION SearchForMaxDelta(ax,bx,cx : EXTENDED;
                           tol : EXTENDED;
                           mach,gamma : EXTENDED;
                           VAR theta : EXTENDED) : EXTENDED;
  CONST
    R = 0.61803399;
    C = 1.0-R;
  VAR
    x0,x1,x2,x3 : EXTENDED;
    f1,f2 : EXTENDED;

BEGIN
  x0:=ax;
  x3:=cx;
  IF Abs(cx-bx) > Abs(bx-ax) THEN
    BEGIN
      x1:=bx;
      x2:=bx+C*(cx-bx)
    END
  ELSE
    BEGIN
      x2:=bx;
      x1:=bx-C*(bx-ax)
    END;
  f1:=-ArcTan(1.0/Eq138(mach, x1, gamma));
  f2:=-ArcTan(1.0/Eq138(mach, x2, gamma));

  WHILE Abs(x3-x0) > tol*(Abs(f1)+Abs(f2)) DO
    BEGIN
       IF f2 < f1 THEN
        BEGIN
          x0:=x1;
          x1:=x2;
          x2:=R*x1+C*x3;
          f1:=f2;
          f2:=-ArcTan(1.0/Eq138(mach,x2,gamma));
        END
      ELSE
        BEGIN
          x3:=x2;
          x2:=x1;
          x1:=R*x2+C*x0;
          f2:=f1;
          f1:=-ArcTan(1.0/Eq138(mach,x1,gamma));
        END;
    END;

  IF f1 < f2 THEN
    BEGIN
      theta:=x1;
      SearchForMaxDelta:=-f1
    END
  ELSE
    BEGIN
      theta:=x2;
      SearchForMaxDelta:=-f2
    END;
END;


PROCEDURE MaxWedgeAngleForAttachedShock(const mach1,gamma : EXTENDED;
                                        var delta,theta : EXTENDED);
  VAR
    ax,bx,cx, tol: EXTENDED;
BEGIN
  ax:=60*PI/180;
  bx:=65*PI/180;
  cx:=85*PI/180;
  tol:=1E-4;
  delta:=SearchForMaxDelta(ax,bx,cx, tol, mach1,gamma, theta)
END;   (* ------------------------- End of Procedure MaxAttachedShockAngle *)

(* ------------------ Oblique Shock -------------------------------- *)

FUNCTION GetThetaFromM1AndM2(const mach1,mach2,gamma : EXTENDED) : EXTENDED;
  CONST
    EPS = 1E-4;
  VAR
//    mach1s,mach1f : EXTENDED;
    thetaOld,thetaNew,mach2Old,mach2New,deriv : EXTENDED;

BEGIN
//  mach1s:=mach1*mach1;
//  mach1f:=mach1s*mach1s;

  thetaOld:=ArcSin(1.0/mach1);
  mach2Old:=mach1;
  thetaNew:= thetaOld + 0.05;
  REPEAT
    mach2New:= Eq132(mach1, thetaNew, gamma);
    deriv:= (mach2New-mach2Old)/(thetaNew-thetaOld);
    mach2Old:=mach2New;
    thetaOld:=thetaNew;
    thetaNew:=thetaOld + (mach2-mach2Old)/deriv;
  UNTIL Abs(mach2-mach2Old) < EPS;

  GetThetaFromM1AndM2:=thetaOld
END;

(* ------------ Fanno Line Inverse Functions ------------- *)


FUNCTION GetMachFromFLDsubsonic(
    CONST fld, gamma : EXTENDED): EXTENDED;
  VAR
    gp,gm : EXTENDED;
    mi,mo : EXTENDED;
    msub : EXTENDED;
    p1,p2 : EXTENDED;
BEGIN
//  result:=-1.0;
  gp:=gamma+1.0;
  gm:=gamma-1.0;

  mi:=1;//  mo:=0.5;

  REPEAT
    p1:=0.5*gp*Sqr(mi)/(1.0+0.5*gm*Sqr(mi));
    p2:=gamma*fld+1.0-0.5*gp*Ln(p1);    msub:=Sqrt(1.0/p2);    mo:=mi;
    mi:=msub
  UNTIL Abs(mi-mo) < 1E-10;
  result:=mi
END;   (* ----- End of Function GetMachFrom FLDsubsonic *)


FUNCTION GetMachFromFLDsupersonic(
    CONST fld, gamma : EXTENDED): EXTENDED;
  VAR
    gp,gm : EXTENDED;
    fldMax : EXTENDED;
    mi,mo : EXTENDED;
    msuper : EXTENDED;
    p1,p2 : EXTENDED;
BEGIN
  result:=-1.0;
  gp:=gamma+1.0;
  gm:=gamma-1.0;

  fldmax:=gp*Ln(gp/gm)/(2.0*gamma)-1.0/gamma;  if fld > fldmax THEN
    BEGIN
      ErrorMsg('fL/D exceeds max. value');
      Exit
    END;

  mi:=1;//  mo:=5;

  REPEAT
    p1:= -2.0*(gamma*fld+1.0-1.0/Sqr(mi))/gp;
    p2:=0.5*gp*Exp(p1)-0.5*gm;    msuper:=Sqrt(1.0/p2);    mo:=mi;
    mi:=msuper
  UNTIL Abs(mi-mo) < 1E-10;
  result:=mi
END;   (* ----- End of Function GetMachFrom FLDsupersonic *)



FUNCTION GetMachFromPtPtStarSubsonic(
    CONST ptpts, gamma : EXTENDED): EXTENDED;
  VAR
    gp,gm : EXTENDED;
    mi,mo : EXTENDED;
    m : EXTENDED;
    p1,p2 : EXTENDED;
BEGIN
//  result:=-1.0;
  gp:=gamma+1.0;
  gm:=gamma-1.0;

  mi:=5.0;//  mo:=1.0;

  REPEAT
    p1:=Power(mi*ptpts, 2.0*gm/gp);
    p2:=2.0*(0.5*p1*gp-1.0)/gm;
    m:=Sqrt(p2);    mo:=mi;
    mi:=m
  UNTIL Abs(mi-mo) < 1E-10;
  result:=mi
END;   (* ----- End of Function GetMachFrom PtPtStarSubsonic *)


FUNCTION GetMachFromPtPtStarSupersonic(
    CONST ptpts, gamma : EXTENDED): EXTENDED;
  VAR
    gp,gm : EXTENDED;
    mi,mo : EXTENDED;
    m : EXTENDED;
    p1,p2 : EXTENDED;
BEGIN
//  result:=-1.0;
  gp:=gamma+1.0;
  gm:=gamma-1.0;

  mi:=5.0;//  mo:=1.0;

  REPEAT
    p1:=Power(mi*ptpts, 2.0*gm/gp);
    p2:=2.0*(0.5*p1*gp-1.0)/gm;
    m:=Sqrt(Abs(p2));    mo:=mi;
    mi:=m
  UNTIL Abs(mi-mo) < 1E-10;
  result:=mi
END;   (* ----- End of Function GetMachFrom PtPtStarSupersonic *)


(* ------------------------------------ End of Unit CalcSubs. *)
END.
