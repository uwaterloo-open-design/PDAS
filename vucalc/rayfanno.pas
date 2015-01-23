UNIT RayFanno;
(* ----------------------------------------------------------------------- *)
(* PURPOSE - Compute the functions for Rayleigh flow                       *)

(* AUTHOR - Ralph L. Carmichael, Public Domain Aeronautical Software       *)
(*          thanks to Bill Mason, David Pratt, xxxxxxx                     *)

(* REVISION HISTORY                                                        *)
(*   DATE  VERS PERSON  STATEMENT OF CHANGES                               *)
(* 30Nov00  5.5   RLC   Original coding                                    *)

INTERFACE

CONST
  RAYFANNO_UNIT_VERSION = '5.6 (10 Dec 2000)';


PROCEDURE AllRayleigh(CONST mach,gamma : EXTENDED;
  VAR tr,pr,tq,vr,pq,rr : EXTENDED);
FUNCTION RayleighP(const mach,gamma : EXTENDED): EXTENDED;             {P/P*}
FUNCTION RayleighT(const mach,gamma : EXTENDED): EXTENDED;             {T/T*}
FUNCTION RayleighTt(const mach,gamma : EXTENDED): EXTENDED;          {Tt/Tt*}
FUNCTION RayleighRho(const mach,gamma : EXTENDED): EXTENDED;       {rho/rho*}
FUNCTION RayleighV(const mach,gamma : EXTENDED): EXTENDED;             {V/V*}
FUNCTION RayleighPt(const mach,gamma : EXTENDED): EXTENDED;          {Pt/Pt*}

FUNCTION RayleighMachFromP(const p,gamma : EXTENDED): EXTENDED;
FUNCTION RayleighMachFromRho(const rho,gamma : EXTENDED): EXTENDED;
FUNCTION RayleighMachFromV(const v,gamma : EXTENDED): EXTENDED;
FUNCTION RayleighMachFromTSubsonic(const t,gamma : EXTENDED): EXTENDED;
FUNCTION RayleighMachFromTSupersonic(const t,gamma : EXTENDED): EXTENDED;
FUNCTION RayleighMachFromTtSubsonic(const tt,gamma : EXTENDED): EXTENDED;
FUNCTION RayleighMachFromTtSupersonic(const tt,gamma : EXTENDED): EXTENDED;
FUNCTION RayleighMachFromPtSubsonic(const pt,gamma : EXTENDED): EXTENDED;
FUNCTION RayleighMachFromPtSupersonic(const pt,gamma : EXTENDED): EXTENDED;

PROCEDURE AllFanno(CONST mach,gamma : EXTENDED;
   VAR tf,pf,pz,ff,vf,fd : EXTENDED);
FUNCTION FannoP(const mach,gamma : EXTENDED): EXTENDED;             {P/P*}
FUNCTION FannoT(const mach,gamma : EXTENDED): EXTENDED;             {T/T*}
FUNCTION FannoRho(const mach,gamma : EXTENDED): EXTENDED;       {rho/rho*}
FUNCTION FannoV(const mach,gamma : EXTENDED): EXTENDED;             {V/V*}
FUNCTION FannoPt(const mach,gamma : EXTENDED): EXTENDED;          {Pt/Pt*}
FUNCTION FannoFLD(const mach,gamma : EXTENDED) : EXTENDED;

FUNCTION FannoMachFromFLDsubsonic(CONST fld, gamma : EXTENDED): EXTENDED;
FUNCTION FannoMachFromFLDsupersonic(CONST fld, gamma : EXTENDED): EXTENDED;
FUNCTION FannoMachFromPtPtStarSubsonic(CONST ptpts, gamma : EXTENDED): EXTENDED;
FUNCTION FannoMachFromPtPtStarSupersonic(CONST ptpts, gamma : EXTENDED): EXTENDED;


IMPLEMENTATION

USES Dialogs,Math,NACA1135;

PROCEDURE ErrorMsg(const s : STRING);
BEGIN
  MessageDlg(s, mtError, [mbOk], 0)
END;   (* -------------------------------------- End of Procedure ErrorMsg *)


FUNCTION RayleighP(const mach,gamma : EXTENDED): EXTENDED;          {Pt/Pt*}
BEGIN
  result:=(gamma+1.0)/(1.0+gamma*mach*mach)
END;  (* ------------------------------ End of Function RayleighP *)

FUNCTION RayleighT(const mach,gamma : EXTENDED): EXTENDED;             {T/T*}
  VAR
    msq : EXTENDED;
BEGIN
  msq:=Sqr(mach);
  result:=Sqr(gamma+1.0)*msq/Sqr(1.0+gamma*msq)
END;  (* ------------------------------ End of Function RayleighT *)

FUNCTION RayleighTt(const mach,gamma : EXTENDED): EXTENDED;          {Tt/Tt*}
  VAR
    denom,msq,num : EXTENDED;
BEGIN
  msq:=Sqr(mach);
  num:=2.0*(gamma+1.0)*msq*(1.0+0.5*(gamma-1.0)*msq);
  denom:=Sqr(1.0+gamma*msq);
  result:=num/denom
END;  (* ------------------------------ End of Function RayleighTt *)

FUNCTION RayleighRho(const mach,gamma : EXTENDED): EXTENDED;       {rho/rho*}
  VAR
    msq : EXTENDED;
BEGIN
  msq:=Sqr(mach);
  result:=(1.0+gamma*msq)/((gamma+1.0)*msq)
END;  (* ---------------------------- End of Function RayleighRho *)

FUNCTION RayleighV(const mach,gamma : EXTENDED): EXTENDED;             {V/V*}
  VAR
    msq : EXTENDED;
BEGIN
  msq:=Sqr(mach);
  result:=(1.0+gamma)*msq/(1.0+gamma*msq)
END;  (* ---------------------------- End of Function RayleighV *)

FUNCTION RayleighPt(const mach,gamma : EXTENDED): EXTENDED;          {Pt/Pt*}
  VAR
    gp,gm : EXTENDED;
    msq : EXTENDED;
    pr,tt : EXTENDED;
BEGIN
  gp:=gamma+1.0; gm:=gamma-1.0; msq:=mach*mach;
  tt:=1.0/(1.0+0.5*gm*msq);
  pr:=gp/(1.0+gamma*msq);
  result:=pr*Power(2.0/(gp*tt), gamma/gm)
END;  (* ----------------------------- End of Function RayleighPt *)

(* ------- Inverse Rayleigh Functions ----------------- *)

PROCEDURE AllRayleigh(CONST mach,gamma : EXTENDED;
  VAR tr,pr,tq,vr,pq,rr : EXTENDED);

  VAR
    gp,gm : EXTENDED;
    msq : EXTENDED;
    tt : EXTENDED;
BEGIN
  gp:=gamma+1.0; gm:=gamma-1.0; msq:=mach*mach;
  tt:=Eq43(mach,gamma);
//  tf:=0.5*gp*tt;
  pr:=gp/(1.0+gamma*msq);
  tr:=pr*pr*msq;
  tq:=2.0*tr/(gp*tt);
  vr:=pr*msq;
  pq:=pr*Power(2.0/(gp*tt), gamma/gm);
  rr:=RayleighRho(mach,gamma);
END;

FUNCTION RayleighMachFromP(const p,gamma : EXTENDED): EXTENDED;
  VAR
    msq : EXTENDED;
BEGIN
  msq:=(gamma+1.0-p)/(gamma*p);
  result:=Sqrt(msq)
END;  (* ---------------------- End of Function RayleighMachFromP *)

FUNCTION RayleighMachFromRho(const rho,gamma : EXTENDED): EXTENDED;
  VAR
    msq : EXTENDED;
BEGIN
  msq:=1.0/(rho*gamma+rho-gamma);
  result:=Sqrt(msq)
END;  (* -------------------- End of Function RayleighMachFromRho *)

FUNCTION RayleighMachFromV(const v,gamma : EXTENDED): EXTENDED;
  VAR
    msq : EXTENDED;
BEGIN
  msq:=v/(gamma+1.0-gamma*v);
  result:=Sqrt(msq)
END;  (* -------------------- End of Function RayleighMachFromV *)

FUNCTION RayleighMachFromTSubsonic(const t,gamma : EXTENDED): EXTENDED;
  CONST
    EPSILON = 1.0E-10;
  VAR
    m1,m2 : EXTENDED;
    mnew : EXTENDED;
    t1,t2 : EXTENDED;
    slope : EXTENDED;
BEGIN
  m1:=1.0; t1:=1.0-t;
  m2:=0.5; t2:=RayleighT(m2,gamma)-t;
  REPEAT
    slope:=(m2-m1)/(t2-t1);
    mnew:=m2-slope*t2;
    m1:=m2; m2:=mnew;
    t1:=t2; t2:=RayleighT(m2,gamma)-t
  UNTIL Abs(m2-m1) < EPSILON;

  result:=m2
END;   (* ----------- End of Function RayleighMachFromTSubsonic *)

FUNCTION RayleighMachFromTSupersonic(const t,gamma : EXTENDED): EXTENDED;
  CONST
    EPSILON = 1.0E-10;
  VAR
    m1,m2 : EXTENDED;
    mnew : EXTENDED;
    t1,t2 : EXTENDED;
    slope : EXTENDED;
BEGIN
  m1:=1.0; t1:=1.0-t;
  m2:=3.0; t2:=RayleighT(m2,gamma)-t;
  REPEAT
    slope:=(m2-m1)/(t2-t1);
    mnew:=m2-slope*t2;
    m1:=m2; m2:=mnew;
    t1:=t2; t2:=RayleighT(m2,gamma)-t
  UNTIL Abs(m2-m1) < EPSILON;

  result:=m2
END;   (* --------- End of Function RayleighMachFromTSupersonic *)

FUNCTION RayleighMachFromTtSubsonic(const tt,gamma : EXTENDED): EXTENDED;
  CONST
    EPSILON = 1.0E-10;
  VAR
    m1,m2 : EXTENDED;
    mnew : EXTENDED;
    t1,t2 : EXTENDED;
    slope : EXTENDED;
BEGIN
  m1:=1.0; t1:=1.0-tt;
  m2:=0.5; t2:=RayleighTt(m2,gamma)-tt;
  REPEAT
    slope:=(m2-m1)/(t2-t1);
    mnew:=m2-slope*t2;
    m1:=m2; m2:=mnew;
    t1:=t2; t2:=RayleighTt(m2,gamma)-tt
  UNTIL Abs(m2-m1) < EPSILON;

  result:=m2
END;   (* ----------- End of Function RayleighMachFromTtSubsonic *)

FUNCTION RayleighMachFromTtSupersonic(const tt,gamma : EXTENDED): EXTENDED;
  CONST
    EPSILON = 1.0E-10;
  VAR
    m1,m2 : EXTENDED;
    mnew : EXTENDED;
    t1,t2 : EXTENDED;
    slope : EXTENDED;
BEGIN
  m1:=1.0; t1:=1.0-tt;
  m2:=3.0; t2:=RayleighTt(m2,gamma)-tt;
  REPEAT
    slope:=(m2-m1)/(t2-t1);
    mnew:=m2-slope*t2;
    m1:=m2; m2:=mnew;
    t1:=t2; t2:=RayleighTt(m2,gamma)-tt
  UNTIL Abs(m2-m1) < EPSILON;

  result:=m2
END;   (* --------- End of Function RayleighMachFromTtSupersonic *)



FUNCTION RayleighMachFromPtSubsonic(const pt,gamma : EXTENDED): EXTENDED;
  CONST
    EPSILON = 1.0E-10;
  VAR
    m1,m2 : EXTENDED;
    mnew : EXTENDED;
    p1,p2 : EXTENDED;
    slope : EXTENDED;
BEGIN
  m1:=1.0; p1:=1.0-pt;
  m2:=0.5; p2:=RayleighPt(m2,gamma)-pt;
  REPEAT
    slope:=(m2-m1)/(p2-p1);
    mnew:=m2-slope*p2;
    m1:=m2; m2:=mnew;
    p1:=p2; p2:=RayleighPt(m2,gamma)-pt
  UNTIL Abs(m2-m1) < EPSILON;

  result:=m2
END;   (* ----------- End of Function RayleighMachFromPtSubsonic *)

FUNCTION RayleighMachFromPtSupersonic(const pt,gamma : EXTENDED): EXTENDED;
  CONST
    EPSILON = 1.0E-10;
  VAR
    m1,m2 : EXTENDED;
    mnew : EXTENDED;
    p1,p2 : EXTENDED;
    slope : EXTENDED;
BEGIN
  m1:=1.0; p1:=1.0-pt;
  m2:=3.0; p2:=RayleighPt(m2,gamma)-pt;
  REPEAT
    slope:=(m2-m1)/(p2-p1);
    mnew:=m2-slope*p2;
    m1:=m2; m2:=mnew;
    p1:=p2; p2:=RayleighPt(m2,gamma)-pt
  UNTIL Abs(m2-m1) < EPSILON;

  result:=m2
END;   (* --------- End of Function RayleighMachFromPtSupersonic *)

(* ------------- Fanno Line Direct Functions --------------- *)

PROCEDURE AllFanno(CONST mach,gamma : EXTENDED;
   VAR tf,pf,pz,ff,vf,fd : EXTENDED);

  VAR
    gp,gm : EXTENDED;
    msq : EXTENDED;
    tt : EXTENDED;
BEGIN
  gp:=gamma+1.0; gm:=gamma-1.0; msq:=mach*mach;
  tt:=Eq43(mach,gamma);
  tf:=0.5*gp*tt;
  pf:=Sqrt(tf)/mach;
  pz:=Power(1.0/tf, 0.5*gp/gm)/mach;
  ff:=(1.0+gamma*msq)/(mach*(Sqrt(2.0*gp/tt)));
  vf:=msq*pf;
  fd:=(1.0-msq)/(gamma*msq) + 0.5*gp*Ln(tf*msq)/gamma;
END;

FUNCTION FannoP(const mach,gamma : EXTENDED): EXTENDED;             {P/P*}
  VAR tf : EXTENDED;
BEGIN
  tf:=0.5*(gamma+1.0)*Eq43(mach,gamma);
  result:=Sqrt(tf)/mach
END;

FUNCTION FannoT(const mach,gamma : EXTENDED): EXTENDED;             {T/T*}
BEGIN
  result:=0.5*(gamma+1.0)*Eq43(mach,gamma)
END;


FUNCTION FannoRho(const mach,gamma : EXTENDED): EXTENDED;       {rho/rho*}
BEGIN
  result:=-1.0
END;

FUNCTION FannoV(const mach,gamma : EXTENDED): EXTENDED;             {V/V*}
BEGIN
  result:=-1.0
END;

FUNCTION FannoPt(const mach,gamma : EXTENDED): EXTENDED;          {Pt/Pt*}
BEGIN
  result:=-1.0
END;

FUNCTION FannoFLD(const mach,gamma : EXTENDED) : EXTENDED;
  VAR
    gp2,msq,tf : EXTENDED;
BEGIN
  gp2:=0.5*(gamma+1.0);
  msq:=Sqr(mach);
  tf:=gp2*Eq43(mach,gamma);
  result:=(1.0-msq)/(gamma*msq) + gp2*Ln(tf*msq)/gamma;
END;



(* ------------ Fanno Line Inverse Functions ------------- *)


FUNCTION FannoMachFromFLDsubsonic(CONST fld, gamma : EXTENDED): EXTENDED;
  VAR
    gp,gm : EXTENDED;
    mi,mo : EXTENDED;
    msub : EXTENDED;
    p1,p2 : EXTENDED;
BEGIN
//  result:=-1.0;
  gp:=gamma+1.0;
  gm:=gamma-1.0;

  mi:=1;//  mo:=0.5;

  REPEAT
    p1:=0.5*gp*Sqr(mi)/(1.0+0.5*gm*Sqr(mi));
    p2:=gamma*fld+1.0-0.5*gp*Ln(p1);    msub:=Sqrt(1.0/p2);    mo:=mi;
    mi:=msub
  UNTIL Abs(mi-mo) < 1E-10;
  result:=mi
END;   (* ----- End of Function FannoMachFrom FLDsubsonic *)


FUNCTION FannoMachFromFLDsupersonic(CONST fld, gamma : EXTENDED): EXTENDED;
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

  fldmax:=gp*Ln(gp/gm)/(2.0*gamma)-1.0/gamma;  if fld > fldmax THEN
    BEGIN
      ErrorMsg('fL/D exceeds max. value');
      Exit
    END;

  mi:=1;  //  mo:=5;

  REPEAT
    p1:= -2.0*(gamma*fld+1.0-1.0/Sqr(mi))/gp;
    p2:=0.5*gp*Exp(p1)-0.5*gm;   msuper:=Sqrt(1.0/p2);    mo:=mi;
    mi:=msuper
  UNTIL Abs(mi-mo) < 1E-10;
  result:=mi
END;   (* ----- End of Function FannoMachFrom FLDsupersonic *)



FUNCTION FannoMachFromPtPtStarSubsonic(CONST ptpts, gamma : EXTENDED): EXTENDED;
  VAR
    gp,gm : EXTENDED;
    mi,mo : EXTENDED;
    m : EXTENDED;
    p1,p2 : EXTENDED;
BEGIN
  gp:=gamma+1.0;
  gm:=gamma-1.0;

  mi:=0.2; //  mo:=1.0;

  REPEAT
    p1:=(gm*Sqr(mi)+2.0)/gp;
    p2:=Power(p1, 0.5*gp/gm);
    m:=p2/ptpts;    mo:=mi;
    mi:=m
  UNTIL Abs(mi-mo) < 1E-10;
  result:=mi
END;   (* ----- End of Function FannoMachFrom PtPtStarSubsonic *)


FUNCTION FannoMachFromPtPtStarSupersonic(CONST ptpts, gamma : EXTENDED): EXTENDED;
  VAR
    gp,gm : EXTENDED;
    mi,mo : EXTENDED;
    m : EXTENDED;
    p1,p2 : EXTENDED;
BEGIN
  gp:=gamma+1.0;
  gm:=gamma-1.0;

  mi:=5.0;  //  mo:=1.0;

  REPEAT
    p1:=Power(mi*ptpts, 2.0*gm/gp);
    p2:=2.0*(0.5*p1*gp-1.0)/gm;
    m:=Sqrt(Abs(p2));    mo:=mi;
    mi:=m
  UNTIL Abs(mi-mo) < 1E-10;
  result:=mi
END;   (* ----- End of Function FannoMachFrom PtPtStarSupersonic *)


    (* -------------------------------------------- End of Unit RayFanno *)
END.

