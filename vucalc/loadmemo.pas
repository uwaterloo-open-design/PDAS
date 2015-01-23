UNIT LoadMemo;
(* AUTHORS - Tom Benson, NASA Glenn Research Center (formerly Lewis)       *)
(*           Ralph L. Carmichael, Public Domain Aeronautical Software      *)
(* REVISION HISTORY                                                        *)
(*   DATE  VERS PERSON  STATEMENT OF CHANGES                               *)
(*    ?     1.0    TB   Original coding (in C for Unix, X-windows, Forms)  *)
(*   1997   2.0   RLC   Recoded in Delphi for Microsoft Windows            *)
(* 21Nov00  6.0   RLC   Gathered all the procedures into a unit            *)

INTERFACE
USES StdCtrls, Sysutils, Math, NACA1135, RayFanno;

PROCEDURE LoadIsenMemo( CONST memo : TMemo;
                        CONST mach,gamma : EXTENDED);
PROCEDURE LoadNormalMemo(CONST memo : TMemo;
       const mach1,gamma : EXTENDED);
PROCEDURE LoadObliqueMemo(CONST memo : TMemo; const mach1,theta,gamma : EXTENDED);

PROCEDURE LoadAtmosMemoSIunits(CONST memo : TMemo;
              const altKm,mach : EXTENDED);
PROCEDURE LoadAtmosMemoUSunits(CONST memo : TMemo;
              const altKm,mach : EXTENDED);

PROCEDURE LoadRayleighMemo(CONST memo : TMemo;
                           CONST mach,gamma : EXTENDED);

PROCEDURE LoadFannoMemo(CONST memo : TMemo;
                        CONST mach,gamma : EXTENDED);


IMPLEMENTATION

uses GetMachFromInput;

CONST
    TZERO      = 288.15;                 { temperature at sealevel, kelvins }
    PZERO      = 101325.0;              { pressure at sealevel, N/sq.m. = Pa}
    RHOZERO    = 1.2250;                    { density at sealevel, kg/cu.m. }
    ASOUNDZERO = 340.294;               { speed of sound at sealevel, m/sec }

    FT2METERS = 0.3048;                   { mult. ft. to get meters (exact) }
    KELVIN2RANKINE = 1.8;                         { mult deg K to get deg R }
    PSF2NSM = 47.880258;                      { mult lb/sq.ft to get N/sq.m }
    SCF2KCM = 515.379;                    { mult slugs/cu.ft to get kg/cu.m }
    FT2KM = 0.0003048;


FUNCTION Formatted(CONST x : EXTENDED) : STRING;

BEGIN
  IF (x >= 0.001) AND (x < 1000.0) THEN
    result:=FloatToStrF(x,ffFixed,20,4)
  ELSE
    result:=FloatToStrF(x,ffExponent,5,4)
END;


PROCEDURE Atmosphere(alt : EXTENDED;    { geometric altitude, in km }
               VAR sigma : EXTENDED;    { density/sea-level density }
               VAR delta : EXTENDED;    { pressure/sea-level pressure }
               VAR theta : EXTENDED);   { temperature/sea-level temperature }

CONST
  TABLESIZE = 8;
  REARTH = 6369.0;         { radius of the Earth (km)}
  GMR = 34.163195;         { gas constant }

TYPE
  TABLEINDEX = 1..TABLESIZE;
  ATMOSTABLE = ARRAY[TABLEINDEX] OF EXTENDED;

CONST
  htab : ATMOSTABLE = (0.0, 11.0, 20.0, 32.0, 47.0, 51.0, 71.0, 84.852);
  ttab : ATMOSTABLE = (288.15, 216.65, 216.65, 228.65, 270.65,
                       270.65, 214.65, 186.87 );
  ptab : ATMOSTABLE = (1.0, 2.2336110E-1, 5.4032950E-2, 8.5666784E-3,
                     1.0945601E-3, 6.6063531E-4, 3.9046834E-5, 3.685010E-6 );
  gtab : ATMOSTABLE = (-6.5, 0.0, 1.0, 2.8, 0.0, -2.8, -2.0, 0.0);

VAR
   i,j,k : WORD;
   h, tgrad, deltah, tbase, tlocal : EXTENDED;
  BEGIN
  h:=alt*REARTH/(alt+REARTH);     { Convert geometric to geopotential altitude }

  i:=1;
  j:=TABLESIZE;
  REPEAT                                   { binary search in ordered table }
    k:=(i+j) DIV 2;
    IF h < htab[k] THEN j:=k ELSE i:=k
  UNTIL j <= i+1;

  tgrad:=gtab[i];
  deltah:=h-htab[i];
  tbase:=ttab[i];

  tlocal:=tbase+tgrad*deltah;
  theta:=tlocal/ttab[1];                                { temperature ratio }

  IF tgrad=0.0 THEN                                        { pressure ratio }
    delta:=ptab[i]*Exp(-GMR*deltah/tbase)
  ELSE
    delta:=ptab[i]*Exp(Ln(tbase/tlocal)*GMR/tgrad);

  sigma:=delta/theta;                                       { density ratio }
END;   (* ----------------------------------- End of Procedure Atmosphere. *)

(* ======================================================================= *)
FUNCTION MetricViscosity(const theta : EXTENDED) : EXTENDED;
  CONST
    TZERO      = 288.15;     { temperature at sealevel, kelvins }
    BETAVISC   = 1.458E-6;   { viscosity term, N sec/(sq.m Sqrt(K) }
    SUTH       = 110.4;      { Sutherland's constant, kelvins }
  VAR
    t : EXTENDED;                             { temperature, Kelvin }
BEGIN                                   { returns viscosity in metric units }
  t:=TZERO*theta;
  MetricViscosity:=BETAVISC*Sqrt(t*t*t)/(t+SUTH)     { Sutherland's formula }
END;   (* -------------------------------- End of Function MetricViscosity *)


PROCEDURE LoadIsenMemo( CONST memo : TMemo;
                        CONST mach,gamma : EXTENDED);
  VAR
    buffer : STRING;
    xx : EXTENDED;
BEGIN
  memo.Clear;
  buffer:='Mach Number = ' + FloatToStrF(mach, ffFixed, 18, 4);
  memo.Lines.Add(buffer);

  buffer:='gamma=' + FloatToStrF(gamma, ffFixed,16,3);
  memo.Lines.Add(buffer);

  xx:=Eq47(mach,gamma);                                     {q/p}
  buffer:='q/p=' + FloatToStrF(xx,ffFixed,18,4);
  memo.Lines.Add(buffer);

  xx:=Eq44(mach,gamma);                                   { p/pt }
  buffer:='p/pt=' + FloatToStrF(xx,ffFixed,18,5);
  memo.Lines.Add(buffer);

  xx:=Eq48(mach,gamma);                                    { q/pt }
  buffer:='q/pt=' + FloatToStrF(xx,ffFixed,18,5);
  memo.Lines.Add(buffer);

  xx:=Eq80(mach,gamma);       { returns A*/A }
  IF xx > 0.0 THEN xx:=1.0/xx;
  buffer:='A/A*=' + FloatToStrF(xx,ffFixed,18,5);
  memo.Lines.Add(buffer);

  xx:=Eq45(mach,gamma);
  buffer:='rho/rho-t=' + FloatToStrF(xx,ffFixed,18,5);
  memo.Lines.Add(buffer);

  xx:=Eq43(mach,gamma);
  buffer:='T/Tt=' + FloatToStrF(xx,ffFixed,18,4);
  memo.Lines.Add(buffer);

  xx:=Sqrt(Eq50(mach,gamma)); { N.B. Eq50 returns SQUARE of ratio }
  buffer:='V/a*=' + FloatToStrF(xx,ffFixed,18,4);
  memo.Lines.Add(buffer);

  IF mach < 1.0 THEN Exit;

  xx:=RadToDeg(Eq171c(mach,gamma));
  buffer:='nu=' + FloatToStrF(xx,ffFixed,18,3) + ' deg.';
  memo.Lines.Add(buffer);

  xx:=RadToDeg(ArcSin(1.0/mach));
  buffer:='mu=' + FloatToStrF(xx,ffFixed,18,3) + ' deg.';
  memo.Lines.Add(buffer);

END;   (* ---------------------------------- End of Procedure LoadIsenMemo *)

PROCEDURE LoadNormalMemo(CONST memo : TMemo;
                         const mach1,gamma : EXTENDED);
  VAR
    buffer : STRING;
    xx : EXTENDED;
BEGIN
  memo.Clear;
  buffer:='Upstream Mach Number= ' + FloatToStrF(mach1, ffFixed, 18, 4);
  memo.Lines.Add(buffer);

  buffer:='gamma=' + FloatToStrF(gamma, ffFixed,16,3);
  memo.Lines.Add(buffer);

  xx:=Sqrt(Eq96(mach1,gamma));                                     { M2 }
  buffer:='Downstream Mach Number= ' + FloatToStrF(xx,ffFixed,18,4);
  memo.Lines.Add(buffer);

  xx:=Eq99(mach1,gamma);                                     { pt2/pt1 }
  buffer:='Total Pressure Ratio= ' + FloatToStrF(xx,ffFixed,18,4);
  memo.Lines.Add(buffer);

  xx:=Eq93(mach1,gamma);                                   {p2/p1}
  buffer:='Static Pressure Ratio= ' + FloatToStrF(xx,ffFixed,18,4);
  memo.Lines.Add(buffer);

  xx:=Eq95(mach1,gamma);                                   { T2/T1}
  buffer:='Static Temperature Ratio= ' + FloatToStrF(xx,ffFixed,18,4);
  memo.Lines.Add(buffer);

  xx:=Eq94(mach1,gamma);                              {rho2/rho1}
  buffer:='Density ratio= ' + FloatToStrF(xx,ffFixed,18,4);
  memo.Lines.Add(buffer);

END;   (* -------------------------------- End of Procedure LoadNormalMemo *)

PROCEDURE LoadObliqueMemo(CONST memo : TMemo;
                          const mach1,theta,gamma : EXTENDED);
  VAR
    buffer : STRING;
    xx : EXTENDED;
BEGIN
  memo.Clear;

  buffer:='Upstream Mach Number = ' + FloatToStrF(mach1, ffFixed, 18, 4);
  memo.Lines.Add(buffer);

  buffer:='gamma=' + FloatToStrF(gamma, ffFixed,16,3);
  memo.Lines.Add(buffer);

  xx:=Eq132(mach1,theta,gamma);                                    { M2 }
  buffer:='Downstream Mach Number=' + FloatToStrF(xx,ffFixed,18,4);
  memo.Lines.Add(buffer);

  xx:=Eq138(mach1,theta,gamma);                            { ramp angle }
  xx:=RadToDeg(ArcTan(1.0/xx));
  buffer:='Wedge angle=' +
    FloatToStrF(xx,ffFixed,18,3) + ' deg.';
  memo.Lines.Add(buffer);


  buffer:='Wave angle=' +
  FloatToStrF(RadToDeg(theta),ffFixed,18,3) + ' deg.';
  memo.Lines.Add(buffer);

  xx:=Eq142(mach1,theta,gamma);                               { pt2/pt1 }
  buffer:='Total Pressure Ratio=' + FloatToStrF(xx,ffFixed,18,4);
  memo.Lines.Add(buffer);

  xx:=Eq128(mach1,theta,gamma);                                  {p2/p1 }
  buffer:='Static Pressure Ratio=' + FloatToStrF(xx,ffFixed,18,4);
  memo.Lines.Add(buffer);

  xx:=Eq130(mach1,theta,gamma);                                   {T2/T1}
  buffer:='Static Temperature Ratio=' + FloatToStrF(xx,ffFixed,18,4);
  memo.Lines.Add(buffer);

  xx:=Eq129(mach1,theta,gamma);                             { rho2/rho1 }
  buffer:='Density ratio=' + FloatToStrF(xx,ffFixed,18,4);
  memo.Lines.Add(buffer);

  buffer:='All ratios are downstream/upstream';
  memo.Lines.Add(buffer)

END;   (* ------------------------------- End of Procedure LoadObliqueMemo *)

PROCEDURE LoadAtmosMemoSIunits(CONST memo : TMemo;
                               const altKm,mach : EXTENDED);
  CONST
    GAMMA = 1.4;

  VAR
    buffer : STRING;   { assemble the output string here }
    sigma,delta,theta : EXTENDED;
    xx : EXTENDED;               { holds current value being computed }
BEGIN
  memo.Clear;
  Atmosphere(altKm, sigma,delta,theta);    { sets sigma,delta,theta }

  buffer:='density/sea-level density= ' + FloatToStrF(sigma, ffExponent, 5,4);
  memo.Lines.Add(buffer);
  buffer:='pressure/sea-level pressure= ' + FloatToStrF(delta, ffExponent, 5,4);
  memo.Lines.Add(buffer);
  buffer:='temperature/sea-level temperature= ' +
                                              FloatToStrF(theta, ffFixed, 18,4);
  memo.Lines.Add(buffer);

  xx:=RHOZERO*sigma;   { density }
  buffer:='density= ' + FloatToStrF(xx, ffExponent, 5,4) + ' kg/cu.m.';
  memo.Lines.Add(buffer);

  xx:=PZERO*delta;   { pressure }
  buffer:='static pressure= ' + FloatToStrF(xx, ffFixed, 18,4) + ' Pa';
  memo.Lines.Add(buffer);

  xx:=TZERO*theta;   { temperature }
  buffer:='static temperature= ' + FloatToStrF(xx, ffFixed, 18,4) + ' K';
  memo.Lines.Add(buffer);

  xx:=ASOUNDZERO*Sqrt(theta);   { speed of sound }
  buffer:='speed of sound= ' + FloatToStrF(xx, ffFixed, 18,1) + ' m/s';
  memo.Lines.Add(buffer);

  xx:=1E6*MetricViscosity(theta);   { viscosity }
  buffer:='viscosity= ' + FloatToStrF(xx, ffFixed, 18,5) + 'E-6 kg/m-s';
  memo.Lines.Add(buffer);

  xx:=MetricViscosity(theta)/(RHOZERO*sigma);   { kinematic viscosity }
  buffer:='kinematic viscosity= ' + FloatToStrF(xx, ffExponent, 8,4) + ' sq.m/s';
  memo.Lines.Add(buffer);

  xx:=ASOUNDZERO*Sqrt(theta)*mach*RHOZERO*sigma/MetricViscosity(theta);   { unit Reynolds number }
  xx:=xx*1E-6;
  buffer:='unit Reynolds number= ' + FloatToStrF(xx, ffFixed, 18,4) + ' million 1/m';
  memo.Lines.Add(buffer);

  xx:=ASOUNDZERO*Sqrt(theta)*mach;   { velocity }
  buffer:='velocity= ' + FloatToStrF(xx, ffFixed, 18,4) + ' m/s';
  memo.Lines.Add(buffer);

  xx:=0.5*GAMMA*PZERO*delta*mach*mach;   { dynamic pressure }
  buffer:='dynamic pressure= ' + FloatToStrF(xx, ffFixed, 18,4) + ' Pa';
  memo.Lines.Add(buffer);

  xx:=TZERO*theta/Eq43(mach,gamma);   { total temperature }
  buffer:='total temperature= ' + FloatToStrF(xx, ffFixed, 18,4) + ' K';
  memo.Lines.Add(buffer);

  xx:=PZERO*delta/Eq44(mach,gamma);   { total pressure }
  buffer:='total pressure= ' + FloatToStrF(xx, ffFixed, 18,4) + ' Pa';
  memo.Lines.Add(buffer);

END;   (* -------------------------------- End of Procedure LoadMemoSIunits *)

PROCEDURE LoadAtmosMemoUSunits(CONST memo : TMemo;
                               const altKm,mach : EXTENDED);
  CONST
    GAMMA = 1.4;

  VAR
    buffer : STRING;
    sigma,delta,theta : EXTENDED;
    xx : EXTENDED;
BEGIN
  memo.Clear;
  Atmosphere(altKm, sigma,delta,theta);

  buffer:='density/sea-level density= ' + FloatToStrF(sigma, ffExponent, 5,4);
  memo.Lines.Add(buffer);
  buffer:='pressure/sea-level pressure= ' + FloatToStrF(delta, ffExponent, 5,4);
  memo.Lines.Add(buffer);
  buffer:='temperature/sea-level temperature= ' + FloatToStrF(theta, ffFixed, 18,4);
  memo.Lines.Add(buffer);

  xx:=RHOZERO*sigma/SCF2KCM;   { density }
  buffer:='density= ' + FloatToStrF(xx, ffExponent, 5,4) + ' slugs/cu.ft.';
  memo.Lines.Add(buffer);

  xx:=PZERO*delta/PSF2NSM;   { pressure }
  buffer:='static pressure= ' + FloatToStrF(xx, ffFixed, 18,4) + ' lbs./sq.ft.';
  memo.Lines.Add(buffer);

  xx:=TZERO*theta*KELVIN2RANKINE;   { temperature }
  buffer:='static temperature= ' + FloatToStrF(xx, ffFixed, 18,1) + ' deg. Rankine';
  memo.Lines.Add(buffer);

  xx:=ASOUNDZERO*Sqrt(theta)/FT2METERS;   { speed of sound }
  buffer:='speed of sound= ' + FloatToStrF(xx, ffFixed, 18,1) + ' ft/sec';
  memo.Lines.Add(buffer);

  xx:=1E6*MetricViscosity(theta)/PSF2NSM;   { viscosity }
  buffer:='viscosity= ' + FloatToStrF(xx, ffFixed, 18,5) + ' E-6 slugs/ft.-sec';
  memo.Lines.Add(buffer);

  xx:=MetricViscosity(theta)/(RHOZERO*sigma);   { kinematic viscosity}
  xx:=xx*SCF2KCM/PSF2NSM;
  buffer:='kinematic viscosity= ' + FloatToStrF(xx, ffExponent, 8,4) + ' sq.ft/sec';
  memo.Lines.Add(buffer);

  xx:=ASOUNDZERO*Sqrt(theta)*mach*RHOZERO*sigma/MetricViscosity(theta);   { unit Reynolds number }
  xx:=xx*1E-6*FT2METERS;
  buffer:='unit Reynolds number= ' + FloatToStrF(xx, ffFixed, 18,4) + ' million 1/ft';
  memo.Lines.Add(buffer);

  xx:=ASOUNDZERO*Sqrt(theta)*mach/FT2METERS;   { velocity }
  buffer:='velocity= ' + FloatToStrF(xx, ffFixed, 18,4) + ' ft/sec';
  memo.Lines.Add(buffer);

  xx:=0.5*GAMMA*PZERO*delta*mach*mach/PSF2NSM;   { dynamic pressure }
  buffer:='dynamic pressure= ' + FloatToStrF(xx, ffFixed, 18,4) + ' lbs/sq.ft.';
  memo.Lines.Add(buffer);

  xx:=TZERO*theta*KELVIN2RANKINE/Eq43(mach,gamma);   { total temperature }
  buffer:='total temperature= ' + FloatToStrF(xx, ffFixed, 18,4) + ' deg Rankine';
  memo.Lines.Add(buffer);

  xx:=PZERO*delta/(PSF2NSM*Eq44(mach,gamma));   { total pressure }
  buffer:='total pressure= ' + FloatToStrF(xx, ffFixed, 18,4) + ' lbs/sq.ft.';
  memo.Lines.Add(buffer);

END;   (* -------------------------------- End of Procedure LoadAtmosMemoUSunits *)



PROCEDURE LoadRayleighMemo(CONST memo : TMemo;
                           CONST mach,gamma : EXTENDED);

  VAR
    buffer : STRING;
//    gm,gp : EXTENDED;
//    msq : EXTENDED;
    pr : EXTENDED;
    pq : EXTENDED;
    rr : EXTENDED;
    tq : EXTENDED;
    tr : EXTENDED;
//    tt : EXTENDED;
    vr : EXTENDED;
BEGIN
  IF mach = 0.0 THEN
    BEGIN
      memo.Clear;
      memo.Lines.Add('Mach number=0');
      memo.Lines.Add('p/p*=gamma+1 (=2.4)');
      memo.Lines.Add('p0/p0*=1.2679');
      memo.Lines.Add('V/V*=0');
      memo.Lines.Add('rho/rho*=infinity');
      Exit
    END;

//  gp:=gamma+1.0; gm:=gamma-1.0; msq:=mach*mach;
//  tt:=Eq43(mach,gamma);
//  tf:=0.5*gp*tt;
//  pr:=gp/(1.0+gamma*msq);
//  tr:=pr*pr*msq;
//  tq:=2.0*tr/(gp*tt);
//  vr:=pr*msq;
//  pq:=pr*Power(2.0/(gp*tt), gamma/gm);
//  rr:=RayleighRho(mach,gamma);

  AllRayleigh(mach,gamma, tr,pr,tq,vr,pq,rr);
  memo.Clear;

  buffer:='Mach number=' + FloatToStrF(mach, ffFixed,16,5);
  memo.Lines.Add(buffer);

  buffer:='gamma=' + FloatToStrF(gamma, ffFixed,16,3);
  memo.Lines.Add(buffer);

  buffer:='T0/T0*='+FloatToStrF(tq,ffFixed,16,6);
  memo.Lines.Add(buffer);

  buffer:='T/T*='+FloatToStrF(tr,ffFixed,16,6);
  memo.Lines.Add(buffer);

  buffer:='p/p*='+FloatToStrF(pr,ffFixed,16,6);
  memo.Lines.Add(buffer);

  buffer:='p0/p0*='+FloatToStrF(pq,ffFixed,16,6);
  memo.Lines.Add(buffer);

  buffer:='V/V*='+FloatToStrF(vr,ffFixed,16,6);
  memo.Lines.Add(buffer);

  buffer:='rho / rho* = '+FloatToStrF(rr,ffFixed,16,6);
  memo.Lines.Add(buffer);

END;   (* --------------------- End of Procedure LoadRayleighMemo *)


PROCEDURE LoadFannoMemo(CONST memo : TMemo;
                        CONST mach,gamma : EXTENDED);

  VAR
    buffer : STRING;
    fd : EXTENDED;
    ff : EXTENDED;
//    gm,gp : EXTENDED;
//    msq : EXTENDED;
    pf : EXTENDED;
    pz : EXTENDED;
    tf : EXTENDED;
//    tt : EXTENDED;
    vf : EXTENDED;
BEGIN
  IF mach = 0.0 THEN
    BEGIN
      memo.Clear;
      memo.Lines.Add('Mach number=0');
      memo.Lines.Add('T/T*=1.2');
      memo.Lines.Add('p/p*=infinity');
      memo.Lines.Add('p0/p0*=infinity');
      memo.Lines.Add('V/V*=0');
      memo.Lines.Add('F/F*=infinity');
      memo.Lines.Add('4fL/D=infinity');
      Exit
    END;

//  gp:=gamma+1.0; gm:=gamma-1.0; msq:=mach*mach;
//  tt:=Eq43(mach,gamma);
//  tf:=0.5*gp*tt;
//  pf:=Sqrt(tf)/mach;
//  pz:=Power(1.0/tf, 0.5*gp/gm)/mach;
//  ff:=(1.0+gamma*msq)/(mach*(Sqrt(2.0*gp/tt)));
//  vf:=msq*pf;
//  fd:=(1.0-msq)/(gamma*msq) + 0.5*gp*Ln(tf*msq)/gamma;

  AllFanno(mach,gamma,tf,pf,pz,ff,vf,fd);
  memo.Clear;

  buffer:='Mach number=' + FloatToStrF(mach, ffFixed,16,5);
  memo.Lines.Add(buffer);

  buffer:='gamma=' + FloatToStrF(gamma, ffFixed,16,3);
  memo.Lines.Add(buffer);

  buffer:='T/T*='+FloatToStrF(tf,ffFixed,16,6);
  memo.Lines.Add(buffer);

  buffer:='p/p*='+FloatToStrF(pf,ffFixed,16,6);
  memo.Lines.Add(buffer);

  buffer:='p0/p0*='+FloatToStrF(pz,ffFixed,16,6);
  memo.Lines.Add(buffer);

  buffer:='V/V*='+FloatToStrF(vf,ffFixed,16,6);
  memo.Lines.Add(buffer);

  buffer:='F/F*='+FloatToStrF(ff,ffFixed,16,6);
  memo.Lines.Add(buffer);

  buffer:='4fL/D='+FloatToStrF(fd,ffFixed,16,6);
  memo.Lines.Add(buffer);

END;   (* --------------------- End of Procedure LoadFannoMemo *)

END.