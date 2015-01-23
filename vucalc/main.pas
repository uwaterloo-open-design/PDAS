unit Main;
(* AUTHORS - Tom Benson, NASA Glenn Research Center (formerly Lewis)       *)
(*           Ralph L. Carmichael, Public Domain Aeronautical Software      *)
(* REVISION HISTORY                                                        *)
(*   DATE  VERS PERSON  STATEMENT OF CHANGES                               *)
(*    ?     1.0    TB   Original coding (in C for Unix, X-windows, Forms)  *)
(*   1997   2.0   RLC   Recoded in Delphi for Microsoft Windows            *)
(* 21Nov00  6.0   RLC   Gathered all the procedures into a unit            *)

interface

uses
  Windows, Messages, SysUtils, Classes, Graphics, Controls, Forms, Dialogs,
  StdCtrls, ComCtrls, ExtCtrls, Math, NACA1135, CalcSubs, LoadMemo;
type
  TForm1 = class(TForm)

    PageControl1: TPageControl;
    TabSheet1: TTabSheet;
    lblGamma1: TLabel;
    lblGamma2: TLabel;
    lblGamma3: TLabel;
    ebGamma: TEdit;
    btnGammaApply: TButton;
    btnGammaQuit: TButton;
    btnGammaHelp: TButton;





    TabSheet2: TTabSheet;
    rgIsen: TRadioGroup;
    gbIsen: TGroupBox;
    ebIsenInput: TEdit;
    btnIsenCompute: TButton;
    btnIsenQuit: TButton;
    btnIsenHelp: TButton;
    memIsen: TMemo;

    TabSheet3: TTabSheet;
    rgNormal: TRadioGroup;
    gbNormal: TGroupBox;
    ebNormalInput: TEdit;
    btnNormalCompute: TButton;
    btnNormalQuit: TButton;
    btnNormalHelp: TButton;
    memNormal: TMemo;

    TabSheet4: TTabSheet;
    rgOblique: TRadioGroup;
    gbObliqueUpstream: TGroupBox;
    ebObliqueM1: TEdit;
    gbObliqueInput: TGroupBox;
    ebObliqueInput: TEdit;
    btnObliqueCompute: TButton;
    btnObliqueQuit: TButton;
    btnObliqueHelp: TButton;
    memOblique: TMemo;

    TabSheet5: TTabSheet;
    btnAtmosCompute: TButton;
    btnAtmosQuit: TButton;
    btnAtmosHelp: TButton;
    memAtmos: TMemo;
    gbAtmosAlt: TGroupBox;
    ebAtmosAlt: TEdit;
    lbAtmosAlt: TLabel;
    gbAtmosMach: TGroupBox;
    ebAtmosMach: TEdit;
    rgAtmos: TRadioGroup;


    TabSheet6: TTabSheet;
    gbRayleighInput: TGroupBox;
    rgRayleigh: TRadioGroup;
    memRayleigh: TMemo;
    btnRayleighQuit: TButton;
    btnRayleighHelp: TButton;
    ebRayleighInput: TEdit;
    btnRayleighCompute: TButton;

    TabSheet7: TTabSheet;
    rgFanno: TRadioGroup;
    gbFannoInput: TGroupBox;
    ebFannoInput: TEdit;
    btnFannoCompute: TButton;
    btnFannoQuit: TButton;
    btnFannoHelp: TButton;
    memFanno: TMemo;
    TabSheet8: TTabSheet;
    Memo1: TMemo;


    procedure btnGammaApplyClick(Sender: TObject);
    procedure btnGammaQuitClick(Sender: TObject);
    procedure btnGammaHelpClick(Sender: TObject);

    procedure btnIsenComputeClick(Sender: TObject);
    procedure btnIsenQuitClick(Sender: TObject);
    procedure btnIsenHelpClick(Sender: TObject);
    procedure rgIsenClick(Sender: TObject);



    procedure btnRayleighQuitClick(Sender: TObject);
    procedure btnRayleighComputeClick(Sender: TObject);
    procedure btnFannoComputeClick(Sender: TObject);
    procedure btnFannoQuitClick(Sender: TObject);
    procedure btnAtmosQuitClick(Sender: TObject);
    procedure btnObliqueQuitClick(Sender: TObject);
    procedure btnNormalQuitClick(Sender: TObject);

    procedure btnNormalComputeClick(Sender: TObject);
    procedure btnObliqueComputeClick(Sender: TObject);
    procedure btnAtmosComputeClick(Sender: TObject);
    procedure rgAtmosClick(Sender: TObject);
    procedure rgObliqueClick(Sender: TObject);
    procedure rgNormalClick(Sender: TObject);
    procedure rgFannoClick(Sender: TObject);
    procedure rgRayleighClick(Sender: TObject);


    procedure btnNormalHelpClick(Sender: TObject);
    procedure btnObliqueHelpClick(Sender: TObject);
    procedure btnAtmosHelpClick(Sender: TObject);
    procedure btnRayleighHelpClick(Sender: TObject);
    procedure btnFannoHelpClick(Sender: TObject);

  private
    { Private declarations }
  public
    { Public declarations }
  end;

var
  Form1: TForm1;

  gamma : EXTENDED = 1.4;  // used in many procedures; initial value is 1.4

implementation

USES GetMachFromInput;

{$R *.DFM}

PROCEDURE ErrorMsg(CONST s : STRING);
BEGIN
  MessageDlg(s, mtError, [mbOK], 0)
END;   (* ------------------------ End of Procedure ErrorMsg *)


(* ---------- Gamma Page Procedures ------------------- *)


procedure TForm1.btnGammaApplyClick(Sender: TObject);
  VAR
    errCode : INTEGER;
    xx : EXTENDED;

begin
  Val(ebGamma.Text, xx, errCode);
  IF errCode <> 0 THEN
    BEGIN
      ErrorMsg('invalid input string for gamma');
      Exit
    END;

  IF xx < 1.0 THEN
      ErrorMsg('gamma must be >= 1.0')
  ELSE
    IF xx > 1.67 THEN
      ErrorMsg('gamma must not exceed 5/3')
    ELSE
      gamma:=xx

end;


procedure TForm1.btnGammaQuitClick(Sender: TObject);
begin
  Application.Terminate
end;

procedure TForm1.btnGammaHelpClick(Sender: TObject);
begin
  ShowMessage('To change gamma, type the new value in the edit box')
end;

(* ---------- Isentropic Flow Procedures ----------------- *)

procedure TForm1.btnIsenComputeClick(Sender: TObject);
  VAR
    errCode : INTEGER;
    inputValue : EXTENDED;
    mach : EXTENDED;
begin
  Val(ebIsenInput.Text,inputValue,errCode);
  IF errCode <> 0 THEN
    BEGIN
      ErrorMsg('Invalid input value in isentropic flow');
      Exit
    END;

  mach:=GetMachFromIsenInput(rgIsen.ItemIndex,
     inputValue, gamma);
  LoadIsenMemo(memIsen,mach,gamma);
  ebIsenInput.SetFocus
end;

procedure TForm1.btnIsenQuitClick(Sender: TObject);
begin
  Application.Terminate
end;

procedure TForm1.btnIsenHelpClick(Sender: TObject);
begin
  ShowMessage('First select the flow parameter you want to '+
    'specify in the Input quantity box. '+
    'Then enter the value of this parameter in the ' +
    'Input value box. Click the Compute button to display the ' +
    'results. Click the Quit button to exit the program.');
end;

procedure TForm1.rgIsenClick(Sender: TObject);
begin
  ebIsenInput.Clear;
  ebIsenInput.SetFocus
end;

(* -------------------- Normal Shock Procedures ------------ *)

procedure TForm1.btnNormalComputeClick(Sender: TObject);
  VAR
    errCode : INTEGER;
    inputValue : EXTENDED;
    mach : EXTENDED;
begin
  Val(ebNormalInput.Text,inputValue,errCode);
  IF errCode <> 0 THEN
    BEGIN
      ErrorMsg('Invalid input value in normal shock');
      Exit
    END;

  mach:=GetMachFromNormalInput(rgNormal.ItemIndex,
     inputValue, gamma);
  LoadNormalMemo(memNormal,mach,gamma);
  ebNormalInput.SetFocus
end;

procedure TForm1.btnNormalQuitClick(Sender: TObject);
begin
  Application.Terminate
end;

procedure TForm1.btnNormalHelpClick(Sender: TObject);
begin
  ShowMessage('First select the flow parameter you want to '+
    'specify in the Input quantity box. '+
    'Then enter the value of this parameter in the ' +
    'Input value box. Click the Compute button to display the ' +
    'results. Click the Quit button to exit the program.');
end;

procedure TForm1.rgNormalClick(Sender: TObject);
begin
  ebNormalInput.Clear;
  ebNormalInput.SetFocus
end;

(* ----------- Oblique Shock Procedures -------------------- *)

procedure TForm1.btnObliqueComputeClick(Sender: TObject);
  VAR
    errCode : INTEGER;
    inputValue : EXTENDED;
    mach1 : EXTENDED;
    theta : EXTENDED;
begin
  Val(ebObliqueM1.Text, mach1,errCode);
  IF errCode <> 0 THEN
    BEGIN
      ErrorMsg('Invalid value for upstream Mach in oblique shock');
      Exit
    END;

  Val(ebObliqueInput.Text,inputValue,errCode);
  IF errCode <> 0 THEN
    BEGIN
      ErrorMsg('Invalid input value in oblique shock');
      Exit
    END;

  theta:=GetThetaFromObliqueInput(rgOblique.ItemIndex,
     inputValue, mach1, gamma);
  LoadObliqueMemo(memOblique,mach1,theta,gamma);
  ebObliqueInput.SetFocus
end;

procedure TForm1.btnObliqueQuitClick(Sender: TObject);
begin
  Application.Terminate
end;

procedure TForm1.btnObliqueHelpClick(Sender: TObject);
begin
  ShowMessage('First select the flow parameter you want to '+
    'specify in the Input quantity box. '+
    'Then enter the value of this parameter in the ' +
    'Input value box. '+
    'Set the value of upstream Mach number.'+
    'Click the Compute button to display the ' +
    'results. Click the Quit button to exit the program.');

end;

procedure TForm1.rgObliqueClick(Sender: TObject);
begin
  ebObliqueInput.Clear;
  ebObliqueInput.SetFocus
end;


(* ----------- Standard Atmosphere Procedures ------------- *)

procedure TForm1.btnAtmosComputeClick(Sender: TObject);
  CONST
    FT2KM = 0.0003048;

  VAR
    alt : EXTENDED;
    altKm : EXTENDED;
    errCode : INTEGER;
    mach : EXTENDED;
begin
  Val(ebAtmosMach.Text, mach, errCode);
  IF errCode <> 0 THEN
    BEGIN
      ErrorMsg('Invalid mach in Atmosphere');
      Exit
    END;
  IF mach < 0.0 THEN
    BEGIN
      ErrorMsg('Mach cannot be < 0');
      Exit
    END;

  Val(ebAtmosAlt.Text,alt,errCode);
  IF errCode <> 0 THEN
    BEGIN
      ErrorMsg('Invalid altitude in Atmosphere');
      Exit
    END;

  IF rgAtmos.ItemIndex = 0 THEN
    BEGIN  {US units}
      altKm:=FT2KM*alt;
      lbAtmosAlt.Caption:='ft.';
      LoadAtmosMemoUSunits(memAtmos,altKm,mach)
    END
  ELSE
    BEGIN   {SI units}
      altKm:=0.001*alt;
      lbAtmosAlt.Caption:='m';
      LoadAtmosMemoSIunits(memAtmos,altKm,mach)
    END;
  ebAtmosAlt.SetFocus
end;

procedure TForm1.rgAtmosClick(Sender: TObject);
  CONST
    FT2METERS = 0.3048;

  VAR
    altFt,altMeters : SINGLE;
    errCode : INTEGER;

begin
 IF rgAtmos.itemIndex = 0 THEN
    IF lbAtmosAlt.Caption = 'm' THEN
      BEGIN
        Val(ebAtmosAlt.Text, altMeters, errCode);
        altFt:=altMeters/FT2Meters;
        ebAtmosAlt.Text:=FloatToStrF(altFt, ffFixed,8,0);
        btnAtmosComputeClick(rgAtmos);
      END;

  IF rgAtmos.itemIndex = 1 THEN
    IF lbAtmosAlt.Caption = 'ft.' THEN
      BEGIN
        Val(ebAtmosAlt.Text, altFt, errCode);
        altMeters:=altFt*FT2METERS;
        ebAtmosAlt.Text:=FloatToStrF(altMeters, ffFixed,8,0);
        btnAtmosComputeClick(rgAtmos);
      END;

end;

procedure TForm1.btnAtmosQuitClick(Sender: TObject);
begin
  Application.Terminate
end;

procedure TForm1.btnAtmosHelpClick(Sender: TObject);
begin
  ShowMessage('First enter appropriate values for altitude and Mach number ' +
    'in the marked boxes. ' +
    'You may select US or SI units by clicking a radio button. ' +
    'Click the Compute button to display the results. ' +
    'If you toggle US/SI, the units change, but the flight condition is ' +
    'the same. ' +
    'Click the Quit button to exit the program.');

end;


(* ----------------- Rayleigh Line Procedures -------------- *)

procedure TForm1.btnRayleighComputeClick(Sender: TObject);
  VAR
    mach : EXTENDED;
    errCode : INTEGER;
    rayleighInput : EXTENDED;

begin
  Val(ebRayleighInput.Text, rayleighInput, errCode);
  IF errCode < 0 THEN
    BEGIN
      ErrorMsg('Invalid input for Rayleigh');
      Exit
    END;
  mach:=GetMachFromRayleighInput(gamma, rayleighInput, rgRayleigh.ItemIndex);
  LoadRayleighMemo(memRayleigh, mach,gamma);
  ebRayleighInput.SetFocus
end;

procedure TForm1.btnRayleighQuitClick(Sender: TObject);
begin
  Application.Terminate
end;

procedure TForm1.btnRayleighHelpClick(Sender: TObject);
begin
  ShowMessage('First select the flow parameter you want to '+
    'specify in the Input quantity box. '+
    'Then enter the value of this parameter in the ' +
    'Input value box. Click the Compute button to display the ' +
    'results. Click the Quit button to exit the program.');

end;

procedure TForm1.rgRayleighClick(Sender: TObject);
begin
  ebRayleighInput.Clear;
  ebRayleighInput.SetFocus
end;


(* --------------------- Fanno Line Procedures ----------- *)

procedure TForm1.btnFannoComputeClick(Sender: TObject);
  VAR
    mach : EXTENDED;
    errCode : INTEGER;
    fannoInput : EXTENDED;

BEGIN
  Val(ebFannoInput.Text, fannoInput, errCode);
  IF errCode  <> 0 THEN
    BEGIN
      ErrorMsg('Invalid input for Fanno');
      Exit
    END;


  mach:=GetMachFromFannoInput(gamma, fannoInput, rgFanno.ItemIndex);
  LoadFannoMemo(memFanno, mach,gamma);
  ebFannoInput.SetFocus
end;

procedure TForm1.btnFannoQuitClick(Sender: TObject);
begin
  Application.Terminate
end;

procedure TForm1.btnFannoHelpClick(Sender: TObject);
begin
  ShowMessage('First select the flow parameter you want to '+
    'specify in the Input quantity box. '+
    'Then enter the value of this parameter in the ' +
    'Input value box. Click the Compute button to display the ' +
    'results. Click the Quit button to exit the program.');

end;

procedure TForm1.rgFannoClick(Sender: TObject);
begin
  ebFannoInput.Clear;
  ebFannoInput.SetFocus
end;


end.
