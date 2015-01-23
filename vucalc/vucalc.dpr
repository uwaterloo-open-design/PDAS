program VuCalc;

uses
  Forms,
  Main in 'Main.pas' {Form1},
  NACA1135 in 'Naca1135.pas',
  CalcSubs in 'Calcsubs.pas',
  LoadMemo in 'loadmemo.pas',
  getMachFromInput in 'GetMachFromInput.pas',
  RayFanno in 'RayFanno.pas';

{$R *.RES}

begin
  Application.Initialize;
  Application.CreateForm(TForm1, Form1);
  Application.Run;
end.
