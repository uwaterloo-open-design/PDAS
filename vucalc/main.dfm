object Form1: TForm1
  Left = 195
  Top = 110
  Caption = '-'
  ClientHeight = 446
  ClientWidth = 632
  Color = clBtnFace
  Font.Charset = DEFAULT_CHARSET
  Font.Color = clWindowText
  Font.Height = -12
  Font.Name = 'Arial'
  Font.Style = []
  OldCreateOrder = False
  Visible = True
  PixelsPerInch = 96
  TextHeight = 15
  object PageControl1: TPageControl
    Left = 0
    Top = 0
    Width = 632
    Height = 446
    ActivePage = TabSheet8
    Align = alClient
    Font.Charset = DEFAULT_CHARSET
    Font.Color = clWindowText
    Font.Height = -13
    Font.Name = 'Arial'
    Font.Style = []
    ParentFont = False
    ParentShowHint = False
    ShowHint = False
    TabOrder = 0
    object TabSheet1: TTabSheet
      Caption = 'gamma'
      object lblGamma1: TLabel
        Left = 28
        Top = 35
        Width = 370
        Height = 16
        Caption = 'The value of gamma, the ratio of specific heats may be set here.'
      end
      object lblGamma2: TLabel
        Left = 28
        Top = 58
        Width = 425
        Height = 16
        Caption = 
          'This value wil be used in all calculations except the standard a' +
          'tmosphere.'
      end
      object lblGamma3: TLabel
        Left = 28
        Top = 80
        Width = 387
        Height = 16
        Caption = 
          'After entering a value  1 <= gamma <= 5/3,  click the Apply butt' +
          'on.'
      end
      object btnGammaHelp: TButton
        Left = 230
        Top = 144
        Width = 64
        Height = 22
        Caption = 'Help'
        TabOrder = 0
        OnClick = btnGammaHelpClick
      end
      object btnGammaQuit: TButton
        Left = 143
        Top = 144
        Width = 65
        Height = 22
        Hint = 'quit VuCalc'
        Cancel = True
        Caption = 'Quit'
        ParentShowHint = False
        ShowHint = True
        TabOrder = 1
        OnClick = btnGammaQuitClick
      end
      object ebGamma: TEdit
        Left = 141
        Top = 106
        Width = 107
        Height = 24
        TabOrder = 2
        Text = '1.4'
      end
      object btnGammaApply: TButton
        Left = 56
        Top = 144
        Width = 65
        Height = 22
        Caption = 'Apply'
        Default = True
        TabOrder = 3
        OnClick = btnGammaApplyClick
      end
    end
    object TabSheet2: TTabSheet
      Caption = 'isentropic flow'
      Font.Charset = DEFAULT_CHARSET
      Font.Color = clWindowText
      Font.Height = -13
      Font.Name = 'Arial'
      Font.Style = []
      ImageIndex = 1
      ParentFont = False
      object rgIsen: TRadioGroup
        Left = 21
        Top = 8
        Width = 185
        Height = 177
        Caption = 'Input quantity'
        ItemIndex = 0
        Items.Strings = (
          'Mach number'
          'q ratio, q / p'
          'pressure ratio, p / pt'
          'area ratio, A / A*  (M<1)'
          'area ratio, A / A*  (M>1)'
          'density ratio, rho/rhot'
          'temperature ratio, T / Tt'
          'velocity ratio, V / a*'
          'Prandtl-Meyer Angle, deg'
          'Mach Angle, deg.')
        TabOrder = 3
        OnClick = rgIsenClick
      end
      object gbIsen: TGroupBox
        Left = 21
        Top = 192
        Width = 185
        Height = 44
        Caption = 'Input Value'
        TabOrder = 0
        object ebIsenInput: TEdit
          Left = 14
          Top = 14
          Width = 128
          Height = 24
          TabOrder = 0
          Text = '0'
        end
      end
      object memIsen: TMemo
        Left = 248
        Top = 14
        Width = 289
        Height = 276
        TabStop = False
        TabOrder = 4
      end
      object btnIsenCompute: TButton
        Left = 21
        Top = 241
        Width = 66
        Height = 21
        Caption = 'Compute'
        Default = True
        TabOrder = 1
        OnClick = btnIsenComputeClick
      end
      object btnIsenQuit: TButton
        Left = 107
        Top = 241
        Width = 64
        Height = 21
        Hint = 'quit VuCalc'
        Cancel = True
        Caption = 'Quit'
        ParentShowHint = False
        ShowHint = True
        TabOrder = 2
        OnClick = btnIsenQuitClick
      end
      object btnIsenHelp: TButton
        Left = 28
        Top = 275
        Width = 65
        Height = 22
        Caption = 'Help'
        TabOrder = 5
        OnClick = btnIsenHelpClick
      end
    end
    object TabSheet3: TTabSheet
      Caption = 'normal shock'
      ImageIndex = 2
      object rgNormal: TRadioGroup
        Left = 8
        Top = 8
        Width = 205
        Height = 153
        Caption = 'Input quantity'
        ItemIndex = 0
        Items.Strings = (
          'Upstream Mach number'
          'Downstream Mach number'
          'Total pressure ratio, pt2 / pt1'
          'Static pressure ratio, p2 / p1'
          'Static temperature ratio, T2 / T1'
          'Static density ratio, rho2 / rho1')
        TabOrder = 4
        OnClick = rgNormalClick
      end
      object gbNormal: TGroupBox
        Left = 8
        Top = 176
        Width = 193
        Height = 41
        Caption = 'Input Value'
        TabOrder = 0
        object ebNormalInput: TEdit
          Left = 24
          Top = 16
          Width = 161
          Height = 24
          TabOrder = 0
          Text = '2'
        end
      end
      object btnNormalCompute: TButton
        Left = 16
        Top = 240
        Width = 75
        Height = 25
        Caption = 'Compute'
        Default = True
        TabOrder = 1
        OnClick = btnNormalComputeClick
      end
      object btnNormalQuit: TButton
        Left = 104
        Top = 240
        Width = 75
        Height = 25
        Hint = 'quit VuCalc'
        Cancel = True
        Caption = 'Quit'
        ParentShowHint = False
        ShowHint = True
        TabOrder = 2
        OnClick = btnNormalQuitClick
      end
      object btnNormalHelp: TButton
        Left = 192
        Top = 240
        Width = 75
        Height = 25
        Caption = 'Help'
        TabOrder = 3
        OnClick = btnNormalHelpClick
      end
      object memNormal: TMemo
        Left = 232
        Top = 16
        Width = 270
        Height = 217
        TabStop = False
        TabOrder = 5
      end
    end
    object TabSheet4: TTabSheet
      Caption = 'oblique shock'
      ImageIndex = 3
      object rgOblique: TRadioGroup
        Left = 8
        Top = 8
        Width = 217
        Height = 177
        Caption = 'Input quantity'
        ItemIndex = 0
        Items.Strings = (
          'ramp angle, deg.'
          'shock angle, deg.'
          'total  pressure ratio, pt2 / pt1'
          'static pressure ratio, p2 / p1'
          'static temperature ratio, T2 / T1'
          'density ratio, rho2/rho1'
          'downstream Mach number')
        TabOrder = 4
        OnClick = rgObliqueClick
      end
      object gbObliqueInput: TGroupBox
        Left = 8
        Top = 193
        Width = 217
        Height = 57
        Caption = 'Input Value'
        TabOrder = 0
        object ebObliqueInput: TEdit
          Left = 24
          Top = 16
          Width = 161
          Height = 24
          TabOrder = 0
          Text = '0'
        end
      end
      object btnObliqueCompute: TButton
        Left = 256
        Top = 293
        Width = 75
        Height = 26
        Caption = 'Compute'
        Default = True
        TabOrder = 1
        OnClick = btnObliqueComputeClick
      end
      object btnObliqueQuit: TButton
        Left = 337
        Top = 293
        Width = 75
        Height = 26
        Hint = 'Quit VuCalc'
        Cancel = True
        Caption = 'Quit'
        ParentShowHint = False
        ShowHint = True
        TabOrder = 2
        OnClick = btnObliqueQuitClick
      end
      object btnObliqueHelp: TButton
        Left = 418
        Top = 293
        Width = 75
        Height = 26
        Caption = 'Help'
        TabOrder = 3
        OnClick = btnObliqueHelpClick
      end
      object memOblique: TMemo
        Left = 264
        Top = 16
        Width = 252
        Height = 267
        TabStop = False
        TabOrder = 6
      end
      object gbObliqueUpstream: TGroupBox
        Left = 7
        Top = 261
        Width = 220
        Height = 58
        Caption = 'Upstream Mach number'
        TabOrder = 5
        object ebObliqueM1: TEdit
          Left = 28
          Top = 21
          Width = 156
          Height = 24
          TabOrder = 0
          Text = '2.0'
        end
      end
    end
    object TabSheet5: TTabSheet
      Caption = 'std. atmosphere'
      ImageIndex = 4
      object btnAtmosCompute: TButton
        Left = 109
        Top = 125
        Width = 75
        Height = 26
        Hint = 'compute for the altitude and Mach shown'
        Caption = 'Compute'
        Default = True
        ParentShowHint = False
        ShowHint = True
        TabOrder = 0
        OnClick = btnAtmosComputeClick
      end
      object btnAtmosQuit: TButton
        Left = 109
        Top = 161
        Width = 75
        Height = 25
        Hint = 'quit VuCalc'
        Cancel = True
        Caption = 'Quit'
        ParentShowHint = False
        ShowHint = True
        TabOrder = 1
        OnClick = btnAtmosQuitClick
      end
      object btnAtmosHelp: TButton
        Left = 109
        Top = 196
        Width = 75
        Height = 25
        Caption = 'Help'
        TabOrder = 2
        OnClick = btnAtmosHelpClick
      end
      object memAtmos: TMemo
        Left = 205
        Top = 14
        Width = 283
        Height = 269
        TabOrder = 3
      end
      object gbAtmosAlt: TGroupBox
        Left = 14
        Top = 7
        Width = 170
        Height = 43
        Caption = 'Altitude'
        TabOrder = 4
        object lbAtmosAlt: TLabel
          Left = 141
          Top = 14
          Width = 11
          Height = 16
          Caption = 'm'
        end
        object ebAtmosAlt: TEdit
          Left = 14
          Top = 14
          Width = 121
          Height = 24
          TabOrder = 0
          Text = '0'
        end
      end
      object gbAtmosMach: TGroupBox
        Left = 14
        Top = 64
        Width = 170
        Height = 43
        Caption = 'Mach number'
        TabOrder = 5
        object ebAtmosMach: TEdit
          Left = 14
          Top = 14
          Width = 114
          Height = 24
          Hint = 'Calculate for the altitude and Mach number in the input boxes.'
          ParentShowHint = False
          ShowHint = True
          TabOrder = 0
          Text = '0'
        end
      end
      object rgAtmos: TRadioGroup
        Left = 14
        Top = 120
        Width = 79
        Height = 79
        Caption = 'Units'
        ItemIndex = 1
        Items.Strings = (
          'US'
          'SI')
        TabOrder = 6
        OnClick = rgAtmosClick
      end
    end
    object TabSheet6: TTabSheet
      Caption = 'Rayleigh'
      ImageIndex = 5
      object gbRayleighInput: TGroupBox
        Left = 15
        Top = 198
        Width = 184
        Height = 64
        Caption = 'Input value'
        TabOrder = 0
        object ebRayleighInput: TEdit
          Left = 22
          Top = 28
          Width = 92
          Height = 24
          TabOrder = 0
          Text = '1'
        end
      end
      object memRayleigh: TMemo
        Left = 247
        Top = 16
        Width = 283
        Height = 253
        TabOrder = 1
      end
      object btnRayleighQuit: TButton
        Left = 86
        Top = 274
        Width = 66
        Height = 22
        Hint = 'quit VuCalc'
        Caption = 'Quit'
        ParentShowHint = False
        ShowHint = True
        TabOrder = 2
        OnClick = btnRayleighQuitClick
      end
      object btnRayleighHelp: TButton
        Left = 163
        Top = 272
        Width = 65
        Height = 21
        Caption = 'Help'
        TabOrder = 3
        OnClick = btnRayleighHelpClick
      end
      object btnRayleighCompute: TButton
        Left = 10
        Top = 276
        Width = 65
        Height = 21
        Hint = 'Compute for the input value of the input quantity'
        Caption = 'Compute'
        Default = True
        ParentShowHint = False
        ShowHint = True
        TabOrder = 4
        OnClick = btnRayleighComputeClick
      end
      object rgRayleigh: TRadioGroup
        Left = 14
        Top = 7
        Width = 177
        Height = 184
        Caption = 'Input quantity'
        ItemIndex = 0
        Items.Strings = (
          'Mach Number'
          'To / To*   (M < 1)'
          'T o/ To*   (M > 1)'
          'T / T*   (M < 1)'
          'T / T*   (M > 1)'
          'p / p*'
          'Po/Po*   (M < 1)'
          'Po/Po*   (M > 1)'
          'V / V*'
          'rho / rho*')
        TabOrder = 5
        OnClick = rgRayleighClick
      end
    end
    object TabSheet7: TTabSheet
      Caption = 'Fanno'
      ImageIndex = 6
      object rgFanno: TRadioGroup
        Left = 14
        Top = 14
        Width = 156
        Height = 185
        Caption = 'Input quantity'
        ItemIndex = 0
        Items.Strings = (
          'Mach number'
          'T / T*'
          'p / p*'
          'pt / pt* (M < 1)'
          'pt / pt* (M > 1)'
          'V / V*'
          '(fLmax/D) (M < 1)'
          '(fLmax/D) (M > 1)')
        TabOrder = 0
        OnClick = rgFannoClick
      end
      object gbFannoInput: TGroupBox
        Left = 14
        Top = 212
        Width = 156
        Height = 43
        Caption = 'Input Value'
        TabOrder = 1
        object ebFannoInput: TEdit
          Left = 7
          Top = 14
          Width = 142
          Height = 24
          TabOrder = 0
          Text = '1'
        end
      end
      object btnFannoCompute: TButton
        Left = 14
        Top = 282
        Width = 66
        Height = 22
        Hint = 'compute, using the input value for the input quantity'
        Caption = 'Compute'
        Default = True
        ParentShowHint = False
        ShowHint = True
        TabOrder = 2
        OnClick = btnFannoComputeClick
      end
      object btnFannoQuit: TButton
        Left = 92
        Top = 282
        Width = 66
        Height = 22
        Cancel = True
        Caption = 'Quit'
        TabOrder = 3
        OnClick = btnFannoQuitClick
      end
      object btnFannoHelp: TButton
        Left = 169
        Top = 282
        Width = 67
        Height = 22
        Caption = 'Help'
        TabOrder = 4
        OnClick = btnFannoHelpClick
      end
      object memFanno: TMemo
        Left = 212
        Top = 21
        Width = 311
        Height = 248
        TabOrder = 5
      end
    end
    object TabSheet8: TTabSheet
      Caption = 'Help'
      ImageIndex = 7
      object Memo1: TMemo
        Left = 40
        Top = 40
        Width = 433
        Height = 281
        Lines.Strings = (
          #39'VuCalc is an aid to making calculations in compressible '
          
            'fluid dynamics, such as isentropic flow, normal shock, oblique s' +
            'hock'
          
            'Rayleigh and Fanno flow. A page allows computation of the standa' +
            'rd'
          'atmosphere at a given altitude and Mach number.'
          ''
          'There are several tabs across the top of the form.'
          'If you click on a tab, a computing form appears for the desired'
          
            'calculation. There will be one or more input boxes to be filled ' +
            'in'
          'and choices that can be made about which input quantity is to be'
          'specified.'
          
            'After making your choices, click the Compute button and the desi' +
            'red'
          'results will appear in a window.'
          
            'At any time, you may leave the program by clicking a Quit button' +
            '.'
          
            'The Help buttons display a message appropriate to the active pag' +
            'e.'
          '  ')
        TabOrder = 0
      end
    end
  end
end
