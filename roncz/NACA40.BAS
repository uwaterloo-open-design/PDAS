10000 'PROGRAM TO GENERATE NACA AIRFOIL COORDINATES
10010 '
10020 'COMPILER METACOMMAND
10030 '$DYNAMIC
10040 '
10050 'WRITTEN BY D. LEDNICER, dave@amiwest.com
10060 ' COPYRIGHT 1995 - THIS SOFTWARE IS DECLARED TO BE FREEWARE
10070 '
10080 VERSN = 4
10090 '
10100 'SET INPUT POINT LIMITS AND OTHER CONSTANTS
10110 PI = 4 * ATN(1)
10120 DEF FNR5 (X) = INT(X * 100000! + .5)
10130 DEF FNX (X) = 9 + INT(X * 300 + .5)
10140 DEF FNY (Y) = 95 - INT(Y * 270 + .5)
10150 NP = 100
10160 MPTS = 400
10170 DTP = 3
10180 QM = 640  'HORIZONTAL SCREEN COORDINATES
10190 RM = 350  'VERTICAL SCREEN COORDINATES
10200 STOL = .0001
10210 N6S = 5  'NUMBER OF 6 SERIES IN DATABASE (63 - 67)
10220 N6T = 8  'NUMBER OF 6 SERIES THICKNESS FORMS FOR EACH SUBSERIES
10230 N6C = 11 'NUMBER OF NACA CAMBERLINES IN DATABASE
10240 N6P = 26 'NUMBER OF 6 SERIES X POINTS IN DATABASE
10250 '
10260 WIDTH 80: SCREEN 0, 0, 0: COLOR 7, 1: KEY OFF: CLS
10270 LOCATE 2, 14: COLOR 0, 7
10280 PRINT USING " THE NACA AIRFOIL ORDINATE GENERATION PROGRAM v#.# "; VERSN
10290 COLOR 7, 1: PRINT : PRINT
10300 VIEW PRINT 5 TO 24
10310 '
10320 INPUT " Do You Want To Reset Constants Used In The Program (Y/N)[N]"; CHGE$
10330 IF CHGE$ = "y" OR CHGE$ = "Y" THEN GOSUB 24400
10340 '
10350 DIM XCC(MPTS), XU(MPTS), YU(MPTS), XL(MPTS), YL(MPTS), YT(MPTS), YC(MPTS)
10360 DIM XDD(MPTS), DYC(MPTS), XQ(2), YQ(2), ALP(14)
10370 DIM P(2, NP * 4), A1(3), B1(3), DT(3), X(MPTS), Y(MPTS), YP(MPTS), DY(2)
10380 DIM N(NP * 4, 4), B(3, NP * 4), L(NP * 4), U(3, NP * 4), C(2, NP * 4), R(2), F(4, 2)
10390 DIM X6(N6P), YT6(40, N6P), T6(N6S), YC6(N6C, N6P), DYC6(N6C, N6P)
10400 DIM ENDSLP(2), DYDX(2), XI(30), YI(30), TF6(1 + N6S * N6T + 2 * N6C, N6P)
10410 '
10420 PRINT
10430 PRINT " You Have The Following Options:"
10440 PRINT "     1 - NACA 4 Digit Airfoil"
10450 PRINT "     2 - NACA 5 Digit Airfoil"
10460 PRINT "     3 - NACA Modified 4 Digit Airfoil"
10470 PRINT "     4 - NACA Modified 5 Digit Airfoil"
10480 PRINT "     5 - NACA 16 Series Airfoil"
10490 PRINT "     6 - NACA 6 Series Airfoil"
10500 PRINT "     7 - Symmetrical Biconvex Airfoil"
10510 PRINT ""
10520 INPUT " Your Choice [1]"; OPT
10530 IF OPT = 0 THEN OPT = 1
10540 '
10550 PRINT ""
10560 PRINT " You Have The Following Options For Coordinate Spacing:"
10570 PRINT "     1 - Equal Spacing"
10580 PRINT "     2 - Half Cosine With Smaller Increments Near 0"
10590 PRINT "     3 - Half Cosine With Smaller Increments Near 1"
10600 PRINT "     4 - Full Cosine"
10610 PRINT "     5 - Read X/C Values From A File"
10620 PRINT "     6 - Manual Entry Of X/C Values"
10630 PRINT ""
10640 INPUT " Your Choice [4]"; SPAC
10650 IF SPAC = 0 THEN SPAC = 4
10660 VIEW PRINT 1 TO 24
10670 CLS : PRINT
10680 '
10690 ' MANUAL ENTRY
10700 '
10710 IF SPAC = 6 THEN
10720   PRINT " You Have Chosen Manual Entry"
10730   INPUT " How Many Points Are You Going To Enter"; NPP
10740   PRINT " Ok, Now Enter The Points, The First Being 0 And The Last Being 1"
10750   PRINT ""
10760   FOR I = 1 TO NPP
10770    PRINT USING " Point Number ### = "; I: INPUT XDD(I)
10780   NEXT I
10790   IF XDD(1) <> 0 THEN PRINT " *** THE FIRST POINT IS NOT 0, TRY AGAIN ***": GOTO 10750
10800   IF XDD(NPP) <> 1 THEN PRINT " *** THE LAST POINT IS NOT 1, TRY AGAIN ***": GOTO 10750
10810 END IF
10820 '
10830 ' FILE INPUT
10840 '
10850 IF SPAC = 5 THEN
10860   PRINT
10870   INPUT " Enter Name Of File Containing Values"; XVAL$
10880   OPEN XVAL$ FOR INPUT AS 3
10890   NPP = 0
10900   IF EOF(3) THEN GOTO 10940
10910   NPP = NPP + 1
10920   INPUT #3, XDD(NPP)
10930   GOTO 10900
10940   CLOSE 3
10950 END IF
10960 '
10970 IF SPAC < 5 THEN
10980   CLS : PRINT
10990   INPUT " How Many Points Do You Want Generated"; NPP
11000 END IF
11010 '
11020 ' EQUAL SPACING
11030 '
11040 IF SPAC = 4 THEN
11050   DELTH = PI / (NPP - 1)
11060   FOR I = 1 TO NPP
11070    XDD(I) = .5 - .5 * COS(DELTH * (I - 1))
11080   NEXT I
11090 END IF
11100 '
11110 ' HALF COSINE AFT PACKED
11120 '
11130 IF SPAC = 3 THEN
11140   DELTH = PI / 2 / (NPP - 1)
11150   FOR I = 1 TO NPP
11160    XDD(I) = COS(PI / 2 - DELTH * (I - 1))
11170   NEXT I
11180 END IF
11190 '
11200 ' HALF COSINE FWD PACKED
11210 '
11220 IF SPAC = 2 THEN
11230   DELTH = PI / 2 / (NPP - 1)
11240   FOR I = 1 TO NPP
11250    XDD(I) = 1 - COS(DELTH * (I - 1))
11260   NEXT I
11270 END IF
11280 '
11290 ' FULL COSINE
11300 '
11310 IF SPAC = 1 THEN
11320   DELTH = 1 / (NPP - 1)
11330   FOR I = 1 TO NPP
11340    XDD(I) = DELTH * (I - 1)
11350   NEXT I
11360 END IF
11370 '
11380 ' REJOIN PATHS
11390 '
11400 DELTH = PI / NP
11410 FOR I = 1 TO NP
11420  XCC(I) = .5 - .5 * COS(DELTH * I)
11430 NEXT I
11440 CLS : PRINT
11450 '
11460 ' NACA 4 DIGIT --------------------------------------------------
11470 '
11480 IF OPT > 1 THEN GOTO 11870
11490 PRINT " You Have Chosen To Create A NACA 4 Digit Airfoil"
11500 PRINT ""
11510 INPUT " Enter The First Digit Of The 4 Digit Designation"; MM
11520 INPUT " Enter The Second Digit Of The 4 Digit Designation"; PP
11530 INPUT " Enter The Last Two Digits Of The 4 Digit Designation"; TOC
11540 PRINT ""
11550 IF MM = 0 THEN
11560   NP = NPP
11570   FOR I = 1 TO NP
11580    XCC(I) = XDD(I)
11590   NEXT I
11600 END IF
11610 MC = MM / 100
11620 PC = PP / 10
11630 TC = TOC / 100
11640 FOR I = 1 TO NP
11650  YT(I) = TC * (1.4845 * SQR(XCC(I)) - .63 * XCC(I) - 1.758 * XCC(I) ^ 2 + 1.4215 * XCC(I) ^ 3 - .5075 * XCC(I) ^ 4)
11660  IF MC = 0 THEN
11670    YC(I) = 0
11680    DYC(I) = 0
11690   ELSE
11700    IF XCC(I) > PC THEN
11710      YC(I) = MC / (1 - PC) ^ 2 * (1 - 2 * PC + 2 * PC * XCC(I) - XCC(I) ^ 2)
11720      DYC(I) = 2 * MC / (1 - PC) ^ 2 * (PC - XCC(I))
11730     ELSE
11740      YC(I) = MC / PC ^ 2 * (2 * PC * XCC(I) - XCC(I) ^ 2)
11750      DYC(I) = 2 * MC / PC ^ 2 * (PC - XCC(I))
11760    END IF
11770  END IF
11780 NEXT I
11790 LER = 1.1019 * TC ^ 2
11800 TEANG = 2 * ATN(1.16925 * TC)
11810 DESIG = MM * 1000 + PP * 100 + TOC
11820 DESIG$ = STR$(DESIG)
11830 GOTO 15450
11840 '
11850 ' NACA 5 DIGIT --------------------------------------------------
11860 '
11870 IF OPT > 2 THEN GOTO 12500
11880 PRINT " You Have Chosen To Create A NACA 5 Digit Airfoil"
11890 PRINT ""
11900 INPUT " Enter The First Digit Of The 5 Digit Designation"; LL
11910 INPUT " Enter The Second Digit Of The 5 Digit Designation"; PP
11920 INPUT " Enter The Third Digit Of The 5 Digit Designation"; QQ
11930 INPUT " Enter The Last Two Digits Of The 5 Digit Designation"; TOC
11940 PRINT ""
11950 IF QQ < 0 OR QQ > 1 THEN
11960   PRINT "  THIRD DIGIT OUT OF RANGE, ASSUMING NEAREST VALID VALUE"
11970 END IF
11980 IF LL = 0 THEN
11990   NP = NPP
12000   FOR I = 1 TO NP
12010    XCC(I) = XDD(I)
12020   NEXT I
12030 END IF
12040 LC = LL / 10
12050 PC = PP / 20
12060 TC = TOC / 100
12070 MC = PC / 2
12080 PCT = SQR(MC ^ 2 * (MC / 3 - 2 * SQR(MC / 3) + 1))
12090 IF ABS(PC - PCT) <= .0001 THEN GOTO 12120
12100 MC = MC * PC / PCT
12110 GOTO 12080
12120 SVAL = 1 - 2 * MC
12130 GOSUB 19320
12140 QC = (3 * MC - 7 * MC ^ 2 + 8 * MC ^ 3 - 4 * MC ^ 4) / SQR(MC * (1 - MC)) - 1.5 * (1 - 2 * MC) * (PI / 2 - THETA)
12150 K1 = 9 * LC / QC
12160 K2K1 = (3 * (MC - PC) ^ 2 - MC ^ 3) / (1 - MC) ^ 3
12170 FOR I = 1 TO NP
12180  YT(I) = TC * (1.4845 * SQR(XCC(I)) - .63 * XCC(I) - 1.758 * XCC(I) ^ 2 + 1.4215 * XCC(I) ^ 3 - .5075 * XCC(I) ^ 4)
12190  IF LL = 0 THEN GOTO 12340
12200  IF QQ > 0 THEN GOTO 12280
12210  IF XCC(I) > MC THEN GOTO 12250
12220  YC(I) = K1 / 6 * (XCC(I) ^ 3 - 3 * MC * XCC(I) ^ 2 + MC ^ 2 * (3 - MC) * XCC(I))
12230  DYC(I) = K1 / 6 * (3 * XCC(I) ^ 2 - 6 * MC * XCC(I) + MC ^ 2 * (3 - MC))
12240  GOTO 12340
12250  YC(I) = K1 / 6 * MC ^ 3 * (1 - XCC(I))
12260  DYC(I) = -K1 / 6 * MC ^ 3
12270  GOTO 12340
12280  IF XCC(I) > MC THEN GOTO 12320
12290  YC(I) = K1 / 6 * ((XCC(I) - MC) ^ 3 - K2K1 * (1 - MC) ^ 3 * XCC(I) - MC ^ 3 * XCC(I) + MC ^ 3)
12300  DYC(I) = K1 / 6 * (3 * (XCC(I) - MC) ^ 2 - K2K1 * (1 - MC) ^ 3 - MC ^ 3)
12310  GOTO 12340
12320  YC(I) = K1 / 6 * (K2K1 * (XCC(I) - MC) ^ 3 - K2K1 * (1 - MC) ^ 3 * XCC(I) - MC ^ 3 * XCC(I) + MC ^ 3)
12330  DYC(I) = K1 / 6 * (3 * K2K1 * (XCC(I) - MC) ^ 2 - K2K1 * (1 - MC) ^ 3 - MC ^ 3)
12340 NEXT I
12350 LER = 1.1019 * TC ^ 2
12360 TEANG = 2 * ATN(1.16925 * TC)
12370 DESIG = LL * 10000 + PP * 1000 + QQ * 100 + TOC
12380 DESIG$ = STR$(DESIG)
12390 'INPUT "  Enter Name Of Camberline Output File"; DEST2$
12400 'OPEN DEST2$ FOR OUTPUT AS 3
12410 'FOR I = 1 TO NP
12420 ' PRINT #3, USING " #.####### #.#######"; XCC(I); YC(I)
12430 'NEXT I
12440 'PRINT ""
12450 'CLOSE 3
12460 GOTO 15450
12470 '
12480 ' NACA MODIFIED 4 DIGIT -----------------------------------------
12490 '
12500 IF OPT > 3 THEN GOTO 13000
12510 PRINT " You Have Chosen To Create A NACA Modified 4 Digit Airfoil"
12520 PRINT ""
12530 INPUT " Enter The First Digit Of The 4 Digit Designation"; MM
12540 INPUT " Enter The Second Digit Of The 4 Digit Designation"; PP
12550 INPUT " Enter The Last Two Digits Of The 4 Digit Designation"; TOC
12560 INPUT " Enter The First Digit Of The Appended Number"; IP
12570 INPUT " Enter The Second Digit Of The Appended Number"; TT
12580 PRINT ""
12590 IF MM = 0 THEN
12600   NP = NPP
12610   FOR I = 1 TO NP
12620    XCC(I) = XDD(I)
12630   NEXT I
12640 END IF
12650 MC = MM / 100
12660 PC = PP / 10
12670 TC = TOC / 100
12680 TP = TT / 10
12690 D1 = (2.24 - 5.42 * TP + 12.3 * TP ^ 2) / 10 / (1 - .878 * TP)
12700 D2 = (.294 - 2 * (1 - TP) * D1) / (1 - TP) ^ 2
12710 D3 = (-.196 + (1 - TP) * D1) / (1 - TP) ^ 3
12720 A0 = .296904 * IP / 6
12730 R1 = (1 - TP) ^ 2 / 5 / (.588 - 2 * D1 * (1 - TP))
12740 AA1 = .3 / TP - 15 * A0 / 8 / SQR(TP) - TP / 10 / R1
12750 A2 = -.3 / TP ^ 2 + 5 * A0 / 4 / TP ^ (1.5) + 1 / 5 / R1
12760 A3 = .1 / TP ^ 3 - .375 * A0 / TP ^ (2.5) - 1 / 10 / R1 / TP
12770 FOR I = 1 TO NP
12780  IF XCC(I) > TP THEN GOTO 12810
12790  YT(I) = 5 * TC * (A0 * SQR(XCC(I)) + AA1 * XCC(I) + A2 * XCC(I) ^ 2 + A3 * XCC(I) ^ 3)
12800  GOTO 12820
12810  YT(I) = 5 * TC * (.002 + D1 * (1 - XCC(I)) + D2 * (1 - XCC(I)) ^ 2 + D3 * (1 - XCC(I)) ^ 3)
12820  IF MM = 0 THEN GOTO 12890
12830  IF XCC(I) > PC THEN GOTO 12870
12840  YC(I) = MC / PC ^ 2 * (2 * PC * XCC(I) - XCC(I) ^ 2)
12850  DYC(I) = 2 * MC / PC ^ 2 * (PC - XCC(I))
12860  GOTO 12890
12870  YC(I) = MC / (1 - PC) ^ 2 * (1 - 2 * PC + 2 * PC * XCC(I) - XCC(I) ^ 2)
12880  DYC(I) = 2 * MC / (1 - PC) ^ 2 * (PC - XCC(I))
12890 NEXT I
12900 LER = 1.1019 * (IP / 6 * TC) ^ 2
12910 IF IP >= 9 THEN LER = 3 * 1.1019 * (TC) ^ 2
12920 TEANG = 2 * ATN(1.16925 * TC)
12930 DESIG = MM * 1000 + PP * 1000 + TOC
12940 SESIG = IP * 10 + TT
12950 DESIG$ = STR$(DESIG) + "-" + STR$(SESIG)
12960 GOTO 15450
12970 '
12980 ' NACA MODIFIED 5 DIGIT -----------------------------------------
12990 '
13000 IF OPT > 4 THEN GOTO 13730
13010 PRINT " You Have Chosen To Create A NACA Modified 5 Digit Airfoil"
13020 PRINT ""
13030 INPUT " Enter The First Digit Of The 5 Digit Designation"; LL
13040 INPUT " Enter The Second Digit Of The 5 Digit Designation"; PP
13050 INPUT " Enter The Third Digit Of The 5 Digit Designation"; QQ
13060 INPUT " Enter The Last Two Digits Of The 5 Digit Designation"; TOC
13070 INPUT " Enter The First Digit Of The Appended Number"; IP
13080 INPUT " Enter The Second Digit Of The Appended Number"; TT
13090 PRINT ""
13100 IF QQ < 0 OR QQ > 1 THEN
13110   PRINT "  THIRD DIGIT OUT OF RANGE, ASSUMING NEAREST VALID VALUE"
13120 END IF
13130 IF LL = 0 THEN
13140   NP = NPP
13150   FOR I = 1 TO NP
13160    XCC(I) = XDD(I)
13170   NEXT I
13180 END IF
13190 LC = LL / 10
13200 PC = PP / 20
13210 TC = TOC / 100
13220 MC = PC / 2
13230 TP = TT / 10
13240 D1 = (2.24 - 5.42 * TP + 12.3 * TP ^ 2) / 10 / (1 - .878 * TP)
13250 D2 = (.294 - 2 * (1 - TP) * D1) / (1 - TP) ^ 2
13260 D3 = (-.196 + (1 - TP) * D1) / (1 - TP) ^ 3
13270 A0 = .296904 * IP / 6
13280 R1 = (1 - TP) ^ 2 / 5 / (.588 - 2 * D1 * (1 - TP))
13290 AA1 = .3 / TP - 15 * A0 / 8 / SQR(TP) - TP / 10 / R1
13300 A2 = -.3 / TP ^ 2 + 5 * A0 / 4 / TP ^ (1.5) + 1 / 5 / R1
13310 A3 = .1 / TP ^ 3 - .375 * A0 / TP ^ (2.5) - 1 / 10 / R1 / TP
13320 PCT = SQR(MC ^ 2 * (MC / 3 - 2 * SQR(MC / 3) + 1))
13330 IF ABS(PC - PCT) <= .0001 THEN GOTO 13360
13340 MC = MC * PC / PCT
13350 GOTO 13320
13360 SVAL = 1 - 2 * MC
13370 GOSUB 19320
13380 QC = (3 * MC - 7 * MC ^ 2 + 8 * MC ^ 3 - 4 * MC ^ 4) / SQR(MC * (1 - MC)) - 1.5 * (1 - 2 * MC) * (PI / 2 - THETA)
13390 K1 = 9 * LC / QC
13400 K2K1 = (3 * (MC - PC) ^ 2 - MC ^ 3) / (1 - MC) ^ 3
13410 FOR I = 1 TO NP
13420  IF XCC(I) > TP THEN GOTO 13450
13430  YT(I) = 5 * TC * (A0 * SQR(XCC(I)) + AA1 * XCC(I) + A2 * XCC(I) ^ 2 + A3 * XCC(I) ^ 3)
13440  GOTO 13460
13450  YT(I) = 5 * TC * (.002 + D1 * (1 - XCC(I)) + D2 * (1 - XCC(I)) ^ 2 + D3 * (1 - XCC(I)) ^ 3)
13460  IF LL = 0 THEN GOTO 13610
13470  IF QQ > 0 THEN GOTO 13550
13480  IF XCC(I) > MC THEN GOTO 13520
13490  YC(I) = K1 / 6 * (XCC(I) ^ 3 - 3 * MC * XCC(I) ^ 2 + MC ^ 2 * (3 - MC) * XCC(I))
13500  DYC(I) = K1 / 6 * (3 * XCC(I) ^ 2 - 6 * MC * XCC(I) + MC ^ 2 * (3 - MC))
13510  GOTO 13610
13520  YC(I) = K1 / 6 * MC ^ 3 * (1 - XCC(I))
13530  DYC(I) = -K1 / 6 * MC ^ 3
13540  GOTO 13610
13550  IF XCC(I) > MC THEN GOTO 13590
13560  YC(I) = K1 / 6 * ((XCC(I) - MC) ^ 3 - K2K1 * (1 - MC) ^ 3 * XCC(I) - MC ^ 3 * XCC(I) + MC ^ 3)
13570  DYC(I) = K1 / 6 * (3 * (XCC(I) - MC) ^ 2 - K2K1 * (1 - MC) ^ 3 - MC ^ 3)
13580  GOTO 13610
13590  YC(I) = K1 / 6 * (K2K1 * (XCC(I) - MC) ^ 3 - K2K1 * (1 - MC) ^ 3 * XCC(I) - MC ^ 3 * XCC(I) + MC ^ 3)
13600  DYC(I) = K1 / 6 * (3 * K2K1 * (XCC(I) - MC) ^ 2 - K2K1 * (1 - MC) ^ 3 - MC ^ 3)
13610 NEXT I
13620 LER = 1.1019 * (IP / 6 * TC) ^ 2
13630 IF IP >= 9 THEN LER = 3 * 1.1019 * (TC) ^ 2
13640 TEANG = 2 * ATN(1.16925 * TC)
13650 DESIG = LL * 10000 + PP * 1000 + QQ * 100 + TOC
13660 DESIG$ = STR$(DESIG)
13670 SESIG = IP * 10 + TT
13680 DESIG$ = STR$(DESIG) + "-" + STR$(SESIG)
13690 GOTO 15450
13700 '
13710 ' NACA 16 SERIES ------------------------------------------------
13720 '
13730 IF OPT > 5 THEN GOTO 14180
13740 PRINT " You Have Chosen To Create A NACA 16 Series Airfoil"
13750 PRINT ""
13760 INPUT " Enter The First Two Digits Of The Designation"; FSTTWO
13770 INPUT " Enter The Third Digit Of The Designation"; MM
13780 INPUT " Enter The Last Two Digits Of The Designation"; TOC
13790 PRINT ""
13800 IF MM = 0 THEN
13810   NP = NPP
13820   FOR I = 1 TO NP
13830    XCC(I) = XDD(I)
13840   NEXT I
13850 END IF
13860 IP = 4
13870 T = 5
13880 MC = MM / 10
13890 TC = TOC / 100
13900 TP = T / 10
13910 D1 = (2.24 - 5.42 * TP + 12.3 * TP ^ 2) / 10 / (1 - .878 * TP)
13920 D2 = (.294 - 2 * (1 - TP) * D1) / (1 - TP) ^ 2
13930 D3 = (-.196 + (1 - TP) * D1) / (1 - TP) ^ 3
13940 A0 = .296904 * IP / 6
13950 R1 = (1 - TP) ^ 2 / 5 / (.588 - 2 * D1 * (1 - TP))
13960 AA1 = .3 / TP - 15 * A0 / 8 / SQR(TP) - TP / 10 / R1
13970 A2 = -.3 / TP ^ 2 + 5 * A0 / 4 / TP ^ (1.5) + 1 / 5 / R1
13980 A3 = .1 / TP ^ 3 - .375 * A0 / TP ^ (2.5) - 1 / 10 / R1 / TP
13990 FOR I = 1 TO NP
14000  IF XCC(I) > TP THEN GOTO 14030
14010  YT(I) = 5 * TC * (A0 * SQR(XCC(I)) + AA1 * XCC(I) + A2 * XCC(I) ^ 2 + A3 * XCC(I) ^ 3)
14020  GOTO 14040
14030  YT(I) = 5 * TC * (.002 + D1 * (1 - XCC(I)) + D2 * (1 - XCC(I)) ^ 2 + D3 * (1 - XCC(I)) ^ 3)
14040  IF MM = 0 THEN GOTO 14080
14050  IF XCC(I) = 0 OR XCC(I) = 1 THEN GOTO 14080
14060  YC(I) = -MC / 4 / PI * ((1 - XCC(I)) * LOG(1 - XCC(I)) + XCC(I) * LOG(XCC(I)))
14070  DYC(I) = -MC / 4 / PI * (-1 / (1 - XCC(I)) - LOG(1 - XCC(I)) + XCC(I) / (1 - XCC(I)) + LOG(XCC(I)) + 1)
14080 NEXT I
14090 LER = 1.1019 * (IP / 6 * TC) ^ 2
14100 IF IP >= 9 THEN LER = 3 * 1.1019 * (TC) ^ 2
14110 TEANG = 2 * ATN(1.16925 * TC)
14120 DESIG = FSTTWO * 1000 + MM * 100 + TOC
14130 DESIG$ = STR$(DESIG)
14140 GOTO 15450
14150 '
14160 ' NACA 6 SERIES ------------------------------------------------
14170 '
14180 IF OPT > 6 THEN GOTO 15220
14190 PRINT " You Have Chosen To Create A NACA 6 Series Airfoil"
14200 PRINT ""
14210 INPUT " Enter The First Digit Of The 6 Series Designation"; LL
14220 INPUT " Enter The Second Digit Of The 6 Series Designation"; PP
14230 INPUT " Enter The Third Digit Of The 6 Series Designation"; MM
14240 INPUT " Enter The Last Two Digits Of The 6 Series Designation"; TOC
14250 IF MM <> 0 THEN
14260   INPUT " Enter Value Of a For Camberline (0-1)[1]"; AA$
14270   IF AA$ = "" THEN AA = 1 ELSE AA = VAL(AA$)
14280 END IF
14290 INPUT " Do You Want The Lower Surface Cove Filled In (Y/N)[N]"; AMOD$
14300 IF AMOD$ = "" OR AMOD$ = "N" OR AMOD$ = "n" THEN AMOD$ = "N" ELSE AMOD$ = "Y"
14310 PRINT
14320 IF PP < 3 OR PP > 7 THEN
14330   PRINT "  YOU HAVE ASKED FOR A 6 SERIES AIRFOIL OUTSIDE THE THICKNESS"
14340   PRINT "   FORM DATABASE.  RESULTS WILL NOT NECESSARILY BE VALID"
14350   PRINT
14360 END IF
14370 GOSUB 22320
14380 '
14390 NP = N6P
14400 MC = MM / 10
14410 TC = TOC / 100
14420 FOR I = 1 TO NP
14430  XCC(I) = X6(I)
14440  IF TC < .06 THEN
14450    JJ = 1
14460    FOR J = 1 TO N6S
14470     T6(J) = YT6(JJ, I) * TC / .06
14480     JJ = JJ + N6T
14490    NEXT J
14500   ELSEIF TC > .21 THEN
14510    JJ = N6T
14520    FOR J = 1 TO N6S
14530     T6(J) = YT6(JJ, I) * TC / .21
14540     JJ = JJ + N6T
14550    NEXT J
14560   ELSE
14570    JJ = 1
14580    FOR J = 1 TO N6S
14590     N = N6T
14600     XI(1) = .06: XI(2) = .08: XI(3) = 9.000001E-02: XI(4) = .1
14610     XI(5) = .12: XI(6) = .15: XI(7) = .18: XI(8) = .21
14620     FOR K = 1 TO N
14630      YI(K) = YT6(JJ + K - 1, I)
14640     NEXT K
14650     XD = TC
14660     GOSUB 23300
14670     T6(J) = YD
14680     JJ = JJ + N6T
14690    NEXT J
14700  END IF
14710  FOR J = 1 TO N6S
14720   XI(J) = J + 2
14730   YI(J) = T6(J)
14740  NEXT J
14750  N = N6S
14760  XD = PP
14770  GOSUB 23300
14780  YT(I) = YD
14790 '
14800  IF MC = 0 THEN
14810    YC(I) = 0
14820    DYC(I) = 0
14830   ELSE
14840    N = N6C
14850    FOR K = 1 TO N
14860     XI(K) = (K - 1) / 10
14870     YI(K) = YC6(K, I)
14880    NEXT K
14890    XD = AA
14900    GOSUB 23300
14910    YC(I) = YD * MC
14920 '
14930    N = N6C
14940    FOR K = 1 TO N
14950     XI(K) = (K - 1) / 10
14960     YI(K) = DYC6(K, I)
14970    NEXT K
14980    XD = AA
14990    GOSUB 23300
15000    DYC(I) = YD * MC
15010 END IF
15020 '
15030 NEXT I
15040 DESIG = LL * 10 + PP
15050 SESIG = MM * 100 + TOC
15060 IF AMOD$ = "N" THEN
15070   DESIG$ = STR$(DESIG) + "-" + STR$(SESIG)
15080  ELSE
15090   DESIG$ = STR$(DESIG) + "A" + STR$(SESIG)
15100 END IF
15110 'INPUT "  Enter Name Of Camberline Output File"; DEST2$
15120 'OPEN DEST2$ FOR OUTPUT AS 3
15130 'FOR I = 1 TO NP
15140 ' PRINT #3, USING " #.####### #.#######"; XCC(I); YC(I)
15150 'NEXT I
15160 'PRINT ""
15170 'CLOSE 3
15180 GOTO 15450
15190 '
15200 ' BICONVEX AIRFOIL ---------------------------------------------
15210 '
15220 IF OPT > 7 THEN GOTO 10300
15230 PRINT " You Have Chosen To Create A Biconvex Airfoil"
15240 PRINT ""
15250 INPUT " Enter The Desired Thickness To Chord Ratio"; TOC
15260 PRINT ""
15270 RAD = (2 + 8 * (TOC / 2) ^ 2) / 16 / (TOC / 2)
15280 NP = NPP
15290 FOR I = 1 TO NP
15300  XCC(I) = XDD(I)
15310 NEXT I
15320 FOR I = 1 TO NP
15330  YT(I) = SQR(RAD ^ 2 - (XCC(I) - .5) ^ 2)
15340  YC(I) = 0
15350  DYC(I) = 0
15360 NEXT I
15370 Y0 = YT(1)
15380 FOR I = 1 TO NP
15390  YT(I) = YT(I) - Y0
15400 NEXT I
15410 DESIG$ = "BICONVEX"
15420 GOTO 15450
15430 '
15440 ' REJOIN PATHS ----------------------------------------------------
15450 '
15460 FOR I = 1 TO NP
15470  THET = ATN(DYC(I))
15480  XU(I) = XCC(I) - YT(I) * SIN(THET)
15490  YU(I) = YC(I) + YT(I) * COS(THET)
15500  XL(I) = XCC(I) + YT(I) * SIN(THET)
15510  YL(I) = YC(I) - YT(I) * COS(THET)
15520 NEXT I
15530 IF AMOD$ = "Y" THEN GOSUB 23920
15540 INPUT " Do You Want To Bypass Point Redistribution (Y/N)[N]"; ANSD$
15550 IF ANSD$ = "Y" OR ANSD$ = "y" THEN
15560   NPP = NP
15570   GOTO 15600
15580 END IF
15590 GOSUB 19750
15600 CLS : LOCATE 12, 18
15610 INPUT " Would You Like To View The Airfoil (Y/N)[Y]"; WOW$
15620 IF WOW$ = "N" OR WOW$ = "n" THEN GOTO 15720
15630 LOCATE 14, 18
15640 INPUT " Connect Points With Lines (Y/N)[N]"; ANSC$
15650 IF ANSC$ = "" OR ANSC$ = "n" OR ANSC$ = "N" THEN ANSC$ = "N" ELSE ANSC$ = "Y"
15660 GOSUB 19390
15670 '
15680 TITLE$ = "NACA " + DESIG$ + " AIRFOIL"
15690 '
15700 ' OUTPUT
15710 '
15720 CLS : PRINT
15730 INPUT " Enter The Filename For Output Of The Coordinates"; DEST$
15740 OPEN DEST$ FOR OUTPUT AS #1
15750 PRINT : PRINT
15760 PRINT "                 You Have The Following OUTPUT Format Options:"
15770 PRINT "                              1 - LED Format"
15780 PRINT "                              2 - ISES Format"
15790 PRINT "                              3 - Eppler Format"
15800 PRINT "                              4 - WIBCO Format"
15810 PRINT "                              5 - VSAERO Format"
15820 PRINT "                              6 - AutoCAD SCR Format"
15830 INPUT "                             Your Choice [1]"; OOPT
15840 ON OOPT + 1 GOSUB 15880, 15880, 16180, 16380, 18600, 18890, 19020
15850 CLOSE #1
15860 END
15870 '-----------------------------------------------------------------
15880 'LED FORMAT OUTPUT
15890 '
15900 IF OPT = 1 THEN
15910   PRINT #1, USING " NACA 4-DIGIT AIRFOIL M = #.#, P = #.#, T/C = ##.#"; MM; PP; TOC
15920  ELSEIF OPT = 2 THEN
15930   PRINT #1, USING " NACA 5-DIGIT AIRFOIL L = #.#, P = #.#, Q = #.#, T/C = ##.#"; LL; PP; QQ; TOC
15940  ELSEIF OPT = 3 THEN
15950   PRINT #1, USING " NACA MODIFIED 4-DIGIT AIRFOIL M = #.#, P = #.#, T/C = ##.# IP = #.#, T = #.#"; MM; PP; TOC; IP; TT
15960  ELSEIF OPT = 4 THEN
15970   PRINT #1, USING " NACA MODIFIED 5-DIGIT AIRFOIL L = #.#, P = #.#, Q = #.#, T/C = ##.# IP = #.#, T = #.#"; LL; PP; QQ; TOC; IP; TT
15980  ELSEIF OPT = 5 THEN
15990   PRINT #1, USING " NACA 16 SERIES AIRFOIL M = #.#, T/C = ##.#"; MM; TOC
16000  ELSEIF OPT = 6 THEN
16010   IF AMOD$ <> "Y" THEN
16020     PRINT #1, USING " NACA 6 SERIES AIRFOIL P = #.#, M = #.#, T/C = ##.#"; PP; MM; TOC
16030    ELSE
16040     PRINT #1, USING " NACA 6 SERIES AIRFOIL P = #.#, M = #.#, T/C = ##.#, COVE FILLED"; PP; MM; TOC
16050   END IF
16060 END IF
16070 PRINT #1, USING "     ###.      ###."; NPP; NPP
16080 PRINT #1, ""
16090 FOR I = 1 TO NPP
16100  PRINT #1, USING " #.####### #.#######"; XU(I); YU(I)
16110 NEXT I
16120 PRINT #1, ""
16130 FOR I = 1 TO NPP
16140  PRINT #1, USING " #.####### #.#######"; XL(I); YL(I)
16150 NEXT I
16160 RETURN
16170 '-----------------------------------------------------------------
16180 'ISES FORMAT OUTPUT
16190 '
16200 CLS : PRINT : PRINT
16210 INPUT " Enter Length Of Computational Domain Upstream (~2)[2]"; CHINL
16220 IF CHINL = 0 THEN CHINL = 2
16230 INPUT " Enter Length Of Computational Domain Downstream (~2)[2]"; CHOUT
16240 IF CHOUT = 0 THEN CHOUT = 2
16250 INPUT " Enter Approximate Mach Number For Analysis [0.2]"; MACH
16260 IF MACH = 0 THEN MACH = .2
16270 CHWID = 6 / SQR(1 - MACH ^ 2)
16280 PRINT #1, TITLE$
16290 PRINT #1, USING "  0.0  0.0 ##.# ##.# ##.#"; CHINL; CHOUT; CHWID
16300 FOR I = NPP TO 1 STEP -1
16310 PRINT #1, USING " #.####### #.#######"; XU(I); YU(I)
16320 NEXT I
16330 FOR I = 2 TO NPP
16340 PRINT #1, USING " #.####### #.#######"; XL(I); YL(I)
16350 NEXT I
16360 RETURN
16370 '-----------------------------------------------------------------
16380 'EPPLER FORMAT OUTPUT
16390 '
16400 PRINT #1, "FXPR"
16410 PRINT #1, TITLE$
16420 PRINT #1, USING "   ##   ##"; NPP; NPP
16430 FORMAT$ = "#.####### #.####### #.####### #.####### #.####### #.####### #.####### #.#######"
16440 FOR I = 1 TO NPP
16450  PRINT #1, USING FORMAT$; XU(I);
16460  GOSUB 18520
16470 NEXT I
16480 GOSUB 18570
16490 FOR I = 1 TO NPP
16500  PRINT #1, USING FORMAT$; YU(I);
16510  GOSUB 18520
16520  NEXT I
16530 GOSUB 18570
16540 FOR I = 1 TO NPP
16550  PRINT #1, USING FORMAT$; XL(I);
16560  GOSUB 18520
16570 NEXT I
16580 GOSUB 18570
16590 FOR I = 1 TO NPP
16600  PRINT #1, USING FORMAT$; YL(I);
16610  GOSUB 18520
16620 NEXT
16630 GOSUB 18570
16640 L% = 3
16650 WHILE ADD$(L%) <> ""
16660   PRINT #1, ADD$(L%)
16670   L% = L% + 1
16680 WEND
16690 '
16700 CLS : PRINT : PRINT
16710 'MACH CARD
16720 INPUT " Enter Mach Number For This Case [0]"; MACHN
16730 MN = INT(MACHN * 1000)
16740 IF MN <> 0 THEN
16750   PRINT #1, USING "MACH       ####"; MN
16760 END IF
16770 'ALFA CARD
16780 PRINT
16790 INPUT " Enter Number Of Angles Of Attack To Be Entered (14 max)[1]"; NALP
16800 IF NALP < 1 THEN GOTO 18480
16810 IF NALP <> 0 THEN
16820   FOR I = 1 TO NALP
16830    PRINT "  Angle Of Attack #"; I; : INPUT ALP(I)
16840   NEXT I
16850 END IF
16860 PRINT #1, "ALFA";
16870 PRINT
16880 PRINT " Is A Printout Of Moment Coefficient Required"
16890 INPUT "                For Each Angle Of Attack (Y/N)[Y]"; MOM$
16900 IF MOM$ <> "N" AND MOM$ <> "n" THEN NUPA = 1 ELSE NUPA = 2
16910 PRINT #1, USING "#"; NUPA;
16920 NUPE = 1
16930 PRINT #1, USING "#"; NUPE;
16940 PRINT
16950 PRINT " Are The Angles Of Attack Relative To:"
16960 PRINT "  1) The Chord Line"
16970 PRINT "  2) The Zero-Lift Line"
16980 INPUT " Option (1/2)[1]"; AREF
16990 IF AREF <> 2 THEN AREF = 1
17000 PRINT
17010 PRINT " What Should Be Printed Out (And Plotted):"
17020 PRINT "  1) Velocity Ratios v/V"
17030 PRINT "  2) Pressure Coefficients Cp"
17040 INPUT " Option (1/2)[2]"; PPR
17050 IF PPR <> 1 THEN PPR = 2
17060 IF AREF = 1 AND PPR = 1 THEN NUPI = 1
17070 IF AREF = 2 AND PPR = 1 THEN NUPI = 0
17080 IF AREF = 1 AND PPR = 2 THEN NUPI = 3
17090 IF AREF = 2 AND PPR = 2 THEN NUPI = 2
17100 PRINT #1, USING "#"; NUPI;
17110 PRINT #1, USING "###"; NALP;
17120 FOR I = 1 TO NALP
17130  PRINT #1, USING "#####"; INT(ALP(I) * 100);
17140 NEXT I
17150 PRINT #1, ""
17160 'DIAG CARD
17170 PRINT
17180 IF PPR = 1 THEN
17190   INPUT " Is A Plot Of The Velocity Distributions Required (Y/N)[N]", PL$
17200  ELSE
17210   INPUT " Is A Plot Of The Cp Distributions Required (Y/N)[N]", PL$
17220 END IF
17230 IF PL$ = "y" OR PL$ = "Y" THEN
17240   PRINT #1, USING "DIAG  #"; PPR
17250 END IF
17260 'RE CARD
17270 PRINT
17280 INPUT " How Many Combinations Of Re And Transition Type Are Required (5 max)[1]"; NCOMBO
17290 IF NCOMBO < 1 THEN NCOMBO = 1
17300 INPUT "  Print Just A Summary Or The Entire Boundary Layer History (1/2)[1]"; PRTO
17310 IF PRTO = 2 THEN
17320   PRINT " At Each Step In The Boundary Layer Print:"
17330   PRINT "  1) Momentum Thickness"
17340   PRINT "  2) Reynolds Number Based Upon Momentum Thickness"
17350   PRINT "  3) Displacement Thickness"
17360   INPUT " Option (1/2/3)[1]"; BLO
17370   IF BLO < 1 OR BLO > 3 THEN BLO = 1
17380 END IF
17390 RST$ = "RE  1"
17400 DUM$ = STR$((BLO + 1) * 10)
17410 NUPE$ = MID$(DUM$, 2, 1)
17420 RST$ = RST$ + NUPE$ + "    "
17430 GOSUB 24240
17440 ITRANS = 0
17450 FOR I = 1 TO NCOMBO
17460  PRINT " For Case #"; I
17470  INPUT "  Enter Reynolds Number"; REY
17480  REYN = INT(REY / 1000)
17490  PRINT "  Enter Transition Type (1/2/3/4)[3]"
17500  INPUT "       (Enter ? To See Options Again)"; TTYPE$
17510  IF TTYPE$ = "?" THEN
17520    GOSUB 24240
17530    GOTO 17490
17540  END IF
17550  TTYPE = VAL(TTYPE$)
17560  IF TTYPE < 1 OR TTYPE > 4 THEN TTYPE = 3
17570  IF TTYPE = 1 THEN
17580    MU = 0
17590  END IF
17600  IF TTYPE = 2 THEN
17610    ITRANS = ITRANS + 1
17620    IF ITRANS > 2 THEN
17630      PRINT "  ONLY TWO FIXED TRANSITION CASES ALLOWED"
17640      ITRANS = 2
17650      GOTO 17490
17660     ELSE
17670      INPUT "  Enter X/C Location of Transition On Upper Surface"; XTU(ITRANS)
17680      INPUT "  Enter X/C Location of Transition On Lower Surface"; XTL(ITRANS)
17690      MU = ITRANS
17700    END IF
17710  END IF
17720  IF TTYPE = 3 THEN
17730    MU = 3
17740  END IF
17750  IF TTYPE = 4 THEN
17760    INPUT "  Enter Roughness Level (0-6, ? For Options)[0]"; ROUGH$
17770    IF ROUGH$ = "?" THEN
17780      GOSUB 24240
17790      GOTO 17760
17800    END IF
17810    ROUGH = VAL(ROUGH$)
17820    MU = INT(ROUGH + 3)
17830  END IF
17840  DUM$ = STR$(MU)
17850  MU$ = " " + MID$(DUM$, 2, 1) + "   "
17860  RST$ = RST$ + MU$
17870  DUM$ = STR$(REYN)
17880  LD = LEN(DUM$)
17890  REY$ = MID$(DUM$, 2, LD - 1)
17900  LR = LEN(REY$)
17910  IF LR <> 5 THEN
17920    FOR J = 1 TO LR
17930     REY$ = " " + REY$
17940    NEXT J
17950  END IF
17960  RST$ = RST$ + REY$
17970 NEXT I
17980 IF ITRANS > 0 THEN
17990   LR = LEN(RST$)
18000   DR = 60 - LR
18010   FOR J = 1 TO DR
18020    RST$ = RST$ + " "
18030   NEXT J
18040   FOR J = 1 TO ITRANS
18050    DUM$ = STR$(INT(XTU(ITRANS) * 10000))
18060    LD = LEN(DUM$)
18070    XT$ = MID$(DUM$, 2, LD - 1)
18080    LD = LEN(XT$)
18090    DL = 5 - LD
18100    FOR K = 1 TO DL
18110     RST$ = RST$ + " "
18120    NEXT K
18130    RST$ = RST$ + XT$
18140    DUM$ = STR$(INT(XTL(ITRANS) * 10000))
18150    LD = LEN(DUM$)
18160    XT$ = MID$(DUM$, 2, LD - 1)
18170    LD = LEN(XT$)
18180    DL = 5 - LD
18190    FOR K = 1 TO DL
18200     RST$ = RST$ + " "
18210    NEXT K
18220    RST$ = RST$ + XT$
18230   NEXT J
18240 END IF
18250 PRINT #1, RST$
18260 PRINT
18270 INPUT " Is A CD-CL Plot Required (Y/N)[Y]"; PL$
18280 IF PL$ <> "n" AND PL$ <> "N" THEN
18290   PRINT #1, "CDCL"
18300 END IF
18310 PRINT
18320 INPUT " Is A Flap Deflection Case Required (Y/N)[N]"; FLP$
18330 IF FLP$ = "Y" OR FLP$ = "y" THEN
18340   INPUT "  Enter X/C Location Of Hinge"; XLOC
18350   INPUT "  Enter Y/C Location Of Hinge"; YLOC
18360   INPUT "  Enter Flap Deflection (deg, Negative Is Upwards)"; DEFL
18370   FST$ = "FLAP      "
18380   DUM$ = STR$(INT((1 - XLOC) * 10000))
18390   FLC$ = LEFT$(DUM$, 5)
18400   DUM$ = STR$(INT(YLOC * 10000))
18410   YLC$ = LEFT$(DUM$, 5)
18420   TRS$ = "  400"
18430   DUM$ = STR$(INT(DEFL * 100))
18440   DFL$ = LEFT$(DUM$, 5)
18450   FST$ = FST$ + FLC$ + YLC$ + TRS$ + DFL$
18460   PRINT #1, RST$
18470 END IF
18480 PRINT #1, "ENDE"
18490 RETURN
18500 '-----------------------------------------------------------------
18510 'EPPLER OUTPUT UTILITY SUBROUTINE
18520 IF I / 8 = INT(I / 8) THEN PRINT #1, ""
18530 RETURN
18540 '-----------------------------------------------------------------
18550 'EPPLER OUTPUT UTILITY SUBROUTINE
18560 I = I - 1
18570 IF I / 8 <> INT(I / 8) THEN PRINT #1, ""
18580 RETURN
18590 '-----------------------------------------------------------------
18600 'WIBCO FORMAT OUTPUT
18610 '
18620 GOSUB 18190
18630 PRINT #1, USING "     ##."; NPP
18640 FORMAT$ = " #.###### #.###### #.###### #.###### #.###### #.######"
18650 FOR I = 1 TO NPP
18660  PRINT #1, USING FORMAT$; XU(I);
18670  GOSUB 18820
18680 NEXT I
18690 GOSUB 18850
18700 FOR I = 1 TO NPP
18710  PRINT #1, USING FORMAT$; YU(I);
18720  GOSUB 18820
18730 NEXT I
18740 GOSUB 18850
18750 FOR I = 1 TO NPP
18760  PRINT #1, USING FORMAT$; YL(I);
18770  GOSUB 18820
18780 NEXT I
18790 GOSUB 18850
18800 RETURN
18810 '
18820 IF I / 7 = INT(I / 7) THEN PRINT #1, ""
18830 RETURN
18840 '
18850 I = I - 1
18860 IF I / 7 <> INT(I / 7) THEN PRINT #1, ""
18870 RETURN
18880 '-----------------------------------------------------------------
18890 'VSAERO FORMAT OUTPUT
18900 '
18910 PRINT #1, "    0.0000    0.0000    0.0000    1.0000    0.0000    0.0000    2    0    0    0"
18920 FOR I = NPP TO 1 STEP -1
18930  PRINT #1, USING " #.####### #.#######"; XL(I); YL(I)
18940 NEXT I
18950 PRINT
18960 FOR I = 2 TO NPP
18970  PRINT #1, USING " #.####### #.#######"; XU(I); YU(I)
18980 NEXT I
18990 PRINT #1, "                                  3    0    0"
19000 RETURN
19010 '-----------------------------------------------------------------
19020 'AUTOCAD SCRIPT FORMAT OUTPUT
19030 PRINT
19040 IF RIGHT$(DEST$, 3) <> "SCR" AND RIGHT$(DEST$, 3) <> "scr" THEN
19050   PRINT " ERROR: Output Filename Must Have A .SCR Extension"
19060   CLOSE 1
19070   INPUT " Enter The Filename For Output Of The Coordinates"; DEST$
19080   OPEN DEST$ FOR OUTPUT AS #1
19090 END IF
19100 '
19110 'TURN OFF LIMITS TO PREVENT LOCKUP
19120 PRINT #1, "LIMITS OFF"
19130 'PLACE NODE POINT AT CHORD LINE LEADING EDGE
19140 PRINT #1, "POINT 0,0"
19150 'DRAW CHORD LINE
19160 PRINT #1, "LINE 0,0 1,0 "
19170 'DRAW UPPER SURFACE STARTING AT 0,0
19180 PRINT #1, "PLINE"
19190 FOR I = 1 TO NPP
19200  PRINT #1, USING "#.#######,#.#######"; XU(I); YU(I)
19210 NEXT I
19220 'DRAW LOWER SURFACE STARTING AT 0,0
19230 PRINT #1, " PLINE"
19240 FOR I = 1 TO NPP
19250  PRINT #1, USING "#.#######,#.#######"; XL(I); YL(I)
19260 NEXT I
19270 'RESTORE LIMITS
19280 PRINT #1, " LIMITS ON"
19290 RETURN
19300 '-----------------------------------------------------------------
19310 'ASIN SUBROUTINE
19320 IF SVAL = 0 THEN THETA = 0
19330 IF SVAL = 1 THEN THETA = PI / 2
19340 IF SVAL = -1 THEN THETA = -PI / 2
19350 IF SVAL = 0 OR SVAL = 1 OR SVAL = -1 THEN GOTO 19370
19360 THETA = ATN(SVAL / SQR(1 - SVAL ^ 2))
19370 RETURN
19380 '-----------------------------------------------------------------
19390 'SCREEN DISPLAY SUBROUTINE
19400 KEY OFF
19410 IF DTP = 1 THEN SCREEN 1
19420 IF DTP > 1 THEN SCREEN 9
19430 COLOR 1, 0: CLS
19440 LINE (.05 * QM, .25 * RM)-(.05 * QM, .75 * RM)
19450 LINE (.05 * QM, .5 * RM)-(.95 * QM, .5 * RM)
19460 COLR = 7
19470 FOR I = 1 TO NPP
19480  XQ(1) = .05 * QM + XU(I) * QM * .9
19490  YQ(1) = .5 * RM - YU(I) * QM * .9
19500  IF ANSC$ = "N" THEN
19510    PSET (XQ(1), YQ(1)), COLR
19520    GOTO 19580
19530  END IF
19540  IF I = 1 THEN GOTO 19560
19550  LINE (XQ(2), YQ(2))-(XQ(1), YQ(1)), COLR
19560  XQ(2) = XQ(1)
19570  YQ(2) = YQ(1)
19580 NEXT I
19590 FOR I = 1 TO NPP
19600  XQ(1) = .05 * QM + XL(I) * QM * .9
19610  YQ(1) = .5 * RM - YL(I) * QM * .9
19620  IF ANSC$ = "N" THEN
19630    PSET (XQ(1), YQ(1)), COLR
19640    GOTO 19700
19650  END IF
19660  IF I = 1 THEN GOTO 19680
19670  LINE (XQ(2), YQ(2))-(XQ(1), YQ(1)), COLR
19680  XQ(2) = XQ(1)
19690  YQ(2) = YQ(1)
19700 NEXT I
19710 IF INKEY$ = "" THEN 19710 ELSE 19720
19720 WIDTH 80: SCREEN 0, 0, 0: COLOR 7, 1: CLS
19730 RETURN
19740 '-----------------------------------------------------------------
19750 'POINT REDISTRIBUTION ROUTINE
19760 'set up constants
19770 DEG.TO.RAD = PI / 180
19780 RAD.TO.DEG = 180 / PI
19790 N.TOT = NP * 2 - 1
19800 NUP = NP
19810 NLO = NP
19820 S = 2
19830 C1 = 1
19840 C2 = 1
19850 COLOR 7, 1: CLS
19860 II = 0
19870 FOR I = NUP TO 1 STEP -1
19880  II = II + 1
19890  P(1, I) = XU(II)
19900  P(2, I) = YU(II)
19910 NEXT I
19920 II = 0
19930 FOR I = NUP TO N.TOT
19940  II = II + 1
19950  P(1, I) = XL(II)
19960  P(2, I) = YL(II)
19970 NEXT I
19980 XTRAIL = P(1, N.TOT)
19990 YTRAIL = P(2, N.TOT)
20000 'compute number of points to add
20010 'print fancy titles, date & time, job identification
20020 N = N.TOT
20030 Z = INT(300 / N.TOT + .5)
20040 LOCATE 12, 19: PRINT "Will Insert"; Z; " Points Between Given Points"
20050 LOCATE 13, 28: PRINT "Calculating Cubic Spline"
20060 N(1, 2) = 1
20070 N(1, 3) = .5
20080 FOR J = 1 TO N - 1
20090  L(J) = SQR((P(1, J + 1) - P(1, J)) ^ 2 + (P(2, J + 1) - P(2, J)) ^ 2)
20100 NEXT J
20110 FOR K = 1 TO 2
20120  B(K, 1) = (3 / (2 * L(1))) * (P(K, 2) - P(K, 1))
20130 NEXT K
20140 N(N, 1) = 2
20150 N(N, 2) = 4
20160 FOR K = 1 TO 2
20170  B(K, N) = (6 / L(N - 1)) * (P(K, N) - P(K, N - 1))
20180 NEXT K
20190 FOR K = 1 TO 2
20200  B(K, N) = (6 / L(N - 1)) * (P(K, N) - P(K, N - 1))
20210 NEXT K
20220 FOR J = 2 TO N - 1
20230  N(J, 1) = L(J)
20240  N(J, 2) = 2 * (L(J) + L(J - 1))
20250  N(J, 3) = L(J - 1)
20260  FOR K = 1 TO 2
20270   B(K, J) = 3 * (L(J - 1) ^ 2 * (P(K, J + 1) - P(K, J)) + L(J) ^ 2 * (P(K, J) - P(K, J - 1)))
20280   B(K, J) = B(K, J) / (L(J) * L(J - 1))
20290  NEXT K
20300 NEXT J
20310 FOR I = 2 TO N
20320  IF N(I, 1) = 0 THEN GOTO 20430
20330  D = N(I - 1, 2) / N(I, 1)
20340  FOR K = 1 TO 3
20350   N(I, K) = N(I, K) * D - N(I - 1, K + 1)
20360   B(K, I) = B(K, I) * D - B(K, I - 1)
20370  NEXT K
20380  Q = N(I, 2)
20390  FOR K = 1 TO 3
20400   N(I, K) = N(I, K) / Q
20410   B(K, I) = B(K, I) / Q
20420  NEXT K
20430 NEXT I
20440 FOR K = 1 TO 2
20450  FOR J = 0 TO N - 1
20460   U(K, N - J) = (B(K, N - J) - N(N - J, 3) * U(K, N + 1 - J)) / N(N - J, 2)
20470  NEXT J
20480 NEXT K
20490 IF FLAG THEN RETURN
20500 CLS : LOCATE 12, 28: PRINT "Performing Interpolation"
20510 I = 1
20520 FOR J = 1 TO N - 1
20530  FOR K = 1 TO 2
20540   F(1, K) = P(K, J)
20550   F(2, K) = U(K, J)
20560   F(3, K) = (3 / L(J) ^ 2) * (P(K, J + 1) - P(K, J)) - (1 / L(J)) * (U(K, J + 1) + 2 * U(K, J))
20570   F(4, K) = (-2 / L(J) ^ 3) * (P(K, J + 1) - P(K, J)) + (1 / L(J) ^ 2) * (U(K, J + 1) + U(K, J))
20580  NEXT K
20590  TINC = L(J) / Z
20600  T = 0
20610  IF J = 1 OR T <> 0 THEN
20620    FOR K = 1 TO 2
20630     R(K) = F(1, K) + F(2, K) * T + F(3, K) * T * T + F(4, K) * T ^ 3
20640    NEXT K
20650    C(1, I) = R(1)
20660    C(2, I) = R(2)
20670    I = I + 1
20680  END IF
20690  IF T = 0 AND ABS(C(1, I - 1) - P(1, J)) > .00001 THEN
20700    C(1, I) = P(1, J)
20710    C(2, I) = P(2, J)
20720    I = I + 1
20730  END IF
20740  T = T + TINC
20750  IF T <= L(J) THEN GOTO 20610
20760 NEXT J
20770 NUM.PTS = I - 1
20780 CLS : LOCATE 12, 27: PRINT "Number Of Points = "; NUM.PTS
20790 'move points to P array, go back to find first derivatives for all points
20800 IF ABS(C(1, NUM.PTS) - XTRAIL) > .00001 THEN
20810   NUM.PTS = NUM.PTS + 1
20820   C(1, NUM.PTS) = XTRAIL
20830   C(2, NUM.PTS) = YTRAIL
20840 END IF
20850 FOR I = 1 TO NUM.PTS
20860  P(1, I) = C(1, I)
20870  P(2, I) = C(2, I)
20880 NEXT I
20890 N = NUM.PTS
20900 FLAG = -1
20910 PRINT TAB(26); "Computing First Derivatives"
20920 GOSUB 20060
20930 'convert to slopes in radians, store as B array
20940 FOR I = 1 TO NUM.PTS
20950  B(1, I) = U(2, I) / U(1, I)
20960 NEXT I
20970 'split airfoil into upper and lower surfaces
20980 MIN = P(1, 1)
20990 FOR I = 1 TO NUM.PTS
21000  IF P(1, I) < MIN THEN
21010    IMIN = I
21020    MIN = P(1, I)
21030  END IF
21040 NEXT I
21050 GOSUB 22020
21060 'move upper surface to N arrays 1&2, move slopes to C1 array
21070 FOR I = IMIN TO 1 STEP -1
21080  N(I, 1) = P(1, I)
21090  N(I, 2) = P(2, I)
21100  C(1, I) = B(1, I)
21110 NEXT I
21120 'move lower surface to N arrays 3&4, move slopes to C2 array
21130 K = 0
21140 FOR I = IMIN TO N
21150  K = K + 1
21160  N(K, 3) = P(1, I)
21170  N(K, 4) = P(2, I)
21180  C(2, K) = B(1, I)
21190 NEXT I
21200 NUM.UP.PTS = IMIN
21210 NUM.LO.PTS = NUM.PTS - IMIN + 1
21220 PRINT
21230 CLS : LOCATE 11, 1
21240 PRINT TAB(26); "Number Of Upper Points ="; NUM.UP.PTS
21250 PRINT TAB(26); "Number Of Lower Points ="; NUM.LO.PTS
21260 FOR I = 1 TO NUM.UP.PTS
21270  K = NUM.UP.PTS - I + 1
21280  X(I) = N(K, 1)
21290  Y(I) = N(K, 2)
21300  YP(I) = C(1, K)
21310 NEXT I
21320 PRINT TAB(24); "Calculating Desired Coordinates"
21330 IF SGN(YP(1)) <> SGN(YP(2)) THEN YP(1) = -YP(1)
21340 N% = NUM.UP.PTS
21350 XU(1) = XDD(1)
21360 YU(1) = 0
21370 FOR I = 2 TO NPP
21380   XU(I) = XDD(I)
21390   XD = XDD(I)
21400  GOSUB 21630
21410  YU(I) = YD
21420 NEXT
21430 FOR I = 1 TO NUM.LO.PTS
21440  X(I) = N(I, 3)
21450  Y(I) = N(I, 4)
21460  YP(I) = C(2, I)
21470 NEXT I
21480 IF SGN(YP(1)) <> SGN(YP(2)) THEN YP(1) = -YP(1)
21490 RESTORE
21500 N% = NUM.LO.PTS
21510 XL(1) = XDD(1)
21520 YL(1) = 0
21530 FOR I = 2 TO NPP
21540   XL(I) = XDD(I)
21550   XD = XDD(I)
21560   GOSUB 21630
21570   YL(I) = YD
21580 NEXT
21590 CLS
21600 RETURN
21610 '-----------------------------------------------------------------
21620 'INTERPOLATION ROUTINE FOR REDISTRIBUTION
21630 LL% = 0
21640 UL% = N% + 1
21650 IF (UL% - LL% = 1) THEN GOTO 21800
21660 IN% = (UL% - LL%) \ 2 + LL%
21670 DX = XD - X(IN%)
21680 IF DX < 0 THEN
21690   UL% = IN%
21700   GOTO 21650
21710 END IF
21720 IF DX = 0 THEN
21730   YD = Y(IN%)
21740   GOTO 21990
21750 END IF
21760 IF DX > 0 THEN
21770   LL% = IN%
21780   GOTO 21650
21790 END IF
21800 IF UL% < 2 THEN UL% = 2
21810 IF UL% > N% THEN UL% = N%
21820 LL% = UL% - 1
21830 SL = (Y(UL%) - Y(LL%)) / (X(UL%) - X(LL%))
21840 YZ = Y(LL%) + SL * (XD - X(LL%))
21850 FOR I% = 1 TO 2
21860 K2% = LL% + I% - 1
21870 DY(I%) = Y(K2%) + YP(K2%) * (XD - X(K2%)) - YZ
21880 NEXT I%
21890 TEST = DY(1) * DY(2)
21900 IF TEST = 0 THEN
21910   YD = YZ
21920   GOTO 21990
21930 END IF
21940 IF TEST > 0 THEN
21950   YD = YZ + TEST / (DY(1) + DY(2))
21960   GOTO 21990
21970 END IF
21980 YD = YZ + (TEST * (XD + XD - X(LL%) - X(UL%))) / ((DY(1) - DY(2)) * (X(UL%) - X(LL%)))
21990 RETURN
22000 '-----------------------------------------------------------------
22010 'ROTATION SUBROUTINE
22020 CLS : LOCATE 12, 30: PRINT "Performing Rotation"
22030 XRP = 1
22040 YRP = 0
22050 YTE = (P(2, 1) + P(2, NUM.PTS)) / 2
22060 DELTA = ATN(((P(2, IMIN) - YTE) - YRP) / (P(1, IMIN) - XRP))
22070 DDELTA = DELTA * 180 / 3.1415927#
22080 LOCATE 13, 28: PRINT USING "Angle = ##.######## Deg"; DDELTA
22090 CDLT = COS(DELTA)
22100 SDLT = SIN(DELTA)
22110 TDLT = TAN(DELTA)
22120 FOR I = 1 TO NUM.PTS
22130  XRS = P(1, I)
22140  P(2, I) = P(2, I) - YTE
22150  P(1, I) = XRP + (P(1, I) - XRP) * CDLT + (P(2, I) - YRP) * SDLT
22160  P(2, I) = YRP + (P(2, I) - YRP) * CDLT - (XRS - XRP) * SDLT
22170  B(1, I) = B(1, I) - TDLT
22180 NEXT I
22190 P(1, NUM.PTS) = (P(1, 1) + P(1, NUM.PTS)) / 2
22200 P(1, 1) = P(1, NUM.PTS)
22210 XRLE = P(1, IMIN)
22220 FOR I = 1 TO NUM.PTS
22230  P(1, I) = P(1, I) - XRLE
22240 NEXT I
22250 XRCRD = P(1, 1)
22260 FOR I = 1 TO NUM.PTS
22270  P(1, I) = P(1, I) / XRCRD
22280  P(2, I) = P(2, I) / XRCRD
22290 NEXT I
22300 RETURN
22310 '-----------------------------------------------------------------
22320 'NACA CONSTANTS
22330 N = 1 + N6S * N6T + 2 * N6C
22340 FOR I = 1 TO N
22350  FOR J = 1 TO N6P
22360   READ TF6(I, J)
22370  NEXT J
22380 NEXT I
22390 FOR I = 1 TO N6P
22400  X6(I) = TF6(1, I) / 100
22410 NEXT I
22420 FOR I = 2 TO 1 + N6S * N6T
22430  FOR J = 1 TO N6P
22440   YT6(I - 1, J) = TF6(I, J) / 100
22450  NEXT J
22460 NEXT I
22470 TFI = 1 + N6S * N6T
22480 FOR I = 1 TO N6C
22490  TFI = TFI + 1
22500  FOR J = 1 TO N6P
22510   YC6(I, J) = TF6(TFI, J) / 100
22520  NEXT J
22530  TFI = TFI + 1
22540  FOR J = 1 TO N6P
22550   DYC6(I, J) = TF6(TFI, J)
22560  NEXT J
22570 NEXT I
22580 RETURN
22590 'NACA X/C LOCATIONS
22600 DATA 0.0000, 0.5000, 0.7500, 1.2500, 2.5000, 5.0000, 7.5000, 10.0000, 15.0000, 20.0000, 25.0000, 30.0000, 35.0000, 40.0000, 45.0000, 50.0000, 55.0000, 60.0000, 65.0000, 70.0000, 75.0000, 80.0000, 85.0000, 90.0000, 95.0000, 100.0000
22610 'NACA THICKNESS FORMS 6-21%, 63-67 SERIES
22620 DATA 0.0000, 0.5030, 0.6090, 0.7710, 1.0570, 1.4620, 1.7660, 2.0100, 2.3860, 2.6560, 2.8410, 2.9540, 3.0000, 2.9710, 2.8770, 2.7230, 2.5170, 2.2670, 1.9820, 1.6700, 1.3420, 1.0080, 0.6830, 0.3830, 0.1380, 0.0000
22630 DATA 0.0000, 0.6683, 0.8089, 1.0244, 1.4073, 1.9517, 2.3588, 2.6858, 3.1896, 3.5501, 3.7965, 3.9459, 3.9999, 3.9579, 3.8254, 3.6133, 3.3338, 2.9984, 2.6173, 2.2003, 1.7616, 1.3171, 0.8870, 0.4963, 0.1777, 0.0000
22640 DATA 0.0000, 0.7490, 0.9060, 1.1510, 1.5820, 2.1960, 2.6550, 3.0240, 3.5910, 3.9970, 4.2750, 4.4420, 4.5000, 4.4470, 4.2960, 4.0560, 3.7390, 3.3580, 2.9280, 2.4580, 1.9660, 1.4710, 0.9900, 0.5500, 0.1960, 0.0000
22650 DATA 0.0000, 0.8290, 1.0040, 1.2750, 1.7560, 2.4400, 2.9500, 3.3620, 3.9940, 4.4450, 4.7530, 4.9380, 5.0000, 4.9380, 4.7660, 4.4960, 4.1400, 3.7150, 3.2340, 2.7120, 2.1660, 1.6180, 1.0880, 0.6040, 0.2140, 0.0000
22660 DATA 0.0000, 0.9850, 1.1940, 1.5190, 2.1020, 2.9250, 3.5420, 4.0390, 4.7990, 5.3420, 5.7120, 5.9300, 6.0000, 5.9200, 5.7040, 5.3700, 4.9350, 4.4200, 3.8400, 3.2100, 2.5560, 1.9020, 1.2740, 0.7070, 0.2500, 0.0000
22670 DATA 0.0000, 1.2040, 1.4620, 1.8780, 2.6100, 3.6480, 4.4270, 5.0550, 6.0110, 6.6930, 7.1550, 7.4210, 7.5000, 7.3860, 7.0990, 6.6650, 6.1080, 5.4530, 4.7210, 3.9340, 3.1190, 2.3100, 1.5410, 0.8520, 0.3000, 0.0000
22680 DATA 0.0000, 1.4040, 1.7130, 2.2170, 3.1040, 4.3620, 5.3080, 6.0680, 7.2250, 8.0480, 8.6000, 8.9130, 9.0000, 8.8450, 8.4820, 7.9420, 7.2560, 6.4550, 5.5670, 4.6220, 3.6500, 2.6910, 1.7870, 0.9850, 0.3480, 0.0000
22690 DATA 0.0000, 1.5830, 1.9370, 2.5270, 3.5770, 5.0650, 6.1820, 7.0800, 8.4410, 9.4100, 10.0530, 10.4120, 10.5000, 10.2980, 9.8540, 9.2060, 8.3900, 7.4410, 6.3960, 5.2900, 4.1600, 3.0540, 2.0210, 1.1130, 0.3920, 0.0000
22700 DATA 0.0000, 0.4940, 0.5960, 0.7540, 1.0240, 1.4050, 1.6920, 1.9280, 2.2980, 2.5720, 2.7720, 2.9070, 2.9810, 2.9950, 2.9190, 2.7750, 2.5750, 2.3310, 2.0500, 1.7400, 1.4120, 1.0720, 0.7370, 0.4230, 0.1570, 0.0000
22710 DATA 0.0000, 0.6580, 0.7940, 1.0050, 1.3650, 1.8750, 2.2590, 2.5740, 3.0690, 3.4370, 3.7040, 3.8840, 3.9790, 3.9920, 3.8830, 3.6840, 3.4110, 3.0810, 2.7040, 2.2910, 1.8540, 1.4040, 0.9610, 0.5500, 0.2060, 0.0000
22720 DATA 0.0000, 0.7390, 0.8920, 1.1280, 1.5330, 2.1090, 2.5430, 2.8980, 3.4550, 3.8680, 4.1700, 4.3730, 4.4790, 4.4900, 4.3640, 4.1360, 3.8260, 3.4520, 3.0260, 2.5610, 2.0690, 1.5640, 1.0690, 0.6110, 0.2270, 0.0000
22730 DATA 0.0000, 0.8200, 0.9890, 1.2500, 1.7010, 2.3430, 2.8260, 3.2210, 3.8420, 4.3020, 4.6390, 4.8640, 4.9800, 4.9880, 4.8430, 4.5860, 4.2380, 3.8200, 3.3450, 2.8270, 2.2810, 1.7220, 1.1760, 0.6710, 0.2480, 0.0000
22740 DATA 0.0000, 0.9780, 1.1790, 1.4900, 2.0350, 2.8100, 3.3940, 3.8710, 4.6200, 5.1730, 5.5760, 5.8440, 5.9780, 5.9810, 5.7980, 5.4800, 5.0560, 4.5480, 3.9740, 3.3500, 2.6950, 2.0290, 1.3820, 0.7860, 0.2880, 0.0000
22750 DATA 0.0000, 1.2080, 1.4560, 1.8420, 2.5280, 3.5040, 4.2400, 4.8420, 5.7850, 6.4800, 6.9850, 7.3190, 7.4820, 7.4730, 7.2240, 6.8100, 6.2660, 5.6200, 4.8950, 4.1130, 3.2960, 2.4720, 1.6770, 0.9500, 0.3460, 0.0000
22760 DATA 0.0000, 1.4280, 1.7200, 2.1770, 3.0050, 4.1860, 5.0760, 5.8030, 6.9420, 7.7820, 8.3910, 8.7890, 8.9790, 8.9520, 8.6300, 8.1140, 7.4450, 6.6580, 5.7820, 4.8420, 3.8660, 2.8880, 1.9510, 1.1010, 0.4000, 0.0000
22770 DATA 0.0000, 1.6460, 1.9850, 2.5170, 3.4850, 4.8710, 5.9150, 6.7690, 8.1080, 9.0950, 9.8070, 10.2690, 10.4810, 10.4310, 10.0300, 9.4040, 8.6070, 7.6780, 6.6490, 5.5490, 4.4160, 3.2870, 2.2130, 1.2450, 0.4490, 0.0000
22780 DATA 0.0000, 0.4760, 0.5740, 0.7170, 0.9560, 1.3100, 1.5890, 1.8240, 2.1970, 2.4820, 2.6970, 2.8520, 2.9520, 2.9980, 2.9830, 2.9000, 2.7410, 2.5180, 2.2460, 1.9350, 1.5940, 1.2330, 0.8650, 0.5100, 0.1950, 0.0000
22790 DATA 0.0000, 0.6270, 0.7560, 0.9450, 1.2670, 1.7450, 2.1180, 2.4320, 2.9310, 3.3120, 3.5990, 3.8050, 3.9380, 3.9980, 3.9740, 3.8570, 3.6380, 3.3370, 2.9710, 2.5530, 2.0960, 1.6170, 1.1310, 0.6640, 0.2520, 0.0000
22800 DATA 0.0000, 0.7000, 0.8450, 1.0580, 1.4210, 1.9610, 2.3830, 2.7360, 3.2990, 3.7270, 4.0500, 4.2820, 4.4310, 4.4960, 4.4690, 4.3360, 4.0860, 3.7430, 3.3280, 2.8560, 2.3420, 1.8050, 1.2600, 0.7380, 0.2800, 0.0000
22810 DATA 0.0000, 0.7720, 0.9320, 1.1690, 1.5740, 2.1770, 2.6470, 3.0400, 3.6660, 4.1430, 4.5030, 4.7600, 4.9240, 4.9960, 4.9630, 4.8120, 4.5300, 4.1460, 3.6820, 3.1560, 2.5840, 1.9870, 1.3850, 0.8100, 0.3060, 0.0000
22820 DATA 0.0000, 0.9230, 1.1090, 1.3870, 1.8750, 2.6060, 3.1720, 3.6470, 4.4020, 4.9750, 5.4060, 5.7160, 5.9120, 5.9970, 5.9490, 5.7570, 5.4120, 4.9430, 4.3810, 3.7430, 3.0590, 2.3450, 1.6300, 0.9470, 0.3560, 0.0000
22830 DATA 0.0000, 1.1240, 1.3560, 1.7020, 2.3240, 3.2450, 3.9590, 4.5550, 5.5040, 6.2230, 6.7640, 7.1520, 7.3960, 7.4980, 7.4270, 7.1680, 6.7200, 6.1180, 5.4030, 4.6000, 3.7440, 2.8580, 1.9770, 1.1440, 0.4280, 0.0000
22840 DATA 0.0000, 1.3370, 1.6080, 2.0140, 2.7510, 3.8660, 4.7330, 5.4570, 6.6060, 7.4760, 8.1290, 8.5950, 8.8860, 8.9990, 8.9010, 8.5680, 8.0080, 7.2670, 6.3950, 5.4260, 4.3960, 3.3380, 2.2950, 1.3190, 0.4900, 0.0000
22850 DATA 0.0000, 1.5220, 1.8380, 2.3010, 3.1540, 4.4720, 5.4980, 6.3520, 7.7000, 8.7200, 9.4870, 10.0360, 10.3750, 10.4990, 10.3660, 9.9520, 9.2770, 8.3900, 7.3600, 6.2240, 5.0240, 3.8000, 2.5980, 1.4840, 0.5460, 0.0000
22860 DATA 0.0000, 0.4610, 0.5540, 0.6930, 0.9180, 1.2570, 1.5240, 1.7520, 2.1190, 2.4010, 2.6180, 2.7820, 2.8990, 2.9710, 3.0000, 2.9850, 2.9250, 2.8150, 2.6110, 2.3160, 1.9530, 1.5430, 1.1070, 0.6650, 0.2620, 0.0000
22870 DATA 0.0000, 0.6100, 0.7350, 0.9190, 1.2190, 1.6730, 2.0310, 2.3350, 2.8260, 3.2010, 3.4900, 3.7090, 3.8650, 3.9620, 4.0000, 3.9780, 3.8960, 3.7400, 3.4590, 3.0620, 2.5740, 2.0270, 1.4470, 0.8640, 0.3380, 0.0000
22880 DATA 0.0000, 0.6870, 0.8240, 1.0300, 1.3680, 1.8800, 2.2830, 2.6260, 3.1780, 3.6010, 3.9270, 4.1730, 4.3480, 4.4570, 4.4990, 4.4750, 4.3810, 4.2040, 3.8820, 3.4280, 2.8770, 2.2630, 1.6110, 0.9610, 0.3740, 0.0000
22890 DATA 0.0000, 0.7590, 0.9130, 1.1410, 1.5160, 2.0870, 2.5360, 2.9170, 3.5300, 4.0010, 4.3630, 4.6360, 4.8320, 4.9530, 5.0000, 4.9710, 4.8650, 4.6650, 4.3020, 3.7870, 3.1760, 2.4940, 1.7730, 1.0540, 0.4080, 0.0000
22900 DATA 0.0000, 0.9060, 1.0870, 1.3580, 1.8080, 2.4960, 3.0370, 3.4960, 4.2340, 4.8010, 5.2380, 5.5680, 5.8030, 5.9470, 6.0000, 5.9650, 5.8360, 5.5880, 5.1390, 4.5150, 3.7670, 2.9440, 2.0830, 1.2340, 0.4740, 0.0000
22910 DATA 0.0000, 1.1220, 1.3430, 1.6750, 2.2350, 3.1000, 3.7810, 4.3580, 5.2860, 5.9950, 6.5430, 6.9560, 7.2500, 7.4300, 7.4950, 7.4500, 7.2830, 6.9590, 6.3720, 5.5760, 4.6320, 3.5980, 2.5300, 1.4890, 0.5660, 0.0000
22920 DATA 0.0000, 1.3230, 1.5710, 1.9520, 2.6460, 3.6900, 4.5130, 5.2100, 6.3330, 7.1880, 7.8480, 8.3460, 8.7010, 8.9180, 8.9980, 8.9420, 8.7330, 8.3230, 7.5800, 6.5970, 5.4510, 4.2060, 2.9340, 1.7140, 0.6460, 0.0000
22930 DATA 0.0000, 1.5250, 1.8040, 2.2400, 3.0450, 4.2690, 5.2330, 6.0520, 7.3690, 8.3760, 9.1530, 9.7380, 10.1540, 10.4070, 10.5000, 10.4340, 10.1860, 9.6920, 8.7930, 7.6100, 6.2510, 4.7960, 3.3240, 1.9240, 0.7170, 0.0000
22940 DATA 0.0000, 0.4879, 0.5847, 0.7351, 0.9888, 1.3236, 1.5808, 1.7909, 2.1262, 2.3856, 2.5886, 2.7464, 2.8649, 2.9451, 2.9902, 2.9985, 2.9686, 2.8954, 2.7713, 2.5854, 2.3063, 1.9056, 1.4274, 0.8994, 0.3585, 0.0000
22950 DATA 0.0000, 0.6462, 0.7740, 0.9739, 1.3118, 1.7595, 2.1036, 2.3844, 2.8325, 3.1791, 3.4503, 3.6611, 3.8196, 3.9267, 3.9869, 3.9980, 3.9579, 3.8592, 3.6917, 3.4395, 3.0581, 2.5134, 1.8738, 1.1706, 0.4637, 0.0000
22960 DATA 0.0000, 0.7240, 0.8669, 1.0915, 1.4715, 1.9759, 2.3637, 2.6800, 3.1848, 3.5752, 3.8807, 4.1182, 4.2967, 4.4174, 4.4853, 4.4977, 4.4525, 4.3410, 4.1514, 3.8652, 3.4308, 2.8123, 2.0917, 1.3013, 0.5139, 0.0000
22970 DATA 0.0000, 0.8008, 0.9586, 1.2078, 1.6297, 2.1911, 2.6229, 2.9748, 3.5365, 3.9709, 4.3108, 4.5750, 4.7737, 4.9081, 4.9837, 4.9975, 4.9471, 4.8227, 4.6109, 4.2901, 3.8013, 3.1081, 2.3061, 1.4288, 0.5628, 0.0000
22980 DATA 0.0000, 0.9507, 1.1374, 1.4354, 1.9412, 2.6173, 3.1376, 3.5615, 4.2377, 4.7605, 5.1696, 5.4876, 5.7269, 5.8889, 5.9802, 5.9971, 5.9365, 5.7860, 5.5290, 5.1373, 4.5356, 3.6902, 2.7246, 1.6744, 0.6565, 0.0000
22990 DATA 0.0000, 1.1651, 1.3925, 1.7625, 2.3940, 3.2443, 3.8994, 4.4325, 5.2826, 5.9396, 6.4537, 6.8536, 7.1547, 7.3588, 7.4747, 7.4966, 7.4209, 7.2307, 6.9048, 6.4020, 5.6210, 4.5399, 3.3270, 2.0209, 0.7880, 0.0000
23000 DATA 0.0000, 1.3646, 1.6289, 2.0689, 2.8258, 3.8530, 4.6458, 5.2904, 6.3174, 7.1106, 7.7315, 8.2148, 8.5791, 8.8268, 8.9683, 8.9963, 8.9060, 8.6761, 8.2797, 7.6591, 6.6858, 5.3628, 3.8991, 2.3435, 0.9097, 0.0000
23010 DATA 0.0000, 1.5470, 1.8438, 2.3506, 3.2327, 4.4397, 5.3736, 6.1323, 7.3397, 8.2720, 9.0017, 9.5703, 9.9994, 10.2922, 10.4609, 10.4966, 10.3924, 10.1232, 9.6546, 8.9093, 7.7296, 6.1596, 4.4414, 2.6443, 1.0229, 0.0000
23020 'NACA CAMBER LINES AND DERIVATIVES
23030 DATA 0.00000, 0.46000, 0.64100, 0.96400, 1.64100, 2.69300, 3.50700, 4.16100, 5.12400, 5.74700, 6.11400, 6.27700, 6.27300, 6.13000, 5.87100, 5.51600, 5.08100, 4.58100, 4.03200, 3.44500, 2.83600, 2.21700, 1.60400, 1.01300, 0.46700, 0.00000
23040 DATA 0.00000, 0.75867, 0.69212, 0.60715, 0.48892, 0.36561, 0.29028, 0.23515, 0.15508, 0.09693, 0.05156, 0.01482, -.01554, -.04086, -.06201, -.07958, -.09395, -.10539, -.11406, -.12003, -.12329, -.12371, -.12099, -.11455, -.10301, -.07958
23050 DATA 0.00000, 0.44000, 0.61600, 0.93300, 1.60800, 2.68900, 3.55100, 4.25300, 5.26100, 5.90500, 6.28200, 6.44900, 6.44300, 6.29600, 6.02900, 5.66400, 5.21800, 4.70600, 4.14200, 3.54100, 2.91600, 2.28100, 1.65200, 1.04500, 0.48200, 0.00000
23060 DATA 0.00000, 0.73441, 0.67479, 0.59896, 0.49366, 0.38235, 0.31067, 0.25057, 0.16087, 0.09981, 0.05281, 0.01498, -.01617, -.04210, -.06373, -.08168, -.09637, -.10806, -.11694, -.12307, -.12644, -.12693, -.12425, -.11781, -.10620, -.08258
23070 DATA 0.00000, 0.41400, 0.58100, 0.88200, 1.53000, 2.58300, 3.44300, 4.16900, 5.31700, 6.11700, 6.57200, 6.77700, 6.78900, 6.64600, 6.37300, 5.99400, 5.52700, 4.98900, 4.39600, 3.76200, 3.10200, 2.43100, 1.76400, 1.11900, 0.51800, 0.00000
23080 DATA 0.00000, 0.69492, 0.64047, 0.57135, 0.47592, 0.37661, 0.31487, 0.26803, 0.19373, 0.12405, 0.06345, 0.02030, -.01418, -.04246, -.06588, -.08522, -.10101, -.11359, -.12317, -.12985, -.13363, -.13440, -.13186, -.12541, -.11361, -.08941
23090 DATA 0.00000, 0.38900, 0.54600, 0.83200, 1.44800, 2.45800, 3.29300, 4.00800, 5.17200, 6.05200, 6.68500, 7.07200, 7.17500, 7.07400, 6.81600, 6.43300, 5.94900, 5.38300, 4.75300, 4.07600, 3.36800, 2.64500, 1.92400, 1.22400, 0.57000, 0.00000
23100 DATA 0.00000, 0.65536, 0.60524, 0.54158, 0.45399, 0.36344, 0.30780, 0.26621, 0.20246, 0.15068, 0.10278, 0.04833, -.00205, -.03710, -.06492, -.08746, -.10567, -.12014, -.13119, -.13901, -.14365, -.14500, -.14279, -.13638, -.12430, -.09907
23110 DATA 0.00000, 0.36600, 0.51400, 0.78400, 1.36700, 2.33000, 3.13100, 3.82400, 4.96800, 5.86200, 6.54600, 7.03900, 7.34300, 7.43900, 7.27500, 6.92900, 6.44900, 5.86400, 5.19900, 4.47500, 3.70900, 2.92200, 2.13200, 1.36100, 0.63600, 0.00000
23120 DATA 0.00000, 0.61759, 0.57105, 0.51210, 0.43106, 0.34764, 0.29671, 0.25892, 0.20185, 0.15682, 0.11733, 0.07988, 0.04136, -.00721, -.05321, -.08380, -.10734, -.12567, -.13962, -.14963, -.15589, -.15837, -.15683, -.15062, -.13816, -.11138
23130 DATA 0.00000, 0.34500, 0.48500, 0.73500, 1.29500, 2.20500, 2.97000, 3.63000, 4.74000, 5.62000, 6.31000, 6.84000, 7.21500, 7.43000, 7.49000, 7.35000, 6.96500, 6.40500, 5.72500, 4.95500, 4.13000, 3.26500, 2.39500, 1.53500, 0.72000, 0.00000
23140 DATA 0.00000, 0.58195, 0.53855, 0.48360, 0.40815, 0.33070, 0.28365, 0.24890, 0.19690, 0.15650, 0.12180, 0.09000, 0.05930, 0.02800, -.00630, -.05305, -.09765, -.12550, -.14570, -.16015, -.16960, -.17435, -.17415, -.16850, -.15565, -.12660
23150 DATA 0.00000, 0.32500, 0.45500, 0.69500, 1.22000, 2.08000, 2.80500, 3.43500, 4.49500, 5.34500, 6.03500, 6.57000, 6.96500, 7.23500, 7.37000, 7.37000, 7.22000, 6.88000, 6.27500, 5.50500, 4.63000, 3.69500, 2.72000, 1.75500, 0.82500, 0.00000
23160 DATA 0.00000, 0.54825, 0.50760, 0.45615, 0.38555, 0.31325, 0.26950, 0.23730, 0.18935, 0.15210, 0.12125, 0.09310, 0.06660, 0.04060, 0.01405, -.01435, -.04700, -.09470, -.14015, -.16595, -.18270, -.19225, -.19515, -.19095, -.17790, -.14550
23170 DATA 0.00000, 0.30500, 0.42500, 0.65500, 1.16000, 1.95500, 2.64500, 3.24000, 4.24500, 5.06000, 5.71500, 6.24000, 6.63500, 6.92500, 7.09500, 7.15500, 7.09000, 6.90000, 6.56500, 6.03000, 5.20500, 4.21500, 3.14000, 2.03500, 0.96500, 0.00000
23180 DATA 0.00000, 0.51620, 0.47795, 0.42960, 0.36325, 0.29545, 0.25450, 0.22445, 0.17995, 0.14595, 0.11740, 0.09200, 0.06840, 0.04570, 0.02315, 0.00000, -.02455, -.05185, -.08475, -.13650, -.18510, -.20855, -.21955, -.21960, -.20725, -.16985
23190 DATA 0.00000, 0.28700, 0.40400, 0.61600, 1.07700, 1.84100, 2.48300, 3.04300, 3.98500, 4.74800, 5.36700, 5.86300, 6.24800, 6.52800, 6.70900, 6.79000, 6.77000, 6.64400, 6.40500, 6.03700, 5.51400, 4.77100, 3.68300, 2.43500, 1.16300, 0.00000
23200 DATA 0.00000, 0.48535, 0.44925, 0.40359, 0.34104, 0.27718, 0.23868, 0.21050, 0.16892, 0.13734, 0.11101, 0.08775, 0.06634, 0.04601, 0.02613, 0.00620, -.01433, -.03611, -.06010, -.08790, -.12311, -.18412, -.23921, -.25583, -.24904, -.20385
23210 DATA 0.00000, 0.26900, 0.37900, 0.57700, 1.00800, 1.72000, 2.31600, 2.83500, 3.70700, 4.41000, 4.98000, 5.43500, 5.78700, 6.04500, 6.21200, 6.29000, 6.27900, 6.17800, 5.98100, 5.68100, 5.26500, 4.71400, 3.98700, 2.98400, 1.50300, 0.00000
23220 DATA 0.00000, 0.45482, 0.42064, 0.37740, 0.31821, 0.25786, 0.22153, 0.19500, 0.15595, 0.12644, 0.10196, 0.08047, 0.06084, 0.04234, 0.02447, 0.00678, -.01111, -.02965, -.04938, -.07103, -.09583, -.12605, -.16727, -.25204, -.31463, -.26086
23230 DATA 0.00000, 0.25000, 0.35000, 0.53500, 0.93000, 1.58000, 2.12000, 2.58500, 3.36500, 3.98000, 4.47500, 4.86000, 5.15000, 5.35500, 5.47500, 5.51500, 5.47500, 5.35500, 5.15000, 4.86000, 4.47500, 3.98000, 3.36500, 2.58500, 1.58000, 0.00000
23240 DATA 0.00000, 0.42120, 0.38875, 0.34770, 0.29155, 0.23430, 0.19995, 0.17485, 0.13805, 0.11030, 0.08745, 0.06745, 0.04925, 0.03225, 0.01595, 0.00000, -.01595, -.03225, -.04925, -.06745, -.08745, -.11030, -.13805, -.17485, -.23430, 0.00000
23250 '
23260 '-----------------------------------------------------------------
23270 'MONOTONE SPLINE SUBROUTINE
23280 'HOST MUST DIMENSION ENDSLP(2), DYDX(2)
23290 'IF ONLY ONE POINT IS IN SET, RETURN THAT VALUE
23300 IF N = 1 THEN
23310  XD = XI(1)
23320  YD = YI(1)
23330  GOTO 23890
23340 END IF
23350 'FIND POINTS IN SET THAT BOUND XD VIA BINARY SEARCH
23360 LOWLIM = 1: UPPLIM = N
23370 WHILE UPPLIM - LOWLIM > 1
23380  TSTLCN = (UPPLIM + LOWLIM) \ 2
23390  IF XI(TSTLCN) > XD THEN
23400    UPPLIM = TSTLCN
23410   ELSE
23420    LOWLIM = TSTLCN
23430  END IF
23440  IF XI(LOWLIM) = XD THEN YD = YI(LOWLIM): GOTO 23890
23450 WEND
23460 IF LOWLIM = N THEN LOWLIM = N - 1
23470 UPPLIM = LOWLIM + 1
23480 'CALC LINEAR INTERPOLANT
23490 LNSLP = (YI(UPPLIM) - YI(LOWLIM)) / (XI(UPPLIM) - XI(LOWLIM))
23500 YLINR = LNSLP * (XD - XI(LOWLIM)) + YI(LOWLIM)
23510 'CALC END SLOPES
23520 PTBFR = LOWLIM
23530 WHCHEND = 1
23540 FOR SGMN = 1 TO 2
23550  IF PTBFR > 1 AND PTBFR < N - 1 THEN
23560    PNT1 = PTBFR - 1
23570   ELSEIF PTBFR > N - 2 THEN
23580    PNT1 = N - 2
23590   ELSE
23600    PNT1 = 1
23610    WHCHEND = 0
23620  END IF
23630  PNT2 = PNT1 + 1
23640  PNT3 = PNT1 + 2
23650  DELTX1 = XI(PNT1) - XI(PNT2)
23660  DELTX2 = XI(PNT3) - XI(PNT2)
23670  DELTY1 = YI(PNT1) - YI(PNT2)
23680  DELTY2 = YI(PNT3) - YI(PNT2)
23690  KNST = (DELTY1 * DELTX2 / DELTX1 - DELTY2 * DELTX1 / DELTX2) / (DELTX2 - DELTX1)
23700  CFCNT = (DELTY1 / DELTX1 - DELTY2 / DELTX2) / (DELTX1 - DELTX2)
23710  ENDSLP(SGMN) = KNST - 2 * CFCNT * XI(PNT2) + 2 * CFCNT * XI(PTBFR)
23720  PTBFR = PTBFR + 1
23730  IF WHCHEND > 0 THEN PNT1 = PNT1 + 1
23740 NEXT SGMN
23750 'FIND DIST END SLOPES PUT POINT ABOVE LINEAR
23760 FOR SGMN = 1 TO 2
23770  PTBFR = LOWLIM + SGMN - 1
23780  DYDX(SGMN) = ENDSLP(SGMN) * (XD - XI(PTBFR)) + YI(PTBFR) - YLINR
23790 NEXT SGMN
23800 'TEST FOR ABOVE, BELOW OR ON LINEAR AND PRODUCE YD
23810 POSN = DYDX(1) * DYDX(2)
23820 IF POSN = 0 THEN
23830   YD = YLINR
23840  ELSEIF POSN > 0 THEN
23850   YD = YLINR + POSN / (DYDX(1) + DYDX(2))
23860  ELSE
23870   YD = YLINR + (POSN * (2 * XD - XI(LOWLIM) - XI(UPPLIM))) / ((DYDX(1) - DYDX(2)) * (XI(UPPLIM) - XI(LOWLIM)))
23880 END IF
23890 RETURN
23900 '-----------------------------------------------------------------
23910 'SUBROUTINE TO FILL IN LOWER SURFACE COVE
23920 LOCATE 11, 30: PRINT "Filling Lower Cove"
23930 STP = .01
23940 N = NP
23950 FOR I = 1 TO N
23960  XI(I) = XL(I)
23970  YI(I) = YL(I)
23980 NEXT I
23990 XD = PP / 10 + STP * 5
24000 GOSUB 23300
24010 XP = XD
24020 YP = YD
24030 XD = XD + .00001
24040 GOSUB 23300
24050 LSLP = (YD - YP) / (XD - XP)
24060 SLP2TE = (YP - YL(N)) / (XP - XL(N))
24070 SLPDIF = ABS(LSLP - SLP2TE)
24080 IF SLPDIF > STOL AND LSLP < SLP2TE THEN
24090   XD = XP + STP
24100   IF XP < XL(N) THEN GOTO 24000
24110  ELSEIF SLPDIF > STOL AND LSLP > SLP2TE THEN
24120   XD = XP - STP
24130   STP = STP * .5
24140   IF STP > .000001 THEN GOTO 24000
24150 END IF
24160 FOR I = 1 TO N
24170  IF XL(I) > XP THEN
24180    YL(I) = (XL(I) - XP) / (XL(N) - XP) * (YL(N) - YP) + YP
24190  END IF
24200 NEXT I
24210 LOCATE 11, 30: PRINT "                  "
24220 RETURN
24230 '-----------------------------------------------------------------
24240 ' TRANSITION AND ROUGHNESS SUMMARY
24250 PRINT
24260 PRINT " Available Transition Options"
24270 PRINT "  1) Transition At Laminar Separation"
24280 PRINT "  2) Fixed Transition"
24290 PRINT "  3) Natural Transtion"
24300 PRINT "  4) Transition With Roughness Factor"
24310 PRINT
24320 PRINT " Available Roughness Levels Are 0 To 6"
24330 PRINT "  0 = Clean Surface"
24340 PRINT "  2 = Dirty Surface"
24350 PRINT "  4 = Excessively Dirty Surface"
24360 PRINT "  6 = Severely Wavy Surface"
24370 PRINT
24380 RETURN
24390 '-----------------------------------------------------------------
24400 'SET INPUT LIMITS SUBROUTINE
24410 PRINT
24420 PRINT " Input Point Limits Are "
24430 PRINT "  1) Number Of Computational Points ="; NP
24440 PRINT "  2) Number Of Points On The Airfoil ="; MPTS
24450 PRINT
24460 INPUT " Input Limit To Be Changed (0=None)[0]"; LCH
24470 ON LCH + 1 GOTO 24550, 24490, 24520
24480 GOTO 24460
24490 PRINT
24500 INPUT "  Enter New Number Of Computational Points"; NP
24510 GOTO 24450
24520 PRINT
24530 INPUT " Enter New Limit On Number Of Points On The Airfoil ="; MPTS
24540 GOTO 24450
24550 PRINT
24560 PRINT " Possible Monitor Types Are:"
24570 PRINT "  1) CGA (RGB)"
24580 PRINT "  2) EGA With 64K Of Memory"
24590 PRINT "  3) EGA With More Than 64K Of Memory"
24600 PRINT "  4) VGA"
24610 PRINT
24620 INPUT " Your Choice [3 or 4]:"; DTP
24630 IF DTP = 0 THEN DTP = 3
24640 IF DTP = 1 THEN QM = 320: RM = 100: NCLR = 3
24650 IF DTP = 2 THEN NCLR = 3
24660 CLS : PRINT
24670 RETURN

