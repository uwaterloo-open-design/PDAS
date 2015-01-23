

VU-CALC - AN INTERACTIVE CALCULATOR FOR COMPRESSIBLE FLOW     \vucalc\readme.txt

The files for this program are in the directory \vucalc on the CD-ROM  and in the
archive file vucalc.zip that may be downloaded from the PDAS web site.

  readme.txt      general description
  main.dfm        the main form
  calcsubs.pas    the code for unit calcsubs
  getMach.pas     the code for unit GetMach
  loadMemo.pas    the code for unit LoadMemo
  main.pas        the code for unit main 
  naca1135.pas    the code for unit naca1135
  rayfanno.pas    the code for unit RayFanno
  calc.c          the original code from NASA VuCalc (for SGI machine)

  vucalc.dpr      the Delphi project file
  vucalc.exe      the executable file for Windows

To use this program, create a directory on your hard disk and copy the 
file vucalc.exe to that directory. The other files are only needed if you 
want to modify the program.

VU-CALC is a calculator that allows you to solve various problems in
compressible flow. In particular, one can solve for isentropic flow,
normal shock, oblique shock, and for characteristics of the standard
atmosphere. You may solve for those quantities that are direct lookup 
and also for those that are a reverse lookup. 

VU-CALC was designed and written by Tom Benson of NASA Lewis Research
Center for SGI workstations (or for any machine that supports the 
library and X-windows). The original code was written in C and used 
the Forms X-windows library. I converted it to run under Microsoft 
Windows using the Delphi programming environment.

At the time it was written, VuCalc was a highly innovative and original
program. Since that time, writing flow calculators has become a very 
popular activity. If you do a search for Compressible Flow, you will 
find a large collection of such calculators. Some, like VuCalc, are 
application programs that run without any web connection. Others are 
web based and use HTML forms to access the user information. The results
are calculated on the server and a new page is presented to the client. 
Still others are java or javascript applets that execute on the user's
machine with no further involvement of the web server.
 

References:

Benson, Thomas J., and Higgs, C. Fred III: X Based Interactive
Computer Graphics Applications for Aerodynamic Design and Education.
Presented at AIAA Aerospace Sciences Meeting, Reno 1996.
