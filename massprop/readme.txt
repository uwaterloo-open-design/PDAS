

MASSPROP- MASS PROPERTIES OF A RIGID STRUCTURE              \massprop\readme.txt

The files for this program are in the directory \massprop on the CD-ROM
  readme.txt      general description
  input.txt       guide for writing the input file
  massprop.f90    the complete source code in modern Fortran
  massprop.exe    the executable for Windows
  massprop.lnx    the executable for Linux
  massprop.mac    the executable for Macintosh OS X (Intel)
  lar12454.txt    the original program description from COSMIC
  original.src    the original copy of the source code (from COSMIC)

Sample cases for this program are:
  example1.inp    input data for example 1
  example1.out    output data for example1.inp
  example2.inp    input data for example 2
  example2.out    output data for example2.inp
  example3.inp    input data for example 3
  example3.out    output data for example3.inp

The reference documents for this program may be accessed
from the web page http://www.pdas.com/massproprefs.html. 

If you need to recompile this program for your machine, use the command
   gfortran  massprop.f90 -o massprop.exe
Linux and Macintosh users may prefer to omit the .exe on the file name.
 
This program first asks for the name of the input file. This must
be a file written to conform to input.txt.  After calculating 
the solution, the program produces a file called massprop.out that contains a 
wealth of information concerning the flow problem to be solved.

