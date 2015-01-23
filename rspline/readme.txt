

RATIONAL SPLINE SUBROUTINES                                  \rspline\readme.txt

The files for this program are in the directory \rspline on the CD-ROM
  readme.txt      this file of general description
  ratspl.f90      the complete source code in modern Fortran
  example1.f90    example program #1
  example2.f90    example program #2
  example1.out    output from running example1
  example2.out    output from running example2
  lar13694.txt    the original program description from COSMIC
  original.src    the original copy of the source code (from COSMIC)

The reference documents for this program may be accessed
from the web page http://www.pdas.com/rsplinerefs.html. 

This is not a program, but a collection of subroutines that may be called by 
a driver program. Two test programs are shown to illustrate the procedures.

To run example1, first compile the sample program with the command
   gfortran example1.f90

This will produce an executable file that you can run to perform the calculations.
On Windows, the file will be called a.exe and on Linux or Macintosh, it will be
called a.out.

