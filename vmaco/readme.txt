

VARIABLE METRIC ALGORITHM FOR CONSTRAINED OPTIMIZATION      \vmaco\readme.txt

The files for this program are in the directory \vmaco on the CD-ROM
  readme.txt      this file of general description
  original.src    the original copy of the source code (from COSMIC)
  vmaco.f90       the source code for the Vmaco module
  test1.f90       the source code for the first sample case
  test1.out       the output produced by running test1
  test2.f90       the source code for the second sample case
  test2.out       the output produced by running test2
  msc21275.txt    the original program description from COSMIC

The reference documents for this program may be accessed
from the web page http://www.pdas.com/vmacorefs.html. 

To compile this program for your machine, use the command
   gfortran  vmaco.f90 -o vmaco.exe
Linux and Macintosh users may prefer to omit the .exe on the file name.

To use this program, make a directory on your hard disk and
copy these files to that directory. Vmaco is not a program, but a
module of subroutines that may be called by a driver program
to search for the minimum. Two test programs are shown to illustrate
the procedures for using Vmaco. To run test1, first compile the
program with the command
   gfortran test1.f90
This will produce a file called a.exe on Windows or a.out on Linux
or Unix or Macintosh. Execute this program to see the results of test1.
