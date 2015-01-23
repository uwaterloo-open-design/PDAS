

GASP - THERMODYNAMIC AND TRANSPORT PROPERTIES                   \gasp\readme.txt
  OF HELIUM, METHANE, NEON, NITROGEN, CARBON MONOXIDE, 
  CARBON DIOXIDE, OXYGEN, AND ARGON
                          
The files for this program are in the directory \gasp on the CD-ROM
  readme.txt      general description
  g.for           the source code of the subroutine package
  sample1.for     a sample program - table of thermo. props.
  sample1.out     expected results from running sample1
  sample2.for     a sample program - table of transport properties
  sample2.out     expected results from running sample2
  problem1.f90    Example problem #1 from the report
  problem1.out    expected results from running problem1
  problem2.f90    Example problem #2 from the report
  problem2.out    expected results from running problem2
  problem3.f90    Example problem #3 from the report
  problem3.out    expected results from running problem3
  problem4.f90    Example problem #4 from the report
  problem4.out    expected results from running problem4
  lew11629.txt    the original program description from COSMIC
  original.src    the original source code (from COSMIC)
  tnd7808.pdf     NASA TN D-7808 by Hendricks, Baron, and Peller

The reference documents for this software may be accessed
from the web page http://www.pdas.com/gasprefs.html.

Note that Gasp is not a program, but a subroutine package to be used
in your own work. The little programs program1.f90, etc. are examples
of the use of Gasp.
  
This program is quite well known and can be incorporated into your
programs if you need thermodynamic or transport properties. The samples 
and problems that were not included with the COSMIC distribution are from
the NASA report TN D-7808. GASP was released in 1971 as Cosmic Program 
11629 from NASA Lewis Research Center.

You might note a bug in the printed report for the output for Problem 2.
The exhaust temperature is 188.437 kelvins, but the printed value is 88.43.
This happened because the output format was F5.2 and no error was noted.