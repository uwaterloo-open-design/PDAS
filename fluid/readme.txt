

FLUID- THERMODYNAMIC AND TRANSPORT PROPERTIES OF FLUIDS        \fluid\readme.txt

The files for this program are in the directory \fluid on the CD-ROM
  readme.txt      general description
  instr.txt       Original instructions from COSMIC
  fluid.f90       the source code of the subroutine package
  ar.f90          tables for argon
  ch4.f90         tables for methane
  co2.f90         tables for carbon dioxide
  dryair.f90      tables for air
  f2.f90          tables for fluorine
  n2.f90          tables for nitrogen
  o2.f90          tables for oxygen
  ph2.f90         tables for parahydrogen
  lew14418.txt    the original program description from COSMIC
  steam.f90       tables for steam
  original.zip    the original source code (10 files) from COSMIC
  sample.for      original sample program from COSMIC 
  samprlc.f90     a variation on SAMPLE 


The reference documents for this software may be accessed
from the web page http://www.pdas.com/fluidrefs.html.

Note that Fluid is not a program, but a subroutine package for you
to include in your own work. The little programs sample.for and
samprlc.f90 are examples of the use of Fluid.

From the COSMIC description...

The accurate computation of the thermodynamic and transport properties of 
fluids is a necessity for many engineering calculations. The FLUID program 
was developed to calculate the thermodynamic and transport properties of 
pure fluids in both the liquid and gas phases. Fluid properties are 
calculated using a simple gas model, empirical corrections, and an efficient 
numerical interpolation scheme. FLUID produces results that are in very good 
agreement with measured values, while being much faster than older more 
complex programs developed for the same purpose.

A Van der Waals equation of state model is used to obtain approximate state 
values. These values are corrected for real gas effects by model correction 
factors obtained from tables based on experimental data. These tables also 
accurately compensate for the special circumstances which arise whenever 
phase conditions occur. Viscosity and thermal conductivity values are 
computed directly from tables. Interpolation within tables is based on 
Lagrange's three point formula. A set of tables must be generated for each 
fluid implemented.

FLUID currently contains tables for nine fluids including dry air and steam. 
The user can add tables for any fluid for which adequate thermal property 
data is available. The FLUID routine is structured so that it may easily be 
incorporated into engineering programs. FLUID was released in 1977 as Cosmic 
Program LEW-14418 from NASA Lewis Research Center.

