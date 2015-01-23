INCLUDE 'g.for'
!+
PROGRAM Problem5
! ---------------------------------------------------------------------------
! PURPOSE - Suppose that an odd set of units such as that given here are,
!   for some unknown reason, required in a user's program and that the
!   user also requires the use of GASP for PVT values.
    
!    QuantityNew system  Equivalent in Sl system
!         of units       of units
!    Length Cubit (cb)    45. 72 cm
!    Mass  Stone (at)     14 lbm or 63 50. 293 18 g
!    Time  Millenium      1000 yr or 3. 1556926X1010 see
!    Energy Megaton       4. 240 16 1
!    Temperature   Electron volt1. 6021xlO- 19 1. 1641940 4 K
!           equivalent (eve)1.380622xlO-22
!    Force Megapoundal    0. 1382549121376 MN
    
!    Converting this system of units so that GASP can calculate PVT values requires the
!    following conversion factors: 15
    
!       To convert from-      To-Multiply by -
    
!    Pressure, megapoundal/eb 2   MN/Tu 20.138254954376 = 0.661406197
!                (or MPa)         (0.4572) 2
!    Temperature, eVe K 1.16419xlO 4
!    Density, st/eb 3 g/cm3       6350.29318 = 3.037955
!                      (45.72) 3
!    Enthalpy, megaton/st J/g     4.2xlO 16 = 6.6138679xlO 12
!                     6350.29318
!    Entropy and specific J/g-K   -          4.2X10 16 5.6810897xlO 8
!    heats, megaton/st-eVe        6350. 293 18 x 1. 16419xio 4
!    
!    Sonic velocity, cb/millenium cm/sec,    45.72     = 1.4488103xlO-9
!                 3.1556926xlO 10
!    
!    Surface tension,      dyne/cm2-,-138254954376xlO 11 3.30239491xlO 9
!    megapoundal/eb                  45.72
!    Thermal conductivity,  W/cm-sec    -4.2xlO 16     2. 500486
!    megaton/cb -eve -millenium         45. 72 x 1. 16419x104 x 3.1556926xlO 10
!    
!     cosity, st/eb-millenium g/cm-sec    6350.29318 _= 4.4014205xlO-9
!              45. 72 x 3. 1556926X10 10
!    
!      15A convenient source for units conversion is
!    Mechtly, E. A.: The International System of Units, Physical Constants and Conversion Factors.
!    NASA SP-7012 (revised), 1973.
!    
!----------------------------------------------------------------------------



USE GaspPdas
IMPLICIT NONE

  WRITE(*,*) "End of Problem 5"
  STOP
END Program Problem5   ! ====================================================
