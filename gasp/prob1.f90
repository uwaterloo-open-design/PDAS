!INCLUDE 'g.for'
INCLUDE 'ga.f90'
!+    
PROGRAM PROBLEM1
! ---------------------------------------------------------------------------    
! PURPOSE - Calculate the heat-transfer coefficient for liquid hydrogen
!   flowing turbulently through a tube for these conditions:
!   (1) Tube diameter, 0.8 cm
!   (2) Average bulk velocity, 60 m/sec
!   (3) Average bulk temperature, 25 K
!   (4) Average bulk pressure, 6 megapascals
!   (5) Tube wall temperature, 450 K

! The heat transfer coefficient is the thermal conductivity times the
!   Nusselt number divided by the tube diameter
!
! The Nusselt number is 0.021 times Reynolds number to the 0.8 power
!   times the Prandtl number to the 0.4 power times SOMETHING
    
!    where
!    Nu  Nusselt number, hd/k
!    Re  Reynolds number, pVd/u
!    Pr  Prandtl number, gC P A
!    h   desired heat-transfer coefficient
!    d   tube diameter
!    k   thermal conductivity
!    P   density
!    ti  dynamic viscosity
!    V   kinematic viscosity
!    Cp  specific heat
!    U   average bulk velocity
!    
!    and the subscripts b, f, and w denote bulk, film, and
!      wall conditions, respectively.


!    THE GASP NAG STATEMENTS:
!    DO NOT USE K AS AN INDEX; IT MEANS THERMAL CONDUCTIVITY
!    DO NOT FORGET YOUR COMMON /PROPTY/ STATEMENT
!    DO NOT FORGET TO CALL SETUP
USE Gasppdas
  IMPLICIT NONE

!  REAL:: mu,mul,muv,k,kl,kv
!  INTEGER:: ku
!  REAL:: dl,dv, hl,hv, s,sl,sv, cv,cvl,cvv
!  REAL:: cp,cpl,cpv, gamma,gammal,gammav, c,cl,cvp
!  REAL:: sigma,excesk,excl,excv
!  COMMON/PROPTY/ku,dl,dv,hl,hv,s,sl,sv,cv,cvl,cvv,cp,cpl,cpv,gamma, &
!    gammal, gammav,c,cl,cvp,mu,mul,muv,k,kl,kv, sigma, excesk, excl,excv

  CHARACTER(LEN=3),PARAMETER:: NAM= "H2 "
  REAL,PARAMETER:: dia = 0.8    ! diameter, cm.
  REAL,PARAMETER:: vb = 60.0    ! bulk velocity, m/sec
  REAL:: tb = 25.0    ! bulk temperature, K
  REAL:: pb = 6.0     ! bulk pressure, megapascals
  REAL:: tw = 450.0   ! wall temperature K

  INTEGER:: kr,ks,kp
  REAL:: anub,anuw,anuf
  REAL:: db,dw,df
  REAL:: hb,hw,hf
  REAL:: ref,prf
  REAL:: tf
  REAL:: hcof

!----------------------------------------------------------------------------

  CALL SETUP(NAM)

  ku=1   ! select units of calculation
  ks=1                        ! compute bulk properties of density, viscosity
  kp=8
  kr=0 
  CALL GasProperties(ks,kp,tb,pb,db,hb,kr)

  anub=mu/db
  WRITE(*,*) "dia=", dia
  WRITE(*,*) "vb=", vb
  WRITE(*,*) "tb=", tb
  WRITE(*,*) "pb=", pb
  WRITE(*,*) "bulk density=", db
  WRITE(*,*) "bulk viscosity=", mu
  WRITE(*,*) "bulk kinematic viscosity=", anub
  WRITE(*,*)

  kr=0                    ! determine wall propertiesof density and viscosity
  ks=1
  kp=8
  CALL GasProperties(ks,kp,tw,pb,dw,hw,kr)
  anuw=mu/dw

  WRITE(*,*) "wall density=", dw
  WRITE(*,*) "wall viscosity=", mu
  WRITE(*,*) "wall kinematic viscosity=", anuw
  WRITE(*,*)

  kr=0  ! determine film properties of density, specific heat, viscosity
  kp=28   !  and thermal conductivity   KP=4+8+16
  ks=1
  TF=0.5*(TW+TB)

  CALL GasProperties(ks,kp,tf,pb,df,hf,kr)

  ref=df*vb*100.0*dia/mu
  prf=cp*mu/k
  anuf=0.021*ref**.8*prf**.4*(1.+0.01453*anuw/anub)
  hcof=anuf*k/dia

  WRITE(*,*) "film density=", df
  WRITE(*,*) "film viscosity=", mu
  WRITE(*,*) "film kinematic viscosity=", mu/df
  WRITE(*,*) "film thermal conductivity=", k
  WRITE(*,*)
  WRITE(*,*) "Reynolds number=", ref
  WRITE(*,*) "Prandtl number=", prf
  WRITE(*,*) "Nusselt factor=", anuf
  WRITE(*,*) "Coeff. of thermal conductivity=", hcof

  WRITE(*,*) "End of Problem 1"
  STOP 
END Program Problem1   ! ====================================================
