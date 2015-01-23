INCLUDE 'ga.f90'
!+
PROGRAM Problem4
! ---------------------------------------------------------------------------
! PURPOSE - Design a liquid-oxygen-cooled, copper-liner-type thrust
!   chamber where it is required to know the variation of enthalpy with
!   pressure at nearly constant temperature. Given the proposed design
!   parameters of an inlet temperature of 172 degrees Rankine and an
!   inlet pressure of 6700 psia, it is necessary to find the partial
!   derivative of enthalpy with respect to pressure at constant
!   temperature over a range in temperature and pressure for use in a
!   solution matrix.

!   The partial derivative is not directly available from GASP;
!   however, the partial derivative of pressure with respect to
!   temperarture at constant density and the partial derivative of
!   pressure with respect to density at constant temperature are both
!   available/ With a little manipulation one can find the
!   desired quantity

!    The derivative can now be found from GASP as follows:
!       (1) The desired units are psia, 0 R, and lbm/ft 3 (KU=3).
!       (2) The input is T and P. The output will be partial of pressure
!           w.r.t temperature at constant density and the partial of
!           pressure w.r.t density at constant temperature.
!        Set KS=1, KP=4, and KR=0
!       (3) The units of the first partial will be ft3/lbm, and
!           must be multiplied by 144/778.161 to return the KU=3 units,
    
!  An alternate procedure is to use the units of the program since the
! quotient is dimensionless
    
!  (1) Assign KU=I, with KP and KR as before.
!  (2) Convert P and T prior to calling GASP. 
!        P1=P/14.696; T1=T/1.8
    
!        CALL GASP(KU,KP,Tl,p1AH,KR)
    
!   (3) Then convert to the desired units.
!          Note that the units of (aH/aP )T are cm 3 /9
!    and J/(g)(MPa) with no further conversion required.
USE GaspPdas
IMPLICIT NONE
    
  REAL:: d,h
      REAL hpt(9)
  INTEGER:: i,j, kp,kr,ks
      REAL,PARAMETER,DIMENSION(5):: p = &
        (/ 4000.0, 5000.0, 6000.0, 6700.0, 7000.0 /)
      REAL,PARAMETER,DIMENSION(9):: t = &
        (/ 172.0, 200.0, 225.0, 250.0, 275.0, 300.0, 325.0, 350.0, 372.0 /)
  REAL:: plocal,tlocal

!      COMMON/PROPTY/ku,dl,dv,hl,hv,s,sl,sv,cv,cvl,cvv,cp,cpl,cpv,gamma, &
!     & gammal,gammav,c,cl,cvp,mu,mul,muv,k,kl,kv,slgma,excesk,excl,excv

!      COMMON/DERIV/pdt,ptv,pdtl,pdtv,ptvl,ptvv

!      COMMON/CON123/ dconv(5), tconv(5), pconv(5)
!      REAL mu, mul, muv, k, kl, kv
  CHARACTER(LEN=3),PARAMETER:: NAMGAS = "O2 "
!----------------------------------------------------------------------------
        CALL SETUP(NAMGAS)
        ku=3
        kp=4
        ks=1
        kr=0
        WRITE(6,*) "TABLE OF (dH/dT) FOR SELECTED TEMPERATURES AND PRESSURES"
        WRITE(6,'(6X,9F8.1)' ) (t(j),j=1,9)
        DO i=1,5
          DO j=1,9
            kr=0
            plocal=p(i)
            tlocal=t(j)
            CALL GasProperties(ks,kp,tlocal,plocal,d,h,kr)
            hpt(j)=(1.-t(j)/d*ptv/pdt*dconv(3)/tconv(3))/d
          END DO
          WRITE(6,'(F6.0,9F8.4)' ) p(i),(hpt(j),j=1,9)
        END DO


! now, an alternate approach...
        ku=1
        kp=4
        ks=1
        kr=0
        WRITE(6,*) "TABLE OF (dH/dT) FOR SELECTED TEMPERATURES AND PRESSURES"
        WRITE(6,'(6X,9F8.1)' ) (t(j),j=1,9)
        DO i=1,5
          DO j=1,9
            kr=0
            plocal=p(i)/14.696
            tlocal=t(j)/1.8
            CALL GasProperties(ks,kp,tlocal,plocal,d,h,kr)
            hpt(j)=(1.-t(j)/d*ptv/pdt)/d
          END DO
          WRITE(6,'(F6.0,9F8.4)' ) p(i),(hpt(j),j=1,9)
        END DO


  WRITE(*,*) "End of Problem 4"
  STOP 
END Program Problem4   ! ====================================================


