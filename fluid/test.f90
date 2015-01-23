! PROGRAM TestFluidMethods
! ------------------------------------------------------------------------------
! PURPOSE - Test the interpolation procedure NTRP and the CUBIC solver.
!
! AUTHOR - Ralph L. Carmichael, Public Domain Aeronautical Software

! REVISION HISTORY
!   DATE  VERS PERSON  STATEMENT OF CHANGES
! 25May10  0.5   RLC   Original coding

INCLUDE 'futil.f90'
INCLUDE 'fillupmod.f90'

!+
MODULE TestProcedures
! ------------------------------------------------------------------------------
! PURPOSE - 

USE FluidUtilities
USE Fillup
IMPLICIT NONE

  CHARACTER(LEN=*),PARAMETER:: VERSION = "Version 0.51 (30 May 2010)"

!-------------------------------------------------------------------------------

CONTAINS

!+
SUBROUTINE EstimateDerivative(x,f,fprime)
! ------------------------------------------------------------------------------
! PURPOSE - Fill fprime with an estimate of df/dx, where x and f are tables and
!  f = {(x(i),f(x(i))) : k=1,n}
  REAL,INTENT(IN),DIMENSION(:)::  x,f
  REAL,INTENT(OUT),DIMENSION(:):: fprime

  INTEGER:: k,n
!-------------------------------------------------------------------------------
  n=SIZE(x)
  fprime(1)=(f(2)-f(1))/(x(2)-x(1))                         ! forward difference
  fprime(n)=(f(n)-f(n-1))/(x(n)-x(n-1))                    ! backward difference

  DO k=2,n-1
    fprime(k)=(f(k+1)-f(k-1))/(x(k+1)-x(k-1))               ! central difference
  END DO
  RETURN
END Subroutine EstimateDerivative   ! ------------------------------------------

!+
SUBROUTINE InterpTester()
! ------------------------------------------------------------------------------

  REAL,PARAMETER,DIMENSION(18):: T = (/ 0.42, 0.45, 0.50, 0.59, 0.68, 0.80, &
    0.92, 0.985, 1.05, 1.15, 1.30, 1.60, 2.0, 2.6, 3.5, 4.6, 6.0, 8.0 /) 
  REAL,PARAMETER,DIMENSION(18):: D = (/ 1.0000E-5, 3.0000E-4, 7.0000E-4, 2.0000E-3, &
    0.01, 0.03, 0.1, 0.3, 0.6, 1.0, &
    1.45, 1.85, 2.20, 2.46, 2.65, 2.79, 2.87, 2.92 /)

  REAL,DIMENSION(18):: FDUM
  REAL,DIMENSION(18,99):: FFDUM

  INTEGER,PARAMETER:: NDIM = 100
  REAL,DIMENSION(NDIM):: x,f,fp
  INTEGER:: ia=0, iout
  INTEGER:: igoto
  INTEGER:: k
  REAL:: fdummy
!-------------------------------------------------------------------------------
  FFDUM=0.0
  CALL FillArray(T(1),t(SIZE(T)),x)
  DO k=1,NDIM
    igoto=1
    ia=0
    CALL NTRP(T,x(k),18,ia,iout,igoto,D,FFDUM,fdummy)
!!!    write(*,*) k,iout
    igoto=4
    CALL NTRP(T,x(k),18,ia,iout,igoto,D,FFDUM,f(k))
!!    write(*,*) k,iout
  END DO
  CALL EstimateDerivative(x,f,fp)
  
  WRITE(*,'(3F12.5)') (x(k),f(k),fp(k),k=1,NDIM)
  RETURN
END Subroutine InterpTester   ! ------------------------------------------------

END Module TestProcedures   ! ==================================================

!+
PROGRAM TestFluidMethods
! ------------------------------------------------------------------------------
USE TestProcedures
IMPLICIT NONE


!-------------------------------------------------------------------------------
  WRITE(*,*) "test"
  WRITE(*,*) VERSION  
  
  CALL InterpTester()
  
  
  STOP
END Program TestFluidMethods   ! ===============================================
