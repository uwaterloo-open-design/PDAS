!+
MODULE IntUtil   ! InterpolationUtilities
! ---------------------------------------------------------------------------
! PURPOSE - Provide various utility functions common to all
!    interpolation computations
!
! AUTHOR - Ralph L. Carmichael, Public Domain Aeronautical Software
!
! REVISION HISTORY
!   DATE  VERS PERSON  STATEMENT OF CHANGES
! 07May97  0.1   RLC   Original coding
! 09May97  0.2   RLC   Added GetCoeff1,GetCoeff2
! 13May97  0.3   RLC   Added PolyInterp
! 27May97  0.4   RLC   Renamed several items
! 02Jun97  0.5   RLC   Added SortAndIndex
! 07Jun97  0.6   RLC   Added TableLookup
! 27Oct97  0.7   RLC   Added SortAndIndexSingle
! 09Feb00  0.7   RLC   Fixed bug in TableLookup
! 23May01  0.75  RLC   Added LookupDecreasing
! 12Mar04  0.8   RLC   Added TableLookupDecreasing
! 20Mar04  0.9   RLC   Added LinearInterpolation1D, 2D, 3D
!
! USAGE NOTES-
!   CheckOrdering - Inspect a one-dimensional array to determine
!     if it is strictly increasing, strictly decreasing, non-decreasing,
!     non-increasing, or none of the above. Many interpolation schemes
!     depend upon this information.
!
!   Lookup - Search a sorted (increasing) array to find the interval
!     bounding a given number. If n is the size of the array a,
!     return 0 if number x is < a(1). Return n if x > a(n).
!     Return i if  a(i) <= x < a(i+1). If x is exactly a(n), return n-1.
!
!   LookupDecreasing - Search a sorted (decreasing) array to find the
!     interval bounding a given number. If n is the size of the array a,
!     return 0 if number x is > a(1). Return n if x < a(n).
!     Return i if  a(i) >= x > a(i+1). If x is exactly a(n), return n-1.
!
!   EvaluatePolynomial - Given an array with the coefficients of a
!     polynomial, evaluate that polynomial at a specified argument.
!     If the coefficient array is a, then a(1) is the constant term,
!     a(2) is the linear term, etc.
!
!   InterpolatePolynomial - Given a set of (x,y) pairs, there is a
!     unique polynomial that interpolates this data. Evaluate this
!     polynomial at a specified value. This routine is a straightforward
!     implementation of Lagrange's algorithm.
!
!   Neville - Performs exactly the same function as InterpolatePolynomial.
!     Uses Neville's algorithm instead of Lagrange. Based very strongly
!     on the routine in Numerical Recipes.
!
!   GetCoeff1 - Compute the coefficients of the cubic polynomial defined
!                 by its value and first derivatives at two points
!
!   GetCoeff2 - Compute the coefficients of the cubic polynomial defined
!                 by its value and second derivatives at two points
!
!   SortAndIndex - Sort an array into increasing order and create an
!      index array that will allow associated arrays to be rearranged
!      with the same reordering.
! 
!   TableLookup - Given a strictly increasing array of x-coordinates and
!     a corresponding array of y-coordinates, an x-coordinate for
!     evaluation and an order of interpolation, determine a value of the
!     interpolated value of x based on a subset of points surrounding the
!     argument. If TableLookup is used to create a function over the range
!     of the original x-array, it will be continuous, but not necessarily
!     smooth.

!   Linterp - Same as TableLookup, but order is 1 (linear)

IMPLICIT NONE

!----------------------------------------------------------------------------

  CHARACTER(LEN=*),PARAMETER:: MODULE_INTUTIL_VERSION = "0.75 (23May01)"
  INTEGER,PARAMETER,PRIVATE:: SP = KIND(1.0)        ! following example of
  INTEGER,PARAMETER,PRIVATE:: DP = KIND(1.0D0)      !  Numerical Recipes

INTERFACE CheckOrdering
  MODULE PROCEDURE CheckOrderingSingle, CheckOrderingDouble
END INTERFACE

INTERFACE Lookup
  MODULE PROCEDURE SingleLookup,DoubleLookup
END INTERFACE

INTERFACE LookupDecreasing
  MODULE PROCEDURE SingleLookupDecreasing,DoubleLookupDecreasing
END INTERFACE

INTERFACE EvaluatePolynomial
  MODULE PROCEDURE EvaluatePolynomialSingle,EvaluatePolynomialDouble
END INTERFACE

INTERFACE InterpolatePolynomial
  MODULE PROCEDURE InterpolatePolynomialSingle,InterpolatePolynomialDouble
END INTERFACE

INTERFACE Neville
  MODULE PROCEDURE NevilleSingle,NevilleDouble
END INTERFACE

INTERFACE GetCoeff1
  MODULE PROCEDURE GetCoeff1Single,GetCoeff1Double
END INTERFACE

INTERFACE GetCoeff2
  MODULE PROCEDURE GetCoeff2Single,GetCoeff2Double
END INTERFACE

INTERFACE LinearInterpolation1D
  MODULE PROCEDURE LinearInterpolation1DSingle,LinearInterpolation1DDouble
END INTERFACE

INTERFACE LinearInterpolation2D
  MODULE PROCEDURE LinearInterpolation2DSingle,LinearInterpolation2DDouble
END INTERFACE

INTERFACE LinearInterpolation3D
  MODULE PROCEDURE LinearInterpolation3DSingle,LinearInterpolation3DDouble
END INTERFACE

INTERFACE TableLookup
  MODULE PROCEDURE TableLookupSingle,TableLookupDouble
END INTERFACE

INTERFACE TableLookupDecreasing
  MODULE PROCEDURE TableLookupDecreasingSingle,TableLookupDecreasingDouble
END INTERFACE

INTERFACE SortAndIndex
  MODULE PROCEDURE SortAndIndexDouble, &
     SortAndIndexSingle, SortAndIndexInteger
END INTERFACE


CONTAINS

!+
FUNCTION CheckOrderingDouble(a) RESULT(k)
! ---------------------------------------------------------------------------
! PURPOSE - Inspect a one-dimensional array array to determine if it is
!  strictly increasing, strictly decreasing, non-decreasing, non-increasing,
!  or none of the above. Many interpolation schemes depend upon this
!  information.

  REAL(DP),INTENT(IN),DIMENSION(:):: a
  INTEGER:: k  ! =1 if a is strictly increasing
               ! =2 if a is strictly decreasing
               ! =3 if a is non-decreasing
               ! =4 if a is non-increasing
               ! =0 otherwise

  INTEGER:: i,n
!----------------------------------------------------------------------------
  n=SIZE(a)
  IF ( n < 2 ) THEN
    k=1        ! call is increasing
    RETURN
  END IF

  DO i=2,n
    IF ( a(i) <= a(i-1) ) EXIT
  END DO
  IF (i==n+1) THEN
    k=1              ! array is strictly increasing
    RETURN
  END IF

  DO i=2,n
    IF ( a(i) >= a(i-1) ) EXIT
  END DO
  IF (i==n+1) THEN
    k=2              ! array is strictly decreasing
    RETURN
  END IF

  DO i=2,n
    IF ( a(i) < a(i-1) ) EXIT
  END DO
  IF (i==n+1) THEN
    k=3              ! array is non-decreasing
    RETURN
  END IF

  DO i=2,n
    IF ( a(i) > a(i-1) ) EXIT
  END DO
  IF (i==n+1) THEN
    k=4              ! array is non-increasing
    RETURN
  END IF

  k=0    ! if it fails all the above tests
  RETURN
END Function CheckOrderingDouble   ! ----------------------------------------

!+
FUNCTION CheckOrderingSingle(a) RESULT(k)
! ---------------------------------------------------------------------------
! PURPOSE - Same as CheckOrderingDouble.
  REAL(SP),INTENT(IN),DIMENSION(:):: a
  INTEGER:: k
  REAL(DP),DIMENSION(SIZE(a)):: ad
!----------------------------------------------------------------------------
  ad=a   ! converts single to double
  k=CheckOrderingDouble(ad)
  RETURN
END Function CheckOrderingSingle   ! ----------------------------------------

!+
FUNCTION SingleLookup(xtab,x) RESULT (i)
! ---------------------------------------------------------------------------
! PURPOSE - Search a sorted (increasing) array to find the interval
!  bounding a given number. If n is the size of the array a,
!  return 0 if number x is < a(1)
!  return n if x > a(n).
!  return i if  a(i) <= x < a(i+1).
!  If x is exactly equal to a(n), return n-1.

  REAL(SP),INTENT(IN),DIMENSION(:)::  xtab    ! input array                        
  REAL(SP),INTENT(IN):: x  ! input number              

  INTEGER:: i  ! index of interval such that xtab(i) <= x < xtab(i+1)
               ! =0 if x < xtab(1) and =n if x > xtab(i)
  INTEGER:: j,k,n
  INTEGER,SAVE:: isave = 1
!----------------------------------------------------------------------------
  n=SIZE(xtab)
  IF (n <= 0) THEN
    i=-1
    RETURN
  END IF

  IF (x < xtab(1)) THEN
    i=0
    RETURN
  END IF

  IF (x > xtab(n)) THEN
    i=n
    RETURN
  END IF

  IF (isave>0 .AND. isave<n) THEN
    IF ((x >= xtab(isave)) .AND. (x < xtab(isave+1)) ) THEN
      i=isave
      RETURN
    END IF
  END IF

  i=1 
  j=SIZE(xtab)
  DO
    IF (j <= i+1) EXIT
    k=(i+j)/2                                     ! integer division
    IF (x < xtab(k)) THEN
      j=k
    ELSE
      i=k
    END IF      
  END DO
  isave=i

  RETURN
END Function SingleLookup   ! ===============================================

!+
FUNCTION DoubleLookup(xtab,x) RESULT (i)
! ---------------------------------------------------------------------------
! PURPOSE - Same as Single Lookup

  REAL(DP),INTENT(IN),DIMENSION(:)::  xtab    ! input array                        
  REAL(DP),INTENT(IN):: x  ! input number              

  INTEGER:: i  ! index of interval such that xtab(i) <= x < xtab(i+1)
               ! =0 if x < xtab(1) and =n if x > xtab(i)
  INTEGER:: j,k,n
  INTEGER,SAVE:: isave = 1
!----------------------------------------------------------------------------
  n=SIZE(xtab)
  IF (n <= 0) THEN
    i=-1
    RETURN
  END IF

  IF (x < xtab(1)) THEN
    i=0
    RETURN
  END IF

  IF (x > xtab(n)) THEN
    i=n
    RETURN
  END IF

  IF (isave>0 .AND. isave<n) THEN
    IF ((x >= xtab(isave)) .AND. (x < xtab(isave+1)) ) THEN
      i=isave
      RETURN
    END IF
  END IF

  i=1 
  j=SIZE(xtab)
  DO
    IF (j <= i+1) EXIT
    k=(i+j)/2                                     ! integer division
    IF (x < xtab(k)) THEN
      j=k
    ELSE
      i=k
    END IF         
  END DO
  isave=i

  RETURN
END Function DoubleLookup   ! -----------------------------------------------

!+
FUNCTION SingleLookupDecreasing(xtab,x) RESULT (i)
! ---------------------------------------------------------------------------
! PURPOSE - Search a sorted (decreasing) array to find the interval
!  bounding a given number. If n is the size of the array xtab,
!  return 0 if number x is > xtab(1)
!  return n if x < xtab(n).
!  return i if  xtab(i) >= x > xtab(i+1).
!  If x is exactly equal to xtab(n), return n-1.

  REAL(SP),INTENT(IN),DIMENSION(:)::  xtab    ! input array                        
  REAL(SP),INTENT(IN):: x  ! input number              

  INTEGER:: i  ! index of interval such that xtab(i) >= x > xtab(i+1)
               ! =0 if x > xtab(1) and =n if x < xtab(n)
  INTEGER:: j,k,n
!----------------------------------------------------------------------------
  n=SIZE(xtab)
  IF (n <= 0) THEN
    i=-1
    RETURN
  END IF

  IF (x > xtab(1)) THEN
    i=0
    RETURN
  END IF

  IF (x < xtab(n)) THEN
    i=n
    RETURN
  END IF

  i=1 
  j=SIZE(xtab)
  DO
    IF (j <= i+1) EXIT
    k=(i+j)/2                                     ! integer division
    IF (x > xtab(k)) THEN
      j=k
    ELSE
      i=k
    END IF      
  END DO

  RETURN
END Function SingleLookupDecreasing   ! =====================================


!+
FUNCTION DoubleLookupDecreasing(xtab,x) RESULT (i)
! ---------------------------------------------------------------------------
! PURPOSE - Search a sorted (decreasing) array to find the interval
!  bounding a given number. If n is the size of the array xtab,
!  return 0 if number x is > xtab(1)
!  return n if x < xtab(n).
!  return i if  xtab(i) >= x > xtab(i+1).
!  If x is exactly equal to xtab(n), return n-1.

  REAL(DP),INTENT(IN),DIMENSION(:)::  xtab    ! input array                        
  REAL(DP),INTENT(IN):: x  ! input number              

  INTEGER:: i  ! index of interval such that xtab(i) >= x > xtab(i+1)
               ! =0 if x > xtab(1) and =n if x < xtab(n)
  INTEGER:: j,k,n
!----------------------------------------------------------------------------
  n=SIZE(xtab)
  IF (n <= 0) THEN
    i=-1
    RETURN
  END IF

  IF (x > xtab(1)) THEN
    i=0
    RETURN
  END IF

  IF (x < xtab(n)) THEN
    i=n
    RETURN
  END IF

  i=1 
  j=SIZE(xtab)
  DO
    IF (j <= i+1) EXIT
    k=(i+j)/2                                     ! integer division
    IF (x > xtab(k)) THEN
      j=k
    ELSE
      i=k
    END IF      
  END DO

  RETURN
END Function DoubleLookupDecreasing   ! =====================================


!+
FUNCTION EvaluatePolynomialSingle(a, x) RESULT(sum)
! ---------------------------------------------------------------------------
! PURPOSE - Evaluate a polynomial with coeff array a

! NOTES - poly = a(1) + a(2)*x + a(3)*x**2 + a(4)*x**3 + ...
      
  REAL(SP),INTENT(IN),DIMENSION(:):: a
  REAL(SP),INTENT(IN):: x

  INTEGER:: i,n
  REAL(SP):: sum
!----------------------------------------------------------------------------
  n=SIZE(a)
  IF (n <= 0) THEN
    sum=0.0
  ELSE
    sum=a(n)
    DO i=n-1,1,-1          ! Horner's rule
      sum=a(i)+x*sum
    END DO
  END IF
  RETURN
END Function EvaluatePolynomialSingle   ! -----------------------------------

!+
FUNCTION EvaluatePolynomialDouble(a, x) RESULT(sum)
! ---------------------------------------------------------------------------
! PURPOSE - Evaluate a polynomial with coeff array a

! NOTES - poly = a(1) + a(2)*x + a(3)*x**2 + a(4)*x**3 + ...
      
  REAL(DP),INTENT(IN),DIMENSION(:):: a
  REAL(DP),INTENT(IN):: x

  INTEGER:: i,n
  REAL(DP):: sum
!----------------------------------------------------------------------------
  n=SIZE(a)
  IF (n <= 0) THEN
    sum=0.0
  ELSE
    sum=a(n)
    DO i=n-1,1,-1
      sum=a(i)+x*sum      ! Horner's rule
    END DO
  END IF
  RETURN
END Function EvaluatePolynomialDouble   ! -----------------------------------

!+
FUNCTION InterpolatePolynomialSingle(x,y,u) RESULT(sum)
! ---------------------------------------------------------------------------
! PURPOSE -Compute the value of the interpolating polynomial thru
!   x- and y-arrays at the x-value of u, using Lagrange's equation.

  REAL(SP),INTENT(IN),DIMENSION(:):: x,y   ! tables of coordinates
  REAL(SP),INTENT(IN):: u   ! value of x-coordinate for interpolation
  REAL(SP):: sum
  INTEGER:: i,j
  REAL(SP):: fact
  REAL(SP),DIMENSION(SIZE(x)):: du
!----------------------------------------------------------------------------
  du(:)=u-x(:)
  sum=0.0
  DO j=1,SIZE(x)
    fact=1.0
    DO i=1,SIZE(x)
      IF (i /= j) fact=fact*du(i)/(x(j)-x(i))
    END DO
    sum=sum+y(j)*fact
  END DO
  RETURN
END Function InterpolatePolynomialSingle   ! --------------------------------

!+
FUNCTION InterpolatePolynomialDouble(x,y,u) RESULT(sum)
! ---------------------------------------------------------------------------
! PURPOSE -Compute the value of the interpolating polynomial thru
!   x- and y-arrays at the x-value of u, using Lagrange's equation.

  REAL(DP),INTENT(IN),DIMENSION(:):: x,y   ! tables of coordinates
  REAL(DP),INTENT(IN):: u   ! value of x-coordinate for interpolation
  REAL(DP):: sum
  INTEGER:: i,j
  REAL(DP):: fact
  REAL(DP),DIMENSION(SIZE(x)):: du
!----------------------------------------------------------------------------
  du(:)=u-x(:)
  sum=0.0
  DO j=1,SIZE(x)
    fact=1.0
    DO i=1,SIZE(x)
      IF (i /= j) fact=fact*du(i)/(x(j)-x(i))
    END DO
   sum=sum+y(j)*fact
  END DO
  RETURN
END Function InterpolatePolynomialDouble   ! --------------------------------

!+
FUNCTION NevilleDouble(x,y,u) RESULT(fu)
! ---------------------------------------------------------------------------
! PURPOSE -Compute the value of the interpolating polynomial thru
!   x- and y-arrays at the x-value of u, using Neville's equation.
!   Based on subroutine polint in Numerical Recipes

  REAL(DP),INTENT(IN),DIMENSION(:):: x,y
  REAL(DP),INTENT(IN):: u
  REAL(DP):: fu

  REAL(DP),DIMENSION(SIZE(x)):: c,d,den,h
  REAL(DP):: dy                                  ! returned in the NR version
  INTEGER:: m,n,near
  INTEGER,DIMENSION(1):: ns      ! so we can use MINLOC, not iminloc
!----------------------------------------------------------------------------
  n=SIZE(x)
  c=y
  d=c
  h=x-u
  ns=MINLOC(ABS(u-x))   ! why not use h
  near=ns(1)
  fu=y(near)
  near=near-1
  DO m=1,n-1
    den(1:n-m)=h(1:n-m)-h(1+m:n)
    den(1:n-m)=(c(2:n-m+1)-d(1:n-m))/den(1:n-m)
    d(1:n-m)=h(1+m:n)*den(1:n-m)   ! updating the c and d arrays
    c(1:n-m)=h(1:n-m)*den(1:n-m)
    IF (2*near < n-m)  THEN
      dy=c(near+1)
    ELSE
      dy=d(near)
      near=near-1
    END IF
    fu=fu+dy
  END DO

  RETURN
END Function NevilleDouble   ! ----------------------------------------------

!+
FUNCTION NevilleSingle(x,y,u) RESULT(fu)
! ---------------------------------------------------------------------------
  REAL(SP),INTENT(IN),DIMENSION(:):: x,y
  REAL(SP),INTENT(IN):: u
  REAL(SP):: fu

  REAL(DP),DIMENSION(SIZE(x)):: xd,yd
  REAL(DP):: xxd,yyd
!----------------------------------------------------------------------------
  xd(:)=x(:)
  yd(:)=y(1:SIZE(x))
  xxd=u
  yyd=NevilleDouble(xd,yd,xxd)
  fu=yyd 
  RETURN
END Function NevilleSingle   ! ----------------------------------------------

!+
SUBROUTINE GetCoeff1Double(x1,x2,f1,f2,fp1,fp2, b)
! ---------------------------------------------------------------------------
! PURPOSE - Compute the coefficients of the cubic polynomial determined by
!    the function value and first derivatives at two points
  REAL(DP),INTENT(IN):: x1,x2
  REAL(DP),INTENT(IN):: f1,f2
  REAL(DP),INTENT(IN):: fp1,fp2
  REAL(DP),INTENT(OUT),DIMENSION(:):: b   ! dimension must be at least 4

  REAL(DP),DIMENSION(4):: a   ! coeff. of the cubic in t=x-x1
  REAL(DP):: dx,dy,dydx
!----------------------------------------------------------------------------
  dx=x2-x1
  dy=f2-f1
  dydx=dy/dx

  a(1)=f1
  a(2)=fp1
  a(3)=(3.0*dydx-2.0*fp1-fp2)/dx
  a(4)=(fp1+fp2-2.0*dydx)/dx**2

  b(1)=a(1) + x1*(-a(2) + x1*(a(3) - x1*a(4)))
  b(2)=a(2) + x1*(-2.0*a(3) + x1*3.0*a(4))
  b(3)=a(3)-3.0*a(4)*x1
  b(4)=a(4)
  RETURN
END Subroutine GetCoeff1Double   ! ------------------------------------------

!+
SUBROUTINE GetCoeff1Single(x1,x2,f1,f2,fp1,fp2, b)
! ---------------------------------------------------------------------------
! PURPOSE - Compute the coefficients of the cubic polynomial determined by
!    the function value and second derivatives at two points
  REAL(SP),INTENT(IN):: x1,x2
  REAL(SP),INTENT(IN):: f1,f2
  REAL(SP),INTENT(IN):: fp1,fp2
  REAL(SP),INTENT(OUT),DIMENSION(:):: b   ! dimension must be at least 4

  REAL(DP):: x1d,x2d, f1d,f2d, fp1d,fp2d
  REAL(DP),DIMENSION(4):: b2  
!----------------------------------------------------------------------------
  x1d=x1
  x2d=x2
  f1d=f1
  f2d=f2
  fp1d=fp1
  fp2d=fp2
  CALL GetCoeff1Double(x1d,x2d, f1d,f2d, fp1d,fp2d, b2)
  b(1:4)=b2(:)
  RETURN
END Subroutine GetCoeff1Single   ! ------------------------------------------

!+
SUBROUTINE GetCoeff2Double(x1,x2,f1,f2,fpp1,fpp2, b)
! ---------------------------------------------------------------------------
! PURPOSE - Compute the coefficients of the cubic polynomial determined by
!    the function value and second derivatives at two points
  REAL(DP),INTENT(IN):: x1,x2
  REAL(DP),INTENT(IN):: f1,f2
  REAL(DP),INTENT(IN):: fpp1,fpp2
  REAL(DP),INTENT(OUT),DIMENSION(:):: b   ! dimension must be at least 4

  REAL(DP),DIMENSION(4):: a   ! coeff. of the cubic in t=x-x1
  REAL(DP):: dx
!----------------------------------------------------------------------------
  dx=x2-x1
  a(1)=f1
  a(3)=0.5*fpp1
  a(4)=(fpp2-fpp1)/(6.0*dx)
  a(2)=(f2-f1)/dx - dx*(a(3)+dx*a(4))

  b(1)=a(1) + x1*(-a(2) + x1*(a(3) - x1*a(4)))
  b(2)=a(2) + x1*(-2.0*a(3) + x1*3.0*a(4))
  b(3)=a(3)-3.0*a(4)*x1
  b(4)=a(4)
  RETURN
END Subroutine GetCoeff2Double   ! ------------------------------------------

!+
SUBROUTINE GetCoeff2Single(x1,x2,f1,f2,fpp1,fpp2, b)
! ---------------------------------------------------------------------------
! PURPOSE - Compute the coefficients of the cubic polynomial determined by
!    the function value and second derivatives at two points
  REAL(SP),INTENT(IN):: x1,x2
  REAL(SP),INTENT(IN):: f1,f2
  REAL(SP),INTENT(IN):: fpp1,fpp2
  REAL(SP),INTENT(OUT),DIMENSION(:):: b   ! dimension must be at least 4

  REAL(DP):: x1d,x2d, f1d,f2d, fpp1d,fpp2d
  REAL(DP),DIMENSION(4):: b2  
!----------------------------------------------------------------------------
  x1d=x1
  x2d=x2
  f1d=f1
  f2d=f2
  fpp1d=fpp1
  fpp2d=fpp2
  CALL GetCoeff2Double(x1d,x2d, f1d,f2d, fpp1d,fpp2d, b2)
  b(1:4)=b2(:)
  RETURN
END Subroutine GetCoeff2Single   ! ------------------------------------------

!+
FUNCTION LinearInterpolation1DSingle(xtable,ftable,x) RESULT(f)
! ---------------------------------------------------------------------------
! PURPOSE - interpolate the tables xtable and ftable at the value of x

  REAL(SP),INTENT(IN),DIMENSION(:):: xtable,ftable
  REAL(SP),INTENT(IN):: x
  REAL(SP):: f
!----------------------------------------------------------------------------
  f=TableLookupSingle(xtable,ftable, 1,x)
  RETURN
END Function LinearInterpolation1DSingle   ! --------------------------------

!+
FUNCTION LinearInterpolation2DSingle(xtable,ytable,ftable,x,y) RESULT(z)
! ---------------------------------------------------------------------------
! PURPOSE - ftable has dimension (SIZE(xtable), SIZE(ytable)). This
!  function interpolates at (x,y) for the bilinear interpolant.

  REAL(SP),INTENT(IN),DIMENSION(:):: xtable,ytable
  REAL(SP),INTENT(IN),DIMENSION(:,:):: ftable
  REAL(SP),INTENT(IN):: x,y

  REAL(SP):: z
  INTEGER:: j
  INTEGER:: ny
  REAL(SP):: z1,z2
!----------------------------------------------------------------------------
  ny=SIZE(ytable)
  j=Lookup(ytable,y)
  IF (j==0) j=1
  IF (j==ny) j=ny-1
  z1=TableLookupSingle(xtable,ftable(:,j),   1,x)
  z2=TableLookupSingle(xtable,ftable(:,j+1), 1,x)
  z=z1 + (y-ytable(j))*(z2-z1)/(ytable(j+1)-ytable(j))
  RETURN
END Function LinearInterpolation2DSingle   ! --------------------------------

!+
FUNCTION LinearInterpolation3DSingle(xtable,ytable,ztable,ftable,x,y,z) &
                                                                   RESULT(f)
! ---------------------------------------------------------------------------
! PURPOSE - ftable has dimension (SIZE(xtable), SIZE(ytable), SIZE(ztable)).
!  This function interpolates at (x,y,z) for the bilinear interpolant.
  REAL(SP),INTENT(IN),DIMENSION(:):: xtable,ytable,ztable
  REAL(SP),INTENT(IN),DIMENSION(:,:,:):: ftable
  REAL(SP),INTENT(IN):: x,y,z

  REAL(SP):: f
  INTEGER:: j
  INTEGER:: nz
  REAL(SP):: f1,f2
!----------------------------------------------------------------------------
  nz=SIZE(ztable)
  j=Lookup(ztable,z)
  IF (j==0) j=1
  IF (j==nz) j=nz-1
  f1=LinearInterpolation2DSingle(xtable,ytable, ftable(:,:,j),   x,y)
  f2=LinearInterpolation2DSingle(xtable,ytable, ftable(:,:,j+1), x,y)
  f=f1 + (z-ztable(j))*(f2-f1)/(ztable(j+1)-ztable(j))
  RETURN
END Function LinearInterpolation3DSingle   ! --------------------------------



!+
FUNCTION LinearInterpolation1DDouble(xtable,ftable,x) RESULT(f)
! ---------------------------------------------------------------------------
! PURPOSE - interpolate the tables xtable and ftable at the value of x

  REAL(DP),INTENT(IN),DIMENSION(:):: xtable,ftable
  REAL(DP),INTENT(IN):: x
  REAL(DP):: f
!----------------------------------------------------------------------------
  f=TableLookupDouble(xtable,ftable, 1,x)
  RETURN
END Function LinearInterpolation1DDouble   ! --------------------------------

!+
FUNCTION LinearInterpolation2DDouble(xtable,ytable,ftable,x,y) RESULT(z)
! ---------------------------------------------------------------------------
! PURPOSE - ftable has dimension (SIZE(xtable), SIZE(ytable)). This
!  function interpolates at (x,y) for the bilinear interpolant.

  REAL(DP),INTENT(IN),DIMENSION(:):: xtable,ytable
  REAL(DP),INTENT(IN),DIMENSION(:,:):: ftable
  REAL(DP),INTENT(IN):: x,y

  REAL(DP):: z
  INTEGER:: j
  INTEGER:: ny
  REAL(DP):: z1,z2
!----------------------------------------------------------------------------
  ny=SIZE(ytable)
  j=Lookup(ytable,y)
  IF (j==0) j=1
  IF (j==ny) j=ny-1
  z1=TableLookupDouble(xtable,ftable(:,j),   1,x)
  z2=TableLookupDouble(xtable,ftable(:,j+1), 1,x)
  z=z1 + (y-ytable(j))*(z2-z1)/(ytable(j+1)-ytable(j))
  RETURN
END Function LinearInterpolation2DDouble   ! --------------------------------

!+
FUNCTION LinearInterpolation3DDouble(xtable,ytable,ztable,ftable,x,y,z) &
                                                                   RESULT(f)
! ---------------------------------------------------------------------------
! PURPOSE - ftable has dimension (SIZE(xtable), SIZE(ytable), SIZE(ztable)).
!  This function interpolates at (x,y,z) for the bilinear interpolant.
  REAL(DP),INTENT(IN),DIMENSION(:):: xtable,ytable,ztable
  REAL(DP),INTENT(IN),DIMENSION(:,:,:):: ftable
  REAL(DP),INTENT(IN):: x,y,z

  REAL(DP):: f
  INTEGER:: j
  INTEGER:: nz
  REAL(DP):: f1,f2
!----------------------------------------------------------------------------
  nz=SIZE(ztable)
  j=Lookup(ztable,z)
  IF (j==0) j=1
  IF (j==nz) j=nz-1
  f1=LinearInterpolation2DDouble(xtable,ytable, ftable(:,:,j),   x,y)
  f2=LinearInterpolation2DDouble(xtable,ytable, ftable(:,:,j+1), x,y)
  f=f1 + (z-ztable(j))*(f2-f1)/(ztable(j+1)-ztable(j))
  RETURN
END Function LinearInterpolation3DDouble   ! --------------------------------



!+
FUNCTION TableLookupSingle(x,y,order,u) RESULT(fu)
! ---------------------------------------------------------------------------
  REAL(SP),INTENT(IN),DIMENSION(:):: x,y
  INTEGER,INTENT(IN):: order   ! order of interpolation
                               ! (1=linear, 2=quadratic...)
  REAL(SP),INTENT(IN):: u
  REAL(SP):: fu
  INTEGER:: j,m
!----------------------------------------------------------------------------
  m=MIN(order+1, SIZE(x))  ! number of points used for interpolating poly
  j=Lookup(x,u)
  j=j-(m/2-1)
  j=MIN(1+SIZE(x)-m, j)                       ! j+m-1 must not exceed SIZE(x)
  j=MAX(1,j)                                  ! j must be positive
                       ! use points j thru j+m-1 for interpolation (m points)
  fu=InterpolatePolynomialSingle(x(j:j+m-1),y(j:j+m-1),u)
  RETURN
END Function TableLookupSingle   ! ------------------------------------------

!+
FUNCTION TableLookupDouble(x,y,order,u) RESULT(fu)
! ---------------------------------------------------------------------------
  REAL(DP),INTENT(IN),DIMENSION(:):: x,y
  INTEGER,INTENT(IN):: order   ! order of interpolation
                               ! (1=linear, 2=quadratic...)
  REAL(DP),INTENT(IN):: u
  REAL(DP):: fu
  INTEGER:: j,m
!----------------------------------------------------------------------------
  m=MIN(order+1, SIZE(x))  ! number of points used for interpolating poly
  j=Lookup(x,u)
  j=j-(m/2-1)
  j=MIN(1+SIZE(x)-m, j)
  j=MAX(1,j)
  fu=InterpolatePolynomialDouble(x(j:j+m-1),y(j:j+m-1),u)
  RETURN
END Function TableLookupDouble   ! ------------------------------------------


!+
FUNCTION TableLookupDecreasingSingle(x,y,order,u) RESULT(fu)
! ---------------------------------------------------------------------------
  REAL(SP),INTENT(IN),DIMENSION(:):: x,y
  INTEGER,INTENT(IN):: order   ! order of interpolation
                               ! (1=linear, 2=quadratic...)
  REAL(SP),INTENT(IN):: u
  REAL(SP):: fu
  INTEGER:: j,m
!----------------------------------------------------------------------------
  m=MIN(order+1, SIZE(x))  ! number of points used for interpolating poly
  j=LookupDecreasing(x,u)
  j=j-(m/2-1)
  j=MIN(1+SIZE(x)-m, j)                       ! j+m-1 must not exceed SIZE(x)
  j=MAX(1,j)                                  ! j must be positive
                       ! use points j thru j+m-1 for interpolation (m points)
  fu=InterpolatePolynomialSingle(x(j:j+m-1),y(j:j+m-1),u)
  RETURN
END Function TableLookupDecreasingSingle   ! ------------------------------------------

!+
FUNCTION TableLookupDecreasingDouble(x,y,order,u) RESULT(fu)
! ---------------------------------------------------------------------------
  REAL(DP),INTENT(IN),DIMENSION(:):: x,y
  INTEGER,INTENT(IN):: order   ! order of interpolation
                               ! (1=linear, 2=quadratic...)
  REAL(DP),INTENT(IN):: u
  REAL(DP):: fu
  INTEGER:: j,m
!----------------------------------------------------------------------------
  m=MIN(order+1, SIZE(x))  ! number of points used for interpolating poly
  j=LookupDecreasing(x,u)
  j=j-(m/2-1)
  j=MIN(1+SIZE(x)-m, j)
  j=MAX(1,j)
  fu=InterpolatePolynomialDouble(x(j:j+m-1),y(j:j+m-1),u)
  RETURN
END Function TableLookupDecreasingDouble   ! ------------------------------------------

!+
SUBROUTINE SortAndIndexDouble(a, ix)
! ---------------------------------------------------------------------------
! PURPOSE - Sort and index a real array
!
! AUTHOR - Ralph L. Carmichael, Public Domain Aeronautical Software
!       with great assistance from the authors of "Numerical Recipes,
!       The Art of Scientific Computing". 
!
! NOTES - This is a replacement for the routine VSRTR in the IMSL library.
!   VSRTR is copyrighted, but this routine is public domain.
!   The routine uses the heapsort algorithm. The coding follows the
!   examples of the routines SORT2 and INDEXX in "Numerical Recipes"
!
!   On input, a contains the vector to be sorted.
!   On input, ix contains an identifying integer for each element
!   of a. Usually, one sets ix=1,2,3,4,....
!   For real generality, the values of ix may be anything that is useful
!   and meaningful to the programmer. Whenever two values of the array a
!   are swapped, the corresponding entries in the array ix are also swapped.
!   After execution, the a-array will be in increasing order and the array
!   ix will hold a history of the transformations
!
  REAL(DP),INTENT(IN OUT),DIMENSION(:):: a
  INTEGER,INTENT(IN OUT),DIMENSION(:):: ix

  INTEGER:: n
  INTEGER:: i,j,k, ir, ixsave
  REAL(DP):: asave
!----------------------------------------------------------------------------
  n=SIZE(a)   !   ix should also be as large
  IF (n <= 1) RETURN

  k=(n/2)+1
  ir=n

  DO
    IF (k > 1) THEN
      k=k-1
      asave=a(k)
      ixsave=ix(k)
    ELSE
      asave=a(ir)
      ixsave=ix(ir)
      a(ir)=a(1)
      ix(ir)=ix(1)
      ir=ir-1
      IF (ir == 1) THEN
        a(1)=asave
        ix(1) = ixsave
        EXIT            ! this is the only way out of here
      END IF
    END IF

    i = k
    j = k+k

    DO
      IF ( j > ir) EXIT   ! gets you out of this loop
      IF (j < ir) THEN
        IF ( a(j) < a(j+1) ) j=j+1
      END IF

      IF ( asave < a(j) ) THEN
        a(i)=a(j)
        ix(i) = ix(j)
        i = j
        j = j+j
      ELSE
        j = ir+1
      END IF
    END DO

    a(i)=asave
    ix(i) = ixsave
  END DO                                  ! only way out is thru RETURN above
  RETURN
END Subroutine SortAndIndexDouble   ! ---------------------------------------

!+
SUBROUTINE SortAndIndexSingle(a, ix)
! ---------------------------------------------------------------------------
! PURPOSE - Sort and index a real(SP) array

  REAL(SP),INTENT(IN OUT),DIMENSION(:):: a
  INTEGER,INTENT(IN OUT),DIMENSION(:):: ix

  REAL(DP),DIMENSION(SIZE(a)):: ad
!----------------------------------------------------------------------------
  ad(:)=a(:)                        ! converts SP to DP
  CALL SortAndIndexDouble(ad,ix)
  a(:)=ad(:)                        ! converts DP to SP
  RETURN
END Subroutine SortAndIndexSingle   ! ---------------------------------------

!+
SUBROUTINE SortAndIndexInteger(a, ix)
! ---------------------------------------------------------------------------
! PURPOSE - Sort and index an integer array
!
! AUTHOR - Ralph L. Carmichael, Public Domain Aeronautical Software
!       with great assistance from the authors of "Numerical Recipes,
!       The Art of Scientific Computing". 
!
! NOTES - This is a replacement for the routine VSRTI in the IMSL library.
!   VSRTI is copyrighted, but this routine is public domain.
!   The routine uses the heapsort algorithm. The coding follows the
!   examples of the routines SORT2 and INDEXX in "Numerical Recipes"
!
!   On input, a contains the vector to be sorted.
!   On input, ix contains an identifying integer for each element
!   of a. Usually, one sets ix=1,2,3,4,....
!   For real generality, the values of ix may be anything that is useful
!   and meaningful to the programmer. Whenever two values of the array a
!   are swapped, the corresponding entries in the array ix are also swapped.
!   After execution, the a-array will be in increasing order and the array
!   ix will hold a history of the transformations

  INTEGER,INTENT(IN OUT),DIMENSION(:):: a  ! just like VSRTR but a is integer
  INTEGER,INTENT(IN OUT),DIMENSION(:):: ix

  INTEGER:: n
  INTEGER:: i,j,k, ir, ixsave
  INTEGER:: asave
!----------------------------------------------------------------------------
  n=SIZE(a)
  IF (n <= 1) RETURN

      k=(n/2)+1
      ir=n

   20 IF (k .GT. 1) THEN
        k=k-1
        asave=a(k)
        ixsave=ix(k)
      ELSE
        asave=a(ir)
        ixsave=ix(ir)
        a(ir)=a(1)
        ix(ir)=ix(1)
        ir=ir-1
        IF (ir .EQ. 1) THEN
          a(1)=asave
          ix(1) = ixsave
          RETURN            ! this is the only way out of here
        END IF
      END IF

      i = k
      j = k+k

   30 IF ( j .LE. ir) THEN
        IF (j .LT. ir) THEN
          IF ( a(j) .LT. a(j+1) ) j=j+1
        END IF

        IF ( asave .LT. a(j) ) THEN
          a(i)=a(j)
          ix(i) = ix(j)
          i = j
          j = j+j
        ELSE
          j = ir+1
        END IF
        GO TO 30
      END IF

      a(i)=asave
      ix(i) = ixsave
      GO TO 20
END Subroutine SortAndIndexInteger   ! --------------------------------------

!+
SUBROUTINE IndexDouble(a, ix)
! ---------------------------------------------------------------------------
! PURPOSE - Index a real array
!
! AUTHOR - Ralph L. Carmichael, Public Domain Aeronautical Software
!       with great assistance from the authors of "Numerical Recipes,
!       The Art of Scientific Computing". 

!   The routine uses the heapsort algorithm. The coding follows the
!   examples of the routines SORT2 and INDEXX in "Numerical Recipes"
!   and the subroutine VSRTR(a, n, ix) from IMSL

!   On input, a contains the vector to be indexed
!   On input, ix contains an identifying integer for each element
!   of a. Usually, one sets ix=1,2,3,4,....
!   For real generality, the values of ix may be anything that is useful
!   and meaningful to the programmer. Whenever two values of the array a
!   would be swapped in a sorting algorithm, the corresponding entries in
!   the array ix will be swapped.
!   After execution, the array ix will hold a history of the transformations
!   If the array ix is pre-loaded with 1,2,3,4,..., then the sorted 
!   order of a will be a(ix(1)), a(ix(2)), a(ix(3)),....
!
  REAL(DP),INTENT(IN),DIMENSION(:):: a
  INTEGER,INTENT(IN OUT),DIMENSION(:):: ix

  INTEGER:: n
  INTEGER:: i,j,k, ir, ixsave
  REAL(DP):: asave
!----------------------------------------------------------------------------
  n=SIZE(a)   !   ix should also be as large
  IF (n <= 1) RETURN

  k=(n/2)+1
  ir=n

  DO
    IF (k > 1) THEN
      k=k-1
      asave=a(k)
      ixsave=ix(k)
    ELSE
      asave=a(ir)
      ixsave=ix(ir)
      ix(ir)=ix(1)
      ir=ir-1
      IF (ir == 1) THEN
        ix(1) = ixsave
        EXIT            ! this is the only way out of here
      END IF
    END IF

    i = k
    j = k+k

    DO
      IF ( j > ir) EXIT   ! gets you out of this loop
      IF (j < ir) THEN
        IF ( a(j) < a(j+1) ) j=j+1
      END IF

      IF ( asave < a(j) ) THEN
        ix(i) = ix(j)
        i = j
        j = j+j
      ELSE
        j = ir+1
      END IF
    END DO

    ix(i) = ixsave
  END DO                                  ! only way out is thru RETURN above
  RETURN
END Subroutine IndexDouble   ! ----------------------------------------------


END Module IntUtil   ! ======================================================


