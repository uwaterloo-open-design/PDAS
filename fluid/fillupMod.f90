! I find that in the programs I develop, one of the most common procedures
! is the filling of an array, where I know the first and last values and it
! is necessary to fill in all of the intermediate values. I have spent many 
! years helping people find elusive errors in their computer programs, and I
! am convinced that more errors are made in filling arrays than any other
! single process (well, maybe writing format statements). This is usually 
! because it is easy to confuse the number of points in the array with the
! number of intervals. It seems that people constantly use n when they should
! use n-1 and vice versa. I have composed the subroutine FillArray to make
! this process a little closer to foolproof. 
! Even those who get their n's and n-1's right usually code this procedure by
! evaluating the distance between points and compute point 2 by adding this
! value to point 1. Then add it to the newly computed point 2 to get point
! 3, and so on, thereby accumulating additional roundoff error at each point.
! The coding given here will have no more roundoff at the end of the loop
! as at the first.
! FillArray is a Fortran90 creation and it has several features that set it
! off from similar procedures written in earlier versions of Fortran:
! 1) The size of the array to be filled is NOT an argument to the subroutine.
!    It is determined by FillArray using the SIZE statement, thereby 
!    eliminating the possibility of overrunning the allocated length of the
!    array.
! 2) Usually, I want uniform spacing of the intermediate points. But, there
!    are times when other spacings are desirable. These are selected with an
!    integer called spacingCode. This argument to the subroutine is coded as
!    an OPTIONAL argument. SpacingCode=2 gives denser spacing near the
!    endpoints and sparser near the middle; =3 is dense at the first endpoint
!    and sparse at the last; =4 is the reverse of 3. Anything but 2,3, or 4
!    gives uniform spacing. And, no argument at all gives uniform spacing.
! 3) Sometimes we work in single precision and other times in double 
!    precision. There are separate versions of FillArray in the module for
!    the different precisions. By the use of MODULE INTERFACE, the user may
!    use the subroutine by the same name, regardless of precision. There is
!    also a mixed interface, where the constants are given in single, while
!    the array is computed as double. Most of my codes use statements like
!       CALL FillArray(0.0, 1.0, a)
!    where a is a double precision array. So, I don't have to write
!       CALL FillArray(0.0D0, 1.0D0, a)
!    or something similar.
! 4) All Fortran programmers learn to use implied do-loops in READ and WRITE
!    statements. Many don't know that an implied do-loop can appear on the
!    right-hand side of the equal sign in an assignment. It is doubtful 
!    whether this construction would lead to any performance improvement, but
!    it does allow you to code the loop in one line rather than three.
! 5) I know the coding is fully compliant with Fortran 90, because it 
!    compiles and executes properly with Essential Lahey Fortran. I write all
!    my new code with Elf90, just to be sure I am not using a language
!    extention that may not work everywhere.

!+
MODULE FillUp
! ---------------------------------------------------------------------------
! PURPOSE - Hold the subroutines FillArraySingle, FillArrayDouble, and
!    FillArrayMixed.
!  Make all available with the name FillArray through a module interface.
! AUTHOR - Ralph L. Carmichael, Public Domain Aeronautical Software
!           ( ralph@pdas.com )   (  http://www.pdas.com  )

IMPLICIT NONE
  
  INTERFACE FillArray
    MODULE PROCEDURE FillArraySingle,FillArrayDouble,FillArrayMixed
  END INTERFACE  

  INTEGER,PRIVATE,PARAMETER:: SP = SELECTED_REAL_KIND(6,30)
  INTEGER,PRIVATE,PARAMETER:: DP = SELECTED_REAL_KIND(10,50)
!----------------------------------------------------------------------------

CONTAINS

!+
SUBROUTINE FillArraySingle(startValue,endValue, array, spacingCode)
! ---------------------------------------------------------------------------
! PURPOSE - fill an array from start to end. The intermediate points are
!    computed according to various spacing rules. Single Precision Version.

  REAL(SP),INTENT(IN):: startValue,endValue
  REAL(SP),INTENT(OUT),DIMENSION(:):: array
  INTEGER,INTENT(IN),OPTIONAL:: spacingCode
                              ! =2 full cosine
                              ! =3 half cosine
                              ! =4 half sine
                              ! anything else (or nothing)= uniform spacing
  INTEGER:: k,n
  REAL(SP),PARAMETER:: HALF = 0.5_SP, ONE = 1.0_SP
  REAL(SP),PARAMETER:: PI=3.14159265_SP, HALFPI=HALF*PI
  REAL(SP),ALLOCATABLE,DIMENSION(:):: temp
!----------------------------------------------------------------------------
  n=SIZE(array)
  IF (n <= 0) RETURN

  array(n)=endValue
  array(1)=startValue
  IF (n <= 2) RETURN

  ALLOCATE(temp(n-2))
  temp= (/ (REAL(k), k=1,n-2) /)         ! implied do-loop on right-hand side
  temp=temp/REAL(n-1)                                 ! very little round-off

  IF (Present(spacingCode)) THEN
    SELECT CASE(spacingCode)
      CASE (2)
        temp=HALF*(ONE-COS(PI*temp))       ! full cosine, dense near both ends
      CASE (3)
        temp=ONE-COS(HALFPI*temp)             ! half cosine, dense near start
      CASE (4)
        temp=SIN(HALFPI*temp)                     ! half sine, dense near end
    END SELECT
  END IF

  array(2:n-1)=startValue + (endValue-startValue)*temp

  DEALLOCATE(temp)
  RETURN
END Subroutine FillArraySingle   ! ------------------------------------------------

!+
SUBROUTINE FillArrayDouble(startValue,endValue, array, spacingCode)
! ---------------------------------------------------------------------------
! PURPOSE - fill an array from start to end. The intermediate points are
!    computed according to various spacing rules. Double Precision Version.

  REAL(DP),INTENT(IN):: startValue,endValue
  REAL(DP),INTENT(OUT),DIMENSION(:):: array
  INTEGER,INTENT(IN),OPTIONAL:: spacingCode
                              ! =2 full cosine
                              ! =3 half cosine
                              ! =4 half sine
                              ! anything else (or nothing)= uniform spacing
  INTEGER:: k,n
  REAL(DP),PARAMETER:: HALF = 0.5_DP, ONE = 1.0_DP
  REAL(DP),PARAMETER:: PI = 3.141592653589793238462643_DP, HALFPI=HALF*PI
  REAL(DP),ALLOCATABLE,DIMENSION(:):: temp
!----------------------------------------------------------------------------
  n=SIZE(array)
  IF (n <= 0) RETURN

  array(n)=endValue
  array(1)=startValue
  IF (n <= 2) RETURN

  ALLOCATE(temp(n-2))
  temp= (/ (REAL(k), k=1,n-2) /)         ! implied do-loop on right-hand side
  temp=temp/REAL(n-1)                                 ! very little round-off

  IF (Present(spacingCode)) THEN
    SELECT CASE(spacingCode)
      CASE (2)
        temp=HALF*(ONE-COS(PI*temp))       ! full cosine, dense near both ends
      CASE (3)
        temp=ONE-COS(HALFPI*temp)             ! half cosine, dense near start
      CASE (4)
        temp=SIN(HALFPI*temp)                     ! half sine, dense near end
    END SELECT
  END IF

  array(2:n-1)=startValue + (endValue-startValue)*temp

  DEALLOCATE(temp)
  RETURN
END Subroutine FillArrayDouble   ! ------------------------------------------

!+
SUBROUTINE FillArrayMixed(startValue,endValue, array, spacingCode)
! ---------------------------------------------------------------------------
! PURPOSE - As above, but startValue and endValue are single precision
!  while the target array is double precision.
  REAL(SP),INTENT(IN):: startValue,endValue
  REAL(DP),INTENT(OUT),DIMENSION(:):: array
  INTEGER,INTENT(IN),OPTIONAL:: spacingCode

  REAL(DP):: startDouble,endDouble
!----------------------------------------------------------------------------
  startDouble=startValue    ! converts SP to DP
  endDouble=endValue
  IF (Present(spacingCode)) THEN
    CALL FillArrayDouble(startDouble,endDouble, array, spacingCode)
  ELSE
    CALL FillArrayDouble(startDouble,endDouble, array)
  END IF
  RETURN
END Subroutine FillArrayMixed   ! -------------------------------------------

END Module FillUp   ! =======================================================
