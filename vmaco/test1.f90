INCLUDE 'vmaco.f90'
!+
PROGRAM DRV3X2
! ------------------------------------------------------------------------------
! PURPOSE - This is a simple driver program for Vmaco. Contained below is a 
!  constrained optimization test case with 3 control variables and 2 
!  constraints. It may be noted that one of the constraints is an equality 
!  and one is an inequality type.  [ m=2, but meq=1 ]

USE VmacoProcedures
IMPLICIT NONE

!D0
!D0   IDENTIFICATION :
!D0         PROGRAMMER -  J.D. FRICK, MDAC/HOUSTON
!D0         VERSION    -  SEPTEMBER 1986

  INTEGER,PARAMETER:: MSIZE = 4 ! size of arrays pertaining to the constraints
  INTEGER,PARAMETER:: NSIZE = 4 ! corresponds to the maximum dimension of the 
                                ! innermost array size from the cnorml matrix
                                ! (i.e. max control variable size). This allows 
                                ! correct communication between the cnorml 
                                ! and cn matrix.

  REAL:: ACC = 1.0E-7 ! controls the final accuracy for convergence. 
                      ! Convergence occurs when the value for the change in
                      ! the objective function plus suitably weighted multipliers
                      ! of the constraint function differs from its previous 
                      ! value by at most acc.

  REAL,DIMENSION(NSIZE):: C ! vector of constraint variables.
  REAL,DIMENSION(NSIZE,MSIZE):: CNORML ! matrix of constraint normals

  REAL:: F   ! value of the objective function (optimized variable)
  REAL,DIMENSION(NSIZE):: G  ! gradient of the objective function

  INTEGER:: iprint = 1  !  flag to indicate amount of printout
                        !   = 0 no printout
                        !   = 1 values of x, c, and f are printed
                        !   = 2 debug printout

  INTEGER:: istat = -1  ! index to indicate the status of the optimization
                        !  algorithm.
                        !  =-1 initial pass through vmaco
                        !  = 0 continue algorithm calculations
                        !  = 1 converge within required acurracy
                        !  > 1 error encountered

  INTEGER:: m  ! total number of constraints (i.e. constraint equations)
  INTEGER,PARAMETER:: MAXFUN = 20 ! maximum number of times Vmaco can be called


  INTEGER:: meq = 1 ! total number of equality constraints
  INTEGER:: n = 3  ! number of control variables
  REAL,DIMENSION(NSIZE):: X ! vector of control variables
!-------------------------------------------------------------------------------
  istat=-1  ! always start Vmaco with istat set to -1
  n=3
  x(1:n) = (/ 1.0, 2.0, 3.0 /) ! initial guess for the control variables
  m=2   ! total number of constraints
  meq=1 ! number of equality constraints
!  iprint=2 ! un-comment this line to see the algorithm in action

  DO   ! this is the computation loop...
    f = x(1)**2 + x(2)**2 + x(3) ! the objective function is calculated
    g(1) = 2.*x(1)  ! the gradient of the objective function are calculated
    g(2) = 2.*x(2)
    g(3) = 1.0

    c(1) = x(1)*x(2) - x(3)  ! the constraints are calculated
    cnorml(1,1) = x(2)       !  and their first partials 
    cnorml(2,1) = x(1)
    cnorml(3,1) = -1.0

    c(2) = x(3) - 1.
    cnorml(1,2) = 0.
    cnorml(2,2) = 0.
    cnorml(3,2) = 1.

    CALL Vmaco (n,m,meq,x,f,g,c,cnorml,MAXFUN,ACC,istat,iprint)
    IF (istat /= 0) EXIT
  END DO

  WRITE(*,*) 'Program stop with istat= ', istat
  STOP
END Program Drv3x2   ! =========================================================
