INCLUDE 'vmaco.f90'
!+
PROGRAM ROSEN
! ------------------------------------------------------------------------------
! PURPOSE - This is a simple driver program for Vmaco. Contained below is the
!  Rosenbrock curved valley test case. It has 2 control variables and no 
!  constraints.
USE VmacoProcedures
IMPLICIT NONE

!D0   IDENTIFICATION :
!D0         PROGRAMMER -  J.D. FRICK, MDAC/HOUSTON
!D0         VERSION    -  SEPTEMBER 1986

  INTEGER,PARAMETER:: MSIZE = 4 ! size of arrays pertaining to the constraints
  INTEGER,PARAMETER:: NSIZE = 4 ! corresponds to the maximum dimension of the 
                                ! innermost array size from the cnorml matrix
                                ! (i.e. max control variable size). This allows 
                                ! correct communication between the cnorml 
                                ! and cn matrix.

  REAL,PARAMETER:: ACC = 1.0E-7 ! controls the final accuracy for convergence. 
    ! Convergence occurs when the value for the change in the objective function
    ! plus suitably weighted multipliers of the constraint function differs from
    ! its previous value by at most ACC.

  REAL,DIMENSION(NSIZE):: c ! vector of constraint variables.
  REAL,DIMENSION(NSIZE,MSIZE):: cnorml ! matrix of constraint normals
  REAL:: f   ! value of the objective function (optimized variable)
  REAL,DIMENSION(NSIZE):: g  ! gradient of the objective function

  INTEGER:: iprint = 1  !  flag to indicate amount of printout
                        !   = 0 no printout
                        !   = 1 values of x, c, and f are printed
                        !   = 2 debug printout

  INTEGER:: istat  ! index to indicate the status of the optimization algorithm.
                   !  =-1 initial pass through vmaco
                   !  = 0 continue algorithm calculations
                   !  = 1 converge within required acurracy
                   !  > 1 error encountered

  INTEGER:: m ! total number of constraints (i.e. constraint equations)
  INTEGER:: maxfun = 100 ! maximum number of times Vmaco can be called
  INTEGER:: meq = 0 ! total number of equality constraints
  INTEGER:: n  ! number of control variables
  REAL,DIMENSION(NSIZE):: x ! vector of control variables
!-------------------------------------------------------------------------------
  istat=-1  ! always start Vmaco with istat set to -1
  n=2  ! the number of control variables
  x(1:n) = (/ -1.0, 1.0 /)   ! initial guess for the control variables
  m=0  ! the number of constraints
  meq=0 ! the number of equality constraints
!  iprint=2 ! un-comment this line to see the algorithm in action

  DO   ! this is the computation loop...
    f = 100*(x(2)-x(1)**2)**2 + (1-x(1))**2   ! objective function
    g(1) = -400*x(1)*(x(2)-x(1)**2) - 2*(1-x(1)) ! gradient of objective function
    g(2) = 200*(x(2)-x(1)**2)

    CALL Vmaco(n,m,meq,x,f,g,c,cnorml,MAXFUN,ACC,istat,iprint)
    IF (istat /= 0) EXIT
  END DO
  WRITE(*,*) 'Program stop with istat= ', istat
  STOP
END Program Rosen   ! ==========================================================
