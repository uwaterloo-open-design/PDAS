! PROGRAM AirTables
! ------------------------------------------------------------------------------
! PURPOSE - Attempt to reproduce the tables from NBS Circular 564 for dry air
!  The tables show thermodynamic and transport properties of dry air as a
!  function of temperature and pressure.

! AUTHORS - Theodore E. Fessler, NASA Lewis Research Center
!           Lt. Mark D. Klem & Margaret P. Proctor, NASA?
!           Ralph L. Carmichael, Public Domain Aeronautical Software

! REVISION HISTORY
!   DATE  VERS PERSON  STATEMENT OF CHANGES
! 18Jun10  0.1   RLC   Original coding



!-------------------------------------------------------------------------------

INCLUDE 'futil.f90'

!+
MODULE AirTablesProcedures
! ------------------------------------------------------------------------------
! PURPOSE - Collect the procedures and global variables used by program AirTables.
USE FluidProperties
USE FluidProcedures
IMPLICIT NONE

  CHARACTER(LEN=*),PARAMETER:: VERSION = ' Version 0.1 (18 June 2010)'
  INTEGER,PARAMETER:: OUT=2, HTML=4   ! DBG is defined as 3 in FluidProcedures

  INTEGER,PARAMETER:: NP = 5   ! length of P_TABLE
  INTEGER,PARAMETER:: NT = 4   ! length of T_TABLE
  REAL,PARAMETER,DIMENSION(NP):: P_TABLE = (/ 0.1, 1.0, 10.0, 40.0, 100.0 /) ! atmospheres
  REAL,PARAMETER,DIMENSION(NT):: T_TABLE = (/ 200.0, 300.0, 500.0, 1000.0 /)  ! kelvins

  REAL,DIMENSION(NT,NP):: compressibility
  REAL,DIMENSION(NT,NP):: sigma   ! density/critical density
  REAL,DIMENSION(NT,NP):: enthalpy  ! H / (R Tc)
  REAL,DIMENSION(NT,NP):: entropy   ! S / R
  REAL,DIMENSION(NT,NP):: sonic
  REAL,DIMENSION(NT,NP):: specificHeat  ! Cv / R
  REAL,DIMENSION(NT,NP):: gamma


!-------------------------------------------------------------------------------

CONTAINS

!+
SUBROUTINE BuildTables()
! ------------------------------------------------------------------------------
! PURPOSE - Compute the compressibility factor at each value of the 
!  two-dimensional space defined by P_TABLE X T_TABLE. This only fills the
!  array, but does not print it

  INTEGER:: j,k
  REAL:: pressure   ! in megapascals
  REAL:: temperature ! in kelvins
  REAL:: rho  ! in grams per cubic centimeter
  REAL,DIMENSION(8):: props
  INTEGER,PARAMETER:: NP = 8
  INTEGER,PARAMETER:: ENTRY = 3 ! input is temperature and pressure
  LOGICAL:: vapor
  INTEGER:: error
  INTEGER:: IGOTO = 1
!------------------------------------------------------------------------------
  DO k=1,NP
    DO j=1,NT
      pressure=0.101325*P_TABLE(k)
      temperature=T_TABLE(j)
      WRITE(DBG,*) '*** BuildTables, k,j=', k,j
      WRITE(DBG,*) 'Before Fluid, temp, pres, rho', temperature, pressure,rho
      CALL Fluid(temperature,pressure,rho,props,NP,ENTRY,vapor,error,IGOTO)
      WRITE(DBG,*) 'After Fluid, temp, pres, rho', temperature, pressure,rho
      compressibility(j,k)=props(1)
      sigma(j,k)=rho/0.001225
      entropy(j,k)=props(2)/c(1)
      enthalpy(j,k)=props(3)/(c(1)*c(2))
      specificHeat(j,k)=props(4)/c(1)
      gamma(j,k)=props(5)/props(4)
      sonic(j,k)=props(6)
    END DO
  END DO

  RETURN
END Subroutine BuildTables   ! -------------------------------------------------

!+
SUBROUTINE PrintTables()
! ------------------------------------------------------------------------------
! PURPOSE - Print all tables with formatting for monospace font.

  CHARACTER(LEN=*),PARAMETER:: FMT1 = '(5X,10F8.2)'
  CHARACTER(LEN=*),PARAMETER:: FMT2 = '(I6,10F8.4)'
  INTEGER:: j
!-------------------------------------------------------------------------------
  WRITE(OUT,*) 'MISC. CONSTANTS', c

  WRITE(OUT,*) '     COMPRESSIBILITY FACTOR FOR AIR'
  WRITE(OUT,FMT1) P_TABLE
  DO j=1,NT
    WRITE(OUT,FMT2) INT(T_TABLE(j)), compressibility(j,:)
  END DO

  WRITE(OUT,'(///A)') '    DENSITY (sigma) OF AIR'
  WRITE(OUT,FMT1) P_TABLE
  DO j=1,NT
    WRITE(OUT,FMT2) INT(T_TABLE(j)), sigma(j,:)
  END DO

  WRITE(OUT,'(///A)') '    ENTROPY (S/R) OF AIR'
  WRITE(OUT,FMT1) P_TABLE
  DO j=1,NT
    WRITE(OUT,FMT2) INT(T_TABLE(j)), entropy(j,:)
  END DO

  WRITE(OUT,'(///A)') '    ENTHALPY (H/(Tc R))  OF AIR'
  WRITE(OUT,FMT1) P_TABLE
  DO j=1,NT
    WRITE(OUT,FMT2) INT(T_TABLE(j)), enthalpy(j,:)
  END DO

  WRITE(OUT,'(///A)') '    GAMMA OF AIR'
  WRITE(OUT,FMT1) P_TABLE
  DO j=1,NT
    WRITE(OUT,FMT2) INT(T_TABLE(j)), gamma(j,:)
  END DO

  WRITE(OUT,'(///A)') '    SONIC VELOCITY OF AIR'
  WRITE(OUT,FMT1) P_TABLE
  DO j=1,NT
    WRITE(OUT,FMT2) INT(T_TABLE(j)), sonic(j,:)
  END DO

  WRITE(OUT,'(///A)') '    SPECIFIC HEAT (Cv/R) OF AIR'
  WRITE(OUT,FMT1) P_TABLE
  DO j=1,NT
    WRITE(OUT,FMT2) INT(T_TABLE(j)), specificHeat(j,:)
  END DO


  RETURN
END SUBROUTINE PrintTables   ! -------------------------------------------------

!+
SUBROUTINE HtmlTables()
! ------------------------------------------------------------------------------
! PURPOSE - Write all of the tables in HTML 5 format for viewing with a browser.

!-------------------------------------------------------------------------------
  WRITE(HTML,*) '<html>'

  WRITE(HTML,*) '</html>'
  RETURN
END Subroutine HtmlTables   ! --------------------------------------------------

END Module AirTablesProcedures   ! =============================================

!+
PROGRAM AirTables
! ------------------------------------------------------------------------------
USE AirTablesProcedures
IMPLICIT NONE


!-------------------------------------------------------------------------------
  WRITE(*,*) 'airtables - reproduce data in NBS Circular 564'
  WRITE(*,*) VERSION
  OPEN(UNIT=DBG,FILE='airtables.dbg',STATUS='REPLACE',ACTION='WRITE')
  OPEN(UNIT=OUT,FILE='airtables.out',STATUS='REPLACE',ACTION='WRITE')
  OPEN(UNIT=HTML,FILE='airtables.html',STATUS='REPLACE',ACTION='WRITE')

  CALL DryAir()

  CALL BuildTables()
  CALL PrintTables()
  CALL HtmlTables()

  WRITE(*,*) 'All tables have been completed.'
  STOP
END Program AirTables   ! ======================================================

