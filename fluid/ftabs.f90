!+
! PROGRAM FluidTables
! ---------------------------------------------------------------------------
! PURPOSE - Prepare new subroutines for the Fluid package. Each subroutine
!  will load the common block /FLUDPC/ with appropriate data
!
! AUTHOR -  Ralph L. Carmichael, Public Domain Aeronautical Software

! REVISION HISTORY
!   DATE  VERS PERSON  STATEMENT OF CHANGES
! 07Jun10  0.5   RLC   Original coding (just O2)
! 10Jun10  0.55  RLC   Added dryair
! 16Jun10  0.6   RLC   Added all remaining fluids



!+
MODULE CommonFludpc
! ------------------------------------------------------------------------------
! PURPOSE - Emulate the old common block /FLUDPC/

!      DIMENSION G(29,4),TT(27),DD(30),GG(27,30,8),SAT(40,8),C(8),M(4) 
!       INTEGER ENTRY,ERROR 
!       DIMENSION PROPS(8) 
!       LOGICAL VAPOR 
!!!       COMMON /FLUDPC/ G,TT,DD,GG,SAT,C,M,NG,N1,N2,N3,NS 
  REAL,DIMENSION(29,4):: g
  REAL,DIMENSION(27):: tt
  REAL,DIMENSION(30):: dd
  REAL,DIMENSION(27,30,8):: gg
  REAL,DIMENSION(40,8):: sat
  REAL,DIMENSION(8):: c
  INTEGER,DIMENSION(4):: m
  INTEGER:: ng
  INTEGER:: n1
  INTEGER:: n2
  INTEGER:: n3
  INTEGER:: ns

!-------------------------------------------------------------------------------

END Module CommonFludpc   !=====================================================

!+
MODULE FluidTablesProcedures
! ------------------------------------------------------------------------------
! PURPOSE - Make the new subroutines
IMPLICIT NONE
!-------------------------------------------------------------------------------

  CHARACTER(LEN=*),PARAMETER:: VERSION = "Version 0.55 (10 June 2010)"
  INTEGER,PARAMETER:: SUB=2, DBG=3
  CHARACTER(LEN=15):: dateTimeStr
  CHARACTER(LEN=*),PARAMETER:: LINE1 = &
    "! ------------------------------------------------------------------------------"
  CHARACTER(LEN=*),PARAMETER:: LINE2 = &
    "!-------------------------------------------------------------------------------"
CONTAINS
 
INCLUDE 'TEMPo2.f90'
INCLUDE 'TEMPdryair.f90'
INCLUDE 'TEMPar.f90'
INCLUDE 'TEMPn2.f90'
INCLUDE 'TEMPph2.f90'
INCLUDE 'TEMPf2.f90'
INCLUDE 'TEMPsteam.f90'
INCLUDE 'TEMPco2.f90'
INCLUDE 'TEMPch4.f90'

!+
FUNCTION GetDateTimeStr() RESULT (s)
!   ----------------------------------------------------------------------------
! PURPOSE - Return a string with the current date and time
  IMPLICIT NONE
  CHARACTER(LEN=*),PARAMETER:: MONTH='JanFebMarAprMayJunJulAugSepOctNovDec'
  CHARACTER(LEN=*),PARAMETER:: FMT='(I2.2,A1,I2.2,I3,A3,I4)'
  CHARACTER(LEN=15):: s
  INTEGER,DIMENSION(8):: v
  INTRINSIC:: DATE_AND_TIME
!-------------------------------------------------------------------------------
  CALL DATE_AND_TIME(VALUES=v)
  WRITE(s,FMT) v(5), ':', v(6), v(3), MONTH(3*v(2)-2:3*v(2)), v(1)
  RETURN
END FUNCTION GetDateTimeStr   ! ================================================

!+
SUBROUTINE MakeSubroutineOxygen()
! ------------------------------------------------------------------------------
! PURPOSE - Make the new oxygen tables
!-------------------------------------------------------------------------------
  WRITE(DBG,*) "Beginning O2"
  OPEN(UNIT=SUB,FILE='newO2.f90',STATUS='REPLACE',ACTION='WRITE')
  WRITE(SUB,'(A)') "!+"
  WRITE(SUB,'(A)') "SUBROUTINE O2()"
  WRITE(SUB,'(A)') LINE1
  WRITE(SUB,'(A)') "! PURPOSE - Make tables for oxygen. Created "//dateTimeStr

  WRITE(SUB,*) 'USE CommonFludPC,ONLY: ng,n1,n2,n3,ns,g,tt,dd,gg,sat,c'

  CALL TempO2(SUB)

  WRITE(SUB,'(A)') 'RETURN'
  WRITE(SUB,'(A)') 'END Subroutine O2   ! -----------------------------------'
  CLOSE(UNIT=SUB)

  WRITE(*,*) "Subroutine Oxygen completed."
  WRITE(DBG,*) "Ending Oxygen"
  RETURN
END Subroutine MakeSubroutineOxygen   ! --------------------------------------------

!+
SUBROUTINE MakeSubroutineDryAir()
! ------------------------------------------------------------------------------
! PURPOSE - Make the new dry air tables
!-------------------------------------------------------------------------------
  WRITE(DBG,*) "Beginning dryair"
  OPEN(UNIT=SUB,FILE='newdryair.f90',STATUS='REPLACE',ACTION='WRITE')
  WRITE(SUB,'(A)') "!+"
  WRITE(SUB,'(A)') "SUBROUTINE DRYAIR()"
  WRITE(SUB,'(A)') LINE1
  WRITE(SUB,'(A)') "! PURPOSE - Make tables for dry air. Created "//dateTimeStr

  WRITE(SUB,*) 'USE CommonFludPC,ONLY: ng,n1,n2,n3,ns,g,tt,dd,gg,sat,c'

  CALL TempDryAir(SUB)

  WRITE(SUB,'(A)') 'RETURN'
  WRITE(SUB,'(A)') 'END Subroutine DRYAIR   ! -----------------------------------'
  CLOSE(UNIT=SUB)

  WRITE(*,*) "Subroutine dryair completed."
  WRITE(DBG,*) "Ending dryair"
  RETURN
END Subroutine MakeSubroutineDryAir   ! --------------------------------------------


!+
SUBROUTINE MakeSubroutineArgon()
! ------------------------------------------------------------------------------
! PURPOSE - Make the new argon tables
!-------------------------------------------------------------------------------
  WRITE(DBG,*) "Beginning argon"
  OPEN(UNIT=SUB,FILE='newar.f90',STATUS='REPLACE',ACTION='WRITE')
  WRITE(SUB,'(A)') "!+"
  WRITE(SUB,'(A)') "SUBROUTINE AR()"
  WRITE(SUB,'(A)') LINE1
  WRITE(SUB,'(A)') "! PURPOSE - Make tables for argon. Created "//dateTimeStr

  WRITE(SUB,*) 'USE CommonFludPC,ONLY: ng,n1,n2,n3,ns,g,tt,dd,gg,sat,c'

  CALL TempAR(SUB)

  WRITE(SUB,'(A)') 'RETURN'
  WRITE(SUB,'(A)') 'END Subroutine AR   ! -----------------------------------'
  CLOSE(UNIT=SUB)

  WRITE(*,*) "Subroutine ar completed."
  WRITE(DBG,*) "Ending argon"
  RETURN
END Subroutine MakeSubroutineArgon   ! --------------------------------------------


!+
SUBROUTINE MakeSubroutineNitrogen()
! ------------------------------------------------------------------------------
! PURPOSE - Make the new nitrogen tables
!-------------------------------------------------------------------------------
  WRITE(DBG,*) "Beginning nitrogen"
  OPEN(UNIT=SUB,FILE='newn2.f90',STATUS='REPLACE',ACTION='WRITE')
  WRITE(SUB,'(A)') "!+"
  WRITE(SUB,'(A)') "SUBROUTINE NITRO()"
  WRITE(SUB,'(A)') LINE1
  WRITE(SUB,'(A)') "! PURPOSE - Make tables for nitrogen. Created "//dateTimeStr

  WRITE(SUB,*) 'USE CommonFludPC,ONLY: ng,n1,n2,n3,ns,g,tt,dd,gg,sat,c'

  CALL TempN2(SUB)

  WRITE(SUB,'(A)') 'RETURN'
  WRITE(SUB,'(A)') 'END Subroutine NITRO   ! -----------------------------------'
  CLOSE(UNIT=SUB)

  WRITE(*,*) "Subroutine nitrogen completed."
  WRITE(DBG,*) "Ending nitrogen"
  RETURN
END Subroutine MakeSubroutineNitrogen   ! --------------------------------------------

!+
SUBROUTINE MakeSubroutineParahydrogen()
! ------------------------------------------------------------------------------
! PURPOSE - Make the new parahydrogen tables
!-------------------------------------------------------------------------------
  WRITE(DBG,*) "Beginning parahydrogen"
  OPEN(UNIT=SUB,FILE='newph2.f90',STATUS='REPLACE',ACTION='WRITE')
  WRITE(SUB,'(A)') "!+"
  WRITE(SUB,'(A)') "SUBROUTINE PH2()"
  WRITE(SUB,'(A)') LINE1
  WRITE(SUB,'(A)') "! PURPOSE - Make tables for parahydrogen. Created "//dateTimeStr

  WRITE(SUB,*) 'USE CommonFludPC,ONLY: ng,n1,n2,n3,ns,g,tt,dd,gg,sat,c'

  CALL TempPH2(SUB)

  WRITE(SUB,'(A)') 'RETURN'
  WRITE(SUB,'(A)') 'END Subroutine PH2   ! -----------------------------------'
  CLOSE(UNIT=SUB)

  WRITE(*,*) "Subroutine parahydrogen completed."
  WRITE(DBG,*) "Ending parahydrogen"
  RETURN
END Subroutine MakeSubroutineParahydrogen   ! --------------------------------------------

!+
SUBROUTINE MakeSubroutineFluorine()
! ------------------------------------------------------------------------------
! PURPOSE - Make the new fluorine tables
!-------------------------------------------------------------------------------
  WRITE(DBG,*) "Beginning fluorine"
  OPEN(UNIT=SUB,FILE='newf2.f90',STATUS='REPLACE',ACTION='WRITE')
  WRITE(SUB,'(A)') "!+"
  WRITE(SUB,'(A)') "SUBROUTINE F2()"
  WRITE(SUB,'(A)') LINE1
  WRITE(SUB,'(A)') "! PURPOSE - Make tables for fluorine. Created "//dateTimeStr

  WRITE(SUB,*) 'USE CommonFludPC,ONLY: ng,n1,n2,n3,ns,g,tt,dd,gg,sat,c'

  CALL TempF2(SUB)

  WRITE(SUB,'(A)') 'RETURN'
  WRITE(SUB,'(A)') 'END Subroutine F2   ! -----------------------------------'
  CLOSE(UNIT=SUB)

  WRITE(*,*) "Subroutine fluorine completed."
  WRITE(DBG,*) "Ending fluorine"
  RETURN
END Subroutine MakeSubroutineFluorine   ! --------------------------------------------

!+
SUBROUTINE MakeSubroutineSteam()
! ------------------------------------------------------------------------------
! PURPOSE - Make the new steam tables
!-------------------------------------------------------------------------------
  WRITE(DBG,*) "Beginning steam"
  OPEN(UNIT=SUB,FILE='newsteam.f90',STATUS='REPLACE',ACTION='WRITE')
  WRITE(SUB,'(A)') "!+"
  WRITE(SUB,'(A)') "SUBROUTINE STEAM()"
  WRITE(SUB,'(A)') LINE1
  WRITE(SUB,'(A)') "! PURPOSE - Make tables for steam. Created "//dateTimeStr

  WRITE(SUB,*) 'USE CommonFludPC,ONLY: ng,n1,n2,n3,ns,g,tt,dd,gg,sat,c'

  CALL TempSTEAM(SUB)

  WRITE(SUB,'(A)') 'RETURN'
  WRITE(SUB,'(A)') 'END Subroutine Steam   ! -----------------------------------'
  CLOSE(UNIT=SUB)

  WRITE(*,*) "Subroutine steam completed."
  WRITE(DBG,*) "Ending steam"
  RETURN
END Subroutine MakeSubroutineSteam   ! --------------------------------------------

!+
SUBROUTINE MakeSubroutineCarbonDioxide()
! ------------------------------------------------------------------------------
! PURPOSE - Make the new CO2 tables
!-------------------------------------------------------------------------------
  WRITE(DBG,*) "Beginning CarbonDioxide"
  OPEN(UNIT=SUB,FILE='newco2.f90',STATUS='REPLACE',ACTION='WRITE')
  WRITE(SUB,'(A)') "!+"
  WRITE(SUB,'(A)') "SUBROUTINE CO2()"
  WRITE(SUB,'(A)') LINE1
  WRITE(SUB,'(A)') "! PURPOSE - Make tables for carbon dioxide. Created "//dateTimeStr

  WRITE(SUB,*) 'USE CommonFludPC,ONLY: ng,n1,n2,n3,ns,g,tt,dd,gg,sat,c'

  CALL TempCO2(SUB)

  WRITE(SUB,'(A)') 'RETURN'
  WRITE(SUB,'(A)') 'END Subroutine CO2   ! -----------------------------------'
  CLOSE(UNIT=SUB)

  WRITE(*,*) "Subroutine CarbonDioxide completed."
  WRITE(DBG,*) "Ending CarbonDioxide"
  RETURN
END Subroutine MakeSubroutineCarbonDioxide   ! --------------------------------------------

!+
SUBROUTINE MakeSubroutineMethane()
! ------------------------------------------------------------------------------
! PURPOSE - Make the new methane tables
!-------------------------------------------------------------------------------
  WRITE(DBG,*) "Beginning Methane"
  OPEN(UNIT=SUB,FILE='newch4.f90',STATUS='REPLACE',ACTION='WRITE')
  WRITE(SUB,'(A)') "!+"
  WRITE(SUB,'(A)') "SUBROUTINE CH4()"
  WRITE(SUB,'(A)') LINE1
  WRITE(SUB,'(A)') "! PURPOSE - Make tables for methane. Created "//dateTimeStr

  WRITE(SUB,*) 'USE CommonFludPC,ONLY: ng,n1,n2,n3,ns,g,tt,dd,gg,sat,c'

  CALL TempCH4(SUB)

  WRITE(SUB,'(A)') 'RETURN'
  WRITE(SUB,'(A)') 'END Subroutine CH4   ! -----------------------------------'
  CLOSE(UNIT=SUB)

  WRITE(*,*) "Subroutine Methane completed."
  WRITE(DBG,*) "Ending Methane"
  RETURN
END Subroutine MakeSubroutineMethane   ! --------------------------------------------


!+
SUBROUTINE DefineRealOneDimension(efu,a,name)
! ------------------------------------------------------------------------------
! PURPOSE - Write the Fortran declaration for a one-dimensional array defined
!  as a parameter. Each real in a is converted into a character variable of
!  length 11 stored in the array named fields. fields is allocated to have the
!  same rank as a.
  INTEGER,INTENT(IN):: efu
  REAL,INTENT(IN),DIMENSION(:):: a
  CHARACTER(LEN=*),INTENT(IN):: name

  INTEGER:: k,n,n1,n2
  INTRINSIC:: SIZE
  CHARACTER(LEN=11),ALLOCATABLE,DIMENSION(:):: fields
  CHARACTER(LEN=*),PARAMETER:: FMT = '(ES10.4E1,",")'
!-------------------------------------------------------------------------------
  n=SIZE(a)
  ALLOCATE(fields(n))   ! the same rank as a

  DO k=1,n
    WRITE(fields(k),FMT) a(k)
  END DO
  fields(n)(11:11)=' '

  WRITE(efu,'(A,I0,A,A,A)') &
    ' REAL,PARAMETER,DIMENSION(', n,'):: ', name, ' = (/ &'

  n1=-7
  DO k=1,n/8
    n1=n1+8      ! 1,9,17,25,33,...
    n2=n1+7      ! 8,16,24,32,40,...
    WRITE(efu,'(2X,8A11,A)') fields(n1:n2), ' &'
  END DO

! see if there are a few left over...
  n1=n1+8
  n2=n
  IF (n2>n1) WRITE(efu,'(T92,"&",T3,8A11)') fields(n1:n2)

  WRITE(efu,*) '    /)'

  DEALLOCATE(fields)
  RETURN
END Subroutine DefineRealOneDimension   ! --------------------------------------


!+
SUBROUTINE DefineRealTwoDimension(efu,a,name)
! ------------------------------------------------------------------------------
! PURPOSE - Write the Fortran declaration for a two-dimensional array defined
!  as a parameter. Each real in a is converted into a character variable of
!  length 11 stored in the array named fields. fields is allocated to have the
!  same rank as a.
  INTEGER,INTENT(IN):: efu
  REAL,INTENT(IN),DIMENSION(:,:):: a
  CHARACTER(LEN=*),INTENT(IN):: name

  INTEGER:: j,k,n,n1,n2,nrows,ncols
  INTRINSIC:: SIZE
  CHARACTER(LEN=11),ALLOCATABLE,DIMENSION(:,:):: fields
  CHARACTER(LEN=*),PARAMETER:: FMT = '(ES10.4E1,",")'
!-------------------------------------------------------------------------------
  nrows=SIZE(a,1)
  ncols=SIZE(a,2)
  ALLOCATE(fields(nrows,ncols))   ! the same rank as a
  WRITE(DBG,*) "fields allocated", nrows,ncols

  DO k=1,ncols
    DO j=1,nrows
      WRITE(fields(j,k),FMT) a(j,k)
    END DO
  END DO
  fields(nrows,ncols)(11:11)=' '

  DO k=1,ncols
    DO j=1,nrows
      WRITE(DBG,*) fields(j,k)
    END DO
  END DO

  WRITE(efu,'(A,I0,A,I0,A,A,A)') &
    ' REAL,PARAMETER,DIMENSION(', nrows,",",ncols,'):: ', name, ' = RESHAPE( (/ &'

  DO k=1,ncols
    n1=-7
    DO j=1,nrows/8
      n1=n1+8      ! 1,9,17,25,33,...
      n2=n1+7      ! 8,16,24,32,40,...
      WRITE(efu,'(2X,8A11,A)') fields(n1:n2,k), ' &'
      write(dbg,*) 'writing fields', n1,n2,k
    END DO

! see if there are a few left over...
    n1=n1+8
    n2=nrows
    IF (n2>n1) WRITE(efu,'(T92,"&",T3,8A11)') fields(n1:n2,k)
  END DO

  WRITE(efu,'(A,I0,A,I0,A)') '    /), (/',nrows, ',', ncols, '/)  )'

  DEALLOCATE(fields)
  RETURN
END Subroutine DefineRealTwoDimension   ! --------------------------------------

END Module FluidTablesProcedures   ! ===========================================


!+
PROGRAM FluidTables
! ------------------------------------------------------------------------------
USE FluidTablesProcedures
IMPLICIT NONE



!-------------------------------------------------------------------------------
  dateTimeStr=GetDateTimeStr()
  WRITE(*,*) "ftabs - FluidTables. Make new gas subroutines for FLUID"
  WRITE(*,*) VERSION
  WRITE(*,*) "run date: "//dateTimeStr

  OPEN(UNIT=DBG,FILE='ftabs.dbg',STATUS='REPLACE',ACTION='WRITE')

  CALL MakeSubroutineOxygen()
  CALL MakeSubroutineDryAir()
  CALL MakeSubroutineArgon()
  CALL MakeSubroutineNitrogen()
  CALL MakeSubroutineParahydrogen()
  CALL MakeSubroutineFluorine()
  CALL MakeSubroutineSteam()
  CALL MakeSubroutineCarbonDioxide()
  CALL MakeSubroutineMethane()

  WRITE(*,*) 'All subroutines created. Normal termination.'

  STOP
END Program FluidTables   ! ====================================================
