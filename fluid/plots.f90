!+
! PROGRAM FluidPlots
! ------------------------------------------------------------------------------
! PURPOSE - Make plots of the data tables used in FLUID.

! AUTHORS - Theodore E. Fessler, NASA Lewis Research Center
!           Lt. Mark D. Klem & Margaret P. Proctor, NASA?
!           Ralph L. Carmichael, Public Domain Aeronautical Software

! REVISION HISTORY
!   DATE  VERS PERSON  STATEMENT OF CHANGES
! 18Jun10  0.1   RLC   Original coding

INCLUDE 'figmods3.f90'

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
  INTEGER,DIMENSION(4):: mIndexValue
  INTEGER:: ng
  INTEGER:: n1
  INTEGER:: n2
  INTEGER:: n3
  INTEGER:: ns

!-------------------------------------------------------------------------------

END Module CommonFludpc   !=====================================================


!+
MODULE FluidPlotsProcedures
! ------------------------------------------------------------------------------
! PURPOSE - Collect the procedures and global variables used by FluidPlots
USE FigUtilities
USE CommonFludPC,ONLY: ng,n1,n2,n3,ns,g,tt,dd,gg,sat,c
IMPLICIT NONE

  CHARACTER(LEN=*),PARAMETER:: VERSION = '0.1 (20 June 2010)'
  CHARACTER(LEN=15):: dateTimeStr
  INTEGER,PARAMETER:: OUT=2, DBG=3, FIG=4, PS=7

!-------------------------------------------------------------------------------

CONTAINS

INCLUDE 'newo2.f90'

!+
SUBROUTINE PlotIdealGasFunctions(fluidName)
! ------------------------------------------------------------------------------
! PURPOSE - Make plots of entropy, internal energy, and Cv as functions of
!  temperature, either reduced temperature or in kelvins.
!  This procedure assumes that the data for the appropriate gas is loaded in 
!  module CommonFludPC
  CHARACTER(LEN=*),INTENT(IN):: fluidName

  REAL,DIMENSION(1000):: xplot,yplot

!-------------------------------------------------------------------------------

! plot entropy vs. reduced temperature

  xplot(1:ng)=g(1:ng,1)
  yplot(1:ng)=g(1:ng,2)

  WRITE(FIG,*) 'VIEWPORT 0.1 0.7  0.1 0.7'
  WRITE(FIG,*) 'WINDOW  0 10 10 30'
  CALL FigGrid(FIG, 10,4, 0,0, 'T/Tc','S/R',  5,12, 10,8)
  WRITE(FIG,*) 'SYMBOL 1'
  CALL FigData(FIG, xplot(1:ng), yplot(1:ng),'entropy vs T/Tc')
  CALL FigPortraitFrame(FIG)



! plot entropy vs. temperature (kelvins)
  xplot(1:ng)=c(2)*xplot(1:ng)
  WRITE(FIG,*) 'VIEWPORT 0.1 0.7  0.1 0.7'
  WRITE(FIG,*) 'WINDOW  0 1500 10 30'
  CALL FigGrid(FIG, 3,4, 0,0, 'T','S/R',  5,12, 6,8)
  WRITE(FIG,*) 'SYMBOL 1'
  CALL FigData(FIG, xplot(1:ng), yplot(1:ng),'entropy vs T/Tc')
  CALL FigPortraitFrame(FIG)


! plot internal energy vs. reduced temperature
  xplot(1:ng)=g(1:ng,1)
  yplot(1:ng)=g(1:ng,3)

  WRITE(FIG,*) 'VIEWPORT 0.1 0.7  0.1 0.7'
  WRITE(FIG,*) 'WINDOW  0 10 10 30'
  CALL FigGrid(FIG, 10,4, 0,0, 'T/Tc','U0',  5,12, 10,8)
  WRITE(FIG,*) 'SYMBOL 1'
  CALL FigData(FIG, xplot(1:ng), yplot(1:ng),'internal energy vs T/Tc')
  CALL FigPortraitFrame(FIG)


! plot iternal energy vs. temperature (kelvins)
  xplot(1:ng)=c(2)*xplot(1:ng)
  WRITE(FIG,*) 'VIEWPORT 0.1 0.7  0.1 0.7'
  WRITE(FIG,*) 'WINDOW  0 1500 10 30'
  CALL FigGrid(FIG, 3,4, 0,0, 'T','U0',  5,12, 6,8)
  WRITE(FIG,*) 'SYMBOL 1'
  CALL FigData(FIG, xplot(1:ng), yplot(1:ng),'entropy vs T/Tc')
  CALL FigPortraitFrame(FIG)

! plot specific heat vs. reduced temperature
  xplot(1:ng)=g(1:ng,1)
  yplot(1:ng)=g(1:ng,4)

  WRITE(FIG,*) 'VIEWPORT 0.1 0.7  0.1 0.7'
  WRITE(FIG,*) 'WINDOW  0 10 2 4'
  CALL FigGrid(FIG, 10,4, 0,1, 'T/Tc','Cv/R',  5,12, 10,10)
  WRITE(FIG,*) 'SYMBOL 1'
  CALL FigData(FIG, xplot(1:ng), yplot(1:ng),'specific heat vs T/Tc')
  CALL FigPortraitFrame(FIG)


! plot specific heat vs. temperature (kelvins)
  xplot(1:ng)=c(2)*xplot(1:ng)
  WRITE(FIG,*) 'VIEWPORT 0.1 0.7  0.1 0.7'
  WRITE(FIG,*) 'WINDOW  0 1500 10 30'
  CALL FigGrid(FIG, 10,4, 0,0, 'T','Cv/R',  5,12, 10,8)
  WRITE(FIG,*) 'SYMBOL 1'
  CALL FigData(FIG, xplot(1:ng), yplot(1:ng),'entropy vs T/Tc')
  CALL FigPortraitFrame(FIG)

  RETURN
END Subroutine PlotIdealGasFunctions   ! ---------------------------------------

!+
SUBROUTINE PlotCorrectionFactors(fluidName)
! ------------------------------------------------------------------------------
! PURPOSE - 
  CHARACTER(LEN=*),INTENT(IN):: fluidName

  INTEGER:: k
  CHARACTER(LEN=20):: plotLabel
  REAL,DIMENSION(1000):: xplot,yplot
!-------------------------------------------------------------------------------

! plot compressibility correction factor
  xplot(1:ng)=tt(1:n1)


  WRITE(FIG,*) 'VIEWPORT 0.1 0.7  0.1 0.7'
  WRITE(FIG,*) 'WINDOW  0 10 0.0 2.0'
  CALL FigGrid(FIG, 10,4, 0,1, 'T/Tc','Z',  5,12, 10,8)
  WRITE(FIG,*) 'SYMBOL 1'
  DO k=1,n2
    WRITE(plotLabel,'(F9.3)') dd(k)
    yplot(1:n1)=gg(1:n1,k,1)    
    CALL FigData(FIG, xplot(1:n1), yplot(1:n1),'compressibility at '//plotLabel)
  END DO
  CALL FigPortraitFrame(FIG)

  RETURN
END Subroutine PlotCorrectionFactors   ! ---------------------------------------

END Module FluidPlotsProcedures   ! ============================================


!+
PROGRAM FluidPlots
! ------------------------------------------------------------------------------
USE FluidPlotsProcedures
USE ParseMetafileToPs,ONLY:InterpretMetafileAsPs
IMPLICIT NONE


!-------------------------------------------------------------------------------
  WRITE(*,*) 'plots - plot various thermodynamic functions'
  WRITE(*,*) VERSION
  dateTimeStr=GetDateTimeStr()

  OPEN(UNIT=FIG,FILE='plots.fig',STATUS='REPLACE',ACTION='WRITE')
  OPEN(UNIT=DBG,FILE='plots.dbg',STATUS='REPLACE',ACTION='WRITE')

  CALL O2()

  CALL PlotCorrectionFactors('Oxygen')
  CALL PlotIdealGasFunctions('Oxygen')
  

  CLOSE(UNIT=FIG)
  OPEN(UNIT=FIG,FILE='plots.fig',STATUS='OLD',ACTION='READ',POSITION='REWIND')
  OPEN(UNIT=PS,FILE='plots.ps',STATUS='REPLACE',ACTION='WRITE')
  CALL InterpretMetafileAsPs(FIG,PS,DBG)
  CLOSE(UNIT=FIG)
  CLOSE(UNIT=PS)
  
  STOP
END Program FluidPlots   ! =====================================================