  
  ! File: tabz.f90
  ! CRE: 2009-10-24
  ! MOD: $Date: 2010/06/05 05:27:19 $; $Revision: 1.1 $
  ! 
  ! Finds parameter theta of the solution to spherical 
  ! collapse problem by solving Eqn (7) of Munoz and Loeb (2008)
  ! using NETLIB routine ZEROIN.  Needs zeroin.f and d1mach.f.
  ! Currently called in pdelta.f90.


MODULE STOTABZ

  USE CONSTANTS
  IMPLICIT NONE
  REAL(KIND = PREC) :: DL, Z
  
END MODULE STOTABZ

FUNCTION THETA(DELTA_L, REDSHIFT)

  USE CONSTANTS; USE STOTABZ; IMPLICIT NONE
  INTERFACE

     FUNCTION FTH(X)

       USE CONSTANTS
       IMPLICIT NONE
       REAL(KIND = PREC), INTENT(IN) :: X
       REAL(KIND = PREC) :: FTH

     END FUNCTION FTH

  END INTERFACE

  REAL(KIND = PREC), INTENT(IN) :: DELTA_L, REDSHIFT
  REAL(KIND = PREC) :: THETA

  REAL(KIND = PREC) :: SOL, AX, BX, TOL, ZEROIN

  DL = DELTA_L
  Z = REDSHIFT

  AX = 1.0_PREC
  BX = 40.0_PREC
  TOL = 0.01_PREC
  THETA = ZEROIN(AX,BX,FTH,TOL)

END FUNCTION THETA

FUNCTION FTH(X)
  
  USE CONSTANTS; USE STOTABZ; IMPLICIT NONE
  REAL(KIND = PREC), INTENT(IN) :: X
  REAL(KIND = PREC) :: FTH

  FTH = X - SIN(X) - (DL*20.0_PREC/((1.0_PREC+Z)&
       &*3.0_PREC*6.0_PREC**(2.0_PREC/3.0_PREC)))&
       &**(3.0_PREC/2.0_PREC)

END FUNCTION FTH
