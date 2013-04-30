SUBROUTINE ZBRAC(FUNC, X1, X2) 

  ! File: zbrac.f90
  ! Cre: 2009-12-28
  ! Mod: $Date: 2010/04/19 05:20:46 $ ($Revision: 1.1 $)
  !
  ! Used to bracket a root of function 
  ! FUNC(X).  Takes in an initial bracket
  ! guess and expands it geometrically
  ! until a bracket is obtained.  See
  ! Numerical Recipes 9.1.

  USE CONSTANTS; IMPLICIT NONE 
  INTEGER, PARAMETER :: NTRY = 1000
  REAL(KIND=PREC), INTENT(INOUT) :: X1, X2
  INTERFACE
     FUNCTION FUNC(X)
       USE CONSTANTS; IMPLICIT NONE
       REAL(KIND=PREC), INTENT(IN) :: X
       REAL(KIND=PREC) :: FUNC
     END FUNCTION FUNC
  END INTERFACE
  REAL(KIND=PREC), PARAMETER :: FACTOR = 1.6
  REAL(KIND=PREC) :: F1, F2 
  INTEGER :: I 

  IF(X1==X2)WRITE(0,'(A34)')'Error: Bad initial range in ZBRAC.'
  F1=FUNC(X1)
  F2=FUNC(X2)
  DO I=1,NTRY
     IF ((F1*F2)<0.0_PREC) RETURN
     IF (ABS(F1)<ABS(F2)) THEN
        X1=X1+FACTOR*(X1-X2)
        F1=FUNC(X1)
     ELSE
        X2=X2+FACTOR*(X2-X1)
        F2=FUNC(X2)
     END IF
  END DO
  WRITE(0,'(A20)')'Error: ZBRAC failed.'

END SUBROUTINE ZBRAC

