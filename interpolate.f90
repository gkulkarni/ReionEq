
SUBROUTINE INTERPOLATE(YARRAY, XARRAY, X, Y)

  ! File: interpolate.f90
  ! Cre: 2010-04-07 
  ! Mod: $Date: 2010/05/03 14:14:38 $ ($Revision: 1.3 $) 
  !
  ! Returns value of Y = Y(X) by linearly interpolating
  ! in the given arrays XARRAY and YARRAY.

  USE CONSTANTS; IMPLICIT NONE
  REAL(KIND=PREC), DIMENSION(:), INTENT(IN) :: YARRAY, XARRAY
  REAL(KIND=PREC), INTENT(IN) :: X
  REAL(KIND=PREC), INTENT(OUT) :: Y 

  REAL(KIND=PREC), DIMENSION(:), ALLOCATABLE :: DIFF 
  INTEGER, DIMENSION(1) :: LOC 
  REAL(KIND=PREC) :: M 

  ! Trap if YARRAY and XARRAY don't match.
  IF (SIZE(XARRAY) /= SIZE(YARRAY)) THEN
     WRITE (0,*) 'INTERPOLATE: Array size mismatch error.'
  END IF

  ! Prepare the locator.
  ALLOCATE(DIFF(SIZE(XARRAY)))
  DIFF = ABS(X-XARRAY)
  LOC = MINLOC(DIFF)

  ! Trap if location the last array cell.
  IF (LOC(1) == SIZE(XARRAY)) THEN 
     WRITE (0,*) 'INTERPOLATE: Array size exceeded.'
  END IF

  ! Linearly interpolate.
  M = (YARRAY(LOC(1)+1)-YARRAY(LOC(1)))/(XARRAY(LOC(1)+1)-XARRAY(LOC(1)))
  Y = YARRAY(LOC(1)) + M*(X-XARRAY(LOC(1)))

  DEALLOCATE(DIFF) 

END SUBROUTINE INTERPOLATE

