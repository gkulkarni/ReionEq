
SUBROUTINE INTERPOLATE2(YARRAY, XARRAY, X, Y)

  ! File: interpolate.f90
  ! Cre: 2010-04-07 
  ! Mod: $Date: 2012/02/13 08:23:11 $ ($Revision: 1.3 $) 
  !
  ! Returns value of Y = Y(X) by linearly interpolating
  ! in the given arrays XARRAY and YARRAY.

  USE CONSTANTS; IMPLICIT NONE
  REAL(KIND=PREC), DIMENSION(:), INTENT(IN) :: YARRAY, XARRAY
  REAL(KIND=PREC), INTENT(IN) :: X
  REAL(KIND=PREC), INTENT(OUT) :: Y 
  LOGICAL :: NEXT 

  REAL(KIND=PREC), DIMENSION(:), ALLOCATABLE :: DIFF 
  INTEGER, DIMENSION(1) :: LOC 
  REAL(KIND=PREC) :: M 

  ! Trap if YARRAY and XARRAY don't match.
  IF (SIZE(XARRAY) /= SIZE(YARRAY)) THEN
     WRITE (0,*) 'INTERPOLATE2: Array size mismatch error.'
     print*, size(xarray), size(yarray) 
  END IF

  ! Prepare the locator.
  ALLOCATE(DIFF(SIZE(XARRAY)))
  DIFF = ABS(X-XARRAY)
  LOC = MINLOC(DIFF)

  ! Trap if location the last array cell.
!!$  IF (LOC(1) == SIZE(XARRAY)) THEN 
!!$     WRITE (0,*) 'INTERPOLATE2: Array size exceeded.'
!!$  END IF
  
  IF(XARRAY(LOC(1))<X .AND. X<XARRAY(LOC(1)+1)) THEN 
     NEXT = .TRUE. 
  ELSE 
     NEXT = .FALSE. 
  END IF
  
  ! Linearly interpolate.
  IF (NEXT) THEN 
     M = (YARRAY(LOC(1)+1)-YARRAY(LOC(1)))/(XARRAY(LOC(1)+1)-XARRAY(LOC(1)))
     Y = YARRAY(LOC(1)) + M*(X-XARRAY(LOC(1)))
  ELSE 
     M = (YARRAY(LOC(1))-YARRAY(LOC(1)-1))/(XARRAY(LOC(1))-XARRAY(LOC(1)-1))
     Y = YARRAY(LOC(1)-1) + M*(X-XARRAY(LOC(1)-1))
  END IF
     
  DEALLOCATE(DIFF) 

END SUBROUTINE INTERPOLATE2



