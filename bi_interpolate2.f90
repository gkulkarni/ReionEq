
SUBROUTINE BI_INTERPOLATE2(XARRAY, YARRAY, ZARRAY, X, Y, Z)

  ! File: bi_interpolate.f90
  ! Cre: 2012-02-07 
  ! Mod: $Date: 2012/02/14 08:07:08 $ ($Revision: 1.2 $) 
  !
  ! Returns value of Z = Z(X,Y) by bilinearly interpolating in the
  ! given arrays XARRAY, YARRAY, and ZARRAY.  I took the formula from
  ! Numerical Recipes.

  USE CONSTANTS; IMPLICIT NONE
  REAL(KIND=PREC), DIMENSION(:), INTENT(IN) :: XARRAY, YARRAY
  REAL(KIND=PREC), DIMENSION(:,:), INTENT(IN) :: ZARRAY
  REAL(KIND=PREC), INTENT(IN) :: X, Y
  REAL(KIND=PREC), INTENT(OUT) :: Z
  
  INTEGER :: XLOC, YLOC  
  REAL(KIND=PREC) :: M, Z1, Z2, Z3, Z4, T, U  
  
  ! Trap if array sizes don't match. 
  IF (SIZE(XARRAY)*SIZE(YARRAY) /= SIZE(ZARRAY)) THEN
     WRITE (0,*) 'BI_INTERPOLATE2: Array size mismatch error.'
  END IF

  CALL LOCATE(XARRAY, X, XLOC) 
  CALL LOCATE(YARRAY, Y, YLOC) 

  ! Following two lines are an ugly hack! Remove them ASAP!
  IF (XLOC==SIZE(XARRAY)) XLOC = XLOC-1
  IF (YLOC==SIZE(YARRAY)) YLOC = YLOC-1 

  Z1 = ZARRAY(XLOC, YLOC) 
  Z2 = ZARRAY(XLOC+1, YLOC) 
  Z3 = ZARRAY(XLOC+1, YLOC+1)
  Z4 = ZARRAY(XLOC, YLOC+1) 
  
  T = (X-XARRAY(XLOC))/(XARRAY(XLOC+1)-XARRAY(XLOC))
  U = (Y-YARRAY(YLOC))/(YARRAY(YLOC+1)-YARRAY(YLOC)) 
  
  Z = (1.0_PREC-T)*(1.0_PREC-U)*Z1 + T*(1.0_PREC-U)*Z2 + T*U*Z3 + (1.0_PREC-T)*U*Z4 
  
CONTAINS 

  SUBROUTINE LOCATE(ARRAY, QUERY, LOCATION) 

    IMPLICIT NONE 
    REAL(KIND=PREC), DIMENSION(:), INTENT(IN) :: ARRAY 
    REAL(KIND=PREC), INTENT(IN) :: QUERY 
    INTEGER, INTENT(OUT) :: LOCATION 
    
    ! Locates QUERY in ARRAY such that QUERY is between
    ! ARRAY(LOCATION) and ARRAY(LOCATION+1).  ARRAY should be ordered.
    ! However, it could be increasing or decreasing.  QUERY = 0 or N
    ! indicated an overflow.

    INTEGER :: N, JL, JU, JM 

    N = SIZE(ARRAY) 
    JL = 1
    JU = N+1 

10  IF (JU-JL > 1) THEN 
       JM = (JU+JL)/2 
       IF ((ARRAY(N)>ARRAY(1)).EQV.(QUERY>ARRAY(JM))) THEN 
          JL = JM 
       ELSE 
          JU = JM 
       END IF
       GOTO 10 
    END IF

    LOCATION = JL 
    
  END SUBROUTINE LOCATE

END SUBROUTINE BI_INTERPOLATE2



