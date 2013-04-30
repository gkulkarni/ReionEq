
FUNCTION GETJMH(RS)

  ! File: getjmh.f90 
  ! Cre: 2009-07-02
  ! Mod: $Date: 2010/04/21 16:43:43 $ ($Revision: 1.1 $)
  ! 
  ! Gives the minimum collapsing mass 
  ! in ionized region.  This mass correspond
  ! to either a T_virial of 10^4 K or to the
  ! local Jeans mass, whichever is
  ! higher.  

  USE CONSTANTS; USE STORAGE; IMPLICIT NONE
  REAL(KIND = PREC), INTENT(IN) :: RS
  REAL(KIND = PREC) :: GETJMH

  REAL(KIND = PREC) :: A, B, RHO, GETJMH1, GETJMH2 
  INTEGER :: LOC 

  RHO = OMEGA_NR*RHO_CRITICAL ! 10^10 M_solar / Mpc^3 

  ! (MSOLKG*1.0E10) does the mass unit 
  ! conversion from kg to 10^10 M_solar. 
  A = 5.0_PREC*KBOLTZ*1.0E4_PREC / (NEWTG*MPROTON*(1.0_PREC+RS)&
       &*MSOLKG*1.0E10_PREC)

  ! (1.0E-8/CMBYMPCCB) does the length 
  ! unit conversion from Mpc to m. 
  B = 3.0_PREC*1.0E-8_PREC / (4.0_PREC*PI*DELTA_VIRIAL&
       &*RHO*CMBYMPCCB) 
  GETJMH1 = SQRT(B) * A ** (3.0_PREC / 2.0_PREC) 

  !---------------------------

  LOC = FLOOR((INITIAL_REDSHIFT-RS)/ABS(DZ)+1)
  IF (LOC < 1) THEN
     GETJMH2 = GETJMH1 
  ELSE
     GETJMH2 = JMHARR(LOC)*OMEGA_NR/OMEGA_B
  END IF

  GETJMH = MAX(GETJMH1, GETJMH2) ! 10^10 M_solar 

END FUNCTION GETJMH

