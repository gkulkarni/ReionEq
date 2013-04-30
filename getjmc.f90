
FUNCTION GETJMC(RS)

  ! File: getjmc.f90 
  ! Cre: 2009-07-02
  ! Mod: $Date: 2010/04/21 16:43:43 $ ($Revision: 1.1 $) 
  ! 
  ! Gives the mass corresponding to a 
  ! T_virial of 10^4 K.  This is the 
  ! atomic cooling bound on the lower 
  ! mass of sources in neutral regions.  

  USE CONSTANTS; IMPLICIT NONE
  REAL(KIND = PREC), INTENT(IN) :: RS
  REAL(KIND = PREC) :: GETJMC
  REAL(KIND = PREC) :: A, B, RHO, T 

  !---------------------------

  RHO = OMEGA_NR*RHO_CRITICAL ! 10^10 M_solar / Mpc^3 
  T = 1.0E4_PREC ! K 

  ! (MSOLKG*1.0E10) does the mass unit 
  ! conversion from kg to 10^10 M_solar. 
  A = 5.0_PREC*KBOLTZ*T / (NEWTG*MPROTON*(1.0_PREC+RS)&
       &*MSOLKG*1.0E10_PREC)

  ! (1.0E-8/CMBYMPCCB) does the length 
  ! unit conversion from Mpc to m. 
  B = 3.0_PREC*1.0E-8_PREC / (4.0_PREC*PI*DELTA_VIRIAL&
       &*RHO*CMBYMPCCB) 

  !---------------------------

  GETJMC = SQRT(B)*A**(3.0_PREC/2.0_PREC) ! 10^10 M_solar 

END FUNCTION GETJMC

