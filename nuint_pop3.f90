FUNCTION NUINT_POP3(NU)

  ! File: nuint_pop3.f90
  ! Cre: 2012-01-05
  ! Mod: $Date: 2012/06/24 15:24:06 $ ($Revision: 1.4 $) 

  USE CONSTANTS; use storage 
  use interfaces, only : fpop3, interpolate2
  IMPLICIT NONE 
  REAL(KIND = PREC), INTENT(IN) :: NU 
  REAL(KIND = PREC) :: NUINT_POP3

  real(kind = prec) :: sgm, m, f3   

  sgm = sqrt(rgnsig**2 + ((deltac/growth-rgnovd)/nu)**2)
  CALL INTERPOLATE2(MSARR, SIGMARR, sgm, m)
  f3 = fpop3(zform, log10(m)+10.0_prec) ! Added 10.0 to convert mass unit to M_solar
  ! WRITE (0,*) '>>>>>', f3, zform, log10(m)   

  NUINT_POP3 = f3*SQRT(2.0_PREC/PI)*EXP(-NU*NU/2.0_PREC)*NU*NU 

END FUNCTION NUINT_POP3

