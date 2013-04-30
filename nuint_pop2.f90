FUNCTION NUINT_POP2(NU)

  ! File: nuint_pop2.f90
  ! Cre: 2012-01-05
  ! Mod: $Date: 2012/02/13 08:23:11 $ ($Revision: 1.3 $) 

  USE CONSTANTS; USE STORAGE
  use interfaces, only : fpop3, interpolate2 
  implicit none 
  REAL(KIND = PREC), INTENT(IN) :: NU 
  REAL(KIND = PREC) :: NUINT_POP2
  real(kind = prec) :: sgm, m, f3   

  sgm = sqrt(rgnsig**2 + ((deltac/growth-rgnovd)/nu)**2)
  CALL INTERPOLATE2(MSARR, SIGMARR, sgm, m)
  f3 = fpop3(zform, log10(m)+10.0_prec) ! Added 10.0 to convert mass unit to M_solar

  NUINT_POP2 = (1.0_prec-f3)*SQRT(2.0_PREC/PI)*EXP(-NU*NU/2.0_PREC)*NU*NU 

END FUNCTION NUINT_POP2

