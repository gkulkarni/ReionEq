
DOUBLE PRECISION FUNCTION fpop3(zcoll,m)

  ! file: fpop3.f90
  !  cre: -- 
  !  mod: $Date: 2012/06/24 15:24:06 $; $Revision: 1.4 $

  ! This is a fitting function to Rafaella Schneider's f_pop3 taken
  ! from Tirth's code.  It is NOT used in reion! It is simply used as
  ! a check once.

  IMPLICIT NONE

  double precision, intent(in) :: zcoll, m ![m]=M_solar (NOT 10^10 M_solar!)
  DOUBLE PRECISION :: bz,m0

  IF (zcoll > 20.d0) THEN
     fpop3=1.d0
  ELSE
     bz=0.24d0*zcoll+3.81d0
     m0=-0.05d0*zcoll+8.23d0
     fpop3=EXP(-bz*(m-m0)**2)
  END IF

  ! fpop3 = 0.0d0

END FUNCTION fpop3

