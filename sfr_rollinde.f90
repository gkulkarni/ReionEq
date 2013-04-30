
! File: sfr_rollinde.f90 
!  Cre: 2012
!  Mod: $Date$ ($Revision$) 

! Implements the SFR using in Rollinde et al. (2009).  This was used
! to check my code by confirming that it reproduces the results.


function sfr_rollinde_pop2(z)

  use constants
  implicit none 
  real(kind = prec), intent(in) :: z 
  real(kind = prec) :: sfr_rollinde_pop2 

  real(kind = prec) :: nu, zm, a, b 

  nu = 0.3_prec
  zm = 2.6_prec
  a = 1.9_prec
  b = 1.2_prec

  sfr_rollinde_pop2 = nu * ((a*exp(b*(z-zm)))/(a-b+b*exp(a*(z-zm)))) ! M_solar yr^-1 Mpc^-3

end function sfr_rollinde_pop2 

function sfr_rollinde_pop3(z)

  use constants
  implicit none 
  real(kind = prec), intent(in) :: z 
  real(kind = prec) :: sfr_rollinde_pop3

  real(kind = prec) :: nu, zm, a, b 

  nu = 0.0016_prec
  zm = 22.8_prec
  a = 4.0_prec
  b = 3.3_prec

  sfr_rollinde_pop3 = nu * ((a*exp(b*(z-zm)))/(a-b+b*exp(a*(z-zm)))) ! M_solar yr^-1 Mpc^-3
  ! sfr_rollinde_pop3 = 0.0_prec ! M_solar yr^-1 Mpc^-3

end function sfr_rollinde_pop3

