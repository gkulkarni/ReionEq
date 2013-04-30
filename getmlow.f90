
function getmlow(rs, temp)

  ! File: getmlow.f90 
  ! Cre: 05-2010
  ! Mod: $Date: 2010/06/05 05:27:19 $, $Revision: 1.1 $ 
  ! 
  ! Calculates mass corresponding to T_vir = 10^4 K and 10^5 K.  This
  ! is used in the simple reionization model of Wyithe and Loeb
  ! (2007), which is implemented in reion-wl.f90.

  use constants; implicit none
  real(kind = prec), intent(in) :: rs, temp 
  real(kind = prec) :: getmlow
  real(kind = prec) :: a, b, rho, t 

  !---------------------------

  rho = omega_nr*rho_critical ! 10^10 m_solar / mpc^3 
  t = temp ! k 

  ! (msolkg*1.0e10) does the mass unit 
  ! conversion from kg to 10^10 m_solar. 
  a = 5.0_prec*kboltz*t / (newtg*mproton*(1.0_prec+rs)&
       &*msolkg*1.0e10_prec)

  ! (1.0e-8/cmbympccb) does the length 
  ! unit conversion from mpc to m. 
  b = 3.0_prec*1.0e-8_prec / (4.0_prec*pi*delta_virial&
       &*rho*cmbympccb) 

  !---------------------------

  getmlow = sqrt(b)*a**(3.0_prec/2.0_prec) ! 10^10 m_solar 

end function getmlow

