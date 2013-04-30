program check_imf

  use constants
  use interfaces, only : imf, imf_pop3
  implicit none 
  real(kind=prec) :: m 

  m = 8.0_prec 
  do 
     if (m > 130.0_prec) exit 
     print *, m, imf(m), imf_pop3(m)
     m = m*1.58_prec 
  end do

end program check_imf

function imf(m) 

  ! Returns the stellar mass function for mass m. Useful, for example,
  ! for the integral in outlow(z).

  use constants 
  implicit none 
  real(kind=prec), intent(in) :: m ! M_solar 
  real(kind=prec) :: imf, norm  

  ! Important: This stellar IMF is normalised such that \int_minf^msup
  ! m*imf(m)*dm = 1 M_solar.  This affects the units: the total number
  ! of stars in a mass interval is (m_star/m_solar)*imf(m)*dm; total
  ! stellar mass in the interval is (m_star/m_solar)*imf(m)*m*dm,
  ! where m_star is the total mass of stars formed.
  norm = (1.0_prec-imf_slope)/(msup**(1.0_prec-imf_slope)-minf**(1.0_prec-imf_slope))
  if ((m<minf).or.(m>msup)) then 
     imf = 0.0_prec 
  else
     imf = norm*(m**(-(1.0_prec+imf_slope))) ! M_solar^-1
  end if

end function imf

function imf_pop3(m) 

  ! Returns the stellar mass function for mass m, for Pop III
  ! stars. Useful, for example, for the integral in outlow(z).

  use constants 
  implicit none 
  real(kind=prec), intent(in) :: m ! M_solar 
  real(kind=prec) :: imf_pop3, norm, msup_pop3, minf_pop3

  minf_pop3 = 35.0_prec
  msup_pop3 = 100.0_prec 

  ! Important: This stellar IMF is normalised such that \int_minf^msup
  ! m*imf(m)*dm = 1 M_solar.  This affects the units: the total number
  ! of stars in a mass interval is (m_star/m_solar)*imf(m)*dm; total
  ! stellar mass in the interval is (m_star/m_solar)*imf(m)*m*dm,
  ! where m_star is the total mass of stars formed.
  norm = (1.0_prec-imf_slope)/(msup_pop3**(1.0_prec-imf_slope)-minf_pop3**(1.0_prec-imf_slope))
  imf_pop3 = norm*(m**(-(1.0_prec+imf_slope))) ! M_solar^-1

  if ((m<minf_pop3).or.(m>msup_pop3)) then 
     imf_pop3 = 0.0_prec 
  else
     imf_pop3 = norm*(m**(-(1.0_prec+imf_slope))) ! M_solar^-1
  end if

end function imf_pop3
