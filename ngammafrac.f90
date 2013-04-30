
function ngammafrac(z) 

  ! File: ngammafrac.f90 
  !  Cre: 2012-07-06
  !  Mod: $Date: 2012/10/23 08:10:11 $ ($Revision: 1.3 $) 

  use constants 
  use storage 
  use interfaces, only : imf, imf_pop3, interpolate2, bi_interpolate2, getmet
  implicit none 
  real(kind=prec), intent(in) :: z 
  real(kind=prec) :: ngammafrac

  real(kind=prec) :: t, m_d, m_up, m1, m2, int1, int2, sum, tmyr      
  integer :: i, j 
  logical :: all_species 

  call interpolate2(tarr, zarr, z, t) ! [t] = yr 
  tmyr = t*1.0e-6_prec ! Myr 
  call interpolate2(pop3_stellar_mass, pop3_stellar_age, tmyr, m_d) ! [md] = M_solar
  m_up = STMASS_UPLIMIT ! M_solar 

  m1 = max(m_d,8.0_prec)
  int1 = integrand(m1) 
  sum = 0.0_prec 
  do 
     m2 = m1*STELLAR_INTEGRAL_MULTIPLIER 
     if (m2 > m_up) exit 
          
     int2 = integrand(m2) 
     sum = sum + 0.5_prec*(m2-m1)*(int1+int2) ! dimensionless
     int1 = int2 
     m1 = m2
  end do

  ngammafrac = sum ! M_solar^-1 (because we divide sum by 1 M_solar.)

contains 

  function integrand(m) 

    implicit none 
    real(kind=prec), intent(in) :: m 
    real(kind=prec) :: integrand
    real(kind=prec) :: mr, st_age, rs, psi2, psi3, mej, nphotons, &
         &lzmet, st_ageyr, integrand2, integrand3, lifetime 
    logical :: m_overflow, m_underflow, zmet_overflow, zmet_underflow 

    call interpolate2(pop3_stellar_age, pop3_stellar_mass, m, st_age) ! [st_age] = Myr
    st_ageyr = st_age*1.0e6_prec ! yr 

    m_overflow = .false. 
    m_underflow = .false. 
    if (m > 500.0_prec) then 
       m_overflow = .true. 
    else if (m < 5.0_prec) then 
       m_underflow = .true. 
    end if

    if (m_overflow) then 
       integrand = pop3_stellar_ngamma(12)*pop3_stellar_age(12)*yrbys*1.0e6_prec
    else if (m_underflow) then 
       integrand = pop3_stellar_ngamma(1)*pop3_stellar_age(1)*yrbys*1.0e6_prec
    else
       call interpolate2(pop3_stellar_age, pop3_stellar_mass, m, lifetime)
       call interpolate2(pop3_stellar_ngamma, pop3_stellar_mass, m, nphotons)
       integrand = nphotons * lifetime * yrbys * 1.0e6_prec ! dimensionless (no. of photons)
    end if

    integrand = imf_pop3(m) * integrand ! M_solar^-1

  end function integrand

end function ngammafrac


