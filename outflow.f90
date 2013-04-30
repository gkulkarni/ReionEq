function outflow(z) 

  ! File: outflow.f90 
  !  Cre: 2012-01-23
  !  Mod: $Date: 2012/10/23 08:10:11 $ ($Revision: 1.9 $) 
  !
  ! Calculates outflow rate from ISM to IGM in the model of Daigne et
  ! al. 2006.  In particular, this code implements their equation 7.
  ! See also accrate.f90.

  use constants 
  use storage 
  use interfaces, only : vesc_sq, getsfr, imf, interpolate2, getmet, &
       &bi_interpolate2, sfr_rollinde_pop3, sfr_rollinde_pop2, imf_pop3, getsfr2, getsfr3 
  implicit none 
  real(kind=prec), intent(in) :: z 
  real(kind=prec) :: outflow 

  real(kind=prec) :: t, m_d, m_up, m1, m2, int1, int2, sum, tmyr, vesq      
  integer :: i, j 

  call interpolate2(tarr, zarr, z, t) ! [t] = yr 
  tmyr = t*1.0e-6_prec ! Myr 
  call interpolate2(stellar_mass, stellar_age, tmyr, m_d) ! [md] = M_solar
  m_up = STMASS_UPLIMIT ! M_solar 

  do i = 1, DATA_COLUMNS 
     ekin(i) = yields(1,i,4) ! log(erg)
     stmass(i) = yields(1,i,2) ! Stellar mass; M_solar 
  end do 

  do i = 1, 4 
     stmet(i) = yields(i,1,1) 
     do j = 1, DATA_COLUMNS 
        e_kin(i,j) = yields(i,j,4) ! log(erg)
     end do
  end do

  m1 = max(m_d,8.0_prec)
  int1 = integrand(m1) 
  sum = 0.0_prec 
  do 
     m2 = m1*STELLAR_INTEGRAL_MULTIPLIER
     if (m2 > m_up) exit 
          
     int2 = integrand(m2) 
     sum = sum + 0.5_prec*(m2-m1)*(int1+int2) ! erg yr^-1 Mpc^-3 
     int1 = int2 
     m1 = m2 
  end do

  vesq = vesc_sq(z) ! (m/s)^2
  ! print *, z, vesq, sum

  ! The factor 1.0e-7 converts from erg to joule.
  outflow = 2.0_prec*sum*outflow_efficiency*1.0e-7_prec&
       &/(vesq*msolkg) ! M_solar yr^-1 Mpc^-3 
  outflow = outflow*1.0e-10_prec ! 10^10 M_solar yr^-1 Mpc^-3 

contains 

  function integrand(m) 

    implicit none 
    real(kind=prec), intent(in) :: m 
    real(kind=prec) :: integrand
    real(kind=prec) :: ekinetic, st_age, rs, psi2, zmet, lzmet, st_ageyr, &
         &psi3, integrand2, integrand3 
    logical :: m_overflow, m_underflow, zmet_overflow, zmet_underflow 

    call interpolate2(stellar_age, stellar_mass, m, st_age) ! [st_age] = Myr
    st_ageyr = st_age*1.0e6_prec ! yr 
    call interpolate2(zarr, tarr, t-st_ageyr, rs) 
    zmet = getmet(rs) 

    m_overflow = .false. 
    m_underflow = .false. 
    if (m > data_mmax) then 
       m_overflow = .true. 
    else if (m < data_mmin) then 
       m_underflow = .true. 
    end if

    zmet_overflow = .false. 
    zmet_underflow = .false. 
    if (zmet > data_zmetmax) then 
       zmet_overflow = .true. 
    else if (zmet < data_zmetmin) then 
       zmet_underflow = .true. 
    end if

    if (m_overflow) then 
       if (zmet_overflow) then 
          call bi_interpolate2(stmet, stmass, e_kin, data_zmetmax, data_mmax, ekinetic)
       else if (zmet_underflow) then 
          call bi_interpolate2(stmet, stmass, e_kin, data_zmetmin, data_mmax, ekinetic)
       else 
          call bi_interpolate2(stmet, stmass, e_kin, zmet, data_mmax, ekinetic)
       end if
    else if (m_underflow) then 
       if (zmet_overflow) then 
          call bi_interpolate2(stmet, stmass, e_kin, data_zmetmax, data_mmin, ekinetic)
       else if (zmet_underflow) then 
          call bi_interpolate2(stmet, stmass, e_kin, data_zmetmin, data_mmin, ekinetic)
       else 
          call bi_interpolate2(stmet, stmass, e_kin, zmet, data_mmin, ekinetic)
       end if
    else
       if (zmet_overflow) then 
          call bi_interpolate2(stmet, stmass, e_kin, data_zmetmax, m, ekinetic)
       else if (zmet_underflow) then 
          call bi_interpolate2(stmet, stmass, e_kin, data_zmetmin, m, ekinetic)
       else 
          call bi_interpolate2(stmet, stmass, e_kin, zmet, m, ekinetic)
       end if
    end if

    psi2 = getsfr2(rs) ! M_solar yr^-1 Mpc^-3
    integrand2 = imf(m)*ekinetic*psi2 ! erg yr^-1 Mpc^-3 M_solar^-1  

    psi3 = getsfr3(rs) ! M_solar yr^-1 Mpc^-3
    integrand3 = imf_pop3(m)*ekinetic*psi3 ! erg yr^-1 Mpc^-3 M_solar^-1  

    if (m > POP2_UPLIMIT) integrand2 = 0.0_prec 

    integrand = integrand2 + integrand3 ! erg yr^-1 Mpc^-3 M_solar^-1  

  end function integrand
  
end function outflow

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
  real(kind=prec) :: imf_pop3, norm

  ! Important: This stellar IMF is normalised such that \int_minf^msup
  ! m*imf(m)*dm = 1 M_solar.  This affects the units: the total number
  ! of stars in a mass interval is (m_star/m_solar)*imf(m)*dm; total
  ! stellar mass in the interval is (m_star/m_solar)*imf(m)*m*dm,
  ! where m_star is the total mass of stars formed.
  norm = (1.0_prec-imf_slope)/(msup_pop3**(1.0_prec-imf_slope)-minf_pop3**(1.0_prec-imf_slope))
  if ((m<minf_pop3).or.(m>msup_pop3)) then 
     imf_pop3 = 0.0_prec 
  else
     imf_pop3 = norm*(m**(-(1.0_prec+imf_slope))) ! M_solar^-1
  end if

end function imf_pop3

function vesc_sq(z) 

  use constants 
  use interfaces, only : getjmh, getjmc, interpolate2, counter, getff, fbint, fescint
  use storage 
  implicit none 
  real(kind=prec), intent(in) :: z 
  real(kind=prec) :: vesc_sq 

  real(kind=prec) :: mmin_i, mmin_ii, sigma_i, sigma_ii, grwth, numin_i, &
       &numin_ii, q, sum_i, sum_ii, dr, nr, error 
  integer :: grloc, inf, last, neval, ier  

  mmin_i = getjmc(z) 
  ! mmin_i = 1.0e-2_prec ! CHECK*****
  call interpolate2(sigmarr, msarr, mmin_i, sigma_i)
  mmin_ii = getjmh(z) 
  ! mmin_ii = 1.0e-2_prec ! CHECK*****
  call interpolate2(sigmarr, msarr, mmin_ii, sigma_ii)
  grloc = counter(z)
  grwth = grfarr(grloc)

  numin_i = deltac/(grwth*sigma_i)
  numin_ii = deltac/(grwth*sigma_ii) 

  q = getff(z) 

  inf = 1; ier = -1  
  call dqagie(fbint,numin_i,inf,abserr,absrel,maxint,sum_i,&
       &error,neval,ier,alist,blist,rlist,elist,iord,last)
  if (ier > 0) write (0,*) 'accrate: Error', ier
  call dqagie(fbint,numin_ii,inf,abserr,absrel,maxint,sum_ii,&
       &error,neval,ier,alist,blist,rlist,elist,iord,last)
  if (ier > 0) write (0,*) 'accrate: Error', ier

  dr = (1.0_prec-q)*sum_i + q*sum_ii
  zvesc = z 

  inf = 1; ier = -1 
  call dqagie(fescint,numin_i,inf,abserr,absrel,maxint,sum_i,&
       &error,neval,ier,alist,blist,rlist,elist,iord,last)
  if (ier > 0) write (0,*) 'accrate: Error', ier
  call dqagie(fescint,numin_ii,inf,abserr,absrel,maxint,sum_ii,&
       &error,neval,ier,alist,blist,rlist,elist,iord,last)
  if (ier > 0) write (0,*) 'accrate: Error', ier

  nr = (1.0_prec-q)*sum_i + q*sum_ii

  vesc_sq = nr/dr ! (m/s)^2 

end function vesc_sq

function fescint(nu) 

  ! Returns integrand for the integral in fescint.f90.
  
  use constants 
  use storage 
  use interfaces, only : counter, interpolate2 
  implicit none 
  real(kind=prec), intent(in) :: nu 
  real(kind=prec) :: fescint 

  real(kind=prec) :: grwth, m, r, sgm, z     
  integer :: grloc 

  z = zvesc; grloc = counter(z); grwth = grfarr(grloc)
  sgm = sqrt(rgnsig**2 + ((deltac/grwth-rgnovd)/nu)**2)
  call interpolate2(msarr, sigmarr, sgm, m) ! [m] = 10^10 M_solar 
  r = ((3.0_prec*m)/(4.0_prec*pi*rho_critical*omega_nr*delta_virial))**(1.0_prec/3.0_prec) ! Mpc
  fescint = sqrt(2.0_prec/pi)*exp(-nu*nu/2.0_prec)*(2.0_prec*newtg*&
       &1.0e12_prec*m*msolkg*cmbympc/r) ! (m/s)^2 [1.0e12 combines
                                        ! msun-10^10msun and cm-m
                                        ! conversion factors.]
end function fescint
