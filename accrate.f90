function accrate(z) 

  ! File: accrate.f90 
  ! Cre: 2012-01-20 
  ! Mod: $Date: 2012/05/09 07:07:12 $ ($Revision: 1.5 $)
  !
  ! Returns average baryonic accretion rate on haloes at redshift z.
  ! Implements Eq. 6 of Daigne et al. 2006 (Ap. J. 647, 773).

  use constants 
  use interfaces, only : getjmh, getjmc, interpolate2, counter, getff, fbint
  use storage 
  implicit none 
  real(kind=prec), intent(in) :: z 
  real(kind=prec) :: accrate 
  real(kind=prec) :: mmin_i, mmin_ii, sigma_i, sigma_ii, grwth, numin_i, &
       &numin_ii, q, sum_i, sum_ii, sum_z1, sum_z2, dfbdz, error 
  integer :: grloc, inf, last, neval, ier  

  ! Calculate f_b,struct at z 
  mmin_i = getjmc(z) 
  call interpolate2(sigmarr, msarr, mmin_i, sigma_i)
  mmin_ii = getjmh(z) 
  call interpolate2(sigmarr, msarr, mmin_ii, sigma_ii)
  grloc = counter(z)
  grwth = grfarr(grloc)

  numin_i = deltac/(grwth*sigma_i)
  numin_ii = deltac/(grwth*sigma_ii) 

  q = getff(z) 

  inf = 1
  call dqagie(fbint,numin_i,inf,abserr,absrel,maxint,sum_i,&
       &error,neval,ier,alist,blist,rlist,elist,iord,last)
  if (ier > 0) write (0,*) 'accrate: Error', ier
  call dqagie(fbint,numin_ii,inf,abserr,absrel,maxint,sum_ii,&
       &error,neval,ier,alist,blist,rlist,elist,iord,last)
  if (ier > 0) write (0,*) 'accrate: Error', ier

  sum_z1 = (1.0_prec-q)*sum_i + q*sum_ii

  ! Calculate f_b,struct at z-dz 
  mmin_i = getjmc(z-dz) 
  call interpolate2(sigmarr, msarr, mmin_i, sigma_i)
  mmin_ii = getjmh(z-dz) 
  call interpolate2(sigmarr, msarr, mmin_ii, sigma_ii)
  grloc = counter(z-dz)
  grwth = grfarr(grloc)

  numin_i = deltac/(grwth*sigma_i)
  numin_ii = deltac/(grwth*sigma_ii) 

  q = getff(z-dz) 

  inf = 1
  call dqagie(fbint,numin_i,inf,abserr,absrel,maxint,sum_i,&
       &error,neval,ier,alist,blist,rlist,elist,iord,last)
  if (ier > 0) write (0,*) 'accrate: Error', ier
  call dqagie(fbint,numin_ii,inf,abserr,absrel,maxint,sum_ii,&
       &error,neval,ier,alist,blist,rlist,elist,iord,last)
  if (ier > 0) write (0,*) 'accrate: Error', ier

  sum_z2 = (1.0_prec-q)*sum_i + q*sum_ii

  ! Calculate the derivative 
  dfbdz = abs((sum_z1-sum_z2)/dz) ! dimensionless 
  accrate = 1.0e-10_prec*dfbdz*sqrt(omega_lambda+omega_nr*(1.0_prec+z)**3)*&
       &(1.0_prec+z)*(omega_b/0.044_prec)*1.2_prec*smallh**3 ! (10^10*M_solar)/(yr*Mpc^3)

end function accrate

function fbint(nu) 

  use constants 
  implicit none 
  real(kind=prec), intent(in) :: nu 
  real(kind=prec) :: fbint 

  ! Returns integrand for the integral in accrate.f90 and outflow.f90.

  fbint = sqrt(2.0_prec/pi)*exp(-nu*nu/2.0_prec)

end function fbint

