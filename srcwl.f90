subroutine sourcewl(source, redshift)

  ! File: srcwl.f90 
  ! Cre: 05-06-2010
  ! Mod: $Date: 2010/06/05 05:27:19 $, $Revision: 1.1 $ 
  !
  ! Calculates source term for the reionization model of Wyithe and
  ! Loeb (2007).  This is used in reion-wl.f90.

  use constants
  use interfaces, only : counter, interpolate, nuint, luminosity, &
       &getjmc, getjmh, getff    
  use storage
  implicit none

  interface 
     function getmlow(rs, temp)
       use constants; implicit none
       real(kind = prec), intent(in) :: rs, temp 
       real(kind = prec) :: getmlow
     end function getmlow
  end interface

  real(kind = prec), intent(in) :: redshift
  real(kind = prec), intent(out) :: source 

  real(kind = prec) :: dzform, err, lum, nh, rgnms, rgnsig, rho,&
       &rho_baryon, sum, t1, t2, zform
  integer :: flag, last, neval

  !---------------------------

  rho = omega_nr*rho_critical ! 10^10 m_solar / mpc^3 
  rho_baryon = omega_b*rho_critical ! 10^10 m_solar / mpc^3 

  ! comoving hydrogen _number_ density; 1.1891e57 is the number of
  ! protons in a m_solar.  note the 10^10 due to our mass units.
  nh = rho_baryon*(1.0_prec+rgnovd)*(1.1891e57_prec*1.0e10_prec) ! mpc^-3 

  rgnms = (4.0_prec/3.0_prec)*pi*rgnsize**3*rho ! 10^10 m_solar
  call interpolate(sigmarr, msarr, rgnms, rgnsig)

  sum = 0.0_prec; dzform = dz 
  zform = initial_redshift + 2.0_prec
  t1 = time_integrand(zform)  
  do 
     zform = zform + dzform 
     if (zform <= redshift) exit 
     t2 = time_integrand(zform) 
     sum = sum + dzform * (t1 + t2) * 0.5_prec
     t1 = t2 
  end do

  source = sum ! yr^-1 mpc^-3 (10^10 msolar)

contains

  function time_integrand(z)

    implicit none
    real(kind = prec), intent(in) :: z
    real(kind = prec) :: time_integrand
    real(kind = prec) :: dmmass, mmax, numax, sum, tform, &
         &t, sigma, dsigmadm, sum2, ovdextra, growth, numinii,&
         &dgrowth, survival_probability, mmini, mminii, sum1,&
         &numini, fii
    integer :: grloc

    mmini = getmlow(z, 1.0e4_prec) ! 10^10_msolar
    mminii = getmlow(z, 1.0e5_prec) ! 10^10_msolar
    
    mmax = rgnms*0.9_prec ! 10^10 m_solar
    grloc = counter(z); growth = grfarr(grloc)
    
    dmmass = mmax ! 10^10 m_solar; total mass 
    call interpolate(sigmarr, msarr, dmmass, sigma)
    call interpolate(dsigmarr, msarr, dmmass, dsigmadm)
    numax = (deltac/growth - rgnovd)/sqrt(sigma**2-rgnsig**2)

    dmmass = mmini ! 10^10 m_solar; total mass 
    call interpolate(sigmarr, msarr, dmmass, sigma)
    call interpolate(dsigmarr, msarr, dmmass, dsigmadm)
    numini = (deltac/growth - rgnovd)/sqrt(sigma**2-rgnsig**2)

    if (mminii > mmax) then 
       numinii = numax 
    else
       dmmass = mminii ! 10^10 m_solar; total mass 
       call interpolate(sigmarr, msarr, dmmass, sigma)
       call interpolate(dsigmarr, msarr, dmmass, dsigmadm)
       numinii = (deltac/growth - rgnovd)/sqrt(sigma**2-rgnsig**2)
    end if

    !-------------------------

    call dqage(nuint,numini,numax,abserr,absrel,key,&
         &maxint,sum1,err,neval,flag,alist,blist,rlist,elist,iord,last)
    call dqage(nuint,numinii,numax,abserr,absrel,key,&
         &maxint,sum2,err,neval,flag,alist,blist,rlist,elist,iord,last)

    fii = getff(z)
    sum = (1.0_prec-fii)*sum1 + fii*sum2 ! dimensionless 

    !------------------------

    call interpolate(tarr, zarr, z, tform) ! [tform] = yr
    call interpolate(tarr, zarr, redshift, t) ! [t] = yr 
    lum = luminosity(z, tform, t) ! yr^-1

    survival_probability=(deltac/grfarr(counter(redshift))-rgnovd)/&
         &(deltac/growth-rgnovd) ! dimensionless 

    dgrowth = dgrfarr(grloc) 
    dgrowth = dgrowth/growth**2 ! dimensionless 

    ovdextra = deltac/(deltac/growth-rgnovd) ! dimensionless 

    time_integrand = sum*lum*survival_probability*dgrowth*&
         &ovdextra*rho_baryon*rho_baryon ! yr^-1 mpc^-3 (10^10 m_solar)
  end function time_integrand

end subroutine sourcewl

