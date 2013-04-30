subroutine calculate_nuintegral(z)

  ! Cre: 2012-05-27 
  ! Mod: $Date: 2012/05/27 16:20:25 $
  
  ! Calculate the nu-integral in the SFR expression.

  use constants
  use storage 
  use interfaces, only : interpolate2, getjmc, getjmh, counter, nuint_pop2, nuint_pop3, getff 
  implicit none
  real(kind=prec), intent(in) :: z 
  real(kind=prec) :: rho, rho_baryon, mmax, dmmass, sigma, dsigmadm, numax, &
       &mmini, mminii, numini, numinii, sum, sum1, sum2, rgnms, fii, err
  integer :: grloc, flag, last, neval

  rho = omega_nr*rho_critical ! 10^10 M_solar / Mpc^3 
  rho_baryon = omega_b*rho_critical ! 10^10 M_solar / Mpc^3 

  rgnms = (4.0_prec/3.0_prec)*pi*rgnsize**3*rho ! 10^10 M_solar
  call interpolate2(sigmarr, msarr, rgnms, rgnsig)

  mmini = getjmc(z)
  mminii = getjmh(z)

  mmax = rgnms*0.9_PREC ! 10^10 M_solar
  grloc = counter(z); growth = grfarr(grloc)

  dmmass = mmax ! 10^10 M_solar; total mass 
  call interpolate2(sigmarr, msarr, dmmass, sigma)
  call interpolate2(dsigmarr, msarr, dmmass, dsigmadm)
  numax = (deltac/growth - rgnovd)/sqrt(sigma**2-rgnsig**2)

  dmmass = mmini ! 10^10 M_solar; total mass 
  call interpolate2(sigmarr, msarr, dmmass, sigma)
  call interpolate2(dsigmarr, msarr, dmmass, dsigmadm)
  numini = (deltac/growth - rgnovd)/sqrt(sigma**2-rgnsig**2)

  if (mminii > mmax) then 
     numinii = numax 
  else
     dmmass = mminii ! 10^10 M_solar; total mass 
     call interpolate2(sigmarr, msarr, dmmass, sigma)
     call interpolate2(dsigmarr, msarr, dmmass, dsigmadm)
     numinii = (deltac/growth - rgnovd)/sqrt(sigma**2-rgnsig**2)
  end if

  call dqage(nuint_pop2,numini,numax,abserr,absrel,key,&
       &maxint,sum1,err,neval,flag,alist,blist,rlist,elist,iord,last)
  call dqage(nuint_pop2,numinii,numax,abserr,absrel,key,&
       &maxint,sum2,err,neval,flag,alist,blist,rlist,elist,iord,last)

  fii = getff(z)
  sum = (1.0_prec-fii)*sum1 + fii*sum2 ! dimensionless 

  nuintegralarr_pop2(countr) = sum

  call dqage(nuint_pop3,numini,numax,abserr,absrel,key,&
       &maxint,sum1,err,neval,flag,alist,blist,rlist,elist,iord,last)
  call dqage(nuint_pop3,numinii,numax,abserr,absrel,key,&
       &maxint,sum2,err,neval,flag,alist,blist,rlist,elist,iord,last)

  sum = (1.0_prec-fii)*sum1 + fii*sum2 ! dimensionless 
  ! write (0,*) z, sum, sum1, sum2 
  ! write (0,*) z, countr
  nuintegralarr_pop3(countr) = sum

end subroutine calculate_nuintegral
