
subroutine nfchk()

  ! File: nfchk.f90 
  ! Cre: 03-06-2010
  ! Mod: $Date: 2012/03/30 12:16:47 $, $Revision: 1.3 $ 
  ! 
  ! Confirms that the time integration of halo formation rate times
  ! survival probability gives the halo mass function.  This check is
  ! important for the calculation that is implemented in lumfn.f90.

  use constants
  use storage 
  use interfaces
  implicit none 

  real(kind=prec) :: m, m_initial, m_final, sum, z, ns1, ns2, redshift,&
       &nofm, rho, rgnms, sigma, dsigmadm, nu, dnudm, &
       &dzc
  integer :: grloc 

  m_initial = 1.0e-4_prec
  m_final = 1.0e0_prec 

  rho = omega_nr*rho_critical ! 10^10 M_solar / Mpc^3 
  rgnms = (4.0_prec/3.0_prec)*pi*rgnsize**3*rho ! 10^10 M_solar
  call interpolate(sigmarr, msarr, rgnms, rgnsig)

  dzc = dz/2.0_prec
  m = m_initial
  redshift = 10.0_prec
  grloc = counter(redshift)
  growth = grfarr(grloc) 
  do 
     if (m > m_final) exit 
     
     ! Integrate number of survived haloes.
     sum = 0.0_prec 
     z = initial_redshift + 20.0_prec
     ns1 = nsurv(m, redshift, z)
     do 
        z = z + dzc 
        if (z < redshift) exit 
        
        ns2 = nsurv(m, redshift, z)
        sum = sum + 0.5_prec*dzc*(ns1+ns2) ! Mpc^-3 (10^10 m_solar)^-1
        ns1 = ns2 
     end do
     
     ! Calculate mass function.
     call interpolate(sigmarr, msarr, m, sigma)
     call interpolate(dsigmarr, msarr, m, dsigmadm)
     nu = (deltac/growth - rgnovd)/sqrt(sigma**2-rgnsig**2)
     dnudm = (-dsigmadm)*(nu**3)*sigma/((deltac/growth-rgnovd)**2)
     nofm = sqrt(2.0_prec/pi)*exp(-nu*nu/2.0_prec)*&
          &(rho/m)*dnudm ! Mpc^-3 (10^10 M_solar)^-1

     ! Print result.
     print *, m, nofm, sum 
     m = m * 1.2_prec 
  end do

contains 

  function nsurv(mass, z, redshift) 

    implicit none 
    real(kind=prec), intent(in) :: mass, z, redshift
    real(kind=prec) :: nsurv 

    real(kind=prec) :: growth_current, &
         &growth_target, dgrowth, prob_survived, nu, sigma, dsigmadm, &
         &dnudm, nofm, nformdot

    integer :: grloc_current, grloc_target 

    ! Calculate survival probability. 
    ! redshift ----------> current
    !    z     ----------> target
    grloc_current = counter(redshift)
    growth_current = grfarr(grloc_current) 

    grloc_target = counter(z)
    growth_target = grfarr(grloc_target)

    dgrowth = dgrfarr(grloc_current)/growth_current**2 ! dimensionless 
    prob_survived=(deltac/growth_target-rgnovd)/&
         &(deltac/growth_current-rgnovd) ! dimensionless 

    ! Calculate halo formation rate.
    call interpolate(sigmarr, msarr, mass, sigma)
    call interpolate(dsigmarr, msarr, mass, dsigmadm)
    nu = (deltac/growth_current - rgnovd)/sqrt(sigma**2-rgnsig**2)

    dnudm = (-dsigmadm)*(nu**3)*sigma/((deltac/growth_current-rgnovd)**2)
    nofm = sqrt(2.0_prec/pi)*exp(-nu*nu/2.0_prec)*(rho/mass)*dnudm
    nformdot = nu**2*deltac*dgrowth*nofm&
         &/(deltac/growth_current-rgnovd) ! Mpc^-3 (10^10 m_solar)^-1

    ! Calculate N(m, z, \tilde z) 
    nsurv = nformdot*prob_survived ! Mpc^-3 (10^10 m_solar)^-1 redshift^-1

  end function nsurv

end subroutine nfchk

