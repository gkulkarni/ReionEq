
subroutine lumfn(region, mab, z, n) 

  ! File: lumfn.f90
  ! Cre: 2010-05-17
  ! Mod: $Date: 2012/01/23 08:51:37 $; $Revision: 1.3 $ 
  !
  ! Calculates luminosity function \Phi in hot or cold region at redshift z.

  use constants 
  use storage 
  use interfaces, only : luminosity, interpolate, counter, getjmc,&
       &getjmh    
  implicit none 

  real(kind=prec), intent(in) :: mab, z ! [l]=yr^-1
  integer, intent(in) :: region 
  ! region = 1 ----> cold region 
  ! region = 2 ----> hot region 
  real(kind=prec), intent(out) :: n 

  real(kind=prec) :: sum, varz, t1, t2, rho, rgnms, l, dzl    

  if ((region /= 1) .AND. (region /= 2)) then 
     write (0, *) 'lumfn: region should be 1 or 2' 
     stop
  end if

  dzl = dz
  l = 10.0_prec**((51.60-mab)/2.5_prec) ! erg/s/Hz 

  rho = omega_nr*rho_critical ! 10^10 M_solar / Mpc^3 
  rgnms = (4.0_prec/3.0_prec)*pi*(rgnsize**3)*rho ! 10^10 m_solar
  call interpolate(sigmarr, msarr, rgnms, rgnsig)

  sum = 0.0_prec 
  varz = initial_redshift + 2.0_prec 
  t1 = nlz(varz) 
  do 
     varz = varz + dzl 
     if (varz <= z) exit 

     t2 = nlz(varz) 
     sum = sum + dzl*(t1+t2)*0.5_prec 
     t1 = t2 
  end do

  n = sum ! Mpc^-3 magnitude^-1

contains 

  function nlz(redshift) 

    implicit none 
    real(kind=prec), intent(in) :: redshift 
    real(kind=prec) :: nlz 

    real(kind=prec) :: t, tform, lum, growth_target, mass, prob_survived,&
         &mmin, sigma, dsigmadm, nu, massfn_survived, nformdot, nofm,&
         &dgrowth, dldm, mmax, dmabdl, growth_current, dnudm, age 
    integer :: grloc_target, grloc_current

    call interpolate(tarr, zarr, z, t) ! [tform] = yr
    call interpolate(tarr, zarr, redshift, tform) ! [t] = yr 
    age = t-tform ! yr 
    call gallum(1.0_prec, age, redshift, lum) ! [lum]=erg/s/Hz/(10^10M_solar)
    mass = l/lum ! 10^10 M_solar

!!$    if (mab < -19.9_prec) then 
!!$       write (0,'(f6.2,3e11.3e2)') mab, mass, age
!!$    end if

    mmax = rgnms*0.9_prec ! 10^10 m_solar

    if (region == 1) then 
       mmin = getjmc(redshift) ! 10^10 m_solar 
    else 
       mmin = getjmh(redshift) ! 10^10 m_solar 
    end if

    if (mass < mmin) then 
       nlz = 0.0_prec 
    else if (mass > mmax) then 
       nlz = 0.0_prec 
    else 
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

       dnudm = (-dsigmadm)*(nu**3)*sigma/(deltac/growth_current-rgnovd)**2
       nofm = sqrt(2.0_prec/pi)*exp(-nu*nu/2.0_prec)*(rho/mass)*dnudm
       nformdot = nu**2*deltac*dgrowth*nofm&
            &/(deltac/growth_current-rgnovd) ! Mpc^-3 (10^10 m_solar)^-1

       ! Calculate N(m, z, \tilde z) 
       massfn_survived = nformdot*prob_survived ! Mpc^-3 (10^10 m_solar)^-1 redshift^-1

       ! Calculate mass derivate of luminosity (dL/dM) 
       dldm = lum ! erg yr^-1 Hz^-1 (10^10 m_solar)^-1

       ! Calculate derivative of AB magnitude wrt luminosity
       dmabdl = 2.5_prec/l ! magnitude (erg s^-1 Hz^-1)^-1 

       nlz = massfn_survived / dldm / dmabdl ! Mpc^-3 redshift^-1 magnitude^-1
    end if

!!$    if(mass<=1.0e-3_prec .or. mass>=1.0e-2_prec) then 
!!$       nlz=0.0_prec
!!$    end if

  end function nlz

end subroutine lumfn

