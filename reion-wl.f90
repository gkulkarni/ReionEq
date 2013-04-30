
program reion_wl 

  ! File : reion-wl.f90 
  ! Cre : 03-05-2010
  ! Mod : $Date: 2010/06/05 05:27:19 $, $Revision: 1.1 $ 
  !
  ! This implements the simple reionization model of Wyithe and Loeb
  ! (2007).  

  use constants
  use storage
  use interfaces, only: sourceov, dtdz, counter 
  implicit none 

  interface 
     subroutine sourcewl(source, redshift)
       use constants; use storage; implicit none
       real(kind = prec), intent(in) :: redshift
       real(kind = prec), intent(out) :: source 
     end subroutine sourcewl
     function getmlow(rs, temp)
       use constants; implicit none
       real(kind = prec), intent(in) :: rs, temp 
       real(kind = prec) :: getmlow
     end function getmlow
  end interface

  real(kind=prec) :: q, z, rho, rho_baryon, nh, clump, ngamma, grf,&
       &dgrf, freq, ndotm, h0t0, t0sec, t0yr, m, sig, dsig, htt, &
       &zoft, lb_age, lb_lambda, lb_lum, reclambda, t, recbrate

  integer :: number_of_lines, i, ncalc

  include 'readin.inc' 

  rgnovd = 4.0_prec 
  rgnsize = 7.746_prec 
  fesc = 1.0_prec 

  rho = omega_nr*rho_critical ! 10^10 M_solar / Mpc^3 
  rho_baryon = omega_b*rho_critical ! 10^10 M_solar / Mpc^3 
  nh = rho_baryon*(1.1891e57_prec*1.0e10_prec) ! Mpc^-3 
  ngamma = 2.82e70_prec ! (10^10 M_solar)^-1 

  ncalc = (final_redshift-initial_redshift)/dz+1 
  allocate(jmharr(ncalc)); allocate(ffarr(ncalc))
  jmharr = 0.0_prec; ffarr = 0.0_prec 

!!$  t = 3.0e4_prec
!!$  reclambda = 2.d0 * thi / t
!!$  recbrate = 2.753d-14 * reclambda**1.500d0 / &
!!$       ( 1.d0 + (reclambda/2.740d0)**0.407d0 )**2.242d0
!!$  print *, recbrate
!!$  stop

  clump = 10.0_prec 

  countr = 1 
  q = 0.0e-20_prec 
  z = initial_redshift
  ffarr(countr) = q 
  do 
     countr = countr + 1 
     if (z < final_redshift) exit 
     q = q + dqdz(z)*dz
     z = z + dz 
     ffarr(countr) = q 
     print *, z, q 
  end do

contains

  function dqdz(z) 

    implicit none 
    real(kind=prec), intent(in) :: z 
    real(kind=prec) :: dqdz 

    real(kind=prec) :: nh_proper, src, recb, dofz 

    call sourcewl(src, z) ! [src] = yr^-1 Mpc^-3 (10^10 M_solar)
    nphdot = fesc*src*ngamma ! yr^-1 Mpc^-3 
    recb = 9.89e-14_prec ! cm^3/s
    dofz = grfarr(counter(z))
    nh_proper = nh*(1.0_prec+rgnovd*dofz)*(1.0_prec+z)**3 ! Mpc^-3
    dqdz =  nphdot*dtdz(z)/nh - &
         &recb*cmbympccb*yrbys*nh_proper*q*clump*dtdz(z) 

  end function dqdz

end program reion_wl

