program alpha_rc

  implicit none 
  integer, parameter :: prec = 8 
  real(kind=prec) :: temph, reclambda, recrate, thi, nh, &
       &rho_baryon, coolrate, erg2j, cmbympccb, yrbys, recrate_b 

  yrbys = 3.154e7_prec
  cmbympccb = 3.4036771916e-74_prec
  erg2j = 1.0e-7_PREC
  rho_baryon = 4.49e-2_prec*2.775e1*(0.71**2) ! 10^10 M_solar / Mpc^3 
  thi = 1.57807E5_prec 
  temph = 10.0_prec 
  do 
     if (temph > 1.0e8_prec) exit 
     reclambda = 2.0_prec*thi/temph 
     recrate = 1.778e-29_prec*temph*reclambda**1.965_prec/&
          &(1.0_prec+(reclambda/0.541_prec)**0.502_prec)&
          &**2.697_prec ! erg cm^3 s^-1

     recrate_b = 3.435e-30_prec * temph * reclambda**1.970d0 / &
         ( 1.0_prec + (reclambda/2.250e0_prec)**0.376_prec )**3.720_prec

     print *, temph, recrate, recrate_b 

     nh = rho_baryon*(1.1891e57_prec*1.0e10_prec) ! Mpc^-3 
     coolrate = (nh**2)*recrate*erg2j*cmbympccb*yrbys ! J Mpc^-3 yr^-1
     
     temph = temph*1.2_prec
  end do
  

end program alpha_rc
