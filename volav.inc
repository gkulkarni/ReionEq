  
  ! Calculate *mean* case-A HII recombination coefficient (alpha_r).
  ! This expression is from Hui and Gnedin (1997).
  reclambda = 2.0_prec*thi/temphva
  alpha_r = 1.269e-13_prec * reclambda**1.503_prec / &
       (1.0_prec+(reclambda/0.522_prec)**0.470_prec)**1.923_prec ! cm^3 s^-1 
  ! Add case-B recombination coeff. to that. 
  alpha_r = alpha_r + 2.753e-14_prec * (reclambda**1.5_prec) / &
       (1.0_prec + (reclambda/2.74_prec)**0.407_prec)**2.242_prec ! cm^3 s^-1 


  !-------------------------

  ! Update mean ionized hydrogen fraction (x_iiva) in ionized region. 
  gamma_rec = nh_proper*alpha_r*cmbympccb*yrbys ! yr^-1 
  gamma_ion = (gpi_pop2*source_pop2 + gpi_pop3*source_pop3)*lmfp*&
       &fesc*(1.0_prec+z)**3*(cmbympc**2)/q ! yr^-1 

  a = gamma_rec*dz*dtdz(z)
  b = 1.0_prec + dz*dtdz(z)*gamma_ion
  c = -(oldxiiva + dz*dtdz(z)*gamma_ion)

  qe = -0.5_prec*(b+sign(dsqrt(b**2-4.0_prec*a*c),b))
  r1 = qe/a; r2 = c/qe 

  oldxiiva = x_iiva 
  if (r1<=1.0_prec .and. r1>=0.0_prec) then
     x_iiva = r1
  else 
     x_iiva = r2 
  end if

  !-------------------------

  ! Calculate *mean* case-A HII recombination cooling rate (recrate).
  ! This expression is from Hui and Gnedin (1997).
  reclambda = 2.0_prec*thi/temphva 
  recrate = 1.778e-29_prec*temphva*reclambda**1.965_prec/&
       &(1.0_prec+(reclambda/0.541_prec)**0.502_prec)&
       &**2.697_prec ! erg cm^3 s^-1
  ! Add case-B HII recombination cooling rate. 
  recrate = recrate + 3.435e-30_prec*temphva*reclambda**1.970d0 / &
       (1.0_prec+(reclambda/2.250e0_prec)**0.376_prec)**3.720_prec ! erg cm^3 s^-1

  ! Update *mean* temperature of ionized regions (temphva).
  nhii = oldxiiva 
  gamma_rc = (nhii**2)*recrate*erg2j*cmbympccb*yrbys ! j mpc^-3 yr^-1
  gamma_totc = gamma_rc*(1.0_prec+z)**6 ! j mpc^-3 yr^-1

  ! Add Compton cooling 
  gamma_totc = gamma_totc + 5.65e-36_prec*((1.0_prec+z)**7)*&
       &(2.726*(1.0_prec+z)-temph)*nhii*erg2j*yrbys*cmbympccb ! J Mpc^-3 yr^-1 

  nhi = nh*(1.0_prec-oldxiiva) 
  gammaph = (gph_pop2*source_pop2 + gph_pop3*source_pop3)*lmfp*&
       &fesc*(1.0_prec+z)**3*(cmbympc**2)/q  ! j/yr 
  gamma_heat = gammaph*nhi*(1.0_prec+z)**3 ! j mpc^-3 yr^-1 

  heatcool = 2.0_prec*(gamma_heat-gamma_totc)/&
       &(3.0_prec*nh_proper*(1.0_prec+x_iiva)*kboltz) ! k yr^-1 
  temphva = (temphva + dz*dtdz(z)*heatcool) / &
       &(1.0_prec + 2.0_prec*hubp(z)*dz*dtdz(z) + &
       &((x_iiva-oldxiiva)/(1.0_prec+x_iiva)))
!!$  temphva = (temphva + dz*dtdz(z)*heatcool) / &
!!$       &(1.0_prec + 2.0_prec*hubp(z)*dz*dtdz(z) + &
!!$       &((x_iiva-x2init)/(1.0_prec+x_iiva)))


