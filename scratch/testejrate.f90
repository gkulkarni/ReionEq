program testejrate

  ! File: testejrate.f90
  !  Cre: 2012-03-09
  !  Mod: $Date: 2012/03/09 08:17:40 $ ($Revision: 1.1 $) 
  !
  ! Tests ejrate.f90. Compiles with `make testejrate'.

  use constants
  use interfaces, only : ejrate 
  use storage 
  implicit none 

  real(kind=prec) :: z, met, fe, c, o, dtrans
  integer :: ncalc, i, j, number_of_lines 
  real(kind=prec) :: a, b, bound, dfm, dgrf, dlt, dlthi, dltlo, &
       &dsig, fmhi, fmlo, fv, f_m, gammaph, gamma_heat, dnlldz, &
       &gamma_ion, gamma_rc, gamma_rec, gamma_totc, grf, h0t0, heatcool, &
       &htt, igmdcrit, jnsln, lmfp, m, ngamma, nh, nhi, gammapi,&
       &nhii, nh_proper, oldfm, q, qe, r, r1, r2, reclambda, recrate, &
       &rho, rho_baryon, sig, source, sum, error, t0sec, t0yr, t1, t2, &
       &tempc, temph, temph_ev, tolzin, x2init, heat, cool, acc, &
       &zoft, freq, ndotm, gph, gpi, tau1, tau2, tau, jnsm, total_age,&
       &temphva, x_iiva, oldxiiva, lb_age, lb_lambda, lb_lum, initialtime,&
       &finaltime, halototalmass, haloformationredshift, t, haloluminosity,&
       &vcirc, dofz, source_pop2, source_pop3, gpi_pop3, gpi_pop2, ejc, &
       &gph_pop2, gph_pop3, st_mass, st_age, yield_record(10), ofl, &
       &ejr, local_hubble_0, metallicity, metmass, m_c, m_o, m_fe, &
       &m_h, fe_abundance, c_abundance, o_abundance, source_dummy, &
       &ofl_factor, ejrate_fe, ejrate_o, ejrate_c 

  include 'readin.inc' 

  ncalc = (final_redshift-initial_redshift)/dz+1 
  allocate(metarr(ncalc))

  open(unit=11, file="abundances.out", status="old", action="read")
  i = 1 
  do
     read(11, *, end=91) z, met, fe, c, o, dtrans
     metarr(i) = met 
     i = i+1 
  end do
91 continue 
  close(11)

  open(unit=11, file="ejrate.chk", status="unknown", action="write")
  z = initial_redshift 
  do 
     z = z + dz 
     if (z < final_redshift) exit 

     ! write (11,*) z, ejrate(z,0), ejrate(z,14), ejrate(z,7), ejrate(z,9) 
     write (11,*) z, ejrate(z,14)
  end do
  close(11)

!!$  do j = 1, 10
!!$     do i = 1, 4
!!$        print *, yields(i,1,1), yields(i,j,15)
!!$     end do
!!$     print *, '----'
!!$  end do

!!$  do i = 1, 4
!!$     print *, yields(i,1,1), yields(i,1,15)
!!$  end do

end program testejrate

