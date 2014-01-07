
PROGRAM REION 

  ! File: main.f90
  ! Cre: 2010-04-07
  ! Mod: $Date: 2012/10/23 08:10:11 $; $Revision: 1.45 $ 

  USE CONSTANTS 
  USE STORAGE 
  USE INTERFACES, ONLY : CLUMPFAC, DTDZ, GAMMA_PH, GAMMA_PI, HUBP, &
       &IGMFPREO, IGMVFRAC, JEANS_LENGTH, PSIG, RTBIS, SIGMAH , &
       &SIGMA_BARYON, SOLDLT, SOLFM, LUMFN, gallum, COUNTER, &
       &SOURCEOV_POP2, SOURCEOV_POP3, GPH_KERNEL_POP2, outflow, &
       &GPI_KERNEL_POP2, GPI_KERNEL_POP3, GPH_KERNEL_POP3, hallum, &
       &ACCRATE, ejrate, sfr_rollinde_pop3, sfr_rollinde_pop2, &
       &interpolate2, getsfr2, getsfr3, getjmc, getjmh, haloyield_nonira, &
       &rtnewt, newtondelta, newtonfm, igmfposto, ejfrac_nonira, outfrac_nonira, &
       &haloyield_species_nonira, counter, ngammafrac, haloyield_species, &
       &hallum2, bi_interpolate2, scale_sed, MajorMergerRate, dzdt 
  IMPLICIT NONE 

  REAL(KIND=PREC) :: A, B, BOUND, C, DFM, DGRF, DLT, DLTHI, DLTLO, &
       &DSIG, FMHI, FMLO, FV, F_M, GAMMAPH, GAMMA_HEAT, DNLLDZ, &
       &GAMMA_ION, GAMMA_RC, GAMMA_REC, GAMMA_TOTC, GRF, H0T0, HEATCOOL, &
       &HTT, IGMDCRIT, JNSLN, LMFP, M, NGAMMA, NH, NHI, GAMMAPI,&
       &NHII, NH_PROPER, OLDFM, Q, QE, R, R1, R2, RECLAMBDA, RECRATE, &
       &RHO, RHO_BARYON, SIG, SOURCE, SUM, ERROR, T0SEC, T0YR, T1, T2, &
       &TEMPC, TEMPH, TEMPH_EV, TOLZIN, X2INIT, Z, HEAT, COOL, acc, &
       &ZOFT, FREQ, NDOTM, GPH, GPI, TAU1, TAU2, TAU, JNSM, TOTAL_AGE,&
       &TEMPHVA, X_IIVA, OLDXIIVA, LB_AGE, LB_LAMBDA, LB_LUM, InitialTime,&
       &FinalTime, HaloTotalMass, HaloFormationRedshift, t, HaloLuminosity,&
       &VCIRC, DOFZ, SOURCE_POP2, SOURCE_POP3, GPI_POP3, GPI_POP2, ejc, &
       &GPH_POP2, GPH_POP3, ST_MASS, ST_AGE, YIELD_RECORD(25), ofl, &
       &ejr, LOCAL_HUBBLE_0, metallicity, metmass, m_c, m_o, m_fe, &
       &m_h, fe_abundance, c_abundance, o_abundance, source_dummy, &
       &dtrans, ofl_factor, ejrate_fe, ejrate_o, ejrate_c, febyh, &
       &dfebyh, zfe, psi, totalmstar, tfe, tnow, age_fe, mstar_fe, &
       &r_local, mhalo, mhalo_high, mhalo_low, mdot, mstardot, mmetaldot, &
       &return_fraction, fgas_in, mgasdot, mstardot_insitu, Omega_mz, &
       &fout, HaloVirialRadius, d, Delta_c, DiscSpinParameter, &
       &DiscScaleLength, HaloCircularVelocity, tdyn, zeta, mCdot, &
       &mFedot, mOdot, mH_halos, nu, f, grw, sgm, ms, n, dsgm, mo, &
       &log_multiplier, BinsPerDecade, dm, sm, smc, smh, mminc, mminh, &
       &sm_pop3, sm_pop2, smc_pop3, smc_pop2, smh_pop3, smh_pop2, &
       &st_ngamma, metal_term1, metal_term2, metal_term3, term1, &
       &term2, global_t, fb_struct, mNdot, mSidot, mZndot, mMgdot,&
       &ejrate_n, ejrate_si, ejrate_mg, ejrate_zn, foo, decr, mcoolgas,&
       &sdt, sdl, HaloVirialTemp, GasCoolingRate, eta, p, hb, fb, ft,&
       &mcooldot, halo_lum, zlim, maglim, limsfrc, limsfrh, limsfr,&
       &temp_sd93, lambda_sd93, t_test, fe_test, l_test, FirstStar_redshift, &
       &FirstStar_time, Zcr_redshift, Zcr_time, gammapi_pop2, gammapi_pop3, &
       &fesc_pop2, oldigmdcrit, oldlmfp, fake_ngamma, hy, hy_gradual, zold, &
       &ejrate_fe_gradual, ejrate_c_gradual, ejrate_o_gradual, ejrate_n_gradual, &
       &ejrate_si_gradual, ejrate_mg_gradual, ejrate_zn_gradual, ejrate_tot_gradual, &
       &burst_enhancement, starburst_duty_cycle 

  real(kind=prec) :: m_igm, m_ism, m_str, xigm_fe, xigm_c, xigm_o, &
       &xism_fe, xism_c, xism_o, dm_ism, dm_igm, dm_str, &
       &xism_tot, xigm_tot, ejrate_tot, m_totz, xigm_n, xigm_mg, xigm_si, &
       &xigm_zn, xism_n, xism_si, xism_mg, xism_zn

  REAL(KIND=PREC) :: MAGFOO, ZFOO, NFOO1, NFOO2, MSTAR, MIGM, MISM 

  INTEGER :: I, J, IER, INF, LAST, NEVAL, NUMBER_OF_LINES, NCALC, &
       &comnum, n_halocalc, glc, yfoo, lim_haloindexc, lim_haloindexh, jc, countr_init

  CHARACTER(100) :: ZLUMARG, RGNOVDARG, FESCARG, RGNSIZEARG, yfooarg 

  LOGICAL :: PREOVERLAP, limnotset_h, limnotset_c

  ! ----------------------------------------

  comnum = command_argument_count()
  if (comnum /= 4) then 
     write (0, '(a)') 'Usage: ./reion RGNOVD RGNSIZE ZLUM FESC'
     write (0, '(a)') 'Example: ./reion 4.0 7.746 6.5 0.3'
     stop
  end if

  CALL GET_COMMAND_ARGUMENT(1, RGNOVDARG) 
  CALL GET_COMMAND_ARGUMENT(2, RGNSIZEARG) 
  CALL GET_COMMAND_ARGUMENT(3, ZLUMARG)
  CALL GET_COMMAND_ARGUMENT(4, FESCARG)

  READ(RGNOVDARG, '(F10.5)') RGNOVD
  READ(RGNSIZEARG, '(F10.5)') RGNSIZE
  READ(ZLUMARG, '(F10.5)') ZLUM
  READ(FESCARG, '(F10.5)') FESC

  ! fesc_pop2 = 0.1
  fesc_pop2 = fesc 

  ! ----------------------------------------

  ! Set initial conditions. 
  INCLUDE 'readin.inc' 

  ! NGAMMA is calculated from Starburst99 model
  ! `reion-generic' by popsyn/ngtot.f90.
  NGAMMA = 2.82E70_PREC ! (10^10 M_solar)^-1 
  NGAMMA_POP2 = NGAMMA 
  ngamma_pop3 = ngammafrac(z)*1.0e10 ! (10^10 M_solar)^-1
  fake_ngamma = scale_sed()
  Pop3_SEDScaleFactor = Ngamma_pop3 / fake_ngamma

  Z = INITIAL_REDSHIFT 
  Q = 1.0E-8_PREC 

  ! Allocate and initialize arrays that will store values of Jeans
  ! mass and filling factor for redshifts at which calculation is 
  ! done.  Length of these arrays will have to be equal to the 
  ! number of such redshift values. 
  NCALC = (FINAL_REDSHIFT-INITIAL_REDSHIFT)/DZ+1 
  ALLOCATE(JMHARR(NCALC)); ALLOCATE(FFARR(NCALC)); ALLOCATE(SFRARR(NCALC))
  allocate(metarr(ncalc)); allocate(febyharr(ncalc)); allocate(zfearr(ncalc))
  allocate(sfrarr_pop2(ncalc)); allocate(sfrarr_pop3(ncalc)) 
  allocate(nuintegralarr_pop2(ncalc)); allocate(nuintegralarr_pop3(ncalc)) 

  JMHARR = 0.0_PREC; FFARR = 0.0_PREC; SFRARR = 0.0_PREC; zfearr = 0.0_prec
  metarr = 0.0_prec; febyharr = 0.0_prec
  sfrarr_pop2 = 0.0_prec; sfrarr_pop3 = 0.0_prec 
  nuintegralarr_pop2 = 0.0_prec; nuintegralarr_pop3 = 0.0_prec 
  COUNTR = 1 
  FFARR(COUNTR) = Q 
  zfearr(countr) = z

  RHO = OMEGA_NR*RHO_CRITICAL ! 10^10 M_solar / Mpc^3 
  RHO_BARYON = OMEGA_B*RHO_CRITICAL ! 10^10 M_solar / Mpc^3 

  ! Comoving hydrogen _number_ density; 1.1891E57 is the number of
  ! protons in a M_solar.  Note the 10^10 due to our mass units.
  DOFZ = GRFARR(COUNTER(Z))
  NH = RHO_BARYON*(1.0_PREC+RGNOVD*DOFZ)*(1.1891E57_PREC*1.0E10_PREC) ! Mpc^-3 
  NH_PROPER = NH*(1.0_PREC+Z)**3 ! Mpc^-3 

  ! TEMPC is temperature of neutral gas at redshift 
  ! INITIAL_REDSHIFT.  Here I have estimated it by 
  ! assuming that CMB and the neutral gas were coupled
  ! till z = 200 after which the gas cooled ~ (1 + z)^2 
  ! while CMB cooled ~ (1 + z).  See Bharadwaj & Ali (2004).
  TEMPC = CMBTEMP*(1.0_PREC+Z)**2 / 201.0_PREC ! K 
  TEMPH = TEMPC ! K

  TEMPHVA = TEMPH ! K 

  JNSLN = JEANS_LENGTH(TEMPHVA,Z) ! Mpc 

  JNSM = (4.0_PREC*PI*RHO_BARYON*JNSLN**3)/3.0_PREC ! 10^10 M_solar
  JMHARR(COUNTR) = JNSM ! 10^10 M_solar 

  ! Normalise matter power spectrum. 
  RCARRY_PSIG = 8.0_PREC/SMALLH ! Mpc 
  BOUND = 0.0_PREC
  INF = 1 
  CALL DQAGIE(PSIG,BOUND,INF,ABSERR,ABSREL,MAXINT,SUM,ERROR,&
       &NEVAL,IER,ALIST,BLIST,RLIST,ELIST,IORD,LAST)
  IF (IER > 0) WRITE (0,*) 'PSPECNORM: Error.' 
  PSPECNORM = (SIGMA_EIGHT**2/SUM)

  ! Calculate the rms linear mass fluctuation in baryons.
  SIGMAB = SIGMA_BARYON(Z, TEMPH) 

  IGMDCRIT = IGMDCRIT_PREO 
  R = CLUMPFAC(IGMDCRIT)
  r_local = r 
  F_M = IGMFPREO()
  FV = IGMVFRAC(IGMDCRIT) 
  LMFP = q**(1.0_prec/3.0_prec)*LMFP0*JNSLN/((1.0_PREC-FV)**(2.0_PREC/3.0_PREC)) ! Mpc 

  X2INIT = 1.2E-5_PREC/(SMALLH*OMEGA_B) ! dimensionless 
  X_II = X2INIT
  OLDXII = X_II 

  X_IIVA = X_II
  OLDXIIVA = OLDXII 

  PREOVERLAP = .TRUE. 
  POSTOCOUNTER = 1 

  GPH_POP2 = GPH_KERNEL_POP2() 
  GPI_POP2 = GPI_KERNEL_POP2()
  GPH_POP3 = GPH_KERNEL_POP3() 
  GPI_POP3 = GPI_KERNEL_POP3()
  TAU = 0.0_PREC 
  TAU1 = TAUINT(FV, Q, X_IIVA, Z) 

  m_str = 0.0_prec 
  m_ism = 1.0e-4_prec
  m_igm = rho_baryon 

  xism_fe = 0.0_prec 
  xism_c = 0.0_prec 
  xism_o = 0.0_prec 
  xism_n = 0.0_prec 
  xism_si = 0.0_prec 
  xism_mg = 0.0_prec 
  xism_zn = 0.0_prec 
  xism_tot = 0.0_prec 

  xigm_fe = 0.0_prec 
  xigm_c = 0.0_prec 
  xigm_o = 0.0_prec 
  xigm_n = 0.0_prec 
  xigm_si = 0.0_prec 
  xigm_zn = 0.0_prec 
  xigm_mg = 0.0_prec 
  xigm_tot = 0.0_prec 

  metarr(1) = 0.0_prec 

  open(unit=31, file='abundances.out', status='unknown', action='write')
  open(unit=32, file='masses.out', status='unknown', action='write')
  open(unit=33, file='rates.out', status='unknown', action='write')
  open(unit=34, file='sfr.out', status='unknown', action='write')
  open(unit=35, file='fracs.out', status='unknown', action='write')
  open(unit=36, file='ejrates.out', status='unknown', action='write')
  open(unit=37, file='mdf.out', status='unknown', action='write')
  open(unit=38, file='reion.out', status='unknown', action='write')
  open(unit=39, file='halos.out', status='unknown', action='write') 
  open(unit=40, file='halos_stars.out', status='unknown', action='write') 
  open(unit=41, file='halos_gas.out', status='unknown', action='write') 
  open(unit=42, file='halos_metals.out', status='unknown', action='write') 
  open(unit=43, file='halos_c.out', status='unknown', action='write') 
  open(unit=44, file='halos_fe.out', status='unknown', action='write') 
  open(unit=45, file='halos_o.out', status='unknown', action='write') 
  open(unit=46, file='halos_zism.out', status='unknown', action='write') 
  open(unit=47, file='mmin.out', status='unknown', action='write')
  open(unit=48, file='strpop.out', status='unknown', action='write')
  open(unit=49, file='halos_aux.out', status='unknown', action='write') 
  open(unit=50, file='avghalos.out', status='unknown', action='write') 
  open(unit=51, file='avghalos_stars.out', status='unknown', action='write') 
  open(unit=52, file='avghalos_gas.out', status='unknown', action='write') 
  open(unit=53, file='avghalos_metals.out', status='unknown', action='write') 
  open(unit=54, file='avghalos_c.out', status='unknown', action='write') 
  open(unit=55, file='avghalos_fe.out', status='unknown', action='write') 
  open(unit=56, file='avghalos_o.out', status='unknown', action='write') 
  open(unit=57, file='avghalos_zism.out', status='unknown', action='write') 
  open(unit=58, file='params.out', status='unknown', action='write')
  open(unit=59, file='halos_N.out', status='unknown', action='write') 
  open(unit=60, file='halos_Si.out', status='unknown', action='write') 
  open(unit=61, file='halos_Zn.out', status='unknown', action='write') 
  open(unit=62, file='avghalos_N.out', status='unknown', action='write') 
  open(unit=63, file='avghalos_Si.out', status='unknown', action='write') 
  open(unit=64, file='avghalos_Zn.out', status='unknown', action='write') 
  open(unit=65, file='nofmc.out', status='unknown', action='write') 
  open(unit=66, file='nofmh.out', status='unknown', action='write') 
  open(unit=67, file='halos_Mg.out', status='unknown', action='write') 
  open(unit=68, file='avghalos_Mg.out', status='unknown', action='write')
  open(unit=73, file='l1500.out', status='unknown', action='write') 
  open(unit=74, file='halos_coolgas.out', status='unknown', action='write') 
  open(unit=75, file='avghalos_coolgas.out', status='unknown', action='write') 
  open(unit=76, file='halosh_sfrcontrib.out', status='unknown', action='write') 
  open(unit=80, file='halosc_sfrcontrib.out', status='unknown', action='write') 
  open(unit=78, file='avghalos_sfrcontrib.out', status='unknown', action='write') 
  open(unit=77, file='halos_mstardot.out', status='unknown', action='write') 
  open(unit=79, file='avghalos_mstardot.out', status='unknown', action='write') 
  open(unit=81, file='halos_Zcrtimehot.out', status='unknown', action='write') 
  open(unit=82, file='halos_Zcrtimecold.out', status='unknown', action='write') 
  open(unit=83, file='halos_Firststartimehot.out', status='unknown', action='write') 
  open(unit=84, file='halos_Firststartimecold.out', status='unknown', action='write') 

  open(unit=69, file='aux_fe.out', status='unknown', action='write') 
  open(unit=70, file='aux_o.out', status='unknown', action='write') 
  open(unit=71, file='aux_si.out', status='unknown', action='write') 
  open(unit=72, file='aux_zn.out', status='unknown', action='write') 

  !------------------------------

  mhalo_high = 0.5e-3_prec ! 10^10 M_solar 
  mhalo_low = 1.0e-6_prec ! 10^10 M_solar 

  BinsPerDecade = 100.0_prec
  log_multiplier = 10.0_prec**(1.0_prec/BinsPerDecade)
  n_halocalc = int(log10(mhalo_high/mhalo_low)/log10(log_multiplier))+1 

  allocate(m_halosc(n_halocalc))
  allocate(mstar_halosc(n_halocalc))
  allocate(mstardot_halosc(n_halocalc))
  allocate(mgas_halosc(n_halocalc)) 
  allocate(mcoolgas_halosc(n_halocalc)) 
  allocate(mmetal_halosc(n_halocalc)) 
  allocate(mC_halosc(n_halocalc)) 
  allocate(mFe_halosc(n_halocalc)) 
  allocate(mO_halosc(n_halocalc)) 
  allocate(mN_halosc(n_halocalc)) 
  allocate(mSi_halosc(n_halocalc)) 
  allocate(mZn_halosc(n_halocalc)) 
  allocate(mMg_halosc(n_halocalc)) 
  allocate(febyh_halosc(n_halocalc)) 
  allocate(cbyh_halosc(n_halocalc)) 
  allocate(obyh_halosc(n_halocalc)) 
  allocate(strpop_halosc(n_halocalc))
  allocate(ismz_halosc(n_halocalc)) 
  allocate(aux_halosc(n_halocalc)) 
  allocate(sfrcontrib_halosc(n_halocalc)) 
  allocate(NotEnrichedYet_cold(n_halocalc))
  allocate(FirstStar_timec(n_halocalc))
  allocate(Zcr_timec(n_halocalc)) 

  allocate(m_halosh(n_halocalc))
  allocate(mstar_halosh(n_halocalc))
  allocate(mstardot_halosh(n_halocalc))
  allocate(mgas_halosh(n_halocalc)) 
  allocate(mcoolgas_halosh(n_halocalc)) 
  allocate(mmetal_halosh(n_halocalc)) 
  allocate(mC_halosh(n_halocalc)) 
  allocate(mFe_halosh(n_halocalc)) 
  allocate(mO_halosh(n_halocalc)) 
  allocate(mN_halosh(n_halocalc)) 
  allocate(mSi_halosh(n_halocalc)) 
  allocate(mZn_halosh(n_halocalc)) 
  allocate(mMg_halosh(n_halocalc)) 
  allocate(febyh_halosh(n_halocalc)) 
  allocate(cbyh_halosh(n_halocalc)) 
  allocate(obyh_halosh(n_halocalc)) 
  allocate(strpop_halosh(n_halocalc))
  allocate(ismz_halosh(n_halocalc)) 
  allocate(aux_halosh(n_halocalc)) 
  allocate(nofmc(n_halocalc))
  allocate(nofmh(n_halocalc)) 
  allocate(halol1500c(n_halocalc))
  allocate(halol1500h(n_halocalc))
  allocate(sfrcontrib_halosh(n_halocalc)) 
  allocate(NotEnrichedYet_hot(n_halocalc))
  allocate(FirstStar_timeh(n_halocalc))
  allocate(Zcr_timeh(n_halocalc)) 

  allocate(tdyn_array(n_halocalc)) 

  febyh_halosc = -10.0_prec
  febyh_halosh = -10.0_prec

  allocate(t_zn(n_halocalc), t_fe(n_halocalc), t_si(n_halocalc), t_o(n_halocalc))

  m = mhalo_low 
  do i = 1, n_halocalc 
     m_halosc(i) = m; m_halosh(i) = m 
     m = m*log_multiplier 

     strpop_halosc(i) = 3; strpop_halosh(i) = 3 

     mmetal_halosc(i) = 0.0_prec; mmetal_halosh(i) = 0.0_prec 
     mgas_halosc(i) = 0.0_prec; mgas_halosh(i) = 0.0_prec 
     mcoolgas_halosc(i) = 0.0_prec; mcoolgas_halosh(i) = 0.0_prec 
     mstar_halosc(i) = 0.0_prec; mstar_halosh(i) = 0.0_prec 
     mstardot_halosc(i) = 0.0_prec; mstardot_halosh(i) = 0.0_prec 
     mC_halosc(i) = 0.0_prec; mC_halosh(i) = 0.0_prec 
     mFe_halosc(i) = 0.0_prec; mFe_halosh(i) = 0.0_prec 
     mO_halosc(i) = 0.0_prec; mO_halosh(i) = 0.0_prec 
     mN_halosc(i) = 0.0_prec; mN_halosh(i) = 0.0_prec 
     mSi_halosc(i) = 0.0_prec; mSi_halosh(i) = 0.0_prec 
     mZn_halosc(i) = 0.0_prec; mZn_halosh(i) = 0.0_prec 
     mMg_halosc(i) = 0.0_prec; mMg_halosh(i) = 0.0_prec 
     ismz_halosc(i) = 0.0_prec; ismz_halosh(i) = 0.0_prec 
     nofmc(i) = 0.0_prec; nofmh(i) = 0.0_prec 
     sfrcontrib_halosc(i) = 0.0_prec; sfrcontrib_halosh(i) = 0.0_prec 

     NotEnrichedYet_cold(i) = .true.; NotEnrichedYet_hot(i) = .true. 
     FirstStar_timec(i) = 0.0_prec; FirstStar_timeh(i) = 0.0_prec 
     Zcr_timec(i) = 0.0_prec; Zcr_timeh(i) = 0.0_prec

     tdyn_array(i) = 0.0_prec 

     t_zn(i) = 0.0_prec 
     t_si(i) = 0.0_prec 
     t_fe(i) = 0.0_prec 
     t_o(i) = 0.0_prec 

     aux_halosc(i) = 0.0_prec 
  end do

  allocate(sfrarr_halocalc_cold(ncalc, n_halocalc))
  allocate(sfrarr_halocalc_hot(ncalc, n_halocalc))
  allocate(sfrarr_halocalc(ncalc, n_halocalc)) ! Used in hallum.f90. 
  allocate(halopop_hot(ncalc, n_halocalc))
  allocate(halopop_cold(ncalc, n_halocalc))

  ! The following four arrays are for implementing delay in
  ! enrichment.
  
  ! In ejrate_array(dim1, dim2), dim2 tracks the various species
  ! according to this mapping (where the numbers in the first column
  ! are defined in ejrate.f90):
  !
  ! 0 -- 1 
  ! 14 -- 2
  ! 7 -- 3 
  ! 9 -- 4 
  ! 8 -- 5 
  ! 11 -- 6 
  ! 10 -- 7 
  ! 15 -- 8 
  ! 16 -- 9 

  allocate(ejrate_array(ncalc, 9)) 
  ejrate_array = 0.0_prec 

  ! The following two arrays are 3-dimensional.  There structure is
  ! as follows: array(dim1, dim2, dim3), where dim1 tracks redshift
  ! steps of code, dim2 tracks halo mass bins, and dim3 tracks hot and
  ! cold regions (dim3=1: cold; dim3=2: hot). 

  allocate(haloyield_nonira_array(ncalc, n_halocalc, 2)) 
  allocate(ejfrac_nonira_array(ncalc, n_halocalc, 2)) 

  haloyield_nonira_array = 0.0_prec 
  ejfrac_nonira_array = 0.0_prec
  
  ! In haloyield_species_nonira_array(dim1, dim2, dim3, dim4), dim3
  ! tracks various species.  Other dimensions are same as
  ! haloyield_nonira.  The species mapping is as follows, where the
  ! numbers in the first column are defined in
  ! haloyield_species_nonira.f90. 
  ! 
  ! 7 -- 1 
  ! 14 -- 2 
  ! 9 -- 3
  ! 8 -- 4 
  ! 11 -- 5 
  ! 15 -- 6 
  ! 10 -- 7 

  allocate(haloyield_species_nonira_array(ncalc, n_halocalc, 7, 2)) 

  haloyield_species_nonira_array = 0.0_prec 

  allocate(countr_old(ncalc))
  countr_old = 0 

  do i = 1, ncalc 
     zold = i*dz + initial_redshift 
     call interpolate2(tarr, zarr, zold, t) ! [t] = yr 
     t = t - Enrich_time_lag2
     call interpolate2(zarr, tarr, t, zold) 
     countr_old(i) = int((zold-initial_redshift)/dz)
  end do

  Zlim = 10.0_prec 
  maglim = -18.0_prec 
  eta = 1.65_prec ! halo density profile index used in mcooldot below.
  starburst_duty_cycle = 0.1 ! See plots/compare_sfr.pro.
  
  write (58,*) 'initial_redshift=', initial_redshift 
  write (58,*) 'final_redshift=', final_redshift 
  write (58,*) 'dz=', dz 
  write (58,*) 'fstar=', fstar
  write (58,*) 'fstar_pop3=', fstar_pop3
  write (58,*) 'fesc=', fesc 
  write (58,*) 'fesc_pop2=', fesc_pop2 
  write (58,*) 'pop2 imf:', minf, msup 
  write (58,*) 'pop3 imf:', minf_pop3, msup_pop3
  write (58,*) 'imf_slope', imf_slope
  write (58,*) 'data_mmin=', data_mmin
  write (58,*) 'data_mmax=', data_mmax 
  write (58,*) 'stmass_uplimit=', stmass_uplimit 
  write (58,*) 'eta=', eta
  write (58,*) 'Enrich_time_lag2=', Enrich_time_lag2


  do 
     countr = countr + 1 
     z = z + dz 
     if (z < final_redshift) exit 

     write (0,*) z
     call flush(0) 
     

     smc = 0.0_prec; smc_pop3 = 0.0_prec; smc_pop2 = 0.0_prec  
     smh = 0.0_prec; smh_pop3 = 0.0_prec; smh_pop2 = 0.0_prec 
     limsfr = 0.0_prec 
     lim_haloindexc = 0 
     lim_haloindexh = 0 

     limnotset_h = .true.
     limnotset_c = .true. 

     ! Calculate magnitude cut-off corresponding to Bouwens. 
     if (z <= zlim) then 
        do i = 2, n_halocalc
           halo_lum = hallum2(i,countr-1,1) ! erg/s/Hz 
           halol1500c(i) = -2.5_prec*log10(halo_lum)+51.60_prec
           if ((halol1500c(i)<maglim) .and. limnotset_c) then 
              lim_haloindexc = i-1
              limnotset_c = .false.
           end if
        end do
     end if

     if (z <= zlim) then 
        do i = 2, n_halocalc
           halo_lum = hallum2(i,countr-1,2) ! erg/s/Hz 
           halol1500h(i) = -2.5_prec*log10(halo_lum)+51.60_prec
           if ((halol1500h(i)<maglim) .and. limnotset_h) then 
              lim_haloindexh = i-1
              limnotset_h = .false.
           end if
        end do
     end if

     ! New source calculation 
     do i = 2, n_halocalc

        call interpolate2(sigmarr, msarr, m_halosc(i), sgm) ! [sgm] = dimensionless 
        call interpolate2(dsigmarr, msarr, m_halosc(i), dsgm) ! [dsgm] = (10^10 M_solar)^-1 
        glc = counter(z)
        grw = grfarr(glc) ! dimensionless 
        nu = deltac/(grw*sgm) ! dimensionless 
        f = sqrt(2.0_prec/pi)*exp(-nu*nu/2.0_prec) ! dimensionless 
        dm = m_halosc(i) - m_halosc(i-1) ! 10^10 M_solar
        n = (rho/m_halosc(i))*f*dsgm*(-deltac/(grw*sgm*sgm))*dm ! Mpc^-3 

        nofmc(i) = n 

        if (strpop_halosc(i) == 3) then 
           smc_pop3 = smc_pop3 + mstardot_halosc(i)*n ! 10^10 M_solar yr^-1 Mpc^-3 
        else 
           smc_pop2 = smc_pop2 + mstardot_halosc(i)*n ! 10^10 M_solar yr^-1 Mpc^-3 
        end if

        sfrcontrib_halosc(i) = mstardot_halosc(i)*n

        if (z<=zlim .and. i>=lim_haloindexc) then 
           limsfrc = limsfrc + mstardot_halosc(i)*n ! 10^10 M_solar yr^-1 Mpc^-3 
        else 
           limsfrc = 0.0_prec 
        end if

        call interpolate2(sigmarr, msarr, m_halosh(i), sgm) ! [sgm] = dimensionless 
        call interpolate2(dsigmarr, msarr, m_halosh(i), dsgm) ! [dsgm] = (10^10 M_solar)^-1 
        glc = counter(z)
        grw = grfarr(glc) ! dimensionless 
        nu = deltac/(grw*sgm) ! dimensionless 
        f = sqrt(2.0_prec/pi)*exp(-nu*nu/2.0_prec) ! dimensionless 
        dm = m_halosh(i) - m_halosh(i-1) ! 10^10 M_solar
        n = (rho/m_halosh(i))*f*dsgm*(-deltac/(grw*sgm*sgm))*dm ! Mpc^-3 

        nofmh(i) = n 

        if (strpop_halosc(i) == 3) then 
           smh_pop3 = smh_pop3 + mstardot_halosh(i)*n ! 10^10 M_solar yr^-1 Mpc^-3 
        else 
           smh_pop2 = smh_pop2 + mstardot_halosh(i)*n ! 10^10 M_solar yr^-1 Mpc^-3 
        end if

        sfrcontrib_halosh(i) = mstardot_halosh(i)*n 

        if (z<=zlim .and. i>=lim_haloindexh) then 
           limsfrh = limsfrh + mstardot_halosh(i)*n ! 10^10 M_solar yr^-1 Mpc^-3 
        else 
           limsfrh = 0.0_prec 
        end if
     end do

     smh = smh_pop3 + smh_pop2 
     smc = smc_pop3 + smc_pop2 
     sm = q*smh + (1.0_prec-q)*smc
     sm_pop3 = q*smh_pop3 + (1.0_prec-q)*smc_pop3
     sm_pop2 = q*smh_pop2 + (1.0_prec-q)*smc_pop2

     limsfrc = limsfrc * 1.0e10_prec ! M_solar yr^-1 Mpc^-3 
     limsfrh = limsfrh * 1.0e10_prec ! M_solar yr^-1 Mpc^-3 
     limsfr = q*limsfrh + (1.0_prec-q)*limsfrc

     !-------------------------

     ! Calculate rate of ionizing photons (nphdot). 
     source = sm 
     source_pop2 = sm_pop2
     source_pop3 = sm_pop3

     if (source_pop3 .ne. source_pop3) source_pop3 = 0.0_prec 
     
     sfrarr(countr-1) = source*1.0e10_prec ! M_solar yr^-1 Mpc^-3
     sfrarr_pop2(countr-1) = source_pop2*1.0e10_prec ! M_solar yr^-1 Mpc^-3
     sfrarr_pop3(countr-1) = source_pop3*1.0e10_prec ! M_solar yr^-1 Mpc^-3
     write (34, '(F4.1,4E11.3E2)') z, sfrarr(countr-1), sfrarr_pop2(countr-1), &
          &sfrarr_pop3(countr-1), limsfr 
     ! nphdot = fesc*(source_pop2*ngamma_pop2 + source_pop3*ngamma_pop3) ! yr^-1 Mpc^-3
     nphdot = fesc_pop2*source_pop2*ngamma_pop2 + fesc*source_pop3*ngamma_pop3 ! yr^-1 Mpc^-3

     !-------------------------

     ! Calculate case-A HII recombination coefficient (alpha_r).  This
     ! expression is from Hui and Gnedin (1997).
     reclambda = 2.0_prec*thi/temph 
     alpha_r = 1.269e-13_prec * (reclambda**1.503_prec) / &
          (1.0_prec + (reclambda/0.522_prec)**0.470_prec)**1.923_prec ! cm^3 s^-1 
     ! Add case-B recombination coeff. to that. 
     alpha_r = alpha_r + 2.753e-14_prec * (reclambda**1.5_prec) / &
          (1.0_prec + (reclambda/2.74_prec)**0.407_prec)**2.242_prec ! cm^3 s^-1 

     !-------------------------

     ! Calculate hydrogen number density (comoving and proper).  The
     ! conversion factor of 1.18e57 converts M_solar to number of
     ! hydrogen atoms. 1.0e10 reduces mass units to M_solar.
     dofz = grfarr(counter(z))
     nh = rho_baryon*(1.0_prec+rgnovd*dofz)*(1.1891e57_prec*1.0e10_prec) ! Mpc^-3 
     nh_proper = nh*(1.0_prec+z)**3 ! Mpc^-3      

     !-------------------------

     ! Update ionized hydrogen fraction (x_ii) in ionized region. 
     gamma_rec = r*nh_proper*alpha_r*cmbympccb*yrbys ! yr^-1 
!!$     gamma_ion = (gpi_pop2*source_pop2+gpi_pop3*source_pop3)*&
!!$          &fesc*lmfp*(1.0_prec+z)**3*(cmbympc**2)/q ! yr^-1 
     gamma_ion = (fesc_pop2*gpi_pop2*source_pop2+fesc*gpi_pop3*source_pop3)*&
          &lmfp*(1.0_prec+z)**3*(cmbympc**2)/q ! yr^-1 

     a = gamma_rec*dz*dtdz(z)
     b = 1.0_prec + dz*dtdz(z)*gamma_ion
     c = -(oldxii + dz*dtdz(z)*gamma_ion)

     qe = -0.5_prec*(b+sign(dsqrt(b**2-4.0_prec*a*c),b))
     r1 = qe/a; r2 = c/qe 

     oldxii = x_ii 
     if (r1<=1.0_prec .and. r1>=0.0_prec) then
        x_ii = r1
     else 
        x_ii = r2 
     end if

     !-------------------------

     ! Calculate case-A HII recombination cooling rate (recrate).
     ! This expression is from Hui and Gnedin (1997).
     reclambda = 2.0_prec*thi/temph 
     recrate = 1.778e-29_prec*temph*reclambda**1.965_prec/&
          &(1.0_prec+(reclambda/0.541_prec)**0.502_prec)&
          &**2.697_prec ! erg cm^3 s^-1
     ! Add case-B HII recombination cooling rate. 
     recrate = recrate + 3.435e-30_prec*temph*reclambda**1.970d0 / &
          (1.0_prec+(reclambda/2.250e0_prec)**0.376_prec)**3.720_prec ! erg cm^3 s^-1

     ! Update temperature of ionized regions (temph) and neutral
     ! regions (tempc).
     nhii = nh*oldxii ! Mpc^-3 
     gamma_rc = (nhii**2)*recrate*erg2j*cmbympccb*yrbys ! J Mpc^-3 yr^-1
     gamma_totc = r_local*gamma_rc*(1.0_prec+z)**6 ! J Mpc^-3 yr^-1

     ! Add Compton cooling 
     gamma_totc = gamma_totc + 5.65e-36_prec*((1.0_prec+z)**7)*&
          &(2.726*(1.0_prec+z)-temph)*nhii*erg2j*yrbys*cmbympccb ! J Mpc^-3 yr^-1 

     nhi = nh*(1.0_prec-oldxii) 
!!$     gammaph = (gph_pop2*source_pop2+gph_pop3*source_pop3)*fesc*&
!!$          &lmfp*(1.0_prec+z)**3*(cmbympc**2)  ! J/yr 
     gammaph = (fesc_pop2*gph_pop2*source_pop2+fesc*gph_pop3*source_pop3)*&
          &lmfp*(1.0_prec+z)**3*(cmbympc**2)  ! J/yr 

     gamma_heat = gammaph*nhi*(1.0_prec+z)**3/q ! J mpc^-3 yr^-1 

     heatcool = 2.0_prec*(gamma_heat-gamma_totc)/&
          &(3.0_prec*nh_proper*(1.0_prec+x_ii)*kboltz) ! K yr^-1 
     temph = (temph + dz*dtdz(z)*heatcool) / &
          &(1.0_prec + 2.0_prec*hubp(z)*dz*dtdz(z)+ &
          &((x_ii-oldxii)/(1.0_prec+oldxii)))
     tempc = tempc+dz*dtempcdz(z) ! K

     mu_MeanMolWt = 1.0_prec / (1.0_prec + x_ii)

     term1 = dtdz(z)*heatcool
     term2 = 2.0_prec*hubp(z)*dz*dtdz(z)+ ((x_ii-oldxii)/(1.0_prec+x_ii))

     !------------------------- 

     ! Minimum circular velocity of haloes that can cool in ionized
     ! regions.
     vcirc = sqrt(2.0_prec*kboltz*temph/mproton)

     ! Calculate corresponding Jeans mass. 
     jnsln = jeans_length(temphva,z) ! Mpc 
     global_t = q*temph+(1.0_prec-q)*tempc
     sigmab = sigma_baryon(z, global_t) 
     jnsm = (4.0_prec*pi*rho_baryon*jnsln**3)/3.0_prec ! 10^10 M_solar
     jmharr(countr) = jnsm

     !------------------------- 

     ! Update filling factor for the pre-overlap phase.
     if (preoverlap) then 
        r = clumpfac(igmdcrit)
        r_local = r 

        oldfm = f_m 
        f_m = igmfpreo()
        fv = igmvfrac(igmdcrit) 
        lmfp = q**(1.0_prec/3.0_prec)*lmfp0*jnsln/((1.0_prec-fv)**(2.0_prec/3.0_prec)) ! Mpc 

        dfm = f_m - oldfm 
        t1 = dz*dtdz(z)*nphdot/(nh*f_m)
        t2 = dz*dtdz(z)*alpha_r*yrbys*cmbympccb*nh*x_ii*r*(1.0_prec+z)**3/f_m
        t2 = t2+dfm/f_m
        q = (q+t1)/(1.0_prec+t2)

        ffarr(countr) = q 
     end if

     !-------------------------

     ! Update filling factor (i.e. F_M) for the post-overlap phase.
     if (q>=1.0_prec) then 
        preoverlap = .false. 
        postocounter = postocounter + 1 

        solfm_z = z 
        solfm_igmdcrit = igmdcrit 
        solfm_ff = q 
        ! fmlo = f_m 
        fmlo = 0.0
        fmhi = 0.99999999999_prec 
        tolzin = 1.0e-15_prec 
        oldigmfrac = f_m
        f_m = rtbis(solfm, fmlo, fmhi, tolzin, 'main-fm') 

        if (f_m < oldigmfrac) f_m = oldigmfrac 

        if (f_m .ne. f_m) f_m = oldigmfrac 

        q = 1.0_prec 
        ffarr(countr) = q 
        soldlt_z = z
        soldlt_igmfrac = f_m 
        dlthi = 1.0e6_prec !igmdcrit+1.0e30_prec
        dltlo = 0.0_prec !igmdcrit-59.0_prec
        tolzin = 1.0e-8_prec
        oldigmdcrit = igmdcrit 
        dlt = rtbis(soldlt, dltlo, dlthi, tolzin, 'main-dlt')
        igmdcrit = dlt 

        if (dlt < oldigmdcrit) dlt = oldigmdcrit 

        if (dlt .ne. dlt) dlt = oldigmdcrit 

        r = clumpfac(igmdcrit)
        fv = igmvfrac(igmdcrit) 
        oldlmfp = lmfp 
        lmfp = q**(1.0_prec/3.0_prec)*lmfp0*jnsln/((1.0_prec-fv)**(2.0_prec/3.0_prec)) ! mpc 
        if (lmfp .ne. lmfp) lmfp = oldlmfp 
     end if

     !-------------------------

     ! Calculate *mean* values of TEMPH and X_II.  (Above values are
     ! *mass averaged*, not mean.)
     include 'volav.inc'

     tau2 = tauint(fv, q, x_iiva, z) 
     tau = tau + 0.5_prec*dz*(tau1+tau2) 
     tau1 = tau2 

     dnlldz = speed_of_light*cmbympc*yrbys/&
          &(sqrt(pi)*lmfp*hubp(z)*(1.0_prec+z)) ! dimensionless 
!!$     gammapi = (gpi_pop2*source_pop2 + gpi_pop3*source_pop3)*fesc*lmfp*&
!!$          &(1.0_prec+z)**3*(cmbympc**2)/yrbys ! s^-1

     gammapi = (fesc_pop2*gpi_pop2*source_pop2 + fesc*gpi_pop3*source_pop3)*lmfp*&
          &(1.0_prec+z)**3*(cmbympc**2)/yrbys ! s^-1

     gammapi_pop2 = gpi_pop2*source_pop2*fesc_pop2*lmfp*(1.0_prec+z)**3*(cmbympc**2)/yrbys ! s^-1
     gammapi_pop3 = gpi_pop3*source_pop3*fesc*lmfp*(1.0_prec+z)**3*(cmbympc**2)/yrbys ! s^-1
     

     !-------------------------

     write (38,'(F4.1,16E11.3E2)') z, q, tau, gammapi, gammapi_pop2, gammapi_pop3, temph, tempc, &
          &(q*temph+(1.0_prec-q)*tempc), x_ii, dnlldz, lmfp, r, igmdcrit, nphdot, &
          &temphva, fv

     !-----------------------------

     ! Evolve various metal species. 

     ofl = outflow(z) ! 10^10 M_solar yr^-1 Mpc^-3 
     acc = accrate(z) ! 10^10 M_solar yr^-1 Mpc^-3
     ejc = ejrate(z,0) ! 10^10 M_solar yr^-1 Mpc^-3 

     ejrate_array(countr-1, 1) = ejc 

     dm_igm = (ofl - acc)*dz*dtdz(z) ! 10^10 M_solar Mpc^-3 
     dm_str = (source - ejc)*dz*dtdz(z) ! 10^10 M_solar Mpc^-3 
     dm_ism = - dm_igm - dm_str ! 10^10 M_solar Mpc^-3 

     m_igm = m_igm + dm_igm ! 10^10 M_solar Mpc^-3 
     m_str = m_str + dm_str ! 10^10 M_solar Mpc^-3 
     m_ism = m_ism + dm_ism ! 10^10 M_solar Mpc^-3 

     fb_struct = (m_str + m_ism ) / (m_igm + m_str + m_ism) 

     xigm_fe = xigm_fe + dz*dtdz(z)*ofl*(xism_fe-xigm_fe)/m_igm ! dimensionless
     xigm_c = xigm_c + dz*dtdz(z)*ofl*(xism_c-xigm_c)/m_igm ! dimensionless
     xigm_o = xigm_o + dz*dtdz(z)*ofl*(xism_o-xigm_o)/m_igm ! dimensionless
     xigm_n = xigm_n + dz*dtdz(z)*ofl*(xism_n-xigm_n)/m_igm ! dimensionless
     xigm_si = xigm_si + dz*dtdz(z)*ofl*(xism_si-xigm_si)/m_igm ! dimensionless
     xigm_mg = xigm_mg + dz*dtdz(z)*ofl*(xism_mg-xigm_mg)/m_igm ! dimensionless
     xigm_zn = xigm_zn + dz*dtdz(z)*ofl*(xism_zn-xigm_zn)/m_igm ! dimensionless
     xigm_tot = xigm_tot + dz*dtdz(z)*ofl*(xism_tot-xigm_tot)/m_igm ! dimensionless

     ! Below, 16=tot, 14=Fe, 7=C, and 9=O.  These are row numbers from
     ! the data tables of yields.  See comment in ejrate.f90.  Note:
     ! All four quantities calculated below are DIMENSIONLESS.

     ejrate_fe = ejrate(z,14) ! 10^10 M_solar yr^-1 Mpc^-3 
     ejrate_array(countr-1, 2) = ejrate_fe
     ejrate_c = ejrate(z,7) ! 10^10 M_solar yr^-1 Mpc^-3 
     ejrate_array(countr-1, 3) = ejrate_c
     ejrate_o = ejrate(z,9) ! 10^10 M_solar yr^-1 Mpc^-3 
     ejrate_array(countr-1, 4) = ejrate_o
     ejrate_n = ejrate(z,8) ! 10^10 M_solar yr^-1 Mpc^-3 
     ejrate_array(countr-1, 5) = ejrate_n
     ejrate_si = ejrate(z,11) ! 10^10 M_solar yr^-1 Mpc^-3 
     ejrate_array(countr-1, 6) = ejrate_si
     ejrate_mg = ejrate(z,10) ! 10^10 M_solar yr^-1 Mpc^-3 
     ejrate_array(countr-1, 7) = ejrate_mg
     ejrate_zn = ejrate(z,15) ! 10^10 M_solar yr^-1 Mpc^-3 
     ejrate_array(countr-1, 8) = ejrate_zn
     ejrate_tot = ejrate(z,16) ! 10^10 M_solar yr^-1 Mpc^-3 
     ejrate_array(countr-1, 9) = ejrate_tot

     if (Enrich_time_lag2 > 0.0_prec) then 
        
        ejrate_fe_gradual = 0.0_prec 
        ejrate_c_gradual = 0.0_prec 
        ejrate_o_gradual = 0.0_prec 
        ejrate_n_gradual = 0.0_prec 
        ejrate_si_gradual = 0.0_prec 
        ejrate_mg_gradual = 0.0_prec 
        ejrate_zn_gradual = 0.0_prec 
        ejrate_tot_gradual = 0.0_prec 
        
        countr_init = countr_old(countr-2)+1
        if (countr_init < 0) countr_init = 1 
        do jc = countr_init, countr-1
           ejrate_fe_gradual = ejrate_fe_gradual + ejrate_array(jc, 2)
           ejrate_c_gradual = ejrate_c_gradual + ejrate_array(jc, 3)
           ejrate_o_gradual = ejrate_o_gradual + ejrate_array(jc, 4)
           ejrate_n_gradual = ejrate_n_gradual + ejrate_array(jc, 5)
           ejrate_si_gradual = ejrate_si_gradual + ejrate_array(jc, 6)
           ejrate_mg_gradual = ejrate_mg_gradual + ejrate_array(jc, 7)
           ejrate_zn_gradual = ejrate_zn_gradual + ejrate_array(jc, 8)
           ejrate_tot_gradual =ejrate_tot_gradual + ejrate_array(jc, 9)
        end do

        ejrate_fe_gradual = ejrate_fe_gradual * dz * dtdz(z) / Enrich_time_lag2 
        ejrate_c_gradual = ejrate_c_gradual * dz * dtdz(z) / Enrich_time_lag2 
        ejrate_o_gradual = ejrate_o_gradual * dz * dtdz(z) / Enrich_time_lag2 
        ejrate_n_gradual = ejrate_n_gradual * dz * dtdz(z) / Enrich_time_lag2 
        ejrate_si_gradual = ejrate_si_gradual * dz * dtdz(z) / Enrich_time_lag2 
        ejrate_mg_gradual = ejrate_mg_gradual * dz * dtdz(z) / Enrich_time_lag2 
        ejrate_zn_gradual = ejrate_zn_gradual * dz * dtdz(z) / Enrich_time_lag2 
        ejrate_tot_gradual =ejrate_tot_gradual * dz * dtdz(z) / Enrich_time_lag2 

     else 

        ejrate_fe_gradual = ejrate_fe
        ejrate_c_gradual = ejrate_c
        ejrate_o_gradual = ejrate_o
        ejrate_n_gradual = ejrate_n
        ejrate_si_gradual = ejrate_si
        ejrate_mg_gradual = ejrate_mg
        ejrate_zn_gradual = ejrate_zn
        ejrate_tot_gradual =ejrate_tot
        
     endif

     xism_fe = xism_fe + dz*dtdz(z)*(acc*(xigm_fe-xism_fe) + &
          &(ejrate_fe_gradual-ejc*xism_fe))/m_ism ! dimensionless
     xism_c = xism_c + dz*dtdz(z)*(acc*(xigm_c-xism_c) + &
          &(ejrate_c_gradual-ejc*xism_c))/m_ism ! dimensionless
     xism_o = xism_o + dz*dtdz(z)*(acc*(xigm_o-xism_o) + &
          &(ejrate_o_gradual-ejc*xism_o))/m_ism ! dimensionless
     xism_n = xism_n + dz*dtdz(z)*(acc*(xigm_n-xism_n) + &
          &(ejrate_n_gradual-ejc*xism_n))/m_ism ! dimensionless
     xism_si = xism_si + dz*dtdz(z)*(acc*(xigm_si-xism_si) + &
          &(ejrate_si_gradual-ejc*xism_si))/m_ism ! dimensionless
     xism_mg = xism_mg + dz*dtdz(z)*(acc*(xigm_mg-xism_mg) + &
          &(ejrate_mg_gradual-ejc*xism_mg))/m_ism ! dimensionless
     xism_zn = xism_zn + dz*dtdz(z)*(acc*(xigm_zn-xism_zn) + &
          &(ejrate_zn_gradual-ejc*xism_zn))/m_ism ! dimensionless
     xism_tot = xism_tot + dz*dtdz(z)*(acc*(xigm_tot-xism_tot) + &
          &(ejrate_tot_gradual-ejc*xism_tot))/m_ism ! dimensionless

     m_h = 0.76*m_ism ! 10^10 M_solar Mpc^-3 
     m_fe = xism_fe*m_ism ! 10^10 M_solar Mpc^-3 
     m_c = xism_c*m_ism ! 10^10 M_solar Mpc^-3 
     m_o = xism_o*m_ism ! 10^10 M_solar Mpc^-3 
     m_totz = xism_tot*m_ism ! 10^10 M_solar Mpc^-3 

     metarr(countr) = xism_tot/0.02_prec ! dimensionless 

     if (m_fe==0.0_prec) then 
        fe_abundance = -10.0_prec 
     else 
        fe_abundance = log10(abs(m_fe/m_h)) + 2.78 
     end if
     febyharr(countr) = fe_abundance
     zfearr(countr) = z 

     if (m_c==0.0_prec) then 
        c_abundance = -10.0_prec 
     else 
        c_abundance = log10(abs(m_c/m_h)) + 2.36
     end if

     if (m_o==0.0_prec) then 
        o_abundance = -10.0_prec 
     else 
        o_abundance = log10(abs(m_o/m_h)) + 1.87
     end if

     dtrans = log10(10.0**c_abundance + 0.3*10.0**o_abundance)

     write (31, '(F4.1,5E11.3E2)') z, metarr(countr), fe_abundance, c_abundance, o_abundance, dtrans
     write (32, '(F4.1,8E11.3E2)') z, fb_struct, m_igm, m_str, m_ism, m_c, m_o, m_fe, m_totz 
     write (33, '(F4.1,3E11.3E2)') z, acc, ofl, ejc
     write (35, '(F4.1,8E11.3E2)') z, xigm_fe, xigm_c, xigm_o, xigm_tot, xism_fe, &
          &xism_c, xism_o, xism_tot  
     write (36, '(F4.1,4E11.3E2)') z, ejc, ejrate_fe, ejrate_o, ejrate_c 

     !-----------------------------

     mminc = getjmc(z)
     mminh = getjmh(z)
     write (47, '(F4.1,3E11.3E2)') z, mminc, mminh, q*mminh+(1.0_prec-q)*mminc

     ! The do-loop below evolves ``cool'' haloes i. e. haloes in H I regions. 
     do i = 1, n_halocalc

        if (ismz_halosc(i) < METALLICITY_POP3TRANS) then 
           strpop_halosc(i) = 3 
           halopop_cold(countr-1, i) = 3 
        else 
           strpop_halosc(i) = 2 
           halopop_cold(countr-1, i) = 2 
           if (NotEnrichedYet_cold(i)) then 
              Zcr_redshift = z 
              call interpolate2(tarr, zarr, Zcr_redshift, Zcr_time) ! yr 
              Zcr_timec(i) = Zcr_time 
              NotEnrichedYet_cold(i) = .false. 
           end if
        end if

        mdot = MassAccretionRate(z, m_halosc(i)*1.0e10_prec) ! 10^10 M_solar yr^-1
        m_halosc(i) = m_halosc(i) + dz*dtdz(z)*mdot ! 10^10 M_solar yr^-1

        !----------------------------

        ! Delta_c is overdensity (~168) of virialized halo.  Definition
        ! used here is from Barkana and Loeb 2001, Phys. Rep. 349,
        ! 125-238; Eqns. 22, 23.
        Omega_mz = Omega_nr * (1.0_prec+z)**3 / (Omega_nr*(1.0_prec+z)**3 + Omega_lambda &
             &+ Omega_k*(1.0_prec+z)**2) ! dimensionless 
        d = Omega_mz - 1.0_prec ! dimensionless 
        Delta_c = 18.0_prec*pi*pi + 82.0_prec*d + 39.0_prec*d*d ! dimensionless

        ! Definition of halo virial radius is from Barkana and Loeb 2001,
        ! Phys. Rep. 349, 125-238; Eqn. 24.
        HaloVirialRadius = 0.784_prec * (m_halosc(i)/(1.0e-2/smallh))**(1.0_prec/3.0_prec) *&
             &((Omega_nr/Omega_mz)*(Delta_c/(18.0_prec*pi*pi)))**(-1.0_prec/3.0_prec) *&
             &(10.0_prec/(1.0_prec+z)) / smallh ! kpc 

        ! Definitions for the DiscSpinParameter and DiscScaleLength,
        ! and the value of DiscSpinParameter, are taken from Krumholz
        ! and Dekel 2012, arXiv:1106.0301.
        DiscSpinParameter = 0.07_prec ! dimensionless 
        DiscScaleLength = 0.05_prec * (DiscSpinParameter/0.1_prec) * HaloVirialRadius ! kpc

        ! HaloCircularVelocity is the circular velocity of the halo at
        ! virial radius.  Definition used here is from Barkana and Loeb
        ! 2001, Phys. Rep. 349, 125-238; Eqn. 25.
        HaloCircularVelocity = 23.4_prec * (m_halosc(i)/(1.0e-2/smallh))**(1.0_prec/3.0_prec) *&
             &((Omega_nr/Omega_mz)*(Delta_c/(18.0_prec*pi*pi)))**(1.0_prec/6.0_prec) *&
             &((1.0_prec+z)/10.0_prec)**0.5_prec ! km/s 

        ! HaloVirialTemp is the virial temperature of the halo.
        ! Definition used here is from Barkana and Loeb 2001,
        ! Phys. Rep. 349, 125-238; Eqn. 26.
        HaloVirialTemp = 1.98e4_prec * (mu_MeanMolWt/0.6_prec) *&
             (m_halosc(i)/(1.0e-2/smallh))**(2.0_prec/3.0_prec) *&
             ((Omega_nr/Omega_mz)*(Delta_c/(18.0_prec*pi*pi)))**(1.0_prec/3.0_prec)*&
             ((1.0_prec+z)/10.0_prec) ! K 

!!$     call interpolate2(sd93_coolrate, sd93_tvir, HaloVirialTemp, GasCoolingRate) ! J m^3 s^-1
        call bi_interpolate2(sd93_temp, sd93_febyh, sd93_lambda, HaloVirialTemp, &
             &febyh_halosc(i), GasCoolingRate) ! J m^3 s^-1
        GasCoolingRate = abs(GasCoolingRate)

        p = (3.0_prec-eta)/eta
        hb = hubp(z) / yrbys ! s^-1
        fb = mgas_halosc(i)/m_halosc(i) 
        ft = 3.0_prec*mproton*kboltz*HaloVirialTemp/(2.0_prec*mu_MeanMolWt*&
             &GasCoolingRate) ! s (kg/m^3) 
        mcooldot = p*mgas_halosc(i)*(((3.0_prec-eta)*fb/(4.0_prec*pi*newtg*ft))**p)*&
             ((10.0_prec*hb)**(3.0_prec/eta)) ! 10^10_Msolar s^-1 
        mcooldot = mcooldot * yrbys ! 10^10_Msolar yr^-1 

        mcoolgas_halosc(i) = dz*dtdz(z)*mcooldot ! 10^10 M_solar 

        ! Expression for galaxy dynamical time is taken from Bouch\'e
        ! et al. 2010, Ap. J. 718, 1001-1018, Eqn. 8.
        tdyn = 2.0e7_prec * (DiscScaleLength/4.0_prec) * (200.0_prec/HaloCircularVelocity) ! yr

        ! tdyn = timedyn(z) ! yr 

        if (strpop_halosc(i) == 3) then 
           mstardot_insitu = fstar_pop3 * mcoolgas_halosc(i) / tdyn ! 10^10 M_solar yr^-1
        else 
           mstardot_insitu = fstar * mcoolgas_halosc(i) / tdyn ! 10^10 M_solar yr^-1
        end if
        sfrarr_halocalc_cold(countr-1,i) = mstardot_insitu*1.0e10_prec ! M_solar yr^-1 

        if (sfrarr_halocalc_cold(countr-2,i) .eq. 0.0_prec) then 
           FirstStar_redshift = z 
           call interpolate2(tarr, zarr, FirstStar_redshift, FirstStar_time) ! yr 
           FirstStar_timec(i) = FirstStar_time 
        end if

        return_fraction = ejfrac_nonira(z,1,i)
        ejfrac_nonira_array(countr-1, i, 1) = return_fraction 

        mstardot = mstardot_insitu - return_fraction ! 10^10 M_solar yr^-1
        burst_enhancement = MajorMergerRate(m_halosc(i)*1.0e10_prec,z) &
             &* dzdt(z) * m_halosc(i) * starburst_duty_cycle ! 10^10 M_solar yr^-1 
        mstardot_halosc(i) = mstardot + abs(burst_enhancement) ! 10^10 M_solar yr^-1
        mstar_halosc(i) = mstar_halosc(i) + dz*dtdz(z)*mstardot ! 10^10 M_solar yr^-1

        !----------------------------

        if (m_halosc(i) < mminc) then 
           fgas_in = 0.0_prec 
        else 
           fgas_in = 0.7_prec 
        end if
        mgasdot = fgas_in * (omega_b/omega_nr) * mdot ! 10^10 M_solar yr^-1
        mgasdot = mgasdot - mstardot ! 10^10 M_solar yr^-1
        fout = outfrac_nonira(m_halosc(i),z,1,i) ! 10^10 M_solar yr^-1
        mgasdot = mgasdot - fout ! 10^10 M_solar yr^-1
        mgas_halosc(i) = mgas_halosc(i) + dz*dtdz(z)*mgasdot ! 10^10 M_solar 

        !----------------------------

        if (mgas_halosc(i) == 0.0_prec) then 

           mmetal_halosc(i) = 0.0_prec
           HaloMetalAbundance = 0.0_prec
           ismz_halosc(i) = HaloMetalAbundance 
           mC_halosc(i) = 0.0_prec
           mFe_halosc(i) = 0.0_prec
           mO_halosc(i) = 0.0_prec

        else

           zeta = 0.9_prec * exp(-m_halosc(i)/3.0e2_prec) ! dimensionless 

           hy = haloyield_nonira(z,1,i) 
           haloyield_nonira_array(countr-1, i, 1) = hy
           if (Enrich_time_lag2 > 0.0_prec) then 
              hy_gradual = 0.0_prec 
              countr_init = countr_old(countr-2)+1
              if (countr_init < 0) countr_init = 1 
              do jc = countr_init, countr-1
                 hy_gradual = hy_gradual + haloyield_nonira_array(jc, i, 1)
              end do
              hy_gradual = hy_gradual * dz * dtdz(z) /  Enrich_time_lag2
           else 
              hy_gradual = hy 
           endif
           mmetaldot = fgas_in*(omega_b/omega_nr)*mdot*xigm_tot - &
                &fout*mmetal_halosc(i)/mgas_halosc(i)  &
                &+ hy_gradual*(1.0_prec-zeta) ! 10^10 M_solar yr^-1

           mmetal_halosc(i) = mmetal_halosc(i) + dz*dtdz(z)*mmetaldot ! 10^10 M_solar 
           HaloMetalAbundance = (mmetal_halosc(i)/mgas_halosc(i))/0.02_prec 
           ismz_halosc(i) = HaloMetalAbundance 
           
           hy = haloyield_species_nonira(z,7,1,i)
           haloyield_species_nonira_array(countr-1, i, 1, 1) = hy 
           if (Enrich_time_lag2 > 0.0_prec) then 
              hy_gradual = 0.0_prec 
              countr_init = countr_old(countr-2)+1
              if (countr_init < 0) countr_init = 1 
              do jc = countr_init, countr-1
                 hy_gradual = hy_gradual + haloyield_species_nonira_array(jc, i, 1, 1)
              end do
              hy_gradual = hy_gradual * dz * dtdz(z) /  Enrich_time_lag2
           else 
              hy_gradual = hy 
           endif
           mCdot = fgas_in*(omega_b/omega_nr)*mdot*xigm_c - &
                &fout*mC_halosc(i)/mgas_halosc(i)  &
                &+ hy_gradual*(1.0_prec-zeta) ! 10^10 M_solar yr^-1

           mC_halosc(i) = mC_halosc(i) + dz*dtdz(z)*mCdot ! 10^10 M_solar 

           hy = haloyield_species_nonira(z,14,1,i)
           haloyield_species_nonira_array(countr-1, i, 2, 1) = hy 
           if (Enrich_time_lag2 > 0.0_prec) then 
              hy_gradual = 0.0_prec 
              countr_init = countr_old(countr-2)+1
              if (countr_init < 0) countr_init = 1 
              do jc = countr_init, countr-1
                 hy_gradual = hy_gradual + haloyield_species_nonira_array(jc, i, 2, 1)
              end do
              hy_gradual = hy_gradual * dz * dtdz(z) /  Enrich_time_lag2
           else 
              hy_gradual = hy 
           endif
           mFedot = fgas_in*(omega_b/omega_nr)*mdot*xigm_fe - &
                &fout*mFe_halosc(i)/mgas_halosc(i)  &
                &+ hy_gradual*(1.0_prec-zeta) ! 10^10 M_solar yr^-1

           mFe_halosc(i) = mFe_halosc(i) + dz*dtdz(z)*mFedot ! 10^10 M_solar  

           hy = haloyield_species_nonira(z,9,1,i)
           haloyield_species_nonira_array(countr-1, i, 3, 1) = hy 
           if (Enrich_time_lag2 > 0.0_prec) then 
              hy_gradual = 0.0_prec 
              countr_init = countr_old(countr-2)+1
              if (countr_init < 0) countr_init = 1 
              do jc = countr_init, countr-1
                 hy_gradual = hy_gradual + haloyield_species_nonira_array(jc, i, 3, 1)
              end do
              hy_gradual = hy_gradual * dz * dtdz(z) /  Enrich_time_lag2
           else 
              hy_gradual = hy 
           endif
           mOdot = fgas_in*(omega_b/omega_nr)*mdot*xigm_o - &
                &fout*mO_halosc(i)/mgas_halosc(i)  &
                &+ hy_gradual*(1.0_prec-zeta) ! 10^10 M_solar yr^-1

           mO_halosc(i) = mO_halosc(i) + dz*dtdz(z)*mOdot ! 10^10 M_solar 

           hy = haloyield_species_nonira(z,8,1,i)
           haloyield_species_nonira_array(countr-1, i, 4, 1) = hy 
           if (Enrich_time_lag2 > 0.0_prec) then 
              hy_gradual = 0.0_prec 
              countr_init = countr_old(countr-2)+1
              if (countr_init < 0) countr_init = 1 
              do jc = countr_init, countr-1
                 hy_gradual = hy_gradual + haloyield_species_nonira_array(jc, i, 4, 1)
              end do
              hy_gradual = hy_gradual * dz * dtdz(z) /  Enrich_time_lag2
           else 
              hy_gradual = hy 
           endif
           mNdot = fgas_in*(omega_b/omega_nr)*mdot*xigm_n - &
                &fout*mN_halosc(i)/mgas_halosc(i)  &
                &+ hy_gradual*(1.0_prec-zeta) ! 10^10 M_solar yr^-1

           mN_halosc(i) = mN_halosc(i) + dz*dtdz(z)*mNdot ! 10^10 M_solar 

           hy = haloyield_species_nonira(z,11,1,i)
           haloyield_species_nonira_array(countr-1, i, 5, 1) = hy 
           if (Enrich_time_lag2 > 0.0_prec) then 
              hy_gradual = 0.0_prec 
              countr_init = countr_old(countr-2)+1
              if (countr_init < 0) countr_init = 1 
              do jc = countr_init, countr-1
                 hy_gradual = hy_gradual + haloyield_species_nonira_array(jc, i, 5, 1)
              end do
              hy_gradual = hy_gradual * dz * dtdz(z) /  Enrich_time_lag2
           else 
              hy_gradual = hy 
           endif
           mSidot = fgas_in*(omega_b/omega_nr)*mdot*xigm_si - &
                &fout*mSi_halosc(i)/mgas_halosc(i)  &
                &+ hy_gradual*(1.0_prec-zeta) ! 10^10 M_solar yr^-1

           mSi_halosc(i) = mSi_halosc(i) + dz*dtdz(z)*mSidot ! 10^10 M_solar 

           hy = haloyield_species_nonira(z,15,1,i)
           haloyield_species_nonira_array(countr-1, i, 6, 1) = hy

           if (Enrich_time_lag2 > 0.0_prec) then 
              hy_gradual = 0.0_prec 
              countr_init = countr_old(countr-2)+1
              if (countr_init < 0) countr_init = 1 
              do jc = countr_init, countr-1
                 hy_gradual = hy_gradual + haloyield_species_nonira_array(jc, i, 6, 1)
              end do
              hy_gradual = hy_gradual * dz * dtdz(z) /  Enrich_time_lag2
           else 
              hy_gradual = hy 
           endif
           mZndot = fgas_in*(omega_b/omega_nr)*mdot*xigm_zn - &
                &fout*mZn_halosc(i)/mgas_halosc(i)  &
                &+ hy_gradual*(1.0_prec-zeta) ! 10^10 M_solar yr^-1
           mZn_halosc(i) = mZn_halosc(i) + dz*dtdz(z)*mZndot ! 10^10 M_solar 

           hy = haloyield_species_nonira(z,10,1,i)
           haloyield_species_nonira_array(countr-1, i, 7, 1) = hy
           if (Enrich_time_lag2 > 0.0_prec) then 
              hy_gradual = 0.0_prec 
              countr_init = countr_old(countr-2)+1
              if (countr_init < 0) countr_init = 1 
              do jc = countr_init, countr-1
                 hy_gradual = hy_gradual + haloyield_species_nonira_array(jc, i, 7, 1)
              end do
              hy_gradual = hy_gradual * dz * dtdz(z) /  Enrich_time_lag2
           else 
              hy_gradual = hy 
           endif
           mMgdot = fgas_in*(omega_b/omega_nr)*mdot*xigm_mg - &
                &fout*mMg_halosc(i)/mgas_halosc(i)  &
                &+ hy_gradual*(1.0_prec-zeta) ! 10^10 M_solar yr^-1

           mMg_halosc(i) = mMg_halosc(i) + dz*dtdz(z)*mMgdot ! 10^10 M_solar 

        end if

        mH_halos = 0.71*mgas_halosc(i) 

        if (mFe_halosc(i)==0.0_prec) then 
           febyh_halosc(i) = -10.0_prec 
        else 
           febyh_halosc(i) = log10(abs(mFe_halosc(i)/mH_halos)) + 2.78 
        end if

        if (mC_halosc(i)==0.0_prec) then 
           cbyh_halosc(i) = -10.0_prec 
        else 
           cbyh_halosc(i) = log10(abs(mC_halosc(i)/mFe_halosc(i))) - 0.41
        end if

        if (mO_halosc(i)==0.0_prec) then 
           obyh_halosc(i) = -10.0_prec 
        else 
           obyh_halosc(i) = log10(abs(mO_halosc(i)/mH_halos)) + 1.87
        end if
     end do

     ! The do-loop below evolves ``hot'' haloes i. e. haloes in H II regions. 

     do i = 1, n_halocalc 

        if (ismz_halosh(i) < METALLICITY_POP3TRANS) then 
           strpop_halosh(i) = 3 
           halopop_hot(countr-1, i) = 3 
        else 
           strpop_halosh(i) = 2 
           halopop_hot(countr-1, i) = 2
           if (NotEnrichedYet_hot(i)) then 
              Zcr_redshift = z 
              call interpolate2(tarr, zarr, Zcr_redshift, Zcr_time) ! yr 
              Zcr_timeh(i) = Zcr_time 
              NotEnrichedYet_hot(i) = .false. 
           end if
        end if

        mdot = MassAccretionRate(z, m_halosh(i)*1.0e10_prec) ! 10^10 M_solar yr^-1
        m_halosh(i) = m_halosh(i) + dz*dtdz(z)*mdot ! 10^10 M_solar yr^-1

        !----------------------------

        ! Delta_c is overdensity (~168) of virialized halo.  Definition
        ! used here is from Barkana and Loeb 2001, Phys. Rep. 349,
        ! 125-238; Eqns. 22, 23.
        Omega_mz = Omega_nr * (1.0_prec+z)**3 / (Omega_nr*(1.0_prec+z)**3 + Omega_lambda &
             &+ Omega_k*(1.0_prec+z)**2) ! dimensionless 
        d = Omega_mz - 1.0_prec ! dimensionless 
        Delta_c = 18.0_prec*pi*pi + 82.0_prec*d + 39.0_prec*d*d ! dimensionless

        ! Definition of halo virial radius is from Barkana and Loeb 2001,
        ! Phys. Rep. 349, 125-238; Eqn. 24.
        HaloVirialRadius = 0.784_prec * (m_halosh(i)/(1.0e-2/smallh))**(1.0_prec/3.0_prec) *&
             &((Omega_nr/Omega_mz)*(Delta_c/(18.0_prec*pi*pi)))**(-1.0_prec/3.0_prec) *&
             &(10.0_prec/(1.0_prec+z)) / smallh ! kpc 

        ! Definitions for the DiscSpinParameter and DiscScaleLength,
        ! and the value of DiscSpinParameter, are taken from Krumholz
        ! and Dekel 2012, arXiv:1106.0301.
        DiscSpinParameter = 0.07_prec ! dimensionless 
        DiscScaleLength = 0.05_prec * (DiscSpinParameter/0.1_prec) * HaloVirialRadius ! kpc

        ! HaloCircularVelocity is the circular velocity of the halo at
        ! virial radius.  Definition used here is from Barkana and Loeb
        ! 2001, Phys. Rep. 349, 125-238; Eqn. 25.
        HaloCircularVelocity = 23.4_prec * (m_halosh(i)/(1.0e-2/smallh))**(1.0_prec/3.0_prec) *&
             &((Omega_nr/Omega_mz)*(Delta_c/(18.0_prec*pi*pi)))**(1.0_prec/6.0_prec) *&
             &((1.0_prec+z)/10.0_prec)**0.5_prec ! km/s 

        ! HaloVirialTemp is the virial temperature of the halo.
        ! Definition used here is from Barkana and Loeb 2001,
        ! Phys. Rep. 349, 125-238; Eqn. 26.
        HaloVirialTemp = 1.98e4_prec * (mu_MeanMolWt/0.6_prec) *&
             (m_halosh(i)/(1.0e-2/smallh))**(2.0_prec/3.0_prec) *&
             ((Omega_nr/Omega_mz)*(Delta_c/(18.0_prec*pi*pi)))**(1.0_prec/3.0_prec)*&
             ((1.0_prec+z)/10.0_prec) ! K 

!!$        call interpolate2(sd93_coolrate, sd93_tvir, HaloVirialTemp, GasCoolingRate) ! J m^3 s^-1
        call bi_interpolate2(sd93_temp, sd93_febyh, sd93_lambda, HaloVirialTemp, &
             &febyh_halosh(i), GasCoolingRate) ! J m^3 s^-1
        GasCoolingRate = abs(GasCoolingRate)

        p = (3.0_prec-eta)/eta
        hb = hubp(z) / yrbys ! s^-1
        fb = mgas_halosh(i)/m_halosh(i) 
        ft = 3.0_prec*mproton*kboltz*HaloVirialTemp/(2.0_prec*mu_MeanMolWt*GasCoolingRate) ! s (kg/m^3) 
        mcooldot = p*mgas_halosh(i)*(((3.0_prec-eta)*fb/(4.0_prec*pi*newtg*ft))**p)*&
             ((10.0_prec*hb)**(3.0_prec/eta)) ! 10^10_Msolar s^-1 
        mcooldot = mcooldot * yrbys ! 10^10_Msolar yr^-1 

        mcoolgas_halosh(i) = dz*dtdz(z)*mcooldot ! 10^10 M_solar 

        ! Expression for galaxy dynamical time is taken from Bouch\'e
        ! et al. 2010, Ap. J. 718, 1001-1018, Eqn. 8.
        tdyn = 2.0e7_prec * (DiscScaleLength/4.0_prec) * (200.0_prec/HaloCircularVelocity) ! yr

        ! tdyn = timedyn(z) ! yr 

        if (strpop_halosh(i) == 3) then 
           mstardot_insitu = fstar_pop3 * mcoolgas_halosh(i) / tdyn ! 10^10 M_solar yr^-1
        else 
           mstardot_insitu = fstar * mcoolgas_halosh(i) / tdyn ! 10^10 M_solar yr^-1
        end if
        sfrarr_halocalc_hot(countr-1,i) = mstardot_insitu*1.0e10_prec ! M_solar yr^-1 

        if (sfrarr_halocalc_hot(countr-2,i) .eq. 0.0_prec) then 
           FirstStar_redshift = z 
           call interpolate2(tarr, zarr, FirstStar_redshift, FirstStar_time) ! yr 
           FirstStar_timeh(i) = FirstStar_time 
        end if

        return_fraction = ejfrac_nonira(z, 2, i) ! 10^10 M_solar yr^-1 
        ejfrac_nonira_array(countr-1, i, 2) = return_fraction 

        mstardot = mstardot_insitu - return_fraction ! 10^10 M_solar yr^-1
        burst_enhancement = MajorMergerRate(m_halosc(i)*1.0e10_prec,z) &
             &* dzdt(z) * m_halosh(i) * starburst_duty_cycle ! 10^10 M_solar yr^-1 
        mstardot_halosh(i) = mstardot + abs(burst_enhancement) ! 10^10 M_solar yr^-1
        mstar_halosh(i) = mstar_halosh(i) + dz*dtdz(z)*mstardot ! 10^10 M_solar yr^-1

        !----------------------------

        if (m_halosh(i) < mminh) then 
           fgas_in = 0.0_prec 
        else 
           fgas_in = 0.7_prec 
        end if
        mgasdot = fgas_in * (omega_b/omega_nr) * mdot ! 10^10 M_solar yr^-1
        mgasdot = mgasdot - mstardot ! 10^10 M_solar yr^-1
        fout = outfrac_nonira(m_halosh(i),z,2,i) ! 10^10 M_solar yr^-1
        mgasdot = mgasdot - fout ! 10^10 M_solar yr^-1
        mgas_halosh(i) = mgas_halosh(i) + dz*dtdz(z)*mgasdot ! 10^10 M_solar 

        !----------------------------

        if (mgas_halosh(i) == 0.0_prec) then 

           mmetal_halosh(i) = 0.0_prec
           HaloMetalAbundance = 0.0_prec
           ismz_halosh(i) = HaloMetalAbundance 
           mC_halosh(i) = 0.0_prec
           mFe_halosh(i) = 0.0_prec
           mO_halosh(i) = 0.0_prec

        else

           zeta = 0.9_prec * exp(-m_halosh(i)/3.0e2_prec) ! dimensionless 

           hy = haloyield_nonira(z,2,i)
           haloyield_nonira_array(countr-1, i, 2) = hy 
                      
           if (Enrich_time_lag2 > 0.0_prec) then 
              hy_gradual = 0.0_prec 
              countr_init = countr_old(countr-2)+1
              if (countr_init < 0) countr_init = 1 
              do jc = countr_init, countr-1
                 hy_gradual = hy_gradual + haloyield_nonira_array(jc, i, 2)
              end do
              hy_gradual = hy_gradual * dz * dtdz(z) /  Enrich_time_lag2
           else 
              hy_gradual = hy 
           endif
           mmetaldot = fgas_in*(omega_b/omega_nr)*mdot*xigm_tot - &
                &fout*mmetal_halosh(i)/mgas_halosh(i)  &
                &+ hy_gradual*(1.0_prec-zeta)
           mmetal_halosh(i) = mmetal_halosh(i) + dz*dtdz(z)*mmetaldot   
           HaloMetalAbundance = (mmetal_halosh(i)/mgas_halosh(i))/0.02_prec 
           ismz_halosh(i) = HaloMetalAbundance 

           hy = haloyield_species_nonira(z,7,2,i)
           haloyield_species_nonira_array(countr-1, i, 1, 2) = hy 
                      
           if (Enrich_time_lag2 > 0.0_prec) then 
              hy_gradual = 0.0_prec 
              countr_init = countr_old(countr-2)+1
              if (countr_init < 0) countr_init = 1 
              do jc = countr_init, countr-1
                 hy_gradual = hy_gradual + haloyield_species_nonira_array(jc, i, 1, 2)
              end do
              hy_gradual = hy_gradual * dz * dtdz(z) /  Enrich_time_lag2
           else 
              hy_gradual = hy 
           endif
           mCdot = fgas_in*(omega_b/omega_nr)*mdot*xigm_c - &
                &fout*mC_halosh(i)/mgas_halosh(i)  &
                &+ hy_gradual*(1.0_prec-zeta) ! 10^10 M_solar yr^-1
           mC_halosh(i) = mC_halosh(i) + dz*dtdz(z)*mCdot ! 10^10 M_solar 

           hy = haloyield_species_nonira(z,14,2,i)
           haloyield_species_nonira_array(countr-1, i, 2, 2) = hy 
           if (Enrich_time_lag2 > 0.0_prec) then 
              hy_gradual = 0.0_prec 
              countr_init = countr_old(countr-2)+1
              if (countr_init < 0) countr_init = 1 
              do jc = countr_init, countr-1
                 hy_gradual = hy_gradual + haloyield_species_nonira_array(jc, i, 2, 2)
              end do
              hy_gradual = hy_gradual * dz * dtdz(z) /  Enrich_time_lag2
           else 
              hy_gradual = hy 
           endif
           mFedot = fgas_in*(omega_b/omega_nr)*mdot*xigm_fe - &
                &fout*mFe_halosh(i)/mgas_halosh(i)  &
                &+ hy_gradual*(1.0_prec-zeta) ! 10^10 M_solar yr^-1

           mFe_halosh(i) = mFe_halosh(i) + dz*dtdz(z)*mFedot ! 10^10 M_solar  

           hy = haloyield_species_nonira(z,9,2,i)
           haloyield_species_nonira_array(countr-1, i, 3, 2) = hy 
           if (Enrich_time_lag2 > 0.0_prec) then 
              hy_gradual = 0.0_prec 
              countr_init = countr_old(countr-2)+1
              if (countr_init < 0) countr_init = 1 
              do jc = countr_init, countr-1
                 hy_gradual = hy_gradual + haloyield_species_nonira_array(jc, i, 3, 2)
              end do
              hy_gradual = hy_gradual * dz * dtdz(z) /  Enrich_time_lag2
           else 
              hy_gradual = hy 
           endif
           mOdot = fgas_in*(omega_b/omega_nr)*mdot*xigm_o - &
                &fout*mO_halosh(i)/mgas_halosh(i)  &
                &+ hy_gradual*(1.0_prec-zeta) ! 10^10 M_solar yr^-1

           mO_halosh(i) = mO_halosh(i) + dz*dtdz(z)*mOdot ! 10^10 M_solar 

           hy = haloyield_species_nonira(z,8,2,i)
           haloyield_species_nonira_array(countr-1, i, 4, 2) = hy 
           if (Enrich_time_lag2 > 0.0_prec) then 
              hy_gradual = 0.0_prec 
              countr_init = countr_old(countr-2)+1
              if (countr_init < 0) countr_init = 1 
              do jc = countr_init, countr-1
                 hy_gradual = hy_gradual + haloyield_species_nonira_array(jc, i, 4, 2)
              end do
              hy_gradual = hy_gradual * dz * dtdz(z) /  Enrich_time_lag2
           else 
              hy_gradual = hy 
           endif
           mNdot = fgas_in*(omega_b/omega_nr)*mdot*xigm_n - &
                &fout*mN_halosh(i)/mgas_halosh(i)  &
                &+ hy_gradual*(1.0_prec-zeta) ! 10^10 M_solar yr^-1

           mN_halosh(i) = mN_halosh(i) + dz*dtdz(z)*mNdot ! 10^10 M_solar 

           hy = haloyield_species_nonira(z,11,2,i)
           haloyield_species_nonira_array(countr-1, i, 5, 2) = hy 
           if (Enrich_time_lag2 > 0.0_prec) then 
              hy_gradual = 0.0_prec 
              countr_init = countr_old(countr-2)+1
              if (countr_init < 0) countr_init = 1 
              do jc = countr_init, countr-1
                 hy_gradual = hy_gradual + haloyield_species_nonira_array(jc, i, 5, 2)
              end do
              hy_gradual = hy_gradual * dz * dtdz(z) /  Enrich_time_lag2
           else 
              hy_gradual = hy 
           endif
           mSidot = fgas_in*(omega_b/omega_nr)*mdot*xigm_si - &
                &fout*mSi_halosh(i)/mgas_halosh(i)  &
                &+ hy_gradual*(1.0_prec-zeta) ! 10^10 M_solar yr^-1

           mSi_halosh(i) = mSi_halosh(i) + dz*dtdz(z)*mSidot ! 10^10 M_solar 

           hy = haloyield_species_nonira(z,15,2,i)
           haloyield_species_nonira_array(countr-1, i, 6, 2) = hy
           if (Enrich_time_lag2 > 0.0_prec) then 
              hy_gradual = 0.0_prec 
              countr_init = countr_old(countr-2)+1
              if (countr_init < 0) countr_init = 1 
              do jc = countr_init, countr-1
                 hy_gradual = hy_gradual + haloyield_species_nonira_array(jc, i, 6, 2)
              end do
              hy_gradual = hy_gradual * dz * dtdz(z) /  Enrich_time_lag2
           else 
              hy_gradual = hy 
           endif
           mZndot = fgas_in*(omega_b/omega_nr)*mdot*xigm_zn - &
                &fout*mZn_halosh(i)/mgas_halosh(i)  &
                &+ hy_gradual*(1.0_prec-zeta) ! 10^10 M_solar yr^-1

           mZn_halosh(i) = mZn_halosh(i) + dz*dtdz(z)*mZndot ! 10^10 M_solar 

           hy = haloyield_species_nonira(z,10,2,i)
           haloyield_species_nonira_array(countr-1, i, 7, 2) = hy
           if (Enrich_time_lag2 > 0.0_prec) then 
              hy_gradual = 0.0_prec 
              countr_init = countr_old(countr-2)+1
              if (countr_init < 0) countr_init = 1 
              do jc = countr_init, countr-1
                 hy_gradual = hy_gradual + haloyield_species_nonira_array(jc, i, 7, 2)
              end do
              hy_gradual = hy_gradual * dz * dtdz(z) /  Enrich_time_lag2
           else 
              hy_gradual = hy 
           endif
           mMgdot = fgas_in*(omega_b/omega_nr)*mdot*xigm_mg - &
                &fout*mMg_halosh(i)/mgas_halosh(i)  &
                &+ hy_gradual*(1.0_prec-zeta) ! 10^10 M_solar yr^-1

           mMg_halosh(i) = mMg_halosh(i) + dz*dtdz(z)*mMgdot ! 10^10 M_solar 

        end if

        mH_halos = 0.71*mgas_halosh(i) 

        if (mFe_halosh(i)==0.0_prec) then 
           febyh_halosh(i) = -10.0_prec 
        else 
           febyh_halosh(i) = log10(abs(mFe_halosh(i)/mH_halos)) + 2.78 
        end if

        if (mC_halosh(i)==0.0_prec) then 
           cbyh_halosh(i) = -10.0_prec 
        else 
           cbyh_halosh(i) = log10(abs(mC_halosh(i)/mFe_halosh(i))) - 0.41
        end if

        if (mO_halosh(i)==0.0_prec) then 
           obyh_halosh(i) = -10.0_prec 
        else 
           obyh_halosh(i) = log10(abs(mO_halosh(i)/mH_halos)) + 1.87
        end if
     end do

     write (39,'(F4.1,270E11.3E2)') z, m_halosh 
     write (40,'(F4.1,270E11.3E2)') z, mstar_halosh 
     write (41,'(F4.1,270E11.3E2)') z, mgas_halosh 
     write (42,'(F4.1,270E11.3E2)') z, mmetal_halosh 
     write (43,'(F4.1,270E11.3E2)') z, mC_halosh 
     write (44,'(F4.1,270E11.3E2)') z, mFe_halosh 
     write (45,'(F4.1,270E11.3E2)') z, mO_halosh 
     write (59,'(F4.1,270E11.3E2)') z, mN_halosh 
     write (60,'(F4.1,270E11.3E2)') z, mSi_halosh 
     write (61,'(F4.1,270E11.3E2)') z, mZn_halosh 
     write (67,'(F4.1,270E11.3E2)') z, mMg_halosh 
     write (46,'(F4.1,270E11.3E2)') z, ismz_halosh 
     write (48,'(F4.1,270I3)') z, strpop_halosh
     write (74,'(F4.1,270E11.3E2)') z, mcoolgas_halosh 

     write (76,'(F4.1,270E11.3E2)') z, sfrcontrib_halosh 
     write (77,'(F4.1,270E11.3E2)') z, mstardot_halosh 
     write (78,'(F4.1,270E11.3E2)') z, q*sfrcontrib_halosh + (1.0_prec-q)*sfrcontrib_halosc 
     write (79,'(F4.1,270E11.3E2)') z, q*mstardot_halosh + (1.0_prec-q)*sfrcontrib_halosc 
     write (80,'(F4.1,270E11.3E2)') z, sfrcontrib_halosc 

     write (49,'(F4.1,270E11.3E2)') z, aux_halosc

     write (50,'(F4.1,270E11.3E2)') z, q*m_halosh + (1.0_prec-q)*m_halosc 
     write (51,'(F4.1,270E11.3E2)') z, q*mstar_halosh + (1.0_prec-q)*mstar_halosc 
     write (52,'(F4.1,270E11.3E2)') z, q*mgas_halosh + (1.0_prec-q)*mgas_halosc 
     write (53,'(F4.1,270E11.3E2)') z, q*mmetal_halosh + (1.0_prec-q)*mmetal_halosc 
     write (54,'(F4.1,270E11.3E2)') z, q*mC_halosh + (1.0_prec-q)*mC_halosc 
     write (55,'(F4.1,270E11.3E2)') z, q*mFe_halosh + (1.0_prec-q)*mFe_halosc 
     write (56,'(F4.1,270E11.3E2)') z, q*mO_halosh + (1.0_prec-q)*mO_halosc 
     write (62,'(F4.1,270E11.3E2)') z, q*mN_halosh + (1.0_prec-q)*mN_halosc 
     write (63,'(F4.1,270E11.3E2)') z, q*mSi_halosh + (1.0_prec-q)*mSi_halosc 
     write (64,'(F4.1,270E11.3E2)') z, q*mZn_halosh + (1.0_prec-q)*mZn_halosc 
     write (68,'(F4.1,270E11.3E2)') z, q*mMg_halosh + (1.0_prec-q)*mMg_halosc 
     write (57,'(F4.1,270E11.3E2)') z, q*ismz_halosh + (1.0_prec-q)*ismz_halosc 
     write (75,'(F4.1,270E11.3E2)') z, q*mcoolgas_halosh + (1.0_prec-q)*mcoolgas_halosc 

     write (65,'(F4.1,270E11.3E2)') z, nofmc
     write (66,'(F4.1,270E11.3E2)') z, nofmh

     write (73, '(F4.1,270E11.3E2)') z, halol1500c

     write (81, '(F4.1,270E11.3E2)') z, Zcr_timeh
     write (82, '(F4.1,270E11.3E2)') z, Zcr_timec
     write (83, '(F4.1,270E11.3E2)') z, FirstStar_timeh
     write (84, '(F4.1,270E11.3E2)') z, FirstStar_timec

  END DO

  close(31)
  close(32)
  close(33)
  close(34)
  close(35)
  close(36)
  close(37)
  close(38)
  close(39)
  close(40)
  close(41)
  close(42)
  close(43)
  close(44)
  close(45)
  close(46)
  close(47)
  close(48)
  close(49)
  close(50)
  close(51)
  close(52)
  close(53)
  close(54)
  close(55)
  close(56)
  close(57)
  close(58)
  close(59)
  close(60)
  close(61)
  close(62)
  close(63)
  close(64)
  close(65)
  close(66)
  close(67)
  close(68)
  close(73)
  close(74)  
  close(75)
  close(76)
  close(77)
  close(78)
  close(79)
  close(80)

  close(81)
  close(82)
  close(83)
  close(84)

  close(69)
  close(70)
  close(71)
  close(72)

  TAU = TAU-FINAL_REDSHIFT*TAU1 
  WRITE(0,*) 'TAU=', TAU 

  !-------------------------

CONTAINS 

  FUNCTION DTEMPCDZ(REDSHIFT)

    ! Redshift derivative of temperature of neutral region.
    IMPLICIT NONE
    REAL(KIND = PREC), INTENT(IN) :: REDSHIFT
    REAL(KIND = PREC) :: DTEMPCDZ

    DTEMPCDZ = -2.0_PREC*HUBP(REDSHIFT)*TEMPC*DTDZ(REDSHIFT) ! K 

  END FUNCTION DTEMPCDZ

  FUNCTION TAUINT(FM, F2, X2, Z) 

    IMPLICIT NONE
    REAL(PREC), INTENT(IN) :: FM, F2, X2, Z
    REAL(PREC) :: TAUINT
    REAL(PREC) :: ECHARGE, MELECTRON, THOMSCROSS,&
         &ZDOT, HUB0_SMALLUNITS, NHCM, NE

    ECHARGE = 4.803204E-10_PREC ! esu
    MELECTRON = 9.109382E-28_PREC ! grams

    HUB0_SMALLUNITS = 3.24175E-18_PREC * SMALLH ! s^-1 
    THOMSCROSS = (8.0_PREC*PI/3.0_PREC)*(ECHARGE*ECHARGE&
         &/(MELECTRON*SPEED_OF_LIGHT**2))**2 ! cm^2  
    ZDOT = -HUB0_SMALLUNITS*(1.0_PREC+Z)*SQRT(OMEGA_NR*(1.0_PREC+Z)**3&
         & +OMEGA_LAMBDA) ! s^-1 
    NHCM = NH*CMBYMPCCB ! cm^-3 

    NE = FM*NHCM*X2*F2 + FM*(1.0_PREC-F2)*NHCM*X2INIT + &
         &(1.0_PREC-FM)*NHCM*X2INIT
    TAUINT = NE*SPEED_OF_LIGHT*THOMSCROSS&
         &*(1.0_PREC+Z)**3/ZDOT ! dimensionless 
  END FUNCTION TAUINT

  function nstar(mu, mstar)

    implicit none 
    real(kind=prec), intent(in) :: mu, mstar
    real(kind=prec) :: nstar

    real(kind=prec) :: norm, ml

    ml = 0.1_prec ! M_solar
    norm = mstar*(1.0_prec-imf_slope)/(msup**(1.0_prec-imf_slope)-minf**(1.0_prec-imf_slope))

    nstar = -(norm/imf_slope)*(mu**(-imf_slope)-ml**(-imf_slope)) 

  end function nstar

  function MassAccretionRate(z, HaloCurrentMass) 

    real(kind=prec), intent(in) :: z, HaloCurrentMass 
    real(kind=prec) :: MassAccretionRate 

    ! Fakhouri et al. 2010, MNRAS 406, 2267-2278; Eqn. 2.
    MassAccretionRate = 46.1_prec*((HaloCurrentMass/1.0e12_prec)**1.1_prec)*&
         &(1.0_prec + 1.11_prec*z)*&
         &sqrt(omega_nr*(1.0_prec+z)**3 + omega_lambda) ! M_solar yr^-1

    MassAccretionRate = MassAccretionRate * 1.0e-10_prec ! 10^10 M_solar yr^-1

  end function MassAccretionRate

  FUNCTION timedyn(rs)

    REAL(KIND = PREC), INTENT(IN) :: rs
    REAL(KIND = PREC) :: timedyn

    REAL(KIND = PREC) :: HDENS, BKGRHO, RHOCS

    RHOCS = 1.879E-29_PREC * SMALLH ** 2 ! g/cm^3
    BKGRHO = RHOCS * OMEGA_NR * (1.0_PREC + rs) ** 3 ! g/cm^3 
    HDENS = DELTA_VIRIAL * BKGRHO ! g/cm^3 
    timedyn = SQRT(3.0_PREC * PI / (16.0_PREC * NEWTG * &
         &HDENS * 1.0E3_PREC)) / YRBYS ! yr
    ! 1.0E3 converts g to kg and cm to m.  

  END FUNCTION timedyn

  function frac_ej(t) 

    implicit none 
    real(kind = prec), intent(in) :: t 
    real(kind = prec) :: frac_ej 
    
    real(kind = prec) :: slope 

    if (t > Enrich_time_lag2) then 
       frac_ej = 0.0_prec 
    else if (t < Enrich_time_lag2) then 
       frac_ej = 0.0_prec 
    else  
       ! This slope is determined from the condition that all mass should
       ! be ejected between time 0 and Enrich_time_lag2, i.e., the
       ! integral of frac_ej(t)dt should be unity over this interval.
       slope = 1.0_prec / Enrich_time_lag2 ! yr^-1
       frac_ej = slope ! yr^-1 
    end if
       
  end function frac_ej
  
END PROGRAM REION

