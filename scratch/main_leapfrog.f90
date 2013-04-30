
PROGRAM REION 

  ! File: main.f90
  ! Cre: 2010-04-07
  ! Mod: $Date: 2012/03/09 08:17:40 $; $Revision: 1.1 $ 

  USE CONSTANTS 
  USE STORAGE 
  USE INTERFACES, ONLY : CLUMPFAC, DTDZ, GAMMA_PH, GAMMA_PI, HUBP, &
       &IGMFPREO, IGMVFRAC, JEANS_LENGTH, PSIG, RTBIS, SIGMAH , &
       &SIGMA_BARYON, SOLDLT, SOLFM, LUMFN, gallum, COUNTER, &
       &SOURCEOV_POP2, SOURCEOV_POP3, GPH_KERNEL_POP2, outflow, &
       &GPI_KERNEL_POP2, GPI_KERNEL_POP3, GPH_KERNEL_POP3, &
       &ACCRATE, ejrate, sfr_rollinde_pop3, sfr_rollinde_pop2  
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
       &GPH_POP2, GPH_POP3, ST_MASS, ST_AGE, YIELD_RECORD(10), ofl, &
       &ejr, LOCAL_HUBBLE_0, metallicity, metmass, m_c, m_o, m_fe, &
       &m_h, fe_abundance, c_abundance, o_abundance, source_dummy, &
       &r1_igm, r2_igm, r1_ism, r2_ism, r1_str, r2_str, r1igm_fe, &
       &r2igm_fe, r1igm_c, r2igm_c, r1igm_o, r2igm_o, r1ism_fe, r1ism_o,&
       &r1ism_c, r2ism_fe, r2ism_c, r2ism_o

  real(kind=prec) :: m_igm, m_ism, m_str, xigm_fe, xigm_c, xigm_o, &
       & xism_fe, xism_c, xism_o, dm_ism, dm_igm, dm_str

  REAL(KIND=PREC) :: MAGFOO, ZFOO, NFOO1, NFOO2, MSTAR, MIGM, MISM 

  INTEGER :: I, J, IER, INF, LAST, NEVAL, NUMBER_OF_LINES, NCALC, comnum
  CHARACTER(100) :: ZLUMARG, RGNOVDARG, FESCARG, RGNSIZEARG 

  LOGICAL :: PREOVERLAP 

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

  ! ----------------------------------------

  ! Set initial conditions. 
  INCLUDE 'readin.inc' 

  ! NGAMMA is calculated from Starburst99 model
  ! `reion-generic' by popsyn/ngtot.f90.
  NGAMMA = 2.82E70_PREC ! (10^10 M_solar)^-1 
  NGAMMA_POP2 = NGAMMA 
  NGAMMA_POP3 = 2.82E71_PREC ! (10^10 M_solar)^-1 

  Z = INITIAL_REDSHIFT 
  Q = 1.0E-15_PREC 

  ! Allocate and initialize arrays that will store values of Jeans
  ! mass and filling factor for redshifts at which calculation is 
  ! done.  Length of these arrays will have to be equal to the 
  ! number of such redshift values. 
  NCALC = (FINAL_REDSHIFT-INITIAL_REDSHIFT)/DZ+1 
  ALLOCATE(JMHARR(NCALC)); ALLOCATE(FFARR(NCALC)); ALLOCATE(SFRARR(NCALC))
  allocate(metarr(ncalc))
  JMHARR = 0.0_PREC; FFARR = 0.0_PREC; SFRARR = 0.0_PREC  
  metarr = 0.0_prec
  COUNTR = 1 
  FFARR(COUNTR) = Q 

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

  JNSLN = JEANS_LENGTH(TEMPH,Z) ! Mpc 
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
  F_M = IGMFPREO()
  FV = IGMVFRAC(IGMDCRIT) 
  LMFP = LMFP0*JNSLN/((1.0_PREC-FV)**(2.0_PREC/3.0_PREC)) ! Mpc 

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
  m_ism = 0.0_prec
  m_igm = rho_baryon 
  ! print *, z, m_igm, m_str, m_ism 

  xism_fe = 0.0_prec 
  xism_c = 0.0_prec 
  xism_o = 0.0_prec 

  xigm_fe = 0.0_prec 
  xigm_c = 0.0_prec 
  xigm_o = 0.0_prec 

  metarr(1) = 0.0_prec 

  open(unit=31, file='abundances.out', status='unknown', action='write')
  open(unit=32, file='masses.out', status='unknown', action='write')
  open(unit=33, file='rates.out', status='unknown', action='write')
  open(unit=34, file='sfr.out', status='unknown', action='write')

  r1_igm = 0.0_prec
  r1_ism = 0.0_prec 
  r1_str = 0.0_prec 

  r1igm_fe = 0.0_prec
  r1igm_c = 0.0_prec
  r1igm_o = 0.0_prec

  r1ism_fe = 0.0_prec
  r1ism_c = 0.0_prec
  r1ism_o = 0.0_prec 

  !------------------------------

  DO 
     COUNTR = COUNTR + 1 
     Z = Z+DZ 
     IF (Z<FINAL_REDSHIFT) EXIT 

     ! CALL SOURCEOV(SOURCE, Z) ! [SOURCE] = yr^-1 Mpc^-3 (10^10 M_solar)
     CALL SOURCEOV_POP2(SOURCE_POP2, Z) ! [SOURCE] = yr^-1 Mpc^-3 (10^10 M_solar)
     CALL SOURCEOV_POP3(SOURCE_POP3, Z) ! [SOURCE] = yr^-1 Mpc^-3 (10^10 M_solar)
     SOURCE = SOURCE_POP2 + SOURCE_POP3 
     ! print *, z, q, source*1.0e10_prec, source_pop2*1.0e10_prec, source_pop3*1.0e10_prec 

     SFRARR(COUNTR-1) = SOURCE*1.0E10_PREC ! M_solar yr^-1 Mpc^-3
     NPHDOT = FESC*(SOURCE_POP2*NGAMMA_POP2 + SOURCE_POP3*NGAMMA_POP3) ! yr^-1 Mpc^-3

     !-------------------------

     RECLAMBDA = 2.0_PREC*THI/TEMPH 
     ALPHA_R = 1.269E-13_PREC * RECLAMBDA**1.503_PREC / &
          (1.0_PREC+(RECLAMBDA/0.522_PREC)**0.470_PREC)&
          &**1.923_PREC ! cm^3 s^-1 
     DOFZ = GRFARR(COUNTER(Z))
     NH = RHO_BARYON*(1.0_PREC+RGNOVD*DOFZ)&
          &*(1.1891E57_PREC*1.0E10_PREC) ! Mpc^-3 
     NH_PROPER = NH*(1.0_PREC+Z)**3 ! Mpc^-3      

     !-------------------------

     ! Update ionized hydrogen fraction in ionized region. 
     ! ALPHA_R is hydrogen recombination coefficient (case A).
     GAMMA_REC = R*NH_PROPER*ALPHA_R*CMBYMPCCB*YRBYS ! yr^-1 
     GAMMA_ION = (GPI_POP2*SOURCE_POP2+GPI_POP3*SOURCE_POP3)*&
          &FESC*LMFP*(1.0_PREC+Z)**3*(CMBYMPC**2)/Q ! yr^-1 

     A = GAMMA_REC*DZ*DTDZ(Z)
     B = 1.0_PREC + DZ*DTDZ(Z)*GAMMA_ION/Q
     C = -(OLDXII + DZ*DTDZ(Z)*GAMMA_ION/Q)

     QE = -0.5_PREC*(B+SIGN(DSQRT(B**2-4.0_PREC*A*C),B))
     R1 = QE/A; R2 = C/QE 

     OLDXII = X_II 
     IF (R1<=1.0_PREC .AND. R1>=0.0_PREC) THEN
        X_II = R1
     ELSE 
        X_II = R2 
     END IF

     !-------------------------

     ! Calculate total cooling rate.  We include only recombination 
     ! cooling (case A).
     NHII = NH*OLDXII ! Mpc^-3 
     RECLAMBDA = 2.0_PREC*THI/TEMPH 

     ! RECRATE is the recombination cooling rate for HII 
     RECRATE = 1.778E-29_PREC*TEMPH*RECLAMBDA**1.965_PREC/&
          &(1.0_PREC+(RECLAMBDA/0.541_PREC)**0.502_PREC)&
          &**2.697_PREC ! erg cm^3 s^-1

     GAMMA_RC = (NHII**2)*RECRATE*ERG2J*CMBYMPCCB*YRBYS ! J Mpc^-3 yr^-1
     GAMMA_TOTC = R*GAMMA_RC!*(1.0_PREC+Z)**6 ! J Mpc^-3 yr^-1

     NHI = NH*(1.0_PREC-OLDXII) 
     GAMMAPH = (GPH_POP2*SOURCE_POP2+GPH_POP3*SOURCE_POP3)*FESC*&
          &LMFP*(1.0_PREC+Z)**3*(CMBYMPC**2)  ! J/yr 

     GAMMA_HEAT = GAMMAPH*NHI*(1.0_PREC+Z)**3/Q ! J Mpc^-3 yr^-1 

     ! Calculate photoheating rate.
     HEATCOOL = 2.0_PREC*(GAMMA_HEAT-GAMMA_TOTC)/&
          &(3.0_PREC*NH_PROPER*KBOLTZ) ! K yr^-1 
     TEMPH = (TEMPH + DZ*DTDZ(Z)*HEATCOOL) / &
          &(1.0_PREC + 2.0_PREC*HUBP(Z)*DZ*DTDZ(Z)+ &
          &((X_II-OLDXII)/(1.0_PREC+OLDXII)))
     VCIRC = SQRT(2.0_PREC*KBOLTZ*TEMPH/MPROTON)
     TEMPC = TEMPC+DZ*DTEMPCDZ(Z) ! K

     !------------------------- 

     JNSLN = JEANS_LENGTH(TEMPH,Z) ! Mpc 
     SIGMAB = SIGMA_BARYON(Z, TEMPH) 
     JNSM = (4.0_PREC*PI*RHO_BARYON*JNSLN**3)/3.0_PREC ! 10^10 M_solar
     JMHARR(COUNTR) = JNSM

     IF (PREOVERLAP) THEN 
        R = CLUMPFAC(IGMDCRIT)
        OLDFM = F_M 
        F_M = IGMFPREO()
        FV = IGMVFRAC(IGMDCRIT) 
        LMFP = LMFP0*JNSLN/((1.0_PREC-FV)**(2.0_PREC/3.0_PREC)) ! Mpc 

        ! Update filling factor 
        DFM = F_M - OLDFM 
        T1 = DZ*DTDZ(Z)*NPHDOT/(NH*F_M)
        T2 = DZ*DTDZ(Z)*ALPHA_R*YRBYS*CMBYMPCCB*NH*X_II*R*(1.0_PREC+Z)**3/F_M
        T2 = T2+DFM/F_M
        Q = (Q+T1)/(1.0_PREC+T2)
        FFARR(COUNTR) = Q 
     END IF

     !-------------------------

     IF (Q>=1.0_PREC) THEN 
        PREOVERLAP = .FALSE. 
        POSTOCOUNTER = POSTOCOUNTER + 1 

        SOLFM_Z = Z 
        SOLFM_IGMDCRIT = IGMDCRIT 
        SOLFM_FF = Q 
        FMLO = F_M 
        FMHI = 0.99999999999_PREC 
        TOLZIN = 1.0E-15_PREC 
        OLDIGMFRAC = F_M
        F_M = RTBIS(SOLFM, FMLO, FMHI, TOLZIN, 'main-fm') 

        Q = 1.0_PREC 
        FFARR(COUNTR) = Q 
        SOLDLT_Z = Z
        SOLDLT_IGMFRAC = F_M 
        DLTHI = IGMDCRIT+1.0E30_PREC
        DLTLO = IGMDCRIT-59.0_PREC
        TOLZIN = 1.0E-4_PREC
        DLT = RTBIS(SOLDLT, DLTLO, DLTHI, TOLZIN, 'main-dlt')
        IGMDCRIT = DLT 

        R = CLUMPFAC(IGMDCRIT)
        FV = IGMVFRAC(IGMDCRIT) 
        LMFP = LMFP0*JNSLN/((1.0_PREC-FV)**(2.0_PREC/3.0_PREC)) ! Mpc 
     END IF

     INCLUDE 'volav.inc'

     TAU2 = TAUINT(FV, Q, X_IIVA, Z) 
     TAU = TAU + 0.5_PREC*DZ*(TAU1+TAU2) 
     TAU1 = TAU2 

     DNLLDZ = SPEED_OF_LIGHT*CMBYMPC*YRBYS/&
          &(SQRT(PI)*LMFP*HUBP(Z)*(1.0_PREC+Z)) ! dimensionless 

     GAMMAPI = (GPI_POP2*SOURCE_POP2 + GPI_POP3*SOURCE_POP3)*FESC*LMFP*&
          &(1.0_PREC+Z)**3*(CMBYMPC**2)/YRBYS ! s^-1

!!$     IF (Z == ZLUM) THEN 
!!$        MAGFOO = -20.0_PREC 
!!$        ZFOO = Z 
!!$        DO 
!!$           IF (MAGFOO > -8.0_PREC) EXIT 
!!$           CALL LUMFN(1, MAGFOO, ZFOO, NFOO1)
!!$           CALL LUMFN(2, MAGFOO, ZFOO, NFOO2)
!!$           PRINT *, Z, MAGFOO, (1.0-Q)*NFOO1+Q*NFOO2
!!$           MAGFOO = MAGFOO+0.1_PREC
!!$        END DO
!!$     END IF

     ofl = outflow(z) ! 10^10 M_solar yr^-1 Mpc^-3 
     acc = accrate(z) ! 10^10*M_solar yr^-1 Mpc^-3
     ejc = ejrate(z,0) ! 10^10 M_solar yr^-1 Mpc^-3 

     r2_igm = ofl-acc
     source_dummy = sfr_rollinde_pop2(z)*1.0e-10_prec ! 10^10 M_solar Mpc^-3
     r2_str = source_dummy-ejc
     r2_ism = - r2_igm - r2_str

!!$     dm_igm = (ofl - acc)*dz*dtdz(z) ! 10^10 M_solar Mpc^-3 
!!$     source_dummy = sfr_rollinde_pop2(z)*1.0e-10_prec ! 10^10 M_solar Mpc^-3
!!$     dm_str = (source_dummy - ejc)*dz*dtdz(z) ! 10^10 M_solar Mpc^-3 
!!$     dm_ism = - dm_igm - dm_str ! 10^10 M_solar Mpc^-3 
!!$     
!!$     m_igm = m_igm + dm_igm ! 10^10 M_solar Mpc^-3 
!!$     m_str = m_str + dm_str ! 10^10 M_solar Mpc^-3 
!!$     m_ism = m_ism + dm_ism ! 10^10 M_solar Mpc^-3 

     m_igm = m_igm + 0.5_prec*dz*dtdz(z)*(r2_igm+r1_igm) ! 10^10 M_solar Mpc^-3 
     m_ism = m_ism + 0.5_prec*dz*dtdz(z)*(r2_ism+r1_ism) ! 10^10 M_solar Mpc^-3 
     m_str = m_str + 0.5_prec*dz*dtdz(z)*(r2_str+r1_str) ! 10^10 M_solar Mpc^-3 

     r1_igm = r2_igm; r1_ism = r2_ism; r1_str = r2_ism 

     r2igm_fe = ofl*(xism_fe-xigm_fe)/m_igm
     r2igm_c = ofl*(xism_c-xigm_c)/m_igm
     r2igm_o = ofl*(xism_o-xigm_o)/m_igm
     
!!$     xigm_fe = xigm_fe + dz*dtdz(z)*ofl*(xism_fe-xigm_fe)/m_igm ! dimensionless
!!$     xigm_c = xigm_fe + dz*dtdz(z)*ofl*(xism_c-xigm_c)/m_igm ! dimensionless
!!$     xigm_o = xigm_fe + dz*dtdz(z)*ofl*(xism_o-xigm_o)/m_igm ! dimensionless

     xigm_fe = xigm_fe + 0.5_prec*dz*dtdz(z)*(r1igm_fe+r2igm_fe) ! dimensionless
     xigm_c = xigm_fe + 0.5_prec*dz*dtdz(z)*(r1igm_c+r2igm_c) ! dimensionless
     xigm_o = xigm_fe + 0.5_prec*dz*dtdz(z)*(r1igm_o+r2igm_o) ! dimensionless

     r1igm_fe = r2igm_fe; r1igm_c = r2igm_c; r1igm_o = r2igm_o 

     r2ism_fe = (acc*(xigm_fe-xism_fe) + (ejrate(z,14)-ejc*xism_fe))/m_ism
     r2ism_c = (acc*(xigm_c-xism_c) + (ejrate(z,7)-ejc*xism_c))/m_ism
     r2ism_o = (acc*(xigm_o-xism_o) + (ejrate(z,9)-ejc*xism_o))/m_ism 

     ! Below, 14=Fe, 7=C, and 9=O.  These are row numbers from the
     ! data tables of yields.  See comment in ejrate.f90.  Note: All
     ! three quantities calculated below are DIMENSIONLESS.
!!$     xism_fe = xism_fe + dz*dtdz(z)*(acc*(xigm_fe-xism_fe) + (ejrate(z,14)-ejc*xism_fe))/m_ism
!!$     xism_c = xism_c + dz*dtdz(z)*(acc*(xigm_c-xism_c) + (ejrate(z,7)-ejc*xism_c))/m_ism
!!$     xism_o = xism_o + dz*dtdz(z)*(acc*(xigm_o-xism_o) + (ejrate(z,9)-ejc*xism_o))/m_ism 

     xism_fe = xism_fe + dz*dtdz(z)*(r1ism_fe+r2ism_fe)
     xism_c = xism_c + dz*dtdz(z)*(r1ism_c+r2ism_c)
     xism_o = xism_o + dz*dtdz(z)*(r1ism_o+r2ism_o)

     r1ism_fe = r2ism_fe; r1ism_c = r2ism_c; r1ism_o = r2ism_o 

     m_h = 0.76*m_ism ! 10^10 M_solar Mpc^-3 
     m_fe = xism_fe*m_ism ! 10^10 M_solar Mpc^-3 
     m_c = xism_c*m_ism ! 10^10 M_solar Mpc^-3 
     m_o = xism_o*m_ism ! 10^10 M_solar Mpc^-3 

     metmass = m_fe+m_c+m_o ! 10^10 M_solar Mpc^-3 
     metallicity = (metmass/m_ism)/0.02_prec ! dimensionless 
     if (m_ism==0.0_prec) metallicity = 0.0_prec 
     ! metarr(countr) = abs(metallicity*1.0e9_prec)
     metarr(countr) = abs(metallicity)

     if (m_fe==0.0_prec) then 
        fe_abundance = 0.0_prec 
     else 
        fe_abundance = log10(abs(m_fe/m_h)) + 5.0 
     end if

     if (m_c==0.0_prec) then 
        c_abundance = 0.0_prec 
     else 
        c_abundance = log10(abs(m_c/m_h)) + 4.5
     end if

     if (m_o==0.0_prec) then 
        o_abundance = 0.0_prec 
     else 
        o_abundance = log10(abs(m_o/m_h)) + 4.0
     end if

     write (31, *) z, metarr(countr), fe_abundance, c_abundance, o_abundance  
     write (32, *) z, m_igm, m_str, m_ism, metmass 
     write (33, *) z, acc, ofl, ejc
     write (34, *) z, sfr_rollinde_pop2(z), sfr_rollinde_pop3(z), sfr_rollinde_pop2(z)+sfr_rollinde_pop3(z)
  END DO

  close(31)
  close(32)
  close(33)
  close(34)

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

END PROGRAM REION

