  
  ! Read population synthesis data.  This is from 
  ! Starburst99 output `reion-generic'.  nnudot.dat 
  ! is written by popsyn/intspectrum.f90.  
  OPEN(11, FILE='/home/girish/reion-eq/nnudot.dat', STATUS='OLD', ACTION='READ') 
  READ(11, *) FREQ, NDOTM  
  NUMBER_OF_LINES = 0 
  DO 
     NUMBER_OF_LINES = NUMBER_OF_LINES + 1 
     READ(11, *, END=911) FREQ, NDOTM 
  END DO
911 CONTINUE 
  REWIND 11 

  ALLOCATE(IONNNUM(NUMBER_OF_LINES)) 
  ALLOCATE(IONFREQ(NUMBER_OF_LINES)) 
  DO I = 1, NUMBER_OF_LINES 
     READ(11,*) FREQ, NDOTM 
     IONFREQ(I) = FREQ ! Hz 
     IONNNUM(I) = 1.0E10_PREC*NDOTM ! /Hz/(10^10 M_solar) 
  END DO
  CLOSE(11) 

  OPEN(11, FILE='/home/girish/reion-eq/dndm_PopIII_salpeter_star', STATUS='OLD', ACTION='READ') 
  READ(11, *) FREQ, NDOTM  
  NUMBER_OF_LINES = 0 
  DO 
     NUMBER_OF_LINES = NUMBER_OF_LINES + 1 
     READ(11, *, END=916) FREQ, NDOTM 
  END DO
916 CONTINUE 
  REWIND 11 

  ALLOCATE(IONNUMPOP3(NUMBER_OF_LINES)) 
  ALLOCATE(IONFREQPOP3(NUMBER_OF_LINES)) 
  DO I = 1, NUMBER_OF_LINES 
     READ(11,*) FREQ, NDOTM 
     IONFREQPOP3(I) = FREQ ! Hz 
     IONNUMPOP3(I) = 1.0E10_PREC*NDOTM ! /Hz/(10^10 M_solar) 
  END DO
  CLOSE(11) 

  !--------------------------

  ! Read in solution to Freidmann equation for the background 
  ! cosmology.  This gives a redshift--cosmic time conversion
  ! table.  age.dat and z.dat are produced using bg_cosmo.f90.
  OPEN(11, FILE='/home/girish/reion-eq/age.dat', STATUS='OLD', ACTION='READ')
  OPEN(12, FILE='/home/girish/reion-eq/z.dat', STATUS='OLD', ACTION='READ') 

  ! Read in age of the universe from age.dat and count number
  ! of lines in z.dat.  
  READ(11, *) H0T0, T0SEC, T0YR
  NUMBER_OF_LINES = 0
  DO
     READ(12, *, END=912) HTT, ZOFT
     NUMBER_OF_LINES = NUMBER_OF_LINES + 1 
  END DO
912 CONTINUE
  REWIND 12

  ! Allocate and fill arrays TARR and ZARR with values of time
  ! and redshift.  
  ALLOCATE(TARR(NUMBER_OF_LINES))
  ALLOCATE(ZARR(NUMBER_OF_LINES))
  DO I = 1, NUMBER_OF_LINES
     READ(12, *) HTT, ZOFT 
     LOCAL_HUBBLE_0 = 1.023E-10_PREC*0.719_PREC ! yr^-1
     TARR(I) = T0YR-HTT/LOCAL_HUBBLE_0 ! yr 
     ZARR(I) = ZOFT
  END DO
  CLOSE(11)
  CLOSE(12)

  !--------------------------

  ! Read in values of mass, sigma and dsigma/dmass from 
  ! tabsig.dat.  tabsig.dat is produced using tabsig.f90. 
  OPEN(11, FILE='/home/girish/reion-eq/tabsig.dat', STATUS='OLD', ACTION='READ') 

  ! Count number of lines in tabsig.dat. 
  NUMBER_OF_LINES = 0
  DO
     READ(11, *, END=913) M, SIG, DSIG 
     NUMBER_OF_LINES = NUMBER_OF_LINES + 1
  END DO
913 CONTINUE
  REWIND 11 

  ! Allocate arrays and store data.  MSARR holds mass, 
  ! SIGMARR holds sigma and DSIGMARR holds the mass 
  ! derivative of sigma.  This stores values of SIGMA 
  ! for mass between 10^4 and 10^16 M_solar.
  ALLOCATE(MSARR(NUMBER_OF_LINES))
  ALLOCATE(SIGMARR(NUMBER_OF_LINES))
  ALLOCATE(DSIGMARR(NUMBER_OF_LINES))
  DO I = 1, NUMBER_OF_LINES
     READ(11, *) M, SIG, DSIG
     MSARR(I) = M ! 10^10 M_solar 
     SIGMARR(I) = SIG ! dimensionless 
     DSIGMARR(I) = DSIG ! (10^10 M_solar)^-1 
  END DO
  CLOSE(11) 

  !--------------------------

  ! Read in values of the growth factor and its redshift
  ! derivative against redshift.  tabgrf.f90 is produced 
  ! using tabgrf.f90.  
  OPEN(11, FILE='/home/girish/reion-eq/tabgrf.dat', STATUS='OLD', ACTION='READ') 

  ! Count number of lines in tabgrf.dat.
  NUMBER_OF_LINES = 0
  DO
     READ(11, *, END=914) Z, GRF, DGRF
     NUMBER_OF_LINES = NUMBER_OF_LINES + 1
  END DO
914 CONTINUE
  REWIND 11

  ! Allocate and fill arrays with data.  RSARR 
  ! holds redshift values, GRFARR holds the 
  ! growth factor and DGRFARR holds its redshift
  ! derivative.
  ALLOCATE(RSARR(NUMBER_OF_LINES))
  ALLOCATE(GRFARR(NUMBER_OF_LINES))
  ALLOCATE(DGRFARR(NUMBER_OF_LINES))
  DO I = 1, NUMBER_OF_LINES
     READ(11, *) Z, GRF, DGRF 
     RSARR(I) = Z ! dimensionless 
     GRFARR(I) = GRF ! dimensionless 
     DGRFARR(I) = DGRF ! dimensionless 
  END DO
  CLOSE(11) 

  !--------------------------

  OPEN(11, FILE='/home/girish/reion-eq/popsyn/l1500-spectrum-lt.dat', STATUS='OLD', ACTION='READ')

  ! Count number of lines in file.
  NUMBER_OF_LINES = 0
  DO
     READ(11, *, END=915) LB_AGE, LB_LAMBDA, LB_LUM
     NUMBER_OF_LINES = NUMBER_OF_LINES + 1
  END DO
915 CONTINUE
  REWIND 11

  ALLOCATE(LBURST(2,NUMBER_OF_LINES))
  DO I = 1, NUMBER_OF_LINES 
     READ(11, *) LB_AGE, LB_LAMBDA, LB_LUM
     LBURST(1,I) = LB_AGE ! yr 
     LBURST(2,I) = (10.0_PREC**LB_LUM)/BURST_TOTAL_MASS ! erg/s/Ang/(10^10 M_solar)
  END DO
  CLOSE(11)

  !--------------------------

  OPEN(11, FILE='/home/girish/reion-eq/data/lifetime.mm89.dat', STATUS='OLD', ACTION='READ')
  ALLOCATE(STELLAR_MASS(24)) 
  ALLOCATE(STELLAR_AGE(24)) 
  DO I = 1, 24 
     READ(11, *) ST_MASS, ST_AGE 
     STELLAR_MASS(I) = ST_MASS ! M_solar (NOT 10^10 M_solar!) 
     STELLAR_AGE(I) = ST_AGE ! Myr (NOT yr!) 
  END DO
  CLOSE(11)

  !--------------------------

  OPEN(11, FILE='/home/girish/reion-eq/data/pop3_ngamma_schaerrer.dat', STATUS='OLD', ACTION='READ')
  ALLOCATE(POP3_STELLAR_MASS(12)) 
  ALLOCATE(POP3_STELLAR_AGE(12)) 
  ALLOCATE(POP3_STELLAR_NGAMMA(12)) 
  DO I = 1, 12
     READ(11, *) ST_MASS, ST_AGE, ST_NGAMMA 
     POP3_STELLAR_MASS(I) = ST_MASS ! M_solar (NOT 10^10 M_solar!) 
     POP3_STELLAR_AGE(I) = ST_AGE ! Myr (NOT yr!) 
     POP3_STELLAR_NGAMMA(I) = ST_NGAMMA ! photons/s 
  END DO
  CLOSE(11)

  !--------------------------

  ALLOCATE(YIELDS(4,DATA_COLUMNS,17))

  YIELDS(1,1,1) = 0.0_prec
  YIELDS(2,1,1) = 0.0001_prec
  YIELDS(3,1,1) = 0.01_prec
  YIELDS(4,1,1) = 0.1_prec

  OPEN(11, FILE='/home/girish/reion-eq/data/yields.wwhw.Z0', STATUS='OLD', ACTION='READ')
  DO I = 2, 17
     READ(11, *) YIELD_RECORD 
     DO J = 1, DATA_COLUMNS 
        YIELDS(1,J,I) = YIELD_RECORD(J) 
     END DO
  END DO
  CLOSE(11)

  OPEN(11, FILE='/home/girish/reion-eq/data/yields.wwhw.Z00001', STATUS='OLD', ACTION='READ')
  DO I = 2, 17
     READ(11, *) YIELD_RECORD 
     DO J = 1, DATA_COLUMNS 
        YIELDS(2,J,I) = YIELD_RECORD(J) 
     END DO
  END DO
  CLOSE(11)

  OPEN(11, FILE='/home/girish/reion-eq/data/yields.wwhw.Z001', STATUS='OLD', ACTION='READ')
  DO I = 2, 17
     READ(11, *) YIELD_RECORD 
     DO J = 1, DATA_COLUMNS 
        YIELDS(3,J,I) = YIELD_RECORD(J) 
     END DO
  END DO
  CLOSE(11)

  OPEN(11, FILE='/home/girish/reion-eq/data/yields.wwhw.Z01', STATUS='OLD', ACTION='READ')
  DO I = 2, 17
     READ(11, *) YIELD_RECORD 
     DO J = 1, DATA_COLUMNS 
        YIELDS(4,J,I) = YIELD_RECORD(J) 
     END DO
  END DO
  CLOSE(11)

  !-------------------------

  allocate(sd93_tvir(181), sd93_coolrate(181)) 
  ! open(11, file='/home/girish/reion-eq/data/SD93_nie_zero.dat', status='old', action='read') 
  open(11, file='/home/girish/reion-eq/data/SD93_nie_metals.dat', status='old', action='read') 
  do i = 1, 181 
     read(11, *) sdt, sdl 
     sd93_tvir(i) = sdt ! K 
     sd93_coolrate(i) = sdl ! J m^3 s^-1 
  end do

  !-------------------------

  ! Read normalised cooling curve data from Sutherland and Dopita for
  ! various metallicities.  Remember, we are using the CIE data.
  allocate(sd93_temp(90), sd93_febyh(8), sd93_lambda(90,8))

  j = 1 
  open(11, file='/home/girish/reion-eq/data/SD93_zero.cie', status='old', action='read') 
  sd93_febyh(j) = -10.0_prec ! Just a very small number (supposed to be -Infinity).
  do i = 1, 90 
     read (11, *) temp_sd93, lambda_sd93 
     sd93_temp(i) = temp_sd93 
     sd93_lambda(i,j) = lambda_sd93 
  end do
  close(11) 

  ! No need to write anything onto temp_sd93 again as that column is
  ! same in all files.
  j = 2 
  open(11, file='/home/girish/reion-eq/data/SD93-30.cie', status='old', action='read') 
  sd93_febyh(j) = -3.0_prec 
  do i = 1, 90 
     read (11, *) temp_sd93, lambda_sd93 
     sd93_lambda(i,j) = lambda_sd93 
  end do
  close(11) 

  j = 3
  open(11, file='/home/girish/reion-eq/data/SD93-20.cie', status='old', action='read') 
  sd93_febyh(j) = -2.0_prec 
  do i = 1, 90 
     read (11, *) temp_sd93, lambda_sd93 
     sd93_lambda(i,j) = lambda_sd93 
  end do
  close(11) 

  j = 4
  open(11, file='/home/girish/reion-eq/data/SD93-15.cie', status='old', action='read') 
  sd93_febyh(j) = -1.5_prec 
  do i = 1, 90 
     read (11, *) temp_sd93, lambda_sd93 
     sd93_lambda(i,j) = lambda_sd93 
  end do
  close(11) 

  j = 5
  open(11, file='/home/girish/reion-eq/data/SD93-10.cie', status='old', action='read') 
  sd93_febyh(j) = -1.0_prec 
  do i = 1, 90 
     read (11, *) temp_sd93, lambda_sd93 
     sd93_lambda(i,j) = lambda_sd93 
  end do
  close(11) 

  j = 6
  open(11, file='/home/girish/reion-eq/data/SD93-05.cie', status='old', action='read') 
  sd93_febyh(j) = -0.5_prec 
  do i = 1, 90 
     read (11, *) temp_sd93, lambda_sd93 
     sd93_lambda(i,j) = lambda_sd93 
  end do
  close(11) 

  j = 7
  open(11, file='/home/girish/reion-eq/data/SD93-00.cie', status='old', action='read') 
  sd93_febyh(j) = 0.0_prec 
  do i = 1, 90 
     read (11, *) temp_sd93, lambda_sd93 
     sd93_lambda(i,j) = lambda_sd93 
  end do
  close(11) 

  j = 8
  open(11, file='/home/girish/reion-eq/data/SD93+05.cie', status='old', action='read') 
  sd93_febyh(j) = 0.5_prec 
  do i = 1, 90 
     read (11, *) temp_sd93, lambda_sd93 
     sd93_lambda(i,j) = lambda_sd93 
  end do
  close(11) 



