
MODULE CONSTANTS

  ! File: constants.f90
  ! Cre: 2010-04-07
  ! Mod: $Date: 2012/10/23 08:10:11 $; $Revision: 1.25 $
  ! 
  ! This file contains all cosmology-related 
  ! global constants that are required for the 
  ! reionization code.  Cosmological parameters 
  ! are WMAP5 best-fit. 
  !
  ! ,----
  ! | REMEMBER TO CHANGE tabgrf.dat AND tabsig.dat 
  ! | IF ANYTHING IS CHANGED HERE.  
  ! `----

  IMPLICIT NONE

  ! Set precision for the whole calculation.  PREC = 4 is usually
  ! single precision, 8 is double and 16 quad.
  INTEGER, PARAMETER :: PREC = 8

  ! These values are from J. Dunkley et al., 
  ! ApJS 180, 306 (2009).
  REAL(PREC), PARAMETER :: OMEGA_B = 4.49E-2_PREC, & ! dimensionless 
       &OMEGA_NR = 0.2646E0_PREC, & ! dimensionless 
       &OMEGA_LAMBDA = 0.734E0_PREC, & ! dimensionless 
       &SPECTRAL_INDEX = 0.963E0_PREC, & ! dimensionless 
       &OPTICAL_DEPTH = 8.8E-2_PREC, & ! dimensionless 
       &SIGMA_EIGHT = 0.801E0_PREC, & ! dimensionless
       &SMALLH = 0.710E0_PREC, & ! dimensionless 
       &OMEGA_K = 1.0E0_PREC-(OMEGA_LAMBDA+OMEGA_NR)  ! dimensionless 
  REAL(PREC) :: RHO_CRITICAL = 2.775E1*SMALLH**2, & ! 10^10 M_solar/Mpc^3
       &HUBBLE_0 = 1.023E-10_PREC*SMALLH ! yr^-1

  !---------------------------

  ! These values are from Wolfram Alpha 
  ! and Arthur N. Cox, Allen's astrophysical
  ! quantities (New York: Springer-Verlag, 
  ! 2000).
  REAL(PREC), PARAMETER :: &
       &YRBYS = 3.154E7_PREC, & ! yr/s conversion factor 
       &CMBYMPCCB = 3.4036771916E-74_PREC, & ! (cm/Mpc)^3 conversion factor 
       &CMBYMPC = 3.24077928965E-25_PREC, & ! cm/Mpc conversion factor 
       &MSOLKG = 1.988435E30_PREC, & ! solar mass; kg
       &MPROTON = 1.672622E-27_PREC, & ! mass of proton; kg
       &KBOLTZ = 1.38065E-23_PREC, & ! Boltzmann constant; J/K 
       &SPEED_OF_LIGHT = 2.998E10_PREC, & ! cm/s
       &NEWTG = 6.673E-11_PREC, & ! gravitational constant; m^3kg^-1s^-2   
       &EPS0 = 2.179E-18_PREC, & ! threshold energy for H ionization; J 
       &NU0 = 3.288E15_PREC, & ! threshold freq for H ionization; s^-1 (Hz)
       &E0 = 13.6_PREC, & ! threshold energy for H ionization; eV
       &DELTA_VIRIAL = 178.0_PREC, &! dimensionless; Virial overdensity
       &DELTAC = 1.686_PREC, &!dimensionless; critical overdensity 
       &PI = 3.14159265358979_PREC, &
       &CMBTEMP = 2.725_PREC, & ! K; CMB temperature today 
       &ERG2J = 1.0E-7_PREC, & ! erg/J conversion factor 
       &EVBYK = 11604.505_PREC, & ! eV/K conversion factor 
       &THI = 1.57807E5, & ! K; H ionization threshold
       &HPLANCK = 6.626069E-34, &! Js; Planck constant 
       &CMBYANG = 1.0E8 ! cm/Angstrom conversion factor 

  !---------------------------

  ! These values change behaviour of code. 
  REAL(PREC), PARAMETER :: &
       &INITIAL_REDSHIFT = 50.0_PREC, & ! when code begins run 
       &DZ = -0.5_PREC, & ! redshift stepsize (also used 
                                ! in places other than main.f90) 
       &FINAL_REDSHIFT = 0.0_PREC, & ! when code stops running 
       &ABSERR = 1.0E-2_PREC, &
       &ABSREL = 1.0E-2_PREC, &
       &TABGRF_ZINITIAL = 0.0_PREC, &
       &TABGRF_DZ = 0.005_PREC, &
       &FSTAR = 0.002E0_PREC, &
       &FSTAR_POP3 = 0.002E0_PREC, &
       &IGMDCRIT_PREO = 60.0_PREC, & ! dimensionless; IGMDCRIT_PREO is the 
                                ! critical density of IGM in the pre-overlap stage.  Only 
                                ! regions with density less than this will be ionized in that
                                ! phase.  See Sec 2.2 of CF05.
       &LMFP0 = 1.7e-3_PREC, & ! Mpc; LMFP0 is a constant factor that 
                                ! multiplies the mean free path to calibrate it with observations.  
                                ! See Sec 5.1 of CF05.
       &ABATIC_GAMMA = 1.0_PREC, & ! dimensionless; adiabatic index for 
                                ! baryons
       &BURST_TOTAL_MASS = 1.0E-4_PREC ! 10^10 M_solar; total mass used 
  ! in starburst population synthesis.

  INTEGER, PARAMETER :: KEY = 6, & ! KEY and MAXINT control DQAGE
       &MAXINT = 200 

  REAL(KIND=PREC) :: ZLUM, & ! Redshift at which luminosity function is to be 
                                ! calculated.
       &RGNOVD, &! dimensionless; RGNOVD is overdensity 
                                ! of the region in which the model is to be solved.  RGNSIZE 
                                ! is the region's comoving Lagrangian size.
       &FESC, & 
       &RGNSIZE ! Mpc (7.746)

  REAL(KIND=PREC) :: MINF = 0.1_PREC, & ! M_solar; Minimum stellar mass (for Pop 2 IMF) 
       MSUP = 100.0_PREC, & ! M_solar; Maximum stellar mass (for Pop 2 IMF) 
       IMF_SLOPE = 1.3_PREC, & ! (Slope of IMF - 1), 1.3 for Salpeter, 1.7 for Scalo
       OUTFLOW_EFFICIENCY = 1.0e-6_PREC, & ! See outflow.f90 
       MINF_POP3 = 100.0_PREC, & ! M_solar; Minimum stellar mass (for Pop 3 IMF) 
       MSUP_POP3 = 260.0_PREC, & ! M_solar; Maximum stellar mass (for Pop 3 IMF) 
!!$       MINF_POP3 = 1.0_PREC, & ! M_solar; Minimum stellar mass (for Pop 3 IMF) 
!!$       MSUP_POP3 = 100.0_PREC, & ! M_solar; Maximum stellar mass (for Pop 3 IMF) 
       DATA_MMIN = 12.07_PREC, & ! M_Solar; Minimum stellar mass in yields data time 
       DATA_MMAX = 266.0_PREC, & ! M_Solar; Maximum stellar mass in yields data time 
     ! DATA_MMAX = 40.22_PREC, & ! M_Solar; Maximum stellar mass in yields data time 
       DATA_ZMETMIN = 0.0_PREC, & ! dimensionless; Minimum Z in yields data 
       DATA_ZMETMAX = 0.01_PREC, & ! dimensionless; Maximum Z in yields data 
       REDSHIFT_POP3TRANS = 60.0_PREC, & ! redshift at which Pop III SFR stops and Pop II begins
       METALLICITY_POP3TRANS = 1.0e-4_prec, & ! metalicity at which Pop III SFR stops and Pop II begins
                                ! STMASS_UPLIMIT = 130.0_prec, & ! M_solar 
       STMASS_UPLIMIT = 260.0_prec, & ! M_solar 
       POP2_UPLIMIT = 130.0_prec, & ! M_solar
       STELLAR_INTEGRAL_MULTIPLIER = 1.25_prec, & 
!!$       Enrich_time_lag = 0.0 ! yr ; used in haloyield_nonira.f90 
       Enrich_time_lag = 0.0e0_prec, & ! yr ; used in haloyield_nonira.f90 
       Enrich_time_lag2 = 1.0e12_prec 

  INTEGER, PARAMETER :: DATA_COLUMNS = 25

END MODULE CONSTANTS

