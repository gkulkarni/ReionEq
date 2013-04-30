SUBROUTINE SOURCEOV_POP2(SOURCE, REDSHIFT)

  ! File: srcov.f90
  ! Cre: 2010-04-07
  ! Mod: $Date: 2012/06/14 13:30:47 $; $Revision: 1.14 $ 
  !  
  ! Calculates the total number of photons per 
  ! hydrogen atom (script-S of CO) produced by sources 
  ! located in a region of size RGNSIZE and overdensity
  ! RGNOVD at redshift REDSHIFT.

  USE CONSTANTS
  USE INTERFACES, ONLY : COUNTER, INTERPOLATE2, NUINT, LUMINOSITY, &
       &GETJMC, GETJMH, GETFF, NUINT_POP2
  USE STORAGE
  IMPLICIT NONE
  REAL(KIND = PREC), INTENT(IN) :: REDSHIFT
  REAL(KIND = PREC), INTENT(OUT) :: SOURCE

  REAL(KIND = PREC) :: DZFORM, ERR, LUM, NH, RGNMS, RHO,&
       &RHO_BARYON, SUM, T1, T2
  INTEGER :: FLAG, LAST, NEVAL

  !---------------------------

  RHO = OMEGA_NR*RHO_CRITICAL ! 10^10 M_solar / Mpc^3 
  RHO_BARYON = OMEGA_B*RHO_CRITICAL ! 10^10 M_solar / Mpc^3 

  SUM = 0.0_PREC; DZFORM = DZ 
  ZFORM = INITIAL_REDSHIFT + 2.0_PREC

  T1 = TIME_INTEGRAND(ZFORM)  
  DO 
     ZFORM = ZFORM + DZFORM 
     IF (ZFORM <= REDSHIFT) EXIT 
     T2 = TIME_INTEGRAND(ZFORM) 

     SUM = SUM + DZFORM * (T1 + T2) * 0.5_PREC
     T1 = T2 
  END DO

  SOURCE = SUM ! yr^-1 Mpc^-3 (10^10 Msolar)

CONTAINS

  FUNCTION TIME_INTEGRAND(Z)

    IMPLICIT NONE
    REAL(KIND = PREC), INTENT(IN) :: Z
    REAL(KIND = PREC) :: TIME_INTEGRAND
    REAL(KIND = PREC) :: DMMASS, MMAX, NUMAX, SUM, TFORM, &
         &T, SIGMA, DSIGMADM, SUM2, OVDEXTRA, NUMINII,&
         &DGROWTH, SURVIVAL_PROBABILITY, MMINI, MMINII, SUM1,&
         &NUMINI, FII
    INTEGER :: GRLOC, loc

    if (z/=zform) write (0,*) Z, ZFORM 

    MMINI = GETJMC(Z)
    MMINII = GETJMH(Z)

    MMAX = RGNMS*0.9_PREC ! 10^10 M_solar
    GRLOC = COUNTER(Z); GROWTH = GRFARR(GRLOC)

    !------------------------

    loc = floor((initial_redshift - z) / abs(dz)) + 1
    if (loc < 1) then
       sum = 0.0_PREC
    else
       sum = nuintegralarr_pop2(loc)
    end if

    !------------------------

    CALL INTERPOLATE2(TARR, ZARR, Z, TFORM) ! [TFORM] = yr
    CALL INTERPOLATE2(TARR, ZARR, REDSHIFT, T) ! [T] = yr 
    LUM = LUMINOSITY(Z, TFORM, T) ! yr^-1

    SURVIVAL_PROBABILITY=(DELTAC/GRFARR(COUNTER(REDSHIFT))-RGNOVD)/&
         &(DELTAC/GROWTH-RGNOVD) ! dimensionless 

    DGROWTH = DGRFARR(GRLOC) 
    DGROWTH = DGROWTH/GROWTH**2 ! dimensionless 

    OVDEXTRA = DELTAC/(DELTAC/GROWTH-RGNOVD) ! dimensionless 

    TIME_INTEGRAND = SUM*LUM*SURVIVAL_PROBABILITY*DGROWTH*&
         &OVDEXTRA*RHO_BARYON*RHO_BARYON ! yr^-1 Mpc^-3 (10^10 M_solar)
  END FUNCTION TIME_INTEGRAND

END SUBROUTINE SOURCEOV_POP2

!----------------------------------------------------------

SUBROUTINE SOURCEOV_POP3(SOURCE, REDSHIFT)

  USE CONSTANTS
  USE INTERFACES, ONLY : COUNTER, INTERPOLATE2, NUINT, LUMINOSITY, &
       &GETJMC, GETJMH, GETFF, NUINT_POP3    
  USE STORAGE
  IMPLICIT NONE
  REAL(KIND = PREC), INTENT(IN) :: REDSHIFT
  REAL(KIND = PREC), INTENT(OUT) :: SOURCE

  REAL(KIND = PREC) :: DZFORM, ERR, LUM, NH, RGNMS, RHO,&
       &RHO_BARYON, SUM, T1, T2
  INTEGER :: FLAG, LAST, NEVAL

  !---------------------------

  RHO = OMEGA_NR*RHO_CRITICAL ! 10^10 M_solar / Mpc^3 
  RHO_BARYON = OMEGA_B*RHO_CRITICAL ! 10^10 M_solar / Mpc^3 

  SUM = 0.0_PREC; DZFORM = DZ 
  ZFORM = INITIAL_REDSHIFT + 2.0_PREC
  T1 = TIME_INTEGRAND(ZFORM)  
  DO 
     ZFORM = ZFORM + DZFORM 
     IF (ZFORM <= REDSHIFT) EXIT 
     T2 = TIME_INTEGRAND(ZFORM) 

     SUM = SUM + DZFORM * (T1 + T2) * 0.5_PREC
     T1 = T2 
  END DO

  SOURCE = SUM ! yr^-1 Mpc^-3 (10^10 Msolar)

CONTAINS

  FUNCTION TIME_INTEGRAND(Z)

    IMPLICIT NONE
    REAL(KIND = PREC), INTENT(IN) :: Z
    REAL(KIND = PREC) :: TIME_INTEGRAND
    REAL(KIND = PREC) :: DMMASS, MMAX, NUMAX, SUM, TFORM, &
         &T, SIGMA, DSIGMADM, SUM2, OVDEXTRA, NUMINII,&
         &DGROWTH, SURVIVAL_PROBABILITY, MMINI, MMINII, SUM1,&
         &NUMINI, FII
    INTEGER :: GRLOC, loc

    if (z/=zform) write (0,*) Z, ZFORM 

    MMINI = GETJMC(Z)
    MMINII = GETJMH(Z)

    MMAX = RGNMS*0.9_PREC ! 10^10 M_solar
    GRLOC = COUNTER(Z); GROWTH = GRFARR(GRLOC)

    !------------------------
    
    loc = floor((initial_redshift - z) / abs(dz)) + 1
    if (loc < 1) then
       sum = 0.0_PREC
    else
       sum = nuintegralarr_pop3(LOC)
    end if

    !------------------------

    CALL INTERPOLATE2(TARR, ZARR, Z, TFORM) ! [TFORM] = yr
    CALL INTERPOLATE2(TARR, ZARR, REDSHIFT, T) ! [T] = yr 
    LUM = LUMINOSITY(Z, TFORM, T) ! yr^-1

    SURVIVAL_PROBABILITY=(DELTAC/GRFARR(COUNTER(REDSHIFT))-RGNOVD)/&
         &(DELTAC/GROWTH-RGNOVD) ! dimensionless 

    DGROWTH = DGRFARR(GRLOC) 
    DGROWTH = DGROWTH/GROWTH**2 ! dimensionless 

    OVDEXTRA = DELTAC/(DELTAC/GROWTH-RGNOVD) ! dimensionless 

    TIME_INTEGRAND = SUM*LUM*SURVIVAL_PROBABILITY*DGROWTH*&
         &OVDEXTRA*RHO_BARYON*RHO_BARYON ! yr^-1 Mpc^-3 (10^10 M_solar)

  END FUNCTION TIME_INTEGRAND

END SUBROUTINE SOURCEOV_POP3

