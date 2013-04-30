
FUNCTION LUMINOSITY(FORMATION_REDSHIFT, TFORM, T)

  ! File: luminosity.f90
  ! Cre: 2010-04-07

  ! Mod: $Date: 2010/05/17 13:23:22 $; $Revision: 1.5 $ 
  !
  ! This function returns value of luminosity/(ngamma*eps0), 
  ! where `luminosity' is energy released at redshift REDSHIFT 
  ! (time T) per unit time from a halo of mass MASS that formed 
  ! at redshift FORMATION_REDSHIFT (time TFORM), `ngamma' is 
  ! energy released per unit time per unit mass (see srcov.f90) 
  ! in stars, and eps0 is 13.6 eV (see constants.f90).  

  USE CONSTANTS; IMPLICIT NONE
  REAL(KIND = PREC), INTENT(IN) :: FORMATION_REDSHIFT, &
       &TFORM, T 
  REAL(KIND = PREC) :: LUMINOSITY

  REAL(KIND = PREC) :: DYNAMICAL_TIME, TIME, FORMATION_TIME,&
       &TKERNEL

  TIME = T ! yr
  FORMATION_TIME = TFORM ! yr 
  DYNAMICAL_TIME = TDYN(FORMATION_REDSHIFT) ! yr 

  TKERNEL = ((TIME - FORMATION_TIME) / DYNAMICAL_TIME**2) * &
       &EXP(-(TIME - FORMATION_TIME) / DYNAMICAL_TIME) ! yr^-1 

  LUMINOSITY = FSTAR*TKERNEL ! yr^-1 

CONTAINS 

  FUNCTION TDYN(Z)

    REAL(KIND = PREC), INTENT(IN) :: Z
    REAL(KIND = PREC) :: TDYN

    REAL(KIND = PREC) :: HDENS, BKGRHO, RHOCS

    RHOCS = 1.879E-29_PREC * SMALLH ** 2 ! g/cm^3
    BKGRHO = RHOCS * OMEGA_NR * (1.0_PREC + Z) ** 3 ! g/cm^3 
    HDENS = DELTA_VIRIAL * BKGRHO ! g/cm^3 
    TDYN = SQRT(3.0_PREC * PI / (16.0_PREC * NEWTG * &
         &HDENS * 1.0E3_PREC)) / YRBYS ! yr
    ! 1.0E3 converts g to kg and cm to m.  

  END FUNCTION TDYN

END FUNCTION LUMINOSITY

