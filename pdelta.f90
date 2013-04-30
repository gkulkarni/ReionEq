PROGRAM PDELTA

  ! File: pdelta.f90
  ! CRE: 2009-10-24
  ! Mod: $Date: 2010/06/05 05:27:19 $; $Revision: 1.1 $
  !
  ! This code gives a value of the comoving Lagrangian 
  ! overdensity DLTLS.  This is then useful for the feedback 
  ! project.  See Makefile for list of dependencies.  
  ! Compiles with `Make pdelta.'

  USE CONSTANTS; USE INTERFACES; USE STORAGE; IMPLICIT NONE
  REAL(KIND = PREC) :: RE, RL, DLTLS, DLTSS, SGSS, SGLS, Z, RHOBAR
  REAL(KIND = PREC) :: T, RS, GRF, DGRF, GROWTH, M, SIG, DSIG, QLS
  REAL(KIND = PREC) :: QSOMASS, QSS, PD, MASS, RHO
  INTEGER :: GRLOC, NUMBER_OF_LINES, I

  INTERFACE

     FUNCTION THETA(DELTA_L, REDSHIFT)

       USE CONSTANTS; USE STOTABZ; IMPLICIT NONE
       INTERFACE

          FUNCTION FTH(X)
            
            USE CONSTANTS;
            IMPLICIT NONE
            REAL(KIND = PREC), INTENT(IN) :: X
            REAL(KIND = PREC) :: FTH

          END FUNCTION FTH

       END INTERFACE

       REAL(KIND = PREC), INTENT(IN) :: DELTA_L, REDSHIFT
       REAL(KIND = PREC) :: THETA

     END FUNCTION THETA

  END INTERFACE

  !-------------------------------

  ! tabgrf.f90 is produced using tabgrf.f90.  
  OPEN(11, FILE = 'tabgrf.dat', STATUS = 'OLD', ACTION = 'READ') 

  NUMBER_OF_LINES = 0
  DO
     READ(11, *, END = 913) Z, GRF, DGRF
     NUMBER_OF_LINES = NUMBER_OF_LINES + 1
  END DO
913 CONTINUE
  REWIND 11

  ALLOCATE(GRFARR(NUMBER_OF_LINES))
  ALLOCATE(DGRFARR(NUMBER_OF_LINES))

  DO I = 1, NUMBER_OF_LINES
     READ(11, *) Z, GRF, DGRF 
     GRFARR(I) = GRF ! dimensionless 
     DGRFARR(I) = DGRF ! dimensionless 
  END DO

  CLOSE(11) 

  !-------------------------------

  ! tabsig.dat is produced using tabsig.f90. 
  OPEN(11, FILE = 'tabsig.dat', STATUS = 'OLD', ACTION = 'READ') 

  NUMBER_OF_LINES = 0
  DO
     READ(11, *, END = 912) M, SIG, DSIG 
     NUMBER_OF_LINES = NUMBER_OF_LINES + 1
  END DO
912 CONTINUE
  REWIND 11

  ALLOCATE(MSARR(NUMBER_OF_LINES))
  ALLOCATE(MDIFF(NUMBER_OF_LINES))
  MDIFF = 0.0_PREC 
  ALLOCATE(SIGMARR(NUMBER_OF_LINES))
  ALLOCATE(DSIGMARR(NUMBER_OF_LINES))
  ALLOCATE(SIGDIFF(NUMBER_OF_LINES))
  SIGDIFF = 0.0_PREC 

  ! This stores values of SIGMA for mass between 10^4 and 10^16 M_solar.  
  DO I = 1, NUMBER_OF_LINES
     READ(11, *) M, SIG, DSIG
     MSARR(I) = M ! 10^10 M_solar 
     SIGMARR(I) = SIG ! dimensionless 
     DSIGMARR(I) = DSIG ! (10^10 M_solar)^-1 
  END DO

  CLOSE(11) 

  !-------------------------------  

  RHO = OMEGA_NR*RHO_CRITICAL ! 10^10 M_solar / Mpc^3 
  RE = 0.673_PREC ! Mpc
  Z = 8.0_PREC
  DLTLS = 4.0_PREC
  QSOMASS = 2.52E1_PREC ! 10^10 M_solar
  DO
     IF(DLTLS > 10.0_PREC) EXIT
     ! Calculate parameter \theta (T).
     T = THETA(DLTLS, Z)
     RS = (20.0_PREC * DLTLS / (3.0_PREC*(6.0_PREC*(T - SIN(T)))&
          &**(2.0_PREC/3.0_PREC))) - 1.0_PREC
     
     GRLOC = COUNTER(Z); GROWTH = GRFARR(GRLOC)
     RL = 10.0_PREC * RE * DLTLS * GROWTH / (3.0_PREC *(1.0_PREC-COS(T))) ! Mpc
     
     !-------------------------------  

     RHOBAR = OMEGA_NR * RHO_CRITICAL ! 10^10 M_solar / Mpc^3 
     M = 4.0_PREC * PI * RL**3 * RHOBAR ! 10^10 M_solar
     MDIFF = ABS(M - MSARR); ALC = MINLOC(MDIFF)
     SGLS = SIGMARR(ALC(1))

     !-------------------------------  
     
     QLS = Q(DLTLS,SGLS)
     
     DLTSS = DELTAC / GROWTH
     M = QSOMASS
     MDIFF = ABS(M - MSARR); ALC = MINLOC(MDIFF)
     SGSS = SIGMARR(ALC(1))

     QSS = Q(DLTSS-DLTLS,SGSS-SGLS)
     PD = QSS * QLS 
     PRINT '(3E23.6)', DLTLS, PD, RL 
     
     DLTLS = DLTLS + 0.01_PREC
  END DO
  
CONTAINS 
  
  FUNCTION Q(D,S)
    
    IMPLICIT NONE
    REAL(KIND = PREC), INTENT(IN) :: D, S
    REAL(KIND = PREC) :: Q 
    
    Q = EXP(-D*D/(2.0_PREC*S))/SQRT(2.0_PREC*PI*S)
    
  END FUNCTION Q
  
  ! Fix RE
  ! Loop over a range of delta values; For each delta
  !    Calculate RL(delta,RE)
  !    Calculare sigma_largescale(RL)
  !    Calculate Q1=Q(delta,sigma_largescale)
  !    Calculate Q2=Q(delta_smallscale-delta, sigma_smallscale-sigma_largescale)
  !    PDELTA is Q1*Q2

END PROGRAM PDELTA
