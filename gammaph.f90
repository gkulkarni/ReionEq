
  ! File: gammaph.f90
  ! Cre: 2010-04-07
  ! Mod: $Date: 2012/09/17 07:25:40 $; $Revision: 1.11 $ 

  ! Calculates several time-independent quantities relevant to
  ! gamma_pi and gamma_ph calculation in reion.

FUNCTION GPH_KERNEL_POP2() 

  USE CONSTANTS; USE STORAGE 
  USE INTERFACES, ONLY : SIGMAH 
  IMPLICIT NONE 

  REAL(KIND=PREC) :: GPH_KERNEL_POP2
  REAL(KIND=PREC) :: SUM, INT1, INT2 
  INTEGER :: I 

  SUM = 0.0_PREC 
  INT1 = IONNNUM(1)*SIGMAH(IONFREQ(1))*&
       &HPLANCK*(IONFREQ(1)-NU0)
  DO I = 2, SIZE(IONFREQ) 
     IF (IONFREQ(I)<NU0) EXIT 
     INT2 = IONNNUM(I)*SIGMAH(IONFREQ(I))*&
          &HPLANCK*(IONFREQ(I)-NU0) 
     SUM = SUM + 0.5*(IONFREQ(I-1)-IONFREQ(I))*(INT1+INT2)
     INT2 = INT1 
  END DO

  GPH_KERNEL_POP2 = SUM 

END FUNCTION GPH_KERNEL_POP2

!------------------------------

FUNCTION GPI_KERNEL_POP2() 

  USE CONSTANTS; USE STORAGE 
  USE INTERFACES, ONLY : SIGMAH 
  IMPLICIT NONE 

  REAL(KIND=PREC) :: GPI_KERNEL_POP2
  REAL(KIND=PREC) :: SUM, INT1, INT2 
  INTEGER :: I 

  SUM = 0.0_PREC
  INT1 = IONNNUM(1)*SIGMAH(IONFREQ(1))
  DO I = 2, SIZE(IONFREQ) 
     IF (IONFREQ(I)<NU0) EXIT 
     INT2 = IONNNUM(I)*SIGMAH(IONFREQ(I))
     SUM = SUM + 0.5*(IONFREQ(I-1)-IONFREQ(I))*(INT1+INT2)
     INT2 = INT1 
  END DO

  GPI_KERNEL_POP2 = SUM 

END FUNCTION GPI_KERNEL_POP2

!------------------------------

FUNCTION GPH_KERNEL_POP3() 

  USE CONSTANTS; USE STORAGE 
  USE INTERFACES, ONLY : SIGMAH 
  IMPLICIT NONE 

  REAL(KIND=PREC) :: GPH_KERNEL_POP3
  REAL(KIND=PREC) :: SUM, INT1, INT2 
  INTEGER :: I 

  SUM = 0.0_PREC 
  INT1 = IONNUMPOP3(1)*SIGMAH(IONFREQPOP3(1))*&
       &HPLANCK*(IONFREQPOP3(1)-NU0)
  DO I = 2, SIZE(IONFREQPOP3) 
     IF (IONFREQPOP3(I)<NU0) EXIT 
     INT2 = IONNUMPOP3(I)*SIGMAH(IONFREQPOP3(I))*&
          &HPLANCK*(IONFREQPOP3(I)-NU0) 
     SUM = SUM + 0.5*(IONFREQPOP3(I-1)-IONFREQPOP3(I))*(INT1+INT2)
     INT2 = INT1 
  END DO

  GPH_KERNEL_POP3 = SUM 

END FUNCTION GPH_KERNEL_POP3

!------------------------------

FUNCTION GPI_KERNEL_POP3() 

  USE CONSTANTS; USE STORAGE 
  USE INTERFACES, ONLY : SIGMAH 
  IMPLICIT NONE 

  REAL(KIND=PREC) :: GPI_KERNEL_POP3
  REAL(KIND=PREC) :: SUM, INT1, INT2 
  INTEGER :: I 

  SUM = 0.0_PREC
  INT1 = IONNUMPOP3(1)*SIGMAH(IONFREQPOP3(1))
  DO I = 2, SIZE(IONFREQPOP3) 
     IF (IONFREQPOP3(I)<NU0) EXIT 
     INT2 = IONNUMPOP3(I)*SIGMAH(IONFREQPOP3(I))
     SUM = SUM + 0.5*(IONFREQPOP3(I-1)-IONFREQPOP3(I))*(INT1+INT2)
     INT2 = INT1 
  END DO

  GPI_KERNEL_POP3 = SUM 

END FUNCTION GPI_KERNEL_POP3

!------------------------------

FUNCTION SIGMAH(NU) 

  USE CONSTANTS 
  IMPLICIT NONE 
  REAL(KIND=PREC), INTENT(IN) :: NU 
  REAL(KIND=PREC) :: SIGMAH 

  REAL(KIND=PREC) :: A0, EPS  

  ! Osterbrock book 
  A0 = 6.3E-18_PREC ! cm^2 

  IF (NU<NU0) THEN 
     SIGMAH = 0.0_PREC 
  ELSE IF ((NU/NU0-1.0_PREC)==0.0_PREC) THEN
     SIGMAH = A0 * (NU0/NU)**4 
  ELSE 
     EPS = SQRT(NU/NU0-1.0_PREC)
     SIGMAH = A0 * (NU0/NU)**4 * EXP(4.0_PREC - 4.0_PREC*ATAN(EPS)/EPS)&
          &/(1.0_PREC-EXP(-2.0_PREC*PI/EPS))
  END IF

END FUNCTION SIGMAH

