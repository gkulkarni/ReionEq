
! File: pigmd.f90
! Cre: 2010-04-07
! Mod: $Date: 2012/09/17 07:25:40 $; $Revision: 1.14 $ 

MODULE PIGMD_PRIVATE 

  USE CONSTANTS; IMPLICIT NONE 
  REAL(KIND=PREC) :: TEMP, BETA, RSHFT 

END MODULE PIGMD_PRIVATE

FUNCTION SIGMA_BARYON(Z, T) 

  USE CONSTANTS; USE INTERFACES, ONLY : SIGBINT, COUNTER
  USE STORAGE; USE PIGMD_PRIVATE; IMPLICIT NONE 
  REAL(KIND=PREC), INTENT(IN) :: Z, T 
  REAL(KIND=PREC) :: SIGMA_BARYON 

  REAL(KIND=PREC) :: BOUND, SUM, ERROR, GRFAC 
  INTEGER :: INF, NEVAL, IER, LAST

  BOUND = 0.0_PREC
  INF = 1 

  TEMP = T
  RSHFT = Z 
  IER = -1 
  CALL DQAGIE(SIGBINT,BOUND,INF,ABSERR,ABSREL,MAXINT,SUM,ERROR,&
       &NEVAL,IER,ALIST,BLIST,RLIST,ELIST,IORD,LAST)
  IF (IER > 0) WRITE (0,*) 'SIGMA_BARYON: Error', IER

  GRFAC = GRFARR(COUNTER(Z))
  SIGMA_BARYON = SQRT(PSPECNORM*SUM)*GRFAC ! dimensionless 

END FUNCTION SIGMA_BARYON

FUNCTION SIGBINT(K)

  USE CONSTANTS; USE STORAGE; USE PIGMD_PRIVATE; 
  USE INTERFACES, ONLY : JEANS_LENGTH, PSPEC
  IMPLICIT NONE 
  REAL(KIND=PREC), INTENT(IN) :: K
  REAL(KIND=PREC) :: SIGBINT 

  REAL(KIND=PREC) :: XB 

  XB = JEANS_LENGTH(TEMP, RSHFT)
  SIGBINT = (PSPEC(K)/((1.0_PREC+(XB*K)**2)**2))&
       &*K**2/(2.0_PREC*PI*PI) ! dimensionless 

END FUNCTION SIGBINT

FUNCTION PSPEC(K)

  USE CONSTANTS; USE STORAGE; IMPLICIT NONE 
  REAL (KIND = PREC), INTENT (IN) :: K
  REAL (KIND = PREC) :: PSPEC

  REAL (KIND = PREC) :: A, B, C, NU, LGAMMA

  ! Calculate the EBW power spectrum P(k).
  LGAMMA = OMEGA_NR*SMALLH
  NU = 1.13_PREC
  A = 6.4_PREC/LGAMMA
  B = 3.0_PREC/LGAMMA
  C = 1.7_PREC/LGAMMA

  PSPEC = (K**SPECTRAL_INDEX)/(1.0_PREC+(A*K+(B*K)**1.5_PREC&
       & + (C*K)**2)**NU)**(2.0_PREC/NU)

END FUNCTION PSPEC

FUNCTION JEANS_LENGTH(TMPTURE, RSHFT) 

  USE CONSTANTS 
  USE STORAGE 
  USE INTERFACES, ONLY : HUBP 
  IMPLICIT NONE
  REAL(KIND=PREC), INTENT(IN) :: TMPTURE, RSHFT 
  REAL(KIND=PREC) :: JEANS_LENGTH

  JEANS_LENGTH = SQRT(2.0_PREC*KBOLTZ*TMPTURE*ABATIC_GAMMA / &
       &(8.0_PREC*PI*NEWTG*MPROTON*RHO_CRITICAL*OMEGA_NR*&
       &(1.0_PREC+RSHFT)*(1.0E10_PREC*MSOLKG)&
       &*(CMBYMPC*100.0_PREC))) ! Mpc 

!!$  JEANS_LENGTH = SQRT(2.0_PREC*KBOLTZ*TMPTURE*ABATIC_GAMMA / &
!!$       &(8.0_PREC*PI*NEWTG*mu_MeanMolWt*MPROTON*RHO_CRITICAL*OMEGA_NR*&
!!$       &(1.0_PREC+RSHFT)*(1.0E10_PREC*MSOLKG)&
!!$       &*(CMBYMPC*100.0_PREC))) ! Mpc 

END FUNCTION JEANS_LENGTH

FUNCTION PSIG(K)

  USE CONSTANTS; USE STORAGE
  USE INTERFACES, ONLY : PSPEC
  IMPLICIT NONE 
  REAL(KIND = PREC), INTENT(IN) :: K
  REAL(KIND = PREC) :: PSIG 

  REAL(KIND = PREC) :: R 

  R = RCARRY_PSIG ! Mpc 
  ! PSIG is the integrand for calculating sigma.
  ! This is for normalization.
  PSIG = 9.0_PREC*K*K*PSPEC(K)*(((SIN(K*R)-K*R*COS(K*R))/((K*R)**3))&
       &**2)/(2.0_PREC*PI*PI)

END FUNCTION PSIG

FUNCTION IGMVFRAC(DLTA)

  USE CONSTANTS; USE STORAGE 
  USE PIGMD_PRIVATE
  USE INTERFACES, ONLY : COUNTER
  IMPLICIT NONE
  REAL(KIND=PREC),INTENT(IN) :: DLTA
  REAL(KIND=PREC) :: IGMVFRAC 

  REAL(KIND=PREC) :: MU, DLTA_V, DOFZ

  ! Calculate value of density at which a transition from lognormal 
  ! to power law PDF occurs; DLTA_V.  See Sec. 3 of CF05.  DLTA_V 
  ! is dimensionless. 
  BETA = -2.3_PREC 
  DLTA_V = (0.5_PREC*(1.0_PREC&
       &-DERF(SIGMAB*(BETA+1.0_PREC)/SQRT(2.0_PREC)))&
       &-EXP(-SIGMAB**2*(BETA+1.0_PREC)**2/2.0_PREC)&
       &/(SIGMAB*(BETA+1.0_PREC)*SQRT(2.0_PREC*PI)))& ! Numerator 
       &/(0.5_PREC*EXP(SIGMAB**2*(BETA+1.5_PREC))&
       &*(1.0_PREC-DERF(SIGMAB*(BETA+2.0_PREC)/SQRT(2.0_PREC)))&
       &-EXP(-SIGMAB**2*(BETA+1.0_PREC)**2/2.0_PREC)&
       &/(SIGMAB*(BETA+2.0_PREC)*SQRT(2.0_PREC*PI))) ! Denominator 
  DOFZ = GRFARR(COUNTER(RSHFT))
  DLTA_V = DLTA_V*(1.0_PREC+RGNOVD*DOFZ) ! For overdense regions.

  ! Calculate parameter mu.
  MU = LOG(DLTA_V)+SIGMAB**2*(BETA+1.0_PREC)

  ! Normalize P(Delta); watch the reciprocal!  
  PIGMD_NORM = 0.5_PREC&
       &+0.5_PREC*DERF((LOG(DLTA_V)-MU)/(SIGMAB*SQRT(2.0_PREC)))&
       &-(EXP(-(LOG(DLTA_V)-MU)**2/(2.0_PREC*SIGMAB**2))&
       &/(SIGMAB*(BETA+1.0_PREC)*SQRT(2.0_PREC*PI)))
  PIGMD_NORM = 1.0_PREC/PIGMD_NORM 

  ! Calculate IGMVFRAC 
  IF (DLTA < DLTA_V) THEN 
     IGMVFRAC = 0.5_PREC*&
          &(1.0_PREC+DERF(LOG(DLTA)-MU)/(SIGMAB*SQRT(2.0_PREC)))
  ELSE
     IGMVFRAC = (0.5_PREC*&
          &(1.0_PREC+DERF((LOG(DLTA_V)-MU)/(SIGMAB*SQRT(2.0_PREC)))))&
          &+(EXP(-(LOG(DLTA_V)-MU)**2/(2.0_PREC*SIGMAB**2))&
          &*((DLTA/DLTA_V)**(BETA+1.0_PREC)-1.0_PREC)&
          &/(SIGMAB*(BETA+1.0_PREC)*SQRT(2.0_PREC*PI)))
  END IF
  IGMVFRAC = PIGMD_NORM*IGMVFRAC 

END FUNCTION IGMVFRAC

FUNCTION IGMFPREO()

  USE CONSTANTS; USE STORAGE 
  USE PIGMD_PRIVATE 
  USE INTERFACES, ONLY : COUNTER 
  IMPLICIT NONE
  REAL(KIND=PREC) :: IGMFPREO

  REAL(KIND=PREC) :: MU, DLTA_V, DOFZ 

  ! Calculate value of density at which a transition from lognormal 
  ! to power law PDF occurs; DLTA_V.  See Sec. 3 of CF05.  DLTA_V 
  ! is dimensionless. 
  BETA = -2.3_PREC 
  DLTA_V = (0.5_PREC*(1.0_PREC&
       &-DERF(SIGMAB*(BETA+1.0_PREC)/SQRT(2.0_PREC)))&
       &-EXP(-SIGMAB**2*(BETA+1.0_PREC)**2/2.0_PREC)&
       &/(SIGMAB*(BETA+1.0_PREC)*SQRT(2.0_PREC*PI)))& ! Numerator 
       &/(0.5_PREC*EXP(SIGMAB**2*(BETA+1.5_PREC))&
       &*(1.0_PREC-DERF(SIGMAB*(BETA+2.0_PREC)/SQRT(2.0_PREC)))&
       &-EXP(-SIGMAB**2*(BETA+1.0_PREC)**2/2.0_PREC)&
       &/(SIGMAB*(BETA+2.0_PREC)*SQRT(2.0_PREC*PI))) ! Denominator 
  DOFZ = GRFARR(COUNTER(RSHFT))
  DLTA_V = DLTA_V*(1.0_PREC+RGNOVD*DOFZ) ! For overdense regions.

  ! Calculate parameter mu.
  MU = LOG(DLTA_V)+SIGMAB**2*(BETA+1.0_PREC)

  ! Normalize P(Delta); watch the reciprocal!  
  PIGMD_NORM = 0.5_PREC&
       &+0.5_PREC*DERF((LOG(DLTA_V)-MU)/(SIGMAB*SQRT(2.0_PREC)))&
       &-(EXP(-(LOG(DLTA_V)-MU)**2/(2.0_PREC*SIGMAB**2))&
       &/(SIGMAB*(BETA+1.0_PREC)*SQRT(2.0_PREC*PI)))
  PIGMD_NORM = 1.0_PREC/PIGMD_NORM 

  ! Return value of IGMFPREO. 
  IF (IGMDCRIT_PREO < DLTA_V) THEN 
     IGMFPREO = 0.5_PREC*EXP(MU+SIGMAB**2/2.0_PREC)*(1.0_PREC+&
          &DERF((LOG(IGMDCRIT_PREO)-MU-SIGMAB**2)/&
          &(SIGMAB*SQRT(2.0_PREC))))
  ELSE
     IGMFPREO = 0.5_PREC*EXP(MU+SIGMAB**2/2.0_PREC)*(1.0_PREC+&
          &DERF((LOG(DLTA_V)-MU-SIGMAB**2)/&
          &(SIGMAB*SQRT(2.0_PREC)))) + &
          &EXP(-(LOG(DLTA_V)-MU)**2/(2.0_PREC*SIGMAB**2))*&
          &(IGMDCRIT_PREO**(BETA+2.0_PREC)-DLTA_V**(BETA+2.0_PREC))/&
          &(SIGMAB*(BETA+2.0_PREC)*DLTA_V**(BETA+1.0_PREC)&
          &*SQRT(2.0_PREC*PI))
  END IF

  IGMFPREO = PIGMD_NORM*IGMFPREO

END FUNCTION IGMFPREO

FUNCTION CLUMPFAC(DLTA)

  USE CONSTANTS; USE STORAGE 
  USE PIGMD_PRIVATE
  USE INTERFACES, ONLY : COUNTER 
  IMPLICIT NONE
  REAL(KIND=PREC),INTENT(IN) :: DLTA
  REAL(KIND=PREC) :: CLUMPFAC

  REAL(KIND=PREC) :: MU, DLTA_V, DOFZ 

  ! Calculate value of density at which a transition from lognormal 
  ! to power law PDF occurs; DLTA_V.  See Sec. 3 of CF05.  DLTA_V 
  ! is dimensionless. 
  BETA = -2.3_PREC 
  DLTA_V = (0.5_PREC*(1.0_PREC&
       &-DERF(SIGMAB*(BETA+1.0_PREC)/SQRT(2.0_PREC)))&
       &-EXP(-SIGMAB**2*(BETA+1.0_PREC)**2/2.0_PREC)&
       &/(SIGMAB*(BETA+1.0_PREC)*SQRT(2.0_PREC*PI)))& ! Numerator 
       &/(0.5_PREC*EXP(SIGMAB**2*(BETA+1.5_PREC))&
       &*(1.0_PREC-DERF(SIGMAB*(BETA+2.0_PREC)/SQRT(2.0_PREC)))&
       &-EXP(-SIGMAB**2*(BETA+1.0_PREC)**2/2.0_PREC)&
       &/(SIGMAB*(BETA+2.0_PREC)*SQRT(2.0_PREC*PI))) ! Denominator 
  DOFZ = GRFARR(COUNTER(RSHFT))
  DLTA_V = DLTA_V*(1.0_PREC+RGNOVD*DOFZ) ! For overdense regions.

  ! Calculate parameter mu.
  MU = LOG(DLTA_V)+SIGMAB**2*(BETA+1.0_PREC)

  ! Normalize P(Delta); watch the reciprocal!  
  PIGMD_NORM = 0.5_PREC&
       &+0.5_PREC*DERF((LOG(DLTA_V)-MU)/(SIGMAB*SQRT(2.0_PREC)))&
       &-(EXP(-(LOG(DLTA_V)-MU)**2/(2.0_PREC*SIGMAB**2))&
       &/(SIGMAB*(BETA+1.0_PREC)*SQRT(2.0_PREC*PI)))
  PIGMD_NORM = 1.0_PREC/PIGMD_NORM 

  ! Evaluate clumping factor. 
  IF (DLTA < DLTA_V) THEN 
     CLUMPFAC = 0.5_PREC*EXP(2.0_PREC*(MU+SIGMAB**2))*&
          &(1.0_PREC+DERF((LOG(DLTA)-MU-2.0_PREC*SIGMAB**2)/&
          &(SIGMAB*SQRT(2.0_PREC))))
  ELSE
     CLUMPFAC = 0.5_PREC*EXP(2.0_PREC*(MU+SIGMAB**2))*&
          &(1.0_PREC+DERF((LOG(DLTA_V)-MU-2.0_PREC*SIGMAB**2)/&
          &(SIGMAB*SQRT(2.0_PREC)))) + &
          &EXP(-(LOG(DLTA_V)-MU)**2/(2.0_PREC*SIGMAB**2))*&
          &(DLTA**(BETA+3.0_PREC)-DLTA_V**(BETA+3.0_PREC))/&
          &(SIGMAB*SQRT(2.0_PREC*PI)*(BETA+3.0_PREC)*&
          &DLTA_V**(BETA+1.0_PREC))
  END IF
  CLUMPFAC = CLUMPFAC*PIGMD_NORM

END FUNCTION CLUMPFAC

FUNCTION SOLDLT(DLT)

  USE CONSTANTS; USE STORAGE 
  USE INTERFACES, ONLY : IGMFPOSTO  
  IMPLICIT NONE 

  REAL(KIND = PREC), INTENT(IN) :: DLT
  REAL(KIND = PREC) :: SOLDLT 

  SOLDLT = SOLDLT_IGMFRAC-IGMFPOSTO(DLT)

END FUNCTION SOLDLT

FUNCTION IGMFPOSTO(DLTA)

  USE CONSTANTS; USE STORAGE 
  USE PIGMD_PRIVATE 
  USE INTERFACES, ONLY : COUNTER 
  IMPLICIT NONE
  REAL(KIND=PREC),INTENT(IN) :: DLTA 
  REAL(KIND=PREC) :: IGMFPOSTO

  REAL(KIND=PREC) :: MU, DLTA_V, DOFZ  

  ! Calculate value of density at which a transition from lognormal 
  ! to power law PDF occurs; DLTA_V.  See Sec. 3 of CF05.  DLTA_V 
  ! is dimensionless. 
  BETA = -2.3_PREC 
  DLTA_V = (0.5_PREC*(1.0_PREC&
       &-DERF(SIGMAB*(BETA+1.0_PREC)/SQRT(2.0_PREC)))&
       &-EXP(-SIGMAB**2*(BETA+1.0_PREC)**2/2.0_PREC)&
       &/(SIGMAB*(BETA+1.0_PREC)*SQRT(2.0_PREC*PI)))& ! Numerator 
       &/(0.5_PREC*EXP(SIGMAB**2*(BETA+1.5_PREC))&
       &*(1.0_PREC-DERF(SIGMAB*(BETA+2.0_PREC)/SQRT(2.0_PREC)))&
       &-EXP(-SIGMAB**2*(BETA+1.0_PREC)**2/2.0_PREC)&
       &/(SIGMAB*(BETA+2.0_PREC)*SQRT(2.0_PREC*PI))) ! Denominator 
  DOFZ = GRFARR(COUNTER(RSHFT))
  DLTA_V = DLTA_V*(1.0_PREC+RGNOVD*DOFZ) ! For overdense regions.

  ! Calculate parameter mu.
  MU = LOG(DLTA_V)+SIGMAB**2*(BETA+1.0_PREC)

  ! Normalize P(Delta); watch the reciprocal!  
  PIGMD_NORM = 0.5_PREC&
       &+0.5_PREC*DERF((LOG(DLTA_V)-MU)/(SIGMAB*SQRT(2.0_PREC)))&
       &-(EXP(-(LOG(DLTA_V)-MU)**2/(2.0_PREC*SIGMAB**2))&
       &/(SIGMAB*(BETA+1.0_PREC)*SQRT(2.0_PREC*PI)))
  PIGMD_NORM = 1.0_PREC/PIGMD_NORM 

  ! Return value of IGMFPOSTO. 
  IF (DLTA < DLTA_V) THEN 
     IGMFPOSTO = 0.5_PREC*EXP(MU+SIGMAB**2/2.0_PREC)*(1.0_PREC+&
          &DERF((LOG(DLTA)-MU-SIGMAB**2)/&
          &(SIGMAB*SQRT(2.0_PREC))))
  ELSE
     IGMFPOSTO = 0.5_PREC*EXP(MU+SIGMAB**2/2.0_PREC)*(1.0_PREC+&
          &DERF((LOG(DLTA_V)-MU-SIGMAB**2)/&
          &(SIGMAB*SQRT(2.0_PREC)))) + &
          &EXP(-(LOG(DLTA_V)-MU)**2/(2.0_PREC*SIGMAB**2))*&
          &(DLTA**(BETA+2.0_PREC)-DLTA_V**(BETA+2.0_PREC))/&
          &(SIGMAB*(BETA+2.0_PREC)*DLTA_V**(BETA+1.0_PREC)&
          &*SQRT(2.0_PREC*PI))
  END IF

  IGMFPOSTO = PIGMD_NORM*IGMFPOSTO

END FUNCTION IGMFPOSTO

FUNCTION SOLFM(FM)

  USE CONSTANTS; USE STORAGE 
  USE INTERFACES, ONLY : SOLDLT, CLUMPFAC, DTDZ, HUBP, RTBIS, IGMVFRAC, &
       &ZBRAC, COUNTER
  IMPLICIT NONE 

  REAL(KIND = PREC), INTENT(IN) :: FM
  REAL(KIND = PREC) :: SOLFM

  REAL(KIND = PREC) :: Z, IGMFRAC, DCRIT, DLTHI, DLTLO, TOLZIN, DLT, &
       &R, OLDFF, FILLFACTOR, RHO, RHO_BARYON, NH, T1, T2, DOFZ 

  Z = SOLFM_Z 
  IGMFRAC = FM 
  DCRIT = SOLFM_IGMDCRIT 

  ! Calculate delta corresponding to F_M.
  SOLDLT_Z = SOLFM_Z
  SOLDLT_IGMFRAC = FM
  DLTHI = 1.0e6_prec  !DCRIT+1.0E30_PREC
  DLTLO = 0.0_prec !DCRIT-59.0_PREC
  TOLZIN = 1.0E-8_PREC
  CALL ZBRAC(SOLDLT,DLTHI,DLTLO) 
  DLT = RTBIS(SOLDLT, DLTLO, DLTHI, TOLZIN, 'solfm-dlt')

  ! "Calculate" f_II. 
  R = CLUMPFAC(DLT) 
  OLDFF = SOLFM_FF 
  FILLFACTOR = 1.0_PREC

  RHO = OMEGA_NR*RHO_CRITICAL ! 10^10 M_solar / Mpc^3 
  RHO_BARYON = OMEGA_B*RHO_CRITICAL ! 10^10 M_solar / Mpc^3 
  DOFZ = GRFARR(COUNTER(Z))
  NH = RHO_BARYON*(1.0_PREC+RGNOVD*DOFZ)*(1.1891E57_PREC*1.0E10_PREC) ! Mpc^-3 
!!$  NH = RHO_BARYON*(1.0_PREC+RGNOVD)*(1.1891E57_PREC*1.0E10_PREC) ! Mpc^-3 

!!$  T1 = OLDFF*(FM-OLDIGMFRAC)/DZ
!!$  T2 = DTDZ(Z)*NPHDOT/NH - DTDZ(Z)*ALPHA_R*YRBYS*CMBYMPCCB*&
!!$       &NH*X_II*R*(1.0_PREC+Z)**3*FILLFACTOR - FM*(FILLFACTOR-OLDFF)/DZ 
!!$  SOLFM = T1-T2 

  T1 = (FM-OLDIGMFRAC)/DZ
  T2 = DTDZ(Z)*NPHDOT/NH - DTDZ(Z)*ALPHA_R*YRBYS*CMBYMPCCB*&
       &NH*X_II*R*(1.0_PREC+Z)**3*FILLFACTOR 
  SOLFM = T1-T2 

END FUNCTION SOLFM

!!$subroutine funcd(x, f, df)
!!$  use constants; implicit none
!!$  real(kind=prec), intent(in) :: x
!!$  real(kind=prec), intent(out) :: f, df 
!!$
!!$  f = 0.0_prec
!!$  df = 0.0_prec 
!!$
!!$end subroutine funcd
!!$
!!$subroutine testf(x, f, df)
!!$  use constants; implicit none
!!$  real(kind=prec), intent(in) :: x
!!$  real(kind=prec), intent(out) :: f, df 
!!$
!!$  f = exp(x)-1.0_prec
!!$  df = exp(x)
!!$
!!$end subroutine testf
!!$
!!$function testf2(x)
!!$  use constants; implicit none
!!$  real(kind=prec), intent(in) :: x
!!$  real(kind=prec) :: testf2
!!$
!!$  testf2 = exp(x)-1.0_prec
!!$
!!$end function testf2

subroutine newtonfm(fm, func, dfunc) 
  use constants
  use storage
  use interfaces, only: dtdz, rtnewt, clumpfac, newtondelta, counter, soldlt, rtbis 
  implicit none
  real(kind=prec), intent(in) :: fm
  real(kind=prec), intent(out) :: func, dfunc
  real(kind=prec) :: z, delta_local, nh, r, rho_baryon, dlta_lo, dlta_high, tolzin

  z = solfm_z 
  nh = rho_baryon*(1.1891e57_prec*1.0e10_prec) ! Mpc^-3 

  soldlt_igmfrac = fm 
  dlta_lo = 0.0_prec
  dlta_high = 1.0e6_prec 
  delta_local = rtnewt(newtondelta, solfm_igmdcrit, dlta_lo, &
       &dlta_high, 1.0e-8_prec, 'delta_in_newtonfm')

!!$  soldlt_z = solfm_z
!!$  soldlt_igmfrac = fm 
!!$  dlta_high = 1.0e6_prec !igmdcrit+1.0e30_prec
!!$  dlta_lo = 0.0_prec !igmdcrit-59.0_prec
!!$  tolzin = 1.0e-4_prec
!!$  delta_local = rtbis(soldlt, dlta_lo, dlta_high, tolzin, 'delta_in_newtonfm')

  r = clumpfac(delta_local)

  func = (fm-oldigmfrac)-dz*dtdz(z)*(nphdot/nh - nh*x_ii*r*alpha_r*cmbympccb*yrbys*(1.0_prec+z)**3)
  dfunc = 1.0_prec + dz*dtdz(z)*(nh*x_ii*alpha_r*cmbympccb*yrbys*delta_local*(1.0_prec+z)**3) 

end subroutine newtonfm

subroutine newtondelta(dlta, func, dfunc)
  use constants
  use storage
  use interfaces, only: igmfposto, pigmd
  implicit none
  real(kind=prec), intent(in) :: dlta
  real(kind=prec), intent(out) :: func, dfunc

  func = soldlt_igmfrac - igmfposto(dlta)
  dfunc = dlta*pigmd(dlta)

end subroutine newtondelta

function pigmd(dlta)

  use constants; use storage
  use pigmd_private
  use interfaces, only: counter
  implicit none
  real(kind=prec),intent(in) :: dlta
  real(kind=prec) :: pigmd

  REAL(KIND=PREC) :: MU, DLTA_V, DOFZ

  ! Calculate value of density at which a transition from lognormal 
  ! to power law PDF occurs; DLTA_V.  See Sec. 3 of CF05.  DLTA_V 
  ! is dimensionless. 
  BETA = -2.3_PREC 
  DLTA_V = (0.5_PREC*(1.0_PREC&
       &-DERF(SIGMAB*(BETA+1.0_PREC)/SQRT(2.0_PREC)))&
       &-EXP(-SIGMAB**2*(BETA+1.0_PREC)**2/2.0_PREC)&
       &/(SIGMAB*(BETA+1.0_PREC)*SQRT(2.0_PREC*PI)))& ! Numerator 
       &/(0.5_PREC*EXP(SIGMAB**2*(BETA+1.5_PREC))&
       &*(1.0_PREC-DERF(SIGMAB*(BETA+2.0_PREC)/SQRT(2.0_PREC)))&
       &-EXP(-SIGMAB**2*(BETA+1.0_PREC)**2/2.0_PREC)&
       &/(SIGMAB*(BETA+2.0_PREC)*SQRT(2.0_PREC*PI))) ! Denominator 
!!$  DOFZ = GRFARR(COUNTER(RSHFT))
!!$  DLTA_V = DLTA_V*(1.0_PREC+RGNOVD*DOFZ) ! For overdense regions.

  ! Calculate parameter mu.
  MU = LOG(DLTA_V)+SIGMAB**2*(BETA+1.0_PREC)

  ! Normalize P(Delta); watch the reciprocal!  
  PIGMD_NORM = 0.5_PREC&
       &+0.5_PREC*DERF((LOG(DLTA_V)-MU)/(SIGMAB*SQRT(2.0_PREC)))&
       &-(EXP(-(LOG(DLTA_V)-MU)**2/(2.0_PREC*SIGMAB**2))&
       &/(SIGMAB*(BETA+1.0_PREC)*SQRT(2.0_PREC*PI)))
  PIGMD_NORM = 1.0_PREC/PIGMD_NORM 

  ! Calculate PIGMD
  IF (DLTA < DLTA_V) THEN 
     PIGMD = EXP(-(LOG(DLTA)-MU)**2/(2.0_PREC*SIGMAB**2))&
          &/(SIGMAB*DLTA*SQRT(2.0_PREC*PI))
  ELSE
     PIGMD = (EXP(-(LOG(DLTA_V)-MU)**2/(2.0_PREC*SIGMAB**2))&
          &/(SIGMAB*DLTA_V*SQRT(2.0_PREC*PI)))&
          &*(DLTA/DLTA_V)**BETA
  END IF
  PIGMD = PIGMD_NORM*PIGMD

end function pigmd
