
MODULE INTERFACES

  ! File: interfaces.f90
  ! Cre: 2010-04-07
  ! Mod: $Date: 2012/08/27 09:39:45 $ ($Revision: 1.22 $)
  ! 
  ! This module contains a single interface block.
  ! It is used in many subprograms of reion to 
  ! make several interfaces explicit.  

  INTERFACE 

     FUNCTION LUMINOSITY(FORMATION_REDSHIFT, TFORM, T)

       USE CONSTANTS; IMPLICIT NONE
       REAL(KIND = PREC), INTENT(IN) :: FORMATION_REDSHIFT, &
            &TFORM, T 
       REAL(KIND = PREC) :: LUMINOSITY

     END FUNCTION LUMINOSITY

     ! -----------------------

     FUNCTION NUINT(NU)

       USE CONSTANTS; IMPLICIT NONE 
       REAL(KIND = PREC), INTENT(IN) :: NU 
       REAL(KIND = PREC) :: NUINT 

     END FUNCTION NUINT

     ! -----------------------

     FUNCTION GETFF(RS)

       USE CONSTANTS; USE STORAGE; IMPLICIT NONE
       REAL(KIND = PREC), INTENT(IN) :: RS
       REAL(KIND = PREC) :: GETFF

     END FUNCTION GETFF

     ! -----------------------

     FUNCTION COUNTER(Z) 

       USE CONSTANTS; IMPLICIT NONE
       REAL(KIND = PREC) :: Z
       INTEGER :: COUNTER

     END FUNCTION COUNTER

     ! -----------------------

     FUNCTION GETJMC(RS)

       USE CONSTANTS; IMPLICIT NONE
       REAL(KIND = PREC), INTENT(IN) :: RS
       REAL(KIND = PREC) :: GETJMC

     END FUNCTION GETJMC

     ! -----------------------

     FUNCTION GETJMH(RS)

       USE CONSTANTS; USE STORAGE; IMPLICIT NONE
       REAL(KIND = PREC), INTENT(IN) :: RS
       REAL(KIND = PREC) :: GETJMH

     END FUNCTION GETJMH

     ! -----------------------

     SUBROUTINE SOURCEOV(SOURCE, REDSHIFT)

       USE CONSTANTS; USE STORAGE; IMPLICIT NONE
       REAL(KIND = PREC), INTENT(IN) :: REDSHIFT
       REAL(KIND = PREC), INTENT(OUT) :: SOURCE 

     END SUBROUTINE SOURCEOV

     ! ---------------------

     FUNCTION PIGMD(DLTA)

       USE CONSTANTS; USE STORAGE; IMPLICIT NONE
       REAL(KIND=PREC),INTENT(IN) :: DLTA
       REAL(KIND=PREC) :: PIGMD

     END FUNCTION PIGMD

     ! ----------------------

     FUNCTION EDENSITY(REDSHIFT)

       USE CONSTANTS; USE STORAGE; IMPLICIT NONE
       REAL(KIND = PREC), INTENT(IN) :: REDSHIFT
       REAL(KIND = PREC) :: EDENSITY

     END FUNCTION EDENSITY

     ! ----------------------

     FUNCTION PSPEC(K)
       USE CONSTANTS; USE STORAGE; IMPLICIT NONE 
       REAL (KIND = PREC), INTENT (IN) :: K
       REAL (KIND = PREC) :: PSPEC
     END FUNCTION PSPEC

     ! ----------------------

     FUNCTION IGMFPREO()

       USE CONSTANTS; USE STORAGE 
       IMPLICIT NONE
       REAL(KIND=PREC) :: IGMFPREO

     END FUNCTION IGMFPREO

     ! ----------------------

     FUNCTION HUBP(Z) 

       USE CONSTANTS; IMPLICIT NONE
       REAL(KIND=PREC), INTENT(IN) :: Z
       REAL(KIND=PREC) :: HUBP

     END FUNCTION HUBP

     ! ----------------------

     FUNCTION DZDT(Z) 

       USE CONSTANTS
       IMPLICIT NONE
       REAL(KIND=PREC), INTENT(IN) :: Z
       REAL(KIND=PREC) :: DZDT

     END FUNCTION DZDT

     ! ----------------------

     FUNCTION DTDZ(Z) 

       USE CONSTANTS
       IMPLICIT NONE
       REAL(KIND=PREC), INTENT(IN) :: Z
       REAL(KIND=PREC) :: DTDZ

     END FUNCTION DTDZ

     ! ----------------------

     FUNCTION SIGMUE(MU)

       USE CONSTANTS; IMPLICIT NONE
       REAL(KIND = PREC), INTENT(IN) :: MU
       REAL(KIND = PREC) :: SIGMUE

     END FUNCTION SIGMUE

     ! ----------------------

     FUNCTION NMU(NU)

       USE CONSTANTS; IMPLICIT NONE 
       REAL(KIND=PREC), INTENT(IN) :: NU
       REAL(KIND=PREC) :: NMU

     END FUNCTION NMU

     ! ----------------------

     FUNCTION SIGMUP(MU)

       USE CONSTANTS; IMPLICIT NONE
       REAL(KIND = PREC), INTENT(IN) :: MU
       REAL(KIND = PREC) :: SIGMUP

     END FUNCTION SIGMUP

     ! ----------------------

     SUBROUTINE INTERPOLATE(YARRAY, XARRAY, X, Y)

       USE CONSTANTS; IMPLICIT NONE
       REAL(KIND=PREC), DIMENSION(:), INTENT(IN) :: YARRAY, XARRAY
       REAL(KIND=PREC), INTENT(IN) :: X
       REAL(KIND=PREC), INTENT(OUT) :: Y 

     END SUBROUTINE INTERPOLATE

     SUBROUTINE INTERPOLATE2(YARRAY, XARRAY, X, Y)

       USE CONSTANTS; IMPLICIT NONE
       REAL(KIND=PREC), DIMENSION(:), INTENT(IN) :: YARRAY, XARRAY
       REAL(KIND=PREC), INTENT(IN) :: X
       REAL(KIND=PREC), INTENT(OUT) :: Y 

     END SUBROUTINE INTERPOLATE2

     ! ----------------------

     FUNCTION CLUMPFAC(DLTA)

       USE CONSTANTS; USE STORAGE 
       IMPLICIT NONE
       REAL(KIND=PREC),INTENT(IN) :: DLTA
       REAL(KIND=PREC) :: CLUMPFAC

     END FUNCTION CLUMPFAC

     ! ----------------------

     FUNCTION IGMDTH(DFM) 

       USE CONSTANTS; IMPLICIT NONE 
       REAL(KIND=PREC), INTENT(IN) :: DFM
       REAL(KIND=PREC) :: IGMDTH 

     END FUNCTION IGMDTH

     ! ----------------------

     FUNCTION SIGMA_BARYON(Z, T) 

       USE CONSTANTS; IMPLICIT NONE 
       REAL(KIND=PREC), INTENT(IN) :: Z, T 
       REAL(KIND=PREC) :: SIGMA_BARYON 

     END FUNCTION SIGMA_BARYON

     ! ----------------------

     FUNCTION SIGBINT(K)

       USE CONSTANTS; IMPLICIT NONE 
       REAL(KIND=PREC), INTENT(IN) :: K
       REAL(KIND=PREC) :: SIGBINT 

     END FUNCTION SIGBINT

     ! ----------------------

     FUNCTION JEANS_LENGTH(TMPTURE, RSHFT) 

       USE CONSTANTS; IMPLICIT NONE
       REAL(KIND=PREC), INTENT(IN) :: TMPTURE, RSHFT 
       REAL(KIND=PREC) :: JEANS_LENGTH

     END FUNCTION JEANS_LENGTH

     ! ----------------------

     FUNCTION PSIG(K)
       USE CONSTANTS; USE STORAGE; IMPLICIT NONE 
       REAL (KIND = PREC), INTENT (IN) :: K
       REAL (KIND = PREC) :: PSIG
     END FUNCTION PSIG

     ! ----------------------

     FUNCTION IGMFPOSTO(DLTA)

       USE CONSTANTS; USE STORAGE 
       IMPLICIT NONE
       REAL(KIND=PREC),INTENT(IN) :: DLTA 
       REAL(KIND=PREC) :: IGMFPOSTO

     END FUNCTION IGMFPOSTO

     ! ----------------------

     FUNCTION DLTA_THR(Z, FM, INIT)

       USE CONSTANTS; USE STORAGE; IMPLICIT NONE 
       REAL(KIND = PREC), INTENT(IN) :: Z, FM, INIT
       REAL(KIND = PREC) :: DLTA_THR 

     END FUNCTION DLTA_THR

     ! ---------------------- 

     FUNCTION SOLDLT(DLT)

       USE CONSTANTS; USE STORAGE 
       IMPLICIT NONE 

       REAL(KIND = PREC), INTENT(IN) :: DLT
       REAL(KIND = PREC) :: SOLDLT 

     END FUNCTION SOLDLT

     ! ---------------------- 

     FUNCTION SOLFM(FM)

       USE CONSTANTS; USE STORAGE 
       IMPLICIT NONE 

       REAL(KIND = PREC), INTENT(IN) :: FM
       REAL(KIND = PREC) :: SOLFM

     END FUNCTION SOLFM

     ! ---------------------- 

     FUNCTION RTBIS(FUNC, XLO, XHI, TOL, CALLER) 

       USE CONSTANTS; IMPLICIT NONE 
       REAL(KIND=PREC), INTENT(IN) :: XLO, XHI, TOL 
       CHARACTER(LEN=*), INTENT(IN) :: CALLER  
       REAL(KIND=PREC) :: RTBIS 
       INTERFACE
          FUNCTION FUNC(X)
            USE CONSTANTS; IMPLICIT NONE
            REAL(KIND=PREC), INTENT(IN) :: X
            REAL(KIND=PREC) :: FUNC
          END FUNCTION FUNC
       END INTERFACE

     END FUNCTION RTBIS

     ! ---------------------- 

     function rtnewt(funcd, first_guess, xlo, xhi, tol, caller) 

       ! find root by bisection method.
       use constants; implicit none 
       real(kind=prec), intent(in) :: first_guess, tol 
       real(kind=prec), intent(inout) :: xlo, xhi 
       character(len=*), intent(in) :: caller 
       real(kind=prec) :: rtnewt 
       interface
          subroutine funcd(x, f, df)
            use constants; implicit none
            real(kind=prec), intent(in) :: x
            real(kind=prec), intent(out) :: f, df 
          end subroutine funcd
       end interface

     end function rtnewt

     ! ---------------------- 

     FUNCTION IGMVFRAC(DLTA)

       USE CONSTANTS; USE STORAGE 
       IMPLICIT NONE
       REAL(KIND=PREC),INTENT(IN) :: DLTA 
       REAL(KIND=PREC) :: IGMVFRAC 

     END FUNCTION IGMVFRAC

     ! ---------------------- 

     FUNCTION GAMMA_PH(Z) 

       USE CONSTANTS 
       IMPLICIT NONE 

       REAL(KIND=PREC), INTENT(IN) :: Z 
       REAL(KIND=PREC) :: GAMMA_PH 

     END FUNCTION GAMMA_PH

     ! ---------------------- 

     FUNCTION SIGMAH(NU) 

       USE CONSTANTS 
       IMPLICIT NONE 
       REAL(KIND=PREC), INTENT(IN) :: NU 
       REAL(KIND=PREC) :: SIGMAH 

     END FUNCTION SIGMAH

     ! ---------------------- 

     FUNCTION GINT(NU) 

       USE CONSTANTS
       IMPLICIT NONE 
       REAL(KIND=PREC), INTENT(IN) :: NU 
       REAL(KIND=PREC) :: GINT 

     END FUNCTION GINT

     ! ---------------------- 

     FUNCTION GAMMA_PI(Z) 

       USE CONSTANTS 
       USE STORAGE 
       IMPLICIT NONE 

       REAL(KIND=PREC), INTENT(IN) :: Z 
       REAL(KIND=PREC) :: GAMMA_PI 

     END FUNCTION GAMMA_PI

     ! ---------------------- 

     FUNCTION GPINT(NU) 

       USE CONSTANTS
       IMPLICIT NONE 
       REAL(KIND=PREC), INTENT(IN) :: NU 
       REAL(KIND=PREC) :: GPINT 

     END FUNCTION GPINT

     ! ---------------------- 

     SUBROUTINE ZBRAC(FUNC, X1, X2) 

       USE CONSTANTS; IMPLICIT NONE 
       INTEGER, PARAMETER :: NTRY = 50
       REAL(KIND=PREC), INTENT(INOUT) :: X1, X2
       INTERFACE
          FUNCTION FUNC(X)
            USE CONSTANTS; IMPLICIT NONE
            REAL(KIND=PREC), INTENT(IN) :: X
            REAL(KIND=PREC) :: FUNC
          END FUNCTION FUNC
       END INTERFACE

     END SUBROUTINE ZBRAC

     ! ---------------------- 

     FUNCTION GPH_KERNEL() 

       USE CONSTANTS; USE STORAGE 
       IMPLICIT NONE 

       REAL(KIND=PREC) :: GPH_KERNEL

     END FUNCTION GPH_KERNEL

     ! ---------------------- 

     FUNCTION GPI_KERNEL() 

       USE CONSTANTS; USE STORAGE 
       IMPLICIT NONE 

       REAL(KIND=PREC) :: GPI_KERNEL

     END FUNCTION GPI_KERNEL

     ! ---------------------- 

     subroutine lumfn(region, mab, z, n) 

       use constants 
       use storage 
       implicit none 

       integer, intent(in) :: region 
       real(kind=prec), intent(in) :: mab, z 
       real(kind=prec), intent(out) :: n

     end subroutine lumfn

     ! ---------------------- 

     subroutine gallum(HaloTotalMass, HaloAge, HaloFormationRedshift, HaloLuminosity)

       use constants 
       implicit none 

       real(kind=prec), intent(in) :: HaloTotalMass, HaloAge, HaloFormationRedshift 
       real(kind=prec), intent(out) :: HaloLuminosity 

     end subroutine gallum

     ! ---------------------- 

     subroutine nfchk()

       use constants
       use storage 
       implicit none 

     end subroutine nfchk

     ! ---------------------- 

     FUNCTION NUINT_POP2(NU)

       USE CONSTANTS; USE STORAGE; IMPLICIT NONE 
       REAL(KIND = PREC), INTENT(IN) :: NU 
       REAL(KIND = PREC) :: NUINT_POP2

     END FUNCTION NUINT_POP2

     ! ---------------------- 

     FUNCTION NUINT_POP3(NU)

       USE CONSTANTS; USE STORAGE; IMPLICIT NONE 
       REAL(KIND = PREC), INTENT(IN) :: NU 
       REAL(KIND = PREC) :: NUINT_POP3

     END FUNCTION NUINT_POP3

     ! ---------------------- 

     DOUBLE PRECISION FUNCTION fpop3(zcoll,m)
       IMPLICIT NONE

       double precision, intent(in) :: zcoll,m
     END FUNCTION fpop3

     ! ---------------------

     SUBROUTINE SOURCEOV_POP3(SOURCE, REDSHIFT)

       USE CONSTANTS; USE STORAGE; IMPLICIT NONE
       REAL(KIND = PREC), INTENT(IN) :: REDSHIFT
       REAL(KIND = PREC), INTENT(OUT) :: SOURCE

     END SUBROUTINE SOURCEOV_POP3

     ! ---------------------

     SUBROUTINE SOURCEOV_POP2(SOURCE, REDSHIFT)

       USE CONSTANTS; USE STORAGE; IMPLICIT NONE
       REAL(KIND = PREC), INTENT(IN) :: REDSHIFT
       REAL(KIND = PREC), INTENT(OUT) :: SOURCE

     END SUBROUTINE SOURCEOV_POP2

     ! ---------------------- 

     FUNCTION GPH_KERNEL_POP2() 

       USE CONSTANTS; USE STORAGE 
       IMPLICIT NONE 

       REAL(KIND=PREC) :: GPH_KERNEL_POP2

     END FUNCTION GPH_KERNEL_POP2

     ! ---------------------- 

     FUNCTION GPH_KERNEL_POP3() 

       USE CONSTANTS; USE STORAGE 
       IMPLICIT NONE 

       REAL(KIND=PREC) :: GPH_KERNEL_POP3

     END FUNCTION GPH_KERNEL_POP3

     ! ---------------------- 

     FUNCTION GPI_KERNEL_POP2() 

       USE CONSTANTS; USE STORAGE 
       IMPLICIT NONE 

       REAL(KIND=PREC) :: GPI_KERNEL_POP2

     END FUNCTION GPI_KERNEL_POP2

     ! ---------------------- 

     FUNCTION GPI_KERNEL_POP3() 

       USE CONSTANTS; USE STORAGE 
       IMPLICIT NONE 

       REAL(KIND=PREC) :: GPI_KERNEL_POP3

     END FUNCTION GPI_KERNEL_POP3

     ! ---------------------- 

     function fbint(nu) 

       use constants 
       implicit none 
       real(kind=prec), intent(in) :: nu 
       real(kind=prec) :: fbint 

     end function fbint

     ! ---------------------- 

     function accrate(z) 

       use constants 
       implicit none 
       real(kind=prec), intent(in) :: z 
       real(kind=prec) :: accrate 

     end function accrate

     ! ---------------------- 

     function imf(m) 

       use constants 
       implicit none 
       real(kind=prec), intent(in) :: m 
       real(kind=prec) :: imf

     end function imf

     ! ---------------------- 

     function imf_pop3(m) 

       use constants 
       implicit none 
       real(kind=prec), intent(in) :: m 
       real(kind=prec) :: imf_pop3

     end function imf_pop3

     ! ---------------------- 

     function fescint(nu) 

       use constants 
       implicit none 
       real(kind=prec), intent(in) :: nu 
       real(kind=prec) :: fescint 

     end function fescint

     ! ---------------------- 

     function vesc_sq(z) 

       use constants 
       implicit none 
       real(kind=prec), intent(in) :: z 
       real(kind=prec) :: vesc_sq 

     end function vesc_sq

     ! ---------------------- 

     function getsfr(rs)

       use constants; 
       implicit none
       real(kind = prec), intent(in) :: rs
       real(kind = prec) :: getsfr

     end function getsfr

     ! ---------------------- 

     function getsfr2(rs)

       use constants; 
       implicit none
       real(kind = prec), intent(in) :: rs
       real(kind = prec) :: getsfr2

     end function getsfr2

     ! ---------------------- 

     function getsfr3(rs)

       use constants; 
       implicit none
       real(kind = prec), intent(in) :: rs
       real(kind = prec) :: getsfr3

     end function getsfr3

     ! ---------------------- 

     function outflow(z) 

       use constants 
       implicit none 
       real(kind=prec), intent(in) :: z 
       real(kind=prec) :: outflow 

     end function outflow

     ! ---------------------- 

     function outfrac(HaloMass, z, pop) 

       use constants 
       implicit none 
       real(kind=prec), intent(in) :: HaloMass, z 
       integer, intent(in) :: pop
       real(kind=prec) :: outfrac 

     end function outfrac

     ! ---------------------- 

     function ejrate(z, species) 

       use constants 
       implicit none 
       real(kind=prec), intent(in) :: z 
       integer, intent(in) :: species 
       real(kind=prec) :: ejrate 

     end function ejrate

     ! ---------------------- 

     function ejfrac(z, pop) 

       use constants 
       implicit none 
       real(kind=prec), intent(in) :: z 
       integer, intent(in) :: pop 
       real(kind=prec) :: ejfrac 

     end function ejfrac

     ! ---------------------- 

     function haloyield(z, pop) 

       use constants 
       implicit none 
       real(kind=prec), intent(in) :: z 
       integer, intent(in) :: pop 
       real(kind=prec) :: haloyield

     end function haloyield

     ! ---------------------- 

     function haloyield_species(z, species, pop) 

       use constants 
       implicit none 
       real(kind=prec), intent(in) :: z 
       integer, intent(in) :: species, pop
       real(kind=prec) :: haloyield_species 

     end function haloyield_species

     ! ---------------------- 

     subroutine bi_interpolate2(xarray, yarray, zarray, x, y, z)
       use constants; implicit none
       real(kind=prec), dimension(:), intent(in) :: xarray, yarray
       real(kind=prec), dimension(:,:) :: zarray
       real(kind=prec), intent(in) :: x, y
       real(kind=prec), intent(out) :: z
     end subroutine bi_interpolate2

     ! ---------------------- 

     function getmet(rs)

       use constants; use storage; implicit none
       real(kind = prec), intent(in) :: rs
       real(kind = prec) :: getmet

     end function getmet

     ! ---------------------- 

     function sfr_rollinde_pop2(z)

       use constants
       implicit none 
       real(kind = prec), intent(in) :: z 
       real(kind = prec) :: sfr_rollinde_pop2 

     end function sfr_rollinde_pop2

     ! ---------------------- 

     function sfr_rollinde_pop3(z)

       use constants
       implicit none 
       real(kind = prec), intent(in) :: z 
       real(kind = prec) :: sfr_rollinde_pop3

     end function sfr_rollinde_pop3

     ! ---------------------- 

     subroutine newtondelta(dlta, func, dfunc)
       use constants
       use storage
       implicit none
       real(kind=prec), intent(in) :: dlta
       real(kind=prec), intent(out) :: func, dfunc
     end subroutine newtondelta

     ! ---------------------- 

     subroutine newtonfm(fm, func, dfunc) 
       use constants
       use storage
       implicit none
       real(kind=prec), intent(in) :: fm
       real(kind=prec), intent(out) :: func, dfunc
     end subroutine newtonfm

     ! ---------------------- 

     subroutine calculate_nuintegral(z)
       use constants 
       implicit none
       real(kind=prec), intent(in) :: z 
     end subroutine calculate_nuintegral

     ! ---------------------- 

     function ngammafrac(z) 
       use constants 
       implicit none 
       real(kind=prec), intent(in) :: z 
       real(kind=prec) :: ngammafrac
     end function ngammafrac

     ! ---------------------- 

     function getsfr_hot(rs, halomassbin)

       use constants; use storage; implicit none
       real(kind = prec), intent(in) :: rs
       integer, intent(in) :: halomassbin
       real(kind = prec) :: getsfr_hot

     end function getsfr_hot

     ! ---------------------- 

     function getsfr_cold(rs, halomassbin)

       use constants; use storage; implicit none
       real(kind = prec), intent(in) :: rs
       integer, intent(in) :: halomassbin 
       real(kind = prec) :: getsfr_cold

     end function getsfr_cold

     ! ---------------------- 

     function getpop_hot(rs, halomassbin)

       use constants; use storage; implicit none
       real(kind = prec), intent(in) :: rs
       integer, intent(in) :: halomassbin
       integer :: getpop_hot

     end function getpop_hot

     ! ---------------------- 

     function getpop_cold(rs, halomassbin)

       use constants; use storage; implicit none
       real(kind = prec), intent(in) :: rs
       integer, intent(in) :: halomassbin
       integer :: getpop_cold

     end function getpop_cold

     ! ---------------------- 

     function ejfrac_nonira(z, hotcold, bin) 

       use constants; use storage; implicit none 
       real(kind=prec), intent(in) :: z 
       integer, intent(in) :: hotcold, bin 
       real(kind=prec) :: ejfrac_nonira

     end function ejfrac_nonira

     ! ---------------------- 

     function outfrac_nonira(HaloMass, z, hotcold, bin) 

       use constants 
       use storage 
       implicit none 
       real(kind=prec), intent(in) :: HaloMass, z ! [HaloMass] = 10^10 M_solar 
       integer, intent(in) :: hotcold, bin
       real(kind=prec) :: outfrac_nonira 

     end function outfrac_nonira

     ! ---------------------- 

     function haloyield_nonira(z, hotcold, bin) 

       use constants 
       use storage 
       implicit none 
       real(kind=prec), intent(in) :: z 
       integer, intent(in) :: hotcold, bin  
       real(kind=prec) :: haloyield_nonira

     end function haloyield_nonira

     ! ---------------------- 

     function haloyield_species_nonira(z, species, hotcold, bin) 

       use constants 
       use storage 
       implicit none 
       real(kind=prec), intent(in) :: z 
       integer, intent(in) :: species, hotcold, bin
       real(kind=prec) :: haloyield_species_nonira

     end function haloyield_species_nonira

     ! ---------------------- 

     function hallum(haloindex, redshiftindex, hotcold) 

       use constants 
       use storage 
       implicit none 
       integer, intent(in) :: haloindex, redshiftindex, hotcold
       real(kind=prec) :: hallum 

     end function hallum

     ! ---------------------- 

     function hallum2(haloindex, redshiftindex, hotcold) 

       use constants 
       use storage 
       implicit none 
       integer, intent(in) :: haloindex, redshiftindex, hotcold
       real(kind=prec) :: hallum2

     end function hallum2

  END INTERFACE

END MODULE INTERFACES

