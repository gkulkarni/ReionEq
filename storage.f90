
MODULE STORAGE

  ! File: storage.f90
  ! Cre: 2010-04-07
  ! Mod: $Date: 2012/10/23 08:10:11 $; $Revision: 1.26 $
  !
  ! Contains several global variables of reion code.

  USE CONSTANTS; IMPLICIT NONE
  REAL(KIND = PREC) :: TFORM, GRFN, DELTA, PIGMD_NORM, mu_MeanMolWt   
  REAL(KIND = PREC), DIMENSION(:), ALLOCATABLE :: SIGMARR, &
       &DSIGMARR, MSARR, GRFARR, DGRFARR, TARR, ZARR, JMHARR,&
       &FFARR, RSARR, IONNNUM, IONFREQ, MDIFF, SFRARR,&
       &SIGDIFF, IONFREQPOP3, IONNUMPOP3, STELLAR_MASS,&
       &STELLAR_AGE, metarr, febyharr, zfearr, sfrarr_pop2,&
       &sfrarr_pop3, nuintegralarr_pop2, mgas_halosc, &
       &nuintegralarr_pop3, m_halosc, mstar_halosc, &
       &mmetal_halosc, mFe_halosc, mC_halosc, mO_halosc,&
       &mN_halosc, mN_halosh, mSi_halosc, mSi_halosh, mZn_halosc, mZn_halosh,&
       &mMg_halosc, mMg_halosh, mcoolgas_halosc, mcoolgas_halosh,& 
       &febyh_halosc, cbyh_halosc, obyh_halosc, mstardot_halosc, &
       &m_halosh, mstar_halosh, mmetal_halosh, mFe_halosh, &
       &mC_halosh, mO_halosh, febyh_halosh, cbyh_halosh, obyh_halosh, &
       &mstardot_halosh, mgas_halosh, pop3_stellar_mass, &
       &pop3_stellar_age, pop3_stellar_ngamma, ismz_halosc, &
       &ismz_halosh, aux_halosc, aux_halosh, nofmc, nofmh, &
       &t_zn, t_si, t_fe, t_o, sd93_tvir, sd93_coolrate, halol1500c, &
       &halol1500h, sfrcontrib_halosc, sfrcontrib_halosh, sd93_temp, sd93_febyh
  integer, dimension(:), allocatable :: strpop_halosc, strpop_halosh 
  INTEGER, DIMENSION(1) :: ALOC, ALC
  REAL(KIND = PREC), DIMENSION(MAXINT) :: ALIST, BLIST, ELIST, &
       &RLIST, ALIST_PS, BLIST_PS, ELIST_PS, RLIST_PS
  INTEGER, DIMENSION(MAXINT) :: IORD, IORD_PS 
  REAL(KIND = PREC) :: ZCARRY_PIGMD, FCARRY, SIGMAB, PSPECNORM,&
       &RCARRY_PSIG, SOLDLT_IGMFRAC, SOLDLT_Z, SOLFM_Z, &
       &SOLFM_IGMDCRIT, OLDIGMFRAC, SOLFM_SOURCE, X_II, &
       &SIGMA_E, SIGMA_P, OLDXII, OLDENERY, NPHDOT, ALPHA_R,&
       &SOLFM_ENERGY, SOLFM_FF, NGAMMA_LF, GROWTH, RGNSIG, &
       &ZFORM, NGAMMA_POP2, NGAMMA_POP3, stmass(DATA_COLUMNS), &
       &ekin(DATA_COLUMNS), remmas(DATA_COLUMNS), spmass(DATA_COLUMNS), zvesc,&
       &stmet(4), e_kin(4,DATA_COLUMNS), rem_mas(4,DATA_COLUMNS), sp_mass(4,DATA_COLUMNS), &
       &HaloMetalAbundance 
  INTEGER :: LCN, COUNTR, POSTOCOUNTER 
  REAL(KIND = PREC), DIMENSION(:,:), ALLOCATABLE :: LBURST, sfrarr_halocalc_hot, &
       &sfrarr_halocalc_cold, sfrarr_halocalc, sd93_lambda
  integer, dimension(:,:), allocatable :: halopop_hot, halopop_cold 
  REAL(KIND=PREC), DIMENSION(:,:,:), ALLOCATABLE :: YIELDS 

END MODULE STORAGE

