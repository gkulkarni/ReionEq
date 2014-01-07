function MajorMergerRate(HaloMass, Redshift) 

  ! File: merger_rate.f90 
  !  Cre: 2014-01-01 
  !  
  ! Returns major merger rate (number of mergers per halo per unit
  ! redshift).  It is obtained by integrating the internal function
  ! MergerRate; see comments therein.  This quantity is plotted in
  ! Figure 2 of Fakhouri et al. 2010 (MNRAS 406, 2267).  
  ! 
  ! Note: (1) MassRatioForMajorMergers decides which merger is major,
  !       (2) Input halo mass must be in Msun.

  use constants
  implicit none 
  real(kind = prec), intent(in) :: HaloMass, Redshift 
  real(kind = prec) :: MajorMergerRate 

  real(kind = prec) :: MassRatioForMajorMergers, dMassRatio, & 
       MassRatio1, MassRatio2, sum, Integrand1, Integrand2

  MassRatioForMajorMergers = 0.3_prec  
  dMassRatio = 0.1_prec 

  MassRatio1 = MassRatioForMajorMergers 
  sum = 0.0_prec 
  Integrand1 = MergerRate(HaloMass, Redshift, MassRatio1) 
  do 
     if (MassRatio1 >= 1.0_prec) exit 
     MassRatio2 = MassRatio1 + dMassRatio 
     Integrand2 = MergerRate(HaloMass, Redshift, MassRatio2)
     sum = sum + 0.5*(Integrand1+Integrand2)*dMassRatio 
     MassRatio1 = MassRatio2 
  end do

  MajorMergerRate = sum ! dimensionless 

contains 

  function MergerRate(mhalo, z, xi) 

    ! mhalo = descendant halo mass (in Msun) 
    !     z = redshift 
    !    xi = progenitor mass ratio (mass of bigger progenitor / mass
    !         of lower progenitor)
    ! 
    ! This function encodes Equation (1) from Fakhouri et al. 2010
    ! (MNRAS 406, 2267).  It returns the merger rate of a halo with
    ! above attributes in units of number of mergers per halo per unit
    ! redshift per unit xi.

    implicit none 
    real(kind = prec), intent(in) :: mhalo, z, xi 
    real(kind = prec) :: MergerRate

    real(kind = prec) :: a, alpha, beta, gamma, eta, xitilde 

    a = 0.0104_prec
    alpha = 0.133_prec 
    beta = -1.995_prec
    gamma = 0.263_prec
    eta = 0.0993_prec 
    xitilde = 9.72e-3_prec

    MergerRate = a * ((mhalo/1.0e12_prec)**alpha) * (xi**beta) * &
         exp((xi/xitilde)**gamma) * ((1.0+z)**eta)

  end function MergerRate

end function MajorMergerRate

