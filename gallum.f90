
subroutine gallum(HaloTotalMass, HaloAge, HaloFormationRedshift, HaloLuminosity)

  ! File: gallum.f90
  ! Cre: 2010-05-17
  ! Mod: $Date: 2010/06/05 05:27:19 $; $Revision: 1.2 $ 
  ! 
  ! Calculate luminosity of a halo with given mass, age and 
  ! formation redshift.  Eqn. 6 from Samui et al. 2007 is used. 

  use constants 
  use storage 
  implicit none 

  real(kind=prec), intent(in) :: HaloTotalMass, HaloAge, HaloFormationRedshift 
  real(kind=prec), intent(out) :: HaloLuminosity 

  real(kind=prec) :: BurstInitialTime, int1, int2, t1, t2, sum, lambda, &
       &dlambdadnu
  integer :: i 

  BurstInitialTime = lburst(1,1)
  int1 = lburst(2,1)*HaloSFR(HaloAge-BurstInitialTime) ! erg/s/Ang
  sum = BurstInitialTime*int1 
  t1 = BurstInitialTime 
  do i = 2, size(lburst,2)
     t2 = lburst(1,i) 
     if (t2 > HaloAge) exit 
     int2 = lburst(2,i)*HaloSFR(HaloAge-t2)
     sum = sum + (t2-t1)*int2
     t1=t2
  end do

  lambda = 1500 ! Ang 
  dlambdadnu = lambda**2/(speed_of_light*cmbyang) ! Ang s 
  HaloLuminosity = sum*dlambdadnu ! erg/s/Hz

contains 

  function HaloSFR(SourceAge)

    implicit none
    real(kind=prec), intent(in) :: SourceAge
    real(kind=prec) :: HaloSFR 
    
    real(kind=prec) :: tdyn, tkernel

    tdyn = DynamicalTime(HaloFormationRedshift) ! yr
    tkernel = (SourceAge/tdyn**2)*exp(-SourceAge/tdyn) ! yr^-1
    HaloSFR = tkernel*HaloTotalMass*fstar*omega_b/omega_nr ! 10^10 M_solar/yr 
    
  end function HaloSFR

  function DynamicalTime(z)

    real(kind = prec), intent(in) :: z
    real(kind = prec) :: DynamicalTime

    real(kind = prec) :: hdens, bkgrho, rhocs

    rhocs = 1.879e-29_prec*smallh**2 ! g/cm^3
    bkgrho = rhocs*omega_nr*(1.0_prec+z)**3 ! g/cm^3 
    hdens = delta_virial*bkgrho ! g/cm^3 
    ! Below, 1.0e3 converts g to kg and cm to m.  
    DynamicalTime = sqrt(3.0_prec*pi/(16.0_prec*newtg*&
         &hdens*1.0e3_prec))/yrbys ! yr

  end function DynamicalTime

end subroutine gallum 

