program compare_cross_sections 

  implicit none 
  real :: nu 
  
  nu = 0.3e15
  do 
     if (nu > 0.3e17) exit 
     print *, nu, sigmah(nu), sigma_hi(nu) 
     nu = nu*1.12
  end do

contains

  function sigma_hi(nu)

    implicit none
    real :: sigma_hi,nu
    real :: sigma_0,nu_0,p,x, hplanck_ev 

    hplanck_ev = 4.1356692e-15
    sigma_0=5.475d-14
    nu_0=4.298d-1/hplanck_ev
    p=2.963
    x=nu/nu_0
    sigma_hi=sigma_0*((x-1.d0)**2)*(x**(0.5*p-5.5))/(1+sqrt(x/32.88d0))**p

  end function sigma_hi

  function sigmah(nu) 

    implicit none 
    real, intent(in) :: nu 
    real :: sigmah 

    real :: a0, eps, nu0, pi 

    ! osterbrock book 
    a0 = 6.3e-18 ! cm^2 
    nu0 = 3.288e15 ! cm^2 
    pi = 3.1615

    if (nu<nu0) then 
       sigmah = 0.0 
    else if ((nu/nu0-1.0)==0.0) then
       sigmah = a0 * (nu0/nu)**4 
    else 
       eps = sqrt(nu/nu0-1.0)
       sigmah = a0 * (nu0/nu)**4 * exp(4.0 - 4.0*atan(eps)/eps)&
            &/(1.0-exp(-2.0*pi/eps))
    end if

  end function sigmah

end program compare_cross_sections
