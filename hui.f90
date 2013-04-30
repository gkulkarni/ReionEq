
! File: hui.f90
!  Cre: 2012-06-13
!  Mod: $Date: 2012/06/14 13:19:58 $
!
! Implemented eq. 10 of Hui and Gnedin 1997. Can be used to check
! effect of spectrum on temperature.  Compiles with `make hui'.

module local_storage 

  use constants
  implicit none 
  real(kind=prec) :: norm 

end module local_storage

module local_interfaces

  interface 
     function sigmah(nu) 
       use constants 
       implicit none 
       real(kind=prec), intent(in) :: nu 
       real(kind=prec) :: sigmah 
     end function sigmah

     function flux(nu) 
       use constants
       use local_storage
       implicit none 
       real(kind=prec), intent(in) :: nu
       real(kind=prec) :: flux 
     end function flux

     function integrand_nr(nu) 
       use constants
       implicit none 
       real(kind=prec), intent(in) :: nu 
       real(kind=prec) :: integrand_nr 
     end function integrand_nr

     function integrand_dr(nu) 
       use constants
       implicit none 
       real(kind=prec), intent(in) :: nu 
       real(kind=prec) :: integrand_dr
     end function integrand_dr

     function integrand_norm_nr(nu) 
       use constants
       implicit none 
       real(kind=prec), intent(in) :: nu 
       real(kind=prec) :: integrand_norm_nr
     end function integrand_norm_nr

     function integrand_norm_dr(nu) 
       use constants
       implicit none 
       real(kind=prec), intent(in) :: nu 
       real(kind=prec) :: integrand_norm_dr
     end function integrand_norm_dr
  end interface

end module local_interfaces

program hui

  use constants 
  use local_storage 
  use local_interfaces
  implicit none 
  integer :: inf, last, neval, ier  
  real(kind=prec) :: error, j_HI, sum_nr, sum_dr, energy, temperature 
  real(kind=prec) :: t1, t2, nu1, nu2, numax 
  real(kind = prec), dimension(maxint) :: alist, blist, elist, rlist 
  integer, dimension(maxint) :: iord

  j_HI = 1.5e-21_prec ! erg cm^-2 sr^-2 Hz^-1 s^-1 
  numax = 1.0e20_prec ! Hz 

  norm = 1.0_prec ! dimensionless 
  nu1 = nu0 
  sum_nr = 0.0_prec
  t1 = integrand_norm_nr(nu1) 
  do 
     nu2 = nu1*1.26_prec 
     if (nu2 > numax) exit 
     t2 = integrand_norm_nr(nu2) 
     sum_nr = sum_nr + 0.5_prec*(nu2-nu1)*(t2+t1) 
     t1 = t2 
     nu1 = nu2 
  end do

  nu1 = nu0 
  sum_dr = 0.0_prec
  t1 = integrand_norm_dr(nu1) 
  do 
     nu2 = nu1*1.26_prec 
     if (nu2 > numax) exit 
     t2 = integrand_norm_dr(nu2) 
     sum_dr = sum_dr + 0.5_prec*(nu2-nu1)*(t2+t1) 
     t1 = t2 
     nu1 = nu2 
  end do

  norm = j_HI/(sum_nr/sum_dr)

  nu1 = nu0 
  sum_nr = 0.0_prec
  t1 = integrand_nr(nu1) 
  do 
     nu2 = nu1*1.26_prec 
     if (nu2 > numax) exit 
     t2 = integrand_nr(nu2) 
     sum_nr = sum_nr + 0.5_prec*(nu2-nu1)*(t2+t1) 
     t1 = t2 
     nu1 = nu2 
  end do

  nu1 = nu0 
  sum_dr = 0.0_prec
  t1 = integrand_dr(nu1) 
  do 
     nu2 = nu1*1.26_prec 
     if (nu2 > numax) exit 
     t2 = integrand_dr(nu2) 
     sum_dr = sum_dr + 0.5_prec*(nu2-nu1)*(t2+t1) 
     t1 = t2 
     nu1 = nu2 
  end do

  energy = sum_nr/sum_dr ! J 
  temperature = energy/(3.0_prec*kboltz) ! K 
  print *, temperature 

end program hui

function sigmah(nu) 

  use constants 
  implicit none 
  real(kind=prec), intent(in) :: nu 
  real(kind=prec) :: sigmah ! cm^2 
  real(kind=prec) :: a0, eps  

  ! Osterbrock book 
  a0 = 6.3e-18_prec ! cm^2 
  if (nu<nu0) then 
     sigmah = 0.0_prec 
  else if ((nu/nu0-1.0_prec)==0.0_prec) then
     sigmah = a0 * (nu0/nu)**4 
  else 
     eps = sqrt(nu/nu0-1.0_prec)
     sigmah = a0 * (nu0/nu)**4 * exp(4.0_prec - 4.0_prec*atan(eps)/eps)&
          &/(1.0_prec-exp(-2.0_prec*pi/eps))
  end if

end function sigmah

function flux(nu) 

  use constants
  use local_storage
  implicit none 
  real(kind=prec), intent(in) :: nu
  real(kind=prec) :: flux 
  real(kind=prec) :: f_reduction, index 

  f_reduction = 0.01_prec 
  index = -1.0_prec 

  if (nu < nu_HeII) then 
     flux = (nu/nu0)**index
  else
     flux = f_reduction*(nu/nu0)**index
  end if

  flux = norm*flux*1.0e-21_prec ! erg cm^-2 sr^-1 s^-1 Hz^-1 

end function flux

function integrand_nr(nu) 

  use constants
  use local_interfaces, only : sigmah, flux 
  implicit none 
  real(kind=prec), intent(in) :: nu 
  real(kind=prec) :: integrand_nr 

  integrand_nr = 4.0_prec*pi*flux(nu)*sigmah(nu)*hplanck*(nu-nu0)/(1.0e7_prec*hplanck*nu) ! J Hz^-1 

end function integrand_nr

function integrand_dr(nu) 

  use constants
  use local_interfaces, only : sigmah, flux
  implicit none 
  real(kind=prec), intent(in) :: nu 
  real(kind=prec) :: integrand_dr

  integrand_dr = 4.0_prec*pi*flux(nu)*sigmah(nu)/(1.0e7_prec*hplanck*nu) ! Hz^-1

end function integrand_dr

function integrand_norm_nr(nu) 

  use constants
  use local_interfaces, only : sigmah, flux
  implicit none 
  real(kind=prec), intent(in) :: nu 
  real(kind=prec) :: integrand_norm_nr

  integrand_norm_nr = 4.0_prec*pi*flux(nu)*sigmah(nu)/nu ! erg Hz^-1

end function integrand_norm_nr

function integrand_norm_dr(nu) 

  use constants
  use local_interfaces, only : sigmah
  implicit none 
  real(kind=prec), intent(in) :: nu 
  real(kind=prec) :: integrand_norm_dr

  integrand_norm_dr = 4.0_prec*pi*sigmah(nu)/nu ! sr cm^2 s

end function integrand_norm_dr

