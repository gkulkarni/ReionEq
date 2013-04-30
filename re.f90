
program re 

  ! File: re.f90 
  ! Cre: 2010-05-20 
  ! Mod: $Date: 2010/06/05 05:27:19 $; $Revision: 1.1 $ 
  ! 
  ! This program calculates the comoving Eulerian distance at redshift
  ! z corresponding to a field of view theta.  It calculates angular
  ! diameter distance for that.  Compiles with `make re'.

  use constants
  use storage 
  implicit none 

  interface 
     function chi_integrand(z) 
       use constants
       implicit none
       real(kind=prec), intent(in) :: z 
       real(kind=prec) :: chi_integrand 
     end function chi_integrand

  end interface

  real(kind=prec) :: err, sum, z, da, h0, theta, x, c, thetain, zin, dl
  integer :: flag, last, neval, comnum 
  character(len=100) :: thetaarg, zarg 

  comnum = command_argument_count()
  if (comnum == 0) then 
     write (0, '(a)') 'usage: ./re theta (in radians) z'
     write (0, '(a)') 'example: ./re 6.596e-4 8.0'
     stop
  end if

  call get_command_argument(1, thetaarg) 
  call get_command_argument(2, zarg) 

  read(thetaarg, '(f)') thetain
  read(zarg, '(f)') zin

  theta = real(thetain,prec) ! radians
  z = real(zin, prec) 

  sum = 0.0_prec 

  flag = -1 
  call dqage(chi_integrand,0.0,z,abserr,absrel,key,&
       &maxint,sum,err,neval,flag,alist,blist,rlist,elist,iord,last)
  if (flag > 0) write (0,*) 're: integration error'

  c = speed_of_light*yrbys*cmbympc ! Mpc/yr
  h0 = smallh*1.023e-10_prec ! yr^-1
  da = c*sum/(h0*(1.0_prec+z)) ! Mpc

  x = da*theta ! Mpc 
  write (*,'(a,e10.3e2)') 'Comoving Eulerian size (Mpc) =', x
  write (*,'(a,e10.3e2)') 'da (Mpc) =', da

  dl = (1.0_prec+z)**2*da 
  write (*,'(a,e10.3e2)') 'dl (Mpc) =', dl
  
end program re 

function chi_integrand(z) 

  use constants
  implicit none
  real(kind=prec), intent(in) :: z 
  real(kind=prec) :: chi_integrand

  chi_integrand = 1.0/sqrt(omega_lambda+omega_nr*(1.0_prec+z)**3)
  
end function chi_integrand

