program thetaofz 

! Cre: 2010
! Mod: $Date$ ($Revision$) 

! This code is used in tabz.f90.  See that file for details. 

  use constants; implicit none 
  interface
     function theta(delta_l, redshift)
       use constants; use stotabz; implicit none
       interface
          function fth(x)
            use constants;
            implicit none
            real(kind = prec), intent(in) :: x
            real(kind = prec) :: fth
          end function fth
       end interface
       real(kind = prec), intent(in) :: delta_l, redshift
       real(kind = prec) :: theta
     end function theta
  end interface

  real(kind=prec) :: deltal, z, th 

  deltal = 7.676_prec
  z = 8.0 

  th = theta(deltal,z)
  print *, th

end program thetaofz
