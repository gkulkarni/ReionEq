
function rtnewt(funcd, first_guess, xlo, xhi, tol, caller) 

  ! find root by bisection method.
  use constants
  use storage
  implicit none 
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

  integer, parameter :: jmax=100 
  real(kind=prec) :: x, dx, dxold, f, df 
  integer :: i 

  dxold = xhi - xlo 
  dx = dxold 

  call funcd(first_guess, f, df) 

  do i = 1, jmax
     if (((x-xhi)*df-f)*((x-xlo)*df-f)>=0.0_prec .or. &
          &abs(2.0_prec*f)>abs(dxold*df)) then
        dxold = dx 
        dx = 0.5_prec*(xhi-xlo)
        x = xlo+dx
     else
        dxold = dx 
        dx = f/df
        x = x-dx
     end if

     if (abs(dx)<tol) exit 
     
     call funcd(x, f, df)
     if (f < 0.0_prec) then 
        xlo = x 
     else
        xhi = x 
     end if
  end do

  rtnewt = x 

END FUNCTION RTNEWT

