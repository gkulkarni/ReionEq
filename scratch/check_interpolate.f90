program check_interpolate

  implicit none 
  interface
     subroutine interpolate2(yarray, xarray, x, y)
       use constants; implicit none
       real(kind=prec), dimension(:), intent(in) :: yarray, xarray
       real(kind=prec), intent(in) :: x
       real(kind=prec), intent(out) :: y 
     end subroutine interpolate2
  end interface

  real(kind=8) :: xarr(4), yarr(4), x, y 
  integer :: i 

  do i = 1, 4 
     xarr(i) = real(i,8)
     yarr(i) = 2.0d0*xarr(i)
  end do

  x = 0.1d0
  call interpolate2(yarr, xarr, x, y) 
  print *, y 

end program
