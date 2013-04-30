program try_biinterpolate

  use constants 
  implicit none 
  interface 
     subroutine bi_interpolate2(xarray, yarray, zarray, x, y, z)
       use constants; implicit none
       real(kind=prec), dimension(:), intent(in) :: xarray, yarray
       real(kind=prec), dimension(:,:) :: zarray
       real(kind=prec), intent(in) :: x, y
       real(kind=prec), intent(out) :: z
     end subroutine bi_interpolate2
  end interface

  real(kind=prec) :: y(3,4), x1a(3), x2a(4), x1, x2, r 
  integer :: i, j 

  do i=1,3
     x1a(i) = real(i,prec)
     do j=1,4 
        y(i,j) = real(i*j,prec) 
        x2a(j) = real(j,prec) 
     end do
  end do

  x1 = 2.0_prec 
  x2 = 3.5_prec 

  call bi_interpolate2(x1a, x2a, y, x1, x2, r)
  print *, r 

end program try_biinterpolate
