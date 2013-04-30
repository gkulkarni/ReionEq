
! File: sample_ns.f90 
!  Cre: 01-2013
!  Mod: $Date$ ($Revision$)

! Produces a random sample from DLA distribution functions created by
! reion. Uses rejction method (see Numerical Recipes).

program sample_ns

  use constants 
  implicit none 
  real(kind=prec) :: rn(1), area_under_comp_fn, cdf_comp_fn_min, cdf_comp_fn_max, xmin, xmax, y, x, &
       &bigf, cdf_min, cdf_max, reference_y 
  real(kind=prec), dimension(:), allocatable :: pdf_x, pdf_y, diff 
  integer :: i, number_of_lines, loc(1), n_points

  open(11, file='plots/pdf_p3.dat', status='old', action='read') 
  read(11, *) x, y 
  number_of_lines = 0 
  do 
     number_of_lines = number_of_lines + 1 
     read(11, *, end=911) x, y 
  end do
911 continue 
  rewind 11 

  allocate(pdf_x(number_of_lines)) 
  allocate(pdf_y(number_of_lines)) 
  allocate(diff(number_of_lines))
  do i = 1, number_of_lines 
     read(11,*) x, y 
     pdf_x(i) = x 
     pdf_y(i) = y 
  end do
  close(11) 
  diff = 0.0_prec 

  area_under_comp_fn = 0.0054_prec  
  cdf_comp_fn_min = 0.001284_prec 
  cdf_comp_fn_max = area_under_comp_fn ! by definition

  xmin = 0.1_prec 
  xmax = 0.7_prec 

  cdf_min = cdf(xmin)
  cdf_max = cdf(xmax) 

  n_points = 0 
  do 
     call random_number(rn) 
     y = rn(1)*(cdf_comp_fn_max-cdf_comp_fn_min) + cdf_comp_fn_min
     x = cdf_root(y) 
     if (x <= 0.0_prec) then 
        write (0,*) 'Error: cdf_root is negative!' 
     end if

     bigf = cdf(x) 

     call random_number(rn) 
     y = rn(1)*cdf_max

     diff = abs(pdf_x - x)
     loc = minloc(diff)
     diff = 0.0_prec 
     reference_y = pdf_y(loc(1)) 
     
     if (y > reference_y) cycle 
     n_points = n_points + 1 

     ! print *, y, x, reference_y 
     print *, x
     if (n_points > 19) exit 
  end do

contains 

  function cdf_root(y) 

    implicit none 
    real(kind=prec), intent(in) :: y 
    real(kind=prec) :: cdf_root 

    real(kind=prec) :: a, b

    a = 0.002857_prec 
    b = 0.006714_prec 
    cdf_root = (-2.0_prec*b + sqrt(4.0_prec*b*b+8.0_prec*a*y)) / (2.0_prec*a) 

  end function cdf_root

  function cdf(x) 

    implicit none 
    real(kind=prec), intent(in) :: x 
    real(kind=prec) :: cdf 

    real(kind=prec) :: a, b 

    a = 0.002857_prec 
    b = 0.006714_prec 
    cdf = a*x + b 

  end function cdf

end program sample_ns
