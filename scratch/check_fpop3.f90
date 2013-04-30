program check_fpop3

  implicit none 
  interface
     double precision function fpop3(zcoll,m)
       implicit none
       double precision, intent(in) :: zcoll,m
     end function fpop3
  end interface

  double precision :: zcoll = 10.0, m 

  m = 1.0d-3
  do 
     if (m > 1.0d2) exit 
     m = m*1.2d0
     print *, m, fpop3(zcoll, log10(m))
  end do

end program check_fpop3
