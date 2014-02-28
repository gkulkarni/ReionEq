
FUNCTION scale_sed()

  use constants 
  use storage 
  implicit none 

  real(kind=prec) :: scale_sed 
  real(kind=prec) :: sum, int1, int2
  integer :: i 

  sum = 0.0_prec 
  int1 = ionnumpop3(1) ! /Hz/(10^10 M_solar) 
  do i = 2, size(ionnumpop3) 
     if (ionfreqpop3(i)<nu0) exit 
     int2 = ionnumpop3(i) 
     sum = sum + 0.5*(ionfreqpop3(i-1)-ionfreqpop3(i))*(int1+int2)
     int2 = int1 
  end do

  scale_sed = sum ! (10^10 M_solar)^-1

END FUNCTION scale_sed
