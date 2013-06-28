program deltat 

use constants 
use storage 
use interfaces, only : dtdz 

real(kind=prec) :: z, dt, dz2 

dz2 = -0.1 

z = initial_redshift 
do 
   z = z + dz2 
   if (z < final_redshift) exit 
   dt = dz2*dtdz(z) ! yr 
   print *, z, dz2, dt 
end do

end program deltat

