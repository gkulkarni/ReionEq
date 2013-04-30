program try 

  implicit none 
  integer :: a(2,2), b(2) 

  a(1,1) = 1 
  a(1,2) = 2
  a(2,1) = 3
  a(2,2) = 4

  b(1) = 5
  b(2) = 6 

  print *, b 
  b = a(1,:) 
  print *, b 

end program try
