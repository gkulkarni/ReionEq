program check_size 

  implicit none 
  integer :: foo(2,2), i, j 

  print *, size(foo)

  do i=1,2 
     do j=1,2
        foo(i,j) = 1.0 
     end do
  end do

  j = size(foo) 

  print *, j

end program check_size
