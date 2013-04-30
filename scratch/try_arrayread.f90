program try 

  implicit none 
  double precision :: foo(10) 
  open(unit=11,file="yields.ww95bis.Z01",action="read")
  do 
     read (11,*,end=911) foo
     print *, foo 
  end do
911 continue 
  close(11)

end program try


