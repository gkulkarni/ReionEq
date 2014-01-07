program test_MergerRate 

  ! File: Test_MergerRate.f90 
  !  Cre: 2014-01-07
  
  ! This code tests code in merger_rate.f90.  In other words, it
  ! reproduces Figure 2 of Fakhouri et al. 2010 (MNRAS 406, 2267).
  ! Compiles with `make test_mergerrate'.

  use constants 
  use interfaces, only : MajorMergerRate
  implicit none 

  real(kind=prec) :: hmass, mrate, redshift 

  redshift = 0.0_prec 
  hmass = 1.0e10_prec 
  do 
     if (hmass > 1.0e15_prec) exit 
     mrate = MajorMergerRate(hmass, redshift) 
     hmass = hmass * 1.6 
     print *, hmass, mrate 
  end do

end program test_MergerRate
