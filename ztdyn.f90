
subroutine ztdyn() 
  
  ! File: ztdyn.f90 
  ! Cre: 05-06-2010
  ! Mod: $Date: 2012/02/09 09:34:38 $, $Revision: 1.2 $
  !
  ! Calculate dynamical time at various redshifts.

  use constants
  use storage 
  use interfaces, only : interpolate2 
  implicit none 
  
  real(kind=prec) :: z, t, td, tn, zn

  z = 12.0_prec 
  do 
     if (z < 6.0_prec) exit 
     call interpolate2(tarr, zarr, z, t) ! [tform] = yr
     td = tdyn(z) 
     tn = t + td 
     call interpolate2(zarr, tarr, tn, zn)
     print *, z, zn, td
     z = z + dz
  end do

contains 

  function tdyn(z)

    real(kind = prec), intent(in) :: z
    real(kind = prec) :: tdyn

    real(kind = prec) :: hdens, bkgrho, rhocs

    rhocs = 1.879e-29_prec * smallh ** 2 ! g/cm^3
    bkgrho = rhocs * omega_nr * (1.0_prec + z) ** 3 ! g/cm^3 
    hdens = delta_virial * bkgrho ! g/cm^3 
    tdyn = sqrt(3.0_prec * pi / (16.0_prec * newtg * &
         &hdens * 1.0e3_prec)) / yrbys ! yr
    ! 1.0e3 converts g to kg and cm to m.  

  end function tdyn

end subroutine ztdyn
