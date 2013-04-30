function getmet(rs)

  use constants; use storage; implicit none
  real(kind = prec), intent(in) :: rs
  real(kind = prec) :: getmet
  integer :: loc 

  loc = floor((initial_redshift - rs) / abs(dz)) + 1
  if (loc < 1) then
     getmet = 0.0_prec
  else 
     getmet = metarr(loc)
  end if

end function getmet

