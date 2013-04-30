function getsfr(rs)

  use constants; use storage; implicit none
  real(kind = prec), intent(in) :: rs
  real(kind = prec) :: getsfr
  integer :: loc 

  loc = floor((initial_redshift - rs) / abs(dz)) + 1
  if (loc < 1) then
     getsfr = 0.0_prec
  else
     getsfr = sfrarr(loc)
  end if

end function getsfr

function getsfr2(rs)

  use constants; use storage; implicit none
  real(kind = prec), intent(in) :: rs
  real(kind = prec) :: getsfr2
  integer :: loc 

  loc = floor((initial_redshift - rs) / abs(dz)) + 1
  if (loc < 1) then
     getsfr2 = 0.0_prec
  else
     getsfr2 = sfrarr_pop2(loc)
  end if

end function getsfr2

function getsfr3(rs)

  use constants; use storage; implicit none
  real(kind = prec), intent(in) :: rs
  real(kind = prec) :: getsfr3
  integer :: loc 

  loc = floor((initial_redshift - rs) / abs(dz)) + 1
  if (loc < 1) then
     getsfr3 = 0.0_prec
  else
     getsfr3 = sfrarr_pop3(loc)
  end if

end function getsfr3

function getsfr_hot(rs, halomassbin)

  use constants; use storage; implicit none
  real(kind = prec), intent(in) :: rs
  integer, intent(in) :: halomassbin
  real(kind = prec) :: getsfr_hot
  integer :: loc 

  loc = floor((initial_redshift - rs) / abs(dz)) + 1
  getsfr_hot = sfrarr_halocalc_hot(loc, halomassbin)

end function getsfr_hot

function getsfr_cold(rs, halomassbin)

  use constants; use storage; implicit none
  real(kind = prec), intent(in) :: rs
  integer, intent(in) :: halomassbin 
  real(kind = prec) :: getsfr_cold
  integer :: loc 

  loc = floor((initial_redshift - rs) / abs(dz)) + 1
  getsfr_cold = sfrarr_halocalc_cold(loc, halomassbin)

end function getsfr_cold

function getpop_hot(rs, halomassbin)

  use constants; use storage; implicit none
  real(kind = prec), intent(in) :: rs
  integer, intent(in) :: halomassbin
  integer :: getpop_hot
  integer :: loc 

  loc = floor((initial_redshift - rs) / abs(dz)) + 1
  getpop_hot = halopop_hot(loc, halomassbin)

end function getpop_hot

function getpop_cold(rs, halomassbin)

  use constants; use storage; implicit none
  real(kind = prec), intent(in) :: rs
  integer, intent(in) :: halomassbin
  integer :: getpop_cold
  integer :: loc 

  loc = floor((initial_redshift - rs) / abs(dz)) + 1
  getpop_cold = halopop_cold(loc, halomassbin)

end function getpop_cold
