
function hallum(haloindex, redshiftindex, hotcold) 

  ! File: halolum.f90 
  !  Cre: 2013-03-26 
  !  Mod: $Date$ ($Revision$) 
  !
  ! Calculates the luminosity of a halo at 1500 Ang at a given time in
  ! the bathtub model.  This luminosity is calculated by adding up the
  ! luminosity of packets of stars of various ages.  For a given
  ! stellar packet of a certain age, the luminosity is calculated by
  ! looking up the array lburst, which comes from
  ! popsyn/l1500-spectrum-lt.  (See readin.inc for lburst.)

  use constants 
  use storage 
  use interfaces, only : interpolate2, dtdz
  implicit none 
  integer, intent(in) :: haloindex, redshiftindex, hotcold
  real(kind=prec) :: hallum 

  real(kind=prec) :: sfr, lum, z, strlum, str_age, haloz, current_time, t, str_mass, dlambdadnu    
  real(kind=prec), dimension(:), allocatable :: lburst_age_array, diff 
  integer :: i, j, loc(1) 

  if ((hotcold /= 1) .AND. (hotcold /= 2)) then 
     write (0, *) 'hallum: hotcold should be 1 or 2' 
     stop
  end if

  ! Copy SFR data
  if (hotcold == 1) then 
     sfrarr_halocalc = sfrarr_halocalc_cold 
     sfrarr_halocalc = sfrarr_halocalc_hot 
  end if

  ! Copy first column of lburst array. 
  allocate(lburst_age_array(size(lburst,2)))
  allocate(diff(size(lburst_age_array)))
  lburst_age_array = lburst(1,:) ! yr 

  ! Calculate halo age. 
  haloz = (redshiftindex-1)*dz + initial_redshift
  call interpolate2(tarr, zarr, haloz, current_time) ! [current_time]=yr

  lum = 0.0_prec 
  do i = 1, redshiftindex 

     ! Look up halo SFR at redshift corresponding to i 
     sfr = sfrarr_halocalc(i,haloindex) ! M_solar/yr 

     ! Calculate age of that packet today.
     z = (i-1)*dz + initial_redshift 
     call interpolate2(tarr, zarr, z, t) ! [str_age]=yr 
     str_age = current_time-t

     ! Look up corresponding luminosity. 
     str_mass = sfr*1.0e-10_prec*dz*dtdz(z) ! 10^10 M_solar
     diff = abs(str_age-lburst_age_array)
     loc = minloc(diff)
     strlum = lburst(2,loc(1))*str_mass ! erg/s/Ang 

     lum = lum + strlum ! erg/s/Ang

  enddo
  
  ! Convert luminosity units to make suitable to magnitude estimation
  dlambdadnu = 1500.0_prec**2/(speed_of_light*cmbyang) ! Ang s 
  hallum = lum*dlambdadnu ! erg/s/Hz

  deallocate(lburst_age_array) 
  deallocate(diff) 

end function hallum


