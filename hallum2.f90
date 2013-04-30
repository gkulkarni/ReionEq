
function hallum2(haloindex, redshiftindex, hotcold) 

  ! File: halolum2.f90 
  !  Cre: 2013-04-23
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
  use interfaces, only : interpolate2
  implicit none 
  integer, intent(in) :: haloindex, redshiftindex, hotcold
  real(kind=prec) :: hallum2 

  real(kind=prec) :: sfr, lum, z, strlum, str_age, haloz, current_time, t, str_mass, &
       &dlambdadnu, lambda, BurstInitialTime, int1, t1, t2, sum, int2 
  real(kind=prec), dimension(:), allocatable :: lburst_age_array, diff 
  integer :: i, j, loc(1) 

  if ((hotcold /= 1) .AND. (hotcold /= 2)) then 
     write (0, *) 'hallum2: hotcold should be 1 or 2' 
     stop
  end if

  ! Copy SFR data
  if (hotcold == 1) then 
     sfrarr_halocalc = sfrarr_halocalc_cold 
     sfrarr_halocalc = sfrarr_halocalc_hot 
  end if

  ! Calculate cosmic time
  haloz = (redshiftindex-1)*dz + initial_redshift
  call interpolate2(tarr, zarr, haloz, current_time) ! [current_time]=yr

  ! Integrate over lburst, weighted by halo SFR
  BurstInitialTime = lburst(1,1)
  int1 = lburst(2,1)*HaloSFR(current_time-BurstInitialTime) ! (erg/s/Ang)/yr
  sum = BurstInitialTime*int1 
  t1 = BurstInitialTime 
  do i = 2, size(lburst,2)
     t2 = lburst(1,i) 
     if (t2 > current_time) exit 
     int2 = lburst(2,i)*HaloSFR(current_time-t2)
     sum = sum + (t2-t1)*int2 ! erg/s/Ang
     t1=t2
  end do

  lambda = 1500 ! Ang 
  dlambdadnu = lambda**2/(speed_of_light*cmbyang) ! Ang s 
  hallum2 = sum*dlambdadnu ! erg/s/Hz

contains

  function HaloSFR(source_time)

    implicit none
    real(kind=prec), intent(in) :: source_time ! yr 
    real(kind=prec) :: HaloSFR 

    real(kind=prec) :: rs
    integer :: sfrindex

    call interpolate2(zarr, tarr, source_time, rs) 
    sfrindex = (rs-initial_redshift)/dz + 1 
    HaloSFR = 1.0e-10_prec*sfrarr_halocalc(sfrindex, haloindex) ! M_solar/yr 

  end function HaloSFR

end function hallum2

