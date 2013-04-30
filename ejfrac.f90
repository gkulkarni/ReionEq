
function ejfrac(z, pop) 

  ! File: ejfrac.f90 
  !  Cre: 2012-06-20
  !  Mod: $Date: 2012/10/23 08:10:11 $ ($Revision: 1.6 $) 

  ! Calculates the return fraction of a given IMF. Assumes IRA. 

  use constants 
  use storage 
  use interfaces, only : imf, imf_pop3, interpolate2, bi_interpolate2, getmet, getsfr_hot, getsfr_cold
  implicit none 
  real(kind=prec), intent(in) :: z 
  integer, intent(in) :: pop
  real(kind=prec) :: ejfrac

  real(kind=prec) :: t, m_d, m_up, m1, m2, int1, int2, sum, tmyr      
  integer :: i, j 
  logical :: all_species 

  call interpolate2(tarr, zarr, z, t) ! [t] = yr 
  tmyr = t*1.0e-6_prec ! Myr 
  call interpolate2(stellar_mass, stellar_age, tmyr, m_d) ! [md] = M_solar
  m_up = STMASS_UPLIMIT ! M_solar 

  do i = 1, 4
     stmet(i) = yields(i,1,1)
     do j = 1, 10
        rem_mas(i,j) = yields(i,j,3) ! Remnant mass; M_solar 
     end do
  end do

  m1 = max(m_d,8.0_prec)
  int1 = integrand(m1) 
  sum = 0.0_prec 
  do 
     m2 = m1*STELLAR_INTEGRAL_MULTIPLIER 
     if (m2 > m_up) exit 

     int2 = integrand(m2) 
     sum = sum + 0.5_prec*(m2-m1)*(int1+int2) ! dimensionless
     int1 = int2 
     m1 = m2
  end do

  ejfrac = sum ! dimensionless

contains 

  function integrand(m) 

    implicit none 
    real(kind=prec), intent(in) :: m 
    real(kind=prec) :: integrand
    real(kind=prec) :: mr, st_age, rs, psi2, psi3, mej, zmet, &
         &lzmet, st_ageyr, integrand2, integrand3 
    logical :: m_overflow, m_underflow, zmet_overflow, zmet_underflow 

    call interpolate2(stellar_age, stellar_mass, m, st_age) ! [st_age] = Myr
    st_ageyr = st_age*1.0e6_prec ! yr 
    call interpolate2(zarr, tarr, t-st_ageyr, rs) 
    zmet = getmet(rs) 

    m_overflow = .false. 
    m_underflow = .false. 
    if (m > data_mmax) then 
       m_overflow = .true. 
    else if (m < data_mmin) then 
       m_underflow = .true. 
    end if

    zmet_overflow = .false. 
    zmet_underflow = .false. 
    if (zmet > data_zmetmax) then 
       zmet_overflow = .true. 
    else if (zmet < data_zmetmin) then 
       zmet_underflow = .true. 
    end if

    if (m_overflow) then 
       if (zmet_overflow) then 
          call bi_interpolate2(stmet, stmass, rem_mas, data_zmetmax, data_mmax, mr)
          mr = (mr/data_mmax)*m
       else if (zmet_underflow) then 
          call bi_interpolate2(stmet, stmass, rem_mas, data_zmetmin, data_mmax, mr)
          mr = (mr/data_mmax)*m
       else 
          call bi_interpolate2(stmet, stmass, rem_mas, zmet, data_mmax, mr)
          mr = (mr/data_mmax)*m
       end if
    else if (m_underflow) then 
       if (zmet_overflow) then 
          call bi_interpolate2(stmet, stmass, rem_mas, data_zmetmax, data_mmin, mr)
          mr = (mr/data_mmin)*m 
       else if (zmet_underflow) then 
          call bi_interpolate2(stmet, stmass, rem_mas, data_zmetmin, data_mmin, mr)
          mr = (mr/data_mmin)*m 
       else 
          call bi_interpolate2(stmet, stmass, rem_mas, zmet, data_mmin, mr)
          mr = (mr/data_mmin)*m 
       end if
    else
       if (zmet_overflow) then 
          call bi_interpolate2(stmet, stmass, rem_mas, data_zmetmax, m, mr)
       else if (zmet_underflow) then 
          call bi_interpolate2(stmet, stmass, rem_mas, data_zmetmin, m, mr)
       else 
          call bi_interpolate2(stmet, stmass, rem_mas, zmet, m, mr)
       end if
    end if

    if (pop==3) then 
       integrand = imf_pop3(m)*(m-mr) ! dimensionless
    else if (pop==2) then 
       if (m < POP2_UPLIMIT) then 
          integrand = imf(m)*(m-mr) ! dimensionless
       else 
          integrand = 0.0_prec 
       end if
    else
       write (0,*) 'outfrac.f90: pop has an illegal value!'
    end if

!!$    if (HaloMetalAbundance < METALLICITY_POP3TRANS) then 
!!$       integrand = imf_pop3(m)*(m-mr) ! dimensionless
!!$    else
!!$       if (m < POP2_UPLIMIT) then 
!!$          integrand = imf(m)*(m-mr) ! dimensionless
!!$       else 
!!$          integrand = 0.0_prec 
!!$       end if
!!$    end if

!!$    if (z > REDSHIFT_POP3TRANS) then 
!!$       integrand = imf_pop3(m)*(m-mr) ! dimensionless
!!$    else 
!!$       integrand = imf(m)*(m-mr) ! dimensionless
!!$    end if

!!$    integrand3 = imf_pop3(m)*(m-mr) ! dimensionless
!!$    integrand2 = imf(m)*(m-mr) ! dimensionless
!!$    ! integrand = integrand2 + integrand3 ! dimensionless
!!$    integrand = integrand2 ! dimensionless

  end function integrand

end function ejfrac

