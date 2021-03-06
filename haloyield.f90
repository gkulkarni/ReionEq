
function haloyield(z, pop) 

  ! File: haloyield.f90 
  !  Cre: 2012-06-21
  !  Mod: $Date: 2012/10/23 08:10:11 $ ($Revision: 1.6 $) 

  use constants 
  use storage 
  use interfaces, only : getsfr, imf, imf_pop3, interpolate2, bi_interpolate2, &
       &getmet, sfr_rollinde_pop3, sfr_rollinde_pop2, getsfr_hot, getsfr_cold
  implicit none 
  real(kind=prec), intent(in) :: z 
  integer, intent(in) :: pop 
  real(kind=prec) :: haloyield

  real(kind=prec) :: t, m_d, m_up, m1, m2, int1, int2, sum, tmyr      
  integer :: i, j, species = 16 ! All species combined (WITHOUT H and He): species=16 
  logical :: all_species 

  call interpolate2(tarr, zarr, z, t) ! [t] = yr 
  tmyr = t*1.0e-6_prec ! Myr 
  call interpolate2(stellar_mass, stellar_age, tmyr, m_d) ! [md] = M_solar
  m_up = STMASS_UPLIMIT ! M_solar 

  sp_mass = 0.0_prec
  do i = 1, 4
     stmet(i) = yields(i,1,1)
     do j = 1, 10
        sp_mass(i,j) = yields(i,j,species+1) ! Species mass; M_solar 
     end do
  end do

  m1 = max(m_d,8.0_prec)
  int1 = integrand(m1) 
  sum = 0.0_prec 
  do 
     m2 = m1*STELLAR_INTEGRAL_MULTIPLIER 
     if (m2 > m_up) exit 
          
     int2 = integrand(m2) 
     sum = sum + 0.5_prec*(m2-m1)*(int1+int2) ! dimensionless (because we divide by 1 M_solar)
     int1 = int2 
     m1 = m2
  end do

  haloyield = sum ! dimensionless 

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
          call bi_interpolate2(stmet, stmass, sp_mass, data_zmetmax, data_mmax, mej)
          mej = (mej/data_mmax)*m
       else if (zmet_underflow) then 
          call bi_interpolate2(stmet, stmass, sp_mass, data_zmetmin, data_mmax, mej)
          mej = (mej/data_mmax)*m
       else 
          call bi_interpolate2(stmet, stmass, sp_mass, zmet, data_mmax, mej)
          mej = (mej/data_mmax)*m
       end if
    else if (m_underflow) then 
       if (zmet_overflow) then 
          call bi_interpolate2(stmet, stmass, sp_mass, data_zmetmax, data_mmin, mej)
          mej = (mej/data_mmin)*m
       else if (zmet_underflow) then 
          call bi_interpolate2(stmet, stmass, sp_mass, data_zmetmin, data_mmin, mej)
          mej = (mej/data_mmin)*m
       else 
          call bi_interpolate2(stmet, stmass, sp_mass, zmet, data_mmin, mej)
          mej = (mej/data_mmin)*m
       end if
    else
       if (zmet_overflow) then 
          call bi_interpolate2(stmet, stmass, sp_mass, data_zmetmax, m, mej)
       else if (zmet_underflow) then 
          call bi_interpolate2(stmet, stmass, sp_mass, data_zmetmin, m, mej)
       else 
          call bi_interpolate2(stmet, stmass, sp_mass, zmet, m, mej)
       end if
    end if

!!$    if (HaloMetalAbundance < METALLICITY_POP3TRANS) then 
!!$       integrand = imf_pop3(m)*mej ! dimensionless
!!$    else
!!$       integrand = imf(m)*mej ! dimensionless
!!$    end if

    if (pop==3) then 
       integrand = imf_pop3(m)*mej ! dimensionless
    else if (pop==2) then 
       if (m < POP2_UPLIMIT) then 
          integrand = imf(m)*mej ! dimensionless
       else 
          integrand = 0.0_prec 
       end if
    else
       write (0,*) 'haloyield.f90: pop has an illegal value!'
    end if


!!$    if (z > REDSHIFT_POP3TRANS) then 
!!$       integrand = imf_pop3(m)*mej ! dimensionless
!!$    else 
!!$       if (m < POP2_UPLIMIT) then 
!!$          integrand = imf(m)*mej ! dimensionless
!!$       else 
!!$          integrand = 0.0_prec 
!!$       end if
!!$    end if

!!$    integrand2 = imf(m)*mej ! dimensionless 
!!$    integrand3 = imf_pop3(m)*mej ! dimensionless
!!$    ! integrand = integrand2 + integrand3 ! dimensionless
!!$    integrand = integrand2 ! dimensionless

  end function integrand

end function haloyield

