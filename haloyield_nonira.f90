function haloyield_nonira(z, hotcold, bin) 

  ! File: haloyield_nonira.f90 
  !  Cre: 2012-06-21
  !  Mod: $Date: 2012/10/23 08:10:11 $ ($Revision: 1.3 $) 

  use constants 
  use storage 
  use interfaces, only : getsfr, imf, imf_pop3, interpolate2, bi_interpolate2, &
       &getmet, sfr_rollinde_pop3, sfr_rollinde_pop2, getsfr_hot, getsfr_cold, getpop_hot, getpop_cold
  implicit none 
  real(kind=prec), intent(in) :: z 
  integer, intent(in) :: hotcold, bin  
  real(kind=prec) :: haloyield_nonira

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
     do j = 1, DATA_COLUMNS
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
     sum = sum + 0.5_prec*(m2-m1)*(int1+int2) ! M_solar/yr (because we divide by 1 M_solar)
     int1 = int2 
     m1 = m2
  end do

  haloyield_nonira = sum * 1.0e-10_prec ! 10^10 M_solar / yr 

contains 

  function integrand(m) 

    implicit none 
    real(kind=prec), intent(in) :: m 
    real(kind=prec) :: integrand
    real(kind=prec) :: mr, st_age, rs, psi2, psi, mej, zmet, &
         &lzmet, st_ageyr, integrand2, integrand3
    integer :: pop
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

    if (hotcold == 1) then 
       ! cold
       pop = getpop_cold(rs, bin)
       psi = getsfr_cold(rs, bin)  

       if (pop == 2) then 
          integrand = imf(m) * psi * mej
          if (m > POP2_UPLIMIT) integrand = 0.0_prec 
       else 
          integrand = imf_pop3(m) * psi * mej
       end if
    else 
       ! hot 
       pop = getpop_hot(rs, bin)
       psi = getsfr_hot(rs, bin)  
       if (pop == 2) then 
          integrand = imf(m) * psi * mej ! M_solar / yr 
          if (m > POP2_UPLIMIT) integrand = 0.0_prec 
       else 
          integrand = imf_pop3(m) * psi * mej
       end if
    end if

  end function integrand

end function haloyield_nonira

