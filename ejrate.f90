
function ejrate(z, species) 

  ! File: ejrate.f90 
  !  Cre: 2012-01-26 
  !  Mod: $Date: 2012/10/23 08:10:11 $ ($Revision: 1.12 $) 

  ! Calculates the *average* returned fraction given an IMF and SFR.
  ! Uses the definition of `e' from Daigne et al. 2004 and 2006.

  use constants 
  use storage 
  use interfaces, only : getsfr, imf, imf_pop3, interpolate2, bi_interpolate2, &
       &getmet, sfr_rollinde_pop3, sfr_rollinde_pop2, getsfr2, getsfr3 
  implicit none 
  real(kind=prec), intent(in) :: z 
  integer, intent(in) :: species 
  real(kind=prec) :: ejrate

  real(kind=prec) :: t, m_d, m_up, m1, m2, int1, int2, sum, tmyr      
  integer :: i, j 
  logical :: all_species 

  ! Fe: species=14
  !  C: species=7
  !  O: species=9
  ! Si: species=11
  !  N: species=8
  ! Zn: species=15
  ! Mg: species=10
  ! All species combined (WITH H and He): species=0
  ! All species combined (WITHOUT H and He): species=16 

  if (species == 0) then 
     all_species = .true. 
  else 
     all_species = .false.
  end if

  call interpolate2(tarr, zarr, z, t) ! [t] = yr 
  tmyr = t*1.0e-6_prec ! Myr 
  call interpolate2(stellar_mass, stellar_age, tmyr, m_d) ! [md] = M_solar
  m_up = STMASS_UPLIMIT ! M_solar 

  do i = 1, DATA_COLUMNS 
     remmas(i) = yields(1,i,3) ! Remnant mass; M_solar 
     stmass(i) = yields(1,i,2) ! Stellar mass; M_solar 
  end do

  sp_mass = 0.0_prec
  do i = 1, 4
     stmet(i) = yields(i,1,1)
     do j = 1, DATA_COLUMNS
        rem_mas(i,j) = yields(i,j,3) ! Remnant mass; M_solar 
        if (.not. all_species) then 
           sp_mass(i,j) = yields(i,j,species+1) ! Species mass; M_solar 
        end if
     end do
  end do

  m1 = max(m_d,8.0_prec)
  int1 = integrand(m1) 
  sum = 0.0_prec 
  do 
     m2 = m1*STELLAR_INTEGRAL_MULTIPLIER 
     if (m2 > m_up) exit 
     int2 = integrand(m2) 
     sum = sum + 0.5_prec*(m2-m1)*(int1+int2) ! yr^-1 Mpc^-3 M_solar
     int1 = int2 
     m1 = m2
  end do

  ejrate = sum*1.0e-10_prec ! 10^10 M_solar yr^-1 Mpc^-3 

contains 

  function integrand(m) 

    implicit none 
    real(kind=prec), intent(in) :: m 
    real(kind=prec) :: integrand
    real(kind=prec) :: mr, st_age, rs, psi2, psi3, mej, zmet, &
         &lzmet, st_ageyr, integrand2, integrand3 
    logical :: m_overflow, m_underflow, zmet_overflow, zmet_underflow 

    if (all_species) then 
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

       psi2 = getsfr2(rs) ! M_solar yr^-1 Mpc^-3
       integrand2 = imf(m)*(m-mr)*psi2 ! yr^-1 Mpc^-3 

       psi3 = getsfr3(rs) ! M_solar yr^-1 Mpc^-3
       integrand3 = imf_pop3(m)*(m-mr)*psi3 ! yr^-1 Mpc^-3 

       if (m > POP2_UPLIMIT) integrand2 = 0.0_prec 

       integrand = integrand2 + integrand3 ! yr^-1 Mpc^-3 
    else 
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

       psi2 = getsfr2(rs) ! M_solar yr^-1 Mpc^-3
       integrand2 = imf(m)*mej*psi2 ! yr^-1 Mpc^-3 

       psi3 = getsfr3(rs) ! M_solar yr^-1 Mpc^-3
       integrand3 = imf_pop3(m)*mej*psi3 ! yr^-1 Mpc^-3 

       if (m > POP2_UPLIMIT) integrand2 = 0.0_prec 

       integrand = integrand2 + integrand3 ! yr^-1 Mpc^-3 
    end if

  end function integrand

end function ejrate

