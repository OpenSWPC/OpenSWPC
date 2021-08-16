module m_prem

  implicit none
  integer, parameter :: DP = 8

contains
 
  !! ------------------------------------------------------------------------------------------ !!
  !> Return the modified PREM 1D velocity structure
  !! 
  !! Compared to the original model of the Diewonski and Anderson (1981), the seawater layer is
  !! substituted to a crustal structure:
  !! rho=2.6 g/cm^3, vp=5.8 km/s, vs=3.2 km/s qmu=600, qkappa = 57823
  !!
  !! This routine returns Q_P and Q_S rather than Q_kappa and Q_mu, by converting using 
  !! a formula of (9.59) & (9.60) of Dahlen and Tromp (1998)
  !--
  subroutine prem(z, rho, vp, vs, qp, qs)

    real(DP), intent(in)  :: z      ! depth
    real(DP), intent(out) :: rho    ! density [g/cm^3]
    real(DP), intent(out) :: vp     ! P wavespeed [km/s]
    real(DP), intent(out) :: vs     ! S wavespped [km/s]
    real(DP), intent(out) :: qp     ! Qp
    real(DP), intent(out) :: qs     ! Qs
    !--
    real(DP)              :: r, x   ! Radius, Normalized Radius
    real(DP) :: qk, qm ! Q_kappa, Q_mu
    real(DP) :: qki, qmi, qpi, qsi
    !----    
  
    r = 6371.0_DP - z
    x = r / 6371.0_DP
  
    if( 0.0 <= r ) then
  
      if      ( r <= 1221.5_DP ) then                          ! Inner Core
        rho = 13.0885_DP              -  8.8381_DP*x*x
        vp  = 11.2622_DP              -  6.3640_DP*x*x
        vs  =  3.6678_DP              -  4.4475_DP*x*x
        qm  = 84.6_DP
        qk  = 1327.7_DP
      else if ( r <= 3480.0_DP ) then                                ! Outer Core
        rho = 12.5815_DP -  1.2638_DP*x -  3.6426_DP*x*x -  5.5281_DP*x*x*x
        vp  = 11.0487_DP -  4.0362_DP*x +  4.8023_DP*x*x - 13.5732_DP*x*x*x
        vs  = 0.0_DP
        qm  = 1d100
        qk  = 57823_DP
      else if ( r <= 3630.0_DP ) then                            ! Lower Mantle 1
        rho =  7.9565_DP -  6.4761_DP*x +  5.5283_DP*x*x -  3.0807_DP*x*x*x
        vp  = 15.3891_DP -  5.3181_DP*x +  5.5242_DP*x*x -  2.5514_DP*x*x*x
        vs  =  6.9254_DP +  1.4672_DP*x -  2.0834_DP*x*x +  0.9783_DP*x*x*x
        qm  = 312_DP
        qk  = 57823_DP
      else if ( r <= 5600.0_DP ) then                            ! Lower Mantle 2
        rho =  7.9565_DP -  6.4761_DP*x +  5.5283_DP*x*x -  3.0807_DP*x*x*x
        vp  = 24.9520_DP - 40.4673_DP*x + 51.4832_DP*x*x - 26.6419_DP*x*x*x
        vs  = 11.1671_DP - 13.7818_DP*x + 17.4575_DP*x*x -  9.2777_DP*x*x*x
        qm  = 312_DP
        qk  = 57823_DP
      else if ( r <= 5701.0_DP ) then                            ! Lower Mantle 3
        rho =  7.9565_DP -  6.4761_DP*x +  5.5283_DP*x*x -  3.0807_DP*x*x*x
        vp  = 29.2766_DP - 23.6027_DP*x +  5.5242_DP*x*x -  2.5514_DP*x*x*x
        vs  = 22.3459_DP - 17.2473_DP*x -  2.0834_DP*x*x +  0.9783_DP*x*x*x
        qm  = 312_DP
        qk  = 57823_DP
      else if ( r <= 5771.0_DP ) then                         ! Transition zone 1
        rho =  5.3197_DP -  1.4836_DP*x
        vp  = 19.0957_DP -  9.8672_DP*x
        vs  =  9.9839_DP -  4.9324_DP*x
        qm  = 143_DP
        qk  = 57823_DP
      else if ( r <= 5971.0_DP ) then                         ! Transition zone 2
        rho = 11.2494_DP -  8.0298_DP*x
        vp  = 39.7027_DP - 32.6166_DP*x
        vs  = 22.3512_DP - 18.5856_DP*x
        qm  = 143_DP
        qk  = 57823_DP
      else if ( r <= 6151.0_DP ) then                         ! Transition zone 3
        rho =  7.1089_DP -  3.8045_DP*x
        vp  = 20.3926_DP - 12.2569_DP*x
        vs  =  8.9496_DP -  4.4597_DP*x
        qm  = 143_DP
        qk  = 57823_DP
      else if ( r <= 6346.6_DP ) then                                       ! LVZ
        rho =  2.6910_DP +  0.6924_DP*x
        vp  =  4.1875_DP +  3.9382_DP*x
        vs  =  2.1519_DP +  2.3481_DP*x
        qm  = 80_DP
        qk  = 57823_DP
      else if ( r <= 6346.6_DP ) then                                       ! LID
        rho =  2.6910_DP +  0.6924_DP*x
        vp  =  4.1875_DP +  3.9382_DP*x
        vs  =  2.1519_DP +  2.3481_DP*x
        qm  =  600_DP
        qk  =  57823_DP
      else if ( r <= 6356.0_DP ) then                                   ! Crust 1
        rho =  2.900_DP
        vp  =  6.800_DP
        vs  =  3.900_DP
        qm  =  600_DP
        qk  =  57823_DP
      else if ( r <= 6368.0_DP ) then                                   ! Crust 2
        rho =  2.600_DP
        vp  =  5.800_DP
        vs  =  3.200_DP
        qm  =  600_DP
        qk  =  57823_DP
      else if ( r <= 6371.0_DP ) then                                   ! Crust 3
        rho =  2.600_DP
        vp  =  5.800_DP
        vs  =  3.200_DP
        qm  =  600_DP
        qk  =  57823_DP
      end if
      
    else 
      rho = 0.0_DP
      vp  = 0.0_DP
      vs  = 0.0_DP                            ! for illegal value
      qm  = 0.0_DP
      qk  = 0.0_DP
    end if

    ! atetnuation
    qmi = 1._dp / qm
    qki = 1._dp / qm
    qpi = ( 1.0_dp - (4.0_dp/3.0_dp)*vs*vs/vp/vp ) * qki &
        + 4.0_dp/3.0_dp * vs*vs/vp/vp * qmi
    qsi = qmi
    qp = 1.0_dp / qpi
    qs = 1.0_dp / qsi

  end subroutine prem

end module m_prem