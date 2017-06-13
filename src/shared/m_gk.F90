!! ----------------------------------------------------------------------------------------------------------------------------- !!
!>
!! Geographyic coordinate transform by Gauss-Krugeur projection
!!
!! @see
!! http://www.gsi.go.jp/common/000062452.pdf
!! http://www.gsi.go.jp/common/000065826.pdf
!! http://surveycalc.gsi.go.jp/sokuchi/surveycalc/algorithm/xy2bl/xy2bl.htm
!! http://surveycalc.gsi.go.jp/sokuchi/surveycalc/algorithm/bl2xy/bl2xy.htm
!!
!! @copyright
!!   Copyright 2013-2017 Takuto Maeda. All rights reserved. This project is released under the MIT license.
!<
!! --
module m_gk

  use m_std
  implicit none
  private
  save

  real(DP) :: alpha(5)
  real(DP) :: beta(5)
  real(DP) :: AA(0:5)
  real(DP) :: delta(6)
  real(DP), parameter :: a = 6378137.0_DP !< Earth radius along the equator
  real(DP), parameter :: F = 298.257222101_DP !< Inverse elliplicity of the Earth (exact value)
  real(DP), parameter :: m0 = 0.9999_DP !< scale factor
  real(DP), parameter :: n = 1.0_DP / ( 2.0_DP*F - 1.0_DP )
  real(DP), parameter :: rho = 3600 * 180 / PI_D
  logical :: is_first = .true.
  public :: gk__lltoxy
  public :: gk__xytoll

  interface gk__lltoxy
    module procedure lltoxy_d, lltoxy_s
  end interface gk__lltoxy
  interface gk__xytoll
    module procedure xytoll_d, xytoll_s
  end interface gk__xytoll

contains

  !! --------------------------------------------------------------------------------------------------------------------------- !!
  subroutine lltoxy_d( lon, lat, lon0, lat0, x, y )

    real(DP), intent(in) :: lon
    real(DP), intent(in) :: lat
    real(DP), intent(in) :: lon0
    real(DP), intent(in) :: lat0
    real(DP), intent(out) :: x
    real(DP), intent(out) :: y
    !! --
    real(DP) :: lam, lam0, phi, phi0
    real(DP) :: e2n
    real(DP) :: lam_c, lam_s
    real(DP) :: tan_chi, cos_chi
    real(DP) :: xi, eta
    real(DP) :: Abar
    integer :: j
    !! ----


    if( is_first ) then
      call set_coef
      is_first = .false.
    end if

    !! ----

    lam  = std__deg2rad(lon)
    lam0 = std__deg2rad(lon0)
    phi  = std__deg2rad(lat)
    phi0 = std__deg2rad(lat0)

    e2n = 2.0_DP*sqrt(n)/(1.0_DP+n)

    lam_c = cos( lam - lam0 )
    lam_s = sin( lam - lam0 )

    tan_chi = sinh( atanh( sin(phi) ) - e2n * atanh( e2n * sin(phi) ) )
    cos_chi = sqrt( 1 + tan_chi**2 )
    xi      = atan( tan_chi / lam_c )
    eta     = atanh( lam_s / cos_chi )

    !! ----

    Abar = m0 * a / ( 1 + n ) * AA(0)

    !! ----

    x = xi
    y = eta

    do j=1, 5
      x = x + alpha(j)*sin(2*j*xi)*cosh(2*j*eta)
      y = y + alpha(j)*cos(2*j*xi)*sinh(2*j*eta)
    end do
    x = Abar * x - S_phi0(phi0)
    y = Abar * y

    x = x / 1000
    y = y / 1000

    !! ----
    !! @todo
    !! 子午線収差角と縮尺係数も計算する？
    !! ----

  end subroutine lltoxy_d
  !! --------------------------------------------------------------------------------------------------------------------------- !!

  !! --------------------------------------------------------------------------------------------------------------------------- !!
  subroutine lltoxy_s( lon, lat, lon0, lat0, x, y )

    real(SP), intent(in) :: lon
    real(SP), intent(in) :: lat
    real(SP), intent(in) :: lon0
    real(SP), intent(in) :: lat0
    real(SP), intent(out) :: x
    real(SP), intent(out) :: y
    !! --
    real(DP) :: xx, yy

    call lltoxy_d( dble(lon), dble(lat), dble(lon0), dble(lat0), xx, yy )
    x = real(xx)
    y = real(yy)

  end subroutine lltoxy_s
  !! --------------------------------------------------------------------------------------------------------------------------- !!


  !! --------------------------------------------------------------------------------------------------------------------------- !!
  subroutine xytoll_d( x, y, lon0, lat0, lon, lat )

    real(DP), intent(in)  :: x
    real(DP), intent(in)  :: y
    real(DP), intent(in)  :: lon0
    real(DP), intent(in)  :: lat0
    real(DP), intent(out) :: lon
    real(DP), intent(out) :: lat
    !! --

    real(DP) :: lam, lam0, phi, phi0
    real(DP) :: Abar
    real(DP) :: xi, eta, xi2, eta2
    real(DP) :: chi
    integer :: j
    if( is_first ) then
      call set_coef
      is_first = .false.
    end if


    lam0 = std__deg2rad(lon0)
    phi0 = std__deg2rad(lat0)

    Abar = m0 * a / ( 1+n) * AA(0)
    xi   = ( x*1000 + S_phi0(phi0) ) / Abar
    eta  = y*1000 / Abar

    xi2 = xi
    eta2 = eta
    do j=1, 5
      xi2  = xi2  - beta(j) * sin(2*j*xi) * cosh(2*j*eta)
      eta2 = eta2 - beta(j) * cos(2*j*xi) * sinh(2*j*eta)
    end do
    chi = asin( sin(xi2)/cosh(eta2) )

    lam = lam0 + atan( sinh(eta2)/cos(xi2) )
    phi = chi
    do j=1, 6
      phi = phi + delta(j)*sin(2*j*chi)
    end do

    lon = std__rad2deg( lam )
    lat = std__rad2deg( phi )

  end subroutine xytoll_d
  !! --------------------------------------------------------------------------------------------------------------------------- !!

  !! --------------------------------------------------------------------------------------------------------------------------- !!
  subroutine xytoll_s( x, y, lon0, lat0, lon, lat )

    real(SP), intent(in)  :: x
    real(SP), intent(in)  :: y
    real(SP), intent(in)  :: lon0
    real(SP), intent(in)  :: lat0
    real(SP), intent(out) :: lon
    real(SP), intent(out) :: lat
    !! --
    real(DP) :: lon2, lat2

    call xytoll_d( dble(x), dble(y), dble(lon0), dble(lat0), lon2, lat2 )
    lon = real(lon2)
    lat = real(lat2)

  end subroutine xytoll_s
  !! --------------------------------------------------------------------------------------------------------------------------- !!


  !! --------------------------------------------------------------------------------------------------------------------------- !!
  real(DP) function S_phi0( phi0 )
    real(DP), intent(in) :: phi0
    integer :: j

    !S_phi0 = AA(0) * phi0 / rho
    S_phi0 = AA(0) * phi0
    do j=1, 5
      S_phi0 = S_phi0 + AA(j) * sin( 2 * j * phi0 )
    end do
    S_phi0 = S_phi0 * ( m0 * a / ( 1 + n ) )
    !        S_phi0 = S_phi0 * ( a / ( 1 + n ) )

  end function S_phi0
  !! --------------------------------------------------------------------------------------------------------------------------- !!


  !! --------------------------------------------------------------------------------------------------------------------------- !!
  !>
  !! Calculate coefficients for series sum of lltoxy
  !<
  !! --
  subroutine set_coef

    !! --

    alpha(1) = (  1/ 2._DP  + (   -2/  3._DP + (  5/ 16._DP + (   41/   180._DP -   127/  288._DP *n ) *n ) *n )*n )* n**1
    alpha(2) =                (   13/48._DP  + ( -3/  5._DP + (  557/  1440._DP +   281/  630._DP *n ) *n ) *n )    * n**2
    alpha(3) =                                 ( 61/240._DP + (- 103/   140._DP + 15061/26880._DP *n ) *n )         * n**3
    alpha(4) =                                                (49561/161280._DP -   179/  168._DP *n )              * n**4
    alpha(5) =                                                                    34729/80640._DP                   * n**5
    beta(1) = (  1/2._DP  + ( -2/ 3._DP + ( 37/ 96._DP + (   -1/   360._DP -    81/   512._DP *n ) *n ) *n ) *n ) *n
    beta(2) = (             (  1/48._DP + (  1/ 15._DP + ( -437/  1440._DP +    46/   105._DP *n ) *n ) *n ) *n ) *n
    beta(3) = (             (             ( 17/480._DP + (  -37/   840._DP -   209/  4480._DP *n ) *n ) *n ) *n ) *n
    beta(4) = (             (             (              ( 4397/161280._DP -    11/   504._DP *n ) *n ) *n ) *n ) *n
    beta(5) = (             (             (              (                    4583/161280._DP *n ) *n ) *n ) *n ) *n

    delta(1) = (2/1._DP + ( -2/3._DP + (  -2/ 1._DP + (  116/ 45._DP + (   26/ 45._DP + (   -2854/  675._DP )*n)*n)*n)*n)*n)*n
    delta(2) = (          (  7/3._DP + (  -8/ 5._DP + ( -227/ 45._DP + ( 2704/315._DP + (    2323/  945._DP )*n)*n)*n)*n)*n)*n
    delta(3) = (          (            (  56/15._DP + ( -136/ 35._DP + (-1262/105._DP + (   73814/ 2835._DP )*n)*n)*n)*n)*n)*n
    delta(4) = (          (            (              ( 4279/630._DP + (- 332/ 35._DP + ( -399572/14175._DP )*n)*n)*n)*n)*n)*n
    delta(5) = (          (            (              (                ( 4174/315._DP + ( -144838/ 6237._DP )*n)*n)*n)*n)*n)*n
    delta(6) = (          (            (              (                (                (  601676/22275._DP )*n)*n)*n)*n)*n)*n

    AA(0) = 1 +  (          1/ 4._DP        + 1/64._DP * n**2 ) * n**2
    AA(1) =   -3/   2._DP * ( 1._DP - 1/ 8._DP * n**2 - 1/64._DP * n**4 ) * n
    AA(2) =   15/  16._DP * ( 1._DP - 1/ 4._DP * n**2                    ) * n**2
    AA(3) =  -35/  48._DP * ( 1._DP - 5/16._DP * n**2                    ) * n**3
    AA(4) =  315/ 512._DP * n**4
    AA(5) = -693/1280._DP * n**5

  end subroutine set_coef
  !!---------------------------------------------------------------------------------------------------------------------------- !!


  !!---------------------------------------------------------------------------------------------------------------------------- !!
  !>
  !! Hyporbolic arc-tangent
  !<
  !! --
  real(DP) function atanh( x )
    real(DP), intent(in) :: x

    atanh = log( ( 1._DP + x ) / ( 1._DP - x ) ) / 2.0_DP

  end function atanh
  !!---------------------------------------------------------------------------------------------------------------------------- !!

end module m_gk
!! ----------------------------------------------------------------------------------------------------------------------------- !!
