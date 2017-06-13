!! ----------------------------------------------------------------------------------------------------------------------------- !!
!>
!! Geographical coordinate <-> Cartesian Coordinate
!!
!! @copyright
!!   Copyright 2013-2017 Takuto Maeda. All rights reserved. This project is released under the MIT license.
!<
!! --
module m_geomap

  use m_std
  use m_gk
  implicit none

  save
  public

  !! ----

contains



  !! --------------------------------------------------------------------------------------------------------------------------- !!
  !>
  !! Geomap change from geographical (lat,lon) to cartesian (x,y) with rotation
  !<
  !! --
  subroutine geomap__g2c( lon, lat, lon0, lat0, phi, x, y )

    !! -- arguments
    real(SP), intent(in)  :: lon
    real(SP), intent(in)  :: lat
    real(SP), intent(in)  :: lon0 !< center longitude
    real(SP), intent(in)  :: lat0 !< center latitude
    real(SP), intent(in)  :: phi  !< clockwise azimuth (deg.). for phi=0, x:north and y:east
    real(SP), intent(out) :: x
    real(SP), intent(out) :: y

    real(SP) :: xx, yy
    real(SP) :: phi_r

    !! ----

    phi_r = std__deg2rad( phi )
    call gk__lltoxy( lon, lat, lon0, lat0, xx, yy )

    x =   cos(phi_r) * xx + sin(phi_r) * yy
    y = - sin(phi_r) * xx + cos(phi_r) * yy

  end subroutine geomap__g2c
  !! --------------------------------------------------------------------------------------------------------------------------- !!


  !! --------------------------------------------------------------------------------------------------------------------------- !!
  !>
  !! Geomap change from geographical (lat,lon) to cartesian (x,y) with rotation
  !<
  !! --
  subroutine geomap__c2g( x, y, lon0, lat0, phi, lon, lat )

    !! -- arguments
    real(SP), intent(in) :: x
    real(SP), intent(in) :: y
    real(SP), intent(in)  :: lon0 !< center longitude
    real(SP), intent(in)  :: lat0 !< center latitude
    real(SP), intent(in)  :: phi !< clockwise azimuth (deg.). for phi=0, x:north and y:east
    real(SP), intent(out)  :: lon
    real(SP), intent(out)  :: lat

    real(SP) :: xx, yy
    real(SP) :: phi_r

    !! ----

    phi_r = std__deg2rad( phi )
    xx = cos(phi_r) * x - sin(phi_r) * y
    yy = sin(phi_r) * x + cos(phi_r) * y

    call gk__xytoll( xx, yy, lon0, lat0, lon, lat )

  end subroutine geomap__c2g
  !! --------------------------------------------------------------------------------------------------------------------------- !!

end module m_geomap
!! ----------------------------------------------------------------------------------------------------------------------------- !!
