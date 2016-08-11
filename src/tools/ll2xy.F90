!! ------------------------------------------------------------------------- !!
!>
!! Projection form longitude & latitude to local Cartesian Coordinate
!! with given center location
!!
!! Usage:
!!   ll2xy.x lon lat lon0 lat0 [phi]
!!     lon lat: location to be projected (in degrees)
!!     lon0, lat0: coordinate center (in degrees)
!!     phi (optional): coordinate rotaiton angle. Default is phi=0
!!   Calculated x and y (in km) coordinate location will be printed in STDOUT
!!
!! @copyright
!!   Copyright 2013-2016 Takuto Maeda. All rights reserved.
!!   This project is released under the MIT license.
!<
!! --
#include "m_debug.h"
program ll2xy
  
  use m_std
  use m_geomap
  use m_system
  use m_debug
  implicit none
  !! --
  real(SP) :: lat0, lon0, lat, lon, x, y, phi
  integer(SP) :: narg
  !! ----
  
  narg = system__iargc()
  if( narg /= 4 .and. narg /= 5 ) call usage_abort
  
  call system__getarg(1, lon )
  call system__getarg(2, lat )
  call system__getarg(3, lon0 )
  call system__getarg(4, lat0  )
  if( narg == 5 ) then
    call system__getarg(5, phi)
  else
    phi = 0
  end if
  
  call assert( -360. <= lon  .and. lon  <= 360 )
  call assert( -90.  <= lat  .and. lat  <= 90  )
  call assert( -360. <= lon0 .and. lon0 <= 360 )
  call assert( -90.  <= lat0 .and. lat0 <= 90  )
  
  call geomap__g2c( lon, lat, lon0, lat0, phi, x, y)
  write(STDOUT,'(2ES20.10)' ) x, y
  
contains
  
  subroutine usage_abort
    
    write(STDERR,'(A)') ' USAGE: ll2xy.x <lon> <lat> <lon0> <lat0> [phi]'
    write(STDERR,'(A)') '        <lon>  : lonitude [deg.]'
    write(STDERR,'(A)') '        <lat>  : latitude  [deg.]'
    write(STDERR,'(A)') '        <lon0> : origin loniitude [deg.]'
    write(STDERR,'(A)') '        <lat0> : origin latitude [deg.]'
    write(STDERR,'(A)') '        [phi]  : rotation angle [deg.] (optional)'
    stop
  end subroutine usage_abort
  
end program ll2xy
!! ------------------------------------------------------------------------- !!
