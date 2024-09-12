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
!!   Calculated x and y (in km) coordinate location will be printed in output_unit
!!
!! @copyright
!!   Copyright 2013-2024 Takuto Maeda. All rights reserved. This project is released under the MIT license.
!<
!! --
#include "../shared/m_debug.h"
program ll2xy

  use iso_fortran_env, only: error_unit, output_unit
  use m_std
  use m_geomap
  use m_system
  use m_getopt
  use m_debug
  use m_version
  implicit none
  !! --
  real(SP) :: lat0, lon0, lat, lon, x, y, phi
  integer(SP) :: narg
  logical :: is_opt1, is_opt2
  !! ----

  call getopt('v', is_opt1)
  call getopt('-version', is_opt2)
  if( is_opt1 .or. is_opt2 ) call version__display('ll2xy')

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
  write(output_unit,'(2ES20.10)' ) x, y

contains

  subroutine usage_abort

    write(error_unit,'(A)') ' USAGE: ll2xy.x <lon> <lat> <lon0> <lat0> [phi]'
    write(error_unit,'(A)') '        <lon>  : lonitude [deg.]'
    write(error_unit,'(A)') '        <lat>  : latitude  [deg.]'
    write(error_unit,'(A)') '        <lon0> : origin loniitude [deg.]'
    write(error_unit,'(A)') '        <lat0> : origin latitude [deg.]'
    write(error_unit,'(A)') '        [phi]  : rotation angle [deg.] (optional)'
    stop
  end subroutine usage_abort

end program ll2xy
!! ------------------------------------------------------------------------- !!
