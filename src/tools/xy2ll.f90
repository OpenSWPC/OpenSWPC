!! ------------------------------------------------------------------------- !!
!>
!! Inverse projection from a local Cartesian coordinate to longitude & latitude
!! with given center location
!!
!! Usage:
!!   xy2ll.x x y  lon0 lat0 [phi]
!!     x y : local coordinate (in km)
!!     lon0, lat0: coordinate center (in degrees)
!!     phi (optional): coordinate rotaiton angle. Default is phi=0
!!   Calculated longitude and latitude (in degrees) will be printed in output_unit
!!
!! Copyright 2013-2024 Takuto Maeda. All rights reserved. This project is released under the MIT license.
!<
!! --
#include "../shared/m_debug.h"
program xy2ll

    use iso_fortran_env, only: error_unit, output_unit
    use m_std
    use m_geomap
    use m_system
    use m_getopt
    use m_version
    use m_debug
    implicit none
  !! --
    real(SP) :: lat0, lon0, lat, lon, x, y, phi
    integer(SP) :: narg
    logical :: is_opt1, is_opt2
  !! ----

    call getopt('v', is_opt1)
    call getopt('-version', is_opt2)
    if (is_opt1 .or. is_opt2) call version__display('xy2ll')

    narg = system__iargc()
    if (narg /= 4 .and. narg /= 5) call usage_abort

    call system__getarg(1, x)
    call system__getarg(2, y)
    call system__getarg(3, lon0)
    call system__getarg(4, lat0)
    if (narg == 5) then
        call system__getarg(5, phi)
    else
        phi = 0
    end if

    call assert(-360. <= lon0 .and. lon0 <= 360)
    call assert(-90. <= lat0 .and. lat0 <= 90)

    call geomap__c2g(x, y, lon0, lat0, phi, lon, lat)
    write (output_unit, '(2ES20.10)') lon, lat

contains

    subroutine usage_abort

        write (error_unit, '(A)') ' USAGE: xy2ll.x <x> <y> <lon0> <lat0> [phi]'
        write (error_unit, '(A)') '        <x>    : x-coordinate [km]'
        write (error_unit, '(A)') '        <y>    : y-coordinate [km]'
        write (error_unit, '(A)') '        <lon0> : origin loniitude [deg.]'
        write (error_unit, '(A)') '        <lat0> : origin latitude [deg.]'
        write (error_unit, '(A)') '        [phi]  : rotation angle [deg.] (optional)'
        stop

    end subroutine usage_abort

end program xy2ll
!! ------------------------------------------------------------------------- !!
