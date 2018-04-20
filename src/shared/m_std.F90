!! ----------------------------------------------------------------------------------------------------------------------------- !!
!>
!! Definition of standard constants, in/out io numbers, precision constants
!!
!! @copyright
!!   Copyright 2013-2018 Takuto Maeda. All rights reserved. This project is released under the MIT license.
!<
!! --
module m_std

  !! -- Declarations
  implicit none
  private


  !! -- Precision Constant
  integer,     parameter, public :: DP = selected_real_kind(13) !< Double Precision
  integer,     parameter, public :: SP = selected_real_kind(5)  !< Single Precision

  !! -- In/Out
  integer,     parameter, public :: STDERR  = 0 !< Standard Error
  integer,     parameter, public :: STDOUT  = 6 !< Standard Output
  integer,     parameter, public :: STDIN   = 5 !< Standard Input

  !! -- Physics/Math Constants
  real(DP),    parameter, public :: PI      = 3.141592653589793238462643383_DP
  real(DP),    parameter, public :: R_EARTH = 6371.0_DP
  complex(DP), parameter, public :: EI      = (0.0_DP,1.0_DP)

  !! -- Precision-specific constant numbers with specific precision
  real(DP),    parameter, public :: PI_D      = 3.141592653589793238462643383_DP
  real(DP),    parameter, public :: R_EARTH_D = 6371.0_DP
  real(SP),    parameter, public :: PI_S      = 3.141592653589793238462643383
  real(SP),    parameter, public :: R_EARTH_S = 6371.0

  !! -- Public Subroutines
  public :: std__getio
  public :: std__genfname
  public :: std__countline
  public :: std__deg2rad
  public :: std__rad2deg


  interface std__deg2rad
    module procedure d2r_s, d2r_d
  end interface std__deg2rad
  interface std__rad2deg
    module procedure r2d_s, r2d_d
  end interface std__rad2deg


  !! -- Private variables
  integer, parameter, private :: io0        = 10     !< initial number for file I/O
  integer, parameter, private :: IOBIG0     = 900    !< Binary data with io >= IOBIG0 are assumed in BIG_ENDIAN data

  !----

  !! --------------------------------------------------------------------------------------------------------------------------- !!
  !>
  !! Return the unusedd unit number measured from io0 constant
  !!
  !! @par Example
  !! - call std__getio( io )
  !! - call std__getio( io, .true. ) !< Big-endian output
  !! - call std__getio( io, 300 ) !< search io number from 300
  !!
  !<
  !--
  interface std__getio
    module procedure getio_0, &  !< default/number specific
        getio_big   !< big_endian

  end interface std__getio
  !! --------------------------------------------------------------------------------------------------------------------------- !!

contains


  !! --------------------------------------------------------------------------------------------------------------------------- !!
  function d2r_d ( deg )

    real(DP) :: d2r_d
    real(DP), intent(in) :: deg

    d2r_d = PI_D / 180.0_DP * deg

  end function d2r_d
  !! --------------------------------------------------------------------------------------------------------------------------- !!
  !! --------------------------------------------------------------------------------------------------------------------------- !!
  function d2r_s ( deg )
    real(SP) :: d2r_s
    real(SP), intent(in) :: deg

    d2r_s = PI_S / 180.0_SP * deg

  end function d2r_s
  !! --------------------------------------------------------------------------------------------------------------------------- !!
  !! --------------------------------------------------------------------------------------------------------------------------- !!
  function r2d_d( rad )
    real(DP) :: r2d_d
    real(DP), intent(in) :: rad

    r2d_d = 180.0_DP / PI_D * rad

  end function r2d_d
  !! --------------------------------------------------------------------------------------------------------------------------- !!
  !! --------------------------------------------------------------------------------------------------------------------------- !!
  function r2d_s( rad )
    real(SP) :: r2d_s
    real(SP), intent(in) :: rad

    r2d_s = 180.0_SP / PI_S * rad

  end function r2d_s
  !! --------------------------------------------------------------------------------------------------------------------------- !!


  !! --------------------------------------------------------------------------------------------------------------------------- !!
  !>
  !!  Generate filename "base.????.ext" which isn't exit in current directory.
  !!
  !! @par Example
  !!   call std__genfname( 'foo','dat', fname )
  !!   call std__getio( fp )
  !!   open ( fp, file=fname ) ! fname = 'foo.0000.dat' if there's no file
  !<
  !! --
  subroutine std__genfname( base, ext, fname )

    !! -- Arguments
    character(len=*), intent(in)  :: base  ! base filnemae
    character(len=*), intent(in)  :: ext   ! extention
    character(len=*), intent(out) :: fname ! output

    integer :: i
    character(4) :: ci
    logical :: isExist

    !! ----

    i=0
    do
      write(ci,'(I4.4)') i
      fname = adjustl(trim(base))//'.'//adjustl(trim(ci))&
          //'.'//adjustl(trim(ext))
      inquire(file=trim(fname), exist=isExist)
      if( isExist ) then
        i=i+1
      else
        exit
      end if
    end do

  end subroutine std__genfname
  !! --------------------------------------------------------------------------------------------------------------------------- !!


  !! --------------------------------------------------------------------------------------------------------------------------- !!
  !>
  !! Count line nubmer n included in the file specified by io
  !<
  !! --
  subroutine std__countline( io, n, comment )

    !! -- Arguments
    integer,   intent(in)  :: io
    integer,   intent(out) :: n
    character, intent(in), optional :: comment

    integer :: stat
    character (256) :: line

    !! ----

    n = 0
    rewind(io)
    do
      read( io, '(A256)', iostat=stat) line
      if( stat /= 0 ) exit

      if( present( comment ) ) then
        if( line(1:1) == comment ) cycle
      end if

      if( trim(line) /= '' ) then
        n = n + 1
      end if

    end do
    rewind( io )

  end subroutine std__countline
  !! --------------------------------------------------------------------------------------------------------------------------- !!


  !! --------------------------------------------------------------------------------------------------------------------------- !!
  !>
  !! Return the unused unit number measured from pre-defined io0 constant
  !! If io00 is given, it seeks the unit number from io00
  !<
  !! --
  subroutine getio_0( io, io00 )

    !! -- Arguments
    integer, intent(out) :: io               !< unit number
    integer, intent(in), optional :: io00    !< start number (optional)
    !+
    logical :: isOpen
    !--

    if( present( io00 ) ) then
      io = io00
    else
      io = io0
    end if

    isOpen = .true.

    do
      inquire( io, opened = isOpen )
      if( .not. isOpen ) exit
      io = io + 1
    end do

  end subroutine getio_0
  !! --------------------------------------------------------------------------------------------------------------------------- !!


  !! --------------------------------------------------------------------------------------------------------------------------- !!
  !>
  !!  Return the unusedd unit number measured from io0 constant
  !!
  !! If is_big = .true., it searchs io number from IOBIG0. Use it together with appropriate compiler settings
  !<
  !! --
  subroutine getio_big( io, is_big )

    !! -- Arguments
    integer, intent(out) :: io      ! unit number
    logical, intent(in)  :: is_big  ! BIG_ENDIAN
    !+
    logical :: isOpen
    !--

    if( is_big) then
      io = IOBIG0
    else
      io = io0
    end if

    isOpen = .true.

    do
      inquire( io, opened = isOpen )
      if( .not. isOpen ) exit
      io = io + 1
    end do

  end subroutine getio_big
  !! --------------------------------------------------------------------------------------------------------------------------- !!

end module m_std
!! ----------------------------------------------------------------------------------------------------------------------------- !!
