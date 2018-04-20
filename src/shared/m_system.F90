!! ----------------------------------------------------------------------------------------------------------------------------- !!
!>
!! Linux/Mac system routines, for processig command line argument, environment variables, and system call
!!
!! @copyright
!!   Copyright 2013-2018 Takuto Maeda. All rights reserved. This project is released under the MIT license.
!<
!! --
module m_system

  !! -- Dependency
  use m_std

  !! -- Declarations
  implicit none
  private

  public :: system__getarg
  public :: system__getenv
  public :: system__call
  public :: system__iargc
  public :: system__expenv

  !! --------------------------------------------------------------------------------------------------------------------------- !!
  !>
  !! Obtain command line arguments for character, integer, single and double precision data
  !!
  !! It uses Fortran2003 statement
  !<
  !! --
  interface system__getarg

    module procedure getarg_a, getarg_i, getarg_f, getarg_d

  end interface system__getarg
  !! --------------------------------------------------------------------------------------------------------------------------- !!

contains

  !! --------------------------------------------------------------------------------------------------------------------------- !!
  !>
  !! Do system-call.
  !!
  !! @note It uses out-of-standard routines, but it works with most of modern fortran compilers.
  !!
  !<
  !! --
  subroutine system__call (cmd)

    character(*), intent(in) :: cmd
    !! ----

    call system( cmd )

  end subroutine system__call
  !! --------------------------------------------------------------------------------------------------------------------------- !!

  !! --------------------------------------------------------------------------------------------------------------------------- !!
  !>
  !! Returns a number of arguments. Fortran2003 wrapper function
  !<
  !! --
  integer function system__iargc()

    system__iargc  = command_argument_count()

  end function system__iargc
  !! --------------------------------------------------------------------------------------------------------------------------- !!

  !! --------------------------------------------------------------------------------------------------------------------------- !!
  !>
  !! Obtain environmental variable "name".
  !<
  !! --
  subroutine system__getenv( name, value )

    !! -- Arguments
    character(*), intent(in)  :: name
    character(*), intent(out) :: value

    !! ----

    call get_environment_variable( name, value )

  end subroutine system__getenv
  !! --------------------------------------------------------------------------------------------------------------------------- !!

  !! --------------------------------------------------------------------------------------------------------------------------- !!
  !>
  !! get i-th command line argument, Fortran2003 wrapper subroutine
  !<
  !! --
  subroutine getarg_a (i, arg)

    !! -- Arguments
    integer,      intent(in)  :: i   ! order of the arguments
    character(*), intent(out) :: arg ! argument

    call get_command_argument( i, arg )

  end subroutine getarg_a
  !! --------------------------------------------------------------------------------------------------------------------------- !!

  !! --------------------------------------------------------------------------------------------------------------------------- !!
  subroutine getarg_i (i, arg)

    !! -- Arguments
    integer, intent(in) :: i
    integer, intent(out) :: arg

    character(256) :: carg
    !! ----

    call getarg_a( i, carg )
    read(carg,*) arg

  end subroutine getarg_i
  !! --------------------------------------------------------------------------------------------------------------------------- !!

  !! --------------------------------------------------------------------------------------------------------------------------- !!
  subroutine getarg_f (i, arg)
    !! -- Arguments
    integer,  intent(in) :: i
    real, intent(out) :: arg

    character(256) :: carg
    !! ----

    call getarg_a( i, carg )
    read(carg,*) arg

  end subroutine getarg_f
  !! --------------------------------------------------------------------------------------------------------------------------- !!

  !! --------------------------------------------------------------------------------------------------------------------------- !!
  subroutine getarg_d (i, arg)

    !! -- Arguments
    integer,  intent(in) :: i
    real(DP), intent(out) :: arg

    character(256) :: carg
    !! ----

    call getarg_a( i, carg )
    read(carg,*) arg

  end subroutine getarg_d
  !! --------------------------------------------------------------------------------------------------------------------------- !!

  !! --------------------------------------------------------------------------------------------------------------------------- !!
  !>
  !! Expand shell environmental variables wrapped in ${...}
  !<
  !! --
  subroutine system__expenv( str )
    character(*), intent(inout) :: str
    character(256) :: str2
    integer :: ikey1, ikey2
    integer :: iptr
    character(256) :: str3

    iptr = 1
    str2 = ''
    do
      ikey1 = scan( str(iptr:), "${" ) + iptr - 1
      if( ikey1==iptr-1 ) exit

      ikey2 = scan( str(iptr:), "}" ) + iptr -1
      str2=trim(str2) // str(iptr:ikey1-1)
      call system__getenv( str(ikey1+2:ikey2-1), str3 )
      str2 = trim(str2) // trim(str3)
      iptr = ikey2+1

    end do
    str2 = trim(str2) // trim(str(iptr:))

    str = trim(str2)

  end subroutine system__expenv
  !! --------------------------------------------------------------------------------------------------------------------------- !!


end module m_system
!! ----------------------------------------------------------------------------------------------------------------------------- !!
