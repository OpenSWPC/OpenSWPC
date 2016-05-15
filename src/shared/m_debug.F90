!! ------------------------------------------------------------------------------------------------------------------------------ !!
!>
!! debug routines.
!! Use with preprocessor macro (#include "m_pdebug.h" at the top of the soruce code) will give richer information.
!!
!! @par usage
!! - call debug(var):     show variable var.
!! - call assert( cond ): abort if cond = .false. . cond must be logical value or condition.
!! - call info( msg ):    write message msg to STDERR.
!!
!! @copyright
!!   Copyright 2013-2016 Takuto Maeda. All rights reserved. This project is released under the MIT license.
!<
!! --
module m_debug

  use m_std
  implicit none
  private

  public :: debug
  public :: info
  public :: assert
  public :: debug__macro
  public :: info__macro
  public :: assert__macro
  public :: debug__void

  !! debug through preprocessor macro
  interface debug__macro
     module procedure debug_c,  debug_i,  debug_r,  debug_d,  debug_l,  debug__void
     module procedure debug_c1, debug_i1, debug_r1, debug_d1, debug_l1
  end interface

  !! regular debug
  interface debug
     module procedure debug_c0, debug_i0, debug_r0, debug_d0, debug_l0
  end interface

  !! assersion with priprosessor macro
  interface assert__macro
     module procedure  assert_1, assert_2
  end interface

  !! regular assertion
  interface assert
     module procedure  assert_0
  end interface



contains

  !! ---------------------------------------------------------------------------------------------------------------------------- !!
  !>
  !! Do Nothing: dummy
  !<
  !! --
  subroutine debug__void()

    return

  end subroutine debug__void
  !! ---------------------------------------------------------------------------------------------------------------------------- !!

  !! ---------------------------------------------------------------------------------------------------------------------------- !!
  !>
  !! Debug output to STDERR.
  !<
  !! --
  subroutine debug_c( var, fname, nline )

    character(*), intent(in) :: var
    character(*), intent(in) :: fname
    integer,      intent(in) :: nline
    !! --
    character(5) :: cline
    !! ----

    write(cline,'(I5)') nline

    write(STDERR,'(A)') '[debug] '//fname//' ('//trim(adjustl(cline))//'):  '//trim(adjustl(var))

  end subroutine debug_c
  !! ---------------------------------------------------------------------------------------------------------------------------- !!

  !! ---------------------------------------------------------------------------------------------------------------------------- !!
  !>
  !! Debug output to STDERR.
  !<
  !! --
  subroutine debug_c1( var, varname, fname, nline )

    character(*), intent(in) :: var
    character(*), intent(in) :: varname
    character(*), intent(in) :: fname
    integer,      intent(in) :: nline
    !! --
    character(5) :: cline
    !! ----

    write(cline,'(I5)') nline

    write(STDERR,'(A)') '[debug] '//fname//' ('//trim(adjustl(cline))//'): '//trim(adjustl(varname)) //' = '//trim(adjustl(var))

  end subroutine debug_c1
  !! ---------------------------------------------------------------------------------------------------------------------------- !!

  !! ---------------------------------------------------------------------------------------------------------------------------- !!
  !>
  !! Debug output to STDERR.
  !<
  !! --
  subroutine debug_c0( var )

    character(*), intent(in) :: var
    !! ----

    write(STDERR,'(A)') '[debug] '//trim(adjustl(var))

  end subroutine debug_c0
  !! ---------------------------------------------------------------------------------------------------------------------------- !!

  !! ---------------------------------------------------------------------------------------------------------------------------- !!
  !>
  !! Debug output to STDERR.
  !<
  !! --
  subroutine debug_i( var, fname, nline )

    integer,      intent(in) :: var
    character(*), intent(in) :: fname
    integer,      intent(in) :: nline
    !! --
    character(10) :: cvar
    !! ----

    write(cvar,'(I10)') var
    call debug_c( cvar, fname, nline )

  end subroutine debug_i
  !! ---------------------------------------------------------------------------------------------------------------------------- !!

  !! ---------------------------------------------------------------------------------------------------------------------------- !!
  !>
  !! Debug output to STDERR.
  !<
  !! --
  subroutine debug_i1( var, varname,  fname, nline )

    integer,      intent(in) :: var
    character(*), intent(in) :: varname
    character(*), intent(in) :: fname
    integer,      intent(in) :: nline
    !! --
    character(10) :: cvar
    !! ----

    write(cvar,'(I10)') var
    call debug_c1( cvar, varname, fname, nline )

  end subroutine debug_i1
  !! ---------------------------------------------------------------------------------------------------------------------------- !!

  !! ---------------------------------------------------------------------------------------------------------------------------- !!
  !>
  !! Debug output to STDERR.
  !<
  !! --
  subroutine debug_i0( var )

    integer,      intent(in) :: var
    !! --
    character(10) :: cvar
    !! ----

    write(cvar,'(I10)') var
    call debug_c0( cvar )

  end subroutine debug_i0
  !! ---------------------------------------------------------------------------------------------------------------------------- !!

  !! ---------------------------------------------------------------------------------------------------------------------------- !!
  !>
  !! Debug output to STDERR.
  !<
  !! --
  subroutine debug_r( var, fname, nline )

    real,     intent(in) :: var
    character(*), intent(in) :: fname
    integer,      intent(in) :: nline
    !! --
    character(15) :: cvar
    !! ----

    if( abs(var) < 10000.) then
       write(cvar,'(F15.5)') var
    else
       write(cvar,'(ES15.5)') var
    end if
    call debug_c( cvar, fname, nline )

  end subroutine debug_r
  !! ---------------------------------------------------------------------------------------------------------------------------- !!

  !! ---------------------------------------------------------------------------------------------------------------------------- !!
  !>
  !! Debug output to STDERR.
  !<
  !! --
  subroutine debug_r1( var, varname, fname, nline )

    real,     intent(in) :: var
    character(*), intent(in) :: varname
    character(*), intent(in) :: fname
    integer,      intent(in) :: nline
    !! --
    character(15) :: cvar
    !! ----

    if( abs(var) < 10000.) then
       write(cvar,'(F15.5)') var
    else
       write(cvar,'(ES15.5)') var
    end if
    call debug_c1( cvar, varname, fname, nline )

  end subroutine debug_r1
  !! ---------------------------------------------------------------------------------------------------------------------------- !!

  !! ---------------------------------------------------------------------------------------------------------------------------- !!
  !>
  !! Debug output to STDERR.
  !<
  !! --
  subroutine debug_r0( var )

    real,     intent(in) :: var
    !! --
    character(15) :: cvar
    !! ----

    if( abs(var) < 10000.) then
       write(cvar,'(F15.5)') var
    else
       write(cvar,'(ES15.5)') var
    end if
    call debug_c0( cvar )

  end subroutine debug_r0
  !! ---------------------------------------------------------------------------------------------------------------------------- !!

  !! ---------------------------------------------------------------------------------------------------------------------------- !!
  !>
  !! Debug output to STDERR.
  !<
  !! --
  subroutine debug_d( var, fname, nline )

    real(DP),     intent(in) :: var
    character(*), intent(in) :: fname
    integer,      intent(in) :: nline
    !! --
    character(15) :: cvar
    !! ----

    if( abs(var) < 10000.) then
       write(cvar,'(F15.5)') var
    else
       write(cvar,'(ES15.5)') var
    end if
    call debug_c( cvar, fname, nline )

  end subroutine debug_d
  !! ---------------------------------------------------------------------------------------------------------------------------- !!
  !! ---------------------------------------------------------------------------------------------------------------------------- !!
  !>
  !! Debug output to STDERR.
  !<
  !! --
  subroutine debug_d1( var, varname, fname, nline )

    real(DP),     intent(in) :: var
    character(*), intent(in) :: varname
    character(*), intent(in) :: fname
    integer,      intent(in) :: nline
    !! --
    character(15) :: cvar
    !! ----

    if( abs(var) < 10000.) then
       write(cvar,'(F15.5)') var
    else
       write(cvar,'(ES15.5)') var
    end if
    call debug_c1( cvar, varname, fname, nline )

  end subroutine debug_d1
  !! ---------------------------------------------------------------------------------------------------------------------------- !!

  !! ---------------------------------------------------------------------------------------------------------------------------- !!
  !>
  !! Debug output to STDERR.
  !<
  !! --
  subroutine debug_d0( var )

    real(DP),     intent(in) :: var
    !! --
    character(15) :: cvar
    !! ----

    if( abs(var) < 10000.) then
       write(cvar,'(F15.5)') var
    else
       write(cvar,'(ES15.5)') var
    end if
    call debug_c0( cvar )

  end subroutine debug_d0
  !! ---------------------------------------------------------------------------------------------------------------------------- !!

  !! ---------------------------------------------------------------------------------------------------------------------------- !!
  !>
  !! Debug output to STDERR.
  !<
  !! --
  subroutine debug_l( var, fname, nline )

    logical,      intent(in) :: var
    character(*), intent(in) :: fname
    integer,      intent(in) :: nline
    !! --
    character(15) :: cvar
    !! ----
    if( var ) then
       cvar = '.true.'
    else
       cvar = '.false.'
    end if
    call debug_c( cvar, fname, nline )

  end subroutine debug_l
  !! ---------------------------------------------------------------------------------------------------------------------------- !!

  !! ---------------------------------------------------------------------------------------------------------------------------- !!
  !>
  !! Debug output to STDERR.
  !<
  !! --
  subroutine debug_l1( var, varname, fname, nline )

    logical,      intent(in) :: var
    character(*), intent(in) :: varname
    character(*), intent(in) :: fname
    integer,      intent(in) :: nline
    !! --
    character(15) :: cvar
    !! ----
    if( var ) then
       cvar = '.true.'
    else
       cvar = '.false.'
    end if
    call debug_c1( cvar, varname, fname, nline )

  end subroutine debug_l1
  !! ---------------------------------------------------------------------------------------------------------------------------- !!

  !! ---------------------------------------------------------------------------------------------------------------------------- !!
  !>
  !! Debug output to STDERR.
  !<
  !! --
  subroutine debug_l0( var )

    logical,      intent(in) :: var
    !! --
    character(15) :: cvar
    !! ----
    if( var ) then
       cvar = '.true.'
    else
       cvar = '.false.'
    end if
    call debug_c0( cvar )

  end subroutine debug_l0
  !! ---------------------------------------------------------------------------------------------------------------------------- !!


  !! --------------------------------------------------------------------------------------------------------------------------- !!
  !>
  !! A simple assetion
  !<
  !! ----
  subroutine assert_0( cond )

    logical, intent(in) :: cond

    !! ----

    if( .not. cond ) then
       write(STDERR,'(A)') '[assert] failed'
       stop
    end if

  end subroutine assert_0
  !! --------------------------------------------------------------------------------------------------------------------------- !!

  !! --------------------------------------------------------------------------------------------------------------------------- !!
  !>
  !! Assertion with filename and line
  !<
  !! ----
  subroutine assert_1( cond, fname, nline )

    logical,      intent(in) :: cond
    character(*), intent(in) :: fname
    integer,      intent(in) :: nline
    !! --
    character(5) :: cl
    !! ----

    if( .not. cond ) then
       write(cl,'(I5)') nline

       write(STDERR,'(A)') '[assert] failed at ' // trim(adjustl(fname)) //'('//trim(adjustl(cl))//')'
       stop
    end if

  end subroutine assert_1
  !! --------------------------------------------------------------------------------------------------------------------------- !!

  !! --------------------------------------------------------------------------------------------------------------------------- !!
  !>
  !! Assertion with condition name, filename and line
  !<
  !! ----
  subroutine assert_2( cond, cname, fname, nline )

    logical,      intent(in) :: cond
    character(*), intent(in) :: cname
    character(*), intent(in) :: fname
    integer,      intent(in) :: nline
    !! ----

    if( .not. cond ) then
       write(STDERR,'(A,I0,A)') '[assert] failed at ' // trim(adjustl(fname)) //'(', nline, '): ' // trim(adjustl(cname))
       stop
    end if

  end subroutine assert_2
  !! --------------------------------------------------------------------------------------------------------------------------- !!


  !! --------------------------------------------------------------------------------------------------------------------------- !!
  subroutine info( msg )

    character(*), intent(in) :: msg

    write(STDERR,*) '[info] ' // trim(adjustl(msg))

  end subroutine info
  !! --------------------------------------------------------------------------------------------------------------------------- !!

  !! --------------------------------------------------------------------------------------------------------------------------- !!
  !>
  !! information message with filename and line number
  !<
  !! --
  subroutine info__macro( msg, fname, nline )

    character(*), intent(in) :: msg
    character(*), intent(in) :: fname
    integer,      intent(in) :: nline
    !! ----

    write(STDERR,'(A,I0,A)') '[info] '// trim(adjustl(fname)) // '(', nline, '): ' // trim(adjustl(msg))

  end subroutine info__macro
  !! --------------------------------------------------------------------------------------------------------------------------- !!


end module m_debug
!! ------------------------------------------------------------------------------------------------------------------------------ !!
