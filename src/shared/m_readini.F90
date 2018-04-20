!! ----------------------------------------------------------------------------------------------------------------------------- !!
!>
!! Read ini-style parameter file
!!
!! @copyright
!!   Copyright 2013-2018 Takuto Maeda. All rights reserved. This project is released under the MIT license.
!<
!! ----
#include "m_debug.h"
module m_readini

  !! -- Dependency
  use m_std
  use m_debug
  use m_system

  !! -- Declarations
  implicit none
  private
  save

  !! -- Public Procedures
  public :: readini

  !--

  !! --------------------------------------------------------------------------------------------------------------------------- !!
  !>
  !! find keyword 'key' from ini-style file, and returns var
  !!
  !!
  !! @par Usage
  !!   call readini( io, key, var, def )
  !!     - io    (integer)   ! file io number
  !!     - key   (character) ! keyword
  !!     - var (any type)    ! output parameter
  !!     - def (any type)    ! default value
  !<
  !! --
  interface readini

    module procedure readini_d, readini_s, readini_i, readini_c, readini_l

  end interface readini
  !! --------------------------------------------------------------------------------------------------------------------------- !!

contains


  !! --------------------------------------------------------------------------------------------------------------------------- !!
  subroutine readini_c( io, key, var, def )

    !! -- Arguments
    integer,       intent(in)  :: io
    character(*),  intent(in)  :: key
    character(*),  intent(in)  :: def
    character(*),  intent(out) :: var

    character(256) :: keyword
    integer        :: ierr
    character(256) :: cline
    integer        :: keylen
    logical        :: isopen
    !! ----

    !! file status
    inquire( io, OPENED=isopen )

    if( .not. isopen ) then

      write(STDERR,'(A)') 'ERROR [readini]: file not open.'
      var = trim(def)
      return

    end if



    !! initialize file I/O location
    rewind(io)

    keyword = trim(adjustl(key))
    keylen  = len_trim(keyword)


    do

      !! get one line
      read(io,'(A)', iostat=ierr) cline

      !! reach to the last line
      if( ierr /= 0 ) then
        call info( 'key ' // trim(keyword) //' is not found.' )
        call info( '    Use default value '//trim(def)//' instead.' )
        var = trim(def)
        return
      end if

      cline = adjustl(cline)


      !! comment line
      if( cline(1:1) == '#' .or. cline(1:1) == '!' ) cycle


      !! find keyword
      if( cline(1:keylen) == trim(keyword) ) then

        cline = adjustl( cline( keylen+1: ) )

        if( cline(1:1) == '=' ) then
          cline = adjustl(cline(2:))
          read(cline,*) var
          exit
        end if
      end if
    end do

    rewind(io)

    !! expand environmental variable
    call system__expenv( var )

  end subroutine readini_c
  !! --------------------------------------------------------------------------------------------------------------------------- !!

  !! --------------------------------------------------------------------------------------------------------------------------- !!
  subroutine readini_d( io, key, var, def )

    integer,      intent(in)  :: io
    character(*), intent(in)  :: key
    real(DP),     intent(in)  :: def
    real(DP),     intent(out) :: var
    !! --
    character(256) :: avar
    character(256) :: adef
    !! ----

    write(adef,*) def
    call readini_c( io, key, avar, adef )
    read( avar,* ) var

  end subroutine readini_d
  !! --------------------------------------------------------------------------------------------------------------------------- !!

  !! --------------------------------------------------------------------------------------------------------------------------- !!
  subroutine readini_s( io, key, var, def )

    integer,      intent(in)  :: io
    character(*), intent(in)  :: key
    real(SP),     intent(in)  :: def
    real(SP),     intent(out) :: var
    !! -
    character(256) :: avar, adef
    !! ----
    write(adef,*) def
    call readini_c( io, key, avar, adef )
    read(avar,*) var

  end subroutine readini_s
  !! --------------------------------------------------------------------------------------------------------------------------- !!


  !! --------------------------------------------------------------------------------------------------------------------------- !!
  subroutine readini_i( io, key, var, def )
    integer,      intent(in)  :: io
    character(*), intent(in)  :: key
    integer,      intent(in)  :: def
    integer,      intent(out) :: var
    !! --
    character(256) :: avar, adef
    !! ----
    write(adef,*) def
    call readini_c( io, key, avar, adef )
    read(avar,*) var

  end subroutine readini_i
  !! --------------------------------------------------------------------------------------------------------------------------- !!

  !! --------------------------------------------------------------------------------------------------------------------------- !!
  subroutine readini_l( io, key, var, def )

    integer,      intent(in)  :: io
    character(*), intent(in)  :: key
    logical,      intent(in)  :: def
    logical,      intent(out) :: var
    !! --
    character(256) :: avar, adef
    !! ----

    write(adef,*) def
    call readini_c( io, key, avar, adef )
    read(avar,*) var

  end subroutine readini_l
  !! --------------------------------------------------------------------------------------------------------------------------- !!

end module m_readini
!! ----------------------------------------------------------------------------------------------------------------------------- !!
