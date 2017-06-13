!! ------------------------------------------------------------------------- !!
!>
!! Snapshot binary
!!
!! @copyright
!!   Copyright 2013-2017 Takuto Maeda. All rights reserved. This project is released under the MIT license.
!<
!! ----
module m_fdsnap

  !! -- Dependency
  use m_std
  use m_daytim
#ifdef _NETCDF
  use netcdf
#endif

  !! -- Declarations
  implicit none
  private
  save

  !! -- Procedures
  public :: fdsnap__open
  public :: fdsnap__readhdr
  public :: fdsnap__writehdr
  public :: fdsnap__checkhdr
  public :: fdsnap__hdr

  !! ----------------------------------------------------------------------- !!
  !>
  !! FDM snapshot file header
  !<
  !!--
  type fdsnap__hdr

    character(8)  :: bintype  !< binary type "STREAMIO" or "UNFORMAT"
    character(8)  :: codetype !< code to use generate snapshot
    integer       :: hdrver   !< version, mostly used for endian check
    character(80) :: title    !< exe title
    integer       :: exedate
    character(2)  :: coordinate !< 'xy', 'xz', 'yz', 'fs', 'ob' etc
    character(2)  :: datatype   !< 'ps', 'v3' ...
    integer       :: ns1, ns2   !< array size
    real(SP)      :: beg1, beg2 !< beggining for each dimensions
    real(SP)      :: ds1, ds2   !< grid width
    real(SP)      :: dt
    integer       :: na1, na2   !< absoption area size
    integer       :: nmed       !< number of medium layers
    integer       :: nsnp
    real(SP)      :: clon, clat, phi
    real(SP)      :: L, w, C0  !< optional: not used

  end type fdsnap__hdr
  !! ----------------------------------------------------------------------- !!


contains

  !! ----------------------------------------------------------------------- !!
  !>
  !! Open fdsnap file with appripriate binary format
  !!
  !! First, it assume that file is Fortran2003 Stream I/O.
  !! If file header indicate it is not stream file, it re-open the same file
  !! as unfomatted binary.
  !<
  !! --
  subroutine fdsnap__open(fname, io, is_exist, snp_type)

    character(*), intent(in)  :: fname
    integer,      intent(out) :: io
    logical,      intent(out) :: is_exist
    character(6), intent(out) :: snp_type
    integer :: ierr
    !! ----

    !! first, try open as netcdf file
#ifdef _NETCDF
    ierr = nf90_open(fname, NF90_NOWRITE, io)
    if(ierr /= NF90_NOERR) then
      snp_type = 'native'
      call native_file_open(fname, io, is_exist)
      return
    else
      is_exist = .true.
      snp_type = 'netcdf'
      return
    endif
#else
    snp_type = 'native'
    call native_file_open(fname, io, is_exist)
    return
#endif

  end subroutine fdsnap__open
  !! ----------------------------------------------------------------------- !!

  !! ----------------------------------------------------------------------- !!
  !>
  !! read native-formatted snapshot file with existence check
  !<
  !! --
  subroutine native_file_open(fname, io, is_exist)

    character(*), intent(in) :: fname
    integer,      intent(out) :: io
    logical,      intent(out) :: is_exist
    character(8) :: bintype

    call std__getio(io, is_big = .true.)
    inquire(file=trim(fname), exist=is_exist)

    if(.not. is_exist) return

    open(io, file=trim(fname), access='stream', action='read', form='unformatted')
    read(io) bintype

    select case (bintype)

    case ("STREAMIO")
      ! DO NOTHING

    case ("UNFORMAT")
      rewind(io)
      close(io)
      open(io, file=trim(fname), form='unformatted', action='read')

    case default
      write(STDERR,'(A)') 'ERROR [fdsnap__open]: unknown binary type '// bintype
      write(STDERR,'(A)') 'ERROR [fdsnap__open]: file is closed'
      close(io)

    end select

    rewind(io)

  end subroutine native_file_open
  !! ----------------------------------------------------------------------- !!

  !! ----------------------------------------------------------------------- !!
  !>
  !! Write header for terminal
  !!
  !<
  !! --
  subroutine fdsnap__checkhdr(io, hdr)

    integer, intent(in) :: io
    type(fdsnap__hdr), intent(in) :: hdr
    integer :: yr, mo, dy, hr, mi, sc

    !! ----

    call daytim__localtime(hdr%exedate, yr, mo, dy, hr, mi, sc)

    write(io, '(A,A)'                         ) "[binary type]   : ", hdr % bintype
    write(io, '(A,A)'                         ) "[code type]     : ", hdr % codetype
    write(io, '(A,I10)'                       ) "[header version]: ", hdr % hdrver
    write(io, '(A,A)'                         ) "[title]         : ", trim(hdr % title)
    write(io, '(A,I10)'                       ) "[date generated]: ", hdr % exedate
    write(io, '(A,I4,A,I2.2,A,I2.2,A,I2.2,A,I2.2,A,I2.2)') "                  ", yr,"-",mo,"-",dy,"T",hr,"-",mi,"-",sc
    write(io, '(A,A)'                         ) "[coordinate]    : ", hdr % coordinate
    write(io, '(A,A)'                         ) "[data type]     : ", hdr % datatype
    write(io, '(A,I10)'                       ) "[ns1]           : ", hdr % ns1
    write(io, '(A,I10)'                       ) "[ns2]           : ", hdr % ns2
    write(io, '(A,F15.5)'                     ) "[beg1]          : ", hdr % beg1
    write(io, '(A,F15.5)'                     ) "[beg2]          : ", hdr % beg2
    write(io, '(A,F15.5)'                     ) "[ds1]           : ", hdr % ds1
    write(io, '(A,F15.5)'                     ) "[ds2]           : ", hdr % ds2
    write(io, '(A,F15.5)'                     ) "[dt]            : ", hdr % dt
    write(io, '(A,I10)'                       ) "[na1]           : ", hdr % na1
    write(io, '(A,I10)'                       ) "[na2]           : ", hdr % na2
    write(io, '(A,I10)'                       ) "[nmed]          : ", hdr % nmed
    write(io, '(A,I10)'                       ) "[nsnp]          : ", hdr % nsnp
    write(io, '(A,F15.5)'                     ) "[clon]          : ", hdr % clon
    write(io, '(A,F15.5)'                     ) "[clat]          : ", hdr % clat
    write(io, '(A,F15.5)'                     ) "[L]             : ", hdr % L
    write(io, '(A,F15.5)'                     ) "[w]             : ", hdr % w
    write(io, '(A,F15.5)'                     ) "[C0]            : ", hdr % C0

  end subroutine fdsnap__checkhdr
  !! ----------------------------------------------------------------------- !!


  !! ----------------------------------------------------------------------- !!
  !>
  !! Write header as native snapshot binary
  !!
  !<
  !! --
  subroutine fdsnap__writehdr(io, hdr)

    integer, intent(in) :: io
    type(fdsnap__hdr), intent(in) :: hdr

    !! ----


    write(io) hdr % bintype
    write(io) hdr % codetype
    write(io) hdr % hdrver
    write(io) hdr % title
    write(io) hdr % exedate
    write(io) hdr % coordinate
    write(io) hdr % datatype
    write(io) hdr % ns1
    write(io) hdr % ns2
    write(io) hdr % beg1
    write(io) hdr % beg2
    write(io) hdr % ds1
    write(io) hdr % ds2
    write(io) hdr % dt
    write(io) hdr % na1
    write(io) hdr % na2
    write(io) hdr % nmed
    write(io) hdr % nsnp
    write(io) hdr % clon
    write(io) hdr % clat
    if (hdr%hdrver == 5) write(io) hdr % phi

    if(hdr%hdrver == 3) then
      write(io) hdr % L
      write(io) hdr % C0
    else if (hdr%hdrver == 4 .or. hdr%hdrver==5) then
      write(io) hdr % L
      write(io) hdr % C0
      write(io) hdr % w
    end if

  end subroutine fdsnap__writehdr
  !! ----------------------------------------------------------------------- !!


  !! ----------------------------------------------------------------------- !!
  !>
  !! Read header part with endian check
  !<
  !! --
  subroutine fdsnap__readhdr(fname, io, snp_type, hdr)
    character(*),      intent(in)  :: fname
    integer,           intent(in)  :: io
    character(6),      intent(in)  :: snp_type
    type(fdsnap__hdr), intent(out) :: hdr

    !! ----

#ifdef _NETCDF
    if(snp_type == 'netcdf') then
      hdr%bintype = 'netcdf'
      call nc_chk(nf90_get_att(io, NF90_GLOBAL, 'codetype',   hdr%codetype))
      call nc_chk(nf90_get_att(io, NF90_GLOBAL, 'hdrver',     hdr%hdrver))
      call nc_chk(nf90_get_att(io, NF90_GLOBAL, 'title',      hdr%title))
      call nc_chk(nf90_get_att(io, NF90_GLOBAL, 'exedate',    hdr%exedate))
      call nc_chk(nf90_get_att(io, NF90_GLOBAL, 'coordinate', hdr%coordinate))
      call nc_chk(nf90_get_att(io, NF90_GLOBAL, 'datatype',   hdr%datatype))
      call nc_chk(nf90_get_att(io, NF90_GLOBAL, 'ns1',        hdr%ns1))
      call nc_chk(nf90_get_att(io, NF90_GLOBAL, 'ns2',        hdr%ns2))
      call nc_chk(nf90_get_att(io, NF90_GLOBAL, 'beg1',       hdr%beg1))
      call nc_chk(nf90_get_att(io, NF90_GLOBAL, 'beg2',       hdr%beg2))
      call nc_chk(nf90_get_att(io, NF90_GLOBAL, 'ds1',        hdr%ds1))
      call nc_chk(nf90_get_att(io, NF90_GLOBAL, 'ds2',        hdr%ds2))
      call nc_chk(nf90_get_att(io, NF90_GLOBAL, 'dt',         hdr%dt))
      call nc_chk(nf90_get_att(io, NF90_GLOBAL, 'na1',        hdr%na1))
      call nc_chk(nf90_get_att(io, NF90_GLOBAL, 'na2',        hdr%na2))
      call nc_chk(nf90_get_att(io, NF90_GLOBAL, 'nmed',       hdr%nmed))
      call nc_chk(nf90_get_att(io, NF90_GLOBAL, 'nsnp',       hdr%nsnp))
      call nc_chk(nf90_get_att(io, NF90_GLOBAL, 'clon',       hdr%clon))
      call nc_chk(nf90_get_att(io, NF90_GLOBAL, 'clat',       hdr%clat))
      call nc_chk(nf90_get_att(io, NF90_GLOBAL, 'phi',        hdr%phi))
    endif
#endif
    if(snp_type == 'native') then
      read(io) hdr % bintype
      read(io) hdr % codetype
      read(io) hdr % hdrver
      call endian_check(io,hdr%hdrver)
      read(io) hdr % title
      read(io) hdr % exedate
      read(io) hdr % coordinate
      read(io) hdr % datatype
      read(io) hdr % ns1
      read(io) hdr % ns2
      read(io) hdr % beg1
      read(io) hdr % beg2
      read(io) hdr % ds1
      read(io) hdr % ds2
      read(io) hdr % dt
      read(io) hdr % na1
      read(io) hdr % na2
      read(io) hdr % nmed
      read(io) hdr % nsnp
      read(io) hdr % clon
      read(io) hdr % clat
      if(hdr%hdrver == 5) read (io) hdr % phi
      if(hdr%hdrver == 3) then
        read(io) hdr % L
        read(io) hdr % C0
      else if (hdr%hdrver == 4 .or. hdr%hdrver==5) then
        read(io) hdr % L
        read(io) hdr % C0
        read(io) hdr % w
      end if
    endif

  contains

    !! ----
    subroutine endian_check(io, hdrver)

      integer, intent(in) :: io
      integer, intent(in) :: hdrver

      if(0 <= hdrver .and. hdrver <= 140101) then
        ! ok
      else
        write(STDERR,'(A,I5,A)') "ERROR [fdsnap__read]: the file " // trim(fname) // " (",io,") is generated in different endian."
        write(STDERR,'(A)'   ) "ERROR [fdsnap__read]: stop reading. close file."
        close(io)
        return
      end if
    end subroutine endian_check
    !! ----

  end subroutine fdsnap__readhdr
  !! ----------------------------------------------------------------------- !!
  !! ------------------------------------------------------------------------ !!
  !>
  !! An internal subroutine to check error in netcdf function calls
  !<
  !! --
  subroutine nc_chk( ierr )

    integer, intent(in) :: ierr
    !! ----

#ifdef _NETCDF
    if( ierr /= NF90_NOERR )  write(STDERR,*) NF90_STRERROR( ierr )
#endif

  end subroutine nc_chk
  !! ------------------------------------------------------------------------ !!

end module m_fdsnap
!! ------------------------------------------------------------------------- !!
