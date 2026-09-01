module m_fdsnap

    !! Snapshot binary
    !!
    !! Copyright 2013-2026 Takuto Maeda. All rights reserved. This project is released under the MIT license.

    use iso_fortran_env, only: error_unit
    use m_std
    use m_daytim
    use netcdf

    implicit none
    private
    save

    public :: fdsnap__open
    public :: fdsnap__readhdr
    public :: fdsnap__checkhdr
    public :: fdsnap__hdr

    type fdsnap__hdr

        !! FDM snapshot file header

        character(8)  :: bintype  !! binary type "STREAMIO" or "UNFORMAT"
        character(8)  :: codetype !! code to use generate snapshot
        integer       :: hdrver   !! version
        character(80) :: title    !! exe title
        integer       :: exedate  !! date of execution
        character(2)  :: coordinate !! 'xy', 'xz', 'yz', 'fs', 'ob' etc
        character(2)  :: datatype   !! 'ps', 'v3' ...
        integer       :: ns1, ns2   !! array size
        real(SP)      :: beg1, beg2 !! beggining for each dimensions
        real(SP)      :: ds1, ds2   !! grid width
        real(SP)      :: dt         !! time step
        integer       :: na1, na2   !! absoption area size
        integer       :: nmed       !! number of medium layers
        integer       :: nsnp
        real(SP)      :: clon, clat, phi
        real(SP)      :: L, w, C0  !! optional: not used

    end type fdsnap__hdr

contains

    subroutine fdsnap__open(fname, io, is_exist)

        !! Open fdsnap file

        character(*), intent(in)  :: fname
        integer, intent(out) :: io
        logical, intent(out) :: is_exist
        integer :: ierr
        
        !! first, try open as netcdf file
        ierr = nf90_open(fname, NF90_NOWRITE, io)
        is_exist = ierr == NF90_NOERR

    end subroutine fdsnap__open

    subroutine fdsnap__checkhdr(io, hdr)

        !! Write header for terminal

        integer, intent(in) :: io
        type(fdsnap__hdr), intent(in) :: hdr

        integer :: yr, mo, dy, hr, mi, sc
        character(100) :: fmt

        call daytim__localtime(hdr%exedate, yr, mo, dy, hr, mi, sc)

        fmt = '(A,I4,A,I2.2,A,I2.2,A,I2.2,A,I2.2,A,I2.2)'

        write (io, '(A,A)'    ) "[code type]     : ", hdr%codetype
        write (io, '(A,I10)'  ) "[header version]: ", hdr%hdrver
        write (io, '(A,A)'    ) "[title]         : ", trim(hdr%title)
        write (io, '(A,I10)'  ) "[date generated]: ", hdr%exedate
        write (io, fmt        ) "                  ", yr, "-", mo, "-", dy, "T", hr, "-", mi, "-", sc
        write (io, '(A,A)'    ) "[coordinate]    : ", hdr%coordinate
        write (io, '(A,A)'    ) "[data type]     : ", hdr%datatype
        write (io, '(A,I10)'  ) "[ns1]           : ", hdr%ns1
        write (io, '(A,I10)'  ) "[ns2]           : ", hdr%ns2
        write (io, '(A,F15.5)') "[beg1]          : ", hdr%beg1
        write (io, '(A,F15.5)') "[beg2]          : ", hdr%beg2
        write (io, '(A,F15.5)') "[ds1]           : ", hdr%ds1
        write (io, '(A,F15.5)') "[ds2]           : ", hdr%ds2
        write (io, '(A,F15.5)') "[dt]            : ", hdr%dt
        write (io, '(A,I10)'  ) "[na1]           : ", hdr%na1
        write (io, '(A,I10)'  ) "[na2]           : ", hdr%na2
        write (io, '(A,I10)'  ) "[nmed]          : ", hdr%nmed
        write (io, '(A,I10)'  ) "[nsnp]          : ", hdr%nsnp
        write (io, '(A,F15.5)') "[clon]          : ", hdr%clon
        write (io, '(A,F15.5)') "[clat]          : ", hdr%clat
        write (io, '(A,F15.5)') "[L]             : ", hdr%L
        write (io, '(A,F15.5)') "[w]             : ", hdr%w
        write (io, '(A,F15.5)') "[C0]            : ", hdr%C0

    end subroutine fdsnap__checkhdr

    subroutine fdsnap__readhdr(io, hdr)

        !! Read header part 
        integer, intent(in)  :: io
        type(fdsnap__hdr), intent(out) :: hdr

        
        call nc_chk(nf90_get_att(io, NF90_GLOBAL, 'codetype', hdr%codetype))
        call nc_chk(nf90_get_att(io, NF90_GLOBAL, 'hdrver', hdr%hdrver))
        call nc_chk(nf90_get_att(io, NF90_GLOBAL, 'title', hdr%title))
        call nc_chk(nf90_get_att(io, NF90_GLOBAL, 'exedate', hdr%exedate))
        call nc_chk(nf90_get_att(io, NF90_GLOBAL, 'coordinate', hdr%coordinate))
        call nc_chk(nf90_get_att(io, NF90_GLOBAL, 'datatype', hdr%datatype))
        call nc_chk(nf90_get_att(io, NF90_GLOBAL, 'ns1', hdr%ns1))
        call nc_chk(nf90_get_att(io, NF90_GLOBAL, 'ns2', hdr%ns2))
        call nc_chk(nf90_get_att(io, NF90_GLOBAL, 'beg1', hdr%beg1))
        call nc_chk(nf90_get_att(io, NF90_GLOBAL, 'beg2', hdr%beg2))
        call nc_chk(nf90_get_att(io, NF90_GLOBAL, 'ds1', hdr%ds1))
        call nc_chk(nf90_get_att(io, NF90_GLOBAL, 'ds2', hdr%ds2))
        call nc_chk(nf90_get_att(io, NF90_GLOBAL, 'dt', hdr%dt))
        call nc_chk(nf90_get_att(io, NF90_GLOBAL, 'na1', hdr%na1))
        call nc_chk(nf90_get_att(io, NF90_GLOBAL, 'na2', hdr%na2))
        call nc_chk(nf90_get_att(io, NF90_GLOBAL, 'nmed', hdr%nmed))
        call nc_chk(nf90_get_att(io, NF90_GLOBAL, 'nsnp', hdr%nsnp))
        call nc_chk(nf90_get_att(io, NF90_GLOBAL, 'clon', hdr%clon))
        call nc_chk(nf90_get_att(io, NF90_GLOBAL, 'clat', hdr%clat))
        call nc_chk(nf90_get_att(io, NF90_GLOBAL, 'phi', hdr%phi))

    end subroutine fdsnap__readhdr

    subroutine nc_chk(ierr)

        !! An internal subroutine to check error in netcdf function calls

        integer, intent(in) :: ierr
        
        if (ierr /= NF90_NOERR) write (error_unit, *) NF90_STRERROR(ierr)

    end subroutine nc_chk
    

end module m_fdsnap
