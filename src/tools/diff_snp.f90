
#include "../shared/m_debug.h"
program diff_snp

    !! Generate differential snapfile from two inputs
    !!
    !! Copyright 2013-2026 Takuto Maeda. All rights reserved. This project is released under the MIT license.

    use iso_fortran_env, only: error_unit
    use m_std
    use m_system
    use m_fdsnap
    use m_daytim
    use m_debug
    use m_version
    use netcdf

    implicit none

    character(256)        :: fn_in1, fn_in2, fn_out
    character(256)        :: vname
    type(fdsnap__hdr)     :: hdr_in1, hdr_in2, hdr_out
    integer               :: io_in1, io_in2, io_out
    integer               :: i, it, nx, ny, nt
    integer               :: ct(3), st(3)
    logical               :: is_exist
    integer, allocatable :: vid(:)
    real(SP), allocatable :: amp1(:, :), amp2(:, :)

    !! Open input files
    call system__getarg(1, fn_in1)
    if (trim(fn_in1) == '-v' .or. trim(fn_in1) == '--version') call version__display('diff_snp')
    call system__getarg(2, fn_in2)
    call system__getarg(3, fn_out)

    call fdsnap__open(fn_in1, io_in1, is_exist)
    if (.not. is_exist) stop
    call fdsnap__open(fn_in2, io_in2, is_exist)
    if (.not. is_exist) stop

    !! Read Header Part
    call fdsnap__readhdr(io_in1, hdr_in1)
    call fdsnap__readhdr(io_in2, hdr_in2)

    !! Check size consistency
    call fdsnap__checkhdr(error_unit, hdr_in1)
    call fdsnap__checkhdr(error_unit, hdr_in2)

    call assert(hdr_in1%ns1 == hdr_in2%ns1 .and. hdr_in1%ns2 == hdr_in2%ns2)
    nx = hdr_in1%ns1
    ny = hdr_in1%ns2

    allocate (amp1(nx, ny), amp2(nx, ny))

    !! diff file

    !! file generation by a simple copy, then modify in what follows
    call execute_command_line('/bin/cp '//trim(fn_in1)//' '//trim(fn_out))
    call nc_chk(nf90_open(fn_out, NF90_WRITE, io_out))

    !! date & time
    call nc_chk(nf90_put_att(io_out, NF90_GLOBAL, 'exedate', hdr_out%exedate))
    call nc_chk(nf90_inquire_dimension(io_in1, 3, vname, nt))

    allocate (vid(hdr_in1%nsnp))
    do i = 1, hdr_in1%nsnp
        call nc_chk(nf90_inquire_variable(io_in1, 3 + hdr_in1%nmed + i, vname))
        call nc_chk(nf90_inq_varid(io_in1, vname, vid(i)))
    end do

    do it = 1, nt
        write (error_unit, *) it
        ct = (/nx, ny, 1/)
        st = (/1, 1, it/)
        do i = 1, hdr_in1%nsnp
            call nc_chk(nf90_get_var(io_in1, vid(i), amp1, start=st, count=ct))
            call nc_chk(nf90_get_var(io_in2, vid(i), amp2, start=st, count=ct))
            call nc_chk(nf90_put_var(io_out, vid(i), amp1 - amp2, start=st, count=ct))
        end do
    end do

    call nc_chk(nf90_close(io_in1))
    call nc_chk(nf90_close(io_in2))
    call nc_chk(nf90_close(io_out))

contains

    subroutine nc_chk(ierr)

        integer, intent(in) :: ierr

        if (ierr /= NF90_NOERR) write (error_unit, *) NF90_STRERROR(ierr)

    end subroutine nc_chk

end program diff_snp
