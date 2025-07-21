#include "../shared/m_debug.h"
program read_snp

    !! Read snap files from output of swpc, and export to figure
    !!
    !! Copyright 2013-2025 Takuto Maeda. All rights reserved. This project is released under the MIT license.

    use iso_fortran_env, only: error_unit
    use m_std
    use m_system
    use m_getopt
    use m_fdsnap
    use m_pnm
    use m_color
    use m_filt2d
    use m_geomap
    use m_debug
    use m_version
    use netcdf

    implicit none

    character(256) :: fn_snp

    type(fdsnap__hdr) :: hdr
    integer       :: io_snp
    logical       :: is_exist
    integer :: nx, ny
    integer :: iskip
    character(6) :: snp_type
    character(8) :: tname
    integer :: nt
    integer :: ierr
    integer, allocatable :: vid(:)
    character(100) :: vname
    integer :: i
    logical :: is_opt1, is_opt2

    call getopt('v', is_opt1)
    call getopt('-version', is_opt2)
    if (is_opt1 .or. is_opt2) call version__display('read_snp')

  !!
  !! Check Input File
  !!
    if (command_argument_count() == 0) then
        call usage_exit()
    end if

    call getopt('i', is_exist, fn_snp)
    if (.not. is_exist) then
        write (error_unit, '(A)') "ERROR [read_snp]: no input file given"
        write (error_unit, *)
        call usage_exit
    end if

    !! existence check
    call fdsnap__open(fn_snp, io_snp, is_exist, snp_type)
    if (.not. is_exist) then
        write (error_unit, '(A)') "ERROR [read_snp]: file "//trim(fn_snp)//" does not exist"
        write (error_unit, *)
        call usage_exit
    end if

    !! Read Header Part
    call fdsnap__readhdr(fn_snp, io_snp, snp_type, hdr)
    nx = hdr%ns1
    ny = hdr%ns2
    if (snp_type == 'netcdf') then
        allocate (vid(hdr%nsnp))
        do i = 1, hdr%nsnp
            call nc_chk(nf90_inquire_variable(io_snp, 3 + hdr%nmed + i, vname))
            call nc_chk(nf90_inq_varid(io_snp, vname, vid(i)))
            call nc_chk(nf90_inquire_dimension(io_snp, 3, tname, nt))
        end do
    end if

    !! header output mode
    call getopt('h', is_exist)
    if (is_exist) call fdsnap__checkhdr(error_unit, hdr)

    !! snap output to ppm file opotion
    call getopt('skip', is_exist, iskip, 0)

    call getopt('ppm', is_exist)
    if (is_exist) call img_output('ppm')

    call getopt('bmp', is_exist)
    if (is_exist) call img_output('bmp')

    call getopt('asc', is_exist)
    if (is_exist) call dat_output('asc')

    call getopt('bin', is_exist)
    if (is_exist) call dat_output('bin')

    if (snp_type == 'native') then
        close (io_snp)
    else
        ierr = nf90_close(io_snp)
    end if

contains

    subroutine usage_exit

        write (error_unit, *)
        write (error_unit, '(A)') ' read_snp.x -i snapshot [-h] '
        write (error_unit, '(A)') '      [-ppm|-bmp] [-pall] [-mul var | -mul1 var -mul2 var ...] '
        write (error_unit, '(A)') '      [-abs] [-bin|-asc] [-skip n] [-lpf ng] [-notim]'
        write (error_unit, *)
        write (error_unit, '(A)') '  -h: display header information to terminal output'
        write (error_unit, '(A)') '  -bmp: output bmp-formatted snapshot figures'
        write (error_unit, '(A)') '  -ppm: output ppm-formatted snapshot figures'
        write (error_unit, '(A)') '  -pall: plot including absorbing boundary area (clipped in default)'
        write (error_unit, '(A)') '  -mul var: scale amplitude by multiplying var for visualization'
        write (error_unit, '(A)') '  -mul1 var, -mul2 var ... : gives scaling factor by each component; default=1000'
        write (error_unit, '(A)') '  -abs: plot absolute value (only for velocity snapshot)'
        write (error_unit, '(A)') '  -bin: export single-precision xyz binary data'
        write (error_unit, '(A)') '  -asc: export xyz ascii data'
        write (error_unit, '(A)') '  -skip n: skip first n snapshots for export'
        write (error_unit, '(A)') '  -notim: do not plot elapsed time on the snapshort figures'
        write (error_unit, '(A)') '  -lpf ng: apply spatial low-pass filter with corner grid-width of ng before figure output' 
        write (error_unit, *)

        stop

    end subroutine usage_exit

    subroutine dat_output(typ)

        !! data file output

        implicit none

        character(3) :: typ ! asc or bin
        character(80), parameter :: D_ODIR = './dat'
        real(SP), allocatable :: den(:, :), rig(:, :), lam(:, :)
        real(SP), allocatable :: vp(:, :), vs(:, :)
        real(SP), allocatable :: amp(:, :, :)
        real(SP) :: xx(nx), yy(ny)
        integer :: it
        integer :: i, j
        integer :: ierr
        real    :: t
        character(80) :: odir = D_ODIR
        character(6)  :: cit    ! time grid number
        character(6)  :: ct     ! elapsed time
        character(3) :: ext
        character :: ci
        integer :: io
        logical :: is_eof
        character(80) :: fmt
        character(256) :: fn_dat
        integer :: vid_rho, vid_lambda, vid_mu
        integer :: start(3), count(3)
        character(3) :: codetype
        !--

        allocate (den(nx, ny), rig(nx, ny), lam(nx, ny))
        allocate (vp(nx, ny), vs(nx, ny))

        !! output directory
        if (trim(hdr%codetype) == 'SWPC_SH') then
            codetype = 'sh'
        else if (trim(hdr%codetype) == 'SWPC_PSV') then
            codetype = 'psv'
        else
            codetype = '3d'
        end if
        odir = trim(hdr%title)
        odir = trim(odir) // '/' // trim(codetype)
        odir = trim(odir) // '/' // trim(hdr%coordinate)
        odir = trim(odir) // '/' // trim(hdr%datatype)

        call execute_command_line('/bin/mkdir -p '//trim(odir)//' > /dev/null 2>&1 ')

        if (typ == 'asc') then
            ext = typ
        else if (typ == 'bin') then
            ext = typ
        else
            write (error_unit, *) 'unknown type'
            stop
        end if

        !! Medium structure
        if (snp_type == 'native') then
            read (io_snp) den, lam, rig
        else
            call nc_chk(nf90_inq_varid(io_snp, 'rho', vid_rho))
            call nc_chk(nf90_inq_varid(io_snp, 'lambda', vid_lambda))
            call nc_chk(nf90_inq_varid(io_snp, 'mu', vid_mu))
            call nc_chk(nf90_get_var(io_snp, vid_rho, den))
            call nc_chk(nf90_get_var(io_snp, vid_lambda, lam))
            call nc_chk(nf90_get_var(io_snp, vid_mu, rig))
        end if

        do j = 1, ny
            do i = 1, nx
                vs(i, j) = sqrt(rig(i, j) / den(i, j))
                vp(i, j) = sqrt((lam(i, j) + 2 * rig(i, j)) / den(i, j))
            end do
        end do

        do i = 1, nx
            xx(i) = hdr%beg1 + (i - 1) * hdr%ds1
        end do
        do j = 1, ny
            yy(j) = hdr%beg2 + (j - 1) * hdr%ds2
        end do

        fn_dat = trim(odir)//'/'//trim(hdr%title)//'.'//trim(hdr%coordinate)//'.'//trim(hdr%datatype)//'.vps.'//ext
        write (error_unit, *) trim(fn_dat)

        if (typ == 'asc') then
            open (newunit=io, file=fn_dat, action='write')

            fmt = "(2F15.5,3ES15.5)"

            do j = 1, ny
                do i = 1, nx
                    write (io, fmt) xx(i), yy(j), vp(i, j), vs(i, j), den(i, j)
                end do
                write (io, *)
            end do

            close (io)
        else
            open (newunit=io, file=fn_dat, action='write', access='stream', form='unformatted')
            do j = 1, ny
                do i = 1, nx
                    write (io) xx(i), yy(j), vp(i, j), vs(i, j), den(i, j)
                end do
            end do
            close (io)
        end if

        allocate (amp(hdr%nsnp, nx, ny))

        !! Shapshot Read/Data
        is_eof = .false.
        it = -iskip
        do
            t = (it + iskip) * hdr%dt

            if (snp_type == 'native') then
                do i = 1, hdr%nsnp
                    read (io_snp, iostat=ierr) amp(i, :, :)
                    if (ierr /= 0) then
                        write (error_unit, *) "EOF detected"
                        is_eof = .true.
                        exit
                    end if
                end do
                if (is_eof) exit
            else
                count = (/nx, ny, 1/)
                start = (/1, 1, it + iskip + 1/)
                do i = 1, hdr%nsnp
                    call nc_chk(nf90_get_var(io_snp, vid(i), amp(i, :, :), start=start, count=count))
                end do

                if (it == nt) then
                    is_eof = .true.
                    exit
                end if
            end if

            if (it >= 0) then
                write (ct, '(F6.1)') t
                write (cit, '(I6.6)') it

                fn_dat = trim(odir)//'/'&
                       //trim(hdr%title)//'.'&
                       //trim(codetype)//'.'&
                       //trim(hdr%coordinate)//'.'&
                       //trim(hdr%datatype)//'.'&
                       //cit//'.'//ext

                write (error_unit, *) trim(fn_dat)

                if (typ == 'asc') then
                    open (newunit=io, file=fn_dat, action='write')

                    write (ci, '(I1)') hdr%nsnp
                    fmt = "(2F15.5,"//ci//"ES15.5)"

                    do j = 1, ny
                        do i = 1, nx
                            write (io, fmt) xx(i), yy(j), amp(:, i, j)
                        end do
                        write (io, *)
                    end do

                    close (io)
                else
                    open (newunit=io, file=fn_dat, action='write', access='stream', form='unformatted')
                    do j = 1, ny
                        do i = 1, nx
                            write (io) xx(i), yy(j), amp(:, i, j)
                        end do
                    end do
                    close (io)
                end if
            end if

            it = it + 1
        end do

        deallocate (amp)

    end subroutine dat_output

    subroutine img_output(typ)

        use m_pnm
        use m_bmp
        use m_stamp
        implicit none

        character(3), intent(in) :: typ ! bmp or ppm

        real, parameter :: D_MUL = 1000.

        real(SP), allocatable :: den(:, :), rig(:, :), lam(:, :)
        real(SP), allocatable :: vp(:, :), vs(:, :), topo(:, :)
        integer, allocatable :: cmed(:, :, :)
        real(SP), allocatable :: medium_bound(:, :)

        real(SP), allocatable :: amp(:, :, :)
        integer, allocatable :: img(:, :, :)
        integer :: nxs, nys
        integer :: is, js
        integer :: it
        integer :: i, j
        integer :: ii, jj
        integer :: ierr
        logical :: is_exist
        logical :: is_eof
        logical :: is_abs
        real    :: amin, amax
        real    :: t
        real    :: wk
        real    :: mul(4)
        character(80) :: odir
        character(6)  :: cit    ! time grid number
        character(6)  :: ct     ! elapsed time
        logical :: is_transpose
        real(SP) :: mula
        type(color__palette) :: cp, cp_hot
        integer :: ctmp(3)
        real(SP) :: mingrid, kmax
        logical :: is_lpf
        real(SP) :: div, rot
        real(SP) :: sobel_edge, sobel_edge_x, sobel_edge_y
        real :: ud, horiz
        integer :: vid_rho, vid_mu, vid_lambda, vid_topo
        integer :: start(3), count(3)
        logical :: no_timemark
        character(3) :: codetype
        !--

        !! Memory allocation
        allocate (den(nx, ny), rig(nx, ny), lam(nx, ny))
        allocate (vp(nx, ny), vs(nx, ny), topo(nx, ny))
        allocate (cmed(3, nx, ny))
        allocate (medium_bound(nx, ny))

        !! output directory
        if (trim(hdr%codetype) == 'SWPC_SH') then
            codetype = 'sh'
        else if (trim(hdr%codetype) == 'SWPC_PSV') then
            codetype = 'psv'
        else
            codetype = '3d'
        end if
        odir = trim(hdr%title)
        odir = trim(odir) // '/' // trim(codetype)
        odir = trim(odir) // '/' // trim(hdr%coordinate)
        odir = trim(odir) // '/' // trim(hdr%datatype)
        
        call execute_command_line('/bin/mkdir -p '//trim(odir)//' > /dev/null 2>&1 ')

        !! lowpass
        call getopt('lpf', is_lpf, mingrid, 2.0)
        kmax = 2 * PI / (mingrid * sqrt(hdr%ds1 * hdr%ds2) / sqrt(2.0))

        !! amplitude scale

        !! independent amplitude weight
        call getopt('mul1', is_exist, mul(1), D_MUL)
        call getopt('mul2', is_exist, mul(2), D_MUL)
        call getopt('mul3', is_exist, mul(3), D_MUL)
        call getopt('mul4', is_exist, mul(4), D_MUL)

        !! ... or equall weight
        call getopt('mul', is_exist, mula)
        if (is_exist) mul(:) = mula

        call getopt('abs', is_abs) ! absolute value

        !! Time mark (default is ON)
        call getopt('notim', no_timemark)

        !! graph size ( if includes absorbing boundary layer )
        call getopt('pall', is_exist)
        if (is_exist) then
            nxs = nx
            nys = ny
            is = 1
            js = 1
        else
            if (hdr%coordinate == 'fs' .or. hdr%coordinate == 'ob' .or. hdr%coordinate == 'xy') then
                nxs = nx - 2 * hdr%na1
                nys = ny - 2 * hdr%na2
                is = hdr%na1 + 1
                js = hdr%na2 + 1
            else
                nxs = nx - 2 * hdr%na1
                nys = ny - 2 * hdr%na2
                is = hdr%na1 + 1
                js = hdr%na2 + 1
            end if
        end if
        if (hdr%coordinate == 'fs' .or. hdr%coordinate == 'ob' .or. hdr%coordinate == 'xy') then
            is_transpose = .true.
        else
            is_transpose = .false.
        end if

        !! Medium structure
        if (snp_type == 'native') then
            read (io_snp) den, lam, rig
        else
            call nc_chk(nf90_inq_varid(io_snp, 'rho', vid_rho))
            call nc_chk(nf90_inq_varid(io_snp, 'lambda', vid_lambda))
            call nc_chk(nf90_inq_varid(io_snp, 'mu', vid_mu))
            call nc_chk(nf90_get_var(io_snp, vid_rho, den))
            call nc_chk(nf90_get_var(io_snp, vid_lambda, lam))
            call nc_chk(nf90_get_var(io_snp, vid_mu, rig))
        end if

        do j = 1, ny
            do i = 1, nx
                vs(i, j) = sqrt(rig(i, j) / den(i, j))
                vp(i, j) = sqrt((lam(i, j) + 2 * rig(i, j)) / den(i, j))
            end do
        end do
        !! horizontal case: read topography
        if (hdr%coordinate == 'fs' .or. hdr%coordinate == 'ob' .or. hdr%coordinate == 'xy') then
            if (snp_type == 'native') then
                read (io_snp) topo
            else
                call nc_chk(nf90_inq_varid(io_snp, 'topo', vid_topo))
                call nc_chk(nf90_get_var(io_snp, vid_topo, topo))
            end if
        end if

        !! Background Coloring
        medium_bound(:, :) = 1
        amin = minval(vp(:, :))
        amax = maxval(vp(:, :)) + 0.01

        if (is_abs) call color__set('hot', cp_hot)

        if (hdr%coordinate == 'fs' .or. hdr%coordinate == 'ob') then

            !! surface map: use topography data for background color

            call color__set('mytopo', cp)

            do j = 1, ny
                do i = 1, nx
                    call color__interpolate(cp, dble(topo(i, j)), cmed(:, i, j))
                end do
            end do

        else

            !! other cases: coloring by density with discontinuous structure detection
            do j = 1, ny
                do i = 1, nx
                    wk = 1.2 - 0.1 * exp(1.2 * ((vp(i, j) - amin) / (amax - amin)))
                    if (vs(i, j) > 0.1) then ! solid
                        cmed(1, i, j) = min(max(int(255 * wk), 0), 255)
                        cmed(2, i, j) = min(max(int(255 * wk), 0), 255)
                        cmed(3, i, j) = min(max(int(220 * wk), 0), 255)
                    else if (vp(i, j) > 0.1) then ! ocean
                        cmed(1, i, j) = int(220 * wk)
                        cmed(2, i, j) = int(235 * wk)
                        cmed(3, i, j) = int(255 * wk)
                    else
                        cmed(:, i, j) = 255 ! air/vaccum
                    end if
                end do
            end do

            !! edge detection (Sobel's edge detection operator)
            do i = 2, nx - 1
                do j = 2, ny - 1

                    sobel_edge_x = den(i + 1, j + 1) + 2 * den(i + 1, j) + den(i + 1, j - 1) &
                                   - den(i - 1, j + 1) - 2 * den(i - 1, j) - den(i - 1, j - 1)
                    sobel_edge_y = den(i + 1, j + 1) + 2 * den(i, j + 1) + den(i - 1, j + 1) &
                                   - den(i + 1, j - 1) - 2 * den(i, j - 1) - den(i - 1, j - 1)
                    sobel_edge = sqrt(sobel_edge_x**2 + sobel_edge_y**2)
                    if (sobel_edge > 0.) then
                        medium_bound(i, j) = max((1 - sobel_edge / 10)**2, 0.)

                    end if
                end do
            end do
            medium_bound(1, :) = medium_bound(2, :)
            medium_bound(nx, :) = medium_bound(nx - 1, :)
            medium_bound(:, 1) = medium_bound(:, 2)
            medium_bound(:, ny) = medium_bound(:, ny - 1)

        end if

        allocate (amp(hdr%nsnp, nx, ny))
        if (is_transpose) then
            allocate (img(3, nys, nxs))
        else
            allocate (img(3, nxs, nys))
        end if

        it = -iskip

        img(1, :, :) = 100
        img(2, :, :) = 100
        img(3, :, :) = 100

        is_eof = .false.
        do
            t = (it + iskip) * hdr%dt

            if (snp_type == 'native') then
                do i = 1, hdr%nsnp
                    read (io_snp, iostat=ierr) amp(i, :, :)
                    if (ierr /= 0) then
                        write (error_unit, *) "EOF detected"
                        is_eof = .true.
                        exit
                    end if

                end do
                if (is_eof) exit
            else
                if (it == nt) then
                    is_eof = .true.
                    exit
                end if
                count = (/nx, ny, 1/)
                start = (/1, 1, it + iskip + 1/)
                do i = 1, hdr%nsnp
                    call nc_chk(nf90_get_var(io_snp, vid(i), amp(i, :, :), start=start, count=count))
                end do
            end if

            if (it < 0) then
                it = it + 1
                cycle
            end if

            if (is_lpf) then
                do i = 1, hdr%nsnp
                    call filt2d__lowpass(nx, ny, hdr%ds1, hdr%ds2, amp(i, :, :), kmax, 2)
                end do
            end if

            write (ct, '(F6.1)') t
            write (cit, '(I6.6)') it

            fn_snp =  trim(odir)//'/'&
                    //trim(hdr%title)//'.' &
                    //trim(codetype) // '.' &
                    //trim(hdr%coordinate)//'.'&
                    //trim(hdr%datatype)//'.'&
                    //cit//'.'//typ

            write (error_unit, *) trim(fn_snp)

            do i = 1, hdr%nsnp
                amp(i, :, :) = mul(i) * abs(amp(i, :, :))
            end do

            !! coloring
            do j = js + 1, ny - js
                do i = is + 1, nx - is

                    if (is_transpose) then
                        ii = j - js + 1
                        jj = nxs - (i - is + 1) + 1

                    else
                        ii = i - is + 1
                        jj = j - js + 1
                    end if

                    if (hdr%datatype == "ps") then

                        ! wave
                        div = abs(amp(1, i, j))
                        rot = sqrt(sum(amp(2:hdr%nsnp, i, j)**2)) ! include psv and 3D

                        img(1, ii, jj) = cmed(1, i, j) - int(255 * rot) / 4
                        img(2, ii, jj) = cmed(2, i, j) - int(255 * div) / 2
                        img(3, ii, jj) = cmed(3, i, j) - int(255 * (div + rot)) / 3

                    else if (hdr%datatype == "v3" .or. hdr%datatype == "u3") then

                        if (is_abs) then

                            !! absolute value plot by hot color palette
                            call color__interpolate(cp_hot, dble(sqrt(sum(amp(:, i, j)**2))), ctmp(:))

                            !! reducing color
                            img(1, ii, jj) = cmed(1, i, j) - (255 - ctmp(1))
                            img(2, ii, jj) = cmed(2, i, j) - (255 - ctmp(2))
                            img(3, ii, jj) = cmed(3, i, j) - (255 - ctmp(3))

                        else

                            ud = abs(amp(3, i, j))
                            horiz = sqrt(amp(1, i, j)**2 + amp(2, i, j)**2)

                            img(1, ii, jj) = cmed(1, i, j) - int(255 * horiz) / 4
                            img(2, ii, jj) = cmed(2, i, j) - int(255 * ud) / 2
                            img(3, ii, jj) = cmed(3, i, j) - int(255 * (ud + horiz)) / 3

                        end if

                    else if (hdr%datatype == "v2" .or. hdr%datatype == "u2") then

                        if (is_abs) then

                            !! absolute value plot by hot color palette
                            call color__interpolate(cp_hot, dble(sqrt(sum(amp(:, i, j)**2))), ctmp(:))

                            !! reducing color
                            img(1, ii, jj) = cmed(1, i, j) - (255 - ctmp(1))
                            img(2, ii, jj) = cmed(2, i, j) - (255 - ctmp(2))
                            img(3, ii, jj) = cmed(3, i, j) - (255 - ctmp(3))

                        else

                            ud = abs(amp(2, i, j))
                            horiz = abs(amp(1, i, j))

                            img(1, ii, jj) = cmed(1, i, j) - int(255 * horiz) / 4
                            img(2, ii, jj) = cmed(2, i, j) - int(255 * ud) / 2
                            img(3, ii, jj) = cmed(3, i, j) - int(255 * (ud + horiz)) / 3

                        end if

                    else if (hdr%datatype == "vy" .or. hdr%datatype == "uy") then

                        ud = 0
                        horiz = abs(amp(1, i, j))

                        img(1, ii, jj) = cmed(1, i, j) - int(255 * horiz) / 4
                        img(2, ii, jj) = cmed(2, i, j) - int(255 * ud) / 2
                        img(3, ii, jj) = cmed(3, i, j) - int(255 * (ud + horiz)) / 3

                    end if

                    ! normalize
                    img(1, ii, jj) = max(20, min(255, img(1, ii, jj)))
                    img(2, ii, jj) = max(20, min(255, img(2, ii, jj)))
                    img(3, ii, jj) = max(20, min(255, img(3, ii, jj)))

                    ! medium boudary
                    img(1, ii, jj) = int(img(1, ii, jj) * medium_bound(i, j))
                    img(2, ii, jj) = int(img(2, ii, jj) * medium_bound(i, j))
                    img(3, ii, jj) = int(img(3, ii, jj) * medium_bound(i, j))

                end do
            end do

            if (is_transpose) then
                !! surrounding line
                img(:, 1:nys, 1:2) = 0
                img(:, 1:nys, nxs - 1:nxs) = 0
                img(:, 1:2, 1:nxs) = 0
                img(:, nys - 1:nys, 1:nxs) = 0

                !! timemark
                if (.not. no_timemark) then
                    call stamp__char("t = "//trim(ct)//' s', 20, nxs - 40, nys, nxs, img, .false.)
                end if

                !! export
                if (typ == 'ppm') then
                    call ppm__write(trim(fn_snp), nys, nxs, img)
                else if (typ == 'bmp') then
                    call bmp__write(trim(fn_snp), nys, nxs, img)
                end if

            else

                !! surrounding line
                img(:, 1:nxs, 1:2) = 0
                img(:, 1:nxs, nys - 1:nys) = 0
                img(:, 1:2, 1:nys) = 0
                img(:, nxs - 1:nxs, 1:nys) = 0

                !! timemark
                if (.not. no_timemark) then
                    call stamp__char("t = "//trim(ct)//'s', 20, nys - 40, nxs, nys, img, .false.)
                end if

                !! export
                if (typ == 'ppm') then
                    call ppm__write(trim(fn_snp), nxs, nys, img)
                else if (typ == 'bmp') then
                    call bmp__write(trim(fn_snp), nxs, nys, img)
                end if

            end if

            it = it + 1
        end do

        deallocate (amp, img)

    end subroutine img_output

    subroutine nc_chk(ierr)

        integer, intent(in) :: ierr

        if (ierr /= NF90_NOERR) write (error_unit, *) NF90_STRERROR(ierr)

    end subroutine nc_chk

end program read_snp
!------------------------------------------------------------------------------!
