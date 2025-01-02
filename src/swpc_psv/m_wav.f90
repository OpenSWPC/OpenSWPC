#include "../shared/m_debug.h"
module m_wav

    !! waveform output
    !!
    !! Copyright 2024 Takuto Maeda. All rights reseaved. This project is released under the MIT license.

    use m_std
    use m_debug
    use m_global
    use m_pwatch
    use m_sac
    use m_readini
    use m_geomap
    use m_tar
    implicit none
    private
    save

    public :: wav__setup
    public :: wav__store
    public :: wav__write

    integer :: ntdec_w
    integer :: ntdec_w_prg
    character(8) :: wav_format

    logical :: sw_wav_v = .false.
    logical :: sw_wav_u = .false.
    logical :: sw_wav_stress = .false.
    logical :: sw_wav_strain = .false.

    integer :: ntw ! number of wave samples
    real(SP), allocatable :: wav_vel(:, :, :), wav_disp(:, :, :)
    real(SP), allocatable :: wav_stress(:, :, :), wav_strain(:, :, :) ! (6,nt,nst)
    type(sac__hdr), allocatable :: sh_vel(:, :), sh_disp(:, :)
    type(sac__hdr), allocatable :: sh_stress(:, :), sh_strain(:, :)
    real(SP), allocatable :: ux(:), uz(:)
    real(SP), allocatable :: exx(:), ezz(:), exz(:)

    integer :: nst = 0
    real(SP), allocatable :: xst(:), zst(:)
    integer, allocatable :: ist(:), kst(:)
    real(SP), allocatable :: stlo(:), stla(:)
    character(8), allocatable :: stnm(:)

    real(MP) :: r40x, r40z, r41x, r41z, r20x, r20z

contains

    subroutine wav__setup(io_prm)

        integer, intent(in) :: io_prm
        !! --
        character(245) :: fn_stloc
        character(2) :: st_format
        !! ----

        call pwatch__on('wav__setup')

        call readini(io_prm, 'ntdec_w', ntdec_w, 10)
        call readini(io_prm, 'sw_wav_v', sw_wav_v, .false.)
        call readini(io_prm, 'sw_wav_u', sw_wav_u, .false.)
        call readini(io_prm, 'sw_wav_stress', sw_wav_stress, .false.)
        call readini(io_prm, 'sw_wav_strain', sw_wav_strain, .false.)
        call readini(io_prm, 'wav_format', wav_format, 'sac')
        call readini(io_prm, 'st_format', st_format, 'xy')
        call readini(io_prm, 'fn_stloc', fn_stloc, '')
        call readini(io_prm, 'ntdec_w_prg', ntdec_w_prg, 0)

        if (.not. (sw_wav_v .or. sw_wav_u .or. sw_wav_stress .or. sw_wav_strain)) then
            nst = 0
            call pwatch__off('wav__setup')
            return
        end if

        ntw = floor(real(nt - 1) / real(ntdec_w) + 1.0)

        !! FDM coefficients
        r40x = 9.0_MP / 8.0_MP / dx
        r40z = 9.0_MP / 8.0_MP / dz
        r41x = 1.0_MP / 24.0_MP / dx
        r41z = 1.0_MP / 24.0_MP / dz
        r20x = 1./dx
        r20z = 1./dz

        call set_stinfo(fn_stloc, st_format)

        call execute_command_line('mkdir -p '//trim(odir)//'/wav > /dev/null 2>&1')

        if (sw_wav_v) then
            allocate (wav_vel(ntw, 2, nst), source=0.0)
            allocate (sh_vel(2, nst))
        end if

        if (sw_wav_u) then
            allocate (wav_disp(ntw, 2, nst), source=0.0)
            allocate (sh_disp(2, nst))
        end if

        if (sw_wav_stress) then
            allocate (wav_stress(ntw, 3, nst), source=0.0)
            allocate (sh_stress(3, nst))
        end if

        if (sw_wav_strain) then
            allocate (wav_strain(ntw, 3, nst), source=0.0)
            allocate (sh_strain(3, nst))
        end if

        call set_sac_header()

        !! FDM coefficients
        r40x = 9.0_MP / 8.0_MP / dx
        r40z = 9.0_MP / 8.0_MP / dz
        r41x = 1.0_MP / 24.0_MP / dx
        r41z = 1.0_MP / 24.0_MP / dz

        call pwatch__off('wav__setup')

    end subroutine wav__setup

    subroutine wav__store(it)

        integer, intent(in) :: it
        !! --
        integer :: n, itw
        real(MP) :: dxVx, dxVz, dzVx, dzVz
        integer :: i, k
        !! ----

        call pwatch__on("wav__store")
        if (nst == 0) then
            call pwatch__off("wav__store")
            return
        end if

        if (it == 1) then
            if (sw_wav_u) then
                allocate (ux(nst), source=0.0)
                allocate (uz(nst), source=0.0)
            end if

            if (sw_wav_strain) then
                allocate (exx(nst), source=0.0)
                allocate (ezz(nst), source=0.0)
                allocate (exz(nst), source=0.0)
            end if
        end if

        if (sw_wav_u) then

            !$omp parallel do private(n,i,k)
            do n = 1, nst
                i = ist(n); k = kst(n)
                ux(n) = ux(n) + real(Vx(k, i) + Vx(k, i - 1)) * 0.5 * dt
                uz(n) = uz(n) - real(Vz(k, i) + Vz(k - 1, i)) * 0.5 * dt ! positive upward
            end do
            !$omp end parallel do

        end if

        if (sw_wav_strain) then

            !$omp parallel do private(n, i, k, dxVx, dzVz, dxVz, dzVx)
            do n = 1, nst
                i = ist(n); k = kst(n)

                dxVx = (Vx(k, i) - Vx(k, i - 1)) * r40x - (Vx(k, i + 1) - Vx(k, i - 2)) * r41x
                dzVz = (Vz(k, i) - Vz(k - 1, i)) * r40z - (Vz(k + 1, i) - Vz(k - 2, i)) * r41z
                dxVz = ((Vz(k, i + 1) - Vz(k, i)) * r40x - (Vz(k, i + 2) - Vz(k, i - 1)) * r41x &
                        + (Vz(k - 1, i + 1) - Vz(k - 1, i)) * r40x - (Vz(k - 1, i + 2) - Vz(k - 1, i - 1)) * r41x &
                        + (Vz(k, i) - Vz(k, i - 1)) * r40x - (Vz(k, i + 1) - Vz(k, i - 2)) * r41x &
                        + (Vz(k - 1, i) - Vz(k - 1, i - 1)) * r40x - (Vz(k - 1, i + 1) - Vz(k - 1, i - 2)) * r41x) / 4.0
                dzVx = ((Vx(k + 1, i) - Vx(k, i)) * r40z - (Vx(k + 2, i) - Vx(k - 1, i)) * r41z &
                        + (Vx(k + 1, i - 1) - Vx(k, i - 1)) * r40z - (Vx(k + 2, i - 1) - Vx(k - 1, i - 1)) * r41z &
                        + (Vx(k, i) - Vx(k - 1, i)) * r40z - (Vx(k + 1, i) - Vx(k - 2, i)) * r41z &
                        + (Vx(k, i - 1) - Vx(k - 1, i - 1)) * r40z - (Vx(k + 1, i - 1) - Vx(k - 2, i - 1)) * r41z) / 4.0

                exx(n) = exx(n) + real(dxVx) * dt
                ezz(n) = ezz(n) + real(dzVz) * dt
                exz(n) = exz(n) + real(dxVz + dzVx) / 2.0 * dt

            end do
            !$omp end parallel do
        end if

        if (mod(it - 1, ntdec_w) == 0) then

            itw = (it - 1) / ntdec_w + 1
            if (sw_wav_v) then

                !$omp parallel do private(n,i,k)
                do n = 1, nst
                    i = ist(n); k = kst(n)
                    wav_vel(itw, 1, n) = real(Vx(k, i) + Vx(k, i - 1)) / 2.0 * M0 * UC * 1e9
                    wav_vel(itw, 2, n) = -real(Vz(k, i) + Vz(k - 1, i)) / 2.0 * M0 * UC * 1e9
                end do
                !$omp end parallel do

            end if

            if (sw_wav_u) then

                !$omp parallel do private(n)
                do n = 1, nst
                    wav_disp(itw, 1, n) = ux(n) * M0 * UC * 1e9
                    wav_disp(itw, 2, n) = uz(n) * M0 * UC * 1e9
                end do
                !$omp end parallel do

            end if

            if (sw_wav_stress) then

                !$omp parallel do private(n, i, k)
                do n = 1, nst
                    i = ist(n); k = kst(n)
                    wav_stress(itw, 1, n) = real(Sxx(k, i)) * M0 * UC * 1e6
                    wav_stress(itw, 2, n) = real(Szz(k, i)) * M0 * UC * 1e6
                    wav_stress(itw, 3, n) = real(Sxz(k, i) + Sxz(k, i - 1) &
                                                 + Sxz(k - 1, i) + Sxz(k - 1, i - 1)) / 4.0 * M0 * UC * 1e6
                end do
                !$omp end parallel do

            end if

            if (sw_wav_strain) then

                !$omp parallel do private(n)
                do n = 1, nst
                    wav_strain(itw, 1, n) = exx(n) * M0 * UC * 1e-3
                    wav_strain(itw, 2, n) = ezz(n) * M0 * UC * 1e-3
                    wav_strain(itw, 3, n) = exz(n) * M0 * UC * 1e-3
                end do
                !$omp end parallel do

            end if

        end if

        if (ntdec_w_prg > 0 .and. mod(it - 1, ntdec_w_prg) == 0) call wav__write()

        call pwatch__off("wav__store")

    end subroutine wav__store

    subroutine wav__write()

        integer :: i, j
        character(256) :: fn_tar
        character(6) :: cid
        integer :: io
        !! ----

        call pwatch__on("wav__write")

        if (nst == 0) then
            call pwatch__off("wav__write")
            return
        end if

        if (wav_format == 'sac') then

            do i = 1, nst

                if (sw_wav_v) then
                    do j = 1, 2
                        call export_wav__sac(sh_vel(j, i), wav_vel(:, j, i))
                    end do
                end if

                if (sw_wav_u) then
                    do j = 1, 2
                        call export_wav__sac(sh_disp(j, i), wav_disp(:, j, i))
                    end do
                end if

                if (sw_wav_stress) then
                    do j = 1, 3
                        call export_wav__sac(sh_stress(j, i), wav_stress(:, j, i))
                    end do
                end if

                if (sw_wav_strain) then
                    do j = 1, 3
                        call export_wav__sac(sh_strain(j, i), wav_strain(:, j, i))
                    end do
                end if

            end do

        else if (wav_format == 'csf') then

            if (sw_wav_v) call export_wav__csf(nst, 2, sh_vel, wav_vel)
            if (sw_wav_u) call export_wav__csf(nst, 2, sh_disp, wav_disp)
            if (sw_wav_stress) call export_wav__csf(nst, 3, sh_stress, wav_stress)
            if (sw_wav_strain) call export_wav__csf(nst, 3, sh_strain, wav_strain)

        else if (trim(wav_format) == 'tar_st' .or. trim(wav_format) == 'tar_node') then

            if (trim(wav_format) == 'tar_node') then
                write (cid, '(I6.6)') myid
                fn_tar = trim(odir)//'/wav/'//trim(title)//'.'//cid//'.sac.tar'
                open(newunit=io, file=fn_tar, action='write', access='stream', status='unknown')
            end if

            do i=1, nst

                if (trim(wav_format) == 'tar_st') then
                    fn_tar = trim(odir)//'/wav/'//trim(title)//'.'//trim(stnm(i))//'.sac.tar'
                    open(newunit=io, file=fn_tar, action='write', access='stream', status='unknown')
                end if

                if (sw_wav_v) then
                    do j = 1, 2
                        call export_wav__tar(io, sh_vel(j, i), wav_vel(:, j, i))
                    end do
                end if

                if (sw_wav_u) then
                    do j = 1, 2
                        call export_wav__tar(io, sh_disp(j, i), wav_disp(:, j, i))
                    end do
                end if

                if (sw_wav_stress) then
                    do j = 1, 3
                        call export_wav__tar(io, sh_stress(j, i), wav_stress(:, j, i))
                    end do
                end if

                if (sw_wav_strain) then
                    do j = 1, 3
                        call export_wav__tar(io, sh_strain(j, i), wav_strain(:, j, i))
                    end do
                end if                

                if(trim(wav_format) == 'tar_st') then
                    call tar__wend(io)
                    close(io)
                end if        

            end do

        end if

        call pwatch__off("wav__write")

    end subroutine wav__write


    subroutine export_wav__tar(io, sh, dat)

        integer, intent(in) :: io
        type(sac__hdr), intent(in) :: sh
        real(SP), intent(in) :: dat(:)

        character(256) :: fn

        fn = trim(title)//'.'//trim(sh%kstnm)//'.'//trim(sh%kcmpnm)//'.sac'
        call sac__wtar(io, trim(fn), sh, dat)

    end subroutine export_wav__tar


    subroutine set_stinfo(fn_stloc, st_format)

        character(*), intent(in) :: fn_stloc
        character(*), intent(in) :: st_format
        !! --
        integer :: io_stlst
        integer :: err
        integer :: nst_g
        character(256) :: abuf
        real(SP) :: xst_g, zst_g, stlo_g, stla_g
        integer  :: ist_g, kst_g
        character(3) :: zsw_g
        character(8) :: stnm_g
        real(SP) :: rdum
        !! ----

        open (newunit=io_stlst, file=trim(fn_stloc), action='read', status='old', iostat=err)

        if (err /= 0) then
            call info('no station location file found')
            nst = 0
            return
        end if

        do
            read (io_stlst, '(a256)', iostat=err) abuf
            if (err /= 0) exit

            abuf = trim(adjustl(abuf))
            if (abuf(1:1) == "#") cycle
            if (trim(adjustl(abuf)) == "") cycle

            select case (st_format)

            case ('xy')

                read (abuf, *) xst_g, rdum, zst_g, stnm_g, zsw_g
                call geomap__c2g(xst_g, 0.0, clon, clat, phi, stlo_g, stla_g)

            case ('ll')
                read (abuf, *) stlo_g, stla_g, zst_g, stnm_g, zsw_g
                call assert(-360. <= stlo_g .and. stlo_g <= 360.)
                call assert(-90. <= stla_g .and. stla_g <= 90.)
                call geomap__g2c(stlo_g, stla_g, clon, clat, phi, xst_g, rdum)

            case default

                call info('unknown st_format: '//st_format)
                call assert(.false.)

            end select

            ist_g = x2i(xst_g, xbeg, real(dx))
            kst_g = z2k(zst_g, zbeg, real(dz))

            if (i2x(1, xbeg, real(dx)) < xst_g .and. xst_g < i2x(nx, xbeg, real(dx)) .and. &
                kbeg < kst_g .and. kst_g < kend) then

                nst_g = nst_g + 1

                if (ibeg <= ist_g .and. ist_g <= iend) then

                    !! station depth setting
                    select case (zsw_g)
                    case ('dep'); kst_g = z2k(zst_g, zbeg, real(dz))
                    case ('fsb'); kst_g = kfs(ist_g) + 1  !! free surface: one-grid below for avoiding vacuum
                    case ('obb'); kst_g = kob(ist_g) + 1   !! ocean column: below seafloor
                    case ('oba'); kst_g = kob(ist_g) - 1   !! ocean column: above seafloor
                    case ('bd0'); kst_g = z2k(bddep(ist_g, 0), zbeg, real(dz))
                    case ('bd1'); kst_g = z2k(bddep(ist_g, 1), zbeg, real(dz))
                    case ('bd2'); kst_g = z2k(bddep(ist_g, 2), zbeg, real(dz))
                    case ('bd3'); kst_g = z2k(bddep(ist_g, 3), zbeg, real(dz))
                    case ('bd4'); kst_g = z2k(bddep(ist_g, 4), zbeg, real(dz))
                    case ('bd5'); kst_g = z2k(bddep(ist_g, 5), zbeg, real(dz))
                    case ('bd6'); kst_g = z2k(bddep(ist_g, 6), zbeg, real(dz))
                    case ('bd7'); kst_g = z2k(bddep(ist_g, 7), zbeg, real(dz))
                    case ('bd8'); kst_g = z2k(bddep(ist_g, 8), zbeg, real(dz))
                    case ('bd9'); kst_g = z2k(bddep(ist_g, 9), zbeg, real(dz))
                    case default
                        call info("unknown zsw type in station file. Assume 'dep'")
                        kst_g = z2k(zst_g, zbeg, real(dz))
                    end select

                    if (kst_g > kend) then
                        call info("station depth exceeds kend at station "//trim(stnm_g))
                        kst_g = kend - 1
                    end if
                    if (kst_g < kbeg) then
                        call info("station depth fall short of kbeg at station"//trim(stnm_g))
                        kst_g = kbeg + 1
                    end if

                    nst = nst + 1
                    if (nst == 1) then
                        allocate (xst(1), source=xst_g)
                        allocate (zst(1), source=zst_g)
                        allocate (stnm(1), source=stnm_g)
                        allocate (stlo(1), source=stlo_g)
                        allocate (stla(1), source=stla_g)
                        allocate (ist(1), source=ist_g)
                        allocate (kst(1), source=kst_g)
                    else
                        call std__extend_array(xst, xst_g)
                        call std__extend_array(zst, zst_g)
                        call std__extend_array(8, stnm, stnm_g)
                        call std__extend_array(stlo, stlo_g)
                        call std__extend_array(stla, stla_g)
                        call std__extend_array(ist, ist_g)
                        call std__extend_array(kst, kst_g)
                    end if
                end if
            else
                if (myid == 0) call info("station "//trim(stnm_g)//" is out of the region")
            end if
        end do

        close (io_stlst)

        if (nst_g == 0) then
            nst = 0
            if (myid == 0) call info("no station is detected. waveform files will not be created")
            return
        end if

    end subroutine set_stinfo

    subroutine set_sac_header()

        integer :: i, j
        real :: mag
        !! ----

        mag = moment_magnitude(M0)

        do i = 1, nst

            if (sw_wav_v) then
                do j = 1, 2
                    call initialize_sac_header(sh_vel(j, i), stnm(i), stlo(i), stla(i), xst(i), zst(i), mag)
                end do
                sh_vel(1, i)%kcmpnm = "Vx"; sh_vel(1, i)%cmpinc = 90.0; sh_vel(1, i)%cmpaz = 0.0 + phi
                sh_vel(2, i)%kcmpnm = "Vz"; sh_vel(2, i)%cmpinc = 90.0; sh_vel(2, i)%cmpaz = 0.0

                sh_vel(:, i)%idep = 7
            end if

            if (sw_wav_u) then
                do j = 1, 2
                    call initialize_sac_header(sh_disp(j, i), stnm(i), stlo(i), stla(i), xst(i), zst(i), mag)
                end do
                sh_disp(1, i)%kcmpnm = "Ux"; sh_disp(1, i)%cmpinc = 90.0; sh_disp(1, i)%cmpaz = 0.0 + phi
                sh_disp(2, i)%kcmpnm = "Uz"; sh_disp(2, i)%cmpinc = 90.0; sh_disp(2, i)%cmpaz = 0.0

                sh_disp(:, i)%idep = 6
            end if

            if (sw_wav_stress) then
                do j = 1, 3
                    call initialize_sac_header(sh_stress(j, i), stnm(i), stlo(i), stla(i), xst(i), zst(i), mag)
                end do

                sh_stress(1, i)%kcmpnm = "Sxx"
                sh_stress(2, i)%kcmpnm = "Szz"
                sh_stress(3, i)%kcmpnm = "Sxz"

                sh_stress(:, i)%idep = 5
            end if

            if (sw_wav_strain) then
                do j = 1, 3
                    call initialize_sac_header(sh_strain(j, i), stnm(i), stlo(i), stla(i), xst(i), zst(i), mag)
                end do

                sh_strain(1, i)%kcmpnm = "Exx"
                sh_strain(2, i)%kcmpnm = "Ezz"
                sh_strain(3, i)%kcmpnm = "Exz"

                sh_strain(:, i)%idep = 5
            end if

        end do
    end subroutine set_sac_header

    subroutine initialize_sac_header(sh, stnm0, stlo0, stla0, xst0, zst0, mag0)

        use m_daytim
        type(sac__hdr), intent(inout) :: sh
        character(*), intent(in) :: stnm0
        real, intent(in) :: stlo0, stla0, xst0, zst0, mag0
        logical, save :: first_call = .true.
        type(sac__hdr), save :: sh0
        !! ----

        if (first_call) then
            call sac__init(sh0)
            sh0%evlo = evlo
            sh0%evla = evla
            sh0%evdp = evdp  !! km-unit after SWPC5.0
            sh0%tim = exedate
            sh0%b = tbeg
            sh0%delta = ntdec_w * dt
            sh0%npts = ntw
            sh0%mag = mag0

            if (bf_mode) then
                sh0%user0 = fx0
                sh0%user2 = fz0
            else
                sh0%user0 = mxx0
                sh0%user2 = mzz0
                sh0%user4 = mxz0
            end if

            sh0%user6 = clon  !< coordinate
            sh0%user7 = clat  !< coordinate
            sh0%user8 = phi
            sh0%o = otim

            call daytim__localtime(exedate, sh0%nzyear, sh0%nzmonth, sh0%nzday, sh0%nzhour, sh0%nzmin, sh0%nzsec)
            call daytim__ymd2jul(sh0%nzyear, sh0%nzmonth, sh0%nzday, sh0%nzjday)
            sh0%nzmsec = 0

            first_call = .false.
        end if

        sh = sh0
        sh%kevnm = trim(adjustl(title(1:16)))
        sh%kstnm = trim(stnm0)
        sh%stlo = stlo0
        sh%stla = stla0
        sh%stdp = zst0 * 1000 ! in meter unit

        sh%lcalda = .false.
        sh%dist = sqrt((sx0 - xst0)**2)
        sh%az = std__rad2deg(atan2(0.0, xst0 - sx0))
        sh%baz = std__rad2deg(atan2(0.0, sx0 - xst0))

    end subroutine initialize_sac_header

    subroutine export_wav__sac(sh, dat)

        type(sac__hdr), intent(in) :: sh
        !! --
        real(SP), intent(in) :: dat(:)
        character(256) :: fn
        !! ----

        fn = trim(odir)//'/wav/'//trim(title)//'.'//trim(sh%kstnm)//'.'//trim(sh%kcmpnm)//'.sac'
        call sac__write(fn, sh, dat, .true.)

    end subroutine export_wav__sac

    subroutine export_wav__csf(nst1, ncmp, sh, dat)

        integer, intent(in) :: nst1, ncmp
        type(sac__hdr), intent(in) :: sh(ncmp, nst1)
        real(SP), intent(in) :: dat(ntw, ncmp, nst1)
        !! --
        character(5) :: cid
        character(256) :: fn
        !! ----

        write (cid, '(I5.5)') myid
        fn = trim(odir)//'/wav/'//trim(title)//'__'//cid//'__.csf'
        call csf__write(fn, nst1 * ncmp, ntw, reshape(sh, (/ncmp * nst1/)), reshape(dat, (/ntw, ncmp * nst1/)))

    end subroutine export_wav__csf

end module m_wav
