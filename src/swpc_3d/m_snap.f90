#include "../shared/m_debug.h"
module m_snap

    !! Snapshot output
    !!
    !! Copyright 2013-2025 Takuto Maeda. All rights reserved. This project is released under the MIT license.

    use iso_fortran_env, only: error_unit
    use m_std
    use m_debug
    use m_global
    use m_pwatch
    use m_fdtool
    use m_daytim
    use m_readini
    use m_geomap
    use mpi
    use netcdf

    implicit none
    private
    save

    public :: snap__setup
    public :: snap__write
    public :: snap__closefiles

    character(8), parameter :: BINARY_TYPE = "STREAMIO"
    character(8), parameter :: CODE_TYPE = "SWPC_3D "   !!< FIXED parameter for file header
    integer, parameter :: HEADER_VERSION = 6
    !! header_version history
    !! 110311 : original
    !! 2      : added binary_type (character(8))
    !! 3      : added Lz and zeta for curvilinear
    !! 4      : added w for curvilinear
    !! 5      : added coordinate rotation angle phi
    !! 6      : seism -> swpc

    type snp
        logical :: sw     ! true for output
        integer :: io     ! file I/O number (or netcdf file id)
        integer :: ionode ! output MPI node (mpi_comm_world)
        integer :: ionode_local ! output MPI node (local comm, if necessary)
        integer :: nsnp
        character(2) :: snaptype ! snapshot type
        character(2) :: coordinate
        integer :: nmed
        integer :: na1, na2 ! absorbing layer
        real    :: ds1, ds2 ! spatial grid width

        !! variables for netcdf mode
        integer :: did_x1, did_x2, did_t ! dimension id for independent vars
        integer :: vid_x1, vid_x2, vid_t ! variable  id for independent vars
        integer :: varid(10)           ! variable  id for dependent vars
        integer :: medid(10)            ! medium array id
        real    :: vmax(10), vmin(10)  ! max/min of dependent vars
        character(10) :: vname(10)
        character(10) :: vunit(10)

    end type snp

    type(snp) :: xy_ps, yz_ps, xz_ps, fs_ps, ob_ps, xy_v, yz_v, xz_v, fs_v, ob_v, xy_u, yz_u, xz_u, fs_u, ob_u

    ! switch
    integer   :: ntdec_s                            !< time step decimation factor: Snap and Waves
    integer   :: idec, jdec, kdec                   !< spatial decimation factor: x, y, z directions

    real(SP) :: z0_xy, x0_yz, y0_xz ! < danmen
    integer  :: k0_xy, i0_yz, j0_xz ! < danmen

    integer :: nxs, nys, nzs !< snapshot grid size
    real(SP), allocatable :: xsnp(:), ysnp(:), zsnp(:)

    ! I/O area in the node
    integer :: is0, is1, js0, js1, ks0, ks1

    ! derivative coefficient
    real(MP) :: r20x, r20y, r20z

    ! snapshot data format
    character(6) :: snp_format ! native or netcdf

    ! displacement snapshot buffer
    real(SP), allocatable :: buf_yz_u(:,:,:), buf_xz_u(:,:,:), buf_xy_u(:,:,:), buf_fs_u(:,:,:), buf_ob_u(:,:,:)
    real(SP), allocatable :: max_ob_v(:,:,:), max_ob_u(:,:,:), max_fs_v(:,:,:), max_fs_u(:,:,:)

    ! cross-section data MPI communicator
    integer :: mpi_comm_xz, mpi_comm_yz
    integer :: myid_xz, myid_yz
    integer :: idy_xz
    integer :: idx_yz

contains

    subroutine snap__setup(io_prm)

        integer, intent(in) :: io_prm

        integer :: i, j, k, ii, jj, kk
        integer :: err
        integer :: idum

        call pwatch__on("snap__setup")

        call readini(io_prm, 'xy_ps%sw', xy_ps%sw, .false.)
        call readini(io_prm, 'yz_ps%sw', yz_ps%sw, .false.)
        call readini(io_prm, 'xz_ps%sw', xz_ps%sw, .false.)
        call readini(io_prm, 'ob_ps%sw', ob_ps%sw, .false.)
        call readini(io_prm, 'fs_ps%sw', fs_ps%sw, .false.)
        call readini(io_prm, 'xy_v%sw', xy_v%sw, .false.)
        call readini(io_prm, 'yz_v%sw', yz_v%sw, .false.)
        call readini(io_prm, 'xz_v%sw', xz_v%sw, .false.)
        call readini(io_prm, 'ob_v%sw', ob_v%sw, .false.)
        call readini(io_prm, 'fs_v%sw', fs_v%sw, .false.)
        call readini(io_prm, 'xy_u%sw', xy_u%sw, .false.)
        call readini(io_prm, 'yz_u%sw', yz_u%sw, .false.)
        call readini(io_prm, 'xz_u%sw', xz_u%sw, .false.)
        call readini(io_prm, 'ob_u%sw', ob_u%sw, .false.)
        call readini(io_prm, 'fs_u%sw', fs_u%sw, .false.)
        call readini(io_prm, 'z0_xy', z0_xy, max(min(10.0, zend), zbeg))
        call readini(io_prm, 'x0_yz', x0_yz, max(min(0.0, xend), xbeg))
        call readini(io_prm, 'y0_xz', y0_xz, max(min(0.0, yend), ybeg))
        call readini(io_prm, 'idec', idec, 1)
        call readini(io_prm, 'jdec', jdec, 1)
        call readini(io_prm, 'kdec', kdec, 1)
        call readini(io_prm, 'ntdec_s', ntdec_s, 10)
        call readini(io_prm, 'snp_format', snp_format, 'native')

        !! snapshot size#2013-0440
        nxs = (nx + (idec / 2)) / idec
        nys = (ny + (jdec / 2)) / jdec
        nzs = (nz + (kdec / 2)) / kdec

        !! coordinate
        allocate (xsnp(nxs), ysnp(nys), zsnp(nzs))
        do i = 1, nxs
            ii = i * idec - (idec / 2)
            xsnp(i) = i2x(ii, xbeg, real(dx))
        end do
        do j = 1, nys
            jj = j * jdec - (jdec / 2)
            ysnp(j) = j2y(jj, ybeg, real(dy))
        end do
        do k = 1, nzs
            kk = k * kdec - (kdec / 2)
            zsnp(k) = k2z(kk, zbeg, real(dz))
        end do

        !! snapshot region covered by the MPI node
        is0 = ceiling((ibeg + idec / 2) / real(idec))
        is1 = floor((iend + idec / 2) / real(idec))
        js0 = ceiling((jbeg + jdec / 2) / real(jdec))
        js1 = floor((jend + jdec / 2) / real(jdec))
        ks0 = ceiling((kbeg + kdec / 2) / real(kdec))
        ks1 = floor((kend + kdec / 2) / real(kdec))

        !! snapshot location grid
        k0_xy = z2k(z0_xy, zbeg, real(dz))
        i0_yz = x2i(x0_yz, xbeg, real(dx))
        j0_xz = y2j(y0_xz, ybeg, real(dy))

        !! snapshot data communicator
        call mpi_comm_split(mpi_comm_world, idy, idx, mpi_comm_xz, err)
        call mpi_comm_split(mpi_comm_world, idx, idy, mpi_comm_yz, err)
        call mpi_comm_rank(mpi_comm_xz, myid_xz, err)
        call mpi_comm_rank(mpi_comm_yz, myid_yz, err)

        !! output node definition: cyclic
        xy_v%ionode = mod(0, nproc)
        xy_u%ionode = mod(1, nproc)
        xy_ps%ionode = mod(2, nproc)
        fs_v%ionode = mod(3, nproc)
        fs_u%ionode = mod(4, nproc)
        fs_ps%ionode = mod(5, nproc)
        ob_v%ionode = mod(6, nproc)
        ob_u%ionode = mod(7, nproc)
        ob_ps%ionode = mod(8, nproc)

        call global__getnode(1, j0_xz, idum, idy_xz)
        call global__getnode(i0_yz, 1, idx_yz, idum)

        xz_v%ionode_local = mod(0, nproc_x)
        xz_u%ionode_local = mod(1, nproc_x)
        xz_ps%ionode_local = mod(2, nproc_x)
        xz_v%ionode = itbl(xz_v%ionode_local, idy_xz)
        xz_u%ionode = itbl(xz_u%ionode_local, idy_xz)
        xz_ps%ionode = itbl(xz_ps%ionode_local, idy_xz)

        yz_v%ionode_local = mod(0, nproc_y)
        yz_u%ionode_local = mod(1, nproc_y)
        yz_ps%ionode_local = mod(2, nproc_y)
        yz_v%ionode = itbl(idx_yz, yz_v%ionode_local)
        yz_u%ionode = itbl(idx_yz, yz_u%ionode_local)
        yz_ps%ionode = itbl(idx_yz, yz_ps%ionode_local)

       !! snapshot settings
        yz_ps%nsnp = 4; yz_v%nsnp = 3; yz_u%nsnp = 3; 
        xz_ps%nsnp = 4; xz_v%nsnp = 3; xz_u%nsnp = 3; 
        xy_ps%nsnp = 4; xy_v%nsnp = 3; xy_u%nsnp = 3; 
        ob_ps%nsnp = 4; ob_v%nsnp = 3; ob_u%nsnp = 3; 
        fs_ps%nsnp = 4; fs_v%nsnp = 3; fs_u%nsnp = 3; 
        ! horizontal medium continas topography, lon & lat
        yz_ps%nmed = 3; yz_v%nmed = 3; yz_u%nmed = 3; 
        xz_ps%nmed = 3; xz_v%nmed = 3; xz_u%nmed = 3; 
        xy_ps%nmed = 6; xy_v%nmed = 6; xy_u%nmed = 6; 
        ob_ps%nmed = 6; ob_v%nmed = 6; ob_u%nmed = 6; 
        fs_ps%nmed = 6; fs_v%nmed = 6; fs_u%nmed = 6; 
        yz_ps%snaptype = 'ps'; yz_v%snaptype = 'v3'; yz_u%snaptype = 'u3'; 
        xz_ps%snaptype = 'ps'; xz_v%snaptype = 'v3'; xz_u%snaptype = 'u3'; 
        xy_ps%snaptype = 'ps'; xy_v%snaptype = 'v3'; xy_u%snaptype = 'u3'; 
        ob_ps%snaptype = 'ps'; ob_v%snaptype = 'v3'; ob_u%snaptype = 'u3'; 
        fs_ps%snaptype = 'ps'; fs_v%snaptype = 'v3'; fs_u%snaptype = 'u3'; 
        yz_ps%vname(1) = 'div'; yz_ps%vname(2) = 'rot_x'; yz_ps%vname(3) = 'rot_y'; yz_ps%vname(4) = 'rot_z'
        xz_ps%vname(1) = 'div'; xz_ps%vname(2) = 'rot_x'; xz_ps%vname(3) = 'rot_y'; xz_ps%vname(4) = 'rot_z'
        xy_ps%vname(1) = 'div'; xy_ps%vname(2) = 'rot_x'; xy_ps%vname(3) = 'rot_y'; xy_ps%vname(4) = 'rot_z'
        fs_ps%vname(1) = 'div'; fs_ps%vname(2) = 'rot_x'; fs_ps%vname(3) = 'rot_y'; fs_ps%vname(4) = 'rot_z'
        ob_ps%vname(1) = 'div'; ob_ps%vname(2) = 'rot_x'; ob_ps%vname(3) = 'rot_y'; ob_ps%vname(4) = 'rot_z'

        yz_v%vname(1) = 'Vx'; yz_v%vname(2) = 'Vy'; yz_v%vname(3) = 'Vz'
        xz_v%vname(1) = 'Vx'; xz_v%vname(2) = 'Vy'; xz_v%vname(3) = 'Vz'
        xy_v%vname(1) = 'Vx'; xy_v%vname(2) = 'Vy'; xy_v%vname(3) = 'Vz'
        fs_v%vname(1) = 'Vx'; fs_v%vname(2) = 'Vy'; fs_v%vname(3) = 'Vz'
        ob_v%vname(1) = 'Vx'; ob_v%vname(2) = 'Vy'; ob_v%vname(3) = 'Vz'

        yz_u%vname(1) = 'Ux'; yz_u%vname(2) = 'Uy'; yz_u%vname(3) = 'Uz'
        xz_u%vname(1) = 'Ux'; xz_u%vname(2) = 'Uy'; xz_u%vname(3) = 'Uz'
        xy_u%vname(1) = 'Ux'; xy_u%vname(2) = 'Uy'; xy_u%vname(3) = 'Uz'
        fs_u%vname(1) = 'Ux'; fs_u%vname(2) = 'Uy'; fs_u%vname(3) = 'Uz'
        ob_u%vname(1) = 'Ux'; ob_u%vname(2) = 'Uy'; ob_u%vname(3) = 'Uz'

        yz_ps%vunit(1) = '1/s'; yz_ps%vunit(2) = '1/s'; yz_ps%vunit(3) = '1/s'; yz_ps%vunit(4) = '1/s'
        xz_ps%vunit(1) = '1/s'; xz_ps%vunit(2) = '1/s'; xz_ps%vunit(3) = '1/s'; xz_ps%vunit(4) = '1/s'
        xy_ps%vunit(1) = '1/s'; xy_ps%vunit(2) = '1/s'; xy_ps%vunit(3) = '1/s'; xy_ps%vunit(4) = '1/s'
        fs_ps%vunit(1) = '1/s'; fs_ps%vunit(2) = '1/s'; fs_ps%vunit(3) = '1/s'; fs_ps%vunit(4) = '1/s'
        ob_ps%vunit(1) = '1/s'; ob_ps%vunit(2) = '1/s'; ob_ps%vunit(3) = '1/s'; ob_ps%vunit(4) = '1/s'

        yz_v%vunit(1) = 'm/s'; yz_v%vunit(2) = 'm/s'; yz_v%vunit(3) = 'm/s'
        xz_v%vunit(1) = 'm/s'; xz_v%vunit(2) = 'm/s'; xz_v%vunit(3) = 'm/s'
        xy_v%vunit(1) = 'm/s'; xy_v%vunit(2) = 'm/s'; xy_v%vunit(3) = 'm/s'
        fs_v%vunit(1) = 'm/s'; fs_v%vunit(2) = 'm/s'; fs_v%vunit(3) = 'm/s'
        ob_v%vunit(1) = 'm/s'; ob_v%vunit(2) = 'm/s'; ob_v%vunit(3) = 'm/s'

        yz_u%vunit(1) = 'm'; yz_u%vunit(2) = 'm'; yz_u%vunit(3) = 'm'
        xz_u%vunit(1) = 'm'; xz_u%vunit(2) = 'm'; xz_u%vunit(3) = 'm'
        xy_u%vunit(1) = 'm'; xy_u%vunit(2) = 'm'; xy_u%vunit(3) = 'm'
        fs_u%vunit(1) = 'm'; fs_u%vunit(2) = 'm'; fs_u%vunit(3) = 'm'
        ob_u%vunit(1) = 'm'; ob_u%vunit(2) = 'm'; ob_u%vunit(3) = 'm'

        yz_ps%coordinate = 'yz'; yz_v%coordinate = 'yz'; yz_u%coordinate = 'yz'; 
        xz_ps%coordinate = 'xz'; xz_v%coordinate = 'xz'; xz_u%coordinate = 'xz'; 
        xy_ps%coordinate = 'xy'; xy_v%coordinate = 'xy'; xy_u%coordinate = 'xy'; 
        ob_ps%coordinate = 'ob'; ob_v%coordinate = 'ob'; ob_u%coordinate = 'ob'; 
        fs_ps%coordinate = 'fs'; fs_v%coordinate = 'fs'; fs_u%coordinate = 'fs'; 
        !! output settings
        if (snp_format == 'native') then

            if (yz_ps%sw) call newfile_yz(trim(odir)//'/'//trim(title)//'.3d.yz.ps.snp', yz_ps)
            if (xz_ps%sw) call newfile_xz(trim(odir)//'/'//trim(title)//'.3d.xz.ps.snp', xz_ps)
            if (xy_ps%sw) call newfile_xy(trim(odir)//'/'//trim(title)//'.3d.xy.ps.snp', xy_ps)
            if (fs_ps%sw) call newfile_xy(trim(odir)//'/'//trim(title)//'.3d.fs.ps.snp', fs_ps)
            if (ob_ps%sw) call newfile_xy(trim(odir)//'/'//trim(title)//'.3d.ob.ps.snp', ob_ps)

            if (yz_v%sw) call newfile_yz(trim(odir)//'/'//trim(title)//'.3d.yz.v.snp', yz_v)
            if (xz_v%sw) call newfile_xz(trim(odir)//'/'//trim(title)//'.3d.xz.v.snp', xz_v)
            if (xy_v%sw) call newfile_xy(trim(odir)//'/'//trim(title)//'.3d.xy.v.snp', xy_v)
            if (fs_v%sw) call newfile_xy(trim(odir)//'/'//trim(title)//'.3d.fs.v.snp', fs_v)
            if (ob_v%sw) call newfile_xy(trim(odir)//'/'//trim(title)//'.3d.ob.v.snp', ob_v)

            if (yz_u%sw) call newfile_yz(trim(odir)//'/'//trim(title)//'.3d.yz.u.snp', yz_u)
            if (xz_u%sw) call newfile_xz(trim(odir)//'/'//trim(title)//'.3d.xz.u.snp', xz_u)
            if (xy_u%sw) call newfile_xy(trim(odir)//'/'//trim(title)//'.3d.xy.u.snp', xy_u)
            if (fs_u%sw) call newfile_xy(trim(odir)//'/'//trim(title)//'.3d.fs.u.snp', fs_u)
            if (ob_u%sw) call newfile_xy(trim(odir)//'/'//trim(title)//'.3d.ob.u.snp', ob_u)

        else

            if (yz_ps%sw) call newfile_yz_nc(trim(odir)//'/'//trim(title)//'.3d.yz.ps.nc', yz_ps)
            if (xz_ps%sw) call newfile_xz_nc(trim(odir)//'/'//trim(title)//'.3d.xz.ps.nc', xz_ps)
            if (xy_ps%sw) call newfile_xy_nc(trim(odir)//'/'//trim(title)//'.3d.xy.ps.nc', xy_ps)
            if (fs_ps%sw) call newfile_xy_nc(trim(odir)//'/'//trim(title)//'.3d.fs.ps.nc', fs_ps)
            if (ob_ps%sw) call newfile_xy_nc(trim(odir)//'/'//trim(title)//'.3d.ob.ps.nc', ob_ps)

            if (yz_v%sw) call newfile_yz_nc(trim(odir)//'/'//trim(title)//'.3d.yz.v.nc', yz_v)
            if (xz_v%sw) call newfile_xz_nc(trim(odir)//'/'//trim(title)//'.3d.xz.v.nc', xz_v)
            if (xy_v%sw) call newfile_xy_nc(trim(odir)//'/'//trim(title)//'.3d.xy.v.nc', xy_v)
            if (fs_v%sw) call newfile_xy_nc(trim(odir)//'/'//trim(title)//'.3d.fs.v.nc', fs_v)
            if (ob_v%sw) call newfile_xy_nc(trim(odir)//'/'//trim(title)//'.3d.ob.v.nc', ob_v)

            if (yz_u%sw) call newfile_yz_nc(trim(odir)//'/'//trim(title)//'.3d.yz.u.nc', yz_u)
            if (xz_u%sw) call newfile_xz_nc(trim(odir)//'/'//trim(title)//'.3d.xz.u.nc', xz_u)
            if (xy_u%sw) call newfile_xy_nc(trim(odir)//'/'//trim(title)//'.3d.xy.u.nc', xy_u)
            if (fs_u%sw) call newfile_xy_nc(trim(odir)//'/'//trim(title)//'.3d.fs.u.nc', fs_u)
            if (ob_u%sw) call newfile_xy_nc(trim(odir)//'/'//trim(title)//'.3d.ob.u.nc', ob_u)

        end if

        !! for taking derivatives
        r20x = 1.0_MP / dx
        r20y = 1.0_MP / dy
        r20z = 1.0_MP / dz

        !! displacement snapshot buffers
        allocate (buf_yz_u(nys, nzs, 3), source=0.0)
        allocate (buf_xz_u(nxs, nzs, 3), source=0.0)
        allocate (buf_xy_u(nxs, nys, 3), source=0.0)
        allocate (buf_fs_u(nxs, nys, 3), source=0.0)
        allocate (buf_ob_u(nxs, nys, 3), source=0.0)

        !! maximum amplitude buffers
        allocate (max_ob_v(nxs, nys, 3), source=0.0)
        allocate (max_ob_u(nxs, nys, 3), source=0.0)
        allocate (max_fs_v(nxs, nys, 3), source=0.0)
        allocate (max_fs_u(nxs, nys, 3), source=0.0)

        !$acc enter data copyin(buf_yz_u, buf_xz_u, buf_xy_u, buf_fs_u, buf_ob_u, &
        !$acc                   max_ob_v, max_ob_u, max_fs_v, max_fs_u, i0_yz, j0_xz, k0_xy)

        call mpi_barrier(mpi_comm_world, err)

        call pwatch__off("snap__setup")

    end subroutine snap__setup

    subroutine newfile_yz(fname, hdr)

        !! Open new file and write header information, and medium parameters for YZ-cross section

        character(*), intent(in)  :: fname
        type(snp), intent(inout) :: hdr

        real(SP) :: buf(nys, nzs, 3)
        integer :: j, k, kk, ii, jj

        ! absorber thickness setting
        hdr%na1 = na / jdec
        hdr%na2 = na / kdec

        if (myid == hdr%ionode) then

            open (newunit=hdr%io, file=trim(fname), access='stream', action='write', status='replace', form='unformatted')
            call write_snp_header(hdr, nys, nzs, ysnp(1:nys), zsnp(1:nzs))

        end if

        buf = 0.0

        ii = i0_yz
        if (ibeg <= ii .and. ii <= iend) then
            do j = js0, js1
                do k = ks0, ks1
                    kk = k * kdec - kdec / 2
                    jj = j * jdec - jdec / 2

                    buf(j, k, 1) = rho(kk, ii, jj)
                    buf(j, k, 2) = lam(kk, ii, jj)
                    buf(j, k, 3) = mu(kk, ii, jj)

                end do
            end do
        end if

        call write_reduce_array2d_r(nys, nzs, hdr%ionode, hdr%io, buf(:, :, 1))
        call write_reduce_array2d_r(nys, nzs, hdr%ionode, hdr%io, buf(:, :, 2))
        call write_reduce_array2d_r(nys, nzs, hdr%ionode, hdr%io, buf(:, :, 3))

    end subroutine newfile_yz

    subroutine newfile_xz(fname, hdr)

        !! Open new file and write header information, and medium parameters for XY-cross section

        character(*), intent(in)  :: fname
        type(snp), intent(inout) :: hdr

        integer :: i, k, kk, ii, jj
        real(SP) :: buf(nxs, nzs, 3)

        ! absorber thickness setting
        hdr%na1 = na / idec
        hdr%na2 = na / kdec

        if (myid == hdr%ionode) then

            open (newunit=hdr%io, file=trim(fname), access='stream', action='write', status='replace', form='unformatted')
            call write_snp_header(hdr, nxs, nzs, xsnp(1:nxs), zsnp(1:nzs))

        end if

        buf = 0.0

        jj = j0_xz
        if (jbeg <= jj .and. jj <= jend) then
            do i = is0, is1
                do k = ks0, ks1

                    ii = i * idec - idec / 2
                    kk = k * kdec - kdec / 2

                    buf(i, k, 1) = rho(kk, ii, jj)
                    buf(i, k, 2) = lam(kk, ii, jj)
                    buf(i, k, 3) = mu(kk, ii, jj)

                end do
            end do
        end if

        call write_reduce_array2d_r(nxs, nzs, hdr%ionode, hdr%io, buf(:, :, 1))
        call write_reduce_array2d_r(nxs, nzs, hdr%ionode, hdr%io, buf(:, :, 2))
        call write_reduce_array2d_r(nxs, nzs, hdr%ionode, hdr%io, buf(:, :, 3))

    end subroutine newfile_xz

    subroutine newfile_xy(fname, hdr)

        !! Open new file and write header information, and medium parameters for XY-cross section

        character(*), intent(in) :: fname
        type(snp), intent(inout) :: hdr

        integer :: i, j, kk, ii, jj
        real(SP) :: buf(nxs, nys, 3)

        ! absorber thickness setting
        hdr%na1 = na / idec
        hdr%na2 = na / jdec

        if (myid == hdr%ionode) then

            open (newunit=hdr%io, file=trim(fname), access='stream', action='write', status='replace', form='unformatted')
            call write_snp_header(hdr, nxs, nys, xsnp(1:nxs), ysnp(1:nys))

        end if

        buf = 0.0

        do j = js0, js1
            do i = is0, is1

                ii = i * idec - idec / 2
                jj = j * jdec - jdec / 2

                if (hdr%coordinate == 'fs') then
                    kk = kfs(ii, jj) + 1
                else if (hdr%coordinate == 'ob') then
                    kk = kob(ii, jj) + 1
                else
                    kk = k0_xy
                end if

                buf(i, j, 1) = rho(kk, ii, jj)
                buf(i, j, 2) = lam(kk, ii, jj)
                buf(i, j, 3) = mu(kk, ii, jj)

            end do
        end do

        call write_reduce_array2d_r(nxs, nys, hdr%ionode, hdr%io, buf(:, :, 1))
        call write_reduce_array2d_r(nxs, nys, hdr%ionode, hdr%io, buf(:, :, 2))
        call write_reduce_array2d_r(nxs, nys, hdr%ionode, hdr%io, buf(:, :, 3))

        do j = js0, js1
            do i = is0, is1

                ii = i * idec - idec / 2
                jj = j * jdec - jdec / 2
                !! topography data
                buf(i, j, 1) = -bddep(ii, jj, 0) * 1000 ! positive upward, in unit of [m]

            end do
        end do
        call write_reduce_array2d_r(nxs, nys, hdr%ionode, hdr%io, buf(:, :, 1))

    end subroutine newfile_xy

    subroutine newfile_yz_nc(fname, hdr)

        !! create snapshot file in netcdf format

        character(*), intent(in)    :: fname
        type(snp), intent(inout) :: hdr

        real(SP), allocatable :: sbuf(:), rbuf1(:), rbuf2(:), rbuf3(:), buf(:, :, :)
        integer :: j, k, ii, jj, kk, err

        if (idx /= idx_yz) return

        if (myid == hdr%ionode) then

            !! initialize
            hdr%vmax = 0.0
            hdr%vmin = 0.0
            hdr%na1 = na / jdec
            hdr%na2 = na / kdec
            hdr%ds1 = jdec * real(dy)
            hdr%ds2 = kdec * real(dz)

            call nc_chk(nf90_create(trim(fname), NF90_CLOBBER, hdr%io))
            call write_nc_header(hdr, nys, nzs, ysnp, zsnp)

        end if

        allocate (buf(nys, nzs, 3), source=0.0)

        ii = i0_yz
        if (ibeg <= ii .and. ii <= iend) then
            do j = js0, js1
                do k = ks0, ks1
                    kk = k * kdec - kdec / 2
                    jj = j * jdec - jdec / 2

                    buf(j, k, 1) = rho(kk, ii, jj)
                    buf(j, k, 2) = lam(kk, ii, jj)
                    buf(j, k, 3) = mu(kk, ii, jj)

                end do
            end do
        end if

        !! medium
        allocate(sbuf(nys * nzs),  source=0.0)
        allocate(rbuf1(nys * nzs), source=0.0) 
        allocate(rbuf2(nys * nzs), source=0.0) 
        allocate(rbuf3(nys * nzs), source=0.0)
        sbuf = reshape(buf(:, :, 1), shape(sbuf))
        call mpi_reduce(sbuf, rbuf1, nys * nzs, MPI_REAL, MPI_SUM, hdr%ionode_local, mpi_comm_yz, err)
        sbuf = reshape(buf(:, :, 2), shape(sbuf))
        call mpi_reduce(sbuf, rbuf2, nys * nzs, MPI_REAL, MPI_SUM, hdr%ionode_local, mpi_comm_yz, err)
        sbuf = reshape(buf(:, :, 3), shape(sbuf))
        call mpi_reduce(sbuf, rbuf3, nys * nzs, MPI_REAL, MPI_SUM, hdr%ionode_local, mpi_comm_yz, err)

        if (myid == hdr%ionode) then
            call nc_chk(nf90_put_att(hdr%io, hdr%medid(1), 'actual_range', (/minval(rbuf1), maxval(rbuf1)/)))
            call nc_chk(nf90_put_att(hdr%io, hdr%medid(2), 'actual_range', (/minval(rbuf2), maxval(rbuf2)/)))
            call nc_chk(nf90_put_att(hdr%io, hdr%medid(3), 'actual_range', (/minval(rbuf3), maxval(rbuf3)/)))

            call nc_chk(nf90_enddef(hdr%io))

            call nc_chk(nf90_put_var(hdr%io, hdr%vid_x1, ysnp))
            call nc_chk(nf90_put_var(hdr%io, hdr%vid_x2, zsnp))
            call nc_chk(nf90_put_var(hdr%io, hdr%medid(1), reshape(rbuf1, shape(buf(:, :, 1)))))
            call nc_chk(nf90_put_var(hdr%io, hdr%medid(2), reshape(rbuf2, shape(buf(:, :, 2)))))
            call nc_chk(nf90_put_var(hdr%io, hdr%medid(3), reshape(rbuf3, shape(buf(:, :, 3)))))

        end if

        deallocate (sbuf, rbuf1, rbuf2, rbuf3, buf)

    end subroutine newfile_yz_nc

    subroutine newfile_xz_nc(fname, hdr)

        !! create snapshot file in netcdf format

        character(*), intent(in)    :: fname
        type(snp), intent(inout) :: hdr

        real(SP), allocatable :: sbuf(:), rbuf1(:), rbuf2(:), rbuf3(:), buf(:, :, :)
        integer :: i, k, ii, jj, kk, err

        if (idy /= idy_xz) return

        if (myid == hdr%ionode) then

          !! initialize
            hdr%vmax = 0.0
            hdr%vmin = 0.0
            hdr%na1 = na / idec
            hdr%na2 = na / kdec
            hdr%ds1 = idec * dx
            hdr%ds2 = kdec * dz

            call nc_chk(nf90_create(trim(fname), NF90_CLOBBER, hdr%io))
            call write_nc_header(hdr, nxs, nzs, xsnp, zsnp)

        end if

        allocate (buf(nxs, nzs, 3), source=0.0)

        jj = j0_xz

        if (jbeg <= jj .and. jj <= jend) then
            do i = is0, is1
                do k = ks0, ks1

                    ii = i * idec - idec / 2
                    kk = k * kdec - kdec / 2

                    buf(i, k, 1) = rho(kk, ii, jj)
                    buf(i, k, 2) = lam(kk, ii, jj)
                    buf(i, k, 3) = mu(kk, ii, jj)

                end do
            end do
        end if

        !! medium
        allocate(sbuf(nxs * nzs), source=0.0)
        allocate(rbuf1(nxs * nzs),source=0.0) 
        allocate(rbuf2(nxs * nzs),source=0.0) 
        allocate(rbuf3(nxs * nzs),source=0.0)
        sbuf = reshape(buf(:, :, 1), shape(sbuf))
        call mpi_reduce(sbuf, rbuf1, nxs * nzs, MPI_REAL, MPI_SUM, hdr%ionode_local, mpi_comm_xz, err)
        sbuf = reshape(buf(:, :, 2), shape(sbuf))
        call mpi_reduce(sbuf, rbuf2, nxs * nzs, MPI_REAL, MPI_SUM, hdr%ionode_local, mpi_comm_xz, err)
        sbuf = reshape(buf(:, :, 3), shape(sbuf))
        call mpi_reduce(sbuf, rbuf3, nxs * nzs, MPI_REAL, MPI_SUM, hdr%ionode_local, mpi_comm_xz, err)

        if (myid == hdr%ionode) then
            call nc_chk(nf90_put_att(hdr%io, hdr%medid(1), 'actual_range', (/minval(rbuf1), maxval(rbuf1)/)))
            call nc_chk(nf90_put_att(hdr%io, hdr%medid(2), 'actual_range', (/minval(rbuf2), maxval(rbuf2)/)))
            call nc_chk(nf90_put_att(hdr%io, hdr%medid(3), 'actual_range', (/minval(rbuf3), maxval(rbuf3)/)))

            call nc_chk(nf90_enddef(hdr%io))

            call nc_chk(nf90_put_var(hdr%io, hdr%vid_x1, xsnp))
            call nc_chk(nf90_put_var(hdr%io, hdr%vid_x2, zsnp))
            call nc_chk(nf90_put_var(hdr%io, hdr%medid(1), reshape(rbuf1, shape(buf(:, :, 1)))))
            call nc_chk(nf90_put_var(hdr%io, hdr%medid(2), reshape(rbuf2, shape(buf(:, :, 1)))))
            call nc_chk(nf90_put_var(hdr%io, hdr%medid(3), reshape(rbuf3, shape(buf(:, :, 1)))))

        end if

        deallocate (sbuf, rbuf1, rbuf2, rbuf3, buf)

    end subroutine newfile_xz_nc

    subroutine write_nc_header(hdr, ns1, ns2, xs1, xs2)

        !! write common netcdf header

        type(snp), intent(inout) :: hdr
        integer, intent(in) :: ns1, ns2
        real(SP), intent(in) :: xs1(ns1), xs2(ns2)
        integer :: i

        if (hdr%coordinate == 'yz') then
            call nc_chk(nf90_def_dim(hdr%io, 'y', nys, hdr%did_x1))
            call nc_chk(nf90_def_dim(hdr%io, 'z', nzs, hdr%did_x2))
            call nc_chk(nf90_def_var(hdr%io, 'y', NF90_REAL, hdr%did_x1, hdr%vid_x1))
            call nc_chk(nf90_def_var(hdr%io, 'z', NF90_REAL, hdr%did_x2, hdr%vid_x2))
            call nc_chk(nf90_put_att(hdr%io, hdr%vid_x1, 'long_name', 'y'))
            call nc_chk(nf90_put_att(hdr%io, hdr%vid_x2, 'long_name', 'z'))
        else if (hdr%coordinate == 'xz') then
            call nc_chk(nf90_def_dim(hdr%io, 'x', nxs, hdr%did_x1))
            call nc_chk(nf90_def_dim(hdr%io, 'z', nzs, hdr%did_x2))
            call nc_chk(nf90_def_var(hdr%io, 'x', NF90_REAL, hdr%did_x1, hdr%vid_x1))
            call nc_chk(nf90_def_var(hdr%io, 'z', NF90_REAL, hdr%did_x2, hdr%vid_x2))
            call nc_chk(nf90_put_att(hdr%io, hdr%vid_x1, 'long_name', 'x'))
            call nc_chk(nf90_put_att(hdr%io, hdr%vid_x2, 'long_name', 'z'))
        else
            call nc_chk(nf90_def_dim(hdr%io, 'x', nxs, hdr%did_x1))
            call nc_chk(nf90_def_dim(hdr%io, 'y', nys, hdr%did_x2))
            call nc_chk(nf90_def_var(hdr%io, 'x', NF90_REAL, hdr%did_x1, hdr%vid_x1))
            call nc_chk(nf90_def_var(hdr%io, 'y', NF90_REAL, hdr%did_x2, hdr%vid_x2))
            call nc_chk(nf90_put_att(hdr%io, hdr%vid_x1, 'long_name', 'x'))
            call nc_chk(nf90_put_att(hdr%io, hdr%vid_x2, 'long_name', 'y'))
        end if

        call nc_chk(nf90_def_dim(hdr%io, 't', NF90_UNLIMITED, hdr%did_t))
        call nc_chk(nf90_def_var(hdr%io, 't', NF90_REAL, hdr%did_t, hdr%vid_t))

        !! medium
        call nc_chk(nf90_def_var(hdr%io, 'rho', NF90_REAL, (/hdr%did_x1, hdr%did_x2/), hdr%medid(1)))
        call nc_chk(nf90_def_var(hdr%io, 'lambda', NF90_REAL, (/hdr%did_x1, hdr%did_x2/), hdr%medid(2)))
        call nc_chk(nf90_def_var(hdr%io, 'mu', NF90_REAL, (/hdr%did_x1, hdr%did_x2/), hdr%medid(3)))

        !! special for horizontal snapshot
        if (hdr%coordinate == 'ob' .or. hdr%coordinate == 'fs' .or. hdr%coordinate == 'xy') then
            call nc_chk(nf90_def_var(hdr%io, 'topo', NF90_REAL, (/hdr%did_x1, hdr%did_x2/), hdr%medid(4)))
            call nc_chk(nf90_def_var(hdr%io, 'lon', NF90_REAL, (/hdr%did_x1, hdr%did_x2/), hdr%medid(5)))
            call nc_chk(nf90_def_var(hdr%io, 'lat', NF90_REAL, (/hdr%did_x1, hdr%did_x2/), hdr%medid(6)))
            call nc_chk(nf90_put_att(hdr%io, hdr%medid(4), 'long_name', 'Topography'))
            call nc_chk(nf90_put_att(hdr%io, hdr%medid(5), 'long_name', 'Longitude'))
            call nc_chk(nf90_put_att(hdr%io, hdr%medid(6), 'long_name', 'Latitude'))
            call nc_chk(nf90_put_att(hdr%io, hdr%medid(4), 'units', 'm'))
            call nc_chk(nf90_put_att(hdr%io, hdr%medid(5), 'units', 'degrees_east'))
            call nc_chk(nf90_put_att(hdr%io, hdr%medid(6), 'units', 'degrees_north'))

            !! for make tools recognize map projection
            call nc_chk(nf90_put_att(hdr%io, hdr%vid_x1, 'standard_name', 'projection_x_coordinate'))
            call nc_chk(nf90_put_att(hdr%io, hdr%vid_x2, 'standard_name', 'projection_y_coordinate'))
            do i = 1, 4
                call nc_chk(nf90_put_att(hdr%io, hdr%medid(i), 'coordinates', 'lat lon'))
            end do
        end if

        !! variables
        do i = 1, hdr%nsnp
            call nc_chk(nf90_def_var(hdr%io, trim(hdr%vname(i)), NF90_REAL, (/hdr%did_x1, hdr%did_x2, hdr%did_t/), hdr%varid(i)))
            call nc_chk(nf90_put_att(hdr%io, hdr%varid(i), 'long_name', trim(hdr%vname(i))))
            call nc_chk(nf90_put_att(hdr%io, hdr%varid(i), 'units', trim(hdr%vunit(i))))

            if (hdr%coordinate == 'ob' .or. hdr%coordinate == 'fs' .or. hdr%coordinate == 'xy') then
                call nc_chk(nf90_put_att(hdr%io, hdr%varid(i), 'coordinates', 'lat lon'))
            end if

        end do

        !! global attribute
        call nc_chk(nf90_put_att(hdr%io, NF90_GLOBAL, 'generated_by', 'SWPC'))
        call nc_chk(nf90_put_att(hdr%io, NF90_GLOBAL, 'codetype', CODE_TYPE))
        call nc_chk(nf90_put_att(hdr%io, NF90_GLOBAL, 'hdrver', HEADER_VERSION))
        call nc_chk(nf90_put_att(hdr%io, NF90_GLOBAL, 'title', trim(title)))
        call nc_chk(nf90_put_att(hdr%io, NF90_GLOBAL, 'exedate', exedate))
        call nc_chk(nf90_put_att(hdr%io, NF90_GLOBAL, 'ns1', ns1))
        call nc_chk(nf90_put_att(hdr%io, NF90_GLOBAL, 'ns2', ns2))
        call nc_chk(nf90_put_att(hdr%io, NF90_GLOBAL, 'beg1', xs1(1)))
        call nc_chk(nf90_put_att(hdr%io, NF90_GLOBAL, 'beg2', xs2(1)))
        call nc_chk(nf90_put_att(hdr%io, NF90_GLOBAL, 'na1', hdr%na1))
        call nc_chk(nf90_put_att(hdr%io, NF90_GLOBAL, 'na2', hdr%na2))
        call nc_chk(nf90_put_att(hdr%io, NF90_GLOBAL, 'ds1', hdr%ds1))
        call nc_chk(nf90_put_att(hdr%io, NF90_GLOBAL, 'ds2', hdr%ds2))
        call nc_chk(nf90_put_att(hdr%io, NF90_GLOBAL, 'nmed', hdr%nmed))
        call nc_chk(nf90_put_att(hdr%io, NF90_GLOBAL, 'nsnp', hdr%nsnp))
        call nc_chk(nf90_put_att(hdr%io, NF90_GLOBAL, 'coordinate', hdr%coordinate))
        call nc_chk(nf90_put_att(hdr%io, NF90_GLOBAL, 'datatype', hdr%snaptype))
        call nc_chk(nf90_put_att(hdr%io, NF90_GLOBAL, 'dt', dt * ntdec_s))

        call nc_chk(nf90_put_att(hdr%io, NF90_GLOBAL, 'evlo', evlo))
        call nc_chk(nf90_put_att(hdr%io, NF90_GLOBAL, 'evla', evla))
        call nc_chk(nf90_put_att(hdr%io, NF90_GLOBAL, 'evdp', evdp))
        call nc_chk(nf90_put_att(hdr%io, NF90_GLOBAL, 'evx', sx0))
        call nc_chk(nf90_put_att(hdr%io, NF90_GLOBAL, 'evy', sy0))

        !! variable attributes
        call nc_chk(nf90_put_att(hdr%io, hdr%vid_t, 'long_name', 't'))
        call nc_chk(nf90_put_att(hdr%io, hdr%vid_x1, 'units', 'km'))
        call nc_chk(nf90_put_att(hdr%io, hdr%vid_x2, 'units', 'km'))
        call nc_chk(nf90_put_att(hdr%io, hdr%vid_t, 'units', 's'))

        call nc_chk(nf90_put_att(hdr%io, hdr%medid(1), 'long_name', 'rho'))
        call nc_chk(nf90_put_att(hdr%io, hdr%medid(2), 'long_name', 'lambda'))
        call nc_chk(nf90_put_att(hdr%io, hdr%medid(3), 'long_name', 'mu'))
        call nc_chk(nf90_put_att(hdr%io, hdr%medid(1), 'units', '10^3 kg/cm^3'))
        call nc_chk(nf90_put_att(hdr%io, hdr%medid(2), 'units', '10^9 Pa'))
        call nc_chk(nf90_put_att(hdr%io, hdr%medid(3), 'units', '10^9 Pa'))

        call nc_chk(nf90_put_att(hdr%io, hdr%vid_x1, 'actual_range', (/xs1(1), xs1(ns1)/)))
        call nc_chk(nf90_put_att(hdr%io, hdr%vid_x2, 'actual_range', (/xs2(1), xs2(ns2)/)))

        call nc_chk(nf90_put_att(hdr%io, NF90_GLOBAL, 'clon', clon))
        call nc_chk(nf90_put_att(hdr%io, NF90_GLOBAL, 'clat', clat))
        call nc_chk(nf90_put_att(hdr%io, NF90_GLOBAL, 'phi', phi))

    end subroutine write_nc_header

    subroutine newfile_xy_nc(fname, hdr)

        !! create snapshot file in netcdf format

        character(*), intent(in)    :: fname
        type(snp), intent(inout) :: hdr

        real(SP), allocatable :: sbuf(:), rbuf1(:), rbuf2(:), rbuf3(:), rbuf4(:), rbuf5(:), rbuf6(:), buf(:, :, :)
        integer :: i, j, ii, jj, kk, err

        if (myid == hdr%ionode) then

            !! initialize
            hdr%vmax = 0.0
            hdr%vmin = 0.0
            hdr%na1 = na / idec
            hdr%na2 = na / jdec
            hdr%ds1 = idec * dx
            hdr%ds2 = jdec * dy

            call nc_chk(nf90_create(trim(fname), NF90_CLOBBER, hdr%io))
            call write_nc_header(hdr, nxs, nys, xsnp, ysnp)

        end if

        allocate (buf(nxs, nys, 6), source=0.0)

        do j = js0, js1
            do i = is0, is1

                ii = i * idec - idec / 2
                jj = j * jdec - jdec / 2

                if (hdr%coordinate == 'fs') then
                    kk = kfs(ii, jj) + 1
                else if (hdr%coordinate == 'ob') then
                    kk = kob(ii, jj) + 1
                else
                    kk = k0_xy
                end if

                buf(i, j, 1) = rho(kk, ii, jj)
                buf(i, j, 2) = lam(kk, ii, jj)
                buf(i, j, 3) = mu(kk, ii, jj)
                buf(i, j, 4) = -bddep(ii, jj, 0) * 1000 ! positive upward, in meter unit

                ! longitude & latitude table
                call geomap__c2g(xsnp(i), ysnp(j), clon, clat, phi, buf(i, j, 5), buf(i, j, 6))
            end do
        end do

        !! medium
        allocate(sbuf (nxs * nys), source=0.0)
        allocate(rbuf1(nxs * nys), source=0.0)
        allocate(rbuf2(nxs * nys), source=0.0)
        allocate(rbuf3(nxs * nys), source=0.0)
        allocate(rbuf4(nxs * nys), source=0.0)
        allocate(rbuf5(nxs * nys), source=0.0)
        allocate(rbuf6(nxs * nys), source=0.0)

        sbuf = reshape(buf(:, :, 1), shape(sbuf))
        call mpi_reduce(sbuf, rbuf1, nxs * nys, MPI_REAL, MPI_SUM, hdr%ionode, mpi_comm_world, err)
        sbuf = reshape(buf(:, :, 2), shape(sbuf))
        call mpi_reduce(sbuf, rbuf2, nxs * nys, MPI_REAL, MPI_SUM, hdr%ionode, mpi_comm_world, err)
        sbuf = reshape(buf(:, :, 3), shape(sbuf))
        call mpi_reduce(sbuf, rbuf3, nxs * nys, MPI_REAL, MPI_SUM, hdr%ionode, mpi_comm_world, err)
        sbuf = reshape(buf(:, :, 4), shape(sbuf))
        call mpi_reduce(sbuf, rbuf4, nxs * nys, MPI_REAL, MPI_SUM, hdr%ionode, mpi_comm_world, err)
        sbuf = reshape(buf(:, :, 5), shape(sbuf))
        call mpi_reduce(sbuf, rbuf5, nxs * nys, MPI_REAL, MPI_SUM, hdr%ionode, mpi_comm_world, err)
        sbuf = reshape(buf(:, :, 6), shape(sbuf))
        call mpi_reduce(sbuf, rbuf6, nxs * nys, MPI_REAL, MPI_SUM, hdr%ionode, mpi_comm_world, err)

        if (myid == hdr%ionode) then

            call nc_chk(nf90_put_att(hdr%io, hdr%medid(1), 'actual_range', (/minval(rbuf1), maxval(rbuf1)/)))
            call nc_chk(nf90_put_att(hdr%io, hdr%medid(2), 'actual_range', (/minval(rbuf2), maxval(rbuf2)/)))
            call nc_chk(nf90_put_att(hdr%io, hdr%medid(3), 'actual_range', (/minval(rbuf3), maxval(rbuf3)/)))
            call nc_chk(nf90_put_att(hdr%io, hdr%medid(4), 'actual_range', (/minval(rbuf4), maxval(rbuf4)/)))
            call nc_chk(nf90_put_att(hdr%io, hdr%medid(5), 'actual_range', (/minval(rbuf5), maxval(rbuf5)/)))
            call nc_chk(nf90_put_att(hdr%io, hdr%medid(6), 'actual_range', (/minval(rbuf6), maxval(rbuf6)/)))

            call nc_chk(nf90_enddef(hdr%io))

            call nc_chk(nf90_put_var(hdr%io, hdr%vid_x1, xsnp))
            call nc_chk(nf90_put_var(hdr%io, hdr%vid_x2, ysnp))
            call nc_chk(nf90_put_var(hdr%io, hdr%medid(1), reshape(rbuf1, shape(buf(:, :, 1)))))
            call nc_chk(nf90_put_var(hdr%io, hdr%medid(2), reshape(rbuf2, shape(buf(:, :, 1)))))
            call nc_chk(nf90_put_var(hdr%io, hdr%medid(3), reshape(rbuf3, shape(buf(:, :, 1)))))
            call nc_chk(nf90_put_var(hdr%io, hdr%medid(4), reshape(rbuf4, shape(buf(:, :, 1)))))
            call nc_chk(nf90_put_var(hdr%io, hdr%medid(5), reshape(rbuf5, shape(buf(:, :, 1)))))
            call nc_chk(nf90_put_var(hdr%io, hdr%medid(6), reshape(rbuf6, shape(buf(:, :, 1)))))

        end if

        deallocate (sbuf, rbuf1, rbuf2, rbuf3, rbuf4, rbuf5, rbuf6, buf)

    end subroutine newfile_xy_nc

    subroutine write_snp_header(hdr, ns1, ns2, xs1, xs2)

        !! write snapshot header in fixed format

        type(snp), intent(in) :: hdr
        integer, intent(in) :: ns1, ns2
        real(SP), intent(in) :: xs1(ns1), xs2(ns2) ! first and second axis
        real(SP) :: dum

        write (hdr%io) BINARY_TYPE
        write (hdr%io) CODE_TYPE
        write (hdr%io) HEADER_VERSION
        write (hdr%io) title
        write (hdr%io) exedate

        !! space grid size
        write (hdr%io) hdr%coordinate                ! coordinate
        write (hdr%io) hdr%snaptype                  ! data type
        write (hdr%io) ns1
        write (hdr%io) ns2                             ! data size
        write (hdr%io) xs1(1)
        write (hdr%io) xs2(1)
        write (hdr%io) xs1(2) - xs1(1)
        write (hdr%io) xs2(2) - xs2(1)

        !! time
        write (hdr%io) dt * ntdec_s                   ! dt

        ! absorb layer
        write (hdr%io) hdr%na1
        write (hdr%io) hdr%na2

        !! number of arrays
        write (hdr%io) hdr%nmed
        write (hdr%io) hdr%nsnp

        !! coordinate
        dum = 1.0
        write (hdr%io) clon
        write (hdr%io) clat
        write (hdr%io) phi

        write (hdr%io) dum
        write (hdr%io) dum
        write (hdr%io) dum

    end subroutine write_snp_header

    subroutine write_reduce_array2d_r(nx1, nx2, ionode, io, array)

        !! write 2d array with MPI

        integer, intent(in) :: nx1, nx2
        integer, intent(in) :: ionode
        integer, intent(in) :: io
        real(SP), intent(in) :: array(nx1, nx2)

        real(SP) :: sbuf(nx1 * nx2)
        real(SP) :: rbuf(nx1 * nx2)
        integer  :: err

        !! prepare send buffer
        sbuf = reshape(array, shape(sbuf))

        !! gather to io_node
        call mpi_reduce(sbuf, rbuf, nx1 * nx2, MPI_REAL, MPI_SUM, ionode, mpi_comm_world, err)

        !! write
        if (myid == ionode) write (io) rbuf

    end subroutine write_reduce_array2d_r


    subroutine snap__write(it)

        !! write snapshot

        integer, intent(in) :: it

        call pwatch__on("snap__write")

        if (xz_ps%sw) call wbuf_xz_ps(it)
        if (xy_ps%sw) call wbuf_xy_ps(it)
        if (fs_ps%sw) call wbuf_fs_ps(it)
        if (ob_ps%sw) call wbuf_ob_ps(it)
        if (yz_ps%sw) call wbuf_yz_ps(it)

        if (yz_v%sw) call wbuf_yz_v(it)
        if (xz_v%sw) call wbuf_xz_v(it)
        if (xy_v%sw) call wbuf_xy_v(it)
        if (fs_v%sw) call wbuf_fs_v(it)
        if (ob_v%sw) call wbuf_ob_v(it)

        if (yz_u%sw) call wbuf_yz_u(it)
        if (xz_u%sw) call wbuf_xz_u(it)
        if (xy_u%sw) call wbuf_xy_u(it)
        if (fs_u%sw) call wbuf_fs_u(it)
        if (ob_u%sw) call wbuf_ob_u(it)

        call pwatch__off("snap__write")

    end subroutine snap__write


    subroutine wbuf_nc(hdr, nvar, nx1, nx2, it0, rbuf)

        type(snp), intent(inout) :: hdr
        integer, intent(in) :: nvar
        integer, intent(in) :: nx1, nx2
        integer, intent(in) :: it0
        real, intent(in) :: rbuf(:)

        integer :: ns, ib, ie
        integer :: stt(3), cnt(3)
        integer :: vid

        ns = nx1 * nx2
        cnt = (/nx1, nx2, 1/)
        stt = (/1, 1, it0 / ntdec_s + 1/)
        do vid = 1, nvar
            ib = (vid - 1) * ns + 1
            ie = ib + ns - 1

            call nc_chk(nf90_put_var(hdr%io, hdr%varid(vid), reshape(rbuf(ib:ie), (/nx1, nx2/)), stt, cnt))
            call nc_chk(nf90_put_var(hdr%io, hdr%vid_t, it0 * dt, start=(/it0 / ntdec_s + 1/)))
            hdr%vmax(vid) = max(hdr%vmax(vid), maxval(rbuf(ib:ie)))
            hdr%vmin(vid) = min(hdr%vmin(vid), minval(rbuf(ib:ie)))
            call nc_chk(nf90_redef(hdr%io))
            call nc_chk(nf90_put_att(hdr%io, hdr%varid(vid), 'actual_range', (/hdr%vmin(vid), hdr%vmax(vid)/)))
            call nc_chk(nf90_enddef(hdr%io))
            !call nc_chk(nf90_sync(hdr%io))
        end do

    end subroutine wbuf_nc    

    subroutine wbuf_yz_ps(it)

        integer, intent(in) :: it
        integer :: i, j, k, l
        integer :: jj, kk
        real(SP) :: div, rot_x, rot_y, rot_z
        real, allocatable, save :: buf(:, :, :), sbuf(:), rbuf(:)
        integer, save :: req
        integer, save :: it0
        integer :: stat(mpi_status_size)
        integer :: err

        if (idx /= idx_yz) return

        if (.not. allocated(buf)) then
            allocate (buf(nys, nzs, 4), source=0.0)
            !$acc enter data copyin(buf)
        end if

        if ((mod(it - 1, ntdec_s) == 0) .or. (it > nt)) then

#ifdef _OPENACC
            !$acc kernels &
            !$acc present(Vx, Vy, Vz, mu, lam, buf, i0_yz, Sxy, Syz, Sxz)
            !$acc loop independent collapse(2)
#else
            !$omp parallel do private( i, j, k, kk, jj, div, rot_x, rot_y ,rot_z)
#endif
            do jj = js0, js1
                do kk = ks0, ks1

                    k = kk * kdec - kdec / 2
                    j = jj * jdec - jdec / 2
                    i = i0_yz

                    div   = real( (Vx(k  ,i  ,j  ) - Vx(k  ,i-1,j  )) * r20x &
                                + (Vy(k  ,i  ,j  ) - Vy(k  ,i  ,j-1)) * r20y &
                                + (Vz(k  ,i  ,j  ) - Vz(k-1,i  ,j  )) * r20z ) 
                    rot_x = real( (Vz(k  ,i  ,j+1) - Vz(k  ,i  ,j  )) * r20y &
                                - (Vy(k+1,i  ,j  ) - Vy(k  ,i  ,j  )) * r20z )
                    rot_y = real( (Vx(k+1,i  ,j  ) - Vx(k  ,i  ,j  )) * r20z &
                                - (Vz(k  ,i+1,j  ) - Vz(k  ,i  ,j  )) * r20x )
                    rot_z = real( (Vy(k  ,i+1,j  ) - Vy(k  ,i  ,j  )) * r20x &
                                - (Vx(k  ,i  ,j+1) - Vx(k  ,i  ,j  )) * r20y )
            
                    !! masking
                    div = div * lam(k,i,j) / abs(lam(k,i,j) + epsilon(1.0))
                    rot_x = rot_x * abs(Syz(k,i,j)) / abs(Syz(k,i,j) + epsilon(1.0))
                    rot_y = rot_y * abs(Sxz(k,i,j)) / abs(Sxz(k,i,j) + epsilon(1.0))
                    rot_z = rot_z * abs(Sxy(k,i,j)) / abs(Sxy(k,i,j) + epsilon(1.0))

                    !! dx, dy, dz have km unit. correction for 1e3 factor.
                    buf(jj,kk,1) = div   * UC * M0 * 1e-3
                    buf(jj,kk,2) = rot_x * UC * M0 * 1e-3
                    buf(jj,kk,3) = rot_y * UC * M0 * 1e-3
                    buf(jj,kk,4) = rot_z * UC * M0 * 1e-3

                end do
            end do
#ifdef _OPENACC
            !$acc end kernels
#else
            !$omp end parallel do
#endif

            if (snp_format == 'native') then
                call write_reduce_array2d_r(nys, nzs, yz_ps%ionode, yz_ps%io, buf(:, :, 1))
                call write_reduce_array2d_r(nys, nzs, yz_ps%ionode, yz_ps%io, buf(:, :, 2))
                call write_reduce_array2d_r(nys, nzs, yz_ps%ionode, yz_ps%io, buf(:, :, 3))
                call write_reduce_array2d_r(nys, nzs, yz_ps%ionode, yz_ps%io, buf(:, :, 4))
            else
                if (.not. allocated(sbuf)) then
                    allocate( sbuf(nys * nzs * 4), source=0.0 )
                    allocate( rbuf(nys * nzs * 4), source=0.0 )
                    !$acc enter data copyin(sbuf)
                else
                    call mpi_wait(req, stat, err)
                    if (myid == yz_ps%ionode) call wbuf_nc(yz_ps, 4, nys, nzs, it0, rbuf)
                end if
                if (it <= nt) then ! except for the last call
                    call pack_3d(nys, nzs, 4, buf, sbuf)

                    !$acc update self(sbuf)
                    call mpi_ireduce(sbuf, rbuf, nys * nzs * 4, mpi_real, mpi_sum, yz_ps%ionode_local, mpi_comm_yz, req, err)

                    it0 = it ! remember
                end if
            end if

        end if

    end subroutine wbuf_yz_ps

    subroutine wbuf_xz_ps(it)

        integer, intent(in) :: it
        integer :: ii, kk
        integer :: i, j, k, l
        real(SP) :: div, rot_x, rot_y, rot_z
        real, allocatable, save :: buf(:, :, :), sbuf(:), rbuf(:)
        integer, save :: req
        integer, save :: it0
        integer :: stat(mpi_status_size)
        integer :: err

        if (idy /= idy_xz) return

        if (.not. allocated(buf)) then
            allocate (buf(nxs, nzs, 4), source=0.0)
            !$acc enter data copyin(buf)
        end if

        if (mod(it - 1, ntdec_s) == 0 .or. (it > nt)) then


#ifdef _OPENACC
            !$acc kernels &
            !$acc present(Vx, Vy, Vz, Sxz, Syz, Sxy, lam, buf, j0_xz)
            !$acc loop independent collapse(2)
#else
            !$omp parallel do private( ii, kk, i, j, k, div, rot_x, rot_y, rot_z)
#endif
            do ii = is0, is1
                do kk = ks0, ks1
                    k = kk * kdec - kdec / 2
                    i = ii * idec - idec / 2
                    j = j0_xz

                    div   = real( (Vx(k  ,i  ,j  ) - Vx(k  ,i-1,j  )) * r20x &
                                + (Vy(k  ,i  ,j  ) - Vy(k  ,i  ,j-1)) * r20y &
                                + (Vz(k  ,i  ,j  ) - Vz(k-1,i  ,j  )) * r20z ) 
                    rot_x = real( (Vz(k  ,i  ,j+1) - Vz(k  ,i  ,j  )) * r20y &
                                - (Vy(k+1,i  ,j  ) - Vy(k  ,i  ,j  )) * r20z )
                    rot_y = real( (Vx(k+1,i  ,j  ) - Vx(k  ,i  ,j  )) * r20z &
                                - (Vz(k  ,i+1,j  ) - Vz(k  ,i  ,j  )) * r20x )
                    rot_z = real( (Vy(k  ,i+1,j  ) - Vy(k  ,i  ,j  )) * r20x &
                                - (Vx(k  ,i  ,j+1) - Vx(k  ,i  ,j  )) * r20y )

                    !! masking
                    div = div * lam(k,i,j) / abs(lam(k,i,j) + epsilon(1.0))
                    rot_x = rot_x * abs(Syz(k,i,j)) / (abs(Syz(k,i,j)) + epsilon(1.0))
                    rot_y = rot_y * abs(Sxz(k,i,j)) / (abs(Sxz(k,i,j)) + epsilon(1.0))
                    rot_z = rot_z * abs(Sxy(k,i,j)) / (abs(Sxy(k,i,j)) + epsilon(1.0))
                                                
                    !! dx, dy, dz have km unit. correction for 1e3 factor.
                    buf(ii, kk, 1) = div   * UC * M0 * 1e-3
                    buf(ii, kk, 2) = rot_x * UC * M0 * 1e-3
                    buf(ii, kk, 3) = rot_y * UC * M0 * 1e-3
                    buf(ii, kk, 4) = rot_z * UC * M0 * 1e-3

                end do
            end do
#ifdef _OPENACC
            !$acc end kernels
#else
            !$omp end parallel do
#endif

            if (snp_format == 'native') then
                call write_reduce_array2d_r(nxs, nzs, xz_ps%ionode, xz_ps%io, buf(:, :, 1))
                call write_reduce_array2d_r(nxs, nzs, xz_ps%ionode, xz_ps%io, buf(:, :, 2))
                call write_reduce_array2d_r(nxs, nzs, xz_ps%ionode, xz_ps%io, buf(:, :, 3))
                call write_reduce_array2d_r(nxs, nzs, xz_ps%ionode, xz_ps%io, buf(:, :, 4))
            else

                if (.not. allocated(sbuf)) then
                    allocate( sbuf(nxs * nzs * 4), source=0.0)
                    allocate( rbuf(nxs * nzs * 4), source=0.0)
                    !$acc enter data copyin(sbuf)
                else
                    call mpi_wait(req, stat, err)
                    if (myid == xz_ps%ionode)  call wbuf_nc(xz_ps, 4, nxs, nzs, it0, rbuf)
                end if
                if (it <= nt) then ! except for the last call
                    call pack_3d(nxs, nzs, 4, buf, sbuf)

                    !$acc update self(sbuf)
                    call mpi_ireduce(sbuf, rbuf, nxs * nzs * 4, mpi_real, mpi_sum, xz_ps%ionode_local, mpi_comm_xz, req, err)

                    it0 = it ! remember
                end if
            end if
        end if

    end subroutine wbuf_xz_ps

    subroutine wbuf_xy_ps(it)

        integer, intent(in) :: it
        integer :: i, j, k, l
        integer :: ii, jj
        real(SP) :: div, rot_x, rot_y, rot_z
        real, allocatable, save :: buf(:, :, :), sbuf(:), rbuf(:)
        integer, save :: req
        integer, save :: it0
        integer :: stat(mpi_status_size)
        integer :: err

        if (.not. allocated(buf)) then
            allocate (buf(nxs, nys, 4), source=0.0)
            !$acc enter data copyin(buf)
        end if

        if (mod(it - 1, ntdec_s) == 0 .or. (it > nt)) then

#ifdef _OPENACC
            !$acc kernels &
            !$acc present(Vx, Vy, Vz, Sxy, Syz, Sxz, lam, buf, k0_xy)
            !$acc loop independent collapse(2)
#else
            !$omp parallel do private( ii, jj, i, j, k, div, rot_x, rot_y, rot_z)
#endif
            do jj = js0, js1
                do ii = is0, is1
                    i = ii * idec - idec / 2
                    j = jj * jdec - jdec / 2
                    k = k0_xy

                    div   = real( (Vx(k  ,i  ,j  ) - Vx(k  ,i-1,j  )) * r20x &
                                + (Vy(k  ,i  ,j  ) - Vy(k  ,i  ,j-1)) * r20y &
                                + (Vz(k  ,i  ,j  ) - Vz(k-1,i  ,j  )) * r20z ) 
                    rot_x = real( (Vz(k  ,i  ,j+1) - Vz(k  ,i  ,j  )) * r20y &
                                - (Vy(k+1,i  ,j  ) - Vy(k  ,i  ,j  )) * r20z )
                    rot_y = real( (Vx(k+1,i  ,j  ) - Vx(k  ,i  ,j  )) * r20z &
                                - (Vz(k  ,i+1,j  ) - Vz(k  ,i  ,j  )) * r20x )
                    rot_z = real( (Vy(k  ,i+1,j  ) - Vy(k  ,i  ,j  )) * r20x &
                                - (Vx(k  ,i  ,j+1) - Vx(k  ,i  ,j  )) * r20y )

                    !! masking
                    div = div * lam(k,i,j) / abs(lam(k,i,j) + epsilon(1.0))
                    rot_x = rot_x * abs(Syz(k,i,j)) / (abs(Syz(k,i,j)) + epsilon(1.0))
                    rot_y = rot_y * abs(Sxz(k,i,j)) / (abs(Sxz(k,i,j)) + epsilon(1.0))
                    rot_z = rot_z * abs(Sxy(k,i,j)) / (abs(Sxy(k,i,j)) + epsilon(1.0))

                    !! dx, dy, dz have km unit. correction for 1e3 factor.
                    buf(ii, jj, 1) = div   * UC * M0 * 1e-3
                    buf(ii, jj, 2) = rot_x * UC * M0 * 1e-3
                    buf(ii, jj, 3) = rot_y * UC * M0 * 1e-3
                    buf(ii, jj, 4) = rot_z * UC * M0 * 1e-3

                end do
            end do
#ifdef _OPENACC
            !$acc end kernels
#else
            !$omp end parallel do
#endif

            if (snp_format == 'native') then
                call write_reduce_array2d_r(nxs, nys, xy_ps%ionode, xy_ps%io, buf(:, :, 1))
                call write_reduce_array2d_r(nxs, nys, xy_ps%ionode, xy_ps%io, buf(:, :, 2))
                call write_reduce_array2d_r(nxs, nys, xy_ps%ionode, xy_ps%io, buf(:, :, 3))
                call write_reduce_array2d_r(nxs, nys, xy_ps%ionode, xy_ps%io, buf(:, :, 4))
            else

                if (.not. allocated(sbuf)) then
                    allocate( sbuf(nxs * nys * 4), source=0.0)
                    allocate( rbuf(nxs * nys * 4), source=0.0)
                    !$acc enter data copyin(sbuf)
                else
                    call mpi_wait(req, stat, err)
                    if (myid == xy_ps%ionode) call wbuf_nc(xy_ps, 4, nxs, nys, it0, rbuf)
                end if
                if (it <= nt) then ! except for the last call
                    call pack_3d(nxs, nys, 4, buf, sbuf)

                    !$acc update self(sbuf)
                    call mpi_ireduce(sbuf, rbuf, nxs * nys * 4, mpi_real, mpi_sum, xy_ps%ionode, mpi_comm_world, req, err)

                    it0 = it ! remember
                end if

            end if

        end if

    end subroutine wbuf_xy_ps

    subroutine wbuf_fs_ps(it)

        integer, intent(in) :: it
        integer :: i, j, k
        integer :: ii, jj
        real(SP) :: div, rot_x, rot_y, rot_z
        real, allocatable, save :: buf(:, :, :), sbuf(:), rbuf(:)
        integer, save :: req
        integer, save :: it0
        integer :: stat(mpi_status_size)
        integer :: err

        if (.not. allocated(buf)) then
            allocate (buf(nxs, nys, 4), source=0.0)
            !$acc enter data copyin(buf)
        end if

        if (mod(it - 1, ntdec_s) == 0 .or. (it > nt)) then

#ifdef _OPENACC
            !$acc kernels &
            !$acc present(Vx, Vy, Vz, Syz, Sxz, Sxy, lam, buf, kfs)
            !$acc loop independent collapse(2)
#else
            !$omp parallel do private( ii, jj, i, j, k, div, rot_x, rot_y, rot_z)
#endif
            do jj = js0, js1
                do ii = is0, is1
                    j = jj * jdec - jdec / 2
                    i = ii * idec - idec / 2
                    k = kfs(i, j) + 1        !! to calculate derivatives in depth, to assure amplitude exist at detpth

                    div   = real( (Vx(k  ,i  ,j  ) - Vx(k  ,i-1,j  )) * r20x &
                                + (Vy(k  ,i  ,j  ) - Vy(k  ,i  ,j-1)) * r20y &
                                + (Vz(k  ,i  ,j  ) - Vz(k-1,i  ,j  )) * r20z ) 
                    rot_x = real( (Vz(k  ,i  ,j+1) - Vz(k  ,i  ,j  )) * r20y &
                                - (Vy(k+1,i  ,j  ) - Vy(k  ,i  ,j  )) * r20z )
                    rot_y = real( (Vx(k+1,i  ,j  ) - Vx(k  ,i  ,j  )) * r20z &
                                - (Vz(k  ,i+1,j  ) - Vz(k  ,i  ,j  )) * r20x )
                    rot_z = real( (Vy(k  ,i+1,j  ) - Vy(k  ,i  ,j  )) * r20x &
                                - (Vx(k  ,i  ,j+1) - Vx(k  ,i  ,j  )) * r20y )
            
                    !! masking
                    div = div * lam(k,i,j) / abs(lam(k,i,j) + epsilon(1.0))
                    rot_x = rot_x * abs(Syz(k,i,j)) / (abs(Syz(k,i,j)) + epsilon(1.0))
                    rot_y = rot_y * abs(Sxz(k,i,j)) / (abs(Sxz(k,i,j)) + epsilon(1.0))
                    rot_z = rot_z * abs(Sxy(k,i,j)) / (abs(Sxy(k,i,j)) + epsilon(1.0))
                                        
                    !! dx, dy, dz have km unit. correction for 1e3 factor.
                    buf(ii, jj, 1) = div   * UC * M0 * 1e-3
                    buf(ii, jj, 2) = rot_x * UC * M0 * 1e-3
                    buf(ii, jj, 3) = rot_y * UC * M0 * 1e-3
                    buf(ii, jj, 4) = rot_z * UC * M0 * 1e-3
                end do
            end do
#ifdef _OPENACC
            !$acc end kernels
#else
            !$omp end parallel do
#endif

            if (snp_format == 'native') then
                call write_reduce_array2d_r(nxs, nys, fs_ps%ionode, fs_ps%io, buf(:, :, 1))
                call write_reduce_array2d_r(nxs, nys, fs_ps%ionode, fs_ps%io, buf(:, :, 2))
                call write_reduce_array2d_r(nxs, nys, fs_ps%ionode, fs_ps%io, buf(:, :, 3))
                call write_reduce_array2d_r(nxs, nys, fs_ps%ionode, fs_ps%io, buf(:, :, 4))
            else

                if (.not. allocated(sbuf)) then
                    allocate( sbuf(nxs * nys * 4), source=0.0)
                    allocate( rbuf(nxs * nys * 4), source=0.0)
                    !$acc enter data copyin(sbuf)
                else
                    call mpi_wait(req, stat, err)
                    if (myid == fs_ps%ionode) call wbuf_nc(fs_ps, 4, nxs, nys, it0, rbuf)
                end if
                if (it <= nt) then ! except for the last call
                    call pack_3d(nxs, nys, 4, buf, sbuf)

                    !$acc update self(sbuf)
                    call mpi_ireduce(sbuf, rbuf, nxs * nys * 4, mpi_real, mpi_sum, fs_ps%ionode, mpi_comm_world, req, err)

                    it0 = it ! remember
                end if

            end if

        end if

    end subroutine wbuf_fs_ps

    subroutine wbuf_ob_ps(it)

        integer, intent(in) :: it
        integer :: i, j, k
        integer :: ii, jj
        real(SP) :: div, rot_x, rot_y, rot_z
        real, allocatable, save :: buf(:, :, :), sbuf(:), rbuf(:)
        integer, save :: req
        integer, save :: it0
        integer :: stat(mpi_status_size)
        integer :: err

        if (.not. allocated(buf)) then
            allocate (buf(nxs, nys, 4), source=0.0)
            !$acc enter data copyin(buf)
        end if

        if (mod(it - 1, ntdec_s) == 0 .or. (it > nt)) then

#ifdef _OPENACC
            !$acc kernels &
            !$acc present(Vx, Vy, Vz, Syz, Sxz, Sxy, lam, buf, kob)
            !$acc loop independent collapse(2)
#else
            !$omp parallel do private( ii, jj, i, j, k, div, rot_x, rot_y, rot_z )
#endif            
            do jj = js0, js1
                do ii = is0, is1
                    j = jj * jdec - jdec / 2
                    i = ii * idec - idec / 2
                    k = kob(i, j) + 1        !! to calculate derivatives in depth, to assure amplitude exist at detpth

                    div   = real( (Vx(k  ,i  ,j  ) - Vx(k  ,i-1,j  )) * r20x &
                                + (Vy(k  ,i  ,j  ) - Vy(k  ,i  ,j-1)) * r20y &
                                + (Vz(k  ,i  ,j  ) - Vz(k-1,i  ,j  )) * r20z ) 
                    rot_x = real( (Vz(k  ,i  ,j+1) - Vz(k  ,i  ,j  )) * r20y &
                                - (Vy(k+1,i  ,j  ) - Vy(k  ,i  ,j  )) * r20z )
                    rot_y = real( (Vx(k+1,i  ,j  ) - Vx(k  ,i  ,j  )) * r20z &
                                - (Vz(k  ,i+1,j  ) - Vz(k  ,i  ,j  )) * r20x )
                    rot_z = real( (Vy(k  ,i+1,j  ) - Vy(k  ,i  ,j  )) * r20x &
                                - (Vx(k  ,i  ,j+1) - Vx(k  ,i  ,j  )) * r20y )
            
                    !! masking
                    div = div * lam(k,i,j) / abs(lam(k,i,j) + epsilon(1.0))
                    rot_x = rot_x * abs(Syz(k,i,j)) / (abs(Syz(k,i,j)) + epsilon(1.0))
                    rot_y = rot_y * abs(Sxz(k,i,j)) / (abs(Sxz(k,i,j)) + epsilon(1.0))
                    rot_z = rot_z * abs(Sxy(k,i,j)) / (abs(Sxy(k,i,j)) + epsilon(1.0))

                    !! dx, dy, dz have km unit. correction for 1e3 factor.
                    buf(ii,jj,1) = div   * UC * M0 * 1e-3
                    buf(ii,jj,2) = rot_x * UC * M0 * 1e-3
                    buf(ii,jj,3) = rot_y * UC * M0 * 1e-3
                    buf(ii,jj,4) = rot_z * UC * M0 * 1e-3

                end do
            end do
#ifdef _OPENACC
            !$acc end kernels
#else
            !$omp end parallel do
#endif

            if (snp_format == 'native') then
                call write_reduce_array2d_r(nxs, nys, ob_ps%ionode, ob_ps%io, buf(:, :, 1))
                call write_reduce_array2d_r(nxs, nys, ob_ps%ionode, ob_ps%io, buf(:, :, 2))
                call write_reduce_array2d_r(nxs, nys, ob_ps%ionode, ob_ps%io, buf(:, :, 3))
                call write_reduce_array2d_r(nxs, nys, ob_ps%ionode, ob_ps%io, buf(:, :, 4))
            else

                if (.not. allocated(sbuf)) then
                    allocate( sbuf(nxs * nys * 4), source=0.0)
                    allocate( rbuf(nxs * nys * 4), source=0.0)
                    !$acc enter data copyin(sbuf)
                else
                    call mpi_wait(req, stat, err)
                    if (myid == ob_ps%ionode) call wbuf_nc(ob_ps, 4, nxs, nys, it0, rbuf)
                end if
                if (it <= nt) then ! except for the last call
                    call pack_3d(nxs, nys, 4, buf, sbuf)

                    !$acc update self(sbuf)
                    call mpi_ireduce(sbuf, rbuf, nxs * nys * 4, mpi_real, mpi_sum, ob_ps%ionode, mpi_comm_world, req, err)

                    it0 = it ! remember
                end if
            end if

        end if

    end subroutine wbuf_ob_ps

    subroutine wbuf_yz_v(it)
        integer, intent(in) :: it
        integer :: i, j, k
        integer :: jj, kk
        real, allocatable, save :: buf(:, :, :), sbuf(:), rbuf(:)
        integer, save :: req
        integer, save :: it0
        integer :: stat(mpi_status_size)
        integer :: err

        if (idx /= idx_yz) return

        if (.not. allocated(buf)) then
            allocate (buf(nys, nzs, 3), source=0.0)
            !$acc enter data copyin(buf)
        end if

        if (mod(it - 1, ntdec_s) == 0 .or. (it > nt)) then


            buf = 0.0

#ifdef _OPENACC
            !$acc kernels &
            !$acc present(Vx, Vy, Vz, buf)
            !$acc loop independent collapse(2)
#else
            !$omp parallel do private( jj, kk, k, j, i )
#endif
            do jj = js0, js1
                do kk = ks0, ks1
                    i = i0_yz
                    j = jj * jdec - jdec / 2
                    k = kk * kdec - kdec / 2

                    buf(jj, kk, 1) = Vx(k, i, j) * UC * M0
                    buf(jj, kk, 2) = Vy(k, i, j) * UC * M0
                    buf(jj, kk, 3) = Vz(k, i, j) * UC * M0

                end do
            end do
#ifdef _OPENACC
            !$acc end kernels
#else
            !$omp end parallel do
#endif

            if (snp_format == 'native') then
                call write_reduce_array2d_r(nys, nzs, yz_v%ionode, yz_v%io, buf(:, :, 1))
                call write_reduce_array2d_r(nys, nzs, yz_v%ionode, yz_v%io, buf(:, :, 2))
                call write_reduce_array2d_r(nys, nzs, yz_v%ionode, yz_v%io, buf(:, :, 3))
            else

                if (.not. allocated(sbuf)) then
                    allocate( sbuf(nys * nzs * 3), source=0.0)
                    allocate( rbuf(nys * nzs * 3), source=0.0)
                    !$acc enter data copyin(sbuf)
                else
                    call mpi_wait(req, stat, err)
                    if (myid == yz_v%ionode) call wbuf_nc(yz_v, 3, nys, nzs, it0, rbuf)
                end if
                if (it <= nt) then ! except for the last call
                    call pack_3d(nys, nzs, 3, buf, sbuf)

                    !$acc update self(sbuf)
                    call mpi_ireduce(sbuf, rbuf, nys * nzs * 3, mpi_real, mpi_sum, yz_v%ionode_local, mpi_comm_yz, req, err)

                    it0 = it ! remember
                end if
            end if

        end if

    end subroutine wbuf_yz_v

    subroutine wbuf_xz_v(it)

        integer, intent(in) :: it
        integer :: i, j, k
        integer :: ii, kk
        real, allocatable, save :: buf(:, :, :), sbuf(:), rbuf(:)
        integer, save :: req
        integer, save :: it0
        integer :: stat(mpi_status_size)
        integer :: err

        if (idy /= idy_xz) return

        if (.not. allocated(buf)) then
            allocate (buf(nxs, nzs, 3), source=0.0)
            !$acc enter data copyin(buf)
        end if

        if (mod(it - 1, ntdec_s) == 0 .or. (it > nt)) then


#ifdef _OPENACC
            !$acc kernels &
            !$acc present(Vx, Vy, Vz, buf)
            !$acc loop independent collapse(2)
#else
            !$omp parallel do private( ii, kk, i, k, j )
#endif
            do ii = is0, is1
                do kk = ks0, ks1
                    i = ii * idec - idec / 2
                    j = j0_xz
                    k = kk * kdec - kdec / 2

                    buf(ii, kk, 1) = Vx(k, i, j) * UC * M0
                    buf(ii, kk, 2) = Vy(k, i, j) * UC * M0
                    buf(ii, kk, 3) = Vz(k, i, j) * UC * M0

                end do
            end do
#ifdef _OPENACC
            !$acc end kernels
#else
            !$omp end parallel do
#endif

            if (snp_format == 'native') then
                call write_reduce_array2d_r(nxs, nzs, xz_v%ionode, xz_v%io, buf(:, :, 1))
                call write_reduce_array2d_r(nxs, nzs, xz_v%ionode, xz_v%io, buf(:, :, 2))
                call write_reduce_array2d_r(nxs, nzs, xz_v%ionode, xz_v%io, buf(:, :, 3))
            else
                if (.not. allocated(sbuf)) then
                    allocate( sbuf(nxs * nzs * 3), source=0.0)
                    allocate( rbuf(nxs * nzs * 3), source=0.0)
                    !$acc enter data copyin(sbuf)
                else
                    call mpi_wait(req, stat, err)
                    if (myid == xz_v%ionode) call wbuf_nc(xz_v, 3, nxs, nzs, it0, rbuf)
                end if
                if (it <= nt) then ! except for the last call
                    call pack_3d(nxs, nzs, 3, buf, sbuf)

                    !$acc update self(sbuf)
                    call mpi_ireduce(sbuf, rbuf, nxs * nzs * 3, mpi_real, mpi_sum, xz_v%ionode_local, mpi_comm_xz, req, err)
                    
                    it0 = it ! remember
                end if
            end if
        end if
    end subroutine wbuf_xz_v

    subroutine wbuf_xy_v(it)
        integer, intent(in) :: it
        integer :: i, j, k
        integer :: ii, jj
        real, allocatable, save :: buf(:, :, :), sbuf(:), rbuf(:)
        integer, save :: req
        integer, save :: it0
        integer :: stat(mpi_status_size)
        integer :: err

        if (.not. allocated(buf)) then
            allocate (buf(nxs, nys, 3), source=0.0)
            !$acc enter data copyin(buf)
        end if

        if (mod(it - 1, ntdec_s) == 0 .or. (it > nt)) then

#ifdef _OPENACC
            !$acc kernels &
            !$acc present(Vx, Vy, Vz, buf)
            !$acc loop independent collapse(2)
#else
            !$omp parallel do private( ii, jj, i, j, k )
#endif
            do jj = js0, js1
                do ii = is0, is1
                    i = ii * idec - idec / 2
                    j = jj * jdec - jdec / 2
                    k = k0_xy

                    buf(ii, jj, 1) = Vx(k, i, j) * UC * M0
                    buf(ii, jj, 2) = Vy(k, i, j) * UC * M0
                    buf(ii, jj, 3) = Vz(k, i, j) * UC * M0

                end do
            end do
#ifdef _OPENACC
            !$acc end kernels
#else
            !$omp end parallel do
#endif

            if (snp_format == 'native') then
                call write_reduce_array2d_r(nxs, nys, xy_v%ionode, xy_v%io, buf(:, :, 1))
                call write_reduce_array2d_r(nxs, nys, xy_v%ionode, xy_v%io, buf(:, :, 2))
                call write_reduce_array2d_r(nxs, nys, xy_v%ionode, xy_v%io, buf(:, :, 3))
            else
                if (.not. allocated(sbuf)) then
                    allocate( sbuf(nxs * nys * 3), source=0.0)
                    allocate( rbuf(nxs * nys * 3), source=0.0)
                    !$acc enter data copyin(sbuf)
                else
                    call mpi_wait(req, stat, err)
                    if (myid == xy_v%ionode) call wbuf_nc(xy_v, 3, nxs, nys, it0, rbuf)
                end if
                if (it <= nt) then ! except for the last call
                    call pack_3d(nxs, nys, 3, buf, sbuf)

                    !$acc update self(sbuf)
                    call mpi_ireduce(sbuf, rbuf, nxs * nys * 3, mpi_real, mpi_sum, xy_v%ionode, mpi_comm_world, req, err)

                    it0 = it ! remember
                end if
            end if
        end if

    end subroutine wbuf_xy_v

    subroutine wbuf_fs_v(it)

        integer, intent(in) :: it
        integer :: i, j, k
        integer :: ii, jj
        real, allocatable, save :: buf(:, :, :), sbuf(:), rbuf(:)
        integer, save :: req
        integer, save :: it0
        integer :: stat(mpi_status_size)
        integer :: err

        if (.not. allocated(buf)) then
            allocate (buf(nxs, nys, 3), source=0.0)
            !$acc enter data copyin(buf)
        end if

#ifdef _OPENACC
        !$acc kernels &
        !$acc present(Vx, Vy, Vz, buf, kfs)
        !$acc loop independent collapse(2)
#else
        !$omp parallel do private( ii, jj, i, j, k )
#endif
        do jj = js0, js1
            do ii = is0, is1
                j = jj * jdec - jdec / 2
                i = ii * idec - idec / 2
                k = kfs(i, j) + 1

                buf(ii, jj, 1) = Vx(k, i, j) * UC * M0
                buf(ii, jj, 2) = Vy(k, i, j) * UC * M0
                buf(ii, jj, 3) = Vz(k, i, j) * UC * M0

            end do
        end do
#ifdef _OPENACC
        !$acc end kernels
#else
        !$omp end parallel do
#endif

#ifdef _OPENACC
        !$acc kernels present(buf, max_fs_v)
        !$acc loop independent collapse(2)
#else
        !$omp parallel do private(ii,jj)
#endif
        do jj = js0, js1
            do ii = is0, is1
                max_fs_v(ii, jj, 1) = max(max_fs_v(ii, jj, 1), abs(buf(ii, jj, 3)))
                max_fs_v(ii, jj, 2) = max(max_fs_v(ii, jj, 2), sqrt(buf(ii, jj, 1)**2 + buf(ii, jj, 2)**2))
                max_fs_v(ii, jj, 3) = sqrt(max_fs_v(ii, jj, 1)**2 + max_fs_v(ii, jj, 2)**2)
            end do
        end do
#ifdef _OPENACC
        !$acc end kernels
#else
        !$omp end parallel do
#endif

        if (mod(it - 1, ntdec_s) == 0 .or. (it > nt)) then

            if (snp_format == 'native') then
                call write_reduce_array2d_r(nxs, nys, fs_v%ionode, fs_v%io, buf(:, :, 1))
                call write_reduce_array2d_r(nxs, nys, fs_v%ionode, fs_v%io, buf(:, :, 2))
                call write_reduce_array2d_r(nxs, nys, fs_v%ionode, fs_v%io, buf(:, :, 3))
            else
                if (.not. allocated(sbuf)) then
                    allocate( sbuf(nxs * nys * 3), source=0.0)
                    allocate( rbuf(nxs * nys * 3), source=0.0)
                    !$acc enter data copyin(sbuf)
                else
                    call mpi_wait(req, stat, err)
                    if (myid == fs_v%ionode) call wbuf_nc(fs_v, 3, nxs, nys, it0, rbuf)
                end if
                if (it <= nt) then ! except for the last call
                    call pack_3d(nxs, nys, 3, buf, sbuf)

                    !$acc update self(sbuf)
                    call mpi_ireduce(sbuf, rbuf, nxs * nys * 3, mpi_real, mpi_sum, fs_v%ionode, mpi_comm_world, req, err)

                    it0 = it ! remember
                end if
            end if
        end if

    end subroutine wbuf_fs_v

    subroutine wbuf_ob_v(it)

        integer, intent(in) :: it
        integer :: i, j, k
        integer :: ii, jj
        real, allocatable, save :: buf(:, :, :), sbuf(:), rbuf(:)
        integer, save :: req
        integer, save :: it0
        integer :: stat(mpi_status_size)
        integer :: err

        if (.not. allocated(buf)) then
            allocate (buf(nxs, nys, 3), source=0.0)
            !$acc enter data copyin(buf)
        end if

#ifdef _OPENACC
        !$acc kernels &
        !$acc present(Vx, Vy, Vz, buf, kob)
        !$acc loop independent collapse(2)
#else
        !$omp parallel do private( ii, jj, i, j, k )
#endif
        do jj = js0, js1
            do ii = is0, is1
                j = jj * jdec - jdec / 2
                i = ii * idec - idec / 2
                k = kob(i, j) + 1

                buf(ii, jj, 1) = Vx(k, i, j) * UC * M0
                buf(ii, jj, 2) = Vy(k, i, j) * UC * M0
                buf(ii, jj, 3) = Vz(k, i, j) * UC * M0

            end do
        end do
#ifdef _OPENACC
        !$acc end kernels
#else
        !$omp end parallel do
#endif

#ifdef _OPENACC
        !$acc kernels present(buf, max_ob_v)
        !$acc loop independent collapse(2)
#else
        !$omp parallel do private(ii,jj)
#endif
        do jj = js0, js1
            do ii = is0, is1
                max_ob_v(ii, jj, 1) = max(max_ob_v(ii, jj, 1), abs(buf(ii, jj, 3)))
                max_ob_v(ii, jj, 2) = max(max_ob_v(ii, jj, 2), sqrt(buf(ii, jj, 1)**2 + buf(ii, jj, 2)**2))
                max_ob_v(ii, jj, 3) = sqrt(max_ob_v(ii, jj, 1)**2 + max_ob_v(ii, jj, 2)**2)
            end do
        end do
#ifdef _OPENACC
        !$acc end kernels
#else
        !$omp end parallel do
#endif

        if (mod(it - 1, ntdec_s) == 0 .or. (it > nt)) then

            if (snp_format == 'native') then
                call write_reduce_array2d_r(nxs, nys, ob_v%ionode, ob_v%io, buf(:, :, 1))
                call write_reduce_array2d_r(nxs, nys, ob_v%ionode, ob_v%io, buf(:, :, 2))
                call write_reduce_array2d_r(nxs, nys, ob_v%ionode, ob_v%io, buf(:, :, 3))
            else
                if (.not. allocated(sbuf)) then
                    allocate( sbuf(nxs * nys * 3), source=0.0)
                    allocate( rbuf(nxs * nys * 3), source=0.0)
                    !$acc enter data copyin(sbuf)
                else
                    call mpi_wait(req, stat, err)
                    if (myid == ob_v%ionode) call wbuf_nc(ob_v, 3, nxs, nys, it0, rbuf)
                end if
                if (it <= nt) then ! except for the last call
                    call pack_3d(nxs, nys, 3, buf, sbuf)

                    !$acc update self(sbuf)
                    call mpi_ireduce(sbuf, rbuf, nxs * nys * 3, mpi_real, mpi_sum, ob_v%ionode, mpi_comm_world, req, err)

                    it0 = it ! remember
                end if
            end if
        end if

    end subroutine wbuf_ob_v

    subroutine wbuf_yz_u(it)
        integer, intent(in) :: it
        integer :: i, j, k
        integer :: jj, kk
        real, allocatable, save :: sbuf(:), rbuf(:)
        integer, save :: req
        integer, save :: it0
        integer :: stat(mpi_status_size)
        integer :: err

        if (idx /= idx_yz) return


#ifdef _OPENACC
        !$acc kernels &
        !$acc present(Vx, Vy, Vz, buf_yz_u)
        !$acc loop independent collapse(2)
#else            
        !$omp parallel do private( jj, kk, k, i, j )
#endif
        do jj = js0, js1
            do kk = ks0, ks1
                i = i0_yz
                k = kk * kdec - kdec / 2
                j = jj * jdec - jdec / 2

                buf_yz_u(jj, kk, 1) = buf_yz_u(jj, kk, 1) + Vx(k, i, j) * UC * M0 * dt
                buf_yz_u(jj, kk, 2) = buf_yz_u(jj, kk, 2) + Vy(k, i, j) * UC * M0 * dt
                buf_yz_u(jj, kk, 3) = buf_yz_u(jj, kk, 3) + Vz(k, i, j) * UC * M0 * dt

            end do
        end do
#ifdef _OPENACC
        !$acc end kernels
#else
        !$omp end parallel do
#endif

        if (mod(it - 1, ntdec_s) == 0 .or. (it > nt)) then
            if (snp_format == 'native') then
                call write_reduce_array2d_r(nys, nzs, yz_u%ionode, yz_u%io, buf_yz_u(:, :, 1))
                call write_reduce_array2d_r(nys, nzs, yz_u%ionode, yz_u%io, buf_yz_u(:, :, 2))
                call write_reduce_array2d_r(nys, nzs, yz_u%ionode, yz_u%io, buf_yz_u(:, :, 3))
            else
                if (.not. allocated(sbuf)) then
                    allocate(sbuf(nys * nzs * 3), source=0.0)
                    allocate(rbuf(nys * nzs * 3), source=0.0)
                    !$acc enter data copyin(sbuf)
                else
                    call mpi_wait(req, stat, err)
                    if (myid == yz_u%ionode) call wbuf_nc(yz_u, 3, nys, nzs, it0, rbuf)
                end if
                if (it <= nt) then ! except for the last call
                    call pack_3d(nys, nzs, 3, buf_yz_u, sbuf)

                    !$acc update self(sbuf)
                    call mpi_ireduce(sbuf, rbuf, nys * nzs * 3, mpi_real, mpi_sum, yz_u%ionode_local, mpi_comm_yz, req, err)

                    it0 = it ! remember
                end if
            end if
        end if

    end subroutine wbuf_yz_u

    subroutine wbuf_xz_u(it)

        integer, intent(in) :: it
        integer :: i, j, k
        integer :: ii, kk
        real, allocatable, save :: sbuf(:), rbuf(:)
        integer, save :: req
        integer, save :: it0
        integer :: stat(mpi_status_size)
        integer :: err

        if (idy /= idy_xz) return

#ifdef _OPENACC
        !$acc kernels &
        !$acc present(Vx, Vy, Vz, buf_xz_u)
        !$acc loop independent collapse(2)
#else            
        !$omp parallel do private( ii, kk, k, i, j )
#endif
        do ii = is0, is1
            do kk = ks0, ks1
                k = kk * kdec - kdec / 2
                i = ii * idec - idec / 2
                j = j0_xz

                buf_xz_u(ii, kk, 1) = buf_xz_u(ii, kk, 1) + Vx(k, i, j) * UC * M0 * dt
                buf_xz_u(ii, kk, 2) = buf_xz_u(ii, kk, 2) + Vy(k, i, j) * UC * M0 * dt
                buf_xz_u(ii, kk, 3) = buf_xz_u(ii, kk, 3) + Vz(k, i, j) * UC * M0 * dt

            end do
        end do
#ifdef _OPENACC
        !$acc end kernels
#else
        !$omp end parallel do
#endif

        if (mod(it - 1, ntdec_s) == 0 .or. (it > nt)) then

            if (snp_format == 'native') then
                call write_reduce_array2d_r(nxs, nzs, xz_u%ionode, xz_u%io, buf_xz_u(:, :, 1))
                call write_reduce_array2d_r(nxs, nzs, xz_u%ionode, xz_u%io, buf_xz_u(:, :, 2))
                call write_reduce_array2d_r(nxs, nzs, xz_u%ionode, xz_u%io, buf_xz_u(:, :, 3))
            else
                if (.not. allocated(sbuf)) then
                    allocate( sbuf(nxs * nzs * 3), source=0.0)
                    allocate( rbuf(nxs * nzs * 3), source=0.0)
                    !$acc enter data copyin(sbuf)
                else
                    call mpi_wait(req, stat, err)
                    if (myid == xz_u%ionode) call wbuf_nc(xz_u, 3, nxs, nzs, it0, rbuf)
                end if
                if (it <= nt) then ! except for the last call
                    call pack_3d(nxs, nzs, 3, buf_xz_u, sbuf)

                    !$acc update self(sbuf)
                    call mpi_ireduce(sbuf, rbuf, nxs * nzs * 3, mpi_real, mpi_sum, xz_u%ionode_local, mpi_comm_xz, req, err)

                    it0 = it ! remember
                end if
            end if

        end if

    end subroutine wbuf_xz_u

    subroutine wbuf_xy_u(it)
        integer, intent(in) :: it
        integer :: i, j, k
        integer :: ii, jj
        real, allocatable, save :: sbuf(:), rbuf(:)
        integer, save :: req
        integer, save :: it0
        integer :: stat(mpi_status_size)
        integer :: err

#ifdef _OPENACC
        !$acc kernels &
        !$acc present(Vx, Vy, Vz, buf_xy_u)
        !$acc loop independent collapse(2)
#else            
        !$omp parallel do private( jj, ii, k, i, j )
#endif
        do jj = js0, js1
            do ii = is0, is1
                j = jj * jdec - jdec / 2
                i = ii * idec - idec / 2
                k = k0_xy

                buf_xy_u(ii, jj, 1) = buf_xy_u(ii, jj, 1) + Vx(k, i, j) * UC * M0 * dt
                buf_xy_u(ii, jj, 2) = buf_xy_u(ii, jj, 2) + Vy(k, i, j) * UC * M0 * dt
                buf_xy_u(ii, jj, 3) = buf_xy_u(ii, jj, 3) + Vz(k, i, j) * UC * M0 * dt

            end do
        end do
#ifdef _OPENACC
        !$acc end kernels
#else
        !$omp end parallel do
#endif

        if (mod(it - 1, ntdec_s) == 0 .or. (it > nt)) then

            if (snp_format == 'native') then
                call write_reduce_array2d_r(nxs, nys, xy_u%ionode, xy_u%io, buf_xy_u(:, :, 1))
                call write_reduce_array2d_r(nxs, nys, xy_u%ionode, xy_u%io, buf_xy_u(:, :, 2))
                call write_reduce_array2d_r(nxs, nys, xy_u%ionode, xy_u%io, buf_xy_u(:, :, 3))
            else
                if (.not. allocated(sbuf)) then
                    allocate( sbuf(nxs * nys * 3), source=0.0)
                    allocate( rbuf(nxs * nys * 3), source=0.0)
                    !$acc enter data copyin(sbuf)
                else
                    call mpi_wait(req, stat, err)
                    if (myid == xy_u%ionode) call wbuf_nc(xy_u, 3, nxs, nys, it0, rbuf)
                end if
                if (it <= nt) then ! except for the last call
                    call pack_3d(nxs, nys, 3, buf_xy_u, sbuf)

                    !$acc update self(sbuf)
                    call mpi_ireduce(sbuf, rbuf, nxs * nys * 3, mpi_real, mpi_sum, xy_u%ionode, mpi_comm_world, req, err)

                    it0 = it ! remember
                end if
            end if

        end if

    end subroutine wbuf_xy_u

    subroutine wbuf_fs_u(it)

        integer, intent(in) :: it
        integer :: i, j, k
        integer :: ii, jj
        real, allocatable, save :: sbuf(:), rbuf(:)
        integer, save :: req
        integer, save :: it0
        integer :: stat(mpi_status_size)
        integer :: err

#ifdef _OPENACC
        !$acc kernels &
        !$acc present(Vx, Vy, Vz, buf_fs_u, kfs)
        !$acc loop independent collapse(2)
#else            
        !$omp parallel do private( jj, ii, k, i, j )
#endif
        do jj = js0, js1
            do ii = is0, is1
                j = jj * jdec - jdec / 2
                i = ii * idec - idec / 2
                k = kfs(i, j) + 1

                buf_fs_u(ii, jj, 1) = buf_fs_u(ii, jj, 1) + Vx(k, i, j) * UC * M0 * dt
                buf_fs_u(ii, jj, 2) = buf_fs_u(ii, jj, 2) + Vy(k, i, j) * UC * M0 * dt
                buf_fs_u(ii, jj, 3) = buf_fs_u(ii, jj, 3) + Vz(k, i, j) * UC * M0 * dt

            end do
        end do
#ifdef _OPENACC
        !$acc end kernels
#else
        !$omp end parallel do
#endif

#ifdef _OPENACC
        !$acc kernels present(buf_fs_u, max_fs_u)
        !$acc loop independent collapse(2)
#else
        !$omp parallel do private(ii,jj)
#endif
        do jj = js0, js1
            do ii = is0, is1
                max_fs_u(ii, jj, 1) = max(max_fs_u(ii, jj, 1), abs(buf_fs_u(ii, jj, 3)))
                max_fs_u(ii, jj, 2) = max(max_fs_u(ii, jj, 2), sqrt(buf_fs_u(ii, jj, 1)**2 + buf_fs_u(ii, jj, 2)**2))
                max_fs_u(ii, jj, 3) = sqrt(max_fs_u(ii, jj, 1)**2 + max_fs_u(ii, jj, 2)**2)
            end do
        end do
#ifdef _OPENACC
        !$acc end kernels
#else
        !$omp end parallel do
#endif

        if (mod(it - 1, ntdec_s) == 0 .or. (it > nt)) then
            if (snp_format == 'native') then
                call write_reduce_array2d_r(nxs, nys, fs_u%ionode, fs_u%io, buf_fs_u(:, :, 1))
                call write_reduce_array2d_r(nxs, nys, fs_u%ionode, fs_u%io, buf_fs_u(:, :, 2))
                call write_reduce_array2d_r(nxs, nys, fs_u%ionode, fs_u%io, buf_fs_u(:, :, 3))
            else
                if (.not. allocated(sbuf)) then
                    allocate( sbuf(nxs * nys * 3), source=0.0)
                    allocate( rbuf(nxs * nys * 3), source=0.0)
                    !$acc enter data copyin(sbuf)
                else
                    call mpi_wait(req, stat, err)
                    if (myid == fs_u%ionode) call wbuf_nc(fs_u, 3, nxs, nys, it0, rbuf)
                end if
                if (it <= nt) then ! except for the last call
                    call pack_3d(nxs, nys, 3, buf_fs_u, sbuf)

                    !$acc update self(sbuf)
                    call mpi_ireduce(sbuf, rbuf, nxs * nys * 3, mpi_real, mpi_sum, fs_u%ionode, mpi_comm_world, req, err)

                    it0 = it ! remember
                end if
            end if
        end if

    end subroutine wbuf_fs_u

    subroutine wbuf_ob_u(it)

        integer, intent(in) :: it
        integer :: i, j, k
        integer :: ii, jj
        real, allocatable, save :: sbuf(:), rbuf(:)
        integer, save :: req
        integer, save :: it0
        integer :: stat(mpi_status_size)
        integer :: err

#ifdef _OPENACC
        !$acc kernels &
        !$acc present(Vx, Vy, Vz, buf_ob_u, kob)
        !$acc loop independent collapse(2)
#else            
        !$omp parallel do private( jj, ii, k, i, j )
#endif
        do jj = js0, js1
            do ii = is0, is1
                j = jj * jdec - jdec / 2
                i = ii * idec - idec / 2
                k = kob(i, j) + 1

                buf_ob_u(ii, jj, 1) = buf_ob_u(ii, jj, 1) + Vx(k, i, j) * UC * M0 * dt
                buf_ob_u(ii, jj, 2) = buf_ob_u(ii, jj, 2) + Vy(k, i, j) * UC * M0 * dt
                buf_ob_u(ii, jj, 3) = buf_ob_u(ii, jj, 3) + Vz(k, i, j) * UC * M0 * dt

            end do
        end do
#ifdef _OPENACC
        !$acc end kernels
#else
        !$omp end parallel do
#endif

#ifdef _OPENACC
        !$acc kernels present(buf_ob_u, max_ob_u)
        !$acc loop independent collapse(2)
#else
        !$omp parallel do private(ii,jj)
#endif
        do jj = js0, js1
            do ii = is0, is1
                max_ob_u(ii, jj, 1) = max(max_ob_u(ii, jj, 1), abs(buf_ob_u(ii, jj, 3)))
                max_ob_u(ii, jj, 2) = max(max_ob_u(ii, jj, 2), sqrt(buf_ob_u(ii, jj, 1)**2 + buf_ob_u(ii, jj, 2)**2))
                max_ob_u(ii, jj, 3) = sqrt(max_ob_u(ii, jj, 1)**2 + max_ob_u(ii, jj, 2)**2)
            end do
        end do
#ifdef _OPENACC
        !$acc end kernels
#else
        !$omp end parallel do
#endif

        if (mod(it - 1, ntdec_s) == 0 .or. (it > nt)) then
            if (snp_format == 'native') then
                call write_reduce_array2d_r(nxs, nys, ob_u%ionode, ob_u%io, buf_ob_u(:, :, 1))
                call write_reduce_array2d_r(nxs, nys, ob_u%ionode, ob_u%io, buf_ob_u(:, :, 2))
                call write_reduce_array2d_r(nxs, nys, ob_u%ionode, ob_u%io, buf_ob_u(:, :, 3))
            else
                if (.not. allocated(sbuf)) then
                    allocate( sbuf(nxs * nys * 3), source=0.0)
                    allocate( rbuf(nxs * nys * 3), source=0.0)
                    !$acc enter data copyin(sbuf)
                else
                    call mpi_wait(req, stat, err)
                    if (myid == ob_u%ionode) call wbuf_nc(ob_u, 3, nxs, nys, it0, rbuf)
                end if
                if (it <= nt) then ! except for the last call
                    call pack_3d(nxs, nys, 3, buf_ob_u, sbuf)

                    !$acc update self(sbuf)
                    call mpi_ireduce(sbuf, rbuf, nxs * nys * 3, mpi_real, mpi_sum, ob_u%ionode, mpi_comm_world, req, err)

                    it0 = it ! remember
                end if
            end if
        end if

    end subroutine wbuf_ob_u

    subroutine close_nc(hdr)

        type(snp), intent(in) :: hdr
        integer :: vid

        call nc_chk(nf90_redef(hdr%io))
        do vid = 1, hdr%nsnp
            call nc_chk(nf90_put_att(hdr%io, hdr%varid(vid), 'actual_range', (/hdr%vmin(vid), hdr%vmax(vid)/)))
        end do
        call nc_chk(nf90_enddef(hdr%io))
        call nc_chk(nf90_sync(hdr%io))
        call nc_chk(nf90_close(hdr%io))

    end subroutine close_nc

    subroutine snap__closefiles

        call pwatch__on('snap__closefiles')

        if (snp_format == 'native') then
            if (yz_ps%sw .and. myid == yz_ps%ionode) close (yz_ps%io)
            if (xz_ps%sw .and. myid == xz_ps%ionode) close (xz_ps%io)
            if (xy_ps%sw .and. myid == xy_ps%ionode) close (xy_ps%io)
            if (fs_ps%sw .and. myid == fs_ps%ionode) close (fs_ps%io)
            if (ob_ps%sw .and. myid == ob_ps%ionode) close (ob_ps%io)

            if (yz_v%sw .and. myid == yz_v%ionode) close (yz_v%io)
            if (xz_v%sw .and. myid == xz_v%ionode) close (xz_v%io)
            if (xy_v%sw .and. myid == xy_v%ionode) close (xy_v%io)
            if (fs_v%sw .and. myid == fs_v%ionode) close (fs_v%io)
            if (ob_v%sw .and. myid == ob_v%ionode) close (ob_v%io)

            if (yz_u%sw .and. myid == yz_u%ionode) close (yz_u%io)
            if (xz_u%sw .and. myid == xz_u%ionode) close (xz_u%io)
            if (xy_u%sw .and. myid == xy_u%ionode) close (xy_u%io)
            if (fs_u%sw .and. myid == fs_u%ionode) close (fs_u%io)
            if (ob_u%sw .and. myid == ob_u%ionode) close (ob_u%io)
        else

            if (yz_ps%sw) call wbuf_yz_ps(nt + 1)
            if (xz_ps%sw) call wbuf_xz_ps(nt + 1)
            if (xy_ps%sw) call wbuf_xy_ps(nt + 1)
            if (fs_ps%sw) call wbuf_fs_ps(nt + 1)
            if (ob_ps%sw) call wbuf_ob_ps(nt + 1)

            if (yz_v%sw) call wbuf_yz_v(nt + 1)
            if (xz_v%sw) call wbuf_xz_v(nt + 1)
            if (xy_v%sw) call wbuf_xy_v(nt + 1)
            if (fs_v%sw) call wbuf_fs_v(nt + 1)
            if (ob_v%sw) call wbuf_ob_v(nt + 1)

            if (yz_u%sw) call wbuf_yz_u(nt + 1)
            if (xz_u%sw) call wbuf_xz_u(nt + 1)
            if (xy_u%sw) call wbuf_xy_u(nt + 1)
            if (fs_u%sw) call wbuf_fs_u(nt + 1)
            if (ob_u%sw) call wbuf_ob_u(nt + 1)

            if (yz_ps%sw .and. myid == yz_ps%ionode) call close_nc(yz_ps)
            if (xz_ps%sw .and. myid == xz_ps%ionode) call close_nc(xz_ps)
            if (xy_ps%sw .and. myid == xy_ps%ionode) call close_nc(xy_ps)
            if (fs_ps%sw .and. myid == fs_ps%ionode) call close_nc(fs_ps)
            if (ob_ps%sw .and. myid == ob_ps%ionode) call close_nc(ob_ps)
            if (yz_v%sw .and. myid == yz_v%ionode) call close_nc(yz_v)
            if (xz_v%sw .and. myid == xz_v%ionode) call close_nc(xz_v)
            if (xy_v%sw .and. myid == xy_v%ionode) call close_nc(xy_v)
            if (fs_v%sw) then
                !$acc update self (max_fs_v)
                call output__put_maxval(fs_v, max_fs_v)
                if (myid == fs_v%ionode) then
                    call close_nc(fs_v)
                end if
            end if
            if (ob_v%sw) then
                !$acc update self (max_ob_v)
                call output__put_maxval(ob_v, max_ob_v)
                if (myid == ob_v%ionode) then
                    call close_nc(ob_v)
                end if
            end if

            if (yz_u%sw .and. myid == yz_u%ionode) call close_nc(yz_u)
            if (xz_u%sw .and. myid == xz_u%ionode) call close_nc(xz_u)
            if (xy_u%sw .and. myid == xy_u%ionode) call close_nc(xy_u)
            if (fs_u%sw) then
                !$acc update self (max_fs_u)
                call output__put_maxval(fs_u, max_fs_u)
                if (myid == fs_u%ionode) then
                    call close_nc(fs_u)
                end if
            end if
            if (ob_u%sw) then
                !$acc update self (max_ob_u)
                call output__put_maxval(ob_u, max_ob_u)
                if (myid == ob_u%ionode) then
                    call close_nc(ob_u)
                end if
            end if

        end if

        call pwatch__off('snap__closefiles')

    end subroutine snap__closefiles

    subroutine output__put_maxval(hdr, maxv)

        type(snp), intent(in) :: hdr
        real(SP), intent(inout) :: maxv(nxs, nys, 3)
        real(SP) :: sbuf(nxs * nys * 3), rbuf(nxs * nys * 3)
        integer :: err
        integer :: vid_V, vid_H, vid_A

        sbuf = reshape(maxv, shape(sbuf))
        call mpi_reduce(sbuf, rbuf, nxs * nys * 3, MPI_REAL, MPI_SUM, hdr%ionode, mpi_comm_world, err)
        maxv = reshape(rbuf, shape(maxv))
        if (myid == hdr%ionode) then

            if (snp_format == 'native') then
        !! pass
            else

                call nc_chk(nf90_redef(hdr%io))
                call nc_chk(nf90_def_var(hdr%io, 'max-V', NF90_REAL, (/hdr%did_x1, hdr%did_x2/), vid_V))
                call nc_chk(nf90_def_var(hdr%io, 'max-H', NF90_REAL, (/hdr%did_x1, hdr%did_x2/), vid_H))
                call nc_chk(nf90_def_var(hdr%io, 'max-A', NF90_REAL, (/hdr%did_x1, hdr%did_x2/), vid_A))
                call nc_chk(nf90_put_att(hdr%io, vid_V, 'long_name', 'Maximum amplitude of the vertical component'))
                call nc_chk(nf90_put_att(hdr%io, vid_H, 'long_name', 'Maximum amplitude of the horizontal components'))
                call nc_chk(nf90_put_att(hdr%io, vid_A, 'long_name', 'Maximum amplitude of the vector motion'))
                call nc_chk(nf90_put_att(hdr%io, vid_V, 'coordinates', 'lat lon'))
                call nc_chk(nf90_put_att(hdr%io, vid_H, 'coordinates', 'lat lon'))
                call nc_chk(nf90_put_att(hdr%io, vid_A, 'coordinates', 'lat lon'))

                if (hdr%snaptype == 'v3') then

                    call nc_chk(nf90_put_att(hdr%io, vid_V, 'units', 'm/s'))
                    call nc_chk(nf90_put_att(hdr%io, vid_H, 'units', 'm/s'))
                    call nc_chk(nf90_put_att(hdr%io, vid_A, 'units', 'm/s'))

                else if (hdr%snaptype == 'u3') then

                    call nc_chk(nf90_put_att(hdr%io, vid_V, 'units', 'm'))
                    call nc_chk(nf90_put_att(hdr%io, vid_H, 'units', 'm'))
                    call nc_chk(nf90_put_att(hdr%io, vid_A, 'units', 'm'))

                end if

                call nc_chk(nf90_put_att(hdr%io, vid_V, 'actual_range', (/minval(maxv(:, :, 1)), maxval(maxv(:, :, 1))/)))
                call nc_chk(nf90_put_att(hdr%io, vid_H, 'actual_range', (/minval(maxv(:, :, 2)), maxval(maxv(:, :, 2))/)))
                call nc_chk(nf90_put_att(hdr%io, vid_A, 'actual_range', (/minval(maxv(:, :, 3)), maxval(maxv(:, :, 3))/)))
                call nc_chk(nf90_enddef(hdr%io))
                call nc_chk(nf90_put_var(hdr%io, vid_V, maxv(:, :, 1)))
                call nc_chk(nf90_put_var(hdr%io, vid_H, maxv(:, :, 2)))
                call nc_chk(nf90_put_var(hdr%io, vid_A, maxv(:, :, 3)))

            end if
        end if

    end subroutine output__put_maxval

    subroutine nc_chk(err)

          !! An internal subroutine to check error in netcdf function calls
        integer, intent(in) :: err

        if (err /= NF90_NOERR) write (error_unit, *) NF90_STRERROR(err)

    end subroutine nc_chk

    subroutine pack_3D(n1, n2, n3, buf3d, buf1d)

        integer,  intent(in)  :: n1, n2, n3
        real(SP), intent(in)  :: buf3d(n1,n2,n3)
        real(SP), intent(out) :: buf1d(n1*n2*n3)
        integer :: i, j, k, idx

        idx = 1
#ifdef _OPENACC
        !$acc kernels present(buf3d, buf1d)
        !$acc loop independent collapse(3) 
#else
        !$omp parallel do private(i, j, k, idx)
#endif
        do k=1, n3
            do j=1, n2
                do i=1, n1
                    idx = (k-1) * n1 * n2 + (j-1) * n1 + i
                    buf1d(idx) = buf3d(i,j,k)
                end do
            end do
        end do
#ifdef _OPENACC
        !$acc end kernels
#else
        !$omp end parallel do
#endif

    end subroutine pack_3D

end module m_snap

