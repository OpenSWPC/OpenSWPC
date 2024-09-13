#include "../shared/m_debug.h"
module m_global

    !! global control parameters, shared arrays and MPI communication
    !!
    !! Copyright 2013-2024 Takuto Maeda. All rights reserved. This project is released under the MIT license.

    use m_std
    use m_debug
    use m_fdtool
    use m_pwatch
    use m_system
    use m_daytim
    use m_readini
    use mpi

    implicit none
    public
    save

    public :: global__setup
    public :: global__setup2
    public :: global__comm_vel
    public :: global__comm_stress
    public :: global__getnode

    real(SP)           :: UC = 10.0**(-15)                           !< Conventional -> SI unit for moment tensor of [Nm]
    integer, parameter :: MP = DP                                    !< DP for mixed precision, SP for pure single precision
    integer, parameter :: NM = 3                                     !< Number of memory variables
    integer, parameter :: NBD = 9                                    !< Number of boundary depths to be memorized

    real(MP), allocatable :: Vx(:, :, :), Vy(:, :, :), Vz(:, :, :)             !<  velocity components
    real(MP), allocatable :: Sxx(:, :, :), Syy(:, :, :), Szz(:, :, :)          !<  normal stress components
    real(MP), allocatable :: Syz(:, :, :), Sxz(:, :, :), Sxy(:, :, :)          !<  shear  stress components
    real(SP), allocatable :: Rxx(:, :, :, :), Ryy(:, :, :, :), Rzz(:, :, :, :) !<  memory variables: normal components
    real(SP), allocatable :: Ryz(:, :, :, :), Rxz(:, :, :, :), Rxy(:, :, :, :) !<  memory variables: shear  components
    real(SP), allocatable :: rho(:, :, :), lam(:, :, :), mu(:, :, :)           !<  density and relaxed moduli
    real(SP), allocatable :: bx(:, :, :), by(:, :, :), bz(:, :, :)             !<  inverse density (buoyancy) at velocity grids
    real(SP), allocatable :: muyz(:, :, :), muxz(:, :, :), muxy(:, :, :)       !<  averaged rigidity at shear stress compoments
    real(SP), allocatable :: taup(:, :, :), taus(:, :, :)                      !<  creep/relax time ratio based on tau-method
    real(SP), allocatable :: ts(:)                                            !<  relaxation time of visco-elastic medium

    logical               :: benchmark_mode                           !<  true for fixed parameter run
    logical               :: pw_mode                                  !< Plane wave mode
    logical               :: bf_mode                                  !< Body force soruce mode
    logical               :: green_mode                               !< Green's function computaiton with reciprocity
!  logical :: fullspace_mode

    character(80)         :: title                                    !<  execution title, used in filename and headers
    integer               :: exedate                                  !<  date and time by seconds from 1970/1/1 0:0:0

    integer               :: nx, ny, nz                               !<  space grid number (global)
    integer               :: nt                                       !<  time grid number
    real(MP)              :: dx, dy, dz                               !<  space grid width
    real(SP)              :: dt                                       !<  time  grid width
    real(SP)              :: xbeg, xend                               !<  global coordinate: x start / end
    real(SP)              :: ybeg, yend                               !<  global coordinate: y start / end
    real(SP)              :: zbeg, zend                               !<  global coordinate: z start / end
    real(SP)              :: tbeg, tend                               !<  beggining and ending elapsed time

    real(SP)              :: vmin                                     !<  minimum velocity
    real(SP)              :: vmax                                     !<  maximum velocity
    real(SP)              :: fmax                                     !<  maximum frequency by the source
    real(SP)              :: fcut                                     !<  cut-off frequency by the source

    integer               :: nproc_x                                  !<  process numbers for x/i - direction
    integer               :: nproc_y                                  !<  process numbers for y/j - direction
    integer               :: nproc                                    !<  total   numbers of process
    integer               :: nxp                                      !<  space grid number in the assigned node
    integer               :: nyp                                      !<  space grid number in the assigned node
    integer               :: myid                                     !<  MPI node number
    integer               :: idx, idy                                 !<  2D horizontal division ID
    integer, allocatable  :: itbl(:, :)                               !<  node layout table
    integer               :: ibeg, iend                               !<  i-region in the node
    integer               :: jbeg, jend                               !<  j-region in the node
    integer               :: kbeg, kend                               !<  k-region in the node
    integer               :: ibeg_m, iend_m                           !<  i- memory allocation area
    integer               :: jbeg_m, jend_m                           !<  j- memory allocation area
    integer               :: kbeg_m, kend_m                           !<  k- memory allocation area
    integer               :: ipad, jpad, kpad                         !<  memory padding size for optimization

    integer               :: na                                       !<  absorber thickness
    integer               :: ibeg_k, iend_k                           !<  i- kernel integration area without absorption band
    integer               :: jbeg_k, jend_k                           !<  j- kernel integration area without absorption band
    integer               :: kbeg_k, kend_k                           !<  k- kernel integration area without absorption band
    integer, allocatable  :: kbeg_a(:, :)                             !<  k>=kbeg_a(i,j) is in absorber region
    character(16)         :: abc_type

    real(SP)              :: M0                                       !<  total moment

    character(256) :: odir                                            !<  output directory

    integer, allocatable :: kfs(:, :)                                 !<  free surface depth grid in the node
    integer, allocatable :: kob(:, :)                                 !<  ocean bottom depth grid in the node
    integer, allocatable :: kfs_top(:, :), kfs_bot(:, :)              !<  region in which 2nd-order FDM is applied for free surface
    integer, allocatable :: kob_top(:, :), kob_bot(:, :)              !<  region in which 2nd-order FDM is applied for ocean bottom
    real(SP), allocatable :: bddep(:, :, :)                           !<  boundary depth in physical coordinate

    real(SP) :: clon                                                  !< center longitude
    real(SP) :: clat                                                  !< center latitude
    real(SP) :: phi                                                   !< azimuth
    real(SP), allocatable :: xc(:), yc(:), zc(:)

    real(SP) :: evlo
    real(SP) :: evla
    real(SP) :: evdp !< unit:km
    real(SP) :: mxx0, myy0, mzz0, myz0, mxz0, mxy0
    real(SP) :: fx0, fy0, fz0
    real(SP) :: otim
    real(SP) :: sx0, sy0

    real(MP), private, allocatable :: sbuf_ip(:), sbuf_im(:)          !<  mpi send buffer for x-dir
    real(MP), private, allocatable :: sbuf_jp(:), sbuf_jm(:)          !<  mpi send buffer for y-dir
    real(MP), private, allocatable :: rbuf_ip(:), rbuf_im(:)          !<  mpi recv buffer for x-dir
    real(MP), private, allocatable :: rbuf_jp(:), rbuf_jm(:)          !<  mpi recv buffer for y-dir

    integer, private :: mpi_precision

    private :: inside_node
    private :: set_mpi_table

contains

    subroutine global__readprm(io_prm)

        integer, intent(in) :: io_prm

        call readini(io_prm, 'benchmark_mode', benchmark_mode, .false.)

    !! read parameters
        call readini(io_prm, 'title', title, 'swpc3d')
        call readini(io_prm, 'nproc_x', nproc_x, 1)
        call readini(io_prm, 'nproc_y', nproc_y, 2)
        call readini(io_prm, 'nx', nx, 256)
        call readini(io_prm, 'ny', ny, 256)
        call readini(io_prm, 'nz', nz, 256)
        call readini(io_prm, 'nt', nt, 1000)
        call readini(io_prm, 'ipad', ipad, 0)
        call readini(io_prm, 'jpad', jpad, 0)
        call readini(io_prm, 'kpad', kpad, 0)
        call readini(io_prm, 'odir', odir, './out')

    !! some parameters are fixed for benchmark mode
        if (benchmark_mode) then
            dx = 0.5
            dy = 0.5
            dz = 0.5
            dt = 0.04
            na = 20
            xbeg = -nx / 2.0 * real(dx) ! force x=0 at center
            ybeg = -ny / 2.0 * real(dy) ! force y=0 at center
            zbeg = -30 * real(dz)
            tbeg = 0.0
            clon = 139.7604
            clat = 35.7182
            phi = 0.0
            abc_type = 'pml'
!      fullspace_mode = .false.
        else !! or read from file for regular run
            call readini(io_prm, 'dx', dx, 0.5_MP)
            call readini(io_prm, 'dy', dy, 0.5_MP)
            call readini(io_prm, 'dz', dz, 0.5_MP)
            call readini(io_prm, 'dt', dt, 0.01)
            call readini(io_prm, 'na', na, 20)
            call readini(io_prm, 'xbeg', xbeg, -nx / 2 * real(dx))
            call readini(io_prm, 'ybeg', ybeg, -ny / 2 * real(dy))
            call readini(io_prm, 'zbeg', zbeg, -30 * real(dz))
            call readini(io_prm, 'tbeg', tbeg, 0.0)
            call readini(io_prm, 'clon', clon, 139.7604)
            call readini(io_prm, 'clat', clat, 35.7182)
            call readini(io_prm, 'phi', phi, 0.0)
            call readini(io_prm, 'abc_type', abc_type, 'pml')
!      call readini( io_prm, 'fullspace_mode', fullspace_mode, .false. )

        end if

    end subroutine global__readprm

    !>
  !! read parameter file, memory allocation, MPI set-ups
    !<
  !!
    subroutine global__setup(io_prm)

        integer, intent(in) :: io_prm
        integer :: err

    !!
    !! MPI status check
    !!
        call mpi_comm_rank(mpi_comm_world, myid, err)
        if (MP == DP) then
            mpi_precision = MPI_DOUBLE_PRECISION
        else
            mpi_precision = MPI_REAL
        end if

    !!
    !! read key parameters
    !!
        call global__readprm(io_prm)

        nproc = nproc_x * nproc_y            !! total number of processes

    !!
    !! obtain date by unixtime: seconds measured from 1970/1/1 0:0:0
    !!
        if (myid == 0) call daytim__getdate(exedate)
        call mpi_bcast(exedate, 1, MPI_INTEGER, 0, mpi_comm_world, err)

    !!
    !! derived parameters
    !!
        xend = xbeg + nx * real(dx)
        yend = ybeg + ny * real(dy)
        zend = zbeg + nz * real(dz)
        tend = ybeg + nt * dt

    end subroutine global__setup

    !>
  !! Common parameter setup, needed only for starting
    !<
    subroutine global__setup2

        integer :: i, j, k
        integer :: err
        integer :: nproc_exe
        integer :: mx, my
        integer :: proc_x, proc_y

        call pwatch__on("global__setup2") !! measure from here

    !! create output directory (if it does not exist)
        call system__call('mkdir -p '//trim(odir)//'> /dev/null 2>&1')

    !! size settings
        call mpi_comm_size(mpi_comm_world, nproc_exe, err)
        call assert(nproc == nproc_exe)

        mx = mod(nx, nproc_x)
        my = mod(ny, nproc_y)
        proc_x = mod(myid, nproc_x)
        proc_y = myid / nproc_x
        if (proc_x <= nproc_x - mx - 1) then
            nxp = (nx - mx) / nproc_x
        else
            nxp = (nx - mx) / nproc_x + 1
        end if
        if (proc_y <= nproc_y - my - 1) then
            nyp = (ny - my) / nproc_y
        else
            nyp = (ny - my) / nproc_y + 1
        end if

    !! MPI coordinate

        allocate (itbl(-1:nproc_x, -1:nproc_y))
        allocate (sbuf_ip(5 * nyp * nz), source=0.0_MP)
        allocate (sbuf_im(5 * nyp * nz), source=0.0_MP)
        allocate (rbuf_ip(5 * nyp * nz), source=0.0_MP)
        allocate (rbuf_im(5 * nyp * nz), source=0.0_MP)
        allocate (sbuf_jp(5 * nxp * nz), source=0.0_MP)
        allocate (sbuf_jm(5 * nxp * nz), source=0.0_MP)
        allocate (rbuf_jp(5 * nxp * nz), source=0.0_MP)
        allocate (rbuf_jm(5 * nxp * nz), source=0.0_MP)

    !! MPI communication table
        call set_mpi_table

    !! computation region in this node (#244)
        if (proc_x <= nproc_x - mx - 1) then
            ibeg = proc_x * (nx - mx) / nproc_x + 1
            iend = (proc_x + 1) * (nx - mx) / nproc_x
        else
            ibeg = proc_x * ((nx - mx) / nproc_x + 1) - (nproc_x - mx) + 1
            iend = (proc_x + 1) * ((nx - mx) / nproc_x + 1) - (nproc_x - mx)
        end if
        if (proc_y <= nproc_y - my - 1) then
            jbeg = proc_y * (ny - my) / nproc_y + 1
            jend = (proc_y + 1) * (ny - my) / nproc_y
        else
            jbeg = proc_y * ((ny - my) / nproc_y + 1) - (nproc_y - my) + 1
            jend = (proc_y + 1) * ((ny - my) / nproc_y + 1) - (nproc_y - my)
        end if

        kbeg = 1
        kend = nz

    !! memory requirements including margin for MPI/boundary conditions
    !! stress glut also requires sleeve area
        ibeg_m = ibeg - 3
        iend_m = iend + 3 + ipad
        jbeg_m = jbeg - 3
        jend_m = jend + 3 + jpad
        kbeg_m = kbeg - 3
        kend_m = kend + 3 + kpad

    !! coordinate setting
        allocate (xc(ibeg_m:iend_m), yc(jbeg_m:jend_m), zc(kbeg_m:kend_m))

        do i = ibeg_m, iend_m
            xc(i) = i2x(i, xbeg, real(dx))
        end do
        do j = jbeg_m, jend_m
            yc(j) = j2y(j, ybeg, real(dy))
        end do
        do k = kbeg_m, kend_m
            zc(k) = k2z(k, zbeg, real(dz))
        end do

    !!
    !! absorbing boundary region                     -+---> i,j(x,y)
        !                                                 |
    !!  +-----+--------------------------+-----+      |
    !!  |     |                          |     |      v k(z)
    !!  |     |                          |     |
    !!  |     |                          |     |
    !!  |     |                          |     |
    !!  |     |  interior region         |     |
    !!  |     |  eveluated by m_kernel   |     |
    !!  |     |                          |     |
    !!  |     |                          |     |
    !!  |     |                          |     |
    !!  |     +--------------------------+     |
    !!  |         exterior region              |
    !!  |         evaluated by m_absorb        |
    !!  +-----+--------------------------+-----+
    !!  1      na                        nx-na+1  nx
    !!  <- na ->

        allocate (kbeg_a(ibeg_m:iend_m, jbeg_m:jend_m))
        do j = jbeg_m, jend_m
            do i = ibeg_m, iend_m
                if (i <= na .or. nx - na + 1 <= i .or. j <= na .or. ny - na + 1 <= j) then
                    kbeg_a(i, j) = kbeg
                else
                    kbeg_a(i, j) = kend - na + 1
                end if
            end do
        end do

    !! Interior Kernel region

    !! initial value
        ibeg_k = ibeg
        iend_k = iend
        jbeg_k = jbeg
        jend_k = jend
        kbeg_k = kbeg
        kend_k = kend

        if (abc_type == 'pml') then

!      if( fullspace_mode ) kbeg_k = na + 1

            if (iend <= na) then; ibeg_k = iend + 1; ! no kernel integration
            else if (ibeg <= na) then; ibeg_k = na + 1; ! pertial kernel
            end if

            if (ibeg >= nx - na + 1) then; iend_k = ibeg - 1; ! no kernel integartion
            else if (iend >= nx - na + 1) then; iend_k = nx - na; 
            end if

            if (jend <= na) then; jbeg_k = jend + 1; ! no kernel integration
            else if (jbeg <= na) then; jbeg_k = na + 1; ! pertial kernel
            end if

            if (jbeg >= ny - na + 1) then; jend_k = jbeg - 1; ! no kernel integartion
            else if (jend >= ny - na + 1) then; jend_k = ny - na; 
            end if
            kend_k = nz - na

        end if

        call pwatch__off("global__setup2") !! measure from here

    end subroutine global__setup2

    !>
  !! Data buffring & communication for velocity vector
    !<
  !!
    subroutine global__comm_vel()

        integer :: isize, jsize
        integer :: err
        integer :: istatus(mpi_status_size, 4)
        integer :: req_i(4), req_j(4)

        if (myid >= nproc) return

        call pwatch__on("global__comm_vel")

    !! unit buffer size
        isize = nyp * nz
        jsize = nxp * nz

        call mpi_irecv(rbuf_ip, 5 * isize, mpi_precision, itbl(idx + 1, idy), 1, mpi_comm_world, req_i(1), err)
        call mpi_irecv(rbuf_im, 4 * isize, mpi_precision, itbl(idx - 1, idy), 2, mpi_comm_world, req_i(2), err)
        call mpi_irecv(rbuf_jp, 5 * jsize, mpi_precision, itbl(idx, idy + 1), 3, mpi_comm_world, req_j(1), err)
        call mpi_irecv(rbuf_jm, 4 * jsize, mpi_precision, itbl(idx, idy - 1), 4, mpi_comm_world, req_j(2), err)

        sbuf_ip(1:2 * isize) = reshape(Vx(1:nz, iend - 1:iend, jbeg:jend), (/2 * isize/))
        sbuf_ip(2 * isize + 1:3 * isize) = reshape(Vy(1:nz, iend:iend, jbeg:jend), (/isize/))
        sbuf_ip(3 * isize + 1:4 * isize) = reshape(Vz(1:nz, iend:iend, jbeg:jend), (/isize/))
        call mpi_isend(sbuf_ip, 4 * isize, mpi_precision, itbl(idx + 1, idy), 2, mpi_comm_world, req_i(3), err)

        sbuf_im(1:isize) = reshape(Vx(1:nz, ibeg:ibeg, jbeg:jend), (/isize/))
        sbuf_im(isize + 1:3 * isize) = reshape(Vy(1:nz, ibeg:ibeg + 1, jbeg:jend), (/2 * isize/))
        sbuf_im(3 * isize + 1:5 * isize) = reshape(Vz(1:nz, ibeg:ibeg + 1, jbeg:jend), (/2 * isize/))
        call mpi_isend(sbuf_im, 5 * isize, mpi_precision, itbl(idx - 1, idy), 1, mpi_comm_world, req_i(4), err)

        sbuf_jp(1:jsize) = reshape(Vx(1:nz, ibeg:iend, jend:jend), (/jsize/))
        sbuf_jp(jsize + 1:3 * jsize) = reshape(Vy(1:nz, ibeg:iend, jend - 1:jend), (/2 * jsize/))
        sbuf_jp(3 * jsize + 1:4 * jsize) = reshape(Vz(1:nz, ibeg:iend, jend:jend), (/jsize/))
        call mpi_isend(sbuf_jp, 4 * jsize, mpi_precision, itbl(idx, idy + 1), 4, mpi_comm_world, req_j(3), err)

        sbuf_jm(1:2 * jsize) = reshape(Vx(1:nz, ibeg:iend, jbeg:jbeg + 1), (/2 * jsize/))
        sbuf_jm(2 * jsize + 1:3 * jsize) = reshape(Vy(1:nz, ibeg:iend, jbeg:jbeg), (/jsize/))
        sbuf_jm(3 * jsize + 1:5 * jsize) = reshape(Vz(1:nz, ibeg:iend, jbeg:jbeg + 1), (/2 * jsize/))
        call mpi_isend(sbuf_jm, 5 * jsize, mpi_precision, itbl(idx, idy - 1), 3, mpi_comm_world, req_j(4), err)

        call mpi_waitall(4, req_i, istatus, err)

        Vx(1:nz, ibeg - 2:ibeg - 1, jbeg:jend) = reshape(rbuf_im(1:2 * isize), (/nz, 2, nyp/))
        Vy(1:nz, ibeg - 1:ibeg - 1, jbeg:jend) = reshape(rbuf_im(2 * isize + 1:3 * isize), (/nz, 1, nyp/))
        Vz(1:nz, ibeg - 1:ibeg - 1, jbeg:jend) = reshape(rbuf_im(3 * isize + 1:4 * isize), (/nz, 1, nyp/))

        Vx(1:nz, iend + 1:iend + 1, jbeg:jend) = reshape(rbuf_ip(1:isize), (/nz, 1, nyp/))
        Vy(1:nz, iend + 1:iend + 2, jbeg:jend) = reshape(rbuf_ip(isize + 1:3 * isize), (/nz, 2, nyp/))
        Vz(1:nz, iend + 1:iend + 2, jbeg:jend) = reshape(rbuf_ip(3 * isize + 1:5 * isize), (/nz, 2, nyp/))

        call mpi_waitall(4, req_j, istatus, err)

        Vx(1:nz, ibeg:iend, jbeg - 1:jbeg - 1) = reshape(rbuf_jm(1:jsize), (/nz, nxp, 1/))
        Vy(1:nz, ibeg:iend, jbeg - 2:jbeg - 1) = reshape(rbuf_jm(jsize + 1:3 * jsize), (/nz, nxp, 2/))
        Vz(1:nz, ibeg:iend, jbeg - 1:jbeg - 1) = reshape(rbuf_jm(3 * jsize + 1:4 * jsize), (/nz, nxp, 1/))

        Vx(1:nz, ibeg:iend, jend + 1:jend + 2) = reshape(rbuf_jp(1:2 * jsize), (/nz, nxp, 2/))
        Vy(1:nz, ibeg:iend, jend + 1:jend + 1) = reshape(rbuf_jp(2 * jsize + 1:3 * jsize), (/nz, nxp, 1/))
        Vz(1:nz, ibeg:iend, jend + 1:jend + 2) = reshape(rbuf_jp(3 * jsize + 1:5 * jsize), (/nz, nxp, 2/))

        call pwatch__off("global__comm_vel")

    end subroutine global__comm_vel

    !>
  !! Data buffring & communication for stress tensor
    !<
  !!
    subroutine global__comm_stress()

        integer :: isize, jsize
        integer :: err
        integer :: istatus(mpi_status_size, 4)
        integer :: req_i(4), req_j(4)

        if (myid >= nproc) return

        call pwatch__on("global__comm_stress")

    !! unit buffer size
        isize = nyp * nz
        jsize = nxp * nz

        call mpi_irecv(rbuf_ip, 4 * isize, mpi_precision, itbl(idx + 1, idy), 5, mpi_comm_world, req_i(1), err)
        call mpi_irecv(rbuf_im, 5 * isize, mpi_precision, itbl(idx - 1, idy), 6, mpi_comm_world, req_i(2), err)
        call mpi_irecv(rbuf_jp, 4 * jsize, mpi_precision, itbl(idx, idy + 1), 7, mpi_comm_world, req_j(1), err)
        call mpi_irecv(rbuf_jm, 5 * jsize, mpi_precision, itbl(idx, idy - 1), 8, mpi_comm_world, req_j(2), err)

        sbuf_ip(1:isize) = reshape(Sxx(1:nz, iend:iend, jbeg:jend), (/isize/))
        sbuf_ip(isize + 1:3 * isize) = reshape(Sxy(1:nz, iend - 1:iend, jbeg:jend), (/2 * isize/))
        sbuf_ip(3 * isize + 1:5 * isize) = reshape(Sxz(1:nz, iend - 1:iend, jbeg:jend), (/2 * isize/))
        call mpi_isend(sbuf_ip, 5 * isize, mpi_precision, itbl(idx + 1, idy), 6, mpi_comm_world, req_i(3), err)

        sbuf_im(1:2 * isize) = reshape(Sxx(1:nz, ibeg:ibeg + 1, jbeg:jend), (/2 * isize/))
        sbuf_im(2 * isize + 1:3 * isize) = reshape(Sxy(1:nz, ibeg:ibeg, jbeg:jend), (/isize/))
        sbuf_im(3 * isize + 1:4 * isize) = reshape(Sxz(1:nz, ibeg:ibeg, jbeg:jend), (/isize/))
        call mpi_isend(sbuf_im, 4 * isize, mpi_precision, itbl(idx - 1, idy), 5, mpi_comm_world, req_i(4), err)

        sbuf_jp(1:jsize) = reshape(Syy(1:nz, ibeg:iend, jend:jend), (/jsize/))
        sbuf_jp(jsize + 1:3 * jsize) = reshape(Sxy(1:nz, ibeg:iend, jend - 1:jend), (/2 * jsize/))
        sbuf_jp(3 * jsize + 1:5 * jsize) = reshape(Syz(1:nz, ibeg:iend, jend - 1:jend), (/2 * jsize/))
        call mpi_isend(sbuf_jp, 5 * jsize, mpi_precision, itbl(idx, idy + 1), 8, mpi_comm_world, req_j(3), err)

        sbuf_jm(1:2 * jsize) = reshape(Syy(1:nz, ibeg:iend, jbeg:jbeg + 1), (/2 * jsize/))
        sbuf_jm(2 * jsize + 1:3 * jsize) = reshape(Sxy(1:nz, ibeg:iend, jbeg:jbeg), (/jsize/))
        sbuf_jm(3 * jsize + 1:4 * jsize) = reshape(Syz(1:nz, ibeg:iend, jbeg:jbeg), (/jsize/))
        call mpi_isend(sbuf_jm, 4 * jsize, mpi_precision, itbl(idx, idy - 1), 7, mpi_comm_world, req_j(4), err)

        call mpi_waitall(4, req_i, istatus, err)

        Sxx(1:nz, ibeg - 1:ibeg - 1, jbeg:jend) = reshape(rbuf_im(1:isize), (/nz, 1, nyp/))
        Sxy(1:nz, ibeg - 2:ibeg - 1, jbeg:jend) = reshape(rbuf_im(isize + 1:3 * isize), (/nz, 2, nyp/))
        Sxz(1:nz, ibeg - 2:ibeg - 1, jbeg:jend) = reshape(rbuf_im(3 * isize + 1:5 * isize), (/nz, 2, nyp/))

        Sxx(1:nz, iend + 1:iend + 2, jbeg:jend) = reshape(rbuf_ip(1:2 * isize), (/nz, 2, nyp/))
        Sxy(1:nz, iend + 1:iend + 1, jbeg:jend) = reshape(rbuf_ip(2 * isize + 1:3 * isize), (/nz, 1, nyp/))
        Sxz(1:nz, iend + 1:iend + 1, jbeg:jend) = reshape(rbuf_ip(3 * isize + 1:4 * isize), (/nz, 1, nyp/))

        call mpi_waitall(4, req_j, istatus, err)

        Syy(1:nz, ibeg:iend, jbeg - 1:jbeg - 1) = reshape(rbuf_jm(1:jsize), (/nz, nxp, 1/))
        Sxy(1:nz, ibeg:iend, jbeg - 2:jbeg - 1) = reshape(rbuf_jm(jsize + 1:3 * jsize), (/nz, nxp, 2/))
        Syz(1:nz, ibeg:iend, jbeg - 2:jbeg - 1) = reshape(rbuf_jm(3 * jsize + 1:5 * jsize), (/nz, nxp, 2/))

        Syy(1:nz, ibeg:iend, jend + 1:jend + 2) = reshape(rbuf_jp(1:2 * jsize), (/nz, nxp, 2/))
        Sxy(1:nz, ibeg:iend, jend + 1:jend + 1) = reshape(rbuf_jp(2 * jsize + 1:3 * jsize), (/nz, nxp, 1/))
        Syz(1:nz, ibeg:iend, jend + 1:jend + 1) = reshape(rbuf_jp(3 * jsize + 1:4 * jsize), (/nz, nxp, 1/))

        call pwatch__off("global__comm_stress")

    end subroutine global__comm_stress

    !>
  !! check if the voxcel location is inside the MPI node
    !<
    logical function inside_node(i, j, k)

        integer, intent(in) :: i, j, k

        if (ibeg <= i .and. i <= iend .and. &
            jbeg <= j .and. j <= jend .and. &
            kbeg <= k .and. k <= kend) then
            inside_node = .true.
        else
            inside_node = .false.
        end if

    end function inside_node

    subroutine set_mpi_table
        integer :: i
        integer :: ii, jj
    !! 2D communication table
        !> @see 2013-0439

        itbl(-1:nproc_x, -1:nproc_y) = MPI_PROC_NULL
        do i = 0, nproc - 1

            ii = mod(i, nproc_x)
            jj = i / nproc_x

            itbl(ii, jj) = i
        end do

    !! location of this process
        idx = mod(myid, nproc_x)
        idy = myid / nproc_x
    end subroutine set_mpi_table

    subroutine global__getnode(i, j, idx_ij, idy_ij)

        integer, intent(in) :: i, j
        integer, intent(out) :: idx_ij, idy_ij
        integer :: mx, my
        integer :: iproc_x, ibeg_node, iend_node
        integer :: iproc_y, jbeg_node, jend_node

        mx = mod(nx, nproc_x)
        my = mod(ny, nproc_y)

        do iproc_x = 0, nproc_x - 1
            if (iproc_x <= nproc_x - mx - 1) then
                ibeg_node = iproc_x * (nx - mx) / nproc_x + 1
                iend_node = (iproc_x + 1) * (nx - mx) / nproc_x
                if (ibeg_node <= i .and. i <= iend_node) then
                    idx_ij = iproc_x
                    exit
                end if
            end if
        end do

        do iproc_y = 0, nproc_y - 1
            if (iproc_y <= nproc_y - my - 1) then
                jbeg_node = iproc_y * (ny - my) / nproc_y + 1
                jend_node = (iproc_y + 1) * (ny - my) / nproc_y
                if (jbeg_node <= j .and. j <= jend_node) then
                    idy_ij = iproc_y
                    exit
                end if
            end if
        end do
    end subroutine global__getnode
end module m_global
!! ------------------------------------------------------------------------------------------------------------------------------ !!
