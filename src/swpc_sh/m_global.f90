#include "../shared/m_debug.h"
module m_global

    ! global control parameters, shared arrays and MPI communication
    !!
    ! Copyright 2013-2025 Takuto Maeda. All rights reserved. This project is released under the MIT license.

    use m_std
    use m_debug
    use m_fdtool
    use m_pwatch
    use m_daytim
    use m_readini
    use mpi
#ifdef _OPENACC
    use openacc
#endif 

    implicit none
    private
    save

    public :: global__setup
    public :: global__setup2
    public :: global__comm_vel
    public :: global__comm_stress

    real(SP), public           :: UC = 10.0**(-12)                   !< Conventional -> SI unit for moment tensor
    integer, parameter, public :: MP = DP                            !< DP for mixed precision, SP for pure single precision
    integer, parameter, public :: NM = 3                             !< Number of memory variables
    integer, parameter, public :: NBD = 9                            !< Number of boundary depths to be memorized

    real(MP), allocatable, public :: Vy(:, :)                        !<  velocity components
    real(MP), allocatable, public :: Syz(:, :), Sxy(:, :)            !<  shear  stress components
    real(SP), allocatable, public :: Ryz(:, :, :), Rxy(:, :, :)      !<  memory variables: shear  components
    real(SP), allocatable, public :: rho(:, :), lam(:, :), mu(:, :)  !<  density, relaxed moduli
    real(SP), allocatable, public :: taup(:, :), taus(:, :)          !<  creep/relax time ratio based on tau-method for attenuation
    real(SP), allocatable, public :: ts(:)                           !<  relaxation time of visco-elastic medium

    character(80), public         :: title                           !<  execution title, used in filename and headers
    integer, public               :: exedate                         !<  date and time by seconds from 1970/1/1 0:0:0

    integer, public               :: nx, nz                          !<  space grid number (global)
    integer, public               :: nt                              !<  time grid number
    integer, public               :: na                              !<  absorber thickness
    real(MP), public              :: dx, dz                          !<  space grid width
    real(SP), public              :: dt                              !<  time  grid width
    real(SP), public              :: xbeg, xend                      !<  global coordinate: x start / end
    real(SP), public              :: zbeg, zend                      !<  global coordinate: z start / end
    real(SP), public              :: tbeg, tend                      !<  beggining and ending elapsed time

    real(SP), public              :: vmin                            !<  minimum velocity
    real(SP), public              :: vmax                            !<  maximum velocity
    real(SP), public              :: fmax                            !<  maximum frequency by the source
    real(SP), public              :: fcut                            !<  cut-off frequency by the source

    integer, public               :: nproc_x                         !<  total   numbers of process
    integer, public               :: nxp                             !<  space grid number in the assigned node
    integer, public               :: myid                            !<  MPI node number
    integer, public               :: idx                             !<  2D horizontal division ID
    integer, allocatable, public  :: itbl(:)                         !<  node layout table
    integer, public               :: ibeg, iend                      !<  i-region in the node
    integer, public               :: kbeg, kend                      !<  k-region in the node
    integer, public               :: ibeg_m, iend_m                  !<  i- memory allocation area
    integer, public               :: kbeg_m, kend_m                  !<  k- memory allocation area
    integer, public               :: ibeg_k, iend_k                  !<  i- kernel integration area without absorption band
    integer, public               :: kbeg_k, kend_k                  !<  k- kernel integration area without absorption band
    integer, allocatable, public  :: kbeg_a(:)                       !<  absorbing boundary region
    integer, public               :: ipad, kpad                      !<  memory padding size for optimization
    character(16), public         :: abc_type                        !<  cerjan or pml

    real(SP), public              :: M0                              !<  total moment

    character(256), public :: odir

    integer, allocatable, public :: kfs(:)                           !<  free surface depth grid in the node
    integer, allocatable, public :: kob(:)                           !<  ocean bottom depth grid in the node
    integer, allocatable, public :: kfs_top(:), kfs_bot(:)           !<  region in which 2nd-order FDM is applied for free surface
    integer, allocatable, public :: kob_top(:), kob_bot(:)           !<  region in which 2nd-order FDM is applied for ocean bottom
    real(SP), allocatable, public :: bddep(:, :)

    real(SP), public :: clon                                         !< center longitude
    real(SP), public :: clat                                         !< center latitude
    real(SP), public :: phi                                          !< azimuth
    real(SP), allocatable, public :: xc(:), zc(:)

    real(SP), public :: evlo
    real(SP), public :: evla
    real(SP), public :: evdp !< unit:km
    real(SP), public :: mxx0, myy0, mzz0, myz0, mxz0, mxy0
    real(SP), public :: fx0, fy0, fz0
    real(SP), public :: otim
    real(SP), public :: sx0, sy0

    logical, public :: benchmark_mode
    logical, public :: pw_mode                                       !< plane wave mode
    logical, public :: bf_mode                                       !< body force mode

    real(MP), private, allocatable :: sbuf_ip(:), sbuf_im(:)         !<  mpi send buffer
    real(MP), private, allocatable :: rbuf_ip(:), rbuf_im(:)         !<  mpi recv buffer
    integer, private :: mpi_precision

  ! fullspace-mdoe
!  logical :: fullspace_mode

    private :: inside_node
    private :: set_mpi_table

contains

    subroutine global__readprm(io_prm)
        integer, intent(in) :: io_prm

        call readini(io_prm, 'benchmark_mode', benchmark_mode, .false.)

        call readini(io_prm, 'title', title, 'swpc_sh')
        call readini(io_prm, 'nproc_x', nproc_x, 1)
        call readini(io_prm, 'nx', nx, 256)
        call readini(io_prm, 'nz', nz, 256)
        call readini(io_prm, 'nt', nt, 1000)
        call readini(io_prm, 'ipad', ipad, 0)
        call readini(io_prm, 'kpad', kpad, 0)
        call readini(io_prm, 'odir', odir, './out')

        if (benchmark_mode) then
            dx = 0.5
            dz = 0.5
            dt = 0.04
            na = 20
            xbeg = -nx / 2 * real(dx)
            zbeg = -30 * real(dz)
            tbeg = 0.0
            clon = 139.7604
            clat = 35.7182
            phi = 0.0
            abc_type = 'pml'
!      fullspace_mode = .false.
        else
            call readini(io_prm, 'dx', dx, 0.5_MP)
            call readini(io_prm, 'dz', dz, 0.5_MP)
            call readini(io_prm, 'dt', dt, 0.01)
            call readini(io_prm, 'na', na, 20)
            call readini(io_prm, 'xbeg', xbeg, -nx / 2 * real(dx))
            call readini(io_prm, 'zbeg', zbeg, -30 * real(dz))
            call readini(io_prm, 'tbeg', tbeg, 0.0)
            call readini(io_prm, 'clon', clon, 139.7604)
            call readini(io_prm, 'clat', clat, 35.7182)
            call readini(io_prm, 'phi', phi, 0.0)
            call readini(io_prm, 'abc_type', abc_type, 'pml')
!      call readini( io_prm, 'fullspace_mode', fullspace_mode, .false.    )
        end if

    end subroutine global__readprm

    subroutine global__setup(io_prm)

        ! read parameter file, memory allocation, MPI set-ups

        integer, intent(in) :: io_prm
        integer :: err

        call pwatch__on("global__setup") ! measure from here

        ! MPI status check
        call mpi_comm_rank(mpi_comm_world, myid, err)
        if (MP == DP) then
            mpi_precision = MPI_DOUBLE_PRECISION
        else
            mpi_precision = MPI_REAL
        end if

        ! read key parameters
        call global__readprm(io_prm)

        ! obtain date by unixtime: seconds measured from 1970/1/1 0:0:0
        if (myid == 0) call daytim__getdate(exedate)
        call mpi_bcast(exedate, 1, MPI_INTEGER, 0, mpi_comm_world, err)

        ! derived parameters
        xend = xbeg + nx * real(dx)
        zend = zbeg + nz * real(dz)
        tend = tbeg + nt * dt

#ifdef _OPENACC
        block
            integer :: ngpus
            ngpus = acc_get_num_devices(acc_device_nvidia)
            call acc_set_device_num(mod(myid, ngpus), acc_device_nvidia)
        end block
#endif                


        call pwatch__off("global__setup") ! measure from here

    end subroutine global__setup

    subroutine global__setup2

        ! Common parameter setup, needed only for starting

        integer :: i, k
        integer :: err, nproc_exe
        integer :: mx, proc_x
        character(256) :: command

        call pwatch__on("global__setup2") ! measure from here


        call mpi_comm_size(mpi_comm_world, nproc_exe, err)
        call assert(nproc_x == nproc_exe)

        mx = mod(nx, nproc_x)
        proc_x = myid
        if (proc_x <= nproc_x - mx + 1) then
            nxp = (nx - mx) / nproc_x
        else
            nxp = (nx - mx) / nproc_x + 1
        end if

        ! MPI coordinate
        allocate (itbl(-1:nproc_x))

        allocate (sbuf_ip(2 * nz), source=0.0_MP)
        allocate (sbuf_im(2 * nz), source=0.0_MP)
        allocate (rbuf_ip(2 * nz), source=0.0_MP)
        allocate (rbuf_im(2 * nz), source=0.0_MP)

        ! MPI communication table
        call set_mpi_table

        ! create output directory (if it does not exist)
        call mpi_barrier(mpi_comm_world, err)
        command = 'if [ ! -d '// trim(odir) // ' ]; then mkdir -p ' &
                 // trim(odir) // ' > /dev/null 2>&1 ; fi'
        do i=0, nproc_exe-1
            if (myid == i) then
                call execute_command_line(trim(command))
            end if
            call mpi_barrier(mpi_comm_world, err)
        end do                

        ! computation region in this node (#244)
        if (proc_x <= nproc_x - mx - 1) then
            ibeg = proc_x * (nx - mx) / nproc_x + 1
            iend = (proc_x + 1) * (nx - mx) / nproc_x
        else
            ibeg = proc_x * ((nx - mx) / nproc_x + 1) - (nproc_x - mx) + 1
            iend = (proc_x + 1) * ((nx - mx) / nproc_x + 1) - (nproc_x - mx)
        end if
        kbeg = 1
        kend = nz

        ! memory requirements including margin for MPI/boundary conditions
        ! stress glut also requires sleeve area
        ibeg_m = ibeg - 3
        iend_m = iend + 3 + ipad
        kbeg_m = kbeg - 3
        kend_m = kend + 3 + kpad

        ! coordinate setting
        allocate (xc(ibeg_m:iend_m), zc(kbeg_m:kend_m))

        do i = ibeg_m, iend_m
            xc(i) = i2x(i, xbeg, real(dx))
        end do
        do k = kbeg_m, kend_m
            zc(k) = k2z(k, zbeg, real(dz))
        end do

        ! Absorbing region definition
        allocate (kbeg_a(ibeg_m:iend_m))
        do i = ibeg_m, iend_m
            if (i <= na .or. nx - na + 1 <= i) then
                kbeg_a(i) = kbeg
            else
                kbeg_a(i) = kend - na + 1
            end if
        end do

        ibeg_k = ibeg
        iend_k = iend
        kbeg_k = kbeg
        kend_k = kend

        if (abc_type == 'pml') then
!      if( fullspace_mode ) kbeg_k = na+1
            if (iend <= na) then ! no kernel integration
                ibeg_k = iend + 1
            else if (ibeg <= na) then ! pertial kernel
                ibeg_k = na + 1
            end if
            if (ibeg >= nx - na + 1) then ! no kernel integartion
                iend_k = ibeg - 1
            else if (iend >= nx - na + 1) then
                iend_k = nx - na
            end if
            kend_k = nz - na
        end if

        !$acc enter data copyin(&
        !$acc sbuf_ip, sbuf_im, rbuf_ip, rbuf_im, itbl)


        call pwatch__off("global__setup2") ! measure from here

    end subroutine global__setup2

    subroutine global__comm_vel()

        ! Data buffring & communication for velocity vector ( 2013-0420, 2013-0421)

        integer :: err, k
        integer :: istatus(mpi_status_size, 4)
        integer :: req(4)

        if (myid >= nproc_x) return

        call pwatch__on("global__comm_vel")

        !$acc host_data use_device(rbuf_ip, rbuf_im)
        call mpi_irecv(rbuf_ip, 2 * nz, mpi_precision, itbl(idx + 1), 1, mpi_comm_world, req(1), err)
        call mpi_irecv(rbuf_im, 1 * nz, mpi_precision, itbl(idx - 1), 2, mpi_comm_world, req(2), err)
        !$acc end host_data

        !$acc kernels present(Vy, sbuf_ip, sbuf_im)
        !$acc loop independent
        do k=1, nz
            sbuf_ip(k) = Vy(k, iend)
            sbuf_im(k) = Vy(k, ibeg)
            sbuf_im(nz+k) = Vy(k, ibeg+1)
        end do
        !$acc end kernels

        !$acc host_data use_device(sbuf_ip, sbuf_im)
        call mpi_isend(sbuf_ip, 1 * nz, mpi_precision, itbl(idx + 1), 2, mpi_comm_world, req(3), err)
        call mpi_isend(sbuf_im, 2 * nz, mpi_precision, itbl(idx - 1), 1, mpi_comm_world, req(4), err)
        !$acc end host_data

        call mpi_waitall(4, req, istatus, err)

        !$acc kernels present(Vy, rbuf_ip, rbuf_im)
        !$acc loop independent
        do k=1, nz
            Vy(k,iend+1) = rbuf_ip(k)
            Vy(k,iend+2) = rbuf_ip(nz+k)
            Vy(k,ibeg-1) = rbuf_im(k)
        end do
        !$acc end kernels

        call pwatch__off("global__comm_vel")

    end subroutine global__comm_vel

    subroutine global__comm_stress()

        ! Data buffring & communication for stress tensor

        integer :: err, k
        integer :: istatus(mpi_status_size, 4)
        integer :: req(4)

        if (myid >= nproc_x) return

        call pwatch__on("global__comm_stress")

        !$acc host_data use_device(rbuf_ip, rbuf_im)
        call mpi_irecv(rbuf_ip, 1 * nz, mpi_precision, itbl(idx + 1), 3, mpi_comm_world, req(1), err)
        call mpi_irecv(rbuf_im, 2 * nz, mpi_precision, itbl(idx - 1), 4, mpi_comm_world, req(2), err)
        !$acc end host_data

        !$acc kernels present(Sxy, sbuf_ip, sbuf_im)
        !$acc loop independent
        do k=1, nz
            sbuf_ip(k) = Sxy(k,iend-1)
            sbuf_ip(nz+k) = Sxy(k,iend)
            sbuf_im(k) = Sxy(k,ibeg)
        end do
        !$acc end kernels

        !$acc host_data use_device(sbuf_ip, sbuf_im)
        call mpi_isend(sbuf_ip, 2 * nz, mpi_precision, itbl(idx + 1), 4, mpi_comm_world, req(3), err)
        call mpi_isend(sbuf_im, 1 * nz, mpi_precision, itbl(idx - 1), 3, mpi_comm_world, req(4), err)
        !$acc end host_data

        call mpi_waitall(4, req, istatus, err)

        ! Resore the data
        !$acc kernels present(Sxy, rbuf_ip, rbuf_im)
        !$acc loop independent
        do k=1, nz
            Sxy(k,iend+1) = rbuf_ip(k)
            Sxy(k,ibeg-2) = rbuf_ip(nz+k)
            Sxy(k,ibeg-1) = rbuf_ip(2*nz+k)
        end do
        !$acc end kernels

        call pwatch__off("global__comm_stress")

    end subroutine global__comm_stress

    logical function inside_node(i, k)

        ! check if the voxcel location is inside the MPI node

        integer, intent(in) :: i, k

        if (ibeg <= i .and. i <= iend .and. &
            kbeg <= k .and. k <= kend) then
            inside_node = .true.
        else
            inside_node = .false.
        end if

    end function inside_node

    subroutine set_mpi_table

        ! 2D communication table (2013-0439)

        integer :: i
        integer :: ii

        itbl(-1:nproc_x) = MPI_PROC_NULL
        do i = 0, nproc_x - 1

            ii = mod(i, nproc_x)
            itbl(ii) = i

        end do

        ! location of this process
        idx = mod(myid, nproc_x)
    end subroutine set_mpi_table

end module m_global
