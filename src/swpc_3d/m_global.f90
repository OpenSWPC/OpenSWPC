#include "../shared/m_debug.h"
module m_global

    ! global control parameters, shared arrays and MPI communication
    ! Copyright 2013-2025 Takuto Maeda. All rights reserved. This project is released under the MIT license.

    use m_std
    use m_debug
    use m_fdtool
    use m_pwatch
    use m_daytim
    use m_readini
    use m_daytim
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
    public :: global__getnode

    real(SP), public :: UC = 10.0**(-15)                                        !< Conventional -> SI unit 
    integer, parameter, public :: MP = DP                                      !< DP for mixed, SP for single precisions
    integer, parameter, public :: NM = 3                                       !< Number of memory variables
    integer, parameter, public :: NBD = 9                                      !< Number of boundary depths to be memorized

    real(MP), allocatable, public :: Vx(:,:,:), Vy(:,:,:), Vz(:,:,:)           !< velocity components
    real(MP), allocatable, public :: Sxx(:,:,:), Syy(:,:,:), Szz(:,:,:)        !< normal stress components
    real(MP), allocatable, public :: Syz(:,:,:), Sxz(:,:,:), Sxy(:,:,:)        !< shear  stress components
    real(SP), allocatable, public :: Rxx(:,:,:,:), Ryy(:,:,:,:), Rzz(:,:,:,:)  !< memory variables: normal components
    real(SP), allocatable, public :: Ryz(:,:,:,:), Rxz(:,:,:,:), Rxy(:,:,:,:)  !< memory variables: shear  components
    real(SP), allocatable, public :: rho(:,:,:), lam(:,:,:), mu(:,:,:)         !< density and relaxed moduli
    real(SP), allocatable, public :: bx(:,:,:), by(:,:,:), bz(:,:,:)           !< inverse density (buoyancy) at velocity grids
    real(SP), allocatable, public :: muyz(:,:,:), muxz(:,:,:), muxy(:,:,:)     !< averaged rigidity at shear stress compoments
    real(SP), allocatable, public :: taup(:,:,:), taus(:,:,:)                  !< creep/relax time ratio based on tau-method
    real(SP), allocatable, public :: ts(:)                                     !< relaxation time of visco-elastic medium

    logical, public :: benchmark_mode                                          !< true for fixed parameter run
    logical, public :: pw_mode                                                 !< Plane wave mode
    logical, public :: bf_mode                                                 !< Body force soruce mode
    logical, public :: green_mode                                              !< Green's function computaiton with reciprocity
    !logical, public :: fullspace_mode

    character(80), public :: title                                             !< execution title, used in filename and headers
    integer, public :: exedate                                                 !< date and time by seconds from 1970/1/1 0:0:0

    integer, public :: nx, ny, nz                                              !< space grid number (global)
    integer, public :: nt                                                      !< time grid number
    real(MP), public :: dx, dy, dz                                             !< space grid width
    real(SP), public :: dt                                                     !< time  grid width
    real(SP), public :: xbeg, xend                                             !< global coordinate: x start / end
    real(SP), public :: ybeg, yend                                             !< global coordinate: y start / end
    real(SP), public :: zbeg, zend                                             !< global coordinate: z start / end
    real(SP), public :: tbeg, tend                                             !< beggining and ending elapsed time

    real(SP), public :: vmin                                                   !< minimum velocity
    real(SP), public :: vmax                                                   !< maximum velocity
    real(SP), public :: fmax                                                   !< maximum frequency by the source
    real(SP), public :: fcut                                                   !< cut-off frequency by the source

    integer, public :: nproc_x                                                 !< process numbers for x/i - direction
    integer, public :: nproc_y                                                 !< process numbers for y/j - direction
    integer, public :: nproc                                                   !< total   numbers of process
    integer, public :: nxp                                                     !< space grid number in the assigned node
    integer, public :: nyp                                                     !< space grid number in the assigned node
    integer, public :: myid                                                    !< MPI node number
    integer, public :: idx, idy                                                !< 2D horizontal division ID
    integer, allocatable, public :: itbl(:,:)                                  !< node layout table
    integer, public :: ibeg, iend                                              !< i-region in the node
    integer, public :: jbeg, jend                                              !< j-region in the node
    integer, public :: kbeg, kend                                              !< k-region in the node
    integer, public :: ibeg_m, iend_m                                          !< i- memory allocation area
    integer, public :: jbeg_m, jend_m                                          !< j- memory allocation area
    integer, public :: kbeg_m, kend_m                                          !< k- memory allocation area
    integer, public :: ipad, jpad, kpad                                        !< memory padding size for optimization

    integer, public :: na                                                      !< absorber thickness
    integer, public :: ibeg_k, iend_k                                          !< i- kernel integration area w/o absorption band
    integer, public :: jbeg_k, jend_k                                          !< j- kernel integration area w/o absorption band
    integer, public :: kbeg_k, kend_k                                          !< k- kernel integration area w/o absorption band
    integer, allocatable, public :: kbeg_a(:,:)                                !< k>=kbeg_a(i,j) is in absorber region
    character(16), public :: abc_type

    real(SP), public :: M0                                                     !< total moment

    character(256), public :: odir                                             !< output directory

    integer, allocatable, public :: kfs(:,:)                                   !< free surface depth grid in the node
    integer, allocatable, public :: kob(:,:)                                   !< ocean bottom depth grid in the node
    integer, allocatable, public :: kfs_top(:,:), kfs_bot(:,:)                 !< region of 2nd-order FDM for free surface
    integer, allocatable, public :: kob_top(:,:), kob_bot(:,:)                 !< region of 2nd-order FDM for ocean bottom
    real(SP), allocatable, public :: bddep(:,:,:)                              !< boundary depth in physical coordinate

    real(SP), public :: clon                                                   !< center longitude
    real(SP), public :: clat                                                   !< center latitude
    real(SP), public :: phi                                                    !< azimuth
    real(SP), allocatable, public :: xc(:), yc(:), zc(:)

    real(SP), public :: evlo
    real(SP), public :: evla
    real(SP), public :: evdp                                                   !< unit:km
    real(SP), public :: mxx0, myy0, mzz0, myz0, mxz0, mxy0
    real(SP), public :: fx0, fy0, fz0
    real(SP), public :: otim
    real(SP), public :: sx0, sy0

    real(MP), private, allocatable :: sbuf_ip(:), sbuf_im(:)                   !< mpi send buffer for x-dir
    real(MP), private, allocatable :: sbuf_jp(:), sbuf_jm(:)                   !< mpi send buffer for y-dir
    real(MP), private, allocatable :: rbuf_ip(:), rbuf_im(:)                   !< mpi recv buffer for x-dir
    real(MP), private, allocatable :: rbuf_jp(:), rbuf_jm(:)                   !< mpi recv buffer for y-dir

    integer, private :: mpi_precision

    private :: inside_node
    private :: set_mpi_table

contains

    subroutine global__readprm(io_prm)

        integer, intent(in) :: io_prm

        call readini(io_prm, 'benchmark_mode', benchmark_mode, .false.)

        ! read parameters
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

        ! some parameters are fixed for benchmark mode
        if (benchmark_mode) then
            dx = 0.5
            dy = 0.5
            dz = 0.5
            dt = 0.04
            na = 20
            xbeg = -nx/2.0*real(dx) ! force x=0 at center
            ybeg = -ny/2.0*real(dy) ! force y=0 at center
            zbeg = -30*real(dz)
            tbeg = 0.0
            clon = 139.7604
            clat = 35.7182
            phi = 0.0
            abc_type = 'pml'
!      fullspace_mode = .false.
        else ! or read from file for regular run
            call readini(io_prm, 'dx', dx, 0.5_mp)
            call readini(io_prm, 'dy', dy, 0.5_mp)
            call readini(io_prm, 'dz', dz, 0.5_mp)
            call readini(io_prm, 'dt', dt, 0.01)
            call readini(io_prm, 'na', na, 20)
            call readini(io_prm, 'xbeg', xbeg, -nx/2*real(dx))
            call readini(io_prm, 'ybeg', ybeg, -ny/2*real(dy))
            call readini(io_prm, 'zbeg', zbeg, -30*real(dz))
            call readini(io_prm, 'tbeg', tbeg, 0.0)
            call readini(io_prm, 'clon', clon, 139.7604)
            call readini(io_prm, 'clat', clat, 35.7182)
            call readini(io_prm, 'phi', phi, 0.0)
            call readini(io_prm, 'abc_type', abc_type, 'pml')
!      call readini( io_prm, 'fullspace_mode', fullspace_mode, .false. )

        end if

    end subroutine global__readprm

    !read parameter file, memory allocation, MPI set-ups
    subroutine global__setup(io_prm)

        integer, intent(in) :: io_prm
        integer :: err

        ! MPI status check
        call mpi_comm_rank(mpi_comm_world, myid, err)
        if (MP == DP) then
            mpi_precision = MPI_DOUBLE_PRECISION
        else
            mpi_precision = MPI_REAL
        end if

        ! read key parameters
        call global__readprm(io_prm)

        nproc = nproc_x*nproc_y ! total number of processes

        ! obtain date by unixtime: seconds measured from 1970/1/1 0:0:0
        if (myid == 0) call daytim__getdate(exedate)
        call mpi_bcast(exedate, 1, MPI_INTEGER, 0, mpi_comm_world, err)

        ! derived parameters
        xend = xbeg+nx*real(dx)
        yend = ybeg+ny*real(dy)
        zend = zbeg+nz*real(dz)
        tend = tbeg+nt*dt

#ifdef _OPENACC
        block
            integer :: ngpus
            ngpus = acc_get_num_devices(acc_device_nvidia)
            call acc_set_device_num(mod(myid, ngpus), acc_device_nvidia)
        end block
#endif                

    end subroutine global__setup

    ! Common parameter setup, needed only for starting
    subroutine global__setup2

        integer :: i, j, k
        integer :: err
        integer :: nproc_exe
        integer :: mx, my
        integer :: proc_x, proc_y
        character(256) :: command

        call pwatch__on("global__setup2") ! measure from here

        ! size settings
        call mpi_comm_size(mpi_comm_world, nproc_exe, err)
        call assert(nproc == nproc_exe)

        mx = mod(nx, nproc_x)
        my = mod(ny, nproc_y)
        proc_x = mod(myid, nproc_x)
        proc_y = myid/nproc_x
        if (proc_x <= nproc_x-mx-1) then
            nxp = (nx-mx)/nproc_x
        else
            nxp = (nx-mx)/nproc_x+1
        end if
        if (proc_y <= nproc_y-my-1) then
            nyp = (ny-my)/nproc_y
        else
            nyp = (ny-my)/nproc_y+1
        end if

        ! MPI coordinate
        allocate (itbl(-1:nproc_x, -1:nproc_y))
        allocate (sbuf_ip(5*nyp*nz), source=0.0_mp)
        allocate (sbuf_im(5*nyp*nz), source=0.0_mp)
        allocate (rbuf_ip(5*nyp*nz), source=0.0_mp)
        allocate (rbuf_im(5*nyp*nz), source=0.0_mp)
        allocate (sbuf_jp(5*nxp*nz), source=0.0_mp)
        allocate (sbuf_jm(5*nxp*nz), source=0.0_mp)
        allocate (rbuf_jp(5*nxp*nz), source=0.0_mp)
        allocate (rbuf_jm(5*nxp*nz), source=0.0_mp)

        ! MPI communication table
        call set_mpi_table

        ! create output directory (if it does not exist)
        call mpi_barrier(mpi_comm_world, err)
        command = 'if [ ! -d '// trim(odir) // ' ]; then mkdir -p ' &
                 // trim(odir) // ' > /dev/null 2>&1 ; fi'
        do i=0, nproc-1
            if (myid == i) then
                call execute_command_line(trim(command))
            end if
            call mpi_barrier(mpi_comm_world, err)
        end do        

        ! computation region in this node (#244)
        if (proc_x <= nproc_x-mx-1) then
            ibeg = proc_x*(nx-mx)/nproc_x+1
            iend = (proc_x+1)*(nx-mx)/nproc_x
        else
            ibeg = proc_x*((nx-mx)/nproc_x+1)-(nproc_x-mx)+1
            iend = (proc_x+1)*((nx-mx)/nproc_x+1)-(nproc_x-mx)
        end if
        if (proc_y <= nproc_y-my-1) then
            jbeg = proc_y*(ny-my)/nproc_y+1
            jend = (proc_y+1)*(ny-my)/nproc_y
        else
            jbeg = proc_y*((ny-my)/nproc_y+1)-(nproc_y-my)+1
            jend = (proc_y+1)*((ny-my)/nproc_y+1)-(nproc_y-my)
        end if

        kbeg = 1
        kend = nz

        ! memory requirements including margin for MPI/boundary conditions
        ! stress glut also requires sleeve area
        ibeg_m = ibeg-3
        iend_m = iend+3+ipad
        jbeg_m = jbeg-3
        jend_m = jend+3+jpad
        kbeg_m = kbeg-3
        kend_m = kend+3+kpad

        ! coordinate setting
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

        ! absorbing boundary region                     -+---> i,j(x,y)
        !                                                 |
        !  +-----+--------------------------+-----+      |
        !  |     |                          |     |      v k(z)
        !  |     |                          |     |
        !  |     |                          |     |
        !  |     |                          |     |
        !  |     |  interior region         |     |
        !  |     |  eveluated by m_kernel   |     |
        !  |     |                          |     |
        !  |     |                          |     |
        !  |     |                          |     |
        !  |     +--------------------------+     |
        !  |         exterior region              |
        !  |         evaluated by m_absorb        |
        !  +-----+--------------------------+-----+
        !  1      na                        nx-na+1  nx
        !  <- na ->

        allocate (kbeg_a(ibeg_m:iend_m, jbeg_m:jend_m))
        do j = jbeg_m, jend_m
            do i = ibeg_m, iend_m
                if (i <= na .or. nx-na+1 <= i .or. j <= na .or. ny-na+1 <= j) then
                    kbeg_a(i, j) = kbeg
                else
                    kbeg_a(i, j) = kend-na+1
                end if
            end do
        end do

        ! Interior Kernel region

        ! initial value
        ibeg_k = ibeg
        iend_k = iend
        jbeg_k = jbeg
        jend_k = jend
        kbeg_k = kbeg
        kend_k = kend

        if (abc_type == 'pml') then

!      if( fullspace_mode ) kbeg_k = na + 1

            if (iend <= na) then; ibeg_k = iend+1; ! no kernel integration
            else if (ibeg <= na) then; ibeg_k = na+1; ! pertial kernel
            end if

            if (ibeg >= nx-na+1) then; iend_k = ibeg-1; ! no kernel integartion
            else if (iend >= nx-na+1) then; iend_k = nx-na; 
            end if

            if (jend <= na) then; jbeg_k = jend+1; ! no kernel integration
            else if (jbeg <= na) then; jbeg_k = na+1; ! pertial kernel
            end if

            if (jbeg >= ny-na+1) then; jend_k = jbeg-1; ! no kernel integartion
            else if (jend >= ny-na+1) then; jend_k = ny-na; 
            end if
            kend_k = nz-na

        end if

        !$acc enter data copyin(&
        !$acc       sbuf_ip, sbuf_im, &
        !$acc       rbuf_ip, rbuf_im, &
        !$acc       sbuf_jp, sbuf_jm, &
        !$acc       rbuf_jp, rbuf_jm, &
        !$acc       itbl(-1:nproc_x, -1:nproc_y))
 

        call pwatch__off("global__setup2") ! measure from here

    end subroutine global__setup2

    ! Data buffring & communication for velocity vector
    subroutine global__comm_vel()

        integer :: isize, jsize
        integer :: err
        integer :: istatus(mpi_status_size, 4)
        integer :: req_i(4), req_j(4)
        integer :: i, j

        if (myid >= nproc) return

        call pwatch__on("global__comm_vel")

        ! unit buffer size
        isize = nyp*nz
        jsize = nxp*nz

        !$acc host_data use_device(rbuf_ip, rbuf_im, rbuf_jp, rbuf_jm)
        call mpi_irecv(rbuf_ip, 5*isize, mpi_precision, itbl(idx+1, idy), 1, mpi_comm_world, req_i(1), err)
        call mpi_irecv(rbuf_im, 4*isize, mpi_precision, itbl(idx-1, idy), 2, mpi_comm_world, req_i(2), err)
        call mpi_irecv(rbuf_jp, 5*jsize, mpi_precision, itbl(idx, idy+1), 3, mpi_comm_world, req_j(1), err)
        call mpi_irecv(rbuf_jm, 4*jsize, mpi_precision, itbl(idx, idy-1), 4, mpi_comm_world, req_j(2), err)
        !$acc end host_data

        !$acc kernels present(Vx, Vy, Vz, sbuf_ip, sbuf_im) async(1)
        !$acc loop independent
        do j=jbeg, jend
            sbuf_ip(0*isize+(j-jbeg)*nz+1:0*isize+(j-jbeg+1)*nz) = Vx(1:nz,iend-1,j)
            sbuf_ip(1*isize+(j-jbeg)*nz+1:1*isize+(j-jbeg+1)*nz) = Vx(1:nz,iend  ,j)
            sbuf_ip(2*isize+(j-jbeg)*nz+1:2*isize+(j-jbeg+1)*nz) = Vy(1:nz,iend  ,j)
            sbuf_ip(3*isize+(j-jbeg)*nz+1:3*isize+(j-jbeg+1)*nz) = Vz(1:nz,iend  ,j)
            
            sbuf_im(0*isize+(j-jbeg)*nz+1:0*isize+(j-jbeg+1)*nz) = Vx(1:nz,ibeg  ,j)
            sbuf_im(1*isize+(j-jbeg)*nz+1:1*isize+(j-jbeg+1)*nz) = Vy(1:nz,ibeg  ,j)
            sbuf_im(2*isize+(j-jbeg)*nz+1:2*isize+(j-jbeg+1)*nz) = Vy(1:nz,ibeg+1,j)
            sbuf_im(3*isize+(j-jbeg)*nz+1:3*isize+(j-jbeg+1)*nz) = Vz(1:nz,ibeg  ,j)
            sbuf_im(4*isize+(j-jbeg)*nz+1:4*isize+(j-jbeg+1)*nz) = Vz(1:nz,ibeg+1,j)
        end do
        !$acc end kernels

        !$acc kernels present(Vx, Vy, Vz, sbuf_jp, sbuf_jm) async(2)
        !$acc loop independent
        do i=ibeg, iend
            sbuf_jp(0*jsize+(i-ibeg)*nz+1:0*jsize+(i-ibeg+1)*nz) = Vx(1:nz,i,jend  )
            sbuf_jp(1*jsize+(i-ibeg)*nz+1:1*jsize+(i-ibeg+1)*nz) = Vy(1:nz,i,jend-1)
            sbuf_jp(2*jsize+(i-ibeg)*nz+1:2*jsize+(i-ibeg+1)*nz) = Vy(1:nz,i,jend  )
            sbuf_jp(3*jsize+(i-ibeg)*nz+1:3*jsize+(i-ibeg+1)*nz) = Vz(1:nz,i,jend  )

            sbuf_jm(0*jsize+(i-ibeg)*nz+1:0*jsize+(i-ibeg+1)*nz) = Vx(1:nz,i,jbeg  )
            sbuf_jm(1*jsize+(i-ibeg)*nz+1:1*jsize+(i-ibeg+1)*nz) = Vx(1:nz,i,jbeg+1)
            sbuf_jm(2*jsize+(i-ibeg)*nz+1:2*jsize+(i-ibeg+1)*nz) = Vy(1:nz,i,jbeg  )
            sbuf_jm(3*jsize+(i-ibeg)*nz+1:3*jsize+(i-ibeg+1)*nz) = Vz(1:nz,i,jbeg  )
            sbuf_jm(4*jsize+(i-ibeg)*nz+1:4*jsize+(i-ibeg+1)*nz) = Vz(1:nz,i,jbeg+1)
        end do
        !$acc end kernels

        !$acc wait

        !$acc host_data use_device(sbuf_ip, sbuf_im, sbuf_jp, sbuf_jm)
        call mpi_isend(sbuf_ip, 4*isize, mpi_precision, itbl(idx+1, idy), 2, mpi_comm_world, req_i(3), err)
        call mpi_isend(sbuf_im, 5*isize, mpi_precision, itbl(idx-1, idy), 1, mpi_comm_world, req_i(4), err)
        call mpi_isend(sbuf_jp, 4*jsize, mpi_precision, itbl(idx, idy+1), 4, mpi_comm_world, req_j(3), err)
        call mpi_isend(sbuf_jm, 5*jsize, mpi_precision, itbl(idx, idy-1), 3, mpi_comm_world, req_j(4), err)
        !$acc end host_data

        call mpi_waitall(4, req_i, istatus, err)
        call mpi_waitall(4, req_j, istatus, err)

        !$acc kernels present(Vx, Vy, Vz, rbuf_ip, rbuf_im) async(1)
        !$acc loop independent
        do j=jbeg, jend
            Vx(1:nz,ibeg-2,j) = rbuf_im(0*isize+(j-jbeg)*nz+1:0*isize+(j-jbeg+1)*nz)
            Vx(1:nz,ibeg-1,j) = rbuf_im(1*isize+(j-jbeg)*nz+1:1*isize+(j-jbeg+1)*nz)
            Vy(1:nz,ibeg-1,j) = rbuf_im(2*isize+(j-jbeg)*nz+1:2*isize+(j-jbeg+1)*nz)
            Vz(1:nz,ibeg-1,j) = rbuf_im(3*isize+(j-jbeg)*nz+1:3*isize+(j-jbeg+1)*nz)

            Vx(1:nz,iend+1,j) = rbuf_ip(0*isize+(j-jbeg)*nz+1:0*isize+(j-jbeg+1)*nz)
            Vy(1:nz,iend+1,j) = rbuf_ip(1*isize+(j-jbeg)*nz+1:1*isize+(j-jbeg+1)*nz)
            Vy(1:nz,iend+2,j) = rbuf_ip(2*isize+(j-jbeg)*nz+1:2*isize+(j-jbeg+1)*nz)
            Vz(1:nz,iend+1,j) = rbuf_ip(3*isize+(j-jbeg)*nz+1:3*isize+(j-jbeg+1)*nz)
            Vz(1:nz,iend+2,j) = rbuf_ip(4*isize+(j-jbeg)*nz+1:4*isize+(j-jbeg+1)*nz)
        end do
        !$acc end kernels

        !$acc kernels present(Vx, Vy, Vz, rbuf_jp, rbuf_jm) async(2)
        !$acc loop independent
        do i=ibeg, iend
            Vx(1:nz,i,jbeg-1) = rbuf_jm(0*jsize+(i-ibeg)*nz+1:0*jsize+(i-ibeg+1)*nz)
            Vy(1:nz,i,jbeg-2) = rbuf_jm(1*jsize+(i-ibeg)*nz+1:1*jsize+(i-ibeg+1)*nz)
            Vy(1:nz,i,jbeg-1) = rbuf_jm(2*jsize+(i-ibeg)*nz+1:2*jsize+(i-ibeg+1)*nz)
            Vz(1:nz,i,jbeg-1) = rbuf_jm(3*jsize+(i-ibeg)*nz+1:3*jsize+(i-ibeg+1)*nz)

            Vx(1:nz,i,jend+1) = rbuf_jp(0*jsize+(i-ibeg)*nz+1:0*jsize+(i-ibeg+1)*nz)
            Vx(1:nz,i,jend+2) = rbuf_jp(1*jsize+(i-ibeg)*nz+1:1*jsize+(i-ibeg+1)*nz)
            Vy(1:nz,i,jend+1) = rbuf_jp(2*jsize+(i-ibeg)*nz+1:2*jsize+(i-ibeg+1)*nz)
            Vz(1:nz,i,jend+1) = rbuf_jp(3*jsize+(i-ibeg)*nz+1:3*jsize+(i-ibeg+1)*nz)
            Vz(1:nz,i,jend+2) = rbuf_jp(4*jsize+(i-ibeg)*nz+1:4*jsize+(i-ibeg+1)*nz)
        end do
        !$acc end kernels

        !$acc wait

        call pwatch__off("global__comm_vel")

    end subroutine global__comm_vel

    !>
    ! Data buffring & communication for stress tensor
    !<
  !!
    subroutine global__comm_stress()

        integer :: isize, jsize
        integer :: err
        integer :: istatus(mpi_status_size, 4)
        integer :: req_i(4), req_j(4)
        integer :: i, j

        if (myid >= nproc) return

        call pwatch__on("global__comm_stress")

        ! unit buffer size
        isize = nyp*nz
        jsize = nxp*nz

        !$acc host_data use_device(rbuf_ip, rbuf_jp, rbuf_im, rbuf_jm)
        call mpi_irecv(rbuf_ip, 4*isize, mpi_precision, itbl(idx+1, idy), 5, mpi_comm_world, req_i(1), err)
        call mpi_irecv(rbuf_im, 5*isize, mpi_precision, itbl(idx-1, idy), 6, mpi_comm_world, req_i(2), err)
        call mpi_irecv(rbuf_jp, 4*jsize, mpi_precision, itbl(idx, idy+1), 7, mpi_comm_world, req_j(1), err)
        call mpi_irecv(rbuf_jm, 5*jsize, mpi_precision, itbl(idx, idy-1), 8, mpi_comm_world, req_j(2), err)
        !$acc end host_data

        !$acc kernels async(1) &
        !$acc present(Sxx, Sxy, Sxz, sbuf_ip, sbuf_im) 
        !$acc loop independent 
        do j=jbeg, jend
            sbuf_ip(0*isize+(j-jbeg)*nz+1:0*isize+(j-jbeg+1)*nz) = Sxx(1:nz,iend  ,j)
            sbuf_ip(1*isize+(j-jbeg)*nz+1:1*isize+(j-jbeg+1)*nz) = Sxy(1:nz,iend-1,j)
            sbuf_ip(2*isize+(j-jbeg)*nz+1:2*isize+(j-jbeg+1)*nz) = Sxy(1:nz,iend  ,j)
            sbuf_ip(3*isize+(j-jbeg)*nz+1:3*isize+(j-jbeg+1)*nz) = Sxz(1:nz,iend-1,j)
            sbuf_ip(4*isize+(j-jbeg)*nz+1:4*isize+(j-jbeg+1)*nz) = Sxz(1:nz,iend  ,j)

            sbuf_im(0*isize+(j-jbeg)*nz+1:0*isize+(j-jbeg+1)*nz) = Sxx(1:nz,ibeg  ,j)
            sbuf_im(1*isize+(j-jbeg)*nz+1:1*isize+(j-jbeg+1)*nz) = Sxx(1:nz,ibeg+1,j)
            sbuf_im(2*isize+(j-jbeg)*nz+1:2*isize+(j-jbeg+1)*nz) = Sxy(1:nz,ibeg  ,j)
            sbuf_im(3*isize+(j-jbeg)*nz+1:3*isize+(j-jbeg+1)*nz) = Sxz(1:nz,ibeg  ,j)
        end do
        !$acc end kernels

        !$acc kernels async(2) &
        !$acc present(Syy, Sxy, Syz, sbuf_jp, sbuf_jm) 
        !$acc loop independent
        do i=ibeg, iend
            sbuf_jp(0*jsize+(i-ibeg)*nz+1:0*jsize+(i-ibeg+1)*nz) = Syy(1:nz,i,jend  )
            sbuf_jp(1*jsize+(i-ibeg)*nz+1:1*jsize+(i-ibeg+1)*nz) = Sxy(1:nz,i,jend-1)
            sbuf_jp(2*jsize+(i-ibeg)*nz+1:2*jsize+(i-ibeg+1)*nz) = Sxy(1:nz,i,jend  )
            sbuf_jp(3*jsize+(i-ibeg)*nz+1:3*jsize+(i-ibeg+1)*nz) = Syz(1:nz,i,jend-1)
            sbuf_jp(4*jsize+(i-ibeg)*nz+1:4*jsize+(i-ibeg+1)*nz) = Syz(1:nz,i,jend  )

            sbuf_jm(0*jsize+(i-ibeg)*nz+1:0*jsize+(i-ibeg+1)*nz) = Syy(1:nz,i,jbeg  )
            sbuf_jm(1*jsize+(i-ibeg)*nz+1:1*jsize+(i-ibeg+1)*nz) = Syy(1:nz,i,jbeg+1)
            sbuf_jm(2*jsize+(i-ibeg)*nz+1:2*jsize+(i-ibeg+1)*nz) = Sxy(1:nz,i,jbeg  )
            sbuf_jm(3*jsize+(i-ibeg)*nz+1:3*jsize+(i-ibeg+1)*nz) = Syz(1:nz,i,jbeg  )
        end do
        !$acc end kernels

        !$acc wait

        !$acc host_data use_device(sbuf_ip, sbuf_jp, sbuf_im, sbuf_jm)
        call mpi_isend(sbuf_ip, 5*isize, mpi_precision, itbl(idx+1, idy), 6, mpi_comm_world, req_i(3), err)
        call mpi_isend(sbuf_im, 4*isize, mpi_precision, itbl(idx-1, idy), 5, mpi_comm_world, req_i(4), err)
        call mpi_isend(sbuf_jp, 5*jsize, mpi_precision, itbl(idx, idy+1), 8, mpi_comm_world, req_j(3), err)
        call mpi_isend(sbuf_jm, 4*jsize, mpi_precision, itbl(idx, idy-1), 7, mpi_comm_world, req_j(4), err)
        !$acc end host_data

        call mpi_waitall(4, req_i, istatus, err)
        call mpi_waitall(4, req_j, istatus, err)

        !$acc kernels async(1) &
        !$acc present(Sxx, Sxy, Sxz, rbuf_ip, rbuf_im) 
        !$acc loop independent
        do j=jbeg, jend
            Sxx(1:nz,ibeg-1,j) = rbuf_im(0*isize+(j-jbeg)*nz+1:0*isize+(j-jbeg+1)*nz)
            Sxy(1:nz,ibeg-2,j) = rbuf_im(1*isize+(j-jbeg)*nz+1:1*isize+(j-jbeg+1)*nz)
            Sxy(1:nz,ibeg-1,j) = rbuf_im(2*isize+(j-jbeg)*nz+1:2*isize+(j-jbeg+1)*nz)
            Sxz(1:nz,ibeg-2,j) = rbuf_im(3*isize+(j-jbeg)*nz+1:3*isize+(j-jbeg+1)*nz)
            Sxz(1:nz,ibeg-1,j) = rbuf_im(4*isize+(j-jbeg)*nz+1:4*isize+(j-jbeg+1)*nz)

            Sxx(1:nz,iend+1,j) = rbuf_ip(0*isize+(j-jbeg)*nz+1:0*isize+(j-jbeg+1)*nz)
            Sxx(1:nz,iend+2,j) = rbuf_ip(1*isize+(j-jbeg)*nz+1:1*isize+(j-jbeg+1)*nz)
            Sxy(1:nz,iend+1,j) = rbuf_ip(2*isize+(j-jbeg)*nz+1:2*isize+(j-jbeg+1)*nz)
            Sxz(1:nz,iend+1,j) = rbuf_ip(3*isize+(j-jbeg)*nz+1:3*isize+(j-jbeg+1)*nz)
        end do
        !$acc end kernels

        !$acc kernels async(2) &
        !$acc present(Syy, Sxy, Syz, rbuf_jp, rbuf_jm) 
        !$acc loop independent
        do i=ibeg, iend
            Syy(1:nz,i,jbeg-1) = rbuf_jm(0*jsize+(i-ibeg)*nz+1:0*jsize+(i-ibeg+1)*nz)
            Sxy(1:nz,i,jbeg-2) = rbuf_jm(1*jsize+(i-ibeg)*nz+1:1*jsize+(i-ibeg+1)*nz)
            Sxy(1:nz,i,jbeg-1) = rbuf_jm(2*jsize+(i-ibeg)*nz+1:2*jsize+(i-ibeg+1)*nz)
            Syz(1:nz,i,jbeg-2) = rbuf_jm(3*jsize+(i-ibeg)*nz+1:3*jsize+(i-ibeg+1)*nz)
            Syz(1:nz,i,jbeg-1) = rbuf_jm(4*jsize+(i-ibeg)*nz+1:4*jsize+(i-ibeg+1)*nz)

            Syy(1:nz,i,jend+1) = rbuf_jp(0*jsize+(i-ibeg)*nz+1:0*jsize+(i-ibeg+1)*nz)
            Syy(1:nz,i,jend+2) = rbuf_jp(1*jsize+(i-ibeg)*nz+1:1*jsize+(i-ibeg+1)*nz)
            Sxy(1:nz,i,jend+1) = rbuf_jp(2*jsize+(i-ibeg)*nz+1:2*jsize+(i-ibeg+1)*nz)
            Syz(1:nz,i,jend+1) = rbuf_jp(3*jsize+(i-ibeg)*nz+1:3*jsize+(i-ibeg+1)*nz)
        end do
        !$acc end kernels

        !$acc wait

        call pwatch__off("global__comm_stress")

    end subroutine global__comm_stress

    ! check if the voxcel location is inside the MPI node
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

        ! 2D communication table
        !> @see 2013-0439
        itbl(-1:nproc_x, -1:nproc_y) = MPI_PROC_NULL
        do i = 0, nproc-1

            ii = mod(i, nproc_x)
            jj = i/nproc_x

            itbl(ii, jj) = i
        end do

        ! location of this process
        idx = mod(myid, nproc_x)
        idy = myid/nproc_x
    end subroutine set_mpi_table

    subroutine global__getnode(i, j, idx_ij, idy_ij)

        integer, intent(in) :: i, j
        integer, intent(out) :: idx_ij, idy_ij
        integer :: mx, my
        integer :: iproc_x, ibeg_node, iend_node
        integer :: iproc_y, jbeg_node, jend_node

        mx = mod(nx, nproc_x)
        my = mod(ny, nproc_y)

        do iproc_x = 0, nproc_x-1
            if (iproc_x <= nproc_x-mx-1) then
                ibeg_node = iproc_x*(nx-mx)/nproc_x+1
                iend_node = (iproc_x+1)*(nx-mx)/nproc_x
                if (ibeg_node <= i .and. i <= iend_node) then
                    idx_ij = iproc_x
                    exit
                end if
            end if
        end do

        do iproc_y = 0, nproc_y-1
            if (iproc_y <= nproc_y-my-1) then
                jbeg_node = iproc_y*(ny-my)/nproc_y+1
                jend_node = (iproc_y+1)*(ny-my)/nproc_y
                if (jbeg_node <= j .and. j <= jend_node) then
                    idy_ij = iproc_y
                    exit
                end if
            end if
        end do
    end subroutine global__getnode
end module m_global
! ------------------------------------------------------------------------------------------------------------------------------ !!
