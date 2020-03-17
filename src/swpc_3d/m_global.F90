!! ----------------------------------------------------------------------------------------------------------------------------- !!
!>
!! global control parameters, shared arrays and MPI communication
!!
!! @copyright
!!   Copyright 2013-2020 Takuto Maeda. All rights reserved. This project is released under the MIT license.
!<
!! ----
#include "m_debug.h"
module m_global

  !! modules
  use m_std
  use m_debug
  use m_fdtool
  use m_pwatch
  use m_system
  use m_daytim
  use m_readini
  use mpi

  !! declarations
  implicit none
  public
  save

  !!
  !! public routine name
  !!
  public :: global__setup
  public :: global__setup2
  public :: global__comm_vel
  public :: global__comm_stress



  !!
  !! fixed parameters
  !!
  integer,  parameter :: Nl = 4                                     !< FDM order ( 4th )
  integer,  parameter :: Nsl = Nl/2                                 !< thickness of "sleeve area"
  real(SP)            :: UC = 10.0**(-15)                           !< Conventional -> SI unit for moment tensor of [Nm]
  integer,  parameter :: MP = DP                                    !< DP for mixed precision, SP for pure single precision
  integer,  parameter :: NM = 3                                     !< Number of memory variables
  integer,  parameter :: NBD = 9                                    !< Number of boundary depths to be memorized
  !!
  !! global arrays
  !!
  real(MP), allocatable :: Vx(:,:,:),    Vy(:,:,:),    Vz(:,:,:)    !<  velocity components
  real(MP), allocatable :: Sxx(:,:,:),   Syy(:,:,:),   Szz(:,:,:)   !<  normal stress components
  real(MP), allocatable :: Syz(:,:,:),   Sxz(:,:,:),   Sxy(:,:,:)   !<  shear  stress components
  real(SP), allocatable :: Rxx(:,:,:,:), Ryy(:,:,:,:), Rzz(:,:,:,:) !<  memory variables: normal components
  real(SP), allocatable :: Ryz(:,:,:,:), Rxz(:,:,:,:), Rxy(:,:,:,:) !<  memory variables: shear  components
  real(SP), allocatable :: rho(:,:,:),   lam(:,:,:),   mu(:,:,:)    !<  density and relaxed moduli
  real(SP), allocatable :: bx(:,:,:),    by(:,:,:),    bz(:,:,:)    !<  inverse density (buoyancy) at velocity grids
  real(SP), allocatable :: muyz(:,:,:),  muxz(:,:,:),  muxy(:,:,:)  !<  averaged rigidity at shear stress compoments
  real(SP), allocatable :: taup(:,:,:),  taus(:,:,:)                !<  creep/relax time ratio based on tau-method for atten.
  real(SP), allocatable :: ts(:)                                    !<  relaxation time of visco-elastic medium

  !!
  !! mode
  !!
  logical               :: benchmark_mode                           !<  true for fixed parameter run

  !!
  !! title, date
  !!
  character(80)         :: title                                    !<  execution title, used in filename and headers
  integer               :: exedate                                  !<  date and time by seconds from 1970/1/1 0:0:0

  !!
  !! global control parameters
  !!
  integer               :: nx, ny, nz                               !<  space grid number (global)
  integer               :: nt                                       !<  time grid number
  real(MP)              :: dx, dy, dz                               !<  space grid width
  real(SP)              :: dt                                       !<  time  grid width
  real(SP)              :: xbeg, xend                               !<  global coordinate: x start / end
  real(SP)              :: ybeg, yend                               !<  global coordinate: y start / end
  real(SP)              :: zbeg, zend                               !<  global coordinate: z start / end
  real(SP)              :: tbeg, tend                               !<  beggining and ending elapsed time

  !!
  !! Medium Info
  !!
  real(SP)              :: vmin                                     !<  minimum velocity
  real(SP)              :: vmax                                     !<  maximum velocity
  real(SP)              :: fmax                                     !<  maximum frequency by the source
  real(SP)              :: fcut                                     !<  cut-off frequency by the source

  !!
  !! MPI domain
  !!
  integer               :: nproc_x                                  !<  process numbers for x/i - direction
  integer               :: nproc_y                                  !<  process numbers for y/j - direction
  integer               :: nproc                                    !<  total   numbers of process
  integer               :: nxp                                      !<  space grid number in the assigned node
  integer               :: nyp                                      !<  space grid number in the assigned node
  integer               :: myid                                     !<  MPI node number
  integer               :: idx, idy                                 !<  2D horizontal division ID
  integer, allocatable  :: itbl(:,:)                                !<  node layout table
  integer               :: ibeg, iend                               !<  i-region in the node
  integer               :: jbeg, jend                               !<  j-region in the node
  integer               :: kbeg, kend                               !<  k-region in the node
  integer               :: ibeg_m, iend_m                           !<  i- memory allocation area
  integer               :: jbeg_m, jend_m                           !<  j- memory allocation area
  integer               :: kbeg_m, kend_m                           !<  k- memory allocation area
  integer               :: ipad, jpad, kpad                         !<  memory padding size for optimization

  !!
  !! Absorbing boundary
  !!
  integer               :: na                                       !<  absorber thickness
  integer               :: ibeg_k, iend_k                           !<  i- kernel integration area without absorption band
  integer               :: jbeg_k, jend_k                           !<  j- kernel integration area without absorption band
  integer               :: kbeg_k, kend_k                           !<  k- kernel integration area without absorption band
  integer, allocatable  :: kbeg_a(:,:)                              !<  k>=kbeg_a(i,j) is in absorber region
  character(16)         :: abc_type

  !!
  !! Source
  !!
  real(SP)              :: M0                                       !<  total moment



  !!
  !! output directory
  !!
  character(256) :: odir

  !!
  !! free surface / ocean bottom boundary
  !!
  integer,  allocatable :: kfs(:,:)                                 !<  free surface depth grid in the node
  integer,  allocatable :: kob(:,:)                                 !<  ocean bottom depth grid in the node
  integer,  allocatable :: kfs_top(:,:), kfs_bot(:,:)               !<  region in which 2nd-order FDM is applied for free surface
  integer,  allocatable :: kob_top(:,:), kob_bot(:,:)               !<  region in which 2nd-order FDM is applied for ocean bottom
  real(SP), allocatable :: bddep(:,:,:)                             !<  boundary depth in physical coordinate

  !!
  !! map coordinate: use dummy value for cartesian problem
  !!
  real(SP) :: clon                                                  !< center longitude
  real(SP) :: clat                                                  !< center latitude
  real(SP) :: phi                                                   !< azimuth
  real(SP), allocatable :: xc(:), yc(:), zc(:)

  !! Source: remember first source grid location in lat/lon form
  real(SP) :: evlo
  real(SP) :: evla
  real(SP) :: evdp !< unit:km
  real(SP) :: mxx0, myy0, mzz0, myz0, mxz0, mxy0
  real(SP) :: fx0, fy0, fz0
  real(SP) :: otim
  real(SP) :: sx0, sy0

  !! Special modes
  logical  :: pw_mode                                               !< Plane wave mode
  logical  :: bf_mode                                               !< Body force soruce mode
  logical  :: green_mode                                            !< Green's function computaiton with reciprocity

  !!
  !! private variables
  !!
  real(MP), private, allocatable :: sbuf_ip(:), sbuf_im(:)          !<  mpi send buffer for x-dir
  real(MP), private, allocatable :: sbuf_jp(:), sbuf_jm(:)          !<  mpi send buffer for y-dir
  real(MP), private, allocatable :: rbuf_ip(:), rbuf_im(:)          !<  mpi recv buffer for x-dir
  real(MP), private, allocatable :: rbuf_jp(:), rbuf_jm(:)          !<  mpi recv buffer for y-dir
  integer :: mpi_precision


  !! ----
  private :: inside_node
  private :: set_mpi_table


contains

  !! ---------------------------------------------------------------------------------------------------------------------------- !!
  subroutine global__readprm( io_prm )

    integer, intent(in) :: io_prm

    call readini( io_prm, 'benchmark_mode', benchmark_mode, .false.  )

    !! read parameters
    call readini( io_prm, 'title',          title,          'swpc3d'       )
    call readini( io_prm, 'nproc_x',        nproc_x,         1              )
    call readini( io_prm, 'nproc_y',        nproc_y,         2              )
    call readini( io_prm, 'nx',             nx,              256            )
    call readini( io_prm, 'ny',             ny,              256            )
    call readini( io_prm, 'nz',             nz,              256            )
    call readini( io_prm, 'nt',             nt,              1000           )
    call readini( io_prm, 'ipad',           ipad,            0              )
    call readini( io_prm, 'jpad',           jpad,            0              )
    call readini( io_prm, 'kpad',           kpad,            0              )
    call readini( io_prm, 'odir',           odir,           './out'         )

    !! some parameters are fixed for benchmark mode
    if( benchmark_mode ) then
      dx = 0.5
      dy = 0.5
      dz = 0.5
      dt = 0.04
      na = 20
      xbeg = -nx/2.0 * dx ! force x=0 at center
      ybeg = -ny/2.0 * dy ! force x=0 at center
      zbeg = - 30 * dz
      tbeg = 0.0
      clon = 139.7604
      clat = 35.7182
      phi  = 0.0
      abc_type = 'pml'
    else !! or read from file for regular run
      call readini( io_prm, 'dx',             dx,              0.5_MP         )
      call readini( io_prm, 'dy',             dy,              0.5_MP         )
      call readini( io_prm, 'dz',             dz,              0.5_MP         )
      call readini( io_prm, 'dt',             dt,              0.01           )
      call readini( io_prm, 'na',             na,              20             )
      call readini( io_prm, 'xbeg',           xbeg,            -nx/2*real(dx) )
      call readini( io_prm, 'ybeg',           ybeg,            -ny/2*real(dy) )
      call readini( io_prm, 'zbeg',           zbeg,            -30*real(dz)   )
      call readini( io_prm, 'tbeg',           tbeg,            0.0            )
      call readini( io_prm, 'clon',           clon,            139.7604       )
      call readini( io_prm, 'clat',           clat,            35.7182        )
      call readini( io_prm, 'phi',            phi,             0.0            )
      call readini( io_prm, 'abc_type',       abc_type,        'pml'          )

    end if


  end subroutine global__readprm
  !! ---------------------------------------------------------------------------------------------------------------------------- !!

  !! ---------------------------------------------------------------------------------------------------------------------------- !!
  !>
  !! read parameter file, memory allocation, MPI set-ups
  !<
  !!
  subroutine global__setup( io_prm )

    integer, intent(in) :: io_prm
    integer :: ierr
    !! ----


    !!
    !! MPI status check
    !!
    call mpi_comm_rank( mpi_comm_world, myid, ierr )
    if( MP == DP ) then
      mpi_precision = MPI_DOUBLE_PRECISION
    else
      mpi_precision = MPI_REAL
    end if


    !!
    !! read key parameters
    !!
    call global__readprm( io_prm )

    nproc = nproc_x * nproc_y            !! total number of processes


    !!
    !! obtain date by unixtime: seconds measured from 1970/1/1 0:0:0
    !!
    if( myid == 0 ) call daytim__getdate( exedate )
    call mpi_bcast( exedate, 1, MPI_INTEGER, 0, mpi_comm_world, ierr )

    !!
    !! derived parameters
    !!
    xend = xbeg + nx * dx
    yend = ybeg + ny * dy
    zend = zbeg + nz * dz
    tend = ybeg + nt * dt



  end subroutine global__setup

  !! --------------------------------------------------------------------------------------------------------------------------- !!

  !! --------------------------------------------------------------------------------------------------------------------------- !!
  !>
  !! Common parameter setup, needed only for starting
  !<
  !! --
  subroutine global__setup2

    integer :: nl3
    integer :: i, j, k
    integer :: ierr
    integer :: nproc_exe

    call pwatch__on( "global__setup2" ) !! measure from here

    !!
    !! output directory create (if it does not exist)
    !!
    call system__call('mkdir -p ' // trim(odir) //'> /dev/null 2>&1' )

    !!
    !! size settings
    !!
    call mpi_comm_size( mpi_comm_world, nproc_exe, ierr )
    call assert( nproc == nproc_exe )

    nxp = ceiling( nx / real(nproc_x) )  !!  nxp-1 <  nx / nproc_x  <= nxp
    nyp = ceiling( ny / real(nproc_y) )  !!  nyp-1 <  ny / nproc_y  <= nyp


    !!
    !! MPI coordinate
    !!

    !! buffer allocation
    !  half of FDM order is used for communication.
    !  Three variables are communicated at once
    nl3 = ( Nl / 2 ) * 3

    allocate( itbl(-1:nproc_x, -1:nproc_y) )
    allocate( sbuf_ip( nyp*nz*nl3 ), sbuf_im( nyp*nz*nl3 ) )
    allocate( rbuf_ip( nyp*nz*nl3 ), rbuf_im( nyp*nz*nl3 ) )
    allocate( sbuf_jp( nxp*nz*nl3 ), sbuf_jm( nxp*nz*nl3 ) )
    allocate( rbuf_jp( nxp*nz*nl3 ), rbuf_jm( nxp*nz*nl3 ) )

    !! initialize buffer
    sbuf_ip(:) = 0.0_MP
    sbuf_im(:) = 0.0_MP
    rbuf_ip(:) = 0.0_MP
    rbuf_im(:) = 0.0_MP
    sbuf_jp(:) = 0.0_MP
    sbuf_jm(:) = 0.0_MP
    rbuf_jp(:) = 0.0_MP
    rbuf_jm(:) = 0.0_MP

    !!
    !! MPI communication table
    !!
    call set_mpi_table


    !!
    !! computation region in this node
    !!
    ibeg = nxp * idx + 1
    iend = min( ibeg + nxp - 1, nx )
    jbeg = nyp * idy + 1
    jend = min( jbeg + nyp - 1, ny )
    kbeg = 1
    kend = nz

    !! re-define node model size
    !! i, j方向の端でサイズが余る場合、モデルサイズを適切に取り直す。
    nxp = iend - ibeg + 1
    nyp = jend - jbeg + 1


    !! memory requirements including margin for MPI/boundary conditions
    !! stress drop also requires sleeve area
    ibeg_m = ibeg - 3
    iend_m = iend + 3 + ipad
    jbeg_m = jbeg - 3
    jend_m = jend + 3 + jpad
    kbeg_m = kbeg - 3
    kend_m = kend + 3 + kpad


    !!
    !! coordinate setting
    !!
    allocate( xc(ibeg_m:iend_m), yc(jbeg_m:jend_m), zc(kbeg_m:kend_m) )

    do i=ibeg_m, iend_m
      xc(i) = i2x( i, xbeg, real(dx) )
    end do
    do j=jbeg_m, jend_m
      yc(j) = j2y( j, ybeg, real(dy) )
    end do
    do k=kbeg_m, kend_m
      zc(k) = k2z( k, zbeg, real(dz) )
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




    allocate( kbeg_a(ibeg_m:iend_m,jbeg_m:jend_m) )
    do j=jbeg_m, jend_m
      do i=ibeg_m, iend_m
        if( i <= na .or. nx-na+1 <= i .or. j <= na .or. ny-na+1 <= j ) then
          kbeg_a(i,j) = kbeg
        else
          kbeg_a(i,j) = kend-na+1
        end if
      end do
    end do

    !!
    !! Interior Kernel region
    !!

    !! initial value
    ibeg_k = ibeg
    iend_k = iend
    jbeg_k = jbeg
    jend_k = jend
    kbeg_k = kbeg
    kend_k = kend

    if( abc_type == 'pml' ) then
      if      ( iend <= na      ) then; ibeg_k = iend+1;  ! no kernel integration
      else if ( ibeg <= na      ) then; ibeg_k = na+1  ;  ! pertial kernel
      end if

      if      ( ibeg >= nx-na+1 ) then; iend_k = ibeg-1;  ! no kernel integartion
      else if ( iend >= nx-na+1 ) then; iend_k = nx-na;
      end if

      if      ( jend <= na      ) then; jbeg_k = jend+1;  ! no kernel integration
      else if ( jbeg <= na      ) then; jbeg_k = na+1;    ! pertial kernel
      end if

      if      ( jbeg >= ny-na+1 ) then; jend_k = jbeg-1;  ! no kernel integartion
      else if ( jend >= ny-na+1 ) then; jend_k = ny-na;
      end if
      kend_k = nz-na
    end if


    call pwatch__off( "global__setup2" ) !! measure from here

  end subroutine global__setup2
  !! ---------------------------------------------------------------------------------------------------------------------------- !!

  !! ---------------------------------------------------------------------------------------------------------------------------- !!
  !>
  !! Data buffring & communication for velocity vector
  !<
  !!
  subroutine global__comm_vel()


    integer :: isize, jsize, s_isize, s_jsize
    integer :: ierr
    integer :: istatus( mpi_status_size, 4 )
    integer :: ireq1(4), ireq2(4)
    integer :: i, j, k, ptr

    if( myid >= nproc ) return

    call pwatch__on( "global__comm_vel" )

    !! unit buffer size
    isize = nyp * Nsl * nz
    jsize = nxp * Nsl * nz
    s_isize = 3 * isize   ! send/recv buffer size
    s_jsize = 3 * jsize   ! send/recv buffer size

    !!
    !! packing buffer: i-direction
    !!
    !$omp parallel do private(j,k,ptr)
    do j=jbeg, jend
#ifdef _ES
      !cdir nodep,nosync
#endif
      do k=kbeg, kend
        ptr = (k-kbeg)*Nsl + (j-jbeg)*Nsl*(kend-kbeg+1) + 1

        sbuf_ip(        ptr:        ptr+Nsl-1) = Vx(k,iend-Nsl+1:iend,j)
        sbuf_ip(  isize+ptr:  isize+ptr+Nsl-1) = Vy(k,iend-Nsl+1:iend,j)
        sbuf_ip(2*isize+ptr:2*isize+ptr+Nsl-1) = Vz(k,iend-Nsl+1:iend,j)

        sbuf_im(        ptr:        ptr+Nsl-1) = Vx(k,ibeg:ibeg+Nsl-1,j)
        sbuf_im(  isize+ptr:  isize+ptr+Nsl-1) = Vy(k,ibeg:ibeg+Nsl-1,j)
        sbuf_im(2*isize+ptr:2*isize+ptr+Nsl-1) = Vz(k,ibeg:ibeg+Nsl-1,j)

      end do
    end do
    !$omp end parallel do

    !!
    !! Issue send & receive orders: i-direction
    !!
    call mpi_isend( sbuf_ip, s_isize, mpi_precision, itbl(idx+1,idy), 1, mpi_comm_world, ireq1(1), ierr )
    call mpi_isend( sbuf_im, s_isize, mpi_precision, itbl(idx-1,idy), 2, mpi_comm_world, ireq1(2), ierr )
    call mpi_irecv( rbuf_ip, s_isize, mpi_precision, itbl(idx+1,idy), 2, mpi_comm_world, ireq1(3), ierr )
    call mpi_irecv( rbuf_im, s_isize, mpi_precision, itbl(idx-1,idy), 1, mpi_comm_world, ireq1(4), ierr )


    !!
    !! packing buffer: j-direction
    !!
    !$omp parallel do private(ptr)
    do i=ibeg,iend
#ifdef _ES
      !cdir nodep,nosync
#endif
      do k=kbeg,kend
        ptr = (k-kbeg)*Nsl + (i-ibeg)*Nsl*(kend-kbeg+1) + 1

        sbuf_jp(        ptr:        ptr+Nsl-1) = Vx(k,i,jend-Nsl+1:jend)
        sbuf_jp(  jsize+ptr:  jsize+ptr+Nsl-1) = Vy(k,i,jend-Nsl+1:jend)
        sbuf_jp(2*jsize+ptr:2*jsize+ptr+Nsl-1) = Vz(k,i,jend-Nsl+1:jend)

        sbuf_jm(        ptr:        ptr+Nsl-1) = Vx(k,i,jbeg:jbeg+Nsl-1)
        sbuf_jm(  jsize+ptr:  jsize+ptr+Nsl-1) = Vy(k,i,jbeg:jbeg+Nsl-1)
        sbuf_jm(2*jsize+ptr:2*jsize+ptr+Nsl-1) = Vz(k,i,jbeg:jbeg+Nsl-1)

      end do
    end do
    !$omp end parallel do

    !!
    !! Issue send & receive orders: j-direction
    !!
    call mpi_isend( sbuf_jp, s_jsize, mpi_precision, itbl(idx,idy+1), 3, mpi_comm_world, ireq2(1), ierr )
    call mpi_isend( sbuf_jm, s_jsize, mpi_precision, itbl(idx,idy-1), 4, mpi_comm_world, ireq2(2), ierr )
    call mpi_irecv( rbuf_jp, s_jsize, mpi_precision, itbl(idx,idy+1), 4, mpi_comm_world, ireq2(3), ierr )
    call mpi_irecv( rbuf_jm, s_jsize, mpi_precision, itbl(idx,idy-1), 3, mpi_comm_world, ireq2(4), ierr )


    !! Terminate mpi data communication
    call mpi_waitall( 4, ireq1, istatus, ierr )

    !!
    !! restoring the data: i-direction
    !!
    !$omp parallel do private(ptr,i,j,k)
    do j=jbeg, jend
      do k=kbeg, kend

        ptr = ( (k-kbeg) + (j-jbeg)*nz ) * Nsl + 1

        do i=1,Nsl
          Vx(k,iend+i,j) = rbuf_ip(        ptr+i-1)
          Vy(k,iend+i,j) = rbuf_ip(  isize+ptr+i-1)
          Vz(k,iend+i,j) = rbuf_ip(2*isize+ptr+i-1)

          Vx(k,ibeg-Nsl+i-1,j) = rbuf_im(        ptr+i-1)
          Vy(k,ibeg-Nsl+i-1,j) = rbuf_im(  isize+ptr+i-1 )
          Vz(k,ibeg-Nsl+i-1,j) = rbuf_im(2*isize+ptr+i-1 )
        end do

      end do
    end do
    !$omp end parallel do

    !! Terminate mpi data communication
    call mpi_waitall( 4, ireq2, istatus, ierr )

    !!
    !! restoring the data: j-direction
    !!
    !$omp parallel do private(ptr,i,j,k)
    do i=ibeg,iend
      do k=kbeg,kend

        ptr = (k-kbeg)*Nsl + (i-ibeg)*Nsl*(kend-kbeg+1) + 1

        do j=1, Nsl
          Vx(k,i,jend+j) = rbuf_jp(ptr+j-1)
          Vy(k,i,jend+j) = rbuf_jp(jsize+ptr+j-1)
          Vz(k,i,jend+j) = rbuf_jp(2*jsize+ptr+j-1)

          Vx(k,i,jbeg-Nsl+j-1) = rbuf_jm(ptr+j-1)
          Vy(k,i,jbeg-Nsl+j-1) = rbuf_jm(jsize+ptr+j-1)
          Vz(k,i,jbeg-Nsl+j-1) = rbuf_jm(2*jsize+ptr+j-1)
        end do

      end do
    end do
    !$omp end parallel do

    call pwatch__off( "global__comm_vel" )

  end subroutine global__comm_vel
  !! ---------------------------------------------------------------------------------------------------------------------------- !!

  !! ---------------------------------------------------------------------------------------------------------------------------- !!
  !>
  !! Data buffring & communication for stress tensor
  !<
  !!
  subroutine global__comm_stress()

    integer :: isize, jsize, s_isize, s_jsize
    integer :: ierr
    integer :: istatus( mpi_status_size, 4 )
    integer :: ireq1(4), ireq2(4)
    integer :: i, j, k, ptr

    if( myid >= nproc ) return

    call pwatch__on( "global__comm_stress" )

    !! unit buffer size
    isize = nyp * Nsl * nz
    jsize = nxp * Nsl * nz
    s_isize = 3 * isize   ! send/recv buffer size
    s_jsize = 3 * jsize   ! send/recv buffer size

    !!
    !! packing buffer: i-direction ( Sxx, Sxy, Sxz )
    !!
    !$omp parallel do private(i,k,ptr)
    do j=jbeg, jend
#ifdef _ES
      !cdir nodep,nosync
#endif
      do k=kbeg, kend
        ptr = (k-kbeg)*Nsl + (j-jbeg)*Nsl*(kend-kbeg+1) + 1

        sbuf_ip(        ptr:        ptr+Nsl-1) = Sxx(k,iend-Nsl+1:iend,j)
        sbuf_ip(  isize+ptr:  isize+ptr+Nsl-1) = Sxy(k,iend-Nsl+1:iend,j)
        sbuf_ip(2*isize+ptr:2*isize+ptr+Nsl-1) = Sxz(k,iend-Nsl+1:iend,j)

        sbuf_im(        ptr:        ptr+Nsl-1) = Sxx(k,ibeg:ibeg+Nsl-1,j)
        sbuf_im(  isize+ptr:  isize+ptr+Nsl-1) = Sxy(k,ibeg:ibeg+Nsl-1,j)
        sbuf_im(2*isize+ptr:2*isize+ptr+Nsl-1) = Sxz(k,ibeg:ibeg+Nsl-1,j)

      end do
    end do
    !$omp end parallel do

    !!
    !! Issue send & receive orders; i-direction
    !!
    !! i-direction
    call mpi_isend( sbuf_ip, s_isize, mpi_precision, itbl(idx+1,idy), 5, mpi_comm_world, ireq1(1), ierr )
    call mpi_isend( sbuf_im, s_isize, mpi_precision, itbl(idx-1,idy), 6, mpi_comm_world, ireq1(2), ierr )
    call mpi_irecv( rbuf_ip, s_isize, mpi_precision, itbl(idx+1,idy), 6, mpi_comm_world, ireq1(3), ierr )
    call mpi_irecv( rbuf_im, s_isize, mpi_precision, itbl(idx-1,idy), 5, mpi_comm_world, ireq1(4), ierr )

    !!
    !! packing buffer: j-direction ( Syy, Syz, Sxy )
    !!
    !$omp parallel do private(ptr)
    do i=ibeg,iend
      !cdir nodep,nosync
      do k=kbeg,kend
        ptr = (k-kbeg)*Nsl + (i-ibeg)*Nsl*(kend-kbeg+1) + 1

        sbuf_jp(        ptr:        ptr+Nsl-1) = Syy(k,i,jend-Nsl+1:jend)
        sbuf_jp(  jsize+ptr:  jsize+ptr+Nsl-1) = Syz(k,i,jend-Nsl+1:jend)
        sbuf_jp(2*jsize+ptr:2*jsize+ptr+Nsl-1) = Sxy(k,i,jend-Nsl+1:jend)

        sbuf_jm(        ptr:        ptr+Nsl-1) = Syy(k,i,jbeg:jbeg+Nsl-1)
        sbuf_jm(  jsize+ptr:  jsize+ptr+Nsl-1) = Syz(k,i,jbeg:jbeg+Nsl-1)
        sbuf_jm(2*jsize+ptr:2*jsize+ptr+Nsl-1) = Sxy(k,i,jbeg:jbeg+Nsl-1)

      end do
    end do
    !$omp end parallel do


    !!
    !! Issue send & receive orders; j-direction
    !!
    call mpi_isend( sbuf_jp, s_jsize, mpi_precision, itbl(idx,idy+1), 7, mpi_comm_world, ireq2(1), ierr )
    call mpi_isend( sbuf_jm, s_jsize, mpi_precision, itbl(idx,idy-1), 8, mpi_comm_world, ireq2(2), ierr )
    call mpi_irecv( rbuf_jp, s_jsize, mpi_precision, itbl(idx,idy+1), 8, mpi_comm_world, ireq2(3), ierr )
    call mpi_irecv( rbuf_jm, s_jsize, mpi_precision, itbl(idx,idy-1), 7, mpi_comm_world, ireq2(4), ierr )


    !! Terminate mpi data communication
    call mpi_waitall( 4, ireq1, istatus, ierr )

    !!
    !! restore the data: i-direction
    !!
    !$omp parallel do private(ptr,i,j,k)
    do j=jbeg, jend
      do k=kbeg, kend

        ptr = (k-kbeg)*Nsl + (j-jbeg)*Nsl*(kend-kbeg+1) + 1

        do i=1, Nsl

          Sxx(k,iend+i,j) = rbuf_ip(        ptr+i-1)
          Sxy(k,iend+i,j) = rbuf_ip(  isize+ptr+i-1)
          Sxz(k,iend+i,j) = rbuf_ip(2*isize+ptr+i-1)

          Sxx(k,ibeg-Nsl+i-1,j) = rbuf_im(       ptr+i-1)
          Sxy(k,ibeg-Nsl+i-1,j) = rbuf_im(  isize+ptr+i-1)
          Sxz(k,ibeg-Nsl+i-1,j) = rbuf_im(2*isize+ptr+i-1)
        end do

      end do
    end do
    !$omp end parallel do


    !! Terminate mpi data communication
    call mpi_waitall( 4, ireq2, istatus, ierr )

    !!
    !! restore the data: j-direction
    !!
    !$omp parallel do private(ptr,i,j,k)
    do i=ibeg,iend
      do k=kbeg,kend

        ptr = (k-kbeg)*Nsl + (i-ibeg)*Nsl*(kend-kbeg+1) + 1

        do j=1, Nsl
          Syy(k,i,jend+j) = rbuf_jp(        ptr+j-1)
          Syz(k,i,jend+j) = rbuf_jp(  jsize+ptr+j-1)
          Sxy(k,i,jend+j) = rbuf_jp(2*jsize+ptr+j-1)

          Syy(k,i,jbeg-Nsl+j-1) = rbuf_jm(        ptr+j-1)
          Syz(k,i,jbeg-Nsl+j-1) = rbuf_jm(  jsize+ptr+j-1)
          Sxy(k,i,jbeg-Nsl+j-1) = rbuf_jm(2*jsize+ptr+j-1)
        end do

      end do
    end do
    !$omp end parallel do

    call pwatch__off( "global__comm_stress" )

  end subroutine global__comm_stress
  !! ---------------------------------------------------------------------------------------------------------------------------- !!



  !! ---------------------------------------------------------------------------------------------------------------------------- !!
  !>
  !! check if the voxcel location is inside the MPI node
  !<
  logical function inside_node ( i, j, k )

    integer, intent(in) :: i, j, k

    !! ----

    if( ibeg <= i .and. i<= iend .and. &
        jbeg <= j .and. j<= jend .and. &
        kbeg <= k .and. k<= kend ) then
      inside_node = .true.
    else
      inside_node = .false.
    end if

  end function inside_node

  !! ---------------------------------------------------------------------------------------------------------------------------- !!

  !! ---------------------------------------------------------------------------------------------------------------------------- !!
  subroutine global__checkpoint( io )

    integer, intent(in) :: io

    write(io) title(1:80)
    write(io) exedate
    write(io) mpi_precision

    write(io) nproc_x
    write(io) nproc_y
    write(io) nproc
    write(io) nx, ny, nz
    write(io) nxp,nyp

    write(io) nt
    write(io) na
    write(io) dx, dy, dz
    write(io) dt
    write(io) xbeg, xend
    write(io) ybeg, yend
    write(io) zbeg, zend
    write(io) tbeg, tend
    write(io) ibeg, iend
    write(io) jbeg, jend
    write(io) kbeg, kend
    write(io) ibeg_m, iend_m
    write(io) jbeg_m, jend_m
    write(io) kbeg_m, kend_m
    write(io) ibeg_k, iend_k
    write(io) jbeg_k, jend_k
    write(io) kbeg_k, kend_k
    write(io) odir(1:256)
    write(io) clon
    write(io) clat
    write(io) phi
    write(io) xc(ibeg_m:iend_m)
    write(io) yc(jbeg_m:jend_m)
    write(io) zc(kbeg_m:kend_m)
    write(io) kbeg_a(ibeg:iend,jbeg:jend)

    deallocate( sbuf_ip, sbuf_im, sbuf_jp, sbuf_jm )
    deallocate( rbuf_ip, rbuf_im, rbuf_jp, rbuf_jm )
    deallocate( xc, yc, zc )
    deallocate( kbeg_a )

  end subroutine global__checkpoint

  subroutine global__restart( io )

    integer, intent(in) :: io
    integer :: nl3
    read(io) title(1:80)
    read(io) exedate
    read(io) mpi_precision

    read(io) nproc_x
    read(io) nproc_y
    read(io) nproc
    read(io) nx, ny, nz
    read(io) nxp,nyp

    nl3 = (Nl/2)*3
    allocate( itbl(-1:nproc_x, -1:nproc_y) )
    allocate( sbuf_ip( nyp*nz*nl3 ), sbuf_im( nyp*nz*nl3 ) )
    allocate( rbuf_ip( nyp*nz*nl3 ), rbuf_im( nyp*nz*nl3 ) )
    allocate( sbuf_jp( nxp*nz*nl3 ), sbuf_jm( nxp*nz*nl3 ) )
    allocate( rbuf_jp( nxp*nz*nl3 ), rbuf_jm( nxp*nz*nl3 ) )
    call set_mpi_table

    read(io) nt
    read(io) na
    read(io) dx, dy, dz
    read(io) dt
    read(io) xbeg, xend
    read(io) ybeg, yend
    read(io) zbeg, zend
    read(io) tbeg, tend
    read(io) ibeg, iend
    read(io) jbeg, jend
    read(io) kbeg, kend
    read(io) ibeg_m, iend_m
    read(io) jbeg_m, jend_m
    read(io) kbeg_m, kend_m
    read(io) ibeg_k, iend_k
    read(io) jbeg_k, jend_k
    read(io) kbeg_k, kend_k
    read(io) odir(1:256)
    read(io) clon
    read(io) clat
    read(io) phi


    allocate( xc(ibeg_m:iend_m), yc(jbeg_m:jend_m), zc(kbeg_m:kend_m) )
    read(io) xc(ibeg_m:iend_m)
    read(io) yc(jbeg_m:jend_m)
    read(io) zc(kbeg_m:kend_m)

    allocate( kbeg_a(ibeg:iend, jbeg:jend) )
    read(io) kbeg_a(ibeg:iend,jbeg:jend)

  end subroutine global__restart

  subroutine set_mpi_table
    integer :: i
    integer :: ii, jj
    !! 2D communication table
    !> @see 2013-0439

    itbl(-1:nproc_x,-1:nproc_y) = MPI_PROC_NULL
    do i=0, nproc-1

      ii = mod( i, nproc_x )
      jj = i / nproc_x

      itbl(ii,jj) = i
    end do


    !! location of this process
    idx = mod( myid, nproc_x )
    idy = myid / nproc_x
  end subroutine set_mpi_table


end module m_global
!! ------------------------------------------------------------------------------------------------------------------------------ !!
