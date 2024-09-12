!! ----------------------------------------------------------------------------------------------------------------------------- !!
!>
!! global control parameters, shared arrays and MPI communication
!!
!! Copyright 2013-2024 Takuto Maeda. All rights reserved. This project is released under the MIT license.
!<
!! ----
#include "../shared/m_debug.h"
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
  real(SP)            :: UC = 10.0**(-12)                           !< Conventional -> SI unit for moment tensor 
  integer,  parameter :: MP = DP                                    !< DP for mixed precision, SP for pure single precision
  integer,  parameter :: NM = 3                                     !< Number of memory variables
  integer,  parameter :: NBD = 9                                    !< Number of boundary depths to be memorized
  !!
  !! global arrays
  !!
  real(MP), allocatable :: Vy(:,:)                                  !<  velocity components
  real(MP), allocatable :: Syz(:,:),   Sxy(:,:)                     !<  shear  stress components
  real(SP), allocatable :: Ryz(:,:,:), Rxy(:,:,:)                   !<  memory variables: shear  components
  real(SP), allocatable :: rho(:,:),   lam(:,:), mu(:,:)            !<  density, relaxed moduli
  real(SP), allocatable :: taup(:,:),  taus(:,:)                    !<  creep/relax time ratio based on tau-method for attenuation
  real(SP), allocatable :: ts(:)                                    !<  relaxation time of visco-elastic medium

  !!
  !! title, date
  !!
  character(80)         :: title                                    !<  execution title, used in filename and headers
  integer               :: exedate                                  !<  date and time by seconds from 1970/1/1 0:0:0

  !!
  !! global control parameters
  !!
  integer               :: nx, nz                                   !<  space grid number (global)
  integer               :: nt                                       !<  time grid number
  integer               :: na                                       !<  absorber thickness
  real(MP)              :: dx, dz                                   !<  space grid width
  real(SP)              :: dt                                       !<  time  grid width
  real(SP)              :: xbeg, xend                               !<  global coordinate: x start / end
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
  integer               :: nproc_x                                  !<  total   numbers of process
  integer               :: nxp                                      !<  space grid number in the assigned node
  integer               :: myid                                     !<  MPI node number
  integer               :: idx                                      !<  2D horizontal division ID
  integer, allocatable  :: itbl(:)                                  !<  node layout table
  integer               :: ibeg, iend                               !<  i-region in the node
  integer               :: kbeg, kend                               !<  k-region in the node
  integer               :: ibeg_m, iend_m                           !<  i- memory allocation area
  integer               :: kbeg_m, kend_m                           !<  k- memory allocation area
  integer               :: ibeg_k, iend_k                           !<  i- kernel integration area without absorption band
  integer               :: kbeg_k, kend_k                           !<  k- kernel integration area without absorption band
  integer, allocatable  :: kbeg_a(:)                                !<  absorbing boundary region
  integer               :: ipad, kpad                               !<  memory padding size for optimization
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
  integer,  allocatable :: kfs(:)                                   !<  free surface depth grid in the node
  integer,  allocatable :: kob(:)                                   !<  ocean bottom depth grid in the node
  integer,  allocatable :: kfs_top(:), kfs_bot(:)                   !<  region in which 2nd-order FDM is applied for free surface
  integer,  allocatable :: kob_top(:), kob_bot(:)                   !<  region in which 2nd-order FDM is applied for ocean bottom
  real(SP), allocatable :: bddep(:,:)

  !!
  !! map coordinate: use dummy value for cartesian problem
  !!
  real(SP) :: clon                                                  !< center longitude
  real(SP) :: clat                                                  !< center latitude
  real(SP) :: phi                                                   !< azimuth
  real(SP), allocatable :: xc(:), zc(:)


  !! Source: remember first source grid location in lat/lon form
  real(SP) :: evlo
  real(SP) :: evla
  real(SP) :: evdp !< unit:km
  real(SP) :: mxx0, myy0, mzz0, myz0, mxz0, mxy0
  real(SP) :: fx0, fy0, fz0
  real(SP) :: otim
  real(SP) :: sx0, sy0
  
  !!
  !! Benchmark
  !!
  logical :: benchmark_mode

  !! special modes
  logical :: pw_mode                                                !< plane wave mode
  logical :: bf_mode                                                !< body force mode


  !!
  !! private variables
  !!
  real(MP), private, allocatable :: sbuf_ip(:), sbuf_im(:)          !<  mpi send buffer
  real(MP), private, allocatable :: rbuf_ip(:), rbuf_im(:)          !<  mpi recv buffer
  integer :: mpi_precision


  !! fullspace-mdoe
!  logical :: fullspace_mode

  !! ----
  private :: inside_node
  private :: set_mpi_table

contains

  !! ---------------------------------------------------------------------------------------------------------------------------- !!
  subroutine global__readprm( io_prm )
    integer, intent(in) :: io_prm

    call readini( io_prm, 'benchmark_mode', benchmark_mode, .false. )

    call readini( io_prm, 'title',   title, 'swpc_sh'  )
    call readini( io_prm, 'nproc_x', nproc_x,   1      )
    call readini( io_prm, 'nx',      nx,      256      )
    call readini( io_prm, 'nz',      nz,      256      )
    call readini( io_prm, 'nt',      nt,      1000     )
    call readini( io_prm, 'ipad',    ipad,    0        )
    call readini( io_prm, 'kpad',    kpad,    0        )
    call readini( io_prm, 'odir',    odir,  './out'    )

    if( benchmark_mode ) then
      dx = 0.5
      dz = 0.5
      dt = 0.04
      na = 20
      xbeg = -nx/2 * real(dx)
      zbeg = -30   * real(dz)
      tbeg = 0.0
      clon = 139.7604
      clat = 35.7182
      phi  = 0.0
      abc_type = 'pml'
!      fullspace_mode = .false.
    else
      call readini( io_prm, 'dx',      dx,      0.5_MP         )
      call readini( io_prm, 'dz',      dz,      0.5_MP         )
      call readini( io_prm, 'dt',      dt,      0.01           )
      call readini( io_prm, 'na',      na,      20             )
      call readini( io_prm, 'xbeg',    xbeg,   -nx/2*real(dx)  )
      call readini( io_prm, 'zbeg',    zbeg,   -30*real(dz)    )
      call readini( io_prm, 'tbeg',    tbeg,    0.0            )
      call readini( io_prm, 'clon',    clon,  139.7604         )
      call readini( io_prm, 'clat',    clat,   35.7182         )
      call readini( io_prm, 'phi',     phi,     0.0            )
      call readini( io_prm, 'abc_type', abc_type, 'pml'        )
!      call readini( io_prm, 'fullspace_mode', fullspace_mode, .false.    )
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
    integer :: err
    !! ----

    call pwatch__on( "global__setup" ) !! measure from here

    !!
    !! MPI status check
    !!
    call mpi_comm_rank( mpi_comm_world, myid, err )
    if( MP == DP ) then
      mpi_precision = MPI_DOUBLE_PRECISION
    else
      mpi_precision = MPI_REAL
    end if


    !!
    !! read key parameters
    !!
    call global__readprm( io_prm )

    !!
    !! obtain date by unixtime: seconds measured from 1970/1/1 0:0:0
    !!
    if(myid == 0) call daytim__getdate( exedate )
    call mpi_bcast( exedate, 1, MPI_INTEGER, 0, mpi_comm_world, err )


    !!
    !! derived parameters
    !!
    xend = xbeg + nx * real(dx)
    zend = zbeg + nz * real(dz)
    tend = tbeg + nt * dt

    call pwatch__off( "global__setup" ) !! measure from here

  end subroutine global__setup

  !! --------------------------------------------------------------------------------------------------------------------------- !!

  !! --------------------------------------------------------------------------------------------------------------------------- !!
  !>
  !! Common parameter setup, needed only for starting
  !<
  !! --
  subroutine global__setup2

    integer :: i, k
    integer :: err, nproc_exe
    integer :: mx, proc_x

    call pwatch__on( "global__setup2" ) !! measure from here

    !!
    !! output directory create (if it does not exist)
    !!
    call system__call('mkdir -p ' // trim(odir) // '> /dev/null 2>&1' )


    call mpi_comm_size( mpi_comm_world, nproc_exe, err )
    call assert( nproc_x == nproc_exe )


    mx = mod(nx, nproc_x)
    proc_x = myid
    if( proc_x <= nproc_x - mx + 1 ) then
      nxp = (nx - mx)/nproc_x
    else
      nxp = (nx - mx)/nproc_x + 1
    end if    


    !! MPI coordinate
    allocate( itbl(-1:nproc_x) )

    allocate(sbuf_ip(2 * nz), source=0.0_MP)
    allocate(sbuf_im(2 * nz), source=0.0_MP)
    allocate(rbuf_ip(2 * nz), source=0.0_MP)
    allocate(rbuf_im(2 * nz), source=0.0_MP)

    !! MPI communication table
    call set_mpi_table

    !! computation region in this node (#244)
    if ( proc_x <= nproc_x -mx - 1 ) then
      ibeg =  proc_x      * (nx - mx) / nproc_x + 1
      iend = (proc_x + 1) * (nx - mx) / nproc_x
    else
      ibeg =  proc_x      * ((nx - mx) / nproc_x + 1 ) - (nproc_x - mx) + 1
      iend = (proc_x + 1) * ((nx - mx) / nproc_x + 1 ) - (nproc_x - mx)
    end if
    kbeg = 1
    kend = nz

    
    !! memory requirements including margin for MPI/boundary conditions
    !! stress drop also requires sleeve area
    ibeg_m = ibeg - 3
    iend_m = iend + 3 + ipad
    kbeg_m = kbeg - 3
    kend_m = kend + 3 + kpad


    !!
    !! coordinate setting
    !!
    allocate( xc(ibeg_m:iend_m), zc(kbeg_m:kend_m) )

    do i=ibeg_m, iend_m
      xc(i) = i2x( i, xbeg, real(dx) )
    end do
    do k=kbeg_m, kend_m
      zc(k) = k2z( k, zbeg, real(dz) )
    end do



    !!
    !! Absorbing region definition
    !!
    allocate( kbeg_a(ibeg_m:iend_m) )
    do i=ibeg_m, iend_m
      if( i <= na .or. nx-na+1 <= i ) then
        kbeg_a(i) = kbeg
      else
        kbeg_a(i) = kend-na+1
      end if
    end do

    ibeg_k = ibeg
    iend_k = iend
    kbeg_k = kbeg
    kend_k = kend

    
    if( abc_type == 'pml' ) then
!      if( fullspace_mode ) kbeg_k = na+1
      if( iend <= na ) then ! no kernel integration
        ibeg_k = iend+1
      else if ( ibeg <= na ) then ! pertial kernel
        ibeg_k = na+1
      end if
      if( ibeg >= nx-na+1 ) then ! no kernel integartion
        iend_k = ibeg-1
      else if ( iend >= nx-na+1 ) then
        iend_k = nx-na
      end if
      kend_k = nz-na
    end if


    call pwatch__off( "global__setup2" ) !! measure from here


  end subroutine global__setup2
  !! ---------------------------------------------------------------------------------------------------------------------------- !!

  !! ---------------------------------------------------------------------------------------------------------------------------- !!
  !>
  !! Data buffring & communication for velocity vector
  !!
  !! @see
  !! 2013-0420, 2013-0421
  !<
  !!
  subroutine global__comm_vel()

    integer :: err
    integer :: istatus( mpi_status_size, 4 )
    integer :: req(4)
    !! ----

    if( myid >= nproc_x ) return

    call pwatch__on( "global__comm_vel" )

    call mpi_irecv(rbuf_ip, 2*nz, mpi_precision, itbl(idx+1), 1, mpi_comm_world, req(1), err)
    call mpi_irecv(rbuf_im, 1*nz, mpi_precision, itbl(idx-1), 2, mpi_comm_world, req(2), err)

    sbuf_ip(1:  nz) = reshape(Vy(1:nz,iend:iend  ), (/  nz/))
    call mpi_isend(sbuf_ip, 1*nz, mpi_precision, itbl(idx+1), 2, mpi_comm_world, req(3), err)
    sbuf_im(1:2*nz) = reshape(Vy(1:nz,ibeg:ibeg+1), (/2*nz/))
    call mpi_isend(sbuf_im, 2*nz, mpi_precision, itbl(idx-1), 1, mpi_comm_world, req(4), err)

    call mpi_waitall( 4, req, istatus, err )

    Vy(1:nz,iend+1:iend+2) = reshape(rbuf_ip(   1:2*nz), (/nz,2/))
    Vy(1:nz,ibeg-1:ibeg-1) = reshape(rbuf_im(   1:  nz), (/nz,1/))

    call pwatch__off( "global__comm_vel" )

  end subroutine global__comm_vel
  !! ---------------------------------------------------------------------------------------------------------------------------- !!

  !! ---------------------------------------------------------------------------------------------------------------------------- !!
  !>
  !! Data buffring & communication for stress tensor
  !<
  !!
  subroutine global__comm_stress()

    integer :: err
    integer :: istatus( mpi_status_size, 4)
    integer :: req(4)
    !! ----

    if( myid >= nproc_x ) return

    call pwatch__on( "global__comm_stress" )

    call mpi_irecv(rbuf_ip, 1*nz, mpi_precision, itbl(idx+1), 3, mpi_comm_world, req(1), err)
    call mpi_irecv(rbuf_im, 2*nz, mpi_precision, itbl(idx-1), 4, mpi_comm_world, req(2), err)

    sbuf_ip(1:2*nz) = reshape(Sxy(1:nz,iend-1:iend), (/2*nz/))
    call mpi_isend(sbuf_ip, 2*nz, mpi_precision, itbl(idx+1), 4, mpi_comm_world, req(3), err)

    sbuf_im(1:  nz) = reshape(Sxy(1:nz,ibeg  :ibeg), (/  nz/))
    call mpi_isend(sbuf_im, 1*nz, mpi_precision, itbl(idx-1), 3, mpi_comm_world, req(4), err)

    call mpi_waitall( 4, req, istatus, err )

    !! Resore the data
    Sxy(kbeg:kend,iend+1:iend+1) = reshape(rbuf_ip(1:  nz), (/nz,1/))
    Sxy(kbeg:kend,ibeg-2:ibeg-1) = reshape(rbuf_im(1:2*nz), (/nz,2/))

    call pwatch__off( "global__comm_stress" )

  end subroutine global__comm_stress
  !! ---------------------------------------------------------------------------------------------------------------------------- !!


  !! ---------------------------------------------------------------------------------------------------------------------------- !!
  !>
  !! check if the voxcel location is inside the MPI node
  !<
  logical function inside_node ( i, k )

    integer, intent(in) :: i, k

    !! ----

    if( ibeg <= i .and. i<= iend .and. &
        kbeg <= k .and. k<= kend ) then
      inside_node = .true.
    else
      inside_node = .false.
    end if

  end function inside_node

  !! ---------------------------------------------------------------------------------------------------------------------------- !!

  subroutine set_mpi_table
    integer :: i
    integer :: ii
    !! 2D communication table
    !> @see 2013-0439

    itbl(-1:nproc_x) = MPI_PROC_NULL
    do i=0, nproc_x-1

      ii = mod( i, nproc_x )
      itbl(ii) = i

    end do

    !! location of this process
    idx = mod( myid, nproc_x )
  end subroutine set_mpi_table


end module m_global
!! ------------------------------------------------------------------------------------------------------------------------------ !!
