!! ----------------------------------------------------------------------------------------------------------------------------- !!
!>
!! terminal/logfile report
!!
!! @copyright
!!   Copyright 2013-2021 Takuto Maeda. All rights reserved. This project is released under the MIT license.
!<
!! ----
#include "m_debug.h"
module m_report

  !! -- Dependency
  use m_std
  use m_debug
  use m_global
  use m_pwatch
  use m_kernel
  use m_readini

  !! -- Declarations
  implicit none
  private
  save

  !! -- procedures
  public :: report__setup
  public :: report__progress
  public :: report__terminate
  public :: report__checkpoint
  public :: report__restart

  integer :: ntdec_r

  !! -- internal parameter
  integer, parameter :: terminal_output_node = 0

  !<< Lapse Time Measurement >>
  integer  :: timcount, timcount0, timprev
  real(SP) :: ttotal

contains

  !! --------------------------------------------------------------------------------------------------------------------------- !!
  !>
  !! Initialization, welcome message to terminal, open logfile
  !<
  !! ----
  subroutine report__setup( io_prm )

    integer, intent(in) :: io_prm

    !! --
    integer :: crate
    real(SP) :: mem_all, mem_node
    real(SP) :: c, r
    !! ----

    call readini( io_prm, 'ntdec_r', ntdec_r, 10 )


    if( myid == terminal_output_node ) then

      write(STDERR,*)
      write(STDERR,'(A)') " ------------------------------------------------------------------------------"
      if( benchmark_mode ) then
        write(STDERR,'(A)') "  SWPC_3D (benchmark mode)                                                    "
      else if ( pw_mode ) then
        write(STDERR,'(A)') "  SWPC_3D (plane wave mode)                                                   "
      else if ( bf_mode ) then
        write(STDERR,'(A)') "  SWPC_3D (body force mode)                                                   "
      else if ( green_mode ) then
        write(STDERR,'(A)') "  SWPC_3D (Green's function mode)                                             "
      else
        write(STDERR,'(A)') "  SWPC_3D                                                                     "
      end if
      write(STDERR,'(A)') " ------------------------------------------------------------------------------"

    end if


    call memory_size( nproc_x, nproc_y, nx, ny, nz, nm, na, mem_all, mem_node )
    call fdm_cond_stability ( real(dx), real(dy), real(dz), vmax, dt,   c )
    call fdm_cond_wavelength( real(dx), real(dy), real(dz), vmin, fmax, r )

    if( myid == terminal_output_node ) then
      write(STDERR,*)
      write(STDERR,'(A,I8,A,I6,A,I6)') "  Grid Size               : ", nx, " x ", ny ," x ", nz
      write(STDERR,'(A,I8,A,I4)'     ) "  MPI Partitioning        : ", nproc_x, " x ", nproc_y
      write(STDERR,'(A,F15.3,A)'     ) "  Total Memory Size       : ", mem_all,  "  [GiB]"
      write(STDERR,'(A,F15.3,A)'     ) "  Node Memory Size        : ", mem_node, "  [GiB]"
      write(STDERR,'(A,F15.3,A)'     ) "  Stability  Condition c  : ", c,        "  (c<1)"
      write(STDERR,'(A,F15.3,A)'     ) "  Wavelength Condition r  : ", r       , "  (r>5-10)"
      write(STDERR,'(A,F15.3,A)'     ) "  Minimum velocity        : ", vmin,     "  [km/s]"
      write(STDERR,'(A,F15.3,A)'     ) "  Maximum velocity        : ", vmax,     "  [km/s]"
      write(STDERR,'(A,F15.3,A)'     ) "  Maximum frequency       : ", fmax,     "  [Hz]"
      write(STDERR,*)
      write(STDERR,'(A)') " ------------------------------------------------------------------------------"
      write(STDERR,*)

      if( r < 5   ) then
        call info( 'wavelength condition is violated! ' )
        call info( 'use smaller grid and/or decrease maximum frequency' )
      end if

      if( c < 0.5 ) then
        call info( 'time step is too small!' )
        call info( 'consider increase time step up to twice' )
      end if

      if( c > 1.0 ) then
        call info( 'stability condition is violated!' )
        call info( 'use smaller time step and/or decrease max velocity' )
        call assert( c <= 1.0 )
      end if

    end if

    !! Initialize elapsed time counter
    if( myid == terminal_output_node ) then

      call system_clock( timcount, crate )
      timcount0 = timcount
      timprev   = timcount
      ttotal = 0

    end if


  end subroutine report__setup
  !! --------------------------------------------------------------------------------------------------------------------------- !!

  !! --------------------------------------------------------------------------------------------------------------------------- !!
  !>
  !! Show progres to the terminal
  !<
  !! ----
  subroutine report__progress( it )

    integer, intent(in) :: it
    real(SP) :: vxm, vym, vzm
    real(SP) :: vxa, vya, vza
    integer  :: ierr
    real(SP) :: etas
    integer  :: etah, etam, etasi
    real(SP) :: tstep
    integer  :: crate, cmax
    real(SP), parameter :: TOL = 1e5
    !! ----

    if( mod(it,ntdec_r) /= 0 ) return

    call pwatch__on("report__progress")

    call kernel__vmax( vxm, vym, vzm )

    call mpi_allreduce( vxm, vxa, 1, MPI_REAL, MPI_MAX, mpi_comm_world, ierr )
    call mpi_allreduce( vym, vya, 1, MPI_REAL, MPI_MAX, mpi_comm_world, ierr )
    call mpi_allreduce( vzm, vza, 1, MPI_REAL, MPI_MAX, mpi_comm_world, ierr )


    !! check numerical divergence
    if( max( vxa, vya, vza ) * UC > TOL ) then
      if( myid == terminal_output_node ) then
        write(STDERR,'(A)') 'numerical divergence detected with max amp =  ', max( vxa, vya, vza ) * UC, '[(m/s)/moment]'
        write(STDERR,'(A)') 'aborting ... '
      end if

      call mpi_finalize( ierr )
      stop
    end if


    if( myid == terminal_output_node ) then

      !! convert to the physical domain
      vxa = vxa * UC * M0
      vya = vya * UC * M0
      vza = vza * UC * M0

      !!
      !! eta count
      !!
      call system_clock( timcount, crate, cmax )
      if( timcount >= timprev ) then
        tstep = real( timcount - timprev ) / real( crate )
      else
        tstep = real( cmax + timcount - timprev ) / real( crate )
      end if

      ttotal = ttotal + tstep

      etas   = real(nt-it)/ real(it) * ttotal

      etah = int( etas/(   60*60) ); etas = etas - etah   *60*60
      etam = int( etas/(      60) ); etas = etas - etam      *60
      etasi = int(etas)
      timprev = timcount

      write(STDERR,'(A,I7.7,  A,F6.3,A,   A,I3.3,A,I2.2,A,I2.2,A, 3(ES9.2,A))') &
          "  it=", it, ",", &
          ttotal/it, " s/loop,", &
          " eta ", etah,":",etam,":",etasi ,", (", vxa, " ",vya," ", vza," )"

    end if

    call pwatch__off("report__progress")

  end subroutine report__progress
  !! --------------------------------------------------------------------------------------------------------------------------- !!

  !! --------------------------------------------------------------------------------------------------------------------------- !!
  subroutine report__terminate

    if( myid == terminal_output_node ) then

      write(STDERR,*)
      write(STDERR,'(A)') " ------------------------------------------------------------------------------"
      write(STDERR,*) ""
      write(STDERR,'(A,F15.3,A)')        "  Total time             : ", ttotal, " s"
      write(STDERR,*)
      write(STDERR,'(A)') " ------------------------------------------------------------------------------"

    end if

  end subroutine report__terminate
  !! --------------------------------------------------------------------------------------------------------------------------- !!

  !! --------------------------------------------------------------------------------------------------------------------------- !!
  subroutine report__checkpoint( io )

    integer, intent(in) :: io

    write( io ) ntdec_r
    write( io ) ttotal

  end subroutine report__checkpoint
  !! --------------------------------------------------------------------------------------------------------------------------- !!

  !! --------------------------------------------------------------------------------------------------------------------------- !!
  subroutine report__restart( io )

    integer, intent(in) :: io
    integer :: crate

    read( io ) ntdec_r
    read( io ) ttotal

    !! Initialize elapsed time counter
    if( myid == terminal_output_node ) then

      call system_clock( timcount, crate )
      timcount0 = timcount
      timprev   = timcount

    end if

  end subroutine report__restart
  !! --------------------------------------------------------------------------------------------------------------------------- !!

end module m_report
!! ----------------------------------------------------------------------------------------------------------------------------- !!
