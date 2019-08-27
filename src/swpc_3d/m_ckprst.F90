!! ----------------------------------------------------------------------------------------------------------------------------- !!
!>
!! Checkpoint and restart
!!
!! @copyright
!!   Copyright 2013-2019 Takuto Maeda. All rights reserved. This project is released under the MIT license.
!<
!! ----
#include "m_debug.h"
module m_ckprst

  use m_std
  use m_debug
  use m_system
  use m_global
  use m_kernel
  use m_medium
  use m_absorb
  use m_report
  use m_source
  use m_output
  use m_pwatch
  use m_readini
  implicit none
  private
  save

  public :: ckprst__checkpoint
  public :: ckprst__restart
  public :: ckprst__setup
  public :: ckprst__filedelete

  character(99) :: ckpdir
  real(SP) :: ckp_time
  integer :: ckp_interval
  integer :: t0, t1 ! elapsed time
  logical :: is_ckp

contains

  !! --------------------------------------------------------------------------------------------------------------------------- !!
  subroutine ckprst__checkpoint(it)

    integer, intent(in) :: it

    character(256) :: fn_ckp
    integer :: io_ckp
    integer :: ierr
    logical :: arrived_checkpoint = .false.
    integer :: t_rate, t_max
    real(SP) :: t_elapse
    !! --

    if( .not. is_ckp ) return
    if( mod(it, ckp_interval) /= 0 ) return

    call pwatch__on( "ckprst__checkpoint" )


    if( myid == 0 ) then
      call system_clock( t1, t_rate, t_max )
      if( t1 < t0 ) then
        t_elapse = ( t_max - t0 + t1 ) / real( t_rate )
      else
        t_elapse = ( t1 - t0 ) / real( t_rate )
      end if

      if( t_elapse > ckp_time ) arrived_checkpoint = .true.
    end if
    call mpi_bcast( arrived_checkpoint, 1, MPI_LOGICAL, 0, mpi_comm_world, ierr )


    if( arrived_checkpoint .and. myid == 0 ) then
      write(STDERR,*)
      write(STDERR,'(A)') "INFO [ckprst__checkpoint]: Checkpoint arrived. Taking snapshot.. "
    end if

    if( arrived_checkpoint ) then

      call output__closefiles()

      call checkpoint_fname( fn_ckp )

      !! Open the checkpoint file with replace mode
      call std__getio( io_ckp )
      open( io_ckp, file=trim(fn_ckp), form='unformatted', action='write', status='replace' )
      write( io_ckp ) it

      call global__checkpoint( io_ckp )
      call kernel__checkpoint( io_ckp )
      call medium__checkpoint( io_ckp )
      call absorb__checkpoint( io_ckp )
      call report__checkpoint( io_ckp )
      call source__checkpoint( io_ckp )
      call output__checkpoint( io_ckp )

      call pwatch__off( "ckprst__checkpoint" )  !! to save pwatch, stop mesurement before checkpoint of pwatch itself
      call pwatch__checkpoint( io_ckp )

      close( io_ckp )

      call mpi_barrier(mpi_comm_world, ierr )

      if( myid == 0 ) then
        write(STDERR,'(A)') "INFO [ckprst__checkpoint]: Finished."
        write(STDERR,*) ""
      end if

      call mpi_finalize(ierr)
      stop
    else
      call pwatch__off( "ckprst__checkpoint" )
    end if


  end subroutine ckprst__checkpoint
  !! --------------------------------------------------------------------------------------------------------------------------- !!

  subroutine ckprst__setup( io_prm )

    integer, intent(in) :: io_prm
    integer :: t_rate, t_max

    !! ----

    call readini( io_prm, 'is_ckp', is_ckp, .false. )
    call readini( io_prm, 'ckpdir', ckpdir, trim(odir)//'/ckp' )
    call readini( io_prm, 'ckp_interval', ckp_interval, 100000000 )
    call readini( io_prm, 'ckp_time', ckp_time, 1e10 )

    if( .not. is_ckp ) return

    call system__call('mkdir -p '// trim(ckpdir) // '> /dev/null 2>&1' )

    !! initialize time
    call system_clock( t0, t_rate, t_max )

  end subroutine ckprst__setup

  !! --------------------------------------------------------------------------------------------------------------------------- !!
  subroutine ckprst__restart( it )

    integer, intent(out) :: it

    integer :: io_ckp
    character(256) :: fn_ckp
    integer :: ierr
    logical :: is_exist
    !! --

    if( .not. is_ckp ) then
      it = 1
      return
    end if


    call checkpoint_fname( fn_ckp )

    inquire( file=trim(fn_ckp), exist=is_exist )
    if( .not. is_exist ) then
      it = 1
      return
    end if


    call std__getio( io_ckp)
    open( io_ckp, file=trim(fn_ckp), form='unformatted', action='read', iostat=ierr, status='old' )

    !! if file does not exist
    if( ierr /= 0 ) then
      it = 1
    else
      read( io_ckp ) it
    end if

    if( it == 1 ) then

      return

    else if( it < 0 ) then

      if( myid == 0 ) then
        write(STDERR,*) 'INFO [ckprst__restart]: Program has already finished in previous execution'
        write(STDERR,*) 'INFO [ckprst__restart]: Terminate .. '
      end if
      call mpi_finalize( ierr )
      stop
    end if


    it = it + 1 ! advance one time step

    if( myid == 0 ) then
      write(STDERR,*)
      write(STDERR,'(A,I5)') "INFO [ckprst__restart]: restart from it = ", it
      write(STDERR,*)
    end if
    call global__restart( io_ckp )
    call kernel__restart( io_ckp )
    call medium__restart( io_ckp )
    call absorb__restart( io_ckp )
    call report__restart( io_ckp )
    call source__restart( io_ckp )
    call output__restart( io_ckp )
    call pwatch__restart( io_ckp )

    close( io_ckp )

    call mpi_barrier( mpi_comm_world, ierr )

  end subroutine ckprst__restart
  !! --------------------------------------------------------------------------------------------------------------------------- !!
  subroutine ckprst__filedelete()

    character(256) :: fn_ckp
    integer :: io
    integer :: ierr
    integer, parameter :: it_minus = -1

    if( .not. is_ckp ) return

    call pwatch__on('ckprst__filedelete')

    call checkpoint_fname( fn_ckp )

    ! replace file
    call std__getio(io)
    open(io, file=fn_ckp, form='unformatted', action='write', iostat=ierr, status='replace')
    write(io) it_minus
    close(io)

    call pwatch__off('ckprst__filedelete')

  end subroutine ckprst__filedelete
  !! --------------------------------------------------------------------------------------------------------------------------- !!
  subroutine checkpoint_fname( fn_ckp )
    character(*), intent(out) :: fn_ckp
    character(6) :: cid

    write(cid,  '(I6.6)') myid

    fn_ckp = trim(ckpdir) //  '/' // trim(title) // '.checkpoint.'// trim(cid)

  end subroutine checkpoint_fname
  !! --------------------------------------------------------------------------------------------------------------------------- !!

end module m_ckprst
!! ----------------------------------------------------------------------------------------------------------------------------- !!
