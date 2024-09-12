!! ----------------------------------------------------------------------------------------------------------------------------- !!
!>
!!  SWPC: Seismic Wave Propagation Code
!!
!! @detail
!!   This software simulate seismic wave propagation, by solving equations of motion with constitutive equations of
!!   elastic/visco-elastic medium by finite difference method (FDM).
!!
!! @copyright
!!   Copyright 2013-2024 Takuto Maeda. All rights reserved. This project is released under the MIT license.
!!
!! @see
!!   - Noguchi et al.     (2016) GJI     doi:10.1093/gji/ggw074
!!   - Maeda et al.       (2014) EPSL    doi:10.1016/j.epsl.2014.04.037
!!   - Maeda et al.       (2013) BSSA    doi:10.1785/0120120118
!!   - Maeda and Furumura (2013) PAGEOPH doi:10.1007/s00024-011-0430-z
!!   - Noguchi et al.     (2013) PAGEOPH doi:10.1007/s00024-011-0412-1
!!   - Furumura et al.    (2008) PAGEOPH doi:10.1007/s00024-008-0318-8
!!   - Furumura ad Chen   (2005) PARCO   doi:10.1016/j.parco.2005.02.003
!<
!! ----------------------------------------------------------------------------------------------------------------------------- !!
#include "../shared/m_debug.h"
program SWPC_SH

  !! -- Dependency
  use m_std
  use m_debug
  use m_global
  use m_kernel
  use m_getopt
  use m_source
  use m_medium
  use m_report
  use m_pwatch
  use m_output
  use m_wav
  use m_absorb
  use m_readini
  use m_version
  use mpi

  !! -- Declarations
  implicit none
  !! --

  character(256) :: fn_prm
  integer :: it
  integer :: ierr
  logical :: is_opt
  logical :: stopwatch_mode
  integer :: io_prm, io_watch
  logical :: strict_mode
  !! ----

  !!
  !! Version 
  !!
  call getopt('v', is_opt);        if( is_opt ) call version__display('swpc_sh')
  call getopt('-version', is_opt); if( is_opt ) call version__display('swpc_sh')

  !!
  !! Launch MPI process
  !!
  call mpi_init( ierr )

  !!
  !! option processing
  !!
  call getopt('i', is_opt, fn_prm, './in/params.nml' )
  open( newunit=io_prm, file=trim(fn_prm), action='read', status='old' )

  call readini( io_prm, 'stopwatch_mode', stopwatch_mode, .true.  )

  call readini( io_prm, 'strict_mode', strict_mode, .false. )
  call readini__strict_mode( strict_mode )

  !!
  !! Read control parameters
  !!
  call global__setup( io_prm )

  !!
  !! stopwatch start
  !!
  call pwatch__setup( stopwatch_mode )
  
  !!
  !! set-up each module
  !!
  call global__setup2( )
  call medium__setup( io_prm )
  call mpi_barrier( mpi_comm_world, ierr )
  call kernel__setup( )
  call source__setup( io_prm )
  call absorb__setup( io_prm )
  call snap__setup( io_prm )
  call wav__setup(io_prm)
  call report__setup( io_prm )

  close( io_prm )


  !! mainloop
  do it = 1, nt

    call report__progress(it)

    call snap__write( it )
    call wav__store( it )

    call kernel__update_stress()

    call absorb__update_stress()
    call source__stressdrop(it)

    call global__comm_stress()

    call kernel__update_vel()
    call source__bodyforce(it)

    call absorb__update_vel()
    call global__comm_vel()


  end do

  call wav__write()
  call snap__closefiles()

  !! ending message
  call report__terminate()

  !!
  !! stopwatch report from 0-th node
  !!
  if( stopwatch_mode ) then
    if( myid == 0 ) then
      open( newunit=io_watch, file=trim( odir ) // '/' // trim( title ) // '.tim', action='write', status='unknown' )
    end if
    call pwatch__report( io_watch, 0 )

  end if

  !!
  !! Program termination
  !!
  call mpi_barrier( mpi_comm_world, ierr ) 
  call mpi_finalize( ierr )

  stop

end program SWPC_SH

!! ----------------------------------------------------------------------------------------------------------------------------- !!
