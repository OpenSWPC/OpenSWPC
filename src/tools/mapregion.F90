!! ----------------------------------------------------------------------------------------------------------------------------- !!
!>
!! Obtains projected area for FDM simulation from input parameter file
!!
!! @copyright
!!   Copyright 2013-2017 Takuto Maeda. All rights reserved. This project is released under the MIT license.
!<
!! ----------------------------------------------------------------------------------------------------------------------------- !!
program mapregion

  use m_std
  use m_getopt
  use m_readini
  use m_fdtool
  use m_geomap

  implicit none

  character(256) :: fn_prm, fn_out
  integer :: io_prm

  integer               :: nx, ny, nz                               !<  space grid number (global)
  integer               :: na                                       !<  absorber thickness
  real(SP)              :: dx, dy, dz                               !<  space grid width
  real(SP)              :: xbeg                               !<  global coordinate: x start / end
  real(SP)              :: ybeg                               !<  global coordinate: y start / end
  real(SP)              :: zbeg                               !<  global coordinate: z start / end

  !!
  !! MPI domain
  !!
  integer               :: nproc_x                                  !<  process numbers for x/i - direction
  integer               :: nproc_y                                  !<  process numbers for y/j - direction
  logical :: is_given
  integer :: ierr
  real(SP) :: clon, clat, phi
  integer :: nn
  real(SP), allocatable :: lon(:), lat(:)
  integer :: i, j, k
  real(SP) :: x, y
  integer :: io
  integer :: nm = 3

  real(SP) :: mem_n, mem_a

  call getopt( 'i', is_given, fn_prm )
  if( .not. is_given ) then
    write(STDERR,*)
    write(STDERR,'(A)') ' mapregion.x [-i input_prm] [-o output (default=stdout)] '
    write(STDERR,*)
    stop
  end if

  call getopt( 'o', is_given, fn_out, '' )

  call std__getio( io_prm )

  open( io_prm, file=trim(fn_prm), action='read', status='old', iostat=ierr )
  if( ierr /= 0 ) then
    write(STDERR,*) 'ERROR [mapregion]: Could not open file ' // trim( fn_prm )
    stop
  end if

  call readini( io_prm, 'nproc_x', nproc_x, 1        )
  call readini( io_prm, 'nproc_y', nproc_y, 1        )
  call readini( io_prm, 'nx',      nx,      256      )
  call readini( io_prm, 'ny',      ny,      256      )
  call readini( io_prm, 'nz',      nz,      256      )
  call readini( io_prm, 'dx',      dx,      0.5      )
  call readini( io_prm, 'dy',      dy,      0.5      )
  call readini( io_prm, 'dz',      dz,      0.5      )
  call readini( io_prm, 'na',      na,      0        )
  call readini( io_prm, 'xbeg',    xbeg,   -nx/2*dx  )
  call readini( io_prm, 'ybeg',    ybeg,   -ny/2*dy  )
  call readini( io_prm, 'zbeg',    zbeg,   -30*dz    )
  call readini( io_prm, 'clon',    clon,  139.7604   )
  call readini( io_prm, 'clat',    clat,   35.7182   )
  call readini( io_prm, 'phi',     phi,     0.0      )

  call memory_size( nproc_x, nproc_y, nx, ny, nz, nm, na, mem_a, mem_n )
  write(STDERR,*) "grid size ", nx, ny, nz
  write(STDERR,*) "kernel memory size ( whole ) ", mem_a, " [GB]"
  write(STDERR,*) "kernel memory size ( node  ) ", mem_n, " [GB]"

  nn = 2*(nx-1)+2*(ny-1) + 1
  allocate( lon(nn), lat(nn) )
  k = 1
  do i=1, nx
    j = 1
    x = xbeg + (i-1)*dx
    y = ybeg + (j-1)*dy
    call geomap__c2g( x, y, clon, clat, phi, lon(k), lat(k) )
    k = k + 1
  end do
  do j=2, ny
    i = nx
    x = xbeg + (i-1)*dx
    y = ybeg + (j-1)*dy
    call geomap__c2g( x, y, clon, clat, phi, lon(k), lat(k) )
    k = k + 1
  end do
  do i=nx-1,1,-1
    j = ny
    x = xbeg + (i-1)*dx
    y = ybeg + (j-1)*dy
    call geomap__c2g( x, y, clon, clat, phi, lon(k), lat(k) )
    k = k + 1
  end do
  do j=ny-1,1,-1
    i = 1
    x = xbeg + (i-1)*dx
    y = ybeg + (j-1)*dy
    call geomap__c2g( x, y, clon, clat, phi, lon(k), lat(k) )
    k = k + 1
  end do

  if( is_given ) then
    call std__getio( io )
    open( io, file=trim(fn_out), action='write', status='unknown' )
    do k=1, nn
      write(io,*) lon(k), lat(k)
    end do
    close( io )
  else
    do k=1, nn
      write(*,*) lon(k), lat(k)
    end do
  end if

end program mapregion
