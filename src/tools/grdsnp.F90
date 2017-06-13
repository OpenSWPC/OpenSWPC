!! ----------------------------------------------------------------------------------------------------------------------------- !!
!>
!! Extract x-y-z data of velocity discontinuity from a specified grd dat
!!
!! @copyright
!!   Copyright 2013-2017 Takuto Maeda. All rights reserved. This project is released under the MIT license.
!<
!! --
#include "m_debug.h"
program grdsnp

  use m_std
  use m_getopt
  use m_bicubic
  use m_fdtool
  use m_geomap
  use m_readini
  use m_debug
  use netcdf


  integer :: i, j
  integer :: nx, ny
  real(SP)  :: dx, dy, dz, clon, clat, phi, xbeg, ybeg, zbeg
  real(SP), allocatable :: xc(:), yc(:), glon(:,:), glat(:,:)
  real(SP), allocatable ::lon(:), lat(:), grddep(:,:)
  real(SP) :: zgrd
  logical :: is_opt
  integer :: kgrd
  integer :: io
  integer :: nlon, nlat
  character(256) :: adum
  type(bicubic__data) :: bcd
  character(256) :: fn_prm, fn_grd
  integer :: ncid, ndim, nvar, xid, yid, zid
  character(80) :: xname, yname, zname
  !! ----

  call getopt( 'i', is_opt, fn_prm ); if( .not. is_opt ) call usage_exit
  call getopt( 'g', is_opt, fn_grd ); if( .not. is_opt ) call usage_exit

  call std__getio(io)
  open(io,file=trim(fn_prm),action='read', status='old')
  call readini( io, 'dx', dx, 0.5 )
  call readini( io, 'dy', dy, 0.5 )
  call readini( io, 'dz',      dz,      0.5        )
  call readini( io, 'nx',      nx,      256           )
  call readini( io, 'ny',      ny,      256           )
  call readini( io, 'clon',    clon,  139.7604        )
  call readini( io, 'clat',    clat,   35.7182        )
  call readini( io, 'phi',     phi,     0.0           )
  call readini( io, 'xbeg',    xbeg,   -nx/2*real(dx) )
  call readini( io, 'ybeg',    ybeg,   -ny/2*real(dy) )
  call readini( io, 'zbeg',    zbeg,   -ny/2*real(dy) )

  !!
  !! Grid locations in geographic coordinate
  !!
  allocate( xc(nx), yc(ny) )
  allocate( glon(nx,ny), glat(nx,ny) )
  do i=1, nx
    xc(i) = i2x( i, xbeg, real(dx) )
  end do
  do j=1, ny
    yc(j) = j2y( j, ybeg, real(dy) )
  end do

  do j=1, ny
    do i=1, nx
      call geomap__c2g( xc(i), yc(j), clon, clat, phi, glon(i,j), glat(i,j) )
    end do
  end do


  !! read netcdf data
  call assert( nf90_open( trim(fn_grd), NF90_NOWRITE, ncid ) == NF90_NOERR  )
  call assert( nf90_inquire( ncid, ndim, nvar ) ==NF90_NOERR )
  call assert( ndim == 2 )  !! assume 2D netcdf file
  nvar = nvar - 2 !! first two variables should be x(lon) and y(lat)

  !! size
  call assert( nf90_inquire_dimension( ncid, 1, len=nlon ) == NF90_NOERR )
  call assert( nf90_inquire_dimension( ncid, 2, len=nlat ) == NF90_NOERR )
  allocate( lon(nlon), lat(nlat) )
  allocate( grddep(nlon,nlat) )

  !! read
  call assert( nf90_inquire_variable( ncid, 1, xname ) == NF90_NOERR )
  call assert( nf90_inquire_variable( ncid, 2, yname ) == NF90_NOERR )
  call assert( nf90_inquire_variable( ncid, 3, zname ) == NF90_NOERR )
  call assert( nf90_inq_varid( ncid, xname, xid )      == NF90_NOERR )
  call assert( nf90_inq_varid( ncid, yname, yid )      == NF90_NOERR )
  call assert( nf90_inq_varid( ncid, zname, zid )      == NF90_NOERR )
  call assert( nf90_get_var( ncid, xid, lon )          == NF90_NOERR )
  call assert( nf90_get_var( ncid, yid, lat )          == NF90_NOERR )
  call assert( nf90_get_var( ncid, zid, grddep )       == NF90_NOERR )

  call assert( nf90_close( ncid ) == NF90_NOERR )


  grddep(:,:) = grddep(:,:) / 1000 ! m -> km

  call bicubic__init( bcd, nlon, nlat, lon(1), lat(1), lon(2)-lon(1), lat(2)-lat(1), grddep )

  do j=1, ny
    do i=1, nx
      call bicubic__interp( bcd, glon(i,j), glat(i,j), zgrd )

      kgrd = z2k( real(zgrd), zbeg, real(dz) ) !! boxel boundary : zc+dz/2

      if( kgrd == 1 ) then
        kgrd = kgrd + 1
      end if
      write(STDOUT,'(2I6,4F12.4,ES12.4,I12)') i, j, xc(i), yc(j), glon(i,j), glat(i,j), zgrd, kgrd
    end do
    write(STDOUT,*)
  end do

  call bicubic__terminate(bcd)

contains

  subroutine usage_exit
    write(STDERR,*) 'grdsnp.x -i input.inf  -g grd_file'
    stop
  end subroutine usage_exit
end program grdsnp
