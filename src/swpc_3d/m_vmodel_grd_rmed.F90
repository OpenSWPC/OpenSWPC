!! ----------------------------------------------------------------------------------------------------------------------------- !!
!>
!! User-routines for defining velocity/attenuation structure: GMT(netcdf) input
!!
!! @copyright
!!   Copyright 2013-2017 Takuto Maeda. All rights reserved. This project is released under the MIT license.
!<
!! ----
#include "m_debug.h"
#ifdef _NETCDF

module m_vmodel_grd_rmed

  use m_std
  use m_debug
  use m_global
  use m_bicubic
  use m_fdtool
  use m_geomap
  use m_system
  use m_readini
  use m_fdtool
  use m_rdrmed
  use mpi
  use netcdf

  implicit none
  private
  save

  public :: vmodel_grd_rmed

contains

  !! --------------------------------------------------------------------------------------------------------------------------- !!
  !>
  !! Define meidum velocity, density and attenuation
  !<
  !! ----
  subroutine vmodel_grd_rmed( io_prm, i0, i1, j0, j1, k0, k1, xc, yc, zc, vcut, rho, lam, mu, Qp, Qs, bd )

    !! -- Arguments
    integer,  intent(in)  :: io_prm
    integer,  intent(in)  :: i0, i1                         !< i-region
    integer,  intent(in)  :: j0, j1                         !< j-region
    integer,  intent(in)  :: k0, k1                         !< k-region
    real(SP), intent(in)  :: xc  ( i0:i1 )                  !< x-coordinate location
    real(SP), intent(in)  :: yc  ( j0:j1 )                  !< y-coordinate location
    real(SP), intent(in)  :: zc  ( k0:k1 )                  !< z-coordinate location
    real(SP), intent(in)  :: vcut                           !< cut-off velocity
    real(SP), intent(out) :: rho ( k0:k1, i0:i1, j0:j1 )    !< mass density [g/cm^3]
    real(SP), intent(out) :: lam ( k0:k1, i0:i1, j0:j1 )    !< Lame's parameter lambda [ (g/cm^3) * (km/s) ]
    real(SP), intent(out) :: mu  ( k0:k1, i0:i1, j0:j1 )    !< Lame's parameter mu     [ (g/cm^3) * (km/s) ]
    real(SP), intent(out) :: qp  ( k0:k1, i0:i1, j0:j1 )    !< P-wave attenuation
    real(SP), intent(out) :: qs  ( k0:k1, i0:i1, j0:j1 )    !< S-wave attenuation
    real(SP), intent(out) :: bd  ( i0:i1, j0:j1, 0:NBD )    !< Boundary depths

    !! --
    type(bicubic__data), allocatable :: bcd(:)
    character(256) :: fn_grdlst
    character(256) :: dir_grd
    logical :: is_ocean
    real(SP) :: rho0, vp0, vs0, qp0, qs0
    real(SP) :: zgrd
    integer :: ngrd
    integer :: ierr
    real(SP), allocatable :: lon(:), lat(:), grddep(:,:), grdbuf(:)
    real(SP), allocatable :: rho1(:), vp1(:), vs1(:), qp1(:), qs1(:)
    integer,  allocatable :: pid(:)
    real(SP) :: glon( i0:i1, j0:j1), glat(i0:i1,j0:j1) !< grid longitude,latitude
    character(256), allocatable :: fn_grd(:)
    integer :: iolst
    integer :: i, j, k, n, kk, l
    character(256) :: adum
    integer :: nlon, nlat
    integer :: node_grd
    real(SP) :: x_AB, x_AE, y_AB, y_AE ! absorber boundary location
    real(SP) :: xx, yy
    logical :: is_flatten
    integer :: idum
    integer :: ktopo
    integer, allocatable :: kgrd(:,:,:)
    character(256), allocatable :: fn_rmed(:)
    real(SP), allocatable :: xi(:,:,:,:)
    integer :: n_rmed
    integer, allocatable :: tbl_rmed(:)
    character(256), allocatable :: fn_rmed2(:)
    character(256) :: dir_rmed
    integer, allocatable :: reflyr(:)
    logical :: is_exist
    real(SP) :: rho2, vp2, vs2
    real(SP) :: vmin, vmax, dh, cc, rhomin
    logical  :: vmax_over, vmin_under, rhomin_under
    integer :: ncid, ndim, nvar, xid, yid, zid
    character(80) :: xname, yname, zname
    !! ----

    call readini( io_prm, 'fn_grdlst_rmed', fn_grdlst, '.' )
    call readini( io_prm, 'dir_grd',  dir_grd, '.' )
    call readini( io_prm, 'node_grd', node_grd, 0 )
    call readini( io_prm, 'is_ocean', is_ocean, .true. )
    call readini( io_prm, 'is_flatten', is_flatten, .false. )

    if( is_flatten ) is_ocean=.true.

    call readini( io_prm, 'dir_grd',  dir_grd, '.' )
    call readini( io_prm, 'dir_rmed',  dir_rmed, '.' )

    call readini( io_prm, 'rhomin', rhomin, 1.0 )

    vmin = vcut

    dh = sqrt(3.) / sqrt( 1./dx**2 + 1./dy**2 + 1./dz**2 )
    cc = 6. / 7. !! assume 4th order
    vmax = 0.95 * cc * dh / dt  ! 0.95 is a safety coefficient

    vmax_over  = .false.
    vmin_under = .false.
    rhomin_under = .false.


    !!
    !! first initialize whole medium by air/ocean
    !!
    do j=j0,j1
      do i=i0,i1

        !! filled by air
        do k=k0, k1

          vp0  = 0.0
          vs0  = 0.0
          rho0 = 0.001
          qp0  = 1.0
          qs0  = 1.0

          rho(k,i,j) = rho0
          lam(k,i,j) = rho0 * ( vp0 * vp0 - 2 * vs0 * vs0 )
          mu (k,i,j) = rho0 * vs0 * vs0
          qp (k,i,j) = qp0
          qs (k,i,j) = qs0

        end do

        !! filled by ocean, if required
        if( is_ocean ) then
          do k=k0, k1

            if( zc(k) < 0 ) cycle

            vp0  = 1.5
            vs0  = 0.0
            rho0 = 1.0
            qp0  = 1000000.0
            qs0  = 1000000.0

            rho(k,i,j) = rho0
            lam(k,i,j) = rho0 * ( vp0 * vp0 - 2 * vs0 * vs0 )
            mu (k,i,j) = rho0 *                   vs0 * vs0
            qp (k,i,j) = qp0
            qs (k,i,j) = qs0

          end do
        end if

      end do
    end do

    !!
    !! Grid locations in geographic coordinate
    !!

    !! boundary value
    x_AB = i2x( na + 1,   xbeg, real(dx) )
    x_AE = i2x( nx - na,  xbeg, real(dx) )
    y_AB = j2y( na + 1,   ybeg, real(dy) )
    y_AE = j2y( ny - na,  ybeg, real(dy) )

    do j=j0, j1
      do i=i0, i1

        xx = min( max( xc(i), x_AB), x_AE )
        yy = min( max( yc(j), y_AB), y_AE )
        call geomap__c2g( xx, yy, clon, clat, phi, glon(i,j), glat(i,j) )
      end do
    end do


    !!
    !! layer list file
    !!
    call std__getio( iolst )
    open( iolst, file=trim( fn_grdlst ), action='read', status='old' )
    call std__countline( iolst, ngrd, '#' )
    allocate( fn_grd(ngrd) )
    allocate( rho1(ngrd), vp1(ngrd), vs1(ngrd), qp1(ngrd), qs1(ngrd), pid(ngrd), fn_rmed(ngrd), reflyr(ngrd) )
    do n=1, ngrd
      read(iolst,*) fn_grd(n), rho1(n), vp1(n), vs1(n), qp1(n), qs1(n), pid(n), fn_rmed(n), reflyr(n)
      fn_grd(n) = trim(adjustl(dir_grd)) // '/' // trim(adjustl(fn_grd(n)))
      fn_rmed(n) = trim(adjustl(dir_rmed)) // '/' // trim(adjustl(fn_rmed(n)))
      call assert( 0<= reflyr(n) .and. reflyr(n) <= ngrd )
    end do
    !!
    !! random media
    !!
    allocate( tbl_rmed(ngrd), fn_rmed2(ngrd) )
    call independent_list( ngrd, fn_rmed, n_rmed, tbl_rmed, fn_rmed2 )

    allocate(xi(k0:k1,i0:i1,j0:j1,n_rmed) )
    do l=1, n_rmed
      inquire( file=trim(fn_rmed2(l)), exist=is_exist )
      if( is_exist ) then
        call rdrmed__3d( i0, i1, j0, j1, k0, k1, fn_rmed2(l), xi(k0:k1,i0:i1,j0:j1,l) )
      else
        xi(k0:k1,i0:i1,j0:j1,l) = 0.0
      end if
    end do


    !! cut-off velocity: filled by medium velocity of the deeper layer
    do n=ngrd-1, 1, -1
      if( (vp1(n) < vcut .or. vs1(n) < vcut) .and. ( vp1(n) > 0 .and. vs1(n) > 0 ) ) then
        vp1 (n) = vp1 (n+1)
        vs1 (n) = vs1 (n+1)
        rho1(n) = rho1(n+1)
        qp1 (n) = qp1 (n+1)
        qs1 (n) = qs1 (n+1)
      end if
    end do

    !! bicubic dataflow
    allocate( bcd(ngrd) )

    !!
    !! read file and interpolate
    !!
    allocate( kgrd( 0:ngrd, i0:i1, j0:j1 ) )
    kgrd(0,i0:i1,j0:j1) = k0 - 1

    ktopo = z2k( 0.0 - real(dz)/2, zbeg, real(dz) )

    do n=1, ngrd
      !! read grd data at id=node_grd

      if( myid == node_grd ) then
        !! netcdf file
        call assert( nf90_open( trim(fn_grd(n)), NF90_NOWRITE, ncid ) == NF90_NOERR  )
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

      end if

      !! Grd data size communication
      call mpi_bcast( nlon, 1, MPI_INTEGER, node_grd, MPI_COMM_WORLD, ierr )
      call mpi_bcast( nlat, 1, MPI_INTEGER, node_grd, MPI_COMM_WORLD, ierr )

      !! memory allocation
      allocate( grdbuf(nlon*nlat) )
      if( myid == node_grd ) then
        !! set up send buffer
        grdbuf = reshape( grddep, shape( grdbuf ) )
      else
        !! memory is already allocated at node_grd by reading file
        allocate( lon(nlon), lat(nlat), grddep(nlon,nlat) )
      end if

      !! send and receive grd data
      call mpi_bcast(lon, nlon, MPI_REAL, node_grd, MPI_COMM_WORLD, ierr )
      call mpi_bcast(lat, nlat, MPI_REAL, node_grd, MPI_COMM_WORLD, ierr )
      call mpi_bcast(grdbuf, nlon*nlat, MPI_REAL, node_grd, MPI_COMM_WORLD, ierr )

      !! reshape received buffer to 2D array
      if( myid /= node_grd ) grddep = reshape( grdbuf, shape( grddep ) )
      deallocate( grdbuf )


      grddep(:,:) = grddep(:,:) / 1000 ! m -> km

      call bicubic__init( bcd(n), nlon, nlat, lon(1), lat(1), lon(2)-lon(1), lat(2)-lat(1), grddep )

      do j=j0, j1
        do i=i0, i1

          call bicubic__interp( bcd(n), glon(i,j), glat(i,j), zgrd )

          if( n == 1 ) bd(i,j,0) = zgrd

          if( is_flatten ) zgrd = zgrd - bd(i,j,0)

          kgrd(n,i,j) = max( z2k( zgrd-real(dz)/2, zbeg, real(dz) ), kgrd(n-1,i,j) )


          !! seafloor correction: (sea column thickness)>=2
          !! the following condition is always false for is_flatten=.true.
          if( n==1 .and. zgrd > 0 ) kgrd(n,i,j) = max( ktopo + 2, kgrd(n,i,j) )

          ! memorize specified discontinuity depth (such as plate boundary)
          if( pid(n) > 0 ) then
            bd(i,j,pid(n)) = zgrd
          end if

        end do
      end do

      call bicubic__terminate( bcd(n) )
      deallocate( lon, lat, grddep )

    end do

    !!
    !! set medium
    !!
    do j=j0, j1
      do i=i0, i1
        do n=1, ngrd
          !! fills deeper structure
          do k=kgrd(n,i,j)+1, k1

            kk = k - kgrd(reflyr(n),i,j) + 1 !< relative depth index
            !! cyclic condition
            if( kk < k0 ) kk = kk + nz
            if( kk > k1 ) kk = kk - k1

            !! background velocity must not exceeds stability condition
            call assert( vp1(n) < vmax )
            call assert( vs1(n) < vmax )

            !! add random medium
            vp2 = vp1(n) * ( 1.0 + xi(kk,i,j,tbl_rmed(n) ) )
            vs2 = vs1(n) * ( 1.0 + xi(kk,i,j,tbl_rmed(n) ) )
            rho2 = rho1(n) * ( 1.0 + 0.8 * xi(kk,i,j,tbl_rmed(n) ) )

            call vcheck(vp2, vs2, rho2, xi(kk,i,j,tbl_rmed(n)), vmin, vmax, rhomin, vmin_under, vmax_over, rhomin_under)

            rho(k,i,j) = rho2
            lam(k,i,j) = rho2 * ( vp2 * vp2 - 2 * vs2 * vs2 )
            mu (k,i,j) = rho2 * vs2 * vs2
            qp (k,i,j) = qp1(n)
            qs (k,i,j) = qs1(n)
          end do

        end do

      end do
    end do

    !! notification for velocity torelance
    if( vmax_over  ) call info( 'Too high velocity due to random media was corrected. ')
    if( vmin_under ) call info( 'Too low  velocity due to random media was corrected. ')
    if( rhomin_under ) call info( 'Too low  density due to random media was corrected. ')



    !!
    !! terminate
    !!
    deallocate( fn_grd )
    deallocate( rho1, vp1, vs1, qp1, qs1 )
    deallocate( pid )
    deallocate( bcd, kgrd )
    deallocate( fn_rmed, fn_rmed2 )
    deallocate( reflyr )
    deallocate( tbl_rmed )
    deallocate( xi )

  end subroutine vmodel_grd_rmed
  !! --------------------------------------------------------------------------------------------------------------------------- !!

end module m_vmodel_grd_rmed
!! ----------------------------------------------------------------------------------------------------------------------------- !!

#else

module m_vmodel_grd_rmed

  !! This is a dummy module for non-netcdf mode

  use m_std
  use m_global
  use mpi

contains
  !! --------------------------------------------------------------------------------------------------------------------------- !!
  subroutine vmodel_grd_rmed( io, i0, i1, j0, j1, k0, k1, xc, yc, zc, vcut, rho, lam, mu, Qp, Qs, bd )

    !! -- Arguments
    integer,  intent(in)  :: io                             !< prm file io
    integer,  intent(in)  :: i0, i1                         !< i-region
    integer,  intent(in)  :: j0, j1                         !< j-region
    integer,  intent(in)  :: k0, k1                         !< k-region
    real(SP), intent(in)  :: xc  ( i0:i1 )                  !< x-coordinate location
    real(SP), intent(in)  :: yc  ( j0:j1 )                  !< y-coordinate location
    real(SP), intent(in)  :: zc  ( k0:k1 )                  !< z-coordinate location
    real(SP), intent(in)  :: vcut                           !< cut-off velocity
    real(SP), intent(out) :: rho ( k0:k1, i0:i1, j0:j1 )    !< mass density [g/cm^3]
    real(SP), intent(out) :: lam ( k0:k1, i0:i1, j0:j1 )    !< Lame's parameter lambda [ (g/cm^3) * (km/s) ]
    real(SP), intent(out) :: mu  ( k0:k1, i0:i1, j0:j1 )    !< Lame's parameter mu     [ (g/cm^3) * (km/s) ]
    real(SP), intent(out) :: qp  ( k0:k1, i0:i1, j0:j1 )    !< P-wave attenuation
    real(SP), intent(out) :: qs  ( k0:k1, i0:i1, j0:j1 )    !< S-wave attenuation
    real(SP), intent(out) :: bd  ( i0:i1, j0:j1, 0:NBD )    !< Boundary depths

    integer :: ierr
    !! --

    rho = 0.0
    lam = 0.0
    mu = 0.0
    qp = 0.0
    qs = 0.0
    bd = 0.0

    if( myid == 0 ) then
      call info('This program is not compiled with netcdf. ')
      call assert( .false. )
    end if

  end subroutine vmodel_grd_rmed
  !! --------------------------------------------------------------------------------------------------------------------------- !!

end module m_vmodel_grd_rmed

#endif
