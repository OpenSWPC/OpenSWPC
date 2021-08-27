!! ----------------------------------------------------------------------------------------------------------------------------- !!
!>
!! User-routines for defining velocity/attenuation structure: Layered medium input
!!
!! @copyright
!!   Copyright 2013-2021 Takuto Maeda. All rights reserved. This project is released under the MIT license.
!<
!! ----
#include "m_debug.h"
#ifdef _NETCDF

module m_vmodel_grd

  use m_std
  use m_debug
  use m_global
  use m_bicubic
  use m_fdtool
  use m_geomap
  use m_system
  use netcdf

  implicit none
  private
  save

  public :: vmodel_grd

contains

  !! --------------------------------------------------------------------------------------------------------------------------- !!
  !>
  !! Define meidum velocity, density and attenuation
  !<
  !! ----
  subroutine vmodel_grd( io_prm, i0, i1, k0, k1, xc, zc, vcut, rho, lam, mu, Qp, Qs, bd )

    !! -- Arguments
    integer,  intent(in)  :: io_prm
    integer,  intent(in)  :: i0, i1                  !< i-region
    integer,  intent(in)  :: k0, k1                  !< k-region
    real(SP), intent(in)  :: xc  ( i0:i1 )           !< x-coordinate location
    real(SP), intent(in)  :: zc  ( k0:k1 )           !< z-coordinate location
    real(SP), intent(in)  :: vcut                    !< cut-off velocity
    real(SP), intent(out) :: rho ( k0:k1, i0:i1 )    !< mass density [g/cm^3]
    real(SP), intent(out) :: lam ( k0:k1, i0:i1 )    !< Lame's parameter lambda [ (g/cm^3) * (km/s) ]
    real(SP), intent(out) :: mu  ( k0:k1, i0:i1 )    !< Lame's parameter mu     [ (g/cm^3) * (km/s) ]
    real(SP), intent(out) :: qp  ( k0:k1, i0:i1 )    !< P-wave attenuation
    real(SP), intent(out) :: qs  ( k0:k1, i0:i1 )    !< S-wave attenuation
    real(SP), intent(out) :: bd  ( i0:i1, 0:NBD )    !< Boundary depths
    !! --

    character(256) :: fn_grdlst, dir_grd
    logical :: is_ocean
    real(SP) :: rho0, vp0, vs0, qp0, qs0
    real(SP) :: zgrd
    integer, allocatable  :: kgrd(:,:)
    integer :: ngrd
    real(DP), allocatable :: lon(:), lat(:), grddep(:,:)
    real(DP) :: dlon, dlat
    real(SP), allocatable :: rho1(:), vp1(:), vs1(:), qp1(:), qs1(:)
    integer,  allocatable :: pid(:)
    real(SP) :: glon(i0:i1), glat(i0:i1) !< grid longitude,latitude
    character(256), allocatable :: fn_grd(:)
    integer :: iolst
    integer :: i, k, n
    integer :: nlon, nlat
    character(256) :: adum
    logical :: is_flatten
    type(bicubic__data), allocatable :: bcd(:)
    integer :: ktopo
    integer :: ncid, ndim, nvar, xid, yid, zid
    character(80) :: xname, yname, zname
    logical :: earth_flattening
    real(SP) :: Cv(k0:k1) ! velocity scaling coefficient for earth_flattening    
    !! ----

    call readini( io_prm, 'fn_grdlst', fn_grdlst, '' )
    call readini( io_prm, 'is_ocean', is_ocean, .true. )
    call readini( io_prm, 'topo_flatten', is_flatten, .false. )
    if( is_flatten ) is_ocean = .true.
    call readini( io_prm, 'dir_grd',  dir_grd, '.' )

   !! earth-flattening transform
    call readini( io_prm, 'earth_flattening', earth_flattening, .false. )
    if( earth_flattening ) then
      do k=k0, k1
        Cv(k) = exp( zc(k) / R_EARTH)
      end do
    else
      Cv(:) = 1.0
    end if

    !!
    !! first initialize with air/ocean
    !!
    do i=i0,i1

      !! filled by air
      do k=k0, k1

        vp0  = 0.0
        vs0  = 0.0
        rho0 = 0.001
        qp0  = 1.0
        qs0  = 1.0

        rho(k,i) = rho0
        lam(k,i) = rho0 * ( vp0 * vp0 - 2 * vs0 * vs0 )
        mu (k,i) = rho0 * vs0 * vs0
        qp (k,i) = qp0
        qs (k,i) = qs0

      end do

      !! filled by ocean, if required
      if( is_ocean ) then
        do k=k0, k1

          if( zc(k) < 0 ) cycle

          vp0  = Cv(k) * 1.5
          vs0  = 0.0
          rho0 = 1.0 * Cv(k)**(-5)
          qp0  = 1000000.0
          qs0  = 1000000.0

          rho(k,i) = rho0
          lam(k,i) = rho0 * ( vp0 * vp0 - 2 * vs0 * vs0 )
          mu (k,i) = rho0 *                   vs0 * vs0
          qp (k,i) = qp0
          qs (k,i) = qs0

        end do
      end if

    end do

    !!
    !! Grid locations in geographic coordinate
    !!
    do i=i0, i1
      call geomap__c2g( xc(i), 0.0, clon, clat, phi, glon(i), glat(i) )
    end do


    !!
    !! layer list file
    !!
    call std__getio( iolst )
    open( iolst, file=trim( fn_grdlst ), action='read', status='old' )
    call std__countline( iolst, ngrd, '#' )
    allocate( fn_grd(ngrd) )
    allocate( rho1(ngrd), vp1(ngrd), vs1(ngrd), qp1(ngrd), qs1(ngrd), pid(ngrd) )
    do n=1, ngrd
      read(iolst,*) fn_grd(n), rho1(n), vp1(n), vs1(n), qp1(n), qs1(n), pid(n)
      fn_grd(n) = trim(adjustl(dir_grd)) // '/' // trim(adjustl(fn_grd(n)))
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

    !! bicubic data
    allocate( bcd(ngrd) )

    allocate( kgrd( 0:ngrd, i0:i1 ) )
    kgrd(0,i0:i1) = k0-1



    ktopo = z2k( 0.0 - real(dz)/2, zbeg, real(dz) )
    !!
    !! read file and interpolate
    !!
    do n=1, ngrd

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


      grddep(:,:) = grddep(:,:) / 1000 ! m -> km

      dlon = (lon(nlon) - lon(1))/(nlon - 1)
      dlat = (lat(nlat) - lat(1))/(nlat - 1)      
      call bicubic__init( bcd(n), nlon, nlat, lon(1), lat(1), dlon, dlat, grddep )

      do i=i0, i1
        call bicubic__interp( bcd(n), glon(i), glat(i), zgrd )

        if( earth_flattening ) then
          zgrd = - R_EARTH * log( (R_EARTH - zgrd) / R_EARTH )
        end if        

        if( n == 1 ) bd(i,0) = zgrd
        if( is_flatten ) zgrd = zgrd - bd(i,0)

        kgrd(n,i) = max( z2k( zgrd-real(dz)/2, zbeg, real(dz) ), kgrd(n-1,i) )

        !! seafloor correction: (sea column thickness)>=2
        !! the following condition is always false for is_flatten=.true.
        if( n==1 .and. zgrd > 0 ) kgrd(n,i) = max( ktopo + 2, kgrd(n,i) )

        ! memorize specified discontinuity depth (such as plate boundary)
        if( pid(n) > 0 ) then
          bd(i,pid(n)) = zgrd
        end if

      end do
      call bicubic__terminate( bcd(n) )
      deallocate( lon, lat, grddep )
    end do

    do i=i0, i1
      do n=1,ngrd
        do k=kgrd(n,i)+1, k1
          rho(k,i) = rho1(n) * Cv(k)**(-5)
          lam(k,i) = rho1(n) * Cv(k)**(-5) * Cv(k)**2 * ( vp1(n) * vp1(n) - 2 * vs1(n) * vs1(n) )
          mu (k,i) = rho1(n) * Cv(k)**(-5) * Cv(k)**2 *                         vs1(n) * vs1(n)
          qp (k,i) = qp1(n)
          qs (k,i) = qs1(n)
        end do

      end do
    end do

    !!
    !! terminate
    !!
    deallocate( fn_grd )
    deallocate( rho1, vp1, vs1, qp1, qs1 )
    deallocate( pid )
    deallocate( kgrd )
    deallocate( bcd )
  end subroutine vmodel_grd

end module m_vmodel_grd

#else

module m_vmodel_grd

  !! This is a dummy module for non-netcdf mode
  use m_std
  use m_global
  use mpi

contains

  !! --------------------------------------------------------------------------------------------------------------------------- !!
  subroutine vmodel_grd( io_prm, i0, i1, k0, k1, xc, zc, vcut, rho, lam, mu, Qp, Qs, bd )

    !! -- Arguments
    integer,  intent(in)  :: io_prm
    integer,  intent(in)  :: i0, i1                  !< i-region
    integer,  intent(in)  :: k0, k1                  !< k-region
    real(SP), intent(in)  :: xc  ( i0:i1 )           !< x-coordinate location
    real(SP), intent(in)  :: zc  ( k0:k1 )           !< z-coordinate location
    real(SP), intent(in)  :: vcut                    !< cut-off velocity
    real(SP), intent(out) :: rho ( k0:k1, i0:i1 )    !< mass density [g/cm^3]
    real(SP), intent(out) :: lam ( k0:k1, i0:i1 )    !< Lame's parameter lambda [ (g/cm^3) * (km/s) ]
    real(SP), intent(out) :: mu  ( k0:k1, i0:i1 )    !< Lame's parameter mu     [ (g/cm^3) * (km/s) ]
    real(SP), intent(out) :: qp  ( k0:k1, i0:i1 )    !< P-wave attenuation
    real(SP), intent(out) :: qs  ( k0:k1, i0:i1 )    !< S-wave attenuation
    real(SP), intent(out) :: bd  ( i0:i1, 0:NBD )    !< Boundary depths

    integer :: ierr
    !! --

    if( myid == 0 ) then
      write(STDERR,'(A)') 'ERROR [vmodel_grd]: This program is not compiled with netcdf. Abort.'
      call mpi_finalize( ierr )
      stop
    end if

  end subroutine vmodel_grd
  !! --------------------------------------------------------------------------------------------------------------------------- !!

end module m_vmodel_grd
!! ----------------------------------------------------------------------------------------------------------------------------- !!

#endif
