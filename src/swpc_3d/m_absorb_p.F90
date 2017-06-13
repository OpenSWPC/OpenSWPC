!! ----------------------------------------------------------------------------------------------------------------------------- !!
!>
!! Absorbing Boundary Condition: ADE-CFS PML based on Zhang and Shen
!!
!! @par PML region definition
!!  +-----+--------------------------+-----+
!!  |     |                          |     |
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
!!
!! @copyright
!!   Copyright 2013-2017 Takuto Maeda. All rights reserved. This project is released under the MIT license.
!!
!<
!! ----
#include "m_debug.h"
module m_absorb_p

  !! -- Dependency
  use m_std
  use m_debug
  use m_global
  use m_fdtool
  use m_readini

  !! -- Declarations
  implicit none
  private
  save

  !! -- Public Procedures
  public :: absorb_p__setup
  public :: absorb_p__update_stress
  public :: absorb_p__update_vel
  public :: absorb_p__checkpoint
  public :: absorb_p__restart
  !! --

  !! Damping profiles
  real(SP), allocatable :: gxc(:,:), gxe(:,:) !< damping profile along x at center/edge of voxel
  real(SP), allocatable :: gyc(:,:), gye(:,:) !< damping profile along x at center/edge of voxel
  real(SP), allocatable :: gzc(:,:), gze(:,:) !< damping profile along x at center/edge of voxel

  !! ADE variables
  real(SP), allocatable :: axVx(:,:,:), ayVx(:,:,:), azVx(:,:,:)
  real(SP), allocatable :: axVy(:,:,:), ayVy(:,:,:), azVy(:,:,:)
  real(SP), allocatable :: axVz(:,:,:), ayVz(:,:,:), azVz(:,:,:)
  real(SP), allocatable :: axSxx(:,:,:), aySxy(:,:,:), azSxz(:,:,:)
  real(SP), allocatable :: axSxy(:,:,:), aySyy(:,:,:), azSyz(:,:,:)
  real(SP), allocatable :: axSxz(:,:,:), aySyz(:,:,:), azSzz(:,:,:)

  real(MP) :: r20x, r20y, r20z

  integer :: kbeg_min

contains

  !! --------------------------------------------------------------------------------------------------------------------------- !!
  !>
  !! Setup
  !! set PML sponge
  !<
  !! ----
  subroutine absorb_p__setup( io_prm )

    integer, intent(in) :: io_prm
    !! --
    integer  :: i, j, k
    real(SP) :: hx, hy, hz
    integer  :: idum
    !! ----

    !!
    !! derivative coefficient
    !!
    r20x = 1.0 / dx
    r20y = 1.0 / dy
    r20z = 1.0 / dz

    !!
    !! damping profile
    !!
    allocate( gxc(4,ibeg:iend), gxe(4,ibeg:iend) )
    allocate( gyc(4,jbeg:jend), gye(4,jbeg:jend) )
    allocate( gzc(4,kbeg:kend), gze(4,kbeg:kend) )

    hx = na * dx
    hy = na * dy
    hz = na * dz
    do i=ibeg, iend
      call damping_profile( xc(i),              hx, xbeg, xend, gxc(:,i) )
      call damping_profile( xc(i)+real(dx)/2.0, hx, xbeg, xend, gxe(:,i) )
    end do
    do j=jbeg, jend
      call damping_profile( yc(j),              hy, ybeg, yend, gyc(:,j) )
      call damping_profile( yc(j)+real(dy)/2.0, hy, ybeg, yend, gye(:,j) )
    end do
    do k=kbeg, kend
      call damping_profile( zc(k),              hz, zbeg, zend, gzc(:,k) )
      call damping_profile( zc(k)+real(dz)/2.0, hz, zbeg, zend, gze(:,k) )
    end do


    !!
    !! PML region definition
    !!

    kbeg_min = minval( kbeg_a(:,:) )

    !! memory allocation
    allocate(  axVx( kbeg_min:kend, ibeg:iend, jbeg:jend ) )
    allocate(  ayVx( kbeg_min:kend, ibeg:iend, jbeg:jend ) )
    allocate(  azVx( kbeg_min:kend, ibeg:iend, jbeg:jend ) )
    allocate(  axVy( kbeg_min:kend, ibeg:iend, jbeg:jend ) )
    allocate(  ayVy( kbeg_min:kend, ibeg:iend, jbeg:jend ) )
    allocate(  azVy( kbeg_min:kend, ibeg:iend, jbeg:jend ) )
    allocate(  axVz( kbeg_min:kend, ibeg:iend, jbeg:jend ) )
    allocate(  ayVz( kbeg_min:kend, ibeg:iend, jbeg:jend ) )
    allocate(  azVz( kbeg_min:kend, ibeg:iend, jbeg:jend ) )
    allocate( axSxx( kbeg_min:kend, ibeg:iend, jbeg:jend ) )
    allocate( aySxy( kbeg_min:kend, ibeg:iend, jbeg:jend ) )
    allocate( azSxz( kbeg_min:kend, ibeg:iend, jbeg:jend ) )
    allocate( axSxy( kbeg_min:kend, ibeg:iend, jbeg:jend ) )
    allocate( aySyy( kbeg_min:kend, ibeg:iend, jbeg:jend ) )
    allocate( azSyz( kbeg_min:kend, ibeg:iend, jbeg:jend ) )
    allocate( axSxz( kbeg_min:kend, ibeg:iend, jbeg:jend ) )
    allocate( aySyz( kbeg_min:kend, ibeg:iend, jbeg:jend ) )
    allocate( azSzz( kbeg_min:kend, ibeg:iend, jbeg:jend ) )

    axVx ( kbeg_min:kend, ibeg:iend, jbeg:jend ) = 0.0
    ayVx ( kbeg_min:kend, ibeg:iend, jbeg:jend ) = 0.0
    azVx ( kbeg_min:kend, ibeg:iend, jbeg:jend ) = 0.0
    axVy ( kbeg_min:kend, ibeg:iend, jbeg:jend ) = 0.0
    ayVy ( kbeg_min:kend, ibeg:iend, jbeg:jend ) = 0.0
    azVy ( kbeg_min:kend, ibeg:iend, jbeg:jend ) = 0.0
    axVz ( kbeg_min:kend, ibeg:iend, jbeg:jend ) = 0.0
    ayVz ( kbeg_min:kend, ibeg:iend, jbeg:jend ) = 0.0
    azVz ( kbeg_min:kend, ibeg:iend, jbeg:jend ) = 0.0
    axSxx( kbeg_min:kend, ibeg:iend, jbeg:jend ) = 0.0
    aySxy( kbeg_min:kend, ibeg:iend, jbeg:jend ) = 0.0
    azSxz( kbeg_min:kend, ibeg:iend, jbeg:jend ) = 0.0
    axSxy( kbeg_min:kend, ibeg:iend, jbeg:jend ) = 0.0
    aySyy( kbeg_min:kend, ibeg:iend, jbeg:jend ) = 0.0
    azSyz( kbeg_min:kend, ibeg:iend, jbeg:jend ) = 0.0
    axSxz( kbeg_min:kend, ibeg:iend, jbeg:jend ) = 0.0
    aySyz( kbeg_min:kend, ibeg:iend, jbeg:jend ) = 0.0
    azSzz( kbeg_min:kend, ibeg:iend, jbeg:jend ) = 0.0


    idum = io_prm

  end subroutine absorb_p__setup
  !! --------------------------------------------------------------------------------------------------------------------------- !!

  !! --------------------------------------------------------------------------------------------------------------------------- !!
  !>
  !! Update velocity component in PML layer
  !<
  !! ----
  subroutine absorb_p__update_vel

    integer :: i, j, k
    real(SP) :: gxc0(4), gxe0(4), gyc0(4), gye0(4), gzc0(4), gze0(4)
    real(MP) :: dxSxx(kbeg:kend), dySxy(kbeg:kend), dzSxz(kbeg:kend)
    real(MP) :: dxSxy(kbeg:kend), dySyy(kbeg:kend), dzSyz(kbeg:kend)
    real(MP) :: dxSxz(kbeg:kend), dySyz(kbeg:kend), dzSzz(kbeg:kend)

    !!
    !! Horizontal zero-derivative boundary (for plane wave mode)
    !!
    if( pw_mode ) then
      if( idx == 0 ) then
        !$omp parallel private(j,k)
        !$omp do schedule(dynamic)
        do j=jbeg, jend
          do k=kbeg_a(1,j), kend
            Sxx(k,0,j) = 2 * Sxx(k,1,j) - Sxx(k,2,j)
            Syy(k,0,j) = 2 * Syy(k,1,j) - Syy(k,2,j)
            Szz(k,0,j) = 2 * Szz(k,1,j) - Szz(k,2,j)
            Syz(k,0,j) = 2 * Syz(k,1,j) - Syz(k,2,j)
            Sxz(k,0,j) = 2 * Sxz(k,1,j) - Sxz(k,2,j)
            Sxy(k,0,j) = 2 * Sxy(k,1,j) - Sxy(k,2,j)
          end do
        end do
        !$omp end do nowait
        !$omp end parallel
      end if

      if( idx == nproc_x -1 ) then
        !$omp parallel private(j,k)
        !$omp do schedule(dynamic)
        do j=jbeg, jend
          do k=kbeg_a(nx,j), kend
            Sxx(k,nx+1,j) = 2 * Sxx(k,nx,j) - Sxx(k,nx-1,j)
            Syy(k,nx+1,j) = 2 * Syy(k,nx,j) - Syy(k,nx-1,j)
            Szz(k,nx+1,j) = 2 * Szz(k,nx,j) - Szz(k,nx-1,j)
            Syz(k,nx+1,j) = 2 * Syz(k,nx,j) - Syz(k,nx-1,j)
            Sxz(k,nx+1,j) = 2 * Sxz(k,nx,j) - Sxz(k,nx-1,j)
            Sxy(k,nx+1,j) = 2 * Sxy(k,nx,j) - Sxy(k,nx-1,j)
          end do
        end do
        !$omp end do nowait
        !$omp end parallel
      end if

      if( idy == 0 ) then
        !$omp parallel private(i,k)
        !$omp do schedule(dynamic)
        do i=ibeg, iend
          do k=kbeg_a(i,1), kend
            Sxx(k,i,0) = 2 * Sxx(k,i,1) - Sxx(k,i,2)
            Syy(k,i,0) = 2 * Syy(k,i,1) - Syy(k,i,2)
            Szz(k,i,0) = 2 * Szz(k,i,1) - Szz(k,i,2)
            Syz(k,i,0) = 2 * Syz(k,i,1) - Syz(k,i,2)
            Sxz(k,i,0) = 2 * Sxz(k,i,1) - Sxz(k,i,2)
            Sxy(k,i,0) = 2 * Sxy(k,i,1) - Sxy(k,i,2)
          end do
        end do
        !$omp end do nowait
        !$omp end parallel
      end if

      if( idy == nproc_y -1 ) then
        !$omp parallel private(i,k)
        !$omp do schedule(dynamic)
        do i=ibeg, iend
          do k=kbeg_a(i,ny), kend
            Sxx(k,i,ny+1) = 2 * Sxx(k,i,ny) - Sxx(k,i,ny-1)
            Syy(k,i,ny+1) = 2 * Syy(k,i,ny) - Syy(k,i,ny-1)
            Szz(k,i,ny+1) = 2 * Szz(k,i,ny) - Szz(k,i,ny-1)
            Syz(k,i,ny+1) = 2 * Syz(k,i,ny) - Syz(k,i,ny-1)
            Sxz(k,i,ny+1) = 2 * Sxz(k,i,ny) - Sxz(k,i,ny-1)
            Sxy(k,i,ny+1) = 2 * Sxy(k,i,ny) - Sxy(k,i,ny-1)
          end do
        end do
        !$omp end do nowait
        !$omp end parallel
      end if
      !$omp barrier
    end if

    !!
    !! time-marching
    !!

    !$omp parallel &
    !$omp private( dxSxx, dySyy, dzSzz, dySyz, dzSyz, dxSxz, dzSxz, dxSxy ,dySxy ) &
    !$omp private( gxc0, gxe0, gyc0, gye0, gzc0, gze0 ) &
    !$omp private( i, j, k )
    !$omp do &
    !$omp schedule(dynamic)
    do j=jbeg, jend

      gyc0(1:4) = gyc(1:4,j)
      gye0(1:4) = gye(1:4,j)

      do i=ibeg, iend

        gxc0(1:4) = gxc(1:4,i)
        gxe0(1:4) = gxe(1:4,i)

        !!
        !! Derivatives
        !!
        do k=kbeg_a(i,j), kend

          dxSxx(k) = (  Sxx(k  ,i+1,j  ) - Sxx(k  ,i  ,j  )  ) * r20x
          dySyy(k) = (  Syy(k  ,i  ,j+1) - Syy(k  ,i  ,j  )  ) * r20y
          dzSzz(k) = (  Szz(k+1,i  ,j  ) - Szz(k  ,i  ,j  )  ) * r20z
          dySyz(k) = (  Syz(k  ,i  ,j  ) - Syz(k  ,i  ,j-1)  ) * r20y
          dzSyz(k) = (  Syz(k  ,i  ,j  ) - Syz(k-1,i  ,j  )  ) * r20z
          dxSxz(k) = (  Sxz(k  ,i  ,j  ) - Sxz(k  ,i-1,j  )  ) * r20x
          dzSxz(k) = (  Sxz(k  ,i  ,j  ) - Sxz(k-1,i  ,j  )  ) * r20z
          dxSxy(k) = (  Sxy(k  ,i  ,j  ) - Sxy(k  ,i-1,j  )  ) * r20x
          dySxy(k) = (  Sxy(k  ,i  ,j  ) - Sxy(k  ,i  ,j-1)  ) * r20y

        end do


        !!
        !! update velocity
        !!
        do k=kbeg_a(i,j), kend

          gzc0(1:4) = gzc(1:4,k)
          gze0(1:4) = gze(1:4,k)

          !!
          !! Velocity Updates
          !!
          Vx(k,i,j) = Vx(k,i,j) &
              + bx(k,i,j) * ( gxe0(1) * dxSxx(k)     + gyc0(1) * dySxy(k)     + gzc0(1) * dzSxz(k)       &
              + gxe0(2) * axSxx(k,i,j) + gyc0(2) * aySxy(k,i,j) + gzc0(2) * azSxz(k,i,j)  ) * dt

          Vy(k,i,j) = Vy(k,i,j) &
              + by(k,i,j) * ( gxc0(1) * dxSxy(k)     + gye0(1) * dySyy(k)     + gzc0(1) * dzSyz(k)       &
              + gxc0(2) * axSxy(k,i,j) + gye0(2) * aySyy(k,i,j) + gzc0(2) * azSyz(k,i,j)  ) * dt

          Vz(k,i,j) = Vz(k,i,j) &
              + bz(k,i,j) * ( gxc0(1) * dxSxz(k)     + gyc0(1) * dySyz(k)     + gze0(1) * dzSzz(k)       &
              + gxc0(2) * axSxz(k,i,j) + gyc0(2) * aySyz(k,i,j) + gze0(2) * azSzz(k,i,j)  ) * dt


          !!
          !! ADE updates
          !!

          axSxx(k,i,j) = gxe0(3) * axSxx(k,i,j) + gxe0(4) * dxSxx(k) * dt
          aySxy(k,i,j) = gyc0(3) * aySxy(k,i,j) + gyc0(4) * dySxy(k) * dt
          azSxz(k,i,j) = gzc0(3) * azSxz(k,i,j) + gzc0(4) * dzSxz(k) * dt

          axSxy(k,i,j) = gxc0(3) * axSxy(k,i,j) + gxc0(4) * dxSxy(k) * dt
          aySyy(k,i,j) = gye0(3) * aySyy(k,i,j) + gye0(4) * dySyy(k) * dt
          azSyz(k,i,j) = gzc0(3) * azSyz(k,i,j) + gzc0(4) * dzSyz(k) * dt

          axSxz(k,i,j) = gxc0(3) * axSxz(k,i,j) + gxc0(4) * dxSxz(k) * dt
          aySyz(k,i,j) = gyc0(3) * aySyz(k,i,j) + gyc0(4) * dySyz(k) * dt
          azSzz(k,i,j) = gze0(3) * azSzz(k,i,j) + gze0(4) * dzSzz(k) * dt

        end do
      end do
    end do
    !$omp end do nowait
    !$omp end parallel

    !$omp barrier


  end subroutine absorb_p__update_vel
  !! --------------------------------------------------------------------------------------------------------------------------- !!

  !! --------------------------------------------------------------------------------------------------------------------------- !!
  subroutine absorb_p__update_stress

    integer :: i, j, k
    real(SP) :: lam2mu_R, lam_R
    real(SP) :: dxVx_ade, dyVy_ade, dzVz_ade
    real(SP) :: gxc0(4), gxe0(4), gyc0(4), gye0(4), gzc0(4), gze0(4)
    real(MP) :: dxVx(kbeg_min:kend),  dyVx(kbeg_min:kend),  dzVx(kbeg_min:kend)
    real(MP) :: dxVy(kbeg_min:kend),  dyVy(kbeg_min:kend),  dzVy(kbeg_min:kend)
    real(MP) :: dxVz(kbeg_min:kend),  dyVz(kbeg_min:kend),  dzVz(kbeg_min:kend)

    !! ----

    !!
    !! Horizontal zero-derivative boundary (for plane wave mode)
    !!
    if( pw_mode ) then
      if( idx == 0 ) then
        !$omp parallel private(j,k)
        !$omp do schedule(dynamic)
        do j=jbeg, jend
          do k=kbeg_a(1,j), kend
            Vx(k,0,j) = 2* Vx(k,1,j)-Vx(k,2,j)
            Vy(k,0,j) = 2* Vy(k,1,j)-Vy(k,2,j)
            Vz(k,0,j) = 2* Vz(k,1,j)-Vz(k,2,j)
          end do
        end do
        !$omp end do nowait
        !$omp end parallel
      end if

      if( idx == nproc_x -1 ) then
        !$omp parallel private(j,k)
        !$omp do schedule(dynamic)
        do j=jbeg, jend
          do k=kbeg_a(nx,j), kend
            Vx(k,nx+1,j) = 2 * Vx(k,nx,j) - Vx(k,nx-1,j)
            Vy(k,nx+1,j) = 2 * Vy(k,nx,j) - Vy(k,nx-1,j)
            Vz(k,nx+1,j) = 2 * Vz(k,nx,j) - Vz(k,nx-1,j)
          end do
        end do
        !$omp end do nowait
        !$omp end parallel
      end if

      if( idy == 0 ) then
        !$omp parallel private(i,k)
        !$omp do schedule(dynamic)
        do i=ibeg, iend
          do k=kbeg_a(i,1), kend
            Vx(k,i,0) = 2 * Vx(k,i,1) - Vx(k,i,2)
            Vy(k,i,0) = 2 * Vy(k,i,1) - Vy(k,i,2)
            Vz(k,i,0) = 2 * Vz(k,i,1) - Vz(k,i,2)
          end do
        end do
        !$omp end do nowait
        !$omp end parallel
      end if

      if( idy == nproc_y -1 ) then
        !$omp parallel private(i,k)
        !$omp do schedule(dynamic)
        do i=ibeg, iend
          do k=kbeg_a(i,ny), kend
            Vx(k,i,ny+1) = 2 * Vx(k,i,ny) - Vx(k,i,ny-1)
            Vy(k,i,ny+1) = 2 * Vy(k,i,ny) - Vy(k,i,ny-1)
            Vz(k,i,ny+1) = 2 * Vz(k,i,ny) - Vz(k,i,ny-1)
          end do
        end do
        !$omp end do nowait
        !$omp end parallel
      end if

      !$omp barrier
    end if


    !!
    !! Time-marching
    !!

    !$omp parallel &
    !$omp private( gxc0, gxe0, gyc0, gye0, gzc0, gze0 ) &
    !$omp private( dxVx, dxVy, dxVz, dyVx, dyVy, dyVz, dzVx, dzVy, dzVz ) &
    !$omp private( lam2mu_R, lam_R ) &
    !$omp private( dxVx_ade, dyVy_ade, dzVz_ade ) &
    !$omp private( i, j, k )
    !$omp do &
    !$omp schedule(dynamic)
    do j=jbeg, jend

      gyc0(1:4) = gyc(1:4,j)
      gye0(1:4) = gye(1:4,j)

      do i=ibeg, iend

        gxc0(1:4) = gxc(1:4,i)
        gxe0(1:4) = gxe(1:4,i)

        !!
        !! Derivatives
        !!
        do k=kbeg_a(i,j), kend

          dxVx(k) = (  Vx(k  ,i  ,j  ) - Vx(k  ,i-1,j  )  ) * r20x
          dxVy(k) = (  Vy(k  ,i+1,j  ) - Vy(k  ,i  ,j  )  ) * r20x
          dxVz(k) = (  Vz(k  ,i+1,j  ) - Vz(k  ,i  ,j  )  ) * r20x
          dyVx(k) = (  Vx(k  ,i  ,j+1) - Vx(k  ,i  ,j  )  ) * r20y
          dyVy(k) = (  Vy(k  ,i  ,j  ) - Vy(k  ,i  ,j-1)  ) * r20y
          dyVz(k) = (  Vz(k  ,i  ,j+1) - Vz(k  ,i  ,j  )  ) * r20y
          dzVx(k) = (  Vx(k+1,i  ,j  ) - Vx(k  ,i  ,j  )  ) * r20z
          dzVy(k) = (  Vy(k+1,i  ,j  ) - Vy(k  ,i  ,j  )  ) * r20z
          dzVz(k) = (  Vz(k  ,i  ,j  ) - Vz(k-1,i  ,j  )  ) * r20z

        end do

        !!
        !! Update Normal Stress
        !!
        do k=kbeg_a(i,j), kend

          gzc0(1:4) = gzc(1:4,k)

          !!
          !! update stress components
          !!
          lam2mu_R = ( lam(k,i,j) + 2 * mu(k,i,j) )
          lam_R    = lam2mu_R - 2*mu(k,i,j)


          !!
          !! Normal Stress
          !!
          dxVx_ade = gxc0(1) * dxVx(k) + gxc0(2) * axVx(k,i,j)
          dyVy_ade = gyc0(1) * dyVy(k) + gyc0(2) * ayVy(k,i,j)
          dzVz_ade = gzc0(1) * dzVz(k) + gzc0(2) * azVz(k,i,j)
          Sxx(k,i,j) = Sxx(k,i,j) + ( lam2mu_R * dxVx_ade + lam_R * ( dyVy_ade + dzVz_ade ) ) * dt
          Syy(k,i,j) = Syy(k,i,j) + ( lam2mu_R * dyVy_ade + lam_R * ( dxVx_ade + dzVz_ade ) ) * dt
          Szz(k,i,j) = Szz(k,i,j) + ( lam2mu_R * dzVz_ade + lam_R * ( dxVx_ade + dyVy_ade ) ) * dt


          !!
          !! ADE
          !!
          axVx(k,i,j) = gxc0(3) * axVx(k,i,j) + gxc0(4) * dxVx(k) * dt
          ayVy(k,i,j) = gyc0(3) * ayVy(k,i,j) + gyc0(4) * dyVy(k) * dt
          azVz(k,i,j) = gzc0(3) * azVz(k,i,j) + gzc0(4) * dzVz(k) * dt

        end do


        !!
        !! Update Shear Stress
        !!
        do k=kbeg_a(i,j), kend

          gze0(1:4) = gze(1:4,k)

          Syz(k,i,j) = Syz(k,i,j) + muyz(k,i,j) * ( gye0(1) * dyVz(k)     + gze0(1) * dzVy(k) &
              + gye0(2) * ayVz(k,i,j) + gze0(2) * azVy(k,i,j) ) * dt

          Sxz(k,i,j) = Sxz(k,i,j) + muxz(k,i,j) * ( gxe0(1) * dxVz(k)     + gze0(1) * dzVx(k) &
              + gxe0(2) * axVz(k,i,j) + gze0(2) * azVx(k,i,j) ) * dt

          Sxy(k,i,j) = Sxy(k,i,j) + muxy(k,i,j) * ( gxe0(1) * dxVy(k)     + gye0(1) * dyVx(k) &
              + gxe0(2) * axVy(k,i,j) + gye0(2) * ayVx(k,i,j) ) * dt

          ayVx(k,i,j) = gye0(3) * ayVx(k,i,j) + gye0(4) * dyVx(k) * dt
          azVx(k,i,j) = gze0(3) * azVx(k,i,j) + gze0(4) * dzVx(k) * dt
          axVy(k,i,j) = gxe0(3) * axVy(k,i,j) + gxe0(4) * dxVy(k) * dt
          azVy(k,i,j) = gze0(3) * azVy(k,i,j) + gze0(4) * dzVy(k) * dt
          axVz(k,i,j) = gxe0(3) * axVz(k,i,j) + gxe0(4) * dxVz(k) * dt
          ayVz(k,i,j) = gye0(3) * ayVz(k,i,j) + gye0(4) * dyVz(k) * dt


        end do
      end do
    end do
    !$omp end do nowait
    !$omp end parallel

    !$omp barrier

  end subroutine absorb_p__update_stress
  !! --------------------------------------------------------------------------------------------------------------------------- !!


  !! --------------------------------------------------------------------------------------------------------------------------- !!
  !>
  !! checkpoint data export
  !<
  !! --
  subroutine absorb_p__checkpoint( io )

    integer, intent(in) :: io
    integer :: j
    !! --

    write(io) r20x, r20y, r20z
    write(io) kbeg_min

    write(io) gxc(:,ibeg:iend), gxe(:,ibeg:iend)
    write(io) gyc(:,jbeg:jend), gye(:,jbeg:jend)
    write(io) gzc(:,kbeg:kend), gze(:,kbeg:kend)
    deallocate( gxc, gxe, gyc, gye, gzc, gze )

    do j=jbeg,jend;  write(io)  axVx(kbeg_min:kend,ibeg:iend,j); end do;  deallocate( axVx )
    do j=jbeg,jend;  write(io)  ayVx(kbeg_min:kend,ibeg:iend,j); end do;  deallocate( ayVx )
    do j=jbeg,jend;  write(io)  azVx(kbeg_min:kend,ibeg:iend,j); end do;  deallocate( azVx )
    do j=jbeg,jend;  write(io)  axVy(kbeg_min:kend,ibeg:iend,j); end do;  deallocate( axVy )
    do j=jbeg,jend;  write(io)  ayVy(kbeg_min:kend,ibeg:iend,j); end do;  deallocate( ayVy )
    do j=jbeg,jend;  write(io)  azVy(kbeg_min:kend,ibeg:iend,j); end do;  deallocate( azVy )
    do j=jbeg,jend;  write(io)  axVz(kbeg_min:kend,ibeg:iend,j); end do;  deallocate( axVz )
    do j=jbeg,jend;  write(io)  ayVz(kbeg_min:kend,ibeg:iend,j); end do;  deallocate( ayVz )
    do j=jbeg,jend;  write(io)  azVz(kbeg_min:kend,ibeg:iend,j); end do;  deallocate( azVz )

    do j=jbeg,jend;  write(io) axSxx(kbeg_min:kend,ibeg:iend,j); end do;  deallocate( axSxx )
    do j=jbeg,jend;  write(io) aySxy(kbeg_min:kend,ibeg:iend,j); end do;  deallocate( aySxy )
    do j=jbeg,jend;  write(io) azSxz(kbeg_min:kend,ibeg:iend,j); end do;  deallocate( azSxz )
    do j=jbeg,jend;  write(io) axSxy(kbeg_min:kend,ibeg:iend,j); end do;  deallocate( axSxy )
    do j=jbeg,jend;  write(io) aySyy(kbeg_min:kend,ibeg:iend,j); end do;  deallocate( aySyy )
    do j=jbeg,jend;  write(io) azSyz(kbeg_min:kend,ibeg:iend,j); end do;  deallocate( azSyz )
    do j=jbeg,jend;  write(io) axSxz(kbeg_min:kend,ibeg:iend,j); end do;  deallocate( axSxz )
    do j=jbeg,jend;  write(io) aySyz(kbeg_min:kend,ibeg:iend,j); end do;  deallocate( aySyz )
    do j=jbeg,jend;  write(io) azSzz(kbeg_min:kend,ibeg:iend,j); end do;  deallocate( azSzz )

  end subroutine absorb_p__checkpoint
  !! --------------------------------------------------------------------------------------------------------------------------- !!
  
  !! --------------------------------------------------------------------------------------------------------------------------- !!
  subroutine absorb_p__restart(io)

    integer, intent(in) :: io
    integer :: j
    
    read(io) r20x, r20y, r20z
    read(io) kbeg_min

    !! memory allocation
    allocate( gxc( 4,ibeg:iend ), gxe( 4,ibeg:iend ) )
    allocate( gyc( 4,jbeg:jend ), gye( 4,jbeg:jend ) )
    allocate( gzc( 4,kbeg:kend ), gze( 4,kbeg:kend ) )

    allocate(  axVx(kbeg_min:kend,ibeg:iend,jbeg:jend) )
    allocate(  ayVx(kbeg_min:kend,ibeg:iend,jbeg:jend) )
    allocate(  azVx(kbeg_min:kend,ibeg:iend,jbeg:jend) )
    allocate(  axVy(kbeg_min:kend,ibeg:iend,jbeg:jend) )
    allocate(  ayVy(kbeg_min:kend,ibeg:iend,jbeg:jend) )
    allocate(  azVy(kbeg_min:kend,ibeg:iend,jbeg:jend) )
    allocate(  axVz(kbeg_min:kend,ibeg:iend,jbeg:jend) )
    allocate(  ayVz(kbeg_min:kend,ibeg:iend,jbeg:jend) )
    allocate(  azVz(kbeg_min:kend,ibeg:iend,jbeg:jend) )
    allocate( axSxx(kbeg_min:kend,ibeg:iend,jbeg:jend) )
    allocate( aySxy(kbeg_min:kend,ibeg:iend,jbeg:jend) )
    allocate( azSxz(kbeg_min:kend,ibeg:iend,jbeg:jend) )
    allocate( axSxy(kbeg_min:kend,ibeg:iend,jbeg:jend) )
    allocate( aySyy(kbeg_min:kend,ibeg:iend,jbeg:jend) )
    allocate( azSyz(kbeg_min:kend,ibeg:iend,jbeg:jend) )
    allocate( axSxz(kbeg_min:kend,ibeg:iend,jbeg:jend) )
    allocate( aySyz(kbeg_min:kend,ibeg:iend,jbeg:jend) )
    allocate( azSzz(kbeg_min:kend,ibeg:iend,jbeg:jend) )

    read(io) gxc( :,ibeg:iend ), gxe( :,ibeg:iend )
    read(io) gyc( :,jbeg:jend ), gye( :,jbeg:jend )
    read(io) gzc( :,kbeg:kend ), gze( :,kbeg:kend )


    do j=jbeg,jend;  read(io)  axVx(kbeg_min:kend,ibeg:iend,j); end do;
    do j=jbeg,jend;  read(io)  ayVx(kbeg_min:kend,ibeg:iend,j); end do;
    do j=jbeg,jend;  read(io)  azVx(kbeg_min:kend,ibeg:iend,j); end do;
    do j=jbeg,jend;  read(io)  axVy(kbeg_min:kend,ibeg:iend,j); end do;
    do j=jbeg,jend;  read(io)  ayVy(kbeg_min:kend,ibeg:iend,j); end do;
    do j=jbeg,jend;  read(io)  azVy(kbeg_min:kend,ibeg:iend,j); end do;
    do j=jbeg,jend;  read(io)  axVz(kbeg_min:kend,ibeg:iend,j); end do;
    do j=jbeg,jend;  read(io)  ayVz(kbeg_min:kend,ibeg:iend,j); end do;
    do j=jbeg,jend;  read(io)  azVz(kbeg_min:kend,ibeg:iend,j); end do;
      
    do j=jbeg,jend;  read(io) axSxx(kbeg_min:kend,ibeg:iend,j); end do;
    do j=jbeg,jend;  read(io) aySxy(kbeg_min:kend,ibeg:iend,j); end do;
    do j=jbeg,jend;  read(io) azSxz(kbeg_min:kend,ibeg:iend,j); end do;
    do j=jbeg,jend;  read(io) axSxy(kbeg_min:kend,ibeg:iend,j); end do;
    do j=jbeg,jend;  read(io) aySyy(kbeg_min:kend,ibeg:iend,j); end do;
    do j=jbeg,jend;  read(io) azSyz(kbeg_min:kend,ibeg:iend,j); end do;
    do j=jbeg,jend;  read(io) axSxz(kbeg_min:kend,ibeg:iend,j); end do;
    do j=jbeg,jend;  read(io) aySyz(kbeg_min:kend,ibeg:iend,j); end do;
    do j=jbeg,jend;  read(io) azSzz(kbeg_min:kend,ibeg:iend,j); end do;

  end subroutine absorb_p__restart
  !! --------------------------------------------------------------------------------------------------------------------------- !!
  !! --------------------------------------------------------------------------------------------------------------------------- !!
  !>
  !! ADE-CFS PML damping factor according to Zhao and Shen
  !!
  subroutine damping_profile( x, H, xbeg, xend, g )
    
    real(SP), intent(in) :: x   !< cartesian coordinate location
    real(SP), intent(in) :: H   !< absorption layer thickness
    real(SP), intent(in) :: xbeg
    real(SP), intent(in) :: xend
    real(SP), intent(out) :: g(4) !< damping prof

    real(SP) :: R0 !! reflection coefficient
    real(SP) :: d0, a0, b0
    integer, parameter :: pd = 1
    integer, parameter :: pa = 1
    integer, parameter :: pb = 2
    real(SP), parameter :: cp = 6.0 !! assumed P-wave velocity
    real :: d, a, b, xx


    R0 = 10**( - ( log10( real(na) ) - 1 ) / log10( 2.0 )  - 3.0 )
    d0 = - ( 1.0 / (2.0*H) ) * ( pd +1 ) * cp * log( R0 )
    b0 = 6.0
    a0 = min( PI * fcut, d0/b0 *0.02 )

    if( x <= xbeg + H ) then
      xx = ( xbeg + H ) - x
    else if ( x >= xend - H ) then
      xx = x - ( xend - H )
    else
      xx = 0.0 !! no absorption
    end if

    d = d0 * abs( xx / H )**pd
    !a = a0 * ( 1.0 - abs( xx / H )**pa )
    b = 1.0 + ( b0 - 1.0 ) * abs( xx / H )**pb
    a = 0.02 * d / b

    g(1) = (  ( 1.0 + ( dt / 2.0 ) * a ) / b     ) / ( 1.0 + ( dt / 2.0 ) * ( a + d / b ) )
    g(2) = ( -1.0 / b                            ) / ( 1.0 + ( dt / 2.0 ) * ( a + d / b ) )
    g(3) = (  1.0 - ( dt / 2.0 ) * ( a + d / b ) ) / ( 1.0 + ( dt / 2.0 ) * ( a + d / b ) )
    g(4) = (  d / b                              ) / ( 1.0 + ( dt / 2.0 ) * ( a + d / b ) )

  end subroutine damping_profile
  !! --------------------------------------------------------------------------------------------------------------------------- !!
end module m_absorb_p
!! ----------------------------------------------------------------------------------------------------------------------------- !!
