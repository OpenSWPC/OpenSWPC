!! ----------------------------------------------------------------------------------------------------------------------------- !!
!>
!! Computation kernel for FDM numerical simulation
!!
!! @copyright
!!   Copyright 2013-2016 Takuto Maeda. All rights reserved. This project is released under the MIT license.
!<
!! ----
#include "m_debug.h"
module m_kernel

  !! -- Dependency
  use m_std
  use m_debug
  use m_global
  use m_pwatch
  use m_medium

  !! -- Declarations
  implicit none
  private
  save

  !! -- Public procedures
  public :: kernel__setup
  public :: kernel__update_vel
  public :: kernel__update_stress
  public :: kernel__vmax
  public :: kernel__checkpoint
  public :: kernel__restart

  !! -- Parameters
  real(MP), parameter   :: C20 = 1.0_MP
  real(MP), parameter   :: C40 = 9.0_MP /  8.0_MP
  real(MP), parameter   :: C41 = 1.0_MP / 24.0_MP

  !! --
  real(MP)              :: r40x, r40y, r40z
  real(MP)              :: r41x, r41y, r41z
  real(MP)              :: r20x, r20y, r20z
  real(SP), allocatable :: c1(:), c2(:), d1(:)
  real(SP) :: d2

contains

  !! --------------------------------------------------------------------------------------------------------------------------- !!
  !>
  !! Setup
  !!
  !! @detail
  !! this routine MUST BE CALLED AFTER medium__setup since it uses viscoelastic function
  !! use ts(:) array, dx, dy, dz from global module
  !<
  !! ----
  subroutine kernel__setup()

    integer :: m

    !! ----

    call pwatch__on("kernel__setup")
    call assert( medium__initialized() )

    r40x = C40 / dx
    r40y = C40 / dy
    r40z = C40 / dz
    r41x = C41 / dx
    r41y = C41 / dy
    r41z = C41 / dz
    r20x = C20 / dx
    r20y = C20 / dy
    r20z = C20 / dz

    call memory_allocate()

    !!
    !! initialize
    !!
    Vx  (       kbeg_m:kend_m, ibeg_m:iend_m, jbeg_m:jend_m ) = 0.0_MP
    Vy  (       kbeg_m:kend_m, ibeg_m:iend_m, jbeg_m:jend_m ) = 0.0_MP
    Vz  (       kbeg_m:kend_m, ibeg_m:iend_m, jbeg_m:jend_m ) = 0.0_MP
    Sxx (       kbeg_m:kend_m, ibeg_m:iend_m, jbeg_m:jend_m ) = 0.0_MP
    Syy (       kbeg_m:kend_m, ibeg_m:iend_m, jbeg_m:jend_m ) = 0.0_MP
    Szz (       kbeg_m:kend_m, ibeg_m:iend_m, jbeg_m:jend_m ) = 0.0_MP
    Syz (       kbeg_m:kend_m, ibeg_m:iend_m, jbeg_m:jend_m ) = 0.0_MP
    Sxz (       kbeg_m:kend_m, ibeg_m:iend_m, jbeg_m:jend_m ) = 0.0_MP
    Sxy (       kbeg_m:kend_m, ibeg_m:iend_m, jbeg_m:jend_m ) = 0.0_MP

    if( nm > 0 ) then
       Rxx ( kbeg_k:kend_k, ibeg_k:iend_k, jbeg_k:jend_k, 1:nm ) = 0.0
       Ryy ( kbeg_k:kend_k, ibeg_k:iend_k, jbeg_k:jend_k, 1:nm ) = 0.0
       Rzz ( kbeg_k:kend_k, ibeg_k:iend_k, jbeg_k:jend_k, 1:nm ) = 0.0
       Ryz ( kbeg_k:kend_k, ibeg_k:iend_k, jbeg_k:jend_k, 1:nm ) = 0.0
       Rxz ( kbeg_k:kend_k, ibeg_k:iend_k, jbeg_k:jend_k, 1:nm ) = 0.0
       Rxy ( kbeg_k:kend_k, ibeg_k:iend_k, jbeg_k:jend_k, 1:nm ) = 0.0

       do m=1, nm
          c1(m) = ( 2 * ts(m) - dt ) / ( 2 * ts(m) + dt )
          c2(m) = ( 2              ) / ( 2 * ts(m) + dt ) / nm
       end do
       d2 = sum( dt / ( 2*ts(:) - dt ) ) / nm
       do m=1, nm
          d1(m) = 2 * ts(m) / ( 2*ts(m) - dt )
       end do

    end if

    call pwatch__off("kernel__setup")

  end subroutine kernel__setup
  !! --------------------------------------------------------------------------------------------------------------------------- !!

  !! --------------------------------------------------------------------------------------------------------------------------- !!
  !>
  !! Update vel for one time step
  !<
  !! ----
  subroutine kernel__update_vel()

    integer :: i, j, k
    real(MP) :: d3Sx3(kbeg_k:kend_k), d3Sy3(kbeg_k:kend_k), d3Sz3(kbeg_k:kend_k)
    !! ----

    call pwatch__on("kernel__update_vel")


    !$omp parallel &
    !$omp private( d3Sx3, d3Sy3, d3Sz3 ) &
    !$omp private( i, j, k )
    !$omp do &
#ifdef _ES
    !$omp schedule(static)
#else
    !$omp schedule(static,1)
#endif
    do j=jbeg_k, jend_k

       do i=ibeg_k, iend_k

          !!
          !! derivateives
          !!

          !! stress derivatives
          do k=kbeg_k, kend_k
             d3Sx3(k) = (  Sxx(k  ,i+1,j  ) - Sxx(k  ,i  ,j  )  ) * r40x  -  (  Sxx(k  ,i+2,j  ) - Sxx(k  ,i-1,j  )  ) * r41x &
                      + (  Sxy(k  ,i  ,j  ) - Sxy(k  ,i  ,j-1)  ) * r40y  -  (  Sxy(k  ,i  ,j+1) - Sxy(k  ,i  ,j-2)  ) * r41y &
                      + (  Sxz(k  ,i  ,j  ) - Sxz(k-1,i  ,j  )  ) * r40z  -  (  Sxz(k+1,i  ,j  ) - Sxz(k-2,i  ,j  )  ) * r41z

             d3Sy3(k) = (  Sxy(k  ,i  ,j  ) - Sxy(k  ,i-1,j  )  ) * r40x  -  (  Sxy(k  ,i+1,j  ) - Sxy(k  ,i-2,j  )  ) * r41x &
                      + (  Syy(k  ,i  ,j+1) - Syy(k  ,i  ,j  )  ) * r40y  -  (  Syy(k  ,i  ,j+2) - Syy(k  ,i  ,j-1)  ) * r41y &
                      + (  Syz(k  ,i  ,j  ) - Syz(k-1,i  ,j  )  ) * r40z  -  (  Syz(k+1,i  ,j  ) - Syz(k-2,i  ,j  )  ) * r41z

             d3Sz3(k) = (  Sxz(k  ,i  ,j  ) - Sxz(k  ,i-1,j  )  ) * r40x  -  (  Sxz(k  ,i+1,j  ) - Sxz(k  ,i-2,j  )  ) * r41x &
                      + (  Syz(k  ,i  ,j  ) - Syz(k  ,i  ,j-1)  ) * r40y  -  (  Syz(k  ,i  ,j+1) - Syz(k  ,i  ,j-2)  ) * r41y &
                      + (  Szz(k+1,i  ,j  ) - Szz(k  ,i  ,j  )  ) * r40z  -  (  Szz(k+2,i  ,j  ) - Szz(k-1,i  ,j  )  ) * r41z
          end do

          !! overwrite around free surface
#ifdef _ES
!CDIR NOVECTOR
#endif
          do k=kfs_top(i,j), kfs_bot(i,j)
             d3Sx3(k) = (  Sxx(k  ,i+1,j  ) - Sxx(k  ,i  ,j  )  ) * r20x  &
                      + (  Sxy(k  ,i  ,j  ) - Sxy(k  ,i  ,j-1)  ) * r20y  &
                      + (  Sxz(k  ,i  ,j  ) - Sxz(k-1,i  ,j  )  ) * r20z

             d3Sy3(k) = (  Sxy(k  ,i  ,j  ) - Sxy(k  ,i-1,j  )  ) * r20x  &
                      + (  Syy(k  ,i  ,j+1) - Syy(k  ,i  ,j  )  ) * r20y  &
                      + (  Syz(k  ,i  ,j  ) - Syz(k-1,i  ,j  )  ) * r20z

             d3Sz3(k) = (  Sxz(k  ,i  ,j  ) - Sxz(k  ,i-1,j  )  ) * r20x  &
                      + (  Syz(k  ,i  ,j  ) - Syz(k  ,i  ,j-1)  ) * r20y  &
                      + (  Szz(k+1,i  ,j  ) - Szz(k  ,i  ,j  )  ) * r20z
          end do

          !! overwrite around seafloor
#ifdef _ES
!CDIR NOVECTOR
#endif
          do k=kob_top(i,j), kob_bot(i,j)
             d3Sx3(k) = (  Sxx(k  ,i+1,j  ) - Sxx(k  ,i  ,j  )  ) * r20x  &
                      + (  Sxy(k  ,i  ,j  ) - Sxy(k  ,i  ,j-1)  ) * r20y  &
                      + (  Sxz(k  ,i  ,j  ) - Sxz(k-1,i  ,j  )  ) * r20z

             d3Sy3(k) = (  Sxy(k  ,i  ,j  ) - Sxy(k  ,i-1,j  )  ) * r20x  &
                      + (  Syy(k  ,i  ,j+1) - Syy(k  ,i  ,j  )  ) * r20y  &
                      + (  Syz(k  ,i  ,j  ) - Syz(k-1,i  ,j  )  ) * r20z

             d3Sz3(k) = (  Sxz(k  ,i  ,j  ) - Sxz(k  ,i-1,j  )  ) * r20x  &
                      + (  Syz(k  ,i  ,j  ) - Syz(k  ,i  ,j-1)  ) * r20y  &
                      + (  Szz(k+1,i  ,j  ) - Szz(k  ,i  ,j  )  ) * r20z
          end do

          !!
          !! update velocity
          !!
          do k=kbeg_k, kend_k

             Vx(k,i,j) = Vx(k,i,j) + bx(k,i,j) * d3Sx3(k) * dt
             Vy(k,i,j) = Vy(k,i,j) + by(k,i,j) * d3Sy3(k) * dt
             Vz(k,i,j) = Vz(k,i,j) + bz(k,i,j) * d3Sz3(k) * dt

          end do
       end do
    end do
    !$omp end do nowait
    !$omp end parallel

    !$omp barrier

    call pwatch__off("kernel__update_vel")


  end subroutine kernel__update_vel
  !! --------------------------------------------------------------------------------------------------------------------------- !!


  !! --------------------------------------------------------------------------------------------------------------------------- !!
  !>
  !! Update stress for one time step
  !!
  !<
  !! ----
  subroutine kernel__update_stress()

    integer :: i, j, k, m, ierr
    real(SP) :: mu2, lam2mu
    real(SP) :: taup1, taus1, taup_plus1, taus_plus1
    real(SP) :: d3v3, dxVx_dyVy, dxVx_dzVz, dyVy_dzVz
    real(SP) :: Rxx_n, Ryy_n, Rzz_n, Ryz_n, Rxz_n, Rxy_n
    real(MP) :: dxVx(kbeg_k:kend_k),  dyVy(kbeg_k:kend_k), dzVz(kbeg_k:kend_k)
    real(MP) :: dxVy_dyVx(kbeg_k:kend_k), dxVz_dzVx(kbeg_k:kend_k), dyVz_dzVy(kbeg_k:kend_k)
    !! ----

    call pwatch__on("kernel__update_stress")

    !$omp parallel  &
    !$omp private( dxVx, dyVy, dzVz ) &
    !$omp private( mu2, lam2mu ) &
    !$omp private( taup1, taus1, taup_plus1, taus_plus1 ) &
    !$omp private( d3v3, dyVy_dzVz, dxVx_dzVz, dxVx_dyVy ) &
    !$omp private( Rxx_n, Ryy_n, Rzz_n ) &
    !$omp private( i, j, k, m )
    !$omp do &
#ifdef _ES
    !$omp schedule(static)
#else
    !$omp schedule(static,1)
#endif
    do j=jbeg_k, jend_k
       do i=ibeg_k, iend_k

          !! Derivatives
          !!

          do k=kbeg_k, kend_k

             dxVx(k) = (  Vx(k  ,i  ,j  ) - Vx(k  ,i-1,j  )  ) * r40x  -  (  Vx(k  ,i+1,j  ) - Vx(k  ,i-2,j  )  ) * r41x
             dyVy(k) = (  Vy(k  ,i  ,j  ) - Vy(k  ,i  ,j-1)  ) * r40y  -  (  Vy(k  ,i  ,j+1) - Vy(k  ,i  ,j-2)  ) * r41y
             dzVz(k) = (  Vz(k  ,i  ,j  ) - Vz(k-1,i  ,j  )  ) * r40z  -  (  Vz(k+1,i  ,j  ) - Vz(k-2,i  ,j  )  ) * r41z

          end do

          !! overwrite around free surface
#ifdef _ES
!CDIR NOVECTOR
#endif
          do k=kfs_top(i,j), kfs_bot(i,j)

             dxVx(k) = (  Vx(k  ,i  ,j  ) - Vx(k  ,i-1,j  )  ) * r20x
             dyVy(k) = (  Vy(k  ,i  ,j  ) - Vy(k  ,i  ,j-1)  ) * r20y
             dzVz(k) = (  Vz(k  ,i  ,j  ) - Vz(k-1,i  ,j  )  ) * r20z

          end do

          !! overwrite around seafloor
#ifdef _ES
!CDIR NOVECTOR
#endif
          do k=kob_top(i,j), kob_bot(i,j)

             dxVx(k) = (  Vx(k  ,i  ,j  ) - Vx(k  ,i-1,j  )  ) * r20x
             dyVy(k) = (  Vy(k  ,i  ,j  ) - Vy(k  ,i  ,j-1)  ) * r20y
             dzVz(k) = (  Vz(k  ,i  ,j  ) - Vz(k-1,i  ,j  )  ) * r20z

          end do


          !!
          !! update memory variables and stress tensors: normal stress components
          !!

#ifdef _FX
!ocl unroll('full')
#endif
          do k=kbeg_k, kend_k

             !!
             !! medium copy
             !!
             mu2    = 2*mu (k,i,j)
             lam2mu = lam(k,i,j) + mu2

             taup1 = taup(k,i,j)
             taus1 = taus(k,i,j)


             !!
             !! update memory variables
             !!

             !! working variables for combinations of velocity derivatives
             d3v3      = dxVx(k) + dyVy(k) + dzVz(k)
             dyVy_dzVz = dyVy(k) + dzVz(k)
             dxVx_dzVz = dxVx(k) + dzVz(k)
             dxVx_dyVy = dxVx(k) + dyVy(k)


             Rxx_n = 0.0
             Ryy_n = 0.0
             Rzz_n = 0.0
             do m=1, nm
                Rxx(k,i,j,m) = c1(m) * Rxx(k,i,j,m) - c2(m) * ( lam2mu * taup1 * d3v3 - mu2 * taus1 * dyVy_dzVz ) * dt
                Ryy(k,i,j,m) = c1(m) * Ryy(k,i,j,m) - c2(m) * ( lam2mu * taup1 * d3v3 - mu2 * taus1 * dxVx_dzVz ) * dt
                Rzz(k,i,j,m) = c1(m) * Rzz(k,i,j,m) - c2(m) * ( lam2mu * taup1 * d3v3 - mu2 * taus1 * dxVx_dyVy ) * dt
                Rxx_n = Rxx_n + d1(m) * Rxx(k,i,j,m)
                Ryy_n = Ryy_n + d1(m) * Ryy(k,i,j,m)
                Rzz_n = Rzz_n + d1(m) * Rzz(k,i,j,m)
             end do

             !!
             !! update stress components
             !!
             taup_plus1 = 1 + taup1 * ( 1 + d2 )
             taus_plus1 = 1 + taus1 * ( 1 + d2 )

             Sxx (k,i,j) = Sxx (k,i,j) + ( lam2mu * taup_plus1 * d3v3 - mu2 * taus_plus1 * dyVy_dzVz + Rxx_n ) * dt
             Syy (k,i,j) = Syy (k,i,j) + ( lam2mu * taup_plus1 * d3v3 - mu2 * taus_plus1 * dxVx_dzVz + Ryy_n ) * dt
             Szz (k,i,j) = Szz (k,i,j) + ( lam2mu * taup_plus1 * d3v3 - mu2 * taus_plus1 * dxVx_dyVy + Rzz_n ) * dt


          end do

       end do
    end do
    !$omp end do nowait
    !$omp end parallel


    !$omp parallel &
    !$omp private( dxVy_dyVx, dxVz_dzVx, dyVz_dzVy ) &
    !$omp private( mu2, taus1, taus_plus1 ) &
    !$omp private( Ryz_n, Rxz_n, Rxy_n ) &
    !$omp private( i, j, k, m )
    !$omp do  &
#ifdef _ES
    !$omp schedule(static)
#else
    !$omp schedule(static,1)
#endif
    do j=jbeg_k, jend_k

       do i=ibeg_k, iend_k

          !!
          !! Derivatives
          !!
          do k=kbeg_k, kend_k

             dxVy_dyVx(k) = (  Vy(k  ,i+1,j  ) - Vy(k  ,i  ,j  )  ) * r40x  -  (  Vy(k  ,i+2,j  ) - Vy(k  ,i-1,j  )  ) * r41x &
                          + (  Vx(k  ,i  ,j+1) - Vx(k  ,i  ,j  )  ) * r40y  -  (  Vx(k  ,i  ,j+2) - Vx(k  ,i  ,j-1)  ) * r41y
             dxVz_dzVx(k) = (  Vz(k  ,i+1,j  ) - Vz(k  ,i  ,j  )  ) * r40x  -  (  Vz(k  ,i+2,j  ) - Vz(k  ,i-1,j  )  ) * r41x &
                          + (  Vx(k+1,i  ,j  ) - Vx(k  ,i  ,j  )  ) * r40z  -  (  Vx(k+2,i  ,j  ) - Vx(k-1,i  ,j  )  ) * r41z
             dyVz_dzVy(k) = (  Vz(k  ,i  ,j+1) - Vz(k  ,i  ,j  )  ) * r40y  -  (  Vz(k  ,i  ,j+2) - Vz(k  ,i  ,j-1)  ) * r41y &
                          + (  Vy(k+1,i  ,j  ) - Vy(k  ,i  ,j  )  ) * r40z  -  (  Vy(k+2,i  ,j  ) - Vy(k-1,i  ,j  )  ) * r41z
          end do

          !! overwrite around free surface
#ifdef _ES
!CDIR NOVECTOR
#endif
          do k=kfs_top(i,j), kfs_bot(i,j)

             dxVy_dyVx(k) = (  Vy(k  ,i+1,j  ) - Vy(k  ,i  ,j  )  ) * r20x  &
                          + (  Vx(k  ,i  ,j+1) - Vx(k  ,i  ,j  )  ) * r20y
             dxVz_dzVx(k) = (  Vz(k  ,i+1,j  ) - Vz(k  ,i  ,j  )  ) * r20x  &
                          + (  Vx(k+1,i  ,j  ) - Vx(k  ,i  ,j  )  ) * r20z
             dyVz_dzVy(k) = (  Vz(k  ,i  ,j+1) - Vz(k  ,i  ,j  )  ) * r20y  &
                          + (  Vy(k+1,i  ,j  ) - Vy(k  ,i  ,j  )  ) * r20z

          end do

          !! overwrite around seafloor
#ifdef _ES
!CDIR NOVECTOR
#endif
          do k=kob_top(i,j), kob_bot(i,j)

             dxVy_dyVx(k) = (  Vy(k  ,i+1,j  ) - Vy(k  ,i  ,j  )  ) * r20x  &
                          + (  Vx(k  ,i  ,j+1) - Vx(k  ,i  ,j  )  ) * r20y
             dxVz_dzVx(k) = (  Vz(k  ,i+1,j  ) - Vz(k  ,i  ,j  )  ) * r20x  &
                          + (  Vx(k+1,i  ,j  ) - Vx(k  ,i  ,j  )  ) * r20z
             dyVz_dzVy(k) = (  Vz(k  ,i  ,j+1) - Vz(k  ,i  ,j  )  ) * r20y  &
                          + (  Vy(k+1,i  ,j  ) - Vy(k  ,i  ,j  )  ) * r20z

          end do


          !!
          !! update memory variables and stress tensors: shear stress components
          !!

#ifdef _FX
!ocl unroll('full')
#endif
          do k=kbeg_k, kend_k

             !!
             !! medium copy
             !!
             mu2    = 2*mu (k,i,j)
             taus1 = taus(k,i,j)


             !!
             !! update memory variables
             !!

             Ryz_n = 0.0
             Rxz_n = 0.0
             Rxy_n = 0.0
             do m=1, nm
                Ryz(k,i,j,m) = c1(m) * Ryz(k,i,j,m) - c2(m) * muyz(k,i,j) * taus1 * dyVz_dzVy(k) * dt
                Rxz(k,i,j,m) = c1(m) * Rxz(k,i,j,m) - c2(m) * muxz(k,i,j) * taus1 * dxVz_dzVx(k) * dt
                Rxy(k,i,j,m) = c1(m) * Rxy(k,i,j,m) - c2(m) * muxy(k,i,j) * taus1 * dxVy_dyVx(k) * dt
                Ryz_n = Ryz_n + d1(m) * Ryz(k,i,j,m)
                Rxz_n = Rxz_n + d1(m) * Rxz(k,i,j,m)
                Rxy_n = Rxy_n + d1(m) * Rxy(k,i,j,m)
             end do

             !!
             !! update stress components
             !!
             taus_plus1 = 1 + taus1 * ( 1 + d2 )

             Syz (k,i,j) = Syz (k,i,j) + ( muyz(k,i,j) * taus_plus1 * dyVz_dzVy(k) + Ryz_n ) * dt
             Sxz (k,i,j) = Sxz (k,i,j) + ( muxz(k,i,j) * taus_plus1 * dxVz_dzVx(k) + Rxz_n ) * dt
             Sxy (k,i,j) = Sxy (k,i,j) + ( muxy(k,i,j) * taus_plus1 * dxVy_dyVx(k) + Rxy_n ) * dt

          end do
       end do
    end do
    !$omp end do nowait
    !$omp end parallel

    !$omp barrier


    call pwatch__off("kernel__update_stress")

  end subroutine kernel__update_stress
  !! --------------------------------------------------------------------------------------------------------------------------- !!


  !! --------------------------------------------------------------------------------------------------------------------------- !!
  !>
  !! maximum value for terminal output
  !<
  !! --------------------------------------------------------------------------------------------------------------------------- !!
  subroutine kernel__vmax( xmax, ymax, zmax )
    real(SP), intent(out) :: xmax, ymax, zmax
    integer :: i, j
    integer, parameter :: margin = 5
    xmax = 0.0
    ymax = 0.0
    zmax = 0.0

    !! avoid nearby the absorbing boundary
    do j=max(na+margin+1,jbeg_k), min(ny-na-margin, jend_k)
       do i=max(na+margin+1,ibeg_k), min(ny-na-margin, iend_k)
          xmax = max( xmax, real( abs( vx(kob(i,j)+1,i,j) ) ) )
          ymax = max( ymax, real( abs( vy(kob(i,j)+1,i,j) ) ) )
          zmax = max( zmax, real( abs( vz(kob(i,j)+1,i,j) ) ) )
       end do
    end do

  end subroutine kernel__vmax
  !! --------------------------------------------------------------------------------------------------------------------------- !!


  !! --------------------------------------------------------------------------------------------------------------------------- !!
  !>
  !! checkpoint data export
  !<
  !! --
  subroutine kernel__checkpoint( io )

    integer, intent(in) :: io
    integer :: i, j
    !! --

    write(io) r40x, r40y, r40z, r41x, r41y, r41z, r20x, r20y, r20z
    do j=jbeg_m,jend_m; write(io)  Vx(kbeg_m:kend_m,ibeg_m:iend_m,j); end do; deallocate(  Vx )
    do j=jbeg_m,jend_m; write(io)  Vy(kbeg_m:kend_m,ibeg_m:iend_m,j); end do; deallocate(  Vy )
    do j=jbeg_m,jend_m; write(io)  Vz(kbeg_m:kend_m,ibeg_m:iend_m,j); end do; deallocate(  Vz )
    do j=jbeg_m,jend_m; write(io) Sxx(kbeg_m:kend_m,ibeg_m:iend_m,j); end do; deallocate( Sxx )
    do j=jbeg_m,jend_m; write(io) Syy(kbeg_m:kend_m,ibeg_m:iend_m,j); end do; deallocate( Syy )
    do j=jbeg_m,jend_m; write(io) Szz(kbeg_m:kend_m,ibeg_m:iend_m,j); end do; deallocate( Szz )
    do j=jbeg_m,jend_m; write(io) Syz(kbeg_m:kend_m,ibeg_m:iend_m,j); end do; deallocate( Syz )
    do j=jbeg_m,jend_m; write(io) Sxz(kbeg_m:kend_m,ibeg_m:iend_m,j); end do; deallocate( Sxz )
    do j=jbeg_m,jend_m; write(io) Sxy(kbeg_m:kend_m,ibeg_m:iend_m,j); end do; deallocate( Sxy )

    if( nm > 0 ) then
       write(io) c1(1:nm)
       write(io) c2(1:nm)

       do j=jbeg_k,jend_k; do i=ibeg_k,iend_k; write(io) Rxx(kbeg_k:kend_k,i,j,1:nm);  end do; end do; deallocate( Rxx )
       do j=jbeg_k,jend_k; do i=ibeg_k,iend_k; write(io) Ryy(kbeg_k:kend_k,i,j,1:nm);  end do; end do; deallocate( Ryy )
       do j=jbeg_k,jend_k; do i=ibeg_k,iend_k; write(io) Rzz(kbeg_k:kend_k,i,j,1:nm);  end do; end do; deallocate( Rzz )
       do j=jbeg_k,jend_k; do i=ibeg_k,iend_k; write(io) Ryz(kbeg_k:kend_k,i,j,1:nm);  end do; end do; deallocate( Ryz )
       do j=jbeg_k,jend_k; do i=ibeg_k,iend_k; write(io) Rxz(kbeg_k:kend_k,i,j,1:nm);  end do; end do; deallocate( Rxz )
       do j=jbeg_k,jend_k; do i=ibeg_k,iend_k; write(io) Rxy(kbeg_k:kend_k,i,j,1:nm);  end do; end do; deallocate( Rxy )

    end if


  end subroutine kernel__checkpoint
  !! --------------------------------------------------------------------------------------------------------------------------- !!


  !! --------------------------------------------------------------------------------------------------------------------------- !!
  !>
  !! checkpoint data inport
  !<
  !! --
  subroutine kernel__restart( io )

    integer, intent(in) :: io
    integer :: i, j
    !! --


    call memory_allocate()

    read(io) r40x, r40y, r40z, r41x, r41y, r41z, r20x, r20y, r20z

    do j=jbeg_m,jend_m; read(io)  Vx(kbeg_m:kend_m,ibeg_m:iend_m,j); end do;
    do j=jbeg_m,jend_m; read(io)  Vy(kbeg_m:kend_m,ibeg_m:iend_m,j); end do;
    do j=jbeg_m,jend_m; read(io)  Vz(kbeg_m:kend_m,ibeg_m:iend_m,j); end do;
    do j=jbeg_m,jend_m; read(io) Sxx(kbeg_m:kend_m,ibeg_m:iend_m,j); end do;
    do j=jbeg_m,jend_m; read(io) Syy(kbeg_m:kend_m,ibeg_m:iend_m,j); end do;
    do j=jbeg_m,jend_m; read(io) Szz(kbeg_m:kend_m,ibeg_m:iend_m,j); end do;
    do j=jbeg_m,jend_m; read(io) Syz(kbeg_m:kend_m,ibeg_m:iend_m,j); end do;
    do j=jbeg_m,jend_m; read(io) Sxz(kbeg_m:kend_m,ibeg_m:iend_m,j); end do;
    do j=jbeg_m,jend_m; read(io) Sxy(kbeg_m:kend_m,ibeg_m:iend_m,j); end do;

    if( nm > 0 ) then
       read(io) c1(1:nm)
       read(io) c2(1:nm)

       do j=jbeg_k,jend_k; do i=ibeg_k,iend_k; read(io) Rxx(kbeg_k:kend_k,i,j,1:nm); end do; end do;
       do j=jbeg_k,jend_k; do i=ibeg_k,iend_k; read(io) Ryy(kbeg_k:kend_k,i,j,1:nm); end do; end do;
       do j=jbeg_k,jend_k; do i=ibeg_k,iend_k; read(io) Rzz(kbeg_k:kend_k,i,j,1:nm); end do; end do;
       do j=jbeg_k,jend_k; do i=ibeg_k,iend_k; read(io) Ryz(kbeg_k:kend_k,i,j,1:nm); end do; end do;
       do j=jbeg_k,jend_k; do i=ibeg_k,iend_k; read(io) Rxz(kbeg_k:kend_k,i,j,1:nm); end do; end do;
       do j=jbeg_k,jend_k; do i=ibeg_k,iend_k; read(io) Rxy(kbeg_k:kend_k,i,j,1:nm); end do; end do;

    end if

  end subroutine kernel__restart
  !! --------------------------------------------------------------------------------------------------------------------------- !!



  !! --------------------------------------------------------------------------------------------------------------------------- !!
  subroutine memory_allocate
    !!
    !! memory allocation
    !!
    allocate(  Vx(       kbeg_m:kend_m, ibeg_m:iend_m, jbeg_m:jend_m ) )
    allocate(  Vy(       kbeg_m:kend_m, ibeg_m:iend_m, jbeg_m:jend_m ) )
    allocate(  Vz(       kbeg_m:kend_m, ibeg_m:iend_m, jbeg_m:jend_m ) )

    allocate( Sxx(       kbeg_m:kend_m, ibeg_m:iend_m, jbeg_m:jend_m ) )
    allocate( Syy(       kbeg_m:kend_m, ibeg_m:iend_m, jbeg_m:jend_m ) )
    allocate( Szz(       kbeg_m:kend_m, ibeg_m:iend_m, jbeg_m:jend_m ) )
    allocate( Syz(       kbeg_m:kend_m, ibeg_m:iend_m, jbeg_m:jend_m ) )
    allocate( Sxz(       kbeg_m:kend_m, ibeg_m:iend_m, jbeg_m:jend_m ) )
    allocate( Sxy(       kbeg_m:kend_m, ibeg_m:iend_m, jbeg_m:jend_m ) )

    if( nm > 0 ) then
       allocate( c1(nm), c2(nm), d1(nm) )
       allocate( Rxx( kbeg_k:kend_k, ibeg_k:iend_k, jbeg_k:jend_k, 1:nm ) )
       allocate( Ryy( kbeg_k:kend_k, ibeg_k:iend_k, jbeg_k:jend_k, 1:nm ) )
       allocate( Rzz( kbeg_k:kend_k, ibeg_k:iend_k, jbeg_k:jend_k, 1:nm ) )
       allocate( Ryz( kbeg_k:kend_k, ibeg_k:iend_k, jbeg_k:jend_k, 1:nm ) )
       allocate( Rxz( kbeg_k:kend_k, ibeg_k:iend_k, jbeg_k:jend_k, 1:nm ) )
       allocate( Rxy( kbeg_k:kend_k, ibeg_k:iend_k, jbeg_k:jend_k, 1:nm ) )
    end if

  end subroutine memory_allocate
  !! --------------------------------------------------------------------------------------------------------------------------- !!


end module m_kernel
!! ----------------------------------------------------------------------------------------------------------------------------- !!
