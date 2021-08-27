!! ----------------------------------------------------------------------------------------------------------------------------- !!
!>
!! Boundary absorber module: Cerjan's Sponge
!!
!! @copyright
!!   Copyright 2013-2021 Takuto Maeda. All rights reserved. This project is released under the MIT license.
!<
!! ----
#include "m_debug.h"
module m_absorb_c

  !! -- Dependency
  use m_std
  use m_debug
  use m_global
  use m_pwatch
  use m_fdtool

  !! -- Declarations
  implicit none
  private
  save

  !! -- Public Procedures
  public :: absorb_c__setup
  public :: absorb_c__update_stress
  public :: absorb_c__update_vel
  public :: absorb_c__checkpoint
  public :: absorb_c__restart

  real(SP), allocatable :: gx_c(:), gy_c(:), gz_c(:)                !<  attenuator for Q and B.C. for voxel center
  real(SP), allocatable :: gx_b(:), gy_b(:), gz_b(:)                !<  attenuator for Q and B.C. for voxel boundary

contains

  !! --------------------------------------------------------------------------------------------------------------------------- !!
  !>
  !! Setup
  !!
  !! @detail
  !! set Cerjan's sponge function for x(i), y(j) and z(k) directions
  !! as
  !!
  !! @see
  !! note: 2013-00446
  !<
  !! ----
  subroutine absorb_c__setup( io_prm )
    integer, intent(in) :: io_prm
    real(SP), parameter :: alpha = 0.09
    real(SP) :: Lx, Lz
    integer :: i, k
    integer :: io0
    !! ----

    !!
    !! memory allocation and initialize
    !!

    allocate( gx_c( ibeg_m:iend_m ), gx_b( ibeg_m:iend_m ) )
    allocate( gz_c( kbeg_m:kend_m ), gz_b( kbeg_m:kend_m ) )
    gx_c(                                 ibeg_m:iend_m ) = 1.0
    gz_c(       kbeg_m:kend_m                           ) = 1.0
    gx_b(                                 ibeg_m:iend_m ) = 1.0
    gz_b(       kbeg_m:kend_m                           ) = 1.0


    Lx = na * dx
    Lz = na * dz


    !!
    !! Calculate attenuator based on distance
    !!
    !$omp parallel do private(i)
    do i=ibeg, iend
      if( i <= na ) then
        gx_c(i) = exp( - alpha * ( 1.0 -  (   i2x(i, 0.0, real(dx)) )          / Lx )**2 )
        gx_b(i) = exp( - alpha * ( 1.0 -  ( ( i2x(i, 0.0, real(dx)) + real(dx)/2 ) ) / Lx )**2 )
      else if( i >= nx-na+1 ) then
        gx_c(i) = exp( - alpha * ( 1.0 -  (   i2x(i, Nx*real(dx), -real(dx)) + real(dx)/2     ) / Lx )**2 )
        gx_b(i) = exp( - alpha * ( 1.0 -  ( ( i2x(i, Nx*real(dx), -real(dx))         ) ) / Lx )**2 )
      else
        gx_c(i) = 1.0
        gx_b(i) = 1.0
      end if
    end do
    !$omp end parallel do
    !$omp parallel do private(k)
    do k=kbeg, kend
      if( k <= na ) then
        ! if( fullspace_mode ) then
        !   gz_c(k) = exp( - alpha * ( 1.0 -  (   k2z(k, 0.0, real(dz)) )          / Lz )**2 )
        !   gz_b(k) = exp( - alpha * ( 1.0 -  ( ( k2z(k, 0.0, real(dz)) + real(dz)/2 ) ) / Lz )**2 )
        ! else
          gz_c(k) = 1.0
          gz_b(k) = 1.0
        ! end if          
      else if( k >= nz-na+1 ) then
        gz_c(k) = exp( - alpha * ( 1.0 -  (   k2z(k, Nz*real(dz), -real(dz))  + real(dz)/2  ) / Lz )**2 )
        gz_b(k) = exp( - alpha * ( 1.0 -  ( ( k2z(k, Nz*real(dz), -real(dz))  ) ) / Lz )**2 )
      else
        gz_c(k) = 1.0
        gz_b(k) = 1.0
      end if
    end do
    !$omp end parallel do

    !! dummy
    io0 = io_prm

  end subroutine absorb_c__setup
  !! --------------------------------------------------------------------------------------------------------------------------- !!

  !! --------------------------------------------------------------------------------------------------------------------------- !!
  subroutine absorb_c__update_stress

    integer :: i, k

    !$omp parallel do private(i,k)
    do i=ibeg, iend
      do k=kbeg, kend
        Syz(k,i) = Syz(k,i) * gx_c(i) * gz_b(k)
        Sxy(k,i) = Sxy(k,i) * gx_b(i) * gz_c(k)
      end do
    end do
    !$omp end parallel do

  end subroutine absorb_c__update_stress
  !! --------------------------------------------------------------------------------------------------------------------------- !!

  !! --------------------------------------------------------------------------------------------------------------------------- !!
  subroutine absorb_c__update_vel

    integer :: i, k

    !$omp parallel do private(i,k)
    do i=ibeg, iend
      do k=kbeg, kend
        Vy(k,i) = Vy(k,i) * gx_c(i) * gz_c(k)
      end do
    end do
    !$omp end parallel do

  end subroutine absorb_c__update_vel
  !! --------------------------------------------------------------------------------------------------------------------------- !!

  subroutine absorb_c__checkpoint( io )
    integer, intent(in) :: io

    write( io ) gx_c( ibeg_m:iend_m ), gx_b( ibeg_m:iend_m )
    write( io ) gz_c( kbeg_m:kend_m ), gz_b( kbeg_m:kend_m )

  end subroutine absorb_c__checkpoint

  subroutine absorb_c__restart( io )
    integer, intent(in) :: io

    allocate( gx_c( ibeg_m:iend_m ), gx_b( ibeg_m:iend_m ) )
    allocate( gz_c( kbeg_m:kend_m ), gz_b( kbeg_m:kend_m ) )
    read( io ) gx_c( ibeg_m:iend_m ), gx_b( ibeg_m:iend_m )
    read( io ) gz_c( kbeg_m:kend_m ), gz_b( kbeg_m:kend_m )

  end subroutine absorb_c__restart

end module m_absorb_c
!! ----------------------------------------------------------------------------------------------------------------------------- !!
