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
!!   Copyright 2013-2020 Takuto Maeda. All rights reserved. This project is released under the MIT license.
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
  real(SP), allocatable :: gzc(:,:), gze(:,:) !< damping profile along x at center/edge of voxel

  !! ADE variables
  real(SP), allocatable :: axVy(:,:), azVy(:,:)
  real(SP), allocatable :: axSxy(:,:)
  real(SP), allocatable :: azSyz(:,:)

  real(SP) :: r20x, r20z

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
    integer  :: i, k
    real(SP) :: hx, hz
    integer  :: idum
    real :: coef
    !! ----

    !!
    !! derivative coefficient
    !!
    r20x = 1.0 / dx
    r20z = 1.0 / dz

    !!
    !! damping profile
    !!
    allocate( gxc(4,ibeg:iend), gxe(4,ibeg:iend) )
    allocate( gzc(4,kbeg:kend), gze(4,kbeg:kend) )

    hx = na * dx
    hz = na * dz

    do i=ibeg, iend
      call damping_profile( xc(i),              hx, xbeg, xend, gxc(:,i) )
      call damping_profile( xc(i)+real(dx)/2.0, hx, xbeg, xend, gxe(:,i) )
    end do
    do k=kbeg, kend
      call damping_profile( zc(k),              hz, zbeg, zend, gzc(:,k) )
      call damping_profile( zc(k)+real(dz)/2.0, hz, zbeg, zend, gze(:,k) )
    end do

    !!
    !! PML region definition
    !!
    kbeg_min = minval( kbeg_a(:) )
    if( fullspace_mode ) kbeg_min = kbeg

    !! memory allocation
    allocate(  axVy( kbeg_min:kend, ibeg:iend ) )
    allocate(  azVy( kbeg_min:kend, ibeg:iend ) )
    allocate( axSxy( kbeg_min:kend, ibeg:iend ) )
    allocate( azSyz( kbeg_min:kend, ibeg:iend ) )

    axVy ( kbeg_min:kend, ibeg:iend ) = 0.0
    azVy ( kbeg_min:kend, ibeg:iend ) = 0.0
    axSxy( kbeg_min:kend, ibeg:iend ) = 0.0
    azSyz( kbeg_min:kend, ibeg:iend ) = 0.0

    idum = io_prm


  end subroutine absorb_p__setup
  !! --------------------------------------------------------------------------------------------------------------------------- !!

  !! --------------------------------------------------------------------------------------------------------------------------- !!
  !>
  !! Update velocity component in PML layer
  !<
  !! ----
  subroutine absorb_p__update_vel

    integer :: i, k
    real(SP) :: gxc0(4), gzc0(4)
    real(SP) :: by
    real(MP) :: dzSyz(kbeg:kend), dxSxy(kbeg:kend)
    !! ----

    !!
    !! Horizontal zero-derivative boundary (for plane wave mode)
    !!
    if( pw_mode ) then
      if( idx == 0 ) then
        !$omp parallel private(k)
        !$omp do schedule(dynamic)
        do k=kbeg_a(1), kend
          Syz(k,0) = 2 * Syz(k,1) - Syz(k,2)
          Sxy(k,0) = 2 * Sxy(k,1) - Sxy(k,2)
        end do
        !$omp end do nowait
        !$omp end parallel
      end if

      if( idx == nproc_x -1 ) then
        !$omp parallel private(k)
        !$omp do schedule(dynamic)
        do k=kbeg_a(nx), kend
          Syz(k,nx+1) = 2 * Syz(k,nx) - Syz(k,nx-1)
          Sxy(k,nx+1) = 2 * Sxy(k,nx) - Sxy(k,nx-1)
        end do
        !$omp end do nowait
        !$omp end parallel
      end if
    end if

    !!
    !! time-marching
    !!
    !$omp parallel &
    !$omp private( gxc0, gzc0, dzSyz, dxSxy, by, i, k )
    !$omp do schedule(dynamic)
    do i=ibeg, iend

      gxc0(1:4) = gxc(1:4,i)

      !!
      !! Derivatives
      !!
      do k=kbeg_a(i), kend

        dzSyz(k) = (  Syz(k  ,i ) - Syz(k-1,i   )  ) * r20z
        dxSxy(k) = (  Sxy(k  ,i ) - Sxy(k  ,i-1 )  ) * r20x

      end do



      !!
      !! update velocity
      !!
      do k=kbeg_a(i), kend

        gzc0(1:4) = gzc(1:4,k)

        !!
        !! Velocity Updates
        !!
        by = 1.0 /  rho(k,i)

        Vy(k,i) = Vy(k,i) &
            + by * ( gxc0(1) * dxSxy(k)   + gzc0(1) * dzSyz(k)       &
            + gxc0(2) * axSxy(k,i) + gzc0(2) * azSyz(k,i)  ) * dt

        !!
        !! ADE updates
        !!
        axSxy(k,i) = gxc0(3) * axSxy(k,i) + gxc0(4) * dxSxy(k) * dt
        azSyz(k,i) = gzc0(3) * azSyz(k,i) + gzc0(4) * dzSyz(k) * dt

      end do
    end do
    !$omp end do nowait
    !$omp end parallel

    if( fullspace_mode ) then
      !$omp parallel &
      !$omp private( gxc0, gzc0, dzSyz, dxSxy, by, i, k )
      !$omp do schedule(dynamic)
      do i=ibeg_k, iend_k

        gxc0(1:4) = gxc(1:4,i)

        !!
        !! Derivatives
        !!
        do k=1,na

          dzSyz(k) = (  Syz(k  ,i ) - Syz(k-1,i   )  ) * r20z
          dxSxy(k) = (  Sxy(k  ,i ) - Sxy(k  ,i-1 )  ) * r20x

        end do

        !!
        !! update velocity
        !!
        do k=1, na

          gzc0(1:4) = gzc(1:4,k)

          !!
          !! Velocity Updates
          !!
          by = 1.0 /  rho(k,i)

          Vy(k,i) = Vy(k,i) &
              + by * ( gxc0(1) * dxSxy(k)   + gzc0(1) * dzSyz(k)       &
              + gxc0(2) * axSxy(k,i) + gzc0(2) * azSyz(k,i)  ) * dt

          !!
          !! ADE updates
          !!
          axSxy(k,i) = gxc0(3) * axSxy(k,i) + gxc0(4) * dxSxy(k) * dt
          azSyz(k,i) = gzc0(3) * azSyz(k,i) + gzc0(4) * dzSyz(k) * dt
          
        end do
      end do
      !$omp end do nowait
      !$omp end parallel
    end if
  

    !$omp barrier

  end subroutine absorb_p__update_vel
  !! --------------------------------------------------------------------------------------------------------------------------- !!

  !! --------------------------------------------------------------------------------------------------------------------------- !!
  subroutine absorb_p__update_stress

    integer :: i, k
    real(SP) :: gxe0(4), gze0(4)
    real(SP) :: nnn, pnn, npn
    real(SP) :: muxy, muyz
    real(SP) :: epsl = epsilon(1.0)
    real(MP) :: dxVy(kbeg:kend), dzVy(kbeg:kend)
    !! ----

    !!
    !! Horizontal zero-derivative boundary (for plane wave mode)
    !!
    if( pw_mode ) then
      if( idx == 0 ) then
        !$omp parallel private(k)
        !$omp do schedule(dynamic)
        do k=kbeg_a(1), kend
          Vy(k,0) = 2* Vy(k,1)-Vy(k,2)
        end do
        !$omp end do nowait
        !$omp end parallel
      end if

      if( idx == nproc_x -1 ) then
        !$omp parallel private(k)
        !$omp do schedule(dynamic)
        do k=kbeg_a(nx), kend
          Vy(k,nx+1) = 2 * Vy(k,nx) - Vy(k,nx-1)
        end do
        !$omp end do nowait
        !$omp end parallel
      end if
    end if

    !!
    !! Time-marching
    !!
    !$omp parallel &
    !$omp private(i, k, gxe0, gze0, dxVy, dzVy, nnn, pnn, npn, muxy, muyz )
    !$omp do schedule(dynamic)
    do i=ibeg, iend

      gxe0(1:4) = gxe(1:4,i)

      !!
      !! Derivatives
      !!
      do k=kbeg_a(i), kend

        dxVy(k) = (  Vy(k  ,i+1) - Vy(k,i)  ) * r20x
        dzVy(k) = (  Vy(k+1,i  ) - Vy(k,i)  ) * r20z

      end do


      !!
      !! Update Shear Stress
      !!
      do k=kbeg_a(i), kend

        gze0(1:4) = gze(1:4,k)


        !!
        !! effective rigidity for shear stress components
        !!

        nnn = mu (k  ,i )
        pnn = mu (k+1,i )
        npn = mu (k,  i+1)
        muxy = 2*nnn*npn / ( nnn + npn + epsl )
        muyz = 2*nnn*pnn / ( nnn + pnn + epsl )

        axVy(k,i) = gxe0(3) * axVy(k,i) + gxe0(4) * dxVy(k) * dt
        azVy(k,i) = gze0(3) * azVy(k,i) + gze0(4) * dzVy(k) * dt
        Syz(k,i) = Syz(k,i) + muyz* ( gze0(1) * dzVy(k) + gze0(2) * azVy(k,i) ) * dt
        Sxy(k,i) = Sxy(k,i) + muxy* ( gxe0(1) * dxVy(k) + gxe0(2) * axVy(k,i) ) * dt

      end do
    end do
    !$omp end do nowait
    !$omp end parallel

    if( fullspace_mode ) then
      !$omp parallel &
      !$omp private(i, k, gxe0, gze0, dxVy, dzVy, nnn, pnn, npn, muxy, muyz )
      !$omp do schedule(dynamic)
      do i=ibeg_k, iend_k

        gxe0(1:4) = gxe(1:4,i)

        !!
        !! Derivatives
        !!
        do k=kbeg, na

          dxVy(k) = (  Vy(k  ,i+1) - Vy(k,i)  ) * r20x
          dzVy(k) = (  Vy(k+1,i  ) - Vy(k,i)  ) * r20z

        end do


        !!
        !! Update Shear Stress
        !!
        do k=kbeg, na

          gze0(1:4) = gze(1:4,k)


          !!
          !! effective rigidity for shear stress components
          !!

          nnn = mu (k  ,i )
          pnn = mu (k+1,i )
          npn = mu (k,  i+1)
          muxy = 2*nnn*npn / ( nnn + npn + epsl )
          muyz = 2*nnn*pnn / ( nnn + pnn + epsl )

          axVy(k,i) = gxe0(3) * axVy(k,i) + gxe0(4) * dxVy(k) * dt
          azVy(k,i) = gze0(3) * azVy(k,i) + gze0(4) * dzVy(k) * dt
          Syz(k,i) = Syz(k,i) + muyz* ( gze0(1) * dzVy(k) + gze0(2) * azVy(k,i) ) * dt
          Sxy(k,i) = Sxy(k,i) + muxy* ( gxe0(1) * dxVy(k) + gxe0(2) * axVy(k,i) ) * dt

        end do
      end do
      !$omp end do nowait
      !$omp end parallel

    end if
    
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
    !! --

    write(io) r20x, r20z
    write(io) kbeg_min
    write(io) gxc(:,ibeg:iend), gxe(:,ibeg:iend)
    write(io) gzc(:,kbeg:kend), gze(:,kbeg:kend)
    deallocate( gxc, gxe, gzc, gze )

    write(io)  axVy(kbeg_min:kend,ibeg:iend); deallocate( axVy )
    write(io)  azVy(kbeg_min:kend,ibeg:iend); deallocate( azVy )
    write(io) axSxy(kbeg_min:kend,ibeg:iend); deallocate( axSxy )
    write(io) azSyz(kbeg_min:kend,ibeg:iend); deallocate( azSyz )

  end subroutine absorb_p__checkpoint
  !! --------------------------------------------------------------------------------------------------------------------------- !!

  !! --------------------------------------------------------------------------------------------------------------------------- !!
  subroutine absorb_p__restart(io)

    integer, intent(in) :: io

    read(io) r20x, r20z
    read(io) kbeg_min

    !! memory allocation
    allocate( gxc( 4,ibeg:iend ), gxe( 4,ibeg:iend ) )
    allocate( gzc( 4,kbeg:kend ), gze( 4,kbeg:kend ) )

    allocate(  axVy(kbeg_min:kend,ibeg:iend) )
    allocate(  azVy(kbeg_min:kend,ibeg:iend) )
    allocate( axSxy(kbeg_min:kend,ibeg:iend) )
    allocate( azSyz(kbeg_min:kend,ibeg:iend) )

    read(io) gxc( :,ibeg:iend ), gxe( :,ibeg:iend )
    read(io) gzc( :,kbeg:kend ), gze( :,kbeg:kend )


    read(io)  axVy(kbeg_min:kend,ibeg:iend)
    read(io)  azVy(kbeg_min:kend,ibeg:iend)

    read(io) axSxy(kbeg_min:kend,ibeg:iend)
    read(io) azSyz(kbeg_min:kend,ibeg:iend)

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
    integer, parameter :: pa = 2
    integer, parameter :: pb = 2
    real(SP), parameter :: cp = 6.0 !! assumed P-wave velocity
    real :: d, a, b, xx

    R0 = 10**( - ( log10( real(na) ) - 1 ) / log10( 2.0 )  - 3.0 )
    d0 = - ( 1.0 / (2.0*H) ) * ( pd +1 ) * cp * log( R0 )
    b0 = 7.0
    a0 = PI * fcut

    if( x <= xbeg + H ) then
      xx = ( xbeg + H ) - x
    else if ( x >= xend - H ) then
      xx = x - ( xend - H )
    else
      xx = 0.0 !! no absorption
    end if

    d = d0 * abs( xx / H )**pd
    a = a0 * ( 1.0 - abs( xx / H )**pa )
    b = 1.0 + ( b0 - 1.0 ) * abs( xx / H )**pb

    g(1) = (  ( 1.0 + ( dt / 2.0 ) * a ) / b     ) / ( 1.0 + ( dt / 2.0 ) * ( a + d / b ) )
    g(2) = ( -1.0 / b                            ) / ( 1.0 + ( dt / 2.0 ) * ( a + d / b ) )
    g(3) = (  1.0 - ( dt / 2.0 ) * ( a + d / b ) ) / ( 1.0 + ( dt / 2.0 ) * ( a + d / b ) )
    g(4) = (  d / b                              ) / ( 1.0 + ( dt / 2.0 ) * ( a + d / b ) )

  end subroutine damping_profile
  !! --------------------------------------------------------------------------------------------------------------------------- !!

end module m_absorb_p
!! ----------------------------------------------------------------------------------------------------------------------------- !!
