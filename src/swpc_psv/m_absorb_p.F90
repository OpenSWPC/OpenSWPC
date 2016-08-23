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
!!   Copyright 2013-2016 Takuto Maeda. All rights reserved. This project is released under the MIT license.
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
  real(SP), allocatable :: axVx(:,:), ayVx(:,:), azVx(:,:)
  real(SP), allocatable :: axVz(:,:), ayVz(:,:), azVz(:,:)
  real(SP), allocatable :: axSxx(:,:), azSxz(:,:)
  real(SP), allocatable :: axSxz(:,:), azSzz(:,:)

  real(SP) :: r20x, r20y, r20z
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
    real(SP) :: hx,hz
    integer  :: idum
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

    !! memory allocation
    allocate(  axVx( kbeg_min:kend, ibeg:iend ) )
    allocate(  azVx( kbeg_min:kend, ibeg:iend ) )
    allocate(  axVz( kbeg_min:kend, ibeg:iend ) )
    allocate(  azVz( kbeg_min:kend, ibeg:iend ) )
    allocate( axSxx( kbeg_min:kend, ibeg:iend ) )
    allocate( azSxz( kbeg_min:kend, ibeg:iend ) )
    allocate( axSxz( kbeg_min:kend, ibeg:iend ) )
    allocate( azSzz( kbeg_min:kend, ibeg:iend ) )

    axVx ( kbeg_min:kend, ibeg:iend ) = 0.0
    azVx ( kbeg_min:kend, ibeg:iend ) = 0.0
    axVz ( kbeg_min:kend, ibeg:iend ) = 0.0
    azVz ( kbeg_min:kend, ibeg:iend ) = 0.0
    axSxx( kbeg_min:kend, ibeg:iend ) = 0.0
    azSxz( kbeg_min:kend, ibeg:iend ) = 0.0
    axSxz( kbeg_min:kend, ibeg:iend ) = 0.0
    azSzz( kbeg_min:kend, ibeg:iend ) = 0.0


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
    real(SP) :: gxc0(4), gxe0(4), gzc0(4), gze0(4)
    real(SP) :: bx, bz
    real(MP) :: dxSxx(kbeg_min:kend), dzSxz(kbeg_min:kend), dxSxz(kbeg_min:kend), dzSzz(kbeg_min:kend)
    !! ----

    !!
    !! Horizontal zero-derivative boundary (for plane wave mode)
    !!
    if( pw_mode ) then
      if( idx == 0 ) then
        !$omp parallel private(k)
        !$omp do schedule(dynamic)
        do k=kbeg_a(1), kend
          Sxx(k,0) = 2 * Sxx(k,1) - Sxx(k,2)
          Szz(k,0) = 2 * Szz(k,1) - Szz(k,2)
          Sxz(k,0) = 2 * Sxz(k,1) - Sxz(k,2)
        end do
        !$omp end do nowait
        !$omp end parallel
      end if
      
      if( idx == nproc_x -1 ) then
        !$omp parallel private(k)
        !$omp do schedule(dynamic)
        do k=kbeg_a(nx), kend
          Sxx(k,nx+1) = 2 * Sxx(k,nx) - Sxx(k,nx-1)
          Szz(k,nx+1) = 2 * Szz(k,nx) - Szz(k,nx-1)
          Sxz(k,nx+1) = 2 * Sxz(k,nx) - Sxz(k,nx-1)
        end do
        !$omp end do nowait
        !$omp end parallel        
      end if
    end if
    
    
    !!
    !! time-marching
    !!
    !$omp parallel &
    !$omp private( dxSxx, dzSzz, dxSxz, dzSxz ) &
    !$omp private( gxc0, gxe0, gzc0, gze0 ) &
    !$omp private( i, k )
    !$omp do &
    !$omp schedule(dynamic)
    do i=ibeg, iend
      
      gxc0(1:4) = gxc(1:4,i)
      gxe0(1:4) = gxe(1:4,i)

      !!
      !! Derivatives
      !!
      do k=kbeg_a(i), kend

        dxSxx(k) = (  Sxx(k  ,i+1) - Sxx(k  ,i  )  ) * r20x
        dzSzz(k) = (  Szz(k+1,i  ) - Szz(k  ,i  )  ) * r20z
        dxSxz(k) = (  Sxz(k  ,i  ) - Sxz(k  ,i-1)  ) * r20x
        dzSxz(k) = (  Sxz(k  ,i  ) - Sxz(k-1,i  )  ) * r20z

      end do


      !!
      !! update velocity
      !!
      do k=kbeg_a(i), kend

        gzc0(1:4) = gzc(1:4,k)
        gze0(1:4) = gze(1:4,k)

        bx = 2.0 / ( rho(k,i) + rho(k,i+1) )
        bz = 2.0 / ( rho(k,i) + rho(k+1,i) )

        !!
        !! Velocity Updates
        !!
        Vx(k,i) = Vx(k,i) &
                + bx * ( gxe0(1) * dxSxx(k)   + gzc0(1) * dzSxz(k)       &
                + gxe0(2) * axSxx(k,i) + gzc0(2) * azSxz(k,i)  ) * dt

        Vz(k,i) = Vz(k,i) &
                + bz * ( gxc0(1) * dxSxz(k)     + gze0(1) * dzSzz(k)       &
                + gxc0(2) * axSxz(k,i) + gze0(2) * azSzz(k,i)  ) * dt


        !!
        !! ADE updates
        !!

        axSxx(k,i) = gxe0(3) * axSxx(k,i) + gxe0(4) * dxSxx(k) * dt
        azSxz(k,i) = gzc0(3) * azSxz(k,i) + gzc0(4) * dzSxz(k) * dt
        axSxz(k,i) = gxc0(3) * axSxz(k,i) + gxc0(4) * dxSxz(k) * dt
        azSzz(k,i) = gze0(3) * azSzz(k,i) + gze0(4) * dzSzz(k) * dt

      end do
    end do
    !$omp end do nowait
    !$omp end parallel

    !$omp barrier

  end subroutine absorb_p__update_vel
  !! --------------------------------------------------------------------------------------------------------------------------- !!

  !! --------------------------------------------------------------------------------------------------------------------------- !!
  subroutine absorb_p__update_stress

    integer :: i, k
    real(SP) :: lam2mu_R, lam_R
    real(SP) :: dxVx_ade, dzVz_ade
    real(SP) :: gxc0(4), gxe0(4), gzc0(4), gze0(4)
    real(SP) :: nnn, pnn, npn, ppn
    real(SP) :: muxz
    real(SP) :: epsl = epsilon(1.0)
    real(MP) :: dxVx(kbeg_min:kend),  dzVx(kbeg_min:kend),  dxVz(kbeg_min:kend),  dzVz(kbeg_min:kend) 

    !! ----

    !!
    !! Horizontal zero-derivative boundary (for plane wave mode)
    !!
    if( pw_mode ) then
      if( idx == 0 ) then
        !$omp parallel private(k)
        !$omp do schedule(dynamic)
        do k=kbeg_a(1), kend
          Vx(k,0) = 2* Vx(k,1)-Vx(k,2)
          Vz(k,0) = 2* Vz(k,1)-Vz(k,2)
        end do
        !$omp end do nowait
        !$omp end parallel
      end if
      
      if( idx == nproc_x -1 ) then
        !$omp parallel private(k)
        !$omp do schedule(dynamic)
        do k=kbeg_a(nx), kend
          Vx(k,nx+1) = 2 * Vx(k,nx) - Vx(k,nx-1)
          Vz(k,nx+1) = 2 * Vz(k,nx) - Vz(k,nx-1)
        end do
        !$omp end do nowait
        !$omp end parallel
      end if
    end if
    
    !!
    !! Time-marching
    !!
    !$omp parallel &
    !$omp private( gxc0, gxe0, gzc0, gze0 ) &
    !$omp private( dxVx, dxVz, dzVx, dzVz ) &
    !$omp private( lam2mu_R, lam_R ) &
    !$omp private( dxVx_ade, dzVz_ade ) &
    !$omp private( i,k ) &
    !$omp private( nnn, pnn, npn, ppn, muxz ) 
    !$omp do &
    !$omp schedule(dynamic)
    do i=ibeg, iend

       gxc0(1:4) = gxc(1:4,i)
       gxe0(1:4) = gxe(1:4,i)

       !!
       !! Derivatives
       !!
       do k=kbeg_a(i), kend

          dxVx(k) = (  Vx(k  ,i  ) - Vx(k  ,i-1)  ) * r20x
          dxVz(k) = (  Vz(k  ,i+1) - Vz(k  ,i  )  ) * r20x
          dzVx(k) = (  Vx(k+1,i  ) - Vx(k  ,i  )  ) * r20z
          dzVz(k) = (  Vz(k  ,i  ) - Vz(k-1,i  )  ) * r20z

       end do

       !!
       !! Update Normal Stress
       !!
       do k=kbeg_a(i), kend

          gzc0(1:4) = gzc(1:4,k)

          !!
          !! update stress components
          !!
          lam2mu_R = ( lam(k,i) + 2 * mu(k,i) )
          lam_R    = lam2mu_R - 2*mu(k,i)


          !!
          !! Normal Stress
          !!
          dxVx_ade = gxc0(1) * dxVx(k) + gxc0(2) * axVx(k,i)
          dzVz_ade = gzc0(1) * dzVz(k) + gzc0(2) * azVz(k,i)
          Sxx(k,i) = Sxx(k,i) + ( lam2mu_R * dxVx_ade + lam_R * dzVz_ade ) * dt
          Szz(k,i) = Szz(k,i) + ( lam2mu_R * dzVz_ade + lam_R * dxVx_ade ) * dt


          !!
          !! ADE
          !!
          axVx(k,i) = gxc0(3) * axVx(k,i) + gxc0(4) * dxVx(k) * dt
          azVz(k,i) = gzc0(3) * azVz(k,i) + gzc0(4) * dzVz(k) * dt

       end do


       !!
       !! Update Shear Stress
       !!
       do k=kbeg_a(i), kend

          gze0(1:4) = gze(1:4,k)

          nnn = mu (k  ,i  )
          pnn = mu (k+1,i  )
          npn = mu (k,  i+1)
          ppn = mu (k+1,i+1)
          muxz = 4*nnn*pnn*npn*ppn / ( nnn*pnn*npn + nnn*pnn*ppn + nnn*npn*ppn + pnn*npn*ppn + epsl)

          Sxz(k,i) = Sxz(k,i) + muxz * ( gxe0(1) * dxVz(k)     + gze0(1) * dzVx(k) &
                                       + gxe0(2) * axVz(k,i) + gze0(2) * azVx(k,i) ) * dt

          azVx(k,i) = gze0(3) * azVx(k,i) + gze0(4) * dzVx(k) * dt
          axVz(k,i) = gxe0(3) * axVz(k,i) + gxe0(4) * dxVz(k) * dt

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
    !! --

    write(io) r20x, r20z
    write(io) kbeg_min

    write(io) gxc(:,ibeg:iend), gxe(:,ibeg:iend)
    write(io) gzc(:,kbeg:kend), gze(:,kbeg:kend)
    deallocate( gxc, gxe, gzc, gze )

    write(io)  axVx(kbeg_min:kend,ibeg:iend);  deallocate( axVx )
    write(io)  azVx(kbeg_min:kend,ibeg:iend);  deallocate( azVx )
    write(io)  axVz(kbeg_min:kend,ibeg:iend);  deallocate( axVz )
    write(io)  azVz(kbeg_min:kend,ibeg:iend);  deallocate( azVz )

    write(io) axSxx(kbeg_min:kend,ibeg:iend);  deallocate( axSxx )
    write(io) azSxz(kbeg_min:kend,ibeg:iend);  deallocate( azSxz )
    write(io) axSxz(kbeg_min:kend,ibeg:iend);  deallocate( axSxz )
    write(io) azSzz(kbeg_min:kend,ibeg:iend);  deallocate( azSzz )

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

    allocate(  axVx(kbeg_min:kend,ibeg:iend) )
    allocate(  azVx(kbeg_min:kend,ibeg:iend) )
    allocate(  axVz(kbeg_min:kend,ibeg:iend) )
    allocate(  azVz(kbeg_min:kend,ibeg:iend) )
    allocate( axSxx(kbeg_min:kend,ibeg:iend) )
    allocate( azSxz(kbeg_min:kend,ibeg:iend) )
    allocate( axSxz(kbeg_min:kend,ibeg:iend) )
    allocate( azSzz(kbeg_min:kend,ibeg:iend) )

    read(io) gxc( :,ibeg:iend ), gxe( :,ibeg:iend )
    read(io) gzc( :,kbeg:kend ), gze( :,kbeg:kend )


    read(io)  axVx(kbeg_min:kend,ibeg:iend)
    read(io)  azVx(kbeg_min:kend,ibeg:iend)
    read(io)  axVz(kbeg_min:kend,ibeg:iend)
    read(io)  azVz(kbeg_min:kend,ibeg:iend)

    read(io) axSxx(kbeg_min:kend,ibeg:iend)
    read(io) azSxz(kbeg_min:kend,ibeg:iend)
    read(io) axSxz(kbeg_min:kend,ibeg:iend)
    read(io) azSzz(kbeg_min:kend,ibeg:iend)

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
    real(SP), save :: d0, a0, b0
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
