!! ----------------------------------------------------------------------------------------------------------------------------- !!
!>
!! Computation kernel for FDM numerical simulation
!!
!! @copyright
!!   Copyright 2013-2021 Takuto Maeda. All rights reserved. This project is released under the MIT license.
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
  real(MP)              :: r40x, r40z
  real(MP)              :: r41x, r41z
  real(MP)              :: r20x, r20z
  real(SP), allocatable :: c1(:), c2(:)

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

    if( .not. medium__initialized() ) then
      write(STDERR,'(A)') 'ERROR [kernel__setup]: call medium__setup() before kernel__setup()'
      stop
    end if

    r40x = C40 / dx
    r40z = C40 / dz
    r41x = C41 / dx
    r41z = C41 / dz
    r20x = C20 / dx
    r20z = C20 / dz

    call memory_allocate()

    !!
    !! initialize
    !!
    Vy  (       kbeg_m:kend_m, ibeg_m:iend_m ) = 0.0_MP
    Syz (       kbeg_m:kend_m, ibeg_m:iend_m ) = 0.0_MP
    Sxy (       kbeg_m:kend_m, ibeg_m:iend_m ) = 0.0_MP
    if( nm > 0 ) then
      Ryz ( 1:nm, kbeg_k:kend_k, ibeg_k:iend_k ) = 0.0
      Rxy ( 1:nm, kbeg_k:kend_k, ibeg_k:iend_k ) = 0.0

      do m=1, nm
        c1(m) = ( 2 * ts(m) - dt ) / ( 2 * ts(m) + dt )
        c2(m) = ( 2              ) / ( 2 * ts(m) + dt ) / nm
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

    integer :: i, k
    real(SP) :: by ! buoyancy
    real(MP) :: dzSyz(kbeg:kend), dxSxy(kbeg:kend)
    !! ----

    call pwatch__on("kernel__update_vel")

    !$omp parallel &
    !$omp private(dzSyz, dxSxy, i, k, by) 
    !$omp do &
#ifdef _ES
    !$omp schedule(static)
#else
    !$omp schedule(static,1)
#endif    
    do i=ibeg_k, iend_k

      !! derivateives
      do k=kbeg_k, kend_k
        dzSyz(k) = (  Syz(k  ,i  ) - Syz(k-1,i  )  ) * r40z  -  (  Syz(k+1,i  ) - Syz(k-2,i  )  ) * r41z
        dxSxy(k) = (  Sxy(k  ,i  ) - Sxy(k  ,i-1)  ) * r40x  -  (  Sxy(k  ,i+1) - Sxy(k  ,i-2)  ) * r41x
      end do

      !! surfaces
#ifdef _ES
      !NEC$ novector
#endif
      do k=kfs_top(i), kfs_bot(i)
        dzSyz(k) = (  Syz(k  ,i  ) - Syz(k-1,i  )  ) * r20z
        dxSxy(k) = (  Sxy(k  ,i  ) - Sxy(k  ,i-1)  ) * r20x
      end do

#ifdef _ES
      !NEC$ novector
#endif
      do k=kob_top(i), kob_bot(i)
        dzSyz(k) = (  Syz(k  ,i  ) - Syz(k-1,i  )  ) * r20z
        dxSxy(k) = (  Sxy(k  ,i  ) - Sxy(k  ,i-1)  ) * r20x
      end do

      !!
      !! update velocity
      !!
      do k=kbeg_k, kend_k

        !!
        !! effective buoyancy
        !!
        by = 1.0 /  rho(k,i)

        !!
        !! update velocity
        !!
        Vy(k,i) =  Vy(k,i) + by * ( dxSxy(k) + dzSyz(k) ) * dt

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

    integer :: i, k, m
    real(SP) :: nnn, pnn, npn
    real(SP) :: mu_xy, mu_yz
    real(SP) :: taus1, taus_plus1
    real(SP) :: Ryz_o, Rxy_o
    real(SP) :: Ryz_n, Rxy_n
    real(SP) :: f_Rxy, f_Ryz
    real(SP) :: epsl = epsilon(1.0)
    real(MP) :: dxVy(kbeg:kend), dzVy(kbeg:kend)
    !! ----

    call pwatch__on("kernel__update_stress")


    !$omp parallel  &
    !$omp private( dxVy, dzVy, nnn, pnn, npn, mu_yz, mu_xy, taus1, taus_plus1, f_Ryz, f_Rxy, Ryz_o, Rxy_o, Ryz_n, Rxy_n )
    !$omp do &
#ifdef _ES
    !$omp schedule(static)
#else
    !$omp schedule(static,1)
#endif
    do i=ibeg_k, iend_k

      !!
      !! Derivatives
      !!
      do k=kbeg_k, kend_k

        dxVy(k) = (  Vy(k  ,i+1) - Vy(k  ,i  )  ) * r40x  -  (  Vy(k  ,i+2) - Vy(k  ,i-1)  ) * r41x
        dzVy(k) = (  Vy(k+1,i  ) - Vy(k  ,i  )  ) * r40z  -  (  Vy(k+2,i  ) - Vy(k-1,i  )  ) * r41z

      end do

      !! free surface
#ifdef _ES
      !NEC$ novector
#endif
      do k=kfs_top(i), kfs_bot(i)

        dxVy(k) = (  Vy(k  ,i+1) - Vy(k  ,i  )  ) * r20x
        dzVy(k) = (  Vy(k+1,i  ) - Vy(k  ,i  )  ) * r20z

      end do

      !! seafloor
#ifdef _ES
      !NEC$ novector
#endif
      do k=kob_top(i), kob_bot(i)

        dxVy(k) = (  Vy(k  ,i+1) - Vy(k  ,i  )  ) * r20x
        dzVy(k) = (  Vy(k+1,i  ) - Vy(k  ,i  )  ) * r20z

      end do


      do k=kbeg_k, kend_k

        !!
        !! effective rigidity for shear stress components
        !!

        nnn = mu (k  ,i )
        pnn = mu (k+1,i )
        npn = mu (k,  i+1)
        mu_xy = 2*nnn*npn / ( nnn + npn + epsl )
        mu_yz = 2*nnn*pnn / ( nnn + pnn + epsl )
        taus1 = taus(k,i)

        !!
        !! update memory variables
        !!

        !! working variables for combinations of velocity derivatives
        f_Ryz =  mu_yz * taus1 * dzVy(k)
        f_Rxy =  mu_xy * taus1 * dxVy(k)


        Ryz_o = 0.0
        Rxy_o = 0.0
        Ryz_n = 0.0
        Rxy_n = 0.0

        if( nm > 0 ) then
          !! previous memory variables
          Ryz_o = sum( Ryz( 1:nm, k,i) )
          Rxy_o = sum( Rxy( 1:nm, k,i) )

          !! Crank-Nicolson Method for avoiding stiff solution
          do m=1, nm
            Ryz(m,k,i) = c1(m) * Ryz(m,k,i) - c2(m) * f_Ryz * dt
            Rxy(m,k,i) = c1(m) * Rxy(m,k,i) - c2(m) * f_Rxy * dt
          end do

          !! new memory variables
          Ryz_n = sum( Ryz( 1:nm, k,i) )
          Rxy_n = sum( Rxy( 1:nm, k,i) )
        end if

        !!
        !! update stress components
        !!
        taus_plus1 = 1 + taus1

        Syz (k,i) = Syz (k,i) + (  mu_yz*taus_plus1*dzVy(k) + ( Ryz_n+Ryz_o )/2 ) * dt
        Sxy (k,i) = Sxy (k,i) + (  mu_xy*taus_plus1*dxVy(k) + ( Rxy_n+Rxy_o )/2 ) * dt

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
  subroutine kernel__vmax( ymax )
    real(SP), intent(out) :: ymax
    integer :: i

    ymax = 0.0
    do i=ibeg_k, iend_k
      ymax = max( ymax, abs( vy(kob(i)+1,i) ) )
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
    integer :: m
    !! --

    write(io) r40x, r40z, r41x, r41z, r20x, r20z

    write(io) Vy (       kbeg_m:kend_m, ibeg_m:iend_m )
    write(io) Syz(       kbeg_m:kend_m, ibeg_m:iend_m )
    write(io) Sxy(       kbeg_m:kend_m, ibeg_m:iend_m )

    do m=1, nm
      write(io) c1(m), c2(m)
      write(io) Ryz( m, kbeg_k:kend_k, ibeg_k:iend_k )
      write(io) Rxy( m, kbeg_k:kend_k, ibeg_k:iend_k )
    end do

  end subroutine kernel__checkpoint
  !! --------------------------------------------------------------------------------------------------------------------------- !!


  !! --------------------------------------------------------------------------------------------------------------------------- !!
  !>
  !! checkpoint data inport
  !<
  !! --
  subroutine kernel__restart( io )

    integer, intent(in) :: io
    integer :: m
    !! --


    call memory_allocate()

    read(io) r40x, r40z, r41x, r41z, r20x, r20z
    read(io) Vy (       kbeg_m:kend_m, ibeg_m:iend_m )
    read(io) Syz(       kbeg_m:kend_m, ibeg_m:iend_m )
    read(io) Sxy(       kbeg_m:kend_m, ibeg_m:iend_m )

    do m=1, nm
      read(io) c1(m), c2(m)
      read(io) Ryz( m, kbeg_k:kend_k, ibeg_k:iend_k )
      read(io) Rxy( m, kbeg_k:kend_k, ibeg_k:iend_k )
    end do

  end subroutine kernel__restart
  !! --------------------------------------------------------------------------------------------------------------------------- !!



  !! --------------------------------------------------------------------------------------------------------------------------- !!
  subroutine memory_allocate
    !!
    !! memory allocation
    !!
    allocate(  Vy(       kbeg_m:kend_m, ibeg_m:iend_m ) )

    allocate( Syz(       kbeg_m:kend_m, ibeg_m:iend_m ) )
    allocate( Sxy(       kbeg_m:kend_m, ibeg_m:iend_m ) )

    if( nm > 0 ) then
      allocate( Ryz( 1:nm, kbeg_k:kend_k, ibeg_k:iend_k ) )
      allocate( Rxy( 1:nm, kbeg_k:kend_k, ibeg_k:iend_k ) )
      allocate( c1(1:nm), c2(1:nm) )
    end if

  end subroutine memory_allocate
  !! --------------------------------------------------------------------------------------------------------------------------- !!


end module m_kernel
!! ----------------------------------------------------------------------------------------------------------------------------- !!
