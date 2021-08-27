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
    Vx  (       kbeg_m:kend_m, ibeg_m:iend_m ) = 0.0_MP
    Vz  (       kbeg_m:kend_m, ibeg_m:iend_m ) = 0.0_MP
    Sxx (       kbeg_m:kend_m, ibeg_m:iend_m ) = 0.0_MP
    Szz (       kbeg_m:kend_m, ibeg_m:iend_m ) = 0.0_MP
    Sxz (       kbeg_m:kend_m, ibeg_m:iend_m ) = 0.0_MP
    if( nm > 0 ) then
      Rxx ( 1:nm, kbeg_m:kend_m, ibeg_m:iend_m ) = 0.0
      Rzz ( 1:nm, kbeg_m:kend_m, ibeg_m:iend_m ) = 0.0
      Rxz ( 1:nm, kbeg_m:kend_m, ibeg_m:iend_m ) = 0.0

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
    real(SP) :: bx, bz ! buoyancy
    real(MP) :: dxSxx(kbeg:kend), dxSxz(kbeg:kend), dzSxz(kbeg:kend), dzSzz(kbeg:kend)
    !! ----

    call pwatch__on("kernel__update_vel")


    !$omp parallel &
    !$omp private(dxSxx, dzSzz, dxSxz, dzSxz ) &
    !$omp private(i,k) &
    !$omp private(bx, bz) 
    !$omp do &
#ifdef _ES
    !$omp schedule(static)
#else
    !$omp schedule(static,1)
#endif
    do i=ibeg_k, iend_k

      !! derivateives
      do k=kbeg_k, kend_k
        dxSxx(k) = (  Sxx(k  ,i+1) - Sxx(k  ,i  )  ) * r40x  -  (  Sxx(k  ,i+2) - Sxx(k  ,i-1)  ) * r41x
        dzSzz(k) = (  Szz(k+1,i  ) - Szz(k  ,i  )  ) * r40z  -  (  Szz(k+2,i  ) - Szz(k-1,i  )  ) * r41z
        dxSxz(k) = (  Sxz(k  ,i  ) - Sxz(k  ,i-1)  ) * r40x  -  (  Sxz(k  ,i+1) - Sxz(k  ,i-2)  ) * r41x
        dzSxz(k) = (  Sxz(k  ,i  ) - Sxz(k-1,i  )  ) * r40z  -  (  Sxz(k+1,i  ) - Sxz(k-2,i  )  ) * r41z
      end do

      !! surfaces
#ifdef _ES
      !NEC$ novector
#endif
      do k=kfs_top(i), kfs_bot(i)
        dxSxx(k) = (  Sxx(k  ,i+1) - Sxx(k  ,i  )  ) * r20x
        dzSzz(k) = (  Szz(k+1,i  ) - Szz(k  ,i  )  ) * r20z
        dxSxz(k) = (  Sxz(k  ,i  ) - Sxz(k  ,i-1)  ) * r20x
        dzSxz(k) = (  Sxz(k  ,i  ) - Sxz(k-1,i  )  ) * r20z
      end do

#ifdef _ES
      !NEC$ novector
#endif
      do k=kob_top(i), kob_bot(i)
        dxSxx(k) = (  Sxx(k  ,i+1) - Sxx(k  ,i  )  ) * r20x
        dzSzz(k) = (  Szz(k+1,i  ) - Szz(k  ,i  )  ) * r20z
        dxSxz(k) = (  Sxz(k  ,i  ) - Sxz(k  ,i-1)  ) * r20x
        dzSxz(k) = (  Sxz(k  ,i  ) - Sxz(k-1,i  )  ) * r20z
      end do


      !! top
      dzSzz(1)    = (  Szz(2  ,i  ) - Szz(1  ,i  )  ) * r20z
      dzSxz(1)    = (  Sxz(1  ,i  ) - 0.0        ) * r20z
      dzSxz(2)    = (  Sxz(2  ,i  ) - Sxz(1,  i  )  ) * r20z

      !! bottom
      dzSzz(nz-1) = (  Szz(nz  ,i  ) - Szz(nz-1,i  )  ) * r20z
      dzSzz(nz)   = (  0.0        - Szz(nz  ,i  )  ) * r20z
      dzSxz(nz)   = (  Sxz(nz  ,i  ) - Sxz(nz-1,i  )  ) * r20z

      !! i-boundary 
#ifdef _ES
      !NEC$ novector
#endif
      if( i == 1 ) then
        do k=kbeg_k, kend_k
          dxSxx(k) = (  Sxx(k  ,2  ) - Sxx(k  ,1  )  ) * r20x
          dxSxz(k) = (  Sxz(k  ,1  ) - 0.0            ) * r20x
        end do
      else if ( i == 2 ) then
        do k=kbeg_k, kend_k
          dxSxz(k) = (  Sxz(k  ,2  ) - Sxz(k  ,1)  ) * r20x
        end do
      else if ( i == nx-1 ) then
        do k=kbeg_k, kend_k
          dxSxx(k) = (  Sxx(k  ,nx) - Sxx(k  ,nx-1)  ) * r20x
        end do
      else if ( i == nx ) then
        do k=kbeg_k, kend_k
          dxSxx(k) = (  0.0        - Sxx(k  ,nx  )  ) * r20x
          dxSxz(k) = (  Sxz(k  ,nx  ) - Sxz(k  ,nx-1)  ) * r20x
        end do
      end if



      !!
      !! update velocity
      !!
      do k=kbeg_k, kend_k

        !!
        !! effective buoyancy
        !!

        bx = 2.0 / ( rho(k,i) + rho(k,i+1) )
        bz = 2.0 / ( rho(k,i) + rho(k+1,i) )

        !!
        !! update velocity
        !!
        Vx(k,i) = Vx(k,i) + bx * ( dxSxx(k) + dzSxz(k) ) * dt
        Vz(k,i) = Vz(k,i) + bz * ( dxSxz(k) + dzSzz(k) ) * dt

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
    real(SP) :: mu2, lam2mu
    real(SP) :: nnn, pnn, npn, ppn
    real(SP) :: mu_xz
    real(SP) :: taup1, taus1, taup_plus1, taus_plus1
    real(SP) :: d2v2, dxVx_dzVz, dxVz_dzVx
    real(SP) :: Rxx_o, Rzz_o, Rxz_o
    real(SP) :: Rxx_n, Rzz_n, Rxz_n
    real(SP) :: f_Rxx, f_Rzz, f_Rxz
    real(SP) :: epsl = epsilon(1.0)
    real(MP) :: dxVx(kbeg:kend), dxVz(kbeg:kend), dzVx(kbeg:kend), dzVz(kbeg:kend)
    !! ----

    call pwatch__on("kernel__update_stress")

    !$omp parallel  &
    !$omp private( dxVx, dxVz, dzVx, dzVz ) &
    !$omp private( mu2, lam2mu ) &
    !$omp private( taup1, taus1, nnn, pnn, npn, ppn, mu_xz, taup_plus1, taus_plus1 ) &
    !$omp private( d2v2, dxVx_dzVz, dxVz_dzVx ) &
    !$omp private( f_Rxx, f_Rzz, f_Rxz ) &
    !$omp private( Rxx_o, Rzz_o, Rxz_o, Rxx_n, Rzz_n, Rxz_n ) &
    !$omp private( i, k, m )
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

        dxVx(k) = (  Vx(k  ,i  ) - Vx(k  ,i-1)  ) * r40x  -  (  Vx(k  ,i+1) - Vx(k  ,i-2)  ) * r41x
        dxVz(k) = (  Vz(k  ,i+1) - Vz(k  ,i  )  ) * r40x  -  (  Vz(k  ,i+2) - Vz(k  ,i-1)  ) * r41x
        dzVx(k) = (  Vx(k+1,i  ) - Vx(k  ,i  )  ) * r40z  -  (  Vx(k+2,i  ) - Vx(k-1,i  )  ) * r41z
        dzVz(k) = (  Vz(k  ,i  ) - Vz(k-1,i  )  ) * r40z  -  (  Vz(k+1,i  ) - Vz(k-2,i  )  ) * r41z

      end do

      !! free surface
#ifdef _ES
      !NEC$ novector
#endif
      do k=kfs_top(i), kfs_bot(i)

        dxVx(k) = (  Vx(k  ,i  ) - Vx(k  ,i-1)  ) * r20x
        dxVz(k) = (  Vz(k  ,i+1) - Vz(k  ,i  )  ) * r20x
        dzVx(k) = (  Vx(k+1,i  ) - Vx(k  ,i  )  ) * r20z
        dzVz(k) = (  Vz(k  ,i  ) - Vz(k-1,i  )  ) * r20z

      end do

      !! seafloor
#ifdef _ES
      !NEC$ novector
#endif
      do k=kob_top(i), kob_bot(i)

        dxVx(k) = (  Vx(k  ,i  ) - Vx(k  ,i-1)  ) * r20x
        dxVz(k) = (  Vz(k  ,i+1) - Vz(k  ,i  )  ) * r20x
        dzVx(k) = (  Vx(k+1,i  ) - Vx(k  ,i  )  ) * r20z
        dzVz(k) = (  Vz(k  ,i  ) - Vz(k-1,i  )  ) * r20z

      end do

      !!
      !! vertical edge
      !!

      !! top
      dzVx(1) = (  Vx(2  ,i  ) - Vx(1  ,i  )  ) * r20z
      dzVz(1) = (  Vz(1  ,i  ) - 0.0       ) * r20z
      dzVz(2) = (  Vz(2  ,i  ) - Vz(1  ,i  )  ) * r20z

      !! bottom
      dzVx(nz-1) = (  Vx(nz,i  ) - Vx(nz-1,i  )  ) * r20z
      dzVx(nz  ) = (  0.0     - Vx(nz  ,i  )  ) * r20z
      dzVz(nz  ) = (  Vz(nz,i  ) - Vz(nz-1,i  )  ) * r20z


      !!
      !! i-edge
      !!
#ifdef _ES
      !NEC$ novector
#endif
      if( i == 1 ) then
        do k=kbeg_k, kend_k
          dxVx(k) = (  Vx(k  ,i  ) - 0.0       ) * r20x
          dxVz(k) = (  Vz(k  ,i+1) - Vz(k  ,i  )  ) * r20x
        end do
      else if ( i == 2 ) then
        do k=kbeg_k, kend_k
          dxVx(k) = (  Vx(k  ,i  ) - Vx(k  ,i-1)  ) * r20x
        end do
      else if ( i == nx-1 ) then
        do k=kbeg_k, kend_k
          dxVz(k) = (  Vz(k  ,i+1) - Vz(k  ,i  )  ) * r20x
        end do
      else if ( i==nx ) then
        do k=kbeg_k, kend_k
          dxVx(k) = (  Vx(k  ,i  ) - Vx(k  ,i-1)  ) * r20x
          dxVz(k) = (  0.0      - Vz(k  ,i  )  ) * r20x
        end do
      end if


      do k=kbeg_k, kend_k

        !!
        !! medium copy
        !!
        mu2    = 2*mu (k,i)
        lam2mu = lam(k,i) + mu2

        taup1 = taup(k,i)
        taus1 = taus(k,i)

        !!
        !! effective rigidity for shear stress components
        !!

        nnn = mu (k  ,i  )
        pnn = mu (k+1,i  )
        npn = mu (k,  i+1)
        ppn = mu (k+1,i+1)
        mu_xz = 4*nnn*pnn*npn*ppn / ( nnn*pnn*npn + nnn*pnn*ppn + nnn*npn*ppn + pnn*npn*ppn + epsl)


        !!
        !! update memory variables
        !!

        !! working variables for combinations of velocity derivatives
        d2v2      = dxVx(k) + dzVz(k)
        dxVx_dzVz = dxVx(k) + dzVz(k)
        dxVz_dzVx = dxVz(k) + dzVx(k)


        f_Rxx = lam2mu * taup1 * d2v2  -  mu2   * taus1 * dzVz(k)
        f_Rzz = lam2mu * taup1 * d2v2  -  mu2   * taus1 * dxVx(k)
        f_Rxz =                           mu_xz * taus1 * dxVz_dzVx


        Rxx_o = 0.0
        Rzz_o = 0.0
        Rxz_o = 0.0
        Rxx_n = 0.0
        Rzz_n = 0.0
        Rxz_n = 0.0
        if( nm > 0 ) then
          !! previous memory variables
          Rxx_o = sum( Rxx( 1:nm, k,i) )
          Rzz_o = sum( Rzz( 1:nm, k,i) )
          Rxz_o = sum( Rxz( 1:nm, k,i) )


          !! Crank-Nicolson Method for avoiding stiff solution
          do m=1, nm
            Rxx(m,k,i) = c1(m) * Rxx(m,k,i) - c2(m) * f_Rxx * dt
            Rzz(m,k,i) = c1(m) * Rzz(m,k,i) - c2(m) * f_Rzz * dt
            Rxz(m,k,i) = c1(m) * Rxz(m,k,i) - c2(m) * f_Rxz * dt
          end do

          !! new memory variables
          Rxx_n = sum( Rxx( 1:nm, k,i) )
          Rzz_n = sum( Rzz( 1:nm, k,i) )
          Rxz_n = sum( Rxz( 1:nm, k,i) )
        end if



        !!
        !! update stress components
        !!
        taup_plus1 = 1 + taup1
        taus_plus1 = 1 + taus1

        Sxx (k,i) = Sxx (k,i) + ( lam2mu*taup_plus1*d2v2 - mu2*taus_plus1*dzVz(k)   + ( Rxx_n+Rxx_o )/2 ) * dt
        Szz (k,i) = Szz (k,i) + ( lam2mu*taup_plus1*d2v2 - mu2*taus_plus1*dxVx(k)   + ( Rzz_n+Rzz_o )/2 ) * dt
        Sxz (k,i) = Sxz (k,i) + (                        mu_xz*taus_plus1*dxVz_dzVx + ( Rxz_n+Rxz_o )/2 ) * dt

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
  subroutine kernel__vmax( xmax, zmax )
    real(SP), intent(out) :: xmax, zmax
    integer :: i

    xmax = 0.0
    zmax = 0.0
    do i=ibeg_k, iend_k
      xmax = max( xmax, abs( vx(kob(i)+1,i) ) )
      zmax = max( zmax, abs( vz(kob(i)+1,i) ) )
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

    write(io) Vx (       kbeg_m:kend_m, ibeg_m:iend_m )
    write(io) Vz (       kbeg_m:kend_m, ibeg_m:iend_m )
    write(io) Sxx(       kbeg_m:kend_m, ibeg_m:iend_m )
    write(io) Szz(       kbeg_m:kend_m, ibeg_m:iend_m )
    write(io) Sxz(       kbeg_m:kend_m, ibeg_m:iend_m )

    do m=1, nm
      write(io) c1(m), c2(m)
      write(io) Rxx( m, kbeg_m:kend_m, ibeg_m:iend_m )
      write(io) Rzz( m, kbeg_m:kend_m, ibeg_m:iend_m )
      write(io) Rxz( m, kbeg_m:kend_m, ibeg_m:iend_m )
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
    read(io) Vx (       kbeg_m:kend_m, ibeg_m:iend_m )
    read(io) Vz (       kbeg_m:kend_m, ibeg_m:iend_m )
    read(io) Sxx(       kbeg_m:kend_m, ibeg_m:iend_m )
    read(io) Szz(       kbeg_m:kend_m, ibeg_m:iend_m )
    read(io) Sxz(       kbeg_m:kend_m, ibeg_m:iend_m )

    do m=1, nm
      read(io) c1(m), c2(m)
      read(io) Rxx( m, kbeg_m:kend_m, ibeg_m:iend_m )
      read(io) Rzz( m, kbeg_m:kend_m, ibeg_m:iend_m )
      read(io) Rxz( m, kbeg_m:kend_m, ibeg_m:iend_m )
    end do


  end subroutine kernel__restart
  !! --------------------------------------------------------------------------------------------------------------------------- !!



  !! --------------------------------------------------------------------------------------------------------------------------- !!
  subroutine memory_allocate
    !!
    !! memory allocation
    !!
    allocate(  Vx(       kbeg_m:kend_m, ibeg_m:iend_m ) )
    allocate(  Vz(       kbeg_m:kend_m, ibeg_m:iend_m ) )

    allocate( Sxx(       kbeg_m:kend_m, ibeg_m:iend_m ) )
    allocate( Szz(       kbeg_m:kend_m, ibeg_m:iend_m ) )
    allocate( Sxz(       kbeg_m:kend_m, ibeg_m:iend_m ) )

    if( nm > 0 ) then
      allocate( Rxx( 1:nm, kbeg_m:kend_m, ibeg_m:iend_m ) )
      allocate( Rzz( 1:nm, kbeg_m:kend_m, ibeg_m:iend_m ) )
      allocate( Rxz( 1:nm, kbeg_m:kend_m, ibeg_m:iend_m ) )
      allocate( c1(1:nm), c2(1:nm) )
    end if


  end subroutine memory_allocate
  !! --------------------------------------------------------------------------------------------------------------------------- !!


end module m_kernel
!! ----------------------------------------------------------------------------------------------------------------------------- !!
