!! ----------------------------------------------------------------------------------------------------------------------------- !!
!>
!! Set-up medium velocity/attenuation structure
!!
!! @copyright
!!   Copyright 2013-2018 Takuto Maeda. All rights reserved. This project is released under the MIT license.
!<
!! ----
#include "m_debug.h"
module m_medium

  use m_std
  use m_debug
  use m_global
  use m_pwatch
  use m_readini
  use m_vmodel_uni
  use m_vmodel_grd
  use m_vmodel_lhm
  use m_vmodel_user
  use m_vmodel_uni_rmed
  use m_vmodel_grd_rmed
  use m_vmodel_lhm_rmed
  implicit none
  private
  save

  public :: medium__setup
  public :: medium__initialized
  public :: medium__checkpoint
  public :: medium__restart

  logical :: init = .false.

contains

  !! --------------------------------------------------------------------------------------------------------------------------- !!
  !>
  !! Obtain elastic/anelastic medium structure
  !<
  !!
  subroutine medium__setup( io_prm )

    integer,      intent(in) :: io_prm

    real(SP) :: fq_min, fq_max                    !< frequency range
    real(SP) :: fq_ref                            !< reference frequency
    integer  :: i, j, k
    real(SP) :: zeta
    real(SP) :: vcut
    character(16) :: vmodel_type
    logical :: is_stabilize_pml

    call pwatch__on("medium__setup")

    !!
    !! allocate memory and initialize
    !!
    call memory_allocate()


    rho ( kbeg_m:kend_m, ibeg_m:iend_m, jbeg_m:jend_m ) = 0.0
    bx  ( kbeg_m:kend_m, ibeg_m:iend_m, jbeg_m:jend_m ) = 0.0
    by  ( kbeg_m:kend_m, ibeg_m:iend_m, jbeg_m:jend_m ) = 0.0
    bz  ( kbeg_m:kend_m, ibeg_m:iend_m, jbeg_m:jend_m ) = 0.0
    lam ( kbeg_m:kend_m, ibeg_m:iend_m, jbeg_m:jend_m ) = 0.0
    mu  ( kbeg_m:kend_m, ibeg_m:iend_m, jbeg_m:jend_m ) = 0.0
    muyz( kbeg_m:kend_m, ibeg_m:iend_m, jbeg_m:jend_m ) = 0.0
    muxz( kbeg_m:kend_m, ibeg_m:iend_m, jbeg_m:jend_m ) = 0.0
    muxy( kbeg_m:kend_m, ibeg_m:iend_m, jbeg_m:jend_m ) = 0.0
    taup( kbeg_m:kend_m, ibeg_m:iend_m, jbeg_m:jend_m ) = 0.0
    taus( kbeg_m:kend_m, ibeg_m:iend_m, jbeg_m:jend_m ) = 0.0


    !!
    !! benchmark mode: fixed medium parameter
    !!
    if( benchmark_mode ) then

      fq_min = 0.05
      fq_max = 5.0
      fq_ref = 1.0
      do k = kbeg_m, kend_m
        if( zc(k) < 0.0 ) then
          rho( k, :, : )  = 0.001
          mu ( k, :, : )  = 0.0
          lam( k, :, : )  = 0.0
        else
          rho(k,:,:)  = 2.7
          mu (k,:,:)  = 2.7 * 3.5*3.5
          lam(k,:,:)  = 2.7 * 3.5*3.5 ! poison solid: lambda = mu
        end if

        !! very large Q value (no attenuation) for benchmark
        taup(k,:,:) = 1e10
        taus(k,:,:) = 1e10
      end do

    else

      !!
      !! read parameters
      !!
      call readini( io_prm, 'fq_min', fq_min, 0.05 )
      call readini( io_prm, 'fq_max', fq_max, 5.00 )
      call readini( io_prm, 'fq_ref', fq_ref, 1.00 )

      call readini( io_prm, 'vmodel_type', vmodel_type, 'uni' )
      call readini( io_prm, 'vcut',  vcut, 0.0 )

      !!
      call pwatch__on("vmodel")
      select case ( trim(vmodel_type) )

      case ( 'user' )
        call vmodel_user( io_prm, ibeg_m, iend_m, jbeg_m, jend_m, kbeg_m, kend_m, xc, yc, zc, vcut, &
            rho, lam, mu, taup, taus, bddep )
      case ( 'uni' )
        call vmodel_uni( io_prm, ibeg_m, iend_m, jbeg_m, jend_m, kbeg_m, kend_m, xc, yc, zc, vcut, &
            rho, lam, mu, taup, taus, bddep )
      case ( 'grd' )
        call vmodel_grd( io_prm, ibeg_m, iend_m, jbeg_m, jend_m, kbeg_m, kend_m, xc, yc, zc, vcut, &
            rho, lam, mu, taup, taus, bddep )
      case ( 'lhm' )
        call vmodel_lhm( io_prm, ibeg_m, iend_m, jbeg_m, jend_m, kbeg_m, kend_m, xc, yc, zc, vcut, &
            rho, lam, mu, taup, taus, bddep )
      case ( 'uni_rmed' )
        call vmodel_uni_rmed( io_prm, ibeg_m, iend_m, jbeg_m, jend_m, kbeg_m, kend_m, xc, yc, zc, vcut, &
            rho, lam, mu, taup, taus, bddep )
      case ( 'grd_rmed' )
        call vmodel_grd_rmed( io_prm, ibeg_m, iend_m, jbeg_m, jend_m, kbeg_m, kend_m, xc, yc, zc, vcut, &
            rho, lam, mu, taup, taus, bddep )
      case ( 'lhm_rmed' )
        call vmodel_lhm_rmed( io_prm, ibeg_m, iend_m, jbeg_m, jend_m, kbeg_m, kend_m, xc, yc, zc, vcut, &
            rho, lam, mu, taup, taus, bddep )

      case default
        call assert( .false. )
      end select
      call pwatch__off("vmodel")

    end if


    !!
    !! homogenize absorber region
    !!
    do i=ibeg_m, na
      !$omp parallel do private(j,k)
      do j=jbeg_m,jend_m
        do k=kbeg_m, kend_m
          rho (k,i,j) = rho (k,na+1,j)
          lam (k,i,j) = lam (k,na+1,j)
          mu  (k,i,j) = mu  (k,na+1,j)
          taup(k,i,j) = taup(k,na+1,j)
          taus(k,i,j) = taus(k,na+1,j)
        end do
      end do
      !$omp end parallel do
    end do
    do i=nx-na+1,iend_m
      !$omp parallel do private(j,k)
      do j=jbeg_m,jend_m
        do k=kbeg_m, kend_m
          rho (k,i,j) = rho (k,nx-na,j)
          lam (k,i,j) = lam (k,nx-na,j)
          mu  (k,i,j) = mu  (k,nx-na,j)
          taup(k,i,j) = taup(k,nx-na,j)
          taus(k,i,j) = taus(k,nx-na,j)
        end do
      end do
      !$omp end parallel do
    end do
    do j=jbeg_m, na
      !$omp parallel do private(i,k)
      do i=ibeg_m,iend_m
        do k=kbeg_m, kend_m
          rho (k,i,j) = rho (k,i,na+1)
          lam (k,i,j) = lam (k,i,na+1)
          mu  (k,i,j) = mu  (k,i,na+1)
          taup(k,i,j) = taup(k,i,na+1)
          taus(k,i,j) = taus(k,i,na+1)
        end do
      end do
      !$omp end parallel do
    end do
    do j=ny-na+1,jend_m
      !$omp parallel do private(i,k)
      do i=ibeg_m,iend_m
        do k=kbeg_m, kend_m
          rho (k,i,j) = rho (k,i,ny-na)
          lam (k,i,j) = lam (k,i,ny-na)
          mu  (k,i,j) = mu  (k,i,ny-na)
          taup(k,i,j) = taup(k,i,ny-na)
          taus(k,i,j) = taus(k,i,ny-na)
        end do
      end do
      !$omp end parallel do
    end do

    !$omp parallel do private(i,j,k)
    do j=jbeg_m,jend_m
      do i=ibeg_m,iend_m
        do k=nz-na+1,kend_m
          rho (k,i,j) = rho (nz-na,i,j)
          lam (k,i,j) = lam (nz-na,i,j)
          mu  (k,i,j) = mu  (nz-na,i,j)
          taup(k,i,j) = taup(nz-na,i,j)
          taus(k,i,j) = taus(nz-na,i,j)
        end do
      end do
    end do
    !$omp end parallel do

    !!
    !! Define visco-elastic medium by tau-method
    !!
    call visco_set_relaxtime  ( nm, ts, fq_min, fq_max )
    zeta = visco_constq_zeta( nm, fq_min, fq_max, ts )
    if( benchmark_mode ) zeta = 0.0 !! no attenuation for benchmark mode

    !!
    !! Redefine taup and taus as relaxation times of P- and S-waves, based on tau-method
    !!
    !$omp parallel do private(i,j,k)
    do j=jbeg_m, jend_m
      do i=ibeg_m, iend_m
        do k=kbeg_m, kend_m
          taup(k,i,j) = nm * zeta / taup(k,i,j)
          taus(k,i,j) = nm * zeta / taus(k,i,j)
        end do
      end do
    end do
    !$omp end parallel do

    !!
    !! Scaling lambda and mu based on physical dispersion & reference frequency
    !!
    call relaxed_medium()

    !!
    !! Free surface & Ocean bottom
    !!
    call surface_detection()

    !!
    !! global min/max of velocity for stability check
    !!
    call velocity_minmax()

    !!
    !! modify low-velocity layers in PML region
    !!
    call readini( io_prm, 'stabilize_pml', is_stabilize_pml, .false. )
    if( is_stabilize_pml ) then
      call stabilize_absorber()
    end if


    !!
    !! Store grid-boundary averaged medium for staggered grid locations
    !!
    call averaged_medium()


    !!
    !! initialized flag
    !!
    init = .true.

    call pwatch__off("medium__setup")


  contains

    !! ------------------------------------------------------------------------------------------------------------------------ !!
    !! scale medium velocity using reference frequency
    !!
    subroutine relaxed_medium()
      
      integer :: i, j, k, im
      real(SP) :: rho_beta2, rho_alpha2
      real(SP) :: chi_mu, chi_lam
      complex(SP) :: cc
      real(SP) :: omega

      !! mu, lam must be re-defined including sleeve area for medium smoothing

      if( nm == 0 ) return

      omega = 2 * PI * fq_ref
      cc = 0.0
      do im=1, nm
        cc = cc + ( EI * omega * ts(im) ) / ( 1.0 - EI * omega * ts(im) )
      end do
      cc = cc / nm

      !$omp parallel do private(rho_beta2, rho_alpha2, chi_mu, chi_lam, i, j, k)
      do j=jbeg_m, jend_m
        do i=ibeg_m, iend_m
          do k=kbeg_m, kend_m

            rho_beta2  =  mu(k,i,j)
            rho_alpha2 = lam(k,i,j) + 2*mu(k,i,j)

            chi_mu  = 1.0 / real( 1.0 / sqrt( 1.0 - cc * taus(k,i,j) ) )
            chi_lam = 1.0 / real( 1.0 / sqrt( 1.0 - cc * taup(k,i,j) ) )

            mu (k,i,j) = rho_beta2  / chi_mu **2
            lam(k,i,j) = rho_alpha2 / chi_lam**2 - 2 * mu(k,i,j)

          end do
        end do
      end do
      !$omp end parallel do

    end subroutine relaxed_medium
    !! ------------------------------------------------------------------------------------------------------------------------ !!

    !! ------------------------------------------------------------------------------------------------------------------------ !!
    !>
    !! Avoid low-velocity layer for stabilize PML absorber
    !<
    !! --
    subroutine stabilize_absorber()

      integer :: i, j, k, k2
      real :: vp, vs, gamma
      real, parameter :: V_DYNAMIC_RANGE = 0.4 ! ratio between maximum and minimum velocity
      real :: vmin_pml
      integer :: LV_THICK = 20 !! minimum thickness of low-velocity layer in grids

      vmin_pml = vmax * V_DYNAMIC_RANGE

      do j=jbeg-1, jend+1
        do i=ibeg-1, iend+1
          k=minval(kbeg_a(i-2:i+2,j-2:j+2))
          do while( k<=kend )

            if( lam(k,i,j) < lam(k-1,i,j) .or. mu(k,i,j) < mu(k-1,i,j) ) then

              do k2=k+1, kend
                if( lam(k2,i,j) > lam(k2-1,i,j) .or. mu(k2,i,j) > mu(k2-1,i,j) ) exit
              end do

              if( k2-k <= LV_THICK ) then


                rho (k,i,j) = rho (k-1,i,j)
                lam (k,i,j) = lam (k-1,i,j)
                mu  (k,i,j) = mu  (k-1,i,j)
                taup(k,i,j) = taup(k-1,i,j)
                taus(k,i,j) = taus(k-1,i,j)
                k = k2 - 1

              end if

            end if

            k = k + 1

          end do
        end do
      end do

      do j=jbeg-1, jend+1
        do i=ibeg-1, iend+1
          do k=minval(kbeg_a(i-2:i+2,j-2:j+2)), kend

            vp = sqrt( (lam(k,i,j) + 2 * mu(k,i,j))/rho(k,i,j) )
            vs = sqrt( mu(k,i,j) / rho(k,i,j) )

            ! skip ocean and air
            if( vs < epsilon(1.0) ) cycle

            gamma = sqrt(3.0)
            if( vs < vmin_pml ) then
              vs = vmin_pml
              vp = vs * gamma

              lam(k,i,j) = rho(k,i,j) * (vp**2 - 2 * vs**2)
              mu (k,i,j) = rho(k,i,j) * (vs**2)
            end if
          end do
        end do
      end do

    end subroutine stabilize_absorber
    !! ------------------------------------------------------------------------------------------------------------------------ !!


    !! ------------------------------------------------------------------------------------------------------------------------ !!
    !! free surface & ocean bottom boundaries
    !!
    subroutine surface_detection

      real(SP) :: epsl = epsilon(1.0)
      integer  :: i, j, k

      !!
      !! initial value
      !! This initial settings implies there is NO free surface boundary in the interior medium
      !!
      kfs(:,:) = kbeg-1
      kob(:,:) = kbeg-1


      !!
      !! Detect free surfaces as material interface
      !! kfs, kob must be defined ibeg-1:iend+2, jbeg-1:jend+2 to detect (kfs|kob)_top/bot .
      !!
      !$omp parallel do private(i,j,k)
      do j=jbeg-1, jend+2
        do i=ibeg-1, iend+2
          do k=kbeg, kend-1

            !! air(ocean)-to-solid boundary
            if( abs(mu (k,i,j)) < epsl .and. abs(mu (k+1,i,j)) > epsl ) then
              kob(i,j) = k
            end if

            !! air-to-solid(ocean) boundary
            if( abs(lam(k,i,j)) < epsl .and. abs(lam(k+1,i,j)) > epsl ) then
              kfs(i,j) = k
            end if

          end do
        end do
      end do
      !$omp end parallel do

      !!
      !! define 2nd-order derivative area #2013-00419
      !! updated (stable) version: 2015-08-18
      !!
      !$omp parallel do private(i,j,k)
      do j=jbeg, jend
        do i=ibeg, iend

          kfs_top(i,j) = max( minval( kfs(i-2:i+3,j-2:j+3) ) - 2, kbeg )
          kfs_bot(i,j) = min( maxval( kfs(i-2:i+3,j-2:j+3) ) + 2, kend )

          kob_top(i,j) = max( minval( kob(i-2:i+3,j-2:j+3) ) - 2, kbeg )
          kob_bot(i,j) = min( maxval( kob(i-2:i+3,j-2:j+3) ) + 2, kend )

        end do
      end do
      !$omp end parallel do

    end subroutine surface_detection
    !! ------------------------------------------------------------------------------------------------------------------------ !!

    !! ------------------------------------------------------------------------------------------------------------------------ !!
    !! maximum & minimum velocities
    !!
    subroutine velocity_minmax()
      real(SP) :: vmin1, vmax1
      real(SP) :: vp, vs
      integer  :: i, j, k
      integer  :: ierr

      !!
      !!
      vmax1 = -1
      vmin1 = 1e30
      !$omp parallel do reduction(max:vmax1) reduction(min:vmin1)
      do j=jbeg, jend
        do i=ibeg, iend
          do k=kfs(i,j)+1, kend
            vp = sqrt( ( lam(k,i,j) + 2 * mu(k,i,j) ) / rho(k,i,j) )
            vs = sqrt(                    mu(k,i,j)   / rho(k,i,j) )

            vmax1 = max( vmax1, vp )
            if( vs < epsilon(1.0) ) cycle
            vmin1 = min( vmin1, vs )

          end do
        end do
      end do
      !$omp end parallel do

      !! obtain global maximum and minimum velocities
      call mpi_allreduce( vmax1, vmax, 1, MPI_REAL, MPI_MAX, mpi_comm_world, ierr )
      call mpi_allreduce( vmin1, vmin, 1, MPI_REAL, MPI_MIN, mpi_comm_world, ierr )

    end subroutine velocity_minmax
    !! ------------------------------------------------------------------------------------------------------------------------ !!



    !! ------------------------------------------------------------------------------------------------------------------------ !!
    subroutine averaged_medium

      real(SP) :: nnn, pnn, npn, ppn, nnp, npp, pnp
      real(SP) :: epsl = epsilon(1.0)
      integer  :: i, j, k

      !$omp parallel do private( nnn,pnn,npn,ppn,nnp,npp,pnp, i,j,k )
      do j=jbeg, jend
        do i=ibeg, iend
          do k=kbeg, kend
            bx(k,i,j) = 2.0 / ( rho(k,i,j) + rho(k,i+1,j) )
            by(k,i,j) = 2.0 / ( rho(k,i,j) + rho(k,i,j+1) )
            bz(k,i,j) = 2.0 / ( rho(k,i,j) + rho(k+1,i,j) )

            nnn = mu(k  ,i  ,j  )
            pnn = mu(k+1,i,  j  )
            npn = mu(k,  i+1,j  )
            ppn = mu(k+1,i+1,j  )
            nnp = mu(k,  i,  j+1)
            npp = mu(k,  i+1,j+1)
            pnp = mu(k+1,i,  j+1)

            muxz(k,i,j) = 4*nnn*pnn*npn*ppn / ( nnn*pnn*npn + nnn*pnn*ppn + nnn*npn*ppn + pnn*npn*ppn + epsl )
            muxy(k,i,j) = 4*nnn*npn*nnp*npp / ( nnn*npn*nnp + nnn*npn*npp + nnn*nnp*npp + npn*nnp*npp + epsl )
            muyz(k,i,j) = 4*nnn*pnn*nnp*pnp / ( nnn*pnn*nnp + nnn*pnn*pnp + nnn*nnp*pnp + pnn*nnp*pnp + epsl )
          end do
        end do
      end do
      !$omp end parallel do
    end subroutine averaged_medium
    !! ------------------------------------------------------------------------------------------------------------------------ !!


  end subroutine medium__setup
  !! --------------------------------------------------------------------------------------------------------------------------- !!

  !! --------------------------------------------------------------------------------------------------------------------------- !!
  !>
  !! Check if medium__setup has already been called
  !<
  !!
  logical function medium__initialized()

    medium__initialized = init

  end function medium__initialized
  !! --------------------------------------------------------------------------------------------------------------------------- !!

  !! --------------------------------------------------------------------------------------------------------------------------- !!
  subroutine memory_allocate()
    !!
    allocate( rho    (  kbeg_m:kend_m, ibeg_m:iend_m, jbeg_m:jend_m ) )
    allocate( bx     (  kbeg_m:kend_m, ibeg_m:iend_m, jbeg_m:jend_m ) )
    allocate( by     (  kbeg_m:kend_m, ibeg_m:iend_m, jbeg_m:jend_m ) )
    allocate( bz     (  kbeg_m:kend_m, ibeg_m:iend_m, jbeg_m:jend_m ) )
    allocate( lam    (  kbeg_m:kend_m, ibeg_m:iend_m, jbeg_m:jend_m ) )
    allocate( mu     (  kbeg_m:kend_m, ibeg_m:iend_m, jbeg_m:jend_m ) )
    allocate( muyz   (  kbeg_m:kend_m, ibeg_m:iend_m, jbeg_m:jend_m ) )
    allocate( muxz   (  kbeg_m:kend_m, ibeg_m:iend_m, jbeg_m:jend_m ) )
    allocate( muxy   (  kbeg_m:kend_m, ibeg_m:iend_m, jbeg_m:jend_m ) )
    allocate( taup   (  kbeg_m:kend_m, ibeg_m:iend_m, jbeg_m:jend_m ) )
    allocate( taus   (  kbeg_m:kend_m, ibeg_m:iend_m, jbeg_m:jend_m ) )
    allocate( kfs    (                 ibeg_m:iend_m, jbeg_m:jend_m ) )
    allocate( kob    (                 ibeg_m:iend_m, jbeg_m:jend_m ) )
    allocate( kfs_top(                 ibeg_m:iend_m, jbeg_m:jend_m ) )
    allocate( kfs_bot(                 ibeg_m:iend_m, jbeg_m:jend_m ) )
    allocate( kob_top(                 ibeg_m:iend_m, jbeg_m:jend_m ) )
    allocate( kob_bot(                 ibeg_m:iend_m, jbeg_m:jend_m ) )
    allocate( bddep  (                 ibeg_m:iend_m, jbeg_m:jend_m, 0:NBD ) ) !! 0: topo, 1-2: discontinuity
    if( nm > 0 ) allocate( ts(  1:nm ) )

  end subroutine memory_allocate
  !! --------------------------------------------------------------------------------------------------------------------------- !!

  !! --------------------------------------------------------------------------------------------------------------------------- !!
  subroutine medium__checkpoint( io )

    integer, intent(in) :: io
    integer :: j
    !! ----

    do j=jbeg_m,jend_m;  write(io)   bx(kbeg_m:kend_m,ibeg_m:iend_m,j); end do;  deallocate(bx)
    do j=jbeg_m,jend_m;  write(io)   by(kbeg_m:kend_m,ibeg_m:iend_m,j); end do;  deallocate(by)
    do j=jbeg_m,jend_m;  write(io)   bz(kbeg_m:kend_m,ibeg_m:iend_m,j); end do;  deallocate(bz)
    do j=jbeg_m,jend_m;  write(io)  lam(kbeg_m:kend_m,ibeg_m:iend_m,j); end do;  deallocate(lam)
    do j=jbeg_m,jend_m;  write(io)   mu(kbeg_m:kend_m,ibeg_m:iend_m,j); end do;  deallocate(mu)
    do j=jbeg_m,jend_m;  write(io) muyz(kbeg_m:kend_m,ibeg_m:iend_m,j); end do;  deallocate(muyz)
    do j=jbeg_m,jend_m;  write(io) muxz(kbeg_m:kend_m,ibeg_m:iend_m,j); end do;  deallocate(muxz)
    do j=jbeg_m,jend_m;  write(io) muxy(kbeg_m:kend_m,ibeg_m:iend_m,j); end do;  deallocate(muxy)
    do j=jbeg_m,jend_m;  write(io) taup(kbeg_m:kend_m,ibeg_m:iend_m,j); end do;  deallocate(taup)
    do j=jbeg_m,jend_m;  write(io) taus(kbeg_m:kend_m,ibeg_m:iend_m,j); end do;  deallocate(taus)

    write(io) kfs    ( ibeg_m:iend_m, jbeg_m:jend_m )
    write(io) kob    ( ibeg_m:iend_m, jbeg_m:jend_m )
    write(io) kfs_top( ibeg_m:iend_m, jbeg_m:jend_m )
    write(io) kfs_bot( ibeg_m:iend_m, jbeg_m:jend_m )
    write(io) kob_top( ibeg_m:iend_m, jbeg_m:jend_m )
    write(io) kob_bot( ibeg_m:iend_m, jbeg_m:jend_m )
    write(io) bddep  ( ibeg_m:iend_m, jbeg_m:jend_m, 0:NBD )
    if( nm > 0 ) write(io) ts(  1:nm )
    
    deallocate( kfs, kob, kfs_top, kfs_bot, kob_top, kob_bot, bddep )
    
  end subroutine medium__checkpoint
  !! --------------------------------------------------------------------------------------------------------------------------- !!

  !! --------------------------------------------------------------------------------------------------------------------------- !!
  subroutine medium__restart( io )
    
    integer, intent(in) :: io
    integer :: j
    !! ----
    
    call memory_allocate()
    do j=jbeg_m,jend_m;  read(io)   bx(kbeg_m:kend_m,ibeg_m:iend_m,j); end do;
    do j=jbeg_m,jend_m;  read(io)   by(kbeg_m:kend_m,ibeg_m:iend_m,j); end do;
    do j=jbeg_m,jend_m;  read(io)   bz(kbeg_m:kend_m,ibeg_m:iend_m,j); end do;
    do j=jbeg_m,jend_m;  read(io)  lam(kbeg_m:kend_m,ibeg_m:iend_m,j); end do;
    do j=jbeg_m,jend_m;  read(io)   mu(kbeg_m:kend_m,ibeg_m:iend_m,j); end do;
    do j=jbeg_m,jend_m;  read(io) muyz(kbeg_m:kend_m,ibeg_m:iend_m,j); end do;
    do j=jbeg_m,jend_m;  read(io) muxz(kbeg_m:kend_m,ibeg_m:iend_m,j); end do;
    do j=jbeg_m,jend_m;  read(io) muxy(kbeg_m:kend_m,ibeg_m:iend_m,j); end do;
    do j=jbeg_m,jend_m;  read(io) taup(kbeg_m:kend_m,ibeg_m:iend_m,j); end do;
    do j=jbeg_m,jend_m;  read(io) taus(kbeg_m:kend_m,ibeg_m:iend_m,j); end do;

    read(io) kfs    ( ibeg_m:iend_m, jbeg_m:jend_m )
    read(io) kob    ( ibeg_m:iend_m, jbeg_m:jend_m )
    read(io) kfs_top( ibeg_m:iend_m, jbeg_m:jend_m )
    read(io) kfs_bot( ibeg_m:iend_m, jbeg_m:jend_m )
    read(io) kob_top( ibeg_m:iend_m, jbeg_m:jend_m )
    read(io) kob_bot( ibeg_m:iend_m, jbeg_m:jend_m )
    read(io) bddep  ( ibeg_m:iend_m, jbeg_m:jend_m, 0:NBD )
    
    if( nm > 0 ) read(io) ts(  1:nm )
    
  end subroutine medium__restart
  !! --------------------------------------------------------------------------------------------------------------------------- !!
  
  
end module m_medium
!! ----------------------------------------------------------------------------------------------------------------------------- !!
  
