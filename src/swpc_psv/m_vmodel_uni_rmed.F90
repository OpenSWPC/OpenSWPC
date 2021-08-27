!! ----------------------------------------------------------------------------------------------------------------------------- !!
!>
!! User-routines for defining velocity/attenuation structure
!!
!! @copyright
!!   Copyright 2013-2021 Takuto Maeda. All rights reserved. This project is released under the MIT license.
!<
!! ----

#include "m_debug.h"
module m_vmodel_uni_rmed

  use m_std
  use m_debug
  use m_readini
  use m_global
  use m_rdrmed
  use m_seawater
  implicit none
  private
  save

  public :: vmodel_uni_rmed

contains

  !! --------------------------------------------------------------------------------------------------------------------------- !!
  !>
  !! Define meidum velocity, density and attenuation
  !<
  !! ----
  subroutine vmodel_uni_rmed( io_prm, i0, i1, k0, k1, xc, zc, vcut, rho, lam, mu, Qp, Qs, bd )

    !! -- Arguments
    integer,  intent(in)  :: io_prm
    integer,  intent(in)  :: i0, i1                  !< i-region
    integer,  intent(in)  :: k0, k1                  !< k-region
    real(SP), intent(in)  :: xc  ( i0:i1 )           !< x-coordinate location
    real(SP), intent(in)  :: zc  ( k0:k1 )           !< z-coordinate location
    real(SP), intent(in)  :: vcut                    !< cut-off minimum velocity
    real(SP), intent(out) :: rho ( k0:k1, i0:i1 )    !< mass density [g/cm^3]
    real(SP), intent(out) :: lam ( k0:k1, i0:i1 )    !< Lame's parameter lambda [ (g/cm^3) * (km/s) ]
    real(SP), intent(out) :: mu  ( k0:k1, i0:i1 )    !< Lame's parameter mu     [ (g/cm^3) * (km/s) ]
    real(SP), intent(out) :: qp  ( k0:k1, i0:i1 )    !< P-wave attenuation
    real(SP), intent(out) :: qs  ( k0:k1, i0:i1 )    !< S-wave attenuation
    real(SP), intent(out) :: bd  ( i0:i1, 0:NBD )    !< Boundary depths
    !! --

    integer  :: i, k
    real(SP) :: vp0, vs0, rho0, qp0, qs0, topo0
    real(SP) :: vp1, vs1
    real(SP) :: vp2, vs2, rho2
    real(SP) :: dum
    character(256) :: fn_rmed0
    real(SP), allocatable :: xi(:,:)
    logical :: is_exist
    character(256) :: dir_rmed
    real(SP) :: vmin, vmax, dh, cc, rhomin
    logical  :: is_vmax_over, is_vmin_under, is_rhomin_under
    logical :: use_munk
    logical :: earth_flattening
    real(SP) :: zs(k0:k1) ! spherical depth for earth_flattening
    real(SP) :: Cv(k0:k1) ! velocity scaling coefficient for earth_flattening    
    !! ----

    call readini( io_prm, 'vp0',    vp0, 5.0 )
    call readini( io_prm, 'vs0',    vs0, vp0/sqrt(3.0) )
    call readini( io_prm, 'rho0',   rho0, 2.7 )
    call readini( io_prm, 'qp0',    qp0, 1000000.0 )
    call readini( io_prm, 'qs0',    qs0, 1000000.0 )
    call readini( io_prm, 'topo0',  topo0, 0.0 )
    call readini( io_prm, 'rhomin', rhomin, 1.0 )

    !! seawater
    call readini( io_prm, 'munk_profile', use_munk, .false. )
    call seawater__init( use_munk )

    !! earth-flattening transformation
    call readini( io_prm, 'earth_flattening', earth_flattening, .false. )
    if( earth_flattening ) then
      do k=k0, k1
        zs(k) = R_EARTH - R_EARTH * exp( - zc(k) / R_EARTH )
        Cv(k) = exp( zc(k) / R_EARTH)
      end do
    else
      zs(:) = zc(:)
      Cv(:) = 1.0
    end if    

    vmin = vcut

    dh = 1. / sqrt( 1./dx**2 + 1./dz**2 )
    cc = 6. / 7. !! assume 4th order
    vmax = cc * dh / dt

    call readini( io_prm, 'dir_rmed', dir_rmed, '' )
    call readini( io_prm, 'fn_rmed0', fn_rmed0, '' )
    fn_rmed0 = trim(dir_rmed) // '/' // trim(fn_rmed0)
    allocate( xi(k0:k1,i0:i1) )
    inquire( file=fn_rmed0, exist=is_exist )
    if( is_exist ) then
      call rdrmed__2d( i0, i1, k0, k1, fn_rmed0, xi )
    else
      call info( 'rmedia file '//trim(fn_rmed0)//' not found' )
      xi(:,:) = 0.0
    end if

    is_vmax_over  = .false.
    is_vmin_under = .false.
    is_rhomin_under = .false.

    do i = i0, i1

      bd(i,0) = topo0

      do k = k0, k1

        if( zs( k ) > bd(i,0) ) then

          !! elastic medium
          rho2 = (1 + 0.8*xi(k,i)) * rho0
          vp2  = (1 +     xi(k,i)) * Cv(k)*vp0
          vs2  = (1 +     xi(k,i)) * Cv(k)*vs0

          call vcheck( vp2, vs2, rho2, xi(k,i), vmin, vmax, rhomin, is_vmin_under, is_vmax_over, is_rhomin_under )

          rho(k,i) = rho2
          mu (k,i) = rho2 * vs2 * vs2
          lam(k,i) = rho2 * ( vp2*vp2 - 2*vs2*vs2 )
          qp (k,i) = qp0
          qs (k,i) = qs0

        else if ( zs (k) > 0.0 ) then

          !! ocean column

          vp1 = Cv(k)*seawater__vel(zc(k))
          vs1 = 0.0

          rho(k,i) = 1.0
          mu (k,i) = rho(k,i) * vs1 * vs1
          lam(k,i) = rho(k,i) * ( vp1*vp1 - 2*vs1*vs1 )
          qp (k,i) = 1000000.0 ! effectively no attenuation in ocean column
          qs (k,i) = 1000000.0

        else

          !! air column

          vp1 = 0.0
          vs1 = 0.0

          rho(k,i) = 0.001
          mu (k,i) = rho(k,i) * vs1 * vs1
          lam(k,i) = rho(k,i) * ( vp1*vp1 - 2*vs1*vs1 )
          qp (k,i) = 10.0 ! artificially strong attenuation in air-column
          qs (k,i) = 10.0 ! artificially strong attenuation in air-column

        end if
      end do
    end do

    !! notification for velocity torelance
    if( is_vmax_over  ) call info( 'Too high velocity due to random media was corrected. ')
    if( is_vmin_under ) call info( 'Too low  velocity due to random media was corrected. ')
    if( is_rhomin_under ) call info( 'Too low  density due to random media was corrected. ')

    ! dummy value
    bd(:,1:NBD) = -9999
    dum = xc(i0)
    dum = zc(k0)
    dum = vcut

  end subroutine vmodel_uni_rmed
  !! --------------------------------------------------------------------------------------------------------------------------- !!

end module m_vmodel_uni_rmed
!! ----------------------------------------------------------------------------------------------------------------------------- !!
