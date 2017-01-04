!! ----------------------------------------------------------------------------------------------------------------------------- !!
!>
!! User-routines for defining velocity/attenuation structure
!!
!! @copyright
!!   Copyright 2013-2016 Takuto Maeda. All rights reserved. This project is released under the MIT license.
!<
!! ----
#include "m_debug.h"
module m_vmodel_uni_rmed

  use m_std
  use m_debug
  use m_readini
  use m_global
  use m_rdrmed
  use m_fdtool
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
  subroutine vmodel_uni_rmed( io_prm, i0, i1, j0, j1, k0, k1, xc, yc, zc, vcut, rho, lam, mu, Qp, Qs, bd )

    !! -- Arguments
    integer,  intent(in)  :: io_prm
    integer,  intent(in)  :: i0, i1                         !< i-region
    integer,  intent(in)  :: j0, j1                         !< j-region
    integer,  intent(in)  :: k0, k1                         !< k-region
    real(SP), intent(in)  :: xc  ( i0:i1 )                  !< x-coordinate location
    real(SP), intent(in)  :: yc  ( j0:j1 )                  !< y-coordinate location
    real(SP), intent(in)  :: zc  ( k0:k1 )                  !< z-coordinate location
    real(SP), intent(in)  :: vcut                           !< cut-off minimum velocity
    real(SP), intent(out) :: rho ( k0:k1, i0:i1, j0:j1 )    !< mass density [g/cm^3]
    real(SP), intent(out) :: lam ( k0:k1, i0:i1, j0:j1 )    !< Lame's parameter lambda [ (g/cm^3) * (km/s) ]
    real(SP), intent(out) :: mu  ( k0:k1, i0:i1, j0:j1 )    !< Lame's parameter mu     [ (g/cm^3) * (km/s) ]
    real(SP), intent(out) :: qp  ( k0:k1, i0:i1, j0:j1 )    !< P-wave attenuation
    real(SP), intent(out) :: qs  ( k0:k1, i0:i1, j0:j1 )    !< S-wave attenuation
    real(SP), intent(out) :: bd  ( i0:i1, j0:j1, 0:NBD )    !< Boundary depths
    !! --

    integer  :: i, j, k
    real(SP) :: vp0, vs0, rho0, qp0, qs0, topo0, rho2, vp2, vs2
    real(SP) :: vp1, vs1
    real(SP) :: dum
    character(256) :: dir_rmed, fn_rmed
    real(SP), allocatable :: xi(:,:,:)
    logical :: is_exist
    real(SP) :: vmin, vmax, dh, cc, rhomin
    logical  :: vmax_over, vmin_under, rhomin_under
    !! ----

    call readini( io_prm, 'vp0',    vp0, 5.0 )
    call readini( io_prm, 'vs0',    vs0, vp0/sqrt(3.0) )
    call readini( io_prm, 'rho0',   rho0, 2.7 )
    call readini( io_prm, 'qp0',    qp0, 1000000.0 )
    call readini( io_prm, 'qs0',    qs0, 1000000.0 )
    call readini( io_prm, 'topo0', topo0, 0.0 )
    call readini( io_prm, 'rhomin', rhomin, 1.0 )

    vmin = vcut

    dh = sqrt(3.) / sqrt( 1./dx**2 + 1./dy**2 + 1./dz**2 )
    cc = 6. / 7. !! assume 4th order
    vmax = 0.95 * cc * dh / dt  ! 0.95 is a safety coefficient

    vmax_over  = .false.
    vmin_under = .false.
    rhomin_under = .false.

    !! random medium
    call readini( io_prm, 'dir_rmed', dir_rmed, '' )
    call readini( io_prm, 'fn_rmed0', fn_rmed, '' )
    fn_rmed = trim(dir_rmed) // '/' // trim(fn_rmed)
    allocate( xi(k0:k1,i0:i1,j0:j1) )
    inquire( file=fn_rmed, exist=is_exist )
    if( is_exist ) then
       call rdrmed__3d( i0, i1, j0, j1, k0, k1, fn_rmed, xi )
    else
       call info( 'rmedia file '//trim(fn_rmed)//' not found' )
       xi(:,:,:) = 0.0
    end if


    do j = j0, j1
       do i = i0, i1

          !! define topography shape here
          bd(i,j,0) = topo0

          do k = k0, k1

             if( zc( k ) > bd(i,j,0) ) then

                rho2 = ( 1.0 + 0.8 * xi(k,i,j) ) * rho0
                vp2  = ( 1.0 +       xi(k,i,j) ) * vp0
                vs2  = ( 1.0 +       xi(k,i,j) ) * vs0

                call vcheck( vp2, vs2, rho2, xi(k,i,j), vmin, vmax, rhomin, vmin_under, vmax_over, rhomin_under )

                !! elastic medium
                rho(k,i,j) = rho2
                mu (k,i,j) = rho(k,i,j) * vs2 * vs2
                lam(k,i,j) = rho(k,i,j) * ( vp2*vp2 - 2*vs2*vs2 )
                qp (k,i,j) = qp0
                qs (k,i,j) = qs0

             else if ( zc (k) > 0.0 ) then

                !! ocean column

                vp1 = 1.5
                vs1 = 0.0

                rho(k,i,j) = 1.0
                mu (k,i,j) = rho(k,i,j) * vs1 * vs1
                lam(k,i,j) = rho(k,i,j) * ( vp1*vp1 - 2*vs1*vs1 )
                qp (k,i,j) = 1000000.0 ! effectively no attenuation in ocean column
                qs (k,i,j) = 1000000.0

             else

                !! air column

                vp1 = 0.0
                vs1 = 0.0

                rho(k,i,j) = 0.001
                mu (k,i,j) = rho(k,i,j) * vs1 * vs1
                lam(k,i,j) = rho(k,i,j) * ( vp1*vp1 - 2*vs1*vs1 )
                qp (k,i,j) = 10.0 ! artificially strong attenuation in air-column
                qs (k,i,j) = 10.0 ! artificially strong attenuation in air-column

             end if
          end do
       end do
    end do

    !! notification for velocity torelance
    if( vmax_over  ) call info( 'Too high velocity due to random media was corrected. ')
    if( vmin_under ) call info( 'Too low  velocity due to random media was corrected. ')
    if( rhomin_under ) call info( 'Too low  density due to random media was corrected. ')

    ! dummy value
    bd(:,:,1:NBD) = -9999

    ! substitute to a dummy variable for avoiding compiler warnings
    dum = xc(i0)
    dum = yc(j0)
    dum = zc(k0)
    dum = vcut

  end subroutine vmodel_uni_rmed
  !! --------------------------------------------------------------------------------------------------------------------------- !!

end module m_vmodel_uni_rmed
!! ----------------------------------------------------------------------------------------------------------------------------- !!
