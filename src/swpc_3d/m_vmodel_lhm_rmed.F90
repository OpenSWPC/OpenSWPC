!! ----------------------------------------------------------------------------------------------------------------------------- !!
!>
!! 1D velocity structure
!!
!! @copyright
!!   Copyright 2013-2017 Takuto Maeda. All rights reserved. This project is released under the MIT license.
!<
!! ----
#include "m_debug.h"
module m_vmodel_lhm_rmed

  use m_std
  use m_debug
  use m_readini
  use m_global
  use m_rdrmed
  use m_fdtool
  implicit none
  private
  save

  public :: vmodel_lhm_rmed

contains

  !! --------------------------------------------------------------------------------------------------------------------------- !!
  !>
  !! Define meidum velocity, density and attenuation
  !<
  !! ----
  subroutine vmodel_lhm_rmed( io_prm, i0, i1, j0, j1, k0, k1, xc, yc, zc, vcut, rho, lam, mu, Qp, Qs, bd )

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
    character(256) :: fn_lhm
    integer  :: i, j, k, l
    real(SP), allocatable, dimension(:) :: vp0, vs0, rho0, qp0, qs0, depth
    real(SP) :: vp1, vs1, rho1, qp1, qs1
    real(SP) :: dum
    integer :: ierr
    integer :: io_vel
    logical :: is_exist
    integer :: nlayer
    character(256) :: adum
    character(256), allocatable :: fn_rmed(:)
    real(SP), allocatable :: xi(:,:,:,:)
    integer :: n_rmed
    integer, allocatable :: tbl_rmed(:)
    character(256), allocatable :: fn_rmed2(:)
    character(256) :: dir_rmed
    real(SP) :: vmin, vmax, dh, cc, rhomin
    logical  :: vmax_over, vmin_under, rhomin_under
    !! ----

    call readini( io_prm, 'fn_lhm_rmed', fn_lhm, '' )
    call readini( io_prm,  'dir_rmed', dir_rmed, '' )
    call readini( io_prm, 'rhomin', rhomin, 1.0 )

    inquire( file=fn_lhm, exist=is_exist )
    call assert( is_exist )

    call std__getio( io_vel )
    open( io_vel, file=fn_lhm, status='old', action='read', iostat=ierr )
    call assert( ierr == 0 )
    call std__countline( io_vel, nlayer, '#' )
    allocate( depth(nlayer), rho0(nlayer), vp0(nlayer), vs0(nlayer), qp0(nlayer), qs0(nlayer), fn_rmed(nlayer) )


    vmin = vcut

    dh = 1. / sqrt( 1./dx**2 + 1./dy**2 + 1./dz**2 )
    cc = 6. / 7. !! assume 4th order
    vmax = 0.95 * cc * dh / dt ! 0.95 is a safety coefficient

    vmax_over  = .false.
    vmin_under = .false.
    rhomin_under = .false.

    l = 0
    do
      read(io_vel,'(A256)', iostat=ierr ) adum
      if( ierr /= 0 ) exit
      adum = trim(adjustl(adum))
      if( trim(adum) == '' ) cycle
      if( adum(1:1) == "#" ) cycle
      l = l + 1
      read(adum,*) depth(l), rho0(l), vp0(l), vs0(l), qp0(l), qs0(l), fn_rmed(l)
    end do
    close( io_vel )

    !! velocity cut-off
    do l = nlayer-1, 1, -1
      if( ( vp0(l) < vcut .or. vs0(l) < vcut ) .and. ( vp0(l) > 0 .and. vs0(l) >0 ) ) then
        vp0(l)  = vp0(l+1)
        vs0(l)  = vs0(l+1)
        rho0(l) = rho0(l+1)
        qp0(l)  = qp0(l+1)
        qs0(l)  = qs0(l+1)
      end if
    end do

    do l=1, nlayer
      fn_rmed(l) = trim( dir_rmed ) // '/' // trim(fn_rmed(l) )
    end do

    !! Read random media
    allocate( tbl_rmed(nlayer), fn_rmed2(nlayer) )
    call independent_list( nlayer, fn_rmed, n_rmed, tbl_rmed, fn_rmed2 )

    allocate(xi(k0:k1,i0:i1,j0:j1,n_rmed) )
    do l=1, n_rmed
      inquire( file=trim(fn_rmed2(l)), exist=is_exist )
      if( is_exist ) then
        call rdrmed__3d( i0, i1, j0, j1, k0, k1, fn_rmed2(l), xi(k0:k1,i0:i1,j0:j1,l) )
      else
        xi(k0:k1,i0:i1,j0:j1,l) = 0.0
      end if
    end do


    !! define topography shape here
    bd(i0:i1,j0:j1,0) = depth(1)

    do k = k0, k1

      !! air column
      if( zc(k) < depth(1) ) then

        vp1 = 0.0
        vs1 = 0.0

        rho(k,i0:i1,j0:j1) = 0.001
        mu (k,i0:i1,j0:j1) = 0.0
        lam(k,i0:i1,j0:j1) = 0.0
        ! give artificially strong attenuation in air-column
        qp (k,i0:i1,j0:j1) = 10.0
        qs (k,i0:i1,j0:j1) = 10.0

        cycle
      end if

      do j=j0, j1
        do i=i0, i1
          !! chose layer
          do l=1, nlayer
            if( zc(k) >= depth(l) ) then
              rho1 = rho0(l) * ( 1 + 0.8*xi(k,i,j,tbl_rmed(l)) )
              vp1  = vp0(l)  * ( 1 +     xi(k,i,j,tbl_rmed(l)) )
              vs1  = vs0(l)  * ( 1 +     xi(k,i,j,tbl_rmed(l)) )

              call vcheck( vp1, vs1, rho1, xi(k,i,j,tbl_rmed(l) ), vmin, vmax, rhomin, vmin_under, vmax_over, rhomin_under )

              qp1  = qp0(l)
              qs1  = qs0(l)
            end if
          end do

          !! set medium parameters
          rho(k,i,j) = rho1
          mu (k,i,j) = rho1 * vs1 * vs1
          lam(k,i,j) = rho1 * ( vp1*vp1 - 2*vs1*vs1 )
          qp (k,i,j) = qp1
          qs (k,i,j) = qs1
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

    deallocate( depth, rho0, vp0, vs0, qp0, qs0 )
    deallocate( xi, tbl_rmed, fn_rmed, fn_rmed2 )

  end subroutine vmodel_lhm_rmed
  !! --------------------------------------------------------------------------------------------------------------------------- !!

end module m_vmodel_lhm_rmed
!! ----------------------------------------------------------------------------------------------------------------------------- !!
