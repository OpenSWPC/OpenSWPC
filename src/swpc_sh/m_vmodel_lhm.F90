!! ----------------------------------------------------------------------------------------------------------------------------- !!
!>
!! 1D velocity structure
!!
!! @copyright
!!   Copyright 2013-2021 Takuto Maeda. All rights reserved. This project is released under the MIT license.
!<
!! ----
#include "m_debug.h"
module m_vmodel_lhm

  use m_std
  use m_debug
  use m_readini
  use m_global
  implicit none
  private
  save

  public :: vmodel_lhm

contains

  !! --------------------------------------------------------------------------------------------------------------------------- !!
  !>
  !! Define meidum velocity, density and attenuation
  !<
  !! ----
  subroutine vmodel_lhm( io_prm, i0, i1, k0, k1, xc, zc, vcut, rho, lam, mu, Qp, Qs, bd )

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
    character(256) :: fn_lhm
    integer  :: k, l
    real(SP), allocatable, dimension(:) :: vp0, vs0, rho0, qp0, qs0, depth
    real(SP) :: vp1, vs1, rho1, qp1, qs1
    real(SP) :: dum
    integer :: ierr
    integer :: io_vel
    logical :: is_exist
    integer :: nlayer
    character(256) :: adum
    logical :: earth_flattening
    real(SP) :: zs(k0:k1) ! spherical depth for earth_flattening
    real(SP) :: Cv(k0:k1) ! velocity scaling coefficient for earth_flattening
    !! ----

    call readini( io_prm, 'fn_lhm', fn_lhm, '' )

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

    inquire( file=fn_lhm, exist=is_exist )
    if( .not. is_exist ) then
      write(STDERR,*) "ERROR [m_vmodel_lhm]: velocity file "//trim(fn_lhm)//" does not exist"
      call mpi_finalize(ierr)
      stop
    end if

    call std__getio( io_vel )
    open( io_vel, file=fn_lhm, status='old', action='read', iostat=ierr )
    call std__countline( io_vel, nlayer, '#' )
    allocate( depth(nlayer), rho0(nlayer), vp0(nlayer), vs0(nlayer), qp0(nlayer), qs0(nlayer) )

    l = 0
    do
      read(io_vel,'(A256)', iostat=ierr ) adum
      if( ierr /= 0 ) exit
      adum = trim(adjustl(adum))
      if( trim(adum) == '' ) cycle
      if( adum(1:1) == "#" ) cycle
      l = l + 1
      read(adum,*) depth(l), rho0(l), vp0(l), vs0(l), qp0(l), qs0(l)
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

    !! define topography shape here
    bd(i0:i1,0) = depth(1)

    do k = k0, k1

      !! air/ocean column
      if( zs(k) < depth(1) ) then

        if( zs(k) < 0.0 ) then 

          vp1 = 0.0
          vs1 = 0.0
          
          rho(k,i0:i1) = 0.001
          mu (k,i0:i1) = 0.0
          lam(k,i0:i1) = 0.0
          ! give artificially strong attenuation in air-column
          qp (k,i0:i1) = 10.0
          qs (k,i0:i1) = 10.0
          
        else

          vp1 = Cv(k) * 1.5
          vs1 = 0.0
          
          rho(k,i0:i1) = 1.0 * Cv(k)**(-5)
          mu (k,i0:i1) = 0.0
          lam(k,i0:i1) = 1.0 * vp1 * vp1

          qp (k,i0:i1) = 1000000.0
          qs (k,i0:i1) = 1000000.0
          
        end if
        
        cycle
      end if

      !! chose layer
      do l=1, nlayer
        if( zs(k) >= depth(l) ) then
          rho1 = rho0(l) * Cv(k)**(-5)
          vp1  = Cv(k) * vp0(l)
          vs1  = Cv(k) * vs0(l)
          qp1  = qp0(l)
          qs1  = qs0(l)
        end if
      end do

      !! set medium parameters
      rho(k,i0:i1) = rho1
      mu (k,i0:i1) = rho1 * vs1 * vs1
      lam(k,i0:i1) = rho1 * ( vp1*vp1 - 2*vs1*vs1 )
      qp (k,i0:i1) = qp1
      qs (k,i0:i1) = qs1

    end do

    ! dummy value
    bd(:,1:NBD) = -9999

    ! substitute to a dummy variable for avoiding compiler warnings
    dum = xc(i0)

    deallocate( depth, rho0, vp0, vs0, qp0, qs0 )

  end subroutine vmodel_lhm
  !! --------------------------------------------------------------------------------------------------------------------------- !!


end module m_vmodel_lhm
!! ----------------------------------------------------------------------------------------------------------------------------- !!
