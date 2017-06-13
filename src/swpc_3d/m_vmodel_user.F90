!! ----------------------------------------------------------------------------------------------------------------------------- !!
!>
!! User-routine for defining velocity/attenuation structure
!!
!! @copyright
!!   Copyright 2013-2017 Takuto Maeda. All rights reserved. This project is released under the MIT license.
!<
!! ----
module m_vmodel_user

  use m_std
  use m_global
  use m_geomap
  use m_fdtool
  use m_readini
  implicit none
  private
  save

  public :: vmodel_user

contains

  !! --------------------------------------------------------------------------------------------------------------------------- !!
  !>
  !! Define meidum velocity, density and attenuation
  !!
  !! This is a user-specific routine to define original veloicty model.
  !!
  !! Input:
  !!    io_prm                          :: I/O number of parameter file (file has been opened already)
  !!    i0,i1, j0,j1, k0,k1             :: model area by indices in i-,j-,and k- directions
  !!    xc(i0:i1), yc(j0:j1), zc(k0:k1) :: Cartesian coordinate location
  !!    vcut                            :: cut-off velocity specified by input parameter
  !!
  !! Output:
  !!    rho(k0:k1, i0:i1, j0:j1)        :: mass density (usually in g/cm^3)
  !!    lam(k0:k1, i0:i1, j0:j1)        :: Lame's parameter (usually in (g/cm^3) * (km/s)^2)
  !!    mu (k0:k1, i0:i1, j0:j1)        :: Lame's parameter (usually in (g/cm^3) * (km/s)^2)
  !!    qp (k0:k1, i0:i1, j0:j1)        :: Attenuation QP
  !!    qs (k0:k1, i0:i1, j0:j1)        :: Attenuation QS
  !!    bd (i0:i1, j0:j1, 0:NBD)        :: Boundary depths
  !!
  !! Note:
  !! bd(:,:,0) are treated as topography shape for output.
  !!    this is only for output and visualization. topography in the simulation will be automatically detected by medium params.
  !! bd(:,:,1:NBD) may contain internal boundary depths. The boundary number can be specified as source depth or station depth. 
  !!
  !<
  !! ----
  subroutine vmodel_user( io_prm, i0, i1, j0, j1, k0, k1, xc, yc, zc, vcut, rho, lam, mu, Qp, Qs, bd )

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
    real(SP), intent(out) :: lam ( k0:k1, i0:i1, j0:j1 )    !< Lame's parameter lambda [ (g/cm^3) * (km/s)**2 ]
    real(SP), intent(out) :: mu  ( k0:k1, i0:i1, j0:j1 )    !< Lame's parameter mu     [ (g/cm^3) * (km/s)**2 ]
    real(SP), intent(out) :: qp  ( k0:k1, i0:i1, j0:j1 )    !< P-wave attenuation
    real(SP), intent(out) :: qs  ( k0:k1, i0:i1, j0:j1 )    !< S-wave attenuation
    real(SP), intent(out) :: bd  ( i0:i1, j0:j1, 0:NBD )    !< Boundary depths
    !! --

    integer  :: i, j, k
    real(SP) :: vp0, vs0, rho0, qp0, qs0, topo0
    real(SP) :: vp1, vs1
    real(SP) :: dum
    !! ----

    !!
    !! The following dummy code is an example how to discribe the routine. 
    !!

    !!
    !! subroutine readini() can access parameters defined in the input file.
    !! Any original parameters can be added in the input file. 
    !!
    call readini( io_prm, 'vp0',    vp0, 5.0 )
    call readini( io_prm, 'vs0',    vs0, vp0/sqrt(3.0) )
    call readini( io_prm, 'rho0',   rho0, 2.7 )
    call readini( io_prm, 'qp0',    qp0, 1000000.0 )
    call readini( io_prm, 'qs0',    qs0, 1000000.0 )
    call readini( io_prm, 'topo0', topo0, 0.0 )

    !!
    !! The medium parameter must be set from given region (i0:i1, j0:j1, k0:k1)
    !! Note that the order of indices is k->i->j, for improving performance
    !! 
    do j = j0, j1
      do i = i0, i1

        !! define topography shape here
        bd(i,j,0) = topo0

        do k = k0, k1

          if( zc( k ) > bd(i,j,0) ) then

            !! elastic medium
            rho(k,i,j) = rho0
            mu (k,i,j) = rho(k,i,j) * vs0 * vs0
            lam(k,i,j) = rho(k,i,j) * ( vp0*vp0 - 2*vs0*vs0 )
            qp (k,i,j) = qp0
            qs (k,i,j) = qs0

          else if ( zc (k) > 0.0 ) then

            !!
            !! ocean column
            !!
            !! The code treat the uppermost layer as ocean column if P-wave velocity is finite and S-wave velocity is zero
            !!
            vp1 = 1.5
            vs1 = 0.0

            rho(k,i,j) = 1.0
            mu (k,i,j) = rho(k,i,j) * vs1 * vs1
            lam(k,i,j) = rho(k,i,j) * ( vp1*vp1 - 2*vs1*vs1 )
            qp (k,i,j) = 1000000.0 ! effectively no attenuation in ocean column
            qs (k,i,j) = 1000000.0

          else

            !!
            !! air column
            !!
            !! The air column must have zero P- & S-wave velocity (i.e., mu=lam=0)
            !! Please use non-zero but very small density (e.g., 0.001) for avoiding zero division with satisfying boundary cond.
            !! Since waves do not penetrate to the air column, qp and qs does not affect. Just set dummy.
            !!
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

    !! dummy value
    bd(:,:,1:NBD) = -9999

    ! substitute to a dummy variable for avoiding compiler warnings
    dum = xc(i0)
    dum = yc(j0)
    dum = zc(k0)
    dum = vcut

  end subroutine vmodel_user
  !! --------------------------------------------------------------------------------------------------------------------------- !!


end module m_vmodel_user
!! ----------------------------------------------------------------------------------------------------------------------------- !!
