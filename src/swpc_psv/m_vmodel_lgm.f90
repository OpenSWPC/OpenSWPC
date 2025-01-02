#include "../shared/m_debug.h"
module m_vmodel_lgm

    !! 1D layered velocity structure with gradient
    !!
    !! Copyright 2013-2025 Takuto Maeda. All rights reserved. This project is released under the MIT license.

    use iso_fortran_env, only: error_unit
    use m_std
    use m_debug
    use m_readini
    use m_global
    use m_seawater
    implicit none
    private
    save

    public :: vmodel_lgm

contains

    subroutine vmodel_lgm(io_prm, i0, i1, k0, k1, xc, zc, vcut, rho, lam, mu, Qp, Qs, bd)

        integer,  intent(in)  :: io_prm
        integer,  intent(in)  :: i0, i1              !< i-region
        integer,  intent(in)  :: k0, k1              !< k-region
        real(SP), intent(in)  :: xc(i0:i1)           !< x-coordinate location
        real(SP), intent(in)  :: zc(k0:k1)           !< z-coordinate location
        real(SP), intent(in)  :: vcut                !< cut-off minimum velocity
        real(SP), intent(out) :: rho(k0:k1, i0:i1)   !< mass density [g/cm^3]
        real(SP), intent(out) :: lam(k0:k1, i0:i1)   !< Lame's parameter lambda [ (g/cm^3) * (km/s) ]
        real(SP), intent(out) :: mu(k0:k1, i0:i1)    !< Lame's parameter mu     [ (g/cm^3) * (km/s) ]
        real(SP), intent(out) :: qp(k0:k1, i0:i1)    !< P-wave attenuation
        real(SP), intent(out) :: qs(k0:k1, i0:i1)    !< S-wave attenuation
        real(SP), intent(out) :: bd(i0:i1, 0:NBD)    !< Boundary depths
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
        logical :: use_munk
        logical :: earth_flattening
        real(SP) :: zs(k0:k1) ! spherical depth for earth_flattening
        real(SP) :: Cv(k0:k1) ! velocity scaling coefficient for earth_flattening

        call readini(io_prm, 'fn_lhm', fn_lhm, '')
        inquire (file=fn_lhm, exist=is_exist)
        call assert(is_exist)

        !! seawater
        call readini(io_prm, 'munk_profile', use_munk, .false.)
        call seawater__init(use_munk)

        call readini(io_prm, 'earth_flattening', earth_flattening, .false.)
        if (earth_flattening) then
            do k = k0, k1
                zs(k) = R_EARTH - R_EARTH * exp(-zc(k) / R_EARTH)
                Cv(k) = exp(zc(k) / R_EARTH)
            end do
        else
            zs(:) = zc(:)
            Cv(:) = 1.0
        end if

        open (newunit=io_vel, file=fn_lhm, status='old', action='read', iostat=ierr)
        call std__countline(io_vel, nlayer, '#')
        allocate (depth(nlayer), rho0(nlayer), vp0(nlayer), vs0(nlayer), qp0(nlayer), qs0(nlayer))

        l = 0
        do
            read (io_vel, '(A256)', iostat=ierr) adum
            if (ierr /= 0) exit
            adum = trim(adjustl(adum))
            if (trim(adum) == '') cycle
            if (adum(1:1) == "#") cycle
            l = l + 1
            read (adum, *) depth(l), rho0(l), vp0(l), vs0(l), qp0(l), qs0(l)
        end do
        close (io_vel)

        !! velocity cut-off
        do l = nlayer - 1, 1, -1
            if ((vp0(l) < vcut .or. vs0(l) < vcut) .and. (vp0(l) > 0 .and. vs0(l) > 0)) then
                vp0 (l) = vp0 (l+1)
                vs0 (l) = vs0 (l+1)
                rho0(l) = rho0(l+1)
                qp0 (l) = qp0 (l+1)
                qs0 (l) = qs0 (l+1)
            end if
        end do

        !! defne topography shape here
        bd(i0:i1, 0) = depth(1)

        do k = k0, k1

            !! air/ocean column
            if (zs(k) < depth(1)) then

                if (zs(k) < 0.0) then

                    rho1 = 0.001
                    vp1  = 0.0
                    vs1  = 0.0
                    qp1  = 10.0   ! give artificially strong attenuation in air-column
                    qs1  = 10.0

                else

                    rho1 = 1.0
                    vp1  = Cv(k) * seawater__vel(zs(k))
                    vs1  = 0.0
                    qp1  = 1000000.0
                    qs1  = 1000000.0

                end if

            else

                ! initialize by the values of the lowermost layer
                rho1 = rho0(nlayer)
                vp1  = Cv(k) * vp0(nlayer)
                vs1  = Cv(k) * vs0(nlayer)
                qp1  = qp0(nlayer)
                qs1  = qs0(nlayer)

                !! chose layer
                do l=1, nlayer-1
                    if ( depth(l) <= zs(k) .and. zs(k) < depth(l+1) ) then
                        rho1 =          rho0(l) + (rho0(l+1) - rho0(l)) * (zs(k) - depth(l)) / (depth(l+1) - depth(l))
                        vp1  = Cv(k) * ( vp0(l) + ( vp0(l+1) -  vp0(l)) * (zs(k) - depth(l)) / (depth(l+1) - depth(l)))
                        vs1  = Cv(k) * ( vs0(l) + ( vs0(l+1) -  vs0(l)) * (zs(k) - depth(l)) / (depth(l+1) - depth(l)))
                        qp1  =           qp0(l) + ( qp0(l+1) -  qp0(l)) * (zs(k) - depth(l)) / (depth(l+1) - depth(l))
                        qs1  =           qs0(l) + ( qs0(l+1) -  qs0(l)) * (zs(k) - depth(l)) / (depth(l+1) - depth(l))
                        exit
                    end if
                end do
            end if

            ! set medium parameters
            rho(k, i0:i1) = rho1
            mu (k, i0:i1) = rho1 * vs1 * vs1
            lam(k, i0:i1) = rho1 * (vp1 * vp1 - 2 * vs1 * vs1)
            qp (k, i0:i1) = qp1
            qs (k, i0:i1) = qs1

        end do

        ! dummy value
        bd(:, 1:NBD) = -9999

        ! substitute to a dummy variable for avoiding compiler warnings
        dum = xc(i0)

        deallocate (depth, rho0, vp0, vs0, qp0, qs0)

    end subroutine vmodel_lgm

end module m_vmodel_lgm
