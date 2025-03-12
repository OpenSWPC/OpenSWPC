#include "../shared/m_debug.h"
module m_kernel

    !! Computation kernel for FDM numerical simulation
    !!
    !! Copyright 2013-2025 Takuto Maeda. All rights reserved. This project is released under the MIT license.

    use iso_fortran_env
    use m_std
    use m_debug
    use m_global
    use m_pwatch
    use m_medium

    implicit none
    private
    save

    public :: kernel__setup
    public :: kernel__update_vel
    public :: kernel__update_stress
    public :: kernel__vmax

    real(SP), allocatable :: c1(:), c2(:), d1(:)
    real(SP) :: d2
    real(MP) :: rc40x, rc41x, rc40z, rc41z
    real(MP) :: rd40x, rd41x, rd40z, rd41z

contains

    subroutine kernel__setup()

        !! Setup
        !!
        !! this routine MUST BE CALLED AFTER medium__setup since it uses viscoelastic function
        !! use ts(:) array, dx, dy, dz from global module

        integer :: m

        call pwatch__on("kernel__setup")

        if (.not. medium__initialized()) then
            write (error_unit, '(A)') 'ERROR [kernel__setup]: call medium__setup() before kernel__setup()'
            stop
        end if

        rc40x = 17.0_MP / 16.0_MP / dx
        rc40z = 17.0_MP / 16.0_MP / dz
        rc41x =  1.0_MP / 48.0_MP / dx
        rc41z =  1.0_MP / 48.0_MP / dz
        rd40x = -1.0_MP / 16.0_MP / dx
        rd40z = -1.0_MP / 16.0_MP / dz
        rd41x = -1.0_MP / 48.0_MP / dx
        rd41z = -1.0_MP / 48.0_MP / dz       

        call memory_allocate()

        if (nm > 0) then

            do m = 1, nm
                c1(m) = (2 * ts(m) - dt) / (2 * ts(m) + dt)
                c2(m) = (2) / (2 * ts(m) + dt) / nm
            end do
            d2 = sum(dt / (2 * ts(:) - dt)) / nm
            do m = 1, nm
                d1(m) = 2 * ts(m) / (2 * ts(m) - dt)
            end do

        end if

        call pwatch__off("kernel__setup")

    end subroutine kernel__setup

    subroutine kernel__update_vel()

        !! Update vel for one time step

        integer :: i, k
        real(SP) :: by ! buoyancy
        real(MP) :: dzSyz, dxSxy
        integer :: isign
        real(MP) :: re40x, re41x, re40z, re41z

        call pwatch__on("kernel__update_vel")

        !$omp parallel &
        !$omp private(dzSyz, dxSxy, i, k, by, isign, re40x, re41x, re40z, re41z)
        !$omp do &
        !$omp schedule(static,1)
        do i = ibeg_k, iend_k

            do k = kbeg_k, kend_k

                isign = sign(1, max((k - kfs_top(i)) * (kfs_bot(i) - k), &
                                    (k - kob_top(i)) * (kob_bot(i) - k)))

                re40x = rc40x + isign * rd40x
                re41x = rc41x + isign * rd41x
                re40z = rc40z + isign * rd40z
                re41z = rc41z + isign * rd41z               

                dzSyz = (Syz(k, i) - Syz(k-1,i  )) * re40z - (Syz(k+1,i  ) - Syz(k-2,i  )) * re41z
                dxSxy = (Sxy(k, i) - Sxy(k  ,i-1)) * re40x - (Sxy(k  ,i+1) - Sxy(k  ,i-2)) * re41x

                by = 1.0 / rho(k, i)

                Vy(k,i) = Vy(k,i) + by * (dxSxy + dzSyz) * dt

            end do
        end do
        !$omp end do nowait
        !$omp end parallel

        !$omp barrier

        call pwatch__off("kernel__update_vel")

    end subroutine kernel__update_vel

    subroutine kernel__update_stress()

        !! Update stress for one time step

        integer :: i, k, m
        real(SP) :: nnn, pnn, npn
        real(SP) :: mu_xy, mu_yz
        real(SP) :: taus1, taus_plus1
        real(SP) :: Ryz_n, Rxy_n
        real(SP) :: f_Rxy, f_Ryz
        real(SP) :: epsl = epsilon(1.0)
        real(MP) :: dxVy, dzVy
        integer :: isign
        real(MP) :: re40x, re41x, re40z, re41z

        call pwatch__on("kernel__update_stress")

        !$omp parallel  &
        !$omp private( dxVy, dzVy, nnn, pnn, npn, mu_yz, mu_xy ) &
        !$omp private( taus1, taus_plus1, f_Ryz, f_Rxy,  Ryz_n, Rxy_n ) &
        !$omp private( re40x, re41x, re40z, re41z, isign) &
        !$omp private( i, k, m)
        !$omp do &
        !$omp schedule(static,1)
        do i = ibeg_k, iend_k

            !ocl unroll('full')
            do k = kbeg_k, kend_k

                isign = sign(1, max((k - kfs_top(i)) * (kfs_bot(i) - k), &
                                    (k - kob_top(i)) * (kob_bot(i) - k)))

                re40x = rc40x + isign * rd40x
                re41x = rc41x + isign * rd41x
                re40z = rc40z + isign * rd40z
                re41z = rc41z + isign * rd41z      

                dxVy = (Vy(k  ,i+1) - Vy(k  ,i  )) * re40x - (Vy(k  ,i+2) - Vy(k  ,i-1)) * re41x
                dzVy = (Vy(k+1,i  ) - Vy(k  ,i  )) * re40z - (Vy(k+2,i  ) - Vy(k-1,i  )) * re41z

                !! effective rigidity for shear stress components
                nnn = mu(k, i)
                pnn = mu(k + 1, i)
                npn = mu(k, i + 1)
                mu_xy = 2 * nnn * npn / (nnn + npn + epsl)
                mu_yz = 2 * nnn * pnn / (nnn + pnn + epsl)
                taus1 = taus(k, i)

                !! update memory variables
                !! working variables for combinations of velocity derivatives
                f_Ryz = mu_yz * taus1 * dzVy
                f_Rxy = mu_xy * taus1 * dxVy

                Ryz_n = 0.0
                Rxy_n = 0.0

                !! Crank-Nicolson Method for avoiding stiff solution
                do m = 1, nm
                    Ryz(m, k, i) = c1(m) * Ryz(m, k, i) - c2(m) * f_Ryz * dt
                    Rxy(m, k, i) = c1(m) * Rxy(m, k, i) - c2(m) * f_Rxy * dt
                    Ryz_n = Ryz_n + d1(m) * Ryz(m,k,i)
                    Rxy_n = Rxy_n + d1(m) * Rxy(m,k,i)
                end do


                !! update stress components
                taus_plus1 = 1 + taus1 * (1 + d2)

                Syz(k,i) = Syz(k,i) + (mu_yz * taus_plus1 * dzVy + Ryz_n ) * dt
                Sxy(k,i) = Sxy(k,i) + (mu_xy * taus_plus1 * dxVy + Rxy_n ) * dt

            end do
        end do
        !$omp end do nowait
        !$omp end parallel

        !$omp barrier

        call pwatch__off("kernel__update_stress")

    end subroutine kernel__update_stress

    subroutine kernel__vmax(ymax)

        !! maximum value for terminal output

        real(SP), intent(out) :: ymax
        integer :: i

        ymax = 0.0
        do i = ibeg_k, iend_k
            ymax = max(ymax, abs(real(vy(kob(i) + 1, i))))
        end do

    end subroutine kernel__vmax

    subroutine memory_allocate

        !! memory allocation

        allocate (Vy(kbeg_m:kend_m, ibeg_m:iend_m), source=0.0_MP)
        allocate (Syz(kbeg_m:kend_m, ibeg_m:iend_m), source=0.0_MP)
        allocate (Sxy(kbeg_m:kend_m, ibeg_m:iend_m), source=0.0_MP)

        if (nm > 0) then
            allocate (Ryz(1:nm, kbeg_k:kend_k, ibeg_k:iend_k), source=0.0)
            allocate (Rxy(1:nm, kbeg_k:kend_k, ibeg_k:iend_k), source=0.0)
            allocate (c1(1:nm), c2(1:nm), d1(1:nm))
        end if

    end subroutine memory_allocate

end module m_kernel
