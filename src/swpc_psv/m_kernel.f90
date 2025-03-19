#include "../shared/m_debug.h"
module m_kernel

    !! Computation kernel for FDM numerical simulation
    !!
    !! Copyright 2013-2025 Takuto Maeda. All rights reserved. This project is released under the MIT license.

    use iso_fortran_env, only: error_unit
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
        real(SP) :: bx, bz ! buoyancy
        real(MP) :: dxSxx, dxSxz, dzSxz, dzSzz
        integer :: isign
        real(MP) :: re40x, re41x, re40z, re41z

        call pwatch__on("kernel__update_vel")

        !$omp parallel &
        !$omp private(dxSxx, dzSzz, dxSxz, dzSxz) &
        !$omp private(i,k) &
        !$omp private(re40x, re41x, re40z, re41z, isign) &
        !$omp private(bx, bz)
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

                dxSxx = (Sxx(k  ,i+1) - Sxx(k  ,i  )) * re40x - (Sxx(k  ,i+2) - Sxx(k  ,i-1)) * re41x
                dzSzz = (Szz(k+1,i  ) - Szz(k  ,i  )) * re40z - (Szz(k+2,i  ) - Szz(k-1,i  )) * re41z
                dxSxz = (Sxz(k  ,i  ) - Sxz(k  ,i-1)) * re40x - (Sxz(k  ,i+1) - Sxz(k  ,i-2)) * re41x
                dzSxz = (Sxz(k  ,i  ) - Sxz(k-1,i  )) * re40z - (Sxz(k+1,i  ) - Sxz(k-2,i  )) * re41z

                !! effective buoyancy
                bx = 2.0 / (rho(k, i) + rho(k, i + 1))
                bz = 2.0 / (rho(k, i) + rho(k + 1, i))

                !! update velocity
                Vx(k,i) = Vx(k,i) + bx * (dxSxx + dzSxz) * dt
                Vz(k,i) = Vz(k,i) + bz * (dxSxz + dzSzz) * dt

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
        real(SP) :: mu2, lam2mu
        real(SP) :: nnn, pnn, npn, ppn
        real(SP) :: mu_xz
        real(SP) :: taup1, taus1, taup_plus1, taus_plus1
        real(SP) :: d2v2, dxVx_dzVz, dxVz_dzVx
        real(SP) :: Rxx_n, Rzz_n, Rxz_n
        real(SP) :: f_Rxx, f_Rzz, f_Rxz
        real(SP) :: epsl = epsilon(1.0)
        real(MP) :: dxVx, dxVz, dzVx, dzVz
        integer :: isign
        real(MP) :: re40x, re41x, re40z, re41z


        call pwatch__on("kernel__update_stress")

        !$omp parallel  &
        !$omp private( dxVx, dxVz ) &
        !$omp private( mu2, lam2mu ) &
        !$omp private( taup1, taus1, taup_plus1, taus_plus1 ) &
        !$omp private( d2v2, dxVx_dzVz ) &
        !$omp private( f_Rxx, f_Rzz ) &
        !$omp private( Rxx_n, Rzz_n ) &
        !$omp private( re40x, re41x, re40z, re41z, isign) &
        !$omp private( i, k, m )
        !$omp do &
        !$omp schedule(static,1)
        do i = ibeg_k, iend_k
            !ocl unroll('full')
            !ocl swp
            !OCL SWP_IREG_RATE(200)
            do k = kbeg_k, kend_k

                isign = sign(1, max((k - kfs_top(i)) * (kfs_bot(i) - k), &
                                    (k - kob_top(i)) * (kob_bot(i) - k)))

                re40x = rc40x + isign * rd40x
                re41x = rc41x + isign * rd41x
                re40z = rc40z + isign * rd40z
                re41z = rc41z + isign * rd41z               

                dxVx = (Vx(k  ,i  ) - Vx(k  ,i-1)) * re40x - (Vx(k  ,i+1) - Vx(k  ,i-2)) * re41x
                dzVz = (Vz(k  ,i  ) - Vz(k-1,i  )) * re40z - (Vz(k+1,i  ) - Vz(k-2,i  )) * re41z

                mu2 = 2 * mu(k, i)
                lam2mu = lam(k, i) + mu2

                taup1 = taup(k, i)
                taus1 = taus(k, i)

                d2v2 = dxVx + dzVz
                dxVx_dzVz = dxVx + dzVz

                f_Rxx = lam2mu * taup1 * d2v2 - mu2 * taus1 * dzVz
                f_Rzz = lam2mu * taup1 * d2v2 - mu2 * taus1 * dxVx

                Rxx_n = 0.0
                Rzz_n = 0.0

                !! Crank-Nicolson Method for avoiding stiff solution
                do m = 1, nm
                    Rxx(m, k, i) = c1(m) * Rxx(m, k, i) - c2(m) * f_Rxx * dt
                    Rzz(m, k, i) = c1(m) * Rzz(m, k, i) - c2(m) * f_Rzz * dt
                    Rxx_n = Rxx_n + d1(m) * Rxx(m,k,i)
                    Rzz_n = Rzz_n + d1(m) * Rzz(m,k,i)
                end do

                !! update stress components
                taup_plus1 = 1 + taup1 * (1 + d2)
                taus_plus1 = 1 + taus1 * (1 + d2)

                Sxx(k,i) = Sxx(k,i) + (lam2mu * taup_plus1 * d2v2 - mu2 * taus_plus1 * dzVz + Rxx_n) * dt
                Szz(k,i) = Szz(k,i) + (lam2mu * taup_plus1 * d2v2 - mu2 * taus_plus1 * dxVx + Rzz_n) * dt

            end do
        end do
        !$omp end do nowait
        !$omp end parallel

        !$omp parallel  &
        !$omp private( dxVz, dzVx ) &
        !$omp private( taus1, nnn, pnn, npn, ppn, mu_xz, taup_plus1, taus_plus1 ) &
        !$omp private( dxVz_dzVx ) &
        !$omp private( f_Rxz ) &
        !$omp private( Rxz_n ) &
        !$omp private( re40x, re41x, re40z, re41z, isign) &
        !$omp private( i, k, m )
        !$omp do &
        !$omp schedule(static,1)
        do i = ibeg_k, iend_k
            !ocl unroll('full')
            !ocl swp
            !OCL SWP_IREG_RATE(200)
            do k = kbeg_k, kend_k

                isign = sign(1, max((k - kfs_top(i)) * (kfs_bot(i) - k), &
                                    (k - kob_top(i)) * (kob_bot(i) - k)))

                re40x = rc40x + isign * rd40x
                re41x = rc41x + isign * rd41x
                re40z = rc40z + isign * rd40z
                re41z = rc41z + isign * rd41z               

                dxVz = (Vz(k  ,i+1) - Vz(k  ,i  )) * re40x - (Vz(k  ,i+2) - Vz(k  ,i-1)) * re41x
                dzVx = (Vx(k+1,i  ) - Vx(k  ,i  )) * re40z - (Vx(k+2,i  ) - Vx(k-1,i  )) * re41z

                taus1 = taus(k, i)

                nnn = mu(k, i)
                pnn = mu(k + 1, i)
                npn = mu(k, i + 1)
                ppn = mu(k + 1, i + 1)
                mu_xz = 4 * nnn * pnn * npn * ppn / (nnn * pnn * npn + nnn * pnn * ppn + nnn * npn * ppn + pnn * npn * ppn + epsl)

                dxVz_dzVx = dxVz + dzVx

                f_Rxz = mu_xz * taus1 * dxVz_dzVx

                Rxz_n = 0.0

                !! Crank-Nicolson Method for avoiding stiff solution
                do m = 1, nm
                    Rxz(m, k, i) = c1(m) * Rxz(m, k, i) - c2(m) * f_Rxz * dt
                    Rxz_n = Rxz_n + d1(m) * Rxz(m,k,i)
                end do

                !! update stress components
                taus_plus1 = 1 + taus1 * (1 + d2)

                Sxz(k,i) = Sxz(k,i) + (mu_xz  * taus_plus1 * dxVz_dzVx + Rxz_n) * dt

            end do
        end do
        !$omp end do nowait
        !$omp end parallel        

        !$omp barrier
        call pwatch__off("kernel__update_stress")

    end subroutine kernel__update_stress

    subroutine kernel__vmax(xmax, zmax)

        !! maximum value for terminal output
        real(SP), intent(out) :: xmax, zmax
        integer :: i

        xmax = 0.0
        zmax = 0.0
        do i = ibeg_k, iend_k
            xmax = max(xmax, abs(real(vx(kob(i) + 1, i))))
            zmax = max(zmax, abs(real(vz(kob(i) + 1, i))))
        end do

    end subroutine kernel__vmax

    subroutine memory_allocate

        !! memory allocation

        allocate (Vx(kbeg_m:kend_m, ibeg_m:iend_m), source=0.0_MP)
        allocate (Vz(kbeg_m:kend_m, ibeg_m:iend_m), source=0.0_MP)
        allocate (Sxx(kbeg_m:kend_m, ibeg_m:iend_m), source=0.0_MP)
        allocate (Szz(kbeg_m:kend_m, ibeg_m:iend_m), source=0.0_MP)
        allocate (Sxz(kbeg_m:kend_m, ibeg_m:iend_m), source=0.0_MP)

        if (nm > 0) then
            allocate (Rxx(1:nm, kbeg_m:kend_m, ibeg_m:iend_m), source=0.0)
            allocate (Rzz(1:nm, kbeg_m:kend_m, ibeg_m:iend_m), source=0.0)
            allocate (Rxz(1:nm, kbeg_m:kend_m, ibeg_m:iend_m), source=0.0)
            allocate (c1(1:nm), c2(1:nm), d1(1:nm))
        end if

    end subroutine memory_allocate

end module m_kernel
