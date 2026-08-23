#include "../shared/m_debug.h"
module m_absorb_p

    !! Absorbing Boundary Condition: ADE-CFS PML based on Zhang and Shen
    !!
    !! Copyright 2013-2026 Takuto Maeda. All rights reserved. This project is released under the MIT license.

    use m_std
    use m_debug
    use m_global
    use m_fdtool
    use m_readini

    implicit none
    private
    save

    public :: absorb_p__setup
    public :: absorb_p__update_stress
    public :: absorb_p__update_vel

    real(SP), allocatable :: gxc(:, :), gxe(:, :) !< damping profile along x at center/edge of voxel
    real(SP), allocatable :: gzc(:, :), gze(:, :) !< damping profile along z at center/edge of voxel

    real(SP), allocatable :: axVx(:), azVx(:)
    real(SP), allocatable :: axVz(:), azVz(:)
    real(SP), allocatable :: axSxx(:), azSxz(:)
    real(SP), allocatable :: axSxz(:), azSzz(:)

    real(SP) :: r20x, r20z
    real(MP) :: rc40x, rc41x, rc40z, rc41z
    real(MP) :: rd40x, rd41x, rd40z, rd41z

contains

    subroutine absorb_p__setup(io_prm)

        !! set PML sponge

        integer, intent(in) :: io_prm
        integer  :: i, k
        real(SP) :: hx, hz
        integer  :: idum

        !! derivative coefficient
        r20x = real(1.0 / dx)
        r20z = real(1.0 / dz)
        rc40x = 17.0_MP / 16.0_MP / dx
        rc40z = 17.0_MP / 16.0_MP / dz
        rc41x =  1.0_MP / 48.0_MP / dx
        rc41z =  1.0_MP / 48.0_MP / dz
        rd40x = -1.0_MP / 16.0_MP / dx
        rd40z = -1.0_MP / 16.0_MP / dz
        rd41x = -1.0_MP / 48.0_MP / dx
        rd41z = -1.0_MP / 48.0_MP / dz       

        !! damping profile
        allocate (gxc(4, ibeg:iend), gxe(4, ibeg:iend))
        allocate (gzc(4, kbeg:kend), gze(4, kbeg:kend))

        hx = real(na * dx)
        hz = real(na * dz)

        do i = ibeg, iend
            call damping_profile(xc(i), hx, xbeg, xend, gxc(:, i))
            call damping_profile(xc(i) + real(dx) / 2.0, hx, xbeg, xend, gxe(:, i))
        end do

        do k = kbeg, kend
            call damping_profile(zc(k), hz, zbeg, zend, gzc(:, k))
            call damping_profile(zc(k) + real(dz) / 2.0, hz, zbeg, zend, gze(:, k))
        end do

        !! memory allocation
        allocate (axVx (n_sponge_cell), source=0.0)
        allocate (azVx (n_sponge_cell), source=0.0)
        allocate (axVz (n_sponge_cell), source=0.0)
        allocate (azVz (n_sponge_cell), source=0.0)
        allocate (axSxx(n_sponge_cell), source=0.0)
        allocate (azSxz(n_sponge_cell), source=0.0)
        allocate (axSxz(n_sponge_cell), source=0.0)
        allocate (azSzz(n_sponge_cell), source=0.0)

        idum = io_prm

        !$acc enter data copyin(axVx, azVx, axVz, azVz, axSxx, azSxz, axSxz, azSzz, &
        !$acc                   gxc, gxe, gzc, gze)

    end subroutine absorb_p__setup

    subroutine absorb_p__update_vel

        integer :: ibox

        if(pw_mode) call set_stress_boundary()

        !! time-marching
        do ibox=1, 4
            if(box(ibox)%ncell > 0) call update_vel_core(box(ibox))
        end do

    end subroutine absorb_p__update_vel

    subroutine set_stress_boundary()

        integer :: k

        if (idx == 0) then
#ifdef _OPENACC
            !$acc kernels present(Sxx, Szz, Sxz)
            !$acc loop independent
#else
            !$omp parallel private(k)
            !$omp do schedule(dynamic)
#endif
            do k = kbeg, kend
                Sxx(k,0) = 2 * Sxx(k,1) - Sxx(k,2)
                Szz(k,0) = 2 * Szz(k,1) - Szz(k,2)
                Sxz(k,0) = 2 * Sxz(k,1) - Sxz(k,2)
            end do
#ifdef _OPENACC
            !$acc end kernels
#else
            !$omp end do nowait
            !$omp end parallel
#endif

        end if

        if (idx == nproc_x - 1) then
#ifdef _OPENACC
            !$acc kernels present(Sxx, Szz, Sxz)
            !$acc loop independent
#else
            !$omp parallel private(k)
            !$omp do schedule(dynamic)
#endif
            do k = kbeg, kend
                Sxx(k,nx+1) = 2 * Sxx(k,nx) - Sxx(k,nx-1)
                Szz(k,nx+1) = 2 * Szz(k,nx) - Szz(k,nx-1)
                Sxz(k,nx+1) = 2 * Sxz(k,nx) - Sxz(k,nx-1)
            end do
#ifdef _OPENACC
            !$acc end kernels
#else
            !$omp end do nowait
            !$omp end parallel
#endif
        end if

    end subroutine set_stress_boundary


    subroutine update_vel_core(bb)

        type(t_box), intent(in) :: bb
        integer :: p
        integer :: i, k
        real(SP) :: bx, bz
        real(MP) :: dxSxx, dzSxz, dxSxz, dzSzz
        integer :: isign
        real(MP) :: re40x, re41x, re40z, re41z

#ifdef _OPENACC
        !$acc kernels &
        !$acc present(Sxx, Szz, Sxz, Vx, Vz, rho, axSxx, azSxz, axSxz, azSzz, &
        !$acc         gxc, gxe, gzc, gze, bb,  kfs_top, kfs_bot, kob_top, kob_bot)
        !$acc loop independent collapse(2)
#else
        !$omp parallel &
        !$omp private( dxSxx, dzSzz, dxSxz, dzSxz ) &
        !$omp private(re40x, re41x, re40z, re41z, isign) &
        !$omp private( i, k, p )
        !$omp do &
        !$omp schedule(dynamic)
#endif
        do i = bb%ib, bb%ie

            do k = bb%kb, bb%ke

                p = bb%offset + (i-bb%ib) * bb%nz + (k - bb%kb + 1)

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

                bx = 2.0 / (rho(k,i) + rho(k,i + 1))
                bz = 2.0 / (rho(k,i) + rho(k + 1, i))

                !! Velocity Updates
                Vx(k,i) = Vx(k,i) &
                        + bx * (gxe(1,i) * dxSxx    + gzc(1,k) * dzSxz &
                              + gxe(2,i) * axSxx(p) + gzc(2,k) * azSxz(p)) * dt

                Vz(k,i) = Vz(k,i) &
                        + bz * (gxc(1,i) * dxSxz    + gze(1,k) * dzSzz &
                              + gxc(2,i) * axSxz(p) + gze(2,k) * azSzz(p)) * dt

                !! ADE updates
                axSxx(p) = real(gxe(3,i) * axSxx(p) + gxe(4,i) * dxSxx * dt)
                azSxz(p) = real(gzc(3,k) * azSxz(p) + gzc(4,k) * dzSxz * dt)
                axSxz(p) = real(gxc(3,i) * axSxz(p) + gxc(4,i) * dxSxz * dt)
                azSzz(p) = real(gze(3,k) * azSzz(p) + gze(4,k) * dzSzz * dt)

            end do
        end do
#ifdef _OPENACC
        !$acc end kernels
#else
        !$omp end do nowait
        !$omp end parallel
        !$omp barrier
#endif

    end subroutine update_vel_core

    subroutine absorb_p__update_stress

        integer :: ibox

        if (pw_mode) call set_vel_boundary()

        do ibox=1, 4
            if(box(ibox)%ncell > 0) call update_stress_core(box(ibox))
        end do

    end subroutine absorb_p__update_stress


    subroutine set_vel_boundary()

        integer :: i, k

        if (idx == 0) then
#ifdef _OPENACC
            !$acc kernels present(Vx, Vz)
            !$acc loop independent
#else
            !$omp parallel private(k)
            !$omp do schedule(dynamic)
#endif
            do k = kbeg, kend
                Vx(k,0) = 2 * Vx(k,1) - Vx(k,2)
                Vz(k,0) = 2 * Vz(k,1) - Vz(k,2)
            end do
#ifdef _OPENACC
            !$acc end kernels
#else
            !$omp end do nowait
            !$omp end parallel
#endif
        end if

        if (idx == nproc_x - 1) then
#ifdef _OPENACC
            !$acc kernels present(Vx, Vz)
            !$acc loop independent
#else
            !$omp parallel private(k)
            !$omp do schedule(dynamic)
#endif
            do k = kbeg, kend
                Vx(k,nx + 1) = 2 * Vx(k,nx) - Vx(k,nx - 1)
                Vz(k,nx + 1) = 2 * Vz(k,nx) - Vz(k,nx - 1)
            end do
#ifdef _OPENACC
            !$acc end kernels
#else
            !$omp end do nowait
            !$omp end parallel
#endif
        end if

    end subroutine set_vel_boundary


    subroutine update_stress_core(bb)

        type(t_box), intent(in) :: bb

        integer :: i, k, p, m
        real(SP) :: mu2, lam2mu
        real(SP) :: dxVx_ade, dzVz_ade, dxVz_ade, dzVx_ade
        real(SP) :: nnn, pnn, npn, ppn
        real(SP) :: mu_xz
        real(SP) :: epsl = epsilon(1.0)
        real(MP) :: dxVx, dzVx, dxVz, dzVz
        real(SP) :: taup1, taus1, taup_plus1, taus_plus1
        real(SP) :: d2v2, dxVx_dzVz, dxVz_dzVx
        real(SP) :: Rxx_n, Rzz_n, Rxz_n
        real(SP) :: f_Rxx, f_Rzz, f_Rxz
        integer :: isign
        real(MP) :: re40x, re41x, re40z, re41z

        !! Time-marching
#ifdef _OPENACC
        !$acc kernels &
        !$acc present(Vx, Vz, Sxx, Szz, axVx, azVz, gxc, gxe, gzc, gze, mu, lam, bb, &
        !$acc         c1, c2, d1, d2, taup, taus, kfs_top, kfs_bot, kob_top, kob_bot, Rxx, Rzz)
        !$acc loop independent collapse(2)
#else
        !$omp parallel &
        !$omp private( dxVx, dxVz) &
        !$omp private( mu2, lam2mu ) &
        !$omp private( dxVx_ade, dzVz_ade) &
        !$omp private( i, k, m, p ) &
        !$omp private( taup1, taus1, taup_plus1, taus_plus1 ) &
        !$omp private( d2v2, dxVx_dzVz) &
        !$omp private( f_Rxx, f_Rzz) &
        !$omp private( Rxx_n, Rzz_n ) &
        !$omp private( re40x, re41x, re40z, re41z, isign)
        !$omp do &
        !$omp schedule(dynamic)
#endif
        do i = bb%ib, bb%ie

            !ocl unroll('full')
            !ocl swp
            !OCL SWP_IREG_RATE(200)
            do k = bb%kb, bb%ke

                p = bb%offset + (i-bb%ib) * bb%nz + (k - bb%kb + 1)

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

                dxVx_ade = real(gxc(1,i) * dxVx + gxc(2,i) * axVx(p))
                dzVz_ade = real(gzc(1,k) * dzVz + gzc(2,k) * azVz(p))

                d2v2 = dxVx_ade + dzVz_ade
                dxVx_dzVz = dxVx_ade + dzVz_ade

                f_Rxx = lam2mu * taup1 * d2v2 - mu2 * taus1 * dzVz_ade
                f_Rzz = lam2mu * taup1 * d2v2 - mu2 * taus1 * dxVx_ade

                Rxx_n = 0.0
                Rzz_n = 0.0

                !$acc loop seq reduction(+:Rxx_n,Rzz_n)
                do m = 1, nm
                    Rxx(m, k, i) = c1(m) * Rxx(m, k, i) - c2(m) * f_Rxx * dt
                    Rzz(m, k, i) = c1(m) * Rzz(m, k, i) - c2(m) * f_Rzz * dt
                    Rxx_n = Rxx_n + d1(m) * Rxx(m,k,i)
                    Rzz_n = Rzz_n + d1(m) * Rzz(m,k,i)
                end do                

                taup_plus1 = 1 + taup1 * (1 + d2)
                taus_plus1 = 1 + taus1 * (1 + d2)

                Sxx(k,i) = Sxx(k,i) + (lam2mu * taup_plus1 * d2v2 - mu2 * taus_plus1 * dzVz_ade + Rxx_n) * dt
                Szz(k,i) = Szz(k,i) + (lam2mu * taup_plus1 * d2v2 - mu2 * taus_plus1 * dxVx_ade + Rzz_n) * dt

                axVx(p) = real(gxc(3,i) * axVx(p) + gxc(4,i) * dxVx * dt)
                azVz(p) = real(gzc(3,k) * azVz(p) + gzc(4,k) * dzVz * dt)

            end do
        end do
#ifdef _OPENACC
        !$acc end kernels
#else        
        !$omp end do nowait
        !$omp end parallel
#endif


#ifdef _OPENACC
        !$acc kernels &
        !$acc present(Vx, Vz, Sxz, azVx, axVz, gxc, gxe, gzc, gze, mu, bb, &
        !$acc         c1, c2, d1, d2, taus, kfs_top, kfs_bot, kob_top, kob_bot, Rxz)
        !$acc loop independent collapse(2)
#else
        !$omp parallel &
        !$omp private( dzVx, dzVz ) &
        !$omp private( dxVz_ade, dzVx_ade ) &
        !$omp private( i, k, m, p ) &
        !$omp private( taus1, taus_plus1, dxVz_dzVx, f_Rxz,  Rxz_n ) &
        !$omp private( re40x, re41x, re40z, re41z, isign) &
        !$omp private( nnn, pnn, npn, ppn, mu_xz) 
        !$omp do &
        !$omp schedule(dynamic)
#endif
        do i = bb%ib, bb%ie

            !ocl unroll('full')
            !ocl swp
            !OCL SWP_IREG_RATE(200)
            do k = bb%kb, bb%ke
 
                p = bb%offset + (i-bb%ib) * bb%nz + (k - bb%kb + 1)

                isign = sign(1, max((k - kfs_top(i)) * (kfs_bot(i) - k), &
                                    (k - kob_top(i)) * (kob_bot(i) - k)))

                re40x = rc40x + isign * rd40x
                re41x = rc41x + isign * rd41x
                re40z = rc40z + isign * rd40z
                re41z = rc41z + isign * rd41z               

                dxVz = (Vz(k  ,i+1) - Vz(k  ,i  )) * re40x - (Vz(k  ,i+2) - Vz(k  ,i-1)) * re41x
                dzVx = (Vx(k+1,i  ) - Vx(k  ,i  )) * re40z - (Vx(k+2,i  ) - Vx(k-1,i  )) * re41z

                taus1 = taus(k, i)

                nnn = mu(k,i)
                pnn = mu(k + 1, i)
                npn = mu(k,i + 1)
                ppn = mu(k + 1, i + 1)
                mu_xz = 4 * nnn * pnn * npn * ppn &
                     / (nnn * pnn * npn + nnn * pnn * ppn + nnn * npn * ppn + pnn * npn * ppn + epsl)

                dxVz_ade = gxe(1,i) * dxVz + gxe(2,i) * axVz(p)
                dzVx_ade = gze(1,k) * dzVx + gze(2,k) * azVx(p)
                dxVz_dzVx = dxVz_ade + dzVx_ade

                f_Rxz = mu_xz * taus1 * dxVz_dzVx

                Rxz_n = 0.0

                !! Crank-Nicolson Method for avoiding stiff solution
                !$acc loop seq reduction(+:Rxz_n)
                do m = 1, nm
                    Rxz(m, k, i) = c1(m) * Rxz(m, k, i) - c2(m) * f_Rxz * dt
                    Rxz_n = Rxz_n + d1(m) * Rxz(m,k,i)
                end do
                
                taus_plus1 = 1 + taus1 * (1 + d2)

                Sxz(k,i) = Sxz(k,i) + (mu_xz  * taus_plus1 * dxVz_dzVx + Rxz_n) * dt

                azVx(p) = real(gze(3,k) * azVx(p) + gze(4,k) * dzVx * dt)
                axVz(p) = real(gxe(3,i) * axVz(p) + gxe(4,i) * dxVz * dt)

            end do
        end do
#ifdef _OPENACC
        !$acc end kernels
#else
        !$omp end do nowait
        !$omp end parallel
        !$omp barrier
#endif

    end subroutine update_stress_core

    subroutine damping_profile(x, H, xb, xe, g)

        !! ADE-CFS PML damping factor according to Zhao and Shen

        real(SP), intent(in) :: x   !! cartesian coordinate location
        real(SP), intent(in) :: H   !! absorption layer thickness
        real(SP), intent(in) :: xb
        real(SP), intent(in) :: xe
        real(SP), intent(out) :: g(4) !< damping prof

        real(SP) :: R0 !! reflection coefficient
        real(SP) :: d0, a0, b0
        integer, parameter :: pd = 2
        integer, parameter :: pa = 1
        integer, parameter :: pb = 2
        real(SP), parameter :: cp = 6.0 !! assumed P-wave velocity
        real :: d, a, b, xx

        R0 = 10**(-(log10(real(na)) - 1) / log10(2.0) - 3.0)
        d0 = -(1.0 / (2.0 * H)) * (pd + 1) * cp * log(R0)
        b0 = 7.0
        a0 = real(PI * fcut)

        if (x <= xb + H) then
            xx = (xb + H) - x
        else if (x >= xe - H) then
            xx = x - (xe - H)
        else
            xx = 0.0 !! no absorption
        end if

        d = d0 * abs(xx / H)**pd
        a = a0 * (1.0 - abs(xx / H)**pa)
        b = 1.0 + (b0 - 1.0) * abs(xx / H)**pb

        g(1) = ((1.0 + (dt / 2.0) * a) / b) / (1.0 + (dt / 2.0) * (a + d / b))
        g(2) = (-1.0 / b) / (1.0 + (dt / 2.0) * (a + d / b))
        g(3) = (1.0 - (dt / 2.0) * (a + d / b)) / (1.0 + (dt / 2.0) * (a + d / b))
        g(4) = (d / b) / (1.0 + (dt / 2.0) * (a + d / b))

    end subroutine damping_profile

end module m_absorb_p
