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
    real(SP), allocatable :: gzc(:, :), gze(:, :) !< damping profile along x at center/edge of voxel

    real(SP), allocatable :: axVy(:), azVy(:)
    real(SP), allocatable :: axSxy(:), azSyz(:)

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

        hx = na * dx
        hz = na * dz

        do i = ibeg, iend
            call damping_profile(xc(i), hx, xbeg, xend, gxc(:, i))
            call damping_profile(xc(i) + real(dx) / 2.0, hx, xbeg, xend, gxe(:, i))
        end do
        do k = kbeg, kend
            call damping_profile(zc(k), hz, zbeg, zend, gzc(:, k))
            call damping_profile(zc(k) + real(dz) / 2.0, hz, zbeg, zend, gze(:, k))
        end do

        !! memory allocation
        allocate (axVy (n_sponge_cell), source=0.0)
        allocate (azVy (n_sponge_cell), source=0.0)
        allocate (axSxy(n_sponge_cell), source=0.0)
        allocate (azSyz(n_sponge_cell), source=0.0)

        idum = io_prm

        !$acc enter data copyin(axVy, azVy, axSxy, azSyz, gxc, gxe, gzc, gze)

    end subroutine absorb_p__setup


    subroutine absorb_p__update_vel
        !! Update velocity component in PML layer
        integer :: ibox

        !! Horizontal zero-derivative boundary (for plane wave mode)
        if (pw_mode) call set_stress_boundary()

        do ibox=1, 4
            if(box(ibox)%ncell > 0) call update_vel_core(box(ibox))
        end do

    end subroutine absorb_p__update_vel

    subroutine set_stress_boundary()

        integer :: i, k

        if (idx == 0) then
#ifdef _OPENACC
            !$acc kernels present(Syz, Sxy)                
            !$acc loop independent
#else
            !$omp parallel private(k)
            !$omp do schedule(dynamic)
#endif
            do k = kbeg, kend
                Syz(k,0) = 2 * Syz(k,1) - Syz(k,2)
                Sxy(k,0) = 2 * Sxy(k,1) - Sxy(k,2)
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
            !$acc kernels present(Sxy, Syz)
            !$acc loop independent
#else
            !$omp parallel private(k)
            !$omp do schedule(dynamic)
#endif
            do k = kbeg, kend
                Syz(k,nx+1) = 2 * Syz(k,nx) - Syz(k,nx-1)
                Sxy(k,nx+1) = 2 * Sxy(k,nx) - Sxy(k,nx-1)
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
        real(SP) :: by
        real(MP) :: dzSyz, dxSxy
        integer :: isign
        real(MP) :: re40x, re41x, re40z, re41z

#ifdef _OPENACC
        !$acc kernels &
        !$acc present(Vy, Sxy, Syz, axSxy, azSyz, rho, gxc, gzc, bb, kfs_top, kfs_bot, kob_top, kob_bot) 
        !$acc loop independent
#else
        !$omp parallel &
        !$omp private( dzSyz, dxSxy, by, i, k, p, isign, re40x, re41x, re40z, re41z)
        !$omp do schedule(dynamic)
#endif
        do i = bb%ib, bb%ie

#ifdef _OPENACC
            !$acc loop vector independent
#endif
            do k = bb%kb, bb%ke

                isign = sign(1, max((k - kfs_top(i)) * (kfs_bot(i) - k), &
                                    (k - kob_top(i)) * (kob_bot(i) - k)))

                re40x = rc40x + isign * rd40x
                re41x = rc41x + isign * rd41x
                re40z = rc40z + isign * rd40z
                re41z = rc41z + isign * rd41z               

                dzSyz = (Syz(k, i) - Syz(k-1,i  )) * re40z - (Syz(k+1,i  ) - Syz(k-2,i  )) * re41z
                dxSxy = (Sxy(k, i) - Sxy(k  ,i-1)) * re40x - (Sxy(k  ,i+1) - Sxy(k  ,i-2)) * re41x                

                p = bb%offset + (i-bb%ib) * bb%nz + (k - bb%kb + 1)

                by = 1.0 / rho(k,i)

                Vy(k,i) = Vy(k,i) &
                        + by * (gxc(1,i) * dxSxy    + gzc(1,k) * dzSyz &
                              + gxc(2,i) * axSxy(p) + gzc(2,k) * azSyz(p)) * dt

                axSxy(p) = gxc(3,i) * axSxy(p) + gxc(4,i) * dxSxy * dt
                azSyz(p) = gzc(3,k) * azSyz(p) + gzc(4,k) * dzSyz * dt

            end do
        end do
#ifdef _OPENACC
        !$acc end kernels
#else
        !$omp end do nowait
        !$omp end parallel
#endif        

    end subroutine update_vel_core

    subroutine absorb_p__update_stress

        integer :: ibox

        !! Horizontal zero-derivative boundary (for plane wave mode)
        if (pw_mode) call set_vel_boundary()

        do ibox=1, 4
            if( box(ibox)%ncell > 0 ) call update_stress_core(box(ibox))
        end do


    end subroutine absorb_p__update_stress

    subroutine set_vel_boundary()

        integer :: i, k

        if (idx == 0) then
#ifdef _OPENACC
            !$acc kernels present(Vy)
            !$acc loop independent
#else
            !$omp parallel private(k)
            !$omp do schedule(dynamic)
#endif
            do k = kbeg, kend
                Vy(k,0) = 2 * Vy(k,1) - Vy(k,2)
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
            !$acc kernels present(Vy)
            !$acc loop independent
#else
            !$omp parallel private(k)
            !$omp do schedule(dynamic)
#endif
            do k = kbeg, kend
                Vy(k,nx+1) = 2 * Vy(k,nx) - Vy(k,nx-1)
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
        integer :: i, k, m
        real(SP) :: nnn, pnn, npn
        real(SP) :: mu_xy, mu_yz
        real(SP) :: epsl = epsilon(1.0)
        real(MP) :: dxVy, dzVy, dxVy_ade, dzVy_ade
        integer :: p, p0
        real(SP) :: taus1, taus_plus1
        real(SP) :: Ryz_n, Rxy_n
        real(SP) :: f_Rxy, f_Ryz
        integer :: isign
        real(MP) :: re40x, re41x, re40z, re41z

#ifdef _OPENACC
        !$acc kernels &
        !$acc present(Vy, Sxy, Syz, gxc, gxe, gzc, gze, axVy, azVy, mu, bb, Ryz, Rxy, c1, c2, d1, d2, taus, kob_top, kfs_top, kob_bot, kfs_bot)
        !$acc loop independent
#else
        !$omp parallel &
        !$omp private(i, k, dxVy, dzVy, dxVy_ade, dzVy_ade, nnn, pnn, npn, mu_xy, mu_yz, p, p0 ) &
        !$omp private( taus1, taus_plus1, f_Ryz, f_Rxy,  Ryz_n, Rxy_n )  &
        !$omp private( re40x, re41x, re40z, re41z, isign) &
        !$omp do schedule(dynamic)
#endif
        do i = bb%ib, bb%ie

            p0 = bb%offset + (i-bb%ib) * bb%nz

            !$acc loop vector independent
            do k = bb%kb, bb%ke

                p = p0 + (k - bb%kb + 1)

                isign = sign(1, max((k - kfs_top(i)) * (kfs_bot(i) - k), &
                                    (k - kob_top(i)) * (kob_bot(i) - k)))

                re40x = rc40x + isign * rd40x
                re41x = rc41x + isign * rd41x
                re40z = rc40z + isign * rd40z
                re41z = rc41z + isign * rd41z      

                dxVy = (Vy(k  ,i+1) - Vy(k  ,i  )) * re40x - (Vy(k  ,i+2) - Vy(k  ,i-1)) * re41x
                dzVy = (Vy(k+1,i  ) - Vy(k  ,i  )) * re40z - (Vy(k+2,i  ) - Vy(k-1,i  )) * re41z

                nnn = mu(k  ,i  )
                pnn = mu(k+1,i  )
                npn = mu(k  ,i+1)
                mu_xy = 2 * nnn * npn / (nnn + npn + epsl)
                mu_yz = 2 * nnn * pnn / (nnn + pnn + epsl)

                dzVy_ade = gze(1,k) * dzVy + gze(2,k) * azVy(p)
                dxVy_ade = gxe(1,i) * dxVy + gxe(2,i) * axVy(p)

                !! update memory variables
                !! working variables for combinations of velocity derivatives
                taus1 = taus(k, i)
                f_Ryz = mu_yz * taus1 * dzVy_ade
                f_Rxy = mu_xy * taus1 * dxVy_ade

                Ryz_n = 0.0
                Rxy_n = 0.0

                !! Crank-Nicolson Method for avoiding stiff solution
                !$acc loop seq reduction(+:Ryz_n,Rxy_n)
                do m = 1, nm
                    Ryz(m,k,i) = c1(m) * Ryz(m,k,i) - c2(m) * f_Ryz * dt
                    Rxy(m,k,i) = c1(m) * Rxy(m,k,i) - c2(m) * f_Rxy * dt
                    Ryz_n = Ryz_n + d1(m) * Ryz(m,k,i)
                    Rxy_n = Rxy_n + d1(m) * Rxy(m,k,i)
                end do

                !! update stress components
                taus_plus1 = 1 + taus1 * (1 + d2)

                Syz(k,i) = Syz(k,i) + (mu_yz * taus_plus1 * dzVy_ade + Ryz_n) * dt
                Sxy(k,i) = Sxy(k,i) + (mu_xy * taus_plus1 * dxVy_ade + Rxy_n) * dt

                axVy(p) = gxe(3,i) * axVy(p) + gxe(4,i) * dxVy * dt
                azVy(p) = gze(3,k) * azVy(p) + gze(4,k) * dzVy * dt

            end do
        end do
#ifdef _OPENACC
        !$acc end kernels
#else
        !$omp end do nowait
        !$omp end parallel
#endif

        
    end subroutine update_stress_core


    subroutine damping_profile(x, H, xbeg, xend, g)

        !! ADE-CFS PML damping factor according to Zhao and Shen

        real(SP), intent(in) :: x   !< cartesian coordinate location
        real(SP), intent(in) :: H   !< absorption layer thickness
        real(SP), intent(in) :: xbeg
        real(SP), intent(in) :: xend
        real(SP), intent(out) :: g(4) !< damping prof

        real(SP) :: R0 !! reflection coefficient
        real(SP) :: d0, a0, b0
        integer, parameter :: pd = 1
        integer, parameter :: pa = 1
        integer, parameter :: pb = 2
        real(SP), parameter :: cp = 3.0 !! assumed S-wave velocity
        real :: d, a, b, xx

        R0 = 10**(-(log10(real(na)) - 1) / log10(2.0) - 3.0)
        d0 = -(1.0 / (2.0 * H)) * (pd + 1) * cp * log(R0)
        b0 = 7.0
        a0 = PI * fcut

        if (x <= xbeg + H) then
            xx = (xbeg + H) - x
        else if (x >= xend - H) then
            xx = x - (xend - H)
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
