#include "../shared/m_debug.h"
module m_absorb_p

    !! Absorbing Boundary Condition: ADE-CFS PML based on Zhang and Shen
    !!
    !! #### PML region definition
    !!
    !! ```text
    !!  +-----+--------------------------+-----+
    !!  |     |                          |     |
    !!  |     |                          |     |
    !!  |     |                          |     |
    !!  |     |                          |     |
    !!  |     |  interior region         |     |
    !!  |     |  eveluated by m_kernel   |     |
    !!  |     |                          |     |
    !!  |     |                          |     |
    !!  |     |                          |     |
    !!  |     +--------------------------+     |
    !!  |         exterior region              |
    !!  |         evaluated by m_absorb        |
    !!  +-----+--------------------------+-----+
    !!  1      na                        nx-na+1  nx
    !!  <- na ->
    !! ```
    !!
    !! Copyright 2013-2025 Takuto Maeda. All rights reserved. This project is released under the MIT license.

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

    real(SP), allocatable :: gxc(:, :), gxe(:, :) !! damping profile along x at center/edge of voxel
    real(SP), allocatable :: gyc(:, :), gye(:, :) !! damping profile along y at center/edge of voxel
    real(SP), allocatable :: gzc(:, :), gze(:, :) !! damping profile along z at center/edge of voxel

    real(SP), allocatable :: axVx(:, :, :), ayVx(:, :, :), azVx(:, :, :)
    real(SP), allocatable :: axVy(:, :, :), ayVy(:, :, :), azVy(:, :, :)
    real(SP), allocatable :: axVz(:, :, :), ayVz(:, :, :), azVz(:, :, :)
    real(SP), allocatable :: axSxx(:, :, :), aySxy(:, :, :), azSxz(:, :, :)
    real(SP), allocatable :: axSxy(:, :, :), aySyy(:, :, :), azSyz(:, :, :)
    real(SP), allocatable :: axSxz(:, :, :), aySyz(:, :, :), azSzz(:, :, :)

    real(MP) :: r20x, r20y, r20z

    integer :: kbeg_min

contains

    subroutine absorb_p__setup(io_prm)

        !! Set PML sponge

        integer, intent(in) :: io_prm
        integer  :: i, j, k
        real(SP) :: hx, hy, hz
        integer  :: idum

        !! derivative coefficient
        r20x = 1.0 / dx
        r20y = 1.0 / dy
        r20z = 1.0 / dz

        !! damping profile
        allocate (gxc(4, ibeg:iend), gxe(4, ibeg:iend))
        allocate (gyc(4, jbeg:jend), gye(4, jbeg:jend))
        allocate (gzc(4, kbeg:kend), gze(4, kbeg:kend))

        hx = na * real(dx)
        hy = na * real(dy)
        hz = na * real(dz)
        do i = ibeg, iend
            call damping_profile(xc(i), hx, xbeg, xend, gxc(:, i))
            call damping_profile(xc(i) + real(dx) / 2.0, hx, xbeg, xend, gxe(:, i))
        end do
        do j = jbeg, jend
            call damping_profile(yc(j), hy, ybeg, yend, gyc(:, j))
            call damping_profile(yc(j) + real(dy) / 2.0, hy, ybeg, yend, gye(:, j))
        end do
        do k = kbeg, kend
            call damping_profile(zc(k), hz, zbeg, zend, gzc(:, k))
            call damping_profile(zc(k) + real(dz) / 2.0, hz, zbeg, zend, gze(:, k))
        end do

        !! PML region definition
        kbeg_min = minval(kbeg_a(:, :))

        !! memory allocation
        allocate (axVx(kbeg_min:kend, ibeg:iend, jbeg:jend), source=0.0)
        allocate (ayVx(kbeg_min:kend, ibeg:iend, jbeg:jend), source=0.0)
        allocate (azVx(kbeg_min:kend, ibeg:iend, jbeg:jend), source=0.0)
        allocate (axVy(kbeg_min:kend, ibeg:iend, jbeg:jend), source=0.0)
        allocate (ayVy(kbeg_min:kend, ibeg:iend, jbeg:jend), source=0.0)
        allocate (azVy(kbeg_min:kend, ibeg:iend, jbeg:jend), source=0.0)
        allocate (axVz(kbeg_min:kend, ibeg:iend, jbeg:jend), source=0.0)
        allocate (ayVz(kbeg_min:kend, ibeg:iend, jbeg:jend), source=0.0)
        allocate (azVz(kbeg_min:kend, ibeg:iend, jbeg:jend), source=0.0)
        allocate (axSxx(kbeg_min:kend, ibeg:iend, jbeg:jend), source=0.0)
        allocate (aySxy(kbeg_min:kend, ibeg:iend, jbeg:jend), source=0.0)
        allocate (azSxz(kbeg_min:kend, ibeg:iend, jbeg:jend), source=0.0)
        allocate (axSxy(kbeg_min:kend, ibeg:iend, jbeg:jend), source=0.0)
        allocate (aySyy(kbeg_min:kend, ibeg:iend, jbeg:jend), source=0.0)
        allocate (azSyz(kbeg_min:kend, ibeg:iend, jbeg:jend), source=0.0)
        allocate (axSxz(kbeg_min:kend, ibeg:iend, jbeg:jend), source=0.0)
        allocate (aySyz(kbeg_min:kend, ibeg:iend, jbeg:jend), source=0.0)
        allocate (azSzz(kbeg_min:kend, ibeg:iend, jbeg:jend), source=0.0)

        idum = io_prm

    end subroutine absorb_p__setup

    subroutine absorb_p__update_vel

        !! Update velocity component in PML layer

        integer :: i, j, k
        real(SP) :: gxc0(4), gxe0(4), gyc0(4), gye0(4), gzc0(4), gze0(4)
        real(MP) :: dxSxx(kbeg:kend), dySxy(kbeg:kend), dzSxz(kbeg:kend)
        real(MP) :: dxSxy(kbeg:kend), dySyy(kbeg:kend), dzSyz(kbeg:kend)
        real(MP) :: dxSxz(kbeg:kend), dySyz(kbeg:kend), dzSzz(kbeg:kend)

        !! Horizontal zero-derivative boundary (for plane wave mode)
        if (pw_mode) then
            if (idx == 0) then
                !$omp parallel private(j,k)
                !$omp do schedule(dynamic)
                do j = jbeg, jend
                    do k = kbeg_a(1, j), kend
                        Sxx(k, 0, j) = 2 * Sxx(k, 1, j) - Sxx(k, 2, j)
                        Syy(k, 0, j) = 2 * Syy(k, 1, j) - Syy(k, 2, j)
                        Szz(k, 0, j) = 2 * Szz(k, 1, j) - Szz(k, 2, j)
                        Syz(k, 0, j) = 2 * Syz(k, 1, j) - Syz(k, 2, j)
                        Sxz(k, 0, j) = 2 * Sxz(k, 1, j) - Sxz(k, 2, j)
                        Sxy(k, 0, j) = 2 * Sxy(k, 1, j) - Sxy(k, 2, j)
                    end do
                end do
                !$omp end do nowait
                !$omp end parallel
            end if

            if (idx == nproc_x - 1) then
                !$omp parallel private(j,k)
                !$omp do schedule(dynamic)
                do j = jbeg, jend
                    do k = kbeg_a(nx, j), kend
                        Sxx(k, nx + 1, j) = 2 * Sxx(k, nx, j) - Sxx(k, nx - 1, j)
                        Syy(k, nx + 1, j) = 2 * Syy(k, nx, j) - Syy(k, nx - 1, j)
                        Szz(k, nx + 1, j) = 2 * Szz(k, nx, j) - Szz(k, nx - 1, j)
                        Syz(k, nx + 1, j) = 2 * Syz(k, nx, j) - Syz(k, nx - 1, j)
                        Sxz(k, nx + 1, j) = 2 * Sxz(k, nx, j) - Sxz(k, nx - 1, j)
                        Sxy(k, nx + 1, j) = 2 * Sxy(k, nx, j) - Sxy(k, nx - 1, j)
                    end do
                end do
                !$omp end do nowait
                !$omp end parallel
            end if

            if (idy == 0) then
                !$omp parallel private(i,k)
                !$omp do schedule(dynamic)
                do i = ibeg, iend
                    do k = kbeg_a(i, 1), kend
                        Sxx(k, i, 0) = 2 * Sxx(k, i, 1) - Sxx(k, i, 2)
                        Syy(k, i, 0) = 2 * Syy(k, i, 1) - Syy(k, i, 2)
                        Szz(k, i, 0) = 2 * Szz(k, i, 1) - Szz(k, i, 2)
                        Syz(k, i, 0) = 2 * Syz(k, i, 1) - Syz(k, i, 2)
                        Sxz(k, i, 0) = 2 * Sxz(k, i, 1) - Sxz(k, i, 2)
                        Sxy(k, i, 0) = 2 * Sxy(k, i, 1) - Sxy(k, i, 2)
                    end do
                end do
                !$omp end do nowait
                !$omp end parallel
            end if

            if (idy == nproc_y - 1) then
                !$omp parallel private(i,k)
                !$omp do schedule(dynamic)
                do i = ibeg, iend
                    do k = kbeg_a(i, ny), kend
                        Sxx(k, i, ny + 1) = 2 * Sxx(k, i, ny) - Sxx(k, i, ny - 1)
                        Syy(k, i, ny + 1) = 2 * Syy(k, i, ny) - Syy(k, i, ny - 1)
                        Szz(k, i, ny + 1) = 2 * Szz(k, i, ny) - Szz(k, i, ny - 1)
                        Syz(k, i, ny + 1) = 2 * Syz(k, i, ny) - Syz(k, i, ny - 1)
                        Sxz(k, i, ny + 1) = 2 * Sxz(k, i, ny) - Sxz(k, i, ny - 1)
                        Sxy(k, i, ny + 1) = 2 * Sxy(k, i, ny) - Sxy(k, i, ny - 1)
                    end do
                end do
                !$omp end do nowait
                !$omp end parallel
            end if
            !$omp barrier
        end if

        !! time-marching

        !$omp parallel &
        !$omp private( dxSxx, dySyy, dzSzz, dySyz, dzSyz, dxSxz, dzSxz, dxSxy ,dySxy ) &
        !$omp private( gxc0, gxe0, gyc0, gye0, gzc0, gze0 ) &
        !$omp private( i, j, k )
        !$omp do &
        !$omp schedule(dynamic)
        do j = jbeg, jend

            gyc0(1:4) = gyc(1:4, j)
            gye0(1:4) = gye(1:4, j)

            do i = ibeg, iend

                gxc0(1:4) = gxc(1:4, i)
                gxe0(1:4) = gxe(1:4, i)

                !! Derivatives
                do k = kbeg_a(i, j), kend

                    dxSxx(k) = (Sxx(k, i + 1, j) - Sxx(k, i, j)) * r20x
                    dySyy(k) = (Syy(k, i, j + 1) - Syy(k, i, j)) * r20y
                    dzSzz(k) = (Szz(k + 1, i, j) - Szz(k, i, j)) * r20z
                    dySyz(k) = (Syz(k, i, j) - Syz(k, i, j - 1)) * r20y
                    dzSyz(k) = (Syz(k, i, j) - Syz(k - 1, i, j)) * r20z
                    dxSxz(k) = (Sxz(k, i, j) - Sxz(k, i - 1, j)) * r20x
                    dzSxz(k) = (Sxz(k, i, j) - Sxz(k - 1, i, j)) * r20z
                    dxSxy(k) = (Sxy(k, i, j) - Sxy(k, i - 1, j)) * r20x
                    dySxy(k) = (Sxy(k, i, j) - Sxy(k, i, j - 1)) * r20y

                end do

                !! update velocity
                do k = kbeg_a(i, j), kend

                    gzc0(1:4) = gzc(1:4, k)
                    gze0(1:4) = gze(1:4, k)

                    !! Velocity Updates
                    Vx(k, i, j) = Vx(k, i, j) + bx(k, i, j) &
                                  * real(gxe0(1) * dxSxx(k) + gyc0(1) * dySxy(k) + gzc0(1) * dzSxz(k) &
                                         + gxe0(2) * axSxx(k, i, j) + gyc0(2) * aySxy(k, i, j) + gzc0(2) * azSxz(k, i, j)) * dt

                    Vy(k, i, j) = Vy(k, i, j) + by(k, i, j) &
                                  * real(gxc0(1) * dxSxy(k) + gye0(1) * dySyy(k) + gzc0(1) * dzSyz(k) &
                                         + gxc0(2) * axSxy(k, i, j) + gye0(2) * aySyy(k, i, j) + gzc0(2) * azSyz(k, i, j)) * dt

                    Vz(k, i, j) = Vz(k, i, j) + bz(k, i, j) &
                                  * real(gxc0(1) * dxSxz(k) + gyc0(1) * dySyz(k) + gze0(1) * dzSzz(k) &
                                         + gxc0(2) * axSxz(k, i, j) + gyc0(2) * aySyz(k, i, j) + gze0(2) * azSzz(k, i, j)) * dt

                    !! ADE updates
                    axSxx(k, i, j) = gxe0(3) * axSxx(k, i, j) + gxe0(4) * real(dxSxx(k)) * dt
                    aySxy(k, i, j) = gyc0(3) * aySxy(k, i, j) + gyc0(4) * real(dySxy(k)) * dt
                    azSxz(k, i, j) = gzc0(3) * azSxz(k, i, j) + gzc0(4) * real(dzSxz(k)) * dt
                    axSxy(k, i, j) = gxc0(3) * axSxy(k, i, j) + gxc0(4) * real(dxSxy(k)) * dt
                    aySyy(k, i, j) = gye0(3) * aySyy(k, i, j) + gye0(4) * real(dySyy(k)) * dt
                    azSyz(k, i, j) = gzc0(3) * azSyz(k, i, j) + gzc0(4) * real(dzSyz(k)) * dt
                    axSxz(k, i, j) = gxc0(3) * axSxz(k, i, j) + gxc0(4) * real(dxSxz(k)) * dt
                    aySyz(k, i, j) = gyc0(3) * aySyz(k, i, j) + gyc0(4) * real(dySyz(k)) * dt
                    azSzz(k, i, j) = gze0(3) * azSzz(k, i, j) + gze0(4) * real(dzSzz(k)) * dt

                end do
            end do
        end do
        !$omp end do nowait
        !$omp end parallel

        !$omp barrier

    end subroutine absorb_p__update_vel

    subroutine absorb_p__update_stress

        integer :: i, j, k
        real(SP) :: lam2mu_R, lam_R
        real(SP) :: dxVx_ade, dyVy_ade, dzVz_ade
        real(SP) :: gxc0(4), gxe0(4), gyc0(4), gye0(4), gzc0(4), gze0(4)
        real(MP) :: dxVx(kbeg_min:kend), dyVx(kbeg_min:kend), dzVx(kbeg_min:kend)
        real(MP) :: dxVy(kbeg_min:kend), dyVy(kbeg_min:kend), dzVy(kbeg_min:kend)
        real(MP) :: dxVz(kbeg_min:kend), dyVz(kbeg_min:kend), dzVz(kbeg_min:kend)

        !! Horizontal zero-derivative boundary (for plane wave mode)
        if (pw_mode) then
            if (idx == 0) then
                !$omp parallel private(j,k)
                !$omp do schedule(dynamic)
                do j = jbeg, jend
                    do k = kbeg_a(1, j), kend
                        Vx(k, 0, j) = 2 * Vx(k, 1, j) - Vx(k, 2, j)
                        Vy(k, 0, j) = 2 * Vy(k, 1, j) - Vy(k, 2, j)
                        Vz(k, 0, j) = 2 * Vz(k, 1, j) - Vz(k, 2, j)
                    end do
                end do
                !$omp end do nowait
                !$omp end parallel
            end if

            if (idx == nproc_x - 1) then
                !$omp parallel private(j,k)
                !$omp do schedule(dynamic)
                do j = jbeg, jend
                    do k = kbeg_a(nx, j), kend
                        Vx(k, nx + 1, j) = 2 * Vx(k, nx, j) - Vx(k, nx - 1, j)
                        Vy(k, nx + 1, j) = 2 * Vy(k, nx, j) - Vy(k, nx - 1, j)
                        Vz(k, nx + 1, j) = 2 * Vz(k, nx, j) - Vz(k, nx - 1, j)
                    end do
                end do
                !$omp end do nowait
                !$omp end parallel
            end if

            if (idy == 0) then
                !$omp parallel private(i,k)
                !$omp do schedule(dynamic)
                do i = ibeg, iend
                    do k = kbeg_a(i, 1), kend
                        Vx(k, i, 0) = 2 * Vx(k, i, 1) - Vx(k, i, 2)
                        Vy(k, i, 0) = 2 * Vy(k, i, 1) - Vy(k, i, 2)
                        Vz(k, i, 0) = 2 * Vz(k, i, 1) - Vz(k, i, 2)
                    end do
                end do
                !$omp end do nowait
                !$omp end parallel
            end if

            if (idy == nproc_y - 1) then
                !$omp parallel private(i,k)
                !$omp do schedule(dynamic)
                do i = ibeg, iend
                    do k = kbeg_a(i, ny), kend
                        Vx(k, i, ny + 1) = 2 * Vx(k, i, ny) - Vx(k, i, ny - 1)
                        Vy(k, i, ny + 1) = 2 * Vy(k, i, ny) - Vy(k, i, ny - 1)
                        Vz(k, i, ny + 1) = 2 * Vz(k, i, ny) - Vz(k, i, ny - 1)
                    end do
                end do
                !$omp end do nowait
                !$omp end parallel
            end if

            !$omp barrier
        end if

        !! Time-marching

        !$omp parallel &
        !$omp private( gxc0, gxe0, gyc0, gye0, gzc0, gze0 ) &
        !$omp private( dxVx, dxVy, dxVz, dyVx, dyVy, dyVz, dzVx, dzVy, dzVz ) &
        !$omp private( lam2mu_R, lam_R ) &
        !$omp private( dxVx_ade, dyVy_ade, dzVz_ade ) &
        !$omp private( i, j, k )
        !$omp do &
        !$omp schedule(dynamic)
        do j = jbeg, jend

            gyc0(1:4) = gyc(1:4, j)
            gye0(1:4) = gye(1:4, j)

            do i = ibeg, iend

                gxc0(1:4) = gxc(1:4, i)
                gxe0(1:4) = gxe(1:4, i)

                !! Derivatives
                do k = kbeg_a(i, j), kend

                    dxVx(k) = (Vx(k, i, j) - Vx(k, i - 1, j)) * r20x
                    dxVy(k) = (Vy(k, i + 1, j) - Vy(k, i, j)) * r20x
                    dxVz(k) = (Vz(k, i + 1, j) - Vz(k, i, j)) * r20x
                    dyVx(k) = (Vx(k, i, j + 1) - Vx(k, i, j)) * r20y
                    dyVy(k) = (Vy(k, i, j) - Vy(k, i, j - 1)) * r20y
                    dyVz(k) = (Vz(k, i, j + 1) - Vz(k, i, j)) * r20y
                    dzVx(k) = (Vx(k + 1, i, j) - Vx(k, i, j)) * r20z
                    dzVy(k) = (Vy(k + 1, i, j) - Vy(k, i, j)) * r20z
                    dzVz(k) = (Vz(k, i, j) - Vz(k - 1, i, j)) * r20z

                end do

                !! Update Normal Stress
                do k = kbeg_a(i, j), kend

                    gzc0(1:4) = gzc(1:4, k)

                    lam2mu_R = (lam(k, i, j) + 2 * mu(k, i, j))
                    lam_R = lam2mu_R - 2 * mu(k, i, j)

                    dxVx_ade = gxc0(1) * real(dxVx(k)) + gxc0(2) * axVx(k, i, j)
                    dyVy_ade = gyc0(1) * real(dyVy(k)) + gyc0(2) * ayVy(k, i, j)
                    dzVz_ade = gzc0(1) * real(dzVz(k)) + gzc0(2) * azVz(k, i, j)

                    Sxx(k, i, j) = Sxx(k, i, j) + (lam2mu_R * dxVx_ade + lam_R * (dyVy_ade + dzVz_ade)) * dt
                    Syy(k, i, j) = Syy(k, i, j) + (lam2mu_R * dyVy_ade + lam_R * (dxVx_ade + dzVz_ade)) * dt
                    Szz(k, i, j) = Szz(k, i, j) + (lam2mu_R * dzVz_ade + lam_R * (dxVx_ade + dyVy_ade)) * dt

                    axVx(k, i, j) = gxc0(3) * axVx(k, i, j) + gxc0(4) * real(dxVx(k)) * dt
                    ayVy(k, i, j) = gyc0(3) * ayVy(k, i, j) + gyc0(4) * real(dyVy(k)) * dt
                    azVz(k, i, j) = gzc0(3) * azVz(k, i, j) + gzc0(4) * real(dzVz(k)) * dt

                end do

                !! Update Shear Stress
                do k = kbeg_a(i, j), kend

                    gze0(1:4) = gze(1:4, k)

                    Syz(k, i, j) = Syz(k, i, j) + muyz(k, i, j) * (gye0(1) * dyVz(k) + gze0(1) * dzVy(k) &
                                                                   + gye0(2) * ayVz(k, i, j) + gze0(2) * azVy(k, i, j)) * dt

                    Sxz(k, i, j) = Sxz(k, i, j) + muxz(k, i, j) * (gxe0(1) * dxVz(k) + gze0(1) * dzVx(k) &
                                                                   + gxe0(2) * axVz(k, i, j) + gze0(2) * azVx(k, i, j)) * dt

                    Sxy(k, i, j) = Sxy(k, i, j) + muxy(k, i, j) * (gxe0(1) * dxVy(k) + gye0(1) * dyVx(k) &
                                                                   + gxe0(2) * axVy(k, i, j) + gye0(2) * ayVx(k, i, j)) * dt

                    ayVx(k, i, j) = gye0(3) * ayVx(k, i, j) + gye0(4) * real(dyVx(k)) * dt
                    azVx(k, i, j) = gze0(3) * azVx(k, i, j) + gze0(4) * real(dzVx(k)) * dt
                    axVy(k, i, j) = gxe0(3) * axVy(k, i, j) + gxe0(4) * real(dxVy(k)) * dt
                    azVy(k, i, j) = gze0(3) * azVy(k, i, j) + gze0(4) * real(dzVy(k)) * dt
                    axVz(k, i, j) = gxe0(3) * axVz(k, i, j) + gxe0(4) * real(dxVz(k)) * dt
                    ayVz(k, i, j) = gye0(3) * ayVz(k, i, j) + gye0(4) * real(dyVz(k)) * dt

                end do
            end do
        end do
        !$omp end do nowait
        !$omp end parallel
        !$omp barrier

    end subroutine absorb_p__update_stress

    subroutine damping_profile(x, H, xbeg0, xend0, g)

        !! ADE-CFS PML damping factor according to Zhao and Shen

        real(SP), intent(in) :: x   !! cartesian coordinate location
        real(SP), intent(in) :: H   !! absorption layer thickness
        real(SP), intent(in) :: xbeg0
        real(SP), intent(in) :: xend0
        real(SP), intent(out) :: g(4) !! damping prof

        real(SP) :: R0 !! reflection coefficient
        real(SP) :: d0, a0, b0
        integer, parameter :: pd = 1
        integer, parameter :: pa = 1
        integer, parameter :: pb = 2
        real(SP), parameter :: cp = 6.0 !! assumed P-wave velocity
        real :: d, a, b, xx

        R0 = 10**(-(log10(real(na)) - 1) / log10(2.0) - 3.0)
        d0 = -(1.0 / (2.0 * H)) * (pd + 1) * cp * log(R0)
        b0 = 7.0
        a0 = real(PI * fcut)

        if (x <= xbeg0 + H) then
            xx = (xbeg0 + H) - x
        else if (x >= xend0 - H) then
            xx = x - (xend0 - H)
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
