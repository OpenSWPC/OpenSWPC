#include "../shared/m_debug.h"
module m_absorb_p

    !! Absorbing Boundary Condition: ADE-CFS PML based on Zhang and Shen
    !!
    !!#### PML region definition
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

        !$acc enter data copyin(axVx, ayVx, azVx, axVy, ayVy, azVy, axVz, ayVz, azVz, &
        !$acc                   axSxx, aySxy, azSxz, axSxy, aySyy, azSyz, axSxz, aySyz, azSzz, &
        !$acc                   gxc, gxe, gyc, gye, gzc, gze)

    end subroutine absorb_p__setup

    subroutine absorb_p__update_vel

        !! Update velocity component in PML layer

        integer :: i, j, k
        real(MP) :: dxSxx, dySxy, dzSxz
        real(MP) :: dxSxy, dySyy, dzSyz
        real(MP) :: dxSxz, dySyz, dzSzz

        !! Horizontal zero-derivative boundary (for plane wave mode)
        if (pw_mode) then
            if (idx == 0) then
#ifdef _OPENACC
                !$acc kernels present(Sxx, Syy, Szz, Syz, Sxz, Sxy) 
                !$acc loop independent collapse(2)
#else
                !$omp parallel private(j,k)
                !$omp do schedule(dynamic)
#endif
                do j = jbeg, jend
                    do k = 1, nz
                        Sxx(k,0,j) = 2 * Sxx(k,1,j) - Sxx(k,2,j)
                        Syy(k,0,j) = 2 * Syy(k,1,j) - Syy(k,2,j)
                        Szz(k,0,j) = 2 * Szz(k,1,j) - Szz(k,2,j)
                        Syz(k,0,j) = 2 * Syz(k,1,j) - Syz(k,2,j)
                        Sxz(k,0,j) = 2 * Sxz(k,1,j) - Sxz(k,2,j)
                        Sxy(k,0,j) = 2 * Sxy(k,1,j) - Sxy(k,2,j)
                    end do
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
                !$acc kernels present(Sxx, Syy, Szz, Syz, Sxz, Sxy) 
                !$acc loop independent collapse(2)
#else
                !$omp parallel private(j,k)
                !$omp do schedule(dynamic)
#endif
                do j = jbeg, jend
                    do k = 1, nz
                        Sxx(k,nx+1,j) = 2 * Sxx(k,nx,j) - Sxx(k,nx-1,j)
                        Syy(k,nx+1,j) = 2 * Syy(k,nx,j) - Syy(k,nx-1,j)
                        Szz(k,nx+1,j) = 2 * Szz(k,nx,j) - Szz(k,nx-1,j)
                        Syz(k,nx+1,j) = 2 * Syz(k,nx,j) - Syz(k,nx-1,j)
                        Sxz(k,nx+1,j) = 2 * Sxz(k,nx,j) - Sxz(k,nx-1,j)
                        Sxy(k,nx+1,j) = 2 * Sxy(k,nx,j) - Sxy(k,nx-1,j)
                    end do
                end do
#ifdef _OPENACC
                !$acc end kernels
#else
                !$omp end do nowait
                !$omp end parallel
#endif
            end if

            if (idy == 0) then
#ifdef _OPENACC
                !$acc kernels present(Sxx, Syy, Szz, Syz, Sxz, Sxy) 
                !$acc loop independent collapse(2)
#else
                !$omp parallel private(i,k)
                !$omp do schedule(dynamic)
#endif
                do i = ibeg, iend
                    do k = 1, nz
                        Sxx(k, i, 0) = 2 * Sxx(k, i, 1) - Sxx(k, i, 2)
                        Syy(k, i, 0) = 2 * Syy(k, i, 1) - Syy(k, i, 2)
                        Szz(k, i, 0) = 2 * Szz(k, i, 1) - Szz(k, i, 2)
                        Syz(k, i, 0) = 2 * Syz(k, i, 1) - Syz(k, i, 2)
                        Sxz(k, i, 0) = 2 * Sxz(k, i, 1) - Sxz(k, i, 2)
                        Sxy(k, i, 0) = 2 * Sxy(k, i, 1) - Sxy(k, i, 2)
                    end do
                end do
#ifdef _OPENACC
                !$acc end kernels
#else
                !$omp end do nowait
                !$omp end parallel
#endif
            end if

            if (idy == nproc_y - 1) then
#ifdef _OPENACC
                !$acc kernels present(Sxx, Syy, Szz, Syz, Sxz, Sxy) 
                !$acc loop independent collapse(2)
#else
                !$omp parallel private(i,k)
                !$omp do schedule(dynamic)
#endif
                do i = ibeg, iend
                    do k = 1, nz
                        Sxx(k,i,ny+1) = 2 * Sxx(k,i,ny) - Sxx(k,i,ny-1)
                        Syy(k,i,ny+1) = 2 * Syy(k,i,ny) - Syy(k,i,ny-1)
                        Szz(k,i,ny+1) = 2 * Szz(k,i,ny) - Szz(k,i,ny-1)
                        Syz(k,i,ny+1) = 2 * Syz(k,i,ny) - Syz(k,i,ny-1)
                        Sxz(k,i,ny+1) = 2 * Sxz(k,i,ny) - Sxz(k,i,ny-1)
                        Sxy(k,i,ny+1) = 2 * Sxy(k,i,ny) - Sxy(k,i,ny-1)
                    end do
                end do
#ifdef _OPENACC
                !$acc end kernels
#else
                !$omp end do nowait
                !$omp end parallel
#endif
            end if
#ifndef _OPENACC
            !$omp barrier
#endif
        end if

        !! time-marching

#ifdef _OPENACC
        !$acc kernels &
        !$acc present(Vx, Vy, Vz, Sxx, Syy, Szz, Syz, Sxz, Sxy, &
        !$acc         axSxx, aySxy, azSxz, axSxy, aySyy, azSyz, axSxz, aySyz, azSzz, &
        !$acc         bx, by, bz, gxc, gxe, gyc, gye, gzc, gze, kbeg_a)
        !$acc loop independent collapse(2)
#else
        !$omp parallel &
        !$omp private( dxSxx, dySyy, dzSzz, dySyz, dzSyz, dxSxz, dzSxz, dxSxy ,dySxy ) &
        !$omp private( i, j, k )
        !$omp do &
        !$omp schedule(dynamic)
#endif
        do j = jbeg, jend
            do i = ibeg, iend

                !$acc loop vector independent
                do k = kbeg_a(i, j), kend

                    dxSxx = (Sxx(k  ,i+1,j) - Sxx(k  ,i  ,j  )) * r20x
                    dySyy = (Syy(k  ,i,j+1) - Syy(k  ,i  ,j  )) * r20y
                    dzSzz = (Szz(k+1,i,  j) - Szz(k  ,i  ,j  )) * r20z
                    dySyz = (Syz(k  ,i  ,j) - Syz(k  ,i  ,j-1)) * r20y
                    dzSyz = (Syz(k  ,i  ,j) - Syz(k-1,i  ,j  )) * r20z
                    dxSxz = (Sxz(k  ,i  ,j) - Sxz(k  ,i-1,j  )) * r20x
                    dzSxz = (Sxz(k  ,i  ,j) - Sxz(k-1,i  ,j  )) * r20z
                    dxSxy = (Sxy(k  ,i  ,j) - Sxy(k  ,i-1,j  )) * r20x
                    dySxy = (Sxy(k  ,i  ,j) - Sxy(k  ,i  ,j-1)) * r20y

                    !! Velocity Updates
                    Vx(k,i,j) = Vx(k,i,j) + bx(k,i,j) &
                              * real(gxe(1,i) * dxSxx + gyc(1,j) * dySxy + gzc(1,k) * dzSxz &
                                   + gxe(2,i) * axSxx(k,i,j) + gyc(2,j) * aySxy(k,i,j) + gzc(2,k) * azSxz(k,i,j)) * dt

                    Vy(k,i,j) = Vy(k,i,j) + by(k,i,j) &
                              * real(gxc(1,i) * dxSxy + gye(1,j) * dySyy + gzc(1,k) * dzSyz &
                                   + gxc(2,i) * axSxy(k,i,j) + gye(2,j) * aySyy(k,i,j) + gzc(2,k) * azSyz(k,i,j)) * dt

                    Vz(k,i,j) = Vz(k,i,j) + bz(k,i,j) &
                              * real(gxc(1,i) * dxSxz + gyc(1,j) * dySyz + gze(1,k) * dzSzz &
                                   + gxc(2,i) * axSxz(k,i,j) + gyc(2,j) * aySyz(k,i,j) + gze(2,k) * azSzz(k,i,j)) * dt

                    !! ADE updates
                    axSxx(k,i,j) = gxe(3,i) * axSxx(k,i,j) + gxe(4,i) * real(dxSxx) * dt
                    aySxy(k,i,j) = gyc(3,j) * aySxy(k,i,j) + gyc(4,j) * real(dySxy) * dt
                    azSxz(k,i,j) = gzc(3,k) * azSxz(k,i,j) + gzc(4,k) * real(dzSxz) * dt
                    axSxy(k,i,j) = gxc(3,i) * axSxy(k,i,j) + gxc(4,i) * real(dxSxy) * dt
                    aySyy(k,i,j) = gye(3,j) * aySyy(k,i,j) + gye(4,j) * real(dySyy) * dt
                    azSyz(k,i,j) = gzc(3,k) * azSyz(k,i,j) + gzc(4,k) * real(dzSyz) * dt
                    axSxz(k,i,j) = gxc(3,i) * axSxz(k,i,j) + gxc(4,i) * real(dxSxz) * dt
                    aySyz(k,i,j) = gyc(3,j) * aySyz(k,i,j) + gyc(4,j) * real(dySyz) * dt
                    azSzz(k,i,j) = gze(3,k) * azSzz(k,i,j) + gze(4,k) * real(dzSzz) * dt

                end do
            end do
        end do
#ifdef _OPENACC
        !$acc end kernels
#else
        !$omp end do nowait
        !$omp end parallel
        !$omp barrier
#endif 

    end subroutine absorb_p__update_vel

    subroutine absorb_p__update_stress

        integer :: i, j, k
        real(SP) :: lam2mu_R, lam_R
        real(SP) :: dxVx_ade, dyVy_ade, dzVz_ade
        real(SP) :: gxc0(4), gxe0(4), gyc0(4), gye0(4), gzc0(4), gze0(4)
        real(MP) :: dxVx, dyVx, dzVx
        real(MP) :: dxVy, dyVy, dzVy
        real(MP) :: dxVz, dyVz, dzVz

        !! Horizontal zero-derivative boundary (for plane wave mode)
        if (pw_mode) then
            if (idx == 0) then
#ifdef _OPENACC
                !$acc kernels present(Vx, Vy, Vz) 
                !$acc loop independent collapse(2)
#else
                !$omp parallel private(j,k)
                !$omp do schedule(dynamic)
#endif
                do j = jbeg, jend
                    do k = 1, nz
                        Vx(k, 0, j) = 2 * Vx(k, 1, j) - Vx(k, 2, j)
                        Vy(k, 0, j) = 2 * Vy(k, 1, j) - Vy(k, 2, j)
                        Vz(k, 0, j) = 2 * Vz(k, 1, j) - Vz(k, 2, j)
                    end do
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
                !$acc kernels present(Vx, Vy, Vz) 
                !$acc loop independent collapse(2)
#else
                !$omp parallel private(j,k)
                !$omp do schedule(dynamic)
#endif
                do j = jbeg, jend
                    do k = 1, nz
                        Vx(k, nx + 1, j) = 2 * Vx(k, nx, j) - Vx(k, nx - 1, j)
                        Vy(k, nx + 1, j) = 2 * Vy(k, nx, j) - Vy(k, nx - 1, j)
                        Vz(k, nx + 1, j) = 2 * Vz(k, nx, j) - Vz(k, nx - 1, j)
                    end do
                end do
#ifdef _OPENACC
                !$acc end kernels
#else
                !$omp end do nowait
                !$omp end parallel
#endif
            end if

            if (idy == 0) then
#ifdef _OPENACC
                !$acc kernels present(Vx, Vy, Vz)
                !$acc loop independent collapse(2)
#else
                !$omp parallel private(i,k)
                !$omp do schedule(dynamic)
#endif
                do i = ibeg, iend
                    do k = 1, nz
                        Vx(k, i, 0) = 2 * Vx(k, i, 1) - Vx(k, i, 2)
                        Vy(k, i, 0) = 2 * Vy(k, i, 1) - Vy(k, i, 2)
                        Vz(k, i, 0) = 2 * Vz(k, i, 1) - Vz(k, i, 2)
                    end do
                end do
#ifdef _OPENACC
                !$acc end kernels
#else
                !$omp end do nowait
                !$omp end parallel
#endif
            end if

            if (idy == nproc_y - 1) then
#ifdef _OPENACC
                !$acc kernels present(Vx, Vy, Vz) 
                !$acc loop independent collapse(2)
#else
                !$omp parallel private(i,k)
                !$omp do schedule(dynamic)
#endif
                do i = ibeg, iend
                    do k = 1, nz
                        Vx(k, i, ny + 1) = 2 * Vx(k, i, ny) - Vx(k, i, ny - 1)
                        Vy(k, i, ny + 1) = 2 * Vy(k, i, ny) - Vy(k, i, ny - 1)
                        Vz(k, i, ny + 1) = 2 * Vz(k, i, ny) - Vz(k, i, ny - 1)
                    end do
                end do
#ifdef _OPENACC
                !$acc end kernels
#else
                !$omp end do nowait
                !$omp end parallel
#endif
            end if
#ifndef _OPENACC
            !$omp barrier
#endif 
        end if

        !! Time-marching

#ifdef _OPENACC
        !$acc kernels &
        !$acc present(Vx, Vy, Vz, Sxx, Syy, Szz, Syz, Sxz, Sxy, &
        !$acc         axVx, ayVx, azVx, axVy, ayVy, azVy, axVz, ayVz, azVz, &
        !$acc         lam, mu, muyz, muxz, muxy, gxc, gxe, gyc, gye, gzc, gze, kbeg_a)
        !$acc loop independent collapse(2)
#else
        !$omp parallel &
        !$omp private( dxVx, dxVy, dxVz, dyVx, dyVy, dyVz, dzVx, dzVy, dzVz ) &
        !$omp private( lam2mu_R, lam_R ) &
        !$omp private( dxVx_ade, dyVy_ade, dzVz_ade ) &
        !$omp private( i, j, k )
        !$omp do &
        !$omp schedule(dynamic)
#endif
        do j = jbeg, jend

            do i = ibeg, iend

                !$acc loop vector independent
                do k = kbeg_a(i, j), kend

                    dxVx = (Vx(k  ,i  ,j  ) - Vx(k  ,i-1,j  )) * r20x
                    dyVy = (Vy(k  ,i  ,j  ) - Vy(k  ,i  ,j-1)) * r20y
                    dzVz = (Vz(k  ,i  ,j  ) - Vz(k-1,i  ,j  )) * r20z
    
                    lam2mu_R = (lam(k,i,j) + 2 * mu(k,i,j))
                    lam_R = lam2mu_R - 2 * mu(k,i,j)

                    dxVx_ade = gxc(1,i) * real(dxVx) + gxc(2,i) * axVx(k,i,j)
                    dyVy_ade = gyc(1,j) * real(dyVy) + gyc(2,j) * ayVy(k,i,j)
                    dzVz_ade = gzc(1,k) * real(dzVz) + gzc(2,k) * azVz(k,i,j)

                    Sxx(k,i,j) = Sxx(k,i,j) + (lam2mu_R * dxVx_ade + lam_R * (dyVy_ade + dzVz_ade)) * dt
                    Syy(k,i,j) = Syy(k,i,j) + (lam2mu_R * dyVy_ade + lam_R * (dxVx_ade + dzVz_ade)) * dt
                    Szz(k,i,j) = Szz(k,i,j) + (lam2mu_R * dzVz_ade + lam_R * (dxVx_ade + dyVy_ade)) * dt

                    axVx(k,i,j) = gxc(3,i) * axVx(k,i,j) + gxc(4,i) * real(dxVx) * dt
                    ayVy(k,i,j) = gyc(3,j) * ayVy(k,i,j) + gyc(4,j) * real(dyVy) * dt
                    azVz(k,i,j) = gzc(3,k) * azVz(k,i,j) + gzc(4,k) * real(dzVz) * dt

                end do

                !$acc loop vector independent
                do k = kbeg_a(i, j), kend

                    dxVy = (Vy(k  ,i+1,j  ) - Vy(k  ,i  ,j  )) * r20x
                    dxVz = (Vz(k  ,i+1,j  ) - Vz(k  ,i  ,j  )) * r20x
                    dyVx = (Vx(k  ,i  ,j+1) - Vx(k  ,i  ,j  )) * r20y
                    dyVz = (Vz(k  ,i  ,j+1) - Vz(k  ,i  ,j  )) * r20y
                    dzVx = (Vx(k+1,i  ,j  ) - Vx(k  ,i  ,j  )) * r20z
                    dzVy = (Vy(k+1,i  ,j  ) - Vy(k  ,i  ,j  )) * r20z


                    Syz(k,i,j) = Syz(k,i,j) + muyz(k,i,j) * (gye(1,j) * dyVz + gze(1,k) * dzVy &
                                                                   + gye(2,j) * ayVz(k,i,j) + gze(2,k) * azVy(k,i,j)) * dt

                    Sxz(k,i,j) = Sxz(k,i,j) + muxz(k,i,j) * (gxe(1,i) * dxVz + gze(1,k) * dzVx &
                                                                   + gxe(2,i) * axVz(k,i,j) + gze(2,k) * azVx(k,i,j)) * dt

                    Sxy(k,i,j) = Sxy(k,i,j) + muxy(k,i,j) * (gxe(1,i) * dxVy + gye(1,j) * dyVx &
                                                                   + gxe(2,i) * axVy(k,i,j) + gye(2,j) * ayVx(k,i,j)) * dt

                    ayVx(k,i,j) = gye(3,j) * ayVx(k,i,j) + gye(4,j) * real(dyVx) * dt
                    azVx(k,i,j) = gze(3,k) * azVx(k,i,j) + gze(4,k) * real(dzVx) * dt
                    axVy(k,i,j) = gxe(3,i) * axVy(k,i,j) + gxe(4,i) * real(dxVy) * dt
                    azVy(k,i,j) = gze(3,k) * azVy(k,i,j) + gze(4,k) * real(dzVy) * dt
                    axVz(k,i,j) = gxe(3,i) * axVz(k,i,j) + gxe(4,i) * real(dxVz) * dt
                    ayVz(k,i,j) = gye(3,j) * ayVz(k,i,j) + gye(4,j) * real(dyVz) * dt

                end do

            end do
        end do
#ifdef _OPENACC
        !$acc end kernels
#else
        !$omp end do nowait
        !$omp end parallel
        !$omp barrier
#endif 

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
