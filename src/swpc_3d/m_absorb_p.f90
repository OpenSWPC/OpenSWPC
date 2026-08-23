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

    real(SP), allocatable :: gxc(:,:), gxe(:,:) !! damping profile along x at center/edge of voxel
    real(SP), allocatable :: gyc(:,:), gye(:,:) !! damping profile along y at center/edge of voxel
    real(SP), allocatable :: gzc(:,:), gze(:,:) !! damping profile along z at center/edge of voxel

    real(SP), allocatable :: axVx(:), ayVx(:), azVx(:)
    real(SP), allocatable :: axVy(:), ayVy(:), azVy(:)
    real(SP), allocatable :: axVz(:), ayVz(:), azVz(:)
    real(SP), allocatable :: axSxx(:), aySxy(:), azSxz(:)
    real(SP), allocatable :: axSxy(:), aySyy(:), azSyz(:)
    real(SP), allocatable :: axSxz(:), aySyz(:), azSzz(:)

    real(MP) :: r20x, r20y, r20z
    real(MP) :: rc40x, rc41x, rc40y, rc41y, rc40z, rc41z
    real(MP) :: rd40x, rd41x, rd40y, rd41y, rd40z, rd41z

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
        rc40x = 17.0_MP / 16.0_MP / dx
        rc40y = 17.0_MP / 16.0_MP / dy
        rc40z = 17.0_MP / 16.0_MP / dz
        rc41x =  1.0_MP / 48.0_MP / dx
        rc41y =  1.0_MP / 48.0_MP / dy
        rc41z =  1.0_MP / 48.0_MP / dz
        rd40x = -1.0_MP / 16.0_MP / dx
        rd40y = -1.0_MP / 16.0_MP / dy
        rd40z = -1.0_MP / 16.0_MP / dz
        rd41x = -1.0_MP / 48.0_MP / dx
        rd41y = -1.0_MP / 48.0_MP / dy
        rd41z = -1.0_MP / 48.0_MP / dz       

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

        !! memory allocation
        allocate (axVx (n_sponge_cell), source=0.0)
        allocate (ayVx (n_sponge_cell), source=0.0)
        allocate (azVx (n_sponge_cell), source=0.0)
        allocate (axVy (n_sponge_cell), source=0.0)
        allocate (ayVy (n_sponge_cell), source=0.0)
        allocate (azVy (n_sponge_cell), source=0.0)
        allocate (axVz (n_sponge_cell), source=0.0)
        allocate (ayVz (n_sponge_cell), source=0.0)
        allocate (azVz (n_sponge_cell), source=0.0)
        allocate (axSxx(n_sponge_cell), source=0.0)
        allocate (aySxy(n_sponge_cell), source=0.0)
        allocate (azSxz(n_sponge_cell), source=0.0)
        allocate (axSxy(n_sponge_cell), source=0.0)
        allocate (aySyy(n_sponge_cell), source=0.0)
        allocate (azSyz(n_sponge_cell), source=0.0)
        allocate (axSxz(n_sponge_cell), source=0.0)
        allocate (aySyz(n_sponge_cell), source=0.0)
        allocate (azSzz(n_sponge_cell), source=0.0)

        idum = io_prm

        !$acc enter data copyin(axVx, ayVx, azVx, axVy, ayVy, azVy, axVz, ayVz, azVz, &
        !$acc                   axSxx, aySxy, azSxz, axSxy, aySyy, azSyz, axSxz, aySyz, azSzz, &
        !$acc                   gxc, gxe, gyc, gye, gzc, gze)

    end subroutine absorb_p__setup


    subroutine absorb_p__update_vel

        integer :: ibox

        if (pw_mode) call set_stress_boundary()

        !! time-marching
        do ibox = 1, 6
            if(box(ibox)%ncell > 0) call update_vel_core(box(ibox))
        end do

    end subroutine absorb_p__update_vel


    subroutine set_stress_boundary()
        !! Horizontal zero-derivative boundary (for plane wave mode)

        integer :: i, j, k

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

    end subroutine set_stress_boundary


    subroutine update_vel_core(bb) 

        type(t_box), intent(in) :: bb
        integer :: p, p0
        integer :: i, j, k
        real(MP) :: dxSxx, dySxy, dzSxz
        real(MP) :: dxSxy, dySyy, dzSyz
        real(MP) :: dxSxz, dySyz, dzSzz
        real(MP) :: d3Sx3, d3Sy3, d3Sz3
        integer :: isign
        real(MP) :: re40x, re41x, re40y, re41y, re40z, re41z

#ifdef _OPENACC
        !$acc kernels &
        !$acc present(Vx, Vy, Vz, Sxx, Syy, Szz, Syz, Sxz, Sxy, &
        !$acc         axSxx, aySxy, azSxz, axSxy, aySyy, azSyz, axSxz, aySyz, azSzz, &
        !$acc         gxc, gxe, gyc, gye, gzc, gze, bb)
        !$acc loop independent collapse(3)
#else
        !$omp parallel &
        !$omp private( dxSxx, dySyy, dzSzz, dySyz, dzSyz, dxSxz, dzSxz, dxSxy ,dySxy ) &
        !$omp private(re40x, re41x, re40y, re41y, re40z, re41z, isign) &
        !$omp private( d3Sx3, d3Sy3, d3Sz3 ) &
        !$omp private( i, j, k, p )
        !$omp do &
        !$omp schedule(dynamic)
#endif
        do j = bb%jb, bb%je
            do i = bb%ib, bb%ie
                do k = bb%kb, bb%ke

                    p = bb%offset + (j - bb%jb) * bb%nx * bb%nz + (i-bb%ib) * bb%nz + (k - bb%kb + 1)

                    isign = sign(1, max((k - kfs_top(i,j)) * (kfs_bot(i,j) - k), &
                                        (k - kob_top(i,j)) * (kob_bot(i,j) - k)))

                    re40x = rc40x + isign * rd40x
                    re41x = rc41x + isign * rd41x
                    re40y = rc40y + isign * rd40y
                    re41y = rc41y + isign * rd41y
                    re40z = rc40z + isign * rd40z
                    re41z = rc41z + isign * rd41z

                    dxSxx = (Sxx(k  ,i+1,j  ) - Sxx(k  ,i  ,j  )) * re40x - (Sxx(k  ,i+2,j  ) - Sxx(k  ,i-1,j  )) * re41x
                    dySxy = (Sxy(k  ,i  ,j  ) - Sxy(k  ,i  ,j-1)) * re40y - (Sxy(k  ,i  ,j+1) - Sxy(k  ,i  ,j-2)) * re41y 
                    dzSxz = (Sxz(k  ,i  ,j  ) - Sxz(k-1,i  ,j  )) * re40z - (Sxz(k+1,i  ,j  ) - Sxz(k-2,i  ,j  )) * re41z
                    dxSxy = (Sxy(k  ,i  ,j  ) - Sxy(k  ,i-1,j  )) * re40x - (Sxy(k  ,i+1,j  ) - Sxy(k  ,i-2,j  )) * re41x
                    dySyy = (Syy(k  ,i  ,j+1) - Syy(k  ,i  ,j  )) * re40y - (Syy(k  ,i  ,j+2) - Syy(k  ,i  ,j-1)) * re41y
                    dzSyz = (Syz(k  ,i  ,j  ) - Syz(k-1,i  ,j  )) * re40z - (Syz(k+1,i  ,j  ) - Syz(k-2,i  ,j  )) * re41z
                    dxSxz = (Sxz(k  ,i  ,j  ) - Sxz(k  ,i-1,j  )) * re40x - (Sxz(k  ,i+1,j  ) - Sxz(k  ,i-2,j  )) * re41x
                    dySyz = (Syz(k  ,i  ,j  ) - Syz(k  ,i  ,j-1)) * re40y - (Syz(k  ,i  ,j+1) - Syz(k  ,i  ,j-2)) * re41y
                    dzSzz = (Szz(k+1,i  ,j  ) - Szz(k  ,i  ,j  )) * re40z - (Szz(k+2,i  ,j  ) - Szz(k-1,i  ,j  )) * re41z

                    d3Sx3 = gxe(1,i) * dxSxx + gxe(2,i) * axSxx(p) &
                          + gyc(1,j) * dySxy + gyc(2,j) * aySxy(p) &
                          + gzc(1,k) * dzSxz + gzc(2,k) * azSxz(p)
                    d3Sy3 = gxc(1,i) * dxSxy + gxc(2,i) * axSxy(p) &
                          + gye(1,j) * dySyy + gye(2,j) * aySyy(p) &
                          + gzc(1,k) * dzSyz + gzc(2,k) * azSyz(p)
                    d3Sz3 = gxc(1,i) * dxSxz + gxc(2,i) * axSxz(p) &
                          + gyc(1,j) * dySyz + gyc(2,j) * aySyz(p) &
                          + gze(1,k) * dzSzz + gze(2,k) * azSzz(p)

                    !! Velocity Updates
                    Vx(k,i,j) = Vx(k,i,j) + 2.0 / (rho(k,i,j) + rho(k,i+1,j)) * d3Sx3 * dt
                    Vy(k,i,j) = Vy(k,i,j) + 2.0 / (rho(k,i,j) + rho(k,i,j+1)) * d3Sy3 * dt
                    Vz(k,i,j) = Vz(k,i,j) + 2.0 / (rho(k,i,j) + rho(k+1,i,j)) * d3Sz3 * dt

                    !! ADE updates
                    axSxx(p) = gxe(3,i) * axSxx(p) + gxe(4,i) * real(dxSxx) * dt
                    aySxy(p) = gyc(3,j) * aySxy(p) + gyc(4,j) * real(dySxy) * dt
                    azSxz(p) = gzc(3,k) * azSxz(p) + gzc(4,k) * real(dzSxz) * dt
                    axSxy(p) = gxc(3,i) * axSxy(p) + gxc(4,i) * real(dxSxy) * dt
                    aySyy(p) = gye(3,j) * aySyy(p) + gye(4,j) * real(dySyy) * dt
                    azSyz(p) = gzc(3,k) * azSyz(p) + gzc(4,k) * real(dzSyz) * dt
                    axSxz(p) = gxc(3,i) * axSxz(p) + gxc(4,i) * real(dxSxz) * dt
                    aySyz(p) = gyc(3,j) * aySyz(p) + gyc(4,j) * real(dySyz) * dt
                    azSzz(p) = gze(3,k) * azSzz(p) + gze(4,k) * real(dzSzz) * dt

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

    end subroutine update_vel_core


    subroutine absorb_p__update_stress

        integer :: ibox

        !! Horizontal zero-derivative boundary (for plane wave mode)
        if (pw_mode) call set_vel_boundary()

        !! Time-marching
        do ibox = 1, 6
            if(box(ibox)%ncell > 0) call update_stress_core(box(ibox))
        end do

    end subroutine absorb_p__update_stress

    subroutine set_vel_boundary()

        integer :: i, j, k

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

    end subroutine set_vel_boundary


    subroutine update_stress_core(bb)

        type(t_box), intent(in) :: bb
        integer :: i, j, k, p0, p, m, isign
        real(SP) :: mu2, lam2mu
        real(SP) :: taup1, taus1, taup_plus1, taus_plus1
        real(SP) :: d3v3, dxVx_dyVy, dxVx_dzVz, dyVy_dzVz
        real(SP) :: dxVx_ade, dyVy_ade, dzVz_ade
        real(SP) :: Rxx_n, Ryy_n, Rzz_n, Ryz_n, Rxz_n, Rxy_n
        real(SP) :: gxc0(4), gxe0(4), gyc0(4), gye0(4), gzc0(4), gze0(4)
        real(SP) :: dxVy_dyVx, dxVz_dzVx, dyVz_dzVy
        real(MP) :: dxVx, dyVx, dzVx
        real(MP) :: dxVy, dyVy, dzVy
        real(MP) :: dxVz, dyVz, dzVz
        real(SP) :: muxz, muyz, muxy
        real(MP) :: re40x, re41x, re40y, re41y, re40z, re41z
        real(SP) :: epsl = epsilon(1.0)

#ifdef _OPENACC
        !$acc kernels &
        !$acc present(Vx, Vy, Vz, Sxx, Syy, Szz, Rxx, Ryy, Rzz, &
        !$acc         mu, lam, taup, taus, c1, c2, d1, d2, &
        !$acc         axVx, ayVy, azVz, &
        !$acc         kfs_top, kfs_bot, kob_top, kob_bot, &
        !$acc         gxc, gxe, gyc, gye, gzc, gze, bb)
        !$acc loop independent collapse(3)
#else
        !$omp parallel &
        !$omp private( dxVx, dyVy, dzVz ) &
        !$omp private( mu2, lam2mu ) &
        !$omp private( taup1, taus1, taup_plus1, taus_plus1 ) &
        !$omp private( d3v3, dyVy_dzVz, dxVx_dzVz, dxVx_dyVy ) &
        !$omp private( Rxx_n, Ryy_n, Rzz_n ) &
        !$omp private( re40x, re41x, re40y, re41y, re40z, re41z, isign) &
        !$omp private( i, j, k, m, p )
        !$omp private( dxVx_ade, dyVy_ade, dzVz_ade ) &
        !$omp do &
        !$omp schedule(dynamic)
#endif
        do j = bb%jb, bb%je
            do i = bb%ib, bb%ie
                !ocl unroll('full')
                !ocl swp
                !OCL SWP_IREG_RATE(200)
                do k = bb%kb, bb%ke

                    p = bb%offset + (j - bb%jb) * bb%nx * bb%nz + (i-bb%ib) * bb%nz + (k - bb%kb + 1)

                    isign = sign(1, max((k - kfs_top(i,j)) * (kfs_bot(i,j) - k), &
                                        (k - kob_top(i,j)) * (kob_bot(i,j) - k)))

                    re40x = rc40x + isign * rd40x
                    re41x = rc41x + isign * rd41x
                    re40y = rc40y + isign * rd40y
                    re41y = rc41y + isign * rd41y
                    re40z = rc40z + isign * rd40z
                    re41z = rc41z + isign * rd41z

                    dxVx = (Vx(k  ,i  ,j  ) - Vx(k  ,i-1,j  )) * re40x - (Vx(k  ,i+1,j  ) - Vx(k  ,i-2,j  )) * re41x
                    dyVy = (Vy(k  ,i  ,j  ) - Vy(k  ,i  ,j-1)) * re40y - (Vy(k  ,i  ,j+1) - Vy(k  ,i  ,j-2)) * re41y
                    dzVz = (Vz(k  ,i  ,j  ) - Vz(k-1,i  ,j  )) * re40z - (Vz(k+1,i  ,j  ) - Vz(k-2,i  ,j  )) * re41z
    
                    mu2 = 2 * mu(k,i,j)
                    lam2mu = lam(k,i,j) + mu2
                    taup1 = taup(k,i,j)
                    taus1 = taus(k,i,j)

                    dxVx_ade = gxc(1,i) * dxVx + gxc(2,i) * axVx(p)
                    dyVy_ade = gyc(1,j) * dyVy + gyc(2,j) * ayVy(p)
                    dzVz_ade = gzc(1,k) * dzVz + gzc(2,k) * azVz(p)

                    d3v3      = dxVx_ade + dyVy_ade + dzVz_ade
                    dyVy_dzVz = dyVy_ade + dzVz_ade
                    dxVx_dzVz = dxVx_ade + dzVz_ade
                    dxVx_dyVy = dxVx_ade + dyVy_ade

                    Rxx_n = 0.0
                    Ryy_n = 0.0
                    Rzz_n = 0.0

                    !$acc loop seq reduction(+:Rxx_n,Ryy_n,Rzz_n)
                    do m=1, nm
                        Rxx(m,k,i,j) = c1(m) * Rxx(m,k,i,j) - c2(m) * (lam2mu * taup1 * d3v3 - mu2 * taus1 * dyVy_dzVz) * dt
                        Ryy(m,k,i,j) = c1(m) * Ryy(m,k,i,j) - c2(m) * (lam2mu * taup1 * d3v3 - mu2 * taus1 * dxVx_dzVz) * dt
                        Rzz(m,k,i,j) = c1(m) * Rzz(m,k,i,j) - c2(m) * (lam2mu * taup1 * d3v3 - mu2 * taus1 * dxVx_dyVy) * dt
                        Rxx_n = Rxx_n + d1(m) * Rxx(m,k,i,j)
                        Ryy_n = Ryy_n + d1(m) * Ryy(m,k,i,j)
                        Rzz_n = Rzz_n + d1(m) * Rzz(m,k,i,j)
                    end do

                    taup_plus1 = 1 + taup1 * ( 1 + d2 )
                    taus_plus1 = 1 + taus1 * ( 1 + d2 )
                    
                    Sxx (k,i,j) = Sxx (k,i,j) + (lam2mu * taup_plus1 * d3v3 - mu2 * taus_plus1 * dyVy_dzVz + Rxx_n) * dt
                    Syy (k,i,j) = Syy (k,i,j) + (lam2mu * taup_plus1 * d3v3 - mu2 * taus_plus1 * dxVx_dzVz + Ryy_n) * dt
                    Szz (k,i,j) = Szz (k,i,j) + (lam2mu * taup_plus1 * d3v3 - mu2 * taus_plus1 * dxVx_dyVy + Rzz_n) * dt

                    axVx(p) = gxc(3,i) * axVx(p) + gxc(4,i) * dxVx * dt
                    ayVy(p) = gyc(3,j) * ayVy(p) + gyc(4,j) * dyVy * dt
                    azVz(p) = gzc(3,k) * azVz(p) + gzc(4,k) * dzVz * dt

                end do
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
        !$acc present(Vx, Vy, Vz, Sxy, Sxz, Syz, Rxy, Rxz, Ryz, mu, &
        !$acc         c1, c2, d1, d2, &
        !$acc         axVx, ayVx, azVx, axVy, ayVy, azVy, axVz, ayVz, azVz, &
        !$acc         kfs_top, kfs_bot, kob_top, kob_bot, &
        !$acc         gxc, gxe, gyc, gye, gzc, gze, bb)
        !$acc loop independent collapse(3)
#else
        !$omp parallel &
        !$omp private( dxVy_dyVx, dxVz_dzVx, dyVz_dzVy, muyz, muxz, muxy) &
        !$omp private( dxVy, dxVz, dyVx, dyVz, dzVx, dzVy) &
        !$omp private( taus1, taus_plus1 ) &
        !$omp private( Ryz_n, Rxz_n, Rxy_n ) &
        !$omp private( re40x, re41x, re40y, re41y, re40z, re41z, isign) &
        !$omp private( i, j, k, m, p ) &
        !$omp do &
        !$omp schedule(dynamic)
#endif
        do j = bb%jb, bb%je
            do i = bb%ib, bb%ie
                !ocl unroll('full')
                !ocl swp
                !OCL SWP_IREG_RATE(300)
                !OCL SWP_FREG_RATE(200)
                do k = bb%kb, bb%ke

                    p = bb%offset + (j - bb%jb) * bb%nx * bb%nz + (i-bb%ib) * bb%nz + (k - bb%kb + 1)

                    isign = sign(1, max((k - kfs_top(i,j)) * (kfs_bot(i,j) - k), &
                                        (k - kob_top(i,j)) * (kob_bot(i,j) - k)))

                    re40x = rc40x + isign * rd40x
                    re41x = rc41x + isign * rd41x
                    re40y = rc40y + isign * rd40y
                    re41y = rc41y + isign * rd41y
                    re40z = rc40z + isign * rd40z
                    re41z = rc41z + isign * rd41z

                    dxVy = (Vy(k  ,i+1,j  ) - Vy(k  ,i  ,j  )) * re40x - (Vy(k  ,i+2,j  ) - Vy(k  ,i-1,j  )) * re41x
                    dyVx = (Vx(k  ,i  ,j+1) - Vx(k  ,i  ,j  )) * re40y - (Vx(k  ,i  ,j+2) - Vx(k  ,i  ,j-1)) * re41y
                    dxVz = (Vz(k  ,i+1,j  ) - Vz(k  ,i  ,j  )) * re40x - (Vz(k  ,i+2,j  ) - Vz(k  ,i-1,j  )) * re41x
                    dzVx = (Vx(k+1,i  ,j  ) - Vx(k  ,i  ,j  )) * re40z - (Vx(k+2,i  ,j  ) - Vx(k-1,i  ,j  )) * re41z
                    dyVz = (Vz(k  ,i  ,j+1) - Vz(k  ,i  ,j  )) * re40y - (Vz(k  ,i  ,j+2) - Vz(k  ,i  ,j-1)) * re41y
                    dzVy = (Vy(k+1,i  ,j  ) - Vy(k  ,i  ,j  )) * re40z - (Vy(k+2,i  ,j  ) - Vy(k-1,i  ,j  )) * re41z

                    dxVy_dyVx = gxe(1,i) * dxVy + gxe(2,i) * axVy(p) &
                              + gye(1,j) * dyVx + gye(2,j) * ayVx(p)
                    dxVz_dzVx = gxe(1,i) * dxVz + gxe(2,i) * axVz(p) &
                              + gze(1,k) * dzVx + gze(2,k) * azVx(p)
                    dyVz_dzVy = gye(1,j) * dyVz + gye(2,j) * ayVz(p) &
                              + gze(1,k) * dzVy + gze(2,k) * azVy(p)

                    muxz = 4 * mu(k  ,i  ,j  ) * mu(k+1,i  ,j  ) * mu(k  ,i+1,j  ) * mu(k+1,i+1,j  ) &
                           / ( mu(k  ,i  ,j  ) * mu(k+1,i  ,j  ) * mu(k  ,i+1,j  ) &
                             + mu(k  ,i  ,j  ) * mu(k+1,i  ,j  ) * mu(k+1,i+1,j  ) &
                             + mu(k  ,i  ,j  ) * mu(k  ,i+1,j  ) * mu(k+1,i+1,j  ) &
                             + mu(k+1,i  ,j  ) * mu(k  ,i+1,j  ) * mu(k+1,i+1,j  ) + epsl)

                    muxy = 4 * mu(k  ,i  ,j  ) * mu(k  ,i+1,j  ) * mu(k  ,i  ,j+1) * mu(k  ,i+1,j+1) &
                           / ( mu(k  ,i  ,j  ) * mu(k  ,i+1,j  ) * mu(k  ,i  ,j+1) &
                             + mu(k  ,i  ,j  ) * mu(k  ,i+1,j  ) * mu(k  ,i+1,j+1) &
                             + mu(k  ,i  ,j  ) * mu(k  ,i  ,j+1) * mu(k  ,i+1,j+1) &
                             + mu(k  ,i+1,j  ) * mu(k  ,i  ,j+1) * mu(k  ,i+1,j+1) + epsl)

                    muyz = 4 * mu(k  ,i  ,j  ) * mu(k+1,i  ,j  ) * mu(k  ,i  ,j+1) * mu(k+1,i  ,j+1) &
                           / ( mu(k  ,i  ,j  ) * mu(k+1,i  ,j  ) * mu(k  ,i  ,j+1) &
                             + mu(k  ,i  ,j  ) * mu(k+1,i  ,j  ) * mu(k+1,i  ,j+1) &
                             + mu(k  ,i  ,j  ) * mu(k  ,i  ,j+1) * mu(k+1,i  ,j+1) &
                             + mu(k+1,i  ,j  ) * mu(k  ,i  ,j+1) * mu(k+1,i  ,j+1) + epsl)

                    !! medium copy
                    taus1 = taus(k, i, j)

                    !! update memory variables
                    Ryz_n = 0.0
                    Rxz_n = 0.0
                    Rxy_n = 0.0
                    !$acc loop seq reduction(+:Rxy_n,Ryz_n,Rxz_n)
                    do m = 1, nm
                        Ryz(m,k,i,j) = c1(m) * Ryz(m,k,i,j) - c2(m) * muyz * taus1 * dyVz_dzVy * dt
                        Rxz(m,k,i,j) = c1(m) * Rxz(m,k,i,j) - c2(m) * muxz * taus1 * dxVz_dzVx * dt
                        Rxy(m,k,i,j) = c1(m) * Rxy(m,k,i,j) - c2(m) * muxy * taus1 * dxVy_dyVx * dt
                        Ryz_n = Ryz_n + d1(m) * Ryz(m,k,i,j)
                        Rxz_n = Rxz_n + d1(m) * Rxz(m,k,i,j)
                        Rxy_n = Rxy_n + d1(m) * Rxy(m,k,i,j)
                    end do

                    !! update stress components
                    taus_plus1 = 1 + taus1 * (1 + d2)

                    Syz(k,i,j) = Syz(k,i,j) + (muyz * taus_plus1 * dyVz_dzVy + Ryz_n) * dt
                    Sxz(k,i,j) = Sxz(k,i,j) + (muxz * taus_plus1 * dxVz_dzVx + Rxz_n) * dt
                    Sxy(k,i,j) = Sxy(k,i,j) + (muxy * taus_plus1 * dxVy_dyVx + Rxy_n) * dt

                    ayVx(p) = gye(3,j) * ayVx(p) + gye(4,j) * dyVx * dt
                    azVx(p) = gze(3,k) * azVx(p) + gze(4,k) * dzVx * dt
                    axVy(p) = gxe(3,i) * axVy(p) + gxe(4,i) * dxVy * dt
                    azVy(p) = gze(3,k) * azVy(p) + gze(4,k) * dzVy * dt
                    axVz(p) = gxe(3,i) * axVz(p) + gxe(4,i) * dxVz * dt
                    ayVz(p) = gye(3,j) * ayVz(p) + gye(4,j) * dyVz * dt

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

    end subroutine update_stress_core


    subroutine damping_profile(x, H, xbeg0, xend0, g)

        !! ADE-CFS PML damping factor according to Zhao and Shen

        real(SP), intent(in) :: x   !! cartesian coordinate location
        real(SP), intent(in) :: H   !! absorption layer thickness
        real(SP), intent(in) :: xbeg0
        real(SP), intent(in) :: xend0
        real(SP), intent(out) :: g(4) !! damping prof

        real(SP) :: R0 !! reflection coefficient
        real(SP) :: d0, a0, b0
        integer, parameter :: pd = 2
        integer, parameter :: pa = 1
        integer, parameter :: pb = 3
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
