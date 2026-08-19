#include "../shared/m_debug.h"
module m_absorb_p

    !! Absorbing Boundary Condition: ADE-CFS PML based on Zhang and Shen
    !!
    !! #### PML region definition
    !!
    !! ```
    !!  <X-Z cross section>
    !!  +-----+--------------------------+-----+ 
    !!  |     |            (3)           |     | 
    !!  |     +--------------------------+     | 
    !!  |     |                          |     | 
    !!  |     |                          |     | 
    !!  |     |                          |     | 
    !!  | (1) |     interior region      | (2) | 
    !!  |     |       (m_kernel)         |     | 
    !!  |     |                          |     | 
    !!  |     |                          |     | 
    !!  |     +--------------------------+     | 
    !!  |     |            (4)           |     | 
    !!  +-----+--------------------------+-----+ 
    !! ```
    !! (thickness of (3) in k-direction is zero for fullspace_mode = .false. )
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
        if (pw_mode) call absorb_p__set_stress_boundary()

        do ibox=1, 4
            if(box(ibox)%ncell > 0) call absorb_p__update_vel_core(box(ibox))
        end do

    end subroutine absorb_p__update_vel

    subroutine absorb_p__set_stress_boundary()

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

    end subroutine absorb_p__set_stress_boundary

    subroutine absorb_p__update_vel_core(bb)

        type(t_box) :: bb
        integer :: p
        integer :: i, k
        real(SP) :: by
        real(MP) :: dzSyz, dxSxy

#ifdef _OPENACC
        !$acc kernels &
        !$acc present(Vy, Sxy, Syz, axSxy, azSyz, rho, gxc, gzc, bb)
        !$acc loop independent
#else
        !$omp parallel &
        !$omp private( dzSyz, dxSxy, by, i, k, p )
        !$omp do schedule(dynamic)
#endif
        do i = bb%ib, bb%ie

#ifdef _OPENACC
            !$acc loop vector independent
#endif
            do k = bb%kb, bb%ke

                p = bb%offset + (i-bb%ib) * bb%nz + (k - bb%kb + 1)

                dzSyz = (Syz(k,i) - Syz(k-1,i)) * r20z
                dxSxy = (Sxy(k,i) - Sxy(k,i-1)) * r20x

                by = 1.0 / rho(k,i)

                Vy(k,i) = Vy(k,i) &
                        + by * (gxc(1,i) * dxSxy    + gzc(1,i) * dzSyz &
                              + gxc(2,i) * axSxy(p) + gzc(2,i) * azSyz(p)) * dt

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

    end subroutine absorb_p__update_vel_core

    subroutine absorb_p__update_stress

        integer :: ibox

        !! Horizontal zero-derivative boundary (for plane wave mode)
        if (pw_mode) call absorb_p__set_vel_boundary()

        do ibox=1, 4
            if( box(ibox)%ncell > 0 ) call absorb_p__update_stress_core(box(ibox))
        end do


    end subroutine absorb_p__update_stress

    subroutine absorb_p__set_vel_boundary()

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
    end subroutine absorb_p__set_vel_boundary


    subroutine absorb_p__update_stress_core(bb)

        type(t_box), intent(in) :: bb
        integer :: i, k
        real(SP) :: nnn, pnn, npn
        real(SP) :: muxy, muyz
        real(SP) :: epsl = epsilon(1.0)
        real(MP) :: dxVy, dzVy
        integer :: p

#ifdef _OPENACC
        !$acc kernels &
        !$acc present(Vy, Sxy, Syz, gxc, gxe, gzc, gze, axVy, azVy, mu, bb)
        !$acc loop independent
#else
        !$omp parallel &
        !$omp private(i, k, dxVy, dzVy, nnn, pnn, npn, muxy, muyz, p )
        !$omp do schedule(dynamic)
#endif
        do i = bb%ib, bb%ie

            !$acc loop vector independent
            do k = bb%kb, bb%ke

                p = bb%offset + (i-bb%ib) * bb%nz + (k - bb%kb + 1)

                dxVy = (Vy(k,i+1) - Vy(k,i)) * r20x
                dzVy = (Vy(k+1,i) - Vy(k,i)) * r20z

                nnn = mu(k  ,i  )
                pnn = mu(k+1,i  )
                npn = mu(k  ,i+1)
                muxy = 2 * nnn * npn / (nnn + npn + epsl)
                muyz = 2 * nnn * pnn / (nnn + pnn + epsl)

                axVy(p) = gxe(3,i) * axVy(p) + gxe(4,i) * dxVy * dt
                azVy(p) = gze(3,k) * azVy(p) + gze(4,k) * dzVy * dt

                Syz(k,i) = Syz(k,i) + muyz * (gze(1,k) * dzVy + gze(2,k) * azVy(p)) * dt
                Sxy(k,i) = Sxy(k,i) + muxy * (gxe(1,i) * dxVy + gxe(2,i) * axVy(p)) * dt

            end do
        end do
#ifdef _OPENACC
        !$acc end kernels
#else
        !$omp end do nowait
        !$omp end parallel
#endif

        
    end subroutine absorb_p__update_stress_core


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
        integer, parameter :: pa = 2
        integer, parameter :: pb = 2
        real(SP), parameter :: cp = 6.0 !! assumed P-wave velocity
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
