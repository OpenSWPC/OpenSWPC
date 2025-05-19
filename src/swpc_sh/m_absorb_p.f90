#include "../shared/m_debug.h"
module m_absorb_p

    !! Absorbing Boundary Condition: ADE-CFS PML based on Zhang and Shen
    !!
    !! #### PML region definition
    !!
    !! ```
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

    real(SP), allocatable :: gxc(:, :), gxe(:, :) !< damping profile along x at center/edge of voxel
    real(SP), allocatable :: gzc(:, :), gze(:, :) !< damping profile along x at center/edge of voxel

    real(SP), allocatable :: axVy(:, :), azVy(:, :)
    real(SP), allocatable :: axSxy(:, :)
    real(SP), allocatable :: azSyz(:, :)

    real(SP) :: r20x, r20z

    integer :: kbeg_min

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

        !! PML region definition
        kbeg_min = minval(kbeg_a(:))
        ! if( fullspace_mode ) kbeg_min = kbeg

        !! memory allocation
        allocate (axVy(kbeg_min:kend, ibeg:iend))
        allocate (azVy(kbeg_min:kend, ibeg:iend))
        allocate (axSxy(kbeg_min:kend, ibeg:iend))
        allocate (azSyz(kbeg_min:kend, ibeg:iend))

        axVy(kbeg_min:kend, ibeg:iend) = 0.0
        azVy(kbeg_min:kend, ibeg:iend) = 0.0
        axSxy(kbeg_min:kend, ibeg:iend) = 0.0
        azSyz(kbeg_min:kend, ibeg:iend) = 0.0

        idum = io_prm

        !$acc enter data copyin(axVy, azVy, axSxy, azSyz, gxc, gxe, gzc, gze)

    end subroutine absorb_p__setup

    subroutine absorb_p__update_vel

        !! Update velocity component in PML layer

        integer :: i, k
        real(SP) :: by
        real(MP) :: dzSyz, dxSxy

        !! Horizontal zero-derivative boundary (for plane wave mode)
        if (pw_mode) then
            if (idx == 0) then
#ifdef _OPENACC
                !$acc kernels present(Syz, Sxy)                
                !$acc loop independent
#else
                !$omp parallel private(k)
                !$omp do schedule(dynamic)
#endif
                do k = kbeg_a(1), kend
                    Syz(k, 0) = 2 * Syz(k, 1) - Syz(k, 2)
                    Sxy(k, 0) = 2 * Sxy(k, 1) - Sxy(k, 2)
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
                do k = kbeg_a(nx), kend
                    Syz(k, nx + 1) = 2 * Syz(k, nx) - Syz(k, nx - 1)
                    Sxy(k, nx + 1) = 2 * Sxy(k, nx) - Sxy(k, nx - 1)
                end do
#ifdef _OPENACC
                !$acc end kernels
#else
                !$omp end do nowait
                !$omp end parallel
#endif
            end if
        end if

        
        !! time-marching
#ifdef _OPENACC
        !$acc kernels &
        !$acc present(Vy, Sxy, Syz, axSxy, azSyz, rho, gxc, gzc, kbeg_a)
        !$acc loop independent
#else
        !$omp parallel &
        !$omp private( dzSyz, dxSxy, by, i, k )
        !$omp do schedule(dynamic)
#endif
        do i = ibeg, iend

#ifdef _OPENACC
            !$acc loop vector independent
#endif
            do k = kbeg_a(i), kend

                dzSyz = (Syz(k,i) - Syz(k-1,i)) * r20z
                dxSxy = (Sxy(k,i) - Sxy(k,i-1)) * r20x

                by = 1.0 / rho(k,i)

                Vy(k,i) = Vy(k,i) &
                        + by * (gxc(1,i) * dxSxy      + gzc(1,i) * dzSyz &
                              + gxc(2,i) * axSxy(k,i) + gzc(2,i) * azSyz(k,i)) * dt

                axSxy(k,i) = gxc(3,i) * axSxy(k,i) + gxc(4,i) * dxSxy * dt
                azSyz(k,i) = gzc(3,k) * azSyz(k,i) + gzc(4,k) * dzSyz * dt

            end do
        end do
#ifdef _OPENACC
        !$acc end kernels
#else
        !$omp end do nowait
        !$omp end parallel
#endif

        ! if( fullspace_mode ) then
        !   ! $omp parallel &
        !   ! $omp private( gxc0, gzc0, dzSyz, dxSxy, by, i, k )
        !   ! $omp do schedule(dynamic)
        !   do i=ibeg_k,iend_k

        !     gxc0(1:4) = gxc(1:4,i)

        !     !!
        !     !! Derivatives
        !     !!
        !     do k=1,na

        !       dzSyz(k) = (  Syz(k  ,i ) - Syz(k-1,i   )  ) * r20z
        !       dxSxy(k) = (  Sxy(k  ,i ) - Sxy(k  ,i-1 )  ) * r20x

        !     end do

        !     !!
        !     !! update velocity
        !     !!
        !     do k=1, na

        !       gzc0(1:4) = gzc(1:4,k)

        !       !!
        !       !! Velocity Updates
        !       !!
        !       by = 1.0 /  rho(k,i)

        !       Vy(k,i) = Vy(k,i) &
        !           + by * ( gxc0(1) * dxSxy(k)   + gzc0(1) * dzSyz(k)       &
        !           + gxc0(2) * axSxy(k,i) + gzc0(2) * azSyz(k,i)  ) * dt

        !       !!
        !       !! ADE updates
        !       !!
        !       axSxy(k,i) = gxc0(3) * axSxy(k,i) + gxc0(4) * dxSxy(k) * dt
        !       azSyz(k,i) = gzc0(3) * azSyz(k,i) + gzc0(4) * dzSyz(k) * dt

        !     end do
        !   end do
        !   ! $omp end do nowait
        !   ! $omp end parallel
        ! end if

        ! ! $omp barrier

    end subroutine absorb_p__update_vel

    subroutine absorb_p__update_stress

        integer :: i, k
        real(SP) :: nnn, pnn, npn
        real(SP) :: muxy, muyz
        real(SP) :: epsl = epsilon(1.0)
        real(MP) :: dxVy, dzVy

        !! Horizontal zero-derivative boundary (for plane wave mode)
        if (pw_mode) then
            if (idx == 0) then
#ifdef _OPENACC
                !$acc kernels present(Vy)
                !$acc loop independent
#else
                !$omp parallel private(k)
                !$omp do schedule(dynamic)
#endif
                do k = kbeg_a(1), kend
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
                do k = kbeg_a(nx), kend
                    Vy(k,nx+1) = 2 * Vy(k,nx) - Vy(k,nx-1)
                end do
#ifdef _OPENACC
                !$acc end kernels
#else
                !$omp end do nowait
                !$omp end parallel
#endif
            end if
        end if

        !! Time-marching
#ifdef _OPENACC
        !$acc kernels &
        !$acc present(Vy, Sxy, Syz, gxc, gxe, gzc, gze, axVy, azVy, mu, kbeg_a)
        !$acc loop independent
#else
        !$omp parallel &
        !$omp private(i, k, dxVy, dzVy, nnn, pnn, npn, muxy, muyz )
        !$omp do schedule(dynamic)
#endif
        do i = ibeg, iend

            !$acc loop vector independent
            do k = kbeg_a(i), kend

                dxVy = (Vy(k,i+1) - Vy(k,i)) * r20x
                dzVy = (Vy(k+1,i) - Vy(k,i)) * r20z

                nnn = mu(k  ,i  )
                pnn = mu(k+1,i  )
                npn = mu(k  ,i+1)
                muxy = 2 * nnn * npn / (nnn + npn + epsl)
                muyz = 2 * nnn * pnn / (nnn + pnn + epsl)

                axVy(k,i) = gxe(3,i) * axVy(k,i) + gxe(4,i) * dxVy * dt
                azVy(k,i) = gze(3,k) * azVy(k,i) + gze(4,k) * dzVy * dt

                Syz(k,i) = Syz(k,i) + muyz * (gze(1,k) * dzVy + gze(2,k) * azVy(k,i)) * dt
                Sxy(k,i) = Sxy(k,i) + muxy * (gxe(1,i) * dxVy + gxe(2,i) * axVy(k,i)) * dt

            end do
        end do
#ifdef _OPENACC
        !$acc end kernels
#else
        !$omp end do nowait
        !$omp end parallel
#endif

        ! if( fullspace_mode ) then
        !   ! $omp parallel &
        !   ! $omp private(i, k, gxe0, gze0, dxVy, dzVy, nnn, pnn, npn, muxy, muyz )
        !   ! $omp do schedule(dynamic)
        !   do i=ibeg_k,iend_k

        !     gxe0(1:4) = gxe(1:4,i)

        !     !!
        !     !! Derivatives
        !     !!
        !     do k=kbeg, na

        !       dxVy(k) = (  Vy(k  ,i+1) - Vy(k,i)  ) * r20x
        !       dzVy(k) = (  Vy(k+1,i  ) - Vy(k,i)  ) * r20z

        !     end do

        !     !!
        !     !! Update Shear Stress
        !     !!
        !     do k=kbeg, na

        !       gze0(1:4) = gze(1:4,k)

        !       !!
        !       !! effective rigidity for shear stress components
        !       !!

        !       nnn = mu (k  ,i )
        !       pnn = mu (k+1,i )
        !       npn = mu (k,  i+1)
        !       muxy = 2*nnn*npn / ( nnn + npn + epsl )
        !       muyz = 2*nnn*pnn / ( nnn + pnn + epsl )

        !       axVy(k,i) = gxe0(3) * axVy(k,i) + gxe0(4) * dxVy(k) * dt
        !       azVy(k,i) = gze0(3) * azVy(k,i) + gze0(4) * dzVy(k) * dt
        !       Syz(k,i) = Syz(k,i) + muyz* ( gze0(1) * dzVy(k) + gze0(2) * azVy(k,i) ) * dt
        !       Sxy(k,i) = Sxy(k,i) + muxy* ( gxe0(1) * dxVy(k) + gxe0(2) * axVy(k,i) ) * dt

        !     end do
        !   end do
        !   ! $omp end do nowait
        !   ! $omp end parallel

        ! end if

        ! ! $omp barrier

    end subroutine absorb_p__update_stress

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
