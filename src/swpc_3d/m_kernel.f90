#include "../shared/m_debug.h"
module m_kernel

    !! Computation kernel for FDM numerical simulation
    !!
    !! Copyright 2013-2025 Takuto Maeda. All rights reserved. This project is released under the MIT license.

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
    real(MP) :: rc40x, rc41x, rc40y, rc41y, rc40z, rc41z
    real(MP) :: rd40x, rd41x, rd40y, rd41y, rd40z, rd41z

contains

    subroutine kernel__setup()

        !! Setup
        !!
        !! this routine MUST BE CALLED AFTER medium__setup since it uses viscoelastic function
        !! use ts(:) array, dx, dy, dz from global module

        integer :: m

        call pwatch__on("kernel__setup")
        call assert(medium__initialized())


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

#ifdef _OPENACC
        !$acc enter data copyin(c1, c2, d1, d2)
#endif

    end subroutine kernel__setup

    subroutine kernel__update_vel()

        !! Update vel for one time step

        integer :: i, j, k
        real(MP) :: d3Sx3, d3Sy3, d3Sz3
        integer :: isign
        real(MP) :: re40x, re41x, re40y, re41y, re40z, re41z

        call pwatch__on("kernel__update_vel")


#ifdef _OPENACC
        !$acc kernels &
        !$acc pcopyin(bx, by, bz, Sxx, Syy, Szz, Syz, Sxz, Sxy, Vx, Vy, Vz, kfs_top, kfs_bot, kob_top, kob_bot)
        !$acc loop independent collapse(3) 
#else
        !$omp parallel &
        !$omp private( d3Sx3, d3Sy3, d3Sz3 ) &
        !$omp private( i, j, k) &
        !$omp private(re40x, re41x, re40y, re41y, re40z, re41z, isign)
        !$omp do &
        !$omp schedule(static,1)
#endif
        do j = jbeg_k, jend_k
            do i = ibeg_k, iend_k
                do k = kbeg_k, kend_k
                
                    isign = sign(1, max((k - kfs_top(i,j)) * (kfs_bot(i,j) - k), &
                                        (k - kob_top(i,j)) * (kob_bot(i,j) - k)))

                    re40x = rc40x + isign * rd40x
                    re41x = rc41x + isign * rd41x
                    re40y = rc40y + isign * rd40y
                    re41y = rc41y + isign * rd41y
                    re40z = rc40z + isign * rd40z
                    re41z = rc41z + isign * rd41z
            
                    d3Sx3 = (Sxx(k  ,i+1,j  ) - Sxx(k  ,i  ,j  )) * re40x - (Sxx(k  ,i+2,j  ) - Sxx(k  ,i-1,j  )) * re41x &
                          + (Sxy(k  ,i  ,j  ) - Sxy(k  ,i  ,j-1)) * re40y - (Sxy(k  ,i  ,j+1) - Sxy(k  ,i  ,j-2)) * re41y &
                          + (Sxz(k  ,i  ,j  ) - Sxz(k-1,i  ,j  )) * re40z - (Sxz(k+1,i  ,j  ) - Sxz(k-2,i  ,j  )) * re41z
                    d3Sy3 = (Sxy(k  ,i  ,j  ) - Sxy(k  ,i-1,j  )) * re40x - (Sxy(k  ,i+1,j  ) - Sxy(k  ,i-2,j  )) * re41x &
                          + (Syy(k  ,i  ,j+1) - Syy(k  ,i  ,j  )) * re40y - (Syy(k  ,i  ,j+2) - Syy(k  ,i  ,j-1)) * re41y &
                          + (Syz(k  ,i  ,j  ) - Syz(k-1,i  ,j  )) * re40z - (Syz(k+1,i  ,j  ) - Syz(k-2,i  ,j  )) * re41z
                    d3Sz3 = (Sxz(k  ,i  ,j  ) - Sxz(k  ,i-1,j  )) * re40x - (Sxz(k  ,i+1,j  ) - Sxz(k  ,i-2,j  )) * re41x &
                          + (Syz(k  ,i  ,j  ) - Syz(k  ,i  ,j-1)) * re40y - (Syz(k  ,i  ,j+1) - Syz(k  ,i  ,j-2)) * re41y &
                          + (Szz(k+1,i  ,j  ) - Szz(k  ,i  ,j  )) * re40z - (Szz(k+2,i  ,j  ) - Szz(k-1,i  ,j  )) * re41z

                    Vx(k, i, j) = Vx(k, i, j) + bx(k, i, j) * d3Sx3 * dt
                    Vy(k, i, j) = Vy(k, i, j) + by(k, i, j) * d3Sy3 * dt
                    Vz(k, i, j) = Vz(k, i, j) + bz(k, i, j) * d3Sz3 * dt

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
         
        call pwatch__off("kernel__update_vel")

    end subroutine kernel__update_vel

    subroutine kernel__update_stress()

        !! Update stress for one time step

        integer :: i, j, k, m
        real(SP) :: mu2, lam2mu
        real(SP) :: taup1, taus1, taup_plus1, taus_plus1
        real(SP) :: d3v3, dxVx_dyVy, dxVx_dzVz, dyVy_dzVz
        real(SP) :: Rxx_n, Ryy_n, Rzz_n, Ryz_n, Rxz_n, Rxy_n
        real(MP) :: dxVx, dyVy, dzVz
        real(MP) :: dxVy_dyVx, dxVz_dzVx, dyVz_dzVy
        integer  :: isign
        real(MP) :: re40x, re41x, re40y, re41y, re40z, re41z

        call pwatch__on("kernel__update_stress")

#ifdef _OPENACC
        !$acc kernels &
        !$acc present(Vx, Vy, Vz,  Sxx, Syy, Szz, Rxx, Ryy, Rzz, &
        !$acc         mu, lam, taup, taus, c1, c2, d1, d2, &
        !$acc         kfs_top, kfs_bot, kob_top, kob_bot)
        !$acc loop independent collapse(3)
#else
        !$omp parallel  &
        !$omp private( dxVx, dyVy, dzVz ) &
        !$omp private( mu2, lam2mu ) &
        !$omp private( taup1, taus1, taup_plus1, taus_plus1 ) &
        !$omp private( d3v3, dyVy_dzVz, dxVx_dzVz, dxVx_dyVy ) &
        !$omp private( Rxx_n, Ryy_n, Rzz_n ) &
        !$omp private( re40x, re41x, re40y, re41y, re40z, re41z, isign) &
        !$omp private( i, j, k, m )
        !$omp do &
        !$omp schedule(static,1)
#endif
        do j = jbeg_k, jend_k
            do i = ibeg_k, iend_k

                !ocl unroll('full')
                !ocl swp
                !OCL SWP_IREG_RATE(200)
                do k = kbeg_k, kend_k


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

                    !! medium copy
                    mu2 = 2 * mu(k, i, j)
                    lam2mu = lam(k, i, j) + mu2

                    taup1 = taup(k, i, j)
                    taus1 = taus(k, i, j)

                    !! update memory variables

                    !! working variables for combinations of velocity derivatives
                    d3v3      = dxVx + dyVy + dzVz
                    dyVy_dzVz = dyVy + dzVz
                    dxVx_dzVz = dxVx + dzVz
                    dxVx_dyVy = dxVx + dyVy
                  
                    Rxx_n = 0.0
                    Ryy_n = 0.0
                    Rzz_n = 0.0

#ifdef _OPENACC
                    !$acc loop seq reduction(+:Rxx_n,Ryy_n,Rzz_n)
#endif
                    do m=1, nm
                      Rxx(m,k,i,j) = c1(m) * Rxx(m,k,i,j) - c2(m) * (lam2mu * taup1 * d3v3 - mu2 * taus1 * dyVy_dzVz) * dt
                      Ryy(m,k,i,j) = c1(m) * Ryy(m,k,i,j) - c2(m) * (lam2mu * taup1 * d3v3 - mu2 * taus1 * dxVx_dzVz) * dt
                      Rzz(m,k,i,j) = c1(m) * Rzz(m,k,i,j) - c2(m) * (lam2mu * taup1 * d3v3 - mu2 * taus1 * dxVx_dyVy) * dt
                      Rxx_n = Rxx_n + d1(m) * Rxx(m,k,i,j)
                      Ryy_n = Ryy_n + d1(m) * Ryy(m,k,i,j)
                      Rzz_n = Rzz_n + d1(m) * Rzz(m,k,i,j)
                    end do

                    !! update stress components

                    taup_plus1 = 1 + taup1 * ( 1 + d2 )
                    taus_plus1 = 1 + taus1 * ( 1 + d2 )
          
                    Sxx (k,i,j) = Sxx (k,i,j) + (lam2mu * taup_plus1 * d3v3 - mu2 * taus_plus1 * dyVy_dzVz + Rxx_n) * dt
                    Syy (k,i,j) = Syy (k,i,j) + (lam2mu * taup_plus1 * d3v3 - mu2 * taus_plus1 * dxVx_dzVz + Ryy_n) * dt
                    Szz (k,i,j) = Szz (k,i,j) + (lam2mu * taup_plus1 * d3v3 - mu2 * taus_plus1 * dxVx_dyVy + Rzz_n) * dt
                    
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
        !$acc present(Vx, Vy, Vz, Sxy, Sxz, Syz, Rxy, Rxz, Ryz, &
        !$acc         muxy, muxz, muyz, taus, c1, c2, d1, d2, &
        !$acc         kfs_top, kfs_bot, kob_top, kob_bot)
        !$acc loop independent collapse(3)
#else
        !$omp parallel &
        !$omp private( dxVy_dyVx, dxVz_dzVx, dyVz_dzVy ) &
        !$omp private( taus1, taus_plus1 ) &
        !$omp private( Ryz_n, Rxz_n, Rxy_n ) &
        !$omp private( re40x, re41x, re40y, re41y, re40z, re41z, isign) &
        !$omp private( i, j, k, m )
        !$omp do  &
        !$omp schedule(static,1)
#endif
        do j = jbeg_k, jend_k
            do i = ibeg_k, iend_k

                !ocl unroll('full')
                !ocl swp
                !OCL SWP_IREG_RATE(300)
                !OCL SWP_FREG_RATE(200)
                do k=kbeg_k, kend_k

                    isign = sign(1, max((k - kfs_top(i,j)) * (kfs_bot(i,j) - k), &
                                        (k - kob_top(i,j)) * (kob_bot(i,j) - k)))

                    re40x = rc40x + isign * rd40x
                    re41x = rc41x + isign * rd41x
                    re40y = rc40y + isign * rd40y
                    re41y = rc41y + isign * rd41y
                    re40z = rc40z + isign * rd40z
                    re41z = rc41z + isign * rd41z

                    dxVy_dyVx = (Vy(k  ,i+1,j  ) - Vy(k  ,i  ,j  )) * re40x - (Vy(k  ,i+2,j  ) - Vy(k  ,i-1,j  )) * re41x &
                              + (Vx(k  ,i  ,j+1) - Vx(k  ,i  ,j  )) * re40y - (Vx(k  ,i  ,j+2) - Vx(k  ,i  ,j-1)) * re41y
                    dxVz_dzVx = (Vz(k  ,i+1,j  ) - Vz(k  ,i  ,j  )) * re40x - (Vz(k  ,i+2,j  ) - Vz(k  ,i-1,j  )) * re41x &
                              + (Vx(k+1,i  ,j  ) - Vx(k  ,i  ,j  )) * re40z - (Vx(k+2,i  ,j  ) - Vx(k-1,i  ,j  )) * re41z
                    dyVz_dzVy = (Vz(k  ,i  ,j+1) - Vz(k  ,i  ,j  )) * re40y - (Vz(k  ,i  ,j+2) - Vz(k  ,i  ,j-1)) * re41y &
                              + (Vy(k+1,i  ,j  ) - Vy(k  ,i  ,j  )) * re40z - (Vy(k+2,i  ,j  ) - Vy(k-1,i  ,j  )) * re41z

                    !! medium copy
                    taus1 = taus(k, i, j)

                    !! update memory variables
                    Ryz_n = 0.0
                    Rxz_n = 0.0
                    Rxy_n = 0.0
#ifdef _OPENACC
                    !$acc loop seq reduction(+:Rxy_n,Ryz_n,Rxz_n)
#endif
                    do m = 1, nm
                        Ryz(m,k,i,j) = c1(m) * Ryz(m,k,i,j) - c2(m) * muyz(k,i,j) * taus1 * dyVz_dzVy * dt
                        Rxz(m,k,i,j) = c1(m) * Rxz(m,k,i,j) - c2(m) * muxz(k,i,j) * taus1 * dxVz_dzVx * dt
                        Rxy(m,k,i,j) = c1(m) * Rxy(m,k,i,j) - c2(m) * muxy(k,i,j) * taus1 * dxVy_dyVx * dt
                        Ryz_n = Ryz_n + d1(m) * Ryz(m,k,i,j)
                        Rxz_n = Rxz_n + d1(m) * Rxz(m,k,i,j)
                        Rxy_n = Rxy_n + d1(m) * Rxy(m,k,i,j)
                    end do

                    !! update stress components
                    taus_plus1 = 1 + taus1 * (1 + d2)

                    Syz (k,i,j) = Syz (k,i,j) + (muyz(k,i,j) * taus_plus1 * dyVz_dzVy + Ryz_n) * dt
                    Sxz (k,i,j) = Sxz (k,i,j) + (muxz(k,i,j) * taus_plus1 * dxVz_dzVx + Rxz_n) * dt
                    Sxy (k,i,j) = Sxy (k,i,j) + (muxy(k,i,j) * taus_plus1 * dxVy_dyVx + Rxy_n) * dt
                              
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

        call pwatch__off("kernel__update_stress")

    end subroutine kernel__update_stress

    subroutine kernel__vmax(xmax, ymax, zmax)

        !! maximum value for terminal output

        real(SP), intent(out) :: xmax, ymax, zmax
        integer :: i, j
        integer, parameter :: margin = 5
        real(SP) :: xmax_local, ymax_local, zmax_local

        !! avoid nearby the absorbing boundary
        !$acc kernels present(Vx, Vy, Vz, kob) copyout(xmax, ymax, zmax)
        xmax = 0.0
        ymax = 0.0
        zmax = 0.0
        !$acc loop reduction(max:xmax, ymax, zmax)
        do j = max(na + margin + 1, jbeg_k), min(ny - na - margin, jend_k)
            do i = max(na + margin + 1, ibeg_k), min(nx - na - margin, iend_k)
                xmax = max(xmax, real(abs(vx(kob(i, j) + 1, i, j))))
                ymax = max(ymax, real(abs(vy(kob(i, j) + 1, i, j))))
                zmax = max(zmax, real(abs(vz(kob(i, j) + 1, i, j))))
            end do
        end do
        !$acc end kernels

    end subroutine kernel__vmax

    subroutine memory_allocate

        !! memory allocation

        allocate (Vx(kbeg_m:kend_m, ibeg_m:iend_m, jbeg_m:jend_m), source=0.0_MP)
        allocate (Vy(kbeg_m:kend_m, ibeg_m:iend_m, jbeg_m:jend_m), source=0.0_MP)
        allocate (Vz(kbeg_m:kend_m, ibeg_m:iend_m, jbeg_m:jend_m), source=0.0_MP)
        allocate (Sxx(kbeg_m:kend_m, ibeg_m:iend_m, jbeg_m:jend_m), source=0.0_MP)
        allocate (Syy(kbeg_m:kend_m, ibeg_m:iend_m, jbeg_m:jend_m), source=0.0_MP)
        allocate (Szz(kbeg_m:kend_m, ibeg_m:iend_m, jbeg_m:jend_m), source=0.0_MP)
        allocate (Syz(kbeg_m:kend_m, ibeg_m:iend_m, jbeg_m:jend_m), source=0.0_MP)
        allocate (Sxz(kbeg_m:kend_m, ibeg_m:iend_m, jbeg_m:jend_m), source=0.0_MP)
        allocate (Sxy(kbeg_m:kend_m, ibeg_m:iend_m, jbeg_m:jend_m), source=0.0_MP)

        if (nm > 0) then
            allocate (c1(nm), c2(nm), d1(nm))
            allocate (Rxx(1:nm, kbeg_k:kend_k, ibeg_k:iend_k, jbeg_k:jend_k), source=0.0)
            allocate (Ryy(1:nm, kbeg_k:kend_k, ibeg_k:iend_k, jbeg_k:jend_k), source=0.0)
            allocate (Rzz(1:nm, kbeg_k:kend_k, ibeg_k:iend_k, jbeg_k:jend_k), source=0.0)
            allocate (Ryz(1:nm, kbeg_k:kend_k, ibeg_k:iend_k, jbeg_k:jend_k), source=0.0)
            allocate (Rxz(1:nm, kbeg_k:kend_k, ibeg_k:iend_k, jbeg_k:jend_k), source=0.0)
            allocate (Rxy(1:nm, kbeg_k:kend_k, ibeg_k:iend_k, jbeg_k:jend_k), source=0.0)
        end if

    end subroutine memory_allocate

end module m_kernel
