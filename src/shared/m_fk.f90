module m_fk

    !! Frequency & Wavenumber FFTs
    !!
    !! Copyright 2013-2024 Takuto Maeda. All rights reserved. This project is released under the MIT license.

    use m_std
    use m_rfft

    private

    public :: fk__t2f     ! time-to-frequency FFT
    public :: fk__f2t     ! frequency-to-time FFT
    public :: fk__x2k_2d  ! space to wavenumber FFT (2D)
    public :: fk__k2x_2d  ! wavenumber to space FFT (2D)
    public :: fk__x2k_3d  ! space to wavenumber FFT (3D)
    public :: fk__k2x_3d  ! wavenumber to space FFT (3D)

    interface fk__t2f
        !! Time-to-frequqency domain Fourier transform using FFT
        module procedure t2f_8, t2f_4
    end interface fk__t2f

    interface fk__f2t
        !! frequqency-to-time domain Fourier transform using FFT. it is an inverse of fk__t2f.
        module procedure f2t_8, f2t_4
    end interface fk__f2t

    interface fk__x2k_2d
        !! Space-to-wavenumber Fourier transform using FFT
        module procedure x2k2d_8, x2k2d_4
    end interface fk__x2k_2d

    interface fk__k2x_2d
        !! Wavenumber-to-space Fourier transform using FFT
        module procedure k2x2d_8, k2x2d_4
    end interface fk__k2x_2d

    interface fk__x2k_3d
        !! Space-to-wavenumber Fourier transform using FFT
        module procedure x2k3d_8, x2k3d_4
    end interface fk__x2k_3d

    interface fk__k2x_3d
        !! Wavenumber-to-space Fourier transform using FFT
        module procedure k2x3d_8, k2x3d_4
    end interface fk__k2x_3d

contains

    subroutine t2f_8(n, dt, r, c, f)

        !! Time-to-frequqency domain Fourier transform using FFT
        !!
        !! Calculate
        !!  c(f) = \int r(t) exp(i 2 \pi f t) dt
        !! at discretized frequencies f(1:n)

        integer,     intent(in)  :: n               !! number of grid (power of 2)
        real(DP),    intent(in)  :: dt              !! sampling period
        real(DP),    intent(in)  :: r(n)            !! time series
        complex(DP), intent(out) :: c(n / 2 + 1)    !! frequency spectrum (unit [r]*s)
        real(DP),    intent(out) :: f(n / 2 + 1)    !! frequency (Hz)
        integer :: i

        call rfft__1drf(n, r, c, 1)

        do i = 1, n / 2 + 1
            f(i) = dble(i - 1) / (dble(n) * dt)
            c(i) = c(i) * dt
        end do

    end subroutine t2f_8


    subroutine t2f_4(n, dt, r, c, f)

        !! Time-to-frequqency domain Fourier transform using FFT
        !!
        !! Calculate
        !!  c(f) = \int r(t) exp(i 2 \pi f t) dt
        !! at discretized frequencies f(1:n)

        integer,     intent(in)  :: n               !! number of grid (power of 2)
        real(SP),    intent(in)  :: dt              !! sampling period
        real(SP),    intent(in)  :: r(n)            !! time series
        complex(SP), intent(out) :: c(n / 2 + 1)    !! frequency spectrum (unit [r]*s)
        real(SP),    intent(out) :: f(n / 2 + 1)    !! frequency (Hz)

        integer :: i

        call rfft__1drf(n, r, c, 1)

        do i = 1, n / 2 + 1
            f(i) = dble(i - 1) / (dble(n) * dt)
            c(i) = c(i) * dt
        end do

    end subroutine t2f_4


    subroutine f2t_8(n, dt, c, r)

        !! Frequency-to-Time domain Fourier transform using FFT. Inverse of t2f
        !!
        !! Calculate
        !!  r(t) = \int c(f) exp(-i 2 \pi f t) df
        !! Result is real value

        integer,     intent(in)  :: n          !! number of grid (power of 2)
        real(DP),    intent(in)  :: dt         !! sampling period
        complex(DP), intent(in)  :: c(n/2 + 1) !! frequency spectrum (unit [r]*s)
        real(DP),    intent(out) :: r(n)       !! time series

        call rfft__1dri(n, c, r, -1)
        r(1:n) = r(1:n) / (n * dt)

    end subroutine f2t_8


    subroutine f2t_4(n, dt, c, r)

        !! Frequency-to-Time domain Fourier transform using FFT. Inverse of t2f
        !!
        !! Calculate
        !!  r(t) = \int c(f) exp(-i 2 \pi f t) df
        !! Result is real value

        integer,     intent(in)  :: n          !! number of grid (power of 2)
        real(SP),    intent(in)  :: dt         !! sampling period
        complex(SP), intent(in)  :: c(n/2 + 1) !! frequency spectrum (unit [r]*s)
        real(SP),    intent(out) :: r(n)       !! time series

        call rfft__1dri(n, c, r, -1)
        r(1:n) = r(1:n) / (n * dt)

    end subroutine f2t_4


    subroutine x2k2d_8(nx, ny, dx, dy, a, b, kx, ky)

        !! Space-to-wavenumber Fourier transform using FFT
        !!
        !! Calculate
        !!  b(kx,ky) = \int a(x,y) exp(- i (kx x + ky y) ) dx dy

        integer,     intent(in)  :: nx, ny
        real(DP),    intent(in)  :: dx, dy
        real(DP),    intent(in)  :: a(nx, ny)
        complex(DP), intent(out) :: b(nx, ny)
        real(DP),    intent(out) :: kx(nx)
        real(DP),    intent(out) :: ky(ny)

        integer :: i, j

        b(1:nx, 1:ny) = dcmplx(a(1:nx, 1:ny))
        call rfft__2d(nx, ny, b, -1)

        b(1:nx, 1:ny) = b(1:nx, 1:ny) * (dx * dy)

        kx(1:nx / 2 + 1) = (/(2 * PI * (i - 1) / dble(nx * dx), i=1, nx / 2 + 1)/)
        ky(1:ny / 2 + 1) = (/(2 * PI * (j - 1) / dble(ny * dy), j=1, ny / 2 + 1)/)
        kx(nx / 2 + 2:nx) = (/(-2 * PI * (nx - i + 1) / dble(nx * dx), i=nx / 2 + 2, nx)/)
        ky(ny / 2 + 2:ny) = (/(-2 * PI * (ny - j + 1) / dble(ny * dy), j=ny / 2 + 2, ny)/)

    end subroutine x2k2d_8


    subroutine k2x2d_8(nx, ny, dx, dy, b, a)

        !! Wavenumber-to-space Fourier transform using FFT
        !!
        !! Calculate
        !!  a(x,y) = \int b(kx,ky) exp( i (kx x + ky y) ) dkx dky / ( 2 \pi )^2
        
        integer,     intent(in)    :: nx, ny
        real(DP),    intent(in)    :: dx, dy
        complex(DP), intent(inout) :: b(nx, ny)
        real(DP),    intent(out)   :: a(nx, ny)

        integer :: i, j

        call rfft__2d(nx, ny, b, 1)
        do j = 1, ny
            do i = 1, nx
                a(i, j) = dble(b(i, j)) / (nx * dx * ny * dy)
            end do
        end do

    end subroutine k2x2d_8


    subroutine x2k2d_4(nx, ny, dx, dy, a, b, kx, ky)

        !! Space-to-wavenumber Fourier transform using FFT
        !!
        !! Calculate
        !!  b(kx,ky) = \int a(x,y) exp(- i (kx x + ky y) ) dx dy

        integer,     intent(in)  :: nx, ny
        real(SP),    intent(in)  :: dx, dy
        real(SP),    intent(in)  :: a(nx, ny)
        complex(SP), intent(out) :: b(nx, ny)
        real(SP),    intent(out) :: kx(nx)
        real(SP),    intent(out) :: ky(ny)

        integer :: i, j

        b(1:nx, 1:ny) = dcmplx(a(1:nx, 1:ny))
        call rfft__2d(nx, ny, b, -1)

        b(1:nx, 1:ny) = b(1:nx, 1:ny) * (dx * dy)

        kx(1:nx / 2 + 1) = (/(2 * PI * (i - 1) / real(nx * dx), i=1, nx / 2 + 1)/)
        ky(1:ny / 2 + 1) = (/(2 * PI * (j - 1) / real(ny * dy), j=1, ny / 2 + 1)/)
        kx(nx / 2 + 2:nx) = (/(-2 * PI * (nx - i + 1) / real(nx * dx), i=nx / 2 + 2, nx)/)
        ky(ny / 2 + 2:ny) = (/(-2 * PI * (ny - j + 1) / real(ny * dy), j=ny / 2 + 2, ny)/)

    end subroutine x2k2d_4


    subroutine k2x2d_4(nx, ny, dx, dy, b, a)

        !! Wavenumber-to-space Fourier transform using FFT
        !!
        !! Calculate
        !!  a(x,y) = \int b(kx,ky) exp( i (kx x + ky y) ) dkx dky / ( 2 \pi )^2

        integer,     intent(in)     :: nx, ny
        real(SP),    intent(in)     :: dx, dy
        complex(SP), intent(inout)  :: b(nx, ny)
        real(SP),    intent(out)    :: a(nx, ny)

        integer :: i, j

        call rfft__2d(nx, ny, b, 1)

        do j = 1, ny
            do i = 1, nx
                a(i, j) = dble(b(i, j)) / (nx * dx * ny * dy)
            end do
        end do

    end subroutine k2x2d_4


    subroutine x2k3d_8(nx, ny, nz, dx, dy, dz, a, b, kx, ky, kz)

        !! Space-to-wavenumber Fourier transform using FFT
        !!
        !! Calculate
        !!  b(kx,ky,kz) = \int a(x,y,z) exp(- i (kx x + ky y + kz z ) ) dx dy dz

        integer,     intent(in)  :: nx, ny, nz
        real(DP),    intent(in)  :: dx, dy, dz
        real(DP),    intent(in)  :: a(nx, ny, nz)
        complex(DP), intent(out) :: b(nx, ny, nz)
        real(DP),    intent(out) :: kx(nx)
        real(DP),    intent(out) :: ky(ny)
        real(DP),    intent(out) :: kz(nz)

        integer :: i, j, k

        b(1:nx, 1:ny, 1:nz) = dcmplx(a(1:nx, 1:ny, 1:nz))
        call rfft__3d(nx, ny, nz, b, -1)

        b(1:nx, 1:ny, 1:nz) = b(1:nx, 1:ny, 1:nz) * (dx * dy * dz)

        kx(1:nx / 2 + 1) = (/(2 * PI * (i - 1) / dble(nx * dx), i=1, nx / 2 + 1)/)
        ky(1:ny / 2 + 1) = (/(2 * PI * (j - 1) / dble(ny * dy), j=1, ny / 2 + 1)/)
        kz(1:nz / 2 + 1) = (/(2 * PI * (k - 1) / dble(nz * dz), k=1, nz / 2 + 1)/)
        kx(nx / 2 + 2:nx) = (/(-2 * PI * (nx - i + 1) / dble(nx * dx), i=nx / 2 + 2, nx)/)
        ky(ny / 2 + 2:ny) = (/(-2 * PI * (ny - j + 1) / dble(ny * dy), j=ny / 2 + 2, ny)/)
        kz(nz / 2 + 2:nz) = (/(-2 * PI * (nz - k + 1) / dble(nz * dz), k=nz / 2 + 2, nz)/)

    end subroutine x2k3d_8


    subroutine x2k3d_4(nx, ny, nz, dx, dy, dz, a, b, kx, ky, kz)

        !! Space-to-wavenumber Fourier transform using FFT
        !!
        !! Calculate
        !!  b(kx,ky,kz) = \int a(x,y,z) exp(- i (kx x + ky y + kz z ) ) dx dy dz

        integer,     intent(in)  :: nx, ny, nz
        real(SP),    intent(in)  :: dx, dy, dz
        real(SP),    intent(in)  :: a(nx, ny, nz)
        complex(SP), intent(out) :: b(nx, ny, nz)
        real(SP),    intent(out) :: kx(nx)
        real(SP),    intent(out) :: ky(ny)
        real(SP),    intent(out) :: kz(nz)

        integer :: i, j, k

        b(1:nx, 1:ny, 1:nz) = dcmplx(a(1:nx, 1:ny, 1:nz))
        call rfft__3d(nx, ny, nz, b, -1)

        b(1:nx, 1:ny, 1:nz) = b(1:nx, 1:ny, 1:nz) * (dx * dy * dz)

        kx(1:nx / 2 + 1) = (/(2 * PI * (i - 1) / real(nx * dx), i=1, nx / 2 + 1)/)
        ky(1:ny / 2 + 1) = (/(2 * PI * (j - 1) / real(ny * dy), j=1, ny / 2 + 1)/)
        kz(1:nz / 2 + 1) = (/(2 * PI * (k - 1) / real(nz * dz), k=1, nz / 2 + 1)/)
        kx(nx / 2 + 2:nx) = (/(-2 * PI * (nx - i + 1) / real(nx * dx), i=nx / 2 + 2, nx)/)
        ky(ny / 2 + 2:ny) = (/(-2 * PI * (ny - j + 1) / real(ny * dy), j=ny / 2 + 2, ny)/)
        kz(nz / 2 + 2:nz) = (/(-2 * PI * (nz - k + 1) / real(nz * dz), k=nz / 2 + 2, nz)/)

    end subroutine x2k3d_4


    subroutine k2x3d_8(nx, ny, nz, dx, dy, dz, b, a)

        !! Space-to-wavenumber Fourier transform using FFT
        !!
        !! Calculate
        !!  a(x,y,z) = \int b(kx,ky,kz) exp( i (kx x + ky y + kz z) ) dkx dky dkz
        !!           / (2 \pi)^3

        integer,     intent(in)     :: nx, ny, nz
        real(DP),    intent(in)     :: dx, dy, dz
        complex(DP), intent(inout)  :: b(nx, ny, nz)
        real(DP),    intent(out)    :: a(nx, ny, nz)

        integer :: i, j, k

        call rfft__3d(nx, ny, nz, b, 1)

        do k = 1, nz
            do j = 1, ny
                do i = 1, nx
                    a(i, j, k) = dble(b(i, j, k)) / (nx * dx * ny * dy * nz * dz)
                end do
            end do
        end do

    end subroutine k2x3d_8


    subroutine k2x3d_4(nx, ny, nz, dx, dy, dz, b, a)

        !! Space-to-wavenumber Fourier transform using FFT
        !!
        !! Calculate
        !!    a(x,y,z) = \int b(kx,ky,kz) exp( i (kx x + ky y + kz z) ) dkx dky dkz
        !!             / (2 \pi)^3

        integer,     intent(in)     :: nx, ny, nz
        real(SP),    intent(in)     :: dx, dy, dz
        complex(SP), intent(inout)  :: b(nx, ny, nz)
        real(SP),    intent(out)    :: a(nx, ny, nz)

        integer :: i, j, k

        call rfft__3d(nx, ny, nz, b, 1)

        do k = 1, nz
            do j = 1, ny
                do i = 1, nx
                    a(i, j, k) = real(b(i, j, k)) / (nx * dx * ny * dy * nz * dz)
                end do
            end do
        end do

    end subroutine k2x3d_4

    
end module m_fk
