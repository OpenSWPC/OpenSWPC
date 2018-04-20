!! ----------------------------------------------------------------------------------------------------------------------------- !!
!>
!! Calculate power spectrum density functions (PSDF) and synthetizing realization of random media
!!
!! @copyright
!!   Copyright 2013-2018 Takuto Maeda. All rights reserved. This project is released under the MIT license.
!<
!! ----
module m_rmedia

  !! -- Dependency
  use m_std
  use m_fk
  use m_gammaf

  !! -- Declaraitons
  implicit none
  private

  !! -- Public procedures
  public :: rmedia__3dgen                !< A realization of 3D random media
  public :: rmedia__3dgen_ani            !< A realization of 3D anisomeric random media
  public :: rmedia__2dgen                !< A realization of 2D random media
  public :: rmedia__2dgen_ani            !< A realization of 2D anisomeric random media
  public :: rmedia__psdf_gauss1d         !< Gaussian PSDF (1D)
  public :: rmedia__psdf_gauss2d         !< Gaussian PSDF (2D)
  public :: rmedia__psdf_gauss3d         !< Gaussian PSDF (3D)
  public :: rmedia__psdf_exp1d           !< Exponential PSDF (1D)
  public :: rmedia__psdf_exp2d           !< Exponential PSDF (2D)
  public :: rmedia__psdf_exp3d           !< Exponential PSDF (3D)
  public :: rmedia__psdf_vonKarman1d     !< von Karman PSDF (1D)
  public :: rmedia__psdf_vonKarman2d     !< von Karman PSDF (2D)
  public :: rmedia__psdf_vonKarman3d     !< von Karman PSDF (3D)
  public :: rmedia__psdf_ani_gauss2d     !< Gaussian Anisomeric PSDF (2D)
  public :: rmedia__psdf_ani_exp2d       !< Gaussian Anisomeric PSDF (2D)
  public :: rmedia__psdf_ani_vonKarman2d !< Gaussian Anisomeric PSDF (2D)
  public :: rmedia__psdf_ani_gauss3d     !< Gaussian PSDF (3D)
  public :: rmedia__psdf_ani_exp3d       !< Exponential PSDF (3D)
  public :: rmedia__psdf_ani_vonKarman3d !< von Karman PSDF (3D)

  !! -- Local Parameters

  integer, parameter :: ptype_gauss     = 1
  integer, parameter :: ptype_exp       = 2
  integer, parameter :: ptype_vonKarman = 3

  !! -----

contains

  !! --------------------------------------------------------------------------------------------------------------------------- !!
  !>
  !! generate 3D random media having specific power spectrum density function
  !!
  !! @par Method
  !!   1. Genrerate random variable n*n array
  !!   2. Calculate spectrum by using 2D fft
  !!   3. Whitening (let spectrum amplitude unity)
  !!   4. Filtering by using specified PSDF
  !!   5. Get random media by ifft
  !!
  !! @par Input parameter
  !!    ptype = 1: Gaussian /  2: Exponential /  3: von Karman
  !!
  !! @see Roth and Korn (1992), doi:10.1111/j.1365-246X.1993.tb01442.x
  !<
  !! ---
  subroutine rmedia__3dgen_ani( ptype, nx, ny, nz, dx, dy, dz, ax, ay, az, epsil, kappa, x, y, z,  media, seed )

    !! -- Arguments
    integer,  intent(in)  :: ptype           !< type of random media
    integer,  intent(in)  :: nx              !< size of x-direction
    integer,  intent(in)  :: ny              !< size of y-direction
    integer,  intent(in)  :: nz              !< size of z-direction
    real(SP), intent(in)  :: dx              !< grid width of x-direction
    real(SP), intent(in)  :: dy              !< grid width of y-direction
    real(SP), intent(in)  :: dz              !< grid width of y-direction
    real(SP), intent(in)  :: ax              !< correlation length of random media
    real(SP), intent(in)  :: ay              !< correlation length of random media
    real(SP), intent(in)  :: az              !< correlation length of random media
    real(SP), intent(in)  :: epsil           !< fractional fluctuation of random media
    real(SP), intent(in)  :: kappa           !< von Karman parameter; not used for gaussian/exponential type
    real(SP), intent(out) :: x(nx)           !< x dataloc array
    real(SP), intent(out) :: y(ny)           !< y dataloc array
    real(SP), intent(out) :: z(nz)           !< z dataloc array
    real(SP), intent(out) :: media(nx,ny,nz) !< random media
    integer,  intent(in), optional :: seed   !< seed number of random variable; Default automatically generate seed from clock
    !! --
    complex(SP), allocatable :: bb(:,:,:)
    real(SP)    :: kx(nx), ky(ny), kz(nz)     ! wavenumber vector
    integer     :: ic
    real(SP)    :: psdf
    integer     :: i, j, k
    real(SP)    :: m
    real(SP)    :: aepsil
    integer, allocatable :: seedval(:)
    integer :: nseed
    !! ----

    allocate( bb(nx,ny,nz) )

    ! initialize random seed by specified value or clock
    call random_seed( size=nseed )
    allocate(seedval(nseed))
    if ( present( seed ) ) then
      seedval = seed
    else
      call system_clock ( count = ic )
      seedval = ic
    end if
    call random_seed( put=seedval )
    deallocate( seedval )

    ! random value
    do k=1, nz
      do j=1, ny
        do i=1, nx
          call random_number( media(i,j,k) )
        end do
      end do
    end do

    call fk__x2k_3d( nx, ny, nz, dx, dy, dz, media, bb, kx, ky, kz )
    do k=1, nz
      do j=1, ny
        do i=1, nx

          !! spectrum whitening
          bb(i,j,k) = bb(i,j,k) / abs( bb(i,j,k) )

          !! obtain PSDF
          m = sqrt( kx(i)**2 + ky(j)**2 + kz(k)**2 )

          select case ( ptype )
          case ( ptype_gauss )
            call rmedia__psdf_ani_gauss3d ( kx(i), ky(j), kz(k), ax, ay, az, epsil, psdf )
          case ( ptype_exp )
            call rmedia__psdf_ani_exp3d ( kx(i), ky(j), kz(k), ax, ay, az, epsil, psdf )
          case ( ptype_vonKarman )
            call rmedia__psdf_ani_vonkarman3d ( kx(i), ky(j), kz(k), ax, ay, az, epsil, kappa, psdf )
          case default
            write(STDERR,'(A)') 'WARNING [rmedia__3dgen]: no such psdf type; assume von Karman'
            call rmedia__psdf_ani_vonkarman3d ( kx(i), ky(j), kz(k), ax, ay, az, epsil, kappa, psdf )
          end select

          ! filetering
          bb(i,j,k) = bb(i,j,k) * sqrt(psdf*nx*ny*nz*dx*dy*dz)
        end do
      end do
    end do

    ! let DC component zero
    bb(1,1,1)=0.0_dp

    ! inverse transform
    call fk__k2x_3d( nx, ny, nz, dx, dy, dz, bb, media )

    ! remove mean
    media = media - sum(media)/dble(nx*ny*nz)

    ! epsilon correction so that  epsilon^2 = sum( media^2 )
    aepsil = 0.0_dp
    do k=1, nz
      do j=1, ny
        do i=1, nx
          aepsil = aepsil + media(i,j,k)**2
        end do
      end do
    end do

    aepsil = sqrt( aepsil / dble( nx*ny*nz ) )

    media = media * epsil/aepsil

    ! axis
    do i=1, nx
      x(i) = (i-1)*dx
    end do
    do j=1, ny
      y(j) = (j-1)*dy
    end do
    do k=1, nz
      z(k) = (k-1)*dz
    end do

    deallocate( bb )

  end subroutine rmedia__3dgen_ani
  !! --------------------------------------------------------------------------------------------------------------------------- !!

  !! --------------------------------------------------------------------------------------------------------------------------- !!
  !>
  !! generate 3D random media having specific power spectrum density function
  !!
  !! @par Method
  !!   1. Genrerate random variable n*n array
  !!   2. Calculate spectrum by using 2D fft
  !!   3. Whitening (let spectrum amplitude unity)
  !!   4. Filtering by using specified PSDF
  !!   5. Get random media by ifft
  !!
  !! @par Input parameter
  !!    ptype = 1: Gaussian /  2: Exponential /  3: von Karman
  !!
  !! @see Roth and Korn (1992), doi:10.1111/j.1365-246X.1993.tb01442.x
  !<
  !! ---
  subroutine rmedia__3dgen( ptype, nx, ny, nz, dx, dy, dz, a, epsil, kappa, x, y, z,  media, seed )

    !! -- Arguments
    integer,  intent(in)  :: ptype           !< type of random media
    integer,  intent(in)  :: nx              !< size of x-direction
    integer,  intent(in)  :: ny              !< size of y-direction
    integer,  intent(in)  :: nz              !< size of z-direction
    real(SP), intent(in)  :: dx              !< grid width of x-direction
    real(SP), intent(in)  :: dy              !< grid width of y-direction
    real(SP), intent(in)  :: dz              !< grid width of y-direction
    real(SP), intent(in)  :: a               !< correlation length of random media
    real(SP), intent(in)  :: epsil           !< fractional fluctuation of random media
    real(SP), intent(in)  :: kappa           !< von Karman parameter; not used for gaussian/exponential type
    real(SP), intent(out) :: x(nx)           !< x dataloc array
    real(SP), intent(out) :: y(ny)           !< y dataloc array
    real(SP), intent(out) :: z(nz)           !< z dataloc array
    real(SP), intent(out) :: media(nx,ny,nz) !< random media
    integer,  intent(in), optional :: seed   !< seed number of random variable; Default automatically generate seed from clock
    !! --
    complex(SP), allocatable :: bb(:,:,:)
    real(SP)    :: kx(nx), ky(ny), kz(nz)     ! wavenumber vector
    integer     :: ic
    real(SP)    :: psdf
    integer     :: i, j, k
    real(SP)    :: m
    real(SP)    :: aepsil
    integer, allocatable :: seedval(:)
    integer :: nseed
    !! ----

    allocate( bb(nx,ny,nz) )

    ! initialize random seed by specified value or clock
    call random_seed( size=nseed )
    allocate(seedval(nseed))
    if ( present( seed ) ) then
      seedval = seed
    else
      call system_clock ( count = ic )
      seedval = ic
    end if
    call random_seed( put=seedval )
    deallocate( seedval )

    ! random value
    do k=1, nz
      do j=1, ny
        do i=1, nx
          call random_number( media(i,j,k) )
        end do
      end do
    end do

    call fk__x2k_3d( nx, ny, nz, dx, dy, dz, media, bb, kx, ky, kz )
    do k=1, nz
      do j=1, ny
        do i=1, nx

          !! spectrum whitening
          bb(i,j,k) = bb(i,j,k) / abs( bb(i,j,k) )

          !! obtain PSDF
          m = sqrt( kx(i)**2 + ky(j)**2 + kz(k)**2 )

          select case ( ptype )
          case ( ptype_gauss )
            call rmedia__psdf_gauss3d ( m, a, epsil, psdf )
          case ( ptype_exp )
            call rmedia__psdf_exp3d ( m, a, epsil, psdf )
          case ( ptype_vonKarman )
            call rmedia__psdf_vonkarman3d ( m, a, epsil, kappa, psdf )
          case default
            write(STDERR,'(A)') 'WARNING [rmedia__3dgen]: no such psdf type; assume von Karman'
            call rmedia__psdf_vonkarman3d ( m, a, epsil, kappa, psdf )
          end select

          ! filetering
          bb(i,j,k) = bb(i,j,k) * sqrt(psdf*nx*ny*nz*dx*dy*dz)
        end do
      end do
    end do

    ! let DC component zero
    bb(1,1,1)=0.0_dp

    ! inverse transform
    call fk__k2x_3d( nx, ny, nz, dx, dy, dz, bb, media )

    ! remove mean
    media = media - sum(media)/dble(nx*ny*nz)

    ! epsilon correction so that  epsilon^2 = sum( media^2 )
    aepsil = 0.0_dp
    do k=1, nz
      do j=1, ny
        do i=1, nx
          aepsil = aepsil + media(i,j,k)**2
        end do
      end do
    end do

    aepsil = sqrt( aepsil / dble( nx*ny*nz ) )

    media = media * epsil/aepsil

    ! axis
    do i=1, nx
      x(i) = (i-1)*dx
    end do
    do j=1, ny
      y(j) = (j-1)*dy
    end do
    do k=1, nz
      z(k) = (k-1)*dz
    end do

    deallocate( bb )

  end subroutine rmedia__3dgen
  !! --------------------------------------------------------------------------------------------------------------------------- !!

  !! --------------------------------------------------------------------------------------------------------------------------- !!
  !>
  !! generate 2D anistropic random media having specific power spectrum density function
  !!
  !! @par Method
  !!   1. Genrerate random variable n*n array
  !!   2. Calculate spectrum by using 2D fft
  !!   3. Whitening (let spectrum amplitude unity)
  !!   4. Filtering by using specified PSDF
  !!   5. Get random media by ifft
  !!
  !! @par Input parameter
  !!    ptype = 1: Gaussian /  2: Exponential /  3: von Karman
  !!
  !! @see Roth and Korn (1992), doi:10.1111/j.1365-246X.1993.tb01442.x
  !<
  !! ---
  subroutine rmedia__2dgen_ani ( ptype, nx, ny, dx, dy, ax, ay, epsil, kappa, x, y, media, seed )

    !! ----
    integer,  intent(in)  :: ptype            !< type of random media
    integer,  intent(in)  :: nx               !< size of x-direction
    integer,  intent(in)  :: ny               !< size of y-direction
    real(SP), intent(in)  :: dx               !< grid width of x-direction
    real(SP), intent(in)  :: dy               !< grid width of y-direction
    real(SP), intent(in)  :: ax               !< correlation length of random media
    real(SP), intent(in)  :: ay               !< correlation length of random media
    real(SP), intent(in)  :: epsil            !< fractional fluctuation of random media
    real(SP), intent(in)  :: kappa            !< von Karman parameter
    real(SP), intent(out) :: x(nx)            !< x dataloc array
    real(SP), intent(out) :: y(ny)            !< y dataloc array
    real(SP), intent(out) :: media(nx,ny)     !< random media
    integer,  intent(in), optional :: seed    !< seed number of random variable
    !+
    real(SP) :: kx(nx), ky(ny)     ! wavenumber vector
    integer  :: ic ! seed
    complex(SP), allocatable :: aa(:,:)
    real(SP) :: psdf
    integer :: i,j
    real(SP) :: aepsil
    integer, allocatable :: seedval(:)
    integer :: nseed
    !! --

    allocate( aa(nx,ny) )

    ! initialize random seed by specified value or clock
    call random_seed( size=nseed )
    allocate(seedval(nseed))
    if ( present( seed ) ) then
      seedval = seed
    else
      call system_clock ( count = ic )
      seedval = ic
    end if
    call random_seed( put=seedval )
    deallocate( seedval )


    ! random value
    do i=1, nx
      do j=1, ny
        call random_number( media(i,j) )
      end do
    end do

    call fk__x2k_2d( nx, ny, dx, dy, media, aa, kx, ky )

    do i=1, nx
      do j=1, ny

        ! whitening
        aa(i,j) = aa(i,j) / abs( aa(i,j) )

        if( ptype == ptype_gauss ) then
          call rmedia__psdf_ani_gauss2d( kx(i), ky(j), ax, ay, epsil, psdf )
        else if( ptype == ptype_exp ) then
          call rmedia__psdf_ani_exp2d( kx(i), ky(j), ax, ay, epsil, psdf )
        else
          call rmedia__psdf_ani_vonKarman2d( kx(i), ky(j), ax, ay, epsil,kappa, psdf )
        end if

        ! filetering
        aa(i,j) = aa(i,j) * sqrt(psdf*nx*ny*dx*dy)
      end do
    end do

    ! let DC component zero
    aa(1,1)=0.0

    ! inverse transform
    call fk__k2x_2d( nx, ny, dx, dy, aa, media )

    ! remove mean
    media = media - sum(media)/dble(nx*ny)

    ! epsilon correction so that  epsilon^2 = sum( media^2 )
    aepsil = 0.0
    do i=1, nx
      do j=1, ny
        aepsil = aepsil + media(i,j)**2
      end do
    end do
    aepsil = sqrt( aepsil / dble( nx*ny ) )

    media = media * epsil/aepsil

    do i=1, nx
      x(i) = (i-1)*dx
    end do
    do j=1, ny
      y(j) = (j-1)*dy
    end do

    deallocate( aa )

  end subroutine rmedia__2dgen_ani
  !! --------------------------------------------------------------------------------------------------------------------------- !!

  !! --------------------------------------------------------------------------------------------------------------------------- !!
  subroutine rmedia__2dgen( ptype, nx, ny, dx, dy, a, epsil, kappa, x, y, media, seed )

    !! -- Arguments
    integer,  intent(in)  :: ptype         !< type of random media
    integer,  intent(in)  :: nx            !< size of x-direction
    integer,  intent(in)  :: ny            !< size of y-direction
    real(SP), intent(in)  :: dx            !< grid width of x-direction
    real(SP), intent(in)  :: dy            !< grid width of y-direction
    real(SP), intent(in)  :: a             !< correlation length of random media
    real(SP), intent(in)  :: epsil         !< fractional fluctuation of random media
    real(SP), intent(in)  :: kappa         !< von Karman parameter
    real(SP), intent(out) :: x(nx)         !< x dataloc array
    real(SP), intent(out) :: y(ny)         !< y dataloc array
    real(SP), intent(out) :: media(nx,ny)  !< random media
    integer, intent(in), optional :: seed  !< seed number of random variable

    !! --
    real(SP) :: kx(nx), ky(ny)     ! wavenumber vector
    integer  :: ic ! seed
    complex(SP), allocatable  :: aa(:,:)
    real(SP) :: psdf
    integer :: i,j
    real(SP) :: m
    real(SP) :: aepsil
    integer, allocatable :: seedval(:)
    integer :: nseed
    !! ----

    allocate( aa(nx,ny) )

    ! initialize random seed by specified value or clock
    call random_seed( size=nseed )
    allocate(seedval(nseed))
    if ( present( seed ) ) then
      seedval = seed
    else
      call system_clock ( count = ic )
      seedval = ic
    end if
    call random_seed( put=seedval )
    deallocate( seedval )

    ! random value
    do i=1, nx
      do j=1, ny
        call random_number(media(i,j))
      end do
    end do

    call fk__x2k_2d( nx, ny, dx, dy, media, aa, kx, ky )

    do i=1, nx
      do j=1, ny

        ! whitening
        aa(i,j) = aa(i,j) / abs( aa(i,j) )

        m = sqrt( kx(i)**2 + ky(j)**2 )

        if( ptype == ptype_gauss ) then
          call rmedia__psdf_gauss2d( m, a, epsil, psdf )
        else if( ptype == ptype_exp ) then
          call rmedia__psdf_exp2d( m, a, epsil, psdf )
        else
          call rmedia__psdf_vonKarman2d( m, a, epsil, kappa, psdf )
        end if

        ! filetering
        aa(i,j) = aa(i,j) * sqrt(psdf*nx*ny*dx*dy)
      end do
    end do

    ! let DC component zero
    aa(1,1)=0.0

    ! inverse transform
    call fk__k2x_2d( nx, ny, dx, dy, aa, media )

    ! remove mean
    media = media - sum(media)/dble(nx*ny)

    ! epsilon correction so that  epsilon^2 = sum( media^2 )
    aepsil = 0.0
    do i=1, nx
      do j=1, ny
        aepsil = aepsil + media(i,j)**2
      end do
    end do
    aepsil = sqrt( aepsil / dble( nx*ny ) )

    media = media * epsil/aepsil
    do i=1, nx
      x(i) = (i-1)*dx
    end do
    do j=1, ny
      y(j) = (j-1)*dy
    end do

    deallocate( aa )

  end subroutine rmedia__2dgen
  !! --------------------------------------------------------------------------------------------------------------------------- !!


  !! --------------------------------------------------------------------------------------------------------------------------- !!
  subroutine rmedia__psdf_exp1d ( m, a, epsil, psdf )

    !! -- Arguments
    real(SP), intent(in)  :: m         ! wavenumber
    real(SP), intent(in)  :: a         ! correlation distance
    real(SP), intent(in)  :: epsil     ! rms amplitude of fluctuation
    real(SP), intent(out) :: psdf      ! psdf value
    !! --

    psdf = ( 2 * epsil*epsil * a  )/( 1 + a*a * m*m  )

  end subroutine rmedia__psdf_exp1d
  !! --------------------------------------------------------------------------------------------------------------------------- !!

  !! --------------------------------------------------------------------------------------------------------------------------- !!
  subroutine rmedia__psdf_exp2d ( m, a, epsil, psdf )

    !! -- Arguments
    real(SP), intent(in)  :: m         ! wavenumber
    real(SP), intent(in)  :: a         ! correlation distance
    real(SP), intent(in)  :: epsil     ! rms amplitude of fluctuation
    real(SP), intent(out) :: psdf      ! psdf value
    !! --

    psdf = ( 2 * epsil*epsil * a*a * PI_D )/( 1 + a*a * m*m  )**1.5

  end subroutine rmedia__psdf_exp2d
  !! --------------------------------------------------------------------------------------------------------------------------- !!

  !! --------------------------------------------------------------------------------------------------------------------------- !!
  subroutine rmedia__psdf_ani_gauss2d( mx, mz, ax, az, epsil, psdf )

    !! -- Arguments
    real(SP), intent(in)  :: mx        ! wavenumber
    real(SP), intent(in)  :: mz        ! wavenumber
    real(SP), intent(in)  :: ax        ! correlation distance
    real(SP), intent(in)  :: az        ! correlation distance
    real(SP), intent(in)  :: epsil     ! rms amplitude of fluctuation
    real(SP), intent(out) :: psdf      ! psdf value
    !! --

    psdf = epsil*epsil * PI_D*ax*az * exp( - ( ax*ax*mx*mx + az*az*mz*mz ) / 4.0_SP )

  end subroutine rmedia__psdf_ani_gauss2d
  !! --------------------------------------------------------------------------------------------------------------------------- !!

  !! --------------------------------------------------------------------------------------------------------------------------- !!
  subroutine rmedia__psdf_ani_exp2d( mx, mz, ax, az, epsil, psdf )

    !! -- Arguments
    real(SP), intent(in)  :: mx        ! wavenumber
    real(SP), intent(in)  :: mz        ! wavenumber
    real(SP), intent(in)  :: ax        ! correlation distance
    real(SP), intent(in)  :: az        ! correlation distance
    real(SP), intent(in)  :: epsil     ! rms amplitude of fluctuation
    real(SP), intent(out) :: psdf      ! psdf value
    !! --
    psdf = ( 2*PI_D * epsil*epsil * ax*az ) / ( 1 + ax*ax*mx*mx+ az*az*mz*mz )**1.5_SP

  end subroutine rmedia__psdf_ani_exp2d
  !! --------------------------------------------------------------------------------------------------------------------------- !!


  !! --------------------------------------------------------------------------------------------------------------------------- !!
  subroutine rmedia__psdf_ani_vonKarman2d( mx, mz, ax, az, epsil, kappa, psdf )

    !! -- Arguments
    real(SP), intent(in)  :: mx        ! wavenumber
    real(SP), intent(in)  :: mz        ! wavenumber
    real(SP), intent(in)  :: ax        ! correlation distance
    real(SP), intent(in)  :: az        ! correlation distance
    real(SP), intent(in)  :: epsil     ! rms amplitude of fluctuation
    real(SP), intent(in)  :: kappa     ! rms amplitude of fluctuation
    real(SP), intent(out) :: psdf      ! psdf value
    !! --

    psdf = ( 4*PI_D * kappa * epsil*epsil * ax*az )  / ( 1 + ax*ax*mx*mx+ az*az*mz*mz )**(1.0_SP + kappa)

  end subroutine rmedia__psdf_ani_vonKarman2d
  !! --------------------------------------------------------------------------------------------------------------------------- !!

  !! --------------------------------------------------------------------------------------------------------------------------- !!
  subroutine rmedia__psdf_exp3d ( m, a, epsil, psdf )

    !! -- Arguments
    real(SP), intent(in)  :: m         ! wavenumber
    real(SP), intent(in)  :: a         ! correlation distance
    real(SP), intent(in)  :: epsil     ! rms amplitude of fluctuation
    real(SP), intent(out) :: psdf      ! psdf value
    real(SP) :: wk
    !! --

    wk = 1.0 + a*a*m*m
    psdf = ( 8 * epsil*epsil * a*a*a * PI_D )/( wk*wk )

  end subroutine rmedia__psdf_exp3d
  !! --------------------------------------------------------------------------------------------------------------------------- !!

  !! --------------------------------------------------------------------------------------------------------------------------- !!
  subroutine rmedia__psdf_gauss1d ( m, a, epsil, psdf )

    !! -- Arguments
    real(SP), intent(in)  :: m         ! wavenumber
    real(SP), intent(in)  :: a         ! correlation distance
    real(SP), intent(in)  :: epsil     ! rms amplitude of fluctuation
    real(SP), intent(out) :: psdf      ! psdf value
    !! --

    psdf = epsil*epsil * a * sqrt(PI_D) * exp ( - m*m * a*a * 0.25 )

  end subroutine rmedia__psdf_gauss1d
  !! --------------------------------------------------------------------------------------------------------------------------- !!

  !! --------------------------------------------------------------------------------------------------------------------------- !!
  subroutine rmedia__psdf_gauss2d ( m, a, epsil, psdf )

    !! -- Arguments
    real(SP), intent(in)  :: m         ! wavenumber
    real(SP), intent(in)  :: a         ! correlation distance
    real(SP), intent(in)  :: epsil     ! rms amplitude of fluctuation
    real(SP), intent(out) :: psdf      ! psdf value
    !! --

    psdf = epsil*epsil * a*a * PI_D * exp ( - m*m * a*a * 0.25 )

  end subroutine rmedia__psdf_gauss2d
  !! --------------------------------------------------------------------------------------------------------------------------- !!

  !! --------------------------------------------------------------------------------------------------------------------------- !!
  subroutine rmedia__psdf_gauss3d ( m, a, epsil, psdf )

    !! -- Arguments
    real(SP), intent(in)  :: m         ! wavenumber
    real(SP), intent(in)  :: a         ! correlation distance
    real(SP), intent(in)  :: epsil     ! rms amplitude of fluctuation
    real(SP), intent(out) :: psdf      ! psdf value
    !! --

    psdf = (epsil*epsil) * (a*a*a) * PI_D**1.5 * exp ( - m*m * a*a * 0.25 )

  end subroutine rmedia__psdf_gauss3d
  !! --------------------------------------------------------------------------------------------------------------------------- !!

  !! --------------------------------------------------------------------------------------------------------------------------- !!
  subroutine rmedia__psdf_vonKarman1d ( m, a, epsil, kappa, psdf )

    !! -- Arguments
    real(SP), intent(in)  :: m         ! wavenumber
    real(SP), intent(in)  :: a         ! correlation distance
    real(SP), intent(in)  :: epsil     ! rms amplitude of fluctuation
    real(SP), intent(in)  :: kappa     ! von karman parameter
    real(SP), intent(out) :: psdf      ! psdf value
    !! --

    psdf = ( 2 * sqrt(PI_D) * epsil*epsil * a * gammaf( kappa + 0.5 ) ) / &
        ( gammaf( kappa ) * ( 1+a*a * m*m )**( kappa + 0.5 ) )

  end subroutine rmedia__psdf_vonKarman1d
  !! --------------------------------------------------------------------------------------------------------------------------- !!

  !! --------------------------------------------------------------------------------------------------------------------------- !!
  subroutine rmedia__psdf_vonKarman2d ( m, a, epsil, kappa, psdf )

    !! -- Arguments
    real(SP), intent(in)  :: m         ! wavenumber
    real(SP), intent(in)  :: a         ! correlation distance
    real(SP), intent(in)  :: epsil     ! rms amplitude of fluctuation
    real(SP), intent(in)  :: kappa     ! von karman parameter
    real(SP), intent(out) :: psdf      ! psdf value
    !! --

    psdf = ( 4 * PI_D * a*a * epsil*epsil * kappa ) /  ( 1 + a*a * m*m )**( 1 + kappa )

  end subroutine rmedia__psdf_vonKarman2d
  !! --------------------------------------------------------------------------------------------------------------------------- !!

  !! --------------------------------------------------------------------------------------------------------------------------- !!
  subroutine rmedia__psdf_vonKarman3d ( m, a, epsil, kappa, psdf )

    !! -- Arguments
    real(SP), intent(in)  :: m         ! wavenumber
    real(SP), intent(in)  :: a         ! correlation distance
    real(SP), intent(in)  :: epsil     ! rms amplitude of fluctuation
    real(SP), intent(in)  :: kappa     ! von karman parameter
    real(SP), intent(out) :: psdf      ! psdf value
    !! --

    psdf = ( 8 * PI_D**(1.5) * &
        epsil*epsil * a*a*a * gammaf( kappa + 1.5)  ) / &
        ( gammaf(kappa) * (1 + a*a * m*m)**( kappa+1.5 ) )

  end subroutine rmedia__psdf_vonKarman3d
  !! --------------------------------------------------------------------------------------------------------------------------- !!

  !! --------------------------------------------------------------------------------------------------------------------------- !!
  subroutine rmedia__psdf_ani_gauss3d ( mx, my, mz, ax, ay, az, epsil, psdf )

    !! -- Arguments
    real(SP), intent(in)  :: mx        ! wavenumber
    real(SP), intent(in)  :: my        ! wavenumber
    real(SP), intent(in)  :: mz        ! wavenumber
    real(SP), intent(in)  :: ax        ! correlation distance
    real(SP), intent(in)  :: ay        ! correlation distance
    real(SP), intent(in)  :: az        ! correlation distance
    real(SP), intent(in)  :: epsil     ! rms amplitude of fluctuation
    real(SP), intent(out) :: psdf      ! psdf value
    !! --

    psdf = (epsil*epsil) * (ax*ay*az) * PI**1.5 * exp ( -( (mx*ax)**2 +(my*ay)**2 + (mz*az)**2 ) * 0.25 )

  end subroutine rmedia__psdf_ani_gauss3d
  !! --------------------------------------------------------------------------------------------------------------------------- !!

  !! --------------------------------------------------------------------------------------------------------------------------- !!
  subroutine rmedia__psdf_ani_exp3d ( mx, my, mz, ax, ay, az, epsil, psdf )

    !! -- Arguments
    real(SP), intent(in)  :: mx        ! wavenumber
    real(SP), intent(in)  :: my        ! wavenumber
    real(SP), intent(in)  :: mz        ! wavenumber
    real(SP), intent(in)  :: ax        ! correlation distance
    real(SP), intent(in)  :: ay        ! correlation distance
    real(SP), intent(in)  :: az        ! correlation distance
    real(SP), intent(in)  :: epsil     ! rms amplitude of fluctuation
    real(SP), intent(out) :: psdf      ! psdf value
    real(SP) :: wk
    !! --

    wk = 1.0 + (mx*ax)**2 + (my*ay)**2 + (mz*az)**2
    psdf = ( 8 * epsil*epsil * ax*ay*az * PI )/( wk*wk )

  end subroutine rmedia__psdf_ani_exp3d
  !! --------------------------------------------------------------------------------------------------------------------------- !!

  !! --------------------------------------------------------------------------------------------------------------------------- !!
  subroutine rmedia__psdf_ani_vonKarman3d ( mx, my, mz, ax, ay, az, epsil, kappa, psdf )

    !! -- Arguments
    real(SP), intent(in)  :: mx        ! wavenumber
    real(SP), intent(in)  :: my        ! wavenumber
    real(SP), intent(in)  :: mz        ! wavenumber
    real(SP), intent(in)  :: ax        ! correlation distance
    real(SP), intent(in)  :: ay        ! correlation distance
    real(SP), intent(in)  :: az        ! correlation distance
    real(SP), intent(in)  :: epsil     ! rms amplitude of fluctuation
    real(SP), intent(in)  :: kappa     ! von karman parameter
    real(SP), intent(out) :: psdf      ! psdf value
    !! --

    psdf = ( 8 * PI_D**(1.5) * &
        epsil*epsil * ax*ay*az * gammaf( kappa + 1.5)  ) / &
        ( gammaf(kappa) * (1 + (mx*ax)**2 + (my*ay)**2 + (mz*az)**2 )**( kappa+1.5 ) )

  end subroutine rmedia__psdf_ani_vonKarman3d
  !! --------------------------------------------------------------------------------------------------------------------------- !!

end module m_rmedia
!! ----------------------------------------------------------------------------------------------------------------------------- !!
