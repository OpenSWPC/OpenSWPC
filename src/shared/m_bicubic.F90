!! ----------------------------------------------------------------------------------------------------------------------------- !!
!>
!! Bicubic interpolation
!!
!! @copyright
!!   Copyright 2013-2020 Takuto Maeda. All rights reserved. This project is released under the MIT license.
!<
!! ----------------------------------------------------------------------------------------------------------------------------- !!
module m_bicubic

  !! Dependency
  use m_std
  implicit none
  private
  save

  type bicubic__data
    real(DP), allocatable :: f(:,:), fx(:,:), fy(:,:), fxy(:,:)
    real(DP) :: x0, y0, dx, dy
    integer  :: nx, ny
    logical  :: is_init  = .false.
    logical  :: is_first = .true.
    integer  :: ii0, jj0
    real(DP) :: aa(0:3,0:3)
  end type bicubic__data

  logical :: is_global_first = .true.

  !! Public Procedures
  public :: bicubic__init
  public :: bicubic__terminate
  public :: bicubic__interp
  public :: bicubic__data


  !! --------------------------------------------------------------------------------------------------------------------------- !!
  interface bicubic__init

    module procedure   bicubic__init_s, bicubic__init_d

  end interface bicubic__init
  !! --------------------------------------------------------------------------------------------------------------------------- !!


  !! --------------------------------------------------------------------------------------------------------------------------- !!
  interface bicubic__interp

    module procedure   bicubic__interp_s, bicubic__interp_d

  end interface bicubic__interp
  !! --------------------------------------------------------------------------------------------------------------------------- !!

contains

  !! --------------------------------------------------------------------------------------------------------------------------- !!
  !>
  !! Set-up new inputdata for bicubic interpolation
  !<
  !! --
  subroutine bicubic__init_s( bd, nxi, nyi, x0i, y0i, dxi, dyi, dati )

    !! Arguments
    type(bicubic__data), intent(inout) :: bd
    integer,             intent(in)    :: nxi, nyi
    real(SP),            intent(in)    :: x0i, y0i
    real(SP),            intent(in)    :: dxi, dyi
    real(SP),            intent(in)    :: dati(nxi,nyi)
    real(DP),            allocatable   :: ddat(:,:)
    !! ----

    allocate( ddat(nxi,nyi) )
    ddat = dble(dati)
    call bicubic__init_d( bd, nxi, nyi, dble(x0i), dble(y0i), dble(dxi), dble(dyi), ddat )
    deallocate( ddat )

  end subroutine bicubic__init_s
  !! --------------------------------------------------------------------------------------------------------------------------- !!


  !! --------------------------------------------------------------------------------------------------------------------------- !!
  !>
  !! Set-up new inputdata for bicubic interpolation
  !<
  !! --
  subroutine bicubic__init_d( bd, nxi, nyi, x0i, y0i, dxi, dyi, dati )

    !! Arguments
    type(bicubic__data), intent(inout) :: bd
    integer,             intent(in)    :: nxi, nyi
    real(DP),            intent(in)    :: x0i, y0i
    real(DP),            intent(in)    :: dxi, dyi
    real(DP),            intent(in)    :: dati(nxi,nyi)
    !! --

    ! deallocate memory for previous dataset
    if( allocated( bd%f ) .or. bd%is_init )  call bicubic__terminate( bd )

    ! private save data
    bd%nx = nxi
    bd%ny = nyi
    bd%x0 = x0i
    bd%y0 = y0i
    bd%dx = dxi
    bd%dy = dyi

    ! memory allocation
    allocate( bd%f(bd%nx,bd%ny), bd%fx(bd%nx,bd%ny), bd%fy(bd%nx,bd%ny), bd%fxy(bd%nx,bd%ny) )

    ! data copy
    bd%f(1:bd%nx,1:bd%ny) = dati(1:bd%nx,1:bd%ny)

    ! pertial derivatives
    call diffx( bd%nx, bd%ny, bd%dx, bd%f,  bd%fx )
    call diffy( bd%nx, bd%ny, bd%dy, bd%f,  bd%fy )
    call diffx( bd%nx, bd%ny, bd%dx, bd%fy, bd%fxy )

    ! scaling
    bd%fx = bd%fx * bd%dx
    bd%fy = bd%fy * bd%dy
    bd%fxy = bd%fxy * bd%dx * bd%dy

    bd%is_init = .true.
    bd%is_first = .true.

  end subroutine bicubic__init_d
  !! --------------------------------------------------------------------------------------------------------------------------- !!

  !! --------------------------------------------------------------------------------------------------------------------------- !!
  !>
  !! Terminate using bicubic module
  !<
  !! --
  subroutine bicubic__terminate( bd )

    type(bicubic__data), intent(inout) :: bd

    !! --

    deallocate( bd%f, bd%fx, bd%fy, bd%fxy )
    bd%is_init = .false.
    bd%is_first = .true.

  end subroutine bicubic__terminate
  !! --------------------------------------------------------------------------------------------------------------------------- !!

  !! --------------------------------------------------------------------------------------------------------------------------- !!
  !>
  !! Get interporated data at xi and yi: single-precision interface
  !<
  !! --
  subroutine bicubic__interp_s( bd, xi, yi, v, default )

    !! - Arguments
    type(bicubic__data), intent(inout) :: bd
    real(SP), intent(in) :: xi
    real(SP), intent(in) :: yi
    real(SP), intent(out) :: v
    real(SP), optional, intent(in) :: default

    real(DP) :: vv

    !! --

    if(present(default)) then
      call bicubic__interp_d( bd, dble(xi), dble(yi), vv, dble(default) )
    else
      call bicubic__interp_d( bd, dble(xi), dble(yi), vv )
    end if

    v = real( vv )

  end subroutine bicubic__interp_s
  !! --------------------------------------------------------------------------------------------------------------------------- !!


  !! --------------------------------------------------------------------------------------------------------------------------- !!
  !>
  !! Get interporated data at xi and yi
  !<
  !! --
  subroutine bicubic__interp_d( bd, xi, yi, v, default )

    !! - Arguments
    type(bicubic__data), intent(inout) :: bd
    real(DP),            intent(in)    :: xi
    real(DP),            intent(in)    :: yi
    real(DP),            intent(out)   :: v
    real(DP), optional,  intent(in)    :: default

    !! --
    real(DP) :: xd, yd
    integer  :: ii, jj
    integer  :: i, j
    real(DP) :: xi2, yi2
    real(DP) :: x0, y0, dx, dy
    integer  :: nx, ny
    real(DP) :: xda(0:3), yda(0:3)

    !! ----

    !! local copy
    x0 = bd%x0
    y0 = bd%y0
    dx = bd%dx
    dy = bd%dy
    nx = bd%nx
    ny = bd%ny

    ii = floor(  ( xi - x0 ) / dx  ) + 1
    jj = floor(  ( yi - y0 ) / dy  ) + 1


    !! outlier
    xi2 = xi
    yi2 = yi

    if( ii < 1 .or. ii > nx-1 .or. jj <= 0 .or. jj > ny-1 ) then
      if( present( default ) ) then
        v = default
        return
      else

        if( ii < 1   ) then
          ii =   1
          xi2 = x0
        end if
        if( ii > nx-1 ) then
          ii = nx-1
          xi2 = x0+(nx-1)*Dx
        end if
        if( jj < 1   ) then
          jj =    1
          yi2 = y0
        end if
        if( jj > ny-1 ) then
          jj = ny-1
          yi2 = y0+(ny-1)*Dy
        end if
      end if
    end if

    if( ( bd%ii0 /= ii ) .or. ( bd%jj0 /= jj ) .or. bd%is_first ) then
      call bicubic__coef( bd, ii, jj )
    end if

    ! in-pixel distance
    xd = ( xi2 - ( x0 + (ii-1)*dx ) ) / dx
    yd = ( yi2 - ( y0 + (jj-1)*dy ) ) / dy

    xda(0) = 1.0_DP
    xda(1) = xd
    xda(2) = xd*xd
    xda(3) = xd*xd*xd

    yda(0) = 1.0_DP
    yda(1) = yd
    yda(2) = yd*yd
    yda(3) = yd*yd*yd

    v = 0.0_DP
    do j=0, 3
      do i=0, 3
        v = v + bd%aa(i,j)*xda(i)*yda(j)
      end do
    end do

    bd%is_first = .false.

    !! store ii, jj value for next subroutine call
    bd%ii0 = ii
    bd%jj0 = jj

  end subroutine bicubic__interp_d
  !! --------------------------------------------------------------------------------------------------------------------------- !!

  !! --------------------------------------------------------------------------------------------------------------------------- !!
  !>
  !! Set-up interpolation table
  !<
  !! --
  subroutine bicubic__coef( bd, ii, jj )

    !! -- Arguments
    type(bicubic__data), intent(inout) :: bd
    integer,             intent(in)    :: ii, jj ! pixel index

    real(DP) :: xx(16)
    real(DP), save :: mat(16,16)
    !! ----


    if( is_global_first ) then
      mat( 1,:) = (/  1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 /)
      mat( 2,:) = (/  0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 /)
      mat( 3,:) = (/ -3, 3, 0, 0,-2,-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 /)
      mat( 4,:) = (/  2,-2, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 /)
      mat( 5,:) = (/  0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0 /)
      mat( 6,:) = (/  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0 /)
      mat( 7,:) = (/  0, 0, 0, 0, 0, 0, 0, 0,-3, 3, 0, 0,-2,-1, 0, 0 /)
      mat( 8,:) = (/  0, 0, 0, 0, 0, 0, 0, 0, 2,-2, 0, 0, 1, 1, 0, 0 /)
      mat( 9,:) = (/ -3, 0, 3, 0, 0, 0, 0, 0,-2, 0,-1, 0, 0, 0, 0, 0 /)
      mat(10,:) = (/  0, 0, 0, 0,-3, 0, 3, 0, 0, 0, 0, 0,-2, 0,-1, 0 /)
      mat(11,:) = (/  9,-9,-9, 9, 6, 3,-6,-3, 6,-6, 3,-3, 4, 2, 2, 1 /)
      mat(12,:) = (/ -6, 6, 6,-6,-3,-3, 3, 3,-4, 4,-2, 2,-2,-2,-1,-1 /)
      mat(13,:) = (/  2, 0,-2, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0 /)
      mat(14,:) = (/  0, 0, 0, 0, 2, 0,-2, 0, 0, 0, 0, 0, 1, 0, 1, 0 /)
      mat(15,:) = (/ -6, 6, 6,-6,-4,-2, 4, 2,-3, 3,-3, 3,-2,-1,-2,-1 /)
      mat(16,:) = (/  4,-4,-4, 4, 2, 2,-2,-2, 2,-2, 2,-2, 1, 1, 1, 1 /)
      is_global_first = .false.
    end if


    xx( 1) = bd%f  (ii,  jj)
    xx( 2) = bd%f  (ii+1,jj  )
    xx( 3) = bd%f  (ii,  jj+1)
    xx( 4) = bd%f  (ii+1,jj+1)

    xx( 5) = bd%fx (ii,  jj  )
    xx( 6) = bd%fx (ii+1,jj  )
    xx( 7) = bd%fx (ii,  jj+1)
    xx( 8) = bd%fx (ii+1,jj+1)

    xx( 9) = bd%fy (ii,  jj  )
    xx(10) = bd%fy (ii+1,jj  )
    xx(11) = bd%fy (ii,  jj+1)
    xx(12) = bd%fy (ii+1,jj+1)

    xx(13) = bd%fxy(ii,  jj  )
    xx(14) = bd%fxy(ii+1,jj  )
    xx(15) = bd%fxy(ii,  jj+1)
    xx(16) = bd%fxy(ii+1,jj+1)


    bd%aa(0:3,0) = matmul( mat( 1: 4,:), xx )
    bd%aa(0:3,1) = matmul( mat( 5: 8,:), xx )
    bd%aa(0:3,2) = matmul( mat( 9:12,:), xx )
    bd%aa(0:3,3) = matmul( mat(13:16,:), xx )


  end subroutine bicubic__coef
  !! --------------------------------------------------------------------------------------------------------------------------- !!


  !! =========================================================================================================================== !!
  !!                                                        private routines
  !! =========================================================================================================================== !!

  !! --------------------------------------------------------------------------------------------------------------------------- !!
  !>
  !! Second-order derivative
  !<
  !! --
  subroutine diffx( nx, ny, dx, ff, ffx )

    integer,  intent(in)  :: nx, ny
    real(DP), intent(in)  :: dx
    real(DP), intent(in)  :: ff(nx,ny)
    real(DP), intent(out) :: ffx(nx,ny)
    !! --
    integer :: i, j
    !! --

    do j=1, ny
      do i=2, nx-1
        ffx(i,j) = ( ff(i+1,j) - ff(i-1,j) ) / ( 2 * dx )
      end do

      ! approximate by one-sided derivatives around the boudnary
      ffx(1,j) = ( ff(2,j) - ff(1,j) ) / dx
      ffx(nx,j) = ( ff(nx,j) - ff(nx-1,j) ) / dx
    end do

  end subroutine diffx
  !! --------------------------------------------------------------------------------------------------------------------------- !!

  !! --------------------------------------------------------------------------------------------------------------------------- !!
  !>
  !! Second-order derivative
  !<
  !! --
  subroutine diffy( nx, ny, dy, ff, ffy )

    integer,  intent(in)  :: nx, ny
    real(DP), intent(in)  :: dy
    real(DP), intent(in)  :: ff(nx,ny)
    real(DP), intent(out) :: ffy(nx,ny)
    !! --
    integer :: i, j
    !! --

    do i=1, nx
      do j=2, ny-1
        ffy(i,j) = ( ff(i,j+1) - ff(i,j-1) ) / ( 2 * dy )
      end do
      !  approximate by one-sided derivatives
      ffy(i,1) = ( ff(i,2) - ff(i,1) ) / dy
      ffy(i,ny) = ( ff(i,ny) - ff(i,ny-1) ) / dy
    end do
  end subroutine diffy
  !! --------------------------------------------------------------------------------------------------------------------------- !!


end module m_bicubic
!! ----------------------------------------------------------------------------------------------------------------------------- !!
