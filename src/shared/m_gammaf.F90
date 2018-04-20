!! ------------------------------------------------------------------------------------------------------------------------------ !!
!>
!! Gamma function
!!
!! @copyright
!!   Copyright 2013-2018 Takuto Maeda. All rights reserved. This project is released under the MIT license.
!<
!! --
module m_gammaf

  use m_std
  implicit none
  private

  public :: gammaf

  real(DP) :: b(16)
  logical, save :: init

  !! generic interface
  interface gammaf
    module procedure gammaf_d, gammaf_s
  end interface gammaf

contains


  !! --------------------------------------------------------------------------------------------------------------------------- !!
  !>
  !! Gamma function, Double precision
  !<
  !!--
  real(DP) function gammaf_d(x)
    real(DP), intent(in) :: x
    if( x < 0 ) then
      gammaf_d = PI / ( sin(PI*x) * exp( gammaf_log( 1-x ) ) )
    else
      gammaf_d = exp( gammaf_log(x) )
    end if

  end function gammaf_d
  !! --------------------------------------------------------------------------------------------------------------------------- !!

  !! --------------------------------------------------------------------------------------------------------------------------- !!
  !>
  !! Gamma function, Single precision with conversion
  !<
  !!--
  real(SP) function gammaf_s(x)
    real(SP), intent(in) :: x

    gammaf_s = real( gammaf_d( dble(x) ) )

  end function gammaf_s
  !! --------------------------------------------------------------------------------------------------------------------------- !!



  !! --------------------------------------------------------------------------------------------------------------------------- !!
  !>
  !! Logarithm of the Gamma function
  !<
  !! --
  real(DP) function gammaf_log( x )

    real(DP), intent(in) :: x
    real(DP) :: v, y
    integer  :: i

    if( .not. init ) then
      b( 2) =     1.0_DP /    6.0_DP
      b( 4) =    -1.0_DP /   30.0_DP
      b( 6) =     1.0_DP /   42.0_DP
      b( 8) =    -1.0_DP /   30.0_DP
      b(10) =     5.0_DP /   66.0_DP
      b(12) =  -691.0_DP / 2730.0_DP
      b(14) =     7.0_DP /    6.0_DP
      b(16) = -3617.0_DP /  510.0_DP
      init = .true.
    end if

    y = x
    v = 1.0_DP
    do while( y < 8 )
      v = v * y
      y = y + 1.0_DP
    end do

    gammaf_log = 0
    do i=1, 8
      gammaf_log = gammaf_log + b(2*i) / ( 2*i * ( 2*i-1 ) * y**(2*i-1) )
    end do
    gammaf_log = gammaf_log + (y-0.5)*log(y) - y + 0.5*log(2*PI) - log(v)

  end function gammaf_log
  !! --------------------------------------------------------------------------------------------------------------------------- !!

end module m_gammaf
!! ----------------------------------------------------------------------------------------------------------------------------- !!
