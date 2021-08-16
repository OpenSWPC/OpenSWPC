module m_seawater
  
  use m_std
  implicit none
  private
  public :: seawater__vel, seawater__init

  real(dp), parameter :: epsil0 = 0.00737_dp
  real(dp), save :: epsil
  logical, save :: init = .false. 
contains

  subroutine seawater__init( use_munk )
    logical, intent(in) :: use_munk
    
    if( use_munk ) then
      epsil = epsil0
    else
      epsil = 0.0_dp  !< force seawater velocity 1500 m/s irrespective to depth
    end if

    init = .true.
  end subroutine seawater__init

  function seawater__vel(z) result(c)

    real :: c
    real, intent(in) :: z
    real(dp) :: zb, zc

    if( .not. init ) call seawater__init(.false.)

    zc = 1300.0_dp
    zb = 2 * ( dble(z) * 1000. - zc ) / zc

    c = real(1.5_dp * ( 1.0_dp + epsil * ( zb - 1.0_dp + exp( - zb )) ) )

  end function seawater__vel

end module m_seawater
