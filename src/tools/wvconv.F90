!! ----------------------------------------------------------------------------------------------------------------------------- !!
!>
!! Convert waveforms of 2D simulation output into pseudo-3D
!!
!! Usage:
!!   wvconv.x in.sac out.sac (V0)
!!     in.sac: input SAC file (must be 2D simulation output)
!!     out.sac: output SAC filename
!!     V0 (optional): assumed average velocity in m/s unit. Default is 4000 m/s
!!
!! @copyright
!!   Copyright 2021 Takuto Maeda. All rights reserved. This project is released under the MIT license.
!<
!! --
program wvconv

    use m_std
    use m_sac
    use m_rfft
    use m_fk
    implicit none

    character(256) :: fn_in, fn_out
    real(dp), allocatable :: dat_sac(:), dat_fft(:), freq(:)
    complex(dp), allocatable :: spec(:)
    type(sac__hdr) :: sh
    integer :: nfft
    integer :: i
    real(dp) :: R
    real(dp), parameter :: Vdef = 4000. ! average velocity ( P or S )
    real(dp) :: V
    complex(dp), parameter  :: zi = ( 0.0_dp, 1.0_dp )
    integer :: p
    character(10) :: adum
    !---

    if( command_argument_count() == 0 ) then
      write(STDERR,*) 'wvconv.x in.sac out.sac (V0)'
      stop
    end if
    

    call get_command_argument(1, fn_in)
    call get_command_argument(2, fn_out)
    if( command_argument_count() == 3 ) then
      call get_command_argument(3, adum)
      read(adum,*) V
    else
      V = Vdef
    end if
  
    call sac__read(fn_in, sh, dat_sac)

    ! fft size
    p = ceiling( log(dble(sh%npts) ) / log(2.0_dp) ) 
    nfft = 2 ** p 
    
    allocate(dat_fft(nfft))
    allocate(spec(nfft/2+1), freq(nfft/2+1))
  
    dat_fft = 0.0
    dat_fft(1:sh%npts) = dat_sac(:)
    call fk__t2f( nfft, sh%delta, dat_fft, spec, freq )

    do i=1, nfft/2+1
      spec(i)= spec(i) * sqrt( - zi * freq(i) ) 
    end do

    call fk__f2t(nfft, sh%delta, spec, dat_fft)

    R = sqrt( sh%dist**2 + sh%evdp**2 ) * 1000. ! in m-unit 
    dat_sac(1:sh%npts) = dat_fft(1:sh%npts) / sqrt( R * V)
    sh%kevnm = trim(adjustl(sh%kevnm)) // '.CONV'
    call sac__write(fn_out, sh, dat_sac, .true. )

end program wvconv
!! ----------------------------------------------------------------------------------------------------------------------------- !!
