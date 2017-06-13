!! ----------------------------------------------------------------------------------------------------------------------------- !!
!>
!! Calculate frequency dependense of Q^(-1) model
!!
!! @copyright
!!   Copyright 2013-2017 Takuto Maeda. All rights reserved. This project is released under the MIT license.
!<
!! --
program qmodel_tau

  use m_std
  use m_getopt
  use m_fdtool
  use m_readini
  implicit none

  integer :: i
  logical :: sw
  integer :: nm
  character(256) :: fn_prm
  real(SP), allocatable :: ts(:) ! tau_sigma
  real(SP) :: q0, f0, f1
  integer :: nf
  integer :: io
  integer :: ierr
  real(SP) :: fq_min, fq_max, fq_ref
  real(SP) :: zeta, tau
  real(SP) :: i1, i2
  integer :: l
  real(SP) :: f, omega
  real(SP) :: qinv
  real(SP), allocatable :: qinv2(:)
  character(20) :: fmt
  real(SP) :: chi, chi_R

  call getopt( 'nm', sw, nm    ); if( .not. sw ) call usage_exit()
  call getopt( 'i',  sw, fn_prm); if( .not. sw ) call usage_exit()
  call getopt( 'q0', sw, q0, 100.0 )
  call getopt( 'f0', sw, f0, 0.001 )
  call getopt( 'f1', sw, f1, 10.0 )
  call getopt( 'nf', sw, nf, 1000 )
  call std__getio( io )

  open( io, file=trim(fn_prm), action='read', iostat=ierr );
  if( ierr /= 0 ) then
    write(STDERR,*) 'input file open error'
    call usage_exit
  end if

  call readini( io, 'fq_min', fq_min, 0.01 )
  call readini( io, 'fq_max', fq_max, 5.00 )
  call readini( io, 'fq_ref', fq_ref, 1.00 )

  allocate( ts(nm) )
  call visco_set_relaxtime( nm, ts, fq_min, fq_max )
  zeta = visco_constq_zeta( nm, fq_min, fq_max, ts )
  tau = nm*zeta / q0
  allocate( qinv2(nm) )

  !! reference modulus
  chi_R = visco_chi( nm, ts, tau, fq_ref )

  write(fmt,'(A,I1.1,A)') 'F12.5,',nm+1,'ES12.5,ES12.5'
  do i=1, nf

    f = f0 * ( f1/f0 )**( real(i-1)/real(nf-1) )
    omega = 2*PI*f

    i1=0
    i2=0
    chi = 0
    do l=1, nm
      i1 = i1 + ( omega * ts(l) ) / ( 1 + omega**2*ts(l)**2 )
      i2 = i2 + (1 + omega**2 * ts(l)**2 * ( 1 + tau ) )/(1+omega**2*ts(l)**2 )
      qinv2(l) = ( omega * ts(l) * tau ) / ( 1 + omega**2*ts(l)**2 )/nm ! 要素Q
    end do
    qinv = tau * i1 / i2
    chi = visco_chi( nm, ts, tau, f )

    write(STDOUT,"("//trim(fmt)//")") f, qinv, qinv2(:), chi / chi_R
  end do


contains

  subroutine usage_exit()

    write(STDERR,*) 'qmodel_tau -nm [nm] -i [prm file] -f0 [min freq (0.001)] -f1 [max freq (10)] -nf [number of grid (1000)] '
    write(STDERR,*) ' [input parameter file (inf format)'
    write(STDERR,*) '   required field : fq_min (minimum freq for constant Q band) (0.05)'
    write(STDERR,*) '                    fq_max (maximum freq for constant Q band) (5.00)'
    write(STDERR,*) '                    fq_ref (reference frequency) (1.00)'

    stop
  end subroutine usage_exit



end program qmodel_tau
