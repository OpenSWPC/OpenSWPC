!! ----------------------------------------------------------------------------------------------------------------------------- !!
!>
!! FDM model parameter condition check
!!
!! @copyright
!!   Copyright 2013-2017 Takuto Maeda. All rights reserved. This project is released under the MIT license.
!<
!! -----
program fdmcond

  implicit none
  real :: ef
  integer :: ng = 7
  integer :: dim
  integer :: stype
  integer :: inq
  real :: dt, vmin, vmax, tr, fmax
  real, allocatable :: dh(:)

  call term_cls()
  write(*,*)
  write(*,'(A)') "----------------------------------------------------------------------"
  write(*,'(A)') "                           FDM CONDITION                              "
  write(*,'(A)') "----------------------------------------------------------------------"

  call fdline(3)
  write(*,'(A)') "   2) 2D"
  write(*,'(A)') "   3) 3D"
  call bkline(3)
  write(*,'(A)', advance="no") " Model Dimension ? --> "
  read(*,*) dim

  if( dim /= 2 .and. dim /= 3 ) then
    write(*,*) "Invalid Dimension. Quit."
    stop
  end if

  allocate( dh(dim) )

  call fdline(5)
  write(*,'(A)') "   1) Triangle"
  write(*,'(A)') "   2) Herrmann"
  write(*,'(A)') "   3) Kupper"
  call bkline(4)
  write(*,'(A)', advance='no') " Source Type ? --> "

  read(*,*) stype

  if( stype == 1 ) then
    ef = 5.4
  else if ( stype == 2 ) then
    ef = 3.2
  else if ( stype == 3 ) then
    ef = 2.3
  else
    write(*,*) "Invalid Source Type. Quit."
    stop
  end if

  call fdline(6)

  write(*,'(A)') "   1) dh   (space grid),  fmax (max freq.),  vmax (max vel.)"
  write(*,'(A)') "   2) dh   (space grid),  Tr   (rise time),  vmax (max vel.)"
  write(*,'(A)') "   3) dh   (space grid),  fmax (max freq.),  dt (time grid)"
  write(*,'(A)') "   4) dh   (space grid),  Tr   (rise time),  dt (time grid)"
  write(*,'(A)') "   5) dh   (space grid),  vmin (min vel.),   vmax (max vel.)"
  write(*,'(A)') "   6) dh   (space grid),  vmin (min vel.),   dt (time grid) "
  write(*,'(A)') "   7) fmax (max freq.) ,  vmax (max vel.),   dt (time grid)"
  write(*,'(A)') "   8) Tr   (rise time) ,  vmax (max vel.),   dt (time grid)"
  write(*,'(A)') "   9) vmin (min vel.)  ,  vmax (max vel.),   dt (time grid)"

  call bkline(10)
  write(*,'(A)', advance='no') " Parameter Combination ? --> "
  read(*,*) inq
  call fdline(11)

  write(*,'(A)') " Assumed Parameters: "
  select case( inq )
  case(1)
    call get_dh ( dh )
    fmax = getv( 'fmax' )
    vmax = getv( 'vmax' )
    call cond_1( dh, fmax, vmax )

  case(2)
    call get_dh ( dh )
    tr = getv( 'Tr' )
    vmax = getv( 'vmax' )
    call cond_2( dh, Tr, vmax )

  case(3)
    call get_dh ( dh )
    fmax = getv( 'fmax' )
    dt   = getv( 'dt' )
    call cond_3( dh, fmax, dt )

  case(4)
    call get_dh ( dh )
    tr  = getv( 'Tr' )
    dt  = getv( 'dt' )
    call cond_4( dh, Tr, dt )

  case(5)
    call get_dh ( dh )
    vmin = getv( 'vmin' )
    vmax = getv( 'vmax' )
    call cond_5( dh, vmin, vmax )

  case(6)
    call get_dh ( dh )
    vmin = getv( 'vmin' )
    dt   = getv( 'dt'  )
    call cond_6( dh, vmin, dt )

  case(7)
    fmax = getv( 'fmax' )
    vmax = getv( 'vmax' )
    dt   = getv( 'dt' )
    call cond_7( fmax, vmax, dt )
  case(8)
    tr   = getv( 'Tr' )
    vmax = getv( 'vmax' )
    dt   = getv( 'dt' )
    call cond_8( tr, vmax, dt )
  case(9)
    vmin = getv( 'vmin' )
    vmax = getv( 'vmax' )
    dt   = getv( 'dt' )
    call cond_9( vmin, vmax, dt )
  end select

  write(*,*)

contains

  real function getv( nm  )
    character(*), intent(in) :: nm
    integer :: nc
    character(5) :: query1
    nc = len_trim(nm)
    query1 = '     '
    query1(1:nc) = trim(adjustl(nm))

    write(*,'(A)',advance='no') '   ' // query1 // '  =   '
    read(*,*) getv
  end function getv


  subroutine fdline(n)
    integer, intent(in) :: n
    integer :: i
    do i=1, n
      write(*,*)
    end do
  end subroutine fdline


  subroutine bkline(n)
    integer, intent(in) :: n
    integer :: i
    do i=1, n
      write(*,'(1X,A)') char(27)//'[2A'
    end do
  end subroutine bkline

  subroutine cond_1( dh, fmax, vmax )
    real, intent(in) :: dh(:), fmax, vmax
    real :: dt, tr, vmin

    dt = 6. / 7. / ( vmax * sqrt( sum( 1/dh(:)**2 )  ) )
    tr = ef / fmax
    vmin = fmax * ng * maxval(dh)

    write(*,*)
    write(*,'(A)') " Derivaed Parameters: "
    write(*,'(A,F9.5)') "   Tr     = ", tr
    write(*,'(A,F9.5)') "   dt    <= ", dt
    write(*,'(A,F9.5)') "   vmin  >= ", vmin

  end subroutine cond_1


  subroutine cond_2( dh, Tr, vmax )

    real, intent(in) :: dh(:), Tr, vmax
    real :: dt, fmax, vmin

    dt = 6. / 7. / ( vmax * sqrt( sum( 1/dh(:)**2 ) ) )
    fmax = ef / tr
    vmin = fmax * ng * maxval(dh)
    write(*,*)
    write(*,'(A)') " Derivaed Parameters: "
    write(*,'(A,F9.5)') "   fmax   = ", fmax
    write(*,'(A,F9.5)') "   dt    <= ", dt
    write(*,'(A,F9.5)') "   vmin  >= ", vmin

  end subroutine cond_2


  subroutine cond_3( dh, fmax, dt )
    real, intent(in) :: dh(:), fmax, dt
    real :: vmax, tr, vmin

    vmax = 6. / 7. / ( dt * sqrt( sum( 1/dh(:)**2 ) ) )
    tr = ef / fmax
    vmin = fmax * ng * maxval(dh)
    write(*,*)
    write(*,'(A)') " Derivaed Parameters: "
    write(*,'(A,F9.5)') "   Tr     = ", Tr
    write(*,'(A,F9.5)') "   vmax  <= ", vmax
    write(*,'(A,F9.5)') "   vmin  >= ", vmin

  end subroutine cond_3


  subroutine cond_4( dh, tr, dt )
    real, intent(in) :: dh(:), tr, dt
    real :: vmax, fmax, vmin

    vmax = 6. / 7. / ( dt * sqrt( sum( 1/dh(:)**2 ) ) )
    fmax = ef / tr
    vmin = fmax * ng * maxval(dh)
    write(*,*)
    write(*,'(A)') " Derivaed Parameters: "
    write(*,'(A,F9.5)') "   fmax    = ", fmax
    write(*,'(A,F9.5)') "   vmax   <= ", vmax
    write(*,'(A,F9.5)') "   vmin   >= ", vmin

  end subroutine cond_4

  subroutine cond_5( dh, vmin, vmax )

    real, intent(in) :: dh(:), vmin, vmax
    real :: dt, fmax, tr

    dt = 6. / 7. / ( vmax * sqrt( sum( 1/dh(:)**2 ) ) )
    fmax = vmin / ( ng * maxval(dh) )
    tr   = ef / fmax

    write(*,*)
    write(*,'(A)') " Derivaed Parameters: "
    write(*,'(A,F9.5)') "   dt    <= ", dt
    write(*,'(A,F9.5)') "   fmax  <= ", fmax
    write(*,'(A,F9.5)') "   Tr    >= ", tr

  end subroutine cond_5

  subroutine cond_6( dh, vmin, dt )
    real, intent(in) :: dh(:), vmin, dt
    real :: vmax, fmax, tr

    vmax = 6. / 7. / ( dt * sqrt( sum( 1/dh(:)**2 ) ) )
    fmax = vmin / ( ng * maxval(dh) )
    tr   = ef / fmax

    write(*,*)
    write(*,'(A)') " Derivaed Parameters: "
    write(*,'(A,F9.5)') "   vmax  <= ", vmax
    write(*,'(A,F9.5)') "   fmax  <= ", fmax
    write(*,'(A,F9.5)') "   Tr    >= ", tr

  end subroutine cond_6


  subroutine cond_7 ( fmax, vmax, dt )
    real, intent(in) :: fmax, vmax, dt
    real :: dh, tr, vmin

    dh = sqrt(real(dim)) * dt * vmax * 7. / 6.
    tr = ef / fmax
    vmin = fmax * ng * dh

    write(*,*)
    write(*,'(A)') " Derivaed Parameters: "
    write(*,'(A,F9.5)') "   Tr     = ", tr
    write(*,'(A,F9.5)') "   dh    <= ", dh
    write(*,'(A,F9.5)') "   vmin  >= ", vmin
    write(*,*)
    write(*,*) "dx=dy=dz=dh has been assumed. vmin value depends dh. "
    write(*,*) "Please decide dx/dy/dz, and retry by giving (dx/dy/dz, dt, tr) or (dx/dy/dz, dt, fmax)"

  end subroutine cond_7

  subroutine cond_8( tr, vmax, dt )
    real, intent(in) :: tr, vmax, dt
    real :: dh, fmax, vmin


    dh = sqrt(real(dim)) * dt * vmax * 7. / 6.
    fmax = ef / tr
    vmin = fmax * ng * dh

    write(*,*)
    write(*,'(A)') " Derivaed Parameters: "
    write(*,'(A,F9.5)') "   fmax   = ", fmax
    write(*,'(A,F9.5)') "   dx    >= ", dh
    if( dim == 3 ) then
      write(*,'(A,F9.5)') "   dy    >= ", dh
    end if
    write(*,'(A,F9.5)') "   dz    >= ", dh
    write(*,'(A,F9.5)') "   vmin  >= ", vmin
    write(*,*)
    if( dim==3 ) write(*,*) "dx=dy=dz=dh has been assumed. vmin depends choice of dh. "
    if( dim==2 ) write(*,*) "dx=dz=dh has been assumed. vmin depends choice of dh. "

  end subroutine cond_8

  subroutine cond_9( vmin, vmax, dt )
    real, intent(in) :: vmin, vmax, dt
    real :: tr, fmax, dh


    dh = sqrt(real(dim)) * dt * vmax * 7. / 6.
    fmax = vmin / ( ng * dh )
    tr = 1 / fmax

    write(*,*)
    write(*,'(A)') " Derivaed Parameters: "
    write(*,'(A,F9.5)') "   dx    >= ", dh
    if( dim == 3 ) then
      write(*,'(A,F9.5)') "   dy    >= ", dh
    end if
    write(*,'(A,F9.5)') "   dz    >= ", dh
    write(*,'(A,F9.5)') "   fmax  <= ", fmax
    write(*,'(A,F9.5)') "   tr    >= ", tr
    write(*,*)
    if( dim==3 ) write(*,*) "Assumed dx=dy=dz. vmin depends choice of dh. "
    if( dim==2 ) write(*,*) "Assumed dx=dz. vmin depends choice of dh. "

  end subroutine cond_9

  subroutine get_dh( dh )
    real, intent(out) :: dh(:)

    if( dim == 2 ) then
      dh(1) = getv( 'dx' )
      dh(2) = getv( 'dz')
    else
      dh(1) = getv( 'dx')
      dh(2) = getv( 'dy')
      dh(3) = getv( 'dz')
    end if
  end subroutine get_dh

  subroutine term_cls

    write(*,'(1X,A)') char(27)//'[2J'
    write(*,'(1X,A)') char(27)//'[500A' ! let the cursor the top of the screen

  end subroutine term_cls


end program fdmcond
