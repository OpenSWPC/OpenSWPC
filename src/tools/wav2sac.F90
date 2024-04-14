!! --------------------------------------------------------------------------------------------- !!
!>
!! Convert waveform binary (.wav) to SAC datafiles
!!
!! @copyright
!!   Copyright 2013-2024 Takuto Maeda. All rights reserved. This project is released under the MIT license.
!<
!! ----
program wav2sac

  use m_std
  use m_sac
  use m_system
  use m_getopt
  use m_version
  implicit none

  character(256) :: fn_wav, fn_sac
  type(sac__hdr), allocatable :: sh(:)
  integer :: ierr
  real(SP), allocatable :: wav(:,:)
  integer :: nch, npts
  integer :: nfile
  integer :: io
  integer :: i, j
  character(80) :: title
  logical :: is_opt1, is_opt2
  !! ----

  call getopt('v', is_opt1)
  call getopt('-version', is_opt2)
  if( is_opt1 .or. is_opt2 ) call version__display('wav2sac')


  if( system__iargc() == 0 ) then
    write(*,*) 'wav2sac.x [wavfiles]'
    stop
  end if

  nfile = system__iargc()

  do i=1, nfile
    call system__getarg(i, fn_wav)
    call std__getio(io, is_big = .true. )
    open( io, file=trim(fn_wav), access='stream', action='read', status='old', iostat=ierr )
    if( ierr /= 0 ) cycle

    write(STDERR,*) "FILE: ", trim(fn_wav)

    do
      read( io, iostat=ierr ) nch, npts, title
      if( ierr /= 0 ) then
        close(io)
        exit
      end if
      allocate( sh(nch), wav(npts,nch) )

      read(io) sh
      read(io) wav

      do j=1, nch
        if( trim(title) == trim(sh(j)%kevnm) ) then
          fn_sac = trim(title) // '__' // trim( sh(j)%kstnm ) // '__' // &
              trim(sh(j)%kcmpnm) // '__.sac'
        else
          fn_sac = trim(title) // '__' // trim(sh(j)%kevnm) // '__' // &
              trim( sh(j)%kstnm ) // '__' // trim(sh(j)%kcmpnm) // '__.sac'
        end if
        write(STDERR,*) "    . ",  trim(fn_sac)
        call sac__write( fn_sac, sh(j), wav(:,j), .true. )
      end do
      deallocate( sh, wav )

    end do

    close(io)

  end do

end program wav2sac
!! --------------------------------------------------------------------------------------------- !!
