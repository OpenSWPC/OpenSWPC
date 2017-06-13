!! ----------------------------------------------------------------------------------------------------------------------------- !!
!>
!! Green's Function Special Mode for SWPC
!!
!! @copyright
!!   Copyright 2013-2017 Takuto Maeda. All rights reserved. This project is released under the MIT license.
!<
!! ----
#include "m_debug.h"
module m_green

  !! -- Declarations
  use m_std
  use m_debug
  use m_global
  use m_wsac
  use m_output
  use m_readini
  use m_geomap
  use m_pwatch
  use m_system
  use mpi
  implicit none
  private
  save

  !! -- Public Procedures
  public :: green__setup
  public :: green__source
  public :: green__store
  public :: green__export

  !!
  !! Green's Function Mode
  !!
  logical :: is_src = .false. !< true if the pseudo source (station loc) is contained in the node

  !! unit conversion coefficient
  real(SP), parameter :: UC_BF    = 10.0**(-12)
  real(SP), parameter :: UC_DERIV = 10.0**(-15)
  real(SP) :: dt_dxyz
  integer  :: isrc, jsrc, ksrc
  real(SP) :: xsrc, ysrc, zsrc, evlo0, evla0
  real(SP) :: fx1, fy1, fz1
  character(3) :: cmp(9)
  character(3) :: wav_format

  integer :: ntw
  integer :: ng

  !! -- Parameters
  real(MP), parameter   :: C20 = 1.0_MP
  real(MP), parameter   :: C40 = 9.0_MP /  8.0_MP
  real(MP), parameter   :: C41 = 1.0_MP / 24.0_MP
  real(MP) :: r40x, r40y, r40z, r41x, r41y, r41z, r20x, r20y, r20z

  !! reuse from m_output.F90
  integer :: ntdec_w

  integer, allocatable :: ig(:), jg(:), kg(:)
  real(SP), allocatable :: xg(:), yg(:), zg(:), long(:), latg(:)
  integer, allocatable :: gid(:)
  type(sac__hdr), allocatable :: sh(:)
  real(SP), allocatable :: gf(:,:)
  character(256), allocatable :: fn(:)
  character(256) :: fn_csf

  real(SP), allocatable :: Ux(:), Uy(:), Uz(:)
  real(SP), allocatable :: dxUx(:), dyUx(:), dzUx(:), dxUy(:), dyUy(:), dzUy(:), dxUz(:), dyUz(:), dzUz(:)
  real(SP), parameter :: green_tbeg = 0.0
  real(SP) :: green_trise
  character :: green_cmp
  character(16) :: stftype
  logical :: green_pbmode
  character(3) :: green_plate
  logical :: green_bforce
  integer :: ncmp
contains


  !! --------------------------------------------------------------------------------------------------------------------------- !!
  !>
  !! Set up Green's function mode
  !!
  !! it must be called after output__setup
  !<
  !! --
  subroutine green__setup(io_prm)

    integer, intent(in) :: io_prm
    !! --
    logical        :: src_somewhere
    integer        :: ierr
    character(8)   :: green_stnm
    integer        :: i, j, k
    real(SP)       :: xg0, yg0, zg0, long0, latg0
    integer        :: ii, jj, kk
    integer        :: gid0
    character(8)   :: cid8
    character(256) :: fn_glst
    integer        :: io
    logical        :: is_exist
    character(3)   :: green_fmt
    real(SP)       :: green_maxdist
    integer        :: ng0
    real(SP)       :: dd
    character(256) :: abuf
    real(SP)       :: evlo1, evla1
    !! ----

    if( benchmark_mode ) then
      green_mode = .false.
      return
    end if

    call pwatch__on( 'green__setup' )

    call readini( io_prm, 'green_mode', green_mode, .false. )


    if( .not. green_mode ) then
      call pwatch__off( 'green__setup' )
      return
    end if

    !! FDM coefficients
    r40x = C40 / dx
    r40y = C40 / dy
    r40z = C40 / dz
    r41x = C41 / dx
    r41y = C41 / dy
    r41z = C41 / dz
    r20x = 1.  / dx
    r20y = 1.  / dy
    r20z = 1.  / dz

    !!
    !! parameter input
    !!

    !! pseudo source location
    call readini( io_prm, 'green_stnm', green_stnm, '' )
    call readini( io_prm, 'green_cmp',  green_cmp,  '' )
    call readini( io_prm, 'green_trise', green_trise, 1.0 )
    call readini( io_prm, 'green_bforce', green_bforce, .false. )
    call readini( io_prm, 'green_maxdist', green_maxdist, 1e30 )    
    call assert( green_maxdist > 0.0 )

    M0 = 1
    fmax = 2.0 / green_trise

    !! green's function output grid layout
    call readini( io_prm, 'fn_glst', fn_glst, '' )
    inquire( file=fn_glst, exist=is_exist )
    call assert( is_exist )

    call readini( io_prm, 'green_fmt', green_fmt, 'xyz' )
    call assert( green_fmt == 'xyz' .or. green_fmt == 'llz' )

    call readini( io_prm, 'ntdec_w',ntdec_w, 10 )
    call readini( io_prm, 'stftype', stftype, 'kupper' )
    call readini( io_prm, 'wav_format', wav_format, 'sac' ) 

    if( trim(adjustl(stftype)) == 'scosine' ) stftype = 'cosine'  !! backward compatibility

    !!
    !! Set up pseudo source, information is obtained from m_output
    !!
    select case (green_cmp)
    case( 'x' );   fx1 = 1.0; fy1 = 0.0; fz1 = 0.0
    case( 'y' );   fx1 = 0.0; fy1 = 1.0; fz1 = 0.0
    case( 'z' );   fx1 = 0.0; fy1 = 0.0; fz1 = 1.0  !! positive upward to fit nature of observation
    case default
      write(STDERR,*) "no matching green_cmp"
      stop
    end select

    dt_dxyz = dt / ( dx*dy*dz )

    !! initialize with unrealistic value
    evlo1 = -12345.0
    evla1 = -12345.0
    call output__station_query( green_stnm, is_src, isrc, jsrc, ksrc, xsrc, ysrc, zsrc, evlo1, evla1 )

    !! Confirm if the computational domain contains the source
    call mpi_allreduce( is_src, src_somewhere, 1, MPI_LOGICAL, MPI_LOR, mpi_comm_world, ierr )
    call assert( src_somewhere )

    !! Distribute evlo & evla. Stored to evlo0 and evla0
    call mpi_allreduce( evlo1, evlo0, 1, MPI_REAL, MPI_MAX, mpi_comm_world, ierr )
    call mpi_allreduce( evla1, evla0, 1, MPI_REAL, MPI_MAX, mpi_comm_world, ierr )

    !!
    !! Read Green's function location table
    !!
    call std__getio(io)
    open(io,file=trim(fn_glst), action='read', status='old', iostat=ierr )
    call assert( ierr == 0 )
    call std__countline( io, ng0, '#' )

    !! First count-up Green's function grid points inside the MPI node
    ng = 0
    ntw = floor( real(nt-1) / real(ntdec_w) + 1.0 )

    do

      read(io,'(A256)', iostat=ierr) abuf
      if( ierr /= 0 ) exit
      abuf = trim(adjustl(abuf))
      if( abuf(1:1) == '#'          ) cycle !! comment line
      if( trim(adjustl(abuf)) == '' ) cycle !! blank line

      if( green_fmt == 'xyz' ) then
        read(abuf,*,iostat=ierr) xg0, yg0, zg0, gid0
        call assert( ierr == 0 )
      else if (green_fmt == 'llz' ) then !! geographic coordinate
        read(abuf,*,iostat=ierr) long0, latg0, zg0, gid0
        call assert( ierr == 0 )
        call geomap__g2c( long0, latg0, clon, clat, phi, xg0, yg0 )
      else
        call assert( .false. )
      end if
      call assert( 0 <= gid0 .and. gid0 <= 99999999 )

      !! horizontal distance
      dd = sqrt( (xg0-xsrc)**2 + (yg0-ysrc)**2 )
      if( dd > green_maxdist ) cycle

      ii = x2i( xg0, xbeg, real(dx) )
      jj = y2j( yg0, ybeg, real(dy) )
      kk = z2k( zg0, zbeg, real(dz) )

      if( ibeg <= ii .and. ii <= iend .and. &
          jbeg <= jj .and. jj <= jend ) then
        if( kob(ii,jj) <= kk .and. kk <= kend ) then !! Onsider only the solid part
          ng = ng + 1
        end if
      end if
    end do

    !! Read the file again with storing the grid location
    rewind(io)
    allocate( ig(ng), jg(ng), kg(ng), xg(ng), yg(ng), zg(ng), gid(ng), long(ng), latg(ng) )
    ng = 0
    do

      read(io,'(A256)', iostat=ierr) abuf
      if( ierr /= 0 ) exit
      abuf = trim(adjustl(abuf))
      if( abuf(1:1) == '#'          ) cycle !! comment line
      if( trim(adjustl(abuf)) == '' ) cycle !! blank line

      if( green_fmt == 'xyz' ) then
        read(abuf,*,iostat=ierr) xg0, yg0, zg0, gid0
        call assert( ierr == 0 )
        call geomap__c2g( xg0, yg0, clon, clat, phi, long0, latg0 )
      else if (green_fmt == 'llz' ) then !! geographic coordinate
        read(abuf,*,iostat=ierr) long0, latg0, zg0, gid0
        call assert( ierr == 0 )
        call geomap__g2c( long0, latg0, clon, clat, phi, xg0, yg0 )
      else
        call assert( .false. )
      end if
      call assert( 0 <= gid0 .and. gid0 <= 99999999 )

      !! horizontal distance
      dd = sqrt( (xg0-xsrc)**2 + (yg0-ysrc)**2 )
      if( dd > green_maxdist ) cycle

      ii = x2i( xg0, xbeg, real(dx) )
      jj = y2j( yg0, ybeg, real(dy) )
      kk = z2k( zg0, zbeg, real(dz) )

      if( ibeg <= ii .and. ii <= iend .and. &
          jbeg <= jj .and. jj <= jend ) then
        if( kob(ii,jj) <= kk .and. kk <= kend ) then !! Consider only the solid part
          ng = ng + 1
          ig(ng) = ii
          jg(ng) = jj
          kg(ng) = kk
          xg(ng) = xg0
          yg(ng) = yg0
          zg(ng) = zg0
          gid(ng) = gid0
          long(ng) = long0
          latg(ng) = latg0
        end if
      end if
    end do
    close(io)


    if( green_bforce ) then
      ncmp = 9
    else
      ncmp = 6
    end if
    allocate( sh(ncmp*ng) )

    do i=1,ng
      do j=1,ncmp
        call sac__init(sh((i-1)*ncmp+j))
      end do
    end do

    !! station location = pseudo event location
    sh(:)%stlo = evlo0
    sh(:)%stla = evla0
    sh(:)%stdp = zsrc * 1000

    if( green_bforce ) then
      cmp(1:9) = (/ 'mxx', 'myy', 'mzz', 'myz', 'mxz', 'mxy', 'fx_', 'fy_', 'fz_' /)
    else
      cmp(1:6) = (/ 'mxx', 'myy', 'mzz', 'myz', 'mxz', 'mxy' /)
    end if

    allocate(gf(ntw,ncmp*ng))
    allocate(fn(ncmp*ng))
    gf(:,:) = 0.0

    call system__call( 'mkdir -p '//trim(odir) // '/green/' // trim(green_stnm) )

    do i=1, ng

      write(cid8,'(I8.8)') gid(i)

      do j=1, ncmp
        k = (i-1)*ncmp + j

        fn(k) = trim(odir)  // '/green/' // trim(green_stnm) // '/' // &
            trim(title) // '__' // cid8 // '__' // trim(green_stnm) // '__' //  green_cmp // '__' // &
            trim(cmp(j)) // '__.sac'
      end do

      select case (green_cmp)
      case( 'x' )
        do j=1, ncmp
          k = (i-1)*ncmp + j
          sh(k)%kcmpnm = 'G_Vx_' // cmp(j)
          sh(k)%cmpinc = 90.0
          sh(k)%cmpaz  =  0.0 + phi
        end do
      case( 'y' )
        do j=1, ncmp
          k = (i-1)*ncmp + j
          sh(k)%kcmpnm = 'G_Vy_' // cmp(j)
          sh(k)%cmpinc = 90.0
          sh(k)%cmpaz  = 90.0 + phi
        end do
      case( 'z' )
        do j=1, ncmp
          k = (i-1)*ncmp + j
          sh(k)%kcmpnm = 'G_Vz_' // cmp(j)
          sh(k)%cmpinc = 0.0
          sh(k)%cmpaz  = 0.0
        end do
      end select

      !! headers (grid location = event location)
      do j=1, ncmp
        k = (i-1)*ncmp + j

        sh(k)%kevnm = cid8
        sh(k)%kstnm = trim(green_stnm)
        sh(k)%evlo = long(i)
        sh(k)%evla = latg(i)
        sh(k)%evdp = zg(i)*1000 ! in m-unit

        sh(k)%tim   = exedate
        if( j<=6 ) then
          sh(k)%idep = 7
        else
          sh(k)%idep = 6
        end if


        call daytim__localtime( sh(k)%tim, &
            sh(k)%nzyear, sh(k)%nzmonth, sh(k)%nzday, sh(k)%nzhour, sh(k)%nzmin, sh(k)%nzsec )
        call daytim__ymd2jul  ( sh(k)%nzyear, sh(k)%nzmonth, sh(k)%nzday, sh(k)%nzjday )

        sh(k)%nzmsec = 0
        sh(k)%b     = tbeg
        sh(k)%delta = ntdec_w * dt
        sh(k)%npts  = ntw

        !! grid location specified by input
        sh(k)%user0 = xg(i)
        sh(k)%user1 = yg(i)
        sh(k)%user2 = zg(i)

        !! true location in grid system
        sh(k)%user3 = dble(i2x( ig(i), xbeg, real(dx) ))
        sh(k)%user4 = dble(j2y( jg(i), ybeg, real(dy) ))
        sh(k)%user5 = dble(k2z( kg(i), zbeg, real(dz) ))

        sh(k)%user6 = clon !< coordinate
        sh(k)%user7 = clat !< coordinate
        sh(k)%user8 = phi  !< coordinate

      end do
    end do

    call pwatch__off( 'green__setup' )


  end subroutine green__setup
  !! --------------------------------------------------------------------------------------------------------------------------- !!

  !! --------------------------------------------------------------------------------------------------------------------------- !!
  subroutine green__store(it)

    integer, intent(in) :: it
    integer :: i
    integer :: ii, jj, kk
    integer :: itw
    real(MP) :: dxVx, dxVy, dxVz, dyVx, dyVy, dyVz, dzVx, dzVy, dzVz
    real(MP) :: dxVy1, dxVy2, dxVy3, dxVy4
    real(MP) :: dxVz1, dxVz2, dxVz3, dxVz4
    real(MP) :: dyVx1, dyVx2, dyVx3, dyVx4
    real(MP) :: dyVz1, dyVz2, dyVz3, dyVz4
    real(MP) :: dzVx1, dzVx2, dzVx3, dzVx4
    real(MP) :: dzVy1, dzVy2, dzVy3, dzVy4
    integer, save :: pol

    if( .not. green_mode ) return

    call pwatch__on( 'green__store' )


    if( it == 1 ) then
      allocate( dxux(ng), dxuy(ng), dxuz(ng), dyux(ng), dyuy(ng),dyuz(ng),dzux(ng),dzuy(ng),dzuz(ng) )

      dxUx(:) = 0.0
      dxUy(:) = 0.0
      dxUz(:) = 0.0
      dyUx(:) = 0.0
      dyUy(:) = 0.0
      dyUz(:) = 0.0
      dzUx(:) = 0.0
      dzUy(:) = 0.0
      dzUz(:) = 0.0


      if( green_bforce ) then
        allocate( ux(ng), uy(ng), uz(ng) )
        ux(:)   = 0.0
        uy(:)   = 0.0
        uz(:)   = 0.0
      end if

      if( green_cmp == 'z' ) then
        pol = -1
      else
        pol = 1
      end if


    end if

    do i=1, ng
      ii = ig(i)
      jj = jg(i)
      kk = kg(i)

      dxVx  = (  Vx(kk  ,ii  ,jj  ) - Vx(kk  ,ii-1,jj  )  ) * r40x  -  (  Vx(kk  ,ii+1,jj  ) - Vx(kk  ,ii-2,jj  )  ) * r41x
      dyVy  = (  Vy(kk  ,ii  ,jj  ) - Vy(kk  ,ii  ,jj-1)  ) * r40y  -  (  Vy(kk  ,ii  ,jj+1) - Vy(kk  ,ii  ,jj-2)  ) * r41y
      dzVz  = (  Vz(kk  ,ii  ,jj  ) - Vz(kk-1,ii  ,jj  )  ) * r40z  -  (  Vz(kk+1,ii  ,jj  ) - Vz(kk-2,ii  ,jj  )  ) * r41z

      dxVy1 = (  Vy(kk  ,ii+1,jj  ) - Vy(kk  ,ii  ,jj  )  ) * r40x  -  (  Vy(kk  ,ii+2,jj  ) - Vy(kk  ,ii-1,jj  )  ) * r41x
      dxVy2 = (  Vy(kk  ,ii+1,jj-1) - Vy(kk  ,ii  ,jj-1)  ) * r40x  -  (  Vy(kk  ,ii+2,jj-1) - Vy(kk  ,ii-1,jj-1)  ) * r41x
      dxVy3 = (  Vy(kk  ,ii  ,jj  ) - Vy(kk  ,ii-1,jj  )  ) * r40x  -  (  Vy(kk  ,ii+1,jj  ) - Vy(kk  ,ii-2,jj  )  ) * r41x
      dxVy4 = (  Vy(kk  ,ii  ,jj-1) - Vy(kk  ,ii-1,jj-1)  ) * r40x  -  (  Vy(kk  ,ii+1,jj-1) - Vy(kk  ,ii-2,jj-1)  ) * r41x

      dxVz1 = (  Vz(kk  ,ii+1,jj  ) - Vz(kk  ,ii  ,jj  )  ) * r40x  -  (  Vz(kk  ,ii+2,jj  ) - Vz(kk  ,ii-1,jj  )  ) * r41x
      dxVz2 = (  Vz(kk-1,ii+1,jj  ) - Vz(kk-1,ii  ,jj  )  ) * r40x  -  (  Vz(kk-1,ii+2,jj  ) - Vz(kk-1,ii-1,jj  )  ) * r41x
      dxVz3 = (  Vz(kk  ,ii  ,jj  ) - Vz(kk  ,ii-1,jj  )  ) * r40x  -  (  Vz(kk  ,ii+1,jj  ) - Vz(kk  ,ii-2,jj  )  ) * r41x
      dxVz4 = (  Vz(kk-1,ii  ,jj  ) - Vz(kk-1,ii-1,jj  )  ) * r40x  -  (  Vz(kk-1,ii+1,jj  ) - Vz(kk-1,ii-2,jj  )  ) * r41x

      dyVx1 = (  Vx(kk  ,ii  ,jj+1) - Vx(kk  ,ii  ,jj  )  ) * r40y  -  (  Vx(kk  ,ii  ,jj+2) - Vx(kk  ,ii  ,jj-1)  ) * r41y
      dyVx2 = (  Vx(kk  ,ii-1,jj+1) - Vx(kk  ,ii-1,jj  )  ) * r40y  -  (  Vx(kk  ,ii-1,jj+2) - Vx(kk  ,ii-1,jj-1)  ) * r41y
      dyVx3 = (  Vx(kk  ,ii  ,jj  ) - Vx(kk  ,ii  ,jj-1)  ) * r40y  -  (  Vx(kk  ,ii  ,jj+1) - Vx(kk  ,ii  ,jj-2)  ) * r41y
      dyVx4 = (  Vx(kk  ,ii-1,jj  ) - Vx(kk  ,ii-1,jj-1)  ) * r40y  -  (  Vx(kk  ,ii-1,jj+1) - Vx(kk  ,ii-1,jj-2)  ) * r41y

      dyVz1 = (  Vz(kk  ,ii  ,jj+1) - Vz(kk  ,ii  ,jj  )  ) * r40y  -  (  Vz(kk  ,ii  ,jj+2) - Vz(kk  ,ii  ,jj-1)  ) * r41y
      dyVz2 = (  Vz(kk-1,ii  ,jj+1) - Vz(kk-1,ii  ,jj  )  ) * r40y  -  (  Vz(kk-1,ii  ,jj+2) - Vz(kk-1,ii  ,jj-1)  ) * r41y
      dyVz3 = (  Vz(kk  ,ii  ,jj  ) - Vz(kk  ,ii  ,jj-1)  ) * r40y  -  (  Vz(kk  ,ii  ,jj+1) - Vz(kk  ,ii  ,jj-2)  ) * r41y
      dyVz4 = (  Vz(kk-1,ii  ,jj  ) - Vz(kk-1,ii  ,jj-1)  ) * r40y  -  (  Vz(kk-1,ii  ,jj+1) - Vz(kk-1,ii  ,jj-2)  ) * r41y

      dzVx1 = (  Vx(kk+1,ii  ,jj  ) - Vx(kk  ,ii  ,jj  )  ) * r40z  -  (  Vx(kk+2,ii  ,jj  ) - Vx(kk-1,ii  ,jj  )  ) * r41z
      dzVx2 = (  Vx(kk+1,ii-1,jj  ) - Vx(kk  ,ii-1,jj  )  ) * r40z  -  (  Vx(kk+2,ii-1,jj  ) - Vx(kk-1,ii-1,jj  )  ) * r41z
      dzVx3 = (  Vx(kk  ,ii  ,jj  ) - Vx(kk-1,ii  ,jj  )  ) * r40z  -  (  Vx(kk+1,ii  ,jj  ) - Vx(kk-2,ii  ,jj  )  ) * r41z
      dzVx4 = (  Vx(kk  ,ii-1,jj  ) - Vx(kk-1,ii-1,jj  )  ) * r40z  -  (  Vx(kk+1,ii-1,jj  ) - Vx(kk-2,ii-1,jj  )  ) * r41z

      dzVy1 = (  Vy(kk+1,ii  ,jj  ) - Vy(kk  ,ii  ,jj  )  ) * r40z  -  (  Vy(kk+2,ii  ,jj  ) - Vy(kk-1,ii  ,jj  )  ) * r41z
      dzVy2 = (  Vy(kk+1,ii  ,jj-1) - Vy(kk  ,ii  ,jj-1)  ) * r40z  -  (  Vy(kk+2,ii  ,jj-1) - Vy(kk-1,ii  ,jj-1)  ) * r41z
      dzVy3 = (  Vy(kk  ,ii  ,jj  ) - Vy(kk-1,ii  ,jj  )  ) * r40z  -  (  Vy(kk+1,ii  ,jj  ) - Vy(kk-2,ii  ,jj  )  ) * r41z
      dzVy4 = (  Vy(kk  ,ii  ,jj-1) - Vy(kk-1,ii  ,jj-1)  ) * r40z  -  (  Vy(kk+1,ii  ,jj-1) - Vy(kk-2,ii  ,jj-1)  ) * r41z


      dxVy = ( dxVy1 + dxVy2 + dxVy3 + dxVy4 ) * 0.25
      dxVz = ( dxVz1 + dxVz2 + dxVz3 + dxVz4 ) * 0.25
      dyVx = ( dyVx1 + dyVx2 + dyVx3 + dyVx4 ) * 0.25
      dyVz = ( dyVz1 + dyVz2 + dyVz3 + dyVz4 ) * 0.25
      dzVx = ( dzVx1 + dzVx2 + dzVx3 + dzVx4 ) * 0.25
      dzVy = ( dzVy1 + dzVy2 + dzVy3 + dzVy4 ) * 0.25

      dxUx(i) = dxUx(i) + dxVx * dt
      dxUy(i) = dxUy(i) + dxVy * dt
      dxUz(i) = dxUz(i) + dxVz * dt
      dyUx(i) = dyUx(i) + dyVx * dt
      dyUy(i) = dyUy(i) + dyVy * dt
      dyUz(i) = dyUz(i) + dyVz * dt
      dzUx(i) = dzUx(i) + dzVx * dt
      dzUy(i) = dzUy(i) + dzVy * dt
      dzUz(i) = dzUz(i) + dzVz * dt

      if( green_bforce ) then
        ux(i) = ux(i) + 0.5 * ( Vx(kk,ii,jj)+Vx(kk,ii-1,jj) ) * dt
        uy(i) = uy(i) + 0.5 * ( Vy(kk,ii,jj)+Vy(kk,ii,jj-1) ) * dt
        uz(i) = uz(i) + 0.5 * ( Vz(kk,ii,jj)+Vz(kk-1,ii,jj) ) * dt
      end if


    end do

    if( mod( it-1, ntdec_w ) == 0 ) then

      itw = (it-1) / ntdec_w + 1

      do i=1, ng

        ii = ig(i)
        jj = jg(i)
        kk = kg(i)

        gf(itw,(i-1)*ncmp+1) =  dxUx(i)             * UC_DERIV * 1e9 ! m/s -> nm/s
        gf(itw,(i-1)*ncmp+2) =  dyUy(i)             * UC_DERIV * 1e9 ! m/s -> nm/s
        gf(itw,(i-1)*ncmp+3) =  dzUz(i)             * UC_DERIV * 1e9 ! m/s -> nm/s
        gf(itw,(i-1)*ncmp+4) = (dyUz(i) + dzUy(i) ) * UC_DERIV * 1e9 ! m/s -> nm/s
        gf(itw,(i-1)*ncmp+5) = (dxUz(i) + dzUx(i) ) * UC_DERIV * 1e9 ! m/s -> nm/s
        gf(itw,(i-1)*ncmp+6) = (dxUy(i) + dyUx(i) ) * UC_DERIV * 1e9 ! m/s -> nm/s

        if( green_bforce ) then
          gf(itw,(i-1)*ncmp+7) = ux(i) * UC_BF * 1e9 ! m -> nm
          gf(itw,(i-1)*ncmp+8) = uy(i) * UC_BF * 1e9 ! m -> nm
          gf(itw,(i-1)*ncmp+9) = uz(i) * UC_BF * 1e9 ! m -> nm
        end if

      end do

    end if

    call pwatch__off( 'green__store' )

  end subroutine green__store
  !! --------------------------------------------------------------------------------------------------------------------------- !!

  !! --------------------------------------------------------------------------------------------------------------------------- !!
  subroutine green__export()

    integer :: i, j, k
    integer :: io
    character(6) :: cmyid
    character(256) :: fn_csf, fn_wav
    !! --

    if( .not. green_mode ) return

    call pwatch__on( 'green__export' )

    !! Make positive upward for z-component
    if( green_cmp == 'z' ) then
      gf(:,:) = -gf(:,:)
    end if

    if( wav_format == 'sac' ) then
      do i=1, ng
        do j=1, ncmp
          k = (i-1)*ncmp + j
          call sac__write( fn(k), sh(k), gf(:,k), .true. )
        end do
      end do
    else if ( wav_format == 'csf' ) then

      if( ng>0 ) then

        write(cmyid,'(I6.6)') myid
        fn_csf = trim(odir)  // '/green/' // trim(sh(1)%kstnm) // '/' // &
            trim(title) //  '__' // trim(sh(1)%kstnm) // '__' // &
            trim(green_cmp) // '__' // trim(cmyid) // '__.csf'
        call csf__write( fn_csf, ng*ncmp, sh(1)%npts, sh, gf, .true. )
      end if

    else if ( wav_format == 'wav' ) then

      if( ng>0 ) then

        write(cmyid,'(I6.6)') myid
        fn_wav = trim(odir)  // '/green/' // trim(sh(1)%kstnm) // '/' // &
            trim(title) //  '__' // trim(sh(1)%kstnm) // '__' // &
            trim(green_cmp) // '__' // trim(cmyid) // '__.wav'

#ifdef _ES
        call std__getio(io, is_big=.true.)
        open(io, file=trim(fn_wav), form='unformatted', action='write', status='replace')
#else
        call std__getio(io, is_big=.true.) 
        open(io, file=trim(fn_wav), access='stream', form='unformatted', action='write', status='replace')
#endif

        write(io) ng*ncmp, ntw, title
        write(io) sh(:)
        write(io) gf(:,:)

        close(io)

      end if


    end if

    call pwatch__off( 'green__export' )

  end subroutine green__export
  !! --------------------------------------------------------------------------------------------------------------------------- !!

  !! --------------------------------------------------------------------------------------------------------------------------- !!
  subroutine green__source( it )

    integer, intent(in) :: it
    real(SP) :: stf
    real(SP) :: t
    real(SP) :: fx, fy, fz

    if( .not. green_mode ) return
    call pwatch__on( 'green__source' )

    if( .not. is_src ) then
      call pwatch__off( 'green__source' )
      return
    end if

    !! velocity needs half-grid shift in time
    t = n2t( it, tbeg, dt ) + dt/2.0

    select case ( trim( stftype ) )
    case ( 'boxcar'   );  stf = boxcar   ( t, green_tbeg, green_trise )
    case ( 'triangle' );  stf = triangle ( t, green_tbeg, green_trise )
    case ( 'herrmann' );  stf = herrmann ( t, green_tbeg, green_trise )
    case ( 'kupper'   );  stf = kupper   ( t, green_tbeg, green_trise )
    case ( 'cosine'   );  stf = cosine   ( t, green_tbeg, green_trise )
    case default;         stf = kupper   ( t, green_tbeg, green_trise )
    end select


    fx = fx1 * dt_dxyz * stf
    fy = fy1 * dt_dxyz * stf
    fz = fz1 * dt_dxyz * stf

    !! Reciprocal density must be avaraged for keeping accuracy
    Vx(ksrc  ,isrc  ,jsrc  ) = Vx(ksrc  ,isrc  ,jsrc  ) + bx(ksrc  ,isrc  ,jsrc  ) * fx / 2
    Vx(ksrc  ,isrc-1,jsrc  ) = Vx(ksrc  ,isrc-1,jsrc  ) + bx(ksrc  ,isrc-1,jsrc  ) * fx / 2
    Vy(ksrc  ,isrc  ,jsrc  ) = Vy(ksrc  ,isrc  ,jsrc  ) + by(ksrc  ,isrc  ,jsrc  ) * fy / 2
    Vy(ksrc  ,isrc  ,jsrc-1) = Vy(ksrc  ,isrc  ,jsrc-1) + by(ksrc  ,isrc  ,jsrc-1) * fy / 2
    Vz(ksrc  ,isrc  ,jsrc  ) = Vz(ksrc  ,isrc  ,jsrc  ) + bz(ksrc  ,isrc  ,jsrc  ) * fz / 2
    Vz(ksrc-1,isrc  ,jsrc  ) = Vz(ksrc-1,isrc  ,jsrc  ) + bz(ksrc-1,isrc  ,jsrc  ) * fz / 2

    call pwatch__off( 'green__source' )

  end subroutine green__source
  !! --------------------------------------------------------------------------------------------------------------------------- !!


end module m_green
!! ----------------------------------------------------------------------------------------------------------------------------- !!
