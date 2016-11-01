!! ----------------------------------------------------------------------------------------------------------------------------- !!
!>
!! Snapshot/waveform output
!!
!! @copyright
!!   Copyright 2013-2016 Takuto Maeda. All rights reserved. This project is released under the MIT license.
!<
!! ----
#include "m_debug.h"
module m_output

  !! -- Dependency
  use m_std
  use m_debug
  use m_global
  use m_pwatch
  use m_fdtool
  use m_daytim
  use m_wsac
  use m_readini
  use m_geomap
  use m_system
#ifdef _NETCDF
  use netcdf
#endif

  !! -- Declarations
  implicit none
  private
  save

  public :: output__setup
  public :: output__write_snap
  public :: output__store_wav
  public :: output__export_wav
  public :: output__checkpoint
  public :: output__restart
  public :: output__closefiles
  public :: output__station_query

  !! -- Internal Parameters
  character(8), parameter :: BINARY_TYPE= "STREAMIO"
  character(8), parameter :: CODE_TYPE  = "SWPC_3D "   !!< FIXED parameter for file header
  integer,      parameter :: HEADER_VERSION = 6
  !! header_version history
  !! 110311 : original
  !! 2      : added binary_type (character(8))
  !! 3      : added Lz and zeta for curvilinear
  !! 4      : added w for curvilinear
  !! 5      : added coordinate rotation angle phi
  !! 6      : seism -> swpc

  !! -- Internal variables
  type snp
    logical :: sw     ! true for output
    integer :: io     ! file I/O number (or netcdf file id)
    integer :: ionode ! output MPI node
    integer :: nsnp
    character(2) :: snaptype ! snapshot type
    character(2) :: coordinate
    integer :: nmed
    integer :: na1, na2 ! absorbing layer
    real    :: ds1, ds2 ! spatial grid width

    !! variables for netcdf mode
    integer :: did_x1, did_x2, did_t ! dimension id for independent vars
    integer :: vid_x1, vid_x2, vid_t ! variable  id for independent vars
    integer :: varid(10)           ! variable  id for dependent vars
    integer :: medid(10)            ! medium array id
    real    :: vmax(10), vmin(10)  ! max/min of dependent vars
    character(10) :: vname(10)
    character(10) :: vunit(10)
  end type snp

  type(snp) :: xy_ps, yz_ps, xz_ps, fs_ps, ob_ps, xy_v, yz_v, xz_v, fs_v, ob_v, xy_u, yz_u, xz_u, fs_u, ob_u
  logical   :: sw_wav, sw_wav_v, sw_wav_u

  !! switch
  integer   :: ntdec_s, ntdec_w                                !< time step decimation factor: Snap and Waves
  integer   :: idec, jdec, kdec                                !< spatial decimation factor: x, y, z directions
  character(2) :: st_format                                       !< station file format
  character(3) :: wav_format                                   !< sac or csf
  
  real(SP) :: z0_xy, x0_yz, y0_xz ! < danmen
  integer  :: k0_xy, i0_yz, j0_xz ! < danmen
  character(256) :: fn_stloc

  integer :: nxs, nys, nzs !< snapshot grid size
  real(SP), allocatable :: xsnp(:), ysnp(:), zsnp(:)

  !! station
  integer :: nst
  real(SP), allocatable :: xst(:), yst(:), zst(:)
  integer,  allocatable :: ist(:), jst(:), kst(:)
  real(SP), allocatable :: stlo(:), stla(:)
  character(8), allocatable :: stnm(:)

  !! waveform
  integer :: ntw ! number of wave samples
  real(SP), allocatable :: vxst(:,:), vyst(:,:), vzst(:,:), uxst(:,:), uyst(:,:), uzst(:,:)
  type(sac__hdr), allocatable :: sh(:,:)
  real(SP), allocatable :: ux(:), uy(:), uz(:)

  !! I/O area in the node
  integer :: is0, is1, js0, js1, ks0, ks1

  !! derivative coefficient
  real(SP) :: r20x, r20y, r20z

  !! snapshot data format
  character(6) :: snp_format ! native or netcdf

  !! displacement snapshot buffer
  real(SP), allocatable :: buf_yz_u(:,:,:), buf_xz_u(:,:,:), buf_xy_u(:,:,:), buf_fs_u(:,:,:), buf_ob_u(:,:,:)

contains


  !! --------------------------------------------------------------------------------------------------------------------------- !!
  !>
  !! Setup
  !!
  !! @see
  !! #2012-41 #32013-00440
  !<
  !! ----
  subroutine output__setup( io_prm )
    integer, intent(in) :: io_prm
    logical :: green_mode
    integer :: i, j, k, ii, jj, kk
    integer :: ierr
    !! ----

    call pwatch__on( "output__setup" )

    !!
    !! read parameter
    !!

    call readini( io_prm, 'xy_ps%sw',   xy_ps%sw,   .false. )
    call readini( io_prm, 'yz_ps%sw',   yz_ps%sw,   .false. )
    call readini( io_prm, 'xz_ps%sw',   xz_ps%sw,   .false. )
    call readini( io_prm, 'ob_ps%sw',   ob_ps%sw,   .false. )
    call readini( io_prm, 'fs_ps%sw',   fs_ps%sw,   .false. )
    call readini( io_prm, 'xy_v%sw',    xy_v%sw,    .false. )
    call readini( io_prm, 'yz_v%sw',    yz_v%sw,    .false. )
    call readini( io_prm, 'xz_v%sw',    xz_v%sw,    .false. )
    call readini( io_prm, 'ob_v%sw',    ob_v%sw,    .false. )
    call readini( io_prm, 'fs_v%sw',    fs_v%sw,    .false. )
    call readini( io_prm, 'xy_u%sw',    xy_u%sw,    .false. )
    call readini( io_prm, 'yz_u%sw',    yz_u%sw,    .false. )
    call readini( io_prm, 'xz_u%sw',    xz_u%sw,    .false. )
    call readini( io_prm, 'ob_u%sw',    ob_u%sw,    .false. )
    call readini( io_prm, 'fs_u%sw',    fs_u%sw,    .false. )
    call readini( io_prm, 'z0_xy',      z0_xy,       max( min(10.0,zend), zbeg) )
    call readini( io_prm, 'x0_yz',      x0_yz,       max( min( 0.0,xend), xbeg) )
    call readini( io_prm, 'y0_xz',      y0_xz,       max( min( 0.0,yend), ybeg) )
    call readini( io_prm, 'idec',       idec,        1       )
    call readini( io_prm, 'jdec',       jdec,        1       )
    call readini( io_prm, 'kdec',       kdec,        1       )
    call readini( io_prm, 'ntdec_s',    ntdec_s,     10      )
    call readini( io_prm, 'ntdec_w',    ntdec_w,     10      )
    call readini( io_prm, 'st_format',  st_format,  'xy'    )
    call readini( io_prm, 'fn_stloc',   fn_stloc,   ''      )
    call readini( io_prm, 'sw_wav_v',   sw_wav_v,   .false. )
    call readini( io_prm, 'sw_wav_u',   sw_wav_u,   .false. )
    call readini( io_prm, 'snp_format', snp_format, 'native' )
    call readini( io_prm, 'wav_format', wav_format, 'sac' ) 

    !! Do not output waveform for Green's function mode
    call readini( io_prm, 'green_mode', green_mode, .false. )

    sw_wav = ( sw_wav_v .or. sw_wav_u )
    if( green_mode ) then
      sw_wav_v = .false.
      sw_wav_u = .false.
      sw_wav   = .true.
    end if


    !!!!
    !!!! snapshot
    !!!!

    !!
    !! snapshot size #2013-0440
    !!
    nxs = ( nx + (idec/2) ) / idec
    nys = ( ny + (jdec/2) ) / jdec
    nzs = ( nz + (kdec/2) ) / kdec

    !!
    !! coordinate
    !!
    allocate( xsnp(nxs), ysnp(nys), zsnp(nzs) )
    do i=1, nxs
      ii =  i*idec - (idec/2)
      xsnp(i) = i2x( ii, xbeg, real(dx) )
    end do
    do j=1, nys
      jj = j*jdec - (jdec/2)
      ysnp(j) = j2y( jj, ybeg, real(dy) )
    end do
    do k=1, nzs
      kk = k*kdec - (kdec/2)
      zsnp(k) = k2z( kk, zbeg, real(dz) )
    end do

    !! snapshot region covered by the MPI node
    is0 = ceiling( ( ibeg + idec/2) / real(idec) )
    is1 = floor  ( ( iend + idec/2) / real(idec) )
    js0 = ceiling( ( jbeg + jdec/2) / real(jdec) )
    js1 = floor  ( ( jend + jdec/2) / real(jdec) )
    ks0 = ceiling( ( kbeg + kdec/2) / real(kdec) )
    ks1 = floor  ( ( kend + kdec/2) / real(kdec) )

    !! snapshot location grid
    k0_xy = z2k( z0_xy, zbeg, real(dz) )
    i0_yz = x2i( x0_yz, xbeg, real(dx) )
    j0_xz = y2j( y0_xz, ybeg, real(dy) )


    !!
    !! output node definition: cyclic
    !!
    yz_ps%ionode = mod( 0, nproc )
    xz_ps%ionode = mod( 1, nproc )
    xy_ps%ionode = mod( 2, nproc )
    fs_ps%ionode = mod( 3, nproc )
    ob_ps%ionode = mod( 4, nproc )
    yz_v%ionode  = mod( 5, nproc )
    xz_v%ionode  = mod( 6, nproc )
    xy_v%ionode  = mod( 7, nproc )
    fs_v%ionode  = mod( 8, nproc )
    ob_v%ionode  = mod( 9, nproc )
    yz_u%ionode  = mod( 10, nproc )
    xz_u%ionode  = mod( 11, nproc )
    xy_u%ionode  = mod( 12, nproc )
    fs_u%ionode  = mod( 13, nproc )
    ob_u%ionode  = mod( 14, nproc )

    !!
    !! snapshot settings
    !!
    yz_ps % nsnp = 4;  yz_v  % nsnp = 3;  yz_u  % nsnp = 3;
    xz_ps % nsnp = 4;  xz_v  % nsnp = 3;  xz_u  % nsnp = 3;
    xy_ps % nsnp = 4;  xy_v  % nsnp = 3;  xy_u  % nsnp = 3;
    ob_ps % nsnp = 4;  ob_v  % nsnp = 3;  ob_u  % nsnp = 3;
    fs_ps % nsnp = 4;  fs_v  % nsnp = 3;  fs_u  % nsnp = 3;

    ! horizontal medium continas topography, lon & lat
    yz_ps % nmed = 3;  yz_v  % nmed = 3;  yz_u  % nmed = 3;
    xz_ps % nmed = 3;  xz_v  % nmed = 3;  xz_u  % nmed = 3;
    xy_ps % nmed = 6;  xy_v  % nmed = 6;  xy_u  % nmed = 6;
    ob_ps % nmed = 6;  ob_v  % nmed = 6;  ob_u  % nmed = 6;
    fs_ps % nmed = 6;  fs_v  % nmed = 6;  fs_u  % nmed = 6;


    yz_ps % snaptype = 'ps';  yz_v % snaptype = 'v3';  yz_u % snaptype = 'u3';
    xz_ps % snaptype = 'ps';  xz_v % snaptype = 'v3';  xz_u % snaptype = 'u3';
    xy_ps % snaptype = 'ps';  xy_v % snaptype = 'v3';  xy_u % snaptype = 'u3';
    ob_ps % snaptype = 'ps';  ob_v % snaptype = 'v3';  ob_u % snaptype = 'u3';
    fs_ps % snaptype = 'ps';  fs_v % snaptype = 'v3';  fs_u % snaptype = 'u3';

    yz_ps % vname(1) = 'div'; yz_ps % vname(2) = 'rot_x';  yz_ps % vname(3) = 'rot_y';  yz_ps % vname(4) = 'rot_z'
    xz_ps % vname(1) = 'div'; xz_ps % vname(2) = 'rot_x';  xz_ps % vname(3) = 'rot_y';  xz_ps % vname(4) = 'rot_z'
    xy_ps % vname(1) = 'div'; xy_ps % vname(2) = 'rot_x';  xy_ps % vname(3) = 'rot_y';  xy_ps % vname(4) = 'rot_z'
    fs_ps % vname(1) = 'div'; fs_ps % vname(2) = 'rot_x';  fs_ps % vname(3) = 'rot_y';  fs_ps % vname(4) = 'rot_z'
    ob_ps % vname(1) = 'div'; ob_ps % vname(2) = 'rot_x';  ob_ps % vname(3) = 'rot_y';  ob_ps % vname(4) = 'rot_z'

    yz_v % vname(1) = 'Vx'; yz_v % vname(2) = 'Vy';  yz_v % vname(3) = 'Vz'
    xz_v % vname(1) = 'Vx'; xz_v % vname(2) = 'Vy';  xz_v % vname(3) = 'Vz'
    xy_v % vname(1) = 'Vx'; xy_v % vname(2) = 'Vy';  xy_v % vname(3) = 'Vz'
    fs_v % vname(1) = 'Vx'; fs_v % vname(2) = 'Vy';  fs_v % vname(3) = 'Vz'
    ob_v % vname(1) = 'Vx'; ob_v % vname(2) = 'Vy';  ob_v % vname(3) = 'Vz'

    yz_u % vname(1) = 'Ux'; yz_u % vname(2) = 'Uy';  yz_u % vname(3) = 'Uz'
    xz_u % vname(1) = 'Ux'; xz_u % vname(2) = 'Uy';  xz_u % vname(3) = 'Uz'
    xy_u % vname(1) = 'Ux'; xy_u % vname(2) = 'Uy';  xy_u % vname(3) = 'Uz'
    fs_u % vname(1) = 'Ux'; fs_u % vname(2) = 'Uy';  fs_u % vname(3) = 'Uz'
    ob_u % vname(1) = 'Ux'; ob_u % vname(2) = 'Uy';  ob_u % vname(3) = 'Uz'

    yz_ps % vunit(1) = '1/s'; yz_ps % vunit(2) = '1/s';  yz_ps % vunit(3) = '1/s';  yz_ps % vunit(4) = '1/s'
    xz_ps % vunit(1) = '1/s'; xz_ps % vunit(2) = '1/s';  xz_ps % vunit(3) = '1/s';  xz_ps % vunit(4) = '1/s'
    xy_ps % vunit(1) = '1/s'; xy_ps % vunit(2) = '1/s';  xy_ps % vunit(3) = '1/s';  xy_ps % vunit(4) = '1/s'
    fs_ps % vunit(1) = '1/s'; fs_ps % vunit(2) = '1/s';  fs_ps % vunit(3) = '1/s';  fs_ps % vunit(4) = '1/s'
    ob_ps % vunit(1) = '1/s'; ob_ps % vunit(2) = '1/s';  ob_ps % vunit(3) = '1/s';  ob_ps % vunit(4) = '1/s'

    yz_v % vunit(1) = 'm/s'; yz_v % vunit(2) = 'm/s';  yz_v % vunit(3) = 'm/s'
    xz_v % vunit(1) = 'm/s'; xz_v % vunit(2) = 'm/s';  xz_v % vunit(3) = 'm/s'
    xy_v % vunit(1) = 'm/s'; xy_v % vunit(2) = 'm/s';  xy_v % vunit(3) = 'm/s'
    fs_v % vunit(1) = 'm/s'; fs_v % vunit(2) = 'm/s';  fs_v % vunit(3) = 'm/s'
    ob_v % vunit(1) = 'm/s'; ob_v % vunit(2) = 'm/s';  ob_v % vunit(3) = 'm/s'

    yz_u % vunit(1) = 'm'; yz_u % vunit(2) = 'm';  yz_u % vunit(3) = 'm'
    xz_u % vunit(1) = 'm'; xz_u % vunit(2) = 'm';  xz_u % vunit(3) = 'm'
    xy_u % vunit(1) = 'm'; xy_u % vunit(2) = 'm';  xy_u % vunit(3) = 'm'
    fs_u % vunit(1) = 'm'; fs_u % vunit(2) = 'm';  fs_u % vunit(3) = 'm'
    ob_u % vunit(1) = 'm'; ob_u % vunit(2) = 'm';  ob_u % vunit(3) = 'm'

    yz_ps % coordinate = 'yz';  yz_v % coordinate = 'yz';  yz_u % coordinate = 'yz';
    xz_ps % coordinate = 'xz';  xz_v % coordinate = 'xz';  xz_u % coordinate = 'xz';
    xy_ps % coordinate = 'xy';  xy_v % coordinate = 'xy';  xy_u % coordinate = 'xy';
    ob_ps % coordinate = 'ob';  ob_v % coordinate = 'ob';  ob_u % coordinate = 'ob';
    fs_ps % coordinate = 'fs';  fs_v % coordinate = 'fs';  fs_u % coordinate = 'fs';

    !!
    !! open
    !!

    !!
    !! output settings
    !!
    if( snp_format == 'native' ) then

      if( yz_ps%sw ) call newfile_yz(   trim(odir) // '/' // trim(title) //'.yz.ps.snp',  yz_ps  )
      if( xz_ps%sw ) call newfile_xz(   trim(odir) // '/' // trim(title) //'.xz.ps.snp',  xz_ps  )
      if( xy_ps%sw ) call newfile_xy(   trim(odir) // '/' // trim(title) //'.xy.ps.snp',  xy_ps  )
      if( fs_ps%sw ) call newfile_xy(   trim(odir) // '/' // trim(title) //'.fs.ps.snp',  fs_ps  )
      if( ob_ps%sw ) call newfile_xy(   trim(odir) // '/' // trim(title) //'.ob.ps.snp',  ob_ps  )

      if( yz_v %sw ) call newfile_yz(   trim(odir) // '/' // trim(title) //'.yz.v.snp',   yz_v )
      if( xz_v %sw ) call newfile_xz(   trim(odir) // '/' // trim(title) //'.xz.v.snp',   xz_v )
      if( xy_v %sw ) call newfile_xy(   trim(odir) // '/' // trim(title) //'.xy.v.snp',   xy_v )
      if( fs_v %sw ) call newfile_xy(   trim(odir) // '/' // trim(title) //'.fs.v.snp',   fs_v )
      if( ob_v %sw ) call newfile_xy(   trim(odir) // '/' // trim(title) //'.ob.v.snp',   ob_v )

      if( yz_u %sw ) call newfile_yz(   trim(odir) // '/' // trim(title) //'.yz.u.snp',   yz_u )
      if( xz_u %sw ) call newfile_xz(   trim(odir) // '/' // trim(title) //'.xz.u.snp',   xz_u )
      if( xy_u %sw ) call newfile_xy(   trim(odir) // '/' // trim(title) //'.xy.u.snp',   xy_u )
      if( fs_u %sw ) call newfile_xy(   trim(odir) // '/' // trim(title) //'.fs.u.snp',   fs_u )
      if( ob_u %sw ) call newfile_xy(   trim(odir) // '/' // trim(title) //'.ob.u.snp',   ob_u )

    else

      if( yz_ps%sw ) call newfile_yz_nc(   trim(odir) // '/' // trim(title) //'.yz.ps.nc',  yz_ps  )
      if( xz_ps%sw ) call newfile_xz_nc(   trim(odir) // '/' // trim(title) //'.xz.ps.nc',  xz_ps  )
      if( xy_ps%sw ) call newfile_xy_nc(   trim(odir) // '/' // trim(title) //'.xy.ps.nc',  xy_ps  )
      if( fs_ps%sw ) call newfile_xy_nc(   trim(odir) // '/' // trim(title) //'.fs.ps.nc',  fs_ps  )
      if( ob_ps%sw ) call newfile_xy_nc(   trim(odir) // '/' // trim(title) //'.ob.ps.nc',  ob_ps  )

      if( yz_v %sw ) call newfile_yz_nc(   trim(odir) // '/' // trim(title) //'.yz.v.nc',   yz_v )
      if( xz_v %sw ) call newfile_xz_nc(   trim(odir) // '/' // trim(title) //'.xz.v.nc',   xz_v )
      if( xy_v %sw ) call newfile_xy_nc(   trim(odir) // '/' // trim(title) //'.xy.v.nc',   xy_v )
      if( fs_v %sw ) call newfile_xy_nc(   trim(odir) // '/' // trim(title) //'.fs.v.nc',   fs_v )
      if( ob_v %sw ) call newfile_xy_nc(   trim(odir) // '/' // trim(title) //'.ob.v.nc',   ob_v )

      if( yz_u %sw ) call newfile_yz_nc(   trim(odir) // '/' // trim(title) //'.yz.u.nc',   yz_u )
      if( xz_u %sw ) call newfile_xz_nc(   trim(odir) // '/' // trim(title) //'.xz.u.nc',   xz_u )
      if( xy_u %sw ) call newfile_xy_nc(   trim(odir) // '/' // trim(title) //'.xy.u.nc',   xy_u )
      if( fs_u %sw ) call newfile_xy_nc(   trim(odir) // '/' // trim(title) //'.fs.u.nc',   fs_u )
      if( ob_u %sw ) call newfile_xy_nc(   trim(odir) // '/' // trim(title) //'.ob.u.nc',   ob_u )

    end if

    !! for taking derivatives
    r20x = 1 / dx
    r20y = 1 / dy
    r20z = 1 / dz

    !! displacement snapshot buffers
    allocate( buf_yz_u(nys,nzs,3), buf_xz_u(nxs,nzs,3), buf_xy_u(nxs,nys,3), buf_fs_u(nxs,nys,3), buf_ob_u(nxs,nys,3) )
    buf_yz_u = 0.0
    buf_xz_u = 0.0
    buf_xy_u = 0.0
    buf_fs_u = 0.0
    buf_ob_u = 0.0

    !!!!
    !!!! waveform
    !!!!
    if( sw_wav ) then

      ntw = floor( real(nt-1)/real(ntdec_w) + 1.0 )

      call read_stinfo()
      call system__call('mkdir -p '//trim(odir)//'/wav > /dev/null 2>&1' )

    end if

    call mpi_barrier( mpi_comm_world, ierr )

    call pwatch__off( "output__setup" )

  end subroutine output__setup
  !! --------------------------------------------------------------------------------------------------------------------------- !!


  !! --------------------------------------------------------------------------------------------------------------------------- !!
  !>
  !<
  !! --
  subroutine output__export_wav()
    integer :: i
    character(256) :: fn1, fn2, fn3, fn4, fn5, fn6
    character(6) :: cid

    
    call pwatch__on("output__export_wav")

    if( nst == 0 ) then
      call pwatch__off("output__export_wav")
      return
    end if
    
    if( wav_format == 'sac' ) then
      do i=1, nst
        
        if( sw_wav_v ) then
          
          fn1 = trim(odir) // '/wav/' // trim(title) // '.' // trim(stnm(i)) // '.Vx.sac'
          fn2 = trim(odir) // '/wav/' // trim(title) // '.' // trim(stnm(i)) // '.Vy.sac'
          fn3 = trim(odir) // '/wav/' // trim(title) // '.' // trim(stnm(i)) // '.Vz.sac'
          
          call sac__write( fn1, sh(1,i), vxst(:,i), .true. )
          call sac__write( fn2, sh(2,i), vyst(:,i), .true. )
          call sac__write( fn3, sh(3,i), vzst(:,i), .true. )
          
        end if
        
        if( sw_wav_u ) then
          fn4 = trim(odir) // '/wav/' // trim(title) // '.' // trim(stnm(i)) // '.Ux.sac'
          fn5 = trim(odir) // '/wav/' // trim(title) // '.' // trim(stnm(i)) // '.Uy.sac'
          fn6 = trim(odir) // '/wav/' // trim(title) // '.' // trim(stnm(i)) // '.Uz.sac'
          call sac__write( fn4, sh(4,i), uxst(:,i), .true. )
          call sac__write( fn5, sh(5,i), uyst(:,i), .true. )
          call sac__write( fn6, sh(6,i), uzst(:,i), .true. )
        end if
        
      end do

    else
      write(cid,'(I6.6)') myid
      
      if( sw_wav_v ) then

        fn1 = trim(odir) // '/wav/' // trim(title) // '.' // trim(cid) // '.Vx.csf'
        fn2 = trim(odir) // '/wav/' // trim(title) // '.' // trim(cid) // '.Vy.csf'
        fn3 = trim(odir) // '/wav/' // trim(title) // '.' // trim(cid) // '.Vz.csf'
        call csf__write( fn1, nst, sh(1,1)%npts, sh(1,:), vxst(:,:), .true.)
        call csf__write( fn2, nst, sh(2,1)%npts, sh(2,:), vyst(:,:), .true.)
        call csf__write( fn3, nst, sh(3,1)%npts, sh(3,:), vzst(:,:), .true.)

      end if

      if( sw_wav_u ) then
        fn4 = trim(odir) // '/wav/' // trim(title) // '.' // trim(cid) // '.Ux.sac'
        fn5 = trim(odir) // '/wav/' // trim(title) // '.' // trim(cid) // '.Uy.sac'
        fn6 = trim(odir) // '/wav/' // trim(title) // '.' // trim(cid) // '.Uz.sac'
        call csf__write( fn4, nst, sh(4,1)%npts, sh(4,:), uxst(:,:), .true. )
        call csf__write( fn5, nst, sh(5,1)%npts, sh(5,:), uyst(:,:), .true. )
        call csf__write( fn6, nst, sh(6,1)%npts, sh(6,:), uzst(:,:), .true. )
      end if
        
    end if
    
    call pwatch__off("output__export_wav")


  end subroutine output__export_wav
  !! --------------------------------------------------------------------------------------------------------------------------- !!


  !! --------------------------------------------------------------------------------------------------------------------------- !!
  !>
  !! Set-up waveform file. Create, header output.
  !<
  !! --
  subroutine read_stinfo( )

    integer                   :: io_stlst
    real(SP), allocatable     :: xst_g(:), yst_g(:), zst_g(:), stlo_g(:), stla_g(:)
    integer, allocatable      :: ist_g(:), jst_g(:), kst_g(:), stid(:)
    character(3), allocatable :: zsw_g(:)
    character(8), allocatable :: stnm_g(:)
    character(256)            :: abuf
    integer                   :: ierr
    integer :: nst_g
    integer :: i, ii, j
    real(SP) :: mw
    !! ----

    !!!!
    !!!! Read station location file
    !!!!

    call std__getio(io_stlst)
    open(io_stlst, file=trim(fn_stloc), action='read', iostat = ierr, status='old' )

    if( ierr /= 0 ) then
      nst_g = 0          ! not exist
      sw_wav = .false.
      call info( "No station file found" )
      return
    else
      call std__countline( io_stlst, nst_g, "#" )
    end if

    !! allocate temp memory
    allocate( xst_g(nst_g), yst_g(nst_g), zst_g(nst_g), stnm_g(nst_g), zsw_g(nst_g), stlo_g(nst_g), stla_g(nst_g) )
    allocate( ist_g(nst_g), jst_g(nst_g), kst_g(nst_g), stid(nst_g) )


    !!
    !! First, store all station location into memory, and select station inside the MPI node
    !!
    nst = 0

    !! once recet total number of station
    !! count-up stations only inside the computation domain
    nst_g = 0
    i = 0

    !! moment magnitude for header part
    mw = moment_magnitude( m0 )

    do
      read(io_stlst,'(A256)', iostat=ierr) abuf
      if( ierr /= 0)  exit
      abuf = trim(adjustl(abuf))
      if( abuf(1:1) == "#" ) cycle ! neglect comment line
      if( trim(adjustl(abuf)) == "" ) cycle ! neglect blank line

      i = i + 1

      select case ( st_format )

      case( 'xy' )

        read(abuf,*,iostat=ierr) xst_g(i), yst_g(i), zst_g(i), stnm_g(i), zsw_g(i)
        call assert( ierr==0 )
        call geomap__c2g( xst_g(i), yst_g(i), clon, clat, phi, stlo_g(i), stla_g(i) )


      case( 'll' )

        read(abuf,*,iostat=ierr) stlo_g(i), stla_g(i), zst_g(i), stnm_g(i), zsw_g(i)
        call assert( ierr==0 )
        call assert( -360. <= stlo_g(i) .and. stlo_g(i) <= 360. )
        call assert(  -90. <= stla_g(i) .and. stla_g(i) <=  90. )
        call geomap__g2c( stlo_g(i), stla_g(i), clon, clat, phi, xst_g(i), yst_g(i) )

      case default

        call assert( .false. )

      end select


      ! digitize
      ist_g(i) = x2i ( xst_g(i), xbeg, real(dx) )
      jst_g(i) = y2j ( yst_g(i), ybeg, real(dy) )
      kst_g(i) = z2k ( zst_g(i), zbeg, real(dz) )


      !! check if the station is in the computational domain
      if(  i2x( 1, xbeg, real(dx) ) < xst_g(i) .and. xst_g(i) < i2x( nx, xbeg, real(dx) ) .and. &
      j2y( 1, ybeg, real(dy) ) < yst_g(i) .and. yst_g(i) < j2y( ny, ybeg, real(dy) ) .and. &
      kbeg       < kst_g(i) .and. kst_g(i) <          kend             )  then


      !! memorize station location
      nst_g = nst_g + 1

      !! MPI region check: memorize station number
      if( ibeg <= ist_g(i) .and. ist_g(i) <= iend .and. jbeg <= jst_g(i) .and. jst_g(i) <= jend ) then

        nst = nst + 1
        stid(nst) = i

      end if

    else
      if( myid == 0 )  call info( "station " // trim( stnm_g(i) ) // " is out of the region" )
    end if
  end do

  close( io_stlst )

  if( nst_g == 0 ) then
    sw_wav = .false.
    if( myid == 0 ) call info( "no station is detected. waveform file will not be created" )
    return
  end if

  if( nst == 0 ) return


  allocate( xst(nst), yst(nst), zst(nst), ist(nst), jst(nst), kst(nst), stnm(nst), stlo(nst), stla(nst) )


  do i = 1, nst

    ii = stid(i)

    xst(i) = xst_g(ii)
    yst(i) = yst_g(ii)
    zst(i) = zst_g(ii)
    ist(i) = ist_g(ii)
    jst(i) = jst_g(ii)
    ist(i) = ist_g(ii)
    stnm(i) = stnm_g(ii)
    stlo(i) = stlo_g(ii)
    stla(i) = stla_g(ii)

    !! station depth setting
    !! this is done only for MPI-node because of the definition of kfs and kob
    select case ( zsw_g(ii) )
    case( 'dep' );    kst(i) = kst_g(ii)
    case( 'fsb' );    kst(i) = kfs( ist(i), jst(i) ) + 1   !! free surface: one-grid below for avoiding vacuum
    case( 'obb' );    kst(i) = kob( ist(i), jst(i) ) + 1   !! ocean column: below seafloor
    case( 'oba' );    kst(i) = kob( ist(i), jst(i) ) - 1   !! ocean column: above seafloor
    case( 'bd0' );    kst(i) = z2k( bddep(ist(i),jst(i),0), zbeg, real(dz) ) !! Boundary interface
    case( 'bd1' );    kst(i) = z2k( bddep(ist(i),jst(i),1), zbeg, real(dz) ) !! Boundary interface
    case( 'bd2' );    kst(i) = z2k( bddep(ist(i),jst(i),2), zbeg, real(dz) ) !! Boundary interface
    case( 'bd3' );    kst(i) = z2k( bddep(ist(i),jst(i),3), zbeg, real(dz) ) !! Boundary interface
    case( 'bd4' );    kst(i) = z2k( bddep(ist(i),jst(i),4), zbeg, real(dz) ) !! Boundary interface
    case( 'bd5' );    kst(i) = z2k( bddep(ist(i),jst(i),5), zbeg, real(dz) ) !! Boundary interface
    case( 'bd6' );    kst(i) = z2k( bddep(ist(i),jst(i),6), zbeg, real(dz) ) !! Boundary interface
    case( 'bd7' );    kst(i) = z2k( bddep(ist(i),jst(i),7), zbeg, real(dz) ) !! Boundary interface
    case( 'bd8' );    kst(i) = z2k( bddep(ist(i),jst(i),8), zbeg, real(dz) ) !! Boundary interface
    case( 'bd9' );    kst(i) = z2k( bddep(ist(i),jst(i),9), zbeg, real(dz) ) !! Boundary interface
    case default
      call info( "unknown zsw type in station file. Assume 'dep'")
      kst(i) = kst_g(ii)
    end select

    !! depth check
    if( kst(i) > kend ) then
      call info( "station depth exceeds kend at station "//trim(stnm(i)) )
      kst(i) = kend - 1
    end if
    if( kst(i) < kbeg ) then
      write(STDERR,*) 'WARNING[output__setup]: station depth exceeds kbeg at station ' // trim(stnm(i))
      kst(i) = kbeg + 1
    end if
       

  end do

  allocate( vxst(ntw,nst), vyst(ntw,nst), vzst(ntw,nst), uxst(ntw,nst), uyst(ntw,nst), uzst(ntw,nst) )
  vxst(:,:) = 0.0
  vyst(:,:) = 0.0
  vzst(:,:) = 0.0
  uxst(:,:) = 0.0
  uyst(:,:) = 0.0
  uzst(:,:) = 0.0

  !!
  !! set-up sac header
  !!
  allocate( sh(6,nst) )
  do i=1, nst

    !! first initialize header type
    do j=1, 6
      call sac__init(sh(j,i))
    end do

    !! common header
    sh(:,i)%evlo    = evlo
    sh(:,i)%evla    = evla
    sh(:,i)%evdp    = evdp*1000
    sh(:,i)%tim     = exedate
    sh(:,i)%b       = tbeg
    sh(:,i)%delta   = ntdec_w * dt
    sh(:,i)%npts    = ntw
    sh(:,i)%mag     = mw
    if( bf_mode ) then
      sh(:,i)%user0   = fx0
      sh(:,i)%user1   = fy0
      sh(:,i)%user2   = fz0
    else
      sh(:,i)%user0   = mxx0
      sh(:,i)%user1   = myy0
      sh(:,i)%user2   = mzz0
      sh(:,i)%user3   = myz0
      sh(:,i)%user4   = mxz0
      sh(:,i)%user5   = mxy0
    end if
    sh(:,i)%user6   = clon !< coordinate
    sh(:,i)%user7   = clat !< coordinate
    sh(:,i)%user8   = phi
    sh(:,i)%o       = otim

    !! exedate -> time
    do j=1, 6
      call daytim__localtime( sh(j,i)%tim, &
      sh(j,i)%nzyear, sh(j,i)%nzmonth, sh(j,i)%nzday, sh(j,i)%nzhour, sh(j,i)%nzmin, sh(j,i)%nzsec )
      call daytim__ymd2jul  ( sh(j,i)%nzyear, sh(j,i)%nzmonth, sh(j,i)%nzday, sh(j,i)%nzjday )
      sh(j,i)%nzmsec = 0
    end do

    !! station dependent
    sh(:,i)%kevnm = trim(adjustl( title(1:16) ))
    sh(:,i)%kstnm = trim(stnm(i))
    sh(:,i)%stlo  = stlo(i)
    sh(:,i)%stla  = stla(i)
    sh(:,i)%stdp  = zst(i)*1000 ! in meter unit

    !! component dependent
    sh(1,i)%kcmpnm = "Vx"
    sh(2,i)%kcmpnm = "Vy"
    sh(3,i)%kcmpnm = "Vz"
    sh(4,i)%kcmpnm = "Ux"
    sh(5,i)%kcmpnm = "Uy"
    sh(6,i)%kcmpnm = "Uz"

    sh(1,i)%cmpinc = 90.0;  sh(1,i)%cmpaz  =  0.0 + phi
    sh(2,i)%cmpinc = 90.0;  sh(2,i)%cmpaz  = 90.0 + phi
    sh(3,i)%cmpinc =  0.0;  sh(3,i)%cmpaz  =  0.0
    sh(4,i)%cmpinc = 90.0;  sh(4,i)%cmpaz  =  0.0 + phi
    sh(5,i)%cmpinc = 90.0;  sh(5,i)%cmpaz  = 90.0 + phi
    sh(6,i)%cmpinc =  0.0;  sh(6,i)%cmpaz  =  0.0

    sh(1:3,i)%idep = 7 ! velocity [nm/s]
    sh(4:6,i)%idep = 6 ! displacement [nm]

  end do



end subroutine read_stinfo
!! --------------------------------------------------------------------------------------------------------------------------- !!




!! --------------------------------------------------------------------------------------------------------------------------- !!
!>
!! Open new file and write header information, and medium parameters for YZ-cross section
!<
!! ----
subroutine newfile_yz ( fname, hdr )

  character(*), intent(in)  :: fname
  type(snp),    intent(inout) :: hdr
  !! ----
  real(SP) :: buf(nys,nzs,3)
  integer :: j, k, kk, ii, jj
  !! ----

  ! absorber thickness setting
  hdr % na1 = na / jdec
  hdr % na2 = na / kdec

  if( myid == hdr % ionode ) then

    call std__getio( hdr%io, is_big=.true. )

#ifdef _ES
    open( hdr%io, file=trim(fname), action='write', status='replace', form='unformatted' )
#else
    open( hdr%io, file=trim(fname), access='stream', action='write', status='replace', form='unformatted' )
#endif
    call write_snp_header( hdr, nys, nzs, ysnp(1:nys), zsnp(1:nzs) )

  end if

  buf = 0.0

  ii = i0_yz
  if( ibeg <= ii .and. ii <= iend ) then
    do j = js0, js1
      do k = ks0, ks1
        kk = k * kdec - kdec/2
        jj = j * jdec - jdec/2

        buf(j,k,1) = rho( kk, ii, jj )
        buf(j,k,2) = lam( kk, ii, jj )
        buf(j,k,3) = mu ( kk, ii, jj )

      end do
    end do
  end if

  call write_reduce_array2d_r( nys, nzs, hdr % ionode, hdr % io, buf(:,:,1) )
  call write_reduce_array2d_r( nys, nzs, hdr % ionode, hdr % io, buf(:,:,2) )
  call write_reduce_array2d_r( nys, nzs, hdr % ionode, hdr % io, buf(:,:,3) )

end subroutine newfile_yz
!! --------------------------------------------------------------------------------------------------------------------------- !!


!! --------------------------------------------------------------------------------------------------------------------------- !!
!>
!! Open new file and write header information, and medium parameters for XY-cross section
!<
!! ----
subroutine newfile_xz ( fname, hdr )

  character(*), intent(in)  :: fname
  type(snp),    intent(inout) :: hdr
  !! --
  integer :: i, k, kk, ii, jj
  real(SP) :: buf(nxs,nzs,3)
  !! ----

  ! absorber thickness setting
  hdr % na1 = na / idec
  hdr % na2 = na / kdec

  if( myid == hdr % ionode ) then

    call std__getio( hdr % io, is_big=.true. )
#ifdef _ES
    open( hdr % io, file=trim(fname), action='write', status='replace', form='unformatted' )
#else
    open( hdr % io, file=trim(fname), access='stream', action='write', status='replace', form='unformatted' )
#endif
    call write_snp_header( hdr, nxs, nzs, xsnp(1:nxs), zsnp(1:nzs) )

  end if

  buf = 0.0

  jj = j0_xz
  if( jbeg <= jj .and. jj <= jend ) then
    do i = is0, is1
      do k = ks0, ks1

        ii = i * idec - idec/2
        kk = k * kdec - kdec/2

        buf(i,k,1) = rho( kk, ii, jj )
        buf(i,k,2) = lam( kk, ii, jj )
        buf(i,k,3) = mu ( kk, ii, jj )

      end do
    end do
  end if

  call write_reduce_array2d_r( nxs, nzs, hdr % ionode, hdr % io, buf(:,:,1) )
  call write_reduce_array2d_r( nxs, nzs, hdr % ionode, hdr % io, buf(:,:,2) )
  call write_reduce_array2d_r( nxs, nzs, hdr % ionode, hdr % io, buf(:,:,3) )

end subroutine newfile_xz
!! --------------------------------------------------------------------------------------------------------------------------- !!

!! --------------------------------------------------------------------------------------------------------------------------- !!
!>
!! Open new file and write header information, and medium parameters for XY-cross section
!<
!! ----
subroutine newfile_xy ( fname, hdr )

  character(*), intent(in) :: fname
  type(snp),    intent(inout) :: hdr
  !! --
  integer :: i, j, kk, ii, jj
  real(SP) :: buf(nxs,nys,3)
  !! ----

  ! absorber thickness setting
  hdr % na1 = na / idec
  hdr % na2 = na / jdec

  if( myid == hdr % ionode ) then

    call std__getio( hdr % io, is_big=.true. )
#ifdef _ES
    open( hdr % io, file=trim(fname), action='write', status='replace', form='unformatted' )
#else
    open( hdr % io, file=trim(fname), access='stream', action='write', status='replace', form='unformatted' )
#endif
    call write_snp_header( hdr, nxs, nys, xsnp(1:nxs), ysnp(1:nys) )

  end if

  buf = 0.0

  do j = js0, js1
    do i = is0, is1

      ii = i * idec - idec/2
      jj = j * jdec - jdec/2


      if( hdr % coordinate == 'fs' ) then
        kk = kfs(ii,jj) + 1
      else if ( hdr % coordinate == 'ob' ) then
        kk = kob(ii,jj) + 1
      else
        kk = k0_xy
      end if

      buf(i,j,1) = rho( kk, ii, jj )
      buf(i,j,2) = lam( kk, ii, jj )
      buf(i,j,3) = mu ( kk, ii, jj )

    end do
  end do

  call write_reduce_array2d_r( nxs, nys, hdr % ionode, hdr % io, buf(:,:,1) )
  call write_reduce_array2d_r( nxs, nys, hdr % ionode, hdr % io, buf(:,:,2) )
  call write_reduce_array2d_r( nxs, nys, hdr % ionode, hdr % io, buf(:,:,3) )

  do j = js0, js1
    do i = is0, is1

      ii = i * idec - idec/2
      jj = j * jdec - jdec/2
      !! topography data
      buf(i,j,1) = - bddep( ii, jj, 0 ) * 1000 ! positive upward, in unit of [m]

    end do
  end do
  call write_reduce_array2d_r( nxs, nys, hdr % ionode, hdr % io, buf(:,:,1) )


end subroutine newfile_xy
!! --------------------------------------------------------------------------------------------------------------------------- !!

!! --------------------------------------------------------------------------------------------------------------------------- !!
!>
!! create snapshot file in netcdf format
!<
!! --
subroutine newfile_yz_nc( fname, hdr )

  character(*), intent(in)    :: fname
  type(snp),    intent(inout) :: hdr
  !!
  real(SP), allocatable :: sbuf(:), rbuf1(:), rbuf2(:), rbuf3(:), buf(:,:,:)
  integer :: j, k, ii, jj, kk, ierr
  !! --

#ifdef _NETCDF

  if( myid == hdr%ionode ) then

    !! initialize
    hdr%vmax = 0.0
    hdr%vmin = 0.0
    hdr % na1 = na / jdec
    hdr % na2 = na / kdec
    hdr % ds1 = jdec * dy
    hdr % ds2 = kdec * dz

    call nc_chk( nf90_create( trim(fname), NF90_CLOBBER, hdr%io ) )
    call write_nc_header( hdr, nys, nzs, ysnp, zsnp )

  end if

  allocate( buf(nys,nzs,3) )

  buf = 0.0
  ii = i0_yz
  if( ibeg <= ii .and. ii <= iend ) then
    do j = js0, js1
      do k = ks0, ks1
        kk = k * kdec - kdec/2
        jj = j * jdec - jdec/2

        buf(j,k,1) = rho( kk, ii, jj )
        buf(j,k,2) = lam( kk, ii, jj )
        buf(j,k,3) = mu ( kk, ii, jj )

      end do
    end do
  end if

  !! medium
  allocate( sbuf(nys*nzs), rbuf1(nys*nzs), rbuf2(nys*nzs), rbuf3(nys*nzs) )
  sbuf = reshape( buf(:,:,1), shape(sbuf) )
  call mpi_reduce( sbuf, rbuf1, nys*nzs, MPI_REAL, MPI_SUM, hdr%ionode, mpi_comm_world, ierr )
  sbuf = reshape( buf(:,:,2), shape(sbuf) )
  call mpi_reduce( sbuf, rbuf2, nys*nzs, MPI_REAL, MPI_SUM, hdr%ionode, mpi_comm_world, ierr )
  sbuf = reshape( buf(:,:,3), shape(sbuf) )
  call mpi_reduce( sbuf, rbuf3, nys*nzs, MPI_REAL, MPI_SUM, hdr%ionode, mpi_comm_world, ierr )

  if( myid == hdr%ionode ) then

    call nc_chk( nf90_put_att( hdr%io, hdr%medid(1), 'actual_range', (/minval(rbuf1), maxval(rbuf1)/) ) )
    call nc_chk( nf90_put_att( hdr%io, hdr%medid(2), 'actual_range', (/minval(rbuf2), maxval(rbuf2)/) ) )
    call nc_chk( nf90_put_att( hdr%io, hdr%medid(3), 'actual_range', (/minval(rbuf3), maxval(rbuf3)/) ) )

    call nc_chk( nf90_enddef( hdr%io ) )

    call nc_chk( nf90_put_var( hdr%io, hdr%vid_x1, ysnp ) )
    call nc_chk( nf90_put_var( hdr%io, hdr%vid_x2, zsnp ) )
    call nc_chk( nf90_put_var( hdr%io, hdr%medid(1), reshape(rbuf1, shape(buf(:,:,1)) ) ) )
    call nc_chk( nf90_put_var( hdr%io, hdr%medid(2), reshape(rbuf2, shape(buf(:,:,1)) ) ) )
    call nc_chk( nf90_put_var( hdr%io, hdr%medid(3), reshape(rbuf3, shape(buf(:,:,1)) ) ) )

  end if

  deallocate( sbuf, rbuf1, rbuf2, rbuf3, buf )

#endif

end subroutine newfile_yz_nc
!! --------------------------------------------------------------------------------------------------------------------------- !!

!! --------------------------------------------------------------------------------------------------------------------------- !!
!>
!! create snapshot file in netcdf format
!<
!! --
subroutine newfile_xz_nc( fname, hdr )

  character(*), intent(in)    :: fname
  type(snp),    intent(inout) :: hdr
  !!
  real(SP), allocatable :: sbuf(:), rbuf1(:), rbuf2(:), rbuf3(:), buf(:,:,:)
  integer :: i, k, ii, jj, kk, ierr
  !! --

#ifdef _NETCDF

  if( myid == hdr%ionode ) then

    !! initialize
    hdr%vmax = 0.0
    hdr%vmin = 0.0
    hdr % na1 = na / idec
    hdr % na2 = na / kdec
    hdr % ds1 = idec * dx
    hdr % ds2 = kdec * dz

    call nc_chk( nf90_create( trim(fname), NF90_CLOBBER, hdr%io ) )
    call write_nc_header( hdr, nxs, nzs, xsnp, zsnp )

  end if

  allocate( buf(nxs,nzs,3) )

  buf = 0.0

  jj = j0_xz

  if( jbeg <= jj .and. jj <= jend ) then
    do i = is0, is1
      do k = ks0, ks1

        ii = i * idec - idec/2
        kk = k * kdec - kdec/2

        buf(i,k,1) = rho( kk, ii, jj )
        buf(i,k,2) = lam( kk, ii, jj )
        buf(i,k,3) = mu ( kk, ii, jj )

      end do
    end do
  end if

  !! medium
  allocate( sbuf(nxs*nzs), rbuf1(nxs*nzs), rbuf2(nxs*nzs), rbuf3(nxs*nzs) )
  sbuf = reshape( buf(:,:,1), shape(sbuf) )
  call mpi_reduce( sbuf, rbuf1, nxs*nzs, MPI_REAL, MPI_SUM, hdr%ionode, mpi_comm_world, ierr )
  sbuf = reshape( buf(:,:,2), shape(sbuf) )
  call mpi_reduce( sbuf, rbuf2, nxs*nzs, MPI_REAL, MPI_SUM, hdr%ionode, mpi_comm_world, ierr )
  sbuf = reshape( buf(:,:,3), shape(sbuf) )
  call mpi_reduce( sbuf, rbuf3, nxs*nzs, MPI_REAL, MPI_SUM, hdr%ionode, mpi_comm_world, ierr )

  if( myid == hdr%ionode ) then
    call nc_chk( nf90_put_att( hdr%io, hdr%medid(1), 'actual_range', (/minval(rbuf1), maxval(rbuf1)/) ) )
    call nc_chk( nf90_put_att( hdr%io, hdr%medid(2), 'actual_range', (/minval(rbuf2), maxval(rbuf2)/) ) )
    call nc_chk( nf90_put_att( hdr%io, hdr%medid(3), 'actual_range', (/minval(rbuf3), maxval(rbuf3)/) ) )

    call nc_chk( nf90_enddef( hdr%io ) )

    call nc_chk( nf90_put_var( hdr%io, hdr%vid_x1, xsnp ) )
    call nc_chk( nf90_put_var( hdr%io, hdr%vid_x2, zsnp ) )
    call nc_chk( nf90_put_var( hdr%io, hdr%medid(1), reshape(rbuf1, shape(buf(:,:,1)) ) ) )
    call nc_chk( nf90_put_var( hdr%io, hdr%medid(2), reshape(rbuf2, shape(buf(:,:,1)) ) ) )
    call nc_chk( nf90_put_var( hdr%io, hdr%medid(3), reshape(rbuf3, shape(buf(:,:,1)) ) ) )

  end if

  deallocate( sbuf, rbuf1, rbuf2, rbuf3, buf )

#endif

end subroutine newfile_xz_nc
!! --------------------------------------------------------------------------------------------------------------------------- !!

!! --------------------------------------------------------------------------------------------------------------------------- !!
!>
!! write common netcdf header
!<
!! ---
subroutine write_nc_header( hdr, ns1, ns2, xs1, xs2 )

  type(snp), intent(inout) :: hdr
  integer, intent(in) :: ns1, ns2
  real(SP), intent(in) :: xs1(ns1), xs2(ns2)
  integer :: i
  !! --
#ifdef _NETCDF
  if( hdr % coordinate == 'yz' ) then
    call nc_chk( nf90_def_dim( hdr%io, 'y', nys, hdr%did_x1 ) )
    call nc_chk( nf90_def_dim( hdr%io, 'z', nzs, hdr%did_x2 ) )
    call nc_chk( nf90_def_var( hdr%io, 'y', NF90_REAL, hdr%did_x1, hdr%vid_x1 ) )
    call nc_chk( nf90_def_var( hdr%io, 'z', NF90_REAL, hdr%did_x2, hdr%vid_x2 ) )
    call nc_chk( nf90_put_att( hdr%io, hdr%vid_x1, 'long_name', 'y' ) )
    call nc_chk( nf90_put_att( hdr%io, hdr%vid_x2, 'long_name', 'z' ) )
  else if( hdr % coordinate == 'xz' ) then
    call nc_chk( nf90_def_dim( hdr%io, 'x', nxs, hdr%did_x1 ) )
    call nc_chk( nf90_def_dim( hdr%io, 'z', nzs, hdr%did_x2 ) )
    call nc_chk( nf90_def_var( hdr%io, 'x', NF90_REAL, hdr%did_x1, hdr%vid_x1 ) )
    call nc_chk( nf90_def_var( hdr%io, 'z', NF90_REAL, hdr%did_x2, hdr%vid_x2 ) )
    call nc_chk( nf90_put_att( hdr%io, hdr%vid_x1, 'long_name', 'x' ) )
    call nc_chk( nf90_put_att( hdr%io, hdr%vid_x2, 'long_name', 'z' ) )
  else
    call nc_chk( nf90_def_dim( hdr%io, 'x', nxs, hdr%did_x1 ) )
    call nc_chk( nf90_def_dim( hdr%io, 'y', nys, hdr%did_x2 ) )
    call nc_chk( nf90_def_var( hdr%io, 'x', NF90_REAL, hdr%did_x1, hdr%vid_x1 ) )
    call nc_chk( nf90_def_var( hdr%io, 'y', NF90_REAL, hdr%did_x2, hdr%vid_x2 ) )
    call nc_chk( nf90_put_att( hdr%io, hdr%vid_x1, 'long_name', 'x' ) )
    call nc_chk( nf90_put_att( hdr%io, hdr%vid_x2, 'long_name', 'y' ) )
  end if

  call nc_chk( nf90_def_dim( hdr%io, 't', NF90_UNLIMITED, hdr%did_t ) )
  call nc_chk( nf90_def_var( hdr%io, 't', NF90_REAL, hdr%did_t, hdr%vid_t ) )

  !! medium
  call nc_chk( nf90_def_var( hdr%io, 'rho',    NF90_REAL, (/hdr%did_x1, hdr%did_x2/), hdr%medid(1) ) )
  call nc_chk( nf90_def_var( hdr%io, 'lambda', NF90_REAL, (/hdr%did_x1, hdr%did_x2/), hdr%medid(2) ) )
  call nc_chk( nf90_def_var( hdr%io, 'mu',     NF90_REAL, (/hdr%did_x1, hdr%did_x2/), hdr%medid(3) ) )

  !! special for horizontal snapshot
  if( hdr % coordinate == 'ob' .or.  hdr % coordinate == 'fs' .or.  hdr % coordinate == 'xy' ) then
    call nc_chk( nf90_def_var( hdr%io, 'topo',   NF90_REAL, (/hdr%did_x1, hdr%did_x2/), hdr%medid(4) ) )
    call nc_chk( nf90_def_var( hdr%io, 'lon',    NF90_REAL, (/hdr%did_x1, hdr%did_x2/), hdr%medid(5) ) )
    call nc_chk( nf90_def_var( hdr%io, 'lat',    NF90_REAL, (/hdr%did_x1, hdr%did_x2/), hdr%medid(6) ) )
    call nc_chk( nf90_put_att( hdr%io, hdr%medid(4), 'long_name', 'Topography'    ) )
    call nc_chk( nf90_put_att( hdr%io, hdr%medid(5), 'long_name', 'Longitude'    ) )
    call nc_chk( nf90_put_att( hdr%io, hdr%medid(6), 'long_name', 'Latitude'     ) )
    call nc_chk( nf90_put_att( hdr%io, hdr%medid(4), 'units', 'm' ) )
    call nc_chk( nf90_put_att( hdr%io, hdr%medid(5), 'units', 'degrees_east' ) )
    call nc_chk( nf90_put_att( hdr%io, hdr%medid(6), 'units', 'degrees_north' ) )

    !! for make tools recognize map projection
    call nc_chk( nf90_put_att( hdr%io, hdr%vid_x1, 'standard_name', 'projection_x_coordinate' ) )
    call nc_chk( nf90_put_att( hdr%io, hdr%vid_x2, 'standard_name', 'projection_y_coordinate' ) )
    do i=1, 4
      call nc_chk( nf90_put_att( hdr%io, hdr%medid(i), 'coordinates', 'lat lon' ) )
    end do
  end if

  !! variables
  do i=1, hdr%nsnp
    call nc_chk( nf90_def_var( hdr%io, trim(hdr%vname(i)), NF90_REAL, (/hdr%did_x1, hdr%did_x2, hdr%did_t/), hdr%varid(i)))
    call nc_chk( nf90_put_att( hdr%io, hdr%varid(i), 'long_name', trim(hdr%vname(i)) ) )
    call nc_chk( nf90_put_att( hdr%io, hdr%varid(i), 'units', trim(hdr%vunit(i)) ) )

    if( hdr % coordinate == 'ob' .or.  hdr % coordinate == 'fs' .or.  hdr % coordinate == 'xy' ) then
      call nc_chk( nf90_put_att( hdr%io, hdr%varid(i), 'coordinates', 'lat lon' ) )
    end if

  end do

  !! global attribute
  call nc_chk( nf90_put_att( hdr%io, NF90_GLOBAL, 'generated_by',   'SWPC' ) )
  call nc_chk( nf90_put_att( hdr%io, NF90_GLOBAL, 'codetype',       CODE_TYPE ) )
  call nc_chk( nf90_put_att( hdr%io, NF90_GLOBAL, 'hdrver',         HEADER_VERSION ) )
  call nc_chk( nf90_put_att( hdr%io, NF90_GLOBAL, 'title',          trim(title) ) )
  call nc_chk( nf90_put_att( hdr%io, NF90_GLOBAL, 'exedate',        exedate     ) )
  call nc_chk( nf90_put_att( hdr%io, NF90_GLOBAL, 'ns1',            ns1 ) )
  call nc_chk( nf90_put_att( hdr%io, NF90_GLOBAL, 'ns2',            ns2 ) )
  call nc_chk( nf90_put_att( hdr%io, NF90_GLOBAL, 'beg1',           xs1(1) ) )
  call nc_chk( nf90_put_att( hdr%io, NF90_GLOBAL, 'beg2',           xs2(1) ) )
  call nc_chk( nf90_put_att( hdr%io, NF90_GLOBAL, 'na1',            hdr%na1 ) )
  call nc_chk( nf90_put_att( hdr%io, NF90_GLOBAL, 'na2',            hdr%na2 ) )
  call nc_chk( nf90_put_att( hdr%io, NF90_GLOBAL, 'ds1',            hdr%ds1 ) )
  call nc_chk( nf90_put_att( hdr%io, NF90_GLOBAL, 'ds2',            hdr%ds2 ) )
  call nc_chk( nf90_put_att( hdr%io, NF90_GLOBAL, 'nmed',           hdr%nmed ) )
  call nc_chk( nf90_put_att( hdr%io, NF90_GLOBAL, 'nsnp',           hdr%nsnp ) )
  call nc_chk( nf90_put_att( hdr%io, NF90_GLOBAL, 'coordinate',     hdr%coordinate ) )
  call nc_chk( nf90_put_att( hdr%io, NF90_GLOBAL, 'datatype',       hdr%snaptype ) )
  call nc_chk( nf90_put_att( hdr%io, NF90_GLOBAL, 'dt',             dt * ntdec_s ) )


  !! variable attributes
  call nc_chk( nf90_put_att( hdr%io, hdr%vid_t, 'long_name', 't' ) )
  call nc_chk( nf90_put_att( hdr%io, hdr%vid_x1, 'units', 'km' ) )
  call nc_chk( nf90_put_att( hdr%io, hdr%vid_x2, 'units', 'km' ) )
  call nc_chk( nf90_put_att( hdr%io, hdr%vid_t, 'units', 's' ) )

  call nc_chk( nf90_put_att( hdr%io, hdr%medid(1), 'long_name', 'rho'    ) )
  call nc_chk( nf90_put_att( hdr%io, hdr%medid(2), 'long_name', 'lambda' ) )
  call nc_chk( nf90_put_att( hdr%io, hdr%medid(3), 'long_name', 'mu'     ) )
  call nc_chk( nf90_put_att( hdr%io, hdr%medid(1), 'units', '10^3 kg/cm^3'    ) )
  call nc_chk( nf90_put_att( hdr%io, hdr%medid(2), 'units', '10^9 Pa' ) )
  call nc_chk( nf90_put_att( hdr%io, hdr%medid(3), 'units', '10^9 Pa'     ) )

  call nc_chk( nf90_put_att( hdr%io, hdr%vid_x1, 'actual_range', (/ xs1(1), xs1(ns1) /) ) )
  call nc_chk( nf90_put_att( hdr%io, hdr%vid_x2, 'actual_range', (/ xs2(1), xs2(ns2) /) ) )

  call nc_chk( nf90_put_att( hdr%io, NF90_GLOBAL, 'clon', clon ) )
  call nc_chk( nf90_put_att( hdr%io, NF90_GLOBAL, 'clat', clat ) )
  call nc_chk( nf90_put_att( hdr%io, NF90_GLOBAL, 'phi',  phi  ) )

#endif

end subroutine write_nc_header
!! --------------------------------------------------------------------------------------------------------------------------- !!

!! --------------------------------------------------------------------------------------------------------------------------- !!
!>
!! create snapshot file in netcdf format
!<
!! --
subroutine newfile_xy_nc( fname, hdr )

  character(*), intent(in)    :: fname
  type(snp),    intent(inout) :: hdr
  !!
  real(SP), allocatable :: sbuf(:), rbuf1(:), rbuf2(:), rbuf3(:), rbuf4(:), rbuf5(:), rbuf6(:), buf(:,:,:)
  integer :: i, j, ii, jj, kk, ierr
  !! --

#ifdef _NETCDF

  if( myid == hdr%ionode ) then

    !! initialize
    hdr%vmax = 0.0
    hdr%vmin = 0.0
    hdr % na1 = na / idec
    hdr % na2 = na / jdec
    hdr % ds1 = idec * dx
    hdr % ds2 = jdec * dy

    call nc_chk( nf90_create( trim(fname), NF90_CLOBBER, hdr%io ) )
    call write_nc_header( hdr, nxs, nys, xsnp, ysnp )

  end if

  allocate( buf(nxs,nys,6) )


  buf = 0.0

  do j = js0, js1
    do i = is0, is1

      ii = i * idec - idec/2
      jj = j * jdec - jdec/2


      if( hdr % coordinate == 'fs' ) then
        kk = kfs(ii,jj) + 1
      else if ( hdr % coordinate == 'ob' ) then
        kk = kob(ii,jj) + 1
      else
        kk = k0_xy
      end if

      buf(i,j,1) = rho( kk, ii, jj )
      buf(i,j,2) = lam( kk, ii, jj )
      buf(i,j,3) = mu ( kk, ii, jj )
      buf(i,j,4) = - bddep( ii, jj, 0 ) * 1000 ! positive upward, in meter unit

      ! longitude & latitude table
      call geomap__c2g( xsnp(i), ysnp(j), clon, clat, phi, buf(i,j,5), buf(i,j,6) )
    end do
  end do

  !! medium
  allocate( sbuf(nxs*nys), rbuf1(nxs*nys), rbuf2(nxs*nys), rbuf3(nxs*nys) )
  allocate( rbuf4(nxs*nys), rbuf5(nxs*nys), rbuf6(nxs*nys) )
  sbuf = reshape( buf(:,:,1), shape(sbuf) )
  call mpi_reduce( sbuf, rbuf1, nxs*nys, MPI_REAL, MPI_SUM, hdr%ionode, mpi_comm_world, ierr )
  sbuf = reshape( buf(:,:,2), shape(sbuf) )
  call mpi_reduce( sbuf, rbuf2, nxs*nys, MPI_REAL, MPI_SUM, hdr%ionode, mpi_comm_world, ierr )
  sbuf = reshape( buf(:,:,3), shape(sbuf) )
  call mpi_reduce( sbuf, rbuf3, nxs*nys, MPI_REAL, MPI_SUM, hdr%ionode, mpi_comm_world, ierr )
  sbuf = reshape( buf(:,:,4), shape(sbuf) )
  call mpi_reduce( sbuf, rbuf4, nxs*nys, MPI_REAL, MPI_SUM, hdr%ionode, mpi_comm_world, ierr )
  sbuf = reshape( buf(:,:,5), shape(sbuf) )
  call mpi_reduce( sbuf, rbuf5, nxs*nys, MPI_REAL, MPI_SUM, hdr%ionode, mpi_comm_world, ierr )
  sbuf = reshape( buf(:,:,6), shape(sbuf) )
  call mpi_reduce( sbuf, rbuf6, nxs*nys, MPI_REAL, MPI_SUM, hdr%ionode, mpi_comm_world, ierr )

  if( myid == hdr%ionode ) then

    call nc_chk( nf90_put_att( hdr%io, hdr%medid(1), 'actual_range', (/minval(rbuf1), maxval(rbuf1)/) ) )
    call nc_chk( nf90_put_att( hdr%io, hdr%medid(2), 'actual_range', (/minval(rbuf2), maxval(rbuf2)/) ) )
    call nc_chk( nf90_put_att( hdr%io, hdr%medid(3), 'actual_range', (/minval(rbuf3), maxval(rbuf3)/) ) )
    call nc_chk( nf90_put_att( hdr%io, hdr%medid(4), 'actual_range', (/minval(rbuf4), maxval(rbuf4)/) ) )
    call nc_chk( nf90_put_att( hdr%io, hdr%medid(5), 'actual_range', (/minval(rbuf5), maxval(rbuf5)/) ) )
    call nc_chk( nf90_put_att( hdr%io, hdr%medid(6), 'actual_range', (/minval(rbuf6), maxval(rbuf6)/) ) )

    call nc_chk( nf90_enddef( hdr%io ) )

    call nc_chk( nf90_put_var( hdr%io, hdr%vid_x1, xsnp ) )
    call nc_chk( nf90_put_var( hdr%io, hdr%vid_x2, ysnp ) )
    call nc_chk( nf90_put_var( hdr%io, hdr%medid(1), reshape(rbuf1, shape(buf(:,:,1)) ) ) )
    call nc_chk( nf90_put_var( hdr%io, hdr%medid(2), reshape(rbuf2, shape(buf(:,:,1)) ) ) )
    call nc_chk( nf90_put_var( hdr%io, hdr%medid(3), reshape(rbuf3, shape(buf(:,:,1)) ) ) )
    call nc_chk( nf90_put_var( hdr%io, hdr%medid(4), reshape(rbuf4, shape(buf(:,:,1)) ) ) )
    call nc_chk( nf90_put_var( hdr%io, hdr%medid(5), reshape(rbuf5, shape(buf(:,:,1)) ) ) )
    call nc_chk( nf90_put_var( hdr%io, hdr%medid(6), reshape(rbuf6, shape(buf(:,:,1)) ) ) )

  end if

  deallocate( sbuf, rbuf1, rbuf2, rbuf3, rbuf4, rbuf5, rbuf6, buf )

#endif

end subroutine newfile_xy_nc
!! --------------------------------------------------------------------------------------------------------------------------- !!

!! --------------------------------------------------------------------------------------------------------------------------- !!
!>
!! write snapshot header in fixed format
!<
!! --------------------------------------------------------------------------------------------------------------------------- !!
subroutine write_snp_header( hdr, ns1, ns2, xs1, xs2 )

  type(snp),    intent(in) :: hdr
  integer,      intent(in) :: ns1, ns2
  real(SP),     intent(in) :: xs1(ns1), xs2(ns2) ! first and second axis
  real(SP) :: dum
  !! ----

  write( hdr % io ) BINARY_TYPE
  write( hdr % io ) CODE_TYPE
  write( hdr % io ) HEADER_VERSION
  write( hdr % io ) title
  write( hdr % io ) exedate

  !! space grid size
  write( hdr % io ) hdr % coordinate                ! coordinate
  write( hdr % io ) hdr % snaptype                  ! data type
  write( hdr % io ) ns1
  write( hdr % io ) ns2                             ! data size
  write( hdr % io ) xs1(1)
  write( hdr % io ) xs2(1)
  write( hdr % io ) xs1(2) - xs1(1)
  write( hdr % io ) xs2(2) - xs2(1)

  !! time
  write( hdr % io ) dt * ntdec_s                   ! dt

  ! absorb layer
  write( hdr % io ) hdr % na1
  write( hdr % io ) hdr % na2

  !! number of arrays
  write( hdr % io ) hdr % nmed
  write( hdr % io ) hdr % nsnp

  !! coordinate
  dum = 1.0
  write( hdr % io ) clon
  write( hdr % io ) clat
  write( hdr % io ) phi

  write( hdr % io ) dum
  write( hdr % io ) dum
  write( hdr % io ) dum

end subroutine write_snp_header
!! --------------------------------------------------------------------------------------------------------------------------- !!

!! --------------------------------------------------------------------------------------------------------------------------- !!
!>
!! write 2d array with MPI
!<
!! --------------------------------------------------------------------------------------------------------------------------- !!
subroutine write_reduce_array2d_r( nx1, nx2, ionode, io, array )

  integer,  intent(in) :: nx1, nx2
  integer,  intent(in) :: ionode
  integer,  intent(in) :: io
  real(SP), intent(in) :: array(nx1,nx2)

  real(SP) :: sbuf(nx1*nx2)
  real(SP) :: rbuf(nx1*nx2)
  integer  :: ierr

  !! prepare send buffer
  sbuf = reshape( array, shape(sbuf) )

  !! gather to io_node
  call mpi_reduce( sbuf, rbuf, nx1*nx2, MPI_REAL, MPI_SUM, ionode, mpi_comm_world, ierr )

  !! write
  if( myid == ionode ) write(io) rbuf

end subroutine write_reduce_array2d_r
!! --------------------------------------------------------------------------------------------------------------------------- !!

!! --------------------------------------------------------------------------------------------------------------------------- !!
!>
!! write 2d array with MPI
!<
!! --------------------------------------------------------------------------------------------------------------------------- !!
subroutine write_reduce_array2d_r_nc( it, vid, nx1, nx2, hdr, array )

  integer,   intent(in) :: it
  integer,   intent(in) :: vid
  integer,   intent(in) :: nx1, nx2
  type(snp), intent(inout) :: hdr
  real(SP),  intent(in) :: array(nx1,nx2)

  real(SP) :: sbuf(nx1*nx2)
  real(SP) :: rbuf(nx1*nx2)
  integer  :: ierr
  integer :: count(3)
  integer :: start(3)

#ifdef _NETCDF
  !! prepare send buffer
  sbuf = reshape( array, shape(sbuf) )

  !! gather to io_node
  call mpi_reduce( sbuf, rbuf, nx1*nx2, MPI_REAL, MPI_SUM, hdr%ionode, mpi_comm_world, ierr )

  !! write
  if( myid == hdr%ionode ) then
    count = (/ nx1, nx2, 1/)
    start = (/ 1, 1, it/ntdec_s+1 /)
    call nc_chk( nf90_put_var( hdr%io, hdr%varid(vid), reshape(rbuf,shape(array)), start, count ))
    call nc_chk( nf90_put_var( hdr%io, hdr%vid_t, it*dt, start=(/ it/ntdec_s+1 /) ) )
    hdr%vmax(vid) = max( hdr%vmax(vid), maxval(rbuf) )
    hdr%vmin(vid) = min( hdr%vmin(vid), minval(rbuf) )
    call nc_chk( nf90_redef( hdr%io ) )
    call nc_chk( nf90_put_att( hdr%io, hdr%varid(vid), 'actual_range', (/hdr%vmin(vid), hdr%vmax(vid)/)) )
    call nc_chk( nf90_enddef( hdr%io ) )
  end if
#endif

end subroutine write_reduce_array2d_r_nc
!! --------------------------------------------------------------------------------------------------------------------------- !!

!! --------------------------------------------------------------------------------------------------------------------------- !!
!>
!! write 1d array with MPI
!<
!! --------------------------------------------------------------------------------------------------------------------------- !!
subroutine write_reduce_array1d_r( nx, ionode, io, array )

  integer,  intent(in) :: nx
  integer,  intent(in) :: ionode
  integer,  intent(in) :: io
  real(SP), intent(in) :: array(nx)

  real(SP) :: rbuf(nx)
  integer  :: ierr

  !! gather to io_node
  call mpi_reduce( array, rbuf, nx, MPI_REAL, MPI_SUM, ionode, mpi_comm_world, ierr )

  !! write
  if( myid == ionode ) write(io) rbuf

end subroutine write_reduce_array1d_r
!! --------------------------------------------------------------------------------------------------------------------------- !!

!! --------------------------------------------------------------------------------------------------------------------------- !!
!>
!! Calculate divergence and abs(rotation) of velocity vector at given location (k,i,j) by 2nd order FDM
!<
!! ----
subroutine divrot( k, i, j, div, rot )
  !! --
  integer,  intent(in)  :: k, i, j
  real(SP), intent(out) :: div, rot(3)
  !!

  real(SP) :: dxVx, dxVy, dxVz, dyVx, dyVy, dyVz, dzVx, dzVy, dzVz

  dxVx = (  Vx(k  ,i  ,j  ) - Vx(k  ,i-1,j  )  ) * r20x
  dxVy = (  Vy(k  ,i+1,j  ) - Vy(k  ,i  ,j  )  ) * r20x
  dxVz = (  Vz(k  ,i+1,j  ) - Vz(k  ,i  ,j  )  ) * r20x
  dyVx = (  Vx(k  ,i  ,j+1) - Vx(k  ,i  ,j  )  ) * r20y
  dyVy = (  Vy(k  ,i  ,j  ) - Vy(k  ,i  ,j-1)  ) * r20y
  dyVz = (  Vz(k  ,i  ,j+1) - Vz(k  ,i  ,j  )  ) * r20y
  dzVx = (  Vx(k+1,i  ,j  ) - Vx(k  ,i  ,j  )  ) * r20z
  dzVy = (  Vy(k+1,i  ,j  ) - Vy(k  ,i  ,j  )  ) * r20z
  dzVz = (  Vz(k  ,i  ,j  ) - Vz(k-1,i  ,j  )  ) * r20z

  div = dxVx + dyVy + dzVz
  rot(1) = dyVz - dzVy  ! x
  rot(2) = dzVx - dxVz  ! y
  rot(3) = dxVy - dyVx  ! z

  !! masking
  div = div * lam(k,i,j) / abs( lam(k,i,j) + epsilon(1.0) )
  rot(1) = rot(1) * muyz(k,i,j) / abs( muyz(k,i,j) + epsilon(1.0) )
  rot(2) = rot(2) * muxz(k,i,j) / abs( muxz(k,i,j) + epsilon(1.0) )
  rot(3) = rot(3) * muxy(k,i,j) / abs( muxy(k,i,j) + epsilon(1.0) )

end subroutine divrot
!! --------------------------------------------------------------------------------------------------------------------------- !!

!! --------------------------------------------------------------------------------------------------------------------------- !!
!>
!! write snapshot
!<
!! ----
subroutine output__write_snap( it )

  integer, intent(in) :: it


  call pwatch__on( "output__write_snap" )


  !! 個別処理

  if( yz_ps%sw ) call wbuf_yz_ps(it)
  if( xz_ps%sw ) call wbuf_xz_ps(it)
  if( xy_ps%sw ) call wbuf_xy_ps(it)
  if( fs_ps%sw ) call wbuf_fs_ps(it)
  if( ob_ps%sw ) call wbuf_ob_ps(it)

  if( yz_v%sw ) call wbuf_yz_v(it)
  if( xz_v%sw ) call wbuf_xz_v(it)
  if( xy_v%sw ) call wbuf_xy_v(it)
  if( fs_v%sw ) call wbuf_fs_v(it)
  if( ob_v%sw ) call wbuf_ob_v(it)

  if( yz_u%sw ) call wbuf_yz_u(it)
  if( xz_u%sw ) call wbuf_xz_u(it)
  if( xy_u%sw ) call wbuf_xy_u(it)
  if( fs_u%sw ) call wbuf_fs_u(it)
  if( ob_u%sw ) call wbuf_ob_u(it)


  call pwatch__off( "output__write_snap" )

end subroutine output__write_snap
!! --------------------------------------------------------------------------------------------------------------------------- !!




!! --------------------------------------------------------------------------------------------------------------------------- !!
subroutine wbuf_yz_ps(it)

  integer, intent(in) :: it
  integer :: i, j, k
  integer :: jj, kk
  real(SP) :: div, rot(3)
  real, allocatable :: buf(:,:,:)
  !! ----

  if( .not. allocated(buf) ) then
    allocate( buf(nys,nzs,4) )
    buf(1:nys,1:nzs,1:4) = 0.0
  end if

  if( mod( it-1, ntdec_s ) == 0 ) then

    i = i0_yz

    if( ibeg <= i .and. i <= iend ) then
      !$omp parallel do private( k, j, kk, jj, div, rot )
      do jj = js0, js1
        do kk= ks0, ks1
          k = kk * kdec - kdec/2
          j = jj * jdec - jdec/2

          call divrot( k, i, j, div, rot )
          !! dx, dy, dz have km unit. correction for 1e3 factor.
          buf(jj,kk,1) = div    * UC * M0 * 1e-3
          buf(jj,kk,2) = rot(1) * UC * M0 * 1e-3
          buf(jj,kk,3) = rot(2) * UC * M0 * 1e-3
          buf(jj,kk,4) = rot(3) * UC * M0 * 1e-3

        end do
      end do
      !$omp end parallel do
    end if

    if( snp_format == 'native' ) then
      call write_reduce_array2d_r( nys, nzs, yz_ps%ionode, yz_ps%io, buf(:,:,1) )
      call write_reduce_array2d_r( nys, nzs, yz_ps%ionode, yz_ps%io, buf(:,:,2) )
      call write_reduce_array2d_r( nys, nzs, yz_ps%ionode, yz_ps%io, buf(:,:,3) )
      call write_reduce_array2d_r( nys, nzs, yz_ps%ionode, yz_ps%io, buf(:,:,4) )
    else
      call write_reduce_array2d_r_nc ( it, 1, nys, nzs, yz_ps, buf(:,:,1) )
      call write_reduce_array2d_r_nc ( it, 2, nys, nzs, yz_ps, buf(:,:,2) )
      call write_reduce_array2d_r_nc ( it, 3, nys, nzs, yz_ps, buf(:,:,3) )
      call write_reduce_array2d_r_nc ( it, 4, nys, nzs, yz_ps, buf(:,:,4) )
    end if

  end if

end subroutine wbuf_yz_ps
!! --------------------------------------------------------------------------------------------------------------------------- !!

!! --------------------------------------------------------------------------------------------------------------------------- !!
subroutine wbuf_xz_ps(it)

  integer, intent(in) :: it
  integer :: i, j, k
  integer :: ii, kk
  real(SP) :: div, rot(3)
  real, allocatable :: buf(:,:,:)
  !! ----

  if( .not. allocated(buf) ) then
    allocate( buf(nxs,nzs,4) )
    buf(1:nxs,1:nzs,1:4) = 0.0
  end if

  if( mod( it-1, ntdec_s ) == 0 ) then

    j = j0_xz

    if( jbeg <= j .and. j <= jend ) then
      !$omp parallel do private( ii, kk, i, k, div, rot )
      do ii = is0, is1
        do kk= ks0, ks1
          k = kk * kdec - kdec/2
          i = ii * idec - idec/2

          call divrot( k, i, j, div, rot )
          !! dx, dy, dz have km unit. correction for 1e3 factor.
          buf(ii,kk,1) = div     * UC * M0 * 1e-3
          buf(ii,kk,2) = rot(1)  * UC * M0 * 1e-3
          buf(ii,kk,3) = rot(2)  * UC * M0 * 1e-3
          buf(ii,kk,4) = rot(3)  * UC * M0 * 1e-3

        end do
      end do
      !$omp end parallel do
    end if

    if( snp_format == 'native' ) then
      call write_reduce_array2d_r( nxs, nzs, xz_ps%ionode, xz_ps%io, buf(:,:,1) )
      call write_reduce_array2d_r( nxs, nzs, xz_ps%ionode, xz_ps%io, buf(:,:,2) )
      call write_reduce_array2d_r( nxs, nzs, xz_ps%ionode, xz_ps%io, buf(:,:,3) )
      call write_reduce_array2d_r( nxs, nzs, xz_ps%ionode, xz_ps%io, buf(:,:,4) )
    else
      call write_reduce_array2d_r_nc( it, 1, nxs, nzs, xz_ps, buf(:,:,1) )
      call write_reduce_array2d_r_nc( it, 2, nxs, nzs, xz_ps, buf(:,:,2) )
      call write_reduce_array2d_r_nc( it, 3, nxs, nzs, xz_ps, buf(:,:,3) )
      call write_reduce_array2d_r_nc( it, 4, nxs, nzs, xz_ps, buf(:,:,4) )
    end if

  end if

end subroutine wbuf_xz_ps
!! --------------------------------------------------------------------------------------------------------------------------- !!

!! --------------------------------------------------------------------------------------------------------------------------- !!
subroutine wbuf_xy_ps(it)

  integer, intent(in) :: it
  integer :: i, j, k
  integer :: ii, jj
  real(SP) :: div, rot(3)
  real, allocatable :: buf(:,:,:)
  !! ----

  if( .not. allocated(buf) ) then
    allocate( buf(nxs,nys,4) )
    buf(1:nxs,1:nys,1:4) = 0.0
  end if

  if( mod( it-1, ntdec_s ) == 0 ) then

    k = k0_xy
    !$omp parallel do private( ii, jj, i, j, div, rot )
    do jj = js0, js1
      do ii = is0, is1
        j = jj * jdec - jdec/2
        i = ii * idec - idec/2

        call divrot( k, i, j, div, rot )
        !! dx, dy, dz have km unit. correction for 1e3 factor.
        buf(ii,jj,1) = div    * UC * M0 * 1e-3
        buf(ii,jj,2) = rot(1) * UC * M0 * 1e-3
        buf(ii,jj,3) = rot(2) * UC * M0 * 1e-3
        buf(ii,jj,4) = rot(3) * UC * M0 * 1e-3

      end do
    end do
    !$omp end parallel do

    if( snp_format == 'native' ) then
      call write_reduce_array2d_r( nxs, nys, xy_ps%ionode, xy_ps%io, buf(:,:,1) )
      call write_reduce_array2d_r( nxs, nys, xy_ps%ionode, xy_ps%io, buf(:,:,2) )
      call write_reduce_array2d_r( nxs, nys, xy_ps%ionode, xy_ps%io, buf(:,:,3) )
      call write_reduce_array2d_r( nxs, nys, xy_ps%ionode, xy_ps%io, buf(:,:,4) )
    else
      call write_reduce_array2d_r_nc( it, 1, nxs, nys, xy_ps, buf(:,:,1) )
      call write_reduce_array2d_r_nc( it, 2, nxs, nys, xy_ps, buf(:,:,2) )
      call write_reduce_array2d_r_nc( it, 3, nxs, nys, xy_ps, buf(:,:,3) )
      call write_reduce_array2d_r_nc( it, 4, nxs, nys, xy_ps, buf(:,:,4) )
    end if

  end if

end subroutine wbuf_xy_ps
!! --------------------------------------------------------------------------------------------------------------------------- !!
!! --------------------------------------------------------------------------------------------------------------------------- !!
subroutine wbuf_fs_ps(it)

  integer, intent(in) :: it
  integer :: i, j, k
  integer :: ii, jj
  real(SP) :: div, rot(3)
  real, allocatable :: buf(:,:,:)
  !! ----

  if( .not. allocated(buf) ) then
    allocate( buf(nxs,nys,4) )
    buf(1:nxs,1:nys,1:4) = 0.0
  end if

  if( mod( it-1, ntdec_s ) == 0 ) then

    !$omp parallel do private( ii, jj, i, j, k, div, rot )
    do jj = js0, js1
      do ii = is0, is1
        j = jj * jdec - jdec/2
        i = ii * idec - idec/2
        k = kfs(i,j) + 1        !! to calculate derivatives in depth, to assure amplitude exist at detpth

        call divrot( k, i, j, div, rot )
        !! dx, dy, dz have km unit. correction for 1e3 factor.
        buf(ii,jj,1) = div    * UC * M0 * 1e-3
        buf(ii,jj,2) = rot(1) * UC * M0 * 1e-3
        buf(ii,jj,3) = rot(2) * UC * M0 * 1e-3
        buf(ii,jj,4) = rot(3) * UC * M0 * 1e-3
      end do
    end do
    !$omp end parallel do

    if( snp_format == 'native' ) then
      call write_reduce_array2d_r( nxs, nys, fs_ps%ionode, fs_ps%io, buf(:,:,1) )
      call write_reduce_array2d_r( nxs, nys, fs_ps%ionode, fs_ps%io, buf(:,:,2) )
      call write_reduce_array2d_r( nxs, nys, fs_ps%ionode, fs_ps%io, buf(:,:,3) )
      call write_reduce_array2d_r( nxs, nys, fs_ps%ionode, fs_ps%io, buf(:,:,4) )
    else
      call write_reduce_array2d_r_nc( it, 1, nxs, nys, fs_ps, buf(:,:,1) )
      call write_reduce_array2d_r_nc( it, 2, nxs, nys, fs_ps, buf(:,:,2) )
      call write_reduce_array2d_r_nc( it, 3, nxs, nys, fs_ps, buf(:,:,3) )
      call write_reduce_array2d_r_nc( it, 4, nxs, nys, fs_ps, buf(:,:,4) )
    end if

  end if

end subroutine wbuf_fs_ps
!! --------------------------------------------------------------------------------------------------------------------------- !!

!! --------------------------------------------------------------------------------------------------------------------------- !!
subroutine wbuf_ob_ps(it)

  integer, intent(in) :: it
  integer :: i, j, k
  integer :: ii, jj
  real(SP) :: div, rot(3)
  real, allocatable :: buf(:,:,:)
  !! ----

  if( .not. allocated(buf) ) then
    allocate( buf(nxs,nys,4) )
    buf(1:nxs,1:nys,1:4) = 0.0
  end if

  if( mod( it-1, ntdec_s ) == 0 ) then

    !$omp parallel do private( ii, jj, i, j, k, div, rot )
    do jj = js0, js1
      do ii = is0, is1
        j = jj * jdec - jdec/2
        i = ii * idec - idec/2
        k = kob(i,j) + 1        !! to calculate derivatives in depth, to assure amplitude exist at detpth

        call divrot( k, i, j, div, rot )
        !! dx, dy, dz have km unit. correction for 1e3 factor.
        buf(ii,jj,1) = div    * UC * M0 * 1e-3
        buf(ii,jj,2) = rot(1) * UC * M0 * 1e-3
        buf(ii,jj,3) = rot(2) * UC * M0 * 1e-3
        buf(ii,jj,4) = rot(3) * UC * M0 * 1e-3

      end do
    end do
    !$omp end parallel do

    if( snp_format == 'native' ) then
      call write_reduce_array2d_r( nxs, nys, ob_ps%ionode, ob_ps%io, buf(:,:,1) )
      call write_reduce_array2d_r( nxs, nys, ob_ps%ionode, ob_ps%io, buf(:,:,2) )
      call write_reduce_array2d_r( nxs, nys, ob_ps%ionode, ob_ps%io, buf(:,:,3) )
      call write_reduce_array2d_r( nxs, nys, ob_ps%ionode, ob_ps%io, buf(:,:,4) )
    else
      call write_reduce_array2d_r_nc( it, 1, nxs, nys, ob_ps, buf(:,:,1) )
      call write_reduce_array2d_r_nc( it, 2, nxs, nys, ob_ps, buf(:,:,2) )
      call write_reduce_array2d_r_nc( it, 3, nxs, nys, ob_ps, buf(:,:,3) )
      call write_reduce_array2d_r_nc( it, 4, nxs, nys, ob_ps, buf(:,:,4) )
    end if

  end if

end subroutine wbuf_ob_ps
!! --------------------------------------------------------------------------------------------------------------------------- !!

!! --------------------------------------------------------------------------------------------------------------------------- !!
subroutine wbuf_yz_v(it)
  integer, intent(in) :: it
  integer :: i, j, k
  integer :: jj, kk
  real(SP), allocatable :: buf(:,:,:)
  !! --

  if( mod( it-1, ntdec_s ) /= 0 ) return

  if( .not. allocated(buf) ) allocate(buf(nys,nzs,3))

  i = i0_yz

  buf = 0.0
  if( ibeg <= i .and. i <= iend ) then
    !$omp parallel do private( jj, kk, k, j )
    do jj = js0, js1
      do kk= ks0, ks1
        k = kk * kdec - kdec/2
        j = jj * jdec - jdec/2

        buf(jj,kk,1) = Vx(k,i,j) * UC * M0
        buf(jj,kk,2) = Vy(k,i,j) * UC * M0
        buf(jj,kk,3) = Vz(k,i,j) * UC * M0

      end do
    end do
    !$omp end parallel do

  end if

  if( snp_format == 'native' ) then
    call write_reduce_array2d_r( nys, nzs, yz_v%ionode, yz_v%io, buf(:,:,1) )
    call write_reduce_array2d_r( nys, nzs, yz_v%ionode, yz_v%io, buf(:,:,2) )
    call write_reduce_array2d_r( nys, nzs, yz_v%ionode, yz_v%io, buf(:,:,3) )
  else
    call write_reduce_array2d_r_nc( it, 1, nys, nzs, yz_v, buf(:,:,1) )
    call write_reduce_array2d_r_nc( it, 2, nys, nzs, yz_v, buf(:,:,2) )
    call write_reduce_array2d_r_nc( it, 3, nys, nzs, yz_v, buf(:,:,3) )
  end if

end subroutine wbuf_yz_v
!! --------------------------------------------------------------------------------------------------------------------------- !!

!! --------------------------------------------------------------------------------------------------------------------------- !!
subroutine wbuf_xz_v(it)

  integer, intent(in) :: it
  integer :: i, j, k
  integer :: ii, kk
  real(SP), allocatable :: buf(:,:,:)

  if( mod( it-1, ntdec_s ) /= 0 ) return
  if( .not. allocated(buf) )  allocate( buf(nxs,nzs,3) )

  j = j0_xz

  buf = 0.0
  if( jbeg <= j .and. j <= jend ) then
    !$omp parallel do private( ii, kk, i, k )
    do ii = is0, is1
      do kk= ks0, ks1
        k = kk * kdec - kdec/2
        i = ii * idec - idec/2

        buf(ii,kk,1) = Vx(k,i,j) * UC * M0
        buf(ii,kk,2) = Vy(k,i,j) * UC * M0
        buf(ii,kk,3) = Vz(k,i,j) * UC * M0

      end do
    end do
    !$omp end parallel do

  end if

  if( snp_format == 'native' ) then
    call write_reduce_array2d_r( nxs, nzs, xz_v%ionode, xz_v%io, buf(:,:,1) )
    call write_reduce_array2d_r( nxs, nzs, xz_v%ionode, xz_v%io, buf(:,:,2) )
    call write_reduce_array2d_r( nxs, nzs, xz_v%ionode, xz_v%io, buf(:,:,3) )
  else
    call write_reduce_array2d_r_nc( it, 1, nxs, nzs, xz_v, buf(:,:,1) )
    call write_reduce_array2d_r_nc( it, 2, nxs, nzs, xz_v, buf(:,:,2) )
    call write_reduce_array2d_r_nc( it, 3, nxs, nzs, xz_v, buf(:,:,3) )
  end if

end subroutine wbuf_xz_v
!! --------------------------------------------------------------------------------------------------------------------------- !!

!! --------------------------------------------------------------------------------------------------------------------------- !!
subroutine wbuf_xy_v(it)
  integer, intent(in) :: it
  integer :: i, j, k
  integer :: ii, jj
  real(SP), allocatable :: buf(:,:,:)
  !! ---

  if( mod( it-1, ntdec_s ) /= 0 ) return
  if( .not. allocated(buf) ) allocate( buf(nx,ny,3) )

  k = k0_xy
  buf = 0.0

  !$omp parallel do private( ii, jj, i, j )
  do jj = js0, js1
    do ii = is0, is1
      j = jj * jdec - jdec/2
      i = ii * idec - idec/2

      buf(ii,jj,1) = Vx(k,i,j) * UC * M0
      buf(ii,jj,2) = Vy(k,i,j) * UC * M0
      buf(ii,jj,3) = Vz(k,i,j) * UC * M0

    end do
  end do
  !$omp end parallel do

  if( snp_format == 'native' ) then
    call write_reduce_array2d_r( nxs, nys, xy_v%ionode, xy_v%io, buf(:,:,1) )
    call write_reduce_array2d_r( nxs, nys, xy_v%ionode, xy_v%io, buf(:,:,2) )
    call write_reduce_array2d_r( nxs, nys, xy_v%ionode, xy_v%io, buf(:,:,3) )
  else
    call write_reduce_array2d_r_nc( it, 1, nxs, nys, xy_v, buf(:,:,1) )
    call write_reduce_array2d_r_nc( it, 2, nxs, nys, xy_v, buf(:,:,2) )
    call write_reduce_array2d_r_nc( it, 3, nxs, nys, xy_v, buf(:,:,3) )
  end if

end subroutine wbuf_xy_v
!! --------------------------------------------------------------------------------------------------------------------------- !!

!! --------------------------------------------------------------------------------------------------------------------------- !!
subroutine wbuf_fs_v(it)

  integer, intent(in) :: it
  integer :: i, j, k
  integer :: ii, jj
  real(SP), allocatable :: buf(:,:,:)

  if( mod( it-1, ntdec_s ) /= 0 ) return
  if( .not. allocated(buf) ) allocate(buf(nxs,nys,3))
  buf = 0.0

  !$omp parallel do private( ii, jj, i, j, k )
  do jj = js0, js1
    do ii = is0, is1
      j = jj * jdec - jdec/2
      i = ii * idec - idec/2
      k = kfs(i,j) + 1

      buf(ii,jj,1) = Vx(k,i,j) * UC * M0
      buf(ii,jj,2) = Vy(k,i,j) * UC * M0
      buf(ii,jj,3) = Vz(k,i,j) * UC * M0

    end do
  end do
  !$omp end parallel do

  if( snp_format == 'native' ) then
    call write_reduce_array2d_r( nxs, nys, fs_v%ionode, fs_v%io, buf(:,:,1) )
    call write_reduce_array2d_r( nxs, nys, fs_v%ionode, fs_v%io, buf(:,:,2) )
    call write_reduce_array2d_r( nxs, nys, fs_v%ionode, fs_v%io, buf(:,:,3) )
  else
    call write_reduce_array2d_r_nc( it, 1, nxs, nys, fs_v, buf(:,:,1) )
    call write_reduce_array2d_r_nc( it, 2, nxs, nys, fs_v, buf(:,:,2) )
    call write_reduce_array2d_r_nc( it, 3, nxs, nys, fs_v, buf(:,:,3) )
  end if

end subroutine wbuf_fs_v
!! --------------------------------------------------------------------------------------------------------------------------- !!

!! --------------------------------------------------------------------------------------------------------------------------- !!
subroutine wbuf_ob_v(it)

  integer, intent(in) :: it
  integer :: i, j, k
  integer :: ii, jj
  real(SP), allocatable :: buf(:,:,:)

  if( mod( it-1, ntdec_s ) /= 0 ) return
  if( .not. allocated(buf) ) allocate( buf(nxs,nys,3) )

  buf = 0.0

  !$omp parallel do private( ii, jj, i, j, k )
  do jj = js0, js1
    do ii = is0, is1
      j = jj * jdec - jdec/2
      i = ii * idec - idec/2
      k = kob(i,j) + 1

      buf(ii,jj,1) = Vx(k,i,j) * UC * M0
      buf(ii,jj,2) = Vy(k,i,j) * UC * M0
      buf(ii,jj,3) = Vz(k,i,j) * UC * M0

    end do
  end do
  !$omp end parallel do

  if( snp_format == 'native' ) then
    call write_reduce_array2d_r( nxs, nys, ob_v%ionode, ob_v%io, buf(:,:,1) )
    call write_reduce_array2d_r( nxs, nys, ob_v%ionode, ob_v%io, buf(:,:,2) )
    call write_reduce_array2d_r( nxs, nys, ob_v%ionode, ob_v%io, buf(:,:,3) )
  else
    call write_reduce_array2d_r_nc( it, 1, nxs, nys, ob_v, buf(:,:,1) )
    call write_reduce_array2d_r_nc( it, 2, nxs, nys, ob_v, buf(:,:,2) )
    call write_reduce_array2d_r_nc( it, 3, nxs, nys, ob_v, buf(:,:,3) )
  end if


end subroutine wbuf_ob_v
!! --------------------------------------------------------------------------------------------------------------------------- !!

!! --------------------------------------------------------------------------------------------------------------------------- !!
subroutine wbuf_yz_u(it)
  integer, intent(in) :: it
  integer :: i, j, k
  integer :: jj, kk
  !! --

  i = i0_yz

  if( ibeg <= i .and. i <= iend ) then
    !$omp parallel do private( jj, kk, k, j )
    do jj = js0, js1
      do kk= ks0, ks1
        k = kk * kdec - kdec/2
        j = jj * jdec - jdec/2

        buf_yz_u(jj,kk,1) = buf_yz_u(jj,kk,1) + Vx(k,i,j) * UC * M0 * dt
        buf_yz_u(jj,kk,2) = buf_yz_u(jj,kk,2) + Vy(k,i,j) * UC * M0 * dt
        buf_yz_u(jj,kk,3) = buf_yz_u(jj,kk,3) + Vz(k,i,j) * UC * M0 * dt

      end do
    end do
    !$omp end parallel do

  end if

  if( mod( it-1, ntdec_s ) == 0 ) then
    if( snp_format == 'native' ) then
      call write_reduce_array2d_r( nys, nzs, yz_u%ionode, yz_u%io, buf_yz_u(:,:,1) )
      call write_reduce_array2d_r( nys, nzs, yz_u%ionode, yz_u%io, buf_yz_u(:,:,2) )
      call write_reduce_array2d_r( nys, nzs, yz_u%ionode, yz_u%io, buf_yz_u(:,:,3) )
    else
      call write_reduce_array2d_r_nc( it, 1, nys, nzs, yz_u, buf_yz_u(:,:,1) )
      call write_reduce_array2d_r_nc( it, 2, nys, nzs, yz_u, buf_yz_u(:,:,2) )
      call write_reduce_array2d_r_nc( it, 3, nys, nzs, yz_u, buf_yz_u(:,:,3) )
    end if
  end if

end subroutine wbuf_yz_u
!! --------------------------------------------------------------------------------------------------------------------------- !!

!! --------------------------------------------------------------------------------------------------------------------------- !!
subroutine wbuf_xz_u(it)

  integer, intent(in) :: it
  integer :: i, j, k
  integer :: ii, kk

  j = j0_xz

  if( jbeg <= j .and. j <= jend ) then
    !$omp parallel do private( ii, kk, i, k )
    do ii = is0, is1
      do kk= ks0, ks1
        k = kk * kdec - kdec/2
        i = ii * idec - idec/2

        buf_xz_u(ii,kk,1) = buf_xz_u(ii,kk,1) + Vx(k,i,j) * UC * M0 * dt
        buf_xz_u(ii,kk,2) = buf_xz_u(ii,kk,2) + Vy(k,i,j) * UC * M0 * dt
        buf_xz_u(ii,kk,3) = buf_xz_u(ii,kk,3) + Vz(k,i,j) * UC * M0 * dt

      end do
    end do
    !$omp end parallel do

  end if

  if( mod( it-1, ntdec_s ) == 0 ) then

    if( snp_format == 'native' ) then
      call write_reduce_array2d_r( nxs, nzs, xz_u%ionode, xz_u%io, buf_xz_u(:,:,1) )
      call write_reduce_array2d_r( nxs, nzs, xz_u%ionode, xz_u%io, buf_xz_u(:,:,2) )
      call write_reduce_array2d_r( nxs, nzs, xz_u%ionode, xz_u%io, buf_xz_u(:,:,3) )
    else
      call write_reduce_array2d_r_nc( it, 1, nxs, nzs, xz_u, buf_xz_u(:,:,1) )
      call write_reduce_array2d_r_nc( it, 2, nxs, nzs, xz_u, buf_xz_u(:,:,2) )
      call write_reduce_array2d_r_nc( it, 3, nxs, nzs, xz_u, buf_xz_u(:,:,3) )
    end if

  end if

end subroutine wbuf_xz_u
!! --------------------------------------------------------------------------------------------------------------------------- !!

!! --------------------------------------------------------------------------------------------------------------------------- !!
subroutine wbuf_xy_u(it)
  integer, intent(in) :: it
  integer :: i, j, k
  integer :: ii, jj
  !! ---

  k = k0_xy

  !$omp parallel do private( ii, jj, i, j )
  do jj = js0, js1
    do ii = is0, is1
      j = jj * jdec - jdec/2
      i = ii * idec - idec/2

      buf_xy_u(ii,jj,1) = buf_xy_u(ii,jj,1) + Vx(k,i,j) * UC * M0 * dt
      buf_xy_u(ii,jj,2) = buf_xy_u(ii,jj,2) + Vy(k,i,j) * UC * M0 * dt
      buf_xy_u(ii,jj,3) = buf_xy_u(ii,jj,3) + Vz(k,i,j) * UC * M0 * dt

    end do
  end do
  !$omp end parallel do

  if( mod( it-1, ntdec_s ) == 0 ) then

    if( snp_format == 'native' ) then
      call write_reduce_array2d_r( nxs, nys, xy_u%ionode, xy_u%io, buf_xy_u(:,:,1) )
      call write_reduce_array2d_r( nxs, nys, xy_u%ionode, xy_u%io, buf_xy_u(:,:,2) )
      call write_reduce_array2d_r( nxs, nys, xy_u%ionode, xy_u%io, buf_xy_u(:,:,3) )
    else
      call write_reduce_array2d_r_nc( it, 1, nxs, nys, xy_u, buf_xy_u(:,:,1) )
      call write_reduce_array2d_r_nc( it, 2, nxs, nys, xy_u, buf_xy_u(:,:,2) )
      call write_reduce_array2d_r_nc( it, 3, nxs, nys, xy_u, buf_xy_u(:,:,3) )
    end if

  end if


end subroutine wbuf_xy_u
!! --------------------------------------------------------------------------------------------------------------------------- !!

!! --------------------------------------------------------------------------------------------------------------------------- !!
subroutine wbuf_fs_u(it)

  integer, intent(in) :: it
  integer :: i, j, k
  integer :: ii, jj


  !$omp parallel do private( ii, jj, i, j, k )
  do jj = js0, js1
    do ii = is0, is1
      j = jj * jdec - jdec/2
      i = ii * idec - idec/2
      k = kfs(i,j) + 1

      buf_fs_u(ii,jj,1) = buf_fs_u(ii,jj,1) + Vx(k,i,j) * UC * M0 * dt
      buf_fs_u(ii,jj,2) = buf_fs_u(ii,jj,2) + Vy(k,i,j) * UC * M0 * dt
      buf_fs_u(ii,jj,3) = buf_fs_u(ii,jj,3) + Vz(k,i,j) * UC * M0 * dt

    end do
  end do
  !$omp end parallel do

  if( mod( it-1, ntdec_s ) == 0 ) then
    if( snp_format == 'native' ) then
      call write_reduce_array2d_r( nxs, nys, fs_u%ionode, fs_u%io, buf_fs_u(:,:,1) )
      call write_reduce_array2d_r( nxs, nys, fs_u%ionode, fs_u%io, buf_fs_u(:,:,2) )
      call write_reduce_array2d_r( nxs, nys, fs_u%ionode, fs_u%io, buf_fs_u(:,:,3) )
    else
      call write_reduce_array2d_r_nc( it, 1, nxs, nys, fs_u, buf_fs_u(:,:,1) )
      call write_reduce_array2d_r_nc( it, 2, nxs, nys, fs_u, buf_fs_u(:,:,2) )
      call write_reduce_array2d_r_nc( it, 3, nxs, nys, fs_u, buf_fs_u(:,:,3) )
    end if
  end if

end subroutine wbuf_fs_u
!! --------------------------------------------------------------------------------------------------------------------------- !!

!! --------------------------------------------------------------------------------------------------------------------------- !!
subroutine wbuf_ob_u(it)

  integer, intent(in) :: it
  integer :: i, j, k
  integer :: ii, jj

  !$omp parallel do private( ii, jj, i, j, k )
  do jj = js0, js1
    do ii = is0, is1
      j = jj * jdec - jdec/2
      i = ii * idec - idec/2
      k = kob(i,j) + 1

      buf_ob_u(ii,jj,1) = buf_ob_u(ii,jj,1) + Vx(k,i,j) * UC * M0 * dt
      buf_ob_u(ii,jj,2) = buf_ob_u(ii,jj,2) + Vy(k,i,j) * UC * M0 * dt
      buf_ob_u(ii,jj,3) = buf_ob_u(ii,jj,3) + Vz(k,i,j) * UC * M0 * dt

    end do
  end do
  !$omp end parallel do

  if( mod( it-1, ntdec_s ) == 0 ) then
    if( snp_format == 'native' ) then
      call write_reduce_array2d_r( nxs, nys, ob_u%ionode, ob_u%io, buf_ob_u(:,:,1) )
      call write_reduce_array2d_r( nxs, nys, ob_u%ionode, ob_u%io, buf_ob_u(:,:,2) )
      call write_reduce_array2d_r( nxs, nys, ob_u%ionode, ob_u%io, buf_ob_u(:,:,3) )
    else
      call write_reduce_array2d_r_nc( it, 1, nxs, nys, ob_u, buf_ob_u(:,:,1) )
      call write_reduce_array2d_r_nc( it, 2, nxs, nys, ob_u, buf_ob_u(:,:,2) )
      call write_reduce_array2d_r_nc( it, 3, nxs, nys, ob_u, buf_ob_u(:,:,3) )
    end if
  end if

end subroutine wbuf_ob_u
!! --------------------------------------------------------------------------------------------------------------------------- !!



!! --------------------------------------------------------------------------------------------------------------------------- !!
!>
!! write waveform
!<
!! ----
subroutine output__store_wav(it)
  integer, intent(in) :: it
  integer :: i, itw

  if( .not. ( sw_wav_v .or. sw_wav_u ) ) return

  call pwatch__on( "output__write_wav" )


  if( it == 1 ) then
    allocate( ux(nst), uy(nst), uz(nst) )
    ux(:) = 0.0
    uy(:) = 0.0
    uz(:) = 0.0
  end if

  !! integrate waveform
  if( sw_wav_u ) then
    !$omp parallel do private(i)
    do i=1, nst
      ux(i) = ux(i) + ( Vx( kst(i), ist(i), jst(i) ) + Vx( kst(i),   ist(i)-1, jst(i)   ) ) * 0.5 * dt
      uy(i) = uy(i) + ( Vy( kst(i), ist(i), jst(i) ) + Vy( kst(i),   ist(i)  , jst(i)-1 ) ) * 0.5 * dt
      uz(i) = uz(i) - ( Vz( kst(i), ist(i), jst(i) ) + Vz( kst(i)-1, ist(i)  , jst(i)   ) ) * 0.5 * dt ! positive upward
    end do
    !$omp end parallel do
  end if


  !! output
  if( mod( it-1, ntdec_w ) == 0 ) then
    itw = (it-1)/ntdec_w + 1
    if( sw_wav_v ) then
      !$omp parallel do private(i)
      do i=1, nst
        vxst(itw,i) =  ( Vx( kst(i), ist(i), jst(i) ) + Vx( kst(i)  , ist(i)-1, jst(i)   ) )*0.5 * M0 * UC * 1e9 !! [nm/s]
        vyst(itw,i) =  ( Vy( kst(i), ist(i), jst(i) ) + Vy( kst(i)  , ist(i)  , jst(i)-1 ) )*0.5 * M0 * UC * 1e9 !! [nm/s]
        vzst(itw,i) = -( Vz( kst(i), ist(i), jst(i) ) + Vz( kst(i)-1, ist(i)  , jst(i)   ) )*0.5 * M0 * UC * 1e9 !! [nm/s]
      end do
      !$omp end parallel do
    end if

    if( sw_wav_u ) then
      !$omp parallel do private(i)
      do i=1, nst
        uxst(itw,i) = ux(i) * M0 * UC * 1e9                          !! [nm]
        uyst(itw,i) = uy(i) * M0 * UC * 1e9                          !! [nm]
        uzst(itw,i) = uz(i) * M0 * UC * 1e9                          !! [nm]
      end do
      !$omp end parallel do
    end if
  end if

  call pwatch__off( "output__write_wav" )

end subroutine output__store_wav
!! --------------------------------------------------------------------------------------------------------------------------- !!

!! --------------------------------------------------------------------------------------------------------------------------- !!
!>
!! Return location and indices of the requested station
!!
!! Used for Green's function computation with reciprocity theorem
!<
!! --
subroutine output__station_query( stnm1, is_exist, ist1, jst1, kst1, xst1, yst1, zst1, stlo1, stla1 )

  character(*), intent(in)  :: stnm1
  logical,      intent(out) :: is_exist
  integer,      intent(out) :: ist1, jst1, kst1 !< station location indices
  real(SP),     intent(out) :: xst1, yst1, zst1 !< station location coordinate
  real(SP),     intent(out) :: stlo1, stla1
  !! --
  integer :: i

  !! ----

  is_exist = .false.
  do i=1, nst
    if( trim(stnm1) == trim(stnm(i)) ) then

      !! copy station information
      ist1  = ist(i)
      jst1  = jst(i)
      kst1  = kst(i)
      xst1  = xst(i)
      yst1  = yst(i)
      zst1  = zst(i)
      stlo1 = stlo(i)
      stla1 = stla(i)
      is_exist = .true.
      exit
    end if
  end do


end subroutine output__station_query
!! --------------------------------------------------------------------------------------------------------------------------- !!

!! --------------------------------------------------------------------------------------------------------------------------- !!
subroutine output__checkpoint( io )
  integer, intent(in) :: io
  !! ----

  write( io ) xy_ps, yz_ps, xz_ps, fs_ps, ob_ps, xy_v, yz_v, xz_v, fs_v, ob_v, xy_u, yz_u, xz_u, fs_u, ob_u
  write( io ) sw_wav, sw_wav_u, sw_wav_v

  write( io ) ntdec_s
  write( io ) idec, jdec, kdec

  write( io ) z0_xy, x0_yz, y0_xz
  write( io ) k0_xy, i0_yz, j0_xz

  write( io ) nxs, nys, nzs
  write( io ) xsnp(1:nxs)
  write( io ) ysnp(1:nys)
  write( io ) zsnp(1:nzs)

  write( io ) is0, is1, js0, js1, ks0, ks1
  write( io ) r20x, r20y, r20z

  write( io ) buf_yz_u(js0:js1,ks0:ks1,1:3)
  write( io ) buf_xz_u(is0:is1,ks0:ks1,1:3)
  write( io ) buf_xy_u(is0:is1,js0:js1,1:3)
  write( io ) buf_fs_u(is0:is1,js0:js1,1:3)
  write( io ) buf_ob_u(is0:is1,js0:js1,1:3)

  write( io ) snp_format

  if( sw_wav ) then
    write( io ) ntdec_w
    write( io ) nst
    write( io ) ntw

    if( nst > 0 ) then
      write( io ) xst(1:nst)
      write( io ) yst(1:nst)
      write( io ) zst(1:nst)
      write( io ) ist(1:nst)
      write( io ) jst(1:nst)
      write( io ) kst(1:nst)
      write( io ) stlo(1:nst)
      write( io ) stla(1:nst)
      write( io ) stnm(1:nst)

      write( io ) vxst(1:ntw,1:nst)
      write( io ) vyst(1:ntw,1:nst)
      write( io ) vzst(1:ntw,1:nst)
      write( io ) uxst(1:ntw,1:nst)
      write( io ) uyst(1:ntw,1:nst)
      write( io ) uzst(1:ntw,1:nst)
      write( io ) sh(1:6,1:nst)
      write( io ) ux(1:nst)
      write( io ) uy(1:nst)
      write( io ) uz(1:nst)
    end if

  end if

end subroutine output__checkpoint
!! --------------------------------------------------------------------------------------------------------------------------- !!

!! --------------------------------------------------------------------------------------------------------------------------- !!
subroutine output__restart( io )
  integer, intent(in) :: io
  !! ----


  read( io ) xy_ps, yz_ps, xz_ps, fs_ps, ob_ps, xy_v, yz_v, xz_v, fs_v, ob_v, xy_u, yz_u, xz_u, fs_u, ob_u
  read( io ) sw_wav, sw_wav_u, sw_wav_v

  read( io ) ntdec_s
  read( io ) idec, jdec, kdec

  read( io ) z0_xy, x0_yz, y0_xz
  read( io ) k0_xy, i0_yz, j0_xz

  read( io ) nxs, nys, nzs

  allocate( xsnp(1:nxs) )
  allocate( ysnp(1:nys) )
  allocate( zsnp(1:nzs) )
  read( io ) xsnp(1:nxs)
  read( io ) ysnp(1:nys)
  read( io ) zsnp(1:nzs)

  read( io ) is0, is1, js0, js1, ks0, ks1
  read( io ) r20x, r20y, r20z

  allocate( buf_yz_u(nys,nzs,3), buf_xz_u(nxs,nzs,3), buf_xy_u(nxs,nys,3), buf_fs_u(nxs,nys,3), buf_ob_u(nxs,nys,3) )
  buf_yz_u = 0.0
  buf_xz_u = 0.0
  buf_xy_u = 0.0
  buf_fs_u = 0.0
  buf_ob_u = 0.0

  read( io ) buf_yz_u(js0:js1,ks0:ks1,1:3)
  read( io ) buf_xz_u(is0:is1,ks0:ks1,1:3)
  read( io ) buf_xy_u(is0:is1,js0:js1,1:3)
  read( io ) buf_fs_u(is0:is1,js0:js1,1:3)
  read( io ) buf_ob_u(is0:is1,js0:js1,1:3)

  read( io ) snp_format

  if( sw_wav ) then
    read( io ) ntdec_w
    read( io ) nst
    read( io ) ntw

    if( nst > 0 ) then
      allocate( xst(nst), yst(nst), zst(nst) )
      allocate( ist(nst), jst(nst), kst(nst) )
      allocate( stnm(nst) )
      allocate( stlo(nst), stla(nst) )

      read( io ) xst(1:nst)
      read( io ) yst(1:nst)
      read( io ) zst(1:nst)
      read( io ) ist(1:nst)
      read( io ) jst(1:nst)
      read( io ) kst(1:nst)
      read( io ) stlo(1:nst)
      read( io ) stla(1:nst)
      read( io ) stnm(1:nst)

      allocate( vxst(ntw,nst), vyst(ntw,nst), vzst(ntw,nst) )
      allocate( uxst(ntw,nst), uyst(ntw,nst), uzst(ntw,nst) )
      allocate( ux(nst), uy(nst), uz(nst) )
      allocate( sh(6,nst) )
      read( io ) vxst(1:ntw,1:nst)
      read( io ) vyst(1:ntw,1:nst)
      read( io ) vzst(1:ntw,1:nst)
      read( io ) uxst(1:ntw,1:nst)
      read( io ) uyst(1:ntw,1:nst)
      read( io ) uzst(1:ntw,1:nst)
      read( io ) sh(1:6,1:nst)
      read( io ) ux(1:nst)
      read( io ) uy(1:nst)
      read( io ) uz(1:nst)
    end if

  end if

  if( snp_format == 'native' ) then

#ifdef _ES
    if( yz_ps%sw .and. myid == yz_ps%ionode ) then
      open( yz_ps%io, file=trim(odir) //'/'// trim(title) //'.yz.ps.snp', &
      status='old', position='append', form='unformatted' )
    end if

    if( xz_ps%sw .and. myid == xz_ps%ionode ) then
      open( xz_ps%io, file=trim(odir) //'/'// trim(title) //'.xz.ps.snp', &
      status='old', position='append', form='unformatted' )
    end if

    if( xy_ps%sw .and. myid == xy_ps%ionode ) then
      open( xy_ps%io, file=trim(odir) //'/'// trim(title) //'.xy.ps.snp', &
      status='old', position='append', form='unformatted' )
    end if

    if( fs_ps%sw .and. myid == fs_ps%ionode ) then
      open( fs_ps%io, file=trim(odir) //'/'// trim(title) //'.fs.ps.snp', &
      status='old', position='append', form='unformatted' )
    end if

    if( ob_ps%sw .and. myid == ob_ps%ionode ) then
      open( ob_ps%io, file=trim(odir) //'/'// trim(title) //'.ob.ps.snp', &
      status='old', position='append', form='unformatted' )
    end if

    if( yz_v%sw .and. myid == yz_v%ionode ) then
      open( yz_v%io, file=trim(odir) //'/'// trim(title) //'.yz.v.snp', &
      status='old', position='append', form='unformatted' )
    end if

    if( xz_v%sw .and. myid == xz_v%ionode ) then
      open( xz_v%io, file=trim(odir) //'/'// trim(title) //'.xz.v.snp', &
      status='old', position='append', form='unformatted' )
    end if

    if( xy_v%sw .and. myid == xy_v%ionode ) then
      open( xy_v%io, file=trim(odir) //'/'// trim(title) //'.xy.v.snp', &
      status='old', position='append', form='unformatted' )
    end if

    if( fs_v%sw .and. myid == fs_v%ionode ) then
      open( fs_v%io, file=trim(odir) //'/'// trim(title) //'.fs.v.snp', &
      status='old', position='append', form='unformatted' )
    end if

    if( ob_v%sw .and. myid == ob_v%ionode ) then
      open( ob_v%io, file=trim(odir) //'/'// trim(title) //'.ob.v.snp', &
      status='old', position='append', form='unformatted' )
    end if
#else
    if( yz_ps%sw .and. myid == yz_ps%ionode ) then
      open( yz_ps%io, file=trim(odir) //'/'// trim(title) //'.yz.ps.snp', &
      status='old', access='stream', position='append', form='unformatted' )
    end if

    if( xz_ps%sw .and. myid == xz_ps%ionode ) then
      open( xz_ps%io, file=trim(odir) //'/'// trim(title) //'.xz.ps.snp', &
      status='old', access='stream', position='append', form='unformatted' )
    end if

    if( xy_ps%sw .and. myid == xy_ps%ionode ) then
      open( xy_ps%io, file=trim(odir) //'/'// trim(title) //'.xy.ps.snp', &
      status='old', access='stream', position='append', form='unformatted' )
    end if

    if( fs_ps%sw .and. myid == fs_ps%ionode ) then
      open( fs_ps%io, file=trim(odir) //'/'// trim(title) //'.fs.ps.snp', &
      status='old', access='stream', position='append', form='unformatted' )
    end if

    if( ob_ps%sw .and. myid == ob_ps%ionode ) then
      open( ob_ps%io, file=trim(odir) //'/'// trim(title) //'.ob.ps.snp', &
      status='old', access='stream', position='append', form='unformatted' )
    end if

    if( yz_v%sw .and. myid == yz_v%ionode ) then
      open( yz_v%io, file=trim(odir) //'/'// trim(title) //'.yz.v.snp', &
      status='old', access='stream', position='append', form='unformatted' )
    end if

    if( xz_v%sw .and. myid == xz_v%ionode ) then
      open( xz_v%io, file=trim(odir) //'/'// trim(title) //'.xz.v.snp', &
      status='old', access='stream', position='append', form='unformatted' )
    end if

    if( xy_v%sw .and. myid == xy_v%ionode ) then
      open( xy_v%io, file=trim(odir) //'/'// trim(title) //'.xy.v.snp', &
      status='old', access='stream', position='append', form='unformatted' )
    end if

    if( fs_v%sw .and. myid == fs_v%ionode ) then
      open( fs_v%io, file=trim(odir) //'/'// trim(title) //'.fs.v.snp', &
      status='old', access='stream', position='append', form='unformatted' )
    end if

    if( ob_v%sw .and. myid == ob_v%ionode ) then
      open( ob_v%io, file=trim(odir) //'/'// trim(title) //'.ob.v.snp', &
      status='old', access='stream', position='append', form='unformatted' )
    end if

    if( yz_u%sw .and. myid == yz_u%ionode ) then
      open( yz_u%io, file=trim(odir) //'/'// trim(title) //'.yz.u.snp', &
      status='old', access='stream', position='append', form='unformatted' )
    end if

    if( xz_u%sw .and. myid == xz_u%ionode ) then
      open( xz_u%io, file=trim(odir) //'/'// trim(title) //'.xz.u.snp', &
      status='old', access='stream', position='append', form='unformatted' )
    end if

    if( xy_u%sw .and. myid == xy_u%ionode ) then
      open( xy_u%io, file=trim(odir) //'/'// trim(title) //'.xy.u.snp', &
      status='old', access='stream', position='append', form='unformatted' )
    end if

    if( fs_u%sw .and. myid == fs_u%ionode ) then
      open( fs_u%io, file=trim(odir) //'/'// trim(title) //'.fs.u.snp', &
      status='old', access='stream', position='append', form='unformatted' )
    end if

    if( ob_u%sw .and. myid == ob_u%ionode ) then
      open( ob_u%io, file=trim(odir) //'/'// trim(title) //'.ob.u.snp', &
      status='old', access='stream', position='append', form='unformatted' )
    end if

#endif
  else

#ifdef _NETCDF
    if( yz_ps%sw .and. myid == yz_ps%ionode ) then
      call nc_chk( nf90_open( trim(odir) //'/'// trim(title) //'.yz.ps.nc', NF90_WRITE, yz_ps%io ) )
    end if

    if( xz_ps%sw .and. myid == xz_ps%ionode ) then
      call nc_chk( nf90_open( trim(odir) //'/'// trim(title) //'.xz.ps.nc', NF90_WRITE, xz_ps%io ) )
    end if

    if( xy_ps%sw .and. myid == xy_ps%ionode ) then
      call nc_chk( nf90_open( trim(odir) //'/'// trim(title) //'.xy.ps.nc', NF90_WRITE, xy_ps%io ) )
    end if

    if( fs_ps%sw .and. myid == fs_ps%ionode ) then
      call nc_chk( nf90_open( trim(odir) //'/'// trim(title) //'.fs.ps.nc', NF90_WRITE, fs_ps%io ) )
    end if

    if( ob_ps%sw .and. myid == ob_ps%ionode ) then
      call nc_chk( nf90_open( trim(odir) //'/'// trim(title) //'.ob.ps.nc', NF90_WRITE, ob_ps%io ) )
    end if


    if( yz_v%sw .and. myid == yz_v%ionode ) then
      call nc_chk( nf90_open( trim(odir) //'/'// trim(title) //'.yz.v.nc', NF90_WRITE, yz_v%io ) )
    end if

    if( xz_v%sw .and. myid == xz_v%ionode ) then
      call nc_chk( nf90_open( trim(odir) //'/'// trim(title) //'.xz.v.nc', NF90_WRITE, xz_v%io ) )
    end if

    if( xy_v%sw .and. myid == xy_v%ionode ) then
      call nc_chk( nf90_open( trim(odir) //'/'// trim(title) //'.xy.v.nc', NF90_WRITE, xy_v%io ) )
    end if

    if( fs_v%sw .and. myid == fs_v%ionode ) then
      call nc_chk( nf90_open( trim(odir) //'/'// trim(title) //'.fs.v.nc', NF90_WRITE, fs_v%io ) )
    end if

    if( ob_v%sw .and. myid == ob_v%ionode ) then
      call nc_chk( nf90_open( trim(odir) //'/'// trim(title) //'.ob.v.nc', NF90_WRITE, ob_v%io ) )
    end if


    if( yz_u%sw .and. myid == yz_u%ionode ) then
      call nc_chk( nf90_open( trim(odir) //'/'// trim(title) //'.yz.u.nc', NF90_WRITE, yz_u%io ) )
    end if

    if( xz_u%sw .and. myid == xz_u%ionode ) then
      call nc_chk( nf90_open( trim(odir) //'/'// trim(title) //'.xz.u.nc', NF90_WRITE, xz_u%io ) )
    end if

    if( xy_u%sw .and. myid == xy_u%ionode ) then
      call nc_chk( nf90_open( trim(odir) //'/'// trim(title) //'.xy.u.nc', NF90_WRITE, xy_u%io ) )
    end if

    if( fs_u%sw .and. myid == fs_u%ionode ) then
      call nc_chk( nf90_open( trim(odir) //'/'// trim(title) //'.fs.u.nc', NF90_WRITE, fs_u%io ) )
    end if

    if( ob_u%sw .and. myid == ob_u%ionode ) then
      call nc_chk( nf90_open( trim(odir) //'/'// trim(title) //'.ob.u.nc', NF90_WRITE, ob_u%io ) )
    end if
#endif
  end if

end subroutine output__restart
!! --------------------------------------------------------------------------------------------------------------------------- !!

subroutine close_nc( hdr )
  type(snp), intent(in) :: hdr
  integer :: vid

#ifdef _NETCDF
  call nc_chk( nf90_redef( hdr%io ) )
  do vid = 1, hdr%nsnp
    call nc_chk( nf90_put_att( hdr%io, hdr%varid(vid), 'actual_range', (/hdr%vmin(vid), hdr%vmax(vid)/)) )
  end do
  call nc_chk( nf90_enddef( hdr%io ) )
  call nc_chk( nf90_close( hdr%io ) )
#endif

end subroutine close_nc


!! --------------------------------------------------------------------------------------------------------------------------- !!
subroutine output__closefiles

  integer :: vid
  !! --

  if( snp_format == 'native' ) then
    if( yz_ps%sw .and. myid == yz_ps%ionode ) close( yz_ps%io )
    if( xz_ps%sw .and. myid == xz_ps%ionode ) close( xz_ps%io )
    if( xy_ps%sw .and. myid == xy_ps%ionode ) close( xy_ps%io )
    if( fs_ps%sw .and. myid == fs_ps%ionode ) close( fs_ps%io )
    if( ob_ps%sw .and. myid == ob_ps%ionode ) close( ob_ps%io )

    if( yz_v%sw .and. myid == yz_v%ionode ) close( yz_v%io )
    if( xz_v%sw .and. myid == xz_v%ionode ) close( xz_v%io )
    if( xy_v%sw .and. myid == xy_v%ionode ) close( xy_v%io )
    if( fs_v%sw .and. myid == fs_v%ionode ) close( fs_v%io )
    if( ob_v%sw .and. myid == ob_v%ionode ) close( ob_v%io )

    if( yz_u%sw .and. myid == yz_u%ionode ) close( yz_u%io )
    if( xz_u%sw .and. myid == xz_u%ionode ) close( xz_u%io )
    if( xy_u%sw .and. myid == xy_u%ionode ) close( xy_u%io )
    if( fs_u%sw .and. myid == fs_u%ionode ) close( fs_u%io )
    if( ob_u%sw .and. myid == ob_u%ionode ) close( ob_u%io )
  else

    if( yz_ps%sw .and. myid == yz_ps%ionode ) call close_nc( yz_ps )
    if( xz_ps%sw .and. myid == xz_ps%ionode ) call close_nc( xz_ps )
    if( xy_ps%sw .and. myid == xy_ps%ionode ) call close_nc( xy_ps )
    if( fs_ps%sw .and. myid == fs_ps%ionode ) call close_nc( fs_ps )
    if( ob_ps%sw .and. myid == ob_ps%ionode ) call close_nc( ob_ps )
    if( yz_v%sw .and. myid == yz_v%ionode )   call close_nc( yz_v )
    if( xz_v%sw .and. myid == xz_v%ionode )   call close_nc( xz_v )
    if( xy_v%sw .and. myid == xy_v%ionode )   call close_nc( xy_v )
    if( fs_v%sw .and. myid == fs_v%ionode )   call close_nc( fs_v )
    if( ob_v%sw .and. myid == ob_v%ionode )   call close_nc( ob_v )
    if( yz_u%sw .and. myid == yz_u%ionode )   call close_nc( yz_u )
    if( xz_u%sw .and. myid == xz_u%ionode )   call close_nc( xz_u )
    if( xy_u%sw .and. myid == xy_u%ionode )   call close_nc( xy_u )
    if( fs_u%sw .and. myid == fs_u%ionode )   call close_nc( fs_u )
    if( ob_u%sw .and. myid == ob_u%ionode )   call close_nc( ob_u )

  end if


end subroutine output__closefiles
!! --------------------------------------------------------------------------------------------------------------------------- !!


!! --------------------------------------------------------------------------------------------------------------------------- !!
!>
!! An internal subroutine to check error in netcdf function calls
!<
!! --
subroutine nc_chk( ierr )

  integer, intent(in) :: ierr
  !! ----

#ifdef _NETCDF
  if( ierr /= NF90_NOERR )  write(STDERR,*) NF90_STRERROR( ierr )
#endif

end subroutine nc_chk
!! --------------------------------------------------------------------------------------------------------------------------- !!

end module m_output
!! ----------------------------------------------------------------------------------------------------------------------------- !!
