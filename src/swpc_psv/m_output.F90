!! ----------------------------------------------------------------------------------------------------------------------------- !!
!>
!! Snapshot/waveform output
!!
!! @copyright
!!   Copyright 2013-2021 Takuto Maeda. All rights reserved. This project is released under the MIT license.
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
  use m_sac
  use m_readini
  use m_geomap
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

  !! -- Internal Parameters
  character(8), parameter :: BINARY_TYPE= "STREAMIO"
  character(8), parameter :: CODE_TYPE  = "SWPC_PSV"   !!< FIXED parameter for file header
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
    character(2) :: coordinate = 'xz'
    integer :: nmed = 3
    integer :: na1, na2 ! absorbing layer

    !! variables for netcdf mode
    integer :: did_x1, did_x2, did_t ! dimension id for independent vars
    integer :: vid_x1, vid_x2, vid_t ! variable  id for independent vars
    integer :: varid(10)           ! variable  id for dependent vars
    integer :: medid(10)            ! medium array id
    real    :: vmax(10), vmin(10)  ! max/min of dependent vars
    character(10) :: vname(10)
    character(10) :: vunit(10)
  end type snp


  type(snp) :: xz_ps, xz_v, xz_u
  logical   :: sw_wav, sw_wav_v, sw_wav_u, sw_wav_stress, sw_wav_strain

  !! switch
  integer   :: ntdec_s, ntdec_w                                !< time step decimation factor: Snap and Waves
  integer   :: idec,  kdec                                !< spatial decimation factor: x, y, z directions
  character(2) :: st_format                                       !< station file format
  character(3) :: wav_format
  character(256) :: fn_stloc

  integer :: nxs, nzs !< snapshot grid size
  real(SP), allocatable :: xsnp(:), zsnp(:)

  !! station
  integer :: nst
  real(SP), allocatable :: xst(:), zst(:)
  integer,  allocatable :: ist(:), kst(:)
  real(SP), allocatable :: stlo(:), stla(:)
  character(8), allocatable :: stnm(:)

  !! waveform
  integer :: ntw ! number of wave samples
  real(SP), allocatable :: wav_vel(:,:,:), wav_disp(:,:,:)
  real(SP), allocatable :: wav_stress(:,:,:), wav_strain(:,:,:)
  type(sac__hdr), allocatable :: sh_vel(:,:), sh_disp(:,:), sh_stress(:,:), sh_strain(:,:)
  real(SP), allocatable :: ux(:), uz(:)
  real(SP), allocatable :: exx(:), ezz(:), exz(:)

  !! I/O area in the node
  integer :: is0, is1, ks0, ks1

  !! derivative coefficient
  real(MP) :: r20x, r20z

  character(6) :: snp_format ! native or netcdf

  !! displacement snapshot buffer
  real(SP), allocatable :: buf_u(:,:,:)

  real(MP) :: r40x, r40z, r41x, r41z  
  logical :: wav_calc_dist

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
    integer :: i, k, ii, kk
    !! ----

    call pwatch__on( "output__setup" )

    !!
    !! read parameter
    !!
    call readini( io_prm, 'xz_ps%sw',  xz_ps%sw,  .false. )
    call readini( io_prm, 'xz_v%sw',   xz_v%sw,  .false. )
    call readini( io_prm, 'xz_u%sw',   xz_u%sw,  .false. )
    call readini( io_prm, 'idec',      idec,      1       )
    call readini( io_prm, 'kdec',      kdec,      1       )
    call readini( io_prm, 'ntdec_s',   ntdec_s,   10      )
    call readini( io_prm, 'ntdec_w',   ntdec_w,   10      )
    call readini( io_prm, 'st_format', st_format, 'xy'    )
    call readini( io_prm, 'fn_stloc',  fn_stloc,  ''      )
    call readini( io_prm, 'sw_wav_v',  sw_wav_v,  .false. )
    call readini( io_prm, 'sw_wav_u',  sw_wav_u,  .false. )
    call readini( io_prm, 'sw_wav_stress',  sw_wav_stress,  .false. )
    call readini( io_prm, 'sw_wav_strain',  sw_wav_strain,  .false. )

    call readini( io_prm, 'snp_format', snp_format, 'native' )
    call readini( io_prm, 'wav_format', wav_format, 'sac' )

    call readini( io_prm, 'wav_calc_dist', wav_calc_dist, .false. )

    sw_wav = ( sw_wav_v .or. sw_wav_u .or. sw_wav_stress .or. sw_wav_strain )

    !!
    !! snapshot
    !!

    !!
    !! snapshot size #2013-0440
    !!
    nxs = ( nx + (idec/2) ) / idec
    nzs = ( nz + (kdec/2) ) / kdec

    !!
    !! coordinate
    !!
    allocate( xsnp(nxs), zsnp(nzs) )
    do i=1, nxs
      ii = i*idec - (idec/2)
      xsnp(i) = i2x( ii, xbeg, real(dx) )
    end do
    do k=1, nzs
      kk = k*kdec - (kdec/2)
      zsnp(k) = k2z( kk, zbeg, real(dz) )
    end do

    !! snapshot region covered by the MPI node
    is0 = ceiling( ( ibeg + idec/2) / real(idec) )
    is1 = floor  ( ( iend + idec/2) / real(idec) )
    ks0 = ceiling( ( kbeg + kdec/2) / real(kdec) )
    ks1 = floor  ( ( kend + kdec/2) / real(kdec) )

    !!
    !! output node definition: cyclic
    !!
    xz_ps%ionode = mod( 1, nproc_x )
    xz_v%ionode = mod( 2, nproc_x )
    xz_u%ionode = mod( 3, nproc_x )

    !! number of snapshots per cycle
    xz_ps % nsnp = 2
    xz_v  % nsnp = 2
    xz_u  % nsnp = 2

    xz_ps % snaptype  = 'ps'
    xz_v  % snaptype  = 'v2'
    xz_u  % snaptype  = 'u2'

    xz_ps % vname(1) = 'divergence'
    xz_ps % vname(2) = 'rotation'
    xz_ps % vunit(1) = '1/s'
    xz_ps % vunit(2) = '1/s'
    xz_v  % vname(1) = 'Vx'
    xz_v  % vname(2) = 'Vz'
    xz_v  % vunit(1) = 'm/s'
    xz_v  % vunit(2) = 'm/s'
    xz_u  % vname(1) = 'Ux'
    xz_u  % vname(2) = 'Uz'
    xz_u  % vunit(1) = 'm'
    xz_u  % vunit(2) = 'm'

    !!
    !! open
    !!

    !!
    !! output settings
    !!
    if( snp_format == 'native' ) then
      if( xz_ps%sw ) call newfile_xz(   trim(odir) // '/' // trim(title) //'.xz.ps.snp',  xz_ps )
      if( xz_v%sw  ) call newfile_xz(   trim(odir) // '/' // trim(title) //'.xz.v.snp',   xz_v )
      if( xz_u%sw  ) call newfile_xz(   trim(odir) // '/' // trim(title) //'.xz.u.snp',   xz_u )
    else
      if( xz_ps%sw ) call newfile_xz_nc( trim(odir) // '/' // trim(title) //'.xz.ps.nc', xz_ps )
      if( xz_v%sw  ) call newfile_xz_nc( trim(odir) // '/' // trim(title) //'.xz.v.nc',  xz_v )
      if( xz_u%sw  ) call newfile_xz_nc( trim(odir) // '/' // trim(title) //'.xz.u.nc',  xz_u )
    end if

    !! always allocate displacement buffer
    allocate(buf_u(nxs,nzs,2))
    buf_u(:,:,:) = 0.0

    !! for taking derivatives
    !! FDM coefficients
    r40x = 9.0_MP /  8.0_MP / dx
    r40z = 9.0_MP /  8.0_MP / dz
    r41x = 1.0_MP / 24.0_MP / dx
    r41z = 1.0_MP / 24.0_MP / dz
    r20x = 1.  / dx
    r20z = 1.  / dz

!!!!
!!!! waveform
!!!!

    if( sw_wav ) then

      ntw = floor( real(nt-1)/real(ntdec_w) + 1.0 )

      call read_stinfo()

    end if


    call pwatch__off( "output__setup" )

  end subroutine output__setup
  !! --------------------------------------------------------------------------------------------------------------------------- !!


  !! --------------------------------------------------------------------------------------------------------------------------- !!
  !>
  !<
  !! --
  subroutine output__export_wav()
    integer :: i, j
    character(6) :: cid
    integer :: io
    character(256) :: fn
    
    call pwatch__on("output__export_wav")

    if( nst>0 ) call system__call('mkdir '//trim(odir)//'/wav > /dev/null 2>&1' )

    if( wav_format == 'sac' ) then
      do i=1, nst

        if( sw_wav_v ) then
          do j=1, 2
            call export_wav__sac(sh_vel(j,i), wav_vel(:,j,i))
          end do
        end if

        if( sw_wav_u ) then
          do j=1, 2
            call export_wav__sac(sh_disp(j,i), wav_disp(:,j,i))
          end do
        end if

        if( sw_wav_stress ) then
          do j=1, 3
            call export_wav__sac( sh_stress(j,i), wav_stress(:,j,i) )
          end do
        end if

        if( sw_wav_strain ) then
          do j=1, 3
            call export_wav__sac( sh_strain(j,i), wav_strain(:,j,i) )
          end do
        end if
        
      end do

    else if ( wav_format == 'csf' ) then

      if( sw_wav_v      ) call export_wav__csf(nst, 3, sh_vel, wav_vel)
      if( sw_wav_u      ) call export_wav__csf(nst, 3, sh_disp, wav_disp)
      if( sw_wav_stress ) call export_wav__csf(nst, 3, sh_stress, wav_stress )
      if( sw_wav_strain ) call export_wav__csf(nst, 3, sh_strain, wav_strain )

    else if ( wav_format == 'wav' ) then

      write(cid,'(I6.6)') myid
      fn = trim(odir) // '/wav/' // trim(title) // '.' // trim(cid) // '.wav'

#ifdef _ES
      call std__getio(io, is_big=.true.)
      open(io, file=trim(fn), form='unformatted', action='write', status='replace')
#else
      call std__getio(io) 
      open(io, file=trim(fn), access='stream', form='unformatted', action='write', status='replace')
#endif

      if( sw_wav_v ) write(io) nst, ntw, title, sh_vel(:,:), wav_vel(:,:,:)
      if( sw_wav_u ) write(io) nst, ntw, title, sh_disp(:,:), wav_disp(:,:,:)
      if( sw_wav_stress ) write(io) nst, ntw, title, sh_stress(:,:), wav_stress(:,:,:)
      if( sw_wav_strain ) write(io) nst, ntw, title, sh_strain(:,:), wav_strain(:,:,:)

      close(io)

    end if

    call pwatch__off("output__export_wav")

  contains
    
    subroutine export_wav__sac( sh, dat )

      type(sac__hdr), intent(in) :: sh
      real(SP), intent(in) :: dat(:)
      character(256) :: fn
      !! --
      
      fn = trim(odir) // '/wav/' // trim(title) // '.' // trim(sh%kstnm) // '.' // trim(sh%kcmpnm) // '.sac'
      call sac__write( fn, sh, dat, .true. )
      
    end subroutine export_wav__sac

    subroutine export_wav__csf(nst, ncmp, sh, dat)
      
      integer, intent(in) :: nst, ncmp
      type(sac__hdr), intent(in) :: sh(ncmp, nst)
      real(SP), intent(in) :: dat(ntw, ncmp, nst)
      character(5) :: cid
      character(256) :: fn

      write(cid,'(I5.5)') myid
      fn = trim(odir) // '/wav/' // trim(title) // '__' // cid // '__.csf'
      call csf__write(fn, nst*ncmp, ntw, reshape(sh,(/ncmp*nst/)), reshape(dat, (/ntw, ncmp*nst/)))

    end subroutine export_wav__csf
    
  end subroutine output__export_wav
  !! --------------------------------------------------------------------------------------------------------------------------- !!


  !! --------------------------------------------------------------------------------------------------------------------------- !!
  !>
  !! Set-up waveform file. Create, header output.
  !<
  !! --
  subroutine read_stinfo( )

    integer                   :: io_stlst
    real(SP), allocatable     :: xst_g(:), zst_g(:), stlo_g(:), stla_g(:)
    integer, allocatable      :: ist_g(:), kst_g(:), stid(:)
    character(3), allocatable :: zsw_g(:)
    character(8), allocatable :: stnm_g(:)
    character(256)            :: abuf
    integer                   :: ierr
    integer :: nst_g
    integer :: i, ii, j
    real(SP) :: rdum
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
      sw_wav_v = .false.
      sw_wav_u = .false.
      sw_wav_stress = .false.
      sw_wav_strain = .false.
      
      if( myid == 0 ) then
        write(STDERR,'(A)') "[INFO] output--read_stinfo: Station file does not exist. WAV file will not be created"
      end if
      return

    else
      call std__countline( io_stlst, nst_g, "#" )
    end if

    !! allocate temp memory
    allocate( xst_g(nst_g), zst_g(nst_g), stnm_g(nst_g), zsw_g(nst_g), stlo_g(nst_g), stla_g(nst_g) )
    allocate( ist_g(nst_g), kst_g(nst_g), stid(nst_g) )


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
      if( adjustl(abuf(1:1)) == "#" ) cycle ! neglect comment line
      if( trim(adjustl(abuf)) == "" ) cycle ! neglect blank line

      i = i + 1

      select case ( st_format )

      case( 'xy' )

        read(abuf,*) xst_g(i), rdum, zst_g(i), stnm_g(i), zsw_g(i)
        call geomap__c2g( xst_g(i), 0.0, clon, clat, phi, stlo_g(i), stla_g(i) )


      case( 'll' )

        read(abuf,*) stlo_g(i), stla_g(i), zst_g(i), stnm_g(i), zsw_g(i)
        call geomap__g2c( stlo_g(i), stla_g(i), clon, clat, phi, xst_g(i), rdum )

      case default

        write(STDERR,'(A)') "ERROR [output__stinfo]: unrecognized station format."
        stop

      end select


      ! digitize
      ist_g(i) = x2i ( xst_g(i), xbeg, real(dx) )
      kst_g(i) = z2k ( zst_g(i), zbeg, real(dz) )


      !! check if the station is in the computational domain
      if(  i2x( 1, xbeg, real(dx) ) < xst_g(i) .and. xst_g(i) < i2x( nx, xbeg, real(dx) ) .and. &
          kbeg             < kst_g(i) .and. kst_g(i) <          kend           )  then


        !! memorize station location
        nst_g = nst_g + 1

        !! MPI region check: memorize station number
        if( ibeg <= ist_g(i) .and. ist_g(i) <= iend ) then

          nst = nst + 1
          stid(nst) = i

        end if

      else
        if( myid == 0 )  write(STDERR,*) "WARNING [output__setup]: station " // trim( stnm_g(i) ) // " is out of the region"
      end if
    end do

    close( io_stlst )

    if( nst_g == 0 ) then
      sw_wav = .false.
      if( myid == 0 ) then
        write(STDERR,*) "[INFO] output--read_stinfo: No station is detected. WAV file will not be created"
      end if
      return
    end if


    if( nst == 0 ) return


    allocate( xst(nst), zst(nst), ist(nst), kst(nst), stnm(nst), stlo(nst), stla(nst) )


    do i = 1, nst

      ii = stid(i)

      xst(i) = xst_g(ii)
      zst(i) = zst_g(ii)
      ist(i) = ist_g(ii)
      ist(i) = ist_g(ii)
      stnm(i) = stnm_g(ii)
      stlo(i) = stlo_g(ii)
      stla(i) = stla_g(ii)

      !! station depth setting
      !! this is done only for MPI-node because of the definition of kfs and kob
      select case ( zsw_g(ii) )
      case( 'dep' );    kst(i) = kst_g(ii)
      case( 'fsb' );    kst(i) = kfs( ist(i) ) + 1   !! free surface: one-grid below for avoiding vacuum
      case( 'obb' );    kst(i) = kob( ist(i) ) + 1   !! ocean column: below seafloor
      case( 'oba' );    kst(i) = kob( ist(i) ) - 1   !! ocean column: above seafloor
      case( 'bd0' );    kst(i) = z2k( bddep(ist(i),0), zbeg, real(dz) ) !! Boundary interface
      case( 'bd1' );    kst(i) = z2k( bddep(ist(i),1), zbeg, real(dz) ) !! Boundary interface
      case( 'bd2' );    kst(i) = z2k( bddep(ist(i),2), zbeg, real(dz) ) !! Boundary interface
      case( 'bd3' );    kst(i) = z2k( bddep(ist(i),3), zbeg, real(dz) ) !! Boundary interface
      case( 'bd4' );    kst(i) = z2k( bddep(ist(i),4), zbeg, real(dz) ) !! Boundary interface
      case( 'bd5' );    kst(i) = z2k( bddep(ist(i),5), zbeg, real(dz) ) !! Boundary interface
      case( 'bd6' );    kst(i) = z2k( bddep(ist(i),6), zbeg, real(dz) ) !! Boundary interface
      case( 'bd7' );    kst(i) = z2k( bddep(ist(i),7), zbeg, real(dz) ) !! Boundary interface
      case( 'bd8' );    kst(i) = z2k( bddep(ist(i),8), zbeg, real(dz) ) !! Boundary interface
      case( 'bd9' );    kst(i) = z2k( bddep(ist(i),9), zbeg, real(dz) ) !! Boundary interface
      case default; kst(i) = kst_g(ii)
      end select

      !! depth check
      if( kst(i) > kend ) then
        write(STDERR,*) 'WARNING[output__setup]: station depth exceeds kend at station ' // trim(stnm(i))
        kst(i) = kend - 1
      end if
      if( kst(i) < kbeg ) then
        write(STDERR,*) 'WARNING[output__setup]: station depth exceeds kbeg at station ' // trim(stnm(i))
        kst(i) = kbeg + 1
      end if

    end do

    if( sw_wav_v ) then
      allocate(wav_vel(ntw,2,nst))
      allocate(sh_vel(2,nst))
      wav_vel(:,:,:) = 0.0
    end if
    
    if( sw_wav_u ) then
      allocate(wav_disp(ntw,2,nst))
      allocate(sh_disp(2,nst))
      wav_vel(:,:,:) = 0.0
    end if
    
      
    if( sw_wav_stress ) then
      allocate( wav_stress(ntw,3,nst) )
      allocate( sh_stress(3,nst) )
      wav_stress(:,:,:) = 0.0
    end if

    if( sw_wav_strain ) then
      allocate( wav_strain(ntw,3,nst) )
      allocate( sh_strain(3,nst) )
      wav_strain(:,:,:) = 0.0
    end if

    !!
    !! set-up sac header
    !!
    do i=1, nst
      
      if( sw_wav_v ) then
        do j=1, 2
          call setup_sac_header( sh_vel(j,i), i )
        end do
        sh_vel(1,i)%kcmpnm = "Vx"
        sh_vel(2,i)%kcmpnm = "Vz"
        sh_vel(:,i)%idep = 7 ! velocity [nm/s]
        sh_vel(1,i)%cmpinc = 90.0;  sh_vel(1,i)%cmpaz  =  0.0 + phi
        sh_vel(2,i)%cmpinc =  0.0;  sh_vel(2,i)%cmpaz  =  0.0

        if( wav_calc_dist ) then
          sh_vel(:,i)%lcalda = .false. 
          sh_vel(:,i)%dist = sqrt( (sx0 - xst(i))**2  )
          sh_vel(:,i)%az = std__rad2deg(atan2(0., xst(i)-sx0))
          sh_vel(:,i)%baz = std__rad2deg(atan2(0., sx0-xst(i)))
        end if
      end if

      if( sw_wav_u ) then
        do j=1, 2
          call setup_sac_header( sh_disp(j,i), i )
          sh_disp(1,i)%kcmpnm = "Ux"
          sh_disp(2,i)%kcmpnm = "Uz"
        end do
        sh_disp(:,i)%idep = 6 ! displacement [nm]
        sh_disp(1,i)%cmpinc = 90.0;  sh_disp(1,i)%cmpaz  =  0.0 + phi
        sh_disp(2,i)%cmpinc =  0.0;  sh_disp(2,i)%cmpaz  =  0.0

        if( wav_calc_dist ) then
          sh_disp(:,i)%lcalda = .false. 
          sh_disp(:,i)%dist = sqrt( (sx0 - xst(i))**2  )
          sh_disp(:,i)%az = std__rad2deg(atan2(0., xst(i)-sx0))
          sh_disp(:,i)%baz = std__rad2deg(atan2(0., sx0-xst(i)))
        end if

      end if

      if( sw_wav_stress ) then
        do j=1, 3
          call setup_sac_header( sh_stress(j,i), i )
        end do
        sh_stress(1,i)%kcmpnm = "Sxx"
        sh_stress(2,i)%kcmpnm = "Szz"
        sh_stress(3,i)%kcmpnm = "Sxz"
        
        sh_stress(:,i)%idep = 5 ! unknown

         if( wav_calc_dist ) then
          sh_stress(:,i)%lcalda = .false. 
          sh_stress(:,i)%dist = sqrt( (sx0 - xst(i))**2  )
          sh_stress(:,i)%az = std__rad2deg(atan2(0., xst(i)-sx0))
          sh_stress(:,i)%baz = std__rad2deg(atan2(0., sx0-xst(i)))
        end if

      end if
      
      if( sw_wav_strain ) then
        do j=1, 3
          call setup_sac_header( sh_strain(j,i), i )
        end do
        sh_strain(1,i)%kcmpnm = "Exx"
        sh_strain(2,i)%kcmpnm = "Ezz"
        sh_strain(3,i)%kcmpnm = "Exz"
        
        sh_strain(:,i)%idep = 5 ! unknown
        
        if( wav_calc_dist ) then
          sh_strain(:,i)%lcalda = .false. 
          sh_strain(:,i)%dist = sqrt( (sx0 - xst(i))**2  )
          sh_strain(:,i)%az = std__rad2deg(atan2(0., xst(i)-sx0))
          sh_strain(:,i)%baz = std__rad2deg(atan2(0., sx0-xst(i)))
        end if
      end if

    end do

  contains
    
    subroutine setup_sac_header( sh, ist )

      type(sac__hdr), intent(out) :: sh
      integer,        intent(in)  :: ist
      !! --

      call sac__init(sh)

      !! common header
      sh%evlo    = evlo
      sh%evla    = evla
      sh%evdp    = evdp  !! evdp changed to km unit from SWPC 5.0
      sh%tim     = exedate
      sh%b       = tbeg
      sh%delta   = ntdec_w * dt
      sh%npts    = ntw
      sh%mag     = mw
      
      if( bf_mode ) then
        sh%user0   = fx0
        sh%user1   = fy0
        sh%user2   = fz0
      else
        sh%user0   = mxx0
        sh%user1   = myy0
        sh%user2   = mzz0
        sh%user3   = myz0
        sh%user4   = mxz0
        sh%user5   = mxy0
      end if

      sh%user6   = clon  !< coordinate
      sh%user7   = clat  !< coordinate
      sh%user8   = phi
      sh%o       = otim

      call daytim__localtime( sh%tim, sh%nzyear, sh%nzmonth, sh%nzday, sh%nzhour, sh%nzmin, sh%nzsec )
      call daytim__ymd2jul  ( sh%nzyear, sh%nzmonth, sh%nzday, sh%nzjday )
      sh%nzmsec = 0

      !! station dependent
      sh%kevnm = trim(adjustl( title(1:16) ))
      sh%kstnm = trim(stnm(ist))
      sh%stlo  = stlo(ist)
      sh%stla  = stla(ist)
      sh%stdp  = zst(ist)*1000 ! in meter unit
      
    end subroutine setup_sac_header
    
  end subroutine read_stinfo
  !! --------------------------------------------------------------------------------------------------------------------------- !!


  !! --------------------------------------------------------------------------------------------------------------------------- !!
  !>
  !! Open new file and write header information, and medium parameters for XZ-cross section
  !<
  !! ----
  subroutine newfile_xz ( fname, hdr )

    character(*), intent(in)    :: fname
    type(snp),    intent(inout) :: hdr
    integer :: io0
    !! ----

    integer :: i, k, kk, ii
    real(SP) :: buf(nxs,nzs,3)


    hdr % na1 = na / idec
    hdr % na2 = na / kdec
    hdr % coordinate = 'xz'

    if( myid == hdr%ionode ) then

      call std__getio( hdr%io, is_big=.true. )
#ifdef _ES
      open( hdr%io, file=trim(fname), action='write', status='replace', form='unformatted' )
#else
      open( hdr%io, file=trim(fname), access='stream', action='write', status='replace', form='unformatted' )
#endif
      call write_snp_header( hdr, nxs, nzs, xsnp(1:nxs), zsnp(1:nzs) )

    end if

    buf = 0.0
    do i = is0, is1
      do k = ks0, ks1

        ii = i * idec - idec/2
        kk = k * kdec - kdec/2

        buf(i,k,1) = rho( kk, ii )
        buf(i,k,2) = lam( kk, ii )
        buf(i,k,3) = mu ( kk, ii )

      end do
    end do

    call write_reduce_array2d_r( nxs, nzs, hdr%ionode, hdr%io, buf(:,:,1) )
    call write_reduce_array2d_r( nxs, nzs, hdr%ionode, hdr%io, buf(:,:,2) )
    call write_reduce_array2d_r( nxs, nzs, hdr%ionode, hdr%io, buf(:,:,3) )

  end subroutine newfile_xz
  !! --------------------------------------------------------------------------------------------------------------------------- !!

  !! --------------------------------------------------------------------------------------------------------------------------- !!
  subroutine newfile_xz_nc ( fname, hdr )
    character(*), intent(in)    :: fname
    type(snp),    intent(inout) :: hdr
    !!
    real, allocatable :: sbuf(:), rbuf1(:), rbuf2(:), rbuf3(:), buf(:,:,:)
    integer :: i, k, ierr, ii, kk
    !! ---

#ifdef _NETCDF

    if( myid == hdr%ionode ) then

      !! initialize
      hdr%vmax = 0.0
      hdr%vmin = 0.0
      hdr%na1  = na / idec
      hdr%na2  = na / kdec

      call nc_chk( nf90_create( trim(fname), NF90_CLOBBER, hdr%io ) )
      call write_nc_header( hdr, nxs, nzs, xsnp, zsnp )
    end if

    allocate( buf(nxs,nzs,3) )
    buf = 0.0
    do i = is0, is1
      do k = ks0, ks1

        ii = i * idec - idec/2
        kk = k * kdec - kdec/2

        buf(i,k,1) = rho( kk, ii )
        buf(i,k,2) = lam( kk, ii )
        buf(i,k,3) = mu ( kk, ii )

      end do
    end do

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
      call nc_chk( nf90_put_var( hdr%io, hdr%medid(1), reshape(rbuf1,shape(buf(:,:,1)) ) ) )
      call nc_chk( nf90_put_var( hdr%io, hdr%medid(2), reshape(rbuf2,shape(buf(:,:,1)) ) ) )
      call nc_chk( nf90_put_var( hdr%io, hdr%medid(3), reshape(rbuf3,shape(buf(:,:,1)) ) ) )

    end if

    deallocate( sbuf, rbuf1, rbuf2, rbuf3, buf )

#endif
  end subroutine newfile_xz_nc
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
    real(SP) :: dum = 0
    !! ----

    write( hdr%io ) BINARY_TYPE
    write( hdr%io ) CODE_TYPE
    write( hdr%io ) HEADER_VERSION
    write( hdr%io ) title
    write( hdr%io ) exedate

    !! space grid size
    write( hdr%io ) hdr%coordinate                ! coordinate
    write( hdr%io ) hdr%snaptype                  ! data type
    write( hdr%io ) ns1
    write( hdr%io ) ns2                       ! data size
    write( hdr%io ) xs1(1)
    write( hdr%io ) xs2(1)
    write( hdr%io ) xs1(2) - xs1(1)
    write( hdr%io ) xs2(2) - xs2(1)

    !! time
    write( hdr%io ) dt * ntdec_s               ! dt

    ! absorb layer
    write( hdr%io ) hdr%na1
    write( hdr%io ) hdr%na2

    !! number of arrays
    write( hdr%io ) hdr%nmed
    write( hdr%io ) hdr%nsnp

    !! coordinate
    write( hdr%io ) clon
    write( hdr%io ) clat
    write( hdr%io ) phi

    write( hdr%io ) dum
    write( hdr%io ) dum
    write( hdr%io ) dum

  end subroutine write_snp_header
  !! --------------------------------------------------------------------------------------------------------------------------- !!
  !! --------------------------------------------------------------------------------------------------------------------------- !!
  !>
  !! write netcdf header
  !<
  !! ---
  subroutine write_nc_header( hdr, ns1, ns2, xs1, xs2 )

    type(snp), intent(inout) :: hdr
    integer, intent(in) :: ns1, ns2
    real(SP), intent(in) :: xs1(ns1), xs2(ns2)
    integer :: i
    !! --

#ifdef _NETCDF
    call nc_chk( nf90_def_dim( hdr%io, 'x', nxs, hdr%did_x1 ) )
    call nc_chk( nf90_def_dim( hdr%io, 'z', nzs, hdr%did_x2 ) )
    call nc_chk( nf90_def_dim( hdr%io, 't', NF90_UNLIMITED, hdr%did_t ) )
    call nc_chk( nf90_def_var( hdr%io, 'x', NF90_REAL, hdr%did_x1, hdr%vid_x1 ) )
    call nc_chk( nf90_def_var( hdr%io, 'z', NF90_REAL, hdr%did_x2, hdr%vid_x2 ) )
    call nc_chk( nf90_def_var( hdr%io, 't', NF90_REAL, hdr%did_t, hdr%vid_t ) )

    !! medium
    call nc_chk( nf90_def_var( hdr%io, 'rho',    NF90_REAL, (/hdr%did_x1, hdr%did_x2/), hdr%medid(1) ) )
    call nc_chk( nf90_def_var( hdr%io, 'lambda', NF90_REAL, (/hdr%did_x1, hdr%did_x2/), hdr%medid(2) ) )
    call nc_chk( nf90_def_var( hdr%io, 'mu',     NF90_REAL, (/hdr%did_x1, hdr%did_x2/), hdr%medid(3) ) )

    !! variables
    do i=1, hdr%nsnp
      call nc_chk( nf90_def_var( hdr%io, trim(hdr%vname(i)), NF90_REAL, (/hdr%did_x1, hdr%did_x2, hdr%did_t/), hdr%varid(i) ) )
      call nc_chk( nf90_put_att( hdr%io, hdr%varid(i), 'long_name', trim(hdr%vname(i)) ) )
      call nc_chk( nf90_put_att( hdr%io, hdr%varid(i), 'units', trim(hdr%vunit(i)) ) )
    end do

    !! global attribute
    call nc_chk( nf90_put_att( hdr%io, NF90_GLOBAL, 'generated_by',   'SWPC' ) )
    call nc_chk( nf90_put_att( hdr%io, NF90_GLOBAL, 'title',          trim(title) ) )
    call nc_chk( nf90_put_att( hdr%io, NF90_GLOBAL, 'exedate',        exedate     ) )
    call nc_chk( nf90_put_att( hdr%io, NF90_GLOBAL, 'hdrver',         HEADER_VERSION ) )
    call nc_chk( nf90_put_att( hdr%io, NF90_GLOBAL, 'codetype',       CODE_TYPE ) )
    call nc_chk( nf90_put_att( hdr%io, NF90_GLOBAL, 'ns1',            ns1 ) )
    call nc_chk( nf90_put_att( hdr%io, NF90_GLOBAL, 'ns2',            ns2 ) )
    call nc_chk( nf90_put_att( hdr%io, NF90_GLOBAL, 'beg1',           xs1(1) ) )
    call nc_chk( nf90_put_att( hdr%io, NF90_GLOBAL, 'beg2',           xs2(1) ) )
    call nc_chk( nf90_put_att( hdr%io, NF90_GLOBAL, 'na1',            hdr%na1 ) )
    call nc_chk( nf90_put_att( hdr%io, NF90_GLOBAL, 'na2',            hdr%na2 ) )
    call nc_chk( nf90_put_att( hdr%io, NF90_GLOBAL, 'ds1',            idec * dx ) )
    call nc_chk( nf90_put_att( hdr%io, NF90_GLOBAL, 'ds2',            kdec * dz ) )
    call nc_chk( nf90_put_att( hdr%io, NF90_GLOBAL, 'nmed',           hdr%nmed ) )
    call nc_chk( nf90_put_att( hdr%io, NF90_GLOBAL, 'nsnp',           hdr%nsnp ) )
    call nc_chk( nf90_put_att( hdr%io, NF90_GLOBAL, 'coordinate',     hdr%coordinate ) )
    call nc_chk( nf90_put_att( hdr%io, NF90_GLOBAL, 'datatype',       hdr%snaptype ) )
    call nc_chk( nf90_put_att( hdr%io, NF90_GLOBAL, 'dt',             dt * ntdec_s ) )
    
    call nc_chk( nf90_put_att( hdr%io, NF90_GLOBAL, 'evlo',           evlo))
    call nc_chk( nf90_put_att( hdr%io, NF90_GLOBAL, 'evla',           evla))
    call nc_chk( nf90_put_att( hdr%io, NF90_GLOBAL, 'evdp',           evdp))
    call nc_chk( nf90_put_att( hdr%io, NF90_GLOBAL, 'evx',            sx0))
    call nc_chk( nf90_put_att( hdr%io, NF90_GLOBAL, 'evy',            sy0))


    !! variable attributes
    call nc_chk( nf90_put_att( hdr%io, hdr%vid_x1, 'long_name', 'x' ) )
    call nc_chk( nf90_put_att( hdr%io, hdr%vid_x2, 'long_name', 'z' ) )
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
      count = (/ nxs, nzs, 1/)
      start = (/ 1, 1, it/ntdec_s+1 /)
      call nc_chk( nf90_put_var( hdr%io, hdr%varid(vid), reshape(rbuf,shape(array)), start, count ))
      call nc_chk( nf90_put_var( hdr%io, hdr%vid_t, it*dt, start=(/ it/ntdec_s+1 /) ) )
      hdr%vmax(vid) = max( hdr%vmax(vid), maxval(rbuf) )
      hdr%vmin(vid) = min( hdr%vmin(vid), minval(rbuf) )
    end if
#endif

  end subroutine write_reduce_array2d_r_nc
  !! --------------------------------------------------------------------------------------------------------------------------- !!

  !! --------------------------------------------------------------------------------------------------------------------------- !!
  !>
  !! Calculate divergence and abs(rotation) of velocity vector at given location (k,i,j) by 2nd order FDM
  !<
  !! ----
  subroutine divrot( k, i, div, rot )
    !! --
    integer,  intent(in)  :: k, i
    real(SP), intent(out) :: div, rot
    !!

    real(SP) :: rot_y
    real(SP) :: dxVx, dxVz, dzVx, dzVz
    real(SP) :: nnn, pnn, npn, ppn, mu_xz


    dxVx = (  Vx(k  ,i  ) - Vx(k  ,i-1)  ) * r20x
    dxVz = (  Vz(k  ,i+1) - Vz(k  ,i  )  ) * r20x
    dzVx = (  Vx(k+1,i  ) - Vx(k  ,i  )  ) * r20z
    dzVz = (  Vz(k  ,i  ) - Vz(k-1,i  )  ) * r20z

    div   = ( dxVx + dzVz )
    rot_y = ( dzVx - dxVz )

    !! masking
    nnn = mu (k  ,i  )
    pnn = mu (k+1,i  )
    npn = mu (k,  i+1)
    ppn = mu (k+1,i+1)
    mu_xz = 4*nnn*pnn*npn*ppn / ( nnn*pnn*npn + nnn*pnn*ppn + nnn*npn*ppn + pnn*npn*ppn + epsilon(1.0) )

    div   = div * lam(k,i) / ( abs( lam(k,i) ) + epsilon(1.0) )
    rot_y = rot_y * mu_xz / abs( mu_xz + epsilon(1.0) )
    rot   = rot_y

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


    if( xz_ps%sw ) call wbuf_xz_ps(it)
    if( xz_v%sw ) call wbuf_xz_v(it)
    if( xz_u%sw ) call wbuf_xz_u(it)

    call pwatch__off( "output__write_snap" )

  end subroutine output__write_snap
  !! --------------------------------------------------------------------------------------------------------------------------- !!




  !! --------------------------------------------------------------------------------------------------------------------------- !!
  subroutine wbuf_xz_ps(it)

    integer, intent(in) :: it
    integer :: i, k
    integer :: ii, kk
    real(SP) :: div, rot
    real(SP) :: buf(nxs,nzs,2)


    if( mod( it-1, ntdec_s ) /= 0 ) return

    buf(:,:,:) = 0.0

    !$omp parallel do private( ii, kk, i, k, div, rot )
    do ii = is0, is1
      do kk= ks0, ks1
        k = kk * kdec - kdec/2
        i = ii * idec - idec/2

        call divrot( k, i, div, rot )
        buf(ii,kk,1) = div
        buf(ii,kk,2) = rot

      end do
    end do
    !$omp end parallel do

    !! dx, dz have km unit. correction for 1e3 factor.
    buf =  buf * UC * M0 * 1e-3

    if( snp_format == 'native' ) then
      call write_reduce_array2d_r( nxs, nzs, xz_ps%ionode, xz_ps%io, buf(:,:,1) )
      call write_reduce_array2d_r( nxs, nzs, xz_ps%ionode, xz_ps%io, buf(:,:,2) )
    else
      call write_reduce_array2d_r_nc( it, 1, nxs, nzs, xz_ps, buf(:,:,1) )
      call write_reduce_array2d_r_nc( it, 2, nxs, nzs, xz_ps, buf(:,:,2) )
    end if

  end subroutine wbuf_xz_ps
  !! --------------------------------------------------------------------------------------------------------------------------- !!

  !! --------------------------------------------------------------------------------------------------------------------------- !!
  subroutine wbuf_xz_v(it)

    integer, intent(in) :: it
    integer :: i, k
    integer :: ii, kk
    real(SP) :: buf(nxs,nzs,2)

    if( mod( it-1, ntdec_s ) /= 0 ) return

    buf(:,:,:) = 0.0

    !$omp parallel do private( ii, kk, i, k )
    do ii = is0, is1
      do kk= ks0, ks1
        k = kk * kdec - kdec/2
        i = ii * idec - idec/2

        buf(ii,kk,1) = Vx(k,i) * UC * M0
        buf(ii,kk,2) = Vz(k,i) * UC * M0

      end do
    end do
    !$omp end parallel do


    if( snp_format == 'native' ) then
      call write_reduce_array2d_r( nxs, nzs, xz_v%ionode, xz_v%io, buf(:,:,1) )
      call write_reduce_array2d_r( nxs, nzs, xz_v%ionode, xz_v%io, buf(:,:,2) )
    else
      call write_reduce_array2d_r_nc( it, 1, nxs, nzs, xz_v, buf(:,:,1) )
      call write_reduce_array2d_r_nc( it, 2, nxs, nzs, xz_v, buf(:,:,2) )
    end if

  end subroutine wbuf_xz_v
  !! --------------------------------------------------------------------------------------------------------------------------- !!

  !! --------------------------------------------------------------------------------------------------------------------------- !!
  subroutine wbuf_xz_u(it)

    integer, intent(in) :: it
    integer :: i, k
    integer :: ii, kk
    !! --

    !$omp parallel do private(ii,kk,i,k)
    do ii = is0, is1
      do kk= ks0, ks1
        k = kk * kdec - kdec/2
        i = ii * idec - idec/2

        buf_u(ii,kk,1) = buf_u(ii,kk,1) + Vx(k,i) * UC * M0 * dt
        buf_u(ii,kk,2) = buf_u(ii,kk,2) + Vz(k,i) * UC * M0 * dt

      end do
    end do
    !$omp end parallel do

    if( mod( it-1, ntdec_s ) == 0 ) then

      if( snp_format == 'native' ) then
        call write_reduce_array2d_r( nxs, nzs, xz_u%ionode, xz_u%io, buf_u(:,:,1) )
        call write_reduce_array2d_r( nxs, nzs, xz_u%ionode, xz_u%io, buf_u(:,:,2) )
      else
        call write_reduce_array2d_r_nc( it, 1, nxs, nzs, xz_u, buf_u(:,:,1) )
        call write_reduce_array2d_r_nc( it, 2, nxs, nzs, xz_u, buf_u(:,:,2) )
      end if

    end if

  end subroutine wbuf_xz_u
  !! --------------------------------------------------------------------------------------------------------------------------- !!

  !! --------------------------------------------------------------------------------------------------------------------------- !!
  !>
  !! write waveform
  !<
  !! ----
  subroutine output__store_wav(it)
    integer, intent(in) :: it
    integer :: i, itw
    real(MP) :: dxVx, dxVz, dzVx, dzVz
    integer :: ii, kk

    if( .not. sw_wav ) return

    call pwatch__on( "output__store_wav" )


    if( it == 1 ) then
      allocate( ux(nst), uz(nst) )
      ux(:) = 0.0
      uz(:) = 0.0

      if( sw_wav_strain ) then
        allocate( exx(nst), ezz(nst), exz(nst) )
        exx(:) = 0.0
        ezz(:) = 0.0
        exz(:) = 0.0
      end if
            
    end if

    !! integrate waveform
    if( sw_wav_u ) then
      !$omp parallel do private(i)
      do i=1, nst
        ux(i) = ux(i) + ( Vx( kst(i), ist(i) ) + Vx( kst(i), ist(i)-1 ) ) * 0.5 * dt
        ! vertical: positive upward for output
        uz(i) = uz(i) - ( Vz( kst(i), ist(i) ) + Vz( kst(i)-1, ist(i) ) ) * 0.5 * dt 
      end do
      !$omp end parallel do
    end if

   if( sw_wav_strain ) then
      !$omp parallel do private(i, ii, kk, dxVx, dzVz, dxVz, dzVx )
      do i=1, nst
        ii = ist(i)
        kk = kst(i)

        dxVx = (Vx(kk  ,ii  ) - Vx(kk  ,ii-1)) * r40x  -  (Vx(kk  ,ii+1) - Vx(kk  ,ii-2)) * r41x
        dzVz = (Vz(kk  ,ii  ) - Vz(kk-1,ii  )) * r40z  -  (Vz(kk+1,ii  ) - Vz(kk-2,ii  )) * r41z

        dxVz = ( (Vz(kk  ,ii+1) - Vz(kk  ,ii  )) * r40x  -  (Vz(kk  ,ii+2) - Vz(kk  ,ii-1)) * r41x &
               + (Vz(kk-1,ii+1) - Vz(kk-1,ii  )) * r40x  -  (Vz(kk-1,ii+2) - Vz(kk-1,ii-1)) * r41x &
               + (Vz(kk  ,ii  ) - Vz(kk  ,ii-1)) * r40x  -  (Vz(kk  ,ii+1) - Vz(kk  ,ii-2)) * r41x &
               + (Vz(kk-1,ii  ) - Vz(kk-1,ii-1)) * r40x  -  (Vz(kk-1,ii+1) - Vz(kk-1,ii-2)) * r41x ) * 0.25
        dzVx = ( (Vx(kk+1,ii  ) - Vx(kk  ,ii  )) * r40z  -  (Vx(kk+2,ii  ) - Vx(kk-1,ii  )) * r41z &
               + (Vx(kk+1,ii-1) - Vx(kk  ,ii-1)) * r40z  -  (Vx(kk+2,ii-1) - Vx(kk-1,ii-1)) * r41z &
               + (Vx(kk  ,ii  ) - Vx(kk-1,ii  )) * r40z  -  (Vx(kk+1,ii  ) - Vx(kk-2,ii  )) * r41z &
               + (Vx(kk  ,ii-1) - Vx(kk-1,ii-1)) * r40z  -  (Vx(kk+1,ii-1) - Vx(kk-2,ii-1)) * r41z ) * 0.25

        exx(i) = exx(i) + dxVx * dt
        ezz(i) = ezz(i) + dzVz * dt
        exz(i) = exz(i) + (dxVz + dzVx) * 0.5 * dt

      end do
      !$omp end parallel do
      
    end if
    
          


    !! output
    if( mod( it-1, ntdec_w ) == 0 ) then
      itw = (it-1)/ntdec_w + 1
      if( sw_wav_v ) then
        !$omp parallel do private(i)
        do i=1, nst
          wav_vel(itw,1,i) =  ( Vx( kst(i), ist(i) ) + Vx( kst(i), ist(i)-1) ) * 0.5 * M0 * UC * 1e9 !! [nm/s]
          wav_vel(itw,2,i) = -( Vz( kst(i), ist(i) ) + Vz( kst(i)-1, ist(i)) ) * 0.5 * M0 * UC * 1e9 !! [nm/s]
        end do
        !$omp end parallel do
      end if

      if( sw_wav_u ) then
        !$omp parallel do private(i)
        do i=1, nst
          wav_disp(itw,1,i) = ux(i) * M0 * UC * 1e9                          !! [nm]
          wav_disp(itw,2,i) = uz(i) * M0 * UC * 1e9                          !! [nm]
        end do
        !$omp end parallel do
      end if

      if( sw_wav_stress ) then
        !$omp parallel do private(i)
        do i=1, nst
          wav_stress(itw,1,i) = Sxx(kst(i),ist(i)) * M0 * UC * 1e6  !! [N/m^2]
          wav_stress(itw,2,i) = Szz(kst(i),ist(i)) * M0 * UC * 1e6
          wav_stress(itw,3,i) = (Sxz(kst(i),  ist(i)) + Sxz(kst(i),  ist(i)-1)  &
                               + Sxz(kst(i)-1,ist(i)) + Sxz(kst(i)-1,ist(i)-1)  ) * 0.25 * M0 * UC * 1e6
        end do
      end if

      if( sw_wav_strain ) then
        !$omp parallel do private(i)
        do i=1, nst
          wav_strain(itw,1,i) = exx(i) * M0 * UC * 1e-3
          wav_strain(itw,2,i) = ezz(i) * M0 * UC * 1e-3
          wav_strain(itw,3,i) = exz(i) * M0 * UC * 1e-3
        end do
        !$omp end parallel do
      end if
      
    end if

    call pwatch__off( "output__store_wav" )

  end subroutine output__store_wav
  !! --------------------------------------------------------------------------------------------------------------------------- !!

  !! --------------------------------------------------------------------------------------------------------------------------- !!
  subroutine output__checkpoint( io )
    integer, intent(in) :: io
    !! ----

    write( io ) xz_ps, xz_v, xz_u
    write( io ) sw_wav, sw_wav_u, sw_wav_v
    write( io ) sw_wav_stress, sw_wav_strain
    write( io ) wav_format
    
    write( io ) ntdec_s
    write( io ) idec, kdec

    write( io ) nxs, nzs
    write( io ) xsnp(1:nxs)
    write( io ) zsnp(1:nzs)

    write( io ) snp_format

    write( io ) is0, is1, ks0, ks1
    write( io ) r20x, r20z

    write( io ) buf_u(is0:is1,ks0:ks1,1:2)

    if( sw_wav ) then
      write( io ) ntdec_w
      write( io ) nst
      write( io ) ntw

      if( nst > 0 ) then
        write( io ) xst(1:nst)
        write( io ) zst(1:nst)
        write( io ) ist(1:nst)
        write( io ) kst(1:nst)
        write( io ) stlo(1:nst)
        write( io ) stla(1:nst)
        write( io ) stnm(1:nst)
        
        if( sw_wav_v  ) then
          write( io ) sh_vel(:,:), wav_vel(:,:,:)
        end if
        
        if( sw_wav_u ) then
          write( io ) sh_disp(:,:), wav_disp(:,:,:)
          write( io ) ux(:), uz(:)
        end if

        if( sw_wav_stress ) then
          write( io ) sh_stress(:,:), wav_stress(:,:,:)
        end if
        
        if( sw_wav_strain ) then
          write( io ) sh_strain(:,:), wav_strain(:,:,:)
          write( io ) exx(:), ezz(:), exz(:)
        end if

      end if

    end if



  end subroutine output__checkpoint
  !! --------------------------------------------------------------------------------------------------------------------------- !!
  !! --------------------------------------------------------------------------------------------------------------------------- !!
  subroutine output__restart( io )
    integer, intent(in) :: io
    !! ----


    read( io ) xz_ps, xz_v, xz_u
    read( io ) sw_wav, sw_wav_u, sw_wav_v
    read( io ) sw_wav_stress, sw_wav_strain
    read( io ) wav_format
    
    read( io ) ntdec_s
    read( io ) idec, kdec

    read( io ) nxs,nzs

    allocate( xsnp(1:nxs) )
    allocate( zsnp(1:nzs) )
    read( io ) xsnp(1:nxs)
    read( io ) zsnp(1:nzs)

    read( io ) snp_format

    read( io ) is0, is1, ks0, ks1
    read( io ) r20x, r20z

    allocate( buf_u(nxs,nzs,2) )
    read( io ) buf_u(is0:is1,ks0:ks1,1:2)

    if( sw_wav ) then
      read( io ) ntdec_w
      read( io ) nst
      read( io ) ntw

      if( nst > 0 ) then
        allocate( xst(nst), zst(nst) )
        allocate( ist(nst), kst(nst) )
        allocate( stnm(nst) )
        allocate( stlo(nst), stla(nst) )

        read( io ) xst(1:nst)
        read( io ) zst(1:nst)
        read( io ) ist(1:nst)
        read( io ) kst(1:nst)
        read( io ) stlo(1:nst)
        read( io ) stla(1:nst)
        read( io ) stnm(1:nst)

        if( sw_wav_v ) then
          allocate( sh_vel(2,nst), wav_vel(ntw,2,nst) )
          read(io) sh_vel, wav_vel
        end if
        
        if( sw_wav_u ) then
          allocate( sh_disp(2,nst), wav_disp(ntw,2,nst), ux(nst), uz(nst) )
          read(io) sh_disp, wav_disp
          read(io) ux, uz
        end if
        
        if( sw_wav_stress  ) then
          allocate( wav_stress(ntw,3,nst), sh_stress(3,nst) )      
          read( io ) sh_stress, wav_stress
        end if

        if( sw_wav_strain ) then
          allocate( wav_strain(ntw,3,nst), sh_strain(3,nst), exx(nst), ezz(nst), exz(nst) )
          read( io ) sh_strain, wav_strain
          read( io ) exx, ezz, exz
        end if

    end if

    end if
    
    
    if( snp_format == 'native' ) then

#ifdef _ES
      if( xz_ps%sw .and. myid == xz_ps%ionode ) then
        open( xz_ps%io, file=trim(odir) //'/'// trim(title) //'.xz.ps.snp', &
            position='append', form='unformatted', status='old' )
      end if

      if( xz_v%sw .and. myid == xz_v%ionode ) then
        open( xz_v%io, file=trim(odir) //'/'// trim(title) //'.xz.v.snp', &
            position='append', form='unformatted', status='old' )
      end if

      if( xz_u%sw .and. myid == xz_u%ionode ) then
        open( xz_u%io, file=trim(odir) //'/'// trim(title) //'.xz.u.snp', &
            position='append', form='unformatted', status='old' )
      end if
#else
      if( xz_ps%sw .and. myid == xz_ps%ionode ) then
        open( xz_ps%io, file=trim(odir) //'/'// trim(title) //'.xz.ps.snp', &
            access='stream', position='append', form='unformatted', status='old' )
      end if

      if( xz_v%sw .and. myid == xz_v%ionode ) then
        open( xz_v%io, file=trim(odir) //'/'// trim(title) //'.xz.v.snp', &
            access='stream', position='append', form='unformatted', status='old' )
      end if

      if( xz_u%sw .and. myid == xz_u%ionode ) then
        open( xz_u%io, file=trim(odir) //'/'// trim(title) //'.xz.u.snp', &
            access='stream', position='append', form='unformatted', status='old' )
      end if
#endif
    else

#ifdef _NETCDF
      if( xz_ps%sw .and. myid == xz_ps%ionode ) then
        call nc_chk( nf90_open( trim(odir) //'/'// trim(title) //'.xz.ps.nc', NF90_WRITE, xz_ps%io ) )
      end if

      if( xz_v%sw .and. myid == xz_v%ionode ) then
        call nc_chk( nf90_open( trim(odir) //'/'// trim(title) //'.xz.v.nc', NF90_WRITE, xz_v%io ) )
      end if

      if( xz_u%sw .and. myid == xz_u%ionode ) then
        call nc_chk( nf90_open( trim(odir) //'/'// trim(title) //'.xz.u.nc', NF90_WRITE, xz_u%io ) )
      end if
#endif
    end if


  end subroutine output__restart
  !! --------------------------------------------------------------------------------------------------------------------------- !!
  subroutine output__closefiles

    integer :: vid

    if( snp_format == 'native' ) then
      if( xz_ps%sw .and. myid == xz_ps%ionode ) close( xz_ps%io )
      if( xz_v%sw .and. myid == xz_v%ionode ) close( xz_v%io )
      if( xz_u%sw .and. myid == xz_u%ionode ) close( xz_u%io )
    else
#ifdef _NETCDF
      if( xz_ps%sw .and. myid == xz_ps%ionode ) then

        ! set max & min for each vars
        call nc_chk( nf90_redef( xz_ps%io ) )
        do vid=1, 2
          call nc_chk( nf90_put_att( xz_ps%io, xz_ps%varid(vid), 'actual_range', (/xz_ps%vmin(vid), xz_ps%vmax(vid)/) ) )
        end do
        call nc_chk( nf90_enddef( xz_ps%io ) )
        call nc_chk( nf90_close( xz_ps%io ) )

      end if

      if( xz_v%sw  .and. myid == xz_v%ionode  ) then

        ! set max & min for each vars
        call nc_chk( nf90_redef( xz_v%io ) )
        do vid=1, 2
          call nc_chk( nf90_put_att( xz_v%io, xz_v%varid(vid), 'actual_range', (/xz_v%vmin(vid), xz_v%vmax(vid)/) ) )
        end do
        call nc_chk( nf90_close( xz_v%io ) )

      end if

      if( xz_u%sw  .and. myid == xz_u%ionode  ) then

        ! set max & min for each vars
        call nc_chk( nf90_redef( xz_u%io ) )
        do vid=1, 2
          call nc_chk( nf90_put_att( xz_u%io, xz_u%varid(vid), 'actual_range', (/xz_u%vmin(vid), xz_u%vmax(vid)/) ) )
        end do
        call nc_chk( nf90_close( xz_u%io ) )

      end if
#endif
    end if

  end subroutine output__closefiles


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
