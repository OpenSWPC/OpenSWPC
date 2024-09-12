!! ----------------------------------------------------------------------------------------------------------------------------- !!
!>
!! Snapshot output
!!
!! @copyright
!!   Copyright 2013-2024 Takuto Maeda. All rights reserved. This project is released under the MIT license.
!<
!! ----
#include "../shared/m_debug.h"
module m_output

  !! -- Dependency
  use iso_fortran_env, only: error_unit
  use m_std
  use m_debug
  use m_global
  use m_pwatch
  use m_fdtool
  use m_daytim
  use m_readini
  use m_geomap
  use netcdf

  !! -- Declarations
  implicit none
  private
  save

  public :: snap__setup
  public :: snap__write
  public :: snap__closefiles

  !! -- Internal Parameters
  character(8), parameter :: BINARY_TYPE= "STREAMIO"
  !! character(8), parameter :: BINARY_TYPE= "UNFORMAT" ! for architectures that does not support fortran2003
  character(8), parameter :: CODE_TYPE  = "SWPC_SH "   !!< FIXED parameter for file header
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
    logical :: sw
    integer :: io
    integer :: ionode
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

  type(snp) :: xz_v, xz_u

  !! switch
  integer   :: ntdec_s                                !< time step decimation factor: Snap
  integer   :: idec, kdec                                !< spatial decimation factor: x, y, z

  integer :: nxs, nzs !< snapshot grid size
  real(SP), allocatable :: xsnp(:), zsnp(:)

  !! I/O area in the node
  integer :: is0, is1, ks0, ks1
  character(6) :: snp_format ! native or netcdf

  real(SP), allocatable :: buf_u(:,:) !! displacement buffer


contains

  !! --------------------------------------------------------------------------------------------------------------------------- !!
  !>
  !! Setup
  !!
  !! @see
  !! #2012-41 #2013-00440
  !<
  !! ----
  subroutine snap__setup( io_prm )

    integer, intent(in) :: io_prm
    integer :: i, k, ii, kk
    !! ----

    call pwatch__on( "snap__setup" )

    !!
    !! read parameter
    !!

    call readini( io_prm, 'xz_v%sw',  xz_v%sw,  .false. )
    call readini( io_prm, 'xz_u%sw',  xz_u%sw,  .false. )
    call readini( io_prm, 'idec',      idec,      1       )
    call readini( io_prm, 'kdec',      kdec,      1       )
    call readini( io_prm, 'ntdec_s',   ntdec_s,   10      )
    call readini( io_prm, 'snp_format', snp_format, 'native' )

    !! snapshot size #2013-0440
    nxs = ( nx + (idec/2) ) / idec
    nzs = ( nz + (kdec/2) ) / kdec

    !! coordinate
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

    !! number of snapshots per cycle
    xz_v  % nsnp = 1
    xz_u  % nsnp = 1

    xz_v  % snaptype  = 'vy'
    xz_u  % snaptype  = 'uy'

    xz_v  % vname(1) = 'Vy'
    xz_v  % vunit(1) = 'm/s'
    xz_u  % vname(1) = 'Uy'
    xz_u  % vunit(1) = 'm'

    !! output node definition: cyclic
    xz_v%ionode = mod( 1, nproc_x )
    xz_u%ionode = mod( 2, nproc_x )

    !! output settings
    if( snp_format == 'native' ) then
      if( xz_v%sw ) call newfile_xz(   trim(odir) // '/' // trim(title) //'.xz.v.snp', xz_v )
      if( xz_u%sw ) call newfile_xz(   trim(odir) // '/' // trim(title) //'.xz.u.snp', xz_u )
    else
      if( xz_v%sw ) call newfile_xz_nc(   trim(odir) // '/' // trim(title) //'.xz.v.nc', xz_v )
      if( xz_u%sw ) call newfile_xz_nc(   trim(odir) // '/' // trim(title) //'.xz.u.nc', xz_u )
    end if

    allocate(buf_u(nxs,nzs))
    buf_u(:,:) = 0.0

    call pwatch__off( "snap__setup" )

  end subroutine snap__setup

  !! --------------------------------------------------------------------------------------------------------------------------- !!
  !>
  !! Open new file and write header information, and medium parameters for XY-cross section
  !<
  !! ----
  subroutine newfile_xz ( fname, hdr )

    character(*), intent(in)    :: fname
    type(snp),    intent(inout) :: hdr
    !! ----

    integer :: i, k, kk, ii
    real :: buf(nxs, nzs,3)

    !! --

    hdr % na1 = na / idec
    hdr % na2 = na / kdec
    hdr % coordinate = 'xz'


    if( myid == hdr%ionode ) then

      open( newunit=hdr%io, file=trim(fname), access='stream', action='write', form='unformatted', status='replace' )
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
    write( hdr%io ) dt * ntdec_s                   ! dt

    ! absorb layer
    write( hdr%io ) hdr%na1
    write( hdr%io ) hdr%na2

    !! number of arrays
    write( hdr%io ) hdr%nmed
    write( hdr%io ) hdr%nsnp

    !! coordinate
    dum = -12345.0
    write( hdr%io ) clon
    write( hdr%io ) clat
    write( hdr%io ) phi

    write( hdr%io ) dum !< dummy
    write( hdr%io ) dum !< dummy
    write( hdr%io ) dum !< dummy

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

  end subroutine write_nc_header
  !! --------------------------------------------------------------------------------------------------------------------------- !!

  !! --------------------------------------------------------------------------------------------------------------------------- !!
  subroutine newfile_xz_nc ( fname, hdr )
    character(*), intent(in)    :: fname
    type(snp),    intent(inout) :: hdr
    !!
    real, allocatable :: sbuf(:), rbuf1(:), rbuf2(:), rbuf3(:), buf(:,:,:)
    integer :: i, k, ierr, ii, kk
    !! ---

    if( myid == hdr%ionode ) then

      !! initialize
      hdr%vmax = 0.0
      hdr%vmin = 0.0
      hdr % na1 = na / idec
      hdr % na2 = na / kdec

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

  end subroutine newfile_xz_nc
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
  !! write snapshot
  !<
  !! ----
  subroutine snap__write( it )

    integer, intent(in) :: it


    call pwatch__on( "snap__write" )

    if( xz_v%sw ) call wbuf_xz_v(it)
    if( xz_u%sw ) call wbuf_xz_u(it)

    call pwatch__off( "snap__write" )

  end subroutine snap__write


  !! --------------------------------------------------------------------------------------------------------------------------- !!
  subroutine wbuf_xz_v(it)
    integer :: it
    integer :: i, k
    integer :: ii, kk
    real(SP), allocatable, save :: buf(:,:), sbuf(:), rbuf(:)
    integer, save :: req
    integer, save :: it0
    integer :: stat(mpi_status_size)
    integer :: err
    !! ----

    if( .not. allocated(buf) ) allocate(buf(nxs,nzs), source=0.0)

    if( (mod( it-1, ntdec_s ) == 0) .or. ( it > nt) ) then
        buf(:,:) = 0.0
        !$omp parallel do private(ii,kk,i,k)
        do ii = is0, is1
            do kk= ks0, ks1
                k = kk * kdec - kdec/2
                i = ii * idec - idec/2

                buf(ii,kk) = Vy(k,i) * UC * M0

            end do
        end do
        !$omp end parallel do

        if( snp_format == 'native' ) then
            call write_reduce_array2d_r( nxs, nzs, xz_v%ionode, xz_v%io, buf )
        else
            if (.not. allocated(sbuf)) then
                allocate(sbuf(nxs*nzs), rbuf(nxs*nzs))
            else
                call mpi_wait(req, stat, err)
                if( myid == xz_v%ionode ) call wbuf_nc(xz_v, 1, nxs, nzs, it0, rbuf)
            end if
            if( it <= nt ) then ! except for the last call
                sbuf = reshape(buf(:,:), (/nxs * nzs/))
                call mpi_ireduce(sbuf, rbuf, nxs * nzs , mpi_real, mpi_sum, xz_v%ionode, mpi_comm_world, req, err)
                it0 = it ! remember
            end if    
        end if
    end if

  end subroutine wbuf_xz_v
  !! --------------------------------------------------------------------------------------------------------------------------- !!

  !! --------------------------------------------------------------------------------------------------------------------------- !!
  subroutine wbuf_xz_u(it)
    integer :: it
    integer :: i, k
    integer :: ii, kk
    real, allocatable, save :: sbuf(:), rbuf(:)
    integer, save :: req
    integer, save :: it0
    integer :: stat(mpi_status_size)
    integer :: err       
    !! ----

    !$omp parallel do private(ii, kk, i, k)
    do ii = is0, is1
      do kk= ks0, ks1
        k = kk * kdec - kdec/2
        i = ii * idec - idec/2

        buf_u(ii,kk) = buf_u(ii,kk) +  Vy(k,i) * UC * M0 * dt

      end do
    end do
    !$omp end parallel do

    if( mod( it-1, ntdec_s ) == 0 .or. ( it > nt) ) then

      if( snp_format == 'native' ) then
        call write_reduce_array2d_r( nxs, nzs, xz_u%ionode, xz_u%io, buf_u )
      else
        if (.not. allocated(sbuf)) then
            allocate(sbuf(nxs*nzs), rbuf(nxs*nzs))
        else
            call mpi_wait(req, stat, err)
            if( myid == xz_u%ionode ) call wbuf_nc(xz_u, 1, nxs, nzs, it0, rbuf)
        end if
        if( it <= nt ) then ! except for the last call
            sbuf = reshape(buf_u(:,:), (/nxs * nzs/))
            call mpi_ireduce(sbuf, rbuf, nxs * nzs, mpi_real, mpi_sum, xz_u%ionode, mpi_comm_world, req, err)
            it0 = it ! remember
        end if

    end if

    end if

  end subroutine wbuf_xz_u

  subroutine wbuf_nc(hdr, nvar, nx1, nx2, it0, rbuf)
    type(snp), intent(inout) :: hdr
    integer,   intent(in) :: nvar
    integer,   intent(in) :: nx1, nx2
    integer,   intent(in) :: it0
    real,      intent(in) :: rbuf(:)

    integer :: ns, ib, ie
    integer :: stt(3), cnt(3)
    integer :: vid
    !! ----

    ns = nx1 * nx2
    cnt = (/ nx1, nx2, 1/)
    stt = (/ 1, 1, it0 / ntdec_s + 1 /)
    do vid=1, nvar
        ib = (vid-1)*ns + 1
        ie = ib + ns - 1

        call nc_chk(nf90_put_var(hdr%io, hdr%varid(vid), reshape(rbuf(ib:ie), (/nx1, nx2/)), stt, cnt))
        call nc_chk( nf90_put_var( hdr%io, hdr%vid_t, it0*dt, start=(/ it0/ntdec_s+1 /) ) )
        hdr%vmax(vid) = max(hdr%vmax(vid), maxval(rbuf(ib:ie)))
        hdr%vmin(vid) = min(hdr%vmin(vid), minval(rbuf(ib:ie)))
        call nc_chk(nf90_redef( hdr%io ))
        call nc_chk(nf90_put_att( hdr%io, hdr%varid(vid), 'actual_range', (/hdr%vmin(vid), hdr%vmax(vid)/)))
        call nc_chk(nf90_enddef( hdr%io ))
!        call nc_chk(nf90_sync(hdr%io))
    end do

  end subroutine wbuf_nc  

  subroutine snap__closefiles

    call pwatch__on('snap__closefiles')

    if( snp_format == 'native' ) then
      if( xz_v%sw .and. myid == xz_v%ionode ) close( xz_v%io )
      if( xz_u%sw .and. myid == xz_u%ionode ) close( xz_u%io )
    else

    if( xz_v%sw   ) call wbuf_xz_v(nt+1)
    if( xz_u%sw   ) call wbuf_xz_u(nt+1)

    if( xz_v%sw .and. myid == xz_v%ionode ) call close_nc(xz_v)
    if( xz_u%sw .and. myid == xz_u%ionode ) call close_nc(xz_u)

    end if

    call pwatch__off('snap__closefiles')

  end subroutine snap__closefiles

  subroutine close_nc( hdr )
    type(snp), intent(in) :: hdr
    integer :: vid

    call nc_chk( nf90_redef( hdr%io ) )
    do vid = 1, hdr%nsnp
      call nc_chk( nf90_put_att( hdr%io, hdr%varid(vid), 'actual_range', (/hdr%vmin(vid), hdr%vmax(vid)/)) )
    end do
    call nc_chk( nf90_enddef( hdr%io ) )
    call nc_chk( nf90_sync( hdr%io ))
    call nc_chk( nf90_close( hdr%io ) )

  end subroutine close_nc  

  !! --------------------------------------------------------------------------------------------------------------------------- !!
  !>
  !! An internal subroutine to check error in netcdf function calls
  !<
  !! --
  subroutine nc_chk( ierr )

    integer, intent(in) :: ierr
    !! ----

    if( ierr /= NF90_NOERR )  write(error_unit,*) NF90_STRERROR( ierr )

  end subroutine nc_chk
  !! --------------------------------------------------------------------------------------------------------------------------- !!

end module m_output
!! ----------------------------------------------------------------------------------------------------------------------------- !!
