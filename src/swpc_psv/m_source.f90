!! ----------------------------------------------------------------------------------------------------------------------------- !!
!>
!! seismic source radiation
!!
!! @copyright
!!   Copyright 2013-2024 Takuto Maeda. All rights reserved. This project is released under the MIT license.
!<
!! --
#include "../shared/m_debug.h"
module m_source

  !! dependency
  use iso_fortran_env, only: error_unit
  use m_std
  use m_debug
  use m_global
  use m_fdtool
  use m_pwatch
  use m_readini
  use m_geomap
  use mpi

  !! declarations
  implicit none
  private
  save

  !! public routines
  public :: source__setup
  public :: source__stressdrop
  public :: source__bodyforce

  !! local variables
  character(16)         :: stftype                          !< type of source time function: used in u_source
  integer               :: stf_file_type                    !< type of stf file format: used in u_source
  integer               :: nsrc                             !< source grid number inside the node
  integer               :: n_stfprm                         !< number of parameters for moment-rate function
  real(SP), allocatable :: srcprm(:,:)                      !< control paramater for moment-rate function at grids
  real(SP), allocatable :: sx(:), sz(:)                     !< source location in distance scale
  integer,  allocatable :: isrc(:), ksrc(:)                 !< source grid location voxel
  real(MP), allocatable :: mxx(:), mzz(:), mxz(:)           !< moment rate
  real(MP), allocatable :: fx(:), fz(:)
  real(MP), allocatable :: mo(:)                            !< moment at grids
  real(MP)              :: dt_dxz
  character(4)          :: sdep_fit                         !< 'bd0', 'bd1', ..., 'bd9'
  logical :: earth_flattening


contains

  !! --------------------------------------------------------------------------------------------------------------------------- !!
  !>
  !! set up module
  !!
  !! @detail
  !! read source grid from user-module u_source, and allocates memory for source grid and time-related parameters
  !! also pre-calculate common variables for stress drop source
  !!
  !<
  subroutine source__setup( io_prm )

    integer, intent(in) :: io_prm
    !! ----

    integer :: nn, i, k
    integer :: nsrc_g
    real(SP), allocatable :: mxx_g(:), mzz_g(:), mxz_g(:), mo_g(:)
    real(SP), allocatable :: fx_g(:), fz_g(:)
    real(SP), allocatable :: sx_g(:), sz_g(:)
    integer,  allocatable :: is_g(:), ks_g(:)
    logical,  allocatable :: inside(:)
    real(SP), allocatable :: srcprm_g(:,:)
    character(99) :: fn_stf
    character(6) :: stf_format
    integer :: io
    character(3) :: sdep0
    !! ----

    call pwatch__on("source__setup")


    !!
    !! plane wave mode
    !!
    call readini(io_prm, 'pw_mode', pw_mode, .false. )
    if( pw_mode ) then
      if( .not. benchmark_mode ) then !! neglect pw_mode for benchmark
        call pw_setup( io_prm )
        nsrc = 0        !! no regular source grid for plane wave mode
        M0 = 1.0 / UC   !! Fictitious scalar moment for output
        call pwatch__off("source__setup")
        return
      end if
    end if

    !! body-force mode?
    call readini(io_prm, 'bf_mode', bf_mode, .false. )

    !!
    !! obtain number of source grid points and number of parameters for specifying soruce time functions
    !!
    call readini( io_prm, 'fn_stf', fn_stf, '' )
    call readini( io_prm, 'stftype', stftype, 'kupper' )
    call readini( io_prm, 'stf_format', stf_format, 'xym0ij' )
    call readini( io_prm, 'sdep_fit', sdep_fit, 'asis' )

    call readini( io_prm, 'earth_flattening', earth_flattening, .false. )


    if( trim(adjustl(stftype)) == 'scosine' ) stftype = 'cosine'  !! backward compatibility

    !! n_stfprm may be changed for future extension.
    select case ( trim(stftype) )
    case ( 'boxcar'   );  n_stfprm = 2
    case ( 'triangle' );  n_stfprm = 2
    case ( 'herrmann' );  n_stfprm = 2
    case ( 'kupper'   );  n_stfprm = 2
    case ( 'cosine'   );  n_stfprm = 2
    case ( 'texp'     );  n_stfprm = 2
    end select



    !! fixed source parameter for benchmarking: no source grid file used
    if( benchmark_mode ) then
      nsrc_g   = 1
      stftype  = 'kupper'
      n_stfprm = 2
      bf_mode = .false.
    else

      !! count source-grid file
      open( newunit=io, file=trim(fn_stf), action='read', status='old')
      call std__countline( io, nsrc_g, "#" )  ! exclude comment line which start by "#"
      close( io )
    end if


    !!
    !! temporal memory allocation
    !!
    allocate( is_g(nsrc_g), ks_g(nsrc_g) )
    allocate( srcprm_g(n_stfprm, nsrc_g) )
    allocate( sx_g(nsrc_g), sz_g(nsrc_g) )
    allocate( inside( nsrc_g ) )
    inside(1:nsrc_g) = .false.

    if( bf_mode ) then
      allocate( fx_g(nsrc_g), fz_g(nsrc_g) )
      call source__grid_bodyforce( fn_stf, stf_format, nsrc_g, n_stfprm, sx_g, sz_g, fx_g, fz_g, srcprm_g )
    else

      allocate( mo_g(nsrc_g), mxx_g(nsrc_g), mzz_g(nsrc_g), mxz_g(nsrc_g) )
      !! fixed source parameter for benchmarking
      if( benchmark_mode ) then
        sx_g(1) = 0.0
        sz_g(1) = 5.0
        mo_g(1) = 1e15
        mxx_g(1) = 1 / sqrt(2.0)
        mzz_g(1) = 1 / sqrt(2.0)
        mxz_g(1) = 0.0
        srcprm_g(1,1) = 0.1
        srcprm_g(2,1) = 2.0

      else

        !!
        !! request all source grid information ( from user-defined routine )
        !!
        call source__grid_moment( fn_stf, stf_format, nsrc_g, n_stfprm, &
            sx_g, sz_g, mo_g, mxx_g, mzz_g, mxz_g, srcprm_g)

      end if
    end if

    ! 
    if (earth_flattening) then
      do k=1, nsrc_g
        sz_g(k) = - R_EARTH * log( ( R_EARTH - sz_g(k) ) / R_EARTH )
      end do
    end if

    ! count up all moment / body force
    if( bf_mode ) then
      M0 = sqrt( sum( fx_g(1:nsrc_g)**2 + fz_g(1:nsrc_g)**2  ) )
      !! unit conversion for body wave source
      UC = UC * 10**3
    else
      M0 = sum( mo_g(1:nsrc_g) )
    end if


    !! check cut-off frequency
    !! currently assumes srcprm(2,:) indicates rise time
    fcut = 0.
    do i=1, nsrc_g
      fcut = max( fcut, 1/srcprm_g(2,i) )
    end do
    fmax = 2 * fcut

    !!
    !! check if the soruce grid is located inside the node
    !!
    nsrc = 0
    do i=1, nsrc_g

      is_g(i) = x2i( sx_g(i), xbeg, real(dx) )
      ks_g(i) = z2k( sz_g(i), zbeg, real(dz) )


      !! Count-up source grid *including sleeve area* so that it works even if the source grid is near MPI node boundary
      if( ibeg - 2 <= is_g(i) .and. is_g(i) <= iend + 3 .and. &
          kbeg - 2 <= ks_g(i) .and. ks_g(i) <= kend + 3      ) then

        inside(i) = .true.

        nsrc = nsrc + 1

      end if
    end do


    !!
    !! memory allocation for source grids inside the node
    !!
    allocate( sx(nsrc), sz(nsrc) )
    allocate( isrc(nsrc), ksrc(nsrc) )
    allocate( srcprm(n_stfprm,nsrc) )

    if( bf_mode ) then
      allocate( fx(nsrc), fz(nsrc) )
    else
      allocate( mo(nsrc), mxx(nsrc), mzz(nsrc), mxz(nsrc) )
    end if


    !!
    !! copy source grid information for the current node
    !!
    nn = 0
    do i=1, nsrc_g
      if( inside(i) ) then

        nn = nn + 1

        isrc(nn) = is_g (i)
        ksrc(nn) = ks_g (i)
        sx  (nn) = sx_g (i)
        sz  (nn) = sz_g (i)

        if( bf_mode ) then
          fx  (nn) = fx_g (i)
          fz  (nn) = fz_g (i)
        else
          mo  (nn) = mo_g (i)
          mxx (nn) = mxx_g(i)
          mzz (nn) = mzz_g(i)
          mxz (nn) = mxz_g(i)
        end if

        srcprm(:,nn) = srcprm_g(:,i)

        !! depth fitting

        do k=0, NBD
          write(sdep0,'(I1.1)') k
          sdep0 = 'bd' // sdep0

          if( trim(sdep_fit) == sdep0 ) then
            sz(nn) = bddep(isrc(nn),  k)
            ksrc(nn) = z2k( sz(nn), zbeg, real(dz) )
          end if
        end do

      end if
    end do

    !!
    !! Confirm that FDM model space contains the srouce grid location
    !!
    do i=1, nsrc
      call assert( xbeg <= sx(i) .and. sx(i) <= xend )
      call assert( zbeg <= sz(i) .and. sz(i) <= zend )
    end do


    !!
    !! release memroy
    !!
    deallocate( sx_g, sz_g )
    deallocate( is_g, ks_g )
    deallocate( inside )
    if( bf_mode ) then
      deallocate( fx_g, fz_g )
    else
      deallocate( mxx_g, mzz_g, mxz_g )
      deallocate( mo_g )
    end if

    !!
    !! Normalization for numerical stability
    !! Total moment M0 will again be multiplied when export the result
    !!
    if( bf_mode ) then
      fx(1:nsrc) = fx(1:nsrc) / M0
      fz(1:nsrc) = fz(1:nsrc) / M0
    else
      mo(1:nsrc) = mo(1:nsrc) / M0
    end if

    !! common grid-related value for stress drip calculation
    dt_dxz = dt / ( dx*dz )


    call pwatch__off("source__setup")

  end subroutine source__setup
  !! --------------------------------------------------------------------------------------------------------------------------- !!


  !! --------------------------------------------------------------------------------------------------------------------------- !!
  !>
  !! returns source grid location, moment and mechanism
  !!
  !! @detail
  !! The source__grid subroutine should return
  !!   - number of source grid ns
  !!   - number of source control parameters ( e.g., rupture start time, rise time ... ) nprm
  !!   - source location sx(ns), sy(ns), sz(ns) in the Cartesian coordinate
  !!   - moment release at source grids mo(ns)
  !!   - moment tensors mij(ns) (six components)
  !!   - source control parameters sprm (nprm,ns)
  !!
  !<
  !! ----
  subroutine source__grid_moment( fn_stf, stf_format, ns, nprm, sx, sz, mo, mxx, mzz, mxz, sprm )

    character(*), intent(in) :: fn_stf
    character(*),  intent(in) :: stf_format
    integer,  intent(in)     :: ns
    integer,  intent(in)     :: nprm
    real(SP), intent(out)    :: sx(ns), sz(ns)
    real(SP), intent(out)    :: mo(ns), mxx(ns), mzz(ns), mxz(ns)
    real(SP), intent(out)    :: sprm(1:nprm,1:ns)
    real(SP) :: strike, dip, rake, lon, lat
    integer :: io
    integer :: i
    character(256) :: adum
    real(SP) :: sy(ns)
    integer :: ierr
    real(SP) :: rdum
    real(SP) :: mw
    real(SP) :: D, S
    integer :: is0, js0, ks0
    real(SP), allocatable :: r0(:)
    real(MP):: M0tmp
    integer :: iex
    !---
        
    open( newunit=io, file = trim(fn_stf), action='read', status='old' )
    i = 0

    do

      read(io,'(A256)', iostat=ierr) adum
      adum = adjustl(adum)
      if( ierr /= 0 ) exit                  ! detect EOF
      if( adjustl(adum(1:1)) == "#" ) cycle ! neglect comment line
      if( trim(adjustl(adum)) == "" ) cycle ! neglect blank line

      i = i+1

      select case( stf_format )

      case( 'xym0ij' )
        read( adum,*,iostat=ierr ) sx(i), sy(i), sz(i), sprm(1,i), sprm(2,i), mo(i), mxx(i), rdum, mzz(i), rdum, mxz(i), rdum
        call assert( ierr == 0 )

      case( 'xym0dc' )
        read( adum,*,iostat=ierr ) sx(i), sy(i), sz(i), sprm(1,i), sprm(2,i), mo(i), strike, dip, rake
        call assert( ierr == 0 )
        call assert( -360. <= strike .and. strike <= 360. )
        call assert(  -90. <= dip    .and. dip    <= 90.  )
        call assert( -180. <= rake   .and. rake   <= 180. )
        ! use strike angle measured from map azimuth
        call sdr2moment( strike-phi, dip, rake, mxx(i), rdum, mzz(i), rdum, mxz(i), rdum )

      case( 'llm0ij' )
        read( adum,*,iostat=ierr ) lon, lat, sz(i), sprm(1,i), sprm(2,i), mo(i), mxx(i), rdum, mzz(i), rdum, mxz(i), rdum
        call assert( ierr == 0 )
        call assert( -360. <= lon .and. lon <= 360 )
        call assert(  -90. <= lat .and. lat <=  90 )
        call geomap__g2c( lon, lat, clon, clat, phi, sx(i), sy(i))

      case( 'llm0dc' )
        read( adum,*,iostat=ierr ) lon, lat, sz(i), sprm(1,i), sprm(2,i), mo(i), strike, dip, rake
        call assert( ierr == 0 )
        call assert( -360. <= lon .and. lon <= 360 )
        call assert(  -90. <= lat .and. lat <=  90 )
        call sdr2moment( strike-phi, dip, rake, mxx(i), rdum, mzz(i), rdum, mxz(i), rdum )
        call geomap__g2c( lon, lat, clon, clat, phi, sx(i), sy(i) )

      case( 'xymwij' )
        read( adum,*,iostat=ierr ) sx(i), sy(i), sz(i), sprm(1,i), sprm(2,i), mw, mxx(i), rdum, mzz(i), rdum, mxz(i), rdum
        call assert( ierr == 0 )
        call assert( mw <= 11. ) !! magnitude
        mo(i) = seismic_moment(mw)

      case( 'xymwdc' )
        read( adum,*,iostat=ierr ) sx(i), sy(i), sz(i), sprm(1,i), sprm(2,i), mw, strike, dip, rake
        call assert( ierr == 0 )
        call assert( -360. <= strike .and. strike <= 360. )
        call assert(  -90. <= dip    .and. dip    <= 90.  )
        call assert( -180. <= rake   .and. rake   <= 180. )
        call assert( mw <= 11. ) !! magnitude
        ! use strike angle measured from map azimuth
        mo(i) = seismic_moment(mw)
        call sdr2moment( strike-phi, dip, rake, mxx(i), rdum, mzz(i), rdum, mxz(i), rdum )

      case( 'llmwij' )
        read( adum,*,iostat=ierr ) lon, lat, sz(i), sprm(1,i), sprm(2,i), mw, mxx(i), rdum, mzz(i), rdum, mxz(i), rdum
        call assert( ierr == 0 )
        call assert( -360. <= lon .and. lon <= 360 )
        call assert(  -90. <= lat .and. lat <=  90 )
        call assert( mw <= 11. ) !! magnitude
        mo(i) = seismic_moment(mw)
        call geomap__g2c( lon, lat, clon, clat, phi, sx(i), sy(i))

      case( 'llmwdc' )
        read( adum,*,iostat=ierr ) lon, lat, sz(i), sprm(1,i), sprm(2,i), mw, strike, dip, rake
        call assert( ierr == 0 )
        call assert( -360. <= lon .and. lon <= 360 )
        call assert(  -90. <= lat .and. lat <=  90 )
        call assert( -360. <= strike .and. strike <= 360. )
        call assert(  -90. <= dip    .and. dip    <= 90.  )
        call assert( -180. <= rake   .and. rake   <= 180. )
        call assert( mw <= 11. ) !! magnitude
        mo(i) = seismic_moment(mw)
        call sdr2moment( strike-phi, dip, rake, mxx(i), rdum, mzz(i), rdum, mxz(i), rdum )
        call geomap__g2c( lon, lat, clon, clat, phi, sx(i), sy(i) )

      case( 'xydsdc' )
        read(adum,*,iostat=ierr) sx(i), sy(i), sz(i), sprm(1,i), sprm(2,i), D, S, strike, dip, rake
        call assert( ierr == 0 )
        call assert( -360. <= strike .and. strike <= 360. )
        call assert(  -90. <= dip    .and. dip    <= 90.  )
        call assert( -180. <= rake   .and. rake   <= 180. )
        call sdr2moment( strike-phi, dip, rake, mxx(i), rdum, mzz(i), rdum, mxz(i), rdum )
        is0 = x2i( sx(i), xbeg, real(dx) )

        if( earth_flattening ) then
          ks0 = z2k( real( - R_EARTH * log( ( R_EARTH - sz(i) )/R_EARTH )), zbeg, real(dz) )
        else
          ks0 = z2k( sz(i), zbeg, real(dz) )
        end if

        if( ibeg - 2 <= is0 .and. is0 <= iend + 3 .and. &
            kbeg - 2 <= ks0 .and. ks0 <= kend + 3      ) then
          mo(i) = (1e9 * mu(ks0,js0)) * D * S
        else
          mo(i) = 0.
        end if

      case( 'lldsdc' )
        read(adum,*,iostat=ierr) lon, lat, sz(i), sprm(1,i), sprm(2,i), D, S, strike, dip, rake
        call assert( ierr == 0 )
        call assert( -360. <= lon .and. lon <= 360 )
        call assert(  -90. <= lat .and. lat <=  90 )
        call assert( -360. <= strike .and. strike <= 360. )
        call assert(  -90. <= dip    .and. dip    <= 90.  )
        call assert( -180. <= rake   .and. rake   <= 180. )

        call sdr2moment( strike-phi, dip, rake, mxx(i), rdum, mzz(i), rdum, mxz(i), rdum )
        call geomap__g2c( lon, lat, clon, clat, phi, sx(i), sy(i) )

        is0 = x2i( sx(i), xbeg, real(dx) )
        if( earth_flattening ) then
          ks0 = z2k( real( - R_EARTH * log( ( R_EARTH - sz(i) )/R_EARTH )), zbeg, real(dz) )
        else
          ks0 = z2k( sz(i), zbeg, real(dz) )
        end if

        if( ibeg - 2 <= is0 .and. is0 <= iend + 3 .and. &
            kbeg - 2 <= ks0 .and. ks0 <= kend + 3      ) then
          mo(i) = (1e9 * mu(ks0,js0)) * D * S
        else
          mo(i) = 0.
        end if    

      case( 'psmeca' )
        read(adum,*,iostat=ierr) lon, lat, sz(i), mzz(i), mxx(i), rdum, mxz(i), rdum, rdum, iex
        
        call geomap__g2c( lon, lat, clon, clat, phi, sx(i), sy(i) )

        ! moment in dyn-cm
        M0tmp = sqrt( mxx(i)**2 + mzz(i)**2 + 2 * ( mxz(i)**2 ))/ sqrt(2.0)
        mo(i) = M0tmp * 10.**(iex)
        
        sprm(1,i) = 0.0
        ! 2 x (empirical half-duration) will be a rise time
        sprm(2,i) = 2 * 1.05 * 1e-8 * mo(i)**(1._dp/3._dp)

        ! convert to N-m unit from Dyn-cm
        mo(i) = mo(i) * 1e-7

        ! scale moment tensor components
        mxx(i) = mxx(i) / M0tmp
        mzz(i) = mzz(i) / M0tmp
        mxz(i) = mxz(i) / M0tmp

      case default
        write(error_unit,*) "ERROR [source__setup]: Invalid source type: "//trim(stf_format)

      end select

      !! remember first record as hypocenter
      if( i==1 ) then
        call geomap__c2g( sx(i), sy(i), clon, clat, phi, evlo, evla )
        sx0 = sx(i)
        sy0 = sy(i)
        evdp = sz(i)
        mxx0 = mxx(i)
        myy0 = -12345.0
        mzz0 = mzz(i)
        myz0 = -12345.0
        mxz0 = mxz(i)
        mxy0 = -12345.0
        otim = sprm(1,i)
      end if

    end do

    close( io )

    if( stf_format == 'lldsdc' .or. stf_format == 'xydsdc' ) then
      allocate(r0(ns))
      call mpi_allreduce(mo, r0, ns, mpi_real, mpi_max, mpi_comm_world, ierr)
      mo(:) = r0(:)
      deallocate(r0)
    end if


  end subroutine source__grid_moment
  !! --------------------------------------------------------------------------------------------------------------------------- !!


  !! --------------------------------------------------------------------------------------------------------------------------- !!
  !>
  !! returns source grid location, moment and mechanism
  !!
  !! @detail
  !! The source__grid subroutine should return
  !!   - number of source grid ns
  !!   - number of source control parameters ( e.g., rupture start time, rise time ... ) nprm
  !!   - source location sx(ns), sy(ns), sz(ns) in the Cartesian coordinate
  !!   - moment release at source grids mo(ns)
  !!   - moment tensors mij(ns) (six components)
  !!   - source control parameters sprm (nprm,ns)
  !!
  !<
  !! ----
  subroutine source__grid_bodyforce( fn_stf, stf_format, ns, nprm, sx, sz, fx, fz, sprm )

    character(*), intent(in) :: fn_stf
    character(*),  intent(in) :: stf_format
    integer,  intent(in)     :: ns
    integer,  intent(in)     :: nprm
    real(SP), intent(out)    :: sx(ns), sz(ns)
    real(SP), intent(out)    :: fx(ns), fz(ns)
    real(SP), intent(out)    :: sprm(1:nprm,1:ns)
    real(SP) :: strike, dip, rake, lon, lat
    integer :: io
    integer :: i
    character(256) :: adum
    integer :: ierr
    real(SP) :: rdum
    real(SP) :: mw
    character(2) :: stf_coord
    !! --

    open( newunit=io, file = trim(fn_stf), action='read', status='old' )
    i = 0

    stf_coord(1:2) = stf_format(1:2)

    do

      read(io,'(A256)', iostat=ierr) adum
      if( ierr /= 0 ) exit                  ! detect EOF
      if( adjustl(adum(1:1)) == "#" ) cycle ! neglect comment line
      if( trim(adjustl(adum)) == "" ) cycle ! neglect blank line

      i = i+1

      select case( stf_coord )

      case( 'xy' )
        read( adum,*,iostat=ierr ) sx(i), rdum, sz(i), sprm(1,i), sprm(2,i), fx(i), rdum, fz(i)
        call assert( ierr == 0 )

      case( 'll' )
        read( adum,*,iostat=ierr ) lon, lat, sz(i), sprm(1,i), sprm(2,i), fx(i), rdum, fz(i)
        call assert( ierr == 0 )
        call assert( -360. <= lon .and. lon <= 360 )
        call assert(  -90. <= lat .and. lat <=  90 )
        call geomap__g2c( lon, lat, clon, clat, phi, sx(i), rdum )

      case default
        write(error_unit,*) "ERROR [source__setup]: Invalid source type: "//trim(stf_format)

      end select

      !! remember first record as hypocenter
      if( i==1 ) then
        call geomap__c2g( sx(i), 0.0, clon, clat, phi, evlo, evla )
        evdp = sz(i)
        otim = sprm(1,i)
        fx0 = fx(i)
        fy0 = -12345.0
        fz0 = fz(i)
      end if


    end do

    close( io )

  end subroutine source__grid_bodyforce
  !! ---------------------------------------------------------------------------------------------------------------------------- !!



  !! --------------------------------------------------------------------------------------------------------------------------- !!
  !>
  !! apply stress drop for source grids
  !<
  !!
  subroutine source__stressdrop( it )

    integer, intent(in) :: it !< time grid number

    !! ----

    real(SP) :: t
    integer  :: ii, kk
    real(MP) :: sdrop
    integer  :: i
    real(MP) :: stime

    !! ----


    if( bf_mode ) return

    call pwatch__on("source__stressdrop")

    t = n2t( it, tbeg, dt )

    do i=1, nsrc

      stime = source__momentrate( t, stftype, n_stfprm, srcprm(:,i) )
      sdrop = mo(i) * stime * dt_dxz

      ii = isrc(i)
      kk = ksrc(i)

      Sxx(kk,ii) = Sxx(kk,ii) - mxx(i)*sdrop
      Szz(kk,ii) = Szz(kk,ii) - mzz(i)*sdrop

      Sxz(kk  ,ii  ) = Sxz(kk  ,ii  ) - mxz(i)*sdrop/4
      Sxz(kk-1,ii  ) = Sxz(kk-1,ii  ) - mxz(i)*sdrop/4
      Sxz(kk  ,ii-1) = Sxz(kk  ,ii-1) - mxz(i)*sdrop/4
      Sxz(kk-1,ii-1) = Sxz(kk-1,ii-1) - mxz(i)*sdrop/4

    end do

    call pwatch__off("source__stressdrop")

  end subroutine source__stressdrop
  !! --------------------------------------------------------------------------------------------------------------------------- !!



  !! --------------------------------------------------------------------------------------------------------------------------- !!
  !>
  !! apply body force for source grids
  !<
  !!
  subroutine source__bodyforce( it )

    integer, intent(in) :: it !< time grid number

    !! ----

    real(SP) :: t
    integer  :: ii, kk
    real(MP) :: sdrop
    integer  :: i
    real(MP) :: stime

    !! ----

    if( .not. bf_mode ) return

    call pwatch__on("source__bodyforce")

    t = n2t( it, tbeg, dt ) + dt / 2

    do i=1, nsrc

      stime = source__momentrate( t, stftype, n_stfprm, srcprm(:,i) )

      ii = isrc(i)
      kk = ksrc(i)


      Vx(kk  ,ii  ) = Vx(kk  ,ii  ) + fx(i) * stime * dt_dxz / 2
      Vx(kk  ,ii-1) = Vx(kk  ,ii-1) + fx(i) * stime * dt_dxz / 2
      Vz(kk  ,ii  ) = Vz(kk  ,ii  ) + fz(i) * stime * dt_dxz / 2
      Vz(kk-1,ii  ) = Vz(kk-1,ii  ) + fz(i) * stime * dt_dxz / 2

    end do

    call pwatch__off("source__bodyforce")

  end subroutine source__bodyforce
  !! --------------------------------------------------------------------------------------------------------------------------- !!



  !!
  !! user defined
  !!

  !! ---------------------------------------------------------------------------------------------------------------------------- !!
  !>
  !! returns number of source grids
  !<
  !! ----
  real(SP) function source__momentrate( t, stftype, nprm, srcprm )

    real(SP), intent(in) :: t
    character(*), intent(in) :: stftype
    integer,  intent(in) :: nprm
    real(SP), intent(in) :: srcprm(1:nprm)

    !! ----

    real(SP) :: tbeg
    real(SP) :: trise

    !! ----

    tbeg   = srcprm(1)
    trise  = srcprm(2)

    select case ( trim(stftype) )
    case ( 'boxcar'   );  source__momentrate = boxcar   ( t, tbeg, trise )
    case ( 'triangle' );  source__momentrate = triangle ( t, tbeg, trise )
    case ( 'herrmann' );  source__momentrate = herrmann ( t, tbeg, trise )
    case ( 'kupper'   );  source__momentrate = kupper   ( t, tbeg, trise )
    case ( 'cosine'   );  source__momentrate = cosine   ( t, tbeg, trise )
    case ( 'texp'     );  source__momentrate = texp     ( t, tbeg, trise )
    case default;         source__momentrate = kupper   ( t, tbeg, trise )
    end select

  end function source__momentrate
  !! ---------------------------------------------------------------------------------------------------------------------------- !!

  !! --------------------------------------------------------------------------------------------------------------------------- !!
  subroutine pw_setup( io_prm )

    integer, intent(in) :: io_prm
    real(SP)  :: pw_ztop
    real(SP)  :: pw_zlen
    character :: pw_ps
    real(SP)  :: pw_strike, pw_dip, pw_rake
    real(SP)  :: vp, vs
    integer   :: i, k
    real(SP)  :: sd, cd, sf, cf, sl, cl, c2d, c2f
    real(SP)  :: x0, z0, x1, z1, la0, mu0
    real(SP)  :: stf_ii, stf_vx, stf_vz, stf_xz
    integer   :: ierr
    !! ----

    !!
    !! parameter input
    !!
    call readini( io_prm, 'pw_ztop', pw_ztop, 1e30 )
    call assert( pw_ztop < zend )

    call readini( io_prm, 'pw_zlen', pw_zlen, -1.)
    call assert( pw_zlen > 0.0 )

    call readini( io_prm, 'pw_ps', pw_ps, '' )
    call assert( pw_ps=='p' .or. pw_ps=='s' .or. pw_ps=='P' .or. pw_ps=='S' )

    call readini( io_prm, 'pw_strike', pw_strike, 0.0 )
    call readini( io_prm, 'pw_dip', pw_dip, 0.0 )
    call readini( io_prm, 'pw_rake', pw_rake, 0. )

    pw_strike = std__deg2rad( 90.0 )
    pw_dip    = std__deg2rad( pw_dip )
    pw_rake   = std__deg2rad( 90.0 )

    call readini(io_prm, 'stftype', stftype, 'kupper' )


    sd  = sin( pw_dip )
    cd  = cos( pw_dip )
    sf  = sin( pw_strike )
    cf  = cos( pw_strike )
    sl  = sin( pw_rake )
    cl  = cos( pw_rake )
    c2d = cos( 2 * pw_dip     )
    c2f = cos( 2 * pw_strike )

    if( pw_ps == 'p' .or. pw_ps == 'P' ) then

      do i=ibeg_m, iend_m
        do k=kbeg_m, kend_m

          vp = sqrt( ( lam(k,i) + 2 * mu(k,i) ) / rho(k,i) )
          if( vp < epsilon(1.0) ) cycle

          x0 = xbeg + (i-0.5) * dx
          z0 = zbeg + (k-0.5) * dz - pw_ztop
          x1 = x0 + dx/2.
          z1 = z0 + dz/2.

          !! source time function evaluated along rotated zeta value at staggered grid points
          stf_ii =  source__momentrate( sd*sf * x0 + cd * z0,            stftype, 2, (/ 0., pw_zlen /) )
          stf_vx =  source__momentrate( sd*sf * x1 + cd * z0 + dt/2.*vp, stftype, 2, (/ 0., pw_zlen /) )
          stf_vz =  source__momentrate( sd*sf * x0 + cd * z1 + dt/2.*vp, stftype, 2, (/ 0., pw_zlen /) )
          stf_xz =  source__momentrate( sd*sf * x1 + cd * z1,            stftype, 2, (/ 0., pw_zlen /) )

          la0 = lam(k,i)
          mu0 = mu(k,i)

          vx (k,i) = - sd * sf * stf_vx
          vz (k,i) = - cd *      stf_vz
          sxx(k,i) = - ( la0 + 2 * mu0 * sd * sd * sf * sf ) * stf_ii / vp
          szz(k,i) = - ( la0 + 2 * mu0 + cd * cd           ) * stf_ii / vp
          sxz(k,i) = - 2 * mu0 * sd * cd * sf * stf_xz / vp

        end do
      end do

    else if ( pw_ps == 's' .or. pw_ps == 'S' ) then

      do i=ibeg_m, iend_m
        do k=kbeg_m, kend_m

          vs = sqrt( mu(k,i) / rho(k,i) )
          if( vs < epsilon(1.0) ) cycle


          x0 = xbeg + (i-0.5) * dx
          z0 = zbeg + (k-0.5) * dz - pw_ztop
          x1 = x0 + dx/2.
          z1 = z0 + dz/2.

          la0 = lam(k,i)
          mu0 = mu(k,i)

          !! source time function evaluated along rotated zeta value at staggered grid points
          stf_ii =  source__momentrate( sd*sf * x0 + cd * z0,            stftype, 2, (/ 0., pw_zlen /) )
          stf_vx =  source__momentrate( sd*sf * x1 + cd * z0 + dt/2.*vs, stftype, 2, (/ 0., pw_zlen /) )
          stf_vz =  source__momentrate( sd*sf * x0 + cd * z1 + dt/2.*vs, stftype, 2, (/ 0., pw_zlen /) )
          stf_xz =  source__momentrate( sd*sf * x1 + cd * z1,            stftype, 2, (/ 0., pw_zlen /) )

          vx (k,i) = ( cl*cf + sl*cd*sf ) * stf_vx
          vz (k,i) =          -sl*sd      * stf_vz
          sxx(k,i) =   2 * mu0 *sd * sf * (cl * cf + sl * cd * sf) * stf_ii / vs
          szz(k,i) = - 2 * mu0 *cd      *            sl * sd       * stf_ii / vs
          sxz(k,i) =   mu0 * (cl * cd * cf  + sl * c2d * sf) * stf_xz / vs

        end do
      end do

    else
      call assert( .false. )
    end if

    !! wavelength condition
    i = x2i( (xbeg+xend)/2, xbeg, real(dx) )
    k = z2k( pw_ztop, zbeg, real(dz) )
    fcut = 0
    fmax = 0
    if( ibeg <= i .and. i <= iend ) then

      if( pw_ps == 'p' .or. pw_ps == 'P' ) then
        vp = sqrt( ( lam(k,i) + 2 * mu(k,i) ) / rho(k,i) )
        fcut = vp / pw_zlen
      else
        vs = sqrt( mu(k,i) / rho(k,i) )
        fcut = vs / pw_zlen
      end if
    end if

    call mpi_allreduce( fcut, fmax, 1, MPI_REAL, MPI_MAX, mpi_comm_world, ierr )
    fcut = fmax
    fmax = fcut * 2.0

  end subroutine pw_setup
  !! --------------------------------------------------------------------------------------------------------------------------- !!


end module m_source
!! ----------------------------------------------------------------------------------------------------------------------------- !!
