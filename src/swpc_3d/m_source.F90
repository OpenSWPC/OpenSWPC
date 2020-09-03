!! ----------------------------------------------------------------------------------------------------------------------------- !!
!>
!! seismic source radiation
!!
!! @copyright
!!   Copyright 2013-2020 Takuto Maeda. All rights reserved. This project is released under the MIT license.
!<
!! ----
#include "m_debug.h"
module m_source

  !! dependency
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
  public :: source__checkpoint
  public :: source__restart

  !! local variables
  character(16)         :: stftype                                        !< type of source time function: used in u_source
  integer               :: stf_file_type                                  !< type of stf file format: used in u_source
  integer               :: nsrc                                           !< source grid number inside the node
  integer               :: n_stfprm                                       !< number of parameters for moment-rate function
  real(SP), allocatable :: srcprm(:,:)                                    !< control paramater for moment-rate function at grids
  real(SP), allocatable :: sx(:), sy(:), sz(:)                            !< source location in distance scale
  integer,  allocatable :: isrc(:), jsrc(:), ksrc(:)                      !< source grid location voxel
  real(MP), allocatable :: mxx(:), myy(:), mzz(:), myz(:), mxz(:), mxy(:) !< moment rate
  real(MP), allocatable :: fx(:), fy(:), fz(:)                            !< body force
  real(MP), allocatable :: mo(:)                                          !< moment at grids
  real(MP)              :: dt_dxyz
  character(4)          :: sdep_fit                                       !< 'bd0', 'bd1', ..., 'bd9'


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
    real(SP), allocatable :: mxx_g(:), myy_g(:), mzz_g(:), myz_g(:), mxz_g(:), mxy_g(:), mo_g(:)
    real(SP), allocatable :: fx_g(:), fy_g(:), fz_g(:)
    real(SP), allocatable :: sx_g(:), sy_g(:), sz_g(:)
    integer,  allocatable :: is_g(:), js_g(:), ks_g(:)
    logical,  allocatable :: inside(:)
    real(SP), allocatable :: srcprm_g(:,:)
    character(99) :: fn_stf
    character(6) :: stf_format
    integer :: io
    logical :: green_mode
    character(3) :: sdep0
    !! ----

    call pwatch__on("source__setup")

    !!
    !! Check other modes
    !!
    call readini(io_prm, 'pw_mode', pw_mode, .false. )
    call readini( io_prm, 'green_mode', green_mode, .false. )
    call readini(io_prm, 'bf_mode', bf_mode, .false. )
    call assert( .not. ( pw_mode .and. green_mode ) )

    !!
    !! plane wave mode
    !!
    if( pw_mode ) then
      if( .not. benchmark_mode ) then !! neglect pw_mode for benchmark
        call pw_setup( io_prm )
        nsrc = 0        !! no regular source grid for plane wave mode
        M0 = 1.0 / UC   !! Fictitious scalar moment for output
        call pwatch__off("source__setup")
        return
      end if
    end if

    !!
    !! Green's Function Mode
    !!
    if( green_mode ) then
      if( .not. benchmark_mode ) then
        nsrc = 0
        pw_mode = .false.
        call pwatch__off("source__setup")
        return
      end if
    end if


    !!
    !! obtain number of source grid points and number of parameters for specifying soruce time functions
    !!
    call readini(io_prm, 'fn_stf', fn_stf, '' )
    call readini(io_prm, 'stftype', stftype, 'kupper' )

    if( trim(adjustl(stftype)) == 'scosine' ) stftype = 'cosine'  !! backward compatibility

    call readini(io_prm, 'stf_format', stf_format, 'xym0ij' )
    call readini(io_prm, 'sdep_fit', sdep_fit, 'asis' )


    !! 将来拡張のため、要素震源時間の時間パラメータ数は可変にしてある
    !! 現時点ではすべて開始時間とライズタイムの2種。
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
      pw_mode = .false.
    else
      !! count source-grid file
      call std__getio(io)
      open(io, file=trim(fn_stf), action='read', status='old')
      call std__countline( io, nsrc_g, "#" )  ! exclude comment line which start by "#"
      close(io)
    end if


    !!
    !! temporal memory allocation
    !!
    allocate( is_g(nsrc_g), js_g(nsrc_g), ks_g(nsrc_g) )
    allocate( srcprm_g(n_stfprm, nsrc_g) )
    allocate( sx_g(nsrc_g),  sy_g(nsrc_g),  sz_g(nsrc_g) )
    allocate( inside( nsrc_g ) )
    inside(1:nsrc_g) = .false.

    if( bf_mode ) then

      allocate( fx_g(nsrc_g), fy_g(nsrc_g), fz_g(nsrc_g) )
      call source__grid_bodyforce( fn_stf, stf_format, nsrc_g, n_stfprm, sx_g, sy_g, sz_g, fx_g, fy_g, fz_g, srcprm_g )

    else

      allocate( mo_g(nsrc_g), mxx_g(nsrc_g), myy_g(nsrc_g), mzz_g(nsrc_g), myz_g(nsrc_g), mxz_g(nsrc_g), mxy_g(nsrc_g) )

      !! fixed source parameter for benchmarking
      if( benchmark_mode ) then
        sx_g(1) = 0.0
        sy_g(1) = 0.0
        sz_g(1) = 5.0
        mo_g(1) = 1e15
        mxx_g(1) = 1 / sqrt(3.0)
        myy_g(1) = 1 / sqrt(3.0)
        mzz_g(1) = 1 / sqrt(3.0)
        myz_g(1) = 0.0
        mxz_g(1) = 0.0
        mxy_g(1) = 0.0
        srcprm_g(1,1) = 0.1 ! rupture start time
        srcprm_g(2,1) = 2.0 ! source time duration
        evlo = clon
        evla = clat
        evdp = sz_g(1)
        mxx0 = mxx_g(1)
        myy0 = myy_g(1)
        mzz0 = mzz_g(1)
        myz0 = myz_g(1)
        mxz0 = mxz_g(1)
        mxy0 = mxy_g(1)
      else
        !!
        !! request all source grid information
        !!
        call source__grid_moment( fn_stf, stf_format, nsrc_g, n_stfprm, &
            sx_g, sy_g, sz_g, mo_g, mxx_g, myy_g, mzz_g, myz_g, mxz_g, mxy_g, srcprm_g )
      end if

    end if


    !! check cut-off frequency
    !! currently assumes srcprm(2,:) indicates rise time
    fcut = 0.
    do i=1, nsrc_g
      fcut = max( fcut, 1/srcprm_g(2,i) )
    end do
    fmax = 2 * fcut

    ! count up all magnitude
    if( bf_mode ) then
      M0 = sqrt( sum( fx_g(1:nsrc_g)**2 + fy_g(1:nsrc_g)**2 + fz_g(1:nsrc_g)**2 ) )

      !! 10^3 for unit conversion of UC for body force
      UC = UC * 10**3
    else
      M0 = sum( mo_g(1:nsrc_g) )
    end if

    !!
    !! check if the soruce grid is located inside the node
    !!
    nsrc = 0
    do i=1, nsrc_g

      is_g(i) = x2i( sx_g(i), xbeg, real(dx) )
      js_g(i) = y2j( sy_g(i), ybeg, real(dy) )
      ks_g(i) = z2k( sz_g(i), zbeg, real(dz) )


      !! Count-up source grid *including sleeve area* so that it works even if the source grid is near MPI node boundary
      if( ibeg - 2 <= is_g(i) .and. is_g(i) <= iend + 3 .and. &
          jbeg - 2 <= js_g(i) .and. js_g(i) <= jend + 3 .and. &
          kbeg - 2 <= ks_g(i) .and. ks_g(i) <= kend + 3      ) then

        inside(i) = .true.

        nsrc = nsrc + 1

      end if
    end do


    !!
    !! memory allocation for source grids inside the node
    !!
    allocate( sx(nsrc), sy(nsrc), sz(nsrc) )
    allocate( isrc(nsrc), jsrc(nsrc), ksrc(nsrc) )
    allocate( srcprm(n_stfprm,nsrc) )
    if( bf_mode ) then
      allocate( fx(nsrc), fy(nsrc), fz(nsrc) )
    else
      allocate( mo(nsrc), mxx(nsrc), myy(nsrc), mzz(nsrc), myz(nsrc), mxz(nsrc), mxy(nsrc) )
    end if


    !!
    !! copy source grid information for the current node
    !!
    nn = 0
    do i=1, nsrc_g
      if( inside(i) ) then

        nn = nn + 1

        isrc(nn) = is_g (i)
        jsrc(nn) = js_g (i)
        ksrc(nn) = ks_g (i)
        sx  (nn) = sx_g (i)
        sy  (nn) = sy_g (i)
        sz  (nn) = sz_g (i)

        if( bf_mode ) then
          fx  (nn) = fx_g (i)
          fy  (nn) = fy_g (i)
          fz  (nn) = fz_g (i)
        else
          mo  (nn) = mo_g (i)
          mxx (nn) = mxx_g(i)
          myy (nn) = myy_g(i)
          mzz (nn) = mzz_g(i)
          myz (nn) = myz_g(i)
          mxz (nn) = mxz_g(i)
          mxy (nn) = mxy_g(i)
        end if

        srcprm(:,nn) = srcprm_g(:,i)

        do k=0, NBD
          write(sdep0,'(I1.1)') k
          sdep0 = 'bd' // sdep0

          if( trim(sdep_fit) == sdep0 ) then
            sz(nn) = bddep(isrc(nn), jsrc(nn), k)
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
      call assert( ybeg <= sy(i) .and. sy(i) <= yend )
      call assert( zbeg <= sz(i) .and. sz(i) <= zend )
    end do

    !!
    !! release memroy
    !!
    if( bf_mode ) then
      deallocate( fx_g, fy_g, fz_g )
    else
      deallocate( mxx_g, myy_g, mzz_g, myz_g, mxz_g, mxy_g )
      deallocate( mo_g )
    end if

    deallocate( sx_g, sy_g, sz_g )
    deallocate( is_g, js_g, ks_g )
    deallocate( inside )

    !!
    !! Normalization for numerical stability
    !! Total moment M0 will again be multiplied when export the result
    !!
    if( bf_mode ) then
      fx(1:nsrc) = fx(1:nsrc) / M0
      fy(1:nsrc) = fy(1:nsrc) / M0
      fz(1:nsrc) = fz(1:nsrc) / M0
    else
      mo(1:nsrc) = mo(1:nsrc) / M0
    end if

    !! common grid-related value for stress drip calculation
    dt_dxyz = dt / ( dx*dy*dz )


    call pwatch__off("source__setup")

  end subroutine source__setup
  !! --------------------------------------------------------------------------------------------------------------------------- !!


  !! --------------------------------------------------------------------------------------------------------------------------- !!
  subroutine pw_setup( io_prm )

    integer, intent(in) :: io_prm
    real(SP)  :: pw_ztop
    real(SP)  :: pw_zlen
    character :: pw_ps
    real(SP)  :: pw_strike, pw_dip, pw_rake
    real(SP)  :: xi1, xi2, vp, vs
    integer   :: i, j, k
    real(SP)  :: sd, cd, sf, cf, sl, cl, c2d, c2f
    real(SP)  :: x0, y0, z0, x1, y1, z1, la0, mu0
    real(SP)  :: stf_ii, stf_vx, stf_vy, stf_vz, stf_yz, stf_xz, stf_xy
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

    pw_strike = std__deg2rad( pw_strike )
    pw_dip    = std__deg2rad( pw_dip )
    pw_rake   = std__deg2rad( pw_rake )

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

      do j=jbeg_m, jend_m
        do i=ibeg_m, iend_m
          do k=kbeg_m, kend_m

            vp = sqrt( ( lam(k,i,j) + 2 * mu(k,i,j) ) / rho(k,i,j) )
            if( vp < epsilon(1.0) ) cycle

            x0 = xbeg + (i-0.5) * dx
            y0 = ybeg + (j-0.5) * dy
            z0 = zbeg + (k-0.5) * dz - pw_ztop
            x1 = x0 + dx/2.
            y1 = y0 + dy/2.
            z1 = z0 + dz/2.

            !! source time function evaluated along rotated zeta value at staggered grid points
            stf_ii =  source__momentrate( sd*sf * x0 - sd*cf * y0 + cd * z0,            stftype, 2, (/ 0., pw_zlen /) )
            stf_vx =  source__momentrate( sd*sf * x1 - sd*cf * y0 + cd * z0 + dt/2.*vp, stftype, 2, (/ 0., pw_zlen /) )
            stf_vy =  source__momentrate( sd*sf * x0 - sd*cf * y1 + cd * z0 + dt/2.*vp, stftype, 2, (/ 0., pw_zlen /) )
            stf_vz =  source__momentrate( sd*sf * x0 - sd*cf * y0 + cd * z1 + dt/2.*vp, stftype, 2, (/ 0., pw_zlen /) )
            stf_yz =  source__momentrate( sd*sf * x0 - sd*cf * y1 + cd * z1,            stftype, 2, (/ 0., pw_zlen /) )
            stf_xz =  source__momentrate( sd*sf * x1 - sd*cf * y0 + cd * z1,            stftype, 2, (/ 0., pw_zlen /) )
            stf_xy =  source__momentrate( sd*sf * x1 - sd*cf * y1 + cd * z0,            stftype, 2, (/ 0., pw_zlen /) )

            la0 = lam(k,i,j)
            mu0 = mu(k,i,j)

            vx (k,i,j) = - sd * sf * stf_vx
            vy (k,i,j) =   sd * cf * stf_vy
            vz (k,i,j) = - cd *      stf_vz
            sxx(k,i,j) = - ( la0 + 2 * mu0 * sd * sd * sf * sf ) * stf_ii / vp
            syy(k,i,j) = - ( la0 + 2 * mu0 * sd * sd * cf * cf ) * stf_ii / vp
            szz(k,i,j) = - ( la0 + 2 * mu0 + cd * cd           ) * stf_ii / vp
            syz(k,i,j) =   2 * mu0 * sd * cd * cf * stf_yz / vp
            sxz(k,i,j) = - 2 * mu0 * sd * cd * sf * stf_xz / vp
            sxy(k,i,j) =   2 * mu0 * sd * cd * sf * stf_xy / vp

          end do
        end do
      end do

    else if ( pw_ps == 's' .or. pw_ps == 'S' ) then

      do j=jbeg_m, jend_m
        do i=ibeg_m, iend_m
          do k=kbeg_m, kend_m

            vs = sqrt( mu(k,i,j) / rho(k,i,j) )
            if( vs < epsilon(1.0) ) cycle


            x0 = xbeg + (i-0.5) * dx
            y0 = ybeg + (j-0.5) * dy
            z0 = zbeg + (k-0.5) * dz - pw_ztop
            x1 = x0 + dx/2.
            y1 = y0 + dy/2.
            z1 = z0 + dz/2.

            la0 = lam(k,i,j)
            mu0 = mu(k,i,j)

            !! source time function evaluated along rotated zeta value at staggered grid points
            stf_ii =  source__momentrate( sd*sf * x0 - sd*cf * y0 + cd * z0,            stftype, 2, (/ 0., pw_zlen /) )
            stf_vx =  source__momentrate( sd*sf * x1 - sd*cf * y0 + cd * z0 + dt/2.*vs, stftype, 2, (/ 0., pw_zlen /) )
            stf_vy =  source__momentrate( sd*sf * x0 - sd*cf * y1 + cd * z0 + dt/2.*vs, stftype, 2, (/ 0., pw_zlen /) )
            stf_vz =  source__momentrate( sd*sf * x0 - sd*cf * y0 + cd * z1 + dt/2.*vs, stftype, 2, (/ 0., pw_zlen /) )
            stf_yz =  source__momentrate( sd*sf * x0 - sd*cf * y1 + cd * z1,            stftype, 2, (/ 0., pw_zlen /) )
            stf_xz =  source__momentrate( sd*sf * x1 - sd*cf * y0 + cd * z1,            stftype, 2, (/ 0., pw_zlen /) )
            stf_xy =  source__momentrate( sd*sf * x1 - sd*cf * y1 + cd * z0,            stftype, 2, (/ 0., pw_zlen /) )

            vx (k,i,j) = ( cl*cf + sl*cd*sf ) * stf_vx
            vy (k,i,j) = ( cl*sf - sl*cd*cf ) * stf_vy
            vz (k,i,j) =          -sl*sd      * stf_vz
            sxx(k,i,j) =   2 * mu0 *sd * sf * (cl * cf + sl * cd * sf) * stf_ii / vs
            syy(k,i,j) = - 2 * mu0 *sd * cf * (cl * sf - sl * cd * cf) * stf_ii / vs
            szz(k,i,j) = - 2 * mu0 *cd      *            sl * sd       * stf_ii / vs
            syz(k,i,j) =   mu0 * (cl * cd * sf  - sl * c2d * cf) * stf_yz / vs
            sxz(k,i,j) =   mu0 * (cl * cd * cf  + sl * c2d * sf) * stf_xz / vs
            sxy(k,i,j) = - mu0 * (cl * sd * c2f + 2 * sl * sd * cd * sf * cf) * stf_xy / vs

          end do
        end do
      end do

    else
      call assert( .false. )
    end if

    !! wavelength condition
    i = x2i( (xbeg+xend)/2, xbeg, real(dx) )
    j = y2j( (ybeg+yend)/2, ybeg, real(dy) )
    k = z2k( pw_ztop,       zbeg, real(dz) )
    fcut = 0.0
    fmax = 0.0
    if( ibeg <= i .and. i <= iend .and. jbeg <= j .and. j<= jend ) then
      if( pw_ps == 'p' .or. pw_ps == 'P' ) then
        vp = sqrt( ( lam(k,i,j) + 2 * mu(k,i,j) ) / rho(k,i,j) )
        fcut = vp / pw_zlen
      else
        vs = sqrt( mu(k,i,j) / rho(k,i,j) )
        fcut = vs / pw_zlen
      end if
    end if
    call mpi_allreduce( fcut, fmax, 1, MPI_REAL, MPI_MAX, mpi_comm_world, ierr )
    fcut = fmax
    fmax = fmax * 2.0

  end subroutine pw_setup
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
  subroutine source__grid_moment( fn_stf, stf_format, ns, nprm, sx, sy, sz, mo, mxx, myy, mzz, myz, mxz, mxy, sprm )

    character(*), intent(in) :: fn_stf
    character(*), intent(in) :: stf_format
    integer,  intent(in)     :: ns
    integer,  intent(in)     :: nprm
    real(SP), intent(out)    :: sx(ns), sy(ns), sz(ns)
    real(SP), intent(out)    :: mo(ns), mxx(ns), myy(ns), mzz(ns), myz(ns), mxz(ns), mxy(ns)
    real(SP), intent(out)    :: sprm(1:nprm,1:ns)
    real(SP) :: strike, dip, rake, lon, lat
    integer :: io
    integer :: i
    character(256) :: adum
    integer :: ierr
    real(SP) :: mw
    real(MP) :: D, S
    integer :: is0, js0, ks0
    real(SP), allocatable :: rdum(:)
    integer :: iex
    real(MP):: M0tmp
    !! ----

    call std__getio( io )
    open( io, file = trim(fn_stf), action='read', status='old' )
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
        read( adum,*,iostat=ierr ) sx(i), sy(i), sz(i), sprm(1,i), sprm(2,i), mo(i), &
            mxx(i), myy(i), mzz(i), myz(i), mxz(i), mxy(i)
        call assert( ierr == 0 )

      case( 'xym0dc' )
        read( adum,*,iostat=ierr ) sx(i), sy(i), sz(i), sprm(1,i), sprm(2,i), mo(i), strike, dip, rake
        call assert( ierr == 0 )
        call assert( -360. <= strike .and. strike <= 360. )
        call assert(  -90. <= dip    .and. dip    <= 90.  )
        call assert( -180. <= rake   .and. rake   <= 180. )
        ! use strike angle measured from map azimuth
        call sdr2moment( strike-phi, dip, rake, mxx(i), myy(i), mzz(i), myz(i), mxz(i), mxy(i) )

      case( 'llm0ij' )
        read( adum,*,iostat=ierr ) lon, lat, sz(i), sprm(1,i), sprm(2,i), mo(i), mxx(i), myy(i), mzz(i), myz(i), mxz(i), mxy(i)
        call assert( ierr == 0 )
        call assert( -360. <= lon .and. lon <= 360 )
        call assert(  -90. <= lat .and. lat <=  90 )
        call geomap__g2c( lon, lat, clon, clat, phi, sx(i), sy(i) )

      case( 'llm0dc' )
        read( adum,*,iostat=ierr ) lon, lat, sz(i), sprm(1,i), sprm(2,i), mo(i), strike, dip, rake
        call assert( ierr == 0 )
        call assert( -360. <= lon .and. lon <= 360 )
        call assert(  -90. <= lat .and. lat <=  90 )
        call sdr2moment( strike-phi, dip, rake, mxx(i), myy(i), mzz(i), myz(i), mxz(i), mxy(i) )
        call geomap__g2c( lon, lat, clon, clat, phi, sx(i), sy(i) )

      case( 'xymwij' )
        read( adum,*,iostat=ierr ) sx(i), sy(i), sz(i), sprm(1,i), sprm(2,i), mw, &
            mxx(i), myy(i), mzz(i), myz(i), mxz(i), mxy(i)
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
        call sdr2moment( strike-phi, dip, rake, mxx(i), myy(i), mzz(i), myz(i), mxz(i), mxy(i) )

      case( 'llmwij' )
        read( adum,*,iostat=ierr ) lon, lat, sz(i), sprm(1,i), sprm(2,i), mw, &
            mxx(i), myy(i), mzz(i), myz(i), mxz(i), mxy(i)
        call assert( ierr == 0 )
        call assert( -360. <= lon .and. lon <= 360 )
        call assert(  -90. <= lat .and. lat <=  90 )
        call assert( mw <= 11. ) !! magnitude
        mo(i) = seismic_moment(mw)
        call geomap__g2c( lon, lat, clon, clat, phi, sx(i), sy(i) )

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
        call sdr2moment( strike-phi, dip, rake, mxx(i), myy(i), mzz(i), myz(i), mxz(i), mxy(i) )
        call geomap__g2c( lon, lat, clon, clat, phi, sx(i), sy(i) )

      case( 'xydsdc' )
        read(adum,*,iostat=ierr) sx(i), sy(i), sz(i), sprm(1,i), sprm(2,i), D, S, strike, dip, rake

        call assert( ierr == 0 )
        call assert( -360. <= strike .and. strike <= 360. )
        call assert(  -90. <= dip    .and. dip    <= 90.  )
        call assert( -180. <= rake   .and. rake   <= 180. )

        ! use strike angle measured from map azimuth
        call sdr2moment( strike-phi, dip, rake, mxx(i), myy(i), mzz(i), myz(i), mxz(i), mxy(i) )
        is0 = x2i( sx(i), xbeg, real(dx) )
        js0 = y2j( sy(i), ybeg, real(dy) )
        ks0 = z2k( sz(i), zbeg, real(dz) )

        if( ibeg - 2 <= is0 .and. is0 <= iend + 3 .and. &
            jbeg - 2 <= js0 .and. js0 <= jend + 3 .and. &
            kbeg - 2 <= ks0 .and. ks0 <= kend + 3      ) then
          mo(i) = (1e9 * mu(ks0,is0,js0)) * D * S
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

        call geomap__g2c( lon, lat, clon, clat, phi, sx(i), sy(i) )

        ! use strike angle measured from map azimuth
        call sdr2moment( strike-phi, dip, rake, mxx(i), myy(i), mzz(i), myz(i), mxz(i), mxy(i) )
        is0 = x2i( sx(i), xbeg, real(dx) )
        js0 = y2j( sy(i), ybeg, real(dy) )
        ks0 = z2k( sz(i), zbeg, real(dz) )

        if( ibeg - 2 <= is0 .and. is0 <= iend + 3 .and. &
            jbeg - 2 <= js0 .and. js0 <= jend + 3 .and. &
            kbeg - 2 <= ks0 .and. ks0 <= kend + 3      ) then
          mo(i) = (1e9 * mu(ks0,is0,js0)) * D * S
        else
          mo(i) = 0.
        end if

      case( 'psmeca' )
        read(adum,*,iostat=ierr) lon, lat, sz(i), mzz(i), mxx(i), myy(i), mxz(i), myz(i), mxy(i), iex
        ! reverse sign
        myz(i) = -myz(i)
        mxy(i) = -mxy(i)
        
        call geomap__g2c( lon, lat, clon, clat, phi, sx(i), sy(i) )

        ! moment in dyn-cm
        M0tmp = sqrt( mxx(i)**2 + myy(i)**2 + mzz(i)**2 + 2 * ( mxz(i)**2 + myz(i)**2 + mxy(i)**2 ))/ sqrt(2.0)
        mo(i) = M0tmp * 10.**(iex)
        
        sprm(1,i) = 0.0
        ! 2 x (empirical half-duration) will be a rise time
        sprm(2,i) = 2 * 1.05 * 1e-8 * mo(i)**(1._dp/3._dp)

        ! convert to N-m unit from Dyn-cm
        mo(i) = mo(i) * 1e-7

        ! scale moment tensor components
        mxx(i) = mxx(i) / M0tmp
        myy(i) = myy(i) / M0tmp
        mzz(i) = mzz(i) / M0tmp
        mxz(i) = mxz(i) / M0tmp
        myz(i) = myz(i) / M0tmp
        mxy(i) = mxy(i) / M0tmp

      case default
        call info( 'invalid source type' )
        call assert( .false. )

      end select

      !! remember first record as hypocenter
      if( i==1 ) then
        call geomap__c2g( sx(i), sy(i), clon, clat, phi, evlo, evla )
        sx0 = sx(i)
        sy0 = sy(i)
        evdp = sz(i)
        mxx0 = mxx(i)
        myy0 = myy(i)
        mzz0 = mzz(i)
        myz0 = myz(i)
        mxz0 = mxz(i)
        mxy0 = mxy(i)
        otim = sprm(1,i)
      end if


    end do

    if( stf_format == 'lldsdc' .or. stf_format == 'xydsdc' ) then
      allocate(rdum(ns))
      call mpi_allreduce(mo, rdum, ns, mpi_real, mpi_max, mpi_comm_world, ierr)
      mo(:) = rdum(:)
      deallocate(rdum)
    end if

    close( io )

  end subroutine source__grid_moment
  !! ---------------------------------------------------------------------------------------------------------------------------- !!



  !! --------------------------------------------------------------------------------------------------------------------------- !!
  !>
  !! returns source grid location, moment and mechanism
  !!
  !! @detail
  !! The source__grid subroutine should return
  !!   - number of source grid ns
  !!   - number of source control parameters ( e.g., rupture start time, rise time ... ) nprm
  !!   - source location sx(ns), sy(ns), sz(ns) in the Cartesian coordinate
  !!   - body force terms
  !!   - source control parameters sprm (nprm,ns)
  !<
  !! ----
  subroutine source__grid_bodyforce( fn_stf, stf_format, ns, nprm, sx, sy, sz, fx, fy, fz, sprm )

    character(*), intent(in) :: fn_stf
    character(*), intent(in) :: stf_format
    integer,  intent(in)     :: ns
    integer,  intent(in)     :: nprm
    real(SP), intent(out)    :: sx(ns), sy(ns), sz(ns)
    real(SP), intent(out)    :: fx(ns), fy(ns), fz(ns)
    real(SP), intent(out)    :: sprm(1:nprm,1:ns)
    real(SP) :: strike, dip, rake, lon, lat
    integer :: io
    integer :: i
    character(256) :: adum
    integer :: ierr
    real(SP) :: mw
    character(2) :: stf_coord
    !! ----

    call std__getio( io )
    open( io, file = trim(fn_stf), action='read', status='old' )
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
        read( adum,*,iostat=ierr ) sx(i), sy(i), sz(i), sprm(1,i), sprm(2,i), fx(i), fy(i), fz(i)
        call assert( ierr == 0 )

      case( 'll' )
        read( adum,*,iostat=ierr ) lon, lat, sz(i), sprm(1,i), sprm(2,i), fx(i), fy(i), fz(i)
        call assert( ierr == 0 )
        call assert( -360. <= lon .and. lon <= 360 )
        call assert(  -90. <= lat .and. lat <=  90 )
        call geomap__g2c( lon, lat, clon, clat, phi, sx(i), sy(i) )

      case default
        call info( 'invalid source type' )
        call assert( .false. )

      end select

      !! remember first record as hypocenter
      if( i==1 ) then
        call geomap__c2g( sx(i), sy(i), clon, clat, phi, evlo, evla )
        evdp = sz(i)
        fx0 = fx(i)
        fy0 = fy(i)
        fz0 = fz(i)
        otim = sprm(1,i)
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
    integer  :: ii, jj, kk
    real(MP) :: sdrop
    integer  :: i
    real(SP) :: stime

    !! ----

    if( bf_mode ) return

    call pwatch__on("source__stressdrop")

    t = n2t( it, tbeg, dt )

    do i=1, nsrc

      stime = source__momentrate( t, stftype, n_stfprm, srcprm(:,i) )
      sdrop = mo(i) * stime * dt_dxyz

      ii = isrc(i)
      jj = jsrc(i)
      kk = ksrc(i)

      Sxx(kk,ii,jj) = Sxx(kk,ii,jj) - mxx(i)*sdrop
      Syy(kk,ii,jj) = Syy(kk,ii,jj) - myy(i)*sdrop
      Szz(kk,ii,jj) = Szz(kk,ii,jj) - mzz(i)*sdrop

      Sxy(kk,ii  ,jj  ) = Sxy(kk,ii  ,jj  ) - mxy(i)*sdrop/4
      Sxy(kk,ii  ,jj-1) = Sxy(kk,ii  ,jj-1) - mxy(i)*sdrop/4
      Sxy(kk,ii-1,jj  ) = Sxy(kk,ii-1,jj  ) - mxy(i)*sdrop/4
      Sxy(kk,ii-1,jj-1) = Sxy(kk,ii-1,jj-1) - mxy(i)*sdrop/4

      Sxz(kk  ,ii  ,jj) = Sxz(kk  ,ii  ,jj) - mxz(i)*sdrop/4
      Sxz(kk-1,ii  ,jj) = Sxz(kk-1,ii  ,jj) - mxz(i)*sdrop/4
      Sxz(kk  ,ii-1,jj) = Sxz(kk  ,ii-1,jj) - mxz(i)*sdrop/4
      Sxz(kk-1,ii-1,jj) = Sxz(kk-1,ii-1,jj) - mxz(i)*sdrop/4

      Syz(kk  ,ii,jj  ) = Syz(kk  ,ii,jj  ) - myz(i)*sdrop/4
      Syz(kk-1,ii,jj  ) = Syz(kk-1,ii,jj  ) - myz(i)*sdrop/4
      Syz(kk  ,ii,jj-1) = Syz(kk  ,ii,jj-1) - myz(i)*sdrop/4
      Syz(kk-1,ii,jj-1) = Syz(kk-1,ii,jj-1) - myz(i)*sdrop/4

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

    integer, intent(in) :: it !< time grid
    !! --
    integer :: i, ii, jj, kk
    real(SP) :: t, stime
    !! ----

    if( .not. bf_mode ) return

    call pwatch__on("source__bodyforce")

    t = n2t( it, tbeg, dt ) + dt / 2

    do i=1, nsrc

      stime = source__momentrate( t, stftype, n_stfprm, srcprm(:,i) )

      ii = isrc(i)
      jj = jsrc(i)
      kk = ksrc(i)

      Vx(kk  ,ii  ,jj  ) = Vx(kk  ,ii  ,jj  ) + bx(kk  ,ii  ,jj  ) * fx(i) * stime * dt_dxyz / 2
      Vx(kk  ,ii-1,jj  ) = Vx(kk  ,ii-1,jj  ) + bx(kk  ,ii-1,jj  ) * fx(i) * stime * dt_dxyz / 2
      Vy(kk  ,ii  ,jj  ) = Vy(kk  ,ii  ,jj  ) + by(kk  ,ii  ,jj  ) * fy(i) * stime * dt_dxyz / 2
      Vy(kk  ,ii  ,jj-1) = Vy(kk  ,ii  ,jj-1) + by(kk  ,ii  ,jj-1) * fy(i) * stime * dt_dxyz / 2
      Vz(kk  ,ii  ,jj  ) = Vz(kk  ,ii  ,jj  ) + bz(kk  ,ii  ,jj  ) * fz(i) * stime * dt_dxyz / 2
      Vz(kk-1,ii  ,jj  ) = Vz(kk-1,ii  ,jj  ) + bz(kk-1,ii  ,jj  ) * fz(i) * stime * dt_dxyz / 2

    end do

    call pwatch__off("source__bodyforce")
  end subroutine source__bodyforce
  !! --------------------------------------------------------------------------------------------------------------------------- !!


  !! --------------------------------------------------------------------------------------------------------------------------- !!
  subroutine source__checkpoint( io )
    integer, intent(in) :: io
    !! ----

    write( io ) stftype
    write( io ) stf_file_type
    write( io ) nsrc
    write( io ) n_stfprm
    write( io ) bf_mode

    if( nsrc > 0 ) then
      write( io ) srcprm(1:n_stfprm,1:nsrc)
      write( io ) sx(1:nsrc)
      write( io ) sy(1:nsrc)
      write( io ) sz(1:nsrc)
      write( io ) isrc(1:nsrc)
      write( io ) jsrc(1:nsrc)
      write( io ) ksrc(1:nsrc)

      if( bf_mode ) then
        write( io ) fx(1:nsrc)
        write( io ) fy(1:nsrc)
        write( io ) fz(1:nsrc)
      else
        write( io ) mxx(1:nsrc)
        write( io ) myy(1:nsrc)
        write( io ) mzz(1:nsrc)
        write( io ) myz(1:nsrc)
        write( io ) mxz(1:nsrc)
        write( io ) mxy(1:nsrc)
        write( io ) mo(1:nsrc)
      end if

    end if

    write( io ) dt_dxyz
    write( io ) M0

    if( bf_mode ) then
      write( io ) fx0, fy0, fz0
    else
      write( io ) mxx0, myy0, mzz0, myz0, mxz0, mxy0
    end if

  end subroutine source__checkpoint
  !! --------------------------------------------------------------------------------------------------------------------------- !!

  !! --------------------------------------------------------------------------------------------------------------------------- !!
  subroutine source__restart( io )
    integer, intent(in) :: io
    !! ----

    read( io ) stftype
    read( io ) stf_file_type
    read( io ) nsrc
    read( io ) n_stfprm
    read( io ) bf_mode

    allocate( sx(nsrc), sy(nsrc), sz(nsrc) )
    allocate( isrc(nsrc), jsrc(nsrc), ksrc(nsrc) )
    allocate( srcprm(n_stfprm,nsrc) )
    if( bf_mode ) then
      allocate( fx(nsrc), fy(nsrc), fz(nsrc) )
    else
      allocate( mo(nsrc), mxx(nsrc), myy(nsrc), mzz(nsrc), myz(nsrc), mxz(nsrc), mxy(nsrc) )
    end if

    if( nsrc > 0 ) then
      read( io ) srcprm(1:n_stfprm,1:nsrc)
      read( io ) sx(1:nsrc)
      read( io ) sy(1:nsrc)
      read( io ) sz(1:nsrc)
      read( io ) isrc(1:nsrc)
      read( io ) jsrc(1:nsrc)
      read( io ) ksrc(1:nsrc)

      if( bf_mode ) then
        read( io ) fx(1:nsrc)
        read( io ) fy(1:nsrc)
        read( io ) fz(1:nsrc)
      else
        read( io ) mxx(1:nsrc)
        read( io ) myy(1:nsrc)
        read( io ) mzz(1:nsrc)
        read( io ) myz(1:nsrc)
        read( io ) mxz(1:nsrc)
        read( io ) mxy(1:nsrc)
        read( io ) mo(1:nsrc)
      end if
    end if

    read( io ) dt_dxyz
    read( io ) M0

    if( bf_mode ) then
      read( io ) fx0, fy0, fz0
    else
      read( io ) mxx0, myy0, mzz0, myz0, mxz0, mxy0
    end if


  end subroutine source__restart
  !! --------------------------------------------------------------------------------------------------------------------------- !!



  !! ---------------------------------------------------------------------------------------------------------------------------- !!
  !>
  !! returns number of source grids
  !<
  !! ----
  real(SP) function source__momentrate( t, stftype, nprm, srcprm )

    real(SP), intent(in) :: t
    character(*),  intent(in) :: stftype
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
  !! --------------------------------------------------------------------------------------------------------------------------- !!

end module m_source
!! ----------------------------------------------------------------------------------------------------------------------------- !!
