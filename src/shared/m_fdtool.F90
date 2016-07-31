!! ----------------------------------------------------------------------------------------------------------------------------- !!
!>
!! individual utility subroutines / functions working without global-variables/parameters
!!
!! @copyright
!!   Copyright 2013-2016 Takuto Maeda. All rights reserved. This project is released under the MIT license.
!<
!! ----
module m_fdtool

  use m_std
  implicit none

  save
  public

contains


  !! --------------------------------------------------------------------------------------------------------------------------- !!
  subroutine vcheck( vp, vs, rho, xi, vmin, vmax, rhomin, is_vmin_under, is_vmax_over, is_rhomin_under)

    real(SP), intent(inout) :: vp, vs, rho
    real(SP), intent(in) :: xi !< velocity purturbation
    real(SP), intent(in) :: vmin, vmax, rhomin
    logical, intent(inout) :: is_vmin_under, is_vmax_over, is_rhomin_under
    real :: gamma, xi2

    !! keep vp/vs ratio
    gamma = vp / vs

    !! velocity check
    if( vp > vmax .or. vs > vmax ) then
       is_vmax_over = .true.

       xi2 = ( 1 + xi ) * vmax / vp - 1
       vp = vmax
       vp = vmax
       vs = vmax / gamma
       rho = rho * ( 1 + 0.8 * xi2 ) / ( 1 + 0.8 * xi )
    end if

    if( vp < vmin .or. vs < vmin ) then
       is_vmin_under = .true.
       vs = vmin
       vp = vmin * gamma
    end if

    if( rho < rhomin ) then
      is_rhomin_under = .true.
      rho = rhomin
    end if

  end subroutine vcheck
  !! --------------------------------------------------------------------------------------------------------------------------- !!

  !! --------------------------------------------------------------------------------------------------------------------------- !!
  !>
  !! gives maximum time step size by means of stable condition, based on 3D 4th-order assumption
  !!
  !<
  !! ----
  subroutine fdm_stable_dt( dx, dy, dz, vmax, dt )
    real(SP), intent(in) :: dx, dy, dz
    real(SP), intent(in) :: vmax
    real(SP), intent(out) :: dt

    real(SP) :: hh ! effective grid size
    real(SP) :: cc ! coefficients from FDM formula

    !! ----

    hh = 1. / sqrt( 1/dx**2 + 1/dy**2 + 1/dz**2 )
    cc = 6. / 7.
    dt = cc * hh / vmax

  end subroutine fdm_stable_dt
  !! --------------------------------------------------------------------------------------------------------------------------- !!


  !! --------------------------------------------------------------------------------------------------------------------------- !!
  !>
  !! Returns stablity condition c: Must be c<=1
  !!
  !<
  !! ----
  subroutine fdm_cond_stability( dx, dy, dz, vmax, dt, c )
    real(SP), intent(in) :: dx, dy, dz
    real(SP), intent(in) :: vmax
    real(SP), intent(in) :: dt
    real(SP), intent(out) :: c
    !!
    real(SP) :: dt2
    !! ----

    call fdm_stable_dt( dx, dy, dz, vmax, dt2 )
    c = dt / dt2

  end subroutine fdm_cond_stability

  !! --------------------------------------------------------------------------------------------------------------------------- !!
  !>
  !! Returns rg: ratio between minimum-wavelength and grid size from given grid size, minimum velocity and maximum frequency
  !! Better to be rg >= 5-10
  !<
  !!
  subroutine fdm_cond_wavelength( dx, dy, dz, vmin, fmax, r )

    real(SP), intent(in) :: dx, dy, dz
    real(SP), intent(in) :: vmin
    real(SP), intent(in) :: fmax
    real(SP), intent(out) :: r

    real(SP) :: lambda_min
    real(SP) :: dh

    dh = max( dx, dy, dz )
    lambda_min = vmin / fmax
    r = lambda_min / dh

  end subroutine fdm_cond_wavelength



  !! --------------------------------------------------------------------------------------------------------------------------- !!
  !>
  !! Memory-size rough estimate of FDM simulation
  !!
  !<
  subroutine memory_size( nproc_x, nproc_y, nx, ny, nz, nm, na, memsize_all, memsize_node )

    integer,  intent(in)  :: nproc_x
    integer,  intent(in)  :: nproc_y
    integer,  intent(in)  :: nx
    integer,  intent(in)  :: ny
    integer,  intent(in)  :: nz
    integer,  intent(in)  :: nm
    integer,  intent(in)  :: na
    real(SP), intent(out) :: memsize_all  !< total memory size in GB
    real(SP), intent(out) :: memsize_node !< average in-node memory size in GB

    !! ----
    real(SP) :: rn_int, rn_abc
    real(SP) :: b_med, b_vel, b_stress, b_mem, b_abc
    integer :: nproc
    real(SP) :: ba_int, ba_abc, ba_com
    integer :: nxpm, nypm, nzpm, nxp, nyp
    real(SP) :: rn_com
    !! ----

    nproc = nproc_x * nproc_y

    ! node inside memory buffer size
    nxp = ceiling( nx / real( nproc_x ) )
    nxpm = nxp + 6
    nyp = ceiling( ny / real( nproc_y ) )
    nypm = nyp + 6
    nzpm = nz  + 6

    ! mesh number ( in real, for avoiding integer overflow )
    rn_com = real( nxpm * nproc_x ) * real( nypm * nproc_y ) * nzpm
    rn_int = real( nx-2*na ) * real( ny-2*na ) * real( nz-na )
    rn_abc = real( nx )      * real( ny )      * real( nz )    - real( nx-2*nxp ) * real( ny-2*nyp ) * real( nz-na )

    ! byte per mesh
    b_med    = 11 * SP     ! rho, lam, mu, taup, taus, muyz, muxz, muyz, bx, by, bz
    b_vel    = 3 * DP      ! Vx, Vy, Vz
    b_stress = 6 * DP      ! Sxx, Syy, Szz, Syz, Sxz, Sxy
    b_mem    = 6 * nm * SP ! Rij

    b_abc = 18 * SP

    ba_com = ( b_med + b_vel + b_stress ) * rn_com
    ba_int = ( b_mem ) * rn_int
    ba_abc = ( b_abc ) * rn_abc


    memsize_all  = ( ba_com + ba_int + ba_abc ) / ( 1024. ) ** 3
    memsize_node = memsize_all / nproc

  end subroutine memory_size
  !! --------------------------------------------------------------------------------------------------------------------------- !!

  !! --------------------------------------------------------------------------------------------------------------------------- !!
  !>
  !! Memory-size estimate of FDM simulation
  !!
  !<
  subroutine memory_size_sh( nproc, nx, nz, nm, na, memsize_all, memsize_node )

    integer,  intent(in)  :: nproc
    integer,  intent(in)  :: nx
    integer,  intent(in)  :: nz
    integer,  intent(in)  :: nm
    integer,  intent(in)  :: na
    real(SP), intent(out) :: memsize_all  !< total memory size in GB
    real(SP), intent(out) :: memsize_node !< average in-node memory size in GB

    !! ----
    real(SP) :: rn_int, rn_abc
    real(SP) :: b_med, b_vel, b_stress, b_mem, b_abc
    real(SP) :: ba_int, ba_abc, ba_com
    integer  :: nxpm, nzpm, nxp
    real(SP) :: rn_com
    !! ----

    ! node inside memory buffer size
    nxp = ceiling( nx / real( nproc ) )
    nxpm = nxp + 6
    nzpm = nz  + 6

    ! mesh number ( in real, for avoiding integer overflow )
    rn_com = real( nxpm * nproc ) * nzpm
    rn_int = real( nx-2*na ) * real( nz-na )
    rn_abc = real( nx )      * real( nz )    - real( nx-2*nxp ) * real( nz-na )

    ! byte per mesh
    b_med    = 5 * SP
    b_vel    = 1 * DP
    b_stress = 2 * DP
    b_mem    = 2 * nm * SP

    b_abc = 2 * DP

    ba_com = ( b_med + b_vel + b_stress ) * rn_com
    ba_int = ( b_mem ) * rn_int
    ba_abc = ( b_abc ) * rn_abc


    memsize_all  = ( ba_com + ba_int + ba_abc ) / ( 1024. ) ** 3
    memsize_node = memsize_all / nproc

  end subroutine memory_size_sh
  !! --------------------------------------------------------------------------------------------------------------------------- !!

  !! --------------------------------------------------------------------------------------------------------------------------- !!
  !>
  !! Memory-size estimate of FDM simulation
  !!
  !<
  subroutine memory_size_psv( nproc, nx, nz, nm, na, memsize_all, memsize_node )

    integer,  intent(in)  :: nproc
    integer,  intent(in)  :: nx
    integer,  intent(in)  :: nz
    integer,  intent(in)  :: nm
    integer,  intent(in)  :: na
    real(SP), intent(out) :: memsize_all  !< total memory size in GB
    real(SP), intent(out) :: memsize_node !< average in-node memory size in GB

    !! ----
    real(SP) :: rn_int, rn_abc
    real(SP) :: b_med, b_vel, b_stress, b_mem, b_abc
    real(SP) :: ba_int, ba_abc, ba_com
    integer  :: nxpm, nzpm, nxp
    real(SP) :: rn_com
    !! ----

    ! node inside memory buffer size
    nxp = ceiling( nx / real( nproc ) )
    nxpm = nxp + 6
    nzpm = nz  + 6

    ! mesh number ( in real, for avoiding integer overflow )
    rn_com = real( nxpm * nproc ) * nzpm
    rn_int = real( nx-2*na ) * real( nz-na )
    rn_abc = real( nx )      * real( nz )    - real( nx-2*nxp ) * real( nz-na )

    ! byte per mesh
    b_med    = 5 * SP
    b_vel    = 2 * DP
    b_stress = 4 * DP
    b_mem    = 4 * nm * SP

    b_abc = 4 * DP

    ba_com = ( b_med + b_vel + b_stress ) * rn_com
    ba_int = ( b_mem ) * rn_int
    ba_abc = ( b_abc ) * rn_abc


    memsize_all  = ( ba_com + ba_int + ba_abc ) / ( 1024. ) ** 3
    memsize_node = memsize_all / nproc

  end subroutine memory_size_psv
  !! --------------------------------------------------------------------------------------------------------------------------- !!

  !! --------------------------------------------------------------------------------------------------------------------------- !!
  !>
  !! Calculate moment-magnituide from moment m0
  !<
  real(SP) function moment_magnitude ( m0 )

    real(SP), intent(in) :: m0

    !! ----

    if( m0 < epsilon(1.0) ) then
      moment_magnitude = -12345.0 !< undefined
    else
      moment_magnitude = (log10( m0 ) - 9.1 ) * 2.0 / 3.0
    end if
    
  end function moment_magnitude
  !! --------------------------------------------------------------------------------------------------------------------------- !!



  !! --------------------------------------------------------------------------------------------------------------------------- !!
  !>
  !! Calculate seismic moment from moment-magnituide mw
  !<
  real(SP) function seismic_moment ( mw )

    real(SP), intent(in) :: mw

    !! ----

    seismic_moment = 10**( 1.5 * mw + 9.1 )

  end function seismic_moment
  !! --------------------------------------------------------------------------------------------------------------------------- !!



  !! --------------------------------------------------------------------------------------------------------------------------- !!
  !>
  !! Calcluate six moment tensor components from strike, dip, and rake [deg] under the double-couple assumptions
  !<
  subroutine sdr2moment( strike, dip, rake, mxx, myy, mzz, myz, mxz, mxy )

    !! Arguments
    real(SP), intent(in)  :: strike, dip, rake
    real(SP), intent(out) :: mxx, myy, mzz
    real(SP), intent(out) :: mxy, myz, mxz

    !! Locals
    real(SP) :: sind, cosd, sin2d, cos2d, cosl, sinl
    real(SP) :: sinf, cosf, sin2f, cos2f

    !! ----

    sind  = sin( std__deg2rad(      dip     ) )
    cosd  = cos( std__deg2rad(      dip     ) )
    sin2d = sin( std__deg2rad(  2 * dip     ) )
    cos2d = cos( std__deg2rad(  2 * dip     ) )
    sinl  = sin( std__deg2rad(      rake    ) )
    cosl  = cos( std__deg2rad(      rake    ) )
    sinf  = sin( std__deg2rad(      strike  ) )
    cosf  = cos( std__deg2rad(      strike  ) )
    sin2f = sin( std__deg2rad(  2 * strike  ) )
    cos2f = cos( std__deg2rad(  2 * strike  ) )

    mxx = - ( sind*cosl*sin2f + sin2d*sinl*sinf*sinf )
    mxy =   ( sind*cosl*cos2f + sin2d*sinl*sin2f/2   )
    mxz = - ( cosd*cosl*cosf  + cos2d*sinl*sinf      )
    myy =   ( sind*cosl*sin2f - sin2d*sinl*cosf*cosf )
    myz = - ( cosd*cosl*sinf  - cos2d*sinl*cosf      )
    mzz =   (                   sin2d*sinl           )

  end subroutine sdr2moment
  !! --------------------------------------------------------------------------------------------------------------------------- !!


  !! --------------------------------------------------------------------------------------------------------------------------- !!
  !>
  !! Single-lobed Kupper function for moment rate
  !!
  !! This function have >0 value among ts <= t <= ts + tr for given start time ts and rise time tr
  !! The amplitude is normalized by \int_0^\infty f(t) dt = 1
  !<
  real(SP) function kupper( t, ts, tr )

    real(SP), intent(in) :: t  !< time
    real(SP), intent(in) :: ts !< rupture start time
    real(SP), intent(in) :: tr !< rise time

    !! ----

    if ( ts <= t .and. t <= ts + tr ) then
       kupper = 3 * Pi * ( sin( Pi*(t-ts)/tr ) )**3 / ( 4 * tr )
    else
       kupper = 0.0
    end if

  end function kupper
  !! --------------------------------------------------------------------------------------------------------------------------- !!


  !! --------------------------------------------------------------------------------------------------------------------------- !!
  !>
  !! t-exp type source time function
  !!
  !! This function have >0 value among ts <= t <= infinity for given start time ts and characteristic time tr
  !! The amplitude is normalized by \int_0^\infty f(t) dt = 1
  !<
  real(SP) function texp( t, ts, tr )

    real(SP), intent(in) :: t  !< time
    real(SP), intent(in) :: ts !< rupture start time
    real(SP), intent(in) :: tr !< characteristic time
    !
    real(SP) :: tt
    !! ----

    if ( ts <= t ) then
       tt = t-ts
       texp = (2*PI)**2 * tt / (tr*tr) * exp( -2*PI*tt / tr )
    else
       texp = 0.0
    end if

  end function texp
  !! --------------------------------------------------------------------------------------------------------------------------- !!


  !! --------------------------------------------------------------------------------------------------------------------------- !!
  !>
  !! Cosine function for moment rate
  !!
  !! This function have >0 value among ts <= t <= ts + tr for given start time ts and rise time tr
  !! The amplitude is normalized by \int_0^\infty f(t) dt = 1
  !<
  real(SP) function cosine( t, ts, tr )

    real(SP), intent(in) :: t  !< time
    real(SP), intent(in) :: ts !< rupture start time
    real(SP), intent(in) :: tr !< rise time

    !! ----

    if ( ts <= t .and. t <= ts + tr ) then
       cosine = ( 1 - cos( 2*PI*(t-ts)/tr ) ) / tr
    else
       cosine = 0.0
    end if

  end function cosine
  !! --------------------------------------------------------------------------------------------------------------------------- !!


  !! --------------------------------------------------------------------------------------------------------------------------- !!
  !>
  !! Box-car function for moment rate
  !!
  !! This function have >0 value among ts <= t <= ts + tr for given start time ts and rise time tr
  !! The amplitude is normalized by \int_0^\infty f(t) dt = 1
  !<
  real(SP) function boxcar( t, ts, tr )

    real(SP), intent(in) :: t   !<  time
    real(SP), intent(in) :: ts  !<  rupture start time
    real(SP), intent(in) :: tr  !<  rise time

    !! ----

    if ( ts <= t .and. t <= ts + tr ) then
       boxcar = 1.0 / tr
    else
       boxcar = 0.0
    end if

  end function boxcar
  !! --------------------------------------------------------------------------------------------------------------------------- !!


  !! --------------------------------------------------------------------------------------------------------------------------- !!
  !>
  !! Triangle function for moment rate
  !!
  !! This function have >0 value among ts <= t <= ts + tr for given start time ts and rise time tr
  !! The amplitude is normalized by \int_0^\infty f(t) dt = 1
  !<
  real(SP) function triangle( t, ts, tr )

    real(SP), intent(in) :: t   !<  time
    real(SP), intent(in) :: ts  !<  rupture start time
    real(SP), intent(in) :: tr  !<  rise time

    !! ----

    if ( ts <= t .and. t <= ts + tr/2) then
       triangle = 4 * ( t-ts ) / ( tr*tr )
    else if ( ts + tr/2 < t .and. t <= ts + tr ) then
       triangle = - 4 * ( t - ts - tr ) / ( tr * tr )
    else
       triangle = 0.0
    end if

  end function triangle
  !! --------------------------------------------------------------------------------------------------------------------------- !!


  !! --------------------------------------------------------------------------------------------------------------------------- !!
  !>
  !! The Herrmann function for moment rate
  !!
  !! This function have >0 value among ts <= t <= ts + tr for given start time ts and rise time tr
  !! The amplitude is normalized by \int_0^\infty f(t) dt = 1
  !<
  real(SP) function herrmann ( t, ts, tr )

    real(SP), intent(in) :: t   !<  time
    real(SP), intent(in) :: ts  !<  rupture start time
    real(SP), intent(in) :: tr  !<  rise time

    real(SP) :: t1
    real(SP) :: t2

    !! ----

    t1 = ts +     tr / 4
    t2 = ts + 3 * tr / 4

    if( ts <= t .and. t < t1 ) then
       herrmann = 16 * ( t - ts )**2 / ( tr**3 )
    else if ( t1 <= t .and. t < t2 ) then
       herrmann = -2 * ( 8*( t*t + tr*ts + ts*ts- t*tr -2*t*ts ) + tr*tr ) / ( tr**3 )
    else if ( t2 <= t .and. t <= ts + tr ) then
       herrmann = 16 * ( ts + tr - t ) **2 / ( tr**3 )
    else
       herrmann = 0.0
    end if

  end function herrmann
  !! --------------------------------------------------------------------------------------------------------------------------- !!


  !! --------------------------------------------------------------------------------------------------------------------------- !!
  !>
  !! Time-integrated box-car function for moment function
  !!
  !! This function have >0 value among ts <= t <= ts + tr for given start time ts and rise time tr
  !! The amplitude will become iboxcar -> 1 as t->infinity
  !<
  real(SP) function iboxcar( t, ts, tr )

    real(SP), intent(in) :: t   !<  time
    real(SP), intent(in) :: ts  !<  rupture start time
    real(SP), intent(in) :: tr  !<  rise time

    !! ----

    if( t < ts ) then
       iboxcar = 0.0
    else if ( ts <= t .and. t <= ts + tr ) then
       iboxcar = (t-ts) / tr
    else
       iboxcar = 1.0
    end if

  end function iboxcar
  !! --------------------------------------------------------------------------------------------------------------------------- !!


  !! --------------------------------------------------------------------------------------------------------------------------- !!
  !>
  !! The time-integrated triangle function for moment function
  !!
  !! This function have >0 value among ts <= t <= ts + tr for given start time ts and rise time tr
  !! The amplitude will become itriangle -> 1 as t->infinity
  !<
  real(SP) function itriangle( t, ts, tr )

    real(SP), intent(in) :: t   !<  time
    real(SP), intent(in) :: ts  !<  rupture start time
    real(SP), intent(in) :: tr  !<  rise time

    !! ----

    if( t-ts < 0 ) then
       itriangle = 0
    else if ( 0 <= t-ts .and. t-ts <= tr/2) then
       itriangle = (2*(t - ts)**2)/tr**2
    else if ( ts + tr/2 < t .and. t <= ts + tr ) then
       itriangle = -1 + (4*(t - ts))/tr - (2*(t - ts)**2)/tr**2
    else
       itriangle = 1.0
    end if

  end function itriangle
  !! --------------------------------------------------------------------------------------------------------------------------- !!


  !! --------------------------------------------------------------------------------------------------------------------------- !!
  !>
  !! The time-integrated Herrmann function for moment function
  !!
  !! This function have >0 value among ts <= t <= ts + tr for given start time ts and rise time tr
  !! The amplitude will become iherrmann -> 1 as t->infinity
  !<
  real(SP) function iherrmann ( t, ts, tr )

    real(SP), intent(in) :: t   !<  time
    real(SP), intent(in) :: ts  !<  rupture start time
    real(SP), intent(in) :: tr  !<  rise time

    real(SP) :: t1
    real(SP) :: t2

    !! ----

    t1 =      tr / 4
    t2 = ts + 3 * tr / 4

    if( t-ts < 0 ) then
       iherrmann = 0
    else if( 0 <= t -ts .and. t -ts < tr / 4 ) then
       iherrmann = (16*(t - ts)**3)/(3*tr**3)
    else if ( tr/4 <= t-ts .and. t-ts < 3*tr/4 ) then
       iherrmann = (tr**3 - 12*tr**2*(t - ts) + 48*tr*(t - ts)**2 - 32*(t - ts)**3)/(6*tr**3)
    else if ( 3*tr/4 <= t-ts .and. t <= ts + tr ) then
       iherrmann = (-13*tr**3 + 48*tr**2*(t - ts) - 48*tr*(t - ts)**2 + 16*(t - ts)**3)/(3*tr**3)
    else
       iherrmann = 1.0
    end if

  end function iherrmann
  !! --------------------------------------------------------------------------------------------------------------------------- !!


  !! --------------------------------------------------------------------------------------------------------------------------- !!
  !>
  !! Time-integrated single-lobed Kupper function for moment function
  !!
  !! This function have >0 value among ts <= t <= ts + tr for given start time ts and rise time tr
  !! The amplitude will become ikupper -> 1 as t->infinity
  !<
  real(SP) function ikupper( t, ts, tr )

    real(SP), intent(in) :: t  !< time
    real(SP), intent(in) :: ts !< rupture start time
    real(SP), intent(in) :: tr !< rise time

    !! ----

    if ( t< ts ) then
       ikupper = 0.0
    else if ( ts <= t .and. t <= ts + tr ) then
       ikupper = (2 + cos(( PI * (t - ts)) / tr )) * sin(( PI *(t - ts) ) / (2*tr) )**4
    else
       ikupper = 1.0
    end if

  end function ikupper
  !! --------------------------------------------------------------------------------------------------------------------------- !!




  !! --------------------------------------------------------------------------------------------------------------------------- !!
  !>
  !! convert coordinate location to voxel number
  !<
  integer function x2i ( x, xbeg, dx )

    real(SP), intent(in) :: x    !< location
    real(SP), intent(in) :: xbeg !< lower boundary value of the coordinate
    real(SP), intent(in) :: dx   !< grid width

    !! ----

    x2i = ceiling( ( x - xbeg )  / dx )

  end function x2i
  !! --------------------------------------------------------------------------------------------------------------------------- !!


  !! --------------------------------------------------------------------------------------------------------------------------- !!
  !>
  !! convert coordinate location to voxel number
  !<
  integer function y2j ( y, ybeg, dy )

    real(SP), intent(in) :: y    !< location
    real(SP), intent(in) :: ybeg !< lower boundary value of the coordinate
    real(SP), intent(in) :: dy   !< grid width

    !! ----

    y2j = ceiling( ( y - ybeg ) / dy )

  end function y2j
  !! --------------------------------------------------------------------------------------------------------------------------- !!


  !! --------------------------------------------------------------------------------------------------------------------------- !!
  !>
  !! convert coordinate location to voxel number
  !<
  integer function z2k ( z, zbeg, dz )

    real(SP), intent(in) :: z    !< location
    real(SP), intent(in) :: zbeg !< lower boundary value of the coordinate
    real(SP), intent(in) :: dz   !< grid width

    !! ----

    z2k = ceiling( ( z - zbeg ) / dz )

  end function z2k
  !! --------------------------------------------------------------------------------------------------------------------------- !!

  !! --------------------------------------------------------------------------------------------------------------------------- !!
  !>
  !! convert voxel number to coordinate location
  !<
  real(SP) function i2x( i, xbeg, dx )

    integer,  intent(in) :: i    !< voxel number
    real(SP), intent(in) :: xbeg !< lower boundary value of the coordinate
    real(SP), intent(in) :: dx   !< grid width

    !! ----

    i2x = xbeg + ( i - 0.5 ) * dx

  end function i2x
  !! --------------------------------------------------------------------------------------------------------------------------- !!

  !! --------------------------------------------------------------------------------------------------------------------------- !!
  !>
  !! convert voxel number to coordinate location
  !<
  real(SP) function j2y( j, ybeg, dy )

    integer,  intent(in) :: j    !< voxel number
    real(SP), intent(in) :: ybeg !< lower boundary value of the coordinate
    real(SP), intent(in) :: dy   !< grid width

    !! ----

    j2y = ybeg + ( j - 0.5 ) * dy

  end function j2y
  !! --------------------------------------------------------------------------------------------------------------------------- !!

  !! --------------------------------------------------------------------------------------------------------------------------- !!
  !>
  !! convert voxel number to coordinate location
  !<
  real(SP) function k2z( k, zbeg, dz )

    integer,  intent(in) :: k    !< voxel number
    real(SP), intent(in) :: zbeg !< lower boundary value of the coordinate
    real(SP), intent(in) :: dz   !< grid width

    !! ----

    k2z = zbeg + ( k - 0.5 ) * dz

  end function k2z
  !! --------------------------------------------------------------------------------------------------------------------------- !!


  !! --------------------------------------------------------------------------------------------------------------------------- !!
  !>
  !! convert time grid number to physical time
  !<
  real(SP) function n2t( n, tbeg, dt )

    integer,  intent(in) :: n    !< time grid number
    real(SP), intent(in) :: tbeg !< start time
    real(SP), intent(in) :: dt   !< grid width

    !! ----

    n2t = tbeg + ( n - 0.5 ) * dt

  end function n2t
  !! --------------------------------------------------------------------------------------------------------------------------- !!


  !! --------------------------------------------------------------------------------------------------------------------------- !!
  !>
  !! Set relaxed time of viscoelastic medium based on the order of visco-elastic medium nm and frequency range (fmin, fmax)
  !!
  !!  - nm = 0 : Do nothing, it is elastic medium
  !!  - nm = 1 : Single viscoelastic medium: tau_sigma is given as inverse of central angular frequency (Hestholm, 1999)
  !!  - nm > 1 : Evenly-spaced tau_sigma in logarithmic space between frequency band fmin and fmax
  !<
  !!
  subroutine visco_set_relaxtime( nm, ts, fmin, fmax )

    integer,  intent(in)  :: nm
    real(SP), intent(in)  :: fmin, fmax
    real(SP), intent(out) :: ts(:)
    real(SP) :: omega_a, omega_b, omega

    integer :: im

    omega_a = 2 * Pi * fmin
    omega_b = 2 * Pi * fmax

    select case (nm)
    case( 0 )  ! elastic: do nothing

    case( 1 )  ! single Zener body: set center frequency

       omega = sqrt( omega_a * omega_b )
       ts(1) = 1.0 / omega

    case default ! generalized zener body

       ! evenly-spaced relaxed time in logarithmic scale from 1/omega_a to 1/omega_b
       do im = 1, nm
          omega = omega_a * ( omega_b / omega_a )**( dble(im-1)/dble(nm-1) )
          ts(im) = 1.0 / omega
       end do

    end select

  end subroutine visco_set_relaxtime
  !! --------------------------------------------------------------------------------------------------------------------------- !!


  !! --------------------------------------------------------------------------------------------------------------------------- !!
  !>
  !! Set chi-function for velocity dispersion as a function of frequency
  !<
  !! ----
  real(SP) function visco_chi( nm, ts, tau, f )

    integer, intent(in)  :: nm      !< order of viscoelastic model \f$ N_m \f$
    real(SP), intent(in) :: ts(nm)  !< relaxation times \f$ \tau_\sigma \f$
    real(SP), intent(in) :: tau     !< tau function of Blanch
    real(SP), intent(in) :: f       !< reference frequency of velocity

    integer  :: im
    real(SP) :: omega
    complex(SP) :: cc
    !! ----

    visco_chi = 0.0
    omega = 2 * Pi * f
    cc = 0.0
    do im=1, nm
       cc = cc + ( 1 - EI * omega * ts(im) * ( 1 + tau ) ) / ( 1 - EI * omega * ts(im) )
    end do
    cc = cc / nm

    visco_chi = 1 / real( cc**(-0.5) )

  end function visco_chi
  !! --------------------------------------------------------------------------------------------------------------------------- !!


  !! --------------------------------------------------------------------------------------------------------------------------- !!
  !>
  !! Set zeta-value from given relax time and order of the viscoelastic medium
  !!
  !!   zeta = \int_{omega_a}^{omega_b} F( omega ) d omega  / int_{omega_a}^{omega_b} F^2( omega ) d omega
  !!   F(omega) = \sum_{im=1}^{nm} omega tau_sigma(l) / ( 1 + omega^2 tau_sigma^2(l)
  !<
  !! ----
  real(SP) function visco_constq_zeta ( nm, fmin, fmax, ts  )

    !! Arguments
    integer, intent(in)   :: nm
    real(SP), intent(in)  :: fmin, fmax
    real(SP), intent(in)  :: ts(:)

    !! --

    real(SP) :: i1(nm), i2(nm,nm), i0(nm)
    real(SP) :: i1sum, i2sum, i0sum
    integer  :: im, km
    real(SP) :: wk1, wk2
    real(SP) :: om_a, om_b
    !! ----

    if( nm==0 ) then
       visco_constq_zeta = 0.0
    else

       !! angular frequency
       om_a = 2 * PI * fmin
       om_b = 2 * PI * fmax

       !! Integrants required for tau-method by Blanch
       do im=1, nm

          i0(im) = ( log( 1.0 + om_b**2 * ts(im)**2 )  -  log(1.0 + om_a**2 * ts(im)**2 ) )  /  ( 2 * ts(im) )

          i1(im) = ( ( atan( om_b*ts(im) ) - om_b*ts(im) / ( 1 + om_b**2 * ts(im)**2 ) ) &
               - ( atan( om_a*ts(im) ) - om_a*ts(im) / ( 1 + om_a**2 * ts(im)**2 ) ) ) / ( 2 * ts(im ) )

       end do

       do im=1, nm-1
          do km=im+1, nm
             wk1 = atan( om_b*ts(im) ) / ts(im) - atan( om_b*ts(km) ) / ts(km)
             wk2 = atan( om_a*ts(im) ) / ts(im) - atan( om_a*ts(km) ) / ts(km)
             i2(km,im) = ts(im) * ts(km) / ( ts(km)**2 - ts(im)**2 ) * ( wk1 - wk2 )
          end do
       end do

       i0sum = sum(i0(:))
       i1sum = sum(i1(:))
       i2sum = 0.0
       do im = 1, nm-1
          do km = im+1, nm
             i2sum = i2sum + i2(km,im)
          end do
       end do

       visco_constq_zeta = i0sum / ( i1sum + 2*i2sum )

    end if

  end function visco_constq_zeta
  !! --------------------------------------------------------------------------------------------------------------------------- !!

  !! --------------------------------------------------------------------------------------------------------------------------- !!
  !>
  !! Obtain list without duplicaiton from input list.
  !! Pointer integer array to the uniq list from the original list also will be returned.
  !<
  subroutine independent_list( nin, list, nind, ptr, list_ind )

    integer,      intent(in)  :: nin            !< number of input list data
    character(*), intent(in)  :: list(nin)      !< input list
    integer,      intent(out) :: nind           !< number of independent lines
    integer,      intent(out) :: ptr(nin)       !< pointer to the independent list
    character(*), intent(out) :: list_ind(nin)  !< independent components. 1..nind will be used.
    !! --
    integer :: i, j
    logical :: is_independent
    !! ----

    !! first data must be an independent component
    list_ind(1) = list(1)
    ptr (1) = 1
    nind    = 1

    do i=2, nin
       is_independent = .true. ! initialize flag

       !! check duplication
       do j=1, nind
          if( trim( list(i) ) == trim( list_ind(j) ) ) then
             is_independent   = .false.
             ptr(i) = j
             exit
          end if
       end do

       !! add independent component
       if( is_independent ) then
          nind   = nind + 1
          ptr(i) = nind
          list_ind(nind) = list(i)
       end if

    end do


  end subroutine independent_list
  !! --------------------------------------------------------------------------------------------------------------------------- !!


end module m_fdtool
!! ----------------------------------------------------------------------------------------------------------------------------- !!
