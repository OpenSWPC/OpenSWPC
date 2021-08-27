
!! ----------------------------------------------------------------------------------------------------------------------------- !!
!>
!! Convert free surface/ocean bottom snapshot data from OpenSWPC to GMT-compatible lon-lat-var grid dataset
!!
!! @Usage
!! fs2grd.x -i swpc-snapshot.nc -R lon0/lon1/lat0/lat1 -dlon delta_lon -dlat delta_lat -v varname
!!
!! @copyright
!! Copyright 2020-2021 Takuto Maeda. All rights reserved. This project is released under the MIT license.
!<
!! --
program fs2grd

  use m_std  
  use m_getopt
  use m_bicubic
  use m_fdsnap
  use m_geomap
  use netcdf
  implicit none

  character(80) :: title
  character(80) :: varname
  integer :: ncid_in
  integer :: nx, ny, nt
  real, allocatable :: lon_in(:,:), lat_in(:,:), Var_in(:,:)
  real, allocatable :: lon_g(:), lat_g(:), xg(:,:), yg(:,:)
  real :: dt
  integer :: nlon, nlat
  real, allocatable :: x_in(:), y_in(:)
  real :: dx, dy
  character(256) :: fn_in
  integer :: var_id
  real :: dlon, dlat
  logical :: is_time_dependent_variable

  call option_processings(nlon, nlat, dlon, dlat, lon_g, lat_g, fn_in, varname)

  is_time_dependent_variable = &
    ( trim(varname) == 'Vx' .or. trim(varname) == 'Vy' .or. trim(varname) == 'Vz' .or. &
      trim(varname) == 'Ux' .or. trim(varname) == 'Uy' .or. trim(varname) == 'Uz' .or. &
      trim(varname) == 'div' .or. trim(varname) == 'rot_x' .or. trim(varname) == 'rot_y' .or. &
      trim(varname) == 'rot_z' )

 
  call open_input_netcdf(ncid_in, var_id, x_in, y_in, xg, yg, title)

  if( is_time_dependent_variable ) then
    call read_export_timedepend
  else
    call read_export_static
  end if
  

contains

  subroutine read_export_timedepend() 

    integer :: it, i, j
    integer :: cnt(3), stt(3)
    real, allocatable :: Var_g(:,:)
    type(bicubic__data) :: bd
    real:: nan
    !--
    nan = -1.0
    nan = sqrt(nan)

    allocate(Var_g(nlon, nlat))

    cnt = (/ nx, ny, 1/)

    do it=1, nt
      stt = (/ 1, 1, it /)
      call nc_chk(nf90_get_var(ncid_in, var_id, var_in, start=stt, count=cnt))
    
      ! interpolate
      call bicubic__init(bd, nx, ny, x_in(1), y_in(1), dx, dy, var_in)
      do j=1, nlat
        do i=1, nlon
          call bicubic__interp( bd, xg(i,j), yg(i,j), var_g(i,j), nan)
        end do
      end do
      call bicubic__terminate(bd)

      call dump(it, Var_g)

    end do

  end subroutine read_export_timedepend


  subroutine read_export_static()

    integer :: i,j 
    real, allocatable :: Var_g(:,:)
    type(bicubic__data) :: bd
    real:: nan

    nan = -1.0
    nan = sqrt(nan)

    allocate(Var_g(nlon, nlat))


   call nc_chk(nf90_get_var(ncid_in, var_id, var_in))
    
    ! interpolate
    call bicubic__init(bd, nx, ny, x_in(1), y_in(1), dx, dy, var_in)
    do j=1, nlat
      do i=1, nlon
        call bicubic__interp( bd, xg(i,j), yg(i,j), var_g(i,j), nan)
      end do
    end do
    call bicubic__terminate(bd)

    call dump(-1, Var_g)

  end subroutine read_export_static

  subroutine open_input_netcdf(ncid_in, var_id, x_in, y_in, xg, yg, title)

    integer, intent(out) :: ncid_in, var_id
    real, allocatable :: x_in(:), y_in(:)
    real, allocatable, intent(out) :: xg(:,:), yg(:,:)
    character(*), intent(out) :: title
    character(80) :: vname
    integer :: i, j
    real :: clon, clat, phi
    integer :: latid, lonid
    real :: xx, yy
    !----


    call nc_chk(nf90_open(fn_in, nf90_nowrite, ncid_in))

    call nc_chk(nf90_inquire_dimension(ncid_in, 1, vname, nx))
    call nc_chk(nf90_inquire_dimension(ncid_in, 2, vname, ny))
    call nc_chk(nf90_inquire_dimension(ncid_in, 3, vname, nt))
    call nc_chk(nf90_get_att(ncid_in, nf90_global, 'title', title))
    call nc_chk(nf90_get_att(ncid_in, nf90_global, 'dt', dt))
    call nc_chk(nf90_get_att(ncid_in, nf90_global, 'ds1', dx))
    call nc_chk(nf90_get_att(ncid_in, nf90_global, 'ds2', dy))
    call nc_chk(nf90_get_att(ncid_in, nf90_global, 'clon', clon))
    call nc_chk(nf90_get_att(ncid_in, nf90_global, 'clat', clat))
    call nc_chk(nf90_get_att(ncid_in, nf90_global, 'phi', phi))

    allocate(lon_in(nx,ny), lat_in(nx,ny), Var_in(nx,ny))
    allocate(x_in(nx), y_in(ny))
    call nc_chk(nf90_inq_varid(ncid_in, 'lon', lonid))
    call nc_chk(nf90_inq_varid(ncid_in, 'lat', latid))
    call nc_chk(nf90_get_var(ncid_in, 1, x_in))
    call nc_chk(nf90_get_var(ncid_in, 2, y_in))

    call nc_chk(nf90_get_var(ncid_in, lonid, lon_in))
    call nc_chk(nf90_get_var(ncid_in, latid, lat_in))

    !! grid locations in the cartesian coordinate
    allocate(xg(nlon,nlat), yg(nlon,nlat))

    do j=1, nlat
      do i=1, nlon
        call geomap__g2c(lon_g(i), lat_g(j), clon, clat, phi, xg(i,j), yg(i,j) )
      end do
    end do

   call nc_chk(nf90_inq_varid(ncid_in, trim(varname), var_id))

    
  end subroutine open_input_netcdf

  subroutine option_processings(nlon, nlat, dlon, dlat, lon_g, lat_g, fn_in, varname)
    
    integer, intent(out) :: nlon, nlat
    real,    intent(out) :: dlon, dlat
    real, allocatable, intent(out) :: lon_g(:), lat_g(:)
    character(*), intent(out) :: fn_in
    character(*), intent(out) :: varname
    !--

    character(80) :: cregion
    character(80), allocatable :: reg(:)
    integer :: nr
    logical :: is_opt
    real :: lon_beg, lon_end, lat_beg, lat_end
    integer :: i, j


    call getopt('i', is_opt, fn_in, ''  )
    if( .not. is_opt ) then 
      write(stderr,*) "no input file given [-i]"
      stop
    end if
        
    call getopt('R', is_opt, cregion, '')
    if( .not. is_opt ) then
      write(stderr,*)  "no region is specified [-R]"
      stop
    end if
    call split(cregion, '/', nr, reg)

    call getopt('dlon', is_opt, dlon, 0.0)
    if( .not. is_opt ) then
      write(stderr,*)  "[-dlon] must be specified"
      stop
    end if
    
    call getopt('dlat', is_opt, dlat, 0.0)
    if( .not. is_opt ) then
      write(stderr,*)  "[-dlat] must be specified"
      stop
    end if

    call getopt('v', is_opt, varname, '')
    if( .not. is_opt ) then
      write(stderr,*) "[-v] must be specified"
      stop
    end if

    if( nr /= 4 ) then
      write(stderr,*) "invalid region specification"
      stop
    end if
    read(reg(1),*) lon_beg
    read(reg(2),*) lon_end
    read(reg(3),*) lat_beg
    read(reg(4),*) lat_end
   
    nlat = floor( ( lat_end - lat_beg ) / dlat + 0.5 ) + 1
    nlon = floor( ( lon_end - lon_beg ) / dlon + 0.5 ) + 1
    allocate(lon_g(nlon), lat_g(nlat))
    do i=1, nlon
      lon_g(i) = lon_beg + (i-1) * dlon
    end do
    do j=1, nlat
      lat_g(j) = lat_beg + (j-1) * dlat
    end do

  end subroutine option_processings
  

  subroutine dump(it, Var)
    
    integer, intent(in) :: it
    real, intent(in) :: Var(nlon, nlat)
    !--
    character(6) :: cit
    character(256) :: fn_out
    integer :: ncid
    integer :: dimid_lon, dimid_lat, varid_lon, varid_lat, varid
    
    if( it>= 0 ) then  
      write(cit,'(I6.6)') it
      fn_out = trim(adjustl(title)) // '.' // trim(varname) // '.' // cit // '.nc' 
    else
      fn_out = trim(adjustl(title)) // '.' // trim(varname) // '.nc' 
    end if

    write(*,*) trim(fn_out)

    call nc_chk(nf90_create(fn_out, nf90_clobber, ncid))
    call nc_chk(nf90_def_dim(ncid, 'lon', nlon, dimid_lon))
    call nc_chk(nf90_def_dim(ncid, 'lat', nlat, dimid_lat))
    call nc_chk(nf90_def_var(ncid, 'lon', nf90_real, dimid_lon, varid_lon))
    call nc_chk(nf90_def_var(ncid, 'lat', nf90_real, dimid_lat, varid_lat))
    call nc_chk(nf90_put_att(ncid, nf90_global, 'title', trim(title)))
    call nc_chk(nf90_put_att(ncid, varid_lon, 'long_name', 'Longitude'))
    call nc_chk(nf90_put_att(ncid, varid_lat, 'long_name', 'Latitude'))
    call nc_chk(nf90_put_att(ncid, varid_lon, 'units', 'degrees_east'))
    call nc_chk(nf90_put_att(ncid, varid_lat, 'units', 'degrees_north'))
    call nc_chk(nf90_put_att(ncid, varid_lon, 'actual_range', (/lon_g(1), lon_g(nlon)/)))
    call nc_chk(nf90_put_att(ncid, varid_lat, 'actual_range', (/lat_g(1), lat_g(nlat)/)))
    
    call nc_chk(nf90_def_var(ncid, trim(varname), nf90_real, (/dimid_lon, dimid_lat/), varid))
    call nc_chk(nf90_put_att(ncid, varid, 'long_name', trim(varname)))
    call nc_chk(nf90_put_att(ncid, varid, 'coordinates', "lat lon"))

    call nc_chk(nf90_enddef(ncid))

    call nc_chk(nf90_put_var(ncid, varid_lon,  lon_g))
    call nc_chk(nf90_put_var(ncid, varid_lat,  lat_g))
    call nc_chk(nf90_put_var(ncid, varid, Var))
    call nc_chk(nf90_close(ncid))

  end subroutine dump

  subroutine nc_chk( ierr )

    integer, intent(in) :: ierr

    if( ierr /= nf90_noerr )  write(stderr,*) nf90_strerror( ierr )

  end subroutine nc_chk  

  subroutine split(buf, sep, n, nm)
    
    character(*), intent(in)  :: buf  !< input buffer
    character(*), intent(in)  :: sep  !< separator
    integer,      intent(out) :: n    !< number of separated fields
    !--
    character(*), allocatable, intent(out) :: nm(:)  !< separated variables
    integer :: i, j, k
    !----
    
    ! count separator
    n = 0
    i = 1
    j = 1
    do
      i = index(buf(j:), sep)
      if( i == 0 ) exit
      j = j + i
      n = n + 1
    end do
    n = n + 1  !! #words = #separator + 1
    allocate(nm(n))

    i = 1
    j = 1
    do k=1, n-1
      i = index(buf(j:), sep)
      nm(k) = trim(adjustl(buf(j:i+j-2)))
      j = j + i
    end do
    nm(n) = trim(adjustl(buf(j:)))

  end subroutine split


end program fs2grd