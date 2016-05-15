!! ----------------------------------------------------------------------------------------------------------------------------- !!
!>
!! Generate random velocity flucuation model in 3D space
!!
!! @copyright
!!   Copyright 2013-2016 Takuto Maeda. All rights reserved. This project is released under the MIT license.
!<
!! --
#include "m_debug.h"
program gen_rmed3d

  use m_std
  use m_getopt
  use m_rmedia
  use m_debug
  use netcdf

  implicit none

  real(SP)              :: ax, ay, az, epsil, kappa
  integer               :: nx, ny, nz
  integer               :: nnx, nny, nnz ! power of 2
  integer               :: px, py, pz
  logical               :: is_opt
  real(SP)              :: dx, dy, dz
  integer               :: ptype
  real(SP), allocatable :: med(:,:,:)
  real(SP), allocatable :: xx(:), yy(:), zz(:)
  integer               :: seed
  logical               :: is_seed
  character(256)        :: fn_out
  integer               :: ncid, vid
  integer :: i, j, k

  !! ----

  !! option processing
  call getopt('nx',    is_opt, nx);    if(.not. is_opt) call usage_exit
  call getopt('ny',    is_opt, ny);    if(.not. is_opt) call usage_exit
  call getopt('nz',    is_opt, nz);    if(.not. is_opt) call usage_exit
  call getopt('dx',    is_opt, dx);    if(.not. is_opt) call usage_exit
  call getopt('dy',    is_opt, dy);    if(.not. is_opt) call usage_exit
  call getopt('dz',    is_opt, dz);    if(.not. is_opt) call usage_exit
  call getopt('ax',    is_opt, ax);    if(.not. is_opt) call usage_exit
  call getopt('ay',    is_opt, ay);    if(.not. is_opt) call usage_exit
  call getopt('az',    is_opt, az);    if(.not. is_opt) call usage_exit
  call getopt('epsil', is_opt, epsil); if(.not. is_opt) call usage_exit
  call getopt('ptype', is_opt, ptype); if(.not. is_opt) call usage_exit
  kappa = 0
  call assert(1<=ptype .and. ptype<=3)

  if(ptype == 3) then
    call getopt('kappa', is_opt, kappa); if(.not. is_opt) call usage_exit
  end if
  call getopt('seed', is_seed, seed, -12345)
  call getopt('o', is_opt, fn_out)

  ! 2^(p-1) < nx <= 2^p
  px = ceiling(log(dble(nx)) / log(2.0_DP))
  py = ceiling(log(dble(ny)) / log(2.0_DP))
  pz = ceiling(log(dble(nz)) / log(2.0_DP))
  nnx = 2**px
  nny = 2**py
  nnz = 2**pz

  if(nnx - nx > 0 .or. nny - ny > 0 .or. nnz - nz > 0) then
    call info ('size is not power of 2')
  end if


  allocate(med(1:nnx, 1:nny, 1:nnz), xx(nnx), yy(nny), zz(nnz))

  !! generate random media (x<->z exchanged for SWPC's array configuration)
  if(is_seed) then
    call rmedia__3dgen_ani(ptype, nnx, nny, nnz, dx, dy, dz, ax, ay, az, epsil, kappa, xx, yy, zz, med, seed)
  else
    call rmedia__3dgen_ani(ptype, nnx, nny, nnz, dx, dy, dz, ax, ay, az, epsil, kappa, xx, yy, zz, med)
  end if

  !! data export
  call export_netcdf()

  deallocate(xx, yy, zz, med)
  !! --------------------------------------------------------------------------------------------------------------------------- !!

contains

  subroutine export_netcdf()

    integer :: ncid
    integer :: dimid_x, dimid_y, dimid_z
    integer :: varid_x, varid_y, varid_z, varid_r
    integer :: dimids(2)

    call nc_chk(nf90_create(fn_out, NF90_CLOBBER, ncid))
    call nc_chk(nf90_def_dim(ncid, 'x', nx, dimid_x))
    call nc_chk(nf90_def_dim(ncid, 'y', ny, dimid_y))
    call nc_chk(nf90_def_dim(ncid, 'z', nz, dimid_z))
    call nc_chk(nf90_def_var(ncid, 'x', NF90_REAL, dimid_x, varid_x))
    call nc_chk(nf90_def_var(ncid, 'y', NF90_REAL, dimid_y, varid_y))
    call nc_chk(nf90_def_var(ncid, 'z', NF90_REAL, dimid_z, varid_z))
    call nc_chk(nf90_put_att(ncid, NF90_GLOBAL, 'title', 'random media'))
    call nc_chk(nf90_put_att(ncid, NF90_GLOBAL, 'ax',ax))
    call nc_chk(nf90_put_att(ncid, NF90_GLOBAL, 'ay',ay))
    call nc_chk(nf90_put_att(ncid, NF90_GLOBAL, 'az',az))
    call nc_chk(nf90_put_att(ncid, NF90_GLOBAL, 'epsilon', epsil))
    select case (ptype)
    case(1)
      call nc_chk(nf90_put_att(ncid, NF90_GLOBAL, 'PSDF type', 'Gaussian'))
    case(2)
      call nc_chk(nf90_put_att(ncid, NF90_GLOBAL, 'PSDF type', 'Exponential'))
    case(3)
      call nc_chk(nf90_put_att(ncid, NF90_GLOBAL, 'PSDF type', 'von Karman'))
      call nc_chk(nf90_put_att(ncid, NF90_GLOBAL, 'kappa',kappa))
    end select

    call nc_chk(nf90_put_att(ncid, NF90_GLOBAL, 'dx', dx))
    call nc_chk(nf90_put_att(ncid, NF90_GLOBAL, 'dy', dy))
    call nc_chk(nf90_put_att(ncid, NF90_GLOBAL, 'dz', dz))

    call nc_chk(nf90_put_att(ncid, varid_x, 'long_name', 'x'))
    call nc_chk(nf90_put_att(ncid, varid_y, 'long_name', 'y'))
    call nc_chk(nf90_put_att(ncid, varid_z, 'long_name', 'z'))
    call nc_chk(nf90_put_att(ncid, varid_x, 'units', 'km'))
    call nc_chk(nf90_put_att(ncid, varid_y, 'units', 'km'))
    call nc_chk(nf90_put_att(ncid, varid_z, 'units', 'km'))
    call nc_chk(nf90_put_att(ncid, varid_x, 'actual_range', (/xx(1), xx(nx)/)))
    call nc_chk(nf90_put_att(ncid, varid_y, 'actual_range', (/yy(1), yy(nx)/)))
    call nc_chk(nf90_put_att(ncid, varid_z, 'actual_range', (/zz(1), zz(nz)/)))
    call nc_chk(nf90_def_var(ncid, "random media", NF90_REAL, (/dimid_x, dimid_y, dimid_z/),  varid_r))
    call nc_chk(nf90_put_att(ncid, varid_r, 'long_name', 'random media'))
    call nc_chk(nf90_put_att(ncid, varid_r, 'units', ''))

    call nc_chk(nf90_enddef(ncid))

    call nc_chk(nf90_put_var(ncid, varid_x, xx(1:nx)))
    call nc_chk(nf90_put_var(ncid, varid_y, yy(1:ny)))
    call nc_chk(nf90_put_var(ncid, varid_z, zz(1:nz)))
    call nc_chk(nf90_put_var(ncid, varid_r, med(1:nx,1:ny,1:nz)))
    call nc_chk(nf90_close(ncid))




  end subroutine export_netcdf
  !! --------------------------------------------------------------------------------------------------------------------------- !!

  !! --------------------------------------------------------------------------------------------------------------------------- !!
  subroutine usage_exit()
    write(STDERR,'(A)') 'gen_rmed3d.x  [-o outfile] [-nx nx] [-ny ny] [-nz nz] [-epsil epsilon] [-kappa kappa] '
    write(STDERR,'(A)') '              [-dx dx] [-dy dy] [-dz dz] [-ax ax] [-ay ay] [-az az]'
    write(STDERR,'(A)') '              [-ptype ptype] (1=Gaussian, 2=Exponential, 3=von Karman]) {-seed seed_number}'
    stop
  end subroutine usage_exit
  !! --------------------------------------------------------------------------------------------------------------------------- !!

  !! --------------------------------------------------------------------------------------------------------------------------- !!
  subroutine nc_chk(ierr)

    integer, intent(in) :: ierr

    if(ierr /= NF90_NOERR) then
      write(STDERR,*) NF90_STRERROR(ierr)
      stop
    end if

  end subroutine nc_chk
  !! --------------------------------------------------------------------------------------------------------------------------- !!

end program gen_rmed3d
!! ----------------------------------------------------------------------------------------------------------------------------- !!
