#include "m_debug.h"
!! ----------------------------------------------------------------------------------------------------------------------------- !!
!>
!! Read random media volume
!!
!! @copyright
!!   Copyright 2013-2016 Takuto Maeda. All rights reserved. This project is released under the MIT license.
!<
!! --
module m_rdrmed

  use m_std
  use m_debug
#ifdef _NETCDF
  use netcdf
#endif
  implicit none

  public

contains


  !! --------------------------------------------------------------------------------------------------------------------------- !!
  subroutine rdrmed__2d( ib, ie, kb, ke, fn_rmed, vol )

    integer,      intent(in)  :: ib, ie
    integer,      intent(in)  :: kb, ke
    character(*), intent(in)  :: fn_rmed
    real(SP),     intent(out) :: vol(kb:ke, ib:ie)
    !! --
    real(SP), allocatable :: hh(:,:)
    integer :: i, k, ii, kk
    integer :: ncid, vid
    character(80) :: xn, zn, vn
    integer :: nxc, nzc !< netcdf volume size
    !! ----

#ifdef _NETCDF
    call debug( fn_rmed )
    call assert( nf90_open( fn_rmed, NF90_NOWRITE, ncid ) == NF90_NOERR )
    !! size
    call assert( nf90_inquire_dimension( ncid, 1, xn, nxc ) == NF90_NOERR )
    call assert( nf90_inquire_dimension( ncid, 2, zn, nzc ) == NF90_NOERR )
    call assert( nf90_inquire_variable ( ncid, 3, vn      ) == NF90_NOERR )
    call assert( nf90_inq_varid        ( ncid, vn,    vid ) == NF90_NOERR )
    call debug( vn )
    call debug( vid )
    call debug( ncid )
    allocate( hh(nxc,nzc) )
    call assert( nf90_get_var( ncid, vid, hh ) == NF90_NOERR )

    do k=kb, min(ke,nzc)

       if( k <= 0 ) then
          kk = k + nzc
       else
          kk = k
       end if

       do i=ib, ie
          ii = mod(i,nxc)
          if( ii<=0 ) ii = ii + nxc
          vol(k,i) = hh(ii,kk)
       end do
    end do

    deallocate( hh )

    !! bottom cyclic part
    do k=nzc+1, ke
       kk = mod(k,nzc)
       if( kk<=0 ) kk = kk + nzc
       vol(k,ib:ie) = vol(kk,ib:ie)
    end do


#endif

  end subroutine rdrmed__2d
  !! --------------------------------------------------------------------------------------------------------------------------- !!


  !! --------------------------------------------------------------------------------------------------------------------------- !!
  subroutine rdrmed__3d( ib, ie, jb, je, kb, ke, fn_rmed, vol )

    integer,      intent(in)  :: ib, ie
    integer,      intent(in)  :: jb, je
    integer,      intent(in)  :: kb, ke
    character(*), intent(in)  :: fn_rmed
    real(SP),     intent(out) :: vol(kb:ke, ib:ie, jb:je)
    !! --
    real(SP), allocatable :: hh(:,:)
    integer :: i, j, k, ii, jj, kk
    integer :: ncid, vid
    character(80) :: xn, yn, zn, vn
    integer :: nxc, nyc, nzc !< netcdf volume size
    integer :: st(3), ct(3)
    !! ----

#ifdef _NETCDF

    call assert( nf90_open( fn_rmed, NF90_NOWRITE, ncid ) == NF90_NOERR )
    !! size
    call assert( nf90_inquire_dimension( ncid, 1, xn, nxc ) == NF90_NOERR )
    call assert( nf90_inquire_dimension( ncid, 2, yn, nyc ) == NF90_NOERR )
    call assert( nf90_inquire_dimension( ncid, 3, zn, nzc ) == NF90_NOERR )
    call assert( nf90_inquire_variable ( ncid, 4, vn      ) == NF90_NOERR )
    call assert( nf90_inq_varid        ( ncid, vn,    vid ) == NF90_NOERR )

    allocate( hh(nxc,nyc) )
    st(1:3) = (/1,1,1/)
    ct(1:3) = (/nxc,nyc,1/)

    do k=kb, min(ke,nzc)

       if( k <= 0 ) then
          kk = k + nzc
       else
          kk = k
       end if

       st(3) = kk
       call assert( nf90_get_var(ncid, vid, hh, start=st, count=ct ) == NF90_NOERR )

       do j=jb, je

          jj = mod(j,nyc)
          if( jj<=0 ) jj = jj + nyc

          do i=ib, ie
             ii = mod(i,nxc)
             if( ii<=0 ) ii = ii + nxc
             vol(k,i,j) = hh(ii,jj)

          end do
       end do

    end do

    deallocate( hh )

    !! bottom cyclic part
    do k=nzc+1, ke
       kk = mod(k,nzc)
       vol(k,ib:ie,jb:je) = vol(kk,ib:ie,jb:je)
    end do

#endif

  end subroutine rdrmed__3d
  !! --------------------------------------------------------------------------------------------------------------------------- !!

end module m_rdrmed
!! ----------------------------------------------------------------------------------------------------------------------------- !!
