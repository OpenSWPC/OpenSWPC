!! ----------------------------------------------------------------------------------------------------------------------------- !!
!>
!! bitmap figure
!!
!! @copyright
!!   Copyright 2013-2016 Takuto Maeda. All rights reserved. This project is released under the MIT license.
!<
!! --
#include "m_debug.h"
module m_bmp

  use m_std
  use m_debug
  implicit none

  private
  public :: bmp__write

contains

  !! --------------------------------------------------------------------------------------------------------------------------- !!
  !>
  !! Create windows bitmap (bmp) formatted color image with size of width x height pixels
  !!
  !! Color code must be given in 256 value (0-255) in the array img(:,:,:) as integer format
  !<
  !! --
  subroutine bmp__write( fname, width, height, img )

    !! -- Arguments
    character(*), intent(in) :: fname
    integer,      intent(in) :: width
    integer,      intent(in) :: height
    integer,      intent(in) :: img( 3, width, height) !< RGB 0-255

    !! BITMAP FILE HEADER
    character(2) :: bfType ='BM' !< Filetype 'BM'
    integer(4)   :: bfSize       !< file size (byte)
    integer(2)   :: bfReserve = 0
    integer(4)   :: bfOffBits = 54

    !! BITMAP INFORMATION HEADER
    integer(4) :: biSize = 40
    integer(2) :: biPlanes = 1
    integer(2) :: biBitCount = 24    !< Color bit number
    integer(4) :: biCompression = 0  !< 0 means no compression
    integer(4) :: biSizeImage
    integer(4) :: biXPixPerMeter = 0 !< Resolusion, not used in this code
    integer(4) :: biYPixPerMeter = 0 !< Resolusion, not used in this code
    integer(4) :: biClrUsed = 0      !< Used color, not used in this code
    integer(4) :: biClrImportant = 0 !< Number of important colors, not used in this code
    integer    :: imgwid
    character(1), allocatable :: aimg(:,:)
    integer    :: io
    integer    :: ierr
    integer    :: i, j

    !! mod( image_width, 4) must be 0: column width including padding
    imgwid = 4 * ceiling( 3 * width / 4. )

    allocate( aimg(imgwid,height) )

    !! bmp image file
    call std__getio( io )
    open(io, file=trim(fname), access='stream', action='write', iostat=ierr, status='unknown', form='unformatted')
    call assert( ierr==0 )

    bfSize = 14 &           ! size of file header
           + 40 &           ! size of windows bitmap information header
           + imgwid*height  ! size of image data

    biSizeImage = 3*width*height

    !! header
    write(io) bfType
    write(io) bfSize
    write(io) bfReserve
    write(io) bfReserve
    write(io) bfOffBits
    write(io) biSize
    write(io) width
    write(io) height
    write(io) biPlanes
    write(io) biBitCount
    write(io) biCompression
    write(io) biSizeImage
    write(io) biXPixPerMeter
    write(io) biYPixPerMeter
    write(io) biClrUsed
    write(io) biClrImportant

    !! image data
    do j=1, height
       do i=1, width
          aimg(3*(i-1)+1,height-j+1) = transfer(img(3,i,j),'a')
          aimg(3*(i-1)+2,height-j+1) = transfer(img(2,i,j),'a')
          aimg(3*i      ,height-j+1) = transfer(img(1,i,j),'a')
       end do
    end do
    aimg(3*width+1:imgwid,height-j+1) = transfer(0,'a')
    write(io) aimg

    close(io)
    deallocate(aimg)

  end subroutine bmp__write
  !! --------------------------------------------------------------------------------------------------------------------------- !!

end module m_bmp
!! ----------------------------------------------------------------------------------------------------------------------------- !!
