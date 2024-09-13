
#include "../shared/m_debug.h"
module m_pnm

    !! Read/Write pnm (color ppm / grayscale pgm ) files
    !!
    !! Copyright 2013-2024 Takuto Maeda. All rights reserved. This project is released under the MIT license.

    use m_std
    use m_debug
    implicit none
    private

    public :: ppm__write
    public :: ppm__read
    public :: pgm__write_8
    public :: pgm__write_16
    public :: pgm__read_8
    public :: pgm__read_16

contains

    subroutine ppm__write(fname, width, height, img)

        !! Create ppm-formatted color image with size of width x height pixels
        !!
        !! Color code must be given in 256 value (0-255) in the array img(:,:,:) as integer format
        !! This routine use Fortran2003 statements

        character(*), intent(in) :: fname
        integer, intent(in) :: width
        integer, intent(in) :: height
        integer, intent(in) :: img(3, width, height) !< RGB 0-255

        character :: aimg(3, width, height)
        integer :: io_ppm
        integer :: i, j, k
        integer :: ierr

        !! Convert one-byte ascii
        do i = 1, height
            do j = 1, width
                do k = 1, 3
                    aimg(k, j, i) = transfer(max(min(img(k, j, i), 255), 0), 'a')
                end do
            end do
        end do

        !! Ascii output of header part, then close file once
        open (newunit=io_ppm, file=fname, iostat=ierr)
        call assert(ierr == 0)
        write (io_ppm, '(A2)') 'P6'
        write (io_ppm, '(2I8)') width, height
        write (io_ppm, '(A3)') '255'
        close (io_ppm)

        !! Re-open the file in the stream i/o mode to add binary data
        open (io_ppm, file=fname, access='stream', position='append', form='unformatted', iostat=ierr)
        call assert(ierr == 0)
        write (io_ppm) aimg
        close (io_ppm)

    end subroutine ppm__write

    subroutine ppm__read(fname, width, height, image)

        !! Read color ppm-formatted image file and stores the color code as integer to array img

        character(*), intent(in)  :: fname
        integer, intent(in)  :: width
        integer, intent(in)  :: height
        integer, intent(out) :: image(3, width, height)

        character, parameter  :: lf = char(10) ! 改行コード
        integer       :: io_ppm
        character(80) :: buf
        character     :: aimage(3, width, height)
        integer       :: i, j, k
        integer       :: ierr
\
        !! Header part: first read long bufffe, then look for line break
        open (newunit=io_ppm, file=fname, access='stream', form='unformatted', iostat=ierr)
        call assert(ierr == 0)
        read (io_ppm) buf
        rewind (io_ppm)
        read (io_ppm) buf(1:index(buf, lf, back=.true.))
        read (io_ppm) aimage
        close (io_ppm)

        !! Convert to integer image
        do k = 1, height
            do j = 1, width
                do i = 1, 3
                    image(i, j, k) = transfer(aimage(i, j, k), width)
                end do
            end do
        end do

    end subroutine ppm__read

    subroutine pgm__write_16(fname, width, height, img)

        !! Create pgm-formatted grayscale image with size of width x height pixels
        !!
        !! Grayscale code must be given in 65536 value (0-65535) in the array img(:,:) as integer format
        
        character(*), intent(in) :: fname
        integer, intent(in) :: width
        integer, intent(in) :: height
        integer, intent(in) :: img(width, height) !< gray color code 0-65535

        character :: aimg(2, width, height)
        integer :: img2(2, width, height)
        integer :: i, j
        integer :: io_pgm
        integer :: ierr
        
        !! data scaling
        do j = 1, height
            do i = 1, width
                img2(1, i, j) = max(min(img(i, j) / 256, 255), 0)
                img2(2, i, j) = max(mod(img(i, j), 256), 0)
            end do
        end do

        !! Convert one-byte character
        do j = 1, height
            do i = 1, width
                aimg(1, i, j) = transfer(img2(1, i, j), 'a')
                aimg(2, i, j) = transfer(img2(2, i, j), 'a')
            end do
        end do

        !! Output header part in formatted text
        open (newunit=io_pgm, file=fname, iostat=ierr)
        call assert(ierr == 0)
        write (io_pgm, '(A2)') 'P5'
        write (io_pgm, '(2I8)') width, height
        write (io_pgm, '(A5)') '65535'

        !! Once close the file
        close (io_pgm)

        !! Re-open the same file in the stream access mode to add binary data
        open (io_pgm, file=fname, access='stream', position='append', form='unformatted', iostat=ierr)
        call assert(ierr == 0)
        write (io_pgm) aimg
        close (io_pgm)

    end subroutine pgm__write_16

    subroutine pgm__read_16(fname, width, height, image)

        !! Read grayscale pgm-formatted image file and stores the color code as integer to array img
        !!
        !! It assume that the data is stored in the Little-endian format

        character(*), intent(in)  :: fname
        integer, intent(in)  :: width
        integer, intent(in)  :: height
        integer, intent(out) :: image(width, height)

        character, parameter  :: lf = char(10) !< Line-break
        integer       :: io_pgm
        character(80) :: buf
        character     :: aimage(2, width, height)
        integer       :: i, j, k
        integer       :: var1, var2
        integer       :: ierr

        !! Read header part: First read long array, then look for line-break
        open (newunit=io_pgm, file=fname, access='stream', form='unformatted', iostat=ierr)
        call assert(ierr == 0)

        read (io_pgm) buf
        rewind (io_pgm)
        i = index(buf(:), lf)
        j = index(buf(i + 1:), lf) + i
        k = index(buf(j + 1:), lf) + j

        !! data portion
        read (io_pgm) aimage
        close (io_pgm)

        !! convert to integer
        do j = 1, height
            do i = 1, width
                var1 = transfer(aimage(1, i, j), width)
                var2 = transfer(aimage(2, i, j), width)
                image(i, j) = var1 * 256 + var2
            end do
        end do

    end subroutine pgm__read_16

    subroutine pgm__write_8(fname, width, height, img)

        !! Create pgm-formatted grayscale image with size of width x height pixels
        !!
        !! Grayscale code must be given in 256 value (0-255) in the array img(:,:) as integer format
        !! This routine use Fortran2003 statements

        character(*), intent(in) :: fname
        integer, intent(in) :: width
        integer, intent(in) :: height
        integer, intent(in) :: img(width, height) !< gray scale 0-255

        character :: aimg(width, height)
        integer :: io_pgm
        integer :: i, j
        integer :: ierr
    
        
        do j = 1, height
            do i = 1, width
                aimg(i, j) = transfer(max(min(img(i, j), 255), 0), 'a')
            end do
        end do

        open (newunit=io_pgm, file=fname, iostat=ierr)
        call assert(ierr == 0)

        write (io_pgm, '(A2)') 'P5'
        write (io_pgm, '(2I8)') width, height
        write (io_pgm, '(A3)') '255'
        close (io_pgm)

        open (io_pgm, file=fname, access='stream', position='append', form='unformatted', iostat=ierr)
        call assert(ierr == 0)
        write (io_pgm) aimg
        close (io_pgm)

    end subroutine pgm__write_8

    subroutine pgm__read_8(fname, width, height, image)

        !! Read grayscale pgm-formatted image file and stores the color code as integer to array img

        character(*), intent(in)  :: fname
        integer, intent(in)  :: width
        integer, intent(in)  :: height
        integer, intent(out) :: image(width, height)

        character, parameter  :: lf = char(10) !< Linebreak
        integer       :: io_pgm
        character(80) :: buf
        character     :: aimage(width, height)
        integer       :: i, j, k
        integer       :: ierr

        open (newunit=io_pgm, file=fname, access='stream', form='unformatted', iostat=ierr)
        call assert(ierr == 0)

        read (io_pgm) buf
        rewind (io_pgm)
        i = index(buf(:), lf)
        j = index(buf(i + 1:), lf) + i
        k = index(buf(j + 1:), lf) + j
        read (io_pgm) buf(1:k)
        read (io_pgm) aimage
        close (io_pgm)

        do j = 1, height
            do i = 1, width
                image(i, j) = transfer(aimage(i, j), width)
            end do
        end do

    end subroutine pgm__read_8

end module m_pnm
