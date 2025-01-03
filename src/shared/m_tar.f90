module m_tar

    !! Tar-archive module
    !!
    !! Copyright 2025 Takuto Maeda, All rights reserved. This project is released under the MIT license. 

    use iso_c_binding, only: c_null_char

    implicit none
    public

    type tar__hdr
        character(100) :: fname
        integer        :: mode
        integer        :: uid
        integer        :: gid
        integer        :: fsize
        integer        :: mtime
        integer        :: checksum
        integer        :: typeflag
        character(100) :: linkname
        character(6)   :: magic
        character(2)   :: version
        character(32)  :: uname
        character(32)  :: gname
        character(8)   :: devmajor
        character(8)   :: devminor
        character(155) :: prefix
    end type tar__hdr

    character(512), private :: null
    logical, save :: initialized = .false.

contains

    subroutine tar__init

        integer :: i
        do i=1, 512
            null(i:i) = c_null_char
        end do

        initialized = .true.

    end subroutine tar__init
    

    subroutine tar__inithdr(th)

        !! Initialize the tar header information.

        type(tar__hdr), intent(out) :: th

        if(.not. initialized) call tar__init()

        th%fname    = c_null_char
        th%mode     = 420
        th%uid      = 0
        th%gid      = 0
        th%fsize    = 0
        th%mtime    = 1735657200
        th%checksum = 0
        th%typeflag = 0
        th%linkname = c_null_char
        th%magic    = 'ustar'
        th%version  = '00'
        th%uname    = 'root'
        th%gname    = 'root'
        th%devmajor = c_null_char
        th%devminor = c_null_char
        th%prefix   = c_null_char

    end subroutine tar__inithdr

    subroutine tar__printhdr(th)

        !! Confirm the header information of the tar file.

        type(tar__hdr), intent(in) :: th

        write (*, '(A, A)') 'fname:     ', trim(th%fname)
        write (*, '(A, I0)') 'mode:      ', th%mode
        write (*, '(A, I0)') 'uid:       ', th%uid
        write (*, '(A, I0)') 'gid:       ', th%gid
        write (*, '(A, I0)') 'fsize:     ', th%fsize
        write (*, '(A, I0)') 'mtime:     ', th%mtime
        write (*, '(A, I0)') 'checksum:  ', th%checksum
        write (*, '(A, I0)') 'typeflag:  ', th%typeflag
        write (*, '(A, A)') 'linkname:  ', trim(th%linkname)
        write (*, '(A, A)') 'magic:     ', trim(th%magic)
        write (*, '(A, A)') 'version:   ', trim(th%version)
        write (*, '(A, A)') 'uname:     ', trim(th%uname)
        write (*, '(A, A)') 'gname:     ', trim(th%gname)
        write (*, '(A, A)') 'devmajor:  ', th%devmajor
        write (*, '(A, A)') 'devminor:  ', th%devminor
        write (*, '(A, A)') 'prefix:    ', trim(th%prefix)

    end subroutine tar__printhdr


    subroutine tar__rfile(io, th, buf)

        integer, intent(in) :: io
        type(tar__hdr), intent(in) :: th
        character(:), allocatable, intent(out) :: buf
        character(512) :: pbuf
        integer :: i, n

        if (.not. allocated(buf)) then
            allocate (character(len=th%fsize) :: buf)
        end if
        n = ceiling(th%fsize/512.0)
        do i = 1, n
            read (io) pbuf
            if (i < n) then
                buf((i - 1)*512 + 1:i*512) = pbuf
            else
                buf((i - 1)*512 + 1:th%fsize) = pbuf(1:th%fsize - (i - 1)*512)
            end if
        end do

    end subroutine tar__rfile


    subroutine tar__rhdr(io, th, end_of_file)

        integer, intent(in) :: io
        type(tar__hdr), intent(out) :: th
        logical, intent(out) :: end_of_file
        integer :: i
        character(512) :: chdr

        read (io) chdr(1:512)
        end_of_file = .true.
        do i = 1, 512
            if (chdr(i:i) /= c_null_char) then
                end_of_file = .false.
            end if
        end do
        if (end_of_file) return

        call cread(chdr(1:100), th%fname)
        read (chdr(101:108), '(O8.8)') th%mode
        read (chdr(109:116), '(O8.8)') th%uid
        read (chdr(117:124), '(O8.8)') th%gid
        read (chdr(125:136), '(O12.12)') th%fsize
        read (chdr(137:148), '(O12.12)') th%mtime
        read (chdr(149:156), '(O8.8)') th%checksum
        read (chdr(157:157), '(I1)') th%typeflag

        call cread(chdr(158:257), th%linkname)
        call cread(chdr(258:265), th%magic)
        if (th%magic /= 'ustar') error stop "Not a tar archive: "

        th%version = chdr(264:265)
        call cread(chdr(266:297), th%uname)
        call cread(chdr(298:329), th%gname)
        call cread(chdr(330:337), th%devmajor)
        call cread(chdr(338:345), th%devminor)
        call cread(chdr(346:500), th%prefix)

    end subroutine tar__rhdr


    subroutine tar__whdr(io, th)

        !! write tar header information to a character block.

        integer, intent(in) :: io 
        type(tar__hdr), intent(inout) :: th

        character(512):: chdr
        integer :: i

        chdr = c_null_char
        chdr(1:1+len_trim(th%fname)) = th%fname // c_null_char
        write (chdr(101:108), '(O8.8)') th%mode
        write (chdr(109:116), '(O8.8)') th%uid
        write (chdr(117:124), '(O8.8)') th%gid
        write (chdr(125:136), '(O12.12)') th%fsize
        write (chdr(137:148), '(O12.12)') th%mtime
        chdr(149:156) = repeat(' ', 8)
        write (chdr(157:157), '(I1)') th%typeflag
        chdr(158:158 + len_trim(th%linkname)) = trim(th%linkname) // c_null_char
        chdr(258:258 + len_trim(th%magic)   ) = trim(th%magic)    // c_null_char
        chdr(266:266 + len_trim(th%uname)   ) = trim(th%uname)    // c_null_char
        chdr(298:298 + len_trim(th%gname)   ) = trim(th%gname)    // c_null_char
        chdr(330:330 + len_trim(th%devmajor)) = trim(th%devmajor) // c_null_char
        chdr(338:338 + len_trim(th%devminor)) = trim(th%devminor) // c_null_char
        chdr(346:346 + len_trim(th%prefix)  ) = trim(th%prefix)   // c_null_char

        !! checksum
        th%checksum = 0
        do i = 1, 512
            th%checksum = th%checksum + ichar(chdr(i:i))
        end do
        write (chdr(149:156), '(O8.8)') th%checksum

        write(io) chdr

    end subroutine tar__whdr


    subroutine tar__wbuf(io, buf)

        integer, intent(in) :: io
        character(*), intent(in) :: buf

        write(io) buf(:)
        call tar__wpad(io, len(buf))

    end subroutine tar__wbuf

    subroutine tar__wpad(io, fsize)

        !! padding zero to the end of the 512-byte block
        
        integer, intent(in) :: io
        integer, intent(in) :: fsize
        integer :: m

        m = mod(fsize, 512)
        write(io) null(m+1:512)

    end subroutine tar__wpad

    subroutine tar__wend(io)

        integer, intent(in) :: io

        write(io) null
        write(io) null

    end subroutine tar__wend


    subroutine cread(buf, c)

        character(*), intent(in) :: buf
        character(*), intent(out) :: c
        integer :: i

        do i = 1, len(buf)
            if (buf(i:i) == c_null_char) then
                c = buf(1:i-1)
                exit
            end if
        end do

    end subroutine cread


    subroutine tar__extract_filename(fullpath, filename)
        implicit none
        character(*), intent(in) :: fullpath
        character(*), intent(out) :: filename
        integer :: last_slash
        integer :: i

        last_slash = 0

        last_slash = max(index(fullpath, '/'), index(fullpath, '\'))

        if (last_slash > 0) then
            filename = fullpath(last_slash + 1:)
        else
            filename = fullpath
        end if

        do i = len_trim(filename) + 1, len(filename)
            filename(i:i) = c_null_char
        end do

    end subroutine tar__extract_filename

end module m_tar
