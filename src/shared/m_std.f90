module m_std

    !! Definition of standard constants, in/out io numbers, precision constants
    !!
    !! Copyright 2013-2025 Takuto Maeda. All rights reserved. This project is released under the MIT license.

    use iso_fortran_env, only: real64, real32
    implicit none
    private

    integer,     parameter, public :: DP = real64  !! Double Precision
    integer,     parameter, public :: SP = real32  !! Single Precision

    real(DP),    parameter, public :: PI        = atan(1.0_DP) * 4
    real(DP),    parameter, public :: R_EARTH   = 6371.0_DP
    complex(DP), parameter, public :: EI        = (0.0_DP, 1.0_DP)

    public :: std__countline
    public :: std__deg2rad
    public :: std__rad2deg
    public :: std__extend_array

    interface std__deg2rad
        module procedure d2r_s, d2r_d
    end interface std__deg2rad

    interface std__rad2deg
        module procedure r2d_s, r2d_d
    end interface std__rad2deg

    interface std__extend_array
        module procedure extend_array_c, extend_array_d, extend_array_f, extend_array_i
    end interface


contains

    subroutine extend_array_f(array, val)

        real(SP), intent(inout), allocatable :: array(:)
        real(SP), intent(in) :: val

        real(SP), allocatable :: tmp(:)
        integer :: n

        n = size(array)

        allocate (tmp(n + 1))
        tmp(1:n) = array(1:n)
        deallocate (array)
        allocate (array(n + 1))
        array(1:n) = tmp(1:n)
        array(n + 1) = val
        deallocate (tmp)

    end subroutine extend_array_f


    subroutine extend_array_d(array, val)

        real(DP), intent(inout), allocatable :: array(:)
        real(DP), intent(in) :: val

        real(DP), allocatable :: tmp(:)
        integer :: n

        n = size(array)

        allocate (tmp(n + 1))
        tmp(1:n) = array(1:n)
        deallocate (array)
        allocate (array(n + 1))
        array(1:n) = tmp(1:n)
        array(n + 1) = val
        deallocate (tmp)

    end subroutine extend_array_d


    subroutine extend_array_i(array, val)
        integer, intent(inout), allocatable :: array(:)
        integer, intent(in) :: val
        integer, allocatable :: tmp(:)
        integer :: n

        n = size(array)

        allocate (tmp(n + 1))
        tmp(1:n) = array(1:n)
        deallocate (array)
        allocate (array(n + 1))
        array(1:n) = tmp(1:n)
        array(n + 1) = val
        deallocate (tmp)

    end subroutine extend_array_i


    subroutine extend_array_c(clen, array, val)

        integer, intent(in) :: clen
        character(clen), intent(inout), allocatable :: array(:)
        character(clen), intent(in) :: val

        character(clen), allocatable :: tmp(:)
        integer :: n

        
        n = size(array)

        allocate (tmp(n + 1))
        tmp(1:n) = array(1:n)
        deallocate (array)
        allocate (array(n + 1))
        array(1:n) = tmp(1:n)
        array(n + 1) = val
        deallocate (tmp)

    end subroutine extend_array_c


    function d2r_d(deg)

        real(DP) :: d2r_d
        real(DP), intent(in) :: deg

        d2r_d = PI / 180.0_DP * deg

    end function d2r_d


    function d2r_s(deg)

        real(SP) :: d2r_s
        real(SP), intent(in) :: deg

        d2r_s = real(PI / 180.0_SP * deg)

    end function d2r_s


    function r2d_d(rad)

        real(DP) :: r2d_d
        real(DP), intent(in) :: rad

        r2d_d = 180.0_DP / PI * rad

    end function r2d_d


    function r2d_s(rad)

        real(SP) :: r2d_s
        real(SP), intent(in) :: rad

        r2d_s = real(180.0_SP / PI * rad)

    end function r2d_s


    subroutine std__countline(io, n, comment)

        !! Count line nubmer n included in the file specified by io

        integer, intent(in)  :: io
        integer, intent(out) :: n
        character, intent(in), optional :: comment

        integer :: stat
        character(256) :: line

        n = 0
        rewind (io)
        do
            read (io, '(A256)', iostat=stat) line
            if (stat /= 0) exit

            if (present(comment)) then
                if (line(1:1) == comment) cycle
            end if

            if (trim(line) /= '') then
                n = n + 1
            end if

        end do
        rewind (io)

    end subroutine std__countline
    

end module m_std
