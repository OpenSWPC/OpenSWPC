module m_system

    !! System routines, for processig command line argument, environment variables, and system call
    !!
    !! Copyright 2013-2024 Takuto Maeda. All rights reserved. This project is released under the MIT license.

    use m_std

    implicit none
    private

    public :: system__getarg
    public :: system__getenv
    public :: system__call
    public :: system__iargc
    public :: system__expenv

    interface system__getarg
        module procedure getarg_a, getarg_i, getarg_f, getarg_d
    end interface system__getarg

contains

    subroutine system__call(cmd)

        !! Do system-call.

        character(*), intent(in) :: cmd

        call execute_command_line(cmd)

    end subroutine system__call

    integer function system__iargc()

        !! Returns a number of arguments. Fortran2003 wrapper function

        system__iargc = command_argument_count()

    end function system__iargc

    subroutine system__getenv(name, value)

        !! Obtain environmental variable "name".

        character(*), intent(in)  :: name
        character(*), intent(out) :: value

        call get_environment_variable(name, value)

    end subroutine system__getenv

    subroutine getarg_a(i, arg)
        !! get i-th command line argument, Fortran2003 wrapper subroutine

        integer, intent(in)  :: i   !! order of the arguments
        character(*), intent(out) :: arg !! argument

        call get_command_argument(i, arg)

    end subroutine getarg_a

    subroutine getarg_i(i, arg)

        integer, intent(in) :: i
        integer, intent(out) :: arg

        character(256) :: carg

        call getarg_a(i, carg)
        read (carg, *) arg

    end subroutine getarg_i

    subroutine getarg_f(i, arg)

        integer, intent(in) :: i
        real, intent(out) :: arg

        character(256) :: carg

        
        call getarg_a(i, carg)
        read (carg, *) arg

    end subroutine getarg_f

    subroutine getarg_d(i, arg)

        integer, intent(in) :: i
        real(DP), intent(out) :: arg

        character(256) :: carg

        call getarg_a(i, carg)
        read (carg, *) arg

    end subroutine getarg_d

    subroutine system__expenv(str)

        !! Expand shell environmental variables wrapped in ${...}

        character(*), intent(inout) :: str
        character(256) :: str2
        integer :: ikey1, ikey2
        integer :: iptr
        character(256) :: str3

        iptr = 1
        str2 = ''
        do
            ikey1 = scan(str(iptr:), "${") + iptr - 1
            if (ikey1 == iptr - 1) exit

            ikey2 = scan(str(iptr:), "}") + iptr - 1
            str2 = trim(str2)//str(iptr:ikey1 - 1)
            call system__getenv(str(ikey1 + 2:ikey2 - 1), str3)
            str2 = trim(str2)//trim(str3)
            iptr = ikey2 + 1

        end do
        str2 = trim(str2)//trim(str(iptr:))

        str = trim(str2)

    end subroutine system__expenv

end module m_system
