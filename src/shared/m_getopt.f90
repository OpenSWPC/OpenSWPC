module m_getopt

    !! Obtain command-line option
    !!
    !! Copyright 2013-2024 Takuto Maeda. All rights reserved. This project is released under the MIT license.

    use iso_fortran_env, only: error_unit
    use m_std
    use m_system

    implicit none
    private
    public :: getopt

    interface getopt

        module procedure getopt_d, & 
            getopt_s, & 
            getopt_i, & 
            getopt_c    

    end interface getopt

contains

    subroutine getopt_c(opt, isExist, val, default)

        !! Obtain command line as characters
        !!
        !! #### Example
        !! For the command "hoge -A -T 01"
        !!   - call getopt( 'A', isA, Aval ) -> isA=.true.,  Aval=''
        !!   - call getopt( 'T', isT, Tval ) -> isT=.true.,  Tval='01'
        !!   - call getopt( 'C', isC, Cval ) -> isC=.false., Cval=''
        
        character(*), intent(in)  :: opt
        logical, intent(out) :: isExist
        character(*), optional, intent(out) :: val
        character(*), optional, intent(in)  :: default

        integer :: narg
        character(256), allocatable :: argv(:)
        integer :: i
        character(256) :: optkey


        narg = system__iargc()
        allocate (argv(1:narg))

        do i = 1, narg
            call system__getarg(i, argv(i) (:))
        end do

        if (present(val)) then
            val = ''
        end if

        isExist = .false.
        optkey = '-'//trim(adjustl(opt))

        do i = 1, narg

            if (trim(optkey) == trim(argv(i))) then

                if (isExist) then
                    write (error_unit, *) 'getopt: ', &
                        'option '//trim(optkey)//' is multiplly defined. '
                end if

                isExist = .true.

                if (present(val)) then

                    if (i == narg) then
                        val = ''
                    else
                        val = argv(i + 1)
                    end if

                end if
            end if
        end do

        deallocate (argv)

        if (.not. isExist .and. present(default)) then
            val = default
        end if

    end subroutine getopt_c

    subroutine getopt_d(opt, isExist, val, default)

        !! Obtain command line option for double precision number

        character(*), intent(in)  :: opt
        logical, intent(out) :: isExist
        real(DP), intent(out) :: val
        real(DP), optional, intent(in)  :: default

        character(1024) :: aval

        call getopt_c(opt, isExist, aval)

        if (.not. isExist) then
            if (present(default)) then
                val = default
            else
                val = -99999.9_DP
            end if
            return
        end if

        read (aval, *) val

    end subroutine getopt_d

    subroutine getopt_s(opt, isExist, val, default)

        !! Obtain command line option for single precision number

        character(*), intent(in)  :: opt
        logical, intent(out) :: isExist
        real, intent(out) :: val
        real, optional, intent(in)  :: default

        character(1024) :: aval

        call getopt_c(opt, isExist, aval)

        if (.not. isExist) then
            if (present(default)) then
                val = default
            else
                val = -99999.9
            end if
            return
        end if

        read (aval, *) val

    end subroutine getopt_s

    subroutine getopt_i(opt, isExist, val, default)

        !! Obtain command line option for integer number

        character(*), intent(in)  :: opt
        logical, intent(out) :: isExist
        integer, intent(out) :: val
        integer, optional, intent(in)  :: default

        character(1024) :: aval

        call getopt_c(opt, isExist, aval)

        if (.not. isExist) then
            if (present(default)) then
                val = default
            else
                val = -99999
            end if
            return
        end if

        read (aval, *) val

    end subroutine getopt_i

end module m_getopt
