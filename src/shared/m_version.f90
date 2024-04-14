!! ----------------------------------------------------------------------------------------------------------------------------- !!
!> Code Version
!!
!! @copyright
!!   Copyright 2013-2024 Takuto Maeda. All rights reserved.
!!   This project is released under the MIT license.
!! ----------------------------------------------------------------------------------------------------------------------------- !!
module m_version

    use iso_fortran_env, only: output_unit
    implicit none
    public 
    !! Do not modify this unless version up
    character(80), public :: version = '5.3.1'

contains

    subroutine version__display(codename)

        character(*), intent(in), optional :: codename

        if (present(codename)) then
            write (output_unit, '(a)') trim(codename)//' (OpenSWPC) version '//VERSION
        else
            write (output_unit, '(a)') 'OpenSWPC version '//VERSION
        end if

        stop

    end subroutine version__display

    subroutine version__get(ver)
        character(*), intent(out) :: ver
        ver = version
    end subroutine version__get

end module m_version