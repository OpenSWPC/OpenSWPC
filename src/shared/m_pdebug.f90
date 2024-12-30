module m_debug

    !! Debug routines under MPI environment.
    !!
    !! Use with preprocessor macro (#include "m_pdebug.h" at the top of the soruce code) will give richer information.
    !! The following routines must be called after mpi_init.
    !!
    !! #### Usage
    !! - call debug(var):     show variable var.
    !! - call assert( cond ): abort if cond = .false. . cond must be logical value or condition.
    !! - call info( msg ):    write message msg to error_unit.
    !!
    !! Copyright 2013-2024 Takuto Maeda. All rights reserved. This project is released under the MIT license.

    use iso_fortran_env, only: error_unit
    use m_std
    use mpi
    implicit none
    private

    public :: debug, debug__macro
    public :: info, info__macro
    public :: assert, assert__macro

    interface debug
        module procedure debug_c0, debug_i0, debug_r0, debug_d0, debug_l0, debug__void
    end interface debug

    interface debug__macro
        module procedure debug_c, debug_i, debug_r, debug_d, debug_l, debug__void
    end interface debug__macro

    logical :: initialized = .false.
    integer :: nproc, myid

contains

    subroutine debug_init
        integer :: ierr
        call mpi_comm_size(mpi_comm_world, nproc, ierr)
        call mpi_comm_rank(mpi_comm_world, myid, ierr)
        initialized = .true.
    end subroutine debug_init


    subroutine debug__void()

        !! Do Nothing: dummy

        return

    end subroutine debug__void


    subroutine debug_c(var, fname, nline)

        character(*), intent(in) :: var
        character(*), intent(in) :: fname
        integer, intent(in) :: nline

        character(5) :: cline

        if (.not. initialized) call debug_init
        write (cline, '(I5)') nline

        write (error_unit, '(A,I0,A)') '[debug{', myid, '}] '//fname//' ('//trim(adjustl(cline))//'):  '//trim(adjustl(var))

    end subroutine debug_c


    subroutine debug_c0(var)

        character(*), intent(in) :: var

        if (.not. initialized) call debug_init
        write (error_unit, '(A,I0,A)') '[debug{', myid, '}] '//trim(adjustl(var))

    end subroutine debug_c0


    subroutine debug_i(var, fname, nline)

        integer, intent(in) :: var
        character(*), intent(in) :: fname
        integer, intent(in) :: nline

        character(10) :: cvar

        write (cvar, '(I10)') var
        call debug_c(cvar, fname, nline)

    end subroutine debug_i


    subroutine debug_i0(var)

        integer, intent(in) :: var

        character(10) :: cvar

        write (cvar, '(I10)') var
        call debug_c0(cvar)

    end subroutine debug_i0


    subroutine debug_r(var, fname, nline)

        real, intent(in) :: var
        character(*), intent(in) :: fname
        integer, intent(in) :: nline

        character(15) :: cvar

        if (abs(var) < 10000.) then
            write (cvar, '(F15.5)') var
        else
            write (cvar, '(ES15.5)') var
        end if
        call debug_c(cvar, fname, nline)

    end subroutine debug_r


    subroutine debug_r0(var)

        real, intent(in) :: var

        character(15) :: cvar

        if (abs(var) < 10000.) then
            write (cvar, '(F15.5)') var
        else
            write (cvar, '(ES15.5)') var
        end if
        call debug_c0(cvar)

    end subroutine debug_r0


    subroutine debug_d(var, fname, nline)

        real(DP), intent(in) :: var
        character(*), intent(in) :: fname
        integer, intent(in) :: nline

        character(15) :: cvar
        
        if (abs(var) < 10000.) then
            write (cvar, '(F15.5)') var
        else
            write (cvar, '(ES15.5)') var
        end if
        call debug_c(cvar, fname, nline)

    end subroutine debug_d


    subroutine debug_d0(var)

        real(DP), intent(in) :: var

        character(15) :: cvar

        if (abs(var) < 10000.) then
            write (cvar, '(F15.5)') var
        else
            write (cvar, '(ES15.5)') var
        end if
        call debug_c0(cvar)

    end subroutine debug_d0


    subroutine debug_l(var, fname, nline)

        logical, intent(in) :: var
        character(*), intent(in) :: fname
        integer, intent(in) :: nline

        character(15) :: cvar

        if (var) then
            cvar = '.true.'
        else
            cvar = '.false.'
        end if
        call debug_c(cvar, fname, nline)

    end subroutine debug_l


    subroutine debug_l0(var)

        logical, intent(in) :: var

        character(15) :: cvar

        if (var) then
            cvar = '.true.'
        else
            cvar = '.false.'
        end if
        call debug_c0(cvar)

    end subroutine debug_l0


    subroutine assert(cond)

        logical, intent(in) :: cond
        integer :: ierr, errcode

        if (.not. initialized) call debug_init
        if (.not. cond) then
            write (error_unit, '(A,I0,A)') '[assert{', myid, '}] failed'
            call mpi_abort(mpi_comm_world, errcode, ierr)
        end if

    end subroutine assert


    subroutine assert__macro(cond, fname, nline)

        logical, intent(in) :: cond
        character(*), intent(in) :: fname
        integer, intent(in) :: nline

        character(5) :: cl
        integer :: ierr, errcode

        if (.not. initialized) call debug_init
        if (.not. cond) then
            write (cl, '(I5)') nline

            write (error_unit, '(A,I0,A)') '[assert{', myid, '}] failed at '//trim(adjustl(fname))//'('//trim(adjustl(cl))//')'
            call mpi_abort(mpi_comm_world, errcode, ierr)
        end if

    end subroutine assert__macro


    subroutine info(msg)

        character(*), intent(in) :: msg

        if (.not. initialized) call debug_init
        write (error_unit, '(A,I0,A)') '[info{', myid, '}] '//trim(adjustl(msg))

    end subroutine info


    subroutine info__macro(msg, fname, nline)

        character(*), intent(in) :: msg
        character(*), intent(in) :: fname
        integer, intent(in) :: nline

        if (.not. initialized) call debug_init
        write (error_unit, '(A,I0,A,I0,A)') '[info{', myid, '}] '//trim(adjustl(fname))//'(', nline, '): '//trim(adjustl(msg))

    end subroutine info__macro

    
end module m_debug
