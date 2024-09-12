#include "../shared/m_debug.h"
module m_report

    !! Terminal/logfile report
    !!
    !! Copyright 2013-2024 Takuto Maeda. All rights reserved. This project is released under the MIT license.

    use iso_fortran_env, only: error_unit
    use m_std
    use m_debug
    use m_global
    use m_pwatch
    use m_kernel
    use m_readini
    use m_version

    implicit none
    private
    save

    public :: report__setup
    public :: report__progress
    public :: report__terminate

    integer :: ntdec_r

    integer, parameter :: terminal_output_node = 0

    !<< Lapse Time Measurement >>
    integer  :: timcount, timcount0, timprev
    real(SP) :: ttotal

contains

    subroutine report__setup(io_prm)

        !! Initialization, welcome message to terminal, open logfile

        integer, intent(in) :: io_prm

        integer :: crate
        real(SP) :: mem_all, mem_node, r, c
        character(256) :: codename, ver

        call readini(io_prm, 'ntdec_r', ntdec_r, 10)

        if (myid == terminal_output_node) then

            call version__get(ver)
            codename = "  SWPC_PSV version "//trim(ver)
            if (benchmark_mode) then; codename = trim(codename)//' (benchmark mode)   '
            else if (pw_mode) then; codename = trim(codename)//' (plane wave mode)  '
            else if (bf_mode) then; codename = trim(codename)//' (body force mode)  '
            end if

            write (error_unit, '(A)')
            write (error_unit, '(A)') " ------------------------------------------------------------------------------"
            write (error_unit, '(A)') trim(codename)
            write (error_unit, '(A)') " ------------------------------------------------------------------------------"
            write (error_unit, '(A)')

        end if

        call memory_size_psv(nproc_x, nx, nz, nm, na, mem_all, mem_node)
        call fdm_cond_stability(real(dx), 1e10, real(dz), vmax, dt, c)
        call fdm_cond_wavelength(real(dx), -1.0, real(dz), vmin, fmax, r)

        if (myid == terminal_output_node) then

            write (error_unit, *)
            write (error_unit, '(A,I8,A,I6)') "  Grid Size               : ", nx, " x ", nz
            write (error_unit, '(A,I15    )') "  MPI Partitioning        : ", nproc_x
            write (error_unit, '(A,F15.3,A)') "  Total Memory Size       : ", mem_all, "  [GiB]"
            write (error_unit, '(A,F15.3,A)') "  Node Memory Size        : ", mem_node, "  [GiB]"
            write (error_unit, '(A,F15.3,A)') "  Stability  Condition c  : ", c, "  (c<1)"
            write (error_unit, '(A,F15.3,A)') "  Wavelength Condition r  : ", r, "  (r>5-10)"
            write (error_unit, '(A,F15.3,A)') "  Minimum velocity        : ", vmin, "  [km/s]"
            write (error_unit, '(A,F15.3,A)') "  Maximum velocity        : ", vmax, "  [km/s]"
            write (error_unit, '(A,F15.3,A)') "  Maximum frequency       : ", fmax, "  [Hz]"
            write (error_unit, *)
            write (error_unit, '(A)') " ------------------------------------------------------------------------------"
            write (error_unit, *)

            if (r < 5) then
                call info('wavelength condition is violated! ')
                call info('use smaller grid and/or decrease maximum frequency')
            end if

            if (c < 0.5) then
                call info('time step is too small!')
                call info('consider increase time step up to twice')
            end if

            if (c > 1.0) then
                call info('stability condition is violated!')
                call info('use smaller time step and/or decrease max velocity')
                call assert(c <= 1.0)
            end if

        end if

        !! Initialize elapsed time counter
        if (myid == terminal_output_node) then

            call system_clock(timcount, crate)
            timcount0 = timcount
            timprev = timcount
            ttotal = 0

        end if

    end subroutine report__setup

    subroutine report__progress(it)

        !! Show progres to the terminal

        integer, intent(in) :: it
        real(SP) :: vxm, vzm
        real(SP) :: vxa, vza
        integer  :: ierr
        real(SP) :: etas
        integer  :: etah, etam, etasi
        real(SP) :: tstep
        integer  :: crate, cmax

        if (mod(it, ntdec_r) /= 0) return

        call pwatch__on("report__progress")

        call kernel__vmax(vxm, vzm)

        call mpi_reduce(vxm, vxa, 1, MPI_REAL, MPI_MAX, terminal_output_node, mpi_comm_world, ierr)
        call mpi_reduce(vzm, vza, 1, MPI_REAL, MPI_MAX, terminal_output_node, mpi_comm_world, ierr)

        if (myid == terminal_output_node) then

            !! to SI unit
            vxa = vxa * UC * M0
            vza = vza * UC * M0

            !! eta count
            call system_clock(timcount, crate, cmax)
            if (timcount >= timprev) then
                tstep = real(timcount - timprev) / real(crate)
            else
                tstep = real(cmax + timcount - timprev) / real(crate)
            end if

            ttotal = ttotal + tstep

            etas = real(nt - it) / real(it) * ttotal

            etah = int(etas / (60 * 60)); etas = etas - etah * 60 * 60
            etam = int(etas / (60)); etas = etas - etam * 60
            etasi = int(etas)
            timprev = timcount

            write (error_unit, '(A,I7.7,  A,F6.3,A,   A,I3.3,A,I2.2,A,I2.2,A, 2(ES9.2,A))') &
                "  it=", it, ",", &
                ttotal / it, " s/loop,", &
                " eta ", etah, ":", etam, ":", etasi, ", (", vxa, " ", vza, " )"

        end if

        call pwatch__off("report__progress")

    end subroutine report__progress

    subroutine report__terminate

        if (myid == terminal_output_node) then

            write (error_unit, *)
            write (error_unit, '(A)') " ------------------------------------------------------------------------------"
            write (error_unit, *) ""
            write (error_unit, '(A,F15.3,A)') "  Total time             : ", ttotal, " s"
            write (error_unit, *)
            write (error_unit, '(A)') " ------------------------------------------------------------------------------"

        end if

    end subroutine report__terminate

end module m_report
