module m_pwatch

    !! Stopwatch module for MPI parallel computation environment
    !!
    !! Copyright 2013-2025 Takuto Maeda. All rights reserved. This project is released under the MIT license.

    use iso_fortran_env, only: error_unit
    use m_std
    use mpi
    implicit none
    private
    save

    public :: pwatch__setup
    public :: pwatch__on
    public :: pwatch__off
    public :: pwatch__report

    integer, parameter :: NBLOCK_MAX = 500 !< maximum number of stopwatch blocks
    integer            :: irank
    logical            :: measure_time
    integer            :: nblock
    integer            :: iblock(NBLOCK_MAX)
    real               :: tim(NBLOCK_MAX)
    real               :: tim_total = 0
    character(255)     :: block_name(NBLOCK_MAX)
    integer            :: c0(NBLOCK_MAX)
    integer            :: count_rate
    integer            :: count_max
    integer            :: nproc !< number of processors
    integer            :: myid  !< my MPI rank

contains

    subroutine pwatch__setup(sw)
      !! Setup

        logical, intent(in) :: sw

        integer :: i, ierr

        measure_time = sw
        irank = 0
        tim_total = 0
        nblock = 0

        !! initializze
        do i = 1, NBLOCK_MAX
            tim(i) = 0
        end do

        c0(:) = 0

        call mpi_comm_size(mpi_comm_world, nproc, ierr)
        call mpi_comm_rank(mpi_comm_world, myid, ierr)

    end subroutine pwatch__setup


    subroutine pwatch__on(name)

        !! stopwatch switch on, labeled by "name"

        character(*), intent(in) :: name

        integer :: i
        integer :: cc, ib

        if (.not. measure_time) return

        irank = irank + 1

        !! table search by name
        iblock(irank) = 0
        do i = 1, nblock
            if (trim(block_name(i)) == trim(name)) then
                iblock(irank) = i
                exit
            end if
        end do
        if (iblock(irank) == 0) then ! generate new block
            nblock = nblock + 1
            iblock(irank) = nblock
            block_name(iblock(irank)) = trim(name)
        end if

        !! measure time
        call system_clock(cc, count_rate, count_max)

        if (irank > 1) then
            ib = iblock(irank - 1)

            !! First stop the time measurement of previous block
            if (cc >= c0(ib)) then
                tim(ib) = tim(ib) + (cc - c0(ib)) / real(count_rate)
            else
                tim(ib) = tim(ib) + (count_max + cc - c0(ib)) / real(count_rate)
            end if

        end if

        !! Then start measuring time of the current block
        ib = iblock(irank)
        c0(ib) = cc

    end subroutine pwatch__on


    subroutine pwatch__off(name)

        !! stopwatch switch off

        character(*), intent(in) :: name
        integer :: cc
        integer :: ib

        if (.not. measure_time) return

        call system_clock(cc, count_rate, count_max)

        !! search block name
        if (trim(block_name(iblock(irank))) == trim(name)) then
            ib = iblock(irank)

            if (cc >= c0(ib)) then
                tim(ib) = tim(ib) + (cc - c0(ib)) / real(count_rate)
            else
                tim(ib) = tim(ib) + (count_max + cc - c0(ib)) / real(count_rate)
            end if

            !! stop if it is the highest rank
            if (irank == 1) then
                return
            else !! check the upper rank
                c0(iblock(irank - 1)) = cc
                irank = irank - 1
            end if
        else
            write (error_unit, '(A)') 'ERROR [pwatch__off]: name '//trim(name) &
                //' does not match with '//trim(block_name(iblock(irank)))
        end if

    end subroutine pwatch__off


    subroutine pwatch__report(io, ionode)

        !! Report the total computation time and their occupation rate to unit io

        integer, intent(in) :: io     !! io file number
        integer, intent(in) :: ionode !! output node

        real, allocatable :: tsum(:), tbuf(:)
        integer :: i, j
        real :: trate(NBLOCK_MAX)
        real(SP), allocatable :: buf(:), tim0(:, :), buf2(:)
        integer :: ierr

        if (.not. measure_time) return

        allocate (tsum(0:nproc - 1), tbuf(0:nproc - 1))
        tbuf(0:nproc - 1) = 0
        tbuf(myid) = sum(tim(1:nblock))
        call mpi_reduce(tbuf(0:nproc - 1), tsum(0:nproc - 1), nproc, MPI_REAL, MPI_SUM, ionode, mpi_comm_world, ierr)

        allocate (tim0(0:nproc - 1, nblock), buf(0:nproc - 1), buf2(0:nproc - 1))
        tim0(:, :) = 0.0
        tim0(myid, 1:nblock) = tim(1:nblock)

        !! MPI data collection
        do i = 1, nblock
            buf = 0.0
            buf(myid) = tim0(myid, i)
            call mpi_reduce(buf(0:nproc - 1), buf2(0:nproc - 1), nproc, MPI_REAL, MPI_SUM, ionode, mpi_comm_world, ierr)
            tim0(0:nproc - 1, i) = buf2(0:nproc - 1)
        end do

        if (myid == ionode) then

            write (io, '(A)') &
                '#   CPU     #ID       Procedure Name         Real Time[s]   Total Time[s]   Occupancy[%]  Total Occp.[%] '
            write (io, '(A)') &
                '# -------+-------+-------------------------+--------------+--------------+--------------+----------------'

            do j = 0, nproc - 1
                trate(1:nblock) = tim0(j, 1:nblock) / tsum(j) * 100.0
                do i = 1, nblock
                    write (io, '(I8.5,I8.5,"    ",A22,4F15.3)') &
                        j, i, block_name(i), tim0(j, i), sum(tim0(j, 1:i)), trate(i), sum(trate(1:i))
                end do
            end do

        end if

    end subroutine pwatch__report

    
end module m_pwatch
