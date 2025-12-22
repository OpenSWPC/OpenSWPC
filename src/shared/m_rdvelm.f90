#include "../shared/m_debug.h"
module m_rdvelm

    !! Read 2D/3D velocity model from netcdf
    !!
    !! Jacob  This project is released under the MIT license.

    use m_std
    use m_debug
    use netcdf
    implicit none

    private
    public :: rdvelm__2d
    public :: rdvelm__3d

contains

    subroutine rdvelm__2d(ib, ie, kb, ke, fn_velm, mu, rho, lambda, Qp, Qs)

        integer, intent(in)  :: ib, ie
        integer, intent(in)  :: kb, ke
        character(*), intent(in)  :: fn_velm
        real(SP), intent(out) :: mu(kb:ke, ib:ie)
        real(SP), intent(out) :: rho(kb:ke, ib:ie)
        real(SP), intent(out) :: lambda(kb:ke, ib:ie)
        real(SP), intent(out) :: Qp(kb:ke, ib:ie)
        real(SP), intent(out) :: Qs(kb:ke, ib:ie)

        real(SP), allocatable :: hhmu(:, :), hhrho(:,:), hhlam(:,:), hhQp(:,:), hhQs(:,:)
        integer :: i, k, ii, kk
        integer :: ncid, muid, rhoid, lamid, qpid, qsid
        character(80) :: xn, zn
        integer :: nxc, nzc !< netcdf volume size

        
        call debug(fn_velm)
        call assert(nf90_open(fn_velm, NF90_NOWRITE, ncid) == NF90_NOERR)
        !! size
        call assert(nf90_inquire_dimension(ncid, 1, xn, nxc) == NF90_NOERR)
        call assert(nf90_inquire_dimension(ncid, 2, zn, nzc) == NF90_NOERR)
        call assert(nf90_inq_varid(ncid, "mu", muid) == NF90_NOERR)
        call assert(nf90_inq_varid(ncid, "rho", rhoid) == NF90_NOERR)
        call assert(nf90_inq_varid(ncid, "lambda", lamid) == NF90_NOERR)
        call assert(nf90_inq_varid(ncid, "Qp", qpid) == NF90_NOERR)
        call assert(nf90_inq_varid(ncid, "Qs", qsid) == NF90_NOERR)
        

        call debug(muid)
        call debug(rhoid)
        call debug(lamid)
        call debug(qpid)
        call debug(qsid)
        call debug(ncid)

        allocate (hhmu(nxc, nzc))
        call assert(nf90_get_var(ncid, muid, hhmu) == NF90_NOERR)
        allocate (hhrho(nxc, nzc))
        call assert(nf90_get_var(ncid, rhoid, hhrho) == NF90_NOERR)
        allocate (hhlam(nxc, nzc))
        call assert(nf90_get_var(ncid, lamid, hhlam) == NF90_NOERR)
        allocate (hhQp(nxc, nzc))
        call assert(nf90_get_var(ncid, qpid, hhQp) == NF90_NOERR)
        allocate (hhQs(nxc, nzc))
        call assert(nf90_get_var(ncid, qsid, hhQs) == NF90_NOERR)

        do k = kb, min(ke, nzc)

            if (k <= 0) then
                kk = k + nzc
            else
                kk = k
            end if

            do i = ib, ie
                ii = mod(i, nxc)
                if (ii <= 0) ii = ii + nxc
                mu(k, i) = hhmu(ii, kk)
                rho(k, i) = hhrho(ii, kk)
                lambda(k, i) = hhlam(ii, kk)
                Qp(k, i) = hhQp(ii, kk)
                Qs(k, i) = hhQs(ii, kk)
            end do
        end do

        deallocate (hhmu)
        deallocate (hhrho)
        deallocate (hhlam)
        deallocate (hhQp)
        deallocate (hhQs)
        

        !! bottom cyclic part
        do k = nzc + 1, ke
            kk = mod(k, nzc)
            if (kk <= 0) kk = kk + nzc
            mu(k, ib:ie) = mu(kk, ib:ie)
            rho(k, ib:ie) = rho(kk, ib:ie)
            lambda(k, ib:ie) = lambda(kk, ib:ie)
            Qp(k, ib:ie) = Qp(kk, ib:ie)
            Qs(k, ib:ie) = Qs(kk, ib:ie)
        end do

    end subroutine rdvelm__2d

    
    subroutine rdvelm__3d(ib, ie, jb, je, kb, ke, fn_velm, mu, rho, lambda, Qp, Qs)

        integer, intent(in)  :: ib, ie
        integer, intent(in)  :: jb, je
        integer, intent(in)  :: kb, ke
        character(*), intent(in)  :: fn_velm
        real(SP), intent(out) :: mu(kb:ke, ib:ie, jb:je)
        real(SP), intent(out) :: rho(kb:ke, ib:ie, jb:je)
        real(SP), intent(out) :: lambda(kb:ke, ib:ie, jb:je)
        real(SP), intent(out) :: Qp(kb:ke, ib:ie, jb:je)
        real(SP), intent(out) :: Qs(kb:ke, ib:ie, jb:je)

        real(SP), allocatable :: hhmu(:, :), hhrho(:,:), hhlam(:,:), hhQp(:,:), hhQs(:,:)
        integer :: i, j, k, ii, jj, kk
        integer :: ncid, muid, rhoid, lamid, qpid, qsid
        character(80) :: xn, yn, zn
        integer :: nxc, nyc, nzc !< netcdf volume size
        integer :: st(3), ct(3)

        call assert(nf90_open(fn_velm, NF90_NOWRITE, ncid) == NF90_NOERR)
        !! size
        call assert(nf90_inquire_dimension(ncid, 1, xn, nxc) == NF90_NOERR)
        call assert(nf90_inquire_dimension(ncid, 2, yn, nyc) == NF90_NOERR)
        call assert(nf90_inquire_dimension(ncid, 3, zn, nzc) == NF90_NOERR)
        call assert(nf90_inq_varid(ncid, "mu", muid) == NF90_NOERR)
        call assert(nf90_inq_varid(ncid, "rho", rhoid) == NF90_NOERR)
        call assert(nf90_inq_varid(ncid, "lambda", lamid) == NF90_NOERR)
        call assert(nf90_inq_varid(ncid, "Qp", qpid) == NF90_NOERR)
        call assert(nf90_inq_varid(ncid, "Qs", qsid) == NF90_NOERR)
        
        allocate (hhmu(nxc, nyc))
        allocate (hhrho(nxc, nyc))
        allocate (hhlam(nxc, nyc))
        allocate (hhQp(nxc, nyc))
        allocate (hhQs(nxc, nyc))
        
        st(1:3) = (/1, 1, 1/)
        ct(1:3) = (/nxc, nyc, 1/)

        do k = kb, min(ke, nzc)

            if (k <= 0) then
                kk = k + nzc
            else
                kk = k
            end if

            st(3) = kk
            call assert(nf90_get_var(ncid, muid,  hhmu,  start=st, count=ct) == NF90_NOERR)
            call assert(nf90_get_var(ncid, rhoid, hhrho, start=st, count=ct) == NF90_NOERR)
            call assert(nf90_get_var(ncid, lamid, hhlam, start=st, count=ct) == NF90_NOERR)
            call assert(nf90_get_var(ncid, qpid,  hhQp,  start=st, count=ct) == NF90_NOERR)
            call assert(nf90_get_var(ncid, qsid,  hhQs,  start=st, count=ct) == NF90_NOERR)

            do j = jb, je

                jj = mod(j, nyc)
                if (jj <= 0) jj = jj + nyc

                do i = ib, ie
                    ii = mod(i, nxc)
                    if (ii <= 0) ii = ii + nxc
                    mu(k, i, j) = hhmu(ii, jj)
                    rho(k, i, j) = hhrho(ii, jj)
                    lambda(k, i, j) = hhlam(ii, jj)
                    Qp(k, i, j) = hhQp(ii, jj)
                    Qs(k, i, j) = hhQs(ii, jj)
                end do
            end do

        end do

        deallocate (hhmu)
        deallocate (hhrho)
        deallocate (hhlam)
        deallocate (hhQp)
        deallocate (hhQs)

        !! bottom cyclic part
        do k = nzc + 1, ke
            kk = mod(k, nzc)
            mu(k, ib:ie, jb:je) = mu(kk, ib:ie, jb:je)
            rho(k, ib:ie, jb:je) = rho(kk, ib:ie, jb:je)
            lambda(k, ib:ie, jb:je) = lambda(kk, ib:ie, jb:je)
            Qp(k, ib:ie, jb:je) = Qp(kk, ib:ie, jb:je)
            Qs(k, ib:ie, jb:je) = Qs(kk, ib:ie, jb:je)
        end do

    end subroutine rdvelm__3d

end module m_rdvelm
